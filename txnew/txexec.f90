module tx_main
  implicit none
  private
  integer(4), dimension(1:2) :: iasg
  real(8), dimension(:), allocatable :: L2
  real(8), dimension(:,:), allocatable :: XN, XP, ASG
  ! Anderson acceleration (NANDER /= 0): scaled iterate and residual of the last
  ! sweep, histories of their differences, workspace, and run statistics.
  real(8), dimension(:,:), allocatable :: and_dx, and_df, and_A
  real(8), dimension(:),   allocatable :: and_x0, and_f0, and_x, and_g, and_b, and_s, and_work
  real(8), dimension(:),   allocatable :: and_sc
  integer(4), dimension(:), allocatable :: and_iwork
  integer(4) :: and_k = 0, and_ncall = 0, and_nrst = 0, and_nfail = 0
  logical :: and_have0 = .false.
  real(8) :: and_condmax = 0.d0, and_cpu = 0.d0
  ! Quality report of the last DGBSVX band solve (MDLPCK = 3, and the MDLPCK = 2
  ! fallback).  Kept here so the diagnostics can print it without threading it
  ! through the call chain.
  real(8) :: solq_rcond = 0.d0, solq_ferr = 0.d0, solq_berr = 0.d0
  real(8) :: solq_berr_row = 0.d0, solq_berr_num = 0.d0, solq_berr_den = 0.d0
  ! How often the non-smooth guards actually bite.  Reported at the end of a RUN
  ! under IDIAG >= 5; zero unless something fired.
  integer(4) :: nsm_minus(1:4) = 0, nsm_tiny = 0, nsm_tiny_nq(1:64) = 0
  integer(4) :: solq_berr_idx = 0
  character(len=1) :: solq_equed = 'N'
  public :: TXEXEC

contains

!***************************************************************
!
!   MAIN ROUTINE
!
!***************************************************************

  subroutine TXEXEC
    use tx_commons, only : IERR, T_TX, AVE_IC
    use tx_interface, only : APTOS

    integer(4) :: NSTR1, NSTR2, NSTR3, NSTR4
    real(4) :: gCTIME1, gCTIME2, gCTIME3
    character(len=10) :: STR1, STR2, STR3, STR4
    real(8) :: t_interval

    if (IERR /= 0) then
       write(6,*) '### ERROR(TXEXEC) : Error should be cleared.'
       return
    end if

    call CPU_TIME(gCTIME1)
    t_interval = T_TX

    ! ***** Core part of simulation *****

    write(6,*) 'Calculating...'
    call TXLOOP

    ! ***********************************

    call CPU_TIME(gCTIME2)
    gCTIME3 = gCTIME2-gCTIME1
    t_interval = T_TX - t_interval

    NSTR1 = 0
    call APTOS(STR1, NSTR1, gCTIME3, 'F2')
    NSTR2 = 0
    call APTOS(STR2, NSTR2, real(t_interval), 'E1')
    if(gCTIME3 /= 0) then
       NSTR3 = 0
       call APTOS(STR3, NSTR3, real(t_interval) / gCTIME3 * 100, 'F3')
       NSTR4 = 0
       call APTOS(STR4, NSTR4, real(AVE_IC),'F2')
       write(6,*) 'CPU = ', STR1(1:NSTR1), ' (sec)   ',  &
            &     'sim time = ', STR2(1:NSTR2), ' (sec)', '  (', STR3(1:NSTR3), '%)  ', &
            &     'ICave = ', STR4(1:NSTR4)
    else
       write(6,*) 'CPU = ', STR1(1:NSTR1), ' (sec)   ',  &
            &     'sim time = ', STR2(1:NSTR2), ' (sec)'
    END IF

    !     ***** Print simulation results *****

    call TXWDAT

  end subroutine TXEXEC

!***************************************************************
!
!   Core routine of simulation
!
!***************************************************************

  subroutine TXLOOP
    use tx_commons, only : T_TX, NTCUM, rIPe, rIPs, NTMAX, IGBDF, NQMAX, NRMAX, X, ICMAX, &
         &                 ErV, PNsV_FIX, PTsV_FIX, &
         &                 ErV_FIX, IERR, LQb1, LQn1, LQn2, LQr1, &
         &                 tiny_cap, EPS, IDIAG, NTSTEP, ICCALC, iaxproj, &
         &                 FSRP, MODEAV, XOLD, bdf_history_ready, bdf_history_dt, &
         &                 NT, DT, rIP, MDLPCK, ICONT, AVE_IC, MODECV, oldmix, oldmixa, Var, RLXN0, MDLN0FL, LQn1, LQn3, NC0CAP, &
         &                 NANDER, NRETRY
    use tx_variables
    use tx_coefficients, only : TXCALA
    use tx_graphic, only : TX_GRAPH_SAVE, TXSTGT, TXSTGV, TXSTGR, TXSTGQ, &
         &                 NGR, NGRM, NGRSTP, NGTSTP, NGVSTP, GY, GT
    use tx_ntv, only : Wnm_spline
    use tx_glob, only : TXGLOB
    use libbnd, only: BANDRD
#if   _LAPACK == 0
    ! for self-compiled lapack
    use f95_lapack, only : GBSV => LA_GBSV
#elif _LAPACK == 1
    ! for intel mkl LAPACK95, 
    !  Note: This module file includes "ptsv" subroutine, 
    !        whose name conflicts with PTsV defined in TASK/TX.  
    use lapack95, only : GBSV
#endif

    real(8), dimension(:,:), allocatable :: BA, BL
    real(8), dimension(:),   allocatable :: BX
    integer(4) :: NR, NQ, IC = 0, IDIV, NTDO, ICSUM, IDIAGL, istat
    integer(4), save :: ntc0cap = 0
    ! --- retry of a non-converged step as two half steps (NRETRY = 1) ---
    integer(4) :: nsub, isub, nretried
    real(8) :: DTfull, Tstart
    real(8), dimension(:,:), allocatable :: Xstart
    ! --- adaptive under-relaxation ---------------------------------------
    !  Below dt ~ 3e-6 the fixed-point map for the radial flows stops
    !  contracting.  It does not blow up: the sweep increment simply sits at
    !  O(1) and stays there -- LQz2 at 1.5 from sweep 2 to sweep 60 at
    !  dt = 3.9e-7 -- so "the increment grew since the last sweep" cannot see
    !  it, and neither can a ratio to the running minimum, which is larger at
    !  the ordinary convergence floor than it is here.
    !
    !  ERRSWEEP is the increment normalised exactly as check_convergence
    !  normalises it, so it is already dimensionless: the largest nodal change
    !  in a sweep against that equation's profile RMS.  Past a few sweeps a
    !  healthy iteration has it far below one.  Measured over every case on
    !  disk that completes -- dt from 1e-4 to 7.8e-7, NRMAX 60 to 360, the
    !  circular equilibrium, a hundred times heavier impurity influx -- the
    !  largest value at sweep DIVIC or beyond is 0.145, while the two cases
    !  that die on a negative density reach 3.1 and 7.7.  DIVBIG sits in that
    !  gap.  One sweep over it is not enough, so the counter decrements on a
    !  good sweep rather than resetting, and damping engages at DIVCNT: with
    !  these settings no case that completes ever reaches 1, and the two that
    !  fail reach 14 and 53.
    !
    !  Normalising against the first sweep was tried and rejected.  At the
    !  start of a run LQm3 -- psi-dot, a divided difference -- can make sweep 1
    !  a thousand times the later ones, and a sustained increment of 7.5 then
    !  reads as 0.008 of it and goes unnoticed.
    real(8) :: errsweep, relax, avnq
    integer(4) :: ndiv, ndamp
    logical :: damping
    integer(4), parameter :: DIVIC = 8, DIVCNT = 6
    real(8), parameter :: DIVBIG = 0.5d0, DIVCAP = 0.9d0
    ! ---------------------------------------------------------------------
    ! *** LAPACK ************************************
    integer(4) :: ierr_la
#if _LAPACK == 2
    integer(4) :: m, kl, n, ku, nrhs, ldBL, ldBX
    integer(4), dimension(1:NQMAX*(NRMAX+1)) :: ipiv
#endif
    ! ***********************************************
    real(8) :: TIME0, DIP, EPSabs
    real(8), dimension(1:NQMAX) :: tiny_array
    character(len=80) :: MSG_NQ
#ifdef TX_PROFILE
    real(8) :: prof_total_t0, prof_total_t1
    real(8) :: prof_ic_t0, prof_ic_t1
    real(8) :: prof_t0, prof_t1
    real(8) :: prof_ic_sum, prof_txcala_sum, prof_txcalb_sum, prof_solver_sum, prof_io_sum
#endif
#ifdef TX_SOLVER_BENCH
    logical, save :: solver_bench_done = .false.
    integer(4) :: solver_bench_unit
    integer(4) :: m_bench, n_bench, kl_bench, ku_bench, ldbl_bench, ldbx_bench
    character(len=*), parameter :: solver_bench_matrix_file = 'solver_bench_matrix.bin'
    character(len=*), parameter :: solver_bench_solution_file = 'solver_bench_solution_ref.bin'
#endif
#ifdef TX_SOLVER_VERIFY
    logical, save :: solver_verify_done = .false.
    real(8), dimension(:,:), allocatable :: bl_verify
    real(8), dimension(:), allocatable :: bx_rhs_verify, bx_ref_verify
    real(8) :: verify_max_abs, verify_max_rel, verify_max_rel_floor, verify_scale, verify_diff
    integer(4) :: verify_i, ierr_verify
#if _LAPACK == 2
    integer(4), dimension(1:NQMAX*(NRMAX+1)) :: ipiv_verify
#endif
#endif

    allocate( BA(1:4*NQMAX-1,1:NQMAX*(NRMAX+1)) &
         &  , BL(1:6*NQMAX-2,1:NQMAX*(NRMAX+1)) &
         &  , BX(1:NQMAX*(NRMAX+1)))
    allocate( XN, XP, ASG, mold=X )
    allocate( L2(1:NQMAX) )

    IDIAGL = mod(IDIAG,10)
    EPSabs = abs(EPS)

    !  Read spline table for neoclassical toroidal viscosity if not loaded when FSRP/=0
    if(FSRP /= 0.d0) call Wnm_spline

    if (MODEAV == 0) then
       IDIV = NTMAX + 1
    else
       IDIV = NTMAX / MODEAV
    end if
    TIME0 = T_TX
    if(NTMAX /= 0) DIP = (rIPe - rIPs) / NTMAX
    ICSUM = 0
    ndamp = 0
    nretried = 0 ; nsub = 1
    and_ncall = 0 ; and_nrst = 0 ; and_nfail = 0 ; and_condmax = 0.d0 ; and_cpu = 0.d0
#ifdef TX_PROFILE
    call CPU_TIME(prof_total_t0)
    prof_ic_sum    = 0.d0
    prof_txcala_sum = 0.d0
    prof_txcalb_sum = 0.d0
    prof_solver_sum = 0.d0
    prof_io_sum     = 0.d0
#endif

    ! Save X -> XP -> XOLD for BDF only when there is no X(t-dt) to lean on.
    ! A RUN starts from scratch; a CONT continues with the history it has.
    if(T_TX == 0.d0) bdf_history_ready = .false.
    ! XOLD belongs to the step size it was taken with.  The weights used here are
    ! the fixed-step BDF2 ones, so carrying the history across a change of DT --
    ! which is what P followed by C does -- is not BDF2 at all.  Measured on the
    ! standard case, stepping to t = 1e-3 as 100 steps of 1e-7 then 79 of 1.25e-5
    ! lands 2.0e-4 from a fine reference when the history is carried over and
    ! 5.3e-6 when it is not, against 1.6e-5 for a uniform 1.25e-5.  Take one
    ! backward-Euler step at the change instead; its O(dt^2) local error, once,
    ! leaves the order alone.
    if(DT /= bdf_history_dt) bdf_history_ready = .false.
    if(IGBDF /= 0 .and. .not. bdf_history_ready) XOLD = X

    L_NTDO:do NTDO = 1, NTMAX
       NT = NTDO
       NTCUM = NTCUM + 1
       T_TX = TIME0 + DT*NT
       rIP  = rIPs  + DIP*NT
       nsub = 1 ; isub = 0

       L_SUB : do
       isub = isub + 1
       if(nsub > 1) then   ! the substeps of a retried step
          T_TX = Tstart + DT * isub
          rIP  = rIPs + DIP * (NT - 1 + real(isub,8) / nsub)
       end if

       ! Create new X := XN
       XN = X

       ! The divergence latch is per time step.
       ndiv    = 0
       damping = .false.
       relax   = oldmix

       ! Negligible order of magnitude for each variable
       do NQ = 1, NQMAX
          tiny_array(NQ) = maxval(XN(:,NQ)) * tiny_cap
       end do

       ! In the following loop, XN is being updated during iteration.
       L_IC : do IC = 1, ICMAX
#ifdef TX_PROFILE
          call CPU_TIME(prof_ic_t0)
#endif
          ! Save past X := XP
          XP = XN

          call TXCALV(XP)
          if(IC <= ICCALC) then
             PNsV_FIX(:,:) = Var(:,:)%n
             PTsV_FIX(:,:) = Var(:,:)%T
             ErV_FIX (:) = ErV (:)
          end if

          call TXCALC(IC)
#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          call TXCALA
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_txcala_sum = prof_txcala_sum + (prof_t1 - prof_t0)
#endif
          ! Get BA or BL, and BX
#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          call TXCALB(BA,BL,BX)
          if(IDIAGL >= 6 .and. MDLPCK /= 0) call report_conditioning(BL,NQMAX,NRMAX,NT,IC)
          if(IDIAGL >= 8 .and. MDLPCK /= 0) call dump_matrix(BL,BX,NQMAX,NRMAX,NT,IC)
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_txcalb_sum = prof_txcalb_sum + (prof_t1 - prof_t0)
#endif
!          call TXGLOB

#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          if(MDLPCK == 0) then
             call BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQMAX-1, IERR)
             if (IERR >= 30000) then
                write(6,'(3(A,I6))') '### ERROR(TXLOOP) : Matrix BA is singular at ',  &
                     &              NT, ' -', IC, ' step. IERR=',IERR
                IERR = 1
                XN = XP
                goto 180
             end if
          else if(MDLPCK == 2) then
#ifdef TX_SOLVER_VERIFY
             if(.not. solver_verify_done) then
                allocate(bl_verify(size(BL,1), size(BL,2)), bx_rhs_verify(size(BX,1)), bx_ref_verify(size(BX,1)))
                bl_verify = BL
                bx_rhs_verify = BX
             end if
#endif
             call TXSOLV_BLOCKTRI(BL, BX, NQMAX, NRMAX, ierr_la)
             if(ierr_la == -2) then
                ! Fallback when BL has non-zero far blocks: the equilibrated band
                ! solve, not the plain one, so the fallback keeps the accuracy the
                ! block path was chosen for.
                call TXSOLV_BAND_EQUIL(BL, size(BL,1), BX, NQMAX, NRMAX, ierr_la, IDIAGL >= 6)
                if(IDIAGL >= 6) call dump_solver_quality(NT,IC,'fallback')
             end if
             if(ierr_la /= 0) then
                write(6,'(3(A,I6))') '### ERROR(TXLOOP) : BLKTRI, NT = ',  &
                     &              NT, ' -', IC, ' step. IERR=',ierr_la
                IERR = 1
                XN = XP
                goto 180
             end if
#ifdef TX_SOLVER_VERIFY
             if(.not. solver_verify_done) then
                bx_ref_verify = bx_rhs_verify
#if   _LAPACK == 2
                m    = NQMAX*(NRMAX+1)
                kl   = 2*NQMAX-1
                n    = m
                ku   = kl
                nrhs = 1
                ldBL = 6*NQMAX-2
                ldBX = NQMAX*(NRMAX+1)
                call LAPACK_DGBSV(n,kl,ku,nrhs,bl_verify,ldBL,ipiv_verify,bx_ref_verify,ldBX,ierr_verify)
#else
                call GBSV(bl_verify,bx_ref_verify,INFO=ierr_verify)
#endif
                if(ierr_verify == 0) then
                   verify_max_abs = 0.d0
                   verify_max_rel = 0.d0
                   verify_max_rel_floor = 0.d0
                   do verify_i = 1, size(BX)
                      verify_diff = abs(BX(verify_i) - bx_ref_verify(verify_i))
                      verify_max_abs = max(verify_max_abs, verify_diff)
                      verify_scale = max(abs(BX(verify_i)), abs(bx_ref_verify(verify_i)))
                      if(verify_scale > 0.d0) then
                         verify_max_rel = max(verify_max_rel, verify_diff / verify_scale)
                         if(verify_scale > 1.d-3) verify_max_rel_floor = max(verify_max_rel_floor, verify_diff / verify_scale)
                      end if
                   end do
                   write(6,'(A,1X,A,1X,ES12.4,1X,A,1X,ES12.4,1X,A,1X,ES12.4,1X,A,1X,I6,1X,A,1X,I6)') &
                        & '## TX_SOLVER_VERIFY', 'max_abs=', verify_max_abs, 'max_rel=', verify_max_rel, &
                        & 'max_rel_floor=', verify_max_rel_floor, 'NT=', NT, 'IC=', IC
                else
                   write(6,'(A,1X,A,1X,I6,1X,A,1X,I6,1X,A,1X,I6)') &
                        & '## TX_SOLVER_VERIFY', 'reference_GBSV_failed ierr=', ierr_verify, &
                        & 'NT=', NT, 'IC=', IC
                end if
                deallocate(bl_verify, bx_rhs_verify, bx_ref_verify)
                solver_verify_done = .true.
             end if
#endif
          else if(MDLPCK == 3) then
             call TXSOLV_BAND_EQUIL(BL, size(BL,1), BX, NQMAX, NRMAX, ierr_la, IDIAGL >= 6)
             if(IDIAGL >= 6) call dump_solver_quality(NT,IC,'dgbsvx  ')
             if(ierr_la /= 0) then
                write(6,'(3(A,I6))') '### ERROR(TXLOOP) : DGBSVX, NT = ',  &
                     &              NT, ' -', IC, ' step. IERR=',ierr_la
                IERR = 1
                XN = XP
                goto 180
             end if
          else
             ! +++ NOTE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             !  These LAPACK subroutines are threaded via BLAS in Intel MKL.
             !  However, with regard to the usage of TASK/TX, say, GBTRF internally calls
             !    the LAPACK routine, DGBTF2, which calls the BLAS routine, DGER, which is
             !    NOT multi-threaded.
             !  Thus, these routines are executed as a single threaded regardless of the
             !    MKL link option.
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef TX_SOLVER_BENCH
             if(.not. solver_bench_done) then
                m_bench = NQMAX*(NRMAX+1)
                n_bench = m_bench
                kl_bench = 2*NQMAX-1
                ku_bench = kl_bench
                ldbl_bench = 6*NQMAX-2
                ldbx_bench = NQMAX*(NRMAX+1)
                open(newunit=solver_bench_unit, file=solver_bench_matrix_file, &
                     & status='replace', form='unformatted', access='stream')
                write(solver_bench_unit) NQMAX, NRMAX, m_bench, n_bench, kl_bench, ku_bench, &
                     & ldbl_bench, ldbx_bench, NT, IC
                write(solver_bench_unit) BL
                write(solver_bench_unit) BX
                close(solver_bench_unit)
             end if
#endif
#if   _LAPACK == 2
             m    = NQMAX*(NRMAX+1)
             kl   = 2*NQMAX-1
!           +++ F77 +++
             n    = m
             ku   = kl
             nrhs = 1
             ldBL = 6*NQMAX-2
             ldBX = NQMAX*(NRMAX+1)
             call LAPACK_DGBSV(n,kl,ku,nrhs,BL,ldBL,ipiv,BX,ldBX,ierr_la) 
!             call DGBTRF(m,n,kl,ku,BL,ldBL,ipiv,ierr_la)
!             call DGBTRS('N',n,kl,ku,nrhs,BL,ldBL,ipiv,BX,ldBX,ierr_la)
!#elif _LAPACK == 1
#else
!           +++ F95 +++
!             call LA_GBSV(BL,BX,INFO=ierr_la) ! for self-compiled LAPACK95
             call GBSV(BL,BX,INFO=ierr_la) ! for intel mkl LAPACK95
!             call GBTRF(BL,kl,m,ipiv,ierr_la)
!             call GBTRS(BL,BX,ipiv,kl,'N',ierr_la)
#endif
             if(ierr_la /= 0) then
                write(6,'(3(A,I6))') '### ERROR(TXLOOP) : GBSV, NT = ',  &
                     &              NT, ' -', IC, ' step. IERR=',ierr_la
                IERR = 1
                XN = XP
                goto 180
             end if
#ifdef TX_SOLVER_BENCH
             if(.not. solver_bench_done) then
                open(newunit=solver_bench_unit, file=solver_bench_solution_file, &
                     & status='replace', form='unformatted', access='stream')
                write(solver_bench_unit) NQMAX, NRMAX, m_bench, n_bench, kl_bench, ku_bench, &
                     & ldbl_bench, ldbx_bench, NT, IC
                write(solver_bench_unit) BX
                close(solver_bench_unit)
                solver_bench_done = .true.
                write(6,'(A,1X,A,1X,A,1X,2(A,I6))') &
                     & '[TX_SOLVER_BENCH] dumped', trim(solver_bench_matrix_file), &
                     & trim(solver_bench_solution_file), 'NT=', NT, ' IC=', IC
             end if
#endif
          end if
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_solver_sum = prof_solver_sum + (prof_t1 - prof_t0)
#endif

          ! Copy calculated variables' vector to variable matrix
          forall (NR = 0:NRMAX, NQ = 1:NQMAX) XN(NR,NQ) = BX(NQMAX * NR + NQ)

          ! Avoid negative values
!          call MINUS_CHECK(XN,LQb1,0)
          call MINUS_CHECK(XN,LQb1,3)
          call MINUS_CHECK(XN,LQn1,1)
          call MINUS_CHECK(XN,LQn2,1)
          IF(FSRP /= 0.d0) call MINUS_CHECK(XN,LQr1,2)
          ! Ignore tiny values
          do NQ = 1, NQMAX
             if(any( abs(XN(:,NQ)) < tiny_array(NQ) .and. XN(:,NQ) /= 0.d0 )) then
                nsm_tiny = nsm_tiny + 1
                if(NQ <= size(nsm_tiny_nq)) nsm_tiny_nq(NQ) = nsm_tiny_nq(NQ) + 1
             end if
             where( abs(XN(:,NQ)) < tiny_array(NQ) ) XN(:,NQ) = 0.d0
          end do
!!$          ! In the case of NBI off after NBI on
!!$          IF(PNBH == 0.d0) call THRESHOLD(XN(:,LQb1),ID)
!!$          IF(ID == 1) then
!!$             X (:,LQb1) = 0.d0
!!$             XP(:,LQb1) = 0.d0
!!$             X (:,LQb3) = 0.d0
!!$             XP(:,LQb3) = 0.d0
!!$             XN(:,LQb3) = 0.d0
!!$             X (:,LQb4) = 0.d0
!!$             XP(:,LQb4) = 0.d0
!!$             XN(:,LQb4) = 0.d0
!!$          END IF

          ! Check negative density or temperature in variable matrix
          call TXCHCK(NT,IC,XN,IERR)
          if (IERR /= 0) then
             X = XN
             call TXCALV(X)
!             call TXCALC(IC)
             call TXGLOB
#ifdef TX_PROFILE
             call CPU_TIME(prof_t0)
#endif
             call TX_GRAPH_SAVE
#ifdef TX_PROFILE
             call CPU_TIME(prof_t1)
             prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif
             if(nsub > 1) DT = DTfull
             return
          end if

!          if(IC == 1) EXIT L_IC

          ! Convergence check
          call check_convergence(IC,IDIV,istat)
          if(IDIAGL >= 5) call dump_iter_floor(NT,IC)
#ifdef TX_PROFILE
          call CPU_TIME(prof_ic_t1)
          prof_ic_sum = prof_ic_sum + (prof_ic_t1 - prof_ic_t0)
#endif
          if(istat == 0) then ! going to next time step
             if(IDIAGL >= 8 .and. MDLPCK /= 0) &
                  & call dump_matrix(BL,BX,NQMAX,NRMAX,NT,IC,final=.true.)
             exit L_IC
          else if(istat == 1) then
             ! Size of this sweep's increment, normalised exactly as
             ! check_convergence normalises it.
             errsweep = 0.d0
             do NQ = 1, NQMAX
                avnq = sqrt(sum(XN(:,NQ)**2) / NRMAX)
                if(avnq > epsilon(1.d0)) &
                     & errsweep = max(errsweep, maxval(abs(XN(:,NQ) - XP(:,NQ))) / avnq)
             end do
             if(oldmixa > 0.d0 .and. IC >= DIVIC) then
                if(errsweep > DIVBIG) then
                   ndiv = ndiv + 1
                else
                   ndiv = max(0, ndiv - 1)
                end if
                if(ndiv >= DIVCNT) then
                   ndiv = 0
                   if(.not. damping) then
                      damping = .true.
                      ndamp = ndamp + 1
                      relax = max(relax, oldmixa)
                   else
                      ! still not contracting with the damping in hand; halve
                      ! the distance to one, which is the only knob left
                      relax = min(DIVCAP, 1.d0 - 0.5d0 * (1.d0 - relax))
                   end if
                   if(IDIAGL >= 3) write(6,'(A,I6,A,I4,2(A,ES10.3))') &
                        & ' ## Picard not contracting: NT =', NT, ', IC =', IC, &
                        & ', increment =', errsweep, ', relaxation now', relax
                end if
             end if

             if(relax > 0.d0) XN = (1.d0 - relax) * XN + relax * XP
             ! Legacy nodal Larsen limiter (MDLN0FL = 1) only: under-relax the
             ! neutral densities alone, RLXN0 being the weight of the previous
             ! sweep.  MDLN0FL = 2 does not need it.
             if(MDLN0FL == 1 .and. RLXN0 > 0.d0) &
                  & XN(:,LQn1:LQn3) = (1.d0 - RLXN0) * XN(:,LQn1:LQn3) + RLXN0 * XP(:,LQn1:LQn3)
             if(NANDER /= 0) call anderson_accelerate(IC, damping)
             cycle L_IC
          end if

       end do L_IC

       ! NRETRY = k > 0: a step (or substep) whose Picard iteration reached ICMAX without
       ! converging is taken again from the start of the step as twice as many substeps,
       ! down to DT/2**k; the first substep is backward Euler (no history at the new size),
       ! the others BDF2.  Afterwards the BDF2 history is X(t - DT) again, so the next step
       ! is an ordinary one.  Converged steps are not touched.
       if(IC >= ICMAX .and. nsub < 2**NRETRY) then
          if(nsub == 1) then
             nretried = nretried + 1
             DTfull = DT ; Tstart = T_TX - DT ; Xstart = X
          else
             X = Xstart
          end if
          nsub = 2 * nsub ; isub = 0
          DT = DTfull / nsub
          if(IDIAGL >= 1) write(6,'(A,I6,A,ES10.3)') ' ## Picard not converged at NT =', NT, &
               & ': retried as steps of', DT
          bdf_history_ready = .false.
          if(IGBDF /= 0) XOLD = X
          cycle L_SUB
       end if
       if(nsub > 1 .and. isub < nsub) then   ! a substep is the start of the next one
          if(IGBDF /= 0) then
             XOLD = X ; bdf_history_ready = .true. ; bdf_history_dt = DT
          end if
          X = XN
          call TXCALV(X,1)
          call TXCALC(IC)
          cycle L_SUB
       end if
       exit L_SUB
       end do L_SUB
       if(nsub > 1) then   ! back to DT; the history saved below is X(t - DT)
          DT = DTfull ; X = Xstart
       end if

       ICSUM = ICSUM + IC

       ! The bound n_C0 = n_b of the source-slaved C0 part (MDLNZD = 1) is the
       ! limit of that model (see TXCALC).  Warn when it is still active in the
       ! last sweep, at its onset and when it clears; IDIAG >= 2 prints the
       ! number of bounded nodes at every step in between.
       if(NC0CAP > 0) then
          ntc0cap = ntc0cap + 1
          if(ntc0cap == 1) then
             write(6,'(A,I6,A,I4,A)') ' ## WARNING: n_C0 = n_b in the last Picard sweep at NT =', &
                  & NT, ',', NC0CAP, ' nodes; the source-slaved C0 model is outside its range'
          else if(IDIAGL >= 2) then
             write(6,'(A,I6,A,I4)') ' ## n_C0 = n_b: NT =', NT, ', nodes =', NC0CAP
          end if
       else if(ntc0cap > 0) then
          write(6,'(A,I6,A,I6,A)') ' ## n_C0 = n_b cleared at NT =', NT, ' after', ntc0cap, ' steps'
          ntc0cap = 0
       end if

#ifdef TX_PROFILE
       call CPU_TIME(prof_t0)
#endif
       if(IDIAGL >= 2) then
          if(MODECV == 0) then
             write(6,'(A2,3(A,X),8X,A,15X,A,12X,A,11X,A)') &
                  & "**","IC","VEQ","VNR","XP","XN","V_ERRMAX","EPS"
             if(IC == ICMAX) then ! NOT converged
                write(6,'(3I4,4ES17.8,A2)') IC,IASG(2),IASG(1)-1,XP(IASG(1)-1,IASG(2)), &
                     &                XN(IASG(1)-1,IASG(2)),ASG(IASG(1)-1,IASG(2)),EPSabs," *"
             else ! converged
                IASG(1:2) = maxloc(ASG(:,:))
                write(6,'(3I4,4ES17.8)') IC,IASG(2),IASG(1)-1,XP(IASG(1)-1,IASG(2)), &
                     &                  XN(IASG(1)-1,IASG(2)),ASG(IASG(1)-1,IASG(2)),EPSabs
             end if
          else
             write(6,'(A,2(A,X),5X,A,7X,A)') &
                  & "**","IC","VEQ","V_ERRMAX","EPS"
             write(6,'(2I4,ES17.8,ES10.2)') &
                  & IC,maxloc(L2)-1,maxval(L2),EPSabs
          end if
       end if
#ifdef TX_PROFILE
       call CPU_TIME(prof_t1)
       prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif

#ifdef TX_PROFILE
       call CPU_TIME(prof_t0)
#endif
       if(IDIAG >= 10) then
          if(MODECV == 0) then
             write(6,'(A,2(X,A,X),5X,A,11X,A)') &
                  & "*****","NQ","VR","V_ERRMAX","EPS"
             do NQ = 1, NQMAX
                IASG(1:2)  = maxloc(ASG(:,NQ:NQ))
                write(MSG_NQ,'(4X,2I4,2ES17.8)') NQ,IASG(1)-1,ASG(IASG(1)-1,NQ),EPSabs
                if( ASG(IASG(1)-1,NQ) > EPSabs) then
                   MSG_NQ = trim(MSG_NQ)//' *'
                else
                   MSG_NQ = trim(MSG_NQ)//'  '
                end if
                write(6,'(A80)') MSG_NQ
             end do
          else
             write(6,'(2(A,X),5X,A,7X,A)') &
                  & "*********","NQ","V_ERRMAX","EPS"
             do NQ = 1, NQMAX
                write(MSG_NQ,'(8X,I4,ES17.8,ES10.2)') NQ,L2(NQ),EPSabs
                if(L2(NQ) > EPSabs) then
                   MSG_NQ = trim(MSG_NQ)//' *'
                else
                   MSG_NQ = trim(MSG_NQ)//'  '
                end if
                write(6,'(A80)') MSG_NQ
             end do
          end if
       end if
#ifdef TX_PROFILE
       call CPU_TIME(prof_t1)
       prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif

#ifdef TX_PROFILE
       call CPU_TIME(prof_t0)
#endif
       if(IDIAGL >= 4 .and. MODECV == 0) then
          do nq = 1, NQMAX
             do nr = 0, NRMAX
                write(6,*) nq,nr,ASG(nr,nq)
             end do
          end do
       end if
#ifdef TX_PROFILE
       call CPU_TIME(prof_t1)
       prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif

       if(istat == 2) then
          IERR = 1
          write(6,'(2(A,I6))') '### ERROR(TXLOOP) : Solutions not converged at ',  &
               &              NT, ' -', IC
          goto 180
       end if

       ! Save past X for BDF
       if(IGBDF /= 0) then
          XOLD = X
          bdf_history_ready = .true.
          bdf_history_dt = DT
       end if

       ! Calculation fully converged
       X = XN

       ! Diagnostic only (iaxproj /= 0): put the near-axis profiles back into the
       ! regular subspace.  Off by default; see axis_regularity_projection.
       if(iaxproj /= 0) call axis_regularity_projection

       ! Calculate mesh and coefficients at the next step
       call TXCALV(X,1) ! Set new values to pres0 and ErV0
!!$       PNsV_FIX(:,:) = Var(:,:)%n
!!$       PTsV_FIX(:,:) = Var(:,:)%T
!!$       ErV_FIX (:) = ErV (:)
       call TXCALC(IC)

#ifdef TX_PROFILE
       call CPU_TIME(prof_t0)
#endif
       if(IDIAGL == 0 .or. IDIAGL == 2) then
          if ((mod(NT, NTSTEP) == 0) .and. (NT /= NTMAX)) &
               & write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3)') NT,T_TX,IC
       else if(IDIAGL > 0) then
          if(NT /= NTMAX) then
             if(IC-1 == ICMAX) then
                write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3,"  *")') NT,T_TX,IC
             else
                write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3)') NT,T_TX,IC
             end if
          end if
       end if
#ifdef TX_PROFILE
       call CPU_TIME(prof_t1)
       prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif

180    if(nsub > 1) DT = DTfull
       if (mod(NT, NGRSTP) == 0) then
#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          call TXSTGR(NGR,GT,GY,NGRM)
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif
       end if

       call TXGLOB
       if (mod(NT, NGTSTP) == 0) then
#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          call TXSTGT(real(T_TX))
          call txstgq !!!temporary
          if(IDIAG < 0) call steady_check
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif
       end if

       if (mod(NT, NGVSTP) == 0) then
#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          call TXSTGV(real(T_TX))
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif
       end if

       if (mod(NT, NTMAX ) == 0) then
#ifdef TX_PROFILE
          call CPU_TIME(prof_t0)
#endif
          call TXSTGQ
#ifdef TX_PROFILE
          call CPU_TIME(prof_t1)
          prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif
       end if

!       call cal_flux

       if (IERR /= 0) exit L_NTDO

    end do L_NTDO

    rIPs = rIPe

#ifdef TX_PROFILE
    call CPU_TIME(prof_t0)
#endif
    if(IC == ICMAX) then
       write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3,"  *")') NT,T_TX,IC
    else
       write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3)') NT,T_TX,IC
    end if
#ifdef TX_PROFILE
    call CPU_TIME(prof_t1)
    prof_io_sum = prof_io_sum + (prof_t1 - prof_t0)
#endif
    if(NTMAX /= 0) then
       if(mod(IDIAG,10) >= 5) then
          write(6,'(A,I8,A,I8)') ' ## time steps where adaptive damping engaged:', ndamp, &
               & '   out of', NTMAX
          if(NANDER /= 0) write(6,'(A,I9,A,I6,A,I6,A,ES10.3,A,F9.2,A)') ' ## Anderson: sweeps', and_ncall, &
               & '   restarts', and_nrst, '   LS failures', and_nfail, '   max cond', and_condmax, &
               & '   CPU', and_cpu, ' s'
          write(6,'(A,4I8,A,I8)') ' ## non-smooth guards fired (MINUS_CHECK ID=0,1,2,3):', &
               & nsm_minus(1:4), '   tiny-zeroing (eq,step) pairs:', nsm_tiny
          if(nsm_tiny > 0) then
             write(6,'(A)',advance='no') ' ##   tiny-zeroing by NQ:'
             do NQ = 1, min(NQMAX,size(nsm_tiny_nq))
                if(nsm_tiny_nq(NQ) > 0) write(6,'(A,I0,A,I0)',advance='no') ' ', NQ, ':', nsm_tiny_nq(NQ)
             end do
             write(6,*)
          end if
       end if

       AVE_IC = real(ICSUM) / NTMAX
    else
       AVE_IC = 0.d0
    end if

#ifdef TX_PROFILE
    call CPU_TIME(prof_total_t1)
    write(6,'(A,ES15.7)') '## TX_PROFILE total_sec=',          prof_total_t1 - prof_total_t0
    write(6,'(A,ES15.7)') '## TX_PROFILE iteration_loop_sec=', prof_ic_sum
    write(6,'(A,ES15.7)') '## TX_PROFILE txcala_sec=',         prof_txcala_sum
    write(6,'(A,ES15.7)') '## TX_PROFILE txcalb_sec=',         prof_txcalb_sum
    write(6,'(A,ES15.7)') '## TX_PROFILE solver_sec=',         prof_solver_sum
    write(6,'(A,ES15.7)') '## TX_PROFILE io_sec=',             prof_io_sum
#endif

    deallocate(BA,BL,BX)
    deallocate(XN,XP,ASG,L2)

  end subroutine TXLOOP

!*********************************************************************************
!
!   Calculate coefficients matrix BA or BL and vector BX for linearlized equation
!
!      BA : coefficient matrix for BANDRD solver
!      BL : coefficient matrix for LAPACK DGBSV or LAPACK95 LA_GBSV solver
!      BX : right-hand-side vector
!
!   ** Simple explanation **
!
!      A(u)u=b(u)+c, where A denotes the matrix, u the vector to be solved
!                          b the coefficients vector of u,
!                          c the additive coefficients vector independent on u
!
!*********************************************************************************

  subroutine TXCALB(BA,BL,BX)

    use tx_commons, only : IGBDF, ADV, MDLPCK, NQMAX, NRMAX, NLC, NLCR, CLC, BLC, ALC, NLCMAX, &
         &              PLC, X, XOLD, bdf_history_ready!, lqb1, pnbv
    real(8), dimension(:,:), intent(inout) :: BA, BL
    real(8), dimension(:), intent(inout) :: BX
    INTEGER(4) :: J, NR, NQ, NC, NC1, IA, IB, IC, NCHvs, NC1Hvs, NC1Hvs2, NC2
    INTEGER(4) :: JA, JB, JC, KL
    real(8), parameter :: C43 = 4.d0/3.d0, C23 = 2.d0/3.d0, C13 = 1.d0/3.d0
    real(8) :: COEF1, COEF2, COEF3, COEF!, suml1,suml2
    
    ! The two-step BDF needs a genuine X(t-dt).  Where there is none, XOLD has
    ! been seeded with X, and using the BDF weights there would solve
    ! M u1 = M u0 + ADV dt (L u1 + S) with ADV = 2/3, advancing the solution by
    ! two thirds of a step.  That single O(dt) error costs the whole run its
    ! second order.  Take one backward-Euler step instead: its O(dt^2) local
    ! error, once, leaves the order alone.
    if(IGBDF /= 0) then
       if(bdf_history_ready) then
          ADV = C23
       else
          ADV = 1.d0
       end if
    end if

    !  NR : number of radial mesh 
    !  NQ : number of equation

    ! *** Left-hand-side banded coefficient matrix, denoted by "A" ***

    !***************************************************************
    !   ALC, BLC and CLC are NQMAX x max(NLCMAX) on each NR,
    !   then IC, IB and IA are NQMAX away from each other,
    !   which are pointers of columns of the matrix BA in terms of 
    !   certain variable in certain term in certain equation.
    !***************************************************************

    ! For BANDRD solver
    if(MDLPCK == 0) then
       BA = 0.d0 ! initialize BA array

!!$       DO NQ = 1, NQMAX
!!$          DO NC = 0, NLCMAX(NQ)
!!$             NCHvs = ( NQMAX - 1 + NC ) / NQMAX ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
!!$             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
!!$             NR = 0
!!$                NC1 = NLCR(NC,NQ,0)
!!$                IF(NC1 /= 0) then
!!$                   IC = NQMAX + (NC1 - 1) - (NQ - 1)
!!$                   IB = IC + NQMAX
!!$                   IA = IB + NQMAX
!!$                   J = NR * NQMAX + NQ
!!$                   BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef
!!$                   BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef
!!$                   BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef
!!$                END IF
!!$             NC1 = NLC(NC,NQ)
!!$             IF(NC1 /= 0) then
!!$                IC = NQMAX + (NC1 - 1) - (NQ - 1)
!!$                IB = IC + NQMAX
!!$                IA = IB + NQMAX
!!$                DO NR = 1, NRMAX-1
!!$                   J = NR * NQMAX + NQ
!!$                   BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef
!!$                   BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef
!!$                   BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef
!!$                END DO
!!$             END IF
!!$             NR = NRMAX
!!$                NC1 = NLCR(NC,NQ,1)
!!$                IF(NC1 /= 0) then
!!$                   IC = NQMAX + (NC1 - 1) - (NQ - 1)
!!$                   IB = IC + NQMAX
!!$                   IA = IB + NQMAX
!!$                   J = NR * NQMAX + NQ
!!$                   BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef
!!$                   BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef
!!$                   BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef
!!$                END IF
!!$          END DO
!!$       END DO

       ! *** Introducing NCHvs makes it possible to eliminate NC branch (NC=0 or else).
       ! *** Introducing NC1Hvs and NC1Hvs2 makes it possible to eliminate IF(NC1 /= 0) statement in do loops

       do NQ = 1, NQMAX
          do NC = 0, NLCMAX(NQ)
             NCHvs = ( NQMAX - 1 + NC ) / NQMAX ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
             ! --- NR = 0
             NC1 = NLCR(NC,NQ,0)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC1Hvs2 = 1 - NC1Hvs
             ! NC1Hvs2 in IC: avoid IC=0 in BL when NR = NQMAX and NC1 = 0
             IC = NQMAX + (NC1 - 1) - (NQ - 1) + NC1Hvs2
             IB = IC + NQMAX
             IA = IB + NQMAX
             NR = 0
                J = NR * NQMAX + NQ
                BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef * NC1Hvs
             ! --- NR = 1 ~ NRMAX
             NC1 = NLC(NC,NQ)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC1Hvs2 = 1 - NC1Hvs
             ! NC1Hvs2 in IC: avoid IC=0 in BL when NR = NQMAX and NC1 = 0
             IC = NQMAX + (NC1 - 1) - (NQ - 1) + NC1Hvs2
             IB = IC + NQMAX
             IA = IB + NQMAX
             do concurrent (NR = 1:NRMAX-1)
                J = NR * NQMAX + NQ
                BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef * NC1Hvs
             end do
             ! --- NR = NRMAX
             NC1 = NLCR(NC,NQ,1)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC1Hvs2 = 1 - NC1Hvs
             ! NC1Hvs2 in IC: avoid IC=0 in BL when NR = NQMAX and NC1 = 0
             IC = NQMAX + (NC1 - 1) - (NQ - 1) + NC1Hvs2
             IB = IC + NQMAX
             IA = IB + NQMAX
             NR = NRMAX
                J = NR * NQMAX + NQ
                BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef * NC1Hvs
          end do
       end do

    ! For LAPACK solver
    else
       BL = 0.d0 ! initialize BL array
       KL = 2 * NQMAX - 1

!!$       do NQ = 1, NQMAX
!!$          do NC = 0, NLCMAX(NQ)
!!$             NCHvs = ( NQMAX - 1 + NC ) / NQMAX ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
!!$             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
!!$             NR = 0
!!$                NC1 = NLCR(NC,NQ,0)
!!$                NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX
!!$                IF(NC1 /= 0) then
!!$                   IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
!!$                   IB = IA + NQMAX
!!$                   JA =(NR + 1) * NQMAX + NC1
!!$                   JB = NR      * NQMAX + NC1
!!$                   BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef
!!$                   BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef
!!$                END IF
!!$             NC1 = NLC(NC,NQ)
!!$             if(NC1 /= 0) then
!!$                IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
!!$                IB = IA + NQMAX
!!$                IC = IB + NQMAX
!!$                do NR = 1, NRMAX-1
!!$                   JA =(NR + 1) * NQMAX + NC1
!!$                   JB = NR      * NQMAX + NC1
!!$                   JC =(NR - 1) * NQMAX + NC1
!!$                   BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef
!!$                   BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef
!!$                   BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef
!!$                end do
!!$             end if
!!$             NR = NRMAX
!!$                NC1 = NLCR(NC,NQ,1)
!!$                if(NC1 /= 0) then
!!$                   IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
!!$                   IB = IA + NQMAX
!!$                   IC = IB + NQMAX
!!$                   JB = NR      * NQMAX + NC1
!!$                   JC =(NR - 1) * NQMAX + NC1
!!$                   BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef
!!$                   BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef
!!$                end if
!!$          end do
!!$       end do

       ! *** Introducing NCHvs makes it possible to eliminate NC branch (NC=0 or else).
       ! *** Introducing NC1Hvs and NC1Hvs2 makes it possible to eliminate IF(NC1 /= 0) statement in do loops

       do NQ = 1, NQMAX
          do NC = 0, NLCMAX(NQ)
             NCHvs = ( NQMAX - 1 + NC ) / NQMAX ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
             ! --- NR = 0
             NC1 = NLCR(NC,NQ,0)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC1Hvs2 = 1 - NC1Hvs
             IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
             IB = IA + NQMAX
             NR = 0
                JA =(NR + 1) * NQMAX + NC1
                ! NC1Hvs2 in JB: avoid JB=0 in BL when NC1 = 0
                JB = NR      * NQMAX + NC1 + NC1Hvs2
                BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef * NC1Hvs
                BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef * NC1Hvs
             ! --- NR = 1 ~ NRMAX
             NC1 = NLC(NC,NQ)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX
             NC1Hvs2 = 1 - NC1Hvs
             ! NC1Hvs2 in IA: avoid overflow of IA in BL when NQ=NQMAX and NC1 = 0
             IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL - NC1Hvs2
             IB = IA + NQMAX
             IC = IB + NQMAX
             do concurrent (NR = 1:NRMAX-1)
                JA =(NR + 1) * NQMAX + NC1
                JB = NR      * NQMAX + NC1
                ! NC1Hvs2 in JC: avoid overflow of JC in BL when NR=NRMAX-1 and NC1 = 0
                JC =(NR - 1) * NQMAX + NC1 + NC1Hvs2
                BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef * NC1Hvs
                BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef * NC1Hvs
             end do
             ! --- NR = NRMAX
             NC1 = NLCR(NC,NQ,1)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX
             NC1Hvs2 = 1 - NC1Hvs
             IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
             IB = IA + NQMAX
             ! NC1Hvs2 in IC: avoid overflow of IC in BL when NQ=NQMAX and NC1 = 0
             IC = IB + NQMAX - NC1Hvs2
             NR = NRMAX
                JB = NR      * NQMAX + NC1
                JC =(NR - 1) * NQMAX + NC1
                BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef * NC1Hvs
          end do
       end do

    end if

    ! *** Right-hand-side vector, denoted by "bu" ***

    BX = 0.d0 ! initialize BX array

    if(IGBDF == 0 .or. .not. bdf_history_ready) then
       COEF1 = 1.d0
       COEF2 = 0.d0
       COEF3 = 1.d0
    else
       COEF1 = C43
       COEF2 = C13
       COEF3 = C23
    end if

    ! In the case of NC=0 i.e. time derivative term effects in any equations
!!$    NC = 0
!!$    NR = 0
!!$       do NQ = 1, NQMAX
!!$          NC1 = NLCR(NC,NQ,0)
!!$          if(NC1 /= 0) then
!!$             BX(NQMAX * NR + NQ) &
!!$                  & = BX(NQMAX * NR + NQ) + BLC(NR,NC,NQ) * X   (NR  ,NC1) * COEF1 &
!!$                  &                       + ALC(NR,NC,NQ) * X   (NR+1,NC1) * COEF1 &
!!$                  &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC1) * COEF2 &
!!$                  &                       - ALC(NR,NC,NQ) * XOLD(NR+1,NC1) * COEF2
!!$          end if
!!$       end do
!!$
!!$    do NQ = 1, NQMAX
!!$       NC1 = NLC(NC,NQ)
!!$       do NR = 1, NRMAX - 1
!!$          if(NC1 /= 0) then
!!$             BX(NQMAX * NR + NQ) &
!!$                  &    = BX(NQMAX * NR + NQ) + CLC(NR,NC,NQ) * X   (NR-1,NC1) * COEF1 &
!!$                  &                          + BLC(NR,NC,NQ) * X   (NR  ,NC1) * COEF1 &
!!$                  &                          + ALC(NR,NC,NQ) * X   (NR+1,NC1) * COEF1 &
!!$                  &                          - CLC(NR,NC,NQ) * XOLD(NR-1,NC1) * COEF2 &
!!$                  &                          - BLC(NR,NC,NQ) * XOLD(NR  ,NC1) * COEF2 &
!!$                  &                          - ALC(NR,NC,NQ) * XOLD(NR+1,NC1) * COEF2
!!$          end if
!!$       end do
!!$    end do
!!$
!!$    NR = NRMAX
!!$       do NQ = 1, NQMAX
!!$          NC1 = NLCR(NC,NQ,1)
!!$          if(NC1 /= 0) then
!!$             BX(NQMAX * NR + NQ) &
!!$                  & = BX(NQMAX * NR + NQ) + CLC(NR,NC,NQ) * X   (NR-1,NC1) * COEF1 &
!!$                  &                       + BLC(NR,NC,NQ) * X   (NR  ,NC1) * COEF1 &
!!$                  &                       - CLC(NR,NC,NQ) * XOLD(NR-1,NC1) * COEF2 &
!!$                  &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC1) * COEF2
!!$          end if
!!$       end do

    NC = 0
    do NQ = 1, NQMAX
       ! --- NR = 0
       NC1 = NLCR(NC,NQ,0)
       NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
       NC2 = NC1 + 1 - NC1Hvs
       NR = 0
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) +(BLC(NR,NC,NQ) * X   (NR  ,NC2) * COEF1 &
               &                       + ALC(NR,NC,NQ) * X   (NR+1,NC2) * COEF1 &
               &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC2) * COEF2 &
               &                       - ALC(NR,NC,NQ) * XOLD(NR+1,NC2) * COEF2) * NC1Hvs
!       write(6,'(3I3,4ES15.7)') NQMAX * NR + NQ,NC1,NC2,BLC(NR,NC,NQ) * X   (NR  ,NC2),ALC(NR,NC,NQ) * X   (NR+1,NC2),BLC(NR,NC,NQ) * XOLD(NR  ,NC2),ALC(NR,NC,NQ) * XOLD(NR+1,NC2)

       ! --- NR = 1 ~ NRMAX
       NC1 = NLC(NC,NQ)
       NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
       NC2 = NC1 + 1 - NC1Hvs
       do concurrent (NR = 1:NRMAX - 1)
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) +(CLC(NR,NC,NQ) * X   (NR-1,NC2) * COEF1 &
               &                       + BLC(NR,NC,NQ) * X   (NR  ,NC2) * COEF1 &
               &                       + ALC(NR,NC,NQ) * X   (NR+1,NC2) * COEF1 &
               &                       - CLC(NR,NC,NQ) * XOLD(NR-1,NC2) * COEF2 &
               &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC2) * COEF2 &
               &                       - ALC(NR,NC,NQ) * XOLD(NR+1,NC2) * COEF2) * NC1Hvs
       end do

       ! --- NR = NRMAX
       NC1 = NLCR(NC,NQ,1)
       NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
       NC2 = NC1 + 1 - NC1Hvs
       NR = NRMAX
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) +(CLC(NR,NC,NQ) * X   (NR-1,NC2) * COEF1 &
               &                       + BLC(NR,NC,NQ) * X   (NR  ,NC2) * COEF1 &
               &                       - CLC(NR,NC,NQ) * XOLD(NR-1,NC2) * COEF2 &
               &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC2) * COEF2) * NC1Hvs
    end do

    ! *** In the case of general term effect in any equations, denoted by "c" ***

    forall (NQ = 1:NQMAX, NR = 0:NRMAX) &
         & BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) + sum(PLC(NR,1:NLCMAX(NQ),NQ)) * COEF3

    ! *** Only the case of not BDF ***
    !  Because "ADV" has no longer original meaning when BDF is used.

    if(IGBDF == 0) then
!!$    do NQ = 1, NQMAX
!!$       NR = 0
!!$          do NC = 1, NLCMAX(NQ)
!!$             NC1 = NLCR(NC,NQ,0)
!!$             if(NC1 /= 0) then
!!$                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
!!$                     &              +(  BLC(NR,NC,NQ) * X(NR  ,NC1) &
!!$                     &                + ALC(NR,NC,NQ) * X(NR+1,NC1)) * (1.d0 - ADV)
!!$             end if
!!$          end do
!!$
!!$       do NR = 1, NRMAX-1
!!$          do NC = 1, NLCMAX(NQ)
!!$             NC1 = NLC(NC,NQ)
!!$             if(NC1 /= 0) then
!!$                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
!!$                     &              +(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
!!$                     &                + BLC(NR,NC,NQ) * X(NR  ,NC1) &
!!$                     &                + ALC(NR,NC,NQ) * X(NR+1,NC1)) * (1.d0 - ADV)
!!$             end if
!!$          end do
!!$       end do
!!$
!!$       NR = NRMAX
!!$          do NC = 1, NLCMAX(NQ)
!!$             NC1 = NLCR(NC,NQ,1)
!!$             if(NC1 /= 0) then
!!$                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
!!$                     &              +(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
!!$                     &                + BLC(NR,NC,NQ) * X(NR  ,NC1)) * (1.d0 - ADV)
!!$             end if
!!$          end do
!!$    end do

    do NQ = 1, NQMAX
       NR = 0
          do NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,0)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC2 = NC1 + 1 - NC1Hvs
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                  &              +(  BLC(NR,NC,NQ) * X(NR  ,NC2) &
                  &                + ALC(NR,NC,NQ) * X(NR+1,NC2)) * (1.d0 - ADV) * NC1Hvs
          end do

       do NC = 1, NLCMAX(NQ)
          NC1 = NLC(NC,NQ)
          NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
          NC2 = NC1 + 1 - NC1Hvs
          do concurrent (NR = 1:NRMAX-1)
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                  &              +(  CLC(NR,NC,NQ) * X(NR-1,NC2) &
                  &                + BLC(NR,NC,NQ) * X(NR  ,NC2) &
                  &                + ALC(NR,NC,NQ) * X(NR+1,NC2)) * (1.d0 - ADV) * NC1Hvs
          end do
       end do

       NR = NRMAX
          do NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,1)
             NC1Hvs = ( NQMAX - 1 + NC1 ) / NQMAX ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC2 = NC1 + 1 - NC1Hvs
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                  &              +(  CLC(NR,NC,NQ) * X(NR-1,NC2) &
                  &                + BLC(NR,NC,NQ) * X(NR  ,NC2)) * (1.d0 - ADV) * NC1Hvs
          end do
    end do
    end if

  end subroutine TXCALB

  subroutine TXSOLV_BAND_EQUIL(BL, ldBL, BX, NQMAXL, NRMAXL, IER, detail)

!   Band solve through DGBSVX with FACT = 'E': row and column equilibration, the
!   factorisation, iterative refinement, and estimates of the condition and of
!   the forward and backward error, in one call.
!
!   The forty equations are assembled in their own physical units -- 1e19 m^-3
!   next to O(1) potentials next to 1e-3 flows -- so the row norms of BL span
!   some seventeen orders of magnitude and the columns eight.  Measured on the
!   standard case, that puts the 1-norm condition of the assembled matrix at
!   ~9e22, of which twelve orders are pure scaling: equilibrating it brings the
!   condition to ~2e11, which is what a second-order FEM operator on this mesh
!   must cost and no more.  Equilibration matters beyond the error bound,
!   because partial pivoting picks the largest entry in a column and would
!   otherwise choose its pivot by the row's physical units rather than by the
!   actual dependence.  Nothing here touches the equations or the variables.

    integer(4), intent(in) :: ldBL, NQMAXL, NRMAXL
    real(8), dimension(ldBL,*), intent(inout) :: BL
    real(8), dimension(:), intent(inout) :: BX
    integer(4), intent(out) :: IER
    logical, intent(in), optional :: detail

    integer(4) :: n, kl, ku, ldafb, i, j, i1, i2
    real(8), dimension(:,:), allocatable, save :: AFB, ABsave
    real(8), dimension(:), allocatable, save :: Req, Ceq, Bwk, Xwk, Bsave, work
    integer(4), dimension(:), allocatable, save :: iwork, ipivl
    real(8) :: ferr1(1), berr1(1), aij, num, den, brow
    logical :: want

    n     = NQMAXL * (NRMAXL + 1)
    kl    = 2 * NQMAXL - 1
    ku    = kl
    ldafb = 2 * kl + ku + 1

    if(allocated(AFB)) then
       if(size(AFB,2) /= n .or. size(AFB,1) /= ldafb) then
          deallocate(AFB, Req, Ceq, Bwk, Xwk, work, iwork, ipivl)
          if(allocated(ABsave)) deallocate(ABsave, Bsave)
       end if
    end if
    if(.not. allocated(AFB)) &
         & allocate(AFB(ldafb,n), Req(n), Ceq(n), Bwk(n), Xwk(n), &
         &          work(3*n), iwork(n), ipivl(n))

    want = .false.
    if(present(detail)) want = detail
    ! DGBSVX overwrites BL with the equilibrated band and BX with diag(R)*BX, so
    ! keep both originals when the backward error is to be attributed to a row.
    ! Only the diagnostics need that, so the copies are not kept otherwise.
    if(want) then
       if(.not. allocated(ABsave)) allocate(ABsave(kl+ku+1,n), Bsave(n))
       ABsave(1:kl+ku+1,1:n) = BL(kl+1:kl+ku+1+kl,1:n)
       Bsave(1:n) = BX(1:n)
    end if

    Bwk(1:n) = BX(1:n)

    ! BL holds the band in the DGBSV layout, i.e. the matrix proper starts at row
    ! kl+1; DGBSVX wants it at row 1, so hand it that row with the same leading
    ! dimension.
    call DGBSVX('E', 'N', n, kl, ku, 1, BL(kl+1,1), ldBL, AFB, ldafb, ipivl, &
         &      solq_equed, Req, Ceq, Bwk, n, Xwk, n, solq_rcond, ferr1, berr1, &
         &      work, iwork, IER)

    solq_ferr = ferr1(1)
    solq_berr = berr1(1)

    ! INFO = n+1 means RCOND is below machine precision.  DGBSVX has still
    ! returned the solution together with FERR and BERR saying what it is worth,
    ! so take it and let the caller's diagnostics speak.
    if(IER == n + 1) IER = 0
    if(IER == 0) BX(1:n) = Xwk(1:n)

    ! Which row owns the componentwise backward error DGBSVX reports.  Same
    ! formula as DGBRFS, on the original unequilibrated system:
    !   berr_i = |b - A x|_i / ( (|A| |x|)_i + |b|_i )
    if(want .and. IER == 0) then
       solq_berr_row = 0.d0 ; solq_berr_idx = 0
       do i = 1, n
          num = Bsave(i) ; den = abs(Bsave(i))
          i1 = max(1,i-kl) ; i2 = min(n,i+ku)
          do j = i1, i2
             aij = ABsave(ku+1+i-j,j)
             num = num - aij * Xwk(j)
             den = den + abs(aij) * abs(Xwk(j))
          end do
          if(den > 0.d0) then
             brow = abs(num) / den
             if(brow > solq_berr_row) then
                solq_berr_row = brow ; solq_berr_idx = i
                solq_berr_num = abs(num) ; solq_berr_den = den
             end if
          end if
       end do
    end if

  end subroutine TXSOLV_BAND_EQUIL

  subroutine TXSOLV_BLOCKTRI(BL, BX, NQMAX, NRMAX, IER)

    real(8), dimension(:,:), intent(in) :: BL
    real(8), dimension(:), intent(inout) :: BX
    integer(4), intent(in) :: NQMAX, NRMAX
    integer(4), intent(out) :: IER

    ! Cap on the refinement passes.  Block elimination on its own is *less*
    ! accurate than band LU here -- pivoting is confined to the diagonal blocks,
    ! and with no refinement the iteration floor on LQm2 is 1.3e-3 against the
    ! band solver's 9.6e-6.  It is the refinement that makes this path the better
    ! one, so the loop below checks its own residual rather than taking a fixed
    ! number of blind passes.
    integer(4), parameter :: NREFINE = 3
    integer(4) :: p, nb, n, ldbl, kl, ku, center
    integer(4) :: k, ir, jc, i, j, row, iter, bi, bj, i_min, i_max
    real(8) :: max_far, resnorm, resprev
    real(8), dimension(:,:,:), allocatable :: DORG, DLU, LBLK, UBLK, WBLK
    real(8), dimension(:,:), allocatable :: RHSBLK, SOLBLK
    integer(4), dimension(:,:), allocatable :: IPIV
    real(8), dimension(:), allocatable :: RHS0, XWORK, AX, RES, DX

    IER = 0
    p = NQMAX
    nb = NRMAX + 1
    n = p * nb
    ldbl = size(BL,1)
    kl = 2 * p - 1
    ku = kl
    center = kl + ku + 1

    if(size(BL,2) < n .or. size(BX,1) < n) then
      IER = -1
      return
    end if

    ! Block-tridiagonal solver assumes only main/upper/lower neighboring blocks.
    ! If far blocks are non-zero, request fallback to the original band solver.
    max_far = 0.d0
    do j = 1, n
       bj = (j - 1) / p + 1
       i_min = max(1, j - ku)
       i_max = min(n, j + kl)
       do i = i_min, i_max
          bi = (i - 1) / p + 1
          if(abs(bi - bj) > 1) then
             row = center + i - j
             max_far = max(max_far, abs(BL(row,j)))
          end if
       end do
    end do
    if(max_far /= 0.d0) then
       IER = -2
       return
    end if

    allocate(DORG(p,p,nb), DLU(p,p,nb), LBLK(p,p,nb-1), UBLK(p,p,nb-1), WBLK(p,p,nb-1))
    allocate(RHSBLK(p,nb), SOLBLK(p,nb), IPIV(p,nb))
    allocate(RHS0(n), XWORK(n), AX(n), RES(n), DX(n))

    DORG = 0.d0
    LBLK = 0.d0
    UBLK = 0.d0

    do k = 1, nb
       do jc = 1, p
          j = (k - 1) * p + jc
          do ir = 1, p
             i = (k - 1) * p + ir
             row = center + i - j
             DORG(ir,jc,k) = BL(row,j)
          end do
       end do
    end do

    do k = 1, nb - 1
       do jc = 1, p
          j = k * p + jc
          do ir = 1, p
             i = (k - 1) * p + ir
             row = center + i - j
             UBLK(ir,jc,k) = BL(row,j)
          end do
       end do
    end do

    do k = 2, nb
       do jc = 1, p
          j = (k - 2) * p + jc
          do ir = 1, p
             i = (k - 1) * p + ir
             row = center + i - j
             LBLK(ir,jc,k-1) = BL(row,j)
          end do
       end do
    end do

    do k = 1, nb
       RHSBLK(:,k) = BX((k - 1) * p + 1:k * p)
    end do
    RHS0 = BX(1:n)

    DLU = DORG
    call TXSOLV_BLOCKTRI_FACTOR_SOLVE(DLU, LBLK, UBLK, WBLK, IPIV, RHSBLK, SOLBLK, IER)
    if(IER /= 0) then
       deallocate(DORG, DLU, LBLK, UBLK, WBLK, RHSBLK, SOLBLK, IPIV, RHS0, XWORK, AX, RES, DX)
       return
    end if

    do k = 1, nb
       XWORK((k - 1) * p + 1:k * p) = SOLBLK(:,k)
    end do

    resprev = huge(0.d0)
    do iter = 1, NREFINE
       call TXSOLV_BLOCKTRI_MATVEC(DORG, LBLK, UBLK, XWORK, AX)
       RES = RHS0 - AX
       resnorm = maxval(abs(RES))
       if(resnorm == 0.d0) exit
       ! A pass that no longer halves the residual has nothing left to give.
       if(iter > 1 .and. resnorm > 0.5d0 * resprev) exit
       call TXSOLV_BLOCKTRI_SOLVE_FACT(DLU, LBLK, WBLK, IPIV, RES, DX, IER)
       if(IER /= 0) then
          deallocate(DORG, DLU, LBLK, UBLK, WBLK, RHSBLK, SOLBLK, IPIV, RHS0, XWORK, AX, RES, DX)
          return
       end if
       XWORK = XWORK + DX
       resprev = resnorm
    end do

    BX(1:n) = XWORK

    deallocate(DORG, DLU, LBLK, UBLK, WBLK, RHSBLK, SOLBLK, IPIV, RHS0, XWORK, AX, RES, DX)

  end subroutine TXSOLV_BLOCKTRI

  subroutine TXSOLV_BLOCKTRI_FACTOR_SOLVE(DLU, LBLK, UBLK, WBLK, IPIV, RHSBLK, SOLBLK, IER)

    real(8), dimension(:,:,:), intent(inout) :: DLU
    real(8), dimension(:,:,:), intent(in) :: LBLK, UBLK
    real(8), dimension(:,:,:), intent(out) :: WBLK
    integer(4), dimension(:,:), intent(out) :: IPIV
    real(8), dimension(:,:), intent(inout) :: RHSBLK
    real(8), dimension(:,:), intent(out) :: SOLBLK
    integer(4), intent(out) :: IER

    integer(4) :: p, nb, k, info

    IER = 0
    p = size(DLU,1)
    nb = size(DLU,3)
    WBLK = 0.d0
    SOLBLK = 0.d0

    do k = 1, nb - 1
       call DGETRF(p, p, DLU(:,:,k), p, IPIV(:,k), info)
       if(info /= 0) then
          IER = 200000 + 1000 * k + info
          return
       end if
       WBLK(:,:,k) = UBLK(:,:,k)
       call DGETRS('N', p, p, DLU(:,:,k), p, IPIV(:,k), WBLK(:,:,k), p, info)
       if(info /= 0) then
          IER = 210000 + 1000 * k + info
          return
       end if

       SOLBLK(:,k) = RHSBLK(:,k)
       call DGETRS('N', p, 1, DLU(:,:,k), p, IPIV(:,k), SOLBLK(:,k), p, info)
       if(info /= 0) then
          IER = 220000 + 1000 * k + info
          return
       end if

       RHSBLK(:,k+1) = RHSBLK(:,k+1) - matmul(LBLK(:,:,k), SOLBLK(:,k))
       DLU(:,:,k+1) = DLU(:,:,k+1) - matmul(LBLK(:,:,k), WBLK(:,:,k))
    end do

    call DGETRF(p, p, DLU(:,:,nb), p, IPIV(:,nb), info)
    if(info /= 0) then
       IER = 200000 + 1000 * nb + info
       return
    end if

    SOLBLK(:,nb) = RHSBLK(:,nb)
    call DGETRS('N', p, 1, DLU(:,:,nb), p, IPIV(:,nb), SOLBLK(:,nb), p, info)
    if(info /= 0) then
       IER = 220000 + 1000 * nb + info
       return
    end if

    do k = nb - 1, 1, -1
       SOLBLK(:,k) = SOLBLK(:,k) - matmul(WBLK(:,:,k), SOLBLK(:,k+1))
    end do

  end subroutine TXSOLV_BLOCKTRI_FACTOR_SOLVE

  subroutine TXSOLV_BLOCKTRI_SOLVE_FACT(DLU, LBLK, WBLK, IPIV, RHSVEC, DXVEC, IER)

    real(8), dimension(:,:,:), intent(in) :: DLU, LBLK, WBLK
    integer(4), dimension(:,:), intent(in) :: IPIV
    real(8), dimension(:), intent(in) :: RHSVEC
    real(8), dimension(:), intent(out) :: DXVEC
    integer(4), intent(out) :: IER

    integer(4) :: p, nb, k, info
    real(8), dimension(:,:), allocatable :: G, Y

    IER = 0
    p = size(DLU,1)
    nb = size(DLU,3)
    allocate(G(p,nb), Y(p,nb))

    do k = 1, nb
       G(:,k) = RHSVEC((k - 1) * p + 1:k * p)
    end do

    do k = 1, nb - 1
       Y(:,k) = G(:,k)
       call DGETRS('N', p, 1, DLU(:,:,k), p, IPIV(:,k), Y(:,k), p, info)
       if(info /= 0) then
          IER = 230000 + 1000 * k + info
          deallocate(G, Y)
          return
       end if
       G(:,k+1) = G(:,k+1) - matmul(LBLK(:,:,k), Y(:,k))
    end do

    Y(:,nb) = G(:,nb)
    call DGETRS('N', p, 1, DLU(:,:,nb), p, IPIV(:,nb), Y(:,nb), p, info)
    if(info /= 0) then
       IER = 230000 + 1000 * nb + info
       deallocate(G, Y)
       return
    end if

    do k = nb - 1, 1, -1
       Y(:,k) = Y(:,k) - matmul(WBLK(:,:,k), Y(:,k+1))
    end do

    do k = 1, nb
       DXVEC((k - 1) * p + 1:k * p) = Y(:,k)
    end do

    deallocate(G, Y)

  end subroutine TXSOLV_BLOCKTRI_SOLVE_FACT

  subroutine TXSOLV_BLOCKTRI_MATVEC(DBLK, LBLK, UBLK, XVEC, AXVEC)

    real(8), dimension(:,:,:), intent(in) :: DBLK, LBLK, UBLK
    real(8), dimension(:), intent(in) :: XVEC
    real(8), dimension(:), intent(out) :: AXVEC

    integer(4) :: p, nb, k
    real(8), dimension(:), allocatable :: xk, axk

    p = size(DBLK,1)
    nb = size(DBLK,3)
    AXVEC = 0.d0
    allocate(xk(p), axk(p))

    do k = 1, nb
       xk = XVEC((k - 1) * p + 1:k * p)
       axk = matmul(DBLK(:,:,k), xk)
       if(k < nb) axk = axk + matmul(UBLK(:,:,k), XVEC(k * p + 1:(k + 1) * p))
       if(k > 1)  axk = axk + matmul(LBLK(:,:,k - 1), XVEC((k - 2) * p + 1:(k - 1) * p))
       AXVEC((k - 1) * p + 1:k * p) = axk
    end do

    deallocate(xk, axk)

  end subroutine TXSOLV_BLOCKTRI_MATVEC

!***************************************************************
!
!   Check negative density or temperature
!
!***************************************************************

  subroutine TXCHCK(NTL,IC,XL,IER)

    use tx_commons, only : NQMAX, NRMAX, LQe1, LQi1, LQz1, LQe5, LQi5, LQz5, LQb1, LQr1
    integer(4), intent(in) :: NTL, IC
    integer(4), intent(inout) :: IER
    real(8), dimension(0:NRMAX,1:NQMAX), intent(in) :: XL
    integer(4) :: NR

    IER = 0

    do NR = 0, NRMAX
       if (XL(NR,LQe1) < 0.d0 .or. XL(NR,LQi1) < 0.d0 .or. XL(NR,LQz1) < 0.d0 .or.&
        &  XL(NR,LQb1) < 0.d0 .or. XL(NR,LQr1) < 0.d0) then
          write(6,'(2A,I4,2(A,I4),A)') '### ERROR(TXLOOP) : Negative density at ', &
               &           'NR =', NR, ', NT=', NTL, ', IC=', IC, '.'
          write(6,'(1X,5(A,ES10.3))')  'ne =', real(XL(NR,LQe1)), &
               &                   ',   ni =', real(XL(NR,LQi1)), &
               &                   ',   nz =', real(XL(NR,LQz1)), &
               &                   ',   nb =', real(XL(NR,LQb1)), &
               &                   ', nbrp =', real(XL(NR,LQr1))
          IER = 1
          return
       end if
    end do

    do NR = 0, NRMAX
       if (XL(NR,LQe5) < 0.d0 .or. XL(NR,LQi5) < 0.d0 .or. XL(NR,LQz5) < 0.d0) then
          write(6,'(2A,I4,2(A,I4),A)') '### ERROR(TXLOOP) : Negative temperature at ', &
               &           'NR =', NR, ', NT=', NTL, ', IC=', IC, '.'
          write(6,'(20X,3(A,ES15.7))') 'Te =', real(XL(NR,LQe5)), &
               &                   ',   Ti =', real(XL(NR,LQi5)), &
               &                   ',   Tz =', real(XL(NR,LQz5))
          IER = 2
          return
       end if
    end do

  end subroutine TXCHCK

!***************************************************************
!
!   Negative value forced to be set to zero for densities
!
!***************************************************************

  subroutine MINUS_CHECK(XL,LQ,ID)
    
    use tx_commons, only : NQMAX, NRMAX, NRA
    integer(4), intent(in) :: LQ, ID
    real(8), dimension(0:NRMAX,1:NQMAX), intent(inout) :: XL
    integer(4) :: NR, NZERO, ind

    ind = 0
    if(ID == 3) then
       if(any(XL(:,LQ) < 0.d0)) nsm_minus(4) = nsm_minus(4) + 1
    else
       if(minval(XL(:,LQ)) < 0.d0) nsm_minus(ID+1) = nsm_minus(ID+1) + 1
    end if
    if(ID == 0) then ! Negative values are set to zero inside rho=1.
       if(minval(XL(:,LQ)) < 0.d0) then
          ind = 1
          do NR = 0, NRMAX
             if(XL(NR,LQ) <= 0.d0) then
                if(NR < NRA) then
                   XL(NR,LQ) = 0.d0
                else
                   NZERO = NR
                   exit
                end if
             end if
          end do

          XL(NZERO:NRMAX,LQ) = 0.d0
       end if
    else if(ID == 1) then ! Searching from outboard
       if(minval(XL(:,LQ)) < 0.d0) then
          ind = 1
          do NR = NRMAX, 0, -1
             if(XL(NR,LQ) <= 0.d0) then
                NZERO = NR
                exit
             end if
          end do

          XL(0:NZERO,LQ) = 0.d0
       end if
    else if(ID == 2) then ! Negative values are set to zero.
       if(minval(XL(:,LQ)) < 0.d0) then
          ind = 1
          where(XL(:,LQ) < 0.d0) XL(:,LQ) = 0.d0
       end if
    else ! Taking absolute
       XL(:,LQ) = abs(XL(:,LQ))
    end if
!    if(ind /= 0 .and. (ID == 1 .or. ID == 2)) write(6,'(3(X,A,I3,X))') "MINUS_CHECK : LQ=",LQ,"ID=",ID,"NR=",NZERO

  end subroutine MINUS_CHECK

!!$  subroutine THRESHOLD(XL,ID)
!!$
!!$    real(8), dimension(:), intent(inout) :: XL
!!$    integer(4), intent(out) :: ID
!!$    integer(4) :: NR, NRL
!!$
!!$    ID = 0
!!$    NRL = 0.5 * NRA
!!$    if(minval(XL(0:NRL)) < 1.D-8 .and. maxval(XL(0:NRL)) > 0.d0) then
!!$       XL(:) = 0.d0
!!$       ID = 1
!!$    end if
!!$
!!$  end subroutine THRESHOLD

!***************************************************************
!
!   Checking convergence
!
!***************************************************************

  subroutine axis_regularity_projection

!   Diagnostic only.  A profile that is smooth across the magnetic axis is an
!   even function of rho, so its expansion starts q0 + q2 rho^2 + q4 rho^4.  The
!   equation set does not enforce that anywhere -- BV_1s is the only quantity
!   whose boundedness would, and its singular part is annihilated by the rho^2
!   factors it meets in every other equation -- so a solution with a finite
!   drho slope at the axis is admissible.  This projects the near-axis part of
!   the chosen profiles back onto that even subspace, to test whether axis
!   regularity is what is missing.  It is not a permanent fix and is off unless
!   iaxproj is set.
!
!   The fit is least squares in the FEM volume weight, which makes the constant
!   basis function exact and therefore conserves the particle (or energy)
!   content of the window identically.

    use tx_commons, only : NRMAX, rho, vv, X, iaxproj, rhoaxproj, IDIAG, NT, &
         &                 LQe1, LQi1, LQz1, LQe5, LQi5, LQz5, LQm1

    integer(4) :: nq(7), nn, k, m, NR, NRW, i, j, info, mode
    integer(4) :: ipiv(3)
    real(8) :: A(3,3), rhs(3), w(0:NRMAX), f(3,0:NRMAX), u, chg

    mode = mod(iaxproj-1,3) + 1
    nn = 0
    if(mode >= 1) then
       nq(nn+1) = LQz1 ; nq(nn+2) = LQz5 ; nn = nn + 2
    end if
    if(mode >= 2) then
       nq(nn+1) = LQe1 ; nq(nn+2) = LQe5
       nq(nn+3) = LQi1 ; nq(nn+4) = LQi5 ; nn = nn + 4
    end if
    if(mode >= 3) then
       nq(nn+1) = LQm1 ; nn = nn + 1
    end if

    if(iaxproj >= 4) then
       ! Local form: make the profile linear in V across the first two elements,
       ! which is the same as saying it is even in rho with no rho^4 term there.
       ! Only the axis node moves, and d q / d V at the axis becomes the slope of
       ! the second element, which the outer solution already bounds.  This is the
       ! operation the permanent fix would perform as a constraint on the trial
       ! space, so the diagnostic and the candidate fix are the same statement.
       do k = 1, nn
          u = X(1,nq(k)) - (X(2,nq(k)) - X(1,nq(k))) * (vv(1) - vv(0)) / (vv(2) - vv(1))
          if(mod(IDIAG,10) >= 3) write(6,'(A,I6,A,I4,A,ES11.3)') &
               & ' ## axis relation: NT =', NT, ', NQ =', nq(k), &
               & ', change =', u - X(0,nq(k))
          X(0,nq(k)) = u
       end do
       return
    end if

    NRW = 0
    do NR = 0, NRMAX
       if(rho(NR) <= rhoaxproj) NRW = NR
    end do
    if(NRW < 4) return                      ! too few nodes to fit three modes

    ! lumped FEM weight in V on the window
    w = 0.d0
    do NR = 0, NRW
       if(NR > 0)     w(NR) = w(NR) + 0.5d0*(vv(NR) - vv(NR-1))
       if(NR < NRMAX) w(NR) = w(NR) + 0.5d0*(vv(NR+1) - vv(NR))
    end do

    do NR = 0, NRW
       u = rho(NR)/rhoaxproj
       f(1,NR) = 1.d0
       f(2,NR) = u*u
       f(3,NR) = u*u*u*u
    end do


    do k = 1, nn
       do i = 1, 3
          do j = 1, 3
             A(i,j) = sum(w(0:NRW)*f(i,0:NRW)*f(j,0:NRW))
          end do
          rhs(i) = sum(w(0:NRW)*f(i,0:NRW)*X(0:NRW,nq(k)))
       end do
       call DGESV(3, 1, A, 3, ipiv, rhs, 3, info)
       if(info /= 0) cycle
       chg = 0.d0
       do NR = 0, NRW
          u = rhs(1)*f(1,NR) + rhs(2)*f(2,NR) + rhs(3)*f(3,NR)
          chg = max(chg, abs(u - X(NR,nq(k))))
          X(NR,nq(k)) = u
       end do
       if(mod(IDIAG,10) >= 3) write(6,'(A,I6,2(A,I4),A,ES11.3)') &
            & ' ## axis projection: NT =', NT, ', NQ =', nq(k), &
            & ', nodes =', NRW+1, ', max change =', chg
    end do

  end subroutine axis_regularity_projection

  subroutine dump_matrix(BL,BX,NQMAXL,NRMAXL,NTL,ICL,final)

!   The assembled system, once, for offline analysis of the near-axis blocks.
!   Writes the band LHS exactly as TXCALB built it, the right-hand side, and the
!   elemental blocks term by term, the source vector PLC and the term-to-unknown
!   maps, so the near-axis rows can be rebuilt offline.  Under IDIAG ones-digit
!   >= 8, on the last time step only: the first sweep goes to tx_matrix.bin and
!   the last sweep to tx_matrix_last.bin, the latter so that the nonlinear
!   residual A(XP) XP - b(XP) of the iterate the sweep loop stopped on can be
!   formed offline.

    use tx_commons, only : NTMAX, BLC, ALC, CLC, PLC, NLC, NLCR, NLCMAX, NCM, &
         &                  ADV, IGBDF, DT, bdf_history_ready, X, vv, rho, ICMAX

    real(8), dimension(:,:), intent(in) :: BL
    real(8), dimension(:), intent(in) :: BX
    integer(4), intent(in) :: NQMAXL, NRMAXL, NTL, ICL
    logical, intent(in), optional :: final

    integer(4) :: n, kl, ku, ldab, io
    real(8) :: coef1
    logical :: islast
    character(len=32) :: fname

    if(NTL /= NTMAX) return
    islast = (ICL == ICMAX)
    if(present(final)) islast = islast .or. final
    if(ICL /= 1 .and. .not. islast) return
    fname = 'tx_matrix.bin'
    if(ICL /= 1 .or. islast) fname = 'tx_matrix_last.bin'

    n    = NQMAXL * (NRMAXL + 1)
    kl   = 2 * NQMAXL - 1
    ku   = kl
    ldab = size(BL,1)
    coef1 = 1.d0
    if(IGBDF /= 0 .and. bdf_history_ready) coef1 = 4.d0 / 3.d0

    open(newunit=io, file=trim(fname), status='replace', &
         & form='unformatted', access='stream')
    write(io) NQMAXL, NRMAXL, n, kl, ku, ldab, NTL, ICL, NCM
    write(io) DT, coef1, ADV
    write(io) BL(1:ldab,1:n)
    write(io) BX(1:n)
    write(io) BLC(0:NRMAXL,0:NCM,1:NQMAXL)
    write(io) ALC(0:NRMAXL,0:NCM,1:NQMAXL)
    write(io) CLC(0:NRMAXL,0:NCM,1:NQMAXL)
    write(io) PLC(0:NRMAXL,1:NCM,1:NQMAXL)
    write(io) NLC(0:NCM,1:NQMAXL)
    write(io) NLCR(0:NCM,0:NQMAXL,0:1)
    write(io) NLCMAX(1:NQMAXL)
    write(io) X(0:NRMAXL,1:NQMAXL)
    write(io) vv(0:NRMAXL), rho(0:NRMAXL)
    write(io) XP(0:NRMAXL,1:NQMAXL)
    close(io)
    write(6,'(3A,2(A,I6))') ' ## ', trim(fname), ' written', '  NT =', NTL, '  IC =', ICL

  end subroutine dump_matrix

  subroutine report_conditioning(BL,NQMAXL,NRMAXL,NTL,ICL)
!   Condition of the assembled band matrix, as it stands and after LAPACK's
!   row/column equilibration.  The equation set mixes 1e19 densities with O(1)
!   potentials and O(1e-3) flows in one matrix; this says what that costs.

    real(8), dimension(:,:), intent(in) :: BL
    integer(4), intent(in) :: NQMAXL, NRMAXL, NTL, ICL

    integer(4) :: n, kl, ku, ldab, i, j, i1, i2, info
    real(8) :: anorm, anorme, rcond0, rcond1, rowcnd, colcnd, amax, csum
    real(8), dimension(:,:), allocatable :: AB
    real(8), dimension(:), allocatable :: R, C, work
    integer(4), dimension(:), allocatable :: ipivl, iwork
    character(len=1) :: equed
    integer(4), save :: iuc
    logical, save :: openedc = .false.

    n    = NQMAXL * (NRMAXL + 1)
    kl   = 2 * NQMAXL - 1
    ku   = kl
    ldab = size(BL,1)

    if(.not. openedc) then
       openedc = .true.
       open(newunit=iuc, file='tx_cond.dat', status='replace', action='write')
       write(iuc,'(A)') '#     NT     IC       anorm      rcond_raw       cond_raw' // &
            & '     rcond_equil     cond_equil   rowcnd   colcnd equed'
    end if

    allocate(AB(ldab,n), R(n), C(n), work(3*n), ipivl(n), iwork(n))

    ! --- as assembled
    AB(:,1:n) = BL(:,1:n)
    anorm = 0.d0
    do j = 1, n
       i1 = max(1,j-ku) ; i2 = min(n,j+kl) ; csum = 0.d0
       do i = i1, i2
          csum = csum + abs(AB(kl+ku+1+i-j,j))
       end do
       anorm = max(anorm,csum)
    end do
    call DGBTRF(n,n,kl,ku,AB,ldab,ipivl,info)
    rcond0 = 0.d0
    if(info == 0) call DGBCON('1',n,kl,ku,AB,ldab,ipivl,anorm,rcond0,work,iwork,info)

    ! --- after row and column equilibration
    AB(:,1:n) = BL(:,1:n)
    call DGBEQU(n,n,kl,ku,AB(kl+1,1),ldab,R,C,rowcnd,colcnd,amax,info)
    equed = 'N'
    if(info == 0) then
       call DLAQGB(n,n,kl,ku,AB(kl+1,1),ldab,R,C,rowcnd,colcnd,amax,equed)
       anorme = 0.d0
       do j = 1, n
          i1 = max(1,j-ku) ; i2 = min(n,j+kl) ; csum = 0.d0
          do i = i1, i2
             csum = csum + abs(AB(kl+ku+1+i-j,j))
          end do
          anorme = max(anorme,csum)
       end do
       call DGBTRF(n,n,kl,ku,AB,ldab,ipivl,info)
       rcond1 = 0.d0
       if(info == 0) call DGBCON('1',n,kl,ku,AB,ldab,ipivl,anorme,rcond1,work,iwork,info)
    else
       rcond1 = 0.d0
    end if

    write(iuc,'(2I7,5ES16.6,2ES9.1,3X,A1)') NTL, ICL, anorm, rcond0, 1.d0/max(rcond0,1.d-300), &
         & rcond1, 1.d0/max(rcond1,1.d-300), rowcnd, colcnd, equed
    flush(iuc)

    deallocate(AB,R,C,work,ipivl,iwork)

  end subroutine report_conditioning

  subroutine dump_solver_quality(NTL,ICL,tag)
!   What DGBSVX reports about the solve it just did: the reciprocal condition of
!   the equilibrated matrix, the estimated forward and backward error of this
!   right-hand side, and which scaling it applied.

    use tx_commons, only : nqmax_for_dump => NQMAX

    integer(4), intent(in) :: NTL, ICL
    character(len=*), intent(in) :: tag
    integer(4), save :: iuq
    logical, save :: openedq = .false.

    if(.not. openedq) then
       openedq = .true.
       open(newunit=iuq, file='tx_solver.dat', status='replace', action='write')
       write(iuq,'(A)') '#     NT     IC  path           rcond            cond' // &
            & '            ferr            berr  equed   berr_row      NQ    NR' // &
            & '     |resid|        |A||x|+|b|'
    end if

    write(iuq,'(2I7,2X,A8,4ES16.6,5X,A1,ES13.4,2I6,2ES16.6)') NTL, ICL, tag, solq_rcond, &
         & 1.d0/max(solq_rcond,1.d-300), solq_ferr, solq_berr, solq_equed, &
         & solq_berr_row, mod(solq_berr_idx-1,nqmax_for_dump)+1, (solq_berr_idx-1)/nqmax_for_dump, &
         & solq_berr_num, solq_berr_den
    flush(iuq)

  end subroutine dump_solver_quality

  subroutine dump_iter_floor(NTL,ICL)
!   Per-sweep record of the Picard increment, normalised exactly as ASG is in
!   check_convergence, for every equation at once.  check_convergence returns at
!   the first node it finds unconverged, so its printout only ever shows the
!   lowest-numbered offender; this writes the whole picture to tx_iter.dat.

    use tx_commons, only : NRMAX, NQMAX, NEMAX, hv, vv

    integer(4), intent(in) :: NTL, ICL
    integer(4) :: NQ, NRL(1:NQMAX)
    real(8) :: AV, AMX(1:NQMAX), AVS(1:NQMAX)
    integer(4), save :: iunit
    logical, save :: opened = .false.

    if(.not. opened) then
       opened = .true.
       open(newunit=iunit, file='tx_iter.dat', status='replace', action='write')
       write(iunit,'(A,I5,3(A,ES12.4))') '# NRMAX=', NRMAX, '  hv(1)=', hv(1), &
            & '  hv(NEMAX)=', hv(NEMAX), '  vv(NRMAX)=', vv(NRMAX)
       write(iunit,'(A)') '# NT IC  then per NQ: NR_of_max  ASG_max  AV(profile RMS)'
    end if

    do NQ = 1, NQMAX
       AV = sqrt(sum(XN(:,NQ)**2) / NRMAX)
       AVS(NQ) = AV
       if(AV < epsilon(1.d0)) then
          NRL(NQ) = 0 ; AMX(NQ) = 0.d0
       else
          NRL(NQ) = sum(maxloc(abs(XN(:,NQ) - XP(:,NQ)))) - 1
          AMX(NQ) = abs(XN(NRL(NQ),NQ) - XP(NRL(NQ),NQ)) / AV
       end if
    end do

    write(iunit,'(2I6,40(I4,2ES12.4))') NTL, ICL, (NRL(NQ), AMX(NQ), AVS(NQ), NQ = 1, NQMAX)
    flush(iunit)

  end subroutine dump_iter_floor

  subroutine anderson_accelerate(IC, noncontracting)

!   Anderson acceleration (type II, mixing beta = 1, depth m = |NANDER|) of the
!   Picard sweeps.  On entry XP is the iterate x_k of this sweep and XN its
!   Picard image g_k = G(x_k), after any relaxation; on exit XN is
!
!     x_{k+1} = g_k - sum_j gamma_j ( dx_j + df_j ),
!     gamma = argmin | f_k - sum_j gamma_j df_j |_2,   f = g - x,
!     dx_j = x_j - x_{j-1},  df_j = f_j - f_{j-1}   (the last m sweeps).
!
!   A fixed point of G is a fixed point of the accelerated iteration, so the
!   converged solution is unchanged.  Each variable is scaled by its profile
!   RMS at the first sweep of the step, as check_convergence scales its
!   increment.  With NANDER < 0 the rows of LQm1-LQm5 do not enter the fit
!   (near convergence their increments are rounding; for psi_t dot, the
!   floor described in check_convergence), but the correction is applied to
!   all variables.  The history starts afresh at every
!   time step, and is dropped whenever the least-squares solve fails or returns
!   a combination larger than GAMCAP, in which case g_k is taken as it is.
!   Once the adaptive relaxation has found the step not contracting
!   (noncontracting = .true.), the rest of the step is left to the relaxation
!   alone: the linear model behind the mixing does not hold there, and the
!   extrapolation drives the iterate away, as in the first steps from the
!   initial profiles at NRMAX = 960.

    use tx_commons, only : NQMAX, NRMAX, IDIAG, NT, NANDER
    integer(4), intent(in) :: IC
    logical, intent(in) :: noncontracting
    real(8), parameter :: GAMCAP = 1.d3, RCOND = 1.d-12
    integer(4) :: N, mdep, NQ, NR, j, m, info, rank, lwork, liwork
    real(8) :: wq(1), t0, t1
    integer(4) :: iwq(1)

    call CPU_TIME(t0)
    N = NQMAX * (NRMAX + 1) ; mdep = abs(NANDER)
    if(allocated(and_x0)) then
       if(size(and_x0) /= N .or. size(and_dx,2) /= mdep) then
          deallocate(and_dx, and_df, and_A, and_x0, and_f0, and_x, and_g, and_b, and_s, &
               &     and_work, and_iwork, and_sc)
       end if
    end if
    if(.not. allocated(and_x0)) then
       allocate(and_dx(N,mdep), and_df(N,mdep), and_A(N,mdep), and_x0(N), and_f0(N), and_x(N), &
            &   and_g(N), and_b(N), and_s(mdep), and_sc(NQMAX))
       call DGELSD(N, mdep, 1, and_A, N, and_b, N, and_s, RCOND, rank, wq, -1, iwq, info)
       lwork = int(wq(1)) ; liwork = max(1, iwq(1))
       allocate(and_work(lwork), and_iwork(liwork))
    end if
    lwork = size(and_work) ; liwork = size(and_iwork)

    if(IC == 1) then
       and_k = 0 ; and_have0 = .false.
       do NQ = 1, NQMAX
          and_sc(NQ) = sqrt(sum(XP(:,NQ)**2) / NRMAX)
          if(and_sc(NQ) <= 0.d0) and_sc(NQ) = 1.d0
       end do
    end if
    if(noncontracting) then
       and_k = 0 ; and_have0 = .false.
       call CPU_TIME(t1)
       and_cpu = and_cpu + (t1 - t0)
       return
    end if

    do NR = 0, NRMAX
       do NQ = 1, NQMAX
          j = NQMAX*NR + NQ
          and_x(j) = XP(NR,NQ) / and_sc(NQ)
          and_g(j) = XN(NR,NQ) / and_sc(NQ)
       end do
    end do

    ! update the history with this sweep: f_k - f_{k-1}, x_k - x_{k-1}
    if(and_have0) then
       if(and_k == mdep) then
          and_dx(:,1:mdep-1) = and_dx(:,2:mdep)
          and_df(:,1:mdep-1) = and_df(:,2:mdep)
          and_k = mdep - 1
       end if
       and_k = and_k + 1
       and_dx(:,and_k) = and_x - and_x0
       and_df(:,and_k) = (and_g - and_x) - and_f0
    end if
    and_x0 = and_x ; and_f0 = and_g - and_x ; and_have0 = .true.

    if(and_k > 0) then
       m = and_k
       and_ncall = and_ncall + 1
       and_A(:,1:m) = and_df(:,1:m) ; and_b = and_f0
       if(NANDER < 0) then
          do NR = 0, NRMAX
             and_A(NQMAX*NR+1:NQMAX*NR+5,1:m) = 0.d0
             and_b(NQMAX*NR+1:NQMAX*NR+5)     = 0.d0
          end do
       end if
       call DGELSD(N, m, 1, and_A, N, and_b, N, and_s, RCOND, rank, and_work, lwork, and_iwork, info)
       if(info == 0 .and. maxval(abs(and_b(1:m))) < GAMCAP) then
          if(and_s(max(rank,1)) > 0.d0) and_condmax = max(and_condmax, and_s(1) / and_s(max(rank,1)))
          do NR = 0, NRMAX
             do NQ = 1, NQMAX
                j = NQMAX*NR + NQ
                XN(NR,NQ) = ( and_g(j) - sum( and_b(1:m) * (and_dx(j,1:m) + and_df(j,1:m)) ) ) &
                     &    * and_sc(NQ)
             end do
          end do
       else
          if(info /= 0) then
             and_nfail = and_nfail + 1
          else
             and_nrst = and_nrst + 1
          end if
          if(mod(IDIAG,10) >= 3) write(6,'(A,I6,A,I4,A,I4,A,ES10.3)') &
               & ' ## Anderson history dropped: NT =', NT, ', IC =', IC, ', info =', info, &
               & ', max|gamma| =', maxval(abs(and_b(1:m)))
          and_k = 0
       end if
    end if
    call CPU_TIME(t1)
    and_cpu = and_cpu + (t1 - t0)

  end subroutine anderson_accelerate

  subroutine check_convergence(IC,IDIV,istat)
    use tx_commons, only : MODECV, NRMAX, NQMAX, EPS, IDIAG, NT, ICMAX, &
         &                 LQm2, LQm5, ADV, DT, rMU0

    integer(4), intent(in)  :: IC, IDIV
    integer(4), intent(out) :: istat
    integer(4) :: NR, NQ, NRAVM, IDISP, IDIAGL
    real(8) :: EPSabs, AV, AVM, ERR1, FLOORm2
    ! The increment of psi_t dot (LQm2) has a rounding floor that no iteration
    ! removes.  psi_t is mostly vacuum flux, so in Ampere's law (the LQm2 row)
    ! the rot B_phi term is a cancellation of terms larger than itself by a
    ! factor ~ NRMAX**2, and psi_t is determined only to a relative
    ! ~ 0.05-0.07 eps_mach NRMAX**2.  Faraday's law (the LQm5 row) then gives
    ! psi_t dot = mu0 (X5 - X5_old) / (ADV DT) with X5 = psi_t / mu0, a
    ! difference quotient that multiplies that rounding by
    ! psi_t / (ADV DT |psi_t dot|), which grows without bound as psi_t
    ! approaches a steady state.  LQm2 stays in the test: its increment is
    ! tested against EPS * RMS(X2), as every variable is, plus this floor,
    ! CRNDm2 eps_mach NRMAX**2 mu0 RMS(X5) / (ADV DT).  The test of every
    ! other variable, and the fixed point, are unchanged.  CRNDm2 is not a
    ! physical tolerance but a safety factor over the floor: 0.3 lies on the
    ! largest floor observed and leaves no margin for rounding that depends on
    ! the compiler, the BLAS or the order of operations; 1 is three to four
    ! times it.
    real(8), parameter :: CRNDm2 = 1.0d0

    IDIAGL = mod(IDIAG,10)
    EPSabs = abs(EPS)
    FLOORm2 = CRNDm2 * epsilon(1.d0) * dble(NRMAX)**2 &
         &  * rMU0 * sqrt(SUM(XN(:,LQm5)**2) / NRMAX) / (ADV * DT)

    if(MODECV == 0) then

       ASG = 0.d0
       L_NQ:do NQ = 1, NQMAX
          ! Calculate maximum local root-mean-square of X and corresponding grid point
          NRAVM = sum(maxloc(abs(XN(:,NQ))))-1
          AVM   = abs(XN(NRAVM,NQ))

          ! Calculate root-mean-square X over the profile (Euclidean norm)
          AV    = sqrt(SUM(XN(:,NQ)**2) / NRMAX)
          if (AV < epsilon(1.d0)) then
             ! It means that the variable for NQ is zero over the profile.
             ASG(:,NQ) = 0.d0
             cycle L_NQ
          end if

          ERR1  = 0.d0
          IDISP = IDIV
          if(NQ == LQm2) then
             ! the part of the increment above the rounding floor
             ASG(:,NQ) = max(0.d0, abs(XN(:,NQ) - XP(:,NQ)) - FLOORm2) / AV
          else
             ASG(:,NQ) = abs(XN(:,NQ) - XP(:,NQ)) / AV
          end if

          L_NR:do NR = 0, NRMAX
             ! Show results
             if (NT == IDISP .and. NR == NRMAX) then
                if (NQ == 1) then
                   write(6,'(1X,A5,A3," =",I3)')'#####','IC',IC
                   write(6,'(1X,A7," =",ES9.2,2X,A7," =",I3)') &
                        & 'EPS    ', EPSabs,'NRMAX  ', NRMAX
                end if
                ! Maximum relative error over the profile
                ERR1 = max(ERR1, ASG(NR,NQ))
                write(6,'((1X,A2," =",I2,2X,A2," =",ES9.2, &
                     & 2X,A5," =",ES9.2,A1,I2,2X,A5," =",ES9.2))') &
                     & 'NQ    ', NQ , &
                     & 'AV    ', AV    ,  'AVMAX ', AVM   ,':',NRAVM, &
                     & 'SCMAX ', ERR1
                IDISP = IDIV + NT
                if (NQ == NQMAX) then ! not converged
                   istat = 1
                   return
                end IF
             else if (NT /= IDISP .and. ASG(NR,NQ) > EPSabs) then
                ! NOT converged
                if(IC /= ICMAX) then
                   if(IDIAGL >= 3) then
                      if(mod(IC,20) == 0 .or. IC == 1) &
                           & write(6,'(A2,3(A,X),8X,A,15X,A,12X,A,11X,A)') &
                           & "**","IC","VEQ","VNR","XP","XN","V_ERRMAX","EPS"
                      write(6,'(3I4,4ES17.8)') IC,NQ,NR,XP(NR,NQ),XN(NR,NQ), &
                           & ASG(NR,NQ),EPSabs
                   end if
                   istat = 1
                   return
                end if
             end if
          end do L_NR

       end do L_NQ
       ! Converged or IC == ICMAX
       IASG(1:2) = maxloc(ASG)

       if(eps > 0.d0) then ! calculation continues even if not converged
          istat = 0
       else ! calculation stops if not converged
          if(istat == 1) then
             istat = 2 ! not converged at IC = ICMAX
          else
             istat = 0 ! converged
          end if
       end if
!       if(EPS > 0.d0) EXIT L_IC

    else
       L_NQ2:do NQ = 1, NQMAX
          AV = abs(sum(XN(:,NQ)))
          if(AV /= 0.d0) then
             ASG(:,NQ) = (XN(:,NQ) - XP(:,NQ))**2
             L2(NQ) = sqrt(sum(ASG(:,NQ))) / AV
          else
             L2(NQ) = 0.d0
          end if
       end do L_NQ2

       ! Converged
       if(maxval(L2) < EPSabs) then
!          if(EPS > 0.d0) EXIT L_IC
          if(eps > 0.d0) then ! calculation continues even if not converged
             istat = 0
             return
          end if
       end if

       if(IC /= ICMAX) then
          if(IDIAGL >= 3) then
             if(mod(IC,20) == 0 .or. IC == 1) &
                  & write(6,'(A,2(A,X),5X,A,7X,A)') &
                  & "**","IC","VEQ","V_ERRMAX","EPS(1)"
             write(6,'(2I4,ES17.8,ES10.2)') &
                  & IC,maxloc(L2)-1,maxval(L2),EPSabs
          end if
       else
          istat = 2
!          EXIT L_IC ! This is used so that IC exceeds ICMAX after "END DO L_IC".
       end if
    end if

  end subroutine CHECK_CONVERGENCE

end module tx_main
