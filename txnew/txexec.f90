module tx_main
  implicit none
  private
  integer(4), dimension(1:2) :: iasg
  real(8), dimension(:), allocatable :: L2
  real(8), dimension(:,:), allocatable :: XN, XP, ASG
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
         &                 tiny_cap, EPS, IDIAG, NTSTEP, &
         &                 FSRP, MODEAV, XOLD, &
         &                 NT, DT, rIP, MDLPCK, ICONT, AVE_IC, MODECV, oldmix, Var
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

    ! Save X -> XP -> XOLD for BDF only at the beginning of the calculation
    if(IGBDF /= 0 .and. (T_TX == 0.d0 .or. ICONT /= 0)) XOLD= X

    L_NTDO:do NTDO = 1, NTMAX
       NT = NTDO
       NTCUM = NTCUM + 1
       T_TX = TIME0 + DT*NT
       rIP  = rIPs  + DIP*NT

       ! Create new X := XN
       XN = X

       ! Negligible order of magnitude for each variable
       do NQ = 1, NQMAX
          tiny_array(NQ) = maxval(XN(:,NQ)) * tiny_cap
       end do

       ! In the following loop, XN is being updated during iteration.
       L_IC : do IC = 1, ICMAX
          ! Save past X := XP
          XP = XN

          call TXCALV(XP)
          if(IC <= 2) then
             PNsV_FIX(:,:) = Var(:,:)%n
             PTsV_FIX(:,:) = Var(:,:)%T
             ErV_FIX (:) = ErV (:)
          end if

          call TXCALC(IC)
          call TXCALA
          ! Get BA or BL, and BX
          call TXCALB(BA,BL,BX)
!          call TXGLOB

          if(MDLPCK == 0) then
             call BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQMAX-1, IERR)
             if (IERR >= 30000) then
                write(6,'(3(A,I6))') '### ERROR(TXLOOP) : Matrix BA is singular at ',  &
                     &              NT, ' -', IC, ' step. IERR=',IERR
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
          end if

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
             call TX_GRAPH_SAVE
             return
          end if

!          if(IC == 1) EXIT L_IC

          ! Convergence check
          call check_convergence(IC,IDIV,istat)
          if(istat == 0) then ! going to next time step
             exit L_IC
          else if(istat == 1) then
              XN = (1.d0 - oldmix) * XN + oldmix * XP
             cycle L_IC
          end if

       end do L_IC

       ICSUM = ICSUM + IC

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
          
       if(IDIAGL >= 4 .and. MODECV == 0) then
          do nq = 1, NQMAX
             do nr = 0, NRMAX
                write(6,*) nq,nr,ASG(nr,nq)
             end do
          end do
       end if

       if(istat == 2) then
          IERR = 1
          write(6,'(2(A,I6))') '### ERROR(TXLOOP) : Solutions not converged at ',  &
               &              NT, ' -', IC
          goto 180
       end if

       ! Save past X for BDF
       if(IGBDF /= 0) XOLD = X

       ! Calculation fully converged
       X = XN

       ! Calculate mesh and coefficients at the next step
       call TXCALV(X,1) ! Set new values to pres0 and ErV0
!!$       PNsV_FIX(:,:) = Var(:,:)%n
!!$       PTsV_FIX(:,:) = Var(:,:)%T
!!$       ErV_FIX (:) = ErV (:)
       call TXCALC(IC)

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

180    if (mod(NT, NGRSTP) == 0) call TXSTGR(NGR,GT,GY,NGRM)

       call TXGLOB
       if (mod(NT, NGTSTP) == 0) then
          call TXSTGT(real(T_TX))
          call txstgq !!!temporary
          if(IDIAG < 0) call steady_check
       end if

       if (mod(NT, NGVSTP) == 0) call TXSTGV(real(T_TX))

       if (mod(NT, NTMAX ) == 0) call TXSTGQ

!       call cal_flux

       if (IERR /= 0) exit L_NTDO

    end do L_NTDO

    rIPs = rIPe

    if(IC == ICMAX) then
       write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3,"  *")') NT,T_TX,IC
    else
       write(6,'(1x,"NT =",I4,"   T =",ES9.2,"   IC =",I3)') NT,T_TX,IC
    end if
    if(NTMAX /= 0) then
       AVE_IC = real(ICSUM) / NTMAX
    else
       AVE_IC = 0.d0
    end if

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
         &              PLC, X, XOLD!, lqb1, pnbv
    real(8), dimension(:,:), intent(inout) :: BA, BL
    real(8), dimension(:), intent(inout) :: BX
    INTEGER(4) :: J, NR, NQ, NC, NC1, IA, IB, IC, NCHvs, NC1Hvs, NC1Hvs2, NC2
    INTEGER(4) :: JA, JB, JC, KL
    real(8), parameter :: C43 = 4.d0/3.d0, C23 = 2.d0/3.d0, C13 = 1.d0/3.d0
    real(8) :: COEF1, COEF2, COEF3, COEF!, suml1,suml2
    
    if(IGBDF /= 0) ADV = C23

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
             do NR = 1, NRMAX-1
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
             do NR = 1, NRMAX-1
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

    if(IGBDF == 0) then
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
       do NR = 1, NRMAX - 1
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
          do NR = 1, NRMAX-1
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

  subroutine check_convergence(IC,IDIV,istat)
    use tx_commons, only : MODECV, NRMAX, NQMAX, EPS, IDIAG, NT, ICMAX

    integer(4), intent(in)  :: IC, IDIV
    integer(4), intent(out) :: istat
    integer(4) :: NR, NQ, NRAVM, IDISP, IDIAGL
    real(8) :: EPSabs, AV, AVM, ERR1

    IDIAGL = mod(IDIAG,10)
    EPSabs = abs(EPS)

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
          ASG(:,NQ) = abs(XN(:,NQ) - XP(:,NQ)) / AV

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
