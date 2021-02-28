!     $Id: txexec.f90,v 1.55 2011/05/02 06:45:52 fukuyama Exp $
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

  SUBROUTINE TXEXEC
    use tx_commons, only : IERR, T_TX, AVE_IC
    use tx_interface, only : APTOS

    INTEGER(4) :: NSTR1, NSTR2, NSTR3, NSTR4
    REAL(4) :: gCTIME1, gCTIME2, gCTIME3
    character(len=10) :: STR1, STR2, STR3, STR4
    real(8) :: t_interval

    IF (IERR /= 0) THEN
       WRITE(6,*) '### ERROR(TXEXEC) : Error should be cleared.'
       RETURN
    END IF

    CALL CPU_TIME(gCTIME1)
    t_interval = T_TX

    ! ***** Core part of simulation *****

    WRITE(6,*) 'Calculating...'
    CALL TXLOOP

    ! ***********************************

    CALL CPU_TIME(gCTIME2)
    gCTIME3 = gCTIME2-gCTIME1
    t_interval = T_TX - t_interval

    NSTR1 = 0
    CALL APTOS(STR1, NSTR1, gCTIME3, 'F2')
    NSTR2 = 0
    CALL APTOS(STR2, NSTR2, REAL(t_interval), 'E1')
    IF(gCTIME3 /= 0) THEN
       NSTR3 = 0
       CALL APTOS(STR3, NSTR3, REAL(t_interval) / gCTIME3 * 100, 'F3')
       NSTR4 = 0
       CALL APTOS(STR4, NSTR4, REAL(AVE_IC),'F2')
       WRITE(6,*) 'CPU = ', STR1(1:NSTR1), ' (sec)   ',  &
            &     'sim time = ', STR2(1:NSTR2), ' (sec)', '  (', STR3(1:NSTR3), '%)  ', &
            &     'ICave = ', STR4(1:NSTR4)
    ELSE
       WRITE(6,*) 'CPU = ', STR1(1:NSTR1), ' (sec)   ',  &
            &     'sim time = ', STR2(1:NSTR2), ' (sec)'
    END IF

    !     ***** Print simulation results *****

    CALL TXWDAT

    RETURN
  END SUBROUTINE TXEXEC

!***************************************************************
!
!   Core routine of simulation
!
!***************************************************************

  SUBROUTINE TXLOOP
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
    USE libbnd
#ifdef lapack95
#ifdef laself
    ! for self-compiled lapack
    use f95_lapack, only : GBSV => LA_GBSV
#else
    ! for intel mkl LAPACK95, 
    !  Note: This module file includes "ptsv" subroutine, 
    !        whose name conflicts with PTsV defined in TASK/TX.  
    use lapack95, only : GBSV
#endif
#endif
    real(8), dimension(:,:), allocatable :: BA, BL
    real(8), dimension(:),   allocatable :: BX
    integer(4) :: NR, NQ, IC = 0, IDIV, NTDO, ICSUM, IDIAGL, istat
    ! *** LAPACK ************************************
    integer(4) :: ierr_la
!    integer(4) :: m, kl, !, n, ku, nrhs, ldBL, ldBX
!    integer(4), dimension(1:NQMAX*(NRMAX+1)) :: ipiv
    ! ***********************************************
    real(8) :: TIME0, DIP, EPSabs
    real(8), dimension(1:NQMAX) :: tiny_array
    character(len=80) :: MSG_NQ
    integer(4):: m,kl,n,ku,nrhs,ldbl,ldbx
    integer(4),DIMENSION(:),ALLOCATABLE:: ipiv

    allocate( BA(1:4*NQMAX-1,1:NQMAX*(NRMAX+1)) &
         &  , BL(1:6*NQMAX-2,1:NQMAX*(NRMAX+1)) &
         &  , BX(1:NQMAX*(NRMAX+1)))
    allocate( XN (0:NRMAX,1:NQMAX) &
         &  , XP (0:NRMAX,1:NQMAX) &
         &  , ASG(0:NRMAX,1:NQMAX), L2(1:NQMAX))

    IDIAGL = MOD(IDIAG,10)
    EPSabs = abs(EPS)

    !  Read spline table for neoclassical toroidal viscosity if not loaded when FSRP/=0
    IF(FSRP /= 0.D0) CALL Wnm_spline

    IF (MODEAV == 0) THEN
       IDIV = NTMAX + 1
    ELSE
       IDIV = NTMAX / MODEAV
    END IF
    TIME0 = T_TX
    IF(NTMAX /= 0) DIP = (rIPe - rIPs) / NTMAX
    ICSUM = 0

    ! Save X -> XP -> XOLD for BDF only at the beginning of the calculation
    IF(IGBDF /= 0 .and. (T_TX == 0.D0 .OR. ICONT /= 0)) XOLD= X

    L_NTDO:DO NTDO = 1, NTMAX
       NT = NTDO
       NTCUM = NTCUM + 1
       T_TX = TIME0 + DT*NT
       rIP  = rIPs  + DIP*NT

       ! Create new X := XN
       XN = X

       ! Negligible order of magnitude for each variable
       DO NQ = 1, NQMAX
          tiny_array(NQ) = maxval(XN(:,NQ)) * tiny_cap
       END DO

       ! In the following loop, XN is being updated during iteration.
       L_IC : DO IC = 1, ICMAX
          ! Save past X := XP
          XP = XN

          CALL TXCALV(XP)
          IF(IC <= 2) THEN
             PNsV_FIX(:,1:2) = Var(:,1:2)%n
             PTsV_FIX(:,1:2) = Var(:,1:2)%T
             ErV_FIX (:) = ErV (:)
          END IF

          CALL TXCALC(IC)
          CALL TXCALA
          ! Get BA or BL, and BX
          CALL TXCALB(BA,BL,BX)
!          CALL TXGLOB

          IF(MDLPCK == 0) THEN
             CALL BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQMAX-1, IERR)
             IF (IERR >= 30000) THEN
                WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : Matrix BA is singular at ',  &
                     &              NT, ' -', IC, ' step. IERR=',IERR
                IERR = 1
                XN = XP
                GOTO 180
             END IF
          ELSE IF(MDLPCK == 1) THEN
! +++ NOTE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  These LAPACK subroutines are threaded via BLAS in Intel MKL.
!  However, with regard to the usage of TASK/TX, say, GBTRF internally calls
!    the LAPACK routine, DGBTF2, which calls the BLAS routine, DGER, which is
!    NOT multi-threaded.
!  Thus, these routines are executed as a single threaded regardless of the
!    MKL link option.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             m    = NQMAX*(NRMAX+1)
             kl   = 2*NQMAX-1
!!           +++ F77 +++
             n    = m
             ku   = kl
             nrhs = 1
             ldBL = 6*NQMAX-2
             ldBX = NQMAX*(NRMAX+1)
             ALLOCATE(ipiv(n))
             CALL LAPACK_DGBSV(n,kl,ku,nrhs,BL,ldBL,ipiv,BX,ldBX,ierr_la) 
             CALL DGBTRF(m,n,kl,ku,BL,ldBL,ipiv,ierr_la)
             CALL DGBTRS('N',n,kl,ku,nrhs,BL,ldBL,ipiv,BX,ldBX,ierr_la)
             DEALLOCATE(ipiv)
#ifdef lapack95
          ELSE
!           +++ F95 +++
!             CALL LA_GBSV(BL,BX,INFO=ierr_la) ! for self-compiled LAPACK95
             CALL GBSV(BL,BX,INFO=ierr_la) ! for intel mkl LAPACK95
!             CALL GBTRF(BL,kl,m,ipiv,ierr_la)
!             CALL GBTRS(BL,BX,ipiv,kl,'N',ierr_la)
             IF(ierr_la /= 0) THEN
                WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : GBSV, NT = ',  &
                     &              NT, ' -', IC, ' step. IERR=',ierr_la
                IERR = 1
                XN = XP
                GOTO 180
             ENDIF
#endif
          END IF

          ! Copy calculated variables' vector to variable matrix
          forall (NR = 0:NRMAX, NQ = 1:NQMAX) XN(NR,NQ) = BX(NQMAX * NR + NQ)

          ! Avoid negative values
!          CALL MINUS_CHECK(XN,LQb1,0)
          CALL MINUS_CHECK(XN,LQb1,3)
          CALL MINUS_CHECK(XN,LQn1,1)
          CALL MINUS_CHECK(XN,LQn2,1)
          IF(FSRP /= 0.D0) CALL MINUS_CHECK(XN,LQr1,2)
          ! Ignore tiny values
          DO NQ = 1, NQMAX
             where( abs(XN(:,NQ)) < tiny_array(NQ) ) XN(:,NQ) = 0.d0
          END DO
!!$          ! In the case of NBI off after NBI on
!!$          IF(PNBH == 0.D0) CALL THRESHOLD(XN(:,LQb1),ID)
!!$          IF(ID == 1) THEN
!!$             X (:,LQb1) = 0.D0
!!$             XP(:,LQb1) = 0.D0
!!$             X (:,LQb3) = 0.D0
!!$             XP(:,LQb3) = 0.D0
!!$             XN(:,LQb3) = 0.D0
!!$             X (:,LQb4) = 0.D0
!!$             XP(:,LQb4) = 0.D0
!!$             XN(:,LQb4) = 0.D0
!!$          END IF

          ! Check negative density or temperature in variable matrix
          CALL TXCHCK(NT,IC,XN,IERR)
          IF (IERR /= 0) THEN
             X = XN
             CALL TXCALV(X)
!             CALL TXCALC(IC)
             CALL TXGLOB
             CALL TX_GRAPH_SAVE
             RETURN
          END IF

!          IF(IC == 1) EXIT L_IC

          ! Convergence check
          call check_convergence(IC,IDIV,istat)
          if(istat == 0) then ! going to next time step
             exit L_IC
          else if(istat == 1) then
              XN = (1.d0 - oldmix) * XN + oldmix * XP
             cycle L_IC
          end if

       END DO L_IC

       ICSUM = ICSUM + IC

       IF(IDIAGL >= 2) THEN
          IF(MODECV == 0) THEN
             WRITE(6,'(A2,3(A,X),8X,A,15X,A,12X,A,11X,A)') &
                  & "**","IC","VEQ","VNR","XP","XN","V_ERRMAX","EPS"
             IF(IC == ICMAX) THEN ! NOT converged
                WRITE(6,'(3I4,1P4E17.8,A2)') IC,IASG(2),IASG(1)-1,XP(IASG(1)-1,IASG(2)), &
                     &                XN(IASG(1)-1,IASG(2)),ASG(IASG(1)-1,IASG(2)),EPSabs," *"
             ELSE ! converged
                IASG(1:2) = MAXLOC(ASG(:,:))
                WRITE(6,'(3I4,1P4E17.8)') IC,IASG(2),IASG(1)-1,XP(IASG(1)-1,IASG(2)), &
                     &                  XN(IASG(1)-1,IASG(2)),ASG(IASG(1)-1,IASG(2)),EPSabs
             END IF
          ELSE
             WRITE(6,'(A,2(A,X),5X,A,7X,A)') &
                  & "**","IC","VEQ","V_ERRMAX","EPS"
             WRITE(6,'(2I4,1PE17.8,1PE10.2)') &
                  & IC,MAXLOC(L2)-1,MAXVAL(L2),EPSabs
          END IF
       END IF

       IF(IDIAG >= 10) THEN
          IF(MODECV == 0) THEN
             WRITE(6,'(A,2(X,A,X),5X,A,11X,A)') &
                  & "*****","NQ","VR","V_ERRMAX","EPS"
             DO NQ = 1, NQMAX
                IASG(1:2)  = MAXLOC(ASG(:,NQ:NQ))
                WRITE(MSG_NQ,'(4X,2I4,1P2E17.8)') NQ,IASG(1)-1,ASG(IASG(1)-1,NQ),EPSabs
                IF( ASG(IASG(1)-1,NQ) > EPSabs) THEN
                   MSG_NQ = trim(MSG_NQ)//' *'
                ELSE
                   MSG_NQ = trim(MSG_NQ)//'  '
                END IF
                WRITE(6,'(A80)') MSG_NQ
             END DO
          ELSE
             WRITE(6,'(2(A,X),5X,A,7X,A)') &
                  & "*********","NQ","V_ERRMAX","EPS"
             DO NQ = 1, NQMAX
                WRITE(MSG_NQ,'(8X,I4,1PE17.8,1PE10.2)') NQ,L2(NQ),EPSabs
                IF(L2(NQ) > EPSabs) THEN
                   MSG_NQ = trim(MSG_NQ)//' *'
                ELSE
                   MSG_NQ = trim(MSG_NQ)//'  '
                END IF
                WRITE(6,'(A80)') MSG_NQ
             END DO
          END IF
       END IF
          
       IF(IDIAGL >= 4 .AND. MODECV == 0) THEN
          do nq = 1, NQMAX
             do nr = 0, NRMAX
                write(6,*) nq,nr,ASG(nr,nq)
             end do
          end do
       END IF

       if(istat == 2) then
          IERR = 1
          WRITE(6,'(2(A,I6))') '### ERROR(TXLOOP) : Solutions not converged at ',  &
               &              NT, ' -', IC
          goto 180
       end if

       ! Save past X for BDF
       IF(IGBDF /= 0) XOLD = X

       ! Calculation fully converged
       X = XN

       ! Calculate mesh and coefficients at the next step
       CALL TXCALV(X,1) ! Set new values to pres0 and ErV0
!!$       PNsV_FIX(:,1) = Var(:,1)%n
!!$       PTsV_FIX(:,1) = Var(:,1)%T
!!$       PNsV_FIX(:,2) = Var(:,2)%n
!!$       PTsV_FIX(:,2) = Var(:,2)%T
!!$       ErV_FIX (:) = ErV (:)
       CALL TXCALC(IC)

       IF(IDIAGL == 0 .OR. IDIAGL == 2) THEN
          IF ((MOD(NT, NTSTEP) == 0) .AND. (NT /= NTMAX)) &
               & WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC
       ELSE IF(IDIAGL > 0) THEN
          IF(NT /= NTMAX) THEN
             IF(IC-1 == ICMAX) THEN
                WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3,"  *")') NT,T_TX,IC
             ELSE
                WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC
             END IF
          END IF
       END IF

180    IF (MOD(NT, NGRSTP) == 0) CALL TXSTGR(NGR,GT,GY,NGRM)

       IF (MOD(NT, NGTSTP) == 0) THEN
          CALL TXGLOB
          CALL TXSTGT(REAL(T_TX))
          call txstgq !!!temporary
          if(IDIAG < 0) call steady_check
       END IF

       IF (MOD(NT, NGVSTP) == 0) CALL TXSTGV(REAL(T_TX))

       IF (MOD(NT, NTMAX ) == 0) CALL TXSTGQ

!       call cal_flux

       IF (IERR /= 0) EXIT L_NTDO

    END DO L_NTDO

    rIPs = rIPe

    IF(IC == ICMAX) THEN
       WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3,"  *")') NT,T_TX,IC
    ELSE
       WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC
    END IF
    IF(NTMAX /= 0) THEN
       AVE_IC = REAL(ICSUM) / NTMAX
    ELSE
       AVE_IC = 0.d0
    END IF

    deallocate(BA,BL,BX)
    deallocate(XN,XP,ASG,L2)

    RETURN
  END SUBROUTINE TXLOOP

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

  SUBROUTINE TXCALB(BA,BL,BX)

    use tx_commons, only : IGBDF, ADV, MDLPCK, NQMAX, NRMAX, NLC, NLCR, CLC, BLC, ALC, NLCMAX, &
         &              PLC, X, XOLD, NCM!, lqb1, pnbv
    real(8), dimension(:,:), intent(inout) :: BA, BL
    real(8), dimension(:), intent(inout) :: BX
    INTEGER(4) :: J, NR, NQ, NC, NC1, IA, IB, IC, NCHvs, NC1Hvs, NC1Hvs2, NC2
    INTEGER(4) :: JA, JB, JC, KL
    real(8), parameter :: C43 = 4.d0/3.d0, C23 = 2.d0/3.d0, C13 = 1.d0/3.d0
    real(8) :: COEF1, COEF2, COEF3, COEF!, suml1,suml2
    
    IF(IGBDF /= 0) ADV = C23

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
    IF(MDLPCK == 0) THEN
       BA = 0.d0 ! initialize BA array

!!$       DO NQ = 1, NQMAX
!!$          DO NC = 0, NLCMAX(NQ)
!!$             NCHvs = real( NCM - 1 + NC ) / NCM ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
!!$             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
!!$             NR = 0
!!$                NC1 = NLCR(NC,NQ,0)
!!$                IF(NC1 /= 0) THEN
!!$                   IC = NQMAX + (NC1 - 1) - (NQ - 1)
!!$                   IB = IC + NQMAX
!!$                   IA = IB + NQMAX
!!$                   J = NR * NQMAX + NQ
!!$                   BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef
!!$                   BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef
!!$                   BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef
!!$                END IF
!!$             NC1 = NLC(NC,NQ)
!!$             IF(NC1 /= 0) THEN
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
!!$                IF(NC1 /= 0) THEN
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

!$omp parallel
!$omp do private(NC,NCHvs,coef,NC1,NC1Hvs,NC1Hvs2,IA,IB,IC,NR,J)
       DO NQ = 1, NQMAX
          DO NC = 0, NLCMAX(NQ)
             NCHvs = ( NCM - 1 + NC ) / NCM ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
             ! --- NR = 0
             NC1 = NLCR(NC,NQ,0)
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
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
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC1Hvs2 = 1 - NC1Hvs
             ! NC1Hvs2 in IC: avoid IC=0 in BL when NR = NQMAX and NC1 = 0
             IC = NQMAX + (NC1 - 1) - (NQ - 1) + NC1Hvs2
             IB = IC + NQMAX
             IA = IB + NQMAX
             DO NR = 1, NRMAX-1
                J = NR * NQMAX + NQ
                BA(IC,J) = BA(IC,J) + CLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IB,J) = BA(IB,J) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BA(IA,J) = BA(IA,J) + ALC(NR,NC,NQ) * coef * NC1Hvs
             END DO
             ! --- NR = NRMAX
             NC1 = NLCR(NC,NQ,1)
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
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
          END DO
       END DO
!$omp end do
!$omp end parallel

    ! For LAPACK solver
    ELSE
       BL = 0.d0 ! initialize BL array
       KL = 2 * NQMAX - 1

!!$       DO NQ = 1, NQMAX
!!$          DO NC = 0, NLCMAX(NQ)
!!$             NCHvs = ( NCM - 1 + NC ) / NCM ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
!!$             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
!!$             NR = 0
!!$                NC1 = NLCR(NC,NQ,0)
!!$                NC1Hvs = ( NCM - 1 + NC1 ) / NCM
!!$                IF(NC1 /= 0) THEN
!!$                   IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
!!$                   IB = IA + NQMAX
!!$                   JA =(NR + 1) * NQMAX + NC1
!!$                   JB = NR      * NQMAX + NC1
!!$                   BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef
!!$                   BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef
!!$                END IF
!!$             NC1 = NLC(NC,NQ)
!!$             IF(NC1 /= 0) THEN
!!$                IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
!!$                IB = IA + NQMAX
!!$                IC = IB + NQMAX
!!$                DO NR = 1, NRMAX-1
!!$                   JA =(NR + 1) * NQMAX + NC1
!!$                   JB = NR      * NQMAX + NC1
!!$                   JC =(NR - 1) * NQMAX + NC1
!!$                   BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef
!!$                   BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef
!!$                   BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef
!!$                END DO
!!$             END IF
!!$             NR = NRMAX
!!$                NC1 = NLCR(NC,NQ,1)
!!$                IF(NC1 /= 0) THEN
!!$                   IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
!!$                   IB = IA + NQMAX
!!$                   IC = IB + NQMAX
!!$                   JB = NR      * NQMAX + NC1
!!$                   JC =(NR - 1) * NQMAX + NC1
!!$                   BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef
!!$                   BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef
!!$                END IF
!!$          END DO
!!$       END DO

       ! *** Introducing NCHvs makes it possible to eliminate NC branch (NC=0 or else).
       ! *** Introducing NC1Hvs and NC1Hvs2 makes it possible to eliminate IF(NC1 /= 0) statement in do loops

!$omp parallel
!$omp do private(NC,NCHvs,coef,NC1,NC1Hvs,NC1Hvs2,IA,IB,IC,JA,JB,JC,NR)
       DO NQ = 1, NQMAX
          DO NC = 0, NLCMAX(NQ)
             NCHvs = ( NCM - 1 + NC ) / NCM ! NCHvs = 0 when NC = 0 ; otherwise NCHvs = 1
             coef  = 1.d0 - ( 1.d0 + adv ) * NCHvs
             ! --- NR = 0
             NC1 = NLCR(NC,NQ,0)
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
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
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM
             NC1Hvs2 = 1 - NC1Hvs
             ! NC1Hvs2 in IA: avoid overflow of IA in BL when NQ=NQMAX and NC1 = 0
             IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL - NC1Hvs2
             IB = IA + NQMAX
             IC = IB + NQMAX
             DO NR = 1, NRMAX-1
                JA =(NR + 1) * NQMAX + NC1
                JB = NR      * NQMAX + NC1
                ! NC1Hvs2 in JC: avoid overflow of JC in BL when NR=NRMAX-1 and NC1 = 0
                JC =(NR - 1) * NQMAX + NC1 + NC1Hvs2
                BL(IA,JA) = BL(IA,JA) + ALC(NR,NC,NQ) * coef * NC1Hvs
                BL(IB,JB) = BL(IB,JB) + BLC(NR,NC,NQ) * coef * NC1Hvs
                BL(IC,JC) = BL(IC,JC) + CLC(NR,NC,NQ) * coef * NC1Hvs
             END DO
             ! --- NR = NRMAX
             NC1 = NLCR(NC,NQ,1)
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM
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
          END DO
       END DO
!$omp end do
!$omp end parallel

    END IF

    ! *** Right-hand-side vector, denoted by "bu" ***

    BX = 0.d0 ! initialize BX array

    IF(IGBDF == 0) THEN
       COEF1 = 1.D0
       COEF2 = 0.D0
       COEF3 = 1.D0
    ELSE
       COEF1 = C43
       COEF2 = C13
       COEF3 = C23
    END IF

    ! In the case of NC=0 i.e. time derivative term effects in any equations
!!$    NC = 0
!!$    NR = 0
!!$       DO NQ = 1, NQMAX
!!$          NC1 = NLCR(NC,NQ,0)
!!$          IF(NC1 /= 0) THEN
!!$             BX(NQMAX * NR + NQ) &
!!$                  & = BX(NQMAX * NR + NQ) + BLC(NR,NC,NQ) * X   (NR  ,NC1) * COEF1 &
!!$                  &                       + ALC(NR,NC,NQ) * X   (NR+1,NC1) * COEF1 &
!!$                  &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC1) * COEF2 &
!!$                  &                       - ALC(NR,NC,NQ) * XOLD(NR+1,NC1) * COEF2
!!$          END IF
!!$       END DO
!!$
!!$    DO NQ = 1, NQMAX
!!$       NC1 = NLC(NC,NQ)
!!$       DO NR = 1, NRMAX - 1
!!$          IF(NC1 /= 0) THEN
!!$             BX(NQMAX * NR + NQ) &
!!$                  &    = BX(NQMAX * NR + NQ) + CLC(NR,NC,NQ) * X   (NR-1,NC1) * COEF1 &
!!$                  &                          + BLC(NR,NC,NQ) * X   (NR  ,NC1) * COEF1 &
!!$                  &                          + ALC(NR,NC,NQ) * X   (NR+1,NC1) * COEF1 &
!!$                  &                          - CLC(NR,NC,NQ) * XOLD(NR-1,NC1) * COEF2 &
!!$                  &                          - BLC(NR,NC,NQ) * XOLD(NR  ,NC1) * COEF2 &
!!$                  &                          - ALC(NR,NC,NQ) * XOLD(NR+1,NC1) * COEF2
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    NR = NRMAX
!!$       DO NQ = 1, NQMAX
!!$          NC1 = NLCR(NC,NQ,1)
!!$          IF(NC1 /= 0) THEN
!!$             BX(NQMAX * NR + NQ) &
!!$                  & = BX(NQMAX * NR + NQ) + CLC(NR,NC,NQ) * X   (NR-1,NC1) * COEF1 &
!!$                  &                       + BLC(NR,NC,NQ) * X   (NR  ,NC1) * COEF1 &
!!$                  &                       - CLC(NR,NC,NQ) * XOLD(NR-1,NC1) * COEF2 &
!!$                  &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC1) * COEF2
!!$          END IF
!!$       END DO

!$omp parallel
    NC = 0
!$omp do private(NC1,NC1Hvs,NC2,NR)
    DO NQ = 1, NQMAX
       ! --- NR = 0
       NC1 = NLCR(NC,NQ,0)
       NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
       NC2 = NC1 + 1 - NC1Hvs
       NR = 0
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) +(BLC(NR,NC,NQ) * X   (NR  ,NC2) * COEF1 &
               &                       + ALC(NR,NC,NQ) * X   (NR+1,NC2) * COEF1 &
               &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC2) * COEF2 &
               &                       - ALC(NR,NC,NQ) * XOLD(NR+1,NC2) * COEF2) * NC1Hvs

       ! --- NR = 1 ~ NRMAX
       NC1 = NLC(NC,NQ)
       NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
       NC2 = NC1 + 1 - NC1Hvs
       DO NR = 1, NRMAX - 1
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) +(CLC(NR,NC,NQ) * X   (NR-1,NC2) * COEF1 &
               &                       + BLC(NR,NC,NQ) * X   (NR  ,NC2) * COEF1 &
               &                       + ALC(NR,NC,NQ) * X   (NR+1,NC2) * COEF1 &
               &                       - CLC(NR,NC,NQ) * XOLD(NR-1,NC2) * COEF2 &
               &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC2) * COEF2 &
               &                       - ALC(NR,NC,NQ) * XOLD(NR+1,NC2) * COEF2) * NC1Hvs
       END DO

       ! --- NR = NRMAX
       NC1 = NLCR(NC,NQ,1)
       NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
       NC2 = NC1 + 1 - NC1Hvs
       NR = NRMAX
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) +(CLC(NR,NC,NQ) * X   (NR-1,NC2) * COEF1 &
               &                       + BLC(NR,NC,NQ) * X   (NR  ,NC2) * COEF1 &
               &                       - CLC(NR,NC,NQ) * XOLD(NR-1,NC2) * COEF2 &
               &                       - BLC(NR,NC,NQ) * XOLD(NR  ,NC2) * COEF2) * NC1Hvs
    END DO
!$omp end do
!$omp end parallel

    ! *** In the case of general term effect in any equations, denoted by "c" ***

    forall (NQ = 1:NQMAX, NR = 0:NRMAX) &
         & BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) + sum(PLC(NR,1:NLCMAX(NQ),NQ)) * COEF3

    ! *** Only the case of not BDF ***
    !  Because "ADV" has no longer original meaning when BDF is used.

    IF(IGBDF == 0) THEN
!!$    DO NQ = 1, NQMAX
!!$       NR = 0
!!$          DO NC = 1, NLCMAX(NQ)
!!$             NC1 = NLCR(NC,NQ,0)
!!$             IF(NC1 /= 0) THEN
!!$                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
!!$                     &              +(  BLC(NR,NC,NQ) * X(NR  ,NC1) &
!!$                     &                + ALC(NR,NC,NQ) * X(NR+1,NC1)) * (1.D0 - ADV)
!!$             END IF
!!$          END DO
!!$
!!$       DO NR = 1, NRMAX-1
!!$          DO NC = 1, NLCMAX(NQ)
!!$             NC1 = NLC(NC,NQ)
!!$             IF(NC1 /= 0) THEN
!!$                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
!!$                     &              +(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
!!$                     &                + BLC(NR,NC,NQ) * X(NR  ,NC1) &
!!$                     &                + ALC(NR,NC,NQ) * X(NR+1,NC1)) * (1.D0 - ADV)
!!$             END IF
!!$          END DO
!!$       END DO
!!$
!!$       NR = NRMAX
!!$          DO NC = 1, NLCMAX(NQ)
!!$             NC1 = NLCR(NC,NQ,1)
!!$             IF(NC1 /= 0) THEN
!!$                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
!!$                     &              +(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
!!$                     &                + BLC(NR,NC,NQ) * X(NR  ,NC1)) * (1.D0 - ADV)
!!$             END IF
!!$          END DO
!!$    END DO

!$omp parallel do private(NR,NC,NC1,NC1Hvs,NC2)
    DO NQ = 1, NQMAX
       NR = 0
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,0)
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC2 = NC1 + 1 - NC1Hvs
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                  &              +(  BLC(NR,NC,NQ) * X(NR  ,NC2) &
                  &                + ALC(NR,NC,NQ) * X(NR+1,NC2)) * (1.D0 - ADV) * NC1Hvs
          END DO

       DO NC = 1, NLCMAX(NQ)
          NC1 = NLC(NC,NQ)
          NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
          NC2 = NC1 + 1 - NC1Hvs
          DO NR = 1, NRMAX-1
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                  &              +(  CLC(NR,NC,NQ) * X(NR-1,NC2) &
                  &                + BLC(NR,NC,NQ) * X(NR  ,NC2) &
                  &                + ALC(NR,NC,NQ) * X(NR+1,NC2)) * (1.D0 - ADV) * NC1Hvs
          END DO
       END DO

       NR = NRMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,1)
             NC1Hvs = ( NCM - 1 + NC1 ) / NCM ! NC1Hvs = 0 when NC1 = 0 ; otherwise NC1Hvs = 1
             NC2 = NC1 + 1 - NC1Hvs
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                  &              +(  CLC(NR,NC,NQ) * X(NR-1,NC2) &
                  &                + BLC(NR,NC,NQ) * X(NR  ,NC2)) * (1.D0 - ADV) * NC1Hvs
          END DO
    END DO
!$omp end parallel do
    END IF

    RETURN
  END SUBROUTINE TXCALB

!***************************************************************
!
!   Check negative density or temperature
!
!***************************************************************

  SUBROUTINE TXCHCK(NTL,IC,XL,IER)

    use tx_commons, only : NQMAX, NRMAX, LQe1, LQi1, LQe5, LQi5, LQb1, LQr1
    INTEGER(4), intent(in) :: NTL, IC
    integer(4), intent(inout) :: IER
    REAL(8), DIMENSION(0:NRMAX,1:NQMAX), intent(in) :: XL
    integer(4) :: NR

    IER = 0

    DO NR = 0, NRMAX
       IF (XL(NR,LQe1) < 0.D0 .OR. XL(NR,LQi1) < 0.D0 .OR. &
        &  XL(NR,LQb1) < 0.D0 .OR. XL(NR,LQr1) < 0.D0) THEN
          WRITE(6,'(2A,I4,2(A,I4),A)') '### ERROR(TXLOOP) : Negative density at ', &
               &           'NR =', NR, ', NT=', NTL, ', IC=', IC, '.'
          WRITE(6,'(1X,4(A,1PE10.3))') 'ne =',  REAL(XL(NR,LQe1)), &
               &                    ',   ni =', REAL(XL(NR,LQi1)), &
               &                    ',   nb =', REAL(XL(NR,LQb1)), &
               &                    ', nbrp =', REAL(XL(NR,LQr1))
          IER = 1
          RETURN
       END IF
    END DO

    DO NR = 0, NRMAX
       IF (XL(NR,LQe5) < 0.D0 .OR. XL(NR,LQi5) < 0.D0) THEN
          WRITE(6,'(2A,I4,2(A,I4),A)') '### ERROR(TXLOOP) : Negative temperature at ', &
               &           'NR =', NR, ', NT=', NTL, ', IC=', IC, '.'
          WRITE(6,'(20X,2(A,1PE15.7))') 'Te =', REAL(XL(NR,LQe5)), &
               &           ',   Ti =', REAL(XL(NR,LQi5))
          IER = 2
          RETURN
       END IF
    END DO

    RETURN
  END SUBROUTINE TXCHCK

!***************************************************************
!
!   Negative value forced to be set to zero for densities
!
!***************************************************************

  SUBROUTINE MINUS_CHECK(XL,LQ,ID)
    
    use tx_commons, only : NQMAX, NRMAX, NRA
    integer(4), intent(in) :: LQ, ID
    real(8), dimension(0:NRMAX,1:NQMAX), intent(inout) :: XL
    integer(4) :: NR, NZERO, ind

    ind = 0
    IF(ID == 0) THEN ! Negative values are set to zero inside rho=1.
       IF(MINVAL(XL(:,LQ)) < 0.D0) THEN
          ind = 1
          DO NR = 0, NRMAX
             IF(XL(NR,LQ) <= 0.D0) THEN
                IF(NR < NRA) THEN
                   XL(NR,LQ) = 0.D0
                ELSE
                   NZERO = NR
                   EXIT
                END IF
             END IF
          END DO

          XL(NZERO:NRMAX,LQ) = 0.D0
       END IF
    ELSE IF(ID == 1) THEN ! Searching from outboard
       IF(MINVAL(XL(:,LQ)) < 0.D0) THEN
          ind = 1
          DO NR = NRMAX, 0, -1
             IF(XL(NR,LQ) <= 0.D0) THEN
                NZERO = NR
                EXIT
             END IF
          END DO

          XL(0:NZERO,LQ) = 0.D0
       END IF
    ELSE IF(ID == 2) THEN ! Negative values are set to zero.
       IF(MINVAL(XL(:,LQ)) < 0.D0) THEN
          ind = 1
          where(XL(:,LQ) < 0.d0) XL(:,LQ) = 0.d0
       END IF
    ELSE ! Taking absolute
       XL(:,LQ) = abs(XL(:,LQ))
    END IF
!    if(ind /= 0 .and. (ID == 1 .or. ID == 2)) write(6,'(3(X,A,I3,X))') "MINUS_CHECK : LQ=",LQ,"ID=",ID,"NR=",NZERO

  END SUBROUTINE MINUS_CHECK

!!$  SUBROUTINE THRESHOLD(XL,ID)
!!$
!!$    real(8), dimension(:), intent(inout) :: XL
!!$    integer(4), intent(out) :: ID
!!$    integer(4) :: NR, NRL
!!$
!!$    ID = 0
!!$    NRL = 0.5 * NRA
!!$    IF(MINVAL(XL(0:NRL)) < 1.D-8 .AND. MAXVAL(XL(0:NRL)) > 0.D0) THEN
!!$       XL(:) = 0.D0
!!$       ID = 1
!!$    END IF
!!$
!!$  END SUBROUTINE THRESHOLD

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

    IDIAGL = MOD(IDIAG,10)
    EPSabs = abs(EPS)

    IF(MODECV == 0) THEN

       ASG = 0.D0
       L_NQ:DO NQ = 1, NQMAX
          ! Calculate maximum local root-mean-square of X and corresponding grid point
          NRAVM = SUM(MAXLOC(abs(XN(:,NQ))))-1
          AVM   = abs(XN(NRAVM,NQ))

          ! Calculate root-mean-square X over the profile (Euclidean norm)
          AV    = SQRT(SUM(XN(:,NQ)**2) / NRMAX)
          IF (AV < epsilon(1.d0)) THEN
             ! It means that the variable for NQ is zero over the profile.
             ASG(:,NQ) = 0.d0
             CYCLE L_NQ
          END IF

          ERR1  = 0.D0
          IDISP = IDIV
          ASG(:,NQ) = ABS(XN(:,NQ) - XP(:,NQ)) / AV

          L_NR:DO NR = 0, NRMAX
             ! Show results
             IF (NT == IDISP .AND. NR == NRMAX) THEN
                IF (NQ == 1) THEN
                   WRITE(6,'(1X,A5,A3," =",I3)')'#####','IC',IC
                   WRITE(6,'(1X,A7," =",1PD9.2,2X,A7," =",I3)') &
                        & 'EPS    ', EPSabs,'NRMAX  ', NRMAX
                END IF
                ! Maximum relative error over the profile
                ERR1 = MAX(ERR1, ASG(NR,NQ))
                WRITE(6,'((1X,A2," =",I2,2X,A2," =",1PD9.2, &
                     & 2X,A5," =",1PD9.2,A1,I2,2X,A5," =",1PD9.2))') &
                     & 'NQ    ', NQ , &
                     & 'AV    ', AV    ,  'AVMAX ', AVM   ,':',NRAVM, &
                     & 'SCMAX ', ERR1
                IDISP = IDIV + NT
                if (NQ == NQMAX) then ! not converged
                   istat = 1
                   return
                end IF
             ELSEIF (NT /= IDISP .AND. ASG(NR,NQ) > EPSabs) THEN
                ! NOT converged
                IF(IC /= ICMAX) THEN
                   IF(IDIAGL >= 3) THEN
                      IF(MOD(IC,20) == 0 .OR. IC == 1) &
                           & WRITE(6,'(A2,3(A,X),8X,A,15X,A,12X,A,11X,A)') &
                           & "**","IC","VEQ","VNR","XP","XN","V_ERRMAX","EPS"
                      WRITE(6,'(3I4,1P4E17.8)') IC,NQ,NR,XP(NR,NQ),XN(NR,NQ), &
                           & ASG(NR,NQ),EPSabs
                   END IF
                   istat = 1
                   return
                END IF
             END IF
          END DO L_NR

       END DO L_NQ
       ! Converged or IC == ICMAX
       IASG(1:2) = MAXLOC(ASG)

       if(eps > 0.d0) then ! calculation continues even if not converged
          istat = 0
       else ! calculation stops if not converged
          if(istat == 1) then
             istat = 2 ! not converged at IC = ICMAX
          else
             istat = 0 ! converged
          end if
       end if
!       if(EPS > 0.D0) EXIT L_IC

    ELSE
       L_NQ2:DO NQ = 1, NQMAX
          AV = abs(sum(XN(:,NQ)))
          IF(AV /= 0.d0) THEN
             ASG(:,NQ) = (XN(:,NQ) - XP(:,NQ))**2
             L2(NQ) = SQRT(sum(ASG(:,NQ))) / AV
          ELSE
             L2(NQ) = 0.D0
          END IF
       END DO L_NQ2

       ! Converged
       IF(MAXVAL(L2) < EPSabs) THEN
!          if(EPS > 0.D0) EXIT L_IC
          if(eps > 0.d0) then ! calculation continues even if not converged
             istat = 0
             return
          end if
       END IF

       IF(IC /= ICMAX) THEN
          IF(IDIAGL >= 3) THEN
             IF(MOD(IC,20) == 0 .OR. IC == 1) &
                  & WRITE(6,'(A,2(A,X),5X,A,7X,A)') &
                  & "**","IC","VEQ","V_ERRMAX","EPS(1)"
             WRITE(6,'(2I4,1PE17.8,1PE10.2)') &
                  & IC,MAXLOC(L2)-1,MAXVAL(L2),EPSabs
          END IF
       ELSE
          istat = 2
!          EXIT L_IC ! This is used so that IC exceeds ICMAX after "END DO L_IC".
       END IF
    END IF

  end subroutine CHECK_CONVERGENCE

end module tx_main
