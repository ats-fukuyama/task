!     $Id$
module tx_main
  implicit none
  private
  public :: TXEXEC

contains

!***************************************************************
!
!   MAIN ROUTINE
!
!***************************************************************

  SUBROUTINE TXEXEC
    use tx_commons, only : IERR, T_TX
    use tx_interface, only : APTOS

    INTEGER(4) :: NDY, NDM, NDD, NTH, NTM, NTS, NSTR1, NSTR2, NSTR3
    REAL(4) :: gCTIME1, gCTIME2, gCTIME3
    character(len=10) :: STR1, STR2, STR3
    INTEGER(4), DIMENSION(1:8) :: TIMES
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
    CALL APTOS(STR2, NSTR2, SNGL(t_interval), 'E1')
    IF(gCTIME3 /= 0) THEN
       NSTR3 = 0
       CALL APTOS(STR3, NSTR3, SNGL(t_interval) / gCTIME3 * 100, 'F5')
       WRITE(6,*) 'CPU = ', STR1(1:NSTR1), ' (sec)   ',  &
            &     'sim time = ', STR2(1:NSTR2), ' (sec)', '  (', STR3(1:NSTR3), '%)'
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
    use tx_commons, only : T_TX, rIPe, rIPs, NTMAX, IGBDF, NQMAX, NRMAX, X, ICMAX, &
         &                 PNeV, PTeV,PNeV_FIX, PTeV_FIX, NQM, IERR, LQb1, LQn1, &
         &                 tiny_cap, EPS, IDIAG, NTSTEP, NGRSTP, NGTSTP, NGVSTP, GT, GY, &
         &                 NGRM, NGYRM, FSRP, fmnq, wnm, umnq, nmnqm, MODEAV, XOLD, &
         &                 NT, DT, rIP, MDLPCK, NGR, MDSOLV
    use tx_variables
    use tx_coefficients, only : TXCALA
    use tx_graphic, only : TX_GRAPH_SAVE, TXSTGT, TXSTGV, TXSTGR, TXSTGQ

    real(8), dimension(:,:), allocatable :: BA, BL
    real(8), dimension(:),   allocatable :: BX, XNvec, FL, FLP, DltXN, DltXP, BAE
    INTEGER(4) :: I, J, NR, NQ, NC, NC1, IC = 0, IDIV, NTDO, IDISP, NRAVM, ID
    INTEGER(4), DIMENSION(1:NQM*(NRMAX+1)) :: IPIV
    REAL(8) :: TIME0, DIP, AVM, ERR1, AV
    real(8), dimension(1:NQMAX) :: tiny_array
    REAL(8), DIMENSION(1:NQM,0:NRMAX) :: XN, XP, ASG
    integer(4), dimension(1:2) :: iasg
    real(8) :: denom, nume, tiny_denom
    real(8), dimension(:,:), allocatable :: XN1, XN2

    allocate(BA(1:4*NQM-1,1:NQM*(NRMAX+1)),BL(1:6*NQM-2,1:NQM*(NRMAX+1)),BX(1:NQM*(NRMAX+1)))
    IF(MDSOLV == 1) THEN
       allocate(XNvec(1:NQM*(NRMAX+1)),FL(1:NQM*(NRMAX+1)),FLP(1:NQM*(NRMAX+1)))
       allocate(BAE(1:(NQM*(NRMAX+1))**2),DltXN(1:NQM*(NRMAX+1)),DltXP(1:NQM*(NRMAX+1)))
    END IF

    !  Read spline table for neoclassical toroidal viscosity if not loaded when FSRP/=0
    IF(FSRP /= 0.D0 .AND. maxval(fmnq) == 0.D0) CALL Wnm_spline(fmnq, wnm, umnq, nmnqm)

    IF (MODEAV == 0) THEN
       IDIV = NTMAX + 1
    ELSE
       IDIV = NTMAX / MODEAV
    END IF
    TIME0 = T_TX
    IF(NTMAX /= 0) DIP = (rIPe - rIPs) / NTMAX

    ! Save X -> XP -> XOLD for BDF only at the beginning of the calculation
    IF(IGBDF /= 0 .and. T_TX == 0.D0) XOLD(1:NQMAX,0:NRMAX) = X(1:NQMAX,0:NRMAX)

    L_NTDO:DO NTDO = 1, NTMAX
       NT = NTDO
       T_TX = TIME0 + DT*NT
       rIP  = rIPs  + DIP*NT

       ! Create new X = XN
       XN(1:NQMAX,0:NRMAX) = X(1:NQMAX,0:NRMAX)

       ! Negligible order of magnitude for each variable
       DO NQ = 1, NQMAX
          tiny_array(NQ) = maxval(XN(NQ,0:NRMAX)) * tiny_cap
       END DO

       L_IC : DO IC = 1, ICMAX
          ! Save past X = XP
          XP(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)
          IF(MDSOLV == 1 .and. IC == 1) THEN
             DO NR = 0, NRMAX
                DO NQ = 1, NQMAX
                   XNvec(NQMAX * NR + NQ) = XN(NQ,NR)
                END DO
             END DO
          END IF

          CALL TXCALV(XP,ID)
          IF(IC == 1) THEN
             PNeV_FIX(0:NRMAX) = PNeV(0:NRMAX)
             PTeV_FIX(0:NRMAX) = PTeV(0:NRMAX)
          END IF

          CALL TXCALC
          CALL TXCALA
          IF(MDSOLV == 1.and. IC /= 1) THEN
             ! Get FL and BX
             CALL TXCALB_SECANT(FL,BX,XNvec)
          ELSE
             ! Get BA or BL, and BX
             CALL TXCALB(BA,BL,BX)
          END IF
          CALL TXGLOB

          IF(MDSOLV /= 1 .or. (MDSOLV == 1 .and. IC == 1)) THEN
             IF(MDLPCK == 0) THEN
                CALL BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQM-1, IERR)
                IF (IERR >= 30000) THEN
                   WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : Matrix BA is singular at ',  &
                        &              NT, ' -', IC, ' step. IERR=',IERR
                   IERR = 1
                   XN(1:NQMAX,0:NRMAX) = XP(1:NQMAX,0:NRMAX)
                   GOTO 180
                END IF
             ELSE
                CALL LAPACK_DGBSV(NQMAX*(NRMAX+1),2*NQMAX-1,2*NQMAX-1,1,BL, &
                     &            6*NQMAX-2,IPIV,BX,NQMAX*(NRMAX+1),IERR)
                IF(IERR /= 0) THEN
                   WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : DGBSV, IERR = ',  &
                        &              NT, ' -', IC, ' step. IERR=',IERR
                   IERR = 1
                   XN(1:NQMAX,0:NRMAX) = XP(1:NQMAX,0:NRMAX)
                   GOTO 180
                ENDIF
             END IF
          ELSE
             IF(MDLPCK == 0) THEN
                stop 'TXLOOP: MDLPCK = 0 cannot be reconciled to MDSOLV = 1. '
             ELSE
                DO J = 1, NQMAX*(NRMAX+1)
                   DO I = 1, NQMAX*(NRMAX+1)
                      IF(DltXP(J) < tiny_cap) THEN
                         BAE((NQMAX*(NRMAX+1))*(J-1)+I) = 0.D0
                      ELSE
                         BAE((NQMAX*(NRMAX+1))*(J-1)+I) = (FL(I) - FLP(I)) / DltXP(J)
                      END IF
!                      IF(I == J) write(6,*) BAE((NQMAX*(NRMAX+1))*(J-1)+I)
                   END DO
                END DO
                DltXN(1:NQMAX*(NRMAX+1)) = - FL(1:NQMAX*(NRMAX+1))
                CALL LAPACK_DGESV(NQMAX*(NRMAX+1),1,BAE,NQMAX*(NRMAX+1), &
                     &            IPIV,DltXN,NQMAX*(NRMAX+1),IERR)
                IF(IERR /= 0) THEN
                   WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : DGESV, IERR = ',  &
                        &              NT, ' -', IC, ' step. IERR=',IERR
                   IERR = 1
                   XN(1:NQMAX,0:NRMAX) = XP(1:NQMAX,0:NRMAX)
                   GOTO 180
                ENDIF
             END IF
          END IF

          IF(MDSOLV /= 1 .or. (MDSOLV == 1 .and. IC == 1)) THEN
             ! Copy calculated variables' vector to variable matrix
             DO NR = 0, NRMAX
                DO NQ = 1, NQMAX
                   XN(NQ,NR) = BX(NQMAX * NR + NQ)
                END DO
             END DO
          ELSE
             DO NR = 0, NRMAX
                DO NQ = 1, NQMAX
                   DltXN(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ)
                   XN(NQ,NR) = XP(NQ,NR) + DltXN(NQMAX * NR + NQ)
                   XNvec(NQMAX * NR + NQ) = XN(NQ,NR)
                END DO
             END DO
          END IF

          ! Avoid negative values
          CALL MINUS_GOES_ZERO(XN,LQb1,0)
          CALL MINUS_GOES_ZERO(XN,LQn1,1)
          ! Ignore tiny values
          DO NQ = 1, NQMAX
             DO NR = 0, NRMAX
                if(abs(XN(NQ,NR)) < tiny_array(NQ)) XN(NQ,NR) = 0.d0
             END DO
          END DO
!!$          ! In the case of NBI off after NBI on
!!$          IF(PNBH == 0.D0) CALL THRESHOLD(XN(LQb1,0:NRMAX),ID)
!!$          IF(ID == 1) THEN
!!$             X (LQb1,0:NRMAX) = 0.D0
!!$             XP(LQb1,0:NRMAX) = 0.D0
!!$             X (LQb3,0:NRMAX) = 0.D0
!!$             XP(LQb3,0:NRMAX) = 0.D0
!!$             XN(LQb3,0:NRMAX) = 0.D0
!!$             X (LQb4,0:NRMAX) = 0.D0
!!$             XP(LQb4,0:NRMAX) = 0.D0
!!$             XN(LQb4,0:NRMAX) = 0.D0
!!$          END IF

          ! Check negative density or temperature in variable matrix
          CALL TXCHCK(NT,IC,XN,IERR)
          IF (IERR /= 0) THEN
             X(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)
             CALL TXCALV(X)
             CALL TXCALC
             CALL TXGLOB
             CALL TX_GRAPH_SAVE
             RETURN
          END IF

          IF(MDSOLV == 1) THEN
             FLP(1:NQMAX*(NRMAX+1)) = FL(1:NQMAX*(NRMAX+1))
             DO NR = 0, NRMAX
                DO NQ = 1, NQMAX
                   XNvec(NQMAX * NR + NQ) = XN(NQ,NR)
                   DltXP(NQMAX * NR + NQ) = XNvec(NQMAX * NR + NQ) - XP(NQ,NR)
                END DO
             END DO
          END IF

          if(MDSOLV == 2) then

             ! Steffensen's method (not completed)

             tiny_denom = tiny_cap * EPS

             allocate(XN1(1:NQMAX,0:NRMAX),XN2(1:NQMAX,0:NRMAX))

             ! Save past X = XN1
             XN1(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)

             CALL TXCALV(XN1,ID)
             CALL TXCALC
             CALL TXCALA
             CALL TXCALB(BA,BL,BX)
             CALL TXGLOB
             IF(MDLPCK == 0) THEN
                CALL BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQM-1, IERR)
                IF (IERR >= 30000) THEN
                   WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : Matrix BA is singular at ',  &
                        &              NT, ' -', IC, ' step. IERR=',IERR
                   IERR = 1
                   XN(1:NQMAX,0:NRMAX) = XN1(1:NQMAX,0:NRMAX)
                   GOTO 180
                END IF
             ELSE
                CALL LAPACK_DGBSV(NQMAX*(NRMAX+1),2*NQMAX-1,2*NQMAX-1,1,BL, &
                     &            6*NQMAX-2,IPIV,BX,NQMAX*(NRMAX+1),IERR)
                IF(IERR /= 0) THEN
                   WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : DGBSV, IERR = ',  &
                        &              NT, ' -', IC, ' step. IERR=',IERR
                   IERR = 1
                   XN(1:NQMAX,0:NRMAX) = XN1(1:NQMAX,0:NRMAX)
                   GOTO 180
                ENDIF
             END IF

             ! Copy calculated variables' vector to variable matrix
             DO NR = 0, NRMAX
                DO NQ = 1, NQMAX
                   XN2(NQ,NR) = BX(NQMAX * NR + NQ)
                END DO
             END DO
             ! Avoid negative values
             CALL MINUS_GOES_ZERO(XN2,LQb1,0)
             CALL MINUS_GOES_ZERO(XN2,LQn1,1)
             ! Ignore tiny values
             DO NQ = 1, NQMAX
                DO NR = 0, NRMAX
                   if(abs(XN2(NQ,NR)) < tiny_array(NQ)) XN2(NQ,NR) = 0.d0
                END DO
             END DO
!             where (abs(xn2) < tiny_cap) xn2 = 0.d0

             do NQ = 1, NQMAX
                do NR = 0, NRMAX
                   if(abs(XN(NQ,NR) - XP(NQ,NR)) > tiny_cap) then
                      denom = XP(NQ,NR) - 2.D0 * XN1(NQ,NR) + XN2(NQ,NR)
                      nume  = (XN1(NQ,NR) - XP(NQ,NR))**2
                      if(abs(denom) < tiny_denom) then
                         XN(NQ,NR) = 0.D0
!                         write(6,*) "PASS: NQ=",NQ,"  NR=",NR
!                         write(6,*) denom, tiny_denom
                      else
                         XN(NQ,NR) = XP(NQ,NR) - nume / denom
                      end if
                   end if
                end do
             end do

             deallocate(XN1,XN2)

          end if

!          IF(IC == 1) EXIT L_IC

          ASG(1:NQMAX,0:NRMAX) = 0.D0
          L_NQ:DO NQ = 1, NQMAX
             ! Calculate maximum local root-mean-square of X and corresponding grid point
             NRAVM = SUM(MAXLOC(SQRT(XN(NQ,0:NRMAX)**2)))-1
             AVM   = SQRT(XN(NQ,NRAVM)**2)

             ! Calculate root-mean-square X over the profile (Euclidean norm)
             AV    = SQRT(SUM(XN(NQ,0:NRMAX)**2)) / NRMAX
             IF (AV == 0.D0) CYCLE L_NQ

             ERR1  = 0.D0
             IDISP = IDIV
             ASG(NQ,0:NRMAX) = ABS(XN(NQ,0:NRMAX) - XP(NQ,0:NRMAX)) / AV
!!$             IF(NQ == LQm2 .OR. NQ == LQi4) ASG(NQ,0:NRMAX) = ASG(NQ,0:NRMAX) * 1.D-2
!!$             IF(NQ == LQm3) ASG(NQ,0:NRMAX) = ASG(NQ,0:NRMAX) * 1.D-1

             L_NR:DO NR = 0, NRMAX
                ! Show results
                IF (NT == IDISP .AND. NR == NRMAX) THEN
                   IF (NQ == 1) THEN
                      WRITE(6,'((1X,A5,A3," =",I3))')'#####','IC',IC
                      WRITE(6,'((1X,A4," =",1PD9.2),2X,A2," =",I3)') &
                           & 'EPS   ', EPS,'NRMAX   ', NRMAX
                   END IF
                   ! Maximum relative error over the profile
                   ERR1 = MAX(ERR1, ASG(NQ,NR))
                   WRITE(6,'((1X,A2," =",I2,2X,A2," =",1PD9.2, &
                        & 2X,A5," =",1PD9.2,A1,I2,2X,A5," =",1PD9.2))') &
                        & 'NQ    ', NQ , &
                        & 'AV    ', AV    ,  'AVMAX ', AVM   ,':',NRAVM, &
                        & 'SCMAX ', ERR1
                   IDISP = IDIV + NT
                   IF (NQ == NQMAX) CYCLE L_IC
                ELSEIF (NT /= IDISP .AND. ASG(NQ,NR) > EPS) THEN
                   IF(IDIAG == 3) THEN
!                      IASG(1:2) = MAXLOC(ASG(1:NQMAX,0:NRMAX))
                      IASG(1:2) = MAXLOC(ASG(1:NQ,0:NR))
                      WRITE(6,'(3I4,1P4E17.8)') IC,IASG(1),IASG(2)-1,XP(IASG(1),IASG(2)-1), &
                           &                  XN(IASG(1),IASG(2)-1),ASG(IASG(1),IASG(2)-1),EPS
                   END IF
                   CYCLE L_IC
                END IF
             END DO L_NR

          END DO L_NQ

          EXIT L_IC
       END DO L_IC

       ! Save past X for BDF
       IF(IGBDF /= 0) XOLD(1:NQMAX,0:NRMAX) = X(1:NQMAX,0:NRMAX)

       IF(IDIAG >= 2) THEN
          IASG(1:2) = MAXLOC(ASG(1:NQMAX,0:NRMAX))
          IF(IC-1 == ICMAX) THEN
             WRITE(6,'(3I4,1P4E17.8,A2)') IC,IASG(1),IASG(2)-1,XP(IASG(1),IASG(2)-1), &
                  &                     XN(IASG(1),IASG(2)-1),ASG(IASG(1),IASG(2)-1),EPS," *"
          ELSE
             WRITE(6,'(3I4,1P4E17.8)') IC,IASG(1),IASG(2)-1,XP(IASG(1),IASG(2)-1), &
                  &                  XN(IASG(1),IASG(2)-1),ASG(IASG(1),IASG(2)-1),EPS
          END IF
       END IF

       ! Calculation fully converged
       X(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)

       ! Calculate mesh and coefficients at the next step
       CALL TXCALV(X)
       CALL TXCALC
!       write(6,'(A4,I4,A1,2X,E20.10)') 'TMP(',NT,')',PTeV(NRA)
!       write(6,'(A4,I4,A1,2X,E20.10)') 'TMP(',NT,')',PTiV(NRA)

       IF(IDIAG == 0 .OR. IDIAG == 2) THEN
          IF ((MOD(NT, NTSTEP) == 0) .AND. (NT /= NTMAX)) &
               & WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC
       ELSE IF(IDIAG > 0) THEN
          IF(NT /= NTMAX) THEN
             IF(IC-1 == ICMAX) THEN
                WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3,"  *")') NT,T_TX,IC
             ELSE
                WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC
             END IF
!             write(6,*) AJT,AJOHT
          END IF
       END IF

180    IF (MOD(NT, NGRSTP) == 0) CALL TXSTGR(NGR,GT,GY,NRMAX,NGRM,NGYRM)

       IF (MOD(NT, NGTSTP) == 0) THEN
          CALL TXGLOB
          CALL TXSTGT(SNGL(T_TX))
          call txstgq !!!temporary
          if(IDIAG < 0) call steady_check
       END IF

       IF (MOD(NT, NGVSTP) == 0) CALL TXSTGV(SNGL(T_TX))

       IF (MOD(NT, NTMAX ) == 0) CALL TXSTGQ

       call cal_flux

       IF (IERR /= 0) EXIT L_NTDO

    END DO L_NTDO

    rIPs = rIPe

    WRITE(6,'(1x ,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC

    deallocate(BA,BL,BX)
    IF(MDSOLV == 1) deallocate(XNvec,FL,FLP,DltXN,DltXP,BAE)

    RETURN
  END SUBROUTINE TXLOOP

!*********************************************************************************
!
!   Calculate coefficients matrix BA or BL and vector BX for linearlized equation
!
!      BA : coefficient matrix for BANDRD solver
!      BL : coefficient matrix for LAPACK DGBSV solver
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

    use tx_commons, only : IGBDF, ADV, MDLPCK, NQMAX, NRMAX, NLCR, CLC, BLC, ALC, NLCMAX, &
         &              PLC, X, XOLD!, lqb1, pnbv
    real(8), dimension(:,:), intent(inout) :: BA, BL
    real(8), dimension(:), intent(inout) :: BX
    INTEGER(4) :: I, J, NR, NQ, NC, NC1, IA, IB, IC
    INTEGER(4) :: JA, JB, JC, KL
    REAL(8) :: C43 = 4.D0/3.D0, C23 = 2.D0/3.D0, C13 = 1.D0/3.D0, COEF1, COEF2, COEF3!, suml1,suml2
    
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
       BA(1:NQMAX*4-1,1:NQMAX*(NRMAX+1)) = 0.D0

       ! Time derivative (NC=0)
       DO NR = 0, NRMAX
          DO NQ = 1, NQMAX
             NC = 0
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 == 0) CYCLE
             IC = NQMAX + (NC1 - 1) - (NQ - 1)
             IB = IC + NQMAX
             IA = IB + NQMAX
             J = NR * NQMAX + NQ
             BA(IC,J) = BA(IC,J) + CLC(NC,NQ,NR)
             BA(IB,J) = BA(IB,J) + BLC(NC,NQ,NR)
             BA(IA,J) = BA(IA,J) + ALC(NC,NQ,NR)
          END DO
       END DO

       ! *** NC /= 0 ***
       DO NR = 0, NRMAX
          DO NQ = 1, NQMAX
             DO NC = 1, NLCMAX(NQ)
                NC1 = NLCR(NC,NQ,NR)
                IF(NC1 == 0) CYCLE
                IC = NQMAX + (NC1 - 1) - (NQ - 1)
                IB = IC + NQMAX
                IA = IB + NQMAX
                J = NR * NQMAX + NQ
                BA(IC,J) = BA(IC,J) - CLC(NC,NQ,NR) * ADV
                BA(IB,J) = BA(IB,J) - BLC(NC,NQ,NR) * ADV
                BA(IA,J) = BA(IA,J) - ALC(NC,NQ,NR) * ADV
             END DO
          END DO
       END DO

    ! For LAPACK DGBSV solver
    ELSE
       BL(1:6*NQMAX-2,1:NQMAX*(NRMAX+1)) = 0.D0
       KL = 2 * NQMAX - 1

       ! Time derivative (NC=0)
       NC = 0
       DO NR = 0, NRMAX
          DO NQ = 1, NQMAX
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
                IB = IA + NQMAX
                IC = IB + NQMAX
                JA =(NR + 1) * NQMAX + NC1
                JB = NR      * NQMAX + NC1
                JC =(NR - 1) * NQMAX + NC1
                IF(NR == 0) THEN
                   BL(IA,JA) = BL(IA,JA) + ALC(NC,NQ,NR)
                   BL(IB,JB) = BL(IB,JB) + BLC(NC,NQ,NR)
                ELSEIF(NR == NRMAX) THEN
                   BL(IB,JB) = BL(IB,JB) + BLC(NC,NQ,NR)
                   BL(IC,JC) = BL(IC,JC) + CLC(NC,NQ,NR)
                ELSE
                   BL(IA,JA) = BL(IA,JA) + ALC(NC,NQ,NR)
                   BL(IB,JB) = BL(IB,JB) + BLC(NC,NQ,NR)
                   BL(IC,JC) = BL(IC,JC) + CLC(NC,NQ,NR)
                END IF
             END IF
          END DO
       END DO

       ! *** NC /= 0 ***
       DO NR = 0, NRMAX
          DO NQ = 1, NQMAX
             DO NC = 1, NLCMAX(NQ)
                NC1 = NLCR(NC,NQ,NR)
                IF(NC1 /= 0) THEN
                   IA = NQMAX - (NC1 - 1) + (NQ - 1) + KL
                   IB = IA + NQMAX
                   IC = IB + NQMAX
                   JA =(NR + 1) * NQMAX + NC1
                   JB = NR      * NQMAX + NC1
                   JC =(NR - 1) * NQMAX + NC1
                   IF(NR == 0) THEN
                      BL(IA,JA) = BL(IA,JA) - ALC(NC,NQ,NR) * ADV
                      BL(IB,JB) = BL(IB,JB) - BLC(NC,NQ,NR) * ADV
                   ELSEIF(NR == NRMAX) THEN
                      BL(IB,JB) = BL(IB,JB) - BLC(NC,NQ,NR) * ADV
                      BL(IC,JC) = BL(IC,JC) - CLC(NC,NQ,NR) * ADV
                   ELSE
                      BL(IA,JA) = BL(IA,JA) - ALC(NC,NQ,NR) * ADV
                      BL(IB,JB) = BL(IB,JB) - BLC(NC,NQ,NR) * ADV
                      BL(IC,JC) = BL(IC,JC) - CLC(NC,NQ,NR) * ADV
                   END IF
                END IF
             END DO
          END DO
       END DO

!!$       ! Ratio between stochastic loss and other losses
!!$       suml1 = 0.d0
!!$       do nr=0,26
!!$          suml1 = suml1 + BLC(7,lqb1,nr) * PNbV(nr)
!!$       end do
!!$       suml2 = 0.d0
!!$       do nr=27,nrmax
!!$          suml2 = suml2 + BLC(7,lqb1,nr) * PNbV(nr)
!!$       end do
!!$       write(6,*) suml1,suml2

    END IF

    ! *** Right-hand-side vector, denoted by "bu" ***

    BX(1:NQMAX*(NRMAX+1)) = 0.D0

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
    NC = 0
    NR = 0
    DO NQ = 1, NQMAX
       NC1 = NLCR(NC,NQ,NR)
       IF(NC1 /= 0) THEN
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) + BLC(NC,NQ,NR) * X(NC1,NR  ) * COEF1 &
               &                       + ALC(NC,NQ,NR) * X(NC1,NR+1) * COEF1 &
               &                       - BLC(NC,NQ,NR) * XOLD(NC1,NR  ) * COEF2 &
               &                       - ALC(NC,NQ,NR) * XOLD(NC1,NR+1) * COEF2
       END IF
    END DO

    DO NR = 1, NRMAX - 1
       DO NQ = 1, NQMAX
          NC1 = NLCR(NC,NQ,NR)
          IF(NC1 /= 0) THEN
             BX(NQMAX * NR + NQ) &
                  &    = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1) * COEF1 &
                  &                          + BLC(NC,NQ,NR) * X(NC1,NR  ) * COEF1 &
                  &                          + ALC(NC,NQ,NR) * X(NC1,NR+1) * COEF1 &
                  &                          - CLC(NC,NQ,NR) * XOLD(NC1,NR-1) * COEF2 &
                  &                          - BLC(NC,NQ,NR) * XOLD(NC1,NR  ) * COEF2 &
                  &                          - ALC(NC,NQ,NR) * XOLD(NC1,NR+1) * COEF2
          END IF
       END DO
    END DO

    NR = NRMAX
    DO NQ = 1, NQMAX
       NC1 = NLCR(NC,NQ,NR)
       IF(NC1 /= 0) THEN
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1) * COEF1 &
               &                       + BLC(NC,NQ,NR) * X(NC1,NR  ) * COEF1 &
               &                       - CLC(NC,NQ,NR) * XOLD(NC1,NR-1) * COEF2 &
               &                       - BLC(NC,NQ,NR) * XOLD(NC1,NR  ) * COEF2
       END IF
    END DO

    ! *** In the case of general term effect in any equations, denoted by "c" ***
    DO NR = 0, NRMAX
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) + PLC(NC,NQ,NR) * COEF3
          END DO
       END DO
    END DO

    ! *** Only the case of not BDF ***
    !  Because "ADV" has no longer original meaning when BDF is used.

    IF(IGBDF == 0) THEN
    NR = 0
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                     &              +(  BLC(NC,NQ,NR) * X(NC1,NR  ) &
                     &                + ALC(NC,NQ,NR) * X(NC1,NR+1)) * (1.D0 - ADV)
             END IF
          END DO
       END DO

    DO NR = 1, NRMAX-1
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                     &              +(  CLC(NC,NQ,NR) * X(NC1,NR-1) &
                     &                + BLC(NC,NQ,NR) * X(NC1,NR  ) &
                     &                + ALC(NC,NQ,NR) * X(NC1,NR+1)) * (1.D0 - ADV)
             END IF
          END DO
       END DO
    END DO

    NR = NRMAX
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                     &              +(  CLC(NC,NQ,NR) * X(NC1,NR-1) &
                     &                + BLC(NC,NQ,NR) * X(NC1,NR  )) * (1.D0 - ADV)
             END IF
          END DO
       END DO
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
    REAL(8), DIMENSION(1:NQMAX,0:NRMAX), intent(in) :: XL
    integer(4) :: NR

    IER = 0

    DO NR = 0, NRMAX
       IF (XL(LQe1,NR) < 0.D0 .OR. XL(LQi1,NR) < 0.D0 .OR. &
        &  XL(LQb1,NR) < 0.D0 .OR. XL(LQr1,NR) < 0.D0) THEN
          WRITE(6,'(2A,I4,2(A,I4),A)') '### ERROR(TXLOOP) : Negative density at ', &
               &           'NR =', NR, ', NT=', NTL, ', IC=', IC, '.'
          WRITE(6,'(1X,4(A,1PE10.3))') 'ne =', SNGL(XL(LQe1,NR)), &
               &                    ',   ni =', SNGL(XL(LQi1,NR)), &
               &                    ',   nb =', SNGL(XL(LQb1,NR)), &
               &                    ', nbrp =', SNGL(XL(LQr1,NR))
          IER = 1
          RETURN
       END IF
    END DO

    DO NR = 0, NRMAX
       IF (XL(LQe5,NR) < 0.D0 .OR. XL(LQi5,NR) < 0.D0) THEN
          WRITE(6,'(2A,I4,2(A,I4),A)') '### ERROR(TXLOOP) : Negative temperature at ', &
               &           'NR =', NR, ', NT=', NTL, ', IC=', IC, '.'
          WRITE(6,'(20X,2(A,1PE15.7))') 'Te =', SNGL(XL(LQe5,NR)), &
               &           ',   Ti =', SNGL(XL(LQi5,NR))
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

  SUBROUTINE MINUS_GOES_ZERO(XL,LQ,ID)
    
    use tx_commons, only : NQMAX, NRMAX, NRA
    integer(4), intent(in) :: LQ, ID
    real(8), dimension(1:NQMAX,0:NRMAX), intent(inout) :: XL
    integer(4) :: NR, NZERO

    IF(ID == 0) THEN
       IF(MINVAL(XL(LQ,0:NRMAX)) < 0.D0) THEN
          DO NR = 0, NRMAX
             IF(XL(LQ,NR) <= 0.D0) THEN
                IF(NR < NRA) THEN
                   XL(LQ,NR) = 0.D0
                ELSE
                   NZERO = NR
                   EXIT
                END IF
             END IF
          END DO

          XL(LQ,NZERO:NRMAX) = 0.D0
       END IF
    ELSE
       IF(MINVAL(XL(LQ,0:NRMAX)) < 0.D0) THEN
          DO NR = NRMAX, 0, -1
             IF(XL(LQ,NR) <= 0.D0) THEN
                NZERO = NR
                EXIT
             END IF
          END DO

          XL(LQ,0:NZERO) = 0.D0
       END IF
    END IF

  END SUBROUTINE MINUS_GOES_ZERO

!!$  SUBROUTINE THRESHOLD(XL,ID)
!!$
!!$    real(8), dimension(0:NRMAX), intent(inout) :: XL
!!$    integer(4), intent(out) :: ID
!!$    integer(4) :: NR, NRL
!!$
!!$    ID = 0
!!$    NRL = 0.5 * NRA
!!$    IF(MINVAL(XL(0:NRL)) < 1.D-8 .AND. MAXVAL(XL(0:NRL)) > 0.D0) THEN
!!$       XL(0:NRMAX) = 0.D0
!!$       ID = 1
!!$    END IF
!!$
!!$  END SUBROUTINE THRESHOLD

!*********************************************************************************
!
!   Calculate coefficients matrix BA or BL and vector BX for linearlized equation
!
!      BA : coefficient matrix for BANDRD solver
!      BL : coefficient matrix for LAPACK DGBSV solver
!      BX : right-hand-side vector
!
!   ** Simple explanation **
!
!      A(u)u=b(u)+c, where A denotes the matrix, u the vector to be solved
!                          b the coefficients vector of u,
!                          c the additive coefficients vector independent on u
!
!*********************************************************************************

  SUBROUTINE TXCALB_SECANT(FL,BX,XNvec)

    use tx_commons, only : IGBDF, ADV, NQMAX, NRMAX, NLCR, CLC, BLC, ALC, NLCMAX, &
         &              PLC, X, XOLD
    real(8), dimension(:), intent(out) :: FL
    real(8), dimension(:), intent(inout) :: BX
    real(8), dimension(:), intent(in) :: XNvec
    INTEGER(4) :: NR, NQ, NC, NC1, IA, IB, IC
    REAL(8) :: C43 = 4.D0/3.D0, C23 = 2.D0/3.D0, C13 = 1.D0/3.D0, COEF1, COEF2, COEF3
    
    IF(IGBDF /= 0) ADV = C23

    !  NR : number of radial mesh 
    !  NQ : number of equation

    ! *** Left-hand-side coefficient matrix, denoted by "A" ***

    !***************************************************************
    !   ALC, BLC and CLC are NQMAX x max(NLCMAX) on each NR,
    !   then IC, IB and IA are NQMAX away from each other,
    !   which are pointers of columns of the matrix BA in terms of 
    !   certain variable in certain term in certain equation.
    !***************************************************************

    FL(1:NQMAX*(NRMAX+1)) = 0.D0

    ! Time derivative (NC=0)
    DO NR = 0, NRMAX
       DO NQ = 1, NQMAX
          NC = 0
          NC1 = NLCR(NC,NQ,NR)
          IF(NC1 == 0) CYCLE
          IC = NQMAX * (NR - 1) + NC1
          IB = IC + NQMAX
          IA = IB + NQMAX
          IF(NR /= 0)     FL(IC) = FL(IC) + CLC(NC,NQ,NR) * XNvec(IC)
                          FL(IB) = FL(IB) + BLC(NC,NQ,NR) * XNvec(IB)
          IF(NR /= NRMAX) FL(IA) = FL(IA) + ALC(NC,NQ,NR) * XNvec(IA)
       END DO
    END DO

    ! *** NC /= 0 ***
    DO NR = 0, NRMAX
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 == 0) CYCLE
             IC = NQMAX * (NR - 1) + NC1
             IB = IC + NQMAX
             IA = IB + NQMAX
             IF(NR /= 0)     FL(IC) = FL(IC) - CLC(NC,NQ,NR) * ADV * XNvec(IC)
                             FL(IB) = FL(IB) - BLC(NC,NQ,NR) * ADV * XNvec(IB)
             IF(NR /= NRMAX) FL(IA) = FL(IA) - ALC(NC,NQ,NR) * ADV * XNvec(IA)
!             write(6,'(7I6)') NR,NQ,NC,NC1,IC,IB,IA
          END DO
       END DO
    END DO

    ! *** Right-hand-side vector, denoted by "bu" ***

    BX(1:NQMAX*(NRMAX+1)) = 0.D0

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
    NC = 0
    NR = 0
    DO NQ = 1, NQMAX
       NC1 = NLCR(NC,NQ,NR)
       IF(NC1 /= 0) THEN
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) + BLC(NC,NQ,NR) * X(NC1,NR  ) * COEF1 &
               &                       + ALC(NC,NQ,NR) * X(NC1,NR+1) * COEF1 &
               &                       - BLC(NC,NQ,NR) * XOLD(NC1,NR  ) * COEF2 &
               &                       - ALC(NC,NQ,NR) * XOLD(NC1,NR+1) * COEF2
       END IF
    END DO

    DO NR = 1, NRMAX - 1
       DO NQ = 1, NQMAX
          NC1 = NLCR(NC,NQ,NR)
          IF(NC1 /= 0) THEN
             BX(NQMAX * NR + NQ) &
                  &    = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1) * COEF1 &
                  &                          + BLC(NC,NQ,NR) * X(NC1,NR  ) * COEF1 &
                  &                          + ALC(NC,NQ,NR) * X(NC1,NR+1) * COEF1 &
                  &                          - CLC(NC,NQ,NR) * XOLD(NC1,NR-1) * COEF2 &
                  &                          - BLC(NC,NQ,NR) * XOLD(NC1,NR  ) * COEF2 &
                  &                          - ALC(NC,NQ,NR) * XOLD(NC1,NR+1) * COEF2
          END IF
       END DO
    END DO

    NR = NRMAX
    DO NQ = 1, NQMAX
       NC1 = NLCR(NC,NQ,NR)
       IF(NC1 /= 0) THEN
          BX(NQMAX * NR + NQ) &
               & = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1) * COEF1 &
               &                       + BLC(NC,NQ,NR) * X(NC1,NR  ) * COEF1 &
               &                       - CLC(NC,NQ,NR) * XOLD(NC1,NR-1) * COEF2 &
               &                       - BLC(NC,NQ,NR) * XOLD(NC1,NR  ) * COEF2
       END IF
    END DO

    ! *** In the case of general term effect in any equations, denoted by "c" ***
    DO NR = 0, NRMAX
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) + PLC(NC,NQ,NR) * COEF3
          END DO
       END DO
    END DO

    ! *** Only the case of not BDF ***
    !  Because "ADV" has no longer original meaning when BDF is used.

    IF(IGBDF == 0) THEN
    NR = 0
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                     &              +(  BLC(NC,NQ,NR) * X(NC1,NR  ) &
                     &                + ALC(NC,NQ,NR) * X(NC1,NR+1)) * (1.D0 - ADV)
             END IF
          END DO
       END DO

    DO NR = 1, NRMAX-1
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                     &              +(  CLC(NC,NQ,NR) * X(NC1,NR-1) &
                     &                + BLC(NC,NQ,NR) * X(NC1,NR  ) &
                     &                + ALC(NC,NQ,NR) * X(NC1,NR+1)) * (1.D0 - ADV)
             END IF
          END DO
       END DO
    END DO

    NR = NRMAX
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             NC1 = NLCR(NC,NQ,NR)
             IF(NC1 /= 0) THEN
                BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) &
                     &              +(  CLC(NC,NQ,NR) * X(NC1,NR-1) &
                     &                + BLC(NC,NQ,NR) * X(NC1,NR  )) * (1.D0 - ADV)
             END IF
          END DO
       END DO
    END IF

    ! Making LHS term of "f(x)=0"
    FL(1:NQMAX*(NRMAX+1)) = FL(1:NQMAX*(NRMAX+1)) - BX(1:NQMAX*(NRMAX+1))

    RETURN
  END SUBROUTINE TXCALB_SECANT

end module tx_main
