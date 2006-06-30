!     $Id$
module main
  use commons
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
    use libraries, only : APRTOS
    use results
    use output_console, only : TXWDAT

    INTEGER :: NDY, NDM, NDD, NTH, NTM, NTS, iRTIME1, iRTIME2, &
         &     iRTIME3, NSTR1, NSTR2, NSTR3
    REAL :: gCTIME1, gCTIME2, gCTIME3
    character(len=10) :: STR1, STR2, STR3
    INTEGER, DIMENSION(1:8) :: TIMES

    IF (IERR /= 0) THEN
       WRITE(6,*) '### ERROR(TXEXEC) : Error should be cleared.'
       RETURN
    END IF

    CALL DATE_AND_TIME(VALUES=TIMES)
    iRTIME1 = TIMES(5) * 60 * 60 + TIMES(6) * 60 + TIMES(7)
    CALL CPU_TIME(gCTIME1)

    ! ***** Core part of simulation *****

    WRITE(6,*) 'Calculating...'
    CALL TXLOOP

    ! ***********************************

    CALL DATE_AND_TIME(VALUES=TIMES)
    iRTIME2 = TIMES(5) * 60 * 60 + TIMES(6) * 60 + TIMES(7)
    iRTIME3 = iRTIME2 - iRTIME1
    if (iRTIME3 < 0) iRTIME3 = iRTIME3 + 24 * 60 * 60
    CALL CPU_TIME(gCTIME2)
    gCTIME3 = gCTIME2-gCTIME1

    NSTR1 = 0
    CALL APRTOS(STR1, NSTR1, gCTIME3, 'F2')
    IF (iRTIME3 == 0) THEN
       WRITE(6,*) 'real =', iRTIME3, '(sec)   ',  &
            &              'CPU = ', STR1(1:NSTR1), '(sec)'
    ELSE
       NSTR2 = 0
       CALL APRTOS(STR2, NSTR2, gCTIME3 / (iRTIME3 + 1) * 100, 'F1')
       NSTR3 = 0
       CALL APRTOS(STR3, NSTR3, gCTIME3 / (iRTIME3 + 0) * 100, 'F1')
       WRITE(6,*) 'real =', iRTIME3, '(sec)   ',  &
            &     'CPU = ', STR1(1:NSTR1), '(sec)   ',  &
            &     '(', STR2(1:NSTR2), '% - ', STR3(1:NSTR3), '%)'
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
    use results
    use variables
    use coefficients, only : TXCALA
    use graphic, only : TXSTGT, TXSTGV, TXSTGR

    INTEGER :: I, J, NR, NQ, NC, NC1, IA, IB, IC, IDIV, NTDO, IDISP, NRAVM
    REAL(8) :: TIME0, DIP, SUML, AVM, ERR1, AV, FCTR
    REAL(8), DIMENSION(NQM,0:NRM) :: XN, XP

    IF (MODEAV == 0) THEN
       IDIV = NTMAX + 1
    ELSE
       IDIV = NTMAX / MODEAV
    END IF
    TIME0 = T_TX
    DIP = (rIPe - rIPs) / NTMAX

    L_NTDO:DO NTDO = 1, NTMAX
       NT = NTDO
       T_TX = TIME0 + DT*NT
       rIP  = rIPs  + DIP*NT

       ! Create new X = XN
       XN(1:NQMAX,0:NRMAX) = X(1:NQMAX,0:NRMAX)

       L_IC : DO IC = 1, ICMAX
          ! Save past X = XP
          XP(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)

          CALL TXCALV(XP)
          CALL TXCALC
          CALL TXCALA
          CALL TXCALB
          CALL TXGLOB
          CALL BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQM-1, IERR)
          IF (IERR >= 30000) THEN
             WRITE(6,'(3(A,I6))') '### ERROR(TXLOOP) : Matrix BA is singular at ',  &
                  &              NT, ' -', IC, ' step. IERR=',IERR
             IERR = 1
             XN(1:NQMAX,0:NRMAX) = XP(1:NQMAX,0:NRMAX)
             GOTO 180
          END IF

          ! Copy calculated variables' vector to variable matrix
          DO NR = 0, NRMAX
             DO NQ = 1, NQMAX
                XN(NQ,NR) = BX(NQMAX * NR + NQ)
             END DO
          END DO
          ! Avoid negative value
          WHERE(XN(16,0:NRMAX) < 0.D0) XN(16,0:NRMAX) = 0.D0
          WHERE(XN(19,0:NRMAX) < 0.D0) XN(19,0:NRMAX) = 0.D0

          ! Check negative density or temperature in variable matrix
          CALL TXCHCK(NT,IC,XN,IERR)
          IF (IERR /= 0) THEN
             X(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)
             CALL TXCALV(X)
             CALL TXCALC
             CALL TXGLOB
             CALL TXSTGT(SNGL(T_TX))
             CALL TXSTGV(SNGL(T_TX))
             CALL TXSTGR
             RETURN
          END IF

!          IF(IC == 1) EXIT L_IC

          L_NQ:DO NQ = 1, NQMAX
             ! Calculate maximum root-mean-square X and grid point at that time
             AVM  = 0.D0
             DO NR = 0, NRMAX
                IF(SQRT(XN(NQ,NR)**2) > AVM) THEN
                   AVM = SQRT(XN(NQ,NR)**2)
                   NRAVM= NR
                END IF
             END DO

             ! Calculate root-mean-square X over the profile
             ERR1  = 0.D0
             IDISP = IDIV
             SUML  = SUM(XN(NQ,0:NRMAX)**2)
             AV    = SQRT(SUML) / NRMAX
             IF (AV == 0.D0) CYCLE L_NQ
             L_NR:DO NR = 0, NRMAX
                ! Maximum relative error over the profile
                ERR1 = MAX(ERR1, ABS(XN(NQ,NR) - XP(NQ,NR)) / AV)
                ! Show results
                IF (NT == IDISP .AND. NR == NRMAX) THEN
                   IF (NQ == 1) THEN
                      WRITE(6,'((1X,A5,A3," =",I3))')'#####','IC',IC
                      WRITE(6,'((1X,A4," =",1PD9.2),2X,A2," =",I3)') &
                           & 'EPS   ', EPS,'NRMAX   ', NRMAX
                   END IF
                   WRITE(6,'((1X,A2," =",I2,2X,A2," =",1PD9.2, &
                        & 2X,A5," =",1PD9.2,A1,I2,2X,A5," =",1PD9.2))') &
                        & 'NQ    ', NQ , &
                        & 'AV    ', AV    ,  'AVMAX ', AVM   ,':',NRAVM, &
                        & 'SCMAX ', ERR1
                   IDISP = IDIV + NT
                   IF (NQ == NQMAX) CYCLE L_IC
                ELSEIF (NT /= IDISP .AND.  &
                     & ABS(XN(NQ,NR) - XP(NQ,NR)) / AV > EPS) THEN
                   CYCLE L_IC
                END IF
             END DO L_NR

          END DO L_NQ
          EXIT L_IC
       END DO L_IC

       ! Calculation fully converged
       X(1:NQMAX,0:NRMAX) = XN(1:NQMAX,0:NRMAX)

       ! Calculate mesh and coefficients at the next step
       CALL TXCALV(X)
       CALL TXCALC

       IF ((MOD(NT, NTSTEP) == 0) .AND. (NT /= NTMAX)) &
            & WRITE(6,'(1x,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC

180    IF (MOD(NT, NGRSTP) == 0) CALL TXSTGR

       IF (MOD(NT, NGTSTP) == 0) THEN
          CALL TXGLOB
          CALL TXSTGT(SNGL(T_TX))
       END IF

       IF (MOD(NT, NGVSTP) == 0) CALL TXSTGV(SNGL(T_TX))

       IF (IERR /= 0) EXIT L_NTDO

    END DO L_NTDO

    rIPs = rIPe

    WRITE(6,'(1x ,"NT =",I4,"   T =",1PD9.2,"   IC =",I3)') NT,T_TX,IC

    RETURN
  END SUBROUTINE TXLOOP

!***************************************************************
!
!   Calculate BA, BX
!      BA : coefficient matrix
!      BX : right-hand-side vector
!
!***************************************************************

  SUBROUTINE TXCALB

    INTEGER :: I, J, NR, NQ, NC, NC1, IA, IB, IC

    !  NR : number of radial mesh 
    !  NQ : number of equation

    ! *** Left-hand-side coefficient matrix ***

    BA(1:NQMAX*4-1,1:NQMAX*(NRMAX+1)) = 0.D0

    !***************************************************************
    !   ALC, BLC and CLC are NQMAX x max(NLCMAX) on each NR,
    !   then IC, IB and IA are NQMAX away from each other,
    !   which are pointers of columns of the matrix BA in terms of 
    !   certain variable in certain term in certain equation.
    !***************************************************************

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
             BA(IC,J) = BA(IC,J) - CLC(NC,NQ,NR)
             BA(IB,J) = BA(IB,J) - BLC(NC,NQ,NR)
             BA(IA,J) = BA(IA,J) - ALC(NC,NQ,NR)
          END DO
       END DO
    END DO

!    BA(40,4) =-1.D0 / 0.15D0!6.6667D0

    ! *** Right-hand-side vector ***

    BX(1:NQMAX*(NRMAX+1)) = 0.D0

    ! In the case of NC=0 i.e. time derivative term effects in any equations
    NC = 0
    NR = 0
    DO NQ = 1, NQMAX
       NC1 = NLCR(NC,NQ,NR)
       BX(NQMAX * NR + NQ) &
            & = BX(NQMAX * NR + NQ) + BLC(NC,NQ,NR) * X(NC1,NR  ) &
            &                       + ALC(NC,NQ,NR) * X(NC1,NR+1)
    END DO

    DO NR = 1, NRMAX - 1
       DO NQ = 1, NQMAX
          NC1 = NLCR(NC,NQ,NR)
          BX(NQMAX * NR + NQ) &
               &    = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1) &
               &                          + BLC(NC,NQ,NR) * X(NC1,NR  ) &
               &                          + ALC(NC,NQ,NR) * X(NC1,NR+1)
       END DO
    END DO

    NR = NRMAX
    DO NQ = 1, NQMAX
       NC1 = NLCR(NC,NQ,NR)
       BX(NQMAX * NR + NQ) &
            & = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1) &
            &                       + BLC(NC,NQ,NR) * X(NC1,NR  )
    END DO

    ! In the case of general term effect in any equations
    DO NR = 0, NRMAX
       DO NQ = 1, NQMAX
          DO NC = 1, NLCMAX(NQ)
             BX(NQMAX * NR + NQ) = BX(NQMAX * NR + NQ) + PLC(NC,NQ,NR)
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE TXCALB

!***************************************************************
!
!   Check negative density or temperature
!
!***************************************************************

  SUBROUTINE TXCHCK(NTL,IC,XL,IER)

    INTEGER :: NTL, IC, IER, NR
    REAL(8), DIMENSION(NQM,0:NRM) :: XL

    IER = 0

    DO NR = 0, NRMAX
       IF (XL(LQe1,NR) < 0.D0 .OR. XL(LQi1,NR) < 0.D0) THEN
          WRITE(6,*) '### ERROR(TXLOOP) : ',  &
               &           'Density becomes negative at ',  &
               &           'NR =', NR, ', ', NTL, ' -', IC, ' step.'
          WRITE(6,*) 'ne =', SNGL(XL(LQe1,NR)), '   ni =', SNGL(XL(LQi1,NR))
          IER = 1
          RETURN
       END IF
    END DO

    DO NR = 0, NRMAX
       IF (XL(LQe5,NR) < 0.D0 .OR. XL(LQi5,NR) < 0.D0) THEN
          WRITE(6,*) '### ERROR(TXLOOP) : ',  &
               &           'Temperature becomes negative at ',  &
               &           'NR =', NR, ', ', NTL, ' -', IC, ' step.'
          WRITE(6,*) 'Te =', SNGL(XL(LQe5,NR)), '   Ti =', SNGL(XL(LQi5,NR))
          IER = 1
          RETURN
       END IF
    END DO

    RETURN
  END SUBROUTINE TXCHCK
end module main

