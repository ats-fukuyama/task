!     $Id$
!***************************************************************
!
!   Write Data
!
!***************************************************************

module output_console
  use commons
  implicit none
  public

contains

  SUBROUTINE TXWDAT
    use physical_constants, only : PI
    use libraries, only : INTG_F

    INTEGER :: NR
    REAL(8) :: rNbar

    !     ***** Volume-averaged density *****

    rNbar = 2.D0 * PI * INTG_F(PNeV) * 1.D20
    rNbar = rNbar / (PI * RB**2)

    WRITE(6,'((1X,A," =",1PD9.2,3(2X,A,"=",1PD9.2)))') &
         &     'Ne(0)',    PNeV(0),  &
         &     'UePhi(0)', (9*X(LQe4,0)-X(LQe4,1))/(8*PNeV(0)) / 1.D3,  &
         &     'UiPhi(0)', (9*X(LQi4,0)-X(LQi4,1))/(8*PNiV(0)) / 1.D3,  &
         &     'N0(RB)',   (9*X(LQn1,NRMAX-1)-X(LQn1,NRMAX-2))/8 * 1.D20,&
         &     'NB(0)',    rLINEAVE(0.D0)   / 1.D20,  &
         &     'NB(0.24)', rLINEAVE(0.24D0) / 1.D20,  &
         &     ' NB(0.6)', rLINEAVE(0.6D0)  / 1.D20,  &
         &     '    PF',   PNeV(0) * 1.D20 / rNbar
    WRITE(6,'(1X,A," =",1PD9.2,2X,A,"=",1PD9.2)') &
         &     'Te(0)',    PTeV(0),  &
         &     'Ti(0)   ', PTiV(0)
    RETURN
  END SUBROUTINE TXWDAT

!***************************************************************
!
!   Write Data 2
!
!***************************************************************

  SUBROUTINE TXWDAT2

    REAL :: gPNeMIN, gPNeMAX, gNB0MIN, gNB0MAX, gUiphMIN, gUiphMAX
    CALL GMNMX1(GTY(0,1),  1, NGT + 1, 1,  gPNeMIN,  gPNeMAX)
    CALL GMNMX1(GTY(0,2),  1, NGT + 1, 1,  gNB0MIN,  gNB0MAX)
    CALL GMNMX1(GTY(0,11), 1, NGT + 1, 1, gUiphMIN, gUiphMAX)
    WRITE(6,'((1X,A,"=",1PD9.2,2(2X,A,"=",1PD9.2)))') &
         &     'MAX(Ne(0))',     gPNeMAX / 1.E20,  &
         &     'MAX(NB(0))',     gNB0MAX / 1.E20,  &
         &     'MAX(UiPhi(0))', gUiphMAX / 1.E3, &
         &     'MIN(Ne(0))',     gPNeMIN / 1.E20,  &
         &     'MIN(NB(0))',     gNB0MIN / 1.E20,  &
         &     'MIN(UiPhi(0))', gUiphMIN / 1.E3

    RETURN
  END SUBROUTINE TXWDAT2

!***********************************************************
!
!  LINE AVERAGE OF rN
!
!***********************************************************

  REAL(8) FUNCTION rLINEAVE(Rho)

    INTEGER :: I, IR, NY = 100
    REAL(8), INTENT(IN) :: Rho
    REAL(8) :: D, DY, Y, RL, SUML

    SUML = 0.D0
    D = Rho * RA
    DY = SQRT(RA*RA - D*D) / NY
    DO I = 0, NY
       Y = DY * I
       RL = SQRT(Y*Y + D*D)
       IR = NINT(RL * NRMAX / RA)
       SUML = SUML + PNeV(IR) * 1.D20 * DY
    END DO
    rLINEAVE = SUML / SQRT(RA**2 - D**2)

    RETURN
  END FUNCTION rLINEAVE
end module output_console

module file_io
  use commons
  implicit none

contains

!***************************************************************
!
!   Save transport data
!
!***************************************************************

  SUBROUTINE TXSAVE
    use libraries, only : TOUPPER

    INTEGER :: IST, NQ, NR, I, IGYT, IGYV
    character(len=100) :: TXFNAM, RCSId
    character(len=1) :: STR
    LOGICAL :: LEX

    RCSId = ' '

    DO 
       WRITE(6,*) '# INPUT : SAVE FILE NAME'
       CALL GUFLSH
       READ(*,'(A100)',IOSTAT=IST) TXFNAM
       IF (IST > 0) CYCLE
       INQUIRE(FILE=TXFNAM,EXIST=LEX)
       IF (LEX) THEN
          WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',  &
               &              'ARE YOU SURE {Y/N} ?'
          CALL GUFLSH
          READ(*,'(A1)') STR
          CALL TOUPPER(STR)
          IF (STR == 'Y') THEN
             OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',FORM='UNFORMATTED')
             IF (IST == 0) THEN
                WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR OUTPUT.'
                EXIT
             ELSEIF (IST > 0) THEN
                WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
             END IF
          END IF
       ELSE
          OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='NEW',FORM='UNFORMATTED')
          IF (IST == 0) THEN
             WRITE(6,*) '# NEW FILE (', TXFNAM, ') IS CREATED FOR OUTPUT.'
             EXIT
          ELSEIF (IST > 0) THEN
             WRITE(6,*) 'XX  NEW FILE OPEN ERROR !, IOSTAT = ', IST
          END IF
       END IF
    END DO

    WRITE(21) SLID
    WRITE(21) RCSId

    WRITE(21) RA,RB,RC,RR,BB
    WRITE(21) PA,PZ,Zeff
    WRITE(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
    WRITE(21) De0,Di0,rMue0,rMui0,WPM0
    WRITE(21) Chie0,Chii0
    WRITE(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
    WRITE(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
    WRITE(21) rLn,rLT
    WRITE(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
    WRITE(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
    WRITE(21) DLT,DT,EPS
    WRITE(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
    WRITE(21) rG1
    WRITE(21) rIPs,rIPe

    WRITE(21) T_TX,TMAX,NT,NQMAX,IERR
    WRITE(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

    WRITE(21) NGT,NGYTM,NGYVM
    WRITE(21) (GTX(I), I=0, NGT)
    WRITE(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
    WRITE(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
    CLOSE(21)
    WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

    RETURN
  END SUBROUTINE TXSAVE

!***************************************************************
!
!   Load transport data
!
!***************************************************************

  SUBROUTINE TXLOAD
    use results
    use output_console
    use variables
    use init_prof
    use coefficients, only : TXCALA

    INTEGER :: IST, NQ, NR, NGYT, NGYV, I, IGYT, IGYV
    character(len=100) ::  TXFNAM, RCSId
    character(len=8) :: LOADSLID
    LOGICAL :: LEX

    ! tmp : NGYT

    DO 
       WRITE(6,*) '# INPUT : LOAD FILE NAME'
       CALL GUFLSH
       READ(*,'(A100)',IOSTAT=IST) TXFNAM
       IF (IST > 0) CYCLE
       INQUIRE(FILE=TXFNAM,EXIST=LEX)
       IF (LEX) THEN
          OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',FORM='UNFORMATTED')
          IF (IST == 0) THEN
             WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR INPUT.'
             EXIT
          ELSEIF (IST > 0) THEN
             WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
          END IF
       ELSE
          WRITE(6,*) 'XX  FILE (', TXFNAM, ') DOES NOT EXIST !'
       END IF
    END DO

    READ(21,IOSTAT=IST) LOADSLID
    IF (IST > 0) THEN
       WRITE(6,*) 'XX READ ERROR in TXLOAD !'
       CLOSE(21)
       RETURN
    END IF
    !  IF(LOADSLID(1:5) == 'tx300') THEN
    READ(21) RCSId

    READ(21) RA,RB,RC,RR,BB
    READ(21) PA,PZ,Zeff
    READ(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
    READ(21) De0,Di0,rMue0,rMui0,WPM0
    READ(21) Chie0,Chii0
    READ(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
    READ(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
    READ(21) rLn,rLT
    READ(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
    READ(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
    READ(21) DLT,DT,EPS
    READ(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
    READ(21) rG1
    READ(21) rIPs,rIPe

    READ(21) T_TX,TMAX,NT,NQMAX,IERR
    READ(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

    READ(21) NGT,NGYT,NGYV
    READ(21) (GTX(I), I=0, NGT)
    READ(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYT)
    READ(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYV)
    !  END IF
    CLOSE(21)
    WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

    IF (NQMAX == 16) THEN
       X(LQm2,0:NRMAX) = 0
       X(LQm5,0:NRMAX) = BB
       X(LQb3,0:NRMAX) = 0
       X(LQn2,0:NRMAX) = 0
    END IF

    NGT=-1
    NGR=-1
    NGVV=-1
    rIP=rIPs

    CALL TXCALM
    CALL TXCALV(X)
    CALL TXCALC
    CALL TXCALA
    CALL TXGLOB
    CALL TXWDAT
    CALL TXWDAT2

    RETURN
  END SUBROUTINE TXLOAD

!***************************************************************
!
!   Save graphic data
!
!***************************************************************

  SUBROUTINE TXGSAV

    INTEGER :: IST, NQ, NR, IGR, I, IGYR, IGYT, IGYV
    character(len=100) :: TXFNAM, RCSId
    character(len=1) :: STR
    LOGICAL :: LEX

    RCSId = ' '

    DO 
       WRITE(6,*) '# INPUT : SAVE FILE NAME'
       CALL GUFLSH
       READ(*,'(A100)',IOSTAT=IST) TXFNAM
       IF (IST > 0) CYCLE
       INQUIRE(FILE=TXFNAM,EXIST=LEX)
       IF (LEX) THEN
          WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',  &
               &              'ARE YOU SURE {Y/N} ?'
          CALL GUFLSH
          READ(*,'(A1)') STR
          IF (STR /= 'Y' .AND. STR /= 'y') THEN
          ELSE
             OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',FORM='UNFORMATTED')
             IF (IST == 0) THEN
                WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR OUTPUT.'
                EXIT
             ELSEIF (IST > 0) THEN
                WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
             END IF
          END IF
       ELSE
          OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='NEW',FORM='UNFORMATTED')
          IF (IST == 0) THEN
             WRITE(6,*) '# NEW FILE (', TXFNAM, ') IS CREATED FOR OUTPUT.'
             EXIT
          ELSEIF (IST > 0) THEN
             WRITE(6,*) 'XX  NEW FILE OPEN ERROR !, IOSTAT = ', IST
          END IF
       END IF
    END DO

    WRITE(21) SLID
    WRITE(21) RCSId

    WRITE(21) RA,RB,RC,RR,BB
    WRITE(21) PA,PZ,Zeff
    WRITE(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
    WRITE(21) De0,Di0,rMue0,rMui0,WPM0
    WRITE(21) Chie0,Chii0
    WRITE(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
    WRITE(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
    WRITE(21) rLn,rLT
    WRITE(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
    WRITE(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
    WRITE(21) DLT,DT,EPS
    WRITE(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
    WRITE(21) rG1
    WRITE(21) rIPs,rIPe

    WRITE(21) T_TX,TMAX,NT,NQMAX,IERR
    WRITE(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

    WRITE(21) NGR,NGYRM
    WRITE(21) NGT,NGYTM
    WRITE(21) NGVV,NGYVM
    WRITE(21) (GT(IGR), IGR=0, NGR)
    WRITE(21)(((GY(I,IGR,IGYR), I=0,NRMAX), IGR=0,NGR), IGYR=1,NGYRM)
    WRITE(21) (GTX(I), I=0, NGT)
    WRITE(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
    WRITE(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
    CLOSE(21)
    WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

    RETURN
  END SUBROUTINE TXGSAV

!***************************************************************
!
!   Load graphic data
!
!***************************************************************

  SUBROUTINE TXGLOD
    use results
    use output_console
    use init_prof
    use variables

    INTEGER :: IST, NQ, NR, NGYR, NGYT, NGYV, IGR, I, IGYR, IGYT, IGYV
    character(len=100) :: TXFNAM, RCSId
    character(len=8) :: LOADSLID
    LOGICAL :: LEX

    ! tmp : NGYT

    DO 
       WRITE(6,*) '# INPUT : LOAD FILE NAME'
       CALL GUFLSH
       READ(*,'(A100)',IOSTAT=IST) TXFNAM
       IF (IST > 0) CYCLE
       INQUIRE(FILE=TXFNAM,EXIST=LEX)
       IF (LEX) THEN
          OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',FORM='UNFORMATTED')
          IF (IST == 0) THEN
             WRITE(6,*) '# OLD FILE (', TXFNAM, ') IS ASSIGNED FOR INPUT.'
             EXIT
          ELSEIF (IST > 0) THEN
             WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
          END IF
       ELSE
          WRITE(6,*) 'XX  FILE (', TXFNAM, ') DOES NOT EXIST !'
       END IF
    END DO

    READ(21,IOSTAT=IST) LOADSLID
    IF (IST > 0) THEN
       WRITE(6,*) 'XX READ ERROR in TXGLOD !'
       CLOSE(21)
       RETURN
    END IF
    !  IF(LOADSLID(1:5) == 'tx300') THEN
    READ(21) RCSId

    READ(21) RA,RB,RC,RR,BB
    READ(21) PA,PZ,Zeff
    READ(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
    READ(21) De0,Di0,rMue0,rMui0,WPM0
    READ(21) Chie0,Chii0
    READ(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
    READ(21) FSCX,FSLC,FSNC,FSLP,FSION,FSD0
    READ(21) rLn,rLT
    READ(21) Eb,RNB,PNBH,rNRF,RRF,PRFH,PNBCD
    READ(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
    READ(21) DLT,DT,EPS
    READ(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
    READ(21) rG1
    READ(21) rIPs,rIPe

    READ(21) T_TX,TMAX,NT,NQMAX,IERR
    READ(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

    READ(21) NGR,NGYR
    READ(21) NGT,NGYT
    READ(21) NGVV,NGYV
    READ(21) (GT(IGR), IGR=0, NGR)
    READ(21) (((GY(I,IGR,IGYR), I=0, NRMAX), IGR=0, NGR), IGYR=1, NGYR)
    READ(21) (GTX(I), I=0, NGT)
    READ(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYT)
    READ(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYV)
    !  END IF
    CLOSE(21)
    WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

    IF (NQMAX == 16) THEN
       X(LQm2,0:NRMAX) = 0
       X(LQm5,0:NRMAX) = BB
       X(LQb3,0:NRMAX) = 0
       X(LQn2,0:NRMAX) = 0
    END IF

    CALL TXCALM
    CALL TXCALV(X)
    CALL TXCALC
    CALL TXGLOB
    CALL TXWDAT
    CALL TXWDAT2

    RETURN
  END SUBROUTINE TXGLOD
end module file_io
