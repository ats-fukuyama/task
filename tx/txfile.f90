!     $Id$
!***************************************************************
!
!   Write Data
!
!***************************************************************

SUBROUTINE TXWDAT
  use tx_commons, only : PI, PNeV, RB, X, LQe4, LQi4, PNiV, LQn1, NRMAX, PTeV, PTiV, WPT
  use tx_interface, only : INTG_F, rLINEAVE

  implicit none
  INTEGER(4) :: NR
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
       &     'NB(0.60)', rLINEAVE(0.6D0)  / 1.D20,  &
       &     'PF    ', PNeV(0) * 1.D20 / rNbar
  WRITE(6,'(1X,A," =",1PD9.2,2(2X,A,"=",1PD9.2))') &
       &     'Te(0)',    PTeV(0),  &
       &     'Ti(0)   ', PTiV(0),  &
       &     'Wst     ', WPT
  RETURN
END SUBROUTINE TXWDAT

!***************************************************************
!
!   Write Data 2
!
!***************************************************************

SUBROUTINE TXWDAT2

  use tx_commons, only : GTY, NGT
  implicit none
  REAL(4) :: gPNeMIN, gPNeMAX, gNB0MIN, gNB0MAX, gUiphMIN, gUiphMAX
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

!***************************************************************
!
!   Write Statistic Data
!
!***************************************************************

subroutine TXSTAT
  use tx_commons, only : VOLAVN, ALI, VLOOP, TAUE1, TAUE2, TAUEP, TAUEH, BETAA, &
       &                 BETAPA, BETAN, Q, ANSAV, rIp, PI, RA, NRA, NRMAX, R, &
       &                 rMui, Chii, UiphV
  implicit none
  integer(4) :: NR, NRL
  real(8) :: RL, rmuil, chiil, uiphl
  real(8) :: aitken2p

  RL = 0.3D0 * RA
  DO NR = 0, NRMAX-1
     IF(R(NR) <= RL.AND.R(NR+1) >= RL) THEN
        NRL = NR
        EXIT
     END IF
  END DO

  rmuil = aitken2p(rl,rmui(nrl),rmui(nrl+1),rmui(nrl+2),r(nrl),r(nrl+1),r(nrl+2))
  chiil = aitken2p(rl,chii(nrl),chii(nrl+1),chii(nrl+2),r(nrl),r(nrl+1),r(nrl+2))
  uiphl = aitken2p(rl,uiphv(nrl),uiphv(nrl+1),uiphv(nrl+2),r(nrl),r(nrl+1),r(nrl+2))

  write(6,'(1X,2(A27,1PD10.3,3X))') "Vol. ave. of neutrality  = ", VOLAVN
  write(6,'(1X,2(A27,1PD10.3,3X))') "Inductance               = ", ALI, &
       &                            "Loop voltage             = ", VLOOP
  write(6,'(1X,2(A27,1PD10.3,3X))') "Confinement time 1       = ", TAUE1, &
       &                            "Confinement time 2       = ", TAUE2
  write(6,'(1X,2(A27,1PD10.3,3X))') "L-mode scaling time      = ", TAUEP, &
       &                            "IPB98(y,2) scaling time  = ", TAUEH
  write(6,'(1X,2(A27,1PD10.3,3X))') "Beta                     = ", BETAA, &
       &                            "Poloidal beta            = ", BETAPA
  write(6,'(1X,2(A27,1PD10.3,3X))') "Normalized beta          = ", BETAN
  write(6,'(1X,2(A27,1PD10.3,3X))') "Line averaged e density  = ", ANSAV(1), &
       &                            "Greenwald density        = ", rIp / (PI * RA**2)
  write(6,'(1X,2(A27,1PD10.3,3X))') "Safety factor on axis    = ", Q(0), &
       &                            "Safety factor at sep.    = ", Q(NRA)
  write(6,'(1X,2(A27,1PD10.3,3X))') "Ion Prandtl num. at 0.3  = ", rmuil/chiil, &
       &                            "Ion tor. velocity at 0.3 = ", uiphl

end subroutine TXSTAT

!***********************************************************
!
!  LINE AVERAGE OF rN
!
!***********************************************************

REAL(8) FUNCTION rLINEAVE(Rho)

  use tx_commons, only : RA, NRMAX, PNeV
  implicit none
  REAL(8), INTENT(IN) :: Rho
  INTEGER(4) :: I, IR, NY = 100
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

!***************************************************************
!
!   Save transport data
!
!***************************************************************

SUBROUTINE TXSAVE
  use tx_commons, only : SLID,RA,RB,RC,RR,BB,PA,PZ,Zeff,PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ, &
       &              De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0,FSDFIX,FSCDBM,FSBOHM,FSPSCL, &
       &              PROFD,FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC,rLn,rLT, &
       &              Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2, &
       &              rNRF,RRF,RRF0,PRFH,PNBCD,PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       &              DT,EPS,ADV,tiny_cap,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,rG1, &
       &              rIPs,rIPe,T_TX,TMAX,NT,NQMAX,IERR,X,NGT,NGYTM,NGYVM,GTX,GVX,NGVV, &
       &              GTY,GVY,NLCMAX,NQM,GQY,NCM,NTCOIL,DltRPn,m_pol,n_tor, &
       &              MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB,MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD, &
       &              PNBMPD,thrp,kappa

  use tx_interface, only : TOUPPER

  implicit none
  INTEGER(4) :: IST, NQ, NR, NC, I, IGYT, IGYV
  character(len=100) :: TXFNAM, RCSId
  character(len=1) :: STR
  LOGICAL :: LEX

  RCSId = ' '

  DO 
     WRITE(6,*) '# INPUT : SAVE FILE NAME'
     CALL GUFLSH
     READ(*,'(A100)',IOSTAT=IST) TXFNAM
     IF (IST > 0) THEN
        CYCLE
     ELSE IF (IST < 0) THEN
        RETURN
     END IF
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
              WRITE(6,*) '# OLD FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)), &
                   &     ' ) IS ASSIGNED FOR OUTPUT.'
              EXIT
           ELSEIF (IST > 0) THEN
              WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
           END IF
        END IF
     ELSE
        OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='NEW',FORM='UNFORMATTED')
        IF (IST == 0) THEN
           WRITE(6,*) '# NEW FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)), &
                &     ' ) IS CREATED FOR OUTPUT.'
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
  WRITE(21) De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0
  WRITE(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
  WRITE(21) FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC
  WRITE(21) rLn,rLT
  WRITE(21) Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  WRITE(21) rNRF,RRF,RRF0,PRFH,PNBCD,PNBMPD
  WRITE(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
  WRITE(21) NTCOIL,DltRPn,kappa,m_pol,n_tor
  WRITE(21) DT,EPS,ADV,tiny_cap
  WRITE(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  WRITE(21) rG1
  WRITE(21) rIPs,rIPe
  WRITE(21) MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB
  WRITE(21) MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD

  WRITE(21) T_TX,TMAX,NT,NQMAX,IERR
  WRITE(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

  WRITE(21) NGT,NGYTM,NGYVM
  WRITE(21) (GTX(I), I=0, NGT)
  WRITE(21) (GVX(I), I=0, NGVV)
  WRITE(21) ((GTY(I,IGYT), I=0, NGT),  IGYT =1, NGYTM)
  WRITE(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  WRITE(21) (NLCMAX(NQ), NQ=1,NQM)
  WRITE(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQM)
  WRITE(21) (thrp(I), I=1, 2*NRMAX)
  CLOSE(21)
  WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

  RETURN
END SUBROUTINE TXSAVE

!***************************************************************
!
!   Load transport data
!
!***************************************************************

SUBROUTINE TXLOAD(IST)
  use tx_commons, only : allocate_txcomm, deallocate_txcomm, &
       &              SLID,RA,RB,RC,RR,BB,PA,PZ,Zeff,PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ, &
       &              De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0,FSDFIX,FSCDBM,FSBOHM,FSPSCL, &
       &              PROFD,FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC,rLn,rLT, &
       &              Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2, &
       &              rNRF,RRF,RRF0,PRFH,PNBCD,PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       &              DT,EPS,ADV,tiny_cap,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,rG1, &
       &              rIPs,rIPe,T_TX,TMAX,NT,NQMAX,IERR,X,NGT,NGYTM,NGYVM,GTX,GVX,NGVV, &
       &              GTY,GVY,NLCMAX,NQM,GQY,NCM,NTCOIL,DltRPn,m_pol,n_tor, &
       &              MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB,MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD, &
       &              PNBMPD,NGR,rIP,thrp,kappa,PTeV_FIX,PNeV_FIX,LQe1,LQe5,TAUE2, &
       &              rMU0,rMUb1,rMUb2
  use tx_variables
  use tx_coefficients, only : TXCALA
  use tx_parameter_control, only : TXPARM_CHECK

  implicit none
  integer(4), intent(out) :: IST
  INTEGER(4) :: NQ, NR, NC, NGYT, NGYV, I, IGYT, IGYV, NGTL
  character(len=100) ::  TXFNAM, RCSId
  character(len=8) :: LOADSLID
  LOGICAL :: LEX

  ! tmp : NGYT

  DO 
     WRITE(6,*) '# INPUT : LOAD FILE NAME'
     CALL GUFLSH
     READ(*,'(A100)',IOSTAT=IST) TXFNAM
     IF (IST > 0) THEN
        CYCLE
     ELSE IF (IST < 0) THEN
        RETURN
     END IF
     INQUIRE(FILE=TXFNAM,EXIST=LEX)
     IF (LEX) THEN
        OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',FORM='UNFORMATTED')
        IF (IST == 0) THEN
           WRITE(6,*) '# OLD FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)),  &
                &     ' ) IS ASSIGNED FOR INPUT.'
           EXIT
        ELSEIF (IST > 0) THEN
           WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
        END IF
     ELSE
        WRITE(6,*) 'XX  FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)), ' ) DOES NOT EXIST !'
     END IF
  END DO

  call allocate_txcomm(ierr)
  if(ierr /= 0) then
     call deallocate_txcomm
     write(6,*) "XX Allocation error : TXGLOD"
  end if

  READ(21,IOSTAT=IST) LOADSLID
  IF (IST > 0) THEN
     WRITE(6,*) 'XX READ ERROR in TXLOAD !'
     CLOSE(21)
     RETURN
  END IF
  !  IF(LOADSLID(1:5) == 'tx449') THEN
  READ(21) RCSId

  READ(21) RA,RB,RC,RR,BB
  READ(21) PA,PZ,Zeff
  READ(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
  READ(21) De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0
  READ(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
  READ(21) FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC
  READ(21) rLn,rLT
  READ(21) Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  READ(21) rNRF,RRF,RRF0,PRFH,PNBCD,PNBMPD
  READ(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
  READ(21) NTCOIL,DltRPn,kappa,m_pol,n_tor
  READ(21) DT,EPS,ADV,tiny_cap
  READ(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  READ(21) rG1
  READ(21) rIPs,rIPe
  READ(21) MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB
  READ(21) MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD

  READ(21) T_TX,TMAX,NT,NQMAX,IERR
  READ(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

  READ(21) NGT,NGYT,NGYV
  READ(21) (GTX(I), I=0, NGT)
  READ(21) (GVX(I), I=0, NGVV)
  READ(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYT)
  READ(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYV)
  READ(21) (NLCMAX(NQ), NQ=1,NQM)
  READ(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQM)
  READ(21) (thrp(I), I=1, 2*NRMAX)
  !  END IF
  CLOSE(21)
  WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  NGTL=NGT
  NGT=-1
  NGR=-1
  NGVV=-1
  rIP=rIPs

  CALL TXPARM_CHECK
  CALL TXCALM

  IF(rMUb1 == rMU0 .and. (PNBHT1 /= 0.D0 .OR. PNBHT2 /= 0.D0 .OR. PNBHP /= 0.D0)) THEN
     rMUb1 = 1.D0
     rMUb2 = rMU0
  END IF

  CALL TXCALV(X)

  PNeV_FIX(0:NRMAX) = X(LQe1,0:NRMAX)
  PTeV_FIX(0:NRMAX) = X(LQe5,0:NRMAX) / X(LQe1,0:NRMAX)

  CALL TXCALC
  CALL TXCALA
  CALL TXGLOB
  CALL TXWDAT
  CALL TXWDAT2

  ! TAUE2 uses data one step before the data was stored.
  ! Then TAUE2 is reconstituted by using the graphic data of TAUE2.
  TAUE2 = DBLE(GTY(NGTL,34))

  RETURN
END SUBROUTINE TXLOAD

!***************************************************************
!
!   Save graphic data
!
!***************************************************************

SUBROUTINE TXGSAV

  use tx_commons, only : SLID,RA,RB,RC,RR,BB,PA,PZ,Zeff,PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ, &
       &              De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0,FSDFIX,FSCDBM,FSBOHM,FSPSCL, &
       &              PROFD,FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC,rLn,rLT, &
       &              Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2, &
       &              rNRF,RRF,RRF0,PRFH,PNBCD,PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       &              DT,EPS,ADV,tiny_cap,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,rG1, &
       &              rIPs,rIPe,T_TX,TMAX,NT,NQMAX,IERR,X,NGT,NGYTM,NGYVM,GTX,GVX,NGVV, &
       &              GTY,GVY,NLCMAX,NQM,GQY,NCM,NGR,NGYRM,GT,GY,NTCOIL,DltRPn,m_pol,n_tor, &
       &              MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB,MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD, &
       &              PNBMPD,thrp,kappa

  implicit none
  INTEGER(4) :: IST, NQ, NR, NC, IGR, I, IGYR, IGYT, IGYV
  character(len=100) :: TXFNAM, RCSId
  character(len=1) :: STR
  LOGICAL :: LEX

  RCSId = ' '

  DO
     WRITE(6,*) '# INPUT : SAVE FILE NAME'
     CALL GUFLSH
     READ(*,'(A100)',IOSTAT=IST) TXFNAM
     IF (IST > 0) THEN
        CYCLE
     ELSE IF(IST < 0) THEN
        RETURN
     END IF
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
              WRITE(6,*) '# OLD FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)), &
                   &     ' ) IS ASSIGNED FOR OUTPUT.'
              EXIT
           ELSEIF (IST > 0) THEN
              WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
           END IF
        END IF
     ELSE
        OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='NEW',FORM='UNFORMATTED')
        IF (IST == 0) THEN
           WRITE(6,*) '# NEW FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)), &
                &     ' ) IS CREATED FOR OUTPUT.'
           EXIT
        ELSEIF (IST > 0) THEN
           WRITE(6,*) 'XX  NEW FILE OPEN ERROR !, IOSTAT = ', IST
        END IF
     END IF
  END DO

!!$    WRITE(21) SLID
!!$    WRITE(21) RCSId

  WRITE(21) RA,RB,RC,RR,BB
  WRITE(21) PA,PZ,Zeff
  WRITE(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
  WRITE(21) De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0
  WRITE(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
  WRITE(21) FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC
  WRITE(21) rLn,rLT
  WRITE(21) Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  WRITE(21) rNRF,RRF,RRF0,PRFH,PNBCD,PNBMPD
  WRITE(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
  WRITE(21) NTCOIL,DltRPn,kappa,m_pol,n_tor
  WRITE(21) DT,EPS,ADV,tiny_cap
  WRITE(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  WRITE(21) rG1
  WRITE(21) rIPs,rIPe
  WRITE(21) MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB
  WRITE(21) MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD

  WRITE(21) T_TX,TMAX,NT,NQMAX,IERR
  WRITE(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

  WRITE(21) NGR,NGYRM
  WRITE(21) NGT,NGYTM
  WRITE(21) NGVV,NGYVM
  WRITE(21) (GT(IGR), IGR=0, NGR)
  WRITE(21)(((GY(I,IGR,IGYR), I=0,NRMAX), IGR=0,NGR), IGYR=1,NGYRM)
  WRITE(21) (GTX(I), I=0, NGT)
  WRITE(21) (GVX(I), I=0, NGVV)
  WRITE(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYTM)
  WRITE(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYVM)
  WRITE(21) (NLCMAX(NQ), NQ=1,NQM)
  WRITE(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQM)
  WRITE(21) (thrp(I), I=1, 2*NRMAX)
  CLOSE(21)
  WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED IN THE FILE.'

  RETURN
END SUBROUTINE TXGSAV

!***************************************************************
!
!   Load graphic data
!
!***************************************************************

SUBROUTINE TXGLOD(IST)

  use tx_commons, only : &
       & allocate_txcomm, deallocate_txcomm, &
       & SLID,RA,RB,RC,RR,BB,PA,PZ,Zeff,PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ, &
       & De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0,FSDFIX,FSCDBM,FSBOHM,FSPSCL, &
       & PROFD,FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC,rLn,rLT, &
       & Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2, &
       & rNRF,RRF,RRF0,PRFH,PNBCD,PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       & DT,EPS,ADV,tiny_cap,NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,rG1, &
       & rIPs,rIPe,T_TX,TMAX,NT,NQMAX,IERR,X,NGT,NGYTM,NGYVM,GTX,GVX,NGVV, &
       & GTY,GVY,NLCMAX,NQM,GQY,NCM,NGR,NGYRM,GT,GY,NTCOIL,DltRPn,m_pol,n_tor, &
       & MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB,MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD, &
       & PNBMPD,rIP,thrp,kappa
  use tx_variables

  implicit none
  integer(4), intent(out) :: IST
  INTEGER(4) :: NQ, NR, NC, NGYR, NGYT, NGYV, IGR, I, IGYR, IGYT, IGYV
  character(len=100) :: TXFNAM, RCSId
  character(len=8) :: LOADSLID
  LOGICAL :: LEX

  ! tmp : NGYT

  DO 
     WRITE(6,*) '# INPUT : LOAD FILE NAME'
     CALL GUFLSH
     READ(*,'(A100)',IOSTAT=IST) TXFNAM
     IF (IST > 0) THEN
        CYCLE
     ELSE IF(IST < 0) THEN
        RETURN
     END IF
     INQUIRE(FILE=TXFNAM,EXIST=LEX)
     IF (LEX) THEN
        OPEN(21,FILE=TXFNAM,IOSTAT=IST,STATUS='OLD',FORM='UNFORMATTED')
        IF (IST == 0) THEN
           WRITE(6,*) '# OLD FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)),&
                &     ' ) IS ASSIGNED FOR INPUT.'
           EXIT
        ELSEIF (IST > 0) THEN
           WRITE(6,*) 'XX  OLD FILE OPEN ERROR !, IOSTAT = ', IST
        END IF
     ELSE
        WRITE(6,*) 'XX  FILE ( ', TXFNAM(1:LEN_TRIM(TXFNAM)), ' ) DOES NOT EXIST !'
     END IF
  END DO

  call allocate_txcomm(ierr)
  if(ierr /= 0) then
     call deallocate_txcomm
     write(6,*) "XX Allocation error : TXGLOD"
  end if

!!$    READ(21,IOSTAT=IST) LOADSLID
!!$    IF (IST > 0) THEN
!!$       WRITE(6,*) 'XX READ ERROR in TXGLOD !'
!!$       CLOSE(21)
!!$       RETURN
!!$    END IF
!!$    !  IF(LOADSLID(1:5) == 'tx449') THEN
!!$    READ(21) RCSId

  READ(21) RA,RB,RC,RR,BB
  READ(21) PA,PZ,Zeff
  READ(21) PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ
  READ(21) De0,Di0,rMue0,rMui0,WPM0,Chie0,Chii0
  READ(21) FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD
  READ(21) FSCX,FSLC,FSRP,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,MDLC
  READ(21) rLn,rLT
  READ(21) Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBH,PNBHP,PNBHT1,PNBHT2
  READ(21) rNRF,RRF,RRF0,PRFH,PNBCD,PNBMPD
  READ(21) PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV
  READ(21) NTCOIL,DltRPn,kappa,m_pol,n_tor
  READ(21) DT,EPS,ADV,tiny_cap
  READ(21) NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP
  READ(21) rG1
  READ(21) rIPs,rIPe
  READ(21) MODEG,MODEAV,MODEGL,MDLPCK,MDLWTB
  READ(21) MDLETA,MDFIXT,IDIAG,IGBDF,MDLNBD

  READ(21) T_TX,TMAX,NT,NQMAX,IERR
  READ(21) ((X(NQ,NR), NQ=1, NQMAX), NR=0, NRMAX)

  READ(21) NGR,NGYR
  READ(21) NGT,NGYT
  READ(21) NGVV,NGYV
  READ(21) (GT(IGR), IGR=0, NGR)
  READ(21) (((GY(I,IGR,IGYR), I=0, NRMAX), IGR=0, NGR), IGYR=1, NGYR)
  READ(21) (GTX(I), I=0, NGT)
  READ(21) (GVX(I), I=0, NGVV)
  READ(21) ((GTY(I,IGYT), I=0, NGT), IGYT =1, NGYT)
  READ(21) ((GVY(I,IGYV), I=0, NGVV), IGYV =1, NGYV)
  READ(21) (NLCMAX(NQ), NQ=1,NQM)
  READ(21) (((GQY(NR,NC,NQ), NR=0, NRMAX), NC=1, NCM), NQ=1, NQM)
  READ(21) (thrp(I), I=1, 2*NRMAX)
  !  END IF
  CLOSE(21)
  WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  rIp = rIps

  CALL TXCALM
!!$    CALL TXCALV(X)
!!$    CALL TXCALC
!!$    CALL TXGLOB
!!$    CALL TXWDAT
!!$    CALL TXWDAT2

  RETURN
END SUBROUTINE TXGLOD
