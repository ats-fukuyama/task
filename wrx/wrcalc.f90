!   wrcalc.f90

MODULE wrcalc

  PRIVATE
  PUBLIC wr_calc

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE WR_CALC(IERR)

    USE wrcomm
    USE wrexecr,ONLY: wr_execr
    USE wrsub,ONLY: wr_cold_rkr0,wrnwtn,wrcale_i
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(4):: TIME1,TIME2
    INTEGER:: NS,NSS,NRAY
    REAL(rkind):: ZA1,ZA2,SINP2,SINT2,ANGZ,ANGPH
    REAL(rkind):: RKR0_1,RKR0_2,EPARA1,EPERP1,EPARA2,EPERP2

    IERR=0

    CALL GUTIME(TIME1)

    CALL WRSETUP(IERR)

!     --- eliminate disp factor for same z/a species ---

    DO NS=1,NSMAX
       NSDP(NS)=1
       DO NSS=1,NS-1
          ZA1=PZ(NS)/PA(NS)
          ZA2=PZ(NSS)/PA(NSS)
          IF(ABS(ZA1-ZA2).LE.1.D-8) NSDP(NS)=0
       ENDDO
    ENDDO

    DO NRAY=1,NRAYMAX
       RF=RFIN(NRAY)
       RPI=RPIN(NRAY)
       ZPI=ZPIN(NRAY)
       PHII=PHIIN(NRAY)
       RKR0=RKRIN(NRAY)
       RNZI=RNZIN(NRAY)
       RNPHII=RNPHIIN(NRAY)
       ANGZ=ANGZIN(NRAY)
       ANGPH=ANGPHIN(NRAY)
       MODEW=MODEWIN(NRAY)
       UUI=UUIN(NRAY)

       IF(MDLWRI.EQ.100)THEN
          WRITE(6,*) &
               '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
       ELSE
          WRITE(6,*) &
               '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
          SINP2=SIN(ANGZ*PI/180.D0)**2
          SINT2=SIN(ANGPH*PI/180.D0)**2
          RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
          RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
          IF(ANGZ.LT.0.D0) RNZI=-RNZI
          IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
       ENDIF

       RAYIN(1,NRAY)=RF
       RAYIN(2,NRAY)=RPI
       RAYIN(3,NRAY)=ZPI
       RAYIN(4,NRAY)=PHII
       RAYIN(5,NRAY)=RKR0
       RAYIN(6,NRAY)=ANGZ
       RAYIN(7,NRAY)=ANGPH
       RAYIN(8,NRAY)=UUI

       RKZI  =2.D6*PI*RF*RNZI  /VC
       RKPHII=2.D6*PI*RF*RNPHII/VC

       SELECT CASE(MODEW)
       CASE(1,2)
          CALL WR_COLD_RKR0(RKR0_1,RKR0_2)
          RKR0=RKR0_1
          CALL WRCALE_I(EPARA1,EPERP1)
          RKR0=RKR0_2
          CALL WRCALE_I(EPARA2,EPERP2)
          WRITE(6,'(A,1P3E12.4)') 'RKR0_1,EPARA,EPWEP=',RKR0_1,EPARA1,EPERP1
          WRITE(6,'(A,1P3E12.4)') 'RKR0_2,EPARA,EPWEP=',RKR0_2,EPARA2,EPERP2
          IF(MODEW.EQ.1) THEN
             RKR0=RKR0_1
          ELSE
             RKR0=RKR0_2
          END IF
       END SELECT

       RKZI  =2.D6*PI*RF*RNZI  /VC
       RKPHII=2.D6*PI*RF*RNPHII/VC
       CALL WRNWTN(IERR)  ! input RKR0,RKZI,RKPHII; output RKRI
       IF(IERR.NE.0) cycle

       CALL wr_execr(NRAY,IERR)
       IF(IERR.NE.0) cycle

       WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
            '    RKRI=  ',RKRI,  '  PABS/PIN=', &
            1.D0-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
    ENDDO

    CALL WRAPWR

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

    RETURN

9000 RETURN
  END SUBROUTINE WR_CALC

!     ***** ABSORBED POWER PROFILE*****

  SUBROUTINE WRSETUP(IERR)

    USE wrcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS,NSS,NRAY
    REAL(rkind):: ZA1,ZA2,ANGZ,ANGPH,SINP2,SINT2

    IERR=0

    IF(MDLWRI.EQ.0) THEN
       WRITE(6,*) &
            '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
       WRITE(6,'(1PE12.4,0P7F9.2)') &
            RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
    ELSEIF(MDLWRI.EQ.1.OR.MDLWRI.EQ.2) THEN
       IF(ABS(RNZI).GT.1.D0) THEN
          ANGZ=0.D0
          ANGPH=0.D0
       ELSE
          ANGZ=ASIN(RNZI)*180.D0/PI
          IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
             ANGZ=0.D0
          ELSE
             ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
          ENDIF
          IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
             ANGPH=0.D0
          ELSE
             ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
          ENDIF
       ENDIF
       IF(MDLWRI.EQ.1) THEN
          WRITE(6,*) &
               '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
       ELSEIF(MDLWRI.EQ.2) THEN
          WRITE(6,*) &
               '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
       ENDIF
    ELSEIF(MDLWRI.EQ.11) THEN
       WRITE(6,*) &
            '# default values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
       WRITE(6,'(1PE12.4,0P7F9.2)') &
            RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
       PHII=0.D0
    ENDIF

!     --- Each ray tracing ---

    DO NRAY=1,NRAYMAX

       IF(MDLWRI.LT.100) THEN

1         WRITE(6,*) '# NRAY = ',NRAY
          IF(MDLWRI.EQ.0) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
             IF(ABS(RNZI).GT.1.D0) THEN
                ANGZ=0.D0
                ANGPH=0.D0
             ELSE
                ANGZ=ASIN(RNZI)*180.D0/PI
                IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
                   ANGZ=0.D0
                ELSE
                   ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
                ENDIF
                IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
                   ANGPH=0.D0
                ELSE
                   ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
                ENDIF
             ENDIF
          ELSEIF(MDLWRI.EQ.1) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1.D0-SINT2)/(1.D0-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1.D0-SINP2)/(1.D0-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
          ELSEIF(MDLWRI.EQ.2) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
             WRITE(6,*) 'XX MDLWRI=2 IS NOT SUPPORTED YET.'
             GOTO 1
          ELSEIF(MDLWRI.EQ.11) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
             IF(ABS(RNZI).GT.1.D0) THEN
                ANGZ=0.D0
                ANGPH=0.D0
             ELSE
                ANGZ=ASIN(RNZI)*180.D0/PI
                IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
                   ANGZ=0.D0
                ELSE
                   ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
                ENDIF
                IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
                   ANGPH=0.D0
                ELSE
                   ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
                ENDIF
             ENDIF
          ENDIF
          RFIN(NRAY)=RF
          RPIN(NRAY)=RPI
          ZPIN(NRAY)=ZPI
          PHIIN(NRAY)=PHII
          RKRIN(NRAY)=RKR0
          RNZIN(NRAY)=RNZI
          RNPHIIN(NRAY)=RNPHII
          ANGZIN(NRAY)=ANGZ
          ANGPHIN(NRAY)=ANGPH
          MODEWIN(NRAY)=MODEW
          UUIN(NRAY)=UUI
       ENDIF
    END DO
    RETURN

9000 CONTINUE
    IERR=10
    RETURN
  END SUBROUTINE WRSETUP

  SUBROUTINE WRAPWR

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind):: PWRRAY(NRDIVMAX,NRAYMAX)
    REAL(rkind):: PWR(NRDIVMAX),AJR(NRDIVMAX)
    INTEGER:: NRAY,NSTP,NR,NRS1,NRS2,NDR,LOCMAX
    REAL(rkind):: DRHO,XL,YL,ZL,RHON1,RHON2,SDR,DELPW
    REAL(rkind):: PWRMAX,RHOMAX,DPWR,DDPWR

!     ----- CALCULATE RADIAL DEPOSITION PROFILE -----

      DRHO=1.D0/(NRDIVMAX-1)
      DO NRAY=1,NRAYMAX
         DO NR=1,NRDIVMAX
            PWRRAY(NR,NRAY)=0.D0
         ENDDO
      ENDDO

      DO NRAY=1,NRAYMAX
         DO NSTP=0,NSTPMAX_NRAY(NRAY)-1
            XL=RAYS(1,NSTP,NRAY)
            YL=RAYS(2,NSTP,NRAY)
            ZL=RAYS(3,NSTP,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON1)
            IF(RHON1.LE.1.D0) THEN
               NRS1=INT(RHON1/DRHO)+1
               XL=RAYS(1,NSTP+1,NRAY)
               YL=RAYS(2,NSTP+1,NRAY)
               ZL=RAYS(3,NSTP+1,NRAY)
               CALL PL_MAG_OLD(XL,YL,ZL,RHON2)
               NRS2=INT(RHON2/DRHO)+1
               NDR=ABS(NRS2-NRS1)
               IF(NDR.EQ.0) THEN
                  PWRRAY(NRS1,NRAY) &
                       =PWRRAY(NRS1,NRAY)+RAYS(8,NSTP+1,NRAY)
               ELSE IF(NRS1.LT.NRS2) THEN
                  SDR=(RHON2-RHON1)/DRHO
                  DELPW=RAYS(8,NSTP+1,NRAY)/SDR
                  PWRRAY(NRS1,NRAY)=PWRRAY(NRS1,NRAY) &
                       +(DBLE(NRS1)-RHON1/DRHO)*DELPW
                  DO NR=NRS1+1,NRS2-1
                     PWRRAY(NR,NRAY)=PWRRAY(NR,NRAY)+DELPW
                  END DO
                  PWRRAY(NRS2,NRAY)=PWRRAY(NRS2,NRAY) &
                       +(RHON2/DRHO-DBLE(NRS2-1))*DELPW
               ELSE
                  SDR=(RHON1-RHON2)/DRHO
                  DELPW=RAYS(8,NSTP+1,NRAY)/SDR
                  PWRRAY(NRS2,NRAY)=PWRRAY(NRS2,NRAY) &
                       +(DBLE(NRS2)-RHON2/DRHO)*DELPW
                  DO NR=NRS2+1,NRS1-1
                     PWRRAY(NR,NRAY)=PWRRAY(NR,NRAY)+DELPW
                  END DO
                  PWRRAY(NRS1,NRAY)=PWRRAY(NRS1,NRAY) &
                       +(RHON1/DRHO-DBLE(NRS1-1))*DELPW
               END IF
            END IF
         END DO
      END DO

!     ----- divide by area -----

      DO NRAY=1,NRAYMAX
         DO NR=1,NRDIVMAX
            PWRRAY(NR,NRAY)=PWRRAY(NR,NRAY) &
                       /(2*PI*(DBLE(NR)-0.5D0)*DRHO*DRHO)
         ENDDO
      ENDDO

!     ----- weight of power for each ray-------

      DO NR=1,NRDIVMAX
         PWR(NR)=0.D0
         DO NRAY=1,NRAYMAX
            PWR(NR)=PWR(NR)+PWRRAY(NR,NRAY)
         ENDDO
      ENDDO

      DO NR=1,NRDIVMAX
         AJR(NR)=0.D0
      ENDDO

!     ----- find location of absorbed power peak -----

      PWRMAX=0.D0
      LOCMAX=0
      DO NR=1,NRDIVMAX
         IF(PWR(NR).GT.PWRMAX) THEN
            PWRMAX=PWR(NR)
            LOCMAX=NR
         ENDIF
      END DO
      IF(LOCMAX.LE.1) THEN
         RHOMAX=0.D0
      ELSE IF(LOCMAX.EQ.NRDIVMAX) THEN
         RHOMAX=1.D0
      ELSE
         DPWR =(PWR(LOCMAX+1)-PWR(LOCMAX-1))/(2.D0*DRHO)
         DDPWR=(PWR(LOCMAX+1)-2*PWR(LOCMAX)+PWR(LOCMAX-1))/DRHO**2
         RHOMAX=(LOCMAX-0.5D0)/(NRDIVMAX-1.D0)-DPWR/DDPWR
         PWRMAX=PWRMAX-DPWR**2/(2.D0*DDPWR)
      ENDIF
      WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
           '    PWRMAX=',PWRMAX,'  AT RHON =',RHOMAX

      RETURN
  END SUBROUTINE WRAPWR

END MODULE wrcalc
