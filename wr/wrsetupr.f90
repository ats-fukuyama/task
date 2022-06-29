!   wrsetupr.f90

MODULE wrsetupr

  PRIVATE
  PUBLIC wr_setup_rays

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE wr_setup_rays(ierr)

    USE wrcomm
    USE wrsub,ONLY: wrnwtn,wrcale_i
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL:: TIME1
    INTEGER:: NS,NSS,NRAY
    REAL(rkind):: ZA1,ZA2,SINP2,SINT2,ANGZ,ANGPH
    EXTERNAL GUTIME

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

    ENDDO

    RETURN
  END SUBROUTINE wr_setup_rays

!     ***** ABSORBED POWER PROFILE*****

  SUBROUTINE WRSETUP(IERR)

    USE wrcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NRAY
    REAL(rkind):: ANGZ,ANGPH,SINP2,SINT2

    IERR=0
    WRITE(6,*) '@@@ point 1: mdlwri=',mdlwri

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


END MODULE wrsetupr
