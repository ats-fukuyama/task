!   obcalc.f90

MODULE obcalc

  PRIVATE
  PUBLIC ob_calc

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE OB_CALC(IERR)

    USE obcomm
    USE obexec,ONLY: ob_exec
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(4):: TIME1,TIME2
    INTEGER:: NS,NSS,NOBT
    REAL(rkind):: ZA1,ZA2,SINP2,SINT2,ANGZ,ANGPH
    REAL(rkind):: RKR0_1,RKR0_2,EPARA,EPERP
    REAL(rkind):: RKR0_11,RKR0_12,RKR0_21,RKR0_22

    IERR=0

    CALL GUTIME(TIME1)

    IF(MDLOBI.LT.100) CALL ob_setup(IERR)

    DO NOBT=1,NOBTMAX
       RF=RFIN(NOBT)
       RPI=RPIN(NOBT)
       ZPI=ZPIN(NOBT)
       PHII=PHIIN(NOBT)
       RKR0=RKRIN(NOBT)
       RNZI=RNZIN(NOBT)
       RNPHII=RNPHIIN(NOBT)
       ANGZ=ANGZIN(NOBT)
       ANGPH=ANGPHIN(NOBT)
       MODEW=MODEWIN(NOBT)
       UUI=UUIN(NOBT)

       IF(MDLOBI.EQ.100)THEN
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

       OBTIN(1,NOBT)=RF
       OBTIN(2,NOBT)=RPI
       OBTIN(3,NOBT)=ZPI
       OBTIN(4,NOBT)=PHII
       OBTIN(5,NOBT)=RKR0
       OBTIN(6,NOBT)=ANGZ
       OBTIN(7,NOBT)=ANGPH
       OBTIN(8,NOBT)=UUI

       CALL ob_exec(NOBT,IERR)
       IF(IERR.NE.0) cycle

    ENDDO

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

    RETURN
  END SUBROUTINE OB_CALC

  ! setup intial condition in various forms
  
  SUBROUTINE ob_setup(IERR)

    USE obcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NOBT
    REAL(rkind):: ANGZ,ANGPH,SINP2,SINT2

    IERR=0

    IF(MDLOBI.EQ.0) THEN
       WRITE(6,*) &
            '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
       WRITE(6,'(1PE12.4,0P7F9.2)') &
            RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
    ELSEIF(MDLOBI.EQ.1.OR.MDLOBI.EQ.2) THEN
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
       IF(MDLOBI.EQ.1) THEN
          WRITE(6,*) &
               '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
       ELSEIF(MDLOBI.EQ.2) THEN
          WRITE(6,*) &
               '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
       ENDIF
    ELSEIF(MDLOBI.EQ.11) THEN
       WRITE(6,*) &
            '# default values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
       WRITE(6,'(1PE12.4,0P7F9.2)') &
            RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
       PHII=0.D0
    ENDIF

!     --- Input for each obt tracing ---

    DO NOBT=1,NOBTMAX

1         WRITE(6,*) '# NOBT = ',NOBT
          IF(MDLOBI.EQ.0) THEN
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
          ELSEIF(MDLOBI.EQ.1) THEN
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
          ELSEIF(MDLOBI.EQ.2) THEN
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
             WRITE(6,*) 'XX MDLOBI=2 IS NOT SUPPORTED YET.'
             GOTO 1
          ELSEIF(MDLOBI.EQ.11) THEN
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
          RFIN(NOBT)=RF
          RPIN(NOBT)=RPI
          ZPIN(NOBT)=ZPI
          PHIIN(NOBT)=PHII
          RKRIN(NOBT)=RKR0
          RNZIN(NOBT)=RNZI
          RNPHIIN(NOBT)=RNPHII
          ANGZIN(NOBT)=ANGZ
          ANGPHIN(NOBT)=ANGPH
          MODEWIN(NOBT)=MODEW
          UUIN(NOBT)=UUI
    END DO
    RETURN

9000 CONTINUE
    IERR=10
    RETURN
  END SUBROUTINE ob_setup


END MODULE obcalc
