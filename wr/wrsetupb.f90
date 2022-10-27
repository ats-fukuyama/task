!   wrsetupb.f

MODULE wrsetupb

  PRIVATE
  PUBLIC wr_setup_beams

CONTAINS

!     ***** Beam tracing module *****

  SUBROUTINE wr_setup_beams(ierr)

    USE wrcomm
    USE wrsub,ONLY: wrcale,wrnwtn
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: NRAY
    REAL(rkind):: ANGZ,ANGPH,SINP2,SINT2
    EXTERNAL GUTIME

    ierr=0

    IF(MDLWRI.EQ.0) THEN
       WRITE(6,*) '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
       WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
       WRITE(6,'(1PE12.4,0P6F10.3)') &
            RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
       WRITE(6,'(12X,1P5E12.4)') &
            RCURVA,RCURVB,RBRADA,RBRADB,UUI
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
          WRITE(6,*) '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
          WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
          WRITE(6,'(1PE12.4,0P6F10.3)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
          WRITE(6,'(12X,1P5E12.4)') &
               RCURVA,RCURVB,RBRADA,RBRADB,UUI
       ELSEIF(MDLWRI.EQ.2) THEN
          WRITE(6,*) '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH'
          WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
          WRITE(6,'(1PE12.4,0P6F10.3)') &
               RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
          WRITE(6,'(12X,1P5E12.4)') &
               RCURVA,RCURVB,RBRADA,RBRADB,UUI
       ENDIF
    ENDIF

    DO NRAY=1,NRAYMAX

       IF(MDLWRI.LT.10) then
1         WRITE(6,*) '# NRAY=' ,NRAY
          IF(MDLWRI.EQ.0) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII, &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
          ELSEIF(MDLWRI.EQ.1) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH, &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
          ELSEIF(MDLWRI.EQ.2) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             WRITE(6,*) 'XX MDLWRI=2 IS NOT SUPPORTED YET.'
             GOTO 1
          ENDIF
       ELSE
          RF=RFIN(NRAY)
          RPI=RPIN(NRAY)
          ZPI=ZPIN(NRAY)
          PHII=PHIIN(NRAY)
          RKR0=RKRIN(NRAY)
          RNZI=RNZIN(NRAY)
          RNPHII=RNPHIIN(NRAY)
          ANGZ=ANGZIN(NRAY)
          ANGPH=ANGPHIN(NRAY)
          UUI=UUIN(NRAY)
          RCURVA=RCURVAIN(NRAY)
          RCURVB=RCURVBIN(NRAY)
          RBRADA=RBRADAIN(NRAY)
          RBRADB=RBRADBIN(NRAY)

          IF(MDLWRI.EQ..10)THEN
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
          ELSE
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
          ENDIF
       ENDIF

       RAYIN(1,NRAY)=RF
       RAYIN(2,NRAY)=RPI
       RAYIN(3,NRAY)=ZPI
       RAYIN(4,NRAY)=PHII
       RAYIN(5,NRAY)=RKR0
!     
       IF(MDLWRI.EQ.0)THEN
          RAYIN(6,NRAY)=RNZI
          RAYIN(7,NRAY)=RNPHII
       ELSE
          RAYIN(6,NRAY)=ANGZ
          RAYIN(7,NRAY)=ANGPH
       ENDIF

       RAYIN(8,NRAY)=UUI

    END DO
    RETURN

9000 CONTINUE
    ierr=9000
    RETURN
  END SUBROUTINE wr_setup_beams
    
END MODULE wrsetupb
