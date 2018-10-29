!   wrcalc.f90

MODULE wrcalc

  PRIVATE
  PUBLIC wr_calc,wrcale,wrnwtn,wrcalk

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE WR_CALC

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind):: Y(NEQ)
    REAL(4):: TIME1,TIME2
    INTEGER:: NS,NSS,NSTP,NRAY,IERR
    REAL(rkind):: ZA1,ZA2,SINP2,SINT2,ANGZ,ANGPH,RKRI,RKZI,RKPHII

    CALL GUTIME(TIME1)

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

!     --- eliminate disp factor for same z/a species ---

    DO NS=1,NSMAX
       NSDP(NS)=1
       DO NSS=1,NS-1
          ZA1=PZ(NS)/PA(NS)
          ZA2=PZ(NSS)/PA(NSS)
          IF(ABS(ZA1-ZA2).LE.1.D-8) NSDP(NS)=0
       ENDDO
    ENDDO

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
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
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
       ENDIF

       RAYIN(1,NRAY)=RF
       RAYIN(2,NRAY)=RPI
       RAYIN(3,NRAY)=ZPI
       RAYIN(4,NRAY)=PHII
       RAYIN(5,NRAY)=RKR0
       RAYIN(6,NRAY)=ANGZ
       RAYIN(7,NRAY)=ANGPH
       RAYIN(8,NRAY)=UUI

       RKRI  = RKR0
       RKZI  =2.D6*PI*RF*RNZI  /VC
       RKPHII=2.D6*PI*RF*RNPHII/VC
       CALL WRNWTN(RKRI,RKZI,RKPHII,IERR)
       WRITE(6,'(A,1PE16.8)') 'RKRI= ',RKRI
       IF(IERR.NE.0) GOTO 1200

       IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
          Y(1)= RPI
          Y(2)= PHII
          Y(3)= ZPI
          Y(4)= RKRI
          Y(5)= RKPHII
          Y(6)= RKZI
       ELSE
          Y(1)= RPI*COS(PHII)
          Y(2)= RPI*SIN(PHII)
          Y(3)= ZPI
          Y(4)= RKRI*COS(PHII)-RKPHII*SIN(PHII)
          Y(5)= RKRI*SIN(PHII)+RKPHII*COS(PHII)
          Y(6)= RKZI
       ENDIF
       Y(7)= UUI
       IF(MDLWRQ.EQ.0) THEN
          CALL WRRKFT(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
       ELSEIF(MDLWRQ.EQ.1) THEN
          CALL WRRKFT_WITHD0(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
       ELSEIF(MDLWRQ.EQ.2) THEN
          CALL WRRKFT_WITHMC(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
       ELSEIF(MDLWRQ.EQ.3) THEN
          CALL WRRKFT_ODE(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
       ELSEIF(MDLWRQ.EQ.4) THEN
          CALL WRSYMP(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
       ELSE
          WRITE(6,*) 'XX WRCALC: unknown MDLWRQ =', MDLWRQ
          RETURN
       ENDIF
       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RAYRB1(NSTP,NRAY)=0.D0
          RAYRB2(NSTP,NRAY)=0.D0
       END DO
       CALL WRCALE(RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY),NRAY)
       WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
            '    RKRI=  ',RKRI,  '  PABS/PIN=', &
            1.D0-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
1200   CONTINUE
    ENDDO

    CALL WRAPWR

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

9000 RETURN
  END SUBROUTINE WR_CALC

!************************************************************************

  SUBROUTINE WRRKFT(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3),F(NEQ)
    INTEGER:: INIT,NSTPLIM,NSTP,I,IOX,ID
    REAL(rkind):: X0,XE,OMG,PSIN,RL,PHIL,ZL,RKRL,Y7,RHON
    REAL(rkind):: RNPHI_IDEI,DELTA,RKPARA,RKPERP

    X0 = 0.D0
    XE = DELS
    INIT = 1
    NSTPLIM=INT(SMAX/DELS)

    NSTP=0
    YN(0,NSTP)=X0
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0
    OMG=2.D6*PI*RF
    IOX=0

    DO NSTP = 1,NSTPLIM
       Y7=Y(7)
       CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=Y7-YM(7)

       CALL PL_MAG_OLD(YM(1),YM(2),YM(3),PSIN)
       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          RL  =YM(1)
          PHIL=ASIN(YM(2)/(2.D0*PI*RR))
          ZL  =YM(3)
          RKRL=YM(4)
       ELSE IF(MODELG.EQ.11) THEN
          RL  =YM(1)
          PHIL=YM(2)
          ZL  =YM(3)
          RKRL=YM(4)
       ELSE
          RL  =SQRT(YM(1)**2+YM(2)**2)
          PHIL=ATAN2(YM(2),YM(1))
          ZL  =YM(3)
          RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
       ENDIF
       RNPHI_IDEI= -YM(4)*VC/OMG*SIN(PHIL) + YM(5)*VC/OMG*COS(PHIL)
       DELTA=DISPXR( YM(1), YM(2), YM(3), YM(4), YM(5), YM(6), OMG)
       RKPARA=YM(4)*BNX+YM(5)*BNY+YM(6)*BNZ
       RKPERP=SQRT((YM(4)*YM(4)+YM(5)*YM(5)+YM(6)*YM(6))-RKPARA**2)

       IF(MDLWRW.GE.1) THEN
          ID=0
          SELECT CASE(MDLWRW)
          CASE(1)
             ID=1
          CASE(2)
             IF(MOD(NSTP-1,10).EQ.0) ID=1
          CASE(3)
             IF(MOD(NSTP-1,100).EQ.0) ID=1
          CASE(4)
             IF(MOD(NSTP-1,1000).EQ.0) ID=1
          CASE(5)
             IF(MOD(NSTP-1,10000).EQ.0) ID=1
          END SELECT
          IF(ID.EQ.1) &
               WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
       ENDIF

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS
       IF(Y(7).LT.UUMIN) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA*RKAP) THEN
             NNSTP = NSTP
             GOTO 11
          ENDIF
       ENDIF
    END DO
    NNSTP=NSTPLIM
!     
11  IF(YN(7,NNSTP).LT.0.D0) THEN
       YN(7,NNSTP)=0.D0
    ENDIF
    IF(MDLWRW.GE.1) THEN
       IF(ID.EQ.0) &
            WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
    ENDIF

    RETURN
  END SUBROUTINE WRRKFT

!************************************************************************

  SUBROUTINE WRRKFT_WITHD0(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3)
    INTEGER:: INIT,NSTPLIM,NSTP,I,IX,ID
    REAL(rkind):: X0,XE,OMG,Y7,DELTA,RHON,RL,PHIL,ZL,RKRL
    REAL(rkind):: RNPHI_IDEI,RKPARA,RKPERP

      X0 = 0.D0
      XE = DELS
      INIT = 1
      NSTPLIM=INT(SMAX/DELS)
      NSTP=0
      YN(0,NSTP)=X0
      DO I=1,7
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(8,NSTP)=0.D0
      OMG=2.D6*PI*RF

      DO NSTP = 1,NSTPLIM
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),OMG)
         IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,YK)
            DO I=1,3
               YM(I+3) = YK(I)
            END DO
         END IF 
         YN(0,NSTP)=XE
         DO I=1,7
            YN(I,NSTP)=YM(I)
         ENDDO
         YN(8,NSTP)=Y7-YM(7)

         CALL PL_MAG_OLD(YM(1),YM(2),YM(3),RHON)
         IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
            RL  =YM(1)
            IF(RR.EQ.0.D0) THEN
               PHIL=YM(2)
            ELSE
               PHIL=ASIN(YM(2)/(2.D0*PI*RR))
            ENDIF
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE IF(MODELG.EQ.11) THEN
            RL  =YM(1)
            PHIL=YM(2)
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE
            RL  =SQRT(YM(1)**2+YM(2)**2)
            PHIL=ATAN2(YM(2),YM(1))
            ZL  =YM(3)
            RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         ENDIF
         RNPHI_IDEI= -YM(4)*VC/OMG*SIN(PHIL) + YM(5)*VC/OMG*COS(PHIL)
         DELTA=DISPXR( YM(1), YM(2), YM(3), YM(4), YM(5), YM(6), OMG )
         RKPARA=YM(4)*BNX+YM(5)*BNY+YM(6)*BNZ
         RKPERP=SQRT((YM(4)*YM(4)+YM(5)*YM(5)+YM(6)*YM(6))-RKPARA**2)

!IDEI	 WRITE(6,*) 'BX=',BNX*BABS,'BY=',BNY*BABS,'BZ=',BNZ*BABS
!	 WRITE(6,6001) YN(0,IT),RL,ZL,PHIL,YN(1,IT), &
!                      YN(2,IT),YN(3,IT),YN(6,IT)*VC/OMG,YN(7,IT), &
!      	               DELTA,RKPARA*VC/OMG,RKPERP*VC/OMG,RKRL, &
!      	               RNPHI_IDEI
! 6001    FORMAT(1H ,1P14E13.5)

         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(NSTP-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(NSTP-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(NSTP-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(NSTP-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) &
                 WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
         ENDIF

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
         IF(MODELG.LE.10) THEN
            CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
            IF(RHON.GT.RB/RA*RKAP) THEN
               NNSTP = NSTP
               GOTO 11
            ENDIF        
         END IF
      END DO
      NNSTP=NSTPLIM
!     
 11   IF(YN(7,NNSTP).LT.0.D0) THEN
         YN(7,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) &
              WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NNSTP)
      ENDIF

      RETURN
  END SUBROUTINE WRRKFT_WITHD0

!************************************************************************

  SUBROUTINE WRRKFT_ODE(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3)
    REAL(rkind):: X0,XE,RL,PHIL,ZL,KL,Y7,RHON,RKRL
    INTEGER:: NSTPLIM,NSTP,I,ID

      X0 = 0.D0
      XE = DELS     
      NSTPLIM=INT(SMAX/DELS)
      NSTP=0
      YN(0,NSTP)=X0
      DO I=1,7
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(8,NSTP)=0.D0

      DO NSTP = 1,NSTPLIM
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         YN(0,NSTP)=XE
         DO I=1,7
            YN(I,NSTP)=YM(I)
         ENDDO
         IF(YM(7).GT.0.D0) THEN
            YN(8,NSTP)=Y7-YM(7)
         ELSE
            YN(8,NSTP)=Y7
         ENDIF

         IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
            RL  =YM(1)
            PHIL=ASIN(YM(2)/(2.D0*PI*RR))
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE
            RL  =SQRT(YM(1)**2+YM(2)**2)
            PHIL=ATAN2(YM(2),YM(1))
            ZL  =YM(3)
            RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         ENDIF

         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(NSTP-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(NSTP-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(NSTP-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(NSTP-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) &
                 WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
         ENDIF

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA*RKAP) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
      END DO
      NNSTP=NSTPLIM
!     
 11   IF(YN(7,NNSTP).LT.0.D0) THEN
         YN(7,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) &
              WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NNSTP)
      ENDIF

      RETURN
  END SUBROUTINE WRRKFT_ODE

!************************************************************************
  
  SUBROUTINE WRRKFT_WITHMC(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3),F(NEQ)
    REAL(rkind):: X0,XE,OMG,RL,RKRL,OXEFF,RHON,PHIL,ZL,Y7
    REAL(rkind):: RNPHI_IDEI,DELTA,RKPARA,RKPERP
    INTEGER:: INIT,NSTPLIM,NSTP,I,IOX,ID

      X0 = 0.D0
      XE = DELS
      INIT = 1
      NSTPLIM=INT(SMAX/DELS)
      NSTP=0
      YN(0,NSTP)=X0
      DO I=1,7
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(8,NSTP)=0.D0
      OMG=2.D6*PI*RF
      IOX=0

      DO NSTP = 1,NSTPLIM
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),OMG)

         IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,YK)
            DO I=1,3
               YM(I+3) = YK(I)
            END DO
         END IF 

!   --- Mode conversion

         RL  =SQRT(YM(1)**2+YM(2)**2)
         RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         IF ( RKRL.GE.0.D0.AND.IOX.EQ.0 ) THEN
            CALL WRMODCONV_OX(IOX,YM,F,OXEFF)
            IF(IOX.GE.100000) THEN
               WRITE(6,*) 'ERROR in WRMODCON_OX routine IOX=',IOX
            ELSE 
               DO I=1,NEQ
                  YM(I) = F(I)
               END DO
            END IF
         ENDIF

         YN(0,NSTP)=XE
         DO I=1,7
            YN(I,NSTP)=YM(I)
         ENDDO
         IF(YM(7).GT.0.D0) THEN
            YN(8,NSTP)=Y7-YM(7)
         ELSE
            YN(8,NSTP)=Y7
         ENDIF

         CALL PL_MAG_OLD(YM(1),YM(2),YM(3),RHON)
         IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
            RL  =YM(1)
            PHIL=ASIN(YM(2)/(2.D0*PI*RR))
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE
            RL  =SQRT(YM(1)**2+YM(2)**2)
            PHIL=ATAN2(YM(2),YM(1))
            ZL  =YM(3)
            RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         ENDIF
         RNPHI_IDEI= -YM(4)*VC/OMG*SIN(PHIL) + YM(5)*VC/OMG*COS(PHIL)
         DELTA=DISPXR( YM(1), YM(2), YM(3), YM(4), YM(5), YM(6), OMG )
         RKPARA=YM(4)*BNX+YM(5)*BNY+YM(6)*BNZ
         RKPERP=SQRT((YM(4)*YM(4)+YM(5)*YM(5)+YM(6)*YM(6))-RKPARA**2)

         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(NSTP-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(NSTP-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(NSTP-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(NSTP-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) &
                 WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
         ENDIF

!IDEI	 WRITE(6,*) 'BX=',BNX*BABS,'BY=',BNY*BABS,'BZ=',BNZ*BABS
!	 WRITE(6,6001) YN(0,IT),RL,ZL,PHIL,YN(1,IT), &
!                      YN(2,IT),YN(3,IT),YN(6,IT)*VC/OMG,YN(7,IT), &
!      	               DELTA,RKPARA*VC/OMG,RKPERP*VC/OMG,RKRL, &
!      	               RNPHI_IDEI
! 6001    FORMAT(1H ,1P14E13.5)

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA*RKAP) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
      END DO
      NNSTP=NSTPLIM
!     
 11   IF(YN(7,NNSTP).LT.0.D0) THEN
         YN(7,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) &
              WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NNSTP)
      ENDIF

      RETURN
  END SUBROUTINE WRRKFT_WITHMC

!************************************************************************

  SUBROUTINE WRMODCONV_OX(IOX, Y, F, OXEFF) 

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IOX
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: F(NEQ)
    REAL(rkind),INTENT(OUT):: OXEFF
    INTEGER:: I
    REAL(rkind):: YK(3)
    REAL(rkind):: OMG,OX_K0,OX_LN,OX_Y,OX_NZOPT,OX_NZ,OX_NY,OX_N0,RHON
    REAL(rkind):: BNX0,BNY0,BNZ0,RL0,Rr_IDEI,RKPARA0,S_O_X
    REAL(rkind):: DELTAB,Y10,Y20,Y30,OX_KC,Y4_OX,Y5_OX,Y6_OX
    REAL(rkind):: Y1_OX,Y2_OX,Y3_OX,DELTA

       OMG=2.D6*PI*RF

       CALL RAMBDA_N_OX(Y(1), Y(2), Y(3), OX_LN)
       CALL REFINDEX(Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)
       OX_K0 = OMG / VC
       WRITE(6,*)'N_OPT=',OX_NZOPT,'NZ=',OX_NZ,'NY=',OX_NY
       WRITE(6,*)'K0=',OX_K0,'N=', &
                    SQRT((Y(4)**2+Y(5)**2+Y(6)**2))*(VC/OMG)	

       OXEFF = ( 2.0*(1.0+OX_Y)*((OX_NZ-OX_NZOPT)**2) + OX_NY**2 )
       OXEFF = EXP(-PI*OX_K0*OX_LN*SQRT(0.5*OX_Y)*OXEFF)
       WRITE(6,*) 'OXEFF=',OXEFF 

       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       CALL PL_PROF_OLD(RHON)	   

       BNX0 = BNX
       BNY0 = BNY
       BNZ0 = BNZ
       RL0  =SQRT(Y(1)**2+Y(2)**2)
       Rr_IDEI = RL0-RR
       RKPARA0=Y(4)*BNX0+Y(5)*BNY0+Y(6)*BNZ0
       S_O_X = 5.0D-5

       DELTAB =1.0D0
       DO I=1,1000000
          IF(I.EQ.1) THEN 
             DELTAB=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), OMG )
             Y10 = (Rr_IDEI/RL0)*Y(1)
             Y20 = (Rr_IDEI/RL0)*Y(2)
             Y30 = Y(3)
             Y10 = Y10 / SQRT(Y10**2+Y20**2+Y30**2)
             Y20 = Y20 / SQRT(Y10**2+Y20**2+Y30**2)
             Y30 = Y30 / SQRT(Y10**2+Y20**2+Y30**2)
             OX_KC = (Y(4)*Y10 + Y(5)*Y20 + Y(6)*Y30)
             Y4_OX = Y(4) - OX_KC * Y10 
             Y5_OX = Y(5) - OX_KC * Y20 
             Y6_OX = Y(6) - OX_KC * Y30
          END IF
          Y1_OX = Y(1) - IOX * S_O_X * Y10*Y(1)
          Y2_OX = Y(2) - IOX * S_O_X * Y20*Y(2)
          Y3_OX = Y(3) - IOX * S_O_X * Y30*Y(3)

          DELTA=DISPXR( Y1_OX, Y2_OX, Y3_OX, Y4_OX, Y5_OX, Y6_OX, OMG )
		
          IF ( DELTA*DELTAB.LT.0D0 ) THEN
             Y(1) =Y1_OX
             Y(2) =Y2_OX
             Y(3) =Y3_OX
             Y(4) =Y4_OX
             Y(5) =Y5_OX
             Y(6) =Y6_OX
             EXIT
          END IF
          IOX=I
       END DO
					    
       DO I =1,NEQ
          F(I) = Y(I)
       END DO

       RETURN
  END SUBROUTINE WRMODCONV_OX

!************************************************************************

  SUBROUTINE  REFINDEX(Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: OX_Y,OX_NZOPT,OX_NZ,OX_NY
    REAL(rkind):: OMG,RHON,OMG_C_OX,RL,Rr_IDEI
    REAL(rkind):: Y10,Y20,Y30,OX_NX,OX_N2

       OMG=2.D6*PI*RF 
		
       CALL PL_MAG_OLD(Y(1), Y(2), Y(3), RHON)
       OMG_C_OX = BABS*AEE/(AME)
       OX_Y = OMG_C_OX / OMG
       OX_NZOPT = SQRT(OX_Y/(1.D0+OX_Y))
       OX_NZ = (Y(4)*BNX + Y(5)*BNY + Y(6)*BNZ)*VC/OMG
	
       RL  =SQRT(Y(1)**2+Y(2)**2)
       Rr_IDEI = RL-RR
       Y10 = (Rr_IDEI/RL)*Y(1)
       Y20 = (Rr_IDEI/RL)*Y(2)
       Y30 = Y(3)
       Y10 = Y10 / SQRT(Y10**2+Y20**2+Y30**2)
       Y20 = Y20 / SQRT(Y10**2+Y20**2+Y30**2)
       Y30 = Y30 / SQRT(Y10**2+Y20**2+Y30**2)
       OX_NX = (Y(4)*Y10 + Y(5)*Y20 + Y(6)*Y30)*VC/OMG
		
       OX_N2 = (Y(4)**2+Y(5)**2+Y(6)**2)*(VC/OMG)*(VC/OMG)		
       OX_NY = OX_N2 - OX_NZ**2 - OX_NX**2
       IF (OX_NY.LT.0.D0) OX_NY=0.0D0
       OX_NY = SQRT(OX_NY)
			
       RETURN
  END SUBROUTINE REFINDEX

!************************************************************************

  SUBROUTINE  RAMBDA_N_OX(X,Y,Z, OX_LN)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y,Z
    REAL(rkind),INTENT(OUT):: OX_LN
    REAL(rkind):: RHON,OX_NE,RL0,Rr_IDEI
    REAL(rkind):: OX_X0,OX_Y0,OX_Z0,D_OX_X0,D_OX_Y0,D_OX_Z0,OX_NE_P,OX_NE_M

      CALL PL_MAG_OLD(X,Y,Z, RHON)
      CALL PL_PROF_OLD(RHON)	  
      OX_NE = RN(1) 
	  
      RL0  =SQRT(X**2+Y**2)
      Rr_IDEI = RL0-RR

      OX_X0 = (Rr_IDEI/RL0)*X
      OX_Y0 = (Rr_IDEI/RL0)*Y
      OX_Z0 = Z
      D_OX_X0 = - (DELDER) * OX_X0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
      D_OX_Y0 = - (DELDER) * OX_Y0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
      D_OX_Z0 = - (DELDER) * OX_Z0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
	  
!	  WRITE(6,*)'X=',X,'Y=',Y,'Z=',Z
!	  WRITE(6,*)'DX=',D_OX_X0,'DY=',D_OX_Y0,'DZ=',D_OX_Z0
	  
      CALL PL_MAG_OLD(X+D_OX_X0, Y+D_OX_Y0, Z+D_OX_Z0, RHON)
      CALL PL_PROF_OLD(RHON)	  
      OX_NE_P = RN(1)	  
      CALL PL_MAG_OLD(X-D_OX_X0, Y-D_OX_Y0, Z-D_OX_Z0, RHON)
      CALL PL_PROF_OLD(RHON)	  
      OX_NE_M = RN(1) 
	  
!	  WRITE(6,*) 'OX_NE_P(E18)=',OX_NE_P,'OX_NE_M(E18)=', OX_NE_M
	  
      OX_LN = OX_NE / ( (OX_NE_P-OX_NE_M)/(2.0D0*DELDER) )
!	  WRITE(6,*) 'NE(E18)=',OX_NE,'OX_LN=',OX_LN

      RETURN
  END SUBROUTINE RAMBDA_N_OX

!************************************************************************

  SUBROUTINE WRSYMP(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: F(NEQ)
    INTEGER:: NSTPLIM,NLPMAX,NSTP,I,ID,NLP,IERR
    REAL(rkind):: EPS,X,RL,PHIL,ZL,RKRL,OMG,RHON,OMG_C,WP2,ERROR,DELTA,Y7

    NSTPLIM=INT(SMAX/DELS)
    NLPMAX=10
    EPS=1.D-6

    NSTP=0
    X=0.D0
    YN(0,NSTP)=X
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0

    DO NSTP = 1,NSTPLIM
       Y7=Y(7)
       CALL SYMPLECTIC(Y,DELS,WRFDRVR,6,NLPMAX,EPS,NLP,ERROR,IERR)
       CALL WRFDRV(0.D0,Y,F)
       X=X+DELS

       YN(0,NSTP)=X
       DO I=1,6
          YN(I,NSTP)=Y(I)
       ENDDO
       Y(7)=Y(7)+F(7)*DELS
       IF(Y(7).GT.0.D0) THEN
          YN(7,NSTP)=Y(7)
          YN(8,NSTP)=-F(7)*DELS
       ELSE
          YN(7,NSTP)=0.D0
          YN(8,NSTP)=Y(7)
       ENDIF

       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          RL  =Y(1)
          PHIL=ASIN(Y(2)/(2.D0*PI*RR))
          ZL  =Y(3)
          RKRL=Y(4)
       ELSE
          RL  =SQRT(Y(1)**2+Y(2)**2)
          PHIL=ATAN2(Y(2),Y(1))
          ZL  =Y(3)
          RKRL=(Y(4)*Y(1)+Y(5)*Y(2))/RL
       ENDIF

       OMG=2.D6*PI*RF
       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       CALL PL_PROF_OLD(RHON)	  
       OMG_C = BABS*AEE/(AME)
       WP2=RN(1)*1.D20*AEE*AEE/(EPS0*AMP*PA(1))
       wp2=sqrt(WP2)

       DELTA=DISPXR(Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),OMG)

       IF(MDLWRW.GE.1) THEN
          ID=0
          SELECT CASE(MDLWRW)
          CASE(1)
             ID=1
          CASE(2)
             IF(MOD(NSTP-1,10).EQ.0) ID=1
          CASE(3)
             IF(MOD(NSTP-1,100).EQ.0) ID=1
          CASE(4)
             IF(MOD(NSTP-1,1000).EQ.0) ID=1
          CASE(5)
             IF(MOD(NSTP-1,10000).EQ.0) ID=1
          END SELECT
          IF(ID.EQ.1) &
               WRITE(6,'(1P7E11.3)') X,RL,PHIL,ZL,RKRL,Y(7),YN(8,NSTP)
       ENDIF

       IF(Y(7).LT.UUMIN) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       IF(RHON.GT.RB/RA*1.2D0) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
    ENDDO
    NNSTP=NSTPLIM

11  IF(YN(7,NNSTP).LT.0.D0) THEN
       YN(7,NNSTP)=0.D0
    ENDIF
    IF(MDLWRW.GE.1) THEN
       IF(ID.EQ.0) &
            WRITE(6,'(1P7E11.3)') X,RL,PHIL,ZL,RKRL,Y(7),YN(8,NNSTP)
    ENDIF

    RETURN
  END SUBROUTINE WRSYMP

!************************************************************************

  SUBROUTINE WRFDRV(X,Y,F) 

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y(7)
    REAL(rkind),INTENT(OUT):: F(7)
    REAL(rkind):: VV,TT,XP,YP,ZP,RKXP,RKYP,RKZP,UU,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ,VDU

      VV=DELDER
      TT=DELDER

      XP=Y(1)
      YP=Y(2)
      ZP=Y(3)
      RKXP=Y(4)
      RKYP=Y(5)
      RKZP=Y(6)
      UU=Y(7)
      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      RXP=MAX(ABS(XP)*VV,TT)
      RYP=MAX(ABS(YP)*VV,TT)
      RZP=MAX(ABS(ZP)*VV,TT)
      RRKXP=MAX(ABS(RKXP)*VV,TT)
      RRKYP=MAX(ABS(RKYP)*VV,TT)
      RRKZP=MAX(ABS(RKZP)*VV,TT)

      DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
      DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
      DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
      DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
      DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
      DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

      IF(DOMG.GT.0.D0) THEN
         DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
      ELSE
         DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
      ENDIF

      VX   =-DKXP/DS
      VY   =-DKYP/DS
      VZ   =-DKZP/DS

      VKX  = DXP/DS
      VKY  = DYP/DS
      VKZ  = DZP/DS

      VDU  =-2.D0*ABS(DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)/DS)

      F(1)=VX
      F(2)=VY
      F(3)=VZ
      F(4)=VKX
      F(5)=VKY
      F(6)=VKZ
      IF(UU.LT.0.D0) THEN
         F(7)=0.D0
      ELSE
         F(7)=VDU*UU 
      ENDIF

!      WRITE(6,'(A,1P6E12.4)') 'XY:',X,Y(1),Y(4),Y(5),Y(6),Y(7)
!      WRITE(6,'(A,1P6E12.4)') 'F :',F(1),F(2),F(3),F(4),F(5),F(6)
      RETURN
  END SUBROUTINE WRFDRV

!************************************************************************

  SUBROUTINE WRFDRVR(Y,F) 

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y(6)
    REAL(rkind),INTENT(OUT):: F(6)
    REAL(rkind):: VV,TT,XP,YP,ZP,RKXP,RKYP,RKZP,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ

      VV=DELDER
      TT=DELDER

      XP=Y(1)
      YP=Y(2)
      ZP=Y(3)
      RKXP=Y(4)
      RKYP=Y(5)
      RKZP=Y(6)
      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      RXP=MAX(ABS(XP)*VV,TT)
      RYP=MAX(ABS(YP)*VV,TT)
      RZP=MAX(ABS(ZP)*VV,TT)
      RRKXP=MAX(ABS(RKXP)*VV,TT)
      RRKYP=MAX(ABS(RKYP)*VV,TT)
      RRKZP=MAX(ABS(RKZP)*VV,TT)

      DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
      DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
      DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
      DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
      DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
      DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

      IF(DOMG.GT.0.D0) THEN
         DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
      ELSE
         DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
      ENDIF

      VX   =-DKXP/DS
      VY   =-DKYP/DS
      VZ   =-DKZP/DS

      VKX  = DXP/DS
      VKY  = DYP/DS
      VKZ  = DZP/DS

      F(1)=VX
      F(2)=VY
      F(3)=VZ
      F(4)=VKX
      F(5)=VKY
      F(6)=VKZ
      RETURN
  END SUBROUTINE WRFDRVR

!************************************************************************

  SUBROUTINE WRNWTN(RKRI,RKZI,RKPHII,IERR)

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: RKRI
    REAL(rkind),INTENT(IN):: RKZI,RKPHII
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ICOUNT
    REAL(rkind):: OMG,S,T,RKR

      IERR=0
      OMG=2.D6*PI*RF
      ICOUNT=0
 10   CONTINUE
      ICOUNT=ICOUNT+1

      S= DISPFN(RKRI,      RKPHII,RKZI,RPI,ZPI,PHII,OMG)
      T=(DISPFN(RKRI+DELKR,RKPHII,RKZI,RPI,ZPI,PHII,OMG) &
        -DISPFN(RKRI-DELKR,RKPHII,RKZI,RPI,ZPI,PHII,OMG))/(2*DELKR)

      RKR=RKRI-S/T
      IF(MDLWRW.EQ.1) &
           WRITE(6,'(A,1P3E12.4)') 'RKR,RKRI,-S/T=',RKR,RKRI,-S/T

      IF(ABS((RKR-RKRI)/RKRI).LE.EPSNW) GOTO 9000
      IF(ICOUNT.GT.LMAXNW) GOTO 8000
      RKRI=RKR
      GOTO 10

 8000 WRITE(6,*) ' WRNWTN: DOES NOT CONVERGE'
      IERR=1000
 9000 RETURN   
  END SUBROUTINE WRNWTN

!************************************************************************

  FUNCTION DISPFN(RKR,RKPHI,RKZ,RP,ZP,PHI,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp,ONLY: cfdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RKR,RKPHI,RKZ,RP,ZP,PHI,OMG
    REAL(rkind):: DISPFN
    INTEGER:: MODELPS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

    CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
    IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
       CKX=DCMPLX(RKR,0.D0)
       CKY=DCMPLX(RKPHI,0.D0)
       CKZ=DCMPLX(RKZ,0.D0)
       X=RP
       Y=2.D0*PI*RR*SIN(PHI)
       Z=ZP
    ELSEIF(MODELG.EQ.11) THEN
       CKX=DCMPLX(RKR,0.D0)
       CKY=DCMPLX(RKPHI,0.D0)
       CKZ=DCMPLX(RKZ,0.D0)
       X=RP
       Y=PHI
       Z=ZP
    ELSE
       CKX=DCMPLX(RKR*COS(PHI)-RKPHI*SIN(PHI),0.D0)
       CKY=DCMPLX(RKR*SIN(PHI)+RKPHI*COS(PHI),0.D0)
       CKZ=DCMPLX(RKZ,0.D0)
       X=RP*COS(PHI)
       Y=RP*SIN(PHI)
       Z=ZP
    ENDIF
        
    DO NS=1,NSMAX
       MODELPS(NS)=MODELP(NS)
       IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
       IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
       IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
    ENDDO

    CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

    CWW=2.D0*PI*1.D6*CRF
    DO NS=1,NSMAX
       CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
       CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
    ENDDO

    DO NS=1,NSMAX
       MODELP(NS)=MODELPS(NS)
    ENDDO

    DISPFN=DBLE(CF)
    RETURN
  END FUNCTION DISPFN

!************************************************************************

  FUNCTION DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp,ONLY: cfdispr
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    REAL(rkind):: DISPXR
    INTEGER:: MODELPS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
!            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
         IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
         IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
      ENDDO

      CF=CFDISPR(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         IF(NSDP(NS).EQ.1) THEN
            CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
         ENDIF
      ENDDO

      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
      ENDDO

      DISPXR=DBLE(CF)
      RETURN
  END FUNCTION DISPXR

!***********************************************************************

  FUNCTION DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp,ONLY: cfdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    REAL(rkind):: DISPXI
    INTEGER:: MODELPS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
         IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
         IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
      ENDDO

      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         IF(NSDP(NS).EQ.1) THEN
            CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
         ENDIF
      ENDDO

      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
      ENDDO

      DISPXI=DIMAG(CF)
      RETURN
  END FUNCTION DISPXI

!***********************************************************************
!     calculate E,K
!***********************************************************************

  SUBROUTINE WRCALE(YN,NSTPMAX_L,NRAY)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(IN):: NSTPMAX_L,NRAY
    COMPLEX(rkind):: CDET(3,3),CDETP(3,3),CDETM(3,3),CDETD(3,3)
    INTEGER:: NSTP,J,I
    REAL(rkind):: OMG,X1,Y1,Z1,EQ,VV,TT,ROMG,UE2,UE,RHON,EA
    COMPLEX(rkind):: CRF,CKX1,CKY1,CKZ1,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ,CRFP,CRFM,CUE2

      OMG=2.D6*PI*RF

      DO NSTP=0,NSTPMAX_L
         CRF =DCMPLX(OMG/(2.D6*PI),0.D0)
         X1  =YN(1,NSTP)
         Y1  =YN(2,NSTP)
         Z1  =YN(3,NSTP)
         CKX1=DCMPLX(YN(4,NSTP),0.D0)
         CKY1=DCMPLX(YN(5,NSTP),0.D0)
         CKZ1=DCMPLX(YN(6,NSTP),0.D0)
         CALL DP_DISP(CRF,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDET)
         CE1=CDET(2,2)*CDET(1,3)-CDET(1,2)*CDET(2,3)
         CE2=CDET(1,1)*CDET(2,3)-CDET(2,1)*CDET(1,3)
         CE3=CDET(2,2)*CDET(1,1)-CDET(1,2)*CDET(2,1)
         CE4=CDET(1,3)*CDET(2,1)-CDET(2,3)*CDET(1,1)
         IF(CE2.EQ.0.OR.CE4.EQ.0)THEN
            CEXY=0
            CEZY=0
         ELSE
            CEXY=CE1/CE2
            CEZY=CE3/CE4
         ENDIF

         EA=SQRT(ABS(CEXY)**2+1.D0+ABS(CEZY)**2)

         CUEX=CEXY/EA
         CUEY=1.D0/EA
         CUEZ=CEZY/EA

         VV=DELDER
         TT=DELDER

         ROMG=MAX(ABS(OMG)*VV,TT)
         CRFP =DCMPLX((OMG+ROMG)/(2.D6*PI),0.D0)
         CALL DP_DISP(CRFP,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETP)
         CRFM =DCMPLX((OMG-ROMG)/(2.D6*PI),0.D0)
         CALL DP_DISP(CRFM,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETM)
         
         DO J=1,3
         DO I=1,3
            CDETD(I,J)=OMG*(CDETP(I,J)-CDETM(I,J))/(2.D0*ROMG) &
                      +CDET(I,J)
         ENDDO
         ENDDO
         CUE2=DCONJG(CUEX)*(CDETD(1,1)*CUEX &
                           +CDETD(1,2)*CUEY &
                           +CDETD(1,3)*CUEZ) &
             +DCONJG(CUEY)*(CDETD(2,1)*CUEX &
                           +CDETD(2,2)*CUEY &
                           +CDETD(2,3)*CUEZ) &
             +DCONJG(CUEZ)*(CDETD(3,1)*CUEX &
                           +CDETD(3,2)*CUEY &
                           +CDETD(3,3)*CUEZ)
         UE2=DBLE(CUE2)
         UE=SQRT(ABS(YN(7,NSTP)/UE2))

         CEXS(NSTP,NRAY)=UE*CUEX
         CEYS(NSTP,NRAY)=UE*CUEY
         CEZS(NSTP,NRAY)=UE*CUEZ
         RKXS(NSTP,NRAY)=DBLE(CKX1)
         RKYS(NSTP,NRAY)=DBLE(CKY1)
         RKZS(NSTP,NRAY)=DBLE(CKZ1)
         RXS(NSTP,NRAY)=X1
         RYS(NSTP,NRAY)=Y1
         RZS(NSTP,NRAY)=Z1

         CALL pl_mag_old(X1,Y1,Z1,RHON)
         BNXS(NSTP,NRAY)=BNX
         BNYS(NSTP,NRAY)=BNY
         BNZS(NSTP,NRAY)=BNZ
         BABSS(NSTP,NRAY)=BABS
      ENDDO

      RETURN
  END SUBROUTINE WRCALE

!********************************** CAL K ********************************

  SUBROUTINE WRCALK(NSTP,NRAY,RKPARA,RKPERP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NSTP,NRAY
    REAL(rkind),INTENT(OUT):: RKPARA,RKPERP
    REAL(rkind):: X,Y,Z,RHON,RKX,RKY,RKZ

      X=RAYS(1,NSTP,NRAY)
      Y=RAYS(2,NSTP,NRAY)
      Z=RAYS(3,NSTP,NRAY)

      CALL PL_MAG_OLD(X,Y,Z,RHON)

      RKX=RAYS(4,NSTP,NRAY)
      RKY=RAYS(5,NSTP,NRAY)
      RKZ=RAYS(6,NSTP,NRAY)

      RKPARA=RKX*BNX+RKY*BNY+RKZ*BNZ
      RKPERP=SQRT(RKX**2+RKY**2+RKZ**2-RKPARA**2)
      RETURN
  END SUBROUTINE WRCALK

!     ***** ABSORBED POWER PROFILE*****

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

      DRHO=1.D0/NRDIVMAX
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


!_ZHENYA_2010_11_1
!************************************************************************

  SUBROUTINE WRMODNWTN(Y_zh, YK) 

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y_zh(NEQ)
    REAL(rkind),INTENT(OUT):: YK(3)
    REAL(rkind):: Y(NEQ)
    INTEGER:: I,IMODNWTN
    REAL(rkind):: OMG,DELTA,XP,YP,ZP,VV,TT
    REAL(rkind):: RKXP,RKYP,RKZP,RRKXP,RRKYP,RRKZP,DKXP,DKYP,DKZP
    REAL(rkind):: DS2,FAC_NWTN,DDELTA

    Y(1:NEQ)=Y_zh(1:NEQ)
    OMG=2.D6*PI*RF
    DELTA=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), OMG )
    XP=Y(1)
    YP=Y(2)
    ZP=Y(3)
    VV=DELDER
    TT=DELDER

    DO I=1,LMAXNW
       RKXP=Y(4)
       RKYP=Y(5)
       RKZP=Y(6)
       RRKXP=MAX(ABS(RKXP)*VV,TT)
       RRKYP=MAX(ABS(RKYP)*VV,TT)
       RRKZP=MAX(ABS(RKZP)*VV,TT)

       DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
            -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
       DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
            -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
       DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
            -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

       DS2 = DKXP**2+DKYP**2+DKZP**2
       FAC_NWTN = 10.D0
       DO IMODNWTN=1,10
          FAC_NWTN = FAC_NWTN/10.D0
          Y(4) = RKXP - (DELTA*DKXP)/DS2*FAC_NWTN
          Y(5) = RKYP - (DELTA*DKYP)/DS2*FAC_NWTN
          Y(6) = RKZP - (DELTA*DKZP)/DS2*FAC_NWTN
          DDELTA = DISPXR(XP,YP,ZP,Y(4),Y(5),Y(6),OMG)
          IF (ABS(DDELTA) .LT. ABS(DELTA)) EXIT
       END DO
       IF (ABS(DDELTA).LT.1.0D-6) EXIT
       DELTA = DDELTA
    END DO

    DO I =1,3
       YK(I) = Y(I+3)
    END DO

    RETURN
  END SUBROUTINE WRMODNWTN
END MODULE wrcalc
