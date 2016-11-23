C     $Id$
C     
C***********************************************************************
C
      SUBROUTINE WRCALC
C
      USE plcomm
      INCLUDE 'wrcomm.inc'
C
      DIMENSION Y(NEQ)
      REAL*4 TIME1,TIME2
      DATA MODEW/0/
C 
      CALL GUTIME(TIME1)
C
C      CALL DPCHEK(IERR)
C      IF(IERR.NE.0) RETURN
C
      IF(MDLWRI.EQ.0) THEN
         WRITE(6,*) 
     &   '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
         WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
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
            WRITE(6,*) 
     &      '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                         RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
         ELSEIF(MDLWRI.EQ.2) THEN
            WRITE(6,*) 
     &      '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                         RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
         ENDIF
      ELSEIF(MDLWRI.EQ.11) THEN
         WRITE(6,*) 
     &   '# default values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
         WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
         PHII=0.D0
      ENDIF
c
c     --- eliminate disp factor for same z/a species ---
C
      DO NS=1,NSMAX
         NSDP(NS)=1
         DO NSS=1,NS-1
            ZA1=PZ(NS)/PA(NS)
            ZA2=PZ(NSS)/PA(NSS)
            IF(ABS(ZA1-ZA2).LE.1.D-8) NSDP(NS)=0
         ENDDO
      ENDDO
C
C     --- Each ray tracing ---
C
      DO NRAY=1,NRAYMX
C
         IF(MDLWRI.LT.100) THEN
C
    1       WRITE(6,*) '# NRAY = ',NRAY
            IF(MDLWRI.EQ.0) THEN
               READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
               WRITE(6,*) 
     &         '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
               WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
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
               READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
               WRITE(6,*) 
     &         '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
               WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
               SINP2=SIN(ANGZ*PI/180.D0)**2
               SINT2=SIN(ANGPH*PI/180.D0)**2
               RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
               RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
               IF(ANGZ.LT.0.D0) RNZI=-RNZI
               IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
            ELSEIF(MDLWRI.EQ.2) THEN
               READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
               WRITE(6,*) 
     &         '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
               WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
               SINP2=SIN(ANGZ*PI/180.D0)**2
               SINT2=SIN(ANGPH*PI/180.D0)**2
               RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
               RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
               IF(ANGZ.LT.0.D0) RNZI=-RNZI
               IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
               WRITE(6,*) 'XX MDLWRI=2 IS NOT SUPPORTED YET.'
               GOTO 1
            ELSEIF(MDLWRI.EQ.11) THEN
               READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
               WRITE(6,*) 
     &         '# initial values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
               WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
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
               WRITE(6,*) 
     &         '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
               WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
            ELSE
               WRITE(6,*) 
     &         '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
               WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
               SINP2=SIN(ANGZ*PI/180.D0)**2
               SINT2=SIN(ANGPH*PI/180.D0)**2
               RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
               RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
               IF(ANGZ.LT.0.D0) RNZI=-RNZI
               IF(ANGPH.LT.0.D0) RNPHII=-RNPHII

               IF(RNZI.EQ.0.D0.AND.RNPHII.EQ.0.D0) THEN
                  ANGZ_=0.D0
                  ANGPH_=0.D0
               ELSE
                  SINP2=RNZI**4  /(RNZI**2+RNPHII**2-RNZI**2*RNPHII**2)
                  SINT2=RNPHII**4/(RNZI**2+RNPHII**2-RNZI**2*RNPHII**2)
                  ANGZ_= 180.D0/PI*ASIN(SQRT(SINP2))
                  ANGPH_=180.D0/PI*ASIN(SQRT(SINT2))
                  IF(RNZI  .LT.0.D0) ANGZ_= -ANGZ_
                  IF(RNPHII.LT.0.D0) ANGPH_=-ANGPH_
               ENDIF
               WRITE(6,'(A,1P2E12.4)') 'ANG:  ',ANGZ,ANGPH
               WRITE(6,'(A,1P2E12.4)') 'RNI:  ',RNZI,RNPHII
               WRITE(6,'(A,1P2E12.4)') 'ANG_: ',ANGZ_,ANGPH_
            ENDIF
         ENDIF
C
         RAYIN(1,NRAY)=RF
         RAYIN(2,NRAY)=RPI
         RAYIN(3,NRAY)=ZPI
         RAYIN(4,NRAY)=PHII
         RAYIN(5,NRAY)=RKR0
         RAYIN(6,NRAY)=ANGZ
         RAYIN(7,NRAY)=ANGPH
         RAYIN(8,NRAY)=UUI
C          
         RKRI  = RKR0
         RKZI  =2.D6*PI*RF*RNZI  /VC
         RKPHII=2.D6*PI*RF*RNPHII/VC
         CALL WRNWTN(RKRI,RKZI,RKPHII,IERR)
         IF(IERR.NE.0) GOTO 1200
C
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
            CALL WRRKFT(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSEIF(MDLWRQ.EQ.1) THEN
            CALL WRRKFT_WITHD0(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSEIF(MDLWRQ.EQ.2) THEN
            CALL WRRKFT_WITHMC(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSEIF(MDLWRQ.EQ.3) THEN
            CALL WRRKFT_ODE(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSEIF(MDLWRQ.EQ.4) THEN
            CALL WRSYMP(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSE
            WRITE(6,*) 'XX WRCALC: unknown MDLWRQ =', MDLWRQ
         ENDIF
         CALL WRCALE(RAYS(0,0,NRAY),NITMAX(NRAY),NRAY)
         WRITE(6,'(A,1PE12.4,A,1PE12.4)') 
     &        '    RKRI=  ',RKRI,  '  PABS/PIN=',
     &        1.D0-RAYS(7,NITMAX(NRAY),NRAY)
 1200    CONTINUE
      ENDDO
C
      CALL WRAPWR
C
      CALL GUTIME(TIME2)
      WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'
C
 9000 RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT(Y,YN,NIT)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      INCLUDE 'wrcomm.inc'  
C
      EXTERNAL WRFDRV
      DIMENSION Y(NEQ),YM(NEQ),YN(0:NEQ,0:NITM),WORK(2,NEQ)
CIDEI
	  DIMENSION YK(3), F(NEQ)
CIDEI	   
      X0 = 0.D0
      XE = DELS
      INIT = 1
      ITMAX=INT(SMAX/DELS)
      IT=0
      YN(0,IT)=X0
      DO I=1,7
         YN(I,IT)=Y(I)
      ENDDO
      YN(8,IT)=0.D0
      OMG=2.D6*PI*RF
      IOX=0
C
      DO 10 IT = 1,ITMAX
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)

         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=YM(I)
         ENDDO
         YN(8,IT)=Y7-YM(7)

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
C
         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(IT-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(IT-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(IT-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(IT-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) 
     &           WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
         ENDIF

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF
         IF(MODELG.LE.10) THEN
            CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
            IF(RHON.GT.RB/RA*RKAP) THEN
               NIT = IT
               GOTO 11
            ENDIF
         ENDIF
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) 
     &        WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT_WITHD0(Y,YN,NIT)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY:PL_MAG_OLD
      INCLUDE 'wrcomm.inc'  
C
      EXTERNAL WRFDRV
      DIMENSION Y(NEQ),YM(NEQ),YN(0:NEQ,0:NITM),WORK(2,NEQ)
CIDEI
	  DIMENSION YK(3), F(NEQ)
CIDEI	   
      X0 = 0.D0
      XE = DELS
      INIT = 1
      ITMAX=INT(SMAX/DELS)
      IT=0
      YN(0,IT)=X0
      DO I=1,7
         YN(I,IT)=Y(I)
      ENDDO
      YN(8,IT)=0.D0
      OMG=2.D6*PI*RF
      IOX=0
C
      DO 10 IT = 1,ITMAX
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
C         write(6,'(A,I5,1p2E12.4)') 
C     &        'WRRKFT_WITHD0: IT,Y,YM=',IT,Y(1),YM(1)

         delta=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),OMG)
CBGNIDEI_WRMOD		 
         IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,YK)
            DO I=1,3
               YM(I+3) = YK(I)
            END DO
         END IF 
CENDIDEI		 

         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=YM(I)
         ENDDO
         YN(8,IT)=Y7-YM(7)

         CALL PL_MAG_OLD(YM(1),YM(2),YM(3),PSIN)
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

CIDEI	 WRITE(6,*) 'BX=',BNX*BABS,'BY=',BNY*BABS,'BZ=',BNZ*BABS
C	 WRITE(6,6001) YN(0,IT),RL,ZL,PHIL,YN(1,IT),
C     &                YN(2,IT),YN(3,IT),YN(6,IT)*VC/OMG,YN(7,IT),
C     &	               DELTA,RKPARA*VC/OMG,RKPERP*VC/OMG,RKRL,
C     &	               RNPHI_IDEI
C 6001    FORMAT(1H ,1P14E13.5)

C
         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(IT-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(IT-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(IT-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(IT-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) 
     &           WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
         ENDIF

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF
         IF(MODELG.LE.10) THEN
            CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
            IF(RHON.GT.RB/RA*RKAP) THEN
               NIT = IT
               GOTO 11
            ENDIF        
         END IF
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) 
     &        WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NIT)
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT_ODE(Y,YN,NIT)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      INCLUDE 'wrcomm.inc'  
C
      EXTERNAL WRFDRV
      DIMENSION Y(NEQ),YM(NEQ),YN(0:NEQ,0:NITM),WORK(2,NEQ),YK(3)
C
      X0 = 0.D0
      XE = DELS     
      ITMAX=INT(SMAX/DELS)
      IT=0
      YN(0,IT)=X0
      DO I=1,7
         YN(I,IT)=Y(I)
      ENDDO
      YN(8,IT)=0.D0
C
      DO 10 IT = 1,ITMAX
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
C
         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=YM(I)
         ENDDO
         IF(YM(7).GT.0.D0) THEN
            YN(8,IT)=Y7-YM(7)
         ELSE
            YN(8,IT)=Y7
         ENDIF
C
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
C
         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(IT-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(IT-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(IT-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(IT-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) 
     &           WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
         ENDIF
C
         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA*RKAP) THEN
            NIT = IT
            GOTO 11
         ENDIF
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) 
     &        WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT_WITHMC(Y,YN,NIT)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      INCLUDE 'wrcomm.inc'  
C
      EXTERNAL WRFDRV
      DIMENSION Y(NEQ),YM(NEQ),YN(0:NEQ,0:NITM),WORK(2,NEQ)
CIDEI
	  DIMENSION YK(3), F(NEQ)
CIDEI	   
      X0 = 0.D0
      XE = DELS
      INIT = 1
      ITMAX=INT(SMAX/DELS)
      IT=0
      YN(0,IT)=X0
      DO I=1,7
         YN(I,IT)=Y(I)
      ENDDO
      YN(8,IT)=0.D0
      OMG=2.D6*PI*RF
      IOX=0
C
      DO 10 IT = 1,ITMAX
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         delta=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),OMG)
CBGNIDEI_WRMOD		 
         IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,YK)
            DO I=1,3
               YM(I+3) = YK(I)
            END DO
         END IF 
CENDIDEI		 

CBGNIDEI_MODCONV_OX
         RL  =SQRT(YM(1)**2+YM(2)**2)
         RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         IF ( RKRL.GE.0D0.AND.IOX.EQ.0 ) THEN
            CALL WRMODCOV_OX(IOX,YM,F)
            IF(IOX.GE.100000) THEN
               WRITE(6,*) 'ERROR in WRMODCON_OX routine IOX=',IOX
            ELSE 
               DO I=1,NEQ
                  YM(I) = F(I)
               END DO
            END IF
         ENDIF
C
         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=YM(I)
         ENDDO
         IF(YM(7).GT.0.D0) THEN
            YN(8,IT)=Y7-YM(7)
         ELSE
            YN(8,IT)=Y7
         ENDIF
C
         CALL PL_MAG_OLD(YM(1),YM(2),YM(3),PSIN)
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
               IF(MOD(IT-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(IT-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(IT-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(IT-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) 
     &           WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
         ENDIF

CIDEI	 WRITE(6,*) 'BX=',BNX*BABS,'BY=',BNY*BABS,'BZ=',BNZ*BABS
C	 WRITE(6,6001) YN(0,IT),RL,ZL,PHIL,YN(1,IT),
C     &                YN(2,IT),YN(3,IT),YN(6,IT)*VC/OMG,YN(7,IT),
C     &	               DELTA,RKPARA*VC/OMG,RKPERP*VC/OMG,RKRL,
C     &	               RNPHI_IDEI
C 6001    FORMAT(1H ,1P14E13.5)

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA*RKAP) THEN
            NIT = IT
            GOTO 11
         ENDIF
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) 
     &        WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRMODCOV_OX(IOX, Y, F, OXEFF) 
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
      INCLUDE 'wrcomm.inc' 
C
      DIMENSION Y(NEQ),F(NEQ), YK(3)
C
		OMG=2.D6*PI*RF
*
		CALL RAMBDA_N_OX(Y(1), Y(2), Y(3), OX_LN)
		CALL REFINDEX(Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)
		OX_K0 = OMG / VC
		WRITE(6,*)'N_OPT=',OX_NZOPT,'NZ=',OX_NZ,'NY=',OX_NY
		WRITE(6,*)'K0=',OX_K0,'N=', 
     &               SQRT((Y(4)**2+Y(5)**2+Y(6)**2))*(VC/OMG)	
		
		OXEFF = ( 2.0*(1.0+OX_Y)*((OX_NZ-OX_NZOPT)**2) 
     &                  + OX_NY**2 )
		OXEFF = EXP(-PI*OX_K0*OX_LN*SQRT(0.5*OX_Y)*OXEFF)
		WRITE(6,*) 'OXEFF=',OXEFF 

		CALL PL_MAG_OLD(Y(1),Y(2),Y(3),PSIN)
		CALL PL_PROF_OLD(PSIN)	   
		BNX0 = BNX
		BNY0 = BNY
		BNZ0 = BNZ
		RL0  =SQRT(Y(1)**2+Y(2)**2)
		Rr_IDEI = RL0-RR
		RKPARA0=Y(4)*BNX0+Y(5)*BNY0+Y(6)*BNZ0
C		WRITE(6,*)'RKPARA0=',RKPARA0
		S_O_X = 5.0D-5
*
CC		0.46, 1.0D-6 EFF=0.33
CC		0.66 5.0D-5 EFF=0.52

		DELTAB =1.0D0
		DO IOX=1,1000000
				IF(IOX.EQ.1) THEN 
				DELTAB=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), OMG )
CCC				WRITE(6,*) 'O=>X LOOP IOX=',IOX,'DELTAB=',DELTAB
CCC				WRITE(6,*) Y(1),Y(2),Y(3)
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
C					Y4_OX = RKPARA0 * BNX0
C					Y5_OX = RKPARA0 * BNY0
C					Y6_OX = RKPARA0 * BNZ0					
				END IF
				Y1_OX = Y(1) - IOX * S_O_X * Y10*Y(1)
				Y2_OX = Y(2) - IOX * S_O_X * Y20*Y(2)
				Y3_OX = Y(3) - IOX * S_O_X * Y30*Y(3)
*
				DELTA=DISPXR( Y1_OX, Y2_OX, Y3_OX, Y4_OX, Y5_OX, Y6_OX, OMG )
*				
C				WRITE(6,*) 'DELTA=',DELTA
CCC				WRITE(6,*) Y1_OX,Y2_OX,Y3_OX,Y4_OX,Y5_OX,Y6_OX
				IF ( DELTA*DELTAB.LT.0D0 ) THEN
C				IF ( DELTA.LT.0.D0.AND.ABS(DELTA).LT.1.0D-6) THEN
C				IF ( (WPE2/OMG/OMG).GT.1.0D0 ) THEN
					Y(1) =Y1_OX
					Y(2) =Y2_OX
					Y(3) =Y3_OX
					Y(4) =Y4_OX
					Y(5) =Y5_OX
					Y(6) =Y6_OX
C					CALL WRMODNWTN(Y,YK)
C					DO I=1,3
C						Y(I+3) = YK(I)
C					END DO
					EXIT
				END IF
		END DO
						    
		DO I =1,NEQ
			F(I) = Y(I)
		END DO  

      RETURN
      END

C
C************************************************************************
C
      SUBROUTINE  REFINDEX(Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      INCLUDE 'wrcomm.inc'

	  DIMENSION Y(NEQ)

		OMG=2.D6*PI*RF 
		
		CALL PL_MAG_OLD(Y(1), Y(2), Y(3), PSIN)
		OMG_C_OX = BABS*AEE/(AME)
		OX_Y = OMG_C_OX / OMG
		OX_NZOPT = SQRT(OX_Y/(1.D0+OX_Y))
		OX_NZ = (Y(4)*BNX + Y(5)*BNY + Y(6)*BNZ)*VC/OMG
C		WRITE(6,*)'OX_NZOPT =',OX_NZOPT,'OX_NZ=',OX_NZ
	
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
      END

C
C************************************************************************
C
      SUBROUTINE  RAMBDA_N_OX(X,Y,Z, OX_LN)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
      INCLUDE 'wrcomm.inc'
C
	  CALL PL_MAG_OLD(X,Y,Z, PSIN)
	  CALL PL_PROF_OLD(PSIN)	  
	  OX_NE = RN(1) 
	  
	  RL0  =SQRT(X**2+Y**2)
	  Rr_IDEI = RL0-RR

	  OX_X0 = (Rr_IDEI/RL0)*X
	  OX_Y0 = (Rr_IDEI/RL0)*Y
	  OX_Z0 = Z
	  D_OX_X0 = - (DELDER) * OX_X0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
	  D_OX_Y0 = - (DELDER) * OX_Y0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
	  D_OX_Z0 = - (DELDER) * OX_Z0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
	  
C	  WRITE(6,*)'X=',X,'Y=',Y,'Z=',Z
C	  WRITE(6,*)'DX=',D_OX_X0,'DY=',D_OX_Y0,'DZ=',D_OX_Z0
	  
	  CALL PL_MAG_OLD(X+D_OX_X0, Y+D_OX_Y0, Z+D_OX_Z0, PSIN)
	  CALL PL_PROF_OLD(PSIN)	  
	  OX_NE_P = RN(1)	  
	  CALL PL_MAG_OLD(X-D_OX_X0, Y-D_OX_Y0, Z-D_OX_Z0, PSIN)
	  CALL PL_PROF_OLD(PSIN)	  
	  OX_NE_M = RN(1) 
	  
C	  WRITE(6,*) 'OX_NE_P(E18)=',OX_NE_P,'OX_NE_M(E18)=', OX_NE_M
	  
	  OX_LN = OX_NE / ( (OX_NE_P-OX_NE_M)/(2.0D0*DELDER) )
C	  WRITE(6,*) 'NE(E18)=',OX_NE,'OX_LN=',OX_LN

      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRSYMP(Y,YN,NIT)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
      INCLUDE 'wrcomm.inc'      
C
      EXTERNAL WRFDRV,WRFDRVR
      DIMENSION Y(NEQ),F(NEQ),YN(0:NEQ,0:NITM)
C
      ITMAX=INT(SMAX/DELS)
      NLPMAX=10
      EPS=1.D-6
C
      IT=0
      X=0.D0
      YN(0,IT)=X
      DO I=1,7
         YN(I,IT)=Y(I)
      ENDDO
      YN(8,IT)=0.D0
C
      DO IT = 1,ITMAX
         Y7=Y(7)
         CALL SYMPLECTIC(Y,DELS,WRFDRVR,6,NLPMAX,EPS,NLP,ERROR,IERR)
         CALL WRFDRV(0.D0,Y,F)
         X=X+DELS
C
         YN(0,IT)=X
         DO I=1,6
            YN(I,IT)=Y(I)
         ENDDO
         Y(7)=Y(7)+F(7)*DELS
         IF(Y(7).GT.0.D0) THEN
            YN(7,IT)=Y(7)
            YN(8,IT)=-F(7)*DELS
         ELSE
            YN(7,IT)=0.D0
            YN(8,IT)=Y(7)
         ENDIF
C
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
C
c_zhenya
         OMG=2.D6*PI*RF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),PSIN)
         CALL PL_PROF_OLD(PSIN)	  
         OMG_C = BABS*AEE/(AME)
         WP2=RN(1)*1.D20*AEE*AEE/(EPS0*AMP*PA(1))
         wp2=sqrt(WP2)
c_zhenya
         delta=DISPXR(Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),OMG)
C
         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(IT-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(IT-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(IT-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(IT-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) 
     &           WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,Y(7),YN(8,IT)
         ENDIF

         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF         
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(MODEG.GT.1) THEN
            IF(RHON.GT.RB/RA*1.2D0) THEN
               NIT = IT
               GOTO 11
            ENDIF
         ENDIF
      ENDDO
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) 
     &        WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,Y(7),YN(8,IT)
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRFDRV(X,Y,F) 
C
      USE plcomm
      INCLUDE 'wrcomm.inc'      
C
      DIMENSION Y(7),F(7)
C
      VV=DELDER
      TT=DELDER
C
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
C
      DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG)
     &     -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
      DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
      DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
      DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
      DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG)
     &     -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
      DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG)
     &     -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)
C
      IF(DOMG.GT.0.D0) THEN
         DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
      ELSE
         DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
      ENDIF
C      WRITE(21,'(1P3E12.4)') X,DOMG,DS
C      CALL GUFLSH
C
C      write(6,'(A,1P3E12.4)') 'WRFDRV: DOMG,DS,DKXP',DOMG,DS,DKXP
      VX   =-DKXP/DS
      VY   =-DKYP/DS
      VZ   =-DKZP/DS
C
      VKX  = DXP/DS
      VKY  = DYP/DS
      VKZ  = DZP/DS
C
C      VDU  =-2*DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)/DS
      VDU  =-2.D0*ABS(DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)/DS)
C
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
C
C      WRITE(6,'(A,1P6E12.4)') 'XY:',X,Y(1),Y(4),Y(5),Y(6),Y(7)
C      WRITE(6,'(A,1P6E12.4)') 'F :',F(1),F(2),F(3),F(4),F(5),F(6)
C      CALL GUFLSH
C      write(6,'(A,1P3E12.4)') 'WRFDRV: X,Y(1),F(1)=',X,Y(1),F(1)
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRFDRVR(Y,F) 
C
      USE plcomm
      INCLUDE 'wrcomm.inc'      
C
      DIMENSION Y(6),F(6)
C
      VV=DELDER
      TT=DELDER
C
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
C
      DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG)
     &     -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
      DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
      DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
      DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG)
     &     -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
      DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG)
     &     -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
      DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG)
     &     -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)
C
      IF(DOMG.GT.0.D0) THEN
         DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
      ELSE
         DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
      ENDIF
C
      VX   =-DKXP/DS
      VY   =-DKYP/DS
      VZ   =-DKZP/DS
C
      VKX  = DXP/DS
      VKY  = DYP/DS
      VKZ  = DZP/DS
C
      F(1)=VX
      F(2)=VY
      F(3)=VZ
      F(4)=VKX
      F(5)=VKY
      F(6)=VKZ
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRNWTN(RKRI,RKZI,RKPHII,IERR)
C
      USE plcomm
      INCLUDE 'wrcomm.inc'
C
      IERR=0
      OMG=2.D6*PI*RF
      ICOUNT=0
 10   CONTINUE
      ICOUNT=ICOUNT+1
C
      S= DISPFN(RKRI,      RKPHII,RKZI,RPI,ZPI,PHII,OMG)
      T=(DISPFN(RKRI+DELKR,RKPHII,RKZI,RPI,ZPI,PHII,OMG)
     &  -DISPFN(RKRI-DELKR,RKPHII,RKZI,RPI,ZPI,PHII,OMG))/(2*DELKR)
C
C      WRITE(6,*) S,T,S/T
      RKR=RKRI-S/T
      IF(MDLWRW.NE.0) 
     &     WRITE(6,'(A,1P3E12.4)') 'RKR,RKRI,-S/T=',RKR,RKRI,-S/T
C
      IF(ABS((RKR-RKRI)/RKRI).LE.EPSNW) GOTO 9000
C      WRITE(6,*) ABS((RKR-RKRI)/RKRI), RKR
      IF(ICOUNT.GT.LMAXNW) GOTO 8000
      RKRI=RKR
      GOTO 10
C
 8000 WRITE(6,*) ' WRNWTN: DOES NOT CONVERGE'
      IERR=1000
 9000 RETURN   
      END
C
C************************************************************************
C
      FUNCTION DISPFN(RKR,RKPHI,RKZ,RP,ZP,PHI,OMG)
C
      USE plcomm
      USE pllocal
      INCLUDE 'wrcomm.inc'
      DIMENSION MODELPS(NSM)
C
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
C            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
         IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
         IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
      ENDDO
C
      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)
C
      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
         CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
      ENDDO
C
      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
      ENDDO
C
      DISPFN=DBLE(CF)
      RETURN
      END
C
C************************************************************************
C
      FUNCTION DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)
C
      USE plcomm
      USE pllocal
      INCLUDE 'wrcomm.inc'
      DIMENSION MODELPS(NSM)
C
      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
C            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
         IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
         IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
      ENDDO
C
C      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)
      CF=CFDISPR(CRF,CKX,CKY,CKZ,X,Y,Z)
C
      CWW=2.D0*PI*1.D6*CRF
C      DO NS=1,NSMAX
      DO NS=1,1
         IF(NSDP(NS).EQ.1) THEN
            CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
         ENDIF
      ENDDO
C
      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
      ENDDO
C
      DISPXR=DBLE(CF)
      RETURN
      END
C
C***********************************************************************
C
      FUNCTION DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)
C
      USE plcomm
      USE pllocal
      INCLUDE 'wrcomm.inc'
C
      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
C            
      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)
C
      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         IF(NSDP(NS).EQ.1) THEN
            CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
         ENDIF
      ENDDO
C
      DISPXI=DIMAG(CF)
      RETURN
      END
C
C***********************************************************************
C     calculate E,K
C***********************************************************************
C
      SUBROUTINE WRCALE(YN,NITMX,NRAY)
C    
      USE plcomm
      INCLUDE 'wrcomm.inc'
      DIMENSION YN(0:NEQ,0:NITM),CDET(3,3)
      DIMENSION CDETP(3,3),CDETM(3,3),CDETD(3,3)
C
      OMG=2.D6*PI*RF
C
      DO NIT=0,NITMX
         CRF =DCMPLX(OMG/(2.D6*PI),0.D0)
         X1  =YN(1,NIT)
         Y1  =YN(2,NIT)
         Z1  =YN(3,NIT)
         CKX1=DCMPLX(YN(4,NIT),0.D0)
         CKY1=DCMPLX(YN(5,NIT),0.D0)
         CKZ1=DCMPLX(YN(6,NIT),0.D0)
         CALL DPDISP(CRF,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDET)
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
C
         EA=SQRT(ABS(CEXY)**2+1.D0+ABS(CEZY)**2)
C
         CUEX=CEXY/EA
         CUEY=1.D0/EA
         CUEZ=CEZY/EA
C
         VV=DELDER
         TT=DELDER
C
         ROMG=MAX(ABS(OMG)*VV,TT)
         CRFP =DCMPLX((OMG+ROMG)/(2.D6*PI),0.D0)
         CALL DPDISP(CRFP,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETP)
         CRFM =DCMPLX((OMG-ROMG)/(2.D6*PI),0.D0)
         CALL DPDISP(CRFM,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETM)
C
         
         DO J=1,3
         DO I=1,3
            CDETD(I,J)=OMG*(CDETP(I,J)-CDETM(I,J))/(2.D0*ROMG)
     &                +CDET(I,J)
         ENDDO
         ENDDO
         CUE2=DCONJG(CUEX)*(CDETD(1,1)*CUEX
     &                     +CDETD(1,2)*CUEY
     &                     +CDETD(1,3)*CUEZ)
     &       +DCONJG(CUEY)*(CDETD(2,1)*CUEX
     &                     +CDETD(2,2)*CUEY
     &                     +CDETD(2,3)*CUEZ)
     &       +DCONJG(CUEZ)*(CDETD(3,1)*CUEX
     &                     +CDETD(3,2)*CUEY
     &                     +CDETD(3,3)*CUEZ)
         UE2=DBLE(CUE2)
         UE=SQRT(ABS(YN(7,NIT)/UE2))
C
C        WRITE(6,'(I4,1P3E12.4)') NIT,YN(7,NIT),UE2,UE
C
         CEXS(NIT,NRAY)=UE*CUEX
         CEYS(NIT,NRAY)=UE*CUEY
         CEZS(NIT,NRAY)=UE*CUEZ
         RKXS(NIT,NRAY)=DBLE(CKX1)
         RKYS(NIT,NRAY)=DBLE(CKY1)
         RKZS(NIT,NRAY)=DBLE(CKZ1)
         RXS(NIT,NRAY)=X1
         RYS(NIT,NRAY)=Y1
         RZS(NIT,NRAY)=Z1
         RAYRB1(NIT,NRAY)=0.d0
         RAYRB2(NIT,NRAY)=0.d0
C
C         IF(NRAYMX.EQ.1) THEN
C            FACTA=0.D0
C         ELSE
C            FACTA=EXP(-1.D0)/(1.D0+EXP(-1.D0)*(NRAYMX))
C         ENDIF
C            FACTB=1.D0-FACTA*(NRAYMX)
C
C         IF(NRAY.EQ.1)THEN
C            CEXS(NIT,NRAY)=CEXS(NIT,NRAY)*FACTB
C            CEYS(NIT,NRAY)=CEYS(NIT,NRAY)*FACTB
C            CEZS(NIT,NRAY)=CEZS(NIT,NRAY)*FACTB
C         ELSE
C            CEXS(NIT,NRAY)=CEXS(NIT,NRAY)*FACTA
C            CEYS(NIT,NRAY)=CEYS(NIT,NRAY)*FACTA
C            CEZS(NIT,NRAY)=CEZS(NIT,NRAY)*FACTA
C         ENDIF
      ENDDO
C
      RETURN
      END
C
C********************************** CAL K ********************************
C
      SUBROUTINE WRCALK(NIT,NRAY,RKPARA,RKPERP)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      INCLUDE 'wrcomm.inc'
C
      X=RAYS(1,NIT,NRAY)
      Y=RAYS(2,NIT,NRAY)
      Z=RAYS(3,NIT,NRAY)
C
      CALL PL_MAG_OLD(X,Y,Z,RHON)
C
      RKX=RAYS(4,NIT,NRAY)
      RKY=RAYS(5,NIT,NRAY)
      RKZ=RAYS(6,NIT,NRAY)
C
      RKA=SQRT(RKX**2+RKY**2+RKZ**2)
      RKPARA=RKX*BNX+RKY*BNY+RKZ*BNZ
      RKPERP=SQRT(RKA**2-RKPARA**2)
      RETURN
      END
C
C     ***** ABSORBED POWER PROFILE*****
C
      SUBROUTINE WRAPWR
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      INCLUDE 'wrcomm.inc'
C
      PARAMETER(NRM=201)
C
      DIMENSION PWRRAY(NRM,NRAYM)
      DIMENSION PWR(NRM,NSM),AJR(NRM)
C
C     ----- CALCULATE RADIAL DEPOSITION PROFILE -----
C
C      CALL PLDATA_GETN(NRMAXPL,NSMAXPL)
C
      NRMAXPL=50
      NSMAXPL=2

      DRHO=1.D0/NRMAXPL
      DO NRAY=1,NRAYMX
         DO NR=1,NRMAXPL
            PWRRAY(NR,NRAY)=0.D0
         ENDDO
      ENDDO

      DO NRAY=1,NRAYMX
         DO IT=0,NITMAX(NRAY)-1
            XL=RAYS(1,IT,NRAY)
            YL=RAYS(2,IT,NRAY)
            ZL=RAYS(3,IT,NRAY)
            CALL PL_MAG_OLD(XL,YL,ZL,RHON1)
            IF(RHON1.LE.1.D0) THEN
               NRS1=INT(RHON1/DRHO)+1
C               WRITE(6,*) RHON1,DRHO,NRS1
               XL=RAYS(1,IT+1,NRAY)
               YL=RAYS(2,IT+1,NRAY)
               ZL=RAYS(3,IT+1,NRAY)
               CALL PL_MAG_OLD(XL,YL,ZL,RHON2)
               NRS2=INT(RHON2/DRHO)+1
C               WRITE(6,*) RHON2,DRHO,NRS2
               NDR=ABS(NRS2-NRS1)
               IF(NDR.EQ.0) THEN
                  PWRRAY(NRS1,NRAY)
     &                 =PWRRAY(NRS1,NRAY)+RAYS(8,IT+1,NRAY)
               ELSE IF(NRS1.LT.NRS2) THEN
                  SDR=(RHON2-RHON1)/DRHO
                  DELP=RAYS(8,IT+1,NRAY)/SDR
                  PWRRAY(NRS1,NRAY)=PWRRAY(NRS1,NRAY)
     &                 +(DBLE(NRS1)-RHON1/DRHO)*DELP
                  DO NR=NRS1+1,NRS2-1
                     PWRRAY(NR,NRAY)=PWRRAY(NR,NRAY)+DELP
                  END DO
                  PWRRAY(NRS2,NRAY)=PWRRAY(NRS2,NRAY)
     &                 +(RHON2/DRHO-DBLE(NRS2-1))*DELP
               ELSE
                  SDR=(RHON1-RHON2)/DRHO
                  DELP=RAYS(8,IT+1,NRAY)/SDR
                  PWRRAY(NRS2,NRAY)=PWRRAY(NRS2,NRAY)
     &                 +(DBLE(NRS2)-RHON2/DRHO)*DELP
                  DO NR=NRS2+1,NRS1-1
                     PWRRAY(NR,NRAY)=PWRRAY(NR,NRAY)+DELP
                  END DO
                  PWRRAY(NRS1,NRAY)=PWRRAY(NRS1,NRAY)
     &                 +(RHON1/DRHO-DBLE(NRS1-1))*DELP
               END IF
            END IF
         END DO
C         WRITE(6,'(5(I3,1PE12.4))') (NR,PWRRAY(NR,NRAY),NR=1,NRMAXPL)
      END DO
C
      DO NR=1,NRMAXPL
         DO NRAY=1,NRAYMX
            PWRRAY(NR,NRAY)=PWRRAY(NR,NRAY)
     &                 /(2*PI*(DBLE(NR)-0.5D0)*DRHO*DRHO)
         ENDDO
C         WRITE(6,'(5(I3,1PE12.4))') (NRZ,GPY(NRZ,NRAY),NRZ=1,NRZMAX)
      ENDDO
C
C     ----- weight of power for each ray-------
C
      DO NR=1,NRMAXPL
         PWR(NR,1)=0.D0
         DO NRAY=1,NRAYMX
            PWR(NR,1)=PWR(NR,1)+PWRRAY(NR,NRAY)
         ENDDO
      ENDDO
C      WRITE(6,'(1P5E12.4)') (PWR(NR,1),NR=1,NRMAXPL)
      DO NS=2,NSMAXPL
         DO NR=1,NRMAXPL
            PWR(NR,NS)=0.D0
         ENDDO
      ENDDO
      DO NR=1,NRMAXPL
         AJR(NR)=0.D0
      ENDDO

      PWRMAX=0.D0
      LOCMAX=0
      DO NR=1,NRMAXPL
         IF(PWR(NR,1).GT.PWRMAX) THEN
            PWRMAX=PWR(NR,1)
            LOCMAX=NR
         ENDIF
      END DO
      IF(LOCMAX.LE.1) THEN
         RHOMAX=0.D0
      ELSE IF(LOCMAX.EQ.NRMAXPL) THEN
         RHOMAX=1.D0
      ELSE
         DPWR =(PWR(LOCMAX+1,1)-PWR(LOCMAX-1,1))/(2.D0*DRHO)
         DDPWR=(PWR(LOCMAX+1,1)-2*PWR(LOCMAX,1)+PWR(LOCMAX-1,1))/DRHO**2
         RHOMAX=(LOCMAX-0.5D0)/(NRMAXPL-1.D0)-DPWR/DDPWR
         PWRMAX=PWRMAX-DPWR**2/(2.D0*DDPWR)
      ENDIF
      WRITE(6,'(A,1PE12.4,A,1PE12.4)') 
     &        '    PWRMAX=',PWRMAX,'  AT RHON =',RHOMAX
C     &            PWR(MIN(LOCMAX+1,NRMAXPL),1),PWRMAX,RHOMAX
C      WRITE(6,'(A,I5,1P5E12.4)') 
C     &     'PWR:',LOCMAX,PWR(MAX(LOCMAX-1,1),1),PWR(LOCMAX,1),
C     &            PWR(MIN(LOCMAX+1,NRMAXPL),1),PWRMAX,RHOMAX
C
C      CALL PLDATA_SETWR(1,'EC1',PWR,AJR)
C
      RETURN
      END


C_ZHENYA_2010_11_1
C************************************************************************
C
      SUBROUTINE WRMODNWTN(Y_zh, YK) 
C
      USE plcomm
      INCLUDE 'wrcomm.inc'      
C
      DIMENSION Y(NEQ),YK(3),Y_zh(NEQ)
C
		Y=Y_zh
		OMG=2.D6*PI*RF
		DELTA=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), OMG )
		XP=Y(1)
		YP=Y(2)
		ZP=Y(3)
		VV=DELDER
		TT=DELDER

		  DO I=1,100
			RKXP=Y(4)
			RKYP=Y(5)
			RKZP=Y(6)
			RRKXP=MAX(ABS(RKXP)*VV,TT)
			RRKYP=MAX(ABS(RKYP)*VV,TT)
			RRKZP=MAX(ABS(RKZP)*VV,TT)
c
			DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG)
     &			-DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
			DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG)
     &			-DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
			DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG)
     &			-DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)
c
			DS2 = DKXP**2+DKYP**2+DKZP**2
			FAC_NWTN = 10.D0
			DO IMODNWTN=1,10
			    FAC_NWTN = FAC_NWTN/10.D0
				Y(4) = RKXP - (DELTA*DKXP)/DS2*FAC_NWTN
				Y(5) = RKYP - (DELTA*DKYP)/DS2*FAC_NWTN
				Y(6) = RKZP - (DELTA*DKZP)/DS2*FAC_NWTN
				DDELTA = DISPXR(XP,YP,ZP,Y(4),Y(5),Y(6),OMG)
c				WRITE (6,*) 'DDELTA', DDELTA
			    IF (ABS(DDELTA) .LT. ABS(DELTA)) EXIT
			END DO
			IF (ABS(DDELTA).LT.1.0D-6) EXIT
			DELTA = DDELTA
		  END DO
		  
		DO I =1,3
			YK(I) = Y(I+3)
		END DO  

      RETURN
      END

