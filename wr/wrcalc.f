C     $Id$
C     
C***********************************************************************
C
      SUBROUTINE WRCALC
C    
      INCLUDE 'wrcomm.inc'
C
      DIMENSION Y(NEQ)
      REAL*4 TIME1,TIME2
      DATA MODEW/0/
C 
      CALL GUTIME(TIME1)
C
      IF(INTYPE.EQ.0) THEN
         WRITE(6,*) 
     &   '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
         WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
      ELSEIF(INTYPE.EQ.1.OR.INTYPE.EQ.2) THEN
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
         IF(INTYPE.EQ.1) THEN
            WRITE(6,*) 
     &      '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                         RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
         ELSEIF(INTYPE.EQ.2) THEN
            WRITE(6,*) 
     &      '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                         RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
         ENDIF
      ENDIF
C
      DO NRAY=1,NRAYMX
C
    1    WRITE(6,*) '# NRAY = ',NRAY
         IF(INTYPE.EQ.0) THEN
            READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
            WRITE(6,*) 
     &      '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
         ELSEIF(INTYPE.EQ.1) THEN
            READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
            WRITE(6,*) 
     &      '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
            SINP2=SIN(ANGZ*PI/180.D0)**2
            SINT2=SIN(ANGPH*PI/180.D0)**2
            RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
            RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
            IF(ANGZ.LT.0.D0) RNZI=-RNZI
            IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
         ELSEIF(INTYPE.EQ.2) THEN
            READ(5,*,ERR=1,END=9000) 
     &                      RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
            WRITE(6,*) 
     &      '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
            WRITE(6,'(1PE12.4,0P7F9.2)') 
     &                      RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
            SINP2=SIN(ANGZ*PI/180.D0)**2
            SINT2=SIN(ANGPH*PI/180.D0)**2
            RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
            RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
            WRITE(6,*) 'XX INTYPE=2 IS NOT SUPPORTED YET.'
            GOTO 1
         ENDIF
C
         RAYIN(1,NRAY)=RF
         RAYIN(2,NRAY)=RPI
         RAYIN(3,NRAY)=ZPI
         RAYIN(4,NRAY)=PHII
         RAYIN(5,NRAY)=RKR0
         RAYIN(6,NRAY)=RNZI
         RAYIN(7,NRAY)=RNPHII
         RAYIN(8,NRAY)=UUI
C          
         RKRI  = RKR0
         RKZI  =2.D6*PI*RF*RNZI  /VC
         RKPHII=2.D6*PI*RF*RNPHII/VC
         CALL WRNWTN(RKRI,RKZI,RKPHII,IERR)
         IF(IERR.NE.0) GOTO 1200
C         
         WRITE(6,*) 'RKRI=',RKRI
C
         IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
            Y(1)= RPI
            Y(2)= 2.D0*PI*RR*SIN(PHII)
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
         IF(IQTYPE.EQ.0) THEN
            CALL WRRKFT(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSEIF(IQTYPE.EQ.1) THEN
            CALL WRRKFT_ODE(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSEIF(IQTYPE.EQ.2) THEN
            CALL WRRKFT_RKF(Y,RAYS(0,0,NRAY),NITMAX(NRAY))
         ELSE
            WRITE(6,*) 'XX WRCALC: unknown IQTYPE =', IQTYPE
         ENDIF
         CALL WRCALE(RAYS(0,0,NRAY),NITMAX(NRAY),NRAY)
         WRITE(6,'(A,F8.4)') 
     &        '# PABS/PIN=',1.D0-RAYS(7,NITMAX(NRAY),NRAY)
 1200    CONTINUE
      ENDDO
C
      CALL GUTIME(TIME2)
      WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'
C
 9000 RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT_ODE(Y,YN,NIT)
C
      INCLUDE 'wrcomm.inc'      
C
      PARAMETER (NMAX=50,KMAXX=200)
      COMMON /DPODE1/ kmax,kount,dxsav,xp(KMAXX),yp(NMAX,KMAXX)
      EXTERNAL WRFDRV,BSSTEP,RKQS
      DIMENSION Y(NEQ),YN(0:NEQ,0:NITM)
C
      X0 = 0.D0
      XE = DELS     
      ITMAX=INT(SMAX/DELS)
      H1=DELS
      HMIN=DELRAY
      KMAX=0
C
      IT=0
      YN(0,IT)=X0
      DO I=1,7
         YN(I,IT)=Y(I)
      ENDDO
      YN(8,IT)=0.D0
C
      DO 10 IT = 1,ITMAX
         Y7=Y(7)
         CALL ODEINT(Y,7,X0,XE,EPSRAY,H1,HMIN,NOK,NBAD,WRFDRV,
     &               BSSTEP)
C     &               RKQS)
         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=Y(I)
         ENDDO
         Y(8)=Y7-Y(7)
         YN(8,IT)=Y(8)
         WRITE(6,6001) XE,SQRT(Y(1)**2+Y(2)**2),ATAN2(Y(2),Y(1)),Y(3),
     &                 ( Y(4)*Y(1)+Y(5)*Y(2))/SQRT(Y(1)**2+Y(2)**2),
     &                 -Y(4)*Y(2)+Y(5)*Y(1),
     &                 Y(7),NBAD
 6001    FORMAT(1P7E11.3,I2)
C         WRITE(25,'(1P8E9.1)') XE,(Y(I),I=1,7)
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF         
         CALL PLMAG(Y(1),Y(2),Y(3),PSIN)
         IF(PSIN.GT.(RB/RA)**2) THEN
            NIT = IT
            GOTO 11
         ENDIF         
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT_RKF(Y,YN,NIT)
C
      INCLUDE 'wrcomm.inc'      
C
      EXTERNAL WRFDRV
      DIMENSION Y(NEQ),YM(NEQ),YN(0:NEQ,0:NITM),ESTERR(NEQ),WORK1(NEQ),
     &           WORK2(NEQ),WORK3(NEQ),WORK4(NEQ,11)
C
      RELERR = EPSRAY
      ABSERR = EPSRAY
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
C
      DO 10 IT = 1,ITMAX
         Y7=Y(7)
         CALL RKF(7,WRFDRV,X0,XE,Y,INIT,RELERR,ABSERR,YM,
     &            ESTERR,NDE,IER,WORK0,WORK1,WORK2,WORK3,WORK4)
         IF (IER .NE. 0) THEN
            WRITE(6,*) 'XX RKF ERROR'
            RETURN
         ENDIF
C         WRITE(6,'(A,1P7E9.1)') 'ESTERR=',ESTERR
         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=YM(I)
         ENDDO
         YN(8,IT)=Y7-YM(7)
         WRITE(6,6001) X0,YN(1,IT),YN(2,IT),YN(3,IT),YN(4,IT),
     &                    YN(7,IT),YN(8,IT)
 6001    FORMAT(1H ,1P7E11.3)
C     &           '    Z=',YN(2,IT),'    PHI=',YN(3,IT),'    R=',YN(1,IT),
C     &           '    KR=',YN(4,IT),'    KZ=',YN(5,IT),
C     &           '    KPHI=',YN(6,IT),' POWER=',YN(8,IT), '  Y(7)=',Y7 
         DO I=1,7
            Y(I)=YN(I,IT)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NIT = IT
            GOTO 11
         ENDIF         
         CALL PLMAG(Y(1),Y(2),Y(3),PSIN)
         IF(PSIN.GT.(RB/RA)**2) THEN
            NIT = IT
            GOTO 11
         ENDIF         
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRRKFT(Y,YN,NIT)
C
      INCLUDE 'wrcomm.inc'      
C
      EXTERNAL WRFDRV
      DIMENSION Y(NEQ),YM(NEQ),YN(0:NEQ,0:NITM),WORK(2,NEQ)
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
         CALL RK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
C
         YN(0,IT)=XE
         DO I=1,7
            YN(I,IT)=YM(I)
         ENDDO
         YN(8,IT)=Y7-YM(7)
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
         WRITE(6,6001) XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,IT)
 6001    FORMAT(1H ,1P7E11.3)
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
         CALL PLMAG(Y(1),Y(2),Y(3),PSIN)
         IF(PSIN.GT.(RB/RA)**2) THEN
            NIT = IT
            GOTO 11
         ENDIF         
 10   CONTINUE
      NIT=ITMAX
C     
 11   IF(YN(7,NIT).LT.0.D0) THEN
         YN(7,NIT)=0.D0
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRFDRV(X,Y,F) 
C
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
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRNWTN(RKRI,RKZI,RKPHII,IERR)
C
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
C      WRITE(6,*) X1
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
      INCLUDE 'wrcomm.inc'
      INCLUDE '../pl/plcom2.inc'
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
         IF(MODELP(NS).GE.10.AND.MODELP(NS).LT.30) MODELP(NS)=0
         IF(MODELP(NS).GE.30.AND.MODELP(NS).LT.50) MODELP(NS)=4
         IF(MODELP(NS).GE.50.AND.MODELP(NS).LT.60) MODELP(NS)=5
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
      INCLUDE 'wrcomm.inc'
      INCLUDE '../pl/plcom2.inc'
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
         IF(MODELP(NS).GE.10.AND.MODELP(NS).LT.30) MODELP(NS)=0
         IF(MODELP(NS).GE.30.AND.MODELP(NS).LT.50) MODELP(NS)=4
         IF(MODELP(NS).GE.50.AND.MODELP(NS).LT.60) MODELP(NS)=5
      ENDDO
C
      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)
C      CF=CFDISPR(CRF,CKX,CKY,CKZ,X,Y,Z)
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
      DISPXR=DBLE(CF)
      RETURN
      END
C
C***********************************************************************
C
      FUNCTION DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)
C
      INCLUDE 'wrcomm.inc'
      INCLUDE '../pl/plcom2.inc'
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
         CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
         CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
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
      ENDDO
C
      RETURN
      END
C
C***********************************************************************
C     save ray data
C***********************************************************************
C
      SUBROUTINE WRSAVE
C
      INCLUDE 'wrcomm.inc'
C      
      CHARACTER*72 KNAM
C      CHARACTER*1 KID
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : WRDATA SAVE FILE NAME : ',KNAMWR
      READ(5,'(A72)',ERR=1,END=900) KNAM
      IF(KNAM(1:2).NE.'/ ') KNAMWR=KNAM
C
      INQUIRE(FILE=KNAMWR,EXIST=LEX)
      IF(LEX) THEN
C         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
C     &              'ARE YOU SURE {Y/N} ?'
C         READ(5,502) KID
C  502    FORMAT(A1)
C         CALL GUCPTL(KID)
C         IF(KID.NE.'Y') GOTO 1
         OPEN(21,FILE=KNAMWR,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMWR,') IS ASSIGNED FOR OUTPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ELSE
         OPEN(21,FILE=KNAMWR,IOSTAT=IST,STATUS='NEW',ERR=20,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (',KNAMWR,') IS CREATED FOR OUTPUT.'
         GOTO 30
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ENDIF
C
   30 CONTINUE
      WRITE(21) RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ
      WRITE(21) NSMAX
      WRITE(21) (PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS),NS=1,NSMAX)
      WRITE(21) (PTPR(NS),PTPP(NS),PTS(NS),PU(NS),PUS(NS),NS=1,NSMAX)
      WRITE(21) PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
      WRITE(21) RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB
      WRITE(21) (MODELP(NS),NDISP1(NS),NDISP2(NS),NS=1,NSMAX)
      WRITE(21) MODELG,MODELN,MODELQ,MODELV
      WRITE(21) KNAMEQ,KNAMFP
C
      WRITE(21) NRAYMX
      DO NRAY=1,NRAYMX
         WRITE(21) NITMAX(NRAY)
         WRITE(21) (RAYIN(I,NRAY),I=1,8)
         WRITE(21) (CEXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (CEYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (CEZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RKXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RKYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RKZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RAYRB1(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(21) (RAYRB2(NIT,NRAY),NIT=0,NITMAX(NRAY))
         DO I=0,8
            WRITE(21) (RAYS(I,NIT,NRAY),NIT=0,NITMAX(NRAY))
         ENDDO
      ENDDO
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
  900 RETURN
      END
C
C********************************** CAL K ********************************
C
      SUBROUTINE WRCALK(NIT,NRAY,RKPARA,RKPERP)
C
      INCLUDE 'wrcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      X=RAYS(1,NIT,NRAY)
      Y=RAYS(2,NIT,NRAY)
      Z=RAYS(3,NIT,NRAY)
C
      CALL PLMAG(X,Y,Z,PSIN)
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
