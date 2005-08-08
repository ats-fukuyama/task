C     $Id$
C     
***********************************************************************
C
      SUBROUTINE WRBEAM
C    
      INCLUDE 'wrcomm.inc'
C
      DIMENSION Y(NBEQ)
      REAL*4 TIME1,TIME2
      DATA MODEW/0/
C 
      CALL GUTIME(TIME1)
C
      IF(INTYPE.EQ.0) THEN
         WRITE(6,*) '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
         WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
         WRITE(6,'(1PE12.4,0P6F10.3)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
         WRITE(6,'(12X,1P5E12.4)') 
     &                      RCURVA,RCURVB,RBRADA,RBRADB,UUI
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
            WRITE(6,*) '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
            WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,',
     &                 'UU'
            WRITE(6,'(1PE12.4,0P6F10.3)') 
     &                         RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
            WRITE(6,'(12X,1P5E12.4)') 
     &                         RCURVA,RCURVB,RBRADA,RBRADB,UUI
         ELSEIF(INTYPE.EQ.2) THEN
            WRITE(6,*) '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH'
            WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,',
     &                 'UU'
            WRITE(6,'(1PE12.4,0P6F10.3)') 
     &                         RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
            WRITE(6,'(12X,1P5E12.4)') 
     &                         RCURVA,RCURVB,RBRADA,RBRADB,UUI
         ENDIF
      ENDIF
C
      DO NRAY=1,NRAYMX
C
 1       WRITE(6,*) '# NRAY=' ,NRAY
      IF(INTYPE.EQ.0) THEN
         READ(5,*,ERR=1,END=9000) 
     &                   RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,
     &                   RCURVA,RCURVB,RBRADA,RBRADB,UUI
         WRITE(6,*) '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
         WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
         WRITE(6,'(1PE12.4,0P6F10.3)') 
     &                      RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
         WRITE(6,'(12X,1P5E12.4)') 
     &                      RCURVA,RCURVB,RBRADA,RBRADB,UUI
      ELSEIF(INTYPE.EQ.1) THEN
         READ(5,*,ERR=1,END=9000) 
     &                   RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,
     &                   RCURVA,RCURVB,RBRADA,RBRADB,UUI
         WRITE(6,*) '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
         WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
         WRITE(6,'(1PE12.4,0P6F10.3)') 
     &                   RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
         WRITE(6,'(12X,1P5E12.4)') 
     &                      RCURVA,RCURVB,RBRADA,RBRADB,UUI
         SINP2=SIN(ANGZ*PI/180.D0)**2
         SINT2=SIN(ANGPH*PI/180.D0)**2
         RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
         RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
         IF(ANGZ.LT.0.D0) RNZI=-RNZI
         IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
      ELSEIF(INTYPE.EQ.2) THEN
         READ(5,*,ERR=1,END=9000) 
     &                     RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
         WRITE(6,*) '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH'
         WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
         WRITE(6,'(1PE12.4,0P6F10.3)') 
     &                     RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
         WRITE(6,'(12X,1P5E12.4)') 
     &                      RCURVA,RCURVB,RBRADA,RBRADB,UUI
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
C     
      IF(INTYPE.EQ.0)THEN
         RAYIN(6,NRAY)=RNZI
         RAYIN(7,NRAY)=RNPHII
      ELSE
         RAYIN(6,NRAY)=ANGZ
         RAYIN(7,NRAY)=ANGPH
      ENDIF
C
      RAYIN(8,NRAY)=UUI
C
      RKRI  = RKR0
      RKZI  =2.D6*PI*RF*RNZI  /VC
      RKPHII=2.D6*PI*RF*RNPHII/VC
      CALL WRNWTN(RKRI,RKZI,RKPHII,IERR)
      IF(IERR.NE.0) GOTO 1200
C
      CALL WRNWTN(RKRI,RKZI,RKPHII,IERR)
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
C
C      WRITE(6,'(1P4E12.4)') RCURVA,RCURVB,RBRADA,RBRADB
C
      CALL WRSETY(Y)
C
C      CALL WRGETY(Y,RCA1,RCB1,RBA1,RBB1,RTH1,RTH2,YY)
C      WRITE(6,'(1P5E12.4)') RCA1,RCB1,RBA1,RBB1,RTH1,RTH2
C      
      Y(19)= UUI
C
      CALL WRRKFTB(Y,RAYB,NITMAXB)
C
C      NRAYMX=1

C      NITMAX(1)=NITMAXB
      NITMAX(NRAY)=NITMAXB
      DO NIT=0,NITMAXB
         RAYS(0,NIT,NRAY)=RAYB( 0,NIT)
         RAYS(1,NIT,NRAY)=RAYB( 1,NIT)
         RAYS(2,NIT,NRAY)=RAYB( 2,NIT)
         RAYS(3,NIT,NRAY)=RAYB( 3,NIT)
         RAYS(4,NIT,NRAY)=RAYB( 4,NIT)
         RAYS(5,NIT,NRAY)=RAYB( 5,NIT)
         RAYS(6,NIT,NRAY)=RAYB( 6,NIT)
         RAYS(7,NIT,NRAY)=RAYB(19,NIT)
         RAYS(8,NIT,NRAY)=RAYB(20,NIT)
         RAYRB1(NIT,NRAY)=RAYB(23,NIT)
         RAYRB2(NIT,NRAY)=RAYB(24,NIT)
      ENDDO
C      CALL WRCALE(RAYS,NITMAX(1),1)
      CALL WRCALE(RAYS(0,0,NRAY),NITMAX(NRAY),NRAY)
      WRITE(6,'(A,F8.4)') 
     &        '# PABS/PIN=',UUI-RAYS(7,NITMAX(NRAY),NRAY)
 1200 CONTINUE
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
      SUBROUTINE WRSETY(Y)
C
      INCLUDE 'wrcomm.inc'      
      INCLUDE '../pl/plcom2.inc'      
C
      DIMENSION Y(NBEQ)
      DIMENSION RUM(3,3),VS(3,3),VP(3,3)
C
      RKABS=SQRT(Y(4)**2+Y(5)**2+Y(6)**2)
      CALL PLMAG(Y(1),Y(2),Y(3),RHON)
C
      RKNX=Y(4)/RKABS
      RKNY=Y(5)/RKABS
      RKNZ=Y(6)/RKABS
      RKBX=RKNY*BNZ-RKNZ*BNY
      RKBY=RKNZ*BNX-RKNX*BNZ
      RKBZ=RKNX*BNY-RKNY*BNX
      RKBABS=SQRT(RKBX**2+RKBY**2+RKBZ**2)
C
      RUM(1,1)=RKNX
      RUM(1,2)=RKNY
      RUM(1,3)=RKNZ
      RUM(2,1)=RKBX/RKBABS
      RUM(2,2)=RKBY/RKBABS
      RUM(2,3)=RKBZ/RKBABS
      RUM(3,1)=RUM(1,2)*RUM(2,3)-RUM(1,3)*RUM(2,2)
      RUM(3,2)=RUM(1,3)*RUM(2,1)-RUM(1,1)*RUM(2,3)
      RUM(3,3)=RUM(1,1)*RUM(2,2)-RUM(1,2)*RUM(2,1)
C
C      RL0=VC/DBLE(CW)
      RL0=1.D0/RKABS
      IF(RCURVA.EQ.0.D0)THEN
         VS22=0.D0
      ELSE
         VS22=1.D0/(RL0*RCURVA)
      ENDIF
      IF(RCURVB.EQ.0.D0)THEN
         VS33=0.D0
      ELSE
         VS33=1.D0/(RL0*RCURVB)
      ENDIF
      DO I=1,3
      DO J=1,3
         VS(I,J)=RUM(2,I)*VS22*RUM(2,J)
     &          +RUM(3,I)*VS33*RUM(3,J)
      ENDDO
      ENDDO
C
      IF(RBRADA.EQ.0.D0)THEN
         VP22=0.D0
      ELSE
         VP22=2.D0/(RBRADA*RBRADA)
      ENDIF
      IF(RBRADB.EQ.0.D0)THEN
         VP33=0.D0
      ELSE
         VP33=2.D0/(RBRADB*RBRADB)
      ENDIF
      DO I=1,3
      DO J=1,3
         VP(I,J)=RUM(2,I)*VP22*RUM(2,J)
     &          +RUM(3,I)*VP33*RUM(3,J)
      ENDDO
      ENDDO
C
      Y( 7)=VS(1,1)
      Y( 8)=VS(1,2)
      Y( 9)=VS(1,3)
      Y(10)=VS(2,2)
      Y(11)=VS(2,3)
      Y(12)=VS(3,3)
C
      Y(13)=VP(1,1)
      Y(14)=VP(1,2)
      Y(15)=VP(1,3)
      Y(16)=VP(2,2)
      Y(17)=VP(2,3)
      Y(18)=VP(3,3)
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRGETY(Y,RC1,RC2,RB1,RB2,RTH1,RTH2,YY)
C
      INCLUDE 'wrcomm.inc'      
      INCLUDE '../pl/plcom2.inc'      
C
      DIMENSION Y(NBEQ),YY(27),DFK(3)
      DIMENSION RUM(3,3),VS(3,3),VP(3,3),WS(3,3),WP(3,3)
C
      RKABS=SQRT(Y(4)**2+Y(5)**2+Y(6)**2)
      CALL PLMAG(Y(1),Y(2),Y(3),RHON)
C
      RKNX=Y(4)/RKABS
      RKNY=Y(5)/RKABS
      RKNZ=Y(6)/RKABS
      RKBX=RKNY*BNZ-RKNZ*BNY
      RKBY=RKNZ*BNX-RKNX*BNZ
      RKBZ=RKNX*BNY-RKNY*BNX
      RKBABS=SQRT(RKBX**2+RKBY**2+RKBZ**2)
C      WRITE(6,*) BNX,BNY,BNZ
C
      RUM(1,1)=RKNX
      RUM(1,2)=RKNY
      RUM(1,3)=RKNZ
      RUM(2,1)=RKBX/RKBABS
      RUM(2,2)=RKBY/RKBABS
      RUM(2,3)=RKBZ/RKBABS
      RUM(3,1)=RUM(1,2)*RUM(2,3)-RUM(1,3)*RUM(2,2)
      RUM(3,2)=RUM(1,3)*RUM(2,1)-RUM(1,1)*RUM(2,3)
      RUM(3,3)=RUM(1,1)*RUM(2,2)-RUM(1,2)*RUM(2,1)
C
      VS(1,1)=Y( 7)
      VS(1,2)=Y( 8)
      VS(1,3)=Y( 9)
      VS(2,1)=Y( 8)
      VS(2,2)=Y(10)
      VS(2,3)=Y(11)
      VS(3,1)=Y( 9)
      VS(3,2)=Y(11)
      VS(3,3)=Y(12)
C
      VP(1,1)=Y(13)
      VP(1,2)=Y(14)
      VP(1,3)=Y(15)
      VP(2,1)=Y(14)
      VP(2,2)=Y(16)
      VP(2,3)=Y(17)
      VP(3,1)=Y(15)
      VP(3,2)=Y(17)
      VP(3,3)=Y(18)
C
      DO I=1,3
      DO J=1,3
         WS(I,J)=0.D0
         WP(I,J)=0.D0
         DO K=1,3
         DO L=1,3
            WS(I,J)=WS(I,J)+RUM(I,K)*VS(K,L)*RUM(J,L)
            WP(I,J)=WP(I,J)+RUM(I,K)*VP(K,L)*RUM(J,L)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      YY(1)=WS(1,1)
      YY(2)=WS(1,2)
      YY(3)=WS(1,3)
      YY(4)=WS(2,2) 
      YY(5)=WS(2,3)
      YY(6)=WS(3,3)
C
      YY(7)=WP(1,1)
      YY(8)=WP(1,2)
      YY(9)=WP(1,3)
      YY(10)=WP(2,2) 
      YY(11)=WP(2,3)
      YY(12)=WP(3,3)
C
      YY(13)=RUM(1,1)
      YY(14)=RUM(1,2)
      YY(15)=RUM(1,3)
      YY(16)=RUM(2,1)
      YY(17)=RUM(2,2)
      YY(18)=RUM(2,3)
      YY(19)=RUM(3,1)
      YY(20)=RUM(3,2)
      YY(21)=RUM(3,3)
C
C     ***********   MODIFICATION BY NISHINA *********
C     calculating eigenvalue
C
C      SD=(WS(2,2)+WS(3,3))**2-4.D0*(WS(2,2)*WS(3,3)-WS(2,3)**2)
      SD=(WS(2,2)-WS(3,3))**2+4.D0*WS(2,3)**2
      RSLAM1=0.5D0*(WS(2,2)+WS(3,3)+SQRT(SD))
      RSLAM2=0.5D0*(WS(2,2)+WS(3,3)-SQRT(SD))
C      BD=(WP(2,2)+WP(3,3))**2-4.D0*(WP(2,2)*WP(3,3)-WP(2,3)**2)
      BD=(WP(2,2)-WP(3,3))**2+4.D0*WP(2,3)**2
      RBLAM1=0.5D0*(WP(2,2)+WP(3,3)+SQRT(BD))
      RBLAM2=0.5D0*(WP(2,2)+WP(3,3)-SQRT(BD))
C
C      WRITE(6,'(A,1P2E12.4)') 'SD,BD        =',SD,BD
C      WRITE(6,'(A,1P2E12.4)') 'RSLAM1,RSLAM2=',RSLAM1,RSLAM2
C      WRITE(6,'(A,1P2E12.4)') 'RBLAM1,RBLAM2=',RBLAM1,RBLAM2
C     
C     calculating curvature radius
C     
      RL0=1.D0/RKABS
      IF(RSLAM1.EQ.0.D0)THEN
         RC1=0.D0
      ELSE
         RC1=1.D0/(RL0*RSLAM1)
      ENDIF
      IF(RSLAM2.EQ.0.D0)THEN
         RC2=0.D0
      ELSE
         RC2=1.D0/(RL0*RSLAM2)
      ENDIF
C
C     calculating beam radius
C
      IF(RBLAM1.EQ.0.D0)THEN
         RB1=0.D0
      ELSE
         RB1=SQRT(2.D0/RBLAM1)
      ENDIF
      IF(RBLAM2.EQ.0.D0)THEN
         RB2=0.D0
      ELSE
         RB2=SQRT(2.D0/RBLAM2)
      ENDIF
C
C      WRITE(6,'(A,1P2E12.4)') 'RC1,RC2      =',RC1,RC2
C      WRITE(6,'(A,1P2E12.4)') 'RB1,RB2      =',RB1,RB2
C 
C     calculating angle
C
      RTH1=ATAN2(-WP(2,3),WP(2,2)-RBLAM1)*180.D0/PI
      RTH2=ATAN2(-WP(2,3),WP(2,2)-RBLAM2)*180.D0/PI
      IF(RTH1.LT.0.D0) RTH1=RTH1+180.D0
      IF(RTH2.LT.0.D0) RTH2=RTH2+180.D0
C
C      WRITE(6,'(A,1P4E12.4)') 'RB1,RB2,RTH1,RTH2=',RB1,RB2,RTH1,RTH2
C
C
C     -----------------CALC VG---------------------------------------
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
C
      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      DKXP=MAX(ABS(RKXP)*VV,TT)
      DKYP=MAX(ABS(RKYP)*VV,TT)
      DKZP=MAX(ABS(RKZP)*VV,TT)   
C
      DOMG=(DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG)
     &     -DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
C
      F000P00=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F000M00=DISPBXR(XP,YP,ZP,RKXP-DKXP,RKYP,RKZP,OMG)
      F0000P0=DISPBXR(XP,YP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F0000M0=DISPBXR(XP,YP,ZP,RKXP,RKYP-DKYP,RKZP,OMG)
      F00000P=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F00000M=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP-DKZP,OMG)
C
      DFK(1)=(F000P00-F000M00)/(2.D0*DKXP)
      DFK(2)=(F0000P0-F0000M0)/(2.D0*DKYP)
      DFK(3)=(F00000P-F00000M)/(2.D0*DKZP)
C
      VGX=DFK(1)/DOMG
      VGY=DFK(2)/DOMG
      VGZ=DFK(3)/DOMG
      RVGABS=SQRT(VGX**2+VGY**2+VGZ**2)
C
      YY(22)=RKXP/RKABS
      YY(23)=RKYP/RKABS
      YY(24)=RKZP/RKABS
      YY(25)=VGX/RVGABS
      YY(26)=VGY/RVGABS
      YY(27)=VGZ/RVGABS
C
      RETURN
      END
C
C***********************************************************************
C************************************************************************
C
      SUBROUTINE WRRKFTB(Y,YN,NIT)
C
      INCLUDE 'wrcomm.inc'      
C
      EXTERNAL WRFDRVB
      DIMENSION Y(NBEQ),YM(NBEQ),YN(0:NBVAR,0:NITM),WORK(NBEQ,2)
      DIMENSION YY(27),F(NBEQ)
C
      ITMAX=INT(SMAX/DELS)
C      WRITE(6,'(A,1P2E12.4,I5)') 'SMAX,DELS,ITMAX=',SMAX,DELS,ITMAX
C
      IT=0
      X0 = 0.D0
      XE = DELS     
      YN(0,IT)=X0
      DO I=1,19
         YN(I,IT)=Y(I)
      ENDDO
      YN(20,IT)=0.D0
      CALL WRGETY(Y,YN(21,IT),YN(22,IT),YN(23,IT),YN(24,IT),
     &                                  YN(25,IT),YN(26,IT),YY)
      DO I=27,53
         YN(I,IT)=0.D0
      ENDDO
C
C         WRITE(6,6002) Y( 7),Y( 8),Y( 9),Y(10),Y(11),Y(12)
C         WRITE(6,6002) Y(13),Y(14),Y(15),Y(16),Y(17),Y(18)
C
      CALL WRFDRVB (X,Y,F)
C
C         CALL WRGETY(Y,RCA,RCB,RBA,RBB,RTH,YY)
C         WRITE(6,6003) RCA,RCB,RBA,RBB,RTH
C
      DO IT = 1,ITMAX
         YPRE=Y(19)
         CALL RK(NBEQ,WRFDRVB,X0,XE,1,Y,YM,WORK)
         YN(0,IT)=XE
         DO I=1,19
           YN(I,IT)=YM(I)
         ENDDO
         YN(20,IT)=YPRE-YM(19)
C
         CALL WRGETY(YM,YN(21,IT),YN(22,IT),YN(23,IT),YN(24,IT),
     &                  YN(25,IT),YN(26,IT),YY)
C
         DO I=27,53
           YN(I,IT)=YY(I-26)
         ENDDO
C
         CALL WRFDRVB (X,Y,F)
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
         WRITE(6,6001) XE,RL,PHIL,ZL,RKRL,YM(19),YN(20,IT)
 6001    FORMAT(1H ,1P7E11.3)
C         WRITE(6,6002) Y( 7),Y( 8),Y( 9),Y(10),Y(11),Y(12)
C         WRITE(6,6002) Y(13),Y(14),Y(15),Y(16),Y(17),Y(18)
C 6002    FORMAT(1H ,11X,1P6E11.3)
         WRITE(6,6003) YN(21,IT),YN(22,IT),YN(23,IT),YN(24,IT),
     &                 YN(25,IT),YN(26,IT)
 6003    FORMAT(1H ,11X,1P6E11.3)
C
         DO I=1,19
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(19).LT.UUMIN) THEN
            NIT = IT
            WRITE(6,*) '--- Absorbed ---'
            GOTO 10
         ENDIF         
         CALL PLMAG(Y(1),Y(2),Y(3),RHON)
C         IF(RHON.GT.RB/RA) THEN
          IF(PSIN.GT.2.D0) THEN
            NIT = IT
            WRITE(6,*) '--- Out of bounds ---'
            GOTO 10
         ENDIF
      ENDDO
      NIT=ITMAX
C     
 10   IF(YN(19,NIT).LT.0.D0) THEN
         YN(19,NIT)=0.D0
      ENDIF
C
      RETURN
      END
C
C************************************************************************
C
      SUBROUTINE WRFDRVB(X,Y,F) 
C
      INCLUDE 'wrcomm.inc'      
C
      DIMENSION Y(NBEQ),F(NBEQ)
      DIMENSION DFR(3),DFK(3),DFRR(3,3),DFKK(3,3),DFRK(3,3)
      DIMENSION S(3,3),P(3,3)
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
      UU=Y(19)
C
      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      DRXP=MAX(ABS(XP)*VV,TT)
      DRYP=MAX(ABS(YP)*VV,TT)
      DRZP=MAX(ABS(ZP)*VV,TT)
      DKXP=MAX(ABS(RKXP)*VV,TT)
      DKYP=MAX(ABS(RKYP)*VV,TT)
      DKZP=MAX(ABS(RKZP)*VV,TT)   
C
      DOMG=(DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG)
     &     -DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      F000000=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)
C
      FP00000=DISPBXR(XP+DRXP,YP,ZP,RKXP,RKYP,RKZP,OMG)
      FM00000=DISPBXR(XP-DRXP,YP,ZP,RKXP,RKYP,RKZP,OMG)
      F0P0000=DISPBXR(XP,YP+DRYP,ZP,RKXP,RKYP,RKZP,OMG)
      F0M0000=DISPBXR(XP,YP-DRYP,ZP,RKXP,RKYP,RKZP,OMG)
      F00P000=DISPBXR(XP,YP,ZP+DRZP,RKXP,RKYP,RKZP,OMG)
      F00M000=DISPBXR(XP,YP,ZP-DRZP,RKXP,RKYP,RKZP,OMG)
      F000P00=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F000M00=DISPBXR(XP,YP,ZP,RKXP-DKXP,RKYP,RKZP,OMG)
      F0000P0=DISPBXR(XP,YP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F0000M0=DISPBXR(XP,YP,ZP,RKXP,RKYP-DKYP,RKZP,OMG)
      F00000P=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F00000M=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP-DKZP,OMG)
C
      FPP0000=DISPBXR(XP+DRXP,YP+DRYP,ZP,RKXP,RKYP,RKZP,OMG)
      FP0P000=DISPBXR(XP+DRXP,YP,ZP+DRZP,RKXP,RKYP,RKZP,OMG)
      FP00P00=DISPBXR(XP+DRXP,YP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      FP000P0=DISPBXR(XP+DRXP,YP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      FP0000P=DISPBXR(XP+DRXP,YP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F0PP000=DISPBXR(XP,YP+DRYP,ZP+DRZP,RKXP,RKYP,RKZP,OMG)
      F0P0P00=DISPBXR(XP,YP+DRYP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F0P00P0=DISPBXR(XP,YP+DRYP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F0P000P=DISPBXR(XP,YP+DRYP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F00PP00=DISPBXR(XP,YP,ZP+DRZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F00P0P0=DISPBXR(XP,YP,ZP+DRZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F00P00P=DISPBXR(XP,YP,ZP+DRZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F000PP0=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP+DKYP,RKZP,OMG)
      F000P0P=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP,RKZP+DKZP,OMG)
      F0000PP=DISPBXR(XP,YP,ZP,RKXP,RKYP+DKYP,RKZP+DKZP,OMG)
C
      DFR(1)=(FP00000-FM00000)/(2.D0*DRXP)
      DFR(2)=(F0P0000-F0M0000)/(2.D0*DRYP)
      DFR(3)=(F00P000-F00M000)/(2.D0*DRZP)
      DFK(1)=(F000P00-F000M00)/(2.D0*DKXP)
      DFK(2)=(F0000P0-F0000M0)/(2.D0*DKYP)
      DFK(3)=(F00000P-F00000M)/(2.D0*DKZP)
C
      DFRR(1,1)=(FP00000-F000000-F000000+FM00000)/(DRXP*DRXP)
      DFRR(1,2)=(FPP0000-FP00000-F0P0000+F000000)/(DRXP*DRYP)
      DFRR(1,3)=(FP0P000-FP00000-F00P000+F000000)/(DRXP*DRZP)
      DFRR(2,1)=DFRR(1,2)
      DFRR(2,2)=(F0P0000-F000000-F000000+F0M0000)/(DRYP*DRYP)
      DFRR(2,3)=(F0PP000-F0P0000-F00P000+F000000)/(DRYP*DRZP)
      DFRR(3,1)=DFRR(1,3)
      DFRR(3,2)=DFRR(2,3)
      DFRR(3,3)=(F00P000-F000000-F000000+F00M000)/(DRZP*DRZP)
C
C      WRITE(6,'(A,1P3E12.4)') 'DRP  :',DRXP,DRYP,DRZP
C      WRITE(6,'(A,1P3E12.4)') 'DKP  :',DKXP,DKYP,DKZP
C      WRITE(6,'(A,1P3E12.4)') 'DFR  :',DFR(1),DFR(2),DFR(3)
C      WRITE(6,'(A,1P3E12.4)') 'DFK  :',DFK(1),DFK(2),DFK(3)
C      WRITE(6,'(A,1P3E12.4)') 'DFRR :',DFRR(1,1),DFRR(1,2),DFRR(1,3)
C      WRITE(6,'(A,1P3E12.4)') 'DFRR :',DFRR(2,1),DFRR(2,2),DFRR(2,3)
C      WRITE(6,'(A,1P6E12.4)') 'DFRR :',DFRR(3,1),DFRR(3,2),DFRR(3,3)
C
      DFKK(1,1)=(F000P00-F000000-F000000+F000M00)/(DKXP*DKXP)
      DFKK(1,2)=(F000PP0-F000P00-F0000P0+F000000)/(DKXP*DKYP)
      DFKK(1,3)=(F000P0P-F000P00-F00000P+F000000)/(DKXP*DKZP)
      DFKK(2,1)=DFKK(1,2)
      DFKK(2,2)=(F0000P0-F000000-F000000+F0000M0)/(DKYP*DKYP)
      DFKK(2,3)=(F0000PP-F0000P0-F00000P+F000000)/(DKYP*DKZP)
      DFKK(3,1)=DFKK(1,3)
      DFKK(3,2)=DFKK(2,3)
      DFKK(3,3)=(F00000P-F000000-F000000+F00000M)/(DKZP*DKZP)
C
C      WRITE(6,'(A,1P3E12.4)') 'DFKK :',DFKK(1,1),DFKK(1,2),DFKK(1,3)
C      WRITE(6,'(A,1P3E12.4)') 'DFKK :',DFKK(2,1),DFKK(2,2),DFKK(2,3)
C      WRITE(6,'(A,1P3E12.4)') 'DFKK :',DFKK(3,1),DFKK(3,2),DFKK(3,3)
C
      DFRK(1,1)=(FP00P00-FP00000-F000P00+F000000)/(DRXP*DKXP)
      DFRK(1,2)=(FP000P0-FP00000-F0000P0+F000000)/(DRXP*DKYP)
      DFRK(1,3)=(FP0000P-FP00000-F00000P+F000000)/(DRXP*DKZP)
      DFRK(2,1)=(F0P0P00-F0P0000-F000P00+F000000)/(DRYP*DKXP)
      DFRK(2,2)=(F0P00P0-F0P0000-F0000P0+F000000)/(DRYP*DKYP)
      DFRK(2,3)=(F0P000P-F0P0000-F00000P+F000000)/(DRYP*DKZP)
      DFRK(3,1)=(F00PP00-F00P000-F000P00+F000000)/(DRZP*DKXP)
      DFRK(3,2)=(F00P0P0-F00P000-F0000P0+F000000)/(DRZP*DKYP)
      DFRK(3,3)=(F00P00P-F00P000-F00000P+F000000)/(DRZP*DKZP)
C
C      WRITE(6,'(A,1P3E12.4)') 'DFRK :',DFRK(1,1),DFRK(1,2),DFRK(1,3)
C      WRITE(6,'(A,1P3E12.4)') 'DFRK :',DFRK(2,1),DFRK(2,2),DFRK(2,3)
C      WRITE(6,'(A,1P3E12.4)') 'DFRK :',DFRK(3,1),DFRK(3,2),DFRK(3,3)
C
      IF(DOMG.GT.0.D0) THEN
         DS=-1.D0/SQRT(DFK(1)**2+DFK(2)**2+DFK(3)**2)
      ELSE
         DS= 1.D0/SQRT(DFK(1)**2+DFK(2)**2+DFK(3)**2)
      ENDIF
C
      S(1,1)=Y( 7)
      S(1,2)=Y( 8)
      S(1,3)=Y( 9)
      S(2,1)=Y( 8)
      S(2,2)=Y(10)
      S(2,3)=Y(11)
      S(3,1)=Y( 9)
      S(3,2)=Y(11)
      S(3,3)=Y(12)
C
      P(1,1)=Y(13)
      P(1,2)=Y(14)
      P(1,3)=Y(15)
      P(2,1)=Y(14)
      P(2,2)=Y(16)
      P(2,3)=Y(17)
      P(3,1)=Y(15)
      P(3,2)=Y(17)
      P(3,3)=Y(18)
C
      F( 1)= DFK(1)*DS
      F( 2)= DFK(2)*DS
      F( 3)= DFK(3)*DS
      F( 4)=-DFR(1)*DS
      F( 5)=-DFR(2)*DS
      F( 6)=-DFR(3)*DS
C
      F( 7)=-DFRR(1,1)*DS
      F( 8)=-DFRR(1,2)*DS
      F( 9)=-DFRR(1,3)*DS
      F(10)=-DFRR(2,2)*DS
      F(11)=-DFRR(2,3)*DS
      F(12)=-DFRR(3,3)*DS
C
      DO I=1,3
         F( 7)=F( 7)-(DFRK(1,I)*S(1,I)+DFRK(1,I)*S(1,I))*DS
         F( 8)=F( 8)-(DFRK(2,I)*S(1,I)+DFRK(1,I)*S(2,I))*DS
         F( 9)=F( 9)-(DFRK(3,I)*S(1,I)+DFRK(1,I)*S(3,I))*DS
         F(10)=F(10)-(DFRK(2,I)*S(2,I)+DFRK(2,I)*S(2,I))*DS
         F(11)=F(11)-(DFRK(3,I)*S(2,I)+DFRK(2,I)*S(3,I))*DS
         F(12)=F(12)-(DFRK(3,I)*S(3,I)+DFRK(3,I)*S(3,I))*DS
      ENDDO
C
      DO I=1,3
      DO J=1,3
         F( 7)=F( 7)
     &        -(DFKK(I,J)*S(1,I)*S(1,J)-DFKK(I,J)*P(1,I)*P(1,J))*DS
         F( 8)=F( 8)
     &        -(DFKK(I,J)*S(1,I)*S(2,J)-DFKK(I,J)*P(1,I)*P(2,J))*DS
         F( 9)=F( 9)
     &        -(DFKK(I,J)*S(1,I)*S(3,J)-DFKK(I,J)*P(1,I)*P(3,J))*DS
         F(10)=F(10)
     &        -(DFKK(I,J)*S(2,I)*S(2,J)-DFKK(I,J)*P(2,I)*P(2,J))*DS
         F(11)=F(11)
     &        -(DFKK(I,J)*S(2,I)*S(3,J)-DFKK(I,J)*P(2,I)*P(3,J))*DS
         F(12)=F(12)
     &        -(DFKK(I,J)*S(3,I)*S(3,J)-DFKK(I,J)*P(3,I)*P(3,J))*DS
      ENDDO
      ENDDO
C
      F(13)=0.D0
      F(14)=0.D0
      F(15)=0.D0
      F(16)=0.D0
      F(17)=0.D0
      F(18)=0.D0
C
      DO I=1,3
         F(13)=F(13)-(DFRK(1,I)*P(1,I)+DFRK(1,I)*P(1,I))*DS
         F(14)=F(14)-(DFRK(2,I)*P(1,I)+DFRK(1,I)*P(2,I))*DS
         F(15)=F(15)-(DFRK(3,I)*P(1,I)+DFRK(1,I)*P(3,I))*DS
         F(16)=F(16)-(DFRK(2,I)*P(2,I)+DFRK(2,I)*P(2,I))*DS
         F(17)=F(17)-(DFRK(3,I)*P(2,I)+DFRK(2,I)*P(3,I))*DS
         F(18)=F(18)-(DFRK(3,I)*P(3,I)+DFRK(3,I)*P(3,I))*DS
      ENDDO
C
      DO I=1,3
      DO J=1,3
         F(13)=F(13)
     &        -(DFKK(I,J)*S(1,I)*P(1,J)+DFKK(J,I)*S(1,I)*P(1,J))*DS
         F(14)=F(14)
     &        -(DFKK(I,J)*S(1,I)*P(2,J)+DFKK(J,I)*S(2,I)*P(1,J))*DS
         F(15)=F(15)
     &        -(DFKK(I,J)*S(1,I)*P(3,J)+DFKK(J,I)*S(3,I)*P(1,J))*DS
         F(16)=F(16)
     &        -(DFKK(I,J)*S(2,I)*P(2,J)+DFKK(J,I)*S(2,I)*P(2,J))*DS
         F(17)=F(17)
     &        -(DFKK(I,J)*S(2,I)*P(3,J)+DFKK(J,I)*S(3,I)*P(2,J))*DS
         F(18)=F(18)
     &        -(DFKK(I,J)*S(3,I)*P(3,J)+DFKK(J,I)*S(3,I)*P(3,J))*DS
      ENDDO
      ENDDO
C
C      VDU  =-2*DISPBXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)*DS
C
      VDU  =-2.D0*ABS(DISPBXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)*DS)
C
      IF(UU.LT.0.D0) THEN
         F(19)=0.D0
      ELSE
         F(19)=VDU*UU 
      ENDIF
C
C      WRITE(6,'(A,1P6E12.4)') 'X:',X
C      WRITE(6,'(A,1P6E12.4)') 'Y :',Y( 1),Y( 2),Y( 3),Y( 4),Y( 5),Y( 6)
C      WRITE(6,'(A,1P6E12.4)') 'Y :',Y( 7),Y( 8),Y( 9),Y(10),Y(11),Y(12)
C      WRITE(6,'(A,1P6E12.4)') 'Y :',Y(13),Y(14),Y(15),Y(16),Y(17),Y(18)
C      WRITE(6,'(A,1P6E12.4)') 'F :',F( 1),F( 2),F( 3),F( 4),F( 5),F( 6)
C      WRITE(6,'(A,1P6E12.4)') 'F :',F( 7),F( 8),F( 9),F(10),F(11),F(12)
C      WRITE(6,'(A,1P6E12.4)') 'F :',F(13),F(14),F(15),F(16),F(17),F(18)
C      CALL GUFLSH
C
      RETURN
      END
C
C************************************************************************
C
      FUNCTION DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)
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
      DISPBXR=DBLE(CF)
      RETURN
      END
C
C***********************************************************************
C
      FUNCTION DISPBXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)
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
      DISPBXI=DIMAG(CF)
      RETURN
      END
