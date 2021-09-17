C     $Id$
C
C     ****** PSI ******
C
      SUBROUTINE WFSPSI(X,Y,PSIN)
C
      INCLUDE 'wfcomm.inc'
C
      SELECT CASE(MODELB)
      CASE(0)
         PSIN=Y*Y/(RA*RA)
      CASE(1)
         PSIN=X*X/(RA*RA)
      CASE(2)
         PSIN=(X*X+Y*Y)/(RA*RA)
      CASE(3,5,6,7)
         CALL WFBPSI(X,Y,PSI)
         CALL WFBPSI(RA,0.D0,PSIA)
         PSIN=PSI/PSIA
      CASE(4)
         CALL WFBPSI(X,Y,PSI)
         CALL WFBPSI(RA,0.D0,PSIA)
         PSIN=PSI**2/PSIA**2
      END SELECT
      RETURN
      END
C
C     ****** MAGNETIC FIELD PROFILE ******
C
      SUBROUTINE WFBMAG(X,Y,BABS,AL)
C
      USE libbes,ONLY: BESINX
      INCLUDE 'wfcomm.inc'
C
      DIMENSION BLO(3),AL(3)
C
      IF(MODELB.EQ.0) THEN
         BLO(1)=BB
         BLO(2)=0.D0
         BLO(3)=0.D0
      ELSEIF(MODELB.EQ.1) THEN
         BLO(1)=0.D0
         BLO(2)=BB
         BLO(3)=0.D0
      ELSEIF(MODELB.EQ.2) THEN
         BLO(1)=0.D0
         BLO(2)=0.D0
         BLO(3)=BB
      ELSEIF(MODELB.EQ.3) THEN
         A0=0.5D0*(1.D0+RMIR)*BB
         A1=0.5D0*(1.D0-RMIR)*BB
         RL=0.5D0*PI*X/ZBB
         BLO(1)=  -A1*SIN(PI*Y/ZBB)*BESINX(1,RL)
         BLO(2)=A0+A1*COS(PI*Y/ZBB)*BESINX(0,RL)
         BLO(3)=0.D0
      ELSEIF(MODELB.EQ.4) THEN
         A0=0.5D0*(1.D0+RMIR)*BB
         A1=0.5D0*(1.D0-RMIR)*BB
         BLO(1)=   A1*SIN(PI*Y/ZBB)*SINH(PI*X/ZBB)
         BLO(2)=A0+A1*COS(PI*Y/ZBB)*COSH(PI*X/ZBB)
         BLO(3)=0.D0
      ELSEIF(MODELB.EQ.5) THEN
         SR =X*X+Y*Y/RKAP**2
         SRA=RA*RA
         IF(SR.LE.SRA) THEN
            QR=Q0+SR*(QA-Q0)/SRA
         ELSE
            QR=QA
         ENDIF
         BLO(3)= BB*RR/(RR+X)
         BLO(1)= Y*BLO(3)/(RR*QR*RKAP)
         BLO(2)=-X*BLO(3)/(RR*QR)
      ELSEIF(MODELB.EQ.6) THEN
         XH1=H1*X
         YH1=H1*Y
         RG=2.D0*BB*HA1/(1.D0+XH1**2+YH1**2)
         BLO(1)=  -RG*YH1*(1.D0+2.D0*XH1**2)
         BLO(2)=  -RG*XH1*(1.D0+2.D0*YH1**2)
         BLO(3)=BB+RG*(XH1**2-YH1**2)
      ELSEIF(MODELB.EQ.7) THEN
         BLO(1)=0.D0
         BLO(2)=0.D0
         BLO(3)=0.D0
         DO NC=1,NCMAX
            CALL WFCOIL(X,Y-ZC(NC),RC(NC),BX,BY)
            BLO(1)=BLO(1)+BC(NC)*BX
            BLO(2)=BLO(2)+BC(NC)*BY
         ENDDO
      ENDIF
C
      BABS=0.D0
      DO 10 I=1,3
         BABS=BABS+BLO(I)*BLO(I)
   10 CONTINUE
C
      IF(BABS.GT.0.D0) THEN
         BABS=SQRT(BABS)
         DO 20 I=1,3
            AL(I)=BLO(I)/BABS
   20    CONTINUE
      ELSE
         DO 30 I=1,3
            AL(I)=0.D0
   30    CONTINUE
      ENDIF
C
      RETURN
      END
C
C     ****** MAGNETIC FLUX PROFILE ******
C
      SUBROUTINE WFBPSI(X,Y,PSI)
C
      USE libbes,ONLY: besinx
      INCLUDE 'wfcomm.inc'
C
      IF(MODELB.EQ.0) THEN
         PSI= BB*Y
      ELSEIF(MODELB.EQ.1) THEN
         PSI=-BB*X
      ELSEIF(MODELB.EQ.2) THEN
         PSI=0.D0
      ELSEIF(MODELB.EQ.3) THEN
         A0=0.5D0*(1.D0+RMIR)*BB
         A1=0.5D0*(1.D0-RMIR)*BB
         RL = PI* X/ZBB
         ZL = PI* Y/ZBB
         PSI = 0.5D0*A0*X *X +A1*(ZBB/PI)**2*COS(ZL)*RL *BESINX(1,RL )
      ELSEIF(MODELB.EQ.4) THEN
         A0=0.5D0*(1.D0+RMIR)*BB
         A1=0.5D0*(1.D0-RMIR)*BB
         RL = PI* X/ZBB
         ZL = PI* Y/ZBB
         PSI = A0*X +A1*(ZBB/PI)*COS(ZL)*SINH(RL )
      ELSEIF(MODELB.EQ.5) THEN
         PSI=(X*X+Y*Y/(RKAP*RKAP))/(RA*RA)
      ELSEIF(MODELB.EQ.6) THEN
         XH1=H1*X
         YH1=H1*Y
         PSI=(XH1**2+YH1**2+2.D0*HA1*(XH1**2-YH1**2))
      ELSEIF(MODELB.EQ.7) THEN
         PSI=0.D0
         DO NC=1,NCMAX
            IF(X.NE.0.D0) THEN
               PSI=PSI+BC(NC)*APSI(X,Y-ZC(NC),RC(NC))
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
C
C     ****** LOOP COIL ******
C
      SUBROUTINE WFCOIL(RL,ZL,RC,BR,BZ)
C
      IMPLICIT REAL*8(A-F,H,O-Z)
C
      IF(RC.LT.0.D0) WRITE(6,*) 'XX WFCOIL:',RL,ZL,RC
      IF(RL.EQ.0.D0) THEN
         BR=0.D0
         BZ=(RC/SQRT(RC**2+ZL**2))**3
      ELSE
         DZ=1.D-6
         DR=MIN(1.D-6,0.5D0*RL)
         BR= (APSI(RL,ZL+DZ,RC)-APSI(RL,ZL-DZ,RC))/(2.D0*RL*DZ)
         BZ=-(APSI(RL+DR,ZL,RC)-APSI(RL-DR,ZL,RC))
     &       /(2.D0*DR*RL)
      ENDIF
      RETURN
      END
C
      FUNCTION APSI(RL,ZL,RC)
C        R*A_psi
      USE libell,ONLY: ELLFC,ELLEC
      IMPLICIT REAL*8(A-F,H,O-Z)
      DATA PI/3.1415926D0/
C
      RX=SQRT(RC**2+RL**2+ZL**2+2.D0*RC*RL)
      RK=SQRT(4.D0*RC*RL)/RX
      IF(RK.LT.0.D0) WRITE(6,*) 'APSI:RK:',RK,RL,ZL,RC
      APSI=(RC/PI)*RX
     &     *((1.D0-0.5D0*RK**2)*ELLFC(RK,IERR1)-ELLEC(RK,IERR2))
      RETURN
      END
C
C     ****** DENSITY & TEMPERATURE PROFILE ******
C
      SUBROUTINE WFSDEN(IN,RN,RTPR,RTPP,RZCL)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION RN(NSM),RT(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
      DIMENSION RNUE(NSM),RNUI(NSM),RNUN(NSM)
C
      IF(NEVOL.EQ.0) THEN
         CALL WFSPSI(XD(IN),YD(IN),PSI)
         IF(MODELP.EQ.0) THEN
            IF(PSI.LT.1.D0) THEN
               FACT=1.D0
            ELSE
               FACT=0.D0
            ENDIF
         ELSEIF(MODELP.EQ.1) THEN
            IF(PSI.LT.1.D0) THEN
               FACT=1.D0-PSI
            ELSE
               FACT=0.D0
            ENDIF
         ELSEIF(MODELP.EQ.2) THEN
            IF(YD(IN).GT.0.D0) THEN
               FACT=-999.D0
            ELSE
               IF(PSI.LT.1.D0) THEN
                  FACT=EXP(YD(IN)/ZBB)*(1.D0-PSI)
               ELSE
                  FACT=0.D0
               ENDIF
            ENDIF
         ELSEIF(MODELP.EQ.3) THEN
            IF(PSI.LT.1.D0) THEN
               Y=YD(IN)
               FACTY=1.D0-Y*Y/(ZBB*ZBB)
               IF(FACTY.GT.0.D0) THEN
                  FACT=(1.D0-PSI)*FACTY
               ELSE
                  FACT=0.D0
               ENDIF
            ELSE
               FACT=0.D0
            ENDIF
         ELSEIF(MODELP.EQ.4) THEN
            FACT=EXP(-(0.16D0-YD(IN))/0.10D0)
         ELSEIF(MODELP.EQ.5) THEN
            IF(PSI.LT.1.D0) THEN
               Y=YD(IN)
               FACTY=1.D0-Y**4/ZBB**4
               IF(FACTY.GT.0.D0) THEN
                  FACT=(1.D0-PSI)*FACTY
               ELSE
                  FACT=0.D0
               ENDIF
            ELSE
               FACT=0.D0
            ENDIF
         ELSEIF(MODELP.EQ.6) THEN
            FACT=1.D0-((XD(IN)-RA)**2+YD(IN)**2)/ZBB**2
            IF(FACT.LT.0.D0) THEN
               FACT=0.D0
            ENDIF
         ELSEIF(MODELP.EQ.7) THEN
            FACTX=1.D0-(XD(IN)/RA)**2
            FACTY=EXP(-(YD(IN)-PNPROFY0)**2/PNPROFYW**2)
            FACT=FACTX*FACTY
            IF(FACT.LT.0.D0) THEN
               FACT=0.D0
            ENDIF
         ELSE
            WRITE(6,*) 'XX WFSDEN: UNKNOWN MODELP: ',MODELP
         ENDIF
C
         DO NS=1,NSMAX
            IF(FACT.EQ.-999.D0) THEN
               RN(NS)=0.D0
            ELSE
               RN(NS)  =(PN(NS)  -PNS(NS))*FACT+PNS(NS)
            ENDIF
            RTPR(NS)=(PTPR(NS)-PTS(NS))*FACT+PTS(NS)
            RTPP(NS)=(PTPP(NS)-PTS(NS))*FACT+PTS(NS)
         ENDDO
         PNE(IN)=RN(1)
         PNI(IN)=RN(2)
         PTE(IN)=RTPR(1)*1.D3
         PTI(IN)=RTPR(2)*1.D3
      ELSE
         RN(1)=PNE(IN)
         RN(2)=PNI(IN)
         RTPR(1)=PTE(IN)*1.D-3
         RTPR(2)=PTI(IN)*1.D-3
         RTPP(1)=PTE(IN)*1.D-3
         RTPP(2)=PTI(IN)*1.D-3
      ENDIF
C
      IF(RN(1).GT.0.D0) THEN
         DO NS=1,NSMAX
            RT(NS)=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
         ENDDO
C
         CALL CALRNU(RN,RT,RNUE,RNUI,RNUN)
C
         DO NS=1,NSMAX
            RNUG(IN,NS)=RNUE(NS)+RNUI(NS)+RNUN(NS)
         ENDDO
C
         DO NS=1,NSMAX
            IF(PZCL(NS).EQ.0.D0) THEN
               RZCL(NS)=RNUG(IN,NS)/(2.D6*PI*RF)
            ELSE
               RZCL(NS)=PZCL(NS)
            ENDIF
         ENDDO
      ELSE
         DO NS=1,NSMAX
            RZCL(NS)=0.D0
         ENDDO
      ENDIF
C
      DO NS=1,NSMAX
         RZCLG(IN,NS)=RZCL(NS)
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE COLLISION FREQUENCY ******
C
      SUBROUTINE CALRNU(RN,RT,RNUE,RNUI,RNUN)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION RN(NSM),RT(NSM),RNUE(NSM),RNUI(NSM),RNUN(NSM)
C
      DO NS=1,NSMAX
         IF(NS.EQ.1) THEN
           IF(RN(1).GT.0) THEN
             RLAMEE= 8.0D0+2.3D0*(LOG10(RT(1))-0.5D0*LOG10(RN(1)))
             RLAMEI= RLAMEE+0.3D0
           ELSE
             RLAMEE= 1.D0
             RLAMEI= 1.D0
           ENDIF
         ELSE
           IF(RN(1).GT.0) THEN
             RLAMII=12.1D0+2.3D0*(LOG10(RT(NS))-0.5D0*LOG10(RN(1)))
           ELSE
             RLAMII= 1.D0
           ENDIF
         ENDIF
      ENDDO
C
      IF(MODELN.EQ.0) THEN
         TE=RT(1)
         VTE=SQRT(2.D0*TE*AEE/AME)
         SNE=0.88D-20
         SVE=SNE*VTE
      ELSE
         CALL ATSIGV(RT(1),SVE,1)
      ENDIF
      SNI=1.D-20
      PNN0=PPN0/(PTN0*AEE)
C
      DO NS=1,NSMAX
         IF(NS.EQ.1) THEN
            TE=RT(1)
            VTE=SQRT(2.D0*TE*AEE/AME)
            RNUEE=RN(1)*RLAMEE
     &           /(1.24D-4*SQRT(TE*1.D-3)**3)
            RNUEI=0.D0
            DO NSI=2,NSMAX
               RNUEI=RNUEI+PZ(NSI)**2*RN(NSI)
            ENDDO
            RNUEI=RNUEI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
            RNUEN=PNN0*SVE
            RNUE(NS)=RNUEE
            RNUI(NS)=RNUEI
            RNUN(NS)=RNUEN
         ELSE
            TI=RT(NS)
            VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
            RNUIE=PZ(NS)**2*RN(1)*RLAMEI
     &           /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
            RNUII=0.D0
            DO NSI=2,NSMAX
               RNUII=RNUII+PZ(NSI)**2*RN(NSI)
            ENDDO
            RNUII=RNUII*RLAMII
     &           /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
            RNUIN=PNN0*SNI*0.88D0*VTI
            RNUE(NS)=RNUIE
            RNUI(NS)=RNUII
            RNUN(NS)=RNUIN
         ENDIF
      ENDDO
C
      RETURN
      END
