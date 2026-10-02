C     $Id$
C
C     ********** WF TRANSPORT SOLVER ( FIRST ORDER ) **********
C
C     ****** INITIALIZATION ******
C
      SUBROUTINE WFEVIN
C
      INCLUDE 'wfcomm.inc'
C
      DO 10 IN=1,NNOD
         PHI(IN)=0.D0
         PSQ(IN)=0.D0
         IF(KNODP(IN).EQ.0) THEN
            CALL WFSPSI(XD(IN),YD(IN),PSI)
            IF(PSI.LT.1.D0) THEN
               IF(MODELP.EQ.1) THEN
                  FACT=1.D0-PSI
               ELSEIF(MODELP.EQ.0) THEN
                  FACT=1.D0
               ELSE
                  WRITE(6,*) 'XX WFEVIN: UNKNOWN MODELP: ',MODELP
               ENDIF
            ELSE
               FACT=0.D0
            ENDIF
            PNE(IN)=(PNE0-PNES)*FACT+PNES
            PTE(IN)=PTE0
            PPE(IN)=PNE(IN)*PTE0
            PNI(IN)=PNE(IN)/PZ(2)
            PTI(IN)=PTI0
            PPI(IN)=PNI(IN)*PTI0
         ELSEIF(KNODP(IN).GE.1.AND.KNODP(IN).LE.9) THEN
            PNE(IN)=PNES
            PTE(IN)=PTES
            PPE(IN)=PNES*PTES
            PNI(IN)=PNES/PZ(2)
            PTI(IN)=PTIS
            PPI(IN)=PNES/PZ(2)*PTIS
         ELSE
            PNE(IN)=0.D0
            PTE(IN)=0.D0
            PPE(IN)=0.D0
            PNI(IN)=0.D0
            PTI(IN)=0.D0
            PPI(IN)=0.D0
         ENDIF
         CEF(1,IN)=(0.D0,0.D0)
         CEF(2,IN)=(0.D0,0.D0)
         CEF(3,IN)=(0.D0,0.D0)
         PWRNS(IN,1)=0.D0
         PWRNS(IN,2)=0.D0
   10 CONTINUE
C
C      PNE(NNOD/2)=10.D0*PNE0
C      PNI(NNOD/2)=10.D0*PNE0/PZ(2)
C
      NGTMAX=0
C
      RETURN
      END
C
C     ***** ONE TIME-STEP ROUTINE *****
C
      SUBROUTINE WFEVOL
C
      INCLUDE 'wfcomm.inc'
C
         DO IN=1,NNOD
            PHIPRE(IN)=PHI(IN)
            PNEPRE(IN)=PNE(IN)
            PPEPRE(IN)=PPE(IN)
            PNIPRE(IN)=PNI(IN)
            PPIPRE(IN)=PPI(IN)
            PSQPRE(IN)=PSQ(IN)
            VDE(IN)=0.D0
            VDI(IN)=0.D0
         ENDDO
C
      DO LOOP=1,MAXIMP
         DO IN=1,NNOD
            PHISAV(IN)=PHI(IN)
            PNESAV(IN)=PNE(IN)
            PPESAV(IN)=PPE(IN)
            PNISAV(IN)=PNI(IN)
            PPISAV(IN)=PPI(IN)
            PSQSAV(IN)=PSQ(IN)
         ENDDO
C
         DO IN=1,NNOD
            PHI(IN)=FACIMP*PHI(IN)+(1.D0-FACIMP)*PHIPRE(IN)
            PNE(IN)=FACIMP*PNE(IN)+(1.D0-FACIMP)*PNEPRE(IN)
            PPE(IN)=FACIMP*PPE(IN)+(1.D0-FACIMP)*PPEPRE(IN)
            PNI(IN)=FACIMP*PNI(IN)+(1.D0-FACIMP)*PNIPRE(IN)
            PPI(IN)=FACIMP*PPI(IN)+(1.D0-FACIMP)*PPIPRE(IN)
            PSQ(IN)=FACIMP*PSQ(IN)+(1.D0-FACIMP)*PSQPRE(IN)
         ENDDO
C
         CALL DVDNOD(IERR)
            IF(IERR.NE.0) GOTO 9000
C
         CALL DVSOLV(IERR)
            IF(IERR.EQ.9002.OR.IERR.EQ.30000) GOTO 9000
            IF(IERR.NE.0) GOTO 9000
C
         CALL DALFLD
C
         DO IN=1,NNOD
            PHIDIF(IN)=PHI(IN)-PHISAV(IN)
            PNEDIF(IN)=PNE(IN)-PNESAV(IN)
            PPEDIF(IN)=PPE(IN)-PPESAV(IN)
            PNIDIF(IN)=PNI(IN)-PNISAV(IN)
            PPIDIF(IN)=PPI(IN)-PPISAV(IN)
            PSQDIF(IN)=PSQ(IN)-PSQSAV(IN)
         ENDDO
C
         PHISUM1=0.D0
         PNESUM1=0.D0
         PPESUM1=0.D0
         PNISUM1=0.D0
         PPISUM1=0.D0
         PSQSUM1=0.D0
         PHISUM2=0.D0
         PNESUM2=0.D0
         PPESUM2=0.D0
         PNISUM2=0.D0
         PPISUM2=0.D0
         PSQSUM2=0.D0
         DO IN=1,NNOD
            PHISUM1=PHISUM1+PHIDIF(IN)**2
            PHISUM2=PHISUM2+PHISAV(IN)**2
            PNESUM1=PNESUM1+PNEDIF(IN)**2
            PNESUM2=PNESUM2+PNESAV(IN)**2
            PPESUM1=PPESUM1+PPEDIF(IN)**2
            PPESUM2=PPESUM2+PPESAV(IN)**2
            PNISUM1=PNISUM1+PNIDIF(IN)**2
            PNISUM2=PNISUM2+PNISAV(IN)**2
            PPISUM1=PPISUM1+PPIDIF(IN)**2
            PPISUM2=PPISUM2+PPISAV(IN)**2
            PSQSUM1=PSQSUM1+PSQDIF(IN)**2
            PSQSUM2=PSQSUM2+PSQSAV(IN)**2
         ENDDO
         PHISUM=PHISUM1/(NNOD*MAX(1.D0,PHISUM2))
         PNESUM=PNESUM1/(NNOD*MAX(1.D0,PNESUM2))
         PPESUM=PPESUM1/(NNOD*MAX(1.D0,PPESUM2))
         PNISUM=PNISUM1/(NNOD*MAX(1.D0,PNISUM2))
         PPISUM=PPISUM1/(NNOD*MAX(1.D0,PPISUM2))
         PSQSUM=PSQSUM1/(NNOD*MAX(1.D0,PSQSUM2))
C
         IF(MAXIMP.GT.1) 
     &      WRITE(6,601) LOOP,PHISUM,PNESUM,PPESUM,
     &                        PNISUM,PPISUM,PSQSUM
  601    FORMAT(' # IMPL LOOP =',I3,1P6E10.2)
C
         SUMMAX=MAX(PHISUM,PNESUM,PPESUM,PNISUM,PPISUM,PSQSUM)
         IF(SUMMAX.LE.EPSIMP) GOTO 9000
      ENDDO
C
 9000 RETURN
      END
C
C     ****** STORE TIME-DEPENDENT VARIABLES ******
C
      SUBROUTINE WFEVST
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION XE(3),YE(3),RL(3)
      REAL*8 A(3),B(3),C(3)
C
      PNETOT=0.D0
      PNITOT=0.D0
      PPETOT=0.D0
      PPITOT=0.D0
      PHITOT=0.D0
      VTOT  =0.D0
C
      DO IEDO=1,NELM
         IE=IEDO
         IF(KELM(IE).EQ.0) THEN
            CALL WFNPOS(IE,XE,YE)
            CALL WFABC(IE,A,B,C,S)
            DO K=1,3
               IF(MODELS.EQ.1) THEN
                  RL(K)=XE(K)
               ELSEIF(MODELS.EQ.2) THEN
                  RL(K)=RR+XE(K)
               ELSE
                  RL(K)=1.D0
               ENDIF
            DO N=1,3
               IN=IELM(N,IE)
               PNETOT=PNETOT+RL(K)*S*PNE(IN)*AIF2(N,K)
               PNITOT=PNITOT+RL(K)*S*PNI(IN)*AIF2(N,K)
               PPETOT=PPETOT+RL(K)*S*PPE(IN)*AIF2(N,K)
               PPITOT=PPITOT+RL(K)*S*PPI(IN)*AIF2(N,K)
               PHITOT=PHITOT+RL(K)*S*PHI(IN)*AIF2(N,K)
               VTOT  =VTOT  +RL(K)*S        *AIF2(N,K)
            ENDDO
            ENDDO
         ENDIF
      ENDDO
C
      PNEAVE=PNETOT/VTOT
      PNIAVE=PNITOT/VTOT
      PPEAVE=PPETOT/VTOT
      PPIAVE=PPITOT/VTOT
      PTEAVE=PPEAVE/PNEAVE
      PTIAVE=PPIAVE/PNIAVE
      PHIAVE=PHITOT/VTOT
C
      PHIMAX=PHI(1)
      PHIMIN=PHI(1)
      PNEMAX=PNE(1)
      PNEMIN=PNE(1)
      PNIMAX=PNE(1)
      PNIMIN=PNE(1)
      PTEMAX=PTE(1)
      PTEMIN=PTE(1)
      PTIMAX=PTI(1)
      PTIMIN=PTI(1)
      DO IN=2,NNOD
         PHIMAX=MAX(PHIMAX,PHI(IN))
         PHIMIN=MIN(PHIMIN,PHI(IN))
         PNEMAX=MAX(PNEMAX,PNE(IN))
         PNEMIN=MIN(PNEMIN,PNE(IN))
         PNIMAX=MAX(PNIMAX,PNI(IN))
         PNIMIN=MIN(PNIMIN,PNI(IN))
         PTEMAX=MAX(PTEMAX,PTE(IN))
         PTEMIN=MIN(PTEMIN,PTE(IN))
         PTIMAX=MAX(PTIMAX,PTI(IN))
         PTIMIN=MIN(PTIMIN,PTI(IN))
      ENDDO
C
      IF(NGTMAX.LE.NGTM) THEN
         NGTMAX=NGTMAX+1
      ENDIF
      TG(NGTMAX)=T
      PHIT(NGTMAX,1)=PHIAVE
      PHIT(NGTMAX,2)=PHIMAX
      PHIT(NGTMAX,3)=PHIMIN
      PNT(NGTMAX,1,1)=PNEAVE
      PNT(NGTMAX,2,1)=PNEMAX
      PNT(NGTMAX,3,1)=PNEMIN
      PNT(NGTMAX,1,2)=PNIAVE
      PNT(NGTMAX,2,2)=PNIMAX
      PNT(NGTMAX,3,2)=PNIMIN
      PTT(NGTMAX,1,1)=PTEAVE
      PTT(NGTMAX,2,1)=PTEMAX
      PTT(NGTMAX,3,1)=PTEMIN
      PTT(NGTMAX,1,2)=PTIAVE
      PTT(NGTMAX,2,2)=PTIMAX
      PTT(NGTMAX,3,2)=PTIMIN
      PABST(NGTMAX,1)=PWRS(1)
      PABST(NGTMAX,2)=PWRS(2)
      PABST(NGTMAX,3)=PWRS(1)+PWRS(2)
C
      WRITE(6,601) T,PNEAVE,PTEAVE,PHIMIN,PHIMAX
  601 FORMAT(1H ,'T=',1PE11.3,'  NAVE,TAVE,PHIMM=',1P4E11.3)
      RETURN
      END
C
C     ******* SET KNOD ARRAY *******
C
      SUBROUTINE DVDNOD(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      DO 10 IE=1,NELM
      DO 10 IN=1,3
         JELM(IN,IE)=IELM(IN,IE)
   10 CONTINUE
C
      IF(MODELW.EQ.0) THEN
         NNDN=0
         NNDT=0
      ELSEIF(MODELW.EQ.1) THEN
         NNDN=2
         NNDT=0
      ELSEIF(MODELW.EQ.2) THEN
         NNDN=2
         NNDT=2
      ELSE
         WRITE(6,*) 'XX INVALID MODELW: ',MODELW
      ENDIF
C
      IL=1
      DO 20 IN=1,NNOD
C     --- PLASMA ---
         IF(KNODP(IN).EQ.0) THEN
            K=5
C     --- PLASMA-WALL ---
         ELSEIF(KNODP(IN).LE.5) THEN
            K=NNDN+NNDT
C     --- PLAMSA-DIELECTRIC / WITHOUT CHARGE ---
         ELSEIF(KNODP(IN).EQ.8) THEN
            K=1+NNDN+NNDT
C     --- PLAMSA-DIELECTRIC / WITH SURFACE CHARGE ---
         ELSEIF(KNODP(IN).EQ.9) THEN
            K=1+NNDN+NNDT+1
C     --- DIELECTRIC ---
         ELSEIF(KNODP(IN).EQ.10) THEN
            K=1
C     --- DIELECTRIC-WALL ---
         ELSEIF(KNODP(IN).LE.15) THEN
            K=0
         ELSE
            WRITE(6,*) 'XX ERROR IN DVDNOD: UNDEFINED KNODP'
         END IF
         JBND(IN)=K
         IBND(IN)=IL
         IL=IL+K
   20 CONTINUE
C
      MLEN=IL-1
      IBND(NNOD+1)=IL
      IF(MLEN.GT.MLENM) GOTO 9000
      IERR=0
      RETURN
C
 9000 WRITE(6,*) 'XX DVDNOD ERROR : MLEN,MLENM =',MLEN,MLENM
      IERR=900
      RETURN
      END
C
C     ****** CALCULATE LOCAL COEFFICIENT MATRIX ******
C
      SUBROUTINE DMCALC(IE)
C
      INCLUDE 'wfcomm.inc'
C
      IF(KELM(IE).NE.0) THEN
         CALL DMCALD(IE)
      ELSE
         CALL DMCALP(IE)
      ENDIF
C
C     +++++ BOUNDARY CONDITIONS +++++
C
C        KNODP  = 0 : not on boundary
C                 1 : conductor boundary with phi = 0
C                 2 : conductor boundary with phi = PHIES
C                 3 : conductor boundary with phi = 0 - PHIES
C                 4 : conductor boundary with phi = PHIES sin omega t
C                 5 : conductor boundary with phi = 0 - PHIES sin omega t
C                 8 : insulator boundary without surface charge
C                 9 : insulator boundary with surface charge
C                10 : not on boundary (without plasma)
C                11 : conductor boundary with phi = 0
C                12 : conductor boundary with phi = PHIES
C                13 : conductor boundary with phi = 0 - PHIES
C                14 : conductor boundary with phi = PHIES sin omega t
C                15 : conductor boundary with phi = 0 - PHIES sin omega t
C
C        MODELW = 0 : fixed density and temperature
C                 1 : fixed density and free temperature
C                 2 : free density and temperature
C

      DO M=1,3
         IM=IELM(M,IE)
         KNODPL=KNODP(IM)
         IF(KNODPL.GE.10) KNODPL=KNODPL-10
         IF(KNODPL.GE.1.AND.KNODPL.LE.5) THEN
            DO N=1,3
               IF(MODELW.EQ.1) THEN
                  DO I=1,6
                     DM(I,1,N,M)=DM(I,2,N,M)
                     DM(I,2,N,M)=DM(I,4,N,M)
                  ENDDO
               ELSEIF (MODELW.EQ.2) THEN
                  DO I=1,6
                     DM(I,1,N,M)=DM(I,2,N,M)
                     DM(I,2,N,M)=DM(I,3,N,M)
                     DM(I,3,N,M)=DM(I,4,N,M)
                     DM(I,4,N,M)=DM(I,5,N,M)
                  ENDDO
               ENDIF
            ENDDO
         ELSEIF(KNODPL.GE.8) THEN
            DO N=1,3
               IF(MODELW.EQ.0) THEN
                  DO I=1,6
                     DM(I,2,N,M)=DM(I,6,N,M)
                  ENDDO
               ELSEIF(MODELW.EQ.1) THEN
                  DO I=1,6
                     DM(I,3,N,M)=DM(I,4,N,M)
                     DM(I,4,N,M)=DM(I,6,N,M)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C
      DO N=1,3
         IN=IELM(N,IE)
         KNODPL=KNODP(IN)
         IF(KNODPL.GE.10) KNODPL=KNODPL-10
         IF(KNODPL.GE.1.AND.KNODPL.LE.5) THEN
            IF(MODELW.EQ.1) THEN
               DO M=1,3
                  DO J=1,6
                     DM(1,J,N,M)=DM(2,J,N,M)
                     DM(2,J,N,M)=DM(4,J,N,M)
                  ENDDO
               ENDDO
               DV(1,N)=DV(2,N)
               DV(2,N)=DV(4,N)
            ELSEIF (MODELW.EQ.2) THEN
               DO M=1,3
                  DO J=1,6
                     DM(1,J,N,M)=DM(2,J,N,M)
                     DM(2,J,N,M)=DM(3,J,N,M)
                     DM(3,J,N,M)=DM(4,J,N,M)
                     DM(4,J,N,M)=DM(5,J,N,M)
                  ENDDO
               ENDDO
               DV(1,N)=DV(2,N)
               DV(2,N)=DV(3,N)
               DV(3,N)=DV(4,N)
               DV(4,N)=DV(5,N)
            ENDIF
         ELSEIF(KNODPL.GE.8) THEN
            IF(MODELW.EQ.0) THEN
               DO M=1,3
                  DO J=1,6
                     DM(2,J,N,M)=DM(6,J,N,M)
                  ENDDO
               ENDDO
               DV(2,N)=DV(6,N)
            ELSEIF(MODELW.EQ.1) THEN
               DO M=1,3
                  DO J=1,6
                     DM(3,J,N,M)=DM(4,J,N,M)
                     DM(4,J,N,M)=DM(6,J,N,M)
                  ENDDO
               ENDDO
               DV(3,N)=DV(4,N)
               DV(4,N)=DV(6,N)
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END
C
C     ****** CALCULATE LOCAL COEFFICIENT MATRIX IN DIELECTRIC ******
C
      SUBROUTINE DMCALD(IE)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION XE(3),YE(3)
      REAL*8 A(3),B(3),C(3)
      DIMENSION RL(3)
     &     
C
      EPSD=EPSDM(KELM(IE))
C
      CALL WFNPOS(IE,XE,YE)
      CALL WFABC(IE,A,B,C,S)
C
      DO 1000 N=1,3
      DO 1000 M=1,3
      DO 1000 I=1,6
      DO 1000 J=1,6
         DM(I,J,N,M)=0.D0
 1000 CONTINUE
C
      DO 1100 N=1,3
      DO 1100 I=1,6
         DV(I,N)=0.D0
 1100 CONTINUE
C
      DO 2000 K=1,3
         IF(MODELS.EQ.1) THEN
            RL(K)=XE(K)
         ELSEIF(MODELS.EQ.2) THEN
            RL(K)=RR+XE(K)
         ELSE
            RL(K)=1.D0
         ENDIF
 2000 CONTINUE
C
      DO 3000 K=1,3
C
      DO 3000 N=1,3
      DO 3000 M=1,3
         DM(1,1,N,M)=DM(1,1,N,M)
     &              +S*RL(K)*AIF1(K)*( B(N)*B(M)
     &                                +C(N)*C(M))*EPSD
 3000 CONTINUE
C
C     +++++ BOUNDARY CONDITIONS +++++
C        KNOD = 10 : not on boundary
C               11 : conductor boundary with phi = 0
C               12 : conductor boundary with phi = PHIES
C               13 : conductor boundary with phi = 0 - PHIES
C               14 : conductor boundary with phi = PHIES sin omega t
C               15 : conductor boundary with phi = 0 - PHIES sin omega t
C
      DO M=1,3
         IM=IELM(M,IE)
         KNODPL=KNODP(IM)
         IF(KNODPL.GE.10) KNODPL=KNODPL-10
         IF(KNODPL.GE.1.AND.KNODPL.LE.5) THEN
            PHIL=WFBPHI(XD(IM),YD(IM),KNODPL)
            DO N=1,3
               DV(1,N)=DV(1,N)-DM(1,1,N,M)*PHIL
            ENDDO
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE GIVEN BOUNDARY POTENTIAL ******
C
      FUNCTION WFBPHI(XNL,YNL,KNODPL)
C
      INCLUDE 'wfcomm.inc'
C
      IF(KNODPL.EQ.1) THEN
         PHIL=0.D0
      ELSE IF(KNODPL.EQ.2) THEN
         PHIL=PHIES
      ELSE IF(KNODPL.EQ.3) THEN
         NV=1
         XMIN=XYVARP(1,NV)
         XMAX=XYVARP(2,NV)
         YMIN=XYVARP(3,NV)
         YMAX=XYVARP(4,NV)
         RLEN0=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
         RLEN1=SQRT((XNL-XMIN)**2+(YNL-YMIN)**2)
         FACT=RLEN1/RLEN0
         PHIL=PHIES*FACT
      ELSE IF(KNODPL.EQ.4) THEN
         PHIL=PHIES*SIN(2.D0*PI*RFES*1.D6*T)
      ELSE IF(KNODPL.EQ.5) THEN
         NV=1
         XMIN=XYVARP(1,NV)
         XMAX=XYVARP(2,NV)
         YMIN=XYVARP(3,NV)
         YMAX=XYVARP(4,NV)
         RLEN0=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
         RLEN1=SQRT((XNL-XMIN)**2+(YNL-YMIN)**2)
         FACT=RLEN1/RLEN0
         PHIL=PHIES*FACT*SIN(2.D0*PI*RFES*1.D6*T)
      ENDIF
      WFBPHI=PHIL
      RETURN
      END
C
C     ****** CALCULATE LOCAL COEFFICIENT MATRIX ******
C
      SUBROUTINE DMCALP(IE)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION XE(3),YE(3)
      DIMENSION AMK(2,2,NSM),ADK(2,2,NSM),AHK(2,2,NSM)
      DIMENSION ANK(MFM,MFM),APK(MFM),AVT(2,NSM)
      DIMENSION VDEL(2),VDIL(2)
      REAL*8 A(3),B(3),C(3)
      DIMENSION RL(3),IBX(2)
      DIMENSION DL(MFM,MFM,3,3),DD(MFM,3,3),DS(MFM,3),V(MFM,3)
C
      CALL WFNPOS(IE,XE,YE)
      CALL WFABC(IE,A,B,C,S)
C
      DO M=1,3
      DO N=1,3
      DO J=1,6
      DO I=1,6
         DL(I,J,N,M)=0.D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO M=1,3
      DO N=1,3
      DO I=1,6
         DD(I,N,M)=0.D0
      ENDDO
      ENDDO
      ENDDO
C
      DO N=1,3
      DO I=1,6
         DS(I,N)=0.D0
      ENDDO
      ENDDO
C
      DTI=1.D0/DT
C
      DO K=1,3
         IF(MODELS.EQ.1) THEN
            RL(K)=XE(K)
         ELSEIF(MODELS.EQ.2) THEN
            RL(K)=RR+XE(K)
         ELSE
            RL(K)=1.D0
         ENDIF
      ENDDO
C
      DO K=1,3
         IK=IELM(K,IE)
         PTEK=PTE(IK)
         PTIK=PTI(IK)
C
         IF(PTEK.LE.0.D0) THEN
            FACTE=0.D0
         ELSE
            FACTE=1.D0/PTEK
         ENDIF
C
         IF(PTIK.LE.0.D0) THEN
            FACTI=0.D0
         ELSE
            FACTI=1.D0/PTIK
         ENDIF
C
         CALL WFCOEF(K,IE,AMK,ADK,AHK,ANK,APK,AVT)
C
      DO N=1,3
      DO M=1,3
         IM=IELM(M,IE)
         FACT1=S*RL(K)*AIF1(K)
         FACT3=S*RL(K)*AIF3(N,M,K)
C
         DD(1,N,M)=0.D0
         DL(1,1,N,M)=DL(1,1,N,M)
     &              +FACT1*(B(N)*B(M)+C(N)*C(M))
         DL(1,2,N,M)=DL(1,2,N,M)
     &              -FACT3*PZ(1)*AEE*1.D20/EPS0
         DL(1,4,N,M)=DL(1,4,N,M)
     &              -FACT3*PZ(2)*AEE*1.D20/EPS0
C
         DDME= FACT1*( B(N)*AMK(1,1,1)*B(M)
     &                +B(N)*AMK(1,2,1)*C(M)
     &                +C(N)*AMK(2,1,1)*B(M)
     &                +C(N)*AMK(2,2,1)*C(M))
         DDDE= FACT1*( B(N)*ADK(1,1,1)*B(M)
     &                +B(N)*ADK(1,2,1)*C(M)
     &                +C(N)*ADK(2,1,1)*B(M)
     &                +C(N)*ADK(2,2,1)*C(M))
         DDHE= FACT1*( B(N)*AHK(1,1,1)*B(M)
     &                +B(N)*AHK(1,2,1)*C(M)
     &                +C(N)*AHK(2,1,1)*B(M)
     &                +C(N)*AHK(2,2,1)*C(M))
C
         DD(2,N,M)  =DD(2,N,M)  +FACT3*DTI
         DL(2,1,N,M)=DL(2,1,N,M)-DDME
         DL(2,2,N,M)=DL(2,2,N,M)+FACT3*ANK(2,2)
C
         IF(MODELT.LE.1) THEN
            DL(2,2,N,M)=DL(2,2,N,M)-DDDE
         ELSE IF(MODELT.EQ.2) THEN
            DL(2,3,N,M)=DL(2,3,N,M)-DDDE*FACTE
         ENDIF
C
         DD(3,N,M)  =DD(3,N,M)     +FACT3*DTI     *1.5D0
         IF(MODELT.EQ.0) THEN
            DL(3,1,N,M)=DL(3,1,N,M)-DDME          *1.5D0*PTEK
            DL(3,2,N,M)=DL(3,2,N,M)+FACT3*ANK(2,2)*1.5D0*PTEK
     &                             -DDDE          *1.5D0*PTEK
         ELSE
            DL(3,1,N,M)=DL(3,1,N,M)-DDME          *2.5D0*PTEK
            DL(3,2,N,M)=DL(3,2,N,M)+FACT3*ANK(3,2)
C
            IF(MODELT.EQ.1) THEN
               DL(3,2,N,M)=DL(3,2,N,M)-DDDE       *2.5D0*PTEK
            ELSE IF(MODELT.EQ.2) THEN
               DL(3,3,N,M)=DL(3,3,N,M)-DDDE       *2.5D0
            ENDIF
C
            DL(3,2,N,M)=DL(3,2,N,M)+DDHE          *1.5D0*PTEK
            DL(3,3,N,M)=DL(3,3,N,M)+FACT3*ANK(3,3)*1.5D0
     &                             -DDHE          *1.5D0
            DL(3,5,N,M)=DL(3,5,N,M)+FACT3*ANK(3,5)*1.5D0
         ENDIF
C
         DDMI= FACT1*( B(N)*AMK(1,1,2)*B(M)
     &                +B(N)*AMK(1,2,2)*C(M)
     &                +C(N)*AMK(2,1,2)*B(M)
     &                +C(N)*AMK(2,2,2)*C(M))
         DDDI= FACT1*( B(N)*ADK(1,1,2)*B(M)
     &                +B(N)*ADK(1,2,2)*C(M)
     &                +C(N)*ADK(2,1,2)*B(M)
     &                +C(N)*ADK(2,2,2)*C(M))
         DDHI= FACT1*( B(N)*AHK(1,1,2)*B(M)
     &                +B(N)*AHK(1,2,2)*C(M)
     &                +C(N)*AHK(2,1,2)*B(M)
     &                +C(N)*AHK(2,2,2)*C(M))
C
         DD(4,N,M)  =DD(4,N,M)  +FACT3*DTI
         DL(4,1,N,M)=DL(4,1,N,M)-DDMI
         DL(4,2,N,M)=DL(4,2,N,M)+FACT3*ANK(4,2)
C
         IF(MODELT.LE.1) THEN
            DL(4,4,N,M)=DL(4,4,N,M)-DDDI
         ELSE
            DL(4,5,N,M)=DL(4,5,N,M)-DDDI*FACTI
         ENDIF
C
         DD(5,N,M)  =DD(5,N,M)     +FACT3*DTI     *1.5D0
         IF(MODELT.EQ.0) THEN
            DL(5,1,N,M)=DL(5,1,N,M)-DDMI          *1.5D0*PTIK
            DL(5,2,N,M)=DL(5,2,N,M)+FACT3*ANK(4,2)*1.5D0*PTIK
            DL(5,4,N,M)=DL(5,4,N,M)-DDDI          *1.5D0*PTIK
         ELSE
            DL(5,1,N,M)=DL(5,1,N,M)-DDMI          *2.5D0*PTIK
            DL(5,3,N,M)=DL(5,3,N,M)+FACT3*ANK(5,3)*1.5D0
C
            IF(MODELT.EQ.1) THEN
               DL(5,4,N,M)=DL(5,4,N,M)-DDDI       *2.5D0*PTIK
            ELSEIF(MODELT.EQ.2) THEN
               DL(5,5,N,M)=DL(5,5,N,M)-DDDI       *2.5D0
            ENDIF
C
            DL(5,4,N,M)=DL(5,4,N,M)+DDHI          *1.5D0*PTIK
            DL(5,5,N,M)=DL(5,5,N,M)+FACT3*ANK(5,5)*1.5D0
     &                             -DDHI          *1.5D0
         ENDIF
C
         DS(2,N)=DS(2,N)   +FACT3*APK(2)
         DS(4,N)=DS(4,N)   +FACT3*APK(4)
         IF(MODELT.EQ.0) THEN
            DS(3,N)=DS(3,N)+FACT3*APK(2)*1.5D0*PTEK
            DS(5,N)=DS(5,N)+FACT3*APK(4)*1.5D0*PTIK
         ELSE
            DS(3,N)=DS(3,N)+FACT3*APK(3)
            DS(5,N)=DS(5,N)+FACT3*APK(5)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      GOTO 9999
C
C     +++++ If two nodes are on boundary, set boundary condition +++++
C
      DO NM=1,3
         INM=IELM(NM,IE)
         IF(KNODP(INM).GE.1) THEN
            NP=MOD(NM,3)+1
            INP=IELM(NP,IE)
            IF(KNODP(INP).GE.1) THEN
C
C     +++++ When the side is boundary +++++
C
               X1=XD(INM)
               Y1=YD(INM)
               X2=XD(INP)
               Y2=YD(INP)
               IBX(1)=NM
               IBX(2)=NP
C
               RX=SQRT((X2-X1)**2+(Y2-Y1)**2)
               ANX= (Y2-Y1)/RX
               ANY=-(X2-X1)/RX
C
               DO NN=1,2
                  VDEL(NN)=0.D0
                  VDIL(NN)=0.D0
               ENDDO
C
               DO KK=1,2
                  K=IBX(KK)
                  CALL WFCOEF(K,IE,AMK,ADK,AHK,ANK,APK,AVT)
                  IK=IELM(K,IE)
                  PTEK=PTE(IK)
                  PTIK=PTI(IK)
                  IF(PTEK.LE.0.D0) THEN
                     FACTE=0.D0
                  ELSE
                     FACTE=1.D0/PTEK
                  ENDIF
                  IF(PTIK.LE.0.D0) THEN
                     FACTI=0.D0
                  ELSE
                     FACTI=1.D0/PTIK
                  ENDIF
C
                  DO NN=1,2
                     N=IBX(NN)
                     IN=IELM(N,IE)
                     FACT2=RX*RL(K)*AIE2(NN,KK)
C
C     +++++ calculation of surface charge sigma +++++
C
                     DO M=1,3
                        VDEL(NN)=VDEL(NN)
     &                      -FACT2*(ANX*AMK(1,1,1)*B(M)
     &                             +ANX*AMK(1,2,1)*C(M)
     &                             +ANY*AMK(2,1,1)*B(M)
     &                             +ANY*AMK(2,2,1)*C(M))
     &                            *PHI(IM)
                        IF(MODELT.LE.1) THEN
                           VDEL(NN)=VDEL(NN)
     &                         -FACT2*(ANX*ADK(1,1,1)*B(M)
     &                                +ANX*ADK(1,2,1)*C(M)
     &                                +ANY*ADK(2,1,1)*B(M)
     &                                +ANY*ADK(2,2,1)*C(M))
     &                               *PNE(IM)
                        ELSE
                           VDEL(NN)=VDEL(NN)
     &                         -FACT2*(ANX*ADK(1,1,1)*B(M)
     &                                +ANX*ADK(1,2,1)*C(M)
     &                                +ANY*ADK(2,1,1)*B(M)
     &                                +ANY*ADK(2,2,1)*C(M))
     &                               *PNE(IM)*PTE(IM)*FACTE
                        ENDIF
C
                        VDIL(NN)=VDIL(NN)
     &                      -FACT2*(ANX*AMK(1,1,2)*B(M)
     &                             +ANX*AMK(1,2,2)*C(M)
     &                             +ANY*AMK(2,1,2)*B(M)
     &                             +ANY*AMK(2,2,2)*C(M))
     &                            *PHI(IM)
                        IF(MODELT.LE.1) THEN
                           VDIL(NN)=VDIL(NN)
     &                         -FACT2*(ANX*ADK(1,1,2)*B(M)
     &                                +ANX*ADK(1,2,2)*C(M)
     &                                +ANY*ADK(2,1,2)*B(M)
     &                                +ANY*ADK(2,2,2)*C(M))
     &                               *PNI(IM)
                        ELSE
                           VDIL(NN)=VDIL(NN)
     &                         -FACT2*(ANX*ADK(1,1,2)*B(M)
     &                                +ANX*ADK(1,2,2)*C(M)
     &                                +ANY*ADK(2,1,2)*B(M)
     &                                +ANY*ADK(2,2,2)*C(M))
     &                               *PNI(IM)*PTI(IM)*FACTI
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NN=1,2
                  N=IBX(NN)
                  IN=IELM(N,IE)
                  VDE(IN)=VDEL(NN)
                  VDI(IN)=VDIL(NN)
                  IF(VDEL(NN).LT.0.D0) VDEL(NN)=0.D0
                  IF(VDIL(NN).LT.0.D0) VDIL(NN)=0.D0
C
C     +++++ calculation of surface charge +++++
C
                  IF(KNODP(IN).EQ.9) THEN
                     DS(6,N)=DS(6,N)+PZ(1)*AEE*1.D20*VDEL(NN)
     &                              +PZ(2)*AEE*1.D20*VDIL(NN)
C
                     DO KK=1,2
                        K=IBX(KK)
                        DO MM=1,2
                           M=IBX(MM)
                           FACT3=RX*RL(K)*AIE3(NN,KK,MM)
                           DD(6,N,M)=DD(6,N,M)+FACT3*DTI
                           DS(6,N)  =DS(6,N)
     &                              +FACT3*PZ(1)*AEE*1.D20*VDEL(NN)
     &                              +FACT3*PZ(2)*AEE*1.D20*VDIL(NN)
C
C     +++++ effect of surface charge +++++
C
                           DL(1,6,N,M)=DL(1,6,N,M)-FACT3/EPS0
C
                        ENDDO
                     ENDDO
                  ENDIF
C
C     +++++ Surface integral of particle and heat fluxes +++++
C
                  IF(MODELW.GE.1) THEN
                     DS(2,N)=DS(2,N)+VDEL(NN)
                     DS(4,N)=DS(4,N)+VDIL(NN)
                  ENDIF
                  IF(MODELW.EQ.2) THEN
                     IF(MODELT.EQ.0) THEN
                        DS(3,N)=DS(3,N)+VDEL(NN)*1.5D0*PTE(IN)
                        DS(5,N)=DS(5,N)+VDIL(NN)*1.5D0*PTI(IN)
                     ELSE
                        DS(3,N)=DS(3,N)+VDEL(NN)*2.5D0*PTE(IN)
                        DS(5,N)=DS(5,N)+VDIL(NN)*2.5D0*PTI(IN)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C
 9999 CONTINUE
C
C     +++++ Assemble to DM and DV +++++
C
      DO M=1,3
         IM=IELM(M,IE)
         V(1,M)=PHIPRE(IM)
         V(2,M)=PNEPRE(IM)
         V(3,M)=PPEPRE(IM)
         V(4,M)=PNIPRE(IM)
         V(5,M)=PPIPRE(IM)
         V(6,M)=PSQPRE(IM)
      ENDDO
C
      DO N=1,3
      DO I=1,6
         DV(I,N)=DS(I,N)
         DO M=1,3
            DO J=1,6
               DM(I,J,N,M)=0.D0
            ENDDO
            DM(I,I,N,M)=DM(I,I,N,M)+DD(I,N,M)
            DV(I,N)    =DV(I,N)    +DD(I,N,M)*V(I,M)
            DO J=1,6
               DM(I,J,N,M)=DM(I,J,N,M)  -FACIMP *DL(I,J,N,M)
               DV(I,N)    =DV(I,N)+(1.D0-FACIMP)*DL(I,J,N,M)*V(J,M)
            ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     +++++ FIXED BOUNDARY CONDITIONS +++++
C
      DO M=1,3
         IM=IELM(M,IE)
C
         IF(KNODP(IM).GE.1.AND.KNODP(IM).LE.5) THEN
            PHIL=WFBPHI(XD(IM),YD(IM),KNODP(IM))
            DO N=1,3
               DO I=1,6
                  DV(I,N)=DV(I,N)-DM(I,1,N,M)*PHIL
               ENDDO
            ENDDO
         ENDIF
C
         IF(KNODP(IM).GE.1) THEN
            IF(MODELW.EQ.0) THEN
               DO N=1,3
                  DO I=1,6
                     DV(I,N)=DV(I,N)-DM(I,2,N,M)*PNES
                     DV(I,N)=DV(I,N)-DM(I,4,N,M)*PNES/PZ(2)
                  ENDDO
               ENDDO
            ENDIF
C
            IF(MODELW.LE.1) THEN
               DO N=1,3
                  DO I=1,6
                     DV(I,N)=DV(I,N)-DM(I,3,N,M)*PNES*PTES
                     DV(I,N)=DV(I,N)-DM(I,5,N,M)*PNES/PZ(2)*PTIS
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE TRANSPORT COEFFICIENTS ******
C
      SUBROUTINE WFCOEF(K,IE,AMK,ADK,AHK,ANK,APK,AVT)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION AMK(2,2,NSM),ADK(2,2,NSM),AHK(2,2,NSM)
      DIMENSION ANK(MFM,MFM),APK(MFM),AVT(2,NSM)
      DIMENSION AL(3)
      DIMENSION RN(NSM),RT(NSM),RNUE(NSM),RNUI(NSM),RNUN(NSM)
      DIMENSION RNU(NSM)
      REAL*8    A(3),B(3),C(3),EL(2)
C
      IN=IELM(K,IE)
C
C     *** SET LOCAL VALUES ***
C
      CALL WFBMAG(XD(IN),YD(IN),BABS,AL)
C
      RN(1)=PNE(IN)
      RT(1)=PTE(IN)
      VTE=SQRT(RT(1)*AEE/AME)
      WCE=-AEE*BABS/AME
C
      RN(2)=PNI(IN)
      RT(2)=PTI(IN)
      VTI=SQRT(RT(2)*AEE/(PA(2)*AMP))
      WCI=PZ(2)*AEE*BABS/(PA(2)*AMP)
C
C     *** CALCULATE STATIC ELECTRIC FIELD ***
C
      CALL WFABC(IE,A,B,C,S)
      EL(1)=0.D0
      EL(2)=0.D0
      DO I=1,3
         INL=IELM(I,IE)
         EL(1)=EL(1)-B(I)*PHI(INL)
         EL(2)=EL(2)-C(I)*PHI(INL)
      ENDDO
      ELPARA=EL(1)*AL(1)+EL(2)*AL(2)
      ELPERP=SQRT(EL(1)**2+EL(2)**2-ELPARA**2)
C
C     *** CALCULATE NEUTRAL DENSITY ***
C
      PNN0=PPN0/(PTN0*AEE)
C
      CALL CALRNU(RN,RT,RNUE,RNUI,RNUN)
C
      RNU(1)=RNUE(1)+RNUI(1)+RNUN(1)
      RNU(2)=RNUE(2)+RNUI(2)+RNUN(2)
C
      CWE=DCMPLX(2.D0*PI*RF*1.D6,RNU(1))
C
C     *** VELOCITY PERTURBATION BY WAVE FIELD ***
C
      CWC=WCE/CWE
      CZ=-CI*AEE/(AME*CWE)
      CX=CZ/(1.D0-CWC*CWC)
      CY=CI*CWC*CX
      CEX=(1.D0-AL(1)*AL(1))*CEF(1,IN)
     &         -AL(1)*AL(2) *CEF(2,IN)
     &         -AL(1)*AL(3) *CEF(3,IN)
      CEY=AL(3)*CEF(2,IN)-AL(2)*CEF(3,IN)
      CEZ=AL(1)*CEF(1,IN)+AL(2)*CEF(2,IN)+AL(3)*CEF(3,IN)
      CVX= CX*CEX+CY*CEY
      CVY=-CY*CEX+CX*CEY
      CVZ= CZ*CEZ
      VE2=ABS(CVX)**2+ABS(CVY)**2+ABS(CVZ)**2
C
C     *** VELOCITY PERTURBATION BY ELECTROSTATIC FIELD ***
C
      UPARAE=-AEE/(AME*RNU(1))
      UPERPE=UPARAE/(1+WCE*WCE/(RNU(1)*RNU(1)))
      VEL2=(UPARAE*ELPARA)**2+(UPERPE*ELPERP)**2
C
C     *** EFFECTIVE ELECTRON TEMPARATURE ***
C
      TEFF=RT(1)+0.5D0*AME*VE2/AEE+0.5D0*AME*VEL2/AEE
C
C     *** IONIZATION ***
C
C        UI     : Ionization Potential [eV]
C        RNUION : Ionization Frequency
C
C      UI=13.6D0   (Hydrogen)
C      UI=15.76D0  (Argon)
C
      UI=13.6D0
      U=UI/TEFF
      IF(U.GT.100.D0) THEN
         RNUION=0.D0
      ELSE
         RNUION=PNN0*1.D-11*SQRT(1.D0/U)*EXP(-U)
     &          /(SQRT(UI)**3*(6.D0+1.D0/U))
      ENDIF
C
C     *** TRANSPORT COEF IN PARA-PERP COORDINATES ***
C
      UPARAE=-AEE/(AME*RNU(1))
      UPERPE=UPARAE/(1+WCE*WCE/(RNU(1)*RNU(1)))
C
      DPARAE=TEFF/RT(1)*VTE*VTE/RNU(1)
      DPERPE=TEFF/RT(1)*VTE*VTE/RNU(1)/(1+WCE*WCE/(RNU(1)*RNU(1)))
C
      HPARAE=3.16D0*TEFF/RT(1)*VTE*VTE/RNU(1)
      HPERPE=4.66D0*TEFF/RT(1)*VTE*VTE/RNU(1)
     &       /(1+WCE*WCE/(RNU(1)*RNU(1)))
C
      UPARAI=PZ(2)*AEE/(PA(2)*AMP*RNU(2))
      UPERPI=UPARAI/(1+WCI*WCI/(RNU(2)*RNU(2)))
C
      DPARAI=VTI*VTI/RNU(2)
      DPERPI=VTI*VTI/RNU(2)/(1+WCI*WCI/(RNU(2)*RNU(2)))
      HPARAI=3.16D0*VTI*VTI/RNU(2)
      HPERPI=4.66D0*VTI*VTI/RNU(2)/(1+WCI*WCI/(RNU(2)*RNU(2)))
C
      VPARAE=VTE
      VPERPE=VTE/(1+WCE*WCE/(RNU(1)*RNU(1)))
      VPARAI=VTI
      VPERPI=VTI/(1+WCI*WCI/(RNU(2)*RNU(2)))
C
C     *** BOHM DIFFUSION ***
C
      IF(BABS.GT.1.D-8) THEN
         DCE=DC*RT(1)/(16.D0*BABS)
         DPERPE=DPERPE+DCE
         HPERPE=HPERPE+DCE
         DCI=DC*RT(2)/(16.D0*BABS)
         DPERPI=DPERPI+DCI
         HPERPI=HPERPI+DCI
      ENDIF
C
C     *** SAVE FOR GRAPHICS ***
C
      RNUG(IN,1)=RNU(1)
      RNUG(IN,2)=RNU(2)
      RZCLG(IN,1)=RNU(1)/(2.D0*PI*RF*1.D6)
      RZCLG(IN,2)=RNU(2)/(2.D0*PI*RF*1.D6)
      UPARAG(IN,1)=UPARAE
      UPERPG(IN,1)=UPERPE
      UPARAG(IN,2)=UPARAI
      UPERPG(IN,2)=UPERPI
      DPARAG(IN,1)=DPARAE
      DPERPG(IN,1)=DPERPE
      DPARAG(IN,2)=DPARAI
      DPERPG(IN,2)=DPERPI
      HPARAG(IN,1)=HPARAE
      HPERPG(IN,1)=HPERPE
      HPARAG(IN,2)=HPARAI
      HPERPG(IN,2)=HPERPI
      RIONG(IN)=RNUION
C      IF(XD(IN).GT.0.1899D0) THEN
C         WRITE(6,*) 'XX HPARA,IN,K,IE=',IN,K,IE,XD(IN),HPARAG(IN,1)
C      ENDIF
C
C     *** TRANSPORT COEF IN REAL COORDINATE ***
C
      URRE=AL(2)*AL(2)* UPERPE+AL(1)*AL(1)*UPARAE
      URZE=AL(1)*AL(2)*(UPARAE-UPERPE)
      UZRE=AL(1)*AL(2)*(UPARAE-UPERPE)
      UZZE=AL(2)*AL(2)* UPARAE+AL(1)*AL(1)*UPERPE
      DRRE=AL(2)*AL(2)* DPERPE+AL(1)*AL(1)*DPARAE
      DRZE=AL(1)*AL(2)*(DPARAE-DPERPE)
      DZRE=AL(1)*AL(2)*(DPARAE-DPERPE)
      DZZE=AL(2)*AL(2)* DPARAE+AL(1)*AL(1)*DPERPE
      HRRE=AL(2)*AL(2)* HPERPE+AL(1)*AL(1)*HPARAE
      HRZE=AL(1)*AL(2)*(HPARAE-HPERPE)
      HZRE=AL(1)*AL(2)*(HPARAE-HPERPE)
      HZZE=AL(2)*AL(2)* HPARAE+AL(1)*AL(1)*HPERPE
C
      URRI=AL(2)*AL(2)* UPERPI+AL(1)*AL(1)*UPARAI
      URZI=AL(1)*AL(2)*(UPARAI-UPERPI)
      UZRI=AL(1)*AL(2)*(UPARAI-UPERPI)
      UZZI=AL(2)*AL(2)* UPARAI+AL(1)*AL(1)*UPERPI
      DRRI=AL(2)*AL(2)* DPERPI+AL(1)*AL(1)*DPARAI
      DRZI=AL(1)*AL(2)*(DPARAI-DPERPI)
      DZRI=AL(1)*AL(2)*(DPARAI-DPERPI)
      DZZI=AL(2)*AL(2)* DPARAI+AL(1)*AL(1)*DPERPI
      HRRI=AL(2)*AL(2)* HPERPI+AL(1)*AL(1)*HPARAI
      HRZI=AL(1)*AL(2)*(HPARAI-HPERPI)
      HZRI=AL(1)*AL(2)*(HPARAI-HPERPI)
      HZZI=AL(2)*AL(2)* HPARAI+AL(1)*AL(1)*HPERPI
C
C     *** SET TRANSPORT COEF MATRIX ***
C
      AMK(1,1,1)=URRE*RN(1)
      AMK(1,2,1)=URZE*RN(1)
      AMK(2,1,1)=UZRE*RN(1)
      AMK(2,2,1)=UZZE*RN(1)
C
      AMK(1,1,2)=URRI*RN(2)
      AMK(1,2,2)=URZI*RN(2)
      AMK(2,1,2)=UZRI*RN(2)
      AMK(2,2,2)=UZZI*RN(2)
C
      ADK(1,1,1)=DRRE
      ADK(1,2,1)=DRZE
      ADK(2,1,1)=DZRE
      ADK(2,2,1)=DZZE
C
      ADK(1,1,2)=DRRI
      ADK(1,2,2)=DRZI
      ADK(2,1,2)=DZRI
      ADK(2,2,2)=DZZI
C
      AHK(1,1,1)=HRRE
      AHK(1,2,1)=HRZE
      AHK(2,1,1)=HZRE
      AHK(2,2,1)=HZZE
C
      AHK(1,1,2)=HRRI
      AHK(1,2,2)=HRZI
      AHK(2,1,2)=HZRI
      AHK(2,2,2)=HZZI
C
      VTRE=SQRT(8.D0/PI)*( AL(1)*VPARAE+AL(2)*VPERPE)
      VTZE=SQRT(8.D0/PI)*(-AL(1)*VPERPE+AL(2)*VPARAE)
      VTRI=SQRT(8.D0/PI)*( AL(1)*VPARAI+AL(2)*VPERPI)
      VTZI=SQRT(8.D0/PI)*(-AL(1)*VPERPI+AL(2)*VPARAI)
C
C     *** secondary electron generation coefficient ***
C
      DGAMSE=0.05D0
      AVT(1,1)=0.25D0*VTRE*(1.D0-DGAMSE)
      AVT(2,1)=0.25D0*VTZE
C
      AVT(1,2)=0.25D0*VTRI
      AVT(2,2)=0.25D0*VTZI
C
C     *** COUPLING TERM ***
C
      DO 1000 I=1,5
      DO 1000 J=1,5
         ANK(I,J)=0.D0
 1000 CONTINUE
C
      ANK(2,2) = RNUION
      ANK(3,2) =-RNUION*UI
      ANK(3,3) =-RNUN(1)-RNUE(2)
      ANK(3,5) =         RNUE(2)*RN(1)/RN(2)
      ANK(4,2) = RNUION/PZ(2)
      ANK(5,3) =         RNUE(2)
      ANK(5,5) =-RNUN(2)-RNUE(2)*RN(1)/RN(2)
C      ANK(5,3) =0.D0
C      ANK(5,5) =0.D0
C
C     *** ABSORBED POWER ***
C
C        *** GIVEN PROFILE ***
C
         X=XD(IN)
         Y=YD(IN)
C
         XX=(X-XGIVEN)**2/RGIVEN**2
         YY=(Y-YGIVEN)**2/RGIVEN**2
         IF(XX.LT.100.D0.AND.YY.LT.100.D0) THEN
            SIONZ=SGIVEN* EXP(-XX-YY)
            PABSZ=PGIVEN* EXP(-XX-YY)
         ELSE
            SIONZ=0.D0
            PABSZ=0.D0
         END IF
C
C        *** JOULE HEATING ***
C
         PABSEL=RN(1)*1.D20*AEE
     &         *(ELPARA*UPARAE*ELPARA
     &          +ELPERP*UPERPE*ELPERP)
C
C        *** TOTAL POWER ***
C
         PABSE=PWRNS(IN,1)
     &        +            PABSZ      *(1.D20*AEE)
     &        +PABSEL
         PABSI=PWRNS(IN,2)
     &        +1.5D0*PTN0*SIONZ/PZ(2)*(1.D20*AEE)
C
C        *** SOURCE TERM ***
C
         APK(1) = 0.D0
         APK(2) = SIONZ
         APK(3) = PABSE/(1.D20*AEE) + 1.5D0*RNUN(1)*RN(1)*PTN0
         APK(4) = SIONZ/PZ(2)
         APK(5) = PABSI/(1.D20*AEE) + 1.5D0*RNUN(2)*RN(2)*PTN0
C         APK(5) = 0.D0
C
      RETURN
      END
C
C     ****** CONVERT FROM DSV TO FIELD VARIABLES ******
C
      SUBROUTINE DALFLD
C
      INCLUDE 'wfcomm.inc'
C
      PNEC=PNE0*1.D-4
C
      DO 10 IN=1,NNOD
         IF(KNODP(IN).EQ.0) THEN
            PHI(IN)=DSV(IBND(IN))
            PNE(IN)=DSV(IBND(IN)+1)
            PPE(IN)=DSV(IBND(IN)+2)
            PNI(IN)=DSV(IBND(IN)+3)
            PPI(IN)=DSV(IBND(IN)+4)
            PSQ(IN)=0.D0
            IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
               PNE(IN)=PNEC
               PNI(IN)=PNEC/PZ(2)
            ENDIF
            PTE(IN)=PPE(IN)/PNE(IN)
            IF(PTE(IN).LT.PTES) PTE(IN)=PTES
            PTI(IN)=PPI(IN)/PNI(IN)
            IF(PTI(IN).LT.PTIS) PTI(IN)=PTIS
            PPE(IN)=PNE(IN)*PTE(IN)
            PPI(IN)=PNI(IN)*PTI(IN)
         ELSE IF(KNODP(IN).GE.1.AND.KNODP(IN).LE.5) THEN
            PHIL=WFBPHI(XD(IN),YD(IN),KNODP(IN))
            PHI(IN)=PHIL
            IF(MODELW.EQ.0) THEN
               PNE(IN)=PNES
               PNI(IN)=PNES/PZ(2)
               PPE(IN)=PTES*PNE(IN)
               PPI(IN)=PTIS*PNI(IN)
               PTE(IN)=PTES
               PTI(IN)=PTIS
            ELSEIF(MODELW.EQ.1) THEN
               PNE(IN)=DSV(IBND(IN))
               PNI(IN)=DSV(IBND(IN)+1)
               IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
                  PNE(IN)=PNEC
                  PNI(IN)=PNEC/PZ(2)
               ENDIF
               PPE(IN)=PTES*PNE(IN)
               PPI(IN)=PTIS*PNI(IN)
               PTE(IN)=PTES
               PTI(IN)=PTIS
            ELSE
               PNE(IN)=DSV(IBND(IN))
               PPE(IN)=DSV(IBND(IN)+1)
               PNI(IN)=DSV(IBND(IN)+2)
               PPI(IN)=DSV(IBND(IN)+3)
               IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
                  PNE(IN)=PNEC
                  PNI(IN)=PNEC/PZ(2)
               ENDIF
               PTE(IN)=PPE(IN)/PNE(IN)
               IF(PTE(IN).LT.PTES) PTE(IN)=PTES
               PTI(IN)=PPI(IN)/PNI(IN)
               IF(PTI(IN).LT.PTIS) PTI(IN)=PTIS
               PPE(IN)=PNE(IN)*PTE(IN)
               PPI(IN)=PNI(IN)*PTI(IN)
            ENDIF
            PSQ(IN)=0.D0
         ELSE IF(KNODP(IN).EQ.8) THEN
            PHI(IN)=DSV(IBND(IN))
            IF(MODELW.EQ.0) THEN
               PNE(IN)=PNES
               PNI(IN)=PNES/PZ(2)
               PPE(IN)=PTES*PNE(IN)
               PPI(IN)=PTIS*PNI(IN)
               PTE(IN)=PTES
               PTI(IN)=PTIS
            ELSEIF(MODELW.EQ.1) THEN
               PNE(IN)=DSV(IBND(IN)+1)
               PNI(IN)=DSV(IBND(IN)+2)
               IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
                  PNE(IN)=PNEC
                  PNI(IN)=PNEC/PZ(2)
               ENDIF
               PPE(IN)=PTES*PNE(IN)
               PPI(IN)=PTIS*PNI(IN)
               PTE(IN)=PTES
               PTI(IN)=PTIS
            ELSE
               PNE(IN)=DSV(IBND(IN)+1)
               PPE(IN)=DSV(IBND(IN)+2)
               PNI(IN)=DSV(IBND(IN)+3)
               PPI(IN)=DSV(IBND(IN)+4)
               IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
                  PNE(IN)=PNEC
                  PNI(IN)=PNEC/PZ(2)
               ENDIF
               PTE(IN)=PPE(IN)/PNE(IN)
               IF(PTE(IN).LT.PTES) PTE(IN)=PTES
               PTI(IN)=PPI(IN)/PNI(IN)
               IF(PTI(IN).LT.PTIS) PTI(IN)=PTIS
               PPE(IN)=PNE(IN)*PTE(IN)
               PPI(IN)=PNI(IN)*PTI(IN)
            ENDIF
            PSQ(IN)=0.D0
         ELSE IF(KNODP(IN).EQ.9) THEN
            PHI(IN)=DSV(IBND(IN))
            IF(MODELW.EQ.0) THEN
               PNE(IN)=PNES
               PNI(IN)=PNES/PZ(2)
               PPE(IN)=PTES*PNE(IN)
               PPI(IN)=PTIS*PNI(IN)
               PSQ(IN)=DSV(IBND(IN)+1)
               PTE(IN)=PTES
               PTI(IN)=PTIS
            ELSEIF(MODELW.EQ.1) THEN
               PNE(IN)=DSV(IBND(IN)+1)
               PNI(IN)=DSV(IBND(IN)+2)
               IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
                  PNE(IN)=PNEC
                  PNI(IN)=PNEC/PZ(2)
               ENDIF
               PPE(IN)=PTES*PNE(IN)
               PPI(IN)=PTIS*PNI(IN)
               PSQ(IN)=DSV(IBND(IN)+3)
               PTE(IN)=PTES
               PTI(IN)=PTIS
            ELSE
               PNE(IN)=DSV(IBND(IN)+1)
               PPE(IN)=DSV(IBND(IN)+2)
               PNI(IN)=DSV(IBND(IN)+3)
               PPI(IN)=DSV(IBND(IN)+4)
               IF(PNE(IN).LT.PNEC.OR.PNI(IN).LT.PNEC/PZ(2)) THEN
                  PNE(IN)=PNEC
                  PNI(IN)=PNEC/PZ(2)
               ENDIF
               PSQ(IN)=DSV(IBND(IN)+5)
               PTE(IN)=PPE(IN)/PNE(IN)
               IF(PTE(IN).LT.PTES) PTE(IN)=PTES
               PTI(IN)=PPI(IN)/PNI(IN)
               IF(PTI(IN).LT.PTIS) PTI(IN)=PTIS
               PPE(IN)=PNE(IN)*PTE(IN)
               PPI(IN)=PNI(IN)*PTI(IN)
            ENDIF
         ELSE IF(KNODP(IN).EQ.10) THEN
            PHI(IN)=DSV(IBND(IN))
            PNE(IN)=0.D0
            PPE(IN)=0.D0
            PNI(IN)=0.D0
            PPI(IN)=0.D0
            PSQ(IN)=0.D0
            PTE(IN)=0.D0
            PTI(IN)=0.D0
         ELSE IF(KNODP(IN).GE.11.AND.KNODP(IN).LE.15) THEN
            PHIL=WFBPHI(XD(IN),YD(IN),KNODP(IN)-10)
            PHI(IN)=PHIL
            PNE(IN)=0.D0
            PPE(IN)=0.D0
            PNI(IN)=0.D0
            PPI(IN)=0.D0
            PSQ(IN)=0.D0
            PTE(IN)=0.D0
            PTI(IN)=0.D0
         ELSE
            WRITE(6,*) 'XX UNKNOWN KNODP: IN,KNODP =',IN,KNODP(IN)
         ENDIF
C
   10 CONTINUE
      RETURN
      END
C
C     ******* FRONTAL ELIMINATION WITH LEAST DIAGONAL PIVOTING *******
C
      SUBROUTINE DVSOLV(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      DATA ND1/23/
C
C     OPEN(ND1,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=MTRACK)
      IBUFF=0
C
      DO 5 I=1,MLEN
         DRV(I)=0.D0
         DSV(I)=0.D0
    5 CONTINUE
C
C     FIND LAST APPEAREANCE OF EACH NODE
C
      DO 10 I=1,NNOD
         NFLG(I)=0
   10 CONTINUE
      DO 20 N=NELM,1,-1
      DO 20 L=1,3
         M=ABS(JELM(L,N))
         IF(NFLG(M).EQ.0) THEN
            NFLG(M)=1
            JELM(L,N)=-JELM(L,N)
         ENDIF
   20 CONTINUE
C
C     PREFRONT
C
      NMAX=MBNDM
      NELL=0
      LCOL=0
      ILEN=0
      IPOS=0
      DO 30 I=1,NMAX
      DO 30 J=1,NMAX
         DEQ(J,I)=0.D0
   30 CONTINUE
C
 1000 NELL=NELL+1
C
      CALL DMCALC(NELL)
C
      KC=0
      DO 40 J=1,3
         NN=JELM(J,NELL)
         M=ABS(NN)
         K=IBND(M)
      DO 40 L=1,JBND(M)
         KC=KC+1
         II=K+L-1
         IF(NN.LT.0) II=-II
         NK(KC)=II
   40 CONTINUE
C
C     SET UP HEADING VECTORS
C
      DO 60 LK=1,KC
         NODE1=NK(LK)
         DO 50 L=1,LCOL
            IF(ABS(NODE1).EQ.ABS(LHED(L))) THEN
               LDEST(LK)=L
               LHED(LDEST(LK))=NODE1
               GOTO 60
            ENDIF
   50    CONTINUE
         LCOL=LCOL+1
         IF(LCOL.GT.NMAX) GOTO 9000
         LDEST(LK)=LCOL
         LHED(LCOL)=NODE1
   60 CONTINUE
C
C     ASSEMBLY
C
      L=0
      DO 70 J=1,3
      DO 70 JJ=1,JBND(ABS(JELM(J,NELL)))
         L=L+1
         LL=LDEST(L)
         K=0
      DO 70 I=1,3
      DO 70 II=1,JBND(ABS(JELM(I,NELL)))
         K=K+1
         KK=LDEST(K)
         DEQ(KK,LL)=DEQ(KK,LL)+DM(II,JJ,I,J)
   70 CONTINUE
C
      K=0
      DO 75 I=1,3
      DO 75 II=1,JBND(ABS(JELM(I,NELL)))
         K=K+1
         KK=LDEST(K)
         LCO=ABS(LHED(KK))
         DRV(LCO)=DRV(LCO)+DV(II,I)
   75 CONTINUE
C
C     FIND OUT WHICH MATRIX ELEMENTS ARE FULLY SUMMED
C
 2000 LC=0
      DO 80 L=1,LCOL
         IF(LHED(L).LT.0) THEN
            LC=LC+1
            LPIV(LC)=L
         ENDIF
   80 CONTINUE
C
      IF(LC.LE.0) GOTO 1000
C
C     SEARCH FOR ABSOLUTE PIVOT
C
      DPIVOT=0.D0
      DO 90 L=1,LC
         LPIVC=LPIV(L)
         DPIVA=DEQ(LPIVC,LPIVC)
         IF(ABS(DPIVA).GE.ABS(DPIVOT)) THEN
            DPIVOT=DPIVA
            LPIVCO=LPIVC
         ENDIF
   90 CONTINUE
C
C     NORMALIZE PIVOTAL ROW
C
      LCO=ABS(LHED(LPIVCO))
      IF(ABS(DPIVOT).LT.1.D-15) THEN
         WRITE(6,601) NELL,DPIVOT
         GOTO 9000
      ENDIF
      DO 100 L=1,LCOL
         DQQ(L)=DEQ(LPIVCO,L)/DPIVOT
  100 CONTINUE
      DRHS=DRV(LCO)/DPIVOT
      DRV(LCO)=DRHS
C
C     ELIMINATE THEN DELETE PIVOTAL ROW AND COLUMN
C
      DO 120 K=1,LPIVCO-1
         DFAC=DEQ(K,LPIVCO)
         KRW=ABS(LHED(K))
         DRV(KRW)=DRV(KRW)-DFAC*DRHS
         DO 110 L=1,LPIVCO-1
            DEQ(K,L)=DEQ(K,L)-DFAC*DQQ(L)
  110    CONTINUE
         DO 120 L=LPIVCO+1,LCOL
            DEQ(K,L-1)=DEQ(K,L)-DFAC*DQQ(L)
  120 CONTINUE
      DO 140 K=LPIVCO+1,LCOL
         DFAC=DEQ(K,LPIVCO)
         KRW=ABS(LHED(K))
         DRV(KRW)=DRV(KRW)-DFAC*DRHS
         DO 130 L=1,LPIVCO-1
            DEQ(K-1,L)=DEQ(K,L)-DFAC*DQQ(L)
  130    CONTINUE
         DO 140 L=LPIVCO+1,LCOL
            DEQ(K-1,L-1)=DEQ(K,L)-DFAC*DQQ(L)
  140 CONTINUE
      DO 150 L=LPIVCO+1,LCOL
         LHED(L-1)=LHED(L)
         DQQ (L-1)=DQQ (L)
  150 CONTINUE
      DO 160 K=1,LCOL
         DEQ(K,LCOL)=0.D0
  160 CONTINUE
      DO 170 L=1,LCOL-1
         DEQ(LCOL,L)=0.D0
  170 CONTINUE
      LCOL=LCOL-1
C
C     WRITE PIVOTAL EQUATION ON DISC
C
      IF(IPOS+LCOL.GT.MBUFM) THEN
         IBUFF=IBUFF+1
         WRITE(ND1,REC=IBUFF) MHED,DBUF
         IPOS=0
      ENDIF
      ILEN=ILEN+1
      MLCO(ILEN)=LCO
      MCOL(ILEN)=LCOL
      MPOS(ILEN)=IPOS
      DO 180 L=1,LCOL
         MHED(IPOS+L)=ABS(LHED(L))
         DBUF(IPOS+L)=DQQ(L)
  180 CONTINUE
      IPOS=IPOS+LCOL
C
C     DETERMINE WHETHER TO ASSEMBLE OR BACKSUBSTITUTE
C
      IF(LCOL.GT.1) GOTO 2000
      LCO=ABS(LHED(1))
      DPIVOT=DEQ(1,1)
      IF(ABS(DPIVOT).LT.1.D-15) THEN
         WRITE(6,601) NELL,DPIVOT
         GOTO 9000
      ENDIF
      DSV(LCO)=DRV(LCO)/DPIVOT
      IF(IBUFF.NE.0) WRITE(6,*) '## BUFFER MAX = ',IBUFF
C
C     BACK SUBSTITUTION
C
      DO 200 ILEN=MLEN-1,1,-1
         IF(IPOS.EQ.0) THEN
            READ(ND1,REC=IBUFF) MHED,DBUF
            IBUFF=IBUFF-1
         ENDIF
         LCO =MLCO(ILEN)
         IPOS=MPOS(ILEN)
         DGASH=0.D0
         DO 190 L=1,MCOL(ILEN)
            DGASH=DGASH-DBUF(IPOS+L)*DSV(MHED(IPOS+L))
  190    CONTINUE
         DSV(LCO)=DRV(LCO)+DGASH
  200 CONTINUE
      IERR=0
      RETURN
C
 9000 IERR=2
      WRITE(6,602) IERR
      RETURN
C
  601 FORMAT(1H ,'## FRONT WARNING : SINGULAR OR ILL CONDITIONED'/
     &       1H ,'       AT NELL = ',I5,'  CPIVOT=',1P2E12.4)
  602 FORMAT(1H ,'## FRONT ERROR : ERR =',I5/
     &       1H ,'       LCOL EXCEEDS NMAX')
      END
