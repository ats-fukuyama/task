C     $Id$
C
C     ********** WF WAVE SOLVER ( A VECTOR, FIRST ORDER ) **********
C
      SUBROUTINE WFWAVE
C
      INCLUDE 'wfcomm.inc'
C
      GTMAIN=0.0
      GTSOLV=0.0
C
      CALL GUTIME(GCPUT0)
C
      WRITE(6,*) '--- WFWPRE started ---'
      CALL WFWPRE(IERR)
         IF(IERR.NE.0) GOTO 9000
C
      WRITE(6,*) '--- CVCALC started ---'
      CALL CVCALC
C
      CALL GUTIME(GCPUT1)
C
      WRITE(6,*) '--- CVSOLV started ---'
      CALL CVSOLV(IERR)
C
      CALL GUTIME(GCPUT2)
         IF(IERR.NE.0) GOTO 9000
C
      CALL CALFLD
      CALL PWRABS
      CALL PWRRAD
      CALL TERMEP
      CALL WFCALB
         CALL GUTIME(GCPUT3)
         GTSOLV=GTSOLV+GCPUT2-GCPUT1
         GTMAIN=GTMAIN+GCPUT3-GCPUT2+GCPUT1-GCPUT0
C
      CALL LPEFLD
         WRITE(6,100) GTMAIN,GTSOLV
  100 FORMAT(' ','****** CPU TIME : MAIN = ',F10.3,' SEC',5X,
     &           ': SOLV = ',F10.3,' SEC ******')
C
 9000 RETURN
      END
C
C     ********** WF WAVE PREPARATION **********
C
      SUBROUTINE WFWPRE(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      IERR=0
C
      CALL LPELMT
C
      WRITE(6,*) '----- SETBDY start ---'
      CALL SETBDY(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,*) '----- SETSID start ---'
      CALL SETSID(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,*) '----- MODANT start ---'
      CALL MODANT(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,*) '----- SETKAN start ---'
      CALL SETKAN(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,*) '----- DEFBND start ---'
      CALL DEFBND(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,*) '----- DEFMBN start ---'
      CALL DEFMBN(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,'(6X,6A10)')  
     &     '     NNMAX','     NEMAX','    NSDMAX',
     &     '    NSFMAX','      MLEN','      MBND'
      WRITE(6,'(A4,2X,6I10)') 
     &     'USED',NNMAX,NEMAX,NSDMAX,NSFMAX,MLEN ,MBND
      WRITE(6,'(A4,2X,6I10)') 
     &     'SIZE',NNM  ,NEM  ,NSDM  ,NSFM  ,MLENM,MBNDM
C
      RETURN
      END
C
C     ******* SET KANOD *******
C
      SUBROUTINE SETKAN(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      IERR=0
      IF(MODELN.EQ.0.OR.MODELN.EQ.4) THEN
         NBMAX=0
         NMMAX=0
         NKMAX=1
         NMKA(1)=0
         DO NE=1,NEMAX
            KAELM(NE)=1
         ENDDO
      ELSEIF(MODELN.EQ.1.OR.MODELN.EQ.5) THEN
         NBMAX=0
      ELSEIF(MODELN.EQ.2.OR.MODELN.EQ.6) THEN
         NMMAX=1
         EPSDM(1)=1.D0
         AMUDM(1)=1.D0
         SIGDM(1)=0.D0
         NKMAX=1
         NMKA(1)=1
         DO NE=1,NEMAX
            KAELM(NE)=1
         ENDDO
      ELSEIF(MODELN.GE.100) THEN
         CALL SETWGB(IERR)
         IF(IERR.NE.0) RETURN
         NMMAX=1
         EPSDM(1)=1.D0
         AMUDM(1)=1.D0
         SIGDM(1)=0.D0
         NKMAX=1
         NMKA(1)=1
         DO NE=1,NEMAX
            KAELM(NE)=1
         ENDDO
      ENDIF
C
      DO NB=1,NBMAX
         WRITE(6,'(3I8,4X,A)') NB,KABDY(NB),NBPMAX(NB),KDBDY(NB)
         IF(KABDY(NB).LT.8) THEN
            DO NBP=1,NBPMAX(NB)
               NE=NENBP(NBP,NB)
               DO IN=1,4
                  NN=NDELM(IN,NE)
                  IF(NN.NE.NDNBP(NBP,NB)) THEN
                     IF(KANOD(NN).EQ.0) THEN
                        WRITE(6,'(A)')
     &                       'XX INCONSISTENT KANOD=0 ON BOUNDARY:'
                        WRITE(6,'(A,6I8)')
     &                       '   NB,NBP,NE,IN,NN,KABDY=',
     &                           NB,NBP,NE,IN,NN,KABDY(NB)
                        WRITE(6,'(A,1P4E12.4)') 
     &                       '   X,Y,Z,R=',XND(NN),YND(NN),ZND(NN),
     &                                  SQRT(XND(NN)**2+YND(NN)**2)
                        IERR=1
                     ELSE
                        IF(KABDY(NB).EQ.4) THEN
C                           WRITE(6,'(A,5I8)') 'NB,NBP,NE,NN,KA=',
C     &                          NB,NBP,NE,NN,KANOD(NN)
                           IF(MODELN.EQ.4.OR.MODELN.EQ.5.OR.
     &                        MODELN.EQ.6.OR.MODELN.EQ.7.OR.
     &                        MODELN.GE.100) THEN
                              KANOD(NN)=-NB
                           ELSE
                              KANOD(NN)=1
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO NBP=1,NBPMAX(NB)
               NN=NDNBP(NBP,NB)
               IF(KANOD(NN).EQ.0) THEN
                  WRITE(6,'(A)')
     &                 'XX INCONSISTENT KANOD=0 ON BOUNDARY:'
                  WRITE(6,'(A,6I8)')
     &                 '   NB,NBP,NE,IN,NN,KABDY=',
     &                     NB,NBP,NE,IN,NN,KABDY(NB)
                  WRITE(6,'(A,1P4E12.4)') 
     &                 '   X,Y,Z,R=',XND(NN),YND(NN),ZND(NN),
     &                            SQRT(XND(NN)**2+YND(NN)**2)
                  IERR=1
               ELSE
                  KANOD(NN)=-NB
               ENDIF
            ENDDO
         ENDIF
C
         DO NSF=1,NSFMAX
            KAN1=KANOD(NDSRF(1,NSF))
            KAN2=KANOD(NDSRF(2,NSF))
            KAN3=KANOD(NDSRF(3,NSF))
            IF(KAN1.EQ.-NB.AND.
     &         KAN2.EQ.-NB.AND.
     &         KAN3.EQ.-NB) THEN
               DO ISD=1,3
                  KASID(ABS(NSDSRF(ISD,NSF)))=-NB
               ENDDO
            ENDIF
         ENDDO
C
         DO NE=1,NEMAX
            DO IN=1,4
               IN1=MOD(IN,4)+1
               IN2=MOD(IN+1,4)+1
               IN3=MOD(IN+2,4)+1
               K1=KANOD(NDELM(IN1,NE))
               K2=KANOD(NDELM(IN2,NE))
               K3=KANOD(NDELM(IN3,NE))
               IF(K1.EQ.-NB.AND.
     &            K2.EQ.-NB.AND.
     &            K3.EQ.-NB) THEN
                  KNELM(IN,NE)=-NB
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** DIELECTRIC TENSOR ******
C
      SUBROUTINE DTENSR(NN,CK)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CK(3,3,NSM)
      DIMENSION WP(NSM),WC(NSM)
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
      DIMENSION AL(3)
C
      WW=2.D0*PI*RF*1.D6
C
      DO NS=1,NSMAX
         WP(NS)=PZ(NS)*PZ(NS)*AEE*AEE*1.D20/(PA(NS)*AMP*EPS0*WW*WW)
         WC(NS)=PZ(NS)*AEE/(PA(NS)*AMP*WW)
      ENDDO
C
      CALL WFSMAG(NN,BABS,AL)
      CALL WFSDEN(NN,RN,RTPR,RTPP,RZCL)
C
C         WRITE(6,*) 'NN,RN(1),RN(2)=',NN,RN(1),RN(2)
C
      DO NS=1,NSMAX
         CWP=WP(NS)*RN(NS)/(1.D0+CI*RZCL(NS))
         CWC=WC(NS)*BABS/(1.D0+CI*RZCL(NS))
         CDT0= CWP/(1.D0-CWC**2)
         CDX0= CI*CWP*CWC/(1.D0-CWC**2)
         CDP0= CWP
C
         CDT=CDT0
         CDP=CDP0-CDT
         CDX=CDX0
C
         FX=AL(1)
         FY=AL(2)
         FZ=AL(3)
         CXX= CDT   +CDP*FX*FX
         CXY= CDX*FZ+CDP*FX*FY
         CXZ=-CDX*FY+CDP*FX*FZ
         CYX=-CDX*FZ+CDP*FY*FX
         CYY= CDT   +CDP*FY*FY
         CYZ= CDX*FX+CDP*FY*FZ
         CZX= CDX*FY+CDP*FZ*FX
         CZY=-CDX*FX+CDP*FZ*FY
         CZZ= CDT   +CDP*FZ*FZ
C
         CK(1,1,NS)=-CXX
         CK(1,2,NS)=-CXY
         CK(1,3,NS)=-CXZ
         CK(2,1,NS)=-CYX
         CK(2,2,NS)=-CYY
         CK(2,3,NS)=-CYZ
         CK(3,1,NS)=-CZX
         CK(3,2,NS)=-CZY
         CK(3,3,NS)=-CZZ
      ENDDO
C
C         WRITE(6,'(A,I8,1P2E12.4)') 'NE,CK(3,3,1)=',NE,CK(3,3,1)
C
      RETURN
      END
C
C     ****** CURRENT COEFFICIENT VECTOR CALCULATION ******
C
      SUBROUTINE CVCALC
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CJ(3),W(4)
      DIMENSION F(4),DF(3,4)
      DIMENSION RWE(3,6),DWE(3,4,6)
C
      RW=2.D0*PI*RF*1.D6
C
      DO NE=1,NEMAX
         DO N=1,7
            CVTOT(N,NE)=0.D0
         ENDDO
      ENDDO
C
      DO NA=1,NAMAX
         PHASE =APH(NA)*PI/180.D0
         CVJ=CI*RW*AMU0*AJ(NA)*EXP(CI*PHASE)
         DO IJ=2,JNUM(NA)
            NE=JELMT(IJ,NA)
            CALL WFABCDX(NE,F,DF,RWE,DWE,V)
            X1=XJ(IJ-1,NA)
            Y1=YJ(IJ-1,NA)
            Z1=ZJ(IJ-1,NA)
            X2=XJ(IJ,NA)
            Y2=YJ(IJ,NA)
            Z2=ZJ(IJ,NA)
            CJ(1)=CVJ*(X2-X1)
            CJ(2)=CVJ*(Y2-Y1)
            CJ(3)=CVJ*(Z2-Z1)
            XM=0.5D0*(X1+X2)
            YM=0.5D0*(Y1+Y2)
            ZM=0.5D0*(Z1+Z2)
            DO IN=1,4
               W(IN)=F(IN)+DF(1,IN)*XM+DF(2,IN)*YM+DF(3,IN)*ZM
            ENDDO
            DO ISD=1,6
               DO IN=1,4
                  CVTOT(ISD,NE)=CVTOT(ISD,NE)
     &                    +(DWE(1,IN,ISD)*W(IN)*CJ(1)
     &                     +DWE(2,IN,ISD)*W(IN)*CJ(2)
     &                     +DWE(3,IN,ISD)*W(IN)*CJ(3))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** LOCAL ELEMENT MATRIX WITH WG******
C
      SUBROUTINE CMCALC(NE,CM,CV)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CM(NMDM,NMDM,7,7),CV(NMDM,7)
      DIMENSION CKL(3,3,4)
      DIMENSION CK(3,3,NSM)
      DIMENSION F(4),DF(3,4)
      DIMENSION RWE(3,6),DWE(3,4,6)
      DIMENSION CEWG(0:NMDM),CBWG(3,0:NMDM,3)
C
      RW=2.D0*PI*RF*1.D6
      WC=RW/VC
      WC2=WC**2
C
      NME=NMKA(KAELM(NE))
      NBE=NBELM(NE)
      IF(NBE.EQ.0) THEN
         ISDMAX=6
         IMDMAX=1
      ELSE
         ISDMAX=7
         IMDMAX=NMBDY(NBE)
      ENDIF
C
      CALL WFABCDX(NE,F,DF,RWE,DWE,V)
C
      DO IM=1,ISDMAX
         DO IN=1,ISDMAX
            DO J=1,IMDMAX
               DO I=1,IMDMAX
                  CM(I,J,IN,IM)=(0.D0,0.D0)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO IN=1,ISDMAX
         DO I=1,IMDMAX
            CV(I,IN)=0.D0
         ENDDO
         CV(1,IN)=CVTOT(IN,NE)
      ENDDO
C
      IF(NME.EQ.0) THEN
         DO K=1,4
            NK=NDELM(K,NE)
            CALL DTENSR(NK,CK)
            DO I=1,3
               DO J=1,3
                  IF(I.EQ.J) THEN
                     CKL(I,J,K)=1.D0
                  ELSE
                     CKL(I,J,K)=0.D0
                  ENDIF
                  DO NS=1,NSMAX
                     CKL(I,J,K)=CKL(I,J,K)+CK(I,J,NS)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         XMU=1.D0
      ELSE
         DO K=1,4
            IF(MODELA.EQ.0) THEN
               ID=0
            ELSE
               NK=NDELM(K,NE)
               IF(MODELA.EQ.1) THEN
                  POS=XND(NK)
               ELSEIF(MODELA.EQ.2) THEN
                  POS=YND(NK)
               ELSEIF(MODELA.EQ.3) THEN
                  POS=ZND(NK)
               ELSEIF(MODELA.EQ.4) THEN
                  POS=SQRT(XND(NK)**2+YND(NK)**2)
               ELSE
                  WRITE(6,*) 'XX CMCALC: UNDEFINED MODELA=',MODELA
               ENDIF
               IF((POSRES-POS)*(POS-POSABS).GE.0.D0) THEN
                  ID=1
               ELSE
                  ID=0
               ENDIF
            ENDIF
C
            IF(ID.EQ.0) THEN
               EPSD=EPSDM(NME)
               SIGD=SIGDM(NME)
               AMUD=AMUDM(NME)
            ELSE
               CEPS=EPSDM(NME)+2.D0*EPSABS*(POS-POSABS)**2
     &              /(ABS(POSRES-POSABS)*(ABS(POSRES-POS)-CI*DLTABS))
               EPSD=DREAL(CEPS)
               SIGD=DREAL(-CI*CEPS)*RW*EPS0
               AMUD=AMUDM(NME)
            ENDIF
C
            DO I=1,3
               DO J=1,3
                  IF(I.EQ.J) THEN
                     CKL(I,J,K)=EPSD+CI*SIGD/(RW*EPS0)
                  ELSE
                     CKL(I,J,K)=0.D0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         XMU=1.D0/AMUD
      ENDIF
C
      DO N=1,6
         DO M=1,6
            CP=RWE(1,N)*RWE(1,M)+RWE(2,N)*RWE(2,M)+RWE(3,N)*RWE(3,M)
C
            CM(1,1,N,M)=CM(1,1,N,M)+XMU*V*CP
C
            DO I=1,3
            DO J=1,3
               DO K2=1,4
                  SUM=0.D0
                  DO K1=1,4
                  DO K3=1,4
                     SUM=SUM+DWE(I,K1,N)*DWE(J,K3,M)*AIG3(K1,K2,K3)
                  ENDDO
                  ENDDO
                  CM(1,1,N,M)=CM(1,1,N,M)-WC2*V*SUM*CKL(I,J,K2)
               ENDDO
            ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C     ***** BOUNDARY SURFACE *****
C
      DO L=1,4
         KN=KNELM(L,NE)
         IF(KN.LT.0) THEN
            NB=-KN
            KA=KABDY(NB)
            CALL WF2ABC(L,NE,S)
C
            IF(KA.GE.8) THEN
               DO K0S=1,3
                  K0=MOD(L+K0S-1,4)+1
                  ND=NDELM(K0,NE)
                  CALL WFBFWG(ND,CBWG(1,0,K0S))
               ENDDO
C
               DO N=1,6
                  DO K1S=1,3
                     K1=MOD(L+K1S-1,4)+1
                     DO K2S=1,3
                        K2=MOD(L+K2S-1,4)+1
                        DO I=1,3
                           CV(1,N)=CV(1,N)
     &                          +CI*WC*XMU*S*DWE(I,K1,N)
     &                          *CBWG(I,0,K2S)*AIF2(K1S,K2S)
                           DO J=1,IMDMAX
                              CM(1,J,N,7)=CM(1,J,N,7)
     &                             -CI*WC*XMU*S*DWE(I,K1,N)
     &                             *CBWG(I,J,K2S)*AIF2(K1S,K2S)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
            ELSEIF(KA.EQ.4) THEN
C               WRITE(6,'(A,3I8)') 'NE,NB,KA=',NE,NB,KA
               DO N=1,6
                  DO M=1,6
                     DO K1S=1,3
                        K1=MOD(L+K1S-1,4)+1
                        DO K2S=1,3
                           K2=MOD(L+K2S-1,4)+1
                           DO I=1,3
                              CM(1,1,N,M)=CM(1,1,N,M)
     &                             -CI*WC*XMU*S*DWE(I,K1,N)*DWE(I,K2,M)
     &                                         *AIF2(K1S,K2S)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C
C     ----- boundary tangential electric field is given -----
C
      DO M=1,6
         NSD=ABS(NSDELM(M,NE))
         KA=KASID(NSD)
         IF(KA.LT.0) THEN
            NB=-KA
            IF(KABDY(NB).GE.8) THEN
               CALL WFEFWG(NSD,CEWG)
               DO N=1,6
                  DO J=1,IMDMAX
                     CM(1,J,N,7)=CM(1,J,N,7)
     &                    +CM(1,1,N,M)*CEWG(J)
                  ENDDO
                  CV(1,N)=CV(1,N)
     &                 -CM(1,1,N,M)*CEWG(0)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C
C     ----- boundary weighting function are coupled -----
C
      DO N=1,6
         NSD=ABS(NSDELM(N,NE))
         KA=KASID(NSD)
         IF(KA.LT.0) THEN
            NB=-KA
            IF(KABDY(NB).GE.8) THEN
               CALL WFEFWG(NSD,CEWG)
               DO M=1,6
                  DO I=1,IMDMAX
                     CM(I,1,7,M)=CM(I,1,7,M)
     &                    +DCMPLX(CEWG(I))*CM(1,1,N,M)
                  ENDDO
               ENDDO
               M=7
               DO I=1,IMDMAX
                  DO J=1,IMDMAX
                     CM(I,J,7,M)=CM(I,J,7,M)
     &                    +DCMPLX(CEWG(I))*CM(1,J,N,M)
                  ENDDO
               ENDDO
               DO I=1,IMDMAX
                  CV(I,7)=CV(I,7)
     &                 +DCMPLX(CEWG(I))*CV(1,N)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ******* ELECTRIC FIELD CALCULATION *******
C
      SUBROUTINE CALFLD
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION F(4),DF(3,4)
      DIMENSION RWE(3,6),DWE(3,4,6)
      DIMENSION CE(3),CEWG(0:NMDM)
C
      DO NB=1,NBMAX
         IF(KABDY(NB).GE.8) THEN
            NN=NDBDY(NB)
            DO L=1,NMBDY(NB)
               CRFL(L,NB)=CSV(IMLEN(NN)+L-1)
               REFL=ABS(CRFL(L,NB))**2
               WRITE(6,'(A,2I5,1P3E12.4)') 
     &              'L,NB,CRFL,REFL=',L,NB,CRFL(L,NB),REFL
            ENDDO
         ENDIF
      ENDDO
C
      DO NSD=1,NSDMAX
         KA=KASID(NSD)
         IF(KA.EQ.0) THEN
            CESD(NSD)=CSV(IMLEN(NSD))
         ELSEIF(KA.EQ.1) THEN
            CESD(NSD)=0.D0
         ELSEIF(KA.EQ.4) THEN
            CESD(NSD)=CSV(IMLEN(NSD))
         ELSEIF(KA.LT.0) THEN
            NB=-KA
            CALL WFEFWG(NSD,CEWG)
            CESD(NSD)=CEWG(0)
            DO L=1,NMBDY(NB)
               CESD(NSD)=CESD(NSD)+CEWG(L)*CRFL(L,NB)
            ENDDO
         ELSE
            WRITE(6,*) 'XX INVALID KASID: KASID(',NSD,')= ',KASID(NSD)
         ENDIF
      ENDDO
C
      DO I=1,3
      DO NN=1,NNMAX
         CEF(I,NN)=0.D0
      ENDDO
      ENDDO
C
      DO NE=1,NEMAX
         CALL WFABCDX(NE,F,DF,RWE,DWE,V)
         DO IN=1,4
            CE(1)=(0.D0,0.D0)
            CE(2)=(0.D0,0.D0)
            CE(3)=(0.D0,0.D0)
            DO ISD=1,6
               NSD=ABS(NSDELM(ISD,NE))
               CE(1)=CE(1)+DWE(1,IN,ISD)*CESD(NSD)
               CE(2)=CE(2)+DWE(2,IN,ISD)*CESD(NSD)
               CE(3)=CE(3)+DWE(3,IN,ISD)*CESD(NSD)
            ENDDO
            NN=NDELM(IN,NE)
            CEF(1,NN)=CEF(1,NN)+CE(1)*0.25D0*VELM(NE)/VNOD(NN)
            CEF(2,NN)=CEF(2,NN)+CE(2)*0.25D0*VELM(NE)/VNOD(NN)
            CEF(3,NN)=CEF(3,NN)+CE(3)*0.25D0*VELM(NE)/VNOD(NN)
         ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** POWER ABSORPTION ******
C
      SUBROUTINE PWRABS
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION PABSL(NSM)
      DIMENSION XE(4),YE(4),ZE(4),CK(3,3,NSM),CKL(3,3,4,NSM)
      DIMENSION F(4),DF(3,4)
      DIMENSION RWE(3,6),DWE(3,4,6)
C
      RW=2.D0*PI*RF*1.D6
      CNST=0.5D0*CI*RW*EPS0
C
      PABST=0.D0
      DO NS=1,NSMAX
         PABSS(NS)=0.D0
      ENDDO
      DO NK=1,NKMAX
         PABSK(NK)=0.D0
      ENDDO
C
      DO NN=1,NNMAX
         PABSTN(NN)=0.D0
         DO NS=1,NSMAX
            PABSSN(NN,NS)=0.D0
         ENDDO
         DO NK=1,NKMAX
            PABSKN(NN,NK)=0.D0
         ENDDO
      ENDDO
C
      DO NEDO=1,NEMAX
         NE=NEDO
         NME=NMKA(KAELM(NE))
         NKE=KAELM(NE)
C
         IF(NME.EQ.0) THEN
            NSMAXL=NSMAX
         ELSE
            NSMAXL=1
         ENDIF
         DO NS=1,NSMAXL
            PABSL(NS)=0.D0
         ENDDO
C
         CALL WFNODE(NE,XE,YE,ZE)
         CALL WFABCDX(NE,F,DF,RWE,DWE,V)
C
         IF(NME.EQ.0) THEN
            DO K=1,4
               NK=NDELM(K,NE)
               CALL DTENSR(NK,CK)
C
               DO I=1,3
                  DO J=1,3
                     DO NS=1,NSMAX
                        CKL(I,J,K,NS)=CK(I,J,NS)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO K=1,4
               IF(MODELA.EQ.0) THEN
                  ID=0
               ELSE
                  NK=NDELM(K,NE)
                  IF(MODELA.EQ.1) THEN
                     POS=XND(NK)
                  ELSEIF(MODELA.EQ.2) THEN
                     POS=YND(NK)
                  ELSEIF(MODELA.EQ.3) THEN
                     POS=ZND(NK)
                  ELSEIF(MODELA.EQ.4) THEN
                     POS=SQRT(XND(NK)**2+YND(NK)**2)
                  ELSE
                     WRITE(6,*) 'XX PWRABS: UNDEFINED MODELA=',MODELA
                  ENDIF
                  IF((POSRES-POS)*(POS-POSABS).GE.0.D0) THEN
                     ID=1
                  ELSE
                     ID=0
                  ENDIF
               ENDIF
C
               IF(ID.EQ.0) THEN
                  EPSD=EPSDM(NME)
                  SIGD=SIGDM(NME)
               ELSE
                  CEPS=EPSDM(NME)+2.D0*EPSABS*(POS-POSABS)**2
     &                 /(ABS(POSRES-POSABS)*(ABS(POSRES-POS)-CI*DLTABS))
                  EPSD=DREAL(CEPS)
                  SIGD=DREAL(-CI*CEPS)*RW*EPS0
               ENDIF
               DO I=1,3
                  DO J=1,3
                     IF(I.EQ.J) THEN
                        CKL(I,J,K,1)=EPSD+CI*SIGD/(RW*EPS0)
                     ELSE
                        CKL(I,J,K,1)=0.D0
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
C
         DO N=1,6
            NSD=ABS(NSDELM(N,NE))
         DO M=1,6
            MSD=ABS(NSDELM(M,NE))
         DO K2=1,4
            DO I=1,3
            DO J=1,3
C
               SUM=0.D0
               DO K1=1,4
               DO K3=1,4
                  SUM=SUM+DWE(I,K1,N)*DWE(J,K3,M)*AIG3(K1,K2,K3)
               ENDDO
               ENDDO
               DO NS=1,NSMAXL
                  PABSL(NS)=PABSL(NS)
     &                     -DBLE(CNST*DCONJG(CESD(NSD))*SUM
     &                       *CKL(I,J,K2,NS)*CESD(MSD))
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
C
         IF(NME.EQ.0) THEN
            SUM=0.D0
            DO NS=1,NSMAX
               SUM      =SUM      +PABSL(NS)*V
               PABSS(NS)=PABSS(NS)+PABSL(NS)*V
               DO N=1,4
                  NN=NDELM(N,NE)
                  PABSTN(NN)    =PABSTN(NN)    +PABSL(NS)*V/4.D0
                  PABSKN(NN,NKE)=PABSKN(NN,NKE)+PABSL(NS)*V/4.D0
                  PABSSN(NN,NS) =PABSSN(NN,NS) +PABSL(NS)*V/4.D0
               ENDDO
            ENDDO
            PABST     =PABST     +SUM
            PABSK(NKE)=PABSK(NKE)+SUM
         ELSE
            PABST     =PABST     +PABSL(1)*V
            PABSK(NKE)=PABSK(NKE)+PABSL(1)*V
            DO N=1,4
               NN=NDELM(N,NE)
               PABSTN(NN)    =PABSTN(NN)    +PABSL(1)*V/4.D0
               PABSKN(NN,NKE)=PABSKN(NN,NKE)+PABSL(1)*V/4.D0
            ENDDO
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ****** RADIATED POWER ******
C
      SUBROUTINE PWRRAD
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CJ(3),W(4)
      DIMENSION F(4),DF(3,4)
      DIMENSION RWE(3,6),DWE(3,4,6)
C
      DO NA=1,NAMAX
         PHASE =APH(NA)*PI/180.D0
         CVJ=AJ(NA)*EXP(CI*PHASE)
         CIMP(NA)=0.D0
         DO IJ=2,JNUM(NA)
            NE=JELMT(IJ,NA)
            CALL WFABCDX(NE,F,DF,RWE,DWE,V)
            X1=XJ(IJ-1,NA)
            Y1=YJ(IJ-1,NA)
            Z1=ZJ(IJ-1,NA)
            X2=XJ(IJ,NA)
            Y2=YJ(IJ,NA)
            Z2=ZJ(IJ,NA)
            CJ(1)=CVJ*(X2-X1)
            CJ(2)=CVJ*(Y2-Y1)
            CJ(3)=CVJ*(Z2-Z1)
            XM=0.5D0*(X1+X2)
            YM=0.5D0*(Y1+Y2)
            ZM=0.5D0*(Z1+Z2)
            DO IN=1,4
               W(IN)=F(IN)+DF(1,IN)*XM+DF(2,IN)*YM+DF(3,IN)*ZM
            ENDDO
            DO ISD=1,6
               NSD=ABS(NSDELM(ISD,NE))
               DO IN=1,4
                  CIMP(NA)=CIMP(NA)
     &                    -0.5D0*DCONJG(CESD(NSD))
     &                     *(DWE(1,IN,ISD)*W(IN)*CJ(1)
     &                      +DWE(2,IN,ISD)*W(IN)*CJ(2)
     &                      +DWE(3,IN,ISD)*W(IN)*CJ(3))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CTIMP=0.D0
      DO NA=1,NAMAX
         CTIMP=CTIMP+CIMP(NA)
      ENDDO
C
      RETURN
      END
C
C     ****** COMPLETE EFIELD AND POWER ******
C
      SUBROUTINE TERMEP
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION EABS(3)
C
      IF(PIN.GE.0.D0.OR.PABST.EQ.0.D0) THEN
         FACT=1.D0
      ELSE
         FACT=-PIN/ABS(PABST)
      ENDIF
      FACTR=SQRT(FACT)
C
      PABST=FACT*PABST
      DO NN=1,NNMAX
         PABSTN(NN)=FACT*PABSTN(NN)/VNOD(NN)
      ENDDO
      DO NS=1,NSMAX
         PABSS(NS)=FACT*PABSS(NS)
         DO NN=1,NNMAX
            PABSSN(NN,NS)=FACT*PABSSN(NN,NS)/VNOD(NN)
         ENDDO
      ENDDO
      DO NK=1,NKMAX
         PABSK(NK)=FACT*PABSK(NK)
         DO NN=1,NNMAX
            PABSKN(NN,NK)=FACT*PABSKN(NN,NK)/VNOD(NN)
         ENDDO
      ENDDO
C
      PNMAX=PABSTN(1)
      DO NN=2,NNMAX
         PNMAX=MAX(PNMAX,PABSTN(NN))
      ENDDO
C
      ETMAX=0.D0
      DO I=1,3
         EMAX(I)=0.D0
      ENDDO
C
      DO NN=1,NNMAX
         DO J=1,3
            CEF(J,NN)=FACTR*CEF(J,NN)
            EABS(J)=ABS(CEF(J,NN))
            EMAX(J)=MAX(EMAX(J),EABS(J))
         ENDDO
         ETOTL =SQRT(EABS(1)*EABS(1)+EABS(2)*EABS(2)+EABS(3)*EABS(3))
         ETMAX =MAX(ETMAX,ETOTL)
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE MAGNETIC FIELD AND POLARIZED FIELD ******
C
      SUBROUTINE WFCALB
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION F(4),DF(3,4)
      DIMENSION RWE(3,6),DWE(3,4,6)
      DIMENSION CB(3),AL(3)
C
      RW=2.D0*PI*RF*1.D6
      COEFB=-CI*VC/RW
C
      DO I=1,3
      DO NN=1,NNMAX
         CBF(I,NN)=0.D0
      ENDDO
      ENDDO
C
      DO NE=1,NEMAX
         CALL WFABCDX(NE,F,DF,RWE,DWE,V)
         CB(1)=(0.D0,0.D0)
         CB(2)=(0.D0,0.D0)
         CB(3)=(0.D0,0.D0)
         DO ISD=1,6
            NSD=ABS(NSDELM(ISD,NE))
            CB(1)=CB(1)+COEFB*RWE(1,ISD)*CESD(NSD)
            CB(2)=CB(2)+COEFB*RWE(2,ISD)*CESD(NSD)
            CB(3)=CB(3)+COEFB*RWE(3,ISD)*CESD(NSD)
         ENDDO
         DO IN=1,4
            NN=NDELM(IN,NE)
            CBF(1,NN)=CBF(1,NN)+CB(1)*0.25D0*VELM(NE)/VNOD(NN)
            CBF(2,NN)=CBF(2,NN)+CB(2)*0.25D0*VELM(NE)/VNOD(NN)
            CBF(3,NN)=CBF(3,NN)+CB(3)*0.25D0*VELM(NE)/VNOD(NN)
         ENDDO
      ENDDO
C
      DO NN=1,NNMAX
         CALL WFSMAG(NN,BABS,AL)
         SUM=SQRT(AL(1)*AL(1)+AL(3)*AL(3))
         IF(SUM.EQ.0.D0) THEN
            CE1=CEF(3,NN)
            CE2=CEF(1,NN)
            CB1=CBF(3,NN)
            CB2=CBF(1,NN)
         ELSE
            CE1=(AL(3)*CEF(1,NN)-AL(1)*CEF(3,NN))/SUM
            CB1=(AL(3)*CBF(1,NN)-AL(1)*CBF(3,NN))/SUM
            CE2=SUM*CEF(2,NN)
     &         -AL(2)*(AL(1)*CEF(1,NN)+AL(3)*CEF(3,NN))/SUM
            CB2=SUM*CBF(2,NN)
     &         -AL(2)*(AL(1)*CBF(1,NN)+AL(3)*CBF(3,NN))/SUM
         ENDIF
         CEP(1,NN)=CE1-CI*CE2
         CEP(2,NN)=CE1+CI*CE2
         CEP(3,NN)=AL(1)*CEF(1,NN)+AL(2)*CEF(2,NN)+AL(3)*CEF(3,NN)
         CBP(1,NN)=CB1-CI*CB2
         CBP(2,NN)=CB1+CI*CB2
         CBP(3,NN)=AL(1)*CBF(1,NN)+AL(2)*CBF(2,NN)+AL(3)*CBF(3,NN)
      ENDDO
C
      FACTOR=0.5D0/(VC*AMU0)
      DO NN=1,NNMAX
         PFV(NN,1)=FACTOR*DBLE((DCONJG(CEF(2,NN))*CBF(3,NN)
     &                         -DCONJG(CEF(3,NN))*CBF(2,NN)))
         PFV(NN,2)=FACTOR*DBLE((DCONJG(CEF(3,NN))*CBF(1,NN)
     &                         -DCONJG(CEF(1,NN))*CBF(3,NN)))
         PFV(NN,3)=FACTOR*DBLE((DCONJG(CEF(1,NN))*CBF(2,NN)
     &                         -DCONJG(CEF(2,NN))*CBF(1,NN)))
C         WRITE(6,'(A,I8,1P3E12.4)') 'NN,PFV=',
C     &        NN,PFV(NN,1),PFV(NN,2),PFV(NN,3)
      ENDDO
C
      RETURN
      END
C
C     ******* OUTPUT FIELD DATA *******
C
      SUBROUTINE LPEFLD
C
      INCLUDE 'wfcomm.inc'
C
      IF(NPRINT.LT.1) RETURN
C
      WRITE(6,110) (EMAX(I),I=1,3),ETMAX,PNMAX
  110 FORMAT(' ','EXMAX  =',1PE12.4
     &      ,3X ,'EYMAX  =',1PE12.4
     &      ,3X ,'EZMAX  =',1PE12.4/
     &       ' ','EMAX   =',1PE12.4
     &      ,3X ,'PNMAX  =',1PE12.4)
C
      WRITE(6,120) DBLE(CTIMP),PABST
  120 FORMAT(' ','RADIATED POWER =',1PE12.4/
     &       ' ','ABSORBED POWER =',1PE12.4)
C
      IF(ABS(PABST).GT.1.D-32) THEN
         DO NS=1,NSMAX
            WRITE(6,125) NS,PABSS(NS)
  125       FORMAT(' ','      PABS(NS=',I2,') =',1PE12.4)
         ENDDO
         DO NK=1,NKMAX
            WRITE(6,126) NK,PABSK(NK)
  126       FORMAT(' ','      PABS(NK=',I2,') =',1PE12.4)
         ENDDO
      ENDIF
C
      IF(NAMAX.GT.0) THEN
         WRITE(6,130)
  130    FORMAT(' ',' I JNUM', ' AJ(I)','  APH(I)','  AWD(I)',
     &                         ' APOS(I)','  XJ(I)','   YJ(I)',
     &          8X,'RADIATED POWER')
      ENDIF
      DO NA=1,NAMAX
         WRITE(6,140) NA,JNUM(NA),AJ(NA),APH(NA),AWD(NA),APOS(NA),
     &                            XJ(1,NA),YJ(1,NA),CIMP(NA)
  140    FORMAT(' ',I3,I3,0PF7.4,F7.2,4F7.4,2X,'(',1P2E12.4,')')
      ENDDO
C
      IF(NPRINT.LT.2) RETURN
C
      WRITE(6,150) (I,KANOD(I),PABSTN(I),(CEF(J,I),J=1,3),I=1,NNMAX)
  150 FORMAT('ABSORBED POWER DENSITY AND ELECTRIC FIELD'/
     &       'NOD ','KB',3X,'POWER',12X,'EX',19X,'EY',19X,'EZ'/
     &      (I4,I2,1PE10.2,3(1X,1P2E10.2)))
C
      RETURN
      END
C
C     ******* OUTPUT ELEMENT DATA *******
C
      SUBROUTINE LPELMT
C
      INCLUDE 'wfcomm.inc'
C
      IF(NPRINT.LT.3) RETURN
C
      WRITE(6,110) NNMAX
  110 FORMAT(/' ','NODE DATA     : #### NNMAX =',I5,' ####'/
     &       ' ',2('  NNMAX',' IMLEN',' KANOD',' INLEN',
     &       9X,'X',14X,'Y',9X))
      WRITE(6,115) (I,IMLEN(I),KANOD(I),INLEN(I),XND(I),YND(I),
     &              I=1,NNMAX)
  115 FORMAT((' ',2(4I6,2X,1P2E15.7,2X)))
C
      WRITE(6,120) NEMAX,(I,(NDELM(J,I),J=1,4),I=1,NEMAX)
  120 FORMAT(/' ','ELEMENT DATA  : #### NEMAX =',I5,' ####'/
     &      (' ',5(I6,'(',4I5,')',2X)))
C
      WRITE(6,125) NEMAX,(I,(ISDELM(J,I),J=1,7),I=1,NEMAX)
  125 FORMAT(/' ','NOP     DATA  : #### NEMAX =',I5,' ####'/
     &      (' ',(I8,'(',7I8,')',2X)))
C
      WRITE(6,130) NBMAX
  130 FORMAT(/' ','BOUNDARY DATA : #### NBMAX =',I5,' ####'/
     &       ' ',3('  NO.',' NDBDY',2X))
      WRITE(6,135) (I,NDBDY(I),I=1,NBMAX)
  135 FORMAT((' ',3(2I5,2X)))
C
      DO NA=1,NAMAX
         WRITE(6,140) NA,JNUM0(NA)
  140    FORMAT(/' ','ORIGINAL ANTENNA DATA : NA =',I5,' JNUM0 =',I5/
     &          ' ',3('  NO.',13X,' XJ0',11X,' YJ0',6X))
         WRITE(6,150) (I,XJ0(I,NA),YJ0(I,NA),I=1,JNUM0(NA))
  150    FORMAT((' ',3(I5,8X,1P2E15.7)))
C
         WRITE(6,154) NA,JNUM(NA)
  154    FORMAT(/' ','MODIFIED ANTENNA DATA : NA =',I5,' JNUM  =',I5/
     &          ' ',3('  NO.',' JELM',8X,' JX ',11X,' JY ',6X))
         WRITE(6,156) (I,JELMT(I,NA),XJ(I,NA),YJ(I,NA),I=1,JNUM(NA))
  156    FORMAT((' ',3(2I5,3X,1P2E15.7)))
      ENDDO
C
      RETURN
      END
