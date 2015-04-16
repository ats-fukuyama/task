C     $Id$
C
C     ********** WF WAVE SOLVER ( A VECTOR, FIRST ORDER ) **********
C
      SUBROUTINE WFWAVE
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFSETW
      CALL LPELMT
C
      GTMAIN=0.0
      GTSOLV=0.0
C
      RKZSAV=RKZ
      CALL SETANT
      CALL INITEP
C
      DO 50 NZDO=1,NZMAX
         NZ=NZDO
            CALL GUTIME(GCPUT0)
         RKZ=RKZF(NZ)
         DO 10 IA=1,NAMAX
            CAJ(IA)=CAJF(NZ,IA)
   10    CONTINUE
         CALL CVDBND(IERR)
            IF(IERR.NE.0) GOTO 9000
         CALL CVCALC
C
            CALL GUTIME(GCPUT1)
         CALL CVSOLV(IERR,NZ)
            CALL GUTIME(GCPUT2)
            IF(IERR.EQ.9002.OR.IERR.EQ.30000) GOTO 9000
            IF(IERR.NE.0) GOTO 9000
C
         CALL CALFLD
         CALL PWRABS(NZ)
         CALL PWRRAD
            CALL GUTIME(GCPUT3)
            GTSOLV=GTSOLV+GCPUT2-GCPUT1
            GTMAIN=GTMAIN+GCPUT3-GCPUT2+GCPUT1-GCPUT0
   50 CONTINUE
C
      RKZ=RKZSAV
      CALL TERMEP
      CALL LPEFLD
         WRITE(6,100) GTMAIN,GTSOLV
  100 FORMAT(1H ,'****** CPU TIME : MAIN = ',F10.3,' SEC',5X,
     &           ': SOLV = ',F10.3,' SEC ******')
C
 9000 RETURN
      END
C
C     ****** SETUP ******
C
      SUBROUTINE WFSETW
C
      USE libbes,ONLY: BESKN
      INCLUDE 'wfcomm.inc'
C
      IF(MODELS.EQ.3) THEN
         RGAMMA=ABS(H1*RRC)
         IF(RGAMMA.LT.1.D-5) THEN
            HA1=0.D0
         ELSE
            HA1=RGAMMA*BESKN(1,2.D0*RGAMMA)+BESKN(2,2.D0*RGAMMA)
         ENDIF
         RKAP=SQRT((1.D0+2.D0*HA1)/(1.D0-2.D0*HA1))
         WRITE(6,'(A,1P2E12.4)') 'HA1,RKAP=',HA1,RKAP
      ENDIF
C
      CALL MODANT(IERR)
C
      DO NB=1,NBDYW
         MB=MOD(NB,NBDYW)+1
         IN=IBDYW(NB)
         IM=IBDYW(MB)
         X1=XD(IN)
         Y1=YD(IN)
         X2=XD(IM)
         Y2=YD(IM)
         RX=SQRT((X2-X1)**2+(Y2-Y1)**2)
         DLLB(NB)=RX
         ANXB(NB)=(Y2-Y1)/RX
         ANYB(NB)=(X1-X2)/RX
      ENDDO
C
      DO IN=1,NNOD
         DLLN(IN)=0.D0
         ANXN(IN)=0.D0
         ANYN(IN)=0.D0
      ENDDO
C
      DO NB=1,NBDYW
         NBN=MOD(NB-2+NBDYW,NBDYW)+1
         NBP=MOD(NB,NBDYW)+1
         INN=IBDYW(NBN)
         IN =IBDYW(NB)
         INP=IBDYW(NBP)
         IF(KNODW(IN).GE.2) THEN
            ANX1=ANXB(NBN)
            ANY1=ANYB(NBN)
            ANX2=ANXB(NB)
            ANY2=ANYB(NB)
            IF(KNODW(INN).GE.2) THEN
               IF(KNODW(INP).GE.2) THEN
                  R=SQRT(2.D0*(1+ANX1*ANX2+ANY1*ANY2))
                  ANXN(IN)=(ANX1+ANX2)/R
                  ANYN(IN)=(ANY1+ANY2)/R
                  DLLN(IN)=0.5D0*(DLLB(NBN)+DLLB(NB))
               ELSE
                  ANXN(IN)=ANX1
                  ANYN(IN)=ANY1
                  DLLN(IN)=0.5D0*DLLB(NBN)
               ENDIF
            ELSE
               IF(KNODW(INP).GE.2) THEN
                  ANXN(IN)=ANX2
                  ANYN(IN)=ANY2
                  DLLN(IN)=0.5D0*DLLB(NB)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      RETURN
      END
C
C     ******* SETING ANTENNA CURRENT DENSITY *******
C
      SUBROUTINE SETANT
C
      USE libfft,ONLY: FFT2L
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CJ(NZM)
      DATA ALOG2/0.6931 47180 55995D0/
C
      IF(NZMAX.EQ.1) THEN
         RKZF(1)=RKZ
         DO 10 IA=1,NAMAX
            PHASE =APH(IA)*PI/180.D0
            CAJF(1,IA)=AJ(IA)*CDEXP(DCMPLX(0.D0,PHASE))
   10    CONTINUE
      ELSE
         LP=INT(LOG(DBLE(NZMAX))/ALOG2 + 0.001D0)
         NZMAX=2**LP
         DZ=2.D0*PI*RZ/NZMAX
         IFFT=1
         DKZ=1.D0/RZ
         RKZF(1)=0.001D0
         RKZF(NZMAX/2+1)=DKZ*DBLE(NZMAX/2)
         DO 20 NZ=1,NZMAX/2-1
            RKZF(NZ+1)=DKZ*NZ
            RKZF(NZMAX-NZ+1)=-DKZ*NZ
   20    CONTINUE
C
         DO 40 IA=1,NAMAX
            DO 30 NZ=1,NZMAX
               CJ(NZ)=(0.D0,0.D0)
   30       CONTINUE
            NZ    =INT(APOS(IA)*PI/(180.D0*DZ)) + 1
            PHASE =APH(IA)*PI/180.D0
            CJ(NZ)=AJ(IA)*CDEXP(DCMPLX(0.D0,PHASE))/DZ
            CALL FFT2L(CJ,CAJF(1,IA),RFFT,LIST,NZMAX,IFFT,1)
   40    CONTINUE
      ENDIF
C
      RETURN
      END
C
C     ******* SET KBND ARRAY *******
C
      SUBROUTINE CVDBND(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      DO 10 IE=1,NELM
      DO 10 IN=1,3
         JELM(IN,IE)=IELM(IN,IE)
   10 CONTINUE
C
C     BOUNDARY CONDITION
C
C        KNODW = 0 : NOT ON BOUNDARY
C                1 : ON AXIS
C                2 : PHI=0
C                3 : PHI=0 ON AXIS
C                4 : PHI=PHIW
C                5 : PHI=PHIW ON AXIS
C                6 : PHI=0..PHIW
C                7 : PHI=0..PHIW ON AXIS
C
      IL=1
      DO 20 IN=1,NNOD
         IF(KNODW(IN).EQ.0) THEN
            K=4
         ELSEIF(MOD(KNODW(IN),2).EQ.0) THEN
            K=1
         ELSEIF(MODELS.EQ.1) THEN
            IF(KNODW(IN).EQ.1) THEN
               IF(NPHI.EQ.0) THEN
                  K=2
               ELSEIF(ABS(NPHI).EQ.1) THEN
                  K=1
               ELSE
                  K=0
               ENDIF
            ELSE
               IF(NPHI.EQ.0) THEN
                  K=1
               ELSEIF(ABS(NPHI).EQ.1) THEN
                  K=0
               ELSE
                  K=0
               ENDIF
            ENDIF
         ELSE
            WRITE(6,*) 'XX ERROR IN DVDNOD: UNDEFINED KNODW'
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
 9000 WRITE(6,*) 'XX CVDBND ERROR : MLEN,MLENM =',MLEN,MLENM
      IERR=900
      RETURN
      END
C
C     ****** DIELECTRIC TENSOR ******
C
      SUBROUTINE DTENSR(IN,CK)
C
      USE libdsp,ONLY: dspfn
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CK(3,3,NSM)
      DIMENSION WP(NSM),WC(NSM),TN(NSM)
      DIMENSION CD(-2:2),CW(-2:2)
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
      DIMENSION AL(3)
      DATA RT2/1.41421356D0/
C
      WW=2.D0*PI*RF*1.D6
      WC2=(WW/VC)**2
      X=XD(IN)
C
      DO 10 IS=1,NSMAX
         WP(IS)=PZ(IS)*PZ(IS)*AEE*AEE*1.D20/(PA(IS)*AMP*EPS0*WW*WW)
         WC(IS)=PZ(IS)*AEE/(PA(IS)*AMP*WW)
         TN(IS)=AEE*1.D3/(PA(IS)*AMP*VC*VC)
   10 CONTINUE
C
         CALL WFBMAG(XD(IN),YD(IN),BABS,AL)
         CALL WFSDEN(IN,RN,RTPR,RTPP,RZCL)
C
         IF(MODELS.EQ.2) THEN
            RKPR=AL(3)*NPHI/(RR+X)
         ELSEIF(MODELS.EQ.1) THEN
            IF(X.EQ.0.D0) THEN
               RKPR=1.D0
            ELSE
               RKPR=AL(3)*NPHI/X
            ENDIF
         ELSE
            RKPR=AL(3)*RKZ
         ENDIF
         IF(ABS(RKPR).LE.1.D-4) RKPR=1.D-4
C
         IF(MODELD.EQ.2) THEN
            R=1.D0
            DO 20 IS=1,NSMAX
               R=R-WP(IS)*RN(IS)/(1.D0+WC(IS)*BABS)
   20       CONTINUE
            RKPP2=WC2*R-RKPR*RKPR
            IF(RKPP2.GT.0.D0) THEN
               RKPP=SQRT(RKPP2)
            ELSE
               RKPP=0.D0
            ENDIF
         ELSE
            RKPP=0.D0
         ENDIF
C
         DO 100 IS=1,NSMAX
            IF(MODELD.EQ.0.OR.
     &         RTPR(IS).LE.0.D0.OR.
     &         RTPP(IS).LE.0.D0) THEN
               CWP=WP(IS)*RN(IS)/DCMPLX(1.D0,RZCL(IS))
               CWC=WC(IS)*BABS/(DCMPLX(1.D0,RZCL(IS)))
               CDT0= CWP/(1.D0-CWC**2)
               CDX0= CI*CWP*CWC/(1.D0-CWC**2)
               CDP0= CWP
               CDT2= 0.D0
               CDX2= 0.D0
               CDM2= 0.D0
               CDP2= 0.D0
            ELSE
               CWP=WP(IS)*RN(IS)
               IF(BB.GT.0.D0) THEN
                  CWC=WC(IS)*BABS
               ELSE
                  CWC=(0.D0,0.D0)
               ENDIF
               TPR=TN(IS)*RTPR(IS)
               TPP=TN(IS)*RTPP(IS)
               APR=1.D0/(RT2*ABS(RKPR)*SQRT(TPR/WC2))
               IF(BB.GT.0.D0) THEN
                  CTP=TPP/(CWC*CWC*WC2)
               ELSE
                  CTP=0.D0
               ENDIF
               RTP=TPP/TPR
               RTR=1.D0-TPR/TPP
C
               DO 30 IC=-2,2
                  CGZ=(1.D0-DBLE(IC)*CWC)*APR
CCC                  CALL ZETA(CGZ,CZ,CDZ,CDDZ)
                  CALL DSPFN(CGZ,CZ,CDZ)
                  CD(IC)=     -(1.D0-DBLE(IC)*CWC*RTR)*RTP*CZ*APR
                  CW(IC)=0.5D0*(1.D0-DBLE(IC)*CWC*RTR)*RTP*CDZ/RKPR**2
   30          CONTINUE
C
               CDT0=CWP*(1.D0-RTP+0.5D0*(CD(1)+CD(-1)))
               CDP0=CWP*WC2*CW(0)/TPP
               CDX0=CWP*(0.D0,0.5D0)*(CD(1)-CD(-1))
               CDT2=CWP*CTP*(-(CD(1)+CD(-1))+(CD(2)+CD(-2)))
               CDM2=CWP*CTP*(4.0D0*CD(0)-2.D0*(CD(1)+CD(-1)))
               IF(BB.GT.0.D0) THEN
                  CDP2=CWP/(CWC*CWC)*(          -2.D0 *CW( 0)
     &                                +(1.D0/CWC-1.D0)*CW( 1)
     &                                +(1.D0/CWC+1.D0)*CW(-1))
               ELSE
                  CDP2=0.D0
               ENDIF
               CDX2=CWP*(0.D0,1.D0)*CTP
     &                    *(-2.0D0*(CD(1)-CD(-1))+(CD(2)-CD(-2)))
            ENDIF
C
            CDT=CDT0+0.5D0*RKPP*RKPP*CDT2
            CDP=CDP0+0.5D0*RKPP*RKPP*CDP2-CDT
            CDX=CDX0+0.5D0*RKPP*RKPP*CDX2
            CDM=     0.5D0*RKPP*RKPP*CDM2
C
            FX=AL(1)
            FY=AL(2)
            FZ=AL(3)
            CXX= CDT   +CDP*FX*FX+0.5D0*(FY*FY+FZ*FZ)*CDM
            CXY= CDX*FZ+CDP*FX*FY-0.5D0*FX*FY*CDM
            CXZ=-CDX*FY+CDP*FX*FZ-0.5D0*FX*FZ*CDM
            CYX=-CDX*FZ+CDP*FY*FX-0.5D0*FY*FX*CDM
            CYY= CDT   +CDP*FY*FY+0.5D0*(FZ*FZ+FX*FX)*CDM
            CYZ= CDX*FX+CDP*FY*FZ-0.5D0*FY*FZ*CDM
            CZX= CDX*FY+CDP*FZ*FX-0.5D0*FZ*FX*CDM
            CZY=-CDX*FX+CDP*FZ*FY-0.5D0*FZ*FY*CDM
            CZZ= CDT   +CDP*FZ*FZ+0.5D0*(FX*FX+FY*FY)*CDM
C
            CK(1,1,IS)=CXX
            CK(1,2,IS)=CXY
            CK(1,3,IS)=CXZ
            CK(2,1,IS)=CYX
            CK(2,2,IS)=CYY
            CK(2,3,IS)=CYZ
            CK(3,1,IS)=CZX
            CK(3,2,IS)=CZY
            CK(3,3,IS)=CZZ
  100    CONTINUE
C
      RETURN
      END
C
C     ****** LOCAL ELEMENT MATRIX ******
C
      SUBROUTINE CMCALC(IE)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CKL(3,3,3),XE(3),YE(3)
      DIMENSION CK(3,3,NSM)
      REAL*8 A(3),B(3),C(3)
      DIMENSION RL(3),RI(3),CKZ(3)
C
      WC2=(2.D0*PI*RF*1.D6/VC)**2
      RW=2.D0*PI*RF*1.D6
C
      CALL WFNPOS(IE,XE,YE)
      CALL WFABC(IE,A,B,C,S)
C
C     +++++ INITIALIE CM AND CV +++++
C
      DO M=1,3
      DO N=1,3
      DO J=1,4
      DO I=1,4
         CM(I,J,N,M)=(0.D0,0.D0)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO N=1,3
      DO I=1,4
         CMV(I,N)=(0.D0,0.D0)
      ENDDO
      ENDDO
C
C     +++++ SET METRIC +++++
C
      DO K=1,3
         IF(MODELS.EQ.1) THEN
            RL(K)=2.D0*PI*XE(K)
            IF(XE(K).LE.0.D0) THEN
               RI(K)=0.D0
            ELSE
               RI(K)=1.D0/XE(K)
            ENDIF
            CKZ(K)=CI*NPHI*RI(K)
         ELSEIF(MODELS.EQ.2) THEN
            RL(K)=2.D0*PI*(RR+XE(K))
            RI(K)=1.D0/(RR+XE(K))
            CKZ(K)=CI*NPHI*RI(K)
         ELSE
            RL(K)=1.D0
            RI(K)=0.D0
            CKZ(K)=CI*RKZ
         ENDIF
      ENDDO
C
C     +++++ SET DIELECTRIC TENSOR +++++
C
      IF(KELM(IE).EQ.0) THEN
         DO K=1,3
            IK=IELM(K,IE)
            CALL DTENSR(IK,CK)
C
            DO I=1,3
            DO J=1,3
               IF(I.EQ.J) THEN
                  CKL(I,J,K)=1.D0
               ELSE
                  CKL(I,J,K)=0.D0
               ENDIF
               DO IS=1,NSMAX
                  CKL(I,J,K)=CKL(I,J,K)-CK(I,J,IS)
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         XMU=1.D0
      ELSE
         IM=KELM(IE)
         EPSD=EPSDM(IM)
         RMUD=RMUDM(IM)
         SIGD=SIGDM(IM)
C         WRITE(6,'(A,2I5,1P3E12.4)') 'IE,IM,EPSD,RNUD,SIGD=',
C     &                                IE,IM,EPSD,RNUD,SIGD
         DO K=1,3
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
         XMU=1.D0/RMUD
      ENDIF
C
C     +++++ MAIN LOOP +++++
C
      DO 1000 N=1,3
      DO 1000 M=1,3
C
         DO 250 K=1,3
            CP=B(N)*B(M)+C(N)*C(M)
C
            CM(1,1,N,M)=CM(1,1,N,M)
     &                 +XMU*S*RL(K)*(CP*AIF1(K)
     &                              -CKZ(K)*CKZ(K)*AIF3(N,M,K)
     &                              +RI(K) *RI(K) *AIF3(N,M,K))
            CM(1,3,N,M)=CM(1,3,N,M)
     &                 +XMU*S*RL(K)*( 2.D0*RI(K)*CKZ(K)*AIF3(N,M,K))
            CM(3,1,N,M)=CM(3,1,N,M)
     &                 +XMU*S*RL(K)*(-2.D0*CKZ(K)*RI(K)*AIF3(N,M,K))
            CM(3,3,N,M)=CM(3,3,N,M)
     &                 +XMU*S*RL(K)*(CP*AIF1(K)
     &                              -CKZ(K)*CKZ(K)*AIF3(N,M,K)
     &                              +RI(K) *RI(K) *AIF3(N,M,K))
            CM(2,2,N,M)=CM(2,2,N,M)
     &                 +XMU*S*RL(K)*(CP*AIF1(K)
     &                              -CKZ(K)*CKZ(K)*AIF3(N,M,K))
  250    CONTINUE
C
         DO 300 K=1,3
         DO 300 I=1,3
         DO 300 J=1,3
            CM(I,J,N,M)=CM(I,J,N,M)
     &                 -WC2*S*RL(K)*CKL(I,J,K)*AIF3(N,M,K)
  300    CONTINUE
C
         DO 400 K=1,3
         DO 400 I=1,3
            CM(I,4,N,M)=CM(I,4,N,M)
     &                 -WC2*S*RL(K)*(+CI*B(M)  *CKL(I,1,K)*AIF2(N,K)
     &                               +CI*C(M)  *CKL(I,2,K)*AIF2(N,K)
     &                               +CI*CKZ(K)*CKL(I,3,K)*AIF3(N,M,K))
  400    CONTINUE
C
         DO 500 K=1,3
         DO 500 J=1,3
            CM(4,J,N,M)=CM(4,J,N,M)
     &                 -WC2*S*RL(K)*(-CI*B(N)  *CKL(1,J,K)*AIF2(M,K)
     &                               -CI*C(N)  *CKL(2,J,K)*AIF2(M,K)
     &                               +CI*CKZ(K)*CKL(3,J,K)*AIF3(N,M,K))
  500    CONTINUE
C
         DO 600 K=1,3
            CM(4,4,N,M)=CM(4,4,N,M)
     &                 -WC2*S*RL(K)*(
     &                      +B(N)  *B(M)  *CKL(1,1,K)*AIF1(K)
     &                      +B(N)  *C(M)  *CKL(1,2,K)*AIF1(K)
     &                      +B(N)  *CKZ(K)*CKL(1,3,K)*AIF2(M,K)
     &                      +C(N)  *B(M)  *CKL(2,1,K)*AIF1(K)
     &                      +C(N)  *C(M)  *CKL(2,2,K)*AIF1(K)
     &                      +C(N)  *CKZ(K)*CKL(2,3,K)*AIF2(M,K)
     &                      -CKZ(K)*B(M)  *CKL(3,1,K)*AIF2(N,K)
     &                      -CKZ(K)*C(M)  *CKL(3,2,K)*AIF2(N,K)
     &                      -CKZ(K)*CKZ(K)*CKL(3,3,K)*AIF3(N,M,K))
  600    CONTINUE
 1000 CONTINUE
C
C     +++++ ON BOUNDARY: only A-NORMAL is unknown +++++
C                        PHI is known and possibly nonzero
C
      DO M=1,3
         IM=IELM(M,IE)
         IF(KNODW(IM).GE.2.AND.MOD(KNODW(IM),2).EQ.0) THEN
            DO N=1,3
               DO I=1,4
                  CM(I,1,N,M)=ANXN(IM)*CM(I,1,N,M)
     &                       +ANYN(IM)*CM(I,2,N,M)
               ENDDO
            ENDDO
            IF(KNODW(IM).EQ.2) THEN
               PHIL=0.D0
            ELSEIF(KNODW(IM).EQ.4) THEN
               PHIL=PHIW
            ELSEIF(KNODW(IM).EQ.6) THEN
               NV=1
               XMIN=XYVARW(1,NV)
               XMAX=XYVARW(2,NV)
               YMIN=XYVARW(3,NV)
               YMAX=XYVARW(4,NV)
               RLEN0=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
               RLEN1=SQRT((XD(IM)-XMIN)**2+(YD(IM)-YMIN)**2)
               FACT=RLEN1/RLEN0
               PHIL=FACT*PHIW
            ENDIF
            DO N=1,3
               DO I=1,4
                  CMV(I,N)=CMV(I,N)-CM(I,4,N,M)*PHIL
               ENDDO
            ENDDO
         ENDIF
      ENDDO
C
      DO N=1,3
         IN=IELM(N,IE)
         IF(KNODW(IN).GE.2.AND.MOD(KNODW(IN),2).EQ.0) THEN
            DO M=1,3
               DO J=1,4
                  CM(1,J,N,M)=ANXN(IN)*CM(1,J,N,M)
     &                       +ANYN(IN)*CM(2,J,N,M)
               ENDDO
            ENDDO
            CMV(1,N)=ANXN(IN)*CMV(1,N)+ANYN(IN)*CMV(2,N)
         ENDIF
      ENDDO
C
C     +++++ ON AXIS: NPHI = 0 : A-Z and PHI are unknown
C                    NPHI = 1 : A-R + i NPHI A-PHI is unknown
C                    NPHI > 1 : ALL variables are known
C
      IF(MODELS.EQ.1) THEN
         DO M=1,3
            IM=IELM(M,IE)
            IF(MOD(KNODW(IM),2).EQ.1) THEN
               IF(KNODW(IM).EQ.1) THEN
                  IF(NPHI.EQ.0) THEN
                     DO N=1,3
                        DO I=1,4
                           CM(I,1,N,M)=CM(I,2,N,M)
                           CM(I,2,N,M)=CM(I,4,N,M)
                        ENDDO
                     ENDDO
                  ELSEIF(ABS(NPHI).EQ.1) THEN
                     DO N=1,3
                        DO I=1,4
                           CM(I,1,N,M)=CM(I,1,N,M)+CI*NPHI*CM(I,3,N,M)
                        ENDDO
                     ENDDO
                  ENDIF
C
               ELSE
                  IF(NPHI.EQ.0) THEN
                     DO N=1,3
                        DO I=1,4
                           CM(I,1,N,M)=CM(I,2,N,M)
                        ENDDO
                     ENDDO
                  ENDIF
                  IF(KNODW(IM).EQ.3) THEN
                     PHIL=0.D0
                  ELSEIF(KNODW(IM).EQ.5) THEN
                     PHIL=PHIW
                  ELSEIF(KNODW(IM).EQ.7) THEN
                     NV=1
                     XMIN=XYVARW(1,NV)
                     XMAX=XYVARW(2,NV)
                     YMIN=XYVARW(3,NV)
                     YMAX=XYVARW(4,NV)
                     RLEN0=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
                     RLEN1=SQRT((XD(IM)-XMIN)**2+(YD(IM)-YMIN)**2)
                     FACT=RLEN1/RLEN0
                     PHIL=FACT*PHIW
                  ENDIF
                  DO N=1,3
                     DO I=1,4
                        CMV(I,N)=CMV(I,N)-CM(I,4,N,M)*PHIL
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
C
         DO N=1,3
            IN=IELM(N,IE)
            IF(MOD(KNODW(IN),2).EQ.1) THEN
               IF(KNODW(IN).EQ.1) THEN
                  IF(NPHI.EQ.0) THEN
                     DO M=1,3
                        DO J=1,4
                           CM(1,J,N,M)=CM(2,J,N,M)
                           CM(2,J,N,M)=CM(4,J,N,M)
                        ENDDO
                     ENDDO
                     CMV(1,N)=CMV(2,N)
                     CMV(2,N)=CMV(4,N)
                  ELSEIF(ABS(NPHI).EQ.1) THEN
                     DO M=1,3
                        DO J=1,4
                           CM(1,J,N,M)=CM(1,J,N,M)-CI*NPHI*CM(3,J,N,M)
                        ENDDO
                     ENDDO
                     CMV(1,N)=CMV(1,N)-CI*NPHI*CMV(3,N)
                  ENDIF
C
               ELSE
                  IF(NPHI.EQ.0) THEN
                     DO M=1,3
                        DO J=1,4
                           CM(1,J,N,M)=CM(2,J,N,M)
                        ENDDO
                     ENDDO
                     CMV(1,N)=CMV(2,N)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDIF
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
      DIMENSION CJ(3)
      REAL*8 A(3),B(3),C(3)
C
      CNST=2.D0*PI*RF*1.D6*RMU0
C
      DO 11 I=1,MLEN
         CRV(I)=(0.D0,0.D0)
   11 CONTINUE
C
      DO 10 NA=1,NAMAX
         IF(NTYPJ0(NA).GE.0.OR.JNUM0(NA).EQ.2) THEN
            IJMAX=JNUM(NA)+1
         ELSE
            IJMAX=JNUM(NA)
         ENDIF
         DO 20 IJ=1,IJMAX
            IF(JNUM0(NA).EQ.1) THEN
               XM=XJ(1,NA)
               YM=YJ(1,NA)
               IE=JAELM(1,NA)
               CALL WFABC(IE,A,B,C,S)
               CVJ=CNST*CAJ(NA)
               CJ(1)=(0.D0,0.D0)
               CJ(2)=(0.D0,0.D0)
               CJ(3)=CVJ
            ELSEIF(NTYPJ0(NA).LT.0.AND.JNUM0(NA).EQ.2) THEN
               IF(IJ.EQ.1) THEN
                  XM=XJ(IJ,NA)
                  YM=YJ(IJ,NA)
                  IE=JAELM(IJ,NA)
                  CVJ= CNST*CAJ(NA)*CI*XM/NPHI
                  CALL WFABC(IE,A,B,C,S)
                  CJ(1)=(0.D0,0.D0)
                  CJ(2)=(0.D0,0.D0)
                  CJ(3)=CVJ
               ELSEIF(IJ.EQ.IJMAX) THEN
                  XM=XJ(IJ-1,NA)
                  YM=YJ(IJ-1,NA)
                  IE=JAELM(IJ-1,NA)
                  CVJ=-CNST*CAJ(NA)*CI*XM/NPHI
                  CALL WFABC(IE,A,B,C,S)
                  CJ(1)=(0.D0,0.D0)
                  CJ(2)=(0.D0,0.D0)
                  CJ(3)=CVJ
               ELSE
                  X1=XJ(IJ-1,NA)
                  Y1=YJ(IJ-1,NA)
                  X2=XJ(IJ,NA)
                  Y2=YJ(IJ,NA)
                  IE=JAELM(IJ,NA)
                  CALL WFABC(IE,A,B,C,S)
                  XM=0.5D0*(X1+X2)
                  YM=0.5D0*(Y1+Y2)
                  CVJ=CNST*CAJ(NA)
                  CJ(1)=CVJ*(X2-X1)
                  CJ(2)=CVJ*(Y2-Y1)
                  CJ(3)=0.D0
               ENDIF
            ELSEIF(NTYPJ0(NA).GE.0) THEN
               PHL=PHJ0(NA)*PI/180.D0
               IF(IJ.EQ.1) THEN
                  XM=XJ(IJ,NA)
                  YM=YJ(IJ,NA)
                  IE=JAELM(IJ,NA)
                  CVJ=CNST*CAJ(NA)*CI*XM/NPHI
                  CALL WFABC(IE,A,B,C,S)
                  CJ(1)=(0.D0,0.D0)
                  CJ(2)=(0.D0,0.D0)
                  CJ(3)=CVJ
               ELSEIF(IJ.EQ.IJMAX) THEN
                  XM=XJ(IJ-1,NA)
                  YM=YJ(IJ-1,NA)
                  IE=JAELM(IJ-1,NA)
                  CVJ=-CNST*CAJ(NA)*CI*XM/NPHI*EXP(CI*PHL)
                  CALL WFABC(IE,A,B,C,S)
                  CJ(1)=(0.D0,0.D0)
                  CJ(2)=(0.D0,0.D0)
                  CJ(3)=CVJ
               ELSE
                  X1=XJ(IJ-1,NA)
                  Y1=YJ(IJ-1,NA)
                  X2=XJ(IJ,NA)
                  Y2=YJ(IJ,NA)
                  IE=JAELM(IJ,NA)
                  CALL WFABC(IE,A,B,C,S)
                  XM=0.5D0*(X1+X2)
                  YM=0.5D0*(Y1+Y2)
                  CVJ=CNST*CAJ(NA)*EXP(CI*PHL*(YM-ZJH1)/(ZJH2-ZJH1))
                  CJ(1)=CVJ*(X2-X1)
                  CJ(2)=CVJ*(Y2-Y1)
                  CJ(3)=-CVJ*PHL*XM/((ZJH2-ZJH1)*NPHI)*(Y2-Y1)
               ENDIF
C               WRITE(6,'(A,I3,1P6E12.4)') 'CV=',IJ,(CJ(I),I=1,3)
            ELSE
               IF(IJ.EQ.1) GOTO 20
               X1=XJ(IJ-1,NA)
               Y1=YJ(IJ-1,NA)
               X2=XJ(IJ,NA)
               Y2=YJ(IJ,NA)
               IE=JAELM(IJ,NA)
               CALL WFABC(IE,A,B,C,S)
               XM=0.5D0*(X1+X2)
               YM=0.5D0*(Y1+Y2)
               CVJ=CNST*CAJ(NA)
               CJ(1)=CVJ*(X2-X1)
               CJ(2)=CVJ*(Y2-Y1)
               CJ(3)=(0.D0,0.D0)
            ENDIF
C
            IF(MODELS.EQ.1) THEN
               RL=2.D0*PI*XM
               CKZ=CI*NPHI/XM
            ELSEIF(MODELS.EQ.2) THEN
               RL=2.D0*PI*(RR+XM)
               CKZ=CI*NPHI/(RR+XM)
            ELSE
               RL=1.D0
               CKZ=CI*RKZ
            ENDIF
C
            DO 50 N=1,3
               WEIGHT=A(N)+XM*B(N)+YM*C(N)
               IN=IELM(N,IE)
               IF(KNODW(IN).EQ.0) THEN
                  CRV(IBND(IN)  )=CRV(IBND(IN)  )
     &                           +WEIGHT*CJ(1)*RL
                  CRV(IBND(IN)+1)=CRV(IBND(IN)+1)
     &                           +WEIGHT*CJ(2)*RL
                  CRV(IBND(IN)+2)=CRV(IBND(IN)+2)
     &                           +WEIGHT*CJ(3)*RL
                  CRV(IBND(IN)+3)=CRV(IBND(IN)+3)
     &                           -CI*B(N)*CJ(1)*RL
     &                           -CI*C(N)*CJ(2)*RL
     &                           +CI*CKZ *CJ(3)*WEIGHT*RL
               ELSE IF(MOD(KNODW(IN),2).EQ.0) THEN
                  CRV(IBND(IN)  )=CRV(IBND(IN)  )
     &                 +ANXN(IN)*WEIGHT*CJ(1)*RL
     &                 +ANYN(IN)*WEIGHT*CJ(2)*RL
               ELSE
C                  WRITE(6,*) 'NA,IJ,IN=',NA,IJ,IN
               ENDIF
   50       CONTINUE
   20    CONTINUE
   10 CONTINUE
      RETURN
      END
C
C     ******* FRONTAL ELIMINATION WITH LEAST DIAGONAL PIVOTING *******
C
      SUBROUTINE CVSOLV(IERR,NZ)
C
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KFNAMB*32
      DATA ND1/23/
C
      KFNAMB='/tmp/wfx-buffer'
      MTRACK=MBUFM*(4+16)
      OPEN(ND1,FILE=KFNAMB,FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=MTRACK)
      IBUFF=0
C
C     FIND LAST APPEAREANCE OF EACH NODE
C
C
      DO 5 I=1,MLEN
         CSV(I)=0.D0
    5 CONTINUE
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
         CEQ(J,I)=0.D0
   30 CONTINUE
C
 1000 NELL=NELL+1
C
      CALL CMCALC(NELL)
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
         CEQ(KK,LL)=CEQ(KK,LL)+CM(II,JJ,I,J)
   70 CONTINUE
C
      K=0
      DO 75 I=1,3
      DO 75 II=1,JBND(ABS(JELM(I,NELL)))
         K=K+1
         KK=LDEST(K)
         LCO=ABS(LHED(KK))
         CRV(LCO)=CRV(LCO)+CMV(II,I)
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
      CPIVOT=0.D0
      DO 90 L=1,LC
         LPIVC=LPIV(L)
         CPIVA=CEQ(LPIVC,LPIVC)
C_MSP    IF(ABS(CPIVA).GE.ABS(CPIVOT)) THEN
         IF(CDABS(CPIVA).GE.CDABS(CPIVOT)) THEN
            CPIVOT=CPIVA
            LPIVCO=LPIVC
         ENDIF
   90 CONTINUE
C
C     NORMALIZE PIVOTAL ROW
C
      LCO=ABS(LHED(LPIVCO))
C_MSP IF(ABS(CPIVOT).LT.1.D-15) WRITE(6,601)
      IF(CDABS(CPIVOT).LT.1.D-15) THEN 
         WRITE(6,601) NELL
         WRITE(6,*) CPIVOT
      ENDIF
      DO 100 L=1,LCOL
         CQQ(L)=CEQ(LPIVCO,L)/CPIVOT
  100 CONTINUE
      CRHS=CRV(LCO)/CPIVOT
      CRV(LCO)=CRHS
C
C     ELIMINATE THEN DELETE PIVOTAL ROW AND COLUMN
C
      DO 120 K=1,LPIVCO-1
         CFAC=CEQ(K,LPIVCO)
         KRW=ABS(LHED(K))
         CRV(KRW)=CRV(KRW)-CFAC*CRHS
         DO 110 L=1,LPIVCO-1
            CEQ(K,L)=CEQ(K,L)-CFAC*CQQ(L)
  110    CONTINUE
         DO 120 L=LPIVCO+1,LCOL
            CEQ(K,L-1)=CEQ(K,L)-CFAC*CQQ(L)
  120 CONTINUE
      DO 140 K=LPIVCO+1,LCOL
         CFAC=CEQ(K,LPIVCO)
         KRW=ABS(LHED(K))
         CRV(KRW)=CRV(KRW)-CFAC*CRHS
         DO 130 L=1,LPIVCO-1
            CEQ(K-1,L)=CEQ(K,L)-CFAC*CQQ(L)
  130    CONTINUE
         DO 140 L=LPIVCO+1,LCOL
            CEQ(K-1,L-1)=CEQ(K,L)-CFAC*CQQ(L)
  140 CONTINUE
      DO 150 L=LPIVCO+1,LCOL
         LHED(L-1)=LHED(L)
         CQQ (L-1)=CQQ (L)
  150 CONTINUE
      DO 160 K=1,LCOL
         CEQ(K,LCOL)=0.D0
  160 CONTINUE
      DO 170 L=1,LCOL-1
         CEQ(LCOL,L)=0.D0
  170 CONTINUE
      LCOL=LCOL-1
C
C     WRITE PIVOTAL EQUATION ON DISC
C
      IF(IPOS+LCOL.GT.MBUFM) THEN
         IBUFF=IBUFF+1
         WRITE(ND1,REC=IBUFF) MHED,CBUF
         IPOS=0
      ENDIF
      ILEN=ILEN+1
      MLCO(ILEN)=LCO
      MCOL(ILEN)=LCOL
      MPOS(ILEN)=IPOS
      DO 180 L=1,LCOL
         MHED(IPOS+L)=ABS(LHED(L))
         CBUF(IPOS+L)=CQQ(L)
  180 CONTINUE
      IPOS=IPOS+LCOL
C
C     DETERMINE WHETHER TO ASSEMBLE OR BACKSUBSTITUTE
C
      IF(LCOL.GT.1) GOTO 2000
      LCO=ABS(LHED(1))
      CPIVOT=CEQ(1,1)
C_MSP IF(ABS(CPIVOT).LT.1.D-15) WRITE(6,601)
      IF(CDABS(CPIVOT).LT.1.D-15) WRITE(6,601)
      CSV(LCO)=CRV(LCO)/CPIVOT
      IF(NZ.EQ.1.AND.IBUFF.NE.0) WRITE(6,*) '## BUFFER MAX = ',IBUFF
C
C     BACK SUBSTITUTION
C
      DO 200 ILEN=MLEN-1,1,-1
         IF(IPOS.EQ.0) THEN
            READ(ND1,REC=IBUFF) MHED,CBUF
            IBUFF=IBUFF-1
         ENDIF
         LCO =MLCO(ILEN)
         IPOS=MPOS(ILEN)
         CGASH=0.D0
         DO 190 L=1,MCOL(ILEN)
            CGASH=CGASH-CBUF(IPOS+L)*CSV(MHED(IPOS+L))
  190    CONTINUE
         CSV(LCO)=CRV(LCO)+CGASH
  200 CONTINUE
      IERR=0
      RETURN
C
 9000 IERR=2
      WRITE(6,602) IERR
      RETURN
C
  601 FORMAT(1H ,'## FRONT WARNING : SINGULAR OR ILL CONDITIONED'/
     &       1H ,'       AT NELL = ',I5)
  602 FORMAT(1H ,'## FRONT ERROR : ERR =',I5/
     &       1H ,'       LCOL EXCEEDS NMAX')
      END
C
C     ******* ELECTRIC FIELD CALCULATION *******
C
      SUBROUTINE CALFLD
C
      INCLUDE 'wfcomm.inc'
      DIMENSION CE(3),CE0(3)
      DIMENSION DXN(8),DYN(8)
C
      DRN=1.D-6
      DSIN8=SIN(0.125D0*PI)
      DCOS8=COS(0.125D0*PI)
      DXN(1)= DRN*DCOS8
      DXN(2)= DRN*DSIN8
      DXN(3)=-DRN*DSIN8
      DXN(4)=-DRN*DCOS8
      DXN(5)=-DRN*DCOS8
      DXN(6)=-DRN*DSIN8
      DXN(7)= DRN*DSIN8
      DXN(8)= DRN*DCOS8
      DYN(1)= DRN*DSIN8
      DYN(2)= DRN*DCOS8
      DYN(3)= DRN*DCOS8
      DYN(4)= DRN*DSIN8
      DYN(5)=-DRN*DSIN8
      DYN(6)=-DRN*DCOS8
      DYN(7)=-DRN*DCOS8
      DYN(8)=-DRN*DSIN8
C
      DO 1 I=1,3
      DO 1 IN=1,NNOD
         CEFK(I,IN)=0.D0
    1 CONTINUE
C
      DO 10 IN=1,NNOD
         IF(KNODW(IN).EQ.0) THEN
            CAF(1,IN)=CSV(IBND(IN))
            CAF(2,IN)=CSV(IBND(IN)+1)
            CAF(3,IN)=CSV(IBND(IN)+2)
            CAF(4,IN)=CSV(IBND(IN)+3)
         ELSEIF(MOD(KNODW(IN),2).EQ.0) THEN
            CAF(1,IN)=ANXN(IN)*CSV(IBND(IN))
            CAF(2,IN)=ANYN(IN)*CSV(IBND(IN))
            CAF(3,IN)=0.D0
            IF(KNODW(IN).EQ.2) THEN
               PHIL=0.D0
            ELSEIF(KNODW(IN).EQ.4) THEN
               PHIL=PHIW
            ELSEIF(KNODW(IN).EQ.6) THEN
               NV=1
               XMIN=XYVARW(1,NV)
               XMAX=XYVARW(2,NV)
               YMIN=XYVARW(3,NV)
               YMAX=XYVARW(4,NV)
               RLEN0=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
               RLEN1=SQRT((XD(IN)-XMIN)**2+(YD(IN)-YMIN)**2)
               FACT=RLEN1/RLEN0
               PHIL=FACT*PHIW
            ENDIF
            CAF(4,IN)=PHIL
         ELSEIF(KNODW(IN).EQ.1) THEN
            IF(NPHI.EQ.0) THEN
               CAF(1,IN)=0.D0
               CAF(2,IN)=CSV(IBND(IN))
               CAF(3,IN)=0.D0
               CAF(4,IN)=CSV(IBND(IN)+1)
            ELSEIF(ABS(NPHI).EQ.1) THEN
               CAF(1,IN)=CSV(IBND(IN))
               CAF(2,IN)=0.D0
               CAF(3,IN)=CI*NPHI*CSV(IBND(IN))
               CAF(4,IN)=0.D0
            ELSE
               CAF(1,IN)=0.D0
               CAF(2,IN)=0.D0
               CAF(3,IN)=0.D0
               CAF(4,IN)=0.D0
            ENDIF
         ELSEIF(MOD(KNODW(IN),2).EQ.1) THEN
            IF(NPHI.EQ.0) THEN
               CAF(1,IN)=0.D0
               CAF(2,IN)=CSV(IBND(IN))
               CAF(3,IN)=0.D0
               IF(KNODW(IN).EQ.3) THEN
                  PHIL=0.D0
               ELSEIF(KNODW(IN).EQ.5) THEN
                  PHIL=PHIW
               ELSEIF(KNODW(IN).EQ.7) THEN
                  NV=1
                  XMIN=XYVARW(1,NV)
                  XMAX=XYVARW(2,NV)
                  YMIN=XYVARW(3,NV)
                  YMAX=XYVARW(4,NV)
                  RLEN0=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
                  RLEN1=SQRT((XD(IN)-XMIN)**2+(YD(IN)-YMIN)**2)
                  FACT=RLEN1/RLEN0
                  PHIL=FACT*PHIW
               ENDIF
               CAF(4,IN)=PHIL
            ELSE
               CAF(1,IN)=0.D0
               CAF(2,IN)=0.D0
               CAF(3,IN)=0.D0
               CAF(4,IN)=0.D0
            ENDIF
         ELSE
            WRITE(6,*) 'XX INVALID KBND: KNODW(',IN,')= ',KNODW(IN)
         ENDIF
   10 CONTINUE
C
      IE=1
      DO 30 IN=1,NNOD
         X0=XD(IN)
         Y0=YD(IN)
         CE0(1)=(0.D0,0.D0)
         CE0(2)=(0.D0,0.D0)
         CE0(3)=(0.D0,0.D0)
         ICOUNT=0
         DO I=1,8
            X=X0+DXN(I)
            Y=Y0+DYN(I)
            IES=IE
            CALL WFFEP(X,Y,IE)
            IF(IE.NE.0) THEN
               ICOUNT=ICOUNT+1
               CALL FIELDE(IE,X,Y,CE)
               CE0(1)=CE0(1)+CE(1)
               CE0(2)=CE0(2)+CE(2)
               CE0(3)=CE0(3)+CE(3)
            ELSE
               IE=IES
            ENDIF
         ENDDO
         IF(ICOUNT.EQ.0) THEN
            WRITE(6,*) 'XX WFFEP ERROR : IN,X,Y = ',IN,X,Y
         ELSE
C            IF(ICOUNT.NE.8) WRITE(6,*) 'IN,ICOUNT=',IN,ICOUNT
            CE0(1)=CE0(1)/DBLE(ICOUNT)
            CE0(2)=CE0(2)/DBLE(ICOUNT)
            CE0(3)=CE0(3)/DBLE(ICOUNT)
         ENDIF
         CEFK(1,IN)=CE0(1)
         CEFK(2,IN)=CE0(2)
         CEFK(3,IN)=CE0(3)
   30 CONTINUE
C
         DO 210 IN=1,NNOD
            CEF(1,IN)=CEF(1,IN)+CEFK(1,IN)
            CEF(2,IN)=CEF(2,IN)+CEFK(2,IN)
            CEF(3,IN)=CEF(3,IN)+CEFK(3,IN)
  210    CONTINUE
C
      RETURN
      END
C
C     ****** POWER ABSORPTION ******
C
      SUBROUTINE PWRABS(NZ)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION EPWRS(NSM),TPWRS(NSM)
      DIMENSION XE(3),YE(3),CK(3,3,NSM)
      REAL*8 A(3),B(3),C(3)
      DIMENSION RL(3),RI(3),CKZ(3),V(3)
C
      CNST=0.5D0*CI*EPS0*2.D0*PI*RF*1.D6
      DR=1.D0/(NRMAX-1)
C
      TPWR=0.D0
      DO NS=1,NSMAX
         TPWRS(NS)=0.D0
      ENDDO
C
      DO IEDO=1,NELM
         IE=IEDO
C
         DO NS=1,NSMAX
            EPWRS(NS)=0.D0
         ENDDO
C
         IF(KELM(IE).EQ.0) THEN
C
         CALL WFNPOS(IE,XE,YE)
         CALL WFABC(IE,A,B,C,S)
C
         DO K=1,3
            IF(MODELS.EQ.1) THEN
               RL(K)=2.D0*PI*XE(K)
               IF(XE(K).LE.0.D0) THEN
                  RI(K)=0.D0
               ELSE
                  RI(K)=1.D0/XE(K)
               ENDIF
               CKZ(K)=CI*NPHI*RI(K)
            ELSEIF(MODELS.EQ.2) THEN
               RL(K)=2.D0*PI*(RR+XE(K))
               RI(K)=1.D0/(RR+XE(K))
               CKZ(K)=CI*NPHI*RI(K)
            ELSE
               RL(K)=1.D0
               RI(K)=0.D0
               CKZ(K)=CI*RKZ
            ENDIF
         ENDDO
C
         DO K=1,3
            IK=IELM(K,IE)
            CALL DTENSR(IK,CK)
         DO NS=1,NSMAX
C
            DO N=1,3
               IN=IELM(N,IE)
            DO M=1,3
               IM=IELM(M,IE)
C
               DO I=1,3
               DO J=1,3
                  EPWRS(NS)=EPWRS(NS)
     &                     +DBLE(CNST*DCONJG(CI*CAF(I,IN))*RL(K)*S
     &                     *CK(I,J,NS)*AIF3(N,M,K)*CI*CAF(J,IM))
               ENDDO
               ENDDO
C
               DO I=1,3
                  EPWRS(NS)=EPWRS(NS)
     &                     +DBLE(CNST*DCONJG(CAF(I,IN))*RL(K)*S
     &                     *(+CI*B(M)  *CK(I,1,NS)*AIF2(N,K) 
     &                       +CI*C(M)  *CK(I,2,NS)*AIF2(N,K)
     &                       +CI*CKZ(K)*CK(I,3,NS)*AIF3(N,M,K))
     &                     *CAF(4,IM))
               ENDDO
C
               DO J=1,3
                  EPWRS(NS)=EPWRS(NS)
     &                     +DBLE(CNST*DCONJG(CAF(4,IN))*RL(K)*S
     &                     *(-CI*B(N)  *CK(1,J,NS)*AIF2(M,K)
     &                       -CI*C(N)  *CK(2,J,NS)*AIF2(M,K)
     &                       +CI*CKZ(K)*CK(3,J,NS)*AIF3(N,M,K))
     &                     *CAF(J,IM))
               ENDDO
C
               EPWRS(NS)=EPWRS(NS)
     &                  +DBLE(CNST*DCONJG(CAF(4,IN))*RL(K)*S
     &                  *(+B(N)  *B(M)  *CK(1,1,NS)*AIF1(K)
     &                    +B(N)  *C(M)  *CK(1,2,NS)*AIF1(K)
     &                    +B(N)  *CKZ(K)*CK(1,3,NS)*AIF2(M,K)
     &                    +C(N)  *B(M)  *CK(2,1,NS)*AIF1(K)
     &                    +C(N)  *C(M)  *CK(2,2,NS)*AIF1(K)
     &                    +C(N)  *CKZ(K)*CK(2,3,NS)*AIF2(M,K)
     &                    -CKZ(K)*B(M)  *CK(3,1,NS)*AIF2(N,K)
     &                    -CKZ(K)*C(M)  *CK(3,2,NS)*AIF2(N,K)
     &                    -CKZ(K)*CKZ(K)*CK(3,3,NS)*AIF3(N,M,K))
     &                  *CAF(4,IM))
            ENDDO
            ENDDO
C
         ENDDO
         ENDDO
C         IF(IE.LE.20) 
C     &   WRITE(6,'(A,I5,1P3E12.4)')
C     &        'IE,EPWRS=',IE,EPWRS(1),EPWRS(2),EPWRS(3)
C
         DO I=1,3
            IN=IELM(I,IE)
            IF(MODELS.EQ.1) THEN
               CALL WFSPSI(XD(IN),YD(IN),PSI)
               IF(PSI.LE.0.D0) THEN
                  RP=0.D0
               ELSE
                  RP=SQRT(PSI)/DR
               ENDIF
            ELSE
               RP=ABS(XD(IN))/DR
            ENDIF
            NR1=INT(RP+1.D0)
            R1=(NR1-RP)*S/3.D0
            R2=(RP-NR1+1.D0)*S/3.D0
            IF(NR1.GT.NRMAX) NR1=NRMAX
            NR2=NR1+1
            IF(NR2.GT.NRMAX) NR2=NRMAX
            DO NS=1,NSMAX
               PWRR(NR1)    =PWRR(NR1)    +EPWRS(NS)*R1
               PWRRS(NR1,NS)=PWRRS(NR1,NS)+EPWRS(NS)*R1
               PWRR(NR2)    =PWRR(NR2)    +EPWRS(NS)*R2
               PWRRS(NR2,NS)=PWRRS(NR2,NS)+EPWRS(NS)*R2
            ENDDO
         ENDDO
C
         DO N=1,3
            V(N)=0.D0
            DO K=1,3
               V(N)=V(N)+RL(K)*S*AIF2(N,K)
            ENDDO
         ENDDO
C
         VTOT=0.D0
         DO N=1,3
            VTOT=VTOT+V(N)
         ENDDO
C         IF(IE.LE.20) 
C     &   WRITE(6,'(A,I5,1P4E12.4)')
C     &        'IE,VTOT,V=',
C     &        IE,VTOT,V(1),V(2),V(3)
C
         EPWR=0.D0
         DO NS=1,NSMAX
            EPWR=EPWR+EPWRS(NS)
         ENDDO
C
         PWRE(IE)=EPWR
         DO NS=1,NSMAX
            PWRES(IE,NS)=EPWRS(NS)
         ENDDO
C
         TPWR   =TPWR   +EPWR
         DO NS=1,NSMAX
            TPWRS(NS)=TPWRS(NS)+EPWRS(NS)
         ENDDO
C
         DO N=1,3
            IN=IELM(N,IE)
            PWRN(IN)    =PWRN(IN)    +EPWR   *V(N)/VTOT
            DO NS=1,NSMAX
               PWRNS(IN,NS)=PWRNS(IN,NS)+EPWRS(NS)*V(N)/VTOT
            ENDDO
         ENDDO
C
      ENDIF
      ENDDO
C
      PWR=PWR+TPWR
      DO NS=1,NSMAX
         PWRS(NS)=PWRS(NS)+TPWRS(NS)
      ENDDO
      PWRF(NZ)=PWR
      DO NS=1,NSMAX
         PWRFS(NZ,NS)=PWRS(NS)
      ENDDO
      RETURN
      END
C
C     ****** RADIATED POWER ******
C
      SUBROUTINE PWRRAD
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CJ(3),CE(3),CIMPK(NAM)
      REAL*8 A(3),B(3),C(3)
C
      DO 100 NA=1,NAMAX
         IF(NTYPJ0(NA).NE.0.OR.JNUM0(NA).EQ.2) THEN
            IJMAX=JNUM(NA)+1
         ELSE
            IJMAX=JNUM(NA)
         ENDIF
         CIMPK(NA)=0.D0
         DO 50 IJ=1,IJMAX
            IF(JNUM(NA).EQ.1) THEN
               XM=XJ(1,NA)
               YM=YJ(1,NA)
               IE=JAELM(1,NA)
               CALL WFABC(IE,A,B,C,S)
               CVJ=CAJ(NA)
               CJ(1)=(0.D0,0.D0)
               CJ(2)=(0.D0,0.D0)
               CJ(3)=CVJ
            ELSEIF(NTYPJ0(NA).NE.0.OR.JNUM0(NA).EQ.2) THEN
               PHL=PHJ0(NA)*PI/180.D0
               IF(IJ.EQ.1) THEN
                  IF(NTYPJ0(NA).EQ.0.OR.
     &               NTYPJ0(NA).EQ.2) THEN
                     XM=XJ(IJ,NA)
                     YM=YJ(IJ,NA)
                     IE=JAELM(IJ,NA)
                     CVJ=CAJ(NA)*CI*XM/NPHI
                     CALL WFABC(IE,A,B,C,S)
                     CJ(1)=(0.D0,0.D0)
                     CJ(2)=(0.D0,0.D0)
                     CJ(3)=CVJ
                  ELSE
                     GOTO 50
                  ENDIF
               ELSEIF(IJ.EQ.IJMAX) THEN
                  IF(NTYPJ0(NA).EQ.0.OR.
     &               NTYPJ0(NA).EQ.1) THEN
                     XM=XJ(IJ-1,NA)
                     YM=YJ(IJ-1,NA)
                     IE=JAELM(IJ-1,NA)
                     CVJ=-CAJ(NA)*CI*XM/NPHI*EXP(CI*PHL)
                     CALL WFABC(IE,A,B,C,S)
                     CJ(1)=(0.D0,0.D0)
                     CJ(2)=(0.D0,0.D0)
                     CJ(3)=CVJ
                  ENDIF
               ELSE
                  X1=XJ(IJ-1,NA)
                  Y1=YJ(IJ-1,NA)
                  X2=XJ(IJ,NA)
                  Y2=YJ(IJ,NA)
                  IE=JAELM(IJ,NA)
                  CALL WFABC(IE,A,B,C,S)
                  XM=0.5D0*(X1+X2)
                  YM=0.5D0*(Y1+Y2)
                  CVJ=CAJ(NA)*EXP(CI*PHL*(YM-ZJH1)/(ZJH2-ZJH1))
                  CJ(1)=CVJ*(X2-X1)
                  CJ(2)=CVJ*(Y2-Y1)
                  CJ(3)=-CVJ*PHL*XM/((ZJH2-ZJH1)*NPHI)*(Y2-Y1)
               ENDIF
            ELSE
               IF(IJ.EQ.1) GOTO 50
               X1=XJ(IJ-1,NA)
               Y1=YJ(IJ-1,NA)
               X2=XJ(IJ,NA)
               Y2=YJ(IJ,NA)
               IE=JAELM(IJ,NA)
               CALL WFABC(IE,A,B,C,S)
               XM=0.5D0*(X1+X2)
               YM=0.5D0*(Y1+Y2)
               CVJ=CAJ(NA)
               CJ(1)=CVJ*(X2-X1)
               CJ(2)=CVJ*(Y2-Y1)
               CJ(3)=(0.D0,0.D0)
            ENDIF
            IF(MODELS.EQ.1) THEN
               RL=2.D0*PI*XM
               CKZ=CI*NPHI/XM
            ELSEIF(MODELS.EQ.2) THEN
               RL=2.D0*PI*(RR+XM)
               CKZ=CI*NPHI/(RR+XM)
            ELSE
               RL=1.D0
               CKZ=CI*RKZ
            ENDIF
            DO 20 N=1,3
               WEIGHT=A(N)+XM*B(N)+YM*C(N)
               IN=IELM(N,IE)
               CE(1)=CI*(WEIGHT*CAF(1,IN)+CI*B(N)*CAF(4,IN))
               CE(2)=CI*(WEIGHT*CAF(2,IN)+CI*C(N)*CAF(4,IN))
               CE(3)=CI*(WEIGHT*CAF(3,IN)+CI*CKZ*WEIGHT*CAF(4,IN))
               DO 30 I=1,3
                  CIMPK(NA)=CIMPK(NA)
     &                 -0.5D0*DCONJG(CE(I))*CJ(I)*RL
   30          CONTINUE
   20       CONTINUE
   50    CONTINUE
  100 CONTINUE
C
      DO NA=1,NAMAX
         CIMP(NA)=CIMP(NA)+CIMPK(NA)
         CTIMP   =CTIMP   +CIMPK(NA)
      ENDDO
C
      DO NA=1,NAMAX
         CIMP(NA)=CIMP(NA)/CAJ(NA)
      ENDDO
C
      RETURN
      END
C
C     ****** INITIALIZE EFIELD AND POWER ******
C
      SUBROUTINE INITEP
C
      INCLUDE 'wfcomm.inc'
C
      DO 10 I=1,3
      DO 10 IN=1,NNOD
         CEF(I,IN)=0.D0
   10 CONTINUE
C
      PWR=0.D0
      DO IN=1,NNOD
         PWRN(IN)=0.D0
      ENDDO

      DO NS=1,NSMAX
         PWRS(NS)=0.D0
         DO IN=1,NNOD
            PWRNS(IN,NS)=0.D0
         ENDDO
      ENDDO
C
      CTIMP=(0.D0,0.D0)
      DO NA=1,NAMAX
         CIMP(NA)=(0.D0,0.D0)
      ENDDO
C
      DO NR=1,NRMAX
         PWRR(NR)=0.D0
         DO NS=1,NSMAX
            PWRRS(NR,NS)=0.D0
         ENDDO
      ENDDO
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
      IF(PIN.EQ.0.D0.OR.ABS(PWR).EQ.0.D0) THEN
         FACT=1.D0
      ELSEIF(PIN.GT.0.D0) THEN
         FACT=PIN/ABS(PWR)
      ELSE
         FACT=-PIN/(ABS(PWR)-PIN)
      ENDIF
      FACTR=SQRT(FACT)
C
      PWR=FACT*PWR
      CTIMP=FACT*CTIMP
C
      DO NS=1,NSMAX
         PWRS(NS)=FACT*PWRS(NS)
      ENDDO
      DO IN=1,NNOD
C         IF(IN.LE.20) WRITE(6,'(A,I5,1P3E12.4)')
C     &        'IN,FACT,PWRN(IN),VNOD(IN)=',IN,FACT,PWRN(IN),VNOD(IN)
         PWRN(IN)=FACT*PWRN(IN)/VNOD(IN)
         DO NS=1,NSMAX
            PWRNS(IN,NS)=FACT*PWRNS(IN,NS)/VNOD(IN)
         ENDDO
      ENDDO
C
      DO NZ=1,NZMAX
         PWRF(NZ)=FACT*PWRF(NZ)
         DO NS=1,NSMAX
            PWRFS(NZ,NS)=FACT*PWRFS(NZ,NS)
         ENDDO
      ENDDO
C
      PNMAX=PWRN(1)
      DO IN=2,NNOD
         PNMAX=MAX(PNMAX,PWRN(IN))
      ENDDO
C
      ETMAX=0.D0
      DO 70 I=1,3
         EMAX(I)=0.D0
   70 CONTINUE
C
      DO 90 IN=1,NNOD
         DO 80 J=1,3
            CEF(J,IN)=FACTR*CEF(J,IN)
            EABS(J)=ABS(CEF(J,IN))
            EMAX(J)=MAX(EMAX(J),EABS(J))
   80    CONTINUE
         ETOTL =SQRT(EABS(1)*EABS(1)+EABS(2)*EABS(2)+EABS(3)*EABS(3))
         ETMAX =MAX(ETMAX,ETOTL)
   90 CONTINUE
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
      DIMENSION P1(NSM),P2(NSM),P3(NSM),P4(NSM)
C
      IF(NPRINT.LT.1) RETURN
C
      WRITE(6,110) (EMAX(I),I=1,3),ETMAX,PNMAX
  110 FORMAT(1H ,'EXMAX  =',1PE12.4
     &      ,3X ,'EYMAX  =',1PE12.4
     &      ,3X ,'EZMAX  =',1PE12.4/
     &       1H ,'EMAX   =',1PE12.4
     &      ,3X ,'PNMAX  =',1PE12.4)
C
      WRITE(6,120) DBLE(CTIMP),PWR
  120 FORMAT(1H ,'RADIATED POWER =',1PE12.4/
     &       1H ,'ABSORBED POWER =',1PE12.4,
     &        5X,' P(1/4)  P(2/4)  P(3/4)  P(4/4)  P(RB)')
C
      NR1=INT(NRMAX*0.25)
      NR2=INT(NRMAX*0.5)
      NR3=INT(NRMAX*0.75)
      NR4=NRMAX-1
C
      DO 5 IS=1,NSMAX
         P=0.D0
         DO 1 NR=1,NR1
            P=P+PWRRS(NR,IS)
    1    CONTINUE
         P1(IS)=P
         DO 2 NR=NR1+1,NR2
            P=P+PWRRS(NR,IS)
    2    CONTINUE
         P2(IS)=P
         DO 3 NR=NR2+1,NR3
            P=P+PWRRS(NR,IS)
    3    CONTINUE
         P3(IS)=P
         DO 4 NR=NR3+1,NR4
            P=P+PWRRS(NR,IS)
    4    CONTINUE
         P4(IS)=P
C
         IF(ABS(PWR).GT.1.D-32) THEN
            WRITE(6,125) IS,PWRS(IS),P1(IS)/PWR,P2(IS)/PWR,
     &                   P3(IS)/PWR,P4(IS)/PWR,PWRS(IS)/PWR
  125       FORMAT(1H ,'      PABS(',I2,') =',1PE12.4,
     &              5X,0P5F8.4)
         ELSE
            WRITE(6,126) IS,PWRS(IS)
  126       FORMAT(1H ,'      PABS(',I2,') =',1PE12.4)
         ENDIF
    5 CONTINUE
C
      WRITE(6,130)
  130 FORMAT(1H ,' I JNUM', '  AJ(I)','  APH(I)','  AWD(I)',
     &                     ' APOS(I)',' XJ(I)','  YJ(I)',
     &           8X,'LOADING IMP.[ohm]')
      DO 10 IA=1,NAMAX
         WRITE(6,140) IA,JNUM(IA),AJ(IA),APH(IA),AWD(IA),APOS(IA),
     &                            XJ(1,IA),YJ(1,IA),CIMP(IA)
  140    FORMAT(1H ,I2,I3,0PF8.2,F8.2,1X,4F7.4,2X,'(',1P2E12.4,')')
C  140    FORMAT(1H ,I3,I3,0PF7.4,F7.2,4F7.4,2X,'(',1P2E12.4,')')
   10 CONTINUE
C
      IF(NZMAX.GT.1) THEN
         WRITE(6,145) (IS,IS=1,NSMAX)
  145    FORMAT(1H ,'   NZ','       RKZ',5X,'   REAL(CFJ)   IMAG(CFJ)',
     &              5X,'        PABS')
         DO 20 NZ=1,NZMAX
            WRITE(6,146) NZ,RKZF(NZ),CAJF(NZ,1),PWRF(NZ)
  146       FORMAT(1H ,I5,F10.4,5X,1P2E12.4,5X,1PE12.4)
   20    CONTINUE
         WRITE(6,147) (IS,IS=1,NSMAX)
  147    FORMAT(1H ,'   NZ',3X,5('     PABS(',I1,')':))
         DO 21 NZ=1,NZMAX
         WRITE(6,148) NZ,(PWRFS(NZ,IS),IS=1,NSMAX)
  148    FORMAT(1H ,I5,3X,5(1PE12.4))
   21    CONTINUE
      ENDIF
C
      IF(NPRINT.LT.2) RETURN
C
      WRITE(6,150) (I,KNODW(I),PWRN(I),(CEF(J,I),J=1,3),I=1,NNOD)
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
      WRITE(6,110) NNOD
  110 FORMAT(1H0,'NODE DATA     : #### NNOD =',I5,' ####'/
     &       1H ,2('  NNOD','  KBND','  IBND','   MDF',
     &       9X,'X',14X,'Y',9X))
      WRITE(6,115) (I,KNODW(I),IBND(I),JBND(I),XD(I),YD(I),
     &              I=1,NNOD)
  115 FORMAT((1H ,2(4I6,2X,1P2E15.7,2X)))
C
      WRITE(6,120) NELM,(I,(IELM(J,I),J=1,3),I=1,NELM)
  120 FORMAT(1H0,'ELEMENT DATA  : #### NELM =',I5,' ####'/
     &      (1H ,5(I6,'(',3I5,')',2X)))
C
      WRITE(6,125) NELM,(I,(JELM(J,I),J=1,3),I=1,NELM)
  125 FORMAT(1H0,'NOP     DATA  : #### NELM =',I5,' ####'/
     &      (1H ,5(I6,'(',3I5,')',2X)))
C
      WRITE(6,130) NBDYW
  130 FORMAT(1H0,'BOUNDARY DATA : #### NBDYW =',I5,' ####'/
     &       1H ,3('  NO.',' IBDYW',2X))
      WRITE(6,135) (I,IBDYW(I),I=1,NBDYW)
  135 FORMAT((1H ,3(2I5,2X)))
C
      DO 10 NA=1,NAMAX
         WRITE(6,140) NA,JNUM0(NA)
  140    FORMAT(1H0,'ORIGINAL ANTENNA DATA : NA =',I5,' JNUM0 =',I5/
     &          1H ,3('  NO.',13X,' XJ0',11X,' YJ0',6X))
         WRITE(6,150) (I,XJ0(I,NA),YJ0(I,NA),I=1,JNUM0(NA))
  150    FORMAT((1H ,3(I5,8X,1P2E15.7)))
C
         WRITE(6,154) NA,JNUM(NA)
  154    FORMAT(1H0,'MODIFIED ANTENNA DATA : NA =',I5,' JNUM  =',I5/
     &          1H ,3('  NO.',' JELM',8X,' JX ',11X,' JY ',6X))
         WRITE(6,156) (I,JAELM(I,NA),XJ(I,NA),YJ(I,NA),I=1,JNUM(NA))
  156    FORMAT((1H ,3(2I5,3X,1P2E15.7)))
   10 CONTINUE
C
      RETURN
      END
