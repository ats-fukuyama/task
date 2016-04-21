C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 2) ******
C
C     ******* LOCAL PLASMA PARAMETERS *******
C
      SUBROUTINE W1DSPQ(ICL)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1PRF1/ PROFB(NXPM),PROFPN(NXPM,ISM),PROFPU(NXPM,ISM),
     &                PROFTR(NXPM,ISM),PROFTP(NXPM,ISM)
      COMMON /W1ZETA/ CZ(NXPM,NHM),CDZ(NXPM,NHM),CDDZ(NXPM,NHM),
     &                GZ(NXPM,NHM)
      PARAMETER (NCLW=2*MATLM+1)
      COMMON /W1QCLM/ CL(3,3,4,NCLM),NCL(NCLW,NXPM,ISM)
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
      COMMON /W1QSST/ XM(NXPM),YX(NXPM),YK(NXPM),
     &                CS0(NXPM),CS1(NXPM),CS2(NXPM),CS3(NXPM)
C
      DIMENSION  X1(4),CSB(3,3,4)
C
      MATL=1
      ICL=0
C
      DO 10 IL=1,4
      DO 10 IA=1,3
      DO 10 IB=1,3
      DO 10 I=1,NCLM
         CL(IA,IB,IL,I)=(0.D0,0.D0)
   10 CONTINUE
      DO 20 IS=1,ISMAX
      DO 20 J=1,2*MATLM+1
      DO 20 I=1,NXP
         NCL(J,I,IS)=0
   20 CONTINUE
C
      RT2  = SQRT ( 2.D0 )
      RW  = 2.D6*PI*RF
C
      DO 1000 IS=1,ISMAX
         FWP = 1.D20*AEE*AEE*PZ(IS)*PZ(IS)/(AMM*PA(IS)*EPS0*RW*RW)
         FWC = AEE*PZ(IS)*BB/(AMM*PA(IS))
         FVT = AEE*1.D3/(AMM*PA(IS))
         DO 100 NX = 1 , NXP
            RKPR=RKZ
            WC = FWC*PROFB(NX)
            UD = SQRT(FVT*PROFPU(NX,IS))
            AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,IS))
C 
*VDIR NOVECTOR
         DO 100 NC=1,2*ABS(IHARM(IS))+1
            NN=NC-ABS(IHARM(IS))-1
            ARG=(RW-NN*WC)/AKPR
            GZ(NX,NC)= ARG
            CZ(NX,NC)= ARG
  100    CONTINUE
C
         DO 200 NC=1,2*ABS(IHARM(IS))+1
            CALL DSPFNA(CZ(1,NC),CDZ(1,NC),CDDZ(1,NC),NXP)
  200    CONTINUE
C
         DO 400 NC=1,2*ABS(IHARM(IS))+1
            NN=NC-ABS(IHARM(IS))-1
            DO 300 NX=1,NXP
               RKPR = RKZ
               WC = FWC*PROFB(NX)
               UD = SQRT(FVT*PROFPU(NX,IS))
               AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,IS))
               RT = PROFTP(NX,IS)/PROFTR(NX,IS)
               XM(NX)=XAM(NX)
               YX(NX)=WC/SQRT(FVT*PROFTP(NX,IS))
               YK(NX)=FWP*PROFPN(NX,IS)*ABS(YX(NX))*RW/AKPR
               CS0(NX)=GZ(NX,NC)*CDZ(NX,NC)
               CS1(NX)=CZ(NX,NC)+0.5D0*(1.D0-RT)*AKPR*CDZ(NX,NC)/RW
               CS2(NX)=(RT+(1.D0-RT)*NN*WC/RW)*CDZ(NX,NC)
     &               /SQRT(2.D0*RT)
               CS3(NX)=(1.D0-1.D0/RT)*CS0(NX)*WC/RW
  300       CONTINUE
         DO 400 I=1,NXP-1
            DELTAX=ABS(XA(I)-XA(I+1))
            IF(MWIDTH.GT.MATLM) MWIDTH=MATLM
            IF(MWIDTH.LT.2)     MWIDTH=2
            MATL=MAX(MATL,MWIDTH)
         DO 400 J=MAX(1,I-MWIDTH+1),MIN(I+MWIDTH-1,NXP-1)
            DXA=XA(I+1)-XA(I)
            DXB=XA(J+1)-XA(J)
            X1(1)=XA(I)-XA(J)
            X1(2)=XA(I+1)-XA(J)
            X1(3)=XA(I)-XA(J+1)
            X1(4)=XA(I+1)-XA(J+1)
            JJ=J-I+MATLM+1
            IF(I.EQ.J) THEN
               X1(1)= 0.D0
               X1(4)= 0.D0
            ELSEIF(I+1.EQ.J) THEN
               X1(2)=0.D0
            ELSEIF(I.EQ.J+1) THEN
               X1(3)=0.D0
            ENDIF
C
            CALL W1QLNI(0.5D0*(XA(I  )+XA(J  )),
     &                  VXA,VKA,CT0A,CT1A,CT2A,CT3A)
            CALL W1QLNI(0.5D0*(XA(I+1)+XA(J  )),
     &                  VXB,VKB,CT0B,CT1B,CT2B,CT3B)
            CALL W1QLNI(0.5D0*(XA(I  )+XA(J+1)),
     &                  VXC,VKC,CT0C,CT1C,CT2C,CT3C)
            CALL W1QLNI(0.5D0*(XA(I+1)+XA(J+1)),
     &                  VXD,VKD,CT0D,CT1D,CT2D,CT3D)
            VX =0.25D0*(VXA +VXB +VXC +VXD )
            VK =0.25D0*(VKA +VKB +VKC +VKD )
            CT0=0.25D0*(CT0A+CT0B+CT0C+CT0D)
            CT1=0.25D0*(CT1A+CT1B+CT1C+CT1D)
            CT2=0.25D0*(CT2A+CT2B+CT2C+CT2D)
            CT3=0.25D0*(CT3A+CT3B+CT3C+CT3D)
            CALL W1QCAL(X1,DXA,DXB,VX,CT0,CT1,CT2,CT3,CSB,NN)
C
            IF(NCL(JJ,I,IS).EQ.0) THEN
               ICL=ICL+1
               NCL(JJ,I,IS)=ICL
               ICLS=ICL
            ELSE
               ICLS=NCL(JJ,I,IS)
            ENDIF
            IF(ICL.GT.NCLM) THEN
               WRITE(6,*) 'XX ICL.GT.NCLM AT IS,I,J = ',IS,',',I,',',J
               ICL=1
               NCL(JJ,I,IS)=ICL
               ICLS=ICL
            ENDIF
         DO 400 IL=1,4
         DO 400 IA=1,3
         DO 400 IB=1,3
            CL(IA,IB,IL,ICLS)=CL(IA,IB,IL,ICLS)
     &                       +VK*CSB(IA,IB,IL)*DXA*DXB
  400    CONTINUE
C
 1000 CONTINUE
C
      RETURN
      END
C     ******* BAND MATRIX COEFFICIENT *******
C
      SUBROUTINE W1BNDQ(IERR)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1BND1/ CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
      COMMON /W1BND2/ CA(NXM)
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1MTRX/ CF(6*MATLM+5,NXM)
      PARAMETER (NCLW=2*MATLM+1)
      COMMON /W1QCLM/ CL(3,3,4,NCLM),NCL(NCLW,NXPM,ISM)
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
C
      DIMENSION DS0(2,2),DS1(2,2),DS2(2,2),DS3(2,2)
      DATA DS0/0.33333 33333 33333D0,0.16666 66666 66667D0,
     &         0.16666 66666 66667D0,0.33333 33333 33333D0/
      DATA DS1/1.D0,0.D0,0.D0,0.D0/
      DATA DS2/-1.D0,1.D0,0.D0,0.D0/
      DATA DS3/1.D0,-1.D0,-1.D0,1.D0/
      DATA CI/(0.D0,1.D0)/
C
      RW=2.D6*PI*RF
C
C
      NSF=3*NXP+4
      RKV=RW/VC
      DO 20 I=1,6*MATLM+5
      DO 20 N=1,NSF
         CF(I,N)=(0.D0,0.D0)
   20 CONTINUE
      DO 30 N=1,NSF
         CA(N)=(0.D0,0.0D0)
   30 CONTINUE
C
      KML=(MATLM-1)*3
      CF(KML+6,1)=CGIN(1,1)
      CF(KML+7,1)=CGIN(2,1)
      CF(KML+9,1)=(-1.D0,0.D0)
      CF(KML+5,2)=CGIN(1,3)
      CF(KML+6,2)=CGIN(2,3)
      CF(KML+9,2)=(-1.D0,0.D0)
      CF(KML+3,4)=CGIN(1,2)
      CF(KML+4,4)=CGIN(2,2)
      CF(KML+2,5)=CGIN(1,4)
      CF(KML+3,5)=CGIN(2,4)
C
      CF(KML+6,NSF-1)=CGOT(1,1)
      CF(KML+7,NSF-1)=CGOT(2,1)
      CF(KML+4,NSF-1)=(-1.D0,0.D0)
      CF(KML+5,NSF  )=CGOT(1,3)
      CF(KML+6,NSF  )=CGOT(2,3)
      CF(KML+4,NSF  )=(-1.D0,0.D0)
      CF(KML+8,NSF-3)=CGOT(1,2)
      CF(KML+9,NSF-3)=CGOT(2,2)
      CF(KML+7,NSF-2)=CGOT(1,4)
      CF(KML+8,NSF-2)=CGOT(2,4)
C
      CF(KML+3,NSF-4)=-1.D0
      CF(KML+6,NSF-4)= 1.D0
C
      CA(1    )=-CGIN(3,1)
      CA(2    )=-CGIN(3,3)
      CA(4    )=-CGIN(3,2)
      CA(5    )=-CGIN(3,4)
      CA(NSF-1)=-CGOT(3,1)
      CA(NSF  )=-CGOT(3,3)
      CA(NSF-3)=-CGOT(3,2)
      CA(NSF-2)=-CGOT(3,4)
C
      DO 6000 I=1,2
      DO 6000 J=1,2
         L=3*(J-I+MATLM+1)-1
         DTT0=DS0(I,J)
         DSS0=DS1(I,J)
         DTS1=DS2(I,J)
         DTS2=DS2(J,I)
         DTTW=DS3(I,J)
C
         DO 5000 NX=1,NXP-1
            RKPR=RKZ
            RNPR=VC*RKPR/RW
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            M=3*(NX+I-1)-1
            CF(L+1,M+1)=CF(L+1,M+1)
     &                 +(1.D0-RNPR*RNPR)*DSS0*DX
            CF(L+1,M+2)=CF(L+1,M+2)
     &                 +(1.D0-RNPR*RNPR)*DTT0*DX
     &                 -DTTW/DX
            CF(L+1,M+3)=CF(L+1,M+3)
     &                 +DTT0*DX
     &                 -DTTW/DX
            CF(L+3,M+1)=CF(L+3,M+1)
     &                 -CI*RNPR*DTS2
            CF(L-1,M+3)=CF(L-1,M+3)
     &                 +CI*RNPR*DTS1
 5000 CONTINUE
 6000 CONTINUE
C
      DO 7000 IS=1,ISMAX
      DO 7000 NX=1,NXP
         M=3*(NX-1)+2
C
      DO 7000 J=MAX(MATLM-NX+2,1)  ,MIN(2*MATLM+1,NXP+1-NX+MATLM)
         L=3*J-1
         ICL=NCL(J,NX,IS)
         IF(ICL.NE.0) THEN
           DO 40 IA=1,3
           DO 40 IB=1,3
             CF(L+IA-IB+1,M+IB  )=CF(L+IA-IB+1,M+IB  )
     &                            +CL(IB,IA,1,ICL)*RKV
             CF(L+IA-IB-2,M+IB+3)=CF(L+IA-IB-2,M+IB+3)
     &                            +CL(IB,IA,2,ICL)*RKV
             CF(L+IA-IB+4,M+IB  )=CF(L+IA-IB+4,M+IB  )
     &                            +CL(IB,IA,3,ICL)*RKV
             CF(L+IA-IB+1,M+IB+3)=CF(L+IA-IB+1,M+IB+3)
     &                            +CL(IB,IA,4,ICL)*RKV
   40      CONTINUE
         ENDIF
C
 7000 CONTINUE
C
C     WRITE(6,999) ((L,M,CF(L,M),L=30,41),M=195,197)
C 999 FORMAT((1H ,3(3X,I3,I4,1P2E14.6)))
C
      DMAT=MATLM-MATL
      DO 8000 IM=1,NSF
      DO 8000 IL=1,6*MATL+5
         CF(IL,IM)=CF(IL+3*DMAT,IM)
 8000 CONTINUE
C
      CALL BANDCD(CF,CA,NSF,6*MATL+5,6*MATLM+5,IND)
         IF(IND.NE.0) WRITE(6,601) IND
      IERR=0
      RETURN
C
  601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
      END
C
C     ******* ELECTROMAGNETIC FIELD IN PLASMA *******
C
      SUBROUTINE W1EPWQ(NZ)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1BND2/ CA(NXM)
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1PWR0/ PABS(NXPM,ISM),FLUX(NXTM)
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      PARAMETER (NCLW=2*MATLM+1)
      COMMON /W1QCLM/ CL(3,3,4,NCLM),NCL(NCLW,NXPM,ISM)
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
C
      DATA CI/(0.D0,1.D0)/
C
      RKV=2.D6*PI*RF/VC
C
      RCE=VC*EPS0
C
      DO 100 IS=1,ISMAX
      DO 100 NX=1,NXP
         PABS(NX,IS)=0.D0
  100 CONTINUE
      DO 110 NX=1,NXP
         FLUX(NX)=0.D0
  110 CONTINUE
C
      DO 200 NX=1,NXP
         CE2DA(NZ,NX,1)=CA(3*NX)
         CE2DA(NZ,NX,2)=CA(3*NX+1)
         CE2DA(NZ,NX,3)=CA(3*NX+2)
  200 CONTINUE
C
      DO 6000 IS=1,ISMAX
      DO 6000 NX=1,NXP-1
         DX=RKV*(XA(NX+1)-XA(NX))
      DO 6000 J=MAX(MATLM-NX+2,1),MIN(2*MATLM+1,NXP+1-NX+MATLM)
         MX=NX+J-MATLM-1
         ICL=NCL(J,NX,IS)
         IF(ICL.NE.0) THEN
            DO 6020 IL=1,4
               CABSL=0.D0
               NN=3*NX-1
               MM=3*MX-1
               IF(IL.EQ.2) NN=NN+3
               IF(IL.EQ.3) MM=MM+3
               IF(IL.EQ.4) NN=NN+3
               IF(IL.EQ.4) MM=MM+3
            DO 6010 IA=1,3
            DO 6010 IB=1,3
               CABSL=CABSL
     &              +CONJG(CA(NN+IA))*CL(IA,IB,IL,ICL)*CA(MM+IB)
 6010       CONTINUE
            PABSL=-CI*RCE*CABSL*RKV
            IF(IL.EQ.1) THEN
               PABS(NX,  IS)=PABS(NX,  IS)+0.5D0*PABSL
               PABS(MX,  IS)=PABS(MX,  IS)+0.5D0*PABSL
            ELSEIF(IL.EQ.2) THEN
               PABS(NX+1,IS)=PABS(NX+1,IS)+0.5D0*PABSL
               PABS(MX,  IS)=PABS(MX,  IS)+0.5D0*PABSL
            ELSEIF(IL.EQ.3) THEN
               PABS(NX,  IS)=PABS(NX,  IS)+0.5D0*PABSL
               PABS(MX+1,IS)=PABS(MX+1,IS)+0.5D0*PABSL
            ELSEIF(IL.EQ.4) THEN
               PABS(NX+1,IS)=PABS(NX+1,IS)+0.5D0*PABSL
               PABS(MX+1,IS)=PABS(MX+1,IS)+0.5D0*PABSL
            ENDIF
 6020       CONTINUE
         ENDIF
 6000 CONTINUE
C
      RETURN
      END
C
C
C      ****** MAKE TABLE OF F,G,H ******
C
      SUBROUTINE W1QTBL(NDMAX,XDMAX,NHARM)
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1QSFG/ SF(NDM,NHM1,3),SG(NDM,NHM1,3),
     &                AF(NDM,NHM1,2),AG(NDM,NHM1,2)
C
      DATA XDMAX1,NDMAX1,NHARM1/0.D0,0,0/
C
      IF(NDMAX.EQ.NDMAX1.AND.
     &   ABS(XDMAX-XDMAX1).LE.1.D-32.AND.
     &   NHARM.EQ.NHARM1) RETURN
      NDMAX1=NDMAX
      XDMAX1=XDMAX
      NHARM1=NHARM
      DXD=XDMAX/NDMAX
      DO 10 NC=1,NHARM+1
         NN=NC-1
      DO 10 N=1,NDMAX+1
         SF(N,NC,1)= W1FNMN((N-1)*DXD,1,NN,1)
         SF(N,NC,2)= W1FNMN((N-1)*DXD,1,NN,3)
         SF(N,NC,3)= W1FNMN((N-1)*DXD,1,NN,5)
         SG(N,NC,1)= W1FNMN((N-1)*DXD,2,NN,0)
         SG(N,NC,2)= W1FNMN((N-1)*DXD,2,NN,2)
         SG(N,NC,3)= W1FNMN((N-1)*DXD,2,NN,4)
         AF(N,NC,1)=-W1FNMN((N-1)*DXD,3,NN,1)
         AF(N,NC,2)=-W1FNMN((N-1)*DXD,3,NN,3)
         AG(N,NC,1)=-W1FNMN((N-1)*DXD,4,NN,0)
         AG(N,NC,2)=-W1FNMN((N-1)*DXD,4,NN,2)
   10 CONTINUE
      RETURN
C
      END
C
C     ****** FUNCTION FOR INTEGRO-DIFFERENTIAL ANALYSIS ******
C
      DOUBLE PRECISION FUNCTION W1FNMN(X,NF1,NN1,NM1)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /W1DEC1/ G1,NF,NN,NM
C
      DATA SP/2.506628275D0/
      DATA PI/3.141592654D0/
      EXTERNAL W1FNFG
C
      NF=NF1
      NN=NN1
      NM=NM1
      IF(ABS(X).LT.1.D-6) THEN
         IF(NF.EQ.0) THEN
            W1FNMN=16.D0/(3.D0*PI*SP)
         ELSEIF(NF.EQ.1) THEN
            IF(NM.EQ.1) W1FNMN=2.D0/(PI*SP*(1-4*NN*NN))
            IF(NM.EQ.3) W1FNMN=3.D0/(2.D0*PI*SP)
     &                        *(1.D0/(1-4*NN*NN)
     &                         -1.D0/(9-4*NN*NN))
            IF(NM.EQ.5) W1FNMN=5.D0/(8.D0*PI*SP)
     &                        *(1.D0/(25-4*NN*NN)
     &                         -3.D0/( 9-4*NN*NN)
     &                         +2.D0/( 1-4*NN*NN))
         ELSEIF(NF.EQ.2) THEN
            IF(NM.EQ.0) W1FNMN=4.D0*NN/(PI*SP*(4.D0*NN*NN-1.D0))
            IF(NM.EQ.2) W1FNMN=1.D0/(2.D0*PI*SP)
     &                        *( 4.D0*NN   /(4*NN*NN-1)
     &                         -(2.D0*NN+1)/((2*NN+1)**2-4)
     &                         -(2.D0*NN-1)/((2*NN-1)**2-4))
            IF(NM.EQ.4) W1FNMN=1.D0/(8.D0*PI*SP)
     &                        *((2.D0*NN+1)/((2*NN+1)**2-16)
     &                         +(2.D0*NN-1)/((2*NN-1)**2-16)
     &                         -4*(2.D0*NN+1)/((2*NN+1)**2-4)
     &                         -4*(2.D0*NN-1)/((2*NN-1)**2-4)
     &                         +3.D0/(2*NN+1)
     &                         +3.D0/(2*NN-1))
         ELSEIF(NF.EQ.3) THEN
            W1FNMN=0.D0
            IF(NM.EQ.1) THEN
               IF(NN.EQ.0) W1FNMN= 0.5D0
               IF(NN.EQ.1) W1FNMN=-0.25D0
            ELSEIF(NM.EQ.3) THEN
               IF(NN.EQ.0) W1FNMN= 0.375D0
               IF(NN.EQ.1) W1FNMN=-0.25D0
               IF(NN.EQ.2) W1FNMN= 0.0625D0
            ELSEIF(NM.EQ.5) THEN
               IF(NN.EQ.0) W1FNMN= 5.D0/16.D0
               IF(NN.EQ.1) W1FNMN=-15.D0/64.D0
               IF(NN.EQ.2) W1FNMN= 3.D0/32.D0
               IF(NN.EQ.3) W1FNMN=-1.D0/64.D0
            ENDIF
         ELSEIF(NF.EQ.4) THEN
            W1FNMN=0.D0
            IF(NM.EQ.0) THEN
               IF(NN.EQ.1) W1FNMN= 0.25D0
            ELSEIF(NM.EQ.2) THEN
               IF(NN.EQ.1) W1FNMN= 0.125D0
               IF(NN.EQ.2) W1FNMN=-0.0625D0
            ELSEIF(NM.EQ.4) THEN
               IF(NN.EQ.1) W1FNMN= 5.D0/64.D0
               IF(NN.EQ.2) W1FNMN=-1.D0/16.D0
               IF(NN.EQ.3) W1FNMN= 1.D0/64.D0
            ENDIF
         ENDIF
      ELSE
         H0=0.5
         EPS=1.D-6
         ILST=0
         G1=X
         CALL DEFT(W1FNFG,CS,ES,H0,EPS,ILST)
         W1FNMN=CS
      ENDIF
      RETURN
      END
C
C     ****** SLAVE FUNTION FOR DE INTEGRATION ******
C
      DOUBLE PRECISION FUNCTION W1FNFG(X,XM,XP)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /W1DEC1/ G1,NF,NN,NM
C
      DATA SP/2.506628275D0/
      DATA PI/3.141592654D0/
C
      DUMMY=X
      DUMMY=XM
      T=0.5D0*PI*XP
      S=SIN(T)
      XX=G1/(2.82843D0*S)
      IF(ABS(XX).LT.10.D0) THEN
         IF(NF.EQ.0) THEN
            W1FNFG = 2.D0*S*S*S*(-XX*ERFC1(XX)+EXP(-XX*XX)/SP)
         ELSEIF(NF.EQ.1) THEN
            W1FNFG = 0.5D0*S**NM*EXP(-XX*XX)*COS(2*NN*T)/SP
         ELSEIF(NF.EQ.2) THEN
            W1FNFG = 0.5D0*COS(T)*S**NM*EXP(-XX*XX)*SIN(2*NN*T)/SP
         ELSEIF(NF.EQ.3) THEN
            W1FNFG = 0.5D0*S**(NM+1)*ERFC1(XX)*COS(2*NN*T)
         ELSEIF(NF.EQ.4) THEN
            W1FNFG = 0.5D0*COS(T)*S**(NM+1)*ERFC1(XX)*SIN(2*NN*T)
         ENDIF
      ELSE
         W1FNFG = 0.D0
      ENDIF
      END
C
C     ******** LINEAR INTERPOLATION  ******************************
C
      SUBROUTINE W1QCAL(GX,DXA,DXB,YX,CT0,CT1,CT2,CT3,CSB,NN)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
      COMMON /W1QSFG/ SF(NDM,NHM1,3),SG(NDM,NHM1,3),
     &                AF(NDM,NHM1,2),AG(NDM,NHM1,2)
C
      DIMENSION  GX(4),CSB(3,3,4),AF0(4,3),AG0(4,3),
     &           AF1(4,2),AF2(4,2),AF3(4,1),AF4(4,1),
     &           AG1(4,2),AG2(4,2),AG3(4,1),AG4(4,1)
      DIMENSION  GV(4),AAF(4),AAG(4),ADF(4),ADG(4),AQF(4),AQG(4),AQH(4)
      DIMENSION  ABF(4),ACF(4),ABG(4),BDG(4)
C
      NNA=ABS(NN)+1
      ESGN=SIGN(1.D0,DBLE(NN))
C
      DO 10 IV=1,4
      GV(IV)=YX*GX(IV)
      VSGN=SIGN(1.D0,GV(IV))
      YSGN=SIGN(1.D0,YX)
      V=ABS(GV(IV))
      IF(V.GT.XDMAX) THEN
         DO 7 IR=1,3
            AF0(IV,IR)=0.D0
            AG0(IV,IR)=0.D0
    7    CONTINUE
         AF1(IV,1)=0.D0
         AG1(IV,1)=0.D0
         AF1(IV,2)=0.D0
         AG1(IV,2)=0.D0
      ELSEIF(ABS(GX(IV)).LT.1.D-70) THEN
         AF0(IV,1)= SF(1,NNA,1)
         AF0(IV,2)= SF(1,NNA,2)
         AF0(IV,3)= SF(1,NNA,3)
         AG0(IV,1)= SG(1,NNA,1)*ESGN
         AG0(IV,2)= SG(1,NNA,2)*ESGN
         AG0(IV,3)= SG(1,NNA,3)*ESGN
         IF(IV.EQ.1.OR.IV.EQ.3) THEN
            AF1(IV,1)= AF(1,NNA,1)*YSGN
            AF1(IV,2)= AF(1,NNA,2)*YSGN
            AG1(IV,1)= AG(1,NNA,1)*YSGN*ESGN
            AG1(IV,2)= AG(1,NNA,2)*YSGN*ESGN
         ELSE
            AF1(IV,1)=-AF(1,NNA,1)*YSGN
            AF1(IV,2)=-AF(1,NNA,2)*YSGN
            AG1(IV,1)=-AG(1,NNA,1)*YSGN*ESGN
            AG1(IV,2)=-AG(1,NNA,2)*YSGN*ESGN
         ENDIF
      ELSE
         PV=V/DXD
         NPV=INT(PV)
         IF(NPV.GE.NDMAX) THEN
            NPV1=NDMAX+1
            NPV2=NDMAX+1
         ELSE
            NPV1=NPV+1
            NPV2=NPV+2
         ENDIF
         PVS=PV-DBLE(NPV)
         PVT=1.D0-PVS
         AF0(IV,1)= PVT*SF(NPV1,NNA,1)+PVS*SF(NPV2,NNA,1)
         AF0(IV,2)= PVT*SF(NPV1,NNA,2)+PVS*SF(NPV2,NNA,2)
         AF0(IV,3)= PVT*SF(NPV1,NNA,3)+PVS*SF(NPV2,NNA,3)
         AG0(IV,1)=(PVT*SG(NPV1,NNA,1)+PVS*SG(NPV2,NNA,1))*ESGN
         AG0(IV,2)=(PVT*SG(NPV1,NNA,2)+PVS*SG(NPV2,NNA,2))*ESGN
         AG0(IV,3)=(PVT*SG(NPV1,NNA,3)+PVS*SG(NPV2,NNA,3))*ESGN
         AF1(IV,1)=(PVT*AF(NPV1,NNA,1)+PVS*AF(NPV2,NNA,1))*VSGN
         AF1(IV,2)=(PVT*AF(NPV1,NNA,2)+PVS*AF(NPV2,NNA,2))*VSGN
         AG1(IV,1)=(PVT*AG(NPV1,NNA,1)+PVS*AG(NPV2,NNA,1))*VSGN*ESGN
         AG1(IV,2)=(PVT*AG(NPV1,NNA,2)+PVS*AG(NPV2,NNA,2))*VSGN*ESGN
      ENDIF
      AF2(IV,1)= 4.D0*AF0(IV,2)+GV(IV)*AF1(IV,1)
      AF2(IV,2)= 4.D0*AF0(IV,3)+GV(IV)*AF1(IV,2)
      AF3(IV,1)=(4.D0*AF1(IV,2)+GV(IV)*AF2(IV,1))/2.D0
      AF4(IV,1)=(4.D0*AF2(IV,2)+GV(IV)*AF3(IV,1))/3.D0
      AG2(IV,1)= 4.D0*AG0(IV,2)+GV(IV)*AG1(IV,1)
      AG2(IV,2)= 4.D0*AG0(IV,3)+GV(IV)*AG1(IV,2)
      AG3(IV,1)=(4.D0*AG1(IV,2)+GV(IV)*AG2(IV,1))/2.D0
      AG4(IV,1)=(4.D0*AG2(IV,2)+GV(IV)*AG3(IV,1))/3.D0
   10 CONTINUE
C
      EI=YX*DXA
      EJ=YX*DXB
      E0=EI
      E1=EI*EJ
      E2=EI*EI*EJ
      E3=EI*EJ*EJ
      E4=EI*EI*EJ*EJ
C
      FA0=AF0(1,1)-AF0(2,1)-AF0(3,1)+AF0(4,1)
      FA1=AF1(1,1)-AF1(2,1)-AF1(3,1)+AF1(4,1)
      FA2=AF2(1,1)-AF2(2,1)-AF2(3,1)+AF2(4,1)
      FA3=AF3(1,1)-AF3(2,1)-AF3(3,1)+AF3(4,1)
      FA4=AF4(1,1)-AF4(2,1)-AF4(3,1)+AF4(4,1)
      GA0=AG0(1,1)-AG0(2,1)-AG0(3,1)+AG0(4,1)
      GA1=AG1(1,1)-AG1(2,1)-AG1(3,1)+AG1(4,1)
      GA2=AG2(1,1)-AG2(2,1)-AG2(3,1)+AG2(4,1)
      GA3=AG3(1,1)-AG3(2,1)-AG3(3,1)+AG3(4,1)
      GA4=AG4(1,1)-AG4(2,1)-AG4(3,1)+AG4(4,1)
      HA0=AF0(1,2)-AF0(2,2)-AF0(3,2)+AF0(4,2)
      HA1=AF1(1,2)-AF1(2,2)-AF1(3,2)+AF1(4,2)
      HA2=AF2(1,2)-AF2(2,2)-AF2(3,2)+AF2(4,2)
C
      FF  =-FA2/E1
      FFS = FA3/E3-(AF2(4,1)-AF2(2,1))/E1
      FFT =-FA3/E2-(AF2(4,1)-AF2(3,1))/E1
      FFST= FA4/E4-(AF3(4,1)-AF3(2,1))/E3
     &            +(AF3(4,1)-AF3(3,1))/E2-AF2(4,1)/E1
      FG  =-GA2/E1
      FGS = GA3/E3-(AG2(4,1)-AG2(2,1))/E1
      FGT =-GA3/E2-(AG2(4,1)-AG2(3,1))/E1
      FGST= GA4/E4-(AG3(4,1)-AG3(2,1))/E3
     &            +(AG3(4,1)-AG3(3,1))/E2-AG2(4,1)/E1
      DF  =-FA1/E1
      DFS = FA2/E3-(AF1(4,1)-AF1(2,1))/E1
      DFT =-FA2/E2-(AF1(4,1)-AF1(3,1))/E1
      DFST= FA3/E4-(AF2(4,1)-AF2(2,1))/E3
     &            +(AF2(4,1)-AF2(3,1))/E2-AF1(4,1)/E1
      DG  =-GA1/E1
      DGS = GA2/E3-(AG1(4,1)-AG1(2,1))/E1
      DGT =-GA2/E2-(AG1(4,1)-AG1(3,1))/E1
      DGST= GA3/E4-(AG2(4,1)-AG2(2,1))/E3
     &            +(AG2(4,1)-AG2(3,1))/E2-AG1(4,1)/E1
      QF  =-FA0/E1
      QFS = FA1/E3-(AF0(4,1)-AF0(2,1))/E1
      QFT =-FA1/E2-(AF0(4,1)-AF0(3,1))/E1
      QFST= FA2/E4-(AF1(4,1)-AF1(2,1))/E3
     &            +(AF1(4,1)-AF1(3,1))/E2-AF0(4,1)/E1
      QG  =-GA0/E1
      QGS = GA1/E3-(AG0(4,1)-AG0(2,1))/E1
      QGT =-GA1/E2-(AG0(4,1)-AG0(3,1))/E1
      QGST= GA2/E4-(AG1(4,1)-AG1(2,1))/E3
     &            +(AG1(4,1)-AG1(3,1))/E2-AG0(4,1)/E1
      QH  =-HA0/E1
      QHS = HA1/E3-(AF0(4,2)-AF0(2,2))/E1
      QHT =-HA1/E2-(AF0(4,2)-AF0(3,2))/E1
      QHST= HA2/E4-(AF1(4,2)-AF1(2,2))/E3
     &            +(AF1(4,2)-AF1(3,2))/E2-AF0(4,2)/E1
C
      IF(ABS(GX(1)).LT.1.D-70) THEN
         FF  =FF  -(AF1(1,1)-AF1(4,1))/E0
         FFS =FFS -(AF1(1,1)-AF1(4,1))/(2.D0*E0)
         FFT =FFT -(AF1(1,1)-AF1(4,1))/(2.D0*E0)
     &            -(AF2(1,1)-AF2(4,1))/E1
         FFST=FFST-(AF1(1,1)-AF1(4,1))/(3.D0*E0)
     &            -(AF2(1,1)-AF2(4,1))/(2.D0*E1)
         FG  =FG  -(AG1(1,1)-AG1(4,1))/E0
         FGS =FGS -(AG1(1,1)-AG1(4,1))/(2.D0*E0)
         FGT =FGT -(AG1(1,1)-AG1(4,1))/(2.D0*E0)
     &            -(AG2(1,1)-AG2(4,1))/E1
         FGST=FGST-(AG1(1,1)-AG1(4,1))/(3.D0*E0)
     &            -(AG2(1,1)-AG2(4,1))/(2.D0*E1)
         DFT =DFT -(AF1(1,1)-AF1(4,1))/E1
         DFST=DFST-(AF1(1,1)-AF1(4,1))/(2.D0*E1)
         DGT =DGT -(AG1(1,1)-AG1(4,1))/E1
         DGST=DGST-(AG1(1,1)-AG1(4,1))/(2.D0*E1)
      ENDIF
C
      AAF(1)=FF-FFS-FFT+FFST
      AAF(2)=FFS-FFST
      AAF(3)=FFT-FFST
      AAF(4)=FFST
      ABF(1)=FF-FFS
      ABF(2)=FFS
      ABF(3)=0.D0
      ABF(4)=0.D0
      ACF(1)=FF-FFT
      ACF(2)=0.D0
      ACF(3)=FFT
      ACF(4)=0.D0
      AAG(1)=FG-FGS-FGT+FGST
      AAG(2)=FGS-FGST
      AAG(3)=FGT-FGST
      AAG(4)=FGST
      ABG(1)=FG
      ABG(2)=0.D0
      ABG(3)=0.D0
      ABG(4)=0.D0
      ADF(1)=DF-DFS-DFT+DFST
      ADF(2)=DFS-DFST
      ADF(3)=DFT-DFST
      ADF(4)=DFST
      ADG(1)=DG-DGS
      ADG(2)=DGS
      ADG(3)=0.D0
      ADG(4)=0.D0
      BDG(1)=DG-DGT
      BDG(2)=0.D0
      BDG(3)=DGT
      BDG(4)=0.D0
      AQF(1)=QF-QFS-QFT+QFST
      AQF(2)=QFS-QFST
      AQF(3)=QFT-QFST
      AQF(4)=QFST
      AQG(1)=QG-QGS-QGT+QGST
      AQG(2)=QGS-QGST
      AQG(3)=QGT-QGST
      AQG(4)=QGST
      AQH(1)=QH-QHS-QHT+QHST
      AQH(2)=QHS-QHST
      AQH(3)=QHT-QHST
      AQH(4)=QHST
C
      DO 25 IU=1,4
        CSB(1,1,IU)=NN*ABG(IU)*CT1
        CSB(1,2,IU)=(0.D0,-1.D0)*NN*ACF(IU)*CT1
        CSB(1,3,IU)= BDG(IU)*CT2
        CSB(2,1,IU)=(0.D0,1.D0)*NN*ABF(IU)*CT1
        CSB(2,2,IU)=(NN*AAG(IU)-2.D0*AQF(IU))*CT1
        CSB(2,3,IU)=-(0.D0,1.D0)*ADF(IU)*CT2
        CSB(3,1,IU)= ADG(IU)*CT2
        CSB(3,2,IU)= (0.D0,1.D0)*ADF(IU)*CT2
        CSB(3,3,IU)=-CT0
     &              *(AAF(IU)+NN*AAG(IU)-2.D0*AQF(IU)+2.D0*AQH(IU))
     &              -CT3*AQG(IU)
   25 CONTINUE
C
      RETURN
      END
C
C     ******** LINEAR INTERPOLATION 2 ******************************
C
      SUBROUTINE W1QLNI(V,VX,VK,CT0,CT1,CT2,CT3)
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1QSST/ XM(NXPM),YX(NXPM),YK(NXPM),
     &                CS0(NXPM),CS1(NXPM),CS2(NXPM),CS3(NXPM)
C
      IF(V.LE.XM(1)) THEN
         VX=YX(1)
         VK=YK(1)
         CT0=CS0(1)
         CT1=CS1(1)
         CT2=CS2(1)
         CT3=CS3(1)
      ELSEIF(V.GE.XM(NXP)) THEN
         VX=YX(NXP)
         VK=YK(NXP)
         CT0=CS0(NXP)
         CT1=CS1(NXP)
         CT2=CS2(NXP)
         CT3=CS3(NXP)
      ELSE
         NHF=INT(NXP*(V-XM(1))/(XM(NXP)-XM(1)))
         IF(NHF.LT.1)   NHF=1
         IF(NHF.GT.NXP) NHF=NXP
         IF(V.LT.XM(NHF)) THEN
            DO 10 I=NHF-1,1,-1
               IF(XM(I).LT.V) THEN
                  NPV=I
                  GOTO 100
               ENDIF
   10       CONTINUE
         ELSE
            DO 20 I=NHF+1,NXP
               IF(XM(I).GT.V) THEN
                  NPV=I-1
                  GOTO 100
               ENDIF
   20       CONTINUE
         ENDIF
  100    PVS=(V-XM(NPV))/(XM(NPV+1)-XM(NPV))
         VX  =(1.D0-PVS)*YX(NPV) +PVS*YX(NPV+1)
         VK  =(1.D0-PVS)*YK(NPV) +PVS*YK(NPV+1)
         CT0 =(1.D0-PVS)*CS0(NPV)+PVS*CS0(NPV+1)
         CT1 =(1.D0-PVS)*CS1(NPV)+PVS*CS1(NPV+1)
         CT2 =(1.D0-PVS)*CS2(NPV)+PVS*CS2(NPV+1)
         CT3 =(1.D0-PVS)*CS3(NPV)+PVS*CS3(NPV+1)
      ENDIF
      RETURN
      END
