C
C     ****** TASK/W1 ******
C
C     ******* DIVISION IN X DIRECTION *******
C
      SUBROUTINE W1SETX(IERR)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1CTRL/ DRF,DRKZ,DXFACT,DXWDTH,APRFPN,APRFTR,APRFTP,
     &                NPRINT,NFILE,NGRAPH,NLOOP,NSYM,NFLR,NALPHA,
     &                NSYS,NDISP,NXABS
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
C
      NXT=NXP+2*NXV
         IF(NXP.GT.NXPM) GO TO 9100
         IF(NXV.GT.NXVM) GO TO 9100
      NXVH=NXV/2
C
      IF(MOD(NMODEL,2).EQ.1) THEN
         NXP1=NXP
      ELSE
         NXP1=NXP-1
      ENDIF
C
      RW=2.D6*PI*RF
      ISACM=0
      DO 10  IS=1,ISMAX
         FWC=AEE*PZ(IS)*BB/(AMM*PA(IS))
         WCH=FWC/(1.D0+RA/RR)
         WCL=FWC/(1.D0-RA/RR)
         IF(ABS(WCH).LT.RW.AND.ABS(WCL).GT.RW) THEN
            XACM=RR*(FWC/RW-1.D0)
            ISACM=IS
         ENDIF
   10 CONTINUE
C
      IF(ISACM.EQ.0.OR.DXFACT.LT.1.D-6) THEN
         DX=2.D0*RA/DBLE(NXP1)
         DO 20 NX=1,NXP1+1
            XA(NX)=DBLE(NX-1)*DX-RA
   20    CONTINUE
      ELSE
         FVT=AEE*1.D3/(AMM*PA(ISACM))
         FWC=AEE*PZ(ISACM)*BB/(AMM*PA(ISACM))
         DXACM=DXWDTH*SQRT(FVT*PTPP(ISACM))/FWC
         WT=2.D0*RA+DXFACT*DXACM
     &             *(ATAN((RA-XACM)/DXACM)-ATAN((-RA-XACM)/DXACM))
         X=-RA
         DO 25 NX=1,NXP1
            XA(NX)=X
            FACT=1.D0+DXFACT/(1.D0+((X-XACM)/DXACM)**2)
            X1=X+0.5D0*WT/(FACT*(NXP1))
            FACT=1.D0+DXFACT/(1.D0+((X1-XACM)/DXACM)**2)
            X=X+WT/(FACT*(NXP1))
   25    CONTINUE
         WRITE(6,601) ISACM,XACM,DXACM,RA-X
         XA(NXP1+1)=RA
      ENDIF
C
      IF(MOD(NMODEL,2).EQ.1) THEN
         DO 27 NX=1,NXP
            XAM(NX)=0.5D0*(XA(NX)+XA(NX+1))
   27    CONTINUE
      ELSE
         DO 28 NX=1,NXP
            XAM(NX)=XA(NX)
   28    CONTINUE
      ENDIF
C
      DX1=(RB-RD)/DBLE(NXVH-1)
      DX2=(RD-RA)/DBLE(NXVH-1)
      DO 30 I=1,NXV
         IF(I.LE.NXVH) THEN
            X=DX1*DBLE(I-1)-RB
         ELSE
            X=DX2*DBLE(I-NXVH-1)-RD
         ENDIF
         NX1=NXP+NXV+I
         XAM(NX1)=X
         NX2=NXP+NXV+1-I
         XAM(NX2)=-X
   30 CONTINUE
      IERR=0
      RETURN
C
 9100 WRITE(6,691) NXP,NXPM,NXV,NXVM
      IERR=10000
      RETURN
  601 FORMAT(1H ,'** DX-ACCUM. : IS,X,XWIDTH,XERR = ',
     &           I3,1P2E12.4,1PE12.4)
  691 FORMAT(1H ,'!! ERROR IN W1SETX   : DIMENSION OVER'/
     &       1H ,'NXP=',I5,5X,'NXPM=',I5/
     &       1H ,'NXV=',I5,5X,'NXVM=',I5/)
      END
C
C     ******* DIVISION IN Z DIRECTION *******
C
      SUBROUTINE W1SETZ(IERR)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1ZDAT/ CJ1(NZPM),CJ2(NZPM),CJ3(NZPM),CJ4(NZPM),
     &                ZA(NZPM),AKZ(NZPM),NZANT1(IAM),NZANT2(IAM),NANT
C
         IF(NZP.GT.NZPM) GO TO 9100
      NZH=NZP/2
C
      DZ=RZ/DBLE(NZP)
      IF(NZP.EQ.1) THEN
         ZA(1)=0.D0
         AKZ(1)=RKZ
         IF(ABS(RKZ).LT.1.D-32) RKZ=0.001D0
      ELSE
         DO 40 NZ=1,NZP
            ZA(NZ)=DZ*DBLE(NZ-1) - 0.5D0*RZ
   40    CONTINUE
         AKZ0=2.D0*PI/RZ
         AKZ(1)   =0.001D0
         AKZ(NZH+1)=DBLE(NZH)*AKZ0
         DO 50 NZ=2,NZH
            AKZ(NZ)      =DBLE(NZ-1)*AKZ0
            AKZ(NZP-NZ+2)=-AKZ(NZ)
   50    CONTINUE
      ENDIF
C
      IERR=0
      RETURN
C
 9100 WRITE(6,691) NZP,NZPM
      IERR=10000
      RETURN
  691 FORMAT(1H ,'!! ERROR IN W1SETZ   : DIMENSION OVER'/
     &       1H ,'NZP=',I5,5X,'NZPM=',I5/)
      END
C
C     ******* SETING ANTENNA CURRENT DENSITY *******
C
      SUBROUTINE W1ANTS
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1ANT1/ AJYL(IAM),AJYH(IAM),AJZL(IAM),AJZH(IAM),
     &                ALYL(IAM),ALYH(IAM),APYL(IAM),APYH(IAM)
      COMMON /W1ZDAT/ CJ1(NZPM),CJ2(NZPM),CJ3(NZPM),CJ4(NZPM),
     &                ZA(NZPM),AKZ(NZPM),NZANT1(IAM),NZANT2(IAM),NANT
C
      DO 10 NZ = 1 , NZP
         CJ1( NZ ) = ( 0.D0 , 0.D0 )
         CJ2( NZ ) = ( 0.D0 , 0.D0 )
         CJ3( NZ ) = ( 0.D0 , 0.D0 )
         CJ4( NZ ) = ( 0.D0 , 0.D0 )
   10 CONTINUE
C
      IF(NZP.EQ.1) THEN
         NZANT1(1)=1
         CJ1(1)   =AJYH(1)/DZ
         NZANT2(1)=1
         CJ2(1)   =AJYL(1)/DZ
         CJ3(1)   =AJZH(1)
         CJ4(1)   =AJZL(1)
      ELSE
      DO 100 IA = 1 , IAMAX
         ANTPOS            = RZ*ALYH(IA)/360.D0 + (RZ + DZ)*.5D0
         NZANT1(IA)        = INT( ANTPOS/DZ ) + 1
         PHASE             = APYH(IA)*PI/180.D0
         CPHASE            = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
         CJ1( NZANT1(IA) ) = AJYH( IA ) / DZ * CPHASE
     &                       + CJ1( NZANT1(IA) )
         CJ3( NZANT1(IA) ) = AJZH( IA ) * CPHASE
     &                       + CJ3( NZANT1(IA) )
C
         ANTPOS            = RZ*ALYL(IA)/360.D0 + (RZ + DZ)*.5D0
         NZANT2(IA)        = INT( ANTPOS/DZ ) + 1
         PHASE             = APYL(IA)*PI/180.D0
         CPHASE            = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
         CJ2( NZANT2(IA) ) = AJYL( IA ) / DZ * CPHASE
     &                       + CJ2( NZANT2(IA) )
         CJ4( NZANT2(IA) ) = AJZL( IA ) * CPHASE
     &                       + CJ4( NZANT2(IA) )
  100 CONTINUE
      ENDIF
C
      RETURN
      END
C
C     ****** INTERFACE FOR FFT ******
C
      SUBROUTINE W1FFTL(CA,N,KEY)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1FFTW/ CT(NZPM),CFFT(NZPM*NZLM),LFFT(NZPM*NZLM)
C
      COMPLEX*16 CA(N)
      DATA NS/0/
C
      IF(N.GT.1) THEN
         LP=LOG(DBLE(N))/LOG(2.D0)+0.5D0
         IF(N.EQ.NS) THEN
            IND=0
         ELSE
            IND=1
            NS=N
         ENDIF
         DO 10 I=1,N
           CT(I)=CA(I)
   10    CONTINUE
         CALL FFT2L(CT,CA,CFFT,LFFT,N/2,IND,KEY+1,LP)
         IF(KEY.EQ.0) CA(N/2+1)=0.D0
      ENDIF
      RETURN
      END
C
C     ****** SUBSTITUTION BY VIRUE OF SYMMETRY IN Z-DIRECTION ******
C
      SUBROUTINE W1SYMS(NZ,NSYM)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1PWR2/ PAKT(NZPM,3),PANTK(NZPM),PAK(NZPM,ISM)
C
      NZ0 = NZP +2 -NZ
C
      IF(NSYM.EQ.1) THEN
         DO 70 NX = 1, NXT
            CE2DA(NZ,NX,1) =  CE2DA(NZ0,NX,1)
            CE2DA(NZ,NX,2) =  CE2DA(NZ0,NX,2)
            CE2DA(NZ,NX,3) = -CE2DA(NZ0,NX,3)
   70    CONTINUE
      ELSEIF(NSYM.EQ.-1) THEN
         DO 80 NX = 1, NXT
            CE2DA(NZ,NX,1) = -CE2DA(NZ0,NX,1)
            CE2DA(NZ,NX,2) = -CE2DA(NZ0,NX,2)
            CE2DA(NZ,NX,3) =  CE2DA(NZ0,NX,3)
   80    CONTINUE
      ENDIF
      PAKT(NZ,1) =  PAKT(NZ0,1)
      PAKT(NZ,2) =  PAKT(NZ0,2)
      PAKT(NZ,3) =  PAKT(NZ0,3)
      PANTK(NZ)  =  PANTK(NZ0)
      DO 90 IS = 1, ISMAX
         PAK(NZ,IS) =  PAK(NZ0,IS)
   90 CONTINUE
C
      RETURN
      END
