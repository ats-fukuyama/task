C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 2) ******
C
C     ******* LOCAL WAVE NUMBER AND POLARIZATION *******
C
      SUBROUTINE W1WKXB
C
      USE w1comm
      IMPLICIT NONE
      INTEGER:: NX,N,NS,KK
      REAL(rkind):: RKV,DX
      COMPLEX(rkind):: CDSP2,CDSP1,CDSP0,CDET,CK
      COMPLEX(rkind):: CDTXX,CDTXY,CDTXZ,CDTYY,CDTYZ,CDTZZ,CDETIP
      COMPLEX(rkind):: CSOL(2,NXPM)
      REAL(rkind),PARAMETER:: EXPARG=80.D0
C
      RKV=2.D6*PI*RF/VC
C
      DO 100 NX = 1 , NXP
         CDSP2 = CD0(1,NX)+CD1(1,NX)**2
         CDSP1 =-CD0(1,NX)*(CD0(3,NX)+CD0(4,NX))-CD0(2,NX)**2
     &          +2.D0*CD0(2,NX)*CD1(1,NX)*CD1(2,NX)
     &          +CD0(1,NX)*CD1(2,NX)**2
     &          -CD0(3,NX)*CD1(1,NX)**2
         CDSP0 = CD0(4,NX)*(CD0(1,NX)*CD0(3,NX)+CD0(2,NX)**2)
C
         CDET  = SQRT(CDSP1**2-4.D0*CDSP2*CDSP0)
         CSOL(1,NX) = (-CDSP1+CDET)/(2.D0*CDSP2)
         CSOL(2,NX) = (-CDSP1-CDET)/(2.D0*CDSP2)
  100 CONTINUE
C
      DO 300 N = 1 , 2
         DO 200 NX = 1 , NXP
            NS = (NX-1)*4 + (N-1)*2 + 1
            KK = (NX-1)*2 + N
            DX = RKV*(XA(NX+1)-XA(NX))
            CK = SQRT(CSOL(N,NX))
            IF(ABS(IMAG(CK)*DX).GT.EXPARG) THEN
               IF(IMAG(CK).GE.0.D0) THEN
                  CK = DCMPLX(DBLE(CK), EXPARG/DX)
               ELSE
                  CK = DCMPLX(DBLE(CK),-EXPARG/DX)
               ENDIF
            ENDIF
            CSKX ( KK ) =   CK
C
            CDTXX = CD0(1,NX) + CSOL(N,NX)*CD2(1,NX)
            CDTXY = CD0(2,NX) + CSOL(N,NX)*CD2(2,NX)
            CDTXZ = CD1(1,NX) * CSKX( KK )
            CDTYY = CD0(3,NX) + CSOL(N,NX)*CD2(3,NX)
            CDTYZ = CD1(2,NX) * CSKX( KK )
            CDTZZ = CD0(4,NX) + CSOL(N,NX)*CD2(4,NX)
C
            CDETIP = 1.D0 / ( CDTXX*CDTYZ + CDTXZ*CDTXY )
            CSPX ( KK ) = (  CDTYY*CDTXZ - CDTXY*CDTYZ ) * CDETIP
            CSPZ ( KK ) =-(  CDTXY*CDTXY + CDTYY*CDTXX ) * CDETIP
C
            CF( 1 , NS+10 ) =  ( 1.D0 , 0.D0 )
            CF( 1 , NS+11 ) =  ( 1.D0 , 0.D0 )
            CF( 2 , NS+10 ) =  (  CD2(2,NX)*CSPX( KK )
     &                          - CD2(3,NX)          )*CSKX( KK )
     &                          - CD1(2,NX)*CSPZ( KK )
            CF( 2 , NS+11 ) = -(  CD2(2,NX)*CSPX( KK )
     &                          - CD2(3,NX)          )*CSKX( KK )
     &                          + CD1(2,NX)*CSPZ( KK )
            CF( 3 , NS+10 ) =   CSPZ( KK )
            CF( 3 , NS+11 ) = - CSPZ( KK )
            CF( 4 , NS+10 ) =   CD1(1,NX)*CSPX( KK )
     &                        + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
            CF( 4 , NS+11 ) =   CD1(1,NX)*CSPX( KK )
     &                        + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
  200    CONTINUE
  300 CONTINUE
      RETURN
      END
C
C     ******* BAND MATRIX COEFFICIENT *******
C
      SUBROUTINE W1BNDB(IERR,NXABSL)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NXABSL
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: NS,NW,NWH,NCF,NSF,NSE,NT1,NF,I,J,II,JJ
      INTEGER:: NS2,NS0,NS1,IMODE,FLUXB,IND,NX
      REAL(rkind):: RKV,X2,X1,DX
      COMPLEX(rkind):: CSKXB,CSPXB,CSPZB,CARG,CPHASE
C
      RKV=2.D6*PI*RF/VC
C
      NS  = 10
      NW  = 3*4 - 1
      NWH = 4 + 4/2
C
      NCF=0
      NSF=0
      NSE=2
      NT1=4
C
      NF=NSF-NCF+NWH
      DO 10 I=2-NF,0
         DO 10 J=1,MIN(NF+I-1,NT1)
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=(0.D0,0.D0)
   10 CONTINUE
      DO 20 I=1,NSE
         DO 20 J=1,NT1
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=-CGIN(I,J)
   20 CONTINUE
      NSF=NSF+NSE
C
      X2  = XA ( 1 )
      NS2 = 4
C
      DO 6000 NX=1,NXP
         X1=X2
         X2=XA(NX+1)
         DX  = RKV * ( X2-X1 )
         NS0 = NT1
         NS1 = NS2
         NS2 = 4
         NSE = NS1
         NT1 = NS1
         IF( NX .EQ. NXP ) NT1 = 4
C
         NF=NSF-NCF+NWH
         DO 5010 I=1,NSE
            DO 5010 J=1,NS0
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=CF(J,NS+I)
 5010    CONTINUE
         DO 5020 I=NSE+1,NW+NS0-NF
            DO 5020 J=MAX(I-NW+NF,1),NS0
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=(0.D0,0.D0)
 5020    CONTINUE
C
         IF(NX.EQ.ABS(NXABSL)) THEN
            DO 5022 I=1,NSE
               IMODE=(I-1)/2+1
               CSKXB=CSKX((NX-1)*2+IMODE)
               CSPXB=CSPX((NX-1)*2+IMODE)
               CSPZB=CSPZ((NX-1)*2+IMODE)
               IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
               IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
               FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB
     &               -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB 
     &                      +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB)
     &                      +CD2(3,NX)
     &                      +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
               IF((FLUXB.LT.0.D0).AND.
     &            ((NXABSL.GT.0).OR.(DBLE(CSKXB).GT.0.D0))) THEN
                  DO 5021 J=1,NS0
                     II=NF+I-J
                     JJ=NCF+J
                     CF(II,JJ)=0.D0
 5021             CONTINUE
               ENDIF
 5022      CONTINUE
         ENDIF
C
         NCF=NCF+NS0
         NF=NSF-NCF+NWH
         DO 5030 I=2-NF,0
            DO 5030 J=1,MIN(NF+I-1,NT1)
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=(0.D0,0.D0)
 5030    CONTINUE
         DO 5040 I=1,NSE
C       ....................................................
C           CPHASE=-CDEXP(DCMPLX(0.D0,DX)*CSKX(NS+I-10))
            IMODE = ( I-1 ) / 2 + 1
            CARG=DCMPLX(0.D0,DX)*CSKX((NX-1)*2+IMODE)
            IF ( MOD(I,2) .EQ. 0 )  THEN
               CPHASE=-CDEXP(-CARG)
            ELSE
               CPHASE=-CDEXP( CARG)
            ENDIF
C       ....................................................
            DO 5040 J=1,NT1
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=CPHASE*CF(J,NS+I)
 5040    CONTINUE
C
         IF(NX.EQ.ABS(NXABSL)) THEN
            DO 5042 I=1,NSE
               IMODE=(I-1)/2+1
               CSKXB=CSKX((NX-1)*3+IMODE)
               CSPXB=CSPX((NX-1)*3+IMODE)
               CSPZB=CSPZ((NX-1)*3+IMODE)
               IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
               IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
               FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB
     &               -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB 
     &                      +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB)
     &                      +CD2(3,NX)
     &                      +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
               IF((FLUXB.GT.0.D0).AND.
     &            ((NXABSL.GT.0).OR.(DBLE(CSKXB).LT.0.D0))) THEN
                  DO 5041 J=1,NT1
                     II=NF+I-J
                     JJ=NCF+J
                     CF(II,JJ)=0.D0
 5041             CONTINUE
               ENDIF
 5042      CONTINUE
         ENDIF
C
         NS=NS+NS1
         NSF=NSF+NSE
C
 6000 CONTINUE
C
      NS0=NT1
      NSE=2
      NF=NSF-NCF+NWH
      DO 6010 I=1,NSE
         DO 6010 J=1,NS0
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=CGOT(I,J)
 6010 CONTINUE
      DO 6020 I=NSE+1,NW+NS0-NF
         DO 6020 J=MAX(I-NW+NF,1),NS0
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=(0.D0,0.D0)
 6020 CONTINUE
C
      NCF=NCF+NS0
      NSF=NSF+NSE
         IF(NCF.NE.NSF) GO TO 9300
C
      DO 6110 J=1,4
         CA(J)=CGIN(3,J)
         CA(NSF-4+J)=-CGOT(3,J)
 6110 CONTINUE
      DO 6120 J=4,NSF-4
         CA(J)=(0.D0,0.D0)
 6120 CONTINUE
C
      CALL BANDCD(CF,CA,NSF,NW,6*MATLM+5,IND)
         IF(IND.NE.0) WRITE(6,601) IND
      IERR=0
      RETURN
C
 9300 WRITE(6,604) NSF,NCF
      IERR=1
C
      RETURN
C
  601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
  604 FORMAT(1H ,'!! ERROR IN CLBAND : NSF,NCF = ',2I5)
      END
C
C     ******* ELECTROMAGNETIC FIELD IN PLASMA *******
C
      SUBROUTINE W1EPWB(NZ)

      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NZ
      INTEGER:: NS,NX,K,KK
      REAL(rkind):: RKV,RCE,DX,PABS0,PABS1,PABSL,FLUX1,FLUX2
      COMPLEX(rkind):: CA1,CA2,CA3,CA4
      COMPLEX(rkind):: CSKX1,CSKX2,CSPX1,CSPX2,CSPZ1,CSPZ2
      COMPLEX(rkind):: CPH1,CPH2,CPH3,CPH4
      COMPLEX(rkind),DIMENSION(NXPM):: CEX,CEY,CEZ,CDEX,CDEY,CDEZ

      RKV=2.D6*PI*RF/VC
C
      RCE=VC*EPS0
C
      DO 100 NS=1,NSMAX
      DO 100 NX=1,NXP
         PABS(NX,NS)=0.D0
  100 CONTINUE
      DO 110 NX=1,NXP
         FLUX(NX)=0.D0
  110 CONTINUE
C
      DO 200 NX = 1 , NXP
         DX = RKV * ( XA( NX+1 ) - XA( NX ) )
         K  = (NX-1)*4 + 2
         KK = (NX-1)*2
C
         CA1=CA(K+1)
         CA2=CA(K+2)
         CA3=CA(K+3)
         CA4=CA(K+4)
         CSKX1=CSKX(KK+1)*CI
         CSKX2=CSKX(KK+2)*CI
         CSPX1=CSPX(KK+1)
         CSPX2=CSPX(KK+2)
         CSPZ1=CSPZ(KK+1)
         CSPZ2=CSPZ(KK+2)
C
         CPH1  = CDEXP( 0.5D0*DX*CSKX1 )
         CPH2  = CDEXP(-0.5D0*DX*CSKX1 )
         CPH3  = CDEXP( 0.5D0*DX*CSKX2 )
         CPH4  = CDEXP(-0.5D0*DX*CSKX2 )
C
         CEX(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPX1
     &            +(CA3*CPH3 +CA4*CPH4 )*CSPX2
         CEY(NX) = (CA1*CPH1 +CA2*CPH2 )
     &            +(CA3*CPH3 +CA4*CPH4 )
         CEZ(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPZ1
     &            +(CA3*CPH3 -CA4*CPH4 )*CSPZ2
         CDEX(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPX1*CSKX1
     &            +(CA3*CPH3 -CA4*CPH4 )*CSPX2*CSKX2
         CDEY(NX) = (CA1*CPH1 -CA2*CPH2 )      *CSKX1
     &            +(CA3*CPH3 -CA4*CPH4 )      *CSKX2
         CDEZ(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPZ1*CSKX1
     &            +(CA3*CPH3 +CA4*CPH4 )*CSPZ2*CSKX2
C
         CE2DA(NZ,NX,1) = CEX(NX)
         CE2DA(NZ,NX,2) = CEY(NX)
         CE2DA(NZ,NX,3) = CEZ(NX)
C
  200    CONTINUE
C
      DO 2000 NS=1,NSMAX
         DO 2000 NX = 1 , NXP
            DX = XA( NX+1 ) - XA( NX )
C
            PABS0 =(CM0(1,NX,NS)* DCONJG(CEX(NX))*CEX(NX)
     &             +CM0(2,NX,NS)* DCONJG(CEX(NX))*CEY(NX)
     &             +CM0(3,NX,NS)* DCONJG(CEY(NX))*CEY(NX)
     &             -CM0(2,NX,NS)* DCONJG(CEY(NX))*CEX(NX)
     &             +CM0(4,NX,NS)* DCONJG(CEZ(NX))*CEZ(NX))*CI
            PABS1 =(CM1(1,NX,NS)*(DCONJG(CEX(NX))*CDEZ(NX))
     &             +CM1(2,NX,NS)*(
     &                           -DCONJG(CDEY(NX))*CEZ(NX))
     &             +CM1(1,NX,NS)*(
     &                           -DCONJG(CDEZ(NX))*CEX(NX))
     &             -CM1(2,NX,NS)*(DCONJG(CEZ(NX))*CDEY(NX)))
C
            PABSL=-RCE*RKV*DX*(PABS0 + PABS1)
            PABS(NX,NS) = PABS(NX,NS) + PABSL
 2000 CONTINUE
C
      DO 3000 NX=1,NXP
            FLUX1 =(
     &             +CD1(2,NX)* DCONJG(CEY(NX))*CEZ(NX)
     &             +CD1(1,NX)* DCONJG(CEZ(NX))*CEX(NX)
     &                                                 )*(-1)
            FLUX2 =((CD2(1,NX))* DCONJG(CEX(NX))*CDEX(NX)
     &             +(CD2(2,NX))* DCONJG(CEX(NX))*CDEY(NX)
     &             +(CD2(3,NX))* DCONJG(CEY(NX))*CDEY(NX)
     &             -(CD2(2,NX))* DCONJG(CEY(NX))*CDEX(NX)
     &             +(CD2(4,NX))* DCONJG(CEZ(NX))*CDEZ(NX))*CI
C
            FLUX(NX)=FLUX(NX)
     &              +RCE*(FLUX1+FLUX2)
 3000 CONTINUE
C
      RETURN
      END
C
C     ******* LOCAL WAVE NUMBER AND POLARIZATION *******
C
      SUBROUTINE W1WKXD
C
      USE w1comm
      IMPLICIT NONE
      INTEGER:: NX,N,NS,KK
      REAL(rkind):: RKV,DX
      COMPLEX(rkind):: CDTXX,CDTXY,CDTXZ,CDTYY,CDTYZ,CDTZZ
      COMPLEX(rkind):: CBXZ,CBYZ,CBZX,CBZY,CDETIP
      COMPLEX(rkind):: CCXX,CCXY,CCYX,CCYY,CCZZ
      COMPLEX(rkind):: CK

      COMPLEX(rkind),DIMENSION(4,NXPM):: CSOL
      ReAL(rkind),PARAMETER:: EXPARG=80.D0
C
      RKV=2.D6*PI*RF/VC
C
      DO 100 NX = 1 , NXP
         CDTXX = CD2( 3,NX )
     &           /( CD2( 1,NX )*CD2( 3,NX ) + CD2( 2,NX )*CD2( 2,NX ) )
         CDTXY =-CD2( 2,NX )
     &           /( CD2( 1,NX )*CD2( 3,NX ) + CD2( 2,NX )*CD2( 2,NX ) )
         CDTYY = CD2( 1,NX )
     &           /( CD2( 1,NX )*CD2( 3,NX ) + CD2( 2,NX )*CD2( 2,NX ) )
         CDTZZ = 1.D0 / CD2( 4 , NX )
C
         CBXZ  =  CD1( 1,NX ) * CDTXX + CD1( 2,NX ) * CDTXY
         CBYZ  = -CD1( 1,NX ) * CDTXY + CD1( 2,NX ) * CDTYY
         CBZX  =  CD1( 1,NX ) * CDTZZ
         CBZY  = -CD1( 2,NX ) * CDTZZ
C
         CCXX  =  CD0( 1,NX ) * CDTXX - CD0( 2,NX ) * CDTXY
         CCXY  =  CD0( 2,NX ) * CDTXX + CD0( 3,NX ) * CDTXY
         CCYX  = -CD0( 1,NX ) * CDTXY - CD0( 2,NX ) * CDTYY
         CCYY  = -CD0( 2,NX ) * CDTXY + CD0( 3,NX ) * CDTYY
         CCZZ  =  CD0( 4,NX ) * CDTZZ
C
         CSOL(1,NX) = ( 1.D0 , 0.D0 )
         CSOL(2,NX) =( CCXX + CCYY + CCZZ - CBXZ * CBZX - CBYZ * CBZY )
         CSOL(3,NX) = ( CCXX * CCYY + CCYY * CCZZ + CCZZ * CCXX
     &               + CCXY * CBZX * CBYZ + CCYX * CBZY * CBXZ
     &               - CCYY * CBXZ * CBZX - CCXX * CBYZ * CBZY
     &               - CCXY * CCYX )
         CSOL(4,NX) = ( CCXX * CCYY - CCXY * CCYX ) * CCZZ
  100 CONTINUE
C
      CALL W1KSOL ( CSOL , NXP , 1.D-14 )
C
      DO 300 N = 1 , 3
         DO 200 NX = 1 , NXP
            NS = (NX-1)*6 + (N-1)*2 + 1
            KK = (NX-1)*3 + N
            DX = RKV*(XA(NX+1)-XA(NX))
            CK = SQRT(CSOL(N,NX))
            IF(ABS(IMAG(CK)*DX).GT.EXPARG) THEN
               IF(IMAG(CK).GE.0.D0) THEN
                  CK = DCMPLX(DBLE(CK), EXPARG/DX)
               ELSE
                  CK = DCMPLX(DBLE(CK),-EXPARG/DX)
               ENDIF
            ENDIF
            CSKX ( KK ) =   CK
C
            CDTXX = CD0(1,NX) + CSOL(N,NX)*CD2(1,NX)
            CDTXY = CD0(2,NX) + CSOL(N,NX)*CD2(2,NX)
            CDTXZ = CD1(1,NX) * CSKX( KK )
            CDTYY = CD0(3,NX) + CSOL(N,NX)*CD2(3,NX)
            CDTYZ = CD1(2,NX) * CSKX( KK )
            CDTZZ = CD0(4,NX) + CSOL(N,NX)*CD2(4,NX)
C
            CDETIP = 1.D0 / ( CDTXX*CDTYZ + CDTXZ*CDTXY )
            CSPX ( KK ) = (  CDTYY*CDTXZ - CDTXY*CDTYZ ) * CDETIP
            CSPZ ( KK ) =-(  CDTXY*CDTXY + CDTYY*CDTXX ) * CDETIP
C
            CF( 1 , NS+10 ) =  ( 1.D0 , 0.D0 )
            CF( 1 , NS+11 ) =  ( 1.D0 , 0.D0 )
            CF( 2 , NS+10 ) =  (  CD2(2,NX)*CSPX( KK )
     &                          - CD2(3,NX)          )*CSKX( KK )
     &                          - CD1(2,NX)*CSPZ( KK )
            CF( 2 , NS+11 ) = -(  CD2(2,NX)*CSPX( KK )
     &                          - CD2(3,NX)          )*CSKX( KK )
     &                          + CD1(2,NX)*CSPZ( KK )
            CF( 3 , NS+10 ) =   CSPZ( KK )
            CF( 3 , NS+11 ) = - CSPZ( KK )
            CF( 4 , NS+10 ) =   CD1(1,NX)*CSPX( KK )
     &                        + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
            CF( 4 , NS+11 ) =   CD1(1,NX)*CSPX( KK )
     &                        + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
            CF( 5 , NS+10 ) =  (  CD2(1,NX)*CSPX( KK )
     &                          + CD2(2,NX)          ) * CSKX( KK )
            CF( 5 , NS+11 ) = -(  CD2(1,NX)*CSPX( KK )
     &                          + CD2(2,NX)          ) * CSKX( KK )
            CF( 6 , NS+10 ) =   CSPX( KK )
            CF( 6 , NS+11 ) =   CSPX( KK )
  200    CONTINUE
  300 CONTINUE
      RETURN
      END
C
C     ****** SOLUTION OF COMPLEX CUBIC EQUATION ******
C
      SUBROUTINE W1KSOL(A,NMAX,EPS)
C
      USE w1comm
      IMPLICIT NONE
      COMPLEX(rkind),DIMENSION(0:3,NXPM),INTENT(INOUT):: A
      INTEGER,INTENT(IN):: NMAX
      REAL(rkind),INTENT(IN):: EPS
      INTEGER:: ICHECK,N,I
      REAL(rkind):: W8,AF1,AF2,AF3,EPSARY(NXPM),EPSMAX,SQEPS
      COMPLEX(rkind):: Z1(NXPM),Z2(NXPM),Z3(NXPM),F1,F2,F3,FF3,CWA,CWB
      INTEGER,PARAMETER:: MAXCNT=50
      INTERFACE
         FUNCTION DCBRT(X)
         USE bpsd_kinds
         IMPLICIT NONE
         REAL(rkind):: X
         REAL(rkind):: DCBRT
         END FUNCTION
      END INTERFACE
C
C     ======( INITIALIZATION )======
      ICHECK = 0
      SQEPS = SQRT( EPS )
      DO 100 N = 1 , NMAX
C     ======( NORMARIZATION )======
         W8 = 4.D0 / (    CDABS( A( 0 , N ) )
     &                  + CDABS( A( 1 , N ) )
     &                  + CDABS( A( 2 , N ) )
     &                  + CDABS( A( 3 , N ) )    )
         A( 0 , N ) = A( 0 , N ) * W8
         A( 1 , N ) = A( 1 , N ) * W8
         A( 2 , N ) = A( 2 , N ) * W8
         A( 3 , N ) = A( 3 , N ) * W8
C     ======( TRIAL SOLUTION )======
C         W8 = .2D0 * DCBRT( CDABS( A( 3 , N )/A( 0 , N ) ) )
         W8 = .2D0 * ABS( A( 3 , N )/A( 0 , N ) )**(1.D0/3.D0)
         Z1( N ) = DCMPLX( 0.D0 ,      W8 )
         Z2( N ) = DCMPLX( - W8 ,      W8 )
         Z3( N ) = DCMPLX( 0.D0 , 2.D0*W8 )
         EPSARY( N ) = 1.E0
  100    CONTINUE
  200 CONTINUE
      DO 300 N = 1 , NMAX
         AF1 = CDABS(((A(0,N)*Z1(N)+A(1,N))*Z1(N)+A(2,N))*Z1(N)+A(3,N))
         AF2 = CDABS(((A(0,N)*Z2(N)+A(1,N))*Z2(N)+A(2,N))*Z2(N)+A(3,N))
         AF3 = CDABS(((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))
         IF ( AF2 .LT. AF1 ) THEN
            CWA     = Z2( N )
            Z2( N ) = Z1( N )
            Z1( N ) = CWA
            ENDIF
         IF ( AF3 .LT. AF2 ) THEN
            IF ( AF3 .LT. AF1 ) THEN
               CWA     = Z1( N )
               CWB     = Z2( N )
               Z1( N ) = Z3( N )
               Z2( N ) = CWA
               Z3( N ) = CWB
            ELSE
               CWA     = Z2( N )
               Z2( N ) = Z3( N )
               Z3( N ) = CWA
               ENDIF
            ENDIF
  300 CONTINUE
C
      I = 1
  400 CONTINUE
        EPSMAX = 0.E0
        DO 500 N = 1 , NMAX
          IF ( EPSARY( N ) .GE. EPS ) THEN
            F1 = ( ( 3.D0*A(0,N)*Z1(N)+2.D0*A(1,N) )*Z1(N)+A(2,N) )
     &         /(((A(0,N)*Z1(N)+A(1,N))*Z1(N)+A(2,N))*Z1(N)+A(3,N))
            F2 = ( ( 3.D0*A(0,N)*Z2(N)+2.D0*A(1,N) )*Z2(N)+A(2,N) )
     &         /(((A(0,N)*Z2(N)+A(1,N))*Z2(N)+A(2,N))*Z2(N)+A(3,N))
            F3 = ( ( 3.D0*A(0,N)*Z3(N)+2.D0*A(1,N) )*Z3(N)+A(2,N) )
     &         /(((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))
            FF3 =
     &          (((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))
C
            CWA = - ( Z2(N)-Z3(N) ) * ( Z3(N)-Z1(N) ) * ( F2-F1 )
     &            / ( ( Z3(N)-Z2(N) )*(F2-F1) + (Z1(N)-Z2(N))*(F3-F2) )
            CWB = 1.D0/F3
            IF ( CDABS(FF3) .LT. SQRT(SQEPS*CDABS(A(3,N)))  ) THEN
              CWB = 2.D0*CWB
              ENDIF
C
            IF ( CDABS(CWA) .LT. CDABS(CWB) ) THEN
              CWA = Z3( N ) - CWA
            ELSE
              CWA = Z3( N ) - CWB
              ENDIF
            Z1( N ) = Z2( N )
            Z2( N ) = Z3( N )
            Z3( N ) = CWA
C
            EPSARY(N) =
     &      CDABS(((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))
     &            /(CDABS( A(0,N)*Z3(N)**3 ) + CDABS( A(1,N)*Z3(N)**2 )
     &             +CDABS( A(2,N)*Z3(N)    ) + CDABS( A(3,N) ) )
            EPSMAX      = DMAX1( EPSARY( N ) , EPSMAX )
            IF ( EPSARY( N ) .LT. EPS ) THEN
              A(2,N) = ( (CWA*A(0,N)+A(1,N))*CWA+A(2,N))/A(0,N)
              A(1,N) =   (CWA*A(0,N)+A(1,N)) / A(0,N)
              A(0,N) =   CWA
              ENDIF
            ENDIF
  500     CONTINUE
          IF ( EPSMAX .LT. EPS ) GOTO 800
          I = I + 1
          IF ( I .GE. MAXCNT ) GOTO 600
          GOTO 400
  600 IF ( ICHECK .EQ. 0 ) THEN
         ICHECK = 1
         DO 700 N = 1 , NMAX
            IF ( EPSARY( N ) .GE. EPS )  THEN
               W8 = .2D0 * DCBRT( CDABS( A( 3 , N )/A( 0 , N ) ) )
               Z1( N ) = DCMPLX( -3.D0 * W8 , 2.D0 * W8 )
               Z2( N ) = DCMPLX( -2.D0 * W8 , 2.D0 * W8 )
               Z3( N ) = DCMPLX( -3.D0 * W8 , 3.D0 * W8 )
         AF1 = CDABS(((A(0,N)*Z1(N)+A(1,N))*Z1(N)+A(2,N))*Z1(N)+A(3,N))
         AF2 = CDABS(((A(0,N)*Z2(N)+A(1,N))*Z2(N)+A(2,N))*Z2(N)+A(3,N))
         AF3 = CDABS(((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))
         IF ( AF2 .LT. AF1 ) THEN
            CWA     = Z2( N )
            Z2( N ) = Z1( N )
            Z1( N ) = CWA
            ENDIF
         IF ( AF3 .LT. AF2 ) THEN
            IF ( AF3 .LT. AF1 ) THEN
               CWA     = Z1( N )
               CWB     = Z2( N )
               Z1( N ) = Z3( N )
               Z2( N ) = CWA
               Z3( N ) = CWB
            ELSE
               CWA     = Z2( N )
               Z2( N ) = Z3( N )
               Z3( N ) = CWA
               ENDIF
            ENDIF
               ENDIF
  700       CONTINUE
         GOTO 200
      ELSE
         WRITE(6,*)
     &     '**** OVER MAX ITERATION COUNT*',MAXCNT,'IN SOBR.W1KSOL ****'
         STOP
      ENDIF
  800 CONTINUE
C     WRITE(6,*) '* ITERATION ',I
      DO 900 N = 1 , NMAX
        AF1 = DREAL( A(1,N) )
        AF2 = DIMAG( A(1,N) )
        IF ( (AF1.GT.1.D32) .OR. (AF2.GT.1.D32) ) THEN
          CWA = A(1,N)*CDSQRT( 1.D0 - 4.D0*A(2,N)/A(1,N)/A(1,N) )
        ELSE
          CWA = CDSQRT( A(1,N)*A(1,N) - 4.D0*A(2,N) )
          ENDIF
        IF ( CDABS(CWA) .LT. 1.D-70 ) THEN
          A(1,N) = -0.5D0*A(1,N)
          A(2,N) = A(1,N)
        ELSE
          F1 = A(1,N)
          F2 = A(2,N)
          IF ( AF1 .GT. 0.D0 ) THEN
            A(1,N) = ( - F1 - CWA ) * 0.5D0
            A(2,N) = - 2.D0 * F2 / ( F1 + CWA )
          ELSE IF ( AF1 .LT. 0.D0 ) THEN
            A(1,N) = 2.D0 * F2 / ( - F1 + CWA )
            A(2,N) = ( - F1 + CWA ) * 0.5D0
          ELSE IF ( AF2 .GE. 0.D0 ) THEN
            A(1,N) = ( - F1 - CWA ) * 0.5D0
            A(2,N) = - 2.D0 * F2 / ( F1 + CWA )
          ELSE
            A(1,N) = + 2.D0 * F2 / ( - F1 + CWA )
            A(2,N) = ( - F1 + CWA ) * 0.5D0
            ENDIF
          ENDIF
 900  CONTINUE
C
      RETURN
      END
C
C     ******* BAND MATRIX COEFFICIENT *******
C
      SUBROUTINE W1BNDD(IERR,NXABSL)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NXABSL
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: NS,NW,NWH,NCF,NSF,NSE,NT1,NF,I,J,II,JJ
      INTEGER:: NS0,NS1,NS2,NX,IMODE,IND
      REAL(rkind):: RKV,X1,X2,DX,FLUXB
      COMPLEX(rkind):: CSKXB,CSPXB,CSPZB,CARG,CPHASE
C
      RKV=2.D6*PI*RF/VC
C
      NS  = 10
      NW  = 3*6 - 1
      NWH = 6 + 6/2
C
      NCF=0
      NSF=0
      NSE=2
      NT1=5
C
      NF=NSF-NCF+NWH
      DO 10 I=2-NF,0
         DO 10 J=1,MIN(NF+I-1,NT1)
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=(0.D0,0.D0)
   10 CONTINUE
      DO 20 I=1,NSE
         DO 20 J=1,NT1
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=-CGIN(I,J)
   20 CONTINUE
      NSF=NSF+NSE
C
      X2  = XA ( 1 )
      NS2 = 6
C
      DO 6000 NX=1,NXP
         X1=X2
         X2=XA(NX+1)
         DX  = RKV * ( X2-X1 )
         NS0 = NT1
         NS1 = NS2
         NS2 = 6
         NSE = NS1
         NT1 = NS1
         IF( NX .EQ. NXP ) NT1 = 5
C
         NF=NSF-NCF+NWH
         DO 5010 I=1,NSE
            DO 5010 J=1,NS0
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=CF(J,NS+I)
 5010    CONTINUE
         DO 5020 I=NSE+1,NW+NS0-NF
            DO 5020 J=MAX(I-NW+NF,1),NS0
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=(0.D0,0.D0)
 5020    CONTINUE
C
         IF(NX.EQ.ABS(NXABSL)) THEN
            DO 5022 I=1,NSE
               IMODE=(I-1)/2+1
               CSKXB=CSKX((NX-1)*3+IMODE)
               CSPXB=CSPX((NX-1)*3+IMODE)
               CSPZB=CSPZ((NX-1)*3+IMODE)
               IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
               IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
               FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB
     &               -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB 
     &                      +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB)
     &                      +CD2(3,NX)
     &                      +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
               IF((FLUXB.LT.0.D0).AND.
     &            ((NXABSL.GT.0).OR.(DBLE(CSKXB).GT.0.D0))) THEN
                  DO 5021 J=1,NS0
                     II=NF+I-J
                     JJ=NCF+J
                     CF(II,JJ)=0.D0
 5021             CONTINUE
               ENDIF
 5022      CONTINUE
         ENDIF
C
         NCF=NCF+NS0
         NF=NSF-NCF+NWH
         DO 5030 I=2-NF,0
            DO 5030 J=1,MIN(NF+I-1,NT1)
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=(0.D0,0.D0)
 5030    CONTINUE
         DO 5040 I=1,NSE
C       ....................................................
C           CPHASE=-CDEXP(DCMPLX(0.D0,DX)*CSKX(NS+I-10))
            IMODE = ( I-1 ) / 2 + 1
            CARG=DCMPLX(0.D0,DX)*CSKX((NX-1)*3+IMODE)
            IF ( MOD(I,2) .EQ. 0 )  THEN
               CPHASE=-CDEXP(-CARG)
            ELSE
               CPHASE=-CDEXP( CARG)
            ENDIF
C       ....................................................
            DO 5040 J=1,NT1
               II=NF+I-J
               JJ=NCF+J
               CF(II,JJ)=CPHASE*CF(J,NS+I)
 5040    CONTINUE
C
         IF(NX.EQ.ABS(NXABSL)) THEN
            DO 5042 I=1,NSE
               IMODE=(I-1)/2+1
               CSKXB=CSKX((NX-1)*3+IMODE)
               CSPXB=CSPX((NX-1)*3+IMODE)
               CSPZB=CSPZ((NX-1)*3+IMODE)
               IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
               IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
               FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB
     &               -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB 
     &                      +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB)
     &                      +CD2(3,NX)
     &                      +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
               IF((FLUXB.GT.0.D0).AND.
     &            ((NXABSL.GT.0).OR.(DBLE(CSKXB).LT.0.D0))) THEN
                  DO 5041 J=1,NT1
                     II=NF+I-J
                     JJ=NCF+J
                     CF(II,JJ)=0.D0
 5041             CONTINUE
               ENDIF
 5042      CONTINUE
         ENDIF
C
         NS=NS+NS1
         NSF=NSF+NSE
C
 6000 CONTINUE
C
      NS0=NT1
      NSE=2
      NF=NSF-NCF+NWH
      DO 6010 I=1,NSE
         DO 6010 J=1,NS0
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=CGOT(I,J)
 6010 CONTINUE
      DO 6020 I=NSE+1,NW+NS0-NF
         DO 6020 J=MAX(I-NW+NF,1),NS0
            II=NF+I-J
            JJ=NCF+J
            CF(II,JJ)=(0.D0,0.D0)
 6020 CONTINUE
C
      NCF=NCF+NS0
      NSF=NSF+NSE
         IF(NCF.NE.NSF) GO TO 9300
C
      DO 6110 J=1,5
         CA(J)=CGIN(3,J)
         CA(NSF-5+J)=-CGOT(3,J)
 6110 CONTINUE
      DO 6120 J=6,NSF-5
         CA(J)=(0.D0,0.D0)
 6120 CONTINUE
C
      CALL BANDCD(CF,CA,NSF,17,6*MATLM+5,IND)
         IF(IND.NE.0) WRITE(6,601) IND
      IERR=0
      RETURN
C
 9300 WRITE(6,604) NSF,NCF
      IERR=1
C
      RETURN
C
  601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
  604 FORMAT(1H ,'!! ERROR IN CLBAND : NSF,NCF = ',2I5)
      END
C
C     ******* ELECTROMAGNETIC FIELD IN PLASMA *******
C
      SUBROUTINE W1EPWD(NZ)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NZ
      INTEGER:: NS,NX,K,KK
      REAL(rkind):: RKV,RCE,DX,PABS0,PABS1,PABS2,PABSL,FLUX1,FLUX2
      COMPLEX(rkind):: CA1,CA2,CA3,CA4,CA5,CA6
      COMPLEX(rkind):: CSKX1,CSKX2,CSKX3,CSPX1,CSPX2,CSPX3
      COMPLEX(rkind):: CSPZ1,CSPZ2,CSPZ3
      COMPLEX(rkind):: CPH1,CPH2,CPH3,CPH4,CPH5,CPH6
      COMPLEX(rkind),DIMENSION(NXPM):: CEX,CEY,CEZ,CDEX,CDEY,CDEZ

      RKV=2.D6*PI*RF/VC
C
      RCE=VC*EPS0
C
      DO 100 NS=1,NSMAX
      DO 100 NX=1,NXP
         PABS(NX,NS)=0.D0
  100 CONTINUE
      DO 110 NX=1,NXP
         FLUX(NX)=0.D0
  110 CONTINUE
C
      DO 200 NX = 1 , NXP
         DX = RKV * ( XA( NX+1 ) - XA( NX ) )
         K  = (NX-1)*6 + 2
         KK = (NX-1)*3
C
         CA1=CA(K+1)
         CA2=CA(K+2)
         CA3=CA(K+3)
         CA4=CA(K+4)
         CA5=CA(K+5)
         CA6=CA(K+6)
         CSKX1=CSKX(KK+1)*CI
         CSKX2=CSKX(KK+2)*CI
         CSKX3=CSKX(KK+3)*CI
         CSPX1=CSPX(KK+1)
         CSPX2=CSPX(KK+2)
         CSPX3=CSPX(KK+3)
         CSPZ1=CSPZ(KK+1)
         CSPZ2=CSPZ(KK+2)
         CSPZ3=CSPZ(KK+3)
C
         CPH1  = CDEXP( 0.5D0*DX*CSKX1 )
         CPH2  = CDEXP(-0.5D0*DX*CSKX1 )
         CPH3  = CDEXP( 0.5D0*DX*CSKX2 )
         CPH4  = CDEXP(-0.5D0*DX*CSKX2 )
         CPH5  = CDEXP( 0.5D0*DX*CSKX3 )
         CPH6  = CDEXP(-0.5D0*DX*CSKX3 )
C
         CEX(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPX1
     &            +(CA3*CPH3 +CA4*CPH4 )*CSPX2
     &            +(CA5*CPH5 +CA6*CPH6 )*CSPX3
         CEY(NX) = (CA1*CPH1 +CA2*CPH2 )
     &            +(CA3*CPH3 +CA4*CPH4 )
     &            +(CA5*CPH5 +CA6*CPH6 )
         CEZ(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPZ1
     &            +(CA3*CPH3 -CA4*CPH4 )*CSPZ2
     &            +(CA5*CPH5 -CA6*CPH6 )*CSPZ3
         CDEX(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPX1*CSKX1
     &            +(CA3*CPH3 -CA4*CPH4 )*CSPX2*CSKX2
     &            +(CA5*CPH5 -CA6*CPH6 )*CSPX3*CSKX3
         CDEY(NX) = (CA1*CPH1 -CA2*CPH2 )      *CSKX1
     &            +(CA3*CPH3 -CA4*CPH4 )      *CSKX2
     &            +(CA5*CPH5 -CA6*CPH6 )      *CSKX3
         CDEZ(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPZ1*CSKX1
     &            +(CA3*CPH3 +CA4*CPH4 )*CSPZ2*CSKX2
     &            +(CA5*CPH5 +CA6*CPH6 )*CSPZ3*CSKX3
C
         CE2DA(NZ,NX,1) = CEX(NX)
         CE2DA(NZ,NX,2) = CEY(NX)
         CE2DA(NZ,NX,3) = CEZ(NX)
C
  200    CONTINUE
C
      DO 2000 NS=1,NSMAX
         DO 2000 NX = 1 , NXP
            DX = XA( NX+1 ) - XA( NX )
C
            PABS0 =(CM0(1,NX,NS)* DCONJG(CEX(NX))*CEX(NX)
     &             +CM0(2,NX,NS)* DCONJG(CEX(NX))*CEY(NX)
     &             +CM0(3,NX,NS)* DCONJG(CEY(NX))*CEY(NX)
     &             -CM0(2,NX,NS)* DCONJG(CEY(NX))*CEX(NX)
     &             +CM0(4,NX,NS)* DCONJG(CEZ(NX))*CEZ(NX))*CI
            PABS1 =(CM1(1,NX,NS)*(DCONJG(CEX(NX))*CDEZ(NX))
     &             +CM1(2,NX,NS)*(
     &                           -DCONJG(CDEY(NX))*CEZ(NX))
     &             +CM1(1,NX,NS)*(
     &                           -DCONJG(CDEZ(NX))*CEX(NX))
     &             -CM1(2,NX,NS)*(DCONJG(CEZ(NX))*CDEY(NX)))
            PABS2 =(CM2(1,NX,NS)* DCONJG(CDEX(NX))*CDEX(NX)
     &             +CM2(2,NX,NS)* DCONJG(CDEX(NX))*CDEY(NX)
     &             +CM2(3,NX,NS)* DCONJG(CDEY(NX))*CDEY(NX)
     &             -CM2(2,NX,NS)* DCONJG(CDEY(NX))*CDEX(NX)
     &             +CM2(4,NX,NS)* DCONJG(CDEZ(NX))*CDEZ(NX))*CI
C
            PABSL=-RCE*RKV*DX*(PABS0 + PABS1 + PABS2)
            PABS(NX,NS) = PABS(NX,NS) + PABSL
 2000 CONTINUE
C
      DO 3000 NX=1,NXP
            FLUX1 =(
     &             +CD1(2,NX)* DCONJG(CEY(NX))*CEZ(NX)
     &             +CD1(1,NX)* DCONJG(CEZ(NX))*CEX(NX)
     &                                                 )*(-1)
            FLUX2 =((CD2(1,NX))* DCONJG(CEX(NX))*CDEX(NX)
     &             +(CD2(2,NX))* DCONJG(CEX(NX))*CDEY(NX)
     &             +(CD2(3,NX))* DCONJG(CEY(NX))*CDEY(NX)
     &             -(CD2(2,NX))* DCONJG(CEY(NX))*CDEX(NX)
     &             +(CD2(4,NX))* DCONJG(CEZ(NX))*CDEZ(NX))*CI
C
            FLUX(NX)=FLUX(NX)
     &              +RCE*(FLUX1+FLUX2)
 3000 CONTINUE
C
      RETURN
      END
