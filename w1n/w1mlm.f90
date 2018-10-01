MODULE w1mlm

!     ******* LOCAL WAVE NUMBER AND POLARIZATION *******

CONTAINS

  SUBROUTINE W1WKXB
    USE w1comm
    IMPLICIT NONE
    COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CSOL
    REAL(rkind),PARAMETER:: EXPARG=80.D0
    INTEGER:: NX,N,NK,KK
    REAL(rkind):: RKV,DX
    COMPLEX(rkind):: CDSP0,CDSP1,CDSP2,CDET,CK
    COMPLEX(rkind):: CDTXX,CDTXY,CDTXZ,CDTYY,CDTYZ,CDTZZ,CDETIP

    MWID = 3*4 - 1
    MLEN = 6*NXPMAX+10

    IF(ALLOCATED(CF)) DEALLOCATE(CF)
    IF(ALLOCATED(CSOL)) DEALLOCATE(CSOL)

    ALLOCATE(CF(MWID,MLEN))
    ALLOCATE(CSOL(2,NXPMAX))

    RKV=2.D6*PI*RF/VC

    DO NX = 1 , NXPMAX
       CDSP2 = CD0(1,NX)+CD1(1,NX)**2
       CDSP1 =-CD0(1,NX)*(CD0(3,NX)+CD0(4,NX))-CD0(2,NX)**2 &
              +2.D0*CD0(2,NX)*CD1(1,NX)*CD1(2,NX) &
              +CD0(1,NX)*CD1(2,NX)**2 &
              -CD0(3,NX)*CD1(1,NX)**2
       CDSP0 = CD0(4,NX)*(CD0(1,NX)*CD0(3,NX)+CD0(2,NX)**2)

       CDET  = SQRT(CDSP1**2-4.D0*CDSP2*CDSP0)
       CSOL(1,NX) = (-CDSP1+CDET)/(2.D0*CDSP2)
       CSOL(2,NX) = (-CDSP1-CDET)/(2.D0*CDSP2)
    END DO

    DO N = 1 , 2
       DO NX = 1 , NXPMAX
          NK = (NX-1)*4 + (N-1)*2 + 1
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

          CDTXX = CD0(1,NX) + CSOL(N,NX)*CD2(1,NX)
          CDTXY = CD0(2,NX) + CSOL(N,NX)*CD2(2,NX)
          CDTXZ = CD1(1,NX) * CSKX( KK )
          CDTYY = CD0(3,NX) + CSOL(N,NX)*CD2(3,NX)
          CDTYZ = CD1(2,NX) * CSKX( KK )
          CDTZZ = CD0(4,NX) + CSOL(N,NX)*CD2(4,NX)

          CDETIP = 1.D0 / ( CDTXX*CDTYZ + CDTXZ*CDTXY )
          CSPX ( KK ) = (  CDTYY*CDTXZ - CDTXY*CDTYZ ) * CDETIP
          CSPZ ( KK ) =-(  CDTXY*CDTXY + CDTYY*CDTXX ) * CDETIP

          CF( 1 , NK+10 ) =  ( 1.D0 , 0.D0 )
          CF( 1 , NK+11 ) =  ( 1.D0 , 0.D0 )
          CF( 2 , NK+10 ) =  (  CD2(2,NX)*CSPX( KK ) &
                              - CD2(3,NX)          )*CSKX( KK ) &
                              - CD1(2,NX)*CSPZ( KK )
          CF( 2 , NK+11 ) = -(  CD2(2,NX)*CSPX( KK ) &
                              - CD2(3,NX)          )*CSKX( KK ) &
                              + CD1(2,NX)*CSPZ( KK )
          CF( 3 , NK+10 ) =   CSPZ( KK )
          CF( 3 , NK+11 ) = - CSPZ( KK )
          CF( 4 , NK+10 ) =   CD1(1,NX)*CSPX( KK ) &
                            + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
          CF( 4 , NK+11 ) =   CD1(1,NX)*CSPX( KK ) &
                            + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
       END DO
    END DO
    DEALLOCATE(CSOL)
    RETURN
  END SUBROUTINE W1WKXB

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BNDB(IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NX,NS,NWH,NCF,NSF,NSE,NT1,NF,I,J,II,JJ
    INTEGER:: NS0,NS1,NS2,IMODE
    REAL(rkind):: RKV,X1,X2,DX,FLUXB
    COMPLEX(rkind):: CSKXB,CSPXB,CSPZB,CARG,CPHASE

    RKV=2.D6*PI*RF/VC

    NS=10
    NWH = 4 + 4/2

    NCF=0
    NSF=0
    NSE=2
    NT1=4

    NF=NSF-NCF+NWH
    DO I=2-NF,0
       DO J=1,MIN(NF+I-1,NT1)
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=(0.D0,0.D0)
       END DO
    END DO
    DO I=1,NSE
       DO J=1,NT1
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=-CGIN(I,J)
       END DO
    END DO
    NSF=NSF+NSE

    X2  = XA ( 1 )
    NS2 = 4

    DO NX=1,NXPMAX
       X1=X2
       X2=XA(NX+1)
       DX  = RKV * ( X2-X1 )
       NS0 = NT1
       NS1 = NS2
       NS2 = 4
       NSE = NS1
       NT1 = NS1
       IF( NX .EQ. NXPMAX ) NT1 = 4

       NF=NSF-NCF+NWH
       DO I=1,NSE
          DO J=1,NS0
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=CF(J,NS+I)
          END DO
       END DO
       DO I=NSE+1,MWID+NS0-NF
          DO J=MAX(I-MWID+NF,1),NS0
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=(0.D0,0.D0)
          END DO
       END DO

       IF(NX.EQ.ABS(NXABS)) THEN
          DO I=1,NSE
             IMODE=(I-1)/2+1
             CSKXB=CSKX((NX-1)*2+IMODE)
             CSPXB=CSPX((NX-1)*2+IMODE)
             CSPZB=CSPZ((NX-1)*2+IMODE)
             IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
             IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
             FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB &
                   -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB &
                          +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB) &
                          +CD2(3,NX) &
                          +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
             IF((FLUXB.LT.0.D0).AND. &
               ((NXABS.GT.0).OR.(DBLE(CSKXB).GT.0.D0))) THEN
                DO J=1,NS0
                   II=NF+I-J
                   JJ=NCF+J
                   CF(II,JJ)=0.D0
                END DO
             END IF
          END DO
       END IF

       NCF=NCF+NS0
       NF=NSF-NCF+NWH
       DO I=2-NF,0
          DO J=1,MIN(NF+I-1,NT1)
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=(0.D0,0.D0)
          END DO
       END DO

       DO I=1,NSE
!       ....................................................
!           CPHASE=-CDEXP(DCMPLX(0.D0,DX)*CSKX(NS+I-10))
          IMODE = ( I-1 ) / 2 + 1
          CARG=DCMPLX(0.D0,DX)*CSKX((NX-1)*2+IMODE)
          IF ( MOD(I,2) .EQ. 0 )  THEN
             CPHASE=-CDEXP(-CARG)
          ELSE
             CPHASE=-CDEXP( CARG)
          ENDIF
!       ....................................................
          DO J=1,NT1
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=CPHASE*CF(J,MWID+I)
          END DO
       END DO

       IF(NX.EQ.ABS(NXABS)) THEN
          DO I=1,NSE
             IMODE=(I-1)/2+1
             CSKXB=CSKX((NX-1)*3+IMODE)
             CSPXB=CSPX((NX-1)*3+IMODE)
             CSPZB=CSPZ((NX-1)*3+IMODE)
             IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
             IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
             FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB &
                   -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB &
                          +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB) &
                          +CD2(3,NX) &
                          +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
             IF((FLUXB.GT.0.D0).AND. &
                ((NXABS.GT.0).OR.(DBLE(CSKXB).LT.0.D0))) THEN
                DO J=1,NT1
                   II=NF+I-J
                   JJ=NCF+J
                   CF(II,JJ)=0.D0
                END DO
             END IF
          END DO
       END IF

       NS=NS+NS1
       NSF=NSF+NSE
    END DO

    NS0=NT1
    NSE=2
    NF=NSF-NCF+NWH
    DO I=1,NSE
       DO J=1,NS0
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=CGOT(I,J)
       END DO
    END DO
    DO I=NSE+1,MWID+NS0-NF
       DO J=MAX(I-MWID+NF,1),NS0
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=(0.D0,0.D0)
       END DO
    END DO

    NCF=NCF+NS0
    NSF=NSF+NSE
    IF(NCF.NE.NSF) GO TO 9300

    DO J=1,4
       CA(J)=CGIN(3,J)
       CA(NSF-4+J)=-CGOT(3,J)
    END DO
    DO J=4,NSF-4
       CA(J)=(0.D0,0.D0)
    END DO

    CALL BANDCD(CF,CA,NSF,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR
    DEALLOCATE(CF)
    RETURN

 9300 CONTINUE
    WRITE(6,604) NSF,NCF
    IERR=1
    RETURN

601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
604 FORMAT(1H ,'!! ERROR IN CLBAND : NSF,NCF = ',2I5)
  END SUBROUTINE W1BNDB

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPWB(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CEX,CEY,CEZ,CDX,CDY,CDZ
    INTEGER:: NS,NX,K,KK
    REAL(rkind)::RKV,RCE,DX,PABS0,PABS1,PABSL,FLUX1,FLUX2
    COMPLEX(rkind):: CA1,CA2,CA3,CA4,CSKX1,CSKX2,CSPX1,CSPX2,CSPZ1,CSPZ2
    COMPLEX(rkind):: CPH1,CPH2,CPH3,CPH4

    ALLOCATE(CEX(NXPMAX),CEY(NXPMAX),CEZ(NXPMAX))
    ALLOCATE(CDX(NXPMAX),CDY(NXPMAX),CDZ(NXPMAX))

    RKV=2.D6*PI*RF/VC
    RCE=VC*EPS0

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABS(NX,NS)=0.D0
       END DO
    END DO
    DO NX=1,NXPMAX
       FLUX(NX)=0.D0
    END DO

    DO NX = 1 , NXPMAX
       DX = RKV * ( XA( NX+1 ) - XA( NX ) )
       K  = (NX-1)*4 + 2
       KK = (NX-1)*2

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

       CPH1  = CDEXP( 0.5D0*DX*CSKX1 )
       CPH2  = CDEXP(-0.5D0*DX*CSKX1 )
       CPH3  = CDEXP( 0.5D0*DX*CSKX2 )
       CPH4  = CDEXP(-0.5D0*DX*CSKX2 )

       CEX(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPX1 &
                +(CA3*CPH3 +CA4*CPH4 )*CSPX2
       CEY(NX) = (CA1*CPH1 +CA2*CPH2 ) &
                +(CA3*CPH3 +CA4*CPH4 )
       CEZ(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPZ1 &
                +(CA3*CPH3 -CA4*CPH4 )*CSPZ2
       CDX(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPX1*CSKX1 &
                +(CA3*CPH3 -CA4*CPH4 )*CSPX2*CSKX2
       CDY(NX) = (CA1*CPH1 -CA2*CPH2 )      *CSKX1 &
                +(CA3*CPH3 -CA4*CPH4 )      *CSKX2
       CDZ(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPZ1*CSKX1 &
                +(CA3*CPH3 +CA4*CPH4 )*CSPZ2*CSKX2

       CE2DA(NZ,NX,1) = CEX(NX)
       CE2DA(NZ,NX,2) = CEY(NX)
       CE2DA(NZ,NX,3) = CEZ(NX)

    END DO

    DO NS=1,NSMAX
       DO NX = 1 , NXPMAX
          DX = XA( NX+1 ) - XA( NX )

          WRITE(6,'(A,I4,1P4E12.4)') 'CM0=',NX,CM0(1,NX,NS),CM0(2,NX,NS)
          WRITE(6,'(A,I4,1P4E12.4)') 'CM0=', 0,CM0(3,NX,NS),CM0(4,NX,NS)
          WRITE(6,'(A,1P6E12.4)') 'CEX=',CEX(NX),CEY(NX),CEZ(NX)
          WRITE(6,'(A,I4,1P4E12.4)') 'CM1=',NX,CM1(1,NX,NS),CM1(2,NX,NS)
          WRITE(6,'(A,1P4E12.4)') 'CDX=',CDY(NX),CDZ(NX)
          PABS0 =(CM0(1,NX,NS)*  DCONJG(CEX(NX))*CEX(NX) &
                 +CM0(2,NX,NS)*  DCONJG(CEX(NX))*CEY(NX) &
                 +CM0(3,NX,NS)*  DCONJG(CEY(NX))*CEY(NX) &
                 -CM0(2,NX,NS)*  DCONJG(CEY(NX))*CEX(NX) &
                 +CM0(4,NX,NS)*  DCONJG(CEZ(NX))*CEZ(NX))*CI
          PABS1 =(CM1(1,NX,NS)*( DCONJG(CEX(NX))*CDZ(NX)) &
                 +CM1(2,NX,NS)*(-DCONJG(CDY(NX))*CEZ(NX)) &
                 +CM1(1,NX,NS)*(      -DCONJG(CDZ(NX))*CEX(NX)) &
                 -CM1(2,NX,NS)*(DCONJG(CEZ(NX))*CDY(NX)))

          PABSL=-RCE*RKV*DX*(PABS0 + PABS1)
          PABS(NX,NS) = PABS(NX,NS) + PABSL
       END DO
    END DO

    DO NX=1,NXPMAX
       FLUX1 =(    +CD1(2,NX)* DCONJG(CEY(NX))*CEZ(NX) &
                   +CD1(1,NX)* DCONJG(CEZ(NX))*CEX(NX))*(-1)
       FLUX2 =((CD2(1,NX))* DCONJG(CEX(NX))*CDX(NX) &
              +(CD2(2,NX))* DCONJG(CEX(NX))*CDY(NX) &
              +(CD2(3,NX))* DCONJG(CEY(NX))*CDY(NX) &
              -(CD2(2,NX))* DCONJG(CEY(NX))*CDX(NX) &
              +(CD2(4,NX))* DCONJG(CEZ(NX))*CDZ(NX))*CI
       FLUX(NX)=FLUX(NX)+RCE*(FLUX1+FLUX2)
    END DO
    RETURN
  END SUBROUTINE W1EPWB

!     ******* LOCAL WAVE NUMBER AND POLARIZATION *******

  SUBROUTINE W1WKXD
    USE w1comm
    IMPLICIT NONE
    COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CSOL
    REAL(rkind),PARAMETER:: EXPARG=80.D0
    INTEGER:: NX,N,NK,KK
    REAL(rkind):: RKV,DX
    COMPLEX(rkind):: CBXZ,CBYZ,CBZX,CBZY,CCXX,CCXY,CCYX,CCYY,CCZZ
    COMPLEX(rkind):: CK,CDTXX,CDTXY,CDTXZ,CDTYY,CDTYZ,CDTZZ,CDETIP

    RKV=2.D6*PI*RF/VC

    MWID=3*6 - 1
    MLEN=6*NXPMAX+10

    IF(ALLOCATED(CF)) DEALLOCATE(CF)
    IF(ALLOCATED(CSOL)) DEALLOCATE(CSOL)

    ALLOCATE(CF(MWID,MLEN))
    ALLOCATE(CSOL(4,NXPMAX))

    DO NX = 1 , NXPMAX
       CDTXX = CD2( 3,NX ) &
               /( CD2( 1,NX )*CD2( 3,NX ) + CD2( 2,NX )*CD2( 2,NX ) )
       CDTXY =-CD2( 2,NX ) &
               /( CD2( 1,NX )*CD2( 3,NX ) + CD2( 2,NX )*CD2( 2,NX ) )
       CDTYY = CD2( 1,NX ) &
               /( CD2( 1,NX )*CD2( 3,NX ) + CD2( 2,NX )*CD2( 2,NX ) )
       CDTZZ = 1.D0 / CD2( 4 , NX )

       CBXZ  =  CD1( 1,NX ) * CDTXX + CD1( 2,NX ) * CDTXY
       CBYZ  = -CD1( 1,NX ) * CDTXY + CD1( 2,NX ) * CDTYY
       CBZX  =  CD1( 1,NX ) * CDTZZ
       CBZY  = -CD1( 2,NX ) * CDTZZ

       CCXX  =  CD0( 1,NX ) * CDTXX - CD0( 2,NX ) * CDTXY
       CCXY  =  CD0( 2,NX ) * CDTXX + CD0( 3,NX ) * CDTXY
       CCYX  = -CD0( 1,NX ) * CDTXY - CD0( 2,NX ) * CDTYY
       CCYY  = -CD0( 2,NX ) * CDTXY + CD0( 3,NX ) * CDTYY
       CCZZ  =  CD0( 4,NX ) * CDTZZ

       CSOL(1,NX) = ( 1.D0 , 0.D0 )
       CSOL(2,NX) =( CCXX + CCYY + CCZZ - CBXZ * CBZX - CBYZ * CBZY )
       CSOL(3,NX) = ( CCXX * CCYY + CCYY * CCZZ + CCZZ * CCXX &
                    + CCXY * CBZX * CBYZ + CCYX * CBZY * CBXZ &
                    - CCYY * CBXZ * CBZX - CCXX * CBYZ * CBZY &
                    - CCXY * CCYX )
       CSOL(4,NX) = ( CCXX * CCYY - CCXY * CCYX ) * CCZZ
    END DO

    CALL W1KSOL ( CSOL , NXPMAX , 1.D-14 )

    DO N = 1 , 3
       DO NX = 1 , NXPMAX
          NK = (NX-1)*6 + (N-1)*2 + 1
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

          CDTXX = CD0(1,NX) + CSOL(N,NX)*CD2(1,NX)
          CDTXY = CD0(2,NX) + CSOL(N,NX)*CD2(2,NX)
          CDTXZ = CD1(1,NX) * CSKX( KK )
          CDTYY = CD0(3,NX) + CSOL(N,NX)*CD2(3,NX)
          CDTYZ = CD1(2,NX) * CSKX( KK )
          CDTZZ = CD0(4,NX) + CSOL(N,NX)*CD2(4,NX)

          CDETIP = 1.D0 / ( CDTXX*CDTYZ + CDTXZ*CDTXY )
          CSPX ( KK ) = (  CDTYY*CDTXZ - CDTXY*CDTYZ ) * CDETIP
          CSPZ ( KK ) =-(  CDTXY*CDTXY + CDTYY*CDTXX ) * CDETIP

          CF( 1 , NK+10 ) =  ( 1.D0 , 0.D0 )
          CF( 1 , NK+11 ) =  ( 1.D0 , 0.D0 )
          CF( 2 , NK+10 ) =  (  CD2(2,NX)*CSPX( KK ) &
                              - CD2(3,NX)          )*CSKX( KK ) &
                              - CD1(2,NX)*CSPZ( KK ) 
          CF( 2 , NK+11 ) = -(  CD2(2,NX)*CSPX( KK ) &
                              - CD2(3,NX)          )*CSKX( KK ) &
                              + CD1(2,NX)*CSPZ( KK )
          CF( 3 , NK+10 ) =   CSPZ( KK )
          CF( 3 , NK+11 ) = - CSPZ( KK )
          CF( 4 , NK+10 ) =   CD1(1,NX)*CSPX( KK ) &
                            + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
          CF( 4 , NK+11 ) =   CD1(1,NX)*CSPX( KK ) &
                            + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
          CF( 5 , NK+10 ) =  (  CD2(1,NX)*CSPX( KK ) &
                              + CD2(2,NX)          ) * CSKX( KK )
          CF( 5 , NK+11 ) = -(  CD2(1,NX)*CSPX( KK ) &
                              + CD2(2,NX)          ) * CSKX( KK )
          CF( 6 , NK+10 ) =   CSPX( KK )
          CF( 6 , NK+11 ) =   CSPX( KK )
       END DO
    END DO
    RETURN
  END SUBROUTINE W1WKXD

!     ****** SOLUTION OF COMPLEX CUBIC EQUATION ******

  SUBROUTINE W1KSOL(A,NMAX,EPS)
    USE w1comm,ONLY: rkind
    IMPLICIT NONE
    INTEGER,PARAMETER:: MAXCNT = 50
    COMPLEX(rkind),INTENT(INOUT):: A (0:3,NMAX)
    INTEGER,INTENT(IN):: NMAX
    REAL(rkind),INTENT(IN):: EPS
    INTEGER:: ICHECK,N,I
    REAL(rkind):: W8,AF1,AF2,AF3,SQEPS
    REAL(rkind):: EPSARY(NMAX),EPSMAX
    COMPLEX(rkind):: Z1(NMAX),Z2(NMAX),Z3(NMAX),F1,F2,F3,FF3,CWA,CWB

!    ======( INITIALIZATION )======
    ICHECK = 0
    SQEPS = SQRT( EPS )

    DO N = 1 , NMAX
!     ======( NORMARIZATION )======
       W8 = 4.D0 / (    CDABS( A( 0 , N ) ) &
                      + CDABS( A( 1 , N ) ) &
                      + CDABS( A( 2 , N ) ) &
                      + CDABS( A( 3 , N ) )    )
       A( 0 , N ) = A( 0 , N ) * W8
       A( 1 , N ) = A( 1 , N ) * W8
       A( 2 , N ) = A( 2 , N ) * W8
       A( 3 , N ) = A( 3 , N ) * W8
!     ======( TRIAL SOLUTION )======
!         W8 = .2D0 * DCBRT( ABS( A( 3 , N )/A( 0 , N ) ) )
       W8 = .2D0 * ABS( A( 3 , N )/A( 0 , N ) )**(1.D0/3.D0)
       Z1( N ) = DCMPLX( 0.D0 ,      W8 )
       Z2( N ) = DCMPLX( - W8 ,      W8 )
       Z3( N ) = DCMPLX( 0.D0 , 2.D0*W8 )
       EPSARY( N ) = 1.E0
    END DO

200 CONTINUE

    DO N = 1 , NMAX
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
    END DO

    I = 1

400 CONTINUE
    EPSMAX = 0.E0
    DO N = 1 , NMAX
       IF ( EPSARY( N ) .GE. EPS ) THEN
          F1 = ( ( 3.D0*A(0,N)*Z1(N)+2.D0*A(1,N) )*Z1(N)+A(2,N) ) &
             /(((A(0,N)*Z1(N)+A(1,N))*Z1(N)+A(2,N))*Z1(N)+A(3,N))
          F2 = ( ( 3.D0*A(0,N)*Z2(N)+2.D0*A(1,N) )*Z2(N)+A(2,N) ) &
             /(((A(0,N)*Z2(N)+A(1,N))*Z2(N)+A(2,N))*Z2(N)+A(3,N))
          F3 = ( ( 3.D0*A(0,N)*Z3(N)+2.D0*A(1,N) )*Z3(N)+A(2,N) ) &
             /(((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))
          FF3 = (((A(0,N)*Z3(N)+A(1,N))*Z3(N)+A(2,N))*Z3(N)+A(3,N))

          CWA = - ( Z2(N)-Z3(N) ) * ( Z3(N)-Z1(N) ) * ( F2-F1 ) &
                / ( ( Z3(N)-Z2(N) )*(F2-F1) + (Z1(N)-Z2(N))*(F3-F2) )
          CWB = 1.D0/F3
          IF ( CDABS(FF3) .LT. SQRT(SQEPS*CDABS(A(3,N)))  ) THEN
             CWB = 2.D0*CWB
          ENDIF

          IF ( CDABS(CWA) .LT. CDABS(CWB) ) THEN
             CWA = Z3( N ) - CWA
          ELSE
             CWA = Z3( N ) - CWB
          ENDIF
          Z1( N ) = Z2( N )
          Z2( N ) = Z3( N )
          Z3( N ) = CWA

          EPSARY(N) = CDABS(((A(0,N)*Z3(N)+A(1,N))*Z3(N) &
                             +A(2,N))*Z3(N)+A(3,N)) &
                    /(CDABS( A(0,N)*Z3(N)**3 ) + CDABS( A(1,N)*Z3(N)**2 ) &
                     +CDABS( A(2,N)*Z3(N)    ) + CDABS( A(3,N) ) )
          EPSMAX = MAX( EPSARY( N ) , EPSMAX )
          IF ( EPSARY( N ) .LT. EPS ) THEN
             A(2,N) = ( (CWA*A(0,N)+A(1,N))*CWA+A(2,N))/A(0,N)
             A(1,N) =   (CWA*A(0,N)+A(1,N)) / A(0,N)
             A(0,N) =   CWA
          ENDIF
       ENDIF
    END DO

    IF ( EPSMAX .LT. EPS ) GOTO 800
    I = I + 1
    IF ( I .GE. MAXCNT ) GOTO 600
    GOTO 400

600 IF ( ICHECK .EQ. 0 ) THEN
       ICHECK = 1
       DO N = 1 , NMAX
          IF ( EPSARY( N ) .GE. EPS )  THEN
             W8 = .2D0 * ( CDABS( A( 3 , N )/A( 0 , N ) ) )**(1.D0/3.D0)
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
       END DO
       GOTO 200
    ELSE
       WRITE(6,*) &
         '**** OVER MAX ITERATION COUNT*',MAXCNT,'IN SOBR.W1KSOL ****'
       STOP
    ENDIF

800 CONTINUE
!     WRITE(6,*) '* ITERATION ',I

    DO N = 1 , NMAX
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
       END IF
    END DO

    RETURN
  END SUBROUTINE W1KSOL

!     ******* BAND MATRIX COEFFICIENT *******

  SUBROUTINE W1BNDD(IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS,NX,NWH,NCF,NSF,NSE,NT1,NF
    INTEGER:: I,J,II,JJ,NS0,NS1,NS2,IMODE
    REAL(rkind):: RKV,X1,X2,DX,FLUXB
    COMPLEX(rkind):: CSKXB,CSPXB,CSPZB,CARG,CPHASE

    RKV=2.D6*PI*RF/VC

    NS  = 10
    MWID=3*6 - 1
    MLEN=6*NXPMAX+5
    ALLOCATE(CF(MWID,MLEN))
    NWH = 6 + 6/2

    NCF=0
    NSF=0
    NSE=2
    NT1=5

    NF=NSF-NCF+NWH
    DO I=2-NF,0
       DO J=1,MIN(NF+I-1,NT1)
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=(0.D0,0.D0)
       END DO
    END DO
    DO I=1,NSE
       DO J=1,NT1
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=-CGIN(I,J)
       END DO
    END DO

    NSF=NSF+NSE

    X2  = XA ( 1 )
    NS2 = 6

    DO NX=1,NXPMAX
       X1=X2
       X2=XA(NX+1)
       DX  = RKV * ( X2-X1 )
       NS0 = NT1
       NS1 = NS2
       NS2 = 6
       NSE = NS1
       NT1 = NS1
       IF( NX .EQ. NXPMAX ) NT1 = 5

       NF=NSF-NCF+NWH
       DO I=1,NSE
          DO J=1,NS0
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=CF(J,NS+I)
          END DO
       END DO
       DO I=NSE+1,MWID+NS0-NF
          DO J=MAX(I-MWID+NF,1),NS0
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=(0.D0,0.D0)
          END DO
       END DO

       IF(NX.EQ.ABS(NXABS)) THEN
          DO I=1,NSE
             IMODE=(I-1)/2+1
             CSKXB=CSKX((NX-1)*3+IMODE)
             CSPXB=CSPX((NX-1)*3+IMODE)
             CSPZB=CSPZ((NX-1)*3+IMODE)
             IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
             IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
             FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB &
                   -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB &
                          +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB) &
                          +CD2(3,NX) &
                          +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
             IF((FLUXB.LT.0.D0).AND. &
                ((NXABS.GT.0).OR.(DBLE(CSKXB).GT.0.D0))) THEN
                DO J=1,NS0
                   II=NF+I-J
                   JJ=NCF+J
                   CF(II,JJ)=0.D0
                END DO
             END IF
          END DO
       END IF

       NCF=NCF+NS0
       NF=NSF-NCF+NWH
       DO I=2-NF,0
          DO J=1,MIN(NF+I-1,NT1)
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=(0.D0,0.D0)
          END DO
       END DO
       DO I=1,NSE
!       ....................................................
!           CPHASE=-CDEXP(DCMPLX(0.D0,DX)*CSKX(NS+I-10))
          IMODE = ( I-1 ) / 2 + 1
          CARG=DCMPLX(0.D0,DX)*CSKX((NX-1)*3+IMODE)
          IF ( MOD(I,2) .EQ. 0 )  THEN
             CPHASE=-CDEXP(-CARG)
          ELSE
             CPHASE=-CDEXP( CARG)
          ENDIF
!       ....................................................
          DO J=1,NT1
             II=NF+I-J
             JJ=NCF+J
             CF(II,JJ)=CPHASE*CF(J,NS+I)
          END DO
       END DO

       IF(NX.EQ.ABS(NXABS)) THEN
          DO I=1,NSE
             IMODE=(I-1)/2+1
             CSKXB=CSKX((NX-1)*3+IMODE)
             CSPXB=CSPX((NX-1)*3+IMODE)
             CSPZB=CSPZ((NX-1)*3+IMODE)
             IF(MOD(I,2).EQ.0) CSKXB=-CSKXB
             IF(MOD(I,2).EQ.0) CSPZB=-CSPZB
             FLUXB=-CD1(2,NX)*CSPZB-CD1(1,NX)*DCONJG(CSPZB)*CSPXB &
                   -CSKXB*(CD2(1,NX)*DCONJG(CSPXB)*CSPXB & 
                          +CD2(2,NX)*(DCONJG(CSPXB)-CSPXB) &
                          +CD2(3,NX) &
                          +CD2(4,NX)*DCONJG(CSPZB)*CSPZB)
             IF((FLUXB.GT.0.D0).AND. &
                ((NXABS.GT.0).OR.(DBLE(CSKXB).LT.0.D0))) THEN
                DO J=1,NT1
                   II=NF+I-J
                   JJ=NCF+J
                   CF(II,JJ)=0.D0
                END DO
             END IF
          END DO
       END IF

       NS=NS+NS1
       NSF=NSF+NSE

    END DO

    NS0=NT1
    NSE=2
    NF=NSF-NCF+NWH
    DO I=1,NSE
       DO J=1,NS0
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=CGOT(I,J)
       END DO
    END DO
    DO I=NSE+1,MWID+NS0-NF
       DO J=MAX(I-MWID+NF,1),NS0
          II=NF+I-J
          JJ=NCF+J
          CF(II,JJ)=(0.D0,0.D0)
       END DO
    END DO

    NCF=NCF+NS0
    NSF=NSF+NSE
    IF(NCF.NE.NSF) GO TO 9300

    DO J=1,5
       CA(J)=CGIN(3,J)
       CA(NSF-5+J)=-CGOT(3,J)
    END DO
    DO J=6,NSF-5
       CA(J)=(0.D0,0.D0)
    END DO

    CALL BANDCD(CF,CA,MLEN,MWID,MWID,IERR)
    IF(IERR.NE.0) WRITE(6,601) IERR
    DEALLOCATE(CF)
    RETURN

9300 CONTINUE
    WRITE(6,604) NSF,NCF
    IERR=1
    RETURN

601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
604 FORMAT(1H ,'!! ERROR IN CLBAND : NSF,NCF = ',2I5)
  END SUBROUTINE W1BNDD

!     ******* ELECTROMAGNETIC FIELD IN PLASMA *******

  SUBROUTINE W1EPWD(NZ)
    USE w1comm 
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CEX,CEY,CEZ,CDX,CDY,CDZ
    INTEGER:: NX,NS,K,KK
    REAL(rkind):: RKV,RCE,DX,PABS0,PABS1,PABS2,PABSL,FLUX1,FLUX2

    COMPLEX(rkind):: CA1,CA2,CA3,CA4,CA5,CA6,CPH1,CPH2,CPH3,CPH4,CPH5,CPH6
    COMPLEX(rkind):: CSKX1,CSKX2,CSKX3,CSPX1,CSPX2,CSPX3,CSPZ1,CSPZ2,CSPZ3

    ALLOCATE(CEX(NXPMAX),CEY(NXPMAX),CEZ(NXPMAX))
    ALLOCATE(CDX(NXPMAX),CDY(NXPMAX),CDZ(NXPMAX))

    RKV=2.D6*PI*RF/VC
    RCE=VC*EPS0

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABS(NX,NS)=0.D0
       END DO
    END DO
    DO NX=1,NXPMAX
       FLUX(NX)=0.D0
    END DO

    DO NX = 1 , NXPMAX
       DX = RKV * ( XA( NX+1 ) - XA( NX ) )
       K  = (NX-1)*6 + 2
       KK = (NX-1)*3

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

       CPH1  = CDEXP( 0.5D0*DX*CSKX1 )
       CPH2  = CDEXP(-0.5D0*DX*CSKX1 )
       CPH3  = CDEXP( 0.5D0*DX*CSKX2 )
       CPH4  = CDEXP(-0.5D0*DX*CSKX2 )
       CPH5  = CDEXP( 0.5D0*DX*CSKX3 )
       CPH6  = CDEXP(-0.5D0*DX*CSKX3 )

       CEX(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPX1 & 
                +(CA3*CPH3 +CA4*CPH4 )*CSPX2 &
                +(CA5*CPH5 +CA6*CPH6 )*CSPX3
       CEY(NX) = (CA1*CPH1 +CA2*CPH2 ) &
                +(CA3*CPH3 +CA4*CPH4 ) &
                +(CA5*CPH5 +CA6*CPH6 )
       CEZ(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPZ1 &
                +(CA3*CPH3 -CA4*CPH4 )*CSPZ2 &
                +(CA5*CPH5 -CA6*CPH6 )*CSPZ3
       CDX(NX) = (CA1*CPH1 -CA2*CPH2 )*CSPX1*CSKX1 &
                +(CA3*CPH3 -CA4*CPH4 )*CSPX2*CSKX2 &
                +(CA5*CPH5 -CA6*CPH6 )*CSPX3*CSKX3
       CDY(NX) = (CA1*CPH1 -CA2*CPH2 )      *CSKX1 &
                +(CA3*CPH3 -CA4*CPH4 )      *CSKX2 &
                +(CA5*CPH5 -CA6*CPH6 )      *CSKX3
       CDZ(NX) = (CA1*CPH1 +CA2*CPH2 )*CSPZ1*CSKX1 &
                +(CA3*CPH3 +CA4*CPH4 )*CSPZ2*CSKX2 &
                +(CA5*CPH5 +CA6*CPH6 )*CSPZ3*CSKX3

       CE2DA(NZ,NX,1) = CEX(NX)
       CE2DA(NZ,NX,2) = CEY(NX)
       CE2DA(NZ,NX,3) = CEZ(NX)
    END DO

    DO NS=1,NSMAX
       DO NX = 1 , NXPMAX
          DX = XA( NX+1 ) - XA( NX )

          PABS0 =(CM0(1,NX,NS)* DCONJG(CEX(NX))*CEX(NX) &
                 +CM0(2,NX,NS)* DCONJG(CEX(NX))*CEY(NX) &
                 +CM0(3,NX,NS)* DCONJG(CEY(NX))*CEY(NX) &
                 -CM0(2,NX,NS)* DCONJG(CEY(NX))*CEX(NX) &
                 +CM0(4,NX,NS)* DCONJG(CEZ(NX))*CEZ(NX))*CI
          PABS1 =(CM1(1,NX,NS)*( DCONJG(CEX(NX))*CDZ(NX)) &
                 +CM1(2,NX,NS)*(-DCONJG(CDY(NX))*CEZ(NX)) &
                 +CM1(1,NX,NS)*(-DCONJG(CDZ(NX))*CEX(NX)) &
                 -CM1(2,NX,NS)*( DCONJG(CEZ(NX))*CDY(NX)))
          PABS2 =(CM2(1,NX,NS)* DCONJG(CDX(NX))*CDX(NX) &
                 +CM2(2,NX,NS)* DCONJG(CDX(NX))*CDY(NX) &
                 +CM2(3,NX,NS)* DCONJG(CDY(NX))*CDY(NX) &
                 -CM2(2,NX,NS)* DCONJG(CDY(NX))*CDX(NX) &
                 +CM2(4,NX,NS)* DCONJG(CDZ(NX))*CDZ(NX))*CI

          PABSL=-RCE*RKV*DX*(PABS0 + PABS1 + PABS2)
          PABS(NX,NS) = PABS(NX,NS) + PABSL
       END DO
    END DO

    DO NX=1,NXPMAX
       FLUX1 =(    +  CD1(2,NX)* DCONJG(CEY(NX))*CEZ(NX) &
                   +  CD1(1,NX)* DCONJG(CEZ(NX))*CEX(NX) )*(-1)
       FLUX2 =((CD2(1,NX))* DCONJG(CEX(NX))*CDX(NX) &
              +(CD2(2,NX))* DCONJG(CEX(NX))*CDY(NX) &
              +(CD2(3,NX))* DCONJG(CEY(NX))*CDY(NX) &
              -(CD2(2,NX))* DCONJG(CEY(NX))*CDX(NX) &
              +(CD2(4,NX))* DCONJG(CEZ(NX))*CDZ(NX))*CI

       FLUX(NX)=FLUX(NX)+RCE*(FLUX1+FLUX2)
    END DO
    RETURN
  END SUBROUTINE W1EPWD
END MODULE w1mlm
