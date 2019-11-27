C     $Id$
C
C     ********** SETUP CALCULATION**********
C
      SUBROUTINE WFSTUP
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFSETZ
      CALL WFFEPI
      CALL WFVNOD
C
      IF(MODELN.NE.0) THEN
         IPA=NINT(PA(2))
         IF(IPA.EQ.1) THEN
            CALL ATINIT('H2  ',1)
         ELSEIF(IPA.EQ.40) THEN
            CALL ATINIT('AR  ',1)
         ELSEIF(IPA.EQ.88) THEN
            CALL ATINIT('CF4 ',1)
         ELSE
            WRITE(6,*) 'XX WFSTUP: UNKNOWN ION SPECIES: IPA =',IPA
            CALL ATINIT('H2  ',1)
         ENDIF
      ENDIF
      RETURN
      END
C
C     ****** SET VNOD ******
C
      SUBROUTINE WFVNOD
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION XE(3),YE(3),RL(3)
      REAL*8 A(3),B(3),C(3)
C
      DO IN=1,NNOD
         VNOD(IN)=0.D0
      ENDDO
C
      DO IEDO=1,NELM
         IE=IEDO
         CALL WFNPOS(IE,XE,YE)
         CALL WFABC(IE,A,B,C,S)
         DO K=1,3
            IF(MODELS.EQ.1) THEN
               RL(K)=2*PI*XE(K)
            ELSEIF(MODELS.EQ.2) THEN
               RL(K)=2*PI*(RR+XE(K))
            ELSE
               RL(K)=1.D0
            ENDIF
            DO N=1,3
               IN=IELM(N,IE)
               VNOD(IN)=VNOD(IN)+RL(K)*S*AIF2(N,K)
            ENDDO
         ENDDO
      ENDDO
C      DO IN=1,20
C         WRITE(6,'(A,I5,1PE12.4)') 'IN,VNOD=',IN,VNOD(IN)
C      ENDDO
      RETURN
      END
C
C     ******* MODIFY ANTENNA DATA *******
C
      SUBROUTINE MODANT(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFSETZ
      CALL WFFEPI
C
      DO 5 NA=1,NAMAX
         CALL WFFEP(XJ0(1,NA),YJ0(1,NA),IES)
         IF(IES.EQ.0) THEN
            IF(JNUM0(NA).EQ.1) GOTO 8500
            NJ1=IBDYW(NBDYW)
            DO 10 IBDO=1,NBDYW
               IB=IBDO
               NK1=IBDYW(IB)
               CALL CROS(XJ0(1,NA),YJ0(1,NA),
     &                   XJ0(2,NA),YJ0(2,NA),
     &                   XD(NJ1), YD(NJ1),
     &                   XD(NK1), YD(NK1),XC,YC,IERR)
C               WRITE(6,'(3I8,1P4E12.4/1P6E12.4)') 
C     &                   NJ1,NK1,IERR,XJ0(1,NA),YJ0(1,NA),
C     &                   XJ0(2,NA),YJ0(2,NA),
C     &                   XD(NJ1),YD(NJ1),
C     &                   XD(NK1),YD(NK1),XC,YC
               IF(IERR.EQ.0) GOTO 1000
               NJ1=NK1
   10       CONTINUE
            GOTO 8000
 1000       CALL EFINDL(IES,NJ1,NK1,IEN)
            IF(IEN.EQ.0) GOTO 8100
            IES=IEN
         ELSE
            XC=XJ0(1,NA)
            YC=YJ0(1,NA)
            NJ1=0
            NK1=0
         ENDIF
C
         N=1
         XJ(N,NA)=XC
         YJ(N,NA)=YC
         JAELM(N,NA)=IES
C         WRITE(6,'(A,2I5,1P2E12.4)') '1: N,NE,XJ,YJ=',
C     &        N,IES,XJ(N,NA),YJ(N,NA)
C
         DO 20 ID=2,JNUM0(NA)
            CALL WFFEP(XJ0(ID,NA),YJ0(ID,NA),IEE)
 3000          N=N+1
               IF(N.GT.JNUMM) GOTO 8200
               IF(IES.EQ.IEE) GOTO 4500
               DO 30 J=1,3
                  K=J+1
                  IF(J.EQ.3) K=1
                  NJ2=IELM(J,IES)
                  NK2=IELM(K,IES)
                  IF(   (NJ1.EQ.NJ2.AND.NK1.EQ.NK2)
     &              .OR.(NJ1.EQ.NK2.AND.NK1.EQ.NJ2)) GOTO 30
                  CALL CROS(XJ(N-1,NA),YJ(N-1,NA),
     &                      XJ0(ID,NA),YJ0(ID,NA),
     &                      XD(NJ2),YD(NJ2),
     &                      XD(NK2),YD(NK2),XC,YC,IERR)
                  IF(IERR.EQ.0) GOTO 4000
   30          CONTINUE
C               IF(JNUM0(NA).EQ.2) GOTO 4500
               GOTO 8300
C
 4000          XJ(N,NA)=XC
               YJ(N,NA)=YC
               JAELM(N,NA)=IES
C               WRITE(6,'(A,2I5,1P2E12.4)') '2: N,NE,XJ,YJ=',
C     &              N,IES,XJ(N,NA),YJ(N,NA)
               NJ1=NJ2
               NK1=NK2
C
               CALL EFINDL(IES,NJ1,NK1,IEN)
               IES=IEN
            IF(IEN.NE.0) GOTO 3000
            IF(ID.EQ.JNUM0(NA).AND.IEE.EQ.0) GOTO 6000
            GOTO 8400
 4500       XJ(N,NA)=XJ0(ID,NA)
            YJ(N,NA)=YJ0(ID,NA)
            JAELM(N,NA)=IES
C            WRITE(6,'(A,2I5,1P2E12.4)') '3: N,NE,XJ,YJ=',
C     &             N,IES,XJ(N,NA),YJ(N,NA)
            NJ1=0
            NK1=0
   20    CONTINUE
C
 6000     JNUM(NA)=N
    5 CONTINUE
      IERR=0
      RETURN
C
 8000 IERR=8000
      WRITE(6,800) IERR
  800 FORMAT(1H ,'## MODANT ERROR : IERR =',I5/
     &       1H ,'           : CAN NOT FIND BOUNDARY POINT')
      RETURN
C
 8100 IERR=8100
      IBS=IB-1
      IF(IBS.EQ.0) IBS=NBDYW
      WRITE(6,810) IERR,IB,IBDYW(IBS),IBDYW(IB)
  810 FORMAT(1H ,'## MODANT ERROR : IERR = ',I5/
     &       1H ,'           : CAN NOT FIND BOUNDARY ELEMENT'/
     &       1H ,'           : IB,IBDYW(IB-1),IBDYW(IB) =',3I7)
      RETURN
C
 8200 WRITE(6,820) N,JNUMM
  820 FORMAT(1H ,'## MODANT ERROR : N.GT.JNUMM '/
     &       1H ,'           : N,JNUMM = ',2I7)
      IERR=8200
      RETURN
C
 8300 IERR=8300
      WRITE(6,830) IERR,IES,ID-1,ID
  830 FORMAT(1H ,'## MODANT ERROR : IERR =',I5/
     &       1H ,'           : CAN NOT FIND A POINT OF INTERSECTION'/
     &       1H ,'           : IES,ID-1,ID =',3I7)
      RETURN
C
 8400 IERR=8400
      WRITE(6,840) IERR,IES,NJ1,NK1
  840 FORMAT(1H ,'## MODANT ERROR : IERR =',I5/
     &       1H ,'           : CAN NOT FIND ELEMENT INCLUDING'/
     &           ' TWO NODES NJ1,NK1 '/
     &       1H ,'           : IES, NJ1, NK1 = ',3I7)
      IERR=8400
      RETURN
C
 8500 IERR=8500
      WRITE(6,850) IERR,NA,IES,JNUM0(NA)
  850 FORMAT(1H ,'## MODANT ERROR : IERR = ',I5/
     &       1H ,'           : NA,IES,JNUM0 = ',3I5)
      RETURN
C
      END
C
C     ******* CALCULATE POINT OF INTERSECTION *******
C
      SUBROUTINE CROS(X1,Y1,X2,Y2,X3,Y3,X4,Y4,XC,YC,IERR)
C
      INCLUDE 'wfcomm.inc'
C
      DATA EPS/1.D-12/
C
      X12=X1-X2
      Y12=Y1-Y2
      X34=X3-X4
      Y34=Y3-Y4
      DELT=Y12*X34-X12*Y34
      AD=ABS(DELT)
      IF(AD.LT.EPS) GOTO 9000
C
      RK=(-Y34*(X4-X2)+X34*(Y4-Y2))/DELT
      RT=(-Y12*(X4-X2)+X12*(Y4-Y2))/DELT
      IF(    RK.LT.0.D0.OR.RK.GT.1.D0
     &   .OR.RT.LT.0.D0.OR.RT.GT.1.D0) GOTO 9100
C
      XC=RK*X12+X2
      YC=RK*Y12+Y2
      IERR=0
      RETURN
C
 9000 IERR=9000
      RETURN
C
 9100 IERR=9100
      RETURN
      END
C
C     ******* FIND ELEMENT INCLUDING NODES NJ,NK *******
C
      SUBROUTINE EFINDL(IES,N1,N2,IE)
C
      INCLUDE 'wfcomm.inc'
C
      IF(IES.LT.0.OR.IES.GT.NELM+1) GOTO 9000
      DO 10 I=1,MAX(NELM-IES,IES)
         IDELT=I
         DO 20 J=1,2
            IDELT=-IDELT
            IE    =IES+IDELT
            IF(IE.GE.1.AND.IE.LE.NELM) THEN
               DO 30 K=1,3
                  NE1=IELM(K,IE)
                  IF(NE1.EQ.N1) THEN
                     DO 40 L=1,2
                        LL=MOD(K+L-1,3)+1
                        NE2=IELM(LL,IE)
                        IF(NE2.EQ.N2) RETURN
   40                CONTINUE
                  ENDIF
   30          CONTINUE
            ENDIF
   20    CONTINUE
   10 CONTINUE
C
      IE=0
      RETURN
C
 9000 IE=0
      RETURN
      END
C
C     ******* INITIALIZE FEP *******
C
      SUBROUTINE WFFEPI
C
      INCLUDE 'wfcomm.inc'
      COMMON /WFFEP1/ FMAX(NYM),FMIN(NYM)
      COMMON /WFFEP2/ NEMAX(NYM),NEMIN(NYM),NEYMAX
C
      DIMENSION       XE(3),YE(3)
C
      NEYMAX   =1
      NEMIN(1)=1
      FMAX(1) =-1.D8
      FMIN(1) = 1.D8
      RGX0    =-1.D8
      DO 10 IEDO=1,NELM
         IE=IEDO
         CALL WFNPOS(IE,XE,YE)
         RGX=(XE(1)+XE(2)+XE(3))/3.D0
         IF(RGX.LE.RGX0) THEN
            NEMAX(NEYMAX)=IE-1
            NEYMAX       =NEYMAX+1
            NEMIN(NEYMAX)=IE
            FMAX(NEYMAX) =-1.D8
            FMIN(NEYMAX) = 1.D8
         ENDIF
         FMAX(NEYMAX)=MAX(FMAX(NEYMAX),YE(1),YE(2),YE(3))
         FMIN(NEYMAX)=MIN(FMIN(NEYMAX),YE(1),YE(2),YE(3))
         RGX0=RGX
         IF(NEYMAX.GT.NYM) GOTO 8000
   10 CONTINUE
      NEMAX(NEYMAX)=NELM
      RETURN
C
 8000 WRITE(6,800) NEYMAX,NYM
  800 FORMAT(1H ,'## FEPINT ERROR : NEYMAX,NYM = ',2I5)
      STOP
      END
C
C     ******* FIND ELEMENT INCLUDING POINT (X,Y) *******
C
      SUBROUTINE WFFEP(X,Y,IEL)
C
      INCLUDE 'wfcomm.inc'
      COMMON /WFFEP1/ FMAX(NYM),FMIN(NYM)
      COMMON /WFFEP2/ NEMAX(NYM),NEMIN(NYM),NEYMAX
C
      DIMENSION     XE(3),YE(3)
C
      DO 10 NY=1,NEYMAX
         IF(Y.LE.FMAX(NY).AND.Y.GE.FMIN(NY)) THEN
            KMIN=NEMIN(NY)
            KMAX=NEMAX(NY)
            DO 20 IE=KMIN,KMAX
               IEL=IE
               CALL WFNPOS(IEL,XE,YE)
               IF(Y.GT.MAX(YE(1),YE(2),YE(3)).OR.
     &            Y.LT.MIN(YE(1),YE(2),YE(3)).OR.
     &            X.GT.MAX(XE(1),XE(2),XE(3)).OR.
     &            X.LT.MIN(XE(1),XE(2),XE(3))) GOTO 20
               DO 30 I=1,3
                  J=I+1
                  IF(I.EQ.3) J=1
                  K=I-1
                  IF(I.EQ.1) K=3
                  VALP=((XE(I)-XE(J))*(Y-YE(J))
     &                 -(YE(I)-YE(J))*(X-XE(J)))
     &                *((XE(I)-XE(J))*(YE(K)-YE(J))
     &                 -(YE(I)-YE(J))*(XE(K)-XE(J)))
                  IF(VALP.LT.0.D0) GOTO 20
C                  IF(VALP.EQ.0.D0) WRITE(6,*) '-- WFFEP: ',IE,X,Y
   30          CONTINUE
               RETURN
   20       CONTINUE
         ENDIF
   10 CONTINUE
      IEL=0
C      WRITE(6,*) '-- WFFEP: ',IEL,X,Y
      RETURN
      END
C
C     ******* A,B,C,D,S  CALCULATION *******
C
      SUBROUTINE WFABC(KE,A,B,C,S)
C
      INCLUDE 'wfcomm.inc'
      REAL*8 C
C
      DIMENSION X(3),Y(3),A(3),B(3),C(3)
C
      DO 10 I=1,3
       X(I)=XD(IELM(I,KE))
       Y(I)=YD(IELM(I,KE))
   10 CONTINUE
C
      D=X(1)*(Y(2)-Y(3))+X(2)*(Y(3)-Y(1))+X(3)*(Y(1)-Y(2))
      DO 20 I=1,3
       J=I+1
       K=I+2
       IF(J.GT.3) J=J-3
       IF(K.GT.3) K=K-3
       A(I)=(X(J)*Y(K)-X(K)*Y(J))/D
       B(I)=(Y(J)-Y(K))/D
       C(I)=(X(K)-X(J))/D
   20 CONTINUE
      IF(D.LT.0.D0) WRITE(6,*) 'XX WFABC: ',KE
      S=0.5D0*D
      RETURN
      END
C
C     ******* CALCULATION OF ELEMENT AREA *******
C
      SUBROUTINE WFSELM(KE,S)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION X(3),Y(3)
C
      DO 10 I=1,3
       X(I)=XD(IELM(I,KE))
       Y(I)=YD(IELM(I,KE))
   10 CONTINUE
C
      D=X(1)*(Y(2)-Y(3))+X(2)*(Y(3)-Y(1))+X(3)*(Y(1)-Y(2))
      S=0.5D0*D
      RETURN
      END
C
C     ******* TOTAL COORDINATE - LOCAL COORDINATE *******
C
      SUBROUTINE WFNPOS(IE,XE,YE)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION       XE(3),YE(3)
C
      DO 10 I=1,3
         M=IELM(I,IE)
         XE(I)=XD(M)
         YE(I)=YD(M)
   10 CONTINUE
      RETURN
      END
C
C     ****** ELECTRIC FIELD AT ELEMENT(IE),POINT(X,Y) ******
C
      SUBROUTINE FIELDE(IE,X,Y,CE)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CE(3)
      REAL*8 A(3),B(3),C(3)
C
      IF(MODELS.EQ.1) THEN
         IF(X.LE.0.D0) THEN
            CKZ=0.D0
         ELSE
            CKZ=CI*NPHI/X
         ENDIF
      ELSEIF(MODELS.EQ.2) THEN
         CKZ=CI*NPHI/(RR+X)
      ELSE
         CKZ=CI*RKZ
      ENDIF
      CALL WFABC(IE,A,B,C,S)
      CE(1)=(0.D0,0.D0)
      CE(2)=(0.D0,0.D0)
      CE(3)=(0.D0,0.D0)
      DO 40 I=1,3
         IN=IELM(I,IE)
         WEIGHT=A(I)+B(I)*X+C(I)*Y
         CE(1)=CE(1)+CI*(WEIGHT*CAF(1,IN)+CI*B(I)*CAF(4,IN))
         CE(2)=CE(2)+CI*(WEIGHT*CAF(2,IN)+CI*C(I)*CAF(4,IN))
         CE(3)=CE(3)+CI*(WEIGHT*CAF(3,IN)+CI*CKZ*CAF(4,IN)*WEIGHT)
   40 CONTINUE
C
      IF(MODELS.EQ.1) THEN
         IF(X.LE.0.D0) THEN
            IF(ABS(NPHI).EQ.1) THEN
               CE(3)=CI*NPHI*CE(1)
            ENDIF
          ENDIF
      ENDIF
C
      RETURN
      END
C
C     ****** Potential at ELEMENT(IE),POINT(X,Y) ******
C
      SUBROUTINE FIELDA(X,Y,IE,NZ,ID,X,Y,CA,CAX,CAY,CAZ)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CA(4),CAX(4),CAY(4),CAZ(4)
      REAL*8 A(3),B(3),C(3)
C
      IF(MODELS.EQ.1) THEN
         IF(X.LE.0.D0) THEN
            CKZ=0.D0
         ELSE
            CKZ=CI*NPHIF(NZ)/X
         ENDIF
      ELSEIF(MODELS.EQ.2) THEN
         CKZ=CI*NPHIF(NZ)/(RR+X)
      ELSE
         CKZ=CI*RKZF(NZ)
      ENDIF
      CALL WFABC(IE,A,B,C,S)
      DO J=1,4
         CA(J)=(0.D0,0.D0)
      END DO
      DO I=1,3
         IN=IELM(I,IE)
         WEIGHT=A(I)+B(I)*X+C(I)*Y
         IF(ID.EQ.0) THEN
            DO J=1,4
               CA(J) =CA(J) +WEIGHT*CAFF(NZ,J,IN)
               CAX(J)=CAX(J)+  B(I)*CAFF(NZ,J,IN)
               CAY(J)=CAY(J)+  C(I)*CAFF(NZ,J,IN)
               CAZ(J)=CAZ(J)+   CKZ*CAFF(NZ,J,IN)
            END DO
         ELSE
            DO J=1,4
               CA(J) =CA(J) +WEIGHT*CAFR(NZ,J,IN)
               CAX(J)=CAX(J)+  B(I)*CAFR(NZ,J,IN)
               CAY(J)=CAY(J)+  C(I)*CAFR(NZ,J,IN)
               CAZ(J)=CAZ(J)+   CKZ*CAFR(NZ,J,IN)
            END DO
         END IF
      END DO
C
      RETURN
      END
C
C     ****** FIELD AT ELEMENT(IE),POINT(X,Y) ******
C
      SUBROUTINE FIELDC(IE,X,Y,CZ,IDM,ID,FR,FI)
C
      INCLUDE 'wfcomm.inc'
C
      REAL*8 A(3),B(3),C(3),RLINT(3)
      DIMENSION CZ(IDM,NNODM)
C
      CALL WFABC(IE,A,B,C,S)
      DO 30 I=1,3
         RLINT(I)=A(I)+B(I)*X+C(I)*Y
   30 CONTINUE
      DO 40 J=1,4
         FR=0.D0
         FI=0.D0
         DO 40 I=1,3
            IN=IELM(I,IE)
            CF=CZ(ID,IN)
            FR=FR+RLINT(I)*DBLE(CF)
            FI=FI+RLINT(I)*DIMAG(CF)
   40 CONTINUE
      RETURN
      END
C
C     ****** FIELD AT ELEMENT(IE),POINT(X,Y) ******
C
      SUBROUTINE FIELDD(IE,X,Y,DZ,ID,F)
C
      INCLUDE 'wfcomm.inc'
C
      REAL*8 A(3),B(3),C(3),RLINT(3)
      DIMENSION DZ(NNODM,*)
C
      CALL WFABC(IE,A,B,C,S)
      DO 30 I=1,3
         RLINT(I)=A(I)+B(I)*X+C(I)*Y
   30 CONTINUE
      DO 40 J=1,4
         F=0.D0
         DO 40 I=1,3
            IN=IELM(I,IE)
            DF=DZ(IN,ID)
            F=F+RLINT(I)*DF
   40 CONTINUE
      RETURN
      END
C
C     ******* CEM INITIALIZE *******
C
      SUBROUTINE SETAIF
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION ID3(3,3),ID2(2,2)
      DATA ID3/1,3*0,1,3*0,1/
      DATA ID2/1,0,0,1/
C
      DO 10 I=1,2
         L1=ID2(1,I)
         L2=ID2(2,I)
         AIE1(I)=AI2(L1,L2)
      DO 10 J=1,2
         L1=ID2(1,I)+ID2(1,J)
         L2=ID2(2,I)+ID2(2,J)
         AIE2(I,J)=AI2(L1,L2)
      DO 10 K=1,2
         L1=ID2(1,I)+ID2(1,J)+ID2(1,K)
         L2=ID2(2,I)+ID2(2,J)+ID2(2,K)
         AIE3(I,J,K)=AI2(L1,L2)
   10 CONTINUE
C
      DO 20 I=1,3
         L1=ID3(1,I)
         L2=ID3(2,I)
         L3=ID3(3,I)
         AIF1(I)=AI3(L1,L2,L3)
      DO 20 J=1,3
         L1=ID3(1,I)+ID3(1,J)
         L2=ID3(2,I)+ID3(2,J)
         L3=ID3(3,I)+ID3(3,J)
         AIF2(I,J)=AI3(L1,L2,L3)
      DO 20 K=1,3
         L1=ID3(1,I)+ID3(1,J)+ID3(1,K)
         L2=ID3(2,I)+ID3(2,J)+ID3(2,K)
         L3=ID3(3,I)+ID3(3,J)+ID3(3,K)
         AIF3(I,J,K)=AI3(L1,L2,L3)
   20 CONTINUE
C
      RETURN
      END
C
C     ******* INTEGRAL OF ELEMENT FUNCTION *******
C
      FUNCTION AI2(L1,L2)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C
      DIMENSION       KAI(0:10)
      DATA KAI/1,1,2,6,24,120,720,5040,40320,362880,3628800/
C
      AI2=DBLE(KAI(L1)*KAI(L2))/DBLE(KAI(L1+L2+1))
      RETURN
      END
C
C     ******* INTEGRAL OF ELEMENT FUNCTION *******
C
      FUNCTION AI3(L1,L2,L3)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C
      DIMENSION       KAI(0:10)
      DATA KAI/1,1,2,6,24,120,720,5040,40320,362880,3628800/
C
      AI3=DBLE(2*KAI(L1)*KAI(L2)*KAI(L3))/DBLE(KAI(L1+L2+L3+2))
      RETURN
      END
