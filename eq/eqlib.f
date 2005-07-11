C     $Id$
C
C     ****** SIMPLE RUNGE-KUTTA METHOD ******
C
      SUBROUTINE EQRK4(X,Y,DYDX,YOUT,H,N,DERIVS)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N,NMAX
      DIMENSION DYDX(N),Y(N),YOUT(N)
      EXTERNAL DERIVS
      PARAMETER (NMAX=50)
      INTEGER I
      REAL*8 H6,HH,XH,DYM(NMAX),DYT(NMAX),YT(NMAX)
C
      IF(N.GT.NMAX) THEN
         WRITE(6,*) 'XX EQRK4: N (number of eqs) .GT. NMAX:'
         WRITE(6,*) '   N,NMAX=',N,NMAX
         STOP
      ENDIF
C
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      DO I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
      ENDDO
      CALL DERIVS(XH,YT,DYT)
      DO I=1,N
        YT(I)=Y(I)+HH*DYT(I)
      ENDDO
      CALL DERIVS(XH,YT,DYM)
      DO I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      ENDDO
      CALL DERIVS(X+H,YT,DYT)
      DO I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
      ENDDO
      RETURN
      END
C
C     ****** TWO-DIMENSIONAL NEWTON METHOD ******
C
      SUBROUTINE NEWTN(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,IER)
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SUB
C
      IER=0
      ITER=0
      IF(ABS(X).GT.1.D0) THEN
         HX=DELT*X
      ELSE
         HX=DELT
      ENDIF
      IF(ABS(Y).GT.1.D0) THEN
         HY=DELT*Y
      ELSE
         HY=DELT
      ENDIF
      CALL SUB(X,   Y,   S0,G0)
      IF(LIST.GT.0) WRITE(6,600) X,Y,S0,G0
      IF(LIST.GT.0) WRITE(6,*) 'HX=',HX,'HY=',HY
      IF(LIST.GT.0) CALL GUFLSH
    1 CALL SUB(X+HX,Y,   SX,GX)
      CALL SUB(X,   Y+HY,SY,GY)
      FXX =(SX-S0)/HX
      FXY1=(GX-G0)/HX
      FYY =(GY-G0)/HY
      FXY2=(SY-S0)/HY
      FXY=0.5D0*(FXY1+FXY2)
      IF(LIST.GT.1) WRITE(6,601) SX,SY,GX,GY
      IF(LIST.GT.1) WRITE(6,601) FXX,FYY,FXY1,FXY2
      IF(LIST.GT.1) CALL GUFLSH
      DF=SQRT(S0*S0+G0*G0)
      DET=FXX*FYY-FXY*FXY
      H11= FYY/DET
      H12=-FXY/DET
      H21=-FXY/DET
      H22= FXX/DET
C
      DX=-(H11*S0+H12*G0)
      DY=-(H21*S0+H22*G0)
      TT=1.0D0
    2 X=X+TT*DX
      Y=Y+TT*DY
      CALL SUB(X,   Y,   SN,GN)
      IF(LIST.GT.0) WRITE(6,600) DX,DY
      IF(LIST.GT.0) WRITE(6,600) X,Y,SN,GN
      IF(LIST.GT.0) CALL GUFLSH
      DFN=SQRT(SN*SN+GN*GN)
      IF(DFN.GT.DF) THEN
         X=X-TT*DX
         Y=Y-TT*DY
         TT=0.5D0*TT
         ITER=ITER+1
         IF(TT.LE.1.D-3) GOTO 8000
         IF(ITER.LE.ILMAX) GOTO 2
      ELSE
         S0=SN
         G0=GN
         DF=DFN
         ITER=ITER+1
         IF(DF.LE.EPS) GO TO 9000
         IF(ITER.LE.ILMAX) GO TO 1
      ENDIF
C
      IER=2
      IF(LIST.GT.0)
     &WRITE(6,*) 'XX NEWTN: LOOP COUNT EXCEEDS UPPER BOUND.'
      GOTO 9000
C
 8000 IER=1
      IF(LIST.GT.0)
     &WRITE(6,*) 'XX NEWTN: DOES NOT CONVERGE.'
      GOTO 9000
C
 9000 XX=X
      YY=Y
      RETURN
  600 FORMAT(" ",6X,'X,Y,FX,FY = ',1P4E15.7)
  601 FORMAT(" ",6X,'FXX,YY,XY = ',1P4E15.7)
      END
C
C   ***************************************
C   **         Power function            **
C   **    0.D0**0.D0 should be 1.D0      **
C   ***************************************
C
      FUNCTION FPOW(X,Y)
C
      REAL*8 X,Y,FPOW
C
      IF(X.EQ.0.D0) THEN
         IF(Y.EQ.0.D0) THEN
            FPOW=1.D0
         ELSE
            FPOW=0.D0
         ENDIF
      ELSE
         FPOW=X**Y
      ENDIF
      RETURN
      END
