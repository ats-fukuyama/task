C     $Id$
C   ************************************************
C   **              ƒXƒvƒ‰ƒCƒ“•âŠÔ                **
C   ************************************************
C
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of coefficients ****
C
      SUBROUTINE SPL1D(X,F,FX,U,NXMAX,ID,IERR)
C
C      INPUT : X(NXMAX)  : COORDINATES
C              F(NXMAX)  : VALUE
C              FX(NXMAX) : EDGE DERIVATIVE FOR ID != 0
C              NXMAX     : NUMBER OF VARIABLES (<NMAX=10001)
C              ID        : 0 : SECOND DERIVATIVES = 0 AT X(1) AND X(NXMAX)
C                          1 : DERIVATIVE FX(1) IS GIVEN
C                          2 : DERIVATIVE FX(NXMAX) IS GIVEN
C                          3 : DERIVATIVES FX(1) AND FX(NXMAX) ARE GIVEN
C      OUTPUT: U(4,NXMAX): SPLINE COEFICIENTS
C              FX(NXMAX) : ESTIMATED DERIVATIVES
C              IERR      : ERROR INDICATOR
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=10001)
C
      DIMENSION X(NXMAX),F(NXMAX)
      DIMENSION U(4,NXMAX)
      DIMENSION FX(NXMAX)
      DIMENSION AX(4,NMAX),BX(NMAX)
C
      IERR=0
      IF(NXMAX.GT.NMAX) GOTO 9001
C
      ID1=MOD(ID,2)
      ID2=MOD(ID/2,2)
C
      IF(ID1.EQ.0) THEN
         DXP=X(2)-X(1)
         AX(1,1)=0.D0
         AX(2,1)=2.D0*DXP
         AX(3,1)=DXP
      ELSE
         AX(1,1)=0.D0
         AX(2,1)=1.D0
         AX(3,1)=0.D0
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         AX(1,NX)=DXP
         AX(2,NX)=2.D0*(DXP+DXM)
         AX(3,NX)=DXM
      ENDDO
      IF(ID2.EQ.0) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         AX(1,NXMAX)=DXM
         AX(2,NXMAX)=2.D0*DXM
         AX(3,NXMAX)=0.D0
      ELSE
         AX(1,NXMAX)=0.D0
         AX(2,NXMAX)=1.D0
         AX(3,NXMAX)=0.D0
      ENDIF
C
      IF(ID1.EQ.0) THEN
         BX(1)=3.D0*(F(2)-F(1))
      ELSE
         BX(1)=FX(1)
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         BX(NX)=3.D0*(DXM*(F(NX+1)-F(NX))/DXP
     &               +DXP*(F(NX)-F(NX-1))/DXM)
      ENDDO
      IF(ID2.EQ.0) THEN
         BX(NXMAX)=3.D0*(F(NXMAX)-F(NXMAX-1))
      ELSE
         BX(NXMAX)=FX(NXMAX)
      ENDIF
C
      CALL TDMSRD(AX,BX,NXMAX,FX,IERR)
      IF(IERR.NE.0) GOTO 9003
C
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)
C
         T11=  F(NX-1)
         T21= FX(NX-1)
         T31=  F(NX  )
         T41= FX(NX  )
C
         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         V31=-3.D0*DX2
         V32=-2.D0*DX1
         V33= 3.D0*DX2
         V34=-1.D0*DX1
         V41= 2.D0*DX3
         V42=      DX2
         V43=-2.D0*DX3
         V44=      DX2
C
         U(1,NX)=    T11
         U(2,NX)=            T21
         U(3,NX)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,NX)=V41*T11+V42*T21+V43*T31+V44*T41
      ENDDO
      RETURN
C
 9001 WRITE(6,*) 'XX SPL1D: NXMAX.GT.NMAX:',NXMAX,NMAX
      IERR=1
      RETURN
 9003 WRITE(6,*) 'XX CSPL1D: TDMSRD ERROR : IERR=',IERR
      IERR=3
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of interpolation  ****
C
      SUBROUTINE SPL1DF(X0,F0,X,U,NXMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX)
      DIMENSION U(4,NXMAX)
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
C
 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2
C
      DX=X0-X(NX-1)
C
      F0= U(1,NX)
     &  + U(2,NX)*DX
     &  + U(3,NX)*DX*DX
     &  + U(4,NX)*DX*DX*DX
C      WRITE(6,'(A,2I5,1P3E12.4)') 
C     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of interpolation  ****
C
      SUBROUTINE SPL1DD(X0,F0,DF0,X,U,NXMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX)
      DIMENSION U(4,NXMAX)
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
C
 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2
C
      DX=X0-X(NX-1)
C
      F0= U(1,NX)
     &  + U(2,NX)*DX
     &  + U(3,NX)*DX*DX
     &  + U(4,NX)*DX*DX*DX
      DF0= U(2,NX)
     &   + U(3,NX)*DX*2
     &   + U(4,NX)*DX*DX*3
      IERR=0
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of U0  ****
C
      SUBROUTINE SPL1DI0(X,U,U0,NXMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX)
      DIMENSION U(4,NXMAX),U0(NXMAX)
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
C
      U0(1)=0.D0
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)
C
         U0(NX)= U0(NX-1)
     &         + U(1,NX)*DX
     &         + U(2,NX)*DX*DX/2.D0
     &         + U(3,NX)*DX*DX*DX/3.D0
     &         + U(4,NX)*DX*DX*DX*DX/4.D0
      ENDDO
C
C     WRITE(6,'(A,2I5,1P3E12.4)') 
C     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of interpolation  ****
C
      SUBROUTINE SPL1DI(X0,F0D,X,U,U0,NXMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX)
      DIMENSION U(4,NXMAX),U0(NXMAX)
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
C
 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2
C
      DX=X0-X(NX-1)
C
      F0D= U0(NX)
     &   + U(1,NX)*DX
     &   + U(2,NX)*DX*DX/2.D0
     &   + U(3,NX)*DX*DX*DX/3.D0
     &   + U(4,NX)*DX*DX*DX*DX/4.D0
C      WRITE(6,'(A,2I5,1P3E12.4)') 
C     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END
C
C     ****** Two-Dimensional Spline Interpolation ******
C       **** Calculation of coefficients ****
C
      SUBROUTINE SPL2D(X,Y,F,FX,FY,FXY,U,NXM,NXMAX,NYMAX,
     &                 IDX,IDY,IERR)
C      INPUT : X(NXMAX)        : COORDINATES
C              Y(NYMAX)        : COORDINATES
C              F(NXM,NYMAX)  : VALUE
C              FX(NXM,NYMAX) : EDGE DERIVATIVES FOR IDX != 0
C              FY(NXM,NYMAX) : EDGE DERIVATIVES FOR IDY != 0
C              FXY(NXM,NYMAX): CORNER DERIVATIVES FOR IDY OR IDY != 0
C              NXM       : ARRAY SIZE
C              NXMAX     : NUMBER OF VARIABLES (<NMAX=10001)
C              NYMAX     : NUMBER OF VARIABLES (<NMAX=10001)
C              IDX       : 0 : SECOND DERIVATIVES = 0 AT X(1) AND X(NXMAX)
C                          1 : DERIVATIVE FX(1) IS GIVEN
C                          2 : DERIVATIVE FX(NXMAX) IS GIVEN
C                          3 : DERIVATIVES FX(1) AXD FX(NXMAX) ARE GIVEN
C              IDY       : 0 : SECOND DERIVATIVES = 0 AT Y(1) AND Y(NYMAX)
C                          1 : DERIVATIVE FY(1) IS GIVEN
C                          2 : DERIVATIVE FY(NYMAX) IS GIVEN
C                          3 : DERIVATIVES FY(1) AXD FY(NYMAX) ARE GIVEN
C      OUTPUT: U(4,4,NXM,NYMAX): SPLINE COEFICIENTS
C              FX(NXM,NYMAX) : ESTIMATED DERIVATIVES
C              FY(NXM,NYMAX) : ESTIMATED DERIVATIVES
C              FXY(NXM,NYMAX): ESTIMATED DERIVATIVES
C              IERR          : ERROR INDICATOR
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=10001)
C
      DIMENSION X(NXMAX),Y(NYMAX),F(NXM,NYMAX)
      DIMENSION U(4,4,NXM,NYMAX)
      DIMENSION FX(NXM,NYMAX),FY(NXM,NYMAX),FXY(NXM,NYMAX)
      DIMENSION AX(4,NMAX),AX0(4,NMAX),BX(NMAX)
      DIMENSION AY(4,NMAX),AY0(4,NMAX),BY(NMAX)
      DIMENSION FXY1(4,NMAX),FXY2(4,NMAX)
C
      IF(NXMAX.GT.NMAX) GOTO 9001
      IF(NYMAX.GT.NMAX) GOTO 9002
C
      IDX1=MOD(IDX,2)
      IDX2=MOD(IDX/2,2)
      IDY1=MOD(IDY,2)
      IDY2=MOD(IDY/2,2)
C
      IF(IDX1.EQ.0) THEN
         DXP=X(2)-X(1)
         AX(1,1)=0.D0
         AX(2,1)=2.D0*DXP
         AX(3,1)=DXP
      ELSE
         AX(1,1)=0.D0
         AX(2,1)=1.D0
         AX(3,1)=0.D0
      ENDIF
      DO 1000 NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         AX(1,NX)=DXP
         AX(2,NX)=2.D0*(DXP+DXM)
         AX(3,NX)=DXM
 1000 CONTINUE
      IF(IDX2.EQ.0) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         AX(1,NXMAX)=DXM
         AX(2,NXMAX)=2.D0*DXM
         AX(3,NXMAX)=0.D0
      ELSE
         AX(1,NXMAX)=0.D0
         AX(2,NXMAX)=1.D0
         AX(3,NXMAX)=0.D0
      ENDIF
      DO NX=1,NXMAX
         AX0(1,NX)=AX(1,NX)
         AX0(2,NX)=AX(2,NX)
         AX0(3,NX)=AX(3,NX)
      ENDDO
C
      IF(IDY1.EQ.0) THEN
         DYP=Y(2)-Y(1)
         AY(1,1)=0.D0
         AY(2,1)=2.D0*DYP
         AY(3,1)=DYP
      ELSE
         AY(1,1)=0.D0
         AY(2,1)=1.D0
         AY(3,1)=0.D0
      ENDIF
      DO NY=2,NYMAX-1
         DYM=Y(NY)-Y(NY-1)
         DYP=Y(NY+1)-Y(NY)
         AY(1,NY)=DYP
         AY(2,NY)=2.D0*(DYP+DYM)
         AY(3,NY)=DYM
      ENDDO
      IF(IDY2.EQ.0) THEN
         DYM=Y(NYMAX)-Y(NYMAX-1)
         AY(1,NYMAX)=DYM
         AY(2,NYMAX)=2.D0*DYM
         AY(3,NYMAX)=0.D0
      ELSE
         AY(1,NYMAX)=0.D0
         AY(2,NYMAX)=1.D0
         AY(3,NYMAX)=0.D0
      ENDIF
      DO NY=1,NYMAX
         AY0(1,NY)=AY(1,NY)
         AY0(2,NY)=AY(2,NY)
         AY0(3,NY)=AY(3,NY)
      ENDDO
C
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            AX(1,NX)=AX0(1,NX)
            AX(2,NX)=AX0(2,NX)
            AX(3,NX)=AX0(3,NX)
         ENDDO
         IF(IDX1.EQ.0) THEN
            BX(1)=3.D0*(F(2,NY)-F(1,NY))
         ELSE
            BX(1)=FX(1,NY)
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(F(NX+1,NY)-F(NX,NY))/DXP
     &                  +DXP*(F(NX,NY)-F(NX-1,NY))/DXM)
         ENDDO
         IF(IDX2.EQ.0) THEN
            BX(NXMAX)=3.D0*(F(NXMAX,NY)-F(NXMAX-1,NY))
         ELSE
            BX(NXMAX)=FX(NXMAX,NY)
         ENDIF
         CALL TDMSRD(AX,BX,NXMAX,FX,IERR)
         IF(IERR.NE.0) GOTO 9003
      ENDDO
C
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            AY(1,NY)=AY0(1,NY)
            AY(2,NY)=AY0(2,NY)
            AY(3,NY)=AY0(3,NY)
         ENDDO
         IF(IDY1.EQ.0) THEN
            BY(1)=3.D0*(F(NX,2)-F(NX,1))
         ELSE
            BY(1)=FY(NX,1)
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(F(NX,NY+1)-F(NX,NY))/DYP
     &                  +DYP*(F(NX,NY)-F(NX,NY-1))/DYM)
         ENDDO
         IF(IDY2.EQ.0) THEN
            BY(NYMAX)=3.D0*(F(NX,NYMAX)-F(NX,NYMAX-1))
         ELSE
            BY(NYMAX)=FY(NX,NYMAX)
         ENDIF
         CALL TDMSRD(AY,BY,NYMAX,FY,IERR)
         IF(IERR.NE.0) GOTO 9004
      ENDDO
C
      DO NY=1,NYMAX,NYMAX-1
         DO NX=1,NXMAX
            AX(1,NX)=AX0(1,NX)
            AX(2,NX)=AX0(2,NX)
            AX(3,NX)=AX0(3,NX)
         ENDDO
         IF(IDX1.EQ.0) THEN
            BX(1)=3.D0*(FY(2,NY)-FY(1,NY))
         ELSE
            BX(1)=FXY(1,NY)
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FY(NX+1,NY)-FY(NX,NY))/DXP
     &                  +DXP*(FY(NX,NY)-FY(NX-1,NY))/DXM)
         ENDDO
         IF(IDX2.EQ.0) THEN
            BX(NXMAX)=3.D0*(FY(NXMAX,NY)-FY(NXMAX-1,NY))
         ELSE
            BX(NXMAX)=FXY(NXMAX,NY)
         ENDIF
         CALL TDMSRD(AX,BX,NXMAX,FXY1,IERR)
         IF(IERR.NE.0) GOTO 9003
      ENDDO
C
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            AY(1,NY)=AY0(1,NY)
            AY(2,NY)=AY0(2,NY)
            AY(3,NY)=AY0(3,NY)
         ENDDO
         IF(IDY1.EQ.0) THEN
            BY(1)=3.D0*(FX(NX,2)-FX(NX,1))
         ELSE
            BY(1)=FXY(NX,1)
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(FX(NX,NY+1)-FX(NX,NY))/DYP
     &                  +DYP*(FX(NX,NY)-FX(NX,NY-1))/DYM)
         ENDDO
         IF(IDY2.EQ.0) THEN
            BY(NYMAX)=3.D0*(FX(NX,NYMAX)-FX(NX,NYMAX-1))
         ELSE
            BY(NYMAX)=FXY(NX,NYMAX)
         ENDIF
         CALL TDMSRD(AY,BY,NYMAX,FXY2,IERR)
         IF(IERR.NE.0) GOTO 9004
      ENDDO
C
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            FXY(NX,NY)=0.5D0*(FXY1(NX,NY)+FXY2(NX,NY))
C            FXY(NX,NY)=FXY1(NX,NY)
C            FXY(NX,NY)=FXY2(NX,NY)
         ENDDO
      ENDDO
C
      DO NX=2,NXMAX
      DO NY=2,NYMAX
         DX=X(NX)-X(NX-1)
         DY=Y(NY)-Y(NY-1)
C
         DY1=1.D0/DY
         DY2=DY1*DY1
         DY3=DY2*DY1
         V13=-3.D0*DY2
         V23=-2.D0*DY1
         V33= 3.D0*DY2
         V43=-1.D0*DY1
         V14= 2.D0*DY3
         V24=      DY2
         V34=-2.D0*DY3
         V44=      DY2
C
         T11=  F(NX-1,NY-1)
         T12=                    FY(NX-1,NY-1)
         T13=  F(NX-1,NY-1)*V13+ FY(NX-1,NY-1)*V23
     &       + F(NX-1,NY  )*V33+ FY(NX-1,NY  )*V43
         T14=  F(NX-1,NY-1)*V14+ FY(NX-1,NY-1)*V24
     &       + F(NX-1,NY  )*V34+ FY(NX-1,NY  )*V44
         T21= FX(NX-1,NY-1)
         T22=                   FXY(NX-1,NY-1)
         T23= FX(NX-1,NY-1)*V13+FXY(NX-1,NY-1)*V23
     &       +FX(NX-1,NY  )*V33+FXY(NX-1,NY  )*V43
         T24= FX(NX-1,NY-1)*V14+FXY(NX-1,NY-1)*V24
     &       +FX(NX-1,NY  )*V34+FXY(NX-1,NY  )*V44
         T31=  F(NX,  NY-1)
         T32=                    FY(NX,  NY-1)
         T33=  F(NX,  NY-1)*V13+ FY(NX,  NY-1)*V23
     &       + F(NX,  NY  )*V33+ FY(NX,  NY  )*V43
         T34=  F(NX,  NY-1)*V14+ FY(NX,  NY-1)*V24
     &       + F(NX,  NY  )*V34+ FY(NX,  NY  )*V44
         T41= FX(NX,  NY-1)
         T42=                   FXY(NX,  NY-1)
         T43= FX(NX,  NY-1)*V13+FXY(NX,  NY-1)*V23
     &       +FX(NX,  NY  )*V33+FXY(NX,  NY  )*V43
         T44= FX(NX,  NY-1)*V14+FXY(NX,  NY-1)*V24
     &       +FX(NX,  NY  )*V34+FXY(NX,  NY  )*V44
C
         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         V31=-3.D0*DX2
         V32=-2.D0*DX1
         V33= 3.D0*DX2
         V34=-1.D0*DX1
         V41= 2.D0*DX3
         V42=      DX2
         V43=-2.D0*DX3
         V44=      DX2
C
         U(1,1,NX,NY)=    T11
         U(2,1,NX,NY)=            T21
         U(3,1,NX,NY)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,1,NX,NY)=V41*T11+V42*T21+V43*T31+V44*T41
         U(1,2,NX,NY)=    T12
         U(2,2,NX,NY)=            T22
         U(3,2,NX,NY)=V31*T12+V32*T22+V33*T32+V34*T42
         U(4,2,NX,NY)=V41*T12+V42*T22+V43*T32+V44*T42
         U(1,3,NX,NY)=    T13
         U(2,3,NX,NY)=            T23
         U(3,3,NX,NY)=V31*T13+V32*T23+V33*T33+V34*T43
         U(4,3,NX,NY)=V41*T13+V42*T23+V43*T33+V44*T43
         U(1,4,NX,NY)=    T14
         U(2,4,NX,NY)=            T24
         U(3,4,NX,NY)=V31*T14+V32*T24+V33*T34+V34*T44
         U(4,4,NX,NY)=V41*T14+V42*T24+V43*T34+V44*T44
      ENDDO
      ENDDO
      IERR=0
      RETURN
C
 9001 WRITE(6,*) 'XX SPL2D: NXMAX.GT.NMAX:',NXMAX,NMAX
      IERR=1
      RETURN
 9002 WRITE(6,*) 'XX SPL2D: NYMAX.GT.NMAX:',NYMAX,NMAX
      IERR=2
      RETURN
 9003 WRITE(6,*) 'XX SPL2D: TDMSRD ERROR: IERR=',IERR
      IERR=3
      RETURN
 9004 WRITE(6,*) 'XX SPL2D: TDMSRD ERROR: IERR=',IERR
      IERR=4
      RETURN
      END
C
C     ****** Two-Dimensional Spline Interpolation ******
C       **** Calculation of interpolation  ****
C
      SUBROUTINE SPL2DF(X0,Y0,F0,X,Y,U,NXM,NXMAX,NYMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX),Y(NYMAX)
      DIMENSION U(4,4,NXM,NYMAX)
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      FX=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FX*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
      FY=1.D0/(Y(NYMAX)-Y(1))
      NY=NINT((Y0-Y(1))*FY*(NYMAX-1))+1
      IF(NY.LT.1) THEN
         IERR=3
         NY=2
      ENDIF
      IF(NY.GT.NYMAX) THEN
         IERR=4
         NY=NYMAX
      ENDIF
C
 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FX.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FX.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2
C
 5005 IF(NY.GE.NYMAX) GOTO 5006
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 5006
         NY=NY+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NY.LE.2) GOTO 5008
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 5008
         NY=NY-1
         GOTO 5007
 5008 CONTINUE
      IF(NY.LT.2)     NY=2
C
      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)
C
C      WRITE(6,'(2I5,1P2E12.4)') NX,NY,X0,Y0
C      IF(NX.EQ.33 .AND. NY.EQ.18) THEN
C         WRITE(6,'(1P4E12.4)') X0,Y0,DX,DY
C         WRITE(6,'(1P4E12.4)') U(1,1,NX,NY),U(1,2,NX,NY),
C     &                         U(1,3,NX,NY),U(1,4,NX,NY)
C         WRITE(6,'(1P4E12.4)') U(2,1,NX,NY),U(2,2,NX,NY),
C     &                         U(2,3,NX,NY),U(2,4,NX,NY)
C         WRITE(6,'(1P4E12.4)') U(3,1,NX,NY),U(3,2,NX,NY),
C     &                         U(3,3,NX,NY),U(3,4,NX,NY)
C         WRITE(6,'(1P4E12.4)') U(4,1,NX,NY),U(4,2,NX,NY),
C     &                         U(4,3,NX,NY),U(4,4,NX,NY)
C         CALL GUFLSH
C      ENDIF
C
C      F0= U(1,1,NX,NY)      +U(1,2,NX,NY)*DY
C     &   +U(1,3,NX,NY)*DY*DY+U(1,4,NX,NY)*DY*DY*DY
C     &  +(U(2,1,NX,NY)      +U(2,2,NX,NY)*DY
C     &   +U(2,3,NX,NY)*DY*DY+U(2,4,NX,NY)*DY*DY*DY)*DX
C     &  +(U(3,1,NX,NY)      +U(3,2,NX,NY)*DY
C     &   +U(3,3,NX,NY)*DY*DY+U(3,4,NX,NY)*DY*DY*DY)*DX*DX
C     &  +(U(4,1,NX,NY)      +U(4,2,NX,NY)*DY
C     &   +U(4,3,NX,NY)*DY*DY+U(4,4,NX,NY)*DY*DY*DY)*DX*DX*DX
C
      F4= ((U(4,4,NX,NY) *DY     +U(4,3,NX,NY))*DY
     &     +U(4,2,NX,NY))*DY     +U(4,1,NX,NY)
      F3= ((U(3,4,NX,NY) *DY     +U(3,3,NX,NY))*DY
     &     +U(3,2,NX,NY))*DY     +U(3,1,NX,NY)
      F2= ((U(2,4,NX,NY) *DY     +U(2,3,NX,NY))*DY
     &     +U(2,2,NX,NY))*DY     +U(2,1,NX,NY)
      F1= ((U(1,4,NX,NY) *DY     +U(1,3,NX,NY))*DY
     &     +U(1,2,NX,NY))*DY     +U(1,1,NX,NY)
      F0= ((F4*DX+F3)*DX+F2)*DX+F1
C
C      F0=((( ((U(4,4,NX,NY) *DY     +U(4,3,NX,NY))*DY
C     &        +U(4,2,NX,NY))*DY     +U(4,1,NX,NY)    )*DX
C     &      +((U(3,4,NX,NY) *DY     +U(3,3,NX,NY))*DY
C     &        +U(3,2,NX,NY))*DY     +U(3,1,NX,NY)    )*DX
C     &      +((U(2,4,NX,NY) *DY     +U(2,3,NX,NY))*DY
C     &        +U(2,2,NX,NY))*DY     +U(2,1,NX,NY)    )*DX
C     &      +((U(1,4,NX,NY) *DY     +U(1,3,NX,NY))*DY
C     &        +U(1,2,NX,NY))*DY     +U(1,1,NX,NY)
C
      IERR=0
      RETURN
      END
C
      SUBROUTINE SPL2DD(X0,Y0,F0,FX0,FY0,X,Y,U,NXM,NXMAX,NYMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX),Y(NYMAX)
      DIMENSION U(4,4,NXM,NYMAX)
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      FX=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FX*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
      FY=1.D0/(Y(NYMAX)-Y(1))
      NY=NINT((Y0-Y(1))*FY*(NYMAX-1))+1
      IF(NY.LT.1) THEN
         IERR=3
         NY=2
      ENDIF
      IF(NY.GT.NYMAX) THEN
         IERR=4
         NY=NYMAX
      ENDIF
C
 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FX.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FX.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2
C
 5005 IF(NY.GE.NYMAX) GOTO 5006
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 5006
         NY=NY+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NY.LE.2) GOTO 5008
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 5008
         NY=NY-1
         GOTO 5007
 5008 CONTINUE
      IF(NY.LT.2)     NY=2
C
      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)
C
C      F0= U(1,1,NX,NY)      +U(1,2,NX,NY)*DY
C     &   +U(1,3,NX,NY)*DY*DY+U(1,4,NX,NY)*DY*DY*DY
C     &  +(U(2,1,NX,NY)      +U(2,2,NX,NY)*DY
C     &   +U(2,3,NX,NY)*DY*DY+U(2,4,NX,NY)*DY*DY*DY)*DX
C     &  +(U(3,1,NX,NY)      +U(3,2,NX,NY)*DY
C     &   +U(3,3,NX,NY)*DY*DY+U(3,4,NX,NY)*DY*DY*DY)*DX*DX
C     &  +(U(4,1,NX,NY)      +U(4,2,NX,NY)*DY
C     &   +U(4,3,NX,NY)*DY*DY+U(4,4,NX,NY)*DY*DY*DY)*DX*DX*DX
C
      F4= ((U(4,4,NX,NY) *DY     +U(4,3,NX,NY))*DY
     &     +U(4,2,NX,NY))*DY     +U(4,1,NX,NY)
      F3= ((U(3,4,NX,NY) *DY     +U(3,3,NX,NY))*DY
     &     +U(3,2,NX,NY))*DY     +U(3,1,NX,NY)
      F2= ((U(2,4,NX,NY) *DY     +U(2,3,NX,NY))*DY
     &     +U(2,2,NX,NY))*DY     +U(2,1,NX,NY)
      F1= ((U(1,4,NX,NY) *DY     +U(1,3,NX,NY))*DY
     &     +U(1,2,NX,NY))*DY     +U(1,1,NX,NY)
      F0= ((F4*DX+F3)*DX+F2)*DX+F1
C
C      FX0=(U(2,1,NX,NY)      +U(2,2,NX,NY)*DY
C     &    +U(2,3,NX,NY)*DY*DY+U(2,4,NX,NY)*DY*DY*DY)
C     &   +(U(3,1,NX,NY)      +U(3,2,NX,NY)*DY
C     &    +U(3,3,NX,NY)*DY*DY+U(3,4,NX,NY)*DY*DY*DY)*DX*2
C     &   +(U(4,1,NX,NY)      +U(4,2,NX,NY)*DY
C     &    +U(4,3,NX,NY)*DY*DY+U(4,4,NX,NY)*DY*DY*DY)*DX*DX*3
C
      FX0= (3*F4*DX+2*F3)*DX+F2
C
C      FY0=                   +U(1,2,NX,NY)
C     &    +U(1,3,NX,NY)*DY*2 +U(1,4,NX,NY)*DY*DY*3
C     &   +(                  +U(2,2,NX,NY)
C     &    +U(2,3,NX,NY)*DY*2 +U(2,4,NX,NY)*DY*DY*3 )*DX
C     &   +(                  +U(3,2,NX,NY)
C     &    +U(3,3,NX,NY)*DY*2 +U(3,4,NX,NY)*DY*DY*3 )*DX*DX
C     &   +(                  +U(4,2,NX,NY)
C     &    +U(4,3,NX,NY)*DY*2 +U(4,4,NX,NY)*DY*DY*3 )*DX*DX*DX
C
      FY4= (3*U(4,4,NX,NY) *DY     +2*U(4,3,NX,NY))*DY
     &       +U(4,2,NX,NY)
      FY3= (3*U(3,4,NX,NY) *DY     +2*U(3,3,NX,NY))*DY
     &       +U(3,2,NX,NY)
      FY2= (3*U(2,4,NX,NY) *DY     +2*U(2,3,NX,NY))*DY
     &       +U(2,2,NX,NY)
      FY1= (3*U(1,4,NX,NY) *DY     +2*U(1,3,NX,NY))*DY
     &       +U(1,2,NX,NY)
      FY0= ((FY4*DX+FY3)*DX+FY2)*DX+FY1
      IERR=0
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of coefficients ****
C
      SUBROUTINE CSPL1D(X,F,FX,U,NXMAX,ID,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=10001)
C
      DIMENSION X(NXMAX),AX(4,NMAX)
      COMPLEX*16 F(NXMAX),U(4,NXMAX),FX(NXMAX),BX(NMAX)
      COMPLEX*16 T11,T21,T31,T41
C
      IF(NXMAX.GT.NMAX) GOTO 9001
C
      ID1=MOD(ID,2)
      ID2=MOD(ID/2,2)
C
      IF(ID1.EQ.0) THEN
         DXP=X(2)-X(1)
         AX(1,1)=0.D0
         AX(2,1)=2.D0*DXP
         AX(3,1)=DXP
      ELSE
         AX(1,1)=0.D0
         AX(2,1)=1.D0
         AX(3,1)=0.D0
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         AX(1,NX)=DXP
         AX(2,NX)=2.D0*(DXP+DXM)
         AX(3,NX)=DXM
      ENDDO
      IF(ID2.EQ.0) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         AX(1,NXMAX)=DXM
         AX(2,NXMAX)=2.D0*DXM
         AX(3,NXMAX)=0.D0
      ELSE
         AX(1,NXMAX)=0.D0
         AX(2,NXMAX)=1.D0
         AX(3,NXMAX)=0.D0
      ENDIF
C
      IF(ID1.EQ.0) THEN
         BX(1)=3.D0*(F(2)-F(1))
      ELSE
         BX(1)=FX(1)
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         BX(NX)=3.D0*(DXM*(F(NX+1)-F(NX))/DXP
     &               +DXP*(F(NX)-F(NX-1))/DXM)
      ENDDO
      IF(ID2.EQ.0) THEN
         BX(NXMAX)=3.D0*(F(NXMAX)-F(NXMAX-1))
      ELSE
         BX(NXMAX)=FX(NXMAX)
      ENDIF
      CALL TDMSCD(AX,BX,NXMAX,FX,IERR)
      IF(IERR.NE.0) GOTO 9003
C
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)
C
         T11=  F(NX-1)
         T21= FX(NX-1)
         T31=  F(NX  )
         T41= FX(NX  )
C
         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         V31=-3.D0*DX2
         V32=-2.D0*DX1
         V33= 3.D0*DX2
         V34=-1.D0*DX1
         V41= 2.D0*DX3
         V42=      DX2
         V43=-2.D0*DX3
         V44=      DX2
C
         U(1,NX)=    T11
         U(2,NX)=            T21
         U(3,NX)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,NX)=V41*T11+V42*T21+V43*T31+V44*T41
      ENDDO
      IERR=0
      RETURN
C
 9001 WRITE(6,*) 'XX CSPL1D: NXMAX.GT.NMAX:',NXMAX,NMAX
      IERR=1
      RETURN
 9003 WRITE(6,*) 'XX CSPL1D: TDMSCD ERROR: IERR=',IERR
      IERR=3
      RETURN
      END
C
C     ****** One-Dimensional Spline Interpolation ******
C       **** Calculation of interpolation  ****
C
      SUBROUTINE CSPL1DF(X0,F0,X,U,NXMAX,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION X(NXMAX)
      COMPLEX*16 U(4,NXMAX),F0
C
      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
C
 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2
C
      DX=X0-X(NX-1)
C
      F0= U(1,NX)
     &  + U(2,NX)*DX
     &  + U(3,NX)*DX*DX
     &  + U(4,NX)*DX*DX*DX
      IERR=0
      RETURN
      END
