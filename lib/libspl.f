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
C                          3 : DERIVATIVES FX(1) AXD FX(NXMAX) ARE GIVEN
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
      DIMENSION AX(4,NMAX),IPX(NMAX),BX(NMAX)
      DIMENSION WK(NMAX)
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
      EPS=1.D-14
      CALL BGLU1(AX,NXMAX,1,1,4,EPS,WK,IPX,IERR)
      IF(IERR.NE.0) GOTO 9003
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
      CALL BGSLV1(AX,NXMAX,1,1,4,BX,IPX)
      DO NX=1,NXMAX
         FX(NX)=BX(NX)
      ENDDO
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
 9001 IERR=1
      WRITE(6,*) 'XX SPL1D: NXMAX.GT.NMAX:',NXMAX,NMAX
      RETURN
 9003 IERR=3
      WRITE(6,*) 'XX SPL1D: BGLU1 ERROR'
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
      DIMENSION AX(4,NMAX),IPX(NMAX),BX(NMAX)
      DIMENSION AY(4,NMAX),IPY(NMAX),BY(NMAX)
      DIMENSION WK(NMAX)
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
C
      EPS=1.D-14
      CALL BGLU1(AX,NXMAX,1,1,4,EPS,WK,IPX,IERR)
      IF(IERR.NE.0) GOTO 9003
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
C
      EPS=1.D-14
      CALL BGLU1(AY,NYMAX,1,1,4,EPS,WK,IPY,IERR)
      IF(IERR.NE.0) GOTO 9004
C
      DO NY=1,NYMAX
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
         CALL BGSLV1(AX,NXMAX,1,1,4,BX,IPX)
         DO NX=1,NXMAX
            FX(NX,NY)=BX(NX)
         ENDDO
      ENDDO
C
      DO NX=1,NXMAX
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
         CALL BGSLV1(AY,NYMAX,1,1,4,BY,IPY)
         DO NY=1,NYMAX
            FY(NX,NY)=BY(NY)
         ENDDO
      ENDDO
C
      DO NY=1,NYMAX,NYMAX-1
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
         CALL BGSLV1(AX,NXMAX,1,1,4,BX,IPX)
         DO NX=1,NXMAX
            FXY(NX,NY)=BX(NX)
         ENDDO
      ENDDO
C
      DO NX=1,NXMAX
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
         CALL BGSLV1(AY,NYMAX,1,1,4,BY,IPY)
         DO NY=1,NYMAX
            FXY(NX,NY)=BY(NY)
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
 9001 IERR=1
      WRITE(6,*) 'XX SPL2D: NXMAX.GT.NMAX:',NXMAX,NMAX
      RETURN
 9002 IERR=2
      WRITE(6,*) 'XX SPL2D: NYMAX.GT.NMAX:',NYMAX,NMAX
      RETURN
 9003 IERR=3
      WRITE(6,*) 'XX SPL2D: BGLU1 ERROR'
      RETURN
 9004 IERR=4
      WRITE(6,*) 'XX SPL2D: BGLU1 ERROR'
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
      DIMENSION X(NXMAX),IPX(NMAX),AX(4,NMAX),WK(NMAX)
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
      EPS=1.D-14
      CALL BGLU1(AX,NXMAX,1,1,4,EPS,WK,IPX,IERR)
      IF(IERR.NE.0) GOTO 9003
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
      CALL CBGSLV1(AX,NXMAX,1,1,4,BX,IPX)
      DO NX=1,NXMAX
         FX(NX)=BX(NX)
      ENDDO
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
 9001 IERR=1
      WRITE(6,*) 'XX CSPL1D: NXMAX.GT.NMAX:',NXMAX,NMAX
      RETURN
 9003 IERR=3
      WRITE(6,*) 'XX CSPL1D: BGLU1 ERROR'
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
C
C     *** LU DECOMPOSITION MATRIX SOLVER ***
C
      SUBROUTINE BGLU1( A, N, ML, MU, MM, EPS, WK, IP, IER )
C
C        BGLU1
C               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
C
C               SOLVES SIMULTANEOUS LINEAR EQUATIONS
C               BY GAUSSIAN ELIMINATION METHOD FOR GENERAL BAND MATRIX.
C
C        INPUT - -
C             A(-ML:MU+ML,N)
C                      R *8  : 2-DIM. ARRAY CONTAINING REAL BAND MATRIX.
C             N        I *4  : ORDER OF MATRIX.
C             ML       I *4  : LOWER BAND WIDTH.
C             MU       I *4  : UPPER BAND WIDTH.
C             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE
C                              MATRIX. ( STANDARD VALUE 3.52D-15 )
C        OUTPUT - -
C             A(-ML:MU+ML,N)
C                            : RESULT OF GAUSSIAN ELIMINATION.
C             IP(N)    I *4  : PIVOT NUMBER.
C             IER      I *4  : = 0,  FOR NORMAL EXECUTION.
C                              = 1,  FOR SINGULAR MATRIX.
C                              = 3,  FOR INVALID ARGUEMENT.
C        WORKING  -
C             WK(N)    R *8  : 1-DIM. ARRAY.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION A(-ML:MU+ML,*), IP(*), WK(*)
      DIMENSION A(-ML:MM-ML-1,*), IP(*), WK(*)
C             LEFT HAND SIDE
      IF( EPS.LT.0.0D0 )  EPS = 3.52D-15
      IF( ( N.LE.0 ).OR.( ML.LE.0 ).OR.( MU.LE.0 ).OR.
     &    ( ML.GE.N ).OR.( MU.GE.N ) )  THEN
         IER = 3
         WRITE(*,*) '  (SUBR. BGLU1)  INVALID ARGUMENT.  ML, MU, N =',
     &              ML, MU, N
         RETURN
      END IF
C
      DO 10 K=1,N
      DO 10 I=1,ML
         A(MU+I,K)=0.D0
   10 CONTINUE
      IER = 0
      DO 100 K = 1, N
C             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN.
        AMAX = ABS(A(0,K))
        IPK = K
        DO 110 I = K+1, MIN(K+ML,N)
          AIK = ABS(A(K-I,I))
          IF( AIK.GT.AMAX )  THEN
             IPK = I
             AMAX = AIK
          END IF
  110   CONTINUE
        IP(K) = IPK
C
        IF( AMAX.GT.EPS )  THEN
           IF( IPK.NE.K )  THEN
              DO 120 J = K, MIN(K+MU+ML,N)
                W = A(J-IPK,IPK)
                A(J-IPK,IPK) = A(J-K,K)
                A(J-K,K) = W
  120         CONTINUE
           END IF
C
           DO 125 J = K+1, MIN(K+MU+ML,N)
  125        WK(J) = A(J-K,K)
C             COMPUTE ALFA AND PERFORM GAUSSIAN ELIMINATION.
           DO 130 I = K+1, MIN(K+ML,N)
             A(K-I,I) = -A(K-I,I)/A(0,K)
             T = A(K-I,I)
             DO 140 J = K+1, MIN(K+MU+ML,N)
  140          A(J-I,I) = A(J-I,I)+T*WK(J)
  130      CONTINUE
C             MATRIX IS SINGULAR.
        ELSE
           IER = 1
           IP(K) = K
           DO 150 I = K+1, MIN(K+ML,N)
  150        A(K-I,I) = 0.0D0
           WRITE(*,*)  '  (SUBR. BGLU1)  MATRIX IS SINGULAR AT K =', K
           RETURN
        END IF
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE BGSLV1( A, N, ML, MU, MM, B, IP )
C
C        BGSLV1
C               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
C
C               SOLVES SIMULTANEOUS LINEAR EQUATIONS
C               BY GAUSSIAN ELIMINATION METHOD FOR GENERAL BAND MATRIX.
C
C        INPUT - -
C             A(-ML:MU+ML,N)
C                      R *8  : RESULT OF GAUSSIAN ELIMINATION.
C             N        I *4  : ORDER OF MATRIX.
C             ML       I *4  : LOWER BAND WIDTH.
C             MU       I *4  : UPPER BAND WIDTH.
C             B(N)     R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND
C                              SIDE VECTOR.
C             IP(N)    I *4  : PIVOT NUMBER.
C        OUTPUT - -
C             B(N)           : SOLUTION.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION A(-ML:MU+ML,*), B(*), IP(*)
      DIMENSION A(-ML:MM-ML-1,*), B(*), IP(*)
C             FORWARD ELIMINATION PROCESS
      DO 100 K = 1, N
        IPK = IP(K)
        IF( IPK.NE.K )  THEN
           W =B(IPK)
           B(IPK) = B(K)
           B(K) = W
        END IF
C             GAUSSIAN ELIMINATION
        T = B(K)
        DO 110 I = K+1, MIN(K+ML,N)
  110     B(I) = B(I)+A(K-I,I)*T
  100 CONTINUE
C             BACKWARD SUBSTITUTION PROCESS
      DO 200 K = N, 1, -1
        S = -B(K)
        DO 210 J = K+1, MIN(K+MU+ML,N)
  210     S = S+A(J-K,K)*B(J)
        B(K) = -S/A(0,K)
  200 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CBGSLV1( A, N, ML, MU, MM, B, IP )
C
C        BGSLV1
C               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
C
C               SOLVES SIMULTANEOUS LINEAR EQUATIONS
C               BY GAUSSIAN ELIMINATION METHOD FOR GENERAL BAND MATRIX.
C
C        INPUT - -
C             A(-ML:MU+ML,N)
C                      R *8  : RESULT OF GAUSSIAN ELIMINATION.
C             N        I *4  : ORDER OF MATRIX.
C             ML       I *4  : LOWER BAND WIDTH.
C             MU       I *4  : UPPER BAND WIDTH.
C             B(N)     R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND
C                              SIDE VECTOR.
C             IP(N)    I *4  : PIVOT NUMBER.
C        OUTPUT - -
C             B(N)           : SOLUTION.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION A(-ML:MU+ML,*), B(*), IP(*)
      DIMENSION IP(*),A(-ML:MM-ML-1,*)
      COMPLEX*16 B(*),W,T,S
C             FORWARD ELIMINATION PROCESS
      DO 100 K = 1, N
        IPK = IP(K)
        IF( IPK.NE.K )  THEN
           W =B(IPK)
           B(IPK) = B(K)
           B(K) = W
        END IF
C             GAUSSIAN ELIMINATION
        T = B(K)
        DO 110 I = K+1, MIN(K+ML,N)
  110     B(I) = B(I)+A(K-I,I)*T
  100 CONTINUE
C             BACKWARD SUBSTITUTION PROCESS
      DO 200 K = N, 1, -1
        S = -B(K)
        DO 210 J = K+1, MIN(K+MU+ML,N)
  210     S = S+A(J-K,K)*B(J)
        B(K) = -S/A(0,K)
  200 CONTINUE
      RETURN
      END
