C     $Id$
C
C
C     ***********************************************************
C
C           CEILING FUNCTION FOR LOG10 PLOT
C
C     ***********************************************************
C
      DOUBLE PRECISION FUNCTION PLOG(X,XMIN,XMAX)
C
      REAL*8 X,XMIN,XMAX
C
      IF(X.LT.XMIN) THEN
         PLOG=LOG10(XMIN)
      ELSEIF(X.GT.XMAX) THEN
         PLOG=LOG10(XMAX)
      ELSE
         PLOG=LOG10(X)
      ENDIF
C
      RETURN
      END
C
C     ****** CONTOUR PLOT : XY, VARIABLE POSITION, PATTERN ******
C
      SUBROUTINE CONTP5(Z,X,Y,NXA,NXMAX,NYMAX,
     &                  ZORG,ZSTEP,NSTEP,IPRD,IPAT,KA)
C
      IMPLICIT LOGICAL(L)
      COMMON /GSEQ00/ NXAC,NXMAXC,NYMAXC
C
      EXTERNAL CONTS5
      DIMENSION Z(NXA,*),KA(2,*),X(NXA,*),Y(NXA,*)
      DATA EPS/1.E-32/
C
      IF(ABS(ZSTEP).LT.EPS) RETURN
C
      NXAC=NXA
      NXMAXC=NXMAX
      NYMAXC=NYMAX
      CALL CONTP0(Z,X,Y,NXA,NXMAX,NYMAX,
     &            ZORG,ZSTEP,NSTEP,IPRD,IPAT,KA,CONTS5)
      RETURN
      END
C
C     ****** CONTOUR PLOT SLAVE ROUTINE : XY VARIABLE ******
C
      SUBROUTINE CONTS5(NAX,NAY,NBX,NBY,U0,UA,UB,X,Y,IPAT,IND)
C
      COMMON /GSEQ00/ NXAC,NXMAXC,NYMAXC
      DIMENSION X(*),Y(*)
C
      CALL INQGDEFIN(PXMIN,PXMAX,PYMIN,PYMAX,
     &                      GXMIN,GXMAX,GYMIN,GYMAX)
      DX=(PXMAX-PXMIN)/(GXMAX-GXMIN)
      DY=(PYMAX-PYMIN)/(GYMAX-GYMIN)
C
      NAXL=ABS(NAX)
      IF(NAXL.GT.NXMAXC) NAXL=1
      NAYL=ABS(NAY)
      IF(NAYL.GT.NYMAXC) NAYL=1
      NA=NXAC*(NAYL-1)+NAXL
C
      NBXL=ABS(NBX)
      IF(NBXL.GT.NXMAXC) NBXL=1
      NBYL=ABS(NBY)
      IF(NBYL.GT.NYMAXC) NBYL=1
      NB=NXAC*(NBYL-1)+NBXL
C
      IF(NAX.GT.0) THEN
         XA=DX*(X(NA)-GXMIN)+PXMIN
         YA=DY*(Y(NA)-GYMIN)+PYMIN
      ELSE
         XA=DX*(0.5*(X(NA)+X(NA+1))-GXMIN)+PXMIN
         YA=DY*(0.5*(Y(NA)+Y(NA+NXAC))-GYMIN)+PYMIN
      ENDIF
      IF(NBX.GT.0) THEN
         XB=DX*(X(NB)-GXMIN)+PXMIN
         YB=DY*(Y(NB)-GYMIN)+PYMIN
      ELSE
         XB=DX*(0.5*(X(NB)+X(NB+1   ))-GXMIN)+PXMIN
         YB=DY*(0.5*(Y(NB)+Y(NB+NXAC))-GYMIN)+PYMIN
      ENDIF
      IF(UB.NE.UA) THEN
         XS=(XB-XA)*(U0-UA)/(UB-UA)+XA
         YS=(YB-YA)*(U0-UA)/(UB-UA)+YA
         IF(IND.EQ.1) THEN
            CALL MOVEPT(XS,YS,IPAT)
         ELSEIF(IND.EQ.-1) THEN
            CALL MOVEPT(XS,YS,-IPAT)
         ELSEIF(IND.EQ.0) THEN
            CALL DRAWPT(XS,YS)
         ENDIF
      ENDIF
C      WRITE(6,*) NAX,NAY,NBX,NBY,XS,YS
      RETURN
      END
