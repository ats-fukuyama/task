C     $Id$
C
C     *********************
C
C     CLIPPING FOR GRAPHICS
C
C     *********************
C
      FUNCTION GCLIP(D)
      REAL*8 D
      IF(ABS(D).LE.1.D-30) THEN
         GCLIP=0.0
      ELSE IF(D.GT. 1.D30) THEN
         GCLIP= 1.E30
      ELSE IF(D.LT.-1.D30) THEN
         GCLIP=-1.E30
      ELSE
         GCLIP=SNGL(D)
      ENDIF
      RETURN
      END
C
C     *****************************
C
C     OPTIMUM NUM LENGTH FOR GVALUE
C
C     *****************************
C
      FUNCTION NGVLEN(GSTEP)
C
      NGX = -INT(LOG10(DBLE(GSTEP*0.11)))
      IF(NGX.LT.-5)THEN
         NGX=-1
      ELSE
         IF(NGX.LT.0) NGX=0
         IF(NGX.GT.5) NGX=-1
      ENDIF
      NGVLEN=NGX
      RETURN
      END
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
      COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
      COMMON /GSEQ00/ NXAC,NXMAXC,NYMAXC
C
      EXTERNAL CONTS5
      DIMENSION Z(NXA,*),KA(2,*),X(NXA,*),Y(NXA,*)
      DATA EPS/1.E-32/
C
      IF(.NOT.LGF.OR.ABS(ZSTEP).LT.EPS) RETURN
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
      IMPLICIT LOGICAL(L)
      COMMON /GSGFXY/ DX,DY,PXS,PYS,PXE,PYE,GXS,GYS,GXE,GYE,LGF
      COMMON /GSEQ00/ NXAC,NXMAXC,NYMAXC
      DIMENSION X(*),Y(*)
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
         XA=DX*(X(NA)-GXS)+PXS
         YA=DY*(Y(NA)-GYS)+PYS
      ELSE
         XA=DX*(0.5*(X(NA)+X(NA+1))-GXS)+PXS
         YA=DY*(0.5*(Y(NA)+Y(NA+NXAC))-GYS)+PYS
      ENDIF
      IF(NBX.GT.0) THEN
         XB=DX*(X(NB)-GXS)+PXS
         YB=DY*(Y(NB)-GYS)+PYS
      ELSE
         XB=DX*(0.5*(X(NB)+X(NB+1))-GXS)+PXS
         YB=DY*(0.5*(Y(NB)+Y(NB+NXAC))-GYS)+PYS
      ENDIF
      XS=(XB-XA)*(U0-UA)/(UB-UA)+XA
      YS=(YB-YA)*(U0-UA)/(UB-UA)+YA
      IF(IND.EQ.1) THEN
         CALL MOVEPT(XS,YS,IPAT)
      ELSEIF(IND.EQ.-1) THEN
         CALL MOVEPT(XS,YS,-IPAT)
      ELSEIF(IND.EQ.0) THEN
         CALL DRAWPT(XS,YS)
      ENDIF
C      WRITE(6,*) NAX,NAY,NBX,NBY,XS,YS
      RETURN
      END
