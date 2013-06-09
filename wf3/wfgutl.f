C     $Id$
C     *********************
C
C     CLIPPING FOR GRAPHICS
C
C     *********************
C
      FUNCTION GCLIP(D)
      REAL*8 D
      IF(ABS(D).GT.1.D-15) THEN
         GCLIP=REAL(D)
      ELSE
         GCLIP=0.0
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
C     *****************************
C
C     INTERPORATE RGB
C
C     *****************************
C
      SUBROUTINE GUSRGB(GL,GRGBL,NRGB,GLA,GRGBLA)
C
      DIMENSION GRGBL(3),GLA(NRGB),GRGBLA(3,NRGB)
C
      DO NDO=2,NRGB
         N=NDO
         IF(GLA(N).GT.GL) GOTO 9
      ENDDO
    9 CONTINUE
C
      GFACT=(GL-GLA(N-1))/(GLA(N)-GLA(N-1))
      GRGBL(1)=GRGBLA(1,N-1)*(1.0-GFACT)+GRGBLA(1,N)*GFACT
      GRGBL(2)=GRGBLA(2,N-1)*(1.0-GFACT)+GRGBLA(2,N)*GFACT
      GRGBL(3)=GRGBLA(3,N-1)*(1.0-GFACT)+GRGBLA(3,N)*GFACT
C
      RETURN
      END
