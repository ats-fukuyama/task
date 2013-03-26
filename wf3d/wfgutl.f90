!     $Id$
!     *********************

!     CLIPPING FOR GRAPHICS

!     *********************

FUNCTION GCLIP(D)

  implicit none
  real(4) :: GCLIP
  real(8) :: D

  IF(ABS(D).GT.1.D-15) THEN
     GCLIP=REAL(D)
  ELSE
     GCLIP=0.0
  ENDIF
  RETURN
END FUNCTION GCLIP

!     *****************************

!     OPTIMUM NUM LENGTH FOR GVALUE

!     *****************************

FUNCTION NGVLEN(GSTEP)

  implicit none
  integer :: NGVLEN,NGX
  real(4) :: GSTEP

  NGX = -INT(LOG10(DBLE(GSTEP*0.11)))
  IF(NGX.LT.-5)THEN
     NGX=-1
  ELSE
     IF(NGX.LT.0) NGX=0
     IF(NGX.GT.5) NGX=-1
  ENDIF
  NGVLEN=NGX
  RETURN
END FUNCTION NGVLEN

!     *****************************

!     INTERPORATE RGB

!     *****************************

!      SUBROUTINE GUSRGB(GL,GRGBL,NRGB,GLA,GRGBLA)

!      DIMENSION GRGBL(3),GLA(NRGB),GRGBLA(3,NRGB)

!      DO NDO=2,NRGB
!         N=NDO
!         IF(GLA(N).GT.GL) GOTO 9
!      ENDDO
!    9 CONTINUE

!      GFACT=(GL-GLA(N-1))/(GLA(N)-GLA(N-1))
!      GRGBL(1)=GRGBLA(1,N-1)*(1.0-GFACT)+GRGBLA(1,N)*GFACT
!      GRGBL(2)=GRGBLA(2,N-1)*(1.0-GFACT)+GRGBLA(2,N)*GFACT
!      GRGBL(3)=GRGBLA(3,N-1)*(1.0-GFACT)+GRGBLA(3,N)*GFACT

!      RETURN
!    END SUBROUTINE GUSRGB
      
