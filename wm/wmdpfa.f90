!
!     ---- Calculate components of dielectric tensor
!     ---- for energetic particles
!     ----                             by futakuchi
 
SUBROUTINE WMDPFA(CX,CFN,COEF,RHOR,CPM,CQM,CRM,MODEFA)

  USE bpsd_constants,ONLY : CI,PI
  IMPLICIT NONE
  COMPLEX(8),INTENT(IN):: CX,CFN,COEF
  REAL(8),INTENT(IN):: RHOR
  COMPLEX(8),INTENT(OUT):: CPM,CQM,CRM
  INTEGER,INTENT(IN):: MODEFA

  COMPLEX(8):: SUM
  COMPLEX(8):: Fun
  REAL(8):: vp, vn
  REAL(8):: dv, dh
  INTEGER::  i

  SUM=(0.D0,0.D0)
  dh=1.D-8
  dv=1.D-2
  vp=REAL(CX)+dh
  vn=REAL(CX)-dh

  DO i=1,200

     SUM=SUM &
        +dv*(Fun(vp,CX)+Fun(vp+dv,CX)+Fun(vn,CX)+Fun(vn-dv,CX))/2.D0

     vp=vp+dv
     vn=vn-dv
  END DO

  dv=5.D-1

  DO i=1,200
     SUM=SUM &
          +dv*(Fun(vp,CX)+Fun(vp+dv,CX)+Fun(vn,CX)+Fun(vn-dv,CX))/2.D0

     vp=vp+dv
     vn=vn-dv
  END DO

  CPM=CFN*COEF*RHOR*RHOR*(CX*(1.D0/2.D0+CX*CX+CX*CX*CX*CX)*SUM &
                         -CX*CX*(3.D0+2.D0*CX)*(PI**0.5D0))*CI
  CQM=CFN*COEF*RHOR*     (CX*CX*(1.D0/2.D0+CX*CX)*SUM &
                         -CX*(2.D0+2.D0*CX*CX)*(PI**0.5D0))*CI
  CRM=CFN*COEF*         ((CX*CX*CX)*SUM &
                         -2.D0*(CX*CX)*(PI**0.5D0))*CI
  RETURN
END SUBROUTINE WMDPFA
    
FUNCTION Fun(v,CX)
  IMPLICIT NONE
  REAL(8),INTENT(IN):: v
  COMPLEX(8),INTENT(IN):: CX
  COMPLEX(8):: Fun
  REAL(8):: EXPV

  IF(ABS(v).LT.5.D0) THEN
     EXPV=EXP(-v*v)
  ELSE
     EXPV=0.D0
  END IF

  Fun=EXPV/(CX-v)

END FUNCTION Fun
