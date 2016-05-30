!
!     ---- Calculate components of dielectric tensor
!     ---- for energetic particles
!     ----                             by futakuchi
 
SUBROUTINE WMDPFA(CX,CFN,COEF,RHOR,CPM,CQM,CRM,MODEFA)
  USE libdsp
  USE bpsd_constants,ONLY : CI,PI
  IMPLICIT NONE
  COMPLEX(8),INTENT(IN):: CX,CFN,COEF
  REAL(8),INTENT(IN):: RHOR
  COMPLEX(8),INTENT(OUT):: CPM,CQM,CRM
  INTEGER,INTENT(IN):: MODEFA

  COMPLEX(8):: SUMf,Zf,DZf,CEX
  COMPLEX(8):: Pdis
  REAL(8):: vpf, vnf
  REAL(8):: dvf, dhf
  INTEGER::  iv

  IF(ABS(CX).GT.5.D0) THEN
     CEX=(0.D0,0.D0)
  ELSE
     CEX=CX*EXP(-(CX*CX))
  ENDIF

  SUMf=(0.D0,0.D0)

  IF(MODEFA.EQ.1) THEN
!---use libdsp-plasma dispersion function Z
     CALL DSPFN(CX,Zf,DZf)
     SUMf=SQRT(PI)*Zf
     CPM=-CFN*COEF*RHOR*RHOR*(CX*(5.D-1+CX*CX+CX*CX*CX*CX)*SUMf &
                            +CX*CX*(1.5D0+CX*CX)*SQRT(PI))*4.D0*CI/PI
     CQM=-CFN*COEF*RHOR*     (CX*CX*(5.D-1+CX*CX)*SUMf &
                            +CX*(1.D0+CX*CX)*SQRT(PI))*4.D0*CI/PI
     CRM=-CFN*COEF*         ((CX*CX*CX)*SUMf &
                            +(CX*CX)*SQRT(PI))*4.D0*CI/PI
!---           
  ELSE
!--suuchisekibunn
     dhf=1.D-8
     dvf=1.D-2
     vpf=REAL(CX)+dhf
     vnf=REAL(CX)-dhf

     DO iv=1,200
        SUMf=SUMf &
            +dvf*(Pdis(vpf,CX) &
                 +Pdis(vpf+dvf,CX) &
                 +Pdis(vnf,CX) &
                 +Pdis(vnf-dvf,CX))/2.D0
        vpf=vpf+dvf
        vnf=vnf-dvf
     END DO
     dvf=1.D-1
     DO iv=1,200
        SUMf=SUMf &
           +dvf*(Pdis(vpf,CX) &
                +Pdis(vpf+dvf,CX) &
                +Pdis(vnf,CX) &
                +Pdis(vnf-dvf,CX))/2.D0
        vpf=vpf+dvf
        vnf=vnf-dvf
     END DO
     dvf=1.D0
     DO iv=1,200
        SUMf=SUMf &
           +dvf*(Pdis(vpf,CX) &
                +Pdis(vpf+dvf,CX) &
                +Pdis(vnf,CX) &
                +Pdis(vnf-dvf,CX))/2.D0
        vpf=vpf+dvf
        vnf=vnf-dvf
     END DO

     CPM=CFN*COEF*RHOR*RHOR*(CX*(5.D-1+CX*CX+CX*CX*CX*CX)*SUMf &
                            -CX*CX*(1.5D0+CX*CX)*SQRT(PI))*5.D-1*CI/PI
     CQM=CFN*COEF*RHOR*     (CX*CX*(5.D-1+CX*CX)*SUMf &
                            -CX*(1.D0+CX*CX)*SQRT(PI))*CI/PI
     CRM=CFN*COEF*          (CX*CX*CX*SUMf &
                            -(CX*CX)*SQRT(PI))*2.D0*CI/PI
 
    IF(AIMAG(CX).EQ.0.D0) THEN
        CPM=CPM &
           +CFN*COEF*RHOR*RHOR*CEX*CX*(5.D-1+CX*CX+CX*CX*CX*CX)*5.D-1/CX
        CQM=CQM &
           +CFN*COEF*RHOR     *CEX*CX*CX*(5.D-1+CX*CX)/CX
        CRM=CRM &
           +CFN*COEF          *CEX*(CX*CX*CX)*2.D0/CX
     ELSEIF(AIMAG(CX).LT.0.D0) THEN
        CPM=CPM &
           +2.D0*CFN*COEF*RHOR*RHOR*CEX*CX*(5.D-1+CX*CX+CX*CX*CX*CX)*5.D-1/CX
        CQM=CQM &
           +2.D0*CFN*COEF*RHOR     *CEX*CX*CX*(5.D-1+CX*CX)/CX
        CRM=CRM &
           +2.D0*CFN*COEF          *CEX*(CX*CX*CX)*2.D0/CX
     ELSE
        CPM=CPM 
        CQM=CQM
        CRM=CRM
     ENDIF
!--
  ENDIF

  RETURN
END SUBROUTINE WMDPFA
    
FUNCTION Pdis(v,CX)
  IMPLICIT NONE
  REAL(8),INTENT(IN):: v
  COMPLEX(8),INTENT(IN):: CX
  COMPLEX(8):: Pdis
  REAL(8):: EXPVV

  IF(ABS(v).LT.5.D0) THEN
     EXPVV=EXP(-v*v)
  ELSE
     EXPVV=0.D0
  END IF

  Pdis=EXPVV/(CX-v)

END FUNCTION Pdis
