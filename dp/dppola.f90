! dppola.f90

MODULE dppola

  PRIVATE
  PUBLIC dp_pola
  
CONTAINS

  ! --- calculate polization vector ---

  SUBROUTINE dp_pola(crf,ckx,cky,ckz,xpos,ypos,zpos,cdet,cepola,err)

    USE dpcomm
    USE dpdisp
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: crf,ckx,cky,ckz
    REAL(rkind),INTENT(IN):: xpos,ypos,zpos
    COMPLEX(rkind),INTENT(OUT):: cdet(3,3),cepola(3)
    REAL(rkind),INTENT(OUT):: err

    COMPLEX(rkind):: cdet11,cdet12,cdet13,cdet21,cdet22,cdet23, &
         cdet31,cdet32,cdet33
    COMPLEX(rkind):: cd1,cd2,cd3,ce1,ce2,ce3
    REAL(rkind):: aca

    CALL DP_DISP(crf,ckx,cky,ckz,xpos,ypos,zpos,cdet)

    cdet11=cdet(1,1)
    cdet12=cdet(1,2)
    cdet13=cdet(1,3)
    cdet21=cdet(2,1)
    cdet22=cdet(2,2)
    cdet23=cdet(2,3)
    cdet31=cdet(3,1)
    cdet32=cdet(3,2)
    cdet33=cdet(3,3)

    cd1=cdet12*cdet23-cdet13*cdet22
    cd2=cdet13*cdet21-cdet11*cdet23
    cd3=cdet11*cdet22-cdet12*cdet21
    IF(ABS(cd1).GE.ABS(cd2)) THEN
       IF(ABS(cd1).GE.ABS(cd3)) THEN
          ce1=(1.D0,0.D0)
          ce2=-( cdet23*cdet11-cdet22*cdet21)/cd1
          ce3=-(-cdet13*cdet11+cdet12*cdet21)/cd1
       ELSE
          ce1=-( cdet22*cdet13-cdet12*cdet23)/cd3
          ce2=-(-cdet21*cdet13+cdet11*cdet23)/cd3
          ce3=(1.D0,0.D0)
       END IF
    ELSE
       IF(ABS(cd2).GE.ABS(cd3)) THEN
          ce1= ( cdet23*cdet12-cdet13*cdet22)/cd2
          ce2=(1.D0,0.D0)
          ce3= (-cdet21*cdet12+cdet11*cdet22)/cd2
       ELSE
          ce1=-( cdet22*cdet13-cdet12*cdet23)/cd3
          ce2=-(-cdet21*cdet13+cdet11*cdet23)/cd3
          ce3=(1.D0,0.D0)
       END IF
    ENDIF
    aca=SQRT(ABS(ce1)**2+ABS(ce2)**2+ABS(ce3)**2)
    cepola(1)=ce1/aca
    cepola(2)=ce2/aca
    cepola(3)=ce3/aca
    err=ABS(cdet31*ce1+cdet32*ce2+cdet33*ce3)
    RETURN
  END SUBROUTINE dp_pola
END MODULE dppola
    
