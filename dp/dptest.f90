! dptest.f90

MODULE dptest

PRIVATE

PUBLIC dp_test

CONTAINS

  !     ***** TASK/DP TEST *****

  SUBROUTINE dp_test
    IMPLICIT NONE
    INTEGER:: id

1   CONTINUE
    WRITE(6,'(A)') '## dp_test: input 1-3 0/9 for end?'
    READ(5,'(I4)',ERR=1,END=9) id
    SELECT CASE(id)
    CASE(1)
       CALL dp_test1
    CASE(2)
       CALL dp_test2
    CASE(0,9)
       GO TO 9
    END SELECT
    GOTO 1

9   CONTINUE
    RETURN
  END SUBROUTINE dp_test
       
  SUBROUTINE dp_test1

    USE dpcomm_local
    USE plprof
    USE plprofw
    USE dptnsr0
    IMPLICIT NONE
    INTEGER:: I
    INTEGER:: NSMAX_SAVE,NRMAX_SAVE,MODELV_SAVE(NSM),MODELP_SAVE(NSM)
    REAL(rkind):: RF0_SAVE,RFI0_SAVE,RKZ0_SAVE,RKX0_SAVE,RL,RHON
    COMPLEX(rkind):: CW,CKPR,CKPP
    COMPLEX(rkind):: CD1(6),CD2(6),CD3(6),CD4(6)
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd

    NSMAX_SAVE=NSMAX
    NRMAX_SAVE=NRMAX_DP
    MODELV_SAVE(1:NSMAX)=MODELV(1:NSMAX)
    MODELP_SAVE(1:NSMAX)=MODELP(1:NSMAX)
    RF0_SAVE=RF0
    RFI0_SAVE=RFI0
    RKZ0_SAVE=RKZ0
    RKX0_SAVE=RKX0
    RL=RR
    
1001 CONTINUE
    WRITE(6,*) '# INPUT: RL,RF0,RFI0,RKZ0,RKX0 ='
    WRITE(6,'(10X,1P5E12.4)') RL,RF0,RFI0,RKZ0,RKX0
    READ(5,*,ERR=1001,END=1002) RL,RF0,RFI0,RKZ0,RKX0
    IF(RF0.LE.0.D0) GOTO 1002
    CW=2.D0*PI*DCMPLX(RF0,RFI0)*1.D6
    CKPR=MAX(RKZ0,1.D-8)
    CKPP=RKX0
    CALL PL_MAG(RL,0.D0,0.D0,mag)
    RHON=mag%rhon
    CALL PL_PROFW(RHON,plfw)
    CALL PL_GRAD(RHON,grd)

    MODELV(1)=0
    MODELP(1)=6
    CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD1)
    MODELV(1)=1
    MODELP(1)=6
    CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD2)
    MODELV(1)=0
    MODELP(1)=9
    CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD3)
    MODELV(1)=3
    MODELP(1)=6
    CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD4)

    WRITE(6,602) 
602 FORMAT(6X,'MODELV=0',13X,'MODELV=1', &
         13X,'MODELV=0',13X,'MODELV=3')
    WRITE(6,603) 
603 FORMAT(6X,'MODELP=6',13X,'MODELP=6', &
         13X,'MODELP=9',13X,'MODELP=6')
    WRITE(6,604) (CD1(I),CD2(I),CD3(I),CD4(I),I=1,6)
604 FORMAT((1PE9.2,1P7E10.2))
    GOTO 1001

1002 CONTINUE
    NSMAX=NSMAX_SAVE
    NRMAX_DP=NRMAX_SAVE
    MODELV(1:NSMAX)=MODELV_SAVE(1:NSMAX)
    MODELP(1:NSMAX)=MODELP_SAVE(1:NSMAX)
    RF0=RF0_SAVE
    RFI0=RFI0_SAVE
    RKZ0=RKZ0_SAVE
    RKX0=RKX0_SAVE
    RETURN
  END SUBROUTINE dp_test1

  SUBROUTINE dp_test2

    USE dpcomm_local
    USE plprof
    USE plprofw
    USE dpdisp
    USE dpdispx
    IMPLICIT NONE
    INTEGER:: I
    INTEGER:: NSMAX_SAVE,NRMAX_SAVE,MODELV_SAVE(NSM),MODELP_SAVE(NSM)
    REAL(rkind):: RF0_SAVE,RFI0_SAVE,RKZ0_SAVE,RKX0_SAVE,RL,RHON
    COMPLEX(rkind):: CW,CKPR,CKPP
    COMPLEX(rkind):: CD1,CD2,CD3,CD4,CD5,CDTNS1(3,3),CDTNS2(3,3),CRF
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd

    NSMAX_SAVE=NSMAX
    NRMAX_SAVE=NRMAX_DP
    MODELV_SAVE(1:NSMAX)=MODELV(1:NSMAX)
    MODELP_SAVE(1:NSMAX)=MODELP(1:NSMAX)
    RF0_SAVE=RF0
    RFI0_SAVE=RFI0
    RKZ0_SAVE=RKZ0
    RKX0_SAVE=RKX0
    RL=RR
    
1001 CONTINUE
    WRITE(6,*) '# INPUT: RL,RF0,RFI0,RKZ0,RKX0 ='
    WRITE(6,'(10X,1P5E12.4)') RL,RF0,RFI0,RKZ0,RKX0
    READ(5,*,ERR=1001,END=1002) RL,RF0,RFI0,RKZ0,RKX0
    IF(RF0.LE.0.D0) GOTO 1002

    CRF=DCMPLX(RF0,RFI0)*1.D6
    CW=2.D0*PI*DCMPLX(RF0,RFI0)*1.D6
    CKPR=MAX(RKZ0,1.D-8)
    CKPP=RKX0
    CALL PL_MAG(RL,0.D0,0.D0,mag)
    RHON=mag%rhon
    CALL PL_PROFW(RHON,plfw)
    CALL PL_GRAD(RHON,grd)

    CD1=CFDISP(CRF,CKPP,(0.D0,0.D0),CKPR,RL,0.D0,0.D0)
    CD2=CF_DISPX(CRF,CKPP,(0.D0,0.D0),CKPR,RL,0.D0,0.D0)
    CD3=CF_DISPX(CRF,CKPP,(0.D0,0.D0),CKPR,RL,0.D0,0.D0,1)
    CALL DP_DTNS(CRF,CKPP,(0.D0,0.D0),CKPR,RL,0.D0,0.D0,CDTNS1)
    CALL DP_DTNS(CRF,CKPP,(0.D0,0.D0),CKPR,RL,0.D0,0.D0,CDTNS2,1)

    WRITE(6,'(A,2ES12.4)') '## DISPX:    ',CD1
    WRITE(6,'(A,2ES12.4)') '## CF_DISPX0:',CD2
    WRITE(6,'(A,2ES12.4)') '## CF_DISPX1:',CD3
    WRITE(6,'(A)'        ) '## DP_DTNS0:'
    WRITE(6,'(A,6ES12.4)') '   ',CDTNS1(1,1),CDTNS1(1,2),CDTNS1(1,3)
    WRITE(6,'(A,6ES12.4)') '   ',CDTNS1(2,1),CDTNS1(2,2),CDTNS1(2,3)
    WRITE(6,'(A,6ES12.4)') '   ',CDTNS1(3,1),CDTNS1(3,2),CDTNS1(3,3)
    WRITE(6,'(A)'        ) '## DP_DTNS1:'
    WRITE(6,'(A,6ES12.4)') '   ',CDTNS2(1,1),CDTNS2(1,2),CDTNS2(1,3)
    WRITE(6,'(A,6ES12.4)') '   ',CDTNS2(2,1),CDTNS2(2,2),CDTNS2(2,3)
    WRITE(6,'(A,6ES12.4)') '   ',CDTNS2(3,1),CDTNS2(3,2),CDTNS2(3,3)

    GOTO 1001

1002 CONTINUE
    NSMAX=NSMAX_SAVE
    NRMAX_DP=NRMAX_SAVE
    MODELV(1:NSMAX)=MODELV_SAVE(1:NSMAX)
    MODELP(1:NSMAX)=MODELP_SAVE(1:NSMAX)
    RF0=RF0_SAVE
    RFI0=RFI0_SAVE
    RKZ0=RKZ0_SAVE
    RKX0=RKX0_SAVE
    RETURN
  END SUBROUTINE dp_test2
END MODULE dptest
