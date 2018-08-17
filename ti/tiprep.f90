! MODULE tiprep

MODULE tiprep

CONTAINS

  SUBROUTINE ti_prep(IERR)

    USE ticomm
    USE plprof
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS,NSA,NZ,NEQ,NR,NV
    REAL(rkind):: RHON
    TYPE(pl_plf_type),DIMENSION(:),ALLOCATABLE:: PLF

    IERR=0

!   *** count NSAMAX: number of active particles ***

    NSA=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(0)
          CONTINUE
       CASE(-1,1,2,5,6)
          NSA=NSA+1
       CASE(10,11,12)
          DO NZ=NZMIN_NS(NS),NZMAX_NS(NS)
             NSA=NSA+1
          END DO
       CASE DEFAULT
          WRITE(6,'(A,2I5)') &
               'XX ti_prep: Undefined ID_NS(NS): ID_NS,NS:', &
               ID_NS(NS),NS
          STOP
       END SELECT
    END DO
    NSAMAX=NSA

!   *** ALLOCATE array for NSAMAX and NRMAX ***

    CALL allocate_ticomm(IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX tiprep: allocate_ticomm ERROR: IERR=',IERR
       STOP
    END IF

!   *** Define NSA variables ***

    NSA=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(0)
          CONTINUE
       CASE(-1,1,2,5,6)
          NSA=NSA+1
          PMA(NSA)=PM(NS)
          PZA(NSA)=PZ(NS)   ! for ID_NS=5,6, PZ will be ealuated later
          NS_NSA(NSA)=NS
          NSA_DN(NSA)=0
          NSA_UP(NSA)=0
       CASE(10,11,12)
          DO NZ=NZMIN_NS(NS),NZMAX_NS(NS)
             NSA=NSA+1
             PMA(NSA)=PM(NS)
             PZA(NSA)=DBLE(NZ)
             NS_NSA(NSA)=NS
             IF(NZ.EQ.NZMIN_NS(NS)) THEN
                NSA_DN(NSA)=0
             ELSE
                NSA_DN(NSA)=NSA-1
             END IF
             IF(NZ.EQ.NZMAX_NS(NS)) THEN
                NSA_UP(NSA)=0
             ELSE
                NSA_UP(NSA)=NSA+1
             END IF
          END DO
       END SELECT
    END DO
    IF(NSA.NE.NSAMAX) THEN
       WRITE(6,*) 'XX ti_prep: INCONSISTENT NSAMAX'
       STOP
    END IF

!   *** Display NSA variables ***

    WRITE(6,*) 'NSA  ','NS   ','ID   ','PMA         ','PZA'
    DO NSA=1,NSAMAX
       WRITE(6,'(3I5,1P2E12.4)') &
            NSA,NS_NSA(NSA),ID_NS(NS_NSA(NSA)),PMA(NSA),PZA(NSA)
    END DO

!   *** Count NEQMAX: Size of equation ***

    NEQ=0
    IF(MODEL_EQB.EQ.1) NEQ=NEQ+1
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(0)
          CONTINUE
       CASE(-1,1,2,5,6)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+1
       CASE(10)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+1
       CASE(11)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+1
       CASE(12)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
       END SELECT
    END DO
    NEQMAX=NEQ
    WRITE(6,'(A,I5)') 'NEQMAX=',NEQMAX

!   *** ALLOCATE array for NEQMAX and NRMAX ***

    CALL allocate_neqmax(IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX tiprep: allocate_neqmax ERROR: IERR=',IERR
       STOP
    END IF

!   *** Define NEQ variables ***

    NEQ=0
    IF(MODEL_EQB.EQ.1) THEN
       NEQ=NEQ+1
       NSA_NEQ(NEQ)=0
       NV_NEQ(NEQ)=1
    END IF
    DO NSA=1,NSAMAX
       NS=NS_NSA(NSA)
       IF(MODEL_EQN.EQ.1) THEN
          SELECT CASE(ID_NS(NS))
          CASE(-1,1,2,5,6,10,11,12)
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=1
          END SELECT
       END IF
       IF(MODEL_EQT.EQ.1) THEN
          SELECT CASE(ID_NS(NS))
          CASE(-1,1,2,5,6,10,11)
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=2
          CASE(12)
             IF(NINT(PZ(NSA)).EQ.NZMIN_NS(NS)) THEN
                NEQ=NEQ+1
                NSA_NEQ(NEQ)=NSA
                NV_NEQ(NEQ)=2
             END IF
          END SELECT
       END IF
       IF(MODEL_EQU.EQ.1) THEN
          SELECT CASE(ID_NS(NS))
          CASE(-1,1,2,5,6,10)
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=3
          CASE(11,12)
             IF(NINT(PZ(NSA)).EQ.NZMIN_NS(NS)) THEN
                NEQ=NEQ+1
                NSA_NEQ(NEQ)=NSA
                NV_NEQ(NEQ)=3
             END IF
          END SELECT
       END IF
    END DO
    IF(NEQ.NE.NEQMAX) THEN
       WRITE(6,*) 'XX ti_prep: INCONSISTENT NEQMAX'
       STOP
    END IF

    NEQ_NVNSA(1:3,0:NSAMAX)=0.D0
    DO NEQ=1,NEQMAX
       NV=NV_NEQ(NEQ)
       NSA=NSA_NEQ(NEQ)
       NEQ_NVNSA(NV,NSA)=NEQ
    END DO

!   *** Display NEQ variables ***

    WRITE(6,*) 'NEQ  ','NS   ','NSA  ','NV   '
    DO NEQ=1,NEQMAX
       WRITE(6,'(4I5)') &
            NEQ,NS_NSA(NSA_NEQ(NEQ)),NSA_NEQ(NEQ),NV_NEQ(NEQ)
    END DO
    WRITE(6,*) 'NS   ','NSA  ','NV   ','NEQ  '
    DO NSA=1,NSAMAX
       DO NV=1,3
          WRITE(6,'(4I5)') NS_NSA(NSA),NSA,NV,NEQ_NVNSA(NV,NSA)
       END DO
    END DO   

!   *** radial mesh ***

    DR=RA/NRMAX

    RG(0)=0.D0
    DO NR=1,NRMAX
       RM(NR)=(NR-0.5D0)*DR
       RG(NR)=NR*DR
    END DO

    DO NR=1,NRMAX
       DVRHO(NR)=4.D0*PI*PI*RR*RKAP*RM(NR)
    END DO

!   *** initial profile ***

    ALLOCATE(PLF(NSMAX))
    DO NR=1,NRMAX
       RHON=RM(NR)/RA
       CALL pl_prof(RHON,plf)
       DO NSA=1,NSAMAX
          NS=NS_NSA(NSA)
          IF(ID_NS(NS).LE.2) THEN
             RNA(NSA,NR)=PLF(NS)%RN
          ELSE
             RNA(NSA,NR)=0.D0
          END IF
          RTA(NSA,NR)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
          RUA(NSA,NR)=PLF(NS)%RU
       END DO
    END DO
    DEALLOCATE(PLF)

    NT=0
    T=0.D0
    NTGT=0
    NTGR=0
    
    RETURN
  END SUBROUTINE ti_prep
END MODULE tiprep

