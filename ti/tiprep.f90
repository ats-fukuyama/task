! MODULE tiprep

MODULE tiprep

CONTAINS

  SUBROUTINE ti_prep(IERR)

    USE ticomm
    USE plprof
    USE tirecord
    1USE libmtx
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS,NSA,NZ,NEQ,NR,NV
    REAL(rkind):: RHON
    TYPE(pl_plf_type),DIMENSION(:),ALLOCATABLE:: PLF

    IERR=0

!   *** count NSA_MAX: number of active particles ***

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
    nsa_max=NSA

!   *** ALLOCATE array for nsa_max and NRMAX ***

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
    IF(NSA.NE.nsa_max) THEN
       WRITE(6,*) 'XX ti_prep: INCONSISTENT nsa_max'
       STOP
    END IF

!   *** Display NSA variables ***

    WRITE(6,*) 'NSA  ','NS   ','ID   ','PMA         ','PZA'
    DO NSA=1,nsa_max
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
    DO NSA=1,nsa_max
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

    NEQ_NVNSA(1:3,0:nsa_max)=0.D0
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
    DO NSA=1,nsa_max
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

    voltot=0.D0
    DO NR=1,NRMAX
       DVRHO(NR)=4.D0*PI*PI*RR*RKAP*RM(NR)
       voltot=voltot+DVRHO(NR)*DR
    END DO
    DVRHOS=4.D0*PI*PI*RR*RKAP*RG(NRMAX)

!   *** initial profile ***

    ALLOCATE(PLF(NSMAX))
    DO NR=1,NRMAX
       RHON=RM(NR)/RA
       CALL pl_prof(RHON,plf)
       DO NSA=1,nsa_max
          NS=NS_NSA(NSA)
          IF(ID_NS(NS).LE.2) THEN
             RNA(NSA,NR)=PLF(NS)%RN
          ELSE
             RNA(NSA,NR)=0.D0
          END IF
          RTA(NSA,NR)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
          RUA(NSA,NR)=PLF(NS)%RU
       END DO
!       WRITE(6,'(A,I5,1P3E12.4)') &
!            'NR,RHON,RNE,RTE=',NR,RHON,RNA(1,NR),RTA(1,NR)
    END DO
    DEALLOCATE(PLF)

    NT=0
    T=0.D0
    ngt_max=0
    ngr_max=0
    IF(ALLOCATED(gt)) DEALLOCATE(gt)
    IF(ALLOCATED(gvt)) DEALLOCATE(gvt)
    IF(ALLOCATED(gvta)) DEALLOCATE(gvta)
    ngt_allocate_max=0
    IF(ALLOCATED(grt)) DEALLOCATE(grt)
    IF(ALLOCATED(gvrt)) DEALLOCATE(gvrt)
    IF(ALLOCATED(gvrta)) DEALLOCATE(gvrta)
    ngr_allocate_max=0

    WRITE(6,*) '--- point 1'
    imax=NRMAX*NEQMAX
    jwidth=4*NEQMAX-1
    WRITE(6,'(A,I5,6I8/)') 'nrank: imax,nrmax=', &
                        nrank, imax,nrmax
    CALL mtx_setup(imax,istart,iend,jwidth)
    NR_START=(istart+NEQMAX-1)/NEQMAX
    NR_END=(iend+NEQMAX-1)/NEQMAX
    WRITE(6,*) '--- point 2'
    CALL mtx_cleanup
    WRITE(6,*) '--- point 3'

    WRITE(6,'(A,I5,6I8/)') 'nrank: imax/s/e,nrmax/s/e=', &
                        nrank, imax,istart,iend,nrmax,nr_start,nr_end

    IF(MOD(NT,NTSTEP ).EQ.0) CALL ti_snap         ! integrate and save
    IF(MOD(NT,NGTSTEP).EQ.0) CALL ti_record_ngt   ! save for time history
    IF(MOD(NT,NGRSTEP).EQ.0) CALL ti_record_ngr   ! save for radial profile

    RETURN
  END SUBROUTINE ti_prep
END MODULE tiprep

