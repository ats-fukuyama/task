! wrprep.f90

MODULE wrprep

  PRIVATE
  PUBLIC wr_prep

CONTAINS
  
!     ***** Set equilibrium profile *****

  SUBROUTINE wr_prep(IERR)

    USE wrcomm_parm
    USE dpprep,ONLY: dp_prep_ns
    USE equnit
    USE plload
    USE plprof
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: INITEQ=0
    INTEGER:: nsu
    REAL(rkind):: dth
    EXTERNAL EQCALQ,EQGETB  !,eqget_rzsu
    INTEGER:: nrr,nrrmax_wr
    REAL(rkind):: rr_wr,drr_wr
    TYPE(pl_mag_type):: mag_wr

    IERR=0

    IF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN ! task/eq and eqdsk
       IF(INITEQ.EQ.0) THEN
          CALL eq_load(MODELG,KNAMEQ,IERR)
          IF(IERR.EQ.0) THEN
             CALL EQCALQ(IERR)
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             CALL pl_rzsu(rsu_wr,zsu_wr,nsumax)
             rmax_eq=rsu_wr(1)
             rmin_eq=rsu_wr(1)
             zmax_eq=zsu_wr(1)
             zmin_eq=zsu_wr(1)
             DO nsu=2,nsumax
                IF(rsu_wr(nsu).GT.rmax_eq) rmax_eq=rsu_wr(nsu)
                IF(rsu_wr(nsu).LT.rmin_eq) rmin_eq=rsu_wr(nsu)
                IF(zsu_wr(nsu).GT.zmax_eq) zmax_eq=zsu_wr(nsu)
                IF(zsu_wr(nsu).LT.zmin_eq) zmin_eq=zsu_wr(nsu)
             END DO
             INITEQ=1
          ELSE
             WRITE(6,*) 'XX EQLOAD: IERR=',IERR
             INITEQ=0
          ENDIF
       ENDIF
    ELSE IF(MODELG.EQ.8) THEN ! task/equ
       IF(INITEQ.EQ.0) THEN
          CALL eq_read(IERR)
          IF(IERR.EQ.0) THEN
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             CALL pl_rzsu(rsu_wr,zsu_wr,nsumax)
             rmax_eq=rsu_wr(1)
             rmin_eq=rsu_wr(1)
             zmax_eq=zsu_wr(1)
             zmin_eq=zsu_wr(1)
             DO nsu=2,nsumax
                IF(rsu_wr(nsu).GT.rmax_eq) rmax_eq=rsu_wr(nsu)
                IF(rsu_wr(nsu).LT.rmin_eq) rmin_eq=rsu_wr(nsu)
                IF(zsu_wr(nsu).GT.zmax_eq) zmax_eq=rsu_wr(nsu)
                IF(zsu_wr(nsu).LT.zmin_eq) zmin_eq=rsu_wr(nsu)
             END DO
             raxis_eq=0.5D0*(rmin_eq+rmax_eq)
             zaxis_eq=0.5D0*(zmin_eq+zmax_eq)
             INITEQ=1
          ELSE
             WRITE(6,*) 'XX EQREAD: IERR=',IERR
             INITEQ=0
          ENDIF
       ENDIF
    ELSE
       INITEQ=0
       raxis_eq=RR
       zaxis_eq=0.D0
       rmax_eq=RR+RA
       rmin_eq=RR-RA
       zmax_eq= RKAP*RA
       zmin_eq=-RKAP*RA
       nsumax=128
       IF(ALLOCATED(rsu_wr)) DEALLOCATE(rsu_wr)
       IF(ALLOCATED(zsu_wr)) DEALLOCATE(zsu_wr)
       ALLOCATE(rsu_wr(nsumax+1),zsu_wr(nsumax+1))
       dth=2.D0*PI/nsumax
       DO nsu=1,nsumax+1
          rsu_wr(nsu)=RR+RA*COS((nsu-1)*dth)
          zsu_wr(nsu)=RKAP*RA*SIN((nsu-1)*dth)
       END DO
    ENDIF
    raxis_eq=0.5D0*(rmin_eq+rmax_eq)
    zaxis_eq=0.5D0*(zmin_eq+zmax_eq)
    IF(bdr_threshold.GT.0.D0) THEN
       rmax_wr=raxis_eq+bdr_threshold*(rmax_eq-raxis_eq)
       rmin_wr=raxis_eq+bdr_threshold*(rmin_eq-raxis_eq)
       zmax_wr=zaxis_eq+bdr_threshold*(zmax_eq-zaxis_eq)
       zmin_wr=zaxis_eq+bdr_threshold*(zmin_eq-zaxis_eq)
    END IF
    IF(rmin_wr.LT.0.D0) rmin_wr=0.D0
    WRITE(6,'(A,4ES12.4)') '_eq:',rmin_eq,rmax_eq,zmin_eq,zmax_eq
    WRITE(6,'(A,2ES12.4)') '_ax:',raxis_eq,zaxis_eq
    WRITE(6,'(A,4ES12.4)') '_wr:',rmin_wr,rmax_wr,zmin_wr,zmax_wr

    CALL pl_load(ierr)
    
    CALL dp_prep_ns(ierr)

    nrrmax_wr=21
    drr_wr=(rmax_wr-rmin_wr)/nrrmax_wr
    DO nrr=1,nrrmax_wr
       rr_wr=rmin_wr+drr_wr*(nrr-1)
       CALL pl_mag(rr_wr,0.D0,0.D0,mag_wr)
       WRITE(6,'(A,I4,6ES12.4)') 'nrr=', &
            nrr,rr_wr,mag_wr%babs,rr_wr*mag_wr%babs, &
            mag_wr%bnx,mag_wr%bny,mag_wr%bnz
    END DO
    RETURN
  END SUBROUTINE wr_prep
END MODULE wrprep
