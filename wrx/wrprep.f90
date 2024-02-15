! wrprep.f90

MODULE wrprep

  PRIVATE
  PUBLIC wr_prep

CONTAINS
  
!     ***** Set equilibrium profile *****

  SUBROUTINE wr_prep(IERR)

    USE wrcomm_parm
    USE dpparm,ONLY: dpprep_local
    USE equnit
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: INITEQ=0
    INTEGER:: nsu
    REAL(rkind):: dth
    EXTERNAL EQCALQ,EQGETB  !,eqget_rzsu
    INTERFACE
       SUBROUTINE eqget_rzsu(rsu,zsu,nsumax)
         USE task_kinds,ONLY: dp
         REAL(dp),ALLOCATABLE,INTENT(OUT):: rsu(:),zsu(:)
         INTEGER,INTENT(OUT):: nsumax
       END SUBROUTINE eqget_rzsu
    END INTERFACE

    IERR=0

    IF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN ! task/eq and eqdsk
       IF(INITEQ.EQ.0) THEN
          CALL eq_load(MODELG,KNAMEQ,IERR)
          IF(IERR.EQ.0) THEN
             WRITE(LINE,'(A,I5)') 'NRMAX =',51
             CALL eq_parm(2,LINE,IERR)
             WRITE(LINE,'(A,I5)') 'NTHMAX=',64
             CALL eq_parm(2,LINE,IERR)
             WRITE(LINE,'(A,I5)') 'NSUMAX=',64
             CALL eq_parm(2,LINE,IERR)
             CALL EQCALQ(IERR)
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             CALL eqget_rzsu(rsu_wr,zsu_wr,nsumax)
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
             CALL eqget_rzsu(rsu_wr,zsu_wr,nsumax)
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
    WRITE(6,'(A,4ES12.4)') 'rmin_wr:',rmin_wr,rmax_wr,zmin_wr,zmax_wr
 !            WRITE(6,'(A,4ES12.4)') '_eq:',rmin_eq,rmax_eq,zmin_eq,zmax_eq
 !            WRITE(6,'(A,2ES12.4)') '_ax:',raxis_eq,zaxis_eq
 !            WRITE(6,'(A,4ES12.4)') '_wr:',rmin_wr,rmax_wr,zmin_wr,zmax_wr

    CALL DPPREP_LOCAL(IERR)

    RETURN
  END SUBROUTINE wr_prep
END MODULE wrprep
