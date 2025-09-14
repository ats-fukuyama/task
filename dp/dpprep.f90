! dpprep.f90

MODULE dpprep

  PRIVATE
  PUBLIC dp_prep_fp
  PUBLIC dp_prep_ns

CONTAINS

!     ***** Setup velocity distribution function *****

  SUBROUTINE dp_prep_fp(NTHMAX_DP_1,NRMAX_DP_1,RMIN_1,RMAX_1,IERR)

    USE dpcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NTHMAX_DP_1,NRMAX_DP_1
    REAL(rkind),INTENT(IN):: RMIN_1,RMAX_1
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NS

    NTHMAX_DP=NTHMAX_DP_1
    NRMAX_DP=NRMAX_DP_1
    DO NS=1,NSMAX
       RHON_MIN(NS)=RMIN_1
       RHON_MAX(NS)=RMAX_1
    END DO
    CALL dp_prep_ns(ierr)
    RETURN
  END SUBROUTINE dp_prep_fp

!     ***** Setup velocity distribution function *****

  SUBROUTINE dp_prep_ns(ierr)

    USE dpcomm
    USE dpfpin
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: nsamax_fp,nsamax_fm,ns,nsa

    IERR=0

    nsamax_fp=0
    nsamax_fm=0
    DO ns=1,nsmax
       IF(modelv(ns).EQ.2.OR.modelv(ns).EQ.4) nsamax_fp=nsamax_fp+1
       IF(modelv(ns).EQ.1.OR.modelv(ns).EQ.3) nsamax_fm=nsamax_fm+1
    END DO

    IF(nsamax_fp.GT.0.AND.nsamax_fm.GT.0) THEN
       WRITE(6,'(A,2I4)') &
            'XX dp_prep_ns: Either fp or fm: nsamax_fp,nsamax_fm=', &
            nsamax_fp,nsamax_fm
       STOP
    END IF

    IF(nsamax_fp.GT.0) THEN
       CALL DPLDFP(IERR)  ! set nsamax_dp and allocated
       IF(IERR.NE.0) RETURN
       IF(nsamax_dp.NE.nsamax_fp) THEN
          WRITE(6,*) 'XX dpprep_fp: nsamax_dp from file is not nsamax_fp'
          STOP
       END IF
    END IF

    IF(nsamax_fm.GT.0) THEN
       nsamax_dp=MAX(nsamax_fm,nsmax)
       CALL dpfp_allocate  ! nsamax_dp=nsamax_fm and allocated
       nsa=0
       DO ns=1,nsmax
          nsa=nsa+1
          ns_nsa_dp(nsa)=ns
          nsa_ns_dp(ns)=nsa
          SELECT CASE(modelv(ns))
          CASE(1)  ! non-relativistic
             CALL DPLDFM(ns,0,ierr)
             IF(ierr.NE.0) THEN
                WRITE(6,*) 'XX dpprep: dpldfm(non rel) error: ierr=',ierr
                STOP
             END IF
          CASE(3) ! relativistic
             CALL DPLDFM(ns,1,ierr)
             IF(ierr.NE.0) THEN
                WRITE(6,*) 'XX dpprep: dpldfm(rel) error: ierr=',ierr
                STOP
             END IF
          END SELECT
       END DO
    ELSE
       nsamax_dp=nsmax
       CALL dpfp_allocate  ! nsamax_dp=nsamax_fm and allocated
    END IF

    IF(nsamax_dp.LT.nsmax) THEN
       DO ns=nsamax_dp+1,nsmax
          nsa=ns
          ns_nsa_dp(nsa)=ns
          nsa_ns_dp(ns)=nsa
       END DO
       nsamax_dp=nsmax
    END IF
    RETURN
  END SUBROUTINE dp_prep_ns
END MODULE dpprep

