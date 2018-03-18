! MODULE plgout

MODULE plgout

  PRIVATE
  PUBLIC pl_gout

CONTAINS

  SUBROUTINE pl_gout

    USE plcomm,ONLY: NSMAX
    USE plprof
    USE libgrf
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: rho
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: rn,rt,ru
    TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
    INTEGER:: nrmax,nr,ns
    REAL(rkind):: drho

    nrmax=101
    drho=1.D0/nrmax

    ALLOCATE(rho(nrmax))
    ALLOCATE(rn(nrmax,nsmax),rt(nrmax,nsmax),ru(nrmax,nsmax))

    DO nr=1,nrmax
       rho(nr)=(nr-1)*drho
       CALL pl_prof(rho(nr),plf)
       DO ns=1,nsmax
          rn(nr,ns)=plf(ns)%rn
          rt(nr,ns)=(plf(ns)%rtpr+2.D0*plf(ns)%rtpp)/3.D0
          ru(nr,ns)=plf(ns)%ru
       END DO
    END DO

    CALL PAGES
    CALL GRD1D(1,rho,rn,nrmax,nrmax,nsmax,'@n vs rho@')
    CALL GRD1D(2,rho,rt,nrmax,nrmax,nsmax,'@T vs rho@')
    CALL GRD1D(3,rho,ru,nrmax,nrmax,nsmax,'@u vs rho@')
    CALL PAGEE

    DEALLOCATE(rn,ru,rt,rho)

    RETURN

  END SUBROUTINE pl_gout
END MODULE plgout
