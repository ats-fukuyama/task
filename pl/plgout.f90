! MODULE plgout

MODULE plgout

  PRIVATE
  PUBLIC pl_gout

CONTAINS

  SUBROUTINE pl_gout

    USE plcomm,ONLY: NSMAX,rkind,pl_prf_type
    USE plprof
    USE libgrf
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: rho
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: rn,rtpr,rtpp,ru
    TYPE(pl_prf_type),DIMENSION(NSMAX):: plf
    INTEGER:: nrmax,nr,ns
    REAL(rkind):: drho

    nrmax=101
    drho=1.D0/nrmax

    ALLOCATE(rho(nrmax))
    ALLOCATE(rn(nrmax,nsmax),ru(nrmax,nsmax))
    ALLOCATE(rtpr(nrmax,nsmax),rtpp(nrmax,nsmax))

    DO nr=1,nrmax
       rho(nr)=(nr-1)*drho
       CALL pl_prof(rho(nr),plf)
       DO ns=1,nsmax
          rn(nr,ns)=plf(ns)%rn
          rtpr(nr,ns)=plf(ns)%rtpr
          rtpp(nr,ns)=plf(ns)%rtpp
          ru(nr,ns)=plf(ns)%ru
       END DO
    END DO

    CALL PAGES
    CALL GRD1D(1,rho,rn,nrmax,nrmax,nsmax,'@n vs rho@')
    CALL GRD1D(2,rho,rtpr,nrmax,nrmax,nsmax,'@T_para vs rho@')
    CALL GRD1D(3,rho,ru,nrmax,nrmax,nsmax,'@u vs rho@')
    CALL GRD1D(4,rho,rtpp,nrmax,nrmax,nsmax,'@T_pepr vs rho@')
    CALL PAGEE

    DEALLOCATE(rn,ru,rtpr,rtpp,rho)

    RETURN

  END SUBROUTINE pl_gout
END MODULE plgout
