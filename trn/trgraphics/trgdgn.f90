MODULE trgdgn

  USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,neqmax,neqrmax, &
       neq_neqr,nsa_neq,nva_neq,rhog
  USE libgrf, ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_diagnostic

  CHARACTER(LEN=30) :: label

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg ! (1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &     ! (0:nrmax,nsamax)
       nrd1g,nrd2g,nrd3g,nrd4g
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &     ! (1:nrmax,nsamax)
       nrd1mg,nrd2mg,nrd3mg,nrd4mg

CONTAINS

  SUBROUTINE tr_gr_diagnostic(k2)
! -------------------------------------------------------------------------
!        Control routine of outputs for diagnostic and debug
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rhom

    CHARACTER(LEN=1),INTENT(IN) :: k2
    INTEGER(ikind) :: ierr,iosts,i2

    CALL tr_gr_diagnostic_alloc(ierr)

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)

    READ(k2,'(I1)',IOSTAT=iosts) i2

    IF(iosts /= 0) THEN
       WRITE(6,*) ' ERROR : Unsupported graoh ID'
       RETURN
    END IF

    SELECT CASE(i2)
    CASE(1)
       CALL tr_gr_dgn1
    CASE(2)
       CALL tr_gr_dgn2
    END SELECT

    RETURN
  END SUBROUTINE tr_gr_diagnostic

! *************************************************************************
  SUBROUTINE tr_gr_dgn1
    USE trcomm, ONLY: nrd1,nrd2,nrd3,nrd4

    CALL tr_gr_dgn_init_nrdg
    CALL tr_gr_dgn_init_nrdmg

    !--- for diagnostic array
!    nrd1mg(1:nrmax,1) = nrd1(1:nrmax)
    nrd1mg(1:nrmax,1) = nrd1(1:nrmax)
    nrd2mg(1:nrmax,1) = nrd2(1:nrmax)
    nrd3mg(1:nrmax,1) = nrd3(1:nrmax)
    nrd4mg(1:nrmax,1) = nrd4(1:nrmax)

    nrd1g(0:nrmax,1) = nrd1(0:nrmax)
    nrd1g(0:nrmax,2) = nrd2(0:nrmax)
    nrd2g(0:nrmax,1) = nrd2(0:nrmax)
    nrd3g(0:nrmax,1) = nrd3(0:nrmax)
    nrd4g(0:nrmax,1) = nrd4(0:nrmax)

    nrd1mg(1:nrmax,2) = nrd2(1:nrmax)

    CALL PAGES
    LABEL = '/diagnostic1 vs rho/'
    CALL GRD1D(1,rhog,nrd1g, nrmax+1, nrmax+1, 2, label, 0)
    LABEL = '/diagnostic2 vs rho/'
    CALL GRD1D(2,rhog,nrd2g, nrmax+1, nrmax+1, 1, label, 0)
    LABEL = '/diagnostic3 vs rho/'
    CALL GRD1D(3,rhog,nrd3g, nrmax+1, nrmax+1, 1, label, 0)
    LABEL = '/diagnostic4 vs rho/'
    CALL GRD1D(4,rhomg,nrd4mg, nrmax, nrmax, 1, label, 0)
    CALL PAGEE    

  END SUBROUTINE tr_gr_dgn1

! *************************************************************************
  SUBROUTINE tr_gr_dgn2
! -------------------------------------------------------------------------
!            Confirmation of 1-D metric quantity
! -------------------------------------------------------------------------
    USE trcomm, ONLY: &
         rjcb,ar1rho,ar2rho,abrho,rmjrho,rmnrho,rkprho,epsrho, &
         pvolrho,psurrho,dvrho,abb1rho,abb2rho,aib2rho,ttrho,  &
         arrho,abvrho,dpdrho, psiprho
    IMPLICIT NONE
    REAL(rkind) :: deriv3
    INTEGER(ikind) :: nr

    CALL tr_gr_dgn_init_nrdg

    DO nr = 0, nrmax
       nrd1g(nr,1) = deriv3(nr,rhog,psiprho,nrmax,0)
    END DO

    CALL PAGES
    label ='/ar1rho vs rho/'
    CALL GRD1D(5,rhog,ar1rho,nrmax+1,nrmax+1,1,label,0)
    label ='/abrho vs rho/'
    CALL GRD1D(6,rhog,abrho,nrmax+1,nrmax+1,1,label,0)
    label ='/pvolrho vs rho/'
    CALL GRD1D(7,rhog,pvolrho,nrmax+1,nrmax+1,1,label,0)
    label ='/dvrho vs rho/'
    CALL GRD1D(8,rhog,dvrho,nrmax+1,nrmax+1,1,label,0)
    label ='/abb1rho vs rho/'
    CALL GRD1D(9,rhog,abb1rho,nrmax+1,nrmax+1,1,label,0)
    label ='/ttrho vs rho/'
    CALL GRD1D(10,rhog,ttrho,nrmax+1,nrmax+1,1,label,0)
    label ='/arrho vs rho/'
    CALL GRD1D(11,rhog,arrho,nrmax+1,nrmax+1,1,label,0)
    label ='/abvrho vs rho/'
    CALL GRD1D(12,rhog,nrd1g,nrmax+1,nrmax+1,1,label,0)
    label ='/dpdrho vs rho/'
    CALL GRD1D(13,rhog,dpdrho,nrmax+1,nrmax+1,1,label,0)
    CALL PAGEE


    RETURN
  END SUBROUTINE tr_gr_dgn2

! *************************************************************************
! *************************************************************************
! *************************************************************************
  SUBROUTINE tr_gr_diagnostic_alloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: nrmax_save,nsamax_save

    IF(nrmax /= nrmax_save .OR. nsamax /= nsamax_save)THEN

       IF(nrmax_save /= 0) CALL tr_gr_diagnostic_dealloc

       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ALLOCATE(nrd1g(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd2g(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd3g(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd4g(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd1mg(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd2mg(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd3mg(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(nrd4mg(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT

          nrmax_save  = nrmax
          nsamax_save = nsamax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_diagnostic_alloc: allocation error: ierr=',ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_diagnostic_alloc

  SUBROUTINE tr_gr_diagnostic_dealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)
    IF(ALLOCATED(nrd1g)) DEALLOCATE(nrd1g)
    IF(ALLOCATED(nrd2g)) DEALLOCATE(nrd2g)
    IF(ALLOCATED(nrd3g)) DEALLOCATE(nrd3g)
    IF(ALLOCATED(nrd4g)) DEALLOCATE(nrd4g)
    IF(ALLOCATED(nrd1mg)) DEALLOCATE(nrd1mg)
    IF(ALLOCATED(nrd2mg)) DEALLOCATE(nrd2mg)
    IF(ALLOCATED(nrd3mg)) DEALLOCATE(nrd3mg)
    IF(ALLOCATED(nrd4mg)) DEALLOCATE(nrd4mg)

    RETURN
  END SUBROUTINE tr_gr_diagnostic_dealloc

! ************************************************************************

  SUBROUTINE tr_gr_dgn_init_nrdg

    nrd1g(0:nrmax,1:nsamax) = 0.d0
    nrd2g(0:nrmax,1:nsamax) = 0.d0
    nrd3g(0:nrmax,1:nsamax) = 0.d0
    nrd4g(0:nrmax,1:nsamax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_dgn_init_nrdg

  SUBROUTINE tr_gr_dgn_init_nrdmg

    nrd1mg(0:nrmax,1:nsamax) = 0.d0
    nrd2mg(0:nrmax,1:nsamax) = 0.d0
    nrd3mg(0:nrmax,1:nsamax) = 0.d0
    nrd4mg(0:nrmax,1:nsamax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_dgn_init_nrdmg

END MODULE trgdgn
