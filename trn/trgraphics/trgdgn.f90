MODULE trgdgn
! -------------------------------------------------------------------------
!  module for debugging or temporary checks of variables
! -------------------------------------------------------------------------
  USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,neq_neqr,nsa_neq,nva_neq,rhog
  USE trgsub,ONLY: tr_gr_time,tr_gr_vnr_alloc,tr_gr_init_vg,tr_gr_init_vm, &
       vg1,vg2,vg3,vg4, vm1,vm2,vm3,vm4, rhomg
  USE libgrf, ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_diagnostic

  CHARACTER(LEN=30) :: label
  INTEGER(ikind)    :: idexp

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

    CALL tr_gr_vnr_alloc(ierr)

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)

    READ(k2,'(I1)',IOSTAT=iosts) i2

    IF(iosts /= 0) THEN
       WRITE(6,*) ' ERROR : Unsupported graoh ID'
       RETURN
    END IF

!    idexp = 0 ! print simulation time on every GSAF page

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

    CALL tr_gr_init_vg
    CALL tr_gr_init_vm

    !--- for diagnostic array
    vg1(0:nrmax,1) = nrd1(0:nrmax)
    vg2(0:nrmax,1) = nrd2(0:nrmax)
    vg3(0:nrmax,1) = nrd3(0:nrmax)
    vg4(0:nrmax,1) = nrd4(0:nrmax)

    vm1(1:nrmax,1) = nrd1(1:nrmax)
    vm2(1:nrmax,1) = nrd2(1:nrmax)
    vm3(1:nrmax,1) = nrd3(1:nrmax)
    vm4(1:nrmax,1) = nrd4(1:nrmax)

    CALL PAGES
    LABEL = '/diagnostic1 vs rho/'
    CALL GRD1D(1,rhog,vg1, nrmax+1, nrmax+1, 1, label, 0)
    LABEL = '/diagnostic2 vs rho/'
    CALL GRD1D(2,rhog,vg2, nrmax+1, nrmax+1, 1, label, 0)
    LABEL = '/diagnostic3 vs rho/'
    CALL GRD1D(3,rhog,vg3, nrmax+1, nrmax+1, 1, label, 0)
    LABEL = '/diagnostic4 vs rho/'
    CALL GRD1D(4,rhomg,vm4, nrmax, nrmax, 1, label, 0)
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
    label ='/rmjrho vs rho/'
    CALL GRD1D(12,rhog,rmjrho,nrmax+1,nrmax+1,1,label,0)
    label ='/dpdrho vs rho/'
    CALL GRD1D(13,rhog,dpdrho,nrmax+1,nrmax+1,1,label,0)
    CALL PAGEE


    RETURN
  END SUBROUTINE tr_gr_dgn2

END MODULE trgdgn
