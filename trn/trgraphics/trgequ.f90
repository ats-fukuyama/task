MODULE trgequ

  USE trcomm,ONLY: ikind,rkind,nrmax,nsamax,rhog
  USE trgsub,ONLY: tr_gr_time, tr_gr1d_rad,                     &
       tr_gr_vnr_alloc, tr_gr_init_vgx,                         &
       vg1,vg2,vg3,vg4,  vm1,vm2,vm3,vm4,  vgx1,vgx2,vgx3,vgx4, &
       vmx1,vmx2,vmx3,vmx4, rhomg
  USE libgrf,ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_equilibrium

  CHARACTER(LEN=30) :: label
  INTEGER(ikind)    :: nr,idexp
  
CONTAINS

  SUBROUTINE tr_gr_equilibrium(k2)
! -------------------------------------------------------------------------
!   control routine of graphic output
!                   for equlibrium information and metric quantities
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rhom
    CHARACTER(LEN=1),INTENT(IN) :: k2
    INTEGER(ikind) :: i2,ierr,iosts

    CALL tr_gr_vnr_alloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)

    READ(k2,'(I1)',IOSTAT=iosts) i2
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    idexp = 0 ! print simulation time on every GSAF page 

    SELECT CASE(i2)
    CASE(1)
       CALL tr_gr_equ1
    CASE(2)
       CALL tr_gr_equ2
    CASE(3)
       CALL tr_gr_equ3
    CASE(4)
       CALL tr_gr_equ4
    CASE(5)
       CALL tr_gr_equ5
    END SELECT

    RETURN
  END SUBROUTINE tr_gr_equilibrium

! ************************************************************************

  SUBROUTINE tr_gr_equ1
    ! show magnetic flux surface (equi-psi surface)
    USE trcomm,ONLY: mdlgmt

    IF(.NOT.(mdlgmt==8 .OR. mdlgmt==9))THEN
       WRITE(6,'(1X,A,I3)') 'Equilibrium code has not been called. MDLGMT = ',MDLGMT
       RETURN
    END IF

    RETURN
  END SUBROUTINE tr_gr_equ1


  SUBROUTINE tr_gr_equ2
!   essential metric quantities for transport equations
    USE trcomm, ONLY: qp,dvrho,ar1rho,ar2rho, pi,rkap,ra,rr,rhog

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = qp(0:nrmax)
    vgx2(0:nrmax,1) = dvrho(0:nrmax)
    vgx3(0:nrmax,1) = ar1rho(0:nrmax)
    vgx4(0:nrmax,1) = ar2rho(0:nrmax)

    CALL PAGES
    label = '@q vs rho@'
    CALL tr_gr1d_rad(1,rhog,vgx1,nrmax+1,1,label,0,FMIN0=0.d0)
    label = "@V' vs rho@"
    CALL tr_gr1d_rad(2,rhog,vgx2,nrmax+1,2,label,0)
    label = '@<|grad rho|> vs rho@'
    CALL tr_gr1d_rad(3,rhog,vgx3,nrmax+1,1,label,0,FMIN0=0.d0)
    label = '@<|grad rho|$+2$=> vs rho@'
    CALL tr_gr1d_rad(4,rhog,vgx4,nrmax+1,1,label,0,FMIN0=0.d0)
    CALL tr_gr_time(idexp)
    CALL PAGEE    

    RETURN
  END SUBROUTINE tr_gr_equ2


  SUBROUTINE tr_gr_equ3
!   essential metric quantities for magnetic diffusion equation
    USE trcomm, ONLY: qp,arrho,ttrho,abb1rho

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = qp(0:nrmax)
    vgx2(0:nrmax,1) = arrho(0:nrmax)
    vgx3(0:nrmax,1) = abb1rho(0:nrmax)
    vgx4(0:nrmax,1) = ttrho(0:nrmax)

    CALL PAGES
    label = '@q vs rho@'
    CALL tr_gr1d_rad(1,rhog,vgx1,nrmax+1,1,label,0,FMIN0=0.d0)
    label = "@<1/R$+2$=>@"
    CALL tr_gr1d_rad(2,rhog,vgx2,nrmax+1,1,label,0)
    label = '@<B> vs rho@'
    CALL tr_gr1d_rad(3,rhog,vgx3,nrmax+1,1,label,0,FMIN0=0.d0)
    label = '@I=R*B vs rho@'
    CALL tr_gr1d_rad(4,rhog,vgx4,nrmax+1,1,label,0,FMIN0=0.d0)
    CALL tr_gr_time(idexp)
    CALL PAGEE    

    RETURN
  END SUBROUTINE tr_gr_equ3


  SUBROUTINE tr_gr_equ4
!   essential metric quantities for values associate with poloidal flux
    USE trcomm, ONLY: jtot,dpdrho,rdpvrho,abvrho

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = jtot(0:nrmax) * 1.d-6
    vgx2(0:nrmax,1) = dpdrho(0:nrmax)
    vgx3(0:nrmax,1) = rdpvrho(0:nrmax)
    vgx4(0:nrmax,1) = abvrho(0:nrmax)

    CALL PAGES
    label = '@j_tot [MA/m$+2$=] vs rho@'
    CALL tr_gr1d_rad(1,rhog,vgx1,nrmax+1,1,label,0)
    label = '@d psi/d rho vs rho@'
    CALL tr_gr1d_rad(2,rhog,vgx2,nrmax+1,1,label,0)
    label = '@d psi/d V vs rho@'
    CALL tr_gr1d_rad(3,rhog,vgx3,nrmax+1,1,label,0)
    label = '@<|grad V|$+2$=/R$+2$=>@'
    CALL tr_gr1d_rad(4,rhog,vgx4,nrmax+1,1,label,0)
    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_equ4


  SUBROUTINE tr_gr_equ5
    USE trcomm, ONLY: pvolrho,psurrho,rmjrho,rmnrho

    vgx1(0:nrmax,1) = pvolrho(0:nrmax)
    vgx2(0:nrmax,1) = psurrho(0:nrmax)
    vgx3(0:nrmax,1) = rmjrho(0:nrmax)
    vgx4(0:nrmax,1) = rmnrho(0:nrmax)

    CALL PAGES
    label = '@V_plasma [m$+3$=] vs rho@'
    CALL tr_gr1d_rad(1,rhog,vgx1,nrmax+1,1,label,0)
    label = '@S_plasma [m$+3$=] vs rho@'
    CALL tr_gr1d_rad(2,rhog,vgx2,nrmax+1,1,label,0)
    label = '@R_major [m] vs rho@'
    CALL tr_gr1d_rad(3,rhog,vgx3,nrmax+1,1,label,0,FMIN0=0.d0)
    label = '@R_minor [m] vs rho@'
    CALL tr_gr1d_rad(4,rhog,vgx4,nrmax+1,1,label,0,FMIN0=0.d0)
    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_equ5

END MODULE trgequ
