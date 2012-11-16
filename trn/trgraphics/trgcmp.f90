MODULE trgcmp
! ***********************************************************************
!    graphic output for comparison of calculation result and exp. data
! ***********************************************************************
  USE trcomm,ONLY: ikind,rkind,nsum,ngt,nrmax,rhog

  USE trgsub,ONLY: tr_gr_time
  USE libgrf,ONLY: grd1d

  IMPLICIT NONE
  PUBLIC tr_gr_comparison

  CHARACTER(LEN=50) :: label
  INTEGER(ikind) :: nr, idexp

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg !(1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,1:nsum)
       vgc1,vgc2,vgc3,vgc4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,1:nsum)
       vmc1,vmc2,vmc3,vmc4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,1:5)
       vgxc1,vgxc2,vgxc3,vgxc4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,1:5)
       vmxc1,vmxc2,vmxc3,vmxc4
  REAL(rkind),DIMENSION(:),ALLOCATABLE   :: gtc !(0:ngt)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:ngt,1:nsum)
       gtc1,gtc2,gtc3,gtc4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:ngt,1:7)
       gtic1,gtic2,gtic3,gtic4

CONTAINS

  SUBROUTINE tr_gr_comparison(k2,k3)
! -----------------------------------------------------------------------
!    Control routine of data comparison outputs
! -----------------------------------------------------------------------
    USE trcomm,ONLY: rhom,gvt

    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,ierr,iosts

    CALL tr_gr_cmp_nralloc(ierr)
    IF(ierr /= 0) RETURN

    CALL tr_gr_cmp_ntalloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)
    gtc(0:ngt)     = gvt(0:ngt,0)
    
    READ(k2,'(I1)',IOSTAT=iosts) i2
    READ(k3,'(I1)',IOSTAT=iosts) i3
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    idexp = 0

    IF(k3 .EQ. ' ')THEN
       SELECT CASE(i2)
       CASE(1)
          CALL tr_gr_cmp1
       CASE(2)
          CALL tr_gr_cmp2
       CASE(3)
          CALL tr_gr_cmp3
       END SELECT

    ELSE IF(i2 == 1)THEN
       SELECT CASE(i3)
       CASE(1)
       CASE(2)
       END SELECT

    ELSE IF(i2 == 2)THEN ! time evolution of 0D variables
       SELECT CASE(i3)
       CASE(1)
          CALL tr_gr_cmp21
       CASE(2)
          CALL tr_gr_cmp22
       END SELECT

    END IF

    RETURN
  END SUBROUTINE tr_gr_comparison

! **********************************************************************

  SUBROUTINE tr_gr_cmp1
! ----------------------------------------------------------------------
!   bulk species density profiles
! ----------------------------------------------------------------------
    USE trcomm, ONLY: nsamax,idnsa,nsa_nsu,nsab_nsaf,rn
    USE trufin, ONLY: rnug
    
    INTEGER(ikind) :: nsa,ns,nsu

    CALL tr_gr_cmp_init_vgxc

    ns = 0
    DO nsu = 1, nsum
       nsa = nsa_nsu(nsu)
       IF(nsa == 0) CYCLE
       IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2) CYCLE
       ns  = ns + 1

       IF(ns==1)THEN
          vgxc1(0:nrmax,1) = rn(nsa,0:nrmax)
          vgxc1(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       ELSE IF(ns==2)THEN
          vgxc2(0:nrmax,1) = rn(nsa,0:nrmax)
          vgxc2(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       ELSE IF(ns==3)THEN
          vgxc3(0:nrmax,1) = rn(nsa,0:nrmax)
          vgxc3(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       ELSE IF(ns==4)THEN
          vgxc4(0:nrmax,1) = rn(nsa,0:nrmax)
          vgxc4(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       END IF
       IF(ns > 4) EXIT
    END DO

    CALL PAGES
    label = '@Ne, Ne(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(1,rhog,vgxc1,nrmax+1,nrmax+1,2,label,0)
    label = '@N1, N1(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(2,rhog,vgxc2,nrmax+1,nrmax+1,2,label,0)
    label = '@N2, N2(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(3,rhog,vgxc3,nrmax+1,nrmax+1,2,label,0)
    label = '@N3, N3(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(4,rhog,vgxc4,nrmax+1,nrmax+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_cmp1

! ***********************************************************************

  SUBROUTINE tr_gr_cmp2
! ----------------------------------------------------------------------
!   bulk species temperature profiles
! ----------------------------------------------------------------------
    USE trcomm, ONLY: nsamax,idnsa,nsa_nsu,nsab_nsaf,rt
    USE trufin, ONLY: rtug
    
    INTEGER(ikind) :: nsa,ns,nsu

    CALL tr_gr_cmp_init_vgxc

    ns = 0
    DO nsu = 1, nsum
       nsa = nsa_nsu(nsu)
       IF(nsa == 0) CYCLE
       IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2) CYCLE
       ns  = ns + 1

       IF(ns==1)THEN
          vgxc1(0:nrmax,1) = rt(nsa,0:nrmax)
          vgxc1(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       ELSE IF(ns==2)THEN
          vgxc2(0:nrmax,1) = rt(nsa,0:nrmax)
          vgxc2(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       ELSE IF(ns==3)THEN
          vgxc3(0:nrmax,1) = rt(nsa,0:nrmax)
          vgxc3(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       ELSE IF(ns==4)THEN
          vgxc4(0:nrmax,1) = rt(nsa,0:nrmax)
          vgxc4(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       END IF
       IF(ns > 4) EXIT
    END DO

    CALL PAGES
    label = '@Te, Te(exp) [keV] vs rho@'
    CALL GRD1D(1,rhog,vgxc1,nrmax+1,nrmax+1,2,label,0)
    label = '@T1, T1(exp) [keV] vs rho@'
    CALL GRD1D(2,rhog,vgxc2,nrmax+1,nrmax+1,2,label,0)
    label = '@T2, T2(exp) [keV] vs rho@'
    CALL GRD1D(3,rhog,vgxc3,nrmax+1,nrmax+1,2,label,0)
    label = '@T3, T3(exp) [keV] vs rho@'
    CALL GRD1D(4,rhog,vgxc4,nrmax+1,nrmax+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_cmp2

  
  SUBROUTINE tr_gr_cmp3
    USE trcomm,ONLY: jtot,joh,jbs_nc,jcd_nb,jcd_ec,jcd_ic,jcd_lh
    USE trufin,ONLY: jtotug,jnbug,jecug,jicug,jlhug,jbsug

    REAL(rkind),DIMENSION(1:nrmax+1) :: johug

    johug(1:nrmax+1) = jtotug(1:nrmax+1) - jbsug(1:nrmax+1) &
       -(jnbug(1:nrmax+1)+jecug(1:nrmax+1)+jicug(1:nrmax+1)*jecug(1:nrmax+1))

    CALL tr_gr_cmp_init_vgxc

    vgxc1(0:nrmax,1) = jtot(0:nrmax)   * 1.d-6
    vgxc1(0:nrmax,2) = - jtotug(1:nrmax+1) * 1.d-6

    vgxc2(0:nrmax,1) = joh(0:nrmax)    * 1.d-6
    vgxc2(0:nrmax,2) = jbs_nc(0:nrmax) * 1.d-6
    vgxc2(0:nrmax,3) = - johug(1:nrmax+1) * 1.d-6
    vgxc2(0:nrmax,4) = - jbsug(1:nrmax+1) * 1.d-6

    vgxc3(0:nrmax,1) = jcd_nb(0:nrmax) * 1.d-6
    vgxc3(0:nrmax,2) = jcd_ec(0:nrmax) * 1.d-6
    vgxc3(0:nrmax,3) = jcd_ic(0:nrmax) * 1.d-6
    vgxc3(0:nrmax,4) = jcd_lh(0:nrmax) * 1.d-6

    vgxc4(0:nrmax,1) = - jnbug(1:nrmax+1) * 1.d-6
    vgxc4(0:nrmax,2) = - jecug(1:nrmax+1) * 1.d-6
    vgxc4(0:nrmax,3) = - jicug(1:nrmax+1) * 1.d-6
    vgxc4(0:nrmax,4) = - jlhug(1:nrmax+1) * 1.d-6

    CALL PAGES
    label = '@jtot (sim,exp) [MA/m^2] vs rho@'
    CALL GRD1D(1,rhog,vgxc1,nrmax+1,nrmax+1,2,label,0)
    label = '@joh,jbs (sim,exp) [MA/m^2] vs rho@'
    CALL GRD1D(2,rhog,vgxc2,nrmax+1,nrmax+1,4,label,0)
    label = '@jcd_nb,ec,ic,lh (sim) [MA/m^2] vs rho@'
    CALL GRD1D(3,rhog,vgxc3,nrmax+1,nrmax+1,4,label,0)
    label = '@jcd_nb,ec,ic,lh (exp) [MA/m^2] vs rho@'
    CALL GRD1D(4,rhog,vgxc4,nrmax+1,nrmax+1,4,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_cmp3


  SUBROUTINE tr_gr_cmp21
    USE trcomm,ONLY: gvt,gvtu

    CALL tr_gr_cmp_init_gtic

    gtic1(0:ngt,1) = gvt(0:ngt,3)   ! rip
    gtic1(0:ngt,2) = gvtu(0:ngt,3)  ! ripu

    gtic2(0:ngt,1) = gvt(0:ngt,1)  ! qp(0)
    gtic2(0:ngt,2) = gvt(0:ngt,2)  ! qp(a)
    gtic2(0:ngt,3) = gvtu(0:ngt,1) ! qpu(0)
    gtic2(0:ngt,4) = gvtu(0:ngt,2) ! qpu(a)

    gtic3(0:ngt,1) = gvt(0:ngt,7)        ! w_th
    gtic3(0:ngt,2) = gvtu(0:ngt,4)*1.d-6 ! w_thu

    gtic4(0:ngt,1) = gvt(0:ngt,8)        ! w_tot
    gtic4(0:ngt,2) = gvtu(0:ngt,5)*1.d-6 ! w_totu   

    CALL PAGES
    label = '@RIP (sim,exp) [MA] vs t@'
    CALL GRD1D(1,gtc,gtic1,ngt+1,ngt+1,2,label,0)
    label = '@q0,qa (sim,exp) vs t@'
    CALL GRD1D(2,gtc,gtic2,ngt+1,ngt+1,4,label,0)
    label = '@Wth (sim,exp) [MJ] vs t@'
    CALL GRD1D(3,gtc,gtic3,ngt+1,ngt+1,2,label,0)
    label = '@Wtot (sim,exp) [MJ] vs t@'
    CALL GRD1D(4,gtc,gtic4,ngt+1,ngt+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_cmp21


  SUBROUTINE tr_gr_cmp22

    RETURN
  END SUBROUTINE tr_gr_cmp22

! ***********************************************************************

  SUBROUTINE tr_gr_cmp_nralloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: nrmax_save=0

    ierr = 0
    IF(nrmax /= nrmax_save)THEN
       IF(nrmax_save /= 0 ) CALL tr_gr_cmp_nrdealloc
       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgc1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgc2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgc3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgc4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmc1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmc2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmc3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmc4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxc1(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxc2(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxc3(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxc4(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxc1(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxc2(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxc3(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxc4(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          nrmax_save = nrmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_nralloc: allocation error: ierr= ', ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_cmp_nralloc
     

  SUBROUTINE tr_gr_cmp_ntalloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: ngt_save = -1

    ierr = 0
    IF(ngt /= ngt_save)THEN
    IF(ngt_save /= -1) CALL tr_gr_cmp_ntdealloc
       DO
          ALLOCATE(gtc(0:ngt),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtc1(0:ngt,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtc2(0:ngt,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtc3(0:ngt,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtc4(0:ngt,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtic1(0:ngt,7),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtic2(0:ngt,7),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtic3(0:ngt,7),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtic4(0:ngt,7),STAT=ierr); IF(ierr /= 0) EXIT
          ngt_save = ngt
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_ntalloc: allocation error: ierr= ', ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_cmp_ntalloc

  SUBROUTINE tr_gr_cmp_nrdealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)
    IF(ALLOCATED(vgc1)) DEALLOCATE(vgc1)
    IF(ALLOCATED(vgc2)) DEALLOCATE(vgc2)
    IF(ALLOCATED(vgc3)) DEALLOCATE(vgc3)
    IF(ALLOCATED(vgc4)) DEALLOCATE(vgc4)
    IF(ALLOCATED(vmc1)) DEALLOCATE(vmc1)
    IF(ALLOCATED(vmc2)) DEALLOCATE(vmc2)
    IF(ALLOCATED(vmc3)) DEALLOCATE(vmc3)
    IF(ALLOCATED(vmc4)) DEALLOCATE(vmc4)
    IF(ALLOCATED(vgxc1)) DEALLOCATE(vgxc1)
    IF(ALLOCATED(vgxc2)) DEALLOCATE(vgxc2)
    IF(ALLOCATED(vgxc3)) DEALLOCATE(vgxc3)
    IF(ALLOCATED(vgxc4)) DEALLOCATE(vgxc4)
    IF(ALLOCATED(vmxc1)) DEALLOCATE(vmxc1)
    IF(ALLOCATED(vmxc2)) DEALLOCATE(vmxc2)
    IF(ALLOCATED(vmxc3)) DEALLOCATE(vmxc3)
    IF(ALLOCATED(vmxc4)) DEALLOCATE(vmxc4)
    
    RETURN
  END SUBROUTINE tr_gr_cmp_nrdealloc

  SUBROUTINE tr_gr_cmp_ntdealloc

    IF(ALLOCATED(gtc)) DEALLOCATE(gtc)
    IF(ALLOCATED(gtc1)) DEALLOCATE(gtc1)
    IF(ALLOCATED(gtc2)) DEALLOCATE(gtc2)
    IF(ALLOCATED(gtc3)) DEALLOCATE(gtc3)
    IF(ALLOCATED(gtc4)) DEALLOCATE(gtc4)
    IF(ALLOCATED(gtic1)) DEALLOCATE(gtic1)
    IF(ALLOCATED(gtic2)) DEALLOCATE(gtic2)
    IF(ALLOCATED(gtic3)) DEALLOCATE(gtic3)
    IF(ALLOCATED(gtic4)) DEALLOCATE(gtic4)

    RETURN
  END SUBROUTINE tr_gr_cmp_ntdealloc

! ***********************************************************************

  SUBROUTINE tr_gr_cmp_init_vgc

    vgc1(0:nrmax,1:nsum) = 0.d0
    vgc2(0:nrmax,1:nsum) = 0.d0
    vgc3(0:nrmax,1:nsum) = 0.d0
    vgc4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_cmp_init_vgc

  SUBROUTINE tr_gr_cmp_init_vmc

    vmc1(1:nrmax,1:nsum) = 0.d0
    vmc2(1:nrmax,1:nsum) = 0.d0
    vmc3(1:nrmax,1:nsum) = 0.d0
    vmc4(1:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_cmp_init_vmc

  SUBROUTINE tr_gr_cmp_init_vgxc

    vgxc1(0:nrmax,1:5) = 0.d0
    vgxc2(0:nrmax,1:5) = 0.d0
    vgxc3(0:nrmax,1:5) = 0.d0
    vgxc4(0:nrmax,1:5) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_cmp_init_vgxc

  SUBROUTINE tr_gr_cmp_init_vmxc

    vmxc1(0:nrmax,1:5) = 0.d0
    vmxc2(0:nrmax,1:5) = 0.d0
    vmxc3(0:nrmax,1:5) = 0.d0
    vmxc4(0:nrmax,1:5) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_cmp_init_vmxc

  SUBROUTINE tr_gr_cmp_init_gtc

    gtc1(0:ngt,1:nsum) = 0.d0
    gtc2(0:ngt,1:nsum) = 0.d0
    gtc3(0:ngt,1:nsum) = 0.d0
    gtc4(0:ngt,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_cmp_init_gtc

  SUBROUTINE tr_gr_cmp_init_gtic

    gtic1(0:ngt,1:7) = 0.d0
    gtic2(0:ngt,1:7) = 0.d0
    gtic3(0:ngt,1:7) = 0.d0
    gtic4(0:ngt,1:7) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_cmp_init_gtic

END MODULE trgcmp
