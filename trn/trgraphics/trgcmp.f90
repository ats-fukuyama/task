MODULE trgcmp
! ***********************************************************************
!    graphic output for comparison of calculation result and exp. data
! ***********************************************************************
  USE trcomm,ONLY: ikind,rkind,nsum,ngt,nrmax,rhog
  USE trgsub,ONLY: &
       tr_gr_time,tr_gr_vnr_alloc,tr_gr_exp_alloc,tr_gr_vnt_alloc, &
       tr_gr_init_vgx,tr_gr_init_vgu,tr_gr_init_gt,tr_gr_init_gti, &
       vgx1,vgx2,vgx3,vgx4, vgu1,vgu2,vgu3,vgu4, gt1,gt2,gt3,gt4,  &
       gti1,gti2,gti3,gti4, rhomg,gt
  USE libgrf,ONLY: grd1d

  IMPLICIT NONE
  PUBLIC tr_gr_comparison

  CHARACTER(LEN=50) :: label
  INTEGER(ikind) :: nr, idexp

CONTAINS

  SUBROUTINE tr_gr_comparison(k2,k3)
! -----------------------------------------------------------------------
!    Control routine of data comparison outputs
! -----------------------------------------------------------------------
    USE trcomm,ONLY: rhom,gvt

    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,ierr,iosts

    CALL tr_gr_vnr_alloc(ierr)
    IF(ierr /= 0) RETURN
    CALL tr_gr_vnt_alloc(ierr)
    IF(ierr /= 0) RETURN
    CALL tr_gr_exp_alloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)
    gt(0:ngt)      = gvt(0:ngt,0)
    
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

    CALL tr_gr_init_vgu

    ns = 0
    DO nsu = 1, nsum
       nsa = nsa_nsu(nsu)
       IF(nsa == 0) CYCLE
       IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2) CYCLE
       ns  = ns + 1

       IF(ns==1)THEN
          vgu1(0:nrmax,1) = rn(nsa,0:nrmax)
          vgu1(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       ELSE IF(ns==2)THEN
          vgu2(0:nrmax,1) = rn(nsa,0:nrmax)
          vgu2(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       ELSE IF(ns==3)THEN
          vgu3(0:nrmax,1) = rn(nsa,0:nrmax)
          vgu3(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       ELSE IF(ns==4)THEN
          vgu4(0:nrmax,1) = rn(nsa,0:nrmax)
          vgu4(0:nrmax,2) = rnug(nsu,1:nrmax+1)
       END IF
       IF(ns > 4) EXIT
    END DO

    CALL PAGES
    label = '@Ne, Ne(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,2,label,0)
    label = '@N1, N1(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,2,label,0)
    label = '@N2, N2(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(3,rhog,vgu3,nrmax+1,nrmax+1,2,label,0)
    label = '@N3, N3(exp) [10$+20$=/m$+3$=] vs rho@'
    CALL GRD1D(4,rhog,vgu4,nrmax+1,nrmax+1,2,label,0)

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

    CALL tr_gr_init_vgu

    ns = 0
    DO nsu = 1, nsum
       nsa = nsa_nsu(nsu)
       IF(nsa == 0) CYCLE
       IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2) CYCLE
       ns  = ns + 1

       IF(ns==1)THEN
          vgu1(0:nrmax,1) = rt(nsa,0:nrmax)
          vgu1(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       ELSE IF(ns==2)THEN
          vgu2(0:nrmax,1) = rt(nsa,0:nrmax)
          vgu2(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       ELSE IF(ns==3)THEN
          vgu3(0:nrmax,1) = rt(nsa,0:nrmax)
          vgu3(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       ELSE IF(ns==4)THEN
          vgu4(0:nrmax,1) = rt(nsa,0:nrmax)
          vgu4(0:nrmax,2) = rtug(nsu,1:nrmax+1)
       END IF
       IF(ns > 4) EXIT
    END DO

    CALL PAGES
    label = '@Te, Te(exp) [keV] vs rho@'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,2,label,0)
    label = '@T1, T1(exp) [keV] vs rho@'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,2,label,0)
    label = '@T2, T2(exp) [keV] vs rho@'
    CALL GRD1D(3,rhog,vgu3,nrmax+1,nrmax+1,2,label,0)
    label = '@T3, T3(exp) [keV] vs rho@'
    CALL GRD1D(4,rhog,vgu4,nrmax+1,nrmax+1,2,label,0)

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

    CALL tr_gr_init_vgx

    vgx1(0:nrmax,1) = jtot(0:nrmax)   * 1.d-6
    vgx1(0:nrmax,2) = - jtotug(1:nrmax+1) * 1.d-6

    vgx2(0:nrmax,1) = joh(0:nrmax)    * 1.d-6
    vgx2(0:nrmax,2) = jbs_nc(0:nrmax) * 1.d-6
    vgx2(0:nrmax,3) = - johug(1:nrmax+1) * 1.d-6
    vgx2(0:nrmax,4) = - jbsug(1:nrmax+1) * 1.d-6

    vgx3(0:nrmax,1) = jcd_nb(0:nrmax) * 1.d-6
    vgx3(0:nrmax,2) = jcd_ec(0:nrmax) * 1.d-6
    vgx3(0:nrmax,3) = jcd_ic(0:nrmax) * 1.d-6
    vgx3(0:nrmax,4) = jcd_lh(0:nrmax) * 1.d-6

    vgx4(0:nrmax,1) = - jnbug(1:nrmax+1) * 1.d-6
    vgx4(0:nrmax,2) = - jecug(1:nrmax+1) * 1.d-6
    vgx4(0:nrmax,3) = - jicug(1:nrmax+1) * 1.d-6
    vgx4(0:nrmax,4) = - jlhug(1:nrmax+1) * 1.d-6

    CALL PAGES
    label = '@jtot (sim,exp) [MA/m^2] vs rho@'
    CALL GRD1D(1,rhog,vgx1,nrmax+1,nrmax+1,2,label,0)
    label = '@joh,jbs (sim,exp) [MA/m^2] vs rho@'
    CALL GRD1D(2,rhog,vgx2,nrmax+1,nrmax+1,4,label,0)
    label = '@jcd_nb,ec,ic,lh (sim) [MA/m^2] vs rho@'
    CALL GRD1D(3,rhog,vgx3,nrmax+1,nrmax+1,4,label,0)
    label = '@jcd_nb,ec,ic,lh (exp) [MA/m^2] vs rho@'
    CALL GRD1D(4,rhog,vgx4,nrmax+1,nrmax+1,4,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_cmp3


  SUBROUTINE tr_gr_cmp21
    USE trcomm,ONLY: gvt,gvtu

    CALL tr_gr_init_gti

    gti1(0:ngt,1) = gvt(0:ngt,3)   ! rip
    gti1(0:ngt,2) = gvtu(0:ngt,3)  ! ripu

    gti2(0:ngt,1) = gvt(0:ngt,1)  ! qp(0)
    gti2(0:ngt,2) = gvt(0:ngt,2)  ! qp(a)
    gti2(0:ngt,3) = gvtu(0:ngt,1) ! qpu(0)
    gti2(0:ngt,4) = gvtu(0:ngt,2) ! qpu(a)

    gti3(0:ngt,1) = gvt(0:ngt,7)        ! w_th
    gti3(0:ngt,2) = gvtu(0:ngt,4)*1.d-6 ! w_thu

    gti4(0:ngt,1) = gvt(0:ngt,8)        ! w_tot
    gti4(0:ngt,2) = gvtu(0:ngt,5)*1.d-6 ! w_totu   

    CALL PAGES
    label = '@RIP (sim,exp) [MA] vs t@'
    CALL GRD1D(1,gt,gti1,ngt+1,ngt+1,2,label,0)
    label = '@q0,qa (sim,exp) vs t@'
    CALL GRD1D(2,gt,gti2,ngt+1,ngt+1,4,label,0)
    label = '@Wth (sim,exp) [MJ] vs t@'
    CALL GRD1D(3,gt,gti3,ngt+1,ngt+1,2,label,0)
    label = '@Wtot (sim,exp) [MJ] vs t@'
    CALL GRD1D(4,gt,gti4,ngt+1,ngt+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_cmp21


  SUBROUTINE tr_gr_cmp22

    RETURN
  END SUBROUTINE tr_gr_cmp22


END MODULE trgcmp
