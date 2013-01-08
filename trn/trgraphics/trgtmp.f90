MODULE trgtmp
! **************************************************************************
!           Temporal evolution outputs
! **************************************************************************
  USE trcomm,ONLY: ikind,rkind,nsamax,neqmax,neqrmax,neq_neqr,nsa_neq,ngt
  USE trgsub,ONLY: tr_gr1d_rad, tr_gr_time,           &
       tr_gr_vnr_alloc, tr_gr_vnrt_alloc,              &
       tr_gr_vnt_alloc, tr_gr_init_gt, tr_gr_init_gti, &
       gt, gt1,gt2,gt3,gt4,  gti1,gti2,gti3,gti4
  USE libgrf,ONLY: grd1d
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_gr_temporal

  CHARACTER(LEN=30) :: label

  INTEGER(ikind) :: nr,nsa,neq,neqr,idexp

CONTAINS

  SUBROUTINE tr_gr_temporal(k2)
! --------------------------------------------------------------------------
!            Control routine of temporal evolution outputs
! --------------------------------------------------------------------------
    USE trcomm, ONLY: gvt

    CHARACTER(LEN=1),INTENT(IN) :: k2
    INTEGER(ikind) :: i2,ierr,iosts

    CALL tr_gr_vnt_alloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    gt(0:ngt) = gvt(0:ngt,0)

    ! control pages
    READ(k2,'(I1)',IOSTAT=iosts) i2
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF

    idexp = 0 ! print simulation time on every GSAF page

    SELECT CASE(i2)
    CASE(1)
       CALL tr_gr_temp1 ! n(0),u(0),T(0),q(0),q(a)
    CASE(2)
       CALL tr_gr_temp2 ! I_pl,W,taue
    CASE(3)
       CALL tr_gr_temp3 ! Pin
    CASE(4)
       CALL tr_gr_temp4 ! n, T
    END SELECT
    
    RETURN
  END SUBROUTINE tr_gr_temporal

! *************************************************************************

  SUBROUTINE tr_gr_temp1
  ! ----- time evolution of (n, u, T, q)-----
    USE trcomm,ONLY: gvt,gvts

    CALL tr_gr_init_gt

    DO nsa=1,nsamax
       gt1(0:ngt,nsa)=gvts(0:ngt,nsa,1)
       gt2(0:ngt,nsa)=gvts(0:ngt,nsa,4)
       gt3(0:ngt,nsa)=gvts(0:ngt,nsa,3)
    END DO

    gt4(0:ngt,1)=gvt(0:ngt,1)
    gt4(0:ngt,2)=gvt(0:ngt,2)

    CALL PAGES
    label = '@n(0) [10^20/m^3] vs t@'
    CALL tr_gr1d_rad(1,gt,gt1,ngt+1,nsamax,label,0,FMIN0=0.d0)
    label = '@p(0) [MPa] vs t@'
    CALL tr_gr1d_rad(2,gt,gt2,ngt+1,nsamax,label,0,FMIN0=0.d0)
    label = '@T(0) [keV] vs t@'
    CALL tr_gr1d_rad(3,gt,gt3,ngt+1,nsamax,label,0,FMIN0=0.d0)
    label = '@q(0),q(a) vs t@'
    CALL tr_gr1d_rad(4,gt,gt4,ngt+1,2,label,0,FMIN0=0.d0)

    CALL tr_gr_time(idexp)
    CALL PAGEE
    
    RETURN
  END SUBROUTINE tr_gr_temp1

! *************************************************************************

  SUBROUTINE tr_gr_temp2
    ! ----- time evolution of (I, W, taue) -----
    USE trcomm, ONLY: gvt

    CALL tr_gr_init_gti

    gti1(0:ngt,1) = gvt(0:ngt,3) ! rip
!    gti1(0:ngt,2) = gvt(0:ngt,4)
!    gti1(0:ngt,3) = gvt(0:ngt,5)

    gti2(0:ngt,1) = gvt(0:ngt,8) ! wp_t

    gti3(0:ngt,1) = gvt(0:ngt,25) ! taue3
    gti3(0:ngt,2) = gvt(0:ngt,26) ! taue89
    gti3(0:ngt,3) = gvt(0:ngt,27) ! taue98

    gti4(0:ngt,1) = gvt(0:ngt,28) ! h89
    gti4(0:ngt,2) = gvt(0:ngt,29) ! h98y2


    CALL PAGES
    label = '@Ipl [MA] vs t@'
    CALL tr_gr1d_rad(1,gt,gti1,ngt+1,1,label,0,FMIN0=0.d0)
    label = '@Wp [MJ] vs t@'
    CALL tr_gr1d_rad(2,gt,gti2,ngt+1,5,label,0,FMIN0=0.d0)
    label = '@tauE,tauE89,tauE98 (H89,H98y2) vs t@'
    CALL tr_gr1d_rad(3,gt,gti3,ngt+1,5,label,0,FMIN0=0.d0)
    label = '@H89,H98y2 vs t@'
    CALL tr_gr1d_rad(4,gt,gti4,ngt+1,5,label,0)
    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_temp2

! *************************************************************************

  SUBROUTINE tr_gr_temp3
    USE trcomm, ONLY: gvt

    CALL tr_gr_init_gti

    gti1(0:ngt,1) = gvt(0:ngt,18) ! pin_t
    gti1(0:ngt,2) = gvt(0:ngt,19) ! poh_t
    gti1(0:ngt,3) = gvt(0:ngt,20) ! pnb_t
    gti1(0:ngt,5) = gvt(0:ngt,21) ! pec_t
    gti1(0:ngt,6) = gvt(0:ngt,22) ! pic_t
    gti1(0:ngt,7) = gvt(0:ngt,23) ! plh_t
    gti1(0:ngt,8) = gvt(0:ngt,24) ! pnf_t

    gti2(0:ngt,1) = gvt(0:ngt,15) ! betap(0)
    gti2(0:ngt,2) = gvt(0:ngt,16) ! betap(nrmax)
    gti2(0:ngt,3) = gvt(0:ngt,17) ! betan

    CALL PAGES
    label = '@Ptot,oh,nb,rf,ec,ic,lh,nf [MW] vs t@'
    CALL tr_gr1d_rad(1,gt,gti1,ngt+1,8,label,0,FMIN0=0.d0)
    label = '@betap(0),betap(nrmax),betan vs t@'
    CALL tr_gr1d_rad(2,gt,gti2,ngt+1,5,label,0,FMIN0=0.d0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

  END SUBROUTINE tr_gr_temp3

  SUBROUTINE tr_gr_temp4
    USE trcomm, ONLY: gvts

    RETURN
  END SUBROUTINE tr_gr_temp4

END MODULE trgtmp
