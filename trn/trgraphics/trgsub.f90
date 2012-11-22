MODULE trgsub
! ------------------------------------------------------------------------
!   common variables and subroutines for graphic output
! ------------------------------------------------------------------------
  USE trcomm, ONLY: ikind,rkind,nrmax,nsamax,ngt,lmaxtr,ntxmax,nsum
  IMPLICIT NONE

  PUBLIC

  ! maximum number of radial profiles for temporal evolution
  INTEGER(ikind),PARAMETER :: nggmax=10

  ! radial profiles
  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg !(1:nrmax) ! for axis
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,nsamax)
       &       vg1,vg2,vg3,vg4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,nsamax)
       &       vm1,vm2,vm3,vm4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,10)
       &       vgx1,vgx2,vgx3,vgx4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,10)
       &       vmx1,vmx2,vmx3,vmx4

  ! temporal evolution of radial profiles
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,0:nggmax)
       &       gg1,gg2,gg3,gg4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,0:nggmax)
       &       gm1,gm2,gm3,gm4

  ! temporal evolution of parameters
  REAL(rkind),DIMENSION(:),ALLOCATABLE   :: gt !(0:ngt) ! for axis
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &  !(0:ngt,nsamax)
       &       gt1,gt2,gt3,gt4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &  !(0:ngt,10)
       &       gti1,gti2,gti3,gti4

  ! numerical profile
  REAL(rkind),DIMENSION(:),ALLOCATABLE   :: lg, lg1 !(1:lmaxtr)

  ! experimental data profiles
  REAL(rkind),DIMENSION(:),ALLOCATABLE   :: rhomgu !(1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &      !(0:nrmax,1:nsum)
       &       vgu1,vgu2,vgu3,vgu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &      !(1:nrmax,1:nsum)
       &       vmu1,vmu2,vmu3,vmu4
  REAL(rkind),DIMENSION(:),ALLOCATABLE   :: gtu    !(1:ntxmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &      !(1:ntxmax,1:nsum)
       &       gtu1,gtu2,gtu3,gtu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &      !(1:ntxmax,10)
       &       gtiu1,gtiu2,gtiu3,gtiu4
  

CONTAINS

  SUBROUTINE tr_gr_time(idexp)
! ----------------------------------------------------------------------
!            Write time on GSAF pages
!
!  idexp = 0 : simulation time
!        = 1 : time of snap shot of experimental data (for mdluf = 1)
!        = 2 : time of snap shot of experimental data (for mdluf = 2)
!               (arbitrary time slice every graphic view)
! ----------------------------------------------------------------------
    USE trcomm, ONLY : t, time_slc, time_snap
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN) :: idexp

    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.0)
    CALL SETFNT(32)
    CALL MOVE(11.8,18.0)
    CALL TEXT('t =',2)

    SELECT CASE(idexp)
    CASE(0)
       CALL NUMBD(t,'(1F7.3)',7)
    CASE(1)
       CALL NUMBD(time_slc,'(1F7.3)',7)
    CASE(2)
       CALL NUMBD(time_snap,'(1F7.3)',7)
    END SELECT

    CALL TEXT(' sec.',4)

    RETURN
  END SUBROUTINE tr_gr_time

! ************************************************************************

  SUBROUTINE tr_gr_vnr_alloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: nrmax_save=0, nsamax_save=0

    ierr = 0
    IF(nrmax /= nrmax_save .OR. nsamax /= nsamax_save)THEN

       IF(nrmax_save /= 0) CALL tr_gr_vnr_dealloc
       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ALLOCATE(vg1(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vg2(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vg3(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vg4(0:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx1(0:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx2(0:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx3(0:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vgx4(0:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx1(1:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx2(1:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx3(1:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vmx4(1:nrmax,10),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm1(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm2(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm3(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(vm4(1:nrmax,nsamax),STAT=ierr); IF(ierr /=0) EXIT

          nrmax_save  = nrmax
          nsamax_save = nsamax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_vnr_alloc: allocation error: ierr= ',ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_vnr_alloc

  SUBROUTINE tr_gr_vnr_dealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)

    IF(ALLOCATED(vg1)) DEALLOCATE(vg1)
    IF(ALLOCATED(vg2)) DEALLOCATE(vg2)
    IF(ALLOCATED(vg3)) DEALLOCATE(vg3)
    IF(ALLOCATED(vg4)) DEALLOCATE(vg4)
    IF(ALLOCATED(vgx1)) DEALLOCATE(vgx1)
    IF(ALLOCATED(vgx2)) DEALLOCATE(vgx2)
    IF(ALLOCATED(vgx3)) DEALLOCATE(vgx3)
    IF(ALLOCATED(vgx4)) DEALLOCATE(vgx4)
    IF(ALLOCATED(vmx1)) DEALLOCATE(vmx1)
    IF(ALLOCATED(vmx2)) DEALLOCATE(vmx2)
    IF(ALLOCATED(vmx3)) DEALLOCATE(vmx3)
    IF(ALLOCATED(vmx4)) DEALLOCATE(vmx4)
    IF(ALLOCATED(vm1)) DEALLOCATE(vm1)
    IF(ALLOCATED(vm2)) DEALLOCATE(vm2)
    IF(ALLOCATED(vm3)) DEALLOCATE(vm3)
    IF(ALLOCATED(vm4)) DEALLOCATE(vm4)

    RETURN
  END SUBROUTINE tr_gr_vnr_dealloc


  SUBROUTINE tr_gr_vnrt_alloc(ierr)
! *** The variables for half integer axis 'rhomg' is allocated
! ***  by the above subrouitne tr_gr_vnr_alloc.
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: nrmax_save=0

    ierr = 0
    IF(nrmax /= nrmax_save)THEN
       IF(nrmax_save /= 0) CALL tr_gr_vnrt_dealloc
       DO
          ALLOCATE(gg1(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gg2(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gg3(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gg4(0:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm1(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm2(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm3(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT
          ALLOCATE(gm4(1:nrmax,0:nggmax),STAT=ierr); IF(ierr /=0) EXIT

          nrmax_save  = nrmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_vnrt_alloc: allocation error: ierr= ',ierr
    END IF

    RETURN    
  END SUBROUTINE tr_gr_vnrt_alloc

  SUBROUTINE tr_gr_vnrt_dealloc

    IF(ALLOCATED(gg1)) DEALLOCATE(gg1)
    IF(ALLOCATED(gg2)) DEALLOCATE(gg2)
    IF(ALLOCATED(gg3)) DEALLOCATE(gg3)
    IF(ALLOCATED(gg4)) DEALLOCATE(gg4)
    IF(ALLOCATED(gm1)) DEALLOCATE(gm1)
    IF(ALLOCATED(gm2)) DEALLOCATE(gm2)
    IF(ALLOCATED(gm3)) DEALLOCATE(gm3)
    IF(ALLOCATED(gm4)) DEALLOCATE(gm4)

    RETURN
  END SUBROUTINE tr_gr_vnrt_dealloc


  SUBROUTINE tr_gr_vnt_alloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: ngt_save, nsamax_save

    ierr = 0
    IF(ngt /= ngt_save .OR. nsamax /= nsamax_save)THEN

       IF(ngt_save /= 0) CALL tr_gr_vnt_dealloc
       DO
          ALLOCATE(gt(0:ngt),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt1(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt2(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt3(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gt4(0:ngt,nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gti1(0:ngt,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gti2(0:ngt,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gti3(0:ngt,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gti4(0:ngt,10),STAT=ierr); IF(ierr /= 0) EXIT

          ngt_save    = ngt
          nsamax_save = nsamax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_vnt_alloc: allocation error: ierr= ',ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_vnt_alloc

  SUBROUTINE tr_gr_vnt_dealloc

    IF(ALLOCATED(gt)) DEALLOCATE(gt)
    IF(ALLOCATED(gt1)) DEALLOCATE(gt1)
    IF(ALLOCATED(gt2)) DEALLOCATE(gt2)
    IF(ALLOCATED(gt3)) DEALLOCATE(gt3)
    IF(ALLOCATED(gt4)) DEALLOCATE(gt4)
    IF(ALLOCATED(gti1)) DEALLOCATE(gti1)
    IF(ALLOCATED(gti2)) DEALLOCATE(gti2)
    IF(ALLOCATED(gti3)) DEALLOCATE(gti3)
    IF(ALLOCATED(gti4)) DEALLOCATE(gti4)

    RETURN
  END SUBROUTINE tr_gr_vnt_dealloc

  
  SUBROUTINE tr_gr_lmt_alloc(ierr)
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: lmaxtr_save = 0

    ierr = 0
    IF(lmaxtr /= lmaxtr_save)THEN

       IF(lmaxtr_save /= 0) CALL tr_gr_lmt_dealloc
       DO
          ALLOCATE(lg(1:lmaxtr),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(lg1(1:lmaxtr),STAT=ierr); IF(ierr /= 0) EXIT

          lmaxtr_save = lmaxtr
          RETURN
       END DO
       WRITE(6,'(A,I5)') ' XX tr_gr_lmt_alloc: allocation error: ierr= ',ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_lmt_alloc

  SUBROUTINE tr_gr_lmt_dealloc

    IF(ALLOCATED(lg)) DEALLOCATE(lg)
    IF(ALLOCATED(lg1)) DEALLOCATE(lg1)

    RETURN
  END SUBROUTINE tr_gr_lmt_dealloc


  SUBROUTINE tr_gr_exp_alloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: nrmax_save=0

    IF(nrmax /= nrmax_save)THEN

       IF(nrmax_save /= 0 ) CALL tr_gr_exp_dealloc
       DO
          ALLOCATE(vgu1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT

          nrmax_save = nrmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_alloc: allocation error: ierr= ', ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_exp_alloc

  SUBROUTINE tr_gr_exp_dealloc

    IF(ALLOCATED(vgu1)) DEALLOCATE(vgu1)
    IF(ALLOCATED(vgu2)) DEALLOCATE(vgu2)
    IF(ALLOCATED(vgu3)) DEALLOCATE(vgu3)
    IF(ALLOCATED(vgu4)) DEALLOCATE(vgu4)
    IF(ALLOCATED(vmu1)) DEALLOCATE(vmu1)
    IF(ALLOCATED(vmu2)) DEALLOCATE(vmu2)
    IF(ALLOCATED(vmu3)) DEALLOCATE(vmu3)
    IF(ALLOCATED(vmu4)) DEALLOCATE(vmu4)

    RETURN
  END SUBROUTINE tr_gr_exp_dealloc


  SUBROUTINE tr_gr_expt_alloc(ierr)
! ***  The deallocation routine is not needed 
! ***   because ntxmax and nsum are not hcange
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: ntalloc_save = 0

    IF(ntalloc_save == 0)THEN
       DO
          ALLOCATE(gtu(1:ntxmax),STAT=ierr); IF(ierr /= 0) EXIT

          ALLOCATE(gtu1(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtu2(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtu3(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtu4(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu1(1:ntxmax,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu2(1:ntxmax,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu3(1:ntxmax,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu4(1:ntxmax,10),STAT=ierr); IF(ierr /= 0) EXIT

          ntalloc_save = 1
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_ntalloc: allocation error: ierr= ', ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_expt_alloc

! ************************************************************************
  
  SUBROUTINE tr_gr_init_vg

    vg1(0:nrmax,1:nsamax) = 0.d0
    vg2(0:nrmax,1:nsamax) = 0.d0
    vg3(0:nrmax,1:nsamax) = 0.d0
    vg4(0:nrmax,1:nsamax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_vg

  SUBROUTINE tr_gr_init_vm

    vm1(1:nrmax,1:nsamax) = 0.d0
    vm2(1:nrmax,1:nsamax) = 0.d0
    vm3(1:nrmax,1:nsamax) = 0.d0
    vm4(1:nrmax,1:nsamax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_vm

  SUBROUTINE tr_gr_init_vgx

    vgx1(0:nrmax,1:10) = 0.d0
    vgx2(0:nrmax,1:10) = 0.d0
    vgx3(0:nrmax,1:10) = 0.d0
    vgx4(0:nrmax,1:10) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_vgx

  SUBROUTINE tr_gr_init_vmx

    vmx1(1:nrmax,1:10) = 0.d0
    vmx2(1:nrmax,1:10) = 0.d0
    vmx3(1:nrmax,1:10) = 0.d0
    vmx4(1:nrmax,1:10) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_vmx


  SUBROUTINE tr_gr_init_gg

    gg1(0:nrmax,0:nggmax) = 0.d0
    gg2(0:nrmax,0:nggmax) = 0.d0
    gg3(0:nrmax,0:nggmax) = 0.d0
    gg4(0:nrmax,0:nggmax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_gg

  SUBROUTINE tr_gr_init_gm

    gm1(1:nrmax,0:nggmax) = 0.d0
    gm2(1:nrmax,0:nggmax) = 0.d0
    gm3(1:nrmax,0:nggmax) = 0.d0
    gm4(1:nrmax,0:nggmax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_gm

  SUBROUTINE tr_gr_init_gt

    gt1(0:ngt,1:nsamax) = 0.d0
    gt2(0:ngt,1:nsamax) = 0.d0
    gt3(0:ngt,1:nsamax) = 0.d0
    gt4(0:ngt,1:nsamax) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_gt

  SUBROUTINE tr_gr_init_gti

    gti1(0:ngt,1:10) = 0.d0
    gti2(0:ngt,1:10) = 0.d0
    gti3(0:ngt,1:10) = 0.d0
    gti4(0:ngt,1:10) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_gti

  SUBROUTINE tr_gr_init_lg

    lg1(1:lmaxtr) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_lg

  SUBROUTINE tr_gr_init_vgu

    vgu1(0:nrmax,1:nsum) = 0.d0
    vgu2(0:nrmax,1:nsum) = 0.d0
    vgu3(0:nrmax,1:nsum) = 0.d0
    vgu4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_vgu

  SUBROUTINE tr_gr_init_vmu

    vmu1(0:nrmax,1:nsum) = 0.d0
    vmu2(0:nrmax,1:nsum) = 0.d0
    vmu3(0:nrmax,1:nsum) = 0.d0
    vmu4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_vmu
  
  SUBROUTINE tr_gr_init_gtu

    gtu1(1:ntxmax,1:nsum) = 0.d0
    gtu2(1:ntxmax,1:nsum) = 0.d0
    gtu3(1:ntxmax,1:nsum) = 0.d0
    gtu4(1:ntxmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_gtu

  SUBROUTINE tr_gr_init_gtiu

    gtiu1(1:ntxmax,1:10) = 0.d0
    gtiu2(1:ntxmax,1:10) = 0.d0
    gtiu3(1:ntxmax,1:10) = 0.d0
    gtiu4(1:ntxmax,1:10) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_init_gtiu

END MODULE trgsub
