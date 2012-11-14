MODULE trgexp
! *************************************************************************
!       graphic output for experimental data
! *************************************************************************
  USE trcomm,ONLY: ikind,rkind,ntum,nrum,nsum,nrmax,rhog,tmu, &
                    mdluf,ntxmax,tlmax
  USE trgsub,ONLY: tr_gr_time
  USE libgrf,ONLY: grd1d

  IMPLICIT NONE
  PUBLIC tr_gr_experiment

  CHARACTER(LEN=50) :: label
  INTEGER(ikind) :: nr, idexp
  INTEGER(ikind) :: ntsl ! time slice point node number for experimental data
  INTEGER(ikind) :: ntxsnap

  REAL(rkind),DIMENSION(:),ALLOCATABLE :: rhomg !(1:nrmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,1:nsum)
       vgu1,vgu2,vgu3,vgu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,1:nsum)
       vmu1,vmu2,vmu3,vmu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(0:nrmax,5)
       vgxu1,vgxu2,vgxu3,vgxu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:nrmax,5)
       vmxu1,vmxu2,vmxu3,vmxu4

  REAL(rkind),DIMENSION(:),ALLOCATABLE   :: gtu !(1:ntxmax)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:ntxmax,1:nsum)
       gtu1,gtu2,gtu3,gtu4
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &   !(1:ntxmax,7)
       gtiu1,gtiu2,gtiu3,gtiu4

CONTAINS

  SUBROUTINE tr_gr_experiment(k2,k3)
! -----------------------------------------------------------------------
!         Control routine of experimental data outputs
! -----------------------------------------------------------------------
    USE trcomm, ONLY: rhom, mdlugt, time_snap, time_slc
    USE trufsub,ONLY: tr_uf_time_slice
    CHARACTER(LEN=1),INTENT(IN) :: k2,k3
    INTEGER(ikind) :: i2,i3,ierr,iosts

    CALL tr_gr_exp_nralloc(ierr)
    IF(ierr /= 0) RETURN

    CALL tr_gr_exp_ntalloc(ierr)
    IF(ierr /= 0) RETURN

    ! set axis
    rhomg(1:nrmax) = rhom(1:nrmax)
    gtu(1:ntxmax)  = tmu(1:ntxmax)

    SELECT CASE(mdluf)
    CASE(1)
       idexp = 1
       ntxsnap = time_slc
    CASE(2)
       IF(k3 == ' ')THEN
          SELECT CASE(mdlugt)
          CASE(0)
             time_snap = -1.d0 ! initialization
          CASE(1)
             CONTINUE
          CASE(2)
             time_snap = tmu(ntxmax)
          END SELECT
          CALL tr_uf_time_slice(time_snap,tmu,ntxmax,ntsl)
          idexp = 2
          ntxsnap = ntsl
       END IF

    CASE DEFAULT
       WRITE(6,*) 'XX Experimental data has not been read. MDLUF= ',mdluf
       RETURN
    END SELECT


    READ(k2,'(I1)',IOSTAT=iosts) i2
    READ(k3,'(I1)',IOSTAT=iosts) i3
    IF(iosts /= 0)THEN
       WRITE(6,*) ' ERROR : Unsupported graph ID'
       RETURN
    END IF


    IF(k3 .EQ. ' ')THEN ! snapshot
       SELECT CASE(i2)
       CASE(1)
          CALL tr_gr_exp1
       CASE(2)
          CALL tr_gr_exp2
       CASE(3)
          CALL tr_gr_exp3
       CASE(4)
          CALL tr_gr_exp4
       CASE(9)
          CALL tr_gr_exp9
       END SELECT

    ELSE IF(i2 == 1)THEN ! time evolution
       SELECT CASE(i3)
       CASE(1)
       CASE(2)
       END SELECT

    ELSE IF(i2 == 2)THEN ! time evolution of 0D variables
       SELECT CASE(i3)
       CASE(1)
          CALL tr_gr_exp21
       CASE(2)
          CALL tr_gr_exp22
       CASE(3)
          CALL tr_gr_exp23
       END SELECT
    END IF   

    RETURN
  END SUBROUTINE tr_gr_experiment

! ************************************************************************
  SUBROUTINE tr_gr_exp1

    USE trcomm, ONLY:rnu,rtu,rpu,qpu
    INTEGER(ikind) :: nsu

    CALL tr_gr_exp_init_vgu

    DO nsu = 1, nsum
       vgu1(0:nrmax,nsu) = rnu(nsu,ntxsnap,1:nrmax+1)
       vgu2(0:nrmax,nsu) = rtu(nsu,ntxsnap,1:nrmax+1)
       vgu3(0:nrmax,nsu) = rpu(nsu,ntxsnap,1:nrmax+1)
    END DO

    vgu4(0:nrmax,1) = qpu(ntxsnap,1:nrmax+1)


    CALL PAGES
    label = '/n(exp) [10$+20$=/m$+3$=] vs rho/'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,nsum,label,0)
    label = '/T(exp) [keV] vs rho/'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,nsum,label,0)
    label = '/p(exp) [Pa] vs rho/'
    CALL GRD1D(3,rhog,vgu3,nrmax+1,nrmax+1,nsum,label,0)
    label = '/q(exp) vs rho/'
    CALL GRD1D(4,rhog,vgu4,nrmax+1,nrmax+1,1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp1

! ************************************************************************
  SUBROUTINE tr_gr_exp2
    ! other profile: zeff ...
!    USE trcomm, ONLY:

    RETURN
  END SUBROUTINE tr_gr_exp2

! ************************************************************************
  SUBROUTINE tr_gr_exp3
    ! heating density (1)
    USE trcomm, ONLY: qnbu,qecu,qicu,qlhu,qohmu,qradu,qfusu

    CALL tr_gr_exp_init_vgu

    ! for electron
    vgu1(0:nrmax,1) = qnbu(1,ntxsnap,1:nrmax+1) *1.d-6
    vgu1(0:nrmax,2) = qecu(1,ntxsnap,1:nrmax+1) *1.d-6
    vgu1(0:nrmax,3) = qicu(1,ntxsnap,1:nrmax+1) *1.d-6
    vgu1(0:nrmax,4) = qlhu(1,ntxsnap,1:nrmax+1) *1.d-6
    ! for ion
    vgu3(0:nrmax,1) = qnbu(2,ntxsnap,1:nrmax+1) *1.d-6
    vgu3(0:nrmax,2) = qecu(2,ntxsnap,1:nrmax+1) *1.d-6
    vgu3(0:nrmax,3) = qicu(2,ntxsnap,1:nrmax+1) *1.d-6
    vgu3(0:nrmax,4) = qlhu(2,ntxsnap,1:nrmax+1) *1.d-6
    ! ohm, rad
    vgu2(0:nrmax,1) = qohmu(ntxsnap,1:nrmax+1) *1.d-6
    vgu2(0:nrmax,2) = qradu(ntxsnap,1:nrmax+1) *1.d-6
    ! fusion
    vgu4(0:nrmax,1) = qfusu(1,ntxsnap,1:nrmax+1) *1.d-6
    vgu4(0:nrmax,2) = qfusu(2,ntxsnap,1:nrmax+1) *1.d-6

    CALL PAGES
    label = '/Pe_nb,ec,ic,lh(exp) [MW/m$+3$=] vs rho/'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,4,label,0)
    label = '/Pe_ohm,rad(exp) [MW/m$+3$=] vs rho/'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,2,label,0)
    label = '/Pi_nb,ec,ic,lh(exp) [MW/m$+3$=] vs rho/'
    CALL GRD1D(3,rhog,vgu3,nrmax+1,nrmax+1,4,label,0)
    label = '/P(e,i)_fus [MW/m$+3$=] vs rho/'
    CALL GRD1D(4,rhog,vgu4,nrmax+1,nrmax+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp3

! ************************************************************************
  SUBROUTINE tr_gr_exp4
    ! current density and particle source density
    USE trcomm, ONLY: jtotu,jbsu,jnbu,jecu,jicu,jlhu,snbu,swallu
    REAL(rkind),DIMENSION(1:ntxmax,1:nrmax+1) :: johu

    CALL tr_gr_exp_init_vgu

    DO nr = 1, nrmax+1
       johu(1:ntxmax,nr) =  jtotu(1:ntxmax,nr) &
                          -(jbsu(1:ntxmax,nr) + jnbu(1:ntxmax,nr) &
                           +jecu(1:ntxmax,nr) + jicu(1:ntxmax,nr) &
                           +jlhu(1:ntxmax,nr))
    END DO

    vgu1(0:nrmax,1) = - jtotu(ntxsnap,1:nrmax+1) * 1.d-6
    vgu1(0:nrmax,2) = - johu(ntxsnap,1:nrmax+1) * 1.d-6
    vgu1(0:nrmax,3) = - jbsu(ntxsnap,1:nrmax+1) * 1.d-6
    vgu1(0:nrmax,4) = -(jnbu(ntxsnap,1:nrmax+1)+jecu(ntxsnap,1:nrmax+1) &
                      + jicu(ntxsnap,1:nrmax+1)+jlhu(ntxsnap,1:nrmax+1))*1.d-6

    vgu2(0:nrmax,1) = - jnbu(ntxsnap,1:nrmax+1) *1.d-6
    vgu2(0:nrmax,2) = - jecu(ntxsnap,1:nrmax+1) *1.d-6
    vgu2(0:nrmax,3) = - jicu(ntxsnap,1:nrmax+1) *1.d-6
    vgu2(0:nrmax,4) = - jlhu(ntxsnap,1:nrmax+1) *1.d-6

    vgu3(0:nrmax,1) = snbu(1,ntxsnap,1:nrmax+1) *1.d-20
    vgu3(0:nrmax,2) = snbu(2,ntxsnap,1:nrmax+1) *1.d-20

    vgu4(0:nrmax,1) = swallu(ntxsnap,1:nrmax+1) *1.d-20


    CALL PAGES
    label = '/jtot,joh,jbs,jcd(exp) [MA] vs rho/'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,4,label,0)
    label = '/j_cd (NB,EC,IC,LH) [MA] vs rho/'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,4,label,0)
    label = '/S_nb(e,i) [10$+20$=/m$+3$= s] vs rho/'
    CALL GRD1D(3,rhog,vgu3,nrmax+1,nrmax+1,2,label,0)
    label = '/Swall [10$+20$=/m$+3$= s] vs rho/'
    CALL GRD1D(4,rhog,vgu4,nrmax+1,nrmax+1,1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp4

! ************************************************************************
  SUBROUTINE tr_gr_exp9
    ! quasi neutrality check of ufile data
    USE trufcalc, ONLY: sumzni,rnuc,rnfuc,zeffruc,ns_mion,ns_mimp
    USE trcomm, ONLY: rnu,zeffru

    INTEGER(ikind) :: nsu
    
    CALL tr_gr_exp_init_vgu
    
    ! original electron density profile
    vgu1(0:nrmax,1) = rnu(1,ntxsnap,1:nrmax+1)
    ! re-calculated profile of sum of ions density
    vgu1(0:nrmax,2) = sumzni(ntxsnap,1:nrmax+1)

    vgu2(0:nrmax,1) = rnuc(ns_mion,ntxsnap,1:nrmax+1)
    vgu2(0:nrmax,2) = rnuc(ns_mimp,ntxsnap,1:nrmax+1)
    vgu2(0:nrmax,3) = rnu(ns_mion,ntxsnap,1:nrmax+1)
    vgu2(0:nrmax,4) = rnu(ns_mimp,ntxsnap,1:nrmax+1)

    DO nsu = 1, nsum-1
       ! original ion density profiles
       vgu3(0:nrmax,nsu) = rnuc(nsu+1,ntxsnap,1:nrmax+1)
    END DO

    ! the original effective charge number profile
    vgu4(0:nrmax,1) = zeffruc(ntxsnap,1:nrmax+1)
    ! re-calculated profile of the effective charge number
    vgu4(0:nrmax,2) = zeffru(ntxsnap,1:nrmax+1)

    CALL PAGES
    label = '/neutrality (re-calc) (n_e and Sum_i (Z_i n_i))/'
    CALL GRD1D(1,rhog,vgu1,nrmax+1,nrmax+1,2,label,0)
    label = '/neutrality (compare) (n_mion,n_mimp)/'
    CALL GRD1D(2,rhog,vgu2,nrmax+1,nrmax+1,4,label,0)
    label = '/n_i(exp) vs rho/'
    CALL GRD1D(3,rhog,vgu3,nrmax+1,nrmax+1,nsum,label,0)
    label = '/Zeff (compare) (org,re-calc) vs rho/'
    CALL GRD1D(4,rhog,vgu4,nrmax+1,nrmax+1,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp9

! ************************************************************************
  SUBROUTINE tr_gr_exp21
    ! values at the axis

    USE trcomm, ONLY: rnu,rtu,rpu,qpu
    INTEGER(ikind) :: nsu

    CALL tr_gr_exp_init_gtu

    DO nsu = 1, nsum
       gtu1(1:ntxmax,nsu) = rnu(nsu,1:ntxmax,1)
       gtu2(1:ntxmax,nsu) = rpu(nsu,1:ntxmax,1)
       gtu3(1:ntxmax,nsu) = rtu(nsu,1:ntxmax,1)
    END DO
    gtu4(1:ntxmax,1) = qpu(1:ntxmax,1)
    gtu4(1:ntxmax,2) = qpu(1:ntxmax,nrmax)

    CALL PAGES
    label = '/n0(exp) [10$+20$=/m$+3$=] vs t/'
    CALL GRD1D(1,gtu,gtu1,ntxmax,ntxmax,nsum,label,0)
    label = '/p0(exp) [Pa] vs t/'
    CALL GRD1D(2,gtu,gtu2,ntxmax,ntxmax,nsum,label,0)
    label = '/T0(exp) [keV] vs t/'
    CALL GRD1D(3,gtu,gtu3,ntxmax,ntxmax,nsum,label,0)
    label = '/q0,qa(exp) vs t/'
    CALL GRD1D(4,gtu,gtu4,ntxmax,ntxmax,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp21

! ************************************************************************
  SUBROUTINE tr_gr_exp22
    USE trcomm,ONLY: ripu,wtotu,wthu

    CALL tr_gr_exp_init_gtiu

    gtiu1(1:ntxmax,1) = - ripu(1:ntxmax) * 1.d-6
    gtiu2(1:ntxmax,1) = wtotu(1:ntxmax) * 1.d-6
    gtiu2(1:ntxmax,2) = wthu(1:ntxmax) * 1.d-6

    CALL PAGES
    label = '/Ipl(exp) [MA] vs t/'
    CALL GRD1D(1,gtu,gtiu1,ntxmax,ntxmax,1,label,0)
    label = '/Wp,Wth [MJ] vs t/'
    CALL GRD1D(2,gtu,gtiu2,ntxmax,ntxmax,2,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp22

! ************************************************************************
  SUBROUTINE tr_gr_exp23
    ! source and sink power
    USE trcomm, ONLY: pnbu,pecu,pibwu,picu,plhu,pohmu,pradu

    CALL tr_gr_exp_init_gtiu

    gtiu1(1:ntxmax,1) = pnbu(1:ntxmax)  * 1.d-6

    gtiu2(1:ntxmax,1) = pecu(1:ntxmax)  * 1.d-6
    gtiu2(1:ntxmax,2) = pibwu(1:ntxmax) * 1.d-6
    gtiu2(1:ntxmax,3) = picu(1:ntxmax)  * 1.d-6
    gtiu2(1:ntxmax,4) = plhu(1:ntxmax)  * 1.d-6

    gtiu3(1:ntxmax,1) = pohmu(1:ntxmax) * 1.d-6
    gtiu4(1:ntxmax,1) = pradu(1:ntxmax) * 1.d-6

    CALL PAGES
    label = '/P(NBI) [MW] vs t/'
    CALL GRD1D(1,gtu,gtiu1,ntxmax,ntxmax,1,label,0)
    label = '/P(EC,IBW,IC,LH) [MW] vs t/'
    CALL GRD1D(2,gtu,gtiu2,ntxmax,ntxmax,7,label,0)
    label = '/P(OHM) [MW] vs t/'
    CALL GRD1D(3,gtu,gtiu3,ntxmax,ntxmax,1,label,0)
    label = '/P(RAD) [MW] vs t/'
    CALL GRD1D(4,gtu,gtiu4,ntxmax,ntxmax,1,label,0)

    CALL tr_gr_time(idexp)
    CALL PAGEE

    RETURN
  END SUBROUTINE tr_gr_exp23

! ************************************************************************
! ************************************************************************
! ************************************************************************

  SUBROUTINE tr_gr_exp_nralloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: nrmax_save=0

    IF(nrmax /= nrmax_save)THEN
       
       IF(nrmax_save /= 0 ) CALL tr_gr_exp_nrdealloc

       DO
          ALLOCATE(rhomg(1:nrmax),STAT=ierr); IF(ierr /= 0) EXIT

          ALLOCATE(vgu1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgu4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu1(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu2(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu3(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmu4(0:nrmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu1(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu2(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu3(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vgxu4(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu1(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu2(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu3(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vmxu4(0:nrmax,5),STAT=ierr); IF(ierr /= 0) EXIT

          nrmax_save = nrmax
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_nralloc: allocation error: ierr= ', ierr

    END IF
    RETURN
  END SUBROUTINE tr_gr_exp_nralloc


  SUBROUTINE tr_gr_exp_ntalloc(ierr)

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind),SAVE :: ntalloc_save = 0

    IF(ntalloc_save == 0)THEN
       DO
          ALLOCATE(gtu(1:ntxmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ALLOCATE(gtu1(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtu2(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtu3(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtu4(1:ntxmax,1:nsum),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu1(1:ntxmax,7),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu2(1:ntxmax,7),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu3(1:ntxmax,7),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gtiu4(1:ntxmax,7),STAT=ierr); IF(ierr /= 0) EXIT
          
          ntalloc_save = 1
          RETURN
       END DO
       WRITE(6,*) ' XX tr_gr_exp_ntalloc: allocation error: ierr= ', ierr
    END IF

    RETURN
  END SUBROUTINE tr_gr_exp_ntalloc


  SUBROUTINE tr_gr_exp_nrdealloc

    IF(ALLOCATED(rhomg)) DEALLOCATE(rhomg)
    IF(ALLOCATED(vgu1)) DEALLOCATE(vgu1)
    IF(ALLOCATED(vgu2)) DEALLOCATE(vgu2)
    IF(ALLOCATED(vgu3)) DEALLOCATE(vgu3)
    IF(ALLOCATED(vgu4)) DEALLOCATE(vgu4)
    IF(ALLOCATED(vmu1)) DEALLOCATE(vmu1)
    IF(ALLOCATED(vmu2)) DEALLOCATE(vmu2)
    IF(ALLOCATED(vmu3)) DEALLOCATE(vmu3)
    IF(ALLOCATED(vmu4)) DEALLOCATE(vmu4)
    IF(ALLOCATED(vgxu1)) DEALLOCATE(vgxu1)
    IF(ALLOCATED(vgxu2)) DEALLOCATE(vgxu2)
    IF(ALLOCATED(vgxu3)) DEALLOCATE(vgxu3)
    IF(ALLOCATED(vgxu4)) DEALLOCATE(vgxu4)
    IF(ALLOCATED(vmxu1)) DEALLOCATE(vmxu1)
    IF(ALLOCATED(vmxu2)) DEALLOCATE(vmxu2)
    IF(ALLOCATED(vmxu3)) DEALLOCATE(vmxu3)
    IF(ALLOCATED(vmxu4)) DEALLOCATE(vmxu4)

    RETURN
  END SUBROUTINE tr_gr_exp_nrdealloc

! ***********************************************************************

  SUBROUTINE tr_gr_exp_init_vgu

    vgu1(0:nrmax,1:nsum) = 0.d0
    vgu2(0:nrmax,1:nsum) = 0.d0
    vgu3(0:nrmax,1:nsum) = 0.d0
    vgu4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_exp_init_vgu

  SUBROUTINE tr_gr_exp_init_vmu

    vmu1(0:nrmax,1:nsum) = 0.d0
    vmu2(0:nrmax,1:nsum) = 0.d0
    vmu3(0:nrmax,1:nsum) = 0.d0
    vmu4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_exp_init_vmu

  SUBROUTINE tr_gr_exp_init_vgxu

    vgxu1(0:nrmax,1:nsum) = 0.d0
    vgxu2(0:nrmax,1:nsum) = 0.d0
    vgxu3(0:nrmax,1:nsum) = 0.d0
    vgxu4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_exp_init_vgxu

  SUBROUTINE tr_gr_exp_init_vmxu

    vmxu1(0:nrmax,1:nsum) = 0.d0
    vmxu2(0:nrmax,1:nsum) = 0.d0
    vmxu3(0:nrmax,1:nsum) = 0.d0
    vmxu4(0:nrmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_exp_init_vmxu

  SUBROUTINE tr_gr_exp_init_gtu

    gtu1(1:ntxmax,1:nsum) = 0.d0
    gtu2(1:ntxmax,1:nsum) = 0.d0
    gtu3(1:ntxmax,1:nsum) = 0.d0
    gtu4(1:ntxmax,1:nsum) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_exp_init_gtu

  SUBROUTINE tr_gr_exp_init_gtiu

    gtiu1(1:ntxmax,1:7) = 0.d0
    gtiu2(1:ntxmax,1:7) = 0.d0
    gtiu3(1:ntxmax,1:7) = 0.d0
    gtiu4(1:ntxmax,1:7) = 0.d0

    RETURN
  END SUBROUTINE tr_gr_exp_init_gtiu


END MODULE trgexp
