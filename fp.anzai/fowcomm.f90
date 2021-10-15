module fowcomm

  USE commpi
  USE fpcomm_parm

  implicit none

  public

  INTEGER:: nthm1          ! number of theta_m grid points
                           !      for 0 <= theta_m <= theta_pnc
  INTEGER:: nthm2          ! number of theta_m grid points
                           !      for theta_pnc <= theta_m <= theta_co_stg
  INTEGER:: nthm3          ! number of theta_m grid points
                           !      for theta_cnt_stg <= theta_m <= pi

  real(rkind),allocatable :: JI(:,:,:,:),& ! dxdydzd(vx)d(vy)d(vz) = JI * dpd(thetam)d(psim)
                             JIR(:,:,:,:)  ! for integrate over velocity space
! COM --------------------------------------------------------------------------------------------------          
  real(rkind),allocatable,dimension(:) :: psim,&                ! maximum poloidal magnetic flux in an orbit, value at half integer grid points
                                          psimg                 ! psim at integer grid points

  real(rkind),allocatable :: rhom_local(:,:,:,:,:), &
                             thetam_local(:,:,:,:,:), &
                             time_loss(:,:,:,:,:)

  real(rkind),allocatable,dimension(:,:,:,:) :: thetam,&              ! pitch angle along orbit in psim, value at half integer grid points
                                                thetamg,&             ! theta_m at integer grid points
                                                thetam_rg,&           ! theta_m for given pm(np), psimg(nr)
                                                thetam_pg             ! theta_m for given pg(np), psim(nr)

! Diffsion coefficients extend FP's difference equation ------------------------------------------------
  real(rkind),allocatable,dimension(:,:,:,:) :: Dppfow, Dptfow, Dprfow,&
                                                Dtpfow, Dttfow, Dtrfow,&
                                                Drpfow, Drtfow, Drrfow,&
                                                Fppfow,  Fthfow,  Frrfow

! equilibrium variables --------------------------------------------------------------------------------
  real(rkind):: psi0                                             ! poloidal flux at the plasma edge
  real(rkind),allocatable,dimension(:) :: Fpsi,&                 ! poloidai current at psim
                                          Bout,&                 ! abs(magnetic field) at psim when psim is outside the axis
                                          Bin,&                  ! abs(magnetic field) at psim when psim is inside the axis
                                          Fpsig,&                ! Fpsi for grid points
                                          Boutg,&                ! Bout for grid points
                                          Bing,&                 ! Bin for grid points
                                          dFdr,&                 ! dFpsi/drhom
                                          dBoutdr,&              ! dBout/drhom
                                          dBindr,&               ! dBin/drhom
                                          dpsimdr,&              ! dpsim/drhom
                                          dFgdr,&                ! dFpsig/drhom
                                          dBoutgdr,&             ! dBoutg/drhom
                                          dBingdr,&              ! dBing/drhom
                                          dpsimgdr,&             ! dpsimg/drhom
                                          safety_factor          ! dpsimg/drhom

  real(rkind),allocatable,dimension(:,:) :: Babs,&                 ! B(psip,thetap)
                                            dBdr,&                 ! dBabs/rhot
                                            dBdthp                 ! dBabs/dtheta_p
  real(rkind),allocatable,dimension(:) :: theta_p                 ! poloidal angle

! use for boundary conditions --------------------------------------------------------------------------
  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc,&         ! theta_m of pinch orbit for given pm(np) and psim(nr)
                                              theta_co_stg,&      ! theta_m of co-stagnation orbit for pm(np) and psim(nr)
                                              theta_cnt_stg,&     ! theta_m of counter-stagnation orbit for pm(np) and psim(nr)
                                              psip_pnc_point      ! poloidal flux of pinch orbit given by I = (pm(np), theta_pnc(p,psi_m), psim(nr))

  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc_pg,&
                                              theta_co_stg_pg,&
                                              theta_cnt_stg_pg,&
                                              psip_pnc_point_pg
                                              
  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc_rg,&
                                              theta_co_stg_rg,&
                                              theta_cnt_stg_rg,&
                                              psip_pnc_point_rg

  real(rkind),allocatable,dimension(:,:,:,:) :: delthm_rg, delthm_pg, delthm

  integer,allocatable,dimension(:) :: nth_stg,&        ! thetamg(nth_stg,(nsa),np,nr,nsa) = theta_stg
                                      nth_pnc          ! thetamg(nth_pnc,(nsa),np,nr,nsa) = theta_pnc

  real(rkind),allocatable,dimension(:,:,:) :: IBCflux_ratio         ! ration between the flux from thetam(nth_pnc-1) to thetam(nth_pnc) and the flux from thetam(nth_pnc-1) to thetam(nth_co_stg)
  integer,allocatable,dimension(:,:,:) ::  nr_pnc_point             ! nr of pinch point of theta_pnc(:,:,:)

  real(rkind),allocatable,dimension(:,:,:) :: rhom_pinch      ! psim of pnc orbit whose pinch point is X stg orbit with np, nr, nsa
  integer,allocatable,dimension(:,:,:) ::  nr_rhom_pinch

  real(rkind),parameter :: NO_PINCH_ORBIT = 19960610d0
  real(rkind),allocatable,dimension(:,:,:,:) :: orbital_loss

! use for bounce average -------------------------------------------------------------------------------
  type :: orbit                                                  ! quantities along orbits culcurated by TASK/OB
    integer :: nstp_max
    real(rkind),allocatable,dimension(:) :: time, &              ! time
                                            psip, &              ! poloidal magnetic flux
                                            Babs, &              ! absolute value of magnetic field
                                            costh,&              ! COS( pitch angle )
                                            sinth,&              ! SIN( pitch angle )
                                            thetap,&             ! poloidal angle
                                            F,&                  ! poloidal current
                                            r,&                  ! minor radius
                                            dpsipdr,&            ! dpsip/drho
                                            dFdr,&               ! dF/drho
                                            dBdr,&
                                            dBdthp
  end type orbit

  type(orbit),allocatable,dimension(:,:,:,:) :: orbit_p,&    ! (np,nth,nr)=(-0.5,  1 ,  1 ) to (npmax+0.5,nthmax,nrmax)
                                                orbit_th,&   ! (np,nth,nr)=(  1 ,-0.5,  1 ) to (npmax,nthmax+0.5,nrmax)
                                                orbit_r,&    ! (np,nth,nr)=(  1 ,  1 ,-0.5) to (npmax,nthmax,nrmax+0.5)
                                                orbit_m      ! (np,nth,nr)=(1,1,1) to (npmax,nthmax,nrmax)

! ------------------------------------------------------------------------------------------------------

contains

  subroutine fow_allocate
    use fpcomm, only:npmax,nthmax,nrmax,nsamax

    IF(ALLOCATED(JI)) CALL fow_deallocate

    allocate(JI          (nthmax  ,npmax  ,nrmax          ,nsamax))
    allocate(JIR         (nthmax  ,npmax  ,nrmax          ,nsamax))
    allocate(rhom_local  (nthmax  ,npmax  ,nrmax  ,nthpmax,nsamax))
    allocate(thetam_local(nthmax  ,npmax  ,nrmax  ,nthpmax,nsamax))
    allocate(time_loss   (nthmax  ,npmax  ,nrmax  ,nthpmax,nsamax))
    allocate(thetam      (nthmax  ,npmax  ,nrmax          ,nsamax))
    allocate(thetamg     (nthmax+1,npmax  ,nrmax          ,nsamax))
    allocate(thetam_pg   (nthmax  ,npmax+1,nrmax          ,nsamax))
    allocate(thetam_rg   (nthmax  ,npmax  ,nrmax+1        ,nsamax))

    allocate(psim    (nrmax  ))
    allocate(Fpsi    (nrmax  ))
    allocate(Bout    (nrmax  ))
    allocate(Bin     (nrmax  ))
    allocate(psimg   (nrmax+1))
    allocate(Fpsig   (nrmax+1))
    allocate(Boutg   (nrmax+1))
    allocate(Bing    (nrmax+1))
    allocate(dFdr    (nrmax  ))
    allocate(dBoutdr (nrmax  ))
    allocate(dBindr  (nrmax  ))
    allocate(dpsimdr (nrmax  ))
    allocate(dFgdr   (nrmax+1))
    allocate(dBoutgdr(nrmax+1))
    allocate(dBingdr (nrmax+1))
    allocate(dpsimgdr(nrmax+1))

    allocate(safety_factor(nrmax))

    allocate(Babs    (nrmax  ,nthpmax))
    allocate(dBdr    (nrmax  ,nthpmax))
    allocate(dBdthp  (nrmax  ,nthpmax))
    allocate(theta_p (        nthpmax))

    allocate(theta_pnc        (npmax  ,nrmax  ,nsamax))
    allocate(theta_co_stg     (npmax  ,nrmax  ,nsamax))
    allocate(theta_cnt_stg    (npmax  ,nrmax  ,nsamax))
    allocate(psip_pnc_point   (npmax  ,nrmax  ,nsamax))
    allocate(theta_pnc_pg     (npmax+1,nrmax  ,nsamax))
    allocate(theta_co_stg_pg  (npmax+1,nrmax  ,nsamax))
    allocate(theta_cnt_stg_pg (npmax+1,nrmax  ,nsamax))
    allocate(psip_pnc_point_pg(npmax+1,nrmax  ,nsamax))
    allocate(theta_pnc_rg     (npmax  ,nrmax+1,nsamax))
    allocate(theta_co_stg_rg  (npmax  ,nrmax+1,nsamax))
    allocate(theta_cnt_stg_rg (npmax  ,nrmax+1,nsamax))
    allocate(psip_pnc_point_rg(npmax  ,nrmax+1,nsamax))

    allocate(delthm       (nthmax  ,npmax  ,nrmax  ,nsamax))
    allocate(delthm_pg    (nthmax  ,npmax+1,nrmax  ,nsamax))
    allocate(delthm_rg    (nthmax  ,npmax  ,nrmax+1,nsamax))
    allocate(IBCflux_ratio(         npmax  ,nrmax  ,nsamax))
    allocate(nr_pnc_point (         npmax  ,nrmax  ,nsamax))
    allocate(nr_rhom_pinch(         npmax  ,nrmax  ,nsamax))
    allocate(orbital_loss (nthmax  ,npmax  ,nrmax  ,nsamax))

    allocate(nth_stg(nsamax))
    allocate(nth_pnc(nsamax))

    allocate(orbit_th(nthmax+1,npmax  ,nrmax  ,nsamax))
    allocate(orbit_p (nthmax  ,npmax+1,nrmax  ,nsamax))
    allocate(orbit_r (nthmax  ,npmax  ,nrmax+1,nsamax))
    allocate(orbit_m (nthmax  ,npmax  ,nrmax  ,nsamax))

    allocate(Dppfow(nthmax  ,npmax+1,nrmax  ,nsamax))
    allocate(Dptfow(nthmax  ,npmax+1,nrmax  ,nsamax))
    allocate(Dprfow(nthmax  ,npmax+1,nrmax  ,nsamax))
    allocate(Dtpfow(nthmax+1,npmax  ,nrmax  ,nsamax))
    allocate(Dttfow(nthmax+1,npmax  ,nrmax  ,nsamax))
    allocate(Dtrfow(nthmax+1,npmax  ,nrmax  ,nsamax))
    allocate(Drpfow(nthmax  ,npmax  ,nrmax+1,nsamax))
    allocate(Drtfow(nthmax  ,npmax  ,nrmax+1,nsamax))
    allocate(Drrfow(nthmax  ,npmax  ,nrmax+1,nsamax))
    allocate(Fppfow(nthmax  ,npmax+1,nrmax  ,nsamax))
    allocate(Fthfow(nthmax+1,npmax  ,nrmax  ,nsamax))
    allocate(Frrfow(nthmax  ,npmax  ,nrmax+1,nsamax))

  end subroutine fow_allocate

  subroutine fow_deallocate

    deallocate(JI)
    deallocate(JIR)
    deallocate(rhom_local)
    deallocate(thetam_local)
    deallocate(time_loss)
    deallocate(thetam)
    deallocate(thetamg)
    deallocate(thetam_pg)
    deallocate(thetam_rg)

    deallocate(psim)
    deallocate(Fpsi)
    deallocate(Bout)
    deallocate(Bin)
    deallocate(psimg)
    deallocate(Fpsig)
    deallocate(Boutg)
    deallocate(Bing)
    deallocate(dFdr)
    deallocate(dBoutdr)
    deallocate(dBindr)
    deallocate(dpsimdr)
    deallocate(dFgdr)
    deallocate(dBoutgdr)
    deallocate(dBingdr)
    deallocate(dpsimgdr)

    deallocate(Babs)
    deallocate(dBdr)
    deallocate(dBdthp)
    deallocate(theta_p)

    deallocate(theta_pnc)
    deallocate(theta_co_stg)
    deallocate(theta_cnt_stg)
    deallocate(psip_pnc_point)
    deallocate(theta_pnc_pg)
    deallocate(theta_co_stg_pg)
    deallocate(theta_cnt_stg_pg)
    deallocate(psip_pnc_point_pg)
    deallocate(theta_pnc_rg)
    deallocate(theta_co_stg_rg)
    deallocate(theta_cnt_stg_rg)
    deallocate(psip_pnc_point_rg)

    deallocate(delthm)
    deallocate(delthm_pg)
    deallocate(delthm_rg)
    deallocate(IBCflux_ratio)
    deallocate(nr_pnc_point)
    deallocate(nr_rhom_pinch)
    deallocate(orbital_loss)

    deallocate(nth_stg)
    deallocate(nth_pnc)

    deallocate(orbit_th)
    deallocate(orbit_p)
    deallocate(orbit_r)
    deallocate(orbit_m)

    deallocate(Dppfow)
    deallocate(Dptfow)
    deallocate(Dprfow)
    deallocate(Dtpfow)
    deallocate(Dttfow)
    deallocate(Dtrfow)
    deallocate(Drpfow)
    deallocate(Drtfow)
    deallocate(Drrfow)
    deallocate(Fppfow)
    deallocate(Fthfow)
    deallocate(Frrfow)

  end subroutine fow_deallocate

end module fowcomm
