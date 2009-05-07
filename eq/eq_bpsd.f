      module eq_bpsd_mod

      use bpsd
      public eq_bpsd_init, eq_bpsd_set, eq_bpsd_get
      private

      type(bpsd_device_type),  private,save :: device
      type(bpsd_equ1D_type),   private,save :: equ1D
      type(bpsd_metric1D_type),private,save :: metric1D
      type(bpsd_species_type), private,save :: species
      type(bpsd_plasmaf_type), private,save :: plasmaf
      logical, private, save :: eq_bpsd_init_flag = .TRUE.
      contains
!=======================================================================
      subroutine eq_bpsd_init(ierr)
!=======================================================================
      INCLUDE '../eq/eqcomq.inc'
!      implicit none
      integer(4) :: ierr
! local variables
      integer(4) :: ns
      real(8) :: pretot, dentot, temave
!=======================================================================
      if(eq_bpsd_init_flag) then
         equ1D%nrmax=0
         metric1D%nrmax=0
         eq_bpsd_init_flag=.FALSE.
      endif

      device%rr   = RR
      device%zz   = 0.d0
      device%ra   = RA
      device%rb   = RB
      device%bb   = BB
      device%ip   = RIP
      device%elip = RKAP
      device%trig = RDLT
      call bpsd_set_data(device,ierr)

      equ1D%time=0.D0
      if(equ1D%nrmax.ne.nrmax) then
         if(associated(equ1D%rho)) deallocate(equ1D%rho)
         if(associated(equ1D%data)) deallocate(equ1D%data)
         equ1D%nrmax=NRMAX
         allocate(equ1D%rho(NRMAX))
         allocate(equ1D%data(NRMAX))
      endif

      metric1D%time=0.D0
      if(metric1D%nrmax.ne.nrmax) then
         if(associated(metric1D%rho)) deallocate(metric1d%rho)
         if(associated(metric1D%data)) deallocate(metric1d%data)
         metric1D%nrmax=NRMAX
         allocate(metric1D%rho(NRMAX))
         allocate(metric1D%data(NRMAX))
      endif

      return
      end subroutine eq_bpsd_init

!=======================================================================
      subroutine eq_bpsd_set(ierr)
!=======================================================================
!     interface eqiulibrium => transport
!         equilibrium grid => transport grid     
!=======================================================================
      INCLUDE '../eq/eqcomq.inc'
      integer(4) :: nr, ierr
      real(8),dimension(nrmax)  :: sa, data, diff, tmp
      real(8),dimension(4,nrmax):: udata
      
! local variables
!=======================================================================

      plasmaf%nrmax=0
      call bpsd_get_data(plasmaf,ierr)
      do nr=1,nrmax
         sa(nr)=rhot(nr)**2 ! = psit(nr)/psit(nrmax)
         data(nr)=1.d0/qps(nr)
      enddo
      call spl1d(sa,data,diff,udata,nrmax,0,ierr)
      if(ierr.ne.0) write(6,*) 'eq_bpsd_set: spl1d: ierr=',ierr
      do nr=1,plasmaf%nrmax
         call spl1df((plasmaf%rho(nr))**2,plasmaf%qinv(nr),
     &               sa,udata,nrmax,ierr)
         if(ierr.ne.0) write(6,*) 'eq_bpsd_set: spl1df: ierr=',ierr
      enddo
      call bpsd_set_data(plasmaf,ierr)

      do nr=1,nrmax
         equ1D%rho(nr)       = rhot(nr) ! = SQRT(psit(nr)/psit(nrmax))
         equ1D%data(nr)%psit = psit(nr)
         equ1D%data(nr)%psip = psip(nr)
         equ1D%data(nr)%ppp  = pps(nr)
         equ1D%data(nr)%piq  = 1.d0/qps(nr)
         equ1D%data(nr)%pip  = tts(nr)/rmu0
         equ1D%data(nr)%pit  = 0.d0
      enddo

      call bpsd_set_equ1D(equ1D,ierr)

      do nr=1,nrmax
         metric1D%rho(nr)          = rhot(nr)
         metric1D%data(nr)%pvol    = fnvps(rhot(nr))    ! vps     on rhot
         metric1D%data(nr)%psur    = fnsps(rhot(nr))    ! sps     on rhot
         metric1D%data(nr)%dvpsit  = fndvdpst(rhot(nr)) ! dvdpsit on rhot
         metric1D%data(nr)%dvpsip  = fndvdpsp(rhot(nr)) ! dvdpsip on rhot
         metric1D%data(nr)%aver2   = fnavrr2 (rhot(nr)) ! averr2  on rhot
         metric1D%data(nr)%aver2i  = fnavir2 (rhot(nr)) ! aveir2  on rhot
         metric1D%data(nr)%aveb2   = fnavbb2 (rhot(nr)) ! avebb2  on rhot
         metric1D%data(nr)%aveb2i  = fnavib2 (rhot(nr)) ! aveib2  on rhot
         metric1D%data(nr)%avegv   = fnavgv  (rhot(nr)) ! avegv   on rhot
         metric1D%data(nr)%avegv2  = fnavgv2 (rhot(nr)) ! avegv2  on rhot
         metric1D%data(nr)%avegvr2 = fnavgvr2(rhot(nr)) ! avegvr2 on rhot
         metric1D%data(nr)%avegr   = fnavgr  (rhot(nr)) ! avegr   on rhot
         metric1D%data(nr)%avegr2  = fnavgr2 (rhot(nr)) ! avegr2  on rhot
         metric1D%data(nr)%avegrr2 = fnavgrr2(rhot(nr)) ! avegrr2 on rhot
         metric1D%data(nr)%avegpp2 = fnavgp2 (rhot(nr)) ! avegp2  on rhot
         metric1D%data(nr)%rr      = fnrrps  (rhot(nr)) ! rrpsi   on rhot
         metric1D%data(nr)%rs      = fnrsps  (rhot(nr)) ! rspsi   on rhot
         metric1D%data(nr)%elip    = fnelpps (rhot(nr)) ! elippsi on rhot
         metric1D%data(nr)%trig    = fntrgps (rhot(nr)) ! trigpsi on rhot
         metric1D%data(nr)%aveb    = fnavbb  (rhot(nr)) ! avebb2  on rhot
      enddo
      call bpsd_set_metric1D(metric1D,ierr)

      end subroutine eq_bpsd_set

!=======================================================================
      subroutine eq_bpsd_get(ierr)
!=======================================================================
!     interface transport => equilibrium
!=======================================================================
      INCLUDE '../eq/eqcomq.inc'
      integer(4) :: ierr

      call bpsd_get_data(device,ierr)

      RR=device%rr
      RA=device%ra
      RB=device%rb
      BB=device%bb
      RIP=device%ip
      RKAP=device%elip
      RDLT=device%trig

      equ1D%time=0.D0
      if(equ1D%nrmax.ne.nrmax) then
         if(associated(equ1D%rho)) deallocate(equ1D%rho)
         if(associated(equ1D%data)) deallocate(equ1D%data)
         equ1D%nrmax=NRMAX
         allocate(equ1D%rho(NRMAX))
         allocate(equ1D%data(NRMAX))
      endif

      metric1D%time=0.D0
      if(metric1D%nrmax.ne.nrmax) then
         if(associated(metric1D%rho)) deallocate(metric1d%rho)
         if(associated(metric1D%data)) deallocate(metric1d%data)
         metric1D%nrmax=NRMAX
         allocate(metric1D%rho(NRMAX))
         allocate(metric1D%data(NRMAX))
      endif
      end subroutine eq_bpsd_get

      end module eq_bpsd_mod
