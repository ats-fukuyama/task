      module eq_bpsd_mod
! interface module with TASK/TR new version
      use bpsd
      public eq_bpsd_init, eq_bpsd_put, eq_bpsd_get
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
!=======================================================================
      if(eq_bpsd_init_flag) then
         equ1D%nrmax=0
         metric1D%nrmax=0
         eq_bpsd_init_flag=.FALSE.
      endif

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
      subroutine eq_bpsd_put(ierr)
!=======================================================================
!     interface eqiulibrium => transport
!          equilibrium grid => transport grid     
!=======================================================================
      INCLUDE '../eq/eqcomq.inc'
      INCLUDE '../eq/eqcom4.inc'
      integer(4) :: nr, ierr
      
! local variables
!=======================================================================

      if(eq_bpsd_init_flag) CALL eq_bpsd_init(ierr)

      device%rr   = RR
      device%zz   = 0.d0
      device%ra   = RA
      device%rb   = RB
      device%bb   = BB
      device%ip   = RIP
      device%elip = RKAP
      device%trig = RDLT
      call bpsd_put_data(device,ierr)

      equ1D%nrmax = nrpmax
      do nr=1,nrmax
         equ1D%rho(nr)       = rhot(nr) ! = SQRT(psit(nr)/psit(nrmax))
         equ1D%data(nr)%psit = psit(nr)
         equ1D%data(nr)%psip = psip(nr)
         equ1D%data(nr)%ppp  = pps(nr)
         equ1D%data(nr)%piq  = 1.d0/qps(nr)
         equ1D%data(nr)%pip  = tts(nr)/rmu0
         equ1D%data(nr)%pit  = 0.d0
      enddo
      call bpsd_put_data(equ1D,ierr)

      metric1D%nrmax = nrpmax
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
         metric1D%data(nr)%aveb    = fnavbb  (rhot(nr)) ! avebb   on rhot
      enddo
      call bpsd_put_data(metric1D,ierr)

      end subroutine eq_bpsd_put

!=======================================================================
      subroutine eq_bpsd_get(ierr)
!=======================================================================
!     interface transport => equilibrium
!=======================================================================
      INCLUDE '../eq/eqcomm.inc'
      INCLUDE '../eq/eqcom4.inc'
      real(8),DIMENSION(NTRM):: ptrrho,qtrrho,deriv
      integer(4) :: ierr,ntr

      ! --- device data ---
      CALL bpsd_get_data(device,ierr)
      RR   = device%rr
      RA   = device%ra
      BB   = device%bb
      RIP  = device%ip
      RKAP = device%elip
      RDLT = device%trig

      ! --- plasmaf data ---
      call bpsd_get_data(plasmaf,ierr)

      ntrmax = plasmaf%nrmax
      DO ntr = 1, ntrmax
         rhotr(ntr)  = plasmaf%rho(ntr)
         psitrx(ntr) = rhotr(ntr)**2
         ptot = 0.d0
         DO ns = 1, plasmaf%nsmax
            ptot = ptot
     &           + plasmaf%data(ntr,ns)%density
     &           * plasmaf%data(ntr,ns)%temperature
     &           * aee
         END DO
         ptrrho(ntr) = ptot*1.D-6
         qtrrho(ntr) = 1.d0/plasmaf%qinv(ntr)
C         write(6,'(A,I5,1P4E12.4)') 'eq_bpsd_get:',ntr,rhotr(ntr),
C     &        psitrx(ntr),ptrrho(ntr),qtrrho(ntr)
      END DO

      deriv=0.d0
      CALL spl1d(psitrx,ptrrho,deriv,uppsi,ntrmax,1,ierr)
      IF(ierr.NE.0) 
     &     WRITE(6,*) 'XX eq_bpsd_get: spl1d ptrrho: ierr=',ierr

      deriv=0.d0
      CALL spl1d(psitrx,qtrrho,deriv,uqpsi,ntrmax,1,ierr)
      IF(ierr.NE.0) 
     &     WRITE(6,*) 'XX eq_bpsd_get: spl1d qtrrho: ierr=',ierr

      mdleqf=9

      END SUBROUTINE eq_bpsd_get

      END MODULE eq_bpsd_mod
