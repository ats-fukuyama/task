c
      module eqpl_mod
c
      use bpsd_mod
      public eqpl_init, eqpl_set
      private
c
      type(bpsd_device_type),private,save  :: device
      type(bpsd_equ1D_type),private,save   :: equ1D
      type(bpsd_metric1D_type),private,save:: metric1D
      type(bpsd_species_type),private,save :: species
      type(bpsd_plasmaf_type),private,save :: plasmaf
      logical, private, save :: eqpl_init_flag = .TRUE.
      contains
c=======================================================================
      subroutine eqpl_init(ierr)
c=======================================================================
      INCLUDE '../eq/eqcomq.inc'
!      implicit none
      integer ierr
! local variables
      integer(4)    ns
      real(8) pretot,dentot,temave
c=======================================================================
      if(eqpl_init_flag) then
         equ1D%nrmax=0
         metric1D%nrmax=0
         eqpl_init_flag=.FALSE.
      endif
c
      device%rr=RR
      device%zz=0.d0
      device%ra=RA
      device%rb=RB
      device%bb=BB
      device%ip=RIP
      device%elip=RKAP
      device%trig=RDLT
      call bpsd_set_data(device,ierr)
c
      equ1D%time=0.D0
      if(equ1D%nrmax.ne.nv) then
         if(allocated(equ1D%s)) deallocate(equ1D%s)
         if(allocated(equ1D%data)) deallocate(equ1D%data)
         equ1D%nrmax=NRMAX
         allocate(equ1D%s(NRMAX))
         allocate(equ1D%data(NRMAX))
      endif
c
      metric1D%time=0.D0
      if(metric1D%nrmax.ne.nv) then
         if(allocated(metric1D%s)) deallocate(metric1d%s)
         if(allocated(metric1D%data)) deallocate(metric1d%data)
         metric1D%nrmax=NRMAX
         allocate(metric1D%s(NRMAX))
         allocate(metric1D%data(NRMAX))
      endif

      return
      end subroutine eqpl_init
c
c=======================================================================
      subroutine eqpl_set(ierr)
c=======================================================================
c     interface eqiulibrium <>transport                            JAERI
c         transport grid -> equilibrium grid     
c=======================================================================
      INCLUDE '../eq/eqcomq.inc'
      integer nr,ierr
      real(8):: s
      real(8),dimension(nrmax):: sa,data,diff
      real(8),dimension(4,nrmax):: udata
      
! local variables
c=======================================================================

      plasmaf%nrmax=0
      call bpsd_get_data(plasmaf,ierr)
      do nr=1,nrmax
         sa(nr)=psit(nr)/psit(nrmax)
         data(nr)=1.d0/qps(nr)
      enddo
      call spl1d(sa,data,diff,udata,nrmax,0,ierr)
      if(ierr.ne.0) write(6,*) 'eqpl_set: spl1d: ierr=',ierr
      do nr=1,plasmaf%nrmax
         call spl1df(plasmaf%s(nr),plasmaf%qinv(nr),sa,udata,nrmax,ierr)
         if(ierr.ne.0) write(6,*) 'eqpl_set: spl1df: ierr=',ierr
      enddo
      call bpsd_set_data(plasmaf,ierr)

      do nr=1,nrmax
         equ1D%s(nr)=psit(nr)/psit(nrmax)
         equ1D%data(nr)%psit=psit(nr)
         equ1D%data(nr)%psip=psip(nr)
         equ1D%data(nr)%ppp=pps(nr)
         equ1D%data(nr)%piq=1.d0/qps(nr)
         equ1D%data(nr)%pip=tts(nr)/rmu0
         equ1D%data(nr)%pit=0.d0
      enddo

      call bpsd_set_equ1D(equ1D,ierr)
c
      do nr=1,nrmax
         s=psit(nr)/psit(nrmax)
         metric1D%s(nr)=s
         metric1D%data(nr)%pvol=vps(nr)
         metric1D%data(nr)%psur=sps(nr)
         metric1D%data(nr)%dvpsit=dvdpsit(nr)
         metric1D%data(nr)%dvpsip=dvdpsip(nr)
         metric1D%data(nr)%aver2= averr2(nr)
         metric1D%data(nr)%aver2i=aveir2(nr)
         metric1D%data(nr)%aveb2= avebb2(nr)
         metric1D%data(nr)%aveb2i=aveib2(nr)
         metric1D%data(nr)%avegv2=avegv2(nr)
         metric1D%data(nr)%avegvr2=avegr2(nr)
         metric1D%data(nr)%avegpp2=avegp2(nr)
         metric1D%data(nr)%rr=rrpsi(nr)
         metric1D%data(nr)%rs=rspsi(nr)
         metric1D%data(nr)%elip=elippsi(nr)
         metric1D%data(nr)%trig=trigpsi(nr)
      enddo
      call bpsd_set_metric1D(metric1D,ierr)
      end subroutine eqpl_set
c
      end module eqpl_mod
