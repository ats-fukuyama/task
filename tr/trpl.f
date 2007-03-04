c
      module trpl_m
      use bpsd_m
      type(bpsd_species_type),private,save :: species
      type(bpsd_equ1D_type),private,save :: equ1D
      type(bpsd_metric1D_type),private,save :: metric1D
      type(bpsd_plasmaf_type),private,save :: plasmaf
      logical, private, save :: trpl_init_flag = .TRUE.
      public
      contains
c
c=======================================================================
      subroutine trpl_init
c=======================================================================
      include 'trcomm.inc'
! local variables
      integer    ns,nr,ierr
c=======================================================================
      if(trpl_init_flag) then
         species%nsmax=0
         equ1D%nrmax=0
         metric1D%nrmax=0
         plasmaf%nsmax=0
         plasmaf%nrmax=0
         trpl_init_flag=.FALSE.
      endif
c
      if(species%nsmax.ne.nsmax) then
         if(allocated(species%data)) then
            deallocate(species%data)
         endif
         species%nsmax=nsmax
         allocate(species%data(species%nsmax))
      endif
c
      do ns=1,species%nsmax
         species%data(ns)%pa=pa(ns)
         species%data(ns)%pz=pz(ns)
         species%data(ns)%pz0=pz(ns)
      enddo
      call bpsd_set_species('species',species,ierr)
c
      if((equ1D%nrmax.ne.nrmax+1)) then
         if(allocated(equ1D%s)) then
            deallocate(equ1D%s)
         endif
         if(allocated(equ1D%data)) then
            deallocate(equ1D%data)
         endif
         equ1D%nrmax=nrmax+1
         allocate(equ1D%s(equ1D%nrmax))
         allocate(equ1D%data(equ1D%nrmax)
      endif
c
      equ1D%time=0.d0
      do nr=1,equ1D%nrmax
         equ1D%s(nr)=rg(nr)**2
      enddo
      do nr=1,equ1D%nrmax
         equ1D%data(nr)%psit=0.d0
         equ1D%data(nr)%psip=0.d0
         equ1D%data(nr)%ppp=0.d0
         equ1D%data(nr)%piq=0.d0
         equ1D%data(nr)%pip=0.d0
         equ1D%data(nr)%pit=0.d0
      enddo
c
      if((metric1D%nrmax.ne.nrmax+1)) then
         if(allocated(metric1D%s)) then
            deallocate(metric1D%s)
         endif
         if(allocated(metric1D%data)) then
            deallocate(metric1D%data)
         endif
         metric1D%nrmax=nrmax+1
         allocate(metric1D%s(metric1D%nrmax))
         allocate(metric1D%data(metric1D%nrmax)
      endif
c
      metric1D%time=0.d0
      do nr=1,metric1D%nrmax
         metric1D%s(nr)=rg(nr)**2
      enddo
      do nr=1,metric1D%nrmax
         metric1D%data(nr)%pvol=0.d0
         metric1D%data(nr)%psur=0.d0
         metric1D%data(nr)%dvpsit=0.d0
         metric1D%data(nr)%dvpsip=0.d0
         metric1D%data(nr)%aver2=0.d0
         metric1D%data(nr)%aver2i=0.d0
         metric1D%data(nr)%aveb2=0.d0
         metric1D%data(nr)%aveb2i=0.d0
         metric1D%data(nr)%avegv2=0.d0
         metric1D%data(nr)%avegvr2=0.d0
         metric1D%data(nr)%avegpp2=0.d0
         metric1D%data(nr)%averr=0.d0
         metric1D%data(nr)%avera=0.d0
         metric1D%data(nr)%avekappa=0.d0
         metric1D%data(nr)%avedelta=0.d0
      enddo
c
      if((plasmaf%nsmax.ne.nsmax).or.
     &   (plasmaf%nrmax.ne.nrmax+1)) then
         if(allocated(plasmaf%s)) then
            deallocate(plasmaf%s)
         endif
         if(allocated(plasmaf%data)) then
            deallocate(plasmaf%data)
         endif
         if(allocated(plasmaf%qinv)) then
            deallocate(plasmaf%qinv)
         endif
         plasmaf%nsmax=nsmax
         plasmaf%nrmax=nrmax+1
         allocate(plasmaf%s(plasmaf%nrmax))
         allocate(plasmaf%data(plasmaf%nrmax,plasmaf%nsmax))
         allocate(plasmaf%qinv(plasmaf%nrmax))
      endif
c
      plasmaf%time=0.d0
      do nr=1,plasmaf%nrmax
         plasmaf%s(nr)=rg(nr)**2
      enddo
      do nr=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=0.d0
            plasmaf%data(nr,ns)%pt=0.d0
            plasmaf%data(nr,ns)%ptpr=0.d0
            plasmaf%data(nr,ns)%ptpp=0.d0
            plasmaf%data(nr,ns)%pu=0.d0
         enddo
         plasmaf%qinv(nr)=0.D0
      enddo
      call bpsd_set_plasmaf('plasmaf',plasmaf,ierr)
      return
      end subroutine trpl_init
c
c=======================================================================
      subroutine trpl_set(ierr)
c=======================================================================
      include 'trcomm.inc'
      integer    ierr
! local variables
      integer    ns,nr
      real*8 temp(nrmp,nsm,3)
c=======================================================================
c
      plasmaf%time=t
      do ns=1,nsmax
         call mesh_convert_mtog(rn(1,ns),temp(1,ns,1),nrmax)
         call mesh_convert_mtog(rt(1,ns),temp(1,ns,2),nrmax)
         call mesh_convert_mtog(ru(1,ns),temp(1,ns,3),nrmax)
      enddo
      do nr=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=temp(nr,ns,1)*1.d20
            plasmaf%data(nr,ns)%pt=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpr=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpp=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%pu=temp(nr,ns,3)
         enddo
      enddo
      do nr=2,plasmaf%nrmax
         plasmaf%qinv(nr)=1.d0/qp(nr-1)
      enddo
      plasmaf%qinv(1)=(4.d0*plasmaf%qinv(2)-plasmaf%qinv(2))/3.d0
c
      call bpsd_set_plasmaf('plasmaf',plasmaf,ierr)
      return
      end subroutine trpl_set
c
c=======================================================================
      subroutine trpl_get(ierr)
c=======================================================================
      include 'trcomm.inc'
      integer    ierr
! local variables
      integer    ns,nr
      real*8 temp(nrmp,nsm,3)
c=======================================================================
c
      call bpsd_get_plasmaf('plasmaf',plasmaf,ierr)
c
      do ns=1,plasmaf%nsmax
         do nr=1,plasmaf%nrmax
            temp(nr,ns,1)=plasmaf%data(nr,ns)%pn*1.d-20
            temp(nr,ns,2)=plasmaf%data(nr,ns)%pt*1.D-3
            temp(nr,ns,3)=plasmaf%data(nr,ns)%pu
         enddo
      enddo
      do nr=2,plasmaf%nrmax
         qp(nr-1)=1.d0/plasmaf%qinv(nr)
      enddo
c
      do ns=1,nsmax
         call mesh_convert_gtom(temp(1,ns,1),rn(1,ns),nrmax)
         call mesh_convert_gtom(temp(1,ns,2),rt(1,ns),nrmax)
         call mesh_convert_gtom(temp(1,ns,3),ru(1,ns),nrmax)
      enddo
c
      call bpsd_get_equ1D('equ1D',equ1D,ierr)
c
c
      call bpsd_get_metric1D('metric1D',metric1D,ierr)
c
      return
      end subroutine trpl_get
c
c     ----- convert half mesh to origin + grid mesh -----
c
      subroutine mesh_convert_mtog(datam,datag,nrmax)
c
      implicit none
      integer nrmax
      real*8 datam(nrmax),datag(nrmax+1)
      integer nr
c
      datag(1)=(9.d0*datam(1)-datam(2))/8.d0
      do nr=2,nrmax
         datag(nr)=(datam(nr-1)+datam(nr))/2.d0
      enddo
      datag(nrmax+1)=(4.d0*datam(nrmax)-datam(nrmax-1))/3.d0
      return
      end subroutine mesh_convert_mtog
c
c     ----- convert origin + grid mesh to half mesh -----
c
      subroutine mesh_convert_gtom(datag,datam,nrmax)
c
      implicit none
      integer nrmax
      real*8 datag(nrmax+1),datam(nrmax)
      real*8 c11,c12,c21,c22,det,a11,a12,a21,a22
      integer nr,ierr
c
      c11= 9.d0/8.d0
      c12=-1.d0/8.d0
      c21= 0.5d0
      c22= 0.5d0
      det=c11*c22-c12*c21
      a11= c22/det
      a12=-c12/det
      a21=-c21/det
      a22= c11/det
      datam(1)=a11*datag(1)+a12*datag(2)
      datam(2)=a21*datag(1)+a22*datag(2)
      do nr=3,nrmax
         datam(nr)=2.d0*datag(nr)-datam(nr-1)
      enddo
      return
      end subroutine mesh_convert_gtom
c
      end module trpl_m
