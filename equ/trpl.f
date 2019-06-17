c
      module trpl_mod
      use bpsd
      type(bpsd_species_type),private,save :: species
      type(bpsd_plasmaf_type),private,save :: plasmaf
      logical, private, save :: trpl_init_flag = .TRUE.
      public
      contains
c
c=======================================================================
      subroutine trpl_init
c=======================================================================
      use trn_mod
      implicit none
! local variables
      integer    ns,nr,ierr
c=======================================================================
      if(trpl_init_flag) then
         species%nsmax=0
         plasmaf%nsmax=0
         plasmaf%nrmax=0
         trpl_init_flag=.FALSE.
      endif
c
      if(species%nsmax.ne.mion+1) then
         if(associated(species%data)) then
            deallocate(species%data)
         endif
         species%nsmax=mion+1
         allocate(species%data(species%nsmax))
      endif
c
      do ns=1,species%nsmax
         species%data(ns)%pa=fmass(ns-1)
         species%data(ns)%pz=fchrg(ns-1)
      enddo
      call bpsd_set_species(species,ierr)
c
      if((plasmaf%nsmax.ne.mion+1).or.
     &   (plasmaf%nrmax.ne.nro)) then
         if(associated(plasmaf%rho)) then
            deallocate(plasmaf%rho)
         endif
         if(associated(plasmaf%data)) then
            deallocate(plasmaf%data)
         endif
         if(associated(plasmaf%qinv)) then
            deallocate(plasmaf%qinv)
         endif
         plasmaf%nsmax=mion+1
         plasmaf%nrmax=nro
         allocate(plasmaf%rho(plasmaf%nrmax))
         allocate(plasmaf%data(plasmaf%nrmax,plasmaf%nsmax))
         allocate(plasmaf%qinv(plasmaf%nrmax))
      endif
c
      plasmaf%time=0.d0
      do nr=1,plasmaf%nrmax
         plasmaf%rho(nr)=ro(nr)
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=0.d0
            plasmaf%data(nr,ns)%pt=0.d0
            plasmaf%data(nr,ns)%ptpr=0.d0
            plasmaf%data(nr,ns)%ptpp=0.d0
            plasmaf%data(nr,ns)%pu=0.d0
         enddo
         plasmaf%qinv(nr)=0.D0
      enddo
      call bpsd_set_plasmaf(plasmaf,ierr)
      return
      end subroutine trpl_init
c
c=======================================================================
      subroutine trpl_set(ierr)
c=======================================================================
      use trn_mod
      implicit none
      integer    ierr
! local variables
      integer    ns,nr
c=======================================================================
c
      plasmaf%time=0.d0
      do nr=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=den(nr,ns-1)
            plasmaf%data(nr,ns)%pt=tem(nr,ns-1)
            plasmaf%data(nr,ns)%ptpr=tem(nr,ns-1)
            plasmaf%data(nr,ns)%ptpp=tem(nr,ns-1)
            plasmaf%data(nr,ns)%pu=0.d0
         enddo
         plasmaf%qinv(nr)=qi(nr)*(2.D0*cnpi)**2
      enddo
      call bpsd_set_plasmaf(plasmaf,ierr)
      return
      end subroutine trpl_set
c
c=======================================================================
      subroutine trpl_get(ierr)
c=======================================================================
      use trn_mod
      implicit none
      integer    ierr
! local variables
      integer    ns,nr
c=======================================================================
c
      call bpsd_get_plasmaf(plasmaf,ierr)
c
      do nr=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            den(nr,ns-1)=plasmaf%data(nr,ns)%pn
            tem(nr,ns-1)=plasmaf%data(nr,ns)%pt
            pre(nr,ns-1)=cnec*tem(nr,ns-1)*den(nr,ns-1)
         enddo
         qi(nr)=plasmaf%qinv(nr)/(2.d0*cnpi)**2
      enddo
      return
      end subroutine trpl_get
c
      end module trpl_mod
