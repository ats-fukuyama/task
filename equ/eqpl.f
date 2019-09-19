c
      module equpl_mod
c
      use bpsd
      use tpxssl_mod
      public eqpl_init, eqpl_prof, eqpl_get, eqpl_put
      private
c
      type(bpsd_device_type),private,save  :: device_eq
      type(bpsd_equ1D_type),private,save   :: equ1D
      type(bpsd_metric1D_type),private,save:: metric1D
      type(bpsd_species_type),private,save :: species
      type(bpsd_plasmaf_type),private,save :: plasmaf
      logical, private, save :: eqpl_init_flag = .TRUE.
      contains
c=======================================================================
      subroutine eqpl_init(ierr)
c=======================================================================
      use aaa_mod
      use com_mod
      use par_mod
      use equ_mod
      use eqv_mod
      use trn_mod
      use eqt_mod
      use vac_mod
      implicit none
      integer ierr
! local variables
      integer(4)    ns
      real(8) pretot,dentot,temave
c=======================================================================
      if(eqpl_init_flag) then
         equ1D%nrmax=0
         metric1D%nrmax=0
         eqpl_init_flag=.FALSE.
         bpsd_debug_flag=.FALSE.
      endif
c
      device_eq%rr=rmaj
      device_eq%zz=zpla
      device_eq%ra=rpla
      device_eq%rb=rpla*1.2d0
      device_eq%bb=btv/rmaj
      device_eq%ip=tcu
      device_eq%elip=elip
      device_eq%trig=trig
      call bpsd_put_data(device_eq,ierr)
c
      equ1D%time=0.D0
      if(equ1D%nrmax.ne.nv) then
         if(associated(equ1D%rho)) deallocate(equ1D%rho)
         if(associated(equ1D%data)) deallocate(equ1D%data)
         equ1D%nrmax=nv
         allocate(equ1D%rho(nv))
         allocate(equ1D%data(nv))
      endif
c
      metric1D%time=0.D0
      if(metric1D%nrmax.ne.nv) then
         if(associated(metric1D%rho)) deallocate(metric1d%rho)
         if(associated(metric1D%data)) deallocate(metric1d%data)
         metric1D%nrmax=nv
         allocate(metric1D%rho(nv))
         allocate(metric1D%data(nv))
      endif

      return
      end subroutine eqpl_init
c
c=======================================================================
      subroutine eqpl_prof(ierr)
c=======================================================================
      use aaa_mod
      use com_mod
      use par_mod
      use equ_mod
      use eqv_mod
      use trn_mod
      use eqt_mod
      use vac_mod
      implicit none
      integer ierr
! local variables
      integer(4)    ns,n
      real(8) pretot,dentot,temave
c=======================================================================

!----- adjust pressure profile -----

      plasmaf%nrmax=0
      call bpsd_get_plasmaf(plasmaf,ierr)
      nt=plasmaf%nrmax
      ntm=nt-1
      do n=1,nt
         ro(n)=plasmaf%rho(n)
         hit(n)=hiv(nv)*plasmaf%rho(n)**2
      enddo
      do n=1,ntm
         roh(n)=0.5*(ro(n)+ro(n+1))
         hih(n)=hiv(nv)*roh(n)**2
      enddo

      call spln(sit,hit,nt,siv,hiv,nv,ww1,0)
      call spln(mut,hit,nt,muv,hiv,nv,ww1,0)
      call spln(nut,hit,nt,nuv,hiv,nv,ww1,0)
      call spln(hdt,hit,nt,hdv,hiv,nv,ww1,0)
      call spln(vlt,hit,nt,vlv,hiv,nv,ww1,0)

c set pressure profiles
      do n=1,nt
         pretot=(mut(n)/cnmu)*hdt(n)**gam
         dentot=0.d0
         do ns=1,plasmaf%nsmax
            dentot=dentot+plasmaf%data(n,ns)%density
         enddo
         temave=pretot/(cnec*dentot)
         do ns=1,plasmaf%nsmax
            plasmaf%data(n,ns)%temperature=temave
            plasmaf%data(n,ns)%temperature_para=temave
            plasmaf%data(n,ns)%temperature_perp=temave
         enddo
      enddo

      do n=1,nt
         plasmaf%qinv(n)=nut(n)*(2.D0*cnpi)**2
      enddo

      call bpsd_put_plasmaf(plasmaf,ierr)

      return
      end subroutine eqpl_prof
c
c=======================================================================
      subroutine eqpl_get(ierr)
c=======================================================================
c     interface eqiulibrium <>transport                            JAERI
c          transport grid -> equilibrium grid
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use eqt_mod
      use trn_mod
      use r2d_mod
      use com_mod
      use imp_mod
      use tpxssl_mod
      implicit none
      integer ierr
! local variables
      integer   i, icend, in, ir, is, it, iturn, iturnm, iz
     >        , izdz, izen, izst, j, jz, n, ns
      real*8    gzte(itdm,0:indmz), muta(itdm),  mutb(itdm)
     >        , mute(itdm,0:indmz), s
     >        , wkt(itdm), wkv(ivdm)
      common/intf/gzte,mute
c=======================================================================
c     transport quantities      trn-grid --> equ-grid    :  hiv(n)
c=======================================================================
      plasmaf%nrmax=0
      call bpsd_get_plasmaf(plasmaf,ierr)

      nt=plasmaf%nrmax
      mion=plasmaf%nsmax-1
      ntm=nt-1
      nro=nv
      do n=1,nt
         ro(n)=plasmaf%rho(n)
         hit(n)=hiv(nv)*plasmaf%rho(n)**2
         qi(n)=plasmaf%qinv(n)/(2.D0*cnpi)**2
         do ns=1,plasmaf%nsmax
            den(n,ns-1)=plasmaf%data(n,ns)%density
            tem(n,ns-1)=plasmaf%data(n,ns)%temperature
            pre(n,ns-1)=cnec*tem(n,ns-1)*den(n,ns-1)
         enddo
      enddo
c
      sit(nt)=0.
      do n=ntm,1,-1
      sit(n)=sit(n+1)-0.5*(qi(n+1)+qi(n))*(hit(n+1)-hit(n))
      enddo
      do n=1,nt
      mut(n)=0.
      do is=0,mion
      gzte(n,is)=den(n,is)/hdt(n)
      mute(n,is)=pre(n,is)/hdt(n)**gam
      mut(n)=mut(n)+mute(n,is)
      enddo
!     if(n.lt.nt)then
!     mutb(n)=prnbi(n)/hdt(n)**gam
!     muta(n)=pralf(n)/hdt(n)**gam
!     else
!     mutb(n)=0.
!     muta(n)=0.
!     endif
!     mut(n)=mut(n)+mutb(n)+muta(n)
      mut(n)=cnmu*mut(n)
      nut(n)=qi(n)
      write(27,'(I5,1p4e12.4)') n,ro(n),mut(n),nut(n),sit(n)
      enddo
c-----------------------------------------------------------------------
      do n=1,nv
      siv(n)=sit(1)+(sit(nt)-sit(1))*float(n-1)/float(nvm)
      enddo
      nnhit=nt
      do n=1,nt
      xxhit(n)=mut(n)
      yyhit(n)=nut(n)
      zzhit(n)=hit(n)
      enddo
c-----------------------------------------------------------------------
      call spln(muv,hiv,nv,mut,hit,nt,wwmu,0)
      call spln(nuv,hiv,nv,nut,hit,nt,wwmu,0)
c-----------------------------------------------------------------------
      return
      end subroutine eqpl_get
c
c=======================================================================
      subroutine eqpl_put(ierr)
c=======================================================================
c     interface eqiulibrium <>transport                            JAERI
c         transport grid -> equilibrium grid     
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use eqt_mod
      use trn_mod
      use r2d_mod
      use com_mod
      use imp_mod
      use tpxssl_mod
      implicit none
      integer ierr
! local variables
      integer   i, icend, in, ir, ire, irs, is, it, iturn, iturnm, iz
     >        , izdz, izen, izst, j, jz, n, izs, ize, ns
      real*8    gzte(itdm,0:indmz), muta(itdm),  mutb(itdm)
     >        , mute(itdm,0:indmz), s
     >        , wkt(itdm), wkv(ivdm)
      common/intf/gzte,mute
c=======================================================================
c     equilibrium quantities      equ-grid --> trn-grid  : ro(n)
c-----------------------------------------------------------------------
c
      call spln(hdt,hit,nt,hdv,hiv,nv,ww1,0)
c
      do n=1,nv
         equ1D%rho(n)=SQRT(hiv(n)/hiv(nv))
         equ1D%data(n)%psit=hiv(n)/(2.D0*cnpi)
         equ1D%data(n)%psip=siv(n)*(2.D0*cnpi)
         equ1D%data(n)%ppp=muv(n)/cnmu*hdv(n)**gam
         equ1D%data(n)%piq=nuv(n)*(2.D0*cnpi)**2
         equ1D%data(n)%pip=hdv(n)/aav(n)/cnmu*(2.D0*cnpi)
         equ1D%data(n)%pit=ckv(n)*sdv(n)/cnmu/(2.D0*cnpi)
      enddo
      call bpsd_put_equ1D(equ1D,ierr)
c
      do n=1,nv
         metric1D%rho(n)=SQRT(hiv(n)/hiv(nv))
         metric1D%data(n)%pvol=vlv(n)
         metric1D%data(n)%psur=arv(n)
         metric1D%data(n)%dvpsit=1.d0/hdv(n)*(2.D0*cnpi)
         metric1D%data(n)%dvpsip=1.d0/sdv(n)/(2.D0*cnpi)
         metric1D%data(n)%aver2=rrv(n)
         metric1D%data(n)%aver2i=aav(n)
         metric1D%data(n)%aveb2=bbv(n)
         metric1D%data(n)%aveb2i=biv(n)
         metric1D%data(n)%avegv2=ssv(n)
         metric1D%data(n)%avegvr2=ckv(n)
         metric1D%data(n)%avegpp2=r2b2(n)*(2.D0*cnpi)**2
         metric1D%data(n)%rr=rtv(n)
         metric1D%data(n)%rs=rpv(n)
         metric1D%data(n)%elip=elv(n)
         metric1D%data(n)%trig=dlv(n)
      enddo
      call bpsd_put_metric1D(metric1D,ierr)
c
         
      call spln(sit,hit,nt,siv,hiv,nv,ww1,0)
      call spln(sdt,hit,nt,sdv,hiv,nv,ww1,0)
      call spln(vlt,hit,nt,vlv,hiv,nv,ww1,0)
      call spln(art,hit,nt,arv,hiv,nv,ww1,0)
      call spln(ckt,hit,nt,ckv,hiv,nv,ww1,0)
      call spln(sst,hit,nt,ssv,hiv,nv,ww1,0)
      call spln(aat,hit,nt,aav,hiv,nv,ww1,0)
      call spln(rrt,hit,nt,rrv,hiv,nv,ww1,0)
c-----------------------------------------------------------------------
c     equilibrium quantities      equ-grid --> trn-grid  : roh(n)
c-----------------------------------------------------------------------
      do n=1,ntm
      roh(n)=0.5*(ro(n)+ro(n+1))
      hih(n)=hiv(nv)*roh(n)**2
      enddo

      call spln(hdh,hih,ntm,hdv,hiv,nv,ww1,0)
      call spln(ckh,hih,ntm,ckv,hiv,nv,ww1,0)
      call spln(aah,hih,ntm,aav,hiv,nv,ww1,0)
      call spln(ssh,hih,ntm,ssv,hiv,nv,ww1,0)
      call spln(vlh,hih,ntm,vlv,hiv,nv,ww1,0)
      call spln(r2b2h,hih,ntm,r2b2,hiv,nv,ww1,0)
c-----
      do i=1,nv
      wkv(i)=rpv(i)**2
      enddo
      call spln(wkt,hih,ntm,wkv,hiv,nv,ww1,0)
      do i=1,ntm
      rph(i)=dsqrt(wkt(i))
      enddo
c-----
      call spln(rth,hih,ntm,rtv,hiv,nv,ww1,0)
      call spln(bbh,hih,ntm,bbv,hiv,nv,ww1,0)
      call spln(bih,hih,ntm,biv,hiv,nv,ww1,0)
c<<TRAPPED PARTICLE FRACTION>>
      do i=1,nv
      wkv(i)=ftr(i)**4.
      enddo
      call spln(wkt,hih,ntm,wkv,hiv,nv,ww1,0)
      do i=1,ntm
      fth(i)=wkt(i)**(1./4.)
      enddo
      do i=1,nv
      wkv(i)=epv(i)**2
      enddo
      call spln(wkt,hih,ntm,wkv,hiv,nv,ww1,0)
      do i=1,ntm
      eph(i)=wkt(i)**(1./2.)
      enddo
c-----------------------------------------------------------------------
c     CALL INTODE
c-----------------------------------------------------------------------
      do n=nt,1,-1
      m=n+1
      if(n.eq.nt)sit(n)=0.
      if(n.lt.nt)sit(n)=sit(m)-0.5*(nut(m)+nut(n))*(hit(m)-hit(n))
      rov(n)=dsqrt(vlt(n)/(2.*cnpi**2*rmaj))
      vro(n)=2.*hit(nt)*ro(n)/hdt(n)
      sdt(n)=nut(n)*hdt(n)
      if(n.eq.1)sro(1)=2.*sro(2)-sro(3)
      if(n.gt.1)sro(n)=-sst(n)*(hdt(n)/(2.*hit(nt)*ro(n)))**2
      if(n.lt.nt)then
      vrh(n)=2.*hit(nt)*roh(n)/hdh(n)
      srh(n)=ssh(n)*(hdh(n)/(2.*hit(nt)*roh(n)))**2
c--
      rovh(n)=dsqrt(vlh(n)/(2.*cnpi**2*rmaj))
      endif
      enddo
c-----
      do is=0,mion
      do n=1,nroblk
      den(n,is)=gzte(n,is)*hdt(n)
      pre(n,is)=mute(n,is)*hdt(n)**gam
      den0(n,is)=den(n,is)
      pre0(n,is)=pre(n,is)
      enddo
      enddo

      do n=1,nroblk
         qi(n)=nut(n)
         qi0(n)=qi(n)
      enddo


      do n=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            tem(n,ns-1)=pre(n,ns-1)/(cnec*den(n,ns-1))
            plasmaf%data(n,ns)%density=den(n,ns-1)
            plasmaf%data(n,ns)%temperature=tem(n,ns-1)
            plasmaf%data(n,ns)%temperature_para=tem(n,ns-1)
            plasmaf%data(n,ns)%temperature_perp=tem(n,ns-1)
            plasmaf%data(n,ns)%velocity_tor=0.d0
         enddo
c         plasmaf%qinv(n)=nut(n)*(2.D0*cnpi)**2
      enddo

      call bpsd_put_plasmaf(plasmaf,ierr)
c-----
c     call spln(cpdt,hit,nt,cpds,hiv,nv,ww1,0)
c=======================================================================
c     2d geometry for transport calculation   :  ro(r,z)
c=======================================================================
      ntr2d=nro
      ntr2db=nroblm
      ntr2dm=nro-1
c-----
      nrr2d=nr
      nzr2d=nz
      nrzr2d=nrz
      drr2d=dr
      dzr2d=dz
      rsr2d=rg(1)
      zsr2d=zg(1)
c=======================================================================
      nsur2d=nsu
      nscl2d=nsu
      do n=1,nsu
      rsur2d(n)=rsu(n)
      zsur2d(n)=zsu(n)
      rscl2d(n)=rsu(n)
      zscl2d(n)=zsu(n)
      enddo
c-----------------------------------------------------------------------
      izs=(zzsmin-zg(1))/dz+1.d0
      ize=(zzsmax-zg(1))/dz+1.d0
      irs=(zrsmin-rg(1))/dr+1.d0
      ire=(zrsmax-rg(1))/dr+1.d0
      it=nro
      do iz=1,nz
      do ir=1,nr
      i=(iz-1)*nr+ir
      rho2d(i)=dfloat(nro+1)
      if(iz.ge.izs.and.iz.le.ize)then
      if(ir.ge.irs.and.ir.le.ire)then
      if(psi(i).le.sit(nro))then
      s=psi(i)*sit(1)/saxis
      do while(s.le.sit(it).and.it.gt.1)
      it=it-1
      enddo
      do while(s.ge.sit(it+1).and.it.lt.nro-1)
      it=it+1
      enddo 
      rho2d(i)=it+(s-sit(it))/(sit(it+1)-sit(it))
      endif
      endif
      endif
      enddo
      enddo
c=======================================================================
      entry int2dp
c=======================================================================
c     for the 2d-calculation of transport process
c=======================================================================
      ionr2d=nion
      jonr2d=mion
      do in=0,jonr2d
      sper2d(in)=specie(in)
      fmsr2d(in)=fmass(in)
      chgr2d(in)=fchrg(in)
      do n=1,ntr2dm
      der2d(n,in)=denh(n,in)
      per2d(n,in)=preh(n,in)
      ter2d(n,in)=temh(n,in)
      enddo
      enddo
      do n=1,ntr2dm
      zef2d(n)=zef(n)
      qir2d(n)=qih(n)
      vlr2d(n)=vlt(n)
      enddo
      vlr2d(ntr2d)=vlt(nro)
c-----------------------------------------------------------------------
      return
      end subroutine eqpl_put
c
      end module equpl_mod
