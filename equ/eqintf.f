c
      module eqintf_mod
      use eqsub_mod
      use tpxssl_mod
      public
      contains
c
c=======================================================================
      subroutine intequ
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
      implicit none
c     include 'BAL'
c     include 'NBI'
c     include 'FUS'
! local variables
      integer   i, icend, in, ir, is, it, iturn, iturnm, iz
     >        , izdz, izen, izst, j, jz, n
      real*8    gzte(itdm,0:indmz), muta(itdm),  mutb(itdm)
     >        , mute(itdm,0:indmz), s
     >        , wkt(itdm), wkv(ivdm)
      common/intf/gzte,mute
c=======================================================================
c     transport quantities      trn-grid --> equ-grid    :  hiv(n)
c=======================================================================
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
      end subroutine intequ
c
c=======================================================================
      subroutine inttrn
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
      implicit none
c     include 'BAL'
c     include 'NBI'
c     include 'FUS'
! local variables
      integer   i, icend, in, ir, ire, irs, is, it, iturn, iturnm, iz
     >        , izdz, izen, izst, j, jz, n, izs, ize
      real*8    gzte(itdm,0:indmz), muta(itdm),  mutb(itdm)
     >        , mute(itdm,0:indmz), s
     >        , wkt(itdm), wkv(ivdm)
      common/intf/gzte,mute
c=======================================================================
c     equilibrium quantities      equ-grid --> trn-grid  : ro(n)
c-----------------------------------------------------------------------
      call spln(sit,hit,nt,siv,hiv,nv,ww1,0)
      call spln(sdt,hit,nt,sdv,hiv,nv,ww1,0)
      call spln(hdt,hit,nt,hdv,hiv,nv,ww1,0)
      call spln(vlt,hit,nt,vlv,hiv,nv,ww1,0)
      call spln(art,hit,nt,arv,hiv,nv,ww1,0)
      call spln(ckt,hit,nt,ckv,hiv,nv,ww1,0)
      call spln(sst,hit,nt,ssv,hiv,nv,ww1,0)
      call spln(aat,hit,nt,aav,hiv,nv,ww1,0)
      call spln(rrt,hit,nt,rrv,hiv,nv,ww1,0)
c-----------------------------------------------------------------------
c     equilibrium quantities      equ-grid --> trn-grid  : roh(n)
c-----------------------------------------------------------------------
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
      end subroutine inttrn
c
      end module eqintf_mod
