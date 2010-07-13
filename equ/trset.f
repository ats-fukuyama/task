c
      module trset_mod
      use tpxssl_mod
      public
      contains
c
c=======================================================================
      subroutine trset
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use eqt_mod
      use trn_mod
      use cnt_mod
      use com_mod
      implicit none
!     include 'TRSET'
!     include 'IMP'
!     include 'PUF'
!     include 'SOR'
! local variables
      integer    is,n
c=======================================================================
c  set up basic parameters of transport calculation
c
c  define transport grid  : ro(n) n=1,nro
c                main plasma : n=1,nroblk
c                scl  plasma : n=noblk+1,nro(=nroblk+nroscl)
c
      nroblk=nv
      nroscl=0
c
      nt=nroblk
      ntm=nt-1
c
      nro=nroblk+nroscl
      nrom=nro-1
      nroblm=nroblk-1
      do n=1,nroblk
      ro(n)=dfloat(n-1)/dfloat(nroblm)
      hit(n)=hiv(nv)*ro(n)**2
      enddo
c
!     if(nroscl.gt.0)then
!     do n=nroblk+1,nro
!     ro(n)=1.+roscl*dfloat(n-nroblk)/dfloat(nroscl)
!     hit(n)=hiv(nv)*ro(n)**2
!     enddo
!     endif
c-----
      do n=1,nrom
      roh(n)=0.5*(ro(n)+ro(n+1))
      hih(n)=hiv(nv)*roh(n)**2
      enddo
c-----
C
      call spln(sit,hit,nt,siv,hiv,nv,ww1,0)
      call spln(mut,hit,nt,muv,hiv,nv,ww1,0)
      call spln(nut,hit,nt,nuv,hiv,nv,ww1,0)
      call spln(hdt,hit,nt,hdv,hiv,nv,ww1,0)
      call spln(vlt,hit,nt,vlv,hiv,nv,ww1,0)
c-----
!     if(nroscl.gt.0)then
!     do n=nroblk+1,nro
!     sit(n)=sit(nv)*ro(n)**2
!     hdt(n)=hdv(nv)*ro(n)**2
!     nut(n)=nuv(nv)*ro(n)**2
!     vlt(n)=vlv(nv)*ro(n)**2
!     enddo
!     endif
c-----------------------------------------------------------------------
c set pressure profiles
      mion=1
      do n=1,nro
      do is=0,mion
      pre(n,is)=(mut(n)/cnmu)*hdt(n)**gam/dfloat(mion+1)
      enddo
      enddo
c-----------------------------------------------------------------------
c set density profiles
      do n=1,nro
      den(n,0)=1.d+20*(1.d0-ro(n)**4)+1.d+18
      do is=1,mion
      den(n,is)=den(n,0)
      enddo
      enddo
c-----------------------------------------------------------------------
c set temperature profiles
      do n=1,nro
      do is=0,mion
      tem(n,is)=pre(n,is)/(cnec*den(n,is))
      enddo
      enddo
c-----------------------------------------------------------------------
c set pressure profiles
      do n=1,nro
      do is=0,mion
      pre(n,is)=cnec*tem(n,is)*den(n,is)
      enddo
      enddo
c-----------------------------------------------------------------------
c set safety profiles
      do n=1,nro
      qi(n)=nut(n)
      enddo
      qi(1)=(4.*qi(2)-qi(3))/3.
c-----------------------------------------------------------------------
c set auxualy quantities
      do n=1,nro
      do is=0,mion
      denc(n,is)=den(n,is)
      temc(n,is)=tem(n,is)
      prec(n,is)=pre(n,is)
      den0(n,is)=den(n,is)
      pre0(n,is)=pre(n,is)
      enddo
c-----
      qic(n)=qi(n)
      qi0(n)=qi(n)
      enddo
      qic(1)=qi(1)
      qi0(1)=qi(1)
c-----------------------------------------------------------------------
      return
      end subroutine trset
c
      end module trset_mod
