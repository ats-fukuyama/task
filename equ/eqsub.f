c
      module eqsub_mod
      use tpxssl_mod
      use eqpfds_mod
      public
      contains
c=======================================================================
      subroutine eqadj
c=======================================================================
c     ADJUST OF VACCUME FIELDS                                     JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use par_mod
      implicit none
! local variables
      integer   i, ill, j, k
      real*8    amat(icvdm2,icvdm2), bmat(icvdm2,icvdm2), err
     >        , rrm, rvnew, savee1, savee2, saver1, saver2
     >        , scale, x
      data savee1, saver1 / 0.0d0, 0.0d0 /
c=======================================================================
      scale=1.d0
      error=0.d0
      do i=1,icvdm1
       do j=1,icvdm2
        amat(i,j)=0.d0
       enddo
       amat(i,i)=1.d0
      enddo
c-----------------------------------------------------------------------
      rrm=(rvac(0)+rvac(1))/2.d0
      error=1.d0-bts/btv
      saver2=saver1
      saver1=rvac(0)
      savee2=savee1
      savee1=error
      rvnew=rvac(0)
      if(iad.eq.2)rvnew=rvac(0)+error*(rrm-rvac(0))
      if(iad.gt.2.and.dabs(savee2-error).gt.1.d-10)
     >            rvnew=rvac(0)-error*(saver2-rvac(0))/(savee2-error)
      x=(rvnew-saver1)/dr
      if(x.gt. 1.d-02)x= 1.d-02
      if(x.lt.-1.d-02)x=-1.d-02
      rvnew=saver1+x*dr
      scale=(rrm-rvnew)/(rrm-saver1)
c-----------------------------------------------------------------------
      do k=0,iabs( msfx )
       rvac(k)=rrm+scale*(rvac(k)-rrm)
       zvac(k)=scale*zvac(k)
      enddo
c-----------------------------------------------------------------------
c     PSI=PSI(P)+PSI(COIL OF I FIXED)
c----------------------------------------------------------------------
      do k=1,icvdm
       if(ivac(k).lt.0)then
        do i=1,nrz
         psi(i)=psi(i)+cvac(k)*svac(i,k)
        enddo
       endif
      enddo
c----------------------------------------------------------------------
      call eqadj0(amat)
c-----------------------------------------------------------------------
      call matslv(icvdm2,icvdm1,icvdm2,amat,bmat,err,ill)
c----------------------------------------------------------------------
c     PSI=PSI+PSI(COIL OF I FREE), SO PSI BECOME TOTAL PSI IN NEXT STEP
c-----------------------------------------------------------------------
      if(ivac(0).ge.0) then
       cvac(0)=amat(1,icvdm2)
       do i=1,nrz
        psi(i)=psi(i)+cvac(0)
       enddo
      endif
c-----
      do k=1,icvdm
       if(ivac(k).ge.0) then
        cvac(k)=amat(k+1,icvdm2)
        do i=1,nrz
         psi(i)=psi(i)+cvac(k)*svac(i,k)
        enddo
       endif
      enddo
c-----------------------------------------------------------------------
      if(iad.eq.1)return
      if(dabs(error).lt.eadmax)iad=-iad
      if(iad.le.iadmax)return
c-----------------------------------------------------------------------
      write(ft06,'(/5x,a29,1pd10.3)')
     >     '***** ERROR:EQADJ------ERROR=',error
      iad=-iad
      return
      end subroutine eqadj
c
c=======================================================================
      subroutine eqadj0(amat)
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use vac_mod
      implicit none
! argumnet
      real*8    amat(icvdm2,icvdm2)
! local variables
      integer   i, i1, i2, ir, iv, iz, j, jkladj, jkloop(30), k
     >        , kadjpt, l, ll, n, nn
      real*8    bmat(icvdm2,icvdm2), psvlk(30,0:icvdm2), work
     >        , works, x, xkloop(30), y, ykloop(30)
      save jkloop, kadjpt, bmat, psvlk, xkloop, ykloop
c-----
      data jkladj/0/
c=======================================================================
      if(jkladj.ge.0)then
       iadjpt=iabs(msfx)+1
       do i=1,iadjpt
        rrloop(i)=rvac(i-1)
        zzloop(i)=zvac(i-1)
        pploop(i)=0.d0
        wfloop(i)=1.d0
       enddo
c-----------------------------------------------------------------------
c      SET THE MATRIX EQ FOR COIL CURRENTS
c          BY USING THE LEAST SQUARE FITTING OF FLUX LOOP DATA
c-----------------------------------------------------------------------
       do j=1,iadjpt
        x=(rrloop(j)-rg(1))/dr+1.d0
        y=(zzloop(j)-zg(1))/dz+1.d0
        ir=x
        iz=y
        i=(iz-1)*nr+ir
        x=x-ir
        y=y-iz
        jkloop(j)=i
        xkloop(j)=x
        ykloop(j)=y
       enddo
       do l= 0, icvdm
        do j=1,iadjpt
         i=jkloop(j)
         i1=i+1
         i2=i+nr
         if(i2.gt.nrz)i2=i
         x=xkloop(j)
         y=ykloop(j)
         if(l.ne.0)then
          psvlk(j,l)=svac(i,l)
     &          + x*(svac(i1,l)-svac(i,l))+y*(svac(i2,l)-svac(i,l))
     &          + x*y*(svac(i2+1,l)+svac(i,l)-svac(i1,l)-svac(i2,l))
         else
          psvlk(j,l)=1.d0
         endif
        enddo
       enddo
       do k=0,icvdm
        if(ivac(k).ge.0)then
         do l= 0, icvdm
          if(ivac(l).ge.0)then
           work = 0.d0
           do j=1,iadjpt
            work = work + wfloop(j)*psvlk(j,k)*psvlk(j,l)
           enddo
           bmat(k+1,l+1) = work
          endif
         enddo
        endif
       enddo
       jkladj=1
      endif
c-----
      do k=0,icvdm
       if(ivac(k).ge.0)then
        do l=0,icvdm
         if(ivac(l).ge.0)then
          amat(k+1,l+1)=bmat(k+1,l+1)
         endif
        enddo
        works = 0.d0
        do j=1,iadjpt
         i=jkloop(j)
         i1=i+1
         i2=i+nr
         if(i2.gt.nrz)i2=i
         x=xkloop(j)
         y=ykloop(j)
         works = works
     &     +wfloop(j)*psvlk(j,k)
     &        *(psi(i)+x*(psi(i1)-psi(i))+y*(psi(i2)-psi(i))
     &                +x*y*(psi(i2+1)+psi(i)-psi(i1)-psi(i))-pploop(j))
        enddo 
        amat(k+1,icvdm2) = - works
        if(k.ge.1.and.ivac(k).ge.0)then
         amat(k+1,k+1)=amat(k+1,k+1)+cvacwg(k)
         amat(k+1,icvdm2)=amat(k+1,icvdm2)+cvacwg(k)*cvacst(k)
        endif
       endif
      enddo
c
      return
      end subroutine eqadj0
c
c=======================================================================
      subroutine eqaxi
c=======================================================================
c     SEARCH MAGNETIC AXIS                                         JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, ir, ix, iy, iz, n, idd,idr,idz, ii
      real*8    dra, drrzz, ds, dsrr, dszz, dza, rbv2, rbvax
      real*8    sigcu,sigbt
c*
      sigcu=1.d0
      sigbt=1.d0
      ieqerr(1)=0
      ir=iraxis
      iz=izaxis
      i=(izaxis-1)*nr+iraxis

      ii=0
      idr=1
      idz=0
      idd=1
      do while(idd.ne.0)
       ii=ii+1
       ir=ir+idr
       iz=iz+idz
       i=(iz-1)*nr+ir
       idr=0
       idz=0
       idd=0
       if(sigcu*psi(i).gt.sigcu*psi(i+1))idr= 1
       if(sigcu*psi(i).gt.sigcu*psi(i-1))idr=-1
       if(sigcu*psi(i).gt.sigcu*psi(i+nr))idz= 1
       if(sigcu*psi(i).gt.sigcu*psi(i-nr))idz=-1
       if(idr.ne.0.or.idz.ne.0)idd=1

       if(ii.gt.1000)idd=0
      enddo
c-----------------------------------------------------------------------
      dra =0.5d0*dr*(psi(i+1)-psi(i-1))/(psi(i+1)-2.d0*psi(i)+psi(i-1))
      dza =0.5d0*dz*(psi(i+nr)-psi(i-nr))
     &    / (psi(i+nr)-2.d0*psi(i)+psi(i-nr))
      iaxis=i
      iraxis=mod(i,nr)
      izaxis=i/nr+1
      if(iudsym.ne.0)izaxis=nzh
      ir = iraxis
      iz = izaxis
      raxis=rg(ir)-dra
      zaxis=zg(iz)-dza
      saxis=psi(i)-0.5d0*dr2i*dra*(psi(i+1)-psi(i-1))
     &            -0.5d0*dz2i*dza*(psi(i+nr)-psi(i-nr))
      siw(1)=saxis
c----
      ds=saxis/DBLE(nvm)
      rbv2=btv**2
      do n=nvm,1,-1
       rbv2=rbv2+(fds(n)+fds(n+1))*ds
      enddo
      rbvax=dsqrt(rbv2)*sigbt
      dsrr=(psi(i+1)-2.d0*psi(i)+psi(i-1))*ddri
      dszz=(psi(i+nr)-2.d0*psi(i)+psi(i-nr))*ddzi
      drrzz =dabs(dsrr*dszz)
      qaxis=rbvax/raxis/dsqrt(drrzz)*sigcu
c-----------------------------------------------------------------------
      return
      end subroutine eqaxi
c
c=======================================================================
      subroutine eqbnd
c=======================================================================
c     LOAD PSI VALUE ON BOUNDARY                                   JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      implicit none
! local variables
      integer   is, l, n, ns
      real*8    cs(isrzdm), ds(isrzdm), flx(isrzdm), rs(isrzdm)
     >        , s, sd, smax, ss, x, xx, zs(isrzdm)
c=======================================================================
      if(nsu.lt.nsumax)then
       ns=nsu-1
       do n=1,ns
        rs(n)=0.5*(rsu(n)+rsu(n+1))
        zs(n)=0.5*(zsu(n)+zsu(n+1))
        cs(n)=dsqrt((rsu(n)-rsu(n+1))**2+(zsu(n)-zsu(n+1))**2)
       enddo
c-----
      else
c-----------------------------------------------------------------------
c*    ----rearrenge surface currents
c-----------------------------------------------------------------------
       xx=sdv(nv)/sdw(nv)
       smax=0.d0
       ds(1)=0.d0
       do n=2,nsu
        ds(n)=dsqrt((rsu(n)-rsu(n-1))**2+(zsu(n)-zsu(n-1))**2)
        smax=smax+ds(n)
       enddo
c-----
       sd=smax/DBLE(nsumax-1)
       ns=1
       is=1
       s=0.d0
       ss=0.d0
       rs(1)=rsu(1)
       zs(1)=zsu(1)
       cs(1)=xx*csu(1)*sd
       do while(ns.lt.nsumax)
        ns=ns+1
        s=sd*DBLE(ns-1)
        if(ns.eq.nsumax)s=smax
        do while(ss.lt.s)
         is=is+1
         ss=ss+ds(is)
         if(is.ge.nsu)then
          is=nsu
          ss=smax
         endif
        enddo
        x=(ss-s)/ds(is)
        rs(ns)=rsu(is)+x*(rsu(is-1)-rsu(is))
        zs(ns)=zsu(is)+x*(zsu(is-1)-zsu(is))
        cs(ns)=csu(is)+x*(csu(is-1)-csu(is))
        cs(ns)=xx*cs(ns)*sd
       enddo
       ns=nsumax
c-----
      endif
c-----------------------------------------------------------------------
c*    set the boundary values of psi
c-----------------------------------------------------------------------
      call flxfun(flx,rbnd,zbnd,irzbnd,rs,zs,cs,ns)

      do l=1,irzbnd
       psi(lrzbnd(l))=flx(l)
      enddo
      return
      end subroutine eqbnd
c
c=======================================================================
      subroutine eqdsk0(io)
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      use cnt_mod
      implicit none
! argument
      integer   io
! local variables
      character fname*19
      integer   i, iieqrd, iieqwt, j, kk, kkread, n
      real*8    rsep, zsep
c=======================================================================
      iieqrd=iabs(ieqrd)
      iieqwt=iabs(ieqwt)
      if(io.eq.2)then
c=======================================================================
      if(iieqwt.le.0)return
      inquire(unit=iieqwt,name=fname)
      rsep=rg(irsep)
      zsep=-zg(izsep)
c-----------------------------------------------------------------------
      write(iieqwt)nr,nz
     >,(rg(i),i=1,nr),(zg(j),j=1,nz)
     >,(psi(i),i=1,nr*nz)
     >,nv,(pds(i),i=1,nv),(fds(i),i=1,nv)
     >   ,(vlv(i),i=1,nv),(qqv(i),i=1,nv),(prv(i),i=1,nv)
     >,btv,ttcu,ttpr,bets,beta,betj
     >,(icp(i),i=1,10),(cp(i),i=1,10)
     >,saxis,raxis,zaxis,ell,trg
     >,nsu,(rsu(i),i=1,nsu),(zsu(i),i=1,nsu),(csu(i),i=1,nsu)
     >,isep,dsep,rsep,zsep
     >,(ivac(i),i=0,5)
     >,(rvac(i),i=0,5),(zvac(i),i=0,5)
     >,(cvac(i),i=0,5),(ncoil(i),i=1,5)
     >,((rcoil(i,j),i=1,50),j=1,5)
     >,((zcoil(i,j),i=1,50),j=1,5)
     >,((ccoil(i,j),i=1,50),j=1,5)
     >,ilimt,(rlimt(i),i=1,20),(zlimt(i),i=1,20)
c-----------------------------------------------------------------------
      jeqwt=jeqwt+1
      write(6,'(//5x,a7,a19,a8,i3,a1,i3,a1)')
     >    'diskio=',fname,': write(',ieqwt,'/',jeqwt,')'
c=======================================================================
      else
c=======================================================================
      if(iieqrd.le.0)return
      inquire(unit=iieqrd,name=fname)
      if(jeqrd.eq.0)jeqrd=10000
      do kk=1,jeqrd
      kkread=0
c-----------------------------------------------------------------------
      read(iieqrd,end=88,err=99)nr,nz
     >,(rg(i),i=1,nr),(zg(j),j=1,nz)
     >,(psi(i),i=1,nr*nz)
     >,nv,(pds(i),i=1,nv),(fds(i),i=1,nv)
     >   ,(vlv(i),i=1,nv),(qqv(i),i=1,nv),(prv(i),i=1,nv)
     >,btv,ttcu,ttpr,bets,beta,betj
     >,(icp(i),i=1,10),(cp(i),i=1,10)
     >,saxis,raxis,zaxis,ell,trg
     >,nsu,(rsu(i),i=1,nsu),(zsu(i),i=1,nsu),(csu(i),i=1,nsu)
     >,isep,dsep,rsep,zsep
     >,(ivac(i),i=0,5)
     >,(rvac(i),i=0,5),(zvac(i),i=0,5)
     >,(cvac(i),i=0,5),(ncoil(i),i=1,5)
     >,((rcoil(i,j),i=1,50),j=1,5)
     >,((zcoil(i,j),i=1,50),j=1,5)
     >,((ccoil(i,j),i=1,50),j=1,5)
     >,ilimt,(rlimt(i),i=1,20),(zlimt(i),i=1,20)
      kkread=kk
      enddo
c-----------------------------------------------------------------------
 99   if(kkread.eq.0) then
         write(6,*) '===== DARAIO/READ ERROR ====='
         stop
      endif
c-----------------------------------------------------------------------
  88  if(kkread.eq.0)then
       backspace iieqrd
       write(6,*)'     dataio/end:backspace'
      endif
c-----------------------------------------------------------------------
      if(kkread.lt.jeqrd)jeqrd=kkread-1
      if(iieqwt.eq.iieqrd.and.jeqwt.le.0)jeqwt=jeqrd
      write(6,'(//5x,a7,a19,a8,i3,a1,i3,a1)')
     >    'diskio=',fname,': read (',ieqrd,'/',jeqrd,')'
c-----------------------------------------------------------------------
      rwmn=rg(1)
      rwmx=rg(nr)
      zwmx=-zg(1)
      call eqgrd
c-----------------------------------------------------------------------
      elip=ell
      trig=trg
      iraxis=(raxis-rg(1))/dr+1.
      izaxis=(zaxis-zg(1))/dz+1.
      iaxis=nr*(izaxis-1)+iraxis
      qaxis=qqv(1)
c-----------------------------------------------------------------------
      isep=-iabs(isep)
      irsep=(rsep-rg(1))/dr
      izsep=(-zsep-zg(1))/dz+1.
c-----------------------------------------------------------------------
      tcu=cnmu0*ttcu
      tcur=cnmu0*ttcu
      ccpl=ttcu
c-----------------------------------------------------------------------
      do n=1,nv
       aav(n)=1./rmaj**2
       hdv(n)=btv*aav(n)
      enddo
c=======================================================================
      endif
c=======================================================================
      return
      end subroutine eqdsk0
c
c=======================================================================
      subroutine eqequ
c=======================================================================
c     SET-UP FREE-BOUNDARY EQUILIBRIUM
c                 FOR GIVEN FUNCTIONAL FORMS OF DP/DS AND DT/DS
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use com_mod
      use cnt_mod
      implicit none
c*
      integer   isetup, jsetup
      real*8    pds0(ivdm), fds0(ivdm)
      common/xxpds0/isetup,jsetup,pds0,fds0
! local variables
      integer   i, irend, izend, kset, n, jst
      real*8    bb, cbetj,  ctcu, curfac, e
     >        , fdsx, fdsx0, pr0, pr1, psi0,  tcu0
     >        , tcu1, tpr, ttcur, x, xvr, xx1, xx2, yvr, yvr0
     >        , zvr, zvr0
c=======================================================================
c*    ITERATION LOOP
c-----------------------------------------------------------------------
      jsetup=0
      yvr=0.d0
      zvr=0.d0
      fdsx=0.d0
      do while(jsetup.ge.0)
      jsetup=jsetup+1
c-----------------------------------------------------------------------
c*    LOAD PDS0 AND FDS0
c-----------------------------------------------------------------------
      call eqpds0
c-----------------------------------------------------------------------
      do n=1,nv
       siv(n)=siw(n)
       sdv(n)=sdw(n)
       pds(n)=pds0(n)
       fds(n)=fds0(n)
      enddo
c-----------------------------------------------------------------------
c*    ADJUST BETA-P                        : ICP(2) > 0
c-----------------------------------------------------------------------
      kset=1
      if(jsetup.eq.1)kset=0
      if(isetup.ge.1)kset=1
      if(icp(2).gt.0.and.kset.eq.1)then
c-- ADJUST PDS(I)
       tpr=0.d0
       pr0=0.d0
       do n=nv-1,1,-1
        pr1=pr0
        pr0=pr1-0.5d0*(siw(n+1)-siw(n))*(pds0(n+1)+pds0(n))
        tpr=tpr+0.5d0*(vlv(n+1)-vlv(n))*(pr0+pr1)
       enddo
       ttcur=cnmu0*tcur
       cbetj=cp(1)*rmaj*ttcur**2/(4.*tpr)
       do n=1,nv
        pds(n)=cbetj*pds0(n)
       enddo
c-----------------------------------------------------------------------
c*    ADJUST FDS(I) : <J/R> IS GIVEN       : ICP(2) = 1
c-----------------------------------------------------------------------
       if(icp(2).eq.1)then
        tcu=0.d0
        do n=nv-1,1,-1
         tcu=tcu-0.5d0*(vlv(n+1)-vlv(n))*(fds0(n+1)+fds0(n))
        enddo
        tcu=tcu/cnpi2
        ctcu=ttcur/tcu
        do n=1,nv
         fds(n)=(ctcu*fds0(n)-pds(n))/aav(n)
        enddo
       endif
c-----------------------------------------------------------------------
c*    ADJUST FDS(I) : <J.B> IS GIVEN       : ICP(2) = 2
c-----------------------------------------------------------------------
       if(icp(2).eq.2)then
        tcu0=0.d0
        tcu1=0.d0
        do i=1,nv
         yvr0=yvr
         zvr0=zvr
         xvr=ckv(i)*sdv(i)**2
         if(rbv(i).le.0.)rbv(i)=btv
         bb=rbv(i)**2*aav(i)+xvr
         yvr=xvr/bb*pds(i)
         zvr=rbv(i)*aav(i)*fds0(i)/bb
         if(i.gt.1)then
          tcu0=tcu0-0.5*(vlv(i)-vlv(i-1))*(yvr0+yvr)
          tcu1=tcu1+0.5*(vlv(i)-vlv(i-1))*(zvr+zvr0)
         endif
        enddo
        tcu0=tcu0/cnpi2
        tcu1=tcu1/cnpi2
        curfac=(ttcur-tcu0)/tcu1
        do i=nv,1,-1
         fdsx0=fdsx
         bb=rbv(i)**2*aav(i)+ckv(i)*sdv(i)**2
         fdsx=-(rbv(i)*pds(i)+curfac*fds0(i))/bb
         if(i.eq.nv)then
          rbv(i)=btv
         else
          rbv(i)=rbv(i+1)-0.5*(fdsx+fdsx0)*(siw(i+1)-siw(i))
         endif
         fds(i)=fdsx*rbv(i)
        enddo
       endif
      endif
c-----------------------------------------------------------------------
c*    STORE PSI FOR CONVERGENCE CHECK
c-----------------------------------------------------------------------
      do i=1,nrz
       rbp(i)=psi(i)
      enddo
c-----------------------------------------------------------------------
c*    LOAD RCU=R*CUR(R,Z)
c-----------------------------------------------------------------------
      call eqrcu
c-----------------------------------------------------------------------
c*    SCALE RCU SUCH THAT TCU=TCUR         : ICP(2) > 0
c-----------------------------------------------------------------------
      if(tcur.gt.0.)then
       x=(cnmu0*tcur)/tcu
       tcu=cnmu0*tcur
       do i=1,nrz
        rcu(i)=x*rcu(i)
       enddo
       do n=1,nv
        pds(n)=x*pds(n)
        fds(n)=x*fds(n)
       enddo
      endif
c-----------------------------------------------------------------------
c*    SCALE CSU SUCH THAT CPL=TCU
c-----------------------------------------------------------------------
      cpl=0.d0
      do n=2,nsu
       cpl = cpl + 0.5d0*(csu(n) + csu(n-1))
     >       *dsqrt((rsu(n)-rsu(n-1))**2+(zsu(n)-zsu(n-1))**2)
      enddo
      cpl = cpl + 0.5d0*(csu(1) + csu(nsu))
     >       *dsqrt((rsu(1)-rsu(nsu))**2+(zsu(1)-zsu(nsu))**2)
      x=tcu/cpl
      do n=1,nsu
       csu(n)=x*csu(n)
      enddo
c-----------------------------------------------------------------------
c*    LOAD PSI ON THE BOUNDARY OF THE COMPUTATIONAL DOMAIN
c-----------------------------------------------------------------------
      call eqbnd
c-----------------------------------------------------------------------
c*    SOLVE L(PSI)=RCU(R,Z)
c-----------------------------------------------------------------------
      call eqpde
c-----------------------------------------------------------------------
c*    ADJUST VACCUME FIELDS
c-----------------------------------------------------------------------
      call eqadj
c-----------------------------------------------------------------------
c*    SEARCH MAGNETIC AXIS
c-----------------------------------------------------------------------
      call eqaxi
c-----------------------------------------------------------------------
c*    SEARCH STAGNATION POINT
c-----------------------------------------------------------------------
      call eqsep
      call eqlim
c-----------------------------------------------------------------------
c*    CONVERGENCE CHECK OF INNER LOOP
c-----------------------------------------------------------------------
      erreq=0.d0
      do i=1,nrz
       e=dabs((psi(i)-rbp(i))/saxis)
       if(e.gt.erreq)erreq=e
      enddo
c-----------------------------------------------------------------------
c*    CALCURATE R*BP
c-----------------------------------------------------------------------
      call eqrbp
c-----------------------------------------------------------------------
c*    TRACE OF PLASMA SURFACE
c-----------------------------------------------------------------------
      if(icp(2).gt.0)then
       call eqlin
      else
       psi0=0.d0
       call eqtrc(psi0,raxis,zaxis,-1,0,irend,izend)
       nsu=nzsu
       do n=1,nzsu
        rsu(n)=zrsu(n)
        zsu(n)=zzsu(n)
        csu(n)=zcsu(n)
       enddo
       aav(nv)=zaav
       sdw(nv)=zsdw
      endif
c-----
      qsurf=1.d0/((2.d0*cnpi)**2)*btv*aav(nv)/sdw(nv)
      ttcur=tcu/cnmu0
c-----------------------------------------------------------------------
c*    OUTPUT
c-----------------------------------------------------------------------
      if(iabs(icp(3)).eq.2.and.xxli.gt.0.)then
       call eqind
       erch(1)=dabs(1.d0-zzli/xxli)
       xx1=zzli
       xx2=xxli
      else
       erch(1)=dabs(1.d0-qaxis/qaxi)
       xx1=qsurf
       xx2=qsur
      endif

      jst=jsetup-100*(jsetup/100)
      if(bavmax.ne.1.)erch(1)=bav
      if(isep.eq.0)then
       write(6,'(1x,2i2,2(1x,1pd8.1),3(1x,0pf6.3))')
     >     isetup,jst,erreq,erch(1),qaxis,xx1,ttcur
      else
       write(6,'(1x,2i2,2(1x,1pd8.1),3(1x,0pf6.3),1x,2i3)')                        
     >     isetup,jst,erreq,erch(1)
     >        ,qaxis,xx1,ttcur,irsep,izsep
      endif
c-----------------------------------------------------------------------
      if(ieqerr(1).ne.0)jsetup=-jsetup
      if(erreq.lt.esetup)jsetup=-jsetup
      if(jsetup.ge.msetup)jsetup=-jsetup
      enddo
c-----------------------------------------------------------------------
      return
      end subroutine eqequ
c     
c=======================================================================
      subroutine eqpds0
c=======================================================================
c     LOAD PDS0 AND FDS0
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      implicit none
c=======================================================================
      integer   isetup, jsetup
      real*8    pds0(ivdm), fds0(ivdm)
      common/xxpds0/isetup,jsetup,pds0,fds0
! local variables
      integer   n
      real*8    cbf, cbp, f, x , v
c=======================================================================
c-----CASE:1
      if(iabs(icp(1)).eq.1)then
       do n=1,nv
        x=DBLE(n-1)/DBLE(nvm)
        siw(n)=saxis*(1.d0-x)
        v=vlv(n)/vlv(nv)
        pds0(n)=(1.d0-cp(4))*(1.d0-v**cp(2))**cp(3)+cp(4)
        fds0(n)=(1.d0-cp(8))*(1.d0-v**cp(6))**cp(7)+cp(8)
       enddo
c-----------------------------------------------------------------------
c-----CASE:2
      elseif(iabs(icp(1)).eq.2)then
       do n=1,nv
        x=DBLE(n-1)/DBLE(nvm)
        siw(n)=saxis*(1.d0-x)
        v=vlv(n)/vlv(nv)
        pds0(n)=1.d0-cp(2)*v**cp(3)-(1.d0-cp(2)+cp(4))*v**cp(5)
        fds0(n)=1.d0-cp(6)*v**cp(7)-(1.d0-cp(6)+cp(8))*v**cp(9)
       enddo
c-----------------------------------------------------------------------
c-----CASE:11    =====  HOLLOW CURRENT PROFILE =====
      elseif(iabs(icp(1)).eq.11)then
       do n=1,nv
        x=DBLE(n-1)/DBLE(nvm)
        siw(n)=saxis*(1.d0-x)
        v=vlv(n)/vlv(nv)
        pds0(n)=(1.d0-cp(4))*(1.d0-v**cp(2))**cp(3)+cp(4)
        fds0(n)=(1.d0-v**cp(6))**cp(7)
     >            *(1.+cp(8)*dexp(-((v-cp(9))/cp(10))**2))
       enddo
c-----------------------------------------------------------------------
      elseif(iabs(icp(1)).eq.9)then
       call eqpfds(pds0,fds0)
      endif
c-----------------------------------------------------------------------
      return
      end subroutine eqpds0
c    
c=======================================================================
      subroutine eqerr(errmsg)
c=======================================================================
c     PRINT PSI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      implicit none
! argument
      character errmsg*(*)
! local variables
      integer   i, n
      real*8    psir(irdm)
c=======================================================================
      write(6,'(//5x,a10,a40)')'==========',errmsg
c----------------------------------------------------------------------
      i=nrzh-nr
      do n=1,nr
       psir(n)=psi(i+n)
      enddo
      call prnts(ft06,'psi ',psir,nr,-siw(1))
c----------------------------------------------------------------------
      write(6,*) '===error stop==='
      stop
      end subroutine eqerr
c   
c=======================================================================
      subroutine prnpsi
c=======================================================================
c     PRINT PSI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      implicit none
! local variables
      integer   i, n
      real*8    psir(irdm)
c=======================================================================
      i=nrzh-nr
      do n=1,nr
       psir(n)=psi(i+n)
      enddo
      call prnts(ft06,'psir',psir,nr,1.d0)
      i=-nr/2
      do n=1,nz
       i=i+nr
       psir(n)=psi(i)
      enddo
      call prnts(ft06,'psiz',psir,nz,1.d0)
      return
      end subroutine prnpsi
c
c=======================================================================
      function flux0(x)
c=======================================================================
      implicit none
! argument
      real*8    x
! local variables
      real*8    flux0
      real*8    a0, a1, a2, b0, b1, b2
     >        , c0, c1, c2, d1, d2, pi
     >        , xe, xk, xl, xx
      data a0,a1,a2/1.3862944d0,0.1119723d0,0.0725296d0/
      data b0,b1,b2/0.5d0      ,0.1213478d0,0.0288729d0/
      data c0,c1,c2/1.0d0      ,0.4630151d0,0.1077812d0/
      data    d1,d2/          0.2452727d0,0.0412496d0/
      data       pi/3.141592654d0/
c=======================================================================
      xx=1.d0-x
      xl=dlog(1.d0/xx)
      xk=a0+xx*(a1+xx*a2)+(b0+xx*(b1+xx*b2))*xl
      xe=c0+xx*(c1+xx*c2)+    xx*(d1+xx*d2) *xl
      flux0=((1.d0-x/2.d0)*xk-xe)/(pi*dsqrt(x))
      return
      end function flux0
c  
c=======================================================================
      subroutine flxtab
c=======================================================================
      implicit none
      integer   ntab, mtab
      real*8    dtab, xtab(2050), ftab(2050)
      common/flx0/ntab,mtab,dtab,xtab,ftab
! local variables
      integer   n
      real*8    r, r0, xn, yn
! function
c      real*8    flux0
c=======================================================================
      ntab=2049
      r=4.1d0
      r0=4.2d0
      xn=((r-r0)/(r+r0))**2
      yn=DBLE(ntab)/DBLE(ntab-1)*dlog(1.d0/xn)
      dtab=float(ntab)/yn
      xtab(1)=0.d0
      ftab(1)=0.d0
      do n=2,ntab
       xtab(n)=1.d0-dexp(-DBLE(n-1)/dtab)
       ftab(n)=flux0(xtab(n))
      enddo
c..
cOMMENT      CALL CHKTAB
c..
      return
      end subroutine flxtab
c
c=======================================================================
      function flux(r,z,r0,z0,c0,n)
c=======================================================================
      implicit none
      integer   ntab, mtab
      real*8    dtab, xtab(2050), ftab(2050)
      common/flx0/ntab,mtab,dtab,xtab,ftab
! argument
      integer   n
      real*8    r, z, r0(n), z0(n), c0(n)
! local variables
      real*8    flux
      integer   i, m
      real*8    d, x, xx
! function
c      real*8    flux0
c=======================================================================
      flux=0.d0
      do m=1,n
       x=r*r0(m)
       xx=4.*x/((r+r0(m))**2+(z-z0(m))**2)
       i=1.d0-dtab*dlog(1.d0-xx)
       if(i.lt.ntab)then
        d=(xx-xtab(i))/(xtab(i+1)-xtab(i))
        flux=flux-c0(m)*dsqrt(x)*(ftab(i)+d*(ftab(i+1)-ftab(i)))
       else
        flux=flux-c0(m)*dsqrt(x)*flux0(xx)
       endif
      enddo
      return
      end function flux
c
c=======================================================================
      subroutine flxfun(flx,r,z,m,r0,z0,c0,n)
c=======================================================================
      implicit none
      integer   ntab, mtab
      real*8    dtab, xtab(2050), ftab(2050)
      common/flx0/ntab,mtab,dtab,xtab,ftab
! argument
      integer   m, n
      real*8    flx(m), r(m), z(m), r0(n), z0(n), c0(n)
! local variables
      integer   i, mm, nn
      real*8    d, x, xx
c=======================================================================
      do mm=1,m
       flx(mm)=0.d0
      enddo
      do nn=1,n
      do mm=1,m
       x=r(mm)*r0(nn)
       xx=4.d0*x/((r(mm)+r0(nn))**2+(z(mm)-z0(nn))**2)
       i=1.d0-dtab*dlog(1.d0-xx)
       i=min0(i,ntab-1)
       d=(xx-xtab(i))/(xtab(i+1)-xtab(i))
       flx(mm)=flx(mm)-c0(nn)*dsqrt(x)*(ftab(i)+d*(ftab(i+1)-ftab(i)))
      enddo
      enddo
      return
      end subroutine flxfun
c     
c=======================================================================
      subroutine eqgrd
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use eqv_mod
      implicit none
! local variables
      integer   ir, iz, l, lo, nszx, istop, iss
c=======================================================================
      nvm=nv-1
c-----------------------------------------------------------------------
      nrm=nr-1
      nzm=nz-1
      nrz=nr*nz
      nzh=(nz-1)/2+1
      nrzh=nr*nzh
c-----
      nnr=nr
      nnz=2*nzh-1
      nnrz=nnr*nnz
c-----
      if(iudsym.eq.1)then
      nnz=nzh
      nnrz=nr*nzh
      endif
c..
      rg( 1)= rwmn
      rg(nr)= rwmx
      zg( 1)=-zwmx
      zg(nzh)=0.d0
      dr=(rwmx-rwmn)/float(nrm)
      dz=zwmx/float(nzh-1)
      drz=dr*dz
      dr2i=1.d0/(2.d0*dr)
      dz2i=1.d0/(2.d0*dz)
      ddri=1.d0/(dr*dr)
      ddzi=1.d0/(dz*dz)
      do ir=1,nr
       rg(ir)=rg(1)+dr*float(ir-1)
       rg2(ir)=rg(ir)**2
      enddo
      do iz=1,nz
       zg(iz)=zg(1)+dz*float(iz-1)
      enddo
c.......................................................................
      nsr=nr
      nsz=nz
      nsrm=nsr-1
      nszm=nsz-1
      nsr2=nsrm/2
      nsz2=nszm/2
      nsrz=nr*(nsz-2)
      sf0=-dz*dz/2.d0
      sf1=(dr/dz)**2
      sf2=rg(1)/dr-1.d0
      sf3=(.5d0+.25d0/(2.d0+sf2))/sf1
      sf4=(.5d0-.25d0/(nrm+sf2))/sf1

      csrz(nsz2)=0.d0
      lo=nsz2
      l=lo/2
      csrz(l)=dsqrt(2.d0+csrz(lo))
      lo=l
      csrz(nszm-l)=-csrz(l)
      l=l+2*lo
      istop=0
      do while(istop.eq.0)
       iss=(2*l/nszm)*(2*lo-3)
       if(iss.lt.0)then
        istop=1
       elseif(iss.eq.0)then
        csrz(l)=(csrz(l+lo)+csrz(l-lo))/csrz(lo)
        csrz(nszm-l)=-csrz(l)
        l=l+2*lo
       else
        l=lo/2
        csrz(l)=dsqrt(2.d0+csrz(lo))
        lo=l
        csrz(nszm-l)=-csrz(l)
        l=l+2*lo
       endif
      enddo
      do l=2,nszm
       csrz(l-1)=1.d0/(2.d0+sf1*(2.d0-csrz(l-1)))
      enddo
c-----------------------------------------------------------------------
      irzbnd=0
      do ir=1,nr
       irzbnd=irzbnd+1
       lrzbnd(irzbnd)=ir
       rbnd(irzbnd)=rg(ir)
       zbnd(irzbnd)=zg( 1)
      enddo
c-----
      if(iudsym.eq.0)then
       do ir=1,nr
        irzbnd=irzbnd+1
        lrzbnd(irzbnd)=nr*(nz-1)+ir
        rbnd(irzbnd)=rg(ir)
        zbnd(irzbnd)=zg(nz)
       enddo
      endif
c-----
      nszx=nz-1
      if(iudsym.ne.0)nszx=nzh
      do iz=2,nszx
       irzbnd=irzbnd+1
       lrzbnd(irzbnd)=nr*(iz-1)+1
       rbnd(irzbnd)=rg( 1)
       zbnd(irzbnd)=zg(iz)
      enddo
c-----
      do iz=2,nszx
       irzbnd=irzbnd+1
       lrzbnd(irzbnd)=nr*iz
       rbnd(irzbnd)=rg(nr)
       zbnd(irzbnd)=zg(iz)
      enddo
c-----------------------------------------------------------------------
      return
      end subroutine eqgrd
c    
c=======================================================================
      subroutine eqind
c=======================================================================
c     CALCULATE INDUCTANCE                                         JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      implicit none
! local variables
      integer   i, ir, iz, jrsep, jzsep, n, nn
     >        , irmin, irmax, izmin, izmax
      real*8    bpav, br, bz, ds, p
     >        , surf, volm, sigcu
     >        , rrmin, rrmax, zzmin, zzmax, xlen, rmaj, fac
c      equivalence(psi(1),sai(1,1))
c=======================================================================
      sigcu=1.d0
      jzsep=0
      if(isep.eq.1.and.iisep.gt.0)jzsep=izsep
      jrsep=nr
      if(isep.eq.2.and.iisep.gt.0)jrsep=irsep
c-----------------------------------------------------------------------
      do i=1,nrz
       rbp(i)=psi(i)-cvac(0)
       do n=1,icvdm
        rbp(i)=rbp(i)-cvac(n)*svac(i,n)
       enddo
      enddo
c=======================================================================
c*    INTEGRATE MAGNETIC ENERGY
      zzlp=0.d0
      zzli=0.d0
      surf=0.d0
      volm=0.d0
c*
      zzmin= 10000.d0
      zzmax=-10000.d0
      rrmin= 1000.d0
      rrmax=-1000.d0
      do n=1,nsu
       if(zzmin.gt.zsu(n))zzmin=zsu(n)
       if(zzmax.lt.zsu(n))zzmax=zsu(n)
       if(rrmin.gt.rsu(n))rrmin=rsu(n)
       if(rrmax.lt.rsu(n))rrmax=rsu(n)
      enddo 
      rmaj=0.5d0*(rrmin+rrmax)
      izmin=(zzmin-zg(1))/dz
      izmax=(zzmax-zg(1))/dz+2.d0
      irmin=(rrmin-rg(1))/dr
      irmax=(rrmax-rg(1))/dr+2.d0
      if(iudsym.eq.1)zzmax=nzh
      do iz=izmin,izmax
       do ir=irmin,irmax
        i=(iz-1)*nr+ir
        if(sigcu*psi(i).lt.0.d0)then
          zzlp=zzlp+rbp(i)*rcu(i)/rg(ir)
          br=dz2i*(psi(i+nr)-psi(i-nr))/rg(ir)
          bz=dr2i*(psi(i+ 1)-psi(i- 1))/rg(ir)
          zzli=zzli+(br*br+bz*bz)*rg(ir)
          volm=volm+rg(ir)
          surf=surf+1.d0
        endif
       enddo
      enddo
c*
      surf=surf*drz
      volm=volm*drz*cnpi2
      zzli=zzli*drz*cnpi2
      zzlp=-cnpi2*zzlp*drz/tcu**2
      zzli=zzli/volm
      bpav=0.d0
      xlen=0.d0
      do n=1,nsu
       nn=n-1
       if(nn.eq.0)nn=nsu
       ds=dsqrt((rsu(n)-rsu(nn))**2+(zsu(n)-zsu(nn))**2)
       bpav=bpav+ds*(csu(n)+csu(nn))/2.d0
       xlen=xlen+ds
      enddo
      zzli=zzli/(bpav/xlen)**2
      return
      end subroutine eqind
c  
c=======================================================================
      subroutine eqlim
c=======================================================================
c     LIMITER                                                      JAERI
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      implicit none
! local variables
      real*8    s, u, v
      integer   i, id, ii, ir, iz, n, nd, istop, nn
      real*8    d, d0, da, dx, p, pvbnd, x, y
      real*8    ra,za,rl,zl,rla,zla,rlz2,zla2,rzla,rzla2,zla2r
      real*8    ux,rla2
c-----------------------------------------------------------------------
      s(m,u,v)=psi(m)+(psi(m+1)-psi(m))*u+(psi(m+nr)-psi(m))*v
     >        +(psi(m+nr+1)+psi(m)-psi(m+nr)-psi(m+1))*u*v
c=======================================================================
      pvbnd=100.
c-----
      if(ivbnd.gt.0)then
       pvbnd=0.d0
       do i=1,ivbnd
        x=(rvbnd(i)-rg(1))/dr+1.d0
        y=(zvbnd(i)-zg(1))/dz+1.d0
        ir=x
        iz=y
        istop=0
        if(ir.le.1.or.ir.ge.nr)istop=1
        if(iz.le.1.or.iz.ge.nz)istop=1
        if(istop.eq.0)then
         x=x-ir
         y=y-iz
         ii=nr*(iz-1)+ir
         p=s(ii,x,y)
         if(p.lt.pvbnd)pvbnd=p
        endif
       enddo
      endif
c=======================================================================
      if(ilimt.gt.0)then
c-----------------------------------------------------------------------
      do i=1,ilimt
       x=(rlimt(i)-rg(1))/dr+1.d0
       y=(zlimt(i)-zg(1))/dz+1.d0
       ir=x
       iz=y
       istop=0
       if(ir.le.1.or.ir.ge.nr)istop=1
       if(iz.le.1.or.iz.ge.nz)istop=1
       if(istop.eq.0)then
        x=x-ir
        y=y-iz
        ii=nr*(iz-1)+ir
        p=s(ii,x,y)
        if(p.lt.0.d0)then
         d=1.d+10
         id=0
         do n=1,nzsu
          ux=(rzla*(zzsu(n)-za)+zla2r+rla2*zrsu(n))/rzla2
          if((rl-ux)*(rl-ra).lt.0.d0)then
           id=1
           nn=n
          endif
         enddo
         if(id.eq.1)then
          if(p.lt.pvbnd)pvbnd=p
         endif
        endif
       endif
      enddo
      endif
c-----------------------------------------------------------------------
      if(pvbnd.lt.0.)then
       do i=1,nrz
        psi(i)=psi(i)-pvbnd
       enddo 
       psep=psep-pvbnd
       saxis=saxis-pvbnd
       call eqlin
      endif
c-----------------------------------------------------------------------
      return
      end subroutine eqlim
c     
c=======================================================================
      subroutine eqlin
c=======================================================================
c     LINE INTEGRALS
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, i1, i2, i3, i4, ir, irst, ist, iz, jr
     >        , k, l, ll, lm, lp, n, kstep, irend, izend, istep, istop
      real*8    aa0, aa1, bb0, bb1, bi0, bi1, bp1, ck0, ck1
     >        , dl, dps, dpsi, ds0, ds1, f, ff, psi0, psi00, r0
     >        , r1,  rr0, rr1
     >        , s1, s2, s3, s4, scale, ss0, ss1, vl0, vl1
     >        , x, z0, z1,  ds, sigcu
c=======================================================================
      sigcu=1.d0
      call eqaxi 
      dpsi=saxis/float(nvm)
      irst=iraxis
      ist= (izaxis-1)*nr + irst
      do while(sigcu*psi(ist+1).le.0.d0.and.irst.lt.nr)
        irst=irst+1
        ist= (izaxis-1)*nr + irst
      enddo
      if(irst.eq.nr-1)then
       write(6,*)'stop at eqlin : irst = nr-1'
       stop
      endif
c-----------------------------------------------------------------------
      istep=0
      do n=nv,2,-1
       psi0=dpsi*DBLE(nv-n)
       ist=(izaxis-1)*nr+irst
       do while((psi0-psi(ist))*(psi0-psi(ist+1)).gt.0.d0)
        irst=irst-1
        ist=ist-1
       enddo 
       istop=0
       do while(istop.eq.0)
        call eqtrc(psi0,rg(irst+1),zaxis,1,0,irend,izend)
        if(irend.eq.0)then
         istep=istep+1
         if(istep.gt.20) then
            write(6,*) '  stop at eqlin : istep > 20'
            stop
         endif
         dps=-0.001d0*saxis
         do i=1,nrz
          psi(i)=psi(i)+dps
         enddo
         saxis=saxis+dps
         psep=psep+dps
         cvac(0)=cvac(0)+dps
         call eqaxi
        else
         istop=1
         if(istep.gt.1)write(6,*)'    surface has been moved : istep=',istep
        endif
       enddo
c-----
       siw(n)=psi0
       vlv(n)=zvlv
       arv(n)=zarv
       ckv(n)=zckv
       ssv(n)=zssv
       aav(n)=zaav
       rrv(n)=zrrv
       sdw(n)=zsdw
       rrv(n)=zrrv
       bbv(n)=zbbv
       biv(n)=zbiv
       rpv(n)=zrpv
       rtv(n)=zrtv
       elv(n)=zelv
       dlv(n)=zdlv

       r2b2(n)=zssv*zsdw**2
c-----
       if(n.eq.nv)then
        nsu=nzsu
        do l=1,nsu
         rsu(l)=zrsu(l)
         zsu(l)=zzsu(l)
         csu(l)=zcsu(l)
        enddo
        rsumax=zrsmax
        rsumin=zrsmin
        zsumax=zzsmax
        zsumin=zzsmin
        rzumax=rzsmax
        rzumin=rzsmin
       endif
      enddo
c..evaluate on the axis
      siw(1)=saxis
      arv(1)=0.d0
      vlv(1)=0.d0
      sdw(1)=-dpsi*(vlv(3)**2-2.*vlv(2)**2)
     >       /(vlv(3)*vlv(2)*(vlv(3)-vlv(2)))
      ckv(1)=0.d0
      ssv(1)=0.d0
      aav(1)=1./(raxis*raxis)
      rrv(1)=raxis*raxis
      bbv(1)=2.*bbv(3)-bbv(2)
      biv(1)=1./bbv(1)
      rpv(1)=0.
      rtv(1)=raxis
      elv(1)=elv(2)
      dlv(1)=dlv(2)

      return
      end subroutine eqlin
c     
c=======================================================================
      subroutine eqout
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      use cnt_mod
      use com_mod
      implicit none
! local variables
      integer   i, ii, ir, irend, ix, iy, izend, k, kc, ksep
     >        , lc, mtime1, n
      real*8    cpu, curr(irdm), psi0, psir(irdm)
     >        , rsep, rzmin
     >        , ssr(ivdm), tcv(ivdm), tpr, tps, x, y, zsep
     >        , asp,betp,elipup_l,trigup_l,elipdw,trigdw,zlen
c=======================================================================
      ttcu=tcu/cnmu0
c-----
      psi0=0.d0
      call eqtrc(psi0,raxis,zaxis,-1,0,irend,izend)
      qaxis=qqv(1)
      qsurf=qqv(nv)
      zlen=0.d0
      do n=2,nsu
       zlen=zlen+dsqrt((zrsu(n)-zrsu(n-1))**2+(zzsu(n)-zzsu(n-1))**2)
      enddo
      rmaj=zrtv
      rpla=zrpv
      asp=zrtv/zrpv
      elip=zelV
      trig=zdlv
c--
      elipup_l=(zzsmax-zaxis)/(zrsmax-zrsmin)*2.d0
      trigup_l=(rmaj-rzsmax)/rpla
      elipdw=(zaxis-zzsmin)/(zrsmax-zrsmin)*2.d0
      trigdw=(rmaj-rzsmin)/rpla
c-----
      qqj=5.d0*rpla**2*btv/rmaj**2/ttcu*(1.d0+elip*2)/2.d0
c-----------------------------------------------------------------------
      psi0=0.05d0*saxis
      call eqtrc(psi0,raxis,zaxis,-1,0,irend,izend)
      q95=1.d0/((2.d0*cnpi)**2)*btv*zaav/zsdw
c-----
      psi0=0.05d0*saxis
      x=nvm*(saxis-psi0)/saxis+1.d0
      ix=x
      if(ix.ge.nv)ix=nvm
      x=x-ix
      q95=qqv(ix)+(qqv(ix+1)-qqv(ix))*x
c-----
      el95=(zzsmax-zzsmin)/(zrsmax-zrsmin)
c-----------------------------------------------------------------------
      tpr=0.d0
      tps=0.d0
      do n=nv,1,-1
       prv(n)=muv(n)*hdv(n)**gam
       qqv(n)=1.d0/(4.d0*cnpi*cnpi*nuv(n))
       rbv(n)=hdv(n)/aav(n)
       tcv(n)=ckv(n)*sdv(n)
       qdv(n)=aav(n)*fds(n)
       rho(n)=dsqrt(arv(n)/cnpi)
       if(n.lt.nv)then
        tpr=tpr+0.5d0*(prv(n+1)+prv(n))*(vlv(n+1)-vlv(n))
        tps=tps+0.5d0*(prv(n+1)**2+prv(n)**2)*(vlv(n+1)-vlv(n))
       endif
      enddo
CCC
      write(6,'(A,1p3E12.4)') 'hdv:',hdv(nv-2),hdv(nv-1),hdv(nv)
      write(6,'(A,1p3E12.4)') 'aav:',aav(nv-2),aav(nv-1),aav(nv)
      write(6,'(A,1p3E12.4)') 'brv:',rbv(nv-2),rbv(nv-1),rbv(nv)
c-----
      prfac=prv(1)*vlv(nv)/tpr
      ttpr=tpr/cnmu0
c-----------------------------------------------------------------------
      call deriv(sha,qqv,vlv,nv)
      do n=1,nv
       sha(n)=2.d0*vlv(n)/qqv(n)*sha(n)
      enddo
      call deriv(cuv,tcv,vlv,nv)
c-----
      do n=1,nv
       prv(n)=prv(n)/cnmu0
      enddo
c-----
      cpl=ckv(nv)*sdv(nv)/cnpi2
      bets=2.d0*dsqrt(tps/vlv(nv))/(rbv(nv)/rmaj)**2
      beta=2.d0*tpr/(vlv(nv)*(rbv(nv)/rmaj)**2)
      betj=4.d0*tpr/(rmaj*cpl**2)
      betp=2.d0*tpr/vlv(nv)/(cpl/zlen)**2
c-----
      ccpl=cpl/cnmu0
      bets=100.d0*bets
      beta=100.d0*beta
c-----
      call eqind
c-----
      x=(rloop-rg(1))/dr+1.d0
      y=(zloop-zg(1))/dz+1.d0
      ix=x
      iy=y
      x=x-ix
      y=y-iy
      i=nsr*(iy-1)+ix
      psloop=psi(i)+x*(psi(i+1)-psi(i))+y*(psi(i+nsr)-psi(i))
     >       +x*y*(psi(i+nsr+1)+psi(i)-psi(i+1)-psi(i+nsr))
     >       -cvac(0)
c=======================================================================
c      call clockm(mtime1)
c      cpu=1.d-03*DBLE(mtime1)-cput
c      cput=1.d-03*DBLE(mtime1)
c-----
      if(ieqout(3).ge.1)then
       write(6,*)
       write(6,'(2x,a22,2x,a7,1pd10.3,a8,i4)')
     >   '===== equilibium =====','error=',erreq
       write(6,'(5(2x,a7,f10.5))')
     >    'beta-j=',betj,'beta-t=',beta,'beta-s=',bets
       write(6,'(5(2x,a7,f10.5))')
     >    'beta-p=',betp ,'tpress=',ttpr ,'p0/pav=',prfac
       write(6,'(5(2x,a7,f10.5))')
     >    'cpl   =',ccpl ,'tcur  =',ttcu ,'btv   =',btv
       write(6,'(5(2x,a7,f10.5))')
     >    'qaxis =',qaxis,'qsurf =',qsurf,'q95   =',q95
     >   ,'qqj   =',qqj 
       write(6,'(5(2x,a7,f10.5))')
     >    'rmaj  =',rmaj ,'rpla  =',rpla  ,'vlv   =',vlv(nv)
     >   ,'rp-av =',rho(nv) 
       write(6,'(5(2x,a7,f10.5))')
     >    'asp   =',asp  ,'rpmax =',rsumax ,'rpmin =',rsumin
       write(6,'(5(2x,a7,f10.5))')
     >    'elip  =',elip ,'trig  =',trig  ,'el95  =',el95
       write(6,'(5(2x,a7,f10.5))')
     >    'elipup=',elipup_l,'trigup=',trigup_l
     >   ,'elipdw=',elipdw,'trigup=',trigdw
       write(6,'(5(2x,a7,f10.5))')
     >    'lp    =',zzlp  ,'li    =',zzli
c     >    'lp    =',zzlp  ,'li    =',zzli ,'lam   =',zzlam 
       write(6,'(2x,a7,f10.5,3(2x,a5,i2,f10.5))')
     >    'surf  =',cvac(0),('cvac-',i,cvac(i),i=1,3)
       write(6,'(4(2x,a5,i2,f10.5))')
     >    ('cvac-',i,cvac(i),i=4,icvdm)
      endif
c-----------------------------------------------------------------------
      if( isep.ne.0.and.ieqout(3).ge.1 )then
       if(isep.gt.0) then
        ksep=isep
        isep=-iabs(isep)
        call eqsep
        isep=ksep
       endif
       rsep=0.d0
       zsep=0.d0
       if(irsep.gt.1.and.irsep.lt.nr) rsep= rg(irsep)
       if(izsep.gt.1.and.izsep.lt.nz) zsep= zg(izsep)
       write(ft06,13)rsep,zsep,rspmx,rspmn,psep
  13  format(
     > 2x,'rsep  =',f10.5,2x,'zsep  =',f10.5/2x,'rspmx ='  ,f10.5
     > ,2x,'rspmin =',f10.5,2x,'psep  =',f10.5)
      endif

      if(ieqout(3).ge.1)then
      write(ft06,'(//)')
      endif
c-----------------------------------------------------------------------
      if(ieqout(1).le.0)return
      if(mod(itime,ieqout(1)).ne.0)return
c-----
      if(ieqout(3).ge.2)then
        write(ft06,'(/1H0)')
        call prnts(ft06,'siv ',siv ,nv,siv( 1))
        call prnts(ft06,'sdw ',sdw ,nv,sdv(nv))
        call prnts(ft06,'sdv ',sdv ,nv,sdv(nv))
        call prnts(ft06,'hiv ',hiv ,nv,hiv(nv))
        call prnts(ft06,'hdv ',hdv ,nv,hdv(nv))
        call prnts(ft06,'ckv ',ckv ,nv,ckv(nv))
        call prnts(ft06,'aav ',aav ,nv,aav(nv))
        call prnts(ft06,'muv ',muv ,nv,muv( 1))
        call prnts(ft06,'nuv ',nuv ,nv,nuv( 1))
        call prnts(ft06,'arv ',arv ,nv,arv(nv))
        call prnts(ft06,'vlv ',vlv ,nv,vlv(nv))
        call prnts(ft06,'ssv ',ssv ,nv,ssv(nv))
        call prnts(ft06,'ssr ',ssr ,nv,1.d0 )
        call prnts(ft06,'rrv ',rrv ,nv,rrv(nv))
        call prnts(ft06,'prv ',prv ,nv,prv( 1))
        call prnts(ft06,'qqv ',qqv ,nv, 1.d0)
        call prnts(ft06,'sha ',sha ,nv, 1.d0)
        call prnts(ft06,'cuv ',cuv ,nv,cuv( 1))
        call prnts(ft06,'pds ',pds ,nv,cuv( 1))
        call prnts(ft06,'qdv ',qdv ,nv,cuv( 1))
        call prnts(ft06,'rbv ',rbv ,nv, btv)
        call prnts(ft06,'epv ',epv ,nv, 1.d0)
        call prnts(ft06,'ftr ',ftr ,nv, 1.d0)
      endif
c-----------------------------------------------------------------------
      if(ieqout(3).ge.3)then
       ii=nrzh-nr
       do ir=1,nr
        ii=ii+1
        psir(ir)=psi(ii)
        curr(ir)=rcu(ii)/rg(ir)/cnmu0
       enddo
       call prnts(ft06,'psi ',psir,nr,-saxis)
       call prnts(ft06,'cur ',curr,nr, 0.d0)
      endif
 999  continue
      return
      end subroutine eqout
c     
c=======================================================================
      subroutine eqpde
c=======================================================================
c     SOLVE PDE BY BUNEMAN TECHNIQUE                               JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      implicit none
! local variables
      integer   i, i0, i1, i10, i2, i20, i3, i30, i4, i40, id, ii, io
     >        , ir, iu, iz, j, j0, j2, jd, jdp, jh, jhp, ji, jo, jt
     >        , jtp, k4, ko, l, li, lo, nm2, nszx
     >        , istp,jstp
      real*8    x, a, as, d(irdm), g(irdm), p(irdm)
c-----------------------------------------------------------------------
      if(iudsym.eq.0)then
       nm2=(nszm-1)*nsr
       nszx=nszm
      else
       nm2=nsz2*nsr
       nszx=nsz2+1
      endif
c-----
      do i=1,irdm
        p(i)=0.d0
        d(i)=0.d0
        g(i)=0.d0
      enddo
c
      do ir=2,nsrm
       i0=ir-nr
       do iz=2,nszx
        i=i0+nr*iz
        psi(i)=sf0*rcu(i)
       enddo
      enddo
c
      i10=1-nr
      i20=2-nr
      i30=nsrm-nr
      i40=nsr-nr
      do iz=2,nszx
       i1=i10+nr*iz
       i2=i20+nr*iz
       i3=i30+nr*iz
       i4=i40+nr*iz
       psi(i2)=psi(i2)+sf3*psi(i1)
       psi(i3)=psi(i3)+sf4*psi(i4)
      enddo
c
      lo=nsz2
      ko=2
      id=1
      istp=0
      do while(istp.eq.0)
       li=2*lo
       k4=2*ko-li/nszm
       jd=nsr*nszm/li
       jh=nsr*(nszm/(2*li))
       jt=jd+jh
       ji=2*jd
       jo=jd*ko
c*
       do j=jo,nm2,ji
        j2=j+2
        iu=j+nsrm
        jtp=jt
        jhp=jh
        jdp=jd
        if(iudsym.ne.0.and.j/nm2.ne.0)then
         jdp=-jd
         jhp=-jh
         jtp=-jt
        endif

        if(k4.eq.4)then
         do i=j2,iu
          x=psi(i)-psi(i+jtp)-psi(i-jt)
          psi(i)=psi(i)-psi(i+jhp)-psi(i-jh)+psi(i+jdp)+psi(i-jd)
          p(i-j)=x+psi(i)
         enddo
        elseif(k4.eq.3)then
         do i=j2,iu
          p(i-j)=2.d0*psi(i)
          psi(i)=psi(i+jdp)+psi(i-jd)
         enddo
        elseif(k4.eq.2)then
         do i=j2,iu
          p(i-j)=2.d0*psi(i)+psi(i+jdp)+psi(i-jd)
          psi(i)=psi(i)-psi(i+jhp)-psi(i-jh)
        enddo
        elseif(k4.eq.1)then
         do i=j2,iu
          p(i-j)=2.d0*psi(i)+psi(i+jdp)+psi(i-jd)
          psi(i)=0.d0
         enddo
        endif
c*
        do l=lo,nszm,li
         a=csrz(l)
         as=a*sf1
         do i=2,nsrm
          p(i)=as*p(i)
          d(i)=a-a/(2.d0*(i+sf2))
          g(i)=2.d0*a-d(i)
         enddo
         g(2)=0.d0
         d(nsrm)=0.d0

         jstp=0
         do while(jstp.eq.0)
          ii=2*id
          io=ii+1
          do i=io,nsrm,ii
           a=1./(1.-d(i)*g(i+id)-g(i)*d(i-id))
           p(i)=a*(p(i)+d(i)*p(i+id)+g(i)*p(i-id))
           d(i)=d(i)*d(i+id)*a
           g(i)=g(i)*g(i-id)*a
          enddo
          id=ii
          if(id-nsr2.ge.0)jstp=1
         enddo

         jstp=0
         do while(jstp.eq.0)
          id=ii/2
          io=id+1
          do i=io,nsrm,ii
           p(i)=p(i)+d(i)*p(i+id)+g(i)*p(i-id)
          enddo
          ii=id
          if(id.le.1)jstp=1
         enddo
        enddo
c*
        do i=j2,iu
         psi(i)=psi(i)+p(i-j)
        enddo
       enddo
c*
       if(ko.eq.1)then
        lo=2*lo
        if(lo.ge.nszm)istp=1
       else
        lo=lo/2
        if(lo.eq.1)ko=1
       endif
      enddo
c-------------------------------------------------------------------
      if(iudsym.ne.0)then
       do iz=1,nzh-1
        i0=nr*(iz-1)
        j0=nr*(nsz-iz)
        do ir=1,nr
         i=i0+ir
         j=j0+ir
         psi(j)=psi(i)
        enddo
       enddo
      endif
c-------------------------------------------------------------------
      return
      end subroutine eqpde
c 
c=======================================================================
      subroutine eqrbp
c=======================================================================
c     LOAD RBP                                                    JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      implicit none
! local variables
      integer   i, ir, iz
      real*8    x, y
c-----------------------------------------------------------------------
      i=nrm
      do iz=2,nzm
       i=i+2
       do ir=2,nrm
        i=i+1
        x=dr2i*(psi(i+ 1)-psi(i- 1))
        y=dz2i*(psi(i+nr)-psi(i-nr))
        rbp(i)=dsqrt(x*x+y*y)
       enddo
      enddo
      return
      end subroutine eqrbp
c 
c=======================================================================
      subroutine eqrcu
c=======================================================================
c     LOAD RCU=R*CU(R,Z)                                           JAERI
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, icend, ir, ire, ire0, irs, irs0, itumax, iturn
     >        , iz, izf, izi, izs, jm, jp, irend,izend, istp
     >        , n,irsmin,irsmax,izsmin,izsmax
      real*8    d0, d1, er, fac, facmax, fm, fp, rcu0, rcu1, rcum
     >        , rcumax, x, psi0
c=======================================================================
      baw=1.-bav
      itumax=2
      facmax=1.
      if(iudsym.eq.1)then
      itumax=1
      facmax=2.
      endif
c-----
      d0=-nvm/saxis
      d1=1.-saxis*d0

      psi0=0.d0
      call eqtrc(psi0,raxis,zaxis,-1,0,irend,izend)
      zrsmax=zrsu(1)
      zrsmin=zrsu(1)
      zzsmax=zzsu(1)
      zzsmin=zzsu(1)
      do n=2,nzsu
       if(zrsmax.lt.zrsu(n))zrsmax=zrsu(n)
       if(zrsmin.gt.zrsu(n))zrsmin=zrsu(n)
       if(zzsmax.lt.zzsu(n))zzsmax=zzsu(n)
       if(zzsmin.gt.zzsu(n))zzsmin=zzsu(n)
      enddo
      irsmin=(zrsmin-rg(1))/dr
      irsmax=(zrsmax-rg(1))/dr+2.d0
      izsmin=(zzsmin-zg(1))/dz
      izsmax=(zzsmax-zg(1))/dz+2.d0

      d0=-DBLE(nvm)/saxis
      d1=1.d0-saxis*d0
      tcu=0.d0
      rcumax=0.d0
      do iz=1,nz
       do ir=1,nr
        i=nr*(iz-1)+ir
        rcu0=rcu(i)
        rcu1=0.d0
        rcu(i)=0.d0
        istp=0
        if(iz.le.izsmin)istp=1
        if(iz.ge.izsmax)istp=1
        if(ir.le.irsmin)istp=1
        if(ir.ge.irsmax)istp=1
        if(psi(i).ge.0.d0)istp=1
        if(istp.eq.0)then
         x=d0*psi(i)+d1
         if(x.le.1.d0)x=1.00001d0
         if(x.ge.DBLE(nv))x=DBLE(nv)-1.d-05
         jm=x
         jp=jm+1
         fp=x-DBLE(jm)
         fm=1.d0-fp
         rcu1=-fm*(rg2(ir)*pds(jm)+fds(jm))
     >           -fp*(rg2(ir)*pds(jp)+fds(jp))
        endif
        rcu(i)=bav*rcu1+baw*rcu0
        tcu=tcu +rcu(i)/rg(ir)
        er=dabs(rcu0-rcu(i))
        if(er.gt.error)error=er
        rcum=dabs(rcu(i))
        if(rcum.gt.rcumax)rcumax=rcum
       enddo
      enddo
      tcu=tcu*drz
      error=error/rcumax
c*
      return
      end subroutine eqrcu
c     
c=======================================================================
      subroutine eqsep
c=======================================================================
c     SEARCH SEPARATRIX                                            JAERI
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use vac_mod
      implicit none
! local variables
c      integer   i, i1, i2, i3, i4, idz, ill, ir, irend, irend0, irsmn
c     >        , irsmx, istep, iz, izend, izend0, izsmn, j, jjsep, k
c     >        , kk
c      real*8    amat(icvdm2,icvdm2), bmat(icvdm2,icvdm2), err, psi0
c     >        , rsep, rvac2, x, xsep, zsep, zvac2
c=======================================================================
c----------------------------------------------------------------------
      integer idz,ir,iz,i,is,jjsep,irsmx,j,irsmn,izsmn
      real*8 xsep,err,x,psi0,rsep,zsep
      integer kk,k,irend0,izend0,irend,izend
      integer istep,istp,jstp
      real*8 sigcu
c=======================================================================
      if(isep.eq.0)return
c=======================================================================
      sigcu=1.d0
      idz= 1
      if(iabs(isep).eq.3)idz=-1
c-----
      ir=iraxis
      iz=izaxis-idz*nzh/5
      i=nr*(iz-1)+ir
      iisep=i
      psep=saxis
c-----
      jstp=0
      istep=0
      do while(jstp.eq.0)
       istp=0
       do while(istp.eq.0)
        if(sigcu*psi(i-1).gt.sigcu*psi(i))then
         istp=1
        else
         if(sigcu*psi(i).lt.sigcu*psep)then
          jstp=1
          istp=1
         else
          i=i-1
          ir=ir-1
          if(ir.le.1)then
           istp=1
           iisep=-1
           jjsep=100
           return
          endif
         endif
        endif
       enddo
c-----
       istp=0
       if(jstp.eq.1)istp=1
       do while(istp.eq.0)
        if(sigcu*psi(i+1).gt.sigcu*psi(i))then
         istp=1
        else
         if(sigcu*psi(i).lt.sigcu*psep)then
          jstp=1
          istp=1
         else
          i=i+1
          ir=ir+1
          if(ir.ge.nr)then
           istp=1
           irsep=1
           izsep=1
           iisep=-2
           jjsep=110
           return
          endif
         endif
        endif
       enddo
c-----
       if(istep.eq.0)then
        psep=psi(i)
        iisep=i
        irsep=ir
        izsep=iz
        i=i-nr*idz
        iz=iz-idz
       else
        if(sigcu*psep.ge.sigcu*psi(i))then
         istp=1
         jstp=1
        else
         psep=psi(i)
         iisep=i
         irsep=ir
         izsep=iz
         i=i-nr*idz
         iz=iz-idz
         if(iz.le. 1.or.iz.ge.nz)then
          irsep=1
          izsep=1
          return
         endif
        endif
       endif
c-----
       if(istep.ge.100)jstp=1
       istep=istep+1
      enddo
c---------------------------------------------------------------------
      i=nr*(izaxis-1)+iraxis-1
      j=iraxis
      istp=0
      do while(istp.eq.0)
       irsmx=j
       i=i+1
       if((psi(i)-psep)*(psi(i-1)-psep).le.0.d0)then
        istp=1
        rspmx=rg(irsmx-1)+dr*(psep-psi(i-1))/(psi(i)-psi(i-1))
       else
        j=j+1
        if(j.ge.nr)istp=1
       endif
      enddo
      if(j.ge.nr)then
       irsmx=nr
       iisep=-5
       jjsep=210
      endif
c-----
      i=nr*(izaxis-1)+iraxis+1
      j=iraxis
      istp=0
      do while(istp.eq.0)
       irsmn=j
       i=i-1
       if((psi(i)-psep)*(psi(i+1)-psep).le.0.d0)then
        rspmn=rg(irsmn)+dr*(psep-psi(i))/(psi(i+1)-psi(i))
        istp=1
       else
        j=j-1
        if(j.le.1)istp=1
       endif
      enddo
      if(j.le.1)then
       irsmn=1
       iisep=-6
       jjsep=230
      endif
c======================================================================
      rsep= rg(irsep)
      zsep= zg(izsep)
c-----------------------------------------------------------------------
      xsep=psep+dsep*(saxis-psep)
c-----------------------------------------------------------------------
       j=izsep
       idz=1
       if(izsep.gt.izaxis)idz=-1  
       istp=0
       jstp=0
       do while(istp.eq.0)
        izsmn=j
        i=nr*(j-1)+irsep
        if((psi(i)-xsep)*(psi(i+nr)-xsep).le.0.d0)then
         istp=1 
         jstp=1
        else
         j=j+idz
         if(idz.eq. 1.and.j.gt.izaxis)istp=1
         if(idz.eq.-1.and.j.lt.izaxis)istp=1
        endif
       enddo
       if(jstp.eq.0)then
        izsmn=izaxis
        iisep=-7
        jjsep=340
        return
       endif
c-----------------------------------------------------------------------
      do i=1,nrz
       psi(i)=psi(i)-xsep
      enddo
      cvac(0)=cvac(0)-xsep
      psep=psep-xsep
      call eqaxi
c-----------------------------------------------------------------------
      istp=0
      k=0
      do while(istp.eq.0) 
       kk=k
       psi0=1.d-03*saxis*DBLE(k)
       call eqtrc(psi0,raxis,zaxis,-1,0,irend,izend)
       if(irend.ge.2)then
        if(irend.le.nr-2)then
         if(izend.ge.2)then
          if(izend.le.nz-2)then
           istp=1
           if(k.gt.0)then
            write(6,'(5x,a33,1pd10.3,3i5)')
     >         '===== surface has been changed : ',psi0,k,irend0,izend0
           endif
          endif
         endif
        endif
       endif
       if(istp.eq.0)then
        irend0=irend
        izend0=izend
        k=k+1
        if(k.gt.10000)istp=-1
       endif
      enddo
      if(istp.eq.-1)then
       write(6,fmt='(1x,i4,1pd10.3,2i4)')k,psi0,irend,izend
      endif
      if(kk.gt.0)then
       do i=1,nrz
        psi(i)=psi(i)-psi0
       enddo
       cvac(0)=cvac(0)-psi0
       psep=psep-psi0
       saxis=saxis-psi0
      endif
c-----
      call eqaxi
c-----------------------------------------------------------------------
      return
      end subroutine eqsep
c
c=======================================================================
      subroutine eqsol
c=======================================================================
c     SET-UP INITIAL GUESS
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      use com_mod
      implicit none
c=======================================================================
      integer   isetup, jsetup
      real*8    pds0(ivdm), fds0(ivdm)
      common/xxpds0/isetup,jsetup,pds0,fds0
! local variables
      integer   i, ir, iz, n
      real*8    x
c=======================================================================
c*    LOAD SOLOVEV EQUILIBRIUM AS INITIAL GUESS
      bav=1.d0
      isetup=0
      if(btol.ne.0.0) btv=rmaj*btol
      if(iudsym.eq.1)zpla=0.d0
      bts=btv
      i=0
      do iz=1,nz
       do ir=1,nr
        i=i+1
        x=((1.d0-qsol)*rg(ir)**2+qsol*asol)*(zg(iz)-zpla)**2/esol
     >                +esol/4.d0*(rg(ir)**2-asol)**2-ssol
        psi(i)=x/ssol
       enddo
      enddo
c*
      saxis=-1.
      iraxis=(dsqrt(asol)-rg(1))/dr+1.d0
      izaxis=(zpla-zg(1))/dz+1.0
      if(iudsym.eq.1)izaxis=nzh
      iaxis=nr*(izaxis-1)+iraxis
      raxis=dsqrt(asol)
      zaxis=zpla
      WRITE(6,'(A,3ES12.4)') 'raxis,zaxis,asol   =',raxis,zaxis,asol
      WRITE(6,'(A,3I12)')    'iraxis,izaxis,iaxis=',iraxis,izaxis,iaxis
c-----
      irsep=1
      izsep=1
      iisep=-1
      if(iabs(isep).eq.1)then
       irsep=iraxis
       izsep=izaxis-10
      endif
      if(iabs(isep).eq.3)then
       irsep=iraxis
       izsep=izaxis+10
      endif
      if(iabs(isep).eq.2)then
       irsep=(iraxis+nr)/2
       izsep=izaxis
      endif
      if(iabs(isep).eq.4)then
       irsep=(iraxis)/2
       izsep=izaxis
      endif
c-----------------------------------------------------------------------
c  DUMMY-LOAD TO AVOID THE DIVIDE CHECK IN SUB:EQLIN
c-----------------------------------------------------------------------
      do n=1,nv
       aav(n)=1./raxis**2
       hdv(n)=btv/raxis**2
      enddo
c-----------------------------------------------------------------------
      return
      end subroutine eqsol
c 
c=======================================================================
      subroutine eqtrc(psi0,rst,zst,irsig,irot,irend,izend)
c=======================================================================
c     LINE INTEGRALS
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      implicit none
! argument
      integer   irsig, irot, irend, izend
      real*8    psi0, rst, zst
! local variables
      integer   i, i1, i2, i3, i4, ir, ist, iz, j, k, ll, lm, lp, n
     >        , istop 
      real*8    aa0, aa1, bp1, ck0, ck1, dl, ds0, ds1, psi00, r0, r1
     >        , rr0, rr1, s1, s2, s3, s4, scale, ss0, ss1, vl0, vl1
     >        , x, z0, z1, bb0, bb1, bi0, bi1 , f, ff
c=======================================================================
      nzsu=0
c-----------------------------------------------------------------------
c..SEARCH STARTING POINT
      ir=(rst-rg(1))/dr+1.
      iz=(zst-zg(1))/dz+1.
      i=nr*(iz-1)+ir
      do while((psi0-psi(i))*(psi0-psi(i+1)).gt.0.)
       i=i-irsig
       ir=ir-irsig
       if(ir.le.1.or.ir.ge.nrm)then
        write(ft06,'(1p3d10.3,i4)')psi0,rst,zst,irsig
        call eqerr('==eqtrc:no starting point==')
       endif
      enddo
      ist=i
c-----------------------------------------------------------------------
c..INITIALIZE LINE INTEGRALS

      n=(saxis-psi0)/saxis*DBLE(nv-1)+1.d0
      f=hdv(n)/aav(n) 
      ff=f*f

      zarv=0.
      zvlv=0.
      zsdw=0.
      zckv=0.
      zssv=0.
      zaav=0.
      zrrv=0.
      zelv=0.
      zdlv=0.
      zbbv=0.
      zbiv=0.
      x=(psi0-psi(i))/(psi(i+1)-psi(i))
      r1=rg(ir)+x*dr
      z1=zg(iz)
      bp1=(rbp(i)+x*(rbp(i+1)-rbp(i)))/r1
      if(dabs(bp1).lt.1.d-08)bp1=1.d-08
      vl1=r1*r1
      ds1=1./bp1
      ck1=bp1
      ss1=vl1*bp1
      aa1=ds1/vl1
      rr1=vl1/bp1
c..TRACE CONTURE
      if(irot.eq.0)then
      k=1
      else
      k=3
      iz=iz+1
      i=i+nr
      endif
c-----
      nzsu=1
      zrsu(1)=r1
      zzsu(1)=z1
      zcsu(1)=bp1
      istop=0
      bb1=0.d0
      bi1=0.d0
      do while(istop.eq.0)
       r0=r1
       z0=z1
       vl0=vl1
       ds0=ds1
       ck0=ck1
       ss0=ss1
       aa0=aa1
       rr0=rr1
       bb0=bb1
       bi0=bi1
       i1=i
       i2=i1+1
       i3=i2-nr
       i4=i1-nr
       s1=psi0-psi(i1)
       s2=psi0-psi(i2)
       s3=psi0-psi(i3)
       s4=psi0-psi(i4)
       if(s1*s2.lt.0.0.and.k.ne.1)then
        x=s1/(s1-s2)
        r1=rg(ir)+dr*x
        z1=zg(iz)
        bp1=(rbp(i1)+x*(rbp(i2)-rbp(i1)))/r1
        i=i+nr
        iz=iz+1
        k=3
       elseif(s2*s3.lt.0.0.and.k.ne.2)then
        x=s2/(s2-s3)
        r1=rg(ir+1)
        z1=zg(iz)-dz*x
        bp1=(rbp(i2)+x*(rbp(i3)-rbp(i2)))/r1
        i=i+1
        ir=ir+1
        k=4
       elseif(s3*s4.lt.0.0.and.k.ne.3)then
        x=s4/(s4-s3)
        r1=rg(ir)+dr*x
        z1=zg(iz-1)
        bp1=(rbp(i4)+x*(rbp(i3)-rbp(i4)))/r1
        i=i-nr
        iz=iz-1
        k=1
       elseif(s4*s1.lt.0.0.and.k.ne.4)then
        x=s1/(s1-s4)
        r1=rg(ir)
        z1=zg(iz)-dz*x
        bp1=(rbp(i1)+x*(rbp(i4)-rbp(i1)))/r1
        i=i-1
        ir=ir-1
        k=2
       else
          write(6,*) ' stop at eqtrc : s1*s2*s3*s4=0 '
          stop
       endif
c
       vl1=r1*r1
       if(dabs(bp1).lt.1.d-08)bp1=1.d-08
       ds1=1./bp1
       ck1=bp1
       ss1=vl1*bp1
       aa1=ds1/vl1
       rr1=vl1/bp1
       bb1=ff/(r1*r1)+bp1*bp1
       bi1=1.d0/bb1
       bb1=bb1/bp1
       bi1=bi1/bp1
       dl=dsqrt((r1-r0)*(r1-r0)+(z1-z0)*(z1-z0))
       zvlv=zvlv+(vl0+vl1)*(z0-z1)/2.0
       zsdw=zsdw+dl*(ds0+ds1)/2.0
       zckv=zckv+dl*(ck0+ck1)/2.0
       zssv=zssv+dl*(ss0+ss1)/2.0
       zaav=zaav+dl*(aa0+aa1)/2.0
       zrrv=zrrv+dl*(rr0+rr1)/2.0
       zbbv=zbbv+dl*(bb0+bb1)/2.0
       zbiv=zbiv+dl*(bi0+bi1)/2.0
       zarv=zarv+(r1-r0)*(z1+z0)/2.0
c-----
       nzsu=nzsu+1
       zrsu(nzsu)=r1
       zzsu(nzsu)=z1
       zcsu(nzsu)=bp1
c-----
       if(ir.le.1.or.ir.ge.nrm.or.iz.le.1.or.iz.ge.nzm)istop=1
       if(iudsym.eq.1.and.iz.gt.nzh)istop=1
       if( i.eq.ist )istop=1
      enddo
c-----------------------------------------------------------------------
      irend=ir
      izend=iz
c-----
      if(iudsym.eq.1)then
       zvlv=2.*zvlv
       zsdw=2.*zsdw
       zckv=2.*zckv
       zssv=2.*zssv
       zaav=2.*zaav
       zrrv=2.*zrrv
       zbbv=2.*zbbv
       zbiv=2.*zbiv
       zarv=2.*zarv
       lp=nzsu
       lm=nzsu
       do ll=2,nzsu
        lp=lp+1
        lm=lm-1
        zrsu(lp)= zrsu(lm)
        zzsu(lp)=-zzsu(lm)
        zcsu(lp)= zcsu(lm)
       enddo
       nzsu=2*nzsu-1
      endif
c-----
      zvlv=cnpi*zvlv
      zsdw=1. /(cnpi2*zsdw)
      zckv=cnpi2*zckv/zsdw
      zssv=cnpi2*zssv/zsdw
      zaav=cnpi2*zaav*zsdw
      zrrv=cnpi2*zrrv*zsdw
      zbbv=cnpi2*zbbv*zsdw
      zbiv=cnpi2*zbiv*zsdw
      zarv=dabs(zarv)
      nzsu = nzsu - 1


      zrsmin=zrsu(1)
      zrsmax=zrsu(1)
      zzsmin=zzsu(1)
      zzsmax=zzsu(1)
      rzsmin=zrsu(1)
      rzsmax=zrsu(1)
      do n=2,nzsu
       if(zrsmin.gt.zrsu(n))zrsmin=zrsu(n)
       if(zrsmax.lt.zrsu(n))zrsmax=zrsu(n)
       if(zzsmin.gt.zzsu(n))then
        zzsmin=zzsu(n)
        rzsmin=zrsu(n)
       endif
       if(zzsmax.lt.zzsu(n))then
        zzsmax=zzsu(n)
        rzsmax=zrsu(n)
       endif
      enddo
      zrpv=(zrsmax-zrsmin)/2.d0
      zrtv=(zrsmax+zrsmin)/2.d0
      zelv=(zzsmax-zzsmin)/(2.d0*zrpv)
c
c      zdlv=(rzsmax+rzsmin)/(2.d0*zrpv)
c
      zdlv=(zrtv-(rzsmax+rzsmin)/2.d0)/zrpv
      return
      end subroutine eqtrc
c     
c=======================================================================
      subroutine eqvac
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, icx(10), ii, ir, iz, j, k, nc
      real*8    r2, r2m, r2m2, r2m3, r2m4, r2m5, rm2, rmp, rmp1, rmp2
     >        , rmp3, rmp4, rmp5, z2
c=======================================================================
c-----LOAD ANALYTIC VACUUM FIELD
      do j=1,6
       icx(j)=0
      enddo

      j=0
      do i=1,icvdm
       if(ncoil(i).eq.0.and.ivac(i).ge.0)then
        if(j.lt.6)then
         j=j+1
         if(iudsym.ne.0.and.j.eq.2)j=3
         icx(j)=i
        else
         ivac(i)=-1
         cvac(i)=0.d0
        endif
       endif
      enddo
c-----------------------------------------------------------------------
      rm2=rmaj*rmaj
      rmp=rmaj*rpla
      rmp1=1.d0/(2.d0*rmp)
      rmp2=1.d0/(4.d0*rmp**2)
      rmp3=1.d0/(8.d0*rmp**3)
      rmp4=1.d0/(5.d0*rmp**4)
      rmp5=4.d0/(7.d0*rmp**5)
      i=0
      do iz=1,nzh
       z2=zg(iz)**2
       do ir=1,nr
        r2=rg2(ir)
        r2m=r2-rm2
        r2m2=r2m*r2m
        r2m3=r2m*r2m2
        r2m4=r2m*r2m3
        r2m5=r2m*r2m3
        i=i+1

        if(icx(1).gt.0)
     >   svac(i,icx(1))=rmp1*r2m
        if(icx(2).gt.0.and.iudsym.eq.0)
     >   svac(i,icx(2))=zg(iz)
        if(icx(3).gt.0)
     >   svac(i,icx(3))=rmp2*(4.*r2*z2-r2m2)
        if(icx(4).gt.0)
     >   svac(i,icx(4))=rmp3*(r2*z2*(8.*z2-12.*r2m)+r2m3)
        if(icx(5).gt.0)
     >   svac(i,icx(5))=rmp4*(r2*z2*(z2*(z2-15./4.*(r2-2./3.*rm2))
     >                +15./8.*r2m2)-5./64.*r2m4)
        if(icx(6).gt.0)
     >   svac(i,icx(6))=rmp5*(r2*z2*(z2*(z2*(z2-7./2.*(2.*r2-rm2))
     >                +35./8.*(2.*r2-rm2)*r2m)-35./16.*r2m3)
     >               +7./128.*r2m5)
        if( iz.ne.nzh ) then
         ii = ( nzh+nzh-1-iz)*nr+ir
         if(icx(1).gt.0)
     >    svac(ii,icx(1)) = svac(i,icx(1))
         if(icx(2).gt.0.and.iudsym.eq.0)
     >    svac(ii,icx(2)) =-svac(i,icx(2))
         if(icx(3).gt.0)
     >    svac(ii,icx(3)) = svac(i,icx(3))
         if(icx(4).gt.0)
     >    svac(ii,icx(4)) = svac(i,icx(4))
         if(icx(5).gt.0)
     >    svac(ii,icx(5)) = svac(i,icx(5))
         if(icx(6).gt.0)
     >    svac(ii,icx(6)) = svac(i,icx(6))
        endif
       enddo
      enddo
c-----------------------------------------------------------------------
c*----SOLVE NUMERICAL VACUUM FIELD PRODUCED BY COIL CURRENTS
      do nc=1,icvdm
      if(ncoil(nc).gt.0)then
      call eqvac0(nc,ncoil(nc),rcoil(1,nc),zcoil(1,nc),ccoil(1,nc))
      endif
      enddo
c-----------------------------------------------------------------------
      return
      end subroutine eqvac
c
c=======================================================================
      subroutine eqvac0(nc,nss,rss,zss,uss)
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use vac_mod
      implicit none
! argument
      integer   nc, nss, istop
      real*8    rss(nss), zss(nss), uss(nss)
! local variables
      integer   i, ip, ir, iz, izz, l, n, nn, ns
      real*8    cs(isrzdm), flx(isrzdm), rgs(irdm), rs(isrzdm)
     >        , zgs(irdm), zs(isrzdm)
c=======================================================================
c*    CLEAR RCU(I)
c-----------------------------------------------------------------------
      do i=1,nrz
       rcu(i)=0.d0
      enddo
c-----------------------------------------------------------------------
c*    LOAD COIL CURRENTS
c-----------------------------------------------------------------------
      ns=0
      do n=1,nss
       if(iudsym.eq.1)zss(n)=-dabs(zss(n))
       ns=ns+1
       rs(ns)=rss(n)
       zs(ns)=zss(n)
       cs(ns)=cnmu0*uss(n)
       istop=0
       if(isvac.ge.1)istop=1
       if(rs(ns).lt.rg(1).or.rs(ns).gt.rg(nr))istop=1
       if(zs(ns).lt.zg(1).or.zs(ns).gt.zg(nz))istop=1
       if(istop.eq.0)then
        ir=(rs(ns)-rg(1))/dr+1.5
        iz=(zs(ns)-zg(1))/dz+1.5
        if(ir.eq. 1)ir=  2
        if(ir.eq.nr)ir=nrm
        if(iz.eq. 1)iz=  2
        if(iz.eq.nz)iz=nzm
        i=(iz-1)*nr+ir
        izz= 2*nzh-iz
        rcu(i)=rg(ir)*cs(ns)/drz
        if(device.ne.devnam)
     *   write(ft06,
     *     '(10x,a33,i1,a1,i2,a8,f10.5,a1,f10.5,a5,f10.5,a1,f10.5,a1)')
     *     '=====  COIL POSITION IS MOVED : #'
     *        ,nc,'/',ns,': FROM (',rs(ns),' ,',zs(ns)
     *        ,') TO(',rg(ir),' ,',zg(iz),')'
        rs(ns) = rg(ir)
        zs(ns) = zg(iz)
       endif
      enddo
c-----
      if(iudsym.eq.1)then
       nn=ns
       do n=1,nn
        ns=ns+1
        rs(ns)= rs(n)
        zs(ns)=-zs(n)
        cs(ns)= cs(n)
       enddo
      endif
c-----------------------------------------------------------------------
c*    LOAD PSI BY USING 'EQPDE'
c-----------------------------------------------------------------------
      if(isvac.eq.0) then
       call flxfun(flx,rbnd,zbnd,irzbnd,rs,zs,cs,ns)
       do l=1,irzbnd
        i=lrzbnd(l)
        psi(i)=flx(l)
       enddo
c-----
       call eqpde
c-----
       do i=1,nrz
        svac(i,nc)=psi(i)
        rcu(i)=0.
       enddo
       return
      endif
c-----------------------------------------------------------------------
c*    LOAD PSI BY GREEN FUNCTION
c-----------------------------------------------------------------------
      if(isvac.eq.1) then
       ip = 0
       do iz=1,nz
        do ir=1,nr
         rgs(ir)=rg(ir)
         zgs(ir)=zg(iz)
        enddo
        call flxfun(flx,rgs,zgs,nsr,rs,zs,cs,ns)
        do ir=1,nsr
         ip = ip + 1
         psi(ip)=flx(ir)
        enddo
       enddo
       do i = 1,nrz
        svac(i,nc) = psi(i)
       enddo
       return
      endif
c-----------------------------------------------------------------------
      return
      end subroutine eqvac0
c
c=======================================================================
      subroutine eqvol
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      use com_mod
      
      implicit none
! local variables
      integer   n
      real*8    dsi, rbv2
c=======================================================================
c-----LINE INTEGRALS
c-----------------------------------------------------------------------
      call eqrbp
      call eqlin
c-----------------------------------------------------------------------
c-----LOAD VALUABLES ON EQUIL GRIDS:SIV(N),N=1,NV
c-----           SIV : POLOIDAT FLUX FUNCTION
c-----------------------------------------------------------------------
      dsi=saxis/DBLE(nvm)
      do n=nv,1,-1
       if(n.eq.nv)then
        prv(n)=0.d0
        rbv2=btv*btv
       else
        prv(n)=prv(n+1)+0.5*(pds(n)+pds(n+1))*dsi
        rbv2=rbv2+(fds(n)+fds(n+1))*dsi
       endif
       rbv(n)=dsqrt(rbv2)
       siv(n)=siw(n)
       sdv(n)=sdw(n)
       hdv(n)=rbv(n)*aav(n)
       write(6,'(I8,1P2E12.4)') n,hdv(n),rbv(n)
       muv(n)=prv(n)/hdv(n)**gam
       nuv(n)=sdv(n)/hdv(n)
       qqv(n)=1.d0/(cnpi2**2*nuv(n))
      enddo
      qqv(1)=qaxis
      nuv(1)=1.d0/(cnpi2**2*qaxis)
      hiv(1)=0.d0
      do n=2,nv
       hiv(n)=hiv(n-1)-0.5*(1./nuv(n)+1./nuv(n-1))*dsi
      enddo
c=======================================================================
      return
      end subroutine eqvol
c
c=======================================================================
      subroutine eqeqv
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      integer   n
      real*8    btv2, dsi, rbv2
c=======================================================================
      call eqrbp
      call eqlin
c=======================================================================
c-----LOAD VALUABLES ON EQUIL GRIDS:SIV(N),N=1,NV
      dsi=saxis/DBLE(nvm)
      btv2=btv*btv
      do n=nv,1,-1
       if(n.eq.nv)then
        prv(n)=2.d0*cnmu*premin
        rbv2=btv2
       else
        prv(n)=prv(n+1)+0.5d0*(pds(n)+pds(n+1))*dsi
        rbv2=rbv2+(fds(n)+fds(n+1))*dsi
       endif
       rbv(n)=dsqrt(rbv2)
       siv(n)=siw(n)
       sdv(n)=sdw(n)
       hdv(n)=rbv(n)*aav(n)
       muv(n)=prv(n)/hdv(n)**gam
       nuv(n)=sdv(n)/hdv(n)
      enddo
      nuv(1)=1.d0/(cnpi2**2*qaxis)
      qqv(1)=1.d0/(cnpi2**2*nuv(1))
      hiv(1)=0.d0
      do n=2,nv
       qqv(n)=1.d0/(cnpi2**2*nuv(n))
       hiv(n)=hiv(n-1)-0.5d0*(1.d0/nuv(n)+1.d0/nuv(n-1))*dsi
      enddo
c-----------------------------------------------------------------------
c-----SAVE EQ-VARIABLES FOR ERROR CHECK
      do n=1,nv
       save(n, 1)=hiv(n)
       save(n, 2)=hdv(n)
       save(n, 3)=siv(n)
       save(n, 4)=sdv(n)
       save(n, 5)=siw(n)
       save(n, 6)=sdw(n)
       save(n, 7)=vlv(n)
       save(n, 8)=ckv(n)
       save(n, 9)=ssv(n)
       save(n,10)=aav(n)
       save(n,11)=rrv(n)
       save(n,12)=pds(n)
       save(n,13)=fds(n)
      enddo
c=======================================================================
      return
      end subroutine eqeqv
c
      end module eqsub_mod
