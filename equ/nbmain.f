      block data nbdata
      implicit none
      include 'AAA'
      include 'NBI'
      include 'NBIC'
c
      data icnnbi,jcnnbi/0,0/
      data dtnbi/0.02/
      data nbiion,fmsnbi,chgnbi/1,1.,1./
      data pnbi,enbi,rnbi,dnbi,sig/20.,100.,0.20,0.20,0./
      data inbi,jnbi,knbi,lnbi/10,10,2,0/
      data mnbi,nnbi/0,0/
      data nbidbg/5*0/
c
c<<injection geometry>>
      data xstnbi /  4.75000,  4.75000,  4.75000,  4.75000 ,4*0./
      data ystnbi /  0.83500,  0.83500,  0.83500,  0.83500 ,4*0./
      data zstnbi /  1.36950,  1.46570, -1.36950, -1.46570 ,4*0./
      data th1nbi / -1.0,     -1.0,     -1.0,     -1.0     ,4*0./
      data th2nbi / 33.6024,  37.3976,  -33.6024, -37.3976 ,4*0./
      data ebnbi  / 100.0,    100.0,    100.0,    100.0    ,4*0./
      data frnbi  /   1.0,      1.0,      1.0,      1.0    ,4*0./
      data wsp1nb /   0.850,    0.850,    0.850,    0.850  ,4*0./
      data wsp2nb /   0.090,    0.090,    0.090,    0.090  ,4*0./
      data wsp3nb /   0.060,    0.060,    0.060,    0.060  ,4*0./
      data nbsrc/4/, nbeng/3/
      end
c
c=======================================================================
      subroutine nbini
c=======================================================================
c     nbi heating
c=======================================================================
      implicit none
      include 'AAA'
      include 'NBI'
      include 'NBIC'
      include 'CNT'
c>>
      real*8    tx(nbidm),ty(nbidm),tz(nbidm)
      integer   ismnb
      common /tvect/tx,ty,tz
      common /smtnb/ismnb
! local variables
      integer   i, ien, ii, ix, iy, j, jen, k, kstepk, l
      real*8    ang1, ang2, cang1, cang2, cos1, cos2, dnbi2, dx, dy
     >        , ratio, rnbi2, rr, sang1, sang2, sin1, sin2, tbeam
     >        , work(100), wtot, x, xst, y, yst, zst
      data kstepk/0/
c>>
c-----------------------------------------------------------------------
      namelist/nbi/dtnbi,nbiion,fmsnbi,chgnbi
     >,pnbi,enbi,rnbi,dnbi,inbi,jnbi,knbi,sig
     >,xstnbi,ystnbi,zstnbi,th1nbi,th2nbi,ebnbi,frnbi,nbsrc,nbeng
     >,wsp1nb,wsp2nb,wsp3nb
     >,mnbi,nnbi
     >,timnbi,pownbi,dtmnbi
     >,nbidbg,icnnbi,jcnnbi
     >,ismnb
c=======================================================================
c::read data
      if(kstepk.le.0)rewind ft05
      kstepk=1
      read(ft05,nbi,end=888,err=999)
      if(nbsrc.gt.4) nbsrc=4
      if(nbeng.gt.3) nbeng=3
c
      write(ft06,601) nbiion,fmsnbi,chgnbi
     >               ,pnbi,enbi,rnbi,dnbi
  601 format(/2x,'namelist &nbi'/
     >/5x,'ion  =',i5,5x,2x,'mass =',f10.4,2x,'charg=',f10.4,
     >/5x,'pnbi =',f10.4,2x,'enbi =',f10.4,2x,'rnbi =',f10.4,
     > 2x,'dnbi =',f10.4)
c
c::species
      do i=1,nbsrc
      wspcnb(i,1)=0.0
      wspcnb(i,2)=0.0
      wspcnb(i,3)=0.0
      wtot=0.0
      do j=1,nbeng
      if(j.eq.1) wspcnb(i,1)=wsp1nb(i)
      if(j.eq.2) wspcnb(i,2)=wsp2nb(i)
      if(j.eq.3) wspcnb(i,3)=wsp3nb(i)
      wtot=wtot+wspcnb(i,j)
      enddo
      if(wtot.eq.0.0 .or. ebnbi(i).eq.0.0 .or. frnbi(i).le.0.0) then
      frnbi(i)=0.0
      else
      do j=1,nbeng
      wspcnb(i,j)=wspcnb(i,j)/wtot
      enddo
      endif
      enddo
c
      wtot=0.0
      do i=1,nbsrc
      wtot=wtot+frnbi(i)
      enddo
      do i=1,nbsrc
      frnbi(i)=frnbi(i)/wtot
      enddo
c
c::energy
      jen=0
      do 140 i=1,nbsrc
      kengnb(i,1)=0
      if(frnbi(i).gt.0.0)then
      if(jen.eq.0) then
      jen=jen+1
      engynb(jen)=ebnbi(i)
      kengnb(i,1)=jen
      else
      do 150 j=1,jen
      ratio=ebnbi(i)/engynb(j)
      if(ratio.ge.0.99 .and. ratio.le.1.01) then
      kengnb(i,1)=j
      go to 140
      endif
  150 continue
      jen=jen+1
      engynb(jen)=ebnbi(i)
      kengnb(i,1)=jen
      endif
      endif
  140 continue
      engynb(0)=0.0
c
c::species
      ien=jen
      do i=2,nbeng
      do j=1,jen
      ien=ien+1
      engynb(ien)=engynb(j)/i
      enddo
      iengnb=ien
      do k=1,nbsrc
      if(kengnb(k,1).ne.0)then
      kengnb(k,i)=jen*(i-1)+kengnb(k,1)
      endif
      enddo
      enddo
c
      write(6,602) (engynb(i),i=1,iengnb)
  602 format(/2x,'  engynb =',6f10.3)
c
      write(ft06,603)
  603 format(/2x,'  injection geometry  '
     >       /10x,'nsrc',7x,'xst',9x,'yst',9x,'zst',9x,'th1'
     >        ,9x,'th2',9x,'eb ',9x,'fr ')
      do i=1,nbsrc
      write(6,604) i,xstnbi(i),ystnbi(i),zstnbi(i)
     >              ,th1nbi(i),th2nbi(i),ebnbi(i),frnbi(i)
      enddo
      write(ft06,605)
  605 format(/2x,'  power ratio '
     >       /10x,'nsrc',7x,'eb1',9x,'eb2',9x,'eb3',12x,'pb1'
     >        ,9x,'pb2',9x,'pb3')
      do i=1,nbsrc
      if(frnbi(i).eq.0.0) then
      write(6,606) i,(0.0,0.0,j=1,nbeng)
      else
      write(6,606) i,(engynb(kengnb(i,j)),j=1,3),
     >(wspcnb(i,j),j=1,3)
      endif
      enddo
  604 format(12x,i1,3x,7(f10.4,2x))
  606 format(12x,i1,3x,3(f10.4,2x),3x,3(f10.4,2x))
c
      write(ft06,12)inbi,jnbi,knbi,mnbi,nnbi
     >          ,(timnbi(j),j=0,5)
     >          ,(pownbi(j),j=0,5)
     >          ,(dtmnbi(j),j=0,5)
     >          ,(nbidbg(j),j=1,5)
  12  format(
     >/5x,'inbi =',i5,7x,'jnbi =',i5,7x,'knbi =',i5,7x,2i5/
     >/5x,'timnbi=',6(1x,f10.3)
     >/5x,'pownbi=',6(1x,f10.3)
     >/5x,'dtmnbi=',6(1x,f10.3)
     >/5x,'dbg  =',5i2/)
c
c
c::beam start point
      do i=1,4
      kst1nb(i)=0
      kst2nb(i)=0
      enddo
      ii=0
      do 400 k=1,nbsrc
      if(frnbi(k).le.0.0) go to 400
      kst1nb(k)=ii+1
      xst=xstnbi(k)
      yst=ystnbi(k)
      zst=zstnbi(k)
      ang1=cnpi/180.d0*th1nbi(k)
      ang2=cnpi/180.d0*th2nbi(k)
      sin1=dsin(ang1)
      cos1=dcos(ang1)
      sin2=dsin(ang2)
      cos2=dcos(ang2)
      rnbi2=rnbi**2
      dnbi2=dnbi**2
      dy=2.*rnbi/dfloat(inbi)
      dx=2.*rnbi/dfloat(jnbi)
      tbeam=0.
      do 499 iy=1,inbi
 499  work(iy)=0.
      do 410 iy=1,inbi
      do 410 ix=1,jnbi
      x=dx*(dfloat(ix)-0.5)-rnbi
      y=dy*(dfloat(iy)-0.5)-rnbi
      rr=x*x+y*y
      if(rr.gt.rnbi2) go to 410
      ii=ii+1
c
      if(ii.gt.nbptl) then
      write(ft06,882) ii,nbptl
  882 format(/2x,'dimension error at sub. nbini   ii.gt.nbptl  ',2i5)
      stop
      endif
c
      fbeam(ii)=dexp(-rr/dnbi2)
      work(iy)=work(iy)+fbeam(ii)
      tbeam=tbeam+fbeam(ii)
c-----------------------------------------------------------------------
      xsbeam(ii)=xst+x*sin1-y*cos1*sin2
      ysbeam(ii)=yst+x*cos1+y*sin1*sin2
      zsbeam(ii)=zst       +y*cos2
c--
      sang1=sin(ang1)
      cang1=cos(ang1)
      sang2=sin(ang2)
      cang2=cos(ang2)
      dxbeam(ii)= cang2*cang1
      dybeam(ii)=-cang2*sang1
      dzbeam(ii)= sang2
 410  continue
      kst2nb(k)=ii
      do 420 l=kst1nb(k),kst2nb(k)
      fbeam(l)=fbeam(l)/tbeam
 420  continue
c
      tx(k)=-cang2*cang1
      ty(k)= cang2*sang1
      tz(k)=-sang2
c
 400  continue
c-----------------------------------------------------------------------
      lnbi=ii
      return
c-----------------------------------------------------------------------
 888  write(6,*)'>>&nbi : end / pnbi = 0<<'
      pnbi=0
      rewind ft05
      return
 999  write(6,*)'>>&nbi : err / pnbi = 0<<'
      pnbi=0
      rewind ft05
      return
      end
c
c=======================================================================
      subroutine nbcrs
c=======================================================================
c<<cross-section for nbi>>
c
      implicit none
      include 'AAA'
      include 'R2D'
      include 'NBI'
      include 'NBIC'
      include 'TRN'
      include 'CRS'
! local variables
      integer   i, ip, k, m, mp, n
      real*8    dlnbi, dnhyd, dnimp, tte, tti, vnbi, x, xm, xrcx
     >        , xrcx0, xrcx1, xrei, xrei0, xrei1, xrii, xrion, xrpi
     >        , xrpi0, xrpi1
c
      dlnbi=drr2d/dfloat(knbi)
      do 50 k=1,iengnb
      enbi=engynb(k)
      vnbi=dsqrt(2.*cnec*1.d+03*enbi/cnmp/fmsnbi)
      xm=dlog(vnbi/vrcrs(1))/dvcrs+1.
      m=xm
      xm=xm-m
      mp=m+1
      do 100 n=1,ntr2dm
      tte=ter2d(n,0)
      tti=ter2d(n,1)
      if(tte.lt.tmcrs(1))tte=tmcrs(1)
      if(tti.lt.tmcrs(1))tti=tmcrs(1)
      x=dlog(tte/tmcrs(1))/dtcrs+1.d0
      i=x
      ip=i+1
      x=x-i
      xrei0=rcei(i,m)*(rcei(ip,m)/rcei(i,m))**x
      xrei1=rcei(i,mp)*(rcei(ip,mp)/rcei(i,mp))**x
      xrei=xrei0*(xrei1/xrei0)**xm
c-----
      x=dlog(tti/tmcrs(1))/dtcrs+1.d0
      i=x
      ip=i+1
      x=x-i
      xrpi0=rcpi(i,m)*(rcpi(ip,m)/rcpi(i,m))**x
      xrpi1=rcpi(i,mp)*(rcpi(ip,mp)/rcpi(i,mp))**x
      xrpi=xrpi0*(xrpi1/xrpi0)**xm
c-----
      xrcx0=rccx(i,m)*(rccx(ip,m)/rccx(i,m))**x
      xrcx1=rccx(i,mp)*(rccx(ip,mp)/rccx(i,mp))**x
      xrcx=xrcx0*(xrcx1/xrcx0)**xm
c-----
      xrii=xrpi+xrcx
c-----
      dnhyd=(fchrg(nion+1)-zef2d(n))/(fchrg(nion+1)-1.d0)
      dnimp=(zef2d(n)-1.d0)/(fchrg(nion+1)-1.d0)/fchrg(nion+1)
c-----
      xrion=xrei+(dnhyd+dnimp*fchrg(nion+1)**1.4)*xrii
c-----
      if(sig.lt.0.)xrion=-sig
      sigvnb(n,k)=der2d(n,0)*xrion*dlnbi
      if(nbidbg(1).gt.0)then
      write(65,fmt='(2x,1p5d10.3)')ter2d(n,0),ter2d(n,1),xrei,xrcx,xrpi
      endif
 100  continue
c
  50  continue
c
      return
      end
c
c=======================================================================
      subroutine nbint
c=======================================================================
c     nbi heating
c=======================================================================
      implicit none
      include 'AAA'
      include 'R2D'
      include 'NBI'
      include 'NBIC'
      include 'TRN'
      include 'CRS'
      include 'CNT'
      real*8    tx(nbidm),ty(nbidm),tz(nbidm)
      real*8    axi(itdm)
      common /tvect/tx,ty,tz
      common /nbixi/axi
! local variables
      integer   i, ien, ii, iii, ir, is, ispmax, isptb(50), ist, it, itt
     >        , iz, j, jeng, k, ksrc, m, meng, msec0, n
      real*8    a, b_r, b_t, b_x, b_y, b_z, db, df, df0, dlnbi, dxi
     >        , dy, dyi, dzi, eb, fac, fb0, ffnbi, fvv(itdm), fvv0
     >        , path, pathl, pb, rr, rvl, s, spcoe(50,4), sum0, sum1
     >        , tht, tnbis, vect_x(1000,nbidm,3), vect_y(1000,nbidm,3)
     >        , vect_z(1000,nbidm,3), x, xa, xi, xl, xxsp(501), y, yi
     >        , yy, yysp(501), z, z0work(10), z1work(10), z2work(10)
     >        , za, zav, zdvl, zi, ztrp(3), zz
      real*8    u, v, sfun
c-----------------------------------------------------------------------
      sfun(m,u,v)=rho2d(m)+(rho2d(m+1)-rho2d(m))*u
     >                    +(rho2d(m+nrr2d)-rho2d(m))*v
     >                    +(rho2d(m+nrr2d+1)+rho2d(m)
     >                     -rho2d(m+nrr2d)-rho2d(m+1))*u*v
c=======================================================================
c
c::clear
      do 110 n=1,ntr2dm
        do 110 i=1,4
          do 110 j=1,3
            fhnbi(n,i,j)=0.d0
 110  continue
c-----------------------------------------------------------------------
      dlnbi=drr2d/dfloat(knbi)
      tnbis=0.d0
c-----------------------------------------------------------------------
      do 150 ksrc=1,nbsrc
        do 170 jeng=1,nbeng
          meng=kengnb(ksrc,jeng)
          if(meng.eq.0) go to 170
          ist=kst1nb(ksrc)
          ien=kst2nb(ksrc)
          if(ist.eq.0 .or. ien.eq.0) go to 170
c
c<<modified>>
          do 239 n=1,ntr2dm
            fvv(n)=0.d0
 239      continue
c>><<
          do 240 i=ist,ien
            xi=xsbeam(i)
            yi=ysbeam(i)
            zi=zsbeam(i)
            dxi=dxbeam(i)
            dyi=dybeam(i)
            dzi=dzbeam(i)
            xl=-dlnbi
            tnbis=tnbis+fbeam(i)
            do 200 iii=1,1000
              xl=xl+dlnbi
              x=xi-xl*dxi
              y=yi-xl*dyi
              z=zi-xl*dzi
              xa=(dsqrt(x*x+y*y)-rsr2d)/drr2d+1.d0
              za=(z-zsr2d)/dzr2d+1.d0
              ir=xa
              iz=za
              if(ir.lt.nrr2d-2.and.iz.gt.2.and.iz.le.nzr2d-1)goto 210
 200        continue
 210        continue
            xa=xa-ir
            za=za-iz
            ii=(iz-1)*nrr2d+ir
            is=0
            path=0.d0
            df=fbeam(i)
            pathl=0.d0
c-----------------------------------------------------------------------
            do 220 j=1,20000
              xl=xl+dlnbi
              x=xi-xl*dxi
              y=yi-xl*dyi
              z=zi-xl*dzi
              xa=(dsqrt(x*x+y*y)-rsr2d)/drr2d+1.d0
              za=(z-zsr2d)/dzr2d+1.d0
              ir=xa
              iz=za
              xa=xa-ir
              za=za-iz
              ii=(iz-1)*nrr2d+ir
              s=sfun(ii,xa,za)
              itt=s
              if(itt.le.0) itt=1
c<<modified>>
              if(nnbi.eq.0)then
                if(itt.le.ntr2dm)then
                  is=1
                  df0=df
                  path=path-sigvnb(itt,meng)
                  pathl=pathl+dlnbi
                  df=fbeam(i)*dexp(path)
                  fhnbi(itt,ksrc,jeng)=fhnbi(itt,ksrc,jeng)+df0-df
c>>
                  rr=dsqrt(x*x+y*y)
                  zz=z
                  call magfld(rr,zz,b_r,b_z,b_t)
c                 b_t=btv/rr
                  if(x.eq.0.d0)then
                    if(y.ge.0)then
                      tht=cnpi/2.d0
                    else
                      tht=3.d0*cnpi/2.d0
                    endif
                  else
                    tht=atan(y/x)
                  endif
                  b_x=b_r*cos(tht)+b_t*sin(tht)
                  b_y=b_r*sin(tht)-b_t*cos(tht)
                  a=dsqrt(b_x**2+b_y**2+b_z**2)
                  vect_x(i,ksrc,jeng)=b_x/a
                  vect_y(i,ksrc,jeng)=b_y/a
                  vect_z(i,ksrc,jeng)=b_z/a
                  axi(itt)=tx(ksrc)*vect_x(i,ksrc,jeng)
     >                    +ty(ksrc)*vect_y(i,ksrc,jeng)
     >                    +tz(ksrc)*vect_z(i,ksrc,jeng)
c>>
                else
                  if(is.gt.0)goto 230
                endif
              else
                s=(s-1.d0)/dfloat(ntr2dm)
                s=s*s*dfloat(ntr2dm)
                it=s+1.d0
                if(it.le.ntr2dm)then
                  is=1
                  df0=df
                  path=path-sigvnb(itt,meng)
                  pathl=pathl+dlnbi
                  df=fbeam(i)*dexp(path)
                  fvv(it)=fvv(it)+df0-df
                else
                  if(is.gt.0)goto 230
                endif
              endif
 220        continue
 230        continue
 240      continue
c<<modified>>
          if(nnbi.ne.0)then
            x=0.d0
            do 241 n=1,ntr2dm
              x=x+fvv(n)
              fvv(n)=x
 241        continue
            m=1
            dy=1.d0/dfloat(ntr2dm)
            y=dy
            fvv0=0.d0
            fb0=0.d0
            do 242 n=1,ntr2dm-1
              x=dfloat(n)/dfloat(ntr2dm)
              x=x*x
 243          if(x.gt.y)then
                fvv0=fvv(m)
                m=m+1
                y=dy*dfloat(m)
                goto 243
              endif
              yy=fvv0+(fvv(m)-fvv0)*(x-(y-dy))/dy
              fhnbi(n,ksrc,jeng)=yy-fb0
              fb0=yy
 242        continue
            fhnbi(ntr2dm,ksrc,jeng)=fvv(ntr2dm)-fb0
          endif
c<<modified-2>>
          if(mnbi.gt.0)then
            do 244 n=1,ntr2d
              xxsp(n)=ro(n)
              if(n.eq.ntr2dm)then
                fhnbi(n,ksrc,jeng)=0.5d0*fhnbi(n-1,ksrc,jeng)
              endif
              yysp(n)=fhnbi(n,ksrc,jeng)
 244        continue
            ispmax=5
            call csplsf(xxsp,yysp,spcoe,sig,2.d-03
     >                                     ,ntr2d,ispmax,10,isptb)
            if(sig.gt.0.d0)then
              sum0=0.d0
              sum1=0.d0
              do 245 n=1,ntr2d
                if(yysp(n).lt.0.d0)yysp(n)=0.d0
                sum0=sum0+yysp(n)
                sum1=sum1+fhnbi(n,ksrc,jeng)
 245          continue
              fac=sum1/sum0
              do 246 n=1,ntr2d
                x=fac*yysp(n)
 246            fhnbi(n,ksrc,jeng)=fac*yysp(n)
            endif
          endif
c-----
  170   continue
  150 continue
c
c::center value
c
      do 400 i=1,nbsrc
        do 410 j=1,nbeng
          zav=(fhnbi(1,i,j)+fhnbi(2,i,j))/(vlr2d(3)-vlr2d(1))
          fhnbi(1,i,j)=zav*(vlr2d(2)-vlr2d(1))
          fhnbi(2,i,j)=zav*(vlr2d(3)-vlr2d(2))
  410   continue
  400 continue
c
c
c::output birth profile
c
      if(nbidbg(2).le.0) return
      write(65,601)
      do 811 jeng=1,nbeng
      z0work(jeng)=0.d0
      z1work(jeng)=0.d0
 811  continue
      do 822 jeng=1,nbeng
      do 821 ksrc=1,nbsrc
      meng=kengnb(ksrc,jeng)
      if(meng.eq.0) go to 822
      eb=1.d+03*engynb(meng)
      pb=1.d+06*pnbi*frnbi(ksrc)*wspcnb(ksrc,jeng)
      db=pb/(cnec*eb)
      z0work(jeng)=z0work(jeng)+db
 821  continue
 822  continue
      do 803 n=1,ntr2dm
      do 812 jeng=1,nbeng
      z2work(jeng)=0.d0
 812  continue
      z2work(nbeng+1)=0.d0
      do 802 jeng=1,nbeng
      do 801 ksrc=1,nbsrc
      meng=kengnb(ksrc,jeng)
      if(meng.eq.0) go to 802
      eb=1.d+03*engynb(meng)
      pb=1.d+06*pnbi*frnbi(ksrc)*wspcnb(ksrc,jeng)
      db=pb/(cnec*eb)
c-----
      ffnbi=fhnbi(n,ksrc,jeng)
      z2work(jeng)=z2work(jeng)+db*ffnbi/(vlr2d(n+1)-vlr2d(n))
      z1work(jeng)=z1work(jeng)+db*ffnbi
 801  continue
      z2work(nbeng+1)=z2work(nbeng+1)+z2work(jeng)
 802  continue
      zdvl=vlr2d(n+1)-vlr2d(n)
      write(65,fmt='(1x,i5,1p10d10.3)')n,(z2work(j),j=1,nbeng+1),zdvl
 803  continue
      do 813 jeng=1,nbeng
 813  z2work(jeng)=z1work(jeng)/z0work(jeng)
      n=0
      write(65,fmt='(1x,i5,1p10d10.3)')n,(z2work(j),j=1,nbeng)
c
 601  format(//2x,80('=')/5x,'beam birth profile')
      do 300 k=1,nbsrc
      write(65,602) k,frnbi(k),kst1nb(k),kst2nb(k)
     >  ,(engynb(kengnb(k,j)),j=1,3)
 602  format(/5x,'ion source =',i2,'  power =',f10.4,'  start point =',
     >        2i8/8x,'  energy =',3f10.4)
      ztrp(1)=0.d0
      ztrp(2)=0.d0
      ztrp(3)=0.d0
      do 310 n=1,ntr2dm
      rvl=(vlr2d(n+1)-vlr2d(n))/vlr2d(ntr2d)
      do 315 j=1,3
      ztrp(j)=ztrp(j)+fhnbi(n,k,j)
 315  continue
      write(65,603) n,ro(n),(fhnbi(n,k,j)/rvl,j=1,3)
 603  format(2x,i3,f10.4,3x,3f10.4)
 310  continue
      write(65,604) (ztrp(j),j=1,3)
 604  format(2x,3x,'    trap  ',3x,3f10.4)
 300  continue
c
      return
      end
c=======================================================================
      subroutine nbmain
c=======================================================================
c     nbi heating
c=======================================================================
      implicit none
      include 'AAA'
      include 'NBI'
      include 'CNT'
      include 'NBIC'
      include 'EQT'
! local variables
      logical*1 calc
      integer   i, is, itnbi, itnbi0, kk, ktnbi, n
      real*8    time0
      save itnbi, itnbi0, ktnbi, time0
c-----------------------------------------------------------------------
      data calc/.false./
c=======================================================================
c((no detail output))----'out=1' at trout
      if(out.eq.1)return
c-----------------------------------------------------------------------
      if(step.eq.0)then
      if(icnnbi.eq.0)then
      call nbini
      endif
      pnbi=0.
      do 100 i=0,10
 100  if(pnbi.lt.pownbi(i))pnbi=pownbi(i)
      if(pnbi.le.0.)calc=.true.
      itnbi=-100
      ktnbi=0
      do 110 kk=1,100
      if(time.ge.timnbi(ktnbi))goto 120
 110  ktnbi=ktnbi+1
 120  continue
      dtnbi=dtmnbi(ktnbi)
      time0=time
      return
      endif
c=======================================================================
      if(calc)return
      if(icnnbi.ge.0)then
c-----
 130  continue
      if(time.gt.timnbi(ktnbi+1))then
      ktnbi=ktnbi+1
      dtnbi=dtmnbi(ktnbi)
      else
      goto 140
      endif
      goto 130
 140  continue
c-----
      pnbi=pownbi(ktnbi)
     >         +(pownbi(ktnbi+1)-pownbi(ktnbi))
     >          *(time-timnbi(ktnbi))/(timnbi(ktnbi+1)-timnbi(ktnbi))
      if(pnbi.le.0.)return
c-----
      itnbi0=itnbi
      itnbi=(time-timnbi(ktnbi))/dtnbi+1
c-----
      if(jcnnbi.eq.1)itnbi0=-100
      if(jcnnbi.eq.2)time0=-time
c-----
      if(itnbi.gt.itnbi0.or.itnbi.le.0)then
      call nbcrs
      call nbint
      call nbcal
      time0=time
c-----
      elseif(time.gt.time0)then
      call nbcal
      time0=time
      endif
c-----
      jcnnbi=0
      endif
c=======================================================================
c<<output source terms for transport equations
c<<    include SOC 
c     tsorc=tsorc+tnbic
c     do 900 n=1,nt
c900  ssorc(n)=ssorc(n)+snbic(n)
c     is=nbiion
c     tsord(is)=tsord(is)+tnbid
c     do 910 n=1,nt
c910  ssord(n,is)=ssord(n,is)+snbid(n,is)
c     do 920 is=0,indmz
c     tsore(is)=tsore(is)+tnbie(is)
c     do 920 n=1,nt
c     ssore(n,is)=ssore(n,is)+snbie(n,is)
c920  continue
      return
      end
c
c
c=======================================================================
      subroutine nbcal
c=======================================================================
c
c     interface of nbmain to 1d f.p. solver
c
c=======================================================================
      implicit none
      include 'AAA'
      include 'EQT'
      include 'TRN'
      include 'CNT'
      include 'NBI'
      include 'NBIC'
c     include 'NTR'
      include 'R2D'
      include 'COM'
c
      real*8    snbie_w(itdm,0:indmz,nbidm,3)
      integer   ismnb
      common /nbs0/ snbie_w
      common /smtnb/ismnb
! local variables
      integer   iistp, is, itfst, ivmax, jeng, ksrc, m, meng, n, nbeam
      real*8    db, dnbm, dnnbi_w(itdm), dtfst, dvl, eb, ebeam(20)
     >        , emax, exmx, fchbm, fmsbm, fnblos, pb, prbm, sbeam(20)
     >        , sbmd, sbme(0:indmz), taul
      save      itfst, ivmax, dtfst, exmx, fnblos, taul
c-----------------------------------------------------------------------
      data iistp/0/
      namelist/fp1dnb/ivmax,fnblos,exmx,taul,dtfst,itfst
c=======================================================================
      if(iistp.eq.0)then
      write(6,*)'                                                   >>'
      write(6,*)'   <<nb heating is evaluated by using the solution >>'
      write(6,*)'   <<      of time dependent 1.d fokker planck eq. >>'
      write(6,*)'                                                   >>'
      iistp=1
      ivmax=201
      fnblos=1.d0
      exmx=1.1d0
      taul=0.001d0
      dtfst=0.001d0
      itfst=5
      timold=-1.d+10
      rewind ft05
      read(ft05,fp1dnb,end=88,err=99)
  99  continue
  88  continue
      write(6,fp1dnb)
      endif
c=======================================================================
      if(time.lt.timold)return
      dtfst=dtime/dfloat(itfst)
      timold=time+dtime
c=======================================================================
c::clear
      tnbid=0.d0
      tpnbi=0.d0
      tdnbi=0.d0
      do 100 is=0,mion
      tnbie(is)=0.d0
      do 100 n=1,ntr2d
      dnnbi(n,is)=0.d0
      snbid(n,is)=0.d0
 100  snbie(n,is)=0.d0
c-----
      tnbic=0.d0
      do 110 n=1,ntr2d
      prnbi(n)=0.d0
 110  snbic(n)=0.d0
c-----------------------------------------------------------------------
      do 230 n=1,ntr2dm
      dvl=vlt(n+1)-vlt(n)
      nbeam=0
      emax=0.d0
      do 200 ksrc=1,nbsrc
      do 200 jeng=1,nbeng
      meng=kengnb(ksrc,jeng)
      if(meng.eq.0) go to 200
      eb=1.d+03*engynb(meng)
      pb=1.d+06*pnbi*frnbi(ksrc)*wspcnb(ksrc,jeng)
      db=pb/(cnec*eb)
      fnbi(n)=fhnbi(n,ksrc,jeng)
      nbeam=nbeam+1
      sbeam(nbeam)=db*fhnbi(n,ksrc,jeng)/dvl
      ebeam(nbeam)=eb
      if(eb.gt.emax)emax=eb
c-----
 200  continue
c-----
      fmsbm=fmass(nbiion)
      fchbm=fchrg(nbiion)
      call nbfp1d(n,ivmax,fnblos,emax,taul,dtfst,itfst
     >           ,nbeam,ebeam,sbeam,fmsbm,fchbm
     >           ,dnbm,prbm,sbme,sbmd)
c-----
      dnnbi(n,nbiion)=dnbm
      prnbi(n)=prbm
      snbid(n,0)=sbmd*dvl
      snbid(n,nbiion)=sbmd*dvl
      do 210 m=0,mion
      snbie(n,m)=sbme(m)*dvl
c
      snbie_w(n,m,1,1)=sbme(m)*dvl
c
  210 continue
c-----
      tdnbi=tdnbi+dnnbi(n,nbiion)*dvl
      tpnbi=tpnbi+prnbi(n)*dvl
      dnnbi(n,0)=dnnbi(n,nbiion)
      tnbid=tnbid+snbid(n,nbiion)
      do 220 is=0,jonr2d
      tnbie(is)=tnbie(is)+snbie(n,is)
  220 continue
  230 continue
c-----
      if(ismnb.eq.1) then
        do 500 is=0,1
          call fitfun(ntm,roh,dnnbi(1,is),dnnbi_w)
          do n=1,ntm
            dnnbi(n,is)=dnnbi_w(n)
          enddo
  500   continue
      endif
   42 format(1x,i3,2(1x,1pe10.3))
   43 format(1x,a3,2(1x,a10))
c-----
c::calc nbcd
      call nbsocc
c-----
      return
      end
c
c=======================================================================
      subroutine nbfp1d(n,ifvmax,flos,emax,taul,dtfst,itfst
     >           ,nbeam,ebeam,sbeam,fmsbm,fchbm
     >           ,dnbm,prbm,sbme,sbmd)
c=======================================================================
c
c     1-d fokker planck solver
c
c=======================================================================
      implicit none
      include 'AAA'
      include 'EQT'
      include 'TRN'
c=======================================================================
      integer   ifpdim
      parameter(ifpdim=401)
      integer   kfvmax
      real*8    fvfst(ifpdim,itdm),vfst(ifpdim),enruni
      common/nbifv0/kfvmax
      common/nbifv1/fvfst,vfst,enruni
! argument
      integer   n, ifvmax, itfst, nbeam
      real*8    flos, emax, taul, dtfst, ebeam(20), sbeam(20)
     >        , fmsbm, fchbm, dnbm, prbm, sbme(0:indmz), sbmd
! local variables
      integer   iist, iistp, ir, is, iv, ivlos, ivsor(20), jv, jvsor
     >        , kvsor, nb
      real*8    aa(ifpdim), bb(ifpdim), cc(ifpdim), dd(ifpdim)
     >        , dif(ifpdim), dnb, dnbl, dts, dv, e, ecrit, elos, enbil
     >        , enbl, enr, ex, ex3, fac(0:indmz), flam, fls(ifpdim)
     >        , flw(ifpdim), se, si, sor(ifpdim), tslow, tthem, vc3
     >        , vcrit, vmax, xfv(ifpdim), xfv0(ifpdim), xl, xsb, zav
      data enbil / 0.0d0 /
      save ivsor, kvsor, dv
c-----------------------------------------------------------------------
      data iistp/0/
c=======================================================================
      if(iistp.eq.0)then
      iistp=1
      kfvmax=ifvmax
      enruni=0.d0
      do 100 nb=1,nbeam
  100 if(enruni.lt.ebeam(nb))enruni=ebeam(nb)
      emax=1.5d0*enruni
      vmax=dsqrt(emax/enruni)
      dv=vmax/dfloat(kfvmax-1)
      do 110 iv=1,kfvmax
      vfst(iv)=dv*dfloat(iv-1)
      do 110 ir=1,nt
      fvfst(iv,ir)=0.d0
 110  continue
c((source terms))
      kvsor=0
      do 130 nb=1,nbeam
      do 120 iv=1,kfvmax
      e=vfst(iv)**2*enruni
      if(e.ge.ebeam(nb))then
      ivsor(nb)=iv
      if(kvsor.lt.iv)kvsor=iv
      goto 130
      endif
 120  continue
 130  continue
      endif
c-----------------------------------------------------------------------
c     clear
c-----------------------------------------------------------------------
      dnbm=0.d0
      prbm=0.d0
      sbmd=0.d0
      do 200 is=0,indmz
      sbme(is)=0.d0
 200  continue
c-----------------------------------------------------------------------
      flam=16.d0
      zav=0.d0
      do 210 is=1,mion
      fac(is)=fchrg(is)**2*denh(n,is)/denh(n,0)*fmsbm/fmass(is)
 210  zav=zav+fac(is)
      do 220 is=1,mion
 220  fac(is)=fac(is)/zav
      ecrit=14.8d0*fmsbm**(1.d0/3.d0)*zav**(2.d0/3.d0)*temh(n,0)
      ex=dsqrt(enruni/ecrit)
      ex3=ex**3
      tslow =0.125d0*fmsbm*(temh(n,0)/1.d+03)**1.5d0
     >     /(fchbm**2*(denh(n,0)/1.d+19)*(flam/16.d0))
      tthem=tslow /3.d0*dlog(1.d0+ex3)
      elos=flos*temh(n,1)
      vcrit=dsqrt(ecrit/enruni)
      vc3=vcrit**3
      ivlos=dsqrt(elos/enruni)/dv+1.d0
c=======================================================================
      do 230 iv=1,kfvmax
 230  xfv(iv)=fvfst(iv,n)
c=======================================================================
      do 800 iist=1,itfst
c=======================================================================
c     temporal xfv
c=======================================================================
      if(tthem.gt.10.d0*dtfst)then
      do 240 iv=1,kfvmax
      xfv0(iv)=xfv(iv)
      sor(iv)=0.d0
      fls(iv)=0.d0
      if(iv.lt.ivlos)fls(iv)=1.d0/taul
      if(iv.gt.1)then
      dif(iv)=temh(n,0)/(2.d0*enruni)*(vfst(iv)**3+vc3)/vfst(iv)
      flw(iv)=vfst(iv)**3+vc3
      endif
 240  continue
      dif(1) = 0.0d0
      flw(1) = 0.0d0
      do 250 nb=1,nbeam
      jvsor=ivsor(nb)
      sor(jvsor)=sor(jvsor)
     >          +sbeam(nb)/(vfst(jvsor)**2*dv)
 250  continue
c-----
      dts=dtfst/tslow
      do 260 iv=2,kfvmax-1
      aa(iv)=-dts*(dif(iv+1)/dv+flw(iv+1)/2.d0)/(2.d0*vfst(iv)**2*dv)
      bb(iv)=1.d0+dts*(dif(iv)+dif(iv))/(2.d0*dv**2*vfst(iv)**2)
     >         +dtfst*fls(iv)/2.d0
      cc(iv)=-dts*(dif(iv-1)/dv-flw(iv-1)/2.d0)/(2.d0*vfst(iv)**2*dv)
      dd(iv)=-aa(iv)*xfv0(iv+1)+(2.d0-bb(iv))*xfv0(iv)
     >       -cc(iv)*xfv0(iv-1)+dtfst*sor(iv)
 260  continue
      xfv(1)=0.d0
      xfv(kfvmax)=0.d0
      call recur(xfv,aa,bb,cc,dd,kfvmax,0)
c=======================================================================
c     steady xfv(iv)
c=======================================================================
      else
      do 270 iv=1,kfvmax
 270  xfv(iv)=0.d0
      do 290 nb=1,nbeam
      jv=ivsor(nb)
      xsb=sbeam(nb)*6.d0/(vfst(jv+1)**3-vfst(jv-1)**3)
     >        *(vfst(jv)**2*dv+dv**3/3.d0)
      xsb=sbeam(nb)*tslow
      do 280 iv=ivlos+1,jv
      xfv(iv)=xfv(iv)+xsb/(vfst(iv)**3+vc3)
 280  continue
      xl=0.d0
      do 285 iv=2,ivlos
      fls(iv)=1.d+05
 285  xl=xl+fls(iv)*vfst(iv)**2*dv
      xl=sbeam(nb)/xl
      do 286 iv=2,ivlos
 286  xfv(iv)=xfv(iv)+xl
 290  continue
      endif
c=======================================================================
 800  continue
c=======================================================================
c     integral values
c=======================================================================
      dnb=0.d0
      enr=0.d0
      dnbl=0.d0
      enbl=0.d0
      se=0.d0
      si=0.d0
      do 300 iv=2,kvsor
      if(iv.gt.ivlos)then
      dnb=dnb+xfv(iv)*vfst(iv)**2*dv
      enr=enr+xfv(iv)*vfst(iv)**4*dv
      se=se-0.5d0*2.d0*(temh(n,0)/(2.d0*enruni*vfst(iv))-1.d0)
     >       *(xfv(iv)*vfst(iv)**4+xfv(iv-1)*vfst(iv-1)**4)*dv
      si=si-0.5d0*2.d0*(temh(n,0)/(2.d0*enruni*vfst(iv))-1.d0)*vc3
     >       *(xfv(iv)*vfst(iv) +xfv(iv-1)*vfst(iv-1) )*dv
      else
      dnbl=dnbl+xfv(iv)*vfst(iv)**2*dv*fls(iv)
      enbl=enbl+xfv(iv)*vfst(iv)**4*dv*fls(iv)
      endif
 300  continue
      enr=enr*cnec*enruni*2.d0/3.d0
      se=se/tslow*cnec*enruni
      si=si/tslow*cnec*enruni
      enbil=enbil*cnec*enruni
c----
      sbme(0)=se
      do 310 is=1,mion
      sbme(is)=fac(is)*si
 310  continue
      sbme(1)=sbme(1)+enbil
      prbm=enr
      dnbm=dnb
      sbmd=dnbl
c=======================================================================
      do 320 iv=1,kfvmax
 320  fvfst(iv,n)=xfv(iv)
c=======================================================================
      return
      end
c*deck nbmain
c=======================================================================
      subroutine nbsocc
c=======================================================================
      implicit none
      include 'AAA'
      include 'TRN'
      include 'NBI'
      include 'NBIC'
      include 'EQT'
      include 'CNT'
      include 'R2D'
c>>
      real*8    axi(itdm)
      real*8    zzflam(itdm,0:indmz,0:indmz)
      real*8    snbie_w(itdm,0:indmz,nbidm,3)
      common /nbixi/axi
      common /lamda/zzflam
      common /nbs0/ snbie_w
! local variables
      integer   is, jeng, ksrc, m, meng, n
      real*8    axih(itdm), eb, epsil, pbm 
     >        , vb, vc, x, xjnbi(itdm), y, zefh(itdm)
     >        , zzflamh(itdm,0:indmz,0:indmz)
      integer   i
      real*8    funcj, funcg
c-----------------------------------------------------------------------
      funcj(x,y)=x**2/(4.d0+3.d0*y+x**2*(x+1.39d0+0.61d0*y**0.7d0))
      funcg(i,x)=dsqrt(x)*(1.55d0+0.85d0/zef(i))
     >                                 -x*(0.2d0+1.55d0/zef(i))
c-----------------------------------------------------------------------
      do 100 n=1,ntr2dm
        m=n+1
        if(m.eq.ntr2d)m=n
        axih(n)=0.5d0*(axi(n)+axi(m))
        zefh(n)=0.5d0*(zef(n)+zef(m))
        zzflam(n,0,1)=30.9-dlog(dsqrt(denh(n,0))/temh(n,0))
        zzflam(m,0,1)=30.9-dlog(dsqrt(denh(m,0))/temh(m,0))
        zzflamh(n,0,1)=0.5d0*(zzflam(n,0,1)+zzflam(m,0,1))
 100  continue
c-----------------------------------------------------------------------
      do 400 n=1,ntr2dm
        xjnbi(n)=0.d0
  400 continue
      do 500 ksrc=1,nbsrc
        do 510 jeng=1,nbeng
          meng=kengnb(ksrc,jeng)
          if(meng.eq.0) go to 500
          eb=1.d+03*engynb(meng)
          vb=dsqrt(2.d0*eb*cnec/(fmsnbi*cnmp))
          do 520 n=1,ntr2dm
            vc=1.42d+6*dsqrt(temh(n,0)*1.d-3)
            x=vb/vc
            y=4.d0*zefh(n)/(5.d0*fmsnbi)
            epsil=rph(n)/rth(n)
            xjnbi(n)=xjnbi(n)+15.82d0*(temh(n,0)/10.d3)
     >                   *axih(n)/(chgnbi*denh(n,0)*1.d-20)
     >                   *(17.d0/zzflamh(n,0,1))*funcj(x,y)
     >                   *(1.d0-(chgnbi/zefh(n))*(1.d0-funcg(n,epsil)))
            pbm=0.d0
            do 530 is=0,jonr2d
              pbm=pbm+snbie_w(n,is,ksrc,jeng)
  530       continue
            snbic(n)=snbic(n)+xjnbi(n)*pbm/(vlt(n+1)-vlt(n))
  520     continue
  510   continue
  500 continue
      do 550 n=1,ntr2dm
        tnbic=tnbic+snbic(n)*(art(n+1)-art(n))
  550 continue
c-----
      return
      end
c
      subroutine magfld(rr,zz,br,bz,bt)
c                        i  i  o  o
c
      implicit none
      include 'AAA'
      include 'EQU'
      include 'EQV'
      include 'GEO'
      include 'R2D'
! argument
      real*8    rr, zz, br, bz, bt
! local variables
      integer   i, ixr, ixz, is
      real*8    xbr1, xbr2, xbr3, xbr4, xbz1, xbz2, xbz3, xbz4, xr, xz
      real*8    si
c
c
      xr=(rr-rg(1))/dr+1.d0
      xz=(zz-zg(1))/dz+1.d0
      ixr=xr
       xr=xr-ixr
      ixz=xz
       xz=xz-ixz
      i  =nzr2d*(ixz-1)+ixr
c
      dr2i=1.d0/(2.d0*drr2d)
      dz2i=1.d0/(2.d0*dzr2d)
c
c<br>::
      xbr1=-(psi(i+nzr2d)-psi(i-nzr2d))*dz2i
      xbr2=-(psi(i+1+nzr2d)-psi(i+1-nzr2d))*dz2i
      xbr3=-(psi(i+2*nzr2d)-psi(i))*dz2i
      xbr4=-(psi(i+1+2*nzr2d)-psi(i+1))*dz2i
      br=xbr1+(xbr2-xbr1)*xr+(xbr3-xbr1)*xz
     >        +(xbr4-xbr3-xbr2+xbr1)*xr*xz
      br=br/rr
c<bz>::
      xbz1=(psi(i+1)-psi(i-1))*dr2i
      xbz2=(psi(i+2)-psi(i))*dr2i
      xbz3=(psi(i+1+nrr2d)-psi(i-1+nrr2d))*dr2i
      xbz4=(psi(i+2+nrr2d)-psi(i+nrr2d))*dr2i
      bz=xbz1+(xbz2-xbz1)*xr+(xbz3-xbz1)*xz
     >        +(xbz4-xbz3-xbz2+xbz1)*xr*xz
      bz=bz/rr
c
      si=psi(i)+xr*(psi(i+1)-psi(i))
     >         +xz*(psi(i+nr)-psi(i))
      is=dfloat(nv-1)*(saxis-si)/saxis+1.d0
      if(is.gt.nv.or.is.lt.1)is=nv
      bt=rbv(is)/rr
c
      return
      end
