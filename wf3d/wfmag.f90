!	  wfmag.f90 MOD. by T. Yokoyama for TASKX/WFX
!     Jan./09/2013	based on wfmag.f90
!     ******* CALCULATION OF GAMMA10 MAGNETIC FIELD ******
!
      subroutine mag(x,y,z,bx,by,bz)
      implicit real*8 (a-h,o-z)
!
      common /cmline/ xline(1200),yline(1200),zline(1200),cline( 200),&
                     idl12( 200),idl23( 200),nline
      common /cmloop/ xloop( 200),yloop( 200),zloop( 200),aloop( 200),&
                     cloop( 200),dzlp( 200),drlp( 200),idzlp( 200),&
                     idrlp( 200),alfa( 200),beta( 200),nloop
      common /cmarc/  xarc( 200),yarc( 200),zarc( 200),aarc( 200),&
                     carc( 200),dzarc( 200),drarc( 200),idzarc( 200),&
                     idrarc( 200),alfarc( 200),betarc( 200),&
                     phai1( 200),phai2( 200),idphac( 200),narc
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200),&
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200),&
                     drhlx( 200),idzhlx( 200),idrhlx( 200),&
                     alfhlx( 200),bethlx( 200),phaix1( 200),&
                     phaix2( 200),idphlx( 200),nhelx
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension bx(1),by(1),bz(1)
!
         xmin = x
         ymin = y
         zmin = z
!
         bx(1) = 0.0d0
         by(1) = 0.0d0
         bz(1) = 0.0d0
!
      if( nline .gt. 0 )  call rect(bx, by, bz)
      if( nloop .gt. 0 )  call loop(bx, by, bz)
      if(  narc .gt. 0 )  call  arc(bx, by, bz)
      if( nhelx .gt. 0 )  call helx(bx, by, bz)
!
      return
      end
      subroutine magprp(kisyu,iwrite)
!
!
!       program name   ***  mag code  ***
!               magnetic field calculation code at plasma research
!               center in the university of tsukuba.
!
!       progammed by i. katanuma       ( august  1980 ).
!
!       main program
!
!
      implicit real*8 (a-h,o-z)
!
      common /cmline/ xline(1200),yline(1200),zline(1200),cline( 200),&
                     idl12( 200),idl23( 200),nline
      common /cmloop/ xloop( 200),yloop( 200),zloop( 200),aloop( 200),&
                     cloop( 200),dzlp( 200),drlp( 200),idzlp( 200),&
                     idrlp( 200),alfa( 200),beta( 200),nloop
      common /cmarc/  xarc( 200),yarc( 200),zarc( 200),aarc( 200),&
                     carc( 200),dzarc( 200),drarc( 200),idzarc( 200),&
                     idrarc( 200),alfarc( 200),betarc( 200),&
                     phai1( 200),phai2( 200),idphac( 200),narc
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200),&
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200),&
                     drhlx( 200),idzhlx( 200),idrhlx( 200),&
                     alfhlx( 200),bethlx( 200),phaix1( 200),&
                     phaix2( 200),idphlx( 200),nhelx
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension xlined(6),ylined(6),zlined(6)
      namelist /data1/ nline,nloop,narc,nhelx,npcoil,ntcoil,nbcoil
      namelist /data2/ xlined,ylined,zlined,clined,idl12d,idl23d
      namelist /data3/ xloopd,yloopd,zloopd,aloopd,cloopd,dzlpd,drlpd,&
                      idzlpd,idrlpd,alfad,betad
      namelist /data4/ xarcd,yarcd,zarcd,aarcd,carcd,dzarcd,drarcd,&
                      idzard,idrard,alfard,betard,phai1d,phai2d,idphad
      namelist /helix/ xhelxd,yhelxd,zhelxd,ahelxd,chelxd,zdhlxd,&
                      rdhlxd,dzhlxd,drhlxd,idzhxd,idrhxd,alfhxd,bethxd,&
                      phai1d,phai2d,idphxd
      namelist /cpcil/ xpcoil,ypcoil,zpcoil,rx,ry,rz,cpcoil,apclx,apcly,&
                      alfap,betap,dxp,dyp,idxp,idyp,idphap
      namelist /ctcil/ xtcoil,ytcoil,ztcoil,rtx,rty,ctcoil,atcoil,&
                      alfat,betat,dyt,dzt,idyt,idzt,idphat
      namelist /cbcil/ xbcoil,ybcoil,zbcoil,rbz,cbcoil,abclx,abcly,&
                      alfab,betab,dxb,dyb,dphab,idxb,idyb,idphab
!*----------------------------------------------------------------------
!*  iwrite .gt. 0 :  coil data are necessary and write coil data.
!*  iwrite .eq. 0 :  coil data are necessary but do not write coil data.
!*  iwrite .lt. 0 :  coil data are read from 'mag.data'.
!*----------------------------------------------------------------------
!
         nline = 0
         nloop = 0
         narc = 0
         nhelx = 0
!
         npcoil = 0
         ntcoil = 0
         nbcoil = 0
!
         ix = 1
         iy = 1
         iz = 1
!
         dx = 1.0d0
         dy = 1.0d0
         dz = 1.0d0
!
!*----------------------------------------------------------------------
      if( iwrite .lt. 0 )  then
!        open(unit=kisyu,file='mag.data',status='old',
!    &        access='sequential',form='formatted',action='read')
!Y         open(unit=kisyu,file='magC1991.data',status='old',
!Y     &        access='sequential',form='formatted')
!!         open(unit=kisyu,file='mag091PU.data',status='old', &
!!             access='sequential',form='formatted')
         CALL FROPEN(kisyu,'mag091PU.data',1,0,'YAMA',IERR)
!Y1         WRITE(6,*) '## YAMA -MAG091PU.DATA- 20040815'
         WRITE(6,*) '## YAMA -MAG091PU.DATA- 20051016'

!Y2         open(unit=kisyu,file='magCR050PU.data',status='old',
!Y2     &        access='sequential',form='formatted')
!Y2         WRITE(6,*) '## YAMA -MAGCR050PU.DATA- 20040818'

!Y3         open(unit=kisyu,file='magCR054PU.data',status='old',
!Y3     &        access='sequential',form='formatted')
!Y3         WRITE(6,*) '## YAMA -MAGCR054PU.DATA- 20040925'
      endif
!*----------------------------------------------------------------------
!
      read(kisyu,data1)
      if( iwrite .gt. 0 )  write(6,56)
      if( iwrite .gt. 0 )  write(6,data1)
!
      if( nline .le. 0 )  go to 15
!
      do 25 i = 1,nline
!
      read(kisyu,data2)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,data2)
!
         xline(i * 6 - 5) = xlined(1)
         xline(i * 6 - 4) = xlined(2)
         xline(i * 6 - 3) = xlined(3)
         xline(i * 6 - 2) = xlined(4)
         xline(i * 6 - 1) = xlined(5)
         xline(i * 6)     = xlined(6)
!
         yline(i * 6 - 5) = ylined(1)
         yline(i * 6 - 4) = ylined(2)
         yline(i * 6 - 3) = ylined(3)
         yline(i * 6 - 2) = ylined(4)
         yline(i * 6 - 1) = ylined(5)
         yline(i * 6)     = ylined(6)
!
         zline(i * 6 - 5) = zlined(1)
         zline(i * 6 - 4) = zlined(2)
         zline(i * 6 - 3) = zlined(3)
         zline(i * 6 - 2) = zlined(4)
         zline(i * 6 - 1) = zlined(5)
         zline(i * 6)     = zlined(6)
!
         cline(i) = clined
         idl12(i) = idl12d
         idl23(i) = idl23d
   25              continue
!
   15 if( nloop .le. 0 )  go to 35
!
      do 45 i = 1,nloop
!
      read(kisyu,data3)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,data3)
!
         xloop(i) = xloopd
         yloop(i) = yloopd
         zloop(i) = zloopd
         aloop(i) = aloopd
         cloop(i) = cloopd
         dzlp(i)  = dzlpd
         drlp(i)  = drlpd
         idzlp(i) = idzlpd
         idrlp(i) = idrlpd
         alfa(i)  = alfad
         beta(i)  = betad
   45              continue
!
   35 if( narc .le. 0 )  go to 55
!
      do 65 i = 1,narc
!
      read(kisyu,data4)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,data4)
!
         xarc(i) = xarcd
         yarc(i) = yarcd
         zarc(i) = zarcd
         aarc(i) = aarcd
         carc(i) = carcd
         dzarc(i) = dzarcd
         drarc(i) = drarcd
         idzarc(i) = idzard
         idrarc(i) = idrard
         alfarc(i) = alfard
         betarc(i) = betard
         phai1(i)  = phai1d
         phai2(i)  = phai2d
         idphac(i) = idphad
   65              continue
!
   55 if( nhelx .le. 0 )  go to 75
!
      do 85 i = 1,nhelx
!
      read(kisyu,helix)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,helix)
!
         xhelx(i) = xhelxd
         yhelx(i) = yhelxd
         zhelx(i) = zhelxd
         ahelx(i) = ahelxd
         chelx(i) = chelxd
         zdhlx(i) = zdhlxd
         rdhlx(i) = rdhlxd
         dzhlx(i) = dzhlxd
         drhlx(i) = drhlxd
         idzhlx(i) = idzhxd
         idrhlx(i) = idrhxd
         alfhlx(i) = alfhxd
         bethlx(i) = bethxd
         phaix1(i) = phai1d
         phaix2(i) = phai2d
         idphlx(i) = idphxd
   85              continue
   75              continue
!
         ii = npcoil
   10 if( ii .le. 0 )  go to 20
      read(kisyu,cpcil)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,cpcil)
      call pcoil(xpcoil,ypcoil,zpcoil,rx,ry,rz,cpcoil,apclx,apcly,&
                alfap,betap,dxp,dyp,idxp,idyp,idphap)
         ii = ii - 1
              go to 10
   20              continue

         ii = ntcoil
   30 if( ii .le. 0 )  go to 40
      read(kisyu,ctcil)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,ctcil)
      call tcoil(xtcoil,ytcoil,ztcoil,rtx,rty,ctcoil,atcoil,&
                alfat,betat,dyt,dzt,idyt,idzt,idphat)
         ii = ii - 1
              go to 30
   40              continue
!
         ii = nbcoil
   50 if( ii .le. 0 )  go to 60
      read(kisyu,cbcil)
      if( iwrite .gt. 0 )  write(6,66)
      if( iwrite .gt. 0 )  write(6,cbcil)
      call bcoil(xbcoil,ybcoil,zbcoil,rbz,cbcoil,abclx,abcly,&
                alfab,betab,dxb,dyb,idxb,idyb,dphab,idphab)
         ii = ii - 1
              go to 50
   60              continue

!
!*----------------------------------------------------------------------
      if( iwrite .lt. 0 )  then
         close(unit=kisyu,status='keep')
      endif
!*----------------------------------------------------------------------
!
!     if( nline .gt. 0 )&
!        write(6,6) (i,(xline(i*6+j-6),j=1,6),(yline(i*6+j-6),j=1,6),&
!                    (zline(i*6+j-6),j=1,6),cline(i),idl12(i),&
!                    idl23(i),i=1,nline)
!     if( nloop .gt. 0 )&
!        write(6,16) (i,xloop(i),yloop(i),zloop(i),aloop(i),cloop(i),&
!                     dzlp(i),drlp(i),idzlp(i),idrlp(i),alfa(i),&
!                     beta(i),i=1,nloop)
!     if(  narc .gt. 0 )&
!        write(6,26) (i,xarc(i),yarc(i),zarc(i),aarc(i),carc(i),&
!                     dzarc(i),drarc(i),idzarc(i),idrarc(i),alfarc(i),&
!                     betarc(i),phai1(i),phai2(i),idphac(i),i=1,narc)
!
         c = 10.0d0
         pi = 3.141592653589793d0
         dpi = pi / 180.0d0
!
      if( nline .le. 0 )  go to 100
!
      do 110 i = 1,nline
         cline(i) = cline(i)/( dfloat(idl12(i) * idl23(i)) * c )
  110          continue
!
  100 if( nloop .le. 0 )  go to 200
!
      do 210 i = 1,nloop
          cloop(i) = cloop(i)/( dfloat(idrlp(i) * idzlp(i)) * c )
          cloop(i) = 2.0d0 * cloop(i)
          alfa(i) = alfa(i) * dpi
          beta(i) = beta(i) * dpi
  210          continue
!
  200 if( narc .le. 0 )  go to 300
!
      do 310 i = 1,narc
          carc(i) = carc(i)/( dfloat(idrarc(i) * idzarc(i)) * c )
          alfarc(i) = alfarc(i) * dpi
          betarc(i) = betarc(i) * dpi
          phai1(i) = phai1(i) * dpi
          phai2(i) = phai2(i) * dpi
  310          continue
!
  300 if( nhelx .le. 0 )  go to 350
!
      do 360 i = 1,nhelx
          chelx(i) = chelx(i)/( dfloat(idrhlx(i) * idzhlx(i)) * c )
          alfhlx(i) = alfhlx(i) * dpi
          bethlx(i) = bethlx(i) * dpi
          phaix1(i) = phaix1(i) * dpi
          phaix2(i) = phaix2(i) * dpi
  360          continue
!
  350          continue
   56 format(1h1/)
   66 format(1h )
!
      return
      end
      subroutine line(bx,by,bz,xx1,yy1,zz1,xx2,yy2,zz2,cur)
      implicit real*8 (a-h,o-z)
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension bx(ix, iy, iz),by(ix, iy, iz),bz(ix, iy, iz)
!
      do 100 k = 1,iz
         zz = zmin + dfloat(k - 1) * dz
         z1 = zz1 - zz
         z2 = zz2 - zz
      do 100 j = 1,iy
         yy = ymin + dfloat(j - 1) * dy
         y1 = yy1 - yy
         y2 = yy2 - yy
         yz = y1 * z2 - y2 * z1
      do 100 i = 1,ix
         xx = xmin + dfloat(i - 1) * dx
         x1 = xx1 - xx
         x2 = xx2 - xx
         zx = z1 * x2 - z2 * x1
         xy = x1 * y2 - x2 * y1
         row1 = x1 * x1 + y1 * y1 + z1 * z1
         row2 = x2 * x2 + y2 * y2 + z2 * z2
         ro1 = dsqrt( row1 )
         ro2 = dsqrt( row2 )
         rod = x1 * x2 + y1 * y2 + z1 * z2
         ro12 = ro1 * ro2
         rod2 = ro12 + rod
      ro12 = dmax1( ro12, 1.0d-8 )
      rod2 = dmax1( rod2, ro12 * 1.0d-8 )
         g = (ro1 + ro2) / ( ro12 * rod2 )
         curr = cur * g
         bx(i, j, k) = bx(i, j, k) + curr * yz
         by(i, j, k) = by(i, j, k) + curr * zx
         bz(i, j, k) = bz(i, j, k) + curr * xy
  100              continue
!
      return
      end
      subroutine rect(bx, by, bz)
      implicit real*8 (a-h,o-z)
      common /cmline/ xline(1200),yline(1200),zline(1200),cline( 200),&
                     idl12( 200),idl23( 200),nline
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension bx(ix, iy, iz),by(ix, iy, iz),bz(ix, iy, iz)
!
!*  in this subroutine  magnetic fields due to straight coils are
!*  calculated.
!
!
         iter = 0
  100    iter = iter + 1
!
         iter6 = iter * 6
!
         x11 = xline(iter6 - 5)
         x12 = xline(iter6 - 4)
         x13 = xline(iter6 - 3)
         x21 = xline(iter6 - 2)
         x22 = xline(iter6 - 1)
         x23 = xline(iter6)
!
         y11 = yline(iter6 - 5)
         y12 = yline(iter6 - 4)
         y13 = yline(iter6 - 3)
         y21 = yline(iter6 - 2)
         y22 = yline(iter6 - 1)
         y23 = yline(iter6)
!
         z11 = zline(iter6 - 5)
         z12 = zline(iter6 - 4)
         z13 = zline(iter6 - 3)
         z21 = zline(iter6 - 2)
         z22 = zline(iter6 - 1)
         z23 = zline(iter6)
!
!     if( idl12(iter).ne.1 .or. idl23(iter).ne.1 )  go to 200
!
!        x12 = x11
!        x13 = x11
!        x22 = x21
!        x23 = x21
!
!        y12 = y11
!        y13 = y11
!        y22 = y21
!        y23 = y21
!
!        z12 = z11
!        z13 = z11
!        z22 = z21
!        z23 = z21
!
! 200              continue
!
         adl1 = 1.0d0 / dfloat( idl23(iter) )
         dx13 = (x13 - x12) * adl1
         dx23 = (x23 - x22) * adl1
         dy13 = (y13 - y12) * adl1
         dy23 = (y23 - y22) * adl1
         dz13 = (z13 - z12) * adl1
         dz23 = (z23 - z22) * adl1
!
         adl2 = 1.0d0 / dfloat( idl12(iter) )
         dx12 = (x12 - x11) * adl2
         dx22 = (x22 - x21) * adl2
         dy12 = (y12 - y11) * adl2
         dy22 = (y22 - y21) * adl2
         dz12 = (z12 - z11) * adl2
         dz22 = (z22 - z21) * adl2
!
         it23 = 0
  300    it23 = it23 + 1
!
         ait23 = dfloat( it23 ) - 0.5d0
!
         xx11 = x11 + ait23 * dx13
         xx21 = x21 + ait23 * dx23
         yy11 = y11 + ait23 * dy13
         yy21 = y21 + ait23 * dy23
         zz11 = z11 + ait23 * dz13
         zz21 = z21 + ait23 * dz23
!
         it12 = 0
  400    it12 = it12 + 1
!
         ait12 = dfloat( it12 ) - 0.5d0
!
         xx1 = xx11 + ait12 * dx12
         xx2 = xx21 + ait12 * dx22
         yy1 = yy11 + ait12 * dy12
         yy2 = yy21 + ait12 * dy22
         zz1 = zz11 + ait12 * dz12
         zz2 = zz21 + ait12 * dz22
!
      call line(bx,by,bz,xx1,yy1,zz1,xx2,yy2,zz2,cline(iter))
!
      if( it12 .lt. idl12(iter) )  go to 400
      if( it23 .lt. idl23(iter) )  go to 300
!
      if( iter .lt. nline )  go to 100
!
      return
      end
      subroutine loop(bx, by, bz)
      implicit real*8 (a-h,o-z)
      common /cmloop/ xloop( 200),yloop( 200),zloop( 200),aloop( 200),&
                     cloop( 200),dzlp( 200),drlp( 200),idzlp( 200),&
                     idrlp( 200),alfa( 200),beta( 200),nloop
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension bx(ix, iy, iz),by(ix, iy, iz),bz(ix, iy, iz)
!
      dimension eic1a(5),eic1b(5),eic2a(5),eic2b(5)
      data eic1a / 0.01451196212d0, 0.03742563713d0, 0.03590092383d0,&
                  0.09666344259d0, 1.38629436112d0 /
      data eic1b / 0.00441787012d0, 0.03328355346d0, 0.06880248576d0,&
                  0.12498593597d0, 0.50000000000d0 /
      data eic2a / 0.01736506451d0, 0.04757383546d0, 0.06260601220d0,&
                  0.44325141463d0, 1.00000000000d0 /
      data eic2b / 0.00526449639d0, 0.04069697526d0, 0.09200180037d0,&
                  0.24998368310d0, 0.00000000000d0 /
!
!*  in this subroutine  magnetic fields due to loop coils are
!*  calculated.
!
!
         iter = 0
  100    iter = iter + 1
!
         cosa = dcos( alfa(iter) )
         sina = dsin( alfa(iter) )
         cosb = dcos( beta(iter) )
         sinb = dsin( beta(iter) )
!
         adr = drlp(iter) / dfloat( idrlp(iter) )
         adz = dzlp(iter) / dfloat( idzlp(iter) )
!
         itz = 0
  200    itz = itz + 1
!
         dd = (dfloat( itz ) - 0.5d0) * adz - 0.5d0 * dzlp(iter)
         xx1 = xloop(iter) + dd * sina * sinb
         yy1 = yloop(iter) - dd * cosa * sinb
         zz1 = zloop(iter) + dd * cosb
!
         itr = 0
  300    itr = itr + 1
!
         ada = aloop(iter) + (dfloat( itr ) - 0.5d0) * adr&
              - 0.5d0 * drlp(iter)
         ada2 = ada * ada
!
      do 400 k = 1,iz
         zz = zmin + dfloat(k - 1) * dz
         z0 = zz - zz1
         zcosb = z0 * cosb
         zsinb = z0 * sinb
      do 400 j = 1,iy
         yy = ymin + dfloat(j - 1) * dy
         y0 = yy - yy1
         ysina = y0 * sina
         ycosa = y0 * cosa
      do 400 i = 1,ix
         xx = xmin + dfloat(i - 1) * dx
         x0 = xx - xx1
         z1 = x0 * sina - ycosa
         x1 = x0 * cosa + ysina
         y1 = zsinb - z1 * cosb
         z1 = zcosb + z1 * sinb
         ro = dsqrt( x1 * x1 + y1 * y1 )
      ro = dmax1( ro , 1.0d-4 )
         aro = 1.0d0 / ro
         rz = ro * ro + z1 * z1
         ar = 2.0d0 * ada * ro
         arz = ada2 + rz
         eta1 = 1.0d0 / ( arz + ar )
         drzr = arz - ar
      if( drzr .lt. 1.0d-8)  drzr = 1.0d-8
!*  start of calculating complete elliptik integral
         eta = drzr * eta1
         etaln = dlog( eta )
!*  first complete elliptik integral
         sum1a = eic1a(1)
         sum1b = eic1b(1)
      do 10 ii = 2,5
         sum1a = sum1a * eta + eic1a(ii)
         sum1b = sum1b * eta + eic1b(ii)
   10              continue
!*  second complete elliptik integral
         sum2a = eic2a(1)
         sum2b = eic2b(1)
      do 20 ii = 2,5
         sum2a = sum2a * eta + eic2a(ii)
         sum2b = sum2b * eta + eic2b(ii)
   20              continue
         fistk = sum1a - sum1b * etaln
         scndk = sum2a - sum2b * etaln
!*  end of calculating complete elliptik integral
         arzm = 1.0d0 / drzr
         arzb = dsqrt( eta1 )
      if( ro .lt. 1.0d-3 )  go to 410
         br = cloop(iter) * z1 * arzb&
             * ( arz * arzm * scndk - fistk ) * aro
              go to 420
  410    br = 0.0d0
  420          continue
         bzz = cloop(iter) * arzb * ( ( ada2 - rz )&
              * arzm * scndk + fistk )
         bxx = br * x1 * aro
         byy = br * y1 * aro
         bb1 = byy * cosb - bzz * sinb
         bx(i, j, k) = bx(i, j, k) + bxx * cosa - bb1 * sina
         by(i, j, k) = by(i, j, k) + bxx * sina + bb1 * cosa
         bz(i, j, k) = bz(i, j, k) + byy * sinb + bzz * cosb
  400         continue
!
      if( itr .lt. idrlp(iter) )  go to 300
      if( itz .lt. idzlp(iter) )  go to 200
!
      if( iter .lt. nloop )  go to 100
!
      return
      end
      subroutine arc(bx, by, bz)
      implicit real*8 (a-h,o-z)
      common /cmarc/  xarc( 200),yarc( 200),zarc( 200),aarc( 200),&
                     carc( 200),dzarc( 200),drarc( 200),idzarc( 200),&
                     idrarc( 200),alfarc( 200),betarc( 200),&
                     phai1( 200),phai2( 200),idphac( 200),narc
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200),&
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200),&
                     drhlx( 200),idzhlx( 200),idrhlx( 200),&
                     alfhlx( 200),bethlx( 200),phaix1( 200),&
                     phaix2( 200),idphlx( 200),nhelx
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension bx(ix, iy, iz),by(ix, iy, iz),bz(ix, iy, iz)
!
!*  in this subroutine  magnetic fields due to arc coils are
!*  calculated.
!
!
         iter = 0
  100    iter = iter + 1
!
         cosa = dcos( alfarc(iter) )
         sina = dsin( alfarc(iter) )
         cosb = dcos( betarc(iter) )
         sinb = dsin( betarc(iter) )
!
         dphai = ( phai2(iter) - phai1(iter) )&
                / dfloat( idphac(iter) )
         adr = drarc(iter) / dfloat( idrarc(iter) )
         adz = dzarc(iter) / dfloat( idzarc(iter) )
!
         itz = 0
  200    itz = itz + 1
!
         dd = (dfloat( itz ) - 0.5d0) * adz - 0.5d0 * dzarc(iter)
         x1 = xarc(iter) + dd * sina * sinb
         y1 = yarc(iter) - dd * cosa * sinb
         z1 = zarc(iter) + dd * cosb
!
         itr = 0
  300    itr = itr + 1
!
         ada = aarc(iter) + (dfloat( itr ) - 0.5d0) * adr&
              - 0.5d0 * drarc(iter)
!
         xx = ada * dcos( phai1(iter) )
         yy = ada * dsin( phai1(iter) )
         xx1 = x1 + xx * cosa - yy * sina * cosb
         yy1 = y1 + xx * sina + yy * cosa * cosb
         zz1 = z1 + yy * sinb
!
         itph = 0
  400    itph = itph + 1
!
         phi = phai1(iter) + dfloat( itph ) * dphai
         xx = ada * dcos( phi )
         yy = ada * dsin( phi )
         xx2 = x1 + xx * cosa - yy * sina * cosb
         yy2 = y1 + xx * sina + yy * cosa * cosb
         zz2 = z1 + yy * sinb
!
      call line(bx,by,bz,xx1,yy1,zz1,xx2,yy2,zz2,carc(iter))
!
         xx1 = xx2
         yy1 = yy2
         zz1 = zz2
!
      if( itph .lt. idphac(iter) )  go to 400
      if( itr  .lt. idrarc(iter) )  go to 300
      if( itz  .lt. idzarc(iter) )  go to 200
!
      if( iter .lt. narc )  go to 100
!
      return
      end
      subroutine helx(bx, by, bz)
      implicit real*8 (a-h,o-z)
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200),&
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200),&
                     drhlx( 200),idzhlx( 200),idrhlx( 200),&
                     alfhlx( 200),bethlx( 200),phaix1( 200),&
                     phaix2( 200),idphlx( 200),nhelx
      common /region/ xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,ix,iy,iz
      dimension bx(ix, iy, iz),by(ix, iy, iz),bz(ix, iy, iz)
!
!*  in this subroutine  magnetic fields due to helical coils are
!*  calculated.
!
!
         iter = 0
  100    iter = iter + 1
!
         cosa = dcos( alfhlx(iter) )
         sina = dsin( alfhlx(iter) )
         cosb = dcos( bethlx(iter) )
         sinb = dsin( bethlx(iter) )
!
         dphai = ( phaix2(iter) - phaix1(iter) )&
                / dfloat( idphlx(iter) )
         zdz = zdhlx(iter) / dfloat( idphlx(iter) )
         rdr = rdhlx(iter) / dfloat( idphlx(iter) )
         adr = drhlx(iter) / dfloat( idrhlx(iter) )
         adz = dzhlx(iter) / dfloat( idzhlx(iter) )
!
         itz = 0
  200    itz = itz + 1
!
         dd = (dfloat( itz ) - 0.5d0) * adz - 0.5d0 * dzhlx(iter)
         x1 = xhelx(iter) + dd * sina * sinb
         y1 = yhelx(iter) - dd * cosa * sinb
         z1 = zhelx(iter) + dd * cosb
!
         itr = 0
  300    itr = itr + 1
!
         ada = ahelx(iter) + (dfloat( itr ) - 0.5d0) * adr&
              - 0.5d0 * drhlx(iter)
!
         xx = ada * dcos( phaix1(iter) )
         yy = ada * dsin( phaix1(iter) )
         xx1 = x1 + xx * cosa - yy * sina * cosb
         yy1 = y1 + xx * sina + yy * cosa * cosb
         zz1 = z1 + yy * sinb
!
         itph = 0
  400    itph = itph + 1
!
         phi = phaix1(iter) + dfloat( itph ) * dphai
         ada1 = ada + dfloat( itph ) * rdr
         dzz = dfloat( itph ) * zdz
         xx = ada1 * dcos( phi )
         yy = ada1 * dsin( phi )
         xx2 = x1 + xx * cosa - yy * sina * cosb + dzz * sina * sinb
         yy2 = y1 + xx * sina + yy * cosa * cosb - dzz * cosa * sinb
         zz2 = z1 + yy * sinb                    + dzz * cosb
!
      call line(bx,by,bz,xx1,yy1,zz1,xx2,yy2,zz2,chelx(iter))
!
         xx1 = xx2
         yy1 = yy2
         zz1 = zz2
!
      if( itph .lt. idphlx(iter) )  go to 400
      if( itr  .lt. idrhlx(iter) )  go to 300
      if( itz  .lt. idzhlx(iter) )  go to 200
!
      if( iter .lt. nhelx )  go to 100
!
      return
      end
      subroutine pcoil(xpcoil,ypcoil,zpcoil,rx,ry,rz,cpcoil,apclx,apcly,&
                      alfap,betap,dxp,dyp,idxp,idyp,idphap)
      implicit real*8 (a-h,o-z)
      common /cmline/ xline(1200),yline(1200),zline(1200),cline( 200),&
                     idl12( 200),idl23( 200),nline
      common /cmarc/  xarc( 200),yarc( 200),zarc( 200),aarc( 200),&
                     carc( 200),dzarc( 200),drarc( 200),idzarc( 200),&
                     idrarc( 200),alfarc( 200),betarc( 200),&
                     phai1( 200),phai2( 200),idphac( 200),narc
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200),&
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200),&
                     drhlx( 200),idzhlx( 200),idrhlx( 200),&
                     alfhlx( 200),bethlx( 200),phaix1( 200),&
                     phaix2( 200),idphlx( 200),nhelx
      dimension xx(6,8),yy(6,8),zz(6,8),rcx(8),rcy(8),rcz(8)
!
         pi = 3.141592653589793d0
         dpi = pi / 180.0d0
!
         xx( 1, 1) = 0.5d0 * rx + apcly - 0.5d0 * dxp
         yy( 1, 1) =-0.5d0 * ry - apclx + 0.5d0 * dyp
         zz( 1, 1) =-0.5d0 * rz
         xx( 2, 1) = xx( 1, 1) + dxp
         yy( 2, 1) = yy( 1, 1)
         zz( 2, 1) = zz( 1, 1)
         xx( 3, 1) = xx( 2, 1)
         yy( 3, 1) = yy( 2, 1) - dyp
         zz( 3, 1) = zz( 1, 1)
!
         xx( 4, 1) = xx( 1, 1)
         yy( 4, 1) = yy( 1, 1)
         zz( 4, 1) = 0.5d0 * rz
         xx( 5, 1) = xx( 2, 1)
         yy( 5, 1) = yy( 2, 1)
         zz( 5, 1) = zz( 4, 1)
         xx( 6, 1) = xx( 3, 1)
         yy( 6, 1) = yy( 3, 1)
         zz( 6, 1) = zz( 4, 1)
!
      do 100 i = 1,6
         xx( i, 2) = xx( i, 1)
         yy( i, 2) = yy( i, 1) + ry + 2.0d0 * apclx
         zz( i, 2) = zz( i, 1)
         xx( i, 3) = xx( i, 1) - rx - 2.0d0 * apcly
         yy( i, 3) = yy( i, 2)
         zz( i, 3) = zz( i, 1)
         xx( i, 4) = xx( i, 3)
         yy( i, 4) = yy( i, 1)
         zz( i, 4) = zz( i, 1)
  100              continue
!
         xx( 1, 5) = 0.5d0 * rx
         yy( 1, 5) =-0.5d0 * ry - apclx + 0.5d0 * dyp
         zz( 1, 5) = 0.5d0 * rz + apcly - 0.5d0 * dxp
         xx( 2, 5) = xx( 1, 5)
         yy( 2, 5) = yy( 1, 5) - dyp
         zz( 2, 5) = zz( 1, 5)
         xx( 3, 5) = xx( 1, 5)
         yy( 3, 5) = yy( 2, 5)
         zz( 3, 5) = zz( 1, 5) + dxp
!
         xx( 4, 5) =-0.5d0 * rx
         yy( 4, 5) = yy( 1, 5)
         zz( 4, 5) = zz( 1, 5)
         xx( 5, 5) = xx( 4, 5)
         yy( 5, 5) = yy( 2, 5)
         zz( 5, 5) = zz( 2, 5)
         xx( 6, 5) = xx( 4, 5)
         yy( 6, 5) = yy( 3, 5)
         zz( 6, 5) = zz( 3, 5)
!
      do 200 i = 1,6
         xx( i, 6) = xx( i, 5)
         yy( i, 6) = yy( i, 5) + ry + 2.0d0 * apclx
         zz( i, 6) = zz( i, 5)
  200              continue
!
         xx( 1, 7) = 0.5d0 * rx + apcly - 0.5d0 * dxp
         yy( 1, 7) =-0.5d0 * ry
         zz( 1, 7) =-0.5d0 * rz - apclx - 0.5d0 * dyp
         xx( 2, 7) = xx( 1, 7) + dxp
         yy( 2, 7) = yy( 1, 7)
         zz( 2, 7) = zz( 1, 7)
         xx( 3, 7) = xx( 2, 7)
         yy( 3, 7) = yy( 1, 7)
         zz( 3, 7) = zz( 2, 7) + dyp
!
         xx( 4, 7) = xx( 1, 7)
         yy( 4, 7) = 0.5d0 * ry
         zz( 4, 7) = zz( 1, 7)
         xx( 5, 7) = xx( 2, 7)
         yy( 5, 7) = yy( 4, 7)
         zz( 5, 7) = zz( 2, 7)
         xx( 6, 7) = xx( 3, 7)
         yy( 6, 7) = yy( 4, 7)
         zz( 6, 7) = zz( 3, 7)
!
      do 300 i = 1,6
         xx( i, 8) = xx( i, 7) - rx - 2.0d0 * apcly
         yy( i, 8) = yy( i, 7)
         zz( i, 8) = zz( i, 7)
  300              continue
!
         rcx( 1) = 0.5d0 * rx + apcly
         rcy( 1) =-0.5d0 * ry
         rcz( 1) =-0.5d0 * rz
         rcx( 2) = rcx( 1)
         rcy( 2) = 0.5d0 * ry
         rcz( 2) = rcz( 1)
         rcx( 3) =-0.5d0 * rx - apcly
         rcy( 3) = rcy( 1)
         rcz( 3) = rcz( 1)
         rcx( 4) = rcx( 3)
         rcy( 4) = rcy( 2)
         rcz( 4) = rcz( 2)
!
         rcx( 5) = 0.5d0 * rx
         rcy( 5) =-0.5d0 * ry - apclx
         rcz( 5) = 0.5d0 * rz
         rcx( 6) = rcx( 5)
         rcy( 6) = 0.5d0 * ry + apclx
         rcz( 6) = rcz( 5)
         rcx( 7) =-0.5d0 * rx
         rcy( 7) = rcy( 5)
         rcz( 7) = rcz( 5)
         rcx( 8) = rcx( 7)
         rcy( 8) = rcy( 6)
         rcz( 8) = rcz( 6)
!
         cosa = dcos( alfap * dpi )
         sina = dsin( alfap * dpi )
         cosb = dcos( betap * dpi )
         sinb = dsin( betap * dpi )
!
      do 400 j = 1,8
         nline = nline + 1
      do 400 i = 1,6
         x1 = xx( i, j) * cosa - yy( i, j) * sina * cosb&
                              + zz( i, j) * sina * sinb
         y1 = xx( i, j) * sina + yy( i, j) * cosa * cosb&
                              - zz( i, j) * cosa * sinb
         z1 = yy( i, j) * sinb + zz( i, j) * cosb
         xline(6 * nline + i - 6) = x1 + xpcoil
         yline(6 * nline + i - 6) = y1 + ypcoil
         zline(6 * nline + i - 6) = z1 + zpcoil
  400              continue
!
      do 500 i = 1,4
         x1 = rcx( i) * cosa - rcy( i) * sina * cosb&
                            + rcz( i) * sina * sinb
         y1 = rcx( i) * sina + rcy( i) * cosa * cosb&
                            - rcz( i) * cosa * sinb
         z1 = rcy( i) * sinb + rcz( i) * cosb
         narc = narc + 1
         xarc(narc) = x1 + xpcoil
         yarc(narc) = y1 + ypcoil
         zarc(narc) = z1 + zpcoil
         aarc(narc) = apclx
         idphac(narc) = idphap
         dzarc(narc) = dxp
         drarc(narc) = dyp
         idzarc(narc) = idxp
         idrarc(narc) = idyp
         alfarc(narc) = 90.0d0 + alfap
         betarc(narc) = 90.0d0
  500              continue
!
      do 600 i = 5,8
         x1 = rcx( i) * cosa - rcy( i) * sina * cosb&
                            + rcz( i) * sina * sinb
         y1 = rcx( i) * sina + rcy( i) * cosa * cosb&
                            - rcz( i) * cosa * sinb
         z1 = rcy( i) * sinb + rcz( i) * cosb
         narc = narc + 1
         xarc(narc) = x1 + xpcoil
         yarc(narc) = y1 + ypcoil
         zarc(narc) = z1 + zpcoil
         aarc(narc) = apcly
         idphac(narc) = idphap
         dzarc(narc) = dyp
         drarc(narc) = dxp
         idzarc(narc) = idyp
         idrarc(narc) = idxp
         alfarc(narc) = 180.0d0 + alfap
         betarc(narc) = 90.0d0 - betap
  600              continue
!
         cline(nline - 7) = - cpcoil
         idl12(nline - 7) = idxp
         idl23(nline - 7) = idyp
!
         cline(nline - 6) = cpcoil
         idl12(nline - 6) = idxp
         idl23(nline - 6) = idyp
!
         cline(nline - 5) = - cpcoil
         idl12(nline - 5) = idxp
         idl23(nline - 5) = idyp
!
         cline(nline - 4) = cpcoil
         idl12(nline - 4) = idxp
         idl23(nline - 4) = idyp
!
         cline(nline - 3) = - cpcoil
         idl12(nline - 3) = idyp
         idl23(nline - 3) = idxp
!
         cline(nline - 2) = cpcoil
         idl12(nline - 2) = idyp
         idl23(nline - 2) = idxp
!
         cline(nline - 1) = cpcoil
         idl12(nline - 1) = idxp
         idl23(nline - 1) = idyp
!
         cline(nline) = - cpcoil
         idl12(nline) = idxp
         idl23(nline) = idyp
!
         carc(narc - 7) = cpcoil
         phai1(narc - 7) = - 180.0d0 + betap
         phai2(narc - 7) = - 90.0d0 + betap
!
         carc(narc - 6) = cpcoil
         phai1(narc - 6) = - 90.0d0 + betap
         phai2(narc - 6) = betap
!
         carc(narc - 5) = - cpcoil
         phai1(narc - 5) = - 180.0d0 + betap
         phai2(narc - 5) = - 90.0d0 + betap
!
         carc(narc - 4) = - cpcoil
         phai1(narc - 4) = - 90.0d0 + betap
         phai2(narc - 4) = betap
!
         carc(narc - 3) = cpcoil
         phai1(narc - 3) = 90.0d0
         phai2(narc - 3) = 180.0d0
!
         carc(narc - 2) = - cpcoil
         phai1(narc - 2) = 90.0d0
         phai2(narc - 2) = 180.0d0
!
         carc(narc - 1) = cpcoil
         phai1(narc - 1) = 0.0d0
         phai2(narc - 1) = 90.0d0
!
         carc(narc) = - cpcoil
         phai1(narc) = 0.0d0
         phai2(narc) = 90.0d0
!
      return
      end
      subroutine tcoil(xtcoil,ytcoil,ztcoil,rtx,rty,ctcoil,atcoil,&
                      alfat,betat,dyt,dzt,idyt,idzt,idphat)
      implicit real * 8 (a-h,o-z)
      common /cmline/ xline(1200),yline(1200),zline(1200),cline( 200),&
                     idl12( 200),idl23( 200),nline
      common /cmarc/  xarc( 200),yarc( 200),zarc( 200),aarc( 200),&
                     carc( 200),dzarc( 200),drarc( 200),idzarc( 200),&
                     idrarc( 200),alfarc( 200),betarc( 200),&
                     phai1( 200),phai2( 200),idphac( 200),narc
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200),&
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200),&
                     drhlx( 200),idzhlx( 200),idrhlx( 200),&
                     alfhlx( 200),bethlx( 200),phaix1( 200),&
                     phaix2( 200),idphlx( 200),nhelx
      dimension xx(6,4),yy(6,4),zz(6,4),rcx(4),rcy(4),rcz(4)
!
         pi = 3.141592653589793d0
         dpi = pi / 180.0d0
!
         xx( 1, 1) = 0.5d0 * rtx + atcoil - 0.5d0 * dyt
         yy( 1, 1) =-0.5d0 * rty
         zz( 1, 1) =-0.5d0 * dzt
         xx( 2, 1) = xx( 1, 1) + dyt
         yy( 2, 1) = yy( 1, 1)
         zz( 2, 1) = zz( 1, 1)
         xx( 3, 1) = xx( 2, 1)
         yy( 3, 1) = yy( 1, 1)
         zz( 3, 1) = 0.5d0 * dzt
!
         xx( 4, 1) = xx( 1, 1)
         yy( 4, 1) = 0.5d0 * rty
         zz( 4, 1) = zz( 1, 1)
         xx( 5, 1) = xx( 2, 1)
         yy( 5, 1) = yy( 4, 1)
         zz( 5, 1) = zz( 2, 1)
         xx( 6, 1) = xx( 3, 1)
         yy( 6, 1) = yy( 4, 1)
         zz( 6, 1) = zz( 3, 1)
!
      do 100 i = 1,6
         xx( i, 2) = xx( i, 1) - rtx - 2.0d0 * atcoil
         yy( i, 2) = yy( i, 1)
         zz( i, 2) = zz( i, 1)
  100              continue
!
         xx( 1, 3) = 0.5d0 * rtx
         yy( 1, 3) = 0.5d0 * rty + atcoil + 0.5d0 * dyt
         zz( 1, 3) =-0.5d0 * dzt
         xx( 2, 3) = xx( 1, 3)
         yy( 2, 3) = yy( 1, 3) - dyt
         zz( 2, 3) = zz( 1, 3)
         xx( 3, 3) = xx( 1, 3)
         yy( 3, 3) = yy( 2, 3)
         zz( 3, 3) = zz( 2, 3) + dzt
!
         xx( 4, 3) =-0.5d0 * rtx
         yy( 4, 3) = yy( 1, 3)
         zz( 4, 3) = zz( 1, 3)
         xx( 5, 3) = xx( 4, 3)
         yy( 5, 3) = yy( 2, 3)
         zz( 5, 3) = zz( 2, 3)
         xx( 6, 3) = xx( 4, 3)
         yy( 6, 3) = yy( 3, 3)
         zz( 6, 3) = zz( 3, 3)
!
      do 200 i = 1,6
         xx( i, 4) = xx( i, 3)
         yy( i, 4) = yy( i, 3) - rty - 2.0d0 * atcoil
         zz( i, 4) = zz( i, 3)
  200              continue
!
         rcx( 1) = 0.5d0 * rtx
         rcy( 1) = 0.5d0 * rty
         rcz( 1) = 0.0d0
         rcx( 2) =-0.5d0 * rtx
         rcy( 2) = rcy( 1)
         rcz( 2) = rcz( 1)
         rcx( 3) = rcx( 2)
         rcy( 3) =-0.5d0 * rty
         rcz( 3) = rcz( 1)
         rcx( 4) = rcx( 1)
         rcy( 4) = rcy( 3)
         rcz( 4) = rcz( 1)
!
         cosa = dcos( alfat * dpi )
         sina = dsin( alfat * dpi )
         cosb = dcos( betat * dpi )
         sinb = dsin( betat * dpi )
!
      do 300 j = 1,4
         nline = nline + 1
         idl12(nline) = idyt
         idl23(nline) = idzt
      do 300 i = 1,6
         x1 = xx( i, j) * cosa - yy( i, j) * sina * cosb &
                              + zz( i, j) * sina * sinb
         y1 = xx( i, j) * sina + yy( i, j) * cosa * cosb &
                              - zz( i, j) * cosa * sinb
         z1 = yy( i, j) * sinb + zz( i, j) * cosb
         xline(6 * nline + i - 6) = x1 + xtcoil
         yline(6 * nline + i - 6) = y1 + ytcoil
         zline(6 * nline + i - 6) = z1 + ztcoil
  300              continue
!
      do 400 i = 1,4
         x1 = rcx( i) * cosa - rcy( i) * sina * cosb &
                            + rcz( i) * sina * sinb
         y1 = rcx( i) * sina + rcy( i) * cosa * cosb &
                            - rcz( i) * cosa * sinb
         z1 = rcy( i) * sinb + rcz( i) * cosb
         narc = narc + 1
         xarc(narc) = x1 + xtcoil
         yarc(narc) = y1 + ytcoil
         zarc(narc) = z1 + ztcoil
         aarc(narc) = atcoil
         carc(narc) = ctcoil
         idphac(narc) = idphat
         dzarc(narc) = dzt
         drarc(narc) = dyt
         idzarc(narc) = idzt
         idrarc(narc) = idyt
         alfarc(narc) = alfat
         betarc(narc) = betat
  400              continue
!
         cline(nline - 3) = ctcoil
         cline(nline - 2) = - ctcoil
         cline(nline - 1) = ctcoil
         cline(nline) = - ctcoil
!
         phai1(narc - 3) = 0.0d0
         phai2(narc - 3) = 90.0d0
!
         phai1(narc - 2) = 90.0d0
         phai2(narc - 2) = 180.0d0
!
         phai1(narc - 1) = 180.0d0
         phai2(narc - 1) = 270.0d0
!
         phai1(narc) = 270.0d0
         phai2(narc) = 360.0d0
!
      return
      end
      subroutine bcoil(xbcoil,ybcoil,zbcoil,rbz,cbcoil,abclx,abcly, &
                      alfab,betab,dxb,dyb,idxb,idyb,dphab,idphab)
      implicit real * 8 (a-h,o-z)
      common /cmline/ xline(1200),yline(1200),zline(1200),cline( 200), &
                     idl12( 200),idl23( 200),nline
      common /cmarc/  xarc( 200),yarc( 200),zarc( 200),aarc( 200), &
                     carc( 200),dzarc( 200),drarc( 200),idzarc( 200), &
                     idrarc( 200),alfarc( 200),betarc( 200), &
                     phai1( 200),phai2( 200),idphac( 200),narc
      common /cmhelx/ xhelx( 200),yhelx( 200),zhelx( 200),ahelx( 200), &
                     chelx( 200),zdhlx( 200),rdhlx( 200),dzhlx( 200), &
                     drhlx( 200),idzhlx( 200),idrhlx( 200), &
                     alfhlx( 200),bethlx( 200),phaix1( 200), &
                     phaix2( 200),idphlx( 200),nhelx
      dimension xx(6,4),yy(6,4),zz(6,4),rcx(4),rcy(4),rcz(4)
!
         pi = 3.141592653589793d0
         dpi = pi / 180.0d0
         tp = 0.5d0 * ( dphab * dpi - pi )
!
         idphyb = idphab
!
         cost = dcos( tp )
         sint = dsin( tp )
!
         xx( 1, 1) = abcly + 0.5d0 * dxb
         yy( 1, 1) = - abclx * cost - 0.5d0 * dyb * cost
         zz( 1, 1) = - 0.5d0 * rbz * cost + 0.5d0 * dyb * sint
         xx( 2, 1) = xx( 1, 1)
         yy( 2, 1) = yy( 1, 1) + dyb * cost
         zz( 2, 1) = zz( 1, 1) - dyb * sint
         xx( 3, 1) = xx( 2, 1) - dxb
         yy( 3, 1) = yy( 2, 1)
         zz( 3, 1) = zz( 2, 1)
!
         xx( 4, 1) = xx( 1, 1)
         yy( 4, 1) = - abclx * cost + rbz * sint - 0.5d0 * dyb * cost
         zz( 4, 1) = 0.5d0 * rbz * cost + 0.5d0 * dyb * sint
         xx( 5, 1) = xx( 4, 1)
         yy( 5, 1) = yy( 4, 1) + dyb * cost
         zz( 5, 1) = zz( 4, 1) - dyb * sint
         xx( 6, 1) = xx( 5, 1) - dxb
         yy( 6, 1) = yy( 5, 1)
         zz( 6, 1) = zz( 5, 1)
!
         xx( 1, 2) = abcly + 0.5d0 * dxb
         yy( 1, 2) = abclx * cost - 0.5d0 * dyb * cost
         zz( 1, 2) = - 0.5d0 * rbz * cost - 0.5d0 * dyb * sint
         xx( 2, 2) = xx( 1, 2)
         yy( 2, 2) = yy( 1, 2) + dyb * cost
         zz( 2, 2) = zz( 1, 2) + dyb * sint
         xx( 3, 2) = xx( 2, 2) - dxb
         yy( 3, 2) = yy( 2, 2)
         zz( 3, 2) = zz( 2, 2)
!
         xx( 4, 2) = xx( 1, 2)
         yy( 4, 2) = abclx * cost - rbz * sint - 0.5d0 * dyb * cost
         zz( 4, 2) = 0.5d0 * rbz * cost - 0.5d0 * dyb * sint
         xx( 5, 2) = xx( 4, 2)
         yy( 5, 2) = yy( 4, 2) + dyb * cost
         zz( 5, 2) = zz( 4, 2) + dyb * sint
         xx( 6, 2) = xx( 5, 2) - dxb
         yy( 6, 2) = yy( 5, 2)
         zz( 6, 2) = zz( 5, 2)
!
      do 100 i = 1, 6
         xx( i, 3) = xx( i, 1) - 2.0d0 * abcly
         yy( i, 3) = yy( i, 1)
         zz( i, 3) = zz( i, 1)
         xx( i, 4) = xx( i, 2) - 2.0d0 * abcly
         yy( i, 4) = yy( i, 2)
         zz( i, 4) = zz( i, 2)
  100              continue
!
         rcx( 1) = abcly
         rcy( 1) = 0.0d0
         rcz( 1) = - abclx * sint - 0.5d0 * rbz * cost
         rcx( 2) = - rcx( 1)
         rcy( 2) = 0.0d0
         rcz( 2) = rcz( 1)
         rcx( 3) = 0.0d0
         rcy( 3) = abclx * cost - rbz * sint
         rcz( 3) = 0.5d0 * rbz * cost
         rcx( 4) = 0.0d0
         rcy( 4) = - rcy( 3)
         rcz( 4) = rcz( 3)
!
         cosa = dcos( alfab * dpi )
         sina = dsin( alfab * dpi )
         cosb = dcos( betab * dpi )
         sinb = dsin( betab * dpi )
!
      do 200 j = 1, 4
         nline = nline + 1
         idl12(nline) = idyb
         idl23(nline) = idxb
      do 200 i = 1, 6
         x1 = xx( i, j) * cosa - yy( i, j) * sina * cosb &
                              + zz( i, j) * sina * sinb
         y1 = xx( i, j) * sina + yy( i, j) * cosa * cosb &
                              - zz( i, j) * cosa * sinb
         z1 = yy( i, j) * sinb + zz( i, j) * cosb
         xline(6 * nline + i - 6) = x1 + xbcoil
         yline(6 * nline + i - 6) = y1 + ybcoil
         zline(6 * nline + i - 6) = z1 + zbcoil
  200              continue
!
      do 300 i = 1,4
         x1 = rcx( i) * cosa - rcy( i) * sina * cosb &
                            + rcz( i) * sina * sinb
         y1 = rcx( i) * sina + rcy( i) * cosa * cosb &
                            - rcz( i) * cosa * sinb
         z1 = rcy( i) * sinb + rcz( i) * cosb
         narc = narc + 1
         xarc(narc) = x1 + xbcoil
         yarc(narc) = y1 + ybcoil
         zarc(narc) = z1 + zbcoil
  300              continue
!
         cline(nline - 3) = - cbcoil
         cline(nline - 2) = cbcoil
         cline(nline - 1) = cbcoil
         cline(nline) = - cbcoil
!
         aarc(narc - 3) = abclx
         carc(narc - 3) = cbcoil
         idphac(narc - 3) = idphab
         dzarc(narc - 3) = dxb
         drarc(narc - 3) = dyb
         idzarc(narc - 3) = idxb
         idrarc(narc - 3) = idyb
         alfarc(narc - 3) = 90.0d0 + alfab
         betarc(narc - 3) = 90.0d0
         phai1(narc - 3) = - 90.0d0 - 0.5d0 * dphab + betab
         phai2(narc - 3) = - 90.0d0 + 0.5d0 * dphab + betab
!
         aarc(narc - 2) = abclx
         carc(narc - 2) = - cbcoil
         idphac(narc - 2) = idphab
         dzarc(narc - 2) = dxb
         drarc(narc - 2) = dyb
         idzarc(narc - 2) = idxb
         idrarc(narc - 2) = idyb
         alfarc(narc - 2) = 90.0d0 + alfab
         betarc(narc - 2) = 90.0d0
         phai1(narc - 2) = - 90.0d0 - 0.5d0 * dphab + betab
         phai2(narc - 2) = - 90.0d0 + 0.5d0 * dphab + betab
!
         aarc(narc - 1) = abcly
         carc(narc - 1) = - cbcoil
         idphac(narc - 1) = idphyb
         dzarc(narc - 1) = dyb
         drarc(narc - 1) = dxb
         idzarc(narc - 1) = idyb
         idrarc(narc - 1) = idxb
         alfarc(narc - 1) = 180.0d0 + alfab
         betarc(narc - 1) = 180.0d0 - 0.5d0 * dphab - betab
         phai1(narc - 1) = 0.0d0
         phai2(narc - 1) = 180.0d0
!
         aarc(narc) = abcly
         carc(narc) = cbcoil
         idphac(narc) = idphyb
         dzarc(narc) = dyb
         drarc(narc) = dxb
         idzarc(narc) = idyb
         idrarc(narc) = idxb
         alfarc(narc) = 180.0d0 + alfab
         betarc(narc) = 0.5d0 * dphab - betab
         phai1(narc) = 0.0d0
         phai2(narc) = 180.0d0
!
      return
      end
