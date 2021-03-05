module mod_eqneo
  implicit none

contains

!=======================================================================
  subroutine eqneo
!=======================================================================
    use tx_commons, only : nrmax, PsitV, mxneo, fmneo, gamneo, rho, pi, sdt, bbrt, Pisq!, epst, q, rr
    use equ_params
!    use neo
    use mod_spln
!-----	    
    integer(4) :: m, n, miv, ii, nrmaxx
    real(8) :: psix, tpfac, eps, ft
    real(8) :: ztpfv, zepsv
!-----	    
    integer(4), dimension(ivdm), save  :: mxnev
    real(8),    dimension(ivdm), save  :: gamnev, fmnerr
    real(8), dimension(10,ivdm), save  :: fmnev
    real(8),    dimension(ivdm)        :: bbav, bbnav, xx
    real(8), dimension(:), allocatable :: yy

!!$    real(8) :: coefmneo, epsl

    nrmaxx = nrmax + 1
    allocate(yy(0:nrmax))
!=======================================================================
!    geometrical factors of neoclassical diffusion on the equal grid
!=======================================================================
    do n = 2, nv
       psix=real((n-1),8)/real((nv-1),8)
       call eqneo0(psix,ztpfv,zepsv,bbav(n),bbnav(n), &
            &      mxnev(n),fmnev(1,n),gamnev(n),fmnerr(n))
    end do
!-----
!!$!   The PS contribution from the lower poloidal mode number overwhelms that 
!!$!     from the higher ones as you go outside of the torus.
!!$!   Hence, cutoff of maximum mode number, mxnev, inside miv is appropriate
!!$!     for saving computation time.
!!$    miv = 2
!!$    do n = 3, nv
!!$       if(mxnev(n) >= mxnev(n-1)) then
!!$          miv=n-1
!!$          exit
!!$       endif
!!$    end do
!!$    do n = 2, miv
!!$       mxnev(n) = mxnev(miv)
!!$    end do
!   fmnev significantly decreases as m increases.
    do n = 2, nv
       loop_m: do m = 2, 10
          if( fmnev(m,n) < 0.d0 .or. fmnev(m,n) > fmnev(m-1,n) ) then
             mxnev(n) = m - 1
             exit loop_m
          end if
       end do loop_m
    end do
!----- 
!   interpolate the values at the magnetic axis
!-----
    mxnev(1)=1
    gamnev(1)=2.d0*gamnev(2)-gamnev(3)
    do m = 1, 10
       fmnev(m,1)=0.d0
!!   Arbitrary assumption to avoid divergence of flows at axis.
!       fmnev(m,1)=fmnev(m,2)
    end do
    bbnav(1)=0.d0
!-----------------------------------------------------------------------
!!$    write(6,*)'hivnorm,fm1,fm5,err,fam,mx'
!!$    do n = 1, nv
!!$       write(6,'(f10.5,4(a1,1pd10.3),a1,i4)') &
!!$            &     hiv(n)/hiv(nv), &
!!$            & ',',fmnev(1,n),',',fmnev(5,n),',',fmnerr(n),',',gamnev(n), &
!!$            & ',',mxnev(n)
!!$       write(6,'(f10.5,1p10e11.3)') hiv(n)/hiv(nv),fmnev(1:10,n)
!!$    end do
!!$    read(5,*)
!----------------------------------------------------- 	 
!
!   gamneo : <n.grad theta> = gamma = 2 pi / (int B/B_p dl_p) ; (B19) in Ref
!                                   = 4 pi**2 / ( <B> dV/dpsi )
!                                   = 4 pi**2 * sdt / bbrt
!
!----------------------------------------------------- 	 
!!    call spln( gamneo, PsitV, nrmaxx, gamnev, hiv, nv, 151 )
    gamneo(0:nrmax) = 4.d0 * Pisq * sdt(0:nrmax) / bbrt(0:nrmax) ! = bthco / bbrt
!-----
    do m = 1, 10
       do n = 1, nv
          xx(n) = fmnev(m,n)
       end do

       call spln( yy, PsitV, nrmaxx, xx, hiv, nv, 151 )
       do n = 0, nrmax
          fmneo(m,n) = yy(n)
       end do
    end do
!-----
    ii=2
    do n = 0, nrmax
       do 
          if(PsitV(n) > hiv(ii)) then
             ii = ii + 1
             if(ii <= nv) cycle
             ii=nv
          end if
          exit
       end do
       mxneo(n) = mxnev(ii)
    end do

!!$    ! Comparison of fmneo and analytic fmneo
!!$    do n = 1, nrmax
!!$       do m = 1, mxneo(n)
!!$          epsl = epst(n)
!!$          coefmneo = 1.d0 - epsl**2
!!$          write(6,*) n,m,fmneo(m,n),real(m,8)* ( (1.d0-sqrt(coefmneo))/epsl)**(2*m) &
!!$                &                 * (1.d0+real(m,8)*sqrt(coefmneo)) &
!!$                &                 / (coefmneo*sqrt(coefmneo)*(q(N)*RR)**2)
!!$       end do
!!$    end do
!-----------------------------------------------------------------------
!!$    call spln( bavcls, PsitV, nrmaxx, bbav,  hiv, nv, 0 )
!!$    call spln( bbnavr, PsitV, nrmaxx, bbnav, hiv, nv, 0 )
!----------------------------------------------------- end intneo ------
!!$    write(6,'(/)')
!!$    write(6,'(/)')
!!$    write(6,*)'geometrical factors of nc transport' 
!!$    write(6,*)'rho,fm1,fm5,gam,gam2,mx'
!!$    do n = 0, nrmax
!!$       write(6,'(f10.5,4(a1,1pd10.3),a1,i4)') &
!!$            &     rho(n),',',fmneo(1,n),',',fmneo(5,n),',',gamneo(n), &
!!$            &     ',',4.d0*pi**2*sdt(n)/bbrt(n),',',mxneo(n)
!!$!       write(6,'(f10.5,1p10e10.3)') rho(n),fmneo(1:10,n)
!!$    end do
!!$    write(6,'(/)')
!!$    read(5,*)
!=======================================================================
    deallocate(yy)
    return
  end subroutine eqneo

!=======================================================================
  subroutine eqneo0(psix,tpf,eps,bbav,bbnavr, &
       &            mxneo,fmneo,gamneo,fmnerr)
!=======================================================================
!     trapped particle fraction : tpf
!         ft=1.-3./4.*<b**2>*integral{0:1/bmax , e*de/<sqrt(1.-e*b)>}
!=======================================================================
    use tx_commons, only : cnpi => PI, q, rr
    use equ_params
    USE libspl1d
    real(8),    intent(in)  :: psix
    integer(4), intent(out) :: mxneo
    real(8),    intent(out) :: tpf,eps,bbav,bbnavr,gamneo,fmnerr
    real(8)   , intent(out), dimension(10) :: fmneo

    integer(4), parameter :: intf=500, maxfmn = 10
    integer(4) :: i,i1,i2,i3,i4,ix,ir,iz,is,ismax,ist,irst,k,j,j0,j1,jr,jz,m,n,nzsu,ierr
    real(8) :: x,psi0,rbt,rbtt,psi00,scale,xs,rs,zs,smax,rbpx &
    &     ,bb,bb0,bp,bp0,bpdl,bpdl2,bbav0,bbav1,rs0,zs0,s1,s2,s3,s4 &
    &     ,drb,dzb,zrmax,zrmin,zbmax,fintx,cgam0,cgam1,fm,xgam &
    &     ,cs,cs1,cs2,cc,cc1,cc2,tetm,cbb,cbm,bbavr,bthavr,suml,sumxx &
    &     ,sl,dsl,clnb,clnb0,clnb1,coebbl,dcoebbl
    real(8), dimension(isrzdm) :: zlsu, zbsu, zrsu, zzsu, zcsu
    real(8), dimension(0:intf) :: fint,fmu
    real(8), dimension(10) :: sumdim
    real(8), dimension(10000) :: coebp,coebb,ds,blint,zthe,ztet
    real(8), dimension(10000) :: zsl,dcoebb,dztet
    real(8), dimension(4,10000) :: ubb,uztet
!=======================================================================
    x=psix*real((nv-1),8)+1.d0
!!$    ix=x
    ix=nint(x)
!!$    x=abs(x-ix)
!!$    y=1.d0-x
    psi0=saxis*(1.d0-psix)
!!$    rbt=y*rbv(ix)+x*rbv(ix+1)
    rbt=rbv(ix)
    rbtt = rbt * rbt
!-----------------------------------------------------------------------
!     search starting point
!-----------------------------------------------------------------------
    psi00 = psi0
    scale = 1.d0
10  continue
    psi0 = scale * psi00
    ir = iraxis - 1
    iz = izaxis
    i = nsr * ( iz - 1 ) + ir
20  continue
    i = i + 1
    ir = ir + 1
    if( ir >=  nsr )  go to 800
    if( ( psi0 - psi(i) ) * ( psi0 - psi(i+1)) > 0.d0 )  go to 20
    ist = i
    irst = ir
!-----------------------------------------------------------------------
!     start to trace the contour
!-----------------------------------------------------------------------
    is = 1
    k = 1
    j0=i
    j1=i+1
    xs = ( psi0 - psi(j0) ) / ( psi(j1) - psi(j0) )
    rs = rg(ir) + dr * xs
    zs = zg(iz)
    ds(1) = 0.d0
    smax = 0.d0
    rbpx=rbp(j0)+xs*(rbp(j1)-rbp(j0))
    bb = sqrt( rbtt + rbpx*rbpx ) / rs
    bp = rbpx / rs
    coebp(1) = bp
    coebb(1) = bb
    bpdl=0.d0
    nzsu=1
    zrsu(1)=rs
    zzsu(1)=zs
    zcsu(1)=bp
    bbav0=0.d0
    bbav1=0.d0
!-----------------------------------------------------------------------
!     search crossing point
!-----------------------------------------------------------------------
30  continue
    rs0=rs
    zs0=zs
    jr = ir
    jz = iz
    is = is + 1
    i1 = i
    i2 = i1 + 1
    i3 = i2 - nsr
    i4 = i1 - nsr
    s1 = psi0 - psi(i1)
    s2 = psi0 - psi(i2)
    s3 = psi0 - psi(i3)
    s4 = psi0 - psi(i4)
    if( s1*s2 < 0.d0 .and. k /= 1 )  go to 40
    if( s2*s3 < 0.d0 .and. k /= 2 )  go to 50
    if( s3*s4 < 0.d0 .and. k /= 3 )  go to 60
    if( s4*s1 < 0.d0 .and. k /= 4 )  go to 70
    scale = scale + 1.d-04
    go to 10
!-----
40  continue
    xs = s1 / ( s1 - s2 )
    j0 = i1
    j1 = i2
    drb = dr
    dzb = 0.d0
    i = i + nsr
    iz = iz + 1
    k = 3
    go to 80
!-----
50  continue
    xs = s2 / ( s2 - s3 )
    j0 = i2
    j1 = i3
    jr = jr + 1
    drb = 0.d0
    dzb = -dz
    i = i + 1
    ir = ir + 1
    k = 4
    go to 80
!-----
60  continue
    xs = s4 / ( s4 - s3 )
    j0 = i4
    j1 = i3
    jz = jz - 1
    drb = dr
    dzb = 0.d0
    i = i - nsr
    iz = iz - 1
    k = 1
    go to 80
!-----
70  continue
    xs = s1 / ( s1 - s4 )
    j0 = i1
    j1 = i4
    drb = 0.d0
    dzb = -dz
    i = i - 1
    ir = ir - 1
    k = 2
!-----------------------------------------------------------------------
!     out of range ?
!-----------------------------------------------------------------------
80  continue
    if( ir <=     1 )  go to 810
    if( ir >= nsr-1 )  go to 810
    if( iz <=     1 )  go to 810
    if( iz >=   nsz )  go to 810
!-----------------------------------------------------------------------
    bp0=bp
    bb0=bb
    rs = rg(jr) + drb * xs
    zs = zg(jz) + dzb * xs
    ds(is) = sqrt((rs-rs0)**2 + (zs-zs0)**2 )
    smax = smax + ds(is)
    rbpx=rbp(j0)+xs*(rbp(j1)-rbp(j0))
    bb = sqrt( rbtt + rbpx*rbpx ) / rs
    bp = rbpx / rs
    coebp(is) = bp
    coebb(is) = bb
    blint(is)=2.d0*ds(is)/(bp+bp0)
    bpdl=bpdl+blint(is) ! int [ dl / Bp ]
    bbav0=bbav0+2.d0/(bb+bb0)*blint(is)
    bbav1=bbav1+(bb+bb0)/2.d0*blint(is)
!-----
    nzsu=nzsu+1
    zrsu(nzsu)=rs
    zzsu(nzsu)=zs
    zcsu(nzsu)=bp
!-----------------------------------------------------------------------
!     end of trace
!-----------------------------------------------------------------------
    if( i /= ist .or. k /= 1 )  go to 30 
    ismax = is
!=======================================================================
!     trapped particle fraction : tpf
!     inverse aspect ratio : eps
!=======================================================================
!*vocl loop,novrec
    do j = 0, intf
       fint(j)=0.d0
       fmu(j)=real(j,8)/real(intf,8)
    end do
    zrmax=0.d0
    zrmin=100000.d0
    do n = 1, nzsu
       if(zrsu(n) > zrmax) zrmax=zrsu(n)
       if(zrsu(n) < zrmin) zrmin=zrsu(n)
    end do
!<<
!*vocl loop,novrec
    do n = 1, nzsu-1
       zlsu(n)=sqrt((zrsu(n+1)-zrsu(n))**2+(zzsu(n+1)-zzsu(n))**2)
       zrsu(n)=(zrsu(n+1)+zrsu(n))/2.d0
       zzsu(n)=(zzsu(n+1)+zzsu(n))/2.d0
       zcsu(n)=(zcsu(n+1)+zcsu(n))/2.d0
       zbsu(n)=sqrt((rbt/zrsu(n))**2+zcsu(n)**2)
    end do
!<<
    zbmax=0.d0
    bbav=0.d0
!*vocl scalar
    do n = 1, nzsu-1
       bbav=bbav+zlsu(n)*zbsu(n)**2/zcsu(n)
       if(zbsu(n) > zbmax) zbmax=zbsu(n)
    end do
    bbav=bbav/zbmax**2
    do n=1,nzsu-1
       zbsu(n)=zbsu(n)/zbmax
!*vocl loop,novrec
       do j = 0, intf
          fint(j)=fint(j)+zlsu(n)/zcsu(n)*sqrt(1.d0-fmu(j)*zbsu(n))
       end do
    end do
    fintx=0.d0
!*vocl scalar
    do j = 1,intf
       fintx=fintx+(fmu(j)**2-fmu(j-1)**2)/(fint(j)+fint(j-1))
    end do
    fintx=3.d0/4.d0*bbav*fintx
!-----
    tpf=1.d0-fintx
    eps=(zrmax-zrmin)/(zrmax+zrmin)
!=======================================================================	 
!     fourier factors for p-s diffusion : fmneo

!     fm = 2 / < b**2 > < b.grad tet>
!         * [ < ( sin m*the ) ( n.grad b ) >
!           * < ( sin m*the ) ( n.grad b ) ( b.grad the ) >     
!           + < ( cos m*the ) ( n.grad b ) >
!           * < ( cos m*the ) ( n.grad b ) ( b.grad the ) > 
!
!        tet = 2*pi*s / smax ! theta, appearing below (B19) in Ref
!        the = gam * int[ ( b / bp ) dl ] ! Theta, (B18) in Ref
!        gam = 2*pi / int[ ( b / bp ) dl ] ! gamma, (B14) in Ref
!        < n.grad the > = gam
!        < b.grad tet > = 2*pi / int[ dl / bp ]
!        < a > = int[ a * dl / bp ] / int[ dl / bp ]
!
!        sum[ fm ] = < ( n.grad b )**2 > / < b**2 >
!
!     Ref: W.A. Houlberg et al PoP 4 (1997) 3230
!=======================================================================	 
    cgam1 = coebb(1) / coebp(1)
    zthe(1) = 0.d0
    ztet(1) = 0.d0
    zsl(1)  = 0.d0
    do is = 2, ismax
       cgam0=cgam1
       cgam1=coebb(is)/coebp(is)
       zthe(is)=zthe(is-1)+ds(is)
       ztet(is)=ztet(is-1)+ds(is)*(cgam0+cgam1)/2.d0 ! int[ ( b / bp ) dl ]
       zsl (is)=zsl (is-1)+ds(is) ! local perimeter
    end do
    xgam=2.d0*cnpi/ztet(ismax) ! gamma
    do is = 1, ismax
       zthe(is)=2.d0*cnpi*zthe(is)/smax ! theta
       ztet(is)=xgam*ztet(is) ! Theta
!       write(6,'(I3,1P3E15.7)') is,zsl(is),zthe(is),ztet(is)
    end do
!-----
    call spl1d(zsl,coebb,dcoebb,ubb  ,ismax,4,ierr)
    call spl1d(zsl,ztet ,dztet ,uztet,ismax,0,ierr)
!
    bpdl2=bpdl*bpdl
    dsl=smax/(ismax-1)
    do m = 1, maxfmn
       fm=dfloat(m)
       sl=0.d0
       cs1=0.d0
       cs2=0.d0
       cc1=0.d0
       cc2=0.d0
!!$       ! Original
!!$       do is = 2, ismax
!!$          tetm=fm*(ztet(is)+ztet(is-1))/2.d0 ! m * Theta
!!$          cs=sin(tetm)       ! sin(m*Theta)
!!$          cc=cos(tetm)       ! cos(m*Theta)
!!$          cbb=2.d0/(coebb(is)+coebb(is-1)) ! 1/B
!!$!          write(6,*) is,m,tetm,cbb
!!$          cbm=coebb(is)-coebb(is-1) ! delta B
!!$          cs1=cs1+cs*cbb*cbm
!!$          cs2=cs2+cs*cbm*xgam
!!$          cc1=cc1+cc*cbb*cbm
!!$          cc2=cc2+cc*cbm*xgam
!!$!          write(6,'(I3,1P6E15.7)') m,tetm,cs,cbb,cbm,cs1
!!$          write(6,'(2I3,1P6E15.7)') m,is,tetm,coebb(is),coebb(is-1),cbm
!!$       end do

       ! Using spline
       do is = 2, ismax
          sl = (real((is-1),8) - 0.5d0) * dsl
          call spl1df(sl,tetm,zsl,uztet,ismax,ierr)
          tetm=fm*tetm    ! m * Theta
          cs=sin(tetm)    ! sin(m*Theta)
          cc=cos(tetm)    ! cos(m*Theta)
          call spl1dd(sl,coebbl,dcoebbl,zsl,ubb,ismax,ierr)
          cbb=1.d0/coebbl ! 1/B
          cbm=dcoebbl     ! dB/dl
          cs1=cs1+cs*cbb*cbm*dsl
          cs2=cs2+cs*cbm*xgam*dsl
          cc1=cc1+cc*cbb*cbm*dsl
          cc2=cc2+cc*cbm*xgam*dsl
!          write(6,'(I3,1P8E15.7)') m,tetm,sl,cs1,cs2,cc1,cc2
!          if(abs(eps-2.24926d-01)<0.001d0) write(101,'(I3,1P8E15.7)') m,tetm,sl,cs1,cs2,cc1,cc2
!          write(6,'(I3,1P8E15.7)') m,tetm,sl,cc1,cc*cbb*cbm*dsl
       end do

!!$       ! Avoid using dB
!!$       clnb1 = log(coebb(1))
!!$       do is = 2, ismax
!!$          clnb0 = clnb1
!!$          clnb1 = log(coebb(is))
!!$          tetm=fm*0.5d0*(ztet(is)+ztet(is-1)) ! m * Theta
!!$          cs=sin(tetm)       ! sin(m*Theta)
!!$          cc=cos(tetm)       ! cos(m*Theta)
!!$          cbb=0.5d0*(coebb(is)+coebb(is-1)) ! B
!!$          clnb=0.5d0*(clnb0+clnb1) ! log(B)
!!$          cs1=cs1-fm*cc*clnb*(ztet(is)-ztet(is-1))
!!$          cs2=cs2-fm*cc*cbb *(ztet(is)-ztet(is-1))*xgam
!!$          cc1=cc1+fm*cs*clnb*(ztet(is)-ztet(is-1))
!!$          cc2=cc2+fm*cs*cbb *(ztet(is)-ztet(is-1))*xgam
!!$!          write(6,'(I3,1P6E15.7)') m,tetm,(ztet(is)-ztet(is-1)),cs1
!!$!          write(6,'(I3,1P6E15.7)') m,tetm,cs1,cs2,cc1,cc2
!!$       end do

!       cs1=cs1/bpdl           ! cf. (3) in Ref
!       cs2=cs2/bpdl           ! cf. (3) in Ref
!       cc1=cc1/bpdl           ! cf. (3) in Ref
!       cc2=cc2/bpdl           ! cf. (3) in Ref
!       fmneo(m)=cs1*cs2+cc1*cc2 ! [ ] in (B9) of Ref
       fmneo(m)=(cs1*cs2+cc1*cc2)/bpdl2 ! [ ] in (B9) of Ref
!       write(101,'(1pe12.5,4x,i4,1x,1p5e13.5)') eps,m,cs1,cs2,cc1,cc2,fmneo(m)
!       write(6,'(1pe12.5,4x,i4,1x,1p5e13.5)') eps,m,cs1/bpdl2,(-1.d0)**(m+1)*((1.d0-sqrt(1.d0-eps**2))/eps)**m/(1.d0+0.5d0*eps**2)/(qqv(ix)*rbt)
!!$       write(6,'(1pe12.5,4x,i4,1x,1p5e13.5)') eps,m,fmneo(m),&
!!$            &         real(m,8)*( (1.D0-SQRT(1.D0-eps**2))/eps)**(2*m) &
!!$            &         *(1.D0+real(m,8)*SQRT(1.D0-eps**2))/((1.D0-eps**2)**1.5D0 &
!!$            &         *(qqv(ix)*rr)**2),&
!!$            &         cs1,cs2
    end do
!------
    bbavr=0.d0
    bbnavr=0.d0
    do is = 2, ismax
       bbavr=bbavr &
          & +0.5d0*(coebb(is)**2+coebb(is-1)**2)*blint(is)
       bbnavr=bbnavr &
          & +0.5d0*(coebp(is)/coebb(is)**2+coebp(is-1)/coebb(is-1)**2) &
          & *(coebb(is)-coebb(is-1))**2/ds(is)
    end do
    bbavr=bbavr/bpdl ! <B^2>
    bbnavr=bbnavr/bpdl        ! <(n.nabla B)^2>
    bthavr=2.d0*cnpi/bpdl     ! <B.nabla theta>
!-----
    gamneo=xgam
    do m = 1, maxfmn
       fmneo(m)=2.d0*fmneo(m)/(bbavr*bthavr) ! Fm, (B9)
!!$       aaa=(1.d0-sqrt(1.d0-eps**2))/eps
!!$       write(101,'(1pd12.5,4x,i4,1x,1p4d12.5)')eps,m,fmneo(m), &
!!$            & dble(m)/(1.d0-eps**2)**1.5d0 & 
!!$            & *(1.d0+dble(m)*sqrt(1.d0-eps**2))*aaa**(2*m) &
!!$            & *(rbv(1)/sqrt(rrv(1))/(qqv(ix)*rbv(ix)))**2 &
!!$            & ,8.d0*aaa**(2*m)*(rbv(1)/sqrt(rrv(1))/(qqv(ix)*rbv(ix)))**2 &
!!$            & ,2.d0*aaa**(2*m)*(rbv(1)/sqrt(rrv(1))/(qqv(ix)*rbv(ix)))**2
    end do
!=======================================================================
!      write(6,'(//)')
!     write(6,*)' ===== eqneo output ====='
!	  write(6,'(a10,f10.5)')' psix    =',psix
!	  write(6,'(a10,f10.5)')' tpf     =',tpf
!	  write(6,'(a10,f10.5)')' eps     =',eps
    suml = 0.d0
    do m = 1, maxfmn
       suml=suml+fmneo(m)
       sumdim(m)=suml*bbavr/bbnavr
       fmnerr=1.d0-sumdim(m)
!       write(6,'(4x,i4,1p4d12.5)')m,sumdim(m),fmneo(m),suml,bbnavr/bbavr
    end do
!    sumxx=0.999d0*sumdim(maxfmn)
    sumxx=(1.d0-1.d-5)*sumdim(maxfmn)
    do m = 1, maxfmn
       mxneo=m
       if(sumdim(m) > sumxx) exit
    end do
!!$    write(6,'(3f10.5,2x,1p4d10.3)')psix,tpf,eps,fmneo(1),fmneo(5),fmnerr,gamneo
!!$    read(5,*)
!=======================================================================
    bbav0=bbav0/bpdl
    bbav1=bbav1/bpdl
    bbav=bbav0-1.d0/bbav1
!=======================================================================	  
    return
800 continue
810 continue
    stop 'error stop'
  end subroutine eqneo0

  ! ********************************************************************
  !      Coefficients for Pfirsch-Schluter viscosity (nccoe)
  ! ********************************************************************

  subroutine wrap_eqneo
    use tx_commons, only : ieqread, NRMAX, gamneo, Pisq, sdt, bbrt, mxneo, epst, fmneo, q, RR
    integer(4) :: NR, i
    real(8) :: epsl, coefmneo, smallvalue = 1.d-4

    if( ieqread >= 2 ) then
       call eqneo
    else
       do NR = 1, NRMAX
          gamneo(NR) = 4.d0 * Pisq * sdt(NR) / bbrt(NR) ! = bthco(NR) / bbrt(NR)
          mxneo(NR) = 3
          epsl = epst(NR)
          coefmneo = 1.d0 - epsl**2
          do i = 1, mxneo(NR)
             fmneo(i,NR) = real(i,8)* ( (1.d0-sqrt(coefmneo))/epsl)**(2*i) &
                  &                 * (1.d0+real(i,8)*sqrt(coefmneo)) &
                  &                 / (coefmneo*sqrt(coefmneo)*(q(NR)*RR)**2)
          end do
          fmneo(mxneo(NR)+1:10,NR) = 0.d0
       end do
    end if
    ! Even at axis, fmneo=0 should be avoided because it would cause K_PS = 0.
    NR = 0
       gamneo(NR) = 4.d0 * Pisq * sdt(NR) / bbrt(NR) ! = bthco(NR) / bbrt(NR)
       mxneo(NR) = 3
!       fmneo(1:10,NR) = 0.d0
       epsl = smallvalue
       coefmneo = 1.d0 - epsl**2
       do i = 1, mxneo(NR)
          fmneo(i,NR) = real(i,8)* ( (1.d0-sqrt(coefmneo))/epsl)**(2*i) &
               &                 * (1.d0+real(i,8)*sqrt(coefmneo)) &
               &                 / (coefmneo*sqrt(coefmneo)*(q(NR)*RR)**2)
       end do
       fmneo(mxneo(NR)+1:10,NR) = 0.d0

  end subroutine wrap_eqneo

end module mod_eqneo
