module mod_trace
  use tx_commons, only : pi, irktrc
  use equ_params, only : rmaxis => raxis, zmaxis => zaxis, btv, vlv, nv, ckv, sdw !(toms760), nw => nsr, nh => nsz, rg, zg, psi
  implicit none
  private

  public :: magtrace, eqmags, psigd, psig
  
  ! For Brent method
  integer(4) :: nzfunc

  integer(4) :: itrace = 0
  real(8), allocatable :: xaa(:), uz(:,:), upsi1d(:,:)
  logical :: lclockwise, lX_lownull

contains

  ! compute quantities for which tracing the surface is required

  subroutine magtrace(rinit,zinit,nn,v)
    use subs
    use mod_num_recipe
    use equ_params, only : tol, miller, ivdm, isrzdm, arv, vlv, sdw, ckv, ssv, aav, rrv &
         &                 , bbv, biv, shv, grbm2v, brv, aiv, rbv, r2b2v, rtv, rpv, epsv &
         &                 , elipv, trigv, ftv, gttiv, dsr, dsz, csu, rsu, zsu, nsu, sigcu &
         &                 , raxis, zaxis, siw, nmax
    use libspl1d, only : spl1d, spl1dd, spl1ddd, spl1df
    real(8), intent(in) :: rinit, zinit
    integer, intent(in) :: nn
    type(miller), intent(inout) :: v
    ! type(miller) includes: Rgeo, Zgeo, rminor, elip, trig, zeta, Rmax, Rmin, zpla

    integer(4) :: n, nround, nrmax, nrmin, nzmin, nzmax, ierr, nbpmin, n1, n2, i, j, ind
    real(8) :: xa(nmax), ya(2,nmax)
    real(8) :: r, z, h, dz, srdz, srzdz, bpmin!, btl, b2l, sumv
    real(8) :: zdpsidr, zdpsidz, elipup, elipdw, trigup, trigdw, zetaup, zetadw
    real(8) :: rzmax, zzmax, rzmin, zzmin, rrmax, rrmin
!    real(8) :: thup, thdw, zdist, zminup, zmindw, thupp, thdwp, zupp, zdwp
    real(8) :: rcup, rcdw, zcup, zcdw, dabup, dacup, dcdup, dabdw, dacdw, dcddw, zdist, zminup, zmindw
    real(8) :: rbup, zbup, rbdw, zbdw, tanup, tandw, rpmax, zmid0
    real(8) :: zsign
    integer(4) :: nup, ndw

    real(8), allocatable :: yar(:), yaz(:), dyar(:), dyaz(:), ur(:,:)
    real(8) :: zl, zlold, rl, rlold, dx, x, zmax, rmax, zmin, rmin, tmp, dtmp
    real(4) :: t1, t2, t3, t4

    real(8) :: dl,dr10,dz10,ds0,ds1,bp0,bp1,r0,r1,z0,z1,vl0,vl1,ck0,ck1,ss0,ss1,aa0,aa1,rr0,rr1 &
         &    ,bb0,bb1,bi0,bi1,bm0,bm1,br0,br1,sh0,sh1,ai0,ai1,dsr0,dsr1,dsz0,dsz1,sdw2

    ! --- trapped particle fraction ---
    ! intf: num. of division in the direction of lambda for trapped particle fraction
    integer(4), parameter :: intf = 100
    real(8) :: fintx, hsq
    integer(4), dimension(:), allocatable :: nsul
    real(8), dimension(:), allocatable :: bmax,fint,flam,dll,zbl,zbpl
    ! --- gttiv ---
    real(8) :: fac
    real(8), dimension(:), allocatable :: rrl, zzl, dsrl, dszl
    ! ---------------------------------

    if( nn == 1 ) then
       ! At the magnetic axis
       arv(1)   = 0.d0

       vlv(1)   = 0.d0
       ! Interpolate sdw(1) by using dpsi/d(sqrt(V))=0 at V=0
       sdw(1)   = ( vlv(3) * sdw(2) - vlv(2) * sdw(3) ) / ( vlv(3) - vlv(2) )
       ckv(1)   = 0.d0
       r2b2v(1) = 0.d0
       ssv(1)   = 0.d0
       aav(1)   = 1.d0 / raxis**2
       rrv(1)   = raxis**2
       bbv(1)   = (rbv(1) / raxis)**2
       biv(1)   = 1.d0 / bbv(1)
       shv(1)   = 0.d0
       grbm2v(1)= 0.d0
       aiv(1)   = 1.d0 / raxis
       brv(1)   = rbv(1) / raxis
       
       rtv(1)   = raxis
       rpv(1)   = 0.d0
       epsv(1)  = 0.d0
!       elipv(1) = 1.d0
       elipv(1) = elipv(2)
       trigv(1) = 0.d0
       
       ftv(1)   = 0.d0

       ! For parameters related to Miller ones

       v%Rgeo   = raxis
       v%Zgeo   = zaxis
       v%rminor = 0.d0
       v%elip   = elipv(2)
       v%trig   = 0.d0
       v%Rmax   = raxis
       v%Rmin   = raxis
       v%zpla   = raxis

       return
    end if

    allocate(bmax(ivdm),fint(0:intf),flam(0:intf))
    allocate(nsul(isrzdm))
    allocate(dll(isrzdm),zbl(isrzdm),zbpl(isrzdm),rrl(isrzdm),zzl(isrzdm),dsrl(isrzdm),dszl(isrzdm))
    do j = 0, intf
       flam(j)  = real( j, 8 ) / real( intf, 8 )
    end do

    ! --- compute all (R,Z) points on the flux surface, anti-clockwise ; xa, ya and nround
    if( nn /= nv ) then
       ind = 0
    else
       ind = -1
    end if
    call eqmags(rinit,zinit,nmax,xa,ya,nround,ind)
    if( ind /= 0 ) then
       write(6,'(X,A,I2)') 'eqmags error: ierr= ', ind
       stop
    end if
    if( lclockwise ) then
       zsign =-1.d0
    else
       zsign = 1.d0
    end if

    ! --- calculate Z_geo and metric(vlv,ckv,sdw)
    srdz  = 0.d0
    srzdz = 0.d0

    arv(nn) = 0.d0
    vlv(nn) = 0.d0
    sdw(nn) = 0.d0
    ckv(nn) = 0.d0
    ssv(nn) = 0.d0
    aav(nn) = 0.d0
    rrv(nn) = 0.d0
    bbv(nn) = 0.d0
    biv(nn) = 0.d0
    shv(nn) = 0.d0
    grbm2v(nn) = 0.d0
    brv(nn) = 0.d0
    aiv(nn) = 0.d0
    bmax(nn) = 0.d0

    bpmin = abs(btv/rmaxis)

    vl1 = ya(1,1) * ya(1,1)
    r1  = ya(1,1)
    z1  = ya(2,1)
    call psigd(r1,z1,zdpsidr,zdpsidz)
    bp1 = sqrt(zdpsidr**2+zdpsidz**2)/r1 ! Bp
    dsr1= zdpsidr
    dsz1= zdpsidz
    ds1 = 1.d0 / bp1
    ck1 = bp1
    ss1 = vl1 * bp1
    aa1 = ds1 / vl1
    rr1 = vl1 / bp1
    bb1 = (rbv(nn) / r1)**2 + bp1**2
    bi1 = ds1 / bb1
    bm1 = sqrt(bb1)
    br1 = bm1 * ds1
    bb1 = bb1 * ds1
    sh1 = r1
    ai1 = ds1 / r1

    rrmax = r1
    rrmin = r1
    zzmax = z1
    zzmin = z1
    rzmax =-1000.d0
    rzmin = 1000.d0

    nsul(nn) = 1
    if(nn == nv) then
       nsu    = 1
       rsu(1) = r1
       zsu(1) = z1
       csu(1) = bp1*sigcu
    endif

    do n = 2, nround
       r     = ya(1,n) ! current R
       z     = ya(2,n) ! current Z
       h     = xa(n) - xa(n-1) ! arc length
       dz    = ya(2,n) - ya(2,n-1)

       bp0 = bp1
       dsr0= dsr1
       dsz0= dsz1
       vl0 = vl1
       r0  = r1
       z0  = z1

       ds0 = ds1
       ck0 = ck1
       ss0 = ss1
       aa0 = aa1
       rr0 = rr1
       bb0 = bb1
       bi0 = bi1
       bm0 = bm1 
       br0 = br1
       sh0 = sh1
       ai0 = ai1

       ! For elevation
       srdz  = srdz  + 0.5d0 * (ya(1,n) + ya(1,n-1)) * dz
       srzdz = srzdz + 0.5d0 * (ya(1,n) + ya(1,n-1)) * 0.5d0 * (ya(2,n) + ya(2,n-1)) * dz
       call psigd(ya(1,n),ya(2,n),zdpsidr,zdpsidz)
       bp1 = sqrt(zdpsidr**2+zdpsidz**2)/r
       if( bp1 < bpmin ) then
          bpmin  = bp1
          nbpmin = n
       endif
!--- calc metric
       vl1 = ya(1,n) * ya(1,n)
       z1  = ya(2,n)
       r1  = ya(1,n)
       z1  = ya(2,n)

       ds1 = 1.d0 / bp1
       ck1 = bp1
       ss1 = vl1 * bp1
       aa1 = ds1 / vl1
       rr1 = vl1 / bp1
       bb1 = (rbv(nn) / r1)**2 + bp1**2
       bi1 = ds1 / bb1
       bm1 = sqrt(bb1)
       br1 = bm1 * ds1
       bb1 = bb1 * ds1
       sh1 = r1
       ai1 = ds1 / r1

       dr10 = r1 - r0
       dz10 =(z1 - z0) * zsign
       dl   = sqrt(dr10 * dr10 + dz10 * dz10)
       arv(nn) = arv(nn) + dz10 * (r0 +r1 ) * 0.5d0
       vlv(nn) = vlv(nn) + dz10 * (vl0+vl1) * 0.5d0
       sdw(nn) = sdw(nn) + dl * (ds0 + ds1) * 0.5d0
       ckv(nn) = ckv(nn) + dl * (ck0 + ck1) * 0.5d0
       ssv(nn) = ssv(nn) + dl * (ss0 + ss1) * 0.5d0
       aav(nn) = aav(nn) + dl * (aa0 + aa1) * 0.5d0
       rrv(nn) = rrv(nn) + dl * (rr0 + rr1) * 0.5d0
       bbv(nn) = bbv(nn) + dl * (bb0 + bb1) * 0.5d0
       biv(nn) = biv(nn) + dl * (bi0 + bi1) * 0.5d0
       shv(nn) = shv(nn) + dl * (sh0 + sh1) * 0.5d0
       grbm2v(nn) = grbm2v(nn) + dl * (  r0**2 * bp0**2 * bi0 &
            &                          + r1**2 * bp1**2 * bi1) * 0.5d0
       aiv(nn) = aiv(nn) + dl * (ai0 + ai1) * 0.5d0
       brv(nn) = brv(nn) + dl * (br0 + br1) * 0.5d0

       ! nsul(nn) is now counting up.
       rrl (nsul(nn)) =         (r0  + r1 ) * 0.5d0 ! R
       zzl (nsul(nn)) =         (z0  + z1 ) * 0.5d0 ! Z
       dll (nsul(nn)) =    dl * (ds0 + ds1) * 0.5d0 ! dlp/Bp
       zbl (nsul(nn)) =         (bm0 + bm1) * 0.5d0 ! B
       zbpl(nsul(nn)) =         (bp0 + bp1) * 0.5d0 ! Bp
       dsrl(nsul(nn)) =         (dsr0+ dsr1)* 0.5d0 ! dpsi/dR
       dszl(nsul(nn)) =         (dsz0+ dsz1)* 0.5d0 ! dpsi/dZ

       bmax(nn) = max(bmax(nn),zbl(nsul(nn)))

!--- geometric factors

       if(r1 > rrmax) rrmax = r1 ! rrmax = R_max
       if(r1 < rrmin) rrmin = r1 ! rrmin = R_min
       if(z1 > zzmax) then
          zzmax = z1 ! zzmax = Z_max
          rzmax = r1 ! rzmax = R at Z_max
       end if
       if(z1 < zzmin) then
          zzmin = z1 ! zzmin = Z_min
          rzmin = r1 ! rzmin = R at Z_min
       end if

       nsul(nn) = nsul(nn) + 1
       if(nn == nv) then
          nsu      = nsu + 1
          rsu(nsu) = r1
          zsu(nsu) = z1
          csu(nsu) = bp1 * sigcu
       endif
    enddo
!--- 
    vlv(nn)    =        pi * vlv(nn)
    sdw(nn)    = 1.d0 / (2.d0 * pi * sdw(nn))
    ckv(nn)    = 2.d0 * pi * ckv(nn) / sdw(nn)
    r2b2v(nn)  = 2.d0 * pi * ssv(nn) * sdw(nn)
    ssv(nn)    = 2.d0 * pi * ssv(nn) / sdw(nn)
    aav(nn)    = 2.d0 * pi * aav(nn) * sdw(nn)
    rrv(nn)    = 2.d0 * pi * rrv(nn) * sdw(nn)
    bbv(nn)    = 2.d0 * pi * bbv(nn) * sdw(nn)
    biv(nn)    = 2.d0 * pi * biv(nn) * sdw(nn)
    shv(nn)    = 2.d0 * pi * shv(nn) * sdw(nn)
    grbm2v(nn) = 2.d0 * pi * grbm2v(nn) * sdw(nn)
    aiv(nn)    = 2.d0 * pi * aiv(nn) * sdw(nn)
    brv(nn)    = 2.d0 * pi * brv(nn) * sdw(nn)
    sdw(nn)    = sdw(nn) * sigcu

    rtv(nn)   = 0.5d0 * (rrmax + rrmin)
    rpv(nn)   = 0.5d0 * (rrmax - rrmin)
    epsv(nn)  = (rrmax - rrmin) / (rrmax + rrmin)
    elipv(nn) = (zzmax - zzmin) / (rrmax - rrmin)
    trigv(nn) = (rtv(nn) - 0.5d0 * (rzmax + rzmin)) / rpv(nn)

    ! *** trapped particle fraction **********************************
    !
    !   ft = 1 - 0.75 <h^2> int_0^1 flam dflam / <sqrt(1 - flam * h)>
    !
    ! ****************************************************************

    nsul(nn) = nsul(nn) - 1
    fint(:) = 0.d0
    do i = 1, nsul(nn)
       h = zbl(i) / bmax(nn) ! h
       do j = 0, intf
          if( 1.d0 - flam(j) * h < 0.d0 ) then
             write(6,'(a,es26.18,a)') &
                  &        '** WARNING  IN SQRT(DX),DX < 0.0(DX=',1.d0-flam(j)*h,').'
          endif
          fint(j) = fint(j) + sqrt( abs( 1.d0 - flam(j) * h ) ) * dll(i)
       enddo
    enddo

    fintx = 0.d0
    do j = 1, intf
       fintx = fintx + ( flam(j)**2 - flam(j-1)**2 ) / ( fint(j) + fint(j-1) )
    enddo
    hsq     = bbv(nn) / bmax(nn)**2 ! <h^2>
    fintx   = 0.75d0 * hsq * fintx / (2.d0 * pi * sdw(nn))
    ftv(nn) = 1.d0 - fintx

    !*** metrics *****************************************************
    !  
    !  gttiv : g^{theta,theta}^{-1} = (e_theta . e_theta)^{-1}
    ! 
    !*****************************************************************

    gttiv(nn) = 0.d0
    do i = 1, nsul(nn)
       fac = dsrl(i) / dszl(i)
       gttiv(nn) = gttiv(nn) + 1.d0 / (rrl(i)**4 *zbpl(i)**2) * dll(i)
    enddo
    gttiv(nn) = (2.d0 * pi * sdw(nn)) * gttiv(nn) * (4.d0 * pi**2 * sdw(nn) / aav(nn))**2

    nsu     = nsu - 1

    v%Zgeo = srzdz / srdz

    deallocate(bmax,fint,flam,nsul,dll,zbl,zbpl,rrl,zzl,dsrl,dszl)

!------------------------------------------------------------------------------------------------------
!    if(nn /= nv) return!if nn ==nv, calc rmaj,rpla etc...
!------------------------------------------------------------------------------------------------------

    ! --- calculate interim Zmax and Zmin and corrensponding R's
    !     using trace points only 

    zzmin = ya(2,1) ; zzmax = ya(2,1)
    do n = 2, nround
       ! call psigd(ya(1,n),ya(2,n),zdpsidr,zdpsidz)
       r     = ya(1,n) ! current R
       z     = ya(2,n) ! current Z
       h     = xa(n) - xa(n-1) ! arc length

       ! linearly interpolate Rmax and Rmin for calculating geometrical center of the surface of interest
       if( (ya(2,n) - v%Zgeo)*(ya(2,n-1) - v%Zgeo) < 0 ) then
          if( ya(1,n) > rmaxis ) then
!             v%Rmax = (ya(1,n) - ya(1,n-1))/(ya(2,n) - ya(2,n-1))*(v%Zgeo - ya(2,n-1))+ya(1,n)
             nrmax = n
          else
!             v%Rmin = (ya(1,n) - ya(1,n-1))/(ya(2,n) - ya(2,n-1))*(v%Zgeo - ya(2,n-1))+ya(1,n)
             nrmin = n
          endif
       endif

       ! For elongation and triangularity
       zzmin = min(zzmin,ya(2,n))
       if( zzmin == ya(2,n) ) nzmin = n
       zzmax = max(zzmax,ya(2,n))
       if( zzmax == ya(2,n) ) nzmax = n
    enddo
!    v%Rgeo   = 0.5d0*(v%Rmax+v%Rmin)
!    v%rminor = 0.5d0*(v%Rmax-v%Rmin)

    ! ==== More accurately estimate (R_max,Z_geo) and (R_min,Z_geo) ====

    ! ******************************************************************    | 
    ! ***                                                            ***    |
    ! ***  More accurately estimate (R_max,Z_geo) and (R_min,Z_geo)  ***  \ | /
    ! ***                                                            ***   \|/
    ! ******************************************************************    V

    allocate(xaa(nround),yar(nround),yaz(nround),dyar(nround),dyaz(nround),ur(4,nround),uz(4,nround))
    ! Spline R and Z over x(n)

    xaa(1:nround) = xa(1:nround)   ! x(n) : arc length from the starting point
    yar(1:nround) = ya(1,1:nround) ! R(n)
    yaz(1:nround) = ya(2,1:nround) ! Z(n)
    call spl1d(xaa,yar,dyar,ur,nround,4,ierr)
    if( ierr /= 0 ) stop 'XX error in spl1d for creating yar(xaa)'
    call spl1d(xaa,yaz,dyaz,uz,nround,4,ierr)
    if( ierr /= 0 ) stop 'XX error in spl1d for creating yaz(xaa)'

    ! --- find maximum R on Z_geo plane

    ! find the position of the arc length x satisfying Z(x)=v%Zgeo in LFS,
    ! i.e., x at which R=Rmax
    call newton1d(xaa(nround),x,spl1dd,xaa,uz,nround,tol,v%Zgeo)
    ! x(n) is of course periodic. If newton1d finds x at which R=Rmax exceeding xaa(nround),
    ! where xaa(nround) is the end of tracing, xaa(nround) should be subtracted from x.
    if( x > xaa(nround) ) x = x - xaa(nround)
    ! get R(x)=Rmax
    call spl1df(x,v%Rmax,xaa,ur,nround,ierr)

    ! --- find minimum R on Z_geo plane

    ! find the position of the arc length x satisfying Z(x)=v%Zgeo in HFS,
    ! i.e., x at which R=Rmin
    !   Rmin must be near the half of tracing, i.e., n = nround/2
    call newton1d(xaa(nround/2),x,spl1dd,xaa,uz,nround,tol,v%Zgeo)
    ! get R(x)=Rmin
    call spl1df(x,v%Rmin,xaa,ur,nround,ierr)

    if( v%Rmin > v%Rmax ) then
       tmp    = v%Rmax
       v%Rmax = v%Rmin
       v%Rmin = tmp
    elseif (  v%Rmax - v%Rmin < epsilon(1.d0) ) then
       stop 'Cannot find Rmax and Rmin'
    endif

    ! ===========================================================================

    v%Rgeo   = 0.5d0*(v%Rmax+v%Rmin)
    v%rminor = 0.5d0*(v%Rmax-v%Rmin)

    ! *************************************************************************    | 
    ! ***                                                                   ***    |
    ! ***  More accurately estimate (R_{Z_max},Z_max) and (R_{Z_min},Zmin)  ***  \ | /
    ! ***                                                                   ***   \|/
    ! *************************************************************************    V

    ! --- find maximum Z and corresponding R

    !  Find x satisfying Z(x) = Zmax by seeking dZ(x)/dx = 0 in the upper half region.
    nzfunc = nround
    if( lclockwise ) then
       ! Zmax must be near the 3/4 of tracing when clockwise.
       ! Hence, nround/4*3 is bounded between nround/8*5 and nround/8*7.
!       n = nround/4*3
       n1 = nround/8*5
       n2 = nround/8*7
    else
       ! Zmax must be near the 1/4 of tracing when anti-clockwise.
       ! Hence, nround/4 is bounded between nround/8 and nround/8*3.
!       n = nround/4
       n1 = nround/8
       n2 = nround/8*3
    end if

!    ! get Zmax
!    call newton1dd(xaa(n),x,zmax,spl1ddd,xaa,uz,nround,tol)

    ! ** Newton method does not work especially near X point. Instead, one out of (1)-(3) is chosen.
    ! (1) use bisection method
!    call bisectiond(0.d0,x,zmax,spl1dd,xaa,uz,nround,xaa(n1),xaa(n2),tol)

    ! (2) use Brent method (computation speed comparable to fine tracing)
!    x = fbrent(dzfunc,xaa(n1),xaa(n2),tol)

    ! (3) use combined method of bisection and newton (faster than Brent method by a factor of two)
    x = rtsafe(dzsub,xaa(n1),xaa(n2),tol)

    ! get Zmax
    call spl1df(x,zmax,xaa,uz,nround,ierr)

    ! get R at Zmax : rmax
    call spl1df(x,rmax,xaa,ur,nround,ierr)

    ! --- find minimum Z and corresponding R

    !  Find x satisfying Z(x) = Zmin by seeking dZ(x)/dx = 0 in the lower half region.
    nzfunc = nround
    if( lclockwise ) then
       ! Zmin must be near the 1/4 of tracing when anti-clockwise.
       ! Hence, nround/4 is bounded between nround/8 and nround/8*3.
!       n = nround/4
       n1 = nround/8
       n2 = nround/8*3
    else
       ! Zmin must be near the 3/4 of tracing when clockwise.
       ! Hence, nround/4*3 is bounded between nround/8*5 and nround/8*7.
!       n = nround/4*3
       n1 = nround/8*5
       n2 = nround/8*7
    end if

!    ! get Zmin
!    call newton1dd(xaa(n),x,zmin,spl1ddd,xaa,uz,nround,tol)

    ! ** Newton method does not work especially near X point. Instead, one out of (1)-(3) is chosen.
    ! (1) use bisection method
!    call bisectiond(0.d0,x,zmin,spl1dd,xaa,uz,nround,xaa(n1),xaa(n2),tol)

    ! (2) use Brent method (computation speed comparable to fine tracing)
!    x = fbrent(dzfunc,xaa(n1),xaa(n2),tol)

    ! (3) use combined method of bisection and newton (faster than Brent method by a factor of two)
    x = rtsafe(dzsub,xaa(n1),xaa(n2),tol)

    ! get Zmin
    call spl1df(x,zmin,xaa,uz,nround,ierr)

    ! get R at Zmin : rmin
    call spl1df(x,rmin,xaa,ur,nround,ierr)

    ! ===========================================================================

!    rzmax  = ya(1,nzmax)
!    rzmin  = ya(1,nzmin)
    rzmax = rmax
    rzmin = rmin
    zzmax = zmax
    zzmin = zmin

    elipup = ( zzmax  - v%Zgeo ) / v%rminor
    elipdw = ( v%Zgeo - zzmin  ) / v%rminor
    trigup = ( v%Rgeo - rzmax  ) / v%rminor
    trigdw = ( v%Rgeo - rzmin  ) / v%rminor

    v%elip = 0.5d0 * ( elipup + elipdw )
    v%trig = 0.5d0 * ( trigup + trigdw )
    v%zpla = 0.5d0 * ( rzmax  + rzmin  )
!    write(6,'(6ES15.7)') v%rminor,v%Rgeo,rzmax,rzmin,trigup,trigdw

    ! Holcomb's squareness [C.T. Holcomb PoP 2009][T. Luce PPCF 2013] defined on the midplane
    ! Aup : (rzmax,zmid0), Adw : (rzmin,zmid0), AR : (rpmax,zmid0)
    ! AZup: (rzmax,zzmax), AZdw: (rzmin,zzmin), Dup: (rpmax,zzmax), Ddw: (rpmax,zzmin)
    ! Cup : (rcup ,zcup ), Cdw : (rcdw ,zcdw ), Bup: (rbup ,zbup ), Bdw: (rbdw ,zbdw )

!!$    rpmax = rinit
!!$    zmid0 = zinit
    rpmax = v%Rmax
    zmid0 = v%Zgeo

    ! upper region
    rcup  = rzmax + ( rpmax - rzmax ) / sqrt(2.d0) ! C
    zcup  = zmid0 + ( zzmax - zmid0 ) / sqrt(2.d0) ! C
    dcdup = distance(rcup,zcup,rpmax,zzmax)    ! CD
    dacup = distance(rzmax,zmid0,rcup,zcup)    ! AC
    tanup = ( zzmax - zmid0 ) / ( rpmax - rzmax )

    ! lower region
    rcdw  = rzmin + ( rpmax - rzmin ) / sqrt(2.d0) ! C
    zcdw  = zmid0 + ( zzmin - zmid0 ) / sqrt(2.d0) ! C
    dcddw = distance(rcdw,zcdw,rpmax,zzmin)    ! CD
    dacdw = distance(rzmin,zmid0,rcdw,zcdw)    ! AC
    tandw = ( zzmin - zmid0 ) / ( rpmax - rzmin )

    zminup = 1.d0
    zmindw = 1.d0
    nup = 0
    ndw = 0
    do n = 2, nround
       if( ya(1,n) > min(rzmax,rzmin) ) then ! LFS only
          if( ya(2,n) > zmid0 ) then ! upper region
             ! Z coordinate at the designated R coordinate, ya(1,n), on line AD
             z = linZ(ya(1,n),rzmax,zmid0,rpmax,zzmax)
             zdist = ya(2,n) - z
             ! find the point where Z, ya(2,n), approaches line AD most.
             if( abs(zdist) < zminup ) then
                zminup = abs(zdist)
                if( zdist >= 0.d0 ) then
                   nup = n
                else
                   nup = n + 1
                endif
             endif
!             write(6,'(2I4,5ES15.7)') n,nup,ya(1,n),ya(2,n),z,zdist,zminup
          else ! lower region
             ! Z coordinate at the designated R coordinate, ya(1,n), on line AD
             z = linZ(ya(1,n),rzmin,zmid0,rpmax,zzmin)
             zdist = ya(2,n) - z
             ! find the point where Z, ya(2,n), approaches line AD most.
             if( abs(zdist) < zmindw ) then
                zmindw = abs(zdist)
                if( zdist >= 0.d0 ) then
                   ndw = n
                else
                   ndw = n + 1
                endif
             endif
!             write(6,'(2I4,5ES15.7)') n,ndw,ya(1,n),ya(2,n),z,zdist,zmindw
          endif
       end if
    enddo
    
!    call intersection(rbup,zbup,ya(1,nup-1),ya(2,nup-1),ya(1,nup),ya(2,nup),rzmax,zmid0,tanup) ! B
!    call intersection(rbdw,zbdw,ya(1,ndw-1),ya(2,ndw-1),ya(1,ndw),ya(2,ndw),rzmin,zmid0,tandw) ! B

    ! ==== more accurately estimating point B for squareness ====

    ! find point B in upper half region by the method chosen among (1)-(3)

    ! (1) use Newton method
    call newton1df(xaa(nup-1),x,spl1dd,xaa,uz,spl1dd,ur,nround,tol,linZ,rzmax,zmid0,rpmax,zzmax)
    call spl1df(x,rbup,xaa,ur,nround,ierr)
    call spl1df(x,zbup,xaa,uz,nround,ierr)

    ! (2) use Brent method (computation speed comparable to fine tracing)
!    r1 = rzmax ; z1 = zmid0 ; r2 = rpmax ; z2 = zzmax ; nzfunc = nround
!    x = fbrent(zfunc,xaa(nup-1),xaa(nup),tol)
!    call spl1df(x,rbup,xaa,ur,nround,ierr)
!    call spl1df(x,zbup,xaa,uz,nround,ierr)

    ! (3) use bisection instead of newton1df (very slow)
!    call bisectionf(linZ,x,spl1df,xaa,uz,spl1df,ur,nround,xaa(nup-1),xaa(nup),tol,rzmax,zmid0,rpmax,zzmax)
!    call spl1df(x,rbup,xaa,ur,nround,ierr)
!    call spl1df(x,zbup,xaa,uz,nround,ierr)

    ! find point B in lower half region by the method chosen among (1) and (2)

    ! (1) use Newton method
    call newton1df(xaa(ndw-1),x,spl1dd,xaa,uz,spl1dd,ur,nround,tol,linZ,rzmin,zmid0,rpmax,zzmin)
    call spl1df(x,rbdw,xaa,ur,nround,ierr)
    call spl1df(x,zbdw,xaa,uz,nround,ierr)

    ! (2) use bisection instead of newton1df (very slow)
!    call bisectionf(linZ,x,spl1df,xaa,uz,spl1df,ur,nround,xaa(ndw-1),xaa(ndw),tol,rzmin,zmid0,rpmax,zzmin)
!    call spl1df(x,rbdw,xaa,ur,nround,ierr)
!    call spl1df(x,zbdw,xaa,uz,nround,ierr)
  
    deallocate(xaa,yar,yaz,dyar,dyaz,ur,uz)
    ! ===========================================================================

    dabup = distance(rzmax,zmid0,rbup,zbup) ! AB in upper region
    dabdw = distance(rzmin,zmid0,rbdw,zbdw) ! AB in lower region

    zetaup = ( dabup - dacup ) / dcdup
    zetadw = ( dabdw - dacdw ) / dcddw
    v%zeta   = 0.5d0 * ( zetaup + zetadw )
!!$    write(6,'(4ES15.7)') dabup,dacup,dcdup,zetaup
!!$    write(6,'(4ES15.7)') dabdw,dacdw,dcddw,zetadw
!!$    write(6,'(4ES15.7)') v%rminor,zetaup,zetadw,v%zeta

  end subroutine magtrace

  !     ***** calculate Z at R="1st argument" by linear interpolation *****

  real(8) function linZ(r,r1,z1,r2,z2)
    real(8), intent(in) :: r, r1, z1, r2, z2

    linZ = (z2-z1)/(r2-r1)*(r-r1)+z1
  end function linZ

  !     ***** calculate distance between (r1,z1) and (r2,z2) *****

  real(8) function distance(r1,z1,r2,z2)
    real(8), intent(in) :: r1, z1, r2, z2

    distance = sqrt((r2-r1)**2+(z2-z1)**2)
    
  end function distance

  !     ***** calculate intersection by linear interpolation *****
!!$  subroutine intersection(r,z,r1,z1,r2,z2,rgeo,zgeo,grad)
!!$    real(8), intent(in)  :: r1,z1,r2,z2,rgeo,zgeo,grad
!!$    real(8), intent(out) :: r,z
!!$
!!$    r = (r2*z1-r1*z2-(r2-r1)*(zgeo-rgeo*grad))/((r2-r1)*grad-z2+z1)
!!$    z = r * grad + zgeo - rgeo * grad
!!$
!!$  end subroutine intersection

  !     ***** for Brent method *****
!  real(8) function zfunc(x)
!    real(8), intent(in) :: x
!
!    integer(4) :: ierr
!    real(8) :: rl, zl, f
!
!    call spl1df(x,rl,xaa,ur,nzfunc,ierr)
!    f = linZ(rl,r1,z1,r2,z2)
!    call spl1df(x,zl,xaa,uz,nzfunc,ierr)
!    if( ierr /= 0 ) stop 'zfunc error'
!    zfunc = f - zl
!
!  end function zfunc

  real(8) function dzfunc(x)
    use libspl1d, ONLY : spl1dd
    real(8), intent(in) :: x

    integer(4) :: ierr
    real(8) :: f, df

    call spl1dd(x,f,df,xaa,uz,nzfunc,ierr)
    if( ierr /= 0 ) stop 'dzfunc error'
    dzfunc = df

  end function dzfunc

  !     ***** for rtsafe method *****

  subroutine dzsub(x,df,ddf)
    use libspl1d, ONLY : spl1ddd
    real(8), intent(in) :: x
    real(8), intent(out) :: df, ddf

    integer(4) :: ierr
    real(8) :: f

    call spl1ddd(x,f,df,ddf,xaa,uz,nzfunc,ierr)
    if( ierr /= 0 ) stop 'dzsub error'

  end subroutine dzsub

  !     ***** INTEGRATE ALONG THE MAGNETIC FIELD LINE *****

  subroutine eqmags(rinit,zinit,nmax,xa,ya,n,ind)
    use equ_params, only : nw => nsr, nh => nsz, rg, zg, upsi, saxis, tol, zaxis, psign
    use subs, only : newtn
    use mod_num_recipe
    use libspl1d, ONLY : spl1d,spl1dd
    use libspl2d, ONLY : spl2df

!     ** Input **
!       RINIT : Initial starting point for tracing
!       ZINIT : Initial starting point for tracing
!       NMAX  : Size of arrays of XA, YA
!     ** Output **
!       XA    : Length from (RINIT,ZINIT) to the current position along the field line
!       YA    : Coordinate of the current position
!       N     : Number of partitions along the magnetic surface
!     ** Input/Output **
!       IND   : (Input) on LCFS if ind == -1 ; (Output) Error indicator

!toms760    use Grid_Interpolation

    integer(4), intent(in)  :: nmax
    real(8),    intent(in)  :: rinit, zinit
    real(8),    intent(out) :: xa(1:nmax), ya(1:2,1:nmax)
    integer(4), intent(out) :: n
    integer(4), intent(inout) :: ind

    integer(4) :: neq, istep, imode, i, j, nstart, istat, iallow_enlarge_stepsize, iosurf, ist
    integer(4) :: istepmax = 4, nimode = 1, idebug = 0
    real(8) :: rkap
    real(8) :: x, h, del, fact, perimeter, h_old, h_org!, psitmp(1), rgl(1), zgl(1)
    real(8), dimension(1:2) :: y, dydx, yout, yerr, yscal = 1.d0
    integer(4) :: nxtmp(1), ierr_newtn, ntmp, ierr_spl2df, ierr_spl1d
    real(8) :: rxtmp, zxtmp, rx, zx, rxmid, zxmid
    real(8) :: psiax, zdpsidr, zdpsidz, zpsi, zrbp, zrbp0
    real(8) :: vertical_criterion = 1.d-4
    real(8), allocatable :: dummy(:), psi1d(:)
    character(len=8) :: file_surf

    neq=2

    if( abs(irktrc) == 1 ) then
       itrace = 0
    else
       itrace = 1
    end if
    if( itrace > 0 ) then
       iallow_enlarge_stepsize = 1
    else
       iallow_enlarge_stepsize = 0
    end if

    ! itrace = 0: use rkck  (Runge-Kutta 4th & 5th)
    !        = 1: use eqrk4 (Runge-Kutta 4th)
    if( itrace == 0 ) then ! rkck
       rkap = 1.5d0
    else                   ! eqrk4
       rkap = 2.d0
    endif
    fact = sqrt(sqrt(2.d0))
    ! perimeter: rough estimate of the perimeter of an ellipse as a flux surface
    if( itrace == 0 ) then ! rkck
       perimeter = pi*(1.d0+rkap)*(rinit-rmaxis)*(1.d0+0.25d0*((1.d0-rkap)/(1.d0+rkap))**2)
       perimeter = 1.1d0 * perimeter ! for redundancy
    else                   ! eqrk4
       perimeter = 2.d0*pi*rkap*(rinit-rmaxis)
    endif
    h = fact * perimeter / nmax
    h_org = h
    istep=0

!      WRITE(6,'(I5,3ES12.4)') 0,H,rmaxis,zmaxis
!      WRITE(6,'(I5,3ES12.4)') NMAX,FACT,PI,RKAP

    lp: do
       x=0.d0
       y(1)=rinit
       y(2)=zinit
    
!       WRITE(6,'(I5,3ES12.4)') 1,X,Y(1),Y(2)

       n = 1
       xa(n)   = x
       ya(1,n) = y(1)
       ya(2,n) = y(2)
!toms760    call rgbi3p(1, nw, nh, rg, zg, psi, 1, y(1), y(2), psitmp(1), ierr)
!toms760    write(6,'(I5,I2,6ES12.4)') n,ierr,x,y(1),y(2),psig(y(1),y(2)),psitmp(1),sibry
!toms760
!toms760    i=28 ; j=15
!toms760    rgl(1) = rg(i)
!toms760    zgl(1) = zg(j)
!toms760    call rgbi3p(1, nw, nh, rg, zg, psi, 1, rgl(1), zgl(1), psitmp(1), ierr)
!toms760    write(6,*) rgl(1),zgl(1)
!toms760    write(6,*) psi(i+nw*(j-1)),psig(rgl(1),zgl(1)),psitmp(1)
!toms760    stop

       ! *** Trace the designated magnetic surface ***

       istat = 1
       ! imode = 0: first half around of a flux surface
       !       = 1: second half around
       imode = 0
       do i = 2, nmax
          call eqderv(x,y,dydx)
          if( itrace == 0 ) then
             call rkck(y,dydx,x,h,yout,yerr,eqderv)
          else
             call eqrk4(x,y,dydx,yout,h,neq,eqderv)
          endif
!          write(6,'(I4,4ES15.7)') i,yout(1),yout(2),yerr(1),yerr(2)
!          write(6,'(2I5,I2,7ES12.4)') i,n,imode,yout(1),yout(2),dydx(1),dydx(2),h,xa(n),x
          if(imode == 0) then ! first half of the surface
             if((yout(2)-zinit)*psign > 0.d0) then
                imode  = 1
                nimode = n
             endif
          else ! second half of the surface
             ! ===== Closed magnetic surface obtained ====
             if((yout(2)-zinit)*psign < 0.d0) then ! trace normal end
                istat = 0
                exit lp
             end if
             ! ===========================================
          endif
          x = x+h
          y(1) = yout(1)
          y(2) = yout(2)

!cont          WRITE(6,'(I5,5ES12.4)') N+1,X,Y(1),Y(2),DYDX(1),DYDX(2)

!toms760       ! ===
!toms760       call rgbi3p(2, nw, nh, rg, zg, psi, 1, y(1), y(2), psitmp(1), ierr)
!toms760       write(6,'(I5,I2,6ES12.4)') N+1,ierr,X,Y(1),Y(2),psig(y(1),y(2)),psitmp(1),sibry
!toms760       ! ===

          n = n+1
          xa(n)   = x
          ya(1,n) = y(1)
          ya(2,n) = y(2)
!          write(6,'(2I5,I2,7ES12.4)') i,n,imode,yout(1),yout(2),dydx(1),dydx(2),h,xa(n),x
       enddo

       ! **************************************************************    | 
       ! ***                                                        ***    |
       ! ***  The tracing points do not close the magnetic surface  ***  \ | /
       ! ***                                                        ***   \|/
       ! **************************************************************    V

       if( iallow_enlarge_stepsize == 1 ) then
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !  For the case that the step size may be too short to enclose the flux surface,
          !   although it is believed that tracing functioned well
          !
          !  Coping plan : ENLARGE THE STEP SIZE.
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          !  In many cases where the tracing fails, the tracing point typically moves vertically,
          !  i.e., in the Z direction.
          !  In such cases, |Z_{i} - Z_{i-1}|, which is the vertical distance at one step, 
          !  becomes quite larger than that of the horizonal distance, |R_{i} - R_{i-1}|.
          !  Try to detect such a situation by checking whether
          !     |( R_{nmax} - R_{nmax-1} )/( Z_{nmax} - Z_{nmax-1} )| 
          !  is small.
          !  If so, enlarge the step size "h" to jump and avoid the point at which the tracing 
          !  point goes the wrong way.

          if( istep <= istepmax ) then
             if(      abs((ya(1,nmax) - ya(1,nmax-1))/(ya(2,nmax) - ya(2,nmax-1))) &
                  & > vertical_criterion ) then
                h_old = h
                h=fact*h
                write(6,'(I2,2(A,E11.4),A)') istep, ": Enlarge the step size from ", h_old &
                     &                            , " to ", h, " and retrace."
                istep=istep+1
                cycle ! Go back and re-trace
             endif

             ! That you are here means that tracing fails (istat = 1) but the magnetic line
             !   may not diverge vertically.
             ! Then, you go to plan B in the following.
             
          else
             write(6,'(I2,A,I2)') istep, " reaches the limit of enlarging the step size, " &
                  &             , istepmax
             write(6,*) "Try to skip the point where some error occurs when tracing."
             ! Exit lp loop
          end if
       end if
       exit

    end do lp

    ! ---------------------------------------------------------------

    !  "ind == -1" means that tracing is executed on the LCFS.
    !  We have to find the X point at least once to determine lclockwise and lX_lownull.

    if( istat /= 0 .or. ind == -1 ) then

       if( ind == -1 ) then
          file_surf = 'lcfs.dat'
       else
          file_surf = 'surf.dat'
       end if

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  For the case that tracing itself does not work properly near the X point
    !
    !  Coping plan : SKIP THE POINT WHERE THE TRACING STARTS TO DIVERGE.
    !
    !  Note that we will go through until step 2. to determine lclockwise and lX_lownull
    !  even if istat == 0, meaning that tracing successfully finished.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    il_skip_X: do 
       !--- 1. Check whether tracing is clockwise or anti-clockwise

       ! Anti-clockwise tracing when saxis is negative (convex downward).
       !      Clockwise tracing when saxis is positive (convex upward).
       ! N.B. saxis is defined as "saxis = simag - sibry"
       !      and saxis is NOT identical to simag in principle.

       if( saxis <= 0.d0 ) then
          lclockwise = .false.
       else
          lclockwise = .true.
       end if

       !--- 2. Find the X point ---

       ! The final trace point (R,Z) at "n" is to large extent away from the actual surface.
       ! At first step, try to detect whether the problem occurred at the upper or lower side 
       !   of the magnetic surface by comparing Z at the final point and Z at the initial point.

       if( ya(2,n) > zinit ) then ! trace error occurs on the upper side
!!$          ! nxtmp is the location where R becomes minimum in the course of tracing.
!!$          nxtmp = minloc(ya(1,1:n))
          nstart = 1
       else                          ! trace error occurs on the lower side
!!$          ! nxtmp is the location where R becomes maximum in the course of tracing.
!!$          nxtmp = maxloc(ya(1,nimode:n))
          nstart = nimode
       endif

       ! Find the point, (rxtmp,zxtmp), which is supposed to be nearest a true X point,
       !  on the trail of the tracing,
       
       ! Find the nearest X point in the course of tracing by searching the minimum of RB_p
       ntmp = 0
       zrbp0 = 100.d0 ! arb. large value, RB_p
       do i = nstart, n
          call psigd(ya(1,i),ya(2,i),zdpsidr,zdpsidz,zpsi)
          zrbp = sqrt(zdpsidr**2 + zdpsidz**2) ! RB_p
!          write(6,'(I4,4ES15.7)') i,ya(1,i),ya(2,i),zrbp,zrbp0
          if( zrbp < zrbp0 ) then
             ntmp = i
             zrbp0 = zrbp
          end if
       end do
!!$       ntmp = nxtmp(1)
       rxtmp = ya(1,ntmp)
       zxtmp = ya(2,ntmp)
       if( ntmp == 1 ) stop 'ntmp should be greater than 1 due to numerical reason.'

       ! Find a true X point (rx, zx) on the 2D psi space

!       write(6,'(A/3ES15.7)') "Prior to NEWTN",rxtmp,zxtmp,psig(rxtmp,zxtmp)
       call newtn( psigd, rxtmp, zxtmp, rx, zx, 1.0d-8, 1.0d-8, 100, idebug, ierr_newtn )
       ! NOTE: rxtmp is equal to rx and zx is equal to zxtmp after newtn.
!       write(6,'(A/3ES15.7)') "Posterior to NEWTN",rx,zx,psig(rx,zx)
       if( ierr_newtn /= 0 ) then
          write(6,'(/X,2A)') ' Note: There may not exist an X point. ' // &
               & 'Please check the shape of the surface by looking into ', file_surf
          open(newunit=iosurf ,file=file_surf,iostat=ist,status='replace',form='formatted')
          do i = 1, n
             write(iosurf,'(3ES15.7)') xa(i), ya(1,i), ya(2,i)
          end do
          close(iosurf)
       end if
       if( zx < zinit ) then
          lX_lownull = .true.
       else
          lX_lownull = .false.
       end if

       ! Exit when "istat == 0", which means tracing was successful.
       ! With istat == 0, we come here only when we are on the LCFS.

       if( istat == 0 .and. ind == -1 ) then
          ind = 0
          exit il_skip_X
       end if

       !----- skipping the region adjacent to the X point -----

       ! renewed starting point, n+1, avoiding X point
       !   n : last point before the X point
!!       write(6,*) "     X",rx,zx,psig(rx,zx)
       if(     (ya(1,ntmp) - rx)*(ya(1,ntmp-1) - rx) >= 0.d0 .and. &
            &  (ya(2,ntmp) - zx)*(ya(2,ntmp-1) - zx) >= 0.d0 ) then
          n = ntmp
       else
          n = ntmp - 1
       endif
!!       write(6,*) "before",ya(1,n),ya(2,n),psig(ya(1,n),ya(2,n)),psig(ya(1,n-1),ya(2,n-1)),psig(rinit,zinit)

       !--- 3. Find (R,Z) on the magnetic surface across the X point

       ! To find R satisfying psi(R,Z) = psi(ya(1,n),Z) on the Z=ya(2,n) surface,
       ! make 1D psi(R) on Z=ya(2,n), i.e., psi(R,ya(2,n)).
       allocate(dummy(1:nw), psi1d(1:nw), upsi1d(1:4,1:nw))
       do i = 1, nw
          call spl2df(rg(i),ya(2,n),psi1d(i),rg,zg,upsi,nw,nw,nh,ierr_spl2df)
          if( ierr_spl2df /= 0 ) stop 'XX SPL2DF ERROR!'
       enddo
       psi1d(:) = psi1d(:) - psig(rinit,zinit)
!       psi1d(:) = psi1d(:) - psig(ya(1,n),ya(2,n))
       call spl1d(rg,psi1d,dummy,upsi1d,nw,0,ierr_spl1d)
       if( ierr_spl1d /= 0 ) stop 'XX SPL1D ERROR!'

       ! (rxtmp,zxtmp) is the point on the same magnetic surface across the X point,
       ! where zxtmp = ya(2,n)
       !  NOTE: +h and -h are added as numerical allowance.
       if     ( lX_lownull .eqv. .true.  .and. lclockwise .eqv. .true.  ) then
          rxtmp = rtsafe(funcpsir,ya(1,n)-h,rx,tol)
       else if( lX_lownull .eqv. .false. .and. lclockwise .eqv. .true.  ) then
          rxtmp = rtsafe(funcpsir,rx,ya(1,n)+h,tol)
       else if( lX_lownull .eqv. .true.  .and. lclockwise .eqv. .false. ) then
          rxtmp = rtsafe(funcpsir,rx,ya(1,n)+h,tol)
       else
          rxtmp = rtsafe(funcpsir,ya(1,n)-h,rx,tol)
       end if
!!$       write(6,*) rxtmp,rx,ya(1,n),ya(1,n)-h
!!$       write(6,*) psig(rxtmp,ya(2,n)),psig(ya(1,n),ya(2,n))
       deallocate(dummy,psi1d,upsi1d)
       zxtmp = ya(2,n)

       !--- 4. Find (R,Z) nearest the X point on the magnetic surface
       !       between the last point before the X point, (ya(1,n),ya(2,n))
       !           and the next point across the X point, (rxtmp  ,zxtmp  ).

       ! rxmid is defined as the center of rxtmp and ya(1,n)
       ! To find Z satisfying psi(rxmid,Z) = psi(rxmid,Z=ya(2,n)) on the Z line,
       ! make 1D psi(Z) on R=rxmid, i.e., psi(rxmid,Z)
       allocate(dummy(1:nh), psi1d(1:nh), upsi1d(1:4,1:nh))
       rxmid = 0.5d0 * (ya(1,n) + rxtmp) ! midpoint
       do i = 1, nh
          call spl2df(rxmid,zg(i),psi1d(i),rg,zg,upsi,nw,nw,nh,ierr_spl2df)
          if( ierr_spl2df /= 0 ) stop 'XX SPL2DF ERROR!'
       enddo
       psi1d(:) = psi1d(:) - psig(rinit,zinit)
!       psi1d(:) = psi1d(:) - psig(ya(1,n),ya(2,n))
       call spl1d(zg,psi1d,dummy,upsi1d,nh,0,ierr_spl1d)
       if( ierr_spl1d /= 0 ) stop 'XX SPL1D ERROR!'

       ! (rxmid,zxmid) is the point nearest the X point on the magnetic surface
       !    rxmid is situated at the center between ya(1,n) and zxtmp.
       !    zxmid is situated at the center between zx      and zxtmp.
       if( lX_lownull ) then
          zxmid = rtsafe(funcpsiz,zx,zxtmp+h,tol)
       else
          zxmid = rtsafe(funcpsiz,zxtmp-h,zx,tol)
       end if
       deallocate(dummy,psi1d,upsi1d)

       !--- 5. Replace the current position with which the trace error began
       !       by (rxmid,zxmid)

       x    = xa(n) + sqrt((rxmid - ya(1,n))**2 + (zxmid - ya(2,n))**2)
!       y(1) = rxmid
!       y(2) = zxmid

       ! Record
       n = n+1
       xa(n)   = x
       ya(1,n) = rxmid
       ya(2,n) = zxmid
!cont       WRITE(6,'(I5,5ES12.4)') N,X,Y(1),Y(2)!,DYDX(1),DYDX(2)

       !--- 6. Hard-code the next position by (rxtmp,zxtmp)

       ! Current position
       y(1) = rxtmp
       y(2) = zxtmp
       x    = xa(n) + sqrt((rxtmp - rxmid)**2 + (zxtmp - zxmid)**2)

       ! Record
       n = n+1
       xa(n)   = x
       ya(1,n) = y(1)
       ya(2,n) = y(2)
!!!       write(6,*) " after",ya(1,n),ya(2,n),psig(ya(1,n),ya(2,n))
!cont       WRITE(6,'(I5,5ES12.4)') N,X,Y(1),Y(2)!,DYDX(1),DYDX(2)
       !-------------------------------------------------------

       !--- 7. Restart tracing
       do i=n+1,nmax
          call eqderv(x,y,dydx)
          if( itrace == 0 ) then
             call rkck(y,dydx,x,h,yout,yerr,eqderv)
          else
             call eqrk4(x,y,dydx,yout,h,neq,eqderv)
          endif
!          write(6,'(2I5,I2,7ES12.4)') i,n,imode,yout(1),yout(2),dydx(1),dydx(2),h,xa(n),x
          if(imode == 0) then
             if((yout(2)-zinit)*psign > 0.d0) then
                imode  = 1
                nimode = n
             endif
          else
          ! ===== Closed magnetic surface obtained ====
             if((yout(2)-zinit)*psign < 0.d0) exit il_skip_X ! trace normal end
          ! ===========================================
          endif
          x = x+h
          y(1) = yout(1)
          y(2) = yout(2)

!cont          WRITE(6,'(I5,5ES12.4)') N+1,X,Y(1),Y(2),DYDX(1),DYDX(2)

          n = n+1
          xa(n)   = x
          ya(1,n) = y(1)
          ya(2,n) = y(2)
!          write(6,'(2I5,I2,7ES12.4)') i,n,imode,yout(1),yout(2),dydx(1),dydx(2),h,xa(n),x
       enddo

       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !   Some tracing error occurs when you are here.
       !   Maybe this is the case where an equilibrium has two null points (double-null)
       !   and tracing fails at the second X point.
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       write(6,*) 'XX EQMAGS: TRACE FAILS. MAYBE AN EQUILIBRIUM HAS TWO NULLS.'
       ind=1
       return

    end do il_skip_X

    end if

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  Reduce the step size by ten times
    !    to exactly get the magnetic line returned to the initial point
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    istat = 1
    j = i
    h=0.1d0*h
    do i = 1, 11
       call eqderv(x,y,dydx)
       if( itrace == 0 ) then
          call rkck(y,dydx,x,h,yout,yerr,eqderv)
       else
          call eqrk4(x,y,dydx,yout,h,neq,eqderv)
       endif
       ! ===== Closed magnetic surface obtained ====
       if((yout(2)-zinit)*psign < 0.d0) then ! trace normal end
          istat = 0
          exit
       end if
       ! ===========================================
       x=x+h
       y(1)=yout(1)
       y(2)=yout(2)
    enddo
    if( istat /= 0 ) then
       write(6,*) 'XX EQMAGS: UNEXPECTED BEHAVIOR'
       ind=2
       return
    end if

    ! *** Final adjustment to exactly get the magnetic line returned to the initial point ***

    del=(zinit-y(2))/(yout(2)-y(2))
    x=x+h*del
    y(1)=y(1)+(yout(1)-y(1))*del
    y(2)=y(2)+(yout(2)-y(2))*del
    n=n+1
    xa(n)=x
    ya(1,n)=y(1)
    ya(2,n)=y(2)
    ind=0

!    write(6,'(2I5,I2,2F12.8,ES15.7)') n,j,istep,rinit,xa(n),h

    return
  end subroutine eqmags

!     ***** DERIVATIVES *****

  subroutine eqderv(x,y,dydx)

    real(8), intent(in) :: x
    real(8), dimension(:), intent(in)  :: y
    real(8), dimension(:), intent(out) :: dydx
    real(8) :: psid, dpsidr, dpsidz

    call psigd(y(1),y(2),dpsidr,dpsidz)

    psid = sqrt(dpsidr**2+dpsidz**2)

    ! clockwise
    dydx(1) =-dpsidz / psid ! dR/dl_p
    dydx(2) = dpsidr / psid ! dZ/dl_p
!    write(6,'(5ES12.4)') x,y(1),y(2),dydx(1),dydx(2)
    return
  end subroutine eqderv

!     ***** INTERPOLATE SUBROUTINE DPSIDR,DPSIDZ(R,Z) *****

  subroutine psigd(r,z,dpsidr,dpsidz,zpsi)
    use equ_params, only : upsi, rg, zg, nw => nsr, nh => nsz
    use libspl2d, ONLY : spl2dd
    real(8), intent(in)  :: r, z
    real(8), intent(out) :: dpsidr, dpsidz
    real(8), intent(out), optional :: zpsi

    integer(4) :: ierr
    real(8) :: zpsil

    call spl2dd(r,z,zpsil,dpsidr,dpsidz,rg,zg,upsi,nw,nw,nh,ierr)
    if( present(zpsi) ) zpsi = zpsil
    if(ierr /= 0) then
       write(6,*) 'XX psigd: SPL2DD ERROR: IERR=',ierr
       write(6,'(A,2ES12.4)') '   R,Z=',r,z
    endif
    return
  end subroutine psigd

!     ***** INTERPOLATE FUNCTION OF PSI(R,Z) *****

  real(8) function psig(r,z)
    use equ_params, only : upsi, rg, zg, nw => nsr, nh => nsz
    use libspl2d, ONLY : spl2df
    real(8), intent(in) :: r, z

    integer(4) :: ierr
    real(8) :: psil

    call spl2df(r,z,psil,rg,zg,upsi,nw,nw,nh,ierr)
    if( ierr /= 0) then
       write(6,*) 'XX PSIG: SPL2DF ERROR: IERR=',IERR
       write(6,'(A,2ES12.4)') '   R,Z=',R,Z
    endif
    psig = psil
    return
  end function psig

  !     ***** for rtsafe method *****

  ! calculate psi and psi' for given R
  subroutine funcpsir(x,fval,fderiv)
    use equ_params, only : rg, nw => nsr
    use libspl1d, ONLY : spl1dd
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: fval, fderiv

    integer(4) :: ierr

    call spl1dd(x,fval,fderiv,rg,upsi1d,nw,ierr)
    if( ierr /= 0 ) stop 'ERROR: funcpsir in trace.f90'

  end subroutine funcpsir

  ! calculate psi and psi' for given Z
  subroutine funcpsiz(x,fval,fderiv)
    use equ_params, only : zg, nh => nsz
    use libspl1d, ONLY : spl1dd
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: fval, fderiv

    integer(4) :: ierr

    call spl1dd(x,fval,fderiv,zg,upsi1d,nh,ierr)
    if( ierr /= 0 ) stop 'ERROR: funcpsiz in trace.f90'

  end subroutine funcpsiz

!     ****** SIMPLE RUNGE-KUTTA METHOD ******

  subroutine eqrk4(x,y,dydx,yout,h,n,derivs)

    implicit none
    integer(4),              intent(in) :: n
    real(8),                 intent(in) :: x, h
    real(8), dimension(1:n), intent(in) :: dydx, y
    real(8), dimension(1:n), intent(out):: yout
    integer(4), parameter      :: nmax = 50
    integer(4)                 :: i
    real(8)                    :: h6, hh, xh
    real(8), dimension(1:nmax) :: dym, dyt, yt

    interface
       subroutine derivs(x,y,dydx)
         real(8), intent(in) :: x
         real(8), dimension(:), intent(in)  :: y
         real(8), dimension(:), intent(out) :: dydx
       end subroutine derivs
    end interface

    if( n > nmax ) then
       write(6,*) 'XX EQRK4: N (number of eqs) > NMAX:'
       write(6,*) '   N,NMAX=',N,NMAX
       stop
    endif

    HH=H*0.5D0
    H6=H/6.D0
    XH=X+HH
    do I=1,N
       YT(I)=Y(I)+HH*DYDX(I)
    enddo
    call derivs(XH,YT,DYT)
    do I=1,N
       YT(I)=Y(I)+HH*DYT(I)
    enddo
    call derivs(XH,YT,DYM)
    do I=1,N
       YT(I)=Y(I)+H*DYM(I)
       DYM(I)=DYT(I)+DYM(I)
    enddo
    call derivs(X+H,YT,DYT)
    do I=1,N
       YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
    enddo
    return
  end subroutine eqrk4

end module mod_trace
