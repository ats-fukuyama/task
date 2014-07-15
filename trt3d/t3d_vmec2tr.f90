  subroutine t3d_vmec2tr(file_name)
    use trcomm, only : dvrho, dvrhog, ar1rho, ar1rhog, ar2rho, ar2rhog, rm, rg, nrmax
    implicit none
    character(len=*), intent(in) :: file_name
    integer(kind=4), parameter :: rkind=8
    integer(kind=4), parameter :: ikind=4
    integer(kind=4), parameter :: dp=rkind
    
    real(kind=rkind),   parameter :: zero=0.0_dp, quat=0.25_dp, half=0.5_dp, one=1.0_dp, two=2.0_dp, c1p5=1.5_dp
    real(kind=rkind),   parameter :: pi=3.141592653589793_dp, twopi=pi+pi
    real(kind=rkind),   parameter :: dmu0=2.0e-7_dp*twopi                    ! magnetic permeability in vacuum [H/m]=[N/A^2]

    !... wout file name
    character(len=80)             :: wout_file_name ! name of wout file from vmec (input)

    !... data from wout_file 
    character(len=4)              :: version        ! version of vmec2000 ; accept '6.90' or '8.40' at present
    integer(kind=ikind)           :: ns             ! the number of radial grids ; j=0:ns
    integer(kind=ikind)           :: nsp1           ! = ns + 1
    integer(kind=ikind)           :: nsm1           ! = ns - 1
    integer(kind=ikind)           :: nfp            ! the number of toroidal field periods
    integer(kind=ikind)           :: mpol           ! the number of poloidal modes ;  0 <= m <= mpol-1
    integer(kind=ikind)           :: ntor           ! the number of toroidal modes ;  -ntor <= n <= ntor
    integer(kind=ikind)           :: mnmax          ! the number of Fourier modes for R, Z, and lambda
    integer(kind=ikind)           :: mnmax2         ! the number of Fourier modes for B and gsqrt ( =mnmax for version 6.9)
    integer(kind=ikind)           :: isgn           ! sign of  gsqrth in vmec2000 (left-hand coordinates -> isgn=-1)
    real(kind=rkind), allocatable :: xm(:)          ! poloidal mode number for R, Z, and lambda ; j=1:mnmax
    real(kind=rkind), allocatable :: xn(:)          ! toroidal mode number for R, Z, and lambda ; j=1:mnmax
    real(kind=rkind), allocatable :: xm2(:)         ! poloidal mode number for B and gsqrt      ; j=1:mnmax2
    real(kind=rkind), allocatable :: xn2(:)         ! toroidal mode number for B and gsqrt      ; j=1:mnmax2
    real(kind=rkind), allocatable :: rmnc(:,:)      ! Rmn on full grids        ; cos,  (1:mnmax,0:ns)
    real(kind=rkind), allocatable :: zmns(:,:)      ! Zmn on full grids        ; sin,  (1:mnmax,0:ns)
    real(kind=rkind), allocatable :: lmnsh(:,:)     ! lambda_mn on half grids  ; sin,  (1:mnmax,0:ns)
    real(kind=rkind), allocatable :: bmnh(:,:)      ! Bmn on half grids        ; cos,  (1:mnmax2,0:ns)
    real(kind=rkind), allocatable :: gmnh(:,:)      ! gsqrt_mn on half grids   ; cos,  (1:mnmax2,0:ns)
    real(kind=rkind), allocatable :: aiotah(:)      ! rotational transform on half grids ; j=0,ns+1
    real(kind=rkind), allocatable :: presh(:)       ! pressure on half grids [Pa]        ; j=0,ns+1
    real(kind=rkind), allocatable :: Itorh(:)       ! twopi*isgn*Itorh/mu0 = net toroidal current[A] on half grids ; j=0,ns+1
    real(kind=rkind), allocatable :: Ipolh(:)       ! twopi*isgn*Ipolh/mu0 = net poloidal current[A] on half grids ; j=0,ns+1
    real(kind=rkind), allocatable :: vprimh(:)      ! (dV/ds)/twopi**2 on half grids, where V is volume [m3] ; j=0,ns+1
    real(kind=rkind)              :: phiedge        ! (total toroidal flux at the edge)/(isgn*twopi) 
    
    !... equilibrium quantities from vmec
    real(kind=rkind), allocatable :: iota_vmec(:)   ! rotational transform from vmec2000 ; j=0,ns
    real(kind=rkind), allocatable :: pprim(:)       ! dp/ds from vmec2000, where p[Pa] ; j=0,ns
    real(kind=rkind), allocatable :: S11(:)         ! S11 on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: S12(:)         ! S12 on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: Bsqav(:)       ! <B^2> [T2] on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: vprim(:)       ! dV/ds [m3] on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: grdssq_av(:)   ! <|grad(s)|^2> on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: Itorf(:)       ! net toroidal current [A] on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: Ipolf(:)       ! net poloidal current [A] on a full mesh ; j=0,ns
    real(kind=rkind), allocatable :: s(:)           ! radial grid points (full mesh) = s ; s is normalized tor. flux ; j=0:ns

    character(len=80)   :: cdummy
    real(kind=rkind)    :: dummy
    integer(kind=ikind) :: idmy, iasym, nbsets, ntheta2, nzeta
    integer(kind=ikind) :: j, mn, mxm, nxn, l, k, iarg
    real(kind=rkind)    :: hs, dtheta, dzeta, arg, dnorm, ohs
    real(kind=rkind)    :: R, Z, R_t, Z_t, R_z, Z_z, l_th, l_zh, g_tt, g_tz, g_zz, l_t, l_z, gsqrth, Bh, gsqrt, Bf
    real(kind=rkind), allocatable :: xtheta(:), xzeta(:)
    real(kind=rkind), allocatable :: tsin(:,:,:), tcos(:,:,:)
    real(kind=rkind), allocatable :: tsin2(:,:,:), tcos2(:,:,:)
    real(kind=rkind), allocatable :: l_th_save(:,:),l_zh_save(:,:),gsqrth_save(:,:),Bh_save(:,:)
    
    real(kind=rkind), dimension(:), allocatable :: dvrho_vmec, ar2rho_vmec, rho_vmec

    !..............................................................................................................................
    !  open(8,file=trim(wout_file_name),form='formatted')
    open(8,file=file_name)
    !  if(ierr.ne.0) return
    
    read(8,'(a15,a4)') cdummy,version
    
    if(version(1:4) == '6.90') then
       read(8,*) (dummy,j=1,7)
       read(8,*) nfp,nsp1,mpol,ntor,mnmax,idmy,idmy,iasym,idmy,idmy
       mnmax2 = mnmax
       if(iasym == 1) then
          print*,'asymmetric equilibrium' ; stop
       endif
       read(8,*) idmy,idmy,idmy,idmy,idmy, idmy
       read(8,'(a)') cdummy
    else if(version(1:4) == '8.40') then
       read(8,*) (dummy,j=1,7)
       read(8,*) nfp,nsp1,mpol,ntor,mnmax,mnmax2,idmy,idmy,iasym,idmy,idmy
       if(iasym == 1) then
          print*,'asymmetric equilibrium' ; stop
       endif
       read(8,*) idmy,idmy,nbsets,idmy,idmy,idmy
       if(nbsets > 0) read(8,*) (dummy,j=1,nbsets)
       read(8,'(a)') cdummy
    else
       print*,'unsupported version : version =',version
       stop
    endif
    
    ns      = nsp1 - 1
    nsm1    = ns - 1
    ntheta2 = 4*mpol+4
    nzeta   = 4*ntor+4
    
    !... allocation of arrays
    allocate(xm(1:mnmax), xn(1:mnmax), xm2(1:mnmax2), xn2(1:mnmax2))
    allocate(rmnc(1:mnmax,0:ns), zmns(1:mnmax,0:ns))
    allocate(lmnsh(1:mnmax,0:ns), bmnh(1:mnmax2,0:ns), gmnh(1:mnmax2,0:ns))
    allocate(aiotah(0:nsp1), presh(0:nsp1), vprimh(0:nsp1), Itorh(0:nsp1), Itorf(0:ns))
    allocate(Ipolh(0:nsp1), Ipolf(0:ns))
    allocate(iota_vmec(0:ns))
    allocate(xtheta(1:ntheta2), xzeta(1:nzeta))
    allocate(tsin(1:mnmax,1:nzeta,1:ntheta2), tcos(1:mnmax,1:nzeta,1:ntheta2))
    allocate(tsin2(1:mnmax2,1:nzeta,1:ntheta2), tcos2(1:mnmax2,1:nzeta,1:ntheta2))
    allocate(pprim(0:ns))
    allocate(l_th_save(1:nzeta,1:ntheta2), l_zh_save(1:nzeta,1:ntheta2))
    allocate(gsqrth_save(1:nzeta,1:ntheta2), Bh_save(1:nzeta,1:ntheta2))
    allocate(S11(0:ns), S12(0:ns), Bsqav(0:ns), vprim(0:ns), grdssq_av(0:ns))
    allocate(s(0:ns))
    
    do j=0,ns
       if(version == '6.90') then
          do mn=1,mnmax                                           
             if(j==0) then
                read(8,*) mxm,nxn
                xm(mn)     = mxm ; xn(mn)     = nxn
                xm2(mn)    = mxm ; xn2(mn)    = nxn
             endif
             read(8,*) rmnc(mn,j),zmns(mn,j),                                            &        ! full mesh
                  &                    lmnsh(mn,j),bmnh(mn,j),gmnh(mn,j),                                &        ! half mesh
                  &                    dummy,dummy,                                                      &        ! half mesh
                  &                    dummy,                                                            &        ! full mesh    
                  &                    dummy,dummy,                                                      &        ! half mesh
                  &                    dummy                                                                      ! full mesh
          enddo
       else if(version == '8.40') then
          do mn=1,mnmax                                           
             if(j==0) then
                read(8,*) mxm,nxn
                xm(mn)     = mxm ; xn(mn)     = nxn
             endif
             read(8,*) rmnc(mn,j),zmns(mn,j),                                            &        ! full mesh
                  &                    lmnsh(mn,j)                                                                ! half mesh
          enddo
          do mn=1,mnmax2                                          
             if(j==0) then
                read(8,*) mxm,nxn
                xm2(mn) = mxm ; xn2(mn) = nxn
             endif
             read(8,*) bmnh(mn,j),gmnh(mn,j),                                            &        ! half mesh
                  &                    dummy,dummy,                                                      &        ! half mesh
                  &                    dummy,                                                            &        ! full mesh    
                  &                    dummy,dummy                                                                ! half mesh
          enddo
       endif
    enddo
    
    if(version == '6.90') then
       read(8,*) (aiotah(j),dummy,presh(j),dummy,phiedge,Itorh(j),Ipolh(j),dummy,vprimh(j),dummy,dummy,dummy,dummy,j=1,ns)
       read(8,*) (dummy,j=1,6)
       read(8,*) isgn
       do j=1,nsm1
          iota_vmec(j) = half*(aiotah(j)+aiotah(j+1))
       enddo
       aiotah(0)    = c1p5*aiotah(1) - half*aiotah(2)
       aiotah(nsp1) = c1p5*aiotah(ns) - half*aiotah(nsm1)
       iota_vmec(0)    = aiotah(0)
       iota_vmec(ns)   = aiotah(nsp1)
       
    else if(version == '8.40') then
       read(8,*) (iota_vmec(j),dummy,dummy,dummy,dummy,dummy, j=0,ns)                                           ! full
       read(8,*) (aiotah(j),dummy,presh(j),dummy,phiedge,Itorh(j),Ipolh(j),vprimh(j),dummy,dummy,j=1,ns)        ! half
       read(8,*) (dummy,j=1,6)
       read(8,*) isgn
       aiotah(0)    = iota_vmec(0)
       aiotah(nsp1) = iota_vmec(ns)
    endif
    
    Itorh(0)     = zero
    Itorh(nsp1)  = c1p5*Itorh(ns)  - half*Itorh(nsm1)
    Ipolh(nsp1)  = c1p5*Ipolh(ns)  - half*Ipolh(nsm1)
    
    do j=1,nsm1
       Itorf(j) = half*(Itorh(j)+Itorh(j+1))
       Ipolf(j) = half*(Ipolh(j)+Ipolh(j+1))
    enddo
    Itorf(0)  = Itorh(0)
    Itorf(ns) = Itorh(nsp1)
    Ipolf(0)  = Ipolh(1)
    Ipolf(ns) = Ipolh(nsp1)
    
    do mn=1,mnmax                                                 
       if(xm(mn) /= zero) then                                           
          rmnc(mn,0) = zero                                              
          zmns(mn,0) = zero                                              
       endif
    enddo
    
    !...................................................................................................................
    
    
    !... set up equilibrium parameters                                                   
    
    hs     = one/ns
    dtheta  = half/(ntheta2-1)                                            
    dzeta   = one/(nzeta*nfp)                                               
    
    do l = 1, ntheta2
       xtheta(l) = (l-1)*dtheta                                         
    enddo
    
    do k = 1, nzeta                                               
       xzeta(k) = dzeta*(k-1)                                          
    enddo
    
    do l = 1, ntheta2
       do k = 1, nzeta
          do mn = 1, mnmax
             arg  = xm(mn)*xtheta(l)-xn(mn)*xzeta(k)
             iarg = arg                                                     
             arg  = twopi * ( arg - iarg )                                  
             tsin(mn,k,l) = sin(arg)                                         
             tcos(mn,k,l) = cos(arg)
          enddo
       enddo
    enddo
    
    do l = 1, ntheta2
       do k = 1, nzeta
          do mn = 1, mnmax2
             arg  = xm2(mn)*xtheta(l)-xn2(mn)*xzeta(k)
             iarg = arg                                                     
             arg  = twopi * ( arg - iarg )                                  
             tsin2(mn,k,l) = sin(arg)                                         
             tcos2(mn,k,l) = cos(arg)
          enddo
       enddo
    enddo
    
    !... get p' on full grid points
    
    do j=1,nsm1
       pprim(j) = (presh(j+1)-presh(j))*ns      ! dp/ds[Pa], unit of presh is [Pa] (in case of vmec2000)
    enddo
    pprim(0)  = pprim(1)
    pprim(ns) = pprim(nsm1)
    !  pprim(ns) = zero                              ! here we assume p'=0 at plasma boundary
    
    dnorm = one/(nzeta*(ntheta2-1))
    ohs   = ns
    
    do l=1,ntheta2
       do k=1,nzeta
          l_th_save(k,l) = zero
          l_zh_save(k,l) = zero
          do mn = 1, mnmax
             l_th_save(k,l) = l_th_save(k,l) + xm(mn)*lmnsh(mn,1)*tcos(mn,k,l)
             l_zh_save(k,l) = l_zh_save(k,l) - xn(mn)*lmnsh(mn,1)*tcos(mn,k,l)
          enddo
          gsqrth_save(k,l) = zero
          Bh_save(k,l)     = zero
          do mn = 1, mnmax2
             gsqrth_save(k,l) = gsqrth_save(k,l) + gmnh(mn,1)*tcos2(mn,k,l)
             Bh_save(k,l)     = Bh_save(k,l)      + bmnh(mn,1)*tcos2(mn,k,l)
          enddo
       enddo
    enddo
    
    do j=1,nsm1
       S11(j) = zero ; S12(j) = zero ; Bsqav(j) = zero ; vprim(j) = zero ; grdssq_av(j) = zero
       do l=1,ntheta2
          do k=1,nzeta
             R = zero ; Z = zero ; R_t = zero ; Z_t = zero ; R_z = zero ; Z_z = zero ; l_th = zero ; l_zh = zero
             do mn = 1, mnmax
                R    = R    + rmnc(mn,j)*tcos(mn,k,l)
                Z    = Z    + zmns(mn,j)*tcos(mn,k,l)
                R_t  = R_t  - xm(mn)*rmnc(mn,j)*tsin(mn,k,l)
                R_z  = R_z  + xn(mn)*rmnc(mn,j)*tsin(mn,k,l)           
                Z_t  = Z_t  + xm(mn)*zmns(mn,j)*tcos(mn,k,l)           
                Z_z  = Z_z  - xn(mn)*zmns(mn,j)*tcos(mn,k,l) 
                l_th = l_th + xm(mn)*lmnsh(mn,j+1)*tcos(mn,k,l)           
                l_zh = l_zh - xn(mn)*lmnsh(mn,j+1)*tcos(mn,k,l) 
             enddo
             g_tt = R_t*R_t + Z_t*Z_t
             g_tz = R_t*R_z + Z_t*Z_z
             g_zz = R_z*R_z + Z_z*Z_z + R*R
             l_t  = half*(l_th_save(k,l)+l_th)
             l_z  = half*(l_zh_save(k,l)+l_zh)
             l_th_save(k,l) = l_th
             l_zh_save(k,l) = l_zh
             
             gsqrth = zero ; Bh = zero
             do mn = 1, mnmax2
                gsqrth = gsqrth + gmnh(mn,j+1)*tcos2(mn,k,l)
                Bh     = Bh     + bmnh(mn,j+1)*tcos2(mn,k,l)
             enddo
             gsqrt = half*(gsqrth_save(k,l)+gsqrth)
             Bf    = half*(Bh_save(k,l)+Bh)
             gsqrth_save(k,l) = gsqrth
             Bh_save(k,l)     = Bh
             
             if(l==1.or.l==ntheta2) then
                S11(j)    = S11(j) + half*g_tt/gsqrt
                S12(j)    = S12(j) + half*(g_tz*(1+l_t)-g_tt*l_z)/gsqrt
                Bsqav(j)  = Bsqav(j) + half*Bf*Bf*gsqrt
                vprim(j)  = vprim(j) + half*gsqrt
                grdssq_av(j) = grdssq_av(j) + half*(g_tt*g_zz-g_tz**2)/gsqrt
             else
                S11(j)    = S11(j) + g_tt/gsqrt
                S12(j)    = S12(j) + (g_tz*(1+l_t)-g_tt*l_z)/gsqrt
                Bsqav(j)  = Bsqav(j) + Bf*Bf*gsqrt
                vprim(j)  = vprim(j) + gsqrt
                grdssq_av(j) = grdssq_av(j) + (g_tt*g_zz-g_tz**2)/gsqrt
             endif
          enddo
       enddo
       S11(j)       = dnorm*S11(j)
       S12(j)       = dnorm*S12(j)
       Bsqav(j)     = Bsqav(j)/vprim(j)
       grdssq_av(j) = grdssq_av(j)/vprim(j)
       vprim(j)     = isgn*dnorm*vprim(j)
       
       !<<check>>
       !    print*,'vprim(j), (vprimh(j)+vprimh(j+1))/2 = ',vprim(j),half*(vprimh(j)+vprimh(j+1))
       
    enddo
    
    !... calculation for the plasma edge
    S11(ns) = zero
    S12(ns) = zero
    Bsqav(ns) = zero
    vprim(ns) = zero
    grdssq_av(ns) = zero
    do l=1,ntheta2
       do k=1,nzeta
          R = zero ; Z = zero ; R_t = zero ; Z_t = zero ; R_z = zero ; Z_z = zero ; l_th = zero ; l_zh = zero
          do mn = 1, mnmax
             R    = R    + rmnc(mn,ns)*tcos(mn,k,l)
             Z    = Z    + zmns(mn,ns)*tcos(mn,k,l)
             R_t  = R_t  - xm(mn)*rmnc(mn,ns)*tsin(mn,k,l)
             R_z  = R_z  + xn(mn)*rmnc(mn,ns)*tsin(mn,k,l)           
             Z_t  = Z_t  + xm(mn)*zmns(mn,ns)*tcos(mn,k,l)           
             Z_z  = Z_z  - xn(mn)*zmns(mn,ns)*tcos(mn,k,l) 
             l_th = l_th + xm(mn)*lmnsh(mn,ns-1)*tcos(mn,k,l)           
             l_zh = l_zh - xn(mn)*lmnsh(mn,ns-1)*tcos(mn,k,l) 
          enddo
          g_tt = R_t*R_t + Z_t*Z_t
          g_tz = R_t*R_z + Z_t*Z_z
          g_zz = R_z*R_z + Z_z*Z_z + R*R
          l_t  = c1p5*l_th_save(k,l) - half*l_th
          l_z  = c1p5*l_zh_save(k,l) - half*l_zh
          
          gsqrth = zero ; Bh = zero
          do mn = 1, mnmax2
             gsqrth = gsqrth + gmnh(mn,ns-1)*tcos2(mn,k,l)
             Bh     = Bh     + bmnh(mn,ns-1)*tcos2(mn,k,l)
          enddo
          gsqrt = c1p5*gsqrth_save(k,l) - half*gsqrth
          Bf    = c1p5*Bh_save(k,l)     - half*Bh
          if(l==1.or.l==ntheta2) then
             S11(ns)    = S11(ns) + half*g_tt/gsqrt
             S12(ns)    = S12(ns) + half*(g_tz*(1+l_t)-g_tt*l_z)/gsqrt
             Bsqav(ns)  = Bsqav(ns) + half*Bf*Bf*gsqrt
             vprim(ns)  = vprim(ns) + half*gsqrt
             grdssq_av(ns) = grdssq_av(ns) + half*(g_tt*g_zz-g_tz*g_tz)/gsqrt
          else
             S11(ns)    = S11(ns) + g_tt/gsqrt
             S12(ns)    = S12(ns) + (g_tz*(1+l_t)-g_tt*l_z)/gsqrt
             Bsqav(ns)  = Bsqav(ns) + Bf*Bf*gsqrt
             vprim(ns)  = vprim(ns) + gsqrt
             grdssq_av(ns) = grdssq_av(ns) + (g_tt*g_zz-g_tz*g_tz)/gsqrt
          endif
       enddo
    enddo
    S11(ns)   = dnorm*S11(ns)
    S12(ns)   = dnorm*S12(ns)
    Bsqav(ns) = Bsqav(ns)/vprim(ns)
    grdssq_av(ns) = grdssq_av(ns)/vprim(ns)
    vprim(ns) = isgn*dnorm*vprim(ns)
    
    S11(0)       = zero
    S12(0)       = zero
    Bsqav(0)     = two*Bsqav(1) - Bsqav(2)
    vprim(0)     = two*vprim(1) - vprim(2)
    grdssq_av(0) = zero
    
    do j = 1, ns
       s(j)     = j*hs                                          
    enddo
    s(0)      = zero  ; s(ns)        = one                                                 
    
    !... force iota_vmec to be (Itor(vmec)/phiedge-S12)/S11
    do j = 1, ns
       iota_vmec(j) = (Itorf(j)/phiedge-S12(j))/S11(j)
    enddo
    iota_vmec(0) = 2*iota_vmec(1) - iota_vmec(2)
    
    !... check consistency
    !      N.B. Unit of Itorf & Ipolf is not [A] and vprim=d(V/twopi**2)/ds at this point.
    
    print*, '* Phi_tor_edge = twopi*isgn*phiedge [Wb] =', twopi*isgn*phiedge
    print*, '* I_tor_edge (in vmec)  [kA]             =',0.001_dp*twopi*isgn/dmu0*Itorf(ns)
    print*
    print*,'      s           S11          S12           -S12/S11      iota(vmec)       Itor[kA]        Itor[kA]'
    print*,'                                                                       phip*(S11*iota+S12)   (vmec)'
    do j = 1,ns
       write(6,'(f10.4,2(e15.5),4(f15.8))') s(j),S11(j),S12(j),-S12(j)/S11(j),iota_vmec(j),&
            &           0.001_dp*twopi*isgn/dmu0*phiedge*(S11(j)*iota_vmec(j)+S12(j)),0.001_dp*twopi*isgn/dmu0*Itorf(j)
    enddo
    print*
    
    print*,'      s        pprim        pprim_check       <B^2>       <B^2>_check      vprim    <|grad(s)|^2>'
    do j = 1,nsm1
       write(6,'(f10.4,6(e15.5))') s(j),pprim(j),&
            &                            -isgn*phiedge*((Ipolh(j+1)-Ipolh(j))+(Itorh(j+1)-Itorh(j))*iota_vmec(j))/hs/dmu0/vprim(j),&
            &                             Bsqav(j),isgn*(Itorf(j)*iota_vmec(j)+Ipolf(j))*phiedge/vprim(j),vprim(j),grdssq_av(j)
    enddo
    print*
    
    !... physical value
    
    do j=0,ns
       Itorf(j) = twopi*isgn/dmu0*Itorf(j)  ! net toroidal current [A]
       Ipolf(j) = twopi*isgn/dmu0*Ipolf(j)  ! net poloidal current [A]
       vprim(j) = twopi**2*vprim(j)         ! dV/ds [m3]
    enddo
    
    deallocate(xtheta, xzeta)
    deallocate(tsin, tcos)
    deallocate(tsin2, tcos2)
    deallocate(l_th_save, l_zh_save)
    deallocate(gsqrth_save, Bh_save)
    
    close(8)
    
    
    !-----------------------------------------------------------------

    allocate(dvrho_vmec(ns), ar2rho_vmec(ns), rho_vmec(ns))
    do j = 1, ns
       dvrho_vmec(j) = vprim(j)*2.0d0*sqrt(s(j))
       ar2rho_vmec(j) = grdssq_av(j)*0.25d0/s(j)
       rho_vmec(j) = sqrt(s(j))
    end do
    
    call t3d_convert_mesh(ar2rhog, rg, nrmax+1, ar2rho_vmec, rho_vmec, ns)
    call t3d_convert_mesh(dvrhog, rg, nrmax+1, dvrho_vmec, rho_vmec, ns)
    
    call t3d_convert_mesh(ar2rho, rm, nrmax, ar2rhog, rg, nrmax+1)
    call t3d_convert_mesh(dvrho, rm, nrmax, dvrhog, rg, nrmax+1)
    
    ar1rhog(0) = sqrt(ar2rhog(0))
    do j = 1, nrmax
       ar1rho(j) = sqrt(ar2rho(j))
       ar1rhog(j) = sqrt(ar2rhog(j))
    end do
    
    deallocate(dvrho_vmec, ar2rho_vmec, rho_vmec)

  end subroutine t3d_vmec2tr


!-------------------------------------------------------------------------------

subroutine t3d_convert_mesh(y, x, imax, y0, x0, imax0)
  implicit none
  integer, intent(in) :: imax, imax0
  real(8), dimension(imax0), intent(in) :: y0, x0
  real(8), dimension(imax), intent(in) :: x
  real(8), dimension(imax), intent(out) :: y
  integer :: i, j
  real(8) :: dz, z

  do i = 1, imax
     z = x(i)/x0(imax0)
     if (x(i) < x0(1)) then
        call t3d_interpol_p(y(i),x(i),x0(1),x0(2),x0(3),y0(1),y0(2),y0(3))
     else if (x(i) > x0(imax0)) then
        dz = (z*x0(imax0)-x0(imax0-1))/(x0(imax0)-x0(imax0-1))
        y(i) = dz*y0(imax0) + (1.0d0-dz)*y0(imax0-1)
     else
        do j = 1, imax0
           if ((x0(j) <= z*x0(imax0)).and.(z*x0(imax0) <= x0(j+1))) then
              dz = (z*x0(imax0)-x0(j))/(x0(j+1)-x0(j))
              y(i) = dz*y0(j+1) + (1.0d0-dz)*y0(j)
              exit
           end if
        end do
     end if
  end do

end subroutine t3d_convert_mesh

!------------------------------------------------------------

subroutine t3d_interpol_p(y, x, x1, x2, x3, y1, y2, y3)
  implicit none
  real(8), intent(in) :: x, x1, x2, x3, y1, y2, y3
  real(8), intent(out) :: y
  
  y = (x-x2)*(x-x3)*y1/((x1-x2)*(x1-x3)) &
       + (x-x1)*(x-x3)*y2/((x2-x1)*(x2-x3)) &
       + (x-x1)*(x-x2)*y3/((x3-x1)*(x3-x2))
  
end subroutine t3d_interpol_p
