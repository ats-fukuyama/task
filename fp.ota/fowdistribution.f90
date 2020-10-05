module fowdistribution
  implicit none
  private
  public :: fow_distribution_u2I, fow_distribution_maxwellian_inCOM, fow_get_ra_from_psip

contains

  subroutine fow_distribution_maxwellian_inCOM(fm_I)
    ! calculate Maxwellian in COM space using the mean psi_p of orbits. fm_I is value on mesh points.
      use fowcomm
      use foworbitclassify
      use fpcomm
      use fpsub
      use fpwrite

      real(rkind),intent(out) :: fm_I(:,:,:,:)
      real(rkind),allocatable :: U_fm_u(:,:), fm_u(:,:,:), fm_ux(:)
      real(rkind) :: mean_psip, mean_ra, dr ! mean_ra [r/a]
      integer :: np, nth, nr, nsa, ierr, ir

      ierr = 0

      allocate(fm_u(npmax,nrmax,nsamax),U_fm_u(4,nrmax),fm_ux(nrmax))

      ! calculate Maxwellian in (p,r/a) space
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            fm_u(np,nr,nsa) = FPMXWL(PM(NP,NSA),NR,NSA)
          end do
        end do
      end do

      do nsa = 1, nsamax
        do nth = 1, nthmax
          do np = 1, nthmax
            do nr = 1, nrmax
              if ( forbitten(np,nth,nr,nsa,[0,0,0]) ) then
                fm_I(np,nth,nr,nsa) = 0.d0
              else
                ! calculate mean psip of an orbit
                mean_psip = sum(orbit_m(np,nth,nr,nsa)%psip)/size(orbit_m(np,nth,nr,nsa)%psip)
                
                ! calculate mean minor radius of an orbit
                call fow_get_ra_from_psip(mean_ra, mean_psip)

                ! calculate fm_I using linear interpolartion
                do ir = 2, nrmax
                  if ( mean_ra <= rm(ir) ) then
                    dr = rm(ir)-mean_ra
                    fm_I(np,nth,nr,nsa) = fm_u(np,ir,nsa)-dr*(fm_u(np,ir,nsa)-fm_u(np,ir-1,nsa))/(rm(ir)-rm(ir-1))
                    exit
                  else if( mean_ra >= rm(nrmax))then
                    fm_I(np,nth,nr,nsa) = fm_u(np,nrmax,nsa)
                    exit
                  else if( mean_ra <= rm(1) )then
                    dr = rm(1)-mean_ra
                    fm_I(np,nth,nr,nsa) = fm_u(np,1,nsa)-dr*(fm_u(np,1,nsa)-fm_u(np,2,nsa))/(rm(1)-rm(2))
                    exit
                  end if
                end do
              end if
            end do
          end do
        end do
      end do

      call fpcsv1D(fm_u(:,2,2)*rnfp0(2)*1.d20,"./csv/fm_u_center.csv")
      call fpcsv1D(fm_u(:,nrmax/2,2)*rnfp0(2)*1.d20,"./csv/fm_u_quarter.csv")
      call fpcsv1D(fm_u(:,nrmax,2)*rnfp0(2)*1.d20,"./csv/fm_u_edge.csv")
      call fpcsv1D(fm_u(1,:,2)*rnfp0(2)*1.d20,"./csv/fm.csv")
      

      ! do nsa = 1, nsamax
      !   do nth = 1, nthmax
      !     do np = 1, nthmax
      !       call SPL1D(rm,fm_u(np,:,nsa),fm_ux,U_fm_u,NRMAX,0,IERR)
      !       do nr = 1, nrmax
      !         if ( forbitten(np,nth,nr,nsa,[0,0,0]) ) then
      !           fm_I(np,nth,nr,nsa) = 0.d0
      !         else
      !           ! calculate mean psip of an orbit
      !           mean_psip = sum(orbit_m(np,nth,nr,nsa)%psip)/size(orbit_m(np,nth,nr,nsa)%psip)
                
      !           ! calculate mean minor radius of an orbit
      !           call fow_get_ra_from_psip(mean_ra, mean_psip)

      !           ! calculate fm_I using spline interpolation
      !           do ir = 2, nrmax
      !             if ( mean_ra <= rm(ir) ) then
      !               dr = mean_ra-rm(ir-1)
      !               fm_I(np,nth,nr,nsa) = U_fm_u(1,ir)&
      !                                   +U_fm_u(2,ir)*dr&
      !                                   +U_fm_u(3,ir)*dr**2&
      !                                   +U_fm_u(4,ir)*dr**3
      !               exit
      !             else if( mean_ra >= rm(nrmax))then
      !               dr = mean_ra-rm(nrmax-1)
      !               fm_I(np,nth,nr,nsa) = U_fm_u(1,nrmax)&
      !                                   +U_fm_u(2,nrmax)*dr&
      !                                   +U_fm_u(3,nrmax)*dr**2&
      !                                   +U_fm_u(4,nrmax)*dr**3
      !               exit
      !             else if( mean_ra <= rm(1) )then
      !               dr = rm(1)-mean_ra
      !               fm_I(np,nth,nr,nsa) = U_fm_u(1,2)&
      !                                   +U_fm_u(2,2)*dr&
      !                                   +U_fm_u(3,2)*dr**2&
      !                                   +U_fm_u(4,2)*dr**3
      !               exit
      !             end if
      !           end do
      !         end if
      !       end do
      !     end do
      !   end do
      ! end do

  end subroutine fow_distribution_maxwellian_inCOM

  subroutine fow_get_ra_from_psip(ra_out, psip_in)
    use fowcomm
    use fpcomm

    real(rkind),intent(out) :: ra_out ! [r/a]
    real(rkind),intent(in) :: psip_in
    real(rkind),allocatable :: fx(:),U(:,:)
    real(rkind) :: dr
    integer :: nr, ierr

    ierr = 0
    DELR=(RMAX-RMIN)/NRMAX

    allocate(U(4,nrmax+1),fx(nrmax+1))

    call SPL1D(rg,psimg,fx,U,NRMAX+1,0,IERR)

    do nr = 2, nrmax+1
      if ( psimg(nr) >= psip_in ) then
        call solve_spl1D(U(:,nr),psip_in,dr)
        ra_out = rg(nr-1)+dr
        exit
      end if
    end do

  end subroutine fow_get_ra_from_psip

  subroutine solve_spl1D(U,a,x)
    ! U4*x**3+U3*x**2+U2*x+U1=a
    implicit none
    real(8),intent(in) :: U(4),a
    real(8),intent(inout):: x
    integer :: i=0,j
    real(8) :: x0,e=1.d10,eps=1.d-12,V(4),d

    V=U
    V(1)=V(1)-a

    do
      d=(V(4)*x**3+V(3)*x**2+V(2)*x+V(1))/(3.d0*V(4)*x**2+2.d0*V(3)*x+V(2))
      if(abs(d)<eps)exit

      x0=x
      x=x0-d

      if(i>e)then
        write(*,*)"Do not converge at subroutine solve_spl1D"
        stop
      end if
    end do

  end subroutine solve_spl1D

  ! subroutine fow_distribution_maxwellian_inCOM(fm_I)
  !   ! calculate Maxwellian in COM space. fm_I is value on mesh points.
  !   use fowcomm
  !   use foworbitclassify
  !   use fpcomm
  !   use fpsub

  !   real(rkind),intent(out) :: fm_I(:,:,:,:)
  !   integer :: np, nth, nr, nsa

  !   do nsa = 1, nsamax
  !     do nr = 1, nrmax
  !       do nth = 1, nthmax
  !         do np = 1, npmax
  !           if ( forbitten(np, nth, nr, nsa, [0,0,0]) ) then
  !             fm_I(np,nth,nr,nsa) = 0.d0
  !           else
  !             fm_I(np,nth,nr,nsa)=FPMXWL(PM(NP,NSA),NR,NSA)
  !           end if
  !         end do
  !       end do
  !     end do
  !   end do
  
  ! end subroutine fow_distribution_maxwellian_inCOM

  subroutine fow_distribution_u2I(fI_out, fu_in, orbit_in)
    ! calculate fI(np,nth,nr,nsa), distribution function in COM space, from fu(np,nth,nr,nsa)
    use fowcomm
    use fpcomm
    use foworbitclassify
    use foworbit
    
    real(rkind),intent(in) :: fu_in(:,:,:,:)
    real(rkind),intent(out) :: fI_out(:,:,:,:)
    type(orbit),intent(in) :: orbit_in(:,:,:,:)
    integer :: np,nth,nr,nsa,mode(3),ierr,nstp,nloop
    real(rkind) :: jacobian_u, pperp, ppara, pabs, cospitch, sinpitch, lorentzfact, pitch, rdot, dt
    real(rkind),allocatable :: jacobian_I(:,:,:,:),momentum(:),pitchAngle(:),minorRadius(:),FX(:,:,:),FY(:,:,:)&
                              ,FZ(:,:,:),FXY(:,:,:),FYZ(:,:,:),FZX(:,:,:),FXYZ(:,:,:),U(:,:,:,:,:,:)

    ierr = 0

    call fow_distribution_get_mode(mode, orbit_in)

    ! calculate J_I from orbit_in
    allocate(jacobian_I(npmax+mode(1),nthmax+mode(2),nrmax+mode(3),nsamax))
    call fow_orbit_jacobian(jacobian_I, orbit_in)

    allocate(momentum(npmax),pitchAngle(nthmax),minorRadius(nrmax),FX(npmax,nthmax,nrmax),FY(npmax,nthmax,nrmax)&
    ,FZ(npmax,nthmax,nrmax),FXY(npmax,nthmax,nrmax),FYZ(npmax,nthmax,nrmax)&
    ,FZX(npmax,nthmax,nrmax),FXYZ(npmax,nthmax,nrmax),U(4,4,4,npmax,nthmax,nrmax))

    ! prepare for SPL3D
    do nth = 1, nthmax
      pitchAngle(nth) = (nth*1.d0-0.5d0)/nthmax*pi
    end do
    do nr = 1, nrmax
      minorRadius(nr) = (nr*1.d0-0.5d0)/nrmax*ra
    end do

    do nsa = 1, nsamax
      ! prepare for SPL3D
      do np = 1, npmax
        momentum(np) = pm(np,nsa)
      end do

      ! calculate spline coefficient for f(np,nth,nr,nsa)
      call SPL3D(momentum,pitchAngle,minorRadius,fu_in(:,:,:,nsa),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,U,&
                npmax,nthmax,npmax,nthmax,nrmax,0,0,0,IERR)
      
      ! execute circular integrations over every orbit
      do nr = 1, nrmax+mode(3)
        do nth = 1, nthmax+mode(2)
          do np = 1, npmax+mode(1)           
            fI_out(np,nth,nr,nsa) = 0.d0

            if ( .not.forbitten(np,nth,nr,nsa,mode) ) then
              nloop = orbit_in(np,nth,nr,nsa)%nstp_max-1
              do nstp = 1, nloop
                lorentzfact = 1.d0&
                /sqrt(1.d0-(((orbit_in(np,nth,nr,nsa)%vpara(nstp))**2+(orbit_in(np,nth,nr,nsa)%vperp(nstp))**2)/vc**2))
                ppara = amfp(nsa)*lorentzfact*orbit_in(np,nth,nr,nsa)%vpara(nstp)
                pperp = amfp(nsa)*lorentzfact*orbit_in(np,nth,nr,nsa)%vperp(nstp)
                pabs = sqrt(ppara**2+pperp**2)
                cospitch = ppara/pabs
                sinpitch = pperp/pabs
                pitch = acos(cospitch)
                jacobian_u = 4*pi*RR*pabs**2*sinpitch
                rdot = sqrt((orbit_in(np,nth,nr,nsa)%psip(nstp+1)-orbit_in(np,nth,nr,nsa)%psip(nstp))**2&
                      +orbit_in(np,nth,nr,nsa)%psip(nstp)**2&
                      *(orbit_in(np,nth,nr,nsa)%theta(nstp+1)-orbit_in(np,nth,nr,nsa)%theta(nstp))**2)
                dt = orbit_in(np,nth,nr,nsa)%time(nstp+1)-orbit_in(np,nth,nr,nsa)%time(nstp)

                fI_out(np,nth,nr,nsa) = fI_out(np,nth,nr,nsa)&
                                        +jacobian_u&
                                        *fow_distribution_spl_fu(U,momentum,pitchAngle,minorRadius&
                                        ,pabs,pitch,orbit_in(np,nth,nr,nsa)%rs(nstp))&
                                        *rdot&
                                        *dt
              end do
              fI_out(np,nth,nr,nsa) = fI_out(np,nth,nr,nsa)/jacobian_I(np,nth,nr,nsa)
            end if

          end do
        end do
      end do
    end do

  end subroutine fow_distribution_u2I

  function fow_distribution_spl_fu(U_in,momentum,pitchAngle,minorRadius,p_in,pitch_in,ra_in) result(fu)
    use fpcomm
    use fowcomm

    real(rkind) :: fu
    real(rkind),intent(in) :: U_in(:,:,:,:,:,:),momentum(:),pitchAngle(:),minorRadius(:),p_in,pitch_in,ra_in
    real(rkind) :: deltap, deltath, deltar
    integer :: i,j,k,np,nth,nr,ip,ith,ir

    do np = 1, npmax
      if ( p_in<=momentum(1) ) then
        ip = 2
        deltap = p_in-momentum(1)
        exit
      else if ( p_in<=momentum(np) ) then
        ip = np
        deltap = p_in-momentum(np-1)
        exit
      else if ( p_in>momentum(npmax) )then
        ip = npmax
        deltap = p_in-momentum(npmax-1)
        exit
      end if
    end do

    do nth = 1, nthmax
      if ( pitch_in<=pitchAngle(1) ) then
        ith = 2
        deltath = pitch_in-pitchAngle(1)
        exit
      else if ( pitch_in<=pitchAngle(nth) ) then
        ith = nth
        deltath = pitch_in-pitchAngle(nth-1)
        exit
      else if ( p_in>pitchAngle(nthmax) )then
        ith = nthmax
        deltath = pitch_in-pitchAngle(nthmax-1)
        exit
      end if
    end do

    do nr = 1, nrmax
      if ( ra_in<=minorRadius(1) ) then
        ir = 2
        deltar = ra_in-minorRadius(1)
        exit
      else if ( ra_in<=minorRadius(nr) ) then
        ir = nr
        deltar = ra_in-minorRadius(nr-1)
        exit
      else if ( ra_in>minorRadius(nrmax) )then
        ir = nrmax
        deltar = ra_in-minorRadius(nrmax-1)
        exit
      end if
    end do

    fu = 0.d0

    do i = 1, 4
      do j = 1, 4
        do k = 1, 4
          fu = fu+U_in(i,j,k,ip,ith,ir)*deltap**(i-1)*deltath**(j-1)*deltar**(k-1)
        end do
      end do
    end do

  end function fow_distribution_spl_fu

  subroutine fow_distribution_get_mode(mode, orbit_in)
    use fowcomm 
    use fpcomm

    integer, intent(out) :: mode(3)
    type(orbit),intent(in) :: orbit_in(:,:,:,:)
    integer :: npm, nthm, nrm

    npm = size(orbit_in,1)
    nthm = size(orbit_in,2)
    nrm = size(orbit_in,3)

    if ( npm == npmax ) then
      mode(1) = 0
    else if ( npm == npmax+1 )then
      mode(1) = 1
    else
      write(*,*)"ERROR : size(orbit_in,1) is not npm = npmax or npmax+1 in fow_distribution_get_mode"
    end if

    if ( nthm == nthmax ) then
      mode(2) = 0
    else if ( nthm == nthmax+1 )then
      mode(2) = 1
    else
      write(*,*)"ERROR : size(orbit_in,2) is not nthm = nthmax or nthmax+1 in fow_distribution_get_mode"
    end if

    if ( nrm == nrmax ) then
      mode(3) = 0
    else if ( nrm == nrmax+1 )then
      mode(3) = 1
    else
      write(*,*)"ERROR : size(orbit_in,3) is not nrm = nrmax or nrmax+1 in fow_distribution_get_mode"
    end if

  end subroutine fow_distribution_get_mode

end module fowdistribution