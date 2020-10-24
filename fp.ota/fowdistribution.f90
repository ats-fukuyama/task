module fowdistribution
  implicit none
  private
  public :: fow_distribution_maxwellian_inCOM, fow_get_ra_from_psip

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
        do np = 1, npmax
          do nth = 1, nthmax
            do nr = 1, nrmax
              if ( forbitten(nth,np,nr,nsa,[0,0,0]) ) then
                fm_I(nth,np,nr,nsa) = 0.d0
              else
                ! calculate mean psip of an orbit
                mean_psip = sum(orbit_m(nth,np,nr,nsa)%psip)/size(orbit_m(nth,np,nr,nsa)%psip)
                ! calculate mean minor radius of an orbit
                call fow_get_ra_from_psip(mean_ra, mean_psip)
                ! calculate fm_I using linear interpolartion
                do ir = 2, nrmax
                  if ( mean_ra <= rm(ir) ) then
                    dr = rm(ir)-mean_ra
                    fm_I(nth,np,nr,nsa) = fm_u(np,ir,nsa)-dr*(fm_u(np,ir,nsa)-fm_u(np,ir-1,nsa))/(rm(ir)-rm(ir-1))
                    exit
                  else if( mean_ra >= rm(nrmax))then
                    fm_I(nth,np,nr,nsa) = fm_u(np,nrmax,nsa)
                    exit
                  else if( mean_ra <= rm(1) )then
                    dr = rm(1)-mean_ra
                    fm_I(nth,np,nr,nsa) = fm_u(np,1,nsa)-dr*(fm_u(np,1,nsa)-fm_u(np,2,nsa))/(rm(1)-rm(2))
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
      !         if ( forbitten(nth,np,nr,nsa,[0,0,0]) ) then
      !           fm_I(nth,np,nr,nsa) = 0.d0
      !         else
      !           ! calculate mean psip of an orbit
      !           mean_psip = sum(orbit_m(nth,np,nr,nsa)%psip)/size(orbit_m(nth,np,nr,nsa)%psip)
                
      !           ! calculate mean minor radius of an orbit
      !           call fow_get_ra_from_psip(mean_ra, mean_psip)

      !           ! calculate fm_I using spline interpolation
      !           do ir = 2, nrmax
      !             if ( mean_ra <= rm(ir) ) then
      !               dr = mean_ra-rm(ir-1)
      !               fm_I(nth,np,nr,nsa) = U_fm_u(1,ir)&
      !                                   +U_fm_u(2,ir)*dr&
      !                                   +U_fm_u(3,ir)*dr**2&
      !                                   +U_fm_u(4,ir)*dr**3
      !               exit
      !             else if( mean_ra >= rm(nrmax))then
      !               dr = mean_ra-rm(nrmax-1)
      !               fm_I(nth,np,nr,nsa) = U_fm_u(1,nrmax)&
      !                                   +U_fm_u(2,nrmax)*dr&
      !                                   +U_fm_u(3,nrmax)*dr**2&
      !                                   +U_fm_u(4,nrmax)*dr**3
      !               exit
      !             else if( mean_ra <= rm(1) )then
      !               dr = rm(1)-mean_ra
      !               fm_I(nth,np,nr,nsa) = U_fm_u(1,2)&
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
  !             fm_I(nth,np,nr,nsa) = 0.d0
  !           else
  !             fm_I(nth,np,nr,nsa)=FPMXWL(PM(NP,NSA),NR,NSA)
  !           end if
  !         end do
  !       end do
  !     end do
  !   end do
  
  ! end subroutine fow_distribution_maxwellian_inCOM


end module fowdistribution