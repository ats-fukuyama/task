!  ***** TASK/PIC PREPARATION *****

MODULE picfield
  PRIVATE
  PUBLIC poissn,fftpic
 
CONTAINS

!***********************************************************************
    subroutine poissn(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)
!***********************************************************************
      implicit none
      complex(8), dimension(nxmaxh1,nymax) :: rhof, phif
      real(8), dimension(nxmaxh1,nymax) :: cform
      integer(4) :: nxmax, nymax, nxmaxh1
      integer(4) :: nx, ny, nymaxh,ipssn
      real(8) :: alxi, alyi, pi, twopi, am, an, amn2, afsp, afsp2

         afsp  = 1.5 !+++ size of finite-size-particle
         afsp2 = afsp * afsp 

      if( ipssn .eq. 0 ) then
         pi    = 3.14159265358979d0
         twopi = 2.d0 * pi 
         alxi   = 1.d0 / dble(nxmax)
         alyi   = 1.d0 / dble(nymax)
         nymaxh    = nymax / 2

         do ny = 1, nymaxh+1
         do nx = 1, nxmaxh1
            am = twopi * dble(nx-1) * alxi
            an = twopi * dble(ny-1) * alyi
            amn2 = am*am + an*an
           
            if( nx .eq. 1 .and. ny .eq. 1 ) then
               cform(nx,ny) = 0.0
            else 
               cform(nx,ny) = 1.d0 / amn2 * exp( - amn2 * afsp2 ) 
            endif

            if( ny .ge. 2 .and. ny .le. nymaxh ) then
               cform(nx,nymax-ny+2) = cform(nx,ny)
            endif
         end do
         end do

      else

         !----- solve poisson equation
         do ny = 1, nymax
         do nx = 1, nxmaxh1
            phif(nx,ny) = cform(nx,ny) * rhof(nx,ny)
         end do
         end do 

      endif
      
    end subroutine poissn

!***********************************************************************
    subroutine fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,a,af,awk,afwk,ifset)
!***********************************************************************
      implicit none
      include 'fftw3.f'
      real(8), dimension(nxmax1,nymax1) :: a
      real(8), dimension(nxmax,nymax) :: awk
      real(8) :: alx, aly 
      complex(8), dimension(nxmaxh1,nymax) :: af, afwk
      integer(4) :: nxmax, nymax, nxmaxh1, nxmax1, nymax1
      integer(4) :: ifset, nx, ny
      !....integer(8), save :: FFTW_ESTIMATE
      integer(8), save :: plan1, plan2

      alx = dble(nxmax)
      aly = dble(nymax)

      if( ifset .eq. 0 ) then
         !----- initialization of fourier transform
         call dfftw_plan_dft_r2c_2d(plan1,nxmax,nymax,awk,afwk,FFTW_ESTIMATE)
         call dfftw_plan_dft_c2r_2d(plan2,nxmax,nymax,afwk,awk,FFTW_ESTIMATE)

      elseif( ifset .eq. -1 ) then

         !----- fourier transform
         do ny = 1, nymax
         do nx = 1, nxmax
            awk(nx,ny) = a(nx,ny)
         end do
         end do

         call dfftw_execute(plan1)

         do ny = 1, nymax
         do nx = 1, nxmaxh1
         af(nx,ny) = afwk(nx,ny) / ( alx * aly )
         end do
         end do

      else

         !----- inverse fourier transform

         do ny = 1, nymax
         do nx = 1, nxmaxh1
         afwk(nx,ny) = af(nx,ny)
         end do
         end do

         call dfftw_execute(plan2)

         do ny = 1, nymax
         do nx = 1, nxmax
            a(nx,ny) = awk(nx,ny)
         end do
         end do

         do ny = 1, nymax
            a(nxmax1,ny) = a(1,ny)
         end do

         do nx = 1, nxmax1
            a(nx,nymax1) = a(nx,1)
         end do

      endif

    end subroutine fftpic

end MODULE picfield
