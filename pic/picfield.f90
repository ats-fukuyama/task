!  ***** TASK/PIC PREPARATION *****

MODULE picfield
  PRIVATE
  PUBLIC poissn,fftpic
 
CONTAINS

!***********************************************************************
    subroutine poissn(nx,ny,nxh1,rhof,phif,cform,ipssn)
!***********************************************************************
      implicit none
      complex(8), dimension(nxh1,ny) :: rhof, phif
      real(8), dimension(nxh1,ny) :: cform
      integer(4) :: nx, ny, nxh1, nyh, i, j, ipssn
      real(8) :: alxi, alyi, pi, twopi, am, an, amn2, afsp, afsp2

         afsp  = 1.5 !+++ size of finite-size-particle
         afsp2 = afsp * afsp 

      if( ipssn .eq. 0 ) then
         pi    = 3.14159265358979d0
         twopi = 2.d0 * pi 
         alxi   = 1.d0 / dble(nx)
         alyi   = 1.d0 / dble(ny)
         nyh    = ny / 2

         do j = 1, nyh+1
         do i = 1, nxh1
            am = twopi * dble(i-1) * alxi
            an = twopi * dble(j-1) * alyi
            amn2 = am*am + an*an
           
            if( i .eq. 1 .and. j .eq. 1 ) then
               cform(i,j) = 0.0
            else 
               cform(i,j) = 1.d0 / amn2 * exp( - amn2 * afsp2 ) 
            endif

            if( j .ge. 2 .and. j .le. nyh ) then
               cform(i,ny-j+2) = cform(i,j)
            endif
         end do
         end do

      else

         !----- solve poisson equation
         do j = 1, ny
         do i = 1, nxh1
            phif(i,j) =  cform(i,j) * rhof(i,j)
         end do
         end do 

      endif
      
    end subroutine poissn

!***********************************************************************
    subroutine fftpic(nx,ny,nxh1,nx1,ny1,a,af,awk,afwk,ifset)
!***********************************************************************
      implicit none
      include 'fftw3.f'
      real(8), dimension(nx1,ny1) :: a
      real(8), dimension(nx,ny) :: awk
      real(8) :: alx, aly 
      complex(8), dimension(nxh1,ny) :: af, afwk
      integer(4) :: nx, ny, nxh1, nx1, ny1, ifset, i, j
      !....integer(8), save :: FFTW_ESTIMATE
      integer(8), save :: plan1, plan2

      alx = dble(nx)
      aly = dble(ny)

      if( ifset .eq. 0 ) then
         !----- initialization of fourier transform
         call dfftw_plan_dft_r2c_2d(plan1,nx,ny,awk,afwk,FFTW_ESTIMATE)
         call dfftw_plan_dft_c2r_2d(plan2,nx,ny,afwk,awk,FFTW_ESTIMATE)

      elseif( ifset .eq. -1 ) then

         !----- fourier transform
         do j = 1, ny
         do i = 1, nx
            awk(i,j) = a(i,j)
         end do
         end do

         call dfftw_execute(plan1)

         do j = 1, ny
         do i = 1, nxh1
         af(i,j) = afwk(i,j) / ( alx * aly )
         end do
         end do

      else

         !----- inverse fourier transform

         do j = 1, ny
         do i = 1, nxh1
         afwk(i,j) = af(i,j)
         end do
         end do

         call dfftw_execute(plan2)

         do j = 1, ny
         do i = 1, nx
            a(i,j) = awk(i,j)
         end do
         end do

         do j = 1, ny
            a(nx1,j) = a(1,j)
         end do

         do i = 1, nx1
            a(i,ny1) = a(i,1)
         end do

      endif

    end subroutine fftpic

end MODULE picfield
