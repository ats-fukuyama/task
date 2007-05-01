!     $Id$

      program femtest

      implicit none
      integer:: mwmax,mlmax
      complex(8),dimension(:,:),allocatable:: fma
      complex(8),dimension(:),allocatable:: fvb,fvx
      real(8),dimension(:),allocatable:: rho,f1,f2

      integer id,nrmax,mw,ml,nr,ierr

      call gsopen
!      id=1
      id=2
      nrmax=11

      call initialize(nrmax)

      call fem_calc(id,nrmax)

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
      enddo
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      write(6,*) '# ierr= ',ierr

      do nr=1,nrmax
         f1(nr)=real(fvx(2*(nr-1)+1))
         f2(nr)=real(fvx(2*(nr-1)+2))
      enddo
      do nr=1,nrmax
         write(6,'(I5,1P3E12.4)') nr,rho(nr),f1(nr),f2(nr)
      enddo

      call pages
      call femgr1d(0,rho,f1,nrmax,'@f1(rho)@')
      call pagee
      call pages
      call femgr1d(0,rho,f2,nrmax,'@f1(rho)@')
      call pagee

      call gsclos

      contains

!----- set profile -----

      subroutine initialize(nrmax)

      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,save:: nrmax_save=0
      real(8):: drho
      integer:: nr

      if(nrmax.ne.nrmax_save) then
         if(allocated(rho)) deallocate(rho)
         if(allocated(f1)) deallocate(f1)
         if(allocated(f2)) deallocate(f2)
         allocate(rho(nrmax))
         allocate(f1(nrmax))
         allocate(f2(nrmax))
      endif
      nrmax_save=nrmax

      drho=1.d0/(nrmax-1)
      do nr=1,nrmax
         rho(nr)=drho*(nr-1)
      enddo
      return
      end subroutine initialize

!----- allocate matrix and vectors -----

      subroutine fem_init(nrmax,nvmax)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: nvmax  ! numbre of variables
      integer,save:: mwmax_save=0,mlmax_save=0

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif
    
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(allocated(fma)) deallocate(fma)
         allocate(fma(mwmax,mlmax))
      endif

      if(mlmax.ne.mlmax_save) then
         if(allocated(fvb)) deallocate(fvb)
         if(allocated(fvx)) deallocate(fvx)
         allocate(fvb(mlmax))
         allocate(fvx(mlmax))
      endif
      mwmax_save=mwmax
      mlmax_save=mlmax
      end subroutine fem_init

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc(id,nrmax)

      use libfem_mod
      implicit none
      integer,intent(in):: id     ! type of model 
                                  !   1 : poisson eq (x)
                                  !   2 : poisson eq (r)
      integer,intent(in):: nrmax  ! number of points including end points
      integer:: nr,ml,mw,mc,nvmax,i,j,k
      real(8):: drho,val1,val2,val3,val4,temp
      real(8),dimension(4):: coef

      select case(id)
      case(1)
         nvmax=2
      case(2)
         nvmax=2
      end select

      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2
      write(6,'(A,3I5)') 'mlmax,mwmax,mc=',mlmax,mwmax,mc

      if(id.eq.1) then
         nvmax=2
         call fem_init(nrmax,nvmax)

         do nr=1,nrmax-1
            drho=rho(nr+1)-rho(nr)
            ml=2*(nr-1)+1

            fma(mc  ,ml  )=fma(mc  ,ml  )+table_hh(5,5)/drho
            fma(mc+1,ml  )=fma(mc+1,ml  )+table_hh(6,5)/drho*drho
            fma(mc+2,ml  )=fma(mc+2,ml  )+table_hh(7,5)/drho
            fma(mc+3,ml  )=fma(mc+3,ml  )+table_hh(8,5)/drho*drho

            fma(mc-1,ml+1)=fma(mc-1,ml+1)+table_hh(5,6)/drho*drho
            fma(mc  ,ml+1)=fma(mc  ,ml+1)+table_hh(6,6)/drho*drho*drho
            fma(mc+1,ml+1)=fma(mc+1,ml+1)+table_hh(7,6)/drho*drho
            fma(mc+2,ml+1)=fma(mc+2,ml+1)+table_hh(8,6)/drho*drho*drho

            fma(mc-2,ml+2)=fma(mc-2,ml+2)+table_hh(5,7)/drho
            fma(mc-1,ml+2)=fma(mc-1,ml+2)+table_hh(6,7)/drho*drho
            fma(mc  ,ml+2)=fma(mc  ,ml+2)+table_hh(7,7)/drho
            fma(mc+1,ml+2)=fma(mc+1,ml+2)+table_hh(8,7)/drho*drho

            fma(mc-3,ml+3)=fma(mc-3,ml+3)+table_hh(5,8)/drho*drho
            fma(mc-2,ml+3)=fma(mc-2,ml+3)+table_hh(6,8)/drho*drho*drho
            fma(mc-1,ml+3)=fma(mc-1,ml+3)+table_hh(7,8)/drho*drho
            fma(mc  ,ml+3)=fma(mc  ,ml+3)+table_hh(8,8)/drho*drho*drho
         enddo
         do mw=1,mwmax
            fma(mw,1) = 0.d0
         enddo
         fma(mc,1)=1.d0
         do mw=1,mwmax
            fma(mw,mlmax-1) = 0.d0
         enddo
         fma(mc,mlmax-1) = 1.d0

         fvb(1)=1.d0
      else if(id.eq.2) then
         nvmax=2
         call fem_init(nrmax,nvmax)

         do nr=1,nrmax-1
            drho=rho(nr+1)-rho(nr)
            ml=2*(nr-1)+1
            val1=(0.1d0+rho(nr)               )**2
            val2=(0.1d0+rho(nr)+drho/3.d0     )**2
            val3=(0.1d0+rho(nr)+drho*2.d0/3.d0)**2
            val4=(0.1d0+rho(nr)+drho          )**2
!            val1=1.d0
!            val2=1.d0
!            val3=1.d0
!            val4=1.d0
            coef(1)=val1
            coef(2)=0.5d0*(-11*val1+18*val2- 9*val3+ 2*val4)
            coef(3)=val4
            coef(4)=0.5d0*( -2*val1+ 9*val2-18*val3+11*val4)
            write(6,'(I5,1P4E12.4)') nr,(coef(i),i=1,4)
            
!            if(nr.eq.1) then
!            do j=5,8
!               do k=5,8
!                  temp=0.d0
!                  do i=1,4
!                     temp=temp+coef(i)*table_hhh(i,j,k)
!                  enddo
!                  write(6,'(2I5,1P4E12.4)') j,k,temp,
!     &                 table_hhh(1,j,k),table_hhh(3,j,k),table_hh(j,k)
!               enddo
!            enddo
!            endif

            do i=1,4
            fma(mc  ,ml  )=fma(mc  ,ml  )
     &                    +coef(i)*table_hhh(i,5,5)/drho
            fma(mc+1,ml  )=fma(mc+1,ml  )
     &                    +coef(i)*table_hhh(i,6,5)/drho*drho
            fma(mc+2,ml  )=fma(mc+2,ml  )
     &                    +coef(i)*table_hhh(i,7,5)/drho
            fma(mc+3,ml  )=fma(mc+3,ml  )
     &                    +coef(i)*table_hhh(i,8,5)/drho*drho

            fma(mc-1,ml+1)=fma(mc-1,ml+1)
     &                    +coef(i)*table_hhh(i,5,6)/drho*drho
            fma(mc  ,ml+1)=fma(mc  ,ml+1)
     &                    +coef(i)*table_hhh(i,6,6)/drho*drho*drho
            fma(mc+1,ml+1)=fma(mc+1,ml+1)
     &                    +coef(i)*table_hhh(i,7,6)/drho*drho
            fma(mc+2,ml+1)=fma(mc+2,ml+1)
     &                    +coef(i)*table_hhh(i,8,6)/drho*drho*drho

            fma(mc-2,ml+2)=fma(mc-2,ml+2)
     &                    +coef(i)*table_hhh(i,5,7)/drho
            fma(mc-1,ml+2)=fma(mc-1,ml+2)
     &                    +coef(i)*table_hhh(i,6,7)/drho*drho
            fma(mc  ,ml+2)=fma(mc  ,ml+2)
     &                    +coef(i)*table_hhh(i,7,7)/drho
            fma(mc+1,ml+2)=fma(mc+1,ml+2)
     &                    +coef(i)*table_hhh(i,8,7)/drho*drho

            fma(mc-3,ml+3)=fma(mc-3,ml+3)
     &                    +coef(i)*table_hhh(i,5,8)/drho*drho
            fma(mc-2,ml+3)=fma(mc-2,ml+3)
     &                    +coef(i)*table_hhh(i,6,8)/drho*drho*drho
            fma(mc-1,ml+3)=fma(mc-1,ml+3)
     &                    +coef(i)*table_hhh(i,7,8)/drho*drho
            fma(mc  ,ml+3)=fma(mc  ,ml+3)
     &                    +coef(i)*table_hhh(i,8,8)/drho*drho*drho
         enddo
         enddo
         do mw=1,mwmax
            fma(mw,1) = 0.d0
         enddo
         fma(mc,1)=1.d0
         do mw=1,mwmax
            fma(mw,mlmax-1) = 0.d0
         enddo
         fma(mc,mlmax-1) = 1.d0

         fvb(1)=1.d0
         fvb(mlmax-1)=1.d0/11.d0
      endif
      return
      end subroutine fem_calc

c
c=======================================================================
      subroutine femgr1d(ngp,x,y,nxmax,str)
c
      implicit none
      integer ngp,nxmax,nx
      real*8, dimension(nxmax) :: x,y
      real*4, dimension(:), allocatable :: gx,gy
      character str*(*)
      real*4 guclip
c
      allocate(gx(nxmax))
      allocate(gy(nxmax))
      do nx=1,nxmax
         gx(nx)=guclip(x(nx))
         gy(nx)=guclip(y(nx))
      end do
c
      call grf1d(ngp,gx,gy,nxmax,nxmax,1,str,0)
c
      deallocate(gy)
      deallocate(gx)
      return
      end subroutine femgr1d
c
c=======================================================================
      subroutine femgr2d(ngp,x,y,z,nxm,nxmax,nymax,str,ntype)
c
      implicit none
      integer ngp,nxm,nxmax,nymax,nx,ny,ntype
      real*8, dimension(nxmax) :: x
      real*8, dimension(nymax) :: y
      real*8, dimension(nxm,nymax) :: z
      real*4, dimension(:), allocatable :: gx,gy
      real*4, dimension(:,:), allocatable :: gz
      integer, dimension(:,:,:), allocatable :: ka
      character str*(*)
      real*4 gxmin,gxmax,gsxmin,gsxmax,gxstep,gxorg
      real*4 gymin,gymax,gsymin,gsymax,gystep,gyorg
      real*4 gzmin,gzmax,gszmin,gszmax,gzstep,gzorg
      real*4 gpx,gpy,gpratio,gpxc,gpyc
      real*4 glx,gly,glratio
      real*4 gp(4),gp3d(6)
      real*4 guclip
c
      allocate(gx(nxmax))
      allocate(gy(nymax))
      allocate(gz(nxmax,nymax))
      allocate(ka(8,nxmax,nymax))
      do nx=1,nxmax
         gx(nx)=guclip(x(nx))
      end do
      do ny=1,nymax
         gy(ny)=guclip(y(ny))
      end do
      do ny=1,nymax
         do nx=1,nxmax
            gz(nx,ny)=guclip(z(nx,ny))
         end do
      end do
C
      call grfut1(gx,nxmax,gxmin,gxmax)
      call grfut1(gy,nymax,gymin,gymax)
      call grfut2(gz,nxm,nxmax,nymax,gzmin,gzmax)
      if(gzmin*gzmax.lt.0.0) then
         gzmax=max(abs(gzmin),abs(gzmax))
         gzmin=-gzmax
      endif
c
      call grfut3(gxmin,gxmax,gsxmin,gsxmax,gxstep,gxorg)
      call grfut3(gymin,gymax,gsymin,gsymax,gystep,gyorg)
      call grfut3(gzmin,gzmax,gszmin,gszmax,gzstep,gzorg)
      gzstep=0.2*gzstep
c
      call grfut4(ngp,gp)
c
      if(ntype.eq.1.or.ntype.eq.2) then
         glx=gsxmax-gsxmin
         gly=gsymax-gsymin
         gpx=gp(2)-gp(1)
         gpy=gp(4)-gp(3)
c
         gpratio=gpy/gpx
         glratio=gly/glx
         if(glratio.gt.gpratio) then
            gpxc=0.5*(gp(1)+gp(2))
            gpx=gpy/glratio
            gp(1)=gpxc-0.5*gpx
            gp(2)=gpxc+0.5*gpx
         elseif(glratio.lt.gpratio) then
            gpyc=0.5*(gp(3)+gp(4))
            gpy=gpx*glratio
            gp(3)=gpyc-0.5*gpy
            gp(4)=gpyc+0.5*gpy
         endif
      endif
c
      if(ntype.eq.1) then
         call grf2dcx(gp,
     &                gsxmin,gsxmax,gxstep,gxorg,
     &                gsymin,gsymax,gystep,gyorg,
     &                gszmin,gszmax,gzstep,gzorg,
     &                gz,nxm,nxmax,nymax,str,ka,0)
c
      elseif(ntype.eq.2) then
         call grf2dax(gp,
     &                gsxmin,gsxmax,gxstep,gxorg,
     &                gsymin,gsymax,gystep,gyorg,
     &                gszmin,gszmax,gzstep,gzorg,
     &                gz,nxm,nxmax,nymax,str,0,0)
c
      elseif(ntype.eq.3) then
         gp3d(1)=10.0*1.5
         gp3d(2)=20.0*1.5
         gp3d(3)=10.0*1.5
C         gp3d(4)=-60.0
         gp3d(4)=-90.0
C         gp3d(5)=65.0
         gp3d(5)=30.0
         gp3d(6)=100.0
         call grf2dbx(gp,
     &                gsxmin,gsxmax,gxstep,gxorg,
     &                gsymin,gsymax,gystep,gyorg,
     &                gszmin,gszmax,gzstep,gzorg,
     &                gz,nxm,nxmax,nymax,
     &                str,gp3d)
      endif
c
      deallocate(ka)
      deallocate(gz)
      deallocate(gy)
      deallocate(gx)
      return
      end subroutine femgr2d

      end program femtest
