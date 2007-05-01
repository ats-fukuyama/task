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

