c
      module eqgout_mod
      public
      contains
c
c=======================================================================
      subroutine eqgout
c=======================================================================
c     interface eqiulibrium <>transport                            JAERI
c          transport grid -> equilibrium grid
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use eqt_mod
      use trn_mod
      USE libchar
      implicit none
c
      integer id,i,idx
      character kid*4,kid1*1

C      write(6,'(A,1P2E12.4)') 'saxis=',saxis
C      DO i=1,nv
C         write(6,'(i5,1P6E12.4)')
C     &        i,vlv(i),siv(i),hiv(i),prv(i),pds(i),qqv(i)
C      END DO
      
    1 continue
      write(6,*) '## input graph kid: Ann,Bnn,Cnn,Dn,Pn,X'
      read(5,'(A4)',ERR=1,END=9000) kid
      kid1=kid(1:1)      
      call toupper(kid1)
C
      if(kid1.eq.'A') then
         read(kid(2:4),'(I3)',ERR=1,END=9000) id
         call pages
            call eqgouta(id,0)
         call pagee
      elseif(kid1.eq.'B') then
         read(kid(2:4),'(I3)',ERR=1,END=9000) id
         call pages
            idx=id/100
            id=mod(id,100)
            call eqgouta(idx*100+4*(id-1)+1,1)
            call eqgouta(idx*100+4*(id-1)+2,2)
            call eqgouta(idx*100+4*(id-1)+3,3)
            call eqgouta(idx*100+4*(id-1)+4,4)
         call pagee
      elseif(kid1.eq.'C') then
         read(kid(2:4),'(I3)',ERR=1,END=9000) id
         call pages
            call eqgouta(    id,1)
            call eqgouta(100+id,2)
            call eqgouta(200+id,3)
            call eqgouta(300+id,4)
         call pagee
      elseif(kid1.eq.'D') then
         read(kid(2:4),'(I3)',ERR=1,END=9000) id
         idx=id/100
         id=mod(id,100)
         call pages
            do i=1,16
               call eqgouta(idx*100+16*(ID-1)+i,13+i)
            enddo 
         call pagee
         call pages
            do i=17,24
               call eqgouta(idx*100+16*(ID-1)+i,-3+i)
            enddo 
         call pagee
      elseif(kid1.eq.'P') then
         read(kid(2:2),'(I1)',ERR=1,END=9000) id
         call pages
         call eqgr2d(0,rg,zg,psi,nr,nr,nz,'@psi(R,Z)@',id)
         call pagee
         call pages
         call eqgr2d(0,rg,zg,rcu,nr,nr,nz,'@rcu(R,Z)@',id)
         call pagee
      elseif(kid1.eq.'X') then
         goto 9000
      endif
      goto 1
c
 9000 return
      end subroutine eqgout
c
c=======================================================================
      subroutine eqgouta(id,ipos)
c=======================================================================
c     interface eqiulibrium <>transport                            JAERI
c          transport grid -> equilibrium grid
c=======================================================================
      use aaa_mod
      use eqv_mod
      use eqt_mod
      use trn_mod
      implicit none
      integer n
      real*8 rovv(ivdm)
c
      integer id,ipos
c
      rovv(1)=0.d0
      do n=2,nv
         rovv(n)=sqrt(hiv(n)/hiv(nv))
      enddo
c
      if(id.eq. 1) call eqgr1d(ipos,hiv,vlv,nv,'@volume(psi_t)@')
      if(id.eq. 2) call eqgr1d(ipos,hiv,arv,nv,'@area(psi_t)@')
      if(id.eq. 3) call eqgr1d(ipos,hiv,hiv,nv,'@psit(psi_t)@')
      if(id.eq. 4) call eqgr1d(ipos,hiv,siv,nv,'@psip(psi_t)@')
      if(id.eq. 5) call eqgr1d(ipos,hiv,sdv,nv,'@dpsip/dV(psi_t)@')
      if(id.eq. 6) call eqgr1d(ipos,hiv,hdv,nv,'@dpsit/dV(psi_t)@')
      if(id.eq. 7) call eqgr1d(ipos,hiv,ckv,nv,
     &                                   '@ave(gradV^2/R^2)(psi_t)@')
      if(id.eq. 8) call eqgr1d(ipos,hiv,ssv,nv,
     &                                   '@ave(grad V^2)(psi_t)@')
      if(id.eq. 9) call eqgr1d(ipos,hiv,aav,nv,'@ave(1/R^2)(psi_t)@')
      if(id.eq.10) call eqgr1d(ipos,hiv,rrv,nv,'@<R^2>(psi_t)@')
      if(id.eq.11) call eqgr1d(ipos,hiv,bbv,nv,'@<B^2>(psi_t)@')
      if(id.eq.12) call eqgr1d(ipos,hiv,biv,nv,'@<1/B^2>(psi_t)@')
      if(id.eq.13) call eqgr1d(ipos,hiv,r2b2,nv,
     &                                   '@<grad_psi_p^2>(psi_t)@')
      if(id.eq.14) call eqgr1d(ipos,hiv,rpv,nv,
     &                                   '@minor_radius(psi_t)@')
      if(id.eq.15) call eqgr1d(ipos,hiv,rtv,nv,
     &                                   '@major_radius(psi_t)@')
      if(id.eq.16) call eqgr1d(ipos,hiv,elv,nv,'@ellipticity(psi_t)@')
      if(id.eq.17) call eqgr1d(ipos,hiv,dlv,nv,
     &                                   '@triangularity(psi_t)@')
      if(id.eq.18) call eqgr1d(ipos,hiv,muv,nv,'@muv(psi_t)@')
      if(id.eq.19) call eqgr1d(ipos,hiv,nuv,nv,'@nuv(psi_t)@')
      if(id.eq.20) call eqgr1d(ipos,hiv,prv,nv,'@prv(psi_t)@')
      if(id.eq.21) call eqgr1d(ipos,hiv,qqv,nv,'@qqv(psi_t)@')
      if(id.eq.22) call eqgr1d(ipos,hiv,sha,nv,'@sha(psi_t)@')
      if(id.eq.23) call eqgr1d(ipos,hiv,cuv,nv,'@cuv(psi_t)@')
      if(id.eq.24) call eqgr1d(ipos,hiv,pds,nv,'@pds(psi_t)@')
      if(id.eq.25) call eqgr1d(ipos,hiv,qdv,nv,'@qdv(psi_t)@')
      if(id.eq.26) call eqgr1d(ipos,hiv,rbv,nv,'@rbv(psi_t)@')
      if(id.eq.27) call eqgr1d(ipos,hiv,epv,nv,'@epv(psi_t)@')
      if(id.eq.28) call eqgr1d(ipos,hiv,ftr,nv,'@ftr(psi_t)@')
c
      if(id.eq.101) call eqgr1d(ipos,rovv,vlv,nv,'@volume(eq_grid)@')
      if(id.eq.102) call eqgr1d(ipos,rovv,arv,nv,'@area(eq_grid)@')
      if(id.eq.103) call eqgr1d(ipos,rovv,hiv,nv,'@psit(eq_grid)@')
      if(id.eq.104) call eqgr1d(ipos,rovv,siv,nv,'@psip(eq_grid)@')
      if(id.eq.105) call eqgr1d(ipos,rovv,sdv,nv,'@dpsip/dV(eq_grid)@')
      if(id.eq.106) call eqgr1d(ipos,rovv,hdv,nv,'@dpsit/dV(eq_grid)@')
      if(id.eq.107) call eqgr1d(ipos,rovv,ckv,nv,
     &                                   '@ave(gradV^2/R^2)(eq_grid)@')
      if(id.eq.108) call eqgr1d(ipos,rovv,ssv,nv,
     &                                   '@ave(grad V^2)(eq_grid)@')
      if(id.eq.109) call eqgr1d(ipos,rovv,aav,nv,
     &                                   '@ave(1/R^2)(eq_grid)@')
      if(id.eq.110) call eqgr1d(ipos,rovv,rrv,nv,'@<R^2>(eq_grid)@')
      if(id.eq.111) call eqgr1d(ipos,rovv,bbv,nv,'@<B^2>(eq_grid)@')
      if(id.eq.112) call eqgr1d(ipos,rovv,biv,nv,'@<1/B^2>(eq_grid)@')
      if(id.eq.113) call eqgr1d(ipos,rovv,r2b2,nv,
     &                                   '@<grad_psi_p^2>(eq_grid)@')
      if(id.eq.114) call eqgr1d(ipos,rovv,rpv,nv,
     &                                   '@minor_radius(eq_grid)@')
      if(id.eq.115) call eqgr1d(ipos,rovv,rtv,nv,
     &                                   '@major_radius(eq_grid)@')
      if(id.eq.116) call eqgr1d(ipos,rovv,elv,nv,
     &                                   '@ellipticity(eq_grid)@')
      if(id.eq.117) call eqgr1d(ipos,rovv,dlv,nv,
     &                                   '@triangularity(eq_grid)@')
c      if(id.eq.118) call eqgr1d(ipos,roh,fth,nv,'@trapped_frac(eq_grid)@')
c      if(id.eq.119) call eqgr1d(ipos,roh,eph,nv,'@trapped_frac^2(eq_grid)@')
c
      if(id.eq.201) call eqgr1d(ipos, ro,vlt,nt,'@volume(tr_grid)@')
      if(id.eq.202) call eqgr1d(ipos, ro,art,nt,'@area(tr_grid)@')
      if(id.eq.203) call eqgr1d(ipos, ro,hit,nt,'@psi_t(tr_grid)@')
      if(id.eq.204) call eqgr1d(ipos, ro,sit,nt,'@psi_p(tr_grid)@')
      if(id.eq.205) call eqgr1d(ipos, ro,sdt,nt,'@dpsi_p/dV(tr_grid)@')
      if(id.eq.206) call eqgr1d(ipos, ro,hdt,nt,'@dpsi_t/dV(tr_grid)@')
      if(id.eq.207) call eqgr1d(ipos, ro,ckt,nt,
     &                                    '@<grad_V^2/R^2>(tr_grid)@')
      if(id.eq.208) call eqgr1d(ipos, ro,sst,nt,'@<grad_V^2>(tr_grid)@')
      if(id.eq.209) call eqgr1d(ipos, ro,aat,nt,'@<1/R^2>(tr_grid)@')
      if(id.eq.210) call eqgr1d(ipos, ro,rrt,nt,'@<R^2>(tr_grid)@')
c      if(id.eq.211) call eqgr1d(ipos, ro,bbt,nt,'@<B^2>(tr_grid)@')
c      if(id.eq.212) call eqgr1d(ipos, ro,bit,nt,'@<1/B^2>(tr_grid)@')
c      if(id.eq.213) call eqgr1d(ipos, ro,r2b2t,nt,'@<grad_psi_p^2>(tr_grid)@')
c      if(id.eq.214) call eqgr1d(ipos, ro,rpt,nt,'@minor_radius(tr_grid)@')
c      if(id.eq.215) call eqgr1d(ipos, ro,rtt,nt,'@major_radius(tr_grid)@')
c      if(id.eq.216) call eqgr1d(ipos, ro,elt,nt,'@ellipticity(tr_grid)@')
c      if(id.eq.217) call eqgr1d(ipos, ro,dlt,nt,'@triangularity(tr_grid)@')
c      if(id.eq.218) call eqgr1d(ipos, ro,ftt,nt,'@trapped_frac(tr_grid)@')
c      if(id.eq.219) call eqgr1d(ipos, ro,ept,nt,'@trapped_frac^2(tr_grid)@')
c
      if(id.eq.301) call eqgr1d(ipos,roh,vlh,ntm,'@volume(tr_h_grid)@')
c      if(id.eq.302) call eqgr1d(ipos,roh,arh,ntm,'@area(tr_h_grid)@')
      if(id.eq.303) call eqgr1d(ipos,roh,hih,ntm,'@psi_t(tr_h_grid)@')
c      if(id.eq.304) call eqgr1d(ipos,roh,sih,ntm,'@psi_p(tr_h_grid)@')
c      if(id.eq.305) call eqgr1d(ipos,roh,sdh,ntm,'@dpsi_p/dV(tr_h_grid)@')
      if(id.eq.306) call eqgr1d(ipos,roh,hdh,ntm,
     &                                   '@dpsi_t/dV(tr_h_grid)@')
      if(id.eq.307) call eqgr1d(ipos,roh,ckh,ntm,
     &                                   '@<grad_V^2/R^2>(tr_h_grid)@')
      if(id.eq.308) call eqgr1d(ipos,roh,ssh,ntm,
     &                                   '@<grad_V^2>(tr_h_grid)@')
      if(id.eq.309) call eqgr1d(ipos,roh,aah,ntm,'@<1/R^2>(tr_h_grid)@')
c      if(id.eq.310) call eqgr1d(ipos,roh,rrh,ntm,'@<R^2>(tr_h_grid)@')
      if(id.eq.311) call eqgr1d(ipos,roh,bbh,ntm,'@<B^2>(tr_h_grid)@')
      if(id.eq.312) call eqgr1d(ipos,roh,bih,ntm,'@<1/B^2>(tr_h_grid)@')
      if(id.eq.313) call eqgr1d(ipos,roh,r2b2h,ntm,
     &                                    '@<grad_psi_p^2>(tr_h_grid)@')
      if(id.eq.314) call eqgr1d(ipos,roh,rph,ntm,
     &                                      '@minor_radius(tr_h_grid)@')
      if(id.eq.315) call eqgr1d(ipos,roh,rth,ntm,
     &                                      '@major_radius(tr_h_grid)@')
      if(id.eq.316) call eqgr1d(ipos,roh,elh,ntm,
     &                                       '@ellipticity(tr_h_grid)@')
      if(id.eq.317) call eqgr1d(ipos,roh,dlh,ntm,
     &                                     '@triangularity(tr_h_grid)@')
      if(id.eq.318) call eqgr1d(ipos,roh,fth,ntm,
     &                                      '@trapped_frac(tr_h_grid)@')
      if(id.eq.319) call eqgr1d(ipos,roh,eph,ntm,
     &                                    '@trapped_frac^2(tr_h_grid)@')
      return
      end subroutine eqgouta
c
c=======================================================================
      subroutine eqgr1d(ngp,x,y,nxmax,str)
c
      USE libgrf,ONLY: grf1d
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
      end subroutine eqgr1d
c
c=======================================================================
      subroutine eqgr2d(ngp,x,y,z,nxm,nxmax,nymax,str,ntype)
!     &                  ngline,
!     &                  ngcoil,rgcoil,zgcoil,
!     &                  ngbound,rgboud,zgboud)
c
      USE aaa_mod,ONLY: icvdm
      USE vac_mod,ONLY: rcoil,zcoil,ccoil,ncoil
      USE libgrf,ONLY: grfut1,grfut2,grfut3,grfut4, 
     &                 grf2dax,grf2dbx,grf2dcx
      implicit none
      integer ngp,nxm,nxmax,nymax,nx,ny,ntype,i,n
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
      real*4 gcross,gxc,gyc
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
      gsymin=gymin
      gsymax=gymax
      gsxmin=gxmin
      gsxmax=gxmax
!      IF(ngline.EQ.0) THEN
         gzstep=0.2*gzstep
!      ELSE
!         gzstep=(gzmax-gzmin)/ngline
!      ENDDO
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
     &                gx,gy,gz,nxm,nxmax,nymax,str,0,0)
c
      elseif(ntype.eq.2) then
         call grf2dax(gp,
     &                gsxmin,gsxmax,gxstep,gxorg,
     &                gsymin,gsymax,gystep,gyorg,
     &                gszmin,gszmax,gzstep,gzorg,
     &                gx,gy,gz,nxm,nxmax,nymax,str,0,0)
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
     &                gx,gy,gz,nxm,nxmax,nymax,
     &                str,gp3d)
      endif

C      gcross=MIN(gxmax-gxmin,gymax-gymin)/20.0
C      CALL SETRGB(0.0,1.0,0.0)
C      DO n=1,icvdm
C         IF(ncoil(n).GT.0) THEN
C            DO i=1,ncoil(n)
C               IF(ABS(ccoil(i,n)).GT.1.E-8) THEN
C                  WRITE(6,'(A,2I5,1P3E12.4)')
C     &                 'COIL:',n,i,ccoil(i,n),rcoil(i,n),zcoil(i,n)
C                  gxc=GUCLIP(rcoil(i,n))
C                  gyc=GUCLIP(zcoil(i,n))
C                  CALL MOVE2D(gxc,gyc-gcross)
C                  CALL DRAW2D(gxc,gyc+gcross)
C                  CALL MOVE2D(gxc-gcross,gyc)
C                  CALL DRAW2D(gxc+gcross,gyc)
C               END IF
C            END DO
C         END IF
C      END DO
            
c
      deallocate(ka)
      deallocate(gz)
      deallocate(gy)
      deallocate(gx)
      return
      end subroutine eqgr2d
c
      end module eqgout_mod
