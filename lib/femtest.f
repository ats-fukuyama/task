!     $Id$

      program femtest

      implicit none
      integer:: mwmax,mlmax
      complex(8),dimension(:,:),allocatable:: fma
      complex(8),dimension(:),allocatable:: fvb,fvx
      real(8),dimension(:),allocatable:: rho
      complex(8),dimension(:),allocatable:: cf1,cf2,cf3

      integer:: id,nrmax,npow,nth,nph,nr
      real(8):: rf,angl

      call gsopen

      write(6,*) ' * id=-2: Laplace eq (cylinder)'
      write(6,*) ' * id=-1: Laplace eq (slab)'
      write(6,*) ' * id= 0: stop'
      write(6,*) ' * id= 1: Maxwell eq (slab) Linear'
      write(6,*) ' * id= 2: Maxwell eq (slab) Linear + edge 0th'
      write(6,*) ' * id= 3: Maxwell eq (slab) Linear + edge 1st'
      write(6,*) ' * id= 4: Maxwell eq (slab) Hermite'
      write(6,*) ' * id= 6: Maxwell eq (cylinder) Linear'
      write(6,*) ' * id= 7: Maxwell eq (cylinder) Linear + edge 0th'
      write(6,*) ' * id= 8: Maxwell eq (cylinder) Linear + edge 1st'
      write(6,*) ' * id= 9: Maxwell eq (cylinder) Hermite'
      write(6,*) ' * id=10: Maxwell eq (cylinder) Hermite+ quadra'
      write(6,*) ' * id=11: Maxwell eq (cylinder) quadra + linear'
      write(6,*) ' * *******************************'

      id=10
      nrmax=11
      npow=1
      nth=0
      nph=0
      rf=1.d0
      angl=0.d0

    1 write(6,'(A,I3,I5,I3,I3,I3,1PE12.4,0PF8.4)') 
     &      '## input: id,nrmax,npow,nth,nph,rf,angl=',
     &                 id,nrmax,npow,nth,nph,rf,angl
      read(5,*,err=1,end=9000) id,nrmax,npow,nth,nph,rf,angl
      if(id.eq.0) goto 9000

c$$$      call fem_exec(id,nrmax,npow,nth,nph,rf,angl)
c$$$      call pages
c$$$      call femgr1dc( 1,rho,cf1,nrmax,'@cf1(rho)@')
c$$$      call femgr1dc( 2,rho,cf2,nrmax,'@cf2(rho)@')
c$$$      call femgr1dc( 3,rho,cf3,nrmax,'@cf3(rho)@')
c$$$      call pagee
c$$$      do nr=1,nrmax
c$$$         write(6,'(I5,1P7E10.2)') nr,rho(nr),cf1(nr),cf2(nr),cf3(nr)
c$$$      enddo

      call pages
      call fem_exec(id,nrmax,npow,nth,nph,rf,angl)
      call femgr1dc( 5,rho,cf1,nrmax,'@cf1(rho)@')
      call femgr1dc( 6,rho,cf2,nrmax,'@cf2(rho)@')
      call femgr1dc( 7,rho,cf3,nrmax,'@cf3(rho)@')
      call fem_exec(id,nrmax,npow,nth+1,nph,rf,angl)
      call femgr1dc( 8,rho,cf1,nrmax,'@cf1(rho)@')
      call femgr1dc( 9,rho,cf2,nrmax,'@cf2(rho)@')
      call femgr1dc(10,rho,cf3,nrmax,'@cf3(rho)@')
      call fem_exec(id,nrmax,npow,nth+2,nph,rf,angl)
      call femgr1dc(11,rho,cf1,nrmax,'@cf1(rho)@')
      call femgr1dc(12,rho,cf2,nrmax,'@cf2(rho)@')
      call femgr1dc(13,rho,cf3,nrmax,'@cf3(rho)@')
      call pagee

      goto 1

 9000 call gsclos
      stop

      contains

      subroutine fem_exec(id,nrmax,npow,nth,nph,rf,angl)

      implicit none
      integer,intent(in):: id,nrmax,npow,nth,nph
      real(8),intent(in):: rf,angl
      integer:: mw,ml,nr,ierr

      select case(id)
         case(-1)
            call fem_calc_l1(nrmax,npow)
         case(-2)
            call fem_calc_l2(nrmax,npow)
         case(1)
            call fem_calc_x1(nrmax,npow,nth,nph,rf,angl)
         case(2)
            call fem_calc_x2(nrmax,npow,nth,nph,rf,angl)
         case(3)
            call fem_calc_x3(nrmax,npow,nth,nph,rf,angl)
         case(4)
            call fem_calc_x4(nrmax,npow,nth,nph,rf,angl)
         case(6)
            call fem_calc_r1(nrmax,npow,nth,nph,rf,angl)
         case(7)
            call fem_calc_r2(nrmax,npow,nth,nph,rf,angl)
         case(8)
            call fem_calc_r3(nrmax,npow,nth,nph,rf,angl)
         case(9)
            call fem_calc_r4(nrmax,npow,nth,nph,rf,angl)
         case(10)
            call fem_calc_r5(nrmax,npow,nth,nph,rf,angl)
         case(11)
            call fem_calc_r6(nrmax,npow,nth,nph,rf,angl)
      end select

c$$$      do ml=1,mlmax
c$$$         write(16,'(1P6E12.4)') (fma(mw,ml),mw=1,mwmax)
c$$$      enddo
c$$$      write(16,'(1P6E12.4)') (fvb(ml),ml=1,mlmax)

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
      enddo
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      if(ierr.ne.0) write(6,*) '# ierr= ',ierr

      if(id.eq.-1.or.id.eq.-2) then
         do nr=1,nrmax
            cf1(nr)=fvx(2*(nr-1)+1)
            cf2(nr)=fvx(2*(nr-1)+2)
            cf3(nr)=0.d0
         enddo

      else if(id.eq.1.or.id.eq.6) then
         do nr=1,nrmax
            cf1(nr)=fvx(3*(nr-1)+1)
            cf2(nr)=fvx(3*(nr-1)+2)
            cf3(nr)=fvx(3*(nr-1)+3)
         enddo

      else if(id.eq.2.or.id.eq.7) then
         do nr=1,nrmax
            cf1(nr)=fvx(3*(nr-1)+1)
            cf2(nr)=fvx(3*(nr-1)+2)
            cf3(nr)=fvx(3*(nr-1)+3)
         enddo

      else if(id.eq.3.or.id.eq.8) then
         do nr=1,nrmax
            cf1(nr)=fvx(4*(nr-1)+1)
            cf2(nr)=fvx(4*(nr-1)+3)
            cf3(nr)=fvx(4*(nr-1)+4)
         enddo

      else if(id.eq.4.or.id.eq.9) then
         do nr=1,nrmax
            cf1(nr)=fvx(6*(nr-1)+1)
            cf2(nr)=fvx(6*(nr-1)+3)
            cf3(nr)=fvx(6*(nr-1)+5)
         enddo

      else if(id.eq.10) then
         do nr=1,nrmax
            cf1(nr)=fvx(6*(nr-1)+1)
            cf2(nr)=fvx(6*(nr-1)+3)
            cf3(nr)=fvx(6*(nr-1)+5)
         enddo
         cf1(nrmax)=cf1(nrmax-1)

      else if(id.eq.11) then
         do nr=1,nrmax
            cf1(nr)=fvx(5*(nr-1)+1)
            cf2(nr)=fvx(5*(nr-1)+2)
            cf3(nr)=fvx(5*(nr-1)+4)
         enddo
         cf1(nrmax)=cf1(nrmax-1)
      endif
      return
      end subroutine fem_exec

!----- set profile -----

      subroutine mesh_init(nrmax,npow)

      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh position
      integer,save:: nrmax_save=0
      real(8):: drho
      integer:: nr

      if(nrmax.ne.nrmax_save) then
         if(allocated(rho)) deallocate(rho)
         if(allocated(cf1)) deallocate(cf1)
         if(allocated(cf2)) deallocate(cf2)
         if(allocated(cf3)) deallocate(cf3)
         allocate(rho(nrmax))
         allocate(cf1(nrmax))
         allocate(cf2(nrmax))
         allocate(cf3(nrmax))
      endif
      nrmax_save=nrmax

      drho=1.d0/(nrmax-1)**npow
      do nr=1,nrmax
         rho(nr)=drho*(nr-1)**npow
      enddo
      return
      end subroutine mesh_init

!----- allocate matrix and vectors -----

      subroutine fem_init(nrmax,nvmax)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: nvmax  ! numbre of variables
      integer,save:: mwmax_save=0,mlmax_save=0
      integer:: mc

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif
    
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(allocated(fma)) deallocate(fma)
         allocate(fma(mwmax,mlmax))
         mc=(mwmax+1)/2
         write(6,'(A,3I5)') 'mlmax,mwmax,mc=',mlmax,mwmax,mc
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

      subroutine fem_calc_l1(nrmax,npow)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer:: nr,ml,mw,mc,nvmax,i,j,k
      real(8):: drho,val1,val2,val3,val4,temp
      real(8),dimension(4):: coef

      call mesh_init(nrmax,npow)

      nvmax=2

      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

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

      return
      end subroutine fem_calc_l1

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_l2(nrmax,npow)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer:: nr,ml,mw,mc,nvmax,i,j,k
      real(8):: drho,val1,val2,val3,val4,temp
      real(8),dimension(4):: coef

      call mesh_init(nrmax,npow)

      nvmax=2

      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

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

      return
      end subroutine fem_calc_l2

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_x1(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=3
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      rkth=nth
      rkph=nph
      factor=rf**2
      
      do inod=1,2
         do k=1,4
            do j=1,3
               do i=1,3
                  fmc(i,j,k,inod)=0.d0
               enddo
            enddo
         enddo
      enddo

      do inod=1,2
         fmc(1,1,1,inod)= factor-rkth**2-rkph**2
         fmc(2,2,1,inod)= factor-rkph**2
         fmc(2,3,1,inod)= rkth*rkph
         fmc(3,2,1,inod)= rkth*rkph
         fmc(3,3,1,inod)= factor-rkth**2
      
         fmc(2,1,2,inod)= ci*rkth
         fmc(3,1,2,inod)= ci*rkph

         fmc(1,2,3,inod)=-ci*rkth
         fmc(1,3,3,inod)=-ci*rkph

         fmc(2,2,4,inod)=-1.d0
         fmc(3,3,4,inod)=-1.d0
      enddo

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)
         ml=3*(nr-1)+1

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
            do i=1,3
               ml=3*(nr-1)+(i-1)+1
               mw=mc+(j-1)-(i-1)
               do inod=1,2

! ----- non derivative terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho

               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho

! ----- dw/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1)
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1)

               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2)
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2)

! ----- dE/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3)
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3)

               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4)
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4)

! ----- dw/dx dE/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho

               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho
               enddo
            enddo
         enddo
      enddo

      do mw=1,mwmax
!         fma(mw,1) = 0.d0
         fma(mw,2) = 0.d0
         fma(mw,3) = 0.d0
      enddo
!      fma(mc,1)=1.d0
      fma(mc,2)=1.d0
      fma(mc,3)=1.d0

      do mw=1,mwmax
!         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
!      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(3*(nr-1)+2)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(3* nr   +2)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(3*(nr-1)+3)=(rho(nr+1)-0.85d0)*angl
            fvb(3* nr   +3)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_x1

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_x2(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=3
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      rkth=nth
      rkph=nph
      factor=rf**2
      
      do inod=1,2
         do k=1,4
            do j=1,3
               do i=1,3
                  fmc(i,j,k,inod)=0.d0
               enddo
            enddo
         enddo
      enddo

      do inod=1,2
         fmc(1,1,1,inod)= factor-rkth**2-rkph**2
         fmc(2,2,1,inod)= factor-rkph**2
         fmc(2,3,1,inod)= rkth*rkph
         fmc(3,2,1,inod)= rkth*rkph
         fmc(3,3,1,inod)= factor-rkth**2
      
         fmc(2,1,2,inod)= ci*rkth
         fmc(3,1,2,inod)= ci*rkph

         fmc(1,2,3,inod)=-ci*rkth
         fmc(1,3,3,inod)=-ci*rkph

         fmc(2,2,4,inod)=-1.d0
         fmc(3,3,4,inod)=-1.d0
      enddo

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
            do i=1,3
               ml=3*(nr-1)+(i-1)+1
               mw=mc+(j-1)-(i-1)
               
! ***** (1,1) *****
            if((i.eq.1).and.(j.eq.1)) then

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,1)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,2)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,1)
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,2)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,3)
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,4)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,1)
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,1)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,2)
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,2)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               enddo

            else

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1)
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1)
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2)
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3)
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3)
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4)
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho
               enddo
            endif

            enddo
         enddo
      enddo

      do mw=1,mwmax
         fma(mw,2) = 0.d0
         fma(mw,3) = 0.d0
      enddo
      fma(mc,2)=1.d0
      fma(mc,3)=1.d0

      do mw=1,mwmax
         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(3*(nr-1)+2)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(3* nr   +2)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(3*(nr-1)+3)=(rho(nr+1)-0.85d0)*angl
            fvb(3* nr   +3)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_x2

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_x3(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer,dimension(3),parameter :: ishift=(/0,2,3/)

      call mesh_init(nrmax,npow)

      nvmax=4
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      rkth=nth
      rkph=nph
      factor=rf**2
      
      do inod=1,2
         do k=1,4
            do j=1,3
               do i=1,3
                  fmc(i,j,k,inod)=0.d0
               enddo
            enddo
         enddo
      enddo

      do inod=1,2
         fmc(1,1,1,inod)= factor-rkth**2-rkph**2
         fmc(2,2,1,inod)= factor-rkph**2
         fmc(2,3,1,inod)= rkth*rkph
         fmc(3,2,1,inod)= rkth*rkph
         fmc(3,3,1,inod)= factor-rkth**2
      
         fmc(2,1,2,inod)= ci*rkth
         fmc(3,1,2,inod)= ci*rkph

         fmc(1,2,3,inod)=-ci*rkth
         fmc(1,3,3,inod)=-ci*rkph

         fmc(2,2,4,inod)=-1.d0
         fmc(3,3,4,inod)=-1.d0
      enddo

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
            do i=1,3
               ml=4*(nr-1)+ishift(i)+1
               mw=mc+ishift(j)-ishift(i)
               
! ***** (1,1) *****
            if((i.eq.1).and.(j.eq.1)) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,1)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,2)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               fma(mw+1,ml  )=fma(mw+1 ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,3)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,3)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,4)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgg(inod,3,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,3,2)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,4,2)/drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgg(inod,3,3)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,4,3)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,3,4)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,4,4)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,1)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgl(inod,3,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,3,3)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,2)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               fma(mw+3,ml+1)=fma(mw+3,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgl(inod,3,2)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,4,2)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,3,4)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,4,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,2)
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw+1 ,ml  )=fma(mw+1,ml  )
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,3)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,3)
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,4)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4)
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,2)
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               fma(mw-3,ml+4)=fma(mw-3,ml+4)
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,3)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,3)
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,4)
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,4)/drho
               enddo

            else

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+4)=fma(mw  ,ml+4)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho

               enddo
            endif

            enddo
         enddo
      enddo

      do mw=1,mwmax
         fma(mw,3) = 0.d0
         fma(mw,4) = 0.d0
      enddo
      fma(mc,3)=1.d0
      fma(mc,4)=1.d0

      do mw=1,mwmax
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(4*(nr-1)+3)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(4* nr   +3)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(4*(nr-1)+4)=(rho(nr+1)-0.85d0)*angl
            fvb(4* nr   +4)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_x3

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_x4(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor
      complex(8),dimension(3,3,4,4):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=6
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      rkth=nth
      rkph=nph
      
      do inod=1,4
         do k=1,4
            do j=1,3
               do i=1,3
                  fmc(i,j,k,inod)=0.d0
               enddo
            enddo
         enddo
      enddo

      factor=rf**2
      do inod=1,4
         fmc(1,1,1,inod)= factor-rkth**2-rkph**2
         fmc(2,2,1,inod)= factor-rkph**2
         fmc(2,3,1,inod)= rkth*rkph
         fmc(3,2,1,inod)= rkth*rkph
         fmc(3,3,1,inod)= factor-rkth**2
      
         fmc(2,1,2,inod)= ci*rkth
         fmc(3,1,2,inod)= ci*rkph

         fmc(1,2,3,inod)=-ci*rkth
         fmc(1,3,3,inod)=-ci*rkph

         fmc(2,2,4,inod)=-1.d0
         fmc(3,3,4,inod)=-1.d0
      enddo

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)
         ml=6*(nr-1)+1

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=0.5d0*(-11*fmc(i,j,k,1)+18*fmc(i,j,k,2)
     &                                - 9*fmc(i,j,k,3)+ 2*fmc(i,j,k,4))
                  fmd(i,j,k,3)=fmc(i,j,k,4)
                  fmd(i,j,k,4)=0.5d0*(- 2*fmc(i,j,k,1)+ 9*fmc(i,j,k,2)
     &                                -18*fmc(i,j,k,3)+11*fmc(i,j,k,4))
               enddo
            enddo
         enddo

         do j=1,3
         do i=1,3
            ml=6*(nr-1)+2*(i-1)+1
            mw=mc+2*(j-1)-2*(i-1)
            do inod=1,4

! ----- non derivative terms

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,1)*drho
            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,1)*drho**2
            fma(mw-6,ml+6)=fma(mw-6,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,1)*drho
            fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,1)*drho**2

            fma(mw+1,ml  )=fma(mw+1,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,2)*drho**2
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,2)*drho**3
            fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,2)*drho**2
            fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,2)*drho**3

            fma(mw+6,ml  )=fma(mw+6,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,3)*drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,3)*drho**2
            fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,3)*drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,3)*drho**2

            fma(mw+7,ml  )=fma(mw+7,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,4)*drho**2
            fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,4)*drho**3
            fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,4)*drho**2
            fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,4)*drho**3

! ----- dw/dx terms

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,1)
            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,1)*drho
            fma(mw-6,ml+6)=fma(mw-6,ml+6)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,1)
            fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,1)*drho

            fma(mw+1,ml  )=fma(mw+1,ml  )
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,2)*drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,2)*drho**2
            fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,2)*drho
            fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,2)*drho**2

            fma(mw+6,ml  )=fma(mw+6,ml  )
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,3)
            fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,3)*drho
            fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,3)
            fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,3)*drho

            fma(mw+7,ml  )=fma(mw+7,ml )
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,4)*drho
            fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,4)*drho**2
            fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,4)*drho
            fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,4)*drho**2

! ----- dE/dx terms

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,5)
            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,5)*drho
            fma(mw-6,ml+6)=fma(mw-6,ml+6)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,5)
            fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,5)*drho

            fma(mw+1,ml  )=fma(mw+1,ml  )
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,6)*drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,6)*drho**2
            fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,6)*drho
            fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,6)*drho**2

            fma(mw+6,ml  )=fma(mw+6,ml  )
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,7)
            fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,7)*drho
            fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,7)
            fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,7)*drho

            fma(mw+7,ml  )=fma(mw+7,ml  )
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,8)*drho
            fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,8)*drho**2
            fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,8)*drho
            fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,8)*drho**2

! ----- dw/dx dE/dx terms

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,5)/drho
            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,5)
            fma(mw-6,ml+6)=fma(mw-6,ml+6)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,5)/drho
            fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,5)

            fma(mw+1,ml  )=fma(mw+1,ml  )
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,6)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,6)*drho
            fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,6)
            fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,6)*drho

            fma(mw+6,ml  )=fma(mw+6,ml  )
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,7)/drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,7)
            fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,7)/drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,7)

            fma(mw+7,ml  )=fma(mw+7,ml  )
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,8)
            fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,8)*drho
            fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,8)
            fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,8)*drho
               enddo
            enddo
         enddo
      enddo

      do mw=1,mwmax
!         fma(mw,1) = 0.d0
         fma(mw,3) = 0.d0
         fma(mw,5) = 0.d0
      enddo
!      fma(mc,1)=1.d0
      fma(mc,3)=1.d0
      fma(mc,5)=1.d0

      do mw=1,mwmax
!         fma(mw,mlmax-5) = 0.d0
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
!      fma(mc,mlmax-5) = 1.d0
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(6*(nr-1)+3)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(6* nr   +3)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(6*(nr-1)+5)=(rho(nr+1)-0.85d0)*angl
            fvb(6* nr   +5)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_x4

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r1(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=3
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         factor=rf**2
         rkph=nph
      
         do inod=1,2
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,2
            if(inod.eq.1) then
               if(nr.eq.1) then
                  rho0=rho(nr+1)/3.d0
               else
                  rho0=rho(nr)
               endif
            else
               rho0=rho(nr+1)
            endif
            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkth**2-rkph**2)
            fmc(1,2,1,inod)= -ci*rkth
            fmc(2,1,1,inod)= +ci*rkth
            fmc(2,2,1,inod)= rho0*(factor-rkph**2)-rkth0
            fmc(2,3,1,inod)= nth*rkph
            fmc(3,2,1,inod)= nth*rkph
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)=       ci*nth
            fmc(2,2,2,inod)= -1.d0
            fmc(3,1,2,inod)= +rho0*ci*rkph

            fmc(1,2,3,inod)= -     ci*nth
            fmc(2,2,3,inod)= -1.d0
            fmc(1,3,3,inod)= -rho0*ci*rkph

            fmc(2,2,4,inod)=-rho0
            fmc(3,3,4,inod)=-rho0
         enddo

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
            do i=1,3
               ml=3*(nr-1)+(i-1)+1
               mw=mc+(j-1)-(i-1)
               do inod=1,2

! ----- non derivative terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho

               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,2) = 0.d0
         enddo
         fma(mc,2)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,2) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc-1,2)=1.d0
         fma(mc,2)=ci*nth
         fma(mc,3)=1.d0
      else
         do mw=1,mwmax
            fma(mw,2) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc,2)=1.d0
         fma(mc,3)=1.d0
      endif


      do mw=1,mwmax
!         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
!      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(3*(nr-1)+2)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(3* nr   +2)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(3*(nr-1)+3)=(rho(nr+1)-0.85d0)*angl
            fvb(3* nr   +3)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_r1

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r2(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=3
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         factor=rf**2
         rkph=nph
      
         do inod=1,2
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,2
            if(inod.eq.1) then
               if(nr.eq.1) then
                  rho0=rho(nr+1)/3.d0
               else
                  rho0=rho(nr)
               endif
            else
               rho0=rho(nr+1)
            endif
            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkth**2-rkph**2)
            fmc(1,2,1,inod)= -ci*rkth
            fmc(2,1,1,inod)= +ci*rkth
            fmc(2,2,1,inod)= rho0*(factor-rkph**2)-rkth0
            fmc(2,3,1,inod)= nth*rkph
            fmc(3,2,1,inod)= nth*rkph
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)=       ci*nth
            fmc(2,2,2,inod)= -1.d0
            fmc(3,1,2,inod)= +rho0*ci*rkph

            fmc(1,2,3,inod)= -     ci*nth
            fmc(2,2,3,inod)= -1.d0
            fmc(1,3,3,inod)= -rho0*ci*rkph

            fmc(2,2,4,inod)=-rho0
            fmc(3,3,4,inod)=-rho0
         enddo

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
            do i=1,3
               ml=3*(nr-1)+(i-1)+1
               mw=mc+(j-1)-(i-1)
               
! ***** (1,1) *****
            if((i.eq.1).and.(j.eq.1)) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,1)
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,2)
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,1)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,2)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,2)
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,2)
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               enddo

            else

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho
               enddo
            endif

            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,2) = 0.d0
         enddo
         fma(mc,2)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,2) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc-1,2)=1.d0
         fma(mc,2)=ci*nth
         fma(mc,3)=1.d0
      else
         do mw=1,mwmax
            fma(mw,2) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc,2)=1.d0
         fma(mc,3)=1.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(3*(nr-1)+2)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(3* nr   +2)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(3*(nr-1)+3)=(rho(nr+1)-0.85d0)*angl
            fvb(3* nr   +3)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_r2

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r3(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer,dimension(3),parameter :: ishift=(/0,2,3/)

      call mesh_init(nrmax,npow)

      nvmax=4
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         factor=rf**2
         rkph=nph
      
         do inod=1,2
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,2
            if(inod.eq.1) then
               if(nr.eq.1) then
!                  rho0=rho(nr+1)/3.d0
                  rho0=rho(nr+1)/3.d0
               else
                  rho0=rho(nr)
               endif
            else
               rho0=rho(nr+1)
            endif
            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkth**2-rkph**2)
            fmc(1,2,1,inod)= -ci*rkth
            fmc(2,1,1,inod)= +ci*rkth
            fmc(2,2,1,inod)= rho0*(factor-rkph**2)-rkth0
            fmc(2,3,1,inod)= nth*rkph
            fmc(3,2,1,inod)= nth*rkph
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)=       ci*nth
            fmc(2,2,2,inod)= -1.d0
            fmc(3,1,2,inod)= +rho0*ci*rkph

            fmc(1,2,3,inod)= -     ci*nth
            fmc(2,2,3,inod)= -1.d0
            fmc(1,3,3,inod)= -rho0*ci*rkph

            fmc(2,2,4,inod)=-rho0
            fmc(3,3,4,inod)=-rho0
         enddo

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
            do i=1,3
               ml=4*(nr-1)+ishift(i)+1
               mw=mc+ishift(j)-ishift(i)
               
! ***** (1,1) *****
            if((i.eq.1).and.(j.eq.1)) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,1)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,2)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               fma(mw+1,ml  )=fma(mw+1 ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,3)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,3)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,4)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgg(inod,3,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,3,2)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,4,2)/drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgg(inod,3,3)*drho
     &                       +fmd(i,j,2,inod)*table_lgg(inod,4,3)
     &                       +fmd(i,j,3,inod)*table_lgg(inod,3,4)
     &                       +fmd(i,j,4,inod)*table_lgg(inod,4,4)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,1)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgl(inod,3,1)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,3,3)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  )
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,2)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               fma(mw+3,ml+1)=fma(mw+3,ml+1)
     &                       +fmd(i,j,1,inod)*table_lgl(inod,3,2)*drho
     &                       +fmd(i,j,2,inod)*table_lgl(inod,4,2)
     &                       +fmd(i,j,3,inod)*table_lgl(inod,3,4)
     &                       +fmd(i,j,4,inod)*table_lgl(inod,4,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,2)
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw+1 ,ml  )=fma(mw+1,ml  )
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,3)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,3)
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,4)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4)
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,2)
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               fma(mw-3,ml+4)=fma(mw-3,ml+4)
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,3)*drho
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,3)
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,4)
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,4)/drho
               enddo

            else

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  )
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+4)=fma(mw  ,ml+4)
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2)
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4)
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho

               enddo
            endif

            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
         enddo
         fma(mc,3)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,4) = 0.d0
         enddo
         fma(mc-2,3)=1.d0
         fma(mc-1,3)=-0.5d0*rho(2)
         fma(mc,3)=ci*nth
         fma(mc,4)=1.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,2) = 0.d0
            fma(mw,3) = 0.d0
            fma(mw,4) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,2)=1.d0
         fma(mc,3)=1.d0
         fma(mc,4)=1.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0


      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(4*(nr-1)+3)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(4* nr   +3)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(4*(nr-1)+4)=(rho(nr+1)-0.85d0)*angl
            fvb(4* nr   +4)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_r3

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r4(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4,4):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=6
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         factor=rf**2
         rkph=nph
      
         do inod=1,4
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,4
            if(inod.eq.1) then
               if(nr.eq.1) then
                  rho0=(8.d0*rho(nr)+rho(nr+1))/9.d0
               else
                  rho0=rho(nr)
               endif
            elseif(inod.eq.2) then
               rho0=(2.d0*rho(nr)+rho(nr+1))/3.d0
            elseif(inod.eq.3) then
               rho0=(rho(nr)+2.d0*rho(nr+1))/3.d0
            else
               rho0=rho(nr+1)
            endif
            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkph**2-rkth**2)
            fmc(1,2,1,inod)= rho0*(-ci*rkth*rkth0)
            fmc(2,1,1,inod)= rho0*(+ci*rkth*rkth0)
            fmc(2,2,1,inod)= rho0*(factor-rkph**2-rkth0**2)
            fmc(2,3,1,inod)= rho0*(rkth*rkph)
            fmc(3,2,1,inod)= rho0*(rkth*rkph)
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)= rho0*( ci*rkth)
            fmc(2,2,2,inod)= rho0*(  -rkth0)
            fmc(3,1,2,inod)= rho0*( ci*rkph)

            fmc(1,2,3,inod)= rho0*(-ci*rkth)
            fmc(1,3,3,inod)= rho0*(-ci*rkph)
            fmc(2,2,3,inod)= rho0*(  -rkth0)

            fmc(1,1,4,inod)= rho0*(1.d-6)
            fmc(2,2,4,inod)= rho0*(-1.d0)
            fmc(3,3,4,inod)= rho0*(-1.d0)
         enddo

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=0.5d0*(-11*fmc(i,j,k,1)+18*fmc(i,j,k,2)
     &                                - 9*fmc(i,j,k,3)+ 2*fmc(i,j,k,4))
                  fmd(i,j,k,3)=fmc(i,j,k,4)
                  fmd(i,j,k,4)=0.5d0*(- 2*fmc(i,j,k,1)+ 9*fmc(i,j,k,2)
     &                                -18*fmc(i,j,k,3)+11*fmc(i,j,k,4))
               enddo
            enddo
         enddo

         do j=1,3
         do i=1,3
            ml=6*(nr-1)+2*(i-1)+1
            mw=mc+2*(j-1)-2*(i-1)
            do inod=1,4

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,1)*drho
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,1)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,5)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,2)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,2)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,6)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,6)
            fma(mw+6,ml  )=fma(mw+6,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,3)*drho
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,3)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,7)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,7)/drho
            fma(mw+7,ml  )=fma(mw+7,ml  )
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,4)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,4)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,8)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,1)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,1)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,5)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,2)*drho**3
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,2)*drho**2
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,6)*drho**2
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,6)*drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,3)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,3)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,7)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,7)
            fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,4)*drho**3
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,4)*drho**2
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,8)*drho**2
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,8)*drho

            fma(mw-6,ml+6)=fma(mw-6,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,1)*drho
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,1)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,5)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,5)/drho
            fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,2)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,2)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,6)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,3)*drho
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,3)
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,7)
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,4)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,4)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,8)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,8)

            fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,1)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,1)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,5)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,5)
            fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,2)*drho**3
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,2)*drho**2
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,6)*drho**2
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,3)*drho**2
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,3)*drho
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,7)*drho
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,4)*drho**3
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,4)*drho**2
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,8)*drho**2
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,8)*drho

               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
         enddo
         fma(mc,3)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc-2,3)=1.d0
         fma(mc,3)=ci*nth
         fma(mc,5)=1.d0
      else
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
c$$$            write(6,'A,I5,1P4E12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
c$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(6*(nr-1)+3)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(6* nr   +3)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(6*(nr-1)+5)=(rho(nr+1)-0.85d0)*angl
            fvb(6* nr   +5)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_r4

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r5(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,mr
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4,4):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer,dimension(3),parameter :: ishift=(/0,2,4/)

      call mesh_init(nrmax,npow)

      nvmax=6
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         factor=rf**2
         rkph=nph
      
         do inod=1,4
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,4
            if(inod.eq.1) then
               if(nr.eq.1) then
!                  rho0=(8.d0*rho(nr)+rho(nr+1))/9.d0
                  rho0=rho(nr+1))/1.d-6
               else
                  rho0=rho(nr)
               endif
            elseif(inod.eq.2) then
               rho0=(2.d0*rho(nr)+rho(nr+1))/3.d0
            elseif(inod.eq.3) then
               rho0=(rho(nr)+2.d0*rho(nr+1))/3.d0
            else
               rho0=rho(nr+1)
            endif

            if(nr.eq.1.and.inod.eq.1) then
            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkph**2-rkth**2)
            fmc(1,2,1,inod)= rho0*(-ci*rkth*rkth0)
            fmc(2,1,1,inod)= rho0*(+ci*rkth*rkth0)
            fmc(2,2,1,inod)= rho0*(factor-rkph**2-rkth0**2)
            fmc(2,3,1,inod)= rho0*(rkth*rkph)
            fmc(3,2,1,inod)= rho0*(rkth*rkph)
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)= rho0*( ci*rkth)
            fmc(2,2,2,inod)= rho0*(  -rkth0)
            fmc(3,1,2,inod)= 0.d0

            fmc(1,2,3,inod)= rho0*(-ci*rkth)
            fmc(1,3,3,inod)= 0.d0
            fmc(2,2,3,inod)= rho0*(  -rkth0)

            fmc(2,2,4,inod)= 0.d0
            fmc(3,3,4,inod)= 0.d0

            else
            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkph**2-rkth**2)
            fmc(1,2,1,inod)= rho0*(-ci*rkth*rkth0)
            fmc(2,1,1,inod)= rho0*(+ci*rkth*rkth0)
            fmc(2,2,1,inod)= rho0*(factor-rkph**2-rkth0**2)
            fmc(2,3,1,inod)= rho0*(rkth*rkph)
            fmc(3,2,1,inod)= rho0*(rkth*rkph)
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)= rho0*( ci*rkth)
            fmc(2,2,2,inod)= rho0*(  -rkth0)
            fmc(3,1,2,inod)= rho0*( ci*rkph)

            fmc(1,2,3,inod)= rho0*(-ci*rkth)
            fmc(1,3,3,inod)= rho0*(-ci*rkph)
            fmc(2,2,3,inod)= rho0*(  -rkth0)

            fmc(2,2,4,inod)= rho0*(-1.d0)
            fmc(3,3,4,inod)= rho0*(-1.d0)
            endif
         enddo

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
!                  fmd(i,j,k,2)=0.5d0*(-11*fmc(i,j,k,1)+18*fmc(i,j,k,2)
!     &                                - 9*fmc(i,j,k,3)+ 2*fmc(i,j,k,4))
                  fmd(i,j,k,2)=1.5d0*(-3*fmc(i,j,k,1)+4*fmc(i,j,k,2)
     &                                -  fmc(i,j,k,3))
                  fmd(i,j,k,3)=fmc(i,j,k,4)
!                  fmd(i,j,k,4)=0.5d0*(- 2*fmc(i,j,k,1)+ 9*fmc(i,j,k,2)
!     &                                -18*fmc(i,j,k,3)+11*fmc(i,j,k,4))
                  fmd(i,j,k,4)=1.5d0*(   fmc(i,j,k,2)
     &                                -4*fmc(i,j,k,3)+3*fmc(i,j,k,4))
               enddo
            enddo
         enddo

         do j=1,3
         do i=1,3
            ml=6*(nr-1)+ishift(i)+1
            mw=mc+ishift(j)-ishift(i)
            mr=6
            do inod=1,4

! (1,1) **********
            if(i.eq.1.and.j.eq.1) then

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_hqq(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,4,1)
     &              +fmd(i,j,3,inod)*table_hqq(inod,1,4)
     &              +fmd(i,j,4,inod)*table_hqq(inod,4,4)/drho

            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,inod)*table_hqq(inod,1,2)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,4,2)
     &              +fmd(i,j,3,inod)*table_hqq(inod,1,5)
     &              +fmd(i,j,4,inod)*table_hqq(inod,4,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_hqq(inod,1,3)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,4,3)
     &              +fmd(i,j,3,inod)*table_hqq(inod,1,6)
     &              +fmd(i,j,4,inod)*table_hqq(inod,4,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hqq(inod,2,1)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,5,1)
     &              +fmd(i,j,3,inod)*table_hqq(inod,2,4)
     &              +fmd(i,j,4,inod)*table_hqq(inod,5,4)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,inod)*table_hqq(inod,2,2)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,5,2)
     &              +fmd(i,j,3,inod)*table_hqq(inod,2,5)
     &              +fmd(i,j,4,inod)*table_hqq(inod,5,5)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hqq(inod,2,3)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,5,3)
     &              +fmd(i,j,3,inod)*table_hqq(inod,2,6)
     &              +fmd(i,j,4,inod)*table_hqq(inod,5,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqq(inod,3,1)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,6,1)
     &              +fmd(i,j,3,inod)*table_hqq(inod,3,4)
     &              +fmd(i,j,4,inod)*table_hqq(inod,6,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqq(inod,3,2)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,6,2)
     &              +fmd(i,j,3,inod)*table_hqq(inod,3,5)
     &              +fmd(i,j,4,inod)*table_hqq(inod,6,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqq(inod,3,3)*drho
     &              +fmd(i,j,2,inod)*table_hqq(inod,6,3)
     &              +fmd(i,j,3,inod)*table_hqq(inod,3,6)
     &              +fmd(i,j,4,inod)*table_hqq(inod,6,6)/drho

! (1,*) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,1)
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,5)
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,2)*drho
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,6)*drho
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,3)*drho
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,3)
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,7)
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  )
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,4)*drho**2
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,4)*drho
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,8)*drho
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,1)*drho
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,1)
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,5)
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,5)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,2)*drho
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,6)*drho
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,6)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,3)*drho
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,3)
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,7)
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,7)/drho
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,4)*drho**2
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,4)*drho
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,8)*drho
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,8)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,1)*drho
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,1)
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,5)
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,2)*drho
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,6)*drho
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,3)*drho
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,3)
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,7)
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,4)*drho**2
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,4)*drho
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,8)*drho
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,8)

! (*,1) **********
            else if(j.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_hhq(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_hhq(inod,5,1)
     &              +fmd(i,j,3,inod)*table_hhq(inod,1,4)
     &              +fmd(i,j,4,inod)*table_hhq(inod,5,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,inod)*table_hhq(inod,1,2)*drho
     &              +fmd(i,j,2,inod)*table_hhq(inod,5,2)
     &              +fmd(i,j,3,inod)*table_hhq(inod,1,5)
     &              +fmd(i,j,4,inod)*table_hhq(inod,5,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_hhq(inod,1,3)*drho
     &              +fmd(i,j,2,inod)*table_hhq(inod,5,3)
     &              +fmd(i,j,3,inod)*table_hhq(inod,1,6)
     &              +fmd(i,j,4,inod)*table_hhq(inod,5,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hhq(inod,2,1)*drho**2
     &              +fmd(i,j,2,inod)*table_hhq(inod,6,1)*drho
     &              +fmd(i,j,3,inod)*table_hhq(inod,2,4)*drho
     &              +fmd(i,j,4,inod)*table_hhq(inod,6,4)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,inod)*table_hhq(inod,2,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hhq(inod,6,2)*drho
     &              +fmd(i,j,3,inod)*table_hhq(inod,2,5)*drho
     &              +fmd(i,j,4,inod)*table_hhq(inod,6,5)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hhq(inod,2,3)*drho**2
     &              +fmd(i,j,2,inod)*table_hhq(inod,6,3)*drho
     &              +fmd(i,j,3,inod)*table_hhq(inod,2,6)*drho
     &              +fmd(i,j,4,inod)*table_hhq(inod,6,6)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhq(inod,3,1)*drho
     &              +fmd(i,j,2,inod)*table_hhq(inod,7,1)
     &              +fmd(i,j,3,inod)*table_hhq(inod,3,4)
     &              +fmd(i,j,4,inod)*table_hhq(inod,7,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhq(inod,3,2)*drho
     &              +fmd(i,j,2,inod)*table_hhq(inod,7,2)
     &              +fmd(i,j,3,inod)*table_hhq(inod,3,5)
     &              +fmd(i,j,4,inod)*table_hhq(inod,7,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhq(inod,3,3)*drho
     &              +fmd(i,j,2,inod)*table_hhq(inod,7,3)
     &              +fmd(i,j,3,inod)*table_hhq(inod,3,6)
     &              +fmd(i,j,4,inod)*table_hhq(inod,7,6)/drho

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhq(inod,4,1)*drho**2
     &              +fmd(i,j,2,inod)*table_hhq(inod,8,1)*drho
     &              +fmd(i,j,3,inod)*table_hhq(inod,4,4)*drho
     &              +fmd(i,j,4,inod)*table_hhq(inod,8,4)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhq(inod,4,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hhq(inod,8,2)*drho
     &              +fmd(i,j,3,inod)*table_hhq(inod,4,5)*drho
     &              +fmd(i,j,4,inod)*table_hhq(inod,8,5)
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhq(inod,4,3)*drho**2
     &              +fmd(i,j,2,inod)*table_hhq(inod,8,3)*drho
     &              +fmd(i,j,3,inod)*table_hhq(inod,4,6)*drho
     &              +fmd(i,j,4,inod)*table_hhq(inod,8,6)


! (*,*) **********
            else
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  )
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,6)*drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,7)
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,8)*drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,8)

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,5)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+mr+1)=fma(mw  ,ml+mr+1)
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,4)*drho**3
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,4)*drho**2
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,8)*drho**2
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,8)*drho
            endif

               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
         enddo
         fma(mc,3)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc-2,3)=1.d0
         fma(mc,3)=ci*nth
         fma(mc,5)=1.d0
      elseif(abs(nth).eq.2) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
      else
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
      endif

      do mw=1,mwmax
!         fma(mw,mlmax-5) = 0.d0
         fma(mw,mlmax-4) = 0.d0
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
!      fma(mc,mlmax-5) = 1.d0
      fma(mc,mlmax-4) = 1.d0
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
            fvb(6*(nr-1)+3)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(6* nr   +3)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(6*(nr-1)+5)=(rho(nr+1)-0.85d0)*angl
            fvb(6* nr   +5)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_r5

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r6(nrmax,npow,nth,nph,rf,angl)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,mr
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4,2):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer,dimension(3),parameter :: ishift=(/0,1,3/)

      call mesh_init(nrmax,npow)

      nvmax=5
     
      call fem_init(nrmax,nvmax)

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         factor=rf**2
         rkph=nph
      
         do inod=1,2
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,2
            if(inod.eq.1) then
               if(nr.eq.1) then
                  rho0=(3.d0*rho(nr)+rho(nr+1))/4.d0
               else
                  rho0=rho(nr)
               endif
            else
               rho0=rho(nr+1)
            endif

            rkth=nth/rho0
            rkth0=1.d0/rho0

            fmc(1,1,1,inod)= rho0*(factor-rkph**2-rkth**2)
            fmc(1,2,1,inod)= rho0*(-ci*rkth*rkth0)
            fmc(2,1,1,inod)= rho0*(+ci*rkth*rkth0)
            fmc(2,2,1,inod)= rho0*(factor-rkph**2-rkth0**2)
            fmc(2,3,1,inod)= rho0*(rkth*rkph)
            fmc(3,2,1,inod)= rho0*(rkth*rkph)
            fmc(3,3,1,inod)= rho0*(factor-rkth**2)
      
            fmc(2,1,2,inod)= rho0*( ci*rkth)
            fmc(2,2,2,inod)= rho0*(  -rkth0)
            fmc(3,1,2,inod)= rho0*( ci*rkph)

            fmc(1,2,3,inod)= rho0*(-ci*rkth)
            fmc(1,3,3,inod)= rho0*(-ci*rkph)
            fmc(2,2,3,inod)= rho0*(  -rkth0)

            fmc(2,2,4,inod)= rho0*(-1.d0)
            fmc(3,3,4,inod)= rho0*(-1.d0)
         enddo

         do k=1,4
            do j=1,3
               do i=1,3
                  fmd(i,j,k,1)=fmc(i,j,k,1)
                  fmd(i,j,k,2)=fmc(i,j,k,2)
               enddo
            enddo
         enddo

         do j=1,3
         do i=1,3
            ml=5*(nr-1)+ishift(i)+1
            mw=mc+ishift(j)-ishift(i)
            mr=5
            do inod=1,2

! (1,1) **********
            if(i.eq.1.and.j.eq.1) then

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_lll(inod,3,1)
     &              +fmd(i,j,3,inod)*table_lll(inod,1,3)
     &              +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
     &              +fmd(i,j,2,inod)*table_lll(inod,3,2)
     &              +fmd(i,j,3,inod)*table_lll(inod,1,4)
     &              +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
     &              +fmd(i,j,2,inod)*table_lll(inod,4,1)
     &              +fmd(i,j,3,inod)*table_lll(inod,2,3)
     &              +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
     &              +fmd(i,j,2,inod)*table_lll(inod,4,2)
     &              +fmd(i,j,3,inod)*table_lll(inod,2,4)
     &              +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho

! (1,*) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_llq(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_llq(inod,3,1)
     &              +fmd(i,j,3,inod)*table_llq(inod,1,4)
     &              +fmd(i,j,4,inod)*table_llq(inod,3,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,inod)*table_llq(inod,1,2)*drho
     &              +fmd(i,j,2,inod)*table_llq(inod,3,2)
     &              +fmd(i,j,3,inod)*table_llq(inod,1,5)
     &              +fmd(i,j,4,inod)*table_llq(inod,3,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_llq(inod,1,3)*drho
     &              +fmd(i,j,2,inod)*table_llq(inod,3,3)
     &              +fmd(i,j,3,inod)*table_llq(inod,1,6)
     &              +fmd(i,j,4,inod)*table_llq(inod,3,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_llq(inod,2,1)*drho
     &              +fmd(i,j,2,inod)*table_llq(inod,4,1)
     &              +fmd(i,j,3,inod)*table_llq(inod,2,4)
     &              +fmd(i,j,4,inod)*table_llq(inod,4,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_llq(inod,2,2)*drho
     &              +fmd(i,j,2,inod)*table_llq(inod,4,2)
     &              +fmd(i,j,3,inod)*table_llq(inod,2,5)
     &              +fmd(i,j,4,inod)*table_llq(inod,4,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_llq(inod,2,3)*drho
     &              +fmd(i,j,2,inod)*table_llq(inod,4,3)
     &              +fmd(i,j,3,inod)*table_llq(inod,2,6)
     &              +fmd(i,j,4,inod)*table_llq(inod,4,6)/drho

! (*,1) **********
            else if(j.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_lql(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_lql(inod,4,1)
     &              +fmd(i,j,3,inod)*table_lql(inod,1,3)
     &              +fmd(i,j,4,inod)*table_lql(inod,4,3)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_lql(inod,1,2)*drho
     &              +fmd(i,j,2,inod)*table_lql(inod,4,2)
     &              +fmd(i,j,3,inod)*table_lql(inod,1,4)
     &              +fmd(i,j,4,inod)*table_lql(inod,4,4)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,inod)*table_lql(inod,2,1)*drho
     &              +fmd(i,j,2,inod)*table_lql(inod,5,1)
     &              +fmd(i,j,3,inod)*table_lql(inod,2,3)
     &              +fmd(i,j,4,inod)*table_lql(inod,5,3)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,inod)*table_lql(inod,2,2)*drho
     &              +fmd(i,j,2,inod)*table_lql(inod,5,2)
     &              +fmd(i,j,3,inod)*table_lql(inod,2,4)
     &              +fmd(i,j,4,inod)*table_lql(inod,5,4)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_lql(inod,3,1)*drho
     &              +fmd(i,j,2,inod)*table_lql(inod,6,1)
     &              +fmd(i,j,3,inod)*table_lql(inod,3,3)
     &              +fmd(i,j,4,inod)*table_lql(inod,6,3)/drho
            fma(mw,ml+mr)=fma(mw,ml+mr)
     &              +fmd(i,j,1,inod)*table_lql(inod,3,2)*drho
     &              +fmd(i,j,2,inod)*table_lql(inod,6,2)
     &              +fmd(i,j,3,inod)*table_lql(inod,3,4)
     &              +fmd(i,j,4,inod)*table_lql(inod,6,4)/drho

! (*,*) **********
            else
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,inod)*table_lqq(inod,1,1)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,4,1)
     &              +fmd(i,j,3,inod)*table_lqq(inod,1,4)
     &              +fmd(i,j,4,inod)*table_lqq(inod,4,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,inod)*table_lqq(inod,1,2)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,4,2)
     &              +fmd(i,j,3,inod)*table_lqq(inod,1,5)
     &              +fmd(i,j,4,inod)*table_lqq(inod,4,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,inod)*table_lqq(inod,1,3)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,4,3)
     &              +fmd(i,j,3,inod)*table_lqq(inod,1,6)
     &              +fmd(i,j,4,inod)*table_lqq(inod,4,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,inod)*table_lqq(inod,2,1)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,5,1)
     &              +fmd(i,j,3,inod)*table_lqq(inod,2,4)
     &              +fmd(i,j,4,inod)*table_lqq(inod,5,4)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,inod)*table_lqq(inod,2,2)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,5,2)
     &              +fmd(i,j,3,inod)*table_lqq(inod,2,5)
     &              +fmd(i,j,4,inod)*table_lqq(inod,5,5)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,inod)*table_lqq(inod,2,3)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,5,3)
     &              +fmd(i,j,3,inod)*table_lqq(inod,2,6)
     &              +fmd(i,j,4,inod)*table_lqq(inod,5,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,inod)*table_lqq(inod,3,1)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,6,1)
     &              +fmd(i,j,3,inod)*table_lqq(inod,3,4)
     &              +fmd(i,j,4,inod)*table_lqq(inod,6,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,inod)*table_lqq(inod,3,2)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,6,2)
     &              +fmd(i,j,3,inod)*table_lqq(inod,3,5)
     &              +fmd(i,j,4,inod)*table_lqq(inod,6,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,inod)*table_lqq(inod,3,3)*drho
     &              +fmd(i,j,2,inod)*table_lqq(inod,6,3)
     &              +fmd(i,j,3,inod)*table_lqq(inod,3,6)
     &              +fmd(i,j,4,inod)*table_lqq(inod,6,6)/drho

            endif

               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,2) = 0.d0
         enddo
         fma(mc,2)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,2) = 0.d0
            fma(mw,4) = 0.d0
         enddo
         fma(mc-1,2)=1.d0
         fma(mc,2)=ci*nth
         fma(mc,4)=1.d0
      else
         do mw=1,mwmax
            fma(mw,2) = 0.d0
            fma(mw,4) = 0.d0
         enddo
         fma(mc,2)=1.d0
         fma(mc,4)=1.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-2) = 0.d0
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-2) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
            fvb(5*(nr-1)+2)=(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(5* nr   +2)=(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(5*(nr-1)+4)=(rho(nr+1)-0.85d0)*angl
            fvb(5* nr   +4)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

c$$$      do ml=1,mlmax
c$$$         do mw=1,mwmax
c$$$            if(abs(fma(mw,ml)).ne.0.d0) then
c$$$               write(6,'(2I5,1P2E12.4)') ml,mw,fma(mw,ml)
c$$$            endif
c$$$         enddo
c$$$      enddo
      return
      end subroutine fem_calc_r6

      include 'femtest-sub.f'

      end program femtest
