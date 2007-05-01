!     $Id$

      program femtest

      implicit none
      integer:: mwmax,mlmax
      complex(8),dimension(:,:),allocatable:: fma
      complex(8),dimension(:),allocatable:: fvb,fvx
      real(8),dimension(:),allocatable:: rho,f1,f2

      integer:: id,nrmax,npow,mw,ml,nr,ierr,nth,nph
      real(8):: rf

      call gsopen

      id=3
      nrmax=11
      npow=1
      nth=0
      nph=0
      rf=1.d0

    1 write(6,'(A,3I5)') '## input: id,nrmax,npow=',id,nrmax,npow,nth,nph,rf
      read(5,*,err=1,end=9000) id,nrmax,npow,nth,nph,rf

      select case(id)
         case(1)
            call fem_calc1(nrmax,npow)
         case(2)
            call fem_calc2(nrmax,npow)
         case(3)
            call fem_calc3(nrmax,npow,nth,nph,rf)
      end select

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
      goto 1

 9000 call gsclos
      stop

      contains

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
         if(allocated(f1)) deallocate(f1)
         if(allocated(f2)) deallocate(f2)
         allocate(rho(nrmax))
         allocate(f1(nrmax))
         allocate(f2(nrmax))
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

      subroutine fem_calc1(nrmax,npow)

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
      write(6,'(A,3I5)') 'mlmax,mwmax,mc=',mlmax,mwmax,mc

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
      end subroutine fem_calc1

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc2(nrmax,npow)

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
      write(6,'(A,3I5)') 'mlmax,mwmax,mc=',mlmax,mwmax,mc

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

      return
      end subroutine fem_calc2

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc3(nrmax,npow,nth,nph,rf)

      use libfem_mod
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph
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

!      write(6,'(A,3I5)') 'mlmax,mwmax,mc=',mlmax,mwmax,mc

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

      do inod=1,4
         fmc(1,1,1,inod)= 1.d0-rkth**2+rkph**2
         fmc(2,2,1,inod)= 1.d0-rkph**2
         fmc(2,3,1,inod)=      rkth*rkph
         fmc(3,2,1,inod)=      rkth*rkph
         fmc(3,3,1,inod)= 1.d0-rkth**2
      
         fmc(2,1,2,inod)= ci*rkth
         fmc(3,1,2,inod)= ci*rkph

         fmc(1,2,3,inod)=-ci*rkth
         fmc(1,3,3,inod)=-ci*rkph

         fmc(2,2,4,inod)= -1.d0
         fmc(3,3,4,inod)= -1.d0
      enddo

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)
         ml=6*(nr-1)+1

         do j=1,3
            do i=1,3
               do k=1,4
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
               mw=mc+2*(j-1)-2*(i-1)

! ----- non derivative terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,k,1)*table_hhh(k,1,1)*drho
               fma(mw+1,ml  )=fma(mw+1,ml  )
     &                       +fmd(i,j,k,1)*table_hhh(k,2,1)*drho*drho
               fma(mw+6,ml  )=fma(mw+6,ml  )
     &                       +fmd(i,j,k,1)*table_hhh(k,3,1)*drho
               fma(mw+7,ml  )=fma(mw+7,ml  )
     &                       +fmd(i,j,k,1)*table_hhh(k,4,1)*drho*drho

               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,k,1)*table_hhh(k,1,2)*drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                       +fmd(i,j,k,1)*table_hhh(k,2,2)*drho*drho
               fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                       +fmd(i,j,k,1)*table_hhh(k,3,2)*drho
               fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                       +fmd(i,j,k,1)*table_hhh(k,4,2)*drho*drho

               fma(mw-6,ml+6)=fma(mw-6,ml+6 )
     &                       +fmd(i,j,k,1)*table_hhh(k,1,3)*drho
               fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                       +fmd(i,j,k,1)*table_hhh(k,2,3)*drho*drho
               fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                       +fmd(i,j,k,1)*table_hhh(k,3,3)*drho
               fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                       +fmd(i,j,k,1)*table_hhh(k,4,3)*drho*drho

               fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                       +fmd(i,j,k,1)*table_hhh(k,1,4)*drho
               fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                       +fmd(i,j,k,1)*table_hhh(k,2,4)*drho*drho
               fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                       +fmd(i,j,k,1)*table_hhh(k,3,4)*drho
               fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                       +fmd(i,j,k,1)*table_hhh(k,4,4)*drho*drho

! ----- dw/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,k,2)*table_hhh(k,5,1)*drho
               fma(mw+1,ml  )=fma(mw+1,ml  )
     &                       +fmd(i,j,k,2)*table_hhh(k,6,1)*drho*drho
               fma(mw+6,ml  )=fma(mw+6,ml  )
     &                       +fmd(i,j,k,2)*table_hhh(k,7,1)*drho
               fma(mw+7,ml  )=fma(mw+7,ml  )
     &                       +fmd(i,j,k,2)*table_hhh(k,8,1)*drho*drho

               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,k,2)*table_hhh(k,5,2)*drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                       +fmd(i,j,k,2)*table_hhh(k,6,2)*drho*drho
               fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                       +fmd(i,j,k,2)*table_hhh(k,7,2)*drho
               fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                       +fmd(i,j,k,2)*table_hhh(k,8,2)*drho*drho

               fma(mw-6,ml+6)=fma(mw-6,ml+6 )
     &                       +fmd(i,j,k,2)*table_hhh(k,5,3)*drho
               fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                       +fmd(i,j,k,2)*table_hhh(k,6,3)*drho*drho
               fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                       +fmd(i,j,k,2)*table_hhh(k,7,3)*drho
               fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                       +fmd(i,j,k,2)*table_hhh(k,8,3)*drho*drho

               fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                       +fmd(i,j,k,2)*table_hhh(k,5,4)*drho
               fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                       +fmd(i,j,k,2)*table_hhh(k,6,4)*drho*drho
               fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                       +fmd(i,j,k,2)*table_hhh(k,7,4)*drho
               fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                       +fmd(i,j,k,2)*table_hhh(k,8,4)*drho*drho

! ----- dE/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,k,3)*table_hhh(k,1,5)*drho
               fma(mw+1,ml  )=fma(mw+1,ml  )
     &                       +fmd(i,j,k,3)*table_hhh(k,2,5)*drho*drho
               fma(mw+6,ml  )=fma(mw+6,ml  )
     &                       +fmd(i,j,k,3)*table_hhh(k,3,5)*drho
               fma(mw+7,ml  )=fma(mw+7,ml  )
     &                       +fmd(i,j,k,3)*table_hhh(k,4,5)*drho*drho

               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,k,3)*table_hhh(k,1,6)*drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                       +fmd(i,j,k,3)*table_hhh(k,2,6)*drho*drho
               fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                       +fmd(i,j,k,3)*table_hhh(k,3,6)*drho
               fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                       +fmd(i,j,k,3)*table_hhh(k,4,6)*drho*drho

               fma(mw-6,ml+6)=fma(mw-6,ml+6 )
     &                       +fmd(i,j,k,3)*table_hhh(k,1,7)*drho
               fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                       +fmd(i,j,k,3)*table_hhh(k,2,7)*drho*drho
               fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                       +fmd(i,j,k,3)*table_hhh(k,3,7)*drho
               fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                       +fmd(i,j,k,3)*table_hhh(k,4,7)*drho*drho

               fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,1,8)*drho
               fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,2,8)*drho*drho
               fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,3,8)*drho
               fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,4,8)*drho*drho

! ----- dw/dx dE/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  )
     &                       +fmd(i,j,k,4)*table_hhh(k,5,5)*drho
               fma(mw+1,ml  )=fma(mw+1,ml  )
     &                       +fmd(i,j,k,4)*table_hhh(k,6,5)*drho*drho
               fma(mw+6,ml  )=fma(mw+6,ml  )
     &                       +fmd(i,j,k,4)*table_hhh(k,7,5)*drho
               fma(mw+7,ml  )=fma(mw+7,ml  )
     &                       +fmd(i,j,k,4)*table_hhh(k,8,5)*drho*drho

               fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &                       +fmd(i,j,k,4)*table_hhh(k,5,6)*drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &                       +fmd(i,j,k,4)*table_hhh(k,6,6)*drho*drho
               fma(mw+5,ml+1)=fma(mw+5,ml+1)
     &                       +fmd(i,j,k,4)*table_hhh(k,7,6)*drho
               fma(mw+6,ml+1)=fma(mw+6,ml+1)
     &                       +fmd(i,j,k,4)*table_hhh(k,8,6)*drho*drho

               fma(mw-6,ml+6)=fma(mw-6,ml+6 )
     &                       +fmd(i,j,k,4)*table_hhh(k,5,7)*drho
               fma(mw-5,ml+6)=fma(mw-5,ml+6)
     &                       +fmd(i,j,k,4)*table_hhh(k,6,7)*drho*drho
               fma(mw  ,ml+6)=fma(mw  ,ml+6)
     &                       +fmd(i,j,k,4)*table_hhh(k,7,7)*drho
               fma(mw+1,ml+6)=fma(mw+1,ml+6)
     &                       +fmd(i,j,k,4)*table_hhh(k,8,7)*drho*drho

               fma(mw-7,ml+7)=fma(mw-7,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,5,8)*drho
               fma(mw-6,ml+7)=fma(mw-6,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,6,8)*drho*drho
               fma(mw-1,ml+7)=fma(mw-1,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,7,8)*drho
               fma(mw  ,ml+7)=fma(mw  ,ml+7)
     &                       +fmd(i,j,k,3)*table_hhh(k,8,8)*drho*drho
            enddo
         enddo
      enddo

      do mw=1,mwmax
         fma(mw,1) = 0.d0
         fma(mw,3) = 0.d0
         fma(mw,5) = 0.d0
      enddo
      fma(mc,1)=1.d0
      fma(mc,3)=1.d0
      fma(mc,5)=1.d0

      do mw=1,mwmax
         fma(mw,mlmax-5) = 0.d0
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
      fma(mc,mlmax-5) = 1.d0
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0

      fvb(mlmax-10)=1.d0

      return
      end subroutine fem_calc3

      include 'femtest-sub.f'

      end program femtest
