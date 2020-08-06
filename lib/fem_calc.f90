MODULE fem_calc

!----- calculate coefficint matrix fma and source vector fvb -----

   USE libtestfem

CONTAINS
   
      subroutine fem_calc_l1(nrmax,npow)

      use libfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer:: nr,ml,mw,mc,nvmax,i,j,k
      real(8):: drho,val1,val2,val3,val4,temp
      real(8),dimension(4):: coef

      call mesh_init(nrmax,npow)

      nvmax=2
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector

      call fem_init

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

      use libfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer:: nr,ml,mw,mc,nvmax,i,j,k
      real(8):: drho,val1,val2,val3,val4,temp
      real(8),dimension(4):: coef

      call mesh_init(nrmax,npow)

      nvmax=2
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector

      call fem_init

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
!                  write(6,'(2I5,4ES12.4)') j,k,temp,
!     &                 table_hhh(1,j,k),table_hhh(3,j,k),table_hh(j,k)
!               enddo
!            enddo
!            endif

            do i=1,4
            fma(mc  ,ml  )=fma(mc  ,ml  ) &
                          +coef(i)*table_hhh(i,5,5)/drho
            fma(mc+1,ml  )=fma(mc+1,ml  ) &
                          +coef(i)*table_hhh(i,6,5)/drho*drho
            fma(mc+2,ml  )=fma(mc+2,ml  ) &
                          +coef(i)*table_hhh(i,7,5)/drho
            fma(mc+3,ml  )=fma(mc+3,ml  ) &
                          +coef(i)*table_hhh(i,8,5)/drho*drho

            fma(mc-1,ml+1)=fma(mc-1,ml+1) &
                          +coef(i)*table_hhh(i,5,6)/drho*drho
            fma(mc  ,ml+1)=fma(mc  ,ml+1) &
                          +coef(i)*table_hhh(i,6,6)/drho*drho*drho
            fma(mc+1,ml+1)=fma(mc+1,ml+1) &
                          +coef(i)*table_hhh(i,7,6)/drho*drho
            fma(mc+2,ml+1)=fma(mc+2,ml+1) &
                          +coef(i)*table_hhh(i,8,6)/drho*drho*drho

            fma(mc-2,ml+2)=fma(mc-2,ml+2) &
                          +coef(i)*table_hhh(i,5,7)/drho
            fma(mc-1,ml+2)=fma(mc-1,ml+2) &
                          +coef(i)*table_hhh(i,6,7)/drho*drho
            fma(mc  ,ml+2)=fma(mc  ,ml+2) &
                          +coef(i)*table_hhh(i,7,7)/drho
            fma(mc+1,ml+2)=fma(mc+1,ml+2) &
                          +coef(i)*table_hhh(i,8,7)/drho*drho

            fma(mc-3,ml+3)=fma(mc-3,ml+3) &
                          +coef(i)*table_hhh(i,5,8)/drho*drho
            fma(mc-2,ml+3)=fma(mc-2,ml+3) &
                          +coef(i)*table_hhh(i,6,8)/drho*drho*drho
            fma(mc-1,ml+3)=fma(mc-1,ml+3) &
                          +coef(i)*table_hhh(i,7,8)/drho*drho
            fma(mc  ,ml+3)=fma(mc  ,ml+3) &
                          +coef(i)*table_hhh(i,8,8)/drho*drho*drho
         enddo
         enddo
         do mw=1,mwmax
            fma(mw,1) = 0.d0
         enddo
         fma(mc,1)=1.d0
         do mw=1,mwmax
            fma(mw,mlmax-5) = 0.d0
         enddo
         fma(mc,mlmax-5) = 1.d0

         fvb(1)=1.d0
         fvb(mlmax-5)=1.d0/11.d0

      return
      end subroutine fem_calc_l2

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_l3(nrmax,npow)

      use libfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer:: nr,ml,mw,mc,nvmax,i,j,k
      real(8):: drho,val1,val2,val3,val4,temp
      real(8),dimension(4):: coef
      real(8):: rho1,rho2,rho3,rho4
      real(8):: moment0,moment1,moment2,moment3,xi,yi
      real(8),parameter:: m0=4.D0
      real(8),parameter:: m2=5.D0/16.D0
      real(8),parameter:: m4=41.D0/1024.D0
      real(8),parameter:: m6=365.D0/65536.D0

      call mesh_init(nrmax,npow)

      nvmax=2
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector

      call fem_init

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo
      mc=(mwmax+1)/2

         do nr=1,nrmax
            rho(nr)=rho(nr)+0.1D0
         end do

         do nr=1,nrmax-1
            drho=rho(nr+1)-rho(nr)
            ml=2*(nr-1)+1

            rho1=(7.D0*rho(nr)+     rho(nr+1))/8.D0
            rho2=(5.D0*rho(nr)+3.D0*rho(nr+1))/8.D0
            rho3=(3.D0*rho(nr)+5.D0*rho(nr+1))/8.D0
            rho4=(     rho(nr)+7.D0*rho(nr+1))/8.D0

!            coef(1)=(0.1d0+rho1)**2
!            coef(2)=(0.1d0+rho2)**2
!            coef(3)=(0.1d0+rho3)**2
!            coef(4)=(0.1d0+rho4)**2
            coef(1)=rho1**2
            coef(2)=rho2**2
            coef(3)=rho3**2
            coef(4)=rho4**2

            moment0=0.D0
            moment1=0.D0
            moment2=0.D0
            moment3=0.D0
            DO i=1,4
               xi=0.125D0*(2*i-1)
               yi=xi-0.5D0
               moment0=moment0+      coef(i)
               moment1=moment1+yi   *coef(i)
               moment2=moment2+yi**2*coef(i)
               moment3=moment3+yi**3*coef(i)
            END DO
            coef(1)=(m4*moment0-m2*moment2)/(m0*m4-m2**2)
            coef(2)=(m6*moment1-m4*moment3)/(m2*m6-m4**2)
            coef(3)=(m2*moment0-m0*moment2)/(m2**2-m0*m4)*2.D0
            coef(4)=(m4*moment1-m2*moment3)/(m4**2-m2*m6)*6.D0
            
            do i=1,4
            fma(mc  ,ml  )=fma(mc  ,ml  ) &
                          +coef(i)*table_hhg(5,5,i)/drho
            fma(mc+1,ml  )=fma(mc+1,ml  ) &
                          +coef(i)*table_hhg(6,5,i)/drho*drho
            fma(mc+2,ml  )=fma(mc+2,ml  ) &
                          +coef(i)*table_hhg(7,5,i)/drho
            fma(mc+3,ml  )=fma(mc+3,ml  ) &
                          +coef(i)*table_hhg(8,5,i)/drho*drho

            fma(mc-1,ml+1)=fma(mc-1,ml+1) &
                          +coef(i)*table_hhg(5,6,i)/drho*drho
            fma(mc  ,ml+1)=fma(mc  ,ml+1) &
                          +coef(i)*table_hhg(6,6,i)/drho*drho*drho
            fma(mc+1,ml+1)=fma(mc+1,ml+1) &
                          +coef(i)*table_hhg(7,6,i)/drho*drho
            fma(mc+2,ml+1)=fma(mc+2,ml+1) &
                          +coef(i)*table_hhg(8,6,i)/drho*drho*drho

            fma(mc-2,ml+2)=fma(mc-2,ml+2) &
                          +coef(i)*table_hhg(5,7,i)/drho
            fma(mc-1,ml+2)=fma(mc-1,ml+2) &
                          +coef(i)*table_hhg(6,7,i)/drho*drho
            fma(mc  ,ml+2)=fma(mc  ,ml+2) &
                          +coef(i)*table_hhg(7,7,i)/drho
            fma(mc+1,ml+2)=fma(mc+1,ml+2) &
                          +coef(i)*table_hhg(8,7,i)/drho*drho

            fma(mc-3,ml+3)=fma(mc-3,ml+3) &
                          +coef(i)*table_hhg(5,8,i)/drho*drho
            fma(mc-2,ml+3)=fma(mc-2,ml+3) &
                          +coef(i)*table_hhg(6,8,i)/drho*drho*drho
            fma(mc-1,ml+3)=fma(mc-1,ml+3) &
                          +coef(i)*table_hhg(7,8,i)/drho*drho
            fma(mc  ,ml+3)=fma(mc  ,ml+3) &
                          +coef(i)*table_hhg(8,8,i)/drho*drho*drho
         enddo
         enddo
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,2) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,2)=1.d0
         do mw=1,mwmax
            fma(mw,mlmax-1) = 0.d0
         enddo
         fma(mc,mlmax-1) = 1.d0

         fvb(1)=1.d0
         fvb(2)=-10.d0
         fvb(mlmax-1)=1.D0/11.D0

      return
      end subroutine fem_calc_l3

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_x1(nrmax,npow,nth,nph,rf,angl)

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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

               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho

               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho

! ----- dw/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,2,inod)*table_lll(inod,3,1)
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,2,inod)*table_lll(inod,4,1)

               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,2,inod)*table_lll(inod,3,2)
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,2,inod)*table_lll(inod,4,2)

! ----- dE/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,3,inod)*table_lll(inod,1,3)
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,3,inod)*table_lll(inod,2,3)

               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,3,inod)*table_lll(inod,1,4)
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,3,inod)*table_lll(inod,2,4)

! ----- dw/dx dE/dx terms

               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho

               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho
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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,2,inod)*table_lgg(inod,2,1)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,3,inod)*table_lgg(inod,1,2)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,2,inod)*table_lgl(inod,2,1)
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,2,inod)*table_lgl(inod,2,2)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,3,inod)*table_lgl(inod,1,3)
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,3,inod)*table_lgl(inod,1,4)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,2,inod)*table_llg(inod,3,1)
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,2,inod)*table_llg(inod,4,1)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,3,inod)*table_llg(inod,1,2)
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,3,inod)*table_llg(inod,2,2)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               enddo

            else

               do inod=1,2
! ----- non derivative terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho
! ----- dw/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,2,inod)*table_lll(inod,3,1)
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,2,inod)*table_lll(inod,4,1)
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,2,inod)*table_lll(inod,3,2)
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,2,inod)*table_lll(inod,4,2)
! ----- dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,3,inod)*table_lll(inod,1,3)
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,3,inod)*table_lll(inod,2,3)
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,3,inod)*table_lll(inod,1,4)
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,3,inod)*table_lll(inod,2,4)
! ----- dw/dx dE/dx terms
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
                             +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
                             +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
                             +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho
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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho &
                             +fmd(i,j,2,inod)*table_lgg(inod,2,1) &
                             +fmd(i,j,3,inod)*table_lgg(inod,1,2) &
                             +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               fma(mw+1,ml  )=fma(mw+1 ,ml  ) &
                             +fmd(i,j,1,inod)*table_lgg(inod,1,3)*drho &
                             +fmd(i,j,2,inod)*table_lgg(inod,2,3) &
                             +fmd(i,j,3,inod)*table_lgg(inod,1,4) &
                             +fmd(i,j,4,inod)*table_lgg(inod,2,4)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1) &
                             +fmd(i,j,1,inod)*table_lgg(inod,3,1)*drho &
                             +fmd(i,j,2,inod)*table_lgg(inod,4,1) &
                             +fmd(i,j,3,inod)*table_lgg(inod,3,2) &
                             +fmd(i,j,4,inod)*table_lgg(inod,4,2)/drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1) &
                             +fmd(i,j,1,inod)*table_lgg(inod,3,3)*drho &
                             +fmd(i,j,2,inod)*table_lgg(inod,4,3) &
                             +fmd(i,j,3,inod)*table_lgg(inod,3,4) &
                             +fmd(i,j,4,inod)*table_lgg(inod,4,4)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                            +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho &
                             +fmd(i,j,2,inod)*table_lgl(inod,2,1) &
                             +fmd(i,j,3,inod)*table_lgl(inod,1,3) &
                             +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1) &
                             +fmd(i,j,1,inod)*table_lgl(inod,3,1)*drho &
                             +fmd(i,j,2,inod)*table_lgl(inod,4,1) &
                             +fmd(i,j,3,inod)*table_lgl(inod,3,3) &
                             +fmd(i,j,4,inod)*table_lgl(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  ) &
                             +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho &
                             +fmd(i,j,2,inod)*table_lgl(inod,2,2) &
                             +fmd(i,j,3,inod)*table_lgl(inod,1,4) &
                             +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               fma(mw+3,ml+1)=fma(mw+3,ml+1) &
                             +fmd(i,j,1,inod)*table_lgl(inod,3,2)*drho &
                             +fmd(i,j,2,inod)*table_lgl(inod,4,2) &
                             +fmd(i,j,3,inod)*table_lgl(inod,3,4) &
                             +fmd(i,j,4,inod)*table_lgl(inod,4,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho &
                             +fmd(i,j,2,inod)*table_llg(inod,3,1) &
                             +fmd(i,j,3,inod)*table_llg(inod,1,2) &
                             +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw+1 ,ml  )=fma(mw+1,ml  ) &
                             +fmd(i,j,1,inod)*table_llg(inod,1,3)*drho &
                             +fmd(i,j,2,inod)*table_llg(inod,3,3) &
                             +fmd(i,j,3,inod)*table_llg(inod,1,4) &
                             +fmd(i,j,4,inod)*table_llg(inod,3,4)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4) &
                             +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho &
                             +fmd(i,j,2,inod)*table_llg(inod,4,1) &
                             +fmd(i,j,3,inod)*table_llg(inod,2,2) &
                             +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               fma(mw-3,ml+4)=fma(mw-3,ml+4) &
                             +fmd(i,j,1,inod)*table_llg(inod,2,3)*drho &
                             +fmd(i,j,2,inod)*table_llg(inod,4,3) &
                             +fmd(i,j,3,inod)*table_llg(inod,2,4) &
                             +fmd(i,j,4,inod)*table_llg(inod,4,4)/drho
               enddo

            else

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
                             +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho &
                             +fmd(i,j,2,inod)*table_lll(inod,3,1) &
                             +fmd(i,j,3,inod)*table_lll(inod,1,3) &
                             +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4) &
                             +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho &
                             +fmd(i,j,2,inod)*table_lll(inod,4,1) &
                             +fmd(i,j,3,inod)*table_lll(inod,2,3) &
                            +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  ) &
                             +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho &
                             +fmd(i,j,2,inod)*table_lll(inod,3,2) &
                             +fmd(i,j,3,inod)*table_lll(inod,1,4) &
                             +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+4)=fma(mw  ,ml+4) &
                             +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho &
                             +fmd(i,j,2,inod)*table_lll(inod,4,2) &
                             +fmd(i,j,3,inod)*table_lll(inod,2,4) &
                             +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho

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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
                  fmd(i,j,k,2)=0.5d0*(-11*fmc(i,j,k,1)+18*fmc(i,j,k,2) &
     &                                - 9*fmc(i,j,k,3)+ 2*fmc(i,j,k,4))
                  fmd(i,j,k,3)=fmc(i,j,k,4)
                  fmd(i,j,k,4)=0.5d0*(- 2*fmc(i,j,k,1)+ 9*fmc(i,j,k,2) &
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

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
                          +fmd(i,j,1,inod)*table_hhh(inod,1,1)*drho
            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
                          +fmd(i,j,1,inod)*table_hhh(inod,2,1)*drho**2
            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
                          +fmd(i,j,1,inod)*table_hhh(inod,3,1)*drho
            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
                          +fmd(i,j,1,inod)*table_hhh(inod,4,1)*drho**2

            fma(mw+1,ml  )=fma(mw+1,ml  ) &
                          +fmd(i,j,1,inod)*table_hhh(inod,1,2)*drho**2
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
                          +fmd(i,j,1,inod)*table_hhh(inod,2,2)*drho**3
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
                          +fmd(i,j,1,inod)*table_hhh(inod,3,2)*drho**2
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
                          +fmd(i,j,1,inod)*table_hhh(inod,4,2)*drho**3

            fma(mw+6,ml  )=fma(mw+6,ml  ) &
                          +fmd(i,j,1,inod)*table_hhh(inod,1,3)*drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
                          +fmd(i,j,1,inod)*table_hhh(inod,2,3)*drho**2
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
                          +fmd(i,j,1,inod)*table_hhh(inod,3,3)*drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
                          +fmd(i,j,1,inod)*table_hhh(inod,4,3)*drho**2

            fma(mw+7,ml  )=fma(mw+7,ml  ) &
                          +fmd(i,j,1,inod)*table_hhh(inod,1,4)*drho**2
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
                          +fmd(i,j,1,inod)*table_hhh(inod,2,4)*drho**3
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
                          +fmd(i,j,1,inod)*table_hhh(inod,3,4)*drho**2
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
                          +fmd(i,j,1,inod)*table_hhh(inod,4,4)*drho**3

! ----- dw/dx terms

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,1)
            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,1)*drho
            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,1)
            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,1)*drho

            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,2)*drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,2)*drho**2
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,2)*drho
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,2)*drho**2

            fma(mw+6,ml  )=fma(mw+6,ml  ) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,3)
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,3)*drho
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,3)
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,3)*drho

            fma(mw+7,ml  )=fma(mw+7,ml ) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,4)*drho
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,4)*drho**2
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,4)*drho
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,4)*drho**2

! ----- dE/dx terms

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,5)
            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,5)*drho
            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,5)
            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,5)*drho

            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,6)*drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,6)*drho**2
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,6)*drho
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,6)*drho**2

            fma(mw+6,ml  )=fma(mw+6,ml  ) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,7)
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,7)*drho
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,7)
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,7)*drho

            fma(mw+7,ml  )=fma(mw+7,ml  ) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,8)*drho
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,8)*drho**2
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,8)*drho
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,8)*drho**2

! ----- dw/dx dE/dx terms

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,5)/drho
            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,5)
            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,5)/drho
            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,5)

            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,6)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,6)*drho
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,6)
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,6)*drho

            fma(mw+6,ml  )=fma(mw+6,ml  ) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,7)/drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,7)
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,7)/drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,7)

            fma(mw+7,ml  )=fma(mw+7,ml  ) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,8)
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,8)*drho
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,8)
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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

               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4) &
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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,1)
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,2)
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,1) &
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,3) &
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,2) &
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,1) &
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,2) &
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho &
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,2) &
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               enddo

            else

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-3,ml+3)=fma(mw-3,ml+3) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+3,ml  )=fma(mw+3,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+3)=fma(mw  ,ml+3) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4) &
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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,1) &
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,2) &
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,2)/drho
               fma(mw+1,ml  )=fma(mw+1 ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgg(inod,1,3)*drho &
     &                       +fmd(i,j,2,inod)*table_lgg(inod,2,3) &
     &                       +fmd(i,j,3,inod)*table_lgg(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_lgg(inod,2,4)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                       +fmd(i,j,1,inod)*table_lgg(inod,3,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lgg(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_lgg(inod,3,2) &
     &                       +fmd(i,j,4,inod)*table_lgg(inod,4,2)/drho
               fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                       +fmd(i,j,1,inod)*table_lgg(inod,3,3)*drho &
     &                       +fmd(i,j,2,inod)*table_lgg(inod,4,3) &
     &                       +fmd(i,j,3,inod)*table_lgg(inod,3,4) &
     &                       +fmd(i,j,4,inod)*table_lgg(inod,4,4)/drho
               enddo

! ***** (1,*) *****
            else if(i.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,1) &
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,3) &
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,3)/drho
               fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                       +fmd(i,j,1,inod)*table_lgl(inod,3,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lgl(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_lgl(inod,3,3) &
     &                       +fmd(i,j,4,inod)*table_lgl(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lgl(inod,1,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lgl(inod,2,2) &
     &                       +fmd(i,j,3,inod)*table_lgl(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_lgl(inod,2,4)/drho
               fma(mw+3,ml+1)=fma(mw+3,ml+1) &
     &                       +fmd(i,j,1,inod)*table_lgl(inod,3,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lgl(inod,4,2) &
     &                       +fmd(i,j,3,inod)*table_lgl(inod,3,4) &
     &                       +fmd(i,j,4,inod)*table_lgl(inod,4,4)/drho
               enddo

! ***** (*,1) *****
            else if(j.eq.1) then

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,1) &
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,2) &
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,2)/drho
               fma(mw+1 ,ml  )=fma(mw+1,ml  ) &
     &                       +fmd(i,j,1,inod)*table_llg(inod,1,3)*drho &
     &                       +fmd(i,j,2,inod)*table_llg(inod,3,3) &
     &                       +fmd(i,j,3,inod)*table_llg(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_llg(inod,3,4)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4) &
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,1)*drho &
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,2) &
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,2)/drho
               fma(mw-3,ml+4)=fma(mw-3,ml+4) &
     &                       +fmd(i,j,1,inod)*table_llg(inod,2,3)*drho &
     &                       +fmd(i,j,2,inod)*table_llg(inod,4,3) &
     &                       +fmd(i,j,3,inod)*table_llg(inod,2,4) &
     &                       +fmd(i,j,4,inod)*table_llg(inod,4,4)/drho
               enddo

            else

               do inod=1,2
               fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,1) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,3) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
               fma(mw-4,ml+4)=fma(mw-4,ml+4) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,1) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,3) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
               fma(mw+4,ml  )=fma(mw+4,ml  ) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,3,2) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,1,4) &
     &                       +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho
               fma(mw  ,ml+4)=fma(mw  ,ml+4) &
     &                       +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho &
     &                       +fmd(i,j,2,inod)*table_lll(inod,4,2) &
     &                       +fmd(i,j,3,inod)*table_lll(inod,2,4) &
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
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
                  fmd(i,j,k,2)=0.5d0*(-11*fmc(i,j,k,1)+18*fmc(i,j,k,2) &
     &                                - 9*fmc(i,j,k,3)+ 2*fmc(i,j,k,4))
                  fmd(i,j,k,3)=fmc(i,j,k,4)
                  fmd(i,j,k,4)=0.5d0*(- 2*fmc(i,j,k,1)+ 9*fmc(i,j,k,2) &
     &                                -18*fmc(i,j,k,3)+11*fmc(i,j,k,4))
               enddo
            enddo
         enddo

         do j=1,3
         do i=1,3
            ml=6*(nr-1)+2*(i-1)+1
            mw=mc+2*(j-1)-2*(i-1)
            do inod=1,4

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,1) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,5) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,6)
            fma(mw+6,ml  )=fma(mw+6,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,3) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,7) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,7)/drho
            fma(mw+7,ml  )=fma(mw+7,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,1,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,5,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,1,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,6)*drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,7)
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,2,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,6,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,2,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,6,8)*drho

            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,1) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,5) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,5)/drho
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,3) &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,7) &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,3,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,7,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,3,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,7,8)

            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,5)
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hhh(inod,4,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhh(inod,8,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhh(inod,4,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhh(inod,8,8)*drho

               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,1)=-ci*nth
         fma(mc+2,1)=1.d0
         fma(mc,5)=1.d0
         fvb(1)=0.d0
         fvb(5)=0.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
         fvb(5)=0.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
                  rho0=rho(nr+1)/1.d-6
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
!                  fmd(i,j,k,2)=0.5d0*(-11*fmc(i,j,k,1)+18*fmc(i,j,k,2) &
!     &                                - 9*fmc(i,j,k,3)+ 2*fmc(i,j,k,4))
                  fmd(i,j,k,2)=1.5d0*(-3*fmc(i,j,k,1)+4*fmc(i,j,k,2) &
     &                                -  fmc(i,j,k,3))
                  fmd(i,j,k,3)=fmc(i,j,k,4)
!                  fmd(i,j,k,4)=0.5d0*(- 2*fmc(i,j,k,1)+ 9*fmc(i,j,k,2) &
!     &                                -18*fmc(i,j,k,3)+11*fmc(i,j,k,4))
                  fmd(i,j,k,4)=1.5d0*(   fmc(i,j,k,2) &
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

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,4,1) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,1,4) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,4,4)/drho

            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,1,2)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,4,2) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,1,5) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,4,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,1,3)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,4,3) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,1,6) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,4,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,2,1)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,5,1) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,2,4) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,5,4)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,2,2)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,5,2) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,2,5) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,5,5)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,2,3)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,5,3) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,2,6) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,5,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,3,1)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,6,1) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,3,4) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,6,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,3,2)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,6,2) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,3,5) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,6,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqq(inod,3,3)*drho &
     &              +fmd(i,j,2,inod)*table_hqq(inod,6,3) &
     &              +fmd(i,j,3,inod)*table_hqq(inod,3,6) &
     &              +fmd(i,j,4,inod)*table_hqq(inod,6,6)/drho

! (1,*) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,1) &
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,5) &
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,2)*drho &
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,6)*drho &
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,3)*drho &
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,3) &
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,7) &
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,1,4)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hqh(inod,4,4)*drho &
     &              +fmd(i,j,3,inod)*table_hqh(inod,1,8)*drho &
     &              +fmd(i,j,4,inod)*table_hqh(inod,4,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,1)*drho &
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,1) &
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,5) &
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,5)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,2)*drho &
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,6)*drho &
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,6)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,3)*drho &
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,3) &
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,7) &
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,7)/drho
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,2,4)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hqh(inod,5,4)*drho &
     &              +fmd(i,j,3,inod)*table_hqh(inod,2,8)*drho &
     &              +fmd(i,j,4,inod)*table_hqh(inod,5,8)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,1)*drho &
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,1) &
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,5) &
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,2)*drho &
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,6)*drho &
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,3)*drho &
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,3) &
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,7) &
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hqh(inod,3,4)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hqh(inod,6,4)*drho &
     &              +fmd(i,j,3,inod)*table_hqh(inod,3,8)*drho &
     &              +fmd(i,j,4,inod)*table_hqh(inod,6,8)

! (*,1) **********
            else if(j.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_hhq(inod,5,1) &
     &              +fmd(i,j,3,inod)*table_hhq(inod,1,4) &
     &              +fmd(i,j,4,inod)*table_hhq(inod,5,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,1,2)*drho &
     &              +fmd(i,j,2,inod)*table_hhq(inod,5,2) &
     &              +fmd(i,j,3,inod)*table_hhq(inod,1,5) &
     &              +fmd(i,j,4,inod)*table_hhq(inod,5,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,1,3)*drho &
     &              +fmd(i,j,2,inod)*table_hhq(inod,5,3) &
     &              +fmd(i,j,3,inod)*table_hhq(inod,1,6) &
     &              +fmd(i,j,4,inod)*table_hhq(inod,5,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,2,1)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhq(inod,6,1)*drho &
     &              +fmd(i,j,3,inod)*table_hhq(inod,2,4)*drho &
     &              +fmd(i,j,4,inod)*table_hhq(inod,6,4)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,2,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhq(inod,6,2)*drho &
     &              +fmd(i,j,3,inod)*table_hhq(inod,2,5)*drho &
     &              +fmd(i,j,4,inod)*table_hhq(inod,6,5)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,2,3)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhq(inod,6,3)*drho &
     &              +fmd(i,j,3,inod)*table_hhq(inod,2,6)*drho &
     &              +fmd(i,j,4,inod)*table_hhq(inod,6,6)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,3,1)*drho &
     &              +fmd(i,j,2,inod)*table_hhq(inod,7,1) &
     &              +fmd(i,j,3,inod)*table_hhq(inod,3,4) &
     &              +fmd(i,j,4,inod)*table_hhq(inod,7,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,3,2)*drho &
     &              +fmd(i,j,2,inod)*table_hhq(inod,7,2) &
     &              +fmd(i,j,3,inod)*table_hhq(inod,3,5) &
     &              +fmd(i,j,4,inod)*table_hhq(inod,7,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,3,3)*drho &
     &              +fmd(i,j,2,inod)*table_hhq(inod,7,3) &
     &              +fmd(i,j,3,inod)*table_hhq(inod,3,6) &
     &              +fmd(i,j,4,inod)*table_hhq(inod,7,6)/drho

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,4,1)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhq(inod,8,1)*drho &
     &              +fmd(i,j,3,inod)*table_hhq(inod,4,4)*drho &
     &              +fmd(i,j,4,inod)*table_hhq(inod,8,4)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,4,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhq(inod,8,2)*drho &
     &              +fmd(i,j,3,inod)*table_hhq(inod,4,5)*drho &
     &              +fmd(i,j,4,inod)*table_hhq(inod,8,5)
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1) &
     &              +fmd(i,j,1,inod)*table_hhq(inod,4,3)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhq(inod,8,3)*drho &
     &              +fmd(i,j,3,inod)*table_hhq(inod,4,6)*drho &
     &              +fmd(i,j,4,inod)*table_hhq(inod,8,6)


! (*,*) **********
            else
            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,1) &
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,5) &
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,2)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,6)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,3)*drho &
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,3) &
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,7) &
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,1,4)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,5,4)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,1,8)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,1)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,1)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,5)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,2)*drho**3 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,2)*drho**2 &
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,6)*drho**2 &
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,6)*drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,3)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,3)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,7)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,7)
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,2,4)*drho**3 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,6,4)*drho**2 &
     &              +fmd(i,j,3,inod)*table_hhh(inod,2,8)*drho**2 &
     &              +fmd(i,j,4,inod)*table_hhh(inod,6,8)*drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,1)*drho &
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,1) &
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,5) &
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,2)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,2)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,6)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,3)*drho &
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,3) &
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,7) &
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,3,4)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,7,4)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,3,8)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,7,8)

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,1)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,1)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,5)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,5)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,2)*drho**3 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,2)*drho**2 &
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,6)*drho**2 &
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1) &
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,3)*drho**2 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,3)*drho &
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,7)*drho &
     &              +fmd(i,j,4,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+mr+1)=fma(mw  ,ml+mr+1) & 
     &              +fmd(i,j,1,inod)*table_hhh(inod,4,4)*drho**3 &
     &              +fmd(i,j,2,inod)*table_hhh(inod,8,4)*drho**2 &
     &              +fmd(i,j,3,inod)*table_hhh(inod,4,8)*drho**2 &
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

      use libfem
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
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_lll(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_lll(inod,3,1) &
     &              +fmd(i,j,3,inod)*table_lll(inod,1,3) &
     &              +fmd(i,j,4,inod)*table_lll(inod,3,3)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_lll(inod,1,2)*drho &
     &              +fmd(i,j,2,inod)*table_lll(inod,3,2) &
     &              +fmd(i,j,3,inod)*table_lll(inod,1,4) &
     &              +fmd(i,j,4,inod)*table_lll(inod,3,4)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lll(inod,2,1)*drho &
     &              +fmd(i,j,2,inod)*table_lll(inod,4,1) &
     &              +fmd(i,j,3,inod)*table_lll(inod,2,3) &
     &              +fmd(i,j,4,inod)*table_lll(inod,4,3)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lll(inod,2,2)*drho &
     &              +fmd(i,j,2,inod)*table_lll(inod,4,2) &
     &              +fmd(i,j,3,inod)*table_lll(inod,2,4) &
     &              +fmd(i,j,4,inod)*table_lll(inod,4,4)/drho

! (1,*) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_llq(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_llq(inod,3,1) &
     &              +fmd(i,j,3,inod)*table_llq(inod,1,4) &
     &              +fmd(i,j,4,inod)*table_llq(inod,3,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_llq(inod,1,2)*drho &
     &              +fmd(i,j,2,inod)*table_llq(inod,3,2) &
     &              +fmd(i,j,3,inod)*table_llq(inod,1,5) &
     &              +fmd(i,j,4,inod)*table_llq(inod,3,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_llq(inod,1,3)*drho &
     &              +fmd(i,j,2,inod)*table_llq(inod,3,3) &
     &              +fmd(i,j,3,inod)*table_llq(inod,1,6) &
     &              +fmd(i,j,4,inod)*table_llq(inod,3,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_llq(inod,2,1)*drho &
     &              +fmd(i,j,2,inod)*table_llq(inod,4,1) &
     &              +fmd(i,j,3,inod)*table_llq(inod,2,4) &
     &              +fmd(i,j,4,inod)*table_llq(inod,4,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_llq(inod,2,2)*drho &
     &              +fmd(i,j,2,inod)*table_llq(inod,4,2) &
     &              +fmd(i,j,3,inod)*table_llq(inod,2,5) &
     &              +fmd(i,j,4,inod)*table_llq(inod,4,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_llq(inod,2,3)*drho &
     &              +fmd(i,j,2,inod)*table_llq(inod,4,3) &
     &              +fmd(i,j,3,inod)*table_llq(inod,2,6) &
     &              +fmd(i,j,4,inod)*table_llq(inod,4,6)/drho

! (*,1) **********
            else if(j.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_lql(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_lql(inod,4,1) &
     &              +fmd(i,j,3,inod)*table_lql(inod,1,3) &
     &              +fmd(i,j,4,inod)*table_lql(inod,4,3)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_lql(inod,1,2)*drho &
     &              +fmd(i,j,2,inod)*table_lql(inod,4,2) &
     &              +fmd(i,j,3,inod)*table_lql(inod,1,4) &
     &              +fmd(i,j,4,inod)*table_lql(inod,4,4)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_lql(inod,2,1)*drho &
     &              +fmd(i,j,2,inod)*table_lql(inod,5,1) &
     &              +fmd(i,j,3,inod)*table_lql(inod,2,3) &
     &              +fmd(i,j,4,inod)*table_lql(inod,5,3)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_lql(inod,2,2)*drho &
     &              +fmd(i,j,2,inod)*table_lql(inod,5,2) &
     &              +fmd(i,j,3,inod)*table_lql(inod,2,4) &
     &              +fmd(i,j,4,inod)*table_lql(inod,5,4)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lql(inod,3,1)*drho &
     &              +fmd(i,j,2,inod)*table_lql(inod,6,1) &
     &              +fmd(i,j,3,inod)*table_lql(inod,3,3) &
     &              +fmd(i,j,4,inod)*table_lql(inod,6,3)/drho
            fma(mw,ml+mr)=fma(mw,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lql(inod,3,2)*drho &
     &              +fmd(i,j,2,inod)*table_lql(inod,6,2) &
     &              +fmd(i,j,3,inod)*table_lql(inod,3,4) &
     &              +fmd(i,j,4,inod)*table_lql(inod,6,4)/drho

! (*,*) **********
            else
            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,1,1)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,4,1) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,1,4) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,4,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,1,2)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,4,2) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,1,5) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,4,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  ) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,1,3)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,4,3) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,1,6) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,4,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,2,1)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,5,1) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,2,4) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,5,4)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,2,2)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,5,2) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,2,5) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,5,5)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,2,3)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,5,3) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,2,6) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,5,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,3,1)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,6,1) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,3,4) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,6,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,3,2)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,6,2) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,3,5) &
     &              +fmd(i,j,4,inod)*table_lqq(inod,6,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr) &
     &              +fmd(i,j,1,inod)*table_lqq(inod,3,3)*drho &
     &              +fmd(i,j,2,inod)*table_lqq(inod,6,3) &
     &              +fmd(i,j,3,inod)*table_lqq(inod,3,6) &
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

!$$$      do ml=1,mlmax
!$$$         do mw=1,mwmax
!$$$            if(abs(fma(mw,ml)).ne.0.d0) then
!$$$               write(6,'(2I5,2ES12.4)') ml,mw,fma(mw,ml)
!$$$            endif
!$$$         enddo
!$$$      enddo
      return
      end subroutine fem_calc_r6

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r7(nrmax,npow,nth,nph,rf,angl)

      use libfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,mr,mu1,mu2,ic2,ic1,ip1,ip2
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      complex(8),dimension(3,3,4):: fmd1,fmd2
      complex(8),dimension(8,8):: fml
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer,dimension(8):: ica =(/2,3,1,2,3,1,2,3/) ! column number
      integer,dimension(8):: ipa =(/1,1,1,3,3,3,5,5/) ! position number

      call mesh_init(nrmax,npow)

      nvmax=8
      mlmax=6*nrmax-4
      mwmax=2*(8+6)-1
      
      call fem_init


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
      
         do k=1,4
            do j=1,3
               do i=1,3
                  fmd1(i,j,k)=0.d0
                  fmd2(i,j,k)=0.d0
               enddo
            enddo
         enddo

         rho0=rho(nr)+0.25D0*drho
         rkth=nth/rho0
         rkth0=1.d0/rho0

         fmd1(1,1,1)= rho0*(factor-rkph**2-rkth**2)
         fmd1(1,2,1)= rho0*(-ci*rkth*rkth0)
         fmd1(2,1,1)= rho0*(+ci*rkth*rkth0)
         fmd1(2,2,1)= rho0*(factor-rkph**2-rkth0**2)
         fmd1(2,3,1)= rho0*(rkth*rkph)
         fmd1(3,2,1)= rho0*(rkth*rkph)
         fmd1(3,3,1)= rho0*(factor-rkth**2)
         
         fmd1(2,1,2)= rho0*( ci*rkth)
         fmd1(2,2,2)= rho0*(  -rkth0)
         fmd1(3,1,2)= rho0*( ci*rkph)
         
         fmd1(1,2,3)= rho0*(-ci*rkth)
         fmd1(1,3,3)= rho0*(-ci*rkph)
         fmd1(2,2,3)= rho0*(  -rkth0)

         fmd1(2,2,4)= rho0*(-1.d0)
         fmd1(3,3,4)= rho0*(-1.d0)

         rho0=rho(nr)+0.75D0*drho
         rkth=nth/rho0
         rkth0=1.d0/rho0

         fmd2(1,1,1)= rho0*(factor-rkph**2-rkth**2)
         fmd2(1,2,1)= rho0*(-ci*rkth*rkth0)
         fmd2(2,1,1)= rho0*(+ci*rkth*rkth0)
         fmd2(2,2,1)= rho0*(factor-rkph**2-rkth0**2)
         fmd2(2,3,1)= rho0*(rkth*rkph)
         fmd2(3,2,1)= rho0*(rkth*rkph)
         fmd2(3,3,1)= rho0*(factor-rkth**2)
      
         fmd2(2,1,2)= rho0*( ci*rkth)
         fmd2(2,2,2)= rho0*(  -rkth0)
         fmd2(3,1,2)= rho0*( ci*rkph)

         fmd2(1,2,3)= rho0*(-ci*rkth)
         fmd2(1,3,3)= rho0*(-ci*rkph)
         fmd2(2,2,3)= rho0*(  -rkth0)

         fmd2(2,2,4)= rho0*(-1.d0)
         fmd2(3,3,4)= rho0*(-1.d0)

         do mu2=1,8
            do mu1=1,8
               fml(mu1,mu2)=0.d0
            end do
         end do

         do mu1=1,8
            ic1=ica(mu1)
            ip1=ipa(mu1)
         do mu2=1,8
            ic2=ica(mu2)
            ip2=ipa(mu2)
            
            if((ic1 == 1) .AND. (ic2 == 1)) THEN
               fml(mu1,mu2)=fml(mu1,mu2) &
              +fmd1(ic1,ic2,1)*table_ppp(1,ip1,  ip2  )*drho &
              +fmd1(ic1,ic2,2)*table_ppp(1,ip1+1,ip2  ) &
              +fmd1(ic1,ic2,3)*table_ppp(1,ip1,  ip2+1) &
              +fmd1(ic1,ic2,4)*table_ppp(1,ip1+1,ip2+1)/drho &
              +fmd2(ic1,ic2,1)*table_ppp(2,ip1,  ip2  )*drho &
              +fmd2(ic1,ic2,2)*table_ppp(2,ip1+1,ip2  ) &
              +fmd2(ic1,ic2,3)*table_ppp(2,ip1,  ip2+1) &
              +fmd2(ic1,ic2,4)*table_ppp(2,ip1+1,ip2+1)/drho

            else if(ic1 == 1) THEN
               fml(mu1,mu2)=fml(mu1,mu2) &
              +fmd1(ic1,ic2,1)*table_ppq(1,ip1,  ip2  )*drho &
              +fmd1(ic1,ic2,2)*table_ppq(1,ip1+1,ip2  ) &
              +fmd1(ic1,ic2,3)*table_ppq(1,ip1,  ip2+1) &
              +fmd1(ic1,ic2,4)*table_ppq(1,ip1+1,ip2+1)/drho &
              +fmd2(ic1,ic2,1)*table_ppq(2,ip1,  ip2  )*drho &
              +fmd2(ic1,ic2,2)*table_ppq(2,ip1+1,ip2  ) &
              +fmd2(ic1,ic2,3)*table_ppq(2,ip1,  ip2+1) &
              +fmd2(ic1,ic2,4)*table_ppq(2,ip1+1,ip2+1)/drho

            else if(ic2 == 1) THEN
               fml(mu1,mu2)=fml(mu1,mu2) &
              +fmd1(ic1,ic2,1)*table_pqp(1,ip1,  ip2  )*drho &
              +fmd1(ic1,ic2,2)*table_pqp(1,ip1+1,ip2  ) &
              +fmd1(ic1,ic2,3)*table_pqp(1,ip1,  ip2+1) &
              +fmd1(ic1,ic2,4)*table_pqp(1,ip1+1,ip2+1)/drho &
              +fmd2(ic1,ic2,1)*table_pqp(2,ip1,  ip2  )*drho &
              +fmd2(ic1,ic2,2)*table_pqp(2,ip1+1,ip2  ) &
              +fmd2(ic1,ic2,3)*table_pqp(2,ip1,  ip2+1) &
              +fmd2(ic1,ic2,4)*table_pqp(2,ip1+1,ip2+1)/drho

            else
               fml(mu1,mu2)=fml(mu1,mu2) &
              +fmd1(ic1,ic2,1)*table_pqq(1,ip1,  ip2  )*drho &
              +fmd1(ic1,ic2,2)*table_pqq(1,ip1+1,ip2  ) &
              +fmd1(ic1,ic2,3)*table_pqq(1,ip1,  ip2+1) &
              +fmd1(ic1,ic2,4)*table_pqq(1,ip1+1,ip2+1)/drho &
              +fmd2(ic1,ic2,1)*table_pqq(2,ip1,  ip2  )*drho &
              +fmd2(ic1,ic2,2)*table_pqq(2,ip1+1,ip2  ) &
              +fmd2(ic1,ic2,3)*table_pqq(2,ip1,  ip2+1) &
              +fmd2(ic1,ic2,4)*table_pqq(2,ip1+1,ip2+1)/drho
            end if
         END DO
         END DO

         ml=6*(nr-1)
         do mu2=1,8
            do mu1=1,8
               fma(mc+mu2-mu1,ml+mu1)=fma(mc+mu2-mu1,ml+mu1) &
                    +fml(mu1,mu2)
            end do
         end do
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
         enddo
         fma(mc,1)=1.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,2) = 0.d0
         enddo
         fma(mc,1)=ci*nth
         fma(mc+2,1)= 4.D0/3.D0
         fma(mc+5,1)=-1.D0/3.D0
         fma(mc,2)=1.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,2) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,2)=1.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-1) = 0.d0
         fma(mw,mlmax  ) = 0.d0
      enddo
      fma(mc,mlmax-1) = 1.d0
      fma(mc,mlmax  ) = 1.d0

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
            fvb(6*(nr-1)+4)=(1.d0-angl)
            fvb(6*(nr-1)+5)=angl
         endif
      enddo

!$$$      do ml=1,mlmax
!$$$         do mw=1,mwmax
!$$$            if(abs(fma(mw,ml)).ne.0.d0) then
!$$$               write(6,'(2I5,2ES12.4)') ml,mw,fma(mw,ml)
!$$$            endif
!$$$         enddo
!$$$      enddo
      return
      end subroutine fem_calc_r7

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r8(nrmax,npow,nth,nph,rf,angl)

      use libfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor,rkth0,rkthp,rkthm,rho0
      complex(8),dimension(3,3,4,4):: fmc,fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=6
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
      
         do inod=1,1
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,1
            rho0=0.5D0*(rho(nr)+rho(nr+1))
            rkth0= nth   /rho0
            rkthp=(nth+1)/rho0
            rkthm=(nth-1)/rho0

            fmd(1,1,1,inod)= rho0*( 0.25D0*rkthp*rkthp+0.5D0*rkph**2-factor)
            fmd(1,2,1,inod)= rho0*( 0.25D0*rkthp*rkthm)
            fmd(1,3,1,inod)= rho0*(-0.50D0*CI*rkth0*rkph)
            fmd(2,1,1,inod)= rho0*( 0.25D0*rkthm*rkthp)
            fmd(2,2,1,inod)= rho0*( 0.25D0*rkthm*rkthm+0.5D0*rkph**2-factor)
            fmd(2,3,1,inod)= rho0*( 0.50D0*CI*rkth0*rkph)
            fmd(3,1,1,inod)= rho0*( 0.50D0*CI*rkth0*rkph)
            fmd(3,2,1,inod)= rho0*(-0.50D0*CI*rkth0*rkph)
            fmd(3,3,1,inod)= rho0*(        rkth0*rkth0              -factor)
      
            fmd(1,1,2,inod)= rho0*( 0.25D0*rkthp)
            fmd(1,2,2,inod)= rho0*( 0.25D0*rkthm)
            fmd(2,1,2,inod)= rho0*(-0.25D0*rkthp)
            fmd(2,2,2,inod)= rho0*(-0.25D0*rkthm)
            fmd(3,1,2,inod)= rho0*(-0.50D0*CI*rkph)
            fmd(3,2,2,inod)= rho0*(-0.50D0*CI*rkph)

            fmd(1,1,3,inod)= rho0*( 0.25D0*rkthp)
            fmd(1,2,3,inod)= rho0*(-0.25D0*rkthp)
            fmd(1,3,3,inod)= rho0*( 0.50D0*CI*rkph)
            fmd(2,1,3,inod)= rho0*( 0.25D0*rkthm)
            fmd(2,2,3,inod)= rho0*(-0.25D0*rkthm)
            fmd(2,3,3,inod)= rho0*( 0.50D0*CI*rkph)

            fmd(1,1,4,inod)= rho0*0.25D0
            fmd(1,2,4,inod)=-rho0*0.25D0
            fmd(2,1,4,inod)=-rho0*0.25D0
            fmd(2,2,4,inod)= rho0*0.25D0
            fmd(3,3,4,inod)= rho0
         enddo

         do j=1,3
         do i=1,3
            ml=6*(nr-1)+2*(i-1)+1
            mw=mc+2*(j-1)-2*(i-1)
            do inod=1,1

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(5,1) &
     &                    +fmd(i,j,3,inod)*table_hh(1,5) &
     &                    +fmd(i,j,4,inod)*table_hh(5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(5,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(1,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(5,6)
            fma(mw+6,ml  )=fma(mw+6,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(5,3) &
     &                    +fmd(i,j,3,inod)*table_hh(1,7) &
     &                    +fmd(i,j,4,inod)*table_hh(5,7)/drho
            fma(mw+7,ml  )=fma(mw+7,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(5,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(1,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(6,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(2,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(6,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(2,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(6,6)*drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(6,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(2,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(6,7)
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(6,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(2,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(6,8)*drho

            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(7,1) &
     &                    +fmd(i,j,3,inod)*table_hh(3,5) &
     &                    +fmd(i,j,4,inod)*table_hh(7,5)/drho
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(7,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(3,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(7,6)
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(7,3) &
     &                    +fmd(i,j,3,inod)*table_hh(3,7) &
     &                    +fmd(i,j,4,inod)*table_hh(7,7)/drho
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(7,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(3,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(7,8)

            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(8,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(4,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(8,5)
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(8,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(4,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(8,6)*drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(8,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(4,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(8,7)
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(8,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(4,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(8,8)*drho

               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
      elseif(nth.eq.1) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,5)=1.d0
         fvb(1)=0.d0
         fvb(5)=0.d0
      elseif(nth.eq.-1) then
         do mw=1,mwmax
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
         fvb(3)=0.d0
         fvb(5)=0.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
         fvb(5)=0.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-5) = 0.d0
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
      fma(mc,mlmax-5) = 1.d0
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0
         fvb(mlmax-5)=0.d0
         fvb(mlmax-3)=0.d0
         fvb(mlmax-1)=0.d0

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-0.85d0,0.85d0-rho(nr)
            fvb(6*(nr-1)+1)= CI*(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(6* nr   +1)= CI*(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(6*(nr-1)+3)=-CI*(rho(nr+1)-0.85d0)*(1.d0-angl)
            fvb(6* nr   +3)=-CI*(0.85d0-rho(nr)  )*(1.d0-angl)
            fvb(6*(nr-1)+5)=(rho(nr+1)-0.85d0)*angl
            fvb(6* nr   +5)=(0.85d0-rho(nr)  )*angl
         endif
      enddo

      return
      end subroutine fem_calc_r8

!----- calculate coefficint matrix fma and source vector fvb -----

      subroutine fem_calc_r9(nrmax,npow,nth,nph,rf,angl)

      use libfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod
      real(8):: drho,rkth,rkph,factor,rkth0,rho0,rd,divj
      complex(8),dimension(3,3,4,4):: fmd
      complex(8),parameter:: ci=(0.d0,1.d0)

      call mesh_init(nrmax,npow)

      nvmax=6
      mwmax=4*nvmax-1         ! width of coefficient matrix
      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector
     
      call fem_init

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
      
         do inod=1,1
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmd(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,1
            rho0=0.5D0*(rho(nr)+rho(nr+1))
            rkth=nth/rho0

            fmd(1,1,1,inod)= rho0*(-rkph**2-rkth**2+factor-1.D0/(rho0**2))
            fmd(1,2,1,inod)= rho0*(-2.D0*ci*rkth/rho0)
            fmd(2,1,1,inod)= rho0*( 2.D0*ci*rkth/rho0)
            fmd(2,2,1,inod)= rho0*(-rkph**2-rkth**2+factor-1.D0/(rho0**2))
            fmd(3,3,1,inod)= rho0*(-rkph**2-rkth**2+factor)
      
            fmd(1,1,4,inod)= rho0*(-1.D0)
            fmd(2,2,4,inod)= rho0*(-1.D0)
            fmd(3,3,4,inod)= rho0*(-1.D0)
         enddo

         do j=1,3
         do i=1,3
            ml=6*(nr-1)+2*(i-1)+1
            mw=mc+2*(j-1)-2*(i-1)
            do inod=1,1

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(5,1) &
     &                    +fmd(i,j,3,inod)*table_hh(1,5) &
     &                    +fmd(i,j,4,inod)*table_hh(5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(5,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(1,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(5,6)
            fma(mw+6,ml  )=fma(mw+6,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(5,3) &
     &                    +fmd(i,j,3,inod)*table_hh(1,7) &
     &                    +fmd(i,j,4,inod)*table_hh(5,7)/drho
            fma(mw+7,ml  )=fma(mw+7,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(5,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(1,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(6,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(2,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(6,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(2,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(6,6)*drho
            fma(mw+5,ml+1)=fma(mw+5,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(6,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(2,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(6,7)
            fma(mw+6,ml+1)=fma(mw+6,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(6,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(2,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(6,8)*drho

            fma(mw-6,ml+6)=fma(mw-6,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(7,1) &
     &                    +fmd(i,j,3,inod)*table_hh(3,5) &
     &                    +fmd(i,j,4,inod)*table_hh(7,5)/drho
            fma(mw-5,ml+6)=fma(mw-5,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(7,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(3,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(7,6)
            fma(mw  ,ml+6)=fma(mw  ,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(7,3) &
     &                    +fmd(i,j,3,inod)*table_hh(3,7) &
     &                    +fmd(i,j,4,inod)*table_hh(7,7)/drho
            fma(mw+1,ml+6)=fma(mw+1,ml+6) &
     &                    +fmd(i,j,1,inod)*table_hh(3,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(7,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(3,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(7,8)

            fma(mw-7,ml+7)=fma(mw-7,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(8,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(4,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(8,5)
            fma(mw-6,ml+7)=fma(mw-6,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(8,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(4,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(8,6)*drho
            fma(mw-1,ml+7)=fma(mw-1,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(8,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(4,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(8,7)
            fma(mw  ,ml+7)=fma(mw  ,ml+7) &
     &                    +fmd(i,j,1,inod)*table_hh(4,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(8,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(4,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(8,8)*drho
               enddo
            enddo
         enddo
      enddo

      if(nth.eq.0) then
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
      elseif(abs(nth).eq.1) then
         do mw=1,mwmax
            fma(mw,1) = 0.D0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,1)=-ci*nth
         fma(mc+2,1)=1.D0
         fma(mc,5)=1.d0
         fvb(1)=0.d0
         fvb(5)=0.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
         fvb(5)=0.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0
      fvb(mlmax-3) = 0.d0
      fvb(mlmax-1) = 0.d0

      rd=0.85D0
      rkth=nth/rd
      rkph=nph
      do nr=1,nrmax-1
         if((rd-rho(nr))*(rho(nr+1)-rd).ge.0.d0) then
!$$$            write(6,'A,I5,4ES12.4') 'Antenna:',nr,rho(nr),rho(nr+1),
!$$$     &           rho(nr+1)-rd,rd-rho(nr)
            divj=rkth*(1.D0-angl)+rkph*angl
            drho=rho(nr+1)-rho(nr)
            fvb(6*(nr-1)+2)=(rho(nr+1)-rd)*(-ci)*divj
            fvb(6* nr   +2)=(rd-rho(nr)  )*(-ci)*divj
            fvb(6*(nr-1)+3)=(rho(nr+1)-rd)*((1.d0-angl)+rkth*divj)
            fvb(6* nr   +3)=(rd-rho(nr)  )*((1.d0-angl)+rkth*divj)
            fvb(6*(nr-1)+5)=(rho(nr+1)-rd)*(angl+rkph*divj)
            fvb(6* nr   +5)=(rd-rho(nr)  )*(angl+rkph*divj)
         endif
      enddo

      return
      end subroutine fem_calc_r9

!----- calculate coefficint matrix fma and source vector fvb -----


END MODULE fem_calc
