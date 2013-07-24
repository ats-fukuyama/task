      subroutine fem_calc_16(nrmax,npow,nth,nph,rf,angl)

      use libfem
      use libtestfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(8),intent(in):: rf     ! wave frequency
      real(8),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nrd
      real(8):: drho,rkth,rkph,factor,rkth0,rho0,rd
      complex(8),dimension(4,4,4,4):: fmd
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer:: ig,l
      complex(8):: moment0,moment1,moment2,moment3,divj
      real(8):: xi,yi,x
      real(8),parameter:: m0=4.D0
      real(8),parameter:: m2=5.D0/16.D0
      real(8),parameter:: m4=41.D0/1024.D0
      real(8),parameter:: m6=365.D0/65536.D0

      call mesh_init(nrmax,npow)

      nvmax=8                ! vector/scalor potentials and their derivatives
      mwmax=4*nvmax-1        ! width of coefficient matrix
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
               do j=1,4
                  do i=1,4
                     fmd(i,j,k,inod)=0.d0
                  enddo
               enddo
            enddo
         enddo

         do inod=1,4
            SELECT CASE(inod)
            CASE(1)
               rho0=(7.D0*rho(nr)+     rho(nr+1))/8.D0
            CASE(2)
               rho0=(5.D0*rho(nr)+3.D0*rho(nr+1))/8.D0
            CASE(3)
               rho0=(3.D0*rho(nr)+5.D0*rho(nr+1))/8.D0
            CASE(4)
               rho0=(     rho(nr)+7.D0*rho(nr+1))/8.D0
            END SELECT

            rkth=nth/rho0

            fmd(1,1,1,inod)= rho0*(-rkph**2-rkth**2+factor-1.D0/(rho0**2))
            fmd(1,2,1,inod)= rho0*(-2.D0*ci*rkth/rho0)
            fmd(1,4,3,inod)= rho0*(ci)
            fmd(2,1,1,inod)= rho0*( 2.D0*ci*rkth/rho0)
            fmd(2,2,1,inod)= rho0*(-rkph**2-rkth**2+factor-1.D0/(rho0**2))
            fmd(2,4,1,inod)= rho0*(-rkth)
            fmd(3,3,1,inod)= rho0*(-rkph**2-rkth**2+factor)
            fmd(3,4,1,inod)= rho0*(-rkph)
      
            fmd(1,1,4,inod)= rho0*(-1.D0)
            fmd(2,2,4,inod)= rho0*(-1.D0)
            fmd(3,3,4,inod)= rho0*(-1.D0)

            fmd(4,1,2,inod)= rho0*(ci)
            fmd(4,2,1,inod)= rho0*( rkth)
            fmd(4,3,1,inod)= rho0*( rkph)
            fmd(4,4,1,inod)= rho0*(-rkph**2-rkth**2)
            fmd(4,4,4,inod)= rho0*(-1.D0)
         enddo

         DO k=1,4
            DO j=1,4
               DO i=1,4
                  moment0=0.D0
                  moment1=0.D0
                  moment2=0.D0
                  moment3=0.D0
                  DO inod=1,4
                     xi=0.125D0*(2*inod-1)
                     yi=xi-0.5D0
                     moment0=moment0+      fmd(i,j,k,inod)
                     moment1=moment1+yi   *fmd(i,j,k,inod)
                     moment2=moment2+yi**2*fmd(i,j,k,inod)
                     moment3=moment3+yi**3*fmd(i,j,k,inod)
                  END DO
                  fmd(i,j,k,1)=(m4*moment0-m2*moment2)/(m0*m4-m2**2)
                  fmd(i,j,k,2)=(m6*moment1-m4*moment3)/(m2*m6-m4**2)
                  fmd(i,j,k,3)=(m2*moment0-m0*moment2)/(m2**2-m0*m4)*2.D0
                  fmd(i,j,k,4)=(m4*moment1-m2*moment3)/(m4**2-m2*m6)*6.D0
               END DO
            END DO
         END DO

         do j=1,4
         do i=1,4
            ml=8*(nr-1)+2*(i-1)+1
            mw=mc+2*(j-1)-2*(i-1)
            do inod=1,4
               ig=2*inod-1

            fma(mw  ,ml  )=fma(mw  ,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhg(1,1,ig)*drho &
     &                    +fmd(i,j,2,inod)*table_hhg(5,1,ig) &
     &                    +fmd(i,j,3,inod)*table_hhg(1,5,ig) &
     &                    +fmd(i,j,4,inod)*table_hhg(5,5,ig)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhg(1,2,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(5,2,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(1,6,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(5,6,ig)
            fma(mw+8,ml  )=fma(mw+8,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhg(1,3,ig)*drho &
     &                    +fmd(i,j,2,inod)*table_hhg(5,3,ig) &
     &                    +fmd(i,j,3,inod)*table_hhg(1,7,ig) &
     &                    +fmd(i,j,4,inod)*table_hhg(5,7,ig)/drho
            fma(mw+9,ml  )=fma(mw+9,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hhg(1,4,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(5,4,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(1,8,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(5,8,ig)

            fma(mw-1,ml+1)=fma(mw-1,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhg(2,1,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(6,1,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(2,5,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(6,5,ig)
            fma(mw  ,ml+1)=fma(mw  ,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhg(2,2,ig)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhg(6,2,ig)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhg(2,6,ig)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhg(6,6,ig)*drho
            fma(mw+7,ml+1)=fma(mw+7,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhg(2,3,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(6,3,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(2,7,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(6,7,ig)
            fma(mw+8,ml+1)=fma(mw+8,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hhg(2,4,ig)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhg(6,4,ig)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhg(2,8,ig)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhg(6,8,ig)*drho

            fma(mw-8,ml+8)=fma(mw-8,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hhg(3,1,ig)*drho &
     &                    +fmd(i,j,2,inod)*table_hhg(7,1,ig) &
     &                    +fmd(i,j,3,inod)*table_hhg(3,5,ig) &
     &                    +fmd(i,j,4,inod)*table_hhg(7,5,ig)/drho
            fma(mw-7,ml+8)=fma(mw-7,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hhg(3,2,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(7,2,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(3,6,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(7,6,ig)
            fma(mw  ,ml+8)=fma(mw  ,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hhg(3,3,ig)*drho &
     &                    +fmd(i,j,2,inod)*table_hhg(7,3,ig) &
     &                    +fmd(i,j,3,inod)*table_hhg(3,7,ig) &
     &                    +fmd(i,j,4,inod)*table_hhg(7,7,ig)/drho
            fma(mw+1,ml+8)=fma(mw+1,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hhg(3,4,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(7,4,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(3,8,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(7,8,ig)

            fma(mw-9,ml+9)=fma(mw-9,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hhg(4,1,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(8,1,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(4,5,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(8,5,ig)
            fma(mw-8,ml+9)=fma(mw-8,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hhg(4,2,ig)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhg(8,2,ig)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhg(4,6,ig)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhg(8,6,ig)*drho
            fma(mw-1,ml+9)=fma(mw-1,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hhg(4,3,ig)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hhg(8,3,ig)*drho &
     &                    +fmd(i,j,3,inod)*table_hhg(4,7,ig)*drho &
     &                    +fmd(i,j,4,inod)*table_hhg(8,7,ig)
            fma(mw  ,ml+9)=fma(mw  ,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hhg(4,4,ig)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hhg(8,4,ig)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hhg(4,8,ig)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hhg(8,8,ig)*drho
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
         do mw=1,mwmax-2
            fma(mw,3) = fma(mw,3)+ci*nth*fma(mw+2,1)
         end do
         fvb(3)=fvb(3)+ci*nth*fvb(1)
         do mw=1,mwmax
            fma(mw,1) = 0.D0
            fma(mw,5) = 0.d0
            fma(mw,8) = 0.d0
         enddo
         fma(mc,1)=-ci*nth
         fma(mc+2,1)=1.D0
         fma(mc+6,1)=0.D0
         fma(mc+7,1)=0.D0
         fma(mc,5)=1.d0
         fma(mc,8)=1.d0
         fvb(1)=0.d0
         fvb(5)=0.d0
         fvb(8)=0.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
            fma(mw,8) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
         fma(mc,8)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
         fvb(5)=0.d0
         fvb(8)=0.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-7) = 0.d0
         fma(mw,mlmax-5) = 0.d0
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
!      fma(mc,mlmax-6) = 1.d0
      fma(mc,mlmax-7) = 1.d0
      fma(mc+1,mlmax-7) = 1.d0
      fma(mc,mlmax-5) = 1.d0
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0
!      fvb(mlmax-6) = 0.d0
      fvb(mlmax-7) = 0.d0
      fvb(mlmax-5) = 0.d0
      fvb(mlmax-3) = 0.d0
      fvb(mlmax-1) = 0.d0

      rd=0.85D0
      rkth=nth/rd
      rkph=nph
      nrd=0
      do nr=1,nrmax-1
         if((rd-rho(nr))*(rho(nr+1)-rd).ge.0.d0) nrd=nr
      end do

      nr=nrd
        x=(rd-rho(nr))/(rho(nr+1)-rho(nr))
        divj=-ci*(nth*(1.D0-angl)/rd +nph*angl)*(rho(nr+1)-rho(nr))

        fvb(8*(nr-1)+1)=divj*(fem_func_h(1.D0,1,2)-fem_func_h(x,1,2))
        fvb(8* nr   +1)=divj*(fem_func_h(1.D0,3,2)-fem_func_h(x,3,2))
        fvb(8*(nr-1)+3)=(1.d0-angl)*fem_func_h(x,1,0)
        fvb(8* nr   +3)=(1.d0-angl)*fem_func_h(x,3,0)
        fvb(8*(nr-1)+5)=angl*fem_func_h(x,1,0)
        fvb(8* nr   +5)=angl*fem_func_h(x,3,0)
      do nr=nrd+1,nrmax-1
!         fvb(8*(nr-1)+1)=fvb(8*(nr-1)+1)+divj*fem_func_h(1.D0,1,2)
!         fvb(8* nr   +1)=fvb(8* nr   +1)+divj*fem_func_h(1.D0,3,2)
      end do

      return
      end subroutine fem_calc_16

