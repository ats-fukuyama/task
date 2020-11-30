      subroutine fem_calc_15(nrmax,npow,nth,nph,rf,angl)

      USE task_kinds,ONLY: dp
      use libfem
      use libtestfem
      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh points
      integer,intent(in):: nth    ! poloidal mode number
      integer,intent(in):: nph    ! toroidal mode number
      real(dp),intent(in):: rf     ! wave frequency
      real(dp),intent(in):: angl   ! antenna angle: 0 perm, 1,para
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nrd
      real(dp):: drho,rkth,rkph,factor,omega,rkth0,rho0,rd,divj
      complex(dp),dimension(4,4,4,1):: fmd
      complex(dp),parameter:: ci=(0.d0,1.d0)

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
         omega=rf
         rkph=nph
      
         do inod=1,1
            do k=1,4
               do j=1,4
                  do i=1,4
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
            fmd(1,4,3,inod)= rho0*(ci)*omega
            fmd(2,1,1,inod)= rho0*( 2.D0*ci*rkth/rho0)
            fmd(2,2,1,inod)= rho0*(-rkph**2-rkth**2+factor-1.D0/(rho0**2))
            fmd(2,4,1,inod)= rho0*(-rkth)*omega
            fmd(3,3,1,inod)= rho0*(-rkph**2-rkth**2+factor)
            fmd(3,4,1,inod)= rho0*(-rkph)*omega
      
            fmd(1,1,4,inod)= rho0*(-1.D0)
            fmd(2,2,4,inod)= rho0*(-1.D0)
            fmd(3,3,4,inod)= rho0*(-1.D0)

            fmd(4,1,2,inod)= rho0*(ci)*omega
            fmd(4,2,1,inod)= rho0*( rkth)*omega
            fmd(4,3,1,inod)= rho0*( rkph)*omega

            fmd(4,4,1,inod)= rho0*(-rkph**2-rkth**2)
            fmd(4,4,4,inod)= rho0*(-1.D0)
         enddo

         do j=1,4
         do i=1,4
            ml=8*(nr-1)+2*(i-1)+1
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
            fma(mw+8,ml  )=fma(mw+8,ml  ) &
     &                    +fmd(i,j,1,inod)*table_hh(1,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(5,3) &
     &                    +fmd(i,j,3,inod)*table_hh(1,7) &
     &                    +fmd(i,j,4,inod)*table_hh(5,7)/drho
            fma(mw+9,ml  )=fma(mw+9,ml  ) &
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
            fma(mw+7,ml+1)=fma(mw+7,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(6,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(2,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(6,7)
            fma(mw+8,ml+1)=fma(mw+8,ml+1) &
     &                    +fmd(i,j,1,inod)*table_hh(2,4)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(6,4)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(2,8)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(6,8)*drho

            fma(mw-8,ml+8)=fma(mw-8,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hh(3,1)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(7,1) &
     &                    +fmd(i,j,3,inod)*table_hh(3,5) &
     &                    +fmd(i,j,4,inod)*table_hh(7,5)/drho
            fma(mw-7,ml+8)=fma(mw-7,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hh(3,2)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(7,2)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(3,6)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(7,6)
            fma(mw  ,ml+8)=fma(mw  ,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hh(3,3)*drho &
     &                    +fmd(i,j,2,inod)*table_hh(7,3) &
     &                    +fmd(i,j,3,inod)*table_hh(3,7) &
     &                    +fmd(i,j,4,inod)*table_hh(7,7)/drho
            fma(mw+1,ml+8)=fma(mw+1,ml+8) &
     &                    +fmd(i,j,1,inod)*table_hh(3,4)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(7,4)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(3,8)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(7,8)

            fma(mw-9,ml+9)=fma(mw-9,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hh(4,1)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(8,1)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(4,5)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(8,5)
            fma(mw-8,ml+9)=fma(mw-8,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hh(4,2)*drho**3 &
     &                    +fmd(i,j,2,inod)*table_hh(8,2)*drho**2 &
     &                    +fmd(i,j,3,inod)*table_hh(4,6)*drho**2 &
     &                    +fmd(i,j,4,inod)*table_hh(8,6)*drho
            fma(mw-1,ml+9)=fma(mw-1,ml+9) &
     &                    +fmd(i,j,1,inod)*table_hh(4,3)*drho**2 &
     &                    +fmd(i,j,2,inod)*table_hh(8,3)*drho &
     &                    +fmd(i,j,3,inod)*table_hh(4,7)*drho &
     &                    +fmd(i,j,4,inod)*table_hh(8,7)
            fma(mw  ,ml+9)=fma(mw  ,ml+9) &
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
         do mw=1,mwmax-2
            fma(mw,3) = fma(mw,3)+ci*nth*fma(mw+2,1)
         end do
         fvb(3)=fvb(3)+ci*nth*fvb(1)
         do mw=1,mwmax
            fma(mw,1) = 0.D0
            fma(mw,5) = 0.d0
            fma(mw,7) = 0.d0
         enddo
         fma(mc,1)=-ci*nth
         fma(mc+2,1)=1.D0
         fma(mc,5)=1.d0
         fma(mc,7)=1.d0
         fvb(1)=0.d0
         fvb(5)=0.d0
         fvb(7)=0.d0
      else
         do mw=1,mwmax
            fma(mw,1) = 0.d0
            fma(mw,3) = 0.d0
            fma(mw,5) = 0.d0
            fma(mw,7) = 0.d0
         enddo
         fma(mc,1)=1.d0
         fma(mc,3)=1.d0
         fma(mc,5)=1.d0
         fma(mc,7)=1.d0
         fvb(1)=0.d0
         fvb(3)=0.d0
         fvb(5)=0.d0
         fvb(7)=0.d0
      endif

      do mw=1,mwmax
         fma(mw,mlmax-7) = 0.d0
         fma(mw,mlmax-5) = 0.d0
         fma(mw,mlmax-3) = 0.d0
         fma(mw,mlmax-1) = 0.d0
      enddo
      fma(mc,mlmax-7) = 1.d0
      fma(mc,mlmax-5) = 1.d0
      fma(mc,mlmax-3) = 1.d0
      fma(mc,mlmax-1) = 1.d0
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
        divj=-ci*(rkth*(1.D0-angl)+rkph*angl)
        fvb(8*(nr-1)+2)=-(rho(nr+1)-rd)*divj/(rho(nr+1)-rho(nr))
        fvb(8*(nr-1)+3)=-(1.d0-angl)
        fvb(8*(nr-1)+5)=-angl
      do nr=nrd+1,nrmax-1
         fvb(8*(nr-1)+2)=-divj
      end do

      return
      end subroutine fem_calc_15
