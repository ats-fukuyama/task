!     $Id$

!---- interface for wm parameter

      subroutine get_wmparm(rr_out,ra_out,rf_out,nth0_out,nph0_out)
      
      include '../wm/wmcomm.inc'
      real(8),intent(out):: rr_out,ra_out,rf_out
      integer,intent(out):: nth0_out,nph0_out
      rr_out=rr
      ra_out=ra
      rf_out=real(crf)
      nth0_out=nth0
      nph0_out=nph0
      return
      end subroutine get_wmparm

      subroutine wmfem(nrmax,nthmax,nphmax,nsmax,rho,cef,cpp,cpa)

      implicit none
      integer,intent(in):: nrmax,nthmax,nphmax,nsmax
      real(8),dimension(nrmax),intent(in):: rho
      complex(8),dimension(3,nthmax,nphmax,nrmax),intent(out):: cef
      complex(8),dimension(nthmax,nphmax,nrmax,nsmax),intent(out):: cpp
      complex(8),dimension(nthmax,nphmax),intent(out):: cpa

      integer:: nbsmax,mlmax,mwmax
      complex(8),dimension(:,:),allocatable:: fma
      complex(8),dimension(:,:,:),allocatable:: fmpsa
      complex(8),dimension(:),allocatable:: fvb,fvx

      nbsmax=nthmax*nphmax  ! size of block matrix
      mlmax=6*nbsmax*nrmax      ! length of coeffient matrix and source vector
      mwmax=4*6*nbsmax-1        ! width of coefficient matrix

      call wmfem_exec

      return

      contains

!---- main routine of wmfem ----

      subroutine wmfem_exec

      integer:: mc,ml,mw,ierr,nr,nth,nph,ns,idx,ivec,idv
      complex(8):: csum

      mc=(mwmax+1)/2    ! diagonal position in mw

!     allocate matrix and vector

      call wmfem_allocate

      if(nrmax.eq.0) return   ! matrix and vector deallocated

!     calculate matrix
!        obtain metric
!        obtain dielectric tensor
!        calculate element matrix
!        FFT
!     calculate RHS vector

      call wmfem_calculate

!     solve matrix

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
      enddo
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      if(ierr.ne.0) write(6,*) '# ierr= ',ierr

!     calculate E field

      do nr=1,nrmax
      do nph=1,nphmax
      do nth=1,nthmax
         idx=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
         cef(1,nth,nph,nr)=fvx(idx+1)
         cef(2,nth,nph,nr)=fvx(idx+3)
         cef(3,nth,nph,nr)=fvx(idx+5)
      enddo
      enddo
      enddo

!     calculate power

      do ns=1,nsmax
      do nr=1,nrmax
      do nph=1,nphmax
      do nth=1,nthmax
         idx=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
         csum=0.d0
         do ivec=1,6
         do mw=1,mwmax
            idv=min(max(1,idx+ivec+mw-mc),mlmax)
            csum=csum+fmpsa(mw,idx+ivec,ns)*fvx(idv)
         enddo
         enddo
         cpp(nth,nph,nr,ns)=conjg(fvx(idx))*csum
      enddo
      enddo
      enddo
      enddo

!     calculate antenna impedance

      do nph=1,nphmax
      do nth=1,nthmax
         cpa(nth,nph)=0.d0
         do nr=1,nrmax
            idx=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
            do ivec=1,6
               cpa(nth,nph)=cpa(nth,nph)
     &              +conjg(fvx(idx+ivec))*fvb(idx+ivec)
            enddo
         enddo
      enddo
      enddo

      return
      end subroutine wmfem_exec

!     ----- allocate -----

      subroutine wmfem_allocate

      implicit none
      integer,save:: mwmax_save=0,mlmax_save=0,nsmax_save=0

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(allocated(fma)) deallocate(fma)
         if(mlmax.ne.0) allocate(fma(mwmax,mlmax))
      endif
      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save).or.
     &   (nsmax.ne.nsmax_save)) then
         if(allocated(fmpsa)) deallocate(fmpsa)
         if(mlmax.ne.0) allocate(fmpsa(mwmax,mlmax,nsmax))
      endif
      if(mlmax.ne.mlmax_save) then
         if(allocated(fvb)) deallocate(fvb)
         if(allocated(fvx)) deallocate(fvx)
         if(mlmax.ne.0) allocate(fvb(mlmax))
         if(mlmax.ne.0) allocate(fvx(mlmax))
      endif
      mwmax_save=mwmax
      mlmax_save=mlmax
      nsmax_save=nsmax
      end subroutine wmfem_allocate


!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate

      implicit none
      complex(8),dimension(:,:,:,:,:),allocatable:: fmc,fmd
      real(8):: drho,rkth,rkph,factor,rkth0,rho0
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nbs,nbsmax,nth,nph
      real(8):: rr,ra,rf
      integer:: nth0,nph0
      complex(8),parameter:: ci=(0.d0,1.d0)
      integer:: id_base=1
      real(8):: angl=0.d0

      call get_wmparm(rr,ra,rf,nth0,nph0)
      factor=rf**2
      nbsmax=nthmax*nphmax
      allocate(fmc(3,3,4,nbsmax,4))
      allocate(fmd(3,3,4,nbsmax,4))

      write(6,*) 'nbsmax=',nbsmax
      write(6,*) 'factor=',factor
      write(6,*) 'nth0=',nth0
      write(6,*) 'nthmax=',nthmax
      write(6,*) 'ra=',ra

      do ml=1,mlmax
         fvb(ml)=0.d0
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
      enddo

      do nr=1,nrmax-1
         drho=rho(nr+1)-rho(nr)

         do inod=1,4
            do nbs=1,nbsmax
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmc(i,j,k,nbs,inod)=0.d0
                     enddo
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
            
            do nph=1,nphmax
               do nth=1,nthmax
                  nbs=nthmax*(nph-1)+nth

                  if(nthmax.eq.1) then
                     rkth=nth0/(ra*rho0)
                  else
                     rkth=(nth-nthmax/2+nth0)/(ra*rho0)
                  endif
                  rkth0=1.d0/rho0
                  if(nphmax.eq.1) then
                     rkph=nph0/rr
                  else
                     rkph=(nph-nphmax/2+nph0)/rr
                  endif

                  fmc(1,1,1,nbs,inod)= rho0*(factor-rkph**2-rkth**2)
                  fmc(1,2,1,nbs,inod)= rho0*(-ci*rkth*rkth0)
                  fmc(2,1,1,nbs,inod)= rho0*(+ci*rkth*rkth0)
                  fmc(2,2,1,nbs,inod)= rho0*(factor-rkph**2-rkth0**2)
                  fmc(2,3,1,nbs,inod)= rho0*(rkth*rkph)
                  fmc(3,2,1,nbs,inod)= rho0*(rkth*rkph)
                  fmc(3,3,1,nbs,inod)= rho0*(factor-rkth**2)
                  
                  fmc(2,1,2,nbs,inod)= rho0*( ci*rkth)
                  fmc(2,2,2,nbs,inod)= rho0*(  -rkth0)
                  fmc(3,1,2,nbs,inod)= rho0*( ci*rkph)
                  
                  fmc(1,2,3,nbs,inod)= rho0*(-ci*rkth)
                  fmc(1,3,3,nbs,inod)= rho0*(-ci*rkph)
                  fmc(2,2,3,nbs,inod)= rho0*(  -rkth0)
                  
                  fmc(1,1,4,nbs,inod)= 0.d0
                  fmc(2,2,4,nbs,inod)= rho0*(-1.d0)
                  fmc(3,3,4,nbs,inod)= rho0*(-1.d0)
               enddo
            enddo
         enddo

! ------ calculate coefficients of basis for profile from four points 

         do nbs=1,nbsmax
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmd(i,j,k,nbs,1)=fmc(i,j,k,nbs,1)
                     fmd(i,j,k,nbs,2)=0.5d0*(-11*fmc(i,j,k,nbs,1)
     &                                       +18*fmc(i,j,k,nbs,2)
     &                                       - 9*fmc(i,j,k,nbs,3)
     &                                       + 2*fmc(i,j,k,nbs,4))
                     fmd(i,j,k,nbs,3)=fmc(i,j,k,nbs,4)
                     fmd(i,j,k,nbs,4)=0.5d0*(- 2*fmc(i,j,k,nbs,1)
     &                                       + 9*fmc(i,j,k,nbs,2)
     &                                       -18*fmc(i,j,k,nbs,3)
     &                                       +11*fmc(i,j,k,nbs,4))
                  enddo
               enddo
            enddo
         enddo

         if(id_base.eq.0) then
            call fem_hhh(nr,fmd,nbsmax,drho)
         else
            call fem_hqq(nr,fmd,nbsmax,drho)
         endif
      enddo

      mc=(mwmax+1)/2

      nr=1
      do nph=1,nphmax
         do nth=1,nthmax
            nbs=nthmax*(nph-1)+nth
            ml=6*nbsmax*(nr-1)+6*(nbs-1)
            
            if(nth-nthmax/2+nth0.eq.0) then
               do mw=1,mwmax
                  fma(mw,ml+3) = 0.d0
               enddo
               fma(mc,ml+3)=1.d0
            elseif(abs(nth-nthmax/2+nth0).eq.1) then
               do mw=1,mwmax
                  fma(mw,ml+3) = 0.d0
                  fma(mw,ml+5) = 0.d0
               enddo
               fma(mc-2,ml+3)=1.d0
               fma(mc  ,ml+3)=ci*(nth-nthmax/2+nth0)
               fma(mc  ,ml+5)=1.d0
            else
               do mw=1,mwmax
                  fma(mw,ml+3) = 0.d0
                  fma(mw,ml+5) = 0.d0
               enddo
               fma(mc,ml+3)=1.d0
               fma(mc,ml+5)=1.d0
            endif
         enddo
      enddo

      nr=nrmax
      do nph=1,nphmax
         do nth=1,nthmax
            nbs=nthmax*(nph-1)+nth
            ml=6*nbsmax*(nr-1)+6*(nbs-1)
            if(id_base.eq.0) then
               do mw=1,mwmax
                  fma(mw,ml+3) = 0.d0
                  fma(mw,ml+5) = 0.d0
               enddo
               fma(mc,ml+3) = 1.d0
               fma(mc,ml+5) = 1.d0
            else
               do mw=1,mwmax
                  fma(mw,ml+1) = 0.d0
                  fma(mw,ml+2) = 0.d0
                  fma(mw,ml+3) = 0.d0
                  fma(mw,ml+5) = 0.d0
               enddo
               fma(mc,ml+1) = 1.d0
               fma(mc,ml+2) = 1.d0
               fma(mc,ml+3) = 1.d0
               fma(mc,ml+5) = 1.d0
            endif
         enddo
      enddo

!------      fvb

      do nr=1,nrmax-1
         if((0.85d0-rho(nr))*(rho(nr+1)-0.85d0).ge.0.d0) then
            do nph=1,nphmax
               do nth=1,nthmax
                  nbs=nthmax*(nph-1)+nth
                  ml=6*nbsmax*(nr-1)+6*(nbs-1)
                  fvb(ml+3)=(rho(nr+1)-0.85d0)*(1.d0-angl)
                  fvb(ml+3)=(0.85d0-rho(nr)  )*(1.d0-angl)
                  fvb(ml+5)=(rho(nr+1)-0.85d0)*angl
                  fvb(ml+5)=(0.85d0-rho(nr)  )*angl
               enddo
            enddo
         endif
      enddo
      deallocate(fmc)
      deallocate(fmd)

!      do ml=1,mlmax
!         write(6,'(1P6E12.4)') (fma(mw,ml),mw=1,mwmax)
!      enddo

      return
      end subroutine wmfem_calculate

!---- FEM cubic hermit ---

      subroutine fem_hhh(nr,fmd,nbsmax,drho)

      use libfem_mod
      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nbsmax,4),intent(in):: fmd
      integer,intent(in):: nbsmax
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nbs,inod,ml,mw

      mr=6*nbsmax  ! line interval between radial points 
      mc=(mwmax+1)/2


         do nbs=1,nbsmax
         do j=1,3
         do i=1,3
            ml=6*nbsmax*(nr-1)+6*(nbs-1)+2*(i-1)+1
            mw=mc+2*(j-1)-6*(nbs-1)-2*(i-1)
            do inod=1,4

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+Mr+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,6)*drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,7)
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,8)*drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,8)

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,5)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+mr+1)=fma(mw  ,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,4)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,4)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,8)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,8)*drho

               enddo
            enddo
         enddo
      enddo

      return
      end subroutine fem_hhh

!---- FEM cubic hermit + quadratic ---

      subroutine fem_hqq(nr,fmd,nbsmax,drho)

      use libfem_mod
      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nbsmax,4),intent(in):: fmd
      integer,intent(in):: nbsmax
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nbs,inod,ml,mw

      mr=6*nbsmax  ! line interval between radial points 
      mc=(mwmax+1)/2

      write(6,*) nr,mr,mc,nbsmax

      write(6,'(1P6E12.4)') fmd
      pause

         do nbs=1,nbsmax
         do j=1,3
         do i=1,3
            ml=6*nbsmax*(nr-1)+6*(nbs-1)+2*(i-1)+1
            mw=mc+2*(j-1)-6*(nbs-1)-2*(i-1)
            write(6,*) 'ml,mw=',ml,mw
            pause
            do inod=1,4

! (1,1) **********
            if(i.eq.1.and.j.eq.1) then

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,1,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,4,1)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,1,4)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,4,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,1,2)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,4,2)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,1,5)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,4,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,1,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,4,3)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,1,6)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,4,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,2,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,5,1)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,2,4)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,5,4)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,2,2)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,5,2)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,2,5)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,5,5)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,2,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,5,3)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,2,6)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,5,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,3,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,6,1)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,3,4)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,6,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,3,2)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,6,2)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,3,5)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,6,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqq(inod,3,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqq(inod,6,3)
     &              +fmd(i,j,3,nbs,inod)*table_hqq(inod,3,6)
     &              +fmd(i,j,4,nbs,inod)*table_hqq(inod,6,6)/drho

! (1,*) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,1,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,4,1)
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,1,5)
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,4,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,4,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,1,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,4,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,1,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,4,3)
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,1,7)
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,4,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,4,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,1,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,4,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,2,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,5,1)
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,2,5)
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,5,5)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,2,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,5,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,2,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,5,6)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,2,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,5,3)
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,2,7)
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,5,7)/drho
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,2,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,5,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,2,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,5,8)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,3,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,6,1)
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,3,5)
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,6,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,6,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,3,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,6,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,3,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,6,3)
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,3,7)
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,6,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hqh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hqh(inod,6,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hqh(inod,3,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hqh(inod,6,8)

! (*,1) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,1,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,5,1)
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,1,4)
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,5,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,1,2)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,5,2)
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,1,5)
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,5,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,1,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,5,3)
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,1,6)
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,5,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,2,1)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,6,1)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,2,4)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,6,4)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,2,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,6,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,2,5)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,6,5)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,2,3)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,6,3)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,2,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,6,6)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,3,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,7,1)
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,3,4)
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,7,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,3,2)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,7,2)
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,3,5)
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,7,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,3,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,7,3)
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,3,6)
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,7,6)/drho

            fma(mw-mw-1,ml+mr+1)=fma(mw-mw-1,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,4,1)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,8,1)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,4,4)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,8,4)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,4,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,8,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,4,5)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,8,5)
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhq(inod,4,3)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhq(inod,8,3)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhq(inod,4,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhq(inod,8,6)

! (*,*) **********
            else
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  )
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,6)*drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,7)
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,6,8)*drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,7,8)

            fma(mw-mw-1,ml+mr+1)=fma(mw-mw-1,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,5)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+mr+1)=fma(mw  ,ml+mr+1)
     &              +fmd(i,j,1,nbs,inod)*table_hhh(inod,4,4)*drho**3
     &              +fmd(i,j,2,nbs,inod)*table_hhh(inod,8,4)*drho**2
     &              +fmd(i,j,3,nbs,inod)*table_hhh(inod,4,8)*drho**2
     &              +fmd(i,j,4,nbs,inod)*table_hhh(inod,8,8)*drho
            endif

               enddo
            enddo
         enddo
      enddo
         write(6,'(1P6E12.4)') (fma(mw,1),mw=1,mwmax)
         write(6,'(1P6E12.4)') (fma(mw,2),mw=1,mwmax)
         write(6,'(1P6E12.4)') (fma(mw,3),mw=1,mwmax)
         write(6,'(1P6E12.4)') (fma(mw,4),mw=1,mwmax)
         write(6,'(1P6E12.4)') (fma(mw,5),mw=1,mwmax)
         write(6,'(1P6E12.4)') (fma(mw,6),mw=1,mwmax)
         pause
       return
      end subroutine fem_hqq

      end subroutine wmfem
