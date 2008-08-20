!     $Id$

!     ***** wmfem main routine *****

      subroutine wmfem(nrmax,nthmax,nphmax,nsmax,rho,cef,cpp,cpa)

      implicit none
      integer,intent(in):: nrmax,nthmax,nphmax,nsmax
      real(8),dimension(nrmax),intent(in):: rho
      complex(8),dimension(3,nthmax,nphmax,nrmax),intent(out):: cef
      complex(8),dimension(nthmax,nphmax,nthmax,nphmax,nrmax,0:nsmax),
     &     intent(out):: cpp
      complex(8),dimension(nthmax,nphmax),intent(out):: cpa
      complex(8),parameter:: ci=(0.d0,1.d0)
      real(8),parameter:: pi = 3.14159265358979D0
      real(8),parameter:: vc = 2.99792458 D8

      complex(8):: crf
      integer:: nth0,nph0
      integer:: idbgwm

      integer:: nfcmax,mlmax,mwmax,nthmax2,nphmax2,nfcmax2
      complex(8),dimension(:,:),allocatable,save:: fma,fmax !(mwmax,mlmax)
      complex(8),dimension(:,:,:),allocatable,save:: fms !(mwmax,mlmax,0:nsmax)
      complex(8),dimension(:),allocatable,save:: fvb,fvx !(mlmax)

      real(8),dimension(:,:,:,:,:),allocatable :: gma,mma 
     &                                          !(3,3,nthmax,nphmax,nrmax)
      real(8),dimension(:,:,:),allocatable:: gj !(nthmax,nphmax,nrmax)

      integer,dimension(:),allocatable,save :: nthnfc,nthfnfc
      integer,dimension(:),allocatable,save :: nphnfc,nphfnfc
      integer,dimension(:),allocatable,save :: nthnfc2,nthfnfc2
      integer,dimension(:),allocatable,save :: nphnfc2,nphfnfc2
!      integer:: nfc,nfc2
!      integer:: nthf1,nphf1,nthf2,nphf2,nthfdiff,nphfdiff

      nfcmax=nthmax*nphmax      ! size of block matrix 
                                !    (number of Fourier components)
      mlmax=6*nfcmax*nrmax      ! length of coeffient matrix and source vector
      mwmax=4*6*nfcmax-1        ! width of coefficient matrix
      nthmax2=2*nthmax          ! number of poloidal modes of coefficients
      nphmax2=2*nphmax          ! number of toroidal modes of coefficients
      nfcmax2=nthmax2*nphmax2   ! size of block matrix of coefficients

      call get_wmparm(crf,nth0,nph0,idbgwm)

!     allocate matrix and vector

      call wmfem_allocate

      if(nrmax.eq.0) return   ! matrix and vector deallocated

!     calculate metric and convaersion tensor

      call wmfem_metric(gma,mma,gj)

!     setup Fourier component index

      call wmfem_setup_index

!     calculate matrix

      call wmfem_calculate

!     solve matrix

      call wmfem_solve

      return

      contains

!     ----- allocate -----

      subroutine wmfem_allocate

      implicit none
      integer,save:: nrmax_save=0,nthmax_save=0,nphmax_save=0
      integer,save:: mwmax_save=0,mlmax_save=0
      integer,save:: nsmax_save=0,nfcmax_save=0

      if((nrmax.ne.nrmax_save).or.(nthmax.ne.nthmax_save)
     &                        .or.(nphmax.ne.nphmax_save)) then
         if(allocated(gma)) deallocate(gma)
         if(allocated(mma)) deallocate(mma)
         if(allocated(gj)) deallocate(gj)
         allocate(gma(3,3,nthmax2,nphmax2,nrmax))
         allocate(mma(3,3,nthmax2,nphmax2,nrmax))
         allocate(gj(nthmax2,nphmax2,nrmax))
      endif

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(allocated(fma)) deallocate(fma)
         if(allocated(fmax)) deallocate(fmax)
         allocate(fma(mwmax,mlmax))
         allocate(fmax(mwmax,mlmax))
      endif
      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save).or.
     &   (nsmax.ne.nsmax_save)) then
         if(allocated(fms)) deallocate(fms)
         allocate(fms(mwmax,mlmax,0:nsmax))
      endif
      if(mlmax.ne.mlmax_save) then
         if(allocated(fvb)) deallocate(fvb)
         if(allocated(fvx)) deallocate(fvx)
         allocate(fvb(mlmax))
         allocate(fvx(mlmax))
      endif
      if(nfcmax.ne.nfcmax_save) then
         if(allocated(nthnfc)) deallocate(nthnfc)
         if(allocated(nphnfc)) deallocate(nphnfc)
         if(allocated(nthfnfc)) deallocate(nthfnfc)
         if(allocated(nphfnfc)) deallocate(nphfnfc)
         if(nfcmax.ne.0) allocate(nthnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nphnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nthfnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nphfnfc(nfcmax))

         if(allocated(nthnfc2)) deallocate(nthnfc2)
         if(allocated(nphnfc2)) deallocate(nphnfc2)
         if(allocated(nthfnfc2)) deallocate(nthfnfc2)
         if(allocated(nphfnfc2)) deallocate(nphfnfc2)
         if(nfcmax2.ne.0) allocate(nthnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nphnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nthfnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nphfnfc2(nfcmax2))
      endif
      mwmax_save=mwmax
      mlmax_save=mlmax
      nsmax_save=nsmax
      nfcmax_save=nfcmax
      end subroutine wmfem_allocate

!---- setup Fourier component inxex ----

      subroutine wmfem_setup_index

      integer:: nth,nph,nfc,nfc2

!     setup an array of mode number

      if(nthmax.eq.1) then
         do nph=1,nphmax
            nfc=nph
            nthnfc(nfc)=1
            nthfnfc(nfc)=nth0
         enddo
      else
         do nph=1,nphmax
         do nth=1,nthmax
            nfc=nthmax*(nph-1)+nth
            nthnfc(nfc)=nth-1
         enddo
         do nth=1,nthmax/2
            nfc=nthmax*(nph-1)+nth
            nthfnfc(nfc)=nth0+nth-1
            nfc=nthmax*(nph-1)+nthmax+1-nth
            nthfnfc(nfc)=nth0-nth
         enddo
         enddo
      endif
      if(nphmax.eq.1) then
         do nth=1,nthmax
            nfc=nth
            nphnfc(nfc)=1
            nphfnfc(nfc)=nph0
         enddo
      else
         do nth=1,nthmax
         do nph=1,nphmax
            nfc=nthmax*(nph-1)+nth
            nphnfc(nfc)=nph
         enddo
         do nph=1,nphmax/2
            nfc=nthmax*(nph-1)+nth
            nphfnfc(nfc)=nph0+nph-1
            nfc=nthmax*(nphmax-nph)+nth
            nphfnfc(nfc)=nph0-nph
         enddo
         enddo
      endif

      if(nthmax.eq.1) then
         do nph=1,nphmax2
            nfc2=nph
            nthnfc2(nfc2)=1
            nthfnfc2(nfc2)=0
         enddo
      else
         do nph=1,nphmax2
         do nth=1,nthmax2
            nfc2=nthmax2*(nph-1)+nth
            nthnfc2(nfc2)=nth-1
         enddo
         do nth=1,nthmax
            nfc2=nthmax2*(nph-1)+nth
            nthfnfc2(nfc2)=nth-1
            nfc2=nthmax2*(nph-1)+nthmax2+1-nth
            nthfnfc2(nfc2)=-nth
         enddo
         enddo
      endif
      if(nphmax.eq.1) then
         do nth=1,nthmax2
            nfc2=nth
            nphnfc2(nfc2)=1
            nphfnfc2(nfc2)=0
         enddo
      else
         do nth=1,nthmax2
         do nph=1,nphmax2
            nfc2=nthmax2*(nph-1)+nth
            nphnfc2(nfc2)=nph
         enddo
         do nph=1,nphmax
            nfc2=nthmax2*(nph-1)+nth
            nphfnfc2(nfc2)=nph-1
            nfc2=nthmax2*(nphmax2-nph)+nth
            nphfnfc2(nfc2)=-nph
         enddo
         enddo
      endif

      return
      end subroutine wmfem_setup_index

!---- Solve matrix equation to solve ----

      subroutine wmfem_solve

      integer:: mc,ml,mw,ierr,nr,nth,nph,ns,nr1,nth1,nph1,i,j
      complex(8):: csum,csums,csum1,csum2

      mc=(mwmax+1)/2    ! diagonal position in mw

!     solve matrix

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
      enddo
      do ml=1,mlmax
         do mw=1,mwmax
            fmax(mw,ml)=fma(mw,ml)
         enddo
      enddo

      write(6,*) 'mlmax,mwmax=',mlmax,mwmax
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      if(ierr.ne.0) write(6,*) '# ierr= ',ierr

!     calculate E field

      do nr=1,nrmax
      do nph=1,nphmax
      do nth=1,nthmax
         ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
         cef(1,nth,nph,nr)=fvx(ml+1)
         cef(2,nth,nph,nr)=fvx(ml+3)
         cef(3,nth,nph,nr)=fvx(ml+5)
      enddo
      enddo
      enddo

!     calculate power

      mc=(mwmax+1)/2
      do ns=0,nsmax
         csums=0.d0
      do nr=1,nrmax-1
      do nph=1,nphmax
      do nth=1,nthmax
         ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
         do nph1=1,nphmax
         do nth1=1,nthmax
            csum=0.d0
            do nr1=nr-1,nr+1
               mw=mc+6*nthmax*nphmax*(nr1-nr)+6*nthmax*(nph1-nph)
     &              +6*(nth1-nth)
               if(ml+mw-mc+1.ge.1.and.ml+mw-mc+6.le.mlmax) then
               do i=1,6
               do j=1,6
                  csum=csum
     &              +conjg(fvx(ml+i))
c     &                      *fms(mw+j-i,ml+i,ns)
     &                           *fvx(ml+mw-mc+j)
               enddo
               enddo
               endif
            enddo
            csums=csums+csum
            cpp(nth,nph,nth1,nph1,nr,ns)=csum
C            write(6,'(6I5,1P2E12.4)') nth,nph,nth1,nph1,nr,ns,
C     &           cpp(nth,nph,nth1,nph1,nr,ns)
         enddo
         enddo
      enddo
      enddo
      enddo
!         write(6,'(A,I5,1P2E12.4)') 'ns,csums=',ns,csums
      enddo

!     calculate antenna impedance

      do nph=1,nphmax
      do nth=1,nthmax
         cpa(nth,nph)=0.d0
         do nr=1,nrmax-1
            ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
            do i=1,6
               cpa(nth,nph)=cpa(nth,nph)+conjg(fvx(ml+i))*fvb(ml+i)
            enddo
         enddo
C         write(6,'(2I5,1P2E12.4)') nth,nph,cpa(nth,nph)
      enddo
      enddo

      return
      end subroutine wmfem_solve

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate

      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmd
      complex(8),dimension(nphmax,nthmax,3):: fvb_nr
      complex(8):: cfactor
      real(8):: drho,rkth,rkph,rkth0,rho0
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,nthf,nphf
      integer:: ns,nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      integer:: id_base=1
      real(8):: angl=0.d0

      cfactor=(2*pi*crf*1.d6)**2/vc**2

      do ns=0,nsmax             ! loop for vacuum and plasma species

         do ml=1,mlmax             ! clear fma
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
         enddo

         do nr=1,nrmax-1        ! loop for elements
            drho=rho(nr+1)-rho(nr)

            if(ns.eq.0) then
               if(idbgwm.eq.0) then
                  call wmfem_calculate_vacuum_g(nr,fmd)
               else if(idbgwm.eq.1) then
c                  call wmfem_calculate_vacuum_l(nr,fmd)
               else if(idbgwm.eq.2) then
                  call wmfem_calculate_vacuum_gc(nr,fmd)
               else if(idbgwm.eq.3) then
                  call wmfem_calculate_vacuum_lc(nr,fmd)
               endif
            else
               call wmfem_calculate_plasma(nr,ns,fmd)
            endif

c$$$            if(ns.eq.0) then
c$$$            do nfc=1,nfcmax
c$$$               nth=nthfc(nfc)
c$$$               nph=nphfc(nfc)
c$$$               write(6,'(5I8)') nr,nth,nph
c$$$               write(6,'(1P2E12.4,2X,1P2E12.4,2X,1P2E12.4)')
c$$$     &              ((fmd(i,j,1,nfc,nfc,1),j=1,3),i=1,3)
c$$$               write(6,'(1P2E12.4,2X,1P2E12.4,2X,1P2E12.4)')
c$$$     &              ((fmd(i,j,2,nfc,nfc,1),j=1,3),i=1,3)
c$$$               write(6,'(1P2E12.4,2X,1P2E12.4,2X,1P2E12.4)')
c$$$     &              ((fmd(i,j,3,nfc,nfc,1),j=1,3),i=1,3)
c$$$               write(6,'(1P2E12.4,2X,1P2E12.4,2X,1P2E12.4)')
c$$$     &              ((fmd(i,j,4,nfc,nfc,1),j=1,3),i=1,3)
c$$$            enddo
c$$$            endif

! ------ calculate coefficients of basis for profile from four points 

            do nfc1=1,nfcmax
               do nfc2=1,nfcmax
                  do k=1,4
                     do j=1,3
                        do i=1,3
                           fmd1=fmd(i,j,k,nfc1,nfc2,1)
                           fmd2=fmd(i,j,k,nfc1,nfc2,2)
                           fmd3=fmd(i,j,k,nfc1,nfc2,3)
                           fmd4=fmd(i,j,k,nfc1,nfc2,4)
                           fmd(i,j,k,nfc1,nfc2,1)=fmd1
                           fmd(i,j,k,nfc1,nfc2,2)
     &                       =0.5d0*(-11*fmd1+18*fmd2- 9*fmd3+ 2*fmd4)
                           fmd(i,j,k,nfc1,nfc2,3)=fmd4
                           fmd(i,j,k,nfc1,nfc2,4)
     &                       =0.5d0*(- 2*fmd1+ 9*fmd2-18*fmd3+11*fmd4)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            
            if(id_base.eq.0) then
               call fem_hhh(nr,fmd,nfcmax,drho)
            else
               call fem_hqq(nr,fmd,nfcmax,drho)
            endif
         enddo

         do ml=1,mlmax
         do mw=1,mwmax
            fms(mw,ml,ns)=fma(mw,ml)
         enddo
         enddo
      enddo

      do ml=1,mlmax
      do mw=1,mwmax
         csum=0.d0
         do ns=0,nsmax
            csum=csum+fms(mw,ml,ns)
         enddo
         fma(mw,ml)=csum
      enddo
      enddo

!----- boundary conditions -----

      mc=(mwmax+1)/2

      nr=1
      do nfc=1,nfcmax
         nthf=nthfnfc(nfc)
         ml=6*nfcmax*(nr-1)+6*(nfc-1)
         if(nthf.eq.0) then
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
            enddo
            fma(mc,ml+3)=1.d0
         elseif(abs(nthf).eq.1) then
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
               fma(mw,ml+5) = 0.d0
            enddo
            fma(mc-2,ml+3)=1.d0
            fma(mc  ,ml+3)=ci*nth
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

      nr=nrmax
      do nfc=1,nfcmax
         ml=6*nfcmax*(nr-1)+6*(nfc-1)
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

!------      fvb

      do ml=1,mlmax
         fvb(ml)=0.d0
      enddo
         
      do nr=1,nrmax-1
         call get_wmfvb(nr,fvb_nr)
         do nfc=1,nfcmax
            nth=nthnfc(nfc)
            nph=nthnfc(nfc)
            ml=6*nfcmax*(nr-1)+6*(nfc-1)
            fvb(ml+1)=fvb_nr(nph,nth,1)
            fvb(ml+3)=fvb_nr(nph,nth,2)
            fvb(ml+5)=fvb_nr(nph,nth,3)
!     write(6,'(4I5,1P2E12.4)') nr,nph,nth,ml+1,fvb(ml+1)
!     write(6,'(4I5,1P2E12.4)') nr,nph,nth,ml+3,fvb(ml+3)
!     write(6,'(4I5,1P2E12.4)') nr,nph,nth,ml+5,fvb(ml+5)
         enddo
      enddo

!      do ml=1,mlmax
!         write(6,'(1P6E12.4)') (fma(mw,ml),mw=1,mwmax)
!      enddo
!      write(6,'(1P6E12.4)') (fvb(ml),ml=1,mlmax)

      return
      end subroutine wmfem_calculate

!----- calculate coefficint matrix fma (E rho,theta,phi )-----

      subroutine wmfem_calculate_vacuum_g(nr,fmd)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8),dimension(3,3,3,3,nthmax2,nphmax2):: fmv1
      complex(8),dimension(3,3,3,nthmax2,nphmax2):: fmv2,fmv3
      complex(8),dimension(3,3,nthmax2,nphmax2):: fmv4
      real(8),dimension(4):: rhol
      real(8):: drho,rkth,rkph,rkth0,rho0
      integer:: ml,mw,mc,nvmax,i,j,k,inod,nfc,nf1,nf2,nth,nph,nthx
      integer:: ns,nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      integer:: nph1,nph2,nphx1,nphx2
      integer:: nphdiff,nphxdiff
      integer:: nth1,nth2,nthx1,nthx2
      integer:: nthdiff,nthxdiff
      integer:: nthf1,nphf1,nthf2,nphf2
      integer:: imn1,imn2

      call wmfem_calculate_vacuum_sub(nr,fmv1,fmv2,fmv3,fmv4)

      do i=1,3
         do j=1,3
            do imn1=1,3
            do imn2=1,3
               do nph=1,nphmax
                  do nth=1,nthmax
                     fv1(nth,nph)=fmv1(i,j,imn1,imn2,nth,nph)
                  enddo
               enddo
               call wmsubfx(fv1,fv1f,nthmax,nphmax)
               do nph=1,nphmax
                  do nth=1,nthmax
                     fmv1(i,j,imn1,imn2,nth,nph)=fv1f(nth,nph)
                  enddo
               enddo
            enddo
            enddo    
            do imn1=1,3
               do nph=1,nphmax
                  do nth=1,nthmax
                     fv1(nth,nph)=fmv2(i,j,imn1,nth,nph)
                  enddo
               enddo
               call wmsubfx(fv1,fv1f,nthmax,nphmax)
               do nph=1,nphmax
                  do nth=1,nthmax
                     fmv2(i,j,imn1,nth,nph)=fv1(nth,nph)
                  enddo
               enddo
            enddo
            do imn1=1,3
               do nph=1,nphmax
                  do nth=1,nthmax
                     fv1(nth,nph)=fmv3(i,j,imn1,nth,nph)
                  enddo
               enddo
               call wmsubfx(fv1,fv1f,nthmax,nphmax)
               do nph=1,nphmax
                  do nth=1,nthmax
                     fmv3(i,j,imn1,nth,nph)=fv1(nth,nph)
                  enddo
               enddo
            enddo
            do nph=1,nphmax
               do nth=1,nthmax
                  fv1(nth,nph)=fmv4(i,j,nth,nph)
               enddo
            enddo
            call wmsubfx(fv1,fv1f,nthmax,nphmax)
            do nph=1,nphmax
               do nth=1,nthmax
                  fmv4(i,j,nth,nph)=fv1(nth,nph)
               enddo
            enddo
         enddo
      enddo

      do nfc1=1,nfcmax          ! Fit to fmd and adjust m and n
         nphf1=nphfnfc(nfc1)
         nthf1=nthfnfc(nfc1)
      do nfc2=1,nfcmax
         nphf2=nphfnfc(nfc2)
         nthf2=nthfnfc(nfc2)

!         nphfdiff=nphf2-nphf1+nphmax
!         nthfdiff=nthf2-nthf1+nthmax

         do inod=1,4
            do j=1,3
            do i=1,3
               fmd(i,j,1,nfc1,nfc2,inod)
     &              =fmv1(i,j,1,1,nthdiff,nphdiff)
     &              +fmv1(i,j,2,1,nthdiff,nphdiff)*nthf1
     &              +fmv1(i,j,3,1,nthdiff,nphdiff)*nphf1
     &              +fmv1(i,j,1,2,nthdiff,nphdiff)      *nthf2
     &              +fmv1(i,j,2,2,nthdiff,nphdiff)*nthf1*nthf2
     &              +fmv1(i,j,3,2,nthdiff,nphdiff)*nphf1*nthf2
     &              +fmv1(i,j,1,3,nthdiff,nphdiff)      *nphf2
     &              +fmv1(i,j,2,3,nthdiff,nphdiff)*nthf1*nphf2
     &              +fmv1(i,j,3,3,nthdiff,nphdiff)*nphf1*nphf2
               fmd(i,j,2,nfc1,nfc2,inod)
     &              =fmv2(i,j,1,nthdiff,nphdiff)
     &              +fmv2(i,j,2,nthdiff,nphdiff)*nthf1
     &              +fmv2(i,j,3,nthdiff,nphdiff)*nphf1
               fmd(i,j,3,nfc1,nfc2,inod)
     &              =fmv3(i,j,1,nthdiff,nphdiff)
     &              +fmv3(i,j,2,nthdiff,nphdiff)      *nthf2
     &              +fmv3(i,j,3,nthdiff,nphdiff)      *nphf2
               fmd(i,j,4,nfc1,nfc2,inod)
     &              =fmv4(i,j,nthdiff,nphdiff)
            enddo
            enddo
         enddo
      enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_g

!----- calculate vacuum matrix fms -----

      subroutine wmfem_calculate_vacuum_sub(nr,fmv1,fmv2,fmv3,fmv4)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,3,3,nthmax,nphmax),intent(out):: fmv1
      complex(8),dimension(3,3,3,nthmax,nphmax),intent(out):: fmv2,fmv3
      complex(8),dimension(3,3,nthmax,nphmax),intent(out):: fmv4
      complex(8),dimension(3,3,3):: cq
      complex(8),dimension(3,3):: cp
      integer:: i,j,k,l,nthm,nthp,nphm,nphp
      integer:: imn1,imn2
      integer:: nph,nth
      real(8):: dph,dth
      complex(8):: csum1,csum2,csum3,csum4,cfactor
      integer:: nrl
      real(8):: drhob,mma3b,mma2b,drhoa,mma3a,mma2a,mma30,mma20
      real(8):: mma3d,mma2d,gjl
      
      cfactor=(2*pi*crf*1.d6)**2/vc**2

C      write(6,'(A,3I5)') 'nr,nthmax,nphmax:',nr,nthmax,nphmax

      do nph=1,nphmax
         if(nph.eq.1) then
            nphm=nphmax
         else
            nphm=nph-1
         endif
         if(nph.eq.nphmax) then
            nphp=1
         else
            nphp=nph+1
         endif
         dph=2*pi/nphmax
      do nth=1,nthmax
         if(nth.eq.1) then
            nthm=nthmax
         else
            nthm=nth-1
         endif
         if(nth.eq.nthmax) then
            nthp=1
         else
            nthp=nth+1
         endif
         dth=2*pi/nthmax

         write(6,'(A,5I5)') 'nr,nph,nths:',nr,nph,nthm,nth,nthp
         write(6,'(1P2E12.4)') gj(nth,nph,nr)
         write(6,'(1P3E12.4)') mma(1,1,nth,nphp,nr),
     &                         mma(1,2,nth,nphp,nr),
     &                         mma(1,3,nth,nphp,nr)
         write(6,'(1P3E12.4)') mma(2,1,nth,nphp,nr),
     &                         mma(2,2,nth,nphp,nr),
     &                         mma(2,3,nth,nphp,nr)
         write(6,'(1P3E12.4)') mma(3,1,nth,nphp,nr),
     &                         mma(3,2,nth,nphp,nr),
     &                         mma(3,3,nth,nphp,nr)

         do j=1,3

            if(nr.eq.1) then
               nrl=3
               drhob=rho(3)
               mma3b=mma(3,j,nth,nph,nrl)
               mma2b=mma(2,j,nth,nph,nrl)
               nrl=2
               drhoa=rho(2)
               mma3a=mma(3,j,nth,nph,nrl)
               mma2a=mma(2,j,nth,nph,nrl)
               nrl=1
               mma30=mma(3,j,nth,nph,nrl)
               mma20=mma(2,j,nth,nph,nrl)
               mma3D=((mma3a-mma30)*drhob**2-(mma3b-mma30)*drhoa**2)
     &              /(drhoa*drhob*(drhob-drhoa))
               mma2D=((mma2a-mma20)*drhob**2-(mma2b-mma20)*drhoa**2)
     &              /(drhoa*drhob*(drhob-drhoa))
               gj(nth,nph,1)=gj(nth,nph,2)/9.d0
               gjl=gj(nth,nph,1)
            elseif(nr.eq.nrmax) then
               nrl=nrmax-2
               drhob=rho(nrl)-rho(nrmax)
               mma3b=mma(3,j,nth,nph,nrl)
               mma2b=mma(2,j,nth,nph,nrl)
               nrl=nrmax-1
               drhoa=rho(nrl)-rho(nrmax)
               mma3a=mma(3,j,nth,nph,nrl)
               mma2a=mma(2,j,nth,nph,nrl)
               nrl=nrmax
               mma30=mma(3,j,nth,nph,nrl)
               mma20=mma(2,j,nth,nph,nrl)
               mma3D=((mma3a-mma30)*drhob**2-(mma3b-mma30)*drhoa**2)
     &              /(drhoa*drhob*(drhob-drhoa))
               mma2D=((mma2a-mma20)*drhob**2-(mma2b-mma20)*drhoa**2)
     &              /(drhoa*drhob*(drhob-drhoa))
               gjl=gj(nth,nph,nr)
            else
               mma3D=(mma(3,j,nth,nph,nr+1)-mma(3,j,nth,nph,nr-1))
     &              /(rho(nr+1)-rho(nr-1))
               mma2D=(mma(2,j,nth,nph,nr+1)-mma(2,j,nth,nph,nr-1))
     &              /(rho(nr+1)-rho(nr-1))
               gjl=gj(nth,nph,nr)
            endif

         cq(1,j,1)=((mma(3,j,nthp,nph,nr)
     &              -mma(3,j,nthm,nph,nr))/dth
     &             -(mma(2,j,nth,nphp,nr)
     &              -mma(2,j,nth,nphm,nr))/dph )/gjl
         cq(1,j,2)=+ci*mma(3,j,nth,nph,nr)/gjl
         cq(1,j,3)=-ci*mma(2,j,nth,nph,nr)/gjl

         cq(2,j,1)=((mma(1,j,nth,nphp,nr)
     &              -mma(1,j,nth,nphm,nr))/dph
     &             -mma3d)/gjl
         cq(2,j,2)=0.d0
         cq(2,j,3)=+ci*mma(1,j,nth,nph,nr)/gjl

         cq(3,j,1)=(mma2d
     &             -(mma(1,j,nthp,nph,nr)
     &              -mma(1,j,nthm,nph,nr))/dth )/gjl
         cq(3,j,2)=-ci*mma(1,j,nth,nph,nr)/gjl
         cq(3,j,3)=0.d0

         cp(1,j)=0.d0
         cp(2,j)=-mma(3,j,nth,nph,nr)/gjl
         cp(3,j)= mma(2,j,nth,nph,nr)/gjl
      enddo

      write(6,*) 'cq(1)'
      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
     &     (cq(i,1,1),cq(i,2,1),cq(i,3,1),i=1,3)
      write(6,*) 'cq(2)'
      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
     &     (cq(i,1,2),cq(i,2,2),cq(i,3,2),i=1,3)
      write(6,*) 'cq(3)'
      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
     &     (cq(i,1,3),cq(i,2,3),cq(i,3,3),i=1,3)
      write(6,*) 'cp'
      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
     &     (cp(i,1),  cp(i,2),  cp(i,3),i=1,3)

      do imn2=1,3
      do imn1=1,3
      do j=1,3
      do i=1,3
         csum1=0.d0
         do k=1,3
         do l=1,3
            csum1=csum1+conjg(cq(k,i,imn1))
     &                       *gma(k,l,nth,nph,nr)
     &                       *cq(l,j,imn2) *gj(nth,nph,nr)
         enddo
         enddo
         if(i.eq.j.and.imn1.eq.1.and.imn2.eq.1) then
            fmv1(i,j,imn1,imn2,nth,nph)=csum1-cfactor*gj(nth,nph,nr)
         else
            fmv1(i,j,imn1,imn2,nth,nph)=csum1
         endif
      enddo
      enddo
      enddo
      enddo

      do imn1=1,3
      do j=1,3
      do i=1,3
         csum2=0.d0
         csum3=0.d0
         do k=1,3
         do l=1,3
            csum2=csum2+conjg(cp(k,i))
     &                       *gma(k,l,nth,nph,nr)
     &                       *cq(l,j,imn1) *gj(nth,nph,nr)
            csum3=csum3+conjg(cq(k,i,imn1))
     &                       *gma(k,l,nth,nph,nr)
     &                       *cp(l,j)  *gj(nth,nph,nr)
         enddo
         enddo
         fmv2(i,j,imn1,nth,nph)=csum2
         fmv3(i,j,imn1,nth,nph)=csum3
      enddo
      enddo
      enddo

      do j=1,3
      do i=1,3
         csum4=0.d0
         do k=1,3
         do l=1,3
            csum4=csum4+conjg(cp(k,i))
     &                       *gma(k,l,nth,nph,nr)
     &                       *cp(l,j) *gj(nth,nph,nr)
         enddo
         enddo
         fmv4(i,j,nth,nph)=csum4
      enddo
      enddo

      enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_sub

!----- calculate coefficint matrix fmd (E cylindrical rho,theta,phi)-----

      subroutine wmfem_calculate_vacuum_gc(nr,fmd)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmc
      complex(8),dimension(nthmax,nphmax):: fv1,fv1f
      real(8),dimension(4):: rhol
      complex(8):: cfactor
      real(8):: drho,rkth,rkph,rkth0,rho0
      integer:: ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,nthx,nphx
      integer:: ns,nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      integer:: nph1,nph2,nph1x,nph2x
      integer:: nphdiff,nphxdiff
      integer:: nth1,nth2,nth1x,nth2x
      integer:: nthdiff,nthxdiff
      integer:: nfcdiff
      real(8):: rr,ra,rb

      call get_wmparm1(rr,ra,rb)
      cfactor=(2*pi*crf*1.d6)**2/vc**2

      do inod=1,4               ! clear fmc
         do nfc2=1,nfcmax
            do nfc1=1,nfcmax
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmc(i,j,k,nfc1,nfc2,inod)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do inod=1,4               ! define local radius for different inod
         if(inod.eq.1) then
            if(nr.eq.1) then
               rhol(inod)=(8.d0*rho(nr)+rho(nr+1))/9.d0
            else
               rhol(inod)=rho(nr)
            endif
         elseif(inod.eq.2) then
            rhol(inod)=(2.d0*rho(nr)+rho(nr+1))/3.d0
         elseif(inod.eq.3) then
            rhol(inod)=(rho(nr)+2.d0*rho(nr+1))/3.d0
         else
            rhol(inod)=rho(nr+1)
         endif
      enddo
            
      do inod=1,4
         rho0=rhol(inod)
         do nfc2=1,nfcmax
            nth=nthfnfc(nfc2)
            nph=nphfnfc(nfc2)
                  
            rkth=nth/(ra*rho0)
            rkth0=1.d0/(ra*rho0)
            rkph=nph/rr

               do nfc1=1,nfcmax
                  fmc(1,1,1,nfc1,nfc2,inod)
     &                 = rho0*(cfactor-rkph**2-rkth**2)
                  fmc(1,2,1,nfc1,nfc2,inod)
     &                 = rho0*(-ci*rkth*rkth0)
                  fmc(2,1,1,nfc1,nfc2,inod)
     &                 = rho0*(+ci*rkth*rkth0)
                  fmc(2,2,1,nfc1,nfc2,inod)
     &                 = rho0*(cfactor-rkph**2-rkth0**2)
                  fmc(2,3,1,nfc1,nfc2,inod)
     &                 = rho0*(rkth*rkph)
                  fmc(3,2,1,nfc1,nfc2,inod)
     &                 = rho0*(rkth*rkph)
                  fmc(3,3,1,nfc1,nfc2,inod)
     &                 = rho0*(cfactor-rkth**2)
                  
                  fmc(2,1,2,nfc1,nfc2,inod)
     &                 = rho0*( ci*rkth)/ra
                  fmc(2,2,2,nfc1,nfc2,inod)
     &                 = rho0*(  -rkth0)/ra
                  fmc(3,1,2,nfc1,nfc2,inod)
     &                 = rho0*( ci*rkph)/ra
                  
                  fmc(1,2,3,nfc1,nfc2,inod)
     &                 = rho0*(-ci*rkth)/ra
                  fmc(1,3,3,nfc1,nfc2,inod)
     &                 = rho0*(-ci*rkph)/ra
                  fmc(2,2,3,nfc1,nfc2,inod)
     &                 = rho0*(  -rkth0)/ra
                  
                  fmc(1,1,4,nfc1,nfc2,inod)
     &                 = 0.d0
                  fmc(2,2,4,nfc1,nfc2,inod)
     &                 = rho0*(-1.d0)/ra**2
                  fmc(3,3,4,nfc1,nfc2,inod)
     &                 = rho0*(-1.d0)/ra**2
            enddo
         enddo
      enddo

      do inod=1,4               ! Fourier transform
         do k=1,4
            do j=1,3
               do i=1,3
                  do nfc2=1,nfcmax
                     do nfc1=1,nfcmax
                        nth=nthnfc(nfc1)
                        nph=nphnfc(nfc1)
                        fv1(nth,nph)=fmc(i,j,k,nfc1,nfc2,inod)
                     enddo
                     call wmsubfx(fv1,fv1f,nthmax,nphmax)
                     do nfc1=1,nfcmax
                        nth=nthnfc(nfc1)
                        nph=nphnfc(nfc1)
                        fmc(i,j,k,nfc1,nfc2,inod)=fv1f(nth,nph)
                     enddo
                  enddo
               enddo    
            enddo
         enddo    
      enddo

      do nph1=1,nphmax          ! Fit to fmd and adjust m and n
         do nth1=1,nthmax
            nfc1=nthmax*(nph1-1)+nth1
            do nph2=1,nphmax
               do nth2=1,nthmax
                  nfc2=nthmax*(nph2-1)+nth2
                  if(nphmax.eq.1) then
                     nphdiff=0
                  else
                     nphdiff=nph1-nph2
                  endif
                  if(nthmax.eq.1) then
                     nthdiff=0
                  else
                     nthdiff=nth1-nth2
                  endif
                  if(nfcmax.eq.1) then
                     nfcdiff=1
                  else
                     nfcdiff=nthmax*nphdiff+nthdiff+nfcmax/2
                  endif
                  if(nfcdiff.ge.1.and.nfcdiff.le.nfcmax) then
                     do inod=1,4
                        do k=1,4
                           do j=1,3
                              do i=1,3
                                 fmd(i,j,k,nfc1,nfc2,inod)
     &                                =fmc(i,j,k,nfcdiff,nfc2,inod)
                              enddo
                           enddo
                        enddo
                     enddo
                  else
                     do inod=1,4
                        do k=1,4
                           do j=1,3
                              do i=1,3
                                 fmd(i,j,k,nfc1,nfc2,inod)=0.d0
                              enddo
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_gc

!----- calculate coefficint matrix fmd (E cylindrical +,-,para)-----

      subroutine wmfem_calculate_vacuum_lc(nr,fmd)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmc
      complex(8),dimension(nthmax,nphmax):: fv1,fv1f
      real(8),dimension(4):: rhol
      complex(8):: cfactor
      real(8):: drho,rkth,rkph,rkth0,rho0
      integer:: ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,nthx
      integer:: ns,nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      integer:: nph1,nph2,nph1x,nph2x
      integer:: nphdiff,nphxdiff
      integer:: nth1,nth2,nth1x,nth2x
      integer:: nthdiff,nthxdiff
      integer:: nfcdiff
      real(8):: rr,ra,rb

      call get_wmparm1(rr,ra,rb)
      cfactor=(2*pi*crf*1.d6)**2/vc**2

      do inod=1,4               ! clear fmc
         do nfc2=1,nfcmax
            do nfc1=1,nfcmax
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmc(i,j,k,nfc1,nfc2,inod)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do inod=1,4               ! define local radius for different inod
         if(inod.eq.1) then
            if(nr.eq.1) then
               rhol(inod)=(8.d0*rho(nr)+rho(nr+1))/9.d0
            else
               rhol(inod)=rho(nr)
            endif
         elseif(inod.eq.2) then
            rhol(inod)=(2.d0*rho(nr)+rho(nr+1))/3.d0
         elseif(inod.eq.3) then
            rhol(inod)=(rho(nr)+2.d0*rho(nr+1))/3.d0
         else
            rhol(inod)=rho(nr+1)
         endif
      enddo
            
      do inod=1,4
         rho0=rhol(inod)
         do nfc2=1,nfcmax
            nth=nthnfc(nfc2)
            nph=nphnfc(nfc2)
                  
            rkth=nth/(ra*rho0)
            rkth0=1.d0/(ra*rho0)
            rkph=nph/rr

               do nfc1=1,nfcmax
                  fmc(1,1,1,nfc1,nfc2,inod)
     &                 = rho0*(cfactor-rkph**2-rkth**2)
                  fmc(1,2,1,nfc1,nfc2,inod)
     &                 = rho0*(-ci*rkth*rkth0)
                  fmc(2,1,1,nfc1,nfc2,inod)
     &                 = rho0*(+ci*rkth*rkth0)
                  fmc(2,2,1,nfc1,nfc2,inod)
     &                 = rho0*(cfactor-rkph**2-rkth0**2)
                  fmc(2,3,1,nfc1,nfc2,inod)
     &                 = rho0*(rkth*rkph)
                  fmc(3,2,1,nfc1,nfc2,inod)
     &                 = rho0*(rkth*rkph)
                  fmc(3,3,1,nfc1,nfc2,inod)
     &                 = rho0*(cfactor-rkth**2)
                  
                  fmc(2,1,2,nfc1,nfc2,inod)
     &                 = rho0*( ci*rkth)/ra
                  fmc(2,2,2,nfc1,nfc2,inod)
     &                 = rho0*(  -rkth0)/ra
                  fmc(3,1,2,nfc1,nfc2,inod)
     &                 = rho0*( ci*rkph)/ra
                  
                  fmc(1,2,3,nfc1,nfc2,inod)
     &                 = rho0*(-ci*rkth)/ra
                  fmc(1,3,3,nfc1,nfc2,inod)
     &                 = rho0*(-ci*rkph)/ra
                  fmc(2,2,3,nfc1,nfc2,inod)
     &                 = rho0*(  -rkth0)/ra
                  
                  fmc(1,1,4,nfc1,nfc2,inod)
     &                 = 0.d0
                  fmc(2,2,4,nfc1,nfc2,inod)
     &                 = rho0*(-1.d0)/ra**2
                  fmc(3,3,4,nfc1,nfc2,inod)
     &                 = rho0*(-1.d0)/ra**2
            enddo
         enddo
      enddo

      do inod=1,4               ! Fourier transform
         do k=1,4
            do j=1,3
               do i=1,3
                  do nfc2=1,nfcmax
                     do nph=1,nphmax
                        do nth=1,nthmax
                           nfc1=nthmax*(nph-1)+nth
                           fv1(nth,nph)=fmc(i,j,k,nfc1,nfc2,inod)
                        enddo
                     enddo
                     call wmsubfx(fv1,fv1f,nthmax,nphmax)
                     do nph=1,nphmax
                        do nth=1,nthmax
                           nfc1=nthmax*(nph-1)+nth
                           fmc(i,j,k,nfc1,nfc2,inod)=fv1f(nth,nph)
                        enddo
                     enddo
                  enddo
               enddo    
            enddo
         enddo    
      enddo

      do nph1=1,nphmax          ! Fit to fmd and adjust m and n
         do nth1=1,nthmax
            nfc1=nthmax*(nph1-1)+nth1
            do nph2=1,nphmax
               do nth2=1,nthmax
                  nfc2=nthmax*(nph2-1)+nth2
                  if(nphmax.eq.1) then
                     nphdiff=0
                  else
                     nphdiff=nph1-nph2
                  endif
                  if(nthmax.eq.1) then
                     nthdiff=0
                  else
                     nthdiff=nth1-nth2
                  endif
                  if(nfcmax.eq.1) then
                     nfcdiff=1
                  else
                     nfcdiff=nthmax*nphdiff+nthdiff+nfcmax/2
                  endif
                  if(nfcdiff.ge.1.and.nfcdiff.le.nfcmax) then
                     do inod=1,4
                        do k=1,4
                           do j=1,3
                              do i=1,3
                                 fmd(i,j,k,nfc1,nfc2,inod)
     &                                =fmc(i,j,k,nfcdiff,nfc2,inod)
                              enddo
                           enddo
                        enddo
                     enddo
                  else
                     do inod=1,4
                        do k=1,4
                           do j=1,3
                              do i=1,3
                                 fmd(i,j,k,nfc1,nfc2,inod)=0.d0
                              enddo
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_lc

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma(nr,ns,fmd)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmc
      complex(8),dimension(3,3,nfcmax,nfcmax):: fms1,fms2
      complex(8),dimension(nthmax,nphmax):: fv1,fv1f
      real(8),dimension(nthmax,nphmax,4):: gjl
      real(8),dimension(4):: rhol
      real(8):: drho,rkth,rkph,rkth0,rho0
      integer:: ml,mw,mc,nvmax,i,j,k,inod,nfc,nf1,nf2,nth,nph,nthx
      integer:: ns,nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      complex(8):: cfactor
      integer:: nph1,nph2,nph1x,nph2x
      integer:: nphdiff,nphxdiff
      integer:: nth1,nth2,nth1x,nth2x
      integer:: nthdiff,nthxdiff
      integer:: nfcdiff

      cfactor=(2*pi*crf*1.d6)**2/vc**2

      do inod=1,4               ! clear fmc
         do nf2=1,nfcmax
            do nf1=1,nfcmax
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmc(i,j,k,nf1,nf2,inod)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do inod=1,4               ! define local radius for different inod
         if(nr.eq.1) then
            do nph=1,nphmax
            do nth=1,nthmax
               gjl(nth,nph,1)=(8.d0*gj(nth,nph,nr)
     &                             +gj(nth,nph,nr+1))/9.d0
            enddo
            enddo
         else
            do nph=1,nphmax
            do nth=1,nthmax
               gjl(nth,nph,1)=gj(nth,nph,nr)
            enddo
            enddo
         endif
         do nph=1,nphmax
         do nth=1,nthmax
            gjl(nth,nph,2)=(2.d0*gj(nth,nph,nr)
     &                          +gj(nth,nph,nr+1))/3.d0
            gjl(nth,nph,3)=(     gj(nth,nph,nr)
     &                     +2.d0*gj(nth,nph,nr+1))/3.d0
            gjl(nth,nph,4)=      gj(nth,nph,nr+1)
         enddo
         enddo
      enddo
            
      call wmfem_disp(nr  ,ns,fms1)
      call wmfem_disp(nr+1,ns,fms2)

      do nfc2=1,nfcmax
         do nfc1=1,nfcmax
            do j=1,3
               do i=1,3
                  fmc(i,j,1,nfc1,nfc2,1)=   fms1(i,j,nfc1,nfc2)
                  fmc(i,j,1,nfc1,nfc2,2)=(2*fms1(i,j,nfc1,nfc2)
     &                                   +  fms2(i,j,nfc1,nfc2))/3.d0
                  fmc(i,j,1,nfc1,nfc2,3)=(  fms1(i,j,nfc1,nfc2)
     &                                   +2*fms2(i,j,nfc1,nfc2))/3.d0
                  fmc(i,j,1,nfc1,nfc2,4)=   fms2(i,j,nfc1,nfc2)
                  do inod=1,4
                     fmc(i,j,1,nfc1,nfc2,inod)
     &                    =cfactor*fmc(i,j,1,nfc1,nfc2,inod)
                     do k=2,4
                        fmc(i,j,k,nfc1,nfc2,inod)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do inod=1,4               ! Fourier transform
         do k=1,4
            do j=1,3
               do i=1,3
                  do nfc2=1,nfcmax
                     do nfc1=1,nfcmax
                        nth=nthnfc(nfc1)
                        nph=nphnfc(nfc1)
                        fv1(nth,nph)=fmc(i,j,k,nfc1,nfc2,inod)
     &                              *gjl(nth,nph,inod)
                     enddo
                     call wmsubfx(fv1,fv1f,nthmax,nphmax)
                     do nfc1=1,nfcmax
                        nth=nthnfc(nfc1)
                        nph=nphnfc(nfc1)
                        fmc(i,j,k,nfc1,nfc2,inod)=fv1f(nth,nph)
                     enddo
                  enddo
               enddo    
            enddo
         enddo    
      enddo

      do nph1=1,nphmax          ! Fit to fmd and adjust m and n
         do nth1=1,nthmax
            nfc1=nthmax*(nph1-1)+nth1
            do nph2=1,nphmax
               do nth2=1,nthmax
                  nfc2=nthmax*(nph2-1)+nth2
                  if(nphmax.eq.1) then
                     nphdiff=0
                  else
                     nphdiff=nph1-nph2
                  endif
                  if(nthmax.eq.1) then
                     nthdiff=0
                  else
                     nthdiff=nth1-nth2
                  endif
                  if(nfcmax.eq.1) then
                     nfcdiff=1
                  else
                     nfcdiff=nthmax*nphdiff+nthdiff+nfcmax/2
                  endif
                  if(nfcdiff.ge.1.and.nfcdiff.le.nfcmax) then
                     do inod=1,4
                        do k=1,4
                           do j=1,3
                              do i=1,3
                                 fmd(i,j,k,nfc1,nfc2,inod)
     &                                =fmc(i,j,k,nfcdiff,nfc2,inod)
                              enddo
                           enddo
                        enddo
                     enddo
                  else
                     do inod=1,4
                        do k=1,4
                           do j=1,3
                              do i=1,3
                                 fmd(i,j,k,nfc1,nfc2,inod)=0.d0
                              enddo
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo
      return

      end subroutine wmfem_calculate_plasma

!---- FEM cubic hermit ---

      subroutine fem_hhh(nr,fmd,nfcmax,drho)

      use libfem_mod
      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(in):: fmd
      integer,intent(in):: nfcmax
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nf1,nf2,inod,ml,mw

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif

      mr=6*nfcmax  ! line interval between radial points 
      mc=(mwmax+1)/2


         do nf1=1,nfcmax
         do nf2=1,nfcmax
         do j=1,3
         do i=1,3
            ml=6*nfcmax*(nr-1)+6*(nf1-1)+2*(i-1)+1
            mw=mc+6*(nf2-nf1)+2*(j-i)
            do inod=1,4

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+Mr+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,6)*drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,7)
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,8)*drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,8)

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,5)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+mr+1)=fma(mw  ,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,4)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,4)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,8)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,8)*drho

            enddo
         enddo
         enddo
         enddo
         enddo

      return
      end subroutine fem_hhh

!---- FEM cubic hermit + quadratic ---

      subroutine fem_hqq(nr,fmd,nfcmax,drho)

      use libfem_mod
      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(in):: fmd
      integer,intent(in):: nfcmax
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nf1,nf2,inod,ml,mw

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif

c$$$      do i=1,3
c$$$      do j=1,3
c$$$      do k=1,4
c$$$      do nf1=1,nfcmax
c$$$      do nf2=1,nfcmax
c$$$      do inod=1,4
c$$$         write(16,'(7I5,A,1P2E12.4)') 
c$$$     &        nr,i,j,k,nf1,nf2,inod,':',fmd(i,j,k,nf1,nf2,inod)
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      enddo

      mr=6*nfcmax  ! line interval between radial points 
      mc=(mwmax+1)/2

         do nf1=1,nfcmax
         do nf2=1,nfcmax
         do j=1,3
         do i=1,3
            ml=6*nfcmax*(nr-1)+6*(nf1-1)+2*(i-1)+1
            mw=mc+6*(nf2-nf1)+2*(j-i)
            do inod=1,4

! (1,1) **********
            if(i.eq.1.and.j.eq.1) then

            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,4,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,1,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,4,4)/drho

            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,1,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,4,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,4,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,4,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,1,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,4,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,2,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,2,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,5,4)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,2,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,5,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,2,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,5,5)/drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,2,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,2,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,5,6)/drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,6,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,3,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,6,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,3,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,6,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,6,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,6,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,3,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,6,6)/drho

! (1,*) **********
            else if(i.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,5)/drho
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,6)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,7)/drho
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,8)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,8)

! (*,1) **********
            else if(j.eq.1) then
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,1,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,5,4)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,1,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,5,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,5,5)/drho
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,1,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,5,6)/drho

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,2,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,6,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,2,4)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,6,4)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,2,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,6,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,2,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,6,5)
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,2,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,6,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,2,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,6,6)

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,7,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,3,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,7,4)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,3,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,7,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,7,5)/drho
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,7,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,3,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,7,6)/drho

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,4,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,8,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,4,4)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,8,4)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,4,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,8,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,4,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,8,5)
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,4,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,8,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,4,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,8,6)

! (*,*) **********
            else
            fma(mw  ,ml  )=fma(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,5)/drho
            fma(mw+1,ml  )=fma(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,6)
            fma(mw+mr,ml  )=fma(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,7)/drho
            fma(mw+mr+1,ml  )=fma(mw+mr+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,8)

            fma(mw-1,ml+1)=fma(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,5)
            fma(mw  ,ml+1)=fma(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,6)*drho
            fma(mw+mr-1,ml+1)=fma(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,7)
            fma(mw+mr,ml+1)=fma(mw+mr,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,8)*drho

            fma(mw-mr,ml+mr)=fma(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,5)/drho
            fma(mw-mr+1,ml+mr)=fma(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,6)
            fma(mw  ,ml+mr)=fma(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,7)/drho
            fma(mw+1,ml+mr)=fma(mw+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,8)

            fma(mw-mr-1,ml+mr+1)=fma(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,5)
            fma(mw-mr,ml+mr+1)=fma(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,6)*drho
            fma(mw-1,ml+mr+1)=fma(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,7)
            fma(mw  ,ml+mr+1)=fma(mw  ,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,4)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,4)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,8)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,8)*drho
            endif

            enddo
         enddo
         enddo
         enddo
         enddo
      return
      end subroutine fem_hqq

      end subroutine wmfem
