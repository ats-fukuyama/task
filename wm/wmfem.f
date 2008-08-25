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

!      real(8),dimension(:,:,:,:,:),allocatable :: gma,mma 
!     &                                          !(3,3,nthmax,nphmax,nrmax)
!      real(8),dimension(:,:,:),allocatable:: gj !(nthmax,nphmax,nrmax)

      integer,dimension(:),allocatable,save :: nthnfc,nthfnfc
      integer,dimension(:),allocatable,save :: nphnfc,nphfnfc
      integer,dimension(:),allocatable,save :: nthnfc2,nthfnfc2
      integer,dimension(:),allocatable,save :: nphnfc2,nphfnfc2

!     ***** define array size  *****

      nfcmax=nthmax*nphmax      ! size of block matrix 
                                !    (number of Fourier components)
      mlmax=6*nfcmax*nrmax      ! length of coeffient matrix and source vector
      mwmax=4*6*nfcmax-1        ! width of coefficient matrix
      nthmax2=2*nthmax          ! number of poloidal modes of coefficients
      nphmax2=2*nphmax          ! number of toroidal modes of coefficients
      nfcmax2=nthmax2*nphmax2   ! size of block matrix of coefficients

!     ***** get additional parameters *****

      call get_wmparm(crf,nth0,nph0,idbgwm)

!     ***** allocate matrix and vector *****

      call wmfem_allocate

      if(nrmax.eq.0) return   ! matrix and vector deallocated

!     ***** setup Fourier component index *****

      call wmfem_setup_index

!     ***** calculate metric and convaersion tensor *****

!      call wmfem_metric(gma,mma,gj)

!     ***** calculate coefficient matrix *****

      call wmfem_calculate

!     ***** solve matrix equation *****

      call wmfem_solve

      return

      contains

!     ---------------------------
!     ----- allocate arrays -----

      subroutine wmfem_allocate

      implicit none
!      integer,save:: nrmax_save=0,nthmax_save=0,nphmax_save=0
      integer,save:: mwmax_save=0,mlmax_save=0
      integer,save:: nsmax_save=0,nfcmax_save=0

!      if((nrmax.ne.nrmax_save).or.(nthmax.ne.nthmax_save)
!     &                        .or.(nphmax.ne.nphmax_save)) then
!         if(allocated(gma)) deallocate(gma)
!         if(allocated(mma)) deallocate(mma)
!         if(allocated(gj)) deallocate(gj)
!         allocate(gma(3,3,nthmax2,nphmax2,nrmax))
!         allocate(mma(3,3,nthmax2,nphmax2,nrmax))
!         allocate(gj(nthmax2,nphmax2,nrmax))
!      endif

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

!     ---------------------------------------
!     ---- setup Fourier component index ----

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

!     -------------------------------=
!     ---- Solve matrix equation  ----

      subroutine wmfem_solve

      integer:: mc,ml,mw,ierr,nr,nth,nph,ns,nr1,nth1,nph1,i,j
      complex(8):: csum,csums,csum1,csum2

      mc=(mwmax+1)/2    ! diagonal position in mw

!     ----- solve matrix -----

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

!     ----- calculate E field -----

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

!     ----- calculate power -----

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
     &                      *fms(mw+j-i,ml+i,ns)
     &                           *fvx(ml+mw-mc+j)
               enddo
               enddo
               endif
            enddo
            csums=csums+csum
            cpp(nth,nph,nth1,nph1,nr,ns)=csum
!            write(6,'(6I5,1P2E12.4)') nth,nph,nth1,nph1,nr,ns,
!     &           cpp(nth,nph,nth1,nph1,nr,ns)
         enddo
         enddo
      enddo
      enddo
      enddo
!         write(6,'(A,I5,1P2E12.4)') 'ns,csums=',ns,csums
      enddo

!     ----- calculate antenna impedance -----

      do nph=1,nphmax
      do nth=1,nthmax
         cpa(nth,nph)=0.d0
         do nr=1,nrmax-1
            ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph-1)+6*(nth-1)
            do i=1,6
               cpa(nth,nph)=cpa(nth,nph)+conjg(fvx(ml+i))*fvb(ml+i)
            enddo
         enddo
!         write(6,'(2I5,1P2E12.4)') nth,nph,cpa(nth,nph)
      enddo
      enddo

      return
      end subroutine wmfem_solve

!     -------------------------------------------
!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate

      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: 
     &     fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax):: 
     &     fmd1,fmd2,fmd3,fmd4
      complex(8),dimension(nphmax,nthmax,3):: fvb_nr
      complex(8):: cfactor
      real(8):: drho,rkth,rkph,rkth0,rho0,rho1,rho2,rho3,rho4
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,nthf,nphf
      integer:: ns,nfc1,nfc2
      complex(8):: csum,f1,f2,f3,f4
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
            rho1=rho(nr)
            rho2=(2.d0*rho(nr)+rho(nr+1))/3.d0
            rho3=(rho(nr)+2.d0*rho(nr+1))/3.d0
            rho4=rho(nr+1)

            if(ns.eq.0) then
               if(idbgwm.eq.0) then
                  call wmfem_calculate_vacuum_g(rho1,fmd1)
                  call wmfem_calculate_vacuum_g(rho2,fmd2)
                  call wmfem_calculate_vacuum_g(rho3,fmd3)
                  call wmfem_calculate_vacuum_g(rho4,fmd4)
               else if(idbgwm.eq.1) then
                  call wmfem_calculate_vacuum_c(rho1,fmd1)
                  call wmfem_calculate_vacuum_c(rho2,fmd2)
                  call wmfem_calculate_vacuum_c(rho3,fmd3)
                  call wmfem_calculate_vacuum_c(rho4,fmd4)
               endif
            else
               call wmfem_calculate_plasma(rho1,ns,fmd1)
               call wmfem_calculate_plasma(rho2,ns,fmd2)
               call wmfem_calculate_plasma(rho3,ns,fmd3)
               call wmfem_calculate_plasma(rho4,ns,fmd4)
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
                           f1=fmd1(i,j,k,nfc1,nfc2)
                           f2=fmd2(i,j,k,nfc1,nfc2)
                           f3=fmd3(i,j,k,nfc1,nfc2)
                           f4=fmd4(i,j,k,nfc1,nfc2)
                           fmd(i,j,k,nfc1,nfc2,1)=f1
                           fmd(i,j,k,nfc1,nfc2,2)
     &                       =0.5d0*(-11*f1+18*f2- 9*f3+ 2*f4)
                           fmd(i,j,k,nfc1,nfc2,3)=f4
                           fmd(i,j,k,nfc1,nfc2,4)
     &                       =0.5d0*(- 2*f1+ 9*f2-18*f3+11*f4)
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

      subroutine wmfem_calculate_vacuum_g(rhol,fmd)

      implicit none
      real(8),intent(in):: rhol
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8),dimension(3,3,3,3,nfcmax2):: fmv1
      complex(8),dimension(3,3,3,nfcmax2):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2):: fmv4
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,nth1,nth2,nthdiff,nph,nph1,nph2,nphdiff
      integer:: imn1,imn2,nfcdiff

      call wmfem_calculate_vacuum_sub(rhol,fmv1,fmv2,fmv3,fmv4)

      do i=1,3
         do j=1,3
            do imn1=1,3
               do imn2=1,3
                  do nph=1,nphmax2
                     do nth=1,nthmax2
                        nfc=nthmax2*(nph-1)+nth
                        fv1(nth,nph)=fmv1(i,j,imn1,imn2,nfc)
                     enddo
                  enddo
                  call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
                  do nph=1,nphmax2
                     do nth=1,nthmax2
                        nfc=nthmax2*(nph-1)+nth
                        fmv1(i,j,imn1,imn2,nfc)=fv1f(nth,nph)
                     enddo
                  enddo
               enddo
            enddo    
            do imn1=1,3
               do nph=1,nphmax2
                  do nth=1,nthmax2
                     nfc=nthmax2*(nph-1)+nth
                     fv1(nth,nph)=fmv2(i,j,imn1,nfc)
                  enddo
               enddo
               call wmsubfx(fv1,fv1f,nthmax,nphmax)
               do nph=1,nphmax2
                  do nth=1,nthmax2
                      nfc=nthmax2*(nph-1)+nth
                     fmv2(i,j,imn1,nfc)=fv1(nth,nph)
                  enddo
               enddo
            enddo
            do imn1=1,3
               do nph=1,nphmax2
                  do nth=1,nthmax2
                     nfc=nthmax2*(nph-1)+nth
                     fv1(nth,nph)=fmv3(i,j,imn1,nfc)
                  enddo
               enddo
               call wmsubfx(fv1,fv1f,nthmax,nphmax)
               do nph=1,nphmax2
                  do nth=1,nthmax2
                     nfc=nthmax2*(nph-1)+nth
                     fmv3(i,j,imn1,nfc)=fv1(nth,nph)
                  enddo
               enddo
            enddo
            do nph=1,nphmax2
               do nth=1,nthmax2
                  nfc=nthmax2*(nph-1)+nth
                  fv1(nth,nph)=fmv4(i,j,nfc)
               enddo
            enddo
            call wmsubfx(fv1,fv1f,nthmax,nphmax)
            do nph=1,nphmax2
               do nth=1,nthmax2
                  nfc=nthmax2*(nph-1)+nth
                  fmv4(i,j,nfc)=fv1(nth,nph)
               enddo
            enddo
         enddo
      enddo

      do nfc1=1,nfcmax          ! Fit to fmd and adjust m and n
         nph1=nphfnfc(nfc1)
         nth1=nthfnfc(nfc1)
         do nfc2=1,nfcmax
            nph2=nphfnfc(nfc2)
            nth2=nthfnfc(nfc2)

            nphdiff=nph1-nph2
            nthdiff=nth1-nth2
            nfcdiff=nthmax*nphdiff+nthdiff+nfcmax

            do j=1,3
               do i=1,3
                  fmd(i,j,1,nfc1,nfc2)
     &              =fmv1(i,j,1,1,nfcdiff)
     &              +fmv1(i,j,2,1,nfcdiff)*nth1
     &              +fmv1(i,j,3,1,nfcdiff)*nph1
     &              +fmv1(i,j,1,2,nfcdiff)     *nth2
     &              +fmv1(i,j,2,2,nfcdiff)*nth1*nth2
     &              +fmv1(i,j,3,2,nfcdiff)*nph1*nth2
     &              +fmv1(i,j,1,3,nfcdiff)     *nph2
     &              +fmv1(i,j,2,3,nfcdiff)*nth1*nph2
     &              +fmv1(i,j,3,3,nfcdiff)*nph1*nph2
               fmd(i,j,2,nfc1,nfc2)
     &              =fmv2(i,j,1,nfcdiff)
     &              +fmv2(i,j,2,nfcdiff)*nth1
     &              +fmv2(i,j,3,nfcdiff)*nph1
               fmd(i,j,3,nfc1,nfc2)
     &              =fmv3(i,j,1,nfcdiff)
     &              +fmv3(i,j,2,nfcdiff)     *nth2
     &              +fmv3(i,j,3,nfcdiff)     *nph2
               fmd(i,j,4,nfc1,nfc2)
     &              =fmv4(i,j,nfcdiff)
            enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_g

!----- calculate vacuum matrix fms -----

      subroutine wmfem_calculate_vacuum_sub(rhol,fmv1,fmv2,fmv3,fmv4)

      implicit none
      real(8),intent(in):: rhol
      complex(8),dimension(3,3,3,3,nfcmax2),intent(out):: fmv1
      complex(8),dimension(3,3,3,nfcmax2),intent(out):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2),intent(out):: fmv4

      real(8),dimension(:,:,:,:),allocatable :: gma,mma,dmma 
     &                                                !(3,3,nthmax,nphmax)
      real(8),dimension(:,:),allocatable:: gj     !(nthmax,nphmax)

      complex(8),dimension(3,3,3):: cq
      complex(8),dimension(3,3):: cp
      integer:: i,j,k,l,nthm,nthp,nphm,nphp
      integer:: imn1,imn2
      integer:: nph,nth
      real(8):: dph,dth
      complex(8):: csum1,csum2,csum3,csum4,cfactor
      integer:: nrl,nfc2
      real(8):: drhob,mma3b,mma2b,drhoa,mma3a,mma2a,mma30,mma20
      real(8):: mma3d,mma2d,gjl
      
      cfactor=(2*pi*crf*1.d6)**2/vc**2

      call wmfem_metrics(rhol,nthmax2,nphmax2,gma,mma,dmma,gj)

      do nph=1,nphmax2
         if(nph.eq.1) then
            nphm=nphmax2
         else
            nphm=nph-1
         endif
         if(nph.eq.nphmax2) then
            nphp=1
         else
            nphp=nph+1
         endif
         dph=2*pi/nphmax2

         do nth=1,nthmax2
            if(nth.eq.1) then
               nthm=nthmax2
            else
               nthm=nth-1
            endif
            if(nth.eq.nthmax2) then
               nthp=1
            else
               nthp=nth+1
            endif
            dth=2*pi/nthmax2
            gjl=gj(nth,nph)
            nfc2=nthmax2*(nph-1)+nth

            do j=1,3
               cq(1,j,1)=((mma(3,j,nthp,nph)
     &                    -mma(3,j,nthm,nph))/dth
     &                   -(mma(2,j,nth,nphp)
     &                    -mma(2,j,nth,nphm))/dph )/gjl
               cq(1,j,2)=+ci*mma(3,j,nth,nph)/gjl
               cq(1,j,3)=-ci*mma(2,j,nth,nph)/gjl

               cq(2,j,1)=((mma(1,j,nth,nphp)
     &                    -mma(1,j,nth,nphm))/dph
     &                    -dmma(3,j,nth,nph))/gjl
               cq(2,j,2)=0.d0
               cq(2,j,3)=+ci*mma(1,j,nth,nph)/gjl

               cq(3,j,1)=(dmma(2,j,nth,nph)
     &                  -(mma(1,j,nthp,nph)
     &                   -mma(1,j,nthm,nph))/dth )/gjl
               cq(3,j,2)=-ci*mma(1,j,nth,nph)/gjl
               cq(3,j,3)=0.d0

               cp(1,j)=0.d0
               cp(2,j)=-mma(3,j,nth,nph)/gjl
               cp(3,j)= mma(2,j,nth,nph)/gjl
            enddo

c$$$      write(6,*) 'cq(1)'
c$$$      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
c$$$     &     (cq(i,1,1),cq(i,2,1),cq(i,3,1),i=1,3)
c$$$      write(6,*) 'cq(2)'
c$$$      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
c$$$     &     (cq(i,1,2),cq(i,2,2),cq(i,3,2),i=1,3)
c$$$      write(6,*) 'cq(3)'
c$$$      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
c$$$     &     (cq(i,1,3),cq(i,2,3),cq(i,3,3),i=1,3)
c$$$      write(6,*) 'cp'
c$$$      write(6,'(1P2E12.4,2X,1P2E12.3,2X,1P2E12.3)') 
c$$$     &     (cp(i,1),  cp(i,2),  cp(i,3),i=1,3)

            do imn2=1,3
               do imn1=1,3
                  do j=1,3
                     do i=1,3
                        csum1=0.d0
                        do k=1,3
                           do l=1,3
                              csum1=csum1+conjg(cq(k,i,imn1))
     &                                   *gma(k,l,nth,nph)
     &                                   *cq(l,j,imn2) *gjl
                           enddo
                        enddo
                        if(i.eq.j.and.imn1.eq.1.and.imn2.eq.1) then
                           fmv1(i,j,imn1,imn2,nfc2)
     &                          =csum1-cfactor*gj(nth,nph)
                        else
                           fmv1(i,j,imn1,imn2,nfc2)=csum1
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
     &                                *gma(k,l,nth,nph)
     &                                *cq(l,j,imn1)*gjl
                           csum3=csum3+conjg(cq(k,i,imn1))
     &                                *gma(k,l,nth,nph)
     &                                *cp(l,j)*gjl
                        enddo
                     enddo
                     fmv2(i,j,imn1,nfc2)=csum2
                     fmv3(i,j,imn1,nfc2)=csum3
                  enddo
               enddo
            enddo
            
            do j=1,3
               do i=1,3
                  csum4=0.d0
                  do k=1,3
                     do l=1,3
                        csum4=csum4+conjg(cp(k,i))
     &                       *gma(k,l,nth,nph)
     &                       *cp(l,j) *gjl
                     enddo
                  enddo
                  fmv4(i,j,nfc2)=csum4
               enddo
            enddo
            
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_sub

!----- calculate coefficint matrix fmd (E cylindrical +,-,para)-----

      subroutine wmfem_calculate_vacuum_c(rhol,fmd)

      implicit none
      real(8),intent(in):: rhol
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8):: cfactor
      real(8):: drho,rkth,rkph,rkth0,rho0
      integer:: ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph
      integer:: ns,nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      integer:: nph1,nph2,nth1,nth2
      integer:: nphdiff,nthdiff,nfcdiff
      real(8):: rr,ra,rb

      call get_wmparm1(rr,ra,rb)
      cfactor=(2*pi*crf*1.d6)**2/vc**2

!     ----- clear fmc -----
      do nfc2=1,nfcmax
         do nfc1=1,nfcmax2
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,nfc1,nfc2)=0.d0
                  enddo
               enddo
            enddo
         enddo
      enddo

      do nfc2=1,nfcmax
         nth=nthnfc(nfc2)
         nph=nphnfc(nfc2)

         if(rhol.ne.0.d0) then
            rkth=nth/(ra*rhol)
            rkth0=1.d0/(ra*rhol)
         else
            rkth=0.d0
            rkth0=0.d0
         endif
         rkph=nph/rr

         do nfc1=1,nfcmax2
            fmc(1,1,1,nfc1,nfc2)= rhol*(cfactor-rkph**2-rkth**2)
            fmc(1,2,1,nfc1,nfc2)= rhol*(-ci*rkth*rkth0)
            fmc(2,1,1,nfc1,nfc2)= rhol*(+ci*rkth*rkth0)
            fmc(2,2,1,nfc1,nfc2)= rhol*(cfactor-rkph**2-rkth0**2)
            fmc(2,3,1,nfc1,nfc2)= rhol*(rkth*rkph)
            fmc(3,2,1,nfc1,nfc2)= rhol*(rkth*rkph)
            fmc(3,3,1,nfc1,nfc2)= rhol*(cfactor-rkth**2)
                  
            fmc(2,1,2,nfc1,nfc2)= rhol*( ci*rkth)/ra
            fmc(2,2,2,nfc1,nfc2)= rhol*(  -rkth0)/ra
            fmc(3,1,2,nfc1,nfc2)= rhol*( ci*rkph)/ra
                  
            fmc(1,2,3,nfc1,nfc2)= rhol*(-ci*rkth)/ra
            fmc(1,3,3,nfc1,nfc2)= rhol*(-ci*rkph)/ra
            fmc(2,2,3,nfc1,nfc2)= rhol*(  -rkth0)/ra
                  
            fmc(1,1,4,nfc1,nfc2)= 0.d0
            fmc(2,2,4,nfc1,nfc2)= rhol*(-1.d0)/ra**2
            fmc(3,3,4,nfc1,nfc2)= rhol*(-1.d0)/ra**2
         enddo
      enddo

      do k=1,4
         do j=1,3
            do i=1,3
               do nfc2=1,nfcmax
                  do nfc1=1,nfcmax2
                     nth=nthnfc2(nfc1)
                     nph=nphnfc2(nfc1)
                     fv1(nth,nph)=fmc(i,j,k,nfc1,nfc2)
                  enddo
                  call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
                  do nfc1=1,nfcmax2
                     nth=nthfnfc2(nfc1)
                     nph=nphfnfc2(nfc1)
                     fmc(i,j,k,nfc1,nfc2)=fv1f(nth,nph)
                  enddo
               enddo
            enddo    
         enddo
      enddo    

!     ----- Fit to fmd and adjust m and n -----

      do nfc2=1,nfcmax
         nth2=nthfnfc(nfc2)
         nph2=nphfnfc(nfc2)
         do nfc1=1,nfcmax
            nth1=nthfnfc(nfc1)
            nph1=nphfnfc(nfc1)

            nphdiff=nph1-nph2
            nthdiff=nth1-nth2
            nfcdiff=nthmax2*nphdiff+nthdiff+nfcmax

            if(nfcdiff.ge.1.and.nfcdiff.le.nfcmax2) then
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmd(i,j,k,nfc1,nfc2)=fmc(i,j,k,nfcdiff,nfc2)
                     enddo
                  enddo
               enddo
            else
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmd(i,j,k,nfc1,nfc2)=0.d0
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_c

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma(rhol,ns,fmd)

      implicit none
      real(8),intent(in):: rhol
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(3,3,nfcmax2,nfcmax):: fms1,fms2
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      real(8),dimension(nthmax2,nphmax2):: gjl
      real(8),dimension(:,:,:,:),allocatable :: gma,mma,dmma 
     &                                                !(3,3,nthmax,nphmax)
      real(8),dimension(:,:),allocatable:: gj         !(nthmax,nphmax)

      real(8),dimension(4):: rhoa
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

      do nf2=1,nfcmax
         do nf1=1,nfcmax
            do k=1,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,nf1,nf2)=0.d0
                  enddo
               enddo
            enddo
         enddo
      enddo

      call wmfem_metrics(rhol,nthmax2,nphmax2,gma,mma,dmma,gj)

      call wmfem_disp(rhol,nthmax2,nphmax2,ns,fms)

      do j=1,3
         do i=1,3
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nphnfc2(nfc1)
                  fv1(nth,nph)=cfactor*fmc(i,j,1,nfc1,nfc2)
     &                                *gjl(nth,nph)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nphnfc2(nfc1)
                  fmc(i,j,1,nfc1,nfc2)=fv1f(nth,nph)
               enddo
            enddo    
         enddo
      enddo    

      do nfc2=1,nfcmax
         do nfc1=1,nfcmax
            do k=2,4
               do j=1,3
                  do i=1,3
                     fmc(i,j,k,nfc1,nfc2)=0.d0
                  end do
               end do
            end do
         end do
      end do

!     ----- Fit to fmd and adjust m and n -----

      do nfc2=1,nfcmax
         nth2=nthfnfc(nfc2)
         nph2=nphfnfc(nfc2)
         do nfc1=1,nfcmax
            nth1=nthfnfc(nfc1)
            nph1=nphfnfc(nfc1)

            nphdiff=nph1-nph2
            nthdiff=nth1-nth2
            nfcdiff=nthmax2*nphdiff+nthdiff+nfcmax

            if(nfcdiff.ge.1.and.nfcdiff.le.nfcmax) then
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmd(i,j,k,nfc1,nfc2)
     &                       =fmc(i,j,k,nfcdiff,nfc2)
                     enddo
                  enddo
               enddo
            else
               do k=1,4
                  do j=1,3
                     do i=1,3
                        fmd(i,j,k,nfc1,nfc2)=0.d0
                     enddo
                  enddo
               enddo
            endif
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
