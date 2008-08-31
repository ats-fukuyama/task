!     $Id$

!     ***** wmfem main routine *****

      subroutine wmfem(nrmax,nthmax,nphmax,nsmax,rhoa,cef,cbf,cpp,cpa)

      implicit none
      integer,intent(in):: nrmax,nthmax,nphmax,nsmax
      real(8),dimension(nrmax),intent(in):: rhoa
      complex(8),dimension(3,nthmax,nphmax,nrmax),intent(out):: cef,cbf
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
      complex(8),dimension(:,:),pointer,save:: fma,fmax !(mwmax,mlmax)
      complex(8),dimension(:,:,:),pointer,save:: fms !(mwmax,mlmax,0:nsmax)
      complex(8),dimension(:),pointer,save:: fvb,fvx !(mlmax)

      integer,dimension(:),pointer,save :: nthnfc,mmnfc
      integer,dimension(:),pointer,save :: nphnfc,nnnfc
      integer,dimension(:),pointer,save :: nthnfc2,mmnfc2
      integer,dimension(:),pointer,save :: nphnfc2,nnnfc2

!     ***** define array size  *****

      nfcmax=nthmax*nphmax      ! size of block matrix 
                                !    (number of Fourier components)
      mlmax=6*nfcmax*nrmax      ! length of coeffient matrix and source vector
      mwmax=4*6*nfcmax-1        ! width of coefficient matrix
      if(nthmax.eq.1) then
         nthmax2=1
      else
         nthmax2=2*nthmax       ! number of poloidal modes of coefficients
      endif
      if(nphmax.eq.1) then
         nphmax2=1
      else
         nphmax2=2*nphmax       ! number of toroidal modes of coefficients
      endif
      nfcmax2=nthmax2*nphmax2   ! size of block matrix of coefficients

!     ***** get additional parameters *****

      call get_wmparm(crf,nth0,nph0,idbgwm)

!     ***** allocate matrix and vector *****

      call wmfem_allocate

      if(nrmax.eq.0) return   ! matrix and vector deallocated

!     ***** setup Fourier component index *****

      call wmfem_setup_index

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

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(associated(fma)) deallocate(fma)
         if(associated(fmax)) deallocate(fmax)
         allocate(fma(mwmax,mlmax))
         allocate(fmax(mwmax,mlmax))
      endif
      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save).or.
     &   (nsmax.ne.nsmax_save)) then
         if(associated(fms)) deallocate(fms)
         allocate(fms(mwmax,mlmax,0:nsmax))
      endif
      if(mlmax.ne.mlmax_save) then
         if(associated(fvb)) deallocate(fvb)
         if(associated(fvx)) deallocate(fvx)
         allocate(fvb(mlmax))
         allocate(fvx(mlmax))
      endif

      if(nfcmax.ne.nfcmax_save) then
         if(associated(nthnfc)) deallocate(nthnfc)
         if(associated(nphnfc)) deallocate(nphnfc)
         if(associated(mmnfc)) deallocate(mmnfc)
         if(associated(nnnfc)) deallocate(nnnfc)
         if(nfcmax.ne.0) allocate(nthnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nphnfc(nfcmax))
         if(nfcmax.ne.0) allocate(mmnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nnnfc(nfcmax))

         if(allocated(nthnfc2)) deallocate(nthnfc2)
         if(allocated(nphnfc2)) deallocate(nphnfc2)
         if(allocated(mmnfc2)) deallocate(mmnfc2)
         if(allocated(nnnfc2)) deallocate(nnnfc2)
         if(nfcmax2.ne.0) allocate(nthnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nphnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(mmnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nnnfc2(nfcmax2))
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
            mmnfc(nfc)=nth0
         enddo
      else
         do nph=1,nphmax
         do nth=1,nthmax
            nfc=nthmax*(nph-1)+nth
            nthnfc(nfc)=nth
         enddo
         do nth=1,nthmax/2
            nfc=nthmax*(nph-1)+nth
            mmnfc(nfc)=nth0+nth-1
            nfc=nthmax*(nph-1)+nthmax+1-nth
            mmnfc(nfc)=nth0-nth
         enddo
         enddo
      endif
      if(nphmax.eq.1) then
         do nth=1,nthmax
            nfc=nth
            nphnfc(nfc)=1
            nnnfc(nfc)=nph0
         enddo
      else
         do nth=1,nthmax
         do nph=1,nphmax
            nfc=nthmax*(nph-1)+nth
            nphnfc(nfc)=nph
         enddo
         do nph=1,nphmax/2
            nfc=nthmax*(nph-1)+nth
            nnnfc(nfc)=nph0+nph-1
            nfc=nthmax*(nphmax-nph)+nth
            nnnfc(nfc)=nph0-nph
         enddo
         enddo
      endif

      if(nthmax.eq.1) then
         do nph=1,nphmax2
            nfc2=nph
            nthnfc2(nfc2)=1
            mmnfc2(nfc2)=0
         enddo
      else
         do nph=1,nphmax2
         do nth=1,nthmax2
            nfc2=nthmax2*(nph-1)+nth
            nthnfc2(nfc2)=nth
         enddo
         do nth=1,nthmax
            nfc2=nthmax2*(nph-1)+nth
            mmnfc2(nfc2)=nth-1
            nfc2=nthmax2*(nph-1)+nthmax2+1-nth
            mmnfc2(nfc2)=-nth
         enddo
         enddo
      endif
      if(nphmax.eq.1) then
         do nth=1,nthmax2
            nfc2=nth
            nphnfc2(nfc2)=1
            nnnfc2(nfc2)=0
         enddo
      else
         do nth=1,nthmax2
         do nph=1,nphmax2
            nfc2=nthmax2*(nph-1)+nth
            nphnfc2(nfc2)=nph
         enddo
         do nph=1,nphmax
            nfc2=nthmax2*(nph-1)+nth
            nnnfc2(nfc2)=nph-1
            nfc2=nthmax2*(nphmax2-nph)+nth
            nnnfc2(nfc2)=-nph
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
      complex(8),dimension(3,nthmax,nphmax):: cbfl

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

!     ----- calculate B field -----

      do nr=1,nrmax
         call wmfem_cbf(nr,cbfl)
         do nph=1,nphmax
            do nth=1,nthmax
               cbf(1,nth,nph,nr)=cbfl(1,nth,nph)
               cbf(2,nth,nph,nr)=cbfl(2,nth,nph)
               cbf(3,nth,nph,nr)=cbfl(3,nth,nph)
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

!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate

      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: 
     &     fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax):: 
     &     fmd1,fmd2,fmd3,fmd4
      complex(8),dimension(nphmax,nthmax,3):: fvb_nr
      real(8):: drho,rkth,rkph,rkth0,rho0,rho1,rho2,rho3,rho4
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,mm,nn
      integer:: ns,nfc1,nfc2
      complex(8):: csum,f1,f2,f3,f4
      integer:: id_base=1
      real(8):: angl=0.d0

      do ns=0,nsmax             ! loop for vacuum and plasma species

         do ml=1,mlmax             ! clear fma
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         enddo
         enddo

         do nr=1,nrmax-1        ! loop for elements
            rho1=rhoa(nr)
            rho2=(2.d0*rhoa(nr)+rhoa(nr+1))/3.d0
            rho3=(rhoa(nr)+2.d0*rhoa(nr+1))/3.d0
            rho4=rhoa(nr+1)
            drho=rhoa(nr+1)-rhoa(nr)

            if(ns.eq.0) then
               if(idbgwm.eq.0) then
                  call wmfem_calculate_vacuum(rho1,fmd1)
                  call wmfem_calculate_vacuum(rho2,fmd2)
                  call wmfem_calculate_vacuum(rho3,fmd3)
                  call wmfem_calculate_vacuum(rho4,fmd4)
               else if(idbgwm.eq.1) then
                  call wmfem_calculate_vacuum_c(rho1,fmd1)
                  call wmfem_calculate_vacuum_c(rho2,fmd2)
                  call wmfem_calculate_vacuum_c(rho3,fmd3)
                  call wmfem_calculate_vacuum_c(rho4,fmd4)
               endif
            else
               if(idbgwm.eq.0) then
                  call wmfem_calculate_plasma(rho1,ns,fmd1)
                  call wmfem_calculate_plasma(rho2,ns,fmd2)
                  call wmfem_calculate_plasma(rho3,ns,fmd3)
                  call wmfem_calculate_plasma(rho4,ns,fmd4)
               else if(idbgwm.eq.1) then
                  call wmfem_calculate_plasma_c(rho1,ns,fmd1)
                  call wmfem_calculate_plasma_c(rho2,ns,fmd2)
                  call wmfem_calculate_plasma_c(rho3,ns,fmd3)
                  call wmfem_calculate_plasma_c(rho4,ns,fmd4)
               endif
            endif

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
         mm=mmnfc(nfc)
         ml=6*nfcmax*(nr-1)+6*(nfc-1)
         if(mm.eq.0) then
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
            enddo
            fma(mc,ml+3)=1.d0
         elseif(abs(mm).eq.1) then
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
               fma(mw,ml+5) = 0.d0
            enddo
            fma(mc-2,ml+3)=1.d0
            fma(mc  ,ml+3)=ci*mm
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
            nph=nphnfc(nfc)
            ml=6*nfcmax*(nr-1)+6*(nfc-1)
            fvb(ml+1)=fvb_nr(nph,nth,1)
            fvb(ml+3)=fvb_nr(nph,nth,2)
            fvb(ml+5)=fvb_nr(nph,nth,3)
         enddo
      enddo

      return
      end subroutine wmfem_calculate

!----- calculate coefficint matrix fma (E rho,theta,phi )-----

      subroutine wmfem_calculate_vacuum(rho,fmd)

      implicit none
      real(8),intent(in):: rho
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8),dimension(3,3,3,3,nfcmax2):: fmv1
      complex(8),dimension(3,3,3,nfcmax2):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2):: fmv4
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,nth1,nth2,nthdiff,nph,nph1,nph2,nphdiff
      integer:: imn,imn1,imn2,nfcdiff

      call wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

      do i=1,3
         do j=1,3
            do imn1=1,3
               do imn2=1,3
                  do nfc2=1,nfcmax2
                     nth=nthnfc2(nfc2)
                     nph=nphnfc2(nfc2)
                     fv1(nth,nph)=fmv1(i,j,imn1,imn2,nfc2)
                  enddo
                  call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
                  do nfc2=1,nfcmax2
                     nth=nthnfc2(nfc2)
                     nph=nphnfc2(nfc2)
                     fmv1(i,j,imn1,imn2,nfc2)=fv1f(nth,nph)
                  enddo
               enddo
            enddo    
            do imn1=1,3
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  fv1(nth,nph)=fmv2(i,j,imn1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  fmv2(i,j,imn1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do imn1=1,3
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  fv1(nth,nph)=fmv3(i,j,imn1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  fmv3(i,j,imn1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nphnfc2(nfc2)
               fv1(nth,nph)=fmv4(i,j,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nphnfc2(nfc2)
               fmv4(i,j,nfc2)=fv1f(nth,nph)
            enddo
         enddo
      enddo

      do nfc1=1,nfcmax          ! Fit to fmd and adjust m and n
         nph1=nnnfc(nfc1)
         nth1=mmnfc(nfc1)
         do nfc2=1,nfcmax
            nph2=nnnfc(nfc2)
            nth2=mmnfc(nfc2)

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
      end subroutine wmfem_calculate_vacuum

!----- calculate vacuum matrix fms -----

      subroutine wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

      implicit none
      real(8),intent(in):: rho
      complex(8),dimension(3,3,3,3,nfcmax2),intent(out):: fmv1
      complex(8),dimension(3,3,3,nfcmax2),intent(out):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2),intent(out):: fmv4

      real(8),dimension(3,3,nthmax2,nphmax2) :: gma,muma,dmuma 
      real(8),dimension(nthmax2,nphmax2):: gja

      complex(8),dimension(3,3,3):: cq
      complex(8),dimension(3,3):: cp
      integer:: i,j,k,l,nthm,nthp,nphm,nphp
      integer:: imn1,imn2
      integer:: nph,nth
      real(8):: dph,dth
      complex(8):: csum1,csum2,csum3,csum4,cfactor
      integer:: nrl,nfc2
      real(8):: drhob,muma3b,muma2b,drhoa,muma3a,muma2a,muma30,muma20
      real(8):: muma3d,muma2d,gj
      
      cfactor=(2*pi*crf*1.d6)**2/vc**2

      call wmfem_tensors(rho,gma,muma,dmuma,gja)

      do nfc2=1,nfcmax2
         nth=nthnfc2(nfc2)
         nph=nphnfc2(nfc2)

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
         gj=gja(nth,nph)

         do j=1,3
            cq(1,j,1)=((muma(3,j,nthp,nph)
     &                 -muma(3,j,nthm,nph))/dth
     &                -(muma(2,j,nth,nphp)
     &                 -muma(2,j,nth,nphm))/dph )/gj
            cq(1,j,2)=+ci*muma(3,j,nth,nph)/gj
            cq(1,j,3)=-ci*muma(2,j,nth,nph)/gj

            cq(2,j,1)=((muma(1,j,nth,nphp)
     &                 -muma(1,j,nth,nphm))/dph
     &                 -dmuma(3,j,nth,nph))/gj
            cq(2,j,2)=0.d0
            cq(2,j,3)=+ci*muma(1,j,nth,nph)/gj

            cq(3,j,1)=(dmuma(2,j,nth,nph)
     &               -(muma(1,j,nthp,nph)
     &                -muma(1,j,nthm,nph))/dth )/gj
            cq(3,j,2)=-ci*muma(1,j,nth,nph)/gj
            cq(3,j,3)=0.d0

            cp(1,j)=0.d0
            cp(2,j)=-muma(3,j,nth,nph)/gj
            cp(3,j)= muma(2,j,nth,nph)/gj
         enddo

         do imn2=1,3
            do imn1=1,3
               do j=1,3
                  do i=1,3
                     csum1=0.d0
                     do k=1,3
                        do l=1,3
                           csum1=csum1+conjg(cq(k,i,imn1))
     &                                *gma(k,l,nth,nph)
     &                                *cq(l,j,imn2) *gj
                        enddo
                     enddo
                     if(i.eq.j.and.imn1.eq.1.and.imn2.eq.1) then
                        fmv1(i,j,imn1,imn2,nfc2)
     &                       =csum1-cfactor*gj
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
     &                             *gma(k,l,nth,nph)
     &                             *cq(l,j,imn1)*gj
                        csum3=csum3+conjg(cq(k,i,imn1))
     &                             *gma(k,l,nth,nph)
     &                             *cp(l,j)*gj
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
     &                    *gma(k,l,nth,nph)
     &                    *cp(l,j) *gj
                  enddo
               enddo
               fmv4(i,j,nfc2)=csum4
            enddo
         enddo

      enddo
      return
      end subroutine wmfem_calculate_vacuum_sub

!----- calculate coefficint matrix fmd (E cylindrical +,-,para)-----

      subroutine wmfem_calculate_vacuum_c(rho,fmd)

      implicit none
      real(8),intent(in):: rho
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

         if(rho.ne.0.d0) then
            rkth=nth/(ra*rho)
            rkth0=1.d0/(ra*rho)
         else
            rkth=0.d0
            rkth0=0.d0
         endif
         rkph=nph/rr

         do nfc1=1,nfcmax2
            fmc(1,1,1,nfc1,nfc2)= rho*(cfactor-rkph**2-rkth**2)
            fmc(1,2,1,nfc1,nfc2)= rho*(-ci*rkth*rkth0)
            fmc(2,1,1,nfc1,nfc2)= rho*(+ci*rkth*rkth0)
            fmc(2,2,1,nfc1,nfc2)= rho*(cfactor-rkph**2-rkth0**2)
            fmc(2,3,1,nfc1,nfc2)= rho*(rkth*rkph)
            fmc(3,2,1,nfc1,nfc2)= rho*(rkth*rkph)
            fmc(3,3,1,nfc1,nfc2)= rho*(cfactor-rkth**2)
                  
            fmc(2,1,2,nfc1,nfc2)= rho*( ci*rkth)/ra
            fmc(2,2,2,nfc1,nfc2)= rho*(  -rkth0)/ra
            fmc(3,1,2,nfc1,nfc2)= rho*( ci*rkph)/ra
                  
            fmc(1,2,3,nfc1,nfc2)= rho*(-ci*rkth)/ra
            fmc(1,3,3,nfc1,nfc2)= rho*(-ci*rkph)/ra
            fmc(2,2,3,nfc1,nfc2)= rho*(  -rkth0)/ra
                  
            fmc(1,1,4,nfc1,nfc2)= 0.d0
            fmc(2,2,4,nfc1,nfc2)= rho*(-1.d0)/ra**2
            fmc(3,3,4,nfc1,nfc2)= rho*(-1.d0)/ra**2
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
                     nth=nthnfc2(nfc1)
                     nph=nphnfc2(nfc1)
                     fmc(i,j,k,nfc1,nfc2)=fv1f(nth,nph)
                  enddo
               enddo
            enddo    
         enddo
      enddo    

!     ----- Fit to fmd and adjust m and n -----

      do nfc2=1,nfcmax
         nth2=mmnfc(nfc2)
         nph2=nnnfc(nfc2)
         do nfc1=1,nfcmax
            nth1=mmnfc(nfc1)
            nph1=nnnfc(nfc1)

            nphdiff=nph1-nph2
            nthdiff=nth1-nth2
            nfcdiff=nthmax2*nphdiff+nthdiff+nfcmax

            do k=1,4
               do j=1,3
                  do i=1,3
                     fmd(i,j,k,nfc1,nfc2)=fmc(i,j,k,nfcdiff,nfc2)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_c

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma(rho,ns,fmd)

      implicit none
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      real(8),dimension(3,3,nthmax2,nphmax2) :: gma,muma,dmuma 
      real(8),dimension(nthmax2,nphmax2):: gja

      complex(8):: cfactor
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nph1,nph2,nphdiff,nth1,nth2,nthdiff,nfcdiff

      cfactor=(2*pi*crf*1.d6)**2/vc**2

      call wmfem_tensors(rho,gma,muma,dmuma,gja)

      call wmfem_disp_tensor(rho,ns,fmc)

      do j=1,3
         do i=1,3
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nphnfc2(nfc1)
                  fv1(nth,nph)=cfactor*fmc(i,j,1,nfc1,nfc2)
     &                                *gja(nth,nph)
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
         nth2=mmnfc(nfc2)
         nph2=nnnfc(nfc2)
         do nfc1=1,nfcmax
            nth1=mmnfc(nfc1)
            nph1=nnnfc(nfc1)

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

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma_c(rho,ns,fmd)

      implicit none
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f

      complex(8):: cfactor
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nph1,nph2,nphdiff,nth1,nth2,nthdiff,nfcdiff

      cfactor=(2*pi*crf*1.d6)**2/vc**2

      call wmfem_disp_tensor(rho,ns,fmc)

      do j=1,3
         do i=1,3
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nphnfc2(nfc1)
                  fv1(nth,nph)=cfactor*fmc(i,j,1,nfc1,nfc2)*rho
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
         nth2=mmnfc(nfc2)
         nph2=nnnfc(nfc2)
         do nfc1=1,nfcmax
            nth1=mmnfc(nfc1)
            nph1=nnnfc(nfc1)

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

      end subroutine wmfem_calculate_plasma_c

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

!     ****** CALCULATE METRIC AND CONVERSION TENSOR ******

      SUBROUTINE wmfem_tensors(rho,gma,muma,dmuma,gja)

      IMPLICIT NONE
      real(8),intent(in):: rho
      real(8),intent(out),dimension(3,3,nthmax2,nphmax2):: 
     &     gma,muma,dmuma
      real(8),intent(out),dimension(nthmax2,nphmax2):: gja
      real(8),dimension(3,3):: gm,gmp,gmm,mum,mump,mumm
      real(8):: th,ph,dth,dph,drhom,drhop,gj,gjp,gjm
      real(8):: babs,bsupth,bsupph,rhol
      integer:: nth,nph,i,j

      dth=2.d0*pi/nthmax2
      dph=2.d0*pi/nphmax2
      IF(rho.EQ.0.d0) THEN
         rhol=1.d-6
         drhom=0.d0
      else
         rhol=rho
         drhom=1.d-6
      endif
      drhop=1.d-6
      DO nph=1,nphmax2
         ph=dph*(nph-1)
         DO nth=1,nthmax2
            th=dth*(nth-1)

            CALL wmfem_metrics(rhol,th,ph,gm,gj)
            CALL wmfem_magnetic(rhol,th,ph,babs,bsupth,bsupph)
            CALL wmfem_rotation_tensor(gm,gj,babs,bsupth,bsupph,mum)

            CALL wmfem_metrics(rhol-drhom,th,ph,gmm,gjm)
            CALL wmfem_magnetic(rhol-drhom,th,ph,babs,bsupth,bsupph)
            CALL wmfem_rotation_tensor(gmm,gjm,babs,bsupth,bsupph,mumm)

            CALL wmfem_metrics(rhol+drhop,th,ph,gmp,gjp)
            CALL wmfem_magnetic(rhol+drhop,th,ph,babs,bsupth,bsupph)
            CALL wmfem_rotation_tensor(gmp,gjp,babs,bsupth,bsupph,mump)

            DO j=1,3
               DO i=1,3
                  gma(i,j,nth,nph)=gm(i,j)
                  muma(i,j,nth,nph)=mum(i,j)
                  dmuma(i,j,nth,nph)
     &                 =(mump(i,j)-mumm(i,j))/(drhom+drhop)
               END DO
            END DO
            gja(nth,nph)=gj

         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE wmfem_tensors

!        ****** CALCULATE CONVERSION TENSOR ******

         SUBROUTINE wmfem_rotation_tensor(gm,gj,babs,bsupth,bsupph,mum)

         IMPLICIT NONE
         real(8),dimension(3,3),intent(in):: gm
         real(8),intent(in):: gj,babs,bsupth,bsupph
         real(8),dimension(3,3),intent(out):: mum
         real(8):: tc2,tc3,gf11,gf12

!        ----- Calculate rotation matrix mum -----

         tc2=bsupth/babs
         tc3=bsupph/babs

!        ----- gf11=gj*SQRT(gm11) -----

         gf11 = SQRT(gm(2,2)*gm(3,3)-gm(2,3)*gm(3,2))
         gf12 =(tc2*(gm(2,3)*gm(1,2)-gm(2,2)*gm(1,3))
     &         +tc3*(gm(3,3)*gm(1,2)-gm(2,3)*gm(1,3)))
            
         mum(1,1)= gj/gf11
         mum(2,1)= 0.D0
         mum(3,1)= 0.D0
         mum(1,2)= gf12/gf11
         mum(2,2)= tc3*gf11
         mum(3,2)=-tc2*gf11
         mum(1,3)= tc2*gm(1,2)+tc3*gm(1,3)
         mum(2,3)= tc2*gm(2,2)+tc3*gm(2,3)
         mum(3,3)= tc2*gm(3,2)+tc3*gm(3,3)

         end subroutine wmfem_rotation_tensor

!     ****** CALCULATE METRIC AND CONVERSION TENSOR ******

         SUBROUTINE wmfem_disp_tensor(rho,ns,fmc)

         IMPLICIT NONE
         real(8),intent(in):: rho
         integer,intent(in):: ns
         complex(8),dimension(3,3,4,nfcmax2,nfcmax),intent(out):: fmc
         complex(8),dimension(3,3):: fms

         real(8):: th,ph,dth,dph
         integer:: nth,nph,mm,nn,i,j,nfc2,nfc

         dth=2.d0*pi/nthmax2
         dph=2.d0*pi/nphmax2
         DO nfc2=1,nfcmax2
            nth=nthnfc2(nfc2)
            nph=nphnfc2(nfc2)
            th=dth*(nth-1)
            ph=dph*(nph-1)
            DO nfc=1,nfcmax
               mm=mmnfc(nfc)
               nn=nnnfc(nfc)
               CALL wmfem_dielectric(rho,th,ph,mm,nn,ns,fms)
               DO j=1,3
                  DO i=1,3
                     fmc(i,j,1,nfc2,nfc)=fms(i,j)
                  END DO
               END DO
            END DO
         END DO
         RETURN
         END SUBROUTINE wmfem_disp_tensor

!     ****** CALCULATE METRIC AND CONVERSION TENSOR ******

         SUBROUTINE wmfem_dielectric(rho,th,ph,mm,nn,ns,fms)

         IMPLICIT NONE
         real(8),intent(in):: rho,th,ph
         integer,intent(in):: mm,nn,ns
         complex(8),dimension(3,3),intent(out):: fms
         complex(8):: cw,ckpara,ckperp
         real(8):: babs,bsupth,bsupph

         cw=2*pi*crf*1.d6
         CALL wmfem_magnetic(rho,th,ph,babs,bsupth,bsupph)
         ckpara=mm*bsupth/babs+nn*bsupph/babs
         ckperp=(0.d0,0.d0)

         CALL dpcalc(cw,ckpara,ckperp,rho,babs,ns,fms)

         END SUBROUTINE wmfem_dielectric

!----- calculate local wave magnetic field -----

      subroutine wmfem_cbf(nr,cbf)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,nthmax,nphmax):: cbf

      real(8),dimension(3,3,nthmax2,nphmax2) :: gma,muma,dmuma 
      real(8),dimension(nthmax2,nphmax2):: gja

      complex(8),dimension(3,3,3):: cq
      complex(8),dimension(3,3):: cp
      complex(8),dimension(3,3,3,nfcmax2):: cqa
      complex(8),dimension(3,3,nfcmax2):: cpa
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      integer:: i,j,k,l,nthm,nthp,nphm,nphp
      integer:: nfc2,nph,nth,imn,ml
      integer:: nfc1,nph1,nth1,nth2,nph2
      integer:: nn1,mm1,nn2,mm2,mmdiff,nndiff,nfcfdiff
      real(8):: rho,dph,dth,gj
      complex(8):: cfactor
      
      cfactor=vc/(2*pi*crf*1.d6)

      rho=rhoa(nr)

      call wmfem_tensors(rho,gma,muma,dmuma,gja)

!     ----- calculation rot coefficients cpa and cqa -----

      do nfc2=1,nfcmax2
         nth=nthnfc2(nfc2)
         nph=nphnfc2(nfc2)
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
         gj=gja(nth,nph)

         do j=1,3
            cq(1,j,1)=((muma(3,j,nthp,nph)
     &                 -muma(3,j,nthm,nph))/dth
     &                -(muma(2,j,nth,nphp)
     &                 -muma(2,j,nth,nphm))/dph )/gj
            cq(1,j,2)=+ci*muma(3,j,nth,nph)/gj
            cq(1,j,3)=-ci*muma(2,j,nth,nph)/gj

            cq(2,j,1)=((muma(1,j,nth,nphp)
     &                 -muma(1,j,nth,nphm))/dph
     &                 -dmuma(3,j,nth,nph))/gj
            cq(2,j,2)=0.d0
            cq(2,j,3)=+ci*muma(1,j,nth,nph)/gj

            cq(3,j,1)=(dmuma(2,j,nth,nph)
     &               -(muma(1,j,nthp,nph)
     &                -muma(1,j,nthm,nph))/dth )/gj
            cq(3,j,2)=-ci*muma(1,j,nth,nph)/gj
            cq(3,j,3)=0.d0

            cp(1,j)=0.d0
            cp(2,j)=-muma(3,j,nth,nph)/gj
            cp(3,j)= muma(2,j,nth,nph)/gj
         enddo

         do imn=1,3
            do i=1,3
               do j=1,3
                  cqa(i,j,imn,nfc2)=0.d0
                  do l=1,3
                     cqa(i,j,imn,nfc2)=cqa(i,j,imn,nfc2)
     &                 +gma(i,l,nth,nph)*cq(l,j,imn)*gj
                  enddo
               enddo
            enddo
         enddo
         do i=1,3
            do j=1,3
               cpa(i,j,nfc2)=0.d0
               do l=1,3
                  cpa(i,j,nfc2)=cpa(i,j,nfc2)
     &                 +gma(i,l,nth,nph)*cp(l,j)*gj
               enddo
            enddo
         enddo
      enddo

!     ----- FFT of cqa and cpa -----
            
      do i=1,3
         do j=1,3
            do imn=1,3
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  fv1(nth,nph)=cqa(i,j,imn,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  cqa(i,j,imn,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nphnfc2(nfc2)
               fv1(nth,nph)=cpa(i,j,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nphnfc2(nfc2)
               cpa(i,j,nfc2)=fv1f(nth,nph)
            enddo
         enddo
      enddo
      
!     ----- calculated cbf -----
      
      do nfc1=1,nfcmax
         nph1=nphnfc(nfc1)
         nth1=nthnfc(nfc1)
         nn1=nnnfc(nfc1)
         mm1=mmnfc(nfc1)
         do nfc2=1,nfcmax
            nph2=nphnfc(nfc2)
            nth2=nthnfc(nfc2)
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)

            nndiff=nn1-nn2
            mmdiff=mm1-mm2
            nfcfdiff=nthmax*nndiff+mmdiff+nfcmax
            ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nph2-1)+6*(nth2-1)

            do i=1,3
               cbf(i,nth1,nph1)=0.d0
               do j=1,3
                  cbf(i,nth1,nph1)=cbf(i,nth1,nph1)
     &                 +cfactor*(
     &                  (cqa(i,j,1,nfcfdiff)
     &                  +cqa(i,j,2,nfcfdiff)*mm2
     &                  +cqa(i,j,3,nfcfdiff)*nn2)*fvx(ml+2*j-1)
     &                 + cpa(i,j,  nfcfdiff)       *fvx(ml+2*j  ))
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_cbf

      end subroutine wmfem
