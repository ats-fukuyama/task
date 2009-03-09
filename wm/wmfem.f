!     $Id$

!     ***** wmfem main routine *****

      subroutine wmfem_main

      use wmfem_com
      implicit none
      integer:: ierr

!     ***** metric setup  *****

         CALL wmfem_setg(ierr)
         IF(IERR.NE.0) RETURN
         CALL wmfem_setj(ierr)
         IF(IERR.NE.0) RETURN
         CALL get_wmfem_size(nrmax,nthmax,nphmax,nsmax,ierr)
         IF(IERR.NE.0) RETURN

!     ***** define array size  *****

      nfcmax=nthmax*nphmax      ! size of block matrix 
                                !    (number of Fourier components)
      mlmax=6*nfcmax*nrmax      ! length of coeffient matrix and source vector
      mwmax=4*6*nfcmax-1        ! width of coefficient matrix

      if(nthmax.eq.1) then
         nthmax2=1
      else
         nthmax2=nthmax*2
      endif
      if(nphmax.eq.1) then
         nphmax2=1
      else
         nphmax2=nphmax*2
      endif
      nfcmax2=nthmax2*nphmax2

!     ***** get additional parameters *****

      call get_wmfem_parm(crf,nth0,nph0,mdlwmf,mdlwmd)

!     ***** allocate matrix and vector *****

      call wmfem_allocate

      if(nrmax.eq.0) return   ! matrix and vector deallocated

!     ***** setup Fourier component index *****

      call wmfem_setup_index

!     ***** calculate coefficient matrix *****

      call wmfem_calculate_matrix

!     ***** solve matrix equation *****

      call wmfem_solve

!     ***** calculate efield *****

      call wmfem_calculate_efield
      call wmfem_efld(cef)

      return

      contains

!     ----- allocate arrays -----

      subroutine wmfem_allocate

      implicit none
      integer,save:: nrmax_save=0,nthmax_save=0,nphmax_save=0
      integer,save:: mwmax_save=0,mlmax_save=0
      integer,save:: nsmax_save=0,nfcmax_save=0
      integer,save:: mdlwmd_save=0 

      if((nrmax.ne.nrmax_save).or.(nthmax.ne.nthmax_save).or. 
     &   (nphmax.ne.nphmax_save)) then
         if(associated(cef)) deallocate(cef)
         allocate(cef(3,nthmax,nphmax,nrmax))
         if(associated(cbf)) deallocate(cbf)
         allocate(cbf(3,nthmax,nphmax,nrmax))
         if(associated(cpp)) deallocate(cpp)
         allocate(cpp(nthmax,nphmax,nthmax2,nphmax2,nrmax,0:nsmax))
         if(associated(cpa)) deallocate(cpa)
         allocate(cpa(nthmax,nphmax))
      endif

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(associated(fma)) deallocate(fma)
         allocate(fma(mwmax,mlmax))
      endif

      if(mdlwmd.ge.1) then
      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save).or.
     &   (nsmax.ne.nsmax_save).or.mdlwmd.ne.mdlwmd_save) then
         if(associated(fms)) deallocate(fms)
         if(mdlwmd.ge.1) then
            allocate(fms(mwmax,12*nfcmax,nrmax,0:nsmax))
         endif 
      endif
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

         if(associated(nthnfc2)) deallocate(nthnfc2)
         if(associated(nphnfc2)) deallocate(nphnfc2)
         if(associated(mmnfc2)) deallocate(mmnfc2)
         if(associated(nnnfc2)) deallocate(nnnfc2)
         if(nfcmax2.ne.0) allocate(nthnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nphnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(mmnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nnnfc2(nfcmax2))
      endif
      mwmax_save=mwmax
      mlmax_save=mlmax
      nsmax_save=nsmax
      nfcmax_save=nfcmax
      mdlwmd_save=mdlwmd
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

      do nph=1,nphmax2
         do nth=1,nthmax2
            nfc2=nthmax2*(nph-1)+nth
            nthnfc2(nfc2)=nth
         enddo
         do nth=1,nthmax
            nfc2=nthmax2*(nph-1)+nth
            mmnfc2(nfc2)=nth0+nth-1
            nfc=nthmax2*(nph-1)+nthmax2+1-nth
            mmnfc2(nfc2)=nth0-nth
         enddo
      enddo
      
      do nth=1,nthmax2
         do nph=1,nphmax2
            nfc2=nthmax2*(nph-1)+nth
            nphnfc2(nfc2)=nph
         enddo
         do nph=1,nphmax
            nfc2=nthmax2*(nph-1)+nth
            nnnfc2(nfc2)=nph0+nph-1
            nfc2=nthmax2*(nphmax2-nph)+nth
            nnnfc2(nfc2)=nph0-nph
         enddo
      enddo

      return
      end subroutine wmfem_setup_index


!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_matrix

      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmd
      complex(8),dimension(3,3,4,nfcmax,nfcmax)::  fmd1,fmd2,fmd3,fmd4
      complex(8),dimension(nphmax,nthmax,3):: fvb_nr
      complex(8),dimension(mwmax,12*nfcmax):: fml
      real(8):: drho,rkth,rkph,rkth0,rho0,rho1,rho2,rho3,rho4
      integer:: nr,ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,mm,nn,mll
      integer:: ns,nfc1,nfc2,ml1,mw1,mr
      complex(8):: csum,f1,f2,f3,f4,cx,cy,cxd,cyd
      integer:: id_base=1
      real(8):: angl=0.d0

      do ml=1,mlmax             ! clear fma
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         end do
      end do

      do ns=0,nsmax             ! loop for vacuum and plasma species

         do nr=1,nrmax-1        ! loop for elements

            call wmfem_calculate_local(nr,ns,fml)

            if(mdlwmd.ge.1) then
               do ml=1,12*nfcmax
                  do mw=1,mwmax
                     fms(mw,ml,nr,ns)=fml(mw,ml)
                  enddo
               enddo
            endif

            do ml=1,12*nfcmax
               do mw=1,mwmax
                  fma(mw,6*nfcmax*(nr-1)+ml)
     &           =fma(mw,6*nfcmax*(nr-1)+ml)+fml(mw,ml)
               enddo
            enddo
         enddo
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

!----- boundary conditions -----

      mc=(mwmax+1)/2
      mr=6*nfcmax

      nr=1
      do nfc=1,nfcmax
         mm=mmnfc(nfc)
         mll=6*(nfc-1)
         ml=6*nfcmax*(nr-1)+mll
         if(mm.eq.0) then
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fms(mw,mll+3,nr,ns)=0.d0
                  enddo
               end if
            enddo
            fma(mc,ml+3)=1.d0
            if(mdlwmd.ge.1) fms(mc,mll+3,nr,0)=1.d0
            fvb(ml+3)=0.d0
         elseif(abs(mm).eq.1) then
            do mw=4,mwmax
               cx= fma(mw  ,ml+1)
               cy= fma(mw-2,ml+3)
               fma(mw  ,ml+1)=cx +ci*mm*cy
               fma(mw-2,ml+3)=cx -ci*mm*cy
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     cx= fms(mw  ,mll+1,nr,ns)
                     cy= fms(mw-2,mll+3,nr,ns)
                     fms(mw  ,mll+1,nr,ns)=cx +ci*mm*cy
                     fms(mw-2,mll+3,nr,ns)=cx -ci*mm*cy
                  enddo
               endif
            enddo
            do mw=-mc+5,mc
               if(ml+mw.ge.1.AND.ml+mw.LE.12*nfcmax) then
                  cx= fma(mc-mw+1,ml+mw)
                  cy= fma(mc-mw+3,ml+mw)
                  fma(mc-mw+1,ml+mw)=cx -ci*mm*cy
                  fma(mc-mw+3,ml+mw)=cx +ci*mm*cy
                  if(mdlwmd.ge.1) then
                     do ns=0,nsmax
                        cx= fms(mc-mw+1,mll+mw,nr,ns)
                        cy= fms(mc-mw+3,mll+mw,nr,ns)
                        fms(mc-mw+1,mll+mw,nr,ns)=cx -ci*mm*cy
                        fms(mc-mw+3,mll+mw,nr,ns)=cx +ci*mm*cy
                        if(ml+mw.GT.mr.and.mc-mw+1-mr.GT.0) then
                           cx= fms(mc-mw+1-mr,mll+mw-mr,nr,ns)
                           cy= fms(mc-mw+3-mr,mll+mw-mr,nr,ns)
                           fms(mc-mw+1-mr,mll+mw-mr,nr+1,ns)=cx-ci*mm*cy
                           fms(mc-mw+3-mr,mll+mw-mr,nr+1,ns)=cx+ci*mm*cy
                        endif
                     end do
                  endif
               endif
            enddo

            do mw=1,mwmax
               fma(mw,ml+1) = 0.d0
               fma(mw,ml+5) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fms(mw,mll+1,nr,ns)=0.d0
                     fms(mw,mll+5,nr,ns)=0.d0
                  end do
               end if
            end do
            fma(mc  ,ml+1)=1.d0
            fma(mc  ,ml+5)=1.d0
            fvb(ml+1)=0.d0
            fvb(ml+5)=0.d0
            if(mdlwmd.ge.1) then
               fms(mc  ,mll+1,nr,0)=1.d0
               fms(mc  ,mll+5,nr,0)=1.d0
            endif
         else
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
               fma(mw,ml+5) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fms(mw,mll+3,nr,ns)=0.d0
                     fms(mw,mll+5,nr,ns)=0.d0
                  enddo
               endif
            enddo
            fma(mc,ml+3)=1.d0
            fma(mc,ml+5)=1.d0
            fvb(ml+3)=0.d0
            fvb(ml+5)=0.d0
            if(mdlwmd.ge.1) then
               fms(mc,mll+3,nr,0)=1.d0
               fms(mc,mll+5,nr,0)=1.d0
            endif
         endif
      enddo

      nr=nrmax
      do nfc=1,nfcmax
         mll=6*(nfc-1)
         ml=6*nfcmax*(nr-1)+mll
         if(id_base.eq.0) then
            do mw=1,mwmax
               fma(mw,ml+3) = 0.d0
               fma(mw,ml+5) = 0.d0
            enddo
            fma(mc,ml+3) = 1.d0
            fma(mc,ml+5) = 1.d0
            fvb(ml+3)=0.d0
            fvb(ml+5)=0.d0
         else
            do mw=1,mwmax
               fma(mw,ml+2) = 0.d0
               fma(mw,ml+3) = 0.d0
               fma(mw,ml+5) = 0.d0
            enddo
            fma(mc,ml+2) = 1.d0
            fma(mc,ml+3) = 1.d0
            fma(mc,ml+5) = 1.d0
            fvb(ml+2)=0.d0
            fvb(ml+3)=0.d0
            fvb(ml+5)=0.d0
         endif
      enddo

      do nr=1,nrmax
         if(nphmax.gt.1) then
            nn=nphmax/2+1
            do mm=1,nthmax
               mll=6*nthmax*(nn-1)+6*(mm-1)
               ml=6*nthmax*nphmax*(nr-1)+mll
               do mw=1,mwmax
                  fma(mw,ml+1)=0.d0
                  fma(mw,ml+3)=0.d0
                  fma(mw,ml+5)=0.d0
                  if(mdlwmd.ge.1) then
                     do ns=0,nsmax
                        fms(mw,mll+1,nr,ns)=0.d0
                        fms(mw,mll+3,nr,ns)=0.d0
                        fms(mw,mll+5,nr,ns)=0.d0
                        if(nr.ne.1) then
                           fms(mw,mll+6*nfcmax+1,nr-1,ns)=0.d0
                           fms(mw,mll+6*nfcmax+3,nr-1,ns)=0.d0
                           fms(mw,mll+6*nfcmax+5,nr-1,ns)=0.d0
                        endif
                     end do
                  endif
               end do
               fma(mc,ml+1)=1.d0
               fma(mc,ml+3)=1.d0
               fma(mc,ml+5)=1.d0
               fvb(ml+1)=0.d0
               fvb(ml+3)=0.d0
               fvb(ml+5)=0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fms(mc,mll+1,nr,ns)=1.d0
                     fms(mc,mll+3,nr,ns)=1.d0
                     fms(mc,mll+5,nr,ns)=1.d0
                     if(nr.ne.1) then
                        fms(mc,mll+6*nfcmax+1,nr-1,ns)=0.d0
                        fms(mc,mll+6*nfcmax+3,nr-1,ns)=0.d0
                        fms(mc,mll+6*nfcmax+5,nr-1,ns)=0.d0
                     endif
                  enddo
               end if
            end do
         endif
         if(nthmax.gt.1) then
            mm=nthmax/2+1
            do nn=1,nphmax
               mll=6*nthmax*(nn-1)+6*(mm-1)
               ml=6*nthmax*nphmax*(nr-1)+mll
               do mw=1,mwmax
                  fma(mw,ml+1)=0.d0
                  fma(mw,ml+3)=0.d0
                  fma(mw,ml+5)=0.d0
                  if(mdlwmd.ge.1) then
                     do ns=0,nsmax
                        fms(mw,mll+1,nr,ns)=0.d0
                        fms(mw,mll+3,nr,ns)=0.d0
                        fms(mw,mll+5,nr,ns)=0.d0
                        if(nr.ne.1) then
                           fms(mw,mll+6*nfcmax+1,nr-1,ns)=0.d0
                           fms(mw,mll+6*nfcmax+3,nr-1,ns)=0.d0
                           fms(mw,mll+6*nfcmax+5,nr-1,ns)=0.d0
                        endif
                     end do
                  endif
               enddo
               fma(mc,ml+1)=1.d0
               fma(mc,ml+3)=1.d0
               fma(mc,ml+5)=1.d0
               fvb(ml+1)=0.d0
               fvb(ml+3)=0.d0
               fvb(ml+5)=0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fms(mc,mll+1,nr,ns)=1.d0
                     fms(mc,mll+3,nr,ns)=1.d0
                     fms(mc,mll+5,nr,ns)=1.d0
                     if(nr.ne.1) then
                        fms(mc,mll+6*nfcmax+1,nr-1,ns)=0.d0
                        fms(mc,mll+6*nfcmax+3,nr-1,ns)=0.d0
                        fms(mc,mll+6*nfcmax+5,nr-1,ns)=0.d0
                     endif
                  end do
               endif
            end do
         endif
      enddo

      return
      end subroutine wmfem_calculate_matrix

!     ***** Solve matrix equation *****

      subroutine wmfem_solve

      integer:: ml,mw,ierr

!     ----- solve matrix -----

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
      end do

!      write(6,*) 'mlmax,mwmax=',mlmax,mwmax
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      if(ierr.ne.0) write(6,*) '# ierr= ',ierr

      return
      end subroutine wmfem_solve

!     ***** Calculate electric field *****

      subroutine wmfem_calculate_efield

      integer:: nr,nn,mm,ml

      nr=1
      do nn=1,nphmax
         do mm=1,nthmax
            ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
            if(abs(mm).eq.1) then
               cef(1,mm,nn,nr)=(fvx(ml+1)+fvx(ml+3))
               cef(2,mm,nn,nr)=(fvx(ml+1)-fvx(ml+3))/(ci*mm)
               cef(3,mm,nn,nr)=fvx(ml+5)
            else
               cef(1,mm,nn,nr)=fvx(ml+1)
               cef(2,mm,nn,nr)=fvx(ml+3)
               cef(3,mm,nn,nr)=fvx(ml+5)
            endif
         enddo
      enddo
      do nr=2,nrmax
         do nn=1,nphmax
            do mm=1,nthmax
               ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
               cef(1,mm,nn,nr)=fvx(ml+1)
               cef(2,mm,nn,nr)=fvx(ml+3)
               cef(3,mm,nn,nr)=fvx(ml+5)
!               write(6,'(3I4,1P6E11.3)') nr,nn,mm,cef(1,mm,nn,nr),
!     &              cef(2,mm,nn,nr),cef(3,mm,nn,nr)
            enddo
         enddo
      enddo

      return
      end subroutine wmfem_calculate_efield
      end subroutine wmfem_main

!     ***** wmfem post routine *****

      subroutine wmfem_post

      use wmfem_com
      implicit none
      integer:: ierr

!     ***** calculate bfield *****

      call wmfem_calculate_bfield
      call wmfem_bfld(cbf)

!     ***** calculate power *****

      call wmfem_calculate_power
      call wmfem_pabs(cpp,cpa)

      return

      contains

!     ***** Calculate magnetic field *****

      subroutine wmfem_calculate_bfield

      use wmfem_com
      implicit none
      integer:: nr,nth,nph
      complex(8),dimension(3,nthmax,nphmax):: cbfl

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

      return
      end subroutine wmfem_calculate_bfield

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
      integer:: nfc1,nph1,nth1
      integer:: nn1,mm1,nn2,mm2,mmdiff,nndiff,nfcdiff
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
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)

            nndiff=nn1-nn2
            if(nndiff.lt.0) nndiff=nndiff+nphmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1
            ml=6*nthmax*nphmax*(nr-1)+6*(nfc2-1)

!          !!!!! on axis should be considered !!!!!

            do i=1,3
               cbf(i,nth1,nph1)=0.d0
               do j=1,3
                  cbf(i,nth1,nph1)=cbf(i,nth1,nph1)
     &                 +cfactor*(
     &                  (cqa(i,j,1,nfcdiff)
     &                  +cqa(i,j,2,nfcdiff)*mm2
     &                  +cqa(i,j,3,nfcdiff)*nn2)*fvx(ml+2*j-1)
     &                 + cpa(i,j,  nfcdiff)     *fvx(ml+2*j  ))
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_cbf

!     ***** Calculate absorbed power *****

      subroutine wmfem_calculate_power

      implicit none
      integer:: mc,ml,mw,nr,ns,mm,nn,i,j
      integer:: ir,nr1,mm1,nn1,mm2,nn2,nndiff,mmdiff
      integer:: i1,ml1,ml2,nfc,nfc1,nfc2,mll,nr0
      real(8):: factor
      complex(8):: csum,csums,cfactor
      complex(8),dimension(mwmax,12*nfcmax):: fml

      mc=(mwmax+1)/2

      do ns=0,nsmax
         do nr=1,nrmax
            do nn2=1,nphmax2
               do mm2=1,nthmax2
                  do nn=1,nphmax
                     do mm=1,nthmax
                        cpp(mm,nn,mm2,nn2,nr,ns)=0.d0
                     end do
                  end do
               end do
            end do
         end do
      end do

      do ns=0,nsmax
         do nr0=1,nrmax-1

            if(mdlwmd.eq.0) then
               call wmfem_calculate_local(nr0,ns,fml)
            else
               do ml=1,12*nfcmax
                  do mw=1,mwmax
                     fml(mw,ml)=fms(mw,ml,nr0,ns)
                  enddo
               enddo
            endif


            do nr=nr0,nr0+1
               do nfc=1,nfcmax
                  nn=nphnfc(nfc)
                  mm=nthnfc(nfc)
                  do i=1,6
                     mll=6*nthmax*nphmax*(nr-nr0)
     &                  +6*nthmax*(nn-1)+6*(mm-1)+i
                     ml =6*nthmax*nphmax*(nr-1)
     &                  +6*nthmax*(nn-1)+6*(mm-1)+i
                     do nr1=nr0,nr0+1
                        do nfc1=1,nfcmax
                           nn1=nphnfc(nfc1)
                           mm1=nthnfc(nfc1)
                           do i1=1,6
                              ml1=6*nthmax*nphmax*(nr1-1)
     &                           +6*nthmax*(nn1-1)+6*(mm1-1)+i1
                              do nfc2=1,nfcmax
                                 nn2=nphnfc(nfc2)
                                 mm2=nthnfc(nfc2)
                                 nndiff=nnnfc(nfc2)-nnnfc(nfc1)
                                 mmdiff=mmnfc(nfc2)-mmnfc(nfc1)
                                 if(nndiff.lt.0) nndiff=nndiff+nphmax2
                                 if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2

                                 mw=6*nthmax*nphmax*(nr1-nr)
     &                             +6*nthmax*(nn2-nn)+6*(mm2-mm)+(i1-i)

                                 cpp(mm,nn,mmdiff+1,nndiff+1,nr0,ns)
     &                          =cpp(mm,nn,mmdiff+1,nndiff+1,nr0,ns)
     &                          -ci*conjg(fvx(ml))
     &                                *fml(mc+mw,mll)*fvx(ml1)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         end do
      end do

c$$$      do ns=0,nsmax
c$$$         do mmdiff=1,nthmax2
c$$$            csum=0.d0
c$$$            do mm=1,nthmax
c$$$               csum=csum+cpp(mm,1,mmdiff,1,1,ns)
c$$$            enddo
c$$$            write(6,'(A,2I5,1P2E12.4)') 
c$$$     &              'ns,mmdiff,cpp=',
c$$$     &               ns,mmdiff,csum
c$$$         enddo
c$$$      enddo
            
      do ns=0,nsmax
         do mmdiff=1,nthmax2
            csum=0.d0
            do mm=1,nthmax
               csum=csum+cpp(mm,1,mmdiff,1,2,ns)
            enddo
c$$$            write(6,'(A,2I5,1P2E12.4)') 
c$$$     &              'ns,mmdiff,cpp=',
c$$$     &               ns,mmdiff,csum
         enddo
      enddo
            
      csums=0.d0
      do ns=0,nsmax
         csum=0.d0
         do nr=1,nrmax
            do nn1=1,nphmax
            do mm1=1,nthmax
               csum=csum+cpp(mm1,nn1,1,1,nr,ns)
c$$$               if(ns.eq.0) then
c$$$                  write(6,'(A,3I5,1P2E12.4)') 
c$$$     &                 'nr,mm,nn,cpp=',nr,mm1,nn1,cpp(mm1,nn1,1,1,nr,ns)
c$$$               endif
            enddo
            enddo
         enddo
         write(6,'(A,I5,1P2E12.4)') 'NS,PABSP=',ns,csum
         csums=csums+csum
      enddo
      write(6,'(A,5X,1P2E12.4)') '   PABSP=',csums

      if(mdlwmd.ge.1) then
      csums=0.d0
      do ns=0,nsmax
         csum=0.d0
         do ml=1,mlmax
            nr=(ml-1)/(6*nfcmax)+1
            mll=ml-6*nfcmax*(nr-1)
            do mw=max(1,mc-ml+1),min(mwmax,mc-ml+mlmax)
               ml1=ml+mw-mc
               csum=csum-ci*conjg(fvx(ml))*fms(mw,mll,nr,ns)*fvx(ml1)
               if(nr.gt.1) then
                  csum=csum-ci*conjg(fvx(ml))
     &                     *fms(mw,mll+6*nfcmax,nr-1,ns)*fvx(ml1)
               endif
            enddo
         enddo
         write(6,'(A,I5,1P2E12.4)') 'NS,PABSM=',ns,csum
         csums=csums+csum
      enddo
      write(6,'(A,5X,1P2E12.4)') '   PABSM=',csums
      endif

!     ----- calculate antenna impedance -----

      do nn=1,nphmax
      do mm=1,nthmax
         cpa(mm,nn)=0.d0
         do nr=1,nrmax
            ml=6*nthmax*nphmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
            do i=1,6
               cpa(mm,nn)=cpa(mm,nn)-ci*conjg(fvx(ml+i))*fvb(ml+i)
            enddo
         enddo
!         write(6,'(2I5,1P2E12.4)') mm,nn,cpa(mm,nn)
      enddo
      enddo
         csum=0.d0
         do ml=1,mlmax
            csum=csum-ci*conjg(fvx(ml))*fvb(ml)
         enddo
         write(6,'(A,5X,1P2E12.4)') '   PRADM=',csum

      return
      end subroutine wmfem_calculate_power

      end subroutine wmfem_post

