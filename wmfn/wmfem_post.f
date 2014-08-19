!     ***** wmfem post routine *****

      subroutine wmfem_post

      use wmfem_comm
      implicit none
      integer:: ierr

!     ***** calculate bfield *****

       !!!!????????!!!!!!!
      print *,"post_seki1"
      call wmfem_calculate_bfield
      print *,"post_seki2"
      call wmfem_bfld
       !!!!????????!!!!!!!

!     ***** calculate power *****

      print *,"post_seki3"
      call wmfem_calculate_power
      print *,"post_seki4"
      call wmfem_pabs

      return

      contains

!     ***** Calculate magnetic field *****

      subroutine wmfem_calculate_bfield

      use wmfem_comm
      implicit none
      integer:: nr,nth,nhh
      complex(8),dimension(3,nthmax,nhhmax):: cbfl

      do nr=1,nrmax
         call wmfem_cbf(nr,cbfl)
         do nhh=1,nhhmax
            do nth=1,nthmax
               cbf(1,nth,nhh,nr)=cbfl(1,nth,nhh)
               cbf(2,nth,nhh,nr)=cbfl(2,nth,nhh)
               cbf(3,nth,nhh,nr)=cbfl(3,nth,nhh)
            enddo
         enddo
      enddo

      return
      end subroutine wmfem_calculate_bfield

!----- calculate local wave magnetic field -----

      subroutine wmfem_cbf(nr,cbf)

      implicit none
      integer,intent(in):: nr
      complex(8),dimension(3,nthmax,nhhmax):: cbf

      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      real(8),dimension(nthmax2,nhhmax2):: gja

      complex(8),dimension(3,3,3):: cq
      complex(8),dimension(3,3):: cp
      complex(8),dimension(3,3,3,nfcmax2):: cqq
      complex(8),dimension(3,3,nfcmax2):: cpp
      complex(8),dimension(nthmax2,nhhmax2):: fv1,fv1f
      integer:: i,j,k,l,nthm,nthp,nhhm,nhhp
      integer:: nfc2,nhh,nth,imn,ml
      integer:: nfc1,nhh1,nth1
      integer:: nn1,mm1,nn2,mm2,mmdiff,nndiff,nfcdiff
      real(8):: rho,dph,dth,gj
      complex(8):: cfactor
      
      cfactor=vc/(2.d0*pi*crf*1.d6)

      rho=rhoa(nr)

!      call wmfem_tensors(rho,gma,muma,dmuma,gja)
      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

!     ----- calculation rot coefficients cpp and cqq -----

      do nfc2=1,nfcmax2
         nth=nthnfc2(nfc2)
         nhh=nhhnfc2(nfc2)
         if(nhh.eq.1) then
            nhhm=nhhmax2
         else
            nhhm=nhh-1
         endif
         if(nhh.eq.nhhmax2) then
            nhhp=1
         else
            nhhp=nhh+1
         endif
         dph=2*pi/nhhmax2

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
         gj=gja(nth,nhh)

         do j=1,3
            cq(1,j,1)=((muma(3,j,nthp,nhh)
     &                 -muma(3,j,nthm,nhh))/dth
     &                -(muma(2,j,nth,nhhp)
     &                 -muma(2,j,nth,nhhm))/dph )/gj
            cq(1,j,2)=+ci*muma(3,j,nth,nhh)/gj
            cq(1,j,3)=-ci*muma(2,j,nth,nhh)/gj

            cq(2,j,1)=((muma(1,j,nth,nhhp)
     &                 -muma(1,j,nth,nhhm))/dph
     &                 -dmuma(3,j,nth,nhh))/gj
            cq(2,j,2)=0.d0
            cq(2,j,3)=+ci*muma(1,j,nth,nhh)/gj

            cq(3,j,1)=(dmuma(2,j,nth,nhh)
     &               -(muma(1,j,nthp,nhh)
     &                -muma(1,j,nthm,nhh))/dth )/gj
            cq(3,j,2)=-ci*muma(1,j,nth,nhh)/gj
            cq(3,j,3)=0.d0

            cp(1,j)=0.d0
            cp(2,j)=-muma(3,j,nth,nhh)/gj
            cp(3,j)= muma(2,j,nth,nhh)/gj
         enddo

         do imn=1,3
            do i=1,3
               do j=1,3
                  cqq(i,j,imn,nfc2)=0.d0
                  do l=1,3
                     cqq(i,j,imn,nfc2)=cqq(i,j,imn,nfc2)
     &                 +gma(i,l,nth,nhh)*cq(l,j,imn)*gj
                  enddo
               enddo
            enddo
         enddo
         do i=1,3
            do j=1,3
               cpp(i,j,nfc2)=0.d0
               do l=1,3
                  cpp(i,j,nfc2)=cpp(i,j,nfc2)
     &                 +gma(i,l,nth,nhh)*cp(l,j)*gj
               enddo
            enddo
         enddo
      enddo

!     ----- FFT of cqq and cpp -----
            
      do i=1,3
         do j=1,3
            do imn=1,3
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nhh=nhhnfc2(nfc2)
                  fv1(nth,nhh)=cqq(i,j,imn,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nhh=nhhnfc2(nfc2)
                  cqq(i,j,imn,nfc2)=fv1f(nth,nhh)
               enddo
            enddo
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nhh=nhhnfc2(nfc2)
               fv1(nth,nhh)=cpp(i,j,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nhh=nhhnfc2(nfc2)
               cpp(i,j,nfc2)=fv1f(nth,nhh)
            enddo
         enddo
      enddo
      
!     ----- calculated cbf -----
      
      do nfc1=1,nfcmax
         nhh1=nhhnfc(nfc1)
         nth1=nthnfc(nfc1)
         nn1=nnnfc(nfc1)
         mm1=mmnfc(nfc1)
         do nfc2=1,nfcmax
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)

            nndiff=nn1-nn2
            if(nndiff.lt.0) nndiff=nndiff+nhhmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1
            ml=6*nthmax*nhhmax*(nr-1)+(nfc2-1)

!          !!!!! on axis should be considered !!!!!

            do i=1,3
               cbf(i,nth1,nhh1)=0.d0
               do j=1,3
                  cbf(i,nth1,nhh1)=cbf(i,nth1,nhh1)
     &                 +cfactor*(
     &                  (cqq(i,j,1,nfcdiff)
     &                  +cqq(i,j,2,nfcdiff)*mm2
     &                  +cqq(i,j,3,nfcdiff)*nn2)*cef(j,nth1,nhh1,nr)
     &                  +cpp(i,j,  nfcdiff)     *cdef(j,nth1,nhh1,nr))
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_cbf

!     ***** Calculate absorbed power *****

      subroutine wmfem_calculate_power

      implicit none
      integer:: mc,mb,mb1,mb2,mw,nr,ns,mm,nn,i,j
      integer:: ir,nr1,mm1,nn1,mm2,nn2,nndiff,mmdiff
      integer:: i1,ml1,nfc,nfc1,nfc2,ml,mll,nr0
      real(8):: factor
      complex(8):: csum,csums,cfactor
      complex(8),dimension(:,:),ALLOCATABLE:: fma_local

      integer ::ml2

!      allocate(fma_local(mbmax,mbmax))
      allocate(fma_local(mwmax,mbmax))

      mc=(mwmax+1)/2

      do ns=0,nsmax
         do nr=1,nrmax
            do nn2=1,nhhmax2
               do mm2=1,nthmax2
                  do nn=1,nhhmax
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

!            if(mdlwmd.eq.0) then
               call wmfem_calculate_local_ef(nr0,ns,fma_local)
!            else
!               do mb2=1,mbmax
!                  do mb1=1,mbmax
!                     fma_local(mb1,mb2)=fma_save(mb1,mb2,nr0,ns)
!                  enddo
!               enddo
!            endif


       !!!!????????!!!!!!!
            do nfc=1,nfcmax
               nn=nhhnfc(nfc)
               mm=nthnfc(nfc)
               do i=1,6
                  mb=nfcmax*(i-1)+nfc
                  ml=6*nfcmax*(nr0-1)+mb
                  do nfc1=1,nfcmax
                     nn1=nhhnfc(nfc1)
                     mm1=nthnfc(nfc1)
                     do i1=1,6
                        mb1=nfcmax*(i1-1)+nfc
                        ml1=6*nfcmax*(nr0-1)+mb1
                        do nfc2=1,nfcmax
                           nn2=nhhnfc(nfc2)
                           mm2=nthnfc(nfc2)
                           nndiff=nnnfc(nfc2)-nnnfc(nfc1)
                           mmdiff=mmnfc(nfc2)-mmnfc(nfc1)
                           if(nndiff.lt.0) nndiff=nndiff+nhhmax2
                           if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
                           cpp(mm,nn,mmdiff+1,nndiff+1,nr0,ns)
      !!!!????????!!!!!!!
     &                          =cpp(mm,nn,mmdiff+1,nndiff+1,nr0,ns)
     &                          -ci*conjg(fvx_ef(ml))
     &                             *fma_local(mc+mb-mb1,mb1)*fvx_ef(ml1)
       !!!!????????!!!!!!!
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
       !!!!????????!!!!!!!

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
            do nn1=1,nhhmax
            do mm1=1,nthmax
               csum=csum+cpp(mm1,nn1,1,1,nr,ns)
            enddo
            enddo
         enddo
         write(6,'(A,I5,1P2E12.4)') 'NS,PABSP=',ns,csum
         csums=csums+csum
      enddo
      write(6,'(A,5X,1P2E12.4)') '   PABSP=',csums

!      if(mdlwmd.ge.1) then
!      csums=0.d0
!      do ns=0,nsmax
!         csum=0.d0
!         do ml=1,mlmax
!            nr=(ml-1)/(6*nfcmax)+1
!            mll=ml-6*nfcmax*(nr-1)
!            do mw=max(1,mc-ml+1),min(mwmax,mc-ml+mlmax)
!               ml1=ml+mw-mc
!               csum=csum-ci*conjg(fvx(ml))*fma_save(mw,mll,nr,ns)
!     &                                    *fvx(ml1)
!               if(nr.gt.1) then
!                  csum=csum-ci*conjg(fvx(ml))
!     &                     *fma_save(mw,mll+6*nfcmax,nr-1,ns)*fvx(ml1)
!               endif
!            enddo
!         enddo
!         write(6,'(A,I5,1P2E12.4)') 'NS,PABSM=',ns,csum
!         csums=csums+csum
!      enddo
!      write(6,'(A,5X,1P2E12.4)') '   PABSM=',csums
!      endif

!     ----- calculate antenna impedance -----

       !!!!????????!!!!!!!
      do nn=1,nhhmax
      do mm=1,nthmax
         cpa(mm,nn)=0.d0
         do nr=1,nrmax-1
            ml=6*nthmax*nhhmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
            ml2=8*nthmax*nhhmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
            do i=1,6
               cpa(mm,nn)=cpa(mm,nn)-ci*conjg(fvx_ef(ml+i))*fvb(ml2+i)
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

      deallocate(fma_local)
       !!!!????????!!!!!!!

      return
      end subroutine wmfem_calculate_power
      
      subroutine wmfem_calculate_local_ef(nr,ns,fml)

      use wmfem_comm
      IMPLICIT NONE
      integer,intent(in):: nr,ns
      complex(8),dimension(mwmax,mbmax),intent(out):: fml
      complex(8),dimension(:,:,:,:,:,:),allocatable :: fmd
      complex(8),dimension(:,:,:,:,:,:),allocatable :: fmd_p
      complex(8):: moment0,moment1,moment2,moment3
      real(8):: rho0,xi,yi
      real(8),parameter:: m0=4.D0
      real(8),parameter:: m2=5.D0/16.D0
      real(8),parameter:: m4=41.D0/1024.D0
      real(8),parameter:: m6=365.D0/65536.D0
      integer ::inod,nfc2,nfc1
      integer :: i,j,k
      real(8) ::drho

      allocate(fmd_p(4,4,4,nfcmax,nfcmax,4))
      fmd_p=0d0
      drho=rhoa(nr+1)-rhoa(nr)
      do inod=1,4
         select case(inod)
           case(1)
              rho0=(7.D0*rhoa(nr)+     rhoa(nr+1))/8.D0
           case(2)
              rho0=(5.D0*rhoa(nr)+3.D0*rhoa(nr+1))/8.D0
           case(3)
              rho0=(3.D0*rhoa(nr)+5.D0*rhoa(nr+1))/8.D0
           case(4)
              rho0=(     rhoa(nr)+7.D0*rhoa(nr+1))/8.D0
        end select

        if(ns.eq.0) then
!          if(mdlwmf.eq.1) then
            call wmfem_calculate_vacuum_ef(rho0,fmd_p(:,:,:,:,:,inod))
!          else if(mdlwmf.eq.2) then
!            call wmfem_calculate_vacuum_c(rho0,fmd_p(:,:,:,:,:,inod))
!          endif
        else
!          if(mdlwmf.eq.1) then
           call wmfem_calculate_plasma_ef(rho0,ns,fmd_p(:,:,:,:,:,inod))
!          else if(mdlwmf.eq.2) then
!            call wmfem_calculate_plasma_c(rho0,ns,fmd_p(:,:,:,:,:,inod))
!          endif
        endif
      enddo
! ------ calculate coefficients of basis for profile from four points 
      allocate(fmd(4,4,4,nfcmax,nfcmax,4))
        fmd=0d0
        do nfc2=1,nfcmax
        do nfc1=1,nfcmax
        do k=1,4
        do j=1,4
           do i=1,4
              moment0=0.D0
              moment1=0.D0
              moment2=0.D0
              moment3=0.D0
              do inod=1,4
                 xi=0.125D0*(2*inod-1)
                 yi=xi-0.5D0
                 moment0=moment0+      fmd_p(i,j,k,nfc1,nfc2,inod)
                 moment1=moment1+yi   *fmd_p(i,j,k,nfc1,nfc2,inod)
                 moment2=moment2+yi**2*fmd_p(i,j,k,nfc1,nfc2,inod)
                 moment3=moment3+yi**3*fmd_p(i,j,k,nfc1,nfc2,inod)
              enddo
                fmd(i,j,k,nfc1,nfc2,1)=
     &                     (m4*moment0-m2*moment2)/(m0*m4-m2**2)
                fmd(i,j,k,nfc1,nfc2,2)=
     &                     (m6*moment1-m4*moment3)/(m2*m6-m4**2)
                fmd(i,j,k,nfc1,nfc2,3)=
     &                     (m2*moment0-m0*moment2)/(m2**2-m0*m4)*2.D0
                fmd(i,j,k,nfc1,nfc2,4)=
     &                     (m4*moment1-m2*moment3)/(m4**2-m2*m6)*6.D0
           enddo
        enddo
        enddo
        enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!not change
       call fem_hhg(fmd,drho,fml)
       deallocate(fmd,fmd_p)
       return
      end subroutine wmfem_calculate_local_ef

      subroutine wmfem_calculate_vacuum_ef(rho,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(4,4,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(nthmax2,nhhmax2):: fv1,fv1f
      complex(8),dimension(3,3,3,3,nfcmax2):: fmv1
      complex(8),dimension(3,3,3,nfcmax2):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2):: fmv4
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
      integer:: imn,imn1,imn2,nfcdiff

      call wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

      fmd=0d0

      do i=1,3
         do j=1,3
            do imn1=1,3
               do imn2=1,3
                  do nfc2=1,nfcmax2
                     nth=nthnfc2(nfc2)
                     nph=nhhnfc2(nfc2)
                     fv1(nth,nph)=fmv1(i,j,imn1,imn2,nfc2)
                  enddo
                  call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
                  do nfc2=1,nfcmax2
                     nth=nthnfc2(nfc2)
                     nph=nhhnfc2(nfc2)
                     fmv1(i,j,imn1,imn2,nfc2)=fv1f(nth,nph)
                  enddo
               enddo
            enddo    
            do imn1=1,3
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fv1(nth,nph)=fmv2(i,j,imn1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fmv2(i,j,imn1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do imn1=1,3
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fv1(nth,nph)=fmv3(i,j,imn1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fmv3(i,j,imn1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nhhnfc2(nfc2)
               fv1(nth,nph)=fmv4(i,j,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nhhnfc2(nfc2)
               fmv4(i,j,nfc2)=fv1f(nth,nph)
            enddo
         enddo
      enddo

      do nfc1=1,nfcmax          ! Fit to fmd and adjust m and n
         nn1=nnnfc(nfc1)
         mm1=mmnfc(nfc1)
         do nfc2=1,nfcmax
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)

            nndiff=nn1-nn2
            if(nndiff.lt.0) nndiff=nndiff+nhhmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1
!            if(nfcdiff.le.0) nfcdiff=nfcdiff+nfcmax2

            do j=1,3
               do i=1,3
                  fmd(i,j,1,nfc1,nfc2)
     &                 =fmv1(i,j,1,1,nfcdiff)
     &                 +fmv1(i,j,2,1,nfcdiff)*mm1
     &                 +fmv1(i,j,3,1,nfcdiff)*nn1
     &                 +fmv1(i,j,1,2,nfcdiff)    *mm2
     &                 +fmv1(i,j,2,2,nfcdiff)*mm1*mm2
     &                 +fmv1(i,j,3,2,nfcdiff)*nn1*mm2
     &                 +fmv1(i,j,1,3,nfcdiff)    *nn2
     &                 +fmv1(i,j,2,3,nfcdiff)*mm1*nn2
     &                 +fmv1(i,j,3,3,nfcdiff)*nn1*nn2
                  fmd(i,j,2,nfc1,nfc2)
     &                 =fmv2(i,j,1,nfcdiff)
     &                 +fmv2(i,j,2,nfcdiff)*mm1
     &                 +fmv2(i,j,3,nfcdiff)*nn1
                  fmd(i,j,3,nfc1,nfc2)
     &                 =fmv3(i,j,1,nfcdiff)
     &                 +fmv3(i,j,2,nfcdiff)    *mm2
     &                 +fmv3(i,j,3,nfcdiff)    *nn2
                  fmd(i,j,4,nfc1,nfc2)
     &                 =fmv4(i,j,nfcdiff)

c$$$                  write(6,'(A,2I5,1P4E12.4)')
c$$$     &                 'fmd:',i,j,fmd(i,j,1,nfc1,nfc2),
c$$$     &                            fmd(i,j,2,nfc1,nfc2)
c$$$                  write(6,'(A,10X,1P4E12.4)')
c$$$     &                 'fmd:',    fmd(i,j,3,nfc1,nfc2),
c$$$     &                            fmd(i,j,4,nfc1,nfc2)
            enddo
            enddo
         enddo
      enddo
      return
      end subroutine wmfem_calculate_vacuum_ef

!----- calculate vacuum matrix fms -----

      subroutine wmfem_calculate_vacuum_ef_sub(rho,fmv1,fmv2,fmv3,fmv4)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(3,3,3,3,nfcmax2),intent(out):: fmv1
      complex(8),dimension(3,3,3,nfcmax2),intent(out):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2),intent(out):: fmv4

      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      real(8),dimension(nthmax2,nhhmax2):: gja

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
      
      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

!      call wmfem_tensors(rho,gma,muma,dmuma,gja)
      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

      fmv1=0d0
      fmv2=0d0
      fmv3=0d0
      fmv4=0d0
      

      do nfc2=1,nfcmax2
         nth=nthnfc2(nfc2)
         nph=nhhnfc2(nfc2)

         if(nph.eq.1) then
            nphm=nhhmax2
         else
            nphm=nph-1
         endif
         if(nph.eq.nhhmax2) then
            nphp=1
         else
            nphp=nph+1
         endif
         dph=2*pi/nhhmax2

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
                           csum1=csum1-conjg(cq(k,i,imn1))
     &                                *gma(k,l,nth,nph)
     &                                *cq(l,j,imn2) *gj
                        enddo
                     enddo
                     if(i.eq.j.and.imn1.eq.1.and.imn2.eq.1) then
                        fmv1(i,j,imn1,imn2,nfc2)
     &                       =csum1+cfactor*gj
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
                        csum2=csum2-conjg(cp(k,i))
     &                             *gma(k,l,nth,nph)
     &                             *cq(l,j,imn1)*gj
                        csum3=csum3-conjg(cq(k,i,imn1))
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
                     csum4=csum4-conjg(cp(k,i))
     &                          *gma(k,l,nth,nph)
     &                          *cp(l,j)*gj
                  enddo
               enddo
               fmv4(i,j,nfc2)=csum4
            enddo
         enddo

      enddo
      return
      end subroutine wmfem_calculate_vacuum_ef_sub

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma_ef(rho,ns,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(4,4,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),allocatable:: fmc
      complex(8),dimension(nthmax2,nhhmax2):: fv1,fv1f
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      real(8),dimension(nthmax2,nhhmax2):: gja

      complex(8):: cfactor
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      integer:: mmadd1,mmadd2,nnadd1,nnadd2
      integer:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

!      call wmfem_tensors(rho,gma,muma,dmuma,gja)
      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

      call wmfem_disp_tensor(rho,ns,fmc)
      
      fmd=0d0
      fmc(:,:,2:4,:,:)=0d0

      do j=1,3
         do i=1,3
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=cfactor*fmc(i,j,1,nfc1,nfc2)
     &                                *gja(nth,nph)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc(i,j,1,nfc1,nfc2)=fv1f(nth,nph)
               enddo
            enddo    
         enddo
      enddo    

      do nfc2=1,nfcmax
         do nfc1=1,nfcmax2
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
         mm2=mmnfc(nfc2)
         nn2=nnnfc(nfc2)
         do nfc1=1,nfcmax
            mm1=mmnfc(nfc1)
            nn1=nnnfc(nfc1)

            nndiff=nn1-nn2
            if(nndiff.lt.0) nndiff=nndiff+nhhmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1

            nnadd1=(nn1+nn2  )/2-nph0
            if(nnadd1.lt.0)      nnadd1=nnadd1+nhhmax
            if(nnadd1.ge.nhhmax) nnadd1=nnadd1-nhhmax
            nnadd2=(nn1+nn2+1)/2-nph0
            if(nnadd2.lt.0)      nnadd2=nnadd2+nhhmax
            if(nnadd2.ge.nhhmax) nnadd2=nnadd2-nhhmax
            mmadd1=(mm1+mm2  )/2-nth0
            if(mmadd1.lt.0)      mmadd1=mmadd1+nthmax
            if(mmadd1.ge.nthmax) mmadd1=mmadd1-nthmax
            mmadd2=(mm1+mm2+1)/2-nth0
            if(mmadd2.lt.0)      mmadd2=mmadd2+nthmax
            if(mmadd2.ge.nthmax) mmadd2=mmadd2-nthmax

            nfcadd1=nthmax*nnadd1+mmadd1+1
            nfcadd2=nthmax*nnadd1+mmadd2+1
            nfcadd3=nthmax*nnadd2+mmadd1+1
            nfcadd4=nthmax*nnadd2+mmadd2+1

            do k=1,4
               do j=1,3
                  do i=1,3
                     fmd(i,j,k,nfc1,nfc2)
     &                    =0.25d0*fmc(i,j,k,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc(i,j,k,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc(i,j,k,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc(i,j,k,nfcdiff,nfcadd4)
                  enddo
               enddo
            enddo
         enddo
      enddo

      deallocate(fmc)

      return

      end subroutine wmfem_calculate_plasma_ef

      end subroutine wmfem_post
