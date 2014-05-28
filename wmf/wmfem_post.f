!     ***** wmfem post routine *****

      subroutine wmfem_post

      use wmfem_comm
      implicit none
      integer:: ierr

!     ***** calculate bfield *****

      call wmfem_calculate_bfield
      call wmfem_bfld

!     ***** calculate power *****

      call wmfem_calculate_power
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

      USE wmcalc,ONLY: wmfem_inverse_tensor,wmfem_tensors
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

      USE wmcalc,ONLY: wmfem_calculate_local
      implicit none
      integer:: mc,mb,mb1,mb2,mw,nr,ns,mm,nn,i,j
      integer:: ir,nr1,mm1,nn1,mm2,nn2,nndiff,mmdiff
      integer:: i1,ml1,nfc,nfc1,nfc2,ml,mll,nr0
      real(8):: factor
      complex(8):: csum,csums,cfactor
      complex(8),dimension(:,:),ALLOCATABLE:: fma_local

      integer ::ml2

      allocate(fma_local(mbmax,mbmax))

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

            if(mdlwmd.eq.0) then
               call wmfem_calculate_local(nr0,ns,fma_local)
            else
               do mb2=1,mbmax
                  do mb1=1,mbmax
                     fma_local(mb1,mb2)=fma_save(mb1,mb2,nr0,ns)
                  enddo
               enddo
            endif


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
     &                          =cpp(mm,nn,mmdiff+1,nndiff+1,nr0,ns)
     &                          -ci*conjg(fvx(ml))
     &                                *fma_local(mb,mb1)*fvx(ml1)
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

      if(mdlwmd.ge.1) then
      csums=0.d0
      do ns=0,nsmax
         csum=0.d0
         do ml=1,mlmax
            nr=(ml-1)/(6*nfcmax)+1
            mll=ml-6*nfcmax*(nr-1)
            do mw=max(1,mc-ml+1),min(mwmax,mc-ml+mlmax)
               ml1=ml+mw-mc
               csum=csum-ci*conjg(fvx(ml))*fma_save(mw,mll,nr,ns)
     &                                    *fvx(ml1)
               if(nr.gt.1) then
                  csum=csum-ci*conjg(fvx(ml))
     &                     *fma_save(mw,mll+6*nfcmax,nr-1,ns)*fvx(ml1)
               endif
            enddo
         enddo
         write(6,'(A,I5,1P2E12.4)') 'NS,PABSM=',ns,csum
         csums=csums+csum
      enddo
      write(6,'(A,5X,1P2E12.4)') '   PABSM=',csums
      endif

!     ----- calculate antenna impedance -----

      do nn=1,nhhmax
      do mm=1,nthmax
         cpa(mm,nn)=0.d0
         do nr=1,nrmax-1
            ml=6*nthmax*nhhmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
            ml2=8*nthmax*nhhmax*(nr-1)+6*nthmax*(nn-1)+6*(mm-1)
            do i=1,6
               cpa(mm,nn)=cpa(mm,nn)-ci*conjg(fvx(ml+i))*fvb(ml2+i)
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
      
      end subroutine wmfem_post
