!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_local(nr,ns,fml)

      USE wmfem_comm
      IMPLICIT NONE
      integer,intent(in):: nr,ns
      complex(8),dimension(mbmax,mbmax),intent(out):: fml
      complex(8),dimension(:,:,:,:,:),ALLOCATABLE :: fmd1,fmd2
      real(8):: drho,rho1,rho2
      integer:: id_base=1

      drho=rhoa(nr+1)-rhoa(nr)
      rho1=rhoa(nr)+0.25D0*drho
      rho2=rhoa(nr)+0.75D0*drho

      allocate(fmd1(3,3,4,nfcmax,nfcmax))
      allocate(fmd2(3,3,4,nfcmax,nfcmax))

      if(ns.eq.0) then
         if(mdlwmf.eq.1) then
            call wmfem_calculate_vacuum(rho1,fmd1)
            call wmfem_calculate_vacuum(rho2,fmd2)
         else if(mdlwmf.eq.2) then
            call wmfem_calculate_vacuum_c(rho1,fmd2)
            call wmfem_calculate_vacuum_c(rho1,fmd2)
         endif
      else
         if(mdlwmf.eq.1) then
            call wmfem_calculate_plasma(rho1,ns,fmd1)
            call wmfem_calculate_plasma(rho2,ns,fmd2)
         else if(mdlwmf.eq.2) then
            call wmfem_calculate_plasma_c(rho1,ns,fmd1)
            call wmfem_calculate_plasma_c(rho2,ns,fmd2)
         endif
      endif

      call fem_ppq(fmd1,fmd2,drho,fml)

      deallocate(fmd1,fmd2)

      return
      end subroutine wmfem_calculate_local

!----- calculate coefficint matrix fma (E rho,theta,phi )-----

      subroutine wmfem_calculate_vacuum(rho,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8),dimension(3,3,3,3,nfcmax2):: fmv1
      complex(8),dimension(3,3,3,nfcmax2):: fmv2,fmv3
      complex(8),dimension(3,3,nfcmax2):: fmv4
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
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
      end subroutine wmfem_calculate_vacuum

!----- calculate vacuum matrix fms -----

      subroutine wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

      USE wmfem_comm
      IMPLICIT NONE
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
      
      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

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
      end subroutine wmfem_calculate_vacuum_sub

!----- calculate coefficint matrix fmd (E cylindrical +,-,para)-----

      subroutine wmfem_calculate_vacuum_c(rho,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),ALLOCATABLE:: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8):: cfactor
      real(8):: rkth,rkph,rkth0,rho0
      real(8):: rkthrh,rkthrh0,rkth2,rkth02,rkth20
      integer:: i,j,k,nn,mm,nth,nph
      integer:: nfc1,nfc2
      complex(8):: csum,fmd1,fmd2,fmd3,fmd4
      integer:: mm1,mm2,nn1,nn2
      integer:: mmdiff,nndiff,nfcdiff
      real(8):: rr,ra,rb

      allocate(fmc(3,3,4,nfcmax2,nfcmax))
      call get_wmfem_parm1(rr,ra,rb)
      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

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
         mm=mmnfc(nfc2)
         nn=nnnfc(nfc2)

!            rkth=mm/(ra*rho)
!            rkth0=1.d0/(ra*rho)
            rkthrh=mm/ra
            rkthrh0=1.d0/ra
         if(rho.gt.1.d-6) then
            rkth=mm/(ra*rho)
            rkth0=1.d0/(ra*rho)
            rkth2=rkth**2
            rkth20=rkth*rkth0
            rkth02=rkth0**2
         else
!            rkth=mm/(ra*rho)
!            rkth0=1.d0/(ra*rho)
            rkth2=0d0
            rkth20=0d0
            rkth02=0d0
         endif
         rkph=nn/rr

         do nfc1=1,nfcmax2
!            fmc(1,1,1,nfc1,nfc2)= rho*(cfactor-rkph**2-rkth**2)
            fmc(1,1,1,nfc1,nfc2)= rho*(cfactor-rkph**2-rkth2)
!            fmc(1,2,1,nfc1,nfc2)= rho*(-ci*rkth*rkth0)
!            fmc(2,1,1,nfc1,nfc2)= rho*(+ci*rkth*rkth0)
            fmc(1,2,1,nfc1,nfc2)= rho*(-ci*rkth20)
            fmc(2,1,1,nfc1,nfc2)= rho*(+ci*rkth20)
!            fmc(2,2,1,nfc1,nfc2)= rho*(cfactor-rkph**2-rkth0**2)
            fmc(2,2,1,nfc1,nfc2)= rho*(cfactor-rkph**2-rkth02)
            fmc(2,3,1,nfc1,nfc2)= (rkthrh*rkph)
            fmc(3,2,1,nfc1,nfc2)= (rkthrh*rkph)
!            fmc(3,3,1,nfc1,nfc2)= rho*(cfactor-rkth**2)
            fmc(3,3,1,nfc1,nfc2)= rho*(cfactor-rkth2)
                  
            fmc(2,1,2,nfc1,nfc2)= ( ci*rkthrh)/ra
            fmc(2,2,2,nfc1,nfc2)= (-rkthrh0)/ra
            fmc(3,1,2,nfc1,nfc2)= rho*( ci*rkph)/ra
                  
            fmc(1,2,3,nfc1,nfc2)= (-ci*rkthrh)/ra
            fmc(1,3,3,nfc1,nfc2)= rho*(-ci*rkph)/ra
            fmc(2,2,3,nfc1,nfc2)= (  -rkthrh0)/ra
                  
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
         mm2=mmnfc(nfc2)
         nn2=nnnfc(nfc2)
         do nfc1=1,nfcmax
            mm1=mmnfc(nfc1)
            nn1=nnnfc(nfc1)

            nndiff=nn1-nn2
            if(nndiff.lt.0) nndiff=nndiff+nphmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1
!            if(nfcdiff.le.0) nfcdiff=nfcdiff+nfcmax2

            do k=1,4
               do j=1,3
                  do i=1,3
                     fmd(i,j,k,nfc1,nfc2)=fmc(i,j,k,nfcdiff,nfc2)
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(fmc)
      return
      end subroutine wmfem_calculate_vacuum_c

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma(rho,ns,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),ALLOCATABLE:: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      real(8),dimension(3,3,nthmax2,nphmax2) :: gma,muma,dmuma 
      real(8),dimension(nthmax2,nphmax2):: gja

      complex(8):: cfactor
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      integer:: mmadd1,mmadd2,nnadd1,nnadd2
      integer:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

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
            if(nndiff.lt.0) nndiff=nndiff+nphmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1

            nnadd1=(nn1+nn2  )/2-nph0
            if(nnadd1.lt.0)      nnadd1=nnadd1+nphmax
            if(nnadd1.ge.nphmax) nnadd1=nnadd1-nphmax
            nnadd2=(nn1+nn2+1)/2-nph0
            if(nnadd2.lt.0)      nnadd2=nnadd2+nphmax
            if(nnadd2.ge.nphmax) nnadd2=nnadd2-nphmax
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

      end subroutine wmfem_calculate_plasma

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma_c(rho,ns,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),ALLOCATABLE:: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f

      complex(8):: cfactor,cx
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      integer:: mmadd1,mmadd2,nnadd1,nnadd2
      integer:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

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
            if(nndiff.lt.0) nndiff=nndiff+nphmax2
            mmdiff=mm1-mm2
            if(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1

            nnadd1=(nn1+nn2  )/2-nph0
            if(nnadd1.lt.0) nnadd1=nnadd1+nphmax
            nnadd2=(nn1+nn2+1)/2-nph0
            if(nnadd2.lt.0) nnadd2=nnadd2+nphmax
            mmadd1=(mm1+mm2  )/2-nth0
            if(mmadd1.lt.0) mmadd1=mmadd1+nthmax
            mmadd2=(mm1+mm2+1)/2-nth0
            if(mmadd2.lt.0) mmadd2=mmadd2+nthmax

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

      end subroutine wmfem_calculate_plasma_c

!     ****** CALCULATE METRIC AND CONVERSION TENSOR ******

      SUBROUTINE wmfem_tensors(rho,gma,muma,dmuma,gja)

      use wmfem_comm
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

!        ****** CALCULATE ROTATION TENSOR ******

      SUBROUTINE wmfem_rotation_tensor(gm,gj,babs,bsupth,bsupph,mum)

      use wmfem_comm
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
     &     +tc3*(gm(3,3)*gm(1,2)-gm(2,3)*gm(1,3)))
            
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

!     ****** CALCULATE dielectric tensor on a mag surface ******

      SUBROUTINE wmfem_disp_tensor(rho,ns,fmc)

      use wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax2,nfcmax),intent(out):: fmc
      complex(8),dimension(3,3):: fml

      real(8):: th,ph,dth,dph
      integer:: nth,nph,mm,nn,i,j,nfc2,nfc
      complex(8):: cx

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
            CALL wmfem_dielectric(rho,th,ph,mm,nn,ns,fml)
            DO j=1,3
               DO i=1,3
                  fmc(i,j,1,nfc2,nfc)=fml(i,j)
                  fmc(i,j,2,nfc2,nfc)=0.d0
                  fmc(i,j,3,nfc2,nfc)=0.d0
                  fmc(i,j,4,nfc2,nfc)=0.d0
               END DO
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE wmfem_disp_tensor

!     ****** CALCULATE dielectric tensor ******

      SUBROUTINE wmfem_dielectric(rho,th,ph,mm,nn,ns,fml)

      use wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho,th,ph
      integer,intent(in):: mm,nn,ns
      complex(8),dimension(3,3),intent(out):: fml
      complex(8):: cw,ckpara,ckperp
      complex(8):: ckppf,ckpps
      real(8):: babs,bsupth,bsupph

      cw=2.d0*pi*crf*1.d6
      CALL wmfem_magnetic(rho,th,ph,babs,bsupth,bsupph)
      ckpara=mm*bsupth/babs+nn*bsupph/babs
      ckperp=(0.d0,0.d0)

c$$$      if(abs(th).lt.1.d-4) then
c$$$         CALL DPCOLD_RKPERP(cw,ckpara,ckppf,ckpps)
c$$$         write(6,'(1P6E12.4)') rho,real(ckpara),ckppf,ckpps 
c$$$      endif

c$$$         if(ns.eq.1.and.abs(th).lt.0.1) 
c$$$     &    write(6,'(a,2i5,1p5e12.4)') 'm,n,kpara:',
c$$$     &         mm,nn,dble(ckpara),rho,babs,bsupth,bsupph

      CALL dpcalc(cw,ckpara,ckperp,rho,babs,ns,fml)

      return

      END SUBROUTINE wmfem_dielectric
