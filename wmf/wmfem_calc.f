!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_local(nr,ns,fml)

      USE wmfem_com
      IMPLICIT NONE
      integer,intent(in):: nr,ns
      complex(8),dimension(mwmax,12*nfcmax),intent(out):: fml
!      complex(8),dimension(3,3,4,nfcmax,nfcmax,4):: fmd
!      complex(8),dimension(3,3,4,nfcmax,nfcmax)::  fmd1,fmd2,fmd3,fmd4
      complex(8),dimension(:,:,:,:,:,:),pointer :: fmd
      complex(8),dimension(:,:,:,:,:),pointer :: fmd1,fmd2,fmd3,fmd4
      complex(8),dimension(nphmax,nthmax,3):: fvb_nr
      real(8):: drho,rkth,rkph,rkth0,rho0,rho1,rho2,rho3,rho4
      integer:: ml,mw,mc,nvmax,i,j,k,inod,nfc,nth,nph,mm,nn,mll
      integer:: nfc1,nfc2,ml1,mw1,mr
      complex(8):: csum,f1,f2,f3,f4,cx,cy,cxd,cyd
      integer:: id_base=1
      real(8):: angl=0.d0

      if(nr.eq.1) then
         rho1=rhoa(nr+1)/9.d0
      else
         rho1=rhoa(nr)+1.D-12
      endif
      rho2=(2.d0*rhoa(nr)+rhoa(nr+1))/3.d0
      rho3=(rhoa(nr)+2.d0*rhoa(nr+1))/3.d0
      rho4=rhoa(nr+1)-1.D-12
      drho=rhoa(nr+1)-rhoa(nr)

      allocate(fmd1(3,3,4,nfcmax,nfcmax),fmd2(3,3,4,nfcmax,nfcmax),
     &         fmd3(3,3,4,nfcmax,nfcmax),fmd4(3,3,4,nfcmax,nfcmax))
      if(ns.eq.0) then
         if(mdlwmf.eq.1) then
            call wmfem_calculate_vacuum(rho1,fmd1)
            call wmfem_calculate_vacuum(rho2,fmd2)
            call wmfem_calculate_vacuum(rho3,fmd3)
            call wmfem_calculate_vacuum(rho4,fmd4)
         else if(mdlwmf.eq.2) then
            call wmfem_calculate_vacuum_c(rho1,fmd1)
            call wmfem_calculate_vacuum_c(rho2,fmd2)
            call wmfem_calculate_vacuum_c(rho3,fmd3)
            call wmfem_calculate_vacuum_c(rho4,fmd4)
         endif
      else
         if(mdlwmf.eq.1) then
            call wmfem_calculate_plasma(rho1,ns,fmd1)
            call wmfem_calculate_plasma(rho2,ns,fmd2)
            call wmfem_calculate_plasma(rho3,ns,fmd3)
            call wmfem_calculate_plasma(rho4,ns,fmd4)
         else if(mdlwmf.eq.2) then
            call wmfem_calculate_plasma_c(rho1,ns,fmd1)
            call wmfem_calculate_plasma_c(rho2,ns,fmd2)
            call wmfem_calculate_plasma_c(rho3,ns,fmd3)
            call wmfem_calculate_plasma_c(rho4,ns,fmd4)
         endif
      endif

! ------ calculate coefficients of basis for profile from four points 

      allocate(fmd(3,3,4,nfcmax,nfcmax,4))
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
     &                    =1.5d0*(-3*f1+ 4*f2- f3)
!     &                    =0.5d0*(-11*f1+18*f2- 9*f3+ 2*f4)
                     fmd(i,j,k,nfc1,nfc2,3)=f4
                     fmd(i,j,k,nfc1,nfc2,4)
     &                    =1.5d0*(f2- 4*f3+ 3*f4)
!     &                    =0.5d0*(- 2*f1+ 9*f2-18*f3+11*f4)
                  enddo
               enddo
            enddo
         enddo
      enddo
            
      if(id_base.eq.0) then
         call fem_hhh(fmd,drho,fml)
      else
         call fem_hqq(fmd,drho,fml)
      endif

      deallocate(fmd,fmd1,fmd2,fmd3,fmd4)

      return
      end subroutine wmfem_calculate_local

!----- calculate coefficint matrix fma (E rho,theta,phi )-----

      subroutine wmfem_calculate_vacuum(rho,fmd)

      USE wmfem_com
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

      USE wmfem_com
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

      USE wmfem_com
      IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),pointer:: fmc
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      complex(8):: cfactor
      real(8):: rkth,rkph,rkth0,rho0
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

         if(rho.ne.0.d0) then
            rkth=mm/(ra*rho)
            rkth0=1.d0/(ra*rho)
         else
            rkth=mm/(ra*1.d-6)
            rkth0=1.d0/(ra*1.d-6)
         endif
         rkph=nn/rr

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

      USE wmfem_com
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),pointer:: fmc
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

      USE wmfem_com
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),pointer:: fmc
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

      use wmfem_com
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

      use wmfem_com
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

      use wmfem_com
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

      use wmfem_com
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

!---- FEM cubic hermit ---

      subroutine fem_hhh(fmd,drho,fml)

      use libfem_mod
      use wmfem_com
      implicit none
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(in):: fmd
      complex(8),dimension(4*6*nfcmax-1,2*6*nfcmax),intent(out):: fml
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nf1,nf2,inod,ml,mw

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif

      mr=6*nfcmax  ! line interval between radial points 
      mc=12*nfcmax ! diagonal column at the top line

      do ml=1,2*6*nfcmax
         do mw=1,4*6*nfcmax-1
            fml(mw,ml)=0.d0
         end do
      end do

         do nf1=1,nfcmax
         do nf2=1,nfcmax
         do j=1,3
         do i=1,3
            ml=6*(nf1-1)+2*(i-1)+1
            mw=mc+6*(nf2-nf1)+2*(j-i)
            do inod=1,4

            fml(mw  ,ml  )=fml(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,5)/drho
            fml(mw+1,ml  )=fml(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,6)
            fml(mw+mr,ml  )=fml(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,7)/drho
            fml(mw+mr+1,ml  )=fml(mw+Mr+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,8)

            fml(mw-1,ml+1)=fml(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,5)
            fml(mw  ,ml+1)=fml(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,6)*drho
            fml(mw+mr-1,ml+1)=fml(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,7)
            fml(mw+mr,ml+1)=fml(mw+mr,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,8)*drho

            fml(mw-mr,ml+mr)=fml(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,5)/drho
            fml(mw-mr+1,ml+mr)=fml(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,6)
            fml(mw  ,ml+mr)=fml(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,7)/drho
            fml(mw+1,ml+mr)=fml(mw+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,8)

            fml(mw-mr-1,ml+mr+1)=fml(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,5)
            fml(mw-mr,ml+mr+1)=fml(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,6)*drho
            fml(mw-1,ml+mr+1)=fml(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,7)
            fml(mw  ,ml+mr+1)=fml(mw  ,ml+mr+1)
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

      subroutine fem_hqq(fmd,drho,fml)

      use libfem_mod
      use wmfem_com
      implicit none
      complex(8),dimension(3,3,4,nfcmax,nfcmax,4),intent(in):: fmd
      complex(8),dimension(4*6*nfcmax-1,2*6*nfcmax),intent(out):: fml
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nf1,nf2,inod,ml,mw

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif

      mr=6*nfcmax  ! line interval between radial points 
      mc=12*nfcmax ! diagonal column at the top line

      do ml=1,2*6*nfcmax
         do mw=1,4*6*nfcmax-1
            fml(mw,ml)=0.d0
         end do
      end do

         do nf1=1,nfcmax
         do nf2=1,nfcmax
         do j=1,3
         do i=1,3
            ml=6*(nf1-1)+2*(i-1)+1
            mw=mc+6*(nf2-nf1)+2*(j-i)
            do inod=1,4

! (1,1) **********
            if(i.eq.1.and.j.eq.1) then

            fml(mw  ,ml  )=fml(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,4,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,1,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,4,4)/drho

            fml(mw+1,ml  )=fml(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,1,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,4,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,4,5)/drho
            fml(mw+mr,ml  )=fml(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,4,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,1,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,4,6)/drho

            fml(mw-1,ml+1)=fml(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,2,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,2,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,5,4)/drho
            fml(mw  ,ml+1)=fml(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,2,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,5,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,2,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,5,5)/drho
            fml(mw+mr-1,ml+1)=fml(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,2,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,2,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,5,6)/drho

            fml(mw-mr,ml+mr)=fml(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,6,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,3,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,6,4)/drho
            fml(mw-mr+1,ml+mr)=fml(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,3,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,6,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,6,5)/drho
            fml(mw  ,ml+mr)=fml(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqq(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqq(inod,6,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqq(inod,3,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqq(inod,6,6)/drho

! (1,*) **********
            else if(i.eq.1) then
            fml(mw  ,ml  )=fml(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,5)/drho
            fml(mw+1,ml  )=fml(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,6)
            fml(mw+mr,ml  )=fml(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,7)/drho
            fml(mw+mr+1,ml  )=fml(mw+mr+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,4,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,1,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,4,8)

            fml(mw-1,ml+1)=fml(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,5)/drho
            fml(mw  ,ml+1)=fml(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,6)
            fml(mw+mr-1,ml+1)=fml(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,7)/drho
            fml(mw+mr,ml+1)=fml(mw+mr,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,2,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,5,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,2,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,5,8)

            fml(mw-mr,ml+mr)=fml(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,5)/drho
            fml(mw-mr+1,ml+mr)=fml(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,6)
            fml(mw  ,ml+mr)=fml(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,7)/drho
            fml(mw+1,ml+mr)=fml(mw+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hqh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hqh(inod,6,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hqh(inod,3,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hqh(inod,6,8)

! (*,1) **********
            else if(j.eq.1) then
            fml(mw  ,ml  )=fml(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,1,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,5,4)/drho
            fml(mw+1,ml  )=fml(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,1,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,5,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,5,5)/drho
            fml(mw+mr,ml  )=fml(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,1,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,5,6)/drho

            fml(mw-1,ml+1)=fml(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,2,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,6,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,2,4)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,6,4)
            fml(mw  ,ml+1)=fml(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,2,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,6,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,2,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,6,5)
            fml(mw+mr-1,ml+1)=fml(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,2,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,6,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,2,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,6,6)

            fml(mw-mr,ml+mr)=fml(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,7,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,3,4)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,7,4)/drho
            fml(mw-mr+1,ml+mr)=fml(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,3,2)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,7,2)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,7,5)/drho
            fml(mw  ,ml+mr)=fml(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,7,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,3,6)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,7,6)/drho

            fml(mw-mr-1,ml+mr+1)=fml(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,4,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,8,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,4,4)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,8,4)
            fml(mw-mr,ml+mr+1)=fml(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,4,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,8,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,4,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,8,5)
            fml(mw-1,ml+mr+1)=fml(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhq(inod,4,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhq(inod,8,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhq(inod,4,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhq(inod,8,6)

! (*,*) **********
            else
            fml(mw  ,ml  )=fml(mw  ,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,5)/drho
            fml(mw+1,ml  )=fml(mw+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,6)
            fml(mw+mr,ml  )=fml(mw+mr,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,7)/drho
            fml(mw+mr+1,ml  )=fml(mw+mr+1,ml  )
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,1,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,5,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,1,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,5,8)

            fml(mw-1,ml+1)=fml(mw-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,5)
            fml(mw  ,ml+1)=fml(mw  ,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,6)*drho
            fml(mw+mr-1,ml+1)=fml(mw+mr-1,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,7)
            fml(mw+mr,ml+1)=fml(mw+mr,ml+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,2,4)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,6,4)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,2,8)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,6,8)*drho

            fml(mw-mr,ml+mr)=fml(mw-mr,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,1)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,1)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,5)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,5)/drho
            fml(mw-mr+1,ml+mr)=fml(mw-mr+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,2)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,2)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,6)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,6)
            fml(mw  ,ml+mr)=fml(mw  ,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,3)*drho
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,3)
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,7)
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,7)/drho
            fml(mw+1,ml+mr)=fml(mw+1,ml+mr)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,3,4)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,7,4)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,3,8)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,7,8)

            fml(mw-mr-1,ml+mr+1)=fml(mw-mr-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,1)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,1)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,5)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,5)
            fml(mw-mr,ml+mr+1)=fml(mw-mr,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,2)*drho**3
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,2)*drho**2
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,6)*drho**2
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,6)*drho
            fml(mw-1,ml+mr+1)=fml(mw-1,ml+mr+1)
     &              +fmd(i,j,1,nf1,nf2,inod)*table_hhh(inod,4,3)*drho**2
     &              +fmd(i,j,2,nf1,nf2,inod)*table_hhh(inod,8,3)*drho
     &              +fmd(i,j,3,nf1,nf2,inod)*table_hhh(inod,4,7)*drho
     &              +fmd(i,j,4,nf1,nf2,inod)*table_hhh(inod,8,7)
            fml(mw  ,ml+mr+1)=fml(mw  ,ml+mr+1)
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
