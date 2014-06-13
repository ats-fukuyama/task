!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_local(nr,ns,fml)

      use wmfem_comm
      IMPLICIT NONE
      integer,intent(in):: nr,ns
      complex(8),dimension(mwmax,mbmax),intent(out):: fml
      complex(8),dimension(:,:,:,:,:,:),pointer :: fmd
      complex(8),dimension(:,:,:,:,:,:),pointer :: fmd_p
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
            call wmfem_calculate_vacuum(rho0,fmd_p(:,:,:,:,:,inod))
!          else if(mdlwmf.eq.2) then
!            call wmfem_calculate_vacuum_c(rho0,fmd_p(:,:,:,:,:,inod))
!          endif
        else
!          if(mdlwmf.eq.1) then
            call wmfem_calculate_plasma(rho0,ns,fmd_p(:,:,:,:,:,inod))
!          else if(mdlwmf.eq.2) then
!            call wmfem_calculate_plasma_c(rho0,ns,fmd_p(:,:,:,:,:,inod))
!          endif
        endif
      enddo

      
! ------ calculate coefficients of basis for profile from four points 
      allocate(fmd(4,4,4,nfcmax,nfcmax,4))
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
!                 print *, fmd_p(i,j,k,nfc1,nfc2,inod)
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
       call fem_hhg(fmd,drho,fml)
       deallocate(fmd,fmd_p)
       return
      end subroutine wmfem_calculate_local


!----- calculate coefficint matrix fma (A  rho,theta,phi )-----

      subroutine wmfem_calculate_vacuum(rho,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(4,4,4,nfcmax,nfcmax),intent(out):: fmd
      complex(8),dimension(4,4,4,nfcmax,nfcmax):: fmd_plasma
      complex(8),dimension(nthmax2,nhhmax2):: fv1,fv1f
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

!      call wmfem_calculate_plasma(rho,0,fmd)
      call wmfem_calculate_plasma(rho,0,fmd_plasma)
      fmd(:,4,:,:,:)=fmd_plasma(:,4,:,:,:)
      fmd(4,1,:,:,:)=fmd_plasma(4,1,:,:,:)
      fmd(4,2,:,:,:)=fmd_plasma(4,2,:,:,:)
      fmd(4,3,:,:,:)=fmd_plasma(4,3,:,:,:)

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
      end subroutine wmfem_calculate_vacuum

!----- calculate vacuum matrix fms -----

      subroutine wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

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
!      real(8),dimension(3,3,nthmax2,nhhmax2) :: muminva 
      real(8),dimension(3,3) :: muminva 

      complex(8),dimension(3,3,3):: cq
      complex(8),dimension(3,3):: cp
      complex(8),dimension(3,3):: cqd
      complex(8),dimension(3):: cpd
      integer:: idiv
      integer:: i,j,k,l,nthm,nthp,nphm,nphp
      integer:: imn1,imn2
      integer:: nph,nth
      real(8):: dph,dth
      complex(8):: csum1,csum2,csum3,csum4
      complex(8):: csumd1,csumd2,csumd3,csumd4
      integer:: nrl,nfc2
      real(8):: drhob,muma3b,muma2b,drhoa,muma3a,muma2a,muma30,muma20
      real(8):: muma3d,muma2d,gj

      complex(8):: cfactor
      complex(8):: cfactor_div

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

      cfactor_div=1d0 !(2.d0*pi*crf*1.d6)
!      cfactor_div=0d0 !(2.d0*pi*crf*1.d6)**2

!!!!

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)
!      call wmfem_inverse_tensor(muma,muminva)



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
     &                 -muma(2,j,nth,nphm))/dph ) /gj
            cq(1,j,2)=+ci*muma(3,j,nth,nph) /gj
            cq(1,j,3)=-ci*muma(2,j,nth,nph) /gj

            cq(2,j,1)=((muma(1,j,nth,nphp)
     &                 -muma(1,j,nth,nphm))/dph
     &                 -dmuma(3,j,nth,nph)) /gj
            cq(2,j,2)=0.d0
            cq(2,j,3)=+ci*muma(1,j,nth,nph) /gj

            cq(3,j,1)=(dmuma(2,j,nth,nph)
     &               -(muma(1,j,nthp,nph)
     &                -muma(1,j,nthm,nph))/dth ) /gj
            cq(3,j,2)=-ci*muma(1,j,nth,nph) /gj
            cq(3,j,3)=0.d0

            cp(1,j)=0.d0
            cp(2,j)=-muma(3,j,nth,nph) /gj
            cp(3,j)= muma(2,j,nth,nph) /gj

            
            cqd(j,1)= (dgjgmuma(j,nth,nph)
     &                +( gja(nthp,nph)*gmuma(2,j,nthp,nph)
     &                  - gja(nthm,nph)*gmuma(2,j,nthm,nph))/dth
     &                +( gja(nth,nphp)*gmuma(3,j,nth,nphp)
     &                  - gja(nth,nphm)*gmuma(3,j,nth,nphm))/dph
     &                )  /gj
            cqd(j,2)= +ci*gmuma(2,j,nth,nph)
            cqd(j,3)= +ci*gmuma(3,j,nth,nph)
                       
            cpd(j)=  gmuma(1,j,nth,nph)
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
     &                                *cq(l,j,imn2)*gj
                        enddo
                     enddo
   
                     csumd1=-conjg(cqd(i,imn1))*cqd(j,imn2) *gj

 
                     if(i.eq.j.and.imn1.eq.1.and.imn2.eq.1) then
                        fmv1(i,j,imn1,imn2,nfc2)
     &                       =csum1+csumd1*cfactor_div +cfactor*gj
!     &                       =csum1+csumd1 +cfactor*gj
!     &                       =csum1-csumd1 +cfactor*gj
!     &                       =csum1 +cfactor*gj
                     else
                        fmv1(i,j,imn1,imn2,nfc2)=csum1 
     &                   + csumd1*cfactor_div
!                        fmv1(i,j,imn1,imn2,nfc2)=csum1 +csumd1
!                        fmv1(i,j,imn1,imn2,nfc2)=csum1 -csumd1
!                        fmv1(i,j,imn1,imn2,nfc2)=csum1
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
     &                             *cq(l,j,imn1)  *gj
                        csum3=csum3-conjg(cq(k,i,imn1))
     &                             *gma(k,l,nth,nph)
     &                             *cp(l,j)  *gj
                     enddo
                  enddo

                  csumd2=-conjg(cpd(i))*cqd(j,imn1)  *gj
                  csumd3=-conjg(cqd(i,imn1))*cpd(j)  *gj

                  fmv2(i,j,imn1,nfc2)=csum2 + csumd2*cfactor_div
                  fmv3(i,j,imn1,nfc2)=csum3 + csumd3*cfactor_div
!                  fmv2(i,j,imn1,nfc2)=csum2 + csumd2
!                  fmv3(i,j,imn1,nfc2)=csum3 + csumd3
!                  fmv2(i,j,imn1,nfc2)=csum2 - csumd2
!                  fmv3(i,j,imn1,nfc2)=csum3 - csumd3
!                  fmv2(i,j,imn1,nfc2)=csum2
!                  fmv3(i,j,imn1,nfc2)=csum3
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
     &                          *cp(l,j)  *gj
                  enddo
               enddo
               csumd4=-conjg(cpd(i))*cpd(j)  *gj

               fmv4(i,j,nfc2)=csum4 + csumd4* cfactor_div
!               fmv4(i,j,nfc2)=csum4 - csumd4
!               fmv4(i,j,nfc2)=csum4
            enddo
         enddo

      enddo
      return
      end subroutine wmfem_calculate_vacuum_sub

!----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_plasma(rho,ns,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(4,4,4,nfcmax,nfcmax),intent(out):: fmd
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),pointer:: fmc
      complex(8),dimension(nthmax2,nhhmax2):: fv1,fv1f

      complex(8),dimension(3,3,nfcmax2,nfcmax)::
     &                                   fmc41,fmc4a1,fmca41
      complex(8),dimension(3,nfcmax2,nfcmax)::
     &                                    fmc42,fmc43,fmc4a2,fmca43
      complex(8),dimension(nfcmax2,nfcmax):: fmc44
      
      complex(8),dimension(3,3) ::fmc41_d,fmc4a1_d,fmca41_d
      complex(8),dimension(3) ::fmc42_d,fmc43_d,fmc4a2_d,fmca43_d
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(nthmax2,nhhmax2):: gja
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
!      real(8),dimension(3,3,nthmax2,nhhmax2) :: muminva 
!      real(8),dimension(3,3) :: muminva 

      complex(8):: cfactor
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      integer:: mmadd1,mmadd2,nnadd1,nnadd2
      integer:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

      if (ns == 0)then
        fmc(:,:,:,:,:)=0d0
        fmc(1,1,:,:,:)=1d0
        fmc(2,2,:,:,:)=1d0
        fmc(3,3,:,:,:)=1d0
      else
        call wmfem_disp_tensor(rho,ns,fmc)
      endif

      call wmfem_calculate_plasma_sub(rho,ns,
     & fmc41,fmc42,fmc43,fmc44,fmca41,fmca43,fmc4a1,fmc4a2)
      do i =1,3 
         do j = 1,3
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=cfactor*fmc(i,j,1,nfc1,nfc2)*gja(nth,nph)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc(i,j,1,nfc1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=fmc41(i,j,nfc1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc41(i,j,nfc1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=fmc4a1(i,j,nfc1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc4a1(i,j,nfc1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=fmca41(i,j,nfc1,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               do nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmca41(i,j,nfc1,nfc2)=fv1f(nth,nph)
               enddo
            enddo
         enddo
         do nfc2=1,nfcmax
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmc42(i,nfc1,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmc42(i,nfc1,nfc2)=fv1f(nth,nph)
            enddo
         enddo
         do nfc2=1,nfcmax
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmc43(i,nfc1,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmc43(i,nfc1,nfc2)=fv1f(nth,nph)
            enddo
         enddo
         do nfc2=1,nfcmax
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmc4a2(i,nfc1,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmc4a2(i,nfc1,nfc2)=fv1f(nth,nph)
            enddo
         enddo
         do nfc2=1,nfcmax
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmca43(i,nfc1,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            do nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmca43(i,nfc1,nfc2)=fv1f(nth,nph)
            enddo
         enddo
      enddo    
      do nfc2=1,nfcmax
         do nfc1=1,nfcmax2
            nth=nthnfc2(nfc1)
            nph=nhhnfc2(nfc1)
            fv1(nth,nph)=fmc44(nfc1,nfc2)
         enddo
         call wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
         do nfc1=1,nfcmax2
            nth=nthnfc2(nfc1)
            nph=nhhnfc2(nfc1)
            fmc44(nfc1,nfc2)=fv1f(nth,nph)
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

            do i =1,3
               do j = 1,3
                  fmc41_d(i,j)
     &                    =0.25d0*fmc41(i,j,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc41(i,j,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc41(i,j,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc41(i,j,nfcdiff,nfcadd4)
                  fmc4a1_d(i,j)
     &                    =0.25d0*fmc4a1(i,j,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc4a1(i,j,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc4a1(i,j,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc4a1(i,j,nfcdiff,nfcadd4)
                  fmca41_d(i,j)
     &                    =0.25d0*fmca41(i,j,nfcdiff,nfcadd1)
     &                    +0.25d0*fmca41(i,j,nfcdiff,nfcadd2)
     &                    +0.25d0*fmca41(i,j,nfcdiff,nfcadd3)
     &                    +0.25d0*fmca41(i,j,nfcdiff,nfcadd4)
               enddo
               fmc42_d(i)
     &                    =0.25d0*fmc42(i,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc42(i,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc42(i,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc42(i,nfcdiff,nfcadd4)
               fmc43_d(i)
     &                    =0.25d0*fmc43(i,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc43(i,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc43(i,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc43(i,nfcdiff,nfcadd4)
               fmc4a2_d(i)
     &                    =0.25d0*fmc4a2(i,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc4a2(i,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc4a2(i,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc4a2(i,nfcdiff,nfcadd4)
               fmca43_d(i)
     &                    =0.25d0*fmca43(i,nfcdiff,nfcadd1)
     &                    +0.25d0*fmca43(i,nfcdiff,nfcadd2)
     &                    +0.25d0*fmca43(i,nfcdiff,nfcadd3)
     &                    +0.25d0*fmca43(i,nfcdiff,nfcadd4)
            enddo
            fmd(4,4,1,nfc1,nfc2)
     &                    =fmc41_d(1,1) 
     &                    +fmc41_d(2,1)*mm1
     &                    +fmc41_d(3,1)*nn1
     &                    +fmc41_d(1,2)     *mm2
     &                    +fmc41_d(2,2)*mm1 *mm2
     &                    +fmc41_d(3,2)*nn1 *mm2
     &                    +fmc41_d(1,3)     *nn2
     &                    +fmc41_d(2,3)*mm1 *nn2
     &                    +fmc41_d(3,3)*nn1 *nn2
            fmd(4,4,2,nfc1,nfc2)
     &                    =fmc42_d(1) 
     &                    +fmc42_d(2)*mm1
     &                    +fmc42_d(3)*nn1
            fmd(4,4,3,nfc1,nfc2)
     &                    =fmc43_d(1) 
     &                    +fmc43_d(2)*mm2
     &                    +fmc43_d(3)*nn2
            fmd(4,4,4,nfc1,nfc2)
     &                    =0.25d0*fmc44(nfcdiff,nfcadd1)
     &                    +0.25d0*fmc44(nfcdiff,nfcadd2)
     &                    +0.25d0*fmc44(nfcdiff,nfcadd3)
     &                    +0.25d0*fmc44(nfcdiff,nfcadd4)
            do i=1,3
               fmd(i,4,1,nfc1,nfc2)
     &                    =fmca41_d(i,1) 
     &                    +fmca41_d(i,2)     *mm2
     &                    +fmca41_d(i,3)     *nn2
               fmd(i,4,3,nfc1,nfc2)
     &                    =fmca43_d(i) 
               fmd(4,i,1,nfc1,nfc2)
     &                    =fmc4a1_d(i,1) 
     &                    +fmc4a1_d(i,2)*mm1
     &                    +fmc4a1_d(i,3)*nn1
               fmd(4,i,2,nfc1,nfc2)
     &                    =fmc4a2_d(i) 

            enddo
            do j=1,3
               do i=1,3
                  fmd(i,j,1,nfc1,nfc2)
     &                    =0.25d0*fmc(i,j,1,nfcdiff,nfcadd1)
     &                    +0.25d0*fmc(i,j,1,nfcdiff,nfcadd2)
     &                    +0.25d0*fmc(i,j,1,nfcdiff,nfcadd3)
     &                    +0.25d0*fmc(i,j,1,nfcdiff,nfcadd4)
                enddo
            enddo
         enddo
      enddo

      return

      end subroutine wmfem_calculate_plasma

      subroutine wmfem_calculate_plasma_sub(rho,ns,
     & fmc41,fmc42,fmc43,fmc44,fmca41,fmca43,fmc4a1,fmc4a2)

      USE wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      integer,intent(in):: ns
      complex(8),dimension(3,3,nfcmax2,nfcmax),intent(out)::
     &                                   fmc41,fmc4a1,fmca41
      complex(8),dimension(3,nfcmax2,nfcmax),intent(out)::
     &                                    fmc42,fmc43,fmc4a2,fmca43
      complex(8),dimension(nfcmax2,nfcmax),intent(out):: fmc44
!      complex(8),dimension(3,3,4,nfcmax2,nfcmax):: fmc
      complex(8),dimension(:,:,:,:,:),pointer:: fmc
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
!      real(8),dimension(3,3,nthmax2,nhhmax2) :: muminva 
      real(8),dimension(3,3) :: muminv 
      real(8),dimension(3,3) :: mum 
      complex(8),dimension(3,3) :: cinv
      real(8),dimension(nthmax2,nhhmax2):: gja

      real(8) ::gj

      complex(8):: cfactor
      complex(8):: cfactora4,cfactor4a,cfactor44,cfactoraa
      complex(8):: csum1,csum2,csum3,csum4
      complex(8):: csuma12,csuma13,csuma2,csuma3,csuma4
      integer:: id,jd,im,jm
      integer:: i,j,k,nfc1,nfc2,nth,nph
      integer:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      integer:: mmadd1,mmadd2,nnadd1,nnadd2
      integer:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactora4=ci*(2.d0*pi*crf*1.d6)/vc
!!!!!!!!!!      cfactora4=ci/vc**2
!!!!!!!!!!       cfactora4=0d0 
       cfactor4a=ci*(2.d0*pi*crf*1.d6)/vc
!!!!!!!!!!      cfactor4a=ci
!!!!!!!!!!       cfactor4a=0d0
!!!!!!!!!!      cfactor4a=(2.d0*pi*crf*1.d6)**2/vc**2
      cfactor44=-1d0
!!!!!!!!!      cfactor44=-1d0/(2.d0*pi*crf*1.d6)
!!!!!!!!!      cfactor44=ci*(2.d0*pi*crf*1.d6)/vc**2
!!!!!!!!!      cfactor44=00

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

!      call wmfem_inverse_tensor(muma,muminva)
      if (ns == 0)then
        fmc(:,:,:,:,:)=0d0
        fmc(1,1,:,:,:)=1d0
        fmc(2,2,:,:,:)=1d0
        fmc(3,3,:,:,:)=1d0
      else
        call wmfem_disp_tensor(rho,ns,fmc)
      endif
      
      do nfc2=1,nfcmax
          do nfc1=1,nfcmax2
              nth=nthnfc2(nfc1)
              nph=nhhnfc2(nfc1)
              gj=gja(nth,nph)
              mum(:,:)=muma(:,:,nth,nph)
              call wmfem_inverse_tensor(mum,muminv)
              do i =1,3
                  cinv(i,1)=muminv(i,1)
                 do j =2,3
                  cinv(i,j)=ci*muminv(i,j)
                 enddo
              enddo
              do im=2,3
                  do jm =2,3
                     csum1=0d0
                     do id =1,3 
                         do jd= 1,3
                         csum1=csum1
     &                   + fmc(id,jd,1,nfc1,nfc2)
     &                   * conjg(cinv(id,jm))*cinv(jd,im)
                         enddo
                     enddo
                     fmc41(im,jm,nfc1,nfc2)=cfactor44 * csum1*gj
                  enddo
              enddo
              do im=2,3
                   csum2=0d0
                   csum3=0d0
                   do id =1,3 
                       do jd= 1,3
                       csum2=csum2
     &                 + fmc(id,jd,2,nfc1,nfc2)
     &                 * conjg(cinv(id,1))*cinv(jd,im)

                       csum3=csum3
     &                 + fmc(id,jd,3,nfc1,nfc2)
     &                 * conjg(cinv(id,im))*cinv(jd,1)
                       enddo
                   enddo
                   fmc42(im,nfc1,nfc2)= cfactor44 * csum2*gj
                   fmc43(im,nfc1,nfc2)= cfactor44 * csum3*gj
              enddo

              csum4=0d0
              do id =1,3 
                  do jd= 1,3
                    csum4=csum4 
     &                   + fmc(id,jd,4,nfc1,nfc2)
     &                   * conjg(cinv(id,1))*cinv(jd,1)
                  enddo
              enddo
              fmc44(nfc1,nfc2)= cfactor44 * csum4*gj

              do i =1,3
                  csuma3=0d0
                  csuma12=0d0
                  csuma13=0d0
                  do id = 1,3
                     csuma3=csuma3 
     &                   + fmc(i,id,2,nfc1,nfc2)*cinv(id,1)
                     csuma12=csuma12 
     &                   + fmc(i,id,1,nfc1,nfc2)*cinv(id,2)
                     csuma13=csuma13 
     &                   + fmc(i,id,1,nfc1,nfc2)*cinv(id,3)
                  enddo
                  fmca43(i,nfc1,nfc2)= cfactora4 * csuma3*gj
                  fmca41(i,2,nfc1,nfc2)= cfactora4 * csuma12*gj
                  fmca41(i,3,nfc1,nfc2)= cfactora4 * csuma13*gj
              enddo

              do i =1,3
                  csuma3=0d0
                  csuma12=0d0
                  csuma13=0d0
                  do id = 1,3
                     csuma3=csuma3 
     &                   + fmc(id,i,3,nfc1,nfc2)*conjg(cinv(id,1))
                     csuma12=csuma12 
     &                   + fmc(id,i,1,nfc1,nfc2)*conjg(cinv(id,2))
                     csuma13=csuma13 
     &                   + fmc(id,i,1,nfc1,nfc2)*conjg(cinv(id,3))
                  enddo
                  fmc4a2(i,nfc1,nfc2)= cfactor4a * csuma3*gj
                  fmc4a1(i,2,nfc1,nfc2)= cfactor4a * csuma12*gj
                  fmc4a1(i,3,nfc1,nfc2)= cfactor4a * csuma13*gj
              enddo
          enddo    
      enddo

      deallocate(fmc)

      return

      end subroutine wmfem_calculate_plasma_sub

!----- calculate coefficint matrix fmd (E cylindrical +,-,para)-----

!     ****** CALCULATE METRIC AND CONVERSION TENSOR ******

!      SUBROUTINE wmfem_tensors(rho,gma,muma,dmuma,gja)
      SUBROUTINE wmfem_tensors(
     &      rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

      use wmfem_comm
      IMPLICIT NONE
      real(8),intent(in):: rho
      real(8),intent(out),dimension(3,3,nthmax2,nhhmax2):: 
     &     gma,muma,dmuma
      real(8),intent(out),dimension(3,3,nthmax2,nhhmax2):: 
     &     gpa,gmuma
      real(8),intent(out),dimension(3,nthmax2,nhhmax2):: 
     &     dgjgmuma
      real(8),intent(out),dimension(nthmax2,nhhmax2):: gja
      real(8),dimension(3,3):: gm,gmp,gmm,mum,mump,mumm
      real(8),dimension(3,3):: gp,gpp,gpm
      real(8),dimension(3,3):: gm_inv,gmp_inv,gmm_inv
      real(8):: th,ph,dth,dph,drhom,drhop,gj,gjp,gjm
      real(8):: babs,bsuprh,bsupth,bsupph,rhol
      integer:: nth,nph,i,j

      dth=2.d0*pi/nthmax2
      dph=2.d0*pi/nhhmax2
      IF(rho.EQ.0.d0) THEN
!         rhol=1.d-6
         rhol=1.d-6
         drhom=0.d0
      else
         rhol=rho
         drhom=1.d-6
      endif
      drhop=1.d-6
      bsuprh=0d0
!      print *,rho
      DO nph=1,nhhmax2
         ph=dph*(nph-1)
         DO nth=1,nthmax2
            th=dth*(nth-1)
            CALL wmfem_metrics(rhol,th,ph,gm,gj)
            CALL wmfem_magnetic(rhol,th,ph,babs,bsupth,bsupph)
            CALL wmfem_rotation_tensor(
     &           gm,gj,babs,bsuprh,bsupth,bsupph,mum)
            CALL wmfem_div_tensor(gm,gj,gp)

            CALL wmfem_metrics(rhol-drhom,th,ph,gmm,gjm)
            CALL wmfem_magnetic(rhol-drhom,th,ph,babs,bsupth,bsupph)
            CALL wmfem_rotation_tensor(
     &           gmm,gjm,babs,bsuprh,bsupth,bsupph,mumm)
            CALL wmfem_div_tensor(gmm,gjm,gpm)

            CALL wmfem_metrics(rhol+drhop,th,ph,gmp,gjp)
            CALL wmfem_magnetic(rhol+drhop,th,ph,babs,bsupth,bsupph)
            CALL wmfem_rotation_tensor(
     &           gmp,gjp,babs,bsuprh,bsupth,bsupph,mump)
            CALL wmfem_div_tensor(gmp,gjp,gpp)

            DO j=1,3
               DO i=1,3
                  gma(i,j,nth,nph)=gm(i,j)
                  gpa(i,j,nth,nph)=gp(i,j)
                  muma(i,j,nth,nph)=mum(i,j)
                  dmuma(i,j,nth,nph)
     &                 =(mump(i,j)-mumm(i,j))/(drhom+drhop)
                  gmuma(i,j,nth,nph)
     &                 =gp(i,1)*mum(1,j)
     &                + gp(i,2)*mum(2,j)+gp(i,3)*mum(3,j)
               END DO
               dgjgmuma(j,nth,nph)
     &              =(gjp*(gpp(1,1)*mump(1,j)
     &                   + gpp(1,2)*mump(2,j)+gpp(1,3)*mump(3,j))
     &            -   gjm*(gpm(1,1)*mumm(1,j)
     &                   + gpm(1,2)*mumm(2,j)+gpm(1,3)*mumm(3,j))
     &               )/(drhom+drhop)

            END DO
            gja(nth,nph)=gj
!            IF(rho.EQ.0.d0) THEN
!              gja(nth,nph)=0d0
!            ENDIF

!                 print *,nph,nhhmax2,nth
!                 print *,'gm'
!                 print *,gp
!                 print *,'gmp'
!                 print *,gpm
!                 print *,'gmm'
!                 print *,gpp
         ENDDO
!                 print *,ph,th
      ENDDO
      RETURN
      END SUBROUTINE wmfem_tensors

!        ****** CALCULATE ROTATION TENSOR ******

      SUBROUTINE wmfem_rotation_tensor(
     &         gm,gj,babs,bsuprh,bsupth,bsupph,mum)

      use wmfem_comm
      IMPLICIT NONE
      real(8),dimension(3,3),intent(in):: gm
      real(8),intent(in):: gj,babs,bsuprh,bsupth,bsupph
      real(8),dimension(3,3),intent(out):: mum
      real(8):: tc1,tc2,tc3,gf11,gf12,gf22,gf32,gf23

!        ----- Calculate rotation matrix mum -----

      tc1=bsuprh/babs
      tc2=bsupth/babs
      tc3=bsupph/babs

!        ----- gf11=gj*SQRT(gm11) -----

      gf11 = SQRT(gm(2,2)*gm(3,3)-gm(2,3)*gm(3,2))
      gf12 =(tc2*(gm(2,3)*gm(1,2)-gm(2,2)*gm(1,3))
     &     +tc3*(gm(3,3)*gm(1,2)-gm(2,3)*gm(1,3)))
      gf22 = gm(1,3)*gm(2,2)-gm(1,2)*gm(2,3)
      gf23 = gm(1,3)*gm(3,2)-gm(1,2)*gm(3,3)

      mum(1,1)= gj/gf11
      mum(2,1)= 0.D0
      mum(3,1)= 0.D0
      mum(1,2)= gf12/gf11
      mum(2,2)= tc1*gf22/gf11 + tc3*gf11
      mum(3,2)= tc1*gf23/gf11 - tc2*gf11
      mum(1,3)= tc1*gm(1,1)+tc2*gm(1,2)+tc3*gm(1,3)
      mum(2,3)= tc1*gm(2,1)+tc2*gm(2,2)+tc3*gm(2,3)
      mum(3,3)= tc1*gm(3,1)+tc2*gm(3,2)+tc3*gm(3,3)

      end subroutine wmfem_rotation_tensor

!        ****** CALCULATE DIVERGENCE TENSOR ******

      SUBROUTINE wmfem_div_tensor(
     &         gm,gj,gp)

      use wmfem_comm
      IMPLICIT NONE
      real(8),dimension(3,3),intent(in):: gm
      real(8),intent(in):: gj
      real(8),dimension(3,3),intent(out):: gp

      gp(1,1) = gm(2,2)*gm(3,3) -gm(3,2)*gm(2,3)
      gp(2,1) = gm(3,2)*gm(1,3) -gm(1,2)*gm(3,3)
      gp(3,1) = gm(1,2)*gm(2,3) -gm(2,2)*gm(1,3)
      gp(1,2) = gm(2,3)*gm(3,1) -gm(3,3)*gm(2,1)
      gp(2,2) = gm(3,3)*gm(1,1) -gm(1,3)*gm(3,1)
      gp(3,2) = gm(1,3)*gm(2,1) -gm(2,3)*gm(1,1)
      gp(1,3) = gm(2,1)*gm(3,2) -gm(3,1)*gm(2,2)
      gp(2,3) = gm(3,1)*gm(1,2) -gm(1,1)*gm(3,2)
      gp(3,3) = gm(1,1)*gm(2,2) -gm(2,1)*gm(1,2)

      gp(:,:) = gp (:,:)/(gj**2)

      end subroutine wmfem_div_tensor

      SUBROUTINE wmfem_inverse_tensor(mum,mum_inv)
      IMPLICIT NONE
      real(8),dimension(3,3),intent(in)  :: mum
      real(8),dimension(3,3),intent(out) :: mum_inv
      real(8) ::detA
       
      detA= mum(1,1)*mum(2,2)*mum(3,3)
     &    + mum(2,1)*mum(3,2)*mum(1,3)
     &    + mum(3,1)*mum(1,2)*mum(2,3)
     &    - mum(1,1)*mum(3,2)*mum(2,3)
     &    - mum(3,1)*mum(2,2)*mum(1,3)
     &    - mum(2,1)*mum(1,2)*mum(3,3)

      mum_inv(1,1) = mum(2,2)*mum(3,3)-mum(2,3)*mum(3,2)
      mum_inv(2,1) = mum(2,3)*mum(3,1)-mum(2,1)*mum(3,3)
      mum_inv(3,1) = mum(2,1)*mum(3,2)-mum(2,2)*mum(3,1)
      mum_inv(1,2) = mum(3,2)*mum(1,3)-mum(3,3)*mum(1,2)
      mum_inv(2,2) = mum(3,3)*mum(1,1)-mum(3,1)*mum(1,3)
      mum_inv(3,2) = mum(3,1)*mum(1,2)-mum(3,2)*mum(1,1)
      mum_inv(1,3) = mum(1,2)*mum(2,3)-mum(1,3)*mum(2,2)
      mum_inv(2,3) = mum(1,3)*mum(2,1)-mum(1,1)*mum(2,3)
      mum_inv(3,3) = mum(1,1)*mum(2,2)-mum(1,2)*mum(2,1)
      mum_inv(:,:) = mum_inv(:,:)/detA

      end subroutine wmfem_inverse_tensor

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
      dph=2.d0*pi/nhhmax2
      DO nfc2=1,nfcmax2
         nth=nthnfc2(nfc2)
         nph=nhhnfc2(nfc2)
         th=dth*(nth-1)
         ph=dph*(nph-1)
         DO nfc=1,nfcmax
            mm=mmnfc(nfc)
            nn=nnnfc(nfc)
            CALL wmfem_dielectric(rho,th,ph,mm,nn,ns,fml)
            DO j=1,3
               DO i=1,3
                  fmc(i,j,1,nfc2,nfc)=fml(i,j)
                  fmc(i,j,2,nfc2,nfc)=fml(i,j)
                  fmc(i,j,3,nfc2,nfc)=fml(i,j)
                  fmc(i,j,4,nfc2,nfc)=fml(i,j)
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

