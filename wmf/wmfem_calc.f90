! $Id$

MODULE wmcalc

PRIVATE
PUBLIC wmfem_calculate_local,wmfem_rotation_tensor,wmfem_inverse_tensor, &
       wmfem_tensors

CONTAINS

!     ----- calculate coefficint matrix fma -----

      SUBROUTINE wmfem_calculate_local(nr,ns,fml)

      USE wmfem_comm
      USE libfem_hhg
      IMPLICIT NONE
      INTEGER,INTENT(IN):: nr,ns
      COMPLEX(8),DIMENSION(mbmax,mbmax),INTENT(OUT):: fml
      COMPLEX(8),DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: fmd,fmd_p
      COMPLEX(8):: moment0,moment1,moment2,moment3
      REAL(8):: rho0,xi,yi
      REAL(8),PARAMETER:: m0=4.D0
      REAL(8),PARAMETER:: m2=5.D0/16.D0
      REAL(8),PARAMETER:: m4=41.D0/1024.D0
      REAL(8),PARAMETER:: m6=365.D0/65536.D0
      INTEGER ::inod,nfc2,nfc1
      INTEGER :: i,j,k
      REAL(8) ::drho

      ALLOCATE(fmd_p(4,4,4,nfcmax,nfcmax,4))
      ALLOCATE(fmd(4,4,4,nfcmax,nfcmax,4))

      drho=rhoa(nr+1)-rhoa(nr)
      DO inod=1,4
         SELECT CASE(inod)
         CASE(1)
            rho0=(7.D0*rhoa(nr)+     rhoa(nr+1))/8.D0
         CASE(2)
            rho0=(5.D0*rhoa(nr)+3.D0*rhoa(nr+1))/8.D0
         CASE(3)
            rho0=(3.D0*rhoa(nr)+5.D0*rhoa(nr+1))/8.D0
         CASE(4)
            rho0=(     rhoa(nr)+7.D0*rhoa(nr+1))/8.D0
         END SELECT

         IF(ns.EQ.0) THEN
            IF(mdlwmf.EQ.1) THEN
               CALL wmfem_calculate_vacuum(rho0,fmd_p(:,:,:,:,:,inod))
            ELSE IF(mdlwmf.EQ.2) THEN
!               CALL wmfem_calculate_vacuum_c(rho0,fmd_p(:,:,:,:,:,inod))
            END IF
         ELSE
            IF(mdlwmf.EQ.1) THEN
               CALL wmfem_calculate_plasma(rho0,ns,fmd_p(:,:,:,:,:,inod))
            ELSE IF(mdlwmf.EQ.2) THEN
!               CALL wmfem_calculate_plasma_c(rho0,ns,fmd_p(:,:,:,:,:,inod))
            END IF
         END IF
      END DO

      
! ------ calculate coefficients of basis for profile from four points 
      DO nfc2=1,nfcmax
      DO nfc1=1,nfcmax
         DO k=1,4
         DO j=1,4
            DO i=1,4
               moment0=0.D0
               moment1=0.D0
               moment2=0.D0
               moment3=0.D0
               DO inod=1,4
                  xi=0.125D0*(2*inod-1)
                  yi=xi-0.5D0
                  moment0=moment0+      fmd_p(i,j,k,nfc1,nfc2,inod)
                  moment1=moment1+yi   *fmd_p(i,j,k,nfc1,nfc2,inod)
                  moment2=moment2+yi**2*fmd_p(i,j,k,nfc1,nfc2,inod)
                  moment3=moment3+yi**3*fmd_p(i,j,k,nfc1,nfc2,inod)
               END DO
               fmd(i,j,k,nfc1,nfc2,1) &
                    =(m4*moment0-m2*moment2)/(m0*m4-m2**2)
               fmd(i,j,k,nfc1,nfc2,2) &
                    =(m6*moment1-m4*moment3)/(m2*m6-m4**2)
               fmd(i,j,k,nfc1,nfc2,3) &
                    =(m2*moment0-m0*moment2)/(m2**2-m0*m4)*2.D0
               fmd(i,j,k,nfc1,nfc2,4) &
                    =(m4*moment1-m2*moment3)/(m4**2-m2*m6)*6.D0
            END DO
         END DO
         END DO
      END DO
      END DO
      CALL fem_hhg(fmd,drho,fml)
      DEALLOCATE(fmd,fmd_p)
      RETURN
    END SUBROUTINE wmfem_calculate_local


!----- calculate coefficint matrix fma (A  rho,theta,phi )-----

      SUBROUTINE wmfem_calculate_vacuum(rho,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      REAL(8),INTENT(IN):: rho
      COMPLEX(8),DIMENSION(4,4,4,nfcmax,nfcmax),INTENT(OUT):: fmd
      COMPLEX(8),DIMENSION(4,4,4,nfcmax,nfcmax):: fmd_plasma
      COMPLEX(8),DIMENSION(nthmax2,nhhmax2):: fv1,fv1f
      COMPLEX(8),DIMENSION(3,3,3,3,nfcmax2):: fmv1
      COMPLEX(8),DIMENSION(3,3,3,nfcmax2):: fmv2,fmv3
      COMPLEX(8),DIMENSION(3,3,nfcmax2):: fmv4
      INTEGER:: i,j,k,nfc,nfc1,nfc2
      INTEGER:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
      INTEGER:: imn,imn1,imn2,nfcdiff

      CALL wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

      DO i=1,3
         DO j=1,3
            DO imn1=1,3
               DO imn2=1,3
                  DO nfc2=1,nfcmax2
                     nth=nthnfc2(nfc2)
                     nph=nhhnfc2(nfc2)
                     fv1(nth,nph)=fmv1(i,j,imn1,imn2,nfc2)
                  END DO
                  CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
                  DO nfc2=1,nfcmax2
                     nth=nthnfc2(nfc2)
                     nph=nhhnfc2(nfc2)
                     fmv1(i,j,imn1,imn2,nfc2)=fv1f(nth,nph)
                  END DO
               END DO
            END DO    
            DO imn1=1,3
               DO nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fv1(nth,nph)=fmv2(i,j,imn1,nfc2)
               END DO
               CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               DO nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fmv2(i,j,imn1,nfc2)=fv1f(nth,nph)
               END DO
            END DO
            DO imn1=1,3
               DO nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fv1(nth,nph)=fmv3(i,j,imn1,nfc2)
               END DO
               CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               DO nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nhhnfc2(nfc2)
                  fmv3(i,j,imn1,nfc2)=fv1f(nth,nph)
               END DO
            END DO
            DO nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nhhnfc2(nfc2)
               fv1(nth,nph)=fmv4(i,j,nfc2)
            END DO
            CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            DO nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nhhnfc2(nfc2)
               fmv4(i,j,nfc2)=fv1f(nth,nph)
            END DO
         END DO
      END DO

!      CALL wmfem_calculate_plasma(rho,0,fmd)
      CALL wmfem_calculate_plasma(rho,0,fmd_plasma)
      fmd(:,4,:,:,:)=fmd_plasma(:,4,:,:,:)
      fmd(4,1,:,:,:)=fmd_plasma(4,1,:,:,:)
      fmd(4,2,:,:,:)=fmd_plasma(4,2,:,:,:)
      fmd(4,3,:,:,:)=fmd_plasma(4,3,:,:,:)

      DO nfc1=1,nfcmax          ! Fit to fmd and adjust m and n
         nn1=nnnfc(nfc1)
         mm1=mmnfc(nfc1)
         DO nfc2=1,nfcmax
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)

            nndiff=nn1-nn2
            IF(nndiff.LT.0) nndiff=nndiff+nhhmax2
            mmdiff=mm1-mm2
            IF(mmdiff.lt.0) mmdiff=mmdiff+nthmax2
            nfcdiff=nthmax2*nndiff+mmdiff+1
!            if(nfcdiff.le.0) nfcdiff=nfcdiff+nfcmax2

            DO j=1,3
               DO i=1,3
                  fmd(i,j,1,nfc1,nfc2) &
                      =fmv1(i,j,1,1,nfcdiff) &
                      +fmv1(i,j,2,1,nfcdiff)*mm1 &
                      +fmv1(i,j,3,1,nfcdiff)*nn1 &
                      +fmv1(i,j,1,2,nfcdiff)    *mm2 &
                      +fmv1(i,j,2,2,nfcdiff)*mm1*mm2 &
                      +fmv1(i,j,3,2,nfcdiff)*nn1*mm2 &
                      +fmv1(i,j,1,3,nfcdiff)    *nn2 &
                      +fmv1(i,j,2,3,nfcdiff)*mm1*nn2 &
                      +fmv1(i,j,3,3,nfcdiff)*nn1*nn2
                  fmd(i,j,2,nfc1,nfc2) &
                      =fmv2(i,j,1,nfcdiff) &
                      +fmv2(i,j,2,nfcdiff)*mm1 &
                      +fmv2(i,j,3,nfcdiff)*nn1
                  fmd(i,j,3,nfc1,nfc2) &
                      =fmv3(i,j,1,nfcdiff) &
                      +fmv3(i,j,2,nfcdiff)    *mm2 &
                      +fmv3(i,j,3,nfcdiff)    *nn2
                 fmd(i,j,4,nfc1,nfc2) &
                      =fmv4(i,j,nfcdiff)

            END DO
            END DO
         END DO
      END DO
      RETURN
      END SUBROUTINE wmfem_calculate_vacuum

!----- calculate vacuum matrix fms -----

      SUBROUTINE wmfem_calculate_vacuum_sub(rho,fmv1,fmv2,fmv3,fmv4)

      USE wmfem_comm
      IMPLICIT NONE
      REAL(8),INTENT(IN):: rho
      COMPLEX(8),DIMENSION(3,3,3,3,nfcmax2),INTENT(OUT):: fmv1
      COMPLEX(8),DIMENSION(3,3,3,nfcmax2),INTENT(OUT):: fmv2,fmv3
      COMPLEX(8),DIMENSION(3,3,nfcmax2),INTENT(OUT):: fmv4

      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma 
      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      REAL(8),DIMENSION(3,nthmax2,nhhmax2) :: dgjgmuma
      REAL(8),DIMENSION(nthmax2,nhhmax2):: gja
!      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: muminva 
      REAL(8),DIMENSION(3,3) :: muminva 

      COMPLEX(8),DIMENSION(3,3,3):: cq
      COMPLEX(8),DIMENSION(3,3):: cp
      COMPLEX(8),DIMENSION(3,3):: cqd
      COMPLEX(8),DIMENSION(3):: cpd
      INTEGER:: idiv
      INTEGER:: i,j,k,l,nthm,nthp,nphm,nphp
      INTEGER:: imn1,imn2
      INTEGER:: nph,nth
      REAL(8):: dph,dth
      COMPLEX(8):: csum1,csum2,csum3,csum4
      COMPLEX(8):: csumd1,csumd2,csumd3,csumd4
      INTEGER:: nrl,nfc2
      REAL(8):: drhob,muma3b,muma2b,drhoa,muma3a,muma2a,muma30,muma20
      REAL(8):: muma3d,muma2d,gj

      COMPLEX(8):: cfactor
      COMPLEX(8):: cfactor_div

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

      cfactor_div=1d0 !(2.d0*pi*crf*1.d6)
!      cfactor_div=0d0 !(2.d0*pi*crf*1.d6)**2

!!!!

      CALL wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)
!      CALL wmfem_inverse_tensor(muma,muminva)



      DO nfc2=1,nfcmax2
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

         DO j=1,3
            cq(1,j,1)=((muma(3,j,nthp,nph) &
                      -muma(3,j,nthm,nph))/dth &
                     -(muma(2,j,nth,nphp) &
                      -muma(2,j,nth,nphm))/dph ) /gj
            cq(1,j,2)=+ci*muma(3,j,nth,nph) /gj
            cq(1,j,3)=-ci*muma(2,j,nth,nph) /gj

            cq(2,j,1)=((muma(1,j,nth,nphp) &
                      -muma(1,j,nth,nphm))/dph &
                      -dmuma(3,j,nth,nph)) /gj
            cq(2,j,2)=0.d0
            cq(2,j,3)=+ci*muma(1,j,nth,nph) /gj

            cq(3,j,1)=(dmuma(2,j,nth,nph) &
                    -(muma(1,j,nthp,nph) &
                     -muma(1,j,nthm,nph))/dth ) /gj
            cq(3,j,2)=-ci*muma(1,j,nth,nph) /gj
            cq(3,j,3)=0.d0

            cp(1,j)=0.d0
            cp(2,j)=-muma(3,j,nth,nph) /gj
            cp(3,j)= muma(2,j,nth,nph) /gj

            
            cqd(j,1)= (dgjgmuma(j,nth,nph) &
                     +( gja(nthp,nph)*gmuma(2,j,nthp,nph) &
                       - gja(nthm,nph)*gmuma(2,j,nthm,nph))/dth &
                     +( gja(nth,nphp)*gmuma(3,j,nth,nphp) &
                       - gja(nth,nphm)*gmuma(3,j,nth,nphm))/dph &
                     )  /gj
            cqd(j,2)= +ci*gmuma(2,j,nth,nph)
            cqd(j,3)= +ci*gmuma(3,j,nth,nph)
                       
            cpd(j)=  gmuma(1,j,nth,nph)
         END DO

         DO imn2=1,3
            DO imn1=1,3
               DO j=1,3
                  DO i=1,3
                     csum1=0.d0
                     DO k=1,3
                        DO l=1,3
                           csum1=csum1-conjg(cq(k,i,imn1)) &
                                     *gma(k,l,nth,nph) &
                                     *cq(l,j,imn2)*gj
                        END DO
                     END DO
   
                     csumd1=-conjg(cqd(i,imn1))*cqd(j,imn2) *gj

 
                     if(i.eq.j.and.imn1.eq.1.and.imn2.eq.1) then
                        fmv1(i,j,imn1,imn2,nfc2) &
                            =csum1+csumd1*cfactor_div +cfactor*gj
!     &                       =csum1+csumd1 +cfactor*gj
!     &                       =csum1-csumd1 +cfactor*gj
!     &                       =csum1 +cfactor*gj
                     else
                        fmv1(i,j,imn1,imn2,nfc2)=csum1 &
                         + csumd1*cfactor_div
!                        fmv1(i,j,imn1,imn2,nfc2)=csum1 +csumd1
!                        fmv1(i,j,imn1,imn2,nfc2)=csum1 -csumd1
!                        fmv1(i,j,imn1,imn2,nfc2)=csum1
                     endif
                  END DO
               END DO
            END DO
         END DO
            
         DO imn1=1,3
            DO j=1,3
               DO i=1,3
                  csum2=0.d0
                  csum3=0.d0
                  DO k=1,3
                     DO l=1,3
                        csum2=csum2-conjg(cp(k,i)) &
                                  *gma(k,l,nth,nph) &
                                  *cq(l,j,imn1)  *gj
                        csum3=csum3-conjg(cq(k,i,imn1)) &
                                  *gma(k,l,nth,nph) &
                                  *cp(l,j)  *gj
                     END DO
                  END DO

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
               END DO
           END DO
         END DO
            
         DO j=1,3
            DO i=1,3
               csum4=0.d0
               DO k=1,3
                  DO l=1,3
                     csum4=csum4-conjg(cp(k,i)) &
                               *gma(k,l,nth,nph) &
                               *cp(l,j)  *gj
                  END DO
               END DO
               csumd4=-conjg(cpd(i))*cpd(j)  *gj

               fmv4(i,j,nfc2)=csum4 + csumd4* cfactor_div
!               fmv4(i,j,nfc2)=csum4 - csumd4
!               fmv4(i,j,nfc2)=csum4
            END DO
         END DO

      END DO
      RETURN
      END SUBROUTINE wmfem_calculate_vacuum_sub

!----- calculate coefficint matrix fma -----

      SUBROUTINE wmfem_calculate_plasma(rho,ns,fmd)

      USE wmfem_comm
      IMPLICIT NONE
      REAL(8),INTENT(IN):: rho
      INTEGER,INTENT(IN):: ns
      COMPLEX(8),DIMENSION(4,4,4,nfcmax,nfcmax),INTENT(OUT):: fmd
!      COMPLEX(8),DIMENSION(3,3,4,nfcmax2,nfcmax):: fmc
      COMPLEX(8),DIMENSION(:,:,:,:,:),ALLOCATABLE:: fmc
      COMPLEX(8),DIMENSION(nthmax2,nhhmax2):: fv1,fv1f

      COMPLEX(8),DIMENSION(3,3,nfcmax2,nfcmax):: &
                                        fmc41,fmc4a1,fmca41
      COMPLEX(8),DIMENSION(3,nfcmax2,nfcmax):: &
                                         fmc42,fmc43,fmc4a2,fmca43
      COMPLEX(8),DIMENSION(nfcmax2,nfcmax):: fmc44
      
      COMPLEX(8),DIMENSION(3,3) ::fmc41_d,fmc4a1_d,fmca41_d
      COMPLEX(8),DIMENSION(3) ::fmc42_d,fmc43_d,fmc4a2_d,fmca43_d
      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      REAL(8),DIMENSION(nthmax2,nhhmax2):: gja
      REAL(8),DIMENSION(3,nthmax2,nhhmax2) :: dgjgmuma
!      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: muminva 
!      REAL(8),DIMENSION(3,3) :: muminva 

      COMPLEX(8):: cfactor
      INTEGER:: i,j,k,nfc1,nfc2,nth,nph
      INTEGER:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      INTEGER:: mmadd1,mmadd2,nnadd1,nnadd2
      INTEGER:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactor=(2.d0*pi*crf*1.d6)**2/vc**2

      CALL wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

      if (ns == 0)then
        fmc(:,:,:,:,:)=0d0
        fmc(1,1,:,:,:)=1d0
        fmc(2,2,:,:,:)=1d0
        fmc(3,3,:,:,:)=1d0
      else
        CALL wmfem_disp_tensor(rho,ns,fmc)
      endif

      CALL wmfem_calculate_plasma_sub(rho,ns, &
              fmc41,fmc42,fmc43,fmc44,fmca41,fmca43,fmc4a1,fmc4a2)
      DO i =1,3 
         DO j = 1,3
            DO nfc2=1,nfcmax
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=cfactor*fmc(i,j,1,nfc1,nfc2)*gja(nth,nph)
               END DO
               CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc(i,j,1,nfc1,nfc2)=fv1f(nth,nph)
               END DO
            END DO
            DO nfc2=1,nfcmax
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=fmc41(i,j,nfc1,nfc2)
               END DO
               CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc41(i,j,nfc1,nfc2)=fv1f(nth,nph)
               END DO
            END DO
            DO nfc2=1,nfcmax
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=fmc4a1(i,j,nfc1,nfc2)
               END DO
               CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmc4a1(i,j,nfc1,nfc2)=fv1f(nth,nph)
               END DO
            END DO
            DO nfc2=1,nfcmax
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fv1(nth,nph)=fmca41(i,j,nfc1,nfc2)
               END DO
               CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
               DO nfc1=1,nfcmax2
                  nth=nthnfc2(nfc1)
                  nph=nhhnfc2(nfc1)
                  fmca41(i,j,nfc1,nfc2)=fv1f(nth,nph)
               END DO
            END DO
         END DO
         DO nfc2=1,nfcmax
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmc42(i,nfc1,nfc2)
            END DO
            CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmc42(i,nfc1,nfc2)=fv1f(nth,nph)
            END DO
         END DO
         DO nfc2=1,nfcmax
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmc43(i,nfc1,nfc2)
            END DO
            CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmc43(i,nfc1,nfc2)=fv1f(nth,nph)
            END DO
         END DO
         DO nfc2=1,nfcmax
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmc4a2(i,nfc1,nfc2)
            END DO
            CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmc4a2(i,nfc1,nfc2)=fv1f(nth,nph)
            END DO
         END DO
         DO nfc2=1,nfcmax
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fv1(nth,nph)=fmca43(i,nfc1,nfc2)
            END DO
            CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
            DO nfc1=1,nfcmax2
               nth=nthnfc2(nfc1)
               nph=nhhnfc2(nfc1)
               fmca43(i,nfc1,nfc2)=fv1f(nth,nph)
            END DO
         END DO
      END DO    
      DO nfc2=1,nfcmax
         DO nfc1=1,nfcmax2
            nth=nthnfc2(nfc1)
            nph=nhhnfc2(nfc1)
            fv1(nth,nph)=fmc44(nfc1,nfc2)
         END DO
         CALL wmsubfx(fv1,fv1f,nthmax2,nhhmax2)
         DO nfc1=1,nfcmax2
            nth=nthnfc2(nfc1)
            nph=nhhnfc2(nfc1)
            fmc44(nfc1,nfc2)=fv1f(nth,nph)
         END DO
      END DO

!     ----- Fit to fmd and adjust m and n -----

      DO nfc2=1,nfcmax
         mm2=mmnfc(nfc2)
         nn2=nnnfc(nfc2)
         DO nfc1=1,nfcmax
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

            DO i =1,3
               DO j = 1,3
                  fmc41_d(i,j) &
                         =0.25d0*fmc41(i,j,nfcdiff,nfcadd1) &
                         +0.25d0*fmc41(i,j,nfcdiff,nfcadd2) &
                         +0.25d0*fmc41(i,j,nfcdiff,nfcadd3) &
                         +0.25d0*fmc41(i,j,nfcdiff,nfcadd4)
                  fmc4a1_d(i,j) &
                         =0.25d0*fmc4a1(i,j,nfcdiff,nfcadd1) &
                         +0.25d0*fmc4a1(i,j,nfcdiff,nfcadd2) &
                         +0.25d0*fmc4a1(i,j,nfcdiff,nfcadd3) &
                         +0.25d0*fmc4a1(i,j,nfcdiff,nfcadd4)
                  fmca41_d(i,j) &
                         =0.25d0*fmca41(i,j,nfcdiff,nfcadd1) &
                         +0.25d0*fmca41(i,j,nfcdiff,nfcadd2) &
                         +0.25d0*fmca41(i,j,nfcdiff,nfcadd3) &
                         +0.25d0*fmca41(i,j,nfcdiff,nfcadd4)
               END DO
               fmc42_d(i) &
                         =0.25d0*fmc42(i,nfcdiff,nfcadd1) &
                         +0.25d0*fmc42(i,nfcdiff,nfcadd2) &
                         +0.25d0*fmc42(i,nfcdiff,nfcadd3) &
                         +0.25d0*fmc42(i,nfcdiff,nfcadd4)
               fmc43_d(i) &
                         =0.25d0*fmc43(i,nfcdiff,nfcadd1) &
                         +0.25d0*fmc43(i,nfcdiff,nfcadd2) &
                         +0.25d0*fmc43(i,nfcdiff,nfcadd3) &
                         +0.25d0*fmc43(i,nfcdiff,nfcadd4)
               fmc4a2_d(i) &
                         =0.25d0*fmc4a2(i,nfcdiff,nfcadd1) &
                         +0.25d0*fmc4a2(i,nfcdiff,nfcadd2) &
                         +0.25d0*fmc4a2(i,nfcdiff,nfcadd3) &
                         +0.25d0*fmc4a2(i,nfcdiff,nfcadd4)
               fmca43_d(i) &
                         =0.25d0*fmca43(i,nfcdiff,nfcadd1) &
                         +0.25d0*fmca43(i,nfcdiff,nfcadd2) &
                         +0.25d0*fmca43(i,nfcdiff,nfcadd3) &
                         +0.25d0*fmca43(i,nfcdiff,nfcadd4)
            END DO
            fmd(4,4,1,nfc1,nfc2) &
                         =fmc41_d(1,1) & 
                         +fmc41_d(2,1)*mm1 &
                         +fmc41_d(3,1)*nn1 &
                         +fmc41_d(1,2)     *mm2 &
                         +fmc41_d(2,2)*mm1 *mm2 &
                         +fmc41_d(3,2)*nn1 *mm2 &
                         +fmc41_d(1,3)     *nn2 &
                         +fmc41_d(2,3)*mm1 *nn2 &
                         +fmc41_d(3,3)*nn1 *nn2
            fmd(4,4,2,nfc1,nfc2) &
                         =fmc42_d(1)  &
                         +fmc42_d(2)*mm1 &
                         +fmc42_d(3)*nn1
            fmd(4,4,3,nfc1,nfc2) &
                         =fmc43_d(1) &
                         +fmc43_d(2)*mm2 &
                         +fmc43_d(3)*nn2
            fmd(4,4,4,nfc1,nfc2) &
                         =0.25d0*fmc44(nfcdiff,nfcadd1) &
                         +0.25d0*fmc44(nfcdiff,nfcadd2) &
                         +0.25d0*fmc44(nfcdiff,nfcadd3) &
                         +0.25d0*fmc44(nfcdiff,nfcadd4)
            DO i=1,3
               fmd(i,4,1,nfc1,nfc2) &
                         =fmca41_d(i,1) &
                         +fmca41_d(i,2)*mm1 &
                         +fmca41_d(i,3)*nn1
               fmd(i,4,3,nfc1,nfc2) &
                         =fmca43_d(i) 
               fmd(4,i,1,nfc1,nfc2) &
                         =fmc4a1_d(i,1) &
                         +fmc4a1_d(i,2)     *mm2 &
                         +fmc4a1_d(i,3)     *nn2
               fmd(4,i,2,nfc1,nfc2) &
                         =fmc4a2_d(i) 

            END DO
            DO j=1,3
               DO i=1,3
                  fmd(i,j,1,nfc1,nfc2) &
                         =0.25d0*fmc(i,j,1,nfcdiff,nfcadd1) &
                         +0.25d0*fmc(i,j,1,nfcdiff,nfcadd2) &
                         +0.25d0*fmc(i,j,1,nfcdiff,nfcadd3) &
                         +0.25d0*fmc(i,j,1,nfcdiff,nfcadd4)
                END DO
            END DO
         END DO
      END DO

      RETURN

      END SUBROUTINE wmfem_calculate_plasma

      SUBROUTINE wmfem_calculate_plasma_sub(rho,ns, &
      fmc41,fmc42,fmc43,fmc44,fmca41,fmca43,fmc4a1,fmc4a2)

      USE wmfem_comm
      IMPLICIT NONE
      REAL(8),INTENT(IN):: rho
      INTEGER,INTENT(IN):: ns
      COMPLEX(8),DIMENSION(3,3,nfcmax2,nfcmax),INTENT(OUT):: &
                                        fmc41,fmc4a1,fmca41
      COMPLEX(8),DIMENSION(3,nfcmax2,nfcmax),INTENT(OUT):: &
                                         fmc42,fmc43,fmc4a2,fmca43
      COMPLEX(8),DIMENSION(nfcmax2,nfcmax),INTENT(OUT):: fmc44
!      COMPLEX(8),DIMENSION(3,3,4,nfcmax2,nfcmax):: fmc
      COMPLEX(8),DIMENSION(:,:,:,:,:),ALLOCATABLE:: fmc
      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      REAL(8),DIMENSION(3,nthmax2,nhhmax2) :: dgjgmuma
!      REAL(8),DIMENSION(3,3,nthmax2,nhhmax2) :: muminva 
      REAL(8),DIMENSION(3,3) :: muminv 
      REAL(8),DIMENSION(3,3) :: mum 
      COMPLEX(8),DIMENSION(3,3) :: cinv
      REAL(8),DIMENSION(nthmax2,nhhmax2):: gja

      REAL(8) ::gj

      COMPLEX(8):: cfactor
      COMPLEX(8):: cfactora4,cfactor4a,cfactor44,cfactoraa
      COMPLEX(8):: csum1,csum2,csum3,csum4
      COMPLEX(8):: csuma12,csuma13,csuma2,csuma3,csuma4
      INTEGER:: id,jd,im,jm
      INTEGER:: i,j,k,nfc1,nfc2,nth,nph
      INTEGER:: nn1,nn2,nndiff,mm1,mm2,mmdiff,nfcdiff,nfcadd
      INTEGER:: mmadd1,mmadd2,nnadd1,nnadd2
      INTEGER:: nfcadd1,nfcadd2,nfcadd3,nfcadd4

      allocate(fmc(3,3,4,nfcmax2,nfcmax))

      cfactora4=ci*(2.d0*pi*crf*1.d6)/vc**2
!      cfactora4=ci/vc**2
!      cfactora4=0d0 
      cfactor4a=ci*(2.d0*pi*crf*1.d6)
!      cfactor4a=ci
!      cfactor4a=0d0
!      cfactor4a=(2.d0*pi*crf*1.d6)**2/vc**2
      cfactor44=-1d0
!      cfactor44=-1d0/(2.d0*pi*crf*1.d6)
!      cfactor44=ci*(2.d0*pi*crf*1.d6)/vc**2
!      cfactor44=00

      CALL wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

!      CALL wmfem_inverse_tensor(muma,muminva)
      if (ns == 0)then
        fmc(:,:,:,:,:)=0d0
        fmc(1,1,:,:,:)=1d0
        fmc(2,2,:,:,:)=1d0
        fmc(3,3,:,:,:)=1d0
      else
        CALL wmfem_disp_tensor(rho,ns,fmc)
      endif
      
      DO nfc2=1,nfcmax
          DO nfc1=1,nfcmax2
              nth=nthnfc2(nfc1)
              nph=nhhnfc2(nfc1)
              gj=gja(nth,nph)
              mum(:,:)=muma(:,:,nth,nph)
              CALL wmfem_inverse_tensor(mum,muminv)
              DO i =1,3
                  cinv(i,1)=muminv(i,1)
                 DO j =2,3
                  cinv(i,j)=ci*muminv(i,j)
                 END DO
              END DO
              DO im=2,3
                  DO jm =2,3
                     csum1=0d0
                     DO id =1,3 
                         DO jd= 1,3
                         csum1=csum1 &
                        + fmc(id,jd,1,nfc1,nfc2) &
                        * conjg(cinv(id,jm))*cinv(jd,im)
                         END DO
                     END DO
                     fmc41(im,jm,nfc1,nfc2)=cfactor44 * csum1*gj
                  END DO
              END DO
              DO im=2,3
                   csum2=0d0
                   csum3=0d0
                   DO id =1,3 
                       DO jd= 1,3
                       csum2=csum2 &
                      + fmc(id,jd,2,nfc1,nfc2) &
                      * conjg(cinv(id,1))*cinv(jd,im)

                       csum3=csum3 &
                      + fmc(id,jd,3,nfc1,nfc2) &
                      * conjg(cinv(id,im))*cinv(jd,1)
                       END DO
                   END DO
                   fmc42(im,nfc1,nfc2)= cfactor44 * csum2*gj
                   fmc43(im,nfc1,nfc2)= cfactor44 * csum3*gj
              END DO

              csum4=0d0
              DO id =1,3 
                  DO jd= 1,3
                    csum4=csum4 &
                        + fmc(id,jd,4,nfc1,nfc2) &
                        * conjg(cinv(id,1))*cinv(jd,1)
                  END DO
              END DO
              fmc44(nfc1,nfc2)= cfactor44 * csum4*gj

              DO i =1,3
                  csuma3=0d0
                  csuma12=0d0
                  csuma13=0d0
                  DO id = 1,3
                     csuma3=csuma3  +fmc(i,id,2,nfc1,nfc2)*cinv(id,1)
                     csuma12=csuma12+fmc(i,id,1,nfc1,nfc2)*cinv(id,2)
                     csuma13=csuma13+fmc(i,id,1,nfc1,nfc2)*cinv(id,3)
                  END DO
                  fmca43(i,nfc1,nfc2)= cfactora4 * csuma3*gj
                  fmca41(i,2,nfc1,nfc2)= cfactora4 * csuma12*gj
                  fmca41(i,3,nfc1,nfc2)= cfactora4 * csuma13*gj
              END DO

              DO i =1,3
                  csuma3=0d0
                  csuma12=0d0
                  csuma13=0d0
                  DO id = 1,3
                     csuma3=csuma3 &
                        + fmc(id,i,3,nfc1,nfc2)*conjg(cinv(id,1))
                    csuma12=csuma12 & 
                        + fmc(id,i,1,nfc1,nfc2)*conjg(cinv(id,2))
                     csuma13=csuma13 &
                        + fmc(id,i,1,nfc1,nfc2)*conjg(cinv(id,3))
                  END DO
                  fmc4a2(i,nfc1,nfc2)= cfactor4a * csuma3*gj
                  fmc4a1(i,2,nfc1,nfc2)= cfactor4a * csuma12*gj
                  fmc4a1(i,3,nfc1,nfc2)= cfactor4a * csuma13*gj
              END DO
          END DO    
      END DO

      deallocate(fmc)

      RETURN

      END SUBROUTINE wmfem_calculate_plasma_sub

!     ****** CALCULATE METRIC AND CONVERSION TENSOR ******

!      SUBROUTINE wmfem_tensors(rho,gma,muma,dmuma,gja)
      SUBROUTINE wmfem_tensors( &
           rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)

      use wmfem_comm
      IMPLICIT NONE
      REAL(8),INTENT(IN):: rho
      REAL(8),INTENT(OUT),DIMENSION(3,3,nthmax2,nhhmax2):: &
          gma,muma,dmuma
      REAL(8),INTENT(OUT),DIMENSION(3,3,nthmax2,nhhmax2):: &
          gpa,gmuma
      REAL(8),INTENT(OUT),DIMENSION(3,nthmax2,nhhmax2):: &
          dgjgmuma
      REAL(8),INTENT(OUT),DIMENSION(nthmax2,nhhmax2):: gja
      REAL(8),DIMENSION(3,3):: gm,gmp,gmm,mum,mump,mumm
      REAL(8),DIMENSION(3,3):: gp,gpp,gpm
      REAL(8),DIMENSION(3,3):: gm_inv,gmp_inv,gmm_inv
      REAL(8):: th,ph,dth,dph,drhom,drhop,gj,gjp,gjm
      REAL(8):: babs,bsuprh,bsupth,bsupph,rhol
      INTEGER:: nth,nph,i,j

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
            CALL wmfem_magnetic(rhol,th,ph,babs,bsuprh,bsupth,bsupph)
            CALL wmfem_rotation_tensor( &
                gm,gj,babs,bsuprh,bsupth,bsupph,mum)
            CALL wmfem_div_tensor(gm,gj,gp)

            CALL wmfem_metrics(rhol-drhom,th,ph,gmm,gjm)
            CALL wmfem_magnetic(rhol-drhom,th,ph, &
                               babs,bsuprh,bsupth,bsupph)
            CALL wmfem_rotation_tensor( &
                gmm,gjm,babs,bsuprh,bsupth,bsupph,mumm)
            CALL wmfem_div_tensor(gmm,gjm,gpm)

            CALL wmfem_metrics(rhol+drhop,th,ph,gmp,gjp)
            CALL wmfem_magnetic(rhol+drhop,th,ph, &
                               babs,bsuprh,bsupth,bsupph)
            CALL wmfem_rotation_tensor( &
                gmp,gjp,babs,bsuprh,bsupth,bsupph,mump)
            CALL wmfem_div_tensor(gmp,gjp,gpp)

            DO j=1,3
               DO i=1,3
                  gma(i,j,nth,nph)=gm(i,j)
                  gpa(i,j,nth,nph)=gp(i,j)
                  muma(i,j,nth,nph)=mum(i,j)
                  dmuma(i,j,nth,nph) &
                      =(mump(i,j)-mumm(i,j))/(drhom+drhop)
                  gmuma(i,j,nth,nph) &
                      =gp(i,1)*mum(1,j) &
                     + gp(i,2)*mum(2,j)+gp(i,3)*mum(3,j)
               END DO
               dgjgmuma(j,nth,nph) &
                   =(gjp*(gpp(1,1)*mump(1,j) &
                        + gpp(1,2)*mump(2,j)+gpp(1,3)*mump(3,j)) &
                 -   gjm*(gpm(1,1)*mumm(1,j) &
                        + gpm(1,2)*mumm(2,j)+gpm(1,3)*mumm(3,j)) &
                    )/(drhom+drhop)

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
         END DO
!                 print *,ph,th
      END DO
      RETURN
      END SUBROUTINE wmfem_tensors

!        ****** CALCULATE ROTATION TENSOR ******

      SUBROUTINE wmfem_rotation_tensor( &
              gm,gj,babs,bsuprh,bsupth,bsupph,mum)

      use wmfem_comm
      IMPLICIT NONE
      REAL(8),DIMENSION(3,3),INTENT(IN):: gm
      REAL(8),INTENT(IN):: gj,babs,bsuprh,bsupth,bsupph
      REAL(8),DIMENSION(3,3),INTENT(OUT):: mum
      REAL(8):: tc1,tc2,tc3,gf11,gf12,gf22,gf32,gf23

!        ----- Calculate rotation matrix mum -----

      tc1=bsuprh/babs
      tc2=bsupth/babs
      tc3=bsupph/babs

!        ----- gf11=gj*SQRT(gm11) -----

      gf11 = SQRT(gm(2,2)*gm(3,3)-gm(2,3)*gm(3,2))
      gf12 =(tc2*(gm(2,3)*gm(1,2)-gm(2,2)*gm(1,3)) &
            +tc3*(gm(3,3)*gm(1,2)-gm(2,3)*gm(1,3)))
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

      END SUBROUTINE wmfem_rotation_tensor

!        ****** CALCULATE DIVERGENCE TENSOR ******

      SUBROUTINE wmfem_div_tensor(gm,gj,gp)

      use wmfem_comm
      IMPLICIT NONE
      REAL(8),DIMENSION(3,3),INTENT(IN):: gm
      REAL(8),INTENT(IN):: gj
      REAL(8),DIMENSION(3,3),INTENT(OUT):: gp

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

      END SUBROUTINE wmfem_div_tensor

      SUBROUTINE wmfem_inverse_tensor(mum,mum_inv)
      IMPLICIT NONE
      REAL(8),DIMENSION(3,3),INTENT(IN)  :: mum
      REAL(8),DIMENSION(3,3),INTENT(OUT) :: mum_inv
      REAL(8) ::detA
       
      detA= mum(1,1)*mum(2,2)*mum(3,3) &
         + mum(2,1)*mum(3,2)*mum(1,3) &
         + mum(3,1)*mum(1,2)*mum(2,3) &
         - mum(1,1)*mum(3,2)*mum(2,3) &
         - mum(3,1)*mum(2,2)*mum(1,3) &
         - mum(2,1)*mum(1,2)*mum(3,3)

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

      END SUBROUTINE wmfem_inverse_tensor

      SUBROUTINE wmfem_disp_tensor(rho,ns,fmc)

      use wmfem_comm
      IMPLICIT NONE
      REAL(8),INTENT(IN):: rho
      INTEGER,INTENT(IN):: ns
      COMPLEX(8),DIMENSION(3,3,4,nfcmax2,nfcmax),INTENT(OUT):: fmc
      COMPLEX(8),DIMENSION(3,3):: fml

      REAL(8):: th,ph,dth,dph
      INTEGER:: nth,nph,mm,nn,i,j,nfc2,nfc
      COMPLEX(8):: cx

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
      use plprofw
      use plprof
      USE dptnsr0
      IMPLICIT NONE
      TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
      TYPE(pl_mag_type):: mag
      TYPE(pl_grd_type),DIMENSION(nsmax):: grd
      REAL(8),INTENT(IN):: rho,th,ph
      INTEGER,INTENT(IN):: mm,nn,ns
      COMPLEX(8),DIMENSION(3,3),INTENT(OUT):: fml
      COMPLEX(8):: cw,ckpara,ckperp
      COMPLEX(8):: ckppf,ckpps
      REAL(8):: babs_l,bsuprh,bsupth,bsupph
      INTEGER:: ns1

      DO ns1=1,nsmax
         grd(ns1)%grdn=0.D0
         grd(ns1)%grdtpr=0.D0
         grd(ns1)%grdtpp=0.D0
         grd(ns1)%grdu=0.D0
      END DO

      cw=2.d0*pi*crf*1.d6

      CALL pl_profw(rho,plfw)

      CALL wmfem_magnetic(rho,th,ph,babs_l,bsuprh,bsupth,bsupph)
      ckpara=mm*bsupth/babs_l+nn*bsupph/babs_l
      ckperp=(0.d0,0.d0)
      mag%babs=babs_l
      mag%bnx=0.D0
      mag%bny=0.D0
      mag%bnz=1.d0
      mag%rhon=rho

!$$$      if(abs(th).lt.1.d-4) then
!$$$         CALL DPCOLD_RKPERP(cw,ckpara,ckppf,ckpps)
!$$$         write(6,'(1P6E12.4)') rho,real(ckpara),ckppf,ckpps 
!$$$      endif

!$$$         if(ns.eq.1.and.abs(th).lt.0.1) 
!$$$     &    write(6,'(a,2i5,1p5e12.4)') 'm,n,kpara:',
!$$$     &         mm,nn,dble(ckpara),rho,babs_l,bsupth,bsupph

      CALL dp_tnsr0(cw,ckpara,ckperp,ns,mag,plfw,grd,fml)

      RETURN

      END SUBROUTINE wmfem_dielectric

    END MODULE wmcalc
