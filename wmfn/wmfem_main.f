
!     ***** wmfem main routine *****

      subroutine wmfem_main

      use wmfem_comm
      implicit none

!     ***** setup Fourier component index *****

      call wmfem_setup_index

!     ***** calculate coefficient matrix *****

      call wmfem_calculate_matrix

!     ***** setup boundary conditions *****

      call wmfem_boundary_condition_axis0
      call wmfem_boundary_condition_wall
      call wmfem_boundary_condition_fourier

!     ***** solve matrix equation *****
      call wmfem_solve

!     ***** calculate efield *****
      call wmfem_calculate_efield0
      call wmfem_efld

      return

      contains

!     ---- setup Fourier component index ----

      subroutine wmfem_setup_index

      integer:: nth,nhh,nfc,nfc2

!     setup an array of mode number

      if(nthmax.eq.1) then
         do nhh=1,nhhmax
            nfc=nhh
            nthnfc(nfc)=1
            mmnfc(nfc)=nth0
         enddo
      else
         do nhh=1,nhhmax
         do nth=1,nthmax
            nfc=nthmax*(nhh-1)+nth
            nthnfc(nfc)=nth
         enddo
         do nth=1,nthmax/2
            nfc=nthmax*(nhh-1)+nth
            mmnfc(nfc)=nth0+nth-1
            nfc=nthmax*(nhh-1)+nthmax+1-nth
            mmnfc(nfc)=nth0-nth
         enddo
         enddo
      endif

      if(nhhmax.eq.1) then
         do nth=1,nthmax
            nfc=nth
            nhhnfc(nfc)=1
            nnnfc(nfc)=nph0
         enddo
      else
         do nth=1,nthmax
         do nhh=1,nhhmax
            nfc=nthmax*(nhh-1)+nth
            nhhnfc(nfc)=nhh
         enddo
         do nhh=1,nhhmax/2
            nfc=nthmax*(nhh-1)+nth
            nnnfc(nfc)=nph0+nhh-1
            nfc=nthmax*(nhhmax-nhh)+nth
            nnnfc(nfc)=nph0-nhh
         enddo
         enddo
      endif

      do nhh=1,nhhmax2
         do nth=1,nthmax2
            nfc2=nthmax2*(nhh-1)+nth
            nthnfc2(nfc2)=nth
         enddo
         do nth=1,nthmax
            nfc2=nthmax2*(nhh-1)+nth
            mmnfc2(nfc2)=nth0+nth-1
            nfc=nthmax2*(nhh-1)+nthmax2+1-nth
            mmnfc2(nfc2)=nth0-nth
         enddo
      enddo
      
      do nth=1,nthmax2
         do nhh=1,nhhmax2
            nfc2=nthmax2*(nhh-1)+nth
            nhhnfc2(nfc2)=nhh
         enddo
         do nhh=1,nhhmax
            nfc2=nthmax2*(nhh-1)+nth
            nnnfc2(nfc2)=nph0+nhh-1
            nfc2=nthmax2*(nhhmax2-nhh)+nth
            nnnfc2(nfc2)=nph0-nhh
         enddo
      enddo
      return
      end subroutine wmfem_setup_index

!     ----- calculate coefficint matrix fma -----

      subroutine wmfem_calculate_matrix
      use libfem

      complex(8),dimension(:,:),allocatable:: fma_local
      complex(8),dimension(:,:,:),allocatable:: fvb1_local,fvb2_local
      integer:: ml,mw,ns,nfc,nth,nhh,nr,mb1,mb2,mm,nn,i
      integer:: mw1

      complex(8),dimension(nfcmax,nfcmax):: 
     &      fmvo_11,fmvo_12,fmvo_13,
     &      fmvo_21,fmvo_22,fmvo_23,
     &      fmvo_31,fmvo_32,fmvo_33
      integer::nfc2
      integer:: ml0,ml_p
      real(8)::DR
      complex(8)::cfactor
      integer :: nrd,ml1
      real(8) :: rd,angl,x
      complex(8) ::divj,div_A
!----- boundary conditions on axis -----

      cfactor=(2.d0*pi*crf*1.d6)

      allocate(fma_local(mwmax,mbmax))
      allocate(fvb1_local(nhhmax,nthmax,3))
      allocate(fvb2_local(nhhmax,nthmax,3))

      do ml=1,mlmax             ! clear global matrix fma
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         end do
      end do

      do ns=0,nsmax             ! loop for vacuum and plasma species

         do nr=1,nrmax-1        ! loop for elements
            print *, ns,nr
            call wmfem_calculate_local(nr,ns,fma_local)

            do mb2=1,mbmax
               do mw1=1,mwmax
                  fma(mw1,8*nfcmax*(nr-1)+mb2)
     &           =fma(mw1,8*nfcmax*(nr-1)+mb2)
     &           +fma_local(mw1,mb2)
               enddo
            enddo
         enddo
      enddo
      do ml=1,mlmax             ! clear RHS vector fvb
         fvb(ml)=0.d0
      enddo
       !!!!????????!!!!!!!

      do nr=1,nrmax-1
         write(6,*) 'nr=',nr
         call get_wmfvb(nr,  fvb1_local)
          call wmfem_ant_sub(rhoa(nr),
     &                 fmvo_11,fmvo_12,fmvo_13,
     &                 fmvo_21,fmvo_22,fmvo_23,
     &                 fmvo_31,fmvo_32,fmvo_33)

         do nfc=1,nfcmax
              nhh=nhhnfc(nfc)
              nth=nthnfc(nfc)
!          print *, fvb1_local(nhh,nth,1)
            do nfc2=1,nfcmax
              nhh=nhhnfc(nfc2)
              nth=nthnfc(nfc2)
              i=1 !
              ml=8*nfcmax*(nr-1)+i+8*(nfc-1)
!              fvb(ml)=fvb(ml)+ fvb1_local(nhh,nth,1)
              fvb(ml)=fvb(ml) 
     &               + fvb1_local(nhh,nth,1)*fmvo_11(nfc,nfc2)
     &               + fvb1_local(nhh,nth,2)*fmvo_12(nfc,nfc2)
     &               + fvb1_local(nhh,nth,3)*fmvo_13(nfc,nfc2)
              i=3 ! 
              ml=8*nfcmax*(nr-1)+i+8*(nfc-1)
!              fvb(ml)=fvb(ml)+fvb1_local(nhh,nth,2)
              fvb(ml)=fvb(ml) 
     &               + fvb1_local(nhh,nth,1)*fmvo_21(nfc,nfc2)
     &               + fvb1_local(nhh,nth,2)*fmvo_22(nfc,nfc2)
     &               + fvb1_local(nhh,nth,3)*fmvo_23(nfc,nfc2)
              i=5 !
              ml=8*nfcmax*(nr-1)+i+8*(nfc-1)
!             fvb(ml)=fvb(ml)+fvb1_local(nhh,nth,3)
              fvb(ml)=fvb(ml) 
     &               + fvb1_local(nhh,nth,1)*fmvo_31(nfc,nfc2)
     &               + fvb1_local(nhh,nth,2)*fmvo_32(nfc,nfc2)
     &               + fvb1_local(nhh,nth,3)*fmvo_33(nfc,nfc2)
            enddo

         enddo
      enddo
      
      deallocate(fma_local) 
      deallocate(fvb1_local) 
      deallocate(fvb2_local) 
      return
      end subroutine wmfem_calculate_matrix

!     ----- setup boundary conditions -----

!     -----     boundary condition on axis: circular -----

      subroutine wmfem_boundary_condition_axis0

      integer:: mc,mr,nr,mm,mll,ml,mw,ns,nfc,nn
      integer::nfc2,mm2,nn2
      integer::mw1o,mw2o,mw3o,mw4o,mw5o,mw6o,mw7o,mw8o

      complex(8):: cx,cy
      integer:: id_base=1

      integer::mlfactor

      complex(8),dimension(nfcmax,nfcmax):: 
     &      fmv_7_1,fmv_7_2,fmv_7_3,fmv_8_1,fmv_8_2,fmv_8_3
      complex(8),dimension(nfcmax,nfcmax):: 
     &      fmv_1,fmv_2,fmv_3,fmv_4,fmv_5,fmv_6
      complex(8)::cfactor
!----- boundary conditions on axis -----

      nr=1
!      cfactor=(2.d0*pi*crf*1.d6)*ci
      cfactor=(2.d0*pi*crf*1.d6)*ci/vc

      call wmfem_nabla_phi_sub(rhoa(1),
     &                  fmv_7_1,fmv_7_2,fmv_7_3,fmv_8_1,fmv_8_2,fmv_8_3)

      call wmfem_boundary_condition_div_sub(rhoa(1),
     &                             fmv_1,fmv_2,fmv_3,fmv_4,fmv_5,fmv_6)
!   --- Etheta = 0 for m=0 ---
      do nfc=1,nfcmax
         mm=mmnfc(nfc)
         nn=nnnfc(nfc)
         if (nn.ne.0)then
            ml=8*(nfc-1) + 7
            do mw=1,mwmax
              fma(mw,ml) = 0.d0
            end do
            fma(mwc,ml)=1d0
            fvb(ml)=0.d0
         endif

         if(mm.eq.0) then 

            do mlfactor =1,2
!            do mlfactor =1,4
!               if (mlfactor == 3)cycle
                ml=8*(nfc-1) + (mlfactor-1)*2 + 1
                do mw=1,mwmax
                   fma(mw,ml) = 0.d0
                enddo
                fma(mwc,ml)=1.d0
                fvb(ml)=0.d0
            enddo

            ml=8*(nfc-1) + 8
            do mw=1,mwmax
              fma(mw,ml) = 0.d0
            end do
            fma(mwc,ml)=1d0
            fvb(ml)=0.d0

         elseif(abs(mm).eq.1) then   ! 

!            ml =8*(nfc-1) + 5
            ml =8*(nfc-1) + 3
            do mw=1,mwmax-2
!            do mw=1,mwmax-4
               fma(mw,ml) = fma(mw,ml) 
     &                    + ci*mm*fma(mw+2,ml-2)
!     &                    + ci*mm*fma(mw+4,ml-4)
            end do
             fvb(ml)=fvb(ml)+ci*mm*fvb(ml-2)
!               fvb(ml)=fvb(ml)+ci*mm*fvb(ml-4)

            ml =8*(nfc-1) + 1
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
            end do

!           A,phi->E
            fma(mwc,ml)=cfactor
!            fma(mwc+4,ml)=cfactor*ci*mm
            fma(mwc+2,ml)=cfactor*ci*mm

            mw7o=mwc+6 -8*nfc
            mw8o=mwc+7 -8*nfc
!            do nfc2=1,nfcmax
!              fma(mw7o+8*nfc2,ml)=    fma(mw7o+8*nfc2,ml)
!     &                            + fmv_7_1(nfc2,nfc) 
!     &                            + fmv_7_2(nfc2,nfc)*ci*mm
!              fma(mw8o+8*nfc2,ml)=    fma(mw8o+8*nfc2,ml)
!     &                            + fmv_8_1(nfc2,nfc) 
!     &                            + fmv_8_2(nfc2,nfc)*ci*mm
!            enddo
            fvb(ml)=0.d0

!            ml =8*(nfc-1) + 4
!            do mw=1,mwmax-2
!               fma(mw,ml) = fma(mw,ml) 
!     &                    + ci*mm*fma(mw+2,ml-2)
!            end do
!            fvb(ml)=fvb(ml)+ci*mm*fvb(ml-2)

!          ml =8*(nfc-1) + 2
!          do mw=1,mwmax
!             fma(mw,ml) = 0.d0
!          end do
!           mw1o=mwc-1 -nfc*8
!           mw2o=mwc+0 -nfc*8
!           mw3o=mwc+1 -nfc*8
!           mw4o=mwc+2 -nfc*8
!           mw5o=mwc+3 -nfc*8
!           mw6o=mwc+4 -nfc*8
!           do nfc2=1,nfcmax
!             fma(mw1o+nfc2*8,ml)=fmv_1(nfc2,nfc)
!             fma(mw2o+nfc2*8,ml)=fmv_2(nfc2,nfc)
!             fma(mw3o+nfc2*8,ml)=fmv_3(nfc2,nfc)
!             fma(mw4o+nfc2*8,ml)=fmv_4(nfc2,nfc)
!             fma(mw5o+nfc2*8,ml)=fmv_5(nfc2,nfc)
!             fma(mw6o+nfc2*8,ml)=fmv_6(nfc2,nfc)
!         enddo
!         fvb(ml)=0.d0


            ml=8*(nfc-1) + 5
            do mw=1,mwmax
              fma(mw,ml) = 0.d0
            end do
            fma(mwc,ml)=1d0
            fvb(ml)=0.d0

            if ( nn .ne. 0)then
              ml=8*(nfc-1) + 8
              do mw=1,mwmax
                fma(mw,ml) = 0.d0
              end do
              fma(mwc,ml)=1d0
              fvb(ml)=0.d0
            endif
         else   
            do mlfactor =1,3
!             if (mlfactor == 3)cycle
              ml=8*(nfc-1) + (mlfactor-1)*2 + 1
              do mw=1,mwmax
                 fma(mw,ml) = 0.d0
              enddo
              fma(mwc,ml)=1.d0
              fvb(ml)=0.d0
            enddo

            ml=8*(nfc-1) + 7
            do mw=1,mwmax
              fma(mw,ml) = 0.d0
            end do
            fma(mwc,ml)=1d0
            fvb(ml)=0.d0


            ml=8*(nfc-1) + 8
            do mw=1,mwmax
              fma(mw,ml) = 0.d0
            end do
            fma(mwc,ml)=1d0
            fvb(ml)=0.d0
         endif
      enddo

      return
      end subroutine wmfem_boundary_condition_axis0

      subroutine wmfem_ant_sub(rho,
     &                 fmvo_11,fmvo_12,fmvo_13,
     &                 fmvo_21,fmvo_22,fmvo_23,
     &                 fmvo_31,fmvo_32,fmvo_33)
       IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(nthmax2,nhhmax2):: fv_11,fv_12, fv_13
      complex(8),dimension(nthmax2,nhhmax2):: fv_21,fv_22, fv_23
      complex(8),dimension(nthmax2,nhhmax2):: fv_31,fv_32, fv_33
      complex(8),dimension(nthmax2,nhhmax2):: fvf_11,fvf_12, fvf_13
      complex(8),dimension(nthmax2,nhhmax2):: fvf_21,fvf_22, fvf_23
      complex(8),dimension(nthmax2,nhhmax2):: fvf_31,fvf_32, fvf_33
      complex(8),dimension(nfcmax2)::
     &            fmv_11,fmv_12,fmv_13,fmv_21,fmv_22,fmv_23,
     &            fmv_31,fmv_32,fmv_33
      complex(8),dimension(nfcmax,nfcmax),intent(out)::
     &                  fmvo_11,fmvo_12,fmvo_13,
     &                  fmvo_21,fmvo_22,fmvo_23,
     &                  fmvo_31,fmvo_32,fmvo_33
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
      integer:: imn,imn1,imn2,nfcdiff

      real(8)::gj,dth,dph
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(nthmax2,nhhmax2):: gja
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      real(8),dimension(3,3) :: muminv 
      real(8),dimension(3,3) :: mum 
      real(8),dimension(3,3) :: gp
      integer::nphm,nphp,nthp,nthm

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)
      fv_11=0d0; fv_12=0d0; fv_13=0d0
      fv_21=0d0; fv_22=0d0; fv_23=0d0
      fv_31=0d0; fv_32=0d0; fv_33=0d0
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
           mum(:,:)=muma(:,:,nth,nph)
           gp(:,:)=gpa(:,:,nth,nph)
           call wmfem_inverse_tensor(mum,muminv)

           fv_11(nth,nph)= muminv(1,1)
           fv_12(nth,nph)= muminv(1,2)
           fv_13(nth,nph)= muminv(1,3)
           fv_21(nth,nph)= muminv(2,1)
           fv_22(nth,nph)= muminv(2,2)
           fv_23(nth,nph)= muminv(2,3)
           fv_31(nth,nph)= muminv(3,1)
           fv_32(nth,nph)= muminv(3,2)
           fv_33(nth,nph)= muminv(3,3)
       enddo
        
        call wmsubfx(fv_11,fvf_11,nthmax2,nhhmax2)
        call wmsubfx(fv_12,fvf_12,nthmax2,nhhmax2)
        call wmsubfx(fv_13,fvf_13,nthmax2,nhhmax2)
        call wmsubfx(fv_21,fvf_21,nthmax2,nhhmax2)
        call wmsubfx(fv_22,fvf_22,nthmax2,nhhmax2)
        call wmsubfx(fv_23,fvf_23,nthmax2,nhhmax2)
        call wmsubfx(fv_31,fvf_31,nthmax2,nhhmax2)
        call wmsubfx(fv_32,fvf_32,nthmax2,nhhmax2)
        call wmsubfx(fv_33,fvf_33,nthmax2,nhhmax2)
 
        do nfc2=1,nfcmax2
           nth=nthnfc2(nfc2)
           nph=nhhnfc2(nfc2)
           fmv_11(nfc2)=fvf_11(nth,nph)
           fmv_12(nfc2)=fvf_12(nth,nph)
           fmv_13(nfc2)=fvf_13(nth,nph)
           fmv_21(nfc2)=fvf_21(nth,nph)
           fmv_22(nfc2)=fvf_22(nth,nph)
           fmv_23(nfc2)=fvf_23(nth,nph)
           fmv_31(nfc2)=fvf_31(nth,nph)
           fmv_32(nfc2)=fvf_32(nth,nph)
           fmv_33(nfc2)=fvf_33(nth,nph)
        enddo


        do nfc1=1,nfcmax
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
              fmvo_11(nfc1,nfc2) = fmv_11(nfcdiff)
              fmvo_12(nfc1,nfc2) = fmv_12(nfcdiff)
              fmvo_13(nfc1,nfc2) = fmv_13(nfcdiff)
              fmvo_21(nfc1,nfc2) = fmv_21(nfcdiff)
              fmvo_22(nfc1,nfc2) = fmv_22(nfcdiff)
              fmvo_23(nfc1,nfc2) = fmv_23(nfcdiff)
              fmvo_31(nfc1,nfc2) = fmv_31(nfcdiff)
              fmvo_32(nfc1,nfc2) = fmv_32(nfcdiff)
              fmvo_33(nfc1,nfc2) = fmv_33(nfcdiff)
           enddo
        enddo
      end subroutine wmfem_ant_sub

      subroutine wmfem_ant_sub_r(rho,
     &                 fmvo_11,fmvo_12,fmvo_13,
     &                 fmvo_21,fmvo_22,fmvo_23,
     &                 fmvo_31,fmvo_32,fmvo_33)
       IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(nthmax2,nhhmax2):: fv_11,fv_12, fv_13
      complex(8),dimension(nthmax2,nhhmax2):: fv_21,fv_22, fv_23
      complex(8),dimension(nthmax2,nhhmax2):: fv_31,fv_32, fv_33
      complex(8),dimension(nthmax2,nhhmax2):: fvf_11,fvf_12, fvf_13
      complex(8),dimension(nthmax2,nhhmax2):: fvf_21,fvf_22, fvf_23
      complex(8),dimension(nthmax2,nhhmax2):: fvf_31,fvf_32, fvf_33
      complex(8),dimension(nfcmax2)::
     &            fmv_11,fmv_12,fmv_13,fmv_21,fmv_22,fmv_23,
     &            fmv_31,fmv_32,fmv_33
      complex(8),dimension(nfcmax,nfcmax),intent(out)::
     &                  fmvo_11,fmvo_12,fmvo_13,
     &                  fmvo_21,fmvo_22,fmvo_23,
     &                  fmvo_31,fmvo_32,fmvo_33
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
      integer:: imn,imn1,imn2,nfcdiff

      real(8)::gj,dth,dph
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(nthmax2,nhhmax2):: gja
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      real(8),dimension(3,3) :: muminv 
      real(8),dimension(3,3) :: mum 
      real(8),dimension(3,3) :: gp
      integer::nphm,nphp,nthp,nthm

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)
      fv_11=0d0; fv_12=0d0; fv_13=0d0
      fv_21=0d0; fv_22=0d0; fv_23=0d0
      fv_31=0d0; fv_32=0d0; fv_33=0d0
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
           mum(:,:)=muma(:,:,nth,nph)
           gp(:,:)=gpa(:,:,nth,nph)

           fv_11(nth,nph)= mum(1,1)
           fv_12(nth,nph)= mum(1,2)
           fv_13(nth,nph)= mum(1,3)
           fv_21(nth,nph)= mum(2,1)
           fv_22(nth,nph)= mum(2,2)
           fv_23(nth,nph)= mum(2,3)
           fv_31(nth,nph)= mum(3,1)
           fv_32(nth,nph)= mum(3,2)
           fv_33(nth,nph)= mum(3,3)
       enddo
        
        call wmsubfx(fv_11,fvf_11,nthmax2,nhhmax2)
        call wmsubfx(fv_12,fvf_12,nthmax2,nhhmax2)
        call wmsubfx(fv_13,fvf_13,nthmax2,nhhmax2)
        call wmsubfx(fv_21,fvf_21,nthmax2,nhhmax2)
        call wmsubfx(fv_22,fvf_22,nthmax2,nhhmax2)
        call wmsubfx(fv_23,fvf_23,nthmax2,nhhmax2)
        call wmsubfx(fv_31,fvf_31,nthmax2,nhhmax2)
        call wmsubfx(fv_32,fvf_32,nthmax2,nhhmax2)
        call wmsubfx(fv_33,fvf_33,nthmax2,nhhmax2)
 
        do nfc2=1,nfcmax2
           nth=nthnfc2(nfc2)
           nph=nhhnfc2(nfc2)
           fmv_11(nfc2)=fvf_11(nth,nph)
           fmv_12(nfc2)=fvf_12(nth,nph)
           fmv_13(nfc2)=fvf_13(nth,nph)
           fmv_21(nfc2)=fvf_21(nth,nph)
           fmv_22(nfc2)=fvf_22(nth,nph)
           fmv_23(nfc2)=fvf_23(nth,nph)
           fmv_31(nfc2)=fvf_31(nth,nph)
           fmv_32(nfc2)=fvf_32(nth,nph)
           fmv_33(nfc2)=fvf_33(nth,nph)
        enddo


        do nfc1=1,nfcmax
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
              fmvo_11(nfc1,nfc2) = fmv_11(nfcdiff)
              fmvo_12(nfc1,nfc2) = fmv_12(nfcdiff)
              fmvo_13(nfc1,nfc2) = fmv_13(nfcdiff)
              fmvo_21(nfc1,nfc2) = fmv_21(nfcdiff)
              fmvo_22(nfc1,nfc2) = fmv_22(nfcdiff)
              fmvo_23(nfc1,nfc2) = fmv_23(nfcdiff)
              fmvo_31(nfc1,nfc2) = fmv_31(nfcdiff)
              fmvo_32(nfc1,nfc2) = fmv_32(nfcdiff)
              fmvo_33(nfc1,nfc2) = fmv_33(nfcdiff)
           enddo
        enddo
      end subroutine wmfem_ant_sub_r

      subroutine wmfem_nabla_phi_sub(rho,
     &                 fmv_7_1,fmv_7_2,fmv_7_3,fmv_8_1,fmv_8_2,fmv_8_3)
       IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(nthmax2,nhhmax2):: fv_7_12, fv_7_13
      complex(8),dimension(nthmax2,nhhmax2):: fvf_7_12,fvf_7_13
      complex(8),dimension(nthmax2,nhhmax2):: fv_7_22, fv_7_23
      complex(8),dimension(nthmax2,nhhmax2):: fvf_7_22,fvf_7_23
      complex(8),dimension(nthmax2,nhhmax2):: fv_7_32, fv_7_33
      complex(8),dimension(nthmax2,nhhmax2):: fvf_7_32,fvf_7_33
      complex(8),dimension(nthmax2,nhhmax2):: fv_8_11,fv_8_21,fv_8_31
      complex(8),dimension(nthmax2,nhhmax2):: fvf_8_11,fvf_8_21,fvf_8_31
      complex(8),dimension(nfcmax2)::
     &            fmv_7_12,fmv_7_13,fmv_7_22,fmv_7_23,fmv_7_32,fmv_7_33,
     &            fmv_8_1d,fmv_8_2d,fmv_8_3d
      complex(8),dimension(nfcmax,nfcmax),intent(out)::
     &                  fmv_7_1,fmv_7_2,fmv_7_3,fmv_8_1,fmv_8_2,fmv_8_3
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
      integer:: imn,imn1,imn2,nfcdiff

      real(8)::gj,dth,dph
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(nthmax2,nhhmax2):: gja
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      real(8),dimension(3,3) :: muminv 
      real(8),dimension(3,3) :: mum 
      real(8),dimension(3,3) :: gp
      integer::nphm,nphp,nthp,nthm

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)
      fv_7_12=0d0; fv_7_13=0d0
      fv_7_22=0d0; fv_7_23=0d0
      fv_7_32=0d0; fv_7_33=0d0
      fv_8_11=0d0; fv_8_21=0d0 ; fv_8_31=0d0
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
           mum(:,:)=muma(:,:,nth,nph)
           gp(:,:)=gpa(:,:,nth,nph)
           call wmfem_inverse_tensor(mum,muminv)

!           print *,mum,muminv
           fv_7_12(nth,nph)= -muminv(1,2)
           fv_7_13(nth,nph)= -muminv(1,3)
           fv_7_22(nth,nph)= -muminv(2,2)
           fv_7_23(nth,nph)= -muminv(2,3)
           fv_7_32(nth,nph)= -muminv(3,2)
           fv_7_33(nth,nph)= -muminv(3,3)
           fv_8_11(nth,nph)= -muminv(1,1)
           fv_8_21(nth,nph)= -muminv(2,1)
           fv_8_31(nth,nph)= -muminv(3,1)
           if (rho < 1d-6)then
              fv_7_12(nth,nph)= 0d0
              fv_8_11(nth,nph)= -muminv(1,1) - muminv(1,2)*1d-6
              fv_7_22(nth,nph)= 0d0
              fv_8_21(nth,nph)= -muminv(2,1) - muminv(2,2)*1d-6
!              fv_7_32(nth,nph)= 0d0
!              fv_8_31(nth,nph)= -muminv(3,1) - muminv(3,2)*1d-6
           endif
       enddo
        
        call wmsubfx(fv_7_12,fvf_7_12,nthmax2,nhhmax2)
        call wmsubfx(fv_7_13,fvf_7_13,nthmax2,nhhmax2)
        call wmsubfx(fv_7_22,fvf_7_22,nthmax2,nhhmax2)
        call wmsubfx(fv_7_23,fvf_7_23,nthmax2,nhhmax2)
        call wmsubfx(fv_7_32,fvf_7_32,nthmax2,nhhmax2)
        call wmsubfx(fv_7_33,fvf_7_33,nthmax2,nhhmax2)
        call wmsubfx(fv_8_11,fvf_8_11,nthmax2,nhhmax2)
        call wmsubfx(fv_8_21,fvf_8_21,nthmax2,nhhmax2)
        call wmsubfx(fv_8_31,fvf_8_31,nthmax2,nhhmax2)
 
        do nfc2=1,nfcmax2
           nth=nthnfc2(nfc2)
           nph=nhhnfc2(nfc2)
           fmv_7_12(nfc2)=fvf_7_12(nth,nph)
           fmv_7_13(nfc2)=fvf_7_13(nth,nph)
           fmv_7_22(nfc2)=fvf_7_22(nth,nph)
           fmv_7_23(nfc2)=fvf_7_23(nth,nph)
           fmv_7_32(nfc2)=fvf_7_32(nth,nph)
           fmv_7_33(nfc2)=fvf_7_33(nth,nph)
           fmv_8_1d(nfc2)=fvf_8_11(nth,nph)
           fmv_8_2d(nfc2)=fvf_8_21(nth,nph)
           fmv_8_3d(nfc2)=fvf_8_31(nth,nph)
        enddo


        do nfc1=1,nfcmax
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
           fmv_7_1(nfc1,nfc2) = fmv_7_12(nfcdiff)*ci*mm2
     &                        + fmv_7_13(nfcdiff)*ci     *nn2
           fmv_7_2(nfc1,nfc2) = fmv_7_22(nfcdiff)*ci*mm2
     &                        + fmv_7_23(nfcdiff)*ci     *nn2
           fmv_7_3(nfc1,nfc2) = fmv_7_32(nfcdiff)*ci*mm2 
     &                        + fmv_7_33(nfcdiff)*ci     *nn2
           fmv_8_1(nfc1,nfc2) = fmv_8_1d(nfcdiff)
           fmv_8_2(nfc1,nfc2) = fmv_8_2d(nfcdiff)
           fmv_8_3(nfc1,nfc2) = fmv_8_3d(nfcdiff)
           enddo
        enddo
      end subroutine wmfem_nabla_phi_sub


!     -----     boundary condition on axis: non-circular -----

      subroutine wmfem_boundary_condition_axis1

      complex(8),dimension(3,3,nthmax2,nhhmax2):: mtxcl
!      complex(8),dimension(3*nfcmax,3*nfcmax):: mtxclx
      complex(8),dimension(:,:),ALLOCATABLE:: mtxclx
      complex(8),dimension(:,:),ALLOCATABLE:: mtx0,mtx1,mtx0i,mtx1i
      complex(8),dimension(:,:),ALLOCATABLE:: vec0,vec1,sol0,sol1
      integer,dimension(nfcmax):: nfc0a,nfc0b,ipos0
      integer,dimension(nfcmax):: nfc1a,nfc1b,ipos1
      real(8):: drho,rkth,rkph,rkth0,rho0,rho1,rho2,rho3,rho4,rhol,thl
      integer:: nfc1,nfc2,nn1,mm1,nn2,mm2,nn3,mm3,i,j,k,nfc,il,jl,kl,mc
      integer:: count0a,count0b,count1a,count1b
      integer:: count1,count2,counta,countb
      integer:: ierr

      allocate(mtxclx(3*nfcmax,3*nfcmax))

      CALL wmeq_get_mtxCL(nthmax2,nhhmax2,mtxcl)

      do nfc1=1,nfcmax
         do nfc2=1,nfcmax
            nn1=nnnfc(nfc1)
            mm1=mmnfc(nfc1)
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)
            nn3=nn1-nn2
            mm3=mm1-mm2
            IF(nn3.LT.0) nn3=nn3+nhhmax2
            IF(mm3.LT.0) mm3=mm3+nthmax2
            do i=1,3
               do j=1,3
                  mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
     &                 =mtxcl(i,j,mm3+1,nn3+1)
               enddo
            enddo
         enddo
      enddo

!      do i=1,3*nfcmax
!         do j=1,nfcmax
!            write(21,'(2I3,1P6E12.4)') i,3*(j-1)+1,
!     &           mtxclx(i,3*(j-1)+1),
!     &           mtxclx(i,3*(j-1)+2),
!     &           mtxclx(i,3*(j-1)+3)
!         enddo
!      enddo

!     ----- generate mm.ne.0 -----

      count0a=0
      count0b=0
      DO nfc=1,nfcmax
         IF(mmnfc(nfc).EQ.0) THEN
            count0b=count0b+1
            nfc0b(count0b)=nfc
            ipos0(nfc)=-count0b
         ELSE
            count0a=count0a+1
            nfc0a(count0a)=nfc
            ipos0(nfc)=count0a
         ENDIF
      ENDDO

      allocate(mtx0(3*count0a,3*count0a))
      allocate(mtx0i(3*count0a,3*count0a))
      allocate(vec0(3*count0a,3*count0b))
      allocate(sol0(3*count0a,3*count0b))


      DO nfc1=1,nfcmax
         IF(ipos0(nfc1).GT.0) THEN
            count1=ipos0(nfc1)
            DO nfc2=1,nfcmax
               IF(ipos0(nfc2).GT.0) THEN
                  count2=ipos0(nfc2)
                  DO i=1,3
                     DO j=1,3
                       mtx0(3*(count1-1)+i,3*(count2-1)+j)
     &                       =mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
                     ENDDO
                  ENDDO
               ELSE
                  count2=-ipos0(nfc2)
                  DO i=1,3
                     DO j=1,3
                       vec0(3*(count1-1)+i,3*(count2-1)+j)
     &                       =mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!      do i=1,3*count0a
!         do j=1,count0a
!            write(21,'(A,2I3,1P6E12.4)') 'M',i,3*(j-1)+1,
!     &           mtx0(i,3*(j-1)+1),
!     &           mtx0(i,3*(j-1)+2),
!     &           mtx0(i,3*(j-1)+3)
!         enddo
!      enddo
!      do j=1,3*count0b
!         do i=1,count0a
!            write(21,'(A,2I3,1P6E12.4)') 'V',3*(i-1)+1,j,
!     &           vec0(3*(i-1)+1,j),
!     &           vec0(3*(i-1)+2,j),
!     &           vec0(3*(i-1)+3,j)
!         enddo
!      enddo

      do i=1,3*count0a
         do j=1,3*count0a
            mtx0i(i,j)=mtx0(i,j)
         enddo
      enddo
      CALL invmcd(mtx0i,3*count0a,3*count0a,ierr)

      do j=1,3*count0b
         do i=1,3*count0a
            sol0(i,j)=(0.d0,0d0)
            do k=1,3*count0a
               sol0(i,j)=sol0(i,j)+mtx0i(i,k)*vec0(k,j)
            enddo
         enddo
      enddo

!      do j=1,3*count0b
!         do i=1,count0a
!            write(21,'(A,2I3,1P6E12.4)') 'S',3*(i-1)+1,j,
!     &           sol0(3*(i-1)+1,j),
!     &           sol0(3*(i-1)+2,j),
!     &           sol0(3*(i-1)+3,j)
!         enddo
!      enddo

!     ----- generate abs(mm).ne.1 -----

      count1a=0
      count1b=0
      DO nfc=1,nfcmax
         IF(ABS(mmnfc(nfc)).EQ.1) THEN
            count1b=count1b+1
            nfc1b(count1b)=nfc
            ipos1(nfc)=-count1b
         ELSE
            count1a=count1a+1
            nfc1a(count1a)=nfc
            ipos1(nfc)=count1a
         ENDIF
      ENDDO

      allocate(mtx1(3*count1a,3*count1a))
      allocate(mtx1i(3*count1a,3*count1a))
      allocate(vec1(3*count1a,3*count1b))
      allocate(sol1(3*count1a,3*count1b))

      DO nfc1=1,nfcmax
         IF(ipos1(nfc1).GT.0) THEN
            count1=ipos1(nfc1)
            DO nfc2=1,nfcmax
               IF(ipos1(nfc2).GT.0) THEN
                  count2=ipos1(nfc2)
                  DO i=1,3
                     DO j=1,3
                        mtx1(3*(count1-1)+i,3*(count2-1)+j)
     &                       =mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
                     ENDDO
                  ENDDO
               ELSE
                  count2=-ipos1(nfc2)
                  DO i=1,3
                     DO j=1,3
                        vec1(3*(count1-1)+i,3*(count2-1)+j)
     &                       =mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!      do i=1,3*count1a
!         do j=1,count1a
!            write(21,'(A,2I3,1P6E12.4)') 'M',i,3*(j-1)+1,
!     &           mtx1(i,3*(j-1)+1),
!     &           mtx1(i,3*(j-1)+2),
!     &           mtx1(i,3*(j-1)+3)
!         enddo
!      enddo
!      do j=1,3*count1b
!         do i=1,count1a
!            write(21,'(A,2I3,1P6E12.4)') 'V',3*(i-1)+1,j,
!     &           vec1(3*(i-1)+1,j),
!     &           vec1(3*(i-1)+2,j),
!     &           vec1(3*(i-1)+3,j)
!         enddo
!      enddo

      do i=1,3*count1a
         do j=1,3*count1a
            mtx1i(i,j)=mtx1(i,j)
         enddo
      enddo
      CALL invmcd(mtx1i,3*count1a,3*count1a,ierr)

      do j=1,3*count1b
         do i=1,3*count1a
            sol1(i,j)=(0.d0,0d0)
            do k=1,3*count1a
               sol1(i,j)=sol1(i,j)+mtx1i(i,k)*vec1(k,j)
            enddo
         enddo
      enddo

!      do j=1,3*count1b
!         do i=1,count1a
!            write(21,'(A,2I3,1P6E12.4)') 'S',3*(i-1)+1,j,
!     &           sol1(3*(i-1)+1,j),
!     &           sol1(3*(i-1)+2,j),
!     &           sol1(3*(i-1)+3,j)
!         enddo
!      enddo

!     ----- modify fma and fmv for E -----

      mc=(mwmax+1)/2

!     ----- add count0a to count0b for E -----
      
      DO countb=1,count0b
         nfc1=nfc0b(countb)
         do counta=1,count0a
            nfc2=nfc0a(counta)
            do i=1,3
               il=6*(nfc1-1)+2*(i-1)+1
               kl=6*(nfc2-1)+2*(i-1)+1
               do j=7,mwmax
!                  write(6,'(A,7I5)') 'nfc1,nfc2,i,il,kl,j,j+il-kl=',
!     &                 nfc1,nfc2,i,il,kl,j,j+il-kl
                  IF(j+il-kl.GE.1.AND.j+il-kl.LE.mwmax) THEN
                     fma(j,il)=fma(j,il)
     &                    +fma(j+il-kl,kl)*conjg(sol0(counta,countb))
                  ENDIF
               enddo
               fvb(il)=fvb(il)+fvb(kl)*conjg(sol0(counta,countb))
            enddo
         enddo
      enddo
               
!     ----- add count1a to count1b for E' -----
      
      DO countb=1,count1b
         nfc1=nfc1b(countb)
         do counta=1,count1a
            nfc2=nfc1a(counta)
            do i=1,3
               il=6*(nfc1-1)+2*(i-1)+2
               kl=6*(nfc2-1)+2*(i-1)+2
               do j=1,mwmax
                  IF(j+il-kl.GE.1.AND.j+il-kl.LE.mwmax) THEN
                  fma(j,il)=fma(j,il)
     &                     +fma(j+il-kl,kl)*conjg(sol1(counta,countb))
                  ENDIF
               enddo
               fvb(il)=fvb(il)+fvb(kl)*conjg(sol1(counta,countb))
            enddo
         enddo
      enddo

!     ----- replace fma for m=0 -----      
            
      DO counta=1,count0a
         nfc1=nfc0a(counta)
         do i=1,3
            il=6*(nfc1-1)+2*(i-1)+1
            do j=1,mwmax
               fma(j,il)=(0.d0,0.d0)
            enddo
            do nfc2=1,nfcmax
               do j=1,3
                  jl=6*(nfc2-1)+2*(j-1)+1
                  fma(jl-il+mc,il)=mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
               enddo
            enddo
            fvb(il)=(0.d0,0.d0)
         enddo
      enddo

!     ----- replace for abs(m)=1 -----      
            
      DO counta=1,count1a
         nfc1=nfc1a(counta)
         do i=1,3
            il=6*(nfc1-1)+2*(i-1)+2
            do j=1,mwmax
               fma(j,il)=(0.d0,0.d0)
            enddo
            do nfc2=1,nfcmax
               do j=1,3
                  jl=6*(nfc2-1)+2*(j-1)+2
                  fma(jl-il+mc,il)=mtxclx(3*(nfc1-1)+i,3*(nfc2-1)+j)
               enddo
            enddo
            fvb(il)=(0.d0,0.d0)
         enddo
      enddo

      deallocate(mtx0,mtx0i,vec0,sol0)
      deallocate(mtx1,mtx1i,vec1,sol1)

!      do i=1,6*nfcmax
!         do j=1,mwmax,3
!            write(21,'(A,2I3,1P6E12.4)') 
!     &           'f',i,j,fma(j,i),fma(j+1,i),fma(j+2,i)
!         enddo
!         write(21,'(A,I3,3X,1P2E12.4)') 
!     &           'v',i,fvb(i)
!      enddo

      deallocate(mtxclx)

      return
      end subroutine wmfem_boundary_condition_axis1

!     -----     boundary condition on wall -----

      subroutine wmfem_boundary_condition_wall

      IMPLICIT NONE
      integer:: nfc,ml,mw,ns,nr

!      complex(8),dimension(nfcmax):: fmv_1
!      complex(8),dimension(nfcmax):: fmv_2,fmv_4,fmv_6
      complex(8),dimension(nfcmax,nfcmax):: fmv_1,fmv_3,fmv_5
      complex(8),dimension(nfcmax,nfcmax):: fmv_2,fmv_4,fmv_6
      integer :: mw1o,mw2o,mw3o,mw4o,mw5o,mw6o
      integer ::nfc2

!      call wmfem_boundary_condition_wall_sub(fmv_1,fmv_2,fmv_4,fmv_6)
      nr=nrmax

      call wmfem_boundary_condition_div_sub(rhoa(nr),
     &                             fmv_1,fmv_2,fmv_3,fmv_4,fmv_5,fmv_6)

      do nfc=1,nfcmax
      
         ml=8*nfcmax*(nr-1)+8*(nfc-1) + 2
         do mw=1,mwmax
           fma(mw,ml) = 0.d0
         enddo

         mw1o=mwc-1 -nfc*8
         mw2o=mwc+0 -nfc*8
         mw3o=mwc+1 -nfc*8
         mw4o=mwc+2 -nfc*8
         mw5o=mwc+3 -nfc*8
         mw6o=mwc+4 -nfc*8
         do nfc2=1,nfcmax
           fma(mw1o+nfc2*8,ml)=fmv_1(nfc2,nfc)
           fma(mw2o+nfc2*8,ml)=fmv_2(nfc2,nfc)
           fma(mw3o+nfc2*8,ml)=fmv_3(nfc2,nfc)
           fma(mw4o+nfc2*8,ml)=fmv_4(nfc2,nfc)
           fma(mw5o+nfc2*8,ml)=fmv_5(nfc2,nfc)
           fma(mw6o+nfc2*8,ml)=fmv_6(nfc2,nfc)
         enddo
         fvb(ml)=0.d0

         ml=8*nfcmax*(nr-1)+8*(nfc-1) + 3
         do mw=1,mwmax
            fma(mw,ml) = 0.d0
         enddo
         fma(mwc,ml)=1.d0
         fvb(ml)=0.d0

        ml=8*nfcmax*(nr-1)+8*(nfc-1) + 5
        do mw=1,mwmax
            fma(mw,ml) = 0.d0
        enddo
        fma(mwc,ml)=1.d0
        fvb(ml)=0.d0

        ml=8*nfcmax*(nr-1)+8*(nfc-1) + 7
        do mw=1,mwmax
            fma(mw,ml) = 0.d0
         enddo
         fma(mwc,ml)=1.d0
         fvb(ml)=0.d0
      enddo
      return
      end subroutine wmfem_boundary_condition_wall
      subroutine wmfem_boundary_condition_div_sub(rho,
     &                             fmv_1,fmv_2,fmv_3,fmv_4,fmv_5,fmv_6)
       IMPLICIT NONE
      real(8),intent(in):: rho
      complex(8),dimension(nthmax2,nhhmax2):: fv_1, fv_1m, fv_1n
      complex(8),dimension(nthmax2,nhhmax2):: fvf_1,fvf_1m,fvf_1n
      complex(8),dimension(nthmax2,nhhmax2):: fv_3, fv_3m, fv_3n
      complex(8),dimension(nthmax2,nhhmax2):: fvf_3,fvf_3m,fvf_3n
      complex(8),dimension(nthmax2,nhhmax2):: fv_5, fv_5m, fv_5n
      complex(8),dimension(nthmax2,nhhmax2):: fvf_5,fvf_5m,fvf_5n
      complex(8),dimension(nthmax2,nhhmax2):: fv_2, fv_4, fv_6
      complex(8),dimension(nthmax2,nhhmax2):: fvf_2,fvf_4,fvf_6
      complex(8),dimension(nfcmax2):: fmv_1d,fmv_3d,fmv_5d
      complex(8),dimension(nfcmax2):: fmv_2d,fmv_4d,fmv_6d
      complex(8),dimension(nfcmax2):: fmv_1m,fmv_1n
      complex(8),dimension(nfcmax2):: fmv_3m,fmv_3n
      complex(8),dimension(nfcmax2):: fmv_5m,fmv_5n
      complex(8),dimension(nfcmax,nfcmax),intent(out)::fmv_1,fmv_3,fmv_5
      complex(8),dimension(nfcmax,nfcmax),intent(out)::fmv_2,fmv_4,fmv_6
      integer:: i,j,k,nfc,nfc1,nfc2
      integer:: nth,mm1,mm2,mmdiff,nph,nn1,nn2,nndiff
      integer:: imn,imn1,imn2,nfcdiff

      real(8)::gj,dth,dph
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gma,muma,dmuma
      real(8),dimension(3,3,nthmax2,nhhmax2) :: gpa,gmuma
      real(8),dimension(nthmax2,nhhmax2):: gja
      real(8),dimension(3,nthmax2,nhhmax2) :: dgjgmuma
      integer::nphm,nphp,nthp,nthm

      fmv_1=0d0
      fmv_3=0d0
      fmv_5=0d0
      fmv_2=0d0
      fmv_4=0d0
      fmv_6=0d0

      call wmfem_tensors(rho,gma,gpa,muma,dmuma,gja,gmuma,dgjgmuma)
      
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
           fv_1(nth,nph)= (dgjgmuma(1,nth,nph)  
     &                   +(gja(nthp,nph)*gmuma(2,1,nthp,nph)
     &                     -gja(nthm,nph)*gmuma(2,1,nthm,nph))/dth
     &                   +(gja(nth,nphp)*gmuma(3,1,nth,nphp)
     &                     -gja(nth,nphm)*gmuma(3,1,nth,nphm))/dph
     &                   )/gj

           fv_1m(nth,nph)= ci*gmuma(2,1,nth,nph)  
           fv_1n(nth,nph)= ci*gmuma(3,1,nth,nph)  

           fv_3(nth,nph)= (dgjgmuma(2,nth,nph)  
     &                   +(gja(nthp,nph)*gmuma(2,2,nthp,nph)
     &                     -gja(nthm,nph)*gmuma(2,2,nthm,nph))/dth
     &                   +(gja(nth,nphp)*gmuma(3,2,nth,nphp)
     &                     -gja(nth,nphm)*gmuma(3,2,nth,nphm))/dph
     &                   )/gj
 
           fv_3m(nth,nph)= ci*gmuma(2,2,nth,nph)  
           fv_3n(nth,nph)= ci*gmuma(3,2,nth,nph)  

           fv_5(nth,nph)= (dgjgmuma(3,nth,nph)  
     &                   +(gja(nthp,nph)*gmuma(2,3,nthp,nph)
     &                     -gja(nthm,nph)*gmuma(2,3,nthm,nph))/dth
     &                   +(gja(nth,nphp)*gmuma(3,3,nth,nphp)
     &                     -gja(nth,nphm)*gmuma(3,3,nth,nphm))/dph
     &                   )/gj
 
           fv_5m(nth,nph)= ci*gmuma(2,3,nth,nph)  
           fv_5n(nth,nph)= ci*gmuma(3,3,nth,nph)  

           fv_2(nth,nph)= gmuma(1,1,nth,nph)  
           fv_4(nth,nph)= gmuma(1,2,nth,nph)  
           fv_6(nth,nph)= gmuma(1,3,nth,nph) 
        enddo
        
        call wmsubfx(fv_1,fvf_1,nthmax2,nhhmax2)
        call wmsubfx(fv_1m,fvf_1m,nthmax2,nhhmax2)
        call wmsubfx(fv_1n,fvf_1n,nthmax2,nhhmax2)

        call wmsubfx(fv_3,fvf_3,nthmax2,nhhmax2)
        call wmsubfx(fv_3m,fvf_3m,nthmax2,nhhmax2)
        call wmsubfx(fv_3n,fvf_3n,nthmax2,nhhmax2)

        call wmsubfx(fv_5,fvf_5,nthmax2,nhhmax2)
        call wmsubfx(fv_5m,fvf_5m,nthmax2,nhhmax2)
        call wmsubfx(fv_5n,fvf_5n,nthmax2,nhhmax2)

        call wmsubfx(fv_2,fvf_2,nthmax2,nhhmax2)
        call wmsubfx(fv_4,fvf_4,nthmax2,nhhmax2)
        call wmsubfx(fv_6,fvf_6,nthmax2,nhhmax2)
 
        do nfc2=1,nfcmax2
           nth=nthnfc2(nfc2)
           nph=nhhnfc2(nfc2)
           fmv_1d(nfc2)=fvf_1(nth,nph)
           fmv_1m(nfc2)=fvf_1m(nth,nph)
           fmv_1n(nfc2)=fvf_1n(nth,nph)

           fmv_3d(nfc2)=fvf_3(nth,nph)
           fmv_3m(nfc2)=fvf_3m(nth,nph)
           fmv_3n(nfc2)=fvf_3n(nth,nph)

           fmv_5d(nfc2)=fvf_5(nth,nph)
           fmv_5m(nfc2)=fvf_5m(nth,nph)
           fmv_5n(nfc2)=fvf_5n(nth,nph)

           fmv_2d(nfc2)=fvf_2(nth,nph)
           fmv_4d(nfc2)=fvf_4(nth,nph)
           fmv_6d(nfc2)=fvf_6(nth,nph)
        enddo

        do nfc1=1,nfcmax
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
           fmv_1(nfc1,nfc2)= fmv_1d(nfcdiff) 
     &                + fmv_1m(nfcdiff)*mm2
     &                + fmv_1n(nfcdiff)     *nn2
           fmv_3(nfc1,nfc2)= fmv_3d(nfcdiff) 
     &                + fmv_3m(nfcdiff)*mm2
     &                + fmv_3n(nfcdiff)     *nn2
           fmv_5(nfc1,nfc2)= fmv_5d(nfcdiff) 
     &                + fmv_5m(nfcdiff)*mm2
     &                + fmv_5n(nfcdiff)     *nn2
           fmv_2(nfc1,nfc2)= fmv_2d(nfcdiff) 
           fmv_4(nfc1,nfc2)= fmv_4d(nfcdiff) 
           fmv_6(nfc1,nfc2)= fmv_6d(nfcdiff) 
           enddo
        enddo
      end subroutine wmfem_boundary_condition_div_sub


!     -----     boundary condition for fourier components -----

      subroutine wmfem_boundary_condition_fourier

!      integer:: mwc,nr,nn,mm,mll,ml,mw,ns,nfc,imax,i
      integer:: nr,nn,mm,mll,ml,mw,ns,nfc,imax,i

      if(nhhmax.gt.1) then
         nn=nhhmax/2+1
         do nr=1,nrmax
!!!!!!!!!!!!???????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            IF(nr.EQ.nrmax) THEN
!               imax=2
!            ELSE
!               imax=6
!            END IF
!!!!!!!!!!!!???????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do mm=1,nthmax
               mll=8*nthmax*(nn-1)+8*(mm-1)
               ml=8*nthmax*nhhmax*(nr-1)+mll
               do i=1,8
!               do i=1,imax
                  do mw=1,mwmax
                     fma(mw,ml+i)=0.d0
                  end do
                  fma(mwc,ml+i)=1.d0
                  fvb(ml+i)=0.d0
               end do
               if(mdlwmd.ge.1) then
                  do mw=1,mwmax
!                     do i=1,imax
                     do i=1,8
                        fma_save(mw,mll+i,nr,0)=fma(mw,ml+i)
                        DO ns=1,nsmax
                           fma_save(mw,mll+i,nr,ns)=0.d0
                        END DO
                     end do
                  end do
               endif
            end do
         end do
      endif

      if(nthmax.gt.1) then
         mm=nthmax/2+1
         do nr=1,nrmax
!!!!!!!!!!!!???????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            IF(nr.EQ.nrmax) THEN
!               imax=2
!            ELSE
!               imax=8
!            END IF
!!!!!!!!!!!!???????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do nn=1,nhhmax
               mll=8*nthmax*(nn-1)+8*(mm-1)
               ml=8*nthmax*nhhmax*(nr-1)+mll
!               do i=1,imax
               do i=1,8
                  do mw=1,mwmax
                     fma(mw,ml+i)=0.d0
                  end do
                  fma(mwc,ml+i)=1.d0
                  fvb(ml+i)=0.d0
               end do
               if(mdlwmd.ge.1) then
                  do mw=1,mwmax
!                     do i=1,imax
                     do i=1,8
                        fma_save(mw,mll+i,nr,0)=fma(mw,ml+i)
                        do ns=1,nsmax
                           fma_save(mw,mll+i,nr,ns)=0.d0
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end if

      return
      end subroutine wmfem_boundary_condition_fourier

!     ***** Solve matrix equation *****

      subroutine wmfem_solve

      integer:: ml,mw,ierr
      integer:: i,j,k

!     ----- solve matrix -----

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
         write(21,*),fvx(ml)
       end do

!       pause
!       do ml=1,mlmax
!          print *,ml,fvx(ml)
!          do mw = 1,mwmax
!          print *,mw,fma(mw,ml)
!          enddo
!       enddo 

!      write(6,*) 'mlmax,mwmax=',mlmax,mwmax
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      if(ierr.ne.0) write(6,*) '# ierr= ',ierr
      do ml=1,mlmax
         write(21,*),fvx(ml)
       end do

      return
      end subroutine wmfem_solve

!     ***** Calculate electric field *****

      subroutine wmfem_calculate_efield

      integer:: nr,nn,mm,ml,nfc,nth,nhh
      real(8):: drho


      nr=1
      drho=rhoa(2)-rhoa(1)
      do nfc=1,nfcmax
         nth=nthnfc(nfc)
         mm=mmnfc(nfc)
         nhh=nhhnfc(nfc)
         nn=nnnfc(nfc)
         ml=8*nfcmax*(nr-1)+nthmax*(nhh-1)+(nth-1)
         IF(abs(mm) == 1) THEN
            cef(1,nth,nhh,nr)=(9.D0*fvx(ml+2*nfcmax+1)
     &                        -1.D0*fvx(ml+5*nfcmax+1))/8.D0
            cdef(1,nth,nhh,nr)=0.D0
            cef(2,nth,nhh,nr)= fvx(ml+1)
            cdef(2,nth,nhh,nr)= 0.D0
         ELSE
            cef(1,nth,nhh,nr)=0.D0
            cdef(1,nth,nhh,nr)=fvx(ml+2*nfcmax+1)/(0.25D0*drho)
            cef(2,nth,nhh,nr)=0.D0
            cdef(2,nth,nhh,nr)=fvx(ml+3*nfcmax+1)/(0.5D0*drho)
         END IF
         IF(abs(mm) == 0) THEN
            cef(3,nth,nhh,nr)= fvx(ml+nfcmax+1)
            cdef(3,nth,nhh,nr)=0.D0
         ELSE
            cef(3,nth,nhh,nr)= 0.D0
            cdef(3,nth,nhh,nr)=fvx(ml+4*nfcmax+1)/(0.5D0*drho)
         END IF
      enddo

      do nr=2,nrmax-1
         drho=rhoa(nr+1)-rhoa(nr)
         do nfc=1,nfcmax
            nth=nthnfc(nfc)
            mm=mmnfc(nfc)
            nhh=nhhnfc(nfc)
            nn=nnnfc(nfc)
            ml=8*nfcmax*(nr-1)+nthmax*(nhh-1)+(nth-1)
            cef(1,nth,nhh,nr)=0.5d0*(fvx(ml+2*nfcmax+1)
     &                              +fvx(ml-  nfcmax+1))
            cdef(1,nth,nhh,nr)=(fvx(ml+2*nfcmax+1)
     &                         -fvx(ml-  nfcmax+1))/(0.5D0*drho)
            cef(2,nth,nhh,nr)=fvx(ml         +1)
            cdef(2,nth,nhh,nr)=(fvx(ml+3*nfcmax+1)
     &                         -fvx(ml-3*nfcmax+1))/drho
            cef(3,nth,nhh,nr)=fvx(ml+  nfcmax+1)
            cdef(3,nth,nhh,nr)=(fvx(ml+4*nfcmax+1)
     &                         -fvx(ml-2*nfcmax+1))/drho
         enddo
      enddo

      nr=nrmax
         drho=rhoa(nr)-rhoa(nr-1)
         do nfc=1,nfcmax
            nth=nthnfc(nfc)
            mm=mmnfc(nfc)
            nhh=nhhnfc(nfc)
            nn=nnnfc(nfc)
            ml=6*nfcmax*(nr-1)+nthmax*(nhh-1)+(nth-1)
            cef(1,nth,nhh,nr)=(3.D0*fvx(ml-  nfcmax+1)
     &                             -fvx(ml-4*nfcmax+1))/2.D0
            cdef(1,nth,nhh,nr)=(fvx(ml-  nfcmax+1)
     &                         -fvx(ml-4*nfcmax+1))/(0.5D0*drho)
            cef(2,nth,nhh,nr)=fvx(ml         +1)
            cdef(2,nth,nhh,nr)=(fvx(ml         +1)
     &                         -fvx(ml-3*nfcmax+1))/(0.5D0*drho)
            cef(3,nth,nhh,nr)=fvx(ml+  nfcmax+1)
            cdef(3,nth,nhh,nr)=(fvx(ml+  nfcmax+1)
     &                         -fvx(ml-2*nfcmax+1))/(0.5D0*drho)
         enddo
      return
      end subroutine wmfem_calculate_efield

      subroutine wmfem_calculate_efield0

      integer:: mc,mr,nr,mm,mll,ml,mw,ns,nfc,nn
      integer:: ml0,ml1,ml2
      integer::nfc2,nn2,mm2
      integer:: nth,nhh
      integer::mlo
      real(8):: drho

      complex(8):: cx,cy
      integer:: id_base=1

      complex(8),dimension(nfcmax,nfcmax)::
     &             fmv_7_1,fmv_7_2,fmv_7_3,fmv_8_1,fmv_8_2,fmv_8_3
      complex(8)::cfactor

      complex(8),dimension(nfcmax,nfcmax):: 
     &      fmv_1,fmv_2,fmv_3,fmv_4,fmv_5,fmv_6
      complex(8),dimension(nfcmax) ::div_A
      complex(8) ::div_A_local
      
      integer ::n_dev

      complex(8),dimension(nfcmax,nfcmax):: 
     &      fmvo_11,fmvo_12,fmvo_13,
     &      fmvo_21,fmvo_22,fmvo_23,
     &      fmvo_31,fmvo_32,fmvo_33
      complex(8),dimension(mlmax)::fvx_ef0
 
      cfactor=(2.d0*pi*crf*1.d6)*ci/vc
!      cfactor=(2.d0*pi*crf)*ci
      fvx_ef=0d0
      div_A=0d0
      do nr=1, nrmax
      call wmfem_nabla_phi_sub(rhoa(nr),
     &            fmv_7_1,fmv_7_2,fmv_7_3,fmv_8_1,fmv_8_2,fmv_8_3)

      call wmfem_boundary_condition_div_sub(rhoa(nr),
     &                             fmv_1,fmv_2,fmv_3,fmv_4,fmv_5,fmv_6)

      div_A=0d0
      mlo=8*nfcmax*(nr-1)
      do nfc=1,nfcmax
           mm=mmnfc(nfc)
           nn=nnnfc(nfc)

           ml0=8*nfcmax*(nr-1)+8*(nfc-1)
           ml1=8*nfcmax*(nr-1)+8*(nfc-1)
           ml2=6*nfcmax*(nr-1)+6*(nfc-1)
           div_A_local=0d0
           do nfc2=1,nfcmax
           div_A_local = div_A_local
     &         + fvx(mlo+8*(nfc2-1)+1)*fmv_1(nfc,nfc2)
     &         + fvx(mlo+8*(nfc2-1)+2)*fmv_2(nfc,nfc2)
     &         + fvx(mlo+8*(nfc2-1)+3)*fmv_3(nfc,nfc2)
     &         + fvx(mlo+8*(nfc2-1)+4)*fmv_4(nfc,nfc2)
     &         + fvx(mlo+8*(nfc2-1)+5)*fmv_5(nfc,nfc2)
     &         + fvx(mlo+8*(nfc2-1)+6)*fmv_6(nfc,nfc2)
           enddo
           div_A(nfc)=div_A(nfc) + div_A_local 
       n_dev=300+nfc
       write(n_dev, '(30es11.3,1x)')real(nr),div_A_local,
     &  fvx(ml1 + 1)*fmv_1(nfc,1),
     &  fvx(ml1 + 2)*fmv_2(nfc,1),
     &  fvx(ml1 + 3)*fmv_3(nfc,1),
     &  fvx(ml1 + 4)*fmv_4(nfc,1),
     &  fvx(ml1 + 5)*fmv_5(nfc,1),
     &  fvx(ml1 + 6)*fmv_6(nfc,1)
!     &   , fvx(ml1         )*fmv_1(nfc)
!     &           + fvx(ml1 + 1*nfcmax)*fmv_2(nfc)
!     &           + fvx(ml1 + 2*nfcmax)*fmv_3(nfc)
!     &           + fvx(ml1 + 3*nfcmax)*fmv_4(nfc)
!     &           + fvx(ml1 + 4*nfcmax)*fmv_5(nfc)
!     &           + fvx(ml1 + 5*nfcmax)*fmv_6(nfc)


       n_dev=200+nfc
       write(n_dev, '(18es15.6,1x)')real(nr),
!     & REAL(fvx(ml1 )*cfactor),IMAG(fvx(ml1 )*cfactor),
!     & REAL(fvx(ml1+1*nfcmax )*cfactor),
!     & IMAG(fvx(ml1+1*nfcmax )*cfactor),
!     & REAL(fvx(ml1+2*nfcmax )*cfactor),
!     & IMAG(fvx(ml1+2*nfcmax )*cfactor),
!     & REAL(fvx(ml1+3*nfcmax )*cfactor),
!     & IMAG(fvx(ml1+3*nfcmax )*cfactor),
!     & REAL(fvx(ml1+4*nfcmax )*cfactor),
!     & IMAG(fvx(ml1+4*nfcmax )*cfactor),
!     & REAL(fvx(ml1+5*nfcmax )*cfactor),
!     & IMAG(fvx(ml1+5*nfcmax )*cfactor),
     & REAL(fvx(ml1+1 )),IMAG(fvx(ml1+1 )),
     & REAL(fvx(ml1+2 )),
     & IMAG(fvx(ml1+2 )),
     & REAL(fvx(ml1+3 )),
     & IMAG(fvx(ml1+3 )),
     & REAL(fvx(ml1+4 )),
     & IMAG(fvx(ml1+4 )),
     & REAL(fvx(ml1+5 )),
     & IMAG(fvx(ml1+5 )),
     & REAL(fvx(ml1+6 )),
     & IMAG(fvx(ml1+6 )),
     & REAL(fvx(ml1+7 )),
     & IMAG(fvx(ml1+7 )),
     & REAL(fvx(ml1+8 )),
     & IMAG(fvx(ml1+8 ))

           ml1=8*nfcmax*(nr-1)+8*(nfc-1)
           ml2=6*nfcmax*(nr-1)+6*(nfc-1)
           fvx_ef(ml2+ 1 )=cfactor*fvx(ml1+ 1) 
           fvx_ef(ml2+ 3 )=cfactor*fvx(ml1+ 3) 
           fvx_ef(ml2+ 5 )=cfactor*fvx(ml1+ 5) 

           do nfc2 = 1, nfcmax
              ml0=8*nfcmax*(nr-1)+8*(nfc2-1)
              fvx_ef(ml2+ 1 )=fvx_ef(ml2+ 1 )
     &                      + fmv_8_1(nfc,nfc2)*fvx(ml0+ 8)
     &                      + fmv_7_1(nfc,nfc2)*fvx(ml0+ 7)
              fvx_ef(ml2+ 3 )=fvx_ef(ml2+ 3 )
     &                    + fmv_8_2(nfc,nfc2)*fvx(ml0+ 8)
     &                    + fmv_7_2(nfc,nfc2)*fvx(ml0+ 7)
              fvx_ef(ml2+ 5 )=fvx_ef(ml2+ 5 )
     &                    + fmv_8_3(nfc,nfc2)*fvx(ml0+ 8)
     &                    + fmv_7_3(nfc,nfc2)*fvx(ml0+ 7)
!              print *,rhoa(nr),fvx(ml0+ 7),fvx(ml0+ 8),
!     &  fmv_8_1(nfc,nfc2)*fvx(ml0+ 8) + fmv_7_1(nfc,nfc2)*fvx(ml0+ 7) ,
!     &  fmv_8_2(nfc,nfc2)*fvx(ml0+ 8) + fmv_7_2(nfc,nfc2)*fvx(ml0+ 7) ,
!     &  fmv_8_3(nfc,nfc2)*fvx(ml0+ 8) + fmv_7_3(nfc,nfc2)*fvx(ml0+ 7) 
           enddo   
       enddo
       enddo
       write(36, '(30es11.3,1x)')
       write(36, '(30es11.3,1x)')

      fvx_ef0=fvx_ef
      do nr=1,nrmax

         call wmfem_ant_sub_r(rhoa(nr),
     &                 fmvo_11,fmvo_12,fmvo_13,
     &                 fmvo_21,fmvo_22,fmvo_23,
     &                 fmvo_31,fmvo_32,fmvo_33)
        fvx_ef=0d0
!!        mlo=6*nfcmax*(nr-1)
         mlo=8*nfcmax*(nr-1)
         do nfc=1,nfcmax
            ml0=6*nfcmax*(nr-1)+6*(nfc-1)
             do nfc2=1,nfcmax
                fvx_ef(ml0+1)= fvx_ef0(ml0+1)
!                fvx_ef(ml0+1)= fvx_ef(ml0+1)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+1)*fmvo_11(nfc,nfc2)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+3)*fmvo_12(nfc,nfc2)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+5)*fmvo_13(nfc,nfc2)

                fvx_ef(ml0+3)= fvx_ef0(ml0+3)
!                fvx_ef(ml0+3)= fvx_ef(ml0+3)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+1)*fmvo_21(nfc,nfc2)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+3)*fmvo_22(nfc,nfc2)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+5)*fmvo_23(nfc,nfc2)

                fvx_ef(ml0+5)= fvx_ef0(ml0+5)
!                fvx_ef(ml0+5)= fvx_ef(ml0+5)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+1)*fmvo_31(nfc,nfc2)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+3)*fmvo_32(nfc,nfc2)
!     &                     + fvx_ef0(mlo+6*(nfc2-1)+5)*fmvo_33(nfc,nfc2)

             enddo
            
            nth=nthnfc(nfc)
            mm=mmnfc(nfc)
            nhh=nhhnfc(nfc)
            nn=nnnfc(nfc)
            ml=6*nfcmax*(nr-1)+6*nthmax*(nhh-1)+6*(nth-1)
            cef(1,nth,nhh,nr)=fvx_ef(ml+1)
            cef(2,nth,nhh,nr)=fvx_ef(ml+3)
            cef(3,nth,nhh,nr)=fvx_ef(ml+5)
!            cef(3,nth,nhh,nr)= fvx(mlo+ 7)
            print *,nr , fvx(mlo+ 7),fvx(mlo+ 8)
         enddo
      enddo
      return
      end subroutine wmfem_calculate_efield0
      end subroutine wmfem_main
