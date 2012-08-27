
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

      call wmfem_calculate_efield
      call wmfem_efld

      return

      contains

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

      complex(8),dimension(:,:),allocatable:: fma_local
      complex(8),dimension(:,:,:),allocatable:: fvb1_local,fvb2_local
      integer:: ml,mw,ns,nfc,nth,nph,nr,mb1,mb2,mm,nn,i

      allocate(fma_local(mbmax,mbmax))
      allocate(fvb1_local(nphmax,nthmax,3))
      allocate(fvb2_local(nphmax,nthmax,3))

      do ml=1,mlmax             ! clear global matrix fma
         do mw=1,mwmax
            fma(mw,ml)=0.d0
         end do
      end do

      do ns=0,nsmax             ! loop for vacuum and plasma species

         do nr=1,nrmax-1        ! loop for elements

            call wmfem_calculate_local(nr,ns,fma_local)

            if(mdlwmd.ge.1) then
               do mb2=1,mbmax
                  do mb1=1,mbmax
                     fma_save(mb1,mb2,nr,ns)=fma_local(mb1,mb2)
                  enddo
               enddo
            endif

            do mb2=1,mbmax
               do mb1=1,mbmax
                  fma(mwc+mb1-mb2,6*nfcmax*(nr-1)+mb2)
     &           =fma(mwc+mb1-mb2,6*nfcmax*(nr-1)+mb2)
     &           +fma_local(mb1,mb2)
               enddo
            enddo
         enddo
      enddo

      do ml=1,mlmax             ! clear RHS vector fvb
         fvb(ml)=0.d0
      enddo
         
      do nr=1,nrmax-1
         call get_wmfvb(nr,  fvb1_local)
         call get_wmfvb(nr+1,fvb2_local)
         do nfc=1,nfcmax
            nph=nphnfc(nfc)
            nth=nthnfc(nfc)
            i=1 ! E_perp (nr)
            ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
            fvb(ml)=fvb(ml)+fvb1_local(nph,nth,2)
            i=2 ! E_para (nr)
            ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
            fvb(ml)=fvb(ml)+fvb1_local(nph,nth,3)
            i=3 ! E_rho (nr+1/4)
            ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
            fvb(ml)=fvb(ml)+0.75D0*fvb1_local(nph,nth,1)
     &                     +0.25D0*fvb2_local(nph,nth,1)
            i=4 ! E_perp (nr+1/2)
            ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
            fvb(ml)=fvb(ml)+0.50D0*fvb1_local(nph,nth,2)
     &                     +0.50D0*fvb2_local(nph,nth,2)
            i=5 ! E_para (nr+1/2)
            ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
            fvb(ml)=fvb(ml)+0.50D0*fvb1_local(nph,nth,3)
     &                     +0.50D0*fvb2_local(nph,nth,3)
            i=6 ! E_rho (nr+3/4)
            ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
            fvb(ml)=fvb(ml)+0.25D0*fvb1_local(nph,nth,1)
     &                     +0.75D0*fvb2_local(nph,nth,1)
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

      complex(8):: cx,cy
      integer:: id_base=1

!----- boundary conditions on axis -----

      nr=1

!   --- Etheta = 0 for m=0 ---

      do nfc=1,nfcmax
         mm=mmnfc(nfc)
         if(mm.eq.0) then   ! E_perp(0)=0.D0
            ml=nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
            enddo
            fma(mwc,ml)=1.d0
            fvb(ml)=0.d0
            if(mdlwmd.ge.1) then
               do mw=1,mwmax
                  fma_save(mw,ml,nr,0)=fma(mw,ml)
                  do ns=1,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               enddo
            end if

         elseif(abs(mm).eq.1) then   ! E+(0)=0.D0, E_para(0)=0.D0
            ml =nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
            end do
            fma(mwc,ml)=1.d0
            fma(mwc+2*nfcmax,ml)= 1.5D0*ci*mm
            fma(mwc+5*nfcmax,ml)=-0.5D0*ci*mm
            fvb(ml)=0.d0
            if(mdlwmd.ge.1) then
               do mw=1,mwmax
                  fma_save(mw,ml,nr,0)=fma(mw,ml)
                  do ns=1,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               enddo
            end if

            ml=nfcmax+nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
            enddo
            fma(mwc,ml)=1.d0
            fvb(ml)=0.d0
            if(mdlwmd.ge.1) then
               do mw=1,mwmax
                  fma_save(ml,ml,nr,0)=fma(mw,ml)
                  do ns=1,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end do
            end if

         else   ! Eperp(0)=0.D0, E_para(0)=0.D0
            ml=nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
            enddo
            fma(mwc,ml)=1.d0
            fvb(ml)=0.d0
            if(mdlwmd.ge.1) then
               do mw=1,mwmax
                  fma_save(ml,ml,nr,0)=fma(mw,ml)
                  do ns=1,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end do
            end if

            ml=nfcmax+nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
            enddo
            fma(mwc,ml)=1.d0
            fvb(ml)=0.d0
            if(mdlwmd.ge.1) then
               do mw=1,mwmax
                  fma_save(ml,ml,nr,0)=fma(mw,ml)
                  do ns=1,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end do
            end if

         endif
      enddo

      return
      end subroutine wmfem_boundary_condition_axis0

!     -----     boundary condition on axis: non-circular -----

      subroutine wmfem_boundary_condition_axis1

      complex(8),dimension(3,3,nthmax2,nphmax2):: mtxcl
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

      CALL wmeq_get_mtxCL(nthmax2,nphmax2,mtxcl)

      do nfc1=1,nfcmax
         do nfc2=1,nfcmax
            nn1=nnnfc(nfc1)
            mm1=mmnfc(nfc1)
            nn2=nnnfc(nfc2)
            mm2=mmnfc(nfc2)
            nn3=nn1-nn2
            mm3=mm1-mm2
            IF(nn3.LT.0) nn3=nn3+nphmax2
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

      integer:: nfc,ml,mw,ns

      do nfc=1,nfcmax
         ml=6*nfcmax*(nrmax-1)+nfc
         do mw=1,mwmax
            fma(mw,ml) = 0.d0
         enddo
         fma(mwc,ml)=1.d0
         fvb(ml)=0.d0
         if(mdlwmd.ge.1) then
            do mw=1,mwmax
               fma_save(ml,ml,nrmax,0)=fma(mw,ml)
               do ns=1,nsmax
                  fma_save(mw,ml,nrmax,ns)=0.d0
               enddo
            end do
         end if

         ml=6*nfcmax*(nrmax-1)+nfcmax+nfc
         do mw=1,mwmax
            fma(mw,ml) = 0.d0
         enddo
         fma(mwc,ml)=1.d0
         fvb(ml)=0.d0
         if(mdlwmd.ge.1) then
            do mw=1,mwmax
               fma_save(ml,ml,nrmax,0)=fma(mw,ml)
               do ns=1,nsmax
                  fma_save(mw,ml,nrmax,ns)=0.d0
               enddo
            end do
         end if
      enddo

      return
      end subroutine wmfem_boundary_condition_wall

!     -----     boundary condition for fourier components -----

      subroutine wmfem_boundary_condition_fourier

      integer:: mwc,nr,nn,mm,mll,ml,mw,ns,nfc,imax,i

      if(nphmax.gt.1) then
         nn=nphmax/2+1
         do nr=1,nrmax
            IF(nr.EQ.nrmax) THEN
               imax=2
            ELSE
               imax=6
            END IF
            do mm=1,nthmax
               mll=6*nthmax*(nn-1)+6*(mm-1)
               ml=6*nthmax*nphmax*(nr-1)+mll
               do i=1,imax
                  do mw=1,mwmax
                     fma(mw,ml+i)=0.d0
                  end do
                  fma(mwc,ml+i)=1.d0
                  fvb(ml+i)=0.d0
               end do
               if(mdlwmd.ge.1) then
                  do mw=1,mwmax
                     do i=1,imax
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
            IF(nr.EQ.nrmax) THEN
               imax=2
            ELSE
               imax=6
            END IF
            do nn=1,nphmax
               mll=6*nthmax*(nn-1)+6*(mm-1)
               ml=6*nthmax*nphmax*(nr-1)+mll
               do i=1,imax
                  do mw=1,mwmax
                     fma(mw,ml+i)=0.d0
                  end do
                  fma(mwc,ml+i)=1.d0
                  fvb(ml+i)=0.d0
               end do
               if(mdlwmd.ge.1) then
                  do mw=1,mwmax
                     do i=1,imax
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

      integer:: nr,nn,mm,ml,nfc,nth,nph
      real(8):: drho

      nr=1
      drho=rhoa(2)-rhoa(1)
      do nfc=1,nfcmax
         nth=nthnfc(nfc)
         mm=mmnfc(nfc)
         nph=nphnfc(nfc)
         nn=nnnfc(nfc)
         ml=6*nfcmax*(nr-1)+nthmax*(nph-1)+(nth-1)
         IF(abs(mm) == 1) THEN
            cef(1,nth,nph,nr)=(9.D0*fvx(ml+2*nfcmax+1)
     &                        -1.D0*fvx(ml+5*nfcmax+1))/8.D0
            cdef(1,nth,nph,nr)=0.D0
            cef(2,nth,nph,nr)= fvx(ml+1)
            cdef(2,nth,nph,nr)= 0.D0
         ELSE
            cef(1,nth,nph,nr)=0.D0
            cdef(1,nth,nph,nr)=fvx(ml+2*nfcmax+1)/(0.25D0*drho)
            cef(2,nth,nph,nr)=0.D0
            cdef(2,nth,nph,nr)=fvx(ml+3*nfcmax+1)/(0.5D0*drho)
         END IF
         IF(abs(mm) == 0) THEN
            cef(3,nth,nph,nr)= fvx(ml+nfcmax+1)
            cdef(3,nth,nph,nr)=0.D0
         ELSE
            cef(3,nth,nph,nr)= 0.D0
            cdef(3,nth,nph,nr)=fvx(ml+4*nfcmax+1)/(0.5D0*drho)
         END IF
      enddo

      do nr=2,nrmax-1
         drho=rhoa(nr+1)-rhoa(nr)
         do nfc=1,nfcmax
            nth=nthnfc(nfc)
            mm=mmnfc(nfc)
            nph=nphnfc(nfc)
            nn=nnnfc(nfc)
            ml=6*nfcmax*(nr-1)+nthmax*(nph-1)+(nth-1)
            cef(1,nth,nph,nr)=0.5d0*(fvx(ml+2*nfcmax+1)
     &                              +fvx(ml-  nfcmax+1))
            cdef(1,nth,nph,nr)=(fvx(ml+2*nfcmax+1)
     &                         -fvx(ml-  nfcmax+1))/(0.5D0*drho)
            cef(2,nth,nph,nr)=fvx(ml         +1)
            cdef(2,nth,nph,nr)=(fvx(ml+3*nfcmax+1)
     &                         -fvx(ml-3*nfcmax+1))/drho
            cef(3,nth,nph,nr)=fvx(ml+  nfcmax+1)
            cdef(3,nth,nph,nr)=(fvx(ml+4*nfcmax+1)
     &                         -fvx(ml-2*nfcmax+1))/drho
         enddo
      enddo

      nr=nrmax
         drho=rhoa(nr)-rhoa(nr-1)
         do nfc=1,nfcmax
            nth=nthnfc(nfc)
            mm=mmnfc(nfc)
            nph=nphnfc(nfc)
            nn=nnnfc(nfc)
            ml=6*nfcmax*(nr-1)+nthmax*(nph-1)+(nth-1)
            cef(1,nth,nph,nr)=(3.D0*fvx(ml-  nfcmax+1)
     &                             -fvx(ml-4*nfcmax+1))/2.D0
            cdef(1,nth,nph,nr)=(fvx(ml-  nfcmax+1)
     &                         -fvx(ml-4*nfcmax+1))/(0.5D0*drho)
            cef(2,nth,nph,nr)=fvx(ml         +1)
            cdef(2,nth,nph,nr)=(fvx(ml         +1)
     &                         -fvx(ml-3*nfcmax+1))/(0.5D0*drho)
            cef(3,nth,nph,nr)=fvx(ml+  nfcmax+1)
            cdef(3,nth,nph,nr)=(fvx(ml+  nfcmax+1)
     &                         -fvx(ml-2*nfcmax+1))/(0.5D0*drho)
         enddo
      return
      end subroutine wmfem_calculate_efield
      end subroutine wmfem_main
