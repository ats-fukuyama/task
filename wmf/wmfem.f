!     $Id$

!     ----- allocate arrays -----

      subroutine wmfem_allocate

      use wmfem_comm
      implicit none
      integer,save:: nrmax_save=0,nthmax_save=0,nphmax_save=0
      integer,save:: mwmax_save=0,mlmax_save=0,mbmax_save=0
      integer,save:: nsmax_save=0,nfcmax_save=0
      integer,save:: mdlwmd_save=0 

      if((nrmax.ne.nrmax_save).or.(nthmax.ne.nthmax_save).or. 
     &   (nphmax.ne.nphmax_save)) then
         if(ALLOCATED(cef)) deallocate(cef)
         allocate(cef(3,nthmax,nphmax,nrmax))
         if(ALLOCATED(cdef)) deallocate(cdef)
         allocate(cdef(3,nthmax,nphmax,nrmax))
         if(ALLOCATED(cbf)) deallocate(cbf)
         allocate(cbf(3,nthmax,nphmax,nrmax))
         if(ALLOCATED(cpp)) deallocate(cpp)
         allocate(cpp(nthmax,nphmax,nthmax2,nphmax2,nrmax,0:nsmax))
         if(ALLOCATED(cpa)) deallocate(cpa)
         allocate(cpa(nthmax,nphmax))
      endif

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(ALLOCATED(fma)) deallocate(fma)
         allocate(fma(mwmax,mlmax))
      endif

      if(mdlwmd.ge.1) then
      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save).or.
     &   (nsmax.ne.nsmax_save).or.mdlwmd.ne.mdlwmd_save) then
         if(ALLOCATED(fma_save)) deallocate(fma_save)
         if(mdlwmd.ge.1) then
            allocate(fma_save(mbmax,mbmax,nrmax,0:nsmax))
         endif 
      endif
      endif

      if(mlmax.ne.mlmax_save) then
         if(ALLOCATED(fvb)) deallocate(fvb)
         if(ALLOCATED(fvx)) deallocate(fvx)
         allocate(fvb(mlmax))
         allocate(fvx(mlmax))
      endif

      if(nfcmax.ne.nfcmax_save) then
         if(ALLOCATED(nthnfc)) deallocate(nthnfc)
         if(ALLOCATED(nphnfc)) deallocate(nphnfc)
         if(ALLOCATED(mmnfc)) deallocate(mmnfc)
         if(ALLOCATED(nnnfc)) deallocate(nnnfc)
         if(nfcmax.ne.0) allocate(nthnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nphnfc(nfcmax))
         if(nfcmax.ne.0) allocate(mmnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nnnfc(nfcmax))

         if(ALLOCATED(nthnfc2)) deallocate(nthnfc2)
         if(ALLOCATED(nphnfc2)) deallocate(nphnfc2)
         if(ALLOCATED(mmnfc2)) deallocate(mmnfc2)
         if(ALLOCATED(nnnfc2)) deallocate(nnnfc2)
         if(nfcmax2.ne.0) allocate(nthnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nphnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(mmnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nnnfc2(nfcmax2))
      endif
      mwmax_save=mwmax
      mlmax_save=mlmax
      nsmax_save=nsmax
      nfcmax_save=nfcmax
      mbmax_save=mbmax
      mdlwmd_save=mdlwmd
      end subroutine wmfem_allocate

!     ***** wmfem pre routine *****

      subroutine wmfem_pre

      use wmfem_comm
      implicit none
      integer:: ierr

!     ***** metric setup  *****

         CALL wmfem_setg(ierr)
         IF(IERR.NE.0) RETURN
         CALL wmfem_setj(ierr)
         IF(IERR.NE.0) RETURN
         CALL get_wmfem_size(nrmax,nthmax,nphmax,nsmax)

!     ***** define array size  *****

      nfcmax=nthmax*nphmax      ! size of block matrix 
                                !    (number of Fourier components)
      mlmax=nfcmax*(6*nrmax-4)  ! length of coeffient matrix and source vector
                                !   Eperp0,Epara0,Erho1/4,
                                !   Eperp1/2,Epara1/2,Erho3/4

      mbmax=nfcmax*8            ! size of block matrix
      mwmax=2*mbmax-1           ! width of coefficient matrix
      mwc=mbmax                 ! position of diagonal coponent

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

      return
      end subroutine wmfem_pre

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
!      call wmfem_boundary_condition_axis1
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

      do ml=1,mlmax             ! clear fma
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

!------      fvb

      do ml=1,mlmax
         fvb(ml)=0.d0
      enddo
         
      do nr=1,nrmax-1
         call get_wmfvb(nr,  fvb1_local)
         call get_wmfvb(nr+1,fvb2_local)
         do i=1,3
            do nfc=1,nfcmax
               nph=nphnfc(nfc)
               nth=nthnfc(nfc)
               ml=6*nfcmax*(nr-1)+nfcmax*(i-1)+nfc
               fvb(ml)=fvb(ml)+fvb1_local(nph,nth,i)
               ml=ml+3*nfcmax
               fvb(ml)=fvb(ml)+fvb2_local(nph,nth,i)
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

      complex(8):: cx,cy
      integer:: id_base=1

!----- boundary conditions on axis -----

      nr=1

!   --- Etheta = 0 for m=0 ---

      do nfc=1,nfcmax
         mm=mmnfc(nfc)
         if(mm.eq.0) then
            ml=nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end if
            enddo
            fma(mwc,ml)=1.d0
            if(mdlwmd.ge.1) then
               fma_save(ml,ml,nr,0)=1.d0
            end if
            fvb(ml)=0.d0

         elseif(abs(mm).eq.1) then
            ml =nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end if
            end do
            fma(mwc,ml)=1.d0
            fma(mwc+2*nfcmax,ml)= 1.5D0*ci*mm
            fma(mwc+5*nfcmax,ml)=-0.5D0*ci*mm
            if(mdlwmd.ge.1) then
               fma_save(ml,ml,nr,0)= 1.d0
               fma_save(ml+2*nfcmax,ml,nr,0)= 1.5D0*ci*mm
               fma_save(ml+5*nfcmax,ml,nr,0)=-0.5D0*ci*mm
            end if
            fvb(ml)=0.d0

            ml=nfcmax+nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end if
            enddo
            fma(mwc,ml)=1.d0
            if(mdlwmd.ge.1) then
               fma_save(ml,ml,nr,0)=1.d0
            end if
            fvb(ml)=0.d0

         else
            ml=nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end if
            enddo
            fma(mwc,ml)=1.d0
            if(mdlwmd.ge.1) then
               fma_save(ml,ml,nr,0)=1.d0
            end if
            fvb(ml)=0.d0

            ml=nfcmax+nfc
            do mw=1,mwmax
               fma(mw,ml) = 0.d0
               if(mdlwmd.ge.1) then
                  do ns=0,nsmax
                     fma_save(mw,ml,nr,ns)=0.d0
                  enddo
               end if
            enddo
            fma(mwc,ml)=1.d0
            if(mdlwmd.ge.1) then
               fma_save(ml,ml,nr,0)=1.d0
            end if
            fvb(ml)=0.d0

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
      integer:: id_base=1

      do nfc=1,nfcmax
         ml=6*nfcmax*(nrmax-1)+nfc
         do mw=1,mwmax
            fma(mw,ml) = 0.d0
            if(mdlwmd.ge.1) then
               do ns=0,nsmax
                  fma_save(mw,ml,nrmax,ns)=0.d0
               enddo
            end if
         enddo
         fma(mwc,ml)=1.d0
         if(mdlwmd.ge.1) then
            fma_save(ml,ml,nrmax,0)=1.d0
         end if
         fvb(ml)=0.d0

         ml=6*nfcmax*(nrmax-1)+nfcmax*nfc
         do mw=1,mwmax
            fma(mw,ml) = 0.d0
            if(mdlwmd.ge.1) then
               do ns=0,nsmax
                  fma_save(mw,ml,nrmax,ns)=0.d0
               enddo
            end if
         enddo
         fma(mwc,ml)=1.d0
         if(mdlwmd.ge.1) then
            fma_save(ml,ml,nrmax,0)=1.d0
         end if
         fvb(ml)=0.d0

      enddo

      return
      end subroutine wmfem_boundary_condition_wall

!     -----     boundary condition for fourier components -----

      subroutine wmfem_boundary_condition_fourier

      integer:: mc,nr,nn,mm,mll,ml,mw,ns,nfc

      mc=(mwmax+1)/2

      if(nphmax.gt.1) then
         do nr=1,nrmax
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
                        fma_save(mw,mll+1,nr,ns)=0.d0
                        fma_save(mw,mll+3,nr,ns)=0.d0
                        fma_save(mw,mll+5,nr,ns)=0.d0
                        if(nr.ne.1) then
                           fma_save(mw,mll+6*nfcmax+1,nr-1,ns)=0.d0
                           fma_save(mw,mll+6*nfcmax+3,nr-1,ns)=0.d0
                           fma_save(mw,mll+6*nfcmax+5,nr-1,ns)=0.d0
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
                     fma_save(mc,mll+1,nr,ns)=1.d0
                     fma_save(mc,mll+3,nr,ns)=1.d0
                     fma_save(mc,mll+5,nr,ns)=1.d0
                     if(nr.ne.1) then
                        fma_save(mc,mll+6*nfcmax+1,nr-1,ns)=0.d0
                        fma_save(mc,mll+6*nfcmax+3,nr-1,ns)=0.d0
                        fma_save(mc,mll+6*nfcmax+5,nr-1,ns)=0.d0
                     endif
                  enddo
               end if
            end do
         end do
      endif

      if(nthmax.gt.1) then
         mm=nthmax/2+1
         do nr=1,nrmax
            do nn=1,nphmax
               mll=6*nthmax*(nn-1)+6*(mm-1)
               ml=6*nthmax*nphmax*(nr-1)+mll
               do mw=1,mwmax
                  fma(mw,ml+1)=0.d0
                  fma(mw,ml+3)=0.d0
                  fma(mw,ml+5)=0.d0
                  if(mdlwmd.ge.1) then
                     do ns=0,nsmax
                        fma_save(mw,mll+1,nr,ns)=0.d0
                        fma_save(mw,mll+3,nr,ns)=0.d0
                        fma_save(mw,mll+5,nr,ns)=0.d0
                        if(nr.ne.1) then
                           fma_save(mw,mll+6*nfcmax+1,nr-1,ns)=0.d0
                           fma_save(mw,mll+6*nfcmax+3,nr-1,ns)=0.d0
                           fma_save(mw,mll+6*nfcmax+5,nr-1,ns)=0.d0
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
                     fma_save(mc,mll+1,nr,ns)=1.d0
                     fma_save(mc,mll+3,nr,ns)=1.d0
                     fma_save(mc,mll+5,nr,ns)=1.d0
                     if(nr.ne.1) then
                        fma_save(mc,mll+6*nfcmax+1,nr-1,ns)=0.d0
                        fma_save(mc,mll+6*nfcmax+3,nr-1,ns)=0.d0
                        fma_save(mc,mll+6*nfcmax+5,nr-1,ns)=0.d0
                     endif
                  end do
               endif
            end do
         enddo
      endif

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
         IF(abs(mm) == 1) THEN
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

!     ***** wmfem_main end *********************************************

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
      complex(8),dimension(3,3,3,nfcmax2):: cqq
      complex(8),dimension(3,3,nfcmax2):: cpp
      complex(8),dimension(nthmax2,nphmax2):: fv1,fv1f
      integer:: i,j,k,l,nthm,nthp,nphm,nphp
      integer:: nfc2,nph,nth,imn,ml
      integer:: nfc1,nph1,nth1
      integer:: nn1,mm1,nn2,mm2,mmdiff,nndiff,nfcdiff
      real(8):: rho,dph,dth,gj
      complex(8):: cfactor
      
      cfactor=vc/(2.d0*pi*crf*1.d6)

      rho=rhoa(nr)

      call wmfem_tensors(rho,gma,muma,dmuma,gja)

!     ----- calculation rot coefficients cpp and cqq -----

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
                  cqq(i,j,imn,nfc2)=0.d0
                  do l=1,3
                     cqq(i,j,imn,nfc2)=cqq(i,j,imn,nfc2)
     &                 +gma(i,l,nth,nph)*cq(l,j,imn)*gj
                  enddo
               enddo
            enddo
         enddo
         do i=1,3
            do j=1,3
               cpp(i,j,nfc2)=0.d0
               do l=1,3
                  cpp(i,j,nfc2)=cpp(i,j,nfc2)
     &                 +gma(i,l,nth,nph)*cp(l,j)*gj
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
                  nph=nphnfc2(nfc2)
                  fv1(nth,nph)=cqq(i,j,imn,nfc2)
               enddo
               call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
               do nfc2=1,nfcmax2
                  nth=nthnfc2(nfc2)
                  nph=nphnfc2(nfc2)
                  cqq(i,j,imn,nfc2)=fv1f(nth,nph)
               enddo
            enddo
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nphnfc2(nfc2)
               fv1(nth,nph)=cpp(i,j,nfc2)
            enddo
            call wmsubfx(fv1,fv1f,nthmax2,nphmax2)
            do nfc2=1,nfcmax2
               nth=nthnfc2(nfc2)
               nph=nphnfc2(nfc2)
               cpp(i,j,nfc2)=fv1f(nth,nph)
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
            ml=6*nthmax*nphmax*(nr-1)+(nfc2-1)

!          !!!!! on axis should be considered !!!!!

            do i=1,3
               cbf(i,nth1,nph1)=0.d0
               do j=1,3
                  cbf(i,nth1,nph1)=cbf(i,nth1,nph1)
     &                 +cfactor*(
     &                  (cqq(i,j,1,nfcdiff)
     &                  +cqq(i,j,2,nfcdiff)*mm2
     &                  +cqq(i,j,3,nfcdiff)*nn2)*cef(j,nth1,nph1,nr)
     &                  +cpp(i,j,  nfcdiff)     *cdef(j,nth1,nph1,nr))
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

      allocate(fma_local(mbmax,mbmax))

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
               call wmfem_calculate_local(nr0,ns,fma_local)
            else
               do mb2=1,mbmax
                  do mb1=1,mbmax
                     fma_local(mb1,mb2)=fma_save(mb1,mb2,nr0,ns)
                  enddo
               enddo
            endif


            do nfc=1,nfcmax
               nn=nphnfc(nfc)
               mm=nthnfc(nfc)
               do i=1,8
                  mb=nfcmax*(i-1)+nfc
                  ml=6*nfcmax*(nr0-1)+mb
                  do nfc1=1,nfcmax
                     nn1=nphnfc(nfc1)
                     mm1=nthnfc(nfc1)
                     do i1=1,8
                        mb1=nfcmax*(i1-1)+nfc
                        ml1=6*nfcmax*(nr0-1)+mb1
                        do nfc2=1,nfcmax
                           nn2=nphnfc(nfc2)
                           mm2=nthnfc(nfc2)
                           nndiff=nnnfc(nfc2)-nnnfc(nfc1)
                           mmdiff=mmnfc(nfc2)-mmnfc(nfc1)
                           if(nndiff.lt.0) nndiff=nndiff+nphmax2
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
            do nn1=1,nphmax
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

      do nn=1,nphmax
      do mm=1,nthmax
         cpa(mm,nn)=0.d0
         do nr=1,nrmax-1
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

      deallocate(fma_local)

      return
      end subroutine wmfem_calculate_power

      end subroutine wmfem_post

!     ***** wmfem_post end *********************************************
