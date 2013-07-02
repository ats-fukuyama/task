!     $Id$

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
         CALL get_wmfem_size(nrmax,nthmax,nhhmax,nsmax)

!     ***** define array size  *****

      nfcmax=nthmax*nhhmax      ! size of block matrix 
                                !    (number of Fourier components)
      mlmax=nfcmax*(6*nrmax-4)  ! length of coeffient matrix and source vector
                                !   E_perp(i)
                                !   E_para(i)
                                !   E_rho (i+1/4)
                                !   E_perp(i+1/2)
                                !   E_para(i+1/2)
                                !   E_rho (i+3/4)

      mbmax=nfcmax*8            ! size of block matrix
                                !   E_perp(i)
                                !   E_para(i)
                                !   E_rho (i+1/4)
                                !   E_perp(i+1/2)
                                !   E_para(i+1/2)
                                !   E_rho (i+3/4)
                                !   E_perp(i+1)
                                !   E_para(i+1)
      mwmax=2*mbmax-1           ! width of coefficient matrix
      mwc=mbmax                 ! position of diagonal coponent

      if(nthmax.eq.1) then
         nthmax2=1
      else
         nthmax2=nthmax*2
      endif
      if(nhhmax.eq.1) then
         nhhmax2=1
      else
         nhhmax2=nhhmax*2
      endif
      nfcmax2=nthmax2*nhhmax2

!     ***** get additional parameters *****

      call get_wmfem_parm(crf,nth0,nph0,mdlwmf,mdlwmd)

      return
      end subroutine wmfem_pre
