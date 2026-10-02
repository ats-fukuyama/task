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
      mlmax=nfcmax*(8*nrmax)    ! length of coeffient matrix and source vector
                                !   A_perp(i)
                                !   A_para(i)
                                !   A_rho (i)
                                !   PHI   (i)
                                !   DA_perp/DRHO(i)
                                !   DA_para/DRHO(i)
                                !   DA_rho/DRHO (i)
                                !   DPHI/DRHO   (i)

      mbmax=nfcmax*16           ! size of block matrix
                                !   A_perp(i)
                                !   A_para(i)
                                !   A_rho (i)
                                !   PHI   (i)
                                !   DA_perp/DRHO(i)
                                !   DA_para/DRHO(i)
                                !   DA_rho/DRHO (i)
                                !   DPHI/DRHO   (i)
                                !   A_perp(i+1)
                                !   A_para(i+1)
                                !   A_rho (i+1)
                                !   PHI   (i+1)
                                !   DA_perp/DRHO(i+1)
                                !   DA_para/DRHO(i+1)
                                !   DA_rho/DRHO (i+1)
                                !   DPHI/DRHO   (i+1)
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
