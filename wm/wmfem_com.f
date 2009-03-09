!     $Id$

      module wmfem_com

      complex(8),parameter:: ci=(0.d0,1.d0)
      real(8),parameter:: pi = 3.14159265358979D0
      real(8),parameter:: vc = 2.99792458 D8

      complex(8):: crf
      integer:: nth0,nph0
      integer:: mdlwmf,mdlwmd

      integer:: nrmax,nthmax,nphmax,nsmax
      real(8),dimension(:),pointer:: rhoa
      integer:: nfcmax,mlmax,mwmax,nthmax2,nphmax2,nfcmax2
      complex(8),dimension(:,:,:,:),pointer:: cef,cbf 
                                !(3,nthmax,nphmax,nrmax)
      complex(8),dimension(:,:,:,:,:,:),pointer:: cpp 
                                !(nthmax,nphmax,nthmax2,nphmax2,nrmax,0:nsmax)
      complex(8),dimension(:,:),pointer:: cpa
                                !(nthmax,nphmax)

      complex(8),dimension(:,:),pointer:: fma !(mwmax,mlmax)
      complex(8),dimension(:,:,:,:),pointer:: fms 
     &                          !(mwmax,12*nfcmax,nrmax,0:nsmax)
      complex(8),dimension(:),pointer:: fvb,fvx !(mlmax)

      integer,dimension(:),pointer :: nthnfc,mmnfc
      integer,dimension(:),pointer :: nphnfc,nnnfc
      integer,dimension(:),pointer :: nthnfc2,mmnfc2
      integer,dimension(:),pointer :: nphnfc2,nnnfc2

      end module wmfem_com

