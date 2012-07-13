!     $Id$

      module wmfem_comm

      complex(8),parameter:: ci=(0.d0,1.d0)
      real(8),parameter:: pi = 3.14159265358979D0
      real(8),parameter:: vc = 2.99792458 D8

      complex(8):: crf
      integer:: nth0,nph0
      integer:: mdlwmf,mdlwmd

      integer:: nrmax,nthmax,nphmax,nsmax
      real(8),dimension(:),ALLOCATABLE:: rhoa
      integer:: nfcmax,nthmax2,nphmax2,nfcmax2
      integer:: mlmax,mwmax,mbmax,mwc
      complex(8),dimension(:,:,:,:),ALLOCATABLE:: cef,cbf 
                                !(3,nthmax,nphmax,nrmax)
      complex(8),dimension(:,:,:,:,:,:),ALLOCATABLE:: cpp 
                                !(nthmax,nphmax,nthmax2,nphmax2,nrmax,0:nsmax)
      complex(8),dimension(:,:),ALLOCATABLE:: cpa
                                !(nthmax,nphmax)

      complex(8),dimension(:,:),ALLOCATABLE:: fma !(mwmax,mlmax)
      complex(8),dimension(:,:,:,:),ALLOCATABLE:: fma_save 
     &                          !(mbmax,mbmax,nrmax,0:nsmax)
      complex(8),dimension(:),ALLOCATABLE:: fvb,fvx !(mlmax)

      integer,dimension(:),ALLOCATABLE :: nthnfc,mmnfc
      integer,dimension(:),ALLOCATABLE :: nphnfc,nnnfc
      integer,dimension(:),ALLOCATABLE :: nthnfc2,mmnfc2
      integer,dimension(:),ALLOCATABLE :: nphnfc2,nnnfc2

      end module wmfem_comm
