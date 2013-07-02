!     $Id$

!     ----- allocate arrays -----

      subroutine wmfem_allocate

      use wmfem_comm
      implicit none
      integer,save:: nrmax_save=0,nthmax_save=0,nhhmax_save=0
      integer,save:: mwmax_save=0,mlmax_save=0,mbmax_save=0
      integer,save:: nsmax_save=0,nfcmax_save=0
      integer,save:: mdlwmd_save=0 

      if((nrmax.ne.nrmax_save).or.(nthmax.ne.nthmax_save).or. 
     &   (nhhmax.ne.nhhmax_save)) then
         if(ALLOCATED(cef)) deallocate(cef)
         allocate(cef(3,nthmax,nhhmax,nrmax))
         if(ALLOCATED(cdef)) deallocate(cdef)
         allocate(cdef(3,nthmax,nhhmax,nrmax))
         if(ALLOCATED(cbf)) deallocate(cbf)
         allocate(cbf(3,nthmax,nhhmax,nrmax))
         if(ALLOCATED(cpp)) deallocate(cpp)
         allocate(cpp(nthmax,nhhmax,nthmax2,nhhmax2,nrmax,0:nsmax))
         if(ALLOCATED(cpa)) deallocate(cpa)
         allocate(cpa(nthmax,nhhmax))
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
         if(ALLOCATED(nhhnfc)) deallocate(nhhnfc)
         if(ALLOCATED(mmnfc)) deallocate(mmnfc)
         if(ALLOCATED(nnnfc)) deallocate(nnnfc)
         if(nfcmax.ne.0) allocate(nthnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nhhnfc(nfcmax))
         if(nfcmax.ne.0) allocate(mmnfc(nfcmax))
         if(nfcmax.ne.0) allocate(nnnfc(nfcmax))

         if(ALLOCATED(nthnfc2)) deallocate(nthnfc2)
         if(ALLOCATED(nhhnfc2)) deallocate(nhhnfc2)
         if(ALLOCATED(mmnfc2)) deallocate(mmnfc2)
         if(ALLOCATED(nnnfc2)) deallocate(nnnfc2)
         if(nfcmax2.ne.0) allocate(nthnfc2(nfcmax2))
         if(nfcmax2.ne.0) allocate(nhhnfc2(nfcmax2))
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
