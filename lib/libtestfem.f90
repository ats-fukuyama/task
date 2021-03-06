! libtestfem.f90

MODULE libtestfem

  USE task_kinds,ONLY: dp
  implicit none
  integer:: mwmax,mlmax
  complex(dp),dimension(:,:),allocatable:: fma
  complex(dp),dimension(:),allocatable:: fvb,fvx
  complex(dp),dimension(:,:),allocatable:: cf1,cf2,cf3
  real(dp),dimension(:),allocatable:: rho

CONTAINS

!----- set profile -----

      subroutine mesh_init(nrmax,npow)

      implicit none
      integer,intent(in):: nrmax  ! number of points including end points
      integer,intent(in):: npow   ! power of mesh position
      integer,save:: nrmax_save=0
      real(dp):: drho
      integer:: nr

      if(nrmax.ne.nrmax_save) then
         if(allocated(rho)) deallocate(rho)
         if(allocated(cf1)) deallocate(cf1)
         if(allocated(cf2)) deallocate(cf2)
         if(allocated(cf3)) deallocate(cf3)
         allocate(rho(nrmax))
         allocate(cf1(nrmax,3))
         allocate(cf2(nrmax,3))
         allocate(cf3(nrmax,3))
      endif
      nrmax_save=nrmax

      drho=1.d0/(nrmax-1)**npow
      do nr=1,nrmax
         rho(nr)=drho*(nr-1)**npow
      enddo
      return
      end subroutine mesh_init

!----- allocate matrix and vectors -----

      subroutine fem_init

      use libfem
      implicit none
      integer,save:: mwmax_save=0,mlmax_save=0
      integer:: mc

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif
    
!      mwmax=4*nvmax-1         ! width of coefficient matrix
!      mlmax=nvmax*nrmax      ! length of coeffient matrix and source vector

      if((mwmax.ne.mwmax_save).or.(mlmax.ne.mlmax_save)) then
         if(allocated(fma)) deallocate(fma)
         allocate(fma(mwmax,mlmax))
         mc=(mwmax+1)/2
         write(6,'(A,3I5)') 'mlmax,mwmax,mc=',mlmax,mwmax,mc
      endif

      if(mlmax.ne.mlmax_save) then
         if(allocated(fvb)) deallocate(fvb)
         if(allocated(fvx)) deallocate(fvx)
         allocate(fvb(mlmax))
         allocate(fvx(mlmax))
      endif
      mwmax_save=mwmax
      mlmax_save=mlmax
      end subroutine fem_init

END MODULE libtestfem

