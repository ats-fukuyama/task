! equtrunit.f90
  

!=======================================================================
!            interface module of "TOPICS-TR"
!                                                         06/08/21
!=======================================================================
      module trunit_mod
      public tr_init, tr_prof, tr_exec
      private
      contains
!=======================================================================
!            initialize tr module
!-----------------------------------------------------------------------
      subroutine tr_init
!
      use trset_mod
      use trpl_mod
      implicit none
!-----------------------------------------------------------------------
!      call trset
      call trpl_init
      return
      end subroutine tr_init
!=======================================================================
!            setup  profile
!-----------------------------------------------------------------------
      subroutine tr_prof
!
      use trset_mod
      use trpl_mod
      implicit none
      integer ierr
!-----------------------------------------------------------------------
      call trset
      call trpl_init
      call trpl_put(ierr)
      return
      end subroutine tr_prof
!=======================================================================
!            calculate transport
!-----------------------------------------------------------------------
      subroutine tr_exec
!
      use tradv_mod
      use trpl_mod
      implicit none
      integer ierr
!-----------------------------------------------------------------------
      call trpl_get(ierr)
      call tradv
      call trpl_put(ierr)
      return
      end subroutine tr_exec
!
      end module trunit_mod
