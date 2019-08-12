c=======================================================================
c            interface module of "TOPICS-TR"
c                                                         06/08/21
c=======================================================================
      module trunit_mod
      public tr_init, tr_prof, tr_exec
      private
      contains
c=======================================================================
c            initialize tr module
c-----------------------------------------------------------------------
      subroutine tr_init
c
      use trset_mod
      use trpl_mod
      implicit none
c-----------------------------------------------------------------------
c      call trset
      call trpl_init
      return
      end subroutine tr_init
c=======================================================================
c            setup  profile
c-----------------------------------------------------------------------
      subroutine tr_prof
c
      use trset_mod
      use trpl_mod
      implicit none
      integer ierr
c-----------------------------------------------------------------------
      call trset
      call trpl_init
      call trpl_put(ierr)
      return
      end subroutine tr_prof
c=======================================================================
c            calculate transport
c-----------------------------------------------------------------------
      subroutine tr_exec
c
      use tradv_mod
      use trpl_mod
      implicit none
      integer ierr
c-----------------------------------------------------------------------
      call trpl_get(ierr)
      call tradv
      call trpl_put(ierr)
      return
      end subroutine tr_exec
c
      end module trunit_mod
