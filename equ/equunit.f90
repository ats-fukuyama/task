! equunit.f90

!=======================================================================
!            interface program of "TOPICS-EQ"
!                                                         06/08/21
!=======================================================================
      module equunit_mod
      use eqinit_mod
      use eqsub_mod
      use eqfct_mod
      use eqset_mod
      use eqinp_mod
      use equpl_mod
      use eqgout_mod
      public equ_init, equ_parm, equ_view, equ_prof, equ_calc, &
             equ_gout, equ_save
      private
      contains
!=======================================================================
!            initialize eq module
!-----------------------------------------------------------------------
      subroutine equ_init
!
      use aaa_mod
      implicit none
!-----------------------------------------------------------------------
      ft05=5
      ft06=6
      out=0
!-----------------------------------------------------------------------
      call flxtab
      call eqinit
      return
      end subroutine equ_init
!=======================================================================
!            set parameters
!-----------------------------------------------------------------------
      subroutine equ_parm(mode,kin,ierr)
!
      implicit none
      character kin*(*)
      integer mode,ierr
      call eqparm(mode,kin,ierr)
      return
      end subroutine equ_parm
!=======================================================================
!            view parameters
!-----------------------------------------------------------------------
      subroutine equ_view
!
      implicit none
      call eqview
      return
      end subroutine equ_view
!=======================================================================
!            setup  profile
!-----------------------------------------------------------------------
      subroutine equ_prof
!
      implicit none
      integer ierr
!-----------------------------------------------------------------------
      call eqflin(18,ierr)
      if(ierr.ne.0) return
      call eqchek(ierr)
      if(ierr.ne.0) return
      call eqgrd
      call eqset
      call eqout
      call eqpl_init(ierr)
      call eqpl_prof(ierr) ! adjust temperature and q profile
      return
      end subroutine equ_prof
!=======================================================================
!            calculate equilibrium
!-----------------------------------------------------------------------
      subroutine equ_calc
!
      implicit none
      integer ierr
!-----------------------------------------------------------------------
!      call intequ
      call eqpl_get(ierr)
      call eqfct
      call eqout
!      call inttrn
      call eqpl_put(ierr)
      return
      end subroutine equ_calc
!=======================================================================
!            graphic output
!-----------------------------------------------------------------------
      subroutine equ_gout
!
      implicit none
      integer ierr
!-----------------------------------------------------------------------
      call eqgout
      return
      end subroutine equ_gout
!=======================================================================
!            file output
!-----------------------------------------------------------------------
      subroutine equ_save
!
      implicit none
      integer ierr
      character(len=80)::  knameq
!-----------------------------------------------------------------------
      knameq='eqdata'
      call equsave(knameq,5,ierr)
      return
      end subroutine equ_save
!=======================================================================
      end module equunit_mod
!=======================================================================
