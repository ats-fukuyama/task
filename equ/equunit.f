c=======================================================================
c            interface program of "TOPICS-EQ"
c                                                         06/08/21
c=======================================================================
      module equunit_mod
      use eqinit_mod
      use eqsub_mod
      use eqfct_mod
      use eqset_mod
      use eqinp_mod
      use equpl_mod
      use eqgout_mod
      public equ_init, equ_parm, equ_view, equ_prof, equ_calc, 
     &       equ_gout, equ_save
      private
      contains
c=======================================================================
c            initialize eq module
c-----------------------------------------------------------------------
      subroutine equ_init
c
      use aaa_mod
      implicit none
c-----------------------------------------------------------------------
      ft05=5
      ft06=6
      out=0
c-----------------------------------------------------------------------
      call flxtab
      call eqinit
      return
      end subroutine equ_init
c=======================================================================
c            set parameters
c-----------------------------------------------------------------------
      subroutine equ_parm(mode,kin,ierr)
c
      implicit none
      character kin*(*)
      integer mode,ierr
      call eqparm(mode,kin,ierr)
      return
      end subroutine equ_parm
c=======================================================================
c            view parameters
c-----------------------------------------------------------------------
      subroutine equ_view
c
      implicit none
      call eqview
      return
      end subroutine equ_view
c=======================================================================
c            setup  profile
c-----------------------------------------------------------------------
      subroutine equ_prof
c
      implicit none
      integer ierr
c-----------------------------------------------------------------------
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
c=======================================================================
c            calculate equilibrium
c-----------------------------------------------------------------------
      subroutine equ_calc
c
      implicit none
      integer ierr
c-----------------------------------------------------------------------
c      call intequ
      call eqpl_get(ierr)
      call eqfct
      call eqout
c      call inttrn
      call eqpl_set(ierr)
      return
      end subroutine equ_calc
c=======================================================================
c            graphic output
c-----------------------------------------------------------------------
      subroutine equ_gout
c
      implicit none
      integer ierr
c-----------------------------------------------------------------------
      call eqgout
      return
      end subroutine equ_gout
c=======================================================================
c            file output
c-----------------------------------------------------------------------
      subroutine equ_save
c
      implicit none
      integer ierr
      character(len=80)::  knameq
c-----------------------------------------------------------------------
      knameq='eqdata'
      call equsave(knameq,5,ierr)
      return
      end subroutine equ_save
c=======================================================================
      end module equunit_mod
c=======================================================================
