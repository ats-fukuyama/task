c=======================================================================
c            interface program of "TOPICS-EQ"
c                                                         06/08/21
c=======================================================================
      module equunit_mod
      public equ_init, equ_parm, equ_view, equ_prof, equ_calc, 
     &       equ_gout, equ_save
      private
      contains
c=======================================================================
c            initialize eq module
c-----------------------------------------------------------------------
      subroutine equ_init
c
      implicit none
c-----------------------------------------------------------------------
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
      ierr=0
      return
      end subroutine equ_parm
c=======================================================================
c            view parameters
c-----------------------------------------------------------------------
      subroutine equ_view
c
      implicit none
      return
      end subroutine equ_view
c=======================================================================
c            setup  profile
c-----------------------------------------------------------------------
      subroutine equ_prof
c
      implicit none
      return
      end subroutine equ_prof
c=======================================================================
c            calculate equilibrium
c-----------------------------------------------------------------------
      subroutine equ_calc
c
      implicit none
      return
      end subroutine equ_calc
c=======================================================================
c            graphic output
c-----------------------------------------------------------------------
      subroutine equ_gout
c
      implicit none
      return
      end subroutine equ_gout
c=======================================================================
c            file output
c-----------------------------------------------------------------------
      subroutine equ_save
c
      implicit none
      return
      end subroutine equ_save
c=======================================================================
      end module equunit_mod
c=======================================================================
