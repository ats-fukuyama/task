c=======================================================================
c            interface program of "EQ"
c                                                         07/04/10
c=======================================================================
      module equnit_mod
      public eq_init,eq_parm,eq_prof,eq_calc,eq_load,eq_gout
      private
      contains
c=======================================================================
c            initialize eq module
c-----------------------------------------------------------------------
      subroutine eq_init
      return
      end subroutine eq_init
c=======================================================================
c            set parameters  profile
c-----------------------------------------------------------------------
      subroutine eq_parm(mode,kin,ierr)
c
      implicit none
      character kin*(*)
      integer mode,ierr
      ierr=0
      return
      end subroutine eq_parm
c=======================================================================
c            setup  profile
c-----------------------------------------------------------------------
      subroutine eq_prof
      return
      end subroutine eq_prof
c=======================================================================
c            calculate equilibrium
c-----------------------------------------------------------------------
      subroutine eq_calc
      return
      end subroutine eq_calc
c=======================================================================
c            load equilibrium file
c-----------------------------------------------------------------------
      subroutine eq_load(modelg1,knameq1,ierr)
c
      integer,intent(in):: modelg1
      character(len=80),intent(in):: knameq1
      integer,intent(out):: ierr
      ierr=0
      return
      end subroutine eq_load
c=======================================================================
c            graphic output
c-----------------------------------------------------------------------
      subroutine eq_gout
      return
      end subroutine eq_gout
c=======================================================================
      end module equnit_mod
c=======================================================================
