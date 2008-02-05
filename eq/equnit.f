c=======================================================================
c            interface program of "TASK-EQ"
c                                                         07/05/15
c=======================================================================
      module equnit_mod
      use eq_bpsd_mod
      public eq_init,eq_parm,eq_prof,eq_calc,eq_load,eq_gout
      private
      contains
c=======================================================================
c            initialize eq module
c-----------------------------------------------------------------------
      subroutine eq_init
c
      implicit none
      call eqinit
      call eqtabl
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
      call eqparm(mode,kin,ierr)
      return
      end subroutine eq_parm
c=======================================================================
c            setup  profile
c-----------------------------------------------------------------------
      subroutine eq_prof
c
      implicit none
c-----------------------------------------------------------------------
      return
      end subroutine eq_prof
c=======================================================================
c            calculate equilibrium
c-----------------------------------------------------------------------
      subroutine eq_calc
c
      implicit none
      integer ierr
c-----------------------------------------------------------------------
      call eqcalc(ierr)
      call eq_bpsd_init(ierr)
      call eqcalq(ierr)
      call eq_bpsd_set(ierr)
      return
      end subroutine eq_calc
c=======================================================================
c            load equilibrium file
c-----------------------------------------------------------------------
      subroutine eq_load(modelg1,knameq1,ierr)
c
      include '../eq/eqcomm.inc'
      integer,intent(in):: modelg1
      character(len=80),intent(in):: knameq1
      integer,intent(out):: ierr
      integer:: nrmax_save,nthmax_save,nsumax_save
c-----------------------------------------------------------------------
      nrmax_save=nrmax
      nthmax_save=nthmax
      nsumax_save=nsumax
      call eqload(modelg1,knameq1,ierr)
      nrmax=nrmax_save
      nthmax=nthmax_save
      nsumax=nsumax_save

      call eq_bpsd_init(ierr)
      call eqcalq(ierr)
      call eq_bpsd_set(ierr)
      return
      end subroutine eq_load
c=======================================================================
c            graphic output
c-----------------------------------------------------------------------
      subroutine eq_gout
c
      implicit none
c-----------------------------------------------------------------------
      call eqgout(0)
      return
      end subroutine eq_gout
c=======================================================================
      end module equnit_mod
c=======================================================================
