!=======================================================================
!            interface program of "TASK-EQ"
!                                                         07/05/15
!=======================================================================
      module equnit_mod
      use eq_bpsd_mod
      public eq_init,eq_parm,eq_prof,eq_calc,eq_load,eq_gout
      private
      contains
!=======================================================================
!            initialize eq module
!-----------------------------------------------------------------------
      subroutine eq_init

      implicit none
      call eqinit
      call eqtabl
      return
      end subroutine eq_init
!=======================================================================
!            set parameters
!-----------------------------------------------------------------------
      subroutine eq_parm(mode,kin,ierr)

      implicit none
      character kin*(*)
      integer mode,ierr
      call eqparm(mode,kin,ierr)
      return
      end subroutine eq_parm
!=======================================================================
!            setup  profile
!-----------------------------------------------------------------------
      subroutine eq_prof

      implicit none
!-----------------------------------------------------------------------
      return
      end subroutine eq_prof
!=======================================================================
!            calculate equilibrium
!-----------------------------------------------------------------------
      subroutine eq_calc

      implicit none
      integer ierr
!-----------------------------------------------------------------------
      call eq_bpsd_get(ierr)
      call eqcalc(ierr)
!      call eq_bpsd_init(ierr)
      call eqcalq(ierr)
      call eq_bpsd_set(ierr)
      return
      end subroutine eq_calc
!=======================================================================
!            load equilibrium file
!-----------------------------------------------------------------------
      subroutine eq_load(modelg1,knameq1,ierr)

      include '../eq/eqcomm.inc'
      integer,intent(in):: modelg1
      character(len=80),intent(in):: knameq1
      integer,intent(out):: ierr
      integer:: nrmax_save,nthmax_save,nsumax_save
!-----------------------------------------------------------------------
      nrmax_save=nrmax
      nthmax_save=nthmax
      nsumax_save=nsumax
      call eqload(modelg1,knameq1,ierr)
      nrmax=nrmax_save
      nthmax=nthmax_save
      nsumax=nsumax_save

      call eq_bpsd_init(ierr)
      if(ierr.ne.0) write(6,*) 'XX eq_bpsd_init: ierr=',ierr
      call eqcalq(ierr)
      if(ierr.ne.0) write(6,*) 'XX eq_calq: ierr=',ierr
      call eq_bpsd_set(ierr)
      if(ierr.ne.0) write(6,*) 'XX eq_bpsd_set: ierr=',ierr
      return
      end subroutine eq_load
!=======================================================================
!            graphic output
!-----------------------------------------------------------------------
      subroutine eq_gout

      implicit none
!-----------------------------------------------------------------------
      call eqgout(0)
      return
      end subroutine eq_gout
!=======================================================================
      end module equnit_mod
!=======================================================================
