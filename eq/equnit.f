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
      call eqcalq(ierr)
      call eq_bpsd_put(ierr)
      return
      end subroutine eq_calc
!=======================================================================
!            load equilibrium file
!-----------------------------------------------------------------------
      SUBROUTINE eq_load(modelg1,knameq1,ierr)

      INCLUDE '../eq/eqcomm.inc'
      INTEGER,INTENT(IN):: modelg1
      CHARACTER(LEN=80),INTENT(IN):: knameq1
      INTEGER,INTENT(OUT):: ierr
      INTEGER:: nrmax_save,nthmax_save,nsumax_save
!-----------------------------------------------------------------------
      nrmax_save=nrmax
      nthmax_save=nthmax
      nsumax_save=nsumax
      CALL eqload(modelg1,knameq1,ierr)
      nrmax=nrmax_save
      nthmax=nthmax_save
      nsumax=nsumax_save
      IF(ierr.NE.0) THEN
         WRITE(6,*) 'XX eq_load: eqload: ierr=',ierr
         RETURN
      ENDIF

      CALL eq_bpsd_init(ierr)
      IF(ierr.NE.0) THEN
         write(6,*) 'XX eq_load: eq_bpsd_init: ierr=',ierr
         RETURN
      ENDIF

      CALL eqcalq(ierr)
      IF(ierr.NE.0) THEN
         WRITE(6,*) 'XX eq_load: eqcalq: ierr=',ierr
         RETURN
      ENDIF

      call eq_bpsd_put(ierr)
      IF(ierr.NE.0) THEN
         WRITE(6,*) 'XX eq_load: eq_bpsd_put: ierr=',ierr
         RETURN
      ENDIF
      RETURN
      END SUBROUTINE eq_load
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
