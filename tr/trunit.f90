!     $Id$
!=======================================================================
!            interface of "task/eq"
!=======================================================================
      module trunit
      use tr_bpsd
      public tr_init,tr_parm,tr_view,tr_prof,tr_exec,tr_load,tr_save, &
      &      tr_gout,tr_fout,tr_term
      private

     contains
!=======================================================================
!            initialize tr component
!-----------------------------------------------------------------------
      subroutine tr_init

      implicit none
      open(7,STATUS='SCRATCH',FORM='FORMATTED')
      call trinit
      return
      end subroutine tr_init
!=======================================================================
!            set parameters
!-----------------------------------------------------------------------
      subroutine tr_parm(mode,kin,ierr)

      implicit none
      character kin*(*)
      integer mode,ierr
      call trparm(mode,kin,ierr)
      return
      end subroutine tr_parm
!=======================================================================
!            view parameters
!-----------------------------------------------------------------------
      subroutine tr_view

      implicit none
      call trview(0)
      return
      end subroutine tr_view
!=======================================================================
!            setup profile
!-----------------------------------------------------------------------
      subroutine tr_prof(ierr)

      implicit none
      integer,intent(out):: ierr
      call trprof             ! initialize profile data
      call tr_bpsd_init       ! initialize tr_bpsd
      call tr_bpsd_set(ierr)  ! set tr_bpsd with initial profile
      return
      end subroutine tr_prof
!=======================================================================
!            calculate transport one time step
!-----------------------------------------------------------------------
      subroutine tr_exec(dt,ierr)

      implicit none
      real(8),intent(in):: dt
      integer,intent(out):: ierr

      call tr_bpsd_get(ierr)
      if(ierr.ne.0) return
      call trexec(dt,ierr)
      if(ierr.ne.0) return
      call tr_bpsd_set(ierr)
      return
      end subroutine tr_exec
!=======================================================================
!            load file
!-----------------------------------------------------------------------
      subroutine tr_load(ierr)

      implicit none
      integer,intent(out):: ierr
      call trload
      call tr_bpsd_init
      call tr_bpsd_set(ierr)  ! set tr_bpsd with initial profile
      return
      end subroutine tr_load
!=======================================================================
!            save file
!-----------------------------------------------------------------------
      subroutine tr_save(ierr)

      implicit none
      integer,intent(out):: ierr
      call trsave(ierr)
      return
      end subroutine tr_save
!=======================================================================
!            graphics
!-----------------------------------------------------------------------
      subroutine tr_gout

      implicit none
      call trgout
      return
      end subroutine tr_gout
!=======================================================================
!            text file
!-----------------------------------------------------------------------
      subroutine tr_fout

      implicit none
      call trfout
      return
      end subroutine tr_fout
!=======================================================================
!            terminate tr component
!-----------------------------------------------------------------------
      subroutine tr_term

      implicit none
      close(7)
      return
      end subroutine tr_term
!=======================================================================
      end module trunit
!=======================================================================
