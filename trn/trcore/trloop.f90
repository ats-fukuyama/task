MODULE trloop

! This module advances time upto ntmax*dt

  PUBLIC tr_loop,tr_save_pvprev
  PRIVATE

CONTAINS

  SUBROUTINE tr_loop

    USE trcomm, ONLY: ikind,rkind,ntmax,t,dt,ntstep,ngtstp
    USE trbpsd, ONLY: tr_bpsd_set,tr_bpsd_get
    USE trstep, ONLY: tr_step
    USE trresult, ONLY: tr_status,tr_calc_global,tr_save_ngt
    IMPLICIT NONE

    REAL(4):: t1,t2
    INTEGER(ikind):: nt,ierr

    CALL GUTIME(t1)

    DO nt = 1, ntmax ! main loop

       t=t+dt
       CALl tr_save_pvprev

!       CALL TASK/EQ
!       CALL tr_bpsd_get(ierr)

       CALL tr_step(ierr); IF(ierr /= 0) GO TO 9000

       CALL tr_bpsd_set(ierr)

       IF(MOD(nt,ntstep) == 0 .OR. &
          MOD(nt,ngtstp) == 0) CALL tr_calc_global
       IF(MOD(nt,ntstep) == 0) CALL tr_status
       IF(MOD(nt,ngtstp) == 0) CALL tr_save_ngt
       
    END DO

9000 CONTINUE

    CALL GUTIME(t2)

    write(6,'(A,F12.3,A)') '## cputime = ',T2 - T1,' (s)'

  END SUBROUTINE tr_loop

! ***** save plasma variables to previous values *****

  SUBROUTINE tr_save_pvprev

    USE trcomm,ONLY: ikind,rkind,nrmax,nsamax, &
         qp,rn,ru,rt,qp_prev,rn_prev,ru_prev,rt_prev
    IMPLICIT NONE

    qp_prev(0:nrmax)=qp(0:nrmax)
    rn_prev(1:nsamax,0:nrmax)=rn(1:nsamax,0:nrmax)
    ru_prev(1:nsamax,0:nrmax)=ru(1:nsamax,0:nrmax)
    rt_prev(1:nsamax,0:nrmax)=rt(1:nsamax,0:nrmax)

    RETURN
  END SUBROUTINE tr_save_pvprev

END MODULE trloop
