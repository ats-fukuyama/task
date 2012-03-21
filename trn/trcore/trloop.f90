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

       ! Interaction with EQ
!       IF(modelg .eq. 3) THEN
!          CALL TASK/EQ
!          CALL tr_bpsd_get(ierr)
!       END IF       

       CALL tr_step(ierr); IF(ierr /= 0) GO TO 9000

       IF(MOD(nt,ntstep) == 0 .OR. &
          MOD(nt,ngtstp) == 0) CALL tr_calc_global
       IF(MOD(nt,ntstep) == 0) CALL tr_status
       IF(MOD(nt,ngtstp) == 0) CALL tr_save_ngt

       CALL tr_bpsd_set(ierr)
       
    END DO

9000 CONTINUE

    CALL GUTIME(t2)

    write(6,'(A,F12.3,A)') '## cputime = ',T2 - T1,' (s)'

  END SUBROUTINE tr_loop

! **************************************************************************

  SUBROUTINE tr_save_pvprev
! ----------------------------------------------------
!      save plasma variables to previous values
! ----------------------------------------------------
    USE trcomm,ONLY: ikind,rkind,nrmax,nsamax, &
         dpdrho,rn,ru,rt,dpdrho_prev,rn_prev,ru_prev,rt_prev
    IMPLICIT NONE

    dpdrho_prev(0:nrmax)=dpdrho(0:nrmax)
    rn_prev(1:nsamax,0:nrmax)=rn(1:nsamax,0:nrmax)
    ru_prev(1:nsamax,0:nrmax)=ru(1:nsamax,0:nrmax)
    rt_prev(1:nsamax,0:nrmax)=rt(1:nsamax,0:nrmax)

    RETURN
  END SUBROUTINE tr_save_pvprev

! **************************************************************************

!!$  SUBROUTINE tr_eq_calc
!!$
!!$    DO
!!$    ! iteration for convergence of q(safety factor)
!!$       CALL eq_prof
!!$       CALL eq_calc
!!$       CALL tr_bpsd_get(ierr)
!!$
!!$       !convergence check of q
!!$    END DO
!!$
!!$    CALL tr_prof_q2dpdrho
!!$
!!$  END SUBROUTINE tr_eq_calc
!!$
!!$  SUBROUTINE tr_prof_q2dpdrho
!!$    USE trcomm, ONLY: rmu0,pi,nrmax,BB,RR,q0,qa,profj1,profj2,knameq, &
!!$         &         rhog,abb1rho,dvrho,ttrho,abrho,arrho,dpdrho,ar1rho,&
!!$         &         jtot,joh,bp,qp, nrd1,nrd2,nrd3,nrd4
!!$
!!$    IMPLICIT NONE
!!$    REAL(rkind) :: dr,prof,factor0,sumfact1
!!$    REAL(rkind) :: FCTR ! function defined in TASK/lib
!!$    INTEGER(ikind) :: nr
!!$
!!$! q --> dpdrho, j
!!$
!!$    dpdrho(
!!$
!!$  END SUBROUTINE tr_prof_q2dpdrho
  

END MODULE trloop
