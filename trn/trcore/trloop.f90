MODULE trloop

! This module advances time upto ntmax*dt

  USE trcomm,ONLY: rkind, ikind, nrmax, neqmax, nsamax
  IMPLICIT NONE

  PUBLIC tr_loop,tr_save_pvprev
  PRIVATE

CONTAINS

  SUBROUTINE tr_loop

    USE trcomm, ONLY: rkev,ntmax,t,dt,ntstep,ngtstp,nwrstp, &
                      mdluf,tlmax,tmax,mdltr_prv,dtr_tb,dtr_tb_prev
    USE trbpsd, ONLY: tr_bpsd_set,tr_bpsd_get
    USE trcalc1, ONLY: tr_calc1
    USE trstep, ONLY: tr_step
    USE trresult, ONLY: tr_status,tr_latest_status, &
                        tr_calc_global,tr_save_ngt,tr_exp_compare
    USE trwrite, ONLY: tr_write_open,tr_write_close,tr_writep_csv,tr_writet_csv
    IMPLICIT NONE

    REAL(4)        :: t1,t2,dt_save
    INTEGER(ikind) :: nt,neq,ierr, fflag

    CALL tr_write_open(ierr)
    CALL tr_writet_csv

    IF(mdluf == 2 .AND. tmax > tlmax)THEN ! fit experimental data
       tmax = tlmax
       WRITE(6,*) ' ## fit the simulation time to the experimental data.'
       WRITE(6,'(A16,F7.2,A4)') ' ##  --> 0.00 - ',tmax,' [s]'
    ELSE IF(t /= 0.d0)THEN
       tmax = tmax + dt*ntmax
    END IF
    fflag = 0 ! flag for end of loop condition 

    CALL GUTIME(t1)

    dt_save = dt
    nt = 0
    time_evolution: DO ! main loop
       IF(ntmax==0) EXIT
       nt = nt + 1

       CALL tr_save_pvprev

       IF(mdltr_prv >= 10 .AND. nt >= 3)THEN
          CALL Pereverzev_step_control(dt)
       END IF
       DO neq = 1, neqmax
          dtr_tb_prev(neq,1:nrmax) = dtr_tb(neq,neq,1:nrmax)
       END DO

       IF(t + dt >= tmax)THEN
          dt = tmax - t
          fflag = 1
       ELSE IF(nt == ntmax - 1)THEN
          ntmax = ntmax + 1
       END IF
       ! incremental addtion
       t = t + dt

       CALL tr_calc1

       ! non-linear iteration
       CALL tr_step(ierr); IF(ierr /= 0) EXIT

       IF(MOD(nt,ntstep) == 0 .OR. &
          MOD(nt,ngtstp) == 0) THEN
                                CALL tr_calc_global
          IF(mdluf > 0)         CALL tr_exp_compare
       END IF

       IF(MOD(nt,ntstep) == 0)  CALL tr_status

       IF(MOD(nt,ngtstp) == 0)  CALL tr_save_ngt
       IF(MOD(nt,nwrstp) == 0)  CALL tr_writet_csv

       CALL tr_bpsd_set(ierr)
       IF(ierr /= 0)THEN
          WRITE(6,*) 'XX tr_loop: Error to set the variables to BPSD interface. NT= ', nt
       END IF

       IF(nt==ntmax .OR. fflag==1)THEN
          dt = dt_save
          EXIT 
       END IF
    END DO time_evolution

    CALL GUTIME(t2)
    WRITE(6,'(A,F12.3,A)') '## cputime = ',T2 - T1,' (s)'

    CALL tr_latest_status
    CALL tr_writep_csv
    CALL tr_write_close(ierr)

    RETURN
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

  SUBROUTINE Pereverzev_step_control(dt)
! -------------------------------------------------------------------------
!
! -------------------------------------------------------------------------
    USE trcomm, ONLY: t,tmax,ntmax,nrmax,neqmax,nsa_neq,nva_neq, &
         idnsa,dtr_tb,dtr_tb_prev,rhog,rhog_fixed,deldtl,cdtprv
    IMPLICIT NONE

    REAL(rkind),   INTENT(INOUT) :: dt

    INTEGER(ikind) :: nsa,nva,neq,nr
    REAL(rkind) :: dt_prv,dt_max,dt_min
    REAL(rkind),DIMENSION(1:neqmax,0:nrmax) :: deldtr

    dt_max = 1.d-2
    dt_min = 1.d-5

    deldtr(1:neqmax,0:nrmax) = 0.d0
    DO neq = 1, neqmax
       nva = nva_neq(neq)
       nsa = nsa_neq(neq)
       IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE

       DO nr = 1, nrmax
          IF(rhog(nr) > rhog_fixed(nva,nsa)) EXIT
          IF(nva==3)THEN
             deldtr(neq,nr) &
                   = ABS((dtr_tb(neq,neq,nr)-dtr_tb_prev(neq,nr)) &
                         /dtr_tb(neq,neq,nr))
          END IF
       END DO
    END DO

    IF(MAXVAL(deldtr(1:neqmax,1:nrmax)) > deldtl)THEN
       ! decrease the time step width
       dt_prv = dt / cdtprv
       IF(dt_prv < dt_min) dt_prv = dt_min
    ELSE
       ! increase the time step width
       dt_prv = cdtprv * dt
       IF(dt_prv > dt_max) dt_prv = dt_max
    END IF

    dt = dt_prv

    RETURN
  END SUBROUTINE Pereverzev_step_control

END MODULE trloop
