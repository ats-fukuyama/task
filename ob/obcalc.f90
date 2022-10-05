!   obcalc.f90

MODULE obcalc

  PRIVATE
  PUBLIC ob_calc

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE ob_calc(ierr)

    USE obcomm
    USE obprep
    USE obexec
    USE oblocal
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL(4):: time1,time2
    REAL(rkind):: y_in(neq_max)
    INTEGER:: nobt
    REAL(rkind):: &
         zetab_pos,thetab_pos,psip_pos,rhopara_pos,ptheta_pos,pzeta_pos, &
         bb_pos,phi_pos,penergy_pos,pcangle_pos,psit_pos,rr_pos,zz_pos, &
         qps_pos
    REAL(rkind):: dummy,fg_pos,fI_pos,pze,v_para,omega_c
    REAL(rkind):: omega_bounce1,omega_bounce2

    ierr=0

    CALL GUTIME(time1)

    ! --- line input of initial parameters ---
    
    SELECT CASE(mdlobi)
    CASE(100)
       WRITE(6,'(A,I5)') '## Initial parameters: nobt_max=',nobt_max
       WRITE(6,'(A)') &
            'nobt  penergy     pcangle     zeta        psipn       theta'
       DO nobt=1,nobt_max
100       CONTINUE
          WRITE(6,'(I4,1P5E12.4)') &
               nobt,penergy_ob_in(nobt),pcangle_ob_in(nobt),zeta_ob_in(nobt), &
               psipn_ob_in(nobt),theta_ob_in(nobt)
          READ(5,*,ERR=100,END=9000) &
               penergy_ob_in(nobt),pcangle_ob_in(nobt),zeta_ob_in(nobt), &
               psipn_ob_in(nobt),theta_ob_in(nobt)
       END DO
    CASE(101)
       WRITE(6,'(A,I5)') '## Initial parameters: nobt_max=',nobt_max
       WRITE(6,'(A)') &
            'nobt  penergy     pcangle     zeta        rr          zz'
       DO nobt=1,nobt_max
200       CONTINUE
          WRITE(6,'(I4,1P5E12.4)') &
               nobt,penergy_ob_in(nobt),pcangle_ob_in(nobt),zeta_ob_in(nobt), &
               rr_ob_in(nobt),zz_ob_in(nobt)
          READ(5,*,ERR=200,END=9000) &
               penergy_ob_in(nobt),pcangle_ob_in(nobt),zeta_ob_in(nobt), &
               rr_ob_in(nobt),zz_ob_in(nobt)
       END DO
    END SELECT

    ! calculate actual initial values ---
    
    SELECT CASE(mdlobi)
    CASE(0,100)
       DO nobt=1,nobt_max
          penergy_pos=penergy_ob_in(nobt)
          pcangle_pos=pcangle_ob_in(nobt)
          zetab_pos=zeta_ob_in(nobt)
          psip_pos=psipa*psipn_ob_in(nobt)
          thetab_pos=theta_ob_in(nobt)

    ! --- calculate constant quantities: peng, pmu, pcangl

          CALL cal_rr_pos(thetab_pos,psip_pos,rr_pos,ierr)
          CALL cal_bb_pos(thetab_pos,psip_pos,bb_pos,ierr)
          CALL cal_qps_pos(psip_pos,qps_pos,dummy,ierr)
          phi_pos=0.d0
          peng=AEE*1.D3*penergy_ob_in(nobt)
          pcangl=pcangle_ob_in(nobt)
          pze=PZ(ns_ob)*AEE
          pmu=(peng-pze*phi_pos)*(1.D0-pcangl**2)/bb_pos
          omega_bounce1=SQRT(pmu*bb_pos/(PA(ns_ob)*AMP) &
               *RA*SQRT(psipn_ob_in(nobt))/(qps_pos**2*rr_pos**3))
          omega_bounce2=pcangl*SQRT(2.D0*peng/(PA(ns_ob)*AMP))/(qps_pos*rr_pos)
          ! WRITE(6,'(A,2ES12.4)') 'omega_bounce=',omega_bounce1,omega_bounce2
          omega_bounce=MAX(omega_bounce1,ABS(omega_bounce2))

    ! --- calculate initial rhopara

          v_para=SQRT(2.D0*(peng-pze*phi_pos)/(PA(ns_ob)*AMP))*pcangl
          omega_c=pze*bb_pos/(PA(ns_ob)*AMP)
          rhopara_pos=v_para/omega_c

          y_in(1)= zetab_pos
          y_in(2)= thetab_pos
          y_in(3)= psip_pos
          y_in(4)= rhopara_pos

          CALL ob_exec(y_in,nobt,ierr)
          IF(ierr.NE.0) cycle
       END DO
    CASE(1,101)
       WRITE(6,'(A,I5)') 'XX obcalc: not-yet supported mdlobi:',mdlobi
    END SELECT
       
    CALL GUTIME(time2)
!     WRITE(6,*) '% CPU TIME = ',time2-time1,' sec'

9000 CONTINUE
    RETURN
  END SUBROUTINE ob_calc


END MODULE obcalc
