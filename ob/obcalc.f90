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
         bb_pos,phi_pos,penergy_pos,pcangle_pos,psit_pos,rr_pos,zz_pos
    REAL(rkind):: dummy,fg_pos,fI_pos,pze,v_para,omega_c

    ierr=0

    CALL GUTIME(time1)

    ! --- line input of initial parameters ---
    
    SELECT CASE(mdlobi)
    CASE(100)
       WRITE(6,'(A,I5)') '## Initial parameters: nobt_max=',nobt_max
       WRITE(6,'(A)') &
            'nobt  penergy_in  pcangle_in  zeta_in     psipn_in    theta_in'
       DO nobt=1,nobt_max
100       CONTINUE
          WRITE(6,'(I4,1P5E12.4)') &
               nobt,penergy_in(nobt),pcangle_in(nobt),zeta_in(nobt), &
               psipn_in(nobt),theta_in(nobt)
          READ(5,*,ERR=100,END=9000) &
               penergy_in(nobt),pcangle_in(nobt),zeta_in(nobt), &
               psipn_in(nobt),theta_in(nobt)
       END DO
    CASE(101)
       WRITE(6,'(A,I5)') '## Initial parameters: nobt_max=',nobt_max
       WRITE(6,'(A)') &
            'nobt  penergy_in  pcangle_in  zeta_in     rr_in       zz_in'
       DO nobt=1,nobt_max
200       CONTINUE
          WRITE(6,'(I4,1P5E12.4)') &
               nobt,penergy_in(nobt),pcangle_in(nobt),zeta_in(nobt), &
               rr_in(nobt),zz_in(nobt)
          READ(5,*,ERR=200,END=9000) &
               penergy_in(nobt),pcangle_in(nobt),zeta_in(nobt), &
               rr_in(nobt),zz_in(nobt)
       END DO
    END SELECT

    ! calculate actual initial values ---
    
    SELECT CASE(mdlobi)
    CASE(0,100)
       DO nobt=1,nobt_max
          penergy_pos=penergy_in(nobt)
          pcangle_pos=pcangle_in(nobt)
          zetab_pos=zeta_in(nobt)
          psip_pos=psipa*psipn_in(nobt)
          thetab_pos=theta_in(nobt)

    ! --- calculate constant quantities: peng, pmu, pcangl

          CALL cal_bb_pos(thetab_pos,psip_pos,bb_pos,ierr)
          phi_pos=0.d0
          peng=AEE*1.D3*penergy_in(nobt)
          pcangl=pcangle_in(nobt)
          pze=PZ(ns_ob)*AEE
          pmu=(peng-pze*phi_pos)*(1.D0-pcangl**2)/bb_pos

    ! --- calculate initial rhopara

          v_para=SQRT(2.D0*(peng-pze*phi_pos)/(PA(ns_ob)*AMP))*pcangl
          omega_c=pze*bb_pos/(PA(ns_ob)*AMP)
          rhopara_pos=v_para/omega_c

          y_in(1)= zetab_pos
          y_in(2)= thetab_pos
          y_in(3)= psip_pos
          y_in(4)= rhopara_pos
       END DO
    CASE(1,101)
       WRITE(6,'(A,I5)') 'XX obcalc: not-yet supported mdlobi:',mdlobi
    END SELECT
       
    DO nobt=1,nobt_max
       CALL ob_exec(y_in,nobt,ierr)
       IF(ierr.NE.0) cycle
    ENDDO

    CALL GUTIME(time2)
    WRITE(6,*) '% CPU TIME = ',time2-time1,' sec'

9000 CONTINUE
    RETURN
  END SUBROUTINE ob_calc


END MODULE obcalc
