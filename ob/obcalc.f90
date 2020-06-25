!   obcalc.f90

MODULE obcalc

  PRIVATE
  PUBLIC ob_calc

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE ob_calc(ierr)

    USE obcomm
    USE obprep
    USE obexec,ONLY: ob_exec
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL(4):: time1,time2
    INTEGER:: nobt
    REAL(rkind):: &
         zetab_pos,thetab_pos,psip_pos,rho_para_pos,ptheta_pos,pzeta_pos, &
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

          CALL cal_rr_pos(thetab_pos,psip_pos,rr_pos,ierr)
          CALL cal_zz_pos(thetab_pos,psip_pos,zz_pos,ierr)
          CALL cal_bb_pos(thetab_pos,psip_pos,bb_pos,ierr)
          CALL cal_psit_pos(psip_pos,psit_pos,dummy,ierr)
          CALL cal_rbps_pos(psip_pos,fg_pos,dummy,ierr)
          CALL cal_ritps_pos(psip_pos,fI_pos,dummy,ierr)

          pze=PZ(ns_ob)*AEE
          v_para=SQRT(2.D0*penergy_pos*AEE*1.D3/(PA(ns_ob)*AMP))*pcangle_pos
          omega_c=pze*bb_pos/(PA(ns_ob)*AMP)
          rho_para_pos=v_para/omega_c

          ptheta_pos=pze*(psit_pos+rho_para_pos*fI_pos)
          pzeta_pos= pze*(rho_para_pos*fg_pos-psip_pos)
          
          time_ob(0,nobt)   = 0.D0
          zetab_ob(0,nobt)  = zetab_pos
          thetab_ob(0,nobt) = thetab_pos
          psip_ob(0,nobt)   = psip_pos
          rhopara_ob(0,nobt)= rho_para_pos
          pzeta_ob(0,nobt)   = pzeta_pos
          ptheta_ob(0,nobt)  = ptheta_pos
          babs_ob(0,nobt)   = bb_pos
          phi_ob(0,nobt)    = 0.D0
          penergy_ob(0,nobt)   = penergy_pos
          pcangle_ob(0,nobt)   = pcangle_pos
          psit_ob(0,nobt)   = psit_pos
          rr_ob(0,nobt)      = rr_pos
          zz_ob(0,nobt)      = zz_pos
          zeta_ob(0,nobt)   = zetab_pos
          rs_ob(0,nobt)     =SQRT((rr_pos-rr_axis)**2+(zz_pos-zz_axis)**2)
          theta_ob(0,nobt)  =ATAN2(zz_pos-zz_axis,rr_pos-rr_axis)
       END DO
    CASE(1,101)
       WRITE(6,'(A,I5)') 'XX obcalc: not-yet supported mdlobi:',mdlobi
    END SELECT
       
    DO nobt=1,nobt_max
       CALL ob_exec(nobt,ierr)
       IF(ierr.NE.0) cycle
    ENDDO

    CALL GUTIME(time2)
    WRITE(6,*) '% CPU TIME = ',time2-time1,' sec'

9000 CONTINUE
    RETURN
  END SUBROUTINE ob_calc


END MODULE obcalc
