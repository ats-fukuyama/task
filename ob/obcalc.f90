!   obcalc.f90

MODULE obcalc

  PRIVATE
  PUBLIC ob_calc

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE ob_calc(ierr)

    USE obcomm
!    USE obexec,ONLY: ob_exec
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL(4):: time1,time2
    INTEGER:: nobt

    ierr=0

    CALL GUTIME(time1)

    ! --- line input of initial parameters ---
    
    SELECT CASE(mdlobi)
    CASE(100)
       WRITE(6,'(A,I5)') '## Initial parameters: nobt_max=',nobt_max
       WRITE(6,'(A)') &
            'nobt  penergy_in  pangle_in   zeta_in     pzeta_in    theta_in'
       DO nobt=1,nobt_max
100       CONTINUE
          WRITE(6,'(I4,1P5E12.4)') &
               nobt,penergy_in(nobt),pangle_in(nobt),zeta_in(nobt), &
               pzeta_in(nobt),theta_in(nobt)
          READ(5,*,ERR=100,END=9000) &
               penergy_in(nobt),pangle_in(nobt),zeta_in(nobt), &
               pzeta_in(nobt),theta_in(nobt)
       END DO
    CASE(101)
       WRITE(6,'(A,I5)') '## Initial parameters: nobt_max=',nobt_max
       WRITE(6,'(A)') &
            'nobt  penergy_in  pangle_in   zeta_in     rr_in       zz_in'
       DO nobt=1,nobt_max
200       CONTINUE
          WRITE(6,'(I4,1P5E12.4)') &
               nobt,penergy_in(nobt),pangle_in(nobt),zeta_in(nobt), &
               rr_in(nobt),zz_in(nobt)
          READ(5,*,ERR=200,END=9000) &
               penergy_in(nobt),pangle_in(nobt),zeta_in(nobt), &
               rr_in(nobt),zz_in(nobt)
       END DO
    END SELECT

    ! calculate actual initial values ---
    
    SELECT CASE(mdlobi)
    CASE(0,100)
       DO nobt=1,nobt_max
          penergy_ob(0,nobt)=penergy_in(nobt)
          pangle_ob(0,nobt)=pangle_in(nobt)
          zetab_ob(0,nobt)=zeta_in(nobt)
          pzeta_ob(0,nobt)=pzeta_in(nobt)
          thetab_ob(0,nobt)=theta_in(nobt)

          ptheta_ob(0,nobt)=0
          psip_ob(0,nobt)=0
          rhopara_ob(0,nobt)=0
          
          babs_ob(0,nobt)=0 ! psip,thetab => babs,phi
          phi_ob(0,nobt)=0
          rr_ob(0,nobt)=0 ! psip,thetab => rr,zz
          zz_ob(0,nobt)=0
          zeta_ob(0,nobt)=zetab_ob(0,nobt)
          rs_ob(0,nobt)=SQRT((rr_ob(0,nobt)-rr_axis)**2 &
                            +(zz_ob(0,nobt)-zz_axis)**2)
          theta_ob(0,nobt)=ATAN2(zz_ob(0,nobt)-zz_axis,rr_ob(0,nobt)-rr_axis)
       END DO
    CASE(1,101)
       DO nobt=1,nobt_max
          penergy_ob(0,nobt)=penergy_in(nobt)
          pangle_ob(0,nobt)=pangle_in(nobt)
          zetab_ob(0,nobt)=zeta_in(nobt)
          rr_ob(0,nobt)=rr_in(nobt)
          zz_ob(0,nobt)=zz_in(nobt)

          psip_ob(0,nobt)=0  ! rr,zz => psip,thetab
          thetab_ob(0,nobt)=0

          pzeta_ob(0,nobt)=0
          ptheta_ob(0,nobt)=0
          rhopara_ob(0,nobt)=0
          
          babs_ob(0,nobt)=0 ! psip,thetab => babs,phi
          phi_ob(0,nobt)=0
          rr_ob(0,nobt)=0
          zz_ob(0,nobt)=0
          zeta_ob(0,nobt)=zetab_ob(0,nobt)
          rs_ob(0,nobt)=SQRT((rr_ob(0,nobt)-rr_axis)**2 &
                            +(zz_ob(0,nobt)-zz_axis)**2)
          theta_ob(0,nobt)=ATAN2(zz_ob(0,nobt)-zz_axis,rr_ob(0,nobt)-rr_axis)
       END DO
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
