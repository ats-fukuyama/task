! wqexec.f90

MODULE wqexec

  PRIVATE
  PUBLIC wq_exec

CONTAINS

  SUBROUTINE wq_exec
    USE wqcomm
    USE wqsolv1
    USE wqsolv2
    USE libsmooth
    IMPLICIT NONE
    INTEGER:: nx,ny,nt,ngt,ngr
    REAL(rkind):: t,pabs_tot,factor_pulse,factor_prof
    REAL(rkind):: y0,yn,x0,xn,ramp_length_temp

    t=t_tot
    nt=nt_tot
    ngt=ngt_max
    ngr=ngr_max
    
    DO nt=1,ntmax

       SELECT CASE(model_source)
       CASE(1:4) ! surface source
          SELECT CASE(model_pulse)
          CASE(0) ! flat
             SELECT CASE(model_ramp)
             CASE(0) ! step
                factor_pulse=1.D0
             CASE(1) ! linear 
                IF(t.LE.ramp_length*period) THEN
                   factor_pulse=t/(ramp_length*period)
                ELSE 
                   factor_pulse=1.D0
                END IF
             CASE(2) ! smooth
                IF(t.LE.ramp_length*period) THEN
                   factor_pulse=smooth_1(t/(ramp_length*period))
                ELSE
                   factor_pulse=1.D0
                END IF
             END SELECT
          CASE(1) ! pulse
             IF(t.LE.pulse_length*period) THEN
                SELECT CASE(model_ramp)
                CASE(0)
                   factor_pulse=1.D0
                CASE(1,2)
                   IF(ramp_length.LE.0.5D0*pulse_length) then
                      ramp_length_temp=0.5D0*pulse_length
                   ELSE
                      ramp_length_temp=ramp_length
                   END IF
                   IF(t.LE.ramp_length*period) THEN
                      SELECT CASE(model_ramp)
                      CASE(1)
                         factor_pulse=t/(ramp_length_temp*period)
                      CASE(2)
                         factor_pulse=smooth_1(t/(ramp_length_temp*period))
                      END SELECT
                   ELSE IF(t.LE.(pulse_length-ramp_length_temp)*period) THEN
                      factor_pulse=1.D0
                   ELSE
                      SELECT CASE(model_ramp)
                      CASE(1)
                         factor_pulse=(pulse_length*period-t) &
                              /(ramp_length_temp*period)
                      CASE(2)
                         factor_pulse=smooth_1((pulse_length*period-t) &
                              /(ramp_length_temp*period))
                      END SELECT
                   END IF
                END SELECT
             ELSE
                factor_pulse=0.D0
             END IF
          END SELECT

          SELECT CASE(model_source)
          CASE(1,3,5,7)
             SELECT CASE(model_source)
             CASE(1,5)
                nx=1
             CASE(3,7)
                nx=nxmax
             END SELECT
             y0=source_position_yn
             DO ny=1,nymax
                yn=yn_ny(ny)
                factor_prof=exp(-(yn-y0)**2/source_width**2)
                SELECT CASE(model_source)
                CASE(1,3)
                   EY(nx,ny) = factor_pulse*factor_prof
                CASE(5,7)
                   EZ(nx,ny) = factor_pulse*factor_prof
                END SELECT
             END DO
          CASE(2,4,6,8)
             SELECT CASE(model_source)
             CASE(2,6)
                ny=1
             CASE(4,8)
                ny=nymax
             END SELECT
             x0=source_position_xn
             DO nx=1,nxmax
                xn=xn_nx(nx)
                factor_prof=exp(-(xn-x0)**2/source_width**2)
                SELECT CASE(model_source)
                CASE(2,4)
                   EY(nx,ny) = factor_pulse*factor_prof
                CASE(6,8)
                   EZ(nx,ny) = factor_pulse*factor_prof
                END SELECT
             END DO
          END SELECT
       END SELECT

       !compute next E

       SELECT CASE(model_solver)
       CASE(1)
          call wq_solv1
       CASE(2)
          call wq_solv2
       CASE DEFAULT
          WRITE(6,*) 'XX wq_exec: unknown model_solver:',model_solver
       END SELECT

       t=t+dt
       nt_tot=nt_tot+1

       ! calculate Absorbed Power when necessary
          
       IF(MOD(nt-1,ntstep).EQ.0.OR. &
          MOD(nt-1,ngtstep).EQ.0.OR. &
          nt.EQ.ntmax)  THEN

          pabs_tot=0.D0

          DO ny=1,nymax
             pabs(    1,ny)=0.D0
             pabs(nxmax,ny)=0.D0
          END DO
          DO nx=1,nxmax
             pabs(nx,    1)=0.D0
             pabs(nx,nymax)=0.D0
          END DO
          
          DO ny=2,nymax-1
             DO nx=2,nxmax-1
                pabs(nx,ny)=0.5d0*real( &
                     conjg(EX(nx,ny))*(CD(1,1,nx,ny)*EX(nx,ny) &
                                      +CD(1,2,nx,ny)*EY(nx,ny) &
                                      +CD(1,3,nx,ny)*EZ(nx,ny))&
                    +conjg(EY(nx,ny))*(CD(2,1,nx,ny)*EX(nx,ny) &
                                      +CD(2,2,nx,ny)*EY(nx,ny) &
                                      +CD(2,3,nx,ny)*EZ(nx,ny))&
                    +conjg(EZ(nx,ny))*(CD(3,1,nx,ny)*EX(nx,ny) &
                                      +CD(3,2,nx,ny)*EY(nx,ny) &
                                      +CD(3,3,nx,ny)*EZ(nx,ny)))
                pabs_tot=pabs_tot+pabs(nx,ny)
             END DO
          END DO
       END IF

       IF(MOD(nt,ntstep).EQ.0.OR.nt.EQ.ntmax) &
          WRITE(6,'(A,I8,2ES12.4)') '## nt,t,pabs_tot:',nt,t,pabs_tot

       IF(MOD(nt,ngtstep).EQ.0.AND.ngt.LT.ngt_m) THEN
          ngt=ngt+1
          t_ngt(ngt)=t
          pabs_tot_ngt(ngt)=pabs_tot
          EX_max_ngt(ngt)=0.D0
          EY_max_ngt(ngt)=0.D0
          EY_max_ngt(ngt)=0.D0
          DO ny=1,nymax
             DO nx=1,nxmax
                EX_max_ngt(ngt)=MAX(EX_max_ngt(ngt),ABS(EX(nx,ny)))
                EY_max_ngt(ngt)=MAX(EY_max_ngt(ngt),ABS(EY(nx,ny)))
                EZ_max_ngt(ngt)=MAX(EZ_max_ngt(ngt),ABS(EZ(nx,ny)))
             END DO
          END DO
!          WRITE(6,'(A,2I6,4ES12.4)') &
!               'nt,ngt,Ex/y/z,pabs=',nt,ngt,EX_max_ngt(ngt), &
!               EY_max_ngt(ngt),EZ_max_ngt(ngt),pabs_tot_ngt(ngt)
       END IF

       IF(MOD(nt,ngrstep).EQ.0.AND.ngr.LT.ngr_m) THEN
          ngr=ngr+1
          t_ngr(ngr)=t
          DO ny=1,nymax
             DO nx=1,nxmax
                EX_save(nx,ny,ngr)=EX(nx,ny)
                EY_save(nx,ny,ngr)=EY(nx,ny)
                EZ_save(nx,ny,ngr)=EZ(nx,ny)
                pabs_save(nx,ny,ngr)=pabs(nx,ny)
             END DO
          END DO
          WRITE(6,'(A,I8,ES12.4,I8)') '## nt,t,ngr:     ',nt,t_ngr(ngr),ngr
       END IF
    END DO
    t_tot=t
    ngt_max=ngt
    ngr_max=ngr

    RETURN

  END SUBROUTINE wq_exec
END MODULE wqexec
