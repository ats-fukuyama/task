! fpcaldrrem.f90

MODULE fpcaldrrem

  ! *******************************************
  !     Radial transport due to EM turbulence
  ! *******************************************

  PRIVATE
  PUBLIC fp_caldrr_em   ! modeld=2

CONTAINS

  SUBROUTINE fp_caldrr_em

    USE fpcomm
    USE fpdebug
    USE plprof
    USE libell
    IMPLICIT NONE
    INTEGER:: NSA,NSSA,NR,NP,NTH,ierr
    
    REAL(rkind):: rhon,ql,rs,rk_para,rk_perp,pfpl,vfpl,rnu_b,v_gradB
    REAL(rkind):: sinml,cosml,v_para,factor,DRR_em_amp_local

    !---- Calculation of EM diffusion coefficient ----

    DO NSA=NSASTART,NSAEND
       NSSA=NS_NSA(NSA)

       DO NR=NRSTART,NRENDWG
          RHON=RG(NR)
          CALL pl_qprf(rhon,ql)
          RS=RA*RG(NR)
          rk_para=1.D0/(ql*RR)
          rk_perp=DRR_em_kdep*rk_para
          DRR_em_amp_local=DRR_em_amp*EXP(-(rhon-DRR_em_r0)**2/DRR_em_rw**2)
          
          DO NP=NPSTART,NPEND
             PFPL=PM(NP,NSSA)*PTFP0(NSA)
             VFPL=PFPL/AMFP(NSA)
             rnu_b=VFPL*SQRT(rs/RR)/(2.D0*ql*RR)
             v_gradB=VFPL*PFPL*rk_perp*DRR_em_amp_local/(2.D0*AEE*BB)

             DO NTH=1,NTHMAX
                SINML=SINM(NTH)
                COSML=COSM(NTH)

!                RAML=SINML**2*(1+EPSRM(NR))
!                RK=SQRT((1.D0-RAML*(1.D0-EPSRM(NR)))/(2.D0*EPSRM(NR)*RAML))
!                IF(RK.LT.1.D0) THEN
!                   factor_th=0.5D0*Pi*SINML**4/(COSML**2*ELLFC(RK,IERR))
!                   IF(ierr.NE.0) THEN
!                      WRITE(6,*) 'XX fp_cal_rd_em: ELLFC error, ierr=',ierr
!                      STOP
!                   END IF
!                ELSE
!                   factor_th=0.D0
!                END IF
!                  WRITE(61,'(4I4,5ES12.4)') NTH,NP,NR,NSA, &
!                       SINML,RAML,RK,ELLFC(RK,IERR),factor_th

                v_para=VFPL*COSML
                factor=rnu_b*SINML/(rk_para**2*v_para**2+rnu_b**2*SINML**2)
                DRR(NTH,NP,NR,NSA)=0.5D0*v_gradB**2*SINML**4*PI**2*factor
!                  WRITE(61,'(4I6,ES12.4)') &
!                       NSA,NR,NP,NTH,DRR(NTH,NP,NR,NSA)
                FRR(NTH,NP,NR,NSA)=0.D0
             END DO ! NTH
          END DO ! NP
       END DO ! NR
    END DO ! NSA

    CALL fp_debug_drr
    RETURN
  END SUBROUTINE fp_caldrr_em

  END MODULE fpcaldrrem
