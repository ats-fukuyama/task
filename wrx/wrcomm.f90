! wrcomm.f90

MODULE wrcomm_parm
  USE commpi
  USE plcomm
  USE dpcomm_parm
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER:: NRAYM=100
  INTEGER,PARAMETER:: idebug_max=99

! --- input parameters ---

  INTEGER:: model_fdrv,model_fdrv_ds
  INTEGER:: NRAYMAX,NSTPMAX,NRSMAX,NRLMAX,LMAXNW
  INTEGER:: mode_beam
  INTEGER:: MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW
  REAL(rkind):: SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW,EPSD0
  REAL(rkind):: pne_threshold,bdr_threshold
  INTEGER:: nsumax
  REAL(rkind):: Rmax_eq,Rmin_eq,Zmin_eq,Zmax_eq,Raxis_eq,Zaxis_eq
  REAL(rkind):: Rmax_wr,Rmin_wr,Zmin_wr,Zmax_wr
  REAL(rkind):: ra_wr
  INTEGER:: nsamax_wr,ns_nsa_wr(nsm),nsa_grf
  REAL(rkind),ALLOCATABLE:: rsu_wr(:),zsu_wr(:)
  INTEGER:: nres_max,nres_type

  REAL(rkind),DIMENSION(NRAYM):: &
       RFIN,RPIN,ZPIN,PHIIN,ANGTIN,ANGPIN,RNPHIN,RNZIN,RNKIN,UUIN, &
       RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN
  INTEGER,DIMENSION(NRAYM):: &
       MODEWIN
  CHARACTER(len=80):: KNAMWRW

  INTEGER,DIMENSION(idebug_max):: idebug_wr

CONTAINS

  SUBROUTINE open_wrcomm_parm
    RETURN
  END SUBROUTINE open_wrcomm_parm

END MODULE wrcomm_parm

MODULE wrcomm
  USE wrcomm_parm
  USE dpcomm
  USE commpi
  IMPLICIT NONE

  INTEGER,PARAMETER:: NEQ=8
  INTEGER,PARAMETER:: NBEQ=19
  INTEGER,PARAMETER:: NBVAR=53

  INTEGER,DIMENSION(:),ALLOCATABLE:: &
       NSTPMAX_NRAY
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       RAYIN,pwr_nsa_nstp
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       RAYS,pwr_nsa_nstp_nray
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       CEXS,CEYS,CEZS
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       RKXS,RKYS,RKZS,RXS,RYS,RZS,BNXS,BNYS,BNZS,BABSS
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       RAYB,RAYRB1,RAYRB2
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: &
       CEXB,CEYB,CEZB
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       RK1B,RP1B
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       RK2B,RP2B
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       RAMPB
  REAL(rkind),ALLOCATABLE:: &
       rs_nstp_nray(:,:),rl_nstp_nray(:,:)
  REAL(rkind),ALLOCATABLE:: &
       pos_nrs(:),pwr_nrs(:,:),pwr_nrs_nray(:,:,:)
  REAL(rkind),ALLOCATABLE:: &
       pos_nrl(:),pwr_nrl(:,:),pwr_nrl_nray(:,:,:)
  REAL(rkind),ALLOCATABLE:: &
       pwr_nsa_nray(:,:),pwr_nsa(:),pwr_nray(:)
  REAL(rkind),ALLOCATABLE:: &
       pos_pwrmax_rs_nray(:,:),pwrmax_rs_nray(:,:)
  REAL(rkind),ALLOCATABLE:: &
       pos_pwrmax_rs(:),pwrmax_rs(:)
  REAL(rkind),ALLOCATABLE:: &
       pos_pwrmax_rl_nray(:,:),pwrmax_rl_nray(:,:)
  REAL(rkind),ALLOCATABLE:: &
       pos_pwrmax_rl(:),pwrmax_rl(:)
  reaL(rkind):: pwr_tot
  INTEGER:: INITIAL_DRV

  TYPE wr_nray_status_type
     INTEGER:: nray,nstp
     REAL(rkind):: RF,RCURVA,RCURVB,RBRADA,RBRADB
  END type wr_nray_status_type
  TYPE(wr_nray_status_type):: wr_nray_status
     
CONTAINS

  SUBROUTINE wr_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: INIT=0
    INTEGER,SAVE:: NRAYMAX_SAVE=0
    INTEGER,SAVE:: NSTPMAX_SAVE=0
    INTEGER,SAVE:: nrsmax_save=0
    INTEGER,SAVE:: nrlmax_save=0
    INTEGER,SAVE:: nsamax_save=0

    IF(INIT.EQ.0) THEN
       INIT=1
    ELSE
       IF((NRAYMAX.EQ.NRAYMAX_SAVE).AND. &
          (NSTPMAX.EQ.NSTPMAX_SAVE).AND. &
          (nrsmax.EQ.nrsmax_save).AND. &
          (nrlmax.EQ.nrlmax_save).AND. &
          (nsamax_wr.EQ.nsamax_save)) RETURN
       CALL wr_deallocate
    END IF

    ALLOCATE(RAYIN(NEQ,NRAYMAX))
    ALLOCATE(NSTPMAX_NRAY(NRAYMAX))
    ALLOCATE(pwr_nsa_nstp(NSAMAX_WR,0:NSTPMAX))
    ALLOCATE(RAYS(0:NEQ,0:NSTPMAX,NRAYMAX))
    ALLOCATE(pwr_nsa_nstp_nray(NSAMAX_WR,0:NSTPMAX,NRAYMAX))

    ALLOCATE(CEXS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(CEYS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(CEZS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RKXS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RKYS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RKZS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RXS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RYS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RZS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(BNXS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(BNYS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(BNZS(0:NSTPMAX,NRAYMAX))
    ALLOCATE(BABSS(0:NSTPMAX,NRAYMAX))

    ALLOCATE(RAYB(0:NBVAR,0:NSTPMAX))
    ALLOCATE(RAYRB1(0:NSTPMAX,NRAYMAX))
    ALLOCATE(RAYRB2(0:NSTPMAX,NRAYMAX))
    ALLOCATE(CEXB(0:NSTPMAX),CEYB(0:NSTPMAX),CEZB(0:NSTPMAX))
    ALLOCATE(RK1B(3,0:NSTPMAX),RP1B(3,0:NSTPMAX))
    ALLOCATE(RK2B(3,3,0:NSTPMAX),RP2B(3,3,0:NSTPMAX))
    ALLOCATE(RAMPB(0:NSTPMAX))

    ALLOCATE(rs_nstp_nray(nstpmax,nraymax),rl_nstp_nray(nstpmax,nraymax))
    ALLOCATE(pos_nrs(nrsmax),pwr_nrs(nsamax_wr,nrsmax))
    ALLOCATE(pwr_nrs_nray(nsamax_wr,nrsmax,nraymax))
    ALLOCATE(pos_nrl(nrlmax),pwr_nrl(nsamax_wr,nrlmax))
    ALLOCATE(pwr_nrl_nray(nsamax_wr,nrlmax,nraymax))
    ALLOCATE(pwr_nsa_nray(nsamax_wr,nraymax))
    ALLOCATE(pwr_nray(nraymax),pwr_nsa(nsamax_wr))
    ALLOCATE(pos_pwrmax_rs_nray(nsamax_wr,nraymax))
    ALLOCATE(pwrmax_rs_nray(nsamax_wr,nraymax))
    ALLOCATE(pos_pwrmax_rl_nray(nsamax_wr,nraymax))
    ALLOCATE(pwrmax_rl_nray(nsamax_wr,nraymax))

  END SUBROUTINE wr_allocate

  SUBROUTINE wr_deallocate
    IMPLICIT NONE

    DEALLOCATE(NSTPMAX_NRAY)
    DEALLOCATE(RAYIN,pwr_nsa_nstp)
    DEALLOCATE(RAYS,pwr_nsa_nstp_nray)
    DEALLOCATE(CEXS,CEYS,CEZS)
    DEALLOCATE(RKXS,RKYS,RKZS,RXS,RYS,RZS,BNXS,BNYS,BNZS,BABSS)
    DEALLOCATE(RAYB,RAYRB1,RAYRB2)
    DEALLOCATE(CEXB,CEYB,CEZB)
    DEALLOCATE(RK1B,RP1B)
    DEALLOCATE(RK2B,RP2B,RAMPB)

    DEALLOCATE(rs_nstp_nray,rl_nstp_nray)
    DEALLOCATE(pos_nrs,pwr_nrs)
    DEALLOCATE(pwr_nrs_nray)
    DEALLOCATE(pos_nrl,pwr_nrl)
    DEALLOCATE(pwr_nrl_nray)
    DEALLOCATE(pwr_nsa_nray,pwr_nray,pwr_nsa)
    DEALLOCATE(pos_pwrmax_rs_nray)
    DEALLOCATE(pwrmax_rs_nray)
    DEALLOCATE(pos_pwrmax_rl_nray)
    DEALLOCATE(pwrmax_rl_nray)

  END SUBROUTINE wr_deallocate
END MODULE wrcomm
