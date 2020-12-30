! wrcomm.f90

MODULE wrcomm_parm
  USE commpi
  USE plcomm
  USE dpcomm_parm
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER:: NRAYM=100

! --- input parameters ---

  INTEGER:: NRAYMAX,NSTPMAX,NRSMAX,NRLMAX,LMAXNW
  INTEGER:: MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW
  REAL(rkind) SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW
  REAL(rkind) RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,UUI,RKRI,RKPHII,RKZI
  REAL(rkind) RCURVA,RCURVB,RBRADA,RBRADB,NRADMX
  INTEGER:: MODEW,nres_max,nres_type

  REAL(rkind),DIMENSION(NRAYM):: &
       RFIN,RPIN,ZPIN,PHIIN,RKRIN,RNZIN,RNPHIIN,ANGZIN,ANGPHIN,UUIN, &
       RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN
  INTEGER,DIMENSION(NRAYM):: &
       MODEWIN

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
       RAYIN
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       RAYS
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
       pos_nrs(:),pwr_nrs(:),pwr_nrs_nray(:,:)
  REAL(rkind),ALLOCATABLE:: &
       pos_nrl(:),pwr_nrl(:),pwr_nrl_nray(:,:)
  REAL(rkind),ALLOCATABLE:: &
       rs_nstp_nray(:,:),rl_nstp_nray(:,:)
  REAL(rkind),ALLOCATABLE:: &
       pos_pwrmax_rs_nray(:),pwrmax_rs_nray(:)
  REAL(rkind):: pos_pwrmax_rs,pwrmax_rs
  REAL(rkind),ALLOCATABLE:: &
       pos_pwrmax_rl_nray(:),pwrmax_rl_nray(:)
  REAL(rkind):: pos_pwrmax_rl,pwrmax_rl
CONTAINS

  SUBROUTINE wr_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: INIT=0
    INTEGER,SAVE:: NRAYMAX_SAVE=0
    INTEGER,SAVE:: NSTPMAX_SAVE=0

    IF(INIT.EQ.0) THEN
       INIT=1
    ELSE
       IF((NRAYMAX.EQ.NRAYMAX_SAVE).AND. &
          (NSTPMAX.EQ.NSTPMAX_SAVE)) RETURN
       CALL wr_deallocate
    END IF

    ALLOCATE(RAYIN(NEQ,NRAYMAX))
    ALLOCATE(NSTPMAX_NRAY(NRAYMAX))
    ALLOCATE(RAYS(0:NEQ,0:NSTPMAX,NRAYMAX))

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
  END SUBROUTINE wr_allocate

  SUBROUTINE wr_deallocate
    IMPLICIT NONE

    DEALLOCATE(NSTPMAX_NRAY)
    DEALLOCATE(RAYIN)
    DEALLOCATE(RAYS)
    DEALLOCATE(CEXS,CEYS,CEZS)
    DEALLOCATE(RKXS,RKYS,RKZS,RXS,RYS,RZS,BNXS,BNYS,BNZS,BABSS)
    DEALLOCATE(RAYB,RAYRB1,RAYRB2)
    DEALLOCATE(CEXB,CEYB,CEZB)
    DEALLOCATE(RK1B,RP1B)
    DEALLOCATE(RK2B,RP2B,RAMPB)
  END SUBROUTINE wr_deallocate
END MODULE wrcomm
