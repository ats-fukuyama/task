! plprofw.f90

MODULE plprofw

  USE bpsd_kinds,ONLY: rkind

  PRIVATE
  PUBLIC pl_prfw_type,pl_profw,pl_profw3d

  TYPE pl_prfw_type       ! local plasma parameter for wave analysis
     real(rkind):: RN,RTPR,RTPP,RUPR,RUPP,RNUC,RZCL
                         ! RN:   number density [10^{20}m^{-3}]
                         ! RTPR: parallel temperature [keV]
                         ! RTPP: perpendicular temperature [keV]
                         ! RUPR: parallel fluid velocity [m/s]
                         ! RUPP: perpendicular fluid velocity [m/s]
                         ! RNUC: collision frequency [1/s]
                         ! RZCL: collision parameter (RNUC/OMEGA)
  END TYPE pl_prfw_type

CONTAINS

  SUBROUTINE pl_profw(rhon,plfw)
    USE plcomm
    USE plcomm_type
    USE plprof
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rhon
    TYPE(pl_prfw_type),DIMENSION(nsmax),INTENT(OUT):: plfw
    TYPE(pl_prf_type),DIMENSION(nsmax):: plf
    REAL(rkind):: ql,bpbt,bnt,bnp
    INTEGER:: ns

    IF(mdlplw.EQ.1) THEN
       CALL pl_qprf(rhon,ql)
       bpbt=rhon*ra/(rr*ql)
       bnt=1.D0/SQRT(1.D0+bpbt**2)
       bnp=bnt*bpbt
    END IF

    CALL pl_prof(rhon,plf)
    DO ns=1,nsmax
       plfw(ns)%rn  =plf(ns)%rn
       plfw(ns)%rtpr=plf(ns)%rtpr
       plfw(ns)%rtpp=plf(ns)%rtpp
       SELECT CASE(mdlplw)
       CASE(0)
          plfw(ns)%rupr=pupr(ns)
          plfw(ns)%rupp=pupp(ns)
       CASE(1)
          plfw(ns)%rupr= plf(ns)%ru*bnt+plf(ns)%rupl*bnp
          plfw(ns)%rupp=-plf(ns)%ru*bnp+plf(ns)%rupl*bnt
       END SELECT
       plfw(ns)%rnuc=plf(ns)%rnuc
       SELECT CASE(model_coll)
       CASE(0)
          plfw(ns)%rzcl=PZCL(NS)
       CASE(1,2)
          plfw(ns)%rzcl=plf(ns)%rnuc/(2.D0*PI*RF_PL*1.D6)
       END SELECT
    END DO
  END SUBROUTINE pl_profw

  SUBROUTINE pl_profw3d(x,y,z,plfw)
    USE plcomm
    USE plprof
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x,y,z
    TYPE(pl_prfw_type),DIMENSION(nsmax),INTENT(OUT):: plfw
    TYPE(pl_prf_type),DIMENSION(nsmax):: plf
    TYPE(pl_mag_type):: mag
    REAL(rkind):: raxis,zaxis,rl,rcost,rsint,bnt,bnr,rs,rsinp,rcosp,bnp
    INTEGER:: ns

    IF(mdlplw.EQ.1) THEN
       CALL pl_mag(x,y,z,mag)
       CALL getaxs(raxis,zaxis)
       rl=SQRT(x**2+y**2)
       rcost= x/rl
       rsint= y/rl
       bnt=-mag%bnx*rsint+mag%bny*rcost
       bnr= mag%bnx*rcost+mag%bny*rsint
       rs   = SQRT((rl-raxis)**2+(z-zaxis)**2)
       rsinp= (z-zaxis)/rs
       rcosp= (rl-raxis)/rs
       bnp=-bnr*rsinp+mag%bnz*rcosp
       IF(bnp.GE.0.D0) THEN
          bnp= SQRT(bnr**2+mag%bnz**2)
       ELSE
          bnp=-SQRT(bnr**2+mag%bnz**2)
       END IF
    END IF

    CALL pl_prof3d(x,y,z,plf)
    DO ns=1,nsmax
       plfw(ns)%rn  =plf(ns)%rn
       plfw(ns)%rtpr=plf(ns)%rtpr
       plfw(ns)%rtpp=plf(ns)%rtpp
       SELECT CASE(mdlplw)
       CASE(0)
          plfw(ns)%rupr=pupr(ns)
          plfw(ns)%rupp=pupp(ns)
       CASE(1)
          plfw(ns)%rupr= plf(ns)%ru*bnt+plf(ns)%rupl*bnp
          plfw(ns)%rupp=-plf(ns)%ru*bnp+plf(ns)%rupl*bnt
       END SELECT
       plfw(ns)%rnuc=plf(ns)%rnuc
       SELECT CASE(model_coll)
       CASE(0)
          plfw(ns)%rzcl=PZCL(NS)
       CASE(1,2)
          plfw(ns)%rzcl=plf(ns)%rnuc/(2.D0*PI*RF_PL*1.D6)
       END SELECT
    END DO
  END SUBROUTINE pl_profw3d

END MODULE plprofw
