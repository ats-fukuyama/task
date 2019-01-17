! MODULE ticoef

MODULE ticoef

  PUBLIC ti_coef

CONTAINS

  SUBROUTINE ti_coef(NR)
    USE ticomm
    USE adpost
    USE tiadas
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind):: RHON,DK_FIXED,rne,rte
    REAL(rkind):: DR_inz,DR_rcb,DR_prb,DR_plt
    INTEGER:: NSA,NEQ,NV,NS,ID,NPA,NPZ,NSA1,ierr

    RHON=RM(NR)/RA

    DO NSA=1,nsa_max
       DO NSA1=1,nsa_max
          CCN(NSA,NSA1,NR)=0.D0
          CCT(NSA,NSA1,NR)=0.D0
          CCU(NSA,NSA1,NR)=0.D0
       END DO
       PRADE(NSA,NR)=0.D0
    END DO

! --- ionization, recombination, and radiation ---
!     !!! presently assuming nsa=1 for electron
!         to be generalized in future 
    
    rne=RNA(1,NR)
    rte=RTA(1,NR)

    DO NSA=1,nsa_max
       NS=NS_NSA(NSA)
       ID=ID_NS(NS)
       SELECT CASE(ID)
          CASE(-1) ! electron
             PZA(NSA)=-1.D0
             PZ2A(NSA)=1.D0
             PRADE(NSA,NR)=0.D0
          CASE(0) ! neutal
             PRADE(NSA,NR)=0.D0
          CASE(1,2) ! ions
             PRADE(NSA,NR)=0.D0
          CASE(5,6) ! ADPOST PZ-variable ions
             NPA=NINT(PA(NS))
             PZA(NSA)=func_adpost(NPA,1,rte)
             PZ2A(NSA)=func_adpost(NPA,2,rte)
             PRADE(NSA,NR)=func_adpost(NPA,3,rte)
          CASE(10,11,12) ! ionization and recombination with OPEN-ADAS
             NPA=NINT(PA(NS))
             NPZ=NINT(PZA(NSA))
             NSA1=NSA_UP(NSA)
             IF(NSA1.NE.0) THEN
                CALL ADAS_scd(NPA,NPZ+1,rne,rte,DR_inz,IERR)
                CCN(NSA,NSA ,NR)=CCN(NSA,NSA ,NR) &           ! ionization
                                -DR_inz*rne            !   to upper
                CALL ADAS_acd(NPA,NPZ+1,rne,rte,DR_rcb,IERR)
                CCN(NSA,NSA1,NR)=CCN(NSA,NSA1,NR) &           ! recombination
                                +DR_rcb*rne            !   from upper
             END IF
             NSA1=NSA_DN(NSA)
             IF(NSA1.NE.0) THEN
                CALL ADAS_acd(NPA,NPZ,rne,rte,DR_rcb,IERR)
                CCN(NSA,NSA ,NR)=CCN(NSA,NSA ,NR) &           ! recombination
                                -DR_rcb*rne            !     to lower
                CALL ADAS_acd(NPA,NPZ,rne,rte,DR_inz,IERR)
                CCN(NSA,NSA1,NR)=CCN(NSA,NSA1,NR) &           ! ionization
                                +DR_inz*rne            !     from lower
             END IF
             CALL ADAS_prb(NPA,NPZ,rne,rte,DR_prb,IERR)       ! recmb/brems pwr
             CALL ADAS_plt(NPA,NPZ,rne,rte,DR_plt,IERR)       ! line ard pwr
             PRADE(NSA,NR)=(10.D0**DR_prb+10.D0**DR_plt)*rne
          END SELECT
       END DO

    DK_FIXED=DK0+(DKS-DK0)*RHON**2

    DO NSA=1,nsa_max
       ADTB(NSA,NR)=AD0*DK_FIXED
       AKTB(NSA,NR)=AK0*DK_FIXED
       AVTB(NSA,NR)=AV0*DK_FIXED*RHON
    END DO

    DO NSA=1,nsa_max
       AKNC(NSA,NR)=0.D0
       ADNC(NSA,NR)=0.D0
       AVNC(NSA,NR)=0.D0
    END DO

    DO NSA=1,nsa_max
       DDN(NSA,NR)=ADTB(NSA,NR)+ADNC(NSA,NR)
       VVN(NSA,NR)=AVTB(NSA,NR)+AVNC(NSA,NR)
       DDT(NSA,NR)=AKTB(NSA,NR)+AKNC(NSA,NR)
       VVT(NSA,NR)=AVTB(NSA,NR)+AVNC(NSA,NR)
       DDU(NSA,NR)=0.D0
       VVU(NSA,NR)=0.D0
    END DO

    RETURN
  END SUBROUTINE ti_coef
END MODULE ticoef
