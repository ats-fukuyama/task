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
    REAL(rkind):: RHON,DR_FIXED,rne,rte
    REAL(rkind):: DR_inz,DR_rcb,DR_prb,DR_plt
    INTEGER:: NSA,NEQ,NV,NS,ID,NPZ,NSA1,ierr

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
             PZA(NSA)=func_adpost(NPA(NS),1,rte)
             PZ2A(NSA)=func_adpost(NPA(NS),2,rte)
             PRADE(NSA,NR)=func_adpost(NPA(NS),3,rte)
          CASE(10,11,12) ! ionization and recombination with OPEN-ADAS
             NPZ=NINT(PZA(NSA))
             NSA1=NSA_UP(NSA)
             IF(NSA1.NE.0) THEN
                CALL ADAS_scd(NPA(NS),NPZ+1,rne,rte,DR_inz,IERR)
                CCN(NSA,NSA ,NR)=CCN(NSA,NSA ,NR) &           ! ionization
                                -DR_inz*rne                   !   to upper
                CALL ADAS_acd(NPA(NS),NPZ+1,rne,rte,DR_rcb,IERR)
                CCN(NSA,NSA1,NR)=CCN(NSA,NSA1,NR) &           ! recombination
                                +DR_rcb*rne                   !   from upper
             END IF
             NSA1=NSA_DN(NSA)
             IF(NSA1.NE.0) THEN
                CALL ADAS_acd(NPA(NS),NPZ,rne,rte,DR_rcb,IERR)
                CCN(NSA,NSA ,NR)=CCN(NSA,NSA ,NR) &           ! recombination
                                -DR_rcb*rne                   !     to lower
                CALL ADAS_scd(NPA(NS),NPZ,rne,rte,DR_inz,IERR)
                CCN(NSA,NSA1,NR)=CCN(NSA,NSA1,NR) &           ! ionization
                                +DR_inz*rne                   !     from lower
             END IF
             CALL ADAS_prb(NPA(NS),NPZ,rne,rte,DR_prb,IERR)   ! recmb/brems pwr
             CALL ADAS_plt(NPA(NS),NPZ,rne,rte,DR_plt,IERR)   ! line ard pwr
             PRADE(NSA,NR)=(10.D0**DR_prb+10.D0**DR_plt)*rne
          END SELECT
       END DO

    DR_FIXED=DR0+(DRS-DR0)*RHON**2

    DO NSA=1,nsa_max
       DNTB(NSA,NR)=DN0*DR_FIXED
       DTTB(NSA,NR)=DT0*DR_FIXED
       DUTB(NSA,NR)=DU0*DR_FIXED
       VDNTB(NSA,NR)=VDN0*DR_FIXED*RHON
       VDTTB(NSA,NR)=VDT0*DR_FIXED*RHON
       VDUTB(NSA,NR)=VDU0*DR_FIXED*RHON
    END DO

    DO NSA=1,nsa_max
       DNNC(NSA,NR)=0.D0
       DTNC(NSA,NR)=0.D0
       DUNC(NSA,NR)=0.D0
       VDNNC(NSA,NR)=0.D0
       VDTNC(NSA,NR)=0.D0
       VDUNC(NSA,NR)=0.D0
    END DO

    DO NSA=1,nsa_max
       DDN(NSA,NR)=DNTB(NSA,NR) +DNNC(NSA,NR)
       VVN(NSA,NR)=VDNTB(NSA,NR)+VDNNC(NSA,NR)
       DDT(NSA,NR)=DTTB(NSA,NR) +DTNC(NSA,NR)
       VVT(NSA,NR)=VDTTB(NSA,NR)+VDTNC(NSA,NR)
       DDU(NSA,NR)=DUTB(NSA,NR) +DUNC(NSA,NR)
       VVU(NSA,NR)=VDUTB(NSA,NR)+VDUNC(NSA,NR)
    END DO

    RETURN
  END SUBROUTINE ti_coef
END MODULE ticoef
