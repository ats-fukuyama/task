! MODULE ticoef

MODULE ticoef

  PUBLIC ti_coef

CONTAINS

  SUBROUTINE ti_coef(NR)
    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind):: RHON,DK_FIXED
    INTEGER:: NSA,NEQ,NV

    RHON=RM(NR)/RA
    DK_FIXED=DK0+(DKS-DK0)*RHON**2

    DO NSA=1,nsa_max
       AKTB(NSA,NR)=AK0*DK_FIXED
       ADTB(NSA,NR)=AD0*DK_FIXED
       AVTB(NSA,NR)=AV0*DK_FIXED*RHON
    END DO

    DO NSA=1,nsa_max
       AKNC(NSA,NR)=0.D0
       ADNC(NSA,NR)=0.D0
       AVNC(NSA,NR)=0.D0
    END DO

    DO NEQ=1,NEQMAX
       NSA=NSA_NEQ(NEQ)
       NV=NV_NEQ(NEQ)
       SELECT CASE(NV)
       CASE(1)
          DD(NEQ,NR)=ADTB(NSA,NR)+ADNC(NSA,NR)
          VV(NEQ,NR)=AVTB(NSA,NR)+AVNC(NSA,NR)
       CASE(2)
          DD(NEQ,NR)=0.D0
          VV(NEQ,NR)=0.D0
       CASE(3)
          DD(NEQ,NR)=AKTB(NSA,NR)+AKNC(NSA,NR)
          VV(NEQ,NR)=AVTB(NSA,NR)+AVNC(NSA,NR)
       END SELECT
    END DO

    RETURN
  END SUBROUTINE ti_coef
END MODULE ticoef
