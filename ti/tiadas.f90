! tiadas.f90

MODULE tiadas

PRIVATE
PUBLIC ADAS_acd,ADAS_scd,ADAS_ccd,ADAS_prb,ADAS_prc,ADAS_plt

!    ADAS class
!                class index    type      GCR data content
!                -----------    ----      ----------------
!                    1          acd       recombination coeffts
!                    2          scd       ionisation coeffts
!                    3          ccd       CX recombination coeffts
!                    4          prb       recomb/brems power coeffts
!                    5          prc       CX power coeffts
!                    6          qcd       base meta. coupl. coeffts
!                    7          xcd       parent meta. coupl. coeffts
!                    8          plt       low level line power coeffts
!                    9          pls       represent. line power coefft
!                   10          zcd       effective charge
!                   11          ycd       effective squared charge
!                   12          ecd       effective ionisation potential

CONTAINS

! --- recombination

  SUBROUTINE ADAS_acd(NZ0,NZ,PN,PT,DR,IERR)

    USE ADF11
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ0,NZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND
    REAL(dp):: DRL

    ND=ND_TABLE(NZ0,1)
    IF(ND.EQ.0) THEN
       IERR=1
       RETURN
    END IF
    IF(PN.LE.0.D0) THEN
       IERR=2
       DR=0.D0
       RETURN
    END IF
    IF(PT.LE.0.D0) THEN
       IERR=3
       DR=0.D0
       RETURN
    END IF

    CALL CALC_ADF11(ND,NZ,LOG10(PN),LOG10(PT),DRL,IERR)
    DR=10.D0**DRL
    RETURN
  END SUBROUTINE ADAS_acd

! --- ionization

  SUBROUTINE ADAS_scd(NZ0,NZ,PN,PT,DR,IERR)

    USE ADF11
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ0,NZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND
    REAL(dp):: DRL

    ND=ND_TABLE(NZ0,2)
    IF(ND.EQ.0) THEN
       IERR=1
       RETURN
    END IF
    IF(PN.LE.0.D0) THEN
       IERR=2
       DR=0.D0
       RETURN
    END IF
    IF(PT.LE.0.D0) THEN
       IERR=3
       DR=0.D0
       RETURN
    END IF

    CALL CALC_ADF11(ND,NZ,LOG10(PN),LOG10(PT),DRL,IERR)
    DR=10.D0**DRL
    RETURN
  END SUBROUTINE ADAS_scd

! --- CX recombination

  SUBROUTINE ADAS_ccd(NZ0,NZ,PN,PT,DR,IERR)

    USE ADF11
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ0,NZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND
    REAL(dp):: DRL

    ND=ND_TABLE(NZ0,3)
    IF(ND.EQ.0) THEN
       IERR=1
       RETURN
    END IF
    IF(PN.LE.0.D0) THEN
       IERR=2
       DR=0.D0
       RETURN
    END IF
    IF(PT.LE.0.D0) THEN
       IERR=3
       DR=0.D0
       RETURN
    END IF

    CALL CALC_ADF11(ND,NZ,LOG10(PN),LOG10(PT),DRL,IERR)
    DR=10.D0**DRL
    RETURN
  END SUBROUTINE ADAS_ccd

! --- recombination/bremsstrahlung power

  SUBROUTINE ADAS_prb(NZ0,NZ,PN,PT,DR,IERR)

    USE ADF11
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ0,NZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND

    ND=ND_TABLE(NZ0,4)
    IF(ND.EQ.0) THEN
       IERR=1
       RETURN
    END IF
    IF(PN.LE.0.D0) THEN
       IERR=2
       DR=0.D0
       RETURN
    END IF
    IF(PT.LE.0.D0) THEN
       IERR=3
       DR=0.D0
       RETURN
    END IF

    CALL CALC_ADF11(ND,NZ,LOG10(PN),LOG10(PT),DR,IERR)
    RETURN
  END SUBROUTINE ADAS_prb

! --- CX power

  SUBROUTINE ADAS_prc(NZ0,NZ,PN,PT,DR,IERR)

    USE ADF11
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ0,NZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND

    ND=ND_TABLE(NZ0,5)
    IF(ND.EQ.0) THEN
       IERR=1
       RETURN
    END IF
    IF(PN.LE.0.D0) THEN
       IERR=2
       DR=0.D0
       RETURN
    END IF
    IF(PT.LE.0.D0) THEN
       IERR=3
       DR=0.D0
       RETURN
    END IF

    CALL CALC_ADF11(ND,NZ,LOG10(PN),LOG10(PT),DR,IERR)
    RETURN
  END SUBROUTINE ADAS_prc

! --- low level line power

  SUBROUTINE ADAS_plt(NZ0,NZ,PN,PT,DR,IERR)

    USE ADF11
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ0,NZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND

    ND=ND_TABLE(NZ0,8)
    IF(ND.EQ.0) THEN
       IERR=1
       DR=0.D0
       RETURN
    END IF
    IF(PN.LE.0.D0) THEN
       IERR=2
       DR=0.D0
       RETURN
    END IF
    IF(PT.LE.0.D0) THEN
       IERR=3
       DR=0.D0
       RETURN
    END IF

    CALL CALC_ADF11(ND,NZ,LOG10(PN),LOG10(PT),DR,IERR)
    RETURN
  END SUBROUTINE ADAS_plt

END MODULE tiadas
