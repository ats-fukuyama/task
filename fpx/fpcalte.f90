module fpcalte

CONTAINS
  SUBROUTINE fp_calte
    USE fpcomm
    USE fpsave
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:),ALLOCATABLE::TAUE,E,SPB
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE::rhot
    INTEGER:: NSA,NR

    ALLOCATE(rhot(NRMAX,NSAMAX)) !thermal energy density
    ALLOCATE(E(NSAMAX))
    ALLOCATE(SPB(NSAMAX)) !sum of power from beam
    ALLOCATE(TAUE(NSAMAX))
    E=0.D0
    SPB=0.D0
    rhot=0.D0

    CALL MOMENT_2ND_ORDER(FNSP,rhot)
    DO NSA=1,NSAMAX
      DO NR=1,NRMAX
        E(NSA)=E(NSA)+rhot(NR,NSA)*VOLR(NR)
      END DO
    END DO
    DO NSA=1,NSAMAX
      DO NR=1,NRMAX
        SPB(NSA)=SPB(NSA)+RSPBL(NR,NSA)*VOLR(NR)
      END DO
    END DO
    DO NSA=1,NSAMAX
      IF(SPB(NSA).EQ.0.D0) THEN
        TAUE(NSA)=0.D0
      ELSE
        TAUE(NSA)=E(NSA)/SPB(NSA)
      END IF
    END DO
    DEALLOCATE(rhot,SPB,TAUE,E)
  END SUBROUTINE fp_calte
END MODULE fpcalte
