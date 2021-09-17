MODULE fpcaldeff

CONTAINS
  SUBROUTINE fp_caldeff
    USE fpcomm
    USE libgrf
    USE fpsave
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE::DNDT,GAMMAP,DEFF,DNDR,PST,GAMMAPR
    REAL(rkind):: DTG2
    INTEGER:: NSA,NR,NP,NTH

    ALLOCATE(DNDT(NRMAX,NSAMAX))
    ALLOCATE(PST(NRMAX,NSAMAX)) !particle source term
    ALLOCATE(GAMMAPR(NRMAX+1,NSAMAX))!particle flux * minor radius
    ALLOCATE(GAMMAP(NRMAX+1,NSAMAX))
    ALLOCATE(DNDR(NRMAX+1,NSAMAX))
    ALLOCATE(DEFF(NRMAX+1,NSAMAX))
    PST=0.D0

    IF(NTG1.EQ.NTG2) THEN
       IF(NTG2.GT.1) THEN
          DTG2=PTG(NTG2)-PTG(NTG2-1)
          DO NSA=1,NSAMAX
             DO NR=1,NRMAX
                DNDT(NR,NSA)=(RNT(NR,NSA,NTG2)-RNT(NR,NSA,NTG2-1))/DTG2
             END DO
          END DO
       ELSE
          DO NSA=1,NSAMAX
             DO NR=1,NRMAX
                DNDT(NR,NSA)=0.D0
             END DO
          END DO
       END IF
    ELSE IF(NTG1.GE.NTG2) THEN
       IF(NTG2.GT.0) THEN
          DTG2=PTG(NTG1)-PTG(NTG2)
          DO NSA=1,NSAMAX
             DO NR=1,NRMAX
                DNDT(NR,NSA)=(RNT(NR,NSA,NTG1)-RNT(NR,NSA,NTG2))/DTG2
             END DO
          END DO
       ELSE
          DO NSA=1,NSAMAX
             DO NR=1,NRMAX
                DNDT(NR,NSA)=0.D0
             END DO
          END DO
       END IF
    END IF
    DO NTH=1,NTHMAX
      DO NP=1,NPMAX
        PST(NR,NSA)=PST(NR,NSA)&
                    +(SPPF(NTH,NP,NR,NSA) &
                    +SPPB(NTH,NP,NR,NSA) &
                    +SPPL(NTH,NP,NR,NSA) &
                    +SPPS(NTH,NP,NR,NSA) &
                    +SPPI(NTH,NP,NR,NSA) &
                    ! +SPPD(NTH,NP,NSA)&
                    +PPL(NTH,NP,NR,NSA)&
                    )*VOLP(NTH,NP,NSA)
      END DO
    END DO
    DO NSA=1,NSAMAX
      GAMMAPR(1,NSA)=0.D0
      GAMMAP(1,NSA)=0.D0
      DO NR=1,NRMAX
        GAMMAPR(NR+1,NSA)=GAMMAPR(NR,NSA)+RA*(NR/NRMAX)*(DNDT(NR,NSA)-PST(NR,NSA))*(RG(NR+1)-RG(NR))
        GAMMAP(NR+1,NSA)=GAMMAPR(NR+1,NSA)/RA*(NR/NRMAX)
      END DO
    END DO
    CALL MOMENT_0TH_ORDER(FNSP,RNSL)
    DO NSA=1,NSAMAX
      DNDR(1,NSA)=0.D0
      DO NR=1,NRMAX
        DNDR(NR+1,NSA)=(RNSL(NR+1,NSA)-RNSL(NR,NSA))/(RG(NR+1)-RG(NR))
      END DO
    END DO
    DO NSA=1,NSAMAX
       DEFF(1,NSA)=0.D0
       DO NR=2,NRMAX+1
          DEFF(NR,NSA)=-GAMMAP(NR,NSA)/DNDR(NR,NSA)
       END DO
    END DO
    CALL PAGES
       CALL GRD1D(1,RM,DNDT,  NRMAX,  NRMAX,  NSAMAX,'@dn/dt vs rho@')
       CALL GRD1D(2,RG,GAMMAP,NRMAX+1,NRMAX+1,NSAMAX,'@Gammap vs rho@')
       CALL GRD1D(3,RG,DNDR,  NRMAX+1,NRMAX+1,NSAMAX,'@dndr vs rho@')
       CALL GRD1D(4,RG,DEFF,  NRMAX+1,NRMAX+1,NSAMAX,'@Deff vs rho@')
    CALL PAGEE
    DEALLOCATE(DNDT,GAMMAP,GAMMAPR,DEFF,DNDR,PST)
  END SUBROUTINE fp_caldeff
END MODULE fpcaldeff
