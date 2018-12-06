MODULE fpcaldeff

CONTAINS
  SUBROUTINE fp_caldeff
    USE fpcomm
    USE libgrf
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE::DNDT,GAMMAP,DEFF,DNDR
    REAL(rkind):: DTG2
    INTEGER:: NSA,NR,NT

    ALLOCATE(DNDT(NRMAX,NSAMAX))
!    DO NSA=1,NSAMAX
!      DO NR=1,NRMAX
!        DO NT=2,NTMAX-1
!          DNDT(NR,NT,NSA)=(RNSL(NR,NT+1,NSA)-RNSL(NR,NT-1,NSA))/(2.D0*DELT)
!        END DO
!        DNDT(NR,1,NSA)=(RNSL(NR,2,NSA)-RNSL(NR,1,NSA))/DELT
!        DNDT(NR,NTMAX,NSA)=(RNSL(NR,NTMAX,NSA)-RNSL(NR,NTMAX-1,NSA))/DELT
!      END DO
!    END DO
    IF(NT.EQ.NTG2) THEN
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
    ELSE IF(NT.GE.NTG2) THEN
       IF(NTG2.GT.0) THEN
          DTG2=PTG(NT)-PTG(NTG2)
          DO NSA=1,NSAMAX
             DO NR=1,NRMAX
                DNDT(NR,NSA)=(RNT(NR,NSA,NT)-RNT(NR,NSA,NTG2))/DTG2
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

    ALLOCATE(GAMMAP(NRMAX+1,NSAMAX))
    DO NSA=1,NSAMAX
       GAMMAP(1,NSA)=0.D0
       DO NR=1,NRMAX
          GAMMAP(NR+1,NSA)=GAMMAP(NR,NSA) &
               +(RSPBL(NR,NSA) &
                +RSPFL(NR,NSA) &
                +RSPSL(NR,NSA) &
                +RSPLL(NR,NSA) &
                +RSPSL_CX(NR,NSA) &
                -DNDT(NR,NSA))*(RG(NR+1)-RG(NR))
       END DO
    END DO
    ALLOCATE(DNDR(NRMAX+1,NSAMAX))
    DO NSA=1,NSAMAX
       DNDR(1,NSA)=0.D0
       DO NR=1,NRMAX
          DNDR(NR+1,NSA)=(RNSL(NR+1,NSA)-RNSL(NR,NSA))/DELR
       END DO
    END DO
    ALLOCATE(DEFF(NRMAX+1,NSAMAX))
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
    DEALLOCATE(DNDT,GAMMAP,DEFF,DNDR)
  END SUBROUTINE fp_caldeff
END MODULE fpcaldeff
