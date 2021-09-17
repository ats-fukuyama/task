MODULE fpcalchieff

CONTAINS
  SUBROUTINE fp_calchieff
    USE fpcomm
    USE fpsave
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE::rhot,rhot_prev,drhotdt,hf,hfr,dtdr,chieff!,ps
    INTEGER::nr,nsa

    ALLOCATE(rhot(nrmax,nsamax))
    ALLOCATE(rhot_prev(nrmax,nsamax))
    ALLOCATE(drhotdt(nrmax,nsamax))
    ! ALLOCATE(ps(nrmax,nsamax))
    ALLOCATE(hf(nrmax+1,nsamax))!heat flux
    ALLOCATE(hfr(nrmax+1,nsamax))!hf * minor radius
    ALLOCATE(dtdr(nrmax+1,nsamax))
    ALLOCATE(chieff(nrmax+1,nsamax))

    CALL MOMENT_0TH_ORDER(FNSP,RNSL)
    CALL MOMENT_2ND_ORDER(FNSP,rhot)
    CALL MOMENT_2ND_ORDER(FNSM,rhot_prev)

    DO nsa=1,nsamax
      DO nr=1,nrmax
        drhotdt(nr,nsa)=(rhot(nr,nsa)-rhot_prev(nr,nsa))/delt
      END DO
    END DO
    ! DO NTH=1,NTHMAX
    !   DO NP=1,NPMAX
    !    PS(NR,NSA)=PS(NR,NSA)&
    !                +(rspbl(NTH,NP,NR,NSA) &
    !                +rspfl(NTH,NP,NR,NSA) &
    !                +rspsl(NTH,NP,NR,NSA) &
    !                +rspll(NTH,NP,NR,NSA) &
    !                +rspsl_cx(NTH,NP,NR,NSA) &
    !                )*VOLP(NTH,NP,NSA)
    !   END DO
    ! END DO
    DO NSA=1,NSAMAX
      hfr(1,NSA)=0.D0
      hf(1,NSA)=0.D0
      DO NR=1,NRMAX
        hfr(NR+1,NSA)=hfr(NR,NSA)+RA*(NR/NRMAX)*(drhotdt(NR,NSA) &
                      -(RSPBL(NR,NSA) &
                      +RSPFL(NR,NSA) &
                      +RSPSL(NR,NSA) &
                      +RSPLL(NR,NSA) &
                      +RSPSL_CX(NR,NSA))) &
                      *(RG(NR+1)-RG(NR))
        hf(NR+1,NSA)=hf(NR+1,NSA)/RA*(NR/NRMAX)
      END DO
    END DO
    DO NSA=1,NSAMAX
      dtdr(1,nsa)=0.d0
      DO NR=1,NRMAX
        dtdr(nr+1,nsa)=((2.d0*rhot(nr+1,nsa))/(3.d0*rnsl(nr+1,nsa))-(2.d0*rhot(nr,nsa))/(3.d0*rnsl(nr,nsa)))/(rg(nr+1)-rg(nr))
      END DO
    END DO
    DO NSA=1,NSAMAX
      chieff(1,NSA)=0.D0
      DO NR=2,NRMAX+1
        chieff(NR,NSA)=-hf(NR,NSA)/(rnsl(nr,nsa)*dtdr(nr,nsa))
      END DO
    END DO

    DEALLOCATE(rhot,rhot_prev,drhotdt,hf,hfr,dtdr,chieff)

  END SUBROUTINE fp_calchieff
END MODULE fpcalchieff
