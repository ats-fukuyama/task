MODULE trexec

! This module solves the matrix equation

  PUBLIC tr_exec
  PRIVATE

CONTAINS

  SUBROUTINE tr_exec

    USE trcomm,ONLY: &
         rkind,ikind,pi,rr,rkap,neqmax,neqrmax,nvmax,nvrmax,nrmax,dt,dr,rg, &
         DFa,Vca,Exa,Pha,elmtx,limtx,rpimtx,r1imtx,r2imtx,r3imtx,rimtx,lhmtx, &
         rhv,xv,xv_new,xv_prev,neqr_neq,id_neq,id_neqnr
    IMPLICIT NONE
    REAL(rkind) :: dh0,dh1,dh2,dh3,dvdrp,dvdrm,dvdr0
    INTEGER(ikind) :: nr,neq,neq1,neqr,neqr1,nvrm,nvrp,ierr,nvm,nvp

    LHMTX(1:4*neqrmax-1,1:nvrmax) = 0.D0
    RHV(1:nvrmax)                 = 0.D0

    DO NR = 1, NRMAX

!       DVDRP = 2.D0*PI*RR*2.D0*PI*RG(NR  )*RKAP
!       DVDRM = 2.D0*PI*RR*2.D0*PI*RG(NR-1)*RKAP
       DVDRP = RG(NR  )
       DVDRM = RG(NR-1)
       DVDR0 = 0.5D0*(DVDRM+DVDRP)

       DH0 = (rg(nr)-rg(nr-1))/12.D0
       DH1 = DT*DH0
       DH2 = DT*DVDR0/(rg(nr)-rg(nr-1))
       DH3 = 0.5D0*DT*DVDR0

       ELMTX(1:2*neqrmax,1:2*neqrmax)=0.D0

       DO NEQ = 1, NEQMAX

          LIMTX(1,1,NEQ) = DH0*(3.D0*DVDRM +      DVDRP)
          LIMTX(2,1,NEQ) = DH0*(     DVDRM +      DVDRP)
          LIMTX(1,2,NEQ) = DH0*(     DVDRM +      DVDRP)
          LIMTX(2,2,NEQ) = DH0*(     DVDRM + 3.D0*DVDRP) 

          RPIMTX(1,1,NEQ) = DH1*(3.D0*DVDRM +      DVDRP)
          RPIMTX(2,1,NEQ) = DH1*(     DVDRM +      DVDRP)
          RPIMTX(1,2,NEQ) = DH1*(     DVDRM +      DVDRP)
          RPIMTX(2,2,NEQ) = DH1*(     DVDRM + 3.D0*DVDRP) 

          DO NEQ1 = 1, NEQMAX
        
             R1IMTX(1,1,NEQ,NEQ1)=   DH2*DFa(NEQ,NEQ1,NR)
             R1IMTX(2,1,NEQ,NEQ1)= - DH2*DFa(NEQ,NEQ1,NR)
             R1IMTX(1,2,NEQ,NEQ1)= - DH2*DFa(NEQ,NEQ1,NR)
             R1IMTX(2,2,NEQ,NEQ1)=   DH2*DFa(NEQ,NEQ1,NR)
           
             R2IMTX(1,1,NEQ,NEQ1)= - DH3*VCa(NEQ,NEQ1,NR)
             R2IMTX(2,1,NEQ,NEQ1)= - DH3*VCa(NEQ,NEQ1,NR)
             R2IMTX(1,2,NEQ,NEQ1)=   DH3*VCa(NEQ,NEQ1,NR)
             R2IMTX(2,2,NEQ,NEQ1)=   DH3*VCa(NEQ,NEQ1,NR)
           
             R3IMTX(1,1,NEQ,NEQ1)= DH1*(3.D0*EXa(NEQ,NEQ1,NR-1)*DVDRM &
                                      +      EXa(NEQ,NEQ1,NR  )*DVDRP)
             R3IMTX(2,1,NEQ,NEQ1)= DH1*(     EXa(NEQ,NEQ1,NR-1)*DVDRM &
                                      +      EXa(NEQ,NEQ1,NR  )*DVDRP)
             R3IMTX(1,2,NEQ,NEQ1)= DH1*(     EXa(NEQ,NEQ1,NR-1)*DVDRM &
                                      +      EXa(NEQ,NEQ1,NR  )*DVDRP)
             R3IMTX(2,2,NEQ,NEQ1)= DH1*(     EXa(NEQ,NEQ1,NR-1)*DVDRM &
                                      + 3.D0*EXa(NEQ,NEQ1,NR  )*DVDRP)

             RIMTX(1,1,NEQ,NEQ1) &
                =R1IMTX(1,1,NEQ,NEQ1)+R2IMTX(1,1,NEQ,NEQ1)-R3IMTX(1,1,NEQ,NEQ1)
             RIMTX(2,1,NEQ,NEQ1) &
                =R1IMTX(2,1,NEQ,NEQ1)+R2IMTX(2,1,NEQ,NEQ1)-R3IMTX(2,1,NEQ,NEQ1)
             RIMTX(1,2,NEQ,NEQ1) &
                =R1IMTX(1,2,NEQ,NEQ1)+R2IMTX(1,2,NEQ,NEQ1)-R3IMTX(1,2,NEQ,NEQ1)
             RIMTX(2,2,NEQ,NEQ1) &
                =R1IMTX(2,2,NEQ,NEQ1)+R2IMTX(2,2,NEQ,NEQ1)-R3IMTX(2,2,NEQ,NEQ1)
          END DO
        
          !--- fixed value ---

          IF(id_neqnr(neq,nr-1) == 2)THEN
             LIMTX(1,1,NEQ) = 1.D0
             LIMTX(1,2,NEQ) = 0.D0

             RPIMTX(1,1,NEQ) = 0.D0
             RPIMTX(1,2,NEQ) = 0.D0

             DO NEQ1=1,NEQMAX
                RIMTX(1,1,NEQ,NEQ1) = 0.D0
                RIMTX(1,2,NEQ,NEQ1) = 0.D0
             END DO
          END IF

          IF(id_neqnr(neq,nr) == 2)THEN
             LIMTX(2,1,NEQ) = 0.D0
             LIMTX(2,2,NEQ) = 1.D0

             RPIMTX(2,1,NEQ) = 0.D0
             RPIMTX(2,2,NEQ) = 0.D0

             DO NEQ1=1,NEQMAX
                RIMTX(2,1,NEQ,NEQ1) = 0.D0
                RIMTX(2,2,NEQ,NEQ1) = 0.D0
             END DO
          END IF

          !--- calculate ELMTX with reduction ---

          IF(id_neq(neq) == 1) THEN
             neqr=neqr_neq(neq)
             DO NEQ1=1,NEQMAX
                IF(id_neq(neq1) == 1) THEN
                   neqr1=neqr_neq(neq1)
                   ELMTX(NEQR        ,NEQR1        ) = RIMTX(1,1,NEQ,NEQ1)
                   ELMTX(NEQR+NEQRMAX,NEQR1        ) = RIMTX(2,1,NEQ,NEQ1)
                   ELMTX(NEQR        ,NEQR1+NEQRMAX) = RIMTX(1,2,NEQ,NEQ1)
                   ELMTX(NEQR+NEQRMAX,NEQR1+NEQRMAX) = RIMTX(2,2,NEQ,NEQ1)
                END IF
             END DO
             ELMTX(NEQR        ,NEQR        ) &
           = ELMTX(NEQR        ,NEQR        )&
                  + LIMTX(1,1,NEQ)
             ELMTX(NEQR+NEQRMAX,NEQR        ) &
           = ELMTX(NEQR+NEQRMAX,NEQR        )&
                  + LIMTX(2,1,NEQ)
             ELMTX(NEQR        ,NEQR+NEQRMAX) &
           = ELMTX(NEQR        ,NEQR+NEQRMAX)&
                  + LIMTX(1,2,NEQ)
             ELMTX(NEQR+NEQRMAX,NEQR+NEQRMAX) &
           = ELMTX(NEQR+NEQRMAX,NEQR+NEQRMAX)&
                  + LIMTX(2,2,NEQ)
          END IF
       END DO

       !=== BUILD THE LEFT HAND BAND MATRIX         

       DO NEQR = 1, 2*NEQRMAX
          DO NEQR1 = 1, 2*NEQRMAX
             LHMTX(2*NEQRMAX+NEQR-NEQR1,(NR-1)*NEQRMAX+NEQR1) &
           = LHMTX(2*NEQRMAX+NEQR-NEQR1,(NR-1)*NEQRMAX+NEQR1) &
                  + ELMTX(NEQR1,NEQR)
          END DO
       END DO

       !=== BUILD THE RHS VECTOR

       NVM=(NR-1)*NEQMAX
       NVP= NR   *NEQMAX
       NVRM=(NR-1)*NEQRMAX
       NVRP= NR   *NEQRMAX
       DO NEQ = 1, NEQMAX
          IF(id_neq(neq) == 1) THEN
             neqr=neqr_neq(neq)
             RHV(NVRM+NEQR) = RHV(NVRM+NEQR)   &
                  + LIMTX(1,1,NEQ) * XV_PREV(NVM+NEQ) &
                  + LIMTX(1,2,NEQ) * XV_PREV(NVP+NEQ) &
                  +RPIMTX(1,1,NEQ) * PHa(NEQ,NR-1) &
                  +RPIMTX(1,2,NEQ) * PHa(NEQ,NR)

             RHV(NVRP+NEQR) = RHV(NVRP+NEQR)   &
                  + LIMTX(2,1,NEQ) * XV_PREV(NVM+NEQ)    &
                  + LIMTX(2,2,NEQ) * XV_PREV(NVP+NEQ)    &
                  +RPIMTX(2,1,NEQ) * PHa(NEQ,NR-1)   &
                  +RPIMTX(2,2,NEQ) * PHa(NEQ,NR)
          END IF
       END DO

    END DO ! NR LOOP END

    CALL BANDRD(lhmtx,rhv,nvrmax,4*neqrmax-1,4*neqrmax-1,ierr)
    IF(ierr /= 0) THEN
       WRITE(6,*) 'XX trexec: BANDRD error =',ierr
       STOP
    ENDIF

    DO NEQ=1,NEQMAX
       IF(id_neq(neq) == 1) THEN
          neqr=neqr_neq(neq)
          DO NR=0,NRMAX
             xv_new(NR*NEQMAX+NEQ) = rhv(NR*NEQRMAX+NEQR)
          END DO
       ELSE
          DO NR=0,NRMAX
             xv_new(NR*NEQMAX+NEQ) = xv(NR*NEQMAX+NEQ)
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE tr_exec
END MODULE trexec
