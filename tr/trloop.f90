!     ***********************************************************

!           MAIN ROUTINE FOR TRANSPORT CALCULATION

!     ***********************************************************

      SUBROUTINE TRLOOP

      USE TRCOMM, ONLY : &
     &   ABRHOG, AJ, AJU, AMM, ANC, ANFE, ANNU, AR1RHOG, AX, AY, AZ, BB, &
     &   DR, DT, DVRHO, DVRHOG, EPSLTR, LDAB, LMAXTR, MDLEQB, MDLEQE, &
     &   MDLEQN, MDLPCK, MDLUF, MDTC, MLM, MODELG, NEQMAX, NFM, NGR, NGRSTP, &
     &   NGT, NGTSTP, NRAMAX, NRM, NRMAX, NROMAX, NSM, NSMAX, NSS, NST, NSV, &
     &   NT, NTEQIT, NTMAX, NTSTEP, NTUM, &
     &   PA, PZ, PZC, PZFE, Q0, QP, RDP, RG, RHOA, RIP, RIPE, RIPS, RIPU, &
     &   RMU0, RN, RR, RT, RU, RW, T, TPRST, TST, TTRHO, TTRHOG, &
     &   VLOOP, VSEC, X, XV, Y, YV, Z, ZV ,NEQMAXM, DIPDT
      USE TRCOM1, ONLY : TMU, TMU1, NTAMAX, NTXMAX, NTXMAX1
      use tr_bpsd,only: tr_bpsd_set, tr_bpsd_get
      use trunit
      use equnit_mod
      use equunit_mod
      IMPLICIT NONE
      INTEGER(4):: I, ICHCK, IDGLOB, IERR, INFO, J, L, LDB, M, MWRMAX, N, NEQ, NEQ1, NEQRMAX, NF, NR, NRHS, NS, NSSN, NSSN1, &
     &             NSTN, NSTN1, NSVN, NSVN1, KL, KU
      REAL(8)   :: AJL, DIP, FACTOR0, FACTORM, FACTORP, FCTR, TSL

      IF(NT.GE.NTMAX) GOTO 9000

      CALL TREVAL(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

 1000 CONTINUE

      CALL tr_exec(DT,IERR)
      IF(IERR.NE.0) GOTO 9000
      NT=NT+1

!     /* Sawtooth Oscillation */
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))
      IF(Q0.LT.1.D0) TST=TST+DT

      IF(TST+0.5D0*DT.GT.TPRST) THEN
         CALL TRSAWT
         TST=0.D0
      ENDIF

!     *** SET GEOMETRY VIA TASK/EQ ***

      IF(NTEQIT.NE.0) THEN
         IF(MOD(NT,NTEQIT).EQ.0) THEN
            if(modelg.eq.8) THEN
               call equ_calc
               call tr_bpsd_get(IERR)
               if(ierr.ne.0) return
            endif
            if(modelg.eq.9) THEN
               call eq_calc
               call tr_bpsd_get(IERR)
               if(ierr.ne.0) return
            endif
         ENDIF
      ENDIF

      CALL TREVAL(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

!     *** READING DATA FROM UFILES FOR NEXT STEP ***

      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) CALL TR_UFREAD
!      IF(MDLUF.EQ.2.AND.MODEP.EQ.3) CALL TR_UFREAD_S
      IF(MDLUF.EQ.2) CALL TR_UFREAD_S

!     ***

      IF(NT.LT.NTMAX) GOTO 1000

 9000 IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         RIPS=RIP
         RIPE=RIP
      ELSE
         RIPS=RIPE
      ENDIF
      RETURN
      END SUBROUTINE TRLOOP
