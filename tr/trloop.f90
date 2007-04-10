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
      use trpl_mod,only: trpl_set, trpl_get
      use equnit_mod
      IMPLICIT NONE
      INTEGER(4):: I, ICHCK, IDGLOB, IERR, INFO, J, L, LDB, M, MWRMAX, N, NEQ, NEQ1, NEQRMAX, NF, NR, NRHS, NS, NSSN, NSSN1, &
     &             NSTN, NSTN1, NSVN, NSVN1, KL, KU
      REAL(8)   :: AJL, DIP, FACTOR0, FACTORM, FACTORP, FCTR, TSL


      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
         DIPDT=0.D0
      ELSE
         RIP=RIPS
         IF(NTMAX.NE.0) DIPDT=(RIPE-RIPS)/(DBLE(NTMAX)*DT)
      ENDIF

      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN
      CALL TRGLOB
      CALL TRSNAP
      IF(NGT.EQ.0) THEN
         CALL TRATOT
         CALL TRATOTN
      ENDIF
      IF(NGR.EQ.0) CALL TRATOG
      IF(NT.GE.NTMAX) GOTO 9000

 1000 CONTINUE

!      write(6,*) 'before trexec:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'before trexec:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      if(modelg.eq.9) then
         call trpl_get(ierr)
         if(ierr.ne.0) goto 9000
      endif

      CALL TREXEC(DT,IERR)
      IF(IERR.NE.0) GOTO 9000
      NT=NT+1

      call trpl_set(ierr)
      if(ierr.ne.0) goto 9000

!     /* Sawtooth Oscillation */
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))
      IF(Q0.LT.1.D0) TST=TST+DT

      IF(TST+0.5D0*DT.GT.TPRST) THEN
         CALL TRSAWT
         TST=0.D0
      ENDIF

!     *** DATA ACQUISITION FOR SHOWING GRAPH AND STATUS ***

      IDGLOB=0
      IF(MOD(NT,NTSTEP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRSNAP
      ENDIF
      IF(MOD(NT,NGTSTP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRATOT
         CALL TRATOTN
      ENDIF
      IF(MOD(NT,NGRSTP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRATOG
      ENDIF

!     *** SET GEOMETRY VIA TASK/EQ ***

      IF(NTEQIT.NE.0) THEN
         IF(MOD(NT,NTEQIT).EQ.0) THEN
            if(modelg.eq.9) THEN
               call eq_calc
               call trpl_get(IERR)
               if(ierr.ne.0) return
            endif
         ENDIF
      ENDIF

      IF(IDGLOB.EQ.0) CALL TRGLOB

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
