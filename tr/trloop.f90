! trloop.f90

MODULE trloop

  PRIVATE
  PUBLIC tr_loop

CONTAINS

!     ***********************************************************

!           MAIN ROUTINE FOR TRANSPORT CALCULATION

!     ***********************************************************

  SUBROUTINE tr_loop

      USE TRCOMM
      USE TRCOM1, ONLY : NTAMAX
      USE trbpsd, ONLY: tr_bpsd_put, tr_bpsd_get
      USE trexec
      USE equnit_mod
      IMPLICIT NONE
      INTEGER(4):: IERR,nr
      REAL(8)   :: FCTR

      IF(NT.GE.NTMAX) GOTO 9000

      CALL tr_eval(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
         DIPDT=0.D0
      ELSE
         RIP=RIPS
         IF(NTMAX.NE.0) DIPDT=(RIPE-RIPS)/(DBLE(NTMAX)*DT)
         write(6,'(A,1P4E12.4)') "**RIP,RIPS,RIPE,DIP=",RIP,RIPS,RIPE,DIPDT
      ENDIF

      call tr_bpsd_get(ierr)
      if(ierr.ne.0) GOTO 9000

 1000 CONTINUE

      CALL tr_exec(DT,IERR)
      IF(IERR.NE.0) GOTO 9000

      DO nr=1,nrmax
         QPINV(nr)=(4.D0*PI**2*RDPVRHOG(nr))/(TTRHOG(nr)*ARRHOG(nr))
      END DO
      call tr_bpsd_put(IERR)
      if(ierr.ne.0) GOTO 9000
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
!               call equ_calc
            endif
            IF(modelg.eq.9) THEN
               call eq_calc
            endif
         ENDIF
      ENDIF
      call tr_bpsd_get(IERR)
      if(ierr.ne.0) return

      CALL tr_eval(NT,IERR)
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
    END SUBROUTINE tr_loop
  END MODULE trloop
