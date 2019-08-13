!     ***********************************************************

!           MAIN ROUTINE FOR TRANSPORT CALCULATION

!     ***********************************************************

      SUBROUTINE TRLOOP

      USE TRCOMM, ONLY : &
           DT, MODELG, NT, NTEQIT, NTMAX, Q0, QP, RG, RIP, RIPE, RIPS, &
           TPRST, TST, DIPDT
      USE TRCOM1, ONLY : NTAMAX
      use tr_bpsd,only: tr_bpsd_put, tr_bpsd_get
      use trunit
      use equnit_mod
      IMPLICIT NONE
      INTEGER(4):: IERR
      REAL(8)   :: FCTR

      IF(NT.GE.NTMAX) GOTO 9000

      CALL TREVAL(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

      RIP=RIPS
      IF(NTMAX.NE.0) DIPDT=(RIPE-RIPS)/(DBLE(NTMAX)*DT)
      write(6,'(A,1P4E12.4)') "**RIP,RIPS,RIPE,DIP=",RIP,RIPS,RIPE,DIPDT

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
!            if(modelg.eq.8) THEN
!               call equ_calc
!               call tr_bpsd_get(IERR)
!               if(ierr.ne.0) return
!            endif
            if(modelg.eq.9) THEN
               call tr_bpsd_put(IERR)
               call eq_calc
               call tr_bpsd_get(IERR)
               if(ierr.ne.0) return
            endif
         ENDIF
      ENDIF

      CALL TREVAL(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

!     ***

      IF(NT.LT.NTMAX) GOTO 1000

 9000 RIPS=RIPE
      RETURN
      END SUBROUTINE TRLOOP
