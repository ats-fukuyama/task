MODULE w1exec

CONTAINS

  SUBROUTINE w1_exec(ierr)

    USE w1comm
    USE w1prepost,ONLY: w1pre,w1post
    USE w1prepostx,ONLY: w1prex,w1postx
    USE w1prep,ONLY: w1syms
    USE w1bcnd,ONLY: w1_bcnd,w1evac
    USE w1disp,ONLY: w1dspa
    USE w1fem,ONLY: w1bnda,w1bndc,w1epwa,w1epwc
    USE w1mlm,ONLY: w1bndb,w1bndd,w1epwb,w1epwd,w1wkxb,w1wkxd
    USE w1intg,ONLY: w1bndq,w1epwq,w1qtbl,w1dspq
    USE w1rslt,ONLY: w1clcd,w1clpw,w1held
    USE w1exec6,ONLY: w1_exec6
    USE w1exec7,ONLY: w1_exec7
    USE w1exec8,ONLY: w1_exec8,w1qtblx
    USE w1exec9,ONLY: w1_exec9
    USE w1exec10,ONLY: w1_exec10
    USE w1exec11,ONLY: w1_exec11
    IMPLICIT NONE
    REAL(rkind):: DXD,RFSAVE,RKSAVE
    INTEGER:: NS,NL,NZP,IERR,IC,ICL,NX

    IERR=0

    WRITE(6,*) '@@@ point 20'
    IF(MOD(NXVMAX,2).NE.0) NXVMAX=NXVMAX+1
!    IF(ABS(RZ).LE.1.D-8) RZ=2.D0*PI*RR

    NHMAX=0
    DO NS=1,NSMAX
       IF(ABS(IHARM(NS)).GT.NHMAX)  NHMAX=ABS(IHARM(NS))
    END DO
    NCMAX=MAX(2*NHMAX+1,5)

    SELECT CASE(NMODEL)
    CASE(0:5)
       CALL w1pre(IERR)
    CASE(6:7,9:10)
       CALL w1prex(IERR)
    CASE(8,11)
       CALL w1prex(IERR)
       DXD=XDMAX/NDMAX
       CALL W1QTBLX
    CASE(12)
       CALL w1pre(IERR)
       DXD=XDMAX/NDMAX
       CALL W1QTBL
    END SELECT

    WRITE(6,*) '@@@ point 21'
!     ******* 2-DIMENSIONAL ANALYSIS *******

    RFSAVE=RF
    RKSAVE=RKZ
    DO NL=1,NLOOP

       IF(NLOOP.NE.1) WRITE(6,630) RF,RKZ
       WRITE(6,*) '@@@ point 22'

!     ******* CALCULATION FOR EACH KZ *******

       DO NZP = 1 , NZMAX
          RKZ   = AKZ(NZP)
          CFJY1 = CJ1(NZP)
          CFJY2 = CJ2(NZP)
          CFJZ1 = CJ3(NZP)
          CFJZ2 = CJ4(NZP)
          CFWG1 = CWG1(NZP)
          CFWG2 = CWG2(NZP)
          CFWG3 = CWG3(NZP)
          CFWG4 = CWG4(NZP)
!          WRITE(6,'(A,1p4E12.4)') 'CFJY=',CFJY1,CFJY2
!          WRITE(6,'(A,1p4E12.4)') 'CFJZ=',CFJZ2,CFJZ2

          WRITE(6,*) '@@@ point 23 NMODEL=',NMODEL
          SELECT CASE(NMODEL)
          CASE(0) ! FEM (no FLR)
             WRITE(6,*) '@@@ point 23'
             CALL W1_BCND
             CALL W1DSPA
             CALL W1BNDA(IERR)
                IF(IERR.NE.0) GOTO 2000
             CALL W1EPWA(NZP)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          CASE(1) ! MLM (no FLR)
             CALL W1_BCND
             CALL W1DSPA
             CALL W1WKXB
             CALL W1BNDB(IERR)
                IF(IERR.NE.0) GOTO 2000
             CALL W1EPWB(NZP)
             CALL W1HELD(4)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          CASE(2) ! FEM (fast wave)
             CALL W1_BCND
             CALL W1DSPA
             CALL W1BNDA(IERR)
                IF(IERR.NE.0) GOTO 2000
             CALL W1EPWA(NZP)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          CASE(3) ! MLM (fast wave)
             CALL W1_BCND
             CALL W1DSPA
             CALL W1WKXB
             CALL W1BNDB(IERR)
                IF(IERR.NE.0) GOTO 2000
             CALL W1EPWB(NZP)
             CALL W1HELD(4)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          CASE(4) ! FEM (2nd order FLR)
             CALL W1_BCND
             CALL W1DSPA
             CALL W1BNDC(IERR)
                IF(IERR.NE.0) GOTO 2000
             CALL W1EPWC(NZP)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          CASE(5) ! MLM  (2nd order FLR)
             WRITE(6,*) '@@@ point 240'
             CALL W1_BCND
             WRITE(6,*) '@@@ point 241'
             CALL W1DSPA
             WRITE(6,*) '@@@ point 242'
             CALL W1WKXD
             WRITE(6,*) '@@@ point 243'
             CALL W1BNDD(IERR)
                IF(IERR.NE.0) GOTO 2000
             CALL W1EPWD(NZP)
             CALL W1HELD(6)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          CASE(6) ! FEM cold collsional (one region, reflecting wall))
             CALL W1_EXEC6(NZP,IERR)
                IF(IERR.NE.0) GOTO 2000
          CASE(7) ! FEM warm differential (one region, reflecting wall))
             CALL W1_EXEC7(NZP,IERR)
                IF(IERR.NE.0) GOTO 2000
          CASE(8) ! FEM hot integro-differential (one region, reflecting wall)
             CALL W1_EXEC8(NZP,IERR)
                IF(IERR.NE.0) GOTO 2000
          CASE(9) ! FEM cold collsional (one region, absorbing wall))
             CALL W1_EXEC9(NZP,IERR)
                IF(IERR.NE.0) GOTO 2000
          CASE(10) ! FEM warm differential (one region, absorbing wall))
             CALL W1_EXEC10(NZP,IERR)
                IF(IERR.NE.0) GOTO 2000
          CASE(11) ! FEM hot integro-differential (one region, absorbing wall)
             CALL W1_EXEC11(NZP,IERR)
                IF(IERR.NE.0) GOTO 2000
          CASE(12) ! old FEM hot integro-differential
             CALL W1_BCND
             CALL W1DSPQ
             CALL W1BNDQ(IERR)
             IF(IERR.NE.0) GOTO 2000
             CALL W1EPWQ(NZP)
             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          END SELECT

       END DO

       SELECT CASE(NMODEL)
       CASE(0:5,12)
          CALL w1post
       CASE(6:11)
          CALL w1postx
       END SELECT

      IF(NLOOP.NE.1) THEN
         RF   = RF   + DRF
         RKZ  = RKZ  + DRKZ
      ENDIF
    END DO

    RF =RFSAVE
    RKZ=RKSAVE

    RETURN

2000 CONTINUE
    RETURN
  616 FORMAT('## INTEGRO-DIFF EQ. : ', &
                 'NLWMAX NCLMAX NCMAX   NDMAX   XDMAX'/ &
             22X,5I8,1P1D12.4)
  630 FORMAT('RF,RKZ=',2ES12.4)
  END SUBROUTINE w1_exec
END MODULE w1exec
