MODULE w1exec

CONTAINS

  SUBROUTINE w1_exec(ierr)

    USE w1comm
    USE w1prof,ONLY: w1_prof
    USE w1prep,ONLY: w1ants,w1fftl,w1setx,w1setz,w1syms
    USE w1bcnd,ONLY: w1_bcnd,w1evac
    USE w1disp,ONLY: w1dspa
    USE w1fem,ONLY: w1bnda,w1bndc,w1epwa,w1epwc
    USE w1mlm,ONLY: w1bndb,w1bndd,w1epwb,w1epwd,w1wkxb,w1wkxd
    USE w1intg,ONLY: w1bndq,w1epwq,w1qtbl,w1dspq
    USE w1rslt,ONLY: w1clcd,w1clpw,w1held,w1pwri,w1pwrs
    USE w1out,ONLY: w1file,w1prnt
    IMPLICIT NONE
    REAL(rkind):: DXD,RFSAVE,RKSAVE
    INTEGER:: NS,NL,NZP,IERR,IC,ICL,NX

    IF(MOD(NXVMAX,2).NE.0) NXVMAX=NXVMAX+1
    IF(ABS(RZ).LE.1.D-8) RZ=2.D0*PI*RR

    NHMAX=0
    DO NS=1,NSMAX
       IF(ABS(IHARM(NS)).GT.NHMAX)  NHMAX=ABS(IHARM(NS))
    END DO
    NCMAX=MAX(2*NHMAX+1,5)

    IF(NMODEL.EQ.6) THEN
       DXD=XDMAX/NDMAX
       CALL W1QTBL
    ENDIF
!     ******* 2-DIMENSIONAL ANALYSIS *******

    RFSAVE=RF
    RKSAVE=RKZ
    DO NL=1,NLOOP
       IF(NLOOP.NE.1) WRITE(6,630) RF,RKZ
       CALL W1SETZ
       CALL W1ANTS
       CALL W1SETX(IERR)
           IF(IERR.NE.0) GOTO 2000
       CALL W1_PROF
       CALL W1PWRI

!     ******* FOURIER TRANSFORM OF ANTENNA CURRENT *******

       CALL W1FFTL(CJ1,NZPMAX,0)
       CALL W1FFTL(CJ2,NZPMAX,0)
       CALL W1FFTL(CJ3,NZPMAX,0)
       CALL W1FFTL(CJ4,NZPMAX,0)

!     ******* CALCULATION FOR EACH KZ *******

       DO NZP = 1 , NZPMAX
          IF(NZP.LE.(NZPMAX/2+1).OR.NSYM.EQ.0) THEN
             RKZ   = AKZ(NZP)
             CFJY1 = CJ1(NZP)
             CFJY2 = CJ2(NZP)
             CFJZ1 = CJ3(NZP)
             CFJZ2 = CJ4(NZP)

             CALL W1_BCND
             SELECT CASE(NMODEL)
             CASE(0) ! FEM (no FLR)
                CALL W1DSPA
                CALL W1BNDA(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWA(NZP)
             CASE(1) ! MLM (no FLR)
                CALL W1DSPA
                CALL W1WKXB
                CALL W1BNDB(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWB(NZP)
                CALL W1HELD(4)
             CASE(2) ! FEM (fast wave)
                CALL W1DSPA
                CALL W1BNDA(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWA(NZP)
             CASE(3) ! MLM (fast wave)
                CALL W1DSPA
                CALL W1WKXB
                CALL W1BNDB(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWB(NZP)
                CALL W1HELD(4)
             CASE(4) ! FEM (2nd order FLR)
                CALL W1DSPA
                CALL W1BNDC(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWC(NZP)
             CASE(5) ! MLM  (2nd order FLR)
                CALL W1DSPA
                CALL W1WKXD
                CALL W1BNDD(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWD(NZP)
                CALL W1HELD(6)
             CASE(6) ! 
                CALL W1DSPQ
                CALL W1BNDQ(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWQ(NZP)
             END SELECT

             CALL W1EVAC(NZP)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP)
          ELSE
             CALL W1SYMS(NZP)
          END IF
       END DO

!     ******* INVERSE FOURIER TRANSFORM *******

       CALL W1FFTL(CJ1,NZPMAX,1)
       CALL W1FFTL(CJ2,NZPMAX,1)
       CALL W1FFTL(CJ3,NZPMAX,1)
       CALL W1FFTL(CJ4,NZPMAX,1)
!
       DO NX=1,NXTMAX
          DO IC=1,3
            CALL W1FFTL(CE2DA(1:NZPMAX,NX,IC),NZPMAX,1)
         END DO
      ENDDO

!     ******* POWER ABSORPTION AND OUTPUT *******

      CALL W1PWRS
      CALL W1PRNT
      CALL W1FILE
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
  630 FORMAT(1H / &
            1H ,'** RF = ',1PD12.4,' (MHZ) , RKZ = ',1PD12.4, &
                ' (/M) **')
  END SUBROUTINE w1_exec
END MODULE w1exec
