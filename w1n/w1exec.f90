MODULE w1exec

CONTAINS

  SUBROUTINE w1_exec(ierr)

    USE w1comm
    IMPLICIT NONE
    REAL(rkind):: DXD,RFSAVE,RKSAVE
    INTEGER:: NS,NL,NZP,IERR,IC,ICL,NX

    IF(MOD(NXVMAX,2).NE.0) NXVMAX=NXVMAX+1
    IF(ABS(RZ).LE.1.D-8) RZ=2.D0*PI*RR
    IF(NDMAX.GT.NDM-1) NDMAX=NDM-1
    NHARM=0
    IF(NMODEL.EQ.6) THEN
       DO NS=1,NSMAX
          IF(ABS(IHARM(NS)).GT.NHARMM) IHARM(NS)=NHARMM
          IF(ABS(IHARM(NS)).GT.NHARM)  NHARM=IHARM(NS)
       END DO
       DXD=XDMAX/NDMAX
       CALL W1QTBL(NDMAX,XDMAX,NHARM)
    ENDIF

!     ******* 2-DIMENSIONAL ANALYSIS *******

    RFSAVE=RF
    RKSAVE=RKZ
    DO NL=1,NLOOP
       IF(NLOOP.NE.1) WRITE(6,630) RF,RKZ
       CALL W1SETZ(IERR)
           IF(IERR.NE.0) GOTO 2000
       CALL W1ANTS
       CALL W1SETX(IERR)
           IF(IERR.NE.0) GOTO 2000
       CALL W1PROF
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
             CALL W1BCND
             IF(NMODEL.LE.5) THEN
                CALL W1DSPA(NALPHA)
             ELSE
                CALL W1DSPQ(ICL)
             ENDIF
             IF(NMODEL.EQ.0.OR.NMODEL.EQ.2) THEN
                CALL W1BNDA(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWA(NZP)
             ELSEIF(NMODEL.EQ.1.OR.NMODEL.EQ.3) THEN
                CALL W1WKXB
                CALL W1BNDB(IERR,NXABS)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWB(NZP)
                CALL W1HELD(4)
             ELSEIF(NMODEL.EQ.4) THEN
                CALL W1BNDC(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWC(NZP)
             ELSEIF(NMODEL.EQ.5) THEN
                CALL W1WKXD
                CALL W1BNDD(IERR,NXABS)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWD(NZP)
                CALL W1HELD(6)
             ELSEIF(NMODEL.EQ.6) THEN
                CALL W1BNDQ(IERR)
                   IF(IERR.NE.0) GOTO 2000
                CALL W1EPWQ(NZP)
             ENDIF
             CALL W1EVAC(NZP,NSYM)
             CALL W1CLCD(NZP)
             CALL W1CLPW(NZP,NSYM)
          ELSE
             CALL W1SYMS(NZP,NSYM)
          END IF
       END DO

       IF(NMODEL.EQ.6) WRITE(6,616) MATL,MATLM,ICL,NCLM,NDMAX,XDMAX

!     ******* INVERSE FOURIER TRANSFORM *******

       CALL W1FFTL(CJ1,NZPMAX,1)
       CALL W1FFTL(CJ2,NZPMAX,1)
       CALL W1FFTL(CJ3,NZPMAX,1)
       CALL W1FFTL(CJ4,NZPMAX,1)
!
       DO NX=1,NXTMAX
          DO IC=1,3
            CALL W1FFTL(CE2DA(1,NX,IC),NZP,1)
         END DO
      ENDDO

!     ******* POWER ABSORPTION AND OUTPUT *******

      CALL W1PWRS
      CALL W1PRNT(NPRINT)
      CALL W1FILE(NFILE)
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
                 'MATL  MATLM  ICL   NCLM  NDMAX XDMAX'/ &
             22X,5I6,1P1D12.4)
  630 FORMAT(1H / &
            1H ,'** RF = ',1PD12.4,' (MHZ) , RKZ = ',1PD12.4, &
                ' (/M) **')
  END SUBROUTINE w1_exec
END MODULE w1exec
