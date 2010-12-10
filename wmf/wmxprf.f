C     $Id$
      SUBROUTINE WMXPRF ( IERR )
C=======================================================================
C
C     Read density and temperature at the magnetic axis
C
C     Input
C       wmcom1.inc
C         NRMAX  : number of the point for calculated
C         NTHMAX : number of poloidal mesh points
C         NSUMAX : number of the coordinates of separatrix surface,
C                  vacuum vessel
C         KNAMEQ : file name of eq data
C         ZEFF   : effective charge number
C         XRHO   : ro data of the calculated point ( sub. WMXRZF )
C
C       plcom1.inc
C         PNS(I) : density on plasma surface           (1.0E20/m3)
C         PTS(I) : temperature on plasma surface       (keV)
C         I=1:electron, I=2:Impurity, I=3:ion
C
C     Output
C       plcom1.inc
C         PN(I)  : density at center                   (1.0E20/m3)
C         PTPR(I): Parallel temperature at center      (keV)
C         PTPP(I): Perpendicular temperature at center (keV)
C         I=1:electron, I=2:Impurity, I=3:ion
C
C=======================================================================
C
      INCLUDE 'wmcomm.inc'
C
      REAL*8 ZEFFSV, PTSSV(NSM), PNSSV(NSM), PNCTR(NSM), PTCTR(NSM)
C
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV
      SAVE ZEFFSV, NSMAXSV, PTSSV, PNSSV
C
      DATA NRMAXSV,NTHMAXSV,NSUMAXSV/0,0,0/
      DATA ZEFFSV / 0.0D0 /
      DATA NSMAXSV /0/
      DATA PTSSV / NSM*0.0D0 /
      DATA PNSSV / NSM*0.0D0 /
C
      IERR=0
C
      IP1 = 1
      IF(NRMAXSV.EQ.NRMAX.AND.
     &   NTHMAXSV.EQ.NTHMAX.AND.
     &   NSUMAXSV.EQ.NSUMAX.AND.
     &   KNAMEQ_SAVE.EQ.KNAMEQ) IP1 = 0
C
      IP2 = 1
      IF(NSMAXSV.EQ.NSMAX.AND.
     &   DABS(ZEFFSV-ZEFF).LT.1.0D-6 ) THEN
         IP2 = 0
         DO NS=1,NSMAX
            IF (DABS(PTSSV(NS)-PTS(NS)).GT.1.D-6 ) IP2 = 1
            IF (DABS(PNSSV(NS)-PNS(NS)).GT.1.D-6 ) IP2 = 1
         ENDDO
      ENDIF
      IF (IP1+IP2.EQ.0) GO TO 9996
C
C      WRITE (6,*) '==========  WMXPRF START  =========='
C
      IERR = 9999
C
C----  Open profile data file and read
C----  PRFNE, PRFTE is data at the point divided equally by rho 
C        defined by toroidal magnetic flux
C
      CALL PLWMXPRF(IERR) ! in pl/plprof.f
C
C----  Read density and temperature at the magnetic axis
C
      NR = 1
      DO NS=1,NSMAX
         CALL WMSPL_PROF(XRHO(NR),NS,PNCTR(NS),PTCTR(NS))
      ENDDO
C
C----  Modification for charge neutrality
C
      VAL=0.D0
      DO NS=2,NSMAX-1
         VAL=VAL+PZ(NS)*PNCTR(NS)
      ENDDO
      PNCTR(NSMAX)=(PNCTR(1)-VAL)/PZ(NSMAX)
C
      DO NS=1,NSMAX
         PN(NS)   = PNCTR(NS) * 1.D-20
         PTPR(NS) = PTCTR(NS) * 1.D-3
         PTPP(NS) = PTCTR(NS) * 1.D-3
      ENDDO
C
      NRMAXSV=NRMAX
      NTHMAXSV=NTHMAX
      NSUMAXSV=NSUMAX
      KNAMEQ_SAVE=KNAMEQ
      ZEFFSV = ZEFF
      NSMAXSV = NSMAX
      DO NS=1,NSMAX
         PTSSV(NS)=PTS(NS)
         PNSSV(NS)=PNS(NS)
      ENDDO
C
      IF(NSMAX.EQ.2) GOTO 9997
C
      IERR = 0
C
      GO TO 9995
 9997 WRITE (6,*) '     *****  NO IMPURITY  *****      '
 9996 WRITE (6,*) '======  WMXPRF ALREADY READ  ======'
      GO TO 9999
 9995 WRITE (6,*) '==========  WMXPRF COMPLETED  =========='
 9999 RETURN
      END
