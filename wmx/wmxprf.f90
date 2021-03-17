! wmxprf.f90

MODULE wmxprf

  PRIVATE
  PUBLIC wm_xprf

CONTAINS
  
  SUBROUTINE wm_xprf ( IERR )

!=======================================================================

!     Read profile data for wm-code.

!     Input
!       wmcom1.inc
!         NRMAX  : number of the point for calculated
!         NTHMAX : number of poloidal mesh points
!         NSUMAX : number of the coordinates of separatrix surface,
!                  vacuum vessel
!         KNAMEQ : file name of eq data
!         ZEFF   : effective charge number
!         XRHO   : ro data of the calculated point ( sub. WMXRZF )

!       plcom1.inc
!         PNS(I) : density on plasma surface           (1.0E20/m3)
!         PTS(I) : temperature on plasma surface       (keV)
!         I=1:electron, I=2:Impurity, I=3:ion

!     Output
!       plcom1.inc
!         PN(I)  : density at center                   (1.0E20/m3)
!         PTPR(I): Parallel temperature at center      (keV)
!         PTPP(I): Perpendicular temperature at center (keV)
!         I=1:electron, I=2:Impurity, I=3:ion

!       wmcom1.inc
!         PN60(N,I) : density at the calculated point     (1/m3)
!         PT60(N,I) : temperature at the calculated point (eV)
!         I=1:electron, I=2:Impurity, I=3:ion

!=======================================================================

    USE wmcomm
    USE plprof,ONLY: pl_wmxprf,wmspl_prof
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER,SAVE:: NRMAXSV=0,NTHMAXSV=0,NSUMAXSV=0,NSMAXSV=0
    REAL(rkind),SAVE:: ZEFFSV=0.D0
    INTEGER:: IP1,IP2,NS,NRMAX1,NR
    REAL(rkind):: VAL

    IERR=0

    IP1 = 1
    IF(NRMAXSV.EQ.NRMAX.AND. &
       NTHMAXSV.EQ.NTHMAX.AND. &
       NSUMAXSV.EQ.NSUMAX.AND. &
       KNAMEQ_SAVE.EQ.KNAMEQ) IP1 = 0

    IP2 = 1
    IF(NSMAXSV.EQ.NSMAX.AND. &
       DABS(ZEFFSV-ZEFF).LT.1.0D-6 ) THEN
       IP2 = 0
       DO NS=1,NSMAX
          IF (DABS(PTSSV(NS)-PTS(NS)).GT.1.D-6 ) IP2 = 1
          IF (DABS(PNSSV(NS)-PNS(NS)).GT.1.D-6 ) IP2 = 1
       ENDDO
    ENDIF
    IF (IP1+IP2.EQ.0) GO TO 9996

!      WRITE (6,*) '==========  WMXPRF START  =========='

    IERR = 9999

!----  Open profile data file and read
!----  PRFNE, PRFTE is data at the point divided equally by rho 
!        defined by toroidal magnetic flux

    CALL PL_WMXPRF(IERR) ! in pl/plprof.f

    NRMAX1 = NRMAX

!----  Set profile data at the point calculated in wm-code.

    DO NR=1,NRMAX1
       DO NS=1,NSMAX
          CALL WMSPL_PROF(XRHO(NR),NS,PN60(NR,NS),PT60(NR,NS))
       ENDDO
    ENDDO

!----  Modification for charge neutrality after spline interpolation

    DO NR=1,NRMAX1
       VAL=0.D0
       DO NS=2,NSMAX-1
          VAL=VAL+PZ(NS)*PN60(NR,NS)
       ENDDO
       PN60(NR,NSMAX)=(PN60(NR,1)-VAL)/PZ(NSMAX)
    ENDDO

    DO NS=1,NSMAX
       PN(NS)   = PN60(1,NS) * 1.D-20
       PTPR(NS) = PT60(1,NS) * 1.D-3
       PTPP(NS) = PT60(1,NS) * 1.D-3
    ENDDO

!----  Debug write

!$$$      WRITE(6,8020)
!$$$      DO N=1,NRMAX1
!$$$         WRITE(6,'(I3,1P7(1XE10.3))') N, XRHO(N), (PN60(N,I),I=1,NXSPC)
!$$$      ENDDO
!$$$      WRITE(6,8030)
!$$$      DO N=1,NRMAX1
!$$$         WRITE(6,'(I3,1P7(1XE10.3))') N, XRHO(N), (PT60(N,I),I=1,NXSPC)
!$$$      ENDDO
!8020 FORMAT(' N ',4X,'XRHO',7X,'PNE',8X,'PNI1', &
!                  7X,'PNI2',7X,'PNI3',7X,'PNI4')
!8030 FORMAT(' N ',4X,'XRHO',7X,'PTE',8X,'PTI1', &
!                  7X,'PTI2',7X,'PTI3',7X,'PTI4')

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

    IF(NSMAX.EQ.2) GOTO 9997

    IERR = 0

    GO TO 9995
9997 WRITE (6,*) '     *****  NO IMPURITY  *****      '
9996 WRITE (6,*) '======  WMXPRF ALREADY READ  ======'
    GO TO 9999
9995 WRITE (6,*) '==========  WMXPRF COMPLETED  =========='
9999 RETURN
  END SUBROUTINE wm_xprf
END MODULE wmxprf
