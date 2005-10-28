C     $Id$
      SUBROUTINE WMXPRF (
     O                    IERR )
C=======================================================================
C
C     This subroutine get profile data for wm-code.
C     Profile data got from tokrd file on nanasvr is for JT-60.
C
C     Input
C       wmcomm.inc
C         NRMAX  : number of the point for calculated
C         NTHMAX : number of poloidal mesh points
C         NSUMAX : number of the coordinates of separatrix surface,
C                  vacuum vessel
C         KNAMEQ : file name of eq data
C         ZEFF   : effective charge number
C         XRHO   : ro data of the calculated point ( sub. WMXRZF )
C         PSIP   : psi data of the calculated point ( sub. EQGETP )
C
C     Output
C       wmcomm.inc
C         PN(I)  : density at center                   (1.0E20/m*3)
C         PNS(I) : density on plasma surface           (1.0E20/m*3)
C         PTPR(I): Parallel temperature at center      (keV)
C         PTPP(I): Perpendicular temperature at center (keV)
C         PTS(I) : temperature on plasma surface       (keV)
C         I=1:electron, I=2:Impurity, I=3:ion
C       wmxprf.inc
C         NPRFI     : calculation mode for ion
C         PN60(N,I) : density at the calculated point     (1/m3)
C         PT60(N,I) : temperature at the calculated point (eV)
C         I=1:electron, I=2:Impurity, I=3:ion
C
C                                     created by Itakura : CSK  00/07/31
C=======================================================================
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'wmxprf.inc'
C
      REAL*8 PRFPSI(NXPRF)  , PRFNE(NXPRF)   , PRFTE(NXPRF)
     >     , PRFTI(NXPRF)
      REAL*8 UPRFNE(4,NXPRF), UPRFTE(4,NXPRF), UPRFTI(4,NXPRF)
      REAL*8 DERIV(NXPRF)
      REAL*8 ZEFFSV, PSIL, PPL, PTSSV(NSM), PNSSV(NSM)
      REAL*8 XPAR1, XPAR2, XPAR
C
      CHARACTER    TRFILE*80, CWK1*10
C
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV
      SAVE ZEFFSV, NSMAXSV, PTSSV, PNSSV
C
      DATA TRFILE / 'trdata' /  ! fixed name
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
      IF (IP1+IP2.EQ.0) GO TO 9997
C
      WRITE (6,*) '==========  WMXPRF START  =========='
C
      IERR = 9999
      NRMAX1 = NRMAX + 1
      XPAR1 = 0.0D0
      XPAR2 = 0.0D0
      DO NS=3,NSMAX
         XPAR1 = XPAR1 + PZ(NS)*PZ(NS)*PNS(NS)
         XPAR2 = XPAR2 + PZ(NS)*PNS(NS)
      ENDDO
      IF (XPAR2.EQ.0.0D0) GO TO 9997
      XPAR = XPAR1/XPAR2
C
C----  Set ion calculation mode
C
      IF (IP2.EQ.1) THEN
        NPRFI = 0
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
        DO WHILE ( NPRFI.NE.1.AND.NPRFI.NE.2 )
          WRITE(6,*) 'SELECT MODE OF CALCULATION FOR ION'
          WRITE(6,*) '  1 : USE PROFILE DATA FORM FILE'
          WRITE(6,*) '  2 : USE DEFAULT METHOD OF WM CODE'
          WRITE(6,*) 'SELECT 1 (DEFAULT) or 2 >>'
          CWK1 = ' '
          READ(5,'(A)',ERR=210) CWK1
          IF (CWK1.EQ.' ') THEN
            NPRFI = 1
          ELSE
            READ(CWK1,*,ERR=210) NPRFI
          ENDIF
       ENDDO
  210  CONTINUE
      ENDIF
C
C----  Open profile data file and read
C----  PRFNE, PRFTE is data at the point divided equally by psi
C
C      CALL HANDLE( IFNO )
C      IF ( IFNO.EQ.0 ) GO TO 9999
C
      IFNO=22
      OPEN ( IFNO, FILE=TRFILE, ERR=9998 )
      READ ( IFNO, '(I3)', END=9999, ERR=9999 ) NPRF
      DO N=1,NPRF
        READ ( IFNO, '(4E14.7)', END=9999, ERR=9999 )
     >                          PRFPSI(N), PRFNE(N), PRFTE(N), PRFTI(N)
      ENDDO
C
C----  Set coefficient for spline
C
      CALL SPL1D(PRFPSI,PRFNE,  DERIV,UPRFNE, NPRF,0,IRC)
      IF (IRC.NE.0) GO TO 9999
      CALL SPL1D(PRFPSI,PRFTE,  DERIV,UPRFTE, NPRF,0,IRC)
      IF (IRC.NE.0) GO TO 9999
      CALL SPL1D(PRFPSI,PRFTI,  DERIV,UPRFTI, NPRF,0,IRC)
      IF (IRC.NE.0) GO TO 9999
C
C----  Set profile data at the point calculated in wm-code.
C----    PSIP  : psi value at the point ( calculated at subroutine eqpsic )
C
      DO NR=1,NRMAX1
        IF (XRHO(NR).GT.1.0D0) THEN
          PN60(NR,1) = 0.0D0
          PT60(NR,1) = PRFTE(NPRF)*PTS(1)
          DO NS=2,NSMAX
             PN60(NR,NS) = 0.0D0
             PT60(NR,NS) = PRFTI(NPRF)*PTS(NS)
          ENDDO
        ELSE
          PSIL=PSIP(NR)
          CALL SPL1DF(PSIL,PPL,PRFPSI,UPRFNE,NPRF,IRC)
          PN60(NR,1)=PPL*PNS(1)
          CALL SPL1DF(PSIL,PPL,PRFPSI,UPRFTE,NPRF,IRC)
          PT60(NR,1)=PPL*PTS(1)
          CALL SPL1DF(PSIL,PPL,PRFPSI,UPRFTI,NPRF,IRC)
          PT60(NR,2)=PPL*PTS(2)
C----  not need PNS(2)
          PN60(NR,2)=(1.D0-(PZ(2)-ZEFF)/(PZ(2)-XPAR))/PZ(2)*PN60(NR,1)
          DO NS=3,NSMAX
             PN60(NR,NS)=(PZ(2)-ZEFF)/(PZ(2)-XPAR)/XPAR2*PN60(NR,1)
     &                  *PNS(NS)
             PT60(NR,NS)=PPL*PTS(NS)
          ENDDO
        ENDIF
      ENDDO
C     
C----  Change the value at center and surface
C
      PN(1)  = PRFNE(1)*1.0D-20*PNS(1)
      PTPR(1) = PRFTE(1)*1.0D-3*PTS(1)
      PTPP(1) = PRFTE(1)*1.0D-3*PTS(1)
C
      IF (NPRFI.EQ.1) THEN
C----  not need PNS(2)
         PN(2)=(1.D0-(PZ(2)-ZEFF)/(PZ(2)-XPAR))/PZ(2)*PN(1)
         PTPR(2) = PRFTI(1)*1.0D-3*PTS(2)
         PTPP(2) = PRFTI(1)*1.0D-3*PTS(2)
C
         DO NS=3,NSMAX
            PN(NS)=(PZ(2)-ZEFF)/(PZ(2)-XPAR)/XPAR2*PN(1)*PNS(NS)
            PTPR(NS) = PRFTI(1)*1.D-3*PTS(NS)
            PTPP(NS) = PRFTI(1)*1.D-3*PTS(NS)
         ENDDO
      ENDIF
C
C----  Debug write
C
C     WRITE(6,8000)
C     DO 8010 N=1,NRMAX1
C       WRITE(6,'(I3,1P7E10.3)') N, XRHO(N), PSIP(N), PRFPSI(N)
C    >                                     , PRFNE(N), PNE(N)
C    >                                     , PRFTI(N), PTI(N)
C8010 CONTINUE
C8000 FORMAT(' N ',3X,'XRHO',6X,'PSIP ',4X,'PRFPSI',5X
C    >      ,'PRFNE',6X,'PNE',6X,'PRFTI',6X,'PTI')
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
      IERR = 0
C
 9999 CONTINUE
      CLOSE( IFNO )
 9998 CONTINUE
      WRITE (6,*) '==========  WMXPRF  END   =========='
 9997 CONTINUE
      RETURN
      END
