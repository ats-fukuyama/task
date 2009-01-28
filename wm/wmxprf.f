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
C
C         PNS(I) : density on plasma surface           (1.0E20/m3)
C         PTS(I) : temperature on plasma surface       (keV)
C         I=1:electron, I=2:Impurity, I=3:ion
C
C     Output
C       wmcomm.inc
C         PN(I)  : density at center                   (1.0E20/m3)
C         PTPR(I): Parallel temperature at center      (keV)
C         PTPP(I): Perpendicular temperature at center (keV)
C         I=1:electron, I=2:Impurity, I=3:ion
C
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
      REAL*8 ZEFFSV, RHON, PPL, PTSSV(NSM), PNSSV(NSM)
C
      CHARACTER    TRFILE*80, CWK1*10
C
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV
      SAVE ZEFFSV, NSMAXSV, PTSSV, PNSSV
C
      DATA TRFILE / 'topics-data' /  ! fixed name
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
      WRITE (6,*) '==========  WMXPRF START  =========='
C
      IERR = 9999
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
 210     CONTINUE
      ENDIF
C
C----  Open profile data file and read
C----  PRFNE, PRFTE is data at the point divided equally by rho 
C        defined by toroidal magnetic flux
C
C      CALL HANDLE( IFNO )
C      IF ( IFNO.EQ.0 ) GO TO 9999
C
      IFNO=22
      OPEN ( IFNO, FILE=TRFILE, ERR=9998 )
      READ ( IFNO, '(I3)', END=9999, ERR=9999 ) NPRF
      DO N=1,NPRF
         READ ( IFNO, '(11E14.7)', END=9999, ERR=9999 )
     >        PRFRHO(N), (PRFN(N,I), I=1,NXSPC),
     >                   (PRFT(N,I), I=1,NXSPC)
      ENDDO
C
C----  Modification for charge neutrality
C
      DO NR=1,NPRF
         VAL=0.D0
         DO NS=2,NSMAX-1
            VAL=VAL+PZ(NS)*PRFN(NR,NS)
         ENDDO
         PRFN(NR,NSMAX)=(PRFN(NR,1)-VAL)/PZ(NSMAX)
      ENDDO
C
C----  Set coefficient for spline
C
      DO NS=1,NSMAX
         CALL SPL1D(PRFRHO,PRFN(1,NS),  DERIV,UPRFN(1,1,NS), NPRF,0,IRC)
         IF (IRC.NE.0) GO TO 9999
         CALL SPL1D(PRFRHO,PRFT(1,NS),  DERIV,UPRFT(1,1,NS), NPRF,0,IRC)
         IF (IRC.NE.0) GO TO 9999
      ENDDO
C
      IF(MDLWMF.EQ.0) THEN
         NRMAX1 = NRMAX + 1
      ELSE
         NRMAX1 = NRMAX
      ENDIF
C
C----  Set profile data at the point calculated in wm-code.
C----    RHON  : rho value at the point
C
      DO NR=1,NRMAX1
         RHON=XRHO(NR)
         DO NS=1,NSMAX
            CALL WMSPL_PROF(RHON,NS,PN60(NR,NS),PT60(NR,NS))
         ENDDO
      ENDDO
C
C----  Modification for charge neutrality after spline interpolation
C
      DO NR=1,NRMAX1
         VAL=0.D0
         DO NS=2,NSMAX-1
            VAL=VAL+PZ(NS)*PN60(NR,NS)
         ENDDO
         PN60(NR,NSMAX)=(PN60(NR,1)-VAL)/PZ(NSMAX)
      ENDDO
C
      IF (NPRFI.EQ.1) THEN
         DO NS=1,NSMAX
            PN(NS)   = PN60(1,NS) * 1.D-20
            PTPR(NS) = PT60(1,NS) * 1.D-3
            PTPP(NS) = PT60(1,NS) * 1.D-3
         ENDDO
      ENDIF
C
C----  Debug write
C
c$$$      WRITE(6,8000)
c$$$      DO N=1,NPRF
c$$$         WRITE(6,'(I3,1P6(1XE10.3))') N,PRFRHO(N),(PRFN(N,I),I=1,NXSPC)
c$$$      ENDDO
c$$$      WRITE(6,8010)
c$$$      DO N=1,NPRF
c$$$         WRITE(6,'(I3,1P6(1XE10.3))') N,PRFRHO(N),(PRFT(N,I),I=1,NXSPC)
c$$$      ENDDO
c$$$      WRITE(6,8020)
c$$$      DO N=1,NRMAX1
c$$$         WRITE(6,'(I3,1P7(1XE10.3))') N, XRHO(N), (PN60(N,I),I=1,NXSPC)
c$$$      ENDDO
c$$$      WRITE(6,8030)
c$$$      DO N=1,NRMAX1
c$$$         WRITE(6,'(I3,1P7(1XE10.3))') N, XRHO(N), (PT60(N,I),I=1,NXSPC)
c$$$      ENDDO
 8000 FORMAT(' N ',3X,'PRFRHO',6X,'PRFNE',6X,'PRFNI1',5X,'PRFNI2',
     >             5X,'PRFNI3',5X,'PRFNI4')
 8010 FORMAT(' N ',3X,'PRFRHO',6X,'PRFTE',6X,'PRFTI1',5X,'PRFTI2',
     >             5X,'PRFTI3',5X,'PRFTI4')
 8020 FORMAT(' N ',4X,'XRHO',7X,'PNE',8X,'PNI1',
     >             7X,'PNI2',7X,'PNI3',7X,'PNI4')
 8030 FORMAT(' N ',4X,'XRHO',7X,'PTE',8X,'PTI1',
     >             7X,'PTI2',7X,'PTI3',7X,'PTI4')
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
 9999 CONTINUE
      CLOSE( IFNO )
 9998 CONTINUE
      WRITE (6,*) '==========  WMXPRF  END   =========='
      GO TO 9995
 9997 WRITE (6,*) '     *****  NO IMPURITY  *****      '
 9996 WRITE (6,*) '======  WMXPRF  ABNORMAL END  ======'
 9995 RETURN
      END
C
C**************************************************
C
C     Interpolation of profile at a given point
C
C**************************************************
C
      SUBROUTINE WMSPL_PROF(Rhol,NS,PNL,PTL)
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'wmxprf.inc'
C
C---- Input
      integer NS
      real*8 Rhol ! Normalized radial mesh
C---- Output
      real*8 PNL, ! density
     &       PTL  ! temperature
C---- Internal
      real*8 PPL
C
C---- The following variables come from "wmxprf.inc".
C        NPRFI,NPRF,
C        PRFRHO,PRFNE,PRFTE,PRFTI,UPRFNE,UPRFTE,UPRFTI,DERIV
C
C---- Carry input parameter "ZEFF" to PLPLOF
      ZEFFWM = ZEFF
C
C----  Set profile data at the point calculated in wm-code.
C
      IF (Rhol.GT.1.0D0) THEN
         PNL = 0.D0
         PTL = PTS(NS)
      ELSE
         CALL SPL1DF(Rhol,PPL,PRFRHO,UPRFN(1,1,NS),NPRF,IRC)
         PNL=PPL
         CALL SPL1DF(Rhol,PPL,PRFRHO,UPRFT(1,1,NS),NPRF,IRC)
         PTL=PPL
      ENDIF
C
      RETURN
      END
