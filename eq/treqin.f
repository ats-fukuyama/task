C     $Id$
C
C   *******************************************
C   **    Initialize EQ interface for TR     **
C   *******************************************
C
      SUBROUTINE TREQIN
C
      CALL EQINIT
      CALL EQPARF
      RETURN
      END
C
C   *******************************************
C   **    EQ interface for TR     **
C   *******************************************
C
C     input:
C
C     NTRMAX        : Maximum array number
C     PSI(NTRMAX)   : Poloidal flux                  (m^2kgs^-2A^-2)
C     PPSI(NTRMAX)  : Pressure                                 (MPa)
C     HJPSI(NTRMAX) : Plasma current density                (MA/m^2)
C     VTPSI(NTRMAX) : Toroidal rotation velocity               (m/s)
C     TPSI(NTRMAX)  : Temperature                              (keV)
C
C     output:
C
C     QPSI(NTRMAX)  : Safety factor
C     VPSI(NTRMAX)  : Plamsa volume
C     SPSI(NTRMAX)  : Plamsa poloidal crosssection area
C
C   ***************************************************************
C
      SUBROUTINE TREQEX(PSIP,PPSI,HJPSI,VTPSI,TPSI,NTRMAX1,
     &                  QPSI,VPSI,SPSI,IERR)
C              
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
      DIMENSION PSIP(NTRMAX1),PPSI(NTRMAX1),HJPSI(NTRMAX1)
      DIMENSION VTPSI(NTRMAX1),TPSI(NTRMAX1)
      DIMENSION QPSI(NTRM),VPSI(NTRM),SPSI(NTRM)
      DIMENSION DERIV(NTRM)
C
      IERR=0
      NTRMAX = NTRMAX1
C
      DO NTR=1,NTRMAX
         PSITR(NTR)=PSIP(NTR)
      ENDDO
C
C     *** Calculate SPLINE coefficients ***
C
      CALL SPL1D(PSITR,PPSI,DERIV,UPPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D PPSI: IERR=',IERR
      CALL SPL1D(PSITR,HJPSI,DERIV,UJPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D HJPSI : IERR=',IERR
      CALL SPL1D(PSITR,VTPSI,DERIV,UVTPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D VTPSI : IERR=',IERR
      CALL SPL1D(PSITR,TPSI,DERIV,UTPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D TPSI : IERR=',IERR
C
      CALL EQMESH
      CALL EQPSIN
      CALL EQLOOP(IERR)
         IF(IERR.NE.0) GOTO 9000
      CALL EQTORZ
      CALL EQSETP
C
      NRMAX1=NRMAX
      NTHMAX1=NTHMAX
      NSUMAX1=NSUMAX
      CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
C
      DO NTR=1,NTRMAX
         PSIN=PSITR(NTR)/PSITR(NTRMAX)
         QPSI(NTR)=FNQPS(PSIN)
         VPSI(NTR)=FNVPS(PSIN)
         SPSI(NTR)=FNSPS(PSIN)
      ENDDO
 9000 RETURN
      END
