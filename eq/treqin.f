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
C     RR            : Major plasma radius (Geometrical center)   (m)
C     RA            : Minor plasma radius (Horizontal plane)     (m)
C     RB            : Minor wall radius                          (m)
C     RKAP          : Elongation
C     RDLT          : Triangularity
C     BB            : Toroidal magnetic field at R=RR            (T)
C     RIP           : Total plasma current                      (MA)
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
      SUBROUTINE TREQEX(RR1,RA1,RB1,RKAP1,RDLT1,BB1,RIP1,
     &                  NTRMAX1,PSIP,PPSI,HJPSI,VTPSI,TPSI,
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
      RR     = RR1
      RA     = RA1
      RB     = RB1
      RKAP   = RKAP1
      RDLT   = RDLT1
      BB     = BB1
      RIP    = RIP1
      NTRMAX = NTRMAX1
C
      DO NTR=1,NTRMAX
         PSITR(NTR)=PSIP(NTR)/PSIP(NTRMAX)
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
      UJPSI0(1)=0.D0
      DO NTR=2,NTRMAX
         PSIN=PSIP(NTR)/PSIP(NTRMAX)-1.D-8
         CALL SPL1DF(PSIN,HJPSIL,PSITR,UJPSI,NTRMAX,IERR)
         CALL SPL1DI(PSIN,HJPSID,PSITR,UJPSI,UJPSI0,NTRMAX,IERR)
         UJPSI0(NTR)=UJPSI0(NTR-1)+HJPSID
C         WRITE(6,'(I5,1P4E12.4)') NTR,PSIN,HJPSIL,HJPSID,UJPSI0(NTR)
      ENDDO
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
