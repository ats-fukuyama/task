C     $Id$
C
C   *******************************************
C   **    Initialize EQ interface for TR     **
C   *******************************************
C
C     input:
C
C     RR            : Major plasma radius (Geometrical center)   (m)
C     RA            : Minor plasma radius (Horizontal plane)     (m)
C     RKAP          : Elongation
C     RDLT          : Triangularity
C     BB            : Toroidal magnetic field at R=RR            (T)
C
C     output:
C
C     IERR          : Error indicator
C
C   ***************************************************************
C
      SUBROUTINE TREQIN(RR1,RA1,RKAP1,RDLT1,BB1,IERR)
C              
      INCLUDE 'eqcomq.inc'
      INCLUDE 'eqcom4.inc'
      CHARACTER KPNAME*80
      SAVE INIT
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         CALL EQINIT
         KPNAME='eqparm'
         CALL EQPARM(1,KPNAME,IERR)
         INIT=1
      ENDIF
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AEE    = 1.60217733D-19
C
      IERR   = 0
      MDLEQF = 8
      MDLEQA = 2
      NPRINT = 1
      NLPMAX = 20
C
      RR     = RR1
      RA     = RA1
      RB     = RA1*1.2D0
      RKAP   = RKAP1
      RDLT   = RDLT1
      BB     = BB1
C
      RETURN
      END
C
C   *******************************************
C   **    EQ interface for TR     **
C   *******************************************
C
C     input:
C
C     NTRMAX         : Maximum array number
C     RHOTR(NTRMAX)  : Normarized radius
C     PRHO(NTRMAX)   : Pressure                                 (MPa)
C     HJRHO(NTRMAX)  : Plasma current density                (MA/m^2)
C     VTRHO(NTRMAX)  : Toroidal rotation velocity               (m/s)
C     TRHO(NTRMAX)   : Temperature                              (keV)
C
C     output:
C
C     QRHO(NTRMAX)   : Safety factor
C     TTRHO(NTRMAX)  : Bphi R
C     DVRHO(NTRMAX)  : dV/drho
C     DSRHO(NTRMAX)  : dS/drho
C     ARHRRHO(NTRMAX): <|nabla rho|^2/R^2>
C     AIR2RHO(NTRMAX): <RR^2/R^2>
C     ARH1RHO(NTRMAX): <|nabla rho|>
C     ARH2RHO(NTRMAX): <|nabla rho|^2>
C     ABB2RHO(NTRMAX): <B^2/BB^2>
C     AIB2RHO(NTRMAX): <BB^2/B^2>
C     ARHBRHO(NTRMAX): <|nabla rho|^2/B^2>
C     EPSRHO(NTRMAX) : (Bmax-Bmin)/(Bmax+Bmin)
C     RMJRHO(NTRMAX) : (Rmax+Rmin)/2
C     RMNRHO(NTRMAX) : (Rmax-Rmin)/2
C     RKAPRHO(NTRMAX): (Zmax-Zmin)/(Rmax-Rmin)
C     IERR           : Error indicator
C
C   ***************************************************************
C
      SUBROUTINE TREQEX(NTRMAX1,RHOTR1,PRHO,HJRHO,VTRHO,TRHO,
     &                  QRHO,TTRHO,DVRHO,DSRHO,ARHRRHO,AIR2RHO,
     &                  ARH1RHO,ARH2RHO,ABB2RHO,AIB2RHO,ARHBRHO,
     &                  EPSRHO,RMJRHO,RMNRHO,RKAPRHO,IERR)
C              
      INCLUDE 'eqcomq.inc'
      INCLUDE 'eqcom4.inc'
C
      DIMENSION PRHO(NTRMAX1),HJRHO(NTRMAX1)
      DIMENSION VTRHO(NTRMAX1),TRHO(NTRMAX1)
      DIMENSION QRHO(NTRMAX1),TTRHO(NTRMAX1)
      DIMENSION DVRHO(NTRMAX1),DSRHO(NTRMAX1)
      DIMENSION ARHRRHO(NTRMAX1),AIR2RHO(NTRMAX1)
      DIMENSION ARH1RHO(NTRMAX1),ARH2RHO(NTRMAX1)
      DIMENSION ABB2RHO(NTRMAX1),AIB2RHO(NTRMAX1)
      DIMENSION ARHBRHO(NTRMAX1),EPSRHO(NTRMAX1)
      DIMENSION RMJRHO(NTRMAX1),RMNRHO(NTRMAX1),RKAPRHO(NTRMAX1)
      DIMENSION RHOTR1(NTRMAX1)
      DIMENSION WORK(NTRM+2),DERIV(NTRM+2),UJPSIX(NTRM+2)
C
      IERR=0
      NTRMAX = NTRMAX1
C
      DO NTR=1,NTRMAX
         RHOTR(NTR)=RHOTR1(NTR)
      ENDDO
      
      PSITRX(1)=0.D0
      DO NTR=2,NTRMAX+1
         PSITRX(NTR)=RHOTR(NTR-1)**2
      ENDDO
      PSITRX(NTRMAX+2)=1.D0
C
C     *** Calculate SPLINE coefficients for PPSI, JPSI, VTPSI, TPSI ***
C
C     <<< PRHO >>>
      WORK(1)=(9.D0*PRHO(1)-PRHO(2))/8.D0
      DO NTR=1,NTRMAX
         WORK(NTR+1)=PRHO(NTR)
      ENDDO
      WORK(NTRMAX+2)=0.D0
      CALL SPL1D(PSITRX,WORK,DERIV,UPPSI,NTRMAX+2,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1D PRHO: IERR=',IERR
C
C     <<< HJRHO >>>
      WORK(1)=(9.D0*HJRHO(1)-HJRHO(2))/8.D0
      DO NTR=1,NTRMAX
         WORK(NTR+1)=HJRHO(NTR)
      ENDDO
      WORK(NTRMAX+2)=0.D0
      CALL SPL1D(PSITRX,WORK,DERIV,UJPSI,NTRMAX+2,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1D HJRHO: IERR=',IERR
C
C     <<< VTRHO >>>
      WORK(1)=(9.D0*VTRHO(1)-VTRHO(2))/8.D0
      DO NTR=1,NTRMAX
         WORK(NTR+1)=VTRHO(NTR)
      ENDDO
      WORK(NTRMAX+2)=0.D0
      CALL SPL1D(PSITRX,WORK,DERIV,UVTPSI,NTRMAX+2,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1D VTRHO: IERR=',IERR
C
C     <<< TRHO >>>
      WORK(1)=(9.D0*TRHO(1)-TRHO(2))/8.D0
      DO NTR=1,NTRMAX
         WORK(NTR+1)=TRHO(NTR)
      ENDDO
      WORK(NTRMAX+2)=0.D0
      CALL SPL1D(PSITRX,WORK,DERIV,UTPSI,NTRMAX+2,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1D TRHO: IERR=',IERR
C
C     <<< HJPSID >>>
      DO NTR=1,NTRMAX+2
         UJPSI0(NTR)=0.D0
      ENDDO
      UJPSIX(1)=0.D0
      DO NTR=2,NTRMAX+2
         PSIL=PSITRX(NTR)
         CALL SPL1DF(PSIL,HJPSIL,PSITRX,UJPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 
     &        'XX TREQEX: SPL1DF HJPSIL: IERR=',IERR
         CALL SPL1DI(PSIL,HJPSID,PSITRX,UJPSI,UJPSI0,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 
     &        'XX TREQEX: SPL1DI HJPSID: IERR=',IERR
         UJPSIX(NTR)=HJPSID
      ENDDO
C
      DO NTR=NTRMAX+1,1,-1
         UJPSI0(NTR)=UJPSI0(NTR+1)-UJPSIX(NTR+1)
      ENDDO
C
C     ***** Solve GS equation for given profile *****
C
      CALL EQCALC(IERR)
         IF(IERR.NE.0) GOTO 9000
C
C     ***** Calculate eqilibrium quantities *****
C
      NRMAX1=NRMAX
      NTHMAX1=64
      CALL EQCALQ(NRMAX1,NTHMAX1,0,IERR)
C
C     ***** Calculate Q, DVRHO and others at given radial position *****
C
      DO NTR=1,NTRMAX
C
C        <<< RHOT >>>
         RHOTL=RHOTR(NTR)
C
C        <<< QRHO >>>
         CALL SPL1DF(RHOTL,QPL,RHOT,UQPS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF QPL: IERR=',IERR
         QRHO(NTR)=QPL
C
C        <<< TTRHO >>>
         CALL SPL1DF(RHOTL,TTL,RHOT,UTTS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF TTL: IERR=',IERR
         TTRHO(NTR)=TTL
C
C        <<< DVRHO >>>
         CALL SPL1DF(RHOTL,VPL,RHOT,UVPS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF VPL: IERR=',IERR
         DVRHO(NTR)=VPL
C
C        <<< DSRHO >>>
         CALL SPL1DF(RHOTL,SPL,RHOT,USPS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF SPL: IERR=',IERR
         DSRHO(NTR)=SPL
C
C        <<< ARHRRHO >>>
         CALL SPL1DF(RHOTL,ARHR,RHOT,UAVERHR,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVERHR: IERR=',IERR
         ARHRRHO(NTR)=ARHR
C
C        <<< AIR2RHO >>>
         CALL SPL1DF(RHOTL,AIR2,RHOT,UAVEIR2,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVEIR2: IERR=',IERR
         AIR2RHO(NTR)=AIR2
C
C        <<< ARH1RHO >>>
         CALL SPL1DF(RHOTL,ARH1,RHOT,UAVERH1,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVERH1: IERR=',IERR
         ARH1RHO(NTR)=ARH1
C
C        <<< ARH2RHO >>>
         CALL SPL1DF(RHOTL,ARH2,RHOT,UAVERH2,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVERH2: IERR=',IERR
         ARH2RHO(NTR)=ARH2
C
C        <<< ABB2RHO >>>
         CALL SPL1DF(RHOTL,ABB2,RHOT,UAVEBB2,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVEBB2: IERR=',IERR
         ABB2RHO(NTR)=ABB2
C
C        <<< AIB2RHO >>>
         CALL SPL1DF(RHOTL,AIB2,RHOT,UAVEIB2,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVEIB2: IERR=',IERR
         AIB2RHO(NTR)=AIB2
C
C        <<< ARHBRHO >>>
         CALL SPL1DF(RHOTL,ARHB,RHOT,UAVERHB,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF AVERHB: IERR=',IERR
         ARHBRHO(NTR)=ARHB
C
C        <<< RRMIN >>>
         CALL SPL1DF(RHOTL,RRMINL,RHOT,URRMIN,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF RRMIN: IERR=',IERR
C
C        <<< RRMAX >>>
         CALL SPL1DF(RHOTL,RRMAXL,RHOT,URRMAX,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF RRMAX: IERR=',IERR
C
C        <<< ZZMIN >>>
         CALL SPL1DF(RHOTL,ZZMINL,RHOT,UZZMIN,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF ZZMIN: IERR=',IERR
C
C        <<< ZZMAX >>>
         CALL SPL1DF(RHOTL,ZZMAXL,RHOT,UZZMAX,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF ZZMAX: IERR=',IERR
C
C        <<< BBMIN >>>
         CALL SPL1DF(RHOTL,BBMINL,RHOT,UBBMIN,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF BBMIN: IERR=',IERR
C
C        <<< BBMAX >>>
         CALL SPL1DF(RHOTL,BBMAXL,RHOT,UBBMAX,NRMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TREQEX: SPL1DF ZZMAX: IERR=',IERR
C
         EPSRHO(NTR)=(BBMAXL-BBMINL)/(BBMAXL+BBMINL)
         RMJRHO(NTR)=0.5D0*(RRMAXL+RRMINL)
         RMNRHO(NTR)=0.5D0*(RRMAXL-RRMINL)
         RKAPRHO(NTR)=(ZZMAXL-ZZMINL)/(RRMAXL-RRMINL)
      ENDDO
 9000 RETURN
      END
