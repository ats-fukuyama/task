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
C     RIP           : Total plasma current                      (MA)
C     NTRMAX        : Maximum array number
C     RHOTR(NTRMAX) : Normalized radial mesh
C     HJRHO(NTRMAX) : Plasma current density                (MA/m^2)
C     MDLUF         :
C
C     output:
C
C     QRHO(NTRMAX+1)  : Safety factor
C     IERR          : Error indicator
C
C   ***************************************************************
C
      SUBROUTINE TREQIN(RR1,RA1,RKAP1,RDLT1,BB1,RIP1,
     &                  NTRMAX1,RHOTR1,HJRHO,QRHO,MDLUF,IERR)
C              
      INCLUDE 'eqcomq.inc'
      INCLUDE 'eqcom4.inc'
      DIMENSION RHOTR1(NTRMAX1),HJRHO(NTRMAX1),QRHO(NTRMAX1+1)
      DIMENSION DERIV(NTRM+2),BP(NTRM+1),DERIVX(NRM)
      DIMENSION ASR(NTRM)
      CHARACTER KFILE*10,KPNAME*80
      SAVE INIT
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         CALL EQINIT
         KPNAME='eqparm'
         CALL EQPARF(KPNAME)
         INIT=1
      ENDIF
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AEE    = 1.60217733D-19
C
      IERR   = 0
      MDLEQF = 5
      RR     = RR1
      RA     = RA1
      RB     = RA1*1.2D0
      RKAP   = RKAP1
      RDLT   = RDLT1
      BB     = BB1
      RIP    = RIP1
      NTRMAX = NTRMAX1
C
C     ***** Define radial array *****
C     ----- RHOTR:  half-integer mesh in TR
C     ----- RHOTRG: integer mesh deduced from RHOTR
C     ----- RHOTRX: half-integer mesh with 0 and 1 for SPLINE
C
      DO NTR=1,NTRMAX-1
         RHOTR(NTR)=RHOTR1(NTR)
         RHOTRG(NTR+1)=0.5D0*(RHOTR1(NTR)+RHOTR1(NTR+1))
         RHOTRX(NTR+1)=RHOTR1(NTR)
      ENDDO
      RHOTRG(1)=0.D0
      RHOTRX(1)=0.D0
      RHOTR(NTRMAX)=RHOTR1(NTRMAX)
      RHOTRG(NTRMAX+1)=1.D0
      RHOTRX(NTRMAX+1)=RHOTR1(NTRMAX)
      RHOTRX(NTRMAX+2)=1.D0
C
C     ***** Calculate interim poloidal magnetic field *****
C
      BP(1)=0.D0
      DO NTR=2,NTRMAX+1
         DRHO=RHOTRG(NTR)-RHOTRG(NTR-1)
         RHOH=0.5D0*(RHOTRG(NTR)+RHOTRG(NTR-1))
         HJRH=HJRHO(NTR-1)/RKAP
         BP(NTR)=(RHOTRG(NTR-1)*BP(NTR-1)+RMU0*RHOH*DRHO*HJRH)
     &           /RHOTRG(NTR)
      ENDDO
      BPA= RMU0*RIP*1.D6/(2.D0*PI*RA*RKAP)
      FACT=BPA/BP(NTRMAX+1)
      DO NTR=2,NTRMAX+1
         BP(NTR)=FACT*BP(NTR)
      ENDDO
      DO NTR=2,NTRMAX+1
         QRHO(NTR)=-RHOTRG(NTR  )*RA*BB/(RR*BP(NTR  ))
      ENDDO
      QRHO(1)=(4.D0*QRHO(2)-QRHO(3))/3.D0
C         DO NTR=1,NTRMAX+1
C            write(6,*) NTR,RHOTRG(NTR),QRHO(NTR)
C         ENDDO
C
C     ***** Calculate interim poloidal flux PSIRHO *****
C
      PSIRHO(1)=0.D0
      DO NTR=2,NTRMAX+1
C         DPSIT=PI*RKAP*(RHOTRG(NTR)**2-RHOTRG(NTR-1)**2)*BB
C         PSIRHO(NTR)=PSIRHO(NTR-1)-2.D0*DPSIT/QRHO(NTR-1)
         DPSIT=PI*RKAP*(RHOTRG(NTR)**2-RHOTRG(NTR-1)**2)*BB*RA**2
         PSIRHO(NTR)=PSIRHO(NTR-1)-DPSIT/QRHO(NTR-1)
C         WRITE(6,'(A,I5,1P3E12.4)') 
C     &        'NTR,RHOTR,QRHO,PSIRHO=',
C     &        NTR,RHOTRG(NTR),QRHO(NTR),PSIRHO(NTR)
      ENDDO
C      PAUSE
      SAXIS=-PSIRHO(NTRMAX+1)
      DO NTR=1,NTRMAX+1
         PSIRHO(NTR)=SAXIS+PSIRHO(NTR)
      ENDDO
C
C     ***** Calculate SPLINE for PSIRHO vs RHO *****
C
      DERIV(1)=0.D0
      CALL SPL1D(RHOTRG,PSIRHO,DERIV,UPSIRHO,NTRMAX+1,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D PSIRHO : IERR=',IERR
C
C     ***** Calculate interim toroidal flux FTS *****
C
      DR=RB/(NRMAX-1)
      DO NR=1,NRMAX
         RHOL=DR*(NR-1)/RA
         IF(RHOL.GT.1.D0) THEN
            FTS(NR)=PI*RKAP*(RA*RHOL)**2*BB
            PSS(NR)=SAXIS-SAXIS*RHOL
         ELSE
            FTS(NR)=PI*RKAP*(RA*RHOL)**2*BB
            CALL SPL1DF(RHOL,PSIL,RHOTRG,UPSIRHO,NRMAX+1,IERR)
            IF(IERR.NE.0) 
     &           WRITE(6,*) 'XX TREQIN: SPL1DF PSS : IERR=',IERR
            PSS(NR)=PSIL
         ENDIF
C         WRITE(6,'(A,I5,1P4E12.4)') 
C     &           'NR,RHOL,PSS,FTS=',NR,RHOL,PSS(NR),FTS(NR)
      ENDDO
C      PAUSE
C
C     ***** Calculate SPLINE for PSS vs FTS and FTS vs PSS *****
C
      CALL SPL1D(FTS,PSS,DERIVX,UFTT,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D PSS: IERR=',IERR
      CALL SPL1D(PSS,FTS,DERIVX,UFTS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN: SPL1D FTS: IERR=',IERR
C
C     ***** Setup 2D mesh and initial psi *****
C
      CALL EQMESH
      CALL EQPSIN
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
C     RIP           : Total plasma current                      (MA)
C     NTRMAX        : Maximum array number
C     PRHO(NTRMAX)  : Pressure                                 (MPa)
C     HJRHO(NTRMAX) : Plasma current density                (MA/m^2)
C     VTRHO(NTRMAX) : Toroidal rotation velocity               (m/s)
C     TRHO(NTRMAX)  : Temperature                              (keV)
C
C     output:
C
C     QRHO(NTRMAX)  : Safety factor
C     TTRHO(NTRMAX) : Bphi R
C     DVRHO(NTRMAX) : dV/drho
C     DSRHO(NTRMAX) : dS/drho
C     ABRHO(NTRMAX) : <(nabla rho)^2/R^2>
C     ARRHO(NTRMAX) : <1/R^2>
C     AR1RHO(NTRMAX): <nabla rho>
C     AR2RHO(NTRMAX): <(nabla rho)^2>
C     EPSRHO(NTRMAX): (Bmax-Bmin)/(Bmax+Bmin)
C     MDLUF         :
C     IERR          : Error indicator
C
C   ***************************************************************
C
      SUBROUTINE TREQEX(RIP1,NTRMAX1,PRHO,HJRHO,VTRHO,TRHO,
     &                  QRHO,TTRHO,DVRHO,DSRHO,ABRHO,ARRHO,
     &                  AR1RHO,AR2RHO,EPSRHO,MDLUF,IERR)
C              
      INCLUDE 'eqcomq.inc'
      INCLUDE 'eqcom4.inc'
      DIMENSION PRHO(NTRMAX1),HJRHO(NTRMAX1)
      DIMENSION VTRHO(NTRMAX1),TRHO(NTRMAX1)
      DIMENSION QRHO(NTRMAX1),TTRHO(NTRMAX1)
      DIMENSION DVRHO(NTRMAX1),DSRHO(NTRMAX1)
      DIMENSION ABRHO(NTRMAX1),ARRHO(NTRMAX1)
      DIMENSION AR1RHO(NTRMAX1),AR2RHO(NTRMAX1),EPSRHO(NTRMAX1)
      DIMENSION WORK(NTRM+2),DERIV(NTRM+2),UJPSIX(NTRM+2)
      DIMENSION ASR(NRM)
      DIMENSION FFQ(NRM),FFTT(NRM),FFVOL(NRM),FFDVOL(NRM)
      DIMENSION FFAREA(NRM),FFRMJ(NRM),FFGR1(NRM),FFGR2(NRM)
      CHARACTER KFILE*10
C
      IERR=0
      RIP    = RIP1
      NTRMAX = NTRMAX1
      MDQ=0
      MDVOL=0
      MDAREA=0
      MDRMJ=0
      MDGR1=0
      MDGR2=0
      MDALL=0
      MDTT =0
C
C     ***** Calculate poloidal flux PSITR/G/X at RHOTR/G/X *****
C
      FTSA=FNFTS(1.D0)
      DO NTR=1,NTRMAX
         FTL=FTSA*RHOTR(NTR)**2
         CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
         PSITR(NTR)=RHOTR(NTR)**2
C         PSITR(NTR)=PSIL/SAXIS
C
         FTL=FTSA*RHOTRG(NTR)**2
         CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
         PSITRG(NTR)=RHOTRG(NTR)**2
         RHOT(NTR)=RHOTRG(NTR)
C
         FTL=FTSA*RHOTRX(NTR)**2
         CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
         PSITRX(NTR)=RHOTRX(NTR)**2
C         PSITRX(NTR)=PSIL/SAXIS
C         WRITE(6,'(A,I5,1P4E12.4)') 
C     &        'NTR,FTL,PSIL,PSIN=',NTR,FTL,
C     &        PSITR(NTR),PSITRG(NTR),PSITRX(NTR)
      ENDDO
         FTL=FTSA*RHOTRG(NTRMAX+1)**2
         CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
         PSITRG(NTRMAX+1)=RHOTRG(NTRMAX+1)**2
C         PSITRG(NTRMAX+1)=PSIL/SAXIS
C
         FTL=FTSA*RHOTRX(NTRMAX+1)**2
         CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
         PSITRX(NTRMAX+1)=RHOTRX(NTRMAX+1)**2
C         PSITRX(NTRMAX+1)=PSIL/SAXIS
C
C         DO NTR=1,NTRMAX+1
C         write(6,*) NTR,FTSA*RHOTRX(NTR+1)**2-FTSA*RHOTRX(NTR)**2
C         ENDDO
         FTL=FTSA*RHOTRX(NTRMAX+2)**2
         CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
         PSITRX(NTRMAX+2)=RHOTRX(NTRMAX+2)**2
C         PSITRX(NTRMAX+2)=PSIL/SAXIS
C
C      PAUSE
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
C      DO NTR=1,NTRMAX+2
C         write(6,*) PSITRX(NTR),WORK(NTR)
C      ENDDO
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
C     <<< HJPSI & HJPSID >>>
      DO NTR=1,NTRMAX+2
         UJPSI0(NTR)=0.D0
      ENDDO
      UJPSIX(1)=0.D0
      DO NTR=2,NTRMAX+2
C         PSIL=PSITRX(NTR)+1.D-8
         PSIL=PSITRX(NTR)
         CALL SPL1DF(PSIL,HJPSIL,PSITRX,UJPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 
     &        'XX TREQEX: SPL1DF HJPSIL: IERR=',IERR
         CALL SPL1DI(PSIL,HJPSID,PSITRX,UJPSI,UJPSI0,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 
     &        'XX TREQEX: SPL1DI HJPSID: IERR=',IERR
         UJPSIX(NTR)=HJPSID
      ENDDO
C      UJPSI0(NTRMAX+2)=0.D0
      DO NTR=NTRMAX+1,1,-1
         UJPSI0(NTR)=UJPSI0(NTR+1)-UJPSIX(NTR+1)
C         WRITE(6,'(A,I5,1P3E12.4)') 
C     &        'NTR,PSIL,HJPSID/0=',
C     &        NTR,PSITRX(NTR),UJPSIX(NTR),UJPSI0(NTR)
      ENDDO
C      PAUSE
C
C     ***** Solve GS equation for given profile *****
C
      CALL EQLOOP(IERR)
         IF(IERR.NE.0) GOTO 9000
      CALL EQTORZ
      CALL EQCALP
C      CALL EQGOUT(1)
C
C     ***** Calculate eqilibrium quantities *****
C
      NRMAX1=NRMAX
      NTHMAX1=64
      CALL EQCALQ(NRMAX1,NTHMAX1,0,IERR)
C
C     ***** Calculate Spline Coefficients *****
C
C
C     ***** Calculate Q, DVRHO and others at given radial position *****
C
      PSITA=FNFTS(1.D0)
      DO NTR=1,NTRMAX
C
C        <<< PSIL >>>
         PSITL=PSITA*RHOTR(NTR)**2
         DPSITDRHO=2.D0*PSITA*RHOTR(NTR)
         CALL SPL1DF(PSITL,PSIL,FTS,UFTT,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF PSIL: IERR=',IERR
C     
C        <<< QRHO >>>
         IF(MDQ.EQ.0) THEN
            CALL SPL1DF(PSIL,QPL,PSS,UQPS,NRMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF QPL: IERR=',IERR
C            WRITE(6,'(I5,1P3E12.4)') NTR,RHOTR(NTR),PSIL,QPL
         ELSEIF(MDQ.EQ.1) THEN
            CALL SPL1DF(PSIL,QPL,ASR,UQEQ,NUFMAX,IERR)
            IF(IERR.NE.0) 
     &           WRITE(6,*) 'XX TREQEX: SPL1DF QPL2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDQ=0
         ENDIF
         QRHO(NTR)=QPL
         DPSIPDRHO=DPSITDRHO/QPL
C         write(6,*) NTR,PSIL,QRHO(NTR)
C
C        <<< TTRHO >>>
         IF(MDTT.EQ.0) THEN
            CALL SPL1DF(PSIL,TTL,PSS,UTTS,NRMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF TTL: IERR=',IERR
         ELSEIF(MDTT.EQ.1.OR.MDTT.EQ.2) THEN
            CALL SPL1DF(PSIL,TTL,ASR,UTTEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF TTL2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDTT=0
         ENDIF
         TTRHO(NTR)=TTL
C         write(6,*) NTR,PSIL,TTL
C
C        <<< DVRHO >>>
         IF(MDVOL.EQ.0) THEN
            CALL SPL1DF(PSIL,VPL,PSS,UVPS,NRMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF VPL: IERR=',IERR
         ELSEIF(MDVOL.EQ.1) THEN
            CALL SPL1DF(PSIL,VPL,ASR,UDVOLEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF VPL2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDVOL=0
         ENDIF
         DVRHO(NTR)=VPL*DPSIPDRHO
C         write(6,*) NTR,PSIL,DVRHO(NTR)
C
C        <<< DSRHO >>>
         IF(MDAREA.EQ.0) THEN
            CALL SPL1DF(PSIL,SPL,PSS,USPS,NRMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF SPL: IERR=',IERR
            DSRHO(NTR)=SPL*DPSIPDRHO
         ELSEIF(MDAREA.EQ.1) THEN
            CALL SPL1DD(PSIL,SPL,DSPL,ASR,UAREAEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF SPL2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDAREA=0
            DSRHO(NTR)=DSPL*DPSIPDRHO
         ENDIF
C
C        <<< ABRHO >>>
         IF(MDRMJ.EQ.0) THEN
            CALL SPL1DF(PSIL,AVBRL,PSS,UAVBR,NRMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVBRL: IERR=',IERR
            ABRHO(NTR)=AVBRL
         ELSEIF(MDRMJ.EQ.1) THEN
            CALL SPL1DF(PSIL,AVRRL,ASR,URMJEQ,NUFMAX,IERR)
            CALL SPL1DF(PSIL,AVR2L,ASR,UGR2EQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVBRL2: IERR=',IERR
            ABRHO(NTR)=(AVRRL*AVR2L)*(PSIL-PSS(1))/DPSIPDRHO**2
         ENDIF
C
C        <<< ARRHO >>>
         IF(MDRMJ.EQ.0) THEN
            CALL SPL1DF(PSIL,AVRRL,PSS,UAVRR,NRMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVRRL: IERR=',IERR
            ARRHO(NTR)=AVRRL
         ELSEIF(MDRMJ.EQ.1) THEN
            CALL SPL1DF(PSIL,AVRRL,ASR,URMJEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVRRL2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDRMJ=0
            ARRHO(NTR)=1.D0/(AVRRL**2)
         ENDIF
C
C        <<< AR1RHO >>>
         IF(MDGR1.EQ.0) THEN
            CALL SPL1DF(PSIL,AVR1L,PSS,UAVR1,NRMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVR1L: IERR=',IERR
            AR1RHO(NTR)=AVR1L
         ELSEIF(MDGR1.EQ.1) THEN
            CALL SPL1DF(PSIL,AVR1L,ASR,UGR1EQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVR1L2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDGR1=0
            AR1RHO(NTR)=AVR1L*SQRT(PSIL-PSS(1))/DPSIPDRHO
         ENDIF
C         write(6,*) NTR,PSIL,AR1RHO(NTR)
C
C        <<< AR2RHO >>>
         IF(MDGR2.EQ.0) THEN
            CALL SPL1DF(PSIL,AVR2L,PSS,UAVR2,NRMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVR2L: IERR=',IERR
            AR2RHO(NTR)=AVR2L
         ELSEIF(MDGR2.EQ.1) THEN
            CALL SPL1DF(PSIL,AVR2L,ASR,UGR2EQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TREQEX: SPL1DF AVR2L2: IERR=',IERR
            IF(NTR.EQ.NTRMAX) MDGR2=0
            AR2RHO(NTR)=AVR2L*(PSIL-PSS(1))/DPSIPDRHO**2
         ENDIF
C         write(6,*) NTR,PSIL,AR2RHO(NTR)
C
C        <<< EPSRHO >>>
         CALL SPL1DF(PSIL,BBMINL,PSS,UBBMIN,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF BBMINL: IERR=',IERR
         CALL SPL1DF(PSIL,BBMAXL,PSS,UBBMAX,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQEX: SPL1DF BBMAXL: IERR=',IERR
         EPSRHO(NTR)=(BBMAXL-BBMINL)/(BBMAXL+BBMINL)
C         WRITE(6,'(I5,1P3E12.4)') NTR,BBMINL,BBMAXL,EPSRHO(NTR)
      ENDDO
C      PAUSE
 9000 RETURN
      END
C
