C     $Id$
C
C     ****** DEFAULT PARAMETERS ******
C
      SUBROUTINE EQINIT
C
      INCLUDE '../eq/eqcomm.inc'
C
C     *** PROFILE PARAMETERS ***
C
C        PP0   : Plasma pressure (main component)              (MPa)
C        PP1   : Plasma pressure (sub component)               (MPa)
C        PP2   : Plasma pressure (increment within ITB)        (MPa)
C        PROFP0: Pressure profile parameter
C        PROFP1: Pressure profile parameter
C        PROFP2: Pressure profile parameter
C
C        PPSI=PP0*(1.D0-PSIN**PROFR0)**PROFP0
C    &       +PP1*(1.D0-PSIN**PROFR1)**PROFP1
C    &       +PP2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFP2
C
C        The third term exists for RHO < RHOITB
C
      PP0    = 0.001D0
      PP1    = 0.0D0
      PP2    = 0.0D0
      PROFP0 = 1.5D0
      PROFP1 = 1.5D0
      PROFP2 = 2.0D0
C
C        PJ0   : Current density at R=RR (main component) : Fixed to 1
C        PJ1   : Current density at R=RR (sub component)       (arb)
C        PJ2   : Current density at R=RR (sub component)       (arb)
C        PROFJ0: Current density profile parameter
C        PROFJ1: Current density profile parameter
C        PROFJ2: Current density profile parameter
C
C      HJPSI=-PJ0*(1.D0-PSIN**PROFR0)**PROFJ0
C     &                *PSIN**(PROFR0-1.D0)
C     &      -PJ1*(1.D0-PSIN**PROFR1)**PROFJ1
C     &                *PSIN**(PROFR1-1.D0)
C     &      -PJ2*(1.D0-PSIN**PROFR2)**PROFJ2
C     &                *PSIN**(PROFR2-1.D0)
C
C        The third term exists for RHO < RHOITB
C
      PJ0    = 1.00D0
      PJ1    = 0.0D0
      PJ2    = 0.0D0
      PROFJ0 = 1.5D0
      PROFJ1 = 1.5D0
      PROFJ2 = 1.5D0
C
C        FF0   : Current density at R=RR (main component) : Fixed to 1
C        FF1   : Current density at R=RR (sub component)       (arb)
C        FF2   : Current density at R=RR (sub component)       (arb)
C        PROFF0: Current density profile parameter
C        PROFF1: Current density profile parameter
C        PROFF2: Current density profile parameter
C
C      FPSI=BB*RR
C     &      +FF0*(1.D0-PSIN**PROFR0)**PROFF0
C     &      +FF1*(1.D0-PSIN**PROFR1)**PROFF1
C     &      +FF2*(1.D0-PSIN**PROFR2)**PROFF2
C
C        The third term exists for RHO < RHOITB
C
      FF0    = 1.0D0
      FF1    = 0.0D0
      FF2    = 0.0D0
      PROFF0 = 1.5D0
      PROFF1 = 1.5D0
      PROFF2 = 1.5D0
C
C        PT0   : Plasma temperature (main component)           (keV)
C        PT1   : Plasma temperature (sub component)            (keV)
C        PT2   : Plasma temperature (increment within ITB)     (keV)
C        PTSEQ : Plasma temperature (at surface)               (keV)
C        PROFT0: Temperature profile parameter
C        PROFT1: Temperature profile parameter
C        PROFT2: Temperature profile parameter
C
C        TPSI=PTSEQ+(PT0-PTSEQ)*(1.D0-PSIN**PROFR0)**PROFT0
C    &       +PT1*(1.D0-PSIN**PROFR1)**PROFT1
C    &       +PT2*(1.D0-PSIN/PSIITB)**PROFR2)**PROFT2
C    &       +PTSEQ
C
C        The third term exits for RHO < RHOITB
C
      PT0    = 1.0D0
      PT1    = 0.0D0
      PT2    = 0.0D0
      PTSEQ  = 0.05D0
      PROFT0 = 1.5D0
      PROFT1 = 1.5D0
      PROFT2 = 2.0D0
C
C        PV0   : Toroidal rotation (main component)              (m/s)
C        PV1   : Toroidal rotation (sub component)               (m/s)
C        PV2   : Toroidal rotation (increment within ITB)        (m/s)
C        PROFV0: Velocity profile parameter
C        PROFV1: Velocity profile parameter
C        PROFV2: Velocity profile parameter
C
C        PVSI=PV0*(1.D0-PSIN**PROFR0)**PROFV0
C    &       +PV1*(1.D0-PSIN**PROFR1)**PROFV1
C    &       +PV2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFV2
C
C        The third term exits for RHO < RHOITB
C
      PV0    = 0.0D0
      PV1    = 0.0D0
      PV2    = 0.0D0
      PROFV0 = 1.5D0
      PROFV1 = 1.5D0
      PROFV2 = 2.0D0
C
C        PN0EQ : Plasma number density(constant)
C
      PN0EQ  = 1.D20
C
C        PROFR0: Profile parameter
C        PROFR1: Profile parameter
C        PROFR2: Profile parameter
C        RHOITB: Normalized radius SQRT(PSI/PSIA) at ITB
C
      PROFR0 = 1.D0
      PROFR1 = 2.D0
      PROFR2 = 2.D0
      RHOITB = 0.5D0
C
C        OTC   : Constant OMEGA**2/TPSI
C        HM    : Constant                                       (Am)
C
      OTC = 0.15D0
      HM  = 1.D6
C
C     *** MESH PARAMETERS ***
C
C        NSGMAX: Number of radial mesh points for Grad-Shafranov eq.
C        NTGMAX: Number of poloidal mesh points for Grad-Shafranov eq.
C        NUGMAX: Number of radial mesh points for flux-average quantities
C        NRGMAX: Number of horizontal mesh points in R-Z plane
C        NZGMAX: Number of vertical mesh points in R-Z plane
C        NPSMAX: Number of flux surfaces
C        NRMAX : Number of radial mesh points for flux coordinates
C        NTHMAX: Number of poloidal mesh points for flux coordinates
C        NSUMAX: Number of boundary points
C        NRVMAX: Number of radial mesh of surface averate
C        NTVMAX: Number of poloidal mesh for surface averate
C
      NSGMAX = 32
      NTGMAX = 32
      NUGMAX = 32
C
      NRGMAX = 33
      NZGMAX = 33
      NPSMAX = 21
C
      NRMAX  = 50
      NTHMAX = 64
      NSUMAX = 65
C
      NRVMAX = 50
      NTVMAX = 200
C
C     *** CONTROL PARAMETERS ***
C
C        EPSEQ  : Convergence criterion for equilibrium
C        NLPMAX : Maximum iteration number of EQ
C        EPSNW  : Convergence criterion for newton method
C        DELNW  : Increment for ferivative in newton method
C        NLPNW  : Maximum iteration number in newton method
C
      EPSEQ  = 1.D-6
      NLPMAX = 20
      EPSNW  = 1.D-2
      DELNW  = 1.D-2
      NLPNW  = 20
C
C        MDLEQF : Profile parameter
C            0: given analytic profile  P,Jtoroidal,T,Vph + Ip
C            1: given analytic profile  P,F               + Ip
C            2: given analytic profile  P,Jparallel       + Ip
C            3: given analytic profile  P,Jparallel
C            4: given analytic profile  P,q
C            5: given spline profile    P,Jtoroial,T,Vph  + Ip
C            6: given spline profile    P,F               + Ip
C            7: given spline profile    P,Jparapllel      + Ip
C            8: given spline profile    P,Jparapllel
C            9: given spline profile    P,q
C
C           10: Free boundary
C
      MDLEQF = 0
C
C        MDLEQA : Rho in P(rho), F(rho), q(rho),...
C            0: SQRT(PSIP/PSIPA)
C            1: SQRT(PSIT/PSITA)
C
      MDLEQA = 0
C
C        MDLEQC : Poloidal coordinate parameter
C            0: Poloidal length coordinate
C            1: Boozer coordinate
C
      MDLEQC = 0
C
C        MDLEQX : Free boundary calculation
C            0: Given PSIB and RIPFC
C            1: PSIB adjusted after loop for given RR,RA,RKAP,RDLT
C            2: PSIB adjusted eqch loop for given RR,RA,RKAP,RDLT
C
      MDLEQX = 0
C
C        NPRINT: Level print out
C            0: no print
C            1: print first and last loop
C            2: print all loop
C
      NPRINT= 0
C
C        RGMIN: Minimum R of coputation region [m]
C        RGMAX: Maxmum  R of coputation region [m]
C        ZGMIN: Minimum Z of coputation region [m]
C        ZGMAX: Maxmum  Z of coputation region [m]
C        ZLIMP: Position of upper X points
C        ZLIMM: Position of lower X points
C
      RGMIN = 1.5D0
      RGMAX = 4.5D0
      ZGMIN =-2.0D0
      ZGMAX = 2.0D0
      ZLIMM =-2.5D0
      ZLIMP = 2.5D0
C
C        PSIB(0:5): Multipole moments of poloidal flux PSIRZ on boundary
C
      PSIB(0) =  2.0D0
      PSIB(1) =  0.5D0
      PSIB(2) =  0.D0
      PSIB(3) =  0.D0
      PSIB(4) =  0.D0
      PSIB(5) =  0.D0
C
C        NPFCMAX : Number of poloidal field coils (PFXs)
C        RIPFC(NPFC) : PFC coil current    [MA]
C        RPFC(NPFC)  : PFC coil position R [m]
C        ZPFC(NPFC)  : PFC coil position Z [m]
C        WPFC(NPFC)  : PFC coil width      [m]
C
      NPFCMAX = 0
      DO NPFC=1,NPFCM
         RIPFC(NPFC) = 0.D0
         RPFC(NPFC)  = 3.D0
         ZPFC(NPFC)  =-1.75D0
         WPFC(NPFC)  = 0.75D0
      ENDDO
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE EQPARM(MODE,KIN,IERR)
C
C     MODE=0 : standard namelinst input
C     MODE=1 : namelist file input
C     MODE=2 : namelist line input
C
C     IERR=0 : normal end
C     IERR=1 : namelist standard input error
C     IERR=2 : namelist file does not exist
C     IERR=3 : namelist file open error
C     IERR=4 : namelist file read error
C     IERR=5 : namelist file abormal end of file
C     IERR=6 : namelist line input error
C     IERR=7 : unknown MODE
C     IERR=10X : input parameter out of range
C
      CHARACTER KIN*(*)
      EXTERNAL EQNLIN,EQPLST
C
    1 CALL TASK_PARM(MODE,'EQ',KIN,EQNLIN,EQPLST,IERR)
      IF(IERR.NE.0) RETURN
C
      CALL EQCHEK(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100
C
      RETURN
      END
C
C     ****** INPUT NAMELIST ******
C
      SUBROUTINE EQNLIN(NID,IST,IERR)
C
      INCLUDE 'eqcomm.inc'
C
      NAMELIST /EQ/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
     &              RHOMIN,QMIN,MODELQ,RHOITB,IDEBUG,KNAMEQ,
     &              PP0,PP1,PP2,PROFP0,PROFP1,PROFP2,
     &              FF0,FF1,FF2,PROFF0,PROFF1,PROFF2,
     &              PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2,
     &              PT0,PT1,PT2,PROFT0,PROFT1,PROFT2,PTSEQ,
     &              PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,PN0EQ,
     &              PROFR0,PROFR1,PROFR2,EPSEQ,NLPMAX,
     &              NSGMAX,NTGMAX,NUGMAX,EPSNW,DELNW,NLPNW,
     &              NRGMAX,NZGMAX,RGMIN,RGMAX,ZGMIN,ZGMAX,ZLIMP,ZLIMM,
     &              NPSMAX,NRVMAX,NTVMAX,
     &              NRMAX,NTHMAX,NSUMAX,
     &              MDLEQF,MDLEQC,MDLEQA,MDLEQX,NPRINT,
     &              PSIB,NPFCMAX,RIPFC,RPFC,ZPFC,WPFC
C
      READ(NID,EQ,IOSTAT=IST,ERR=9800,END=9900)
      IERR=0
      RETURN
C
 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE EQPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &EQ : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ'/
     &       9X,'RHOMIN,QMIN,MODELQ,RHOITB,IDEBUG,KNAMEQ'/
     &       9X,'PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2'/
     &       9X,'PP0,PP1,PP2,PROFP0,PROFP1,PROFP2'/
     &       9X,'FF0,FF1,FF2,PROFF0,PROFF1,PROFF2'/
     &       9X,'PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2'/
     &       9X,'PT0,PT1,PT2,PROFT0,PROFT1,PROFT2,PTSEQ'/
     &       9X,'PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,PN0EQ,HM'/
     &       9X,'PROFR0,PROFR1,PROFR2'/
     &       9X,'NSGMAX,NTGMAX,NUGMAX,NRGMAX,NZGMAX,NPSMAX'/
     &       9X,'NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX'/
     &       9X,'MDLEQF,MDLEQC,MDLEQA,MDLEQX,NPRINT'/
     &       9X,'EPSEQ,NLPMAX,EPSNW,DELNW,NLPNW'/
     &       9X,'RGMIN,RGMAX,RZMIN,RZMAX,ZLIMP,ZLIMM'/
     &       9X,'PSIB,NPFCMAX,RIPFC,RPFC,ZPFC,WPFC')
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE EQCHEK(IERR)
C
      INCLUDE 'eqcomm.inc'
      INCLUDE 'eqcom2.inc'
      INCLUDE 'eqcom3.inc'
C
      IERR=0
C
      IF(NSGMAX.GT.NSGM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NSGMAX.GT.NSGM: ',NSGMAX,NSGM
         IERR=1
      ENDIF
      IF(NTGMAX.GT.NTGM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NTGMAX.GT.NTGM: ',NTGMAX,NTGM
         IERR=2
      ENDIF
      IF(NUGMAX.GT.NUGM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NUGMAX.GT.NUGM: ',NUGMAX,NUGM
         IERR=2
      ENDIF
      IF(NRGMAX.GT.NRGM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NRGMAX.GT.NRGM: ',NRGMAX,NRGM
         IERR=3
      ENDIF
      IF(NZGMAX.GT.NZGM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NZGMAX.GT.NZGM: ',NRGMAX,NRGM
         IERR=4
      ENDIF
      IF(NPSMAX.GT.NPSM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NPSMAX.GT.NPSM: ',NPSMAX,NPSM
         IERR=5
      ENDIF
      IF(NRMAX.GT.NRM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NRMAX.GT.NRM: ',NRMAX,NRM
         IERR=6
      ENDIF
      IF(NTHMAX.GT.NTHM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NTHMAX.GT.NTHM: ',NTHMAX,NTHM
         IERR=7
      ENDIF
      IF(NSUMAX.GT.NSUM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NSUMAX.GT.NSUM: ',NSUMAX,NSUM
         IERR=8
      ENDIF
      IF(NRVMAX.GT.NRVM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NRVMAX.GT.NRVM: ',NRVMAX,NRVM
         IERR=9
      ENDIF
      IF(NTVMAX.GT.NTVM) THEN
         WRITE(6,'(A,I8,I8)') 'XX EQCHEK: NTVMAX.GT.NTVM: ',NTVMAX,NTVM
         IERR=10
      ENDIF
      RETURN
      END
C
C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE EQVIEW
C
      INCLUDE 'eqcomm.inc'
C
      WRITE(6,601) 'RR    ',RR,
     &             'RA    ',RA,
     &             'RKAP  ',RKAP,
     &             'RDLT  ',RDLT
      WRITE(6,601) 'BB    ',BB,
     &             'RIP   ',RIP,
     &             'Q0    ',Q0,
     &             'QA    ',QA
      WRITE(6,601) 'RHOMIN',RHOMIN,
     &             'QMIN  ',QMIN,
     &             'RB    ',RB,
     &             'PROFJ ',PROFJ
      WRITE(6,601) 'RGMIN ',RGMIN,
     &             'RGMAX ',RGMAX,
     &             'ZGMIN ',ZGMIN,
     &             'ZGMAX ',ZGMAX
      WRITE(6,601) 'ZLIMP ',ZLIMP,
     &             'ZLIMM ',ZLIMM
      WRITE(6,601) 'PP0   ',PP0,
     &             'PROFP0',PROFP0,
     &             'PJ0   ',PJ0,
     &             'PROFJ0',PROFJ0
      WRITE(6,601) 'PP1   ',PP1,
     &             'PROFP1',PROFP1,
     &             'PJ1   ',PJ1,
     &             'PROFJ1',PROFJ1
      WRITE(6,601) 'FF0   ',FF0,
     &             'PROFF0',PROFF0
      WRITE(6,601) 'FF1   ',FF1,
     &             'PROFF1',PROFF1
      WRITE(6,601) 'FF2   ',FF2,
     &             'PROFF2',PROFF2
      WRITE(6,601) 'PT0   ',PT0,
     &             'PROFT0',PROFT0,
     &             'PV0   ',PV0,
     &             'PROFV0',PROFV0
      WRITE(6,601) 'PT1   ',PT1,
     &             'PROFT1',PROFT1,
     &             'PV1   ',PV1,
     &             'PROFV1',PROFV1
      WRITE(6,601) 'PT2   ',PT2,
     &             'PROFT2',PROFT2,
     &             'PV2   ',PV2,
     &             'PROFV2',PROFV2
      WRITE(6,601) 'PTSEQ ',PTSEQ,
     &             'PN0EQ ',PN0EQ
      WRITE(6,601) 'PROFR0',PROFR0,
     &             'PROFR1',PROFR1,
     &             'PROFR2',PROFR2,
     &             'RHOITB',RHOITB
      WRITE(6,601) 'EPSEQ ',EPSEQ,
     &             'EPSNW ',EPSNW,
     &             'DELNW ',DELNW
      IF(MDLEQF.GE.10.AND.MDLEQF.LT.20) THEN
         WRITE(6,601) 'PSIB:0',PSIB(0),
     &                'PSIB:1',PSIB(1),
     &                'PSIB:2',PSIB(2),
     &                'PSIB:3',PSIB(3)
         WRITE(6,604) 
     &        (NPFC,RIPFC(NPFC),RPFC(NPFC),ZPFC(NPFC),WPFC(NPFC),
     &         NPFC=1,NPFCMAX)
      ENDIF
      WRITE(6,602) 'NSGMAX',NSGMAX,
     &             'NTGMAX',NTGMAX,
     &             'NUGMAX',NUGMAX
      WRITE(6,602) 'NRGMAX',NRGMAX,
     &             'NZGMAX',NZGMAX,
     &             'NRVMAX',NRVMAX,
     &             'NTVMAX',NTVMAX
      WRITE(6,602) 'NRMAX ',NRMAX,
     &             'NTHMAX',NTHMAX,
     &             'NPSMAX',NPSMAX,
     &             'NSUMAX',NSUMAX
      WRITE(6,602) 'MDLEQF',MDLEQF,
     &             'MDLEQC',MDLEQC,
     &             'MDLEQA',MDLEQA,
     &             'MODELQ',MODELQ
      WRITE(6,602) 'MDLEQX',MDLEQX,
     &             'NPRINT',NPRINT,
     &             'NLPMAX',NLPMAX,
     &             'NLPNW ',NLPNW
C
      RETURN
  601 FORMAT(4(A6,'=',1PE11.2:2X))
  602 FORMAT(4(A6,'=',I7:6X))
C  603 FORMAT(4(A6,'=',I7:6X))
  604 FORMAT(' NPFC ','RIPFC ',6X,'RPFC  ',6X,'ZPFC  ',6X,'WPFC'/
     &       (I6,1P4E12.4))
      END
