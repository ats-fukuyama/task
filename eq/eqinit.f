C     $Id$
C
C     ****** DEFAULT PARAMETERS ******
C
      SUBROUTINE EQINIT
C
      INCLUDE 'eqcomm.inc'
C
C
C     *** CONSTANTS ****
C
C        PI    : Pi
C        RMU0  : Permeability of free space
C        AMP   : Proton mass
C        AEE   : Electron charge
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AMP    = 1.6726231D-27
      AEE    = 1.60217733D-19
C
C     *** CONFIGURATION PARAMETERS ***
C
C        RR    : Plasma major radius                             (m)
C        RA    : Plasma minor radius                             (m)
C        RB    : Wall minor radius                               (m)
C        RKAP  : Plasma shape elongation
C        RDLT  : Plasma shape triangularity 
C        BB    : Magnetic field at center                        (T)
C        RIP   : Plasma current                                 (MA)
C
      RR     = 3.D0
      RA     = 1.D0
      RB     = RA*1.1D0
      RKAP   = 1.6D0
      RDLT   = 0.25D0
      BB     = 3.D0
      RIP    = 3.D0
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
C        QQ0   : Safety factor on axis for QQ1=QQ2=0
C        QQS   : Safety factor on surface
C        QQ0   : Safety factor 
C        QQ1   : Safety factor (sub component)
C        QQ2   : Safety factor (increment within ITB)
C        PROFQ0: Safety factor profile parameter
C        PROFQ1: Safety factor profile parameter
C        PROFP2: Pressure profile parameter
C
C        QPSI=QQS
C    &       +(QQ0-QQS)*(1.D0-PSIN**PROFR0)**PROFQ0
C    &       +QQ1*(1.D0-PSIN**PROFR1)**PROFQ1
C    &       +QQ2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFQ2
C
C        The third term exists for RHO < RHOITB
C
      QQ0    = 1.D0
      QQS    = 3.D0
      QQ1    = 0.0D0
      QQ2    = 0.0D0
      PROFQ0 = 1.0D0
      PROFQ1 = 1.0D0
      PROFQ2 = 1.0D0
C
C        PT0   : Plasma temperature (main component)           (keV)
C        PT1   : Plasma temperature (sub component)            (keV)
C        PT2   : Plasma temperature (increment within ITB)     (keV)
C        PTS   : Plasma temperature (at surface)               (keV)
C        PROFT0: Temperature profile parameter
C        PROFT1: Temperature profile parameter
C        PROFT2: Temperature profile parameter
C
C        TPSI=PTS+(PT0-PTS)*(1.D0-PSIN**PROFR0)**PROFT0
C    &       +PT1*(1.D0-PSIN**PROFR1)**PROFT1
C    &       +PT2*(1.D0-PSIN/PSIITB)**PROFR2)**PROFT2
C    &       +PTS
C
C        The third term exits for RHO < RHOITB
C
      PT0    = 1.0D0
      PT1    = 0.0D0
      PT2    = 0.0D0
      PTS    = 0.05D0
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
C        PN0 : Plasma number density(constant)
C
      PN0    = 1.D20
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
C        NRGMAX: Number of horizontal mesh points in R-Z plane
C        NZGMAX: Number of vertical mesh points in R-Z plane
C        NPSMAX: Number of flux surfaces
C        NRMAX : Number of radial mesh points for flux coordinates
C        NTHMAX: Number of poloidal mesh points for flux coordinates
C        NSUMAX: Number of boundary points
C
      NSGMAX = 32
      NTGMAX = 32
C
      NRGMAX = 32
      NZGMAX = 32
      NPSMAX = 21
C
      NRMAX  = 50
      NTHMAX = 64
      NSUMAX = 65
C
C     *** CONTROL PARAMETERS ***
C
C        EPSEQ : Convergence criterion for equilibrium
C
      EPSEQ  = 1.D-6
C
C        MDLEQF : Profile parameter
C            0: given analytic profile  P,Jtoroidal,T,Vph
C            1: given analytic profile  P,F,T,Vph
C            2: given analytic profile  P,Jparallel,T,Vph
C            3: given analytic profile  P,q,T,Vph
C            5: given spline profile  P,Jtoroial,T,Vph
C            6: given spline profile  P,F,T,Vph
C            7: given spline profile  P,Jparapllel,T,Vph
C            8: given spline profile  P,q,T,Vph
C
      MDLEQF = 0
C
C        MDLEQC : Poloidal coordinate parameter
C            0: Poloidal length coordinate
C            1: Boozer coordinate
C
      MDLEQC = 0
C
C        NPRINT: Level print out
C
      NPRINT= 0
C
C     *** FILE NAME ***
C
C        KNAMEQ: Filename of equilibrium data
C
      KNAMEQ = 'eqdata'
C
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE EQPARM(MODE,LINE,IERR)
C
C     MODE=0 : normal namelinst input
C     MODE=1 : namelist file input
C     MODE=2 : namelist line input
C
      INCLUDE 'eqcomm.inc'
C
      LOGICAL LEX
      CHARACTER KPNAME*80,KNAME*90,LINE*80,KID*1
C
      NAMELIST /EQ/ RR,BB,RIP,
     &              RA,RKAP,RDLT,RB,
     &              PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2,
     &              PP0,PP1,PP2,PROFP0,PROFP1,PROFP2,
     &              FF0,FF1,FF2,PROFF0,PROFF1,PROFF2,
     &              QQ0,QQ1,QQ2,PROFQ0,PROFQ1,PROFQ2,QQS,
     &              PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2,
     &              PT0,PT1,PT2,PROFT0,PROFT1,PROFT2,PTS,
     &              PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,
     &              PROFR0,PROFR1,PROFR2,RHOITB,EPSEQ,
     &              NSGMAX,NTGMAX,
     &              NRGMAX,NZGMAX,
     &              NPSMAX,KNAMEQ,
     &              NRMAX,NTHMAX,NSUMAX,
     &              MDLEQF,MDLEQC,NPRINT
C
      IF(MODE.EQ.0) THEN
    1    CONTINUE
         WRITE(6,*) '# INPUT &EQ :'
         READ(5,EQ,ERR=2,END=3)
         KID=' '
         GOTO 4
C
    2    CALL EQPLST
         GOTO 1
C
    3    IERR=1
    4    CONTINUE
C
      ELSEIF(MODE.EQ.1) THEN
C
         INQUIRE(FILE=LINE,EXIST=LEX,ERR=9800)
         IF(.NOT.LEX) THEN
            IERR=2
            RETURN
         ENDIF
         OPEN(25,FILE=LINE,IOSTAT=IST,STATUS='OLD',ERR=9100)
         READ(25,EQ,IOSTAT=IST,ERR=9800,END=9900)
         CLOSE(25)
         CALL KTRIM(LINE,KL)
         WRITE(6,*) 
     &     '## FILE (',LINE(1:KL),') IS ASSIGNED FOR PARM INPUT'
C
      ELSE(MODE.EQ.2) THEN
         KNAME=' &EQ '//LINE//' &END'
         WRITE(7,'(A90)') KNAME
         REWIND(7)
         READ(7,EQ,ERR=8,END=8)
         WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
         GOTO 9
    8    CALL EQPLST
         IERR=3
         RETURN
    9    REWIND(7)
      ELSE
         WRITE(6,*) 'XX EQPARM : UNKNOWN MODE =',MODE
         IERR=4
         RETURN
      ENDIF
C
      IERR=0
C
C     PARAMETER CHECK SHOULD BE HERE
C
      IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1
C
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      IERR=5
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR : IOSTAT = ',IST
      IERR=6
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      IERR=7
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
  601 FORMAT(' ','# &EQ : RR,BB,RIP,RA,RKAP,RDLT,RB'/
     &       9X,'PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2'/
     &       9X,'PP0,PP1,PP2,PROFP0,PROFP1,PROFP2'/
     &       9X,'FF0,FF1,FF2,PROFF0,PROFF1,PROFF2'/
     &       9X,'QQ0,QQ1,QQ2,PROFQ0,PROFQ1,PROFQ2,QQS'/
     &       9X,'PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2'/
     &       9X,'PT0,PT1,PT2,PROFT0,PROFT1,PROFT2,PTS'/
     &       9X,'PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,HM'/
     &       9X,'PROFR0,PROFR1,PROFR2,RHOITB,EPSEQ,'/
     &       9X,'NSGMAX,NTGMAX,NRGMAX,NZGMAX,NPSMAX'/
     &       9X,'NRMAX,NTHMAX,NSUMAX,KNAMEQ'/
     &       9X,'MDLEQF,MDLEQC,NPRINT')
      END
C
C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE EQVIEW
C
      INCLUDE 'eqcomm.inc'
C
      WRITE(6,601) 'RR    ',RR,
     &             'BB    ',BB,
     &             'RIP   ',RIP,
     &             'EPSEQ ',EPSEQ
      WRITE(6,601) 'RA    ',RA,
     &             'RKAP  ',RKAP,
     &             'RDLT  ',RDLT,
     &             'RB    ',RB
      WRITE(6,601) 'PP0   ',PP0,
     &             'PROFP0',PROFP0,
     &             'PJ0   ',PJ0,
     &             'PROFJ0',PROFJ0
      WRITE(6,601) 'PP1   ',PP1,
     &             'PROFP1',PROFP1,
     &             'PJ1   ',PJ1,
     &             'PROFJ1',PROFJ1
      WRITE(6,601) 'FF0   ',FF0,
     &             'PROFF0',PROFF0,
     &             'QQ0   ',QQ0,
     &             'PROFQ0',PROFQ0
      WRITE(6,601) 'FF1   ',FF1,
     &             'PROFF1',PROFF1,
     &             'QQ1   ',QQ1,
     &             'PROFQ1',PROFQ1
      WRITE(6,601) 'FF2   ',FF2,
     &             'PROFF2',PROFF2,
     &             'QQ2   ',QQ2,
     &             'PROFQ2',PROFQ2
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
      WRITE(6,601) 'PTS   ',PTS,
     &             'QQS   ',QQS,
     &             'PN0   ',PN0
      WRITE(6,601) 'PROFR0',PROFR0,
     &             'PROFR1',PROFR1,
     &             'PROFR2',PROFR2,
     &             'RHOITB',RHOITB
      WRITE(6,602) 'NSGMAX',NSGMAX,
     &             'NTGMAX',NTGMAX,
     &             'NRGMAX',NRGMAX,
     &             'NZGMAX',NZGMAX
      WRITE(6,602) 'NRMAX ',NRMAX,
     &             'NTHMAX',NTHMAX,
     &             'NPSMAX',NPSMAX,
     &             'NSUMAX',NSUMAX
      WRITE(6,602) 'MDLEQF',MDLEQF,
     &             'MDLEQC',MDLEQC,
     &             'NPRINT',NPRINT
C
      RETURN
  601 FORMAT(4(A6,'=',1PE11.2:2X))
  602 FORMAT(4(A6,'=',I7:6X))
      END
