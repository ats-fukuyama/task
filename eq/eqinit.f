C     $Id$
C
C     ****** DEFAULT PARAMETERS ******
C
      SUBROUTINE EQINIT
C
      INCLUDE 'eqcomm.h'
C
C
C     *** CONSTANTS ****
C
C        PI    : Pi
C        RMU0  : Permeability of free space
C        BLTZ  : Boltzmann constant
C        AN0   : Avogadro's number
C        AMP   : Proton mass
C        RGAS  : Gas constant
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AEE    = 1.60217733D-19
      AMP    = 1.6726231D-27
      BLTZ   = 1.38066D-23
      AN0    = 6.0221367D23
      RGAS   = AN0*BLTZ
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
C        PPSI=PP0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFP0
C    &       +PP1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFP1
C    &       +PP2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFP2
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
C      HJPSI=-PJ0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFJ0
C     &                *(1.D0-PSIN)**(PROFR0-1.D0)
C     &      -PJ1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFJ1
C     &                *(1.D0-PSIN)**(PROFR1-1.D0)
C     &      -PJ2*(1.D0-(1.D0-PSIN)**PROFR2)**PROFJ2
C     &               *(1.D0-PSIN)**(PROFR2-1.D0)
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
C        PH0   : Plasma pressure (main component)              (MPa)
C        PH1   : Plasma pressure (sub component)               (MPa)
C        PH2   : Plasma pressure (increment within ITB)        (MPa)
C        PROFH0: Pressure profile parameter
C        PROFH1: Pressure profile parameter
C        PROFH2: Pressure profile parameter
C
C       PTPSI=PH0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFH0
C    &       +PH1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFH1
C    &       +PH2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFH2
C
C        The third term exists for RHO < RHOITB
C
      PH0    = 0.001D0
      PH1    = 0.0D0
      PH2    = 0.0D0
      PROFH0 = 1.5D0
      PROFH1 = 1.5D0
      PROFH2 = 2.0D0
C
C        PT0   : Plasma temperature (main component)           (keV)
C        PT1   : Plasma temperature (sub component)            (keV)
C        PT2   : Plasma temperature (increment within ITB)     (keV)
C        PTS   : Plasma temperature (at surface)              (keV)
C        PROFT0: Temperature profile parameter
C        PROFT1: Temperature profile parameter
C        PROFT2: Temperature profile parameter
C
C         TSI=PT0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFT0
C    &       +PT1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFT1
C    &       +PT2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFT2
C    &       +PTS
C
C        The third term exits for RHO < RHOITB
C
      PT0    = 1.0D0
      PP1    = 0.0D0
      PP2    = 0.0D0
      PTS    = 0.05D0
      PROFP0 = 1.5D0
      PROFP1 = 1.5D0
      PROFP2 = 2.0D0
C
C        PV0   : Toroidal rotation (main component)              (m/s)
C        PV1   : Toroidal rotation (sub component)               (m/s)
C        PV2   : Toroidal rotation (increment within ITB)        (m/s)
C        PROFV0: Pressure profile parameter
C        PROFV1: Pressure profile parameter
C        PROFV2: Pressure profile parameter
C
C        PVSI=PV0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFV0
C    &       +PV1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFV1
C    &       +PV2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFV2
C
C        The third term exits for RHO < RHOITB
C
      PV0    = 1.D4
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
      NTHMAX = 16
      NSUMAX = 41
C
C     *** CONTROL PARAMETERS ***
C
C        EPSEQ : Convergence criterion for equilibrium
C
      EPSEQ  = 1.D-6
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
      SUBROUTINE EQPARM
C
      INCLUDE 'eqcomm.h'
C
      LOGICAL LEX
      CHARACTER KPNAME*32
C
      NAMELIST /EQ/ RR,BB,RIP,
     &              RA,RKAP,RDLT,RB,
     &              PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2,
     &              PP0,PP1,PP2,PROFP0,PROFP1,PROFP2,
     &              PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2,
     &              PH0,PH1,PH2,PROFH0,PROFH1,PROFH2,
     &              PT0,PT1,PT2,PROFT0,PROFT1,PROFT2,PTS,
     &              PV0,PV1,PV2,PROFV0,PROFV1,PROFV2,
     &              PROFR0,PROFR1,PROFR2,RHOITB,
     &              NSGMAX,NTGMAX,
     &              NRGMAX,NZGMAX,
     &              EPSEQ,
     &              NPSMAX,KNAMEQ,
     &              NRMAX,NTHMAX,NSUMAX
C
    1 WRITE(6,*) '# INPUT &eq :'
      READ(5,EQ,ERR=1,END=9000)
C
 9000 RETURN
C
C
      ENTRY EQPARF
C
      KPNAME='eqparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX,ERR=9800)
      IF(LEX) THEN
         OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
         READ(25,EQ,ERR=9800,END=9900)
         CLOSE(25)
         WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
      ENDIF
      GOTO 9000
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      END
C
C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE EQVIEW
C
      INCLUDE 'eqcomm.h'
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
      WRITE(6,601) 'PP2   ',PP2,
     &             'PROFP2',PROFP2,
     &             'PJ2   ',PJ2,
     &             'PROFJ2',PROFJ2
      WRITE(6,601) 'PH0   ',PH0,
     &             'PROFH0',PROFH0,
     &             'PT0   ',PT0,
     &             'PROFT0',PROFT0
      WRITE(6,601) 'PH1   ',PH1,
     &             'PROFH1',PROFH1,
     &             'PT1   ',PT1,
     &             'PROFT1',PROFT1
      WRITE(6,601) 'PH2   ',PH2,
     &             'PROFH2',PROFH2,
     &             'PT2   ',PT2,
     &             'PROFT2',PROFT2
      WRITE(6,601) 'PV0   ',PV0,
     &             'PROFV0',PROFV0,
     &             'PTS   ',PTS
      WRITE(6,601) 'PV1   ',PV1,
     &             'PROFV1',PROFV1,
     &             'PN0   ',PN0
      WRITE(6,601) 'PV2   ',PV2,
     &             'PROFV2',PROFV2
      WRITE(6,601) 'PROFR0',PROFR0,
     &             'PROFR1',PROFR1,
     &             'PROFR2',PROFR2,
     &             'RHOITB',RHOITB
      WRITE(6,602) 'NSGMAX ',NSGMAX,
     &             'NTGMAX',NTGMAX,
     &             'NRGMAX',NRGMAX,
     &             'NZGMAX',NZGMAX
      WRITE(6,602) 'NRMAX ',NRMAX,
     &             'NTHMAX',NTHMAX,
     &             'NPSMAX',NPSMAX,
     &             'NSUMAX',NSUMAX
C
      RETURN
  601 FORMAT(4(A6,'=',1PE11.2:2X))
  602 FORMAT(4(A6,'=',I7:6X))
      END
