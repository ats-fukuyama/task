C     $Id$
C
C     ****** INITIALIZE INPUT PARAMETERS ******
C
      SUBROUTINE PLINIT
C
      INCLUDE '../pl/plcomm.h'
C
C     ======( PHYSICAL CONSTANTS )======
C
C        PI    : Pi
C        AEE   : Elementaty charge
C        AME   : Electron mass
C        AMP   : Proton mass
C        VC    : Speed of light in vacuum
C        RMU0  : Permeability of free space
C        EPS0  : Permittivity of free space
C
      PI    = 2.D0*ACOS(0.D0)
      AEE   = 1.6021892  D-19
      AME   = 9.109534   D-31
      AMP   = 1.6726485  D-27
      VC    = 2.99792458 D  8
      RMU0  = 4.D0*PI*1.D-7
      EPS0  = 1.D0/(VC*VC*RMU0)
C
C     ======( DEVICE PARAMETERS )======
C
C        RR    : Plasma major radius                             (m)
C        RA    : Plasma minor radius                             (m)
C        RB    : Wall minor radius                               (m)
C        RKAP  : Plasma shape elongation
C        RDEL  : Plasma shape triangularity *
C        BB    : Magnetic field at center                        (T)
C        Q0    : Safety factor at center
C        QA    : Safety factor on plasma surface
C        RIP   : Plasma current                                 (MA)
C        PROFJ : Curren density profile parameter (power of (1 - rho^2))
C
      RR    = 3.D0
      RA    = 1.D0
      RB    = 1.2D0
      RKAP  = 1.D0
      RDLT  = 0.D0
C
      BB    = 3.D0
      Q0    = 1.D0
      QA    = 3.D0
      RIP   = 3.D0
      PROFJ = 2.D0
C
C     ======( PLASMA PARAMETERS )======
C
C        NSMAX : Number of particle species
C        PA    : Mass number
C        PZ    : Charge number
C        PN    : Density at center                     (1.0E20/m**3)
C        PNS   : Density on plasma surface             (1.0E20/m**3)
C        PZCL  : Ratio of collision frequency to wave frequency
C        PTPR  : Parallel temperature at center                (keV)
C        PTPP  : Perpendicular temperature at center           (keV)
C        PTS   : Temperature on surface                        (keV)
C        PU    : Toroidal rotation velocity at center          (m/s)
C        PUS   : Toroidal rotation velocity on surface         (m/s)
C        PNITB : Density increment at ITB              (1.0E20/Mm*3)
C        PTITB : Temperature increment at ITB                  (keV)
C        PUITB : Toroidal rotation velocity increment at ITB   (m/s)
C
      NSMAX = MIN(2,NSM)
C
         PA(1)   = AME/AMP
         PZ(1)   =-1.0D0
         PN(1)   = 1.0D0
         PNS(1)  = 0.0D0
         PZCL(1) = 0.00D0
         PTPR(1) = 5.0D0
         PTPP(1) = 5.0D0
         PTS(1)  = 0.05D0
         PU(1)   = 0.D0
         PUS(1)  = 0.D0
         PNITB(1)= 0.D0
         PTITB(1)= 0.D0
         PUITB(1)= 0.D0
C
      IF(NSM.GE.2) THEN
         PA(2)   = 1.0D0
         PZ(2)   = 1.0D0
         PN(2)   = 1.0D0
         PNS(2)  = 0.0D0
         PZCL(2) = 0.00D0
         PTPR(2) = 5.0D0
         PTPP(2) = 5.0D0
         PTS(2)  = 0.05D0
         PU(2)   = 0.D0
         PUS(2)  = 0.D0
         PNITB(2)= 0.D0
         PTITB(2)= 0.D0
         PUITB(2)= 0.D0
      ENDIF
C
      DO NS=3,NSM
         PA(NS)   = 1.0D0
         PZ(NS)   = 1.0D0
         PN(NS)   = 0.0D0
         PNS(NS)  = 0.0D0
         PZCL(NS) = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.0D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
      ENDDO
C
C     ======( PROFILE PARAMETERS )======
C
C
C        PROFN1: Density profile parameter (power of rho)
C        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
C        PROFT1: Temperature profile parameter (power of rho)
C        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
C        PROFU1: Rotation profile parameter (power of rho)
C        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))
C
      PROFN1= 2.D0
      PROFN2= 0.5D0
      PROFT1= 2.D0
      PROFT2= 1.D0
      PROFU1= 2.D0
      PROFU2= 1.D0
C
C     ======( MODEL PARAMETERS )======
C
C        MODELG: Control plasma geometry model
C                   0: Slab geometry
C                   1: Cylindrical geometry
C                   2: Toroidal geometry
C                   3: TASK/EQ output geometry
C                   4: VMEC output geometry
C        MODELN: Control plasma profile
C                   0: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS
C                   1: PT calculated from TASK/EQ pressure profile
C                   2: PN*PT proportional to TASK/EQ pressure profile
C                   9: Read from file by means of WMXPRF routine
C        MODELQ: Control safety factor profile (for MODELG=0,1,2)
C                   0: Parabolic q profile (Q0,QA,RHOMIN,RHOITB)
C                   1: Given current profile (RIP,PROFJ)
C
      MODELG= 2
      MODELN= 0
      MODELQ= 0
C
C        RHOMIN: rho at minimum q (0 for positive shear)
C        QMIN  : q minimum for reversed shear
C        RHOITB: rho at ITB (0 for no ITB)
C        RHOEDG: rho at EDGE for smoothing (1 for no smooth)
C
      RHOMIN = 0.D0
      QMIN   = 1.5D0
      RHOITB = 0.D0
      RHOEDG = 1.D0
C
C     ======( MODEL PARAMETERS )======
C
C        KNAMEQ: Filename of equilibrium data
C        KNAMWR: Filename of ray tracing data
C        KNAMFP: Filename of Fokker-Planck data
C
      KNAMEQ = 'eqdata'
      KNAMWR = 'wrdata'
      KNAMFP = 'fpdata'
C
      IDEBUG = 0
C
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE PLPARM(IERR)
C
      INCLUDE '../pl/plcomm.h'
C
      LOGICAL LEX
      CHARACTER KPNAME*32
      NAMELIST /PL/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,
     &              MODELG,MODELN,MODELQ,
     &              KNAMEQ,KNAMWR,KNAMFP,IDEBUG
C
    1 WRITE(6,*) '## INPUT : &PL'
      READ(5,PL,ERR=1,END=9000)
      CALL PLCHEK(IERR)
C
 9000 RETURN
C
      ENTRY PLPARF(IERR)
C
      KPNAME='plparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(LEX) THEN
         OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
         READ(25,PL,ERR=9800,END=9900)
         CLOSE(25)
         WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
      ENDIF
      CALL PLCHEK(IERR)
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      END
C
C     ****** INPUT PARAMETER CHECK******
C
      SUBROUTINE PLCHEK(IERR)
C
      INCLUDE '../pl/plcomm.h'
C
      IERR=0
C
      IF(NSMAX.GT.NSM) THEN
         WRITE(6,*) 'XX PLPARM: NSMAX.GT.NSM: NSMAX,NSM=',NSMAX,NSM
         IERR=1
      ENDIF
C
      IF(MODELG.NE.3) THEN
         IF((MODELG.LT.0).OR.(MODELG.GT.2)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELG: MODELG=',MODELG
            IERR=1
         ENDIF
         IF((MODELN.NE.0).AND.(MODELN.NE.9)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ELSE
         IF((MODELN.NE.0).AND.(MODELN.NE.1).AND.(MODELN.NE.2).AND.
     &      (MODELN.NE.9)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1).AND.(MODELG.NE.3)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ENDIF
C
      RHOE0=0.D0
      RHOES=0.D0
      DO NS=1,NSMAX
         RHOE0=RHOE0+PZ(NS)*PN(NS)
         RHOES=RHOES+PZ(NS)*PNS(NS)
      ENDDO
      IF(ABS(RHOE0).GT.1.D-10) THEN
         WRITE(6,*) 'XX PLPARM: CHARGE NEUTRALITY ERROR AT CENTER'
         IERR=1
      ENDIF
      IF(ABS(RHOES).GT.1.D-10) THEN
         WRITE(6,*) 'XX PLPARM: CHARGE NEUTRALITY ERROR AT SURFACE'
         IERR=1
      ENDIF
C
      RETURN
      END
C
C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE PLVIEW
C
      INCLUDE '../pl/plcomm.h'
C
      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    ,
     &             'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'RKAP  ',RKAP  ,'RDLT  ',RDLT  ,
     &             'Q0    ',Q0    ,'QA    ',QA
      WRITE(6,601) 'RIP   ',RIP   ,'PROFJ ',PROFJ ,
     &             'PROFN1',PROFN1,'PROFN2',PROFN2
      WRITE(6,601) 'PROFT1',PROFT1,'PROFT2',PROFT2,
     &             'PROFU1',PROFU1,'PROFU2',PROFU2
      WRITE(6,601) 'RHOEDG',RHOEDG
      WRITE(6,604) 'MODELG',MODELG,'MODELN',MODELN,
     &             'MODELQ',MODELQ
C
      WRITE(6,100)
      DO NS=1,NSMAX
         WRITE(6,110) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS)
      ENDDO
      WRITE(6,120)
      DO NS=1,NSMAX
         WRITE(6,130) NS,PTPR(NS),PTPP(NS),PTS(NS),PU(NS),PUS(NS)
      ENDDO
      IF(RHOITB.GT.0.D0) THEN
         WRITE(6,140)
         DO NS=1,NSMAX
           WRITE(6,150) NS,PNITB(NS),PTITB(NS),PUITB(NS)
         ENDDO
      ENDIF
      RETURN
C
  100 FORMAT(1H ,'NS    PA          PZ          PN          ',
     &           'PNS         PZCL')
  110 FORMAT(1H ,I2,' ',1P5E12.4)
  120 FORMAT(1H ,'NS    PTPR        PTPP        PTS         ',
     &           'PU          PUS')
  130 FORMAT(1H ,I2,' ',1P5E12.4)                               
  140 FORMAT(1H ,'NS    PNITB       PTITB       PUITB')
  150 FORMAT(1H ,I2,' ',1P3E12.4)                               
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  604 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END
