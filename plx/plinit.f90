!     $Id$

  MODULE plinit

  CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

    SUBROUTINE pl_init

      USE plcomm
      IMPLICIT NONE
      INTEGER:: NS

      CALL pl_allocate_ns

!     ======( DEVICE PARAMETERS )======

!        RR    : Plasma major radius                             (m)
!        RA    : Plasma minor radius                             (m)
!        RB    : Wall minor radius                               (m)
!        RKAP  : Plasma shape elongation
!        RDEL  : Plasma shape triangularity *
!        BB    : Magnetic field at center                        (T)
!        Q0    : Safety factor at center
!        QA    : Safety factor on plasma surface
!        RIP   : Plasma current                                 (MA)
!        PROFJ : Curren density profile parameter (power of (1 - rho^2))

      RR    = 3.D0
      RA    = 1.D0
      RB    = 1.2D0
      RKAP  = 1.D0
      RDLT  = 0.D0

      BB    = 3.D0
      Q0    = 1.D0
      QA    = 3.D0
      RIP   = 3.D0
      PROFJ = 2.D0

!     ======( PLASMA PARAMETERS )======

!        NSMAX : Number of particle species
!        PA    : Mass number
!        PZ    : Charge number
!        PZ0   : Atomic number (-1 for electron)
!        PN    : Density at center                     (1.0E20/m**3)
!        PNS   : Density on plasma surface             (1.0E20/m**3)
!        PTPR  : Parallel temperature at center                (keV)
!        PTPP  : Perpendicular temperature at center           (keV)
!        PTS   : Temperature on surface                        (keV)
!        PU    : Toroidal rotation velocity at center          (m/s)
!        PUS   : Toroidal rotation velocity on surface         (m/s)
!        PNITB : Density increment at ITB              (1.0E20/Mm*3)
!        PTITB : Temperature increment at ITB                  (keV)
!        PUITB : Toroidal rotation velocity increment at ITB   (m/s)
!        PZCL  : normalized collision frequency

!        KIDNS : index of particle species
!        IDION :  1 : fast ion particle
!                 0 : else

      NSMAX = 2                  ! Default number of particle species

         ! electron
         NS = 1

         KIDNS(NS)= 'e'
         IDION(NS)= 0.0D0
         PA(NS)   = AME/AMP
         PZ(NS)   =-1.0D0
         PZ0(NS)  =-1.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PZCL(NS) = 0.D0

      IF(NSM.GE.2) THEN
<<<<<<< plinit.f90
         ! *** hydrogen ***
!!$         NS = 2
!!$
!!$         KIDNS(NS)= 'H'
!!$         IDION(NS)= 0.0D0
!!$         PA(NS)   = 1.0D0
!!$         PZ(NS)   = 1.0D0
!!$         PZ0(NS)  = 1.0D0
!!$         PN(NS)   = 1.0D0
!!$         PNS(NS)  = 0.0D0
!!$         PTPR(NS) = 5.0D0
!!$         PTPP(NS) = 5.0D0
!!$         PTS(NS)  = 0.05D0
!!$         PU(NS)   = 0.D0
!!$         PUS(NS)  = 0.D0
!!$         PNITB(NS)= 0.D0
!!$         PTITB(NS)= 0.D0
!!$         PUITB(NS)= 0.D0

=======
>>>>>>> 1.7
         ! *** deuterium ***
         NS = 2

         KIDNS(NS)= 'D'
         IDION(NS)= 0.0D0
         PA(NS)   = 2.0D0
         PZ(NS)   = 1.0D0
         PZ0(NS)  = 1.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
<<<<<<< plinit.f90

         ! *** tritium ***
!!$         NS = 3
!!$
!!$         KIDNS(NS)= 'T'
!!$         IDION(NS)= 0.0D0
!!$         PA(NS)   = 3.0D0
!!$         PZ(NS)   = 1.0D0
!!$         PZ0(NS)  = 1.0D0
!!$         PN(NS)   = 1.0D0
!!$         PNS(NS)  = 0.0D0
!!$         PTPR(NS) = 5.0D0
!!$         PTPP(NS) = 5.0D0
!!$         PTS(NS)  = 0.05D0
!!$         PU(NS)   = 0.D0
!!$         PUS(NS)  = 0.D0
!!$         PNITB(NS)= 0.D0
!!$         PTITB(NS)= 0.D0
!!$         PUITB(NS)= 0.D0

!!$         ! *** helium ***
!!$         NS =
!!$
!!$         KIDNS(NS)= 'A'
!!$         IDION(NS)= 0.0D0
!!$         PA(NS)   = 4.0D0
!!$         PZ(NS)   = 2.0D0
!!$         PZ0(NS)  = 2.0D0
!!$         PN(NS)   = 1.0D0
!!$         PNS(NS)  = 0.0D0
!!$         PTPR(NS) = 5.0D0
!!$         PTPP(NS) = 5.0D0
!!$         PTS(NS)  = 0.05D0
!!$         PU(NS)   = 0.D0
!!$         PUS(NS)  = 0.D0
!!$         PNITB(NS)= 0.D0
!!$         PTITB(NS)= 0.D0
!!$         PUITB(NS)= 0.D0

         ! *** hydrogen (fast) ***
!!$         NS = 5
!!$
!!$         KIDNS(NS)= 'H'
!!$         IDION(NS)= 1.0D0
!!$         PA(NS)   = PA(2)
!!$         PZ(NS)   = PZ(2)
!!$         PZ0(NS)  = PZ0(2)
!!$         PN(NS)   = 0.0001D0
!!$         PNS(NS)  = 0.00005D0
!!$         PTPR(NS) = 50.D0
!!$         PTPP(NS) = 50.D0
!!$         PTS(NS)  = 10.D0
!!$         PU(NS)   = 0.D0
!!$         PUS(NS)  = 0.D0
!!$         PNITB(NS)= 0.D0
!!$         PTITB(NS)= 0.D0
!!$         PUITB(NS)= 0.D0

         ! *** deuterium (fast) ***
         NS = 3

         KIDNS(NS)= 'D'
         IDION(NS)= 1.0D0
         PA(NS)   = 2.0D0
         PZ(NS)   = 1.0D0
         PZ0(NS)  = 1.0D0
         PN(NS)   = 0.D0
         PNS(NS)  = 0.D0
         PTPR(NS) = 50.D0
         PTPP(NS) = 50.D0
         PTS(NS)  = 10.D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0

!!$         ! *** helium (fast,alpha) ***
!!$         NS = 5
!!$
!!$         KIDNS(NS)= 'A'
!!$         IDION(NS)= 1.0D0
!!$         PA(NS)   = PA(5)
!!$         PZ(NS)   = PZ(5)
!!$         PZ0(NS)  = PZ0(5)
!!$         PN(NS)   = 0.0D0
!!$         PNS(NS)  = 0.0D0
!!$         PTPR(NS) = 50.D0
!!$         PTPP(NS) = 50.D0
!!$         PTS(NS)  = 50.D0
!!$         PU(NS)   = 0.D0
!!$         PUS(NS)  = 0.D0
!!$         PNITB(NS)= 0.D0
!!$         PTITB(NS)= 0.D0
!!$         PUITB(NS)= 0.D0

         ! *** carbon ***
         NS = 4
         
         KIDNS(NS)= 'C'
         IDION(NS)= 0.0D0
         PA(NS)   = 12.d0
         PZ(NS)   = 6.d0
         PZ0(NS)  = 6.d0
         PN(NS)   = 0.D0
         PNS(NS)  = 0.D0
         PTPR(NS) = 1.D0
         PTPP(NS) = 1.D0
         PTS(NS)  = 0.1D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
=======
         PZCL(NS) = 0.D0
>>>>>>> 1.7

         ! *** dummy ***
      DO NS = 3, NSM
         KIDNS(NS)= ' '
         IDION(NS)= 0.0D0
         PA(NS)   = 1.0D0
         PZ(NS)   = 1.0D0
         PZ0(NS)  = 1.0D0
         PN(NS)   = 0.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.0D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PZCL(NS) = 0.D0
      ENDDO

      ENDIF

!     ======( PROFILE PARAMETERS )======


!        PROFN1: Density profile parameter (power of rho)
!        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
!        PROFT1: Temperature profile parameter (power of rho)
!        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
!        PROFU1: Rotation profile parameter (power of rho)
!        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))

      PROFN1= 2.D0
      PROFN2= 0.5D0
      PROFT1= 2.D0
      PROFT2= 1.D0
      PROFU1= 2.D0
      PROFU2= 1.D0

!     ======( MODEL PARAMETERS )======

!        MODELG: Control plasma geometry model
!                   0: Slab geometry
!                   1: Cylindrical geometry
!                   2: Toroidal geometry
!                   3: TASK/EQ output geometry
!                   4: VMEC output geometry
!                   5: EQDSK output geometry
!                   6: Boozer output geometry
!                   7: new VMEC output geometry
!                   8: call TOPICS/EQU
!                   9: call TASK/EQ
!        MODELN: Control plasma profile
!                   0: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; 0 in SOL
!                   1: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PNS in SOL
!                   7: Read from file by means of WMDPRF routine (DIII-D)
!                   8: Read from file by means of WMXPRF routine (JT-60)
!                   9: Read from file KNAMTR (TASK/TR)
!        MODELQ: Control safety factor profile (for MODELG=0,1,2)
!                   0: Parabolic q profile (Q0,QA,RHOMIN,RHOITB)
!                   1: Given current profile (RIP,PROFJ)

      MODELG= 2
      MODELN= 0
      MODELQ= 0

!        RHOMIN: rho at minimum q (0 for positive shear)
!        QMIN  : q minimum for reversed shear
!        RHOITB: rho at ITB (0 for no ITB)
!        RHOEDG: rho at EDGE for smoothing (1 for no smooth)

      RHOMIN = 0.D0
      QMIN   = 1.5D0
      RHOITB = 0.D0
      RHOEDG = 1.D0

!     ======( GRAPHIC PARAMETERS )======

!        RHOGMN: minimum rho in radial profile
!        RHOGMX: maximum rho in radial profile

      RHOGMN = 0.D0
      RHOGMX = 1.D0

!     ======( IO FILE NAMES )======

!        KNAMEQ: Filename of equilibrium data
!        KNAMWR: Filename of ray tracing data
!        KNAMWM: Filename of full wave data
!        KNAMFP: Filename of Fokker-Planck data
!        KNAMFO: Filename of File output
!        KNAMPF: Filename of profile data

      KNAMEQ = 'eqdata'
      KNAMWR = 'wrdata'
      KNAMWM = 'wmdata'
      KNAMFP = 'fpdata'
      KNAMFO = 'fodata'
      KNAMPF = 'pfdata'

!     ======( FILE IO MODES )======

!        MODEFR: File name interaction at reading file
!                 0 : WITHOUT PROMPT
!                 1 : WITH FILE NAME INPUT
!        MODEFW: File name interaction at writing file
!                 0 : WITHOUT PROMPT, ALWAYS OVERWRITE
!                 1 : WITHOUT PROMPT, CONFIRM, IF FILE EXISTS
!                 2 : WITHOUT PROMPT, ASK NEW NAME, IF FILE EXISTS
!                 3 : WITHOUT PROMPT, ERROR, IF FILE EXISTS
!                 4 : WITH FILE NAME INPUT, ALWAYS OVERWRITE
!                 5 : WITH FILE NAME INPUT, CONFIRM, IF FILE EXISTS
!                 6 : WITH FILE NAME INPUT, ASK NEW NAME, IF FILE EXISTS
!                 7 : WITH FILE NAME INPUT, ERROR, IF FILE EXISTS

!     ======( Default values )======

      MODEFR = 1
      MODEFW = 5

      IDEBUG = 0

      RETURN
    end subroutine pl_init

!     ****** INPUT PARAMETERS ******

    SUBROUTINE pl_parm(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

      IMPLICIT NONE
      INTEGER,INTENT(IN):: mode
      CHARACTER(LEN=*),INTENT(IN)::  kin
      INTEGER,INTENT(OUT):: ierr

    1 CALL TASK_PARM(MODE,'PL',KIN,pl_nlin,pl_plst,IERR)
      IF(IERR.NE.0) RETURN

      CALl pl_check(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
    END subroutine pl_parm

!     ****** INPUT NAMELIST ******

    SUBROUTINE pl_nlin(NID,IST,IERR)

      use plcomm, only:RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
                    NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
                    PNITB,PTITB,PUITB,PZCL, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    RHOMIN,QMIN,RHOITB,RHOEDG, &
                    MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
                    KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
                    MODEFR,MODEFW,IDEBUG, &
                    rkind,pl_allocate_ns

      implicit none
      integer,intent(in) :: NID
      integer,intent(out) :: IST,IERR
      integer,parameter:: NSM=100
      integer:: NS

      NAMELIST /PL/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
                    NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
                    MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
                    KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
                    MODEFR,MODEFW,IDEBUG

      READ(NID,PL,IOSTAT=IST,ERR=9800,END=9900)

      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
    END SUBROUTINE pl_nlin

!     ***** INPUT PARAMETER LIST *****

    SUBROUTINE pl_plst

      implicit none
      WRITE(6,601)
      RETURN

  601 FORMAT(' ','# &EQ : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'MODELG,MODELN,MODELQ,RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,IDEBUG'/ &
             9X,'MODEFW,MODEFR')
    END SUBROUTINE pl_plst

!     ****** CHECK INPUT PARAMETER ******

    SUBROUTINE pl_check(IERR)

      use plcomm
      implicit none
      integer:: IERR,NS
      real(rkind):: RHOE0,RHOES

      IERR=0

      IF(MODELG.NE.3) THEN
         IF((MODELG.LT.0).OR.(MODELG.GT.2)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELG: MODELG=',MODELG
            IERR=1
         ENDIF
         IF((MODELN.NE.0).AND.(MODELN.NE.9)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ELSE
         IF((MODELN.NE.0).AND.(MODELN.NE.1).AND.(MODELN.NE.2).AND. &
            (MODELN.NE.9)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1).AND.(MODELG.NE.3)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ENDIF

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

      RETURN
    END SUBROUTINE pl_check

!     ****** SHOW PARAMETERS ******

    SUBROUTINE pl_view

      use plcomm
      implicit none
      integer:: NS

      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    , &
                   'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'RKAP  ',RKAP  ,'RDLT  ',RDLT  , &
                   'Q0    ',Q0    ,'QA    ',QA
      WRITE(6,601) 'RIP   ',RIP   ,'PROFJ ',PROFJ , &
                   'PROFN1',PROFN1,'PROFN2',PROFN2
      WRITE(6,601) 'PROFT1',PROFT1,'PROFT2',PROFT2, &
                   'PROFU1',PROFU1,'PROFU2',PROFU2
      WRITE(6,601) 'RHOEDG',RHOEDG,'RHOGMN',RHOGMN, &
                   'RHOGMX',RHOGMX
      WRITE(6,604) 'MODELG',MODELG,'MODELN',MODELN, &
                   'MODELQ',MODELQ
      WRITE(6,604) 'MODEFR',MODEFR,'MODEFW',MODEFW

      WRITE(6,100)
      DO NS=1,NSMAX
         WRITE(6,110) NS,PA(NS),PZ(NS),PZ0(NS),PN(NS),PNS(NS)
      ENDDO
      WRITE(6,120)
      DO NS=1,NSMAX
         WRITE(6,130) NS,PTPR(NS),PTPP(NS),PTS(NS),PU(NS),PUS(NS)
      ENDDO
      IF(RHOITB.GT.0.D0) THEN
         WRITE(6,140)
         DO NS=1,NSMAX
           WRITE(6,150) NS,PNITB(NS),PTITB(NS),PUITB(NS),PZCL(NS)
         ENDDO
      ENDIF
      RETURN

  100 FORMAT(1H ,'NS    PA          PZ          PZ0         ', &
                 'PN          PNS')
  110 FORMAT(1H ,I2,' ',1P5E12.3)
  120 FORMAT(1H ,'NS    PTPR        PTPP        PTS         ', &
                 'PU          PUS')
  130 FORMAT(1H ,I2,' ',1P5E12.3)                               
  140 FORMAT(1H ,'NS    PNITB       PTITB       PUITB'      , &
                 'PZCL')
  150 FORMAT(1H ,I2,' ',1P4E12.4)                               
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
             2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  604 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
             2X,A6,'=',I7,4X  :2X,A6,'=',I7)
    END SUBROUTINE pl_view

  END MODULE plinit
