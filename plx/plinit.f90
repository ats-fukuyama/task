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
!        RHOITB: rho at ITB (0 for no ITB)
!        PNITB : Density increment at ITB              (1.0E20/Mm*3)
!        PTITB : Temperature increment at ITB                  (keV)
!        PUITB : Toroidal rotation velocity increment at ITB   (m/s)
!        PZCL  : normalized collision frequency

!        KID_NS: index of particle species
!        ID_NS : -1 : electron
!                 0 : neutral
!                 1 : ion
!                 2 : fast ion

      NSMAX = 2                  ! Default number of particle species

!     *** electron ***
         NS = 1

         KID_NS(NS)= ' e'
         ID_NS(NS) = -1
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
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PZCL(NS) = 0.D0

!     *** deuteron ***
         NS = 2

         KID_NS(NS)= ' D'
         ID_NS(NS) = 1
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
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0

!     *** triton ***
         NS = 3

         KID_NS(NS)= ' T'
         ID_NS(NS) = 1
         PA(NS)   = 3.0D0
         PZ(NS)   = 1.0D0
         PZ0(NS)  = 1.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0

!     *** Helium ion ***
         NS = 4

         KID_NS(NS)= 'He'
         ID_NS(NS) = 1
         PA(NS)   = 4.0D0
         PZ(NS)   = 2.0D0
         PZ0(NS)  = 2.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0

         ! *** dummy ***
      DO NS = 5, NSM
         KID_NS(NS)= ' H'
         ID_NS(NS)= 1
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
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PZCL(NS) = 0.D0
      ENDDO

!     *** PLANE  PARAMETERS ***
  r_corner(1)=0.D0
  r_corner(2)=RA
  r_corner(3)=0.D0
  z_corner(1)=0.D0
  z_corner(2)=0.D0
  z_corner(3)=2.D0*PI*RR

  br_corner(1:3)=0.D0
  bz_corner(1:3)=0.D0
  bt_corner(1:3)=BB

  DO ns=1,nsm
     pn_corner(1,ns)=pn(ns)
     pn_corner(2,ns)=pns(ns)
     pn_corner(3,ns)=pn(ns)
     ptpr_corner(1,ns)=ptpr(ns)
     ptpr_corner(2,ns)=pts(ns)
     ptpr_corner(3,ns)=ptpr(ns)
     ptpp_corner(1,ns)=ptpp(ns)
     ptpp_corner(2,ns)=pts(ns)
     ptpp_corner(3,ns)=ptpp(ns)
  END DO

!     ======( PROFILE PARAMETERS )======


!        PROFN1: Density profile parameter (power of rho)
!        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
!        PROFT1: Temperature profile parameter (power of rho)
!        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
!        PROFU1: Rotation profile parameter (power of rho)
!        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))

  DO NS=1,NSM
     PROFN1(NS)= 2.D0
     PROFN2(NS)= 0.5D0
     PROFT1(NS)= 2.D0
     PROFT2(NS)= 1.D0
     PROFU1(NS)= 2.D0
     PROFU2(NS)= 1.D0
  END DO

!     ======( MODEL PARAMETERS )======

!        MODELG: Control plasma geometry model
!              0: XYZ Slab geometry
!                 MODELB: 
!                      0: Translational mirror
!                      1: Translational mirror with curent rod (NCMAX.GT.0)
!                      2: XY 2D plane profile (defined at corners)
!                               (MODELGX: 0:linear, 1,2:parabolic)
!                      3: XY 2D plane Read 2D mag file 
!              1: RZphi Cylindrical geometry
!                 MODEL B:
!                      0: Axisymmetric mirror
!                      1: Axisymmetric mirror with circular coils (NCMAX.GT.0)
!                      2: RZ 2D toroidal profile (defined at corners)
!                          (MODELGX: 0:linear,1,2:parabolic)
!                      3: RZ 2D toroidal Read 2D mag file 
!              2: RZphi Toroidal geometry
!              3: RZphi Read TASK/EQ output geometry
!              4: RZphi Read VMEC output geometry
!              5: RZphi Read EQDSK output geometry
!              6: RZphi Read Boozer output geometry
!              7: RZphi Read new VMEC output geometry
!              8: RZphi call TASK/EQU
!              9: RZphi call TASK/EQ
!             10: reserved for GAMMA-10
!             11: Straight helical geometry
!        MODELN: Control plasma profile
!                   0: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; 0 in SOL
!                   1: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PNS in SOL
!                   7: Read from file through WMDPRF (DIII-D)
!                   8: Read from file through WMXPRF (JT-60)
!                   9: Read from file KNAMTR (TASK/TR)
!                  12: Read from 2D nT file
!                  14: Read from 2D nT file
!        MODELQ: Control safety factor profile (for MODELG=1,2)
!                   0: Parabolic q profile (Q0,QA,RHOMIN)
!                   1: Given current profile (RIP,PROFJ)
!        MODEL_PROF: profile parameter
!                   0: PROFX1(NS)=PROFX1(1),PROFX2(NS)=PROFX2(1): compatibility
!                   1: PROFX1(NS),PROFX2(NS): defined separately
!        MODEL_NPROF: neutral profile parameter
!                   0: Flat profile
!                   1: Flat only in plasma, 0 outside
!                   2: (1-psi) dependence, 0 outside

      MODELG= 2
      MODELB= 0
      MODELN= 0
      MODELQ= 0
      MODEL_PROF=0
      MODEL_NPROF=0

!        RHOMIN: rho at minimum q (0 for positive shear)
!        QMIN  : q minimum for reversed shear
!        RHOITB: rho at ITB (0 for no ITB)
!        RHOEDG: rho at EDGE for smoothing (1 for no smooth)

      RHOMIN = 0.D0
      QMIN   = 1.5D0
      RHOEDG = 1.D0

!        PPN0: Neutral pressure [Pa] 1 Torr = 1 mmHg = 133.322 Pa
!        PTN0: Neutral temperature [eV]
!        RF_PL: wave frequency [MHz], usually set in wave code

      PPN0 = 3.0D0
      PTN0 = 0.03D0
      RF_PL = 1.D6    ! 1THz to reduce pzcl

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
                    RHOITB,PNITB,PTITB,PUITB,PZCL, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    r_corner,z_corner, &
                    br_corner,bz_corner,bt_corner, &
                    pn_corner,ptpr_corner,ptpp_corner, &
                    RHOMIN,QMIN,RHOEDG,PPN0,PTN0,RF_PL, &
                    MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
                    RHOGMN,RHOGMX, &
                    KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
                    MODEFR,MODEFW,IDEBUG, &
                    rkind,pl_allocate_ns

      implicit none
      integer,intent(in) :: NID
      integer,intent(out) :: IST,IERR
      integer,parameter:: NSM=100
      INTEGER:: NS

      NAMELIST /PL/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
                    NSMAX,PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL, &
                    r_corner,z_corner, &
                    br_corner,bz_corner,bt_corner, &
                    pn_corner,ptpr_corner,ptpp_corner, &
                    PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                    RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
                    PPN0,PTN0,RF_PL, &
                    MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
                    RHOGMN,RHOGMX, &
                    KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
                    MODEFR,MODEFW,IDEBUG

      READ(NID,PL,IOSTAT=IST,ERR=9800,END=9900)

      IF(MODEL_PROF.EQ.0) THEN
         DO NS=2,NSMAX
            PROFN1(NS)=PROFN1(1)
            PROFN2(NS)=PROFN2(1)
            PROFT1(NS)=PROFT1(1)
            PROFT2(NS)=PROFT2(1)
            PROFU1(NS)=PROFU1(1)
            PROFU2(NS)=PROFU2(1)
         END DO
      END IF

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

  601 FORMAT(' ','# &PL : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'r_corner,z_corner,br_corner,bz_corner,bt_corner,'/ &
             9X,'pn_corner,ptpr_corner,ptpp_corner,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'PPN0,PTN0,RFCL,'/ &
             9X,'MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG')
    END SUBROUTINE pl_plst

!     ****** CHECK INPUT PARAMETER ******

    SUBROUTINE pl_check(IERR)

      use plcomm
      implicit none
      integer:: IERR,NS
      real(rkind):: RHOE0,RHOES

      IERR=0

      IF(MODELG.NE.3) THEN
         IF((MODELG.LT.0).OR.(MODELG.GT.11)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELG: MODELG=',MODELG
            IERR=1
         ENDIF
         IF((MODELN.LT.0).OR.(MODELN.GT.14)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1)) THEN
            WRITE(6,*) 'XX pl_check: INVALID MODELQ: MODELQ=',MODELQ
            IERR=1
         ENDIF
      ELSE
         IF((MODELN.LT.0).OR.(MODELN.GT.14)) THEN
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
      integer:: NS,i

      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    , &
                   'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'RKAP  ',RKAP  ,'RDLT  ',RDLT  , &
                   'Q0    ',Q0    ,'QA    ',QA
      WRITE(6,601) 'RIP   ',RIP   ,'PROFJ ',PROFJ
      WRITE(6,601) 'RHOEDG',RHOEDG,'RHOGMN',RHOGMN, &
                   'RHOGMX',RHOGMX
      WRITE(6,604) 'MODELG',MODELG,'MODELB',MODELB, &
                   'MODELN',MODELN,'MODELQ',MODELQ
      WRITE(6,604) 'MODEFR',MODEFR,'MODEFW',MODEFW
      WRITE(6,'(A,I5)') 'MODEL_PROF  =',MODEL_PROF
      WRITE(6,'(A,I5)') 'MODEL_NPROF =',MODEL_NPROF

      WRITE(6,100)
      DO NS=1,NSMAX
         WRITE(6,110) NS,PA(NS),PZ(NS),PZ0(NS),PN(NS),PNS(NS)
      ENDDO
      WRITE(6,120)
      DO NS=1,NSMAX
         WRITE(6,130) NS,PTPR(NS),PTPP(NS),PTS(NS),PU(NS),PUS(NS)
      ENDDO

      WRITE(6,140)
      DO NS=1,NSMAX
         WRITE(6,150) NS,RHOITB(NS),PNITB(NS),PTITB(NS),PUITB(NS),PZCL(NS)
      END DO

      WRITE(6,160)
      DO NS=1,NSMAX
         WRITE(6,170) NS,PROFN1(NS),PROFN2(NS),PROFT1(NS),PROFT2(NS), &
                         PROFU1(NS),PROFU2(NS)
      END DO

      WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  , &
                   'RF_PL ',RF_PL

      IF(MODELG.EQ.11.OR.MODELG.EQ.13) THEN
         WRITE(6,606) 'r_corner:  ',(r_corner(i),i=1,3)
         WRITE(6,606) 'z_corner:  ',(z_corner(i),i=1,3)
         WRITE(6,606) 'br_corner: ',(br_corner(i),i=1,3)
         WRITE(6,606) 'bz_corner: ',(bz_corner(i),i=1,3)
         WRITE(6,606) 'bt_corner: ',(bt_corner(i),i=1,3)
         DO ns=1,NSMAX
            WRITE(6,'(A,I3)') 'ns=',ns
            WRITE(6,606) 'pn_corner: ',(pn_corner(i,ns),i=1,3)
            WRITE(6,606) 'ptpr_corner:',(ptpr_corner(i,ns),i=1,3)
            WRITE(6,606) 'ptpp_corner:',(ptpp_corner(i,ns),i=1,3)
         END DO
      END IF

      RETURN

  100 FORMAT(' ','NS    PA          PZ          PZ0         ', &
                       'PN          PNS')
  110 FORMAT(' ',I2,' ',1P5E12.4)
  120 FORMAT(' ','NS    PTPR        PTPP        PTS         ', &
                       'PU          PUS')
  130 FORMAT(' ',I2,' ',1P5E12.4)                               
  140 FORMAT(' ','NS    RHOITB      PNITB       PTITB       ', &
                       'PUITB       PZCL')
  150 FORMAT(' ',I2,' ',1P5E12.4)                               
  160 FORMAT(' ','NS    PROFN1      PROFN2      PROFT1      ', &
                       'PROFT2      PROFU1      PROFU2')
  170 FORMAT(' ',I2,' ',1P6E12.4)                               
  601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
             2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
             2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  606 FORMAT(' ',A12,1P3E12.4)
    END SUBROUTINE pl_view

  END MODULE plinit
