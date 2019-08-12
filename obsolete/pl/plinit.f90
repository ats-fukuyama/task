!     ****** INITIALIZE INPUT PARAMETERS ******

      SUBROUTINE PLINIT

!      INCLUDE '../pl/plcomm.inc'
!      INCLUDE '../pl/plcnst.inc'
      USE PLCNS1, ONLY : AME, AMP
      USE PLDAT0, ONLY : NRMAXPL, NSMAXPL
      USE PLDBG1
      USE PLNAM1
      USE PLPRG1
      USE PLPRM1
      USE PLPRM2
      USE PLPRM3
      USE PLPRM4
      USE PLPRM5
      USE PLPRM6
      USE PLPRM7
      USE PLPRM8
      USE PLPRM9

      IMPLICIT NONE
      INTEGER(4)  :: NS

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
!        PN    : Density at center                     (1.0E20/m**3)
!        PNS   : Density on plasma surface             (1.0E20/m**3)
!        PZCL  : Ratio of collision frequency to wave frequency
!        PTPR  : Parallel temperature at center                (keV)
!        PTPP  : Perpendicular temperature at center           (keV)
!        PTS   : Temperature on surface                        (keV)
!        PU    : Toroidal rotation velocity at center          (m/s)
!        PUS   : Toroidal rotation velocity on surface         (m/s)
!        PNITB : Density increment at ITB              (1.0E20/Mm*3)
!        PTITB : Temperature increment at ITB                  (keV)
!        PUITB : Toroidal rotation velocity increment at ITB   (m/s)

      NSMAX = MIN(2,NSM)

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

!     ======( MODEL PARAMETERS )======

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

      NRMAXPL= 100
      NSMAXPL= NSMAX

      IDEBUG = 0

      RETURN
      END

!     ****** INPUT PARAMETERS ******

      SUBROUTINE PLPARM(MODE,KIN,IERR)

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
      EXTERNAL PLNLIN,PLPLST
      INTEGER(4),      INTENT(IN)  :: MODE
      CHARACTER(LEN=*),INTENT(IN):: KIN
      INTEGER(4),      INTENT(OUT) :: IERR
!      CHARACTER KIN*(*)

    1 CALL TASK_PARM(MODE,'PL',KIN,PLNLIN,PLPLST,IERR)
      IF(IERR.NE.0) RETURN

      CALl PLCHEK(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
      END

!     ****** INPUT NAMELIST ******

      SUBROUTINE PLNLIN(NID,IST,IERR)

!      INCLUDE '../pl/plcomm.inc'
      USE PLDBG1
      USE PLNAM1, ONLY : KNAMEQ, KNAMFO, KNAMFP, KNAMWR
      USE PLPRG1
      USE PLPRM1
      USE PLPRM2
      USE PLPRM3
      USE PLPRM4
      USE PLPRM5
      USE PLPRM6
      USE PLPRM7
      USE PLPRM8
      USE PLPRM9
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)  :: NID
      INTEGER(4), INTENT(OUT) :: IST, IERR

      NAMELIST /PL/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS, &
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
     &              MODELG,MODELN,MODELQ,RHOGMN,RHOGMX, &
     &              KNAMEQ,KNAMWR,KNAMFP,KNAMFO,IDEBUG

      READ(NID,PL,IOSTAT=IST,ERR=9800,END=9900)
      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
      END

!     ***** INPUT PARAMETER LIST *****

      SUBROUTINE PLPLST

      IMPLICIT NONE
      WRITE(6,601)
      RETURN

  601 FORMAT(' ','# &EQ : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
     &       9X,'NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,'/   &
     &       9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/      &
     &       9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/    &
     &       9X,'MODELG,MODELN,MODELQ,RHOGMN,RHOGMX,'/             &
     &       9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,IDEBUG')
      END

!     ****** CHECK INPUT PARAMETER ******

      SUBROUTINE PLCHEK(IERR)

!      INCLUDE '../pl/plcomm.inc'
      USE PLCOM0
      USE PLPRM3
      USE PLPRM4, ONLY : PN, PNS, PZ
      USE PLPRM9
      IMPLICIT NONE
      INTEGER(4), INTENT(OUT) :: IERR
      REAL(8)     :: RHOE0, RHOES
      INTEGER(4)  :: NS

      IERR=0

      IF(NSMAX.GT.NSM) THEN
         WRITE(6,*) 'XX PLPARM: NSMAX.GT.NSM: NSMAX,NSM=',NSMAX,NSM
         IERR=1
      ENDIF

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
         IF((MODELN.NE.0).AND.(MODELN.NE.1).AND.(MODELN.NE.2).AND.(MODELN.NE.9)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELN: MODELN=',MODELN
            IERR=1
         ENDIF
         IF((MODELQ.NE.0).AND.(MODELQ.NE.1).AND.(MODELG.NE.3)) THEN
            WRITE(6,*) 'XX PLPARM: INVALID MODELQ: MODELQ=',MODELQ
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
      END

!     ****** SHOW PARAMETERS ******

      SUBROUTINE PLVIEW

!      INCLUDE '../pl/plcomm.inc'
      USE PLPRG1
      USE PLPRM1
      USE PLPRM2
      USE PLPRM3
      USE PLPRM4
      USE PLPRM5
      USE PLPRM6
      USE PLPRM7, ONLY : RHOEDG
      USE PLPRM8
      USE PLPRM9

      IMPLICIT NONE
      INTEGER(4)  :: NS


      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    ,'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'RKAP  ',RKAP  ,'RDLT  ',RDLT  ,'Q0    ',Q0    ,'QA    ',QA
      WRITE(6,601) 'RIP   ',RIP   ,'PROFJ ',PROFJ ,'PROFN1',PROFN1,'PROFN2',PROFN2
      WRITE(6,601) 'PROFT1',PROFT1,'PROFT2',PROFT2,'PROFU1',PROFU1,'PROFU2',PROFU2
      WRITE(6,601) 'RHOEDG',RHOEDG,'RHOGMN',RHOGMN,'RHOGMX',RHOGMX
      WRITE(6,604) 'MODELG',MODELG,'MODELN',MODELN,'MODELQ',MODELQ

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

  100 FORMAT(1H ,'NS    PA          PZ          PN          PNS         PZCL')
  110 FORMAT(1H ,I2,' ',1P5E12.4)
  120 FORMAT(1H ,'NS    PTPR        PTPP        PTS         PU          PUS')
  130 FORMAT(1H ,I2,' ',1P5E12.4)
  140 FORMAT(1H ,'NS    PNITB       PTITB       PUITB')
  150 FORMAT(1H ,I2,' ',1P3E12.4)
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  604 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END
