C     $Id$
C
C     *********** INPUT PARAMETER FROM NAMELIST /DP/ ***********
C
C     MODELP: TYPE OF ANALYTIC DIELECTRIC TENSOR
C                 0 : COLLISIONLESS COLD MODEL
C                 1 : COLLISIONAL COLD MODEL
C                 2 : IDEAL MHD MODEL
C                 3 : RESISTIVE MHD MODEL
C                 4 : KINETIC MODEL WITHOUT FLR
C                 5 : KINETIC MODEL WITH FLR
C                 6 : KINETIC MODEL WITH RELATIVISTIC EFFECTS (test)
C                 7 : GYROKINETIC MODEL (coming)
C                 8 : GYROKINETIC MODEL (coming)
C                 9 : LOCAL MODEL (MODELP locally specified by MODELPR)
C                -1 : (WM) MHD plasma
C                -2 : (WM) Cold plasma
C                -3 : (WM) Hot plasma (No FLR)
C                -4 : (WM) Hot plasma (Cold FLR)
C                -5 : (WM) Hot plasma (FLR)
C
C              0- 9 : PROPAGATION  = GIVEN MODEL
C                     POLARIZATION = GIVEN MODEL
C                     ABSORPTION   = GIVEN MODEL
C
C             10-19 : PROPAGATION  = COLD
C                     POLARIZATION = GIVEN MODEL
C                     ABSORPTION   = GIVEN MODEL
C
C             20-29 : PROPAGATION  = COLD
C                     POLARIZATION = COLD
C                     ABSORPTION   = GIVEN MODEL
C
C             30-39 : PROPAGATION  = KINETIC
C                     POLARIZATION = KINETIC
C                     ABSORPTION   = GIVEN MODEL
C
C     MODELV : NUMERICAL MODEL (*: not yet implemented)
C              0 : ANALYTIC MODEL
C              1 : KINETIC: ANALYTIC MAXWELLIAN DISTRIBUTION
C              2 : KINETIC: READ FPDATA DISTRIBUTION
C              3 : KINETIC: ANALYTIC MAXWELLIAN DISTRIBUTUION (RELATIVISTIC)
C              4 : KINETIC: READ FPDATA DISTRIBUTION (RELATIVISTIC)
C              5*: GYROKINETIC: ANALYTIC MAXWELLIAN DISTRIBUTION
C              6*: GYROKINETIC: READ FPDATA DISTRIBUTION
C              7*: GYROKINETIC: ANALYTIC MAXWELLIAN DISTRIBUTUION (REL.)
C              8*: GYROKINETIC: READ FPDATA DISTRIBUTION (REL.)
C              9 : LOCAL MODEL (MODELV locally specified by MODELVR)
C
C     NDISP1: MINIMUM HARMONIC NUMBER (VALID ONLY FOR MODELP>=5)
C     NDISP2: MAXMUM  HARMONIC NUMBER (VALID ONLY FOR MODELP>=5)
C
C     RF0,RFI0,RKX0,RKY0,RKZ0 : STANDARD PARAMETER FOR ROOT FINDING
C     RF1,RF2                 : SCAN RANGE OF REAL FREQUENCY (MHZ)
C     RFI1,RFI2               : SCAN RANGE OF IMAGINARY FREQUENCY (MHZ)
C     RKX1,RKX2               : SCAN RANGE OF WAVE NUMBER KX (1/M)
C     RKY1,RKY2               : SCAN RANGE OF WAVE NUMBER KY (1/M)
C     RKZ1,RKZ2               : SCAN RANGE OF WAVE NUMBER KZ (1/M)
C     X1,X2 Å@Å@              : SCAN RANGE OF POSITION X (M)
C
C     NXMAX  : NUMBER OF SCAN POINTS
C     EPSRT  : CONVERGENCE CRITERION OF ROOT FINDING
C     LMAXRT : MAXIMUM ITERATION COUNT OF ROOT FINDING
C
C     ****** INITIALIZE INPUT PARAMETERS ******
C
      SUBROUTINE DPINIT
C
      INCLUDE 'dpcomm.inc'
C
         MODELP(1)= 5
         NDISP1(1)=-2
         NDISP2(1)= 2
         modelv(1)= 0
C
      IF(NSM.GE.2) THEN
         MODELP(2)= 0
         NDISP1(2)=-2
         NDISP2(2)= 2
         MODELV(2)= 0
      ENDIF
C
      DO NS=3,NSM
         MODELP(NS)= 0
         NDISP1(NS)=-2
         NDISP2(NS)= 2
         modelv(NS)= 0
      ENDDO
C
      RF0    = 160.D3
      RFI0   =   0.D0
      RKX0   = 800.D0
      RKY0   = 160.D0
      RKZ0   =   0.D0
      RX0    = RR
      RY0    = 0.D0
      RZ0    = 0.D0
C
      RF1    = 80000.D0
      RF2    = 16000.D0
      RFI1   =-5.D0
      RFI2   = 5.D0
      RKX1   = 0.D0
      RKX2   = 1600.D0
      RKY1   = 0.D0
      RKY2   = 1600.D0
      RKZ1   = 0.D0
      RKZ2   = 1600.D0
      RX1    = RR+0.D0
      RX2    = RR+0.5D0
C
      NXMAX  = 21
      EPSRT  = 1.D-8
      LMAXRT = 20
C
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE DPPARM(MODE,KIN,IERR)
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
      EXTERNAL DPNLIN,DPPLST
      CHARACTER KIN*(*)
C
    1 CALL TASK_PARM(MODE,'DP',KIN,DPNLIN,DPPLST,IERR)
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
      SUBROUTINE DPNLIN(NID,IST,IERR)
C
      INCLUDE 'dpcomm.inc'
C
C
      NAMELIST /DP/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,
     &              MODELG,MODELN,MODELQ,
     &              KNAMEQ,KNAMWR,KNAMFP,MODEFR,MODEFW,IDEBUG,
     &              MODELP,NDISP1,NDISP2,
     &              RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,
     &              RF1,RFI1,RKX1,RKY1,RKZ1,RX1,
     &              RF2,RFI2,RKX2,RKY2,RKZ2,RX2,
     &              NXMAX,EPSRT,LMAXRT,
     &              MODELV
C
      READ(NID,DP,IOSTAT=IST,ERR=9800,END=9900)
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
      SUBROUTINE DPPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &DP : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/
     &       9X,'NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,'/
     &       9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/
     &       9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,'/
     &       9X,'MODELG,MODELN,MODELQ,'/
     &       9X,'KNAMEQ,KNAMWR,KNAMFP,IDEBUG,MODEFR,MODEFW,'/
     &       9X,'MODELP,NDISP1,NDISP2,'/
     &       9X,'RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,'/
     &       9X,'RF1,RFI1,RKX1,RKY1,RKZ1,RX1,'/
     &       9X,'RF2,RFI2,RKX2,RKY2,RKZ2,RX2,'/
     &       9X,'NXMAX,EPSRT,LMAXRT,'/
     &       9X,'MODELV')
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE DPCHEK(IERR)
C
      INCLUDE 'dpcomm.inc'

      DATA INITFP/0/
C
      DO NS=1,NSMAX
         IF((MODELV(NS).EQ.2.OR.MODELV(NS).EQ.4)) THEN
            IF(INITFP.EQ.0) THEN
               write(6,*) '----- DPLDFP ----- NS=',NS
               CALL DPLDFP
               INITFP=1
            ENDIF
         ELSE
            RHON_MIN=0.D0
            RHON_MAX=1.D0
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE DPVIEW
C
      INCLUDE 'dpcomm.inc'
C
      WRITE(6,100)
      DO NS=1,NSMAX
        WRITE(6,110) NS,MODELP(NS),MODELV(NS),NDISP1(NS),NDISP2(NS)
      ENDDO
C
      RETURN
C
  100 FORMAT(1H ,'NS    MODELP  MODELV  NDISP1  NDISP2')
  110 FORMAT(1H ,I2,' ',4I8)                               
C  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
C     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END
