C     $Id$
C
C     *********** INPUT PARAMETER FROM NAMELIST /DP/ ***********
C
C     MODELP: TYPE OF DIELECTRIC TENSOR
C                 0 : COLD PLASMA WITH COLLISION
C                 1 : WARM PLASMA WITH COLLISION
C                 2 : HOT PLASMA, HARMONICS NUMBER -1,0,+1
C                 3 : HOT PLASMA, HARMONICS NUMBER -2,-1,0,+1,2
C                 4 : HOT PLASMA, HARMONICS NUMBER FROM NDISP1 TO NDISP2
C                 5 : HOT PLASMA, RELATIVISTIC (TO BE IMPLEMENTED)
C                 6 : HOT PLASMA, FM(P)
C                 7 : HOT PLASMA, RM(P), WITH RELATIVISTIC EFFECT
C                 8 : HOT PLASMA, F(P),  
C                 9 : HOT PLASMA, F(P), WITH RELATIVISTIC EFFECT
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
C     MODELF : 0 : ANALYTIC MAXWELLIAN DISTRIBUTION
C              1 : READ FPDATA DISTRIBUTION
C
C
C     ****** INITIALIZE INPUT PARAMETERS ******
C
      SUBROUTINE DPINIT
C
      INCLUDE 'dpcomm.h'
C
         MODELP(1)= 5
         NDISP1(1)=-2
         NDISP2(1)= 2
C
      IF(NSM.GE.2) THEN
         MODELP(2)=0
         NDISP1(2)=-3
         NDISP2(2)= 3
      ENDIF
C
      DO NS=3,NSM
         MODELP(NS)=0
         NDISP1(NS)=-2
         NDISP2(NS)= 2
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
      SUBROUTINE DPPARM
C
      INCLUDE 'dpcomm.h'
C
      LOGICAL LEX
      CHARACTER KPNAME*32
      NAMELIST /DP/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,
     &              MODELG,MODELN,MODELQ,
     &              KNAMEQ,KNAMWR,KNAMFP,IDEBUG,
     &              MODELP,NDISP1,NDISP2,
     &              RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,
     &              RF1,RFI1,RKX1,RKY1,RKZ1,RX1,
     &              RF2,RFI2,RKX2,RKY2,RKZ2,RX2,
     &              NXMAX,EPSRT,LMAXRT,
     &              MODELF
      DATA INITEQ,INITFP/0,0/
C
    1 WRITE(6,*) '## INPUT : &DP'
      READ(5,DP,ERR=1,END=9000)
C
 3000 IF(MODELG.EQ.3) THEN
         IF(INITEQ.EQ.0) THEN
            CALL EQLOAD(MODELG,KNAMEQ,IERR)
            IF(IERR.EQ.0) THEN
               CALL EQSETP
               CALL EQPSIC(51,32,64)
               CALL EQGETB(BB,RR,RIP,RA,RKAP,RDEL,RB)
            ENDIF
            INITEQ=1
         ENDIF
      ELSE
         INITEQ=0
      ENDIF
C
      IF(MODELF.EQ.1) THEN
         IF(INITFP.EQ.0) THEN
            CALL DPLDFP
            INITFP=1
         ENDIF
      ELSE
         INITFP=0
      ENDIF
C
 9000 RETURN
C
C
      ENTRY DPPARF
C
      KPNAME='dpparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(LEX) THEN
         OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
         READ(25,DP,ERR=9800,END=9900)
         CLOSE(25)
         WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
      ENDIF
      GOTO 3000
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
      SUBROUTINE DPVIEW
C
      INCLUDE 'dpcomm.h'
C
      WRITE(6,100)
      DO I=1,NSMAX
        WRITE(6,110) I,MODELP(I),NDISP1(I),NDISP2(I)
      ENDDO
C
      WRITE(6,602) 'MODELF',MODELF
      RETURN
C
  100 FORMAT(1H ,'NS    MODELP  NDISP1  NDISP2')
  110 FORMAT(1H ,I2,' ',3I8)                               
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END


