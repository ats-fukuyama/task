C     $Id$
C
C     *********** INPUT PARAMETER FROM NAMELIST /WR/ ***********
C
C     RF     : WAVE FREQUENCY FOR RAY TRACING [MHZ]
C     RPI    : INITIAL MAJOR RADIUS R [M]
C     ZPI    : INITIAL VERTICAL POSITION Z [M]
C     PHII   : INITIAL TOROIDAL ANGLE [RADIAN]
C     RNZI   : INITIAL VERTICAL REFRACTIVE INDEX
C     RNPHII : INITIAL TOROIDAL REFRACTIVE INDEX
C     RKR0   : SPECULATED INITIAL RADIAL WAVE NUMBER [1/M]
C     UUI    : INITIAL WAVE ENERGY
C
C     SMAX   : MAXIMUM RAY LENGTH
C     DELS   : INCREMENTAL LENGTH OF RAY
C     UUMIN  : MINIMUM POWER TO ADVANCE RAY
C     NRAYMX : NUMBER OF RAYS
C
C     EPSRAY : CONVERGENCE CRITEIRION IN RAY TRACING
C     DELRAY : MINIMUM STEP SIZE IN RAY TRACING
C     DELDER : STEP SIZE TO CALCULATE DERIVATIVES IN RAY TRACING
C
C     DELKR  : STEP SIZE TO ESTIMATE D/DKR IN NEWTON METHOD
C     EPSNW  : CONVERGENCE CRITEIRION IN NEWTON METHOD
C     LMAXNW : MAXIMUM ITERATION COUNT IN NEWTON METHOD
C
C     NRZMAX : NUMBER OF RADIAL DIVISION FOR ABSORBED POWER
C
C     INTYPE : INPUT TYPE OF WAVE PARAMETERS
C              0 : RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU
C              1 : RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU
C              2 : RF,RP,ZP,PHI,MODE,ANGZ,ANGPH,UU
C
C     RCURVA : INITIAL WAVE-FRONT CURVATURE RADIUS (0 for Plane wave)
C     RCURVB : INITIAL WAVE-FRONT CURVATURE RADIUS (0 for Plane wave)
C                 RCURVA perp to k and B
C                 RCURVB perp to k and in kxB plane
C     RBRADA  : INITIAL BEAM RADIUS
C     RBRADB  : INITIAL BEAM RADIUS
C                 RBRADA perp to k and B
C                 RBRADB perp to k and in kxB plane
C     NRADMX  : NUMBER OF RADIAL DIVISION IN BEAM TRACING
C
C     ****** INITIALIZE INPUT PARAMETERS ******
C
      SUBROUTINE WRINIT
C
      INCLUDE 'wrcomm.h'
C
      RF     = 170.D3
      RPI    = 3.95D0
      ZPI    = 0.D0
      PHII   = 0.D0
      RNZI   = 0.D0
      RNPHII = 0.5D0
      RKR0   = -1000.D0
      UUI    = 1.D0
C
      SMAX   = 1.0D0
      DELS   = 0.05D0
      UUMIN  = 1.D-4
      NRAYMX = 1
C
      DELKR  = 1.D0
      EPSNW  = 1.D-6
      LMAXNW = 100
C
      EPSRAY = 1.D-4
      DELRAY = 1.D-3
      DELDER = 1.D-4
C
      NRZMAX = 1000
      INTYPE = 0
C
      RCURVA = 0.D0
      RCURVB = 0.D0
      RBRADA = 0.03D0
      RBRADB = 0.03D0
      NRADMX= 100
C
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE WRPARM
C
      INCLUDE 'wrcomm.h'
C
      LOGICAL LEX
      CHARACTER KPNAME*32
      NAMELIST /WR/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
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
     &              MODELV,INTYPE,
     &              RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,UUI,
     &              SMAX,DELS,UUMIN,EPSNW,DELKR,EPSRAY,DELRAY,DELDER,
     &              LMAXNW,NRZMAX,NRAYMX,KNAMWR,
     &              RCURVA,RCURVB,RBRADA,RBRADB,NRADMX
      DATA INITEQ,INITFP/0,0/
C
    1 WRITE(6,*) '## INPUT : &WR'
      READ(5,WR,ERR=1,END=9000)
C
 3000 IERR=0
      IF(MODELG.EQ.3) THEN
         IF(INITEQ.EQ.0) THEN
            CALL EQLOAD(1,KNAMEQ,IERR)
            IF(IERR.EQ.0) THEN
               CALL EQSETP
               CALL EQPSIC(51,32,64,IERR)
               CALL EQGETB(BB,RR,RIP,RA,RKAP,RDEL,RB)
            ENDIF
            INITEQ=1
         ENDIF
      ELSE
         INITEQ=0
      ENDIF
      IF(IERR.EQ.1) RETURN
C
      IF(MODELV.EQ.1) THEN
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
      ENTRY WRPARF
C
      KPNAME='wrparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(LEX) THEN
         OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
         READ(25,WR,ERR=9800,END=9900)
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
      SUBROUTINE WRVIEW
C
      INCLUDE 'wrcomm.h'
C
      WRITE(6,601) 'SMAX  ',SMAX  ,'DELS  ',DELS  ,
     &             'UUMIN ',UUMIN
      WRITE(6,601) 'EPSNW ',EPSNW ,'DELKR ',DELKR
      WRITE(6,601) 'EPSRAY',EPSRAY,'DELRAY',DELRAY,
     &             'DELDER',DELDER
      WRITE(6,602) 'NRAYMX',NRAYMX,'LMAXNW',LMAXNW,
     &             'NRZMAX',NRZMAX,'NRADMX',NRADMX
      WRITE(6,602) 'INTYPE',INTYPE
      RETURN
C
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END
