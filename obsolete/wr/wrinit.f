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
C     MODEW  : 0: RKR0, 1:RKR0_1, 2:RKR0_2
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
C     MDLWRI : INPUT TYPE OF WAVE PARAMETERS
C              0 : RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU
C              1 : RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU
C             11 : RF,RP,ZP,RKR0,RNZ,RNPHI,UU
C            100 : RFIN,RPIN,ZPIN,PHIIN,RKRIN,RNZIN,RNPHIIN,UUIN: namelist
C            101 : RFIN,RPIN,ZPIN,PHIIN,RKRIN,ANGZIN,ANGPHIN,UUIN: namelist
C
C     MDLWRG : TYPE OF GRAPHICS
C              0 : FULL TORUS, FULL RADIUS FOR DEPOSITION
C              1 : PARTIAL TORUS, FULL RADIUS FOR DEPOSITION
C              2 : FULL TORUS, PARTIAL RADIUS FOR DEPOSITION
C              3 : PARTIAL TORUS, PARTIAL RADIUS FOR DEPOSITION
C             11 : 2D plane
C
C     MDLWRQ : TYPE OF DIFFERENTIAL EQUATION IN RAY TRACING
C              0 : RUNGE-KUTTA, FIXED STEPSIZE
C              1 : RUNGE-KUTTA, FIXED STEPSIZE, with k_X adjust to satisfy D=0
C              2 : RUNGE-KUTTA, FIXED STEPSIZE, with mode conversion
C              3 : RUNGE-KUTTA, VARIABLE STEPSIZE
C              4 : SYMPLECTIC METHOD, FIXED STEPSIZE (not completed)
C
C     MDLWRW : Level of PRINT OUTPUT
C              0 : NO output
C             -1 : Write initial kr calculation
C              1 : Write data every step
C              2 : Write data every 10 steps
C              3 : Write data every 100 step
C              4 : Write data every 1000 step
C              5 : Write data every 10000 step
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
      USE plcomm
      INCLUDE 'wrcomm.inc'
C
      RF     = 170.D3
      RPI    = 3.95D0
      ZPI    = 0.D0
      PHII   = 0.D0
      RNZI   = 0.D0
      RNPHII = 0.5D0
      RKR0   = -1000.D0
      UUI    = 1.D0
      MODEW  = 0
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
      NRZMAX = 100
      MDLWRI = 0
      MDLWRG = 0
      MDLWRQ = 1
      MDLWRW = 0
C
      RCURVA = 0.D0
      RCURVB = 0.D0
      RBRADA = 0.03D0
      RBRADB = 0.03D0
      NRADMX = 100
C
      DO NRAY=1,NRAYM
         RFIN(NRAY)     = 170.D3
         RPIN(NRAY)    = 3.95D0
         ZPIN(NRAY)    = 0.D0
         PHIIN(NRAY)   = 0.D0
         RNZIN(NRAY)   = 0.D0
         RNPHIIN(NRAY) = 0.5D0
         RKRIN(NRAY)   = -1000.D0
         UUIN(NRAY)    = 1.D0
         MODEWIN(NRAY) = 0
         RCURVAIN(NRAY)= 0.D0
         RCURVBIN(NRAY)= 0.D0
         RBRADAIN(NRAY)= 0.03D0
         RBRADBIN(NRAY)= 0.03D0
      ENDDO

      MODEFW=0
      MODEFR=0
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE WRPARM(MODE,KIN,IERR)
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
      EXTERNAL WRNLIN,WRPLST
      CHARACTER KIN*(*)
C
    1 CALL TASK_PARM(MODE,'WR',KIN,WRNLIN,WRPLST,IERR)
      IF(IERR.NE.0) RETURN
C
      CALl WRCHEK(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100
C
      RETURN
      END
C
C     ****** INPUT NAMELIST ******
C
      SUBROUTINE WRNLIN(NID,IST,IERR)
C
      INCLUDE 'wrcomm.inc'
C
      NAMELIST /WR/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,
     &              r_corner,z_corner,
     &              br_corner,bz_corner,bt_corner,
     &              pn_corner,ptpr_corner,ptpp_corner,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,
     &              PPN0,PTN0,RFPL,
     &              MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,
     &              RHOGMN,RHOGMX,
     &              KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2,
     &              MODELP,NDISP1,NDISP2,
     &              MODELV,MDLWRI,MDLWRG,RHOGMN,RHOGMX,MDLWRQ,MDLWRW,
     &              RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,UUI,MODEW,
     &              SMAX,DELS,UUMIN,EPSNW,DELKR,EPSRAY,DELRAY,DELDER,
     &              LMAXNW,NRZMAX,NRAYMX,KNAMEQ,KNAMWR,
     &              RCURVA,RCURVB,RBRADA,RBRADB,NRADMX,
     &              RFIN,RPIN,ZPIN,PHIIN,RNZIN,RNPHIIN,RKRIN,UUIN,
     &              MODEWIN,
     &              ANGZIN,ANGPHIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN,
     &              MODEFW,MODEfR,IDEBUG

C
      READ(NID,WR,IOSTAT=IST,ERR=9800,END=9900)

      IF(MODEL_PROF.EQ.0) THEN
         DO NS=1,NSMAX
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
C
 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE WRPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &WR : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/
     &       9X,'NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,'/
     &       9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/
     &       9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/
     &       9X,'MODELG,MODELN,MODELQ,'/
     &       9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/
     &       9X,'MODELP,NDISP1,NDISP2,'/
     &       9X,'RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,'/
     &       9X,'RF1,RFI1,RKX1,RKY1,RKZ1,RX1,'/
     &       9X,'RF2,RFI2,RKX2,RKY2,RKZ2,RX2,'/
     &       9X,'NXMAX,EPSRT,LMAXRT,'/
     &       9X,'MODELV,MDLWRI,MDLWRG,RHOGMN,RHOGMX,MDLWRQ,MDLWRW,'/
     &       9X,'RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,UUI,MODEW,'/
     &       9X,'RFIN,RPIN,ZPIN,PHIIN,RNZIN,RNPHIIN,RKR0IN,UUIN,'/
     &       9X,'ANGZIN,ANGPHIN,MODEWIN,'/
     &       9X,'SMAX,DELS,UUMIN,EPSNW,DELKR,EPSRAY,DELRAY,DELDER,'/
     &       9X,'LMAXNW,NRZMAX,NRAYMX,KNAMWR,'/
     &       9X,'RCURVA,RCURVB,RBRADA,RBRADB,NRADMX,'/
     &       9X,'MODEFW,MODEFR,IDEBUG')
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE WRCHEK(IERR)
C
      INCLUDE 'wrcomm.inc'
C  
      CHARACTER(LEN=80):: LINE
      DATA INITEQ,INITFP/0,0/
C
      IERR=0
C
      IF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN
         IF(INITEQ.EQ.0) THEN
            CALL EQLOAD(MODELG,KNAMEQ,IERR)
            IF(IERR.EQ.0) THEN
               write(LINE,'(A,I5)') 'nrmax =',51
               call eqparm(2,line,ierr)
               write(line,'(a,i5)') 'nthmax=',64
               call eqparm(2,line,ierr)
               write(line,'(a,i5)') 'nsumax=',64
               call eqparm(2,line,ierr)
               call eqcalq(ierr)
               call eqgetb(bb,rr,rip,ra,rkap,rdlt,rb)
               initeq=1
            else
               write(6,*) 'xx eqload: ierr=',ierr
               initeq=0
            endif
         ENDIF
      ELSE IF(MODELG.EQ.8) THEN
         IF(INITEQ.EQ.0) THEN
            CALL EQREAD(IERR)
            IF(IERR == 0) THEN
               CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
               INITEQ=1
            ELSE
               WRITE(6,*) 'XX EQREAD: IERR=',IERR
               INITEQ=0
            ENDIF
         ENDIF
      ELSE
         INITEQ=0
      ENDIF
C
      DO NS=1,NSMAX
         IF(MODELV(NS).EQ.1) THEN
            IF(INITFP.EQ.0) THEN
               CALL DPLDFP
               INITFP=1
            ENDIF
         ELSE
            INITFP=0
         ENDIF
      ENDDO
C
      RETURN
      END

C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE WRVIEW
C
      INCLUDE 'wrcomm.inc'
C
      WRITE(6,601) 'SMAX  ',SMAX  ,'DELS  ',DELS  ,
     &             'UUMIN ',UUMIN
      WRITE(6,601) 'EPSNW ',EPSNW ,'DELKR ',DELKR
      WRITE(6,601) 'EPSRAY',EPSRAY,'DELRAY',DELRAY,
     &             'DELDER',DELDER
      WRITE(6,602) 'NRAYMX',NRAYMX,'LMAXNW',LMAXNW,
     &             'NRZMAX',NRZMAX,'NRADMX',NRADMX
      WRITE(6,602) 'MDLWRI',MDLWRI,'MDLWRG',MDLWRG,
     &             'MDLWRQ',MDLWRQ,'MDLWRW',MDLWRW
      RETURN
C
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END
