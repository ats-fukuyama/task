C     $Id$
C
C ***************************
C     INITIAL PARAMETERS
C ***************************
C
      SUBROUTINE FPINIT
C
      INCLUDE 'fpcomm.inc'
C
C-----------------------------------------------------------------------
C     ZEFF  : effective ion charge
C     NSFP  : particle species of distribution
C
      ZEFF  = 2.D0
      NSFP  = 1
C
C-----------------------------------------------------------------------
C     DRR0  : radial diffusion coefficient (m^2/s)
C     E0    : toroidal electric field (V/m)
C     R1    : radial position for NRMAX=1 (r/a)
C     DELR1 : radial spacing for NRMAX=1 (r/a)
C     RMIN  : minimum radius for NRMAX<>1 (r/a)
C     RMIN  : maximum radius for NRMAX<>1 (r/a)
C
      DRR0  = 0.D0
      E0    = 0.D0
      R1    = 0.1D0*RR
      DELR1 = 0.1D0
      RMIN  = 0.05D0
      RMAX  = 0.3D0
C
C-----------------------------------------------------------------------
C     DEC   : Electron Cyclotron wave diffusion coefficient (normalized)
C     PEC1  : EC wave spectrum N para center
C     PEC2  : EC wave spectrum N para width
C     RFEC  : EC wave frequency / electron cyclotron grequency at r = 0
C     DELYEC: EC wave beam vertical  width (m)
C     DLH   : Lower Hybrid wave diffusion coefficient (normalized)
C     PLH1  : LH wave spectrum (velocity minimum/N para center)
C     PLH2  : LH wave spectrum (velocity maximum/N para width)
C     RLH   : LH wave minimum accessible minor radius (m)
C     DFW   : Fast Wave diffusion coefficient (normalized)
C     PFW1  : FW spectrum (velocity minimum/N para center)
C     PFW2  : FW spectrum (velocity maximum/N para width)
C     RFW   : FW mimum accessible minor radius (m)
C
      DEC   = 0.D0
      PEC1  = 0.5D0
      PEC2  = 0.1D0
      RFEC  = 0.8D0
      DELYEC= 0.1D0
      DLH   = 0.D0
C      PLH1  = 4.D0
C      PLH2  = 7.D0
      PLH1  = 1.92D0
      PLH2  = 0.3D0
      RLH   = 0.D0
      DFW   = 0.D0
C      PFW1  = 4.D0
C      PFW2  = 7.D0
      PFW1  = 0.D0
      PFW2  = 1.1D0
      RFW   = 0.D0
C
C-----------------------------------------------------------------------
C     RFDW  : wave frequency [MHz]
C     DELNPR: width of toroidal mode number
C     NCMIN : minimum order of cyclotron harmonics
C     NCMAX : maximum order of cyclotron harmonics
C     CEWR  : radial component of wave electric field [V/m]
C     CEWTH : poloidal component of wave electric field [V/m]
C     CEWPH : toroidal component of wave electric field [V/m]
C     RKWR  : radial component of wave number [1/m]
C     RKWTH : poloidal component of wave number [1/m]
C     RKWPH : toroidal component of wave number [1/m]
C     REWY  : vertical position of ray [r/a]
C     DREWY : vertical half-width of ray [r/a]
C
      RFDW  = 170.D3
      DELNPR= 0.05D0
      NCMIN = -3
      NCMAX =  3
      CEWR  = (0.D0,0.D0)
      CEWTH = (0.D0,0.D0)
      CEWPH = (1.D3,0.D0)
      RKWR  = 0.D0
      RKWTH = 0.D0
      RKWPH = 1.D3
      REWY  = 0.D0
      DREWY = 0.1D0
C
C-----------------------------------------------------------------------
C     PMAX  : maximum momentum (normailzed by central thermal momentum)
C     DELT  : time step size (s)
C     RIMPL : implicit computation parameter
C     EPSM  : convergence limit in matrix equation solver
C     EPSE  : convergence limit in electric field prediction
C     LMAXE : maximum loop count in electric field prediction
C     EPSDE : convergence limit in double-exponential integration
C     H0DE  : initial step size in double-exponential integration
C     NGLINE: maximum number of contour lines
C
CCC      PMAX  = 15.D0
      PMAX  = 7.D0
      DELT  = 1.D-2
      RIMPL = 1.D0
      EPSM  = 1.D-8
      EPSE  = 1.D-4
      LMAXE = 10
      EPSDE = 1.D-8
      H0DE  = 0.25D0
      NGLINE= 30
C
C-----------------------------------------------------------------------
C     NPMAX : momentum magnitude division number
C     NTHMAX: momentum angle division number
C     NRMAX : radial division number
C     NAVMAX: wave diuffusion bounce average division number
C
      NPMAX = 50
      NTHMAX= 50
      NRMAX = 1
      NAVMAX= 100
C
C-----------------------------------------------------------------------
C     NTMAX : maximum time step count
C     NTSTP1: time step interval for storing radial profile  data
C     NTSTP2: time step interval for storing global data
C     NTSTPC: time step interval for recalculating coefficients
C
      NTMAX = 10
      NTSTP1= 1
      NTSTP2= 1
      NTSTPC= 1000
C
C-----------------------------------------------------------------------
C     MODELE: 0 for fixed electric field
C             1 for predicting electric field
C     MODELA: 0 without bounce average
C             1 with bounce average
C     MODELR: 0 without relativistic effect
C             1 with relativistic effect
C     MODELC: 0 for linear collision operator
C             1 for nonlinear collision operator
C     MODELW: 0 for given diffusion coefficient model
C             1 for given wave electric field model
C             2 for wave electric field calculated by WR(without beam radius)
C             3 for wave electric field calcurated by WR(with  beam radius)
C
      MODELE= 0
      MODELA= 1
      MODELR= 1
      MODELC= 1
      MODELW= 0
C
C-----------------------------------------------------------------------
C     LLMAX : dimension of legendre polynomials's calculation
C
      LLMAX = 2
C
C-----------------------------------------------------------------------
C     PWAVE : input power
C     LMAXNWR: max loop count in newton method to find ray position
C     EPSNWR: convergence criterion in newton method to find ray position
C
      PWAVE = 1.0D0
      LMAXNWR=100
      EPSNWR=1.D-6
C      
      RETURN
      END
C
C ***********************
C     PARAMETER INPUT
C ***********************
C
      SUBROUTINE FPPARM(KID)
C
      INCLUDE 'fpcomm.inc'
C
      NAMELIST /FP/ BB,RR,RA,RB,RKAP,RDLT,QA,Q0,RIP,PROFJ,
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,
     &              MODELG,MODELN,MODELQ,RHOGMN,RHOGMX,
     &              ZEFF,NSFP,DRR0,E0,R1,DELR1,RMIN,RMAX,
     &              DEC,PEC1,PEC2,RFEC,DELYEC,
     &              DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW,
     &              PMAX,RIMPL,EPSM,EPSE,EPSDE,H0DE,LMAXE,
     &              NPMAX,NTHMAX,NRMAX,NAVMAX,NGLINE,
     &              DELT,NTMAX,NTSTP1,NTSTP2,NTSTPC,
     &              MODELE,MODELA,MODELC,MODELW,MODELR,LLMAX,
     &              RFDW,DELNPR,NCMIN,NCMAX,
     &              CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH,REWY,DREWY,
     &              EPSNWR,LMAXNWR,PWAVE,DELCRI,NTHWAV
      LOGICAL LEX
      CHARACTER KPNAME*72,KLINE*70,KNAME*80,KID*1
C
      MODE=0
    1    CONTINUE
         WRITE(6,*) '## INPUT &FP : '
         READ(5,FP,ERR=2,END=3)
         KID=' '
         GOTO 4
C
    2    CALL FPPLST
      GOTO 1
C
    3 KID='Q'
    4 GOTO 3000
C
      ENTRY FPPARL(KLINE)
C
      MODE=1
      KNAME=' &FP '//KLINE//' &END'
      WRITE(33,'(A80)') KNAME
      REWIND(33)
      READ(33,FP,ERR=8,END=8)
      WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
      GOTO 9
    8 CALL FPPLST
    9 REWIND(33)
      GOTO 3000
C
      ENTRY FPPARF
C
      MODE=2
      KPNAME='fpparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(.NOT.LEX) RETURN
C
      OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
      READ(25,FP,IOSTAT=IST,ERR=9800,END=9900)
      WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
C
 3000 IERROR=0
      IF(NRMAX.GT.NRM) THEN
         WRITE(6,*) 'FP : ERROR : NRMAX.GE.NRM'
         IERROR=1
      ENDIF
      IF(NPMAX.GT.NPM) THEN
         WRITE(6,*) 'FP : ERROR : NPMAX.GE.NPM'
         IERROR=1
      ENDIF
      IF(NTHMAX.GT.NTHM) THEN
         WRITE(6,*) 'FP : ERROR : NTHMAX.GE.NTHM'
         IERROR=1
      ENDIF
      IF(MOD(NTHMAX,2).NE.0) THEN
         WRITE(6,*) 'FP : ERROR : NTHMAX MUST BE EVEN'
         IERROR=1
      ENDIF
C
      IF(IERROR.NE.0.AND.MODE.EQ.0) GOTO 1
C
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR : IOSTAT = ',IST
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE FPPLST
C
      WRITE(6,*) '# &FP : BB,RR,RA,RB,RKAP,RDLT,QA,Q0,RIP,PROFJ,'
      WRITE(6,*) '      NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,'
      WRITE(6,*) '      PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,*) '      RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,'
      WRITE(6,*) '      MODELG,MODELN,MODELQ,RHOGMN,RHOGMX,'
      WRITE(6,*) '      ZEFF,NSFP,DRR0,E0,R1,DELR1,RMIN,RMAX,'
      WRITE(6,*) '      DEC,PEC1,PEC2,RFEC,DELYEC,'
      WRITE(6,*) '      DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW,'
      WRITE(6,*) '      PMAX,RIMPL,EPSM,EPSE,EPSDE,H0DE,LMAXE,'
      WRITE(6,*) '      NPMAX,NTHMAX,NRMAX,NAVMAX,NGLINE,'
      WRITE(6,*) '      MODELE,MODELA,MODELC,MODELW,MODELR,LLMAX'
      WRITE(6,*) '      RFDW,DELNPR,NCMIN,NCMAX,'
      WRITE(6,*) '      CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH,REWY,DREWY,'
      WRITE(6,*) '      EPSNWR,LMAXNWR,PWAVE,DELCRI,NTHWAV'
      RETURN
      END
C
C ***********************
C     PARAMETER VIEW
C ***********************
C
      SUBROUTINE FPVIEW
C
      INCLUDE 'fpcomm.inc'
C
      WRITE(6,600) 'E0    ',E0    ,'DRR0  ',DRR0  ,
     &             'ZEFF  ',ZEFF
C
      WRITE(6,600) 'R1    ',R1    ,'DELR1 ',DELR1 ,
     &             'RMIN  ',RMIN  ,'RMAX  ',RMAX
C
      IF(MODELW.EQ.0) THEN
C
         WRITE(6,600) 'DEC   ',DEC   ,'RFEC  ',RFEC  ,
     &                'PEC1  ',PEC1  ,'PEC2  ',PEC2
         WRITE(6,600) 'DELYEC',DELYEC
C
         WRITE(6,600) 'DLH   ',DLH   ,'RLH   ',RLH   ,
     &                'PLH1  ',PLH1  ,'PLH2  ',PLH2
C
         WRITE(6,600) 'DFW   ',DFW   ,'RFW   ',RFW   ,
     &                'PFW1  ',PFW1  ,'PFW2  ',PFW2
C
      ELSEIF(MODELW.EQ.1) THEN
         WRITE(6,602) 'RFDW  ',RFDW  ,'DELNPR',DELNPR,
     &                'NCMIN ',NCMIN ,'NCMAX ',NCMAX
C
         WRITE(6,600) 'CEWR/R',DBLE(CEWR) ,'CEWR/I',DIMAG(CEWR),
     &                'CEWTHR',DBLE(CEWTH),'CEWTHI',DIMAG(CEWTH)
C
         WRITE(6,600) 'CEWPHR',DBLE(CEWPH),'CEWPHI',DIMAG(CEWPH),
     &                'REWY  ',REWY  ,'DREWY ',DREWY
C
         WRITE(6,600) 'RKWR  ',RKWR  ,'RKWTH ',RKWTH,
     &                'RKWPH ',RKWPH
C
      ELSEIF(MODELW.EQ.2) THEN
         WRITE(6,602) 'RFDW  ',RFDW  ,'DELNPR',DELNPR,
     &                'NCMIN ',NCMIN ,'NCMAX ',NCMAX
C
         WRITE(6,601) 'PWAVE ',PWAVE ,'DELYEC',DELYEC,
     &                'EPSNWR',EPSNWR,'LMAXNW',LMAXNWR
C
      ENDIF
C
      WRITE(6,600) 'PMAX  ',PMAX  ,'DELT  ',DELT  ,
     &             'RIMPL ',RIMPL ,'EPSM  ',EPSM
C
      WRITE(6,601) 'EPSDE ',EPSDE ,'H0DE  ',H0DE  ,
     &             'EPSE  ',EPSE  ,'LMAXE ',LMAXE

      WRITE(6,604) 'LLMAX ',LLMAX ,'NGLINE',NGLINE,
     &             'NSFP  ',NSFP
C
      WRITE(6,604) 'NPMAX ',NPMAX ,'NTHMAX',NTHMAX,
     &             'NRMAX ',NRMAX ,'NAVMAX',NAVMAX
C
      WRITE(6,604) 'NTMAX ',NTMAX ,'NTSTP1',NTSTP1,
     &             'NTSTP2',NTSTP2,'NTSTPC',NTSTPC
C
      WRITE(6,604) 'MODELE',MODELE,'MODELA',MODELA ,
     &             'MODELC',MODELC,'MODELW',MODELW
C
      IF(MODELE.EQ.0)THEN
         WRITE(6,*) 'FIXED ELECTRIC FIELD'
      ELSE IF(MODELE.EQ.1)THEN
         WRITE(6,*) 'CONSISTENT ELECTRIC FIELD'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELE: MODELE =',MODELE
      ENDIF
C
      IF(MODELR.EQ.0)THEN
         WRITE(6,*) 'NONRELATIVISTIC'
      ELSE IF(MODELR.EQ.1)THEN
         WRITE(6,*) 'RELATIVISTIC'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELR: MODELR =',MODELR
      ENDIF
C
      IF(MODELA.EQ.0)THEN
         WRITE(6,*) 'NOT BOUNCE AVERAGED'
      ELSE IF(MODELA.EQ.1)THEN
         WRITE(6,*) 'BOUNCE AVERAGED'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELA: MODELA =',MODELA
      ENDIF
C
      IF(MODELC.EQ.0)THEN
         WRITE(6,*) 'MAXWELLIAN COLLISION OPERATOR'
      ELSE IF(MODELC.GE.1)THEN
         WRITE(6,*) 'NONLINEAR COLLISION OPERATOR'
      ELSE IF(MODELC.EQ.11)THEN
         WRITE(6,*) 'LINEAR COLLISION OPERATOR WITH ION SCATTERING'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELC: MODELC =',MODELC
      END IF
C
      IF(MODELW.EQ.0)THEN
         WRITE(6,*) 'GIVEN WAVE DIFFUSION COEFFICIENTS'
      ELSE IF(MODELW.EQ.1)THEN
         WRITE(6,*) 'GIVEN WAVE AMPLITUDE'
      ELSE IF(MODELW.EQ.2)THEN
         WRITE(6,*) 'RAY TRACING WAVE DATA'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELW: MODELW =',MODELW
      END IF
C
      IF(MODELG.EQ.2)THEN
         WRITE(6,*) 'GIVEN PLASMA GEOMETRY'
      ELSE IF(MODELG.EQ.3)THEN
         WRITE(6,*) 'MHD EQUILIBRIUM FROM TASK/EQ'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELG: MODELG =',MODELG
      END IF
C
      RETURN
C
  600 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',I7)
  602 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  604 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
      END
