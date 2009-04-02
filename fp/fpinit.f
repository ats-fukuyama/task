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
      ZEFF  = 1.D0
      NSBEAM = 0
      NSFPMI = 1
      NSFPMA = 1
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
      RMIN  = 0.1D0
      RMAX  = 0.4D0
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
C     FACTWM: Numerical factor for wave amplitude
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
      FACTWM= 1.D0
C
C-----------------------------------------------------------------------
C     SPTOT(NSA) : Total number of particle source per second  [/s]
C     SPR0(NSA)  : Radial center of particle source            [m]
C     SPRW(NSA)  : Radial width of particle source             [m]
C     SPENG(NSA) : Energy of particle source                   [eV]
C     SPANG(NSA) : Pitch angle of particle source              [degree]
C
      DO NSA=1,NSAM
         SPTOT(NSA)= 0.D0
         SPR0(NSA) = 0.D0
         SPRW(NSA) = 0.D0
         SPENG(NSA)= 1.D6
         SPANG(NSA)= 0.D0
      ENDDO
C
C-----------------------------------------------------------------------
C     PMAX  : maximum momentum (normailzed by central thermal momentum)
C     DELT  : time step size (s)
C     RIMPL : implicit computation parameter
C     EPSM  : convergence limit in matrix equation solver
C     EPSE  : convergence limit in electric field prediction
C     LMAXE : maximum loop count in electric field prediction
C     EPSFP : convergence limit in nonlinear collision operator iteration
C     LMAXFP: maximum loop count in nonlinear collision operator iteration
C     EPSDE : convergence limit in double-exponential integration
C     H0DE  : initial step size in double-exponential integration
C
CCC      PMAX  = 15.D0
      PMAX  = 7.D0
      DELT  = 1.D-2
      RIMPL = 1.D0
      EPSM  = 1.D-8
      EPSE  = 1.D-4
      LMAXE = 10
      EPSFP = 1.D-7
      LMAXFP= 10
      EPSDE = 1.D-8
      H0DE  = 0.25D0
C
C-----------------------------------------------------------------------
C     PGMAX:  maximum p in graphics (if not 0)
C     RGMIN:  minimum rho in graphics (if not 0)
C     RGMAX:  maximum rho in graphics (if not 1)
C     NGLINE: maximum number of contour lines
C     NGRAPH: graphic mode: 0 for file out, 1 for contour plot

      PGMAX=10.D0
      RGMIN=0.D0
      RGMAX=1.D0
      NGLINE= 25
      NGRAPH=1
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
C     MODELC: 0 or 1: linear collision operator and (const./variable) T
C             2 or 3 : NL for same species and L for different species
C             4 : nonlinear collision operator
C             5 : linear coll. operator for different species (for debug)
C             6 : nonlinear coll. operator for different species (for debug)
C            -1 : linear collision operator for same with ion scattering
C            -2 : nonlinear collision operator for same with ion scattering
C     MODELW(NSA) :
C             0 : given diffusion coefficient model
C             1 : wave electric field calculated by WR(ray tracing)
C             2 : wave electric field calculated by WR(beam tracing)
C             3 : given wave electric field model
C             4 : wave electric field calculated by WM
C     MODELS(NSA) :
C             0 : no particle source
C             1 : beam particle source
C             2 : 3.5 MeV alpha particle source
C
      MODELE= 0
      MODELA= 1
      MODELR= 0
      MODELC= 0
      DO NSA=1,NSAM
         MODELW(NSA)= 0
         MODELS(NSA)= 0
      ENDDO
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
C-----------------------------------------------------------------------
C     IDBGFP : debug graphic control parameter 
C          1 : Legendre polynomials
C          2 : fpl, M_l, N_l
C          4 : psy, phy and their derivatives
C          8 : dcpp, dctt, fcp
C
      IDBGFP=0
C
      RETURN
      END
C
C ***********************
C     PARAMETER INPUT
C ***********************
C
      SUBROUTINE FPPARM(MODE,KIN,IERR)
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
      CHARACTER KIN*(*)
      EXTERNAL FPNLIN,FPPLST
C
    1 CALL TASK_PARM(MODE,'FP',KIN,FPNLIN,FPPLST,IERR)
      IF(IERR.NE.0) RETURN
C
      CALL FPCHEK(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100
C
      RETURN
      END
C
C     ****** INPUT NAMELIST ******
C
      SUBROUTINE FPNLIN(NID,IST,IERR)
C
      INCLUDE 'fpcomm.inc'
C
      NAMELIST /FP/ BB,RR,RA,RB,RKAP,RDLT,QA,Q0,RIP,PROFJ,
     &              NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,
     &              MODELG,MODELN,MODELQ,
     &              ZEFF,DRR0,E0,R1,DELR1,RMIN,RMAX,
     &              DEC,PEC1,PEC2,RFEC,DELYEC,
     &              DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW,
     &              SPTOT,SPR0,SPRW,SPENG,SPANG,
     &              PMAX,RIMPL,EPSM,EPSE,EPSFP,EPSDE,H0DE,LMAXE,LMAXFP,
     &              NPMAX,NTHMAX,NRMAX,NAVMAX,
     &              PGMAX,RGMIN,RGMAX,NGLINE,NGRAPH,LLMAX,
     &              DELT,NTMAX,NTSTP1,NTSTP2,NTSTPC,
     &              MODELE,MODELA,MODELC,MODELW,MODELS,MODELR,
     &              RFDW,DELNPR,NCMIN,NCMAX,FACTWM,
     &              CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH,REWY,DREWY,
     &              EPSNWR,LMAXNWR,PWAVE,DELCRI,NTHWAV,IDBGFP,
     &              KNAMEQ,KNAMWR,KNAMWM,KNAMFP,MODEFR,MODEFW,
     &              NSFPMI,NSFPMA
C
      READ(NID,FP,IOSTAT=IST,ERR=9800,END=9900)
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
      SUBROUTINE FPPLST
C
      WRITE(6,*) '# &FP : BB,RR,RA,RB,RKAP,RDLT,QA,Q0,RIP,PROFJ,'
      WRITE(6,*) '      NSMAX,PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,'
      WRITE(6,*) '      PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,*) '      RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,'
      WRITE(6,*) '      MODELG,MODELN,MODELQ,'
      WRITE(6,*) '      ZEFF,DRR0,E0,R1,DELR1,RMIN,RMAX,'
      WRITE(6,*) '      DEC,PEC1,PEC2,RFEC,DELYEC,'
      WRITE(6,*) '      DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW,'
      WRITE(6,*) '      SPTOT,SPR0,SPRW,SPENG,SPANG,'
      WRITE(6,*) '      EPSM,EPSE,EPSFP,EPSDE,H0DE,LMAXE,LMAXFP'
      WRITE(6,*) '      NPMAX,NTHMAX,NRMAX,NAVMAX,'
      WRITE(6,*) '      PGMAX,RGMIN.RGMAX,NGLINE,NGRAPH,LLMAX'
      WRITE(6,*) '      MODELE,MODELA,MODELC,MODELW,MODELS,MODELR,'
      WRITE(6,*) '      RFDW,DELNPR,PMAX,RIMPL,NCMIN,NCMAX,FACTWM,'
      WRITE(6,*) '      CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH,REWY,DREWY,'
      WRITE(6,*) '      EPSNWR,LMAXNWR,PWAVE,DELCRI,NTHWAV,IDBGFP,'
      WRITE(6,*) '      KNAMEQ,KNAMWR,KNAMWM,KNAMFP,MODEFR,MODEFW'
      RETURN
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE FPCHEK(IERR)
C
      INCLUDE 'fpcomm.inc'
C
      IERR=0
C
      IF(NRMAX+1.GT.NRM) THEN
         WRITE(6,*) 'FP : ERROR : NRMAX+1.GE.NRM'
         IERR=1
      ENDIF
      IF(NPMAX+1.GT.NPM) THEN
         WRITE(6,*) 'FP : ERROR : NPMAX+1.GE.NPM'
         IERR=2
      ENDIF
      IF(NTHMAX+1.GT.NTHM) THEN
         WRITE(6,*) 'FP : ERROR : NTHMAX+1.GE.NTHM'
         IERR=3
      ENDIF
      IF(MOD(NTHMAX,2).NE.0) THEN
         WRITE(6,*) 'FP : ERROR : NTHMAX MUST BE EVEN'
         IERR=4
      ENDIF
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
      DO NSA=NSFPMI,NSFPMA
         WRITE(6,'(A,I2)') 'NSA = ',NSA

      IF(MODELW(NSA).EQ.0) THEN
         WRITE(6,600) 'DEC   ',DEC   ,'RFEC  ',RFEC  ,
     &                'PEC1  ',PEC1  ,'PEC2  ',PEC2
         WRITE(6,600) 'DELYEC',DELYEC
         WRITE(6,600) 'DLH   ',DLH   ,'RLH   ',RLH   ,
     &                'PLH1  ',PLH1  ,'PLH2  ',PLH2
         WRITE(6,600) 'DFW   ',DFW   ,'RFW   ',RFW   ,
     &                'PFW1  ',PFW1  ,'PFW2  ',PFW2
C
      ELSEIF(MODELW(NSA).EQ.1) THEN
         WRITE(6,602) 'RFDW  ',RFDW  ,'DELNPR',DELNPR,
     &                'NCMIN ',NCMIN ,'NCMAX ',NCMAX
         WRITE(6,601) 'PWAVE ',PWAVE ,'DELYEC',DELYEC,
     &                'EPSNWR',EPSNWR,'LMAXNW',LMAXNWR
C
      ELSEIF(MODELW(NSA).EQ.2) THEN
         WRITE(6,602) 'RFDW  ',RFDW  ,'DELNPR',DELNPR,
     &                'NCMIN ',NCMIN ,'NCMAX ',NCMAX
         WRITE(6,601) 'PWAVE ',PWAVE ,'DELYEC',DELYEC,
     &                'EPSNWR',EPSNWR,'LMAXNW',LMAXNWR
C
      ELSEIF(MODELW(NSA).EQ.3) THEN
         WRITE(6,602) 'RFDW  ',RFDW  ,'DELNPR',DELNPR,
     &                'NCMIN ',NCMIN ,'NCMAX ',NCMAX
         WRITE(6,600) 'CEWR/R',DBLE(CEWR) ,'CEWR/I',DIMAG(CEWR),
     &                'CEWTHR',DBLE(CEWTH),'CEWTHI',DIMAG(CEWTH)
         WRITE(6,600) 'CEWPHR',DBLE(CEWPH),'CEWPHI',DIMAG(CEWPH),
     &                'REWY  ',REWY  ,'DREWY ',DREWY
         WRITE(6,600) 'RKWR  ',RKWR  ,'RKWTH ',RKWTH,
     &                'RKWPH ',RKWPH
C
      ELSEIF(MODELW(NSA).EQ.4) THEN
         WRITE(6,602) 'RFDW  ',RFDW  ,'DELNPR',DELNPR,
     &                'NCMIN ',NCMIN ,'NCMAX ',NCMAX
         WRITE(6,600) 'FACTWM',FACTWM
      ENDIF
C
      IF(MODELS(NSA).EQ.1) THEN
         WRITE(6,600) 'SPTOT ',SPTOT(NSA) ,'SPR0  ',SPR0(NSA),
     &                'SPRW  ',SPRW(NSA)
         WRITE(6,600) 'SPENG ',SPENG(NSA) ,'SPANG ',SPANG(NSA)
      ELSEIF(MODELS(NSA).EQ.2) THEN
         WRITE(6,600) 'SPTOT ',SPTOT(NSA) ,'SPR0  ',SPR0(NSA),
     &                'SPRW  ',SPRW(NSA)
      ENDIF
      ENDDO
C
      WRITE(6,600) 'PMAX  ',PMAX  ,'DELT  ',DELT  ,
     &             'RIMPL ',RIMPL ,'EPSM  ',EPSM
C
      WRITE(6,600) 'EPSDE ',EPSDE ,'H0DE  ',H0DE  ,
     &             'EPSE  ',EPSE  ,'EPSFP ',EPSFP

      WRITE(6,604) 'LMAXE ',LMAXE ,'LMAXFP',LMAXFP

      WRITE(6,601) 'PGMAX ',PGMAX ,'RGMIN ',RGMIN ,
     &             'RGMAX ',RGMAX

      WRITE(6,604) 'LLMAX ',LLMAX ,'NGLINE',NGLINE,
     &             'IDBGFP',IDBGFP,'NGRAPH',NGRAPH
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
      WRITE(6,604) 'MODEFR',MODEFR,'MODEFW',MODEFW
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
         WRITE(6,*) 'LINEAR COLLISION OPERATOR & CONST. T'
      ELSE IF(MODELC.eq.1)THEN
         WRITE(6,*) 'LINEAR COLLISION OPERATOR & VARIABLE. T'
      ELSE IF(MODELC.EQ.2)THEN
         WRITE(6,*) 
     &    'NONLINEAR COLLISION OPERATOR FOR LIKE PARTILCES & CONST. T'
      ELSE IF(MODELC.eq.3)THEN
         WRITE(6,*) 
     &  'NONLINEAR COLLISION OPERATOR FOR LIKE PARTILCES & VARIABLE T'
      ELSE IF(MODELC.eq.4)THEN
         WRITE(6,*) 'NONLINEAR COLLISION OPERATOR'
      ELSE IF(MODELC.EQ.5)THEN
         WRITE(6,*) 'NONLINEAR COLLISION OPERATOR'
      ELSE IF(MODELC.EQ.-1)THEN
         WRITE(6,*) 'LINEAR COLLISION OPERATOR WITH ION SCATTERING'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELC: MODELC =',MODELC
      END IF
C
      DO NSA=NSFPMI,NSFPMA
         WRITE(6,'(A,I2)') 'NSA = ',NSA

      IF(MODELW(NSA).EQ.0)THEN
         WRITE(6,*) 'GIVEN WAVE DIFFUSION COEFFICIENTS'
      ELSE IF(MODELW(NSA).EQ.1)THEN
         WRITE(6,*) 'RAY TRACING WAVE DATA'
      ELSE IF(MODELW(NSA).EQ.2)THEN
         WRITE(6,*) 'BEAM TRACING WAVE DATA'
      ELSE IF(MODELW(NSA).EQ.3)THEN
         WRITE(6,*) 'GIVEN WAVE AMPLITUDE'
      ELSE IF(MODELW(NSA).EQ.4)THEN
         WRITE(6,*) 'FULL WAVE DATA'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELW: MODELW =',MODELW
      END IF

      IF(MODELS(NSA).EQ.1) THEN
         WRITE(6,*) 'BEAM PARTICLE SOURCE'
      ELSEIF(MODELS(NSA).EQ.2) THEN
         WRITE(6,*) 'ALPHA PARTICLE SOURCE'
      ENDIF
      ENDDO

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
