C     $Id$
C
C     ****** INITIALIZE CONSTANTS & PARAMETERS ******
C
      SUBROUTINE WMINIT
C
      INCLUDE 'wmcomm.h'
C
C     *** CONSTANTS ****
C
C        CI    : Imaginary unit
C        PI    : Pi
C        AEE   : Elementaty charge
C        AME   : Electron mass
C        AMP   : Proton mass
C        VC    : Speed of light in vacuum
C        RMU0  : Permeability of free space
C        EPS0  : Permittivity of free space
C
      CI     = (0.D0,1.D0)
      PI     = ACOS(0.D0)*2.D0
      AEE    = 1.60217733D-19
      AME    = 9.1093897D-31
      AMP    = 1.6726231D-27
      VC     = 2.997924580D8
      RMU0   = 4.D0*PI*1.D-7
      EPS0   = 1.D0/(VC*VC*RMU0)
C
C     *** CONFIGURATION PARAMETERS ***
C
C        BB    : Magnetic field at center                        (T)
C        RR    : Plasma major radius                             (m)
C        RA    : Plasma minor radius                             (m)
C        RB    : Wall minor radius                               (m)
C        Q0    : Safety factor at center
C        QA    : Safety factor on plasma surface
C        RKAP  : Plasma shape elongation
C        RDEL  : Plasma shape triangularity *
C
      BB     = 3.00D0
      RR     = 3.00D0
      RA     = 1.00D0
      RB     = 1.20D0
      Q0     = 1.00D0
      QA     = 3.00D0
      RKAP   = 1.00D0
      RDEL   = 0.00D0
C
C     *** PLASMA PARAMETERS ***
C
C        NSMAX : Number of particle species
C        PA    : Mass number
C        PZ    : Charge number
C        PN    : Density at center                      (1.0E20/m^3)
C        PNS   : Density on plasma surface              (1.0E20/m^3)
C        PTPR  : Parallel temperature at center                (keV)
C        PTPP  : Perpendicular temperature at center           (keV)
C        PTS   : Temperature on surface                        (keV)
C        PZCL  : Ratio of collision frequency to wave frequency
C        PU    : Toroidal rotation velocity at center          (m/s)
C        PUS   : Toroidal rotation velocity on surface         (m/s)
C        PNITB : Density increment at ITB               (1.0E20/m^3)
C        PTITB : Temperature increment at ITB                  (keV)
C        PUITB : Toroidal rotation velocity increment at ITB   (m/s)
C
C        ZEFF  : Effective Z (\sum n Z^2 / \sum n Z)
C        PROFN1: Density profile parameter (power of rho)
C        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
C        PROFT1: Temperature profile parameter (power of rho)
C        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
C        PROFU1: Rotation profile parameter (power of rho)
C        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))
C
      NSMAX  = 3
C
C     < E >
C
      PA(1)   =  AME/AMP
      PZ(1)   = -1.0D0
      PN(1)   =  1.D0
      PNS(1)  =  0.1D0
      PZCL(1) =  3.D-2
      PTPR(1) =  3.D0
      PTPP(1) =  3.D0
      PTS(1)  =  0.3D0
      PU(1)   =  0.D0
      PUS(1)  =  0.D0
      PNITB(1)=  0.D0
      PTITB(1)=  0.D0
      PUITB(1)=  0.D0
C
C     < D >
C
      PA(2)   =  2.0D0
      PZ(2)   =  1.0D0
      PN(2)   =  0.9D0
      PNS(2)  =  0.09D0
      PZCL(2) =  3.0D-2
      PTPR(2) =  3.D0
      PTPP(2) =  3.D0
      PTS(2)  =  0.3D0
      PU(2)   =  0.D0
      PUS(2)  =  0.D0
      PNITB(2)=  0.D0
      PTITB(2)=  0.D0
      PUITB(2)=  0.D0
C
C     < H >
C
      IF(NSM.GE.3) THEN
         PA(3)   =  1.0D0
         PZ(3)   =  1.0D0
         PN(3)   =  0.1D0
         PNS(3)  =  0.01D0
         PZCL(3) =  3.0D-2
         PTPR(3) =  3.D0
         PTPP(3) =  3.D0
         PTS(3)  =  0.3D0
         PU(3)   =  0.D0
         PUS(3)  =  0.D0
         PNITB(3)=  0.D0
         PTITB(3)=  0.D0
         PUITB(3)=  0.D0
      ENDIF
C
C     < He >
C
      IF(NSM.GE.4) THEN
         NS=4
         PA(NS)   =  4.0D0
         PZ(NS)   =  2.0D0
         PN(NS)   =  0.5D0
         PNS(NS)  =  0.05D0
         PZCL(NS) =  3.0D-2
         PTPR(NS) =  3.D0
         PTPP(NS) =  3.D0
         PTS(NS)  =  1.D0
         PU(NS)   =  0.D0
         PUS(NS)  =  0.D0
         PNITB(NS)=  0.D0
         PTITB(NS)=  0.D0
         PUITB(NS)=  0.D0
      ENDIF
C
      ZEFF  = 2.D0
      PROFN1= 2.D0
      PROFN2= 0.5D0
      PROFT1= 2.D0
      PROFT2= 1.D0
      PROFU1= 2.D0
      PROFU2= 1.D0
C
C     **** ALPHA PARTICLE PARAMETERS ****
C
C        PNA   : Alpha density at center               (1.0E20/Mm*3)
C        PNAL  : Density scale length                            (m)
C        PTA   : Effective temperature                         (keV)
C
      PNA  = 0.02D0
      PNAL = 0.5D0
      PTA  = 3.5D3
C
C     *** WAVE PARAMETER ***
C
C     CRF   : Wave frequency                            (MHz)
C     RD    : Antenna minor radius                      (m)
C     BETAJ : Antenna current profile parameter
C     NTH0  : Central value of poloidal mode number
C     NPH0  : Central value of toroidal mode number
C     NHC   : Number of helical coils
C     PRFIN : Input Power (0 for given antenna current) (W)
C
      CRF    = (50.0D0,0.D0)
      RD     = 1.1D0
      BETAJ  = 0.D0
      NTH0   = 0
      NPH0   = 8
      NHC    = 10
      PRFIN  = 0.D0
C
C     *** ANTENNA PARAMETERS ***
C
C        NAMAX : Number of antennae
C        AJ    : Antenna current density                       (A/m)
C        APH   : Antenna phase                              (degree)
C        THJ1  : Start poloidal angle of antenna            (degree)
C        THJ2  : End poloidal angle of antenna              (degree)
C        PHJ1  : Start toroidal angle of antenna            (degree)
C        PHJ2  : End toroidal angle of antenna              (degree)
C
      NAMAX=1
      DO NA=1,NAM
         AJ(NA)   = 1.D0
         APH(NA)  = 0.D0
         THJ1(NA) =-45.D0
         THJ2(NA) = 45.D0
         PHJ1(NA) = 0.D0
         PHJ2(NA) = 0.D0
      ENDDO
C
C     *** MESH PARAMETERS ***
C
C        NRMAX : Number of radial mesh points
C        NTHMAX: Number of poloidal mesh points
C        NPHMAX : Number of toroidal mesh points
C
      NRMAX   = 50
      NTHMAX  = 1
      NPHMAX  = 1
C
C     *** CONTROL PARAMETERS ***
C
C        NPRINT: Control print output
C                   0: No print out
C                   1: Minimum print out (without input data)
C                   2: Minimum print out (with input data)
C                   3: Standard print out
C                   4: More print out
C        NGRAPH: Control graphic output
C                   0: No graphic out
C                   1: Standard graphic out (2D: Coutour)
C                   2: Standard graphic out (2D: Paint)
C                   3: Standard graphic out (2D: Bird's eye)
C        MODELG: Control plasma geometry model
C                   0: Cylindrical geometry
C                   1: Toroidal geometry
C                   2: Straight helical geometry
C                   3: TASK/EQ output geometry
C                   4: VMEC output geometry
C                   5: EQDSK output geometry
C        MODELJ: Control antenna current model
C                   0: Real antenna
C                   1: Real antenna
C                   2: Poloidal current
C                   3: Toroidal current
C                  2X: Vacuum eigen mode, poloidal current
C                  3X: Vacuum eigen mode, toroidal current
C        MODELP: Control plasma dielectric tensor model
C                   0: Vacuum
C                   1: MHD plasma
C                   2: Cold plasma
C                   3: Hot plasma (No FLR)
C                   4: Hot plasma (Cold FLR)
C                   5: Hot plasma (FLR)
C        MODELN: Control plasma profile
C                   0: Calculate from PN,PNS,PTPR,PTPP,PTS,PU,PUS
C                   8: Read from file by means of WMDPRF routine (DIII-D)
C                   9: Read from file by means of WMXPRF routine (JT-60)
C        MODELA: Control alpha particle contribution
C                   0: No alpha effect
C                   1: Precession of alpha particles
C                   2: Precession of electrons
C                   3: Precession of both alpha particles and electrons
C                   4: Calculate alpha particle density using slowing down
C        MODELK: Control mode number cutoff
C                   0: No cutoff
C                   1: With cutoff (this should not be used)
C        MODELM: Control matrix solver
C                   0: BANDCD
C                   1: BCGCDB
C                   2: CGSCDB
C                   3: BCGSTAB
C                   4: BANDCDM
C                   5: BCGCDBM
C                   6: CGSCDBM
C                   7: BSTABCDBM
C                   8: BANDCDBM
C                   9: BCGCDBMA
C                  10: CGSCDBMA
C                  11: BSTABCDBMA
C                  12: BANDCDB
C        MODELW: Control writting a data of absorped power
C                   0: Not writting
C                   1: Writting
C
C        RHOMIN: rho at minimum q for reversed shear
C        QMIN  : q minimum for reversed shear
C        RHOITB: rho at ITB
C
      NPRINT = 2
      NGRAPH = 1
      MODELG = 1
      MODELJ = 0
      MODELP = 4
      MODELN = 0
      MODELA = 0
      MODELK = 0
      MODELM = 2
      MODELW = 0
C
      RHOMIN = 0.D0
      QMIN   = 1.5D0
      RHOITB = 0.D0
C
C     *** FILE NAME ***
C
C        KNAMEQ: Filename of equilibrium data
C        KNAMPF: Filename of profile data
C
      KNAMEQ = 'eqdata'
      KNAMPF = 'pfdata'
C
C     *** EIGEN VALUE PARAMETERS ***
C
C        FRMIN : Minimum real part of frequency in amplitude scan
C        FRMAX : Maximum real part of frequency in amplitude scan
C        FIMIN : Minimum imag part of frequency in amplitude scan
C        FIMAX : Maximum imag part of frequency in amplitude scan
C        FI0   : Imag part of frequency in 1D amplitude scan
C
C        NGFMAX: Number of real freq mesh in 1D amplitude scan
C        NGXMAX: Number of real freq mesh in 2D amplitude scan
C        NGYMAX: Number of imag freq mesh in 2D amplitude scan
C
C        SCMIN : Minimum value in parameter scan
C        SCMAX : Maximum value in parameter scan
C        NSCMAX: Number of mesh in parameter scan
C
C        LISTEG: Listing in parameter scan
C
C        FRINI : Initial real part of frequency in Newton method
C        FIINI : Initial imag part of frequency in Newton method
C
C        DLTNW : Step size in evaluating derivatives in Newton method
C        EPSNW : Convergence criterion in Newton method
C        LMAXNW: Maximum iteration count in Newton method
C        LISTNW: Listing in Newton method
C        MODENW: Type of Newton method
C
C        NCONT : Number of contour lines
C        ILN1  : Line type of lower contours
C        IBL1  : Line boldness of lower contours
C        ICL1  : Line color of lower contours
C        ILN2  : Line type of higher contours
C        IBL2  : Line boldness of higher contours
C        ICL2  : Line color of higher contours
C
      FRMIN = 0.1D0
      FRMAX = 1.D0
      FIMIN =-0.1D0
      FIMAX = 0.1D0
      FI0   = 0.D0
C
      FRINI = DBLE(CRF)
      FIINI = DIMAG(CRF)
C
      NGFMAX= 11
      NGXMAX= 11
      NGYMAX= 11
C
      SCMIN = 0.1D0
      SCMAX = 1.D0
      NSCMAX= 11
C
      LISTEG= 1
C
      DLTNW = 1.D-6
      EPSNW = 1.D-6
      LMAXNW= 10
      LISTNW= 1
      MODENW= 0
C
C     *** ALFVEN FREQUENCY PARAMETERS ***
C
C        WAEMIN : Minimum frequency in Alfven frequency scan
C        WAEMAX : Maximum frequency in Alfven frequency scan
C
      WAEMIN = 0.001D0
      WAEMAX = 0.200D0
C
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE WMPLST
C
      WRITE(6,*) '## INPUT &WM : BB,RR,RA,RB,Q0,QA,RKAP,RDEL,'
      WRITE(6,*) '               PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,'
      WRITE(6,*) '               PROFN1,PROFN2,PROFT1,PROFT2,ZEFF,'
      WRITE(6,*) '               NSMAX,PNA,PNAL,PTA,RF,RFI,RD,BETAJ,'
      WRITE(6,*) '               AJ,APH,THJ1,THJ2,PHJ1,PHJ2,NAMAX,'
      WRITE(6,*) '               NRMAX,NTHMAX,NPHMAX,NTH0,NPH0,NHC,'
      WRITE(6,*) '               MODELG,MODELJ,MODELP,MODELA,MODELN,'
      WRITE(6,*) '               MODELM,MODELW,KNAMEQ,KNAMPF,'
      WRITE(6,*) '               NPRINT,NGRAPH,PRFIN,'
      WRITE(6,*) '               FRMIN,FRMAX,FIMIN,FIMAX,FI0,'
      WRITE(6,*) '               FRINI,FIINI,NGFMAX,NGXMAX,NGYMAX,'
      WRITE(6,*) '               SCMIN,SCMAX,NSCMAX,LISTEG,'
      WRITE(6,*) '               DLTNW,EPSNW,LMAXNW,LISTNW,MODENW,'
      WRITE(6,*) '               RHOMIN,QMIN,PU,PUS,PROFU1,PROFU2'
      WRITE(6,*) '               RHOITB,PNITB,PTITB,PUITB'
      WRITE(6,*) '               WAEMIN,WAEMAX,KNAMEQ,KNAMPF'
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE WMPARM(KID)
C
      INCLUDE 'wmcomm.h'
C
      LOGICAL LEX
      CHARACTER KPNAME*6,KLINE*70,KNAME*80,KID*1
C
      NAMELIST /WM/ BB,RR,RA,RB,Q0,QA,RKAP,RDEL,
     &              PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,NSMAX,
     &              PROFN1,PROFN2,PROFT1,PROFT2,ZEFF,
     &              PNA,PNAL,PTA,
     &              RF,RFI,RD,BETAJ,AJ,APH,THJ1,THJ2,PHJ1,PHJ2,NAMAX,
     &              NRMAX,NTHMAX,NPHMAX,NTH0,NPH0,NHC,
     &              NPRINT,NGRAPH,MODELG,MODELJ,MODELP,MODELN,MODELA,
     &              MODELM,MODELW,KNAMEQ,KNAMPF,
     &              FRMIN,FRMAX,FIMIN,FIMAX,FI0,FRINI,FIINI,
     &              NGFMAX,NGXMAX,NGYMAX,SCMIN,SCMAX,NSCMAX,LISTEG,
     &              DLTNW,EPSNW,LMAXNW,LISTNW,MODENW,
     &              RHOMIN,QMIN,PU,PUS,PROFU1,PROFU2,
     &              RHOITB,PNITB,PTITB,PUITB,WAEMIN,WAEMAX,
     &              KNAMEQ,KNAMPF,PRFIN
C
      MODE=0
 1000 CONTINUE
      RF=DBLE(CRF)
      RFI=DIMAG(CRF)
    1    CONTINUE
         WRITE(6,*) '# INPUT : &WM'
         READ(5,WM,ERR=2,END=3)
         KID=' '
         GOTO 4
C
    2    CALL WMPLST
      GOTO 1
C
    3 KID='Q'
    4 CONTINUE
      CRF=DCMPLX(RF,RFI)
      GOTO 3000
C
      ENTRY WMPARL(KLINE)
C
      MODE=1
      RF=DBLE(CRF)
      RFI=DIMAG(CRF)
      KNAME=' &WM '//KLINE//' &END'
      WRITE(7,'(A80)') KNAME
      REWIND(7)
      READ(7,WM,ERR=8,END=8)
      WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
      GOTO 9
    8 CALL WMPLST
    9 REWIND(7)
      CRF=DCMPLX(RF,RFI)
      GOTO 3000
C
      ENTRY WMPARF
C
      MODE=2
      KPNAME='wmparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(.NOT.LEX) RETURN
C
      RF=DBLE(CRF)
      RFI=DIMAG(CRF)
      OPEN(7,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
      READ(7,WM,ERR=9800,END=9900)
      CLOSE(7)
      WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
      CRF=DCMPLX(RF,RFI)
C
 3000 IERR=0
C
      IF(NSMAX.LT.0.OR.NSMAX.GT.NSM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NSMAX'
         WRITE(6,*) '                  NSMAX,NSM =',NSMAX,NSM
         IERR=1
      ENDIF
C
      IF((NAMAX.LT.1).OR.(NAMAX.GT.NAM)) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NAMAX'
         WRITE(6,*) '                  NAMAX,NAM =',NAMAX,NAM         
         IERR=1
      ENDIF
C
      IF(NTHMAX.LT.1.OR.NTHMAX.GT.MDM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTHMAX'
         WRITE(6,*) '                  NTHMAX,MDM =',NTHMAX,MDM
         IERR=1
      ELSE
         MDP=NINT(LOG(DBLE(NTHMAX))/LOG(2.D0))
         IF(2**MDP.NE.NTHMAX) THEN
            WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTHMAX'
            WRITE(6,*) '                  NTHMAX,MDM =',NTHMAX,MDM
            IERR=1
         ENDIF
      ENDIF
C
      IF(NPHMAX.LT.1.OR.NPHMAX.GT.NDM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NPHMAX'
         WRITE(6,*) '                  NPHMAX,NDM =',NPHMAX,NDM
         IERR=1
      ELSE
         NDP=NINT(LOG(DBLE(NPHMAX))/LOG(2.D0))
         IF(2**NDP.NE.NPHMAX) THEN
            WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NPHMAX'
            WRITE(6,*) '                  NPHMAX,NDM =',NPHMAX,NDM
            IERR=1
         ENDIF
      ENDIF
C
      IF(NRMAX.LT.1.OR.NRMAX+1.GT.NRM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX,NRM =',NRMAX,NRM
         IERR=1
      ENDIF
C
      IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1000
C
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      RETURN
C
      END
C
C     ****** DISPLAY INPUT DATA ******
C
      SUBROUTINE WMVIEW
C
      INCLUDE 'wmcomm.h'
C
      IF(NPRINT.LT.2) RETURN
C
      IF(MODELG.EQ.0) THEN
         WRITE(6,*) '## 0: CYLINDRICAL ##'
      ELSE IF(MODELG.EQ.1) THEN 
         WRITE(6,*) '## 1: TOROIDAL ##'
      ELSE IF(MODELG.EQ.2) THEN 
         WRITE(6,*) '## 2: STRAIGHT HELICAL ##'
      ELSE IF(MODELG.EQ.3) THEN 
         WRITE(6,*) '## 3: TASK/EQ ##'
      ELSE IF(MODELG.EQ.4) THEN 
         WRITE(6,*) '## 4: VMEC ##'
      ELSE IF(MODELG.EQ.5) THEN 
         WRITE(6,*) '## 5: EQDSK ##'
      END IF
C
      IF(MODELJ.EQ.0) THEN
         WRITE(6,*) '## 0: REAL ANTENNA ##'
      ELSE IF(MODELJ.EQ.1) THEN 
         WRITE(6,*) '## 1: REAL ANTENNA ##'
      ELSE IF(MODELJ.EQ.2) THEN 
         WRITE(6,*) '## 2: POLOIDAL MODE ##'
      ELSE IF(MODELJ.EQ.3) THEN 
         WRITE(6,*) '## 3: TOROIDAL MODE ##'
      ELSE
         WRITE(6,*) '##',MODELJ,': VACUUM EIGEN MODE ##'
      END IF
C
      IF(MODELP.EQ.0) THEN
         WRITE(6,*) '## 0: VACUUM ##'
      ELSE IF(MODELP.EQ.1) THEN 
         WRITE(6,*) '## 1: MHD PLASMA ##'
      ELSE IF(MODELP.EQ.2) THEN 
         WRITE(6,*) '## 2: COLD PLASMA ##'
      ELSE IF(MODELP.EQ.3) THEN 
         WRITE(6,*) '## 3: HOT PLASMA (NO FLR) ##'
      ELSE IF(MODELP.EQ.4) THEN 
         WRITE(6,*) '## 4: HOT PLASMA (COLD FLR) ##'
      ELSE IF(MODELP.EQ.5) THEN 
         WRITE(6,*) '## 5: HOT PLASMA (FLR) ##'
      END IF
C
      IF(MODELN.EQ.8) THEN
         WRITE(6,*) '## 8: READ PROFILE DATA : WMDPRF ##'
      ELSEIF(MODELN.EQ.9) THEN
         WRITE(6,*) '## 9: READ PROFILE DATA : WMXPRF ##'
      ENDIF
C
      IF(MOD(MODELA,4).EQ.1) THEN
         WRITE(6,*) '## 1: ALPHA PARTICLE EFFECT ##'
      ELSE IF(MOD(MODELA,4).EQ.2) THEN
         WRITE(6,*) '## 2: ELECTRON BETA EFFECT ##'
      ELSE IF(MOD(MODELA,4).EQ.3) THEN
         WRITE(6,*) '## 3: ALPHA PARTICLE AND ELECTRON BETA EFFECTS ##'
      ENDIF
      IF(MODELA.GE.4) THEN
         WRITE(6,*) '## 4: ALPHA PARTICLE DENSITY CALCULATED ##'
      ENDIF
C
C      IF(MODELK.EQ.0) THEN
C         WRITE(6,*) '## 0: NO MODE NUMBER CUTOFF ##'
C      ELSE IF(MODELK.EQ.1) THEN
C         WRITE(6,*) '## 1: WITH MODE NUMBER CUTOFF ##'
C      ENDIF
C
      IF(MODELM.EQ.0) THEN
         WRITE(6,*) '## 0: BANDCD ##'
      ELSE IF(MODELM.EQ.1) THEN
         WRITE(6,*) '## 1: BCGCDB ##'
      ELSE IF(MODELM.EQ.2) THEN
         WRITE(6,*) '## 2: CGSCDB ##'
      ELSE IF(MODELM.EQ.3) THEN
         WRITE(6,*) '## 3: BSTABCDB ##'
      ELSE IF(MODELM.EQ.4) THEN
         WRITE(6,*) '## 4: BANDCDM ##'
      ELSE IF(MODELM.EQ.5) THEN
         WRITE(6,*) '## 5: BCGCDBM ##'
      ELSE IF(MODELM.EQ.6) THEN
         WRITE(6,*) '## 6: CGSCDBM ##'
      ELSE IF(MODELM.EQ.7) THEN
         WRITE(6,*) '## 7: BSTABCDBM ##'
      ELSE IF(MODELM.EQ.8) THEN
         WRITE(6,*) '## 8: BANDCDBM ##'
      ELSE IF(MODELM.EQ.9) THEN
         WRITE(6,*) '## 9: BCGCDBMA ##'
      ELSE IF(MODELM.EQ.10) THEN
         WRITE(6,*) '## 10: CGSCDBMA ##'
      ELSE IF(MODELM.EQ.11) THEN
         WRITE(6,*) '## 11: BSTABCDBMA ##'
      ELSE IF(MODELM.EQ.12) THEN
         WRITE(6,*) '## 12: BANDCDB'
      ENDIF
C
C      IF(MODELW.EQ.0) THEN
C         WRITE(6,*) '## 0: NOT WRITTING PABS. DATA ##'
C      ELSE IF(MODELW.EQ.1) THEN
C         WRITE(6,*) '## 1: WRITTING PABS. DATA ##'
C      ENDIF
C
      RF =DBLE(CRF)
      RFI=DIMAG(CRF)
      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    ,
     &             'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'Q0    ',Q0    ,'QA    ',QA    ,
     &             'RKAP  ',RKAP  ,'RDEL  ',RDEL
      WRITE(6,601) 'PROFN1',PROFN1,'PROFN2',PROFN2,
     &             'PROFT1',PROFT1,'PROFT2',PROFT2
      WRITE(6,601) 'ZEFF  ',ZEFF  ,'PNA   ',PNA   ,
     &             'PNAL  ',PNAL  ,'PTA   ',PTA
      WRITE(6,601) 'PROFU1',PROFU1,'PROFU2',PROFU2,
     &             'RHOMIN',RHOMIN,'QMIN  ',QMIN
      WRITE(6,601) 'RHOITB',RHOITB,'PRFIN ',PRFIN
      WRITE(6,601) 'RF    ',RF    ,'RFI   ',RFI   ,
     &             'RD    ',RD    ,'BETAJ ',BETAJ
      WRITE(6,602) 'NRMAX ',NRMAX ,'NTHMAX',NTHMAX,
     &             'NPHMAX',NPHMAX
      WRITE(6,602) 'NTH0  ',NTH0  ,'NPH0  ',NPH0  ,
     &             'NHC   ',NHC
C
      IF(MODELP.EQ.1.OR.MODELP.EQ.2) THEN
         WRITE(6,691)
         DO NS=1,NSMAX
            WRITE(6,610) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS)
         ENDDO
      ELSE IF(MODELP.GT.2) THEN
         WRITE(6,692)
         DO NS=1,NSMAX
            WRITE(6,611) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),
     &                      PTPR(NS),PTPP(NS),PTS(NS)
         ENDDO
         DO NS=1,NSMAX
            WRITE(6,612) NS,PU(NS),PUS(NS),
     &                   PNITB(NS),PTITB(NS),PUITB(NS)
         ENDDO
      ENDIF
C
      WRITE(6,693)
      DO NA=1,NAMAX
         WRITE(6,610) NA,AJ(NA),APH(NA),THJ1(NA),THJ2(NA),
     &                                  PHJ1(NA),PHJ2(NA)
      ENDDO
      RETURN
C
  601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  610 FORMAT(' ',I1,6(1PE11.3))
  611 FORMAT(' ',I1,7(1PE11.3))
  612 FORMAT(' ',I1,22X,5(1PE11.3))
  691 FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS',
     &                      8X,'PZCL')
  692 FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS',
     &                      8X,'PTPR',7X,'PTPP',7X,'PTS'/
     &       ' ','        ',11X,    9X,'PU',9X,'PUS',
     &                      7X,'PNITB',6X,'PTITB',5X,'PUITB')
  693 FORMAT(' ','NA    AJ',9X,'APH',8X,'THJ1',7X,'THJ2',
     &                      7X,'PHJ1',7X,'PHJ2')
      END
C
C     ****** BROADCAST PARAMETERS ******
C
      SUBROUTINE WMPRBC
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION IPARA(22),DPARA(28)
C
      IF(MYRANK.EQ.0) THEN
         RF=DBLE(CRF)
         RFI=DIMAG(CRF)
         IPARA(1) =NSMAX
         IPARA(2) =NAMAX
         IPARA(3) =NRMAX
         IPARA(4) =NTHMAX
         IPARA(5) =NPHMAX
         IPARA(6) =NTH0
         IPARA(7) =NPH0
         IPARA(8) =NHC
         IPARA(9) =NPRINT
         IPARA(10)=NGRAPH
         IPARA(11)=MODELG
         IPARA(12)=MODELJ
         IPARA(13)=MODELP
         IPARA(14)=MODELN
         IPARA(15)=MODELA
         IPARA(16)=MODELK
         IPARA(17)=MODELM
         IPARA(18)=MODELW
         IPARA(19)=LISTEG
         IPARA(20)=LMAXNW
         IPARA(21)=LISTNW
         IPARA(22)=MODENW
C
         DPARA(1) =BB
         DPARA(2) =RR
         DPARA(3) =RA
         DPARA(4) =RB
         DPARA(5) =Q0
         DPARA(6) =QA
         DPARA(7) =RKAP
         DPARA(8) =RDEL
         DPARA(9) =PROFN1
         DPARA(10)=PROFN2
         DPARA(11)=PROFT1
         DPARA(12)=PROFT2
         DPARA(13)=ZEFF
         DPARA(14)=PNA
         DPARA(15)=PNAL
         DPARA(16)=PTA
         DPARA(17)=RF
         DPARA(18)=RFI
         DPARA(19)=RD
         DPARA(20)=BETAJ
         DPARA(21)=DLTNW
         DPARA(22)=EPSNW
         DPARA(23)=RHOMIN
         DPARA(24)=QMIN
         DPARA(25)=RHOITB
         DPARA(26)=PROFU1
         DPARA(27)=PROFU2
         DPARA(28)=PRFIN
      ENDIF
C
      CALL MPBCIN(IPARA,22)
      CALL MPBCDN(DPARA,28)
C
      IF(MYRANK.NE.0) THEN
         NSMAX =IPARA(1) 
         NAMAX =IPARA(2) 
         NRMAX =IPARA(3) 
         NTHMAX=IPARA(4) 
         NPHMAX=IPARA(5) 
         NTH0  =IPARA(6) 
         NPH0  =IPARA(7) 
         NHC   =IPARA(8) 
         NPRINT=IPARA(9) 
         NGRAPH=IPARA(10)
         MODELG=IPARA(11)
         MODELJ=IPARA(12)
         MODELP=IPARA(13)
         MODELN=IPARA(14)
         MODELA=IPARA(15)
         MODELK=IPARA(16)
         MODELM=IPARA(17)
         MODELW=IPARA(18)
         LISTEG=IPARA(19)
         LMAXNW=IPARA(20)
         LISTNW=IPARA(21)
         MODENW=IPARA(22)
C       
         BB    =DPARA(1) 
         RR    =DPARA(2) 
         RA    =DPARA(3) 
         RB    =DPARA(4) 
         Q0    =DPARA(5) 
         QA    =DPARA(6) 
         RKAP  =DPARA(7) 
         RDEL  =DPARA(8) 
         PROFN1=DPARA(9) 
         PROFN2=DPARA(10)
         PROFT1=DPARA(11)
         PROFT2=DPARA(12)
         ZEFF  =DPARA(13)
         PNA   =DPARA(14)
         PNAL  =DPARA(15)
         PTA   =DPARA(16)
         RF    =DPARA(17)
         RFI   =DPARA(18)
         RD    =DPARA(19)
         BETAJ =DPARA(20)
         DLTNW =DPARA(21)
         EPSNW =DPARA(22)
         RHOMIN=DPARA(23)
         QMIN  =DPARA(24)
         RHOITB=DPARA(25)
         PROFU1=DPARA(26)
         PROFU2=DPARA(27)
         PRFIN =DPARA(28)
         CRF=DCMPLX(RF,RFI)
      ENDIF
C
      CALL MPBCDN(PA,NSMAX)
      CALL MPBCDN(PZ,NSMAX)
      CALL MPBCDN(PN,NSMAX)
      CALL MPBCDN(PNS,NSMAX)
      CALL MPBCDN(PZCL,NSMAX)
      CALL MPBCDN(PTPR,NSMAX)
      CALL MPBCDN(PTPP,NSMAX)
      CALL MPBCDN(PTS,NSMAX)
      CALL MPBCDN(PU,NSMAX)
      CALL MPBCDN(PUS,NSMAX)
      CALL MPBCDN(PNITB,NSMAX)
      CALL MPBCDN(PTITB,NSMAX)
      CALL MPBCDN(PUITB,NSMAX)
      CALL MPBCDN(AJ,NAMAX)
      CALL MPBCDN(APH,NAMAX)
      CALL MPBCDN(THJ1,NAMAX)
      CALL MPBCDN(THJ2,NAMAX)
      CALL MPBCDN(PHJ1,NAMAX)
      CALL MPBCDN(PHJ2,NAMAX)
      CALL MPBCKN(KNAMEQ,80)
      CALL MPBCKN(KNAMPF,80)
C
      RETURN
      END
