C     $Id$
C
C     ***************************************************************
C
C        Set constants and initial parameters
C
C     ***************************************************************
C
      SUBROUTINE TXINIT
C
      INCLUDE 'txcomm.inc'
C
C     ***** Version ID *****
C     SLID is used to identify data file.
      SLID = 'tx211.00'
C
C     ***** Physical constants *****
C
C     Electron charge (C)
      AEE  = 1.60217733D-19
C
C     Electron mass (kg)
      AME  = 9.1093897D-31
C
C     Proton mass (kg)
      AMP  = 1.6726231D-27
C
C     Light velocity (m/s)
      VC   = 2.99792458D8
C
C     Pi
      PI   = 2 * ASIN(1.D0)
C
C     mu0 (H/m)
      rMU0 = 4 * PI * 1.D-7
C
C     epsilon0 (F/m)
      EPS0 = 1.D0 / (rMU0 * VC**2)
C
C     Conversion factor from keV to joule
      rKEV = 1.D3 * AEE
C
C     ***** Configuration parameters *****
C
C     Plasma minor radius (m)
      RA = 0.35D0
C
C     Wall radius (m)
      RB = 0.4D0
C
C     Plasma major radius (m)
      RR = 1.3D0
C
C     Toroidal magnetic field (T)
      BB = 1.3D0
C
C     Plasma current start (MA)
      rIPs= 0.150D0
C
C     Plasma current end (MA)
      rIPe= 0.150D0
C
C     ***** Plasma components *****
C
C     Atomic number of ion
      PA = 1.D0
C
C     Charge number of ion
      PZ = 1.D0
C
C     Effective charge
      Zeff = 2.D0
C
C     ***** Initial plasma parameters *****
C
C     Initial electron density at r = 0 (10^20 m^-3)
      PN0 = 0.4D0
C
C     Initial electron density at r = a (10^20 m^-3)
      PNa = 0.05D0
C
C     Initial electron temperature at r = 0 (keV)
      PTe0 = 700.D-3
C
C     Initial electron temperature at r = a (keV)
      PTea =  50.D-3
C
C     Initial ion temperature  at r = 0 (keV)
      PTi0 = 700.D-3
C
C     Initial ion temperature  at r = a (keV)
      PTia =  50.D-3
C
C     Initail current profile parameter
      PROFJ = 2.D0
C
C     ***** Particle diffusivity and viscosity parameters *****
C
C     Electron-driven diffusion parameter
      De0 = 0.D0
C
C     Ion-driven diffusion parameter
      Di0 = 1.D0
C
C     Electron viscosity parameter
      rMue0 = 3.D0
C
C     Ion viscosity parameter
      rMui0 = 3.D0
C
C     Drift frequency parameter (omega/omega*e)
      WPM0 = 0.D0
C
C    ***** Thermal diffusivity parameters *****
C
C     Electron thermal diffusivity parameter (Chie/D)
C       0 for fixed temperature profiel
      Chie0 = 0.D0
C
C     Ion thermal diffusivity parameter (Chie/D)
C       0 for fixed temperature profiel
      Chii0 = 0.D0
C
C     ***** Turbulent transport control parameters *****
C
C     Fixed transport coefficient parameter
      FSDFIX = 1.0D0
C
C     CDBM transport coefficient parameter
C       (Current diffusive ballooning mode)
      FSCDBM = 0.D0
C
C     Bohm transport coefficient parameter in SOL
      FSBOHM = 0.D0
C
C     Pseud-classical transport coefficient parameter in SOL
      FSPSCL = 0.D0
C
C     Diffusion coefficient profile parameter (D(r=a)/D(r=0))
      PROFD = 10.D0
C
C     ***** Other transport parameters *****
C
C     Charge exchange parameter
      FSCX = 1.D0
C
C     Orbit loss parameter
CCCC      FSLC = 1.D0
      FSLC = 0.D0
C
C     Toroidal Neoclassical viscosity parameter
      FSNC = 1.D0
C
C     Helical Neoclassical viscosity parameter
      FSHL = 0.D0
C
C     Particle loss to divertor parameter
      FSLP = 1.D0
C
C     Ionization parameter
      FSION = 1.D0
C
C     Neutral diffusion factor
      FSD0 = 1.D0
C
C     Factor of E x B rotation shear
      rG1 = 24.D0
C
C     ***** initial parameters *****
C
C     Initial Density scale length in SOL (m)
      rLn = 0.03D0
C
C     Initail Temperature scale length in SOL (m)
      rLT = 0.030D0
C
C     ***** Heating parameters *****
C
C     NBIl beam energy (keV)
      Eb = 32
C
C     Heating radius of NBI heating (m)
      RNB = 0.175D0
C
C     NBI input power (MW)
      PNBH = 0.D0
C
C     NBI current drive parameter
      PNBCD= 0.D0
C
C     Refractive index of RF waves
      rNRF = 0.D0
C
C     Heating radius of RF heating (m)
      RRF = 0.175D0
C
C     RF input power (MW)
      PRFH = 0.D0
C
C     ***** Neutral parameters *****
C
C     Initial Neutral density (10^20 m^-3)
      PN0s = 1.D-8
C
C     Neutral thermal velocity (m/s)
      V0 = 1.5D3
C
C     Recycling rate in SOL
      rGamm0 = 0.8D0
C
C     Gas-puff particle flux (10^20 1/s)
      rGASPF = 0.1D0
C
C     Electron density in diverter region (Minimum density in SOL)
      PNeDIV = 1.D-2
C
C     Electron temperature in diverter region (Minimum Te in SOL)
      PTeDIV = 1.D-2
C
C     Ion temperature in diverter region (Minimum Ti in SOL)
      PTiDIV = 1.D-2
C
C     ***** Numerical paramters *****
C
C     Implicitness parameter (Not used now)
      DLT = 1.0D0
C
C     Time step size(s)
CCCC      DT = 1.D-4
      DT = 1.D-3
C
C     Convergence parameter
      EPS = 1.D99
C
C     Iteration
      ICMAX=2
C
C     ***** Mesh number parameters *****
C
C     Radial step number
      NRMAX = 40
C
C     Number of time step
      NTMAX = 10
C
C     Time step interval between print output
      NTSTEP = 10
C
C     Time step interval between lines in f(r) graph
      NGRSTP = 1
C
C     Time step interval between points in f(t) graph
      NGTSTP = 1
C
C     Time step interval between points in f(t) graph
      NGVSTP = 1
C
C     Mode of Graph
C     1 : for Display
C     2 : for Print Out
      MODEG = 2
C
C     MODE of Graph Line
C     0 : Change Line Color (Last Color Fixed)
C     1 : Change Line Color and Style
C     2 : Change Line Color, Style and Mark
C     3 : Change Line Color, Style and Mark (With Legend)
      MODEl=1
C
C     Mode of AV
C     0 : OFF
C     n : Number of Display
      MODEAV = 0
C
      DO I = 1, NGYRM
         gDIV(I) = 1.E0
      ENDDO
      gDIV(1)  = 1.E20
      gDIV(2)  = 1.E14
      gDIV(4)  = 1.E3
      gDIV(5)  = 1.E3
      gDIV(7)  = 1.E3
      gDIV(8)  = 1.E3
      gDIV(9)  = 1.E3
      gDIV(16) = 1.E14
      gDIV(18) = 1.E6
C
C     Radius where density increase by command DEL
      DelR = 0.175D0
C
C     Amount of increase of density by command DEL
      DelN = 5.D-1
C
C     Helical ripple amplitude at r=a, Bhelical/Btoroidal, linear to r/a
      EpsH = 0.1D0
C
C     Helical pitch number
      NCphi = 10
C
C     Safety factor for helical
      Q0 = 3.D0
      QA = 2.D0
C
      NGR=-1
C
      RETURN
      END
C
C     ***************************************************************
C
C        Change input parameters
C
C     ***************************************************************
C
      SUBROUTINE TXPARM(KID)
C
      INCLUDE 'txcomm.inc'
C
      NAMELIST /TX/
     & RA,RB,RR,BB,
     & PA,PZ,Zeff,
     & PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,
     & De0,Di0,rMue0,rMui0,WPM0,
     & Chie0,Chii0,
     & FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD,
     & FSCX,FSLC,FSNC,FSLP,FSION,FSD0,
     & rLn,rLT,
     & Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH,
     & PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,
     & DLT,DT,EPS,ICMAX,
     & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,
     & DelR,DelN,
     & rG1,EpsH,FSHL,NCPHI,Q0,QA,
     & rIPs,rIPe,
     & MODEG, gDIV, MODEAV, MODEl,
     & FSHL,EpsH, NCphi
C
      LOGICAL LEX
      CHARACTER KPNAME*80,KLINE*80,KNAME*90,KID*1
C
      MODE=0
    1    CONTINUE
         WRITE(6,*) '# INPUT &TX :'
         READ(5,TX,ERR=2,END=3)
         KID=' '
         GOTO 4
C
    2    CALL TXPLST
      GOTO 1
C
    3 KID='Q'
    4 GOTO 3000
C
      ENTRY TXPARL(KLINE)
C
      MODE=1
      KNAME=' &TX '//KLINE//' &END'
      WRITE(7,'(A90)') KNAME
      REWIND(7)
      READ(7,TX,ERR=8,END=8)
      WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
      GOTO 9
    8 CALL TXPLST
    9 REWIND(7)
      GOTO 3000
C
      ENTRY TXPARF(KPNAME)
C
      MODE=2
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(.NOT.LEX) RETURN
C
      OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
      READ(25,TX,IOSTAT=IST,ERR=9800,END=9900)
      CALL KTRIM(KPNAME,KL)
      WRITE(6,*) 
     &     '## FILE (',KPNAME(1:KL),') IS ASSIGNED FOR PARM INPUT'
C
 3000 IERR=0
C
C     ERROR CHECK
C
      IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1
C
      RETURN

 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      RETURN

      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE TXPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &TX : RA,RB,RR,BB,PA,PZ,Zeff,'/
     &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,'/
     &       ' ',8X,'De0,Di0,rMue0,rMui0,WPM0,'/
     &       ' ',8X,'Chie0,Chii0,'/
     &       ' ',8X,'FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD,'/
     &       ' ',8X,'FSCX,FSLC,FSNC,FSLP,FSION,FSD0,'/
     &       ' ',8X,'rLn,rLT,'/
     &       ' ',8X,'Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH,'/
     &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,'/
     &       ' ',8X,'DLT,DT,EPS,ICMAX,'/
     &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,'/
     &       ' ',8X,'DelR,DelN,'/
     &       ' ',8X,'rG1,EpsH,FSHL,NCPHI,Q0,QA,'/
     &       ' ',8X,'rIPs,rIPe,'/
     &       ' ',8X,'MODEG, gDIV, MODEAV, MODEl,'/
     &       ' ',8X,'FSHL,EpsH, NCphi')
      END
C
C     ***************************************************************
C
C        View input parameters
C
C     ***************************************************************
C
      SUBROUTINE TXVIEW
C
      INCLUDE 'txcomm.inc'
C
      WRITE(6,'((1H ,A6,2H =,1PD9.2,3(2X,A6,2H =,1PD9.2)))')
     &   'RA    ', RA    ,  'RB    ', RB    , 
     &   'RR    ', RR    ,  'BB    ', BB    , 
     &   'PA    ', PA    ,  'PZ    ', PZ    , 
     &   'PN0   ', PN0   ,  'PNa   ', PNa   , 
     &   'PTe0  ', PTe0  ,  'PTea  ', PTea  , 
     &   'PTi0  ', PTi0  ,  'PTia  ', PTia  , 
     &   'rIP   ', rIP   ,  'Zeff  ', Zeff  , 
     &   'PROFJ ', PROFJ ,
     &   'De0   ', De0   ,  'Di0   ', Di0   , 
     &   'rMue0 ', rMue0 ,  'rMui0 ', rMui0 , 
     &   'WPM0  ', WPM0  ,  'PROFD ', PROFD , 
     &   'Chie0 ', Chie0 ,  'Chii0 ', Chii0 , 
     &   'FSDFIX', FSDFIX,  'FSCDBM', FSCDBM, 
     &   'FSBOHM', FSBOHM,  'FSPSCL', FSPSCL, 
     &   'FSCX  ', FSCX  ,  'FSLC  ', FSLC  , 
     &   'FSNC  ', FSNC  ,  'FSLP  ', FSLP  , 
     &   'FSION ', FSION ,  'FSD0  ', FSD0  , 
     &   'rLn   ', rLn   ,  'rLT   ', rLT   , 
     &   'Eb    ', Eb    ,  'RNB   ', RNB   , 
     &   'PNBH  ', PNBH  ,  'rNRF  ', rNRF  , 
     &   'RRF   ', RRF   ,  'PRFH  ', PRFH  , 
     &   'rGamm0', rGamm0,  'V0    ', V0    , 
     &   'rGASPF', rGASPF,  'PNeDIV', PNeDIV, 
     &   'PTeDIV', PTeDIV,  'PTiDIV', PTiDIV, 
     &   'PN0s  ', PN0s  ,  'DLT   ', DLT   , 
     &   'EPS   ', EPS   ,  'DT    ', DT    ,
     &   'rG1   ', rG1   ,  'Zeff  ', Zeff  , 
     &   'rIPs  ', rIPs  ,  'rIPe  ', rIPe  ,
     &   'FSHL  ', FSHL  ,  'EpsH  ', EpsH  ,
     &   'Q0    ', Q0    ,  'QA    ', QA
      WRITE(6,'((1H ,A6,2H =,I5,3(6X,A6,2H =,I5)))')
     &   'NRMAX ', NRMAX , 
     &   'NTMAX ', NTMAX ,  'NTSTEP', NTSTEP, 
     &   'NGRSTP', NGRSTP,  'NGTSTP', NGTSTP, 
     &   'NGVSTP', NGVSTP,  'ICMAX ', ICMAX ,
     &   'MODEG ', MODEG ,  'MODEAV', MODEAV,
     &   'MODEl ', MODEl ,  'NCPHI ', NCPHI
C
      RETURN
      END
C
C     ***************************************************************
C
C        Initialize profiles
C
C     ***************************************************************
C
      SUBROUTINE TXPROF
C
      INCLUDE 'txcomm.inc'
C
      CALL TXCALM
      CALL TXPRFG
C
C  Initialize variable vector
C
      DO NR = 0, NRMAX
         DO NQ = 1, NQMAX
            X(NQ,NR) = 0
         ENDDO
      ENDDO
C
C  Half integer mesh variables
C
      DO NR = 0, NRMAX - 1
         RL=RHI(NR)
         IF (RL .LT. RA) THEN
            PROF  = 1.D0 - (RL / RA)**2
            PROFT = PROF**2
            X(LQe1,NR)  = (PN0 - PNa) * PROF + PNa
            X(LQi1,NR)  = X(LQe1,NR) / PZ
            X(LQe5,NR) = ((PTe0 - PTea) * PROFT + PTea)*X(LQe1,NR)
            X(LQi5,NR) = ((PTi0 - PTia) * PROFT + PTia)*X(LQi1,NR)
         ELSE
            X(LQe1,NR)  = PNa * EXP(-(RL-RA) / rLn)
            X(LQi1,NR)  = X(LQe1,NR) / PZ
CCCC            X(LQe5,NR) = PTea * EXP(-(RL-RA) / rLT)
CCCC            X(LQi5,NR) = PTia * EXP(-(RL-RA) / rLT)
            X(LQe5,NR) = PTea*X(LQe1,NR)
            X(LQi5,NR) = PTia*X(LQi1,NR)
         ENDIF
         X(LQn1,NR) = PN0s
         X(LQn2,NR) = 0.D0
         X(LQm5,NR) = BB
C
         IF((1.D0-(RHI(NR)/RA)**PROFJ).LE.0.D0) THEN
            PROF=0.D0    
         ELSE             
            PROF= (1.D0-(RHI(NR)/RA)**PROFJ)
         ENDIF             
         IF(FSHL.EQ.0.D0) THEN
            X(LQe4,NR) = - rIps * 1.D6 / (AEE * PI * RA**2 * 1.D20)
     &                     * (PROFJ + 1) * PROF**PROFJ
            AJOH(NR)= PROF
         ELSE
            X(LQe4,NR) = 0.D0
            AJOH(NR)= 0.D0
         ENDIF
      ENDDO
C
C Integer mesh variables
C
      IF(FSHL.EQ.0.D0) THEN
         X(LQm4,0) = 0.D0
         DO NR = 1, NRMAX
            RL=R(NR)
            IF (RL .LT. RA) THEN
               PROF = 1.D0 - (RL / RA)**2
               X(LQm4,NR) = rMU0 * rIps * 1.D6 / (2 * PI * RL)
     &              * (1 - PROF**(PROFJ+1))
            ELSE
               X(LQm4,NR) = rMU0 * rIps * 1.D6 / (2 * PI * RL)
            ENDIF
         ENDDO
      ELSE
         X(LQm4,0) = 0.D0
         DO NR = 1, NRMAX
            RL=R(NR)
            QL=(Q0-QA)*(1-(RL/RA)**2)+QA
            X(LQm4,NR) = BB*RL/(QL*RR)
         ENDDO
      ENDIF
C
      IF(FSHL.EQ.0.D0) THEN
         DO NR=0,NRMAX
            AJV(NR)=0.D0
         ENDDO
      ELSE
         DO NR=0,NRMAX
            RL=R(NR)
            RIP1=2.D0*PI*RL*X(LQm4,NR)/rMU0
            RL=R(NR+1)
            RIP2=2.D0*PI*RL*X(LQm4,NR+1)/rMU0
            RL=RHI(NR)
            AJV(NR)=(RIP2-RIP1)/(2.D0*PI*RL*DR)
         ENDDO
      ENDIF
C
      DO NR = 0, NRMAX - 1
         RL=RHI(NR)
         IF (RL .LT. RA) THEN
            PROF = (1.D0 - (RL / RA)**2)**PROFJ
         ELSE
            PROF = 0.D0
         ENDIF
         PNeHI(NR)=X(LQe1,NR)
         PTeHI(NR)=X(LQe5,NR)
         rLnLam = + 15 - LOG(ABS(PNeHI(NR))) / 2 + LOG(ABS(PTeHI(NR)))
         ETA =  SQRT(AME) * Zeff  * AEE**2 * rLnLam
     &              / (3 * (2 * PI)**1.5D0
     &                   * EPS0**2
     &                   * (ABS(PTeHI(NR)) * rKeV)**1.5D0)
         rJP = rIps * 1.D6 / (PI * RA**2) * (PROFJ + 1) * PROF
         IF(FSHL.EQ.0.D0) THEN
            X(LQm3,NR) = ETA * rJP
         ELSE
            X(LQm3,NR) = 0.D0
         ENDIF
      ENDDO
C
      TIME=0.D0
      NGT=-1
      NGR=-1
      NGVV=-1
      rIP=rIPs
C
      CALL TXCALV(X)
      CALL TXCALC
      CALL TXCALA
      CALL TXGLOB
      CALL TXSTGT(SNGL(TIME))
      CALL TXSTGV(SNGL(TIME))
      CALL TXSTGR
C
      RETURN
      END
