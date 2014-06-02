MODULE T2COMM
  
  USE T2CNST, ONLY: rkind, ikind, i0lmaxm, i0spcsm, twopi
  
  IMPLICIT NONE

  !---- input parameters added by AF
  INTEGER(ikind):: &
       ntmax,ntstep,nconvmax, &
       nt0dmax,nt0dstep, &
       nt2dmax,nt2dstep, &
       idfile,idprint,idplot,idmode,idebug
  REAL(rkind):: &
       Dt,time_init,eps_conv

  !---- global parameters added by AF
  INTEGER(ikind):: nrhomax ! number of variables in rho
  INTEGER(ikind):: nchimax ! number of variables in chi
  INTEGER(ikind):: nequmax ! number of equations to be solved
  INTEGER(ikind):: nv0dmax ! number of global variables to be saved
  INTEGER(ikind):: nv2dmax ! number of profile varibales to be saved
  REAL(   rkind):: &
       time_t2               ! global time
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       vv0d,vv2d             ! storage for 0D and 2D data
  
  !C---------------------------------------------------------
  !C
  !C
  !C  GLOBAL PARAMETERS
  !C 
  !C
  !C---------------------------------------------------------
  
  INTEGER(ikind)::i0dbg
  
  INTEGER(ikind) ::&
       i0vmax, & !C NUMBER OF DEPENDENT VARIABLES
       i0wmax, & !C NUMBER OF WORKING VARIABLES FOR DIFFERENTIAL
       i0qmax, & !C NUMBER OF INTEGRAL POINTS PAR DIRECTION
       i0dmax, & !C NUMBER OF DIMENSIONS
       i0smax, & !C NUMBER OF SPECIES
       i0pmax, & !C MAXMUM NUMBER OF PICARD ITERATION LOOP 
       i0nmax, & !C NUMBER OF NODES PAR ELEMENT
       i0mmax, & !C NUMBER OF NODES W   OVERLAP IN DOMAIN
       i0xmax, & !C NUMBER OF NODES W/O OVERLAP W INTERPOLATION IN DOMAIN 
       i0bmax, & !C NUMBER OF NODES W/O OVERLAP AND INTERPOLATION IN DOMAIN
       i0rmax, & !C NUMBER OF NODES IN RADIAL DIRECTION (FOR 1D)
       i0emax, & !C NUMBER OF ELEMENTS IN DOMAIN
       i0hmax, & !C NUMBER OF INTERPOLATION NODES IN DOMAIN
       i0lmax, & !C NUMBER OF SUB DOMAINS 
       i0amax, & !C NUMBER OF NON-ZERO COMPONENTS OF MATRIX (CRS-METHOD)
       i0nrmx, & !C ARRAY SIZE OF I1NIDR (CRS-METHOD) 
       i0ermx, & !C ARRAY SIZE OF I1EIDR (CRS-METHOD) 
       i0ecmx, & !C ARRAY SIZE OF I1EIDC (CRS-METHOD) 
       i0pdiv_number

  !---------------------------------------------------------
  !
  !
  !       GLOBAL PARAMETERS (NEW)
  ! 
  !
  !---------------------------------------------------------
  INTEGER(ikind),SAVE::&
       NNMAX,&
       NQMAX,&
       NDMAX,&
       NSMAX,&
       NVMAX,&
       NRMAX,&
       NKMAX,& 
       NEMAX,&
       NBMAX,&
       NXMAX,&
       NMMAX,&
       NHMAX,&
       NAMAX,&
       NNRMX,&
       NERMX,&
       NECMX,&
       NLMAX,&
       NPMAX,&
       NPMIN,&
       NAVMX,& !    i0avmax = i0amax*i0vmax*i0vmax
       NBVMX   !    i0bvmax = i0bmax*i0vmax

  INTEGER(ikind)::&
       i0mfcs, & !C INDICATOR FOR COORDINATE SYSTEM (1: torus coordinate)
       i0supg, & !C INDICATOR FOR SUPG METHOD (0: w/o SUPG, 1: w SUPG)
       i0wstp, & !C INDICATOR FOR RESULT OUTPUT TIMING
       i0solv, & !C INDICATOR FOR SOLVED DEPENDENT VARIAVLES
       i0anom,&  !C INDICATOR FOR ANOMALOUS TRANSPORT
                 !C      0: w/o anomalous transport 
                 !C      1: w   anomalous transport
       i0cchk    !C INDICATOR FOR COEFFICIENT CHECK (1: coefs check)
  INTEGER(ikind)::&
       i0bvmax,& !C VECTOR SIZE OF b FOR MTXP (Ax=b)
       i0avmax   !C NUMBER OF NONZERO COMPONENT OF A FOR MTXP (Ax=b)

  REAL(   rkind)::&
       d0mfcst, d0mffct,& !C NORMALIZATION CONSTANTS FOR \psi'
       d0btcst, d0btfct,& !C NORMALIZATION CONSTANTS FOR I
       d0etcst, d0etfct,& !C NORMALIZATION CONSTANTS FOR E_{\zeta}
       d0epcst, d0epfct,& !C NORMALIZATION CONSTANTS FOR \bar{E}_{\chi }
       d0ercst, d0erfct,& !C NORMALIZATION CONSTANTS FOR E_{\rho }
       d0nncst, d0nnfct,& !C NORMALIZATION CONSTANTS FOR n_{a}
       d0frcst, d0frfct,& !C NORMALIZATION CONSTANTS FOR n_{a}\bar{u}_{a}^{\rho}
       d0fbcst, d0fbfct,& !C NORMALIZATION CONSTANTS FOR n_{a}u_{a\para}
       d0ftcst, d0ftfct,& !C NORMALIZATION CONSTANTS FOR n_{a}u_{a\zeta}
       d0fpcst, d0fpfct,& !C NORMALIZATION CONSTANTS FOR n_{a}u_{a\zeta}
       d0ppcst, d0ppfct,& !C NORMALIZATION CONSTANTS FOR p_{a}
       d0qrcst, d0qrfct,& !C NORMALIZATION CONSTANTS FOR \bar{Q}_{a}^{\rho}
       d0qbcst, d0qbfct,& !C NORMALIZATION CONSTANTS FOR Q_{a\para}
       d0qtcst, d0qtfct,& !C NORMALIZATION CONSTANTS FOR Q_{a\zeta}
       d0qpcst, d0qpfct,& !C NORMALIZATION CONSTANTS FOR Q_{a}^{\chi}
       d0ubcst, & !C NORMALIZATION CONSTANT FOR u_{a\para}
       d0wbcst, & !C NORMALIZATION CONSTANT FOR w_{a\para}
       d0rmjr,  & !C MAJOR RADIUS (R_{0} [m])
       d0rmnr,  & !C MINOR RADIUS (a     [m])
       d0iar,   & !C INVERSE ASPECT RATIO (a/R_{0})
       d0eps      !C CONVERGENCE CRITERION FOR PICARD ITERATION

  REAL(   rkind),SAVE::&
       BpNF,&
       BtNF,&
       EtNF,&
       EpNF,&
       ErNF,&
       NnNF,&
       FrNF,&
       FbNF,&
       FtNF,&
       FpNF,&
       PpNF,&
       QrNF,&
       QbNF,&
       QtNF,&
       QpNF
  REAL(   rkind),SAVE::&
       EqFaraday,&
       EqAmpere,&
       EqGauss,&
       EqConti,&
       EqMotion,&
       EqEnergy,&
       EqHFlux
       

  !-------------------------------------------------------------------
  !
  !       DEFINITION OF GLOBAL VARIABLES  FOR T2INTG
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------

  !
  ! INTEGRATION ARRAYS BY GAUSSIAN INTEGRATION
  !
  
  !  MassScaIntgPG: INTEGRATION ARRAY FOR MASS SCALAR SUBMATRIX
  !  AdveVecIntgPG: INTEGRATION ARRAY FOR ADVE VECTOR SUBMATRIX
  !  AdveTenIntgPG: INTEGRATION ARRAY FOR ADVE TENSOR SUBMATRIX
  !  DiffTenIntgPG: INTEGRATION ARRAY FOR DIFF TENSOR SUBMATRIX
  !  GradVecIntgPG: INTEGRATION ARRAY FOR GRAD VECTOR SUBMATRIX
  !  GradTenIntgPG: INTEGRATION ARRAY FOR GRAD TENSOR SUBMATRIX
  !  ExciScaIntgPG: INTEGRATION ARRAY FOR EXCI SCALAR SUBMATRIX
  !  ExciVecIntgPG: INTEGRATION ARRAY FOR EXCI VECTOR SUBMATRIX
  !  ExciTenIntgPG: INTEGRATION ARRAY FOR EXCI TENSOR SUBMATRIX
  !  SourScaIntgPG: INTEGRATION ARRAY FOR SOUR SCALAR SUBMATRIX

  REAL(   rkind),ALLOCATABLE,SAVE::&
       MassScaIntgPG(:,:,:        ), AdveVecIntgPG(:,:,:,:      ),&
       AdveTenIntgPG(:,:,:,:,:,:  ), DiffTenIntgPG(:,:,:,:,:    ),&
       GradVecIntgPG(:,:,:,:      ), GradTenIntgPG(:,:,:,:,:,:  ),&
       ExciScaIntgPG(:,:,:        ), ExciVecIntgPG(:,:,:,:,:    ),&
       ExciTenIntgPG(:,:,:,:,:,:,:), SourScaIntgPG(:,:          )


  !-------------------------------------------------------------------
  !
  !       DEFINITION OF GLOBAL VARIABLES  FOR T2VGRA
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------
  LOGICAL,SAVE::&
       UsePotentialDescription,&
       UseNormalization,&
       UseAnomalousTransportFT,&
       UseAnomalousTransportGT,&
       !
       SolveField,&
       SolveElectron,&
       SolveIons,&
       SolveDensity,&
       SolveFlux,&
       SolvePressure,&
       SolveHeatFlux,&
       !
       LockPoloidalMageticFieldOnAxis,&
       LockToroidalMageticFieldOnAxis,&
       LockRadialElectricFieldOnAxis,&
       LockPoloidalElectricFieldOnAxis,&
       LockToroidalElectricFieldOnAxis,&
       LockDensityOnAxis,&
       LockRaidalFluxOnAxis,&
       LockParallelFluxOnAxis,&
       LockToroidalFluxOnAxis,&
       LockPoroidalFluxOnAxis,&
       LockPressureOnAxis,&
       LockRaidalHeatFluxOnAxis,&
       LockParallelHeatFluxOnAxis,&
       LockToroidalHeatFluxOnAxis,&
       LockPoroidalHeatFluxOnAxis,&
       !
       LockPoloidalMageticFieldOnWall,&
       LockToroidalMageticFieldOnWall,&
       LockRadialElectricFieldOnWall,&
       LockPoloidalElectricFieldOnWall,&
       LockToroidalElectricFieldOnWall,&
       LockDensityOnWall,&
       LockRaidalFluxOnWall,&
       LockParallelFluxOnWall,&
       LockToroidalFluxOnWall,&
       LockPoroidalFluxOnWall,&
       LockPressureOnWall,&
       LockRaidalHeatFluxOnWall,&
       LockParallelHeatFluxOnWall,&
       LockToroidalHeatFluxOnWall,&
       LockPoroidalHeatFluxOnWall

  LOGICAL,SAVE,ALLOCATABLE::&
       HaveMassScaCoef(:,:    ),HaveAdveVecCoef(:,:    ),&
       HaveAdveTenCoef(:,:    ),HaveDiffTenCoef(:,:    ),&
       HaveGradVecCoef(:,:    ),HaveGradTenCoef(:,:    ),&
       HaveExciScaCoef(:,:    ),HaveExciVecCoef(:,:    ),&
       HaveExciTenCoef(:,:    ),HaveSourScaCoef(:,:    ),&
       !
       HaveAdveTenKval(:,:,:  ),HaveGradTenKval(:,:,:  ),&
       HaveExciVecKval(:,:,:  ),HaveExciTenKval(:,:,:,:),&
       !
       HaveMat(:,:), LockEqs(:),  LockAxi(:),  LockWal(:)
  
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2NGRA
  !C
  !C------------------------------------------------------------------ 

  !C
  !C ARRAYS ALLOCATED   BY T2NGRA_ALLOCATE
  !C        DEALLOCATED BY T2NGRA_DEALLOCATE
  !C

  !C I1MMAX : NUMBER OF NODES IN EACH SUBDOMAIN
  !C I1BMAX : NUMBER OF NODES IN EACH SUBDOMAIN (W/O OVERLAP)
  !C I1EMAX : NUMBER OF ELEMENTS IN EACH SUBDOMAIN
  !C I1MLEL : MESH LEVEL OF EACH SUBDOMAIN
  !C I1RDN1 : NUMBER OF PARTITION IN RADIAL   DIREACTON IN EACH SUBDOMAIN 
  !C I1RDN2 : NUMBER OF PARTITION IN RADIAL   DIREACTON IN EACH SUBDOMAIN 
  !C          (W/O OVERLAP)
  !C I1PDN1 : NUMBER OF PARTITION IN POLOIDAL DIREACTON IN EACH SUBDOMAIN
  !C I1PDN2 : NUMBER OF PARTITION IN POLOIDAL DIREACTON IN EACH SUBDOMAIN
  !C          (W/O OVERLAP)
  !C

  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::& 
       i1mmax,i1bmax,i1emax,i1mlel,i1rdn1,i1pdn1,i1pdn2

  !C
  !C ARRAYS ALLOCATED   BY T2COMM_ALLOCATE
  !C        DEALLOCATED BY T2COMM_DEALLOCATE
  !C
  
  !C
  !C FOR COMPRESSED ROW STORAGE METHOD FOR NODE-NODE RERATIONS
  !C
  !C I1NDIR: CUMULTATIVE NUMBER OF NODE-NODE CONNECTIVITY + 1 
  !C         UP TO EACH ROW IN MATRIX A
  !C I1NDIC: NODE NUMBER 
  !C
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::& 
       i1nidr,&
       i1nidc

  !C
  !C FOR COMPRESSED ROW STORAGE METHOD FOR NODE-ELEMENT RERATIONS
  !C
  !C I1EDIR: CUMULTATIVE NUMBER OF NODE-ELEMENT CONNECTIVITY + 1
  !C         UP TO EACH NODE 
  !C I1EDIC: ELEMENT NUMBER 
  !C
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::&
       i1eidr, i1eidc

  !C I2CRT: INFORMATION NODE PROJECTION MAP
  !C        1: FOR COEFFICIENT CALCULATION
  !C        2: FOR 2D
  !C        3: FOR 1D
  !C
  !C I3ENR: INFORMATION OF ELEMENT-NODE CONNECTIVITY 
  !C        1: FOR COEFFICIENT CALCULATION
  !C        2: FOR 2D 
  !C        3: FOR 1D-2D 
  !C        4: FOR 1D
  !C
  !C I2HBC: INFORMATION OF BOUND NODE AND BINDING NODES FOR h-MESH
  !C
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:,:,:)::i3enr 
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:,:  )::i2crt,i2hbc
  
  !C
  !C FOR 1D PROBLEM
  !C
  !C I1MC1D: NODE   NUMBERS   FOR DATA STOCK OF 1D VALUES
  !C D1MC1D: RADIAL POSITIONS FOR DATA STOCK OF 1D VALUES
  !C 

  INTEGER(ikind),ALLOCATABLE,DIMENSION(:):: i1mc1d 
  REAL(   rkind),ALLOCATABLE,DIMENSION(:):: d1mc1d



  !C
  !C D1MSIZ: SIZE OF ELEMENT
  !C D1RSIZ: LENGTH IN RADIAL   DIRECTION OF ELEMENT 
  !C D1RSIZ: LENGTH IN POLOIDAL DIRECTION OF ELEMENT 
  !C
  REAL(rkind),ALLOCATABLE,DIMENSION(:)::d1msiz,d1rsiz,d1psiz
  
  
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::&
       i1dbc2, i1mfc1 
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2mfc1

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2MFCS
  !C 
  !C------------------------------------------------------------------
  
  !C
  !C D2UG  : FRAME MOVING VELOCITY VECTOR
  !C D2JM1 : METRICS OF MAGNETIC SURFACE COORDINATES (will be changed) 
  !C D2RZM : RZ COORDINATES w   OVERLAP
  !C D2RZX : RZ COORDINATES w/o OVERLAP w INTERPOLATION 
  !C
  REAL(   rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2ug,d2rzm,d2rzx,d2jm1,d2mtrc

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2PROF
  !C
  !C------------------------------------------------------------------
  INTEGER(ikind)::&
       i0nm,i0nn,i0tm,i0tn
  REAL(   rkind)::&
       d0qc,d0qs,d0bc,d0rw
  REAL(   rkind),DIMENSION(1:i0spcsm)::&
       d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,d1pa,d1pz
  INTEGER(ikind),DIMENSION(0:i0lmaxm+1)::&
       i1mlvl
  INTEGER(ikind),DIMENSION(-1:i0lmaxm)::i1rdn2
  REAL(   rkind),DIMENSION(0:i0lmaxm)::&
       d1rec
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2STEP
  !C 
  !C------------------------------------------------------------------

  !C
  !C     Ax=b 
  !C D2XVEC[         1:I0VMAX,1:I0XMAX]: x
  !C D2BVEC[         1:I0VMAX,1:I0BMAX]: b
  !C D2XVEC_BEFOR[  1:I0VMAX,1:I0XMAX] : FOR PICARD ITERATION 
  !C D2XVEC_AFTER[  1:I0VMAX,1:I0XMAX] : FOR PICARD ITERATION 
  !C I1NLCT: NUMBER   OF ITERATION 
  !C D1RSDL: RESIDUAL OF ITERATION 
  REAL(   rkind),SAVE,ALLOCATABLE::&
       d2xvec(:,:),XvecIn(:,:),XvecOut(:,:)
  INTEGER(ikind),SAVE,ALLOCATABLE::&
       i1nlct(:)
  REAL(   rkind),SAVE,ALLOCATABLE::&
       d1rsdl(:)
  
  !------------------------------------------------------------------
  !
  !                         FOR T2CALV
  !
  !------------------------------------------------------------------ 
  
  REAL(   rkind),DIMENSION(:,:,:        ),ALLOCATABLE::d3ms
  REAL(   rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE::d4av
  REAL(   rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE::d6at
  REAL(   rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE::d5dt
  REAL(   rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE::d4gv
  REAL(   rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE::d6gt
  REAL(   rkind),DIMENSION(:,:,:        ),ALLOCATABLE::d3es
  REAL(   rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE::d5ev
  REAL(   rkind),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE::d7et
  REAL(   rkind),DIMENSION(:,:,:        ),ALLOCATABLE::d3ss

  REAL(   rkind),DIMENSION(:,:          ),ALLOCATABLE::d2ws

  REAL(   rkind),DIMENSION(:),ALLOCATABLE::&
       d1ee,d1mm,d1nn,d1ni,d1pp,d1pi,d1tt,d1ti,&
       d1ur,d1up,d1ut,d1ub,d1u2,&
       d1qr,d1qp,d1qt,d1qb,d1wb,d1wt,d1wp,d1vt,&
       d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,d1hex
  REAL(   rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2x,d2y,d2z,d2bcf
  REAL(   rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4

  REAL(   rkind)::&
       d0ct1_anom,d0ct2_anom
  REAL(   rkind),DIMENSION(:),ALLOCATABLE::&
       d1cx1_anom,d1cx2_anom
  
  INTEGER(ikind)::i0xa

  !------------------------------------------------------------------
  !
  !                         FOR T2CALV
  !
  !------------------------------------------------------------------ 

  REAL(   rkind),SAVE::&
       GRt,&
       G11xCo,&
       G12Co,&
       G22Co,&
       G33Co,&
       G11Ct,&
       G12Ct,&
       G22xCt,&
       G33Ct,&  
       UgrCt,&
       UgpCt,&
       R_rz,&
       R_mc

  REAL(   rkind),SAVE::&
       BtCo,& ! covariant     toroidal magnetic field [T*m]
       BtCt,& ! contravariant toroidal magnetic field [T/m]
       BtSq,& ! squared       toroidal megnetic field [T^2]
       BpCo,& ! covariant     poloidal magnetic field [T*m]
       BpCt,& ! contravariant poloidal magnetic field [T/m]
       BpSq,& ! squared       poloidal megnetic field [T^2]
       Bb,&   !                        megnetic field [T  ]
       BbSq,& ! squared                megnetic field [T^2]
       RR,RA
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       Pa,Pz,&
       Ee,&   ! electric charge                        [C]
       Mm,&   ! particle mass                          [kg]
       Nn,&   ! number density                         [ /m^3]
       Vv,&   ! thermal velocity                       [m/s     ]
       FrCt,& ! contravariant radial particle flux     [ /m^3*s ]
       Fb,&   ! parallel particle flux                 [ /m^2*s ]
       FtCo,& ! covariant toroildal particle flux      [ /m  *s ]
       FpCt,& ! contravariant poloidal particle flux   [ /m^3*s ]
       UrCt,& ! contravariant radial flow velocity     [   /s   ]
       Ub,&   ! parallel flow velocity                 [  m/s   ]
       UtCo,& ! covariant toroildal flow velocity      [m^2/s   ]
       UpCt,& ! contravariant flow velocity            [   /s   ]
       UuSq,& ! squared flow velocity                  [m^2/s^2 ]
       Pp,&   ! pressure                               [ J/m^3  ]
       QrCt,& ! contravariant radial total heat flux   [ J/m^3*s]
       Qb,&   ! parallel total heat flux               [ J/m^2*s]
       QtCo,& ! covariant toroildal total heat flux    [ J/m  *s]
       QpCt,& ! contravariant poloidal total heat flux [ J/m^3*s]
       WrCt,& ! contravariant radial energy velocity   [   /s   ]
       Wb,&   ! parallel energy flow velocity          [  m/s   ]
       WtCo,& ! covariant toroildal flow velocity      [m^2/s   ]
       WpCt,& ! contravariant flow velocity            [   /s   ]
       Tt     ! temperature                            [J       ]

  REAL(   rkind),SAVE,DIMENSION(:,:),ALLOCATABLE::&
       X,&    ! X_ab = Vth_b/Vth_a
       Y,&    ! Y_ab = m_a/m_b
       Z,&    ! Z_ab = T_a/T_b
       BaseNu ! base collision frequency [Hz]

  REAL(   rkind),SAVE,DIMENSION(:,:),ALLOCATABLE::&
       L11,L12,L21,L22,Lx1,Lx2,Lx3,Lx4
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       Hex ! Heat Exchange Rate [Hz]
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       Mu1,Mu2,Mu3,Mux1,Mux2,Mux3,Mux4
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       FtAnom1,&
       FtAnom2,&
       FtAnom3,&
       FtAnom4,&
       FtAnom5,&
       GtAnom1,&
       GtAnom2,&
       GtAnom3,&
       GtAnom4,&
       GtAnom5
       
  INTEGER(ikind),SAVE::&
       I_xa

  REAL(   rkind),SAVE,ALLOCATABLE::&
       KnownVar(:,:),&              ! d2ws
       !
       MassScaCoef(:,:,:        ),& ! d3ms
       AdveVecCoef(:,:,:,:      ),& ! d4av
       AdveTenCoef(:,:,:,:,:,:  ),& ! d6at
       DiffTenCoef(:,:,:,:,:    ),& ! d5dt
       GradVecCoef(:,:,:,:      ),& ! d4gv
       GradTenCoef(:,:,:,:,:,:  ),& ! d6gt
       ExciScaCoef(:,:,:        ),& ! d3es
       ExciVecCoef(:,:,:,:,:    ),& ! d5ev
       ExciTenCoef(:,:,:,:,:,:,:),& ! d7et
       SourScaCoef(:,:,:        )   ! d3ss 

  !C------------------------------------------------------------------
  !C
  !C                         FOR T2CONV
  !C
  !C------------------------------------------------------------------
  !C 
  !C d2xout[1   ,ixidi]: Poloidal Magnetic Field  [        T]
  !C d2xout[2   ,ixidi]: Toroidal Magnetic Field  [        T]
  !C d2xout[3   ,ixidi]: Toroidal Electric Field  [      V/m]
  !C d2xout[4   ,ixidi]: Poloidal Electric Field  [     mV/m]
  !C d2xout[5   ,ixidi]: Radial   Electric Field  [     mV/m]
  !C d2xout[8N-2,ixidi]: Particle Density         [10^20/m^3]
  !C d2xout[8N-1,ixidi]: Radial   Flow Velocity   [      m/s]
  !C d2xout[8N  ,ixidi]: Parallel Flow Velocity   [      m/s]
  !C d2xout[8N+1,ixidi]: Toroidal Rotation        [      m/s]
  !C d2xout[8N+2,ixidi]: Temperature              [      keV]
  !C d2xout[8N+3,ixidi]: Radial   Total Heat Flow [  keV*m/s]
  !C d2xout[8N+4,ixidi]: Parallel Total Heat Flow [  keV*m/s]
  !C d2xout[8N+5,ixidi]: Toroidal Total Heat Flow [  keV*m/s]
  !C
  !C
  !C
  !C
  
  REAL(   rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2xout

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2WRIT
  !C
  !C------------------------------------------------------------------
 
  CHARACTER(10)::c10rname
  INTEGER(ikind)::i0fnum
  REAL(   rkind),ALLOCATABLE,DIMENSION(:  )::&
       d1bp3,d1bt3,d1er3,d1ep3,d1et3
  REAL(   rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2n3, d2fr3,d2fb3,d2ft3,&
       d2p3, d2qr3,d2qb3,d2qt3

  !C
  !C TMP
  !C
  REAL(   rkind),ALLOCATABLE,DIMENSION(:)::&
       d1jm1,d1jm2,d1jm3,d1jm4,d1jm5
CONTAINS

  !C
  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2NGRA_ALLOCATE
    INTEGER(ikind),SAVE::i0lmax_save=0
    INTEGER(ikind)     :: i0err
    IF(i0lmax.NE.i0lmax_save) THEN
       
       CALL T2NGRA_DEALLOCATE
       
       DO 
          ALLOCATE(i1mmax( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1bmax( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1emax( 0:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1rdn1(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1pdn1(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1pdn2(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1msiz( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1rsiz( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1psiz( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          
          i0lmax_save=i0lmax
          
          WRITE(6,'(A)') '-- T2NGRD_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')'XX T2NGRD_ALLOCATE: ALLOCATION ERROR: ECODE=',i0err
       STOP
       
    END IF
    
    RETURN
  END SUBROUTINE T2NGRA_ALLOCATE
  
  SUBROUTINE T2NGRA_DEALLOCATE

    IF(ALLOCATED(i1mmax)) DEALLOCATE(i1mmax)
    IF(ALLOCATED(i1bmax)) DEALLOCATE(i1bmax)
    IF(ALLOCATED(i1emax)) DEALLOCATE(i1emax)
    IF(ALLOCATED(i1rdn1)) DEALLOCATE(i1rdn1)
    IF(ALLOCATED(i1pdn1)) DEALLOCATE(i1pdn1)
    IF(ALLOCATED(i1pdn2)) DEALLOCATE(i1pdn2)
    IF(ALLOCATED(d1msiz)) DEALLOCATE(d1msiz)
    IF(ALLOCATED(d1rsiz)) DEALLOCATE(d1rsiz)
    IF(ALLOCATED(d1psiz)) DEALLOCATE(d1psiz)

    RETURN

  END SUBROUTINE T2NGRA_DEALLOCATE

  !-------------------------------------------------------------------
  !
  !
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2COMM_ALLOCATE
    
    CALL T2COMM_ADHOC
    
    CALL T2COMM_ALLOCATE_INTG
    CALL T2COMM_ALLOCATE_VGRA
    
    CALL T2COMM_ALLOCATE_REMAINS
    
    RETURN

  END SUBROUTINE T2COMM_ALLOCATE
  
  SUBROUTINE T2COMM_DEALLOCATE

    CALL T2COMM_DEALLOCATE_INTG
    CALL T2COMM_DEALLOCATE_VGRA

    CALL T2COMM_DEALLOCATE_REMAINS

    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE

  SUBROUTINE T2COMM_ALLOCATE_REMAINS
    
    INTEGER(ikind),SAVE::&
         i0vmax_save=0,i0nmax_save=0,i0smax_save=0,i0dmax_save=0,&
         i0emax_save=0,i0hmax_save=0,i0qmax_save=0,&
         i0mmax_save=0,i0xmax_save=0,i0bmax_save=0,i0amax_save=0,&
         i0nrmx_save=0,i0ermx_save=0,i0ecmx_save=0
    
    INTEGER(ikind),SAVE::&
         nv0dmax_save=0,nt0dmax_save=0, &
         nv2dmax_save=0,nt2dmax_save=0
    
    INTEGER(ikind)     :: i0err

    IF(  (i0smax .NE.i0smax_save ).OR.&
         (i0nmax .NE.i0nmax_save ).OR.&
         (i0dmax .NE.i0dmax_save ).OR.&
         (i0qmax .NE.i0qmax_save ).OR.&
         (i0mmax .NE.i0mmax_save ).OR.&
         (i0amax .NE.i0amax_save ).OR.&
         (i0xmax .NE.i0xmax_save ).OR.&
         (i0bmax .NE.i0bmax_save ).OR.&
         (i0emax .NE.i0emax_save ).OR.&
         (i0nrmx .NE.i0nrmx_save ).OR.&
         (i0ermx .NE.i0ermx_save ).OR.&
         (i0ecmx .NE.i0ecmx_save ).OR.&
         (i0hmax .NE.i0hmax_save ).OR.&
         (i0vmax .NE.i0vmax_save ).OR.&
         (nv0dmax.NE.nv0dmax_save).OR.&
         (nt0dmax.NE.nt0dmax_save).OR.&
         (nv2dmax.NE.nv2dmax_save).OR.&
         (nt2dmax.NE.nt2dmax_save))THEN
       
       CALL T2COMM_DEALLOCATE_REMAINS

       DO
                    
          !C
          !C T2NGRA
          !C

          ALLOCATE(i1nidr(    1:i0nrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1nidc(    1:i0amax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1eidr(    1:i0ermx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1eidc(    1:i0ecmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mlel(    1:i0emax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2crt( 1:3,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1dbc2(    1:i0bmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mfc1(    1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mc1d(    1:i0rmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1mc1d(    1:i0rmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i3enr( 1:i0nmax,1:4,1:i0emax),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(d2mfc1(1:2,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          IF(i0hmax.NE.0)&
               ALLOCATE(i2hbc(1:2,1:i0hmax),STAT=i0err);&
               IF(i0err.NE.0)EXIT
                    
          !C
          !C T2MFCS
          !C
          
          ALLOCATE(d2ug( 1:2,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jm1(1:9,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzm(1:2,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzx(1:2,1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ! TMP
          ALLOCATE(d2mtrc(1:2,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT

          !C
          !C T2CALV
          !C

          ALLOCATE(d3ms(1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4av(1:i0dmax,1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d6at(1:i0dmax,1:i0dmax,1:i0wmax,1:i0vmax,1:i0vmax,&
               1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d5dt(1:i0dmax,1:i0dmax,1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4gv(1:i0dmax,1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d6gt(1:i0dmax,1:i0dmax,1:i0wmax,1:i0vmax,1:i0vmax,&
               1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3es(1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d5ev(1:i0dmax,1:i0wmax,1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d7et(1:i0dmax,1:i0dmax,1:i0wmax,1:i0wmax,&
               1:i0vmax,1:i0vmax,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3ss(1:i0vmax,1:i0vmax,1:i0mmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2ws(1:i0wmax,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d1ee(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1mm(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nn(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ur(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1up(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ut(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ub(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1u2(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1pp(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qr(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qp(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qt(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qb(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1wb(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1wt(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1wp(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1tt(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1vt(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ni(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1pi(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ti(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT

          ALLOCATE(d1nvcc1(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nvcc2(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nvcc3(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nvcc4(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1hex(  1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf1(1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf2(1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf3(1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf4(1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2x(    1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2y(    1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2z(    1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2bcf(  1:i0smax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d1cx1_anom(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1cx2_anom(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT

          !C FOR T2WRIT
          ALLOCATE(d1bp3(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1bt3(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1er3(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ep3(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1et3(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2n3( 1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2fr3(1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2fb3(1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2ft3(1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2p3( 1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2qr3(1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2qb3(1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2qt3(1:i0xmax,1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT

          !C
          !C T2STEP
          !C
          ALLOCATE(i1nlct(     1:10000),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1rsdl(     1:10000),STAT=i0err);IF(i0err.NE.0)EXIT
          !C
          !C T2CONV
          !C
          ALLOCATE(d2xout(1:i0vmax,1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          
          !C TMP

          ALLOCATE(d1jm1(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm2(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm3(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm4(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm5(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(vv0d(nv0dmax,nt0dmax),STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(vv2d(nv2dmax,nt2dmax),STAT=i0err); IF(i0err.NE.0) EXIT

          i0smax_save = i0smax
          i0nmax_save = i0nmax
          i0dmax_save = i0dmax
          i0qmax_save = i0qmax
          i0mmax_save = i0mmax
          i0bmax_save = i0bmax
          i0xmax_save = i0xmax
          i0amax_save = i0amax
          i0emax_save = i0emax
          i0nrmx_save = i0nrmx
          i0ermx_save = i0ermx
          i0ecmx_save = i0ecmx
          i0hmax_save = i0hmax 
          i0vmax_save = i0vmax

          nv0dmax_save = nv0dmax
          nt0dmax_save = nt0dmax
          nv2dmax_save = nv2dmax
          nt2dmax_save = nt2dmax
          
          WRITE(6,'(A)') '-- T2COMM_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)') 'XX T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',i0err
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2COMM_ALLOCATE_REMAINS
  
  SUBROUTINE T2COMM_DEALLOCATE_REMAINS
        
    !C
    !C T2NGRA
    !C

    IF(ALLOCATED(i1nidr)) DEALLOCATE(i1nidr)
    IF(ALLOCATED(i1nidc)) DEALLOCATE(i1nidc)
    IF(ALLOCATED(i1eidr)) DEALLOCATE(i1eidr)
    IF(ALLOCATED(i1eidc)) DEALLOCATE(i1eidc)
    IF(ALLOCATED(i1mlel)) DEALLOCATE(i1mlel)
    IF(ALLOCATED(i2crt )) DEALLOCATE(i2crt )
    IF(ALLOCATED(i1dbc2)) DEALLOCATE(i1dbc2)
    IF(ALLOCATED(i1mfc1)) DEALLOCATE(i1mfc1)
    IF(ALLOCATED(d2mfc1)) DEALLOCATE(d2mfc1)
    IF(ALLOCATED(i3enr )) DEALLOCATE(i3enr )
    IF(ALLOCATED(i2hbc )) DEALLOCATE(i2hbc )
    IF(ALLOCATED(i1mc1d)) DEALLOCATE(i1mc1d)
    IF(ALLOCATED(d1mc1d)) DEALLOCATE(d1mc1d)    

    
    !C
    !C T2MFCS
    !C

    IF(ALLOCATED(d2ug )) DEALLOCATE(d2ug )
    IF(ALLOCATED(d2jm1)) DEALLOCATE(d2jm1)
    IF(ALLOCATED(d2rzm)) DEALLOCATE(d2rzm)
    IF(ALLOCATED(d2rzx)) DEALLOCATE(d2rzx)
    IF(ALLOCATED(d2mtrc)) DEALLOCATE(d2mtrc)
    !C
    !C FOR T2CALV
    !C

    IF(ALLOCATED(d3ms)) DEALLOCATE(d3ms)
    IF(ALLOCATED(d4av)) DEALLOCATE(d4av)
    IF(ALLOCATED(d6at)) DEALLOCATE(d6at)
    IF(ALLOCATED(d5dt)) DEALLOCATE(d5dt)
    IF(ALLOCATED(d4gv)) DEALLOCATE(d4gv)
    IF(ALLOCATED(d6gt)) DEALLOCATE(d6gt)
    IF(ALLOCATED(d3es)) DEALLOCATE(d3es)
    IF(ALLOCATED(d5ev)) DEALLOCATE(d5ev)
    IF(ALLOCATED(d7et)) DEALLOCATE(d7et)
    IF(ALLOCATED(d3ss)) DEALLOCATE(d3ss)
    IF(ALLOCATED(d2ws)) DEALLOCATE(d2ws)
    
    IF(ALLOCATED(d1ee)) DEALLOCATE(d1ee)
    IF(ALLOCATED(d1mm)) DEALLOCATE(d1mm)
    IF(ALLOCATED(d1nn)) DEALLOCATE(d1nn)
    IF(ALLOCATED(d1ur)) DEALLOCATE(d1ur)
    IF(ALLOCATED(d1up)) DEALLOCATE(d1up)
    IF(ALLOCATED(d1ut)) DEALLOCATE(d1ut)
    IF(ALLOCATED(d1ub)) DEALLOCATE(d1ub)
    IF(ALLOCATED(d1u2)) DEALLOCATE(d1u2)
    IF(ALLOCATED(d1pp)) DEALLOCATE(d1pp)
    IF(ALLOCATED(d1qr)) DEALLOCATE(d1qr)
    IF(ALLOCATED(d1qp)) DEALLOCATE(d1qp)
    IF(ALLOCATED(d1qt)) DEALLOCATE(d1qt)
    IF(ALLOCATED(d1qb)) DEALLOCATE(d1qb)
    IF(ALLOCATED(d1wb)) DEALLOCATE(d1wb)
    IF(ALLOCATED(d1wt)) DEALLOCATE(d1wt)
    IF(ALLOCATED(d1wp)) DEALLOCATE(d1wp)
    IF(ALLOCATED(d1tt)) DEALLOCATE(d1tt)
    IF(ALLOCATED(d1vt)) DEALLOCATE(d1vt)
    IF(ALLOCATED(d1ni)) DEALLOCATE(d1ni)
    IF(ALLOCATED(d1pi)) DEALLOCATE(d1pi)
    IF(ALLOCATED(d1ti)) DEALLOCATE(d1ti)

    IF(ALLOCATED(d1hex)) DEALLOCATE(d1hex)
    IF(ALLOCATED(d2bcf)) DEALLOCATE(d2bcf)

    IF(ALLOCATED(d1nvcc1)) DEALLOCATE(d1nvcc1)
    IF(ALLOCATED(d1nvcc2)) DEALLOCATE(d1nvcc2)
    IF(ALLOCATED(d1nvcc3)) DEALLOCATE(d1nvcc3)
    IF(ALLOCATED(d1nvcc4)) DEALLOCATE(d1nvcc4)

    IF(ALLOCATED(d2nfcf1)) DEALLOCATE(d2nfcf1)
    IF(ALLOCATED(d2nfcf2)) DEALLOCATE(d2nfcf2)
    IF(ALLOCATED(d2nfcf3)) DEALLOCATE(d2nfcf3)
    IF(ALLOCATED(d2nfcf4)) DEALLOCATE(d2nfcf4)

    IF(ALLOCATED(d2x)) DEALLOCATE(d2x)
    IF(ALLOCATED(d2y)) DEALLOCATE(d2y)
    IF(ALLOCATED(d2z)) DEALLOCATE(d2z)

    IF(ALLOCATED(d1cx1_anom)) DEALLOCATE(d1cx1_anom)
    IF(ALLOCATED(d1cx2_anom)) DEALLOCATE(d1cx2_anom)
    
    !C
    !C T2STEP
    !C
    
    IF(ALLOCATED(i1nlct))       DEALLOCATE(i1nlct) 
    IF(ALLOCATED(d1rsdl))       DEALLOCATE(d1rsdl)
    
    IF(ALLOCATED(d2xout))       DEALLOCATE(d2xout)
    
    !C
    !C  FOR T2WRIT
    !C
    
    IF(ALLOCATED(d1bp3)) DEALLOCATE(d1bp3)
    IF(ALLOCATED(d1bt3)) DEALLOCATE(d1bt3)
    IF(ALLOCATED(d1er3)) DEALLOCATE(d1er3)
    IF(ALLOCATED(d1ep3)) DEALLOCATE(d1ep3)
    IF(ALLOCATED(d1et3)) DEALLOCATE(d1et3)
    IF(ALLOCATED(d2n3 )) DEALLOCATE(d2n3 )
    IF(ALLOCATED(d2fr3)) DEALLOCATE(d2fr3)
    IF(ALLOCATED(d2fb3)) DEALLOCATE(d2fb3)
    IF(ALLOCATED(d2ft3)) DEALLOCATE(d2ft3)
    IF(ALLOCATED(d2p3 )) DEALLOCATE(d2p3 )
    IF(ALLOCATED(d2qr3)) DEALLOCATE(d2qr3)
    IF(ALLOCATED(d2qb3)) DEALLOCATE(d2qb3)
    IF(ALLOCATED(d2qt3)) DEALLOCATE(d2qt3)    


    IF(ALLOCATED(d1jm1  )) DEALLOCATE(d1jm1  )
    IF(ALLOCATED(d1jm2  )) DEALLOCATE(d1jm2  )
    IF(ALLOCATED(d1jm3  )) DEALLOCATE(d1jm3  )
    IF(ALLOCATED(d1jm4  )) DEALLOCATE(d1jm4  )
    IF(ALLOCATED(d1jm5  )) DEALLOCATE(d1jm5  )

    !C
    !C added by A.FUKUYAMA
    !C

    IF(ALLOCATED(vv0d   )) DEALLOCATE(vv0d   )
    IF(ALLOCATED(vv2d   )) DEALLOCATE(vv2d   )

    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE_REMAINS

  SUBROUTINE T2COMM_ADHOC
    
    NNMAX = i0nmax
    NDMAX = i0dmax
    NQMAX = i0qmax
    NSMAX = i0smax
    NVMAX = i0vmax
    NRMAX = i0rmax
    NKMAX = i0wmax
    NMMAX = i0mmax
    NXMAX = i0xmax
    NBMAX = i0bmax
    NEMAX = i0emax
    NHMAX = i0hmax
    NLMAX = i0lmax
    NAMAX = i0amax
    NNRMX = i0nrmx

    RETURN
    
  END SUBROUTINE T2COMM_ADHOC

  !-------------------------------------------------------------------
  !
  !       ALLOCATOR OF GLOBAL VARIABLES  FOR T2INTG
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2COMM_ALLOCATE_INTG
    
    INTEGER(ikind),SAVE::&
         NNMAX_save=0,NDMAX_save=0,NQMAX_save=0
    
    INTEGER(ikind):: ierr
    
    IF(  (NNMAX .NE.NNMAX_save ).OR.&
         (NDMAX .NE.NDMAX_save ).OR.&
         (NQMAX .NE.NQMAX_save ))THEN
       
       CALL T2COMM_DEALLOCATE_INTG
       
       DO
          ! for PG-FEM
          ALLOCATE(MassScaIntgPG(1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(AdveVecIntgPG(1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(AdveTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(DiffTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(GradVecIntgPG(1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(GradTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(ExciScaIntgPG(1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(ExciVecIntgPG(1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(ExciTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(SourScaIntgPG(1:NNMAX,1:NNMAX),&
               STAT=ierr); IF(ierr.NE.0) EXIT
          
          ! for SUPG-FEM

          !ALLOCATE(massScaIntgSUPG(1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=ierr); IF(i0err.NE.0) EXIT
          !ALLOCATE(adveVecIntgSUPG(1:NDMAX,1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=ierr); IF(i0err.NE.0) EXIT
          !ALLOCATE(sourScaIntgSUPG(1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=ierr); IF(i0err.NE.0) EXIT
          
          NNMAX_save = NNMAX
          NDMAX_save = NDMAX
          NQMAX_save = NQMAX
          
          WRITE(6,'(A)') '-- T2INTG_PUBLIC_ALLOCATE: SUCCESSED'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)') 'XX T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2COMM_ALLOCATE_INTG
  
  !-------------------------------------------------------------------
  !
  !       DEALLOCATOR OF GLOBAL VARIABLES  FOR T2INTG
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------
  SUBROUTINE T2COMM_DEALLOCATE_INTG
    
    ! for PG-FEM

    IF(ALLOCATED(massScaIntgPG)) DEALLOCATE(massScaIntgPG)
    IF(ALLOCATED(adveVecIntgPG)) DEALLOCATE(adveVecIntgPG)
    IF(ALLOCATED(adveTenIntgPG)) DEALLOCATE(adveTenIntgPG)
    IF(ALLOCATED(diffTenIntgPG)) DEALLOCATE(diffTenIntgPG)
    IF(ALLOCATED(gradVecIntgPG)) DEALLOCATE(gradVecIntgPG)
    IF(ALLOCATED(gradTenIntgPG)) DEALLOCATE(gradTenIntgPG)
    IF(ALLOCATED(exciScaIntgPG)) DEALLOCATE(exciScaIntgPG)
    IF(ALLOCATED(exciVecIntgPG)) DEALLOCATE(exciVecIntgPG)
    IF(ALLOCATED(exciTenIntgPG)) DEALLOCATE(exciTenIntgPG)
    IF(ALLOCATED(sourScaIntgPG)) DEALLOCATE(sourScaIntgPG)
    
    !C for SUPG-FEM
    
    !IF(ALLOCATED(massScaIntgSUPG)) DEALLOCATE(massScaIntgSUPG)
    !IF(ALLOCATED(adveVecIntgSUPG)) DEALLOCATE(adveVecIntgSUPG)
    !IF(ALLOCATED(sourScaIntgSUPG)) DEALLOCATE(sourScaIntgSUPG)
   
    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE_INTG

  !-------------------------------------------------------------------
  !
  !       ALLOCATOR OF GLOBAL VARIABLES  FOR T2VGRA
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2COMM_ALLOCATE_VGRA
    
    INTEGER(ikind),SAVE::&
         NVMAX_save=0,NKMAX_save=0
    
    INTEGER(ikind):: ierr
    
    IF(  (NVMAX.NE.NVMAX_save).OR.&
         (NKMAX.NE.NKMAX_save))THEN
       
       CALL T2COMM_DEALLOCATE_VGRA
       
       DO
          
          ALLOCATE(HaveMassScaCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveAdveVecCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveAdveTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveDiffTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradVecCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciScaCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciVecCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveSourScaCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
    
          ALLOCATE(HaveAdveTenKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradTenKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciVecKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciTenKval(1:NKMAX,1:NKMAX,1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
    
          ALLOCATE(HaveMat(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockEqs(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockAxi(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockWal(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          
          NVMAX_save = NVMAX
          NKMAX_save = NKMAX 
                 
          WRITE(6,'(A)') '-- T2VGRA_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')&
            '--T2VGRA_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2COMM_ALLOCATE_VGRA

  !-------------------------------------------------------------------
  !
  !       DEALLOCATOR OF GLOBAL VARIABLES  FOR T2VGRA
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2COMM_DEALLOCATE_VGRA
    
    IF(ALLOCATED(HaveMassScaCoef))  DEALLOCATE(HaveMassScaCoef)
    IF(ALLOCATED(HaveAdveVecCoef))  DEALLOCATE(HaveAdveVecCoef)
    IF(ALLOCATED(HaveAdveTenCoef))  DEALLOCATE(HaveAdveTenCoef)
    IF(ALLOCATED(HaveDiffTenCoef))  DEALLOCATE(HaveDiffTenCoef)
    IF(ALLOCATED(HaveGradVecCoef))  DEALLOCATE(HaveGradVecCoef)
    IF(ALLOCATED(HaveGradTenCoef))  DEALLOCATE(HaveGradTenCoef)
    IF(ALLOCATED(HaveExciScaCoef))  DEALLOCATE(HaveExciScaCoef)
    IF(ALLOCATED(HaveExciVecCoef))  DEALLOCATE(HaveExciVecCoef)
    IF(ALLOCATED(HaveExciTenCoef))  DEALLOCATE(HaveExciTenCoef)
    IF(ALLOCATED(HaveSourScaCoef))  DEALLOCATE(HaveSourScaCoef)
    
    IF(ALLOCATED(HaveAdveTenKval))  DEALLOCATE(HaveAdveTenKval)
    IF(ALLOCATED(HaveGradTenKval))  DEALLOCATE(HaveGradTenKval)
    IF(ALLOCATED(HaveExciVecKval))  DEALLOCATE(HaveExciVecKval)
    IF(ALLOCATED(HaveExciTenKval))  DEALLOCATE(HaveExciTenKval)
    
    IF(ALLOCATED(HaveMat))          DEALLOCATE(HaveMat)
    IF(ALLOCATED(LockEqs))          DEALLOCATE(LockEqs)
    IF(ALLOCATED(LockAxi))          DEALLOCATE(LockAxi)
    IF(ALLOCATED(LockWal))          DEALLOCATE(LockWal)
    
    RETURN
    
  END SUBROUTINE T2COMM_DEALLOCATE_VGRA

END MODULE T2COMM
