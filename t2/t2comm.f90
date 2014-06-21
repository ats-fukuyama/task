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
  
  !---------------------------------------------------------
  !
  !
  !       GLOBAL PARAMETERS (NEW)
  ! 
  !
  !---------------------------------------------------------

  !
  ! INPUT PARAMETERS
  !
  INTEGER(ikind),SAVE::&
       NNMAX,& ! NUMBER OF NODES PAR ELEMENT
       NQMAX,& ! NUMBER OF INTEGRAL POINTS PAR DIRECTION
       NDMAX,& ! NUMBER OF DIMENSIONS
       NSMAX,& ! NUMBER OF SPECIES
       NPMIN,& ! POLOIDAL DIVISION NUMBER IN LOWEST LEVEL MESH 
       NLMAX   ! NUMBER OF SUB DOMAINS 

  LOGICAL,SAVE::&
       UsePotentialDescription,&
       UseNormalization,&
       UseSUPGFEM,&
       UseCoefficientCheck,&
       UseAnomalousTransportFT,&
       UseAnomalousTransportGT,&
       !
       SolveElectron,SolveIons,&
       SolveBp,SolveBt,SolveEt,SolveEp,SolveEr,&
       SolveNn,SolveFr,SolveFb,SolveFt,SolveFp,&
       SolvePp,SolveQr,SolveQb,SolveQt,SolveQp,&
       !
       LockBpOnAxis,LockBtOnAxis,LockEtOnAxis,LockEpOnAxis,LockErOnAxis,&
       LockNnOnAxis,LockFrOnAxis,LockFbOnAxis,LockFtOnAxis,LockFpOnAxis,&
       LockPpOnAxis,LockQrOnAxis,LockQbOnAxis,LockQtOnAxis,LockQpOnAxis,&
       !
       LockBpOnWall,LockBtOnWall,LockEtOnWall,LockEpOnWall,LockErOnWall,&
       LockNnOnWall,LockFrOnWall,LockFbOnWall,LockFtOnWall,LockFpOnWall,&
       LockPpOnWall,LockQrOnWall,LockQbOnWall,LockQtOnWall,LockQpOnWall,&
       TestMS,TestAV,TestAT,TestDT,TestGV,TestGT,TestES,TestEV,TestET,TestSS,&
       TestLEQ,TestLAX,TestLWL
  
  INTEGER(ikind),SAVE::&
       CoordinateSwitch,TestCase,EqSet
  CHARACTER(10)::c10rname
  INTEGER(ikind)::i0fnum
  
  INTEGER(ikind)::&
       i0nm,i0nn,i0tm,i0tn
  REAL(   rkind)::&
       d0qc,d0qs,d0bc,d0rw,RR,RA
  REAL(   rkind),DIMENSION(1:i0spcsm)::&
       d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,Pa,Pz
  INTEGER(ikind),DIMENSION(0:i0lmaxm+1)::&
       i1mlvl
  INTEGER(ikind),DIMENSION(-1:i0lmaxm)::&
       i1rdn2  ! NUMBER OF PARTITION IN RADIAL DIREACTON
               !      IN EACH SUBDOMAIN 
  REAL(   rkind),DIMENSION(0:i0lmaxm)::&
       d1rec
  
  !
  ! global parameters calculated in TASK/T2
  !
  INTEGER(ikind),SAVE::&
       NVMAX,& ! NUMBER OF DEPENDENT VARIABLES
       NVFMX,& ! NUMBER OF DEPENDENT VARIABLES (1D)
       NRMAX,& ! NUMBER OF NODES IN RADIAL DIRECTION (FOR 1D)
       NKMAX,& ! NUMBER OF WORKING VARIABLES FOR DIFFERENTIAL
       NEMAX,& ! NUMBER OF ELEMENTS IN DOMAIN
       NBMAX,& ! NUMBER OF NODES W/O OVERLAP AND INTERPOLATION IN DOMAIN
       NXMAX,& ! NUMBER OF NODES W/O OVERLAP W INTERPOLATION IN DOMAIN 
       NMMAX,& ! NUMBER OF NODES W   OVERLAP IN DOMAIN
       NHMAX,& ! NUMBER OF INTERPOLATION NODES IN DOMAIN
       NAMAX,& ! NUMBER OF NON-ZERO COMPONENTS OF MATRIX (CRS-METHOD)
       NNRMX,& ! ARRAY SIZE OF NodeRowCRS (CRS-METHOD) 
       NERMX,& ! ARRAY SIZE OF I1EIDR     (CRS-METHOD) 
       NECMX,& ! ARRAY SIZE OF I1EIDC     (CRS-METHOD) 
       NAVMX,& ! i0avmax = i0amax*i0vmax*i0vmax
       NBVMX   ! i0bvmax = i0bmax*i0vmax
  INTEGER(ikind)::&
       i0wstp  ! INDICATOR FOR RESULT OUTPUT TIMING    
     
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
       QpNF,&
       EqBpNF,&
       EqBtNF,&
       EqEtNF,&
       EqEpNF,&
       EqErNF,&
       EqNnNF,&
       EqFrNF,&
       EqFbNF,&
       EqFtNF,&
       EqFpNF,&
       EqPpNF,&
       EqQrNF,&
       EqQbNF,&
       EqQtNF,&
       EqQpNF
  
  INTEGER(ikind),SAVE::&
       StartEqs,EndEqs,StartAxi,EndAxi,StartWal,EndWal
  !-------------------------------------------------------------------
  !
  !       DEFINITION OF GLOBAL VARIABLES 
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------
  REAL(   rkind),SAVE,DIMENSION(:,:),ALLOCATABLE::&
       Xvec  ! [1:NVMAX,1:NXMAX]: 
  
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
  

  !------------------------------------------------------------------
  !
  !       DEFINITION OF GLOBAL VARIABLES FOR T2NGRA
  !
  !                 LAST UPDATE 2014-05-27 H.Seto
  !
  !------------------------------------------------------------------ 
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::& 
       i1emax,&     ! NUMBER OF ELEMENTS IN EACH SUBDOMAIN
       i1rdn1,&     ! NUMBER OF PARTITION IN RADIAL   DIREACTON
                    !      IN EACH SUBDOMAIN 
       i1pdn1,&     ! NUMBER OF PARTITION IN POLOIDAL DIREACTON
                    !      IN EACH SUBDOMAIN
       i1pdn2       ! NUMBER OF PARTITION IN POLOIDAL DIREACTON
                    !      IN EACH SUBDOMAIN WITHOUT OVERLAP
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::& 
       NodeRowCRS,& ! CUMULTATIVE NUMBER OF NODE-NODE CONNECTIVITY+1 
                    !     UP TO EACH ROW IN MATRIX A
       NodeColCRS,& ! NODE NUMBER 
       NodeDiaCRS   ! Position of DIagonal component in MATRIX A 
  
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::&
       i1eidr,&     ! CUMULTATIVE NUMBER OF NODE-ELEMENT CONNECTIVITY + 1
                    !      UP TO EACH NODE 
       i1eidc       ! ELEMENT NUMBER 
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:,:,:)::&
       ElementNodeGraph ! INFORMATION OF ELEMENT-NODE CONNECTIVITY 
                        !     1: FOR COEFFICIENT CALCULATION
                        !     2: FOR 2D 
                        !     3: FOR 1D-2D 
                        !     4: FOR 1D
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:,:  )::&
       i2crt,&         ! INFORMATION NODE PROJECTION MAP
                       !      1: FOR COEFFICIENT CALCULATION
                       !      2: FOR 2D
                       !      3: FOR 1D
       HangedNodeTable ! INFORMATION OF BOUND NODE AND BINDING NODES
    INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::&
       i1mc1d ! NODE   NUMBERS   FOR DATA STOCK OF 1D VALUES
  REAL(   rkind),ALLOCATABLE,DIMENSION(:)::&
       d1mc1d ! RADIAL POSITIONS FOR DATA STOCK OF 1D VALUES
  INTEGER(ikind),ALLOCATABLE,DIMENSION(:)::&
       i1mfc1 
  
  !-------------------------------------------------------------------
  !
  !       DEFINITION OF GLOBAL VARIABLES FOR T2VGRA
  !
  !                 LAST UPDATE 2014-05-27 H.Seto
  ! 
  !-------------------------------------------------------------------
  
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
  !------------------------------------------------------------------
  !
  !                          FOR T2PREP
  ! 
  !------------------------------------------------------------------
  REAL(   rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2ug,&  ! FRAME MOVING VELOCITY VECTOR  with OVERLAP
       d2rzm,& ! RZ COORDINATES                with OVERLAP
       Metric  ! Metrics                       with OVERLAP 
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:)::&
       GlobalCrd
  !------------------------------------------------------------------
  !
  !                          FOR T2STEP
  ! 
  !------------------------------------------------------------------
  !    Amat*Xvec = Bvec
  REAL(   rkind),SAVE,DIMENSION(:,:),ALLOCATABLE::&
       XvecIn,& ! [1:NVMAX,1:NXMAX]: FOR PICARD ITERATION
       XvecOut  ! [1:NVMAX,1:NXMAX]: FOR PICARD ITERATION 
  INTEGER(ikind),SAVE,ALLOCATABLE::&
       i1nlct(:)     ! [1:10000]: num. of NL-its. in each time step
  REAL(   rkind),SAVE,ALLOCATABLE::&
       d1rsdl(:)     ! [1:10000]: residuals       in each time step
    
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
       Bb,  & !                        megnetic field [T  ]
       BbSq   ! squared                megnetic field [T^2]
   
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       Ee,&   ! electric charge                        [C]
       Mm,&   ! particle mass                          [kg]
       Nn,&   ! number density                         [ /m^3]
       Vv,&   ! thermal velocity                       [m/s     ]
       FrCt,& ! contravariant radial particle flux     [ /m^3*s ]
       Fb,&   ! parallel particle flux                 [ /m^2*s ]
       FtCo,& ! covariant toroildal particle flux      [ /m^3*s ]
       FtCt,& ! contravariant toroildal particle flux  [ /m  *s ]
       FpCt,& ! contravariant poloidal particle flux   [ /m^3*s ]
       FpCo,& ! covariant poloidal particle flux       [ /m  *s ]
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
       X,&      ! X_ab = Vth_b/Vth_a
       Y,&      ! Y_ab = m_a/m_b
       Z,&      ! Z_ab = T_a/T_b
       BaseNu,& ! base collision frequency [Hz]
       Nu,&     !      collision frequency [Hz]
       Nuh      ! heat exchange  frequency [Hz]

  REAL(   rkind),SAVE,DIMENSION(:,:),ALLOCATABLE::&
       L11,L12,L21,L22
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       Hex ! Heat Exchange Rate [Hz]
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       Mu1,Mu2,Mu3
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

  ! additional intermidiate variables 
  REAL(   rkind),SAVE::&
       & BNCXb1,BNCXb2,BNCXb3,BNCXb4,BNCXb5,BNCXb6,&
       & BNCXt1,BNCXt2,BNCXt3,&
       & BNCPp1,BNCPp2,BNCPp3,BNCPp4,BNCPp5,BNCPp6,&
       & BNCPp7,BNCPp8,BNCPp9,&
       & BNCQb1,BNCQb2,BNCQb3,BNCQb4,&
       & BNCQt1,BNCQt2,BNCQt3,BNCQt4
       
  REAL(   rkind),SAVE,DIMENSION(:),ALLOCATABLE::&
       & CNCV01,CNCV02,CNCV03,CNCV04,CNCV05,&
       & CNCV06,CNCV07,CNCV08,CNCV09,CNCV10,&
       & CNCV11,CNCV12,CNCV13

  REAL(   rkind),SAVE,DIMENSION(:,:),ALLOCATABLE::&
       & CNCF01,CNCF02,CNCF03,CNCF04

  ! 
  !  T2COEF
  !
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

  !------------------------------------------------------------------
  !
  !                         FOR T2CONV
  !
  !------------------------------------------------------------------
  ! 
  ! d2xout[1   ,ixidi]: Poloidal Magnetic Field  [        T]
  ! d2xout[2   ,ixidi]: Toroidal Magnetic Field  [        T]
  ! d2xout[3   ,ixidi]: Toroidal Electric Field  [      V/m]
  ! d2xout[4   ,ixidi]: Poloidal Electric Field  [     mV/m]
  ! d2xout[5   ,ixidi]: Radial   Electric Field  [     mV/m]
  ! d2xout[10*N-4,ixidi]: Particle Density       [1.D20/m^3]
  ! d2xout[10*N-3,ixidi]: Radial   Flow Velocity [      m/s]
  ! d2xout[10*N-2,ixidi]: Parallel Flow Velocity [      m/s]
  ! d2xout[10*N-1,ixidi]: Toroidal Rotation      [      m/s]
  ! d2xout[10*N  ,ixidi]: Poloidal Rotation      [      m/s]
  ! d2xout[10*N+1,ixidi]: Temperature            [      keV]
  ! d2xout[10*N+2,ixidi]: Radial   Heat Flow     [  keV*m/s]
  ! d2xout[10*N+3,ixidi]: Parallel Heat Flow     [  keV*m/s]
  ! d2xout[10*N+4,ixidi]: Toroidal Heat Flow     [  keV*m/s]
  ! d2xout[10*N+5,ixidi]: Poloidal Heat Flow     [  keV*m/s]
  
  REAL(   rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2xout
  
CONTAINS
  
  SUBROUTINE T2NGRA_ALLOCATE
    INTEGER(ikind),SAVE::NLMAX_save=0
    INTEGER(ikind)     :: ierr
    IF(NLMAX.NE.NLMAX_save) THEN
       
       CALL T2NGRA_DEALLOCATE
       
       DO 
          !ALLOCATE(i1mmax( 1:NLMAX  ),STAT=ierr);IF(ierr.NE.0)EXIT
          !ALLOCATE(i1bmax( 1:NLMAX  ),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1emax( 0:NLMAX  ),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1rdn1(-1:NLMAX  ),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1pdn1(-1:NLMAX  ),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1pdn2(-1:NLMAX  ),STAT=ierr);IF(ierr.NE.0)EXIT
          
          NLMAX_save=NLMAX
          
          WRITE(6,'(A)') '-- T2NGRD_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')'XX T2NGRD_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2NGRA_ALLOCATE
  
  SUBROUTINE T2NGRA_DEALLOCATE
    
    !IF(ALLOCATED(i1mmax)) DEALLOCATE(i1mmax)
    !IF(ALLOCATED(i1bmax)) DEALLOCATE(i1bmax)
    IF(ALLOCATED(i1emax)) DEALLOCATE(i1emax)
    IF(ALLOCATED(i1rdn1)) DEALLOCATE(i1rdn1)
    IF(ALLOCATED(i1pdn1)) DEALLOCATE(i1pdn1)
    IF(ALLOCATED(i1pdn2)) DEALLOCATE(i1pdn2)

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
    
    CALL T2COMM_ALLOCATE_NGRA
    CALL T2COMM_ALLOCATE_INTG
    CALL T2COMM_ALLOCATE_VGRA
    CALL T2COMM_ALLOCATE_CALV
    CALL T2COMM_ALLOCATE_COEF
    CALL T2COMM_ALLOCATE_REMAINS
    
    RETURN

  END SUBROUTINE T2COMM_ALLOCATE
  
  SUBROUTINE T2COMM_DEALLOCATE
    
    CALL T2COMM_DEALLOCATE_NGRA
    CALL T2COMM_DEALLOCATE_INTG
    CALL T2COMM_DEALLOCATE_VGRA
    CALL T2COMM_DEALLOCATE_CALV
    CALL T2COMM_DEALLOCATE_COEF
    CALL T2COMM_DEALLOCATE_REMAINS
    
    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE
  
  SUBROUTINE T2COMM_ALLOCATE_REMAINS
    
    INTEGER(ikind),SAVE::&
         NVMAX_save=0,NXMAX_save=0,NMMAX_save=0,&
         nv0dmax_save=0,nt0dmax_save=0, &
         nv2dmax_save=0,nt2dmax_save=0
    
    INTEGER(ikind)     :: ierr
    
    IF(  (NMMAX  .NE.NMMAX_save ) .OR.&
         (NVMAX  .NE.NVMAX_save ) .OR.&
         (NXMAX  .NE.NXMAX_save ) .OR.&
         (nv0dmax.NE.nv0dmax_save).OR.&
         (nt0dmax.NE.nt0dmax_save).OR.&
         (nv2dmax.NE.nv2dmax_save).OR.&
         (nt2dmax.NE.nt2dmax_save))THEN
       
       CALL T2COMM_DEALLOCATE_REMAINS
       
       DO         
          ALLOCATE(d2ug( 1:2,1:NMMAX) ,STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(d2rzm(1:2,1:NMMAX) ,STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Metric(1:9,1:NMMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ! T2STEP
          ALLOCATE(Xvec(   1:NVMAX,1:NXMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(XvecIn( 1:NVMAX,1:NXMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(XvecOut(1:NVMAX,1:NXMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1nlct(     1:10000),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(d1rsdl(     1:10000),STAT=ierr);IF(ierr.NE.0)EXIT
          ! T2CONV
          ALLOCATE(d2xout(1:NVMAX,1:NXMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ! added by A.Fukuyama
          ALLOCATE(vv0d(nv0dmax,nt0dmax),STAT=ierr); IF(ierr.NE.0) EXIT
          ALLOCATE(vv2d(nv2dmax,nt2dmax),STAT=ierr); IF(ierr.NE.0) EXIT

          NMMAX_save = NMMAX
          NXMAX_save = NXMAX
          NVMAX_save = NVMAX
         
          nv0dmax_save = nv0dmax
          nt0dmax_save = nt0dmax
          nv2dmax_save = nv2dmax
          nt2dmax_save = nt2dmax
          
          WRITE(6,'(A)') '-- T2COMM_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)') 'XX T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2COMM_ALLOCATE_REMAINS
  
  SUBROUTINE T2COMM_DEALLOCATE_REMAINS     
    ! T2PREP
    IF(ALLOCATED(d2ug  ))  DEALLOCATE(d2ug  )
    IF(ALLOCATED(d2rzm ))  DEALLOCATE(d2rzm )    
    IF(ALLOCATED(Metric))  DEALLOCATE(Metric)
    ! T2STEP
    IF(ALLOCATED(i1nlct))  DEALLOCATE(i1nlct) 
    IF(ALLOCATED(d1rsdl))  DEALLOCATE(d1rsdl)
    ! T2CONV
    IF(ALLOCATED(Xvec   )) DEALLOCATE(Xvec   )
    IF(ALLOCATED(XvecIn )) DEALLOCATE(XvecIn )
    IF(ALLOCATED(XvecOut)) DEALLOCATE(XvecOut)
    IF(ALLOCATED(d2xout )) DEALLOCATE(d2xout )
    ! added by A.FUKUYAMA
    IF(ALLOCATED(vv0d   )) DEALLOCATE(vv0d   )
    IF(ALLOCATED(vv2d   )) DEALLOCATE(vv2d   )

    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE_REMAINS
  
  !-------------------------------------------------------------------
  !
  !       ALLOCATOR OF GLOBAL VARIABLES  FOR T2NGRA
  !
  !                                     LAST UPDATE 2014-06-05
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2COMM_ALLOCATE_NGRA
    
    INTEGER(ikind),SAVE::&
         NNMAX_save=0,NDMAX_save=0,NEMAX_save=0,NHMAX_save=0,&
         NMMAX_save=0,NAMAX_save=0,NNRMX_save=0,NERMX_save=0,&
         NECMX_save=0,NBMAX_save=0
    
    INTEGER(ikind):: ierr
    
    IF(  (NNMAX.NE.NNMAX_save).OR.&
         (NDMAX.NE.NDMAX_save).OR.&
         (NMMAX.NE.NMMAX_save).OR.&
         (NBMAX.NE.NBMAX_save).OR.&
         (NAMAX.NE.NAMAX_save).OR.&
         (NEMAX.NE.NEMAX_save).OR.&
         (NNRMX.NE.NNRMX_save).OR.&
         (NERMX.NE.NERMX_save).OR.&
         (NECMX.NE.NECMX_save).OR.&
         (NHMAX.NE.NHMAX_save))THEN
       
       CALL T2COMM_DEALLOCATE_NGRA

       DO
          ALLOCATE(NodeRowCRS(1:NNRMX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(NodeColCRS(1:NAMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(NodeDiaCRS(1:NBMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ElementNodeGraph(1:NNMAX,1:4,1:NEMAX),&
               &                       STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GlobalCrd(1:NDMAX,1:NMMAX),&
               &                       STAT=ierr);IF(ierr.NE.0)EXIT
          IF(NHMAX.NE.0) ALLOCATE(HangedNodeTable(1:2,1:NHMAX),&
               &                       STAT=ierr);IF(ierr.NE.0)EXIT

          ALLOCATE(i1eidr(    1:NERMX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1eidc(    1:NECMX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i2crt( 1:3,1:NMMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1mfc1(    1:NMMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(i1mc1d(    1:NRMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(d1mc1d(    1:NRMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          
          NNMAX_save = NNMAX
          NDMAX_save = NDMAX
          NMMAX_save = NMMAX
          NAMAX_save = NAMAX
          NEMAX_save = NEMAX
          NNRMX_save = NNRMX
          NERMX_save = NERMX
          NECMX_save = NECMX
          NHMAX_save = NHMAX
          
          WRITE(6,'(A)') '-- T2COMM_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)') 'XX T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2COMM_ALLOCATE_NGRA
  
  !-------------------------------------------------------------------
  !
  !      DEALLOCATOR OF GLOBAL VARIABLES  FOR T2NGRA
  !
  !                                     LAST UPDATE 2014-06-05
  ! 
  !-------------------------------------------------------------------
  SUBROUTINE T2COMM_DEALLOCATE_NGRA
    
    IF(ALLOCATED(NodeRowCRS)) DEALLOCATE(NodeRowCRS)
    IF(ALLOCATED(NodeColCRS)) DEALLOCATE(NodeColCRS)
    IF(ALLOCATED(NodeDiaCRS)) DEALLOCATE(NodeDiaCRS)
    IF(ALLOCATED(GlobalCrd)) DEALLOCATE(GlobalCrd)
    IF(ALLOCATED(ElementNodeGraph)) DEALLOCATE(ElementNodeGraph)
    IF(ALLOCATED(HangedNodeTable )) DEALLOCATE(HangedNodeTable )
    IF(ALLOCATED(i1eidr)) DEALLOCATE(i1eidr)
    IF(ALLOCATED(i1eidc)) DEALLOCATE(i1eidc)
    IF(ALLOCATED(i2crt )) DEALLOCATE(i2crt )
    IF(ALLOCATED(i1mfc1)) DEALLOCATE(i1mfc1)
    IF(ALLOCATED(i1mc1d)) DEALLOCATE(i1mc1d)
    IF(ALLOCATED(d1mc1d)) DEALLOCATE(d1mc1d)    
    
    RETURN
    
  END SUBROUTINE T2COMM_DEALLOCATE_NGRA
  
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
          !     STAT=ierr); IF(ierr.NE.0) EXIT
          !ALLOCATE(adveVecIntgSUPG(1:NDMAX,1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=ierr); IF(ierr.NE.0) EXIT
          !ALLOCATE(sourScaIntgSUPG(1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=ierr); IF(ierr.NE.0) EXIT
          
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

    IF(ALLOCATED(MassScaIntgPG)) DEALLOCATE(MassScaIntgPG)
    IF(ALLOCATED(AdveVecIntgPG)) DEALLOCATE(AdveVecIntgPG)
    IF(ALLOCATED(AdveTenIntgPG)) DEALLOCATE(AdveTenIntgPG)
    IF(ALLOCATED(DiffTenIntgPG)) DEALLOCATE(DiffTenIntgPG)
    IF(ALLOCATED(GradVecIntgPG)) DEALLOCATE(GradVecIntgPG)
    IF(ALLOCATED(GradTenIntgPG)) DEALLOCATE(GradTenIntgPG)
    IF(ALLOCATED(ExciScaIntgPG)) DEALLOCATE(ExciScaIntgPG)
    IF(ALLOCATED(ExciVecIntgPG)) DEALLOCATE(ExciVecIntgPG)
    IF(ALLOCATED(ExciTenIntgPG)) DEALLOCATE(ExciTenIntgPG)
    IF(ALLOCATED(SourScaIntgPG)) DEALLOCATE(SourScaIntgPG)
    
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
    If(ALLOCATED(HaveExciVecCoef))  DEALLOCATE(HaveExciVecCoef)
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

  !-------------------------------------------------------------------
  !
  !       ALLOCATOR OF GLOBAL VARIABLES  FOR T2CALV
  !
  !                                     LAST UPDATE 2014-06-05
  ! 
  !------------------------------------------------------------------  
  SUBROUTINE T2COMM_ALLOCATE_CALV
    INTEGER(ikind),SAVE::&
         NSMAX_save=0
    INTEGER(ikind):: ierr
    
    IF(NSMAX.NE.NSMAX_save)THEN
       
       CALL T2COMM_DEALLOCATE_COEF
       
       DO
          ALLOCATE(Ee(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Mm(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Nn(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Vv(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FrCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Fb(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FtCo(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FtCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FpCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FpCo(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(UrCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Ub(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(UtCo(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(UpCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(UuSq(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Pp(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(QrCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Qb(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(QtCo(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(QpCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(WrCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Wb(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(WtCo(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(WpCt(1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Tt(  1:NSMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          
          ALLOCATE(X(      1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Y(      1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Z(      1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(BaseNu( 1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Nu(     1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Nuh(    1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT

          ALLOCATE(L11(    1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(L12(    1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(L21(    1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(L22(    1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT

          ALLOCATE(Hex(    1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Mu1(    1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Mu2(    1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(Mu3(    1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT

          ALLOCATE(FtAnom1(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FtAnom2(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FtAnom3(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FtAnom4(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(FtAnom5(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GtAnom1(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GtAnom2(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GtAnom3(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GtAnom4(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GtAnom5(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          
          ALLOCATE(CNCV01(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV02(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV03(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV04(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV05(1:NSMAX),&
               &                STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV06(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV07(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV08(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV09(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV10(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV11(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV12(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCV13(1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT

          ALLOCATE(CNCF01(1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCF02(1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCF03(1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(CNCF04(1:NSMAX,1:NSMAX),&
               &                 STAT=ierr);IF(ierr.NE.0)EXIT
          
          NSMAX_save = NSMAX
          
          WRITE(6,'(A)') '-- T2VGRA_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')&
            '--T2VGRA_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    ENDIF
        
  END SUBROUTINE T2COMM_ALLOCATE_CALV

  !-------------------------------------------------------------------
  !
  !       DEALLOCATOR OF GLOBAL VARIABLES  FOR T2CALV
  !
  !                                     LAST UPDATE 2014-06-05
  ! 
  !------------------------------------------------------------------  
  SUBROUTINE T2COMM_DEALLOCATE_CALV

    IF(ALLOCATED(Ee  ))    DEALLOCATE(Ee  )
    IF(ALLOCATED(Mm  ))    DEALLOCATE(Mm  )
    IF(ALLOCATED(Nn  ))    DEALLOCATE(Nn  )
    IF(ALLOCATED(Vv  ))    DEALLOCATE(Vv  )
    IF(ALLOCATED(FrCt))    DEALLOCATE(FrCt)
    IF(ALLOCATED(Fb  ))    DEALLOCATE(Fb  )
    IF(ALLOCATED(FtCo))    DEALLOCATE(FtCo)
    IF(ALLOCATED(FtCt))    DEALLOCATE(FtCt)
    IF(ALLOCATED(FpCo))    DEALLOCATE(FpCo)
    IF(ALLOCATED(FpCt))    DEALLOCATE(FpCt)
    IF(ALLOCATED(UrCt))    DEALLOCATE(UrCt)
    IF(ALLOCATED(Ub  ))    DEALLOCATE(Ub  )
    IF(ALLOCATED(UtCo))    DEALLOCATE(UtCo)
    IF(ALLOCATED(UpCt))    DEALLOCATE(UpCt)
    IF(ALLOCATED(UuSq))    DEALLOCATE(UuSq)
    IF(ALLOCATED(Pp  ))    DEALLOCATE(Pp  )
    IF(ALLOCATED(QrCt))    DEALLOCATE(QrCt)
    IF(ALLOCATED(Qb  ))    DEALLOCATE(Qb  )
    IF(ALLOCATED(QtCo))    DEALLOCATE(QtCo)
    IF(ALLOCATED(QpCt))    DEALLOCATE(QpCt)
    IF(ALLOCATED(WrCt))    DEALLOCATE(WrCt)
    IF(ALLOCATED(Wb  ))    DEALLOCATE(Wb  )
    IF(ALLOCATED(WtCo))    DEALLOCATE(WtCo)
    IF(ALLOCATED(WpCt))    DEALLOCATE(WpCt)
    IF(ALLOCATED(Tt  ))    DEALLOCATE(Tt  )
    
    IF(ALLOCATED(X     ))  DEALLOCATE(X     )
    IF(ALLOCATED(Y     ))  DEALLOCATE(Y     )
    IF(ALLOCATED(Z     ))  DEALLOCATE(Z     )
    IF(ALLOCATED(BaseNu))  DEALLOCATE(BaseNu)
    IF(ALLOCATED(Nu    ))  DEALLOCATE(Nu    )
    IF(ALLOCATED(Nuh   ))  DEALLOCATE(Nuh   )
          
    IF(ALLOCATED(L11))     DEALLOCATE(L11)
    IF(ALLOCATED(L12))     DEALLOCATE(L12)
    IF(ALLOCATED(L21))     DEALLOCATE(L21)
    IF(ALLOCATED(L22))     DEALLOCATE(L22)
    
    IF(ALLOCATED(Hex ))    DEALLOCATE(Hex )
    IF(ALLOCATED(Mu1 ))    DEALLOCATE(Mu1 )
    IF(ALLOCATED(Mu2 ))    DEALLOCATE(Mu2 )
    IF(ALLOCATED(Mu3 ))    DEALLOCATE(Mu3 )

    IF(ALLOCATED(FtAnom1)) DEALLOCATE(FtAnom1)
    IF(ALLOCATED(FtAnom2)) DEALLOCATE(FtAnom2)
    IF(ALLOCATED(FtAnom3)) DEALLOCATE(FtAnom3)
    IF(ALLOCATED(FtAnom4)) DEALLOCATE(FtAnom4)
    IF(ALLOCATED(FtAnom5)) DEALLOCATE(FtAnom5)
    IF(ALLOCATED(GtAnom1)) DEALLOCATE(GtAnom1)
    IF(ALLOCATED(GtAnom2)) DEALLOCATE(GtAnom2)
    IF(ALLOCATED(GtAnom3)) DEALLOCATE(GtAnom3)
    IF(ALLOCATED(GtAnom4)) DEALLOCATE(GtAnom4)
    IF(ALLOCATED(GtAnom5)) DEALLOCATE(GtAnom5)
    
    IF(ALLOCATED(CNCV01)) DEALLOCATE(CNCV01)
    IF(ALLOCATED(CNCV02)) DEALLOCATE(CNCV02)
    IF(ALLOCATED(CNCV03)) DEALLOCATE(CNCV03)
    IF(ALLOCATED(CNCV04)) DEALLOCATE(CNCV04)
    IF(ALLOCATED(CNCV05)) DEALLOCATE(CNCV05)
    IF(ALLOCATED(CNCV06)) DEALLOCATE(CNCV06)
    IF(ALLOCATED(CNCV07)) DEALLOCATE(CNCV07)
    IF(ALLOCATED(CNCV08)) DEALLOCATE(CNCV08)
    IF(ALLOCATED(CNCV09)) DEALLOCATE(CNCV09)
    IF(ALLOCATED(CNCV10)) DEALLOCATE(CNCV10)
    IF(ALLOCATED(CNCV11)) DEALLOCATE(CNCV11)
    IF(ALLOCATED(CNCV12)) DEALLOCATE(CNCV12)
    IF(ALLOCATED(CNCV13)) DEALLOCATE(CNCV13)

    IF(ALLOCATED(CNCF01)) DEALLOCATE(CNCF01)
    IF(ALLOCATED(CNCF02)) DEALLOCATE(CNCF02)
    IF(ALLOCATED(CNCF03)) DEALLOCATE(CNCF03)
    IF(ALLOCATED(CNCF04)) DEALLOCATE(CNCF04)

    RETURN
  
  END SUBROUTINE T2COMM_DEALLOCATE_CALV
  
  !-------------------------------------------------------------------
  !
  !       DEALLOCATOR OF GLOBAL VARIABLES  FOR T2COEF
  !
  !                                     LAST UPDATE 2014-06-05
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2COMM_ALLOCATE_COEF
    
    INTEGER(ikind),SAVE::&
         NVMAX_save=0,NDMAX_save=0,NMMAX_save=0,NKMAX_save=0
    INTEGER(ikind):: ierr
    
    IF(  (NDMAX.NE.NDMAX_save ).OR.&
         (NMMAX.NE.NMMAX_save ).OR.&
         (NKMAX.NE.NKMAX_save ).OR.&
         (NVMAX.NE.NVMAX_save ))THEN
       
       CALL T2COMM_DEALLOCATE_COEF
       
       DO
          
          ALLOCATE(KnownVar(   1:NKMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(MassScaCoef(1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(AdveVecCoef(1:NDMAX,1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(AdveTenCoef(1:NDMAX,1:NDMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(DiffTenCoef(1:NDMAX,1:NDMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GradVecCoef(1:NDMAX,1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GradTenCoef(1:NDMAX,1:NDMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ExciScaCoef(1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ExciVecCoef(1:NDMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ExciTenCoef(1:NDMAX,1:NDMAX,1:NKMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(SourScaCoef(1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          
          NKMAX_save = NKMAX
          NVMAX_save = NVMAX
          NDMAX_save = NDMAX
          NMMAX_save = NMMAX
          
          WRITE(6,'(A)') '-- T2COMM_ALLOCATE_COEF: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')&
            '--T2COMM_ALLOCATE_COEF: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2COMM_ALLOCATE_COEF

  !-------------------------------------------------------------------
  !
  !       DEALLOCATOR OF GLOBAL VARIABLES  FOR T2COEF
  !
  !                                     LAST UPDATE 2014-06-05
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2COMM_DEALLOCATE_COEF

    IF(ALLOCATED(KnownVar))         DEALLOCATE(KnownVar)
    !
    IF(ALLOCATED(MassScaCoef))      DEALLOCATE(MassScaCoef)
    IF(ALLOCATED(AdveVecCoef))      DEALLOCATE(AdveVecCoef)
    IF(ALLOCATED(AdveTenCoef))      DEALLOCATE(AdveTenCoef)
    IF(ALLOCATED(DiffTenCoef))      DEALLOCATE(DiffTenCoef)
    IF(ALLOCATED(GradVecCoef))      DEALLOCATE(GradVecCoef)
    IF(ALLOCATED(GradTenCoef))      DEALLOCATE(GradTenCoef)
    IF(ALLOCATED(ExciScaCoef))      DEALLOCATE(ExciScaCoef)
    IF(ALLOCATED(ExciVecCoef))      DEALLOCATE(ExciVecCoef)
    IF(ALLOCATED(ExciTenCoef))      DEALLOCATE(ExciTenCoef)
    IF(ALLOCATED(SourScaCoef))      DEALLOCATE(SourScaCoef)
    
    RETURN
    
  END SUBROUTINE T2COMM_DEALLOCATE_COEF

END MODULE T2COMM
