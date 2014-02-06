MODULE T2COMM
  
  USE T2CNST, ONLY: i0rkind, i0ikind, i0lmaxm, i0spcsm, twopi
  
  IMPLICIT NONE

  !---- input parameters added by AF
  INTEGER(i0ikind):: &
       ntmax,ntstep,nconvmax, &
       nt0dmax,nt0dstep, &
       nt2dmax,nt2dstep, &
       idfile,idprint,idplot,idmode,idebug
  REAL(i0rkind):: &
       dt,time_init,eps_conv

  !---- global parameters added by AF
  INTEGER(i0ikind):: nrhomax ! number of variables in rho
  INTEGER(i0ikind):: nchimax ! number of variables in chi
  INTEGER(i0ikind):: nequmax ! number of equations to be solved
  INTEGER(i0ikind):: nv0dmax ! number of global variables to be saved
  INTEGER(i0ikind):: nv2dmax ! number of profile varibales to be saved
  REAL(   i0rkind):: &
       time_t2                            ! global time
  REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: &
       vv0d,vv2d                          ! storage for 0D and 2D data
  
  !C---------------------------------------------------------
  !C
  !C
  !C  GLOBAL PARAMETERS
  !C 
  !C
  !C---------------------------------------------------------
  INTEGER(i0ikind)::i0dbg

  INTEGER(i0ikind) ::&
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

  INTEGER(i0ikind)::&
       i0mfcs, & !C INDICATOR FOR COORDINATE SYSTEM (1: torus coordinate)
       i0sflg, & !C INDICATOR FOR SUPG METHOD (0: w/o SUPG, 1: w SUPG)
       i0wstp    !C INDICATOR FOR RESULT OUTPUT TIMING
  INTEGER(i0ikind)::&
       i0bvmax,& !C VECTOR SIZE OF b FOR MTXP (Ax=b)
       i0avmax   !C NUMBER OF NONZERO COMPONENT OF A FOR MTXP (Ax=b)

  REAL(   i0rkind)::&
       d0mfcst, & !C NORMALIZATION CONSTANT FOR \psi'
       d0btcst, & !C NORMALIZATION CONSTANT FOR I
       d0etcst, & !C NORMALIZATION CONSTANT FOR E_{\zeta}
       d0epcst, & !C NORMALIZATION CONSTANT FOR E_{\chi }
       d0ercst, & !C NORMALIZATION CONSTANT FOR E_{\rho }
       d0nncst, & !C NORMALIZATION CONSTANT FOR n_{a}
       d0frcst, & !C NORMALIZATION CONSTANT FOR n_{a}u_{a}^{\rho}
       d0fbcst, & !C NORMALIZATION CONSTANT FOR n_{a}u_{a\para}
       d0ftcst, & !C NORMALIZATION CONSTANT FOR n_{a}u_{a\zeta}
       d0ppcst, & !C NORMALIZATION CONSTANT FOR p_{a}
       d0qrcst, & !C NORMALIZATION CONSTANT FOR Q_{a}^{\rho }
       d0qbcst, & !C NORMALIZATION CONSTANT FOR Q_{a\para}
       d0qtcst, & !C NORMALIZATION CONSTANT FOR Q_{a\zeta}
       d0rmjr,  & !C MAJOR RADIUS (R_{0} [m])
       d0rmnr,  & !C MINOR RADIUS (a     [m])
       d0iar,   & !C INVERSE ASPECT RATIO (a/R_{0})
       d0eps      !C CONVERGENCE CRITERION FOR PICARD ITERATION

  !C--------------------------------------------------------
  !C
  !C                 FOR T2INTG: MODIFIED 2014-01-29
  !C 
  !C--------------------------------------------------------


  !C
  !C WORKING ARRAY FOR GAUSSIAN INTEGRATION
  !C

  !C D1ABSC: ABSCISSAS FOR GAUSS INTEGRATION
  !C D1WFCT: WHEIGHT FACTOR (1D)
  !C D2WFCT: WHEIGHT FACTOR (2D)
  !C D4IFNC: INTERPOLATION FUNCTIONS 
  
  REAL(i0rkind),DIMENSION(:      ),ALLOCATABLE:: d1absc    
  REAL(i0rkind),DIMENSION(:      ),ALLOCATABLE:: d1wfct
  REAL(i0rkind),DIMENSION(:,:    ),ALLOCATABLE:: d2wfct
  REAL(i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: d4ifnc
  
  !C
  !C INTEGRATION ARRAYS BY GAUSSIAN INTEGRATION
  !C

  !C
  !C FOR W/O SUPG 
  !C
  !C D3IMSN: INTEGRATION ARRAY FOR MASS       SCALAR SUBMATRIX
  !C D4IAVN: INTEGRATION ARRAY FOR ADVECTION  VECTOR SUBMATRIX
  !C D6IATN: INTEGRATION ARRAY FOR ADVECTION  TENSOR SUBMATRIX
  !C D5IDTN: INTEGRATION ARRAY FOR DIFFUSION  TENSOR SUBMATRIX
  !C D4IGVN: INTEGRATION ARRAY FOR GRADIENT   VECTOR SUBMATRIX
  !C D6IGTN: INTEGRATION ARRAY FOR GRADIENT   TENSOR SUBMATRIX
  !C D3IESN: INTEGRATION ARRAY FOR EXCITATION SCALAR SUBMATRIX
  !C D5IEVN: INTEGRATION ARRAY FOR EXCITATION VECTOR SUBMATRIX
  !C D7IETN: INTEGRATION ARRAY FOR EXCITATION TENSOR SUBMATRIX
  !C D2ISSN: INTEGRATION ARRAY FOR SOURCE     SCALAR SUBMATRIX
  
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE:: d3imsn
  REAL(   i0rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE:: d4iavn
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE:: d6iatn
  REAL(   i0rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE:: d5idtn
  REAL(   i0rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE:: d4igvn
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE:: d6igtn
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE:: d3iesn
  REAL(   i0rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE:: d5ievn
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE:: d7ietn
  REAL(   i0rkind),DIMENSION(:,:          ),ALLOCATABLE:: d2issn
  
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

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::& 
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
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::& 
       i1nidr, i1nidc

  !C
  !C FOR COMPRESSED ROW STORAGE METHOD FOR NODE-ELEMENT RERATIONS
  !C
  !C I1EDIR: CUMULTATIVE NUMBER OF NODE-ELEMENT CONNECTIVITY + 1
  !C         UP TO EACH NODE 
  !C I1EDIC: ELEMENT NUMBER 
  !C
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::&
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
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:,:)::i3enr 
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:  )::i2crt,i2hbc
  
  !C
  !C FOR 1D PROBLEM
  !C
  !C I1MC1D: NODE   NUMBERS   FOR DATA STOCK OF 1D VALUES
  !C D1MC1D: RADIAL POSITIONS FOR DATA STOCK OF 1D VALUES
  !C 
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:):: i1mc1d 
  REAL(   i0rkind),ALLOCATABLE,DIMENSION(:):: d1mc1d



  !C
  !C D1MSIZ: SIZE OF ELEMENT
  !C D1RSIZ: LENGTH IN RADIAL   DIRECTION OF ELEMENT 
  !C D1RSIZ: LENGTH IN POLOIDAL DIRECTION OF ELEMENT 
  !C
  REAL(i0rkind),ALLOCATABLE,DIMENSION(:)::d1msiz,d1rsiz,d1psiz
  
  
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::&
       i1dbc2, i1mfc1 
  REAL(i0rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2mfc1

  !C------------------------------------------------------------------
  !C  
  !C                    FOR T2VGRA: MODIFIED 2014-01-28
  !C
  !C------------------------------------------------------------------
  
  !C
  !C DEPENDENT VARIABLE GRAPHS: I2**VT[I0VIDI,I0VIDJ]
  !C       IVIDI >>
  !C                1     : \psi'
  !C                2     : I
  !C                3     : E_{\zeta}
  !C                4     : E_{\chi }
  !C                5     : E_{\rho }
  !C                8*N-2 : n_{a}
  !C                8*N-1 : n_{a}u_{a}^{\rho}
  !C                8*N   : n_{a}u_{a\para}
  !C                8*N+1 : n_{a}u_{a\zeta}
  !C                8*N+2 : p_{a}
  !C                8*N+3 : Q_{a}^{\rho}
  !C                8*N+4 : Q_{a\para}
  !C                8*N+5 : Q_{a\zeta}
  !C
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:)::&
       i2msvt,i2avvt,i2atvt,i2dtvt,i2gvvt,i2gtvt,&
       i2esvt,i2evvt,i2etvt,i2ssvt,i2vvvt

  !C
  !C DEPENDENT VARIABLE & DIFFERENATIAL GRAPHS: 
  !C 
  !C I3**WT[I0WIDI,       I0VIDI,I0VIDJ]
  !C I4**WT[I0WIDI,I0WIDJ,I0VIDI,I0VIDJ]
  !C       IWIDI >>
  !C                1     : B 
  !C                2     : R 
  !C                2*N+1 : u_{a\para}
  !C                2*N+2 : Q_{a\para}/p_{a}
  !C
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:,:  )::i3atwt,i3gtwt,i3evwt
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:,:,:)::i4etwt

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
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE:: d2ug,d2rzm,d2rzx,d2jm1

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2PROF
  !C
  !C------------------------------------------------------------------
  INTEGER(i0ikind)::&
       i0m0,i0n0
  REAL(   i0rkind)::&
       d0qc,d0qs,d0bc,d0rw
  REAL(   i0rkind),DIMENSION(1:i0spcsm)::&
       d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,d1pa,d1pz
  INTEGER(i0ikind),DIMENSION(0:i0lmaxm+1)::&
       i1mlvl
  INTEGER(i0ikind),DIMENSION(-1:i0lmaxm)::i1rdn2
  REAL(   i0rkind),DIMENSION(0:i0lmaxm)::&
       d1rec
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2STEP
  !C 
  !C------------------------------------------------------------------

  !C
  !C     Ax=b 
  !C D3AMAT[1:I0VMAX,1:I0VMAX,1:I0AMAX]: A (CRS-METHOD)
  !C D2XVEC[         1:I0VMAX,1:I0XMAX]: x
  !C D2BVEC[         1:I0VMAX,1:I0BMAX]: b
  !C D2XVEC_BEFOR[  1:I0VMAX,1:I0XMAX] : FOR PICARD ITERATION 
  !C D2XVEC_AFTER[  1:I0VMAX,1:I0XMAX] : FOR PICARD ITERATION 
  REAL(   i0rkind),DIMENSION(:,:,:),ALLOCATABLE::d3amat
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2xvec,d2bvec,d2xvec_befor,d2xvec_after
  !C 
  !C I1NLCT: NUMBER   OF ITERATION 
  !C D1RSDL: RESIDUAL OF ITERATION 
  !C
  INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1nlct
  REAL(   i0rkind),DIMENSION(:),ALLOCATABLE::d1rsdl
  
  !C------------------------------------------------------------------
  !C
  !C                         FOR T2CALV
  !C
  !C------------------------------------------------------------------ 
  
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE::d3ms
  REAL(   i0rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE::d4av
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE::d6at
  REAL(   i0rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE::d5dt
  REAL(   i0rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE::d4gv
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE::d6gt
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE::d3es
  REAL(   i0rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE::d5ev
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE::d7et
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE::d3ss

  REAL(   i0rkind),DIMENSION(:,:          ),ALLOCATABLE::d2ws

  REAL(   i0rkind),DIMENSION(:),ALLOCATABLE::&
       d1ee,d1mm,d1nn,d1ni,d1pp,d1pi,d1tt,d1ti,&
       d1ur,d1up,d1ut,d1ub,d1u2,&
       d1qr,d1qp,d1qt,d1qb,d1vb,d1vt,&
       d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,d1hex
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2x,d2y,d2z,d2bcf
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4
  INTEGER(i0ikind)::i0xa
  !C------------------------------------------------------------------
  !C
  !C                         FOR T2EXEC
  !C
  !C------------------------------------------------------------------
  
  INTEGER(i0ikind),DIMENSION(:,:),ALLOCATABLE::&
       i2enr0
  REAL(   i0rkind)::&
       d0jdmp
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2jmpm
  REAL(   i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE::&
       d4smat
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2svec
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2kwns
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2wrks

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
  
  REAL(   i0rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2xout

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2WRIT
  !C
  !C------------------------------------------------------------------
 
  CHARACTER(10)::c10rname
  INTEGER(i0ikind)::i0fnum
  REAL(   i0rkind),ALLOCATABLE,DIMENSION(:  )::&
       d1bp3,d1bt3,d1er3,d1ep3,d1et3
  REAL(   i0rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2n3, d2fr3,d2fb3,d2ft3,&
       d2p3, d2qr3,d2qb3,d2qt3

  !C
  !C TMP
  !C
  REAL(   i0rkind),ALLOCATABLE,DIMENSION(:)::&
       d1jm1,d1jm2,d1jm3,d1jm4,d1jm5
CONTAINS

  !C
  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2NGRA_ALLOCATE
    INTEGER(i0ikind),SAVE::i0lmax_save=0
    INTEGER(i0ikind)     :: i0err
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

  SUBROUTINE T2COMM_ALLOCATE
    
    INTEGER(i0ikind),SAVE::&
         i0vmax_save=0,i0nmax_save=0,i0smax_save=0,i0dmax_save=0,&
         i0emax_save=0,i0hmax_save=0,i0qmax_save=0,&
         i0mmax_save=0,i0xmax_save=0,i0bmax_save=0,i0amax_save=0,&
         i0nrmx_save=0,i0ermx_save=0,i0ecmx_save=0
    
    INTEGER(i0ikind),SAVE::&
         nv0dmax_save=0,nt0dmax_save=0, &
         nv2dmax_save=0,nt2dmax_save=0
    
    INTEGER(i0ikind)     :: i0err

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
       
       CALL T2COMM_DEALLOCATE

       DO
          
          !C
          !C T2INTG
          !C
         
          ALLOCATE(d3imsn(1:i0nmax,1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4iavn(1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d6iatn(1:i0dmax,1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax,&
               1:i0nmax),STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d5idtn(1:i0dmax,1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4igvn(1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d6igtn(1:i0dmax,1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax,&
               1:i0nmax),STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d3iesn(1:i0nmax,1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d5ievn(1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d7ietn(1:i0dmax,1:i0dmax,1:i0nmax,1:i0nmax,1:i0nmax,&
               1:i0nmax,1:i0nmax),STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d2issn(1:i0nmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          
          ALLOCATE(d1wfct(1:i0qmax),STAT=i0err);IF(i0err.NE.0) EXIT
          ALLOCATE(d1absc(1:i0qmax),STAT=i0err);IF(i0err.NE.0) EXIT
          ALLOCATE(d2wfct(1:i0qmax,1:i0qmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4ifnc(1:i0qmax,1:i0qmax,0:i0dmax,1:i0nmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT

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
          !C T2VGRA
          !C
          
          ALLOCATE(i2msvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2avvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2atvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2dtvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2gvvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2gtvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2esvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2evvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2etvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2ssvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2vvvt(1:i0vmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(i3atwt(1:i0wmax,         1:i0vmax,1:i0vmax), &
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i3gtwt(1:i0wmax,         1:i0vmax,1:i0vmax), &
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i3evwt(1:i0wmax,         1:i0vmax,1:i0vmax), &
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i4etwt(1:i0wmax,1:i0wmax,1:i0vmax,1:i0vmax), &
               STAT=i0err);IF(i0err.NE.0)EXIT
          
          !C
          !C T2MFCS
          !C
          
          ALLOCATE(d2ug( 1:2,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jm1(1:5,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzm(1:2,1:i0mmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzx(1:2,1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT

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
          ALLOCATE(d1vb(1:i0smax),STAT=i0err);IF(i0err.NE.0)EXIT
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
          ALLOCATE(d3amat(1:i0vmax,1:i0vmax,1:i0amax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2bvec(      1:i0vmax,1:i0bmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2xvec(      1:i0vmax,1:i0xmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2xvec_befor(1:i0vmax,1:i0xmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2xvec_after(1:i0vmax,1:i0xmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT

          !C
          !C T2CONV
          !C
          ALLOCATE(d2xout(1:i0vmax,1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          
          !C
          !C T2EXEC
          !C
          ALLOCATE(i2enr0(1:i0nmax,1:4     ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2wrks(1:i0nmax,1:i0wmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jmpm(1:i0dmax,1:i0dmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2kwns(1:i0nmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4smat(1:i0nmax,1:i0nmax,1:i0vmax,1:i0vmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2svec(1:i0nmax,1:i0vmax),STAT=i0err);IF(i0err.NE.0)EXIT
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
    
  END SUBROUTINE T2COMM_ALLOCATE
  
  SUBROUTINE T2COMM_DEALLOCATE
    
    !C
    !C T2INTG
    !C

    IF(ALLOCATED(d3imsn)) DEALLOCATE(d3imsn)
    IF(ALLOCATED(d4iavn)) DEALLOCATE(d4iavn)
    IF(ALLOCATED(d6iatn)) DEALLOCATE(d6iatn)
    IF(ALLOCATED(d5idtn)) DEALLOCATE(d5idtn)
    IF(ALLOCATED(d4igvn)) DEALLOCATE(d4igvn)
    IF(ALLOCATED(d6igtn)) DEALLOCATE(d6igtn)
    IF(ALLOCATED(d3iesn)) DEALLOCATE(d3iesn)
    IF(ALLOCATED(d5ievn)) DEALLOCATE(d5ievn)
    IF(ALLOCATED(d7ietn)) DEALLOCATE(d7ietn)
    IF(ALLOCATED(d2issn)) DEALLOCATE(d2issn)

    IF(ALLOCATED(d1wfct)) DEALLOCATE(d1wfct)
    IF(ALLOCATED(d1absc)) DEALLOCATE(d1absc)
    IF(ALLOCATED(d2wfct)) DEALLOCATE(d2wfct)
    IF(ALLOCATED(d4ifnc)) DEALLOCATE(d4ifnc)
    
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
    
!    IF(ALLOCATED(i1mfc1)) DEALLOCATE(i1mfc1)
!    IF(ALLOCATED(d2mfc1)) DEALLOCATE(d2mfc1)
    

    !C
    !C I2VGRA
    !C

    IF(ALLOCATED(i2msvt)) DEALLOCATE(i2msvt)
    IF(ALLOCATED(i2avvt)) DEALLOCATE(i2avvt)
    IF(ALLOCATED(i2atvt)) DEALLOCATE(i2atvt)
    IF(ALLOCATED(i2dtvt)) DEALLOCATE(i2dtvt)
    IF(ALLOCATED(i2gvvt)) DEALLOCATE(i2gvvt)
    IF(ALLOCATED(i2gtvt)) DEALLOCATE(i2gtvt)
    IF(ALLOCATED(i2esvt)) DEALLOCATE(i2esvt)
    IF(ALLOCATED(i2evvt)) DEALLOCATE(i2evvt)
    IF(ALLOCATED(i2etvt)) DEALLOCATE(i2etvt)
    IF(ALLOCATED(i2ssvt)) DEALLOCATE(i2ssvt)
    IF(ALLOCATED(i2vvvt)) DEALLOCATE(i2vvvt)
    
    IF(ALLOCATED(i3atwt)) DEALLOCATE(i3atwt)
    IF(ALLOCATED(i3gtwt)) DEALLOCATE(i3gtwt)
    IF(ALLOCATED(i3evwt)) DEALLOCATE(i3evwt)
    IF(ALLOCATED(i4etwt)) DEALLOCATE(i4etwt)
    
    !C
    !C T2MFCS
    !C

    IF(ALLOCATED(d2ug )) DEALLOCATE(d2ug )
    IF(ALLOCATED(d2jm1)) DEALLOCATE(d2jm1)
    IF(ALLOCATED(d2rzm)) DEALLOCATE(d2rzm)
    IF(ALLOCATED(d2rzx)) DEALLOCATE(d2rzx)

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
    IF(ALLOCATED(d1vb)) DEALLOCATE(d1vb)
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

    
    !C
    !C T2STEP
    !C
    
    IF(ALLOCATED(i1nlct))       DEALLOCATE(i1nlct) 
    IF(ALLOCATED(d1rsdl))       DEALLOCATE(d1rsdl)
    IF(ALLOCATED(d3amat))       DEALLOCATE(d3amat)
    IF(ALLOCATED(d2bvec))       DEALLOCATE(d2bvec)
    IF(ALLOCATED(d2xvec))       DEALLOCATE(d2xvec)
    IF(ALLOCATED(d2xvec_befor)) DEALLOCATE(d2xvec_befor)
    IF(ALLOCATED(d2xvec_after)) DEALLOCATE(d2xvec_after)

    IF(ALLOCATED(d2xout))       DEALLOCATE(d2xout)
    !C
    !C T2EXEC
    !C

    IF(ALLOCATED(d4smat)) DEALLOCATE(d4smat)
    IF(ALLOCATED(d2svec)) DEALLOCATE(d2svec)

    IF(ALLOCATED(i2enr0 )) DEALLOCATE(i2enr0 )
    IF(ALLOCATED(d2wrks )) DEALLOCATE(d2wrks )
    IF(ALLOCATED(d2jmpm )) DEALLOCATE(d2jmpm )
    IF(ALLOCATED(d2kwns )) DEALLOCATE(d2kwns )

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

  END SUBROUTINE T2COMM_DEALLOCATE
  
END MODULE T2COMM
