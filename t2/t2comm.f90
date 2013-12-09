MODULE T2COMM
  
  USE T2CNST, ONLY: i0rkind, i0ikind, i0lmaxm, i0spcsm
  
  IMPLICIT NONE
  !C------------------------------------------------------------------
  !C
  !C
  !C
  !C
  !C------------------------------------------------------------------
  REAL(   i0rkind)::&
       d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
       d0nncst,d0frcst,d0fbcst,d0ftcst,&
       d0ppcst,d0qrcst,d0qbcst,d0qtcst,d0iar
  
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2INTG: MODIFIED 2013-09-07
  !C 
  !C------------------------------------------------------------------
  !C INTEGRATION ARRAYS BY GAUSSIAN QUADRATURE
  !C
  
  INTEGER(i0ikind)::i0dbg
  INTEGER(i0ikind)::i0amax0 ! NUMBER OF ABSCISSAS PAR DIRECTION
  
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE:: d3imsn0
  REAL(   i0rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE:: d4iavn0
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE:: d6iatn0
  REAL(   i0rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE:: d5idtn0
  REAL(   i0rkind),DIMENSION(:,:,:,:      ),ALLOCATABLE:: d4igvn0
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:  ),ALLOCATABLE:: d6igtn0
  REAL(   i0rkind),DIMENSION(:,:,:        ),ALLOCATABLE:: d3iesn0
  REAL(   i0rkind),DIMENSION(:,:,:,:,:    ),ALLOCATABLE:: d5ievn0
  REAL(   i0rkind),DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE:: d7ietn0
  REAL(   i0rkind),DIMENSION(:,:          ),ALLOCATABLE:: d2issn0

  !C
  !C WORKING ARRAY FOR GAUSSIAN INTEGRATION
  !C  D4IFNC: INTERPOLATION FUNCTIONS 
  !C  D2WFCT: WHEIGHT FACTOR
     
  REAL(i0rkind),DIMENSION(:      ),ALLOCATABLE:: d1wfct0
  REAL(i0rkind),DIMENSION(:      ),ALLOCATABLE:: d1absc0
  REAL(i0rkind),DIMENSION(:,:    ),ALLOCATABLE:: d2wfct0
  REAL(i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: d4ifnc0

  
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2NGRA
  !C
  !C------------------------------------------------------------------
  INTEGER(i0ikind) ::&
       i0nmax0,& !C NUMBER OF NODES PAR ELEMENT
       i0nmax1,& 
       i0nmax2,& 
       i0nmax3,&
       i0nmax4,&
       i0emax ,&
       i0hmax ,&
       i0stm2 ,&
       i0lmax, &
       i0nrmx, &
       i0ncmx, &
       i0ermx, &
       i0ecmx, &
       i0dmax0,& !C NUMBER OF DIMENSIONS
       i0dmax1,&
       i0tmax ,& 
       i0sflg ,&
       i0spcs ,&
       i0pdiv_number
  
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::& 
       i1nmax1,&
       i1nmax2,&
       i1emax, &
       i1mlel, &
       i1rdn1, &
       i1pdn1, &
       i1pdn2, &
       i1nidr, &
       i1nidc, &
       i1nidr1, & !C NODEGRAPH FOR 1D CONFIGURATION
       i1nidc1, & !C NODEGRAPH FOR 1D CONFIGURATION
       !i1nidr2, & !C NODEGRAPH FOR 2D CONFIGURATION
       !i1nidc2, & !C NODEGRAPH FOR 2D CONFIGURATION
       i1eidr, &
       i1eidc, &
       i1dbc2, &
       i1mfc1, &
       i1mfc4

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:)::& 
       i2crt,&
       i2hbc
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:,:)::&
       i3enr 

  REAL(i0rkind),ALLOCATABLE,DIMENSION(:)::&
       d1msiz,&
       d1rsiz,&
       d1psiz,&
       d1mfc4
  REAL(i0rkind),ALLOCATABLE,DIMENSION(:,:)::&
       d2mfc1
  
  !C------------------------------------------------------------------
  !C  
  !C                    FOR T2VGRA: MODIFIED 2013-09-07
  !C
  !C------------------------------------------------------------------
  INTEGER(i0ikind)::&
       i0vmax,&
       i0vgrmx,i0msrmx,i0avrmx,i0atrmx,i0dtrmx,&
       i0gvrmx,i0gtrmx,i0esrmx,i0evrmx,i0etrmx,i0ssrmx,&
       i0vgcmx,i0mscmx,i0avcmx,i0atcmx,i0dtcmx,&
       i0gvcmx,i0gtcmx,i0escmx,i0evcmx,i0etcmx,i0sscmx

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:)::&
       i2vtbl

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::&
       i1vgidr,i1msidr,i1avidr,i1atidr,i1dtidr,&
       i1gvidr,i1gtidr,i1esidr,i1evidr,i1etidr,i1ssidr,&
       i1vgidc,i1msidc,i1avidc,i1atidc,i1dtidc,&
       i1gvidc,i1gtidc,i1esidc,i1evidc,i1etidc,i1ssidc
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:)::&
       i2atidc,i2gtidc,i2evidc,i2etidc
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::&
       i1atws,i1gtws,i1evws

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:)::&
       i2etws

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2MFCS
  !C 
  !C------------------------------------------------------------------
  INTEGER(i0ikind)::&
       i0mfcs
  REAL(   i0rkind)::&
       d0rmjr,&
       d0rmnr
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2ug,  &
       d2rzc1,&
       d2rzc3,&
       d2jm1

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
  INTEGER(i0ikind),DIMENSION(-1:i0lmaxm)::&
       i1rdn2
  REAL(   i0rkind),DIMENSION(0:i0lmaxm)::&
       d1rec
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2STEP
  !C 
  !C------------------------------------------------------------------
  INTEGER(i0ikind)::&
       i0pmax,&
       i0xmax,&
       i0cmax,&
       i0bmax,&
       i0nlct
  REAL(   i0rkind)::&
       d0eps
  INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::&
       i1nlct
  REAL(   i0rkind),DIMENSION(:),ALLOCATABLE::&
       d1rsdl,&
       d1gsm,&
       d1grv,&
       d1guv,&
       d1guv_befor,&
       d1guv_after

  !C------------------------------------------------------------------
  !C
  !C                          FOR T2LOOP
  !C 
  !C------------------------------------------------------------------
  INTEGER(i0ikind)::&
       i0tstp,&
       i0wstp
  REAL(   i0rkind)::&
       d0time,&
       d0tstp,&
       d0tmax

  !C------------------------------------------------------------------
  !C
  !C                         FOR T2CALV
  !C
  !C------------------------------------------------------------------ 

  REAL(   i0rkind),DIMENSION(:,:    ),ALLOCATABLE::d2ms
  REAL(   i0rkind),DIMENSION(:,:,:  ),ALLOCATABLE::d3av
  REAL(   i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE::d4at
  REAL(   i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE::d4dt
  REAL(   i0rkind),DIMENSION(:,:,:  ),ALLOCATABLE::d3gv
  REAL(   i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE::d4gt
  REAL(   i0rkind),DIMENSION(:,:    ),ALLOCATABLE::d2es
  REAL(   i0rkind),DIMENSION(:,:,:  ),ALLOCATABLE::d3ev
  REAL(   i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE::d4et
  REAL(   i0rkind),DIMENSION(:,:    ),ALLOCATABLE::d2ss

  REAL(   i0rkind),DIMENSION(:  ),ALLOCATABLE::&
       d1e, d1m, d1n, d1p, d1t,&
       d1ur,d1up,d1ut,d1ub,d1u2,d1qb,d1qr,d1qp,d1qt,d1qbp,&
       d1vt,d1c1,d1c2,d1c3,d1c4,d1hx,d1ni,d1pi
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2x,d2y,d2z,d2bcf
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2f1,d2f2,d2f3,d2f4
  
  INTEGER(i0ikind)::i0xa
  !C------------------------------------------------------------------
  !C
  !C                         FOR T2EXEC
  !C
  !C------------------------------------------------------------------
  
  INTEGER(i0ikind)::&
       i0eid,i0lid,i0supg,i0vida,i0vidb,&
       i0msid,i0avid,i0atid,i0dtid,i0gvid,i0gtid,&
       i0esid,i0evid,i0etid,i0ssid,i0wmax
  
  INTEGER(i0ikind),DIMENSION(:,:),ALLOCATABLE::&
       i2enr0
  REAL(   i0rkind)::&
       d0jdmp
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2jmpm
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2smat0
  REAL(   i0rkind),DIMENSION(:  ),ALLOCATABLE::&
       d1svec0
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2knv0
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2ws0,d2ws

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
  
    SUBROUTINE T2NGRA_ALLOCATE_1
    INTEGER(i0ikind),SAVE::i0lmax_save=0,i0spcs_save=0
    INTEGER(i0ikind)     :: i0err
    IF((i0lmax.NE.i0lmax_save).OR.&
       (i0spcs.NE.i0spcs_save))THEN
       
       IF(i0lmax_save.NE.0) CALL T2NGRA_DEALLOCATE_1
       
       DO 
          ALLOCATE(i1nmax1(1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1nmax2(1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1emax( 0:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
!          ALLOCATE(i1mlvl( 0:i0lmax+1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1rdn1(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1pdn1(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
!          ALLOCATE(i1rdn2(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1pdn2(-1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
!          ALLOCATE(d1rec(  0:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1msiz( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1rsiz( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1psiz( 1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          
          i0lmax_save=i0lmax

          RETURN
       
       ENDDO
       
       WRITE(6,*)'T2GRID_ALLOCATE_1: ALLOCATION ERROR: ECODE=',i0err
       STOP
       
    END IF
    
    RETURN
  END SUBROUTINE T2NGRA_ALLOCATE_1

  SUBROUTINE T2NGRA_DEALLOCATE_1

    IF(ALLOCATED(i1nmax1)) DEALLOCATE(i1nmax1)
    IF(ALLOCATED(i1nmax2)) DEALLOCATE(i1nmax2)
    IF(ALLOCATED(i1emax )) DEALLOCATE(i1emax )
!    IF(ALLOCATED(i1mlvl )) DEALLOCATE(i1mlvl )
    IF(ALLOCATED(i1rdn1 )) DEALLOCATE(i1rdn1 )
    IF(ALLOCATED(i1pdn1 )) DEALLOCATE(i1pdn1 )
!    IF(ALLOCATED(i1rdn2 )) DEALLOCATE(i1rdn2 )
    IF(ALLOCATED(i1pdn2 )) DEALLOCATE(i1pdn2 )
!    IF(ALLOCATED(d1rec  )) DEALLOCATE(d1rec  )
    IF(ALLOCATED(d1msiz )) DEALLOCATE(d1msiz )
    IF(ALLOCATED(d1rsiz )) DEALLOCATE(d1rsiz )
    IF(ALLOCATED(d1psiz )) DEALLOCATE(d1psiz )

    RETURN

  END SUBROUTINE T2NGRA_DEALLOCATE_1

  SUBROUTINE T2COMM_ALLOCATE
    INTEGER(i0ikind),SAVE::&
         i0nmax0_save=0,i0nmax1_save=0,i0nmax2_save=0,&
         i0nmax3_save=0,i0dmax0_save=0,i0amax0_save=0,&
         i0emax_save =0,&
         i0nrmx_save =0,i0ncmx_save =0,&
         i0ermx_save =0,i0ecmx_save =0,&
         i0hmax_save =0,i0vmax_save =0,&
         i0vgrmx_save=0,i0msrmx_save=0,i0avrmx_save=0,&
         i0atrmx_save=0,i0dtrmx_save=0,i0gvrmx_save=0,&
         i0gtrmx_save=0,i0esrmx_save=0,i0evrmx_save=0,&
         i0etrmx_save=0,i0ssrmx_save=0,i0vgcmx_save=0,&
         i0mscmx_save=0,i0avcmx_save=0,i0atcmx_save=0,&
         i0dtcmx_save=0,i0gvcmx_save=0,i0gtcmx_save=0,&
         i0escmx_save=0,i0evcmx_save=0,i0etcmx_save=0,&
         i0sscmx_save=0,i0spcs_save =0,i0tmax_save =0,&
         i0xmax_save =0,i0bmax_save =0,i0cmax_save =0

    INTEGER(i0ikind)     :: i0err

    IF(  (i0nmax0.NE.i0nmax0_save).OR.&
         (i0dmax0.NE.i0dmax0_save).OR.&
         (i0amax0.NE.i0amax0_save).OR.&
         (i0nmax1.NE.i0nmax1_save).OR.&
         (i0nmax2.NE.i0nmax2_save).OR.&
         (i0nmax3.NE.i0nmax3_save).OR.&
         (i0emax .NE.i0emax_save ).OR.&
         (i0nrmx .NE.i0nrmx_save ).OR.&
         (i0ncmx .NE.i0ncmx_save ).OR.&
         (i0ermx .NE.i0ermx_save ).OR.&
         (i0ecmx .NE.i0ecmx_save ).OR.&
         (i0hmax .NE.i0hmax_save ).OR.&
         (i0vmax .NE.i0vmax_save ).OR.&
         (i0vgrmx.NE.i0vgrmx_save).OR.&
         (i0msrmx.NE.i0msrmx_save).OR.&
         (i0avrmx.NE.i0avrmx_save).OR.&
         (i0atrmx.NE.i0atrmx_save).OR.&
         (i0dtrmx.NE.i0dtrmx_save).OR.&
         (i0gvrmx.NE.i0gvrmx_save).OR.&
         (i0gtrmx.NE.i0gtrmx_save).OR.&
         (i0esrmx.NE.i0esrmx_save).OR.&
         (i0evrmx.NE.i0evrmx_save).OR.&
         (i0etrmx.NE.i0etrmx_save).OR.&
         (i0ssrmx.NE.i0ssrmx_save).OR.&
         (i0vgcmx.NE.i0vgcmx_save).OR.&
         (i0mscmx.NE.i0mscmx_save).OR.&
         (i0avcmx.NE.i0avcmx_save).OR.&
         (i0atcmx.NE.i0atcmx_save).OR.&
         (i0dtcmx.NE.i0dtcmx_save).OR.&
         (i0gvcmx.NE.i0gvcmx_save).OR.&
         (i0gtcmx.NE.i0gtcmx_save).OR.&
         (i0escmx.NE.i0escmx_save).OR.&
         (i0evcmx.NE.i0evcmx_save).OR.&
         (i0etcmx.NE.i0etcmx_save).OR.&
         (i0sscmx.NE.i0sscmx_save).OR.&
         (i0spcs .NE.i0spcs_save ).OR.&
         (i0tmax .NE.i0tmax_save ).OR.&
         (i0cmax .NE.i0cmax_save ).OR.&
         (i0bmax .NE.i0bmax_save ).OR.&
         (i0xmax .NE.i0xmax_save ))THEN
       
       IF(i0spcs_save.NE.0) CALL T2COMM_DEALLOCATE

       DO
          !C T2INTG
          ALLOCATE(d3imsn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4iavn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d6iatn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0,1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d5idtn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0,1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4igvn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d6igtn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0,1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d3iesn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d5ievn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d7ietn0(&
               1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
               1:i0dmax0,1:i0dmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d2issn0(&
               1:i0nmax0,1:i0nmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          
          ALLOCATE(d1wfct0(1:i0amax0),STAT=i0err);IF(i0err.NE.0) EXIT
          ALLOCATE(d1absc0(1:i0amax0),STAT=i0err);IF(i0err.NE.0) EXIT
          ALLOCATE(d2wfct0(&
               1:i0amax0,1:i0amax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4ifnc0(&
               1:i0nmax0,0:i0dmax0,&
               1:i0amax0,1:i0amax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT

          !C T2NGRA
          ALLOCATE(i1nidr( 1:i0nrmx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1nidc( 1:i0ncmx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1eidr( 1:i0ermx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1eidc( 1:i0ecmx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mlel( 1:i0emax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2crt(  1:i0nmax1,1:2),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(i1dbc2( 1:i0nmax2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mfc1( 1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mfc4( 1:i0nmax4),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1mfc4( 1:i0nmax4),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i3enr(  1:i0emax,1:4,1:i0nmax0),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(d2mfc1( 1:i0nmax1,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          IF(i0hmax.NE.0)&
               ALLOCATE(i2hbc(1:i0hmax,1:2),STAT=i0err);&
               IF(i0err.NE.0)EXIT

          !C T2VGRA
          ALLOCATE(i2vtbl( 1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i1vgidr(1:i0vgrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1msidr(1:i0msrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1avidr(1:i0avrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1atidr(1:i0atrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1dtidr(1:i0dtrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1gvidr(1:i0gvrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1gtidr(1:i0gtrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1esidr(1:i0esrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1evidr(1:i0evrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1etidr(1:i0etrmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1ssidr(1:i0ssrmx),STAT=i0err);IF(i0err.NE.0)EXIT

          ALLOCATE(i1vgidc(1:i0vgcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1msidc(1:i0mscmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1avidc(1:i0avcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1atidc(1:i0atcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1dtidc(1:i0dtcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1gvidc(1:i0gvcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1gtidc(1:i0gtcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1esidc(1:i0escmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1evidc(1:i0evcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1etidc(1:i0etcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1ssidc(1:i0sscmx),STAT=i0err);IF(i0err.NE.0)EXIT
       
          ALLOCATE(i2atidc(1:i0atcmx,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2gtidc(1:i0gtcmx,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2evidc(1:i0evcmx,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2etidc(1:i0evcmx,1:3),STAT=i0err);IF(i0err.NE.0)EXIT

          ALLOCATE(i1atws(    1:i0atcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1gtws(    1:i0gtcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1evws(    1:i0evcmx),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2etws(1:2,1:i0etcmx),STAT=i0err);IF(i0err.NE.0)EXIT

          !C T2MFCS
          ALLOCATE(d2ug(  1:2,1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jm1( 1:5,1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzc1(1:i0nmax1,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzc3(1:i0nmax3,1:2),STAT=i0err);IF(i0err.NE.0)EXIT

          !C T2CALV
          ALLOCATE(d2ms(                    1:i0mscmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3av(1:i0dmax0,          1:i0avcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4at(1:i0dmax0,1:i0dmax0,1:i0atcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4dt(1:i0dmax0,1:i0dmax0,1:i0dtcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3gv(1:i0dmax0,          1:i0gvcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4gt(1:i0dmax0,1:i0dmax0,1:i0gtcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2es(                    1:i0escmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3ev(1:i0dmax0,          1:i0evcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4et(1:i0dmax0,1:i0dmax0,1:i0etcmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2ss(                    1:i0mscmx,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d1e( 1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1m( 1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1n( 1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ur(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1up(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ut(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ub(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1u2(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1p( 1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qb(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qr(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qp(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qt(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qbp(1:i0spcs        ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1t( 1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1vt(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ni(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1pi(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1c1(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1c2(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1c3(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1c4(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1hx(1:i0spcs         ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2f1(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2f2(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2f3(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2f4(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2x(  1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2y(  1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2z(  1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2bcf(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT

          ALLOCATE(d1bp3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1bt3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1er3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ep3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1et3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2n3( 1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2fr3(1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2fb3(1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2ft3(1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2p3( 1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2qr3(1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2qb3(1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2qt3(1:i0nmax3,1:i0spcs),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1nlct(     1:10000),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1rsdl(     1:10000),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1gsm(      1:i0cmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1grv(      1:i0bmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1guv(      1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1guv_befor(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1guv_after(1:i0xmax),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2ws(  1:i0wmax, 1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2enr0( 1:4,      1:i0nmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2ws0(  1:i0wmax, 1:i0nmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jmpm( 1:i0dmax0,1:i0dmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2knv0( 1:i0nmax0,1:i0vmax ),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2smat0(1:i0nmax0,1:i0nmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1svec0(1:i0nmax0),STAT=i0err);IF(i0err.NE.0)EXIT

          !C TMP
          ALLOCATE(d1jm1(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm2(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm4(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm5(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT

          i0spcs_save  = i0spcs
          i0nmax0_save = i0nmax0
          i0dmax0_save = i0dmax0
          i0amax0_save = i0amax0
          i0nmax1_save = i0nmax1
          i0nmax2_save = i0nmax2
          i0nmax3_save = i0nmax3
          i0emax_save  = i0emax
          i0nrmx_save  = i0nrmx
          i0ncmx_save  = i0ncmx
          i0hmax_save  = i0hmax 
          i0vmax_save  = i0vmax
          i0vgrmx_save = i0vgrmx
          i0msrmx_save = i0msrmx
          i0avrmx_save = i0avrmx
          i0atrmx_save = i0atrmx
          i0dtrmx_save = i0dtrmx
          i0gvrmx_save = i0gvrmx
          i0gtrmx_save = i0gtrmx
          i0esrmx_save = i0esrmx
          i0etrmx_save = i0evrmx
          i0ssrmx_save = i0ssrmx
          i0vgcmx_save = i0vgcmx
          i0mscmx_save = i0mscmx
          i0avcmx_save = i0avcmx
          i0atcmx_save = i0atcmx
          i0dtcmx_save = i0dtcmx
          i0gvcmx_save = i0gvcmx
          i0gtcmx_save = i0gtcmx
          i0escmx_save = i0escmx
          i0evcmx_save = i0evcmx
          i0etcmx_save = i0etcmx
          i0sscmx_save = i0sscmx
          i0tmax_save = i0tmax
          i0cmax_save = i0cmax
          i0bmax_save = i0bmax
          i0xmax_save = i0xmax

          RETURN

       ENDDO

       WRITE(6,*)'T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',i0err
       STOP

    END IF

    RETURN

  END SUBROUTINE T2COMM_ALLOCATE
  
  SUBROUTINE T2COMM_DEALLOCATE

    IF(ALLOCATED(i1nmax1)) DEALLOCATE(i1nmax1)
    IF(ALLOCATED(i1nmax2)) DEALLOCATE(i1nmax2)
    IF(ALLOCATED(i1emax )) DEALLOCATE(i1emax )
!    IF(ALLOCATED(i1mlvl )) DEALLOCATE(i1mlvl )
    IF(ALLOCATED(i1rdn1 )) DEALLOCATE(i1rdn1 )
    IF(ALLOCATED(i1pdn1 )) DEALLOCATE(i1pdn1 )
!    IF(ALLOCATED(i1rdn2 )) DEALLOCATE(i1rdn2 )
    IF(ALLOCATED(i1pdn2 )) DEALLOCATE(i1pdn2 )
!    IF(ALLOCATED(d1rec  )) DEALLOCATE(d1rec  )
    IF(ALLOCATED(d1msiz )) DEALLOCATE(d1msiz )
    IF(ALLOCATED(d1rsiz )) DEALLOCATE(d1rsiz )
    IF(ALLOCATED(d1psiz )) DEALLOCATE(d1psiz )
    
    IF(ALLOCATED(i1nidr )) DEALLOCATE(i1nidr )
    IF(ALLOCATED(i1nidc )) DEALLOCATE(i1nidc )
    IF(ALLOCATED(i1eidr )) DEALLOCATE(i1eidr )
    IF(ALLOCATED(i1eidc )) DEALLOCATE(i1eidc )
    IF(ALLOCATED(i2crt  )) DEALLOCATE(i2crt  )
    IF(ALLOCATED(i1dbc2 )) DEALLOCATE(i1dbc2 )
    IF(ALLOCATED(i1mfc1 )) DEALLOCATE(i1mfc1 )
    IF(ALLOCATED(i1mfc4 )) DEALLOCATE(i1mfc4 )
    IF(ALLOCATED(d1mfc4 )) DEALLOCATE(d1mfc4 )
    IF(ALLOCATED(i3enr  )) DEALLOCATE(i3enr  )
    IF(ALLOCATED(i1mlel )) DEALLOCATE(i1mlel )
    IF(ALLOCATED(d2ug   )) DEALLOCATE(d2ug   )
    IF(ALLOCATED(d2mfc1 )) DEALLOCATE(d2mfc1 )
    IF(ALLOCATED(d2rzc1 )) DEALLOCATE(d2rzc1 )
    IF(ALLOCATED(d2rzc3 )) DEALLOCATE(d2rzc3 )
    IF(ALLOCATED(i2hbc  )) DEALLOCATE(i2hbc  )
    
    IF(ALLOCATED(i2vtbl )) DEALLOCATE(i2vtbl )
    IF(ALLOCATED(i1vgidr)) DEALLOCATE(i1vgidr)
    IF(ALLOCATED(i1msidr)) DEALLOCATE(i1msidr)
    IF(ALLOCATED(i1avidr)) DEALLOCATE(i1avidr)
    IF(ALLOCATED(i1atidr)) DEALLOCATE(i1atidr)
    IF(ALLOCATED(i1dtidr)) DEALLOCATE(i1dtidr)
    IF(ALLOCATED(i1gvidr)) DEALLOCATE(i1gvidr)
    IF(ALLOCATED(i1gtidr)) DEALLOCATE(i1gtidr)
    IF(ALLOCATED(i1esidr)) DEALLOCATE(i1esidr)
    IF(ALLOCATED(i1evidr)) DEALLOCATE(i1evidr)
    IF(ALLOCATED(i1etidr)) DEALLOCATE(i1etidr)
    IF(ALLOCATED(i1ssidr)) DEALLOCATE(i1ssidr)
    
    IF(ALLOCATED(i1vgidc)) DEALLOCATE(i1vgidc)
    IF(ALLOCATED(i1msidc)) DEALLOCATE(i1msidc)
    IF(ALLOCATED(i1avidc)) DEALLOCATE(i1avidc)
    IF(ALLOCATED(i1atidc)) DEALLOCATE(i1atidc)
    IF(ALLOCATED(i1dtidc)) DEALLOCATE(i1dtidc)
    IF(ALLOCATED(i1gvidc)) DEALLOCATE(i1gvidc)
    IF(ALLOCATED(i1gtidc)) DEALLOCATE(i1gtidc)
    IF(ALLOCATED(i1esidc)) DEALLOCATE(i1esidc)
    IF(ALLOCATED(i1evidc)) DEALLOCATE(i1evidc)
    IF(ALLOCATED(i1etidc)) DEALLOCATE(i1etidc)
    IF(ALLOCATED(i1ssidc)) DEALLOCATE(i1ssidc)

    IF(ALLOCATED(i2atidc)) DEALLOCATE(i2atidc)
    IF(ALLOCATED(i2gtidc)) DEALLOCATE(i2gtidc)
    IF(ALLOCATED(i2evidc)) DEALLOCATE(i2evidc)
    IF(ALLOCATED(i2etidc)) DEALLOCATE(i2etidc)
    
    IF(ALLOCATED(i1atws )) DEALLOCATE(i1atws )
    IF(ALLOCATED(i1gtws )) DEALLOCATE(i1gtws )
    IF(ALLOCATED(i1evws )) DEALLOCATE(i1evws )
    IF(ALLOCATED(i2etws )) DEALLOCATE(i2etws )
    
    IF(ALLOCATED(d2jm1  )) DEALLOCATE(d2jm1  )
    IF(ALLOCATED(d2rzc1 )) DEALLOCATE(d2rzc1 )
    IF(ALLOCATED(d2rzc3 )) DEALLOCATE(d2rzc3 )
    
    IF(ALLOCATED(d3imsn0)) DEALLOCATE(d3imsn0)
    IF(ALLOCATED(d4iavn0)) DEALLOCATE(d4iavn0)
    IF(ALLOCATED(d6iatn0)) DEALLOCATE(d6iatn0)
    IF(ALLOCATED(d5idtn0)) DEALLOCATE(d5idtn0)
    IF(ALLOCATED(d4igvn0)) DEALLOCATE(d4igvn0)
    IF(ALLOCATED(d6igtn0)) DEALLOCATE(d6igtn0)
    IF(ALLOCATED(d3iesn0)) DEALLOCATE(d3iesn0)
    IF(ALLOCATED(d5ievn0)) DEALLOCATE(d5ievn0)
    IF(ALLOCATED(d7ietn0)) DEALLOCATE(d7ietn0)
    IF(ALLOCATED(d2issn0)) DEALLOCATE(d2issn0)
    IF(ALLOCATED(d1wfct0)) DEALLOCATE(d1wfct0)
    IF(ALLOCATED(d1absc0)) DEALLOCATE(d1absc0)
    IF(ALLOCATED(d2wfct0)) DEALLOCATE(d2wfct0)
    IF(ALLOCATED(d4ifnc0)) DEALLOCATE(d4ifnc0)

    IF(ALLOCATED(d2ms)) DEALLOCATE(d2ms)
    IF(ALLOCATED(d3av)) DEALLOCATE(d3av)
    IF(ALLOCATED(d4at)) DEALLOCATE(d4at)
    IF(ALLOCATED(d4dt)) DEALLOCATE(d4dt)
    IF(ALLOCATED(d3gv)) DEALLOCATE(d3gv)
    IF(ALLOCATED(d4gt)) DEALLOCATE(d4gt)
    IF(ALLOCATED(d2es)) DEALLOCATE(d2es)
    IF(ALLOCATED(d3ev)) DEALLOCATE(d3ev)
    IF(ALLOCATED(d4et)) DEALLOCATE(d4et)
    IF(ALLOCATED(d2ss)) DEALLOCATE(d2ss)

    IF(ALLOCATED(d1e )) DEALLOCATE(d1e )
    IF(ALLOCATED(d1m )) DEALLOCATE(d1m )
    IF(ALLOCATED(d1n )) DEALLOCATE(d1n )
    IF(ALLOCATED(d1ni)) DEALLOCATE(d1ni)
    IF(ALLOCATED(d1pi)) DEALLOCATE(d1pi)
    IF(ALLOCATED(d1ur)) DEALLOCATE(d1ur)
    IF(ALLOCATED(d1up)) DEALLOCATE(d1up)
    IF(ALLOCATED(d1ut)) DEALLOCATE(d1ut)
    IF(ALLOCATED(d1u2)) DEALLOCATE(d1u2)
    IF(ALLOCATED(d1p )) DEALLOCATE(d1p )
    IF(ALLOCATED(d1qb)) DEALLOCATE(d1qb)
    IF(ALLOCATED(d1qr)) DEALLOCATE(d1qr)
    IF(ALLOCATED(d1qp)) DEALLOCATE(d1qp)
    IF(ALLOCATED(d1qt)) DEALLOCATE(d1qt)
    IF(ALLOCATED(d1qbp)) DEALLOCATE(d1qbp)
    IF(ALLOCATED(d1t )) DEALLOCATE(d1t )
    IF(ALLOCATED(d1vt)) DEALLOCATE(d1vt)
    IF(ALLOCATED(d1c1)) DEALLOCATE(d1c1)
    IF(ALLOCATED(d1c2)) DEALLOCATE(d1c2)
    IF(ALLOCATED(d1c3)) DEALLOCATE(d1c3)
    IF(ALLOCATED(d1c4)) DEALLOCATE(d1c4)
    IF(ALLOCATED(d1hx)) DEALLOCATE(d1hx)
    IF(ALLOCATED(d2f1)) DEALLOCATE(d2f1)
    IF(ALLOCATED(d2f2)) DEALLOCATE(d2f2)
    IF(ALLOCATED(d2f3)) DEALLOCATE(d2f3)
    IF(ALLOCATED(d2f4)) DEALLOCATE(d2f4)
    IF(ALLOCATED(d2x )) DEALLOCATE(d2x )
    IF(ALLOCATED(d2y )) DEALLOCATE(d2y )
    IF(ALLOCATED(d2z )) DEALLOCATE(d2z )
    IF(ALLOCATED(d2bcf)) DEALLOCATE(d2bcf)
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
    IF(ALLOCATED(i1nlct     )) DEALLOCATE(i1nlct     ) 
    IF(ALLOCATED(d1rsdl     )) DEALLOCATE(d1rsdl     )
    IF(ALLOCATED(d1gsm      )) DEALLOCATE(d1gsm      )
    IF(ALLOCATED(d1grv      )) DEALLOCATE(d1grv      )
    IF(ALLOCATED(d1guv      )) DEALLOCATE(d1guv      )
    IF(ALLOCATED(d1guv_befor)) DEALLOCATE(d1guv_befor)
    IF(ALLOCATED(d1guv_after)) DEALLOCATE(d1guv_after)
    IF(ALLOCATED(d2ws       )) DEALLOCATE(d2ws       )

    IF(ALLOCATED(i2enr0 )) DEALLOCATE(i2enr0 )
    IF(ALLOCATED(d2ws0  )) DEALLOCATE(d2ws0  )
    IF(ALLOCATED(d2jmpm )) DEALLOCATE(d2jmpm )
    IF(ALLOCATED(d2knv0 )) DEALLOCATE(d2knv0 )
    IF(ALLOCATED(d2smat0)) DEALLOCATE(d2smat0)
    IF(ALLOCATED(d1svec0)) DEALLOCATE(d1svec0)

   
    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE
  
END MODULE T2COMM
