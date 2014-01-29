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
  REAL(i0rkind):: &
       time_t2                            ! global time
  REAL(i0rkind),DIMENSION(:,:),ALLOCATABLE:: &
       vv0d,vv2d                          ! storage for 0D and 2D data

  !C---------------------------------------------------------
  !C
  !C
  !C
  !C
  !C---------------------------------------------------------
  REAL(   i0rkind)::&
       d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
       d0nncst,d0frcst,d0fbcst,d0ftcst,&
       d0ppcst,d0qrcst,d0qbcst,d0qtcst,d0iar
  
  !C--------------------------------------------------------
  !C
  !C                 FOR T2INTG: MODIFIED 2014-01-29
  !C 
  !C--------------------------------------------------------
  !C INTEGRATION ARRAYS BY GAUSSIAN QUADRATURE
  !C
  
  INTEGER(i0ikind)::i0dbg
  INTEGER(i0ikind)::i0amax ! NUMBER OF ABSCISSAS PAR DIRECTION
  
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
  
  
  !C------------------------------------------------------------------
  !C
  !C                          FOR T2NGRA
  !C
  !C------------------------------------------------------------------
  INTEGER(i0ikind) ::&
       i0nmax ,&
       i0nmax0,& !C NUMBER OF NODES PAR ELEMENT
       i0nmax1,& 
       i0nmax2,& 
       i0nmax3,&
       i0nmax4,&
       i0emax ,& ! total number of elements
       i0hmax ,&
       i0stm2 ,&
       i0lmax, &
       i0nrmx, &
       i0ncmx, &
       i0ermx, &
       i0ecmx, &
       i0dmax0,& !C NUMBER OF DIMENSIONS
       i0dmax1,&
       i0sflg ,&
       i0spcs ,&
       i0smax,& !XXXXXXXXXXXXXXXXXXXXX
       i0pdiv_number
  
  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:)::& 
       i1nmax1,&
       i1nmax2,&
       i1emax, & ! number of elements in the level
       i1mlel, &
       i1rdn1, &
       i1pdn1, &
       i1pdn2, &
       i1nidr, &
       i1nidc, &
       !i1nidr1, & !C NODEGRAPH FOR 1D CONFIGURATION
       !i1nidc1, & !C NODEGRAPH FOR 1D CONFIGURATION
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
  !C                    FOR T2VGRA: MODIFIED 2014-01-28
  !C
  !C------------------------------------------------------------------
  INTEGER(i0ikind)::&
       i0vmax

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:)::&
       i2msvt,i2avvt,i2atvt,i2dtvt,i2gvvt,i2gtvt,&
       i2esvt,i2evvt,i2etvt,i2ssvt,i2vvvt

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:,:)::&
       i3atwt,i3gtwt,i3evwt

  INTEGER(i0ikind),ALLOCATABLE,DIMENSION(:,:,:,:)::&
       i4etwt

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
!       i0tstp,&
       i0wstp
!  REAL(   i0rkind)::&
!       d0time,&
!       d0tstp,&
!       d0tmax

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

  REAL(   i0rkind),DIMENSION(:),ALLOCATABLE::&
       d1ee,d1mm,d1nn,d1ni,d1pp,d1pi,d1tt,d1ti,&
       d1ur,d1up,d1ut,d1ub,d1u2,&
       d1qr,d1qp,d1qt,d1qb,d1vb,d1vt,&
       d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,d1hex
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2x,d2y,d2z,d2bcf
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4
  INTEGER(i0ikind)::&
       i0wmax,i0dmax
  INTEGER(i0ikind)::i0xa
  !C------------------------------------------------------------------
  !C
  !C                         FOR T2EXEC
  !C
  !C------------------------------------------------------------------
  
  INTEGER(i0ikind)::&
       i0eid,i0lid,i0supg
  
  INTEGER(i0ikind),DIMENSION(:,:),ALLOCATABLE::&
       i2enr0
  REAL(   i0rkind)::&
       d0jdmp
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2jmpm
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2smat0,d2smat

  REAL(   i0rkind),DIMENSION(:,:,:,:),ALLOCATABLE::&
       d4smat
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2svec
  REAL(   i0rkind),DIMENSION(:  ),ALLOCATABLE::&
       d1svec0,d1svec
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2kwns,d2knv
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2wrks,d2ws
  REAL(   i0rkind),DIMENSION(:,:,:),ALLOCATABLE::&
       d3gmat
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
  
    SUBROUTINE T2NGRA_ALLOCATE
    INTEGER(i0ikind),SAVE::i0lmax_save=0
    INTEGER(i0ikind)     :: i0err
    IF(i0lmax.NE.i0lmax_save) THEN
       
       CALL T2NGRA_DEALLOCATE
       
       DO 
          ALLOCATE(i1nmax1(1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1nmax2(1:i0lmax  ),STAT=i0err);IF(i0err.NE.0)EXIT
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

    IF(ALLOCATED(i1nmax1)) DEALLOCATE(i1nmax1)
    IF(ALLOCATED(i1nmax2)) DEALLOCATE(i1nmax2)
    IF(ALLOCATED(i1emax )) DEALLOCATE(i1emax )
    IF(ALLOCATED(i1rdn1 )) DEALLOCATE(i1rdn1 )
    IF(ALLOCATED(i1pdn1 )) DEALLOCATE(i1pdn1 )
    IF(ALLOCATED(i1pdn2 )) DEALLOCATE(i1pdn2 )
    IF(ALLOCATED(d1msiz )) DEALLOCATE(d1msiz )
    IF(ALLOCATED(d1rsiz )) DEALLOCATE(d1rsiz )
    IF(ALLOCATED(d1psiz )) DEALLOCATE(d1psiz )

    RETURN

  END SUBROUTINE T2NGRA_DEALLOCATE

  SUBROUTINE T2COMM_ALLOCATE
    INTEGER(i0ikind),SAVE::&
         i0nmax0_save=0,i0nmax1_save=0,i0nmax2_save=0,&
         i0nmax3_save=0,i0dmax0_save=0,i0amax_save=0,&
         i0emax_save =0,&
         i0nrmx_save =0,i0ncmx_save =0,&
         i0ermx_save =0,i0ecmx_save =0,&
         i0hmax_save =0,i0vmax_save =0,&
         i0spcs_save =0,i0xmax_save =0,&
         i0bmax_save =0,i0cmax_save =0
    INTEGER(i0ikind),SAVE::&
         nv0dmax_save=0,nt0dmax_save=0, &
         nv2dmax_save=0,nt2dmax_save=0

    INTEGER(i0ikind)     :: i0err

    IF(  (i0spcs .NE.i0spcs_save ).OR.&
         (i0nmax0.NE.i0nmax0_save).OR.&
         (i0dmax0.NE.i0dmax0_save).OR.&
         (i0amax.NE.i0amax_save).OR.&
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
         (i0cmax .NE.i0cmax_save ).OR.&
         (i0bmax .NE.i0bmax_save ).OR.&
         (i0xmax .NE.i0xmax_save ).OR.&
         (nv0dmax.NE.nv0dmax_save).OR.&
         (nt0dmax.NE.nt0dmax_save).OR.&
         (nv2dmax.NE.nv2dmax_save).OR.&
         (nt2dmax.NE.nt2dmax_save))THEN
       
       CALL T2COMM_DEALLOCATE

       DO
          !C T2INTG
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
          ALLOCATE(d2issn(1:i0nmax0,1:i0nmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          
          ALLOCATE(d1wfct(1:i0amax),STAT=i0err);IF(i0err.NE.0) EXIT
          ALLOCATE(d1absc(1:i0amax),STAT=i0err);IF(i0err.NE.0) EXIT
          ALLOCATE(d2wfct(1:i0amax,1:i0amax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(d4ifnc(1:i0amax,1:i0amax,0:i0dmax,1:i0nmax0),&
               STAT=i0err); IF(i0err.NE.0) EXIT

          !C
          !C T2NGRA
          !C

          ALLOCATE(i1nidr( 1:i0nrmx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1nidc( 1:i0ncmx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1eidr( 1:i0ermx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1eidc( 1:i0ecmx  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mlel( 1:i0emax  ),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i2crt(  1:i0nmax1,1:3),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(i1dbc2( 1:i0nmax2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mfc1( 1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i1mfc4( 1:i0nmax4),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1mfc4( 1:i0nmax4),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(i3enr(  1:i0nmax,1:4,1:i0emax),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(d2mfc1( 1:i0nmax1,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          IF(i0hmax.NE.0)&
               ALLOCATE(i2hbc(1:i0hmax,1:2),STAT=i0err);&
               IF(i0err.NE.0)EXIT

          !C
          !C T2VGRA
          !C

          ALLOCATE(i2msvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2avvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2atvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2dtvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2gvvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2gtvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2esvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2evvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2etvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2ssvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          ALLOCATE(i2vvvt(1:i0vmax,1:i0vmax),STAT=i0err)&
               ;IF(i0err.NE.0)EXIT
          
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

          ALLOCATE(d2ug(  1:2,1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jm1( 1:5,1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzc1(1:i0nmax1,1:2),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2rzc3(1:i0nmax3,1:2),STAT=i0err);IF(i0err.NE.0)EXIT

          !C
          !C T2CALV
          !C

          ALLOCATE(d3ms(1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4av(1:i0dmax,1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d6at(1:i0dmax,1:i0dmax,1:i0wmax,1:i0vmax,1:i0vmax,&
               1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d5dt(1:i0dmax,1:i0dmax,1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d4gv(1:i0dmax,1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d6gt(1:i0dmax,1:i0dmax,1:i0wmax,1:i0vmax,1:i0vmax,&
               1:i0nmax1),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3es(1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d5ev(1:i0dmax,1:i0wmax,1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d7et(1:i0dmax,1:i0dmax,1:i0wmax,1:i0wmax,&
               1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3ss(1:i0vmax,1:i0vmax,1:i0nmax1),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d1ee(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1mm(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nn(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ur(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1up(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ut(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ub(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1u2(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1pp(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qr(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qp(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qt(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1qb(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1vb(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1tt(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1vt(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ni(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1pi(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1ti(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT

          ALLOCATE(d1nvcc1(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nvcc2(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nvcc3(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1nvcc4(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d1hex(1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d2nfcf1(1:i0spcs,1:i0spcs),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf2(1:i0spcs,1:i0spcs),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf3(1:i0spcs,1:i0spcs),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          ALLOCATE(d2nfcf4(1:i0spcs,1:i0spcs),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          
          ALLOCATE(d2x(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2y(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2z(1:i0spcs,1:i0spcs),STAT=i0err);IF(i0err.NE.0)EXIT
          
          ALLOCATE(d2bcf(1:i0spcs,1:i0spcs),STAT=i0err);&
               IF(i0err.NE.0)EXIT
          
          !C FOR T2WRIT
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
          ALLOCATE(i2enr0( 1:i0nmax,1:4),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2wrks(  1:i0wmax, 1:i0nmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2jmpm( 1:i0dmax0,1:i0dmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2kwns( 1:i0nmax0,1:i0vmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2smat0(1:i0nmax0,1:i0nmax0),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1svec0(1:i0nmax0),STAT=i0err);IF(i0err.NE.0)EXIT

          !C T2EXEC NEW
          ALLOCATE(d4smat(1:i0nmax,1:i0nmax,1:i0vmax,1:i0vmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d2svec(1:i0nmax,1:i0vmax),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d3gmat(1:i0vmax,1:i0vmax,1:i0ncmx),&
               STAT=i0err);IF(i0err.NE.0)EXIT
          !C TMP
          ALLOCATE(d1jm1(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm2(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm3(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm4(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT
          ALLOCATE(d1jm5(1:i0nmax3),STAT=i0err);IF(i0err.NE.0)EXIT

          ALLOCATE(vv0d(nv0dmax,nt0dmax),STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(vv2d(nv2dmax,nt2dmax),STAT=i0err); IF(i0err.NE.0) EXIT

          i0spcs_save  = i0spcs
          i0nmax0_save = i0nmax0
          i0dmax0_save = i0dmax0
          i0amax_save  = i0amax
          i0nmax1_save = i0nmax1
          i0nmax2_save = i0nmax2
          i0nmax3_save = i0nmax3
          i0emax_save  = i0emax
          i0nrmx_save  = i0nrmx
          i0ncmx_save  = i0ncmx
          i0ermx_save  = i0ermx
          i0ecmx_save  = i0ecmx
          i0hmax_save  = i0hmax 
          i0vmax_save  = i0vmax
          i0cmax_save = i0cmax
          i0bmax_save = i0bmax
          i0xmax_save = i0xmax

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
    !C FOR T2INTG
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

    IF(ALLOCATED(i1nidr)) DEALLOCATE(i1nidr)
    IF(ALLOCATED(i1nidc)) DEALLOCATE(i1nidc)
    IF(ALLOCATED(i1eidr)) DEALLOCATE(i1eidr)
    IF(ALLOCATED(i1eidc)) DEALLOCATE(i1eidc)
    IF(ALLOCATED(i1mlel)) DEALLOCATE(i1mlel)
    IF(ALLOCATED(i2crt )) DEALLOCATE(i2crt )
    IF(ALLOCATED(i1dbc2)) DEALLOCATE(i1dbc2)
    IF(ALLOCATED(i1mfc1)) DEALLOCATE(i1mfc1)
    IF(ALLOCATED(i1mfc4)) DEALLOCATE(i1mfc4)
    IF(ALLOCATED(d1mfc4)) DEALLOCATE(d1mfc4)
    IF(ALLOCATED(i3enr )) DEALLOCATE(i3enr )
    IF(ALLOCATED(d2mfc1)) DEALLOCATE(d2mfc1)
    IF(ALLOCATED(i2hbc )) DEALLOCATE(i2hbc )


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
    !C
    !C

    IF(ALLOCATED(d2ug   )) DEALLOCATE(d2ug   )
    IF(ALLOCATED(d2jm1  )) DEALLOCATE(d2jm1  )
    IF(ALLOCATED(d2rzc1 )) DEALLOCATE(d2rzc1 )
    IF(ALLOCATED(d2rzc3 )) DEALLOCATE(d2rzc3 )

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

    IF(ALLOCATED(d2x )) DEALLOCATE(d2x )
    IF(ALLOCATED(d2y )) DEALLOCATE(d2y )
    IF(ALLOCATED(d2z )) DEALLOCATE(d2z )

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
    IF(ALLOCATED(i1nlct     )) DEALLOCATE(i1nlct     ) 
    IF(ALLOCATED(d1rsdl     )) DEALLOCATE(d1rsdl     )
    IF(ALLOCATED(d1gsm      )) DEALLOCATE(d1gsm      )
    IF(ALLOCATED(d1grv      )) DEALLOCATE(d1grv      )
    IF(ALLOCATED(d1guv      )) DEALLOCATE(d1guv      )
    IF(ALLOCATED(d1guv_befor)) DEALLOCATE(d1guv_befor)
    IF(ALLOCATED(d1guv_after)) DEALLOCATE(d1guv_after)
    IF(ALLOCATED(d2ws       )) DEALLOCATE(d2ws       )
    IF(ALLOCATED(i2enr0 )) DEALLOCATE(i2enr0 )
    IF(ALLOCATED(d2wrks )) DEALLOCATE(d2wrks )
    IF(ALLOCATED(d2jmpm )) DEALLOCATE(d2jmpm )
    IF(ALLOCATED(d2kwns )) DEALLOCATE(d2kwns )
    IF(ALLOCATED(d2smat0)) DEALLOCATE(d2smat0)
    IF(ALLOCATED(d1svec0)) DEALLOCATE(d1svec0)

    IF(ALLOCATED(d4smat)) DEALLOCATE(d4smat)
    IF(ALLOCATED(d2svec)) DEALLOCATE(d2svec)

    IF(ALLOCATED(d1jm1  )) DEALLOCATE(d1jm1  )
    IF(ALLOCATED(d1jm2  )) DEALLOCATE(d1jm2  )
    IF(ALLOCATED(d1jm3  )) DEALLOCATE(d1jm3  )
    IF(ALLOCATED(d1jm4  )) DEALLOCATE(d1jm4  )
    IF(ALLOCATED(d1jm5  )) DEALLOCATE(d1jm5  )

    IF(ALLOCATED(vv0d   )) DEALLOCATE(vv0d   )
    IF(ALLOCATED(vv2d   )) DEALLOCATE(vv2d   )

    RETURN

  END SUBROUTINE T2COMM_DEALLOCATE
  
END MODULE T2COMM
