MODULE trufsub
! -----------------------------------------------------------------------
! Interface for reading and interpolating data from experimental database
!  (Selecting MDSplus or UFILE)
! -----------------------------------------------------------------------

  USE trcomm,ONLY: ikind, rkind, ntum, nrum
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf1d,tr_uf2d,tr_uftl_check

CONTAINS

  SUBROUTINE tr_uf1d(kfid,tl,fout,ntxmax,mdlxp,id_bin,errout,ierr)
! ----------------------------------------------------------------------
!   *** Reading interface of 1D experimental data ***
!
!   << CAUTION >>
!   Before calling this subroutine, 
!    'ufile_init'(TASK/lib/libufile/ufinit.f90) must have been called.
!
!   * ntum  : maximum size of time point data
!
!   input:
!     kfid  : Parameter index of database
!     mdlxp : switch variable which select UFILE or MDSplus
!    errout : Switch for writing inquire error message to standard outoput
!              = 0: write, = 1(else): suppress
!
!   output:
!     tl(ntum) : time point vector
!     f1(ntum) : value vector corresponding to time point vector
!     ntxmax : numver of time points
!     ierr   : error notifier
! ----------------------------------------------------------------------
    USE ufread,ONLY: ufread_1d_time

    ! interim way: This variable declaration are only for MDSplus library
    USE trcomm,ONLY: kuf_dev,kuf_dcg

    IMPLICIT NONE
    
    CHARACTER(10), INTENT(IN)    :: kfid
    INTEGER(ikind),INTENT(IN)    :: mdlxp,id_bin,errout

    REAL(rkind),DIMENSION(ntum),INTENT(OUT) :: tl,fout
    INTEGER(ikind),             INTENT(OUT) :: ntxmax, ierr

    REAL(rkind),DIMENSION(ntum) :: f1
    INTEGER(ikind) :: ntsl

    fout(1:ntum) = 0.d0
    
    IF(mdlxp .EQ. 0)THEN
       CALL ufread_1d_time(kfid,TL,F1,ntxmax,ntum,id_bin,errout,ierr)
    ELSE
       ! TASK/lib mdsplus.f
       CALL IPDB_MDS1(kuf_dev,kuf_dcg,kfid,ntum,TL,F1,ntxmax,ierr)
    END IF

    IF(ierr /= 0) THEN
       IF(errout ==0)THEN
          WRITE(6,'(A,A10,A)') '## tr_uf1d: Failed to read "',KFID,'" file.'
       END IF
       fout(1:ntum) = 0.d0
       RETURN
    ENDIF

!   Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
    tl(2:ntxmax) = tl(2:ntxmax) - tl(1)
    tl(1)        = 0.d0

    fout(1:ntxmax) = f1(1:ntxmax)

    RETURN
  END SUBROUTINE tr_uf1d

! *************************************************************************

  SUBROUTINE tr_uf2d(kfid,tl,fout,ntxmax,nrmax,rhog,rhom, &
                             mdlxp,id_bin,id_mesh,id_deriv,errout,ierr)
! ----------------------------------------------------------------------
!   *** Reading interface of 1D experimental data ***
!
!   << CAUTION >>
!   Before calling this subroutine, 
!    'ufile_init'(TASK/lib/libufile/ufinit.f90) must have been called.
!
!   * ntum  : maximum size of time point data
!   * nrum  : maximum size of radial point data
!
!  < input >
!     kfid  : Parameter index of database
!     nrmax : index number of node at the edge(separatrix)
!     rhog  : radial mesh points (integer)
!     rhom  : radial mesh points (half integer)
!     mdlxp : switch variable which select UFILE or MDSplus
!    errout : Switch for writing inquire error message to standard outoput
!              = 0: write, = 1(else): suppress
!
!   +++ for interpolation;
!   id_mesh : mesh selector in spline interpolation of radial profile
!      = 0  : integer mesh
!      = 1  : half integer mesh
!
!  id_deriv : assumption selector for spline interpolation
!       = 1 : derivative deriv(0) = 0
!       = 2 : derivative deriv(0) = 0 and deriv(nrmax) = 0
!   +++
!
!  < output >
!     tl(ntum)          : time point vector
!     fout(ntum,nrum)   : value vector corresponding to time point vector
!     ntxmax : numver of time points
!     ierr   : error notifier
! ----------------------------------------------------------------------
    USE ufread,ONLY: ufread_2d_time

    ! interim way: This variable declaration are only for MDSplus library
    USE trcomm,ONLY: kuf_dev,kuf_dcg

    IMPLICIT NONE
    CHARACTER(10), INTENT(IN) :: kfid
    INTEGER(ikind),INTENT(IN) :: nrmax,mdlxp,id_bin,errout
    INTEGER(ikind),INTENT(IN) :: id_mesh,id_deriv
    REAL(rkind),DIMENSION(1:nrmax+1),       INTENT(IN)    :: rhog,rhom

    INTEGER(ikind),                       INTENT(OUT)   :: ntxmax,ierr
    REAL(rkind),DIMENSION(ntum),          INTENT(OUT)   :: tl
    REAL(rkind),DIMENSION(ntum,nrum),     INTENT(OUT)   :: fout

    INTEGER(ikind) :: nrlmax,ntx,ntsl,id
    REAL(rkind),DIMENSION(nrum)      :: deriv,rl,f1
    REAL(rkind),DIMENSION(1:nrmax+1)   :: fint
    REAL(rkind),DIMENSION(4,nrum)    :: u
    REAL(rkind),DIMENSION(ntum,nrum) :: f2

    rl(1:nrum)          = 0.d0
    tl(1:nrum)          = 0.d0
    f2(1:ntum,1:nrum)   = 0.d0
    fout(1:ntum,1:nrum) = 0.d0

    IF(mdlxp .EQ. 0) THEN
       CALL ufread_2d_time(kfid,rl,tl,f2,nrlmax,ntxmax,nrum,ntum, &
                           id_bin,errout,ierr)
    ELSE
       CALL IPDB_MDS2(kuf_dev,kuf_dcg,kfid,nrum,ntum, &
                      rl,tl,f2,nrlmax,ntxmax,ierr)
    END IF

    IF(ierr == 1) THEN
       fout(1:ntum,1:nrum) = 0.d0
       RETURN
    ELSE IF(ierr >= 2) THEN
       IF(errout == 0)THEN
          WRITE(6,'(A15,A10,A14)') '## tr_uf2d: Failed to read "',KFID,'" file.'
       END IF
       fout(1:ntum,1:nrum) = 0.d0
       RETURN
    ENDIF

! <<<  Calculate values suitable for TASK/TR mesh using spline >>>

    ! preparation for subroutine SPL1D
    IF(id_deriv == 1)THEN
       deriv(1) = 0.d0
       id       = 1
    ELSE IF(id_deriv == 2)THEN
       deriv(1)        = 0.d0
       deriv(nrlmax)   = 0.d0
       id              = 3
    END IF

    ! Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
    tl(2:ntxmax) = tl(2:ntxmax) - tl(1)
    tl(1)        = 0.d0

    DO ntx = 1, ntxmax
       f1(1:nrlmax) = f2(ntx,1:nrlmax)
       CALL tr_uf2d_interpolate(kfid,rl,f1,fint,u,deriv,nrlmax,nrmax, &
            rhog,rhom,id,id_mesh,ierr)
       fout(ntx,1:nrmax+1) = fint(1:nrmax+1)
    END DO

    RETURN
  END SUBROUTINE tr_uf2d

! *************************************************************************

  SUBROUTINE tr_uftl_check(kfid,tl,ntxmax,tlmax,dt,tl_check,tl_save)
! -------------------------------------------------------------------------
!   ***   Consistency check of UFILE data   ***
!
!   This subroutine checks the number of time point data of each variable
!    data file and the equivalence of it.
! -------------------------------------------------------------------------

    IMPLICIT NONE
    CHARACTER(10),INTENT(IN) :: kfid
    INTEGER(ikind),INTENT(IN)    :: ntxmax, tl_check
    REAL(rkind),                INTENT(IN)    :: dt
    REAL(rkind),DIMENSION(ntum),INTENT(IN)    :: tl
    REAL(rkind),                INTENT(INOUT) :: tlmax
    INTEGER(ikind),INTENT(INOUT) :: tl_save

    IF(tl_save.NE.0) THEN
       IF(tl_check == 0 .AND. tlmax == 0.D0) THEN
          CONTINUE
       ELSE
          IF(tlmax.NE.tl(ntxmax)) THEN
             WRITE(6,*) 'XX tr_uftl_check: ',KFID,'UFILE HAS AN ERROR!'
             WRITE(6,*) 'XX tlmax= ',tlmax,' tl(ntxmax)= ',tl(ntxmax)
             STOP
          END IF
       END IF
    END IF

    tlmax  = tl(ntxmax)

    tl_save = 1

    RETURN
  END SUBROUTINE tr_uftl_check

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_uf2d_interpolate(kfid,rl,f2,fint,u,deriv,nrlmax,nrmax, &
                                 rhog,rhom,id_grd,id_mesh,ierr)
! ------------------------------------------------------------------------
!  data interpolation: UFILE mesh --> TASK/TR mesh
!
!   < input >
!  kfid    : the name of variable data of the UFILE
!  rl      : the radial mesh points data of the UFILE
!  f2      : the functional value of the UFILE
!  deriv   : the derivative value for spline interpolation
!  nrlmax  : the number of radial mesh points of UFILE
!  nrmax   : the number of radial mesh points of TASK/TR
!  rhog    : the radial mesh points data of TAKS/TR (integer mesh)
!  rhom    : the radial mesh points data of TAKS/TR (half integer mesh)
!  id_grd  : the boundary condition selector for spline interpolation
!  id_mesh : the selector of mesh type (=0: integer, else=half integer)
!
!   < output >
!  fint    : the interpolated functional data
!  u       : the spline coefficients used in the interpolation
!  ierr    : error identifier
! ------------------------------------------------------------------------
    
    IMPLICIT NONE
    CHARACTER(10),                 INTENT(IN) :: kfid
    INTEGER(ikind),                INTENT(IN) :: nrlmax,nrmax,id_grd,id_mesh
    REAL(rkind),DIMENSION(nrum),   INTENT(IN) :: rl, f2
    REAL(rkind),DIMENSION(1:nrmax+1),INTENT(IN) :: rhog,rhom
    REAL(rkind),DIMENSION(nrum),INTENT(INOUT) :: deriv

    REAL(rkind),DIMENSION(1:nrmax+1),INTENT(OUT) :: fint
    REAL(rkind),DIMENSION(4,nrum), INTENT(OUT) :: u
    INTEGER(ikind),                INTENT(OUT) :: ierr

    INTEGER(ikind) :: nr
    REAL(rkind)    :: rsl, f0

    CALL SPL1D(rl,f2,deriv,u,nrlmax,id_grd,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX tr_uf2d_interpolate: SPL1D ERROR. IERR= ',IERR
       WRITE(6,*) 'XX KFID= ',kfid
    END IF

    fint(1:nrmax+1) = 0.d0
    DO nr = 1, nrmax+1
       IF(id_mesh.EQ.0)THEN
          rsl = rhog(nr)
       ELSE IF(nr.NE.nrmax+1 .AND. id_mesh.EQ.1)THEN
          rsl = rhom(nr)
       END IF
       CALL SPL1DF(rsl,f0,rl(1:nrlmax),u(1:4,1:nrlmax),nrlmax,ierr)
       IF(ierr.NE.0)THEN
          WRITE(6,*) 'XX tr_uf2d_interpolate: SPL1DF ERROR. IERR= ',ierr
          WRITE(6,*) 'XX KFID= ',KFID, 'NR= ',NR
       END IF
       
       fint(nr) = f0
    END DO

    RETURN
  END SUBROUTINE tr_uf2d_interpolate

END MODULE trufsub
