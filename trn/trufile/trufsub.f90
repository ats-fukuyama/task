MODULE trufsub
! -----------------------------------------------------------------------
! Interface for reading and interpolating data from experimental database
!  (Selecting MDSplus or UFILE)
! -----------------------------------------------------------------------

  USE trcomm,ONLY: ikind, rkind, ntum, nrum
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf1d,tr_uf2d,tr_uftl_check,tr_uf_time_slice

CONTAINS

  SUBROUTINE tr_uf1d(kfid,tl,fout,ntxmax,mdlxp,id_bin,time_slc,id_time, &
                     errout,ierr)
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
!  time_slc : designated time for slicing at a certain time
!
!   id_time : steady state simulaion or time evolution simulation
!       = 1 : steady state (fout(1:ntxmax) have all the same functional values.)
!       = 2 : time evolution
!
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
    INTEGER(ikind),INTENT(IN)    :: mdlxp,id_bin,errout,id_time
    REAL(rkind),   INTENT(INOUT) :: time_slc

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

    ! ***   steady state simulation   ***
    IF(id_time == 1) THEN
       CALL tr_uf_time_slice(time_slc,tl,ntxmax,ntsl)
       fout(1:ntxmax) = f1(ntsl)
    ELSE IF(id_time == 2) THEN
       fout(1:ntxmax) = f1(1:ntxmax)
    END IF

    RETURN
  END SUBROUTINE tr_uf1d

! *************************************************************************

  SUBROUTINE tr_uf2d(kfid,tl,fout,ntxmax,nrmax,rhog,rhom, &
               mdlxp,id_bin,time_slc,id_time,id_mesh,id_deriv,errout,ierr)
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
!  time_slc : designated time for slicing at a certain time
!    errout : Switch for writing inquire error message to standard outoput
!              = 0: write, = 1(else): suppress
!
!   id_time : steady state simulaion or time evolution simulation
!       = 1 : steady state (fout(1:ntxmax) have all the same functional values.)
!       = 2 : time evolution
!
!   +++ for interpolation;
!   id_mesh : mesh selector in spline interpolation of radial profile
!       = 0 : integer mesh
!       = 1 : half integer mesh
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
    INTEGER(ikind),INTENT(IN) :: id_time,id_mesh,id_deriv
    REAL(rkind),DIMENSION(1:nrmax+1),       INTENT(IN)    :: rhog,rhom
    REAL(rkind),                          INTENT(INOUT) :: time_slc

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

    ! ***   steady state simulation   ***
    IF(id_time == 1) THEN
       CALL tr_uf_time_slice(time_slc,tl,ntxmax,ntsl)
       f1(1:nrlmax) = f2(ntsl,1:nrlmax)

       CALL tr_uf2d_interpolate(kfid,rl,f1,fint,u,deriv,nrlmax,nrmax, &
                                 rhog,rhom,id,id_mesh,ierr)
       DO ntx = 1, ntxmax
          fout(ntx,1:nrmax+1) = fint(1:nrmax+1)
       END DO

    ! ***   time evolution simulation   ***
    ELSE IF(id_time == 2) THEN
       
       DO ntx = 1, ntxmax
          f1(1:nrlmax) = f2(ntx,1:nrlmax)
          CALL tr_uf2d_interpolate(kfid,rl,f1,fint,u,deriv,nrlmax,nrmax, &
                                   rhog,rhom,id,id_mesh,ierr)
          fout(ntx,1:nrmax+1) = fint(1:nrmax+1)
       END DO 
    END IF

    RETURN
  END SUBROUTINE tr_uf2d

! *************************************************************************

  SUBROUTINE tr_uftl_check(kfid,tl,ntxmax,tlmax,ntlmax,dt, &
                           tl_check,tl_save)
! -------------------------------------------------------------------------
!   Consistency check of UFILE data
!    check the number of time point data
! -------------------------------------------------------------------------

    IMPLICIT NONE
    CHARACTER(10),INTENT(IN) :: kfid
    INTEGER(ikind),INTENT(IN)    :: ntxmax, tl_check
    REAL(rkind),                INTENT(IN)    :: dt
    REAL(rkind),DIMENSION(ntum),INTENT(IN)    :: tl
    REAL(rkind),                INTENT(INOUT) :: tlmax
    INTEGER(ikind),INTENT(INOUT) :: tl_save

    INTEGER(ikind),INTENT(OUT)   :: ntlmax

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
    ntlmax = INT(DINT(tl(ntxmax)*1.d2)*1.d-2/dt)

    tl_save = 1

    RETURN
  END SUBROUTINE tr_uftl_check

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_uf_time_slice(time_slc,tl,ntxmax,ntsl)
! ------------------------------------------------------------------------
!   acquire the profile of a 2D variable at the certain time point
! ------------------------------------------------------------------------

    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN)  :: ntxmax
    INTEGER(ikind),INTENT(OUT) :: ntsl
    REAL(rkind),                INTENT(INOUT) :: time_slc
    REAL(rkind),DIMENSION(ntum),INTENT(IN)    :: tl

    INTEGER(ikind) :: ioerr, ntx, ntx_min
    REAL(rkind) :: tl_min,tl_min_old


    IF(ntxmax.NE.1)THEN
       DO
          IF(time_slc.LT.tl(1) .OR. time_slc.GT.tl(ntxmax))THEN
             WRITE(6,'(A,F9.5,A,F9.5,A)')             &
             &  '# Input arbitrary time between: ',   &
             &   tl(1),' sec. - ',tl(ntxmax),' sec.'
             READ(5,*,IOSTAT=ioerr) time_slc
             IF(ioerr.NE.0 .OR. time_slc.LT.tl(1)) CYCLE
          END IF
          EXIT
       END DO

       IF(time_slc.GT.tl(ntxmax))THEN
          time_slc = tl(ntxmax)
          WRITE(6,'(A,F9.5,A,A,F9.5,A)')                    &
          &    ' Designated time: ',time_slc,' sec.',       &
          &    ' has been replaced by ',tl(ntxmax),' sec.'
          ntsl = ntxmax
          RETURN
       ELSE IF(time_slc==tl(ntxmax))THEN
          ntsl = ntxmax
          RETURN
       END IF

       tl_min = tl(ntxmax)
       DO ntx = 1, ntxmax
          IF(ABS(tl(ntx)-time_slc) .LE. 1.d-5)THEN
             ntsl = ntx
             EXIT
          END IF

          tl_min_old = tl_min
          tl_min     = MIN(ABS(tl(ntx)-time_slc), tl_min)
          IF(tl_min_old==tl_min .OR. ntx==ntxmax)THEN
             ntx_min  = ntx - 1
             IF(ntx==ntxmax) ntx_min = ntxmax
             WRITE(6,'(A,F9.5,A,A,F9.5,A)')                    &
             &    ' Designated time: ',time_slc,' sec.',       &
             &    ' has been replaced by ',tl(ntx_min),' sec.'

             time_slc = tl(ntx_min)
             ntsl     = ntx_min
             EXIT
          END IF             
       END DO
    END IF
       
    RETURN
  END SUBROUTINE tr_uf_time_slice

! *************************************************************************

  SUBROUTINE tr_uf2d_interpolate(kfid,rl,f1,fint,u,deriv,nrlmax,nrmax, &
                                 rhog,rhom,id,id_mesh,ierr)
    !  UFILE mesh --> TASK/TR mesh
    
    IMPLICIT NONE
    CHARACTER(10),                 INTENT(IN) :: kfid
    INTEGER(ikind),                INTENT(IN) :: nrlmax,nrmax, id,id_mesh
    REAL(rkind),DIMENSION(nrum),   INTENT(IN) :: rl, f1
    REAL(rkind),DIMENSION(1:nrmax+1),INTENT(IN) :: rhog,rhom
    REAL(rkind),DIMENSION(nrum),INTENT(INOUT) :: deriv

    REAL(rkind),DIMENSION(1:nrmax+1),INTENT(OUT) :: fint
    REAL(rkind),DIMENSION(4,nrum), INTENT(OUT) :: u
    INTEGER(ikind),                INTENT(OUT) :: ierr

    INTEGER(ikind) :: nr
    REAL(rkind)    :: rsl, f0

    CALL SPL1D(rl,f1,deriv,u,nrlmax,id,ierr)
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
       CALL SPL1DF(rsl,f0,rl,u,nrlmax,ierr)
       IF(ierr.NE.0)THEN
          WRITE(6,*) 'XX tr_uf2d_interpolate: SPL1DF ERROR. IERR= ',ierr
          WRITE(6,*) 'XX KFID= ',KFID, 'NR= ',NR
       END IF
       
       fint(nr) = f0
    END DO

    RETURN
  END SUBROUTINE tr_uf2d_interpolate

! *************************************************************************

!!$  SUBROUTINE tr_uf2d_edge(kfid,rhogmax,rhoa,rl,u,nrlmax,pv,pva,id_rhoa)
!!$    ! edge extrapolation
!!$    !
!!$    ! CAUTION: This routine must be called
!!$    !                       after 'tr_uf2d_interpolation' is called.
!!$
!!$    IMPLICIT NONE
!!$    CHARACTER(10),INTENT(IN) :: kfid
!!$    INTEGER(ikind),INTENT(IN) :: nrlmax,id_rhoa
!!$    REAL(rkind),                  INTENT(IN)  :: rhoa,rhogmax
!!$    REAL(rkind),DIMENSION(nrum),  INTENT(IN)  :: rl
!!$    REAL(rkind),DIMENSION(4,nrum),INTENT(IN)  :: u
!!$    REAL(rkind),                  INTENT(OUT) :: pv, pva
!!$
!!$    INTEGER(ikind) :: ierr
!!$    REAL(rkind)    :: f0
!!$
!!$    ierr = 0
!!$
!!$    pv  = 0.d0
!!$    pva = 0.d0
!!$
!!$    ! edge value
!!$    IF(id_rhoa.EQ.0) THEN
!!$       RETURN
!!$
!!$    ELSE IF(id_rhoa.EQ.1) THEN
!!$       CALL SPL1DF(rhogmax,f0,rl,u,nrlmax,ierr)
!!$       IF(ierr == 0) pv  = f0
!!$
!!$    ELSE IF(id_rhoa.EQ.2 .AND. rhoa.NE.1.d0) THEN
!!$       CALL SPL1DF(rhoa,f0,rl,u,nrlmax,ierr)
!!$       IF(ierr == 0) pva = f0
!!$    END IF
!!$
!!$    IF(ierr.NE.0)THEN
!!$       WRITE(6,*) 'XX tr_uf_edge: SPL1DF ERROR. IERR= ',ierr
!!$       WRITE(6,*) 'XX KFID= ',KFID, 'NR --> edge value '
!!$    END IF
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_uf2d_edge

END MODULE trufsub
