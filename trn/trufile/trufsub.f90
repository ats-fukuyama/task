MODULE trufsub

  USE trcomm,ONLY: ikind, rkind, ntum, nrum

  PRIVATE
  PUBLIC tr_uf_check_impurity,tr_uf1d,tr_uf2d,tr_uftl_check

CONTAINS

  SUBROUTINE tr_uf_check_impurity(kdirx,kuf_dcg,kuf_dev,mdni,mdslct,nmchk)
! ------------------------------------------------------------------------ 
!           CHECKING WHETHER IMPURITY EXISTS
! ------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(LEN=80),INTENT(IN) :: kdirx,kuf_dcg,kuf_dev
    INTEGER(ikind),INTENT(IN)    :: mdni
    INTEGER(ikind),INTENT(OUT)   :: mdslct,nmchk

    CHARACTER(LEN=80)            :: KDIRR2, KFILE
    CHARACTER(LEN=20)            :: KFID
    INTEGER(ikind)               :: IKDIRX,IKNDCG,IKNDEV,KL2,ierr
    LOGICAL LEX


    IKNDEV= len_trim(KUF_DEV)
    IKNDCG= len_trim(KUF_DCG)

    IKDIRX= len_trim(KDIRX)
    KDIRR2= KDIRX(1:IKDIRX)//KUF_DEV(1:IKNDEV)//'2d'//KUF_DCG(1:IKNDCG)//'.'
    KL2   = len_trim(KDIRR2)

    IF(MDNI.LT.0.OR.MDNI.GT.3) THEN
       MDNI=0
       WRITE(6,*)'Warning: tr_uf_check_impurity: '
       WRITE(6,*)'   Parameter "MDNI" is out of range, and set to zero.'
    END IF
    MDSLCT=0

    DO
       KFID='ZEFFR'
       KFILE=KDIRR2(1:KL2)//KFID
       INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=ierr)
       IF(ierr /= 0) EXIT
       IF(LEX) MDSLCT=MDSLCT+1

       KFID='NM1'
       KFILE=KDIRR2(1:KL2)//KFID
       INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=ierr)
       IF(ierr /= 0) EXIT
       IF(LEX) THEN
          NMCHK=1
          MDSLCT=MDSLCT+2
          
          KFID='NM2'
          KFILE=KDIRR2(1:KL2)//KFID
          INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=ierr)
          IF(ierr /= 0) EXIT
          IF(LEX) NMCHK=2
          
          KFID='NM3'
          KFILE=KDIRR2(1:KL2)//KFID
          INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=ierr)
          IF(ierr /= 0) EXIT
          IF(LEX) NMCHK=3
       ENDIF
       
       KFID='NIMP'
       KFILE=KDIRR2(1:KL2)//KFID
       INQUIRE(FILE=KFILE,EXIST=LEX,IOSTAT=ierr)
       IF(ierr /= 0) EXIT
       IF(LEX) MDSLCT=MDSLCT+4

       EXIT
    END DO

!

    SELECT CASE(mdslct)
    CASE(0)
    CASE(1)
    CASE(2)
    CASE(3)
    CASE(4)
    CASE(5)
    CASE(6)
    CASE(7)
    END SELECT


    RETURN
  END SUBROUTINE tr_uf_check_impurity

! **********************************************************************

  SUBROUTINE tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tl,f1,ntxmax,mdlxp,ierr)
! ntxmax : number of time step of experimental data

    IMPLICIT NONE
    
    CHARACTER(10),INTENT(IN) :: kfid
    CHARACTER(80),INTENT(IN) :: kdirx,kuf_dev,kuf_dcg

    REAL(rkind),DIMENSION(ntum),INTENT(OUT) :: tl,f1
    
    INTEGER(ikind),INTENT(IN) :: mdlxp
    INTEGER(ikind),INTENT(OUT) :: ntxmax, ierr

    
    IF(mdlxp .EQ. 0)THEN
       CALL UFREAD_TIME(kdirx,kuf_dev,kuf_dcg,kfid, &
                        TL,F1,ntxmax,ntum,ierr)
    ELSE
       ! TASK/lib mdsplus.f
       CALL IPDB_MDS1(kuf_dev,kuf_dcg,kfid,ntum,TL,F1,ntxmax,ierr)
    END IF

    IF(ierr.NE.0) THEN
       WRITE(6,'(A13,A10,A14)') '## tr_uf1d: NO "',KFID,'" FILE EXISTS.'
       RETURN
    ENDIF

!   Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
    tl(2:ntxmax) = tl(2:ntxmax) - tl(1)
    tl(1)        = 0.d0

    RETURN
  END SUBROUTINE tr_uf1d

! *************************************************************************

  SUBROUTINE tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,   &
                     pv,pva,tl,fout,ntxmax,rhoa,nrmax,mdlxp, &
                     time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)

!   id_time : steady state simulaion or time evolution simulation
!       = 1 : steady state
!       = 2 : time evolution
!
!   id_mesh : mesh selector in spline interpolation of radial profile
!       = 0 : integer mesh
!       = 1 : half integer mesh
!
!  id_deriv : assumption selector for spline interpolation
!       = 1 : derivative deriv(0) = 0
!       = 2 : derivative deriv(0) = 0 and deriv(nrmax) = 0
!
!   id_rhoa : whether peripheral values (rho = rhoa) are caluculated or not
!       = 0 : not calculated
!       = 1 : calculated ( rho = rhog(nrmax) )
!       = 2 : calculated ( rho = rhoa        )             


    IMPLICIT NONE
    CHARACTER(10),INTENT(IN) :: kfid
    CHARACTER(80),INTENT(IN) :: kdirx,kuf_dev,kuf_dcg

    INTEGER(ikind),INTENT(IN)  :: nrmax,mdlxp
    INTEGER(ikind),INTENT(IN)  :: id_time,id_mesh,id_deriv,id_rhoa
    INTEGER(ikind),INTENT(OUT) :: ntxmax,ierr

    REAL(rkind),                          INTENT(IN)    :: rhoa
    REAL(rkind),                          INTENT(INOUT) :: time_slc
    REAL(rkind),DIMENSION(0:nrmax),       INTENT(IN)    :: rhog,rhom
    REAL(rkind),DIMENSION(ntum),          INTENT(OUT)   :: tl,pv,pva
    REAL(rkind),DIMENSION(ntum,nrum),     INTENT(OUT)   :: fout

    INTEGER(ikind) :: nrlmax,ntx,ntsl,id
    REAL(rkind)    :: pvs, pvas
    REAL(rkind),DIMENSION(nrum)         :: deriv,rl,f1
    REAL(rkind),DIMENSION(0:nrmax)      :: fint
    REAL(rkind),DIMENSION(4,nrum)       :: u
    REAL(rkind),DIMENSION(ntum,nrum) :: f2

    rl(1:nrum)          = 0.d0
    tl(1:nrum)          = 0.d0
    f2(1:ntum,1:nrum)   = 0.d0
    fout(1:ntum,1:nrum) = 0.d0
    pv(1:ntum)  = 0.d0
    pva(1:ntum) = 0.d0

    IF(mdlxp .EQ. 0) THEN
       CALL UFREAD2_TIME(kdirx,kuf_dev,kuf_dcg,kfid,rl,tl,f2, &
                         nrlmax,ntxmax,nrum,ntum,ierr)
    ELSE
       CALL IPDB_MDS2(kuf_dev,kuf_dcg,kfid,nrum,ntum, &
                      rl,tl,f2,nrlmax,ntxmax,ierr)
    END IF

    IF(ierr == 1) THEN
       fout(1:ntum,1:nrum) = 0.d0
       RETURN
    ELSE IF(ierr >= 2) THEN
       fout(1:ntum,1:nrum) = 0.d0
       WRITE(6,'(A15,A10,A14)') '## tr_uf2d: NO "',KFID,'" FILE EXISTS.'
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
       CALL tr_uf2d_time_slice(time_slc,tl,ntxmax,ntsl)
       f1(1:nrlmax) = f2(ntsl,1:nrlmax)

       CALL tr_uf2d_interpolate(kfid,rl,f1,fint,deriv,nrlmax,nrmax, &
                                 rhog,rhom,id,id_mesh,ierr)
       fout(ntxmax,1:nrmax+1) = fint(0:nrmax)

       IF(id_rhoa /= 0)THEN
          CALL tr_uf2d_edge(kfid,rhog(nrmax),rhoa,rl,u,nrlmax,pvs,pvas, &
                            id_rhoa,id_mesh)
          pv(1)  = pvs
          pva(1) = pvas
       END IF

    ! ***   time evolution simulation   ***
    ELSE IF(id_time == 2) THEN
       
       DO ntx = 1, ntxmax
          f1(1:nrlmax) = f2(ntx,1:nrlmax)
          CALL tr_uf2d_interpolate(kfid,rl,f1,fint,deriv,nrlmax,nrmax, &
                                   rhog,rhom,id,id_mesh,ierr)
          fout(ntx,1:nrmax+1) = fint(0:nrmax)
          
          IF (id_rhoa /= 0)THEN
             CALL tr_uf2d_edge(kfid,rhog(nrmax),rhoa,rl,u,nrlmax,pvs,pvas, &
                               id_rhoa,id_mesh)
             pv(ntx)  = pvs
             pva(ntx) = pvas
          END IF
       END DO 

    END IF ! id_time

 RETURN
END SUBROUTINE tr_uf2d

! *************************************************************************

  SUBROUTINE tr_uftl_check(kfid,tl,ntxmax,tlmax,ntlmax,dt, &
                           uftl_check,uftl_save)
! -------------------------------------------------------------------------
!   Consistency check of UFILE data
!    check the number of time data
! -------------------------------------------------------------------------

    IMPLICIT NONE
    CHARACTER(10),INTENT(IN) :: kfid
    REAL(rkind),                INTENT(IN)    :: dt
    REAL(rkind),DIMENSION(ntum),INTENT(IN)    :: tl
    REAL(rkind),                INTENT(INOUT) :: tlmax
    INTEGER(ikind),INTENT(IN)    :: ntxmax, uftl_check
    INTEGER(ikind),INTENT(INOUT) :: uftl_save
    INTEGER(ikind),INTENT(OUT)   :: ntlmax

    IF(uftl_save.NE.0) THEN
       IF(uftl_check.GT.0 .AND. tlmax.EQ.0.D0) THEN
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

    uftl_save = 1

    RETURN
  END SUBROUTINE tr_uftl_check

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_uf2d_time_slice(time_slc,tl,ntxmax,ntsl)

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
             IF(ioerr.NE.0 .OR. &
             &      (time_slc.LT.tl(1) .OR. time_slc.GT.tl(ntxmax))) CYCLE
          END IF
          EXIT
       END DO

       tl_min = tl(ntxmax)
       DO ntx = 1, ntxmax
          IF(ABS(tl(ntx)-time_slc) .LE. 1.d-5)THEN
             ntsl = ntx
             EXIT
          END IF

          tl_min_old = tl_min
          tl_min     = MIN(ABS(tl(ntx)-time_slc), tl_min)
          IF(tl_min_old.EQ.tl_min)THEN
             ntx_min  = ntx - 1
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
  END SUBROUTINE tr_uf2d_time_slice

! *************************************************************************

  SUBROUTINE tr_uf2d_interpolate(kfid,rl,f1,fint,deriv,nrlmax,nrmax, &
                                 rhog,rhom,id,id_mesh,ierr)
!  UFILE mesh --> TASK/TR mesh
    
    IMPLICIT NONE
    CHARACTER(10),              INTENT(IN)     :: kfid
    INTEGER(ikind),             INTENT(IN)     :: nrlmax,nrmax, id,id_mesh
    INTEGER(ikind),             INTENT(OUT)    :: ierr
    REAL(rkind),DIMENSION(nrum),INTENT(IN)     :: rl, f1
    REAL(rkind),DIMENSION(nrum),INTENT(INOUT)  :: deriv
    REAL(rkind),DIMENSION(0:nrmax),INTENT(IN)  :: rhog,rhom
    REAL(rkind),DIMENSION(0:nrmax),INTENT(OUT) :: fint

    INTEGER(ikind) :: nr
    REAL(rkind)    :: rsl, f0
    REAL(rkind),DIMENSION(4,nrum) :: u


    CALL SPL1D(rl,f1,deriv,u,nrlmax,id,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX tr_uf2d_interpolate: SPL1D ERROR. IERR= ',IERR
       WRITE(6,*) 'XX KFID= ',kfid
    END IF

    fint(0:nrmax) = 0.d0
    DO nr = 0, nrmax
       IF(id_mesh.EQ.0)THEN
          rsl = rhog(nr)
       ELSE IF(nr.NE.0 .AND. id_mesh.EQ.1)THEN
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

  SUBROUTINE tr_uf2d_edge(kfid,rhogmax,rhoa,rl,u,nrlmax,pv,pva, &
                          id_rhoa,id_mesh)
    ! edge extrapolation

    IMPLICIT NONE
    CHARACTER(10),INTENT(IN) :: kfid
    INTEGER(ikind),INTENT(IN) :: nrlmax,id_rhoa,id_mesh
    REAL(rkind),                  INTENT(IN)  :: rhoa,rhogmax
    REAL(rkind),DIMENSION(nrum),  INTENT(IN)  :: rl
    REAL(rkind),DIMENSION(4,nrum),INTENT(IN)  :: u
    REAL(rkind),                  INTENT(OUT) :: pv, pva

    INTEGER(ikind) :: ierr
    REAL(rkind)    :: f0

    ierr = 0

    pv  = 0.d0
    pva = 0.d0

    ! edge value
    IF(id_rhoa.EQ.0) THEN
       RETURN

    ELSE IF(id_rhoa.EQ.1 .AND. id_mesh.EQ.1) THEN !for variables on half mesh
       CALL SPL1DF(rhogmax,f0,rl,u,nrlmax,ierr)
       IF(ierr == 0) pv  = f0

    ELSE IF(id_rhoa.EQ.2 .AND. rhoa.NE.1.d0) THEN
       CALL SPL1DF(rhoa,f0,rl,u,nrlmax,ierr)
       IF(ierr == 0) pva = f0
    END IF

    IF(ierr.NE.0)THEN
       WRITE(6,*) 'XX tr_uf_edge: SPL1DF ERROR. IERR= ',ierr
       WRITE(6,*) 'XX KFID= ',KFID, 'NR --> edge value '
    END IF

    RETURN
  END SUBROUTINE tr_uf2d_edge

END MODULE trufsub
