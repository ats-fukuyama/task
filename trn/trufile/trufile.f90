MODULE trufile

  USE trcomm, ONLY: ikind, rkind, ntum,nrum,nsum

  PRIVATE
  PUBLIC tr_ufile

CONTAINS

  SUBROUTINE tr_ufile

    USE trcomm, ONLY: kdirx,kuf_dcg,kuf_dev,mdni,mdslct,nmchk,&
         nrmax,ntmax,mdluf
    IMPLICIT NONE

    ! initialization of ufile variables

    
    CALL tr_uf_check_impurity(kdirx,kuf_dcg,kuf_dev,mdni,mdslct,nmchk)

    SELECT CASE(mdluf)
    CASE(1)
       CALL tr_ufile_steady
!!$    CASE(2)
!!$       CALL tr_ufile_time
!!$    CASE(3)
!!$       CALL tr_ufile_topics
    END SELECT

    RETURN
  END SUBROUTINE tr_ufile

! *************************************************************************

  SUBROUTINE tr_ufile_steady

    USE trcomm, ONLY: nrmax,nsamax,mdlxp,kdirx,kuf_dir,kuf_dev,kuf_dcg, &
         ntxmax,tlmax,ntlmax,dt,rhog,rhom,rhoa, &
         time_slc,tmu,rru,rau,ripu,bbu,rkapu,phiau,rtu,ptsu,rnu,pnsu, &
         nrd1

    USE trufsub, ONLY: tr_uf1d,tr_uf2d,tr_uftl_check

    IMPLICIT NONE
    REAL(rkind),DIMENSION(ntum) :: pv, pva
    REAL(rkind),DIMENSION(ntum,nrum) :: f2out
    INTEGER(ikind) :: nsu
    INTEGER(ikind) :: ierr, uftl_save, uftl_check
    INTEGER(ikind) :: id_time,id_mesh,id_deriv,id_rhoa
    CHARACTER(10)  :: kfid

    uftl_check = 1 ! default

    uftl_save = 0
    tlmax     = 0.d0

    ! locating experimental data directory (TASK/lib)
    IF(mdlxp == 0)THEN
       CALL ufile_interface(kdirx,kuf_dir,kuf_dev,kuf_dcg,0)
    ELSE IF(mdlxp == 1)THEN
       CALL IPDB_OPEN(kuf_dev,kuf_dcg)
    ELSE
       WRITE(6,*) 'XX "MDLXP" must be 0(UFILE) or 1(MDSplus). MDLXP= ',mdlxp
       RETURN
    END IF


!    CALL tr_ufile_get ??

!   **************************   1D data   **************************

    kfid = 'RGEO' ! --> RRU
    CALL tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tmu,rru,ntxmax,mdlxp,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

!    write(*,*) 'RRU =',(RRU(i), i = 1,ntxmax)
!    write(*,*) 'TMU =',(TMU(i), i = 1,ntxmax)

    kfid = 'AMIN' ! --> RAU
    CALL tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tmu,rau,ntxmax,mdlxp,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

!    write(*,*) 'RAU =',(RAU(i), i = 1,ntxmax)
!    write(*,*) 'TMU =',(TMU(i), i = 1,ntxmax)

    kfid = 'IP' ! --> RIPU
    CALL tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tmu,ripu,ntxmax,mdlxp,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'BT' ! --> BBU
    CALL tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tmu,bbu,ntxmax,mdlxp,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

!      IF(KUFDEV.EQ.'tftr' .AND. &
!         (KUFDCG.EQ.'50862' .OR. KUFDCG.EQ.'50921'.OR.KUFDCG.EQ.'52527')) THEN
!         RKAPU(1:NTXMAX1)=1.D0

    kfid = 'KAPPA' ! --> RKAPU
    CALL tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tmu,rau,ntxmax,mdlxp,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'PHIA' ! --> PHIAU
    CALL tr_uf1d(kfid,kdirx,kuf_dev,kuf_dcg,tmu,rau,ntxmax,mdlxp,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


!   **************************   2D data   *****************************
    id_time = 1 ! steady state simulation

    id_mesh  = 0 ! integer mesh
    id_deriv = 1 ! derivative deriv(0) = 0
    id_rhoa  = 1 ! for the time being
    kfid = 'TE'   ! --> RTU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    rtu(1,ntxmax,1:nrmax+1) = f2out(ntxmax,1:nrmax+1) * 1.d-3
    ptsu(1,1) = pv(1)
    ! ??? rhoa value ???
    IF(rhoa.EQ.2.d0)THEN
       CONTINUE
    END IF

    kfid = 'TI'   ! --> RTU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    DO nsu = 2, nsum ! ??? ns ?
       rtu(nsu,ntxmax,1:nrmax+1) = f2out(ntxmax,1:nrmax+1) * 1.d-3
       ptsu(nsu,1) = pv(1)
       IF(rhoa.EQ.2.d0)THEN
          CONTINUE
       END IF
    END DO


    kfid = 'NE'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    rnu(1,ntxmax,1:nrmax+1) = f2out(ntxmax,1:nrmax+1) * 1.d-20


    kfid = 'NFAST'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


    kfid = 'NM1'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


    kfid = 'NM2'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


    kfid = 'NM3'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


    kfid = 'ZEFFR'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


    kfid = 'NIMP'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'PBEAM'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'Q'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'CURTOT'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'CURNBI'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'BPOL'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'QNBIE'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'QNBII'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'QICRHE'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'PICRHI'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'QECH'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'QRAD'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'QOHM'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'VROT'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


! ***********************   Geometry factors   ************************

    kfid = 'RMAJOR'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'RMINOR'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'GRHO1'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'GRHO2'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'VOLUME'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)

    kfid = 'SURF'   ! --> RNU
    CALL tr_uf2d(kfid,kdirx,kuf_dev,kuf_dcg,rhog,rhom,         &
                 pv,pva,tmu,f2out,ntxmax,rhoa,nrmax,mdlxp,     &
                 time_slc,id_time,id_mesh,id_deriv,id_rhoa,ierr)
    CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,ntlmax,dt,uftl_check,uftl_save)


    ! read error handling
    ! substitution, interpolate or extrapolate, and so on.

    RETURN
  END SUBROUTINE tr_ufile_steady

! -------------------------------------------------------------------------

!!$  SUBROUTINE tr_ufile_time
!!$
!!$    ! locating UFILE directory (TASK/lib)
!!$    CALL ufile_interface(kdirx,kufdir,kufdev,kufdcg,0)
!!$
!!$!    CALL tr_ufile_get ??
!!$
!!$    ! substitution, interpolate or extrapolate, and so on.
!!$
!!$  END SUBROUTINE tr_ufile_time
!!$
!!$! -------------------------------------------------------------------------
!!$
!!$  SUBROUTINE tr_ufile_topics
!!$
!!$    ! locating UFILE directory (TASK/lib)
!!$    CALL ufile_interface(kdirx,kufdir,kufdev,kufdcg,0)
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_ufile_topics

! ************************************************************************
! ************************************************************************

  SUBROUTINE tr_uf_complete
!     MDNI is a parameter that can control which data is used
!     to determine bulk density, impurity density or effective
!     charge number among those data.
!     If all the data above do not exist, MDNI is set to zero
!     automatically regardless of original MDNI.
!           0 : NSMAX=2, ne=ni
!           1 : calculate nimp and ni profiles from NE, ZIMP and ZEFFR
!           2 : calculate nimp and zeff profiles from NE, ZIMP and NM1
!           3 : calculate zeff and ni profiles from NE, ZIMP and NIMP

    INTEGER(ikind),INTENT(IN) :: mdslct
    INTEGER(ikind),INTENT(INOUT) :: mdni

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

    SELECT CASE(mdni)
    CASE(0)
    CASE(1)
    CASE(2)
    CASE(3)
    END SELECT

    RETURN
  END SUBROUTINE tr_uf_complete

END MODULE trufile
