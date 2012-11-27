MODULE trufile
! ------------------------------------------------------------------------
!   substitution of experimental data into TASK/TR variables
! ------------------------------------------------------------------------

  USE trcomm, ONLY: ikind,rkind,ntum,nrum,nsum, &
       nrmax,dt,rhog,rhom,mdlxp,ndmax,ntxmax,ntlmax,tlmax,time_slc,ufid_bin

  USE truf0d,ONLY: idnm,idnfast,idnmaz,idnfastaz,idzeff
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_ufile

  INTEGER(ikind) :: tlcheck, tlsave
  REAL(rkind),DIMENSION(1:nrum) :: tmu_save

CONTAINS

  SUBROUTINE tr_ufile(ierr)
    USE ufinit, ONLY: ufile_init
    USE truf0d, ONLY: tr_uf0d,tr_ufile_0d_view,tr_uf_set_table
    USE trufcalc,ONLY: tr_uf_nicomplete
    USE trcomm, ONLY: kuf_dir,kuf_dcg,kuf_dev,nrmax,ntmax,mdluf,tmu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind) :: id_mesh,id_deriv

    ! locating experimental data directory and initialization (TASK/lib)
    IF(mdlxp == 0)THEN
       CALL ufile_init(kuf_dir,kuf_dev,kuf_dcg,ierr)
    ELSE IF(mdlxp == 1)THEN
       CALL IPDB_OPEN(kuf_dev,kuf_dcg)
    ELSE
       WRITE(6,*) 'XX "MDLXP" must be 0(UFILE) or 1(MDSplus). MDLXP= ',mdlxp
       RETURN
    END IF

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX Preparation for reading UFILE failed.'
       WRITE(6,*) 'XX Stop reading UFILEs. IERR= ',ierr
       RETURN
    END IF

    tlcheck = 1 !  = 1: consistency check ON (default)
                !  = 0: consistency check OFF
    tlsave  = 0

    ntxmax = 0
    ntlmax = 0
    tlmax  = 0.d0

    id_mesh  = 0 ! integer mesh
    id_deriv = 1 ! the derivative at the axis: deriv(0) = 0


    SELECT CASE(mdluf)
    CASE(1:3)
       CONTINUE
    CASE(4)
       ! read UFILE data from TOPICS
       WRITE(6,*) 'XX tr_ufile: Unsuppoted flag, for now. MDLUF= ',mdluf
       WRITE(6,*) 'XX skipped to read experimental data.'
       RETURN

    CASE DEFAULT
       WRITE(6,*) 'XX tr_ufile: Unsuppoted flag. MDLUF= ',mdluf
       WRITE(6,*) 'XX skipped to read experimental data.'
       RETURN       
    END SELECT


    CALL tr_uf0d(mdlxp,ndmax,ierr)

    CALL tr_ufget_global(ierr)

    CALL tr_ufget_profile(id_mesh,id_deriv,ierr)

    CALL tr_ufget_source(id_mesh,id_deriv,ierr)

    CALL tr_ufget_geometric(id_mesh,id_deriv,ierr)

    ! load time array
    tmu(1:ntum) = tmu_save(1:ntum)


    CALL tr_uf_nicomplete

    CALL tr_ufile_0d_view
    CALL tr_uf_set_table

    ! read error handling
    ! substitution, interpolate or extrapolate, and so on.

    WRITE(6,*) ! spacing

    RETURN
  END SUBROUTINE tr_ufile

! ************************************************************************
! ************************************************************************

  SUBROUTINE tr_ufget_global(ierr)
! ------------------------------------------------------------------------
!   *** acquire the global variables from experimental data ***
!   ***  and substitute them into TASK/TR variables         ***
! ------------------------------------------------------------------------
    USE trufsub,ONLY: tr_uf1d,tr_uftl_check
    USE trcomm,ONLY: tmu,rru,rau,ripu,bbu,rkapu,rdltu,phiau,zeffu,wthu,wtotu

    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT) :: ierr

    CHARACTER(LEN=10) :: kfid
    INTEGER(ikind)    :: errout

    errout = 0 ! write inquire error message to standard output

    ! ***  1D data  ***

    kfid = 'RGEO' ! --> rru
    CALL tr_uf1d(kfid,tmu,rru,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    ! save time array
    tmu_save(1:ntum) = tmu(1:ntum)

    kfid = 'AMIN' ! --> rau
    CALL tr_uf1d(kfid,tmu,rau,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'IP' ! --> ripu
    CALL tr_uf1d(kfid,tmu,ripu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'BT' ! --> bbu
    CALL tr_uf1d(kfid,tmu,bbu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

!      IF(KUFDEV.EQ.'tftr' .AND. &
!         (KUFDCG.EQ.'50862' .OR. KUFDCG.EQ.'50921'.OR.KUFDCG.EQ.'52527')) THEN
!         RKAPU(1:NTXMAX1)=1.D0

    kfid = 'KAPPA' ! --> rkapu
    CALL tr_uf1d(kfid,tmu,rkapu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'DELTA' ! --> rdltu
    CALL tr_uf1d(kfid,tmu,rdltu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'PHIA' ! --> phiau
    CALL tr_uf1d(kfid,tmu,phiau,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'ZEFF' ! --> zeffu
    CALL tr_uf1d(kfid,tmu,zeffu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'WTH' ! --> wthu
    CALL tr_uf1d(kfid,tmu,wthu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'WTOT' ! --> wtotu
    CALL tr_uf1d(kfid,tmu,wtotu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    RETURN
  END SUBROUTINE tr_ufget_global

! ************************************************************************

  SUBROUTINE tr_ufget_profile(id_mesh,id_deriv,ierr)
! ------------------------------------------------------------------------
!
! -----------------------------------------------------------------------
    USE trufsub,ONLY: tr_uf2d,tr_uftl_check
    USE trcomm,ONLY: rkev,tmu,pau,pzu,rtu,rnu,rnfu,rpu,zeffru,zeffu,qpu,bpu,wrotu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: id_mesh,id_deriv
    INTEGER(ikind),INTENT(OUT) :: ierr

    CHARACTER(LEN=10)                :: kfid
    CHARACTER(LEN=1)                 :: knum
    INTEGER(ikind)                   :: nsu, nsi, num, ntx, nr, errout
    REAL(rkind),DIMENSION(ntum,nrum) :: f2out

    errout = 0 ! write inquire error message to standard output

    rtu(1:nsum,1:ntum,1:nrum) = 0.d0
    rnu(1:nsum,1:ntum,1:nrum) = 0.d0
    zeffru(1:ntum,1:nrum)     = 0.d0

    ! ***  2D data  ***
    kfid = 'TE'   ! --> RTU
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    rtu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1) * 1.d-3

    kfid = 'TI'   ! --> RTU
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    DO nsi = 2, nsum !
       rtu(nsi,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1) * 1.d-3
    END DO

    kfid = 'NE'   ! --> RNU
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    rnu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1) * 1.d-20

    ! *** density of particles including impurity ***
    DO nsi = 2, nsum
       WRITE(knum,'(I1)') nsi-1

       IF(idnfast(nsi))THEN
          kfid = 'NFAST'//knum   ! --> RNU
          CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom, &
                                mdlxp,ufid_bin,id_mesh,id_deriv,errout,ierr)
          IF(ierr == 0)THEN
             CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
          END IF
          rnfu(nsi,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1) * 1.d-20
       END IF

       IF(idnm(nsi))THEN
          kfid = 'NM'//knum   ! --> RNU
          CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom, &
                                mdlxp,ufid_bin,id_mesh,id_deriv,errout,ierr)
          IF(ierr == 0)THEN
             CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
          END IF          
          rnu(nsi,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1) * 1.d-20

       END IF
    END DO

    IF(idzeff(1))THEN
       kfid = 'ZEFFR'   ! --> zeffru
       CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom, &
                                mdlxp,ufid_bin,id_mesh,id_deriv,errout,ierr)
       IF(ierr == 0)THEN
          CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
       END IF
       zeffru(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    ELSE IF(idzeff(2))THEN
       DO ntx = 1, ntxmax
          zeffru(ntx,1:nrmax+1) = zeffu(ntx)
       END DO
       WRITE(6,*) '## tr_ufget_profile: set the flat profile of ZEFF(1d) to ZEFFR(2d) due to lacking of ZEFFR data.'
    END IF
    
    ! correction of negative profile due to the interpolation
    FORALL(nsu=1:nsum,ntx=1:ntxmax,nr=1:nrmax+1,rtu(nsu,ntx,nr)<0.d0)
       rtu(nsu,ntx,nr) = 1.d-3
    END FORALL
    FORALL(nsu=1:nsum,ntx=1:ntxmax,nr=1:nrmax+1,rnu(nsu,ntx,nr)<0.d0)
       rnu(nsu,ntx,nr) = 1.d-17
    END FORALL

    ! **************************************************************

    kfid = 'Q'   ! --> qpu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qpu(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    errout = 1 ! following variables is not essential

    kfid = 'BPOL'   ! --> bpu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    bpu(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'VROT'   ! --> wrotu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    wrotu(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    ! associated values
    DO nsu = 1, nsum
       DO ntx = 1, ntxmax
          rpu(nsu,ntx,1:nrmax+1) = rnu(nsu,ntx,1:nrmax+1)*1.d20  &
                                  *rtu(nsu,ntx,1:nrmax+1)*rkev
       END DO
    END DO


    RETURN
  END SUBROUTINE tr_ufget_profile

! ************************************************************************

  SUBROUTINE tr_ufget_source(id_mesh,id_deriv,ierr)
! ------------------------------------------------------------------------
!   *** acquire the source and sink variables from experimental data ***
!   ***  and substitute them into TASK/TR variables                  ***
! ------------------------------------------------------------------------
    USE trufsub,ONLY: tr_uf1d,tr_uf2d,tr_uftl_check
    USE trcomm,ONLY: tmu,pnbu,pecu,pibwu,picu,plhu,pohmu,pradu,    &
         jtotu,jnbu,jecu,jicu,jlhu,jbsu,qnbu,qecu,qibwu,qicu,qlhu, &
         qfusu,qohmu,qradu,snbu,swallu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: id_mesh,id_deriv
    INTEGER(ikind),INTENT(OUT) :: ierr

    CHARACTER(LEN=10)                :: kfid
    INTEGER(ikind)                   :: errout
    REAL(rkind),DIMENSION(ntum,nrum) :: f2out

    errout = 1 ! NOT write inquire error message to standard output

    ! ***  1d data  ***

    kfid = 'PNBI'   ! --> pnbu
    CALL tr_uf1d(kfid,tmu,pnbu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'PECH'   ! --> pecu
    CALL tr_uf1d(kfid,tmu,pecu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'PIBW'   ! --> pibwu
    CALL tr_uf1d(kfid,tmu,pibwu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'PICRH'  ! --> picu
    CALL tr_uf1d(kfid,tmu,picu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'PLH'    ! --> plhu
    CALL tr_uf1d(kfid,tmu,plhu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'POHM'   ! --> pohmu
    CALL tr_uf1d(kfid,tmu,pohmu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'PRAD'   ! --> pradu
    CALL tr_uf1d(kfid,tmu,pradu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF


    ! ***  2D data  ***
    ! current density
    kfid = 'CURTOT'   ! --> jtotu
    CALL tr_uf2d(kfid,tmu,jtotu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'CURNBI'   ! --> jnbu
    CALL tr_uf2d(kfid,tmu,jnbu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'CURECH'   ! --> jecu
    CALL tr_uf2d(kfid,tmu,jecu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'CURICRH'   ! --> jicu
    CALL tr_uf2d(kfid,tmu,jicu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'CURLH'   ! --> jlhu
    CALL tr_uf2d(kfid,tmu,jlhu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'CURBS'   ! --> jbsu
    CALL tr_uf2d(kfid,tmu,jbsu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF


    ! heating density
    kfid = 'QNBIE'   ! --> qnbu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qnbu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QNBII'   ! --> qnbu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qnbu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'QECHE'   ! --> qecu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qecu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QECHI'   ! --> qecu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qecu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'QIBWHE'   ! --> qibwu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qibwu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QIBWHI'   ! --> qibwu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qibwu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'QICRHE'   ! --> qicu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qicu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QICRHI'   ! --> qicu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qicu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'QLHE'   ! --> qecu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qlhu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QLHI'   ! --> qecu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qlhu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'QFUSIE'   ! --> qfusu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qfusu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QFUSII'   ! --> qfusu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qfusu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'QOHM'   ! --> qohmu
    CALL tr_uf2d(kfid,tmu,qohmu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'QRAD'   ! --> qradu
    CALL tr_uf2d(kfid,tmu,qradu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF


    kfid = 'SNBIE'   ! --> snbu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    snbu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'SNBII'   ! --> snbu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    snbu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    kfid = 'SWALL'   ! --> swallu
    CALL tr_uf2d(kfid,swallu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF


    RETURN
  END SUBROUTINE tr_ufget_source

! ************************************************************************

  SUBROUTINE tr_ufget_geometric(id_mesh,id_deriv,ierr)
! ------------------------------------------------------------------------
! *** acquire the geometric quantity variables from experimental data ***
! ***  and substitute them into TASK/TR variables                     ***
! ------------------------------------------------------------------------
    USE trufsub,ONLY: tr_uf2d,tr_uftl_check
    USE trcomm,ONLY: tmu,pvolu,psuru,rmjrhou,rmnrhou,ar1rhou,ar2rhou, &
                     rkprhou,dvrhou, rru,bbu
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: id_mesh,id_deriv
    INTEGER(ikind),INTENT(OUT) :: ierr

    CHARACTER(LEN=10)                :: kfid
    INTEGER(ikind)                   :: nr, errout
    REAL(rkind),DIMENSION(ntum,nrum) :: f2out

    errout = 0 ! write inquire error message to standard output

    ! ***  Geometric factors  ***

    kfid = 'RMAJOR'   ! --> rmjrhou
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    rmjrhou(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'RMINOR'   ! --> rmnrhou
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    rmnrhou(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'KAPPAR'   ! --> rkprhou
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    rkprhou(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'GRHO1'   ! --> ar1rhou
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    ar1rhou(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'GRHO2'   ! --> ar2rhou
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    ar2rhou(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'VOLUME'   ! --> pvolu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    pvolu(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'SURF'   ! --> psuru
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    psuru(1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


    ! asscociated values
    ! d V/d rho
    dvrhou(1:ntxmax,2:nrmax+1) = psuru(1:ntxmax,2:nrmax+1) &
                              /ar1rhou(1:ntxmax,2:nrmax+1)
    dvrhou(1:ntxmax,1) = 0.d0 ! at the magnetic axis

    RETURN
  END SUBROUTINE tr_ufget_geometric

END MODULE trufile
