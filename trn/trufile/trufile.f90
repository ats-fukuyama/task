MODULE trufile
! ------------------------------------------------------------------------
!   substitution of experimental data into TASK/TR variables
! ------------------------------------------------------------------------

  USE trcomm, ONLY: ikind,rkind,ntum,nrum,nsum, &
       nrmax,dt,rhog,rhom,mdlxp,ndmax,ntxmax,tlmax,time_slc,ufid_bin

  USE truf0d,ONLY: idnm,idnfast,idnmaz,idnfastaz,idzeff
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_ufile

  INTEGER(ikind) :: tlcheck, tlsave
  REAL(rkind),DIMENSION(1:ntum) :: tmu_save

CONTAINS

  SUBROUTINE tr_ufile(ierr)
    USE ufinit,  ONLY: ufile_init
    USE trufsub, ONLY: tr_uf_time_slice
    USE truf0d,  ONLY: tr_uf0d,tr_ufile_0d_view,tr_uf_set_table
    USE trufcalc,ONLY: tr_uf_nicomplete
    USE trcomm,  ONLY: kuf_dir,kuf_dcg,kuf_dev,nrmax,ntmax,tmu, &
                       mdluf,mdlgmt,mdlsrc,mdlglb
    IMPLICIT NONE

    INTEGER(ikind),INTENT(OUT) :: ierr ! -1: fatal error
    INTEGER(ikind) :: id_mesh,id_deriv, dumntx

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
       ierr = -1
       RETURN
    END IF

    tlcheck = 1 !  = 1: consistency check ON (default)
                !  = 0: consistency check OFF
    tlsave  = 0

    ntxmax = 0
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

    IF(mdlglb==6 .OR. mdlglb==7) CALL tr_ufget_global(ierr)

    IF(mdlsrc==6 .OR. mdlsrc==7) CALL tr_ufget_source(id_mesh,id_deriv,ierr)

    IF(mdlgmt==6 .OR. mdlgmt==7) CALL tr_ufget_geometry(id_mesh,id_deriv,ierr)

    CALL tr_ufget_profile(id_mesh,id_deriv,ierr)

    CALL tr_ufget_field(id_mesh,id_deriv,ierr)


    ! load time array; save 'tmu' array in reading (KFID = NE)
    tmu(1:ntum) = tmu_save(1:ntum)
    IF(mdluf==1)THEN
       ! get the initial time from the exp. data at an arbitrary moment
       CALL tr_uf_time_slice(time_slc,tmu,ntxmax,dumntx)

    ELSE IF(mdluf==2 .AND. ntxmax==1)THEN
       mdluf = 1
       IF(mdlglb == 7) mdlglb = 6
       IF(mdlsrc == 7) mdlsrc = 6
       IF(mdlgmt == 7) mdlgmt = 6

       WRITE(6,*) ! spacing
       WRITE(6,*) '## tr_ufile: experimental data has only one point time data.'
       WRITE(6,'(1X,A,I2)') '## MDLUF is replaced. MDLUF= ',mdluf
    END IF

    
    WRITE(6,*) ! spacing
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
    USE truf0d,ONLY: tr_uf0dget_global
    USE trcomm,ONLY: tmu,rru,rau,bbu,rkapu,rdltu,phiau,wthu,wtotu

    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT) :: ierr

    CHARACTER(LEN=10) :: kfid
    INTEGER(ikind)    :: errout

!    errout    = 0 !     write inquire error message to standard output
    errout    = 1 ! not write inquire error message to standard output

    ! ***  1D data  ***

    kfid = 'RGEO' ! --> rru
    CALL tr_uf1d(kfid,tmu,rru,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       CALL tr_uf0dget_global(kfid,rru)
    END IF

    kfid = 'AMIN' ! --> rau
    CALL tr_uf1d(kfid,tmu,rau,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       CALL tr_uf0dget_global(kfid,rau)
    END IF

    kfid = 'BT' ! --> bbu
    CALL tr_uf1d(kfid,tmu,bbu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       CALL tr_uf0dget_global(kfid,bbu)
    END IF

!      IF(KUFDEV.EQ.'tftr' .AND. &
!         (KUFDCG.EQ.'50862' .OR. KUFDCG.EQ.'50921'.OR.KUFDCG.EQ.'52527')) THEN
!         RKAPU(1:NTXMAX1)=1.D0

    kfid = 'KAPPA' ! --> rkapu
    CALL tr_uf1d(kfid,tmu,rkapu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       CALL tr_uf0dget_global(kfid,rkapu)
    END IF

    kfid = 'DELTA' ! --> rdltu
    CALL tr_uf1d(kfid,tmu,rdltu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       CALL tr_uf0dget_global(kfid,rdltu)
    END IF

    errout = 1

    kfid = 'PHIA' ! --> phiau
    CALL tr_uf1d(kfid,tmu,phiau,ntxmax,mdlxp,ufid_bin,errout,ierr)
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
    USE trufsub,ONLY: tr_uf1d,tr_uf2d,tr_uftl_check
    USE trcomm,ONLY: rkev,tmu,pau,pzu,rtu,rnu,rnfu,rpu,zeffru,zeffu, &
                     qpu,jtotu,bpu,wrotu
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
    kfid = 'NE'   ! --> RNU
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    rnu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1) * 1.d-20

    ! save time array
    tmu_save(1:ntum) = tmu(1:ntum)

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
       CALL tr_uf2d(kfid,tmu,zeffru,ntxmax,nrmax,rhog,rhom, &
                                mdlxp,ufid_bin,id_mesh,id_deriv,errout,ierr)
       IF(ierr == 0)THEN
          CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
       END IF

    ELSE IF(idzeff(2))THEN
       kfid = 'ZEFF' ! --> zeffu
       CALL tr_uf1d(kfid,tmu,zeffu,ntxmax,mdlxp,ufid_bin,errout,ierr)
       IF(ierr == 0)THEN
          CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
       END IF

       DO ntx = 1, ntxmax
          zeffru(ntx,1:nrmax+1) = zeffu(ntx)
       END DO
       WRITE(6,*) &
    ' ## tr_ufget_profile: set ZEFF(1d) to ZEFFR(2d) due to lacking of ZEFFR.'
    END IF

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
    
    ! correction of negative profile due to the interpolation
    FORALL(nsu=1:nsum,ntx=1:ntxmax,nr=1:nrmax+1,rtu(nsu,ntx,nr)<0.d0)
       rtu(nsu,ntx,nr) = 1.d-3
    END FORALL
    FORALL(nsu=1:nsum,ntx=1:ntxmax,nr=1:nrmax+1,rnu(nsu,ntx,nr)<0.d0)
       rnu(nsu,ntx,nr) = 1.d-17
    END FORALL

    ! **************************************************************

    errout = 1

    kfid = 'VROT'   ! --> wrotu
    CALL tr_uf2d(kfid,tmu,wrotu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

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

  SUBROUTINE tr_ufget_field(id_mesh,id_deriv,ierr)
! ------------------------------------------------------------------------
!
! ------------------------------------------------------------------------
    USE trufsub,ONLY: tr_uf1d,tr_uf2d,tr_uftl_check
    USE truf0d,ONLY: tr_uf0dget_global
    USE trcomm,ONLY: tmu,ripu,qpu,jtotu,bpu,jnbu,jecu,jicu,jlhu,jbsu, &
                     mdlijq,mdltr_jbs
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN)  :: id_mesh,id_deriv
    INTEGER(ikind),INTENT(OUT) :: ierr

    CHARACTER(LEN=10)                :: kfid
    INTEGER(ikind)                   :: errout
    REAL(rkind),DIMENSION(ntum,nrum) :: f2out

    errout = 0 ! write inquire error message to standard output

    kfid = 'Q'   ! --> qpu
    CALL tr_uf2d(kfid,tmu,qpu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       IF(MOD(mdlijq,2)==0)THEN
          mdlijq = mdlijq - 1
          WRITE(6,*) 'XX tr_ufget_field: Q profile is not exist.'
          WRITE(6,*) '## MDLIJQ is replaced. MDLIJQ= ', mdlijq
       END IF
    END IF

    kfid = 'CURTOT'   ! --> jtotu
    CALL tr_uf2d(kfid,tmu,jtotu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       IF(MOD(mdlijq,2)==1)THEN
          mdlijq = mdlijq + 1
          WRITE(6,*) 'XX tr_ufget_field: CURTOT profile is not exist.'
          WRITE(6,*) '## MDLIJQ is replaced. MDLIJQ= ', mdlijq
       END IF       
    END IF

    kfid = 'IP' ! --> ripu
    CALL tr_uf1d(kfid,tmu,ripu,ntxmax,mdlxp,ufid_bin,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    ELSE
       IF(ntxmax==1)THEN
          CALL tr_uf0dget_global(kfid,ripu)
       ELSE IF(mdlijq <= 2)THEN
          mdlijq = mdlijq + 2
          WRITE(6,*) 'XX tr_ufget_field: IP evolution data is not exist.'
          WRITE(6,*) '## MDLIJQ is replaced. MDLIJQ= ', mdlijq
       END IF
    END IF


    errout = 1 ! following variables is not essential

    kfid = 'BPOL'   ! --> bpu
    CALL tr_uf2d(kfid,tmu,bpu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF


    ! driven current density
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
    ELSE IF(mdltr_jbs==9)THEN
          WRITE(6,*) 'XX tr_ufget_field: bootstrap current data does not exist.'
          mdltr_jbs = 2 ! default: Sauter model
          WRITE(6,'(1X,A39,I3)') 'XX  "mdltr_jbs" is reset to mdltr_jbs =',mdltr_jbs
    END IF

    RETURN
  END SUBROUTINE tr_ufget_field

! ************************************************************************

  SUBROUTINE tr_ufget_source(id_mesh,id_deriv,ierr)
! ------------------------------------------------------------------------
!   *** acquire the source and sink variables from experimental data ***
!   ***  and substitute them into TASK/TR variables                  ***
! ------------------------------------------------------------------------
    USE trufsub,ONLY: tr_uf1d,tr_uf2d,tr_uftl_check
    USE trcomm,ONLY: tmu,pnbu,pecu,pibwu,picu,plhu,pohmu,pradu,         &
         qnbu,qecu,qibwu,qicu,qlhu,qfusu,qohmu,qradu,qwallu,snbu,swallu
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


    kfid = 'QWALLE'   ! --> qwallu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qwallu(1,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)

    kfid = 'QWALLI'   ! --> qwallu
    CALL tr_uf2d(kfid,tmu,f2out,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF
    qwallu(2,1:ntxmax,1:nrmax+1) = f2out(1:ntxmax,1:nrmax+1)


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

  SUBROUTINE tr_ufget_geometry(id_mesh,id_deriv,ierr)
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

    REAL(rkind),DIMENSION(1:nrum) :: temp
    REAL(rkind) :: deriv4,deriv3
    INTEGER(ikind) :: ntx

    errout = 0 ! write inquire error message to standard output

    ! ***  Geometric factors  ***

    kfid = 'RMAJOR'   ! --> rmjrhou
    CALL tr_uf2d(kfid,tmu,rmjrhou,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'RMINOR'   ! --> rmnrhou
    CALL tr_uf2d(kfid,tmu,rmnrhou,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'KAPPAR'   ! --> rkprhou
    CALL tr_uf2d(kfid,tmu,rkprhou,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'GRHO1'   ! --> ar1rhou
    CALL tr_uf2d(kfid,tmu,ar1rhou,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'GRHO2'   ! --> ar2rhou
    CALL tr_uf2d(kfid,tmu,ar2rhou,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'VOLUME'   ! --> pvolu
    CALL tr_uf2d(kfid,tmu,pvolu,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    kfid = 'SURF'   ! --> psuru
    CALL tr_uf2d(kfid,tmu,psuru,ntxmax,nrmax,rhog,rhom,mdlxp,ufid_bin,   &
                 id_mesh,id_deriv,errout,ierr)
    IF(ierr == 0)THEN
       CALL tr_uftl_check(kfid,tmu,ntxmax,tlmax,dt,tlcheck,tlsave)
    END IF

    ! asscociated values
    ! d V/d rho
!!$    dvrhou(1:ntxmax,2:nrmax+1) = psuru(1:ntxmax,2:nrmax+1) &
!!$                              /ar1rhou(1:ntxmax,2:nrmax+1)
!!$    dvrhou(1:ntxmax,1) = 0.d0 ! at the magnetic axis

    DO ntx = 1, ntxmax
       DO nr = 1,nrmax+1
          temp(1:nrum) = pvolu(ntx,1:nrum)
          dvrhou(ntx,nr) = deriv3(nr,rhog(0:nrmax),temp(1:nrmax+1),nrmax+1,1)
       END DO
       dvrhou(ntx,1) = 0.d0 ! at the magnetic axis
    END DO

    RETURN
  END SUBROUTINE tr_ufget_geometry

END MODULE trufile
