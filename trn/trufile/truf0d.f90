MODULE truf0d
  USE uflist, ONLY: uf0d
  USE trcomm, ONLY: ikind, rkind, nsum,ntum,nrum,ntxmax

  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf0d, &
         idnm,idnfast,idnmaz,idnfastaz,idzeff, toknam,shotnum

  LOGICAL,DIMENSION(1:9) :: idnm,idnfast,idnmaz,idnfastaz
  LOGICAL,DIMENSION(1:2) :: idzeff

  CHARACTER(LEN=15) :: toknam,shotnum

CONTAINS

  SUBROUTINE tr_uf0d(mdlxp,ndmax,ierr)
! ----------------------------------------------------------------------
!   This subroutine acquires the 0D data only from UFILE format for now.
!    ( MDSplus is not supportted. )
! ----------------------------------------------------------------------
    USE ufread,ONLY: ufread_0d
    USE trcomm,ONLY: zeffu

    INTEGER(ikind),INTENT(IN)  :: mdlxp
    INTEGER(ikind),INTENT(OUT) :: ndmax,ierr


    IF(mdlxp == 0)THEN
       CALL ufread_0d(ndmax,ierr)
       IF(ierr /= 0)THEN
          WRITE(6,'(A13,A10,A14)') '## tr_uf0d: Failed to read 0d.dat file.'
          ierr = 1
          RETURN
       END IF

       toknam  = uf0d(158)%fc0
       shotnum = uf0d(139)%fc0

       WRITE(6,*) '# DEVICE: ',TRIM(toknam),'     SHOT: ',TRIM(shotnum)

       CALL tr_uf_check_species

       CALL tr_uf_get_papz

       zeffu(1:ntxmax) = 0.d0
       IF(idzeff(1))  zeffu(1:ntxmax) = uf0d(175)%fr0

    ENDIF

    RETURN
  END SUBROUTINE tr_uf0d

! *********************************************************************

  SUBROUTINE tr_uf_check_species
! ------------------------------------------------------------------------ 
!   *** Checking whether impurity exists ***
!
!   This subrouitne must be called after 'tr_uf0d' is called.
! ------------------------------------------------------------------------
    USE ufinit, ONLY: ufile_inquire
    IMPLICIT NONE

    INTEGER(ikind),PARAMETER :: pamax = 100, pzmax = 50

    CHARACTER(LEN=100) :: kfile, kfileb ! dummy for 'ufile_inquire'
    CHARACTER(LEN=10)  :: kfid, kufdim
    CHARACTER(LEN=1)   :: knum
    INTEGER(ikind)     :: num, nsi, id_bin, errout, ierr

    id_bin = 2 ! confirm the existence of only ASCII file
    errout = 1 ! suppress the error message in checking the existence of UFILE

    idnm(1:9)      = .FALSE.
    idnfast(1:9)   = .FALSE.
    idnmaz(1:9)    = .FALSE.
    idnfastaz(1:9) = .FALSE.
    idzeff(1:2)    = .FALSE.

    kufdim = '2d'
    DO nsi = 1, 9
       WRITE(knum,'(I1)') nsi

       ! NMn
       kfid = 'NM'//knum
       CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr)
       IF(ierr == 0) idnm(nsi) = .TRUE.

       ! NFASTn
       kfid = 'NFAST'//knum
       CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr)
       IF(ierr == 0) idnfast(nsi) = .TRUE.
    END DO

    kfid = 'ZEFFR'
    CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr)
    IF(ierr == 0) idzeff(1) = .TRUE.


    ierr = 0

    DO nsi = 1, 9
       num = (nsi-1)*2 + 68 ! see uf0d list: NFAST1A ~ NFAST9Z
       IF(uf0d(num)%lex .EQV. .TRUE. .AND. uf0d(num+1)%lex .EQV. .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnfastaz(nsi) = .TRUE.
          END IF
       END IF

       num = (nsi-1)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
       IF(uf0d(num)%lex .EQV. .TRUE. .AND. uf0d(num+1)%lex .EQV. .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnmaz(nsi) = .TRUE.
          END IF
       END IF

       ! The species can not be idetified in the following case.
       IF(idnfast(nsi).EQV. .TRUE. .AND. idnfastaz(nsi).EQV. .FALSE.)THEN
          ierr = ierr + 10**(nsi-1)
       END IF
       IF(idnm(nsi).EQV. .TRUE. .AND. idnmaz(nsi).EQV. .FALSE.)THEN
          ierr = ierr + 2*10**(nsi-1)
       END IF
    END DO
    
    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX tr_uf_check_species: NO consistency in species data.'
       WRITE(6,*) 'XX IERR= ',ierr
       RETURN
    END IF

    ! 0d.dat: kfid = ZEFF
    IF(uf0d(175)%lex .EQV. .TRUE. .AND. &
         0 < uf0d(175)%fr0 .AND. uf0d(175)%fr0 < pzmax)THEN
       idzeff(2) = .TRUE.
    END IF

    RETURN
  END SUBROUTINE tr_uf_check_species

! **********************************************************************

  SUBROUTINE tr_uf_get_papz
! ----------------------------------------------------------------------
!   set the atomic number and the charge number of species which are
!    contained in experimental data
! ----------------------------------------------------------------------
    USE trcomm, ONLY: ame,amp,pau,pzu,pafu,pzfu

    IMPLICIT NONE
    CHARACTER(LEN=10) :: kfid
    CHARACTER(LEN=1)  :: knum
    INTEGER(ikind)    :: nsi, num
    REAL(rkind)       :: pans, pzns


    pau(1) = ame/amp
    pzu(1) = -1.d0
    
    WRITE(6,*) ! spacing
    ! bulk ions
    WRITE(6,'(A33)',ADVANCE="NO") ' ## UFILE data of ions contains: '
    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          WRITE(knum,'(I1)') nsi
          kfid = 'NM'//knum
          WRITE(6,'(A3)',ADVANCE='NO') TRIM(kfid)
          num = (nsi-1)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
          
          pans = uf0d(num)%fr0
          pzns = uf0d(num+1)%fr0
          CALL tr_uf_identify_ions(pans,pzns)
          pau(nsi) = pans
          pzu(nsi) = pzns
       END IF
    END DO

    pafu(1) = 0.d0
    pzfu(1) = 0.d0

    WRITE(6,*) ! breaking line
    ! fast ions
    WRITE(6,'(A33)',ADVANCE="NO") ' ## UFILE data of fast ions contains: '
    DO nsi = 2, nsum
       IF(idnfast(nsi))THEN
          WRITE(knum,'(I1)') nsi
          kfid = 'NFAST'//knum
          WRITE(6,'(A3)',ADVANCE='NO') TRIM(kfid)
          num = (nsi-1)*2 + 68 ! see uf0d list: NFAST1A ~ NFAST9Z
          
          pans = uf0d(num)%fr0
          pzns = uf0d(num+1)%fr0
          CALL tr_uf_identify_ions(pans,pzns)
          pafu(nsi) = pans
          pzfu(nsi) = pzns       
       END IF
    END DO
    WRITE(6,*) ! spacing
       
    RETURN
  END SUBROUTINE tr_uf_get_papz

! ************************************************************************

  SUBROUTINE tr_uf_identify_ions(pans,pzns)
! -----------------------------------------------------------------------
!   Identify the atomic number and the charge number of ions
! -----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT) :: pans, pzns

    IF(0.9d0 < pzns .AND. pzns < 1.1d0)THEN
       pzns = 1.d0 ! hydrogenic ions

       IF(0.9d0 < pans .AND. pans < 1.1d0)THEN
          pans = 1.d0 ! H
          WRITE(6,'(A4)',ADVANCE="NO") ':H, '

       ELSE IF(1.9d0 < pans .AND. pans < 2.1d0)THEN
          pans = 2.d0 ! D
          WRITE(6,'(A4)',ADVANCE="NO") ':D, '

       ELSE IF(2.9d0 < pans .AND. pans < 3.1d0)THEN
          pans = 3.d0 ! T
          WRITE(6,'(A4)',ADVANCE="NO") ':T, '
       END IF

    ELSE IF(1.9d0 < pzns .AND. pzns < 2.1d0)THEN
       pans = 4.d0
       pzns = 2.d0 ! He
       WRITE(6,'(A5)',ADVANCE="NO") ':He, '

    ELSE IF(2.9d0 < pzns .AND. pzns < 3.1d0)THEN
       pans = 6.d0
       pzns = 3.d0 ! Li             
       WRITE(6,'(A5)',ADVANCE="NO") ':Li, '

    ELSE IF(3.9d0 < pzns .AND. pzns < 4.1d0)THEN
       pans = 8.d0
       pzns = 4.d0 ! Be
       WRITE(6,'(A5)',ADVANCE="NO") ':Be, '

    ELSE IF(4.9d0 < pzns .AND. pzns < 5.1d0)THEN
       pans = 10.d0
       pzns = 5.d0 ! B
       WRITE(6,'(A4)',ADVANCE="NO") ':B, '

    ELSE IF(5.9d0 < pzns .AND. pzns < 6.1d0)THEN
       pans = 12.d0
       pzns = 6.d0 ! C
       WRITE(6,'(A4)',ADVANCE="NO") ':C, '

    ELSE
       ! other impurity
       WRITE(6,'(A8,F9.3,A4,F9.3,A1)',ADVANCE="NO") &
                                    ': IMP(A=',pans, ', Z=',pzns,')'
    END IF

    RETURN
  END SUBROUTINE tr_uf_identify_ions

END MODULE truf0d
