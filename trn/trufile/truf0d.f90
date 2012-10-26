MODULE truf0d
  USE uflist, ONLY: uf0d

  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf0d,idnm,idnfast,idnmaz,idnfastaz,idzeff

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

       CALL tr_uf_get_species
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
       IF(uf0d(num)%lex == .TRUE. .AND. uf0d(num+1)%lex == .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnfastaz(nsi) = .TRUE.
          END IF
       END IF

       num = (nsi-1)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
       IF(uf0d(num)%lex == .TRUE. .AND. uf0d(num+1)%lex == .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnmaz(nsi) = .TRUE.
          END IF
       END IF

       ! The species can not be idetified in the following case.
       IF(idnfast(nsi)==.TRUE. .AND. idnfastaz(nsi)==.FALSE.)THEN
          ierr = ierr + 10**(nsi-1)
       END IF
       IF(idnm(nsi)==.TRUE. .AND. idnmaz(nsi)==.FALSE.)THEN
          ierr = ierr + 2*10**(nsi-1)
       END IF
    END DO
    
    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX tr_uf_check_species: NO consistency in species data.'
       WRITE(6,*) 'XX IERR= ',ierr
       RETURN
    END IF

    ! 0d.dat: kfid = ZEFF
    IF(uf0d(175)%lex == .TRUE. .AND. &
         0 < uf0d(175)%fr0 .AND. uf0d(175)%fr0 < pzmax)THEN
       idzeff(2) = .TRUE.
    END IF

    RETURN
  END SUBROUTINE tr_uf_check_species

! **********************************************************************

  SUBROUTINE tr_uf_get_species!(idnm,idnfast,idnmaz,idnfastaz,idzeff)
    USE trcomm, ONLY: pau,pzu
    IMPLICIT NONE

    CHARACTER(LEN=10) :: kfid
    CHARACTER(LEN=1)  :: knum
    INTEGER(ikind)    :: nsi, num

    pau(1:nsum) = 0.d0
    pzu(1:nsum) = 0.d0
    zeffu(1:ntxmax) = 0.d0

    pau(1) = ame/amp
    pzu(1) = -1.d0

    WRITE(6,*) ! spacing
    WRITE(6,'(A33)',ADVANCE="NO") ' ## UFILE data of ions contains: '
    DO nsi = 2, nsum
       WRITE(knum,'(I1)') ns
       kfid = 'NM'//knum

       num = (nsi-1)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
       IF(idnm(nsi))THEN
          IF(0.9d0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < 1.1d0)THEN
             pzu(nsi) = 1.d0 ! hydrogenic ions

             IF(0.9d0 < uf0d(num)%fr0 .AND. uf0d(num)%fr0 < 1.1d0)THEN
                pau(nsi) = 1.d0 ! H
                WRITE(6,'(I1,A4)',ADVANCE="NO") nsi,'.H, '

             ELSE IF(1.9d0 < uf0d(num)%fr0 .AND. uf0d(num)%fr0 < 2.1d0)THEN
                pau(nsi) = 2.d0 ! D
                WRITE(6,'(I1,A4)',ADVANCE="NO") nsi,'.D, '

             ELSE IF(2.9d0 < uf0d(num)%fr0 .AND. uf0d(num)%fr0 < 3.1d0)THEN
                pau(nsi) = 3.d0 ! T
                WRITE(6,'(I1,A4)',ADVANCE="NO") nsi,'.T, '
             END IF

          ELSE IF(1.9d0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < 2.1d0)THEN
             pau(nsi) = 4.d0
             pzu(nsi) = 2.d0 ! He
             WRITE(6,'(I1,A5)',ADVANCE="NO") nsi,'.He, '

          ELSE IF(2.9d0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < 3.1d0)THEN
             pau(nsi) = 6.d0
             pzu(nsi) = 3.d0 ! Li             
             WRITE(6,'(I1,A5)',ADVANCE="NO") nsi,'.Li, '

          ELSE IF(3.9d0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < 4.1d0)THEN
             pau(nsi) = 8.d0
             pzu(nsi) = 4.d0 ! Be
             WRITE(6,'(I1,A5)',ADVANCE="NO") nsi,'.Be, '

          ELSE IF(4.9d0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < 5.1d0)THEN
             pau(nsi) = 10.d0
             pzu(nsi) = 5.d0 ! B
             WRITE(6,'(I1,A4)',ADVANCE="NO") nsi,'.B, '

          ELSE IF(5.9d0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < 6.1d0)THEN
             pau(nsi) = 12.d0
             pzu(nsi) = 6.d0 ! C
             WRITE(6,'(I1,A4)',ADVANCE="NO") nsi,'.C, '

          ELSE
             ! other impurity
             pau(nsi) = uf0d(num)%fr0
             pzu(nsi) = uf0d(num+1)%fr0
             WRITE(6,'(I1,A7,F9.3,A4,F9.3,A1)',ADVANCE="NO") &
                  nsi,'. IMP(A=',pau(nsi), ', Z=',pzu(nsi),')'
          END IF
       END IF

    END DO  

    ELSE IF(idzeff(1))THEN
       zeffu(1:ntxmax) = uf0d(175)%fr0
    END IF

    WRITE(6,*) ! spacing

    RETURN
  END SUBROUTINE tr_uf_get_species

END MODULE truf0d
