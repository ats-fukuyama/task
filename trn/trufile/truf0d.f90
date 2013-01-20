MODULE truf0d
  USE uflist, ONLY: uf0d
  USE trcomm, ONLY: ikind, rkind, nsum,ntum,nrum,ntxmax

  IMPLICIT NONE

  PUBLIC

  LOGICAL,DIMENSION(1:nsum) :: idnm,idnfast,idnmaz,idnfastaz
  LOGICAL,DIMENSION(1:2)    :: idzeff
  LOGICAL                   :: idpgasa,idpgasz,idpimpa,idpimpz

  CHARACTER(LEN=15) :: toknam,shotnum,auxheat,phase,itb,itbtype
  REAL(rkind) :: pgasa,pgasz,pimpa,pimpz, itbtime

CONTAINS

  SUBROUTINE tr_uf0d(mdlxp,ndmax,ierr)
! ----------------------------------------------------------------------
!   This subroutine acquires the 0D data only from UFILE format for now.
!    ( MDSplus is not supportted. )
! ----------------------------------------------------------------------
    USE ufread,ONLY: ufread_0d

    INTEGER(ikind),INTENT(IN)  :: mdlxp
    INTEGER(ikind),INTENT(OUT) :: ndmax,ierr

    INTEGER(ikind) :: ntx

    IF(mdlxp == 0)THEN
       CALL ufread_0d(ndmax,ierr)
       IF(ierr /= 0)THEN
          WRITE(6,'(A,I4)') '## tr_uf0d: Failed to read 0d.dat file. IERR= ',ierr
          ierr = 1
          RETURN
       END IF

       CALL tr_uf_check_species

       CALL tr_uf_get_papz

    ENDIF

    RETURN
  END SUBROUTINE tr_uf0d

! ***********************************************************************
! ***********************************************************************

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

    idnm(1:nsum)      = .FALSE.
    idnfast(1:nsum)   = .FALSE.
    idnmaz(1:nsum)    = .FALSE.
    idnfastaz(1:nsum) = .FALSE.
    idzeff(1:2)       = .FALSE.

    kufdim = '2d'

    kfid = 'NE'
    CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr)
    IF(ierr == 0)THEN
       idnm(1)   = .TRUE.
       idnmaz(1) = .TRUE.
    END IF

    DO nsi = 2, nsum
       WRITE(knum,'(I1)') nsi-1

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


    kufdim = '1d'
    kfid   = 'ZEFF'
    CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr)
    IF(ierr == 0) idzeff(2) = .TRUE.



    ierr = 0

    DO nsi = 2, nsum
       num = (nsi-2)*2 + 68 ! see uf0d list: NFAST1A ~ NFAST9Z
       IF(uf0d(num)%lex .EQV. .TRUE. .AND. uf0d(num+1)%lex .EQV. .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnfastaz(nsi) = .TRUE.
          END IF
       END IF

       num = (nsi-2)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
       IF(uf0d(num)%lex .EQV. .TRUE. .AND. uf0d(num+1)%lex .EQV. .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnmaz(nsi) = .TRUE.
          END IF
       END IF

       ! The species can not be idetified in the following case.
       IF(idnfast(nsi).EQV. .TRUE. .AND. idnfastaz(nsi).EQV. .FALSE.)THEN
          ierr = ierr + 10**(nsi-2)
       END IF
       IF(idnm(nsi).EQV. .TRUE. .AND. idnmaz(nsi).EQV. .FALSE.)THEN
          ierr = ierr + 2*10**(nsi-2)
       END IF
    END DO
    
    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX tr_uf_check_species: NO consistency in species data.'
       WRITE(6,*) 'XX IERR= ',ierr
       RETURN
    END IF


    RETURN
  END SUBROUTINE tr_uf_check_species

! **********************************************************************

  SUBROUTINE tr_uf_get_papz
! ----------------------------------------------------------------------
!   set the atomic number and the charge number of species which are
!    contained in experimental data
!
!   setup the conversion table of (nsu --> nsa)
! ----------------------------------------------------------------------
    USE trcomm, ONLY: ame,amp,pau,pzu,pafu,pzfu

    IMPLICIT NONE
    INTEGER(ikind)    :: nsi, num
    REAL(rkind)       :: pans, pzns

    ! electron
    pau(1) = ame/amp
    pzu(1) = -1.d0
    
    ! bulk ions
    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          num = (nsi-2)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
          pans = uf0d(num)%fr0
          pzns = uf0d(num+1)%fr0
          CALL tr_uf_identify_ions(pans,pzns)
          pau(nsi) = pans
          pzu(nsi) = pzns
       END IF
    END DO

    pafu(1) = 0.d0
    pzfu(1) = 0.d0

    ! fast ions
    DO nsi = 2, nsum
       IF(idnfast(nsi))THEN
          num  = (nsi-2)*2 + 68 ! see uf0d list: NFAST1A ~ NFAST9Z
          pans = uf0d(num)%fr0
          pzns = uf0d(num+1)%fr0
          CALL tr_uf_identify_ions(pans,pzns)
          pafu(nsi) = pans
          pzfu(nsi) = pzns       
       END IF
    END DO

    idpgasa = uf0d(107)%lex
    idpgasz = uf0d(108)%lex
    idpimpa = uf0d(112)%lex
    idpimpz = uf0d(113)%lex

    pgasa = uf0d(107)%fr0 ! 'mean' mass number of the plasma ions
    pgasz = uf0d(108)%fr0 ! 'mean' charge number of the plasma ions
    pimpa = uf0d(112)%fr0 ! mass number of the plasma main impurity
    pimpz = uf0d(113)%fr0 ! charge number of the plasma main impurity

       
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

       ELSE IF(1.9d0 < pans .AND. pans < 2.1d0)THEN
          pans = 2.d0 ! D

       ELSE IF(2.9d0 < pans .AND. pans < 3.1d0)THEN
          pans = 3.d0 ! T
       END IF

    ELSE IF(1.9d0 < pzns .AND. pzns < 2.1d0)THEN
       pans = 4.d0
       pzns = 2.d0 ! He

    ELSE IF(2.9d0 < pzns .AND. pzns < 3.1d0)THEN
       pans = 6.d0
       pzns = 3.d0 ! Li             

    ELSE IF(3.9d0 < pzns .AND. pzns < 4.1d0)THEN
       pans = 8.d0
       pzns = 4.d0 ! Be

    ELSE IF(4.9d0 < pzns .AND. pzns < 5.1d0)THEN
       pans = 10.d0
       pzns = 5.d0 ! B

    ELSE IF(5.9d0 < pzns .AND. pzns < 6.1d0)THEN
       pans = 12.d0
       pzns = 6.d0 ! C

    ELSE
       ! other impurity
       CONTINUE
    END IF

    RETURN
  END SUBROUTINE tr_uf_identify_ions

! ************************************************************************
! ************************************************************************

  SUBROUTINE tr_uf_set_table
! ------------------------------------------------------------------------
!   set the conversion table: [nsu --> nsa], [nsfu --> nsa]
! ------------------------------------------------------------------------
    USE trcomm, ONLY: idion,nsamax,nsa_nsu,nsa_nsfu,ns_nsa, &
                      pa,pz,pau,pzu,pafu,pzfu
    IMPLICIT NONE
    INTEGER(ikind) :: ns,nsa,nsi

    nsa_nsu(1:nsum)  = 0
    nsa_nsfu(1:nsum) = 0

    nsa_nsu(1)  = 1 ! electron
    nsa_nsfu(1) = 0 ! dummy

    exp: DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          sim: DO nsa = 2, nsamax
             ns = ns_nsa(nsa)
             IF(idion(ns) == 0.d0)THEN ! bulk ions
                IF(pzu(nsi) == pz(ns) .AND. pau(nsi) == pa(ns))THEN
                   nsa_nsu(nsi) = nsa
                   EXIT
                END IF
             END IF

             IF(nsa==nsamax)THEN
                WRITE(6,'(1X,A25,I1,A34)') &
       '## tr_uf_set_table:    NM',nsi-1,'  is discarded in the calculation.'
             END IF
          END DO sim
       END IF
    END DO exp

    expf: DO nsi = 2, nsum
       IF(idnfast(nsi))THEN
          simf: DO nsa = 2, nsamax
             ns = ns_nsa(nsa)
             IF(idion(ns) /= 0.d0)THEN ! fast ions
                IF(pzu(nsi) == pz(ns) .AND. pau(nsi) == pa(ns))THEN
                   nsa_nsfu(nsi) = nsa
                   EXIT
                END IF
             END IF

             IF(nsa==nsamax)THEN
                WRITE(6,'(1X,A25,I1,A34)') &
       '## tr_uf_set_table: NFAST',nsi-1,'  is discarded in the calculation.'
             END IF
          END DO simf
       END IF
    END DO expf

    RETURN
  END SUBROUTINE tr_uf_set_table

! ************************************************************************

  SUBROUTINE tr_ufile_0d_view
    USE trcomm, ONLY: mdlxp,pau,pzu,pafu,pzfu
    IMPLICIT NONE

    CHARACTER(LEN=10) :: kfid
    CHARACTER(LEN=1)  :: knum
    CHARACTER(LEN=30) :: fmt_aaaa,fmt_arar,fmt_aa
    INTEGER(ikind)    :: nsi

    fmt_aaaa = '(1X,A12,A15,A12,A15)'
    fmt_arar = '(1X,A11,ES10.2,3X,A11,ES10.2)'
    fmt_aa   = '(1X,A22,A15)'

    IF(mdlxp == 0)THEN ! UFILE
       shotnum = uf0d(139)%fc0
       toknam  = uf0d(158)%fc0
       auxheat = uf0d(3)%fc0
       phase   = uf0d(109)%fc0
       itb     = uf0d(53)%fc0
       itbtime = uf0d(54)%fr0
       itbtype = uf0d(55)%fc0

       WRITE(6,*) ! spacing
       WRITE(6,'(A55)') '# UFILE information ----------------------------------#'
       WRITE(6,fmt_aaaa) 'DEVICE    : ',toknam,'SHOT      : ',shotnum

       WRITE(6,'(1X,A12)',ADVANCE='NO') 'Bulk ions : '
       DO nsi = 2, nsum
          IF(idnm(nsi))THEN
             WRITE(knum,'(I1)') nsi-1
             kfid = 'NM'//knum
             WRITE(6,'(A6)',ADVANCE='NO') TRIM(kfid)
          
             IF(pzu(nsi) == 1.d0)THEN
                IF(pau(nsi)==1.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : H,  '
                IF(pau(nsi)==2.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : D,  '
                IF(pau(nsi)==3.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : T,  '
             ELSE IF(pzu(nsi) == 2.d0)THEN
                IF(pau(nsi)==4.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : He, '
             ELSE IF(pzu(nsi) == 3.d0)THEN
                IF(pau(nsi)==6.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : Li, '
             ELSE IF(pzu(nsi) == 4.d0)THEN
                IF(pau(nsi)==8.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : Be, '
             ELSE IF(pzu(nsi) == 5.d0)THEN
                IF(pau(nsi)==10.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : B,  '
             ELSE IF(pzu(nsi) == 6.d0)THEN
                IF(pau(nsi)==12.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : C,  '
             ELSE
                WRITE(6,'(A9,ES9.2,A4,ES9.2,A1)',ADVANCE="NO") &
                     ' : IMP(A=',pau(nsi), ', Z=',pzu(nsi),')'
             END IF
          END IF
       END DO
       WRITE(6,*) ! line break
       WRITE(6,'(1X,A12)',ADVANCE='NO') 'Fast ions : '
       DO nsi = 2, nsum
          IF(idnfast(nsi))THEN
             WRITE(knum,'(I1)') nsi-1
             kfid = 'NFAST'//knum
             WRITE(6,'(A6)',ADVANCE='NO') TRIM(kfid)
             
             IF(pzfu(nsi) == 1.d0)THEN
                IF(pafu(nsi)==1.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : H,  '
                IF(pafu(nsi)==2.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : D,  '
                IF(pafu(nsi)==3.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : T,  '
             ELSE IF(pzfu(nsi) == 2.d0)THEN
                IF(pafu(nsi)==4.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : He, '
             ELSE IF(pzfu(nsi) == 3.d0)THEN
                IF(pafu(nsi)==6.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : Li, '
             ELSE IF(pzfu(nsi) == 4.d0)THEN
                IF(pafu(nsi)==8.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : Be, '
             ELSE IF(pzfu(nsi) == 5.d0)THEN
                IF(pafu(nsi)==10.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : B,  '
             ELSE IF(pzfu(nsi) == 6.d0)THEN
                IF(pafu(nsi)==12.d0) WRITE(6,'(A7)',ADVANCE='NO') ' : C,  '
             ELSE
                WRITE(6,'(A9,ES9.2,A4,ES9.2,A1)',ADVANCE="NO") &
                     ' : IMP(A=',pafu(nsi), ', Z=',pzfu(nsi),')'
             END IF
          END IF
       END DO
       WRITE(6,*) ! line break
       
       WRITE(6,fmt_arar) 'PGASA    = ',pgasa,  'PGASZ    = ',pgasz
       WRITE(6,fmt_arar) 'PIMPA    = ',pimpa,  'PIMPZ    = ',pimpz
       WRITE(6,fmt_aa) 'Auxiliary heating   : ',auxheat
       WRITE(6,fmt_aa) 'Plasma phase        : ',phase
       WRITE(6,'(1X,A22,A15,A1,A15)') &
                       'ITB condition, type : ',itb,',',itbtype
       WRITE(6,'(1X,A22,ES10.2)')     &
                       'ITB triggering time : ',itbtime
       WRITE(6,'(A55)') '#-----------------------------------------------------#'
       WRITE(6,*) ! spacing
       
    END IF

    RETURN
  END SUBROUTINE tr_ufile_0d_view

END MODULE truf0d
