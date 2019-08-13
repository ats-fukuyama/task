MODULE truf0d
  USE uflist, ONLY: uf0d
  USE trcomm, ONLY: ikind, rkind, nsum,ntum,nrum

  IMPLICIT NONE

  PUBLIC
  PRIVATE tr_uf_check_species

  LOGICAL,DIMENSION(1:nsum) :: idnm,idnfast,idnmaz,idnfastaz
  LOGICAL,DIMENSION(1:2)    :: idzeff
  LOGICAL                   :: idpgasa,idpgasz,idpimpa,idpimpz

  CHARACTER(LEN=15) :: toknam,shotnum,auxheat,phase,itb,itbtype
  REAL(rkind)       :: itbtime

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

       CALL tr_uf_get_mainpapz

       CALL tr_uf_check_species

       CALL tr_uf_get_nspapz

       CALL tr_uf_main_ions(ierr)

    ENDIF

    RETURN
  END SUBROUTINE tr_uf0d

! ***********************************************************************
! ***********************************************************************

  SUBROUTINE tr_uf_get_mainpapz
! ------------------------------------------------------------------------
!   - Set the atomic number and the charge number of main ion and impurity
!     ** default main ion: deuterium
!     ** default main impurity: carbon
! ------------------------------------------------------------------------
    USE trcomm, ONLY: pa_mion,pz_mion,pa_mimp,pz_mimp
    IMPLICIT NONE

    INTEGER(ikind) :: ierr

    idpgasa = uf0d(107)%lex
    idpgasz = uf0d(108)%lex
    idpimpa = uf0d(112)%lex
    idpimpz = uf0d(113)%lex

    pa_mion = uf0d(107)%fr0 ! PGASA: 'mean' mass number of the plasma ions
    pz_mion = uf0d(108)%fr0 ! PGASZ: 'mean' charge number of the plasma ions
    pa_mimp = uf0d(112)%fr0 ! PIMPA: mass number of the plasma main impurity
    pz_mimp = uf0d(113)%fr0 ! PIMPZ: charge number of the plasma main impurity

    ! Determination of main ion which is used for correction of profile.
    IF(idpgasa .EQ. .FALSE. .OR. idpgasz .EQ. .FALSE.)THEN
       pa_mion = 2.d0
       pz_mion = 1.d0 ! Deuterium: defalut

    ELSE
       CALL tr_uf_identify_ions(pa_mion,pz_mion,ierr)
       IF(ierr /= 0)THEN
          WRITE(6,'(1X,A,A5,A2,A5,A)') &
               '++ Warning: ',uf0d(107)%kfid,', ',uf0d(108)%kfid,' in 0d.dat file'
          IF(ierr == -1) STOP
       END IF

    END IF

    ! Determination of main impurity which is used for correction of profile.
    IF(idpimpa .EQ. .FALSE. .OR. idpimpz .EQ. .FALSE.)THEN
       pa_mimp = 12.d0
       pz_mimp = 6.d0 ! Carbon: default

    ELSE
       CALL tr_uf_identify_ions(pa_mimp,pz_mimp,ierr)
       IF(ierr /= 0)THEN
          WRITE(6,'(1X,A,A5,A2,A5,A)') &
               '++ Warning: ',uf0d(112)%kfid,', ',uf0d(113)%kfid,' in 0d.dat file'
          IF(ierr == -1) STOP
       END IF

    END IF

    RETURN
  END SUBROUTINE tr_uf_get_mainpapz

! ***********************************************************************

  SUBROUTINE tr_uf_check_species
! ------------------------------------------------------------------------ 
!   *** Checking whether impurity exists ***
!
!   This subrouitne must be called after 'tr_uf0d' is called.
! ------------------------------------------------------------------------
    USE trcomm, ONLY: pa_mion, pz_mion
    USE ufinit, ONLY: ufile_inquire
    IMPLICIT NONE

    INTEGER(ikind),PARAMETER :: pamax = 100, pzmax = 50

    CHARACTER(LEN=100) :: kfile, kfileb ! dummy for 'ufile_inquire'
    CHARACTER(LEN=10)  :: kfid, kufdim
    CHARACTER(LEN=1)   :: knum
    INTEGER(ikind)     :: num,nsi,id_bin,errout, ierr1,ierr2,nsi1,nsi2

    ierr1 = 0
    ierr2 = 0

    id_bin = 2 ! confirm the existence of only ASCII file
    errout = 1 ! suppress the error message in checking the existence of UFILE

    idnm(1:nsum)      = .FALSE.
    idnfast(1:nsum)   = .FALSE.
    idnmaz(1:nsum)    = .FALSE.
    idnfastaz(1:nsum) = .FALSE.
    idzeff(1:2)       = .FALSE.

    kufdim = '2d'

    kfid = 'NE'
    CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr1)
    IF(ierr1 == 0)THEN
       idnm(1)   = .TRUE.
       idnmaz(1) = .TRUE.
    END IF

    DO nsi = 2, nsum
       WRITE(knum,'(I1)') nsi-1

       ! NMn
       kfid = 'NM'//knum
       CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr1)
       IF(ierr1 == 0) idnm(nsi) = .TRUE.

       ! NFASTn
       kfid = 'NFAST'//knum
       CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr1)
       IF(ierr1 == 0) idnfast(nsi) = .TRUE.
    END DO

    kfid = 'ZEFFR'
    CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr1)
    IF(ierr1 == 0) idzeff(1) = .TRUE.


    kufdim = '1d'
    kfid   = 'ZEFF'
    CALL ufile_inquire(kufdim,kfid,kfile,kfileb,id_bin,errout,ierr1)
    IF(ierr1 == 0) idzeff(2) = .TRUE.


    ierr1 = 0
    ierr2 = 0

    DO nsi = 2, nsum
       num = (nsi-2)*2 + 68 ! see uf0d list: NFAST1A ~ NFAST9Z
       IF(uf0d(num)%lex .EQ. .TRUE. .AND. uf0d(num+1)%lex .EQ. .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnfastaz(nsi) = .TRUE.
          END IF
       END IF

       num = (nsi-2)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
       IF(uf0d(num)%lex .EQ. .TRUE. .AND. uf0d(num+1)%lex .EQ. .TRUE.)THEN
          IF(0 < uf0d(num)%fr0   .AND.   uf0d(num)%fr0 < pamax .AND.  &
             0 < uf0d(num+1)%fr0 .AND. uf0d(num+1)%fr0 < pzmax)THEN
             idnmaz(nsi) = .TRUE.
          END IF
       END IF

       ! The species can not be idetified in the following case.
       IF(idnfast(nsi).EQ. .TRUE. .AND. idnfastaz(nsi).EQ. .FALSE.)THEN
          ierr1 = ierr1 + 1
          nsi1  = nsi
       END IF
       IF(idnm(nsi).EQ. .TRUE. .AND. idnmaz(nsi).EQ. .FALSE.)THEN
          ierr2 = ierr2 + 1
          nsi2  = nsi
       END IF
    END DO

    ! complete atomic and charge numbers when the lacking data 
    !  is only one species in the presence of PGASA and PGASZ data.
    IF(ierr1 == 1 )THEN
       idnfastaz(nsi1) = .TRUE.
       num = (nsi1-2)*2 + 68

       uf0d(num)%lex   = .TRUE.
       uf0d(num+1)%lex = .TRUE.
       uf0d(num)%fr0   = pa_mion
       uf0d(num+1)%fr0 = pz_mion
       WRITE(6,*) ! spacing
       WRITE(6,*) '## tr_uf_check_species: beam ion information does not exist.'
       WRITE(6,'(1X,A21,I1,A2,ES9.2,A7,I1,A2,ES9.2)') &
               '## completed as NFAST',nsi1-1,'A=',pa_mion, &
                            '  NFAST',nsi1-1,'Z=',pz_mion
       ierr1 = 0
    END IF

    IF(ierr2 == 1 )THEN
       idnmaz(nsi2) = .TRUE.
       num = (nsi2-2)*2 + 86

       uf0d(num)%lex   = .TRUE.
       uf0d(num+1)%lex = .TRUE.
       uf0d(num)%fr0   = pa_mion
       uf0d(num+1)%fr0 = pz_mion
       WRITE(6,*) ! spacing
       WRITE(6,*) '## tr_uf_check_species: main ion information does not exist.'
       WRITE(6,'(1X,A18,I1,A2,ES9.2,A4,I1,A2,ES9.2)')   &
               '## completed as NM',nsi2-1,'A=',pa_mion, &
                            '  NM',nsi2-1,'Z=',pz_mion
       ierr2 = 0
    END IF

    IF(ierr1 /= 0 .OR. ierr2 /= 0)THEN
       WRITE(6,*) 'XX tr_uf_check_species: NO consistency in species data.'
       WRITE(6,'(1X,A11,I2,A10,I2)') 'XX IERR1 = ',ierr1,'  IERR2 = ',ierr2
       RETURN
    END IF


    RETURN
  END SUBROUTINE tr_uf_check_species

! **********************************************************************

  SUBROUTINE tr_uf_get_nspapz
! ----------------------------------------------------------------------
!   - Set the atomic number and the charge number of species which are
!      contained in experimental data
! ----------------------------------------------------------------------
    USE trcomm, ONLY: ame,amp,pau,pzu,pafu,pzfu
    IMPLICIT NONE
    INTEGER(ikind)    :: nsi, num, ierr
    REAL(rkind)       :: pans, pzns

    ierr = 0

    ! electron
    pau(1) = ame/amp
    pzu(1) = -1.d0
    
    ! bulk ions
    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          num = (nsi-2)*2 + 86 ! see uf0d list: NM1A ~ NM9Z
          pans = uf0d(num)%fr0
          pzns = uf0d(num+1)%fr0
          CALL tr_uf_identify_ions(pans,pzns,ierr)
          IF(ierr /= 0)THEN
             WRITE(6,'(1X,A,A4,A2,A4,A)') &
             '++ Warning: ',uf0d(num)%kfid,', ',uf0d(num+1)%kfid,' in 0d.dat file'
             IF(ierr == -1) STOP
          END IF
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
          CALL tr_uf_identify_ions(pans,pzns,ierr)
          IF(ierr /= 0)THEN
             WRITE(6,'(1X,A,A7,A2,A7,A)') &
             '++ Warning: ',uf0d(num)%kfid,', ',uf0d(num+1)%kfid,' in 0d.dat file'
             IF(ierr == -1) STOP
          END IF
          pafu(nsi) = pans
          pzfu(nsi) = pzns       
       END IF
    END DO
       
    RETURN
  END SUBROUTINE tr_uf_get_nspapz

! ************************************************************************

  SUBROUTINE tr_uf_main_ions(ierr)
! ------------------------------------------------------------------------
!   check of the existence of main ion and impurity profile data
! ------------------------------------------------------------------------
    USE trcomm, ONLY: mdlni,pau,pzu,pzfu,nsu_mion,nsu_mimp,nsu_fion, &
                      pa_mion,pz_mion,pa_mimp,pz_mimp
    IMPLICIT NONE

    INTEGER(ikind),INTENT(OUT)   :: ierr
    INTEGER(ikind) :: nsi,id_mion,id_mimp,id_fion

    ierr   = 0

    id_mion = .FALSE.
    id_mimp = .FALSE.
    id_fion = .FALSE.

    nsu_mion = 0
    nsu_mimp = 0
    nsu_fion = 0

    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          IF(pzu(nsi) == pz_mion .AND. pau(nsi) == pa_mion)THEN
             id_mion = .TRUE.
             nsu_mion = nsi
          ELSE IF(pzu(nsi) == pz_mimp .AND. pau(nsi) == pa_mimp)THEN
             id_mimp = .TRUE.
             nsu_mimp = nsi
          END IF
       END IF
       IF(idnfast(nsi) .EQ. .TRUE.)THEN
          id_fion = .TRUE.
          nsu_fion = nsi ! for now, just identifier of the existence
       END IF
    END DO

    ! the case that the main ion density data can not be found.
    IF(id_mion .EQ. .FALSE.)THEN
       DO nsi = 2, nsum
          IF(idnm(nsi) .EQ. .FALSE.)THEN
             idnmaz(nsi) = .TRUE.

             nsu_mion = nsi
             pzu(nsu_mion) = pz_mion
             pau(nsu_mion) = pa_mion
             EXIT
          END IF
          IF(nsi==nsum)THEN
             WRITE(6,*) 'XX tr_uf_main_ions: No blank in UFILE array.'
             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
             ierr = 1
             RETURN
          ENDIF
       END DO
    END IF

    ! the case that the main impurity density data can not be found.
    IF(id_mimp .EQ. .FALSE.)THEN
       DO nsi = 3, nsum ! nsi = 2: impurity data should not be stored.
          IF(idnm(nsi) .EQ. .FALSE.)THEN
             idnmaz(nsi) = .TRUE.

             nsu_mimp = nsi
             pzu(nsu_mimp) = pz_mimp
             pau(nsu_mimp) = pa_mimp
             EXIT
          END IF
          IF(nsi==nsum)THEN
             WRITE(6,*) 'XX tr_uf_main_ions: No blank in UFILE array.'
             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
             ierr = 1
             RETURN
          ENDIF
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_uf_main_ions

! ************************************************************************

  SUBROUTINE tr_uf_identify_ions(pans,pzns,ierr)
! -----------------------------------------------------------------------
!   Identify the atomic number and the charge number of ions
!
!   ierr =  0 : no error
!        =  1 : high charged impurity is detected, and is reset as defalut
!               impurity element (carbon) for now.
!        = -1 : fatal error; element cannot be identified
! -----------------------------------------------------------------------
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT)  :: pans, pzns
    INTEGER(ikind),INTENT(OUT) :: ierr

    REAL(rkind) :: eps

    eps = 0.2d0 ! allowance to identify the element
    ierr = 0

    IF(1.d0-eps < pzns .AND. pzns < 1.d0+eps)THEN
       pzns = 1.d0 ! hydrogenic ions

       IF(1.d0-eps < pans .AND. pans < 1.d0+eps)THEN
          pans = 1.d0 ! H

       ELSE IF(2.d0-eps < pans .AND. pans < 2.d0+eps)THEN
          pans = 2.d0 ! D

       ELSE IF(3.d0-eps < pans .AND. pans < 3.d0+eps)THEN
          pans = 3.d0 ! T
       END IF

    ELSE IF(2.d0-eps < pzns .AND. pzns < 2.d0+eps)THEN
       pzns = 2.d0

       IF(3.d0-eps < pans .AND. pans < 3.d0+eps)THEN
          pans = 3.d0 ! He3
       ELSE IF(4.d0-eps < pans .AND. pans < 4.d0+eps)THEN
          pans = 4.d0 ! He4
       END IF

    ELSE IF(3.d0-eps < pzns .AND. pzns < 3.d0+eps)THEN
       pzns = 3.d0 ! Li             
       pans = 6.d0

    ELSE IF(4.d0-eps < pzns .AND. pzns < 4.d0+eps)THEN
       pzns = 4.d0 ! Be
       pans = 8.d0

    ELSE IF(5.d0-eps < pzns .AND. pzns < 5.d0+eps)THEN
       pzns = 5.d0 ! B
       pans = 10.d0

    ELSE IF(6.d0-eps < pzns .AND. pzns < 6.d0+eps)THEN
       pzns = 6.d0 ! C
       pans = 12.d0

    ELSE IF( pzns > 6.d0)THEN
       ! other impurity; for now, set to defalut impurity: carbon
       ierr = 1   ! correction
       WRITE(6,*) ! line break
       WRITE(6,*) 'XX tr_uf_identify_ions: unsupported impurity element.'
       WRITE(6,'(1X,A,ES9.2,A8,ES9.2)') &
                  'XX unknown element: PZ = ',pzns, '   PA = ',pans
       pzns = 6.d0 ! C
       pans = 12.d0
       WRITE(6,*) '## tr_uf_identify_ions: set to defalut impurity values.'
       WRITE(6,'(1X,A,ES9.2,A8,ES9.2)') &
                  '## replaced element: PZ = ',pzns, '   PA = ',pans

    ELSE
       ierr = -1  ! fatal error
       WRITE(6,*) ! line break
       WRITE(6,*) 'XX tr_uf_identify_ions: cannot identify the element.'
       WRITE(6,'(1X,A,ES9.2,A8,ES9.2)') &
                  'XX unknown element: PZ = ',pzns, '   PA = ',pans
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
       IF(idnmaz(nsi))THEN
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
       IF(idnfastaz(nsi))THEN
          simf: DO nsa = 2, nsamax
             ns = ns_nsa(nsa)
             IF(idion(ns) /= 0.d0)THEN ! fast ions
                IF(pzfu(nsi) == pz(ns) .AND. pafu(nsi) == pa(ns))THEN
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

  SUBROUTINE tr_uf0dget_global(kfid,glbu)
!-------------------------------------------------------------------------
!   alternative routine for get global data in the absence of 1d.dat file
!-------------------------------------------------------------------------
    USE trcomm, ONLY: ntum,mdlxp
    IMPLICIT NONE

    CHARACTER(LEN=10),            INTENT(IN)  :: kfid
    REAL(rkind),DIMENSION(1:ntum),INTENT(OUT) :: glbu

    IF(mdlxp == 0)THEN ! UFILE
       IF(TRIM(kfid) == 'AMIN')THEN
          glbu(1:ntum)   = uf0d(1)%fr0
       ELSE IF(TRIM(kfid) == 'RGEO')THEN
          glbu(1:ntum)   = uf0d(128)%fr0
       ELSE IF(TRIM(kfid) == 'BT')THEN
          glbu(1:ntum)   = uf0d(15)%fr0
       ELSE IF(TRIM(kfid) == 'KAPPA')THEN
          glbu(1:ntum) = uf0d(56)%fr0
       ELSE IF(TRIM(kfid) == 'DELTA')THEN
          glbu(1:ntum) = uf0d(19)%fr0
       ELSE IF(TRIM(kfid) == 'IP')THEN
          glbu(1:ntum) = uf0d(51)%fr0
       END IF
    END IF

    RETURN
  END SUBROUTINE tr_uf0dget_global

! ************************************************************************

  SUBROUTINE tr_ufile_0d_view
    USE trcomm, ONLY: mdlxp,pau,pzu,pafu,pzfu, &
                      pa_mion,pz_mion,pa_mimp,pz_mimp
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
          IF(idnmaz(nsi))THEN
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
          IF(idnfastaz(nsi))THEN
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
       
       WRITE(6,fmt_arar) 'PGASA    = ',pa_mion,  'PGASZ    = ',pz_mion
       WRITE(6,fmt_arar) 'PIMPA    = ',pa_mimp,  'PIMPZ    = ',pz_mimp
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
