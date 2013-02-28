MODULE trufcalc
! ------------------------------------------------------------------------
!   Correction of experimental data, and calculation of associated values
! ------------------------------------------------------------------------

  USE trcomm, ONLY: ikind,rkind,nrum,ntum,nsum,ntxmax,nrmax
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf_nicomplete, &
         sumzni,rnuc,rnfuc,zeffruc ! for graphic output

  ! the profiles before correction by neutrality for graphic output
  REAL(rkind),DIMENSION(1:ntum,1:nrum)        :: zeffruc, sumzni
  REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum) :: rnuc,rnfuc

CONTAINS

  SUBROUTINE tr_uf_nicomplete
! -----------------------------------------------------------------------
!   Correct the profile of main ion, main impurity and effective charge
!    number(Zeff) based on the assumption of charge neutrality.
! -----------------------------------------------------------------------
    USE trcomm, ONLY: pzu,pzfu,rnu,rnfu,zeffru,mdlni, &
         nsu_mion,nsu_mimp,pa_mion,pz_mion,pa_mimp,pz_mimp
    IMPLICIT NONE

    LOGICAL        :: id_mion, id_mimp, test
    INTEGER(ikind) :: id,nsu,nsi,ierr

    ierr = 0

    rnuc(1:nsum,1:ntum,1:nrum)  = 0.d0
    rnfuc(1:nsum,1:ntum,1:nrum) = 0.d0
    zeffruc(1:ntum,1:nrum)      = 0.d0

    ! save original profiles
    DO nsu = 1, nsum
       rnuc(nsu,1:ntxmax,1:nrmax+1)  = rnu(nsu,1:ntxmax,1:nrmax+1)
       rnfuc(nsu,1:ntxmax,1:nrmax+1) = rnfu(nsu,1:ntxmax,1:nrmax+1)
    END DO
    zeffruc(1:ntxmax,1:nrmax+1) = zeffru(1:ntxmax,1:nrmax+1)

    ! check of the existence of main ion and impurity, and correct MDLNI
    CALL tr_uf_nicheck(ierr)

    IF(ierr == 0)THEN
       ! correction of profiles based on charge neutrality
       CALL tr_uf_neutrality(rnu,pzu,rnfu,pzfu,nsu_mion,nsu_mimp,&
                             zeffru,mdlni)
    END IF

    RETURN
  END SUBROUTINE tr_uf_nicomplete


! **********************************************************************

  SUBROUTINE tr_uf_nicheck(ierr)
! ----------------------------------------------------------------------
!   check the consistency of experimental data and MDLNI
! ----------------------------------------------------------------------
    USE truf0d, ONLY: idnm,idnfast,idnmaz,idnfastaz,idzeff
    USE trcomm, ONLY: mdlni,pau,pzu,pzfu,pa_mion,pz_mion,pa_mimp,pz_mimp, &
                      nsu_mion,nsu_mimp,nsu_fion
    IMPLICIT NONE

    INTEGER(ikind),INTENT(OUT)   :: ierr
    INTEGER(ikind) :: nsi

    ierr = 0

    SELECT CASE(mdlni)
    CASE(1) ! complete n_i and n_imp from Zeff, n_e (and n_bulk,n_fast)
       DO
          IF((idzeff(1).EQV..TRUE.) .OR. (idzeff(2).EQV..TRUE.))THEN
             IF(idnm(nsu_mion).EQV..FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp).EQV..FALSE.) idnm(nsu_mimp) = .TRUE.
             EXIT
          END IF

          WRITE(6,*) 'XX tr_uf_nicheck: Lack of ZEFF data.'
          IF(idnm(nsu_mion))THEN
             mdlni = 2
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(idnm(nsu_mimp))THEN
             mdlni = 3
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 1
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
             IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
          END IF
          EXIT
       END DO

    CASE(2) ! complete n_imp and Zeff from n_e, n_i (and n_bulk,n_fast)
       IF(idnm(nsu_mimp).EQV..FALSE.) idnm(nsu_mimp) = .TRUE.

       IF(idnm(nsu_mion).EQV..FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_i data.'
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             mdlni = 1
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(idnm(nsu_mimp))THEN
             mdlni = 3
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 1
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
             IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
          END IF
       END IF

    CASE(3) ! complete n_i and Zeff   from n_e, n_imp (and n_bulk,n_fast)
       IF(idnm(nsu_mion).EQV..FALSE.) idnm(nsu_mion) = .TRUE.

       IF(idnm(nsu_mimp) .EQV. .FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_imp data.'
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             mdlni = 1
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(idnm(nsu_mion))THEN
             mdlni = 2
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 1
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
             IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
          END IF
       END IF

    CASE(4) ! complete n_i and n_imp from Zeff, n_e (and n_bulk)
       DO
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             IF(idnm(nsu_mion).EQV..FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp).EQV..FALSE.) idnm(nsu_mimp) = .TRUE.
             EXIT
          END IF

          WRITE(6,*) 'XX tr_uf_nicheck: Lack of ZEFF data.'
          IF(idnm(nsu_mion))THEN
             mdlni = 5
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(idnm(nsu_mimp))THEN
             mdlni = 6
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 1
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
             IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
          END IF
          EXIT
       END DO

    CASE(5) ! complete n_imp and Zeff from n_e, n_i (and n_bulk)
       IF(idnm(nsu_mimp).EQV..FALSE.) idnm(nsu_mimp) = .TRUE.

       IF(idnm(nsu_mion).EQV..FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_i data.'
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             mdlni = 4
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(idnm(nsu_mimp))THEN
             mdlni = 6
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 1
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
             IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
          END IF
       END IF

    CASE(6) ! complete n_i and Zeff   from n_e, n_imp (and n_bulk)
       IF(idnm(nsu_mion).EQV..FALSE.) idnm(nsu_mion) = .TRUE.

       IF(idnm(nsu_mimp) .EQV. .FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_imp data.'
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             mdlni = 4
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(idnm(nsu_mion))THEN
             mdlni = 5
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 1
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
             IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
             IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
          END IF
       END IF

    CASE(8) ! complete n_i and Zeff using the assumption n_e = n_i + n_f
       IF(idnm(nsu_mion).EQV..FALSE.) idnm(nsu_mion) = .TRUE.

       IF(idnfast(nsu_fion) .EQV. .FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_fast data.'
!          WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!          ierr = 1
          mdlni = 9
          WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
          IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.
       END IF
          
    CASE(9) ! complete n_i and Zeff using the assumption n_e = n_i
       IF(idnm(nsu_mion) .EQV. .FALSE.) idnm(nsu_mion) = .TRUE.
       IF(idnm(nsu_mimp) .EQV. .TRUE.)  idnm(nsu_mimp) = .FALSE.

    CASE DEFAULT
       WRITE(6,*) 'XX tr_uf_nicheck: the value of "MDLNI" is invalid. MDLNI= ',mdlni
       WRITE(6,*) 'XX Stopped to correction of profiles by charge neutrality.'
       ierr = -1
    END SELECT


    RETURN
  END SUBROUTINE tr_uf_nicheck

! **********************************************************************

  SUBROUTINE tr_uf_neutrality(rnc,pzc,rnfc,pzfc,ns_ion,ns_imp,zeffrc,id)
! ----------------------------------------------------------------------
!  < NOTATION >
!  n_e    : particle density of electron
!  n_i    : particle density of main ion
!  n_bulk : particle density of other ions 
!                               (not including main ion and impurity)
!  n_imp  : particle density of impurity
!
!  ns_ion : the subscript of arrays for main ion
!  ns_imp : the subscript of arrays for main impurity
!
!  rnc(nsu=1) must have electron density profile.
!  * rnfc(nsu=1) is dummy.
!
!  id --> mdlni
!  id = 1 : complete n_i and n_imp  from Zeff, n_e  (and n_bulk, n_fast)
!  id = 2 : complete n_imp and Zeff from n_e, n_i   (and n_bulk, n_fast)
!  id = 3 : complete n_i and Zeff   from n_e, n_imp (and n_bulk, n_fast)
!  id = 4 : complete n_i and n_imp  from Zeff, n_e  (and n_bulk)
!  id = 5 : complete n_imp and Zeff from n_e, n_i   (and n_bulk)
!  id = 6 : complete n_i and Zeff   from n_e, n_imp (and n_bulk)
!  id = 8 : complete n_i and Zeff using the assumption n_e = n_i (+ n_fast)
!  id = 9 : complete n_i and Zeff using the assumption n_e = n_i
! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(ikind),               INTENT(IN)  :: ns_ion,ns_imp,id
    REAL(rkind),DIMENSION(1:nsum),INTENT(IN)  :: pzc,pzfc
    REAL(rkind),DIMENSION(1:ntum,1:nrum),       INTENT(INOUT) :: zeffrc
    REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum),INTENT(INOUT) :: rnfc
    REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum),INTENT(INOUT) :: rnc

    INTEGER(ikind) :: nsi,nsu,ntx,nr

    sumzni(1:ntum,1:nrum)  = 0.d0

    SELECT CASE(id)
    CASE(1) ! Zeff, n_e (and n_bulk, n_fast) --> n_i, n_imp
       DO ntx = 1, ntxmax
          rnc(ns_ion,ntx,1:nrmax+1) = (pzc(ns_imp)-zeffrc(ntx,1:nrmax+1)) &
                                      *rnc(1,ntx,1:nrmax+1)
          rnc(ns_imp,ntx,1:nrmax+1) = (zeffrc(ntx,1:nrmax+1)-pzc(ns_ion)) &
                                      *rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_ion .AND. nsi /= ns_imp)THEN
                rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)     &
                     - (pzc(ns_imp)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)

                rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)     &
                     - (pzc(nsi)-pzc(ns_ion))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
             rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)        &
                  - (pzc(ns_imp)-pzfc(nsi))*pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)

             rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)        &
                  - (pzfc(nsi)-pzc(ns_ion))*pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)
          END DO

          rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)           &
                                     /(pzc(ns_ion)*(pzc(ns_imp)-pzc(ns_ion)))

          rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)           &
                                     /(pzc(ns_imp)*(pzc(ns_imp)-pzc(ns_ion)))
       END DO

    ! CASE(2) and CASE(3) are symmetry w.r.t 'ns_ion' or 'ns_imp'
    CASE(2) ! n_e, n_i (and n_bulk, n_fast) --> n_imp, Zeff
       DO ntx = 1, ntxmax
          rnc(ns_imp,ntx,1:nrmax+1) = rnc(1,ntx,1:nrmax+1)

          zeffrc(ntx,1:nrmax+1) = pzc(ns_imp)*rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_imp)THEN
                rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1) &
                                 - pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)

                zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)        &
                     - (pzc(ns_imp)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
             rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)   &
                             - pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)

             zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)           &
                  - (pzc(ns_imp)-pzfc(nsi))*pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)
          END DO
          rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)/pzc(ns_imp)

          zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)/rnc(1,ntx,1:nrmax+1)
       END DO

    CASE(3) ! n_e, n_imp (and n_bulk, n_fast) --> n_i, Zeff
       DO ntx = 1, ntxmax
          rnc(ns_ion,ntx,1:nrmax+1) = rnc(1,ntx,1:nrmax+1)

          zeffrc(ntx,1:nrmax+1) = pzc(ns_ion)*rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_ion)THEN
                rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1) &
                                  -pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)

                zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)         &
                     - (pzc(ns_ion)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
             rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)    &
                              -pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)

             zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)            &
                  - (pzc(ns_ion)-pzfc(nsi))*pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)
          END DO
          rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)/pzc(ns_ion)

          zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)/rnc(1,ntx,1:nrmax+1)
       END DO

    CASE(4) ! Zeff, n_e (and n_bulk) --> n_i, n_imp
       DO ntx = 1, ntxmax
          rnc(ns_ion,ntx,1:nrmax+1) = (pzc(ns_imp)-zeffrc(ntx,1:nrmax+1)) &
                                      *rnc(1,ntx,1:nrmax+1)
          rnc(ns_imp,ntx,1:nrmax+1) = (zeffrc(ntx,1:nrmax+1)-pzc(ns_ion)) &
                                      *rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_ion .AND. nsi /= ns_imp)THEN
                rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)     &
                     - (pzc(ns_imp)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)

                rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)     &
                     - (pzc(nsi)-pzc(ns_ion))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
          END DO

          rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)           &
                                     /(pzc(ns_ion)*(pzc(ns_imp)-pzc(ns_ion)))

          rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)           &
                                     /(pzc(ns_imp)*(pzc(ns_imp)-pzc(ns_ion)))
       END DO

    ! CASE(5) and CASE(6) are symmetry w.r.t 'ns_ion' or 'ns_imp'
    CASE(5) ! n_e, n_i (and n_bulk) --> n_imp, Zeff
       DO ntx = 1, ntxmax
          rnc(ns_imp,ntx,1:nrmax+1) = rnc(1,ntx,1:nrmax+1)

          zeffrc(ntx,1:nrmax+1) = pzc(ns_imp)*rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_imp)THEN
                rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1) &
                                 - pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)

                zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)        &
                     - (pzc(ns_imp)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
          END DO
          rnc(ns_imp,ntx,1:nrmax+1) = rnc(ns_imp,ntx,1:nrmax+1)/pzc(ns_imp)

          zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)/rnc(1,ntx,1:nrmax+1)
       END DO

    CASE(6) ! n_e, n_imp (and n_bulk) --> n_i, Zeff
       DO ntx = 1, ntxmax
          rnc(ns_ion,ntx,1:nrmax+1) = rnc(1,ntx,1:nrmax+1)

          zeffrc(ntx,1:nrmax+1) = pzc(ns_ion)*rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_ion)THEN
                rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1) &
                                  -pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)

                zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)         &
                     - (pzc(ns_ion)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
          END DO
          rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)/pzc(ns_ion)

          zeffrc(ntx,1:nrmax+1) = zeffrc(ntx,1:nrmax+1)/rnc(1,ntx,1:nrmax+1)
       END DO


    CASE(8) ! n_e, n_fast --> n_i, Zeff (n_e = n_i + n_fast; not including impurity)
       rnc(ns_ion,1:ntxmax,1:nrmax+1) = rnc(1,1:ntxmax,1:nrmax+1)

       zeffrc(1:ntxmax,1:nrmax+1) = pzc(ns_ion)*rnc(1,1:ntxmax,1:nrmax+1)
       DO nsi = 2, nsum
          rnc(ns_ion,1:ntxmax,1:nrmax+1) = rnc(ns_ion,1:ntxmax,1:nrmax+1) &
                                  -pzfc(nsi)*rnfc(nsi,1:ntxmax,1:nrmax+1)

          zeffrc(1:ntxmax,1:nrmax+1) = zeffrc(1:ntxmax,1:nrmax+1)         &
             - (pzc(ns_ion)-pzfc(nsi))*pzfc(nsi)*rnfc(nsi,1:ntxmax,1:nrmax+1)
       END DO

       zeffrc(1:ntxmax,1:nrmax+1) = zeffrc(1:ntxmax,1:nrmax+1)            &
                                    /rnc(1,1:ntxmax,1:nrmax+1)

    CASE(9) ! n_e --> n_i, Zeff (n_e = n_i)
       rnc(ns_ion,1:ntxmax,1:nrmax+1) = rnc(1,1:ntxmax,1:nrmax+1)
       zeffrc(1:ntxmax,1:nrmax+1)     = pzc(ns_ion)

    END SELECT

    ! correction for negative density
    FORALL(nsu=1:nsum,ntx=1:ntum,nr=1:nrmax+1,rnc(nsu,ntx,nr) < 0.d0)
       rnc(nsu,ntx,nr) = 1.d-6
    END FORALL
    FORALL(nsu=1:nsum,ntx=1:ntum,nr=1:nrmax+1,rnfc(nsu,ntx,nr) < 0.d0)
       rnfc(nsu,ntx,nr) = 1.d-9
    END FORALL

    ! check for magnitude of correction
    DO nsi = 2, nsum ! ion only: SUM_i(Z_i * n_i)
       sumzni(1:ntxmax,1:nrmax+1) = sumzni(1:ntxmax,1:nrmax+1)     &
                         + pzc(nsi)*rnc(nsi,1:ntxmax,1:nrmax+1)    &
                         + pzfc(nsi)*rnfc(nsi,1:ntxmax,1:nrmax+1)
    END DO

    RETURN
  END SUBROUTINE tr_uf_neutrality

END MODULE trufcalc
