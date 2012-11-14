MODULE trufcalc
! ------------------------------------------------------------------------
!   Correction of experimental data, and calculation of associated values
! ------------------------------------------------------------------------

  USE trcomm, ONLY: ikind,rkind,nrum,ntum,nsum,ntxmax,nrmax
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf_nicomplete, &
       sumzni,rnuc,rnfuc,zeffruc,ns_mion,ns_mimp

  REAL(rkind),DIMENSION(1:ntum,1:nrum)        :: zeffruc, sumzni
  REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum) :: rnuc,rnfuc

  INTEGER(ikind) :: ns_mion, ns_mimp

CONTAINS

  SUBROUTINE tr_uf_nicomplete
! -----------------------------------------------------------------------
!   Correct the profile of main ion, main impurity and effective charge
!    number(Zeff) based on the assumption of charge neutrality.
! -----------------------------------------------------------------------
    USE trcomm, ONLY: pzu,pzfu,rnu,rnfu,zeffru,mdlni

    IMPLICIT NONE

    LOGICAL        :: id_mion, id_mimp, test
    INTEGER(ikind) :: id,nsu,nsi,ierr
    REAL(rkind)    :: pa_mion,pz_mion,pa_mimp,pz_mimp

    rnuc(1:nsum,1:ntum,1:nrum)  = 0.d0
    rnfuc(1:nsum,1:ntum,1:nrum) = 0.d0
    zeffruc(1:ntum,1:nrum)      = 0.d0

    ! save original profiles
    DO nsu = 1, nsum
       rnuc(nsu,1:ntxmax,1:nrmax+1)  = rnu(nsu,1:ntxmax,1:nrmax+1)
       rnfuc(nsu,1:ntxmax,1:nrmax+1) = rnfu(nsu,1:ntxmax,1:nrmax+1)
    END DO
    zeffruc(1:ntxmax,1:nrmax+1) = zeffru(1:ntxmax,1:nrmax+1)

    ! identify the main ion and impurity (their 'pa' and 'pz')
    CALL tr_uf_id_main(pa_mion,pz_mion,pa_mimp,pz_mimp)

    ! check of the existence of main ion and impurity, and correct MDLNI
    CALL tr_uf_nicheck(pa_mion,pz_mion,pa_mimp,pz_mimp,ns_mion,ns_mimp,ierr)

    IF(ierr == 0)THEN
       ! correction of profiles based on charge neutrality
       CALL tr_uf_neutrality(rnu,pzu,rnfu,pzfu,ns_mion,ns_mimp,&
                             zeffru,mdlni)
    END IF

    RETURN
  END SUBROUTINE tr_uf_nicomplete


! **********************************************************************

  SUBROUTINE tr_uf_id_main(pa_mion,pz_mion,pa_mimp,pz_mimp)
! ----------------------------------------------------------------------
!   Identify the main ion and main impurity which are re-caluculated
!    to satisfy the charge neutrality.
!   ** default main ion: deuterium
!   ** default main impurity: carbon
! ----------------------------------------------------------------------
    USE truf0d, ONLY: idnmaz,idnfastaz,idpgasa,idpgasz,idpimpa,idpimpz, &
                      pgasa,pgasz,pimpa,pimpz
    IMPLICIT NONE

!    INTEGER(ikind),INTENT(OUT) :: ierr
    REAL(rkind),   INTENT(OUT) :: pa_mion,pz_mion,pa_mimp,pz_mimp


    pa_mion = pgasa ! 'mean' mass number of the plasma ions
    pz_mion = pgasz ! 'mean' charge number of the plasma ions

    pa_mimp = pimpa ! mass number of the plasma main impurity
    pz_mimp = pimpz ! charge number of the plasma main impurity

    ! Determination of main ion which is used for correction of profile.
    IF(idpgasa .EQV. .FALSE. .OR. idpgasz .EQV. .FALSE.)THEN
       pa_mion = 2.d0
       pz_mion = 1.d0 ! Deuterium: defalut
    ELSE
       IF(pa_mion < 1.5d0)THEN
          pa_mion = 1.d0 ! Light Hydrogen
       ELSE IF(1.5d0 <= pa_mion .AND. pa_mion < 2.5d0)THEN
          pa_mion = 2.d0 ! Deuterium
       ELSE IF(2.5d0 <= pa_mion .AND. pa_mion < 3.5d0)THEN
          pa_mion = 3.d0 ! Tritium
       ELSE IF(pa_mion < 4.d0)THEN
          pa_mion = 4.d0 ! Herium
       END IF
       
       IF(pz_mion < 1.5d0)THEN
          pz_mion = 1.d0 ! Hydrogenic ion
       ELSE
          pz_mion = 2.d0 ! Herium
       END IF

       IF(pa_mion <= 3.d0 .AND. pz_mion == 2.d0)THEN
          WRITE(6,*) 'XX tr_uf_comp_check: Main ion identification error.'
          WRITE(6,*) 'PGASA= ',pgasa, 'PGASZ= ',pgasz
          pa_mion = 2.d0
          pz_mion = 1.d0 ! Deuterium: defalut
          WRITE(6,*) '## tr_uf_comp_check: set to default values.'
          WRITE(6,*) 'PGASA= ',pa_mion, 'PGASZ= ',pz_mion
!          ierr = 1
       END IF
    END IF

    ! Determination of main impurity which is used for correction of profile.
    IF(idpimpa .EQV. .FALSE. .OR. idpimpz .EQV. .FALSE.)THEN
       pa_mimp = 12.d0
       pz_mimp = 6.d0 ! Carbon: default
    ELSE
       DO
          IF(pa_mimp < 9.d0)THEN
             pa_mimp = 8.d0 ! Beryllium
          ELSE IF(9.d0 <= pa_mimp .AND. pa_mimp < 11.d0)THEN
             pa_mimp = 10.d0 ! Boron
          ELSE IF(11.d0 <= pa_mimp .AND. pa_mimp < 13.d0)THEN
             pa_mimp = 12.d0 ! Carbon
          ELSE
             WRITE(6,*) 'XX tr_uf_comp_check: Unsupported impurity element.'
             WRITE(6,*) 'PIMPA= ',pimpa, 'PIMPZ= ',pimpz
             pa_mimp = 12.d0
             pz_mimp = 6.d0 ! Carbon: default
             WRITE(6,*) '## tr_uf_comp_check: set to default values.'
             WRITE(6,*) 'PIMPA= ',pa_mimp, 'PIMPZ= ',pz_mimp
!             ierr = 2
             EXIT
          END IF
    
          IF(pz_mimp < 4.5d0)THEN
             pz_mimp = 4.d0 ! Beryllium
          ELSE IF(4.5d0 <= pz_mimp .AND. pz_mimp < 5.5d0)THEN
             pz_mimp = 5.d0 ! Boron
          ELSE IF(5.5d0 <= pz_mimp .AND. pz_mimp < 6.5d0)THEN
             pz_mimp = 6.d0 ! Carbon
          ELSE
             WRITE(6,*) 'XX tr_uf_comp_check: Unsupported impurity element.'
             WRITE(6,*) 'PIMPA= ',pimpa, 'PIMPZ= ',pimpz
             pa_mimp = 12.d0
             pz_mimp = 6.d0 ! Carbon: default
             WRITE(6,*) '## tr_uf_comp_check: set to default values.'
             WRITE(6,*) 'PIMPA= ',pa_mimp, 'PIMPZ= ',pz_mimp
!             ierr = 2
             EXIT
          END IF
          EXIT
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_uf_id_main

! **********************************************************************

  SUBROUTINE tr_uf_nicheck(pa_mion,pz_mion,pa_mimp,pz_mimp,ns_ion,ns_imp,ierr)
! ----------------------------------------------------------------------
!   check of the existence of main ion and impurity, and correct MDLNI
! ----------------------------------------------------------------------
    USE truf0d, ONLY: idnm,idnfast,idnmaz,idnfastaz,idzeff
    USE trcomm, ONLY: pau,pzu,pzfu, mdlni
    IMPLICIT NONE

    REAL(rkind),   INTENT(IN)    :: pa_mion,pz_mion,pa_mimp,pz_mimp
    INTEGER(ikind),INTENT(OUT)   :: ns_ion, ns_imp, ierr

    LOGICAL        :: id_ion,id_imp
    INTEGER(ikind) :: nsi

    id_ion = .FALSE.
    id_imp = .FALSE.
    ierr   = 0

    ns_ion = 0
    ns_imp = 0
    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          IF(pzu(nsi) == pz_mion .AND. pau(nsi) == pa_mion)THEN
             id_ion = .TRUE.
             ns_ion = nsi
          ELSE IF(pzu(nsi) == pz_mimp .AND. pau(nsi) == pa_mimp)THEN
             id_imp = .TRUE.
             ns_imp = nsi
          END IF
       END IF
    END DO

    ! the case that the main ion density data can not be found.
    IF(ns_ion == 0)THEN
       DO nsi = 2, nsum
          IF(idnm(nsi) .EQV. .FALSE.)THEN
             idnm(nsi)   = .TRUE.
             idnmaz(nsi) = .TRUE.

             ns_ion = nsi
             pzu(ns_ion) = pz_mion
             pau(ns_ion) = pa_mion
             EXIT
          END IF
          IF(nsi==nsum)THEN
             WRITE(6,*) 'XX tr_uf_nicheck: No blank in UFILE array.'
             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
             ierr = 1
             RETURN
          ENDIF
       END DO
    END IF

    ! the case that the main impurity density data can not be found.
    IF(ns_imp == 0)THEN
       DO nsi = 3, nsum ! nsi = 2: impurity data should not be store.
          IF(idnm(nsi) .EQV. .FALSE.)THEN
             idnm(nsi)   = .TRUE.
             idnmaz(nsi) = .TRUE.

             ns_imp = nsi
             pzu(ns_imp) = pz_mimp
             pau(ns_imp) = pa_mimp
             EXIT
          END IF
          IF(nsi==nsum)THEN
             WRITE(6,*) 'XX tr_uf_nicheck: No blank in UFILE array.'
             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
             ierr = 1
             RETURN
          ENDIF
       END DO
    END IF


    SELECT CASE(mdlni)
    CASE(1) ! complete n_i and n_imp from Zeff, n_e (and n_bulk)
       IF(idzeff(1).EQV..FALSE. .AND. idzeff(2).EQV..FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of ZEFF data.'
          IF(id_ion)THEN
             mdlni = 2
             WRITE(6,*) '## "MDNI" is reset to MDLNI= ',mdlni
          ELSE IF(id_imp)THEN
             mdlni = 3
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 2
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          END IF
       END IF

    CASE(2) ! complete n_imp and Zeff from n_e, n_i (and n_bulk)
       IF(id_ion .EQV. .FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_i data.'
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             mdlni = 1
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(id_imp)THEN
             mdlni = 3
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 2
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          END IF
       END IF

    CASE(3) ! complete n_i and Zeff   from n_e, n_imp (and n_bulk)
       IF(id_imp .EQV. .FALSE.)THEN
          WRITE(6,*) 'XX tr_uf_nicheck: Lack of n_imp data.'
          IF(idzeff(1).EQV..TRUE. .OR. idzeff(2).EQV..TRUE.)THEN
             mdlni = 1
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE IF(id_ion)THEN
             mdlni = 2
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          ELSE
             WRITE(6,*) 'XX tr_uf_nicheck: Lack of all essential data, Zeff, n_i, n_imp.'
!             WRITE(6,*) 'XX Stopped to correction of profiles by neutrality.'
!             ierr = 2
             mdlni = 9
             WRITE(6,*) '## "MDLNI" is reset to MDLNI= ',mdlni
          END IF
       END IF

    CASE DEFAULT
       WRITE(6,*) 'XX tr_uf_nicheck: the value of "MDLNI" is invalid. MDLNI= ',mdlni
       WRITE(6,*) 'XX Stopped to correction of profiles by charge neutrality.'
       ierr = 2
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
!  id = 1 : complete n_i and n_imp  from Zeff, n_e (and n_bulk)
!  id = 2 : complete n_imp and Zeff from n_e, n_i (and n_bulk)
!  id = 3 : complete n_i and Zeff   from n_e, n_imp (and n_bulk)
!  id = 9 : complete n_i adn Zeff using the assumption n_e = n_i
! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(ikind),               INTENT(IN)  :: ns_ion,ns_imp,id
    REAL(rkind),DIMENSION(1:nsum),INTENT(IN)  :: pzc,pzfc
    REAL(rkind),DIMENSION(1:ntum,1:nrum),       INTENT(INOUT) :: zeffrc
    REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum),INTENT(IN)    :: rnfc
    REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum),INTENT(INOUT) :: rnc

    INTEGER(ikind) :: nsi, ntx

    sumzni(1:ntum,1:nrum)  = 0.d0

    SELECT CASE(id)
    CASE(1) ! Zeff, n_e (and n_bulk) --> n_i, n_imp
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
    CASE(2) ! n_e, n_i (and n_bulk) --> n_imp, Zeff
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

    CASE(3) ! n_e, n_imp (and n_bulk) --> n_i, Zeff
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

    CASE(9) ! n_e --> n_i, Zeff (n_e = n_i)
       rnc(ns_ion,1:ntxmax,1:nrmax+1) = rnc(1,1:ntxmax,1:nrmax+1)
       write(6,*) rnc(ns_ion,1,1:nrmax+1)
       zeffrc(1:ntxmax,1:nrmax+1)     = pzc(ns_ion)

    END SELECT

    ! check for magnitude of correction
    DO nsi = 2, nsum ! ion only: SUM_i(Z_i * n_i)
       sumzni(1:ntxmax,1:nrmax+1) = sumzni(1:ntxmax,1:nrmax+1)     &
                         + pzc(nsi)*rnc(nsi,1:ntxmax,1:nrmax+1)    &
                         + pzfc(nsi)*rnfc(nsi,1:ntxmax,1:nrmax+1)
    END DO

    RETURN
  END SUBROUTINE tr_uf_neutrality

END MODULE trufcalc
