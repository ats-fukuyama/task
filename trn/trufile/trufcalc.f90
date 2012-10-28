MODULE trufcalc
! ------------------------------------------------------------------------
!   Correction of experimental data, and calculation of associated values
! ------------------------------------------------------------------------

  USE trcomm, ONLY: ikind,rkind,nrum,ntum,nsum,ntxmax,nrmax
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf_complete, &
       sumzni,rnuc,rnfuc,ns_mion,ns_mimp

  REAL(rkind),DIMENSION(1:ntum,1:nrum)        :: zeffruc, sumzni,sumznin
  REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum) :: rnuc,rnfuc

  INTEGER(ikind) :: ns_mion, ns_mimp

CONTAINS

  SUBROUTINE tr_uf_complete
    USE uflist, ONLY: uf0d
    USE truf0d, ONLY: idnm,idnfast,idnmaz,idnfastaz,idzeff
    USE trcomm, ONLY: pau,pzu,pzfu,rnu,rnfu,zeffru, ntxmax,nrmax

    IMPLICIT NONE

    INTEGER(ikind) :: id,nsu,nsi,ierr
    REAL(rkind) :: pa_mion,pz_mion,pa_mimp,pz_mimp

    rnuc(1:nsum,1:ntum,1:nrum)  = 0.d0
    rnfuc(1:nsum,1:ntum,1:nrum) = 0.d0
    zeffruc(1:ntum,1:nrum)      = 0.d0
    
    DO nsu = 1, nsum
       rnuc(nsu,1:ntxmax,1:nrmax+1)  = rnu(nsu,1:ntxmax,1:nrmax+1)
       rnfuc(nsu,1:ntxmax,1:nrmax+1) = rnfu(nsu,1:ntxmax,1:nrmax+1)
    END DO
    zeffruc(1:ntxmax,1:nrmax+1) = zeffru(1:ntxmax,1:nrmax+1)

    ! Identify the main ion and impurity (their 'pa' and 'pz')
    CALL tr_uf_main_ion(pa_mion,pz_mion,pa_mimp,pz_mimp)

    ns_mion = 0
    ns_mimp = 0
    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          IF(pzu(nsu) == pz_mion .AND. pau(nsu) == pa_mion)THEN
             ns_mion = nsi
          ELSE IF(pzu(nsu) == pz_mimp .AND. pau(nsu) == pa_mimp)THEN
             ns_mimp = nsi
          END IF
       END IF
    END DO

    ! the case that the main ion density data can not be found.
    IF(ns_mion == 0)THEN
       DO nsi = 2, nsum
          IF(idnm(nsi) .EQV. .FALSE.)THEN
             idnm(nsi)   = .TRUE.
             idnmaz(nsi) = .TRUE.

             ns_mion = nsi
             pzu(ns_mion) = pz_mion
             pau(ns_mion) = pa_mion
             EXIT
          END IF
          WRITE(6,*) 'XX tr_uf_complete: No blank in UFILE array.'
          ierr = 1
          RETURN
       END DO
    END IF

    ! the case that the main impurity density data can not be found.
    IF(ns_mimp == 0)THEN
       DO nsi = 3, nsum ! nsi = 2: impurity data should not be store.
          IF(idnm(nsi) .EQV. .FALSE.)THEN
             idnm(nsi)   = .TRUE.
             idnmaz(nsi) = .TRUE.

             ns_mimp = nsi
             pzu(ns_mimp) = pz_mimp
             pau(ns_mimp) = pa_mimp
             EXIT
          END IF
       END DO
    END IF
             

    id = 1
    CALL tr_uf_neutrality_check(rnuc,pzu,rnfuc,pzfu,ns_mion,ns_mimp, &
         zeffruc,id)

    RETURN
  END SUBROUTINE tr_uf_complete

! **********************************************************************

  SUBROUTINE tr_uf_main_ion(pa_mion,pz_mion,pa_mimp,pz_mimp)
! ----------------------------------------------------------------------
!   Identify the main ion and main impurity which are re-caluculated
!    to satisfy the charge quasi-neutrality.
!   ** default main ion: deuterium
!   ** default main impurity: carbon
! ----------------------------------------------------------------------
    USE uflist, ONLY: uf0d
    IMPLICIT NONE

!    INTEGER(ikind),INTENT(OUT) :: ierr
    REAL(rkind),   INTENT(OUT) :: pa_mion,pz_mion,pa_mimp,pz_mimp


    pa_mion = uf0d(107)%fr0 ! 'mean' mass number of the plasma ions
    pz_mion = uf0d(108)%fr0 ! 'mean' charge number of the plasma ions

    pa_mimp = uf0d(112)%fr0 ! mass number of the plasma main impurity
    pz_mimp = uf0d(113)%fr0 ! charge number of the plasma main impurity

    WRITE(6,*) 'tr_uf_main_ion:'
    WRITE(6,*) 'PGASA= ',pa_mion, 'PGASZ= ',pz_mion
    WRITE(6,*) 'PIMPA= ',pa_mimp, 'PIMPZ= ',pz_mimp


    ! Determination of main ion which is used for correction of profile.
    IF(uf0d(107)%lex .EQV. .FALSE. .OR. uf0d(108)%lex .EQV. .FALSE.)THEN
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
          WRITE(6,*) 'PGASA= ',uf0d(107)%fr0, 'PGASZ= ',uf0d(108)%fr0
          pa_mion = 2.d0
          pz_mion = 1.d0 ! Deuterium: defalut
!          ierr = 1
          RETURN
       END IF
    END IF

    ! Determination of main impurity which is used for correction of profile.
    IF(uf0d(112)%lex .EQV. .FALSE. .OR. uf0d(113)%lex .EQV. .FALSE.)THEN
       pa_mimp = 12.d0
       pz_mimp = 6.d0 ! Carbon: default
    ELSE
       IF(pa_mimp < 9.d0)THEN
          pa_mimp = 8.d0 ! Beryllium
       ELSE IF(9.d0 <= pa_mimp .AND. pa_mimp < 11.d0)THEN
          pa_mimp = 10.d0 ! Boron
       ELSE IF(11.d0 <= pa_mimp .AND. pa_mimp < 12.d0)THEN
          pa_mimp = 12.d0 ! Carbon
       ELSE
          WRITE(6,*) 'XX tr_uf_comp_check: Unsupported impurity element.'
          WRITE(6,*) 'PIMPA= ',uf0d(112)%fr0, 'PIMPZ= ',uf0d(113)%fr0
          pa_mimp = 12.d0
          pz_mimp = 6.d0 ! Carbon: default
!          ierr = 2
          RETURN
       END IF
    
       IF(pz_mimp < 4.5d0)THEN
          pz_mimp = 4.d0 ! Beryllium
       ELSE IF(4.5d0 <= pz_mimp .AND. pz_mimp < 5.5d0)THEN
          pz_mimp = 5.d0 ! Boron
       ELSE IF(5.5d0 <= pz_mimp .AND. pz_mimp < 6.5d0)THEN
          pz_mimp = 6.d0 ! Carbon
       ELSE
          WRITE(6,*) 'XX tr_uf_comp_check: Unsupported impurity element.'
          WRITE(6,*) 'PIMPA= ',uf0d(112)%fr0, 'PIMPZ= ',uf0d(113)%fr0
          pa_mimp = 12.d0
          pz_mimp = 6.d0 ! Carbon: default
!          ierr = 2
          RETURN
       END IF
    END IF

    RETURN
  END SUBROUTINE tr_uf_main_ion


!!$  SUBROUTINE tr_uf_comp_table
!!$! ------------------------------------------------------------------------
!!$!   ???
!!$! ------------------------------------------------------------------------
!!$    USE truf0d, ONLY: idnmaz
!!$    IMPLICIT NONE
!!$
!!$    INTEGER(ikind) :: nsi,nsc,nslex
!!$    INTEGER(ikind),DIMENSION(1:9) :: nsi_nslex
!!$
!!$    nslex = 0
!!$    DO nsi = 2, nsum-1 ! NM1 ~ NM8
!!$       IF(idnmaz(nsi-1)) nslex = nslex + 1
!!$       IF(nslex /= nsi) THEN
!!$          ! ion only
!!$          DO nsc = nsi, nsum-1
!!$             rnu(nsc,1:ntxmax,1:nrmax+1) = rnu(nsc+1,1:ntxmax,1:nrmax+1)
!!$          END DO
!!$          rnu(nsum,1:ntxmax,1:nrmax+1) = 0.d0
!!$
!!$          nslex = nsi
!!$       END IF
!!$    END DO
!!$    
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_uf_comp_table

! **********************************************************************

  SUBROUTINE tr_uf_neutrality_check(rnc,pzc,rnfc,pzfc,ns_ion,ns_imp, &
                                    zeffrc,id)
! ----------------------------------------------------------------------
!  < NOTATION >
!  n_e    : particle density of electron
!  n_i    : particle density of main ion
!  n_bulk : particle density of other ions (not impurity)
!  n_imp  : particle density of impurity
!
!  id = 1 : complete n_i and n_imp  from Zeff, n_e (and n_bulk)
!  id = 2 : complete n_imp and Zeff from n_e, n_i (and n_bulk)
!  id = 3 : complete n_i and Zeff   from n_e, n_imp (and n_bulk)
! ----------------------------------------------------------------------
    USE trcomm,ONLY: nrmax,ntxmax,rnu

    IMPLICIT NONE

    INTEGER(ikind),               INTENT(IN)  :: ns_ion,ns_imp,id
    REAL(rkind),DIMENSION(1:nsum),INTENT(IN)  :: pzc,pzfc
    REAL(rkind),DIMENSION(1:ntum,1:nrum),       INTENT(IN)    :: zeffrc
    REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum),INTENT(IN)    :: rnfc
    REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum),INTENT(INOUT) :: rnc

    INTEGER(ikind) :: nsi, ntx

    sumzni(1:ntum,1:nrum)  = 0.d0
    sumznin(1:ntum,1:nrum) = 0.d0

    SELECT CASE(id)
    CASE(1)
       DO ntx = 1, ntxmax
          rnc(ns_ion,ntx,1:nrmax+1) = (pzc(ns_imp)-zeffrc(ntx,1:nrmax+1)) &
                                      *rnc(1,ntx,1:nrmax+1)
          DO nsi = 2, nsum
             IF(nsi /= ns_ion .AND. nsi /= ns_imp)THEN
                rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)     &
                     - (pzc(ns_imp)-pzc(nsi))*pzc(nsi)*rnc(nsi,ntx,1:nrmax+1)
             END IF
             rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)        &
                  - (pzc(ns_imp)-pzfc(nsi))*pzfc(nsi)*rnfc(nsi,ntx,1:nrmax+1)
          END DO

          rnc(ns_ion,ntx,1:nrmax+1) = rnc(ns_ion,ntx,1:nrmax+1)           &
                                     /(pzc(ns_ion)*(pzc(ns_imp)-pzc(ns_ion)))
       END DO
    CASE(2)
       CONTINUE
    CASE(3)
       CONTINUE
    END SELECT


    DO nsi = 2, nsum ! ion only: SUM_i(Z_i * n_i)
       sumzni(1:ntxmax,1:nrmax+1) = sumzni(1:ntxmax,1:nrmax+1)     &
                         + pzc(nsi)*rnc(nsi,1:ntxmax,1:nrmax+1)    &
                         + pzfc(nsi)*rnfc(nsi,1:ntxmax,1:nrmax+1)
       write(6,*) nsi+1, pzc(nsi+1), pzfc(nsi+1)
    END DO


    RETURN
  END SUBROUTINE tr_uf_neutrality_check

END MODULE trufcalc
