MODULE trufcalc

  USE trcomm, ONLY: ikind,rkind,nrum,ntum,nsum,ntxmax,nrmax
  IMPLICIT NONE

  PRIVATE
  PUBLIC tr_uf_complete, &
       sumzni,rniu,rnfiu

  REAL(rkind),DIMENSION(1:ntum,1:nrum)        :: sumzni
  REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum) :: rnuc,rnfuc

CONTAINS

  SUBROUTINE tr_uf_complete
    USE uflist, ONLY: uf0d
    USE trufsub, ONLY: idnm,idnfast,idnmaz,idnfastaz,idzeff
    USE trcomm, ONLY: pau,pzu,rnu,rnfu,ntxmax,nrmax

    IMPLICIT NONE

    INTEGER(ikind) :: nsu,nsi,ns_mion,ns_mimp,ierr

    
    DO nsu = 1, nsum
       rnuc(nsu,1:ntxmax,1:nrmax+1)  = rnu(nsu,1:ntxmax,1:nrmax+1)
       rnfuc(nsu,1:ntxmax,1:nrmax+1) = rnfu(nsu,1:ntxmax,1:nrmax+1)
    END DO

    CALL tr_uf_comp_check


    ns_mion = 0
    ns_mimp = 0
    DO nsi = 2, nsum
       IF(idnm(nsi))THEN
          IF(pzu(nsu) == pz_mion .AND. pau(nsu) == pa_mion)THEN
             ns_mion = nsu
          ELSE IF(pzu(nsu) == pz_mimp .AND. pau(nsu) == pa_mimp)THEN
             ns_mimp = nsu
          END IF
       END IF
    END DO
    

    id = 0
    CALL tr_uf_neutrality_check(pzi,id,ierr)

    RETURN
  END SUBROUTINE tr_uf_complete

! **********************************************************************

  SUBROUTINE tr_uf_comp_check(pa_mion,pz_mion,pa_mimp,pz_mimp,ierr)
    USE uflist, ONLY: uf0d
    IMPLICIT NONE

    INTEGER(ikind),INTENT(OUT) :: ierr
    REAL(rkind),   INTENT(OUT) :: pa_mion,pz_mion,pa_mimp,pz_mimp

    ierr = 0

    pa_mion = uf0d(107)%fr0 ! 'mean' mass number of the plasma ions
    pz_mion = uf0d(108)%fr0 ! 'mean' charge number of the plasma ions

    pa_mimp = uf0d(112)%fr0 ! mass number of the plasma main impurity
    pz_mimp = uf0d(113)%fr0 ! charge number of the plasma main impurity

    WRITE(6,*) 'PGASA= ',pa_mion, 'PGASZ= ',pz_mion
    WRITE(6,*) 'PIMPA= ',pa_mimp, 'PIMPZ= ',pz_mimp


    ! Determination of main ion which is used for correction of profile.
    IF(uf0d(107)%lex /= .TRUE. .OR. uf0d(108)%lex /= .TRUE.)THEN
       pa_mion = 2.d0
       pz_mion = 1.d0 ! Deuterium: defalut
    ELSE IF
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
          ierr = 1
          RETURN
       END IF
    END IF

    ! Determination of main impurity which is used for correction of profile.
    IF(uf0d(112)%lex /= .TRUE. .OR. uf0d(113)%lex /= .TRUE.)THEN
       pa_mimp = 12.d0
       pz_mimp = 6.d0 ! Carbon: default
    ELSE IF
       IF(pa_mimp < 9.d0)THEN
          pa_mimp = 8.d0 ! Beryllium
       ELSE IF(9.d0 <= pa_mimp .AND. pa_mimp < 11.d0)THEN
          pa_mimp = 10.d0 ! Boron
       ELSE IF(11.d0 <= pa_mimp .AND. pa_mimp < 12.d0)THEN
          pa_mimp = 12.d0 ! Carbon
       ELSE
          WRITE(6,*) 'XX tr_uf_comp_check: Unsupported impurity element.'
          WRITE(6,*) 'PIMPA= ',uf0d(112)%fr0, 'PIMPZ= ',uf0d(113)%fr0
          ierr = 1
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
          ierr = 1
          RETURN
       END IF
    END IF

    RETURN
  END SUBROUTINE tr_uf_comp_check


  SUBROUTINE tr_uf_comp_table
    USE trufile, ONLY: idnmaz
    IMPLICIT NONE

    INTEGER(ikind) :: nsi,nsc,nslex
    INTEGER(ikind),DIMENSION(1:9) :: nsi_nslex

    nslex = 0
    DO nsi = 2, nsum-1 ! NM1 ~ NM8
       IF(idnmaz(nsi-1)) nslex = nslex + 1
       IF(nslex /= nsi) THEN
          ! ion only
          DO nsc = nsi, nsum-1
             rnu(nsc,1:ntxmax,1:nrmax+1) = rnu(nsc+1,1:ntxmax,1:nrmax+1)
          END DO
          rnu(nsum,1:ntxmax,1:nrmax+1) = 0.d0

          nslex = nsi
       END IF
    END DO
    

    RETURN
  END SUBROUTINE tr_uf_comp_table

!!$  SUBROUTINE tr_uf_comp_nimp
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_uf_comp_nimp
!!$
!!$  SUBROUTINE tr_uf_comp_zeff
!!$
!!$    RETURN
!!$  END SUBROUTINE tr_uf_comp_zeff

! **********************************************************************

  SUBROUTINE tr_uf_neutrality_check(pzi,id,ierr)
! ----------------------------------------------------------------------
!
!  id = 1 : from zeff make n_x
!  id = 2 : from n_i make n_x and zeff
! ----------------------------------------------------------------------
    USE trcomm,ONLY: nrmax,ntxmax,rnu

    IMPLICIT NONE

    INTEGER(ikind),             INTENT(IN)  :: id
    REAL(rkind),                INTENT(IN)  :: pzi
    INTEGER(ikind),             INTENT(OUT) :: ierr

    INTEGER(ikind) :: nsi

    ierr = 0
    sumzni(1:ntum,1:nrum) = 0.d0

    DO nsi = 1, 9 ! ion only: SUM_i(Z_i * n_i)
       sumzni(1:ntxmax,1:nrmax+1) = sumzni(1:ntxmax,1:nrmax+1)           &
                         + pziu(nsi+1)*rniu(nsi+1,1:ntxmax,1:nrmax+1)    &
                         + pzfiu(nsi+1)*rnfiu(nsi+1,1:ntxmax,1:nrmax+1)
       write(6,*) nsi+1, pziu(nsi+1), pzfiu(nsi+1)
    END DO

    RETURN
  END SUBROUTINE tr_uf_neutrality_check

END MODULE trufcalc
