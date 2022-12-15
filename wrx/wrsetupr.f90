!   wrsetupr.f90

MODULE wrsetupr

  PRIVATE
  PUBLIC wr_setup_rays

CONTAINS

!   ***** setup initial ray tracing parameters *****
!                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY),
!                  ANGTIN(NRAY),ANGPIN(NRAY),
!                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)

  SUBROUTINE wr_setup_rays(ierr)

    USE wrcomm
    USE wrsub,ONLY: wrcale_i
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nray
    REAL(rkind):: deg
    EXTERNAL GUTIME

    ierr=0
    deg=PI/180.D0

    ! --- interactive input for ray tracing ---

    IF(mdlwri.GT.100) THEN
       SELECT CASE(mdlwri)
       CASE(101,102,103,121,122,123) ! interaactive ANGT-ANGP input
          WRITE(6,'(A,I8)') '## nraymax=',nraymax
          DO NRAY=1,NRAYMAX
10           CONTINUE
             WRITE(6,'(A,I8)') '# nray=',nray
             WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,ANGT,ANGP,MODEW,RNK,UU'
             WRITE(6,'(ES12.4,5F9.2,I4,2F9.2)') &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  ANGTIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
             READ(5,*,ERR=10,END=9000) &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  ANGTIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
          END DO
       CASE(111,112,113) ! interactive Nphi-ANGP input
          WRITE(6,'(A,I8)') '## nraymax=',nraymax
          DO nray=1,nraymax
20           CONTINUE
             WRITE(6,'(A,I8)') '# nray=',nray
             WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,RNPH,ANGP,MODEW,RNK,UU'
             WRITE(6,'(ES12.4,5F9.2,I4,2F9.2)') &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
             READ(5,*,ERR=20,END=9000) &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
          END DO
       CASE DEFAULT
          WRITE(6,'(A,I4)') 'XX wr_setup_rays: UNKNOWN mdlwri:',mdlwri
          STOP
       END SELECT
    END IF

    ! --- conversion from rnph -> angt ---
    
    SELECT CASE(mdlwri)
    CASE(11,111) ! RNPH->ANGT: poloidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/(RNKIN(NRAY)*COS(ANGTIN(NRAY)*deg))) &
                      /deg
       END DO
    CASE(12,112) ! RNPH->ANGT: toroidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(13,113) ! RNPH->ANGT: intuitive
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    END SELECT

    RETURN
9000 CONTINUE
    ierr=1
    RETURN
  END SUBROUTINE wr_setup_rays

END MODULE wrsetupr
