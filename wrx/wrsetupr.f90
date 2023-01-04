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
    REAL(rkind):: deg,arg
    EXTERNAL GUTIME

    ierr=0
    deg=PI/180.D0

    ! --- interactive input for ray tracing ---

    IF(mdlwri.GT.100) THEN
       SELECT CASE(mdlwri)
       CASE(101,102,103,111,112,113) ! interaactive ANGT-ANGP input
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
       CASE(121,122,123,131,132,133) ! interactive Nphi-ANGP input
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
       CASE(141,142,143,151,152,153) ! interactive Nphi-Nz input
          WRITE(6,'(A,I8)') '## nraymax=',nraymax
          DO nray=1,nraymax
30           CONTINUE
             WRITE(6,'(A,I8)') '# nray=',nray
             WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,RNPH,RNZ,MODEW,RNK,UU'
             WRITE(6,'(ES12.4,5F9.2,I4,2F9.2)') &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),RNZIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
             READ(5,*,ERR=20,END=9000) &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),RNZIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
          END DO
       CASE DEFAULT
          WRITE(6,'(A,I4)') 'XX wr_setup_rays: UNKNOWN mdlwri:',mdlwri
          STOP
       END SELECT
    END IF

    ! --- conversion from rnz -> angp ---

    SELECT CASE(mdlwri)
    CASE(41,141) ! RNZ->ANGP positive: poloidal first
       DO nray=1,nraymax
          ANGPIN(NRAY)= ASIN(RNZIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(51,151) ! RNZ->ANGP negative: poloidal first
       DO nray=1,nraymax
          ANGPIN(NRAY)=-ASIN(RNZIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(42,142) ! RNZ->ANGP positive: toroidal first
       DO nray=1,nraymax
          ANGPIN(NRAY)= ASIN(RNZIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGTIN(NRAY)*deg)))/deg
       END DO
    CASE(52,152) ! RNZ->ANGP negative: toroidal first
       DO nray=1,nraymax
          ANGPIN(NRAY)=-ASIN(RNZIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGTIN(NRAY)*deg)))/deg
       END DO
    END SELECT

    ! --- conversion from rnph -> angt ---
    
    SELECT CASE(mdlwri)
    CASE(21,41,51,121,141,151) ! RNPH->ANGT: poloidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGPIN(NRAY)*deg)))/deg
       END DO
    CASE(31,131) ! RNPH->ANGT: absolute angle:  poloidal first
       DO nray=1,nraymax
          ARG=180.D0-ASIN(RNPHIN(NRAY) &
               /(RNKIN(NRAY)*COS((180.D0-ANGPIN(NRAY))*deg)))/deg
       END DO
    CASE(22,42,52,122,142,152) ! RNPH->ANGT: toroidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(32,132) ! RNPH->ANGT: absolue angel: toroidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=180.D0-ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(23,43,53,123,143,153) ! RNPH->ANGT: intuitive
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(33,133) ! RNPH->ANGT: absolute angle: intuitive
       DO nray=1,nraymax
          ANGTIN(NRAY)=180.D0-ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    END SELECT

    RETURN
9000 CONTINUE
    ierr=1
    RETURN
  END SUBROUTINE wr_setup_rays

END MODULE wrsetupr
