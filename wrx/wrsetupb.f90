!   wrsetupb.f

MODULE wrsetupb

  PRIVATE
  PUBLIC wr_setup_beams

CONTAINS

!     ***** Beam tracing module *****

  SUBROUTINE wr_setup_beams(ierr)

    USE wrcomm
    USE wrsub,ONLY: wrcale
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nray
    REAL(rkind):: deg,factort
    EXTERNAL GUTIME

    ierr=0
    deg=PI/180.D0

    ! --- interactive input for beam tracing ---
    
    SELECT CASE(mdlwri)
    CASE(0,1) ! interaactive ANG input
       WRITE(6,'(A,I8)') '## nraymax=',nraymax
       DO nray=1,nraymax
10        CONTINUE
          WRITE(6,'(A,I8)') '# nray=',nray
          WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,ANGT,ANGP,RNK,UU,MODEW'
          WRITE(6,'(A)')    '#        RCURVA,RCURVB,RBRADA,RBRADB'
          WRITE(6,'(ES12.4,7F9.2,I4)') &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               ANGTIN(NRAY),ANGPIN(NRAY),RNKIN(NRAY),UUIN(NRAY),MODEWIN(NRAY)
          WRITE(6,'(12X,4F9.2)') &
               RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
          READ(5,*,ERR=10,END=9000) &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               ANGTIN(NRAY),ANGPIN(NRAY),RNKIN(NRAY),UUIN(NRAY),MODEWIN(NRAY),&
               RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
       END DO
    CASE(2,3) ! interactive RN input
       WRITE(6,'(A,I8)') '## nraymax=',nraymax
       DO nray=1,nraymax
20        CONTINUE
          WRITE(6,'(A,I8)') '# nray=',nray
          WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,RNPN,ANGP,RNK,UU,MODEW'
          WRITE(6,'(A)')    '#        RCURVA,RCURVB,RBRADA,RBRADB'
          WRITE(6,'(ES12.4,7F9.2,I4)') &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               RNPHIN(NRAY),ANGPIN(NRAY),RNKIN(NRAY),UUIN(NRAY),MODEWIN(NRAY)
          WRITE(6,'(12X,4F9.2)') &
               RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
          READ(5,*,ERR=20,END=9000) &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               RNPHIN(NRAY),ANGPIN(NRAY),RNKIN(NRAY),UUIN(NRAY),MODEWIN(NRAY),&
               RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
       END DO
    END SELECT

    ! --- conversion between ang <-> rnk ---
    
    SELECT CASE(mdlwri)
    CASE(0,100) ! ANG->RNPH: poloidal first
       DO nray=1,nraymax
          RNPHIN(NRAY)=RNKIN(NRAY)*SIN(ANGTIN(NRAY)*deg)*COS(ANGPIN(NRAY)*deg)
       END DO
    CASE(1,101) ! ANG->RNPH: toroidal first
       DO nray=1,nraymax
          RNPHIN(NRAY)=RNKIN(NRAY)*SIN(ANGTIN(NRAY)*deg)
       END DO
    CASE(2,102) ! RNPH->ANG: poloidal first
       DO nray=1,nraymax
          factort=RNPHIN(NRAY)/(RNKIN(NRAY)*COS(ANGPIN(NRAY)*deg))
          IF(factort.GT.1.D0) THEN
             RNKIN(NRAY)=RNPHIN(NRAY)/COS(ANGPIN(NRAY)*deg)
             factort=1.D0
          END IF
          ANGTIN(NRAY)=ASIN(factort)/deg
       END DO
    CASE(3,103) ! RNPH->ANG: toroidal first
       DO nray=1,nraymax
          factort=RNPHIN(NRAY)/RNKIN(NRAY)
          IF(factort.GT.1.D0) THEN
             RNKIN(NRAY)=RNPHIN(NRAY)
             factort=1.D0
          END IF
          ANGTIN(NRAY)=ASIN(factort)/deg
       END DO
       
    END SELECT
    RETURN

9000 CONTINUE
    ierr=9000
    RETURN
  END SUBROUTINE wr_setup_beams
    
END MODULE wrsetupb
