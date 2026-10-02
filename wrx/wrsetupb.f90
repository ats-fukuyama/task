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
    
    IF(mdlwri.GT.100) THEN
       SELECT CASE(mdlwri/10)
       CASE(10,11) ! interaactive ANGT-ANGP input
          WRITE(6,'(A,I8)') '## nraymax=',nraymax
          DO NRAY=1,NRAYMAX
10           CONTINUE
             WRITE(6,'(A,I8)') '# nray=',nray
             WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,ANGT,ANGP,MODEW,RNK,UU'
             WRITE(6,'(ES12.4,5F9.2,I4,2F9.2)') &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  ANGTIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
             WRITE(6,'(12X,4F9.2)') &
                  RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
             READ(5,*,ERR=10,END=9000) &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  ANGTIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY), &
                  RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
          END DO
       CASE(12,13) ! interactive Nph-ANGP input
          WRITE(6,'(A,I8)') '## nraymax=',nraymax
          DO nray=1,nraymax
20           CONTINUE
             WRITE(6,'(A,I8)') '# nray=',nray
             WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,RNPH,ANGP,MODEW,RNK,UU'
             WRITE(6,'(A)')    '#        RCURVA,RCURVB,RBRADA,RBRADB'
             WRITE(6,'(ES12.4,5F9.2,I4,2F9.2)') &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
             WRITE(6,'(12X,4F9.2)') &
                  RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
             READ(5,*,ERR=20,END=9000) &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),ANGPIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY), &
                  RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
          END DO
       CASE(14,15) ! interactive Nph-Nz input
          WRITE(6,'(A,I8)') '## nraymax=',nraymax
          DO nray=1,nraymax
30           CONTINUE
             WRITE(6,'(A,I8)') '# nray=',nray
             WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,RNPH,RNZ,MODEW,RNK,UU'
             WRITE(6,'(A)')    '#        RCURVA,RCURVB,RBRADA,RBRADB'
             WRITE(6,'(ES12.4,5F9.2,I4,2F9.2)') &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),RNZIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
             WRITE(6,'(12X,4F9.2)') &
                  RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
             READ(5,*,ERR=30,END=9000) &
                  RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
                  RNPHIN(NRAY),RNZIN(NRAY), &
                  MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY), &
                  RCURVAIN(NRAY),RCURVBIN(NRAY),RBRADAIN(NRAY),RBRADBIN(NRAY)
          END DO
       CASE DEFAULT
          WRITE(6,'(A,I4)') 'XX wr_setup_rays: UNKNOWN mdlwri:',mdlwri
          STOP
       END SELECT
    END IF
    
    ! --- conversion from rnz -> angp ---

    SELECT CASE(mdlwri)
    CASE(41,43,141,143) ! RNZ->ANGP positive: poloidal first, intuitive,  LFS
       DO nray=1,nraymax
          ANGPIN(NRAY)= ASIN(RNZIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(42,142) ! RNZ->ANGP positive: toroidal first LFS
       DO nray=1,nraymax
          ANGPIN(NRAY)= ASIN(RNZIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGTIN(NRAY)*deg)))/deg
       END DO
    CASE(51,53,151,153) ! RNZ->ANGP negative: poloidal first, intuitive,  HFS
       DO nray=1,nraymax
          ANGPIN(NRAY)=180.D0+ASIN(RNZIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(52,152) ! RNZ->ANGP negative: toroidal first HFS
       DO nray=1,nraymax
          ANGPIN(NRAY)=180.D0+ASIN(RNZIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGTIN(NRAY)*deg)))/deg
       END DO
    END SELECT

    ! --- conversion from rnph -> angt ---
    
    SELECT CASE(MOD(mdlwri,100))
    CASE(21) ! RNPH->ANGT: poloidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGPIN(NRAY)*deg)))/deg
       END DO
    CASE(22,23) ! RNPH->ANGT: toroidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(31) ! RNPH->ANGT: poloidal first: absolute angle
       DO nray=1,nraymax
          ANGTIN(NRAY)=180.D0-ASIN(RNPHIN(NRAY) &
               /(RNKIN(NRAY)*COS(180.D0-ANGPIN(NRAY)*deg))/deg
       END DO
    CASE(32,33) ! RNPH->ANGT: absolue angel: toroidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=180.D0-ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    CASE(41,51) ! RNPH->ANGT: poloidal first
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY) &
               /(RNKIN(NRAY)*COS(ANGPIN(NRAY)*deg)))/deg
       END DO
    CASE(42,43,52,53) ! RNPH->ANGT: intuitive
       DO nray=1,nraymax
          ANGTIN(NRAY)=ASIN(RNPHIN(NRAY)/RNKIN(NRAY))/deg
       END DO
    END SELECT

    RETURN
9000 CONTINUE
    ierr=1
    RETURN
  END SUBROUTINE wr_setup_beams
    
END MODULE wrsetupb
