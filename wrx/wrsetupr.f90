!   wrsetupr.f90

MODULE wrsetupr

  PRIVATE
  PUBLIC wr_setup_rays

CONTAINS

!   ***** setup initial ray tracing parameters *****

  SUBROUTINE wr_setup_rays(ierr)

    USE wrcomm
    USE wrsub,ONLY: wrnwtn,wrcale_i
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: nray
    REAL(rkind):: deg,factorp,factort
    EXTERNAL GUTIME

    ierr=0
    deg=PI/180.D0

    ! --- interactive input for ray tracing ---
    
    SELECT CASE(mdlwri)
    CASE(0,1) ! interaactive ANG input
       WRITE(6,'(A,I8)') '## nraymax=',nraymax
       DO NRAY=1,NRAYMAX
10        CONTINUE
          WRITE(6,'(A,I8)') '# nray=',nray
          WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,ANGP,ANGT,RNK,UU,MODEW'
          WRITE(6,'(ES12.4,7F9.2,I4)') &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               ANGPIN(NRAY),ANGTIN(NRAY),RNKIN(NRAY),UUIN(NRAY),MODEWIN(NRAY)
          READ(5,*,ERR=10,END=9000) &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               ANGPIN(NRAY),ANGTIN(NRAY),RNKIN(NRAY),UUIN(NRAY),MODEWIN(NRAY)
       END DO
    CASE(2,3) ! interactive RN input
       WRITE(6,'(A,I8)') '## nraymax=',nraymax
       DO nray=1,nraymax
20        CONTINUE
          WRITE(6,'(A,I8)') '# nray=',nray
          WRITE(6,'(A)')    '# input: RF,RP,ZP,PHI,RNKP,RNKT,RNK,UU,MODEW'
          WRITE(6,'(ES12.4,5F10.3,I4,F10.3)') &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               RNKPIN(NRAY),RNKTIN(NRAY),MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
          READ(5,*,ERR=20,END=9000) &
               RFIN(NRAY),RPIN(NRAY),ZPIN(NRAY),PHIIN(NRAY), &
               RNKPIN(NRAY),RNKTIN(NRAY),MODEWIN(NRAY),RNKIN(NRAY),UUIN(NRAY)
       END DO
    END SELECT

    ! --- conversion between ang <-> rnk ---
    
    SELECT CASE(mdlwri)
    CASE(0,100) ! ANG->RN: poloidal first
       DO nray=1,nraymax
          RNKPIN(NRAY)=RNKIN(NRAY)*SIN(ANGPIN(NRAY)*deg)
          RNKTIN(NRAY)=RNKIN(NRAY)*SIN(ANGTIN(NRAY)*deg)*COS(ANGPIN(NRAY)*deg)
       END DO
    CASE(1,101) ! ANG->RN: toroidal first
       DO nray=1,nraymax
          RNKTIN(NRAY)=RNKIN(NRAY)*SIN(ANGTIN(NRAY)*deg)
          RNKPIN(NRAY)=RNKIN(NRAY)*SIN(ANGPIN(NRAY)*deg)*COS(ANGPIN(NRAY)*deg)
       END DO
    CASE(2,102) ! RN->ANG: poloidal first
       DO nray=1,nraymax
          factorp=RNKPIN(NRAY)/RNKIN(NRAY)
          IF(ABS(factorp).GT.1.D0) MODEWIN(NRAY)=-1
          factort=RNKIN(NRAY)**2-RNKPIN(NRAY)**2
          IF(RNKTIN(NRAY)**2.GT.factort) MODEWIN(NRAY)=-2
          IF(MODEWIN(NRAY).GE.0) THEN
             ANGPIN(NRAY)=ASIN(factorp)
             ANGTIN(NRAY)=ASIN(RNKTIN(NRAY)/SQRT(factort))
          ELSE
             WRITE(6,'(A)') &
                'XX wr_setup_rays: nray,modewin,factorp,factort,rnkp,rnkt,rnk'
             WRITE(6,'(A,2I4,5ES12.4)') &
                  '      ',nray,modewin(nray),factorp,factort, &
                  rnkpin(nray),rnktin(nray),rnkin(nray)
          END IF
       END DO
    CASE(3,103) ! RN->ANG: toroidal first
       DO nray=1,nraymax
          factort=RNKTIN(NRAY)/RNKIN(NRAY)
          IF(ABS(factort).GT.1.D0) MODEWIN(NRAY)=-3
          factorp=RNKIN(NRAY)**2-RNKTIN(NRAY)**2
          IF(RNKPIN(NRAY)**2.GT.factorp) MODEWIN(NRAY)=-4
          IF(MODEWIN(NRAY).GE.0) THEN
             ANGTIN(NRAY)=ASIN(factort)
             ANGPIN(NRAY)=ASIN(RNKPIN(NRAY)/SQRT(factorp))
          ELSE
             WRITE(6,'(A)') &
                'XX wr_setup_rays: nray,modewin,factort,factorp,rnkt,rnkp,rnk'
             WRITE(6,'(A,2I4,5ES12.4)') &
                  '      ',nray,modewin(nray),factort,factorp, &
                  rnktin(nray),rnkpin(nray),rnkin(nray)
          END IF
       END DO
       
    END SELECT
    RETURN
9000 CONTINUE
    ierr=9000
    RETURN
  END SUBROUTINE wr_setup_rays

END MODULE wrsetupr
