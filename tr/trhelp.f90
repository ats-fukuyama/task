!     ***********************************************************

!           HELP MESSAGE CONTROL SUBROUTINE

!     ***********************************************************

      SUBROUTINE TRHELP(KID)

      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: KID

      IF(KID.EQ.'M') CALL TRHLPM
      IF(KID.EQ.'G') CALL TRHLPG
      IF(KID.EQ.'W') CALL TRHLPW

      RETURN
      END SUBROUTINE TRHELP

!     ***********************************************************

!           HELP MESSAGE FOR MENU

!     ***********************************************************

      SUBROUTINE TRHLPM

      USE libchar
      IMPLICIT NONE
      INTEGER      :: IEND
      CHARACTER(LEN=1):: KID

    1 WRITE(6,601)
  601 FORMAT( &
     &' ','# INPUT FIRST CHARACTER OF MENU ITEMS:'/                     &
     &' ',' (R)UN       : START NEW TRANSPORT CALCULATION'/             &
     &' ','               AFTER PARAMETER INPUT.'/                      &
     &' ',' (C)ONTINUE  : CONTINUE CALCULATION WITH SAME PARAMETERS.'/  &
     &' ',' (P)ARAMETER : CHANGE PARAMETERS AND CONTINUE CALCULATION.'/ &
     &' ',' (G)RAPH     : GRAPHIC PRESENTATION OF TRANSPORT DATA.'/     &
     &' ',' (W)RITE     : PRINT DETAIL OF TRANSPORT DATA.')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,602)
  602 FORMAT( &
     &' ',' (S)AVE      : SAVE TRANSPORT DATA OF PRESENT CALCULATION'/ &
     &' ','               INTO FILE.'/                                 &
     &' ',' (L)OAD      : LOAD TRANSPORT DATA FROM FILE AND CONTINUE'/ &
     &' ','               OLD CALCULATION AFTER PARAMETER INPUT.'/     &
     &' ',' (H)ELP      : SHOW HELP MESSAGE WHICH YOU ARE READING.'/   &
     &' ',' (E)ND       : TERMINATE TRANSPORT CODE.'/                  &
     &' ','# IF YOU NEED LIST OF INPUT PARAMETERS THROUGH NAMELIST ',  &
     &    '&TR, INPUT P.'/                                             &
     &' ','   ANY OTHER CHARACTER BRINGS YOU BACK TO INITIAL MENU.')
      READ(5,'(A1)',ERR=1,END=900) KID
      CALL toupper(KID)
      IF(KID.EQ.'P') CALL TRHLPP
  900 RETURN
      END SUBROUTINE TRHLPM

!     ***********************************************************

!           HELP MESSAGE FOR GRAPHIC MENU

!     ***********************************************************

      SUBROUTINE TRHLPG

      IMPLICIT NONE
      INTEGER :: IEND

      WRITE(6,601)
  601 FORMAT( &
     &' ','          ***  GRAPHIC INFORMATION  ***'/                &
     &' ','# INPUT GRAPH TYPE OR FIRST CHARACTER OF MENU ITEMS:'/   &
     &' ',' (R1)    : RADIAL PROFILE OF N AND T.'/                  &
     &' ',' (R2)    : RADIAL PROFILE OF W AND BETA.'/               &
     &' ',' (R3)    : RADIAL PROFILE OF POWER.'/                    &
     &' ',' (R4)    : RADIAL PROFILE OF Q,EZ,J,S.'/                 &
     &' ',' (R5)    : RADIAL PROFILE OF ZEFF,IMPURITY.'/            &
     &' ',' (R6)    : RADIAL PROFILE OF PIN,SSIN,PELLET.'/          &
     &' ',' (R7)    : RADIAL PROFILE OF ETA,AK.'/                   &
     &' ',' (R8)    : RADIAL PROFILE OF DEDW,DIDW,ETAC,ETAI,LN,LT'/ &
     &' ',' (R9)    : RADIAL PROFILE OF AD,AV,AVK,TAUB,TAUF.')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,602)
  602 FORMAT( &
     &' ',' (G1/P1) : RADIAL/TEMPORAL PROFILE OF NE,ND,NT,NA.'/     &
     &' ',' (G2/P2) : RADIAL/TEMPORAL PROFILE OF TE,TD,TT,TA.'/     &
     &' ',' (G3/P3) : RADIAL/TEMPORAL PROFILE OF Q,JTOT,EZOH,JOH.'/ &
     &' ',' (G4/P4) : RADIAL/TEMPORAL PROFILE OF JNB,JBS,WB,WF.'/   &
     &' ',' (G5/P5) : RADIAL/TEMPORAL PROFILE OF PIN,POH,PNB,PNF.'/ &
     &' ',' (T1)    : TIME EVOLUTION OF T0,I,PIN,POUT.'/            &
     &' ',' (T2)    : TIME EVOLUTION OF W,T0,TAV,ANE.'/             &
     &' ',' (T3)    : TIME EVOLUTION OF TAUE,QF,BETA.'/             &
     &' ',' (T4)    : TIME EVOLUTION OF VLOOP,LI,Q0,RQ1.')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,603)
  603 FORMAT( &
     &' ',' (C)     : CLEAR GRAPHIC DATA GVR AND GVT.'/             &
     &' ',' (S)     : SAVE GRAPHIC DATA.'/                          &
     &' ',' (L)     : LOAD GRAPHIC DATA.'/                          &
     &' ',' (Q)     : CHANGE GRAPH SCALE INQUIRE MODE.'/            &
     &' ',' (H)ELP  : SHOW HELP MESSAGE WHICH YOU ARE READING.'/    &
     &' ',' (E)ND   : GO BACK TO MAIN MENU.')

      RETURN
      END SUBROUTINE TRHLPG

!     ***********************************************************

!           HELP MESSAGE FOR WRITE MENU

!     ***********************************************************

      SUBROUTINE TRHLPW

      IMPLICIT NONE

      WRITE(6,601)
  601 FORMAT( &
     &' ','          ***  TR INFORMATION  ***'/                     &
     &' ','# INPUT PRINT TYPE OR FIRST CHARACTER OF MENU ITEMS:'/   &
     &' ',' (1)    : DETAIL OF GLOBAL DATA.'/                       &
     &' ',' (2)    : Q PROFILE.'/                                   &
     &' ',' (3)    : USED CPU TIME.'/                               &
     &' ',' (4)    : MIN MAX OF GVT DATA.'/                         &
     &' ',' (5)    : SHOW GRAPHIC INFORMATION (NGR,NGT)'/           &
     &' ',' (H)ELP : SHOW HELP MESSAGE WHICH YOU ARE READING.'/     &
     &' ',' (E)ND  : GO BACK TO MAIN MENU.')

      RETURN
      END SUBROUTINE TRHLPW

!     ***********************************************************

!           LIST OF INPUT PARAMETERS

!     ***********************************************************

      SUBROUTINE TRHLPP

      IMPLICIT NONE
      INTEGER :: IEND

      WRITE(6,601)
  601 FORMAT( &
     &' ','# INPUT PARAMETERS THROUGH NAMELIST &TR ##'/             &
     &' ','      ==== DEVICE PARAMETERS ===='/                      &
     &' ','  RR     : PLASMA MAJOR RADIUS (M)'/                     &
     &' ','  RA     : PLASMA MINOR RADIUS (M)'/                     &
     &' ','  RKAP   : ELIPTICITY OF POLOIDAL CROSS SECTION'/        &
     &' ','  BB     : TOROIDAL MAGNETIC FIELD ON PLASMA AXIS (T)'/  &
     &' ','  RIPS   : INITIAL VALUE OF PLASMA CURRENT (MA)'/        &
     &' ','  RIPE   : FINAL VALUE OF PLASMA CURRENT (MA)')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,603)
  603 FORMAT( &
     &' ','      ==== PLASMA PARAMETERS ===='/                          &
     &' ','  NSMAX  : NUMBER OF PARTICLE SPECIES (NS=1:ELECTRON)'/      &
     &' ','  PA(IS) : ATOMIC NUMBER'/                                   &
     &' ','  PZ(IS) : CHARGE NUMBER'/                                   &
     &' ','  PN(IS) : INITIAL NUMBER DENSITY ON AXIS (1.E20 M**-3)'/    &
     &' ','  PNS(IS): INITIAL NUMBER DENSITY ON SURFACE (1.E20 M**-3)'/ &
     &' ','  PT(IS) : INITIAL TEMPERATURE ON AXIS (KEV)'/               &
     &' ','  PTS(IS): INITIAL TEMPERATURE ON SURFACE (KEV)')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,604)
  604 FORMAT( &
     &' ','  PROFN*: PROFILE PARAMETER OF INITIAL DENSITY'/         &
     &' ','  PROFT*: PROFILE PARAMETER OF INITIAL TEMPERATURE'/     &
     &' ','  PROFU*: PROFILE PARAMETER OF INITIAL CURRENT DENSITY'/ &
     &' ','  PROFJ*: PROFILE PARAMETER OF NEUTRAL DENSITY')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,605)
  605 FORMAT( &
     &' ','      ==== IMPURITY PARAMETERS ===='/                        &
     &' ','  PNC    : CARBON DENSITY FACTOR'/                           &
     &' ','  PNFE   : IRON DENSITY FACTOR'/                             &
     &' ','              COMPARED WITH ITER PHYSICS DESIGN GUIDELINE'/  &
     &' ','  PNNU   : NEUTRAL NUMBER DENSITY ON AXIS (1.E20 M**-3)'/    &
     &' ','  PNNUS  :                        ON SURFACE (1.E20 M**-3)'/ &
     &' ','      ==== TRANSPORT PARAMETERS ===='/                       &
     &' ','  AD0    : INWARD PARTICLE DIFFUSION COEFFICIENT (M2/S)'/    &
     &' ','  CNP    : COEFFICIENT FOR NEOCLASICAL PARTICLE DIFFUSION'/  &
     &' ','  CNH    : COEFFICIENT FOR NEOCLASICAL HEAT DIFFUSION'/      &
     &' ','  CNP    : COEFFICIENT FOR TURBULENT PARTICLE DIFFUSION'/    &
     &' ','  CNH    : COEFFICIENT FOR TURBULENT HEAT DIFFUSION'/        &
     &' ','  CDW    : COEFFICIENTS OF DW MODEL')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,606)
  606 FORMAT( &
     &' ','  ALP(1) : RAIDUS REDUCTION FACTOR'/                     &
     &' ','  ALP(2) : MASS WEIGHTING FACTOR OF AD'/                 &
     &' ','  ALP(3) : CHARGE WEIGHTING FACTOR OF ADL'/              &
     &' ','  MDLKAI : DW TRANSPORT MODEL TYPE'/                     &
     &' ','           0:GYRO-BOHM  1:ETAC=1  2:ROMANELLI 3 :ITOH'/  &
     &' ','  MDLETA : RESISTIVITY MODEL TYPE'/                      &
     &' ','           0:CLASSICAL 1:NEOCLASSICAL'/                  &
     &' ','  MDLAD  : PARTICLE DIFFUSION MODEL TYPE'/               &
     &' ','           0:NO PARTICL TRANSPORT 1:CONSTANT D'/         &
     &' ','  MDLAVK : HEAT PINCH MODEL TYPE'/                       &
     &' ','           0:NO HEAT PINCH')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,607)
  607 FORMAT( &
     &' ','      ==== CONTROL PARAMETERS ===='/                     &
     &' ','  NRMAX  : NUMBER OF RADIAL MESH POINTS'/                &
     &' ','  DT     : SIZE OF TIME STEP'/                           &
     &' ','  NTMAX  : NUMBER OF TIME STEP'/                         &
     &' ','  NTSTEP : INTERVAL OF SNAP DATA PRINT'/                 &
     &' ','  NGRSTP : INTERVAL OF RADIAL PROFILE SAVE'/             &
     &' ','  NGTSTP : INTERVAL OF TIME EVOLUTION SAVE'/             &
     &' ','  EPSLTR : CONVERGENCE CRITERION OF ITERATION'/          &
     &' ','  LMAXTR : MAXIMUM COUNT OF ITERATION')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,608)
  608 FORMAT( &
     &' ','      ==== SAWTOOTH PARAMETERS ===='/                    &
     &' ','  TPRST  : SAWTOOTH PERIOD (S)'/                         &
     &' ','  MDLST  : SAWTOOTH MODEL TYPE'/                         &
     &' ','           0:OFF  1:ON'/                                 &
     &' ','      ==== FUSION REACTION PARAMETERS ===='/             &
     &' ','  MDLNF  : FUSION REACTION MODEL TYPE'/                  &
     &' ','           0:OFF  1:ON'/                                 &
     &' ','      ==== NBI HEATING PARAMETERS ===='/                 &
     &' ','  PNBTOT : NBI TOTAL INPUT POWER (MW)'/                  &
     &' ','  PNBR0  : RADIAL POSITION OF NBI POWER DEPOSITION (M)'/ &
     &' ','  PNBRW  : RADIAL WIDTH OF NBI POWER DEPOSITION (M)'/    &
     &' ','  PNBVY  : VERTICAL POSITION OF NBI (M)'/                &
     &' ','  PNBENG : NBI BEAM ENERGY (KEV)'/                       &
     &' ','  PNBRT  : TANGENTIAL RADIUS OF NBI BEAM (M)'/           &
     &' ','  MDLNB  : NBI MODEL TYPE'/                              &
     &' ','           0:OFF  1:GAUSSIAN  2:PENCIL BEAM')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,609)
  609 FORMAT( &
     &' ','      ==== RF HEATING PARAMETERS ===='/                  &
     &' ','  PRFTOT : RF INPUT POWER (MW)'/                         &
     &' ','  PRFR0  : RADIAL POSITION OF RF POWER DEPOSITION (M)'/  &
     &' ','  PRFRW  : RADIAL WIDTH OF RF POWER DEPOSITION (M)'/     &
     &' ','  DIVRF  : DIVISION OF RF POWER TO PARTICLE SPECIES'/    &
     &' ','      ==== CURRENT DRIVE PARAMETERS ===='/               &
     &' ','  PNBCD  : CURRENT DRIVE EFFICIENCY FACTOR OF NBI'/      &
     &' ','  PRFCD  : CURRENT DRIVE EFFICIENCY FACTOR OF RF'/       &
     &' ','  PBSCD  : BOOTSTRAP CURRENT DRIVE FACTOR'/              &
     &' ','  MDLCD  : CURRENT DRIVE MODE TYPE'/                     &
     &' ','           0:CONSTANT IP  1:NO OHMIC')
      CALL TRHLPQ(IEND)
      IF(IEND.EQ.1) RETURN

      WRITE(6,610)
  610 FORMAT( &
     &' ','       ==== PELLET PARAMETERS ===='/                     &
     &' ','  MDLPEL : PELLET INJECTION MODEL TYPE'/                 &
     &' ','           0:OFF  1:GAUSSIAN  2:NAKAMURA  3:HO'/         &
     &' ','  PELTOT : TOTAL NUMBER OF PARTICLES IN PELLET'/         &
     &' ','  PELR0  : RADIAL POSITION OF PELLET DEPOSITION (M)'/    &
     &' ','  PELRW  : RADIAL WIDTH OF PELLET DEPOSITION (M)'/       &
     &' ','  PELTIM : TIME FOR PELLET TO BE INJECTED'/              &
     &' ','  PELRAD : RADIUS OF PELLET (M)'/                        &
     &' ','  PELVEL : PELLET INJECTION VELOCITY (M/S)'/             &
     &' ','  PELPAT : PARTICLE RATIO IN PELLET')

      RETURN
      END SUBROUTINE TRHLPP

!     ***********************************************************

!           INQUIRY IN HELP MESSAGE

!     ***********************************************************

      SUBROUTINE TRHLPQ(IEND)

      USE libchar
      IMPLICIT NONE
      INTEGER, INTENT(OUT):: IEND
      CHARACTER(LEN=1) :: KID

      IEND=0
    1 WRITE(6,*) '# INPUT "E" KEY TO GO BACK TO MENU, OTHERWISE CONTINUE HELP.'
      READ(5,'(A1)',ERR=1,END=900) KID
      CALL toupper(KID)
      IF(KID.EQ.'E') IEND=1
  900 RETURN
      END SUBROUTINE TRHLPQ
