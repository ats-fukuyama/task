!!! Miscellaneous libraries related to the physics or physical model

!***************************************************************
!
!   Coefficient function of CDBM model
!
!***************************************************************

pure REAL(8) FUNCTION TRCOFS(S,ALFA,RKCV)

  implicit none
  real(8), intent(in) :: S, ALFA, RKCV
  real(8) :: SA, FS1, FS2

  IF(ALFA > 0.D0) THEN
     SA = S - ALFA
     IF(SA > 0.D0) THEN
        FS1 = (1.D0 + 9.0D0 * SQRT(2.D0) * SA**2.5D0) &
             &  / (SQRT(2.D0) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2 + 2.0D0 * SA**3))
     ELSE
        FS1 = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SA) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2))
     ENDIF
     IF(RKCV > 0.D0) THEN
        FS2 = SQRT(RKCV)**3 / S**2
     ELSE
        FS2 = 0.D0
     ENDIF
  ELSE
     SA = ALFA - S
     IF(SA > 0.D0) THEN
        FS1 = (1.D0 + 9.0D0 * SQRT(2.D0) * SA**2.5D0) &
             &  / (SQRT(2.D0) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2 + 2.0D0 * SA**3))
     ELSE
        FS1 = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SA) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2))
     ENDIF
     IF(RKCV < 0.D0) THEN
        FS2 = SQRT(-RKCV)**3 / S**2
     ELSE
        FS2 = 0.D0
     ENDIF
  ENDIF
  TRCOFS = MAX(FS1,FS2)

END FUNCTION TRCOFS

!***************************************************************
!
!   Correction factor for resistivity
!     (Hirshman and Sigmar, (1981), Eq. (7.36))
!
!***************************************************************

pure REAL(8) FUNCTION CORR(X)
  ! X is the effective charge number
  real(8), intent(in) :: X

  CORR = (1.D0 + (1.198D0 + 0.222D0 * X) * X) * X &
  &    / (1.D0 + (2.966D0 + 0.753D0 * X) * X)

END FUNCTION CORR
