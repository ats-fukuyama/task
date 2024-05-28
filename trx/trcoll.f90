!     ***********************************************************

!           COULOMB LOGARITHM

!     ***********************************************************

      FUNCTION coulomb_log(NS1,NS2,ANEL,TL)

!     ANEL : electron density [10^20 /m^3]
!     TL   : electron or ion temperature [keV]
!            in case of ion-ion collision, TL becomes ion temp.

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      INTEGER:: NS1,NS2
      REAL(rkind)   :: ANEL,TL,coulomb_log

      ! Coulomb log: Tokamaks 2Ed. p.661
      
      IF(NS1.EQ.1.AND.NS2.EQ.1) THEN
         coulomb_log=14.9D0-0.5D0*LOG(ANEL)+LOG(TL)
      ELSE
         IF(NS1.EQ.1.OR.NS2.EQ.1) THEN
            coulomb_log=15.2D0-0.5D0*LOG(ANEL)+LOG(TL)
         ELSE
            coulomb_log=17.3D0-0.5D0*LOG(ANEL)+1.5D0*LOG(TL)
         ENDIF
      ENDIF

      RETURN
      END FUNCTION coulomb_log

!     ***********************************************************

!           COLLISION TIME

!     ***********************************************************

!     between electrons and ions

      FUNCTION FTAUE(ANEL,ANIL,TEL,ZL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TEL  : electron temperature [kev]
!     ZL   : ion charge number

      USE TRCOMM, ONLY : AEE, AME, EPS0, PI, PZ, RKEV, rkind
      IMPLICIT NONE
      REAL(rkind) :: ANEL, ANIL, TEL, ZL, FTAUE
      REAL(rkind) :: COEF, coulomb_log

      IF(ABS(ANIL).LE.1.D-8) THEN
         FTAUE=1.D8
         RETURN
      END IF
      COEF = 6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)/(AEE**4*1.D20)
      IF(ZL-PZ(2).LE.1.D-7) THEN
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANIL*ZL**2*coulomb_log(1,2,ANEL,TEL))
      ELSE
!     If the plasma contains impurities, we need to consider the
!     effective charge number instead of ion charge number.
!     From the definition of Zeff=sum(n_iZ_i^2)/n_e,
!     n_iZ_i^2 is replaced by n_eZ_eff at the denominator of tau_e.
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANEL*ZL*coulomb_log(1,2,ANEL,TEL))
      ENDIF

      RETURN
      END FUNCTION FTAUE

!     between ions and ions

      FUNCTION FTAUI(ANEL,ANIL,TIL,ZL,PAL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TIL  : ion temperature [kev]
!     ZL   : ion charge number
!     PAL  : ion atomic number

      USE TRCOMM, ONLY : AEE, AMP, EPS0, PI, RKEV, rkind
      IMPLICIT NONE
      REAL(rkind):: ANEL, ANIL, PAL, TIL, ZL, FTAUI
      REAL(rkind):: COEF, coulomb_log

      IF(ABS(ANIL).LE.1.D-8) THEN
         FTAUI=1.D8
         RETURN
      END IF
      COEF = 12.D0*PI*SQRT(PI)*EPS0**2*SQRT(PAL*AMP)/(AEE**4*1.D20)
      FTAUI = COEF*(TIL*RKEV)**1.5D0/(ANIL*ZL**4*coulomb_log(2,2,ANEL,TIL))

      RETURN
      END FUNCTION FTAUI
