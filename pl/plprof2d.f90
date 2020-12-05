MODULE plprof2d

PRIVATE

PUBLIC plsmag11,plsmag13,plsden11,plsden13,plspsi

CONTAINS

!     ****** set psi ******

SUBROUTINE PLSPSI(RL,ZL,PSI)

  use plcomm,ONLY: RA,rkind
  implicit none
  real(rkind),intent(in) :: RL,ZL
  real(rkind),intent(out):: PSI

  PSI=(RL*RL+ZL*ZL)/(RA*RA)

  RETURN
END SUBROUTINE PLSPSI

SUBROUTINE PLSMAG11(R,Z,BABS,AL)

  use plcomm
  USE plload
  implicit none
  integer :: I
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: BABS,AL(3)
  real(rkind) :: rfactor,zfactor,br,bz,bt

  rfactor=(r-r_corner(1))/(r_corner(2)-r_corner(1))
  zfactor=(z-z_corner(1))/(z_corner(3)-z_corner(1))

  br=br_corner(1)+(br_corner(2)-br_corner(1))*rfactor &
                 +(br_corner(3)-br_corner(1))*zfactor
  bt=bt_corner(1)+(bt_corner(2)-bt_corner(1))*rfactor &
                 +(bt_corner(3)-bt_corner(1))*zfactor
  bz=bz_corner(1)+(bz_corner(2)-bz_corner(1))*rfactor &
                 +(bz_corner(3)-bz_corner(1))*zfactor
  babs=SQRT(br*br+bt*bt+bz*bz)
  al(1)=br/babs
  al(2)=bt/babs
  al(3)=bz/babs
  RETURN
END SUBROUTINE PLSMAG11

SUBROUTINE PLSMAG13(R,Z,BABS,AL)

  use PLcomm
  implicit none
  integer :: I
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: BABS,AL(3)
  real(rkind) :: BLO(3),LR,LZ
  real(rkind) :: L,Q

! L : distance from the center of plasma
! Q : safety factor

! --- initialize ---
  LR = R-RR
  LZ = Z
  L  = sqrt(LR*LR+LZ*LZ)

  Q0 = 1.d0
  QA = 3.d0
  Q  = Q0+(QA-Q0)*(L/RA)**2

  BLO(1) =-BB*LZ/(Q*RR)  ! r   direction
  BLO(2) = BB*RR/(RR+LR) ! phi direction
  BLO(3) = BB*LR/(Q*RR)  ! z   direction
     
  BABS=0.d0
  DO I=1,3
     BABS=BABS+BLO(I)*BLO(I)
  ENDDO
  BABS=SQRT(BABS)

  DO I=1,3
     if(BABS.eq.0.d0) then
        AL(I)=0.0
     else
        AL(I)=BLO(I)/BABS
     end if
  ENDDO
  
  RETURN
END SUBROUTINE PLSMAG13

!     ****** set density & collision frequency ******

SUBROUTINE PLSDEN11(R,Z,RN,RTPR,RTPP)

  USE plcomm
  USE plcomm_type
  implicit none
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM)
  real(rkind) :: rfactor,zfactor
  INTEGER :: ns

  ! --- set FACT ---

  rfactor=(r-r_corner(1))/(r_corner(2)-r_corner(1))
  zfactor=(z-z_corner(1))/(z_corner(3)-z_corner(1))

  ! --- set DENSITY

  SELECT CASE(MODELN)
  CASE(0)
     DO ns=1,nsmax
        rn(ns)=pn_corner(1,ns) &
              +(pn_corner(2,ns)-pn_corner(1,ns))*rfactor &
              +(pn_corner(3,ns)-pn_corner(1,ns))*zfactor
        rtpr(ns)=ptpr_corner(1,ns) &
                +(ptpr_corner(2,ns)-ptpr_corner(1,ns))*rfactor &
                +(ptpr_corner(3,ns)-ptpr_corner(1,ns))*zfactor
        rtpp(ns)=ptpp_corner(1,ns) &
                +(ptpp_corner(2,ns)-ptpp_corner(1,ns))*rfactor &
                +(ptpp_corner(3,ns)-ptpp_corner(1,ns))*zfactor
     END DO
  CASE(1)
     DO ns=1,nsmax
        rn(ns)=pn_corner(1,ns) &
              +(pn_corner(2,ns)-pn_corner(1,ns))*rfactor**2 &
              +(pn_corner(3,ns)-pn_corner(1,ns))*zfactor**2
        rtpr(ns)=ptpr_corner(1,ns) &
                +(ptpr_corner(2,ns)-ptpr_corner(1,ns))*rfactor**2 &
                +(ptpr_corner(3,ns)-ptpr_corner(1,ns))*zfactor**2
        rtpp(ns)=ptpp_corner(1,ns) &
                +(ptpp_corner(2,ns)-ptpp_corner(1,ns))*rfactor**2 &
                +(ptpp_corner(3,ns)-ptpp_corner(1,ns))*zfactor**2
     END DO
  END SELECT
  RETURN
END SUBROUTINE PLSDEN11

SUBROUTINE PLSDEN13(R,Z,RN,RTPR,RTPP)

  USE plcomm
  USE plcomm_type
  implicit none
  integer :: NS,NSI
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM)
  real(rkind) :: LR,LZ,FACT,PSI

  ! --- set FACT ---

  SELECT CASE(MODEL_NPROF)
  CASE(0)
     FACT=1.D0
  CASE(1)
     LR=R-RR
     LZ=Z
     call PLSPSI(LR,LZ,PSI)
     if(PSI.lt.1.D0) then
        FACT=1.D0
     else
        FACT=0.D0
     endif
  CASE(2)
     LR=R-RR
     LZ=Z
     call PLSPSI(LR,LZ,PSI)
     if(PSI.lt.1.D0) then
        FACT=1.D0-PSI
     else
        FACT=0.D0
     endif
  END SELECT

  ! --- set density at NODE NN ---

  DO NS=1,NSMAX
     RN(NS)  =(PN(NS)-PNS(NS))*FACT+PNS(NS)
     RTPR(NS)=PTPR(NS)
     RTPP(NS)=PTPP(NS)
  ENDDO
  RETURN
END SUBROUTINE PLSDEN13

END MODULE plprof2d
