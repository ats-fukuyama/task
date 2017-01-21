MODULE plprof2d

PRIVATE

PUBLIC plsmag11,plsmag13,plsden11,plsden13 

CONTAINS

!     ****** set psi ******

SUBROUTINE PLSPSI(RL,ZL,PSI)

  use plcomm,ONLY: RA
  implicit none
  real(8),intent(in) :: RL,ZL
  real(8),intent(out):: PSI

  PSI=(RL*RL+ZL*ZL)/(RA*RA)

  RETURN
END SUBROUTINE PLSPSI

SUBROUTINE PLCOLL(rn,rtpr,rtpp,rzcl)
  USE plcomm
  IMPLICIT NONE
  REAL(8),DIMENSION(NSM):: rn,rtpr,rtpp,rzcl
  REAL(8):: TE,TI,RNTI,RNZI,RLAMEE,RLAMEI,RLAMII,SN,PNN0
  REAL(8):: VTE,RNUEE,RNUEI,RNUEN,RNUE
  REAL(8):: VTI,RNUIE,RNUII,RNUIN,RNUI
  INTEGER:: ns

  ! --- set collision frequency ---
  
  IF(rn(1).GT.0.D0) THEN
     TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
     RNTI=0.D0
     RNZI=0.D0
     DO ns=2,nsmax
        RNTI=RNTI+RN(ns)*(RTPR(ns)+2.D0*RTPP(ns))*1.D3/3.D0
        RNZI=RNZI+RN(ns)*PZ(ns)**2
     END DO
     TI=RNTI/rn(1)
     RLAMEE= 8.0D0+2.3D0*(LOG10(TE)-0.5D0*LOG10(RN(1)))
     RLAMEI= RLAMEE+0.3D0
     RLAMII=12.1D0+2.3D0*(LOG10(TI)-0.5D0*LOG10(RN(1)))

     SN=1.D-20 ! tytpical ionizatioin crosssection
     PNN0=PPN0/(PTN0*AEE) ! neutral density

     DO NS=1,NSMAX
        IF(PZCL(NS).EQ.0) THEN
           IF(NS.EQ.1) THEN
              VTE=SQRT(2.D0*TE*AEE/AME)
              RNUEE=RN(1)*RLAMEE/(1.24D-4*SQRT(TE*1.D-3)**3)
              RNUEI=RNZI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
              RNUEN=PNN0*SN*0.88D0*VTE
              RNUE=RNUEE+RNUEI+RNUEN
              RZCL(NS)=RNUE/(2.D6*PI*RF_PL)
           ELSE
              TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
              VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
              RNUIE=PZ(NS)**2*RN(1)*RLAMEI &
                       /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
              RNUII=RNZI*RLAMII &
                       /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
              RNUIN=PNN0*SN*0.88D0*VTI
              RNUI=RNUIE+RNUII+RNUIN
!              RZCL(NS)=RNUI/(2.D6*PI*RF_PL)
              RZCL(NS)=0.D0
           ENDIF
        ELSE
           RZCL(NS)=PZCL(NS)
        ENDIF
     ENDDO
  ELSE
     DO NS=1,NSMAX
        RZCL(NS)=0.D0
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE PLCOLL


SUBROUTINE PLSMAG11(R,Z,BABS,AL)

  use PLcomm
  USE plload
  implicit none
  integer :: I
  real(8),intent(in) :: R,Z
  real(8),intent(out):: BABS,AL(3)
  real(8) :: rfactor,zfactor,br,bz,bt

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
  real(8),intent(in) :: R,Z
  real(8),intent(out):: BABS,AL(3)
  real(8) :: BLO(3),LR,LZ
  real(8) :: L,Q

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

SUBROUTINE PLSDEN11(R,Z,RN,RTPR,RTPP,RZCL)

  use plcomm
  implicit none
  real(8),intent(in) :: R,Z
  real(8),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  real(8) :: rfactor,zfactor
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
!     write(6,'(1P3E12.4)') rfactor,zfactor
!  DO NS=1,NSMAX
!     write(6,'(1P3E12.4)') pn_corner(1,NS),pn_corner(2,NS),pn_corner(3,NS)
!     write(6,'(1P3E12.4)') ptpr_corner(1,NS),ptpr_corner(2,NS),&
!                           ptpr_corner(3,NS)
!     write(6,'(1P3E12.4)') ptpp_corner(1,NS),ptpp_corner(2,NS),&
!                           ptpp_corner(3,NS)
!     write(6,'(1P3E12.4)') R,Z,rn(ns)
!     write(6,'(1P3E12.4)') rtpr(NS),rtpp(NS),rzcl(NS)
!  END DO
!     STOP

  CALL PLCOLL(rn,rtpr,rtpp,rzcl)

  RETURN
END SUBROUTINE PLSDEN11

SUBROUTINE PLSDEN13(R,Z,RN,RTPR,RTPP,RZCL)

  use plcomm
  implicit none
  integer :: NS,NSI
  real(8),intent(in) :: R,Z
  real(8),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  real(8) :: LR,LZ
  real(8) :: TE,TI,RLAMEE,RLAMEI,RLAMII,SN,PNN0,VTE,RNUEE,RNUEI,RNUEN
  real(8) :: RNUE,VTI,RNUIE,RNUII,RNUIN,RNUI,FACT,PSI

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

  do NS=1,NSMAX
     RN(NS)  =(PN(NS)-PNS(NS))*FACT+PNS(NS)
     RTPR(NS)=PTPR(NS)
     RTPP(NS)=PTPP(NS)
  enddo

  ! --- set collision frequency ---
  
  CALL PLCOLL(rn,rtpr,rtpp,rzcl)

!  WRITE(6,*) 'ZND= ',ZND(IN)
!  WRITE(6,*) 'RN = ',RN(1),RN(2)
!  WRITE(6,*) 'RT = ',RTPR(1),RTPR(2)
!  WRITE(6,*) 'RZ = ',RZCL(1),RZCL(2)
!  WRITE(6,*) 'E  = ',RNUE,RNUEE,RNUEI,RNUEN
!  WRITE(6,*) 'I  = ',RNUI,RNUIE,RNUII,RNUIN
!  STOP

  RETURN
END SUBROUTINE PLSDEN13

END MODULE plprof2d
