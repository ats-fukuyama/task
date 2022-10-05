!     $Id: wfprof.f90,v 1.6 2011/12/15 19:06:55 maruyama Exp $

!     ****** set psi ******

SUBROUTINE WFBPSI(R,Z,PSI)

  USE wfcomm,ONLY: rkind,PI,ZBB,RMIR,BB,MODELG,MODELB,NCOILMAX, &
       RCOIL,ZCOIL,BCOIL,Hpitch1,HA1,rkind
  USE libbes,ONLY: BESINX
  IMPLICIT NONE
  REAL(rkind),INTENT(IN) :: R,Z
  REAL(rkind),INTENT(OUT):: PSI
  INTEGER:: NCOIL
  REAL(rkind):: A0,A1,RL,ZL,WFPSIC,XH1,YH1,PSIG

  SELECT CASE(MODELG)
  CASE(0)  ! slab
     SELECT CASE(MODELB)
     CASE(0)
        A0=0.5D0*(1.D0+RMIR)*BB
        A1=0.5D0*(1.D0-RMIR)*BB
        RL = PI* R/ZBB
        ZL = PI* Z/ZBB
        PSI = A0*R +A1*(ZBB/PI)*COS(ZL)*SINH(RL )
     CASE(1)
        PSI=0.D0
        DO NCOIL=1,NCOILMAX
           IF(R.NE.0.D0) THEN
              PSI=PSI+BCOIL(NCOIL) &
                      *LOG(RCOIL(NCOIL)) &
                      /LOG(SQRT((R-RCOIL(NCOIL))**2+(Z-ZCOIL(NCOIL))**2))
           ENDIF
        ENDDO
     END SELECT
  CASE(1)  ! cylindrical
     SELECT CASE(MODELB)
     CASE(0)
        A0=0.5D0*(1.D0+RMIR)*BB
        A1=0.5D0*(1.D0-RMIR)*BB
        RL = PI* R/ZBB
        ZL = PI* Z/ZBB
        PSI = 0.5D0*A0*R *R +A1*(ZBB/PI)**2*COS(ZL)*RL *BESINX(1,RL)
     CASE(1)
         PSI=0.D0
         DO NCOIL=1,NCOILMAX
            IF(R.NE.0.D0) THEN
               PSI=PSI+BCOIL(NCOIL)*WFPSIC(R,Z-ZCOIL(NCOIL),RCOIL(NCOIL))
            ENDIF
         ENDDO
     END SELECT
  CASE(3,5,6,8,9)   ! toroidal
     PSI=PSIG(R,Z)
  CASE(11)   ! straight helical
     XH1=Hpitch1*R
     YH1=Hpitch1*Z
     PSI=(XH1**2+YH1**2+2.D0*HA1*(XH1**2-YH1**2))
  END SELECT
  RETURN
END SUBROUTINE WFBPSI

!   --- psi function for circular coil ---

  FUNCTION WFPSIC(RL,ZL,RC)

!        R*A_psi

    USE wfcomm,ONLY: rkind,PI
    USE libell,ONLY: ELLFC,ELLEC
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RL,ZL,RC
    REAL(rkind):: WFPSIC
    REAL(rkind):: RX,RK
    INTEGER:: IERR1,IERR2

    RX=SQRT(RC**2+RL**2+ZL**2+2.D0*RC*RL)
    RK=SQRT(4.D0*RC*RL)/RX
    WFPSIC=(RC/PI)*RX &
          *((1.D0-0.5D0*RK**2)*ELLFC(RK,IERR1)-ELLEC(RK,IERR2))
    RETURN
  END FUNCTION WFPSIC


!     ****** set magnetic field ******

SUBROUTINE WFSMAG(R,Z,BABS,AL)

  USE wfcomm,ONLY: MODELG,rkind
  USE plload
  implicit none
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: BABS,AL(3)
  REAL(rkind):: BR,BZ,BT
  integer:: IERR

  SELECT CASE(MODELG)
  CASE(0)
     CALL WFSMAG0(R,Z,BABS,AL)
  CASE(1,2)
     CALL WFSMAG2(R,Z,BABS,AL)
  CASE(12)
     CALL pl_read_p2Dmag(R,Z,BR,BZ,BT,IERR)
     BABS=SQRT(BR**2+BZ**2+BT**2)
     AL(1)=BR/BABS
     AL(2)=BT/BABS
     AL(3)=BZ/BABS
  END SELECT
  RETURN
END SUBROUTINE WFSMAG

SUBROUTINE WFSMAG0(R,Z,BABS,AL)

  use wfcomm
  USE plload
  implicit none
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
END SUBROUTINE WFSMAG0

SUBROUTINE WFSMAG2(R,Z,BABS,AL)

  use wfcomm
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

! --- set B field at NODE NN ---
  if (MODELB.eq.0) then
     BLO(1) =-BB*LZ/(Q*RR)  ! r   direction
     BLO(2) = BB*RR/(RR+LR) ! phi direction
     BLO(3) = BB*LR/(Q*RR)  ! z   direction
  else
     BLO=0.d0
  end if
     
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
END SUBROUTINE WFSMAG2

!     ****** set density & collision frequency ******

SUBROUTINE WFSDEN(R,Z,RN,RTPR,RTPP,RZCL)

  use wfcomm,ONLY: modelg,nsm,nsmax,mdamp,rkind, &
       model_coll_enhance,factor_coll_enhance, &
       xpos_coll_enhance,xwidth_coll_enhance, &
       ypos_coll_enhance,ywidth_coll_enhance
       
  use plload
  implicit none
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  REAL(rkind):: RU(NSM),factor
  INTEGER:: NSMAXL,IERR,NS

  SELECT CASE(MODELG)
  CASE(0)
     CALL WFSDEN0(R,Z,RN,RTPR,RTPP,RZCL)
  CASE(1,2)
     CALL WFSDEN2(R,Z,RN,RTPR,RTPP,RZCL)
  CASE(12)
     CALL pl_read_p2D(R,Z,RN,RTPR,RTPP,RU,IERR)
     IF(mdamp.NE.0) THEN
        RN(NSMAX)=0.D0
        RTPR(NSMAX)=1.D0
        RTPP(NSMAX)=1.D0
        RU(NSMAX)=0.D0
     END IF
     CALL WFCOLL(rn,rtpr,rtpp,rzcl,0)
     SELECT CASE(model_coll_enhance)
     CASE(1)
        factor=factor_coll_enhance &
             *(1.D0+EXP(-(R-xpos_coll_enhance)**2/xwidth_coll_enhance**2))
     CASE(2)
        factor=factor_coll_enhance &
             *(1.D0+EXP(-(Z-ypos_coll_enhance)**2/ywidth_coll_enhance**2))
     CASE DEFAULT
        factor=1.D0
     END SELECT
     IF(mdamp.EQ.0) THEN
        DO NS=1,NSMAX
           RZCL(NS)=RZCL(NS)*factor
        END DO
     ELSE
        DO NS=1,NSMAX-1
           RZCL(NS)=RZCL(NS)*factor
        END DO
     END IF
  END SELECT
  RETURN
END SUBROUTINE WFSDEN

SUBROUTINE WFSDEN0(R,Z,RN,RTPR,RTPP,RZCL)

  use wfcomm
  implicit none
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
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
  CALL WFCOLL(rn,rtpr,rtpp,rzcl,0)

  RETURN
END SUBROUTINE WFSDEN0

SUBROUTINE WFSDEN2(R,Z,RN,RTPR,RTPP,RZCL)

  use wfcomm
  use plprof2d
  implicit none
  integer :: NS
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  real(rkind) :: LR,LZ
  real(rkind) :: FACT,PSI

  ! --- set FACT ---

  if(MODELP.eq.0) then
     FACT=1.D0
  else
     LR=R-RR
     LZ=Z
     call PLSPSI(LR,LZ,PSI)
     if(PSI.lt.1.D0) then
        if(MODELP.eq.1) then
           FACT=1.D0
        elseif(MODELP.eq.2) then
           FACT=1.D0-PSI
        else
           write(6,*) 'XX WDSDEN: UNKNOWN MODELP = ',MODELP
        endif
     else
        FACT=0.D0
     endif
  end if

  ! --- set density at NODE NN ---

  do NS=1,NSMAX
     RN(NS)  =(PN(NS)-PNS(NS))*FACT+PNS(NS)
     RTPR(NS)=PTPR(NS)
     RTPP(NS)=PTPP(NS)
  enddo

  ! --- set collision frequency ---
  
  CALL WFCOLL(rn,rtpr,rtpp,rzcl,0)

!  WRITE(6,*) 'ZND= ',ZND(IN)
!  WRITE(6,*) 'RN = ',RN(1),RN(2)
!  WRITE(6,*) 'RT = ',RTPR(1),RTPR(2)
!  WRITE(6,*) 'RZ = ',RZCL(1),RZCL(2)
!  WRITE(6,*) 'E  = ',RNUE,RNUEE,RNUEI,RNUEN
!  WRITE(6,*) 'I  = ',RNUI,RNUIE,RNUII,RNUIN
!  STOP

  RETURN
END SUBROUTINE WFSDEN2

SUBROUTINE WFCOLL(rn,rtpr,rtpp,rzcl,id)
  USE wfcomm
  IMPLICIT NONE
  REAL(rkind),INTENT(IN):: rn(NSM),rtpr(NSM),rtpp(NSM)
  REAL(rkind),INTENT(OUT):: rzcl(NSM)
  INTEGER,INTENT(IN):: id  ! id=0 without output, id=1 with output
  REAL(rkind):: TE,TI,RNTI,RNZI,RLAMEE,RLAMEI,RLAMII,SN,PNN0
  REAL(rkind):: VTE,RNUEE,RNUEI,RNUEN,RNUE
  REAL(rkind):: VTI,RNUIE,RNUII,RNUIN,RNUI
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
              RZCL(NS)=RNUE/(2.D6*PI*RF)
              IF(ID.NE.0) THEN
                 WRITE(6,'(A,1P3E12.4)') &
                      'PPN0,PTN0,PNN0    =',PPN0,PTN0,PNN0
                 WRITE(6,'(A,I12,1P4E12.4)') &
                      'NS,RN,PTPR,PTPP,TE=',NS,RN(1),PTPR(1),PTPP(1),TE
                 WRITE(6,'(A,1P5E12.4)') &
                      'RNUEE/I/N/TOT/RZCL=', &
                      RNUEE,RNUEI,RNUEN,RNUE,RZCL(1)
              END IF
           ELSE
              TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
              VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
              RNUIE=PZ(NS)**2*RN(1)*RLAMEI &
                       /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
              RNUII=RNZI*RLAMII &
                       /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
              RNUIN=PNN0*SN*0.88D0*VTI
              RNUI=RNUIE+RNUII+RNUIN
              RZCL(NS)=RNUI/(2.D6*PI*RF)
              IF(ID.NE.0) THEN
                 WRITE(6,'(A,I12,1P4E12.4)') &
                      'NS,RN,RTPR,RTPP,TI=', &
                      NS,RN(NS),PTPR(NS),PTPP(NS),TI
                 WRITE(6,'(A,1P5E12.4)') &
                      'RNUIE/I/N/TOT/RZCL=', &
                      RNUIE,RNUII,RNUIN,RNUI,RZCL(NS)
              END IF
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
END SUBROUTINE WFCOLL
