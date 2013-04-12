!     $Id$

!     ****** PSI ******

SUBROUTINE WFSPSI(X,Y,Z,PSI)

  use wfcomm
  implicit none
  real(8) :: PSI,X,Y,Z
  real(8) :: DFZ,DFZMIN,XPT,YPT,ZPT
  integer:: J

  IF(MODELB.EQ.0.OR.&
 &   MODELB.EQ.1.OR.&
 &   MODELB.EQ.2) THEN
     PSI=(X*X+Y*Y)/(RA*RA)
  ELSE
!
! ----- Add. By YOKOYAMA Mar./05/2013 ----
!
!      XCM0 = 0.D0
!      YCM0 = 0.D0
!      ZCM0 = 0.D0
!         CALL MAG(XCM0,YCM0,ZCM0,BXG0,BYG0,BZG0)
!      BX0 = BXG0/1.D4
!      BY0 = BYG0/1.D4
!      BZ0 = BZG0/1.D4
!      BA0 = SQRT(BX0*BX0+BY0*BY0+BZ0*BZ0)      
!      PSIA = RA*RA*BZ0
!
!      XCM = 0.D0
!      YCM = 0.D0
!      ZCM = Z*1.D2
!         CALL MAG(XCM,YCM,ZCM,BXG,BYG,BZG)
!      BX = BXG/1.D4
!      BY = BYG/1.D4
!      BZ = BZG/1.D4
!      BA = SQRT(BX*BX+BY*BY+BZ*BZ)
!      PSI = (X*X+Y*Y)*BZ
!
!      PSI = PSI/PSIA
!
!     セントラル部中央に存在するリミタ(半径18cm)で切られる磁力管を境界として，
!     密度分布を作成する．磁力管は楕円で近似する．
!     ZPT: Z座標
!     XPT: 磁力管のX軸との交点(X座標)
!     YPT: 磁力管のY軸との交点(Y座標)
!
!     与えられたZ座標に，最も近いZ座標を持つ磁力管断面を探す．
!     FLZ,FLX,FLYは cm 単位
!     ZPT,XPT,YPTは  m 単位

      J=1
      ZPT = FLZ(J)/1.D2
      DFZMIN = ABS(Z-ZPT)
      XPT = FLX(J)/1.D2
      YPT = FLY(J)/1.D2
      DO J=2,NGFLIN
         ZPT = FLZ(J)/1.D2
         DFZ = ABS(Z-ZPT)
         IF(DFZ.LT.DFZMIN) THEN
            DFZMIN = DFZ
            XPT = FLX(J)/1.D2
            YPT = FLY(J)/1.D2
         ENDIF
      ENDDO
!
      PSI = X*X/(XPT*XPT) + Y*Y/(YPT*YPT)
      IF(PSI.LT.1.D0) THEN
         FACTC = 1.D0 - PSI
!         FACTA = RA*RA/(XPT*YPT)*(1.D0-PSI)
         FACTA = 1.D0 - PSI
      ELSE
         FACTC = 0.D0
         FACTA = 0.D0
      ENDIF
   ENDIF
!
!
!
  RETURN
END SUBROUTINE WFSPSI

!     ****** MAGNETIC FIELD PROFILE ******

SUBROUTINE WFSMAG(IN,BABS,AL)

  use wfcomm
  implicit none
  integer :: IN,I
  real(8) :: BLO(3),AL(3),XC,YC,ZC,BABS
  real(8) :: XT,YT,ZT,XCM,YCM,ZCM,BX,BY,BZ
  real(8) :: R,Q,PHI,AC(3),BLOT(3)

  XC=XND(IN)
  YC=YND(IN)
  ZC=ZND(IN)
  call RCTORT(XC,YC,ZC,XT,YT,ZT)

  R=sqrt(XT*XT+YT*YT)
  Q0=1.d0
  QA=3.d0
  Q = Q0+(QA-Q0)*(R/RA)**2

  IF(MODELB.EQ.0) THEN
     BLO(1)=BB
     BLO(2)=0.D0
     BLO(3)=0.D0
  ELSEIF(MODELB.EQ.1) THEN
     BLO(1)=0.D0
     BLO(2)=BB
     BLO(3)=0.D0
  ELSEIF(MODELB.EQ.2) THEN
     BLO(1)=0.D0
     BLO(2)=0.D0
     BLO(3)=BB
  ELSEIF(MODELB.eq.3) THEN
     if(MODELG.eq.1) then
        BLOT(1)=-BB*YT/(Q*RR)
        BLOT(2)= BB*XT/(Q*RR)
        BLOT(3)= BB*RR/(RR+XT)
        CALL ATTOAC(BLOT(1),BLOT(2),BLOT(3),ZT,BLO(1),BLO(2),BLO(3))
     else
        BLO(1)=-BB*YT/(Q*RR)
        BLO(2)= BB*XT/(Q*RR)
        BLO(3)= BB*RR/(RR+XT)
     end if
!
! ----- Add. By YOKOYAMA Mar./05/2013 ----  
  ELSEIF(MODELB.eq.4) THEN
!        "m" --> "cm"
         XCM=XC*1.0D2
         YCM=YC*1.0D2
         ZCM=ZC*1.0D2
!
!        Magnetic Field of GAMMA 10
         CALL MAG(XCM,YCM,ZCM,BX,BY,BZ)
!        "Gauss" --> "Tesla"
         BLO(1)=BX/1.0D4
         BLO(2)=BY/1.0D4
         BLO(3)=BZ/1.0D4
  ENDIF
!
! ----- Mar./05/2013 -----
! 
  BABS=0.D0
  DO I=1,3
     BABS=BABS+BLO(I)*BLO(I)
  ENDDO
  BABS=SQRT(BABS)

  DO I=1,3
     if(BABS.eq.0.0) then
        AL(I)=0.0
     else
        AL(I)=BLO(I)/BABS
     end if
  ENDDO
  
  RETURN
END SUBROUTINE WFSMAG

!     ****** DENSITY & TEMPERATURE PROFILE ******

SUBROUTINE WFSDEN(IN,RN,RTPR,RTPP,RZCL)

  use wfcomm
  implicit none
  integer :: IN,NS,NSI,J
  real(8) :: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM),FACT,PSI,Z
  real(8) :: TE,TI,RLAMEE,RLAMEI,RLAMII,SN,PNN0,VTE,RNUEE,RNUEI,RNUEN
  real(8) :: RNUE,VTI,RNUIE,RNUII,RNUIN,RNUI
  real(8) :: XT,YT,ZT
  real(8) :: DFZ,DFZMIN,XPT,YPT,ZPT,WAZ

  IF(MODELP.EQ.0) THEN
     FACT=1.D0
  ELSE
     if(MODELG.eq.1) then
        CALL RCTORT(XND(IN),YND(IN),ZND(IN),XT,YT,ZT)
        CALL WFSPSI(XT,YT,ZT,PSI)
     else
        CALL WFSPSI(XND(IN),YND(IN),ZND(IN),PSI)
     end if
     IF(PSI.LT.1.D0) THEN
        IF(MODELP.EQ.1) THEN
           FACT=1.D0
        ELSEIF(MODELP.EQ.2) THEN
! ----- Add. By YOKOYAMA Mar./05/2013 ----  
               IF(PSI.LT.1.D0) THEN
                  FACT=1.D0-PSI
               ELSE
                  FACT=0.D0
               ENDIF
! ----- Mar./05/2013 -----
        ELSEIF(MODELP.EQ.3) THEN
           FACT=EXP(-(ZPMAX-ZND(IN))/0.10D0)
        ELSE
           WRITE(6,*) 'XX WDSDEN: UNKNOWN MODELP = ',MODELP
        ENDIF
     ELSE
        FACT=0.D0
     ENDIF
! ----- Mod. By YOKOYAMA Mar./05/2013 ----  
!     Z=ZND(IN)
!     IF(Z.LT.ZPMIN.OR.Z.GT.ZPMAX) THEN
!        FACT=0.D0
!     ENDIF
! ----- Mar./05/2013 -----
  ENDIF
  
  DO NS=1,NSMAX
!
! ----- Mod. & Add. By YOKOYAMA Mar./05/2013 ----  
!     RN(NS)  =(PN(NS)  -PNS(NS))*FACT+PNS(NS)
!     与えられたZ座標に，最も近いZ座標を持つ磁力管断面を探す．
!     FLZ,FLX,FLYは cm 単位
!     ZPT,XPT,YPTは  m 単位
      DFZMIN = 1.D2
      DO J=1,NGFLIN
         ZPT = FLZ(J)/1.D2
         DFZ = ABS(ZND(IN)-ZPT)
         IF(DFZ.LT.DFZMIN) THEN
            DFZMIN = DFZ
            XPT = FLX(J)/1.D2
            YPT = FLY(J)/1.D2
         ENDIF
      ENDDO
!
!    Density Profile in G-10 Central Cell
!
         IF(ABS(ZND(IN)).LT.2.8D0) THEN
            RN(NS)&
     &      = (PN(NS)-PNS(NS))*FACTC&
!  ----  Z方向分布
     &        *(0.52D0*EXP(-5.D0*(DABS(ZND(IN))/2.7D0)**7)+0.48D0)&
!  ----  X方向分布
     &        *EXP(-1.D0*(DABS(XND(IN))/XPT)**3)&
!  ----  Y方向分布
     &        *EXP(-1.D0*(DABS(YND(IN))/YPT)**3)&
!  ----  周辺密度
     &        +PNS(NS)
!
!    Density Profile in Transition and Anshor Cell
!
         ELSE
            WAZ = DABS(ZND(IN)-5.2D0)
!            WAB = 1.D0/0.3128231D0
            RN(NS)&
     &      =(PN(NS)-PNS(NS))&
!  ----  Z方向分布
     &       *(&
!     &            FACTA/4.98414D0
!     &            *((WAB-0.48D0)*EXP(-5.D0*(WAZ/0.4D0)**3)+0.48D0)
     &            FACTA*(1.D0-(WAZ/0.8D0)**2)*EXP(-1.D0*(WAZ/1.D0)**2)&
     &         )&
!  ----  X方向分布
     &       *EXP(-1.D0*(DABS(XND(IN))/(XPT*1.D0))**2)&
!  ----  Y方向分布
     &       *EXP(-1.D0*(DABS(YND(IN))/(YPT*1.D0))**2)&
!  ----  周辺密度
     &       + PNS(NS)
         ENDIF
!          IF((XND(IN).EQ.0.D0).and.
!     &       (YND(IN).EQ.0.D0)) THEN
!                WRITE(6,'(4F20.7)') ZND(IN),XND(IN),YND(IN),RN(NS)
!          ENDIF
!  ----
     RTPR(NS)=PTPR(NS)
     RTPP(NS)=PTPP(NS)
  ENDDO
  
  IF(RN(1).GT.0.D0) THEN
     
     TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
     TI=(RTPR(2)+2.D0*RTPP(2))*1.D3/3.D0
     RLAMEE= 8.0D0+2.3D0*(LOG10(TE)-0.5D0*LOG10(RN(1)))
     RLAMEI= RLAMEE+0.3D0
     RLAMII=12.1D0+2.3D0*(LOG10(TI)-0.5D0*LOG10(RN(1)))
     SN=1.D-20
     PNN0=PPN0/(PTN0*AEE)
     
     DO NS=1,NSMAX
        IF(PZCL(NS).EQ.0) THEN
           IF(NS.EQ.1) THEN
              TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
              VTE=SQRT(2.D0*TE*AEE/AME)
              RNUEE=RN(1)*RLAMEE &
                   &                 /(1.24D-4*SQRT(TE*1.D-3)**3)
              RNUEI=0.D0
              DO NSI=2,NSMAX
                 RNUEI=RNUEI+PZ(NSI)**2*RN(NSI)
              ENDDO
              RNUEI=RNUEI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
              RNUEN=PNN0*SN*0.88D0*VTE
              RNUE=RNUEE+RNUEI+RNUEN
              RZCL(NS)=RNUE/(2.D6*PI*RF)
           ELSE
              TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
              VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
              RNUIE=PZ(NS)**2*RN(1)*RLAMEI &
                   &                 /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
              RNUII=0.D0
              DO NSI=2,NSMAX
                 RNUII=RNUII+PZ(NSI)**2*RN(NSI)
              ENDDO
              RNUII=RNUII*RLAMII &
                   &                 /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
              RNUIN=PNN0*SN*0.88D0*VTI
              RNUI=RNUIE+RNUII+RNUIN
              RZCL(NS)=RNUI/(2.D6*PI*RF)
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
  
!  WRITE(6,*) 'ZND= ',ZND(IN)
!  WRITE(6,*) 'RN = ',RN(1),RN(2)
!  WRITE(6,*) 'RT = ',RTPR(1),RTPR(2)
!  WRITE(6,*) 'RZ = ',RZCL(1),RZCL(2)
!  WRITE(6,*) 'E  = ',RNUE,RNUEE,RNUEI,RNUEN
!  WRITE(6,*) 'I  = ',RNUI,RNUIE,RNUII,RNUIN
!  STOP

  RETURN
END SUBROUTINE WFSDEN
