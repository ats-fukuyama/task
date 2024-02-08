! wrfdrv.f90

!  --- slave routine for ray tracing ---

!  Y(1)=X
!  Y(2)=Y
!  Y(3)=Z
!  Y(4)=RKX
!  Y(5)=RKY
!  Y(6)=RKZ
!  Y(7)=W

MODULE wrfdrv

  PRIVATE
  PUBLIC wr_fdrv
  PUBLIC wr_fdrvr

CONTAINS

  SUBROUTINE wr_fdrv(x,y,f)

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x,y(NEQ)
    REAL(rkind),INTENT(OUT):: f(7)

    SELECT CASE(model_fdrv)
    CASE(1)
       CALL wr_fdrv1(x,y,f)
    CASE(2)
       CALL wr_fdrv2(x,y,f)
    CASE(3)
       CALL wr_fdrv3(x,y,f)
    END SELECT
    RETURN
  END SUBROUTINE wr_fdrv
       
  ! DRXP=ABS(RXP)*DELDER
  ! DRKXP=ABS(RKXP)*DELDER

  SUBROUTINE wr_fdrv1(X,Y,F)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,DISPXI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y(NEQ)
    REAL(rkind),INTENT(OUT):: F(7)
    REAL(rkind):: VV,TT,RF,XP,YP,ZP,RKXP,RKYP,RKZP,UU,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ,VDU
    REAL(rkind):: DUMMY

      VV=DELDER
      TT=DELDER
      DUMMY=X

      RF=wr_nray_status%RF
      XP=Y(1)
      YP=Y(2)
      ZP=Y(3)
      RKXP=Y(4)
      RKYP=Y(5)
      RKZP=Y(6)
      UU=Y(7)
      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      RXP=MAX(ABS(XP)*VV,TT)
      RYP=MAX(ABS(YP)*VV,TT)
      RZP=MAX(ABS(ZP)*VV,TT)
      RRKXP=MAX(ABS(RKXP)*VV,TT)
      RRKYP=MAX(ABS(RKYP)*VV,TT)
      RRKZP=MAX(ABS(RKZP)*VV,TT)

      DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
      DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
      DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
      DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
      DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
      DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

      SELECT CASE(model_fdrv_ds)
      CASE(0)
         IF(DOMG.GT.0.D0) THEN
            DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
         ELSE
            DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
         ENDIF
      CASE(1)
         DS=DOMG
      END SELECT

      VX   =-DKXP/DS
      VY   =-DKYP/DS
      VZ   =-DKZP/DS

      VKX  = DXP/DS
      VKY  = DYP/DS
      VKZ  = DZP/DS

      VDU  =-2.D0*ABS(DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)/DS)

      F(1)=VX
      F(2)=VY
      F(3)=VZ
      F(4)=VKX
      F(5)=VKY
      F(6)=VKZ
      IF(UU.LT.0.D0) THEN
         F(7)=0.D0
      ELSE
         F(7)=VDU*UU 
      ENDIF

      IF(idebug_wr(12).NE.0) THEN
         WRITE(6,'(A)') '*** idebug_wr(12): wrfdrv'
         WRITE(6,'(A,3ES12.4)') 'x7:',X,Y(7),F(7)
         WRITE(6,'(A,6ES12.4)') 'y :',Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
         WRITE(6,'(A,6ES12.4)') 'f :',F(1),F(2),F(3),F(4),F(5),F(6)
      END IF
      RETURN
  END SUBROUTINE wr_fdrv1

  ! DRXP=dels*DELDER
  ! DRKXP=(RKXP/RXP)*dels*DELDER

  SUBROUTINE wr_fdrv2(X,Y,F)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,DISPXI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y(NEQ)
    REAL(rkind),INTENT(OUT):: F(7)
    REAL(rkind):: VV,TT,RF,XP,YP,ZP,RKXP,RKYP,RKZP,UU,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ,VDU
    REAL(rkind):: DUMMY,AVR,AVKR

    VV=DELDER
    TT=DELDER
    DUMMY=X

    RF=wr_nray_status%RF
    XP=Y(1)
    YP=Y(2)
    ZP=Y(3)
    RKXP=Y(4)
    RKYP=Y(5)
    RKZP=Y(6)
    UU=Y(7)
    
    OMG=2.D6*PI*RF
    ROMG=MAX(ABS(OMG)*VV,TT)
    AVR=dels
    RXP=MAX(AVR*VV,TT)
    RYP=MAX(AVR*VV,TT)
    RZP=MAX(AVR*VV,TT)
    AVKR=dels*SQRT(RKXP**2+RKYP**2+RKZP**2)/SQRT(XP**2+YP**2+ZP**2)
    RRKXP=MAX(AVKR*VV,TT)
    RRKYP=MAX(AVKR*VV,TT)
    RRKZP=MAX(AVKR*VV,TT)

    DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
         -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
    DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
    DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
    DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
    DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
    DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
         -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
    DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
         -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

    SELECT CASE(model_fdrv_ds)
    CASE(0)
       IF(DOMG.GT.0.D0) THEN
          DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
       ELSE
          DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
       ENDIF
    CASE(1)
       DS=DOMG
    END SELECT

    VX   =-DKXP/DS
    VY   =-DKYP/DS
    VZ   =-DKZP/DS

    VKX  = DXP/DS
    VKY  = DYP/DS
    VKZ  = DZP/DS

    VDU  =-2.D0*ABS(DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)/DS)

    F(1)=VX
    F(2)=VY
    F(3)=VZ
    F(4)=VKX
    F(5)=VKY
    F(6)=VKZ
    IF(UU.LT.0.D0) THEN
       F(7)=0.D0
    ELSE
       F(7)=VDU*UU 
    ENDIF

    IF(idebug_wr(12).NE.0) THEN
       WRITE(26,'(A)') '*** idebug_wr(12): wrfdrv'
       WRITE(26,'(A,4ES12.4)') 'x7ds:',X,Y(7),F(7),DS
       WRITE(26,'(A,6ES12.4)') 'y   :',Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
       WRITE(26,'(A,6ES12.4)') 'f   :',F(1),F(2),F(3),F(4),F(5),F(6)
    END IF
    RETURN
  END SUBROUTINE wr_fdrv2

  ! DRXP=dels*DELDER
  ! DRKXP=(RKXP/RXP)*dels*DELDER

  SUBROUTINE wr_fdrv3(X,Y,F)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,DISPXI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y(NEQ)
    REAL(rkind),INTENT(OUT):: F(7)
    REAL(rkind):: VV,TT,RF,XP,YP,ZP,RKXP,RKYP,RKZP,UU,OMG,ROMG
    REAL(rkind):: DROMG,DRXP,DRYP,DRZP,DRKXP,DRKYP,DRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ,VDU
    REAL(rkind):: DUMMY,delk
    REAL(rkind),SAVE:: DXP_save,DYP_save,DZP_save
    REAL(rkind),SAVE:: DKXP_save,DKYP_save,DKZP_save

    VV=DELDER
    TT=DELDER
    DUMMY=X

    RF=wr_nray_status%RF
    XP=Y(1)
    YP=Y(2)
    ZP=Y(3)
    RKXP=Y(4)
    RKYP=Y(5)
    RKZP=Y(6)
    UU=Y(7)

    OMG=2.D6*PI*RF
    DROMG=MAX(ABS(OMG)*VV,TT)

    IF(INITIAL_DRV.EQ.0) THEN
       DRXP=MAX(dels*VV,TT)
       DRYP=MAX(dels*VV,TT)
       DRZP=MAX(dels*VV,TT)
       delk=SQRT(RKXP**2+RKYP**2+RKZP**2)
       DRKXP=MAX(delk*VV,TT)
       DRKYP=MAX(delk*VV,TT)
       DRKZP=MAX(delk*VV,TT)
       INITIAL_DRV=1
    ELSE
       DRXP=MAX(DXP_save*VV,TT)
       DRYP=MAX(DYP_save*VV,TT)
       DRZP=MAX(DZP_save*VV,TT)
       DRKXP=MAX(DKXP_save*VV,TT)
       DRKYP=MAX(DKYP_save*VV,TT)
       DRKZP=MAX(DKZP_save*VV,TT)
    END IF

    DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+DROMG) &
         -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-DROMG))/(2.D0*DROMG)
    DXP =(DISPXR(XP+DRXP,YP,ZP,RKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP-DRXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*DRXP)
    DYP =(DISPXR(XP,YP+DRYP,ZP,RKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP,YP-DRYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*DRYP)
    DZP =(DISPXR(XP,YP,ZP+DRZP,RKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP,YP,ZP-DRZP,RKXP,RKYP,RKZP,OMG))/(2.D0*DRZP)
    DKXP=(DISPXR(XP,YP,ZP,RKXP+DRKXP,RKYP,RKZP,OMG) &
         -DISPXR(XP,YP,ZP,RKXP-DRKXP,RKYP,RKZP,OMG))/(2.D0*DRKXP)
    DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+DRKYP,RKZP,OMG) &
         -DISPXR(XP,YP,ZP,RKXP,RKYP-DRKYP,RKZP,OMG))/(2.D0*DRKYP)
    DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+DRKZP,OMG) &
         -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-DRKZP,OMG))/(2.D0*DRKZP)

    SELECT CASE(model_fdrv_ds)
    CASE(0)
       IF(DOMG.GT.0.D0) THEN
          DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
       ELSE
          DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
       ENDIF
    CASE(1)
       DS=DOMG
    END SELECT

    VX   =-DKXP/DS
    VY   =-DKYP/DS
    VZ   =-DKZP/DS

    VKX  = DXP/DS
    VKY  = DYP/DS
    VKZ  = DZP/DS

    VDU  =-2.D0*ABS(DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)/DS)

    F(1)=VX
    F(2)=VY
    F(3)=VZ
    F(4)=VKX
    F(5)=VKY
    F(6)=VKZ
    IF(UU.LT.0.D0) THEN
       F(7)=0.D0
    ELSE
       F(7)=VDU*UU 
    ENDIF

    DXP_save=VX*DS
    DYP_save=VY*DS
    DZP_save=VZ*DS
    DKXP_save=VKX*DS
    DKYP_save=VKY*DS
    DKZP_save=VKZ*DS

    IF(idebug_wr(12).NE.0) THEN
       WRITE(26,'(A)') '*** idebug_wr(12): wrfdrv'
       WRITE(26,'(A,4ES12.4)') 'x7ds:',X,Y(7),F(7),DS
       WRITE(26,'(A,6ES12.4)') 'y   :',Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
       WRITE(26,'(A,6ES12.4)') 'f   :',F(1),F(2),F(3),F(4),F(5),F(6)
    END IF
    RETURN
  END SUBROUTINE wr_fdrv3

  !  --- slave routine for symplectic method ---

  SUBROUTINE wr_fdrvr(Y,F) 

    USE wrcomm
    USE wrsub,ONLY: DISPXR
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y(6)
    REAL(rkind),INTENT(OUT):: F(6)
    REAL(rkind):: VV,TT,RF,XP,YP,ZP,RKXP,RKYP,RKZP,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ

      VV=DELDER
      TT=DELDER

      RF=wr_nray_status%RF
      XP=Y(1)
      YP=Y(2)
      ZP=Y(3)
      RKXP=Y(4)
      RKYP=Y(5)
      RKZP=Y(6)
      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      RXP=MAX(ABS(XP)*VV,TT)
      RYP=MAX(ABS(YP)*VV,TT)
      RZP=MAX(ABS(ZP)*VV,TT)
      RRKXP=MAX(ABS(RKXP)*VV,TT)
      RRKYP=MAX(ABS(RKYP)*VV,TT)
      RRKZP=MAX(ABS(RKZP)*VV,TT)

      DOMG=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      DXP =(DISPXR(XP+RXP,YP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP-RXP,YP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RXP)
      DYP =(DISPXR(XP,YP+RYP,ZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP-RYP,ZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RYP)
      DZP =(DISPXR(XP,YP,ZP+RZP,RKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP-RZP,RKXP,RKYP,RKZP,OMG))/(2.D0*RZP)
      DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
      DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
      DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
           -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

      IF(DOMG.GT.0.D0) THEN
         DS= SQRT(DKXP**2+DKYP**2+DKZP**2)
      ELSE
         DS=-SQRT(DKXP**2+DKYP**2+DKZP**2)
      ENDIF

      VX   =-DKXP/DS
      VY   =-DKYP/DS
      VZ   =-DKZP/DS

      VKX  = DXP/DS
      VKY  = DYP/DS
      VKZ  = DZP/DS

      F(1)=VX
      F(2)=VY
      F(3)=VZ
      F(4)=VKX
      F(5)=VKY
      F(6)=VKZ
      RETURN
  END SUBROUTINE wr_fdrvr

END MODULE wrfdrv
