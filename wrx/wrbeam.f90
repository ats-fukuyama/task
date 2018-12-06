!   wrbeam.f

MODULE wrbeam

CONTAINS

!     ***** Beam tracing module *****

  SUBROUTINE WR_BEAM

    USE wrcomm
    USE wrsub,ONLY: wrcale,wrnwtn
    IMPLICIT NONE
    REAL(rkind):: Y(NBEQ)
    REAL(4):: TIME1,TIME2
    INTEGER:: IERR,NRAY,NSTP
    REAL(rkind):: ANGZ,ANGPH,SINP2,SINT2

    CALL GUTIME(TIME1)

    IF(MDLWRI.EQ.0) THEN
       WRITE(6,*) '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
       WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
       WRITE(6,'(1PE12.4,0P6F10.3)') &
            RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
       WRITE(6,'(12X,1P5E12.4)') &
            RCURVA,RCURVB,RBRADA,RBRADB,UUI
    ELSEIF(MDLWRI.EQ.1.OR.MDLWRI.EQ.2) THEN
       IF(ABS(RNZI).GT.1.D0) THEN
          ANGZ=0.D0
          ANGPH=0.D0
       ELSE
          ANGZ=ASIN(RNZI)*180.D0/PI
          IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
             ANGZ=0.D0
          ELSE
             ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
          ENDIF
          IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
             ANGPH=0.D0
          ELSE
             ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
          ENDIF
       ENDIF
       IF(MDLWRI.EQ.1) THEN
          WRITE(6,*) '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
          WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
          WRITE(6,'(1PE12.4,0P6F10.3)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
          WRITE(6,'(12X,1P5E12.4)') &
               RCURVA,RCURVB,RBRADA,RBRADB,UUI
       ELSEIF(MDLWRI.EQ.2) THEN
          WRITE(6,*) '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH'
          WRITE(6,*) '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
          WRITE(6,'(1PE12.4,0P6F10.3)') &
               RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
          WRITE(6,'(12X,1P5E12.4)') &
               RCURVA,RCURVB,RBRADA,RBRADB,UUI
       ENDIF
    ENDIF

    DO NRAY=1,NRAYMAX

       IF(MDLWRI.LT.10) then
1         WRITE(6,*) '# NRAY=' ,NRAY
          IF(MDLWRI.EQ.0) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII, &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
          ELSEIF(MDLWRI.EQ.1) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH, &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
          ELSEIF(MDLWRI.EQ.2) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             WRITE(6,*) 'XX MDLWRI=2 IS NOT SUPPORTED YET.'
             GOTO 1
          ENDIF
       ELSE
          RF=RFIN(NRAY)
          RPI=RPIN(NRAY)
          ZPI=ZPIN(NRAY)
          PHII=PHIIN(NRAY)
          RKR0=RKRIN(NRAY)
          RNZI=RNZIN(NRAY)
          RNPHII=RNPHIIN(NRAY)
          ANGZ=ANGZIN(NRAY)
          ANGPH=ANGPHIN(NRAY)
          UUI=UUIN(NRAY)
          RCURVA=RCURVAIN(NRAY)
          RCURVB=RCURVBIN(NRAY)
          RBRADA=RBRADAIN(NRAY)
          RBRADB=RBRADBIN(NRAY)

          IF(MDLWRI.EQ..10)THEN
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
          ELSE
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH'
             WRITE(6,*) &
                  '#                 RCURVA,RCURVB,RBRADA,RBRADB,UU'
             WRITE(6,'(1PE12.4,0P6F10.3)') &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH
             WRITE(6,'(12X,1P5E12.4)') &
                  RCURVA,RCURVB,RBRADA,RBRADB,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
          ENDIF
       ENDIF

       RAYIN(1,NRAY)=RF
       RAYIN(2,NRAY)=RPI
       RAYIN(3,NRAY)=ZPI
       RAYIN(4,NRAY)=PHII
       RAYIN(5,NRAY)=RKR0
!     
       IF(MDLWRI.EQ.0)THEN
          RAYIN(6,NRAY)=RNZI
          RAYIN(7,NRAY)=RNPHII
       ELSE
          RAYIN(6,NRAY)=ANGZ
          RAYIN(7,NRAY)=ANGPH
       ENDIF

       RAYIN(8,NRAY)=UUI

       RKZI  =2.D6*PI*RF*RNZI  /VC
       RKPHII=2.D6*PI*RF*RNPHII/VC
       CALL WRNWTN(IERR)
       IF(IERR.NE.0) GOTO 1200

       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          Y(1)= RPI
          Y(2)= 2.D0*PI*RR*SIN(PHII)
          Y(3)= ZPI
          Y(4)= RKRI
          Y(5)= RKPHII
          Y(6)= RKZI
       ELSE
          Y(1)= RPI*COS(PHII)
          Y(2)= RPI*SIN(PHII)
          Y(3)= ZPI
          Y(4)= RKRI*COS(PHII)-RKPHII*SIN(PHII)
          Y(5)= RKRI*SIN(PHII)+RKPHII*COS(PHII)
          Y(6)= RKZI
       ENDIF

       CALL WRSETY(Y)

       Y(19)= UUI

       CALL WRRKFTB(Y,RAYB,NSTPMAX_NRAY(NRAY))

       DO NSTP=0,NSTPMAX_NRAY(NRAY)
          RAYS(0,NSTP,NRAY)=RAYB( 0,NSTP)
          RAYS(1,NSTP,NRAY)=RAYB( 1,NSTP)
          RAYS(2,NSTP,NRAY)=RAYB( 2,NSTP)
          RAYS(3,NSTP,NRAY)=RAYB( 3,NSTP)
          RAYS(4,NSTP,NRAY)=RAYB( 4,NSTP)
          RAYS(5,NSTP,NRAY)=RAYB( 5,NSTP)
          RAYS(6,NSTP,NRAY)=RAYB( 6,NSTP)
          RAYS(7,NSTP,NRAY)=RAYB(19,NSTP)
          RAYS(8,NSTP,NRAY)=RAYB(20,NSTP)
          RAYRB1(NSTP,NRAY)=RAYB(23,NSTP)
          RAYRB2(NSTP,NRAY)=RAYB(24,NSTP)
       ENDDO
       CALL WRCALE(RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY),NRAY)
       WRITE(6,'(A,F8.4)') &
            '# PABS/PIN=',UUI-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
1200   CONTINUE
    END DO

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

9000 RETURN
  END SUBROUTINE WR_BEAM

! ***** set uo beam tracing initial condition *****

  SUBROUTINE WRSETY(Y)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NBEQ)
    REAL(rkind):: RUM(3,3),VS(3,3),VP(3,3)
    INTEGER:: I,J
    REAL(rkind):: RKABS,RHON,RKNX,RKNY,RKNZ,RKBX,RKBY,RKBZ,RKBABS
    REAL(rkind):: RL0,VS22,VS33,VP22,VP33

    RKABS=SQRT(Y(4)**2+Y(5)**2+Y(6)**2)
    CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)

    RKNX=Y(4)/RKABS
    RKNY=Y(5)/RKABS
    RKNZ=Y(6)/RKABS
    RKBX=RKNY*BNZ-RKNZ*BNY
    RKBY=RKNZ*BNX-RKNX*BNZ
    RKBZ=RKNX*BNY-RKNY*BNX
    RKBABS=SQRT(RKBX**2+RKBY**2+RKBZ**2)

    RUM(1,1)=RKNX
    RUM(1,2)=RKNY
    RUM(1,3)=RKNZ
    RUM(2,1)=RKBX/RKBABS
    RUM(2,2)=RKBY/RKBABS
    RUM(2,3)=RKBZ/RKBABS
    RUM(3,1)=RUM(1,2)*RUM(2,3)-RUM(1,3)*RUM(2,2)
    RUM(3,2)=RUM(1,3)*RUM(2,1)-RUM(1,1)*RUM(2,3)
    RUM(3,3)=RUM(1,1)*RUM(2,2)-RUM(1,2)*RUM(2,1)

    RL0=1.D0/RKABS
    IF(RCURVA.EQ.0.D0)THEN
       VS22=0.D0
    ELSE
       VS22=1.D0/(RL0*RCURVA)
    ENDIF
    IF(RCURVB.EQ.0.D0)THEN
       VS33=0.D0
    ELSE
       VS33=1.D0/(RL0*RCURVB)
    ENDIF
    DO I=1,3
       DO J=1,3
          VS(I,J)=RUM(2,I)*VS22*RUM(2,J) &
                 +RUM(3,I)*VS33*RUM(3,J)
       ENDDO
    ENDDO

    IF(RBRADA.EQ.0.D0)THEN
       VP22=0.D0
    ELSE
       VP22=2.D0/(RBRADA*RBRADA)
    ENDIF
    IF(RBRADB.EQ.0.D0)THEN
       VP33=0.D0
    ELSE
       VP33=2.D0/(RBRADB*RBRADB)
    ENDIF
    DO I=1,3
       DO J=1,3
          VP(I,J)=RUM(2,I)*VP22*RUM(2,J) &
                 +RUM(3,I)*VP33*RUM(3,J)
       ENDDO
    ENDDO

    Y( 7)=VS(1,1)
    Y( 8)=VS(1,2)
    Y( 9)=VS(1,3)
    Y(10)=VS(2,2)
    Y(11)=VS(2,3)
    Y(12)=VS(3,3)

    Y(13)=VP(1,1)
    Y(14)=VP(1,2)
    Y(15)=VP(1,3)
    Y(16)=VP(2,2)
    Y(17)=VP(2,3)
    Y(18)=VP(3,3)

    RETURN
  END SUBROUTINE WRSETY

!************************************************************************

  SUBROUTINE WRGETY(Y,RC1,RC2,RB1,RB2,RTH1,RTH2,YY)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y(NBEQ)
    REAL(rkind),INTENT(OUT):: YY(27)
    REAL(rkind),INTENT(OUT):: RC1,RC2,RB1,RB2,RTH1,RTH2
    REAL(rkind):: DFK(3),RUM(3,3),VS(3,3),VP(3,3),WS(3,3),WP(3,3)
    INTEGER:: I,J,K,L
    REAL(rkind):: RHON,RKNX,RKNY,RKNZ,RKBX,RKBY,RKBZ,RKABS,RKBABS
    REAL(rkind):: SD,RSLAM1,RSLAM2,BD,RBLAM1,RBLAM2,RL0
    REAL(rkind):: VV,TT,XP,YP,ZP,RKXP,RKYP,RKZP
    REAL(rkind):: OMG,ROMG,DKXP,DKYP,DKZP,DOMG
    REAL(rkind):: F000P00,F000M00,F0000P0,F0000M0,F00000P,F00000M
    REAL(rkind):: VGX,VGY,VGZ,RVGABS

    CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)

    RKNX=Y(4)/RKABS
    RKNY=Y(5)/RKABS
    RKNZ=Y(6)/RKABS
    RKBX=RKNY*BNZ-RKNZ*BNY
    RKBY=RKNZ*BNX-RKNX*BNZ
    RKBZ=RKNX*BNY-RKNY*BNX
    RKBABS=SQRT(RKBX**2+RKBY**2+RKBZ**2)

    RUM(1,1)=RKNX
    RUM(1,2)=RKNY
    RUM(1,3)=RKNZ
    RUM(2,1)=RKBX/RKBABS
    RUM(2,2)=RKBY/RKBABS
    RUM(2,3)=RKBZ/RKBABS
    RUM(3,1)=RUM(1,2)*RUM(2,3)-RUM(1,3)*RUM(2,2)
    RUM(3,2)=RUM(1,3)*RUM(2,1)-RUM(1,1)*RUM(2,3)
    RUM(3,3)=RUM(1,1)*RUM(2,2)-RUM(1,2)*RUM(2,1)

    VS(1,1)=Y( 7)
    VS(1,2)=Y( 8)
    VS(1,3)=Y( 9)
    VS(2,1)=Y( 8)
    VS(2,2)=Y(10)
    VS(2,3)=Y(11)
    VS(3,1)=Y( 9)
    VS(3,2)=Y(11)
    VS(3,3)=Y(12)

    VP(1,1)=Y(13)
    VP(1,2)=Y(14)
    VP(1,3)=Y(15)
    VP(2,1)=Y(14)
    VP(2,2)=Y(16)
    VP(2,3)=Y(17)
    VP(3,1)=Y(15)
    VP(3,2)=Y(17)
    VP(3,3)=Y(18)

    DO I=1,3
       DO J=1,3
          WS(I,J)=0.D0
          WP(I,J)=0.D0
          DO K=1,3
             DO L=1,3
                WS(I,J)=WS(I,J)+RUM(I,K)*VS(K,L)*RUM(J,L)
                WP(I,J)=WP(I,J)+RUM(I,K)*VP(K,L)*RUM(J,L)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    YY(1)=WS(1,1)
    YY(2)=WS(1,2)
    YY(3)=WS(1,3)
    YY(4)=WS(2,2) 
    YY(5)=WS(2,3)
    YY(6)=WS(3,3)

    YY(7)=WP(1,1)
    YY(8)=WP(1,2)
    YY(9)=WP(1,3)
    YY(10)=WP(2,2) 
    YY(11)=WP(2,3)
    YY(12)=WP(3,3)

    YY(13)=RUM(1,1)
    YY(14)=RUM(1,2)
    YY(15)=RUM(1,3)
    YY(16)=RUM(2,1)
    YY(17)=RUM(2,2)
    YY(18)=RUM(2,3)
    YY(19)=RUM(3,1)
    YY(20)=RUM(3,2)
    YY(21)=RUM(3,3)

!     ***********   MODIFICATION BY NISHINA *********
!     --- calculate eigenvalue ---

      SD=(WS(2,2)-WS(3,3))**2+4.D0*WS(2,3)**2
      RSLAM1=0.5D0*(WS(2,2)+WS(3,3)+SQRT(SD))
      RSLAM2=0.5D0*(WS(2,2)+WS(3,3)-SQRT(SD))
      BD=(WP(2,2)-WP(3,3))**2+4.D0*WP(2,3)**2
      RBLAM1=0.5D0*(WP(2,2)+WP(3,3)+SQRT(BD))
      RBLAM2=0.5D0*(WP(2,2)+WP(3,3)-SQRT(BD))

!      WRITE(6,'(A,1P2E12.4)') 'SD,BD        =',SD,BD
!      WRITE(6,'(A,1P2E12.4)') 'RSLAM1,RSLAM2=',RSLAM1,RSLAM2
!      WRITE(6,'(A,1P2E12.4)') 'RBLAM1,RBLAM2=',RBLAM1,RBLAM2

!     calculate curvature radius

      RL0=1.D0/RKABS
      IF(RSLAM1.EQ.0.D0)THEN
         RC1=0.D0
      ELSE
         RC1=1.D0/(RL0*RSLAM1)
      ENDIF
      IF(RSLAM2.EQ.0.D0)THEN
         RC2=0.D0
      ELSE
         RC2=1.D0/(RL0*RSLAM2)
      ENDIF

!     calculating beam radius

      IF(RBLAM1.LE.0.D0)THEN
         RB1=0.D0
      ELSE
         RB1=SQRT(2.D0/RBLAM1)
      ENDIF
      IF(RBLAM2.LE.0.D0)THEN
         RB2=0.D0
      ELSE
         RB2=SQRT(2.D0/RBLAM2)
      ENDIF

!      WRITE(6,'(A,1P2E12.4)') 'RC1,RC2      =',RC1,RC2
!      WRITE(6,'(A,1P2E12.4)') 'RB1,RB2      =',RB1,RB2

!     calculate angle

      RTH1=ATAN2(-WP(2,3),WP(2,2)-RBLAM1)*180.D0/PI
      RTH2=ATAN2(-WP(2,3),WP(2,2)-RBLAM2)*180.D0/PI
      IF(RTH1.LT.0.D0) RTH1=RTH1+180.D0
      IF(RTH2.LT.0.D0) RTH2=RTH2+180.D0

!      WRITE(6,'(A,1P4E12.4)') 'RB1,RB2,RTH1,RTH2=',RB1,RB2,RTH1,RTH2


!     -----------------CALC VG---------------------------------------

      VV=DELDER
      TT=DELDER

      XP=Y(1)
      YP=Y(2)
      ZP=Y(3)
      RKXP=Y(4)
      RKYP=Y(5)
      RKZP=Y(6)

      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      DKXP=MAX(ABS(RKXP)*VV,TT)
      DKYP=MAX(ABS(RKYP)*VV,TT)
      DKZP=MAX(ABS(RKZP)*VV,TT)   

      DOMG=(DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
           -DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)

      F000P00=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F000M00=DISPBXR(XP,YP,ZP,RKXP-DKXP,RKYP,RKZP,OMG)
      F0000P0=DISPBXR(XP,YP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F0000M0=DISPBXR(XP,YP,ZP,RKXP,RKYP-DKYP,RKZP,OMG)
      F00000P=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F00000M=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP-DKZP,OMG)

      DFK(1)=(F000P00-F000M00)/(2.D0*DKXP)
      DFK(2)=(F0000P0-F0000M0)/(2.D0*DKYP)
      DFK(3)=(F00000P-F00000M)/(2.D0*DKZP)

      VGX=DFK(1)/DOMG
      VGY=DFK(2)/DOMG
      VGZ=DFK(3)/DOMG
      RVGABS=SQRT(VGX**2+VGY**2+VGZ**2)

      YY(22)=RKXP/RKABS
      YY(23)=RKYP/RKABS
      YY(24)=RKZP/RKABS
      YY(25)=VGX/RVGABS
      YY(26)=VGY/RVGABS
      YY(27)=VGZ/RVGABS

      RETURN
  END SUBROUTINE WRGETY

!************************************************************************

  SUBROUTINE WRRKFTB(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NBEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NBVAR,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NBEQ),WORK(NBEQ,2),YY(27),F(NBEQ)
    INTEGER:: NSTPLIM,NSTP,I,ID
    REAL(rkind):: X0,XE,YPRE,RL,PHIL,ZL,RKRL,X,RHON

      NSTPLIM=INT(SMAX/DELS)

      NSTP=0
      X0 = 0.D0
      XE = DELS     
      YN(0,NSTP)=X0
      DO I=1,NBEQ
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(20,NSTP)=0.D0
      CALL WRGETY(Y,YN(21,NSTP),YN(22,NSTP),YN(23,NSTP),YN(24,NSTP), &
                                            YN(25,NSTP),YN(26,NSTP),YY)
      DO I=27,53
         YN(I,NSTP)=0.D0
      ENDDO

      CALL WRFDRVB (X,Y,F)

      DO NSTP=1,NSTPLIM
         YPRE=Y(19)
         CALL ODERK(NBEQ,WRFDRVB,X0,XE,10,Y,YM,WORK)
         YN(0,NSTP)=XE
         DO I=1,19
           YN(I,NSTP)=YM(I)
         ENDDO
         YN(20,NSTP)=YPRE-YM(19)

         CALL WRGETY(YM,YN(21,NSTP),YN(22,NSTP),YN(23,NSTP),YN(24,NSTP), &
                                                YN(25,NSTP),YN(26,NSTP),YY)

         DO I=27,53
           YN(I,NSTP)=YY(I-26)
         ENDDO

         CALL WRFDRVB (X,Y,F)

         IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
            RL  =YM(1)
            PHIL=ASIN(YM(2)/(2.D0*PI*RR))
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE
            RL  =SQRT(YM(1)**2+YM(2)**2)
            PHIL=ATAN2(YM(2),YM(1))
            ZL  =YM(3)
            RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         ENDIF

         IF(MDLWRW.GE.1) THEN
            ID=0
            SELECT CASE(MDLWRW)
            CASE(1)
               ID=1
            CASE(2)
               IF(MOD(NSTP-1,10).EQ.0) ID=1
            CASE(3)
               IF(MOD(NSTP-1,100).EQ.0) ID=1
            CASE(4)
               IF(MOD(NSTP-1,1000).EQ.0) ID=1
            CASE(5)
               IF(MOD(NSTP-1,10000).EQ.0) ID=1
            END SELECT
            IF(ID.EQ.1) THEN
               WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(19),YN(20,NSTP)
               WRITE(6,'(11X,1P6E11.3)') YN(21,NSTP),YN(22,NSTP),YN(23,NSTP), &
                                         YN(24,NSTP),YN(25,NSTP),YN(26,NSTP)
            END IF
         ENDIF

         DO I=1,19
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(19).LT.UUMIN) THEN
            NNSTP = NSTP
            WRITE(6,*) '--- Absorbed ---'
            GOTO 10
         ENDIF         
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA) THEN
            NNSTP = NSTP
            WRITE(6,*) '--- Out of bounds ---'
            GOTO 10
         ENDIF
      ENDDO
      NNSTP=NSTPLIM
!     
 10   IF(YN(19,NNSTP).LT.0.D0) THEN
         YN(19,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) THEN
            WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(19),YN(20,NNSTP)
            WRITE(6,'(11X,1P6E11.3)') YN(21,NNSTP),YN(22,NNSTP),YN(23,NNSTP), &
                                      YN(24,NNSTP),YN(25,NNSTP),YN(26,NNSTP)
         END IF
      ENDIF

      RETURN
  END SUBROUTINE WRRKFTB

!************************************************************************

  SUBROUTINE WRFDRVB(X,Y,F) 

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y(NBEQ)
    REAL(rkind),INTENT(OUT):: F(NBEQ)
    REAL(rkind):: DFR(3),DFK(3),DFRR(3,3),DFKK(3,3),DFRK(3,3),S(3,3),P(3,3)
    INTEGER:: I,J
    REAL(rkind):: VV,TT,XP,YP,ZP,RKXP,RKYP,RKZP,UU
    REAL(rkind):: OMG,ROMG,DRXP,DRYP,DRZP,DKXP,DKYP,DKZP,DOMG,DS,VDU
    REAL(rkind):: F000000,FP00000,FM00000,F0P0000,F0M0000,F00P000,F00M000
    REAL(rkind):: F000P00,F000M00,F0000P0,F0000M0,F00000P,F00000M
    REAL(rkind):: FPP0000,FP0P000,FP00P00,FP000P0,FP0000P
    REAL(rkind):: F0PP000,F0P0P00,F0P00P0,F0P000P
    REAL(rkind):: F00PP00,F00P0P0,F00P00P
    REAL(rkind):: F000PP0,F000P0P
    REAL(rkind):: F0000PP

      VV=DELDER
      TT=DELDER

      XP=Y(1)
      YP=Y(2)
      ZP=Y(3)
      RKXP=Y(4)
      RKYP=Y(5)
      RKZP=Y(6)
      UU=Y(19)

      OMG=2.D6*PI*RF
      ROMG=MAX(ABS(OMG)*VV,TT)
      DRXP=MAX(ABS(XP)*VV,TT)
      DRYP=MAX(ABS(YP)*VV,TT)
      DRZP=MAX(ABS(ZP)*VV,TT)
      DKXP=MAX(ABS(RKXP)*VV,TT)
      DKYP=MAX(ABS(RKYP)*VV,TT)
      DKZP=MAX(ABS(RKZP)*VV,TT)   

      DOMG=(DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG+ROMG) &
           -DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG-ROMG))/(2.D0*ROMG)
      F000000=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

      FP00000=DISPBXR(XP+DRXP,YP,ZP,RKXP,RKYP,RKZP,OMG)
      FM00000=DISPBXR(XP-DRXP,YP,ZP,RKXP,RKYP,RKZP,OMG)
      F0P0000=DISPBXR(XP,YP+DRYP,ZP,RKXP,RKYP,RKZP,OMG)
      F0M0000=DISPBXR(XP,YP-DRYP,ZP,RKXP,RKYP,RKZP,OMG)
      F00P000=DISPBXR(XP,YP,ZP+DRZP,RKXP,RKYP,RKZP,OMG)
      F00M000=DISPBXR(XP,YP,ZP-DRZP,RKXP,RKYP,RKZP,OMG)
      F000P00=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F000M00=DISPBXR(XP,YP,ZP,RKXP-DKXP,RKYP,RKZP,OMG)
      F0000P0=DISPBXR(XP,YP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F0000M0=DISPBXR(XP,YP,ZP,RKXP,RKYP-DKYP,RKZP,OMG)
      F00000P=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F00000M=DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP-DKZP,OMG)

      FPP0000=DISPBXR(XP+DRXP,YP+DRYP,ZP,RKXP,RKYP,RKZP,OMG)
      FP0P000=DISPBXR(XP+DRXP,YP,ZP+DRZP,RKXP,RKYP,RKZP,OMG)
      FP00P00=DISPBXR(XP+DRXP,YP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      FP000P0=DISPBXR(XP+DRXP,YP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      FP0000P=DISPBXR(XP+DRXP,YP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F0PP000=DISPBXR(XP,YP+DRYP,ZP+DRZP,RKXP,RKYP,RKZP,OMG)
      F0P0P00=DISPBXR(XP,YP+DRYP,ZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F0P00P0=DISPBXR(XP,YP+DRYP,ZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F0P000P=DISPBXR(XP,YP+DRYP,ZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F00PP00=DISPBXR(XP,YP,ZP+DRZP,RKXP+DKXP,RKYP,RKZP,OMG)
      F00P0P0=DISPBXR(XP,YP,ZP+DRZP,RKXP,RKYP+DKYP,RKZP,OMG)
      F00P00P=DISPBXR(XP,YP,ZP+DRZP,RKXP,RKYP,RKZP+DKZP,OMG)
      F000PP0=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP+DKYP,RKZP,OMG)
      F000P0P=DISPBXR(XP,YP,ZP,RKXP+DKXP,RKYP,RKZP+DKZP,OMG)
      F0000PP=DISPBXR(XP,YP,ZP,RKXP,RKYP+DKYP,RKZP+DKZP,OMG)

      DFR(1)=(FP00000-FM00000)/(2.D0*DRXP)
      DFR(2)=(F0P0000-F0M0000)/(2.D0*DRYP)
      DFR(3)=(F00P000-F00M000)/(2.D0*DRZP)
      DFK(1)=(F000P00-F000M00)/(2.D0*DKXP)
      DFK(2)=(F0000P0-F0000M0)/(2.D0*DKYP)
      DFK(3)=(F00000P-F00000M)/(2.D0*DKZP)

      DFRR(1,1)=(FP00000-F000000-F000000+FM00000)/(DRXP*DRXP)
      DFRR(1,2)=(FPP0000-FP00000-F0P0000+F000000)/(DRXP*DRYP)
      DFRR(1,3)=(FP0P000-FP00000-F00P000+F000000)/(DRXP*DRZP)
      DFRR(2,1)=DFRR(1,2)
      DFRR(2,2)=(F0P0000-F000000-F000000+F0M0000)/(DRYP*DRYP)
      DFRR(2,3)=(F0PP000-F0P0000-F00P000+F000000)/(DRYP*DRZP)
      DFRR(3,1)=DFRR(1,3)
      DFRR(3,2)=DFRR(2,3)
      DFRR(3,3)=(F00P000-F000000-F000000+F00M000)/(DRZP*DRZP)

!      WRITE(6,'(A,1P3E12.4)') 'DRP  :',DRXP,DRYP,DRZP
!      WRITE(6,'(A,1P3E12.4)') 'DKP  :',DKXP,DKYP,DKZP
!      WRITE(6,'(A,1P3E12.4)') 'DFR  :',DFR(1),DFR(2),DFR(3)
!      WRITE(6,'(A,1P3E12.4)') 'DFK  :',DFK(1),DFK(2),DFK(3)
!      WRITE(6,'(A,1P3E12.4)') 'DFRR :',DFRR(1,1),DFRR(1,2),DFRR(1,3)
!      WRITE(6,'(A,1P3E12.4)') 'DFRR :',DFRR(2,1),DFRR(2,2),DFRR(2,3)
!      WRITE(6,'(A,1P6E12.4)') 'DFRR :',DFRR(3,1),DFRR(3,2),DFRR(3,3)

      DFKK(1,1)=(F000P00-F000000-F000000+F000M00)/(DKXP*DKXP)
      DFKK(1,2)=(F000PP0-F000P00-F0000P0+F000000)/(DKXP*DKYP)
      DFKK(1,3)=(F000P0P-F000P00-F00000P+F000000)/(DKXP*DKZP)
      DFKK(2,1)=DFKK(1,2)
      DFKK(2,2)=(F0000P0-F000000-F000000+F0000M0)/(DKYP*DKYP)
      DFKK(2,3)=(F0000PP-F0000P0-F00000P+F000000)/(DKYP*DKZP)
      DFKK(3,1)=DFKK(1,3)
      DFKK(3,2)=DFKK(2,3)
      DFKK(3,3)=(F00000P-F000000-F000000+F00000M)/(DKZP*DKZP)

!      WRITE(6,'(A,1P3E12.4)') 'DFKK :',DFKK(1,1),DFKK(1,2),DFKK(1,3)
!      WRITE(6,'(A,1P3E12.4)') 'DFKK :',DFKK(2,1),DFKK(2,2),DFKK(2,3)
!      WRITE(6,'(A,1P3E12.4)') 'DFKK :',DFKK(3,1),DFKK(3,2),DFKK(3,3)

      DFRK(1,1)=(FP00P00-FP00000-F000P00+F000000)/(DRXP*DKXP)
      DFRK(1,2)=(FP000P0-FP00000-F0000P0+F000000)/(DRXP*DKYP)
      DFRK(1,3)=(FP0000P-FP00000-F00000P+F000000)/(DRXP*DKZP)
      DFRK(2,1)=(F0P0P00-F0P0000-F000P00+F000000)/(DRYP*DKXP)
      DFRK(2,2)=(F0P00P0-F0P0000-F0000P0+F000000)/(DRYP*DKYP)
      DFRK(2,3)=(F0P000P-F0P0000-F00000P+F000000)/(DRYP*DKZP)
      DFRK(3,1)=(F00PP00-F00P000-F000P00+F000000)/(DRZP*DKXP)
      DFRK(3,2)=(F00P0P0-F00P000-F0000P0+F000000)/(DRZP*DKYP)
      DFRK(3,3)=(F00P00P-F00P000-F00000P+F000000)/(DRZP*DKZP)

!      WRITE(6,'(A,1P3E12.4)') 'DFRK :',DFRK(1,1),DFRK(1,2),DFRK(1,3)
!      WRITE(6,'(A,1P3E12.4)') 'DFRK :',DFRK(2,1),DFRK(2,2),DFRK(2,3)
!      WRITE(6,'(A,1P3E12.4)') 'DFRK :',DFRK(3,1),DFRK(3,2),DFRK(3,3)

      IF(DOMG.GT.0.D0) THEN
         DS=-1.D0/SQRT(DFK(1)**2+DFK(2)**2+DFK(3)**2)
      ELSE
         DS= 1.D0/SQRT(DFK(1)**2+DFK(2)**2+DFK(3)**2)
      ENDIF

      S(1,1)=Y( 7)
      S(1,2)=Y( 8)
      S(1,3)=Y( 9)
      S(2,1)=Y( 8)
      S(2,2)=Y(10)
      S(2,3)=Y(11)
      S(3,1)=Y( 9)
      S(3,2)=Y(11)
      S(3,3)=Y(12)

      P(1,1)=Y(13)
      P(1,2)=Y(14)
      P(1,3)=Y(15)
      P(2,1)=Y(14)
      P(2,2)=Y(16)
      P(2,3)=Y(17)
      P(3,1)=Y(15)
      P(3,2)=Y(17)
      P(3,3)=Y(18)

      F( 1)= DFK(1)*DS
      F( 2)= DFK(2)*DS
      F( 3)= DFK(3)*DS
      F( 4)=-DFR(1)*DS
      F( 5)=-DFR(2)*DS
      F( 6)=-DFR(3)*DS

      F( 7)=-DFRR(1,1)*DS
      F( 8)=-DFRR(1,2)*DS
      F( 9)=-DFRR(1,3)*DS
      F(10)=-DFRR(2,2)*DS
      F(11)=-DFRR(2,3)*DS
      F(12)=-DFRR(3,3)*DS

      DO I=1,3
         F( 7)=F( 7)-(DFRK(1,I)*S(1,I)+DFRK(1,I)*S(1,I))*DS
         F( 8)=F( 8)-(DFRK(2,I)*S(1,I)+DFRK(1,I)*S(2,I))*DS
         F( 9)=F( 9)-(DFRK(3,I)*S(1,I)+DFRK(1,I)*S(3,I))*DS
         F(10)=F(10)-(DFRK(2,I)*S(2,I)+DFRK(2,I)*S(2,I))*DS
         F(11)=F(11)-(DFRK(3,I)*S(2,I)+DFRK(2,I)*S(3,I))*DS
         F(12)=F(12)-(DFRK(3,I)*S(3,I)+DFRK(3,I)*S(3,I))*DS
      ENDDO

      DO I=1,3
      DO J=1,3
         F( 7)=F( 7) &
              -(DFKK(I,J)*S(1,I)*S(1,J)-DFKK(I,J)*P(1,I)*P(1,J))*DS
         F( 8)=F( 8) &
              -(DFKK(I,J)*S(1,I)*S(2,J)-DFKK(I,J)*P(1,I)*P(2,J))*DS
         F( 9)=F( 9) &
              -(DFKK(I,J)*S(1,I)*S(3,J)-DFKK(I,J)*P(1,I)*P(3,J))*DS
         F(10)=F(10) &
              -(DFKK(I,J)*S(2,I)*S(2,J)-DFKK(I,J)*P(2,I)*P(2,J))*DS
         F(11)=F(11) &
              -(DFKK(I,J)*S(2,I)*S(3,J)-DFKK(I,J)*P(2,I)*P(3,J))*DS
         F(12)=F(12) &
              -(DFKK(I,J)*S(3,I)*S(3,J)-DFKK(I,J)*P(3,I)*P(3,J))*DS
      ENDDO
      ENDDO

      F(13)=0.D0
      F(14)=0.D0
      F(15)=0.D0
      F(16)=0.D0
      F(17)=0.D0
      F(18)=0.D0

      DO I=1,3
         F(13)=F(13)-(DFRK(1,I)*P(1,I)+DFRK(1,I)*P(1,I))*DS
         F(14)=F(14)-(DFRK(2,I)*P(1,I)+DFRK(1,I)*P(2,I))*DS
         F(15)=F(15)-(DFRK(3,I)*P(1,I)+DFRK(1,I)*P(3,I))*DS
         F(16)=F(16)-(DFRK(2,I)*P(2,I)+DFRK(2,I)*P(2,I))*DS
         F(17)=F(17)-(DFRK(3,I)*P(2,I)+DFRK(2,I)*P(3,I))*DS
         F(18)=F(18)-(DFRK(3,I)*P(3,I)+DFRK(3,I)*P(3,I))*DS
      ENDDO

      DO I=1,3
      DO J=1,3
         F(13)=F(13) &
              -(DFKK(I,J)*S(1,I)*P(1,J)+DFKK(J,I)*S(1,I)*P(1,J))*DS
         F(14)=F(14) &
              -(DFKK(I,J)*S(1,I)*P(2,J)+DFKK(J,I)*S(2,I)*P(1,J))*DS
         F(15)=F(15) &
              -(DFKK(I,J)*S(1,I)*P(3,J)+DFKK(J,I)*S(3,I)*P(1,J))*DS
         F(16)=F(16) &
              -(DFKK(I,J)*S(2,I)*P(2,J)+DFKK(J,I)*S(2,I)*P(2,J))*DS
         F(17)=F(17) &
              -(DFKK(I,J)*S(2,I)*P(3,J)+DFKK(J,I)*S(3,I)*P(2,J))*DS
         F(18)=F(18) &
              -(DFKK(I,J)*S(3,I)*P(3,J)+DFKK(J,I)*S(3,I)*P(3,J))*DS
      ENDDO
      ENDDO

!      VDU  =-2*DISPBXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)*DS
      VDU  =-2.D0*ABS(DISPBXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)*DS)

      IF(UU.LT.0.D0) THEN
         F(19)=0.D0
      ELSE
         F(19)=VDU*UU 
      ENDIF

!      WRITE(6,'(A,1P6E12.4)') 'X:',X
!      WRITE(6,'(A,1P6E12.4)') 'Y :',Y( 1),Y( 2),Y( 3),Y( 4),Y( 5),Y( 6)
!      WRITE(6,'(A,1P6E12.4)') 'Y :',Y( 7),Y( 8),Y( 9),Y(10),Y(11),Y(12)
!      WRITE(6,'(A,1P6E12.4)') 'Y :',Y(13),Y(14),Y(15),Y(16),Y(17),Y(18)
!      WRITE(6,'(A,1P6E12.4)') 'F :',F( 1),F( 2),F( 3),F( 4),F( 5),F( 6)
!      WRITE(6,'(A,1P6E12.4)') 'F :',F( 7),F( 8),F( 9),F(10),F(11),F(12)
!      WRITE(6,'(A,1P6E12.4)') 'F :',F(13),F(14),F(15),F(16),F(17),F(18)
!      CALL GUFLSH

      RETURN
  END SUBROUTINE WRFDRVB

! ***** Dispersion relation *****

  FUNCTION DISPBXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    REAL(rkind):: DISPBXR
    REAL(rkind):: MODELPS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CWW,CWC,CF

      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         IF(MODELP(NS).GE.10.AND.MODELP(NS).LT.30) MODELP(NS)=0
         IF(MODELP(NS).GE.30.AND.MODELP(NS).LT.50) MODELP(NS)=4
         IF(MODELP(NS).GE.50.AND.MODELP(NS).LT.60) MODELP(NS)=5
      ENDDO

      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
         CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
      ENDDO

      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
      ENDDO

      DISPBXR=DBLE(CF)
      RETURN
  END FUNCTION DISPBXR

!***********************************************************************

  FUNCTION DISPBXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    REAL(rkind):: DISPBXI
    REAL(rkind):: MODELPS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CWW,CWC,CF

      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         IF(MODELP(NS).GE.10.AND.MODELP(NS).LT.30) MODELP(NS)=0
         IF(MODELP(NS).GE.30.AND.MODELP(NS).LT.50) MODELP(NS)=4
         IF(MODELP(NS).GE.50.AND.MODELP(NS).LT.60) MODELP(NS)=5
      ENDDO

      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
         CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
      ENDDO

      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
      ENDDO

      DISPBXI=DIMAG(CF)
      RETURN
  END FUNCTION DISPBXI
END MODULE wrbeam
