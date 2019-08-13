!   wrexecr.f90

MODULE wrexecr

  PRIVATE
  PUBLIC wr_execr

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE wr_execr(NRAY,IERR)

    USE wrcomm
    USE wrsub,ONLY: wrcale
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NRAY
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: Y(NEQ)
    INTEGER:: NSTP

    IERR=0
    IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
       Y(1)= RPI
       Y(2)= PHII
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
    Y(7)= UUI
    IF(MDLWRQ.EQ.0) THEN
       CALL WRRKFT(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.1) THEN
       CALL WRRKFT_WITHD0(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.2) THEN
       CALL WRRKFT_WITHMC(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.3) THEN
       CALL WRRKFT_ODE(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.4) THEN
       CALL WRSYMP(Y,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSE
       WRITE(6,*) 'XX WRCALC: unknown MDLWRQ =', MDLWRQ
       IERR=1
       RETURN
    ENDIF
    DO NSTP=0,NSTPMAX_NRAY(NRAY)
       RAYRB1(NSTP,NRAY)=0.D0
       RAYRB2(NSTP,NRAY)=0.D0
    END DO

    CALL WRCALE(RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY),NRAY)

    RETURN
  END SUBROUTINE wr_execr

!  --- original Runge-Kutta method ---

  SUBROUTINE WRRKFT(Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ)
    INTEGER:: INIT,NSTPLIM,NSTP,I,IOX,ID
    REAL(rkind):: X0,XE,OMG,PSIN,RL,PHIL,ZL,RKRL,Y7,RHON
    REAL(rkind):: RNPHI_IDEI,DELTA,RKPARA,RKPERP

    X0 = 0.D0
    XE = DELS
    INIT = 1
    NSTPLIM=INT(SMAX/DELS)

    NSTP=0
    YN(0,NSTP)=X0
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0
    OMG=2.D6*PI*RF
    IOX=0

    DO NSTP = 1,NSTPLIM
       Y7=Y(7)
       CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=Y7-YM(7)

       CALL PL_MAG_OLD(YM(1),YM(2),YM(3),PSIN)
       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          RL  =YM(1)
          PHIL=ASIN(YM(2)/(2.D0*PI*RR))
          ZL  =YM(3)
          RKRL=YM(4)
       ELSE IF(MODELG.EQ.11) THEN
          RL  =YM(1)
          PHIL=YM(2)
          ZL  =YM(3)
          RKRL=YM(4)
       ELSE
          RL  =SQRT(YM(1)**2+YM(2)**2)
          PHIL=ATAN2(YM(2),YM(1))
          ZL  =YM(3)
          RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
       ENDIF
       RNPHI_IDEI= -YM(4)*VC/OMG*SIN(PHIL) + YM(5)*VC/OMG*COS(PHIL)
       DELTA=DISPXR( YM(1), YM(2), YM(3), YM(4), YM(5), YM(6), OMG)
       RKPARA=YM(4)*BNX+YM(5)*BNY+YM(6)*BNZ
       RKPERP=SQRT((YM(4)*YM(4)+YM(5)*YM(5)+YM(6)*YM(6))-RKPARA**2)

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
          IF(ID.EQ.1) &
               WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
       ENDIF

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS
       IF(Y(7).LT.UUMIN) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA*RKAP) THEN
             NNSTP = NSTP
             GOTO 11
          ENDIF
       ENDIF
    END DO
    NNSTP=NSTPLIM
!     
11  IF(YN(7,NNSTP).LT.0.D0) THEN
       YN(7,NNSTP)=0.D0
    ENDIF
    IF(MDLWRW.GE.1) THEN
       IF(ID.EQ.0) &
            WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
    ENDIF

    RETURN
  END SUBROUTINE WRRKFT

!  --- original Runge-Kutta method with correction for D=0 ---

  SUBROUTINE WRRKFT_WITHD0(Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,WRMODNWTN
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3)
    INTEGER:: INIT,NSTPLIM,NSTP,I,ID
    REAL(rkind):: X0,XE,OMG,Y7,DELTA,RHON,RL,PHIL,ZL,RKRL
    REAL(rkind):: RNPHI_IDEI,RKPARA,RKPERP

      X0 = 0.D0
      XE = DELS
      INIT = 1
      NSTPLIM=INT(SMAX/DELS)
      NSTP=0
      YN(0,NSTP)=X0
      DO I=1,7
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(8,NSTP)=0.D0
      OMG=2.D6*PI*RF

      DO NSTP = 1,NSTPLIM
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),OMG)
         IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,YK)
            DO I=1,3
               YM(I+3) = YK(I)
            END DO
         END IF 
         YN(0,NSTP)=XE
         DO I=1,7
            YN(I,NSTP)=YM(I)
         ENDDO
         YN(8,NSTP)=Y7-YM(7)

         CALL PL_MAG_OLD(YM(1),YM(2),YM(3),RHON)
         IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
            RL  =YM(1)
            IF(RR.EQ.0.D0) THEN
               PHIL=YM(2)
            ELSE
               PHIL=ASIN(YM(2)/(2.D0*PI*RR))
            ENDIF
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE IF(MODELG.EQ.11) THEN
            RL  =YM(1)
            PHIL=YM(2)
            ZL  =YM(3)
            RKRL=YM(4)
         ELSE
            RL  =SQRT(YM(1)**2+YM(2)**2)
            PHIL=ATAN2(YM(2),YM(1))
            ZL  =YM(3)
            RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         ENDIF
         RNPHI_IDEI= -YM(4)*VC/OMG*SIN(PHIL) + YM(5)*VC/OMG*COS(PHIL)
         DELTA=DISPXR( YM(1), YM(2), YM(3), YM(4), YM(5), YM(6), OMG )
         RKPARA=YM(4)*BNX+YM(5)*BNY+YM(6)*BNZ
         RKPERP=SQRT((YM(4)*YM(4)+YM(5)*YM(5)+YM(6)*YM(6))-RKPARA**2)

!IDEI	 WRITE(6,*) 'BX=',BNX*BABS,'BY=',BNY*BABS,'BZ=',BNZ*BABS
!	 WRITE(6,6001) YN(0,IT),RL,ZL,PHIL,YN(1,IT), &
!                      YN(2,IT),YN(3,IT),YN(6,IT)*VC/OMG,YN(7,IT), &
!      	               DELTA,RKPARA*VC/OMG,RKPERP*VC/OMG,RKRL, &
!      	               RNPHI_IDEI
! 6001    FORMAT(1H ,1P14E13.5)

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
            IF(ID.EQ.1) &
                 WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
         ENDIF

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
         IF(MODELG.LE.10) THEN
            CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
            IF(RHON.GT.RB/RA*RKAP) THEN
               NNSTP = NSTP
               GOTO 11
            ENDIF        
         END IF
      END DO
      NNSTP=NSTPLIM
!     
 11   IF(YN(7,NNSTP).LT.0.D0) THEN
         YN(7,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) &
              WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NNSTP)
      ENDIF

      RETURN
  END SUBROUTINE WRRKFT_WITHD0

!  --- Runge-Kutta method using ODE library ---

  SUBROUTINE WRRKFT_ODE(Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ)
    REAL(rkind):: X0,XE,RL,PHIL,ZL,Y7,RHON,RKRL
    INTEGER:: NSTPLIM,NSTP,I,ID

      X0 = 0.D0
      XE = DELS     
      NSTPLIM=INT(SMAX/DELS)
      NSTP=0
      YN(0,NSTP)=X0
      DO I=1,7
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(8,NSTP)=0.D0

      DO NSTP = 1,NSTPLIM
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         YN(0,NSTP)=XE
         DO I=1,7
            YN(I,NSTP)=YM(I)
         ENDDO
         IF(YM(7).GT.0.D0) THEN
            YN(8,NSTP)=Y7-YM(7)
         ELSE
            YN(8,NSTP)=Y7
         ENDIF

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
            IF(ID.EQ.1) &
                 WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
         ENDIF

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA*RKAP) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
      END DO
      NNSTP=NSTPLIM
!     
 11   IF(YN(7,NNSTP).LT.0.D0) THEN
         YN(7,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) &
              WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NNSTP)
      ENDIF

      RETURN
  END SUBROUTINE WRRKFT_ODE

!  --- Runge-Kutta method with tunneling of cutoff-resonant layer ---

  SUBROUTINE WRRKFT_WITHMC(Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,WRMODCONV,WRMODNWTN
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3),F(NEQ)
    REAL(rkind):: X0,XE,OMG,RL,RKRL,OXEFF,RHON,PHIL,ZL,Y7
    REAL(rkind):: RNPHI_IDEI,DELTA,RKPARA,RKPERP
    INTEGER:: INIT,NSTPLIM,NSTP,I,IOX,ID

      X0 = 0.D0
      XE = DELS
      INIT = 1
      NSTPLIM=INT(SMAX/DELS)
      NSTP=0
      YN(0,NSTP)=X0
      DO I=1,7
         YN(I,NSTP)=Y(I)
      ENDDO
      YN(8,NSTP)=0.D0
      OMG=2.D6*PI*RF
      IOX=0

      DO NSTP = 1,NSTPLIM
         Y7=Y(7)
         CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
         DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),OMG)

         IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,YK)
            DO I=1,3
               YM(I+3) = YK(I)
            END DO
         END IF 

!   --- Mode conversion

         RL  =SQRT(YM(1)**2+YM(2)**2)
         RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
         IF ( RKRL.GE.0.D0.AND.IOX.EQ.0 ) THEN
            CALL WRMODCONV(IOX,YM,F,OXEFF)
            IF(IOX.GE.100000) THEN
               WRITE(6,*) 'ERROR in WRMODCON_OX routine IOX=',IOX
            ELSE 
               DO I=1,NEQ
                  YM(I) = F(I)
               END DO
            END IF
         ENDIF

         YN(0,NSTP)=XE
         DO I=1,7
            YN(I,NSTP)=YM(I)
         ENDDO
         IF(YM(7).GT.0.D0) THEN
            YN(8,NSTP)=Y7-YM(7)
         ELSE
            YN(8,NSTP)=Y7
         ENDIF

         CALL PL_MAG_OLD(YM(1),YM(2),YM(3),RHON)
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
         RNPHI_IDEI= -YM(4)*VC/OMG*SIN(PHIL) + YM(5)*VC/OMG*COS(PHIL)
         DELTA=DISPXR( YM(1), YM(2), YM(3), YM(4), YM(5), YM(6), OMG )
         RKPARA=YM(4)*BNX+YM(5)*BNY+YM(6)*BNZ
         RKPERP=SQRT((YM(4)*YM(4)+YM(5)*YM(5)+YM(6)*YM(6))-RKPARA**2)

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
            IF(ID.EQ.1) &
                 WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NSTP)
         ENDIF

!IDEI	 WRITE(6,*) 'BX=',BNX*BABS,'BY=',BNY*BABS,'BZ=',BNZ*BABS
!	 WRITE(6,6001) YN(0,IT),RL,ZL,PHIL,YN(1,IT), &
!                      YN(2,IT),YN(3,IT),YN(6,IT)*VC/OMG,YN(7,IT), &
!      	               DELTA,RKPARA*VC/OMG,RKPERP*VC/OMG,RKRL, &
!      	               RNPHI_IDEI
! 6001    FORMAT(1H ,1P14E13.5)

         DO I=1,7
            Y(I)=YM(I)
         ENDDO
         X0=XE
         XE=X0+DELS
         IF(Y(7).LT.UUMIN) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
         CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
         IF(RHON.GT.RB/RA*RKAP) THEN
            NNSTP = NSTP
            GOTO 11
         ENDIF
      END DO
      NNSTP=NSTPLIM
!     
 11   IF(YN(7,NNSTP).LT.0.D0) THEN
         YN(7,NNSTP)=0.D0
      ENDIF
      IF(MDLWRW.GE.1) THEN
         IF(ID.EQ.0) &
              WRITE(6,'(1P7E11.3)') XE,RL,PHIL,ZL,RKRL,YM(7),YN(8,NNSTP)
      ENDIF

      RETURN
  END SUBROUTINE WRRKFT_WITHMC

! --- Symplectic method (not completed) ---

  SUBROUTINE WRSYMP(Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: F(NEQ)
    INTEGER:: NSTPLIM,NLPMAX,NSTP,I,ID,NLP,IERR
    REAL(rkind):: EPS,X,RL,PHIL,ZL,RKRL,OMG,RHON,OMG_C,WP2,ERROR,DELTA,Y7

    NSTPLIM=INT(SMAX/DELS)
    NLPMAX=10
    EPS=1.D-6

    NSTP=0
    X=0.D0
    YN(0,NSTP)=X
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0

    DO NSTP = 1,NSTPLIM
       Y7=Y(7)
       CALL SYMPLECTIC(Y,DELS,WRFDRVR,6,NLPMAX,EPS,NLP,ERROR,IERR)
       CALL WRFDRV(0.D0,Y,F)
       X=X+DELS

       YN(0,NSTP)=X
       DO I=1,6
          YN(I,NSTP)=Y(I)
       ENDDO
       Y(7)=Y(7)+F(7)*DELS
       IF(Y(7).GT.0.D0) THEN
          YN(7,NSTP)=Y(7)
          YN(8,NSTP)=-F(7)*DELS
       ELSE
          YN(7,NSTP)=0.D0
          YN(8,NSTP)=Y(7)
       ENDIF

       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          RL  =Y(1)
          PHIL=ASIN(Y(2)/(2.D0*PI*RR))
          ZL  =Y(3)
          RKRL=Y(4)
       ELSE
          RL  =SQRT(Y(1)**2+Y(2)**2)
          PHIL=ATAN2(Y(2),Y(1))
          ZL  =Y(3)
          RKRL=(Y(4)*Y(1)+Y(5)*Y(2))/RL
       ENDIF

       OMG=2.D6*PI*RF
       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       CALL PL_PROF_OLD(RHON)
       OMG_C = BABS*AEE/(AME)
       WP2=RN(1)*1.D20*AEE*AEE/(EPS0*AMP*PA(1))
       wp2=sqrt(WP2)

       DELTA=DISPXR(Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),OMG)

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
          IF(ID.EQ.1) &
               WRITE(6,'(1P7E11.3)') X,RL,PHIL,ZL,RKRL,Y(7),YN(8,NSTP)
       ENDIF

       IF(Y(7).LT.UUMIN) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       IF(RHON.GT.RB/RA*1.2D0) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
    ENDDO
    NNSTP=NSTPLIM

11  IF(YN(7,NNSTP).LT.0.D0) THEN
       YN(7,NNSTP)=0.D0
    ENDIF
    IF(MDLWRW.GE.1) THEN
       IF(ID.EQ.0) &
            WRITE(6,'(1P7E11.3)') X,RL,PHIL,ZL,RKRL,Y(7),YN(8,NNSTP)
    ENDIF

    RETURN
  END SUBROUTINE WRSYMP

!  --- slave routine for ray tracing ---

!  Y(1)=X
!  Y(2)=Y
!  Y(3)=Z
!  Y(4)=RKX
!  Y(5)=RKY
!  Y(6)=RKZ
!  Y(7)=W

  SUBROUTINE WRFDRV(X,Y,F)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,DISPXI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y(7)
    REAL(rkind),INTENT(OUT):: F(7)
    REAL(rkind):: VV,TT,XP,YP,ZP,RKXP,RKYP,RKZP,UU,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ,VDU
    REAL(rkind):: DUMMY

      VV=DELDER
      TT=DELDER
      DUMMY=X

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

!      WRITE(6,'(A,1P6E12.4)') 'XY:',X,Y(1),Y(4),Y(5),Y(6),Y(7)
!      WRITE(6,'(A,1P6E12.4)') 'F :',F(1),F(2),F(3),F(4),F(5),F(6)
      RETURN
  END SUBROUTINE WRFDRV

!  --- slave routine for symplectic method ---

  SUBROUTINE WRFDRVR(Y,F) 

    USE wrcomm
    USE wrsub,ONLY: DISPXR
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y(6)
    REAL(rkind),INTENT(OUT):: F(6)
    REAL(rkind):: VV,TT,XP,YP,ZP,RKXP,RKYP,RKZP,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS
    REAL(rkind):: VX,VY,VZ,VKX,VKY,VKZ

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
  END SUBROUTINE WRFDRVR

END MODULE wrexecr
