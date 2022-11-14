! wrexecr.f90

MODULE wrexecr

  PRIVATE
  PUBLIC wr_exec_rays
  PUBLIC wr_calc_pwr
  
CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE wr_exec_rays(ierr)

    USE wrcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL:: time1,time2
    REAL(rkind):: RK,PABSN
    INTEGER:: NRAY,nstp

    CALL GUTIME(TIME1)
    DO NRAY=1,NRAYMAX
       CALL wr_setup_start_point(NRAY,RAYS(0,0,NRAY),nstp,IERR)
       CALL wr_exec_single_ray(NRAY,RAYS(0,0,NRAY),nstp,IERR)
       IF(IERR.NE.0) cycle
       RK=SQRT(RAYS(4,NSTPMAX_NRAY(NRAY),NRAY)**2 &
              +RAYS(5,NSTPMAX_NRAY(NRAY),NRAY)**2 &
              +RAYS(6,NSTPMAX_NRAY(NRAY),NRAY)**2)
       PABSN=1.D0-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
       WRITE(6,'(A,I3,A,ES12.4,A,ES12.4)') &
            '    NRAY=',NRAY,'  RK=  ',RK,  '  PABS/PIN=', PABSN
    ENDDO

    CALL wr_calc_pwr

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

    RETURN
  END SUBROUTINE wr_exec_rays

  ! *** single ray tracing ***

  SUBROUTINE wr_exec_single_ray(NRAY,YN,nstp,IERR)
    USE wrcomm
    USE wrsub,ONLY: wrcale,wrcale_xyz,wr_cold_rkperp,wr_newton
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NRAY
    INTEGER,INTENT(INOUT):: nstp
    REAL(rkind),INTENT(IN):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: YA(NEQ)
    REAL(rkind):: RF,S,XP,YP,ZP,RKX,RKY,RKZ,UU,omega,rkv,rk

    IERR=0
    
    RF  =RFIN(NRAY)
    wr_nray_status%RF=RF
    omega=2.D6*PI*RF
    rkv=omega/VC

    S   =YN(0,nstp)
    XP  =YN(1,nstp)
    YP  =YN(2,nstp)
    ZP  =YN(3,nstp)
    RKX =YN(4,nstp)
    RKY =YN(5,nstp)
    RKZ =YN(6,nstp)
    UU  =YN(7,nstp)

    IF(idebug_wr(8).NE.0) THEN
       rk=SQRT(rkx**2+rky**2+rkz**2)
       WRITE(6,'(A,A,I4,I8)') '*** idebug_wr(8): wr_exec_single_ray: ', &
            'nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp       =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz    =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnkx,rnky,rnkz =',RKX/rkv,RKY/rkv,RKZ/rkv
       WRITE(6,'(A,2ES12.4)') '   UU,S,rk,rn     =',UU,S,rk,rk/rkv
    END IF

    YA(1)= XP
    YA(2)= YP
    YA(3)= ZP
    YA(4)= RKX
    YA(5)= RKY
    YA(6)= RKZ
    YA(7)= UU
    
    IF(MDLWRQ.EQ.0) THEN
       CALL WRRKFT(nstp,YA,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.1) THEN
       CALL WRRKFT_WITHD0(nstp,YA,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.2) THEN
       CALL WRRKFT_WITHMC(nstp,YA,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.3) THEN
       CALL WRRKFT_RKF(nstp,YA,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.4) THEN
       CALL WRSYMP(nstp,YA,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSEIF(MDLWRQ.EQ.5) THEN
       CALL WRRKFT_ODE(nstp,YA,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY))
    ELSE
       WRITE(6,*) 'XX WRCALC: unknown MDLWRQ =', MDLWRQ
       IERR=1
       RETURN
    ENDIF
    DO NSTP=0,NSTPMAX_NRAY(NRAY)
       RAYRB1(NSTP,NRAY)=0.D0
       RAYRB2(NSTP,NRAY)=0.D0
    END DO

    CALL WRCALE(RF,RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY),NRAY)

    RETURN
  END SUBROUTINE

!  --- original Runge-Kutta method ---

  SUBROUTINE WRRKFT(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    USE librk
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_in
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I
    REAL(rkind):: X0,XE,RHON,PW,R

    X0 = 0.D0
    XE = DELS
    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_in
    YN(0,NSTP)=X0
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_in+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS

       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.RR+1.1D0*RB.OR. &
          R.LT.RR-1.1D0*RB.OR. &
          Y(3).GT. RKAP*1.1D0*RB.OR. &
          Y(3).LT.-RKAP*1.1D0*RB))) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             NNSTP = NSTP
             GOTO 11
          ENDIF
       ENDIF
    END DO
    NNSTP=NSTPLIM

11  IF(YN(7,NNSTP).LT.0.D0) YN(7,NNSTP)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT

!  --- original Runge-Kutta method with correction for D=0 ---

  SUBROUTINE WRRKFT_WITHD0(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,WRMODNWTN
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    USE librk
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_in
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I
    REAL(rkind):: X0,XE,PW,DELTA,RHON,R,RF,omega

    RF=wr_nray_status%RF
    omega=2.D6*PI*RF
    nstp=nstp_in

    X0 = 0.D0
    XE = DELS

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(11): wrrkft_withd0: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,xe,smax,dels =',X0,XE,SMAX,DELS
    END IF
       
    NSTPLIM=MIN(INT(SMAX/DELS),NSTPMAX)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    YN(0,NSTP_in)=X0
    DO I=1,7
       YN(I,NSTP_in)=Y(I)
    ENDDO
    YN(8,NSTP_in)=0.D0

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(11): wrrkft_withd0: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,y1,y2,y3 =',X0,Y(1),Y(2),Y(3)
       WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',Y(4),Y(5),Y(6),Y(7)
    END IF

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_in+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)

       IF(idebug_wr(11).NE.0) THEN
          WRITE(6,'(A,I8)') '*** idebug_wr(11): wrrkft_withd0: nstp=',nstp
          WRITE(6,'(A,4ES12.4)') '      x0,y1,y2,y3 =',X0,Y(1),Y(2),Y(3)
          WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',Y(4),Y(5),Y(6),Y(7)
          WRITE(6,'(A,4ES12.4)') '      xe,y1,y2,y3 =',XE,YM(1),YM(2),YM(3)
          WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',YM(4),YM(5),YM(6),YM(7)
       END IF

       DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),omega)
       IF (ABS(DELTA).GT.1.0D-6) THEN
          CALL WRMODNWTN(YM,omega,YK)
          DO I=1,3
             YM(I+3) = YK(I)
          END DO
       END IF

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS
       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.RR+1.1D0*RB.OR. &
          R.LT.RR-1.1D0*RB.OR. &
          Y(3).GT. RKAP*1.1D0*RB.OR. &
          Y(3).LT.-RKAP*1.1D0*RB))) THEN
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
 
11  IF(YN(7,NNSTP).LT.0.D0) YN(7,NNSTP)=0.D0
    CALL wr_write_line(NNSTP,X0,YM,YN(8,NNSTP))

    RETURN
  END SUBROUTINE WRRKFT_WITHD0

!  --- Runge-Kutta method using ODE library ---

  SUBROUTINE WRRKFT_ODE(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    USE librk
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_in
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I
    REAL(rkind):: X0,XE,PW,RHON,R

    X0 = 0.D0
    XE = DELS     
    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)
    
    NSTP=nstp_in
    YN(0,NSTP)=X0
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_in+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
       
       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS

       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.RR+1.1D0*RB.OR. &
          R.LT.RR-1.1D0*RB.OR. &
          Y(3).GT. RKAP*1.1D0*RB.OR. &
          Y(3).LT.-RKAP*1.1D0*RB))) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             NNSTP = NSTP
             GOTO 11
          ENDIF
       END IF
    END DO
    NNSTP=NSTPLIM

11  IF(YN(7,NNSTP).LT.0.D0) YN(7,NNSTP)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_ODE

!  --- Runge-Kutta method with tunneling of cutoff-resonant layer ---

  SUBROUTINE WRRKFT_WITHMC(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,WRMODCONV,WRMODNWTN
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    USE librk
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_in
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3),F(NEQ)
    REAL(rkind):: X0,XE,RF,omega,OXEFF,RHON,PW,RL,RKRL,DELTA,R
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I,IOX

    RF=wr_nray_status%RF
    omega=2.D6*PI*RF
    
    X0 = 0.D0
    XE = DELS
    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)
    
    NSTP=nstp_in
    YN(0,NSTP)=X0
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0
    IOX=0

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_in+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,WRFDRV,X0,XE,1,Y,YM,WORK)
       
       DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),omega)

       IF (ABS(DELTA).GT.1.0D-6) THEN
          CALL WRMODNWTN(YM,omega,YK)
          DO I=1,3
             YM(I+3) = YK(I)
          END DO
       END IF

!   --- Mode conversion

       RL  =SQRT(YM(1)**2+YM(2)**2)
       RKRL=(YM(4)*YM(1)+YM(5)*YM(2))/RL
       IF ( RKRL.GE.0.D0.AND.IOX.EQ.0 ) THEN
          CALL WRMODCONV(IOX,YM,F,OXEFF,omega)
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
       YN(8,NSTP)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS
       
       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.RR+1.1D0*RB.OR. &
          R.LT.RR-1.1D0*RB.OR. &
          Y(3).GT. RKAP*1.1D0*RB.OR. &
          Y(3).LT.-RKAP*1.1D0*RB))) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       IF(RHON.GT.RB/RA) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
    END DO
    NNSTP=NSTPLIM

11  IF(YN(7,NNSTP).LT.0.D0) YN(7,NNSTP)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_WITHMC

! --- Auto-step-size Runge-Kutta-F method ---

  SUBROUTINE WRRKFT_RKF(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE librkf
    USE plprof,ONLY: pl_mag_old
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_in
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,INIT,NDE,IER,I
    REAL(rkind):: RELERR,ABSERR,X0,XE,WORK0,PW,YM(NEQ),RHON,R
    REAL(rkind):: ESTERR(NEQ),WORK1(NEQ),WORK2(NEQ),WORK3(NEQ),WORK4(NEQ,11)

    RELERR = EPSRAY
    ABSERR = EPSRAY
    INIT = 1

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(11): wrrkft_rkf: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      X0,xe,smax,dels =',X0,XE,SMAX,DELS
    END IF
       
    X0 = 0.D0
    XE = DELS     
    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_in
    YN(0,NSTP)=X0
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(12): wrrkft_rkf: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,y1,y2,y3 =',X0,Y(1),Y(2),Y(3)
       WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',Y(4),Y(5),Y(6),Y(7)
    END IF

       CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_in+1,NSTPLIM
       PW=Y(7)
       CALL RKF(7,WRFDRV,X0,XE,Y,INIT,RELERR,ABSERR,YM, &
                ESTERR,NDE,IER,WORK0,WORK1,WORK2,WORK3,WORK4)
       IF (IER .NE. 0) THEN
          WRITE(6,'(A,2I6)') 'XX wrrkft_rkf: NSTP,IER=',NSTP,IER
          RETURN
       ENDIF

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS
       
       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.RR+1.1D0*RB.OR. &
          R.LT.RR-1.1D0*RB.OR. &
          Y(3).GT. RKAP*1.1D0*RB.OR. &
          Y(3).LT.-RKAP*1.1D0*RB))) THEN
          NNSTP = NSTP
          GOTO 8000
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL pl_mag_old(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             NNSTP = NSTP
             GOTO 8000
          END IF
       ENDIF
    END DO
    NNSTP=NSTPLIM
     
8000 CONTINUE
    IF(YN(7,NNSTP).LT.0.D0) YN(7,NNSTP)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_RKF

! --- Symplectic method (not completed) ---

  SUBROUTINE WRSYMP(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    USE libsympl
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_in
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: NNSTP
    REAL(rkind):: X,F(NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NLPMAX,NSTP,I,NLP,IERR
    REAL(rkind):: EPS,PW,ERROR,RHON,R

    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)
    NLPMAX=10
    EPS=1.D-6

    NSTP=nstp_in
    X=0.D0
    YN(0,NSTP)=X
    DO I=1,7
       YN(I,NSTP)=Y(I)
    ENDDO
    YN(8,NSTP)=0.D0

    CALL wr_write_line(NSTP,X,Y,YN(8,NSTP))

    DO NSTP = nstp_in+1,NSTPLIM
       PW=Y(7)
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

       CALL wr_write_line(NSTP,X,Y,YN(8,NSTP))

       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.RR+1.1D0*RB.OR. &
          R.LT.RR-1.1D0*RB.OR. &
          Y(3).GT. RKAP*1.1D0*RB.OR. &
          Y(3).LT.-RKAP*1.1D0*RB))) THEN
          NNSTP = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             NNSTP = NSTP
             GOTO 11
          ENDIF
       END IF
    END DO
    NNSTP=NSTPLIM

11  IF(YN(7,NNSTP).LT.0.D0) YN(7,NNSTP)=0.D0
    CALL wr_write_line(NSTP,X,Y,YN(8,NSTP))
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

      IF(idebug_wr(12).NE.0) THEN
         WRITE(6,'(A)') '*** idebug_wr(12): wrfdrv'
         WRITE(6,'(A,3ES12.4)') 'x7:',X,Y(7),F(7)
         WRITE(6,'(A,6ES12.4)') 'y :',Y(1),Y(2),Y(3),Y(4),Y(5),Y(6)
         WRITE(6,'(A,6ES12.4)') 'f :',F(1),F(2),F(3),F(4),F(5),F(6)
      END IF
      RETURN
  END SUBROUTINE WRFDRV

!  --- slave routine for symplectic method ---

  SUBROUTINE WRFDRVR(Y,F) 

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
  END SUBROUTINE WRFDRVR

  ! --- write line ---
  
  SUBROUTINE wr_write_line(NSTP,X,Y,PABS)
    USE wrcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NSTP
    REAL(rkind),INTENT(IN):: X,Y(NEQ),PABS
    REAL(rkind):: RL,PHIL,ZL,RKRL
    INTEGER:: ID
    INTEGER,SAVE:: NSTP_SAVE=-1

    IF(MDLWRW.EQ.0) RETURN

    IF(NSTP.EQ.0) THEN
       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          WRITE(6,'(A,A)') &
               '      S          X        ANG          Z    ', &
               '     RX          W       PABS'
       ELSE IF(MODELG.EQ.11) THEN
          WRITE(6,'(A,A)') &
               '      S          X          Y          Z    ', &
               '     RX          W       PABS'
       ELSE
          WRITE(6,'(A,A)') &
               '      S          R        PHI          Z    ', &
               '    RKR          W       PABS'
       ENDIF
    END IF

    IF(NSTP.EQ.0) THEN
       ID=1
    ELSE
       ID=0
       SELECT CASE(MDLWRW)
       CASE(1)
          ID=1
       CASE(2)
          IF(MOD(NSTP,10).EQ.0) ID=1
       CASE(3)
          IF(MOD(NSTP,100).EQ.0) ID=1
       CASE(4)
          IF(MOD(NSTP,1000).EQ.0) ID=1
       CASE(5)
          IF(MOD(NSTP,10000).EQ.0) ID=1
       END SELECT
       IF(NSTP.EQ.NSTP_SAVE) ID=0
    END IF
    
    IF(ID.EQ.1) THEN
       IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
          RL  =Y(1)
          PHIL=ASIN(Y(2)/(2.D0*PI*RR))
          ZL  =Y(3)
          RKRL=Y(4)
       ELSE IF(MODELG.EQ.11) THEN
          RL  =Y(1)
          PHIL=Y(2)
          ZL  =Y(3)
          RKRL=Y(4)
       ELSE
          RL  =SQRT(Y(1)**2+Y(2)**2)
          PHIL=ATAN2(Y(2),Y(1))
          ZL  =Y(3)
          RKRL=(Y(4)*Y(1)+Y(5)*Y(2))/RL
       ENDIF
       WRITE(6,'(1P7E11.3)') X,RL,PHIL,ZL,RKRL,Y(7),PABS
       NSTP_SAVE=NSTP
    END IF
  END SUBROUTINE wr_write_line

  ! --- calculation of radial profile of power absorption ---

  SUBROUTINE wr_calc_pwr

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: pl_mag_old,pl_rzsu
    USE libgrf
    IMPLICIT NONE
    INTEGER,PARAMETER:: nsum=201
    REAL(rkind),ALLOCATABLE:: rsu(:),zsu(:)
    INTEGER:: nrs,nray,nstp,nrs1,nrs2,ndrs,locmax
    INTEGER:: nrl1,nrl2,ndrl,nsumax,nsu,nrl
    REAL(rkind):: drs,xl,yl,zl,rs1,rs2,sdrs,delpwr,pwrmax,dpwr,ddpwr
    REAL(rkind):: rlmin,rlmax,drl,rl1,rl2,sdrl
    INTEGER,SAVE:: nrsmax_save=0,nrlmax_save=0,nraymax_save=0,nstpmax_save=0

!   ----- evaluate plasma major radius range -----

    ALLOCATE(rsu(nsum),zsu(nsum))
    CALL pl_rzsu(rsu,zsu,nsum,nsumax)
    rlmin=rsu(1)
    rlmax=rsu(1)
    DO nsu=2,nsumax
       rlmin=MIN(rlmin,rsu(nsu))
       rlmax=MAX(rlmax,rsu(nsu))
    END DO
    DEALLOCATE(rsu,zsu)

    !   --- allocate variables for power deposition profile ---

    IF(nstpmax.NE.nstpmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(rs_nstp_nray)) DEALLOCATE(rs_nstp_nray)
       IF(ALLOCATED(rl_nstp_nray)) DEALLOCATE(rl_nstp_nray)
       ALLOCATE(rs_nstp_nray(nstpmax,nraymax),rl_nstp_nray(nstpmax,nraymax))
    END IF
    IF(nrsmax.NE.nrsmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_nrs)) DEALLOCATE(pos_nrs)
       IF(ALLOCATED(pwr_nrs)) DEALLOCATE(pwr_nrs)
       IF(ALLOCATED(pwr_nrs_nray)) DEALLOCATE(pwr_nrs_nray)
       ALLOCATE(pos_nrs(nrsmax),pwr_nrs(nrsmax),pwr_nrs_nray(nrsmax,nraymax))
    END IF
    IF(nrlmax.NE.nrlmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_nrl)) DEALLOCATE(pos_nrl)
       IF(ALLOCATED(pwr_nrl)) DEALLOCATE(pwr_nrl)
       IF(ALLOCATED(pwr_nrl_nray)) DEALLOCATE(pwr_nrl_nray)
       ALLOCATE(pos_nrl(nrlmax),pwr_nrl(nrlmax),pwr_nrl_nray(nrlmax,nraymax))
    END IF
    IF(nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_pwrmax_rs_nray)) DEALLOCATE(pos_pwrmax_rs_nray)
       IF(ALLOCATED(pwrmax_rs_nray)) DEALLOCATE(pwrmax_rs_nray)
       IF(ALLOCATED(pos_pwrmax_rl_nray)) DEALLOCATE(pos_pwrmax_rl_nray)
       IF(ALLOCATED(pwrmax_rl_nray)) DEALLOCATE(pwrmax_rl_nray)
       ALLOCATE(pos_pwrmax_rs_nray(nraymax),pwrmax_rs_nray(nraymax))
       ALLOCATE(pos_pwrmax_rl_nray(nraymax),pwrmax_rl_nray(nraymax))
    END IF
    nrsmax_save=nrsmax
    nrlmax_save=nrlmax
    nraymax_save=nraymax
    nstpmax_save=nstpmax

!     ----- Setup for RADIAL DEPOSITION PROFILE (Minor radius) -----

    drs=1.D0/nrsmax
    DO nrs=1,nrsmax
       pos_nrs(nrs)=(DBLE(nrs)-0.5D0)*drs
    ENDDO
    DO nray=1,nraymax
       DO nrs=1,nrsmax
          pwr_nrs_nray(nrs,nray)=0.D0
       ENDDO
    ENDDO
    DO nrs=1,nrsmax
       pwr_nrs(nrs)=0.D0
    ENDDO

!     ----- Setup for RADIAL DEPOSITION PROFILE (Major radius) -----

    drl=(rlmax-rlmin)/(nrlmax-1)
    DO nrl=1,nrlmax
       pos_nrl(nrl)=rlmin+(nrl-1)*drl
    ENDDO
    DO nray=1,nraymax
       DO nrl=1,nrlmax
          pwr_nrl_nray(nrl,nray)=0.D0
       ENDDO
    ENDDO
    DO nrl=1,nrlmax
       pwr_nrl(nrl)=0.D0
    ENDDO

!   --- calculate power deposition density ----

    DO nray=1,nraymax
       DO nstp=0,nstpmax_nray(nray)-1

          ! --- minor radius porfile ---
          
          xl=rays(1,nstp,nray)
          yl=rays(2,nstp,nray)
          zl=rays(3,nstp,nray)
          CALL pl_mag_old(xl,yl,zl,rs1)
          xl=rays(1,nstp+1,nray)
          yl=rays(2,nstp+1,nray)
          zl=rays(3,nstp+1,nray)
          CALL pl_mag_old(xl,yl,zl,rs2)
          IF(rs1.LE.1.D0.OR.rs1.LE.1.D0) THEN
             nrs1=INT(rs1/drs)+1
             nrs2=INT(rs2/drs)+1
             IF(nrs1.GT.nrsmax) THEN
                nrs1=nrsmax
                IF(nrs2.GT.nrsmax) EXIT
             ENDIF
             IF(nrs2.GT.nrsmax) nrs2=nrsmax
                   
             ndrs=ABS(nrs2-nrs1)
             IF(ndrs.EQ.0) THEN
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +rays(8,nstp+1,nray)
             ELSE IF(nrs1.lt.nrs2) THEN
                sdrs=(rs2-rs1)/drs
                delpwr=rays(8,nstp+1,nray)/sdrs
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +(DBLE(nrs1)-rs1/drs)*delpwr
                DO nrs=nrs1+1,nrs2-1
                   pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray)+delpwr
                ENDDO
                pwr_nrs_nray(nrs2,nray)=pwr_nrs_nray(nrs2,nray) &
                     +(rs2/drs-DBLE(nrs2-1))*delpwr
             ELSE
                sdrs=(rs1-rs2)/drs
                delpwr=rays(8,nstp+1,nray)/sdrs
                pwr_nrs_nray(nrs2,nray)=pwr_nrs_nray(nrs2,nray) &
                     +(DBLE(nrs2)-rs2/drs)*delpwr
                DO nrs=nrs2+1,nrs1-1
                   pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray)+delpwr
                ENDDO
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +(rs1/drs-DBLE(nrs1-1))*delpwr
             END IF
          ENDIF

          ! --- major radius profile ---

          xl=rays(1,nstp,nray)
          yl=rays(2,nstp,nray)
          rl1=SQRT(xl**2+yl**2)
          nrl1=INT((rl1-rlmin)/drl)+1
          xl=rays(1,nstp+1,nray)
          yl=rays(2,nstp+1,nray)
          rl2=SQRT(xl**2+yl**2)
          nrl2=INT((rl2-rlmin)/drl)+1

          IF((nrl1.GE.1.and.nrl1.LE.nrlmax).OR. &
             (nrl2.GE.1.and.nrl2.LE.nrlmax)) THEN
             IF(nrl1.LT.1) nrl1=1
             IF(nrl1.GT.nrlmax) nrl1=nrlmax
             IF(nrl2.LT.1) nrl2=1
             IF(nrl2.GT.nrlmax) nrl2=nrlmax
             ndrl=ABS(nrl2-nrl1)
             IF(ndrl.EQ.0) THEN
                pwr_nrl_nray(nrl1,nray)=pwr_nrl_nray(nrl1,nray) &
                     +rays(8,nstp+1,nray)
             ELSE IF(nrl1.LT.nrl2) THEN
                sdrl=(rl2-rl1)/drl
                delpwr=rays(8,nstp+1,nray)/sdrl
                pwr_nrl_nray(nrl1,nray)=pwr_nrl_nray(nrl1,nray) &
                     +(DBLE(nrl1)-(rl1-rlmin)/drl)*delpwr
                DO nrl=nrl1+1,nrl2-1
                   pwr_nrl_nray(nrl,nray) &
                        =pwr_nrl_nray(nrl,nray)+delpwr
                END DO
                pwr_nrl_nray(nrl2,nray)=pwr_nrl_nray(nrl2,nray) &
                     +((rl2-rlmin)/drl-dble(nrl2-1))*delpwr
             ELSE
                sdrl=(rl1-rl2)/drl
                delpwr=rays(8,nstp+1,nray)/sdrl
                pwr_nrl_nray(nrl2,nray)=pwr_nrl_nray(nrl2,nray) &
                     +(DBLE(nrl2)-(rl2-rlmin)/drl)*delpwr
                DO nrl=nrl2+1,nrl1-1
                   pwr_nrl_nray(nrl,nray) &
                        =pwr_nrl_nray(nrl,nray)+delpwr
                END DO
                pwr_nrl_nray(nrl1,nray)=pwr_nrl_nray(nrl1,nray) &
                     +((rl1-rlmin)/drl-DBLE(nrl1-1))*delpwr
             END IF
          END IF
       ENDDO

       ! --- power divided by division area ---

       DO nrs=1,nrsmax
          pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray) &
               /(2.D0*PI*(DBLE(nrs)-0.5D0)*drs*drs)
       ENDDO
       DO nrl=1,nrlmax
          pwr_nrl_nray(nrl,nray)=pwr_nrl_nray(nrl,nray) &
               /(2.D0*PI*(DBLE(nrl)-0.5D0)*drl*drl)
       ENDDO
    ENDDO

!     ----- sum of power for each ray-------

    DO nrs=1,nrsmax
       pwr_nrs(nrs)=0.D0
       DO nray=1,nraymax
          pwr_nrs(nrs)=pwr_nrs(nrs)+pwr_nrs_nray(nrs,nray)
       ENDDO
    ENDDO

    DO nrl=1,nrlmax
       pwr_nrl(nrl)=0.D0
       DO nray=1,nraymax
          pwr_nrl(nrl)=pwr_nrl(nrl)+pwr_nrl_nray(nrl,nray)
       ENDDO
    ENDDO

!     ----- find location of absorbed power peak -----

    DO nray=1,nraymax
       pwrmax=0.D0
       locmax=0
       DO nrs=1,nrsmax
          IF(pwr_nrs_nray(nrs,nray).GT.pwrmax) THEN
             pwrmax=pwr_nrs_nray(nrs,nray)
             locmax=nrs
          ENDIF
       END DO
       IF(locmax.LE.1) THEN
          pos_pwrmax_rs_nray(nray)=pos_nrs(1)
       ELSE IF(locmax.GE.nrsmax) THEN
          pos_pwrmax_rs_nray(nray)=pos_nrs(nrsmax)
       ELSE
          dpwr =(pwr_nrs_nray(locmax+1,nray) &
                -pwr_nrs_nray(locmax-1,nray))/(2.D0*drs)
          ddpwr=(pwr_nrs_nray(locmax+1,nray) &
              -2*pwr_nrs_nray(locmax  ,nray) &
                +pwr_nrs_nray(locmax-1,nray))/drs**2
          pos_pwrmax_rs_nray(nray)=(locmax-0.5D0)/(nrsmax-1.D0)-dpwr/ddpwr
          pwrmax_rs_nray(nray)=pwrmax-dpwr**2/(2.D0*ddpwr)
       ENDIF
    END DO
   
    pwrmax=0.D0
    locmax=0
    DO nrl=1,nrlmax
       IF(pwr_nrl(nrl).GT.pwrmax) THEN
          pwrmax=pwr_nrl(nrl)
          locmax=nrl
       ENDIF
    END DO
    IF(locmax.LE.1) THEN
       pos_pwrmax_rl=pos_nrl(1)
    ELSE IF(locmax.GE.nrlmax) THEN
       pos_pwrmax_rl=pos_nrl(nrlmax)
    ELSE
       dpwr =(pwr_nrl(locmax+1) &
             -pwr_nrl(locmax-1))/(2.D0*drl)
       ddpwr=(pwr_nrl(locmax+1) &
           -2*pwr_nrl(locmax  ) &
             +pwr_nrl(locmax-1))/drl**2
       pos_pwrmax_rl=(locmax-0.5D0)/(nrlmax-1.D0)-dpwr/ddpwr
       pwrmax_rl=pwrmax-dpwr**2/(2.D0*ddpwr)
    ENDIF
   
    WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
         '    PWRMAX=',pwrmax_rl,'  AT RL =',pos_pwrmax_rl

!    CALL PAGES
!    CALL grd1d(1,pos_nrs,pwr_nrs_nray,nrsmax,nrsmax,nraymax, &
!         '@pwr-nrs vs. pos-nrs@')
!    CALL grd1d(2,pos_nrl,pwr_nrl_nray,nrlmax,nrlmax,nraymax, &
!         '@pwr-nrl vs. pos-nrl@')
!    CALL PAGEE

    RETURN
  END SUBROUTINE wr_calc_pwr

END MODULE wrexecr
