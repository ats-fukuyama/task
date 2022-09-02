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
    INTEGER:: NRAY

    CALL GUTIME(TIME1)
    DO NRAY=1,NRAYMAX
       CALL wr_exec_single_ray(NRAY,IERR)
       IF(IERR.NE.0) cycle

       WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
            '    RKRI=  ',RKRI,  '  PABS/PIN=', &
            1.D0-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
    ENDDO

    CALL wr_calc_pwr

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

    RETURN
  END SUBROUTINE wr_exec_rays

!   ***** Ray tracing module *****

  SUBROUTINE wr_exec_single_ray(NRAY,IERR)

    USE wrcomm
    USE wrsub,ONLY: wrcale,wrcale_xyz,wr_cold_rkperp,wrnwtn
    USE plprof,ONLY: pl_mag_type,pl_mag,pl_prf_type,pl_prof
    USE plprofw,ONLY: pl_prfw_type,pl_profw
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NRAY
    INTEGER,INTENT(OUT):: IERR
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_prf_type),DIMENSION(nsmax):: plf
    REAL(rkind):: YA(NEQ)
    REAL(rkind):: YN(0:NEQ,0:NSTPMAX)
    REAL(rkind):: ANGPH,ANGZ,EPARA,EPERP
    REAL(rkind):: R,Z,PHI
    REAL(rkind):: rhon,rk,rkpara,rkperp_1,rkperp_2,rkphi,rkr
    REAL(rkind):: rkr_1,rkr_2,rkrz_old,rkrz_new_1,rkrz_new_2
    REAL(rkind):: rkx,rky,rkz,rkz_1,rkz_2,rnr,rnz,rnphi,rnr2
    REAL(rkind):: x,y,s,sinp2,sint2
    REAL(rkind):: W(4)
    INTEGER:: mode,nstp,id,i

    IERR=0

    ! --- setup common values and initial values ---

    RF=RFIN(NRAY)
    RPI=RPIN(NRAY)
    ZPI=ZPIN(NRAY)
    PHII=PHIIN(NRAY)
    RKR0=RKRIN(NRAY)
    RNZI=RNZIN(NRAY)
    RNPHII=RNPHIIN(NRAY)
    ANGZ=ANGZIN(NRAY)
    ANGPH=ANGPHIN(NRAY)
    MODEW=MODEWIN(NRAY)
    UUI=UUIN(NRAY)

    IF(MDLWRI.NE.100)THEN
       SINP2=SIN(ANGZ*PI/180.D0)**2
       SINT2=SIN(ANGPH*PI/180.D0)**2
       RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
       RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
       IF(ANGZ.LT.0.D0) RNZI=-RNZI
       IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
    ENDIF
    ! --- initialize variable array RAYIN ---

    RAYIN(1,NRAY)=RF
    RAYIN(2,NRAY)=RPI
    RAYIN(3,NRAY)=ZPI
    RAYIN(4,NRAY)=PHII
    RAYIN(5,NRAY)=RKR0
    RAYIN(6,NRAY)=ANGZ
    RAYIN(7,NRAY)=ANGPH
    RAYIN(8,NRAY)=UUI

    RKR  = RKR0
    RKZ  = 2.D6*PI*RF*RNZI  /VC
    RKPHI= 2.D6*PI*RF*RNPHII/VC
    
    R=RPI
    Z=ZPI
    PHI=PHII
    
    mode=0
    nstp=0
    s=0.D0

    ! --- mode=0: vacuum region ---

    DO WHILE(mode.EQ.0)
       X=R*COS(PHI)
       Y=R*SIN(PHI)
       CALL pl_mag(X,Y,Z,mag)
       rhon=mag%rhon
       CALL pl_prof(rhon,plf)
       CALL pl_profw(rhon,plfw)

!       WRITE(6,'(A,3ES12.4)') 'pne:',plfw(1)%rn,pne_threshold,rhon

!       WRITE(6,'(A)') '## wr_exec_single_ray: '
!       WRITE(6,'(A,3ES12.4)') '   R,Z,PHI      =',R,Z,PHI
!       WRITE(6,'(A,3ES12.4)') '   X,Y,Z        =',X,Y,Z
!       WRITE(6,'(A,3ES12.4)') '   RKX,RKY,RKZ  =',RKX,RKY,RKZ

       WRITE(6,'(A,2ES12.4)') 'rn,pne_th:',plfw(1)%rn,pne_threshold
       IF(plfw(1)%rn.LE.pne_threshold) THEN

          ! --- if in vacuume, advance by ds ---
          
          RNR2=1.D0-RNZ**2-RNPHI**2
          IF(RNR2.LT.0.D0) THEN
             mode=4
             WRITE(6,'(A)') 'XX wr_exec_single_ray: RNRI2 is negative'
             WRITE(6,'(A,3ES12.4)') '   RNZ,RNPHI,RNR2=',RNZ,RNPHI,RNR2
             ierr=401
             RETURN
          END IF
          RNR=SIGN(SQRT(RNR2),RKR)
          RKR=2.D6*PI*RF*RNR/VC
          RKX=RKR*COS(PHI)-RKPHI*SIN(PHI)
          RKY=RKR*SIN(PHI)+RKPHI*COS(PHI)
          RKZ=RKZI
          RK=SQRT(RKX**2+RKY**2+RKZ**2)
          X=X+DELS*RKX/RK
          Y=Y+DELS*RKY/RK
          Z=Z+DELS*RKZ/RK
          R=SQRT(X**2+Y**2)
          PHI=ATAN2(Y,X)
          s=s+dels
          nstp=nstp+1
          YN(0,nstp)= s
          YN(1,nstp)= R
          YN(2,nstp)= PHI
          YN(3,nstp)= Z
          YN(4,nstp)= RKR
          YN(5,nstp)= RKPHI
          YN(6,nstp)= RKZ
          YN(7,nstp)= UUI
          nstp=nstp+1

!          WRITE(6,'(A,I4,6ES12.4)') 'RK:',nstp,R,PHI,Z,RKR,RKPHI,RKZ
!          WRITE(6,'(A,2ES12.4)') 'R: ',R,RMAX_WR
!          WRITE(6,'(A,2ES12.4)') 'R: ',R,RMIN_WR
!          WRITE(6,'(A,2ES12.4)') 'Z: ',Z,ZMAX_WR
!          WRITE(6,'(A,2ES12.4)') 'Z: ',Z,ZMIN_WR
!          WRITE(6,'(A,2ES12.4)') 'S: ',S,SMAX

          IF(R.GT.RMAX_WR.OR. &
             R.LT.RMIN_WR.OR. &
             Z.GT.ZMAX_WR.OR. &
             Z.LT.ZMIN_WR.OR. &
             S.GT.SMAX) THEN

             ! --- If ray goes out of calculation region, mode=2 ---

             mode=2
             WRITE(6,'(A,2I6,2ES12.4)') &
                  'wr_exec_ray_single: nray,nstp,R,Z=',NRAY,nstp,R,Z
             ierr=102
             RETURN
          END IF

       ELSE

          ! --- If ray is in plasma, mode=1 ---
          
          mode=1
       END IF
    END DO

    ! --- calculate wave number vector of two EM modes ---
    
    ! --- solve cold dispersion for given k_para ---
       
    RKX=RKR*COS(PHI)-RKPHI*SIN(PHI)
    RKY=RKR*SIN(PHI)+RKPHI*COS(PHI)
    RKZ=RKZ
    RKPARA=RKX*mag%bnx+RKY*mag%bny+RKZ*mag%bnz
    CALL WR_COLD_RKPERP(R,Z,PHI,RKPARA,RKPERP_1,RKPERP_2)
    WRITE(6,'(A,1P3E12.4)') &
         'RKRPARA,RKPERP_1,RKPERP_2=',RKPARA,RKPERP_1,RKPERP_2
       
    ! --- new k_r_new and k_z_new is determined from k_perp
    !         and the direction of old k_r and k_z
    !             k_rz_old^2=k_r_old^2+k_z_old^2
    !             k_rz_new^2=k_perp^2+k_para^2-k_phi\2
    !             k_r_new=(k_rz_new/k_rz_old) k_r_oldf
    !             k_z_new=(k_rz_new/k_rz_old) k_z_oldf
       
    RKRZ_old=SQRT(RKR**2+RKZ**2)
    RKRZ_new_1=RKPERP_1**2-RKPARA**2-RKPHI**2
    RKRZ_new_2=RKPERP_2**2-RKPARA**2-RKPHI**2
    IF(RKRZ_new_1.GE.0.D0) THEN
       RKRZ_new_1=SQRT(RKRZ_new_1)
       RKR_1=(RKRZ_new_1/RKRZ_old)*RKR
       RKZ_1=(RKRZ_new_1/RKRZ_old)*RKZ
       CALL WRCALE_XYZ(X,Y,Z,RKR_1,RKZ_1,RKPHI,EPARA,EPERP)
       WRITE(6,'(A,1P4E12.4)') &
            'RKR_1,RKZ_1,EPARA,EPERP=',RKR_1,RKZ_1,EPARA,EPERP
    ELSE
       RKR_1=RKR
       RKZ_1=RKZ
    END IF
    
    IF(RKRZ_new_2.GE.0.D0) THEN
       RKRZ_new_2=SQRT(RKRZ_new_2)
       RKR_2=(RKRZ_new_2/RKRZ_old)*RKR
       RKZ_2=(RKRZ_new_2/RKRZ_old)*RKZ
       CALL WRCALE_XYZ(X,Y,Z,RKR_2,RKZ_2,RKPHI,EPARA,EPERP)
       WRITE(6,'(A,1P4E12.4)') &
            'RKR_2,RKZ_2,EPARA,EPERP=',RKR_2,RKZ_2,EPARA,EPERP
    ELSE
       RKR_2=RKR
       RKZ_2=RKZ
    END IF

    SELECT CASE(MODEW)
    CASE(0)
       W(1)=-ABS( RKR_1-RKR0)
       W(2)=-ABS( RKR_2-RKR0)
       W(3)=-ABS(-RKR_1-RKR0)
       W(4)=-ABS(-RKR_2-RKR0)
       ID=MAXLOC(W,1)
       WRITE(6,'(A,5ES12.4)')    'RKR0:',RKR0,RKR_1,RKR_2,-RKR_1,-RKR_2
       WRITE(6,'(A,I6,4ES12.4)') 'ID,W:',ID,(W(I),I=1,4)
       SELECT CASE(ID)
       CASE(1)
          RKR=RKR_1
       CASE(2)
          RKR=RKR_2
       CASE(3)
          RKR=-RKR_1
       CASE(4)
          RKR=-RKR_2
       END SELECT
    CASE(1)
       RKR=RKR_1
    CASE(2)
       RKR=RKR_2
    CASE(3)
       RKR=-RKR_1
    CASE(4)
       RKR=-RKR_2
    END SELECT
       
    ! --- rkr0 is given for MODEW=0 ---
    CALL WRNWTN(IERR)  ! input RKR0,RKZI,RKPHII; output RKRI
    IF(IERR.NE.0) THEN
       ierr=300
       RETURN
    END IF

    WRITE(6,'(A,3ES12.4)') 'WRNWTN results: RKR0,RKRI,RKR=',RKR0,RKRI,RKR


    IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
       YA(1)= R
       YA(2)= PHI
       YA(3)= Z
       YA(4)= RKR
       YA(5)= RKPHI
       YA(6)= RKZ
    ELSE
       YA(1)= R*COS(PHI)
       YA(2)= R*SIN(PHI)
       YA(3)= Z
       YA(4)= RKR*COS(PHI)-RKPHI*SIN(PHI)
       YA(5)= RKR*SIN(PHI)+RKPHI*COS(PHI)
       YA(6)= RKZ
    ENDIF
    YA(7)= UUI
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

    CALL WRCALE(RAYS(0,0,NRAY),NSTPMAX_NRAY(NRAY),NRAY)

    RETURN
  END SUBROUTINE

!  --- original Runge-Kutta method ---

  SUBROUTINE WRRKFT(nstp_in,Y,YN,NNSTP)

    USE wrcomm
    USE wrsub,ONLY: DISPXR
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
    REAL(rkind):: X0,XE,OMG,PW,DELTA,RHON,R

    OMG=2.D6*PI*RF

    X0 = 0.D0
    XE = DELS
    NSTPLIM=MIN(INT(SMAX/DELS),NSTPMAX)
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
    REAL(rkind):: X0,XE,OMG,OXEFF,RHON,PW,RL,RKRL,DELTA,R
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I,IOX

    OMG=2.D6*PI*RF
    
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
    USE wrsub,ONLY: DISPXR
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
