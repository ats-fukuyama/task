! wrexecr.f90

MODULE wrexecr

  USE wrcomm,ONLY: rkind
  REAL(rkind):: omega,rkv,rnv
  INTEGER:: nray_exec,nstp_exec
  
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
       nray_exec=nray
       omega=2.D6*PI*RFIN(nray)
       rkv=omega/VC
       rnv=VC/omega
       WRITE(6,'(A,I4,3ES12.4)') 'nray,omega,rkv,rnv=',nray,omega,rkv,rnv
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

!   ***** setup start point *****

  SUBROUTINE wr_setup_start_point(NRAY,YN,nstp,IERR)

    USE wrcomm
    USE wrsub,ONLY: wrcale,wrcale_xyz,wr_cold_rkperp,wr_newton
    USE plprof,ONLY: pl_mag_type,pl_mag,pl_prf_type,pl_prof
    USE plprofw,ONLY: pl_prfw_type,pl_profw
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NRAY
    INTEGER,INTENT(INOUT):: nstp
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: IERR
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_prf_type),DIMENSION(nsmax):: plf
    REAL(rkind):: RF,RP,ZP,PHI,RNPH,ANGT,ANGP,RNK,UU
    INTEGER:: MODEW,mode,icount
    REAL(rkind):: XP,YP,s,deg,factor,omega_pe2,rne
    REAL(rkind):: rhon,rk,rk_new,rkpara,rkperp,rkperp_1,rkperp_2
    REAL(rkind):: rk_x,rk_y,rk_z,dXP,dYP,dZP
    REAL(rkind):: b_x,b_y,b_z,b_R,b_phi
    REAL(rkind):: ut_x,ut_y,ut_z,ut_R,ut_phi
    REAL(rkind):: un_x,un_y,un_z,un_R,un_phi
    REAL(rkind):: rk_b,rk_n,rk_t,rk_R,rk_phi
    REAL(rkind):: rk_x1,rk_y1,rk_z1,rk_R1,rk_phi1
    REAL(rkind):: rk_x2,rk_y2,rk_z2,rk_R2,rk_phi2
    REAL(rkind):: alpha1,alpha2,diff1,diff2

    IERR=0
    deg=PI/180.D0

    ! --- setup common values and initial values ---

    RF=RFIN(NRAY)
    RP=RPIN(NRAY)
    ZP=ZPIN(NRAY)
    PHI=PHIIN(NRAY)
    RNPH=RNPHIN(NRAY)
    ANGT=ANGTIN(NRAY)
    ANGP=ANGPIN(NRAY)
    RNK=RNKIN(NRAY)
    UU=UUIN(NRAY)
    MODEW=MODEWIN(NRAY)

    ! --- save initial vaariables in RAYIN ---

    RAYIN(1,NRAY)=RF
    RAYIN(2,NRAY)=RP
    RAYIN(3,NRAY)=ZP
    RAYIN(4,NRAY)=PHI
    RAYIN(5,NRAY)=RNK
    RAYIN(6,NRAY)=ANGP
    RAYIN(7,NRAY)=ANGT
    RAYIN(8,NRAY)=UU

    ! --- initial setup ---
    
    mode=0            ! status:  0 : vacuum, 1: plasma, 2: started
                      !         11 : out of region, 12: over count 
    s=0.D0            ! initial ray length
    
    ! --- initial position and wave number ---

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    rk=rkv*rnk
    WRITE(6,'(A,3ES12.4)') 'rk,rkv,rnk=',rk,rkv,rnk
    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rk_r  = -rk*COS(angp*deg)*COS(angt*deg)
       rk_phi=  rk*COS(angp*deg)*SIN(angt*deg)
       rk_z  =  rk*SIN(angp*deg)
    CASE(1,3,101,103)
       rk_r  = -rk*COS(angp*deg)*COS(angt*deg)
       rk_phi=  rk              *SIN(angt*deg)
       rk_z  =  rk*SIN(angp*deg)*COS(angt*deg)
    END SELECT
    rk_x=rk_r*COS(phi)-rk_phi*SIN(phi)
    rk_y=rk_r*SIN(phi)+rk_phi*COS(phi)

    WRITE(6,'(A,5ES12.4)') &
         'rk0:',rk_x,rk_y,rk_z,rk_r,rk_phi
    WRITE(6,'(A,5ES12.4)') &
         'rn0:',rk_x*rnv,rk_y*rnv,rk_z*rnv,rk_r*rnv,rk_phi*rnv
    
    ! --- initial save ---

    nstp=0
    YN(0,nstp)= s
    YN(1,nstp)= XP
    YN(2,nstp)= YP
    YN(3,nstp)= ZP
    YN(4,nstp)= RK_X
    YN(5,nstp)= RK_Y
    YN(6,nstp)= RK_Z
    YN(7,nstp)= UU
    
    ! --- set magnetic field and minor radius at the start point ---

    CALL pl_mag(XP,YP,ZP,mag)
    rhon=mag%rhon
    CALL pl_prof(rhon,plf)
    CALL pl_profw(rhon,plfw)

    ! --- set electron density ---

    rne=plfw(1)%rn
    omega_pe2=rne*1.D20*AEE*AEE/(AME*EPS0)
    factor=omega_pe2/omega**2

    IF(idebug_wr(1).NE.0) THEN
       WRITE(6,'(A,A,I4,I8)') '*** idebug_wr(1): wr_setup_start_point: ', &
            'nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   rf,rnk,rk      =',rf,rnk,rk
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp       =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkr,rkph,rkz   =',rk_r,rk_phi,rk_z
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz    =',rk_x,rk_y,rk_z
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz    =',rk_x/rkv,rk_y/rkv,rk_z/rkv
       WRITE(6,'(A,2ES12.4)') '   uu,s           =',UU,S
       WRITE(6,'(A,3ES12.4)') '   rhon,rne,factor=',rhon,rne,factor
    END IF
    
    ! --- If normalized density is below the threshold ---
    
    IF(factor.LE.pne_threshold) THEN

       ! --- advance straight ray in vacuum ---

       dXP=dels*rk_x/rk
       dYP=dels*rk_y/rk
       dZP=dels*rk_z/rk

       DO WHILE(mode.EQ.0)
          XP=XP+dXP
          YP=YP+dYP
          ZP=ZP+dZP
          RP=SQRT(XP*XP+YP*YP)
          PHI=ATAN2(YP,XP)
          s=s+dels

          nstp=nstp+1
          YN(0,nstp)= s
          YN(1,nstp)= XP
          YN(2,nstp)= YP
          YN(3,nstp)= ZP
          YN(4,nstp)= RK_X
          YN(5,nstp)= RK_Y
          YN(6,nstp)= RK_Z
          YN(7,nstp)= UU
          
          ! --- If the ray is out of region, exit mode=11 ---

          IF(RP.GT.RMAX_WR.OR. &
             RP.LT.RMIN_WR.OR. &
             ZP.GT.ZMAX_WR.OR. &
             ZP.LT.ZMIN_WR.OR. &
             S.GT.SMAX) THEN

             WRITE(6,'(A,2I6,2ES12.4)') &
                  'wr_exec_ray_single: nray,nstp,R,Z=',NRAY,nstp,RP,ZP
             WRITE(6,'(A,I4,6ES12.4)') 'RK:',nstp,RP,PHI,ZP,RK_R,RK_PHI,RK_Z
             WRITE(6,'(A,2ES12.4)') 'R: ',RP,RMAX_WR
             WRITE(6,'(A,2ES12.4)') 'R: ',RP,RMIN_WR
             WRITE(6,'(A,2ES12.4)') 'Z: ',ZP,ZMAX_WR
             WRITE(6,'(A,2ES12.4)') 'Z: ',ZP,ZMIN_WR
             WRITE(6,'(A,2ES12.4)') 'S: ',S,SMAX
             mode=11
             EXIT
          END IF

          ! --- set magnetic field and minor radius at the new point ---

          CALL pl_mag(XP,YP,ZP,mag)
          rhon=mag%rhon
          CALL pl_prof(rhon,plf)
          CALL pl_profw(rhon,plfw)

          ! --- set electron density ---

          rne=plfw(1)%rn
          omega_pe2=rne*1.D20*AEE*AEE/(AME*EPS0)
          factor=omega_pe2/omega**2

          IF(idebug_wr(2).NE.0) THEN
             WRITE(6,'(A,A)') '*** idebug_wr(2): wr_setup_start_point: ', &
                  'xp,yp,zp,rhon,rne,factor='
             WRITE(6,'(6ES12.4)') xp,yp,zp,rhon,rne,factor
          END IF
          
          ! --- If normalized density is above the threshold, exit mode=1---

          IF(factor.GT.pne_threshold) THEN
             mode=1
             EXIT
          END IF

       END DO
    ELSE
       mode=1
    END IF

    ! --- solve cold dispersion for given k_para ---
       
    rk=rkv*rnk
    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rk_r  =-rk*COS(angp*deg)*COS(angt*deg)
       rk_phi= rk*COS(angp*deg)*SIN(angt*deg)
       rk_z  = rk*SIN(angp*deg)
    CASE(1,3,101,103)
       rk_r  =-rk*COS(angp*deg)*COS(angt*deg)
       rk_phi= rk              *SIN(angt*deg)
       rk_z  = rk*SIN(angp*deg)*COS(angt*deg)
    END SELECT
    rk_x=rk_r*COS(phi)-rk_phi*SIN(phi)
    rk_y=rk_r*SIN(phi)+rk_phi*COS(phi)
    WRITE(6,'(A,5ES12.4)') &
         'rk1:',rk_x,rk_y,rk_z,rk_r,rk_phi
    WRITE(6,'(A,5ES12.4)') &
         'rn1:',rk_x*rnv,rk_y*rnv,rk_z*rnv,rk_r*rnv,rk_phi*rnv

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    CALL pl_mag(XP,YP,ZP,mag)

    rkpara=rk_x*mag%bnx+rk_y*mag%bny+rk_z*mag%bnz
    CALL WR_COLD_RKPERP(omega,RP,ZP,PHI,RKPARA,RKPERP_1,RKPERP_2)
    
    WRITE(6,'(A,1P3E12.4)') 'R,PHI,Z=',RP,phi,ZP
    WRITE(6,'(A,1P3E12.4)') &
         'RKRPARA,RKPERP_1,RKPERP_2=',RKPARA,RKPERP_1,RKPERP_2
    WRITE(6,'(A,1P3E12.4)') &
         'RKPARA,RNPARA=',RKPARA,RKPARA*VC/(2.D6*PI*RF)

    ! unit vectors

    ! parallel vector

    b_X=mag%bnx
    b_Y=mag%bny
    b_Z=mag%bnz

    b_R  = b_X*COS(phi)+b_Y*SIN(phi)
    b_phi=-b_X*SIN(phi)+b_Y*COS(phi)

    ! normal vector

    un_R  = b_Z/SQRT(b_Z**2+b_R**2)
    un_phi= 0.D0
    un_Z  =-b_R/SQRT(b_Z**2+b_R**2)
    
    un_X=un_R*COS(phi)-un_phi*SIN(phi)
    un_Y=un_R*SIN(phi)+un_phi*COS(phi)

    ! tangential vector

    ut_X=un_Y*b_Z-un_Z*b_Y
    ut_Y=un_Z*b_X-un_X*b_Z
    ut_Z=un_X*b_Y-un_Y*b_X

    ut_R  = ut_X*COS(phi)+ut_Y*SIN(phi)
    ut_phi=-ut_X*SIN(phi)+ut_Y*COS(phi)

    ! normal, parallel, tangential  component of wave vector

    rk_n=rk_x*un_X+rk_y*un_Y+rk_z*un_Z
    rk_b=rk_x*b_X +rk_y*b_Y +rk_z*b_Z
    rk_t=rk_x*ut_X+rk_y*ut_Y+rk_z*ut_Z

    ! new wave number vector #1

    diff1=rkperp_1**2-rk_t**2
    IF(diff1.GT.0.D0) THEN
       alpha1=SQRT(diff1/rk_n**2)
    ELSE
       alpha1=0.D0
    END IF
    rk_x1=rk_b*b_x+rk_t*ut_x+alpha1*rk_n*un_x
    rk_y1=rk_b*b_y+rk_t*ut_y+alpha1*rk_n*un_y
    rk_z1=rk_b*b_z+rk_t*ut_z+alpha1*rk_n*un_z
    
    ! new wave number vector #2

    diff2=rkperp_2**2-rk_t**2
    IF(diff2.GT.0.D0) THEN
       alpha2=SQRT(diff2/rk_n**2)
    ELSE
       alpha2=0.D0
    END IF
    rk_x2=rk_b*b_x+rk_t*ut_x+alpha2*rk_n*un_x
    rk_y2=rk_b*b_y+rk_t*ut_y+alpha2*rk_n*un_y
    rk_z2=rk_b*b_z+rk_t*ut_z+alpha2*rk_n*un_z
    
    rk_R1  = rk_x1*COS(phi)+rk_y1*SIN(phi)
    rk_phi1=-rk_x1*SIN(phi)+rk_y1*COS(phi)
    rk_R2  = rk_x2*COS(phi)+rk_y2*SIN(phi)
    rk_phi2=-rk_x2*SIN(phi)+rk_y2*COS(phi)

    SELECT CASE(MODEW)
    CASE(1)
       rk_x=rk_x1
       rk_y=rk_y1
       rk_z=rk_z1
       rk_R=rk_R1
       rk_phi=rk_phi1
    CASE(2)
       rk_x=rk_x2
       rk_y=rk_y2
       rk_z=rk_z2
       rk_R=rk_R2
       rk_phi=rk_phi2
    CASE DEFAULT
       WRITE(6,'(A,I4)') 'XX wr_exec_single_ray: MODEW is not 1 nor 2:', MODEW
       STOP
    END SELECT
    rkpara=rk_x*mag%bnx+rk_y*mag%bny+rk_z*mag%bnz

    WRITE(6,'(A,1ES12.4)') 'rkpara   :',rkpara
    WRITE(6,'(A,5ES12.4)') 'b_xyzrp  :',b_x,b_y,b_z,b_R,b_phi
    WRITE(6,'(A,5ES12.4)') 'ut_xyzrp :',ut_x,ut_y,ut_z,ut_R,ut_phi
    WRITE(6,'(A,5ES12.4)') 'un_xyzrp :',un_x,un_y,un_z,un_R,un_phi
    WRITE(6,'(A,3ES12.4)') 'rk_bnt   :',rk_b,rk_n,rk_t
    WRITE(6,'(A,3ES12.4)') 'rn_bnt   :',rk_b*rnv,rk_n*rnv,rk_t*rnv
    WRITE(6,'(A,2ES12.4)') 'alpha_12 :',alpha1,alpha2
    WRITE(6,'(A,5ES12.4)') 'rk1_xyzrp:',rk_x1,rk_y1,rk_z1,rk_R1,rk_phi1
    WRITE(6,'(A,5ES12.4)') 'rn1_xyzrp:',rk_x1*rnv,rk_y1*rnv,rk_z1*rnv, &
                                        rk_R1*rnv,rk_phi1*rnv
    WRITE(6,'(A,5ES12.4)') 'rk2_xyzrp:',rk_x2,rk_y2,rk_z2,rk_R2,rk_phi2
    WRITE(6,'(A,5ES12.4)') 'rn2_xyzrp:',rk_x2*rnv,rk_y2*rnv,rk_z2*rnv, &
                                        rk_R2*rnv,rk_phi2*rnv
    WRITE(6,'(A,5ES12.4)') 'rk_xyzrp: ',rk_x,rk_y,rk_z,rk_R,rk_phi
    WRITE(6,'(A,5ES12.4)') 'rn_xyzrp: ',rk_x*rnv,rk_y*rnv,rk_z*rnv, &
                                        rk_R*rnv,rk_phi*rnv

    rk_x=rk_r*COS(phi)-rk_phi*SIN(phi)
    rk_y=rk_r*SIN(phi)+rk_phi*COS(phi)

    YN(0,nstp)= s
    YN(1,nstp)= XP
    YN(2,nstp)= YP
    YN(3,nstp)= ZP
    YN(4,nstp)= rk_x
    YN(5,nstp)= rk_y
    YN(6,nstp)= rk_z
    YN(7,nstp)= UU
          
    RETURN
  END SUBROUTINE wr_setup_start_point

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
