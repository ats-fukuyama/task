! wrexecr.f90

MODULE wrexecr

  USE wrcomm,ONLY: rkind
  REAL(rkind):: omega,rkv,rnv
  INTEGER:: nray_exec
  
  PRIVATE
  PUBLIC wr_exec_rays
  
CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE wr_exec_rays(ierr)

    USE wrcomm
    USE wrcalpwr
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL:: time1,time2
    REAL(rkind):: RK,PABSN
    INTEGER:: NRAY,nstp,nsa

    nsamax_dp=nsamax_wr
    CALL GUTIME(TIME1)
    DO NRAY=1,NRAYMAX
       nray_exec=nray
       omega=2.D6*PI*RFIN(nray)
       rkv=omega/VC
       rnv=VC/omega
       WRITE(6,'(A,I4,3ES12.4)') 'nray,omega,rkv,rnv=',nray,omega,rkv,rnv

       CALL wr_setup_start_point(NRAY,RAYS(0,0,NRAY),nstp,IERR)
       nstpmax_nray(nray)=nstp
       IF(IERR.NE.0) CYCLE
       CALL wr_exec_single_ray(NRAY,RAYS(0,0,NRAY),nstp,IERR)
       nstpmax_nray(nray)=nstp
       IF(IERR.NE.0) CYCLE

       DO nsa=1,nsamax_wr
          DO nstp=0,nstpmax_nray(nray)
             pwr_nsa_nstp_nray(nsa,nstp,nray)=pwr_nsa_nstp(nsa,nstp)
          END DO
       END DO
       
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
    USE wrsub,ONLY: wrcale,wrcale_xyz,wr_cold_rkperp,wr_newton, &
         wr_write_line,wr_cal_ep
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
    REAL(rkind):: RF,RP,ZP,PHI,ANGT,ANGP,RNK,UU
    INTEGER:: MODEW,mode
    REAL(rkind):: XP,YP,s,deg,factor,omega_pe2,rne,arg,err
    REAL(rkind):: rhon,rkpara,rkperp_1,rkperp_2
    REAL(rkind):: rk,rk_x,rk_y,rk_z,dXP,dYP,dZP
    REAL(rkind):: ub_x,ub_y,ub_z,ub_R,ub_phi
    REAL(rkind):: ut_x,ut_y,ut_z,ut_R,ut_phi
    REAL(rkind):: un_x,un_y,un_z,un_R,un_phi
    REAL(rkind):: rk_b,rk_n,rk_t,rk_R,rk_phi
    REAL(rkind):: rk_x1,rk_y1,rk_z1,rk_R1,rk_phi1
    REAL(rkind):: rk_x2,rk_y2,rk_z2,rk_R2,rk_phi2
    REAL(rkind):: alpha_1,alpha_2,diff_1,diff_2
    COMPLEX(rkind):: cepola(3),cenorm(3)
    
    IERR=0
    deg=PI/180.D0

    ! --- setup common values and initial values ---

    RF=RFIN(NRAY)
    RP=RPIN(NRAY)
    ZP=ZPIN(NRAY)
    PHI=PHIIN(NRAY)
    SELECT CASE(mdlwri)
       CASE(11,12,13,31,32,33,111,112,113,131,132,133)
          ANGT=180.D0-ANGTIN(NRAY)
          ANGP=180.D0-ANGPIN(NRAY)
       CASE DEFAULT
          ANGT=ANGTIN(NRAY)
          ANGP=ANGPIN(NRAY)
       END SELECT
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

    SELECT CASE(MOD(mdlwri,10))
    CASE(1)
       rk_r  =-rk*COS(angp*deg)*COS(angt*deg)
       rk_phi= rk*COS(angp*deg)*SIN(angt*deg)
       rk_z  = rk*SIN(angp*deg)
    CASE(2)
       rk_r  =-rk*COS(angp*deg)*COS(angt*deg)
       rk_phi= rk              *SIN(angt*deg)
       rk_z  = rk*SIN(angp*deg)*COS(angt*deg)
    CASE(3)
       arg=1.D0-SIN(angt*deg)**2-SIN(angp*deg)**2
       IF(arg.GT.0.D0) THEN
          rk_r= -rk*SQRT(arg)
       ELSE
          rk_r=0.D0
       END IF
       rk_phi= rk              *SIN(angt*deg)
       rk_z  = rk              *SIN(angp*deg)
    END SELECT
    rk_x=rk_r*COS(phi)-rk_phi*SIN(phi)
    rk_y=rk_r*SIN(phi)+rk_phi*COS(phi)

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
    YN(8,nstp)= 0.D0
    CALL wr_write_line(NSTP,YN(0,NSTP),YN(1:NEQ,NSTP),YN(8,NSTP))
    
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
       WRITE(6,'(A,A,I4,I8,I4)') '*** idebug_wr(1): wr_setup_start_point: ', &
            'nray,nstp,mdlwri=',nray,nstp,mdlwri
       WRITE(6,'(A,3ES12.4)') '   rf,rnk,rk      =',rf,rnk,rk
       WRITE(6,'(A,3ES12.4)') '   rp,phi,zp      =',rp,phi,zp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp       =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rk,angt,angp   =',rk,angt,angp
       WRITE(6,'(A,3ES12.4)') '   rkr,rkph,rkz   =',rk_r,rk_phi,rk_z
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz    =',rk_x,rk_y,rk_z
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz    =',rk_x*rnv,rk_y*rnv,rk_z*rnv
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
          YN(8,nstp)= 0.D0
          CALL wr_write_line(NSTP,YN(0,NSTP),YN(1:NEQ,NSTP),YN(8,NSTP))

          
          ! --- If the ray is out of region, exit mode=11 ---

          IF(RP.GT.RMAX_WR.OR. &
             RP.LT.RMIN_WR.OR. &
             ZP.GT.ZMAX_WR.OR. &
             ZP.LT.ZMIN_WR.OR. &
             S.GT.SMAX) THEN

             WRITE(6,'(A,2I6,2ES12.4)') &
                  'wr_exec_ray_single: nray,nstp,R,Z=',NRAY,nstp,RP,ZP
             WRITE(6,'(A,I4,6ES12.4)') 'RK:',nstp,RP,PHI,ZP,RK_R,RK_PHI,RK_Z
             WRITE(6,'(A,3ES12.4)') 'R,min,max: ',RP,RMIN_WR,RMAX_WR
             WRITE(6,'(A,3ES12.4)') 'Z,min,max: ',ZP,ZMIN_WR,ZMAX_WR
             WRITE(6,'(A,2ES12.4)') 'S,max:     ',S,SMAX
             IERR=2
             RETURN
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
       
    rkpara=rk_x*mag%bnx+rk_y*mag%bny+rk_z*mag%bnz
    CALL WR_COLD_RKPERP(omega,RP,ZP,PHI,RKPARA,RKPERP_1,RKPERP_2)
    
    IF(idebug_wr(3).NE.0) THEN
       WRITE(6,'(A)') '*** idebug_wr(3): wr_setup_start_point: '
       WRITE(6,'(A,3ES12.4)') 'rhon,rne,wpe2/w2    =',rhon,rne,factor
       WRITE(6,'(A,3ES12.4)') 'R,PHI,Z=',RP,phi,ZP
       WRITE(6,'(A,3ES12.4)') &
            'RKPARA,RKPERP_1,RKPERP_2=',RKPARA,RKPERP_1,RKPERP_2
       WRITE(6,'(A,3ES12.4)') &
            'RNPARA,RNPERP_1,RNPERP_2=',RKPARA*rnv,RKPERP_1*rnv,RKPERP_2*rnv
    END IF

    ! unit vectors

    ! parallel vector

    ub_X=mag%bnx
    ub_Y=mag%bny
    ub_Z=mag%bnz

    ub_R  = ub_X*COS(phi)+ub_Y*SIN(phi)
    ub_phi=-ub_X*SIN(phi)+ub_Y*COS(phi)

    ! normal vector (assuming toroidal symmetry)

    un_R  = ub_Z/SQRT(ub_Z**2+ub_R**2)
    un_phi= 0.D0
    un_Z  =-ub_R/SQRT(ub_Z**2+ub_R**2)
    
    un_X=un_R*COS(phi)
    un_Y=un_R*SIN(phi)

    ! tangential vector (u_t = u_n x u_b)

    ut_X=un_Y*ub_Z-un_Z*ub_Y
    ut_Y=un_Z*ub_X-un_X*ub_Z
    ut_Z=un_X*ub_Y-un_Y*ub_X

    ut_R  = ut_X*COS(phi)+ut_Y*SIN(phi)
    ut_phi=-ut_X*SIN(phi)+ut_Y*COS(phi)

    ! normal, parallel, tangential  component of wave vector

    rk_n=rk_x*un_X+rk_y*un_Y+rk_z*un_Z
    rk_b=rk_x*ub_X+rk_y*ub_Y+rk_z*ub_Z
    rk_t=rk_x*ut_X+rk_y*ut_Y+rk_z*ut_Z

    ! new wave number vector #1

    diff_1=rkperp_1**2-rk_t**2
    IF(rk_n.NE.0.D0) THEN
       IF(diff_1.GT.0.D0) THEN
          alpha_1=SQRT(diff_1/rk_n**2)
       ELSE
          alpha_1=0.D0
       END IF
    ELSE
       alpha_1=0.D0
    END IF
    rk_x1=rk_b*ub_x+rk_t*ut_x+alpha_1*rk_n*un_x
    rk_y1=rk_b*ub_y+rk_t*ut_y+alpha_1*rk_n*un_y
    rk_z1=rk_b*ub_z+rk_t*ut_z+alpha_1*rk_n*un_z
    
    ! new wave number vector #2

    diff_2=rkperp_2**2-rk_t**2
    IF(rk_n.NE.0.D0) THEN
       IF(diff_2.GT.0.D0) THEN
          alpha_2=SQRT(diff_2/rk_n**2)
       ELSE
          alpha_2=0.D0
       END IF
    ELSE
       alpha_2=0.D0
    END IF
    rk_x2=rk_b*ub_x+rk_t*ut_x+alpha_2*rk_n*un_x
    rk_y2=rk_b*ub_y+rk_t*ut_y+alpha_2*rk_n*un_y
    rk_z2=rk_b*ub_z+rk_t*ut_z+alpha_2*rk_n*un_z
    
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

    IF(idebug_wr(4).NE.0) THEN
       WRITE(6,'(A)') '*** idebug_wr(4): wr_setup_start_point: '
       WRITE(6,'(A,3ES12.4)') 'rkpara,pp:',rkpara,rkperp_1,rkperp_2
       WRITE(6,'(A,5ES12.4)') 'ub_xyzrp :',ub_x,ub_y,ub_z,ub_R,ub_phi
       WRITE(6,'(A,5ES12.4)') 'ut_xyzrp :',ut_x,ut_y,ut_z,ut_R,ut_phi
       WRITE(6,'(A,5ES12.4)') 'un_xyzrp :',un_x,un_y,un_z,un_R,un_phi
       WRITE(6,'(A,3ES12.4)') 'rk_bnt   :',rk_b,rk_n,rk_t
       WRITE(6,'(A,3ES12.4)') 'rn_bnt   :',rk_b*rnv,rk_n*rnv,rk_t*rnv
       WRITE(6,'(A,2ES12.4)') 'diff_12  :',diff_1,diff_2
       WRITE(6,'(A,2ES12.4)') 'alpha_12 :',alpha_1,alpha_2
       WRITE(6,'(A,5ES12.4)') 'rk1_xyzrp:',rk_x1,rk_y1,rk_z1,rk_R1,rk_phi1
       WRITE(6,'(A,5ES12.4)') 'rn1_xyzrp:',rk_x1*rnv,rk_y1*rnv,rk_z1*rnv, &
                                           rk_R1*rnv,rk_phi1*rnv
       WRITE(6,'(A,5ES12.4)') 'rk2_xyzrp:',rk_x2,rk_y2,rk_z2,rk_R2,rk_phi2
       WRITE(6,'(A,5ES12.4)') 'rn2_xyzrp:',rk_x2*rnv,rk_y2*rnv,rk_z2*rnv, &
                                           rk_R2*rnv,rk_phi2*rnv
       WRITE(6,'(A,5ES12.4)') 'rk_xyzrp: ',rk_x,rk_y,rk_z,rk_R,rk_phi
       WRITE(6,'(A,5ES12.4)') 'rn_xyzrp: ',rk_x*rnv,rk_y*rnv,rk_z*rnv, &
                                           rk_R*rnv,rk_phi*rnv
       WRITE(6,'(A,3ES12.4)') 'rk_bnt:   ', &
            rk_x*ub_X+rk_y*ub_Y+rk_z*ub_Z, &
            rk_x*un_X+rk_y*un_Y+rk_z*un_Z, &
            rk_x*ut_X+rk_y*ut_Y+rk_z*ut_Z
       WRITE(6,'(A,3ES12.4)') 'rn_bnt:   ', &
            (rk_x*ub_X+rk_y*ub_Y+rk_z*ub_Z)*rnv, &
            (rk_x*un_X+rk_y*un_Y+rk_z*un_Z)*rnv, &
            (rk_x*ut_X+rk_y*ut_Y+rk_z*ut_Z)*rnv
    END IF

    YN(0,nstp)= s
    YN(1,nstp)= XP
    YN(2,nstp)= YP
    YN(3,nstp)= ZP
    YN(4,nstp)= rk_x
    YN(5,nstp)= rk_y
    YN(6,nstp)= rk_z
    YN(7,nstp)= UU
    YN(8,nstp)= 0.D0

    IF(idebug_wr(7).NE.0) THEN
       CALL wr_cal_ep(nstp,nray,cepola,cenorm,err)
       WRITE(6,'(A,2i6)') 'polarization: nstp,nray=',nstp,nray
       WRITE(6,'(A,6ES12.4)') 'CEXYZ:',cepola(1),cepola(2),cepola(3)
       WRITE(6,'(A,6ES12.4)') 'CEOXP:',cenorm(1),cenorm(2),cenorm(3)
       WRITE(6,'(A,6ES12.4)') &
            'CEOXP-ABS:',ABS(cenorm(1)),ABS(cenorm(2)),ABS(cenorm(3))
    END IF
    RETURN
  END SUBROUTINE wr_setup_start_point

  ! *** single ray tracing ***

  SUBROUTINE wr_exec_single_ray(NRAY,YN,nstp,IERR)
    USE wrcomm
    USE wroxb,ONLY: WRRKFT_RKF_OXB
    USE wrsub,ONLY: wrcale,wrcale_xyz,wr_cold_rkperp,wr_newton,wr_write_line
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NRAY
    INTEGER,INTENT(INOUT):: nstp
    REAL(rkind),INTENT(IN):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind):: RF,S,XP,YP,ZP,RKX,RKY,RKZ,UU,omega,rkv,rk
    INTEGER:: nstp1,nstp_end

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
    INITIAL_DRV=0

    IF(idebug_wr(8).NE.0) THEN
       rk=SQRT(rkx**2+rky**2+rkz**2)
       WRITE(6,'(A,A,I4,I8)') '*** idebug_wr(8): wr_exec_single_ray: ', &
            'nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp       =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz    =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnkx,rnky,rnkz =',RKX/rkv,RKY/rkv,RKZ/rkv
       WRITE(6,'(A,4ES12.4)') '   UU,S,rk,rn     =',UU,S,rk,rk/rkv
    END IF

    IF(MDLWRQ.EQ.0) THEN
       CALL WRRKFT(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSEIF(MDLWRQ.EQ.1) THEN
       CALL WRRKFT_WITHD0(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSEIF(MDLWRQ.EQ.2) THEN
       CALL WRRKFT_WITHMC(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSEIF(MDLWRQ.EQ.3) THEN
       CALL WRRKFT_RKF(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSEIF(MDLWRQ.EQ.4) THEN
       CALL WRSYMP(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSEIF(MDLWRQ.EQ.5) THEN
       CALL WRRKFT_ODE(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSEIF(MDLWRQ.EQ.6) THEN
       CALL WRRKFT_RKF_OXB(nstp,RAYS(:,:,NRAY),nstp_end)
    ELSE
       WRITE(6,*) 'XX WRCALC: unknown MDLWRQ =', MDLWRQ
       IERR=1
       nstp_end=nstp
       RETURN
    ENDIF

    DO NSTP1=0,nstp_end
       RAYRB1(NSTP1,NRAY)=0.D0
       RAYRB2(NSTP1,NRAY)=0.D0
    END DO

    CALL WRCALE(RF,RAYS(:,:,NRAY),nstp_end,NRAY)

    nstp=nstp_end

    RETURN
  END SUBROUTINE wr_exec_single_ray

!  --- original Runge-Kutta method ---

  SUBROUTINE WRRKFT(nstp_start,YN,nstp_end)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    USE wrfdrv
    USE librk
    USE wrsub,ONLY: wr_write_line
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ)
    REAL(rkind):: YM(NEQ),WORK(2,NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I
    REAL(rkind):: X0,XE,RHON,PW,R

    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_start
    X0 = YN(0,NSTP)
    XE = X0+DELS
    DO I=1,7
       Y(I)=YN(I,NSTP)
    ENDDO

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,wr_fdrv,X0,XE,1,Y,YM,WORK)

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=PW-YM(7)
       pwr_nsa_nstp(1,nstp)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       
       X0=XE
       XE=X0+DELS

       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             nstp_end = NSTP
             GOTO 11
          ENDIF
       ENDIF
    END DO
    nstp_end=NSTPLIM

11  IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT

!  --- original Runge-Kutta method with correction for D=0 ---

  SUBROUTINE WRRKFT_WITHD0(nstp_start,YN,nstp_end)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,WRMODNWTN,wr_write_line

    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    USE wrfdrv
    USE librk
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ)
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I
    REAL(rkind):: X0,XE,PW,DELTA,RHON,R,RF,omega

    RF=wr_nray_status%RF
    omega=2.D6*PI*RF

    NSTPLIM=MIN(INT(SMAX/DELS),NSTPMAX)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    nstp=nstp_start

    X0 = YN(0,NSTP)
    XE = X0+DELS
    DO I=1,7
       Y(I)=YN(I,nstp)
    ENDDO

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(11): wrrkft_withd0: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,xe,smax,dels =',X0,XE,SMAX,DELS
    END IF
       
    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(11): wrrkft_withd0: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,y1,y2,y3 =',X0,Y(1),Y(2),Y(3)
       WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',Y(4),Y(5),Y(6),Y(7)
    END IF

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,wr_fdrv,X0,XE,1,Y,YM,WORK)

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
       pwr_nsa_nstp(1,nstp)=PW-YM(7)
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
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA*RKAP) THEN
             nstp_end = NSTP
             GOTO 11
          ENDIF
       END IF
    END DO
    nstp_end=NSTPLIM
 
11  IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(nstp_end,X0,YM,YN(8,nstp_end))

    RETURN
  END SUBROUTINE WRRKFT_WITHD0

!  --- Runge-Kutta method using ODE library ---

  SUBROUTINE WRRKFT_ODE(nstp_start,YN,nstp_end)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    USE wrfdrv
    USE librk
    USE wrsub,ONLY: wr_write_line
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ)
    REAL(rkind):: YM(NEQ),WORK(2,NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I
    REAL(rkind):: X0,XE,PW,RHON,R

    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_start
    X0 = YN(0,NSTP)
    XE = X0+DELS     
    DO I=1,7
       Y(I)=YN(I,NSTP)
    ENDDO

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,wr_fdrv,X0,XE,1,Y,YM,WORK)
       
       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       pwr_nsa_nstp(1,nstp)=PW-YM(7)
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
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             nstp_end = NSTP
             GOTO 11
          ENDIF
       END IF
    END DO
    nstp_end=NSTPLIM

11  IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_ODE

!  --- Runge-Kutta method with tunneling of cutoff-resonant layer ---

  SUBROUTINE WRRKFT_WITHMC(nstp_start,YN,nstp_end)

    USE wrcomm
    USE wrsub,ONLY: DISPXR,WRMODCONV,WRMODNWTN,wr_write_line

    USE pllocal
    USE plprof,ONLY:PL_MAG_OLD
    USE wrfdrv
    USE librk
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ)
    REAL(rkind):: YM(NEQ),WORK(2,NEQ),YK(3),F(NEQ)
    REAL(rkind):: X0,XE,RF,omega,OXEFF,RHON,PW,RL,RKRL,DELTA,R
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,I,IOX

    RF=wr_nray_status%RF
    omega=2.D6*PI*RF
    
    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_start
    X0 = YN(0,NSTP)
    XE = X0+DELS
    DO I=1,7
       Y(I)=YN(I,NSTP)
    ENDDO

    IOX=0

    CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL ODERK(7,wr_fdrv,X0,XE,1,Y,YM,WORK)
       
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
       pwr_nsa_nstp(1,nstp)=PW-YM(7)
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
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 11
       ENDIF
       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       IF(RHON.GT.RB/RA) THEN
          nstp_end = NSTP
          GOTO 11
       ENDIF
    END DO
    nstp_end=NSTPLIM

11  IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_WITHMC

! --- Auto-step-size Runge-Kutta-F method ---

  SUBROUTINE WRRKFT_RKF(nstp_start,YN,nstp_end)

    USE wrcomm
    USE plprof,ONLY: pl_mag_old
    USE wrfdrv
    USE librkf
    USE wrsub,ONLY: wr_write_line
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,INIT,NDE,IER,I
    REAL(rkind):: RELERR,ABSERR,X0,XE,WORK0,PW,YM(NEQ),RHON,R
    REAL(rkind):: ESTERR(NEQ),WORK1(NEQ),WORK2(NEQ),WORK3(NEQ),WORK4(NEQ,11)

    RELERR = EPSRAY
    ABSERR = EPSRAY
    INIT = 1

    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_start
    X0 = YN(0,NSTP)
    XE = X0+DELS     
    DO I=1,7
       Y(I)=YN(I,NSTP)
    ENDDO

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(12): wrrkft_rkf: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,y1,y2,y3 =',X0,Y(1),Y(2),Y(3)
       WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',Y(4),Y(5),Y(6),Y(7)
    END IF

       CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL RKF(7,wr_fdrv,X0,XE,Y,INIT,RELERR,ABSERR,YM, &
                ESTERR,NDE,IER,WORK0,WORK1,WORK2,WORK3,WORK4)
       IF (IER .NE. 0) THEN
          WRITE(6,'(A,2I6)') 'XX wrrkft_rkf: NSTP,IER=',NSTP,IER
          nstp_end=nstp-1
          RETURN
       ENDIF

       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       pwr_nsa_nstp(1,nstp)=PW-YM(7)
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
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 8000
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL pl_mag_old(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             nstp_end = NSTP
             GOTO 8000
          END IF
       ENDIF
    END DO
    nstp_end=NSTPLIM
    WRITE(6,*) nstplim
     
8000 CONTINUE
    IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_RKF

! --- Symplectic method (not completed) ---

  SUBROUTINE WRSYMP(nstp_start,YN,nstp_end)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    USE wrfdrv
    USE libsympl
    USE wrsub,ONLY: wr_write_line
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ)
    REAL(rkind):: X,F(NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NLPMAX,NSTP,I,NLP,IERR
    REAL(rkind):: EPS,PW,ERROR,RHON,R

    NLPMAX=10
    EPS=1.D-6

    NSTPLIM=INT(SMAX/DELS)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_start
    X=YN(0,NSTP)
    DO I=1,7
       Y(I)=YN(I,NSTP)
    ENDDO

    CALL wr_write_line(NSTP,X,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL SYMPLECTIC(Y,DELS,wr_fdrvr,6,NLPMAX,EPS,NLP,ERROR,IERR)
       CALL wr_fdrv(0.D0,Y,F)
       X=X+DELS

       YN(0,NSTP)=X
       DO I=1,6
          YN(I,NSTP)=Y(I)
       ENDDO
       Y(7)=Y(7)+F(7)*DELS
       IF(Y(7).GT.0.D0) THEN
          YN(7,NSTP)=Y(7)
          YN(8,NSTP)=-F(7)*DELS
          pwr_nsa_nstp(1,nstp)=-F(7)*DELS
       ELSE
          YN(7,NSTP)=0.D0
          YN(8,NSTP)=Y(7)
          pwr_nsa_nstp(1,nstp)=-F(7)*DELS
       ENDIF

       CALL wr_write_line(NSTP,X,Y,YN(8,NSTP))

       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 11
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             nstp_end = NSTP
             GOTO 11
          ENDIF
       END IF
    END DO
    nstp_end=NSTPLIM

11  IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(NSTP,X,Y,YN(8,NSTP))
    RETURN
  END SUBROUTINE WRSYMP


END MODULE wrexecr
