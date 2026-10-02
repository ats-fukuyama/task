! wrstartpoint.f90

!   ***** Singel ray tracing module *****

  SUBROUTINE wr_setup_start_point(NRAY,YN,nstp,IERR)

    USE wrcomm
    USE wrsub,ONLY: wrcale,wrcale_xyz,wr_cold_rkperp,wr_newton
    USE plprof,ONLY: pl_mag_type,pl_mag,pl_prf_type,pl_prof
    USE plprofw,ONLY: pl_prfw_type,pl_profw
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: NRAY,nstp
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: IERR
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_prf_type),DIMENSION(nsmax):: plf
    REAL(rkind):: RF,RP,ZP,PHI,ANGT,ANGP,RNK,UU
    INTEGER:: MODEW,mode,icount
    REAL(rkind):: XP,YP,omega,rkv,s,deg,factor,omega_pe2,rne
    REAL(rkind):: rhon,rk,rk_new,rkpara,rkperp,rkperp1,rkperp2
    REAL(rkind):: rkr,rkph,rkz,rkx,rky,dXP,dYP,dZP
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
    SELECT CASE(MOD(mdlwri,10))
    CASE(1)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk*COS(angp*deg)*SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)
    CASE(2)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk              *SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)*COS(angt*deg)
    CASE(3)
       arg=1.D0-SIN(angt*deg)**2-SIN(angp*deg)**2
       IF(arg.GT.0.D0) THEN
          rkr= -rk*SQRT(arg)
       ELSE
          rkr=0.D0
       END IF
       rkph= rk              *SIN(angt*deg)
       rkz=  rk              *SIN(angp*deg)
    END SELECT
    rkx=rkr*COS(phi)-rkph*SIN(phi)
    rky=rkr*SIN(phi)+rkph*COS(phi)

    ! --- initial save ---

    nstp=0
    YN(0,nstp)= s
    YN(1,nstp)= XP
    YN(2,nstp)= YP
    YN(3,nstp)= ZP
    YN(4,nstp)= RKX
    YN(5,nstp)= RKY
    YN(6,nstp)= RKZ
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
       WRITE(6,'(A,3ES12.4)') '   rkr,rkph,rkz   =',RKR,RKPH,RKZ
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz    =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz    =',RKX/rkv,RKY/rkv,RKZ/rkv
       WRITE(6,'(A,2ES12.4)') '   uu,s           =',UU,S
       WRITE(6,'(A,3ES12.4)') '   rhon,rne,factor=',rhon,rne,factor
    END IF
    
    ! --- If normalized density is below the threshold ---
    
    IF(factor.LE.pne_threshold) THEN

       ! --- advance straight ray in vacuum ---

       dXP=dels*rkx/rk
       dYP=dels*rky/rk
       dZP=dels*rkz/rk

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
          YN(4,nstp)= RKX
          YN(5,nstp)= RKY
          YN(6,nstp)= RKZ
          YN(7,nstp)= UU
          
          ! --- If the ray is out of region, exit mode=11 ---

          IF(RP.GT.RMAX_WR.OR. &
             RP.LT.RMIN_WR.OR. &
             ZP.GT.ZMAX_WR.OR. &
             ZP.LT.ZMIN_WR.OR. &
             S.GT.SMAX) THEN

             WRITE(6,'(A,2I6,2ES12.4)') &
                  'wr_exec_ray_single: nray,nstp,R,Z=',NRAY,nstp,R,Z
             WRITE(6,'(A,I4,6ES12.4)') 'RK:',nstp,R,PHI,Z,RKR,RKPHI,RKZ
             WRITE(6,'(A,2ES12.4)') 'R: ',R,RMAX_WR
             WRITE(6,'(A,2ES12.4)') 'R: ',R,RMIN_WR
             WRITE(6,'(A,2ES12.4)') 'Z: ',Z,ZMAX_WR
             WRITE(6,'(A,2ES12.4)') 'Z: ',Z,ZMIN_WR
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
    SELECT CASE(MOD(mdlwri,10))
    CASE(1)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk*COS(angp*deg)*SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)
    CASE(2)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk              *SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)*COS(angt*deg)
    CASE(3)
       arg=1.D0-SIN(angt*deg)**2-SIN(angp*deg)**2
       IF(arg.GT.0.D0) THEN
          rkr= -rk*SQRT(arg)
       ELSE
          rkr=0.D0
       END IF
       rkph= rk              *SIN(angt*deg)
       rkz=  rk              *SIN(angp*deg)
    END SELECT
    RKX=RKR*COS(PHI)-RKPHI*SIN(PHI)
    RKY=RKR*SIN(PHI)+RKPHI*COS(PHI)

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    CALL pl_mag(XP,YP,ZP,mag)

    rkpara=rkx*mag%bnx+rky*mag%bny+rkz*mag%bnz
    CALL WR_COLD_RKPERP(RP,ZP,PHI,RKPARA,RKPERP_1,RKPERP_2)
    
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
    ut_Y=un_X*b_X-un_X*b_Z
    ut_Z=un_Z*b_Y-un_Y*b_X

    ut_R  = ut_X*COS(phi)+ut_Y*SIN(phi)
    ut_phi=-ut_X*SIN(phi)+ut_Y*COS(phi)

    ! normal, parallel, tangential  component of wave vector

    rk_n=RKX*un_X+RKY*un_Y+RKZ*un_Z
    rk_b=RKX*b_X +RKY*b_Y +RKZ*b_Z
    rk_t=RKX*ut_X+RKY*ut_Y+RKZ*ut_Z

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
       rkx=rk_x1
       rky=rk_y1
       rkz=rk_z1
    CASE(2)
       rkx=rk_x2
       rky=rk_y2
       rkz=rk_z2
    CASE DEFAULT
       WRITE(6,'(A,I4)') 'XX wr_exec_single_ray: MODEW is not 1 nor 2:', MODEW
       STOP
    END SELECT
    rkpara=rkx*mag%bnx+rky*mag%bny+rkz*mag%bnz

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
    WRITE(6,'(A,5ES12.4)') 'rk_xyzrp: ',rkx,rky,rkz,rk_R,rk_phi
    WRITE(6,'(A,5ES12.4)') 'rn_xyzrp: ',rkx*rnv,rky*rnv,rkz*rnv, &
                                        rk_R*rnv,rk_phi*rnv

    rkx=rk_r*COS(phi)-rk_ph*SIN(phi)
    rky=rk_r*SIN(phi)+rk_ph*COS(phi)

    YN(0,nstp)= s
    YN(1,nstp)= XP
    YN(2,nstp)= YP
    YN(3,nstp)= ZP
    YN(4,nstp)= RKX
    YN(5,nstp)= RKY
    YN(6,nstp)= RKZ
    YN(7,nstp)= UU
          
    RETURN
  END SUBROUTINE wr_setup_start_point
