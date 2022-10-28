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
    REAL(rkind):: RF,RP,ZP,PHI,RNK,RNKP,RNKT,ANGP,ANGT,UU
    INTEGER:: MODEW,mode,icount
    REAL(rkind):: XP,YP,omega,rkv,s,deg,factor,omega_pe2,rne
    REAL(rkind):: rhon,rk,rk_new,rkpara,rkperp,rkperp1,rkperp2
    REAL(rkind):: rkr,rkph,rkz,rkx,rky,dXP,dYP,dZP

    IERR=0
    deg=PI/180.D0

    ! --- setup common values and initial values ---

    RF=RFIN(NRAY)
    RP=RPIN(NRAY)
    ZP=ZPIN(NRAY)
    PHI=PHIIN(NRAY)
    RNK=RNKIN(NRAY)
    RNKP=RNKPIN(NRAY)
    RNKT=RNKTIN(NRAY)
    ANGP=ANGPIN(NRAY)
    ANGT=ANGTIN(NRAY)
    MODEW=MODEWIN(NRAY)
    UU=UUIN(NRAY)

    IF(idebug_wr(1).NE.0) THEN
       WRITE(6,'(A,I4)')         '*** idebug_wr(1): nray=',nray
       WRITE(6,'(A,3ES12.4)')    'rp,zp,phi=         ',rp,zp,phi
       WRITE(6,'(A,3ES12.4)')    'rnk,rnkp,rnkt=     ',rnk,rnkp,rnkt
       WRITE(6,'(A,3ES12.4)')    'rnk,angp,angt=     ',rnk,angp,angt
       WRITE(6,'(A,2ES12.4,I4)') 'rf,uu,modew=       ',rf,uu,modew
    END IF

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
    
    omega=2.D6*PI*RF  ! angular frequency
    rkv=omega/VC      ! vacuum wave number
    mode=0            ! status:  0 : vacuum, 1: plasma, 2: started
                      !         11 : out of region, 12: over count 
    s=0.D0            ! initial ray length
    
    ! --- initial position and wave number ---

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    rk=rkv*rnk
    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk*COS(angp*deg)*SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)
    CASE(1,3,101,103)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk*SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)*COS(angt*deg)
    END SELECT
    SELECT CASE(modew)
    CASE(2,3)
       rkr= -rkr
       rkph=-rkph
       rkz= -rkz
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
    
    IF(idebug_wr(2).NE.0) THEN
       WRITE(6,'(A,2I4)')     '*** idebug_wr(2): nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp    =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkr,rkph,rkz=',RKR,RKPH,RKZ
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz =',RKX/rkv,RKY/rkv,RKZ/rkv
       WRITE(6,'(A,2ES12.4)') '   UU,S           =',UU,S
    END IF

    ! --- set magnetic field and minor radius at the start point ---

    CALL pl_mag(XP,YP,ZP,mag)
    rhon=mag%rhon
    CALL pl_prof(rhon,plf)
    CALL pl_profw(rhon,plfw)

    ! --- set electron density ---

    rne=plfw(1)%rn
    omega_pe2=rne*1.D20*AEE*AEE/(AME*EPS0)
    factor=omega_pe2/omega**2

    IF(idebug_wr(3).NE.0) THEN
       WRITE(6,'(A)') '*** idebug_wr(3): xp,yp,zp,rhon,rne,factor='
       WRITE(6,'(6ES12.4)') xp,yp,zp,rhon,rne,factor
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
          
          IF(idebug_wr(4).NE.0) THEN
             WRITE(6,'(A,2I4)') '*** idebug_wr(4): nray,nstp=',nray,nstp
             WRITE(6,'(A,4ES12.4)') '   xp,yp,zp,s =',XP,YP,ZP,S
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

          IF(idebug_wr(5).NE.0) THEN
             WRITE(6,'(A)') '*** idebug_wr(5): xp,yp,zp,rhon,rne,factor='
             WRITE(6,'(6ES12.4)') xp,yp,zp,rhon,rne,factor
          END IF
          
          ! --- If normalized density is above the threshold, exit mode=1---

          IF(factor.GT.pne_threshold) THEN
             mode=1
             EXIT
          END IF

          ! --- If the ray is out of region, exit mode=11 ---

          IF(RP.GT.RMAX_WR.OR. &
             RP.LT.RMIN_WR.OR. &
             ZP.GT.ZMAX_WR.OR. &
             ZP.LT.ZMIN_WR.OR. &
             S.GT.SMAX) THEN
             mode=11
             EXIT
          END IF
       END DO

       ! --- If the ray is out of region, return with ierr=1011 ---

       IF(mode.EQ.11) THEN
          ierr=1011
          WRITE(6,'(A)') &
               'XX wr_exec_single_ray: out of range: nray,X,Y,Z,R,phi:'
          WRITE(6,'(A,I6,5ES12.4)') &
               '      ',nray,XP,YP,ZP,RP,phi
          WRITE(6,'(A,I6,5ES12.4)') &
               '      ',nray,RMAX_WR,RMIN_WR,ZMAX_WR,ZMIN_WR,S
          RETURN
       END IF

    ELSE
       mode=1
    END IF

    ! --- Solve cold plasma dispersion relation ---

    rk=rkv*rnk
    icount=0
    DO WHILE(mode.EQ.1)
       SELECT CASE(mdlwri)
       CASE(0,2,100,102)
          rkr= -rk*COS(angp*deg)*COS(angt*deg)
          rkph= rk*COS(angp*deg)*SIN(angt*deg)
          rkz=  rk*SIN(angp*deg)
       CASE(1,3,101,103)
          rkr= -rk*COS(angp*deg)*COS(angt*deg)
          rkph= rk*SIN(angt*deg)
          rkz=  rk*SIN(angp*deg)*COS(angt*deg)
       END SELECT
       SELECT CASE(modew)
       CASE(2,3)
          rkr =-rkr
          rkph=-rkph
          rkz =-rkz
       END SELECT
       rkx=rkr*COS(phi)-rkph*SIN(phi)
       rky=rkr*SIN(phi)+rkph*COS(phi)

       IF(idebug_wr(6).NE.0) THEN
          WRITE(6,'(A,2I4)')     '*** idebug_wr(6): nray,nstp=',nray,nstp
          WRITE(6,'(A,4ES12.4)') '   xp,yp,zp,s  =',XP,YP,ZP,S
          WRITE(6,'(A,3ES12.4)') '   rkr,rkph,rkz=',RKR,RKPH,RKZ
          WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz =',RKX,RKY,RKZ
          WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz =',RKX/rkv,RKY/rkv,RKZ/rkv
       END IF

       ! --- calculate perpendicular wave numbers for parallel wave number ---
       !        kperp_1**2 (slow wave) > kperp_2**2 (fast wave)
       
       rkpara=rkx*mag%bnx+rky*mag%bny+rkz*mag%bnz
       CALL wr_cold_rkperp(omega,RP,ZP,PHI,rkpara,rkperp1,rkperp2)

       IF(idebug_wr(7).NE.0) THEN
          WRITE(6,'(A,2I4)') '*** idebug_wr(7): nray,nstp=',nray,nstp
          WRITE(6,'(A)')     '  rkpara,rkperp1,rkperp2,rnpara,rnperp1,rnperp2='
          WRITE(6,'(6ES12.4)')  &
               rkpara,rkperp1,rkperp2,rkpara/rkv,rkperp1/rkv,rkperp2/rkv
          WRITE(6,'(4ES12.4)')  &
               SQRT(rkpara**2+rkperp1**2),SQRT(rkpara**2+rkperp2**2), &
               SQRT(rkpara**2+rkperp1**2)/rkv,SQRT(rkpara**2+rkperp2**2)/rkv
       END IF

       SELECT CASE(modew)
       CASE(0,2) ! slow wave
          rkperp=rkperp1
       CASE(1,3) ! fast wave
          rkperp=rkperp2
       CASE DEFAULT
          WRITE(6,'(A,2I4)') &
               'XX wr_exec_single_ray: unknown modew: nray,modewin(nray)=', &
               nray,modew
          ierr=1012
          RETURN
       END SELECT

       rk_new=SQRT(rkpara**2+rkperp**2)

       ! --- if wave number (absolute value) converged, EXIT ---
       
       IF(ABS(rk_new-rk).LT.epsnw) THEN
          mode=2
          IF(idebug_wr(8).NE.0) THEN
             WRITE(6,'(A,2I4)')  '*** idebug_wr(8): nray,nstp=',nray,nstp
             WRITE(6,'(A)')         '    rk,rk_new,residual='
             WRITE(6,'(A,3ES12.4)') '   ',rk,rk_new,ABS(rk_new-rk)
          END IF
          rk=rk_new
          EXIT
       END IF

       ! --- if loop count exceeds lmaxnw, EXIT with ierr=1013 ---
       
       IF(icount.GT.lmaxnw) THEN
          WRITE(6,'(A,2I4)') &
               'XX wr_exec_single_ray: icount exceed lmaxnw: nray,icount=', &
               nray,icount
          ierr=1013
          RETURN
       END IF

       icount=icount+1
       rk=rk_new
    END DO
          
    ! --- confirm wave number by newton ---
    
    CALL wr_newton(omega,RP,ZP,PHI,angp,angt,rk,rk_new,ierr)
    IF(ierr.NE.0) THEN
       ierr=2000+ierr
       RETURN
    END IF

    IF(idebug_wr(9).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(9): nray,nstp=',nray,nstp
       WRITE(6,'(A,2ES12.4)') '      rk,rk_new=',rk,rk_new
    END IF

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk*COS(angp*deg)*SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)
    CASE(1,3,101,103)
       rkr= -rk*COS(angp*deg)*COS(angt*deg)
       rkph= rk*SIN(angt*deg)
       rkz=  rk*SIN(angp*deg)*COS(angt*deg)
    END SELECT
    SELECT CASE(modew)
    CASE(2,3)
       rkr =-rkr
       rkph=-rkph
       rkz =-rkz
    END SELECT
    rkx=rkr*COS(phi)-rkph*SIN(phi)
    rky=rkr*SIN(phi)+rkph*COS(phi)

    YN(0,nstp)= s
    YN(1,nstp)= XP
    YN(2,nstp)= YP
    YN(3,nstp)= ZP
    YN(4,nstp)= RKX
    YN(5,nstp)= RKY
    YN(6,nstp)= RKZ
    YN(7,nstp)= UU
          
    IF(idebug_wr(10).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(10): nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp    =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz =',RKX/rkv,RKY/rkv,RKZ/rkv
       WRITE(6,'(A,2ES12.4)') '   UU,S           =',UU,S
    END IF

    RETURN
  END SUBROUTINE wr_setup_start_point
