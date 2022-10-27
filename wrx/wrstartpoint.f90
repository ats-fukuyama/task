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
    REAL(rkind):: XP,YP,omega,rkv,s
    REAL(rkind):: rhon,rk,rk_new,rkpara,rkperp,rkperp1,rkperp2
    REAL(rkind):: rkr,rkph,rkz,rkx,rky,dXP,dYP,dZP

    IERR=0

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
       WRITE(6,'(A,2ES12.4,I4)') 'rf,uu,modew=       ',rf,uu,modew
    END IF

    ! --- save RAYIN array ---

    RAYIN(1,NRAY)=RF
    RAYIN(2,NRAY)=RP
    RAYIN(3,NRAY)=ZP
    RAYIN(4,NRAY)=PHI
    RAYIN(5,NRAY)=RNK
    RAYIN(6,NRAY)=ANGP
    RAYIN(7,NRAY)=ANGT
    RAYIN(8,NRAY)=UU

    ! --- initial setup ---
    
    omega=2.D6*PI*RF
    rkv=omega/VC
    mode=0 !  0 : vacuum, 1: plasma, 2: started
           ! 11 : out of region, 12: over count 
    s=0.D0
    
    ! --- initial position and wave number ---

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    rk=rkv*rnk
    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rkr= -rk*COS(angp)*COS(angt)
       rkph=rk*COS(angp)*SIN(angt)
       rkz=  rk*SIN(angp)
    CASE(1,3,101,103)
       rkr= -rk*COS(angp)*COS(angt)
       rkph=rk*COS(angp)*SIN(angt)
       rkz=  rk*SIN(angp)*COS(angt)
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
    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,2I4)')     '*** idebug_wr(11): nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp    =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkr,rkoh,rkz=',RKR,RKPH,RKZ
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz =',RKR/rkv,RKPH/rkv,RKZ/rkv
       WRITE(6,'(A,2ES12.4)') '   UU,S           =',UU,S
    END IF


    ! --- setup start point ---

    CALL pl_mag(XP,YP,ZP,mag)
    rhon=mag%rhon
    CALL pl_prof(rhon,plf)
    CALL pl_profw(rhon,plfw)

    ! --- check electron density ---

    IF(plfw(1)%rn.LE.pne_threshold) THEN
       IF(idebug_wr(2).NE.0) THEN
          WRITE(6,'(A)') '*** idebug_wr(2): xp,yp,zp,rho,rn,pne='
          WRITE(6,'(6ES12.4)') xp,yp,zp,rhon,plfw(1)%rn,pne_threshold
       END IF

       ! --- start straight ray in vacuum ---

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
          
          IF(idebug_wr(11).NE.0) THEN
             WRITE(6,'(A,2I4)') '*** idebug_wr(11): nray,nstp=',nray,nstp
             WRITE(6,'(A,3ES12.4)') '   xp,yp,zp    =',XP,YP,ZP
             WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz =',RKX,RKY,RKZ
             WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz =',RKR/rkv,RKPH/rkv,RKZ/rkv
             WRITE(6,'(A,2ES12.4)') '   UU,S           =',UU,S
          END IF

          CALL pl_mag(XP,YP,ZP,mag)
          rhon=mag%rhon
          CALL pl_prof(rhon,plf)
          CALL pl_profw(rhon,plfw)

          IF(idebug_wr(3).NE.0) THEN
             WRITE(6,'(A,2I4)') '*** idebug_wr(3): nray,nstp=',nray,nstp
             WRITE(6,'(A)')     '  xp,yp,zp,rho,pn,pne='
             WRITE(6,'(6ES12.4)') xp,yp,zp,rhon,plfw(1)%rn,pne_threshold
          END IF

          ! --- check electron density ---

          IF(plfw(1)%rn.GT.pne_threshold) THEN
             mode=1
             EXIT
          END IF
          IF(RP.GT.RMAX_WR.OR. &
             RP.LT.RMIN_WR.OR. &
             ZP.GT.ZMAX_WR.OR. &
             ZP.LT.ZMIN_WR.OR. &
             S.GT.SMAX) THEN
             mode=11
             EXIT
          END IF
       END DO

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
       rk_new=rk
       mode=1
    END IF

    ! --- Solve cold plasma dispersion relation

    WRITE(6,*) 'mode=',mode
    rk=rkv*rnk
    icount=0
    DO WHILE(mode.EQ.1)
       SELECT CASE(mdlwri)
       CASE(0,2,100,102)
          rkr= -rk*COS(angp)*COS(angt)
          rkph=rk*COS(angp)*SIN(angt)
          rkz=  rk*SIN(angp)
       CASE(1,3,101,103)
          rkr= -rk*COS(angp)*COS(angt)
          rkph=rk*COS(angp)*SIN(angt)
          rkz=  rk*SIN(angp)*COS(angt)
       END SELECT
       rkx=rkr*COS(phi)-rkph*SIN(phi)
       rky=rkr*SIN(phi)+rkph*COS(phi)
       rkpara=rkx*mag%bnx+rky*mag%bny+rkz*mag%bnz

       IF(idebug_wr(4).NE.0) THEN
          WRITE(6,'(A,2I4)') '*** idebug_wr(4): nray,nstp=',nray,nstp
          WRITE(6,'(A)')     '  rkr,rkph,rkz,rkx,rky,rkpara='
          WRITE(6,'(6ES12.4)')  rkr,rkph,rkz,rkx,rky,rkpara
       END IF

       CALL wr_cold_rkperp(omega,RP,ZP,PHI,rkpara,rkperp1,rkperp2)

       IF(idebug_wr(5).NE.0) THEN
          WRITE(6,'(A,2I4)') '*** idebug_wr(5): nray,nstp=',nray,nstp
          WRITE(6,'(A)')     '  rkpara,rkperp1,rkperp2,rnpara,rnperp1,rnperp2='
          WRITE(6,'(6ES12.4)')  &
               rkpara,rkperp1,rkperp2,rkpara/rkv,rkperp1/rkv,rkperp2/rkv
       END IF

       SELECT CASE(modew)
       CASE(0,2) ! fast wave
          rkperp=rkperp1
       CASE(1,3)
          rkperp=rkperp2
       CASE DEFAULT
          WRITE(6,'(A,2I4)') &
               'XX wr_exec_single_ray: unknown modew: nray,modewin(nray)=', &
               nray,modew
          ierr=1012
          RETURN
       END SELECT

       rk_new=SQRT(rkpara**2+rkperp**2)
       IF(ABS(rk_new-rk).LT.epsnw) THEN
          mode=2
          IF(idebug_wr(6).NE.0) THEN
             WRITE(6,'(A,2I4)')  '*** idebug_wr(6): nray,nstp=',nray,nstp
             WRITE(6,'(A)')         '    rk,rk_new,residual='
             WRITE(6,'(A,3ES12.4)') '   ',rk,rk_new,ABS(rk_new-rk)
          END IF
          rk=rk_new
          EXIT
       END IF

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
          
    IF(idebug_wr(7).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(7): nray,nstp=',nray,nstp
       WRITE(6,'(A,ES12.4)') '      rk=',rk
    END IF
    ! --- confirm wave number by newton ---
    
    CALL wr_newton(omega,RP,ZP,PHI,angp,angt,rk,rk_new,ierr)
    IF(ierr.NE.0) THEN
       ierr=2000+ierr
       RETURN
    END IF

    IF(idebug_wr(10).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(10): nray,nstp=',nray,nstp
       WRITE(6,'(A,2ES12.4)') '      rk,rk_new=',rk,rk_new
    END IF

    XP=RP*COS(PHI)
    YP=RP*SIN(PHI)
    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rkr= -rk*COS(angp)*COS(angt)
       rkph=rk*COS(angp)*SIN(angt)
       rkz=  rk*SIN(angp)
    CASE(1,3,101,103)
       rkr= -rk*COS(angp)*COS(angt)
       rkph=rk*COS(angp)*SIN(angt)
       rkz=  rk*SIN(angp)*COS(angt)
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
          
    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(11): nray,nstp=',nray,nstp
       WRITE(6,'(A,3ES12.4)') '   xp,yp,zp    =',XP,YP,ZP
       WRITE(6,'(A,3ES12.4)') '   rkx,rky,rkz =',RKX,RKY,RKZ
       WRITE(6,'(A,3ES12.4)') '   rnx,rny,rnz =',RKR/rkv,RKPH/rkv,RKZ/rkv
       WRITE(6,'(A,2ES12.4)') '   UU,S           =',UU,S
    END IF

    RETURN
  END SUBROUTINE wr_setup_start_point
