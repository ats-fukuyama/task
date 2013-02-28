MODULE trprf
! ***********************************************************
!      RF HEATING AND CURRENT DRIVE (GAUSSIAN PROFILE)
! ***********************************************************
  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_power_rf

CONTAINS

  SUBROUTINE tr_power_rf
! -------------------------------------------------------------------------
! < indice >
!    tot: total power
!    rw : radial width of deposition
!    cd : current drive factor
!    r0 : radial position of deposition
!   toe : power partition to electron
!   npr : parallel refractive index
! -------------------------------------------------------------------------
    USE trcomm, ONLY : &
         ame,vc,rkev,nrmax,dvrho,epsrho,ra,rmnrho,rhog, &
         rn,rt,z_eff,jcd_ec,jcd_lh,jcd_ic,pec,plh,pic,  &
         peccd, pecnpr, pecr0, pecrw, pectoe, pectot,   &
         piccd, picnpr, picr0, picrw, pictoe, pictot,   &
         plhcd, plhnpr, plhr0, plhrw, plhtoe, plhtot

    IMPLICIT NONE

    REAL(rkind)   :: &
         efcdec, efcdic, efcdlh, &! current drive efficiency
         fact, pec0, pecl, pic0, picl, plh0, plhl,  &
         plhr0l, rlnlmd, sumec, sumic, sumlh,       &
         pecm,pecp,plhm,plhp,picm,picp,             &
         vphec, vphic, vphlh,    &! phase velocity
         vte, vtep, drhog
    INTEGER(ikind):: nr
    

    IF(pectot + plhtot + pictot <= 0.d0) RETURN

    IF(plhr0 < 0.d0) THEN
       vphlh = vc / (plhnpr*ABS(plhr0))
       DO nr = nrmax, 0, -1
          vte = SQRT(ABS(rt(1,nr))*rkev/ame)

          IF(vte > vphlh) THEN
             IF(nr == nrmax) THEN
                plhr0l = ra
             ELSE
                vtep   = SQRT(ABS(rt(1,nr+1))*rkev/ame)
                fact   = (vte-vphlh) / (vte-vtep)
                plhr0l = (fact*rhog(nr+1) + (1.d0-fact)*rhog(nr))*ra
             ENDIF
             GOTO 6
          ENDIF

       ENDDO
       plhr0l = 0.d0
6      CONTINUE
       !         WRITE(6,*) '*** PLHR0L = ',PLHR0L
    ELSE
       plhr0l = plhr0
       vphlh  = vc/plhnpr
    ENDIF


    sumec = 0.d0
    sumlh = 0.d0
    sumic = 0.d0

    DO nr = 1, nrmax
       drhog = rhog(nr) - rhog(nr-1)
       pecm = DEXP(-((rmnrho(nr-1)-pecr0 )/pecrw)**2)*dvrho(nr-1)
       pecp = DEXP(-((rmnrho(nr  )-pecr0 )/pecrw)**2)*dvrho(nr  )
       plhm = DEXP(-((rmnrho(nr-1)-plhr0l)/plhrw)**2)*dvrho(nr-1)
       plhp = DEXP(-((rmnrho(nr  )-plhr0l)/plhrw)**2)*dvrho(nr  )
       picm = DEXP(-((rmnrho(nr-1)-picr0 )/picrw)**2)*dvrho(nr-1)
       picp = DEXP(-((rmnrho(nr  )-picr0 )/picrw)**2)*dvrho(nr  )

       sumec = sumec + 0.5d0*(pecm + pecp)*drhog
       sumlh = sumlh + 0.5d0*(plhm + plhp)*drhog
       sumic = sumic + 0.5d0*(picm + picp)*drhog
    ENDDO

    pec0 = pectot*1.d6 / sumec
    plh0 = plhtot*1.d6 / sumlh
    pic0 = pictot*1.d6 / sumic

!      IF(ABS(PLHNPR).LE.1.D0) THEN
!         NLH=PLHR0/DR+1.D0
!         RTLH=((RM(NLH+1)-PLHR0)*RT(NLH,1)
!     &        +(PLHR0-RM(NLH))*RT(NLH+1,1))/DR
!         VTLH=SQRT(RTLH*RKEV/AME)
!         VPHLH=VTLH/ABS(PLHNPR)
!         IF(VPHLH.GT.VC) VPHLH=VC
!         IF(PLHNPR.LT.0.D0) VPHLH=-VPHLH
!         WRITE(6,*) '** PLHNPR = ',VC/VPHLH
!      ELSE
!         VPHLH=VC/PLHNPR
!      ENDIF

    radial : DO nr = 0, nrmax
       ! gaussian profile
       pecl = pec0 * DEXP(-((rmnrho(nr)-pecr0) /pecrw)**2)
       plhl = plh0 * DEXP(-((rmnrho(nr)-plhr0l)/plhrw)**2)
       picl = pic0 * DEXP(-((rmnrho(nr)-picr0) /picrw)**2)

       ! power to electron
       pec(1,nr) = pectoe * pecl
       plh(1,nr) = plhtoe * plhl
       pic(1,nr) = pictoe * picl

       ! power to ion
       pec(2,nr) = (1.d0 - pectoe) * pecl
       plh(2,nr) = (1.d0 - plhtoe) * plhl
       pic(2,nr) = (1.d0 - pictoe) * picl
       
       rlnlmd = 16.1D0 - 1.15D0*LOG10(rn(1,nr)) + 2.30d0*LOG10(rt(1,nr))
       vte    = SQRT(ABS(rt(1,nr))*rkev/ame)

       ! current drive efficiency (FUNCTION 'trcdef' is defined below.)
       IF(peccd /= 0.d0) THEN
          vphec = vc / (vte*pecnpr)
          IF(pecnpr <= 1.d0) THEN
             efcdec = 0.D0
          ELSE
             efcdec = trcdef(vphec,z_eff(nr),0.d0,epsrho(nr),0)
          ENDIF
       ELSE
          efcdec = 0.d0
       ENDIF

       IF(plhcd /= 0.d0) THEN
          IF(ABS(vphlh) > vc) THEN
             efcdlh = 0.d0
          ELSE
             efcdlh = trcdef(vphlh/vte,z_eff(nr),0.d0,epsrho(nr),0)
          ENDIF
       ELSE
          efcdlh = 0.d0
       ENDIF

       IF(piccd /= 0.d0) THEN
          vphic = vc / (vte*picnpr)
          IF(picnpr <= 1.d0) THEN
             efcdic = 0.d0
          ELSE
             efcdic = trcdef(vphic,z_eff(nr),0.d0,epsrho(nr),1)
          ENDIF
       ELSE
          efcdic = 0.d0
       ENDIF
!       write(6,*) 'rlnlmd, peccd, pectoe, efcdec, pecl'
!       write(6,*) rlnlmd,peccd,pectoe,efcdec,pecl
       jcd_ec(nr) = 0.384D0 * rt(1,nr)/(rn(1,nr)*rlnlmd)  &
                            * (peccd*pectoe*efcdec*pecl)
       jcd_lh(nr) = 0.384d0 * rt(1,nr)/(rn(1,nr)*rlnlmd)  &
                            * (plhcd*plhtoe*efcdlh*plhl)
       jcd_ic(nr) = 0.384d0 * rt(1,nr)/(rn(1,nr)*rlnlmd)  &
                            * (piccd*pictoe*efcdic*picl)
    ENDDO radial

    RETURN
  END SUBROUTINE tr_power_rf


  REAL(8) FUNCTION trcdef(WT,Z,XR,YR,ID)
! ------------------------------------------------------------------------
!      *** CURRENT DRIVE EFFICIENCY ***
!
!      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
!      Z  = ZEFF   : EFFECTIVE Z
!      XR = X / RR : NORMALIZED X
!      YR = Y / RR : NORMALIZED Y
!      ID : 0 : LANDAU DAMPING
!           1 : TTMP
! ------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind) :: ID
    REAL(rkind) :: WT,Z,XR,YR
    REAL(rkind) :: R,A,C,D,W,RM,RC,EFF0,EFF1,EFF2,EFF3,Y1,Y2,YT,ARG

    R=SQRT(XR*XR+YR*YR)
    IF(ID.EQ.0) THEN
       D=3.D0/Z
       C=3.83D0
       A=0.D0
       RM=1.38D0
       RC=0.389D0
    ELSE
       D=11.91D0/(0.678D0+Z)
       C=4.13D0
       A=12.3D0
       RM=2.48D0
       RC=0.0987D0
    ENDIF
    IF(WT.LE.1.D-20) THEN
       W=1.D-20
    ELSE
       W=WT
    ENDIF
    EFF0=D/W+C/Z**0.707D0+4.D0*W*W/(5.D0+Z)
    EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)
    
    Y2=(R+XR)/(1.D0+R)
    IF(Y2.LT.0.D0) Y2=0.D0
    Y1=SQRT(Y2)
    EFF2=1.D0+A*(Y1/W)**3
    
    IF(Y2.LE.1.D-20) THEN
       YT=(1.D0-Y2)*WT*WT/1.D-60
    ELSE
       YT=(1.D0-Y2)*WT*WT/Y2
    ENDIF
    IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
       ARG=(RC*YT)**RM
       IF(ARG.LE.100.D0) THEN
          EFF3=1.D0-MIN(EXP(-ARG),1.D0)
       ELSE
          EFF3=1.D0
       ENDIF
    ELSE
       EFF3=1.D0
    ENDIF
    
    TRCDEF=EFF0*EFF1*EFF2*EFF3
    RETURN
  END FUNCTION trcdef

END MODULE trprf
