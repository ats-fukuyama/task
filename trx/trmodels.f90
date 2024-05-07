!     $Id$

MODULE trmodels
  PRIVATE
  PUBLIC mbgb_driver,mmm95_driver,mmm71_driver

CONTAINS

! *** Mixed Bohm Gyro Bohm transport model driver ***

  SUBROUTINE mbgb_driver(RSL,RRL,ANEL,TEL,TIL,DTEL,DNEL,PAL,PZL, &
                         QPL,BBL,WEXBL,SL, &
                         ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB, &
                         ierr)
    USE trcomm, ONLY: &
         rkind,RT,NRMAX,MDLKAI,RKEV,AMP
    USE mixed_Bohm_gyro_Bohm, ONLY: mixed_model
    IMPLICIT NONE
    REAL(rkind),INTENT(IN) :: &
         RSL,RRL,ANEL,TEL,TIL,DTEL,DNEL,PAL,PZL,QPL,BBL,WEXBL,SL
    REAL(rkind),INTENT(OUT) :: ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB
    INTEGER,INTENT(OUT):: ierr
    REAL(rkind),DIMENSION(1):: &
            rminor,rmajor,tekev,tikev,q,btor,aimass,charge,wexbs, &
            grdte,grdne,shear, &
            chi_i_mix,themix,thdmix,thigb,thegb,thibohm,thebohm
    REAL(rkind):: zte_p8,zte_edge,VTI,GAMMA0,EXBfactor,SHRfactor,factor
    INTEGER:: NR8,NPOINTS,lflowshear

    rminor(1)=RSL             ! minor radius [m]
    rmajor(1)=RRL             ! major radius [m]
    tekev(1)=TEL              ! T_e [keV]
    tikev(1)=TIL              ! T_i [keV]
    q(1)=QPL                  ! safety-factor
    btor(1)=BBL               ! toroidal magnetic field [T]
    aimass(1)=PAL             ! average ion mass [AMU]
    charge(1)=PZL             ! charge number of main thermal ions
    wexbs(1)=WEXBL            ! ExB shearing rate [rad/s]
    grdte(1)=-RRL*DTEL/TEL    ! -R ( d T_e / d r ) / T_e
    grdne(1)=-RRL*DNEL/ANEL   ! -R ( d n_e / d r ) / n_e
    shear(1)=SL               !  r ( d q   / d r ) / q
    NR8=NINT(NRMAX*0.8d0)
    zte_p8=RT(NR8,1)          ! T_e(0.8a)
    zte_edge=RT(NRMAX,1)      ! T_e(a)
    npoints=1
    SELECT CASE(MDLKAI)
    CASE(140)
       lflowshear=0
       SHRfactor=1.D0
       EXBfactor=1.D0
    CASE(141)
       lflowshear=1
       SHRfactor=1.D0
       EXBfactor=1.D0
    CASE(142)
       lflowshear=0
       VTI=SQRT(2.D0*TIL*RKEV/(PAL*AMP))
       GAMMA0=VTI/(QPL*RRL)
       SHRfactor=1.D0
       EXBfactor=1.D0/(1.D0+(WEXBL/GAMMA0)**2)
    CASE(143)
       lflowshear=0
       VTI=SQRT(2.D0*TIL*RKEV/(PAL*AMP))
       GAMMA0=VTI/(QPL*RRL)
       SHRfactor=1.D0/MAX(1.D0,(SL-0.5d0)**2)
       EXBfactor=1.D0/(1.D0+(WEXBL/GAMMA0)**2)
    CASE DEFAULT
       lflowshear=0
       SHRfactor=1.D0
       EXBfactor=1.D0
    END SELECT

    call mixed_model ( &
         rminor,  rmajor,  tekev,   tikev,   q, &
         btor,    aimass,  charge,  wexbs, &
         grdte,   grdne,   shear, &
         zte_p8, zte_edge, npoints, &
         chi_i_mix,  themix,   thdmix, &
         thigb,   thegb,    thibohm, thebohm, &
         ierr, lflowshear)

    factor=EXBfactor*SHRfactor
    ADFFI=thdmix(1)   *factor
    ACHIE=themix(1)   *factor
    ACHII=chi_i_mix(1)*factor
    ACHIEB=thebohm(1) *factor
    ACHIIB=thibohm(1) *factor
    ACHIEGB=thegb(1)  *factor
    ACHIIGB=thigb(1)  *factor
    RETURN
  END SUBROUTINE mbgb_driver

! *** Multi Mode Model 1995 driver ***

  SUBROUTINE mmm95_driver(NR,CHIIW,DIFHW,CHIEW,DIFZW, &
                             CHIIRB,DIFHRB,CHIERB,DIFZRB, &
                             CHIIKB,DIFHKB,CHIEKB,DIFZKB,IERR)
    USE trcomm, ONLY: &
         rkind,PM,PZ,RN,PNSS,RT,PTS,RA,RR,NRMAX,MDLKAI,RM,RG,RNF, &
         PZC,ANC,PZFE,ANFE,QP,BB,WEXBP,RKPRHO,S,NSMAX,NFMAX
    USE modmmm95
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind),INTENT(OUT):: CHIIW,DIFHW,CHIEW,DIFZW
    REAL(rkind),INTENT(OUT):: CHIIRB,DIFHRB,CHIERB,DIFZRB
    REAL(rkind),INTENT(OUT):: CHIIKB,DIFHKB,CHIEKB,DIFZKB
    INTEGER,INTENT(OUT):: IERR
    INTEGER,PARAMETER:: matdim=5  ! dimension of transport matricies
    INTEGER,PARAMETER:: npoints=1 ! number of radial points
    REAL,DIMENSION(npoints):: &
         rminor,  rmajor,   elong, &
         dense,   densh,    densimp,  densfe, &
         xzeff,   tekev,    tikev,    q,       btor, &
         avezimp, amassimp, amasshyd, aimass,  wexbs, &
         grdne,   grdni,    grdnh,    grdnz,   grdte,   grdti,  grdq
    REAL,DIMENSION(npoints):: &
         thiig,   thdig,    theig,    thzig, &
         thirb,   thdrb,    therb,    thzrb, &
         thikb,   thdkb,    thekb,    thzkb
    REAL,DIMENSION(matdim,npoints):: &
         gamma, omega, velthi, vflux
    REAL,DIMENSION(matdim,matdim,npoints):: &
         difthi
    INTEGER:: nprout, lprint, lsuper, lreset
    INTEGER:: nerr
    INTEGER:: lswitch(8)
    REAL:: cswitch(25), fig(4), frb(4), fkb(4)
    INTEGER:: NS,NF
    REAL(rkind),DIMENSION(NSMAX):: RNL,RNTL,DRNL,DRNTL
    REAL(rkind),DIMENSION(NFMAX):: RNFL,DRNFL
    REAL(rkind):: PZCL,PZFEL,AMCL,AMFEL
    REAL(rkind):: PNEL,PNTEL,PNHL,PNZL,ANCL,ANFEL,PNFL,PTEL,PTIL
    REAL(rkind):: PNIL,PNTIL,ZEFFL,ZNIL,ZZL,AMZL,AMHL,AMIL,QPL,BBL,WEXBL,RSL
    REAL(rkind):: DNEDR,DNIDR,DNHDR,DNZDR,DTEDR,DTIDR

    IF(NR.EQ.NRMAX) THEN
       DO NS=1,NSMAX
          RNL(NS)=PNSS(NS)
          RNTL(NS)=PNSS(NS)*PTS(NS)
          DRNL(NS)=(PNSS(NS)-RN(NRMAX,NS))/(RA-RM(NRMAX))
          DRNTL(NS)=(PNSS(NS)*PTS(NS) &
                    -RN(NRMAX,NS)*RT(NRMAX,NS))/(RA-RM(NRMAX))
       END DO
       DO NF=1,NFMAX
          RNFL(NF)=0.D0
          DRNFL(NF)=(0.D0-RNF(NRMAX,NF))/(RA-RM(NRMAX))
       END DO
       PZCL=PZC(NR)
       PZFEL=PZFE(NR)
       ANCL=ANC(NR)
       ANFEL=ANFE(NR)
    ELSE
       DO NS=1,NSMAX
          RNL(NS)= 0.5D0*(RN(NR+1,NS)+RN(NR,  NS))
          RNTL(NS)=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS) &
                         +RN(NR,  NS)*RT(NR,  NS))
          DRNL(NS)=      (RN(NR+1,NS)-RN(NR,  NS))/(RM(NR+1)-RM(NR))
          DRNTL(NS)=     (RN(NR+1,NS)*RT(NR+1,NS) &
                         -RN(NR,  NS)*RT(NR,  NS))/(RM(NR+1)-RM(NR))
       END DO
       DO NF=1,NFMAX
          RNFL(NF)=0.5D0*(RNF(NR+1,NF)+RNF(NR,NF))
          DRNFL(NF)=     (RNF(NR+1,NF)-RNF(NR,NF))/(RM(NR+1)-RM(NR))
       END DO
       PZCL =0.5D0*(PZC(NR) +PZC(NR+1) )
       PZFEL=0.5D0*(PZFE(NR)+PZFE(NR+1))
       ANCL =0.5D0*(ANC(NR) +ANC(NR+1) )
       ANFEL=0.5D0*(ANFE(NR)+ANFE(NR+1))
    END IF

    PNEL=0.D0
    PNTEL=0.D0
    PNHL=0.D0
    PNZL=0.D0
    PNIL=0.D0
    PNTIL=0.D0
    ZEFFL=0.D0
    ZNIL=0.D0
    ZZL=0.D0
    AMZL=0.D0
    AMHL=0.D0
    AMIL=0.D0

    DNEDR=0.D0
    DNIDR=0.D0
    DNHDR=0.D0
    DNZDR=0.D0
    DTEDR=0.D0
    DTIDR=0.D0

    DO NS=1,NSMAX
       IF(PM(NS).LE.0.01D0) THEN
          PNEL=PNEL+RNL(NS)
          PNTEL=PNTEL+RNTL(NS)
          DNEDR=DNEDR+DRNL(NS)
          DTEDR=DTEDR+DRNTL(NS)
       ELSE
          IF(PM(NS).LE.3.D0.AND.PZ(NS).EQ.1.D0) THEN  
             ! 3He+ and 3H+ are not distinguished; waiting for the use of NPA
             PNHL=PNHL+RNL(NS)
             AMHL=AMHL+PM(NS)*RNL(NS)
             DNHDR=DNHDR+DRNL(NS)
          ELSE
             PNZL=PNZL+RNL(NS)
             ZZL =ZZL +PZ(NS)*RNL(NS)
             AMZL=AMZL+PM(NS)*RNL(NS)
             DNZDR=DNZDR+DRNL(NS)
          END IF
          ZEFFL=ZEFFL+PZ(NS)*PZ(NS)*RNL(NS)
          ZNIL =ZNIL +PZ(NS)*RNL(NS)
          PNIL =PNIL +RNL(NS)
          PNTIL=PNTIL+RNTL(NS)
          AMIL =AMIL +PM(NS)*RNL(NS)
          DNIDR=DNIDR+DRNL(NS)
          DTIDR=DTIDR+DRNTL(NS)
       END IF
    END DO
    PNFL=RNFL(1)+2.D0*RNFL(2)  ! NF=1: H or D, NF=2: He

    ZEFFL=ZEFFL+PZCL**2*ANCL+PZFEL**2*ANFEL
    ZNIL =ZNIL +PZCL   *ANCL+PZFEL   *ANFEL
    ZEFFL=ZEFFL/ZNIL

    PTEL=PNTEL/PNEL
    PTIL=PNTIL/PNIL
    QPL=QP(NR)
    BBL=BB
    RSL=RA*RG(NR)

    AMHL=AMHL/PNHL
    AMIL=AMIL/PNIL
    IF(PNZL.EQ.0.D0) THEN
       ZZL=ZEFFL
       AMZL=AMIL
    ELSE
       ZZL=ZZL/PNZL
       AMZL=AMZL/PNZL
    END IF
    SELECT CASE(MDLKAI)
    CASE(150)
       WEXBL=0.D0
    CASE(151:159)
       WEXBL=WEXBP(NR)
    END SELECT
    DTEDR=DTEDR/PNEL
    DTIDR=DTIDR/PNIL

!   WRITE(6,'(A,I5,1P6E12.4)') '1: ',NR,RSL,RR,RKPRHO(NR),QPL,BBL
!   WRITE(6,'(A,I5,1P6E12.4)') '2: ',NR,PNEL,PNHL,PNZL,PNFL,PTEL,PTIL
!   WRITE(6,'(A,I5,1P6E12.4)') '3: ',NR,ZEFFL,ZZL,AMZL,AMHL,AMIL,WEXBL
!   WRITE(6,'(A,I5,1P6E12.4)') '4: ',NR,DNEDR,DNIDR,DNHDR,DNZDR,DTEDR,DTIDR

    rminor(1)=RSL      ! minor radius (half-width) of zone boundary [m]
    rmajor(1)=RR       ! major radius to geometric center of zone bndry [m]
    elong(1)=RKPRHO(NR)! local elongation of zone boundary

    dense(1)=PNEL*1.D20! electron density [m^-3]
    densh(1)=PNHL*1.D20! sum over thermal hydrogenic ion densities [m^-3]
    IF(PNZL.LE.1.D-8) THEN
       densimp(1)=1.D12 ! nonzero dummy impurity density
    ELSE
       densimp(1)=PNZL*1.D20! sum over impurity ion densities [m^-3]
    END IF
    densfe(1)=PNFL*1.D20 ! electron density from fast (non-thermal) ions [m^-3]
    xzeff(1)=ZEFFL     ! Z_eff
    tekev(1)=PTEL      ! T_e (electron temperature) [keV] 
    tikev(1)=PTIL      ! T_i (temperature of thermal ions) [keV]
    q(1)=QPL           ! magnetic q-value
    btor(1)=BBL        ! ( R B_tor ) / rmajor(jz)  [tesla]
    avezimp(1)=ZZL     ! ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp )
    amassimp(1)=AMZL   ! ( sum_imp n_imp M_imp ) / ( sum_imp n_imp )
    amasshyd(1)=AMHL   ! ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd )
    aimass(1)=AMIL     ! ( sum_i n_i M_i ) / ( sum_i n_i )
    wexbs(1)=WEXBL     ! ExB shearing rate in [rad/s]
    grdne(1)=-RR*DNEDR/PNEL ! -R ( d n_e / d r ) / n_e
    grdni(1)=-RR*DNIDR/PNIL ! -R ( d n_i / d r ) / n_i
    grdnh(1)=-RR*DNHDR/PNHL ! -R ( d n_h / d r ) / n_h
    IF(PNZL.EQ.0.D0) THEN
       grdnz(1)=grdnh(1)
    ELSE
       grdnz(1)=-RR*DNZDR/PNZL ! -R ( d Z n_Z / d r ) / ( Z n_Z )
    END IF
    grdte(1)=-RR*DTEDR/PTEL ! -R ( d T_e / d r ) / T_e
    grdti(1)=-RR*DTIDR/PTIL ! -R ( d T_i / d r ) / T_i
    grdq(1)=RR/RSL*S(NR)    !  R ( d q   / d r ) / q  related to magnetic shear

    nprout=21          ! output unit number for long printout
    lprint=0           ! controls the amount of printout (0 => no printout)
    lsuper=0           ! 0 for non supershot 
    lreset=0           ! 0 to use internal settings for lswitch, cswitch and f**

    CALL mmm95 ( &
         rminor,  rmajor,   elong &
         , dense,   densh,    densimp,  densfe &
         , xzeff,   tekev,    tikev,    q,       btor &
         , avezimp, amassimp, amasshyd, aimass,  wexbs &
         , grdne,   grdni,    grdnh,    grdnz,   grdte,   grdti,  grdq &
         , thiig,   thdig,    theig,    thzig &
         , thirb,   thdrb,    therb,    thzrb &
         , thikb,   thdkb,    thekb,    thzkb &
         , gamma,   omega,    difthi,   velthi,  vflux &
         , matdim,  npoints,  nprout,   lprint,  nerr &
         , lsuper,  lreset,   lswitch,  cswitch, fig,    frb,     fkb)

    IERR=nerr        ! status code returned; 0 = OK; .ne. 0 indicates error
    CHIIW=thiig(1)  ! ion thermal diffusivity from the Weiland model
    DIFHW=thdig(1)  ! hydrogenic ion diffusivity from the Weiland model
    CHIEW=theig(1)  ! elelctron thermal diffusivity from the Weiland model
    DIFZW=thzig(1)  ! impurity ion diffusivity from the Weiland model

    CHIIRB=thirb(1) ! ion thermal diffusivity from resistive ballooning modes
    DIFHRB=thdrb(1) ! hydrogenic ion diffusivity from resistive ballooning md.
    CHIERB=therb(1) ! elelctron thermal diffusivity from resistive ball. modes
    DIFZRB=thzrb(1) ! impurity ion diffusivity from resistive ballooning modes

    CHIIKB=thikb(1) ! ion thermal diffusivity from kinetic ballooning modes
    DIFHKB=thdkb(1) ! hydrogenic ion diffusivity from kinetic ballooning modes
    CHIEKB=thekb(1) ! elelctron thermal diffusivity from kinetic ball. modes
    DIFZKB=thzkb(1) ! impurity ion diffusivity from kinetic ballooning modes

    IF(IERR.NE.0) THEN
       WRITE(6,'(A,I5)') 'XX mmm95_driver: IERR=',IERR
    END IF

    RETURN
  END SUBROUTINE mmm95_driver

! *** Multi Mode Model version 7.1 driver ***

  SUBROUTINE mmm71_driver(NR,CHII,DIFH,CHIE,DIFZ,VIST,VISP, &
                             CHIIW,DIFHW,CHIEW,CHIIB,DIFHB,CHIEB,CHIEG, &
                             IERR)
    USE trcomm, ONLY: &
         rkind,PM,PZ,RN,PNSS,RT,PTS,RA,RR,NRMAX,MDLKAI,RM,RG,RNF, &
         PZC,ANC,PZFE,ANFE,QP,BB,WEXBP,RKPRHO,S,NSMAX,NFMAX
    USE modmmm7_1
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind),INTENT(OUT):: CHII,DIFH,CHIE,DIFZ,VIST,VISP
    REAL(rkind),INTENT(OUT):: CHIIW,DIFHW,CHIEW,CHIIB,DIFHB,CHIEB,CHIEG
    INTEGER,INTENT(OUT):: IERR

    INTEGER,PARAMETER:: npoints=1  ! Number of values in all of the 1-D arrays
!    INTEGER,PARAMETER:: MMM_NCH=6  ! Maximum number of transport channels
    Real(rkind), Dimension(npoints) :: &
         rmin, rmaj, elong, ne, nh, nz, nf, zeff, &
         te, ti, q, btor, zimp, aimp, ahyd, aimass, wexbs, &
         gne, gni, gnh, gnz, gte, gti, gq, &
         gvtor, vtor, gvpol, vpol, gvpar, vpar
    Real(rkind), Dimension(npoints) :: & ![m^2/s]
         xti, xdi, xte, xdz, xvt, xvp
    Real(rkind), Dimension(npoints) :: & ![m^2/s]
         xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG
    Real(rkind), Dimension(4,npoints) :: &
         gammaW20, omegaW20
    Real(rkind), Dimension(npoints) :: &
         gammaDBM, omegaDBM
    Real(rkind), Dimension(MMM_NCH,npoints) :: &
         vconv, vflux
    Integer:: lprint,nprout,nerr
    Real(rkind):: cmodel(3),cswitch(MAXNOPT,4)
    INTEGER:: lswitch(MAXNOPT,4)

    INTEGER:: NS,NFL,jz
    REAL(rkind),DIMENSION(NSMAX):: RNL,RNTL,DRNL,DRNTL
    REAL(rkind),DIMENSION(NFMAX):: RNFL,DRNFL
    REAL(rkind):: PZCL,PZFEL,AMCL,AMFEL
    REAL(rkind):: PNEL,PNTEL,PNHL,PNZL,ANCL,ANFEL,PNFL,PTEL,PTIL
    REAL(rkind):: PNIL,PNTIL,ZEFFL,ZNIL,ZZL,AMZL,AMHL,AMIL,QPL,BBL,WEXBL,RSL
    REAL(rkind):: DNEDR,DNIDR,DNHDR,DNZDR,DTEDR,DTIDR

    IF(NR.EQ.NRMAX) THEN
       DO NS=1,NSMAX
          RNL(NS)=PNSS(NS)
          RNTL(NS)=PNSS(NS)*PTS(NS)
          DRNL(NS)=(PNSS(NS)-RN(NRMAX,NS))/(RA-RM(NRMAX))
          DRNTL(NS)=(PNSS(NS)*PTS(NS) &
                    -RN(NRMAX,NS)*RT(NRMAX,NS))/(RA-RM(NRMAX))
       END DO
       DO NFL=1,NFMAX
          RNFL(NFL)=0.D0
          DRNFL(NFL)=(0.D0-RNF(NRMAX,NFL))/(RA-RM(NRMAX))
       END DO
       PZCL=PZC(NR)
       PZFEL=PZFE(NR)
       ANCL=ANC(NR)
       ANFEL=ANFE(NR)
    ELSE
       DO NS=1,NSMAX
          RNL(NS)= 0.5D0*(RN(NR+1,NS)+RN(NR,  NS))
          RNTL(NS)=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS) &
                         +RN(NR,  NS)*RT(NR,  NS))
          DRNL(NS)=      (RN(NR+1,NS)-RN(NR,  NS))/(RM(NR+1)-RM(NR))
          DRNTL(NS)=     (RN(NR+1,NS)*RT(NR+1,NS) &
                         -RN(NR,  NS)*RT(NR,  NS))/(RM(NR+1)-RM(NR))
       END DO
       DO NFL=1,NFMAX
          RNFL(NFL)=0.5D0*(RNF(NR+1,NFL)+RNF(NR,NFL))
          DRNFL(NFL)=     (RNF(NR+1,NFL)-RNF(NR,NFL))/(RM(NR+1)-RM(NR))
       END DO
       PZCL =0.5D0*(PZC(NR) +PZC(NR+1) )
       PZFEL=0.5D0*(PZFE(NR)+PZFE(NR+1))
       ANCL =0.5D0*(ANC(NR) +ANC(NR+1) )
       ANFEL=0.5D0*(ANFE(NR)+ANFE(NR+1))
    END IF

    PNEL=0.D0
    PNTEL=0.D0
    PNHL=0.D0
    PNZL=0.D0
    PNIL=0.D0
    PNTIL=0.D0
    ZEFFL=0.D0
    ZNIL=0.D0
    ZZL=0.D0
    AMZL=0.D0
    AMHL=0.D0
    AMIL=0.D0

    DNEDR=0.D0
    DNIDR=0.D0
    DNHDR=0.D0
    DNZDR=0.D0
    DTEDR=0.D0
    DTIDR=0.D0

    DO NS=1,NSMAX
       IF(PM(NS).LE.0.01D0) THEN
          PNEL=PNEL+RNL(NS)
          PNTEL=PNTEL+RNTL(NS)
          DNEDR=DNEDR+DRNL(NS)
          DTEDR=DTEDR+DRNTL(NS)
       ELSE
          IF(PM(NS).LE.3.D0.AND.PZ(NS).EQ.1.D0) THEN  
             ! 3He+ and 3H+ are not distinguished; waiting for the use of NPA
             PNHL=PNHL+RNL(NS)
             AMHL=AMHL+PM(NS)*RNL(NS)
             DNHDR=DNHDR+DRNL(NS)
          ELSE
             PNZL=PNZL+RNL(NS)
             ZZL =ZZL +PZ(NS)*RNL(NS)
             AMZL=AMZL+PM(NS)*RNL(NS)
             DNZDR=DNZDR+DRNL(NS)
          END IF
          ZEFFL=ZEFFL+PZ(NS)*PZ(NS)*RNL(NS)
          ZNIL =ZNIL +PZ(NS)*RNL(NS)
          PNIL =PNIL +RNL(NS)
          PNTIL=PNTIL+RNTL(NS)
          AMIL =AMIL +PM(NS)*RNL(NS)
          DNIDR=DNIDR+DRNL(NS)
          DTIDR=DTIDR+DRNTL(NS)
       END IF
    END DO
    PNFL=RNFL(1)+2.D0*RNFL(2)  ! NF=1: H or D, NF=2: He

    ZEFFL=ZEFFL+PZCL**2*ANCL+PZFEL**2*ANFEL
    ZNIL =ZNIL +PZCL   *ANCL+PZFEL   *ANFEL
    ZEFFL=ZEFFL/ZNIL

    PTEL=PNTEL/PNEL
    PTIL=PNTIL/PNIL
    QPL=QP(NR)
    BBL=BB
    RSL=RA*RG(NR)

    AMHL=AMHL/PNHL
    AMIL=AMIL/PNIL
    IF(PNZL.EQ.0.D0) THEN
       ZZL=ZEFFL
       AMZL=AMIL
    ELSE
       ZZL=ZZL/PNZL
       AMZL=AMZL/PNZL
    END IF
    SELECT CASE(MDLKAI)
    CASE(160)
       WEXBL=0.D0
    CASE(161:169)
       WEXBL=WEXBP(NR)
    END SELECT
    DTEDR=DTEDR/PNEL
    DTIDR=DTIDR/PNIL

!   WRITE(6,'(A,I5,1P6E12.4)') '1: ',NR,RSL,RR,RKPRHO(NR),QPL,BBL
!   WRITE(6,'(A,I5,1P6E12.4)') '2: ',NR,PNEL,PNHL,PNZL,PNFL,PTEL,PTIL
!   WRITE(6,'(A,I5,1P6E12.4)') '3: ',NR,ZEFFL,ZZL,AMZL,AMHL,AMIL,WEXBL
!   WRITE(6,'(A,I5,1P6E12.4)') '4: ',NR,DNEDR,DNIDR,DNHDR,DNZDR,DTEDR,DTIDR

   jz=1 
   rmin(jz)=RSL ! half-width of the magnetic surface [m]
   rmaj(jz)=RR ! major radius to geometric center of the magnetic surface [m] 
   elong(jz)=RKPRHO(NR) ! local elongation

   ne(jz)=PNEL*1.D20  ! electron density [m^-3]
   nh(jz)=PNHL*1.D20  ! hydrogenic thermal particle density [m^-3]
   IF(PNZL.EQ.0.D0) THEN
      nz(jz)=1.D12    ! dummy impurity density
   ELSE
      nz(jz)=PNZL*1.D20  ! impurity ion density [m^-3]
   END IF
   nf(jz)=PNFL*1.D20  ! density from fast (non-thermal) ions [m^-3]

   zeff(jz)=ZEFFL     ! mean charge, Z_effective
   te(jz)=PTEL        ! T_e (electron temperature) [keV]
   ti(jz)=PTIL        ! T_i (temperature of thermal ions) [keV]
   q(jz)=QPL          ! magnetic q-value
   btor(jz)=BBL       ! ( R B_tor ) / rmaj(jz)  [Tesla]

   zimp(jz)=ZZL       ! mean charge of impurities
!                      = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
!                      sum_imp = sum over impurity ions with charge state Z_imp
   aimp(jz)=AMZL      ! mean atomic mass of impurities
!                      = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
!                      sum_imp = sum over impurity ions, each with mass M_imp
   ahyd(jz)=AMHL      ! mean atomic mass of hydrogen ions
!                      = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
!                      sum_hyd = sum over hydrogenic ions, each with mass M_hyd
   aimass(jz)=AMIL    ! mean atomic mass of thermal ions
!                      = ( sum_i n_i M_i ) / ( sum_i n_i ) where
!                      sum_i = sum over all ions, each with mass M_i
!
   wexbs(jz)=WEXBL    ! ExB shearing rate [rad/s].
!                       See K.H. Burrell, Phys. of Plasmas, 4, 1499 (1997).
!
!   All of the following normalized gradients are at radial points.
!   r = half-width, R = major radius to center of flux surface
!
   gne(jz)=-RR*DNEDR/PNEL  ! -R ( d n_e / d r ) / n_e
   gni(jz)=-RR*DNIDR/PNIL  ! -R ( d n_i / d r ) / n_i
   gnh(jz)=-RR*DNHDR/PNHL  ! -R ( d n_h / d r ) / n_h
   IF(PNZL.EQ.0.D0) THEN
      gnz(jz)=gnh(jz)
   ELSE
      gnz(jz)=-RR*DNZDR/PNZL ! -R ( d Z n_Z / d r ) / ( Z n_Z )
   END IF
   gte(jz)=-RR*DTEDR/PTEL  ! -R ( d T_e / d r ) / T_e
   gti(jz)=-RR*DTIDR/PTIL  ! -R ( d T_i / d r ) / T_i
   gq(jz)=RR/RSL*S(NR)     !  R ( d q   / d r ) / q

   gvtor(jz)=0.D0
   vtor(jz)=0.D0
   gvpol(jz)=0.D0
   vpol(jz)=0.D0
   gvpar(jz)=0.D0
   vpar(jz)=0.D0
!
! where:
!   n_i  = thermal ion density (sum over hydrogenic and impurity)
!   n_h  = thermal hydrogenic density (sum over hydrogenic species)
!   n_Z  = thermal impurity density,  Z = average impurity charge
!                     summed over all impurities

   cmodel(1)=1.D0   ! Weiland20
   cmodel(2)=1.D0   ! DRIBM
   cmodel(3)=0.D0   ! ETG (ETG requires multi-points)
   cswitch(1:MAXNOPT,1:4)=0.D0
    
   lprint=1  ! Verbose level
   nprout=21 ! Output unit number for long printout

call mmm7_1( &
   rmin   = rmin,   rmaj   = rmaj,   rmaj0  = rmaj(1),    &
   elong  = elong,  ne     = ne,     nh     = nh,         &
   nz     = nz,     nf     = nf,     zeff   = zeff,       &
   te     = te,     ti     = ti,     q      = q,          &
   btor   = btor,   zimp   = zimp,   aimp   = aimp,       &
   ahyd   = ahyd,   aimass = aimass, wexbs  = wexbs,      &
   gne    = gne,    gni    = gni,    gnh    = gnh,        &
   gnz    = gnz,    gte    = gte,    gti    = gti,        &
   gq     = gq,                                           &
   gvtor  = gvtor,  vtor   = vtor,   gvpol  = gvpol,      &
   vpol   = vpol,   gvpar  = gvpar,  vpar   = vpar,       &
   xti    = xti,    xdi    = xdi,    xte    = xte,        &
   xdz    = xdz,    xvt    = xvt,    xvp    = xvp,        &
   xtiW20 = xtiW20, xdiW20 = xdiW20, xteW20 = xteW20,     &
   xtiDBM = xtiDBM, xdiDBM = xdiDBM, xteDBM = xteDBM,     &
   xteETG = xteETG,                                       &
   gammaW20 = gammaW20, omegaW20 = omegaW20,              &
   gammaDBM = gammaDBM, omegaDBM = omegaDBM,              &
   npoints = npoints,                                     &
   lprint  = lprint, nprout  = nprout, nerr    = nerr,   &
   vconv   = vconv,  vflux   = vflux ,                    &
   cmodel  = cmodel, cswitch = cswitch, lswitch = lswitch)

   IERR=nerr

! Diffusivity profiles, which are the sum of all of the component
! diffusivities in the mmm7.1 model:
!
   CHII=xti(jz) ! Effective ion thermal diffusivity
   DIFH=xdi(jz) ! Effective hydrogenic ion diffusivity
   CHIE=xte(jz) ! Effective electron thermal diffusivity
   DIFZ=xdz(jz) ! Impurity ion diffusivity from the Weiland model
   VIST=xvt(jz) ! Toroidal momentum diffusivity from the Weiland model
   VISP=xvp(jz) ! Poloidal momentum diffusivity from the Weiland model

! The following component output arrays give the separate contribution from
! each internal model. Note that the momentum diffusivities are only provided
! by the Weiland model. Generally, these arrays are used for diagnostic
! output only.

!
! The component output profiles are optional. They should be used for the
! diagnostic output purpose.
!
   CHIIW=xtiW20(jz) ! Ion thermal diffusivity from the Weiland (W20) component
   DIFHW=xdiW20(jz) ! Hydrogenic ion particle diffusivity from the Weiland comp.
   CHIEW=xteW20(jz) ! Electron thermal diffusivity from the Weiland component
   CHIIB=xtiDBM(jz) ! Ion thermal diffusivity from the DRIBM component
   DIFHB=xdiDBM(jz) ! Hydrogenic ion diffusivity from the DRIBM component
   CHIEB=xteDBM(jz) ! Electron thermal diffusivity from the DRIBM component
   CHIEG=xteETG(jz) ! Electron thermal diffusivity from the Horton ETG component

! The following are growth rates and mode frequencies from the
! Weiland model for drift modes such as ITG and TEM.
! These arrays are intended for diagnostic output.
!
!   GAMMAL=gamma(jz) ! growth rate for the dominating mode at point jz ( 1/sec )
!   OMEGAL=omega(jz) ! frequency for dominating mode jm at point jz ( rad/sec )

!   vconv, &! Convective velocities [m/s]
!   vflux    ! Flux matrix

    IF(IERR.NE.0) THEN
       WRITE(6,'(A,I5)') 'XX mmm71_driver: IERR=',IERR
    END IF

    RETURN
  END SUBROUTINE mmm71_driver

END MODULE trmodels
