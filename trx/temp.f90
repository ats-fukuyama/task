   INTEGER,PARAMETER:: npoints ! Number of values in all of the 1-D arrays
   INTEGER,PARAMETER:: MMM_NCH=6  ! Maximum number of transport channels
   Real(8), Dimension(npoints) :: &
        rmin, rmaj, elong, ne, nh, nz, nf, zeff, &
        te, ti, q, btor, zimp, aimp, ahyd, aimass, wexbs, &
        gne, gni, gnh, gnz, gte, gti, gq
   Real(8), Dimension(npoints) :: & ![m^2/s]
        xti, xdi, xte, xdz, xvt, xvp
   Real(8), Dimension(npoints) :: & ![m^2/s]
        xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG
   Real(R8), Dimension(MMM_NCH,npoins) :: &
        vconv, vflux
   Integer:: &
        lprint,nprout,nerr

    IF(NR.EQ.NRMAX) THEN
       DO NS=1,NSMAX
          RNL(NS)=PNSS(NS)
          RNTL(NS)=PNSS(NS)*PTS(NS)
          DRNL(NS)=(PNSS(NS)-RN(NRMAX,NS))/(RA-RM(NRMAX))
          DRNTL(NS)=(PNSS(NS)*PTS(NS) &
                    -RN(NRMAX,NS)*RT(NRMAX,NS))/(RA-RM(NRMAX))
       END DO
       DO NF=1,NFM
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
       DO NF=1,NFM
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
       IF(PA(NS).LE.0.01D0) THEN
          PNEL=PNEL+RNL(NS)
          PNTEL=PNTEL+RNTL(NS)
          DNEDR=DNEDR+DRNL(NS)
          DTEDR=DTEDR+DRNTL(NS)
       ELSE
          IF(PA(NS).LE.3.D0.AND.PZ(NS).EQ.1.D0) THEN  
             ! 3He+ and 3H+ are not distinguished; waiting for the use of PZ0
             PNHL=PNHL+RNL(NS)
             AMHL=AMHL+PA(NS)*RNL(NS)
             DNHDR=DNHDR+DRNL(NS)
          ELSE
             PNZL=PNZL+RNL(NS)
             ZZL =ZZL +PZ(NS)*RNL(NS)
             AMZL=AMZL+PA(NS)*RNL(NS)
             DNZDR=DNZDR+DRNL(NS)
          END IF
          ZEFFL=ZEFFL+PZ(NS)*PZ(NS)*RNL(NS)
          ZNIL =ZNIL +PZ(NS)*RNL(NS)
          PNIL =PNIL +RNL(NS)
          PNTIL=PNTIL+RNTL(NS)
          AMIL =AMIL +PA(NS)*RNL(NS)
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

   jz=1 rmin(jz)=RSL ! half-width of the magnetic surface [m]
   rmaj(jz)=RR ! major radius to geometric center of the magnetic
   surface [m] elong(jz)=RKPRHO(NR) ! local elongation

   ne(jz)=PNEL*1.D20  ! electron density [m^-3]
   nh(jz)=PNHL*1,D20  ! hydrogenic thermal particle density [m^-3]
   nz(jz)=PNZL*1,D20  ! impurity ion density [m^-3]
   nf(jz)=PNFL*1,D20  ! density from fast (non-thermal) ions [m^-3]

   zeff(jz)=ZEFFL     ! mean charge, Z_effective
   te(jz)=PTEL        ! T_e (electron temperature) [keV]
   ti(jz)=PTIL        ! T_i (temperature of thermal ions) [keV]
   q(jz)=QPL          ! magnetic q-value
   btor(jz)BBL        ! ( R B_tor ) / rmaj(jz)  [Tesla]

   zimp(jz)=ZZL       ! mean charge of impurities
!                       = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
!                       sum_imp = sum over impurity ions with charge state Z_imp
   aimp(jz)=AMZL      ! mean atomic mass of impurities
!                       = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
!                       sum_imp = sum over impurity ions, each with mass M_imp
   ahyd(jz)=AMHL      ! mean atomic mass of hydrogen ions
!                       = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
!                       sum_hyd = sum over hydrogenic ions, each with mass M_hyd
   aimass(jz)=AMIL    ! mean atomic mass of thermal ions
!                       = ( sum_i n_i M_i ) / ( sum_i n_i ) where
!                       sum_i = sum over all ions, each with mass M_i
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
!
! where:
!   n_i  = thermal ion density (sum over hydrogenic and impurity)
!   n_h  = thermal hydrogenic density (sum over hydrogenic species)
!   n_Z  = thermal impurity density,  Z = average impurity charge
!                     summed over all impurities

   lprint=0  ! Verbose level
   nprout=21 ! Output unit number for long printout

   call mmm7_1( &
   rmin   = rmin,   rmaj   = rmaj,   & !rmaj0  = rmaj(1),    &
   elong  = elong,  ne     = ne,     nh     = nh,         &
   nz     = nz,     nf     = nf,     zeff   = zeff,       &
   te     = te,     ti     = ti,     q      = q,          &
   btor   = btor,   zimp   = zimp,   aimp   = aimp,       &
   ahyd   = ahyd,   aimass = aimass, wexbs  = wexbs,      &
   gne    = gne,    gni    = gni,    gnh    = gnh,        &
   gnz    = gnz,    gte    = gte,    gti    = gti,        &
   gq     = gq,                                           &
!   gvtor  = gvtor,  vtor   = vtor,   gvpol  = gvpol,      &
!   vpol   = vpol,   gvpar  = gvpar,  vpar   = vpar,       &
   xti    = xti,    xdi    = xdi,    xte    = xte,        &
   xdz    = xdz,    xvt    = xvt,    xvp    = xvp,        &
   xtiW20 = xtiW20, xdiW20 = xdiW20, xteW20 = xteW20,     &
   xtiDBM = xtiDBM, xdiDBM = xdiDBM, xteDBM = xteDBM,     &
   xteETG = xteETG,                                       &
   gammaW20 = gammaW20, omegaW20 = omegaW20,              &
   gammaDBM = gammaDBM, omegaDBM = omegaDBM,              &
   npoints = npoints,                                     &
   lprint  = lprint, nprout  = nprout,                    &
   vconv   = vconv,  vflux   = vflux ,                    &
!   cmodel  = cmodel, cswitch = cswitch, lswitch = lswitch, &
   nerr = nerr)

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


