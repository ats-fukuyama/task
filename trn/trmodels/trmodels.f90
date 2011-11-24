MODULE trmodels
  PRIVATE
  PUBLIC mbgb_driver

CONTAINS

  SUBROUTINE mbgb_driver(RSL,RRL,ANEL,TEL,TIL,DTEL,DNEL,PAL,PZL, &
                         QPL,BBL,WEXBL,SL, &
                         ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB, &
                         ierr)
    USE trcomm, ONLY: &
         RT,NRMAX,mdld,aee,amp
    USE mixed_Bohm_gyro_Bohm, ONLY: mixed_model
    IMPLICIT NONE
    REAL(8),INTENT(IN) :: &
         RSL,RRL,ANEL,TEL,TIL,DTEL,DNEL,PAL,PZL,QPL,BBL,WEXBL,SL
    REAL(8),INTENT(OUT) :: ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB
    INTEGER,INTENT(OUT):: ierr
    REAL(8),DIMENSION(1):: &
            rminor,rmajor,tekev,tikev,q,btor,aimass,charge,wexbs, &
            grdte,grdne,shear, &
            chi_i_mix,themix,thdmix,thigb,thegb,thibohm,thebohm
    REAL(8):: zte_p8,zte_edge,VTI,GAMMA0,EXBfactor,SHRfactor,factor
    INTEGER:: NR8,NPOINTS,lflowshear

    REAL(8):: rkev,amm
    integer:: mdlkai

    rkev=aee*1.D-3
    amm=amp
    mdlkai=mdld

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
       VTI=SQRT(2.D0*TIL*RKEV/(PAL*AMM))
       GAMMA0=VTI/(QPL*RRL)
       SHRfactor=1.D0
       EXBfactor=1.D0/(1.D0+(WEXBL/GAMMA0)**2)
    CASE(143)
       lflowshear=0
       VTI=SQRT(2.D0*TIL*RKEV/(PAL*AMM))
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

END MODULE trmodels
