!     $Id$

MODULE trmodels
  PRIVATE
  PUBLIC mbgb_driver

CONTAINS

  SUBROUTINE mbgb_driver(RSL,RRL,ANEL,TEL,TIL,DTEL,DNEL,PAL,PZL, &
                         QPL,BBL,WEXBL,SL, &
                         ADFFI,ACHIE,ACHII,ACHIEB,ACHIIB,ACHIEGB,ACHIIGB, &
                         ierr)
    USE trcomm, ONLY: &
         RT,NRMAX,MDLKAI
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
    REAL(8):: zte_p8,zte_edge
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
    IF(MDLKAI.EQ.140) THEN
       lflowshear=0
    ELSE
       lflowshear=1
    ENDIF

    call mixed_model ( &
         rminor,  rmajor,  tekev,   tikev,   q, &
         btor,    aimass,  charge,  wexbs, &
         grdte,   grdne,   shear, &
         zte_p8, zte_edge, npoints, &
         chi_i_mix,  themix,   thdmix, &
         thigb,   thegb,    thibohm, thebohm, &
         ierr, lflowshear)
           
    ADFFI=thdmix(1)
    ACHIE=themix(1)
    ACHII=chi_i_mix(1)
    ACHIEB=thebohm(1)
    ACHIIB=thibohm(1)
    ACHIEGB=thegb(1)
    ACHIIGB=thigb(1)
    RETURN
  END SUBROUTINE mbgb_driver

END MODULE trmodels
