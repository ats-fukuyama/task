MODULE trcalv
!---------------------------------------------------------------------------
!
!   This module calculates variables to evaluate transport coefficients
!    and other effects. 
!
!   The variables declared below should be refered 
!    from ONLY 'trn/trmodels' directory.
!
!---------------------------------------------------------------------------

  USE trcomm, ONLY: rkind,ikind,nrmax,nsamax
  IMPLICIT NONE

  REAL(rkind),DIMENSION(:,:),ALLOCATABLE ::&
       rp,       &!the pressure of each species (nT)
       rp_d       !the deriv. of pressure of each species (dnT/dr) 
       
  REAL(rkind),DIMENSION(:),ALLOCATABLE ::&
       rp_tot,   &! the total pressure
       rp_totd,  &! the deriv. of total pressure
       rp_add,   &! the additional pressure
       rp_beam,  &! the beam pressure
!
       rt_e,     &! the electron temperature
       rt_em,    &! the electron temperature (half-mesh)
       rt_ed,    &! the deriv. of electron temperature
       rt_ecl,   &! the scale length of electron temperature 
       rt_i,     &! the effective hydrogenic ion temperature
       rt_im,    &! the effective hydrogenic ion temperature (half-mesh)
       rt_id,    &! the deriv. of effective hydrogenic ion temperature 
       rt_icl,   &! the scale length of hydrogenic ion temperature 
!
       rn_e,     &! the electron density
       rn_em,    &! the electron density (half-mesh)
       rn_ed,    &! the deriv. of electron density 
       rn_ecl,   &! the scale length of electron density 
       rn_i,     &! the sum of ion density
       rn_im,    &! the sum of ion density (half-mesh)
       rn_id,    &! the deriv. of ion density 
       rn_icl,   &! the scale length of ion density 
       qp_m,     &! safety factor (half-mesh)
       qp_d,     &! the deriv. of safety factor 
!       
       mshear,   &! magnetic shear            r/q * (dq/dr)
       mshear_cl,&! magnetic shear length  R*q**2/(r*dq/dr)
       mcurv,    &! magnetic curvature
       vexb,     &! ExBt velocity [m/s]
       dvexbdr,  &! ExBt velocity gradient [1/s]
       vexbp,    &! ExBp velocity [m/s]
       dvexbpdr, &! ExBp velocity gradient [1/s]
!       wexb,     &! ExBt shearing rate [rad/s]
       wexbp,    &! ExBp shearing rate [rad/s]
!       v_alf,    &! Alfven wave velocity
       v_se       ! speed of sound for electron

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       z_eff      ! Z_eff: effective charge

  ! for NCLASS output
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       chi_ncp,  &! coef. of the pres. grad. for matrix expression 
                  !  of heat flux
       chi_nct,  &! coef. of the temp. grad. for matrix expression 
                  !  of heat flux
       d_ncp,    &! coef. of the pres. grad. for matrix expression
                  !  of particle flux
       d_nct      ! coef. of the temp. grad. for martix expression
                  !  of particle flux
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       gfls,     &! radial particle flux components of s [/(m^2 s)]
       qfls       ! radial heat conduction flux components of s [W/m^2
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       fls_tot    ! total radial flux of s (heat and particle)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE::&
       vebs,     &! <E.B> particle convection velocity of s [m/s]
       qebs,     &! <E.B> heat convection velocity of s [m/s]
       dia_gdnc,  &! diagonal diffusivity [m^2/s]
       dia_gvnc,  &! diagonal convection driven by off-diagonal part [m/s]
       cjbs_p,   &! <J_bs.B> driven by unit p'/p of s [A*T/1.d-20*m^2]
       cjbs_t     ! <J_bs.B> driven by unit T'/T of s [A*T/1.d-20*m^2]

CONTAINS
!==========================================================================
  SUBROUTINE tr_calc_variables
!--------------------------------------------------------------------------
!   Calculate variables and some effects used in calculate transport 
!--------------------------------------------------------------------------
    USE trcomm, ONLY: ikind,nrmax,RR,rhog,rmnrho,BB,dpdrho,bp,er &
         ,nrd1,nrd2,nrd3,nrd4

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: deriv3 ! the 2nd order accuracy derivative in TASK/lib

    CALL tr_calv_fundamental

    CALL tr_er_field

    CALL tr_zeff

    vexbp(0:nrmax)    = 0.d0
    dvexbpdr(0:nrmax) = 0.d0
    wexbp(0:nrmax)    = 0.d0

    ! on grid; these variables are not defined at rho = 0
    DO nr = 1, nrmax
       ! ExB velocity : er [V/m] / BB [T] -> [m/s]
       vexb(nr)  = -er(nr) / BB
       vexbp(nr) = -er(nr) / (RR*bp(nr)) ! [1/s]
       
       !:the 2nd order accuracy derivative of V_exb
       dvexbdr(nr)  = deriv3(nr,rmnrho,vexb,nrmax,0)
       dvexbpdr(nr) = deriv3(nr,rmnrho,vexbp,nrmax,0)

       ! poloidal rotational shear [1/s]
!      ExB Rotation shear
!       "Effects of {ExB} velocity shear and magnetic shear
!           on turbulence and transport in magnetic confinement devices"
!       [Phys. of Plasmas, 4, 1499 (1997)]          
       wexbp(nr) = (RR*bp(nr))**2/(dpdrho(nr)*sqrt(BB**2+bp(nr)**2)) &
                   *dvexbpdr(nr)
       ! wexbp(nr) = RR*bp(nr)*dvexbpdr(nr)/BB
       
       ! Doppler shift
       ! AGMP(NR) = QP(NR)/EPS*WEXB(NR)

       ! sound speed for electron
       v_se  = 0.d0
       ! pressure gradient for MHD instability

    END DO

!    nrd1(0:nrmax) = er(0:nrmax)
!    nrd2(0:nrmax) = vexbp(0:nrmax)
!    nrd3(0:nrmax) = dvexbpdr(0:nrmax)
!    nrd4(0:nrmax) = bp(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc_variables

!==========================================================================

  SUBROUTINE tr_calv_fundamental
!-------------------------------------------------------------------------
!   The derivatives are calculated on the each grid.
!-------------------------------------------------------------------------
    USE trcomm, ONLY: rkev,pa,pz,pz0,idnsa,ar1rho,RR,BB,rhog,rmnrho, &
         rt,rn,bp,qp, nrd1,nrd2

    IMPLICIT NONE

    ! --- control and internal variables
    INTEGER(ikind) :: nr,ns,nsa
    REAL(rkind) :: rt_isum
    REAL(rkind) :: deriv3
    REAL(rkind),DIMENSION(0:nrmax) :: nr_array

!    *****
!    ! This part will be implemented in introducing the NBI heating modules.
!    *****
!    additional pressure due to NBI
!
!    IF(SUMPBM.EQ.0.D0) THEN
!       PADD(1:NRMAX)=0.D0
!    ELSE
!       ! density [10^20 m^-3] * temperature [keV]
!       PADD(1:NRMAX)=PBM(1:NRMAX)*1.D-20/RKEV-RNF(1:NRMAX,1)*RT(1:NRMAX,2)
!                                                             RTF(nr, 1)
!       ! padd = experimental - calculating result (beam pressure)
!    ENDIF


    ! initilization (zero at NR = 0) ----
    rn_i(0:nrmax)    = 0.d0
    rp_tot(0:nrmax)  = 0.d0
    rp_totd(0:nrmax) = 0.d0
    rp_add(0:nrmax)  = 0.d0
    rp_beam(0:nrmax) = 0.d0

    rp_d(1:nsamax,0:nrmax)   = 0.d0
    rn_ed(0:nrmax)  = 0.d0
    rn_ecl(0:nrmax) = 0.d0
    rn_id(0:nrmax)  = 0.d0
    rn_icl(0:nrmax) = 0.d0
    rt_ed(0:nrmax)  = 0.d0
    rt_ecl(0:nrmax) = 0.d0
    rt_id(0:nrmax)  = 0.d0
    rt_icl(0:nrmax) = 0.d0
    qp_d(0:nrmax)   = 0.d0
    mshear(0:nrmax) = 0.d0

    rt_e(0) = rt(1,0)
    rn_e(0) = rn(1,0)

    ! '_i' means 'hydrogennic ions'. 
    ! for now, in following, nsa are used as loop counters.
    rt_isum = 0.d0
    DO nsa = 1, nsamax
       rp(nsa,0) = rn(nsa,0)*1.d20 * rt(nsa,0)*rkev
       rp_tot(0) = rp_tot(0) + rp(nsa,0)       

       IF(idnsa(nsa) == 1) THEN ! ion
          rn_i(0) = rn_i(0) + rn(nsa,0)
          rt_isum = rt_isum + rn(nsa,0)*rt(nsa,0)
       END IF
    END DO
       rt_i(0) = rt_isum / rn_i(0)
    ! ---------------

    DO nr = 1, nrmax
       rt_isum = 0.d0
       
       rt_e(nr) = rt(1,nr)
       rn_e(nr) = rn(1,nr)
       DO nsa = 1, nsamax
          ! the pressure of each species
          rp(nsa,nr)  = rn(nsa,nr)*1.d20 * rt(nsa,nr)*rkev
          rp_tot(nr)  = rp_tot(nr) + rp(nsa,nr)

          IF(idnsa(nsa) == 1) THEN ! ion
             rn_i(nr) = rn_i(nr) + rn(nsa,nr)
             rt_isum  = rt_isum + rn(nsa,nr)*rt(nsa,nr)
          END IF            
       END DO         
       rt_i(nr)    = rt_isum / rn_i(nr)

       ! safety factor
       qp_m(nr)  = 0.5d0*(qp(nr)+qp(nr-1))
       qp_d(nr) = deriv3(nr,rmnrho,qp,nrmax,0)

       ! magnetic curvature
       mcurv = 0.d0
       
    END DO

    ! derivatives with respect to 'r' ( d XX/dr = <|grad rho|>*d XX/d rho)
    ! ** deriv3: second order accuracy on grid
    !          ; at nr=0, d XX/dr = 0 and 
    ! need consideration for PADD
    DO nr = 1, nrmax
       DO nsa = 1, nsamax
          nr_array(0:nrmax) = rp(nsa,0:nrmax)
          rp_d(nsa,nr) = deriv3(nr,rmnrho,nr_array,nrmax,0)
       END DO

       rn_ed(nr)   = deriv3(nr,rmnrho,rn_e,nrmax,0)
       rn_id(nr)   = deriv3(nr,rmnrho,rn_i,nrmax,0)
       rt_ed(nr)   = deriv3(nr,rmnrho,rt_e,nrmax,0)
       rt_id(nr)   = deriv3(nr,rmnrho,rt_i,nrmax,0)
       rp_totd(nr) = deriv3(nr,rmnrho,rp_tot,nrmax,0)

       ! scale length ( (d XX/d rho)/XX ) on grid
       rt_ecl(nr)  = rt_ed(nr) / rt_e(nr)
       rt_icl(nr)  = rt_id(nr) / rt_i(nr)
       rn_ecl(nr)  = rn_ed(nr) / rn_e(nr)
       rn_icl(nr)  = rn_id(nr) / rn_i(nr)

       ! mean atomic mass of thermal ions [AMU]

       ! impurity and hydrogen density, density weighted charge/atomic number

       ! magnetic shear
       mshear(nr) =  rmnrho(nr)/qp(nr) * qp_d(nr)
    END DO


  END SUBROUTINE tr_calv_fundamental


  SUBROUTINE tr_er_field
!---------------------------------------------------------------------------
!        Radial Electric Field (on grid)
!---------------------------------------------------------------------------
    USE trcomm, ONLY: nrmax,aee,bb,bp,rg,rt,rn,rg,rm,er, &
         pa,pz, &
         mdler,vtor,vpol

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAl(rkind) :: term_dp,term_dpg

    DO nr=0,nrmax

!    *****
!    ! This part will be implemented in introducing the NBI heating modules.
!    *****
!       DRL=RJCB(NR)/DR
!       IF(NR.EQ.NRMAX) THEN
!          ! second speices (main ion) only
!          DPD = DERIV3P(PNSS(2)*PTS(2), &
!                        RN(NR  ,2)*RT(NR  ,2)-PADD(NR  ), &
!                        RN(NR-1,2)*RT(NR-1,2)-PADD(NR-1), &
!                        RHOG(NR),RHOM(NR),RHOM(NR-1))
!          TERM_DP = DPD*RKEV/(PZ(2)*AEE*PNSS(2))
!       ELSE
       ! pressure gradient of bulk D(main ion)
!       DPD =(RN(NR+1,2)*RT(NR+1,2)-PADD(NR+1)-(RN(NR  ,2)*RT(NR  ,2)-PADD(NR  )))*DRL
!          TERM_DP = DPD*RKEV/(PZ(2)*AEE*0.5D0*(RN(NR+1,2)+RN(NR,2)))
!       ENDIF

! DPD -> rp_d(nsa,nr)

       ! second species (main ion: D(H)) only
       term_dp = rp_d(2,nr) / (pz(2)*aee*rn(2,nr)*1.d20)

       ! toroidal and poloidal rotation velocity <- from experiments for now
!       vtor = 0.d0
!       vpol = 0.d0

       IF(mdler.EQ.0) THEN
          ! pressure gradient only (nabla p)
          er(nr)  = term_dp
       ELSE IF(mdler.EQ.1) THEN
          ! nabla p + toroidal rotation (V_tor)
          er(nr) = term_dp+vtor(nr)*bp(nr)
       ELSE IF(mdler.EQ.2) THEN
          ! nabla p + V_tor + poloidal rotation (V_pol) *** typical ER ***
          er(nr) = term_dp+vtor(nr)*bp(nr)-vpol(nr)*bb
!       ELSEIF(MDLER.EQ.3) THEN
!          !     Waltz definition
!          EPS = EPSRHO(NR)
!          F_UNTRAP = 1.D0-1.46D0*SQRT(EPS)+0.46D0*EPS**1.5D0
!          ALPHA_NEO = 1.D0-0.8839D0*F_UNTRAP/(0.3477D0+0.4058D0*F_UNTRAP)
!          IF(NR.EQ.NRMAX) THEN
!             TEL = PTS(1)
!             TIL = PTS(2)
!             RLNI = -DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2),RHOG(NR),RHOM(NR),R\
!             HOM(NR-1))/PNSS(2)
!             RLTI = -DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RH\
!             OM(NR-1))/PTS(2)
!          ELSE
!             TEL = 0.5D0*(RT(NR,1)+RT(NR+1,1))
!             TIL = 0.5D0*(RT(NR,2)+RT(NR+1,2))
!             RLNI = -(LOG(RN(NR+1,2))-LOG(RN(NR,2)))*DRL
!             RLTI = -(LOG(RT(NR+1,2))-LOG(RT(NR,2)))*DRL
!          ENDIF
!          CS = SQRT(TEL*RKEV/(PA(2)*AMM))
!          RHO_S = CS*PA(2)*AMM/(PZ(2)*AEE*BB)
!          ER(NR) =-BB*( (TIL/TEL)*RHO_S*CS*(RLNI+ALPHA_NEO*RLTI)-EPS/QP(NR)*\
!          VTOR(NR))
       ENDIF

    ENDDO

    RETURN

  END SUBROUTINE tr_er_field


  SUBROUTINE tr_zeff
!---------------------------------------------------------------------------
!           Caluculate Z_eff (effective charge)
!---------------------------------------------------------------------------
!!$    USE TRCOMM, ONLY : ANC, ANFE, MDLEQN, MDLUF, NRMAX, PZ, PZC, PZFE, RN, R\
!!$    NF, RT, ZEFF
    USE trcomm, ONLY: nrmax,idnsa,ns_nsa,pz,rn,mdluf
    IMPLICIT NONE
    INTEGER(4):: nr,nsa,ns

    z_eff =0.d0

    ! not include impurities
    IF(mdluf == 0)THEN
       DO nr = 0, nrmax
          DO nsa = 1, nsamax
             IF(idnsa(nsa) == 1)THEN
                ns=ns_nsa(nsa)
                z_eff(nr) = z_eff(nr) + pz(ns)*rn(nsa,nr)
             END IF
          END DO
       END DO
    END IF

!!$    IF(MDLUF.EQ.0) THEN
!!$       DO NR=1,NRMAX
!!$          TE =RT(NR,1)
!!$          ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +PZ(\
!!$4         )  *PZ(4)  *RN(NR,4) &
!!$               &     &                +TRZEC(TE)**2   *ANC (NR) +TRZEFE(TE)**2  *ANFE(NR))/RN(\
!!$          NR,1)
!!$       ENDDO
!!$    ELSE
!!$       IF(MDLEQN.EQ.0) THEN ! fixed density
!!$          DO NR=1,NRMAX
!!$             ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +\
!!$             PZ(2)  *PZ(2)  *RNF(NR,1))/RN(NR,1)
!!$          ENDDO
!!$       ENDIF
!!$    ENDIF
!!$
!!$    DO NR=1,NRMAX
!!$       TE=RT(NR,1)
!!$       PZC(NR)=TRZEC(TE)
!!$       PZFE(NR)=TRZEFE(TE)
!!$    ENDDO

    RETURN
  END SUBROUTINE tr_zeff

!============================================================================
  SUBROUTINE tr_calv_nr_alloc

    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind),SAVE:: nsamax_save
    INTEGER(ikind)     :: ierr


    IF(nrmax /= nrmax_save .OR. &
       nsamax /= nsamax_save) THEN

       IF(nrmax_save /= 0) CALL tr_calv_nr_dealloc

       DO
          ALLOCATE(rp(nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(rp_d(nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          
          ALLOCATE(rp_tot(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(rp_totd(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT  
          ALLOCATE(rp_add(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(rp_beam(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT  
          
          ALLOCATE(rt_e(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rt_em(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rt_ed(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rt_ecl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(rt_i(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rt_im(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rt_id(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rt_icl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          
          ALLOCATE(rn_e(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rn_em(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rn_ed(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rn_ecl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rn_i(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rn_im(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rn_id(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(rn_icl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT        
          ALLOCATE(qp_m(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT      
          ALLOCATE(qp_d(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT      
          
          ALLOCATE(mshear(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(mcurv(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(vexb(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(dvexbdr(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT 
          ALLOCATE(vexbp(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(dvexbpdr(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT 
          ALLOCATE(wexbp(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT       
!          ALLOCATE(v_alf(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(v_se(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ALLOCATE(z_eff(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ! for NCLASS
          ALLOCATE(chi_ncp(1:nsamax,1:nsamax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(chi_nct(1:nsamax,1:nsamax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(d_ncp(1:nsamax,1:nsamax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(d_nct(1:nsamax,1:nsamax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(gfls(5,1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(qfls(5,1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(fls_tot(3,1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(vebs(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(qebs(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(dia_gdnc(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(dia_gvnc(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(cjbs_p(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(cjbs_t(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          nrmax_save  = nrmax
          nsamax_save = nsamax
          RETURN
       END DO

       WRITE(6,*) 'XX tr_calv_nr_alloc: allocation error: ierr=',ierr
       STOP
   END IF

   RETURN
  END SUBROUTINE tr_calv_nr_alloc

  SUBROUTINE tr_calv_nr_dealloc

    IF(ALLOCATED(rp)) DEALLOCATE(rp)
    IF(ALLOCATED(rp_d)) DEALLOCATE(rp_d)
!    
    IF(ALLOCATED(rp_tot)) DEALLOCATE(rp_tot)
    IF(ALLOCATED(rp_totd)) DEALLOCATE(rp_totd)
    IF(ALLOCATED(rp_add)) DEALLOCATE(rp_add)
    IF(ALLOCATED(rp_beam)) DEALLOCATE(rp_beam)
!
    IF(ALLOCATED(rt_e)) DEALLOCATE(rt_e)
    IF(ALLOCATED(rt_em)) DEALLOCATE(rt_em)
    IF(ALLOCATED(rt_ed)) DEALLOCATE(rt_ed)
    IF(ALLOCATED(rt_ecl)) DEALLOCATE(rt_ecl)
    IF(ALLOCATED(rt_i)) DEALLOCATE(rt_i)
    IF(ALLOCATED(rt_im)) DEALLOCATE(rt_im)
    IF(ALLOCATED(rt_id)) DEALLOCATE(rt_id)
    IF(ALLOCATED(rt_icl)) DEALLOCATE(rt_icl)
!
    IF(ALLOCATED(rn_e)) DEALLOCATE(rn_e)
    IF(ALLOCATED(rn_em)) DEALLOCATE(rn_em)
    IF(ALLOCATED(rn_ed)) DEALLOCATE(rn_ed)
    IF(ALLOCATED(rn_ecl)) DEALLOCATE(rn_ecl)
    IF(ALLOCATED(rn_i)) DEALLOCATE(rn_i)
    IF(ALLOCATED(rn_im)) DEALLOCATE(rn_im)
    IF(ALLOCATED(rn_id)) DEALLOCATE(rn_id)
    IF(ALLOCATED(rn_icl)) DEALLOCATE(rn_icl)
    IF(ALLOCATED(qp_m)) DEALLOCATE(qp_m)
    IF(ALLOCATED(qp_d)) DEALLOCATE(qp_d)
!       
    IF(ALLOCATED(mshear)) DEALLOCATE(mshear)
    IF(ALLOCATED(mcurv)) DEALLOCATE(mcurv)
    IF(ALLOCATED(vexb)) DEALLOCATE(vexb)
    IF(ALLOCATED(dvexbdr)) DEALLOCATE(dvexbdr)
    IF(ALLOCATED(vexbp)) DEALLOCATE(vexbp)
    IF(ALLOCATED(dvexbpdr)) DEALLOCATE(dvexbpdr)
    IF(ALLOCATED(wexbp)) DEALLOCATE(wexbp)
!    IF(ALLOCATED(v_alf)) DEALLOCATE(v_alf)
    IF(ALLOCATED(v_se)) DEALLOCATE(v_se)

    IF(ALLOCATED(z_eff)) DEALLOCATE(z_eff)

    ! for NCLASS
    IF(ALLOCATED(chi_ncp)) DEALLOCATE(chi_ncp)
    IF(ALLOCATED(chi_nct)) DEALLOCATE(chi_nct)
    IF(ALLOCATED(d_ncp)) DEALLOCATE(d_ncp)
    IF(ALLOCATED(d_nct)) DEALLOCATE(d_nct)
    IF(ALLOCATED(gfls)) DEALLOCATE(gfls)
    IF(ALLOCATED(qfls)) DEALLOCATE(qfls)
    IF(ALLOCATED(fls_tot)) DEALLOCATE(fls_tot)
    IF(ALLOCATED(vebs)) DEALLOCATE(vebs)
    IF(ALLOCATED(qebs)) DEALLOCATE(qebs)
    IF(ALLOCATED(dia_gdnc)) DEALLOCATE(dia_gdnc)
    IF(ALLOCATED(dia_gvnc)) DEALLOCATE(dia_gvnc)
    IF(ALLOCATED(cjbs_p)) DEALLOCATE(cjbs_p)
    IF(ALLOCATED(cjbs_t)) DEALLOCATE(cjbs_t)

    RETURN
  END SUBROUTINE tr_calv_nr_dealloc

END MODULE trcalv

