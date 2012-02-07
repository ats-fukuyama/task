MODULE trcalv
!---------------------------------------------------------------------------
!
!   This routine calculates variables to evaluate transport coefficients
!    and other effects. 
!
!   The variables declared below should be refered 
!    from only 'trn/trmodels' directory.
!
!---------------------------------------------------------------------------

  USE trcomm, ONLY: rkind,ikind,nrmax,nsamax
  IMPLICIT NONE

  REAL(rkind),DIMENSION(:,:),ALLOCATABLE ::&
       rp,       &!the pressure of each species (nT)
       rp_d       !the deriv. of pressure of each species (dnT/dr) (half-mesh)
       
  REAL(rkind),DIMENSION(:),ALLOCATABLE ::&
       rp_tot,   &! the total pressure
       rp_totd,  &! the deriv. of total pressure
       rp_add,   &! the additional pressure
       rp_beam,  &! the beam pressure
!
       rt_e,     &! the electron temperature
       rt_em,    &! the electron temperature (half-mesh)
       rt_ed,    &! the deriv. of electron temperature (half-mesh)
       rt_ecl,   &! the scale length of electron temperature (half-mesh)
       rt_i,     &! the effective ion temperature
       rt_im,    &! the effective ion temperature (half-mesh)
       rt_id,    &! the deriv. of effective ion temperature (half-mesh)
       rt_icl,   &! the scale length of ion temperature (half-mesh)
!
       rn_e,     &! the electron density
       rn_em,    &! the electron density (half-mesh)
       rn_ed,    &! the deriv. of electron density (half-mesh)
       rn_ecl,   &! the scale length of electron density (half-mesh)
       rn_i,     &! the sum of ion density
       rn_im,    &! the sum of ion density (half-mesh)
       rn_id,    &! the deriv. of ion density (half-mesh)
       qp_m,     &! safety factor (half-mesh)
       qp_d,     &! the deriv. of safety factor (half-mesh)
!       
       mshear,   &! magnetic shear            r/q * (dq/dr)
       mshear_cl,&! magnetic shear length  R*q**2/(r*dq/dr)
       mcurv,    &! magnetic curvature
       vexb,     &! ExB velocity [m/s]
       vexbp,    &! poloidal ExB velocity [m/s]
       dvexbpdr, &! poloidal ExB velocity gradient [1/s]
       wexbp,    &! ExB shearing rate [rad/s]
!       v_alf,    &! Alfven wave velocity
       v_se       ! speed of sound for electron

CONTAINS

  SUBROUTINE tr_calv_nr_alloc

    INTEGER(ikind),SAVE:: nrmax_save
    INTEGER(ikind),SAVE:: nsamax_save
    INTEGER(ikind)     :: ierr


    IF(nrmax /= nrmax_save .OR. &
       nsamax /= nsamax_save) THEN

       IF(nrmax_save /= 0) CALL tr_calv_nr_dealloc

       ALLOCATE(rp(nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000
       ALLOCATE(rp_d(nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000
       
       ALLOCATE(rp_tot(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000   
       ALLOCATE(rp_totd(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000  
       ALLOCATE(rp_add(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000   
       ALLOCATE(rp_beam(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000  

       ALLOCATE(rt_e(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rt_em(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rt_ed(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(rt_ecl(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000   
       ALLOCATE(rt_i(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rt_im(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rt_id(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(rt_icl(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000   

       ALLOCATE(rn_e(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rn_em(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rn_ed(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(rn_ecl(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(rn_i(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rn_im(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(rn_id(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(qp_m(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000      
       ALLOCATE(qp_d(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000      
       
       ALLOCATE(mshear(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000   
       ALLOCATE(mcurv(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(vexb(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000     
       ALLOCATE(vexbp(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(dvexbpdr(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000 
       ALLOCATE(wexbp(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000       
!       ALLOCATE(v_alf(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000    
       ALLOCATE(v_se(0:nrmax),STAT=ierr); IF(ierr /=0 ) GOTO 9000

       nrmax_save  = nrmax
       nsamax_save = nsamax
   END IF
   RETURN
9000 WRITE(6,*) 'XX tr_calv_nr_alloc: allocation error: ierr=',ierr
   STOP
  END SUBROUTINE tr_calv_nr_alloc

  SUBROUTINE tr_calv_nr_dealloc

    IF(ALLOCATED(rp)) DEALLOCATE(rp)
    IF(ALLOCATED(rp_d)) DEALLOCATE(rp_d)
    
    IF(ALLOCATED(rp_tot)) DEALLOCATE(rp_tot)
    IF(ALLOCATED(rp_totd)) DEALLOCATE(rp_totd)
    IF(ALLOCATED(rp_add)) DEALLOCATE(rp_add)
    IF(ALLOCATED(rp_beam)) DEALLOCATE(rp_beam)
!
    IF(ALLOCATED(rt_e)) DEALLOCATE(rt_e)
    IF(ALLOCATED(rt_em)) DEALLOCATE(rt_em)
    IF(ALLOCATED(rt_ed)) DEALLOCATE(rt_ed)
    IF(ALLOCATED(rt_ecl)) DEALLOCATE(rt_icl)
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
    IF(ALLOCATED(qp_m)) DEALLOCATE(qp_m)
    IF(ALLOCATED(qp_d)) DEALLOCATE(qp_d)
!       
    IF(ALLOCATED(mshear)) DEALLOCATE(mshear)
    IF(ALLOCATED(mcurv)) DEALLOCATE(mcurv)
    IF(ALLOCATED(vexb)) DEALLOCATE(vexb)
    IF(ALLOCATED(vexbp)) DEALLOCATE(vexbp)
    IF(ALLOCATED(dvexbpdr)) DEALLOCATE(dvexbpdr)
    IF(ALLOCATED(wexbp)) DEALLOCATE(wexbp)
!    IF(ALLOCATED(v_alf)) DEALLOCATE(v_alf)
    IF(ALLOCATED(v_se)) DEALLOCATE(v_se)

    RETURN
  END SUBROUTINE tr_calv_nr_dealloc

!==========================================================================

  SUBROUTINE tr_calc_variables
    USE trcomm, ONLY: rkev,pa,pz,pz0,idnsa,ar1rho,RR,BB,rg,rhog,rm, &
         rt,rn,bp,qp,er

    IMPLICIT NONE

    ! --- control and internal variables
    INTEGER(ikind) :: nr,ns,nsa
    REAL(rkind) :: rt_isum, dr_norm
    REAL(rkind) :: deriv3 ! the 2nd order accuracy derivative in TASK/lib

    rp_add  = 0.d0
    rp_beam = 0.d0
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

    
    rp_tot  = 0.d0
    rp_totd = 0.d0
    rt_isum = 0.d0

    ! ---  NR = 0 ----
       rt_e(0) = rt(1,0)
       rn_e(0) = rn(1,0)

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
       dr_norm = ar1rho(nr)/(rhog(nr)-rhog(nr-1))
       
       rt_e(nr) = rt(1,nr)
       rn_e(nr) = rn(1,nr)
       DO nsa = 1, nsamax
          ! the pressure of each species
          rp(nsa,nr)  = rn(nsa,nr)*1.d20 * rt(nsa,nr)*rkev
          rp_d(nsa,nr) = (rp(nsa,nr)-rp(nsa,nr-1)) * dr_norm

          rp_tot(nr)  = rp_tot(nr) + rp(nsa,nr)

          IF(idnsa(nsa) == 1) THEN ! ion
             rn_i(nr) = rn_i(nr) + rn(nsa,nr)
             rt_isum  = rt_isum + rn(nsa,0)*rt(nsa,0)
          END IF            
       END DO         
       rt_i(nr)    = rt_isum / rn_i(nr)

       ! half-mesh
       rn_em(nr)   = 0.5d0*(rn_e(nr)+rn_e(nr-1))
       rn_im(nr)   = 0.5d0*(rn_i(nr)+rn_i(nr-1))          
       rt_em(nr)   = 0.5d0*(rt_e(nr)+rt_e(nr-1))
       rt_im(nr)   = 0.5d0*(rt_i(nr)+rt_i(nr-1))

       ! derivatives
       rn_ed(nr)   = (rn_e(nr)-rn_e(nr-1)) * dr_norm
       rn_id(nr)   = (rn_i(nr)-rn_i(nr-1)) * dr_norm
       rt_ed(nr)   = (rt_e(nr)-rt_e(nr-1)) * dr_norm
       rt_id(nr)   = (rt_i(nr)-rt_i(nr-1)) * dr_norm
       rp_totd(nr) = (rp_tot(nr)-rp_tot(nr-1)) * dr_norm

       ! scale length
       rt_icl(nr)  = rt_im(nr) / rt_id(nr)
       rt_ecl(nr)  = rt_em(nr) / rt_ed(nr)
       rn_ecl(nr)  = rn_em(nr) / rn_ed(nr)

       ! safety factor
       qp_m(nr) = 0.5d0*(qp(nr)+qp(nr-1))
       qp_d(nr) = (qp(nr)-qp(nr-1)) * dr_norm

       ! magnetic shear
       mshear(nr) =  rm(nr)/qp_m(nr) * qp_d(nr)

       ! magnetic curvature
       mcurv = 0.d0
       
       ! ExB velocity : er [V/m] / BB [T] -> [m/s]
       CALL tr_er_field
       vexb(nr)  = -er(nr) / BB
       vexbp(nr) = -er(nr) / (RR*bp(nr))

       !:the 2nd order accuracy derivative of V_exb
       dvexbpdr(nr) = deriv3(nr,rg,vexbp,nrmax,1)

       ! rotational shear (poloidal) ??? the dimension is NOT [rad/s]
       wexbp(nr) = RR*bp(nr)*dvexbpdr(nr)/BB
    END DO

    ! Doppler shift
    ! AGMP(NR) = QP(NR)/EPS*WEXB(NR)

    ! sound speed for electron
    v_se  = 0.d0
    ! pressure gradient for MHD instability

  END SUBROUTINE tr_calc_variables


!---------------------------------------------------------------------------
!        Radial Electric Field
!---------------------------------------------------------------------------

  SUBROUTINE tr_er_field
    USE trcomm, ONLY: nrmax,aee,bb,bp,rg,rt,rn,rg,rm,er, &
         pa,pz, &
         mdler,vtor,vpol

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAl(rkind) :: term_dp


    DO nr=1,nrmax

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

       ! second species (main ion: D) only  for now (half-mesh)
       term_dp = rp_d(2,nr) / (pz(2)*aee*0.5d0*(rn(2,nr)+rn(2,nr-1))*1.d20)

       ! toroidal and poloidal rotation velocity <- from experiments for now
       vtor = 0.d0
       vpol = 0.d0

       IF(mdler.EQ.0) THEN
          ! pressure gradient only (nabla p)
          er(nr) = term_dp
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

END MODULE trcalv

