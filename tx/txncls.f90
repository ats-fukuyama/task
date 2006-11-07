module nclass_mod
  use commons
  implicit none
  private
  public :: TX_NCLASS

contains

!***********************************************************
!
!       NCLASS
!
!***********************************************************

  SUBROUTINE TX_NCLASS(NR,NueNC,NuiNC,ETAL,JBSL,IER)
!***********************************************************************
!TX_NCLASS calculates various paramters and arrays for NCLASS.
!Please note that type declarations of all variables except "INTEGER" 
!  in NCLASS subroutine are "REAL(*4)" or "SINGLE" but not "REAL*8" 
!  or "DOUBLE".
!Input:
!  k_order - order of v moments to be solved (-)
!          =2 u and q
!          =3 u, q, and u2
!          =else error
!  k_potato - option to include potato orbits (-)
!           =0 off
!           =else on
!  m_i - number of isotopes (1<m_i<mx_mi+1)
!  m_z - highest charge state of all species (0<m_z<mx_mz+1)
!  c_den - density cutoff below which species is ignored (/m**3)
!  c_potb - kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl - q(0)*R(0) (m)
!  p_b2 - <B**2> (T**2)
!  p_bm2 - <1/B**2> (/T**2)
!  p_eb - <E.B> (V*T/m)
!  p_fhat - mu_0*F/(dPsi/dr) (rho/m)
!  p_fm(3) - poloidal moments of geometric factor for PS viscosity (-)
!  p_ft - trapped particle fraction (-)
!  p_grbm2 - <grad(rho)**2/BB**2> (rho**2/m**2/T**2)
!  p_grphi - radial electric field Phi' (V/rho)
!  p_gr2phi - radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_ngrth - <n.grad(Theta)> (1/m)
!  amu_i(i) - atomic mass number of i (-)
!  grt_i(i) - temperature gradient of i (keV/rho)
!  temp_i(i) - temperature of i (keV)
!  den_iz(i,z) - density of i,z (/m**3)
!  fex_iz(3,i,z) - moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z) - pressure gradient of i,z (keV/m**3/rho)
!Using Output:
!  p_bsjb - <J_bs.B> (A*T/m**2)
!  p_etap - parallel electrical resistivity (Ohm*m)
!  p_exjb - <J_ex.B> current response to fex_iz (A*T/m**2)
!  bsjbp_s(s) - <J_bs.B> driven by unit p'/p of s (A*T*rho/m**3)
!  bsjbt_s(s) - <J_bs.B> driven by unit T'/T of s (A*T*rho/m**3)
!  gfl_s(m,s) - radial particle flux comps of s (rho/m**3/s)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  qfl_s(m,s) - radial heat conduction flux comps of s (W*rho/m**3)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  veb_s(s) - <E.B> particle convection velocity of s (rho/s)
!  qeb_s(s) - <E.B> heat convection velocity of s (rho/s)
!  chip_ss(s1,s2) - heat cond coefficient of s2 on p'/p of s1 (rho**2/s)
!  chit_ss(s1,s2) - heat cond coefficient of s2 on T'/T of s1 (rho**2/s)
!  dp_ss(s1,s2) - diffusion coefficient of s2 on p'/p of s1 (rho**2/s)
!  dt_ss(s1,s2) - diffusion coefficient of s2 on T'/T of s1 (rho**2/s)
!  iflag - warning and error flag
!        =-4 warning: no viscosity
!        =-3 warning: no banana viscosity
!        =-2 warning: no Pfirsch-Schluter viscosity
!        =-1 warning: no potato orbit viscosity
!        = 0 no warnings or errors
!        = 1 error: order of v moments to be solved must be 2 or 3
!        = 2 error: number of species must be 1<m_i<mx_mi+1
!        = 3 error: number of species must be 0<m_z<mx_mz+1
!        = 4 error: number of species must be 1<m_s<mx_ms+1
!        = 5 error: inversion of flow matrix failed
!        = 6 error: trapped fraction must be 0.0.le.p_ft.le.1.0
!***********************************************************************
    use physical_constants, only : AME
    INCLUDE 'nclass/pamx_mi.inc'
    INCLUDE 'nclass/pamx_ms.inc'
    INCLUDE 'nclass/pamx_mz.inc'
    INCLUDE 'txncls.inc'
    INTEGER, INTENT(IN)  :: NR
    INTEGER, INTENT(OUT) :: IER
    REAL(8), INTENT(OUT) :: NueNC, NuiNC, ETAL, JBSL
    INTEGER :: i, k_out, k_v, ier_check, im, iz
    REAL(4) :: a0, bt0, e0, p_eps, p_q, q0l, r0
    REAL(8) :: EpsL, BBL, PZMAX, p_gr2phi1, p_gr2phi2, p_gr2phi3, &
         &                       p_fhat1, p_fhat2, p_fhat3, &
         &     btot, uthai, VPOL(0:NRMAX)
    REAL(8) :: DERIV3, AITKEN2P

    !     *** Initialization ***

    amu_i (1:mx_mi) = 0.0
    grt_i (1:mx_mi) = 0.0
    temp_i(1:mx_mi) = 0.0
    den_iz(1:mx_mi,1:mx_mz) = 0.0
    grp_iz(1:mx_mi,1:mx_mz) = 0.0
    fex_iz(1:3,1:mx_mi,1:mx_mz) = 0.0

    IER = 0

    !  k_out-option for output to nout (-)
    !       = 0 no display
    !        +1 warnings
    !        +2 errors
    !        +4 results
    !  k_v-option for neoclassical v_tor,v_pol,v_para,v_perp
    !       =1 output
    !       =else no output
    k_out    = 2
    k_v      = 1
    k_order  = 2
    k_potato = 1
    m_i      = 2
    PZMAX    = PZ
    m_z      = INT(PZMAX)
    c_den    = 1.E10
    c_potb   = SNGL(1.D0*BB/(2.D0*Q0**2))
    c_potl   = SNGL(Q0*RR)

    amu_i(1) = SNGL(AME/AMI)
    amu_i(2) = SNGL(PA)

    BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
    EpsL = R(NR) / RR
    p_b2   = SNGL(BphV(NR)**2)
    p_bm2  = SNGL(1.D0 / BphV(NR)**2)
    p_eb   = SNGL(EphV(NR)*BphV(NR)+EthV(NR)*BthV(NR)) 
    IF(NR == 0) THEN
       p_fhat1 = BphV(NR+1)/(RA*BthV(NR+1))
       p_fhat2 = BphV(NR+2)/(RA*BthV(NR+2))
       p_fhat3 = BphV(NR+3)/(RA*BthV(NR+3))
       p_fhat  = SNGL(AITKEN2P(R(0),p_fhat1,p_fhat2,p_fhat3,R(1),R(2),R(3)))
    ELSE
       !        p_fhat = SNGL(2.D0*PI*BphV(NR)/(RA*BthV(NR)))
       p_fhat = SNGL(BphV(NR)/(RA*BthV(NR)))
    END IF

    !  Approximation inverse aspect ratio at the magnetix axis
    IF(NR == 0) EpsL = 0.1D0*R(NR+1)/RR
    p_fm(1:3) = 0.0
    IF(SNGL(EpsL) > 0.0) THEN
       DO i=1,3
          p_fm(i)=SNGL(DBLE(i)*((1.D0-SQRT(1.D0-EpsL**2))/EpsL)**(2*i) &
               &                 *(1.D0+DBLE(i)*SQRT(1.D0-EpsL**2))/((1.D0-EpsL**2) &
               &                 **1.5D0*(Q(NR)*RR)**2))
       ENDDO
    ENDIF
    p_ft=SNGL(1.46D0 * SQRT(EpsL) - 0.46 * EpsL * SQRT(EpsL))

    p_grbm2=SNGL(1.D0/RA) * p_bm2
    p_grphi=SNGL(ErV(NR))
    IF(NR == 0) THEN
       p_gr2phi1 = -RA**2*DERIV3(NR+1,R,ErV,NRMAX,NRM,0) &
            &      +RA*ErV(NR+1)*DERIV3(NR+1,R,BthV,NRMAX,NRM,0)/BthV(NR+1)
       p_gr2phi2 = -RA**2*DERIV3(NR+2,R,ErV,NRMAX,NRM,0) &
            &      +RA*ErV(NR+2)*DERIV3(NR+2,R,BthV,NRMAX,NRM,0)/BthV(NR+2)
       p_gr2phi3 = -RA**2*DERIV3(NR+3,R,ErV,NRMAX,NRM,0) &
            &      +RA*ErV(NR+3)*DERIV3(NR+3,R,BthV,NRMAX,NRM,0)/BthV(NR+3)
       p_gr2phi = SNGL(AITKEN2P(R(0),p_gr2phi1,p_gr2phi2,p_gr2phi3,R(1),R(2),R(3)))
    ELSE
       p_gr2phi=SNGL(-RA**2*DERIV3(NR,R,ErV,NRMAX,NRM,0) &
            &        +RA*ErV(NR)*DERIV3(NR,R,BthV,NRMAX,NRM,0)/BthV(NR))
    END IF
    p_ngrth=SNGL(BphV(NR)/(RR*Q(NR)*BBL))
    temp_i(1) = SNGL(PTeV(NR))
    temp_i(2) = SNGL(PTiV(NR))
    grt_i(1)  = SNGL(RA * DERIV3(NR,R,PTeV,NRMAX,NRM,0))
    grt_i(2)  = SNGL(RA * DERIV3(NR,R,PTiV,NRMAX,NRM,0))
    den_iz(1,1)       = SNGL(PNeV(NR)) * 1.E20
    den_iz(2,INT(PZ)) = SNGL(PNiV(NR)) * 1.E20
    grp_iz(1) = SNGL(RA * DERIV3(NR,R,X(LQe5,0:NRMAX),NRMAX,NRM,0)) * 1.E20
    grp_iz(2) = SNGL(RA * DERIV3(NR,R,X(LQi5,0:NRMAX),NRMAX,NRM,0)) * 1.E20
    fex_iz(1:3,1:mx_mi,1:mx_mz) = 0.0

    CALL NCLASS( &
!     Input
         &      k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2, &
         &      p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
         &      p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz, &
!     Output
         &      m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii, &
         &      capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s, &
         &      sqz_s,upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s, &
         &      chip_ss,chit_ss,dp_ss,dt_ss,iflag)

    IF(k_out >0 .and.k_out <= 7) THEN
       p_eps = SNGL(R(NR)/RR)
       p_q   = SNGL(Q(NR))
       r0    = SNGL(RR)
       a0    = SNGL(RA)
       e0    = SNGL(1.D0)
       bt0   = SNGL(BB)
       q0l   = SNGL(Q0)
       CALL NCLASS_CHECK(6,NR, &
            &        k_out,k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2, &
            &        p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
            &        p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz, &
            &        p_eps,p_q,r0,a0,e0,bt0,q0l, &
            &        m_s,jm_s,jz_s,p_bsjb, &
            &        p_etap,p_exjb,calm_i,caln_ii, &
            &        capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s, &
            &        sqz_s,upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s, &
            &        chip_ss,chit_ss,dp_ss,dt_ss,iflag,ier_check)
       IER = ier_check
    ENDIF

    !     *** Takeover Parameters ***

    IF(k_potato == 0) THEN
       IF(iflag /= -1) THEN
          WRITE(6,*) "XX iflag=",iflag
          IER=1
          RETURN
       ENDIF
    ELSE
       IF(iflag /= 0) WRITE(6,*) "XX iflag=",iflag
    ENDIF

    JBSL =-(  DBLE(bsjbt_s(1)) *(RA * DERIV3(NR,R,PTeV,NRMAX,NRM,0) / PTeV(NR)) &
         &  + DBLE(bsjbp_s(1)) *(RA * DERIV3(NR,R,X(LQe5,0:NRMAX),NRMAX,NRM,0) &
         &                     / X(LQe5,NR)) &
         &  + DBLE(bsjbt_s(2)) *(RA * DERIV3(NR,R,PTiV,NRMAX,NRM,0) / PTiV(NR)) &
         &  + DBLE(bsjbp_s(2)) *(RA * DERIV3(NR,R,X(LQi5,0:NRMAX),NRMAX,NRM,0) &
         &                     / X(LQi5,NR))) / BphV(NR)

    ETAL = DBLE(p_etap)

    IF(NR == 0) THEN
       NueNC = 0.D0
       NuiNC = 0.D0
    ELSE
       NueNC = p_b2 * ymu_s(1,1,1) / (PNeV(NR) * 1.D20 * AME * BthV(NR)**2)
       NuiNC = p_b2 * ymu_s(1,1,2) / (PNiV(NR) * 1.D20 * AMI * BthV(NR)**2)
    END IF
    
!!$    AJBSNC(NR)=DBLE(p_bsjb)/BB
!!$    ETANC(NR) =DBLE(p_etap)
!!$    AJEXNC(NR)=DBLE(p_exjb)/BB
!!$
!!$    !     ADNCG : Diagonal diffusivity for graphic use only
!!$    !     AVNCG : Off-diagonal part driven pinch and neoclassical pinch
!!$    !             for graphic use only
!!$
!!$    DO NS=1,NSM
!!$       CJBSP(NR,NS)=DBLE(bsjbp_s(NS))
!!$       CJBST(NR,NS)=DBLE(bsjbt_s(NS))
!!$       ADNCG(NR,NS)=DBLE(dn_s(NS))/AR2RHO(NR)
!!$       AVNCG(NR,NS)=DBLE(vn_s(NS)+veb_s(NS))/AR1RHO(NR)
!!$       DO NS1=1,NSLMAX
!!$          AKNCP(NR,NS,NS1)=DBLE(chip_ss(NS,NS1))/AR2RHO(NR)
!!$          AKNCT(NR,NS,NS1)=DBLE(chit_ss(NS,NS1))/AR2RHO(NR)
!!$          ADNCP(NR,NS,NS1)=DBLE(dp_ss(NS,NS1))/AR2RHO(NR)
!!$          ADNCT(NR,NS,NS1)=DBLE(dt_ss(NS,NS1))/AR2RHO(NR)
!!$       ENDDO
!!$       DO NM=1,5
!!$          RGFLS(NR,NM,NS)=DBLE(gfl_s(NM,NS))*1.D-20/AR1RHO(NR)
!!$          RQFLS(NR,NM,NS)=DBLE(qfl_s(NM,NS))*1.D-20/AR1RHO(NR)
!!$       ENDDO
!!$       AVKNC(NR,NS)=DBLE(qeb_s(NS))/AR1RHO(NR)
!!$       AVNC (NR,NS)=DBLE(veb_s(NS))/AR1RHO(NR)
!!$    ENDDO
!!$    IF(MDLEQZ.NE.0) THEN
!!$       DO NSZ=1,NSZMAX
!!$          NS =NSM+NSZ
!!$          NSN=NSLMAX+NSZ
!!$          CJBSP(NR,NS)=DBLE(bsjbp_s(NSN))
!!$          CJBST(NR,NS)=DBLE(bsjbt_s(NSN))
!!$          ADNCG(NR,NS)=DBLE(dn_s(NSN))/AR2RHO(NR)
!!$          AVNCG(NR,NS)=DBLE(vn_s(NSN)+veb_s(NSN))/AR1RHO(NR)
!!$          DO NS1=1,NSLMAX
!!$             AKNCP(NR,NS,NS1)=DBLE(chip_ss(NSN,NS1))/AR2RHO(NR)
!!$             AKNCT(NR,NS,NS1)=DBLE(chit_ss(NSN,NS1))/AR2RHO(NR)
!!$             ADNCP(NR,NS,NS1)=DBLE(dp_ss(NSN,NS1))/AR2RHO(NR)
!!$             ADNCT(NR,NS,NS1)=DBLE(dt_ss(NSN,NS1))/AR2RHO(NR)
!!$          ENDDO
!!$          DO NM=1,5
!!$             RGFLS(NR,NM,NS)=DBLE(gfl_s(NM,NSN))*1.D-20/AR1RHO(NR)
!!$             RQFLS(NR,NM,NS)=DBLE(qfl_s(NM,NSN))*1.D-20/AR1RHO(NR)
!!$          ENDDO
!!$          AVKNC(NR,NS)=DBLE(qeb_s(NSN))/AR1RHO(NR)
!!$          AVNC (NR,NS)=DBLE(veb_s(NSN))/AR1RHO(NR)
!!$       ENDDO
!!$    ENDIF
!!$    !
!!$    !     /* neoclassical bulk ion toroidal and poloidal velocities */
!!$    !
!!$    IF(k_v /= 0) THEN
       btot=SQRT(BphV(NR)**2+BthV(NR)**2)*BphV(NR)/ABS(BphV(NR))
       DO i=1,m_s
          uthai=utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i)
          IF(DBLE(amu_i(jm_s(i))) == PA .and. jz_s(i) == INT(PZ)) then
             !     Poloidal
          VPOL(NR)=DBLE(uthai*BthV(NR))
!!$             !     Parallel
!!$             VPAR(NR)=DBLE(BthV(NR)/btot*VPOL(NR)+BphV(NR)/btot*VTOR(NR))
!!$             !     Perpendicular
!!$             VPRP(NR)=DBLE(BphV(NR)/btot*VPOL(NR)-BthV(NR)/btot*VTOR(NR))
          ENDIF
       ENDDO
!       IF(NT==500) write(6,*) R(NR),VPOL(NR)
!!$    ENDIF

    RETURN
  end SUBROUTINE TX_NCLASS

!***********************************************************
!
!    PARAMETER AND CONSISTENCY CHECK
!
!***********************************************************

  SUBROUTINE NCLASS_CHECK(nout,nr,k_out, &
       &            k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2, &
       &            p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
       &            p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz, &
       &            p_eps,p_q,r0,a0,e0,bt0,q0, &
       &            m_s,jm_s,jz_s,p_bsjb, &
       &            p_etap,p_exjb,calm_i,caln_ii, &
       &            capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s, &
       &            sqz_s,upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s, &
       &            chip_ss,chit_ss,dp_ss,dt_ss,iflag,ier)

    INCLUDE 'nclass/pamx_mi.inc'
    INCLUDE 'nclass/pamx_ms.inc'
    INCLUDE 'nclass/pamx_mz.inc'
    INCLUDE 'txncls.inc'
    !Declaration of local variables
    CHARACTER :: label*120
    INTEGER, intent(in)  :: nr, nout, k_out
    integer, intent(out) :: ier
    INTEGER   :: i, im, iz, iza, j, jm, jza, k, l
    INTEGER   :: idum(8)
    REAL      :: a0, bt0, bpol, btor, btot, e0, p_eps, p_q, ppr, q0, r0, uthai
    REAL      :: dq_s(mx_ms), vq_s(mx_ms)
    REAL      :: z_coulomb, z_electronmass, z_j7kv, z_mu0, z_pi, z_protonmass
    REAL      :: dum(8), edum(8), rdum(8)
    !Declaration of functions
    REAL      ::  RARRAY_SUM

    !Physical and conversion constants
    z_coulomb=1.6022e-19
    z_electronmass=9.1095e-31
    z_j7kv=1.6022e-16
    z_mu0=1.2566e-06
    z_pi=ACOS(-1.0)
    z_protonmass=1.6726e-27

    ier = 0

    IF(mod(k_out,2) == 1) THEN
       IF(iflag /= 0) WRITE(nout,'(A3,I3)') "NR=",NR
       !Check warning flags
       IF(iflag == -1) THEN
          label='WARNING:NCLASS-no potato orbit viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ELSEIF(iflag == -2) THEN
          label='WARNING:NCLASS-Pfirsch-Schluter viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ELSEIF(iflag == -3) THEN
          label='WARNING:NCLASS-no banana viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ELSEIF(iflag == -4) THEN
          label='WARNING:NCLASS-no viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ENDIF
    END IF
    !Check error flags
    IF(iflag > 0) THEN
       IF(int(k_out/2) == 1 .or. int(k_out/2) == 3) THEN
          IF(iflag == 1) THEN
             label='ERROR:NCLASS-k_order must be 2 or 3, k_order='
             idum(1)=k_order
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 2) THEN
             label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
             idum(1)=m_i
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 3) THEN
             label='ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
             idum(1)=m_z
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 4) THEN
             label='ERROR:NCLASS-require 0<m_s<mx_ms, m_s='
             idum(1)=m_s
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 5) THEN
             label='ERROR:NCLASS-inversion of flow matrix failed'
             CALL WRITE_LINE(nout,label,0,0)
          ELSEIF(iflag == 6) THEN
             label='ERROR:NCLASS-trapped fraction must be 0.0<p_ft<1.0'
             CALL WRITE_LINE(nout,label,0,0)
          ENDIF
       ENDIF
       WRITE(6,*) 'XX NCLASS_CHECK: non zero iflag: ',iflag
       ier = 1
       return
    ENDIF
    !Check for optional output
    IF(k_out >= 4) THEN
       !  Species identification
       label='     *** Species Identification ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Isotope     Species      Charge        Mass'// &
            &'     Density Temperature  Chg Factor'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -           -     Coulomb         AMU'// &
            &'       /m**3         keV           -'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          idum(1)=im
          idum(2)=i
          idum(3)=iz
          rdum(1)=amu_i(im)
          rdum(2)=den_iz(im,iza)
          rdum(3)=temp_i(im)
          rdum(4)=xi_s(i)
          CALL WRITE_IR(nout,3,idum,4,rdum,2)
       ENDDO
       !  Friction coefficients
       label='     *** Friction Coefficients ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO im=1,m_i
          DO jm=1,m_i
             CALL WRITE_LINE(nout,' ',0,0)            
             label='  Isotopes ='
             idum(1)=im
             idum(2)=jm
             CALL WRITE_LINE_IR(nout,label,2,idum,0,rdum,0)           
             !           Mkl
             label='  Mkl(-) ='
             CALL WRITE_LINE(nout,label,0,0)           
             DO k=1,k_order
                DO l=1,k_order
                   rdum(l)=capm_ii(k,l,im,jm)
                ENDDO
                CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
             ENDDO
             !           Nkl
             label='  Nkl(-) ='
             CALL WRITE_LINE(nout,label,0,0)           
             DO k=1,k_order
                DO l=1,k_order
                   rdum(l)=capn_ii(k,l,im,jm)
                ENDDO
                CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
             ENDDO
          ENDDO
       ENDDO
       !  Reduced friction coefficients
       label='     *** Reduced Friction Coefficients ***'
       CALL WRITE_LINE(nout,label,2,0)
       !       cal(M)
       DO im=1,m_i
          CALL WRITE_LINE(nout,' ',0,0)
          label='  calMij (kg/m**3/s) for Isotope ='
          idum(1)=im
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)           
          DO j=1,k_order
             DO k=1,k_order
                rdum(k)=calm_i(j,k,im)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !       cal(N)
       DO im=1,m_i
          CALL RARRAY_ZERO(3,dum)
          DO jm=1,m_i
             CALL WRITE_LINE(nout,' ',0,0)
             label='  calNij (kg/m**3/s) for Isotopes ='
             idum(1)=im
             idum(2)=jm
             CALL WRITE_LINE_IR(nout,label,2,idum,0,rdum,0)           
             DO k=1,k_order
                dum(k)=dum(k)-caln_ii(k,1,im,jm)
                DO l=1,k_order
                   rdum(l)=caln_ii(k,l,im,jm)
                ENDDO
                CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
             ENDDO
          ENDDO
          !         Momentum check
          CALL WRITE_LINE(nout,' ',0,0)
          label='Momentum check for Isotope ='
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          label='  sum_b[-calNk1] (kg/m**3/s) ='
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,dum,2)           
          label='  calMk1 (kg/m**3/s)         ='
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,calm_i(1,1,im),2)           
          label='  Difference ='
          DO k=1,k_order
             dum(k+3)=dum(k)-calm_i(k,1,im)
          ENDDO
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,dum(4),2)           
       ENDDO
       !  Normalized viscosities
       label='     *** Normalized Viscosities ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  muij (kg/m**3/s) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO k=1,k_order
             DO l=1,k_order
                rdum(l)=ymu_s(k,l,i)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !  Normalized parallel flows
       label='     *** Normalized Parallel Flows ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  upar_ij (T*m/s) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO j=1,k_order
             DO k=1,k_order
                rdum(k)=upar_s(j,k,i)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !  Normalized poloidal flows
       label='     *** Normalized Poloidal Flows ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  utheta_ij (m/s/T) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO j=1,k_order
             DO k=1,k_order
                rdum(k)=utheta_s(j,k,i)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !  Radial particle fluxes
       label='     *** Radial Particle Fluxes ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species          BP          PS          CL'// &
            &'       <E.B>         src       total'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -     /m**2/s     /m**2/s     /m**2/s'// &
            &'     /m**2/s     /m**2/s     /m**2/s'
       CALL WRITE_LINE(nout,label,0,0)
       CALL RARRAY_ZERO(6,dum)
       CALL RARRAY_ZERO(6,edum)
       !       Load into rdum the five flux components and total
       !       Load into dum the z-weighted + charge (ion) components
       !       Load into edum the z-weighted - charge (electron) components 
       DO i=1,m_s
          iz=jz_s(i)
          idum(1)=i
          CALL RARRAY_COPY(5,gfl_s(1,i),1,rdum,1)
          rdum(6)=RARRAY_SUM(5,rdum,1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          DO k=1,5
             IF(iz > 0) THEN
                dum(k)=dum(k)+iz*gfl_s(k,i)
             ELSE
                edum(k)=edum(k)+iz*gfl_s(k,i)
             ENDIF
          ENDDO
          IF(iz > 0) THEN
             dum(6)=dum(6)+iz*rdum(6)
          ELSE
             edum(6)=edum(6)+iz*rdum(6)
          ENDIF
       ENDDO
       !  Ambipolarity check
       label='     *** Ambipolarity Check ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='  Sum_i Z_i*Gamma_i ='
       CALL WRITE_LINE(nout,label,0,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,dum,2)           
       label='  Sum_e Z_e*Gamma_e ='
       CALL WRITE_LINE(nout,label,0,0)
       CALL WRITE_IR(nout,1,idum,6,edum,2)           
       label='  Difference ='
       CALL WRITE_LINE(nout,label,0,0)
       DO k=1,6
          rdum(k)=edum(k)+dum(k)
       ENDDO
       CALL WRITE_IR(nout,1,idum,6,rdum,2)
       !  Particle transport is returned in three forms
       !  Particle flux consistency check
       label='     *** Particle Flux Consistency Check ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species       gfl_s   dn_s,vn_s   dp_s,dt_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -     /m**2/s     /m**2/s     /m**2/s'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          !    Gamma total is sum of components
          rdum(1)=RARRAY_SUM(4,gfl_s(1,i),1)
          !    Gamma total is diffusive plus pinch
          rdum(2)=-dn_s(i)*(grp_iz(im,iza)-den_iz(im,iza)*grt_i(im)) &
               & /temp_i(im)+den_iz(im,iza)*(vn_s(i)+veb_s(i))
          !    Gamma total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
             jm=jm_s(j)
             jza=IABS(jz_s(j))
             rdum(3)=rdum(3)-dt_ss(j,i)*grt_i(jm)/temp_i(jm) &
                  & -dp_ss(j,i)*grp_iz(jm,jza)/den_iz(jm,jza)/temp_i(jm)
          ENDDO
          rdum(3)=(rdum(3)+veb_s(i))*den_iz(im,iza)
          CALL WRITE_IR(nout,1,idum,3,rdum,2)
       ENDDO
       !  Particle diffusion, velocity
       label='     *** Particle Diffusion, Velocity ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species         D_s         V_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s         m/s'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          rdum(1)=dn_s(i)
          rdum(2)=vn_s(i)+veb_s(i)+gfl_s(5,i)/den_iz(im,iza)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
       ENDDO
       !  Particle diffusivity matrices
       !    On p'/p           
       label='     *** Particle Diffusivity Matrices ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    on dp/dr ... by Species'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,dp_ss(1,i),2)
       ENDDO
       !    On T'/T
       label='     Species    on dT/dr ... by Species'
       CALL WRITE_LINE(nout,label,1,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,dt_ss(1,i),2)
       ENDDO
       !  Radial conduction fluxes
       label='     *** Radial Conduction Fluxes ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species          BP          PS          CL'// &
            &'       <E.B>         src       total'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      w/m**2      w/m**2      w/m**2'// &
            &'      w/m**2      w/m**2      w/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       CALL RARRAY_ZERO(6,dum)
       CALL RARRAY_ZERO(6,edum)
       DO i=1,m_s
          iz=jz_s(i)
          idum(1)=i
          CALL RARRAY_COPY(5,qfl_s(1,i),1,rdum,1)
          rdum(6)=RARRAY_SUM(5,qfl_s(1,i),1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          DO k=1,5
             IF(iz > 0) THEN
                dum(k)=dum(k)+qfl_s(k,i)
             ELSE
                edum(k)=edum(k)+qfl_s(k,i)
             ENDIF
          ENDDO
          IF(iz > 0) THEN
             dum(6)=dum(6)+rdum(6)
          ELSE
             edum(6)=edum(6)+rdum(6)
          ENDIF
       ENDDO
       !     Total electron radial conduction flux
       label='  Sum of electron conduction fluxes ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,edum,2)
       !     Total ion radial conduction flux
       label='  Sum of ion conduction fluxes ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,dum,2)
       !  Heat conduction is returned in two forms, add diagonal form
       !  Heat flux consistency check
       label='     *** Heat Flux Consistency Check ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species       qfl_s   dq_s,vq_s,  cp_s,ct_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      w/m**2      w/m**2      w/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          !    Conduction total is sum of components
          rdum(1)=RARRAY_SUM(5,qfl_s(1,i),1)
          !    Diagonal conductivity plus convective velocity
          dq_s(i)=chit_ss(i,i)+chip_ss(i,i)
          vq_s(i)=rdum(1)/(den_iz(im,iza)*z_j7kv*temp_i(im)) &
               &              +dq_s(i)*grt_i(im)/temp_i(im)
          rdum(2)=den_iz(im,iza)*z_j7kv*(-dq_s(i)*grt_i(im) &
               &              +vq_s(i)*temp_i(im))
          !    Conduction total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
             jm=jm_s(j)
             jza=IABS(jz_s(j))
             rdum(3)=rdum(3)-chit_ss(j,i)*grt_i(jm)/temp_i(jm) &
                  &         -chip_ss(j,i)*grp_iz(jm,jza) &
                  &                      /den_iz(jm,jza)/temp_i(jm)
          ENDDO
          rdum(3)=(rdum(3)+qeb_s(i))*den_iz(im,iza)*temp_i(im)*z_j7kv
          CALL WRITE_IR(nout,1,idum,3,rdum,2)
       ENDDO
       !  Heat conduction, velocity
       label='     *** Heat Conduction, Velocity ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species         X_s         V_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s         m/s'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          rdum(1)=dq_s(i)
          rdum(2)=vq_s(i)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
       ENDDO
       !     Effective ion heat conduction, velocity
       label='  Effective ion conduction, velocity ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL RARRAY_ZERO(4,rdum)
       DO i=1,m_s
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          IF(iz > 0) THEN
             rdum(1)=rdum(1)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza)
             rdum(2)=rdum(2)+RARRAY_SUM(5,qfl_s(1,i),1)/temp_i(im)/z_j7kv
             rdum(3)=rdum(1)+den_iz(im,iza)
             rdum(4)=rdum(4)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza) &
                  &                      *grt_i(im)/temp_i(im)
          ENDIF
       ENDDO
       rdum(1)=rdum(1)/rdum(3)
       rdum(2)=(rdum(2)+rdum(4))/rdum(3)
       CALL WRITE_IR(nout,1,idum,2,rdum,2)
       !  Thermal diffusivity matrices
       !    On p'/p           
       label='     *** Thermal Diffusivity Matrices ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    on dp/dr ... by Species'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,chip_ss(1,i),2)
       ENDDO
       !    On T'/T
       label='     Species    on dT/dr ... by Species'
       CALL WRITE_LINE(nout,label,1,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,chit_ss(1,i),2)
       ENDDO
       !  Radial energy (conduction+convection) fluxes
       label='     *** Radial Energy (Cond+Conv) Fluxes ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    ban-plat          PS   classical'// &
            &'       <E.B>  extern src       total'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      w/m**2      w/m**2      w/m**2'// &
            &'      w/m**2      w/m**2      w/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       CALL RARRAY_ZERO(6,dum)
       DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          DO k=1,5
             rdum(k)=qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
          ENDDO
          rdum(6)=RARRAY_SUM(5,rdum,1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          IF(iz > 0) THEN
             DO k=1,5
                dum(k)=dum(k)+qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
             ENDDO
             dum(6)=dum(6)+rdum(6)
          ENDIF
       ENDDO
       !     Total ion radial energy flux
       label='  Sum of ion energy fluxes ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,dum,2)           
       !  Bootstrap current is returned in two forms
       !  Bootstrap current consistency check
       label='     *** Bootstrap Current Consistency Check ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='      p_bsjb   dp_s,dt_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='    A*T/m**2    A*T/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       !    Total bootstrap current
       rdum(1)=p_bsjb
       !    Bootstrap current summed over T' and p' components
       rdum(2)=0.0
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          rdum(2)=rdum(2)-bsjbt_s(i)*grt_i(im)/temp_i(im) &
               & -bsjbp_s(i)*grp_iz(im,iza) &
               &               /den_iz(im,iza)/temp_i(im)
       ENDDO
       CALL WRITE_IR(nout,0,idum,2,rdum,2)
       !  Bootstrap current arrays
       !    On p'/p           
       label='     *** Bootstrap Current Arrays ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    on dp/dr    on dT/dr'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -    A*T/m**2    A*T/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          rdum(1)=bsjbp_s(i)
          rdum(2)=bsjbt_s(i)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
       ENDDO
       !  Current response to external source
       label='     *** Bootstrap and External Source Current ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='   Bootstrap    External       Total'
       CALL WRITE_LINE(nout,label,0,0)
       label='    A*T/m**2    A*T/m**2    A*T/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       rdum(1)=p_bsjb
       rdum(2)=p_exjb
       rdum(3)=rdum(1)+rdum(2)
       CALL WRITE_IR(nout,0,idum,3,rdum,2)
       !  Flow velocities on outside midplane
       label='     *** Flow Velocities on Outside Midplane ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species       v-tor       v-pol       v-par'// &
            &'      v-perp'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -         m/s         m/s         m/s'// &
            &'         m/s'
       CALL WRITE_LINE(nout,label,0,0)                        
       btor=bt0/(1.0+p_eps)
       bpol=btor/p_fhat
       btot=SQRT(btor**2+bpol**2)*btor/ABS(btor)
       DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          ppr=p_fhat*grp_iz(im,iza)*z_j7kv &
               & /(z_coulomb*iz*den_iz(im,iza))+p_fhat*p_grphi
          uthai=utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i)
          !         Toroidal
          rdum(1)=uthai*btor-ppr/btor
          !         Poloidal
          rdum(2)=uthai*bpol
          !         Parallel
          rdum(3)=uthai*btot-ppr/btot
          !         Perpendicular
          rdum(4)=ppr*bpol/btot/btor
          CALL WRITE_IR(nout,1,idum,4,rdum,2)           
       ENDDO
       !  Miscellaneous parameters
       label='     *** Miscellaneous Parameters ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='<J_bs.B>/Bt0 (A/m**2) ='
       rdum(1)=p_bsjb/bt0
       CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
       label='<J_ex.B>/Bt0 (A/m**2) ='
       rdum(1)=p_exjb/bt0
       CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
       label='eta parallel (Ohm*m) ='
       rdum(1)=p_etap
       CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
    ENDIF

    RETURN
  end SUBROUTINE NCLASS_CHECK
end module nclass_mod
