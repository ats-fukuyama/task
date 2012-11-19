MODULE trncls

  USE trcomm, ONLY: ikind,rkind,nrmax,nsamax

  PRIVATE
  PUBLIC tr_nclass,tr_ncls_allocate, &
       chi_ncp,chi_nct,d_ncp,d_nct,gfls,qfls,fls_tot,vebs,qebs,dia_gdnc, &
       dia_gvnc,cjbs_p,cjbs_t,eta_ncls,jbs_ncls,jex_ncls,vpol_ncls,       &
       vpar_ncls,vprp_ncls

  ! ----- for NCLASS output -----                                          
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
       qfls       ! radial heat conduction flux components of s [W/(m^2 s)]
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       fls_tot    ! total radial flux of s (heat and particle)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE::&
       vebs,     &! <E.B> particle convection velocity of s [m/s]
       qebs,     &! <E.B> heat convection velocity of s [m/s]
       dia_gdnc, &! diagonal diffusivity [m^2/s]
       dia_gvnc, &! diagonal convection driven by off-diagonal part [m/s]
       cjbs_p,   &! <J_bs.B> driven by unit p'/p of s [A*T/1.d-20*m^2]
       cjbs_t     ! <J_bs.B> driven by unit T'/T of s [A*T/1.d-20*m^2]
  REAL(rkind),DIMENSION(:),ALLOCATABLE::&
       eta_ncls, &!
       jbs_ncls, &!
       jex_ncls, &!
       vpol_ncls,&!
       vpar_ncls,&!
       vprp_ncls  !

CONTAINS

! **********************************************************************
!
!  The interface between TASK/TR(trcalv) and NCLASS model(NCLASS)
!
! **********************************************************************

  SUBROUTINE tr_nclass(ierr)
!***********************************************************************
!TR_NCLASS CALCULATES various parameters and arrays for NCLASS.
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
    USE trcomm, ONLY : &
         rkev,nsm,nrmax,nsamax,nsabmax,ns_nsa,nsab_nsa,idnsa,idion,pa,pz, &
         RR,ra,BB,rkap,abb1rho,abb2rho,aib2rho,ttrho,ar1rho,ar2rho,epsrho,&
         rhog,rhom,rt,rn,rp,dpdrho,pts,pns,qp,q0,bp,joh,er,vtor,eta
!    PADD,! additional pressure due to NBI
!    MDLTPF,! Trapped particle fraction model

    USE trcalnc, ONLY: ftpf

    IMPLICIT NONE
    INCLUDE 'nclass/pamx_mi.inc'
    INCLUDE 'nclass/pamx_ms.inc'
    INCLUDE 'nclass/pamx_mz.inc'
!    INCLUDE 'trncls.inc'

    ! Declaration of input to NCLASS
    INTEGER(4) :: k_order, k_potato, m_i, m_z
    REAL(4) :: &
         c_den,    c_potb,   c_potl,   p_b2,     p_bm2,    p_eb,   &
         p_fhat,   p_ft,   p_grbm2,  p_ngrth,  p_grphi,  p_gr2phi
    REAL(4),DIMENSION(3)            ::  p_fm
    REAL(4),DIMENSION(mx_mi)        ::  amu_i,    grt_i,    temp_i
    REAL(4),DIMENSION(mx_mi,mx_mz)  ::  den_iz,   grp_iz
    REAL(4),DIMENSION(3,mx_mi,mx_mz)::  fex_iz

    ! Declaration of output from NCLASS
    INTEGER(4)::                        iflag,    m_s
    INTEGER(4),DIMENSION(mx_ms)::       jm_s,     jz_s
    REAL(4)::                           p_bsjb,   p_etap,    p_exjb
    REAL(4),DIMENSION(3,3,mx_mi)::      calm_i
    REAL(4),DIMENSION(3,3,mx_mi,mx_mi)::caln_ii,  capm_ii,  capn_ii
    REAL(4),DIMENSION(mx_ms)::&
         bsjbp_s,  bsjbt_s,  dn_s,  sqz_s,  vn_s,  veb_s,  qeb_s,  xi_s
    REAL(4),DIMENSION(5,mx_ms)::        gfl_s,    qfl_s
    REAL(4),DIMENSION(3,3,mx_ms)::      upar_s,   utheta_s, ymu_s
    REAL(4),DIMENSION(mx_ms,mx_ms)::    chip_ss,  chit_ss,  dp_ss,   dt_ss
    
    INTEGER(4),INTENT(OUT):: IERR
    INTEGER(4)::  i,iz,k_out,k_v,na,nm,nr,ns,ns1,nsn,nsz
    REAL(4)   ::  a0,bt0,e0,p_eps,p_q,q0l,r0
    REAL(8)   ::  bpol,btor,btot,uthai
    REAL(8)   ::  aitken2p,deriv4,deriv3,FCTR

    ! internal parameters for tr_nclass
    INTEGER(4) :: mdltpf ! interim definition of switch variables
    INTEGER(4) :: nsa,nsa1,nk,nsabnum, nsab, nsab1
    REAL(8),DIMENSION(0:nrmax) :: eropsi,nr_array,drhog
    REAL(8) :: abb1rhom,abb2rhom,aib2rhom, &
         ttrhom,dpdrhom,ar1rhom,ar2rhom,epsrhom,qpm,bpm,etam,johm
    REAL(8) :: cm,cp

    mdltpf=0

    ! *** Initialization ***
    amu_i(1:mx_mi)  = 0.0
    grt_i(1:mx_mi)  = 0.0
    temp_i(1:mx_mi) = 0.0
    den_iz(1:mx_mi,1:mx_mz)     = 0.0
    grp_iz(1:mx_mi,1:mx_mz)     = 0.0
    fex_iz(1:3,1:mx_mi,1:mx_mz) = 0.0
    
    ierr = 0

!  ----- < Output option switches > ------------------------
!   k_out: option for output to nout (-)
!         = 1    : errors only
!         = 2    : errors and results
!         = else : no output
!
!   k_v  : option for neoclassical v_tor,v_pol,v_para,v_perp
!         = 1    : output
!         = else : no output
! ----------------------------------------------------------
    k_out    = 1
    k_v      = 1
    
    k_order  = 2
    k_potato = 1
      
    !*** mdleqz --> nsamax
    ! number of isotopes (1 < m_i < mx_mi+1)
    m_i = nsabmax

    !*** mdleqz,pzmax --> pz(nsamax)
    ! highest charge state of all species (0 < m_z < mx_mz+1)
    m_z = INT(MAXVAL(pz(1:nsamax)))

    c_den  = SNGL(1.E10)
    c_potb = SNGL(rkap*BB/(2.d0*q0**2))
    c_potl = SNGL(q0*RR)

    ! atomic mass number of i (-)
    amu_i(1:nsamax) = SNGL(PA(1:nsamax))

    ! ----- preparation for calculatio in nr loop -----
    eropsi(1:nrmax) = er(1:nrmax)                      &
         /(0.5d0*(ar1rho(0:nrmax-1)+ar1rho(1:nrmax)))  &
         /(0.5d0*(dpdrho(0:nrmax-1)+dpdrho(1:nrmax)))

      
    DO nr = 1, nrmax
       drhog(nr) = rhog(nr)-rhog(nr-1)
       abb1rhom  = 0.5d0*(abb1rho(nr) + abb1rho(nr-1))
       abb2rhom  = 0.5d0*(abb2rho(nr) + abb2rho(nr-1))
       aib2rhom  = 0.5d0*(aib2rho(nr) + aib2rho(nr-1))
       ar1rhom   = 0.5d0*(ar1rho(nr)  +  ar1rho(nr-1))
       ar2rhom   = 0.5d0*(ar2rho(nr)  +  ar2rho(nr-1))
       ttrhom    = 0.5d0*(ttrho(nr)   +   ttrho(nr-1))
       dpdrhom   = 0.5d0*(dpdrho(nr)  +  dpdrho(nr-1))
       epsrhom   = 0.5d0*(epsrho(nr)  +  epsrho(nr-1))
       qpm       = 0.5d0*(qp(nr)      +      qp(nr-1))
       bpm       = 0.5d0*(bp(nr)      +      bp(nr-1))
       johm      = 0.5d0*(joh(nr)     +     joh(nr-1))
       etam      = eta_ncls(nr)

       p_b2      = SNGL(abb2rhom)
       p_bm2     = SNGL(aib2rhom)
       p_fhat    = SNGL(ttrhom/dpdrhom)
       p_fm(1:3) = 0.0
       IF(SNGL(epsrhom).gt.0.0) THEN
          DO i = 1, 3
             ! poloidal moments of geometic factor for PS viscosity (-)
             p_fm(i) = SNGL( & ! ??????????
                  DBLE(i)                                        &
                  *((1.D0-SQRT(1.D0-epsrhom**2))/epsrhom)**(2*i) &
                  *(1.D0+DBLE(i)*SQRT(1.D0-epsrhom**2))          &
                  /((1.D0-epsrhom**2)**1.5D0*(qpm*RR)**2)        &
                  )
          ENDDO
       ENDIF
       ! trapped particle fraction
       p_ft     = SNGL(ftpf(mdltpf,epsrhom))
       p_grbm2  = SNGL(ar2rhom*aib2rhom)
       ! radial electric field phi' (V/rho)
       p_grphi  = SNGL(er(nr)/ar1rhom)
       p_gr2phi = SNGL(dpdrhom*deriv3(nr,rhom(1:nrmax),eropsi(1:nrmax),nrmax,1))
!       p_gr2phi=0.0
       p_ngrth  = SNGL(bpm/(BB*rhom(nr))) ! ?????
       
       nsabnum = 0
       DO nsa = 1, nsamax ! only for bulk species
          ns = ns_nsa(nsa)
          IF(pz(ns) == 0.d0 .OR.  idion(ns) /= 0.d0) CYCLE
          nsabnum = nsabnum + 1
               
          ! temperature of i (keV)
          temp_i(nsabnum) = SNGL(0.5d0*(rt(nsa,nr)+rt(nsa,nr-1)))
          ! temperature gradient of i (keV/rho)
          grt_i(nsabnum)  = SNGL((rt(nsa,nr)-rt(nsa,nr-1))/drhog(nr))
          ! density of i,z (/m**3)
          den_iz(nsabnum,INT(ABS(pz(ns)))) &
                     = SNGL(0.5d0*(rn(nsa,nr)+rn(nsa,nr-1)))*1.E20
          ! pressure gradient of i,z (keV/m**3 /rho)
          ! ***** PADD should be included *****
          grp_iz(nsabnum,INT(ABS(pz(ns)))) &
                     = SNGL((rp(nsa,nr)-rp(nsa,nr-1))/(drhog(nr)*rkev))
          DO i = 1, 3
             ! moments of external parallel force on i,z (T*j/m**3)
             fex_iz(i,nsabnum,INT(ABS(pz(ns)))) = 0.0
          ENDDO
       ENDDO
       ! <E.B> (V*T/m)  *** E_para = eta_para * J_oh(para) ***
       p_eb = SNGL(etam*johm*abb1rhom)


       CALL NCLASS( &
           ! << Input >>
           k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,    &
           p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
           p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,      &
           ! << Output >>
           m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,    &
           capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,     &
           sqz_s,upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,    &
           chip_ss,chit_ss,dp_ss,dt_ss,iflag)

       IF(k_out.eq.1 .or. k_out.eq.2) THEN ! for error output
          p_eps = SNGL(epsrho(nr))
          p_q   = SNGL(qp(nr))
          r0    = SNGL(RR)
          a0    = SNGL(ra)
          e0    = SNGL(rkap)
          bt0   = SNGL(BB)
          q0l   = SNGL(q0)

          CALL NCLASS_CHECK(6,nr,                                &
              k_out,k_order,m_i,m_z,p_fhat,p_grphi,              &
              amu_i,grt_i,temp_i,den_iz,grp_iz,p_eps,bt0,        &
              m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii, &
              capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,  &
              upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,       &
              chip_ss,chit_ss,dp_ss,dt_ss,iflag)
         
       ENDIF


       ! *****  Takeover Parameters *****
       IF(k_potato.EQ.0) THEN
          IF(iflag.NE.-1) THEN
             WRITE(6,*) "XXX NCLASS error: iflag =",iflag
             ierr=1
             RETURN
          ENDIF
       ELSE
          IF(iflag.NE. 0) WRITE(6,*) "XXX NCLASS error: iflag =",iflag
       ENDIF

       ! *** interim substitution 
       !  (Values on half grids -> arrays on grids)
       eta_ncls(nr) = DBLE(p_etap)
       jbs_ncls(nr) = DBLE(p_bsjb)/abb1rhom
       jex_ncls(nr) = DBLE(p_exjb)/abb1rhom

       ! For now, contribution of the flux driven by the external force
       !  gfl_s(5,i)/n_i is not included.
       DO nsa = 1, nsabnum
          nsab = nsab_nsa(nsa)
          IF(nsab /= 0)THEN
             cjbs_p(nsab,nr) = DBLE(bsjbp_s(nsa))
             cjbs_t(nsab,nr) = DBLE(bsjbt_s(nsa))

             ! flux form
             DO nk  = 1, 5
                gfls(nk,nsab,nr) = DBLE(gfl_s(nk,nsa))*1.d-20/ar1rhom
                qfls(nk,nsab,nr) = DBLE(qfl_s(nk,nsa))*1.d-20/ar1rhom
             ENDDO
          
             ! ADNCG : Diagonal diffusivity for graphic use only
             dia_gdnc(nsab,nr) = DBLE(dn_s(nsa))/ar2rhom

             ! AVNCG : Off-diagonal part driven pinch and neoclassical pinch
             !          for graphic use only
             dia_gvnc(nsab,nr) = DBLE(vn_s(nsa)+veb_s(nsa))/ar1rhom
!             dia_gvnc(nsa,nr) = DBLE(vn_s(nsa))/ar1rhom

             vebs(nsab,nr) = DBLE(veb_s(nsa))/ar1rhom
             qebs(nsab,nr) = DBLE(qeb_s(nsa))/ar1rhom

             ! complete matrix form as coefficients of 
             !  the temperature gradients and the pressure gradients.
             DO nsa1 = 1, nsamax
                nsab1 = nsab_nsa(nsa1)
                IF(nsab1 /= 0)THEN
                   chi_ncp(nsab,nsab1,nr) = DBLE(chip_ss(nsa,nsa1)) / ar2rhom
                   chi_nct(nsab,nsab1,nr) = DBLE(chit_ss(nsa,nsa1)) / ar2rhom

                   d_ncp(nsab,nsab1,nr)   = DBLE(  dp_ss(nsa,nsa1)) / ar2rhom
                   d_nct(nsab,nsab1,nr)   = DBLE(  dt_ss(nsa,nsa1)) / ar2rhom
                END IF
             ENDDO
          END IF
       ENDDO

       ! neoclassical bulk ion toroidal and poloidal velocities
       IF(k_v.ne.0) THEN
          btor = bt0 / (1.0 + 0.5d0*p_eps**2)
          bpol = ar1rhom/r0 * dpdrhom
          btot = SQRT(btor**2+bpol**2)*btor / ABS(btor)

          DO i = 1, m_s
             ! poloidal flow of s
             uthai = utheta_s(1,1,i) + utheta_s(1,2,i) + utheta_s(1,3,i)
             
             IF(DBLE(amu_i(jm_s(i))).EQ.PA(2) &
                       .AND. jz_s(i).EQ.INT(PZ(2))) THEN
                ! poloidal
                vpol_ncls(nr)=DBLE(uthai*bpol)
                ! parallel
                vpar_ncls(nr)=DBLE(bpol/btot*vpol_ncls(nr)+btor/btot*vtor(nr))
                ! perpendicular
                vprp_ncls(nr)=DBLE(btor/btot*vpol_ncls(nr)-bpol/btot*vtor(nr))
             ENDIF
          ENDDO
       ENDIF

    ENDDO ! end of nr loop

    RETURN
  END SUBROUTINE tr_nclass

!     ***********************************************************

!            PARAMETER AND CONSISTENCY CHECK

!     ***********************************************************

  SUBROUTINE NCLASS_CHECK(nout,nr,                                   &
     &            k_out,k_order,m_i,m_z,p_fhat,p_grphi,              &
     &            amu_i,grt_i,temp_i,den_iz,grp_iz,p_eps,bt0,        &
     &            m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii, &
     &            capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,  &
     &            upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,       &
     &            chip_ss,chit_ss,dp_ss,dt_ss,iflag)

      IMPLICIT NONE
      INCLUDE 'nclass/pamx_mi.inc'
      INCLUDE 'nclass/pamx_ms.inc'
      INCLUDE 'nclass/pamx_mz.inc'
!      INCLUDE 'trncls.inc'
!Declaration of input to NCLASS
      INTEGER(4)::                        k_order,  m_i,      m_z
      REAL(4)::                           p_fhat,   p_grphi
      REAL(4),DIMENSION(mx_mi)        ::  amu_i,    grt_i,    temp_i
      REAL(4),DIMENSION(mx_mi,mx_mz)  ::  den_iz,   grp_iz

!Declaration of output from NCLASS
      INTEGER(4)::                        iflag,    m_s
      INTEGER(4),DIMENSION(mx_ms)::       jm_s,     jz_s
      REAL(4)::                           p_bsjb,   p_etap,    p_exjb
      REAL(4),DIMENSION(3,3,mx_mi)::      calm_i
      REAL(4),DIMENSION(3,3,mx_mi,mx_mi)::caln_ii,  capm_ii,  capn_ii
      REAL(4),DIMENSION(mx_ms)::          bsjbp_s,  bsjbt_s,  dn_s,     vn_s,     veb_s,    qeb_s,    xi_s
      REAL(4),DIMENSION(5,mx_ms)::        gfl_s,    qfl_s
      REAL(4),DIMENSION(3,3,mx_ms)::      upar_s,   utheta_s, ymu_s
      REAL(4),DIMENSION(mx_ms,mx_ms)::    chip_ss,  chit_ss,  dp_ss,    dt_ss

!Declaration of local variables
      CHARACTER(LEN=120):: label
      INTEGER        nr
      INTEGER        i,                       im, &
     &               iz,                      iza, &
     &               j,                       jm, &
     &               jza,                     k, &
     &               k_out,                   l, &
     &               nout
      INTEGER        idum(8)
      REAL           bt0, &
     &               bpol,                    btor, &
     &               btot, &
     &               p_eps, &
     &               ppr, &
     &               uthai
      REAL           dq_s(mx_ms),             vq_s(mx_ms)
      REAL           z_coulomb,               z_electronmass, &
     &               z_j7kv,                  z_mu0, &
     &               z_pi,                    z_protonmass
      REAL           dum(8),                  edum(8), &
     &               rdum(8)
!Declaration of functions
      REAL           RARRAY_SUM
!Physical and conversion constants
      z_coulomb=1.6022e-19
      z_electronmass=9.1095e-31
      z_j7kv=1.6022e-16
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
      z_protonmass=1.6726e-27

      IF(iflag.ne.0) WRITE(nout,'(A3,I3)') "NR=",nr
!Check warning flags
      IF(iflag.eq.-1) THEN
        label='WARNING:NCLASS-no potato orbit viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ELSEIF(iflag.eq.-2) THEN
        label='WARNING:NCLASS-Pfirsch-Schluter viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ELSEIF(iflag.eq.-3) THEN
        label='WARNING:NCLASS-no banana viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ELSEIF(iflag.eq.-4) THEN
        label='WARNING:NCLASS-no viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ENDIF
!Check error flags
      IF(iflag.gt.0) THEN
        IF(k_out.gt.0) THEN
          IF(iflag.eq.1) THEN
            label='ERROR:NCLASS-k_order must be 2 or 3, k_order='
            idum(1)=k_order
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.2) THEN
            label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
            idum(1)=m_i
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.3) THEN
            label='ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
            idum(1)=m_z
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.4) THEN
            label='ERROR:NCLASS-require 1<m_s<mx_ms, m_s='
            idum(1)=m_s
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.5) THEN
            label='ERROR:NCLASS-inversion of flow matrix failed'
            CALL WRITE_LINE(nout,label,0,0)
          ELSEIF(iflag.eq.6) THEN
            label='ERROR:NCLASS-trapped fraction must be 0.0<p_ft<1.0'
            CALL WRITE_LINE(nout,label,0,0)
          ENDIF
        ENDIF
        WRITE(6,*) 'XX NCLASS_CHECK: non zero iflag',iflag
        STOP
      ENDIF
!Check for optional output
      IF(k_out.gt.1) THEN
!  Species identification
        label='     *** Species Identification ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='     Isotope     Species      Charge        Mass'// &
     &        '     Density Temperature  Chg Factor'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -           -     Coulomb         AMU'// &
     &        '       /m**3         keV           -'
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
     &        '       <E.B>         src       total'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -     /m**2/s     /m**2/s     /m**2/s'// &
     &        '     /m**2/s     /m**2/s     /m**2/s'
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
            IF(iz.gt.0) THEN
              dum(k)=dum(k)+iz*gfl_s(k,i)
            ELSE
              edum(k)=edum(k)+iz*gfl_s(k,i)
            ENDIF
          ENDDO
          IF(iz.gt.0) THEN
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
     &            /temp_i(im)+den_iz(im,iza)*(vn_s(i)+veb_s(i))
!    Gamma total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
            jm=jm_s(j)
            jza=IABS(jz_s(j))
            rdum(3)=rdum(3)-dt_ss(j,i)*grt_i(jm)/temp_i(jm) &
     &              -dp_ss(j,i)*grp_iz(jm,jza)/den_iz(jm,jza)/temp_i(jm)
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
     &        '       <E.B>         src       total'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      w/m**2      w/m**2      w/m**2'// &
     &        '      w/m**2      w/m**2      w/m**2'
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
            IF(iz.gt.0) THEN
              dum(k)=dum(k)+qfl_s(k,i)
            ELSE
              edum(k)=edum(k)+qfl_s(k,i)
            ENDIF
          ENDDO
          IF(iz.gt.0) THEN
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
          vq_s(i)=rdum(1)/(den_iz(im,iza)*z_j7kv*temp_i(im))+dq_s(i)*grt_i(im)/temp_i(im)
          rdum(2)=den_iz(im,iza)*z_j7kv*(-dq_s(i)*grt_i(im) +vq_s(i)*temp_i(im))
!    Conduction total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
            jm=jm_s(j)
            jza=IABS(jz_s(j))
            rdum(3)=rdum(3)-chit_ss(j,i)*grt_i(jm)/temp_i(jm) &
     &                     -chip_ss(j,i)*grp_iz(jm,jza)/den_iz(jm,jza)/temp_i(jm)
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
          IF(iz.gt.0) THEN
            rdum(1)=rdum(1)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza)
            rdum(2)=rdum(2)+RARRAY_SUM(5,qfl_s(1,i),1)/temp_i(im)/z_j7kv
            rdum(3)=rdum(1)+den_iz(im,iza)
            rdum(4)=rdum(4)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza)*grt_i(im)/temp_i(im)
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
     &        '       <E.B>  extern src       total'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      w/m**2      w/m**2      w/m**2'// &
     &        '      w/m**2      w/m**2      w/m**2'
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
          IF(iz.gt.0) THEN
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
     &                   -bsjbp_s(i)*grp_iz(im,iza)/den_iz(im,iza)/temp_i(im)
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
     &        '      v-perp'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -         m/s         m/s         m/s'// &
     &        '         m/s'
        CALL WRITE_LINE(nout,label,0,0)
        btor=bt0/(1.0+p_eps)
        bpol=btor/p_fhat
        btot=SQRT(btor**2+bpol**2)*btor/ABS(btor)
        DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          ppr=p_fhat*grp_iz(im,iza)*z_j7kv/(z_coulomb*iz*den_iz(im,iza))+p_fhat*p_grphi
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
      END SUBROUTINE NCLASS_CHECK

! ************************************************************************

  SUBROUTINE tr_ncls_allocate

    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind),SAVE:: nsamax_save=0
    INTEGER(ikind)     :: ierr

    IF(nrmax /= nrmax_save .OR. nsamax /= nsamax_save)THEN

       IF(nrmax_save /= 0) CALL tr_ncls_deallocate

       nr_allocation: DO
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
          ALLOCATE(eta_ncls(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jbs_ncls(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jex_ncls(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vpol_ncls(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vpar_ncls(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vprp_ncls(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT

          nrmax_save  = nrmax
          nsamax_save = nsamax
          RETURN
       END DO nr_allocation

       WRITE(6,*) 'XX tr_nr_allocate: allocation error: ierr=',ierr
       STOP
    END IF

    RETURN
  END SUBROUTINE tr_ncls_allocate

  SUBROUTINE tr_ncls_deallocate

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
    IF(ALLOCATED(eta_ncls)) DEALLOCATE(eta_ncls)
    IF(ALLOCATED(jbs_ncls)) DEALLOCATE(jbs_ncls)
    IF(ALLOCATED(jex_ncls)) DEALLOCATE(jex_ncls)
    IF(ALLOCATED(vpol_ncls)) DEALLOCATE(vpol_ncls)
    IF(ALLOCATED(vpar_ncls)) DEALLOCATE(vpar_ncls)
    IF(ALLOCATED(vprp_ncls)) DEALLOCATE(vprp_ncls)

    RETURN
  END SUBROUTINE tr_ncls_deallocate

END MODULE trncls
