C  
C     ***********************************************************
C
C            NCLASS
C
C     ***********************************************************
C
      SUBROUTINE TR_NCLASS
!***********************************************************************
!TR_NCLASS calculates various paramters and arrays for TR_NCLASS_DRIVER.
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
!  c_den - density cutoff below which species is ignored (/m**3)
!  default - switch parameter
!          =0 use default parameters described by nclass_pt_dr.f
!          =1 use parameters calculated here (strongly recommended)
!  amu_i(i) - atomic mass number of i (-)
!  p_grphi - radial electric field Phi' (V/rho)
!  p_gr2phi - radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  grt_i(i) - temperature gradient of i (keV/rho)
!  temp_i(i) - temperature of i (keV)
!  den_iz(i,z) - density of i,z (/m**3)
!  fex_iz(3,i,z) - moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z) - pressure gradient of i,z (keV/m**3/rho)
!  p_eps - inverse aspect ratio (-) given by EPSRHO(NR)
!  p_q - safety factor (-) given by QP(NR)
!  r0 - major radius (m) given by RR
!  a0 - minor radius (m) given by RA
!  e0 - axial elongation (-) given by RKAP
!  bt0 - axial toroidal magnetic field (T) given by BB
!  q0 - axial safety factor (-) given by Q0
!***********************************************************************
      INCLUDE 'trcomm.h'
      INCLUDE 'nclass/pamx_mi.inc'
      INCLUDE 'nclass/pamx_ms.inc'
      INCLUDE 'nclass/pamx_mz.inc'
      INCLUDE 'trncls.inc'
C
      k_order=2
      k_potato=0
      c_den=1.E10
      default=1
C
      DO NS=1,NSMAX
         amu_i(NS)=SNGL(PA(NS))
      ENDDO
C
      DO NR=1,NRMAX
         p_grphi=0.E0
         p_gr2phi=0.E0
         IF(NR.EQ.NRMAX) THEN
            DO NS=1,NSMAX
               temp_i(NS)=SNGL(PTS(NS))
               grt_i(NS)=SNGL(2.D0*(PTS(NS)-RN(NR,NS))*AR1RHO(NR)/DR)
               den_iz(NS,INT(ABS(PZ(NS))))=SNGL(PNSS(NS))*1.E20
               grp_iz(NS,INT(ABS(PZ(NS))))
     &              =SNGL(2.D0*(PNSS(NS)*PTS(NS)-RN(NR,NS)*RT(NR,NS))
     &              *AR1RHO(NR)/DR)*1.E20
               DO NA=1,3
                  fex_iz(NA,NS,INT(ABS(PZ(NS))))=0.E0
               ENDDO
            ENDDO
         ELSE
            DO NS=1,NSMAX
               temp_i(NS)=SNGL(0.5D0*(RT(NR+1,NS)+RT(NR,NS)))
               grt_i(NS)=SNGL((RT(NR+1,NS)-RT(NR,NS))*AR1RHO(NR)/DR)
               den_iz(NS,INT(ABS(PZ(NS))))=SNGL(0.5D0*(RN(NR+1,NS)
     &                                                +RN(NR,NS)))*1.E20
               grp_iz(NS,INT(ABS(PZ(NS))))
     &             =SNGL((RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))
     &              *AR1RHO(NR)/DR)*1.E20
               DO NA=1,3
                  fex_iz(NA,NS,INT(ABS(PZ(NS))))=0.E0
               ENDDO
            ENDDO
         ENDIF
C
         CALL TR_NCLASS_DRIVER(
C     Input
     &                         SNGL(EPSRHO(NR)),SNGL(QP(NR)),SNGL(RR),
     &                         SNGL(RA),SNGL(RKAP),SNGL(BB),SNGL(Q0),
     &                         k_order,k_potato,c_den,p_grphi,
     &                         p_gr2phi,amu_i,grt_i,temp_i,
     &                         den_iz,fex_iz,grp_iz,default,
C     Output
     &                         iflag,m_s,p_bsjb,p_etap,p_exjb,
     &                         calm_i,caln_ii,capm_ii,capn_ii,
     &                         bsjbp_s,bsjbt_s,chip_ss,chit_ss,dn_s,
     &                         dp_ss,dt_ss,gfl_s,jm_s,jz_s,qeb_s,
     &                         qfl_s,sqz_s,upar_s,utheta_s,vn_s,
     &                         veb_s,xi_s,ymu_s)
C
C     *** Takeover Parameters ***
C
         IF(k_potato.EQ.0) THEN
            IF(iflag.NE.-1) WRITE(6,*) "XX iflag=",iflag
         ELSE
            IF(iflag.NE. 0) WRITE(6,*) "XX iflag=",iflag
         ENDIF
         AJBSNC(NR)=DBLE(p_bsjb)/BB
         ETANC(NR) =DBLE(p_etap)
         AJEXNC(NR)=DBLE(p_exjb)/BB
C
         DO NS=1,NSMAX
            CJBSP(NR,NS)=DBLE(bsjbp_s(NS))/AR1RHO(NR)
            CJBST(NR,NS)=DBLE(bsjbt_s(NS))/AR1RHO(NR)
C            if(ns.eq.2) write(6,*) NR,CJBSP(NR,NS),CJBST(NR,NS)
            DO NS1=1,NSMAX
               AKNCP(NR,NS,NS1)=DBLE(chip_ss(NS,NS1))/AR2RHO(NR)
               AKNCT(NR,NS,NS1)=DBLE(chit_ss(NS,NS1))/AR2RHO(NR)
               ADNCP(NR,NS,NS1)=DBLE(dp_ss(NS,NS1))/AR2RHO(NR)
               ADNCT(NR,NS,NS1)=DBLE(dt_ss(NS,NS1))/AR2RHO(NR)
c$$$               if(ns.eq.ns1.and.ns.eq.1) ! or ns.eq.2
c$$$     &              write(6,*) NR,ADNCP(NR,NS,NS1)+ADNCT(NR,NS,NS1),
c$$$     &                            AKNCP(NR,NS,NS1)+AKNCT(NR,NS,NS1)
            ENDDO
            DO NM=1,5
               RGFLS(NR,NM,NS)=DBLE(gfl_s(NM,NS))*1.D-20/AR1RHO(NR)
               RQFLS(NR,NM,NS)=DBLE(qfl_s(NM,NS))*1.D-20/AR1RHO(NR)
            ENDDO
            AVKNCS(NR,NS)=DBLE(qeb_s(NS))/AR1RHO(NR)
            AVNCS (NR,NS)=DBLE(veb_s(NS))/AR1RHO(NR)
         ENDDO
C
      ENDDO
C
C      STOP
      RETURN
      END
C  
C     ***********************************************************
C
C           NCLASS DRIVER
C
C     ***********************************************************
C
      SUBROUTINE TR_NCLASS_DRIVER(
C     Input
     &                            p_eps,p_q,r0,a0,e0,bt0,q0,
     &                            k_order,k_potato,c_den,p_grphi,
     &                            p_gr2phi,amu_i,grt_i,temp_i,
     &                            den_iz,fex_iz,grp_iz,default,
C     Output
     &                            iflag,m_s,p_bsjb,p_etap,p_exjb,
     &                            calm_i,caln_ii,capm_ii,capn_ii,
     &                            bsjbp_s,bsjbt_s,chip_ss,chit_ss,dn_s,
     &                            dp_ss,dt_ss,gfl_s,jm_s,jz_s,qeb_s,
     &                            qfl_s,sqz_s,upar_s,utheta_s,vn_s,
     &                            veb_s,xi_s,ymu_s)
C
      IMPLICIT NONE
      INCLUDE 'nclass/pamx_mi.inc'
      INCLUDE 'nclass/pamx_ms.inc'
      INCLUDE 'nclass/pamx_mz.inc'
      INTEGER        default
!Declaration of input to NCLASS
      INTEGER        k_order,                 k_potato
      INTEGER        m_i,                     m_z
      REAL           c_den,                   c_potb,
     #               c_potl
      REAL           p_b2,                    p_bm2,
     #               p_eb,                    p_fhat,
     #               p_fm(3),                 p_ft,
     #               p_grbm2,                 p_grphi,
     #               p_gr2phi,                p_ngrth
      REAL           amu_i(mx_mi),            grt_i(mx_mi),
     #               temp_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz),
     #               grp_iz(mx_mi,mx_mz)
!Declaration of output from NCLASS
      INTEGER        iflag,                   m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           p_bsjb,                  p_etap,
     #               p_exjb
      REAL           calm_i(3,3,mx_mi)
      REAL           caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),
     #               capn_ii(3,3,mx_mi,mx_mi)
      REAL           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     #               dn_s(mx_ms),             gfl_s(5,mx_ms),
     #               qfl_s(5,mx_ms),          sqz_s(mx_ms),
     #               upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms),
     #               vn_s(mx_ms),             veb_s(mx_ms),
     #               qeb_s(mx_ms),            xi_s(mx_ms),
     #               ymu_s(3,3,mx_ms)
      REAL           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     #               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
      INTEGER        i,                       im,
     #               iz,                      iza,
     #               j,                       jm,
     #               jza,                     k,
     #               l
      REAL           a0,                      bt0,
     #               bpol,                    btor,
     #               btot,                    ds,
     #               e0,                      p_eps,
     #               p_q,                     ppr,
     #               q0,                      r0
      REAL           z_coulomb,               z_electronmass,
     #               z_j7kv,                  z_mu0,
     #               z_pi,                    z_protonmass
!Declaration of functions
      REAL           RARRAY_SUM
!Physical and conversion constants
      z_coulomb=1.6022e-19
      z_electronmass=9.1095e-31
      z_j7kv=1.6022e-16
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
      z_protonmass=1.6726e-27
C
      IF(default.EQ.0) THEN
!Set default data
!  k_order-order of v moments to be solved (-)
!         =2 u and p_q
!         =3 u, p_q, and u2
!         =else error
      k_order=2
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
      k_potato=0
!  c_den-density cutoff below which species is ignored (/m**3)
      c_den=1.0e10
!  p_eps-inverse aspect ratio (-)
      p_eps=0.1
!  p_grphi-radial electric field Phi' (V/rho)
      p_grphi=-5.0e3
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
      p_gr2phi=-2.0e4
!  p_q-safety factor (-)
      p_q=-3.0
!  r0-major radius (m)
      r0=3.0
!  a0-minor radius, scale length for rho (m)
      a0=1.0
!  e0-axial elongation (-)
      e0=1.0
!  bt0-axial toroidal field (T)
      bt0=-3.0
!  q0-axial safety factor (-)
      q0=-3.0
      ENDIF
C
!Calculate input for NCLASS
!  m_i-number of isotopes (1<mi<mx_mi+1)
!  m_z-highest charge state of all species (0<mz<mx_mz+1)
      m_i=0
      m_z=0
	DO i=1,mx_mi
        ds=0.0
	  DO j=1,mx_mz
	    IF(den_iz(i,j).gt.c_den) THEN
            ds=den_iz(i,j)
	      IF(j.gt.m_z) m_z=j
	    ENDIF
	  ENDDO
	  IF(ds.gt.c_den) m_i=m_i+1
	ENDDO
!  grt_i(i)-temperature gradient of i (keV/rho)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
      IF(default.EQ.0) THEN
      DO i=1,m_i
        grt_i(i)=-temp_i(i)/a0/2.0
        DO j=1,m_z
          grp_iz(i,j)=-temp_i(i)*den_iz(i,j)/a0
        ENDDO   
      ENDDO
      ENDIF
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
      c_potb=e0*bt0/2/q0**2
!  c_potl-q(0)*R(0) (m)
      c_potl=q0*r0
!  p_b2-<B**2> (T**2)
      p_b2=bt0**2*(1.0+0.5*p_eps**2)
!  p_bm2-<1/B**2> (/T**2)
      p_bm2=(1.0+1.5*p_eps**2)/bt0**2
!  p_eb-<E.B> (V*T/m)
      p_eb=0.1*bt0/(2.0*z_pi*r0)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
      p_fhat=p_q/p_eps
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
      DO i=1,3
        p_fm(i)=0.0
	ENDDO
      IF(p_eps.gt.0.0) THEN
        DO i=1,3
          p_fm(i)=i*((1.0-SQRT(1.0-p_eps**2))/p_eps)**(2.0*i)
     #            *(1.0+i*SQRT(1.0-p_eps**2))/((1.0-p_eps**2)**1.5
     #            *(p_q*r0)**2)
        ENDDO   
      ENDIF
!  p_ft-trapped fraction (-)
      p_ft=1.46*SQRT(p_eps)
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
      p_grbm2=1.0/bt0**2
!  p_ngrth-<n.grad(Theta)> (1/m)
      p_ngrth=1.0/(p_q*r0)
!Evaluate friction and viscosity coefficients and flows
      CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     #            p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi,
     #            p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,m_s,
     #            jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii,
     #            capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     #            utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss,
     #            dp_ss,dt_ss,iflag)
C
      RETURN
      END
