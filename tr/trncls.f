C  
C     ***********************************************************
C
C            NCLASS
C
C     ***********************************************************
C
      SUBROUTINE TR_NCLASS
!***********************************************************************
!TR_NCLASS calculates various paramters and arrays for NCLASS.
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
      INCLUDE 'trcomm.h'
      INCLUDE 'nclass/pamx_mi.inc'
      INCLUDE 'nclass/pamx_ms.inc'
      INCLUDE 'nclass/pamx_mz.inc'
      INCLUDE 'trncls.inc'
      REAL p_grstr(NRM),p_gr2str(NRM),p_grrm(NRM)
C
      k_order=2
      k_potato=0
      m_i=NSMAX
      DO NS=1,NSMAX
         PZMAX=PZ(NS)
         IF(PZ(NS).GE.PZMAX) PZMAX=PZ(NS)
      ENDDO
      m_z=INT(PZMAX)
      c_den=1.E10
      c_potb=SNGL(RKAP*BB/(2.D0*Q0**2))
      c_potl=SNGL(Q0*RR)
C
      DO NS=1,NSMAX
         amu_i(NS)=SNGL(PA(NS))
      ENDDO
C
      IF(NSMAX.EQ.2) THEN
         NS=NSMAX
         DO NR=1,NRMAX-1
            p_grphi=SNGL((RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))
     &             /DR)*SNGL(RKEV)
     &             /(SNGL(PZ(NS)*AEE*0.5D0*(RN(NR+1,NS)+RN(NR,NS))))
            p_grstr(NR)=p_grphi
         ENDDO
         NR=NRMAX
         p_grphi=SNGL(2.D0*(PNSS(NS)*PTS(NS)-RN(NR,NS)*RT(NR,NS))
     &          /DR)*SNGL(RKEV)
     &          /(SNGL(PZ(NS)*AEE*PNSS(NS)))
         p_grstr(NR)=p_grphi
c$$$         p_gr2str(1)=(p_grstr(2)-0.E0)/DR
c$$$         DO NR=2,NRMAX-1
c$$$            p_gr2str(NR)=(p_grstr(NR+1)-p_grstr(NR-1))/DR
c$$$         ENDDO
c$$$         p_gr2str(NRMAX)=(0.E0-p_grstr(NRMAX-1))/DR
C
         DO NR=2,NRMAX
            p_grrm(NR)=0.5E0*(p_grstr(NR  )/SNGL(BP(NR  ))
     &                       +p_grstr(NR-1)/SNGL(BP(NR-1)))
         ENDDO
         p_grrm(1)=2.E0*p_grrm(2)-p_grrm(3)
         DO NR=1,NRMAX-1
            p_gr2str(NR)=SNGL(BP(NR))*(p_grrm(NR+1)-p_grrm(NR))/SNGL(DR)
         ENDDO
         p_gr2str(NRMAX)=2.E0*p_gr2str(NRMAX-1)-p_gr2str(NRMAX-2)
      ELSE
         DO NR=1,NRMAX
            p_grstr(NR)=0.E0
            p_gr2str(NR)=0.E0
         ENDDO
      ENDIF
C
      DO NR=1,NRMAX
         EPS=EPSRHO(NR)
C
         p_b2=SNGL(BB**2*(1.D0+0.5D0*EPS**2))
         p_bm2=SNGL((1.D0+1.5D0*EPS**2)/BB**2)
         p_fhat=SNGL(QP(NR)/EPS*AR1RHO(NR))
         DO i=1,3
            p_fm(i)=0.0
         ENDDO
         IF(SNGL(EPS).gt.0.0) THEN
            DO i=1,3
               p_fm(i)=SNGL(DBLE(i)*((1.D0-SQRT(1.D0-EPS**2))/EPS)
     &                 **(2*i)
     &                 *(1.D0+DBLE(i)*SQRT(1.D0-EPS**2))/((1.D0-EPS**2)
     &                 **1.5*(QP(NR)*RR)**2))
            ENDDO   
         ENDIF
         p_ft1=SNGL(1.D0-(1.D0-EPS)**2.D0/(DSQRT(1.D0-EPS**2)
     &        *(1.D0+1.46D0*DSQRT(EPS))))
         p_ft2=SNGL(1.46D0*SQRT(EPS))
         p_ft3=SNGL(0.75D0*(1.5D0*SQRT(EPS))
     &        +0.25D0*(3.D0*SQRT(2.D0)/PI*SQRT(EPS)))
         p_ft=p_ft2
C         p_grbm2=SNGL(AR2RHO(NR)/BB**2)
C         p_grbm2=SNGL(AR2RHO(NR))*p_bm2
         p_grbm2=SNGL(AR2RHO(NR))/p_b2
         p_grphi=p_grstr(NR)
         p_gr2phi=p_gr2str(NR)
         p_ngrth=1.0/SNGL(QP(NR)*RR)
         IF(NR.EQ.NRMAX) THEN
            DO NS=1,NSMAX
               temp_i(NS)=SNGL(PTS(NS))
               grt_i(NS)=SNGL(2.D0*(PTS(NS)-RN(NR,NS))/DR)
               den_iz(NS,INT(ABS(PZ(NS))))=SNGL(PNSS(NS))*1.E20
               grp_iz(NS,INT(ABS(PZ(NS))))
     &              =SNGL(2.D0*(PNSS(NS)*PTS(NS)-RN(NR,NS)*RT(NR,NS))
     &              /DR)*1.E20
               DO NA=1,3
                  fex_iz(NA,NS,INT(ABS(PZ(NS))))=0.E0
               ENDDO
               p_eb=SNGL((1.5D0*ETA(NR  )*AJOH(NR  )
     &                   -0.5D0*ETA(NR-1)*AJOH(NR-1))*BB)
            ENDDO
         ELSE
            DO NS=1,NSMAX
               temp_i(NS)=SNGL(0.5D0*(RT(NR+1,NS)+RT(NR,NS)))
               grt_i(NS)=SNGL((RT(NR+1,NS)-RT(NR,NS))/DR)
               den_iz(NS,INT(ABS(PZ(NS))))=SNGL(0.5D0*(RN(NR+1,NS)
     &                                                +RN(NR,NS)))*1.E20
               grp_iz(NS,INT(ABS(PZ(NS))))
     &             =SNGL((RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))
     &              /DR)*1.E20
               DO NA=1,3
                  fex_iz(NA,NS,INT(ABS(PZ(NS))))=0.E0
               ENDDO
               p_eb=SNGL(0.5D0*(ETA(NR+1)*AJOH(NR+1)+ETA(NR)*AJOH(NR))
     &                   *BB)
            ENDDO
         ENDIF
C
      CALL NCLASS(
C     Input
     &            k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     &            p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi,
     &            p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,
C     Output
     &            m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,
     &            capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,
     &            sqz_s,upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,
     &            chip_ss,chit_ss,dp_ss,dt_ss,iflag)
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
            DO NS1=1,NSMAX
               AKNCP(NR,NS,NS1)=DBLE(chip_ss(NS,NS1))/AR2RHO(NR)
               AKNCT(NR,NS,NS1)=DBLE(chit_ss(NS,NS1))/AR2RHO(NR)
               ADNCP(NR,NS,NS1)=DBLE(dp_ss(NS,NS1))/AR2RHO(NR)
               ADNCT(NR,NS,NS1)=DBLE(dt_ss(NS,NS1))/AR2RHO(NR)
c$$$               if(ns.eq.ns1.and.ns.eq.1) ! or ns.eq.2
c$$$     &              write(6,*) NR,ADNCP(NR,NS,NS1)+ADNCT(NR,NS,NS1),
c$$$     &                            AKNCP(NR,NS,NS1)+AKNCT(NR,NS,NS1)
            ENDDO
C            DO NS1=3,4
C               write(6,'(3I5,1P2E20.7)') 
C     &              NR,NS,NS1,dp_ss(NS,NS1),dt_ss(NS,NS1)
C            ENDDO
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
      RETURN
      END
