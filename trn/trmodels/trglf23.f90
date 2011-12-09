MODULE trglf23

  PRIVATE

  PUBLIC tr_glf23

CONTAINS

!     ***********************************************************

!            GLF23 Model

!     ***********************************************************

      SUBROUTINE tr_glf23(dtr_tb,vtr_tb)

!   *************************************************************
!     In case of jshoot=0 and sometimes jmm=0, zeroth arguments of
!         te_m, ti_m, ne_m and ni_m
!     require finite values
!     typically the same ones as first arguments,
!     and those of
!         rho and rmin_exp
!     require zero avoiding numerical error due to the absence
!     of the value. However, the value calculated by using
!     zeroth value of rmin_exp are not used in this case.

!     Input parameters except pressure gradients : half grid
!     Output transport coefficients              : grid
!   *************************************************************

      USE TRCOMM, ONLY : &
           BB, MDLD, &
           NRMAX, NSAMAX, PA, PI, PZ, Q0, QP, RA, RG, RKAP, &
           RN, RR, RT, RMU0, aee
      IMPLICIT NONE
!     Inputs
      INTEGER(4),DIMENSION(5)  :: itport_pt
      REAL(8),DIMENSION(0:Nrmax):: &
           te_m, ti_m, rne_m, rni_m, rns_m, angrotp_exp, egamma_exp, &
           rgamma_p_exp, vphi_m, vpar_m, vper_m, zeff_exp, rho, &
           rgradrho_exp, rgradrhosq_exp, rmin_exp, rmaj_exp, q_exp, &
           shat_exp, alpha_exp, elong_exp

!     Outputs
      REAL(8),DIMENSION(0:nrmax)   :: &
           diff_m, chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, &
           egamma_m, rgamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m
      REAL(8),DIMENSION(0:NRMAX,10)::egamma_d

!     Auxiliaries
      REAL(8),DIMENSION(0:NRMAX):: &
           zpte_m, zpti_m, zpne_m, zpni_m, qe0, qi0, qn0, ddnn, ddne, ddni, &
           chien, chiee, chiei, chiin,chiie,chiii
      INTEGER(4):: &
           i_delay, idengrad, iglf, igrad, irotstab, j, jm, jmaxm, jmm, &
           jshoot, leigen, mode, nbt_flag, nr, nroot, ns, ns1
      REAL(8)   :: &
           alpha_e, amassgas_exp, amassimp_exp, arho_exp, bt_exp, chietem, &
           chiitim, deltat, diffnem, etaparm, etaperm, etaphim, exchm, &
           fctr, qe, qi, qn, rmajor_exp, x_alpha, zimp_exp, zpne_in, &
           zpni_in, zpte_in, zpti_in

      REAL(8),DIMENSION(3*nsamax,nrmax),INTENT(OUT):: dtr_tb
      REAL(8),DIMENSION(3*nsamax,3*nsamax,nrmax),INTENT(OUT):: vtr_tb
      REAL(8),DIMENSION(nrmax):: ar1rho,ar2rho,rkprho,rmjrho,rmnrho,dr,rm
      REAL(8),DIMENSION(nrmax):: wexb,wrot,vpar,vprp,vtor,zeff,alpha,s_hm
      REAL(8)   :: &
           CDH,ppp,ppm,dpp,phia,rkev
      INTEGER:: nsm,nsmax,mddw,mdluf,mdleqn,mdleqt,mdleoi,nsa,mdlkai

      mdlkai=mdld
      rkev=aee*1.D3
      phia=0.d0
      cdh=1.D0
      nsm=nsamax
      nsmax=nsamax
      mddw=0
      mdluf=0
      mdleqn=0
      mdleqt=1
      mdleoi=0
      DO nr=1,nrmax
         dr(nr)=rg(nr)-rg(nr-1)
         rm(nr)=0.5D0*(rg(nr)+rg(nr-1))
         ar1rho(nr)=1.D0/(sqrt(rkap)*ra)
         ar2rho(nr)=1.D0/(sqrt(rkap)*ra)**2
         rkprho(nr)=rkap
         rmjrho(nr)=rr
         rmnrho(nr)=ra*rg(nr)
         ppp=0.D0
         ppm=0.D0
         do nsa=1,nsamax
            ppp=ppp+rn(nsa,nr  )*rt(nsa,nr  )
            ppm=ppm+rn(nsa,nr-1)*rt(nsa,nr-1)
         end do
         dpp=(ppp-ppm)/dr(nr)
         alpha(nr)=-2.D0*RMU0*QP(NR)**2*RR/BB**2*(DPP*1.D20*RKEV)
         wexb(nr)=0.d0
         wrot(nr)=0.d0
         vpar(nr)=0.d0
         vprp(nr)=0.d0
         vtor(nr)=0.d0
         zeff(nr)=2.d0
      ENDDO

      DO nr=1,nrmax
         s_hm(nr)= (rg(nr)+rg(nr-1))*(qp(nr)-qp(nr-1)) &
                 /((rg(nr)-rg(nr-1))*(qp(nr)+qp(nr-1)))
      END DO

      MDDW=1
!     INPUTS
      leigen=1        ! 1 for tomsqz, 0 for cgg solver
      IF(MDLUF.NE.0.AND.NSMAX.GT.2) THEN
         nroot=12   ! num. of equations, 8 for default, 12 for imp.
      ELSE
         nroot=8
      ENDIF
      iglf=1          ! 0 for original model, 1 for retuned model
      jshoot=0        ! 0 for time-dep code, 1 for shooting code
!   In case of jshoot=0,
!      maximum argument of array is important;
!      values of zeroth argument are not important, but avoiding
!      Inf error due to logarithm calculation some finite values
!      need to be stored;
!      NR=1 to NRMAX corresponds to jm=1 to jmaxm in callglf2d.f.
!
      jmm=0           ! jmm=0 does full grid from jm=1,jmaxm-1
      jmaxm=NRMAX     ! maximum num. of grid points
      itport_pt(1)=MDLEQN  ! density transport
      IF(MDLEQT.NE.0) THEN
         IF(NSMAX.EQ.1) THEN
            IF(MDLEOI.EQ.1) THEN
               itport_pt(2)=MDLEQT ! electron transport
               itport_pt(3)=0      ! ion transport
            ELSEIF(MDLEOI.EQ.2) THEN
               itport_pt(2)=0
               itport_pt(3)=MDLEQT
            ELSE
               itport_pt(2)=MDLEQT
               itport_pt(3)=MDLEQT
            ENDIF
         ELSE
            itport_pt(2)=MDLEQT
            itport_pt(3)=MDLEQT
         ENDIF
      ELSE
         itport_pt(2)=MDLEQT
         itport_pt(3)=MDLEQT
      ENDIF
      itport_pt(4)=0  ! v_phi transport
      itport_pt(5)=0  ! v_theta transport
      irotstab=1      ! 1 uses internally computed ExB shear

      jm=0
         te_m(jm) =RT (1,jm+1)
         ti_m(jm) =RT (2,jm+1)
         rne_m(jm)=RN (1,jm+1)*10.D0
         rni_m(jm)=RN (2,jm+1)*10.D0
         rns_m(jm)=0.D0
      DO jm=1,jmaxm
         te_m(jm) =RT (1,jm)        ! Te [keV] ! halfmesh
         ti_m(jm) =RT (2,jm)        ! Ti [keV]
         rne_m(jm)=RN (1,jm)*10.D0  ! Ne [^19m^-3]
         rni_m(jm)=RN (2,jm)*10.D0  ! Ni [^19m^-3]
         rns_m(jm)=0.D0
      ENDDO

      IF(MDLUF.NE.0.AND.NSMAX.GT.2) THEN
         idengrad=3   ! compute simple dilution (3=actual dilution)
      ELSE
         idengrad=2
      ENDIF

      jm=1
         angrotp_exp(jm)  =FCTR(RG(jm),RG(jm+1),WROT(jm),WROT(jm+1))
         egamma_exp(jm)  = FCTR(RG(jm),RG(jm+1),WEXB(jm),WEXB(jm+1))
         rgamma_p_exp(jm)= 0.D0
         vphi_m(jm)      = FCTR(RG(jm),RG(jm+1),VTOR(jm),VTOR(jm+1))
         vpar_m(jm)      = FCTR(RG(jm),RG(jm+1),VPAR(jm),VPAR(jm+1))
         vper_m(jm)      = FCTR(RG(jm),RG(jm+1),VPAR(jm),VPAR(jm+1))
         zeff_exp(jm)    = ZEFF(jm)
      DO jm=2,jmaxm
!     exp. toroidal ang. vel. [1/s]
         angrotp_exp(jm) = 0.5D0*(WROT(jm-1)+WROT(jm))
!     exp. ExB shearing rate: valid if irotstab=0
         egamma_exp(jm)  = 0.5D0*(WEXB(jm-1)+WEXB(jm))
!     exp. para. vel. shearing rate
!       : valid if itport_pt(4)=-1 or irotstab=0
         rgamma_p_exp(jm)= 0.D0
!     toroidal velocity [m/s]
!       : valid if finite itport_pt(4) and itport_pt(5)=0
         vphi_m(jm)      = 0.5D0*(VTOR(jm-1)+VTOR(jm))
!     parallel velocity [m/s]: valid if finite itport_pt(4&5)
         vpar_m(jm)      = 0.5D0*(VPAR(jm-1)+VPAR(jm))
!     perpendicular velocity [m/s]: valid if finite itport_pt(5)
         vper_m(jm)      = 0.5D0*(VPRP(jm-1)+VPRP(jm))
!     effective charge
         zeff_exp(jm)    = ZEFF(jm)
      ENDDO

      bt_exp=BB      ! toroidal field [T]
      nbt_flag=1     ! >0 for Beff, Bt otherwise

!     normalized flux surface; 0 < rho < 1.
      rho(0)=0.D0
      DO jm=1,jmaxm
         rho(jm)=RM(jm) ! norm. toroidal flux surf. label
      ENDDO

      IF(PHIA.EQ.0.D0) THEN
         arho_exp=SQRT(RKAP)*RA ! rho at last closed flux surface [m]
      ELSE
         arho_exp=SQRT(PHIA/(PI*bt_exp))
      ENDIF

      rmin_exp(0)=0.D0
      CALL AITKEN(RM(1),rmin_exp(1),RG,RMNRHO,1,NRMAX)
      CALL AITKEN(RM(1),rmaj_exp(1),RG,RMJRHO,1,NRMAX)
      rgradrho_exp(1)  =AR1RHO(1)*arho_exp
      rgradrhosq_exp(1)=AR2RHO(1)*arho_exp**2
      DO jm=2,jmaxm
         rgradrho_exp(jm)  =AR1RHO(jm)*arho_exp
         rgradrhosq_exp(jm)=AR2RHO(jm)*arho_exp**2
         rmin_exp(jm)      =0.5D0*(RMNRHO(jm-1)+RMNRHO(jm))
                                       ! local minor radius [m]
         rmaj_exp(jm)      =0.5D0*(RMJRHO(jm-1)+RMJRHO(jm))
                                       ! local major radius [m]
      ENDDO
      rmajor_exp=RR  ! geometrical major radius of magnetix axis [m]

      zimp_exp=PZ(3)         ! Zimp; finite data is necessary
      amassimp_exp=PA(3)     ! Aimp; finite data is necessary

      q_exp(1)=0.5D0*(Q0+QP(1))
      DO jm=2,jmaxm
         q_exp(jm)=0.5D0*(QP(jm-1)+QP(jm))  ! safety factor
      ENDDO

      shat_exp (1)=S_HM(1)
      alpha_exp(1)=FCTR(RG(1),RG(2),ALPHA(1),ALPHA(2))
      elong_exp(1)=RKPRHO(1)
      DO jm=2,jmaxm
         shat_exp (jm)=S_HM(jm)      ! magnetic shear
         alpha_exp(jm)=0.5D0*(ALPHA(jm-1)+ALPHA(jm)) ! MHD alpha
         elong_exp(jm)=RKPRHO(jm)    ! local elongation
      ENDDO

      amassgas_exp=PA(2) ! atomic num. of working gas
      alpha_e=1.D0       ! ExB shear stabilization (0=off,>0=on)
      x_alpha=1.D0       ! alpha stabilization (0=off,>0=on)
      i_delay=0          ! default(usually recommended)

      IF(MDLKAI.EQ.60) THEN
!     +++ Normal type +++

         igrad=0         ! compute gradients (1=input gradients)
         zpte_in=0.D0    ! 1/Lte (necessary if igrad and jmm != 0)
         zpti_in=0.D0    ! 1/Lti (necessary if igrad and jmm != 0)
         zpne_in=0.D0    ! 1/Lne (necessary if igrad and jmm != 0)
         zpni_in=0.D0    ! 1/Lni (necessary if igrad and jmm != 0)

         call callglf2d( leigen, nroot, iglf, jshoot, jmm, jmaxm, itport_pt, &
              irotstab, te_m, ti_m, rne_m, rni_m, rns_m, igrad, idengrad, &
              zpte_in, zpti_in, zpne_in, zpni_in, angrotp_exp, egamma_exp, &
              rgamma_p_exp, vphi_m, vpar_m, vper_m, zeff_exp, bt_exp, &
              nbt_flag, rho, arho_exp, rgradrho_exp, rgradrhosq_exp, &
              rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp, &
              q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, &
              x_alpha, i_delay, diffnem, chietem, chiitim, etaphim, etaparm, &
              etaperm, exchm, diff_m, chie_m, chii_m, etaphi_m, etapar_m, &
              etaper_m, exch_m, egamma_m, egamma_d, rgamma_p_m, anrate_m, &
              anrate2_m, anfreq_m, anfreq2_m )

         dtr_tb(1:3*nsamax,1:3*nsamax,1:nrmax)=0.D0
         vtr_tb(1:3*nsamax,1:nrmax)=0.D0
         DO NR=1,NRMAX
            DO nsa=1,nsamax
               ns=ns_nsa(nsa)
               base=3*(nsa-1)
               IF(pz0(ns) < 0.D0) THEN ! for electron
                  dtr_tb(base+1,base+1,nr)=MAX(diff_m(NR),0.D0)
                  dtr_tb(base+3,base+3,nr)=MAX(chie_m(NR),0.D0)
               ELSE IF
                  IF(pz(ns) /= 0.d0) THEN ! for ion
                     dtr_tb(base+1,base+1,nr)=MAX(diff_m(NR),0.D0)
                     dtr_tb(base+3,base+3,nr)=MAX(chii_m(NR),0.d0)
                  END IF
               END IF
            END DO
         END DO

      ELSEIF(MDLKAI.EQ.61) THEN
!     +++ D-V (diffusion convection) method +++

         igrad=1  ! compute gradients (1=input gradients)
         deltat=0.03D0
         MODE=0   ! 0: Kinsey type, 1: Angioni type
         do j=1,jmaxm-1
            zpte_m(j)=-(LOG(RT(j+1,1))-LOG(RT(j,1)))/DR(j) ! grid
            zpti_m(j)=-(LOG(RT(j+1,2))-LOG(RT(j,2)))/DR(j) ! grid
            zpne_m(j)=-(LOG(RN(j+1,1))-LOG(RN(j,1)))/DR(j) ! grid
            zpni_m(j)=-(LOG(RN(j+1,2))-LOG(RN(j,2)))/DR(j) ! grid
            zpte_in=zpte_m(j)
            zpti_in=zpti_m(j)
            zpne_in=zpne_m(j)
            zpni_in=zpni_m(j)

            jmm=j
            call callglf2d( leigen, nroot, iglf, jshoot, jmm, jmaxm, itport_pt &
     & , irotstab, te_m, ti_m, rne_m, rni_m, rns_m, igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in &
     & , angrotp_exp, egamma_exp, rgamma_p_exp, vphi_m, vpar_m, vper_m, zeff_exp, bt_exp, nbt_flag, rho &
     & , arho_exp, rgradrho_exp, rgradrhosq_exp, rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp &
     & , q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha, i_delay &
     & , diffnem, chietem, chiitim, etaphim, etaparm, etaperm, exchm, diff_m, chie_m, chii_m, etaphi_m &
     & , etapar_m, etaper_m, exch_m, egamma_m, egamma_d, rgamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )

            qe0(j)=chietem*zpte_in
            qi0(j)=chiitim*zpti_in
            qn0(j)=diffnem*zpni_in

!     1/Lt1=1/Lt0+DELTAt*1/Lt0

            zpte_in=zpte_m(j)*(1.D0+deltat)
            zpti_in=zpti_m(j)*(1.D0+deltat)
            zpne_in=zpne_m(j)*(1.D0+deltat)
            zpni_in=zpni_m(j)*(1.D0+deltat)

            call callglf2d( leigen, nroot, iglf, jshoot, jmm, jmaxm, itport_pt &
     & , irotstab, te_m, ti_m, rne_m, rni_m, rns_m, igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in &
     & , angrotp_exp, egamma_exp, rgamma_p_exp, vphi_m, vpar_m, vper_m, zeff_exp, bt_exp, nbt_flag, rho &
     & , arho_exp, rgradrho_exp, rgradrhosq_exp, rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp &
     & , q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha, i_delay &
     & , diffnem, chietem, chiitim, etaphim, etaparm, etaperm, exchm, diff_m, chie_m, chii_m, etaphi_m &
     & , etapar_m, etaper_m, exch_m, egamma_m, egamma_d, rgamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )

            qe=chietem*zpte_in
            qi=chiitim*zpti_in
            qn=diffnem*zpni_in

            chien(j)=(qe-qe0(j))/(zpni_in-zpni_m(j))
            chiee(j)=(qe-qe0(j))/(zpte_in-zpte_m(j))
            chiei(j)=(qe-qe0(j))/(zpti_in-zpti_m(j))
            chiin(j)=(qi-qi0(j))/(zpni_in-zpni_m(j))
            chiie(j)=(qi-qi0(j))/(zpte_in-zpte_m(j))
            chiii(j)=(qi-qi0(j))/(zpti_in-zpti_m(j))
            ddnn (j)=(qn-qn0(j))/(zpni_in-zpni_m(j))
            ddne (j)=(qn-qn0(j))/(zpte_in-zpte_m(j))
            ddni (j)=(qn-qn0(j))/(zpti_in-zpti_m(j))

            IF(MODE.EQ.0) THEN
               AKDW (j,1)  =chiee(j)
               AKDWP(j,1,1)=chiee(j)
               AVKDW(j,1)  =(qe0(j)/zpte_m(j)-chiee(j))*zpte_m(j)
               IF(chiee(j).LT.0.D0) THEN
                  AKDW(j,1)   =0.D0
                  AKDWP(j,1,1)=0.D0
                  AVKDW(j,1)  =qe0(j)
               ENDIF
               DO NS=2,NSM
                  AKDW (j,NS)   =chiii(j)
                  AKDWP(j,NS,NS)=chiii(j)
                  AVKDW(j,NS)   =(qi0(j)/zpti_m(j)-chiii(j))*zpti_m(j)
                  IF(chiii(j).LT.0.D0) THEN
                     AKDW (j,NS)   =0.D0
                     AKDWP(j,NS,NS)=0.D0
                     AVKDW(j,NS)   =qi0(j)
                  ENDIF
               ENDDO
            ENDIF
         enddo
         qe0(jmaxm)=0.D0
         qi0(jmaxm)=0.D0
         qn0(jmaxm)=0.D0

!     Let AVDW enable by turning on CDH for anomalous particle convection
         CDH=1.D0

         dtr_tb(1:3*nsamax,1:3*nsamax,1:nrmax)=0.D0
         vtr_tb(1:3*nsamax,1:nrmax)=0.D0

         IF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
            DO nsa=1,nsamax
               ns=ns_nsa(nsa)
               base=3*(nsa-1)
               IF(pz0(ns) < 0.D0) THEN ! for electron
                  vtr_tb(base+1,nr)=qn0(NR)-ddnn (NR)*zpni_m(NR) &
                                   -ddne (NR)*zpte_m(NR) &
                                   -ddni (NR)*zpti_m(NR)
                  vtr_tb(base+3,nr)=qe0(NR)-chien(NR)*zpni_m(NR) &
                                   -chiee(NR)*zpte_m(NR) &
                                   -chiei(NR)*zpti_m(NR)
                  DO nsa1=1,nsamax
                     ns1=ns_nsa(nsa1)
                     base1=3*(nsa1-1)
                     IF(pz0(ns1) < 0.D0) THEN ! for electron
                        dtr_tb(base+1,base1+1,nr)=MAX(ddnn(NR)+ddne(NR),0.D0)
                        dtr_tb(base+1,base1+3,nr)=MAX(         ddne(NR),0.D0)
                        dtr_tb(base+3,base1+1,nr)=MAX(chien(NR)+chiee(NR),0.D0)
                        dtr_tb(base+3,base1+3,nr)=MAX(          chiee(NR),0.D0)
                     ELSE IF(pz(ns1) /= 0.d0) THEN ! for ion
                        dtr_tb(base+1,base1+1,nr)=MAX(ddnn(NR)+ddne(NR),0.D0)
                        dtr_tb(base+1,base1+3,nr)=MAX(         ddne(NR),0.D0)
                        dtr_tb(base+3,base1+1,nr)=MAX(chien(NR)+chiee(NR),0.D0)
                        dtr_tb(base+3,base1+3,nr)=MAX(          chiee(NR),0.D0)
                     END IF
                  END DO
               ELSE IF(pz(ns) /= 0.d0) THEN ! for ion
                  vtr_tb(base+1,nr)=qn0(NR)-ddnn (NR)*zpni_m(NR) &
                                           -ddne (NR)*zpte_m(NR) &
                                           -ddni (NR)*zpti_m(NR)
                  vtr_tb(base+3,nr)=qi0(NR)-chiin(NR)*zpni_m(NR) &
                                           -chiie(NR)*zpte_m(NR) &
                                           -chiii(NR)*zpti_m(NR)
                  DO nsa1=1,nsamax
                     IF(pz0(ns1) < 0.D0) THEN ! for electron
                        dtr_tb(base+1,base1+1,nr)=MAX(ddnn(NR)+ddni(NR),0.D0)
                        dtr_tb(base+1,base1+3,nr)=MAX(         ddni(NR),0.D0)
                        dtr_tb(base+3,base1+1,nr)=MAX(chien(NR)+chiei(NR),0.D0)
                        dtr_tb(base+3,base1+3,nr)=MAX(          chiei(NR),0.D0)
                     ELSE IF(pz(ns1) /= 0.d0) THEN ! for ion
                        dtr_tb(base+1,base1+1,nr)=MAX(ddnn(NR)+ddni(NR),0.D0)
                        dtr_tb(base+1,base1+3,nr)=MAX(         ddni(NR),0.D0)
                        dtr_tb(base+3,base1+1,nr)=MAX(chiin(NR)+chiii(NR),0.D0)
                        dtr_tb(base+3,base1+3,nr)=MAX(          chiii(NR),0.D0)
                     END IF
                  END DO
               END IF
            END DO
         END DO

         IF(MODE.EQ.1) THEN
            DO NR=1,NRMAX-1
               ADDWD(NR,1,1)=ddnn(NR)
               ADDWP(NR,1,1)=ddne(NR)
               AKDWD(NR,1,1)=chien(NR)
               AKDWP(NR,1,1)=chiee(NR)
               DO NS=2,NSM
                  ADDWD(NR,NS,1)=ddnn(NR)
                  ADDWP(NR,NS,1)=ddni(NR)
                  AKDWD(NR,NS,1)=chien(NR)
                  AKDWP(NR,NS,1)=chiei(NR)
               ENDDO
               DO NS1=2,NSM
                  ADDWD(NR,1,NS1)=ddnn(NR)
                  ADDWP(NR,1,NS1)=ddne(NR)
                  AKDWD(NR,1,NS1)=chiin(NR)
                  AKDWP(NR,1,NS1)=chiie(NR)
                  DO NS=2,NSM
                     ADDWD(NR,NS,NS1)=ddnn(NR)
                     ADDWP(NR,NS,NS1)=ddni(NR)
                     AKDWD(NR,NS,NS1)=chiin(NR)
                     AKDWP(NR,NS,NS1)=chiii(NR)
                  ENDDO
               ENDDO
               DO NS=1,NSM
                  ADDW(NR,NS)=ADDWD(NR,NS,NS)
                  AKDW(NR,NS)=AKDWP(NR,NS,NS)
               ENDDO
            ENDDO
            NR=NRMAX
               ADDWD(NR,1:NSM,1:NSM)=ADDWD(NR-1,1:NSM,1:NSM)
               ADDWP(NR,1:NSM,1:NSM)=ADDWP(NR-1,1:NSM,1:NSM)
               AKDWD(NR,1:NSM,1:NSM)=AKDWD(NR-1,1:NSM,1:NSM)
               AKDWP(NR,1:NSM,1:NSM)=AKDWP(NR-1,1:NSM,1:NSM)
               DO NS=1,NSM
                  ADDW(NR,NS)=ADDWD(NR,NS,NS)
                  AKDW(NR,NS)=AKDWP(NR,NS,NS)
               ENDDO

            DO NR=1,NRMAX
               AVKDW(NR, 1)=qe0(NR)-chien(NR)*zpni_m(NR) &
                                   -chiee(NR)*zpte_m(NR) &
                                   -chiei(NR)*zpti_m(NR)
               AVDW (NR, 1)=qn0(NR)-ddnn (NR)*zpni_m(NR) &
                                   -ddne (NR)*zpte_m(NR) &
                                   -ddni (NR)*zpti_m(NR)
               DO NS=2,NSM
                  AVKDW(NR,NS)=qi0(NR)-chiin(NR)*zpni_m(NR) &
                                      -chiie(NR)*zpte_m(NR) &
                                      -chiii(NR)*zpti_m(NR)
                  AVDW (NR,NS)=qn0(NR)-ddnn (NR)*zpni_m(NR) &
                                      -ddne (NR)*zpte_m(NR) &
                                      -ddni (NR)*zpti_m(NR)
               ENDDO
!               write(6,'(I3,6F12.7)') NR,chiin(NR),zpni_m(NR)
!     &              ,chiie(NR),zpte_m(NR),chiii(NR),zpti_m(NR)
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE tr_glf23
    END MODULE trglf23
    
