!     ***********************************************************

!            GLF23 Model

!     ***********************************************************

      SUBROUTINE GLF23_DRIVER(S_HM)

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

      USE TRCOMM
      USE libitp
      IMPLICIT NONE
      INCLUDE '../trmlib/glf/trglf.inc'
      REAL(rkind),DIMENSION(NRMAX), INTENT(IN):: S_HM
      INTEGER:: &
           i_delay, idengrad, iglf, igrad, irotstab, j, jm, jmaxm, jmm, &
           jshoot, leigen, mode, nbt_flag, nr, nroot, ns, ns1
      REAL(rkind)   :: &
           alpha_e, amassgas_exp, amassimp_exp, arho_exp, bt_exp, chietem, &
           chiitim, deltat, diffnem, etaparm, etaperm, etaphim, exchm, &
           qe, qi, qn, rmajor_exp, x_alpha, zimp_exp, zpne_in, &
           zpni_in, zpte_in, zpti_in

      MDLDW=1
!     INPUTS
      leigen=1        ! 1 for tomsqz, 0 for cgg solver
      nroot=8
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
         te_m(jm) =RT (jm+1,1)
         ti_m(jm) =RT (jm+1,2)
         rne_m(jm)=RN (jm+1,1)*10.D0
         rni_m(jm)=RN (jm+1,2)*10.D0
         rns_m(jm)=RNF(jm+1,1)*10.D0
      DO jm=1,jmaxm
         te_m(jm) =RT (jm,1)        ! Te [keV] ! halfmesh
         ti_m(jm) =RT (jm,2)        ! Ti [keV]
         rne_m(jm)=RN (jm,1)*10.D0  ! Ne [^19m^-3]
         rni_m(jm)=RN (jm,2)*10.D0  ! Ni [^19m^-3]
         rns_m(jm)=RNF(jm,1)*10.D0  ! Nf [^19m^-3]
      ENDDO

      idengrad=2

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

      arho_exp=SQRT(RKAP)*RA ! rho at last closed flux surface [m]

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

         call callglf2d( leigen, nroot, iglf, jshoot, jmm, jmaxm, itport_pt &
     & , irotstab, te_m, ti_m, rne_m, rni_m, rns_m, igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in &
     & , angrotp_exp, egamma_exp, rgamma_p_exp, vphi_m, vpar_m, vper_m, zeff_exp, bt_exp, nbt_flag, rho &
     & , arho_exp, rgradrho_exp, rgradrhosq_exp, rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp &
     & , q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp, alpha_e, x_alpha, i_delay &
     & , diffnem, chietem, chiitim, etaphim, etaparm, etaperm, exchm, diff_m, chie_m, chii_m, etaphi_m, etapar_m &
     & , etaper_m, exch_m, egamma_m, egamma_d, rgamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m )

         DO NR=1,NRMAX-1
            AKDW(NR,1)=chie_m(NR)
            AKDW(NR,2)=chii_m(NR)
            AKDW(NR,3)=chii_m(NR)
            AKDW(NR,4)=chii_m(NR)
            ADDW(NR,1)=diff_m(NR)
            ADDW(NR,2)=diff_m(NR)
            ADDW(NR,3)=diff_m(NR)
            ADDW(NR,4)=diff_m(NR)
!            write(6,*) "chii_m(",NR,")=",chii_m(NR)
         ENDDO
         NR=NRMAX
            AKDW(NR,1)=chie_m(NR-1)
            AKDW(NR,2)=chii_m(NR-1)
            AKDW(NR,3)=chii_m(NR-1)
            AKDW(NR,4)=chii_m(NR-1)
            ADDW(NR,1)=diff_m(NR-1)
            ADDW(NR,2)=diff_m(NR-1)
            ADDW(NR,3)=diff_m(NR-1)
            ADDW(NR,4)=diff_m(NR-1)
!            write(6,*) "chii_m(",NR,")=",chii_m(NR)

         DO NR=1,NRMAX
            DO NS=1,4
               IF(AKDW(NR,NS).LT.0.D0) THEN
                  AKDW(NR,NS)=0.D0
               ENDIF
            ENDDO
         ENDDO

      ELSEIF(MDLKAI.EQ.61) THEN
!     +++ D-V (diffusion convection) method +++

         igrad=1  ! compute gradients (1=input gradients)
         deltat=0.03D0
         MODE=0   ! 0: Kinsey type, 1: Angioni type
         do j=1,jmaxm-1
            zpte_m(j)=-(LOG(RT(j+1,1))-LOG(RT(j,1)))/DR ! grid
            zpti_m(j)=-(LOG(RT(j+1,2))-LOG(RT(j,2)))/DR ! grid
            zpne_m(j)=-(LOG(RN(j+1,1))-LOG(RN(j,1)))/DR ! grid
            zpni_m(j)=-(LOG(RN(j+1,2))-LOG(RN(j,2)))/DR ! grid
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
      END SUBROUTINE GLF23_DRIVER

!     ***********************************************************

!            Weiland Model

!     ***********************************************************

      SUBROUTINE WEILAND_DRIVER

!***********************************************************************
!  <INPUT>
!     ENL    : 2 Lne/Lb (Lb=R) (epsilon_ne)
!     EIL    : Lni/Lti (eta_i)
!     EEL    : Lne/Lte (eta_e)
!     TAUL   : Te/Ti (tau_i)
!     FLL    : Finite Larmor Radius parameter: (k_perp*rhos)**2=0.1
!     FTL    : trapped particle fraction
!     BQL    : fraction of an impurity
!     EQL    : Lnq/Ltq (q represents an impurity) (eta_q)
!     ENQL   : 2 Lnq/Lb (epsilon_nq)
!     ZL     : Z value for an impurity
!     BETAEL : electron beta
!     AZL    : impurity mass number
!     COLL   : factor multiplying collisionality (0 or 1)
!     ELL    : factor multiplying electromagnetic effects (0 or 1)
!     TEL    : electron temperature [keV]
!     TAUZL  : Te/Tq (tau_q)
!     RA     : minor radius (a)[m]
!     QL     : safety factor
!     SL     : magnetic shear
!     EPS    : local inverse aspect ratio (r/R)
!     EPSA   : inverse aspect ratio (a/R)
!     RNL    : electron density [10^19(20?) m^-3]
!     RLIST  : controlling printout in disp9t (0: off)
!     RNEQL  : number of equations (NDISP)
!     RKAP   : ellipticity
!     RIWL   : controlling printout in DIFFTD (0: off)
!     RISBL  : If ISB=1, we use the strong ballooning approximation
!              (GAV=1).
!     BB     : toroidal magnetic field [T]
!     SEARCH : The way the code chooses between ion and electron
!              eigenvalues for the use in the eigenfunction (WZ).
!              1 : Only eigenvalues with negative real part are used.
!                 (Real(WZ)<0)
!              2 : Eigenvalues with positive real part are used unless
!                  there are eigenvalues with negative real part.
!                 (Real(WZ)>0)
!              3 : The fastest growing mode with positive real part is
!                  used for eigenvalues with positive real part
!                  if there is such a root.
!     PMA    : mass number for an main ion
!     RGKL   : GRHO1/(a*GRHO2)
!     WEXBL  : ExB shearing rate
!     ROTL   : factor multiplying WEXBL (0 or 1)
!     NR     : radial grid number
!     IST    : 1   : iterations starting from an analytical
!                    approximation
!                    You should use only the first time step.
!              else: iterations starting from the value stored
!                    in WZJ(IK) which is the eigenvalue from the
!                    previous time step

!  <OUTPUT>
!     CHIL(5) : ion thermal transport coefficients for Ti, Te, Ne, Tq
!               and Nq equations [m^2/s (same as above)]
!     CHEL(5) : electron thermal transport coefficients
!     DL(5)   : ion particle transport coefficients
!     CHQL(5) : impurity thermal transport coefficients
!     DQL(5)  : impurity particle transport coefficients
!     SCHI    : effective ion thermal transport coefficient
!     SCHE    : effective electron thermal transport coefficient
!     SD      : effective ion particle transport coefficient
!     SCHQ    : effective impurity thermal transport coefficient
!     SDQ     : effective impurity particle transport coefficient

!***********************************************************************

      USE TRCOMM, ONLY : &
           AR1RHOG, AR2RHOG, BB, DR, EPSRHO, MDLDW, MDLKAI, MDLTPF, &
           NRMAX, NT, PA, PNSS, PTS, PZ, QP, RA, RHOG, RHOM, RJCB, RKAP, &
           RKEV, RMU0, RN, RR, RT, WEXB, S, rkind
      USE libitp
      IMPLICIT NONE
      INTEGER:: ist, nr
      REAL(rkind)   :: &
           azl, betael, bql, coll, drl, eel, eil, ell, enl, enql, eps, &
           epsa, eql, fll, fls, ftl, ftpf, pma, ql, rgkl, risbl, riwl, &
           rlist, rneql, rnl, rotl, sche, schi, schq, sd, sdq, search, &
           shat, sl, slbl, slnel, slnil, slnql, sltel, sltil, sltql, taul, &
           tauzl, tel, wexbl, zl
      REAL(rkind),DIMENSION(5):: CHEL, CHIL, CHQL, DL, DQL

      MDLDW=1
      IF(NT.EQ.0) THEN
         IST=1
      ELSE
         IST=0
      ENDIF
      ZL    = PZ(3)
      AZL   = PA(3)
      COLL  = 1.D0
      ELL   = 1.D0
      RLIST = 1.D0
      RNEQL = 9.D0
      RIWL  = 2.D0
      RISBL = 2.D0
      SEARCH= 2.D0
      PMA   = PA(2)
      ROTL  = 1.D0
      EPSA  = RA/RR
      DO NR=1,NRMAX-1
         DRL   = RJCB(NR)/DR
         EPS   = EPSRHO(NR)
         SLNEL =-0.5D0*(RN(NR+1,1)+RN(NR,1))/((RN(NR+1,1)-RN(NR,1))*DRL)
         SLNIL =-0.5D0*(RN(NR+1,2)+RN(NR,2))/((RN(NR+1,2)-RN(NR,2))*DRL)
         SLNQL =-0.5D0*(RN(NR+1,3)+RN(NR,3))/((RN(NR+1,3)-RN(NR,3))*DRL)
         SLTEL =-0.5D0*(RT(NR+1,1)+RT(NR,1))/((RT(NR+1,1)-RT(NR,1))*DRL)
         SLTIL =-0.5D0*(RT(NR+1,2)+RT(NR,2))/((RT(NR+1,2)-RT(NR,2))*DRL)
         SLTQL =-0.5D0*(RT(NR+1,3)+RT(NR,3))/((RT(NR+1,3)-RT(NR,3))*DRL)
         SLBL  = RR
         ENL   = 2.D0*(SLNEL/SLBL )
         EIL   =       SLNIL/SLTIL
         EEL   =       SLNEL/SLTEL
         TAUL  = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,2)+RT(NR,2))
         FLL   = 1.D-1
         FTL   = FTPF(MDLTPF,EPS)
         BQL   = (RN(NR+1,3)+RN(NR,3))/(RN(NR+1,1)+RN(NR,1))
         EQL   =       SLNQL/SLTQL
         ENQL  = 2.D0*(SLNQL/SLBL )
         BETAEL= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))*RKEV*1.D20 &
                 /(BB**2/(2.D0*RMU0))
         TEL   = 0.5D0*(RT(NR+1,1)+RT(NR,1))
         TAUZL = (RT(NR+1,1)+RT(NR,1))/(RT(NR+1,3)+RT(NR,3))
         QL    = QP(NR)
         SL    = S(NR)
         RNL   = 0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1

         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = WEXB(NR)
         IF(MDLKAI.EQ.64) THEN
            SHAT  = SQRT(2.D0*SL-1.D0+RKAP**2*(SL-1.D0)**2)
            FLS   = (0.7D0+2.4D0/(7.14D0*QL*SHAT+0.1D0))*FLL
            FLL   = 2.D0*FLS/(1.D0+1.D0/TAUL)
         ENDIF

!         COEF = PZ(2)**2*AEE**4*1.D20
!     &         /(6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*RKEV**1.5D0)
!         RLAMB =15.2D0-0.5D0*DLOG(RNL)+DLOG(TEL)
!         VEI  = COEF*RNL*RLAMB/TEL**1.5D0
!         write(6,*) NR,COEF,VEI

         CALL TR_WEILAND_BRIDGE (ENL,EIL,EEL,TAUL,FLL,FTL,BQL,EQL,ENQL,ZL, &
              BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,QL,SL,EPS,EPSA,RNL,RLIST, &
              RNEQL,RKAP,RIWL,RISBL,BB,SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST, &
              CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)

!         write(6,'(I3,5F15.7)') NR,CHIL(2),CHEL(2),DL(2),CHQL(2),DQL(2)

         CALL WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL) !,SCHI,SCHE,SD,SCHQ,SDQ)
!       write(6,'(I3,5F15.7)') NR,CHEL(1),CHEL(2),CHEL(3),CHEL(4),CHEL(5)
!       write(6,'(I3,5F15.7)') NR,CHIL(1),CHIL(2),CHIL(3),CHIL(4),CHIL(5)
!       write(6,'(I3,5F15.7)') NR,CHQL(1),CHQL(2),CHQL(3),CHQL(4),CHQL(5)
!       write(6,'(I3,5F15.7)') NR,DL(1),DL(2),DL(3),DL(4),DL(5)
!         if(nr.eq.4) write(6,'(I4,5F15.7)') NT,DL(3),CHQL(4)
      ENDDO

      NR=NRMAX
         EPS   = EPSRHO(NR)
         SLNEL =-PNSS(1)/DERIV3P(PNSS(1),RN(NR,1),RN(NR-1,1), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLNIL =-PNSS(2)/DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLNQL =-PNSS(3)/DERIV3P(PNSS(3),RN(NR,3),RN(NR-1,3), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLTEL =-PTS (1)/DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLTIL =-PTS (2)/DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLTQL =-PTS (3)/DERIV3P(PTS(3),RT(NR,3),RT(NR-1,3), &
                                 RHOG(NR),RHOM(NR),RHOM(NR-1))
         SLBL  = RR
         ENL   = 2.D0*(SLNEL/SLBL )
         EIL   =       SLNIL/SLTIL
         EEL   =       SLNEL/SLTEL
         TAUL  = PTS(1)/PTS(2)
         FLL   = 1.D-1
         FTL   = FTPF(MDLTPF,EPS)
         BQL   = PNSS(3)/PNSS(1)
         EQL   =       SLNQL/SLTQL
         ENQL  = 2.D0*(SLNQL/SLBL )
         BETAEL= PNSS(1)*PTS(1)*RKEV*1.D20/(BB**2/(2.D0*RMU0))
         TEL   = PTS(1)
         TAUZL = PTS(1)/PTS(3)
         QL    = QP(NR)
         SL    = S(NR)
         RNL   = PNSS(1)*1.D1

         RGKL  = AR1RHOG(NR)/(RA*AR2RHOG(NR))
         WEXBL = WEXB(NR)
         IF(MDLKAI.EQ.64) THEN
            SHAT  = SQRT(2.D0*SL-1.D0+RKAP**2*(SL-1.D0)**2)
            FLS   = (0.7D0+2.4D0/(7.14D0*QL*SHAT+0.1D0))*FLL
            FLL   = 2.D0*FLS/(1.D0+1.D0/TAUL)
         ENDIF

         CALL TR_WEILAND_BRIDGE (ENL,EIL,EEL,TAUL,FLL,FTL,BQL,EQL,ENQL,ZL, &
              BETAEL,AZL,COLL,ELL,TEL,TAUZL,RA,QL,SL,EPS,EPSA,RNL,RLIST, &
              RNEQL,RKAP,RIWL,RISBL,BB,SEARCH,PMA,RGKL,WEXBL,ROTL,NR,IST, &
              CHIL,CHEL,DL,CHQL,DQL,SCHI,SCHE,SD,SCHQ,SDQ)

         CALL WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL) !,SCHI,SCHE,SD,SCHQ,SDQ)

      RETURN
      END SUBROUTINE WEILAND_DRIVER

!     ****************************************************************

      SUBROUTINE WEILAND_COEF(NR,CHIL,CHEL,DL,CHQL,DQL) !,SCHI,SCHE,SD,SCHQ,SDQ)

      USE TRCOMM, ONLY : &
           ADDW, ADDWD, ADDWP, AKDW, AKDWD, AKDWP, MDLWLD, NSM, PA, rkind
      IMPLICIT NONE
      INTEGER  :: NR
!      REAL(rkind)     :: SCHE, SCHI, SCHQ, SD, SDQ
      REAL(rkind),DIMENSION(5):: CHEL, CHIL, CHQL, DL, DQL
      INTEGER:: ns, ns1


!     The diagonal value of the transport coefficient matrix
!     is set to be zero if it becomes negative.

      IF(DL  (3).LT.0.D0) DL  (3)=0.D0
      IF(CHEL(2).LT.0.D0) CHEL(2)=0.D0
      IF(CHIL(1).LT.0.D0) CHIL(1)=0.D0
      IF(DQL (5).LT.0.D0) DQL (5)=0.D0
      IF(CHQL(4).LT.0.D0) CHQL(4)=0.D0

!     It is assumed that De=Di in the followings.

      IF(PA(3).EQ.3.D0) THEN
!         AKDW(NR,1)=SCHE
         AKDW(NR,1)=CHEL(2)
         DO NS=2,NSM
!            AKDW(NR,NS)=SCHI
            AKDW(NR,NS)=CHIL(1)
         ENDDO
         DO NS=1,NSM
!            ADDW(NR,NS)=SD
            ADDW(NR,NS)=DL(3)
         ENDDO
         IF(MDLWLD.NE.0) THEN
            ADDWD(NR,1,1)=DL(3)
            ADDWP(NR,1,1)=CHEL(3)
            DO NS=2,NSM
               ADDWD(NR,NS,1)=DL(3)
               ADDWP(NR,NS,1)=CHIL(3)
            ENDDO
            AKDWD(NR,1,1)=DL(2)
            AKDWP(NR,1,1)=CHEL(2)
            DO NS=2,NSM
               AKDWD(NR,NS,1)=DL(2)
               AKDWP(NR,NS,1)=CHIL(2)
            ENDDO
            DO NS1=2,NSM
               ADDWD(NR,1,NS1)=DL(3)
               ADDWP(NR,1,NS1)=CHEL(3)
               DO NS=2,NSM
                  ADDWD(NR,NS,NS1)=DL(3)
                  ADDWP(NR,NS,NS1)=CHIL(3)
               ENDDO
               AKDWD(NR,1,NS1)=DL(1)
               AKDWP(NR,1,NS1)=CHEL(1)
               DO NS=2,NSM
                  AKDWD(NR,NS,NS1)=DL(1)
                  AKDWP(NR,NS,NS1)=CHIL(1)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         AKDW(NR,1)=CHEL(2)
         AKDW(NR,2)=CHIL(1)
         AKDW(NR,3)=CHQL(4)
         AKDW(NR,4)=CHQL(4)
         ADDW(NR,1)=DL(3)
         ADDW(NR,2)=DL(3)
!         ADDW(NR,2)=0.D0
         ADDW(NR,3)=DQL(5)
         ADDW(NR,4)=DQL(5)
         IF(MDLWLD.NE.0) THEN
            ADDWD(NR,1,1)=DL(3)
            ADDWP(NR,1,1)=CHEL(3)
            ADDWD(NR,2,1)=DL(3)
            ADDWP(NR,2,1)=CHIL(3)
            DO NS=3,NSM
               ADDWD(NR,NS,1)=DQL(3)
               ADDWP(NR,NS,1)=CHQL(3)
            ENDDO
            AKDWD(NR,1,1)=DL(2)
            AKDWP(NR,1,1)=CHEL(2)
            AKDWD(NR,2,1)=DL(2)
            AKDWP(NR,2,1)=CHIL(2)
            DO NS=3,NSM
               AKDWD(NR,NS,1)=DQL(2)
               AKDWP(NR,NS,1)=CHQL(2)
            ENDDO
            ADDWD(NR,1,2)=DL(3)
            ADDWP(NR,1,2)=CHEL(3)
            ADDWD(NR,2,2)=DL(3)
            ADDWP(NR,2,2)=CHIL(3)
            DO NS=3,NSM
               ADDWD(NR,NS,2)=DQL(3)
               ADDWP(NR,NS,2)=CHQL(3)
            ENDDO
            AKDWD(NR,1,2)=DL(1)
            AKDWP(NR,1,2)=CHEL(1)
            AKDWD(NR,2,2)=DL(1)
            AKDWP(NR,2,2)=CHIL(1)
            DO NS=3,NSM
               AKDWD(NR,NS,2)=DQL(1)
               AKDWP(NR,NS,2)=CHQL(1)
            ENDDO
            DO NS1=3,NSM
               ADDWD(NR,1,NS1)=DL(5)
               ADDWP(NR,1,NS1)=CHEL(5)
               ADDWD(NR,2,NS1)=DL(5)
               ADDWP(NR,2,NS1)=CHIL(5)
               DO NS=3,NSM
                  ADDWD(NR,NS,NS1)=DQL(5)
                  ADDWP(NR,NS,NS1)=CHQL(5)
               ENDDO
               AKDWD(NR,1,NS1)=DL(4)
               AKDWP(NR,1,NS1)=CHEL(4)
               AKDWD(NR,2,NS1)=DL(4)
               AKDWP(NR,2,NS1)=CHIL(4)
               DO NS=3,NSM
                  AKDWD(NR,NS,NS1)=DQL(4)
                  AKDWP(NR,NS,NS1)=CHQL(4)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE WEILAND_COEF

!     ****************************************************************

      SUBROUTINE TR_WEILAND_BRIDGE (EN,EI,EE,TAU,FL,FT,BQ,EQ,ENQ,Z,BETAE,AZ, &
           COL,EL,TE,TAUZ,PR,Q,S,EPSR,E,N,LIST,RNEQ,KAPPA,RIW,RISB,BTOR, &
           SEARCH,MA,GKL,WEXBL,ROTL,IKL,ISTL,CHI,CHE,D,CHQ,DQ,SCHIL,SCHEL, &
           SDL,SCHQL,SDQL)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),   INTENT(IN):: &
           TE,BTOR,RIW,RISB,EI,TAU,FL,EN,FT,EE,BQ,EQ,ENQ,Z,TAUZ,BETAE,MA,E, &
           Q,S,KAPPA,PR,N,RNEQ,EPSR,AZ,COL,EL,LIST,SEARCH,GKL,WEXBL,ROTL
      INTEGER,INTENT(IN):: IKL,ISTL
      REAL(rkind),DIMENSION(5),INTENT(INOUT):: CHI,CHE,D,CHQ,DQ
      REAL(rkind),INTENT(OUT)  :: SCHIL,SCHEL,SDL,SCHQL,SDQL
      INTEGER:: ICP, IR, NDISP
      REAL(rkind)   :: &
           A, B, C, DS, DTOT, ENN, G, H, LAMB, LNH, NI, NQ, R, TAUI, TVR, &
           VEI, WI, WIS, WR, WRS, WST, ZEFF, ZX
      REAL(rkind),DIMENSION(5,100):: U
      COMPLEX(rkind):: W
      COMPLEX(rkind),DIMENSION(10):: RP

      INTEGER:: &
           I,IW,IX,IK,IST,ITS,ITL,ITERA,ITC,ISB,LPRINTIN,NDIM,NEQ,IRET, &
           SEARCHMODE
      REAL(rkind)   :: &
           THRD,STR,XIH,RFL,H1,FTR,GQ,BF,ZE,CS,KPPA,RAV,ALP,ALF,KPC,D1,SI, &
           KIQ,KXQ,WDE,EPS,SCHI,SCHE,SD,SCHQ,SDQ,EA,GK,ETE,ETI,ETQ,AZL,ALAF, &
           GAV,VEF,BTA,EM,EIC,EEC,ENC,TAUC,FLC,FTC,EQC,ENQC,BETAEC,TAUZC,QC, &
           SC,ENHC,ZFS,KPS,CHIC,WEXB,ROT,WZIMAX,TOL,SHPE,SCHEF,DEF,LTH,LTE, &
           LN,LTQ,LNQ
      REAL(rkind),DIMENSION(5):: HPT
      REAL(rkind),DIMENSION(32):: CETAIN
      REAL(rkind),DIMENSION(10,10):: ZVR,ZVI
      COMPLEX(rkind):: HQ,WZ,WZP
      COMPLEX(rkind),DIMENSION(10):: ZZ
      COMMON/IRET/ IRET
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
      COMMON/EFFDIFF/ SCHI,SCHE,SD,SCHQ,SDQ
!      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/GRAD/ ENC,ENHC,EIC,EEC
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BF,GQ,SI,ZE,TAUZC,ENQC,EQC,KIQ,KXQ,ZFS
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUC,AZL,FTC,FLC,LPRINTIN,NDIM
      COMMON/BETAE/ BETAEC,QC,SC,CS,EM
      COMMON/TEST/ ALP,ALF,KPC,KPS
      COMMON/HQ/ H1,HQ,GAV,WZ,WZP,WZIMAX
      COMMON/SHAFRS/ ALAF
      COMMON/KAPPA/ KPPA,RAV
      COMMON/IK/ IK,IST,ITC,ITL,ITS,ITERA,TOL
      COMMON/TP/ EA,HPT
      COMMON/WROT/ WEXB,ROT
      COMMON/GRKVOT/ GK
      COMMON/ZV/ ZVR,ZVI
      COMMON/ZZ/ ZZ
      COMMON/NEQ/ NEQ
      COMMON/EM/ SHPE,CHIC,SCHEF,DEF
      COMMON/LT/ LTH,LTE,LN,LTQ,LNQ
      COMMON/ISB/ ISB
      COMMON/SEARCHMODE/ SEARCHMODE
      COMMON/WDE/ WDE ! added


!     ICP : controlling printout in TR_WEILAND_BRIDGE (0: off)
      ICP=0
      IF(ICP.NE.0) THEN
         WRITE(6,*) '//////////////////////////////////////////////////'
      ENDIF

      GK   = GKL     ! GRHO1/(a*GRHO2)
      WEXB = WEXBL   ! EXB shearing rate according to Hahm and Burrell
      ROT  = ROTL    ! factor that multiplies the EXB shearing rate
      IK   = IKL     ! the space profile index for WZL(IK)
      IST  = ISTL    ! 1:analytical form used, other:previous value used

      ZFS = 0.D0     ! the product of charge and fraction (to Ne) of fast particles

      THRD=1.D0/3.D0 ! fixed coefficient (THiRD)
      TVR=2.D0*THRD  ! fixed coefficient (Two thiRd)
      FTR=5.D0/3.D0  ! fixed coefficient (Five ThiRd)
      STR=7.D0/3.D0  ! fixed coefficient (Seven ThiRd)
      BTA=1.5D0      ! fixed, parameter used in the collision model

      EPS=EPSR       ! local inverse aspect ratio
      SEARCHMODE=INT(SEARCH) ! the way to search eigenvalues
      AZL=AZ         ! Aq (atomic number of impurity)
      EM=EL          ! factor that multiplies electromagnetic effects
      ZE=Z           ! Zq
      BF=BQ          ! nq/ne

      G=1.D0-Z*BQ    ! transforming factor from ne to ni
      ENN=1.D0-Z*BQ*EN/ENQ ! (Lne/Lni)*(ni/ne)
      IF(ABS(ENN).GE.0.001D0) GO TO 89
      IF(ENN.LT.0.D0) ENN=-0.001D0
      IF(ENN.GT.0.D0) ENN= 0.001D0
   89 CONTINUE
      ENC=EN         ! 2*Lne/R
      ENHC=G*EN/ENN  ! 2*Lni/R
      ENQC=ENQ       ! 2*Lnq/R
!      ENQC=EN        ! 2*Lne/R (original definition)

      EIC=EI         ! Lni/Lti
      EEC=EE         ! Lne/Lte
      EQC=EQ         ! Lnq/Ltq
      TAUC=TAU       ! Te/Ti
      TAUZC=TAUZ     ! Te/Tq

      FLC=DSQRT(FL)  ! FLR parameter (sqrt((Kperp*rhos)**2=0.1))
      FTC=FT         ! trapped particle fraction
      BETAEC=BETAE   ! electron beta
      QC=Q           ! safety factor
      SC=S           ! magnetic shear
      NEQ=INT(RNEQ)  ! number of equations (NDISP)
      ISB=INT(RISB)  ! ballooning parameter
      NDISP=NEQ      ! number of equations (NDISP)
      NDIM=5         ! dimension of transport matrix
      KPPA=KAPPA     ! elongation
      EA=E           ! inverse aspect ratio
      IW=INT(RIW)    ! controls printout in DIFFTD
      R=PR/E         ! major radius
      D1=6.462D0*DSQRT(MA)/(R*BTOR**2) ! machine dependent parameter

      IX=IK          ! the space profile index for U(J,IX)

      DO I=1,31
         CETAIN(I)=0.D0 ! control vector
      ENDDO
!      CETAIN(32)=1.D-15
      CETAIN(32)=0.001  ! accuracy parameter in the NAG routine
      LPRINTIN=INT(LIST)

      TAUI=1.D0/TAU  ! Ti/Te
      RFL=SQRT(FL)   ! FLR parameter (same with FLC)
      GQ=1.D0-Z*BQ   ! transforming factor from ne to ni (same with G)
      NQ=BQ*N        ! nq
      NI=GQ*N        ! ni
      ZEFF=(NI+Z*Z*NQ)/N ! Zeff

      LN =0.5D0*R*ENC/PR  ! Lne/a
      LNH=0.5D0*R*ENHC/PR ! Lni/a
      LNQ=0.5D0*R*ENQC/PR ! Lnq/a
      LTE=LN/EEC          ! Lte/a
      LTH=LNH/EIC         ! Lti/a
      LTQ=LNQ/EQC         ! Ltq/a

      ETI=ENHC/EIC   ! 2*Lti/R
      ETE=ENC/EEC    ! 2*Lte/R
      ETQ=ENQC/EQC   ! 2*Ltq/R

      KIQ=ENQC/ENHC  ! Lnq/Lni
      KXQ=ENC/ENQC   ! Lne/Lnq

      CS=3.095D5*DSQRT(TE/MA) ! sound speed
      IF(ICP.NE.0) THEN
         WRITE(*,00126) EN,EI,EE,FL,TAU
00126    FORMAT(2X,'EN=',F8.3,' EI=',F8.3,' EE=',F8.3,' FL=',F8.3, ' TAU=',G12.4)
         WRITE(*,00299) ENQ,ENHC
00299    FORMAT(2X,'ENQ=',G11.3,' ENH=',G11.3)
         WRITE(*,00129) BETAE
00129    FORMAT(2X,'BETAE=',G11.3)
         WRITE(*,00131) Q,S,CS
00131    FORMAT(2X,'q=',G12.4,' S=',G12.4,' CS=',G12.4)
         WRITE(*,00132) FT
00132    FORMAT(2X,'FT=',G12.4)
         WRITE(*,00134) ALAF
00134    FORMAT(2X,' Ballooning alpha =',G11.3)
      ENDIF

      WST=DSQRT(FL)*CS/(PR*ABS(LN))  ! diamagnetic frequency (omega_star)
      WDE=ABS(EN)*WST           ! curvature drift frequency (omega_drift_electron)
      LAMB=15.95D0-DLOG(DSQRT(N)/TE) ! coulomb logarithm (lambda)
      VEI=9.19D2*NI*LAMB/TE**1.5D0   ! collisionality (nu_electron_ion)
      VEF=COL*VEI/(EPS*WDE)     ! effective electron ion collision frequency for trapped electrons, normalized by omega_de
      WEXB=WEXB/WDE ! ExB shearing rate should be normalized with WDE

      IF(ICP.NE.0) THEN
         WRITE(*,00130) WST,VEF,VEI,ZEFF,COL
00130    FORMAT(2X,'WST=',G11.3,' VEF=',G11.3,' VEI=',G11.3, ' ZEFF=',G11.3,' COL=',G11.3)

         H=0.5D0*ABS(S)/q       ! not used
         WRITE(*,00133) H
00133    FORMAT(2X,'H=',G12.4)
      ENDIF
!   -----------------------------------------------

      IF(ICP.NE.0) THEN
         A=1.D0-EN*(1.D0+10.D0/(3.D0*TAU))-FL*(1.D0+EI+5.D0*EN/3.D0)/TAU
         B=EI-7.D0/3.D0+5.D0*EN*(1.D0+1.D0/TAU)/3.D0 +5.D0*FL*(1.D0+EI)/(3.D0*TAU)
         B=B*EN/TAU
         C=A/(2.D0*(1.D0+FL))
         DS=C*C-B/(1.D0+FL)
         IF(DS.LT.0.D0) GOTO 140
         WR=C+SQRT(DS)          ! not used
         WI=0.D0                ! not used
         GO TO 160
 140     WR=C                   ! not used
         WI=SQRT(-DS)           ! not used
         WRITE(*,170) WR,WI
 160     CONTINUE
 170     FORMAT(2X,'WR=',F7.3,' WI=',F7.3)
      ENDIF

      ITC=1      ! 1:iteration, other:previous eigenvalues used
      ITL=80     ! Maximum number of iterations
      TOL=0.01D0 ! relative error for convergence
!      IST=1

!      IF(IX.GE.48) write(6,'(I3,2F15.7)') IX,HQ
      CALL disp9t(NDISP,ZZ)

      IF(ICP.NE.0) THEN
         WRITE(*,174) ISB
 174     FORMAT(' ISB=',I5)
         WRITE(*,175) ITC,ITS,ITERA
 175     FORMAT('  ITC=',I5,' ITS=',I5,' ITER=',I5)
         WRITE(*,177) WZ,WZP
 177     FORMAT('  WZ=',2G11.3,' WZP=',2G11.3)
      ENDIF

      IR=0  ! the number of unstable roots (given by the following)

      DO 00199 I=1,NDISP
      ZX=DIMAG(ZZ(I))
      IF(ZX.LE. 0.001) GOTO 00199
      IR=IR+1
      RP(IR)=ZZ(I) ! unstable roots found by disp9t
00199 CONTINUE

      IF(ICP.NE.0) THEN
         WRITE(*,310) IR
 310     FORMAT(2X,' IR=',I4)
         WRITE(*,00134) ALAF
      ENDIF

      DO 0200 I=1,IR
      W=RP(I)
      WR=EN*DREAL(W) ! real part of unstable roots
      WI=EN*DIMAG(W) ! imaginary part of unstable roots
      IF(ICP.NE.0) WRITE(*,311) WR,WI,I
      WRS=WDE*DREAL(W)
      WIS=WDE*DIMAG(W)
      IF(ICP.NE.0) WRITE(*,321) WRS,WIS
 0200 CONTINUE
  311 FORMAT(//,2X,'WR=',G11.3,' WI=',G11.3,' I=',I5)
  321 FORMAT(' WRS=',G11.3,' WIS=',G11.3)

      WZ=EN*WZ
      HQ=EN*HQ
      IF(ICP.NE.0) THEN
         WRITE(*,00128) ALF,ALP,WZ,KAPPA
00128    FORMAT(/,2X,' ALF=',G11.3,' ALP=',G11.3,' WZ=',2G11.3,/,' KAPPA=',G11.3)
         WRITE(*,312) H1,HQ,GAV,RAV
 312     FORMAT(//,2X,'H1=',G12.4,' HQ=',2G12.4,' GAV=',G12.4,/,' RAV=',G12.4)
      ENDIF
      U(1,IX)=TAUI*TE
      U(2,IX)=TE
      U(3,IX)=N
      U(4,IX)=TE/TAUZ
      U(5,IX)=BQ*N

      CALL DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DQ)

      IF(ICP.NE.0) THEN
         WRITE(*,330) SCHI,SCHE,SD,SCHQ,SDQ
 330     FORMAT(/,2X,'CHIEFF=',G11.3,' CHEEFF=',G11.3,' DEFF=',G11.3, ' CHQEFF=',G11.3,' DQEFF=',G11.3)

!         WRITE(*,331) CHE(2),SCHEF
! 331     FORMAT(' CHE(2)=',G11.3,' SCHEF=',G11.3)
         WRITE(*,331) CHQ(5),SCHQ
 331     FORMAT(' CHQ(5)=',G11.3,' SCHQ=',G11.3)
         DTOT=SD+DEF
         WRITE(*,332) D(3),DEF,DTOT
 332     FORMAT(' D(3)=',G11.3,' DEF=',G11.3,' DTOT=',G11.3)
         WRITE(*,335) XIH,SHPE,CHIC,SCHEF
 335     FORMAT('  XIH=',G11.3,' SHPE=',G11.3,' CHIC=',G11.3,' SCHEF=',G11.3)
      ENDIF
!00150 CONTINUE
      SCHIL = SCHI
      SCHEL = SCHE
      SDL   = SD
      SCHQL = SCHQ
      SDQL  = SDQ

      RETURN
      END SUBROUTINE TR_WEILAND_BRIDGE


!     ***********************************************************

!            IFS/PPPL Model

!     ***********************************************************

      SUBROUTINE IFSPPPL_DRIVER(NSTM,NRMAX,RN,RR,DR,RJCB,RHOG,RHOM,QP, &
           S,EPSRHO,RKPRHOG,RT,BB,AMP,AME,PNSS,PTS,RNFL,RBEEDG,NSMAX, &
           AR1RHOG,AR2RHOG,AKDW)

      USE trcomm,ONLY: rkind
      USE libitp
      IMPLICIT NONE

      INTEGER,INTENT(IN):: NSTM,NRMAX,NSMAX
      REAL(rkind)   ,INTENT(IN):: RR,DR,BB,AMP,AME,RBEEDG
      REAL(rkind),DIMENSION(NRMAX,NSMAX),INTENT(IN):: RN, RT
      REAL(rkind),DIMENSION(NRMAX),INTENT(IN):: &
           RJCB,RHOG,RHOM,QP,S,EPSRHO,RKPRHOG,RNFL,AR1RHOG,AR2RHOG
      REAL(rkind),DIMENSION(NSMAX)    ,INTENT(IN):: PNSS,PTS
      REAL(rkind),DIMENSION(NRMAX,NSTM),INTENT(OUT)::AKDW
      integer :: ii, ierr, nr
      integer,parameter :: ipin=7,iptmp=8,ipout=9,screen=6
      integer,dimension(32):: switches
      real :: &
           btesla, chie, chii, gnu, grhoi, gtau, gvti, ne19, omegaexb, &
           rmajor, tekev, tikev, zchie1, zchie2, zchii1, zchii2, zchiicyc, &
           zeps, zkappa, znbne, zncne, znine, zq, zrln, zrlt, zrlt1, zrlt2, &
           zshat


      ierr=0

      do ii = 1, 32
         switches(ii) = 0
      end do
      switches(1)  = 0 ! 0: it produces no diagnostic output
      switches(2)  = 0 ! 0: it uses inputs ne19, tekev, tilev and btesla
                       ! 1: it uses inputs grhoi, gvti, gnu and gtau
      switches(3)  = 0 ! 0: the 1995 model, 1: the 1994 model
      switches(4)  = 1 ! 0: use gnu as given
                       ! 1: the definition in Dorland's IFS/PPPL routine
      switches(5)  = 1 ! 0: won't relax the restrictions on znu
                       ! 1: allows znu to be larger than 10.0
      switches(30) = 0
      switches(31) = 0
      switches(32) = 0

      DO NR=1,NRMAX-1
         znine  = SNGL( (RN(NR+1,2)+RN(NR,2)) /(RN(NR+1,1)+RN(NR,1)))
         zncne  = 0.0
         znbne  = SNGL( (RNFL(NR+1 )+RNFL(NR ))/(RN (NR+1,1)+RN (NR,1)))
         zrlt   =-SNGL(RR/(0.5D0*(RT(NR+1,2)+RT(NR,2))) &
                           *(RT(NR+1,2)-RT(NR,2))/DR*RJCB(NR))
         zrln   =-SNGL(RR/(0.5D0*(RN(NR+1,1)+RN(NR,1))) &
                           *(RN(NR+1,1)-RN(NR,1))/DR*RJCB(NR))
         zq     = SNGL(QP(NR))
         zshat  = SNGL(S(NR))
         zeps   = SNGL(EPSRHO(NR))
         zkappa = SNGL(RKPRHOG(NR))
         gnu    = SNGL((AME/AMP)*1.5625D-15*RN(NR,2)*1D20/RT(NR,1)**1.5D0)
!         gnu    = 2.1*rmajor*ne19/(tekev**1.5 * tikev**0.5)
         gtau   = SNGL( (RT(NR+1,2)+RT(NR,2))/(RT(NR+1,1)+RT(NR,1)))

         ne19   = SNGL(0.5D0*(RN(NR+1,1)+RN(NR,1))*1.D1)
         tekev  = SNGL(0.5D0*(RT(NR+1,1)+RT(NR,1)))
         tikev  = SNGL(0.5D0*(RT(NR+1,2)+RT(NR,2)))
         rmajor = SNGL(RR)
         btesla = SNGL(BB)
         grhoi  = SNGL(6.46D-3*SQRT(0.5D0*(RT(NR+1,2)+RT(NR,2)))/BB)
         gvti   = SNGL(2.19D5*SQRT(0.5D0*(RT(NR+1,2)+RT(NR,2))))

         CALL IFSPPPL( znine, zncne, znbne, zrlt, zrln, zq, zshat, zeps, &
              zkappa, omegaexb, ne19, tekev, tikev, rmajor, btesla, &
              switches, grhoi, gvti, gnu, gtau, chii, chie, zchiicyc, &
              zchii1, zchii2, zchie1, zchie2, zrlt1, zrlt2, ierr )
!         IF(IERR.NE.0) THEN
!            WRITE(6,*) 'XX IFS/PPPL : ERROR IERR=',IERR
!            STOP
!         ENDIF

         AKDW(NR,1) = DBLE(chie)*AR1RHOG(NR)/AR2RHOG(NR)
         AKDW(NR,2) = DBLE(chii)*AR1RHOG(NR)/AR2RHOG(NR)
      ENDDO

      NR=NRMAX
         znine  = SNGL(PNSS(2)/PNSS(1))
         zncne  = 0.0
         znbne  = SNGL(RBEEDG)
         zrlt   =-SNGL(RR/PTS(2)*DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2), &
                                         RHOG(NR),RHOM(NR),RHOM(NR-1)))
         zrln   =-SNGL(RR/PTS(1)*DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1), &
                                         RHOG(NR),RHOM(NR),RHOM(NR-1)))
         zq     = SNGL(QP(NR))
         zshat  = SNGL(S(NR))
         zeps   = SNGL(EPSRHO(NR))
         zkappa = SNGL(RKPRHOG(NR))
         gnu    = SNGL((AME/AMP)*1.5625D-15*RN(NR,2)*1D20/RT(NR,1)**1.5D0)
         gtau   = SNGL(PTS(2)/PTS(1))

         ne19   = SNGL(PNSS(1)*1.D1)
         tekev  = SNGL(PTS(1))
         tikev  = SNGL(PTS(2))
         rmajor = SNGL(RR)
         btesla = SNGL(BB)
         grhoi  = SNGL(6.46D-3*SQRT(PTS(2))/BB)
         gvti   = SNGL(2.19D5*SQRT(PTS(2)))

         CALL IFSPPPL( znine, zncne, znbne, zrlt, zrln,zq, zshat, zeps, &
              zkappa, omegaexb, ne19, tekev, tikev, rmajor, btesla, &
              switches, grhoi, gvti, gnu, gtau,chii, chie, zchiicyc, &
              zchii1, zchii2, zchie1, zchie2,zrlt1, zrlt2, ierr )
!         IF(IERR.NE.0) THEN
!            WRITE(6,*) 'XX IFS/PPPL : ERROR IERR=',IERR
!            STOP
!         ENDIF

         AKDW(NR,1) = DBLE(chie)*AR1RHOG(NR)/AR2RHOG(NR)
         AKDW(NR,2) = DBLE(chii)*AR1RHOG(NR)/AR2RHOG(NR)

      RETURN
      END SUBROUTINE IFSPPPL_DRIVER
