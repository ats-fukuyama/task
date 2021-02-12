!     $Id: fpnfrr.f90,v 1.20 2013/01/14 16:48:26 fukuyama Exp $

! **********************************
!      Nuclear Fusion Reaction Rate
! **********************************

      MODULE fpnfrr

      USE fpcomm
      USE fpinit

      PUBLIC NF_REACTION_COEF
      PUBLIC NF_REACTION_RATE
      PUBLIC NF_REACTION_COEF_BT
      PUBLIC NF_REACTION_RATE_BT 
      PUBLIC NF_REACTION_RATE_TT
      PUBLIC SIGMA_E_Bosch
      PUBLIC CHECK_NF_RATE_ON_STDIO

      PRIVATE

      REAL(8),dimension(1:5,1:6):: Duane_Coef
      DATA Duane_Coef/ &
             46.097D0,  372.D0,4.360D-4,1.220D0,  0.D0, & ! D-D reaction 1 ID=1
             47.880D0,  482.D0,3.080D-4,1.177D0,  0.D0, & ! D-D reaction 2 ID=2
             45.950D0,50200.D0,1.368D-2,1.076D0,409.D0, & ! D-T reaction   ID=3
             89.270D0,25900.D0,3.980D-3,1.297D0,647.D0, & ! D-He3 reaction ID=4
             38.390D0,  448.D0,1.020D-3,2.090D0,  0.D0, & ! T-T reaction   ID=5
            123.100D0,11250.D0,    0.D0,   0.D0,  0.D0/   ! T-He3 reaction ID=6

      REAL(8),dimension(1:5,1:6):: Bosch_Coef_A
      DATA Bosch_Coef_A/ &
           5.5576D4, 2.1054D2, -3.2638D-2, 1.4987D-6, 1.8181D-10, & ! D-D reaction 1 ID=1
           5.3701D4, 3.3027D2, -1.2706D-1, 2.9327D-5, -2.5151D-9, & ! D-D reaction 2 ID=2
            6.927D4,  7.454D8,    2.050D6,  5.2002D4,       0.D0, & ! D-T reaction   ID=3
           5.7501D6, 2.5226D3,   4.5566D1,      0.D0,       0.D0, & ! D-He3 reaction ID=4
               0.D0,     0.D0,       0.D0,      0.D0,       0.D0, & ! T-T reaction   ID=5
               0.D0,     0.D0,       0.D0,      0.D0,       0.D0/   ! T-He3 reaction ID=6

      REAL(8),dimension(1:4,1:6):: Bosch_Coef_B
      DATA Bosch_Coef_B/ &
                 0.D0,       0.D0,      0.D0,      0.D0, & ! D-D reaction 1 ID=1
                 0.D0,       0.D0,      0.D0,      0.D0, & ! D-D reaction 2 ID=2
               6.38D1,   -9.95D-1,  6.981D-5,  1.728D-4, & ! D-T reaction   ID=3
           -3.1995D-3, -8.5530D-6, 5.9014D-8,      0.D0, & ! D-He3 reaction ID=4
                 0.D0,       0.D0,      0.D0,      0.D0, & ! T-T reaction   ID=5
                 0.D0,       0.D0,      0.D0,      0.D0/   ! T-He3 reaction ID=6

      REAL(8),dimension(1:7,1:6):: Bosch_Coef_C
      DATA Bosch_Coef_C/ &
          5.65718D-12, 3.41267D-3, 1.99167D-3,       0.D0, 1.05060D-5,       0.D0, 0.D0, & ! D-D reaction 1 ID=1
          5.43360D-12, 5.85778D-3, 7.68222D-3,       0.D0,-2.96400D-6,       0.D0, 0.D0, & ! D-D reaction 2 ID=2
           1.17302D-9, 1.51361D-2, 7.51886D-2, 4.60643D-3, 1.35000D-2,-1.06750D-4,1.366D-5, & ! D-T reaction   ID=3
          5.51036D-10, 6.41918D-3,-2.02896D-3, -1.9108D-5, 1.35776D-4,       0.D0, 0.D0, & ! D-He3 reaction ID=4
                 0.D0,       0.D0,      0.D0,        0.D0,       0.D0,       0.D0, 0.D0, & ! T-T reaction   ID=5
                 0.D0,       0.D0,      0.D0,        0.D0,       0.D0,       0.D0, 0.D0/   ! T-He3 reaction ID=6

      REAL(8),dimension(1:6):: Rest_mass_E ! ID=1-6
      DATA Rest_mass_E/ & ! in keV
           937814.D0, 937814.D0, 1124656.D0, 1124572.D0, 0.D0, 0.D0/

      REAL(8),dimension(1:6):: Gamov_const ! ID=1-6
      DATA Gamov_const/ &
           31.3970D0, 31.3970D0, 34.3827D0, 68.7508D0, 0.D0, 0.D0/ !correct
!           31.3970D0, 31.790D0, 34.3827D0, 68.7508D0, 0.D0, 0.D0/  ! wrong

      INTEGER:: NSA_H,NSA_D,NSA_T,NSA_HE3,NSA_HE4,NSB_D,NSB_T,NSB_HE3
      INTEGER,PARAMETER:: NEMAX=101 ! Number of energy mesh

      REAL(8),DIMENSION(4,4,NEMAX,NEMAX,6):: USV2D
!      INTEGER:: NSB1_NF(6),NSB2_NF(6)

      contains
!-------------------------------------------------
      SUBROUTINE NF_bt_CROSS_SECTION(ID,energy,temp_i,sigmav_bt)
! The first term in Mikkelsen ( NF vol29 eq. (7) )
! Cross section is chosen from Bosch (NF vol.32 eq. 9)
      IMPLICIT NONE

      INTEGER,intent(in):: ID
      REAL(8),intent(in):: energy, temp_i
      REAL(8),intent(out):: sigmav_bt ! [m^3/t]
      REAL(8):: E_star, v_star, capR, capD, E_star_in_keV
      REAL(8):: beta, v_th, v_b, gamma, mu, capB
      REAL(8):: E_b, m_b, m_i, T_i
! E_b and T_i should be in keV or kJ

      E_b=energy ! [keV]
      m_b=2.D0*AMP
      T_i=temp_i ! [keV]
      m_i=2.D0*AMP
      
      IF(energy.gt.0.D0)THEN
         ! define variables
         mu = m_i*m_b/(m_i+m_b)
         v_th = sqrt(2.D0*T_i/m_i*AEE*1.D3)
         v_b = sqrt(2.D0*E_b/m_b*AEE*1.D3)
         capB= SQRT(2.D0/mu)*Gamov_const(2)*SQRT(AEE*1.D3)
         beta = (v_th/v_b)**2*capB/(2.D0*v_b)
         capR = 0.5D0*beta + 1.D0/27.D0
         capD = sqrt( beta*(0.25D0*beta+1.D0/27.D0) )
         v_star = v_b*sqrt(1.D0/3.D0 + (capR+capD)**(1.D0/3.D0) + (capR-capD)**(1.D0/3.D0))
         E_star = 0.5D0*mu*v_star**2
         E_star_in_keV = E_star*1.D-3/AEE
         gamma = 1.D0+2.D0*(v_star-v_b)/v_star
         
         sigmav_bt = SIGMA_E_Bosch(E_star_in_keV,ID)*v_star &
              *EXP(-((v_star-v_b)/v_th)**2)*(v_star/v_b)/(sqrt(gamma))
      ELSE
         sigmav_bt =0.D0
      END IF

      IF(sigmav_bt.lt.0.D0) WRITE(*,'(A,3E14.6)') "negative sigmav_bt", energy, temp_i, sigmav_bt
      IF(sigmav_bt.ne.sigmav_bt) WRITE(*,'(A,99E14.6)') "NaN sigmav_bt", energy, temp_i, sigmav_bt,&
           SIGMA_E_Bosch(E_star_in_keV,ID), v_star, v_th, v_b , gamma

      END SUBROUTINE NF_bt_CROSS_SECTION
!-------------------------------------------------
!---- thermal-thermal reactivity (According to Bosch NF1992)
      FUNCTION SIGMAV_TT(in_T, ID)

      IMPLICIT NONE
      REAL(8):: SIGMAV_TT       ! [m^3/s], [cm^3/s]=10^-6[m^3/s]
      REAL(8),INTENT(IN):: in_T ! ion Temperature in keV
      INTEGER,INTENT(IN):: ID
      REAL(8):: theta, xi

      theta = in_T/(&
           1.D0- in_T*(Bosch_Coef_C(2,ID)+in_T*(Bosch_Coef_C(4,ID)+in_T*Bosch_Coef_C(6,ID)))&
          /(1.D0+in_T*(Bosch_Coef_C(3,ID)+in_T*(Bosch_Coef_C(5,ID)+in_T*Bosch_Coef_C(7,ID) ) ) ) &
           )

      xi = (Gamov_const(ID)**2/(4.D0*theta))**(1.D0/3.D0)
      SIGMAV_TT=Bosch_Coef_C(1,ID)*theta*SQRT(xi/(Rest_mass_E(ID)*in_T**3))*EXP(-3.D0*xi)*1.D-6

      END FUNCTION SIGMAV_TT
!-------------------------------------------------
!---- Astrophysical S Function (According to Bosch NF1992) ----
      FUNCTION S_FUNCTION(E,ID)
      IMPLICIT NONE
      REAL(8):: S_FUNCTION
      REAL(8),INTENT(IN):: E ! Energy in CM frame[keV]
      INTEGER,INTENT(IN):: ID
      REAL(8):: denominator, numerator
      INTEGER:: i

      IF(E.LE.0.D0) THEN
         S_FUNCTION=0.D0
      ELSE IF(Gamov_const(ID)/SQRT(E).GT.200.D0) THEN
         S_FUNCTION=0.D0
      ELSE
         numerator=0.D0
         denominator=0.D0
         DO i=5,1,-1
            numerator = numerator + Bosch_Coef_A(i,ID)*E**(i-1)
         END DO
         DO i=4,1,-1
            denominator = denominator + Bosch_Coef_B(i,ID)*E**(i)
         END DO
         denominator = denominator + 1.D0
         S_FUNCTION=numerator/denominator
      ENDIF
      RETURN
      END FUNCTION S_FUNCTION

!---- Crosssection of Nuclear Fusion Reaction (According to Bosch) ----

      FUNCTION SIGMA_E_Bosch(E,ID)
      IMPLICIT NONE
      REAL(8):: SIGMA_E_Bosch      ! [m^2] (1 barn = 10^{-28} m^2)
      REAL(8),INTENT(IN):: E ! Energy in [keV]
      INTEGER,INTENT(IN):: ID

      IF(E.LE.0.D0) THEN
         SIGMA_E_Bosch=0.D0
      ELSE IF(Duane_Coef(1,ID)/SQRT(E).GT.200.D0) THEN
         SIGMA_E_Bosch=0.D0
      ELSE
         SIGMA_E_Bosch=S_FUNCTION(E,ID)/(E*EXP(Gamov_const(ID)/SQRT(E)))*1.D-31
      ENDIF
      RETURN
      END FUNCTION SIGMA_E_BOSCH
      
!---- Crosssection of Nuclear Fusion Reaction (According to NRL Formulary) ----

      FUNCTION SIGMA_E(E,ID)
      IMPLICIT NONE
      REAL(8):: SIGMA_E      ! [m^2] (1 barn = 10^{-28} m^2)
      REAL(8),INTENT(IN):: E ! Energy in [keV]
      INTEGER,INTENT(IN):: ID

!      write(6,'(1P3E12.4)') E, &
!                            (Duane_Coef(4,ID)-Duane_Coef(3,ID)*E)**2+1.D0, &
!                             Duane_Coef(1,ID)/SQRT(E)

      IF(E.LE.0.D0) THEN
         SIGMA_E=0.D0
      ELSE IF(Duane_Coef(1,ID)/SQRT(E).GT.200.D0) THEN
         SIGMA_E=0.D0
      ELSE
         SIGMA_E=(Duane_Coef(5,ID) &
                +Duane_Coef(2,ID) &
                 /((Duane_Coef(4,ID)-Duane_Coef(3,ID)*E)**2+1.D0))  &
                /(E*(exp(Duane_Coef(1,ID)/SQRT(E))-1.D0))*1.D-28
      ENDIF
      RETURN
      END FUNCTION SIGMA_E

!---- Sigma v as a function of phi: multiplied by SQRT(RMASS/2) ----

      FUNCTION SIGMAV_PHI(E0L,E1L,PHI,ID)
      IMPLICIT NONE
      REAL(8):: SIGMAV_PHI
      REAL(8),INTENT(IN):: E0L,E1L,PHI
      INTEGER,INTENT(IN):: ID
      REAL(8):: E

      E=E0L-E1L*COS(PHI)
!      write(*,*) "E0L",E0L
      IF(E.LE.0.D0) THEN
         SIGMAV_PHI=0.D0
      ELSE
         IF(MODEL_NF_CS.le.0)THEN
            SIGMAV_PHI=SIGMA_E(E,ID)*SQRT(E)
         ELSE
            SIGMAV_PHI=SIGMA_E_Bosch(E,ID)*SQRT(E)
         END IF
      ENDIF
      RETURN
      END FUNCTION SIGMAV_PHI

!---- gyro averaged Sigma v : multiplied by SQRT(RMASS/2) ----

!      FUNCTION SIGMAV_E(E0L,E1L,ID)
      FUNCTION SIGMAV_E(E0L,E1L,ID)
      IMPLICIT NONE
      REAL(8):: SIGMAV_E
      REAL(8),INTENT(IN):: E0L,E1L
      INTEGER,INTENT(IN):: ID
      INTEGER,PARAMETER:: NPHIMAX=50 ! Number of phi mesh
      INTEGER:: NPHI
      REAL(8):: PHI,RSUM,DPHI

      DPHI=PI/NPHIMAX
      RSUM=0.D0
       DO NPHI=1,NPHIMAX
         PHI=DPHI*(NPHI-0.5D0)
!         PHI=0.D0 ! temporal
         RSUM=RSUM+SIGMAV_PHI(E0L,E1L,PHI,ID)
      ENDDO
      SIGMAV_E=RSUM/DBLE(NPHIMAX)   ! (1/(2*pi))*2*(\pi/nphimax)
      RETURN
      END FUNCTION SIGMAV_E

!---- calculation of relative kinetic energy ----

      SUBROUTINE RELATIVE_ENERGY(NTH1,NP1,NTH2,NP2,NSB1,NSB2,E0L,E1L)

      USE fpcomm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NP1, NTH1, NP2, NTH2, NSB1, NSB2
      REAL(8),INTENT(OUT):: E0L, E1L ! Energy in keV
      REAL(8):: VPARA1, VPARA2, VPERP1, VPERP2, V_MS, VA_MS
      REAL(8):: RED_MASS
      INTEGER:: NS1, NS2

      NS1=NS_NSB(NSB1)
      NS2=NS_NSB(NSB2)

      VPARA1=PM(NP1,NS1)*COSM(NTH1)*VTFD0(NSB1)
      VPARA2=PM(NP2,NS2)*COSM(NTH2)*VTFD0(NSB2)
      VPERP1=PM(NP1,NS1)*SINM(NTH1)*VTFD0(NSB1)
      VPERP2=PM(NP2,NS2)*SINM(NTH2)*VTFD0(NSB2)
!        MEAN SQUARE VELOCITY
      V_MS=VPERP1**2+VPERP2**2+(VPARA1-VPARA2)**2
      VA_MS=VPERP1**2+VPARA1**2
!        REDUCED MASS
      IF(MODEL_NF_CS.eq.0)THEN
         RED_MASS=AMFD(NSB1)
      ELSEIF(MODEL_NF_CS.eq.1)THEN
         RED_MASS=AMFD(NSB1)
      ELSEIF(MODEL_NF_CS.eq.2)THEN
         RED_MASS=AMFD(NSB1)*AMFD(NSB2)/(AMFD(NSB1)+AMFD(NSB2))
      END IF
!      MASS_A=AMFD(NSB1) ! incident particle (NSB1) mass and relative v

      E0L=0.5D0*RED_MASS*V_MS/(AEE*1.D3)
      E1L=RED_MASS*VPERP1*VPERP2/(AEE*1.D3)
!      EAL=0.5D0*RED_MASS*VA_MS/(AEE*1.D3)
!      EAL=0.5D0*MASS_A*V_MS/(AEE*1.D3)
      RETURN
      END SUBROUTINE RELATIVE_ENERGY
!---------------------------------------------------------
      SUBROUTINE DEFINE_NF_IDs

      USE fpcomm
      IMPLICIT NONE
      INTEGER:: NSA, NS, NSB, ID

      NSA_H=0
      NSA_D=0
      NSA_T=0
      NSA_HE3=0
      NSA_HE4=0
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         IF(PA(NS).EQ.1.D0.AND.PZ(NS).EQ.1.D0) NSA_H=NSA
         IF(PA(NS).EQ.2.D0.AND.PZ(NS).EQ.1.D0) NSA_D=NSA
         IF(PA(NS).EQ.3.D0.AND.PZ(NS).EQ.1.D0) NSA_T=NSA
         IF(PA(NS).EQ.3.D0.AND.PZ(NS).EQ.2.D0) NSA_HE3=NSA
         IF(PA(NS).EQ.4.D0.AND.PZ(NS).EQ.2.D0) NSA_HE4=NSA
      ENDDO
      NSB_D=0
      NSB_T=0
      NSB_HE3=0
      DO NSB=1,NSBMAX
         NS=NS_NSB(NSB)
         IF(PA(NS).EQ.2.D0.AND.PZ(NS).EQ.1.D0) NSB_D=NSB
         IF(PA(NS).EQ.3.D0.AND.PZ(NS).EQ.1.D0) NSB_T=NSB
         IF(PA(NS).EQ.3.D0.AND.PZ(NS).EQ.2.D0) NSB_HE3=NSB
      ENDDO

!     ---- identify reacting and produced species ----

      DO ID=1,NF_IDMAX
         NSA1_NF(ID)=0
         NSA2_NF(ID)=0
         SELECT CASE(ID)
         CASE(1)
            IF(NSB_D.NE.0) THEN
               NSB1_NF(ID)=NSB_D
               NSB2_NF(ID)=NSB_D
               NSA1_NF(ID)=NSA_T
               ENG1_NF(ID)=1.01D6
               NSA2_NF(ID)=NSA_H
               ENG2_NF(ID)=3.02D6
            ENDIF
         CASE(2)
            IF(NSB_D.NE.0) THEN
               NSB1_NF(ID)=NSB_D
               NSB2_NF(ID)=NSB_D
               NSA1_NF(ID)=NSA_HE3
               ENG1_NF(ID)=0.82D6
            ENDIF
         CASE(3)
            IF(NSB_D.NE.0.AND.NSB_T.NE.0) THEN
               NSB1_NF(ID)=NSB_D
               NSB2_NF(ID)=NSB_T
               NSA1_NF(ID)=NSA_HE4
               ENG1_NF(ID)=3.5D6
            ENDIF
         CASE(4)
            IF(NSB_D.NE.0.AND.NSB_HE3.NE.0) THEN
               NSB1_NF(ID)=NSB_D
               NSB2_NF(ID)=NSB_HE3
               NSA1_NF(ID)=NSA_HE4
               ENG1_NF(ID)=3.5D6
               NSA2_NF(ID)=NSA_H
               ENG2_NF(ID)=14.7D6
            ENDIF
         CASE(5)
            IF(NSB_T.NE.0) THEN
               NSB1_NF(ID)=NSB_T
               NSB2_NF(ID)=NSB_T
               NSA1_NF(ID)=NSA_HE4
               ENG1_NF(ID)=1.25D6   ! This value is an isotropic estimation
            ENDIF
         CASE(6)
            IF(NSB_T.NE.0.AND.NSB_HE3.NE.0) THEN ! One reaction is implemented
               NSB1_NF(ID)=NSB_T
               NSB2_NF(ID)=NSB_HE3
               NSA1_NF(ID)=NSA_HE4
               ENG1_NF(ID)=4.8D6
               NSA2_NF(ID)=NSA_D
               ENG2_NF(ID)=9.5D6
            ENDIF
         END SELECT
      END DO

      END SUBROUTINE DEFINE_NF_IDs
!---------------------------------------------------------
!---- Calculation of NF reaction rate coefficient ----
      SUBROUTINE NF_REACTION_COEF_BT ! SIGMAV_NF_BT is defined. available for DD only

      USE fpcomm
      IMPLICIT NONE
      INTEGER:: ID, NS1, NS2, NP, NR, NSB1, NSB2
      double precision:: RED_MASS, E_beam, sigmav_bt, temp_i

      sigmav_nf_bt(:,:,:)=0.D0
      CALL DEFINE_NF_IDs

      DO ID=1,2
         NSB1=NSB1_NF(ID)
         NSB2=NSB2_NF(ID)
         IF(NSB1.ne.0)THEN
            NS1=NS_NSB(NSB1)
         ELSE
            NS1=0
         END IF
         IF(NSB2.ne.0)THEN
            NS2=NS_NSB(NSB2)
         ELSE
            NS2=0
         END IF
         RED_MASS=AMFD(NSB1)*AMFD(NSB2)/(AMFD(NSB1)+AMFD(NSB2))

         IF(NSB1.NE.0.AND.NSB2.NE.0) THEN
            DO NR=NRSTART, NREND
               temp_i=RT_TEMP(NR,NS1)
               DO NP=NPSTART, NPEND
!                  E_beam=0.5D0*VTFD0(NSB1)**2*PM(NP,NSB1)**2*RED_MASS/(AEE*1.D3) ! keV
                  E_beam=0.5D0*VTFD0(NSB1)**2*PM(NP,NSB1)**2*AMFD(NSB1)/(AEE*1.D3) ! keV
                  CALL NF_bt_CROSS_SECTION(ID,E_beam,temp_i,sigmav_bt)
                  SIGMAV_NF_BT(NP,NR,ID)=sigmav_bt
!                  WRITE(*,'(A,2I5, 3E14.6)') "TEST1", NR, NP, temp_i, E_beam, sigmav_bt
               END DO
            END DO
         END IF
      END DO

      END SUBROUTINE NF_REACTION_COEF_BT
!*********************************************************
      SUBROUTINE NF_REACTION_RATE_BT(NR,ID)

      USE fpcomm
      USE libmpi
      USE fpfunc
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,ID
      INTEGER:: NSB1, NSB2, NSA1, NSA2, NS1, NS2, NTH, NP, VLOC, NS, NSA, NSA_
      double precision:: double_count, fact, RSUM, RSUM_SUM, F_maxwell, f_beam, cs_tt

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)
      NS1=NS_NSB(NSB1)
      NS2=NS_NSB(NSB2)

      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         IF(PA(NS).EQ.2.D0.AND.PZ(NS).EQ.1.D0) NSA_=NSA
      ENDDO

      IF(NSB1*NSB2.ne.0)THEN
         double_count=1.D0
         FACT = RNFD0(NSB2)*1.D20*double_count*RN_BULK(NR,NS1) !double count is not required for bt
         IF(ID.eq.1.or.ID.eq.2.or.ID.eq.5) double_count=0.5D0 ! used only for tt
         
         RSUM=0.D0
         RSUM_SUM=0.D0
         DO NP=NPSTART, NPEND
            IF(MODEL_EX_READ_Tn.eq.0)THEN
               f_maxwell = FPMXWL(PM(NP,NS1),NR,NS1)
            ELSE
               f_maxwell = FPMXWL_EXP(PM(NP,NS1),NR,NS1)               
            END IF
            DO NTH=1, NTHMAX
               f_beam = FNSB(NTH,NP,NR,NSB1)-f_maxwell
               RSUM = RSUM + FACT*f_beam*VOLP(NTH,NP,NS1) &
                    *RLAMDAG(NTH,NR)*RFSADG(NR)*SIGMAV_NF_BT(NP,NR,ID)
            END DO
         END DO
         CALL mtx_set_communicator(comm_np)
         CALL mtx_allreduce1_real8(RSUM,3,RSUM_sum,vloc)
         CALL mtx_reset_communicator
!         RATE_NF(NR,ID) = RSUM_SUM
         RATE_NF_BT(NR,ID) = RSUM_SUM

!         cs_tt = SIGMAV_TT(RT_TEMP(NR,NS1),ID)*1.D20*double_count*RN_BULK(NR,NS1)*RN_BULK(NR,NS2)
!         IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NRANK.eq.0)THEN
!            WRITE(6,'(A,E14.6,3I4)') '|-NF_REACTION_RATE:', TIMEFP, nrank, comm_nr%rank, comm_nsa%rank
!            WRITE(6,'(A,I5,A,2I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
!                 ,ID,':  ',NSB1,NSB2,' -> ',NSA1,NSA2
!            WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,     RATE_NF, RATE_NF_TT"
!            WRITE(6,'("  ",4I5,1P3E14.6)') ID,NR,NSB1,NSB2, RSUM_SUM/FACT, RSUM_SUM, cs_tt
!         END IF
!         RATE_NF_BB(NR,ID) = 0.D0
!         RATE_NF_TT(NR,ID) = cs_tt

      END IF

      END SUBROUTINE NF_REACTION_RATE_BT
!*********************************************************
      SUBROUTINE NF_REACTION_RATE_TT(NR,ID) ! IF MODELS>=1

      USE fpcomm
      USE libmpi
      USE fpfunc
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,ID
      INTEGER:: NSB1, NSB2, NS1, NS2
      double precision:: cs_tt, double_count

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NS1=NS_NSB(NSB1)
      NS2=NS_NSB(NSB2)

      double_count=1.D0
      IF(ID.eq.1.or.ID.eq.2.or.ID.eq.5) double_count=0.5D0

      cs_tt = double_count*SIGMAV_TT(RT_TEMP(NR,NS1),ID)*1.D20*RN_BULK(NR,NS1)*RN_BULK(NR,NS2)
      RATE_NF_TT(NR,ID) = cs_tt

      END SUBROUTINE NF_REACTION_RATE_TT
!*********************************************************
!---- Calculation of NF reaction rate coefficient ----
      SUBROUTINE NF_REACTION_COEF ! SIGMAV_NF is defined

      USE fpcomm
      IMPLICIT NONE
!      INTEGER,PARAMETER:: NEMAX=101 ! Number of energy mesh
      REAL(8),DIMENSION(NEMAX):: E0A,E1A
      REAL(8),DIMENSION(NEMAX,NEMAX):: SIGMAVA,FX,FY,FXY
      REAL(8),DIMENSION(4,4,NEMAX,NEMAX):: USV
      INTEGER:: NSB1,NSB2,ID,NE0,NE1,NTH1,NTH2,NP1,NP2,IERR,L,K,NS1,NS2
      REAL(8):: V1MAX,V2MAX,E0MAX,E1MAX,DELE0,DELE1,E0L,E1L,F0,RED_MASS

!---- identify possible fusion reactions ----
!     ---- identify related particle species ----

      CALL DEFINE_NF_IDs

      DO ID=1,NF_IDMAX
         NSB1=NSB1_NF(ID)
         NSB2=NSB2_NF(ID)
         IF(NSB1.ne.0)THEN
            NS1=NS_NSB(NSB1)
         ELSE
            NS1=0
         END IF
         IF(NSB2.ne.0)THEN
            NS2=NS_NSB(NSB2)
         ELSE
            NS2=0
         END IF

!     ---- identify whether reaction can be described or not ----

         IF(NSB1.NE.0.AND.NSB2.NE.0) THEN

!     ---- calculate reaction rate table for spline ----
!          - energy mesh is equally spaced in SQRT(E), not E

            V1MAX=VTFD0(NSB1)*PMAX(NS1)
            V2MAX=VTFD0(NSB2)*PMAX(NS2)
            IF(MODEL_NF_CS.eq.0)THEN
               RED_MASS=AMFD(NSB1)
               E0MAX= 0.5D0*AMFD(NSB1)*(V1MAX+V2MAX)**2 /(AEE*1.D3)
               E1MAX=       AMFD(NSB1)* V1MAX*V2MAX     /(AEE*1.D3)
            ELSEIF(MODEL_NF_CS.eq.2)THEN
               RED_MASS=AMFD(NSB1)*AMFD(NSB2)/(AMFD(NSB1)+AMFD(NSB2))
               E0MAX=0.5D0*RED_MASS*(V1MAX+V2MAX)**2 /(AEE*1.D3)!*2.D0
               E1MAX=      RED_MASS* V1MAX*V2MAX     /(AEE*1.D3)!*2.D0
            END IF
            DELE0=(SQRT(E0MAX)-1.D-6)/(NEMAX-1)
            DELE1= SQRT(E1MAX)/(NEMAX-1)
            DO NE0=1,NEMAX
               E0A(NE0)=(DELE0*(NE0-1))**2+1.D-6
            ENDDO
            DO NE1=1,NEMAX
               E1A(NE1)=(DELE1*(NE1-1))**2
            ENDDO
            DO NE0=1,NEMAX
               DO NE1=1,NEMAX
                  SIGMAVA(NE0,NE1)=SIGMAV_E(E0A(NE0),E1A(NE1),ID) &
                       *SQRT(AEE*1.D3*2.D0/RED_MASS) ! SQER(E) is not in keV
               ENDDO
            ENDDO
             
!     ---- calculate spline coefficients ----
            CALL SPL2D(E0A,E1A,SIGMAVA,FX,FY,FXY,USV,NEMAX,NEMAX,NEMAX,0,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX NF_REACTION_COEF: SPL2D: IERR=',IERR

!     ---- calculate reaction rate on 4D momentum mesh ----

            IF(ABS(MODELS).EQ.2) THEN
               DO NP2=1,NPMAX
                  DO NTH2=1,NTHMAX
                     DO NP1=NPSTART,NPEND
                        DO NTH1=1,NTHMAX
                           CALL RELATIVE_ENERGY(NTH1,NP1,NTH2,NP2,NSB1,NSB2,E0L,E1L)
                           CALL SPL2DF(E0L,ABS(E1L), &
                                SIGMAV_NF(NTH1,NP1,NTH2,NP2,ID), &
                                E0A,E1A,USV,NEMAX,NEMAX,NEMAX,IERR)
                           IF(IERR.NE.0) WRITE(6,*) &
                                'XX NF_REACTION_COEF: SPL2DF: IERR=',IERR
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF ! MODELS=2

            IF(MODELS.EQ.3) THEN
               DO NP2=1,NPMAX
                  DO K=0,LLMAX_NF
                     DO NP1=NPSTART,NPEND
                        DO L=0,LLMAX_NF
                           SIGMAV_LG(L,NP1,K,NP2,ID)=0.D0
                        END DO
                     END DO
                  END DO
               END DO
               DO NP2=1,NPMAX
                  DO NTH2=1,NTHMAX
                     DO NP1=NPSTART,NPEND
                        DO NTH1=1,NTHMAX
                           CALL RELATIVE_ENERGY(NTH1,NP1,NTH2,NP2,NSB1,NSB2,E0L,E1L)
                           CALL SPL2DF(E0L,ABS(E1L),F0, &
                                E0A,E1A,USV,NEMAX,NEMAX,NEMAX,IERR)
                           IF(IERR.NE.0) WRITE(6,*) &
                                'XX NF_REACTION_COEF_LG: SPL2DF: IERR=',IERR
                           DO L=0,LLMAX_NF
                              DO K=0,LLMAX_NF
                                 SIGMAV_LG(L,NP1,K,NP2,ID)= &
                                      SIGMAV_LG(L,NP1,K,NP2,ID) &
                                      +F0*PL_NF(L,NTH1)*PL_NF(K,NTH2) &
                                      *SINM(NTH1)*SINM(NTH2)*DELTH**2 &
                                      *0.5D0**2*(2.D0*L+1.D0)*(2.D0*K+1.D0)
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF ! MODELS=3

         END IF ! NSB1.NE.0.AND.NSB2.NE.0
      ENDDO ! ID

      RETURN
      END SUBROUTINE NF_REACTION_COEF
!*********************************************************
!===========================================================
!     CALCULATION OF INTEGRAL SIGMA*V*f(p_a)*f(p_b)
!===========================================================
!     out RATE_NF, _BB, _TT, _D1, _D2. MODELS=-2 does not use _D1, _D2.
      SUBROUTINE NF_REACTION_RATE(NR,ID)

      USE fpcomm
      USE libmpi
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,ID
      INTEGER:: NSB1, NSB2, NSA1, NSA2, NP1, NP2, NTH1, NTH2, VLOC, NS1, NS2
      REAL(8):: FACT, RSUM2, FACT1, FACT2, FACT3
      real(8):: double_count, RSUM3, RSUM_B2, RSUM_sum, RSUM_TT
      real(8),dimension(NTHMAX,NPMAX):: FNSB_B2_VOLP
      real(8),dimension(NTHMAX,NPSTART:NPEND):: FNSB_temp
      real:: gut1, gut2, gut3, gut4, gut5, gut6
      double precision:: cs_tt

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)
      NS1=NS_NSB(NSB1)
      NS2=NS_NSB(NSB2)

      IF(NSB1*NSB2.ne.0)THEN
         double_count=1.D0
         IF(ID.eq.1.or.ID.eq.2.or.ID.eq.5) double_count=0.5D0

         FACT = RNFD0(NSB1)*RNFD0(NSB2)*1.D20*double_count

         DO NP1=NPSTART,NPEND
            DO NTH1=1,NTHMAX
               RATE_NF_D1(NTH1,NP1,NR,ID)=0.D0
               RATE_NF_D2(NTH1,NP1,NR,ID)=0.D0
            END DO
         END DO

!     comm fnsb for nsb2
         call GUTIME(gut1)
         DO NP1=NPSTART, NPEND
            DO NTH1=1,NTHMAX
               FNSB_temp(nth1,np1)=FNSB(nth1,np1,nr,nsb2)*VOLP(NTH1,NP1,NS2)*RLAMDAG(NTH1,NR)*RFSADG(NR)
            END DO
         END DO
         CALL mtx_set_communicator(comm_np)
         CALL mtx_allgather_real8(FNSB_temp,NTHMAX*(NPEND-NPSTART+1),FNSB_B2_VOLP)
         call GUTIME(gut2)

         RSUM3=0.D0
         RSUM_TT=0.D0
         RSUM_sum=0.D0
         DO NP2=1,NPMAX
            DO NTH2=1,NTHMAX
               FACT2 = FNSB_B2_VOLP(NTH2,NP2)
               DO NP1=NPSTART,NPEND
                  DO NTH1=1,NTHMAX
                     FACT3 = SIGMAV_NF(NTH1,NP1,NTH2,NP2,ID) * FACT
                     
                     RATE_NF_D1(NTH1,NP1,NR,ID) = RATE_NF_D1(NTH1,NP1,NR,ID) &
                          +                     FNSB(NTH1,NP1,NR,NSB1) &
                          * FACT2 &
                          * FACT3
                     
                     IF(PM(NP1,NS1).ge.pmax_bb(NS1).and.PM(NP2,NS2).ge.pmax_bb(NS2))THEN
                        FACT1 = VOLP(NTH1,NP1,NS1)*FNSB(NTH1,NP1,NR,NSB1)*RLAMDAG(NTH1,NR)*RFSADG(NR)
                        RSUM3 = RSUM3 + FACT3*FACT1*FACT2
                     ELSEIF(PM(NP1,NS1).le.pmax_tt(NS1).and.PM(NP2,NS2).le.pmax_tt(NS2))THEN
                        FACT1 = VOLP(NTH1,NP1,NS1)*FNSB(NTH1,NP1,NR,NSB1)*RLAMDAG(NTH1,NR)*RFSADG(NR)
                        RSUM_TT = RSUM_TT + FACT3*FACT1*FACT2
                     END IF
                  END DO
               END DO
            END DO
         END DO

         CALL mtx_allreduce1_real8(RSUM3,3,RSUM_sum,vloc) ! integrate np1, nth1
         RSUM3=RSUM_sum
         call GUTIME(gut3)
!

!     comm fnsb for nsb2
         DO NP1=NPSTART, NPEND
            DO NTH1=1,NTHMAX
               FNSB_temp(nth1,np1)=FNSB(nth1,np1,nr,nsb2)
            END DO
         END DO
         CALL mtx_set_communicator(comm_np)
         CALL mtx_allgather_real8(FNSB_temp,NTHMAX*(NPEND-NPSTART+1),FNSB_B2_VOLP)

         call GUTIME(gut4)

         DO NP2=1,NPMAX ! not devide
            DO NTH2=1,NTHMAX
               RSUM_B2=0.D0
               RSUM_sum=0.D0
               DO NP1=NPSTART,NPEND
                  DO NTH1=1,NTHMAX
                     FACT1 = VOLP(NTH1,NP1,NS1)*FNSB(NTH1,NP1,NR,NSB1)*RLAMDAG(NTH1,NR)*RFSADG(NR)
                     FACT3 = SIGMAV_NF(NTH1,NP1,NTH2,NP2,ID) * FACT
                     
                     RSUM_B2 = RSUM_B2 &
                          + FACT1 &
                          *                     FNSB_B2_VOLP(NTH2,NP2) &
                          * FACT3
                  END DO
               END DO
               CALL mtx_allreduce1_real8(RSUM_B2,3,RSUM_sum,vloc) ! integrate np1, nth1
               
               RATE_NF_D2(NTH2,NP2,NR,ID) = RSUM_sum
            END DO
         END DO
         call GUTIME(gut5)
         
!     int RATE_NF_D1 = RATE_NF
         RSUM2=0.D0
         DO NP1=NPSTART,NPEND
            DO NTH1=1,NTHMAX
               RSUM2 = RSUM2 &
                    + RATE_NF_D1(NTH1,NP1,NR,ID) *VOLP(NTH1,NP1,NS1)*RLAMDAG(NTH1,NR)*RFSADG(NR)
            END DO
         END DO
         CALL mtx_allreduce1_real8(RSUM2,3,RSUM_sum,vloc) ! integrate np1, nth1
         RSUM2=RSUM_sum

         CALL mtx_reset_communicator
         call GUTIME(gut6)

!         IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NRANK.eq.0)THEN
!            WRITE(6,'(A,E14.6,3I4)') '|-NF_REACTION_RATE:', TIMEFP, nrank, comm_nr%rank, comm_nsa%rank
!            WRITE(6,'(A,I5,A,2I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
!                 ,ID,':  ',NSB1,NSB2,' -> ',NSA1,NSA2
!            WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,      ENG1_NF,      RATE_NF,   RATE_NF_BB"
!            WRITE(6,'("  ",4I5,1PE14.6,1PE12.4, 1P4E14.6)') ID,NR,NSB1,NSB2,RSUM2/FACT,ENG1_NF(ID), RSUM2, RSUM3
!         END IF
         RATE_NF(NR,ID) = RSUM2
         IF(MODELS_bt.eq.2)THEN
            RATE_NF_BT(NR,ID) = RSUM2 - RSUM3 - RATE_NF_TT(NR,ID)
         END IF

      END IF

      end SUBROUTINE NF_REACTION_RATE
!===========================================================
      SUBROUTINE CHECK_NF_RATE_ON_STDIO(NR,ID)

      IMPLICIT NONE
      integer,intent(in):: NR, ID
      integer:: NSB1, NSB2, NSA1, NSA2, NS1, NS2
      double precision:: FACT, double_count
      
      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)
      NS1=NS_NSB(NSB1)
      NS2=NS_NSB(NSB2)

      IF(NSB1*NSB2.ne.0)THEN
         double_count=1.D0
         IF(ID.eq.1.or.ID.eq.2.or.ID.eq.5) double_count=0.5D0
         FACT = RNFD0(NSB1)*RNFD0(NSB2)*1.D20*double_count

         IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NRANK.eq.0)THEN
            WRITE(6,'(A,E14.6,3I4)') '|-NF_REACTION_RATE:', TIMEFP, nrank, comm_nr%rank, comm_nsa%rank
            WRITE(6,'(A,I5,A,2I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
                 ,ID,':  ',NSB1,NSB2,' -> ',NSA1,NSA2
            WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,      ENG1_NF,      RATE_NF,   RATE_NF_BT,   RATE_NF_TT"
            WRITE(6,'("  ",4I5,1PE14.6,1PE12.4, 1P4E14.6)') ID,NR,NSB1,NSB2,RATE_NF(NR,ID)/FACT,&
                 ENG1_NF(ID), RATE_NF(NR,ID), RATE_NF_BT(NR,ID),RATE_NF_TT(NR,ID)
         END IF
      END IF

      END SUBROUTINE CHECK_NF_RATE_ON_STDIO
!===========================================================

      end MODULE fpnfrr
