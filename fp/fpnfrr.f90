!     $Id$

! **********************************
!      Nuclear Fusion Reaction Rate
! **********************************

      MODULE fpnfrr

      USE fpcomm
      USE fpinit

      PUBLIC NF_REACTION_COEF
      PUBLIC NF_REACTION_RATE

      PRIVATE

      REAL(8),dimension(1:5,1:6):: Duane_Coef
      DATA Duane_Coef/ &
             46.097D0,  372.D0,4.360D-4,1.220D0,  0.D0, & ! D-D reaction 1 ID=1
             47.880D0,  482.D0,3.080D-4,1.177D0,  0.D0, & ! D-D reaction 2 ID=2
             45.950D0,50200.D0,1.368D-2,1.076D0,409.D0, & ! D-T reaction   ID=3
             89.270D0,25900.D0,3.980D-3,1.297D0,647.D0, & ! D-He3 reaction ID=4
             38.390D0,  448.D0,1.020D-3,2.090D0,  0.D0, & ! T-T reaction   ID=5
            123.100D0,11250.D0,    0.D0,   0.D0,  0.D0/   ! T-He3 reaction ID=6
      INTEGER:: NSA_H,NSA_D,NSA_T,NSA_HE3,NSA_HE4,NSB_D,NSB_T,NSB_HE3
!      INTEGER:: NSB1_NF(6),NSB2_NF(6)

      contains

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
         SIGMAV_PHI=SIGMA_E(E,ID)*SQRT(E)
      ENDIF
      RETURN
      END FUNCTION SIGMAV_PHI

!---- gyro averaged Sigma v : multiplied by SQRT(RMASS/2) ----

      FUNCTION SIGMAV_E(E0L,E1L,ID)
      IMPLICIT NONE
      REAL(8):: SIGMAV_E
      REAL(8),INTENT(IN):: E0L,E1L
      INTEGER,INTENT(IN):: ID
      INTEGER,PARAMETER:: NPHIMAX=50 ! Number of phi mesh
      INTEGER:: NPHI
      REAL(8):: PHI,E,RSUM,DPHI

      DPHI=PI/NPHIMAX
      RSUM=0.D0
       DO NPHI=1,NPHIMAX
         PHI=DPHI*(NPHI-0.5D0)
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
      REAL(8):: VPARA1, VPARA2, VPERP1, VPERP2, V_MS
      REAL(8):: RED_MASS

      VPARA1=PM(NP1,NSB1)*COSM(NTH1)*VTFD0(NSB1)
      VPARA2=PM(NP2,NSB2)*COSM(NTH2)*VTFD0(NSB2)
      VPERP1=PM(NP1,NSB1)*SINM(NTH1)*VTFD0(NSB1)
      VPERP2=PM(NP2,NSB2)*SINM(NTH2)*VTFD0(NSB2)
!        MEAN SQUARE VELOCITY
      V_MS=VPERP1**2+VPERP2**2+(VPARA1-VPARA2)**2
!        REDUCED MASS
!      RED_MASS=AMFD(NSB1)*AMFD(NSB2)/(AMFD(NSB1)+AMFD(NSB2))
      RED_MASS=AMFD(NSB1)

      E0L=0.5D0*RED_MASS*V_MS/(AEE*1.D3)
      E1L=RED_MASS*VPERP1*VPERP2/(AEE*1.D3)
      RETURN
      END SUBROUTINE RELATIVE_ENERGY

!---- Calculation of NF reaction rate coefficient ----

      SUBROUTINE NF_REACTION_COEF

      USE fpcomm
      IMPLICIT NONE
      INTEGER,PARAMETER:: NEMAX=101 ! Number of energy mesh
      REAL(8),DIMENSION(NEMAX):: E0A,E1A
      REAL(8),DIMENSION(NEMAX,NEMAX):: SIGMAVA,FX,FY,FXY
      REAL(8),DIMENSION(4,4,NEMAX,NEMAX):: USV
      INTEGER:: NSA,NSB,NS,NSB1,NSB2,ID,NE0,NE1,NTH1,NTH2,NP1,NP2,IERR
      REAL(8):: RED_MASS,V1MAX,V2MAX,E0MAX,E1MAX,DELE0,DELE1,SUM,E0L,E1L
      REAL(8):: E3

!---- identify possible fusion reactions ----
!     ---- identify related particle species ----

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

      DO ID=1,6
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

         NSB1=NSB1_NF(ID)
         NSB2=NSB2_NF(ID)

!     ---- identify whether reaction can be described or not ----

         IF(NSB1.NE.0.AND.NSB2.NE.0.AND.NSA1_NF(ID).NE.0) THEN
         IF(MODELS.EQ.2) THEN

!         WRITE(6,*) 'NF_REACTION_COEF:'
!         WRITE(6,'(A,5I5)') 'ID,NSB1,NSB2,NSA1_NF(ID),NSA2_NF(ID)=', &
!                             ID,NSB1,NSB2,NSA1_NF(ID),NSA2_NF(ID)

         RED_MASS=AMFD(NSB1)*AMFD(NSB2)/(AMFD(NSB1)+AMFD(NSB2))

!     ---- calculate reaction rate table for spline ----
!          - energy mesh is equally spaced in SQRT(E), not E

         V1MAX=VTFD0(NSB1)*PMAX(NSB1)
         V2MAX=VTFD0(NSB2)*PMAX(NSB2)
!         WRITE(*,*) V1MAX**2*AMFD(NSB1)/(AEE*1.D3)
!         E0MAX=0.5D0*RED_MASS*(V1MAX+V2MAX)**2 /(AEE*1.D3)*2.D0
!         E1MAX=      RED_MASS* V1MAX*V2MAX     /(AEE*1.D3)*2.D0
         E0MAX= 0.5D0*AMFD(NSB1)*(V1MAX+V2MAX)**2 /(AEE*1.D3)
         E1MAX=       AMFD(NSB1)* V1MAX*V2MAX     /(AEE*1.D3)
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
                               *SQRT(2.D0*AEE*1.D3/AMFD(NSB1))
!                               *SQRT(2.D0*AEE*1.D3/RED_MASS) &!/SQRT(2.D0) &
!                               *SQRT(RED_MASS/AMFD(NSB1))
            ENDDO
         ENDDO

!         CALL PAGES
!         CALL GDD2D(0,E0A,E1A,SIGMAVA,NEMAX,NEMAX,NEMAX,'/sigmav/',0,0,0)
!         CALL PAGEE

!     ---- calculate spline coefficients ----

         CALL SPL2D(E0A,E1A,SIGMAVA,FX,FY,FXY,USV,NEMAX,NEMAX,NEMAX,0,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX NF_REACTION_COEF: SPL2D: IERR=',IERR

!     ---- calculate reaction rate on 4D momentum mesh ----

         DO NP1=1,NPMAX
            DO NTH1=1,NTHMAX
               DO NP2=1,NPMAX
                  DO NTH2=1,NTHMAX
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

 !        DO NP1=1,NPMAX
 !        DO NTH1=1,NTHMAX
 !        DO NP2=1,NPMAX
 !        DO NTH2=1,NTHMAX
 !        E3=(PTFP0(NSB1)*PM(NP1,NSB1))**2/AMFP(NSB1)/(AEE*1.D3)
 !        SIGMAV_NF(NTH1,NP1,NTH2,NP2,ID)=SIGMA_E(E3,ID)*SQRT(E3) &
 !                              *SQRT(AEE*1.D3/AMFP(NSB1))
 !        END DO
 !        END DO
 !        END DO
 !        END DO
 !        write(*,*)E3


         ENDIF
         ENDIF


      ENDDO

      RETURN
      END SUBROUTINE NF_REACTION_COEF

!===========================================================
!     CALCULATION OF INTEGRAL SIGMA*V*f(p_a)*f(p_b)
!===========================================================

      SUBROUTINE NF_REACTION_RATE(NR,ID)

      USE fpcomm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,ID
      INTEGER:: NSB1, NSB2, NSA1, NSA2, NP1, NP2, NTH1, NTH2
      REAL(8):: RSUM, FACT, RSUM2, FACT1, FACT2, FACT3

      NSB1=NSB1_NF(ID)
      NSB2=NSB2_NF(ID)
      NSA1=NSA1_NF(ID)
      NSA2=NSA2_NF(ID)

      FACT = RNFD0(NSB1)*RNFD0(NSB2)*1.D20

      DO NP1=1,NPMAX
         DO NTH1=1,NTHMAX
            RATE_NF_D1(NTH1,NP1,NR,ID)=0.D0
            RATE_NF_D2(NTH1,NP1,NR,ID)=0.D0
         END DO
      END DO

      DO NP2=1,NPMAX
      DO NTH2=1,NTHMAX
         RSUM=0.D0
         FACT2 = VOLP(NTH2,NP2,NSB2)*FNS(NTH2,NP2,NR,NSB2)*RLAMDA(NTH2,NR)/RFSADG(NR)
         DO NP1=1,NPMAX
         DO NTH1=1,NTHMAX
            FACT1 = VOLP(NTH1,NP1,NSB1)*FNS(NTH1,NP1,NR,NSB1)*RLAMDA(NTH1,NR)/RFSADG(NR)
            FACT3 = SIGMAV_NF(NTH1,NP1,NTH2,NP2,ID) * FACT

            RATE_NF_D1(NTH1,NP1,NR,ID) = RATE_NF_D1(NTH1,NP1,NR,ID) &
                 +                     FNS(NTH1,NP1,NR,NSB1) &
                 * FACT2 &
                 * FACT3
            
            RATE_NF_D2(NTH2,NP2,NR,ID) = RATE_NF_D2(NTH2,NP2,NR,ID) &
                 + FACT1 &
                 *                     FNS(NTH2,NP2,NR,NSB2) &
                 * FACT3
         END DO
         END DO
      END DO
      END DO

!     int RATE_NF_D1 = RATE_NF
      RSUM2=0.D0
      DO NP1=1,NPMAX
      DO NTH1=1,NTHMAX
         RSUM2 = RSUM2 &
              + RATE_NF_D1(NTH1,NP1,NR,ID) *VOLP(NTH1,NP1,NSB1) &
              * RLAMDA(NTH1,NR)/RFSADG(NR)
      END DO
      END DO

      IF((N_IMPL.eq.0.or.N_IMPL.gt.LMAXFP).and.NR.eq.1)THEN
         WRITE(6,*) '|-NF_REACTION_RATE:'
         WRITE(6,'(A,3I5,A,2I5)') '   |-ID,NSB1,NSB2 -> NSA1,NSA2=' &
              ,ID,NSB1,NSB2,' -> ',NSA1,NSA2
         WRITE(6, *) "  |-ID,  NR,NSB1,NSB2,  <sigma*v>,      ENG1_NF"
         WRITE(6,'("  ",4I5,1PE14.6,1PE12.4)') ID,NR,NSB1,NSB2,RSUM2/FACT,ENG1_NF(ID)
      END IF
      RATE_NF(NR,ID) = RSUM2

      end SUBROUTINE NF_REACTION_RATE
!===========================================================

      end MODULE fpnfrr

