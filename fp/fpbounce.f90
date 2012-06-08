!     $Id$

! *****************************
!     PREPARATION OF PARAMETERS FOR BOUNCE AVERAGE
! *****************************
      MODULE fpbounce

      USE fpcomm
      USE fpinit
      USE equnit_mod

      INTEGER,parameter:: mmax=100

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_BOUNCE_PARAM(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH

!     ON RM(NR)
      IF(NR.ne.NRMAX+1)THEN
         CALL SET_ETAMG(NR)

!         CALL SET_RLAMDA_DE(NR)
!         CALL SET_RLAMDA_ELL(NR)
         CALL SET_RLAMDA(NR) ! NAVMAX
!         CALL SET_RLAMDA_TPB(NR) ! Kileen
!         CALL SET_RLAMDA_TPB2(NR) ! Kileen with correction
         CALL SET_RLAMDA_TPB3(NR) ! NAVMAX
      END IF

!     ON RG(NR)
      IF(NR.eq.1)THEN
         DO NTH=1,NTHMAX
            RLAMDA_G(NTH,NR)=1.D0
         END DO
      ELSEIF(NR.eq.NRMAX+1)THEN
         CALL SET_ETAMG_GMAX(NR)

!         CALL SET_RLAMDA_DE_GMAX(NR)
!         CALL SET_RLAMDA_ELL_GMAX(NR)
         CALL SET_RLAMDA_GMAX(NR)

!         CALL SET_RLAMDA_TPB_GMAX(NR)
!         CALL SET_RLAMDA_TPB2_G(NR) ! Killeen with correction
         CALL SET_RLAMDA_TPB3_G(NR)
      ELSE
         CALL SET_ETAMG_G(NR)

!         CALL SET_RLAMDA_DE_G(NR)
!         CALL SET_RLAMDA_ELL_G(NR)
         CALL SET_RLAMDA_G(NR)

!         CALL SET_RLAMDA_TPB_G(NR) ! Killeen
!         CALL SET_RLAMDA_TPB2_G(NR) ! Killeen with correction
         CALL SET_RLAMDA_TPB3_G(NR)
      END IF

      END SUBROUTINE SET_BOUNCE_PARAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET BOUNCE POINT phi_b/2=ETAM ON RM(NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      DOUBLE PRECISION:: EPSL, FACT, A1

      EPSL=EPSRM2(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      
      DO NTH=1,ITL(NR)
         ETAM(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITL(NR)+1,ITU(NR)-1
         A1=FACT*COSM(NTH)**2
         ETAM(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITU(NR),NTHMAX
         ETAM(NTH,NR)=PI*0.5D0
      ENDDO
!
      DO NTH=1,ITL(NR)
         ETAG(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITL(NR)+1,ITU(NR)
         A1=FACT*COSG(NTH)**2
         ETAG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITU(NR)+1,NTHMAX+1
         ETAG(NTH,NR)=PI*0.5D0
      ENDDO

      END SUBROUTINE SET_ETAMG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET BOUNCE POINT phi_b/2=ETAM ON RG(NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG_G(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      DOUBLE PRECISION:: EPSL, FACT, A1

      EPSL=EPSRG2(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      
      DO NTH=1,ITLG(NR)
         ETAM_G(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITLG(NR)+1,ITUG(NR)-1
         A1=FACT*COSM(NTH)**2
         ETAM_G(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITUG(NR),NTHMAX
         ETAM_G(NTH,NR)=PI*0.5D0
      ENDDO
!
      DO NTH=1,ITLG(NR)
         ETAG_G(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITLG(NR)+1,ITUG(NR)
         A1=FACT*COSG(NTH)**2
         ETAG_G(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITUG(NR)+1,NTHMAX+1
         ETAG_G(NTH,NR)=PI*0.5D0
      ENDDO

      END SUBROUTINE SET_ETAMG_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOUNDARY LAYER THM = THETA_B
      SUBROUTINE SET_RLAMDA_TPB(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub,sinb
      DOUBLE PRECISION:: k1, k2, sum, K_1_x
      DOUBLE PRECISION,dimension(0:m):: alpha, beta
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

!     define DEL2T
      mub = SQRT(2.D0*EPSRM2(NR) / (1.D0 + EPSRM2(NR) ) )
      sinb= SQRT( (1.D0-EPSRM2(NR))/(1.D0+EPSRM2(NR)) )
      DEL2T1 = COSG(ITL(NR))/mub -1.D0
      DEL2T2 =-COSG(ITL(NR)+1)/mub +1.D0
      DEL2T = 0.5D0 *( DEL2T1 + DEL2T2 )

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM = 0.D0
      DO i = 0, m
         SUM = SUM + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
!     
      K_1_x  = mub*sinb*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
!
      RLAMDA(ITL(NR),NR) = 4.D0*QLM(NR)*RR*(K_1_x - sum)*DEL2T/(sinb*DELTH)
      RLAMDA(ITU(NR),NR) = RLAMDA(ITL(NR),NR)

      END SUBROUTINE SET_RLAMDA_TPB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB_G(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub
      DOUBLE PRECISION:: k1, k2, sum, K_1_x
      DOUBLE PRECISION,dimension(0:m):: alpha, beta
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b
      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

!     define DEL2T
      mub = SQRT(2.D0*EPSRG2(NR) / (1.D0 + EPSRG2(NR) ) )
      DEL2T1 = COSG(ITLG(NR))/mub -1.D0
      DEL2T2 =-COSG(ITLG(NR)+1)/mub +1.D0
      DEL2T = 0.5D0 *( DEL2T1 + DEL2T2 )

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM = 0.D0
      DO i = 0, m
         SUM = SUM + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
!     
      K_1_x  = mub*SQRT(1.D0-mub**2)*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
!
      RLAMDA_G(ITLG(NR),NR) = 2.D0*QLG(NR)*RR*(K_1_x - sum)/SINM(ITLG(NR))
      RLAMDA_G(ITUG(NR),NR) = RLAMDA_G(ITLG(NR),NR)

      END SUBROUTINE SET_RLAMDA_TPB_G

!-----------------------------
      SUBROUTINE RECURRENCE_ALPHA(m,alpha)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: m
      DOUBLE PRECISION,dimension(0:m),INTENT(OUT):: alpha
      INTEGER:: i
      DOUBLE PRECISION,parameter:: a0=1.D0, a1=-0.5D0
      
      alpha(0)=a0
      alpha(1)=a1

      DO i = 2, m
         alpha(i) = alpha(i-1)*(2.D0*i-3.D0)/(2.D0*i)
      END DO

      END SUBROUTINE RECURRENCE_ALPHA
!-----------------------------
      SUBROUTINE RECURRENCE_BETA(m,beta)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: m
      DOUBLE PRECISION,dimension(0:m),INTENT(OUT):: beta
      INTEGER:: i
      DOUBLE PRECISION,parameter:: b0=0.D0, b1=1.D0
      
      beta(0)=b0
      beta(1)=b1

      DO i = 2, m
         beta(i) = ( 2.D0*(2.D0*i-2.D0)*beta(i-1) - (2.D0*i-3.D0)*beta(i-2) ) &
              /(2.D0*i-1.D0)
      END DO

      END SUBROUTINE RECURRENCE_BETA
!
!----------------------------
!
!     USE FOR PASSING REGION ONLY ON RM(NR)
      SUBROUTINE SET_RLAMDA_DE(NR)

      USE libde,ONLY: DEFT
!      USE libmtx
!      USE plprof

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      DOUBLE PRECISION:: EPSL, FACT, A1, RINT0, ES, mub
!
      IF(ITL(NR).ne.1)THEN
         NRX=NR
         DO NTH = 1, ITL(NR)-1
            NTHX = NTH
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE)
            RLAMDA(NTH,NR) = RINT0 * ( QLM(NR)*RR ) * 2.D0
         END DO
!     symmetry
         DO NTH=1,ITL(NR)-1
            RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
         ENDDO
      END IF

      END SUBROUTINE SET_RLAMDA_DE
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) )
      mu02 = COSM(NTHX)**2
      A0 = ETAM(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DEM(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) )
      mu02 = COSG(NTHX)**2
      A0 = ETAG(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DEM = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DEM
!
!----------------------------
!
!     USE FOR PASSING REGION ONLY ON RG(NR)
      SUBROUTINE SET_RLAMDA_DE_G(NR)

      USE libde,ONLY: DEFT

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      DOUBLE PRECISION:: EPSL, FACT, A1, RINT0, ES
!
      IF(ITLG(NR).ne.1)THEN
         NRX=NR
         DO NTH = 1, ITLG(NR)-1
            NTHX = NTH
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_G)
            RLAMDA_G(NTH,NR) = RINT0 * ( QLG(NR)*RR ) * 2.D0
         END DO

!     symmetry
         DO NTH=1,ITLG(NR)-1
            RLAMDA_G(NTHMAX-NTH+1,NR)=RLAMDA_G(NTH,NR)
         ENDDO
      END IF

    END SUBROUTINE SET_RLAMDA_DE_G
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_G(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRG2(NRX) / (1.D0 + EPSRG2(NRX) )
      mu02 = COSM(NTHX)**2
      A0 = ETAM_G(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_G = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_G
!
!-----------------------------
!     
!    USE FOR TRAPPED REGION ONLY ON RM(NR)
     SUBROUTINE SET_RLAMDA_ELL(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: z, sum, temp, mub, a1,a2
      DOUBLE PRECISION,dimension(0:m):: alpha, j2m

      mub=SQRT(2.D0*EPSRM2(NR)/(1.D0+EPSRM2(NR)))
      IF(ITL(NR).ne.NTHMAX/2)THEN
         DO NTH = ITL(NR)+1, NTHMAX/2
            z = COSM(NTH)/mub
            call RECURRENCE_ALPHA(m,alpha)
            call RECURRENCE_J(m,z,j2m)
            SUM = 0.D0
            DO i = 0, m
               SUM = SUM + alpha(i)*( mub )**(2*i)*j2m(i)
            END DO
            RLAMDA(NTH,NR) = SUM * ( QLM(NR)*RR ) * 2.D0
         END DO
! symmetry
         DO NTH=ITL(NR)+1,NTHMAX/2
            RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
         ENDDO
      END IF


      END SUBROUTINE SET_RLAMDA_ELL
!
!-----------------------------
!     ON RG(NR)
      SUBROUTINE SET_RLAMDA_ELL_G(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: z, sum, mub, a1, a2
      DOUBLE PRECISION,dimension(0:m):: alpha, j2m

      mub=SQRT(2.D0*EPSRG2(NR)/(1.D0+EPSRG2(NR)))
!      mub = COSM(ITLG(NR))
      IF(ITLG(NR).ne.NTHMAX/2)THEN
         DO NTH = ITLG(NR)+1, NTHMAX/2
            z = COSM(NTH)/mub
            call RECURRENCE_ALPHA(m,alpha)
            call RECURRENCE_J(m,z,j2m)
            SUM = 0.D0
            DO i = 0, m
               SUM = SUM + alpha(i)*( mub )**(2*i)*j2m(i)
            END DO
            RLAMDA_G(NTH,NR) = SUM * ( QLG(NR)*RR ) * 2.D0
         END DO
!     BOUNDARY
!         RLAMDA_G(ITLG(NR),NR) = 0.5D0* &
!              ( RLAMDA_G(ITLG(NR)-1,NR) + RLAMDA_G(ITLG(NR)+1,NR) )
!     symmetry
         DO NTH=ITLG(NR)+1,NTHMAX/2
            RLAMDA_G(NTHMAX-NTH+1,NR)=RLAMDA_G(NTH,NR)
         ENDDO
      END IF

      END SUBROUTINE SET_RLAMDA_ELL_G
!-----------------------------
      SUBROUTINE RECURRENCE_J(m,z,j2m)

      USE LIBELL,only : ellfc, ellec
      IMPLICIT NONE
      INTEGER,INTENT(IN):: m
      INTEGER:: IERR, i
      DOUBLE PRECISION,INTENT(IN):: z
      DOUBLE PRECISION,dimension(0:m),INTENT(OUT):: j2m
      DOUBLE PRECISION:: J0, J2, ellK, ellE

      IF(z.ge.1.D0)THEN
!         WRITE(*,*) "ELL1", 1/z**2
         ellK = ELLFC(1/z**2,IERR)
         ellE = ELLEC(1/z**2,IERR)
!         WRITE(*,*) "ELL12", ellK, ellE
         J0 = ellK
         J2 = z**2*( ellK - ellE )
      ELSEIF(z.lt.1.D0.and.z.gt.0)THEN
!         WRITE(*,*) "ELL2", z**2
         ellK = ELLFC(z**2,IERR)
         ellE = ELLEC(z**2,IERR)
!         WRITE(*,*) "ELL22", ellK, ellE
         J0 = z* ellK
         J2 = z* ( ellK - ellE )
      ENDIF
      j2m(0) = J0
      j2m(1) = J2
      DO i = 2, m
         j2m(i) = ( (2.D0*i - 2.D0)*(1.D0 + z**2)*j2m(i-1) &
              - (2.D0*m - 3.D0)*z**2*j2m(i-2) ) &
              /(2.D0*i-1.D0)
      END DO

      END SUBROUTINE RECURRENCE_J
!
!-----------------------------
!!    CALCULATE FLUX SURFACE AVERAGE oint 1/psi ds
      SUBROUTINE SET_RFSAD(NR)

      IMPLICIT NONE
      INTEGER:: NR

      RFSADG(NR)=PI*QLM(NR)*RR/(1.D0+EPSRM2(NR))
      RFSAD_GG(NR)=PI*QLG(NR)*RR/(1.D0+EPSRG2(NR))
         
      END SUBROUTINE SET_RFSAD
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA(NR)

      USE libde,ONLY: DEFT
      USE libmtx
      USE plprof
!      USE fpbroadcast
!      USE fpwrin
!      USE fpwmin

      IMPLICIT NONE
      integer:: NR, NTH, NG
      real(kind8):: EPSL, FACT, A1, RINT0, ES, DELH
      real(kind8):: SUML, ETAL, X, PSIB, PCOS, RINT2, SUML2, RINTL

      NRX=NR
      DO NTH=1,NTHMAX/2
!         NTHX=NTH
         IF(NTH.ne.ITL(NR)) THEN
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            suml=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRM2(NR)*COS(ETAL)*RR
               PSIB=(1.D0+EPSRM2(NR))/(1.D0+X/RR)
               IF(COSM(NTH).ge.0.D0)THEN
                  PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               ELSE
                  PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
               END IF
               suml=suml+1.D0/PCOS
            END DO
            RINT0 = suml*DELH
            RLAMDA(NTH,NR)=RINT0*ABS(COSM(NTH)) * ( QLM(NR)*RR )! lambda = v_0*|cos\theta_0|*tau_B
!         RLAMDA(NTH,NR)=RINT0*ABS(COSM(NTH)) ! lambda = v_0*|cos\theta_0|*tau_B / qR
!            CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A)
!            RLAMDC(NTH,NR)=RINT2/(PI*(1.D0+EPSRM(NR))*(COSG(NTH)))
            RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
!            RLAMDC(NTHMAX-NTH+2,NR)=RLAMDC(NTH,NR)
         ENDIF
      ENDDO
!      RLAMDC(NTHMAX/2+1,NR)=0.D0

      END SUBROUTINE SET_RLAMDA
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA_G(NR)

      USE libde,ONLY: DEFT

      IMPLICIT NONE
      integer:: NR, NTH, NG
      real(kind8):: EPSL, FACT, A1, RINT0, ES, DELH
      real(kind8):: SUML, ETAL, X, PSIB, PCOS, RINT2, SUML2, RINTL

      NRX=NR
      DO NTH=1,NTHMAX/2
         NTHX=NTH
         IF(NTH.ne.ITLG(NR)) THEN
            DELH=2.D0*ETAM_G(NTH,NR)/NAVMAX
            suml=0.D0
            SUML2=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRG2(NR)*COS(ETAL)*RR
               PSIB=(1.D0+EPSRG2(NR))/(1.D0+X/RR)
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               suml=suml+1.D0/PCOS
            END DO
            RINT0 = suml *DELH
!            CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A_G)
            RLAMDA_G(NTH,NR)=RINT0*ABS(COSM(NTH)) * ( QLG(NR)*RR )! lambda = v_0*|cos\theta_0|*tau_B
!            RLAMDC_G(NTH,NR)=RINT2/(PI*(1.D0+EPSRG(NR))*(COSG(NTH)))

            RLAMDA_G(NTHMAX-NTH+1,NR)=RLAMDA_G(NTH,NR)
!            RLAMDC_G(NTHMAX-NTH+2,NR)=RLAMDC_G(NTH,NR)
         ENDIF
      ENDDO
!      RLAMDC_G(NTHMAX/2+1,NR)=0.D0

      END SUBROUTINE SET_RLAMDA_G
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA_GMAX(NR)

      IMPLICIT NONE
      integer:: NR, NTH, NG
      real(kind8):: EPSL, FACT, A1, RINT0, ES, DELH
      real(kind8):: SUML, ETAL, X, PSIB, PCOS, RINT2, SUML2, RINTL

      NRX=NR
      DO NTH=1,NTHMAX/2
         NTHX=NTH
         IF(NTH.ne.ITLG(NR)) THEN
            DELH=2.D0*ETAM_GG(NTH,NR)/NAVMAX
            suml=0.D0 
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRG2(NR)*COS(ETAL)*RR
               PSIB=(1.D0+EPSRG2(NR))/(1.D0+X/RR)
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               suml=suml + 1.D0/PCOS
            END DO
            RINT0 = suml*DELH
            RLAMDA_GG(NTH,NR)=RINT0*ABS(COSM(NTH)) * ( QLG(NR)*RR )! lambda = v_0*|cos\theta_0|*tau_B
            RLAMDA_GG(NTHMAX-NTH+1,NR)=RLAMDA_GG(NTH,NR)
         ENDIF
      ENDDO

      END SUBROUTINE SET_RLAMDA_GMAX

! ***************************************************************

!                       SET OF INTEGRAND

! ***************************************************************

! ============================================================

      REAL*8 FUNCTION  FPFN0U(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRM(NRX))*SINM(NTHX)**2
      FPFN0U=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN0U

! ============================================================

      REAL*8 FUNCTION  FPFN0U_G(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM_G(NTHX,NRX)
      A1=1.D0+EPSRG(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRG(NRX))*SINM(NTHX)**2
      FPFN0U_G=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN0U_G

! ============================================================

      REAL*8 FUNCTION FPFN0T(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=2.D0*EPSRM(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN0T=A0*SQRT(A1/A2)

      RETURN
      END FUNCTION FPFN0T

! ============================================================

      REAL*8 FUNCTION FPFN1A(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      FPFN1A=A0*SQRT(A1)

      RETURN
      END FUNCTION FPFN1A

! ============================================================

      REAL*8 FUNCTION FPFN2A(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSRM(NRX))*SING(NTHX)**2
      FPFN2A=A0*SQRT(A1*A2)

      RETURN
      END FUNCTION FPFN2A
! ============================================================

      REAL*8 FUNCTION FPFN2A_G(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG_G(NTHX,NRX)
      A1=1.D0+EPSRG(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSRG(NRX))*SING(NTHX)**2
      FPFN2A_G=A0*SQRT(A1*A2)

      RETURN
      END FUNCTION FPFN2A_G
!--------------------------------
      REAL*8 FUNCTION  FPFN2U(X,XM,XP)
                           
      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRM(NRX))*SING(NTHX)**2
      FPFN2U=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN2U
!------------------------------
      REAL*8 FUNCTION FPFN2T(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=2.D0*EPSRM(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN2T=A0*SQRT(A1/A2)

      RETURN
      END FUNCTION FPFN2T

!-------------------------------------------------

      REAL*8 FUNCTION FPFN2A_GG(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG_G_GL(NTHX,NRX)
      A1=1.D0+EPSRG(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSRG(NRX))*SING(NTHX)**2
      FPFN2A_GG=A0*SQRT(A1*A2)

      RETURN
      END FUNCTION FPFN2A_GG
! ============================================================
      REAL*8 FUNCTION  FPFN0U_GG(X,XM,XP)

      IMPLICIT NONE
      real(8),INTENT(IN):: X, XM, XP
      real(8):: XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM_GG(NTHX,NRX)
      A1=1.D0+EPSRG(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRG(NRX))*SINM(NTHX)**2
      FPFN0U_GG=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN0U_GG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG_GMAX(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      DOUBLE PRECISION:: EPSL, FACT, A1

      EPSL=EPSRG2(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      
      DO NTH=1,ITLG(NR)
         ETAM_GG(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITLG(NR)+1,ITUG(NR)-1
         A1=FACT*COSM(NTH)**2
         ETAM_GG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITUG(NR),NTHMAX
         ETAM_GG(NTH,NR)=PI*0.5D0
      ENDDO
!
!      DO NTH=1,ITLG(NR)
!         ETAG_GG(NTH,NR)=PI*0.5D0
!      ENDDO
!      DO NTH=ITLG(NR)+1,ITUG(NR)
!         A1=FACT*COSG(NTH)**2
!         ETAG_GG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
!      ENDDO
!      DO NTH=ITUG(NR)+1,NTHMAX+1
!         ETAG_GG(NTH,NR)=PI*0.5D0
!      ENDDO

      END SUBROUTINE SET_ETAMG_GMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB_GMAX(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub
      DOUBLE PRECISION:: k1, k2, sum, K_1_x
      DOUBLE PRECISION,dimension(0:m):: alpha, beta
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b
      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

!     define DEL2T
      mub = SQRT(2.D0*EPSRG2(NR) / (1.D0 + EPSRG2(NR) ) )
      DEL2T1 = COSG(ITLG(NR))/mub -1.D0
      DEL2T2 =-COSG(ITLG(NR)+1)/mub +1.D0
      DEL2T = 0.5D0 *( DEL2T1 + DEL2T2 )

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM = 0.D0
      DO i = 0, m
         SUM = SUM + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
!     
      K_1_x  = mub*SQRT(1.D0-mub**2)*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
!
      RLAMDA_GG(ITLG(NR),NR) = 2.D0*QLG(NR)*RR*(K_1_x - sum)/SINM(ITLG(NR))
      RLAMDA_GG(ITUG(NR),NR) = RLAMDA_GG(ITLG(NR),NR)

      END SUBROUTINE SET_RLAMDA_TPB_GMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     USE FOR PASSING REGION ONLY ON RG(NR)
      SUBROUTINE SET_RLAMDA_DE_GMAX(NR)

      USE libde,ONLY: DEFT
      USE libmtx
      USE plprof

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      DOUBLE PRECISION:: EPSL, FACT, A1, RINT0, ES
!
      IF(ITLG(NR).ne.1)THEN
         NRX=NR
         DO NTH = 1, ITLG(NR)-1
            NTHX = NTH
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_GMAX)
            RLAMDA_GG(NTH,NR) = RINT0 * ( QLG(NR)*RR ) * 2.D0
         END DO

!     symmetry
         DO NTH=1,ITLG(NR)-1
            RLAMDA_GG(NTHMAX-NTH+1,NR)=RLAMDA_GG(NTH,NR)
         ENDDO
      END IF

      END SUBROUTINE SET_RLAMDA_DE_GMAX
!-----------------------------
!     ON RG(NR)
      SUBROUTINE SET_RLAMDA_ELL_GMAX(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: z, sum, mub, a1, a2
      DOUBLE PRECISION,dimension(0:m):: alpha, j2m

      mub=SQRT(2.D0*EPSRG2(NR)/(1.D0+EPSRG2(NR)))
      IF(ITLG(NR).ne.NTHMAX/2)THEN
         DO NTH = ITLG(NR)+1, NTHMAX/2
            z = COSM(NTH)/mub
            call RECURRENCE_ALPHA(m,alpha)
            call RECURRENCE_J(m,z,j2m)
            SUM = 0.D0
            DO i = 0, m
               SUM = SUM + alpha(i)*( mub )**(2*i)*j2m(i)
            END DO
            RLAMDA_GG(NTH,NR) = SUM * ( QLG(NR)*RR ) * 2.D0
         END DO
!     symmetry
         DO NTH=ITLG(NR)+1,NTHMAX/2
            RLAMDA_GG(NTHMAX-NTH+1,NR)=RLAMDA_GG(NTH,NR)
         ENDDO
      END IF

      END SUBROUTINE SET_RLAMDA_ELL_GMAX
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_GMAX(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRG2(NRX) / (1.D0 + EPSRG2(NRX) )
      mu02 = COSM(NTHX)**2
      A0 = ETAM_GG(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_GMAX = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_GMAX
! ============================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB2(NR)

      USE libde, only: DEFT
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub,sinb, DELTHB, z, RINT0, ES, ram_bl
      DOUBLE PRECISION:: k1, k2, sum_bl, sum_ex, K_1_x, L_bl, L_ex, L_b, THEX, DELTHEX
      DOUBLE PRECISION,dimension(0:m):: alpha, beta, j2m
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRM2(NR) / (1.D0 + EPSRM2(NR) ) )
      sinb= SQRT( (1.D0-EPSRM2(NR))/(1.D0+EPSRM2(NR)) )

      IF(COSM(ITL(NR)).le.mub)THEN ! ITL = TRAPPED
!!     define DEL2T
         DEL2T = COSG(ITL(NR))/mub -1.D0
         DELTHB =( ACOS(mub) - THG(ITL(NR)) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITL(NR)+1) - DELTHEX*0.5D0
!
         z = COS(THEX)/mub
         call RECURRENCE_ALPHA(m,alpha)
         call RECURRENCE_J(m,z,j2m)
         SUM_EX = 0.D0
         DO i = 0, m
            SUM_EX = SUM_EX + alpha(i)*( mub )**(2*i)*j2m(i)
         END DO
         L_ex = SUM_EX*( QLM(NR)*RR )*2.D0 * SIN(THEX)*DELTHEX
      ELSE ! ITL = PASSING
         DEL2T = -COSG(ITL(NR)+1)/mub +1.D0
         DELTHB = ( THG(ITL(NR)+1) - ACOS(mub) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITL(NR)) + DELTHEX*0.5D0
!
         NTHX = ITL(NR)
         NRX = NR
         CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_TPB)
         L_ex = RINT0 * ( QLM(NR)*RR ) * 2.D0 * SIN(THEX)*DELTHEX
      END IF

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM_BL = 0.D0
      DO i = 0, m
         SUM_BL = SUM_BL + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
      K_1_x  = mub*sinb*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
      L_b = 8.D0*QLM(NR)*RR*DEL2T*( K_1_x - SUM_BL ) * 0.5D0
!     L_b means lambda_B * SIN(TEHTA_B) *DELTHB

!      IF(COSM(ITL(NR)).le.mub)THEN ! ITL = TRAPPED
         RLAMDA(ITL(NR),NR) = ( L_b + L_ex )/( SINM(ITL(NR))*DELTH )
         RLAMDA(ITU(NR),NR) = RLAMDA(ITL(NR),NR)
!      END IF

      END SUBROUTINE SET_RLAMDA_TPB2
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_TPB(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02, DELTHB, THEX, mub

      mub = SQRT(2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) ) )
      DELTHB = ( THG(ITL(NRX)+1) - ACOS(mub) )*2.D0
      THEX = THG(ITL(NRX)) + (DELTH-DELTHB)*0.5D0

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) )
      mu02 = COS(THEX)**2
      A0 = ETAM(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_TPB = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_TPB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB2_G(NR)

      USE libde, only: DEFT
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub,sinb, DELTHB, z, RINT0, ES, ram_bl
      DOUBLE PRECISION:: k1, k2, sum_bl, sum_ex, K_1_x, L_bl, L_ex, L_b, THEX, DELTHEX
      DOUBLE PRECISION,dimension(0:m):: alpha, beta, j2m
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRG2(NR) / (1.D0 + EPSRG2(NR) ) )
      sinb= SQRT( (1.D0-EPSRG2(NR))/(1.D0+EPSRG2(NR)) )

      IF(COSM(ITLG(NR)).le.mub)THEN ! ITL = TRAPPED
!!     define DEL2T
         DEL2T = COSG(ITLG(NR))/mub -1.D0
         DELTHB =( ACOS(mub) - THG(ITLG(NR)) )*2.D0

!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITLG(NR)+1) - DELTHEX*0.5D0
!
         z = COS(THEX)/mub
         call RECURRENCE_ALPHA(m,alpha)
         call RECURRENCE_J(m,z,j2m)
         SUM_EX = 0.D0
         DO i = 0, m
            SUM_EX = SUM_EX + alpha(i)*( mub )**(2*i)*j2m(i)
         END DO
         L_ex = SUM_EX*( QLG(NR)*RR )*2.D0 * SIN(THEX)*DELTHEX
      ELSE ! ITL = PASSING
         DEL2T = -COSG(ITLG(NR)+1)/mub +1.D0
         DELTHB = ( THG(ITLG(NR)+1) - ACOS(mub) )*2.D0

!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITLG(NR)) + DELTHEX*0.5D0

         NTHX = ITLG(NR)
         NRX = NR
         IF(NR.eq.NRMAX+1)THEN
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_TPBGMAX)
         ELSE
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_TPBG)
         END IF

         L_ex = RINT0 * ( QLG(NR)*RR ) * 2.D0 * SIN(THEX)*DELTHEX
      END IF

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM_BL = 0.D0
      DO i = 0, m
         SUM_BL = SUM_BL + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
      K_1_x  = mub*sinb*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
      L_b = 8.D0*QLG(NR)*RR*DEL2T*( K_1_x - SUM_BL ) * 0.5D0
!     L_b means lambda_B * SIN(TEHTA_B) *DELTHB

      IF(NR.ne.NRMAX+1)THEN
         RLAMDA_G(ITLG(NR),NR) = ( L_b + L_ex )/( SINM(ITLG(NR))*DELTH )
         RLAMDA_G(ITUG(NR),NR) = RLAMDA_G(ITLG(NR),NR)
      ELSEIF(NR.eq.NRMAX+1)THEN
         RLAMDA_GG(ITLG(NR),NR) = ( L_b + L_ex )/( SINM(ITLG(NR))*DELTH )
         RLAMDA_GG(ITUG(NR),NR) = RLAMDA_GG(ITLG(NR),NR)
      END IF

!      END IF

      END SUBROUTINE SET_RLAMDA_TPB2_G
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_TPBG(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02, DELTHB, THEX, mub

      mub = SQRT(2.D0*EPSRG2(NRX) / (1.D0 + EPSRG2(NRX) ) )
      DELTHB = ( THG(ITLG(NRX)+1) - ACOS(mub) )*2.D0
      THEX = THG(ITLG(NRX)) + (DELTH-DELTHB)*0.5D0

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRG2(NRX) / (1.D0 + EPSRG2(NRX) )
      mu02 = COS(THEX)**2
      A0 = ETAM_G(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_TPBG = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_TPBG
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_TPBGMAX(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02, DELTHB, THEX, mub

      mub = SQRT(2.D0*EPSRG2(NRX) / (1.D0 + EPSRG2(NRX) ) )
      DELTHB = ( THG(ITLG(NRX)+1) - ACOS(mub) )*2.D0
      THEX = THG(ITLG(NRX)) + (DELTH-DELTHB)*0.5D0

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRG2(NRX) / (1.D0 + EPSRG2(NRX) )
      mu02 = COS(THEX)**2
      A0 = ETAM_GG(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_TPBGMAX = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_TPBGMAX

!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_TPBL(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02, DELTHB, THEX, mub

      mub = SQRT(2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) ) )
!      DELTHB = ( THG(ITL(NRX)+1) - ACOS(mub) )*2.D0
!      THEX = THG(ITL(NRX)) + (DELTH-DELTHB)*0.5D0

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) )
      mu02 = COS(THG(ITL(NRX)))**2
      A0 = ETAM(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_TPBL = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_TPBL
!---------------------------------
      DOUBLE PRECISION FUNCTION FPF_RLAMDA_DE_TPBH(X,XM,XP)

      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(IN):: X, XM, XP
      DOUBLE PRECISION:: XX, A0, A1, A2
      DOUBLE PRECISION:: mub2, mu02, DELTHB, THEX, mub

      mub = SQRT(2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) ) )
      DELTHB = ( THG(ITL(NRX)+1) - ACOS(mub) )*2.D0
      THEX = THG(ITL(NRX)) + (DELTH-DELTHB)

      XX=X
      XX=XM
      mub2 = 2.D0*EPSRM2(NRX) / (1.D0 + EPSRM2(NRX) )
      mu02 = COS(THEX)**2
      A0 = ETAM(NTHX,NRX) * 0.5D0
      A1 = 1.D0-mub2*SIN(A0*XP)**2
      A2 = 1.D0-(mub2/mu02)*SIN(A0*XP)**2

      FPF_RLAMDA_DE_TPBH = A0*SQRT(A1/A2)!*0.5D0

      RETURN

      END FUNCTION FPF_RLAMDA_DE_TPBH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB3(NR)

      USE libde, only: DEFT
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i, NG
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub,sinb, DELTHB, z, RINT0, ES, ram_bl
      DOUBLE PRECISION:: k1, k2, sum_bl, sum_ex, K_1_x, L_bl, L_ex, L_b, THEX, DELTHEX
      DOUBLE PRECISION,dimension(0:m):: alpha, beta, j2m
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b
      DOUBLE PRECISION:: phib, DELH, suml, ETAL, PSIB, X, PCOS, FACT

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRM2(NR) / (1.D0 + EPSRM2(NR) ) )
      sinb= SQRT( (1.D0-EPSRM2(NR))/(1.D0+EPSRM2(NR)) )

      IF(COSM(ITL(NR)).le.mub)THEN ! ITL = TRAPPED
!!     define DEL2T
         DEL2T = COSG(ITL(NR))/mub -1.D0
         DELTHB =( ACOS(mub) - THG(ITL(NR)) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITL(NR)+1) - DELTHEX*0.5D0
!
         FACT = (1.D0+EPSRM2(NR)) / EPSRM2(NR)
         phib = ACOS(1.D0- FACT*COS(THEX)**2 )
         DELH=phib/NAVMAX
         suml=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            X=EPSRM2(NR)*COS(ETAL)*RR
            PSIB=(1.D0+EPSRM2(NR))/(1.D0+X/RR)
            PCOS=SQRT(1.D0-PSIB*SIN(THEX)**2)
            suml=suml+1.D0/PCOS
         END DO
         RINT0 = suml*DELH
         L_ex = RINT0*ABS(COS(THEX)) * ( QLM(NR)*RR ) * SIN(THEX)*DELTHEX
      ELSE ! ITL = PASSING
         DEL2T = -COSG(ITL(NR)+1)/mub +1.D0
         DELTHB = ( THG(ITL(NR)+1) - ACOS(mub) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITL(NR)) + DELTHEX*0.5D0
!
!         NTHX = ITL(NR)
!         NRX = NR
!         CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_TPB)
!         L_ex = RINT0 * ( QLM(NR)*RR ) * 2.D0 * SIN(THEX)*DELTHEX
         DELH=PI/NAVMAX
         suml=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            X=EPSRM2(NR)*COS(ETAL)*RR
            PSIB=(1.D0+EPSRM2(NR))/(1.D0+X/RR)
            PCOS=SQRT(1.D0-PSIB*SIN(THEX)**2)
            suml=suml+1.D0/PCOS
         END DO
         RINT0 = suml*DELH
         L_ex = RINT0*ABS(COS(THEX)) * ( QLM(NR)*RR ) * SIN(THEX)*DELTHEX
      END IF

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM_BL = 0.D0
      DO i = 0, m
         SUM_BL = SUM_BL + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
      K_1_x  = mub*sinb*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
      L_b = 8.D0*QLM(NR)*RR*DEL2T*( K_1_x - SUM_BL ) * 0.5D0
!     L_b means lambda_B * SIN(TEHTA_B) *DELTHB

!      IF(COSM(ITL(NR)).le.mub)THEN ! ITL = TRAPPED
         RLAMDA(ITL(NR),NR) = ( L_b + L_ex )/( SINM(ITL(NR))*DELTH )
         RLAMDA(ITU(NR),NR) = RLAMDA(ITL(NR),NR)
!      END IF

      END SUBROUTINE SET_RLAMDA_TPB3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB3_G(NR)

      USE libde, only: DEFT
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i, NG
      DOUBLE PRECISION:: DEL2T1, DEL2T2, DEL2T, mub,sinb, DELTHB, z, RINT0, ES, ram_bl
      DOUBLE PRECISION:: k1, k2, sum_bl, sum_ex, K_1_x, L_bl, L_ex, L_b, THEX, DELTHEX
      DOUBLE PRECISION,dimension(0:m):: alpha, beta, j2m
      DOUBLE PRECISION,DIMENSION(0:4):: ell_a, ell_b
      DOUBLE PRECISION:: phib, DELH, suml, ETAL, PSIB, X, PCOS, FACT

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRG2(NR) / (1.D0 + EPSRG2(NR) ) )
      sinb= SQRT( (1.D0-EPSRG2(NR))/(1.D0+EPSRG2(NR)) )

      IF(COSM(ITLG(NR)).le.mub)THEN ! ITL = TRAPPED
!!     define DEL2T
         DEL2T = COSG(ITLG(NR))/mub -1.D0
         DELTHB =( ACOS(mub) - THG(ITLG(NR)) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITLG(NR)+1) - DELTHEX*0.5D0
!
         FACT = (1.D0+EPSRG2(NR)) / EPSRG2(NR)
         phib = ACOS(1.D0- FACT*COS(THEX)**2 )
         DELH=phib/NAVMAX
         suml=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            X=EPSRG2(NR)*COS(ETAL)*RR
            PSIB=(1.D0+EPSRG2(NR))/(1.D0+X/RR)
            PCOS=SQRT(1.D0-PSIB*SIN(THEX)**2)
            suml=suml+1.D0/PCOS
         END DO
         RINT0 = suml*DELH
         L_ex = RINT0*ABS(COS(THEX)) * ( QLG(NR)*RR ) * SIN(THEX)*DELTHEX
      ELSE ! ITL = PASSING
         DEL2T = -COSG(ITLG(NR)+1)/mub +1.D0
         DELTHB = ( THG(ITLG(NR)+1) - ACOS(mub) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITLG(NR)) + DELTHEX*0.5D0
!
!         NTHX = ITLG(NR)
!         NRX = NR
!         CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPF_RLAMDA_DE_TPBG)
!         L_ex = RINT0 * ( QLG(NR)*RR ) * 2.D0 * SIN(THEX)*DELTHEX
         DELH=PI/NAVMAX
         suml=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            X=EPSRG2(NR)*COS(ETAL)*RR
            PSIB=(1.D0+EPSRG2(NR))/(1.D0+X/RR)
            PCOS=SQRT(1.D0-PSIB*SIN(THEX)**2)
            suml=suml+1.D0/PCOS
         END DO
         RINT0 = suml*DELH
         L_ex = RINT0*ABS(COS(THEX)) * ( QLG(NR)*RR ) * SIN(THEX)*DELTHEX
      END IF

!     obtain k1, k2, the coefficients of K(1-x)
      k1 = 0.D0
      k2 = 0.D0
      DO i = 0, 4
         k1 = k1 + ell_a(i) * (2.D0*DEL2T)**i
         k2 = k2 + ell_b(i) * (2.D0*DEL2T)**i
      END DO
!     calculate summention
      call RECURRENCE_ALPHA(m,alpha)
      call RECURRENCE_BETA(m,beta)
      SUM_BL = 0.D0
      DO i = 0, m
         SUM_BL = SUM_BL + alpha(i)*( mub )**(2*i+1)*beta(i)
      END DO
      K_1_x  = mub*sinb*(k1-k2*( LOG(2.D0*DEL2T)-1.D0 ) )
      L_b = 8.D0*QLG(NR)*RR*DEL2T*( K_1_x - SUM_BL ) * 0.5D0
!     L_b means lambda_B * SIN(TEHTA_B) *DELTHB

      IF(NR.ne.NRMAX+1)THEN
         RLAMDA_G(ITLG(NR),NR) = ( L_b + L_ex )/( SINM(ITLG(NR))*DELTH )
         RLAMDA_G(ITUG(NR),NR) = RLAMDA_G(ITLG(NR),NR)
      ELSE
         RLAMDA_GG(ITLG(NR),NR) = ( L_b + L_ex )/( SINM(ITLG(NR))*DELTH )
         RLAMDA_GG(ITUG(NR),NR) = RLAMDA_GG(ITLG(NR),NR)
      END IF

      END SUBROUTINE SET_RLAMDA_TPB3_G
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE FPBOUNCE
