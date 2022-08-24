!     $Id: fpbounce.f90,v 1.11 2013/01/11 07:06:52 nuga Exp $

! *****************************
!     PREPARATION OF PARAMETERS FOR BOUNCE AVERAGE
! *****************************
      MODULE fpbounce

      USE fpcomm
      USE fpinit
      USE equnit

      INTEGER,parameter:: mmax=500

      PRIVATE
      PUBLIC set_bounce_param
      PUBLIC set_rfsad

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_BOUNCE_PARAM

      USE libmtx 
      IMPLICIT NONE
      INTEGER:: NR, NTH
      REAL(rkind):: FACT

!     ON RM(NR)
      DO NR=NRSTART, NREND
         CALL SET_ETAMG(NR) ! obtain the poloidal angle of bounce points
         CALL SET_RLAMDA(NR) ! NAVMAX 
!multiply A_chi0
         A_chi0(NR) = 2.D0*PI*RA**2*RM(NR)*(RR+RA*RM(NR))/(QLM(NR)*RR)
         FACT = RA**2*RM(NR)*( RR+RA*RM(NR) )*2.D0 ! *2.D0 from oint 
         DO NTH=1, NTHMAX
            RLAMDA(NTH,NR)=RLAMDA(NTH,NR)*FACT
         END DO
         CALL SET_RLAMDA_TPB_FROM_DENS(NR) ! RLAMDA(ITL,NR) is set to satisfy init dens
      END DO
! on RM(NRMAX)+DELR
      CALL SET_ETAMG_NRMAXP1
      CALL SET_RLAMDA_NRMAXP1
      FACT = RA**2*(RM(NRMAX)+DELR)*( RR+RA*(RM(NRMAX)+DELR) )*2.D0 ! *2.D0 from oint 
      DO NTH=1, NTHMAX
         RLAMDA_NRMAXP1(NTH)=RLAMDA_NRMAXP1(NTH)*FACT
      END DO
      CALL SET_RLAMDA_TPB_FROM_DENS_NRMAXP1

!     ON RG(NR)
      DO NR=NRSTART, NRENDWG
         IF(NR.eq.1)THEN
            DO NTH=1,NTHMAX
               RLAMDA_RG(NTH,NR)=0.D0
            END DO
         ELSE
            CALL SET_ETAMG_RG(NR)
            CALL SET_RLAMDA_RG(NR)
         END IF
         IF(NR.ne.1)THEN
            DO NTH=1, NTHMAX
               RLAMDA_RG(NTH,NR)=RLAMDA_RG(NTH,NR)*RA**2*RG(NR)*( RR+RA*RG(NR) )*2.D0
            END DO
         END IF
         IF(NR.ne.1) CALL SET_RLAMDA_RG_TPB_FROM_DENS(NR) ! RLAMDA(ITL,NR) is set to satisfy init dens
      END DO

!     RLAMDA means lambda*A_chi0
!     lambda can be obtained by RLAMDA/A_chi0
!     Flux surface averaged dens. can be obtained by RLAMDA*RFSAD
!     RFSAD means 2*pi/V'

      END SUBROUTINE SET_BOUNCE_PARAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET BOUNCE POINT chi_b/2=ETAM ON RM(NR)
!     ETAM: ON THM, ETAG: ON THG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      REAL(RKIND):: EPSL, FACT, A1, pitch

      EPSL=EPSRM2(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)

      pitch = ACOS(SQRT(2.D0*EPSRM2(NR)/(1.D0+EPSRM2(NR))))
      
      DO NTH=1,ITL(NR)-1
         ETAM(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITL(NR),ITU(NR)
         A1=FACT*COSM(NTH)**2
         ETAM(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITU(NR)+1,NTHMAX
         ETAM(NTH,NR)=PI*0.5D0
      ENDDO
!
      DO NTH=1,ITL(NR)-1
         ETAG(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITL(NR)+1,ITU(NR)
         A1=FACT*COSG(NTH)**2
         ETAG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITU(NR)+2,NTHMAX+1
         ETAG(NTH,NR)=PI*0.5D0
      ENDDO

      IF(THG(ITL(NR)).ge.pitch)THEN ! trapped
         A1=FACT*COSG(ITL(NR))**2
         ETAG(ITL(NR),NR)  =0.5D0*DACOS(1.D0-2.D0*A1)
         ETAG(ITU(NR)+1,NR)=ETAG(ITL(NR),NR)
      ELSE
         ETAG(ITL(NR),NR)=PI*0.5D0
         ETAG(ITU(NR)+1,NR)=PI*0.5D0
      END IF

      END SUBROUTINE SET_ETAMG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG_NRMAXP1

      IMPLICIT NONE
      INTEGER:: NTH
      REAL(RKIND):: EPSL, FACT, A1, pitch

      EPSL=EPSRM2(NRMAX+1)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)

      pitch = ACOS(SQRT(2.D0*EPSRM2(NRMAX+1)/(1.D0+EPSRM2(NRMAX+1))))
      
      DO NTH=1,ITL(NRMAX+1)-1
         ETAM_NRMAXP1(NTH)=PI*0.5D0
      ENDDO
      DO NTH=ITL(NRMAX+1),ITU(NRMAX+1)
         A1=FACT*COSM(NTH)**2
         ETAM_NRMAXP1(NTH)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITU(NRMAX+1)+1,NTHMAX
         ETAM_NRMAXP1(NTH)=PI*0.5D0
      ENDDO
!
      DO NTH=1,ITL(NRMAX+1)-1
         ETAG_NRMAXP1(NTH)=PI*0.5D0
      ENDDO
      DO NTH=ITL(NRMAX+1)+1,ITU(NRMAX+1)
         A1=FACT*COSG(NTH)**2
         ETAG_NRMAXP1(NTH)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITU(NRMAX+1)+2,NTHMAX+1
         ETAG_NRMAXP1(NTH)=PI*0.5D0
      ENDDO

      IF(THG(ITL(NRMAX+1)).ge.pitch)THEN ! trapped
         A1=FACT*COSG(ITL(NRMAX+1))**2
         ETAG_NRMAXP1(ITL(NRMAX+1))  =0.5D0*DACOS(1.D0-2.D0*A1)
         ETAG_NRMAXP1(ITU(NRMAX+1)+1)=ETAG_NRMAXP1(ITL(NRMAX+1))
      ELSE
         ETAG_NRMAXP1(ITL(NRMAX+1))=PI*0.5D0
         ETAG_NRMAXP1(ITU(NRMAX+1)+1)=PI*0.5D0
      END IF

      END SUBROUTINE SET_ETAMG_NRMAXP1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET BOUNCE POINT chi_b/2=ETAM ON RG(NR)
!     _G means ON RG GRID
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG_RG(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      REAL(RKIND):: EPSL, FACT, A1, pitch

      EPSL=EPSRG2(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      pitch = ACOS(SQRT(2.D0*EPSRG2(NR)/(1.D0+EPSRG2(NR))))
      
      DO NTH=1,ITLG(NR)-1
         ETAM_RG(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITLG(NR),ITUG(NR)
         A1=FACT*COSM(NTH)**2
         ETAM_RG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITUG(NR)+1,NTHMAX
         ETAM_RG(NTH,NR)=PI*0.5D0
      ENDDO
!
      DO NTH=1,ITLG(NR)-1
         ETAG_RG(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITLG(NR)+1,ITUG(NR)
         A1=FACT*COSG(NTH)**2
         ETAG_RG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITUG(NR)+2,NTHMAX+1
         ETAG_RG(NTH,NR)=PI*0.5D0
      ENDDO

      IF(THG(ITL(NR)).ge.pitch)THEN ! trapped
         A1=FACT*COSG(ITL(NR))**2
         ETAG_RG(ITLG(NR),NR)  =0.5D0*DACOS(1.D0-2.D0*A1)
         ETAG_RG(ITUG(NR)+1,NR)=ETAG_RG(ITL(NR),NR)
      ELSE
         ETAG_RG(ITLG(NR),NR)=PI*0.5D0
         ETAG_RG(ITUG(NR)+1,NR)=PI*0.5D0
      END IF

      END SUBROUTINE SET_ETAMG_RG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     BOUNDARY LAYER THM = THETA_B
!     P.129 and 100 Killeen 
      SUBROUTINE SET_RLAMDA_TPB(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: i
      REAL(RKIND):: DEL2T1, DEL2T2, DEL2T, mub,sinb
      REAL(RKIND):: k1, k2, sum, K_1_x
      REAL(RKIND),dimension(0:m):: alpha, beta
      REAL(RKIND),DIMENSION(0:4):: ell_a, ell_b

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
!      RLAMDA(ITL(NR),NR) = 4.D0*QLM(NR)*RR*(K_1_x - sum)*DEL2T/(sinb*DELTH)
!      RLAMDA(ITL(NR),NR) = 2.D0*(K_1_x - sum)*DEL2T/(sinb*DELTH)
      RLAMDA(ITL(NR),NR) = 4.D0*(K_1_x - sum)*DEL2T/(sinb*DELTH)
      RLAMDA(ITU(NR),NR) = RLAMDA(ITL(NR),NR)

      END SUBROUTINE SET_RLAMDA_TPB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB_G(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: i
      REAL(RKIND):: DEL2T1, DEL2T2, DEL2T, mub
      REAL(RKIND):: k1, k2, sum, K_1_x
      REAL(RKIND),dimension(0:m):: alpha, beta
      REAL(RKIND),DIMENSION(0:4):: ell_a, ell_b
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
      RLAMDA_RG(ITLG(NR),NR) = 2.D0*QLG(NR)*RR*(K_1_x - sum)/SINM(ITLG(NR))
      RLAMDA_RG(ITUG(NR),NR) = RLAMDA_RG(ITLG(NR),NR)

      END SUBROUTINE SET_RLAMDA_TPB_G

!-----------------------------
      SUBROUTINE RECURRENCE_ALPHA(m,alpha)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: m
      REAL(RKIND),dimension(0:mmax),INTENT(OUT):: alpha
      INTEGER:: i
      REAL(RKIND),parameter:: a0=1.D0, a1=-0.5D0
      
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
      REAL(RKIND),dimension(0:m),INTENT(OUT):: beta
      INTEGER:: i
      REAL(RKIND),parameter:: b0=0.D0, b1=1.D0
      
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
!    USE FOR TRAPPED REGION ONLY ON RM(NR)
     SUBROUTINE SET_RLAMDA_ELL(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      REAL(RKIND):: z, sum, mub, diff, pre_sum
      REAL(RKIND),dimension(0:m):: alpha, j2m

      mub=SQRT(2.D0*EPSRM2(NR)/(1.D0+EPSRM2(NR)))
         DO NTH = 1, NTHMAX/2
            IF(NTH.ne.ITL(NR))THEN
               z = COSM(NTH)/mub
               call RECURRENCE_ALPHA(m,alpha)
               call RECURRENCE_J(m,z,j2m)
               SUM = 0.D0
               i=0
               diff=1.D0
               DO WHILE(diff.ge.1.e-6)
                  PRE_SUM = SUM
                  SUM = SUM + alpha(i)*( mub )**(2*i)*j2m(i)
                  diff = ABS(SUM - PRE_SUM)
                  IF(NR.eq.1)THEN
                     IF(NTH.eq.ITL(NR)-1.or.NTH.eq.ITL(NR)-2)THEN
                        WRITE(*,'(A,2I5,10E14.6)') "TEST", NTH, i, alpha(i), j2m(i), ( mub )**(2*i), SUM, PRE_SUM, DIFF
                     END IF
                  END IF
                  i = i + 1
                  IF(i.eq.m) diff=0.D0
               END DO
!               RLAMDA(NTH,NR) = SUM * ( QLM(NR)*RR ) * 2.D0
!               RLAMDA(NTH,NR) = SUM * 4.DO*PI
               RLAMDA(NTH,NR) = SUM * 2.D0 ! not back and forth, one way of oint
               RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
            END IF
         END DO


      END SUBROUTINE SET_RLAMDA_ELL
!
!-----------------------------
!     ON RG(NR)
      SUBROUTINE SET_RLAMDA_ELL_G(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      REAL(RKIND):: z, sum, mub
      REAL(RKIND),dimension(0:m):: alpha, j2m

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
            RLAMDA_RG(NTH,NR) = SUM * ( QLG(NR)*RR ) * 2.D0
         END DO
!     BOUNDARY
!         RLAMDA_RG(ITLG(NR),NR) = 0.5D0* &
!              ( RLAMDA_RG(ITLG(NR)-1,NR) + RLAMDA_RG(ITLG(NR)+1,NR) )
!     symmetry
         DO NTH=ITLG(NR)+1,NTHMAX/2
            RLAMDA_RG(NTHMAX-NTH+1,NR)=RLAMDA_RG(NTH,NR)
         ENDDO
      END IF

      END SUBROUTINE SET_RLAMDA_ELL_G
!-----------------------------
      SUBROUTINE RECURRENCE_J(m,z,j2m)

      USE LIBELL,only : ellfc, ellec
      IMPLICIT NONE
      INTEGER,INTENT(IN):: m
      INTEGER:: IERR, i
      REAL(RKIND),INTENT(IN):: z
      REAL(RKIND),dimension(0:mmax),INTENT(OUT):: j2m
      REAL(RKIND):: J0, J2, ellK, ellE

      IF(z.ge.1.D0)THEN ! passing
!         WRITE(*,*) "ELL1", 1/z**2
         ellK = ELLFC(1/z**2,IERR)
         ellE = ELLEC(1/z**2,IERR)
!         WRITE(*,*) "ELL12", ellK, ellE
         J0 = ellK
         J2 = z**2*( ellK - ellE )
      ELSEIF(z.lt.1.D0.and.z.gt.0)THEN ! trapped
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
              - (2.D0*i - 3.D0)*z**2*j2m(i-2) ) &
              /(2.D0*i-1.D0)
      END DO

      END SUBROUTINE RECURRENCE_J
!
!-----------------------------
!!    CALCULATE FLUX SURFACE AVERAGE oint 1/psi ds
      SUBROUTINE SET_RFSAD(NR)

      IMPLICIT NONE
      INTEGER:: NR

      RFSADG(NR)=1.D0
      RFSADG_RG(NR)=1.D0 ! RFSAD GLOBAL on R_GRID

      IF(NR.ne.NRMAX+1)THEN
!         RFSADG(NR)=(1.D0+EPSRM2(NR))/(2.D0*PI*RA**2*RM(NR)*(RR+RA*RM(NR)) )
         RFSADG(NR)=2.D0*PI/( 2.D0*PI*RR *2.D0*PI*RA**2*RM(NR) )
         Line_Element(NR) = RR*QLM(NR)/(2.D0*PI)
      ELSE
         RFSADG(NR)=2.D0*PI/( 2.D0*PI*RR *2.D0*PI*RA**2*(RM(NRMAX)+DELR) )
         Line_Element(NR) = RR*QLM(NR)/(2.D0*PI)         
      END IF

      IF(NR.eq.1)THEN
         RFSADG_RG(NR)=(1.D0+EPSRM2(NR))/(2.D0*PI*RA**2)
      ELSE
         RFSADG_RG(NR)=(1.D0+EPSRG2(NR))/(2.D0*PI*RA**2*RG(NR)*(RR+RA*RG(NR)) )
      END IF

      END SUBROUTINE SET_RFSAD
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA(NR)

      USE libmtx
      USE plprof

      IMPLICIT NONE
      integer:: NR, NTH, NG
      REAL(rkind):: RINT0, DELH
      REAL(rkind):: SUML, ETAL, PSIB, PCOS

      DO NTH=1,NTHMAX/2
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            suml=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               PSIB=(1.D0+EPSRM2(NR))/(1.D0+EPSRM2(NR)*COS(ETAL))
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               suml=suml+1.D0/PCOS
            END DO
            RINT0 = suml*DELH
            RLAMDA(NTH,NR)=RINT0*ABS(COSM(NTH)) ! lambda = |cos\theta_0| oint 1/|cos\theta|
            RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
      ENDDO

      END SUBROUTINE SET_RLAMDA
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA_NRMAXP1

      USE libmtx
      USE plprof

      IMPLICIT NONE
      integer:: NTH, NG
      REAL(rkind):: RINT0, DELH
      REAL(rkind):: SUML, ETAL, PSIB, PCOS

      DO NTH=1,NTHMAX/2
            DELH=2.D0*ETAM_NRMAXP1(NTH)/NAVMAX
            suml=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               PSIB=(1.D0+EPSRM2(NRMAX+1))/(1.D0+EPSRM2(NRMAX+1)*COS(ETAL))
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               suml=suml+1.D0/PCOS
            END DO
            RINT0 = suml*DELH
            RLAMDA_NRMAXP1(NTH)=RINT0*ABS(COSM(NTH)) ! lambda = |cos\theta_0| oint 1/|cos\theta|
            RLAMDA_NRMAXP1(NTHMAX-NTH+1)=RLAMDA_NRMAXP1(NTH)
      ENDDO

      END SUBROUTINE SET_RLAMDA_NRMAXP1
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA_RG(NR)

      IMPLICIT NONE
      integer:: NR, NTH, NG
      REAL(rkind):: RINT0, DELH
      REAL(rkind):: SUML, ETAL, X, PSIB, PCOS

      DO NTH=1,NTHMAX/2
!         IF(NTH.ne.ITLG(NR)) THEN
            DELH=2.D0*ETAM_RG(NTH,NR)/NAVMAX
            suml=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRG2(NR)*COS(ETAL)
               PSIB=(1.D0+EPSRG2(NR))/(1.D0+X)
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               suml=suml+1.D0/PCOS
            END DO
            RINT0 = suml *DELH
            RLAMDA_RG(NTH,NR)=RINT0*ABS(COSM(NTH)) 
            RLAMDA_RG(NTHMAX-NTH+1,NR)=RLAMDA_RG(NTH,NR)
!         ENDIF
      ENDDO

      END SUBROUTINE SET_RLAMDA_RG
!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA_RGMAX(NR)

      IMPLICIT NONE
      integer:: NR, NTH, NG
      REAL(rkind):: RINT0, DELH
      REAL(rkind):: SUML, ETAL, X, PSIB, PCOS

      DO NTH=1,NTHMAX/2
            DELH=2.D0*ETAMG_RG(NTH,NR)/NAVMAX
            suml=0.D0 
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRG2(NR)*COS(ETAL)
               PSIB=(1.D0+EPSRG2(NR))/(1.D0+X)
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               suml=suml + 1.D0/PCOS
            END DO
            RINT0 = suml*DELH
            RLAMDAG_RG(NTH,NR)=RINT0*ABS(COSM(NTH))! lambda = v_0*|cos\theta_0|*tau_B
            RLAMDAG_RG(NTHMAX-NTH+1,NR)=RLAMDAG_RG(NTH,NR)
      ENDDO

      END SUBROUTINE SET_RLAMDA_RGMAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_ETAMG_RGMAX(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER:: NTH
      REAL(RKIND):: EPSL, FACT, A1

      EPSL=EPSRG2(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      
      DO NTH=1,ITLG(NR)-1
         ETAMG_RG(NTH,NR)=PI*0.5D0
      ENDDO
      DO NTH=ITLG(NR),ITUG(NR)
         A1=FACT*COSM(NTH)**2
         ETAMG_RG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      DO NTH=ITUG(NR)+1,NTHMAX
         ETAMG_RG(NTH,NR)=PI*0.5D0
      ENDDO

      END SUBROUTINE SET_ETAMG_RGMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB_GMAX(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: i
      REAL(RKIND):: DEL2T1, DEL2T2, DEL2T, mub
      REAL(RKIND):: k1, k2, sum, K_1_x
      REAL(RKIND),dimension(0:m):: alpha, beta
      REAL(RKIND),DIMENSION(0:4):: ell_a, ell_b
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
      RLAMDAG_RG(ITLG(NR),NR) = 2.D0*QLG(NR)*RR*(K_1_x - sum)/SINM(ITLG(NR))
      RLAMDAG_RG(ITUG(NR),NR) = RLAMDAG_RG(ITLG(NR),NR)

      END SUBROUTINE SET_RLAMDA_TPB_GMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------
!     ON RG(NR)
      SUBROUTINE SET_RLAMDA_ELL_GMAX(NR)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: NTH, i
      REAL(RKIND):: z, sum, mub
      REAL(RKIND),dimension(0:m):: alpha, j2m

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
            RLAMDAG_RG(NTH,NR) = SUM * ( QLG(NR)*RR ) * 2.D0
         END DO
!     symmetry
         DO NTH=ITLG(NR)+1,NTHMAX/2
            RLAMDAG_RG(NTHMAX-NTH+1,NR)=RLAMDAG_RG(NTH,NR)
         ENDDO
      END IF

      END SUBROUTINE SET_RLAMDA_ELL_GMAX
!---------------------------------
! ============================================================
      SUBROUTINE SET_RLAMDA_TPB4(NR)

      USE libde, only: DEFT
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: i
      REAL(RKIND):: DEL2T, mub,sinb, DELTHB, z
      REAL(RKIND):: k1, k2, sum_bl, sum_ex, K_1_x, L_ex, L_b, THEX, DELTHEX
      REAL(RKIND):: diff, PRE_SUM
      REAL(RKIND),dimension(0:m):: alpha, beta, j2m
      REAL(RKIND),DIMENSION(0:4):: ell_a, ell_b

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRM2(NR) / (1.D0 + EPSRM2(NR) ) )
      sinb= SQRT( (1.D0-EPSRM2(NR))/(1.D0+EPSRM2(NR)) )

      IF(COSM(ITL(NR)).le.mub)THEN ! THM(ITL) = TRAPPED
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
         i=0
         diff=1.D0
         DO WHILE(diff.ge.1.e-6)
            PRE_SUM=SUM_EX
            SUM_EX = SUM_EX + alpha(i)*( mub )**(2*i)*j2m(i)
            diff = ABS(SUM_EX - PRE_SUM) 
            i = i + 1
            IF(i.eq.m) diff=0.D0
         END DO
         L_ex = SUM_EX*( QLM(NR)*RR )*2.D0 * SIN(THEX)*DELTHEX
      ELSE ! ITL = PASSING
         DEL2T = -COSG(ITL(NR)+1)/mub +1.D0
         DELTHB = ( THG(ITL(NR)+1) - ACOS(mub) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITL(NR)) + DELTHEX*0.5D0
!
         z = COS(THEX)/mub
         call RECURRENCE_ALPHA(m,alpha)
         call RECURRENCE_J(m,z,j2m)
         SUM_EX = 0.D0
         i=0
         diff=1.D0
         DO WHILE(diff.ge.1.e-6)
            PRE_SUM=SUM_EX
            SUM_EX = SUM_EX + alpha(i)*( mub )**(2*i)*j2m(i)
            diff = ABS(SUM_EX - PRE_SUM) 
            i = i + 1
            IF(i.eq.m) diff=0.D0
         END DO
         L_ex = SUM_EX*( QLM(NR)*RR )*2.D0 * SIN(THEX)*DELTHEX
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
         RLAMDA(ITL(NR),NR) = ( L_b + L_ex )/( SINM(ITL(NR))*DELTH )/(QLM(NR)*RR)
         RLAMDA(ITU(NR),NR) = RLAMDA(ITL(NR),NR)
!      END IF

      END SUBROUTINE SET_RLAMDA_TPB4
!---------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB3(NR)

      USE libde, only: DEFT
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      INTEGER,parameter:: m=mmax
      INTEGER:: i, NG
      REAL(RKIND):: DEL2T1, DEL2T2, DEL2T, mub,sinb, DELTHB, RINT0
      REAL(RKIND):: k1, k2, sum_bl, K_1_x, L_ex, L_b, THEX, DELTHEX
      REAL(RKIND),dimension(0:m):: alpha, beta
      REAL(RKIND),DIMENSION(0:4):: ell_a, ell_b
      REAL(RKIND):: chib, DELH, suml, ETAL, PSIB, X, PCOS, FACT

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRM2(NR) / (1.D0 + EPSRM2(NR) ) )
      sinb= SQRT( (1.D0-EPSRM2(NR))/(1.D0+EPSRM2(NR)) )

      IF(COSM(ITL(NR)).lt.mub)THEN ! ITL = TRAPPED
!!     define DEL2T
         DEL2T = COSG(ITL(NR))/mub -1.D0
         DELTHB =( ACOS(mub) - THG(ITL(NR)) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITL(NR)+1) - DELTHEX*0.5D0
!
         FACT = (1.D0+EPSRM2(NR)) / EPSRM2(NR)
         chib = ACOS(1.D0- FACT*COS(THEX)**2 )
         DELH=chib/NAVMAX
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
      ELSEIF(COSM(ITL(NR)).gt.mub)THEN ! ITL = PASSING
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
      ELSE ! (COSM(ITL(NR)).eq.mub)THEN
         DEL2T1 = COSG(ITL(NR))/mub -1.D0
         DEL2T2 =-COSG(ITL(NR)+1)/mub +1.D0
         DEL2T = 0.5D0 *( DEL2T1 + DEL2T2 )

         L_ex=0.D0
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
      INTEGER:: i, NG
      REAL(RKIND):: DEL2T1, DEL2T2, DEL2T, mub,sinb, DELTHB, RINT0
      REAL(RKIND):: k1, k2, sum_bl, K_1_x, L_ex, L_b, THEX, DELTHEX
      REAL(RKIND),dimension(0:m):: alpha, beta
      REAL(RKIND),DIMENSION(0:4):: ell_a, ell_b
      REAL(RKIND):: chib, DELH, suml, ETAL, PSIB, X, PCOS, FACT

      DATA ell_a/ 1.38629436112D0, 0.09666344259D0, 0.03590092383D0, &
           0.03742563713D0, 0.01451196212D0/
      DATA ell_b/ 0.5D0, 0.12498593597D0, 0.06880248576D0, &
           0.03328355346D0, 0.00441787012D0/

      mub = SQRT(2.D0*EPSRG2(NR) / (1.D0 + EPSRG2(NR) ) )
      sinb= SQRT( (1.D0-EPSRG2(NR))/(1.D0+EPSRG2(NR)) )

      IF(COSM(ITLG(NR)).lt.mub)THEN ! ITL = TRAPPED
!!     define DEL2T
         DEL2T = COSG(ITLG(NR))/mub -1.D0
         DELTHB =( ACOS(mub) - THG(ITLG(NR)) )*2.D0
!     obtain L_ex
         DELTHEX = DELTH - DELTHB
         THEX = THG(ITLG(NR)+1) - DELTHEX*0.5D0
!
         FACT = (1.D0+EPSRG2(NR)) / EPSRG2(NR)
         chib = ACOS(1.D0- FACT*COS(THEX)**2 )
         DELH=chib/NAVMAX
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
      ELSEIF(COSM(ITL(NR)).gt.mub)THEN ! ITL = PASSING
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
      ELSE ! (COSM(ITL(NR)).eq.mub)THEN
         DEL2T1 = COSG(ITL(NR))/mub -1.D0
         DEL2T2 =-COSG(ITL(NR)+1)/mub +1.D0
         DEL2T = 0.5D0 *( DEL2T1 + DEL2T2 )

         L_ex=0.D0
      END IF

!      IF(NR.eq.NRMAX+1) WRITE(*,'(A,3E14.6,I4)') "BOUNCE", DEL2T, COSG(ITLG(NR)), mub, ITLG(NR)
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
         RLAMDA_RG(ITLG(NR),NR) = ( L_b + L_ex )/( SINM(ITLG(NR))*DELTH )
         RLAMDA_RG(ITUG(NR),NR) = RLAMDA_RG(ITLG(NR),NR)
      ELSE
         RLAMDAG_RG(ITLG(NR),NR) = ( L_b + L_ex )/( SINM(ITLG(NR))*DELTH )
         RLAMDAG_RG(ITUG(NR),NR) = RLAMDAG_RG(ITLG(NR),NR)
      END IF

      END SUBROUTINE SET_RLAMDA_TPB3_G
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB_FROM_DENS(NR)

      IMPLICIT NONE
      INTEGER:: NTH, ITLB, NS
      INTEGER,intent(in):: NR
      REAL(rkind):: RSUM1, RSUM2, RSUM3

!      NS=NS_NSA(NSASTART)
      NS=1
      RSUM1=0.D0
      RSUM2=0.D0
      RSUM3=0.D0
      ITLB=ITL_judge(NR)
!      ITLB=ITL(NR)-1

      DO NTH=1, NTHMAX/2
         RSUM1 = RSUM1 + VOLP(NTH,NPSTART,NS)
         IF(NTH.ne.ITLB)THEN
            RSUM2 = RSUM2 + VOLP(NTH,NPSTART,NS)*RLAMDA(NTH,NR)*RFSADG(NR)
         ELSE
            RSUM3 = VOLP(NTH,NPSTART,NS)*RFSADG(NR)
         END IF
      END DO

      RATIO_NAVMAX(NR) = ( (RSUM1 - RSUM2)/RSUM3 ) / RLAMDA(ITLB,NR)
!      RATIO_NAVMAX(NR) = ( (RSUM1 - RSUM2)/RSUM3 ) / RLAMDA(ITL_judge(NR),NR)
      
      RLAMDA(ITLB,NR) = (RSUM1 - RSUM2)/RSUM3
      RLAMDA(NTHMAX-ITLB+1,NR)=RLAMDA(ITLB,NR)
!      RLAMDA(ITL(NR)-1,NR) = (RSUM1 - RSUM2)/RSUM3
!      RLAMDA(NTHMAX-ITL(NR)+2,NR)=RLAMDA(ITL(NR)-1,NR)

      END SUBROUTINE SET_RLAMDA_TPB_FROM_DENS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_TPB_FROM_DENS_NRMAXP1

      IMPLICIT NONE
      INTEGER:: NTH, ITLB, NS
      REAL(rkind):: RSUM1, RSUM2, RSUM3
      
!      NS=NS_NSA(NSASTART)
      NS=1
      RSUM1=0.D0
      RSUM2=0.D0
      RSUM3=0.D0
      ITLB=ITL_judge(NRMAX+1)

      DO NTH=1, NTHMAX/2
         RSUM1 = RSUM1 + VOLP(NTH,NPSTART,NS)
         IF(NTH.ne.ITLB)THEN
            RSUM2 = RSUM2 + VOLP(NTH,NPSTART,NS)*RLAMDA_NRMAXP1(NTH)*RFSADG(NRMAX+1)
         ELSE
            RSUM3 = VOLP(NTH,NPSTART,NS)*RFSADG(NRMAX+1)
         END IF
      END DO
      
      RLAMDA_NRMAXP1(ITLB) = (RSUM1 - RSUM2)/RSUM3
      RLAMDA_NRMAXP1(NTHMAX-ITLB+1)=RLAMDA_NRMAXP1(ITLB)

      END SUBROUTINE SET_RLAMDA_TPB_FROM_DENS_NRMAXP1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SET_RLAMDA_RG_TPB_FROM_DENS(NR)

      IMPLICIT NONE
      INTEGER:: NTH, ITLB, NS
      INTEGER,intent(in):: NR
      REAL(rkind):: RSUM1, RSUM2, RSUM3

!      NS=NS_NSA(NSASTART)
      NS=1
      RSUM1=0.D0
      RSUM2=0.D0
      RSUM3=0.D0
      ITLB=ITLG_judge(NR)
!      ITLB=ITLG(NR)-1

      DO NTH=1, NTHMAX/2
         RSUM1 = RSUM1 +  VOLP(NTH,NPSTART,NS)
         IF(NTH.ne.ITLB)THEN
            RSUM2 = RSUM2 + VOLP(NTH,NPSTART,NS)*RLAMDA_RG(NTH,NR)*RFSADG_RG(NR)
         ELSE
            RSUM3 = VOLP(NTH,NPSTART,NS)*RFSADG_RG(NR)
         END IF
      END DO
      
!      RLAMDA_RG(ITLG_judge(NR),NR) = (RSUM1 - RSUM2)/RSUM3
!      RLAMDA_RG(NTHMAX-ITLG_judge(NR)+1,NR)=RLAMDA_RG(ITLG_judge(NR),NR)
      RLAMDA_RG(ITLB,NR) = (RSUM1 - RSUM2)/RSUM3
      RLAMDA_RG(NTHMAX-ITLB+1,NR)=RLAMDA_RG(ITLB,NR)

      END SUBROUTINE SET_RLAMDA_RG_TPB_FROM_DENS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE FPBOUNCE
