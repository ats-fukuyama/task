MODULE fpcalr

! ****************************************
!     Radial transport
! ****************************************

CONTAINS

  SUBROUTINE FP_CALR(NSA)

      USE fpcomm
      USE plprof
      USE libmpi,ONLY: mtx_set_communicator, mtx_reset_communicator
      USE fpmpi,ONLY: p_theta_integration
      IMPLICIT NONE
      integer:: NSA, NSBA, NS, NR, NTH, NP, NG
      real(kind8):: RHON, RTFPL, FACTR, FACTP, SV
      real(kind8):: PSIB, PCOS, X, ETAL, sumd, sumf, DELH
      real(kind8):: DNDR, NEDGE, FACT, DINT_D, DINT_F, WRL
      real(kind8):: sum1, temp1, SRHODP, SRHODM, SRHOFP, SRHOFM
      real(kind8):: WRH, DFDR_D, DFDR_F, F_R2, DFDR_R2, F_R1, DFDR_R1
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      double precision:: densm, densp, rgama
      INTEGER:: ISW_D

      ISW_D=MOD(MODELD,10)
      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

      DO NR=NRSTART,NRENDWG
         RHON=RG(NR)
         IF(NR.ne.NRMAX+1)THEN
            RTFPL= RTFP_G(NR,NSA)/RTFP0(NSA)
         ELSE
            RTFPL=RTFPS(NSA)/RTFP0(NSA)
         END IF
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
!------------- SET P DEPENDENCE
               SELECT CASE(ISW_D)
               CASE(0) ! no transport
                  FACTP=0.D0
                  FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
               CASE(1) ! no p dependence
                  FACTP=1.D0
                  FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
               CASE(2) ! depend on 1/p_perp
                  FACTP=SQRT(RTFPL)/SQRT(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
!                  FACTP=1.D0/SQRT(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
                  FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
               CASE(3) ! depend on 1/sqrt{p_perp}
                  FACTP=SQRT(RTFPL)/SQRT(RTFPL+PM(NP,NSBA)*SINM(NTH))
!                  FACTP=1.D0/SQRT(RTFPL+PM(NP,NSBA)*SINM(NTH))
                  FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
               CASE(4) ! depend on 1/p^2
                  FACTP=RTFPL/(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
!                  FACTP=1.D0/(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
                  FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
               CASE(5) ! stochastic delta B /B
                  RGAMA=SQRT(1.D0+THETA0(NSBA)*PM(NP,NSBA)**2)
                  FACTP=PM(NP,NSBA)*ABS(COSM(NTH))/RGAMA
                  FACTR=PI*RR*QLM(NR)*deltaB_B**2 *PTFP0(NSBA)/AMFP(NSBA)
               CASE(6) ! depend on 1/p and H-mode like (constant(1.D0) when rho >= 0.9)
                  FACTP=SQRT(RTFPL)/SQRT(RTFPL+PM(NP,NSBA)**2*SINM(NTH)**2)
                  IF(RHON.GE.0.9D0) THEN
                     FACTR=1.D0
                  ELSE
                     FACTR= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
                  END IF
               END SELECT
                  DRR(NTH,NP,NR,NSA)= FACTR*FACTP/(RA**2)
            ENDDO
         ENDDO

! ------ SET PINCH TERM
         DINT_D=0.D0
         DINT_F=1.D0
         IF(MODELD.ge.10)THEN
            DINT_F=0.D0
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  IF(TIMEFP.eq.0.and.N_IMPL.eq.0)THEN
                     IF(NR.eq.1)THEN
                        WRL=0.25D0 ! not necessary
                     ELSE
                        WRL=(4.D0*RG(NR)+DELR)/(8.D0*RG(NR))                     
                     END IF
                  ELSE
                     WRL=WEIGHR(NTH,NP,NR,NSA)
                  END IF

                  IF(NR.eq.1)THEN
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FS0(NTH,NP,NSA) ) / DELR *2.D0*0
                     F_R1 = FS0(NTH,NP,NSA)
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  ELSEIF(NR.eq.NRMAX+1)THEN ! FS2 = F at rho=1+delR/2
                     DFDR_R1 = ( FS2(NTH,NP,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     IF(MODELD_boundary.eq.0)THEN
                        F_R1 = ( (1.D0-WRL)*FS2(NTH,NP,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) )
                     ELSEIF(MODELD_boundary.eq.1)THEN
                        F_R1 = FS1(NTH,NP,NSA)
                     END IF
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  ELSE
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     F_R1 = ( (1.D0-WRL)*FNSP(NTH,NP,NR,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) ) 
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  END IF
                  DINT_D = DINT_D + VOLP(NTH,NP,NSBA)*SRHODM
                  DINT_F = DINT_F + VOLP(NTH,NP,NSBA)*SRHOFM
               END DO
            END DO
! integration
            CALL mtx_set_communicator(comm_np) 
            CALL p_theta_integration(DINT_D) 
            CALL p_theta_integration(DINT_F) 
            CALL mtx_reset_communicator 
!         WRITE(*,'(A,2I3,2E14.6)') "DINT=", NSA,NR,DINT_D, DINT_F
         ENDIF

         FACT=DINT_D/DINT_F
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               FRR(NTH,NP,NR,NSA) = FACT * DRR(NTH,NP,NR,NSA)
            END DO
         END DO

!         IF(NR.eq.NREND+1.and.NSA.eq.1) WRITE(*,'(I4,3E16.8)') NR, DRR(1,1,NR,1), FRR(1,1,NR,1),WRL
!         IF(NPSTART.eq.1.and.NSA.eq.1) WRITE(*,'(I4,3E16.8)') NR, DRR(1,1,NR,1), FRR(1,1,NR,1), WRL

! boudary condition
!         IF(NR.eq.NRMAX+1)THEN
!            DO NP=NPSTART, NPEND
!               DO NTH=1,NTHMAX
!                  DRR(NTH,NP,NR,NSA)= 0.D0
!                  FRR(NTH,NP,NR,NSA)= 0.D0
!               END DO
!            END DO
!         END IF

      ENDDO ! NR

      RETURN
      END SUBROUTINE FP_CALR

    END MODULE fpcalr
