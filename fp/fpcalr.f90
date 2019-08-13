MODULE fpcalr

! ****************************************
!     Radial transport
! ****************************************

CONTAINS

  SUBROUTINE FP_CALR

      USE fpcomm
      USE plprof
      USE cdbmfp_mod
      USE libmpi,ONLY: mtx_set_communicator, mtx_reset_communicator
      USE fpmpi,ONLY: p_theta_integration
      IMPLICIT NONE
      integer:: NSA, NSBA, NS, NR, NTH, NP, NG, NSB, NR1, NR2
      real(kind8):: RHON, RTFPL, FACTR, FACTP, SV, PPP, PPM
      real(kind8):: PSIB, PCOS, X, ETAL, sumd, sumf, DELH
      real(kind8):: DNDR, NEDGE, FACT, DINT_D, DINT_F, WRL
      real(kind8):: sum1, temp1, SRHODP, SRHODM, SRHOFP, SRHOFM
      real(kind8):: WRH, DFDR_D, DFDR_F, F_R2, DFDR_R2, F_R1, DFDR_R1
      REAL(kind8):: SHEAR,PNEL,RHONI,DPDR,DVEXBDR,CALF,CKAP,CEXB,FSZ,FEZ
      REAL(kind8),DIMENSION(1:NRMAX+1):: CHI_CDBM
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      double precision:: densm, densp, rgama
      INTEGER:: ISW_D

!---- Calculation of CDBM diffusion coefficient ----

      IF(MODELD_RDEP.EQ.2) THEN
         IF(nrank.EQ.0) THEN
            NR1=1
            NR2=NRMAX+1
!            write(18,'(A,1PE12.4/A)') 'T =',TIMEFP,&
!                     'NR,RS/QLM,SHEAR,PNEL,RHONI,DPDR,CHI_CDBM'
         ELSE
            NR1=NRSTART
            NR2=NRENDWG
         END IF
         DO NR=NR1,NR2
            RHON=RG(NR)

            IF(NR.EQ.1) THEN        ! magnetic shear s=(r/q)(dq/dr)
               SHEAR=0.D0
            ELSE IF(NR.EQ.NRMAX+1) THEN  !QLM(NRMAX+1)=QLG(rhon=1)
               SHEAR=(RHON/QLM(NR))*(QLM(NR)-QLM(NR-1)) &
                     *0.5D0/(RG(NR)-RM(NR-1))
            ELSE
               SHEAR=(RHON/(RM(NR)-RM(NR-1)))*(QLM(NR)-QLM(NR-1)) &
                     *2.D0/(QLM(NR)+QLM(NR-1))
            END IF

            PNEL=0.D0               ! electron density
            DO NS=1,NSMAX
               IF(ID_NS(NS).EQ.-1) THEN
                  IF(NR.EQ.1) THEN
                     PNEL=PNEL+RN_TEMP(NR,NS)*1.D20
                  ELSE IF(NR.EQ.NRMAX+1) THEN
                     PNEL=PNEL+RN_TEMP(NR-1,NS)*1.D20
                  ELSE
                     PNEL=PNEL+0.5D0*(RN_TEMP(NR-1,NS)+RN_TEMP(NR,NS))*1.D20
                  END IF
               END IF
            END DO

            RHONI=0.D0              ! ion mass density
            DO NS=1,NSMAX
               IF(ID_NS(NS).EQ.1) THEN
                  IF(NR.EQ.1) THEN
                     RHONI=RHONI+PA(NS)*AMP*RN_TEMP(NR,NS)*1.D20
                  ELSE IF(NR.EQ.NRMAX+1) THEN
                     RHONI=RHONI+PA(NS)*AMP*RN_TEMP(NR-1,NS)*1.D20
                  ELSE
                     RHONI=RHONI+PA(NS)*AMP &
                                *0.5D0*(RN_TEMP(NR-1,NS)+RN_TEMP(NR,NS))*1.D20
                  END IF
               END IF
            END DO

            DPDR=0.D0               ! pressure gradient
            PPP=0.D0
            PPM=0.D0
            DO NS=1,NSMAX
               IF(ID_NS(NS).EQ.1.OR.ID_NS(NS).EQ.-1) THEN
                  IF(NR.EQ.1) THEN
                     PPP=PPP+RHONI+RN_TEMP(NR,NS)*1.D20*RT_TEMP(NR,NS)*RKEV
                     PPM=PPM+RHONI+RN_TEMP(NR,NS)*1.D20*RT_TEMP(NR,NS)*RKEV
                  ELSE IF(NR.EQ.NRMAX) THEN
                     PPP=PPP+RHONI+RN_TEMP(NR  ,NS)*1.D20*RT_TEMP(NR  ,NS)*RKEV
                     PPM=PPM+RHONI+RN_TEMP(NR-1,NS)*1.D20*RT_TEMP(NR-1,NS)*RKEV
                  ELSE IF(NR.EQ.NRMAX+1) THEN
                     PPP=PPP+RHONI+RN_TEMP(NR-1,NS)*1.D20*RT_TEMP(NR-1,NS)*RKEV
                     PPM=PPM+RHONI+RN_TEMP(NR-2,NS)*1.D20*RT_TEMP(NR-2,NS)*RKEV
                  ELSE
                     PPP=PPP+RHONI+RN_TEMP(NR+1,NS)*1.D20*RT_TEMP(NR+1,NS)*RKEV
                     PPM=PPM+RHONI+RN_TEMP(NR-1,NS)*1.D20*RT_TEMP(NR-1,NS)*RKEV
                  END IF
               END IF
            END DO
            IF(NR.EQ.1) THEN
               DPDR=0.D0
            ELSE IF(NR.EQ.NRMAX) THEN
               DPDR=(PPP-PPM)/(RA*(RM(NR  )-RM(NR-1)))
            ELSE IF(NR.EQ.NRMAX+1) THEN
               DPDR=(PPP-PPM)/(RA*(RM(NR-1)-RM(NR-2)))
            ELSE
               DPDR=(PPP-PPM)/(RA*(RM(NR+1)-RM(NR-1)))
            END IF

            dvexbdr=0.D0
            CALF=1.D0
            CKAP=1.D0
            CEXB=0.D0
!            FSZ=1.D0     ! option
!            CURVZ=0.D0   ! option
!            FEZ=0.D0     ! option

            CALL CDBMFP(BB,RR,RA*RHON,RKAP,QLM(NR),SHEAR,PNEL,RHONI,DPDR, &
                      DVEXBDR,CALF,CKAP,CEXB,MODELD_CDBM,CHI_CDBM(NR))
!            IF(nrank.EQ.0) THEN
!               write(18,'(I2,1P7E11.3)') &
!                 NR,RA*RHON,QLM(NR),SHEAR,PNEL,RHONI,DPDR,CHI_CDBM(NR)
!            END IF
         END DO
      END IF

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         NSBA=NSB_NSA(NSA)

      DO NR=NRSTART,NRENDWG
         RHON=RG(NR)
         IF(NR.EQ.1) THEN
            NR1=NR
            NR2=NR
         ELSE IF(NR.EQ.NRMAX+1) THEN
            NR1=NR-1
            NR2=NR-1
         ELSE
            NR1=NR-1
            NR2=NR
         ENDIF

         IF(NR.EQ.1) THEN
            RTFPL= RTFP_G(NR,NSA)/RTFP0(NSA)
         ELSE
            RTFPL=RTFPS(NSA)/RTFP0(NSA)
         END IF

!------------- SET R DEPENDENCE
         SELECT CASE(MODELD_RDEP)
         CASE(0)
            FACTR=(DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
         CASE(1)
            FACTR=PI*RR*QLM(NR)*deltaB_B**2 *PTFP0(NS)/AMFP(NS)
         CASE(2)
            FACTR=FACTOR_CDBM*CHI_CDBM(NR)
         CASE DEFAULT
            WRITE(6,*) 'XX FPCALR: Undefined MODELD_RDEP'
            STOP
         END SELECT

!------------- SET EDGE VALUE
         SELECT CASE(MODELD_EDGE)
         CASE(1)
            IF(RHON.GT.RHO_EDGE) FACTR=DRR_EDGE
         CASE(2)
            IF(RHON.GT.RHO_EDGE) FACTR=FACTOR_DRR_EDGE*FACTR
         END SELECT
            
         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX

!------------- SET P DEPENDENCE
               SELECT CASE(MODELD_PDEP)
               CASE(0) ! no p dependence
                  FACTP=1.D0
               CASE(1) ! proportional to 1/sqrt{p_perp}
                  FACTP=(RTFPL/(RTFPL+PM(NP,NS)**2*SINM(NTH)**2))**0.25D0
               CASE(2) ! proportional to 1/p_perp
                  FACTP=SQRT(RTFPL/(RTFPL+PM(NP,NS)**2*SINM(NTH)**2))
               CASE(3) ! proportional to 1/p_perp^2
                  FACTP=RTFPL/(RTFPL+PM(NP,NS)**2*SINM(NTH)**2)
               CASE(4) ! stochastic delta B /B; relativistic
                  RGAMA=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                  FACTP=PM(NP,NS)*ABS(COSM(NTH))/RGAMA
               END SELECT
               DRR(NTH,NP,NR,NSA)= FACTR*FACTP/RA**2*RLAMDA_RG(NTH,NR)   ! normalization for rhon
            ENDDO
         ENDDO

! ------ SET PINCH TERM

         SELECT CASE(MODELD_PINCH)
         CASE(0) ! no pinch
            FACT=0.D0
         CASE(1) ! no radial particle transport (particle flux = 0)
            DINT_D=0.D0
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
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FS0(NTH,NP,NSA) ) &
                               / DELR *2.D0*0
                     F_R1 = FS0(NTH,NP,NSA)
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  ELSEIF(NR.eq.NRMAX+1)THEN ! FS2 = F at rho=1+delR/2
                     DFDR_R1 = ( FS2(NTH,NP,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                     IF(MODELD_boundary.eq.0)THEN
                        F_R1 = ( (1.D0-WRL)*FS2(NTH,NP,NSA) &
                             + WRL*FNSP(NTH,NP,NR-1,NSA) )
                     ELSEIF(MODELD_boundary.eq.1)THEN
                        F_R1 = FS1(NTH,NP,NSA)
                     END IF
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  ELSE
                     DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP,NR-1,NSA) ) &
                               / DELR
                     F_R1 = ( (1.D0-WRL)*FNSP(NTH,NP,NR,NSA) &
                            + WRL*FNSP(NTH,NP,NR-1,NSA) ) 
                     SRHODM=DFDR_R1 * DRR(NTH,NP,NR,NSA)
                     SRHOFM=F_R1    * DRR(NTH,NP,NR,NSA)
                  END IF
                  DINT_D = DINT_D + VOLP(NTH,NP,NS)*SRHODM
                  DINT_F = DINT_F + VOLP(NTH,NP,NS)*SRHOFM
               END DO
            END DO
! integration
            CALL mtx_set_communicator(comm_np) 
            CALL p_theta_integration(DINT_D) 
            CALL p_theta_integration(DINT_F) 
            CALL mtx_reset_communicator 

!         WRITE(*,'(A,2I3,2E14.6)') "DINT=", NSA,NR,DINT_D, DINT_F

            FACT=DINT_D/DINT_F

         CASE(2) ! pinch for Gaussian profile
            FACT=-2.D0*FACTOR_PINCH*RHON/RA
         END SELECT

         DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               FRR(NTH,NP,NR,NSA) = FACT * DRR(NTH,NP,NR,NSA)
            END DO
         END DO

      END DO ! NR
      END DO ! NSA

      RETURN
      END SUBROUTINE FP_CALR

    END MODULE fpcalr
