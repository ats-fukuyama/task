!
! *************************
!     SAVE DATA ROUTINE
! *************************
!
      MODULE fpcaleind

      USE fpcomm
      USE libmpi
      USE plprof
      USE equnit_mod

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FP_CALE_IND(NSA)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      INTEGER:: NR, NTH, NP
      double precision:: rv, E_IND

      CALL Ip_r
      CALL UPDATE_PSIP_P ! poloidal flux at present step

      DO NR=NRSTART,NREND
         rv = EPSRM2(NR)*RR
         DO NP=1,NPMAX+1
            DO NTH=1,NTHMAX
               E_IND=-( PSIPM_P(NTH,NR)-PSIPM_M(NTH,NR) )/ (2.D0*PI*(RR+rv)*DELT)
               FEPP_IND(NTH,NP,NR,NSA) = AEFP(NSA)*E_IND/PTFP0(NSA)*COSM(NTH)
            END DO
         END DO
     END DO

      DO NR=NRSTART,NREND
         rv = EPSRM2(NR)*RR
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX+1
               E_IND=-( PSIPG_P(NTH,NR)-PSIPG_M(NTH,NR) )/ (2.D0*PI*(RR+rv)*DELT)
               FETH_IND(NTH,NP,NR,NSA) = -AEFP(NSA)*E_IND/PTFP0(NSA)*SING(NTH)
            END DO
         END DO
      END DO

      IF(MODELA.ne.0)THEN
         DO NR=NRSTART,NREND
            DO NP=1,NPMAX+1
               DO NTH=ITL(NR)+1,ITU(NR)-1
                  FEPP_IND(NTH,NP,NR,NSA)=0.D0
               END DO
               FEPP_IND(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
                    *( FEPP_IND(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +FEPP_IND(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) & 
                    +FEPP_IND(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                    +FEPP_IND(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR))
               FEPP_IND(ITU(NR),NP,NR,NSA)=FEPP_IND(ITL(NR),NP,NR,NSA) 
            END DO
!
            DO NP=1,NPMAX
               DO NTH=ITL(NR)+1,ITU(NR)
                  FETH_IND(NTH,NP,NR,NSA)=0.D0
               END DO
            END DO
!
         END DO
      END IF

      END SUBROUTINE FP_CALE_IND
!------------------------------------
      SUBROUTINE UPDATE_PSIP_M

      IMPLICIT NONE
      integer:: NR,NSA,NTH

      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            PSIPM_M(NTH,NR)=PSIPM_P(NTH,NR)
         END DO
         DO NTH=1,NTHMAX+1
            PSIPG_M(NTH,NR)=PSIPG_P(NTH,NR)
         END DO
      END DO

      END SUBROUTINE UPDATE_PSIP_M
!------------------------------------
      SUBROUTINE UPDATE_PSIP_P

      IMPLICIT NONE
      INTEGER:: NR, NTH
      double precision:: PSIP_P, E_IND

      DO NR=NRSTART, NREND
         IF(MODELA.eq.0)THEN
            CALL POLOIDAL_FLUX_IP(NR,PSIP_P)
            DO NTH=1,NTHMAX
               PSIPM_P(NTH,NR)=PSIP_P
            END DO
         ELSE
            DO NTH=1,NTHMAX
               CALL BOUNCE_AV_PSIP(0,NR,NTH,PSIP_P)
               PSIPM_P(NTH,NR)=PSIP_P
            END DO
         END IF
         E_IND=-( PSIPM_P(1,NR)-PSIPM_M(1,NR) )/ (2.D0*PI*RR*(1.D0+EPSRM2(NR))*DELT)
!         IF(NRANK.eq.0) WRITE(*,'(2I4,4E16.8)') N_IMPL, NR, E_IND, &
!              PSIPM_P(1,NR), PSIPM_M(1,NR), RIPP(1,1) 
      END DO
!
      DO NR=NRSTART, NREND
         IF(MODELA.eq.0)THEN
            CALL POLOIDAL_FLUX_IP(NR,PSIP_P)
            DO NTH=1,NTHMAX+1
               PSIPG_P(NTH,NR)=PSIP_P
            END DO
         ELSE
            DO NTH=1,NTHMAX+1
               CALL BOUNCE_AV_PSIP(1,NR,NTH,PSIP_P)
               PSIPG_P(NTH,NR)=PSIP_P
            END DO
         END IF
      END DO

      END SUBROUTINE UPDATE_PSIP_P
!------------------------------------
      SUBROUTINE POLOIDAL_FLUX_IP(NR,PSIP_P)
        
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      double precision,intent(out):: PSIP_P
      INTEGER:: NSA, NR2
      double precision:: rv, SUM, DELrho
      double precision,dimension(NRMAX):: RIP
      double precision:: PSIP0, PSIP1

      DO NR2=1,NRMAX
         RIP(NR2)=0.D0
      END DO
      DO NSA=1,NSAMAX
         DO NR2=1,NRMAX
            RIP(NR2)=RIP(NR2)+RIPP(NR2,NSA)
         END DO
      END DO

      SUM=0.D0
      DO NR2=1,NRMAX 
         rv = EPSRM2(NR2)*RR
         DELrho = (EPSRM2(NR2+1)-EPSRM2(NR2) ) * RR
         SUM = SUM + RMU0*RIP(NR2)*(RR/rv - 1.D0) * DELrho
      END DO
      rv = EPSRM2(NRMAX)*RR
      PSIP0 = SUM + &
           RMU0*RIP(NRMAX)*( RR*LOG(RR/rv) - (RR-rv) )

      SUM=0.D0
      DO NR2=1,NR
         rv = EPSRM2(NR2)*RR
         DELrho = (EPSRM2(NR2+1)-EPSRM2(NR2) ) * RR
         SUM = SUM + RMU0*RIP(NR2)*(RR/rv + 1.D0) * DELrho
      END DO
      PSIP1 = SUM
!      PSIP0=0.D0

      PSIP_P = PSIP0 - PSIP1

      END SUBROUTINE POLOIDAL_FLUX_IP
!------------------------------------
      SUBROUTINE BOUNCE_AV_PSIP(ISW,NR,NTH,PSIP_P)

      IMPLICIT NONE
      integer,intent(IN):: ISW, NR
      double precision,intent(out):: PSIP_P
      integer:: NSA, NP, NTH, NG
      double precision:: SUM, DELH, ETAL, PSI_PHI

      IF(ISW.eq.0)THEN
         DELH=2.D0*ETAM(NTH,NR)/NAVMAX 
      ELSEIF(ISW.eq.1)THEN
         DELH=2.D0*ETAG(NTH,NR)/NAVMAX 
      END IF

      SUM=0.D0
      DO NG=1,NAVMAX
         ETAL=DELH*(NG-0.5D0)
         CALL POLOIDAL_FLUX_IP_PHI(NR,ETAL,PSI_PHI)
         SUM= SUM + PSI_PHI
      END DO
      PSIP_P = SUM*DELH 

      END SUBROUTINE BOUNCE_AV_PSIP
!------------------------------------
      SUBROUTINE POLOIDAL_FLUX_IP_PHI(NR,ETAL,PSI_PHI)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR
      double precision,intent(in)::ETAL
      double precision,intent(out):: PSI_PHI
      INTEGER:: NSA, nx, NR2
      double precision:: SUM, PSIP0, PSIP1
      double precision,dimension(NRMAX):: RIP
      double precision:: zz,rv,xmax,xx,DELxx
      integer,parameter:: NXMAX=100
      double precision:: XIP

      RIP(:)=0.D0
      DO NSA=1,NSAMAX
         DO NR2=1,NRMAX
            RIP(NR2)=RIP(NR2)+RIPP(NR2,NSA)
         END DO
      END DO
      IF(ETAL.le.PI*0.5D0)THEN
         xmax=0.D0
      ELSE
         xmax=EPSRM2(NR)*RR*ABS(COS(ETAL))
      END IF
      zz = EPSRM2(NR)*RR*SIN(ETAL)

      DELxx = (RR-xmax)/NXMAX
      SUM=0.D0
      DO nx=1,nxmax
         xx = xmax + DELxx * (nx-0.5D0)
         rv = SQRT(xx**2 + zz*2)
         CALL IP_INTERPOLATION(rv,RIP,XIP)
         SUM = SUM + RMU0*XIP*xx/(rv**2)*(RR-xx)*DELxx
      END DO
      PSIP0=SUM

      SUM=0.D0
      IF(ETAL.le.PI*0.5D0)THEN
         xmax=EPSRM2(NR)*RR*ABS(COS(ETAL))
         DELxx = xmax/NXMAX
         DO nx=1,NXMAX
            xx = xmax + DELxx * (nx-0.5D0)
            rv = SQRT(xx**2 + zz*2)
            CALL IP_INTERPOLATION(rv,RIP,XIP)
            SUM = SUM + RMU0*XIP*xx/(rv**2)*(RR+xx)*DELxx
         END DO
         PSIP1=SUM
      ELSE
         PSIP1=0.D0
      END IF
      
      PSI_PHI = PSIP0 - PSIP1

      END SUBROUTINE POLOIDAL_FLUX_IP_PHI
!------------------------------------
      SUBROUTINE IP_INTERPOLATION(rv,RIP,XIP)

      IMPLICIT NONE
      double precision,INTENT(IN):: rv
      double precision,dimension(NRMAX),INTENT(IN):: RIP
      double precision,INTENT(out):: XIP
      integer:: NR, NRL, NRU
      double precision:: rv0, rv1, rv2

      NRL=0
      NRU=NRMAX
      rv0 = EPSRM2(NRMAX)*RR
      IF(rv0.le.rv)THEN
         XIP=RIP(NRMAX)
      ELSE
         DO NR=1,NRMAX
            rv0 = EPSRM2(NR)*RR
            IF(rv0.le.rv)THEN
               NRL=NR
            END IF
         END DO
         DO NR=NRMAX,1,-1
            rv0 = EPSRM2(NR)*RR
            IF(rv.le.rv0)THEN
               NRU=NR
            END IF
         END DO

!         WRITE(*,'(A,2I4,3E16.8)') "NRL,NRU ",NRL,NRU,EPSRM2(NRL)*RR, rv, EPSRM2(NRU)*RR
         IF(NRL.eq.0)THEN
            rv1=0.D0
            rv2=EPSRM2(NRU)*RR
            XIP=( RIP(NRU) )/( rv2-rv1 )*( rv-rv1 )
         ELSE
            rv1=EPSRM2(NRL)*RR
            rv2=EPSRM2(NRU)*RR
            XIP=RIP(NRL)+( RIP(NRU)-RIP(NRL) )/( rv2-rv1 )*( rv-rv1 )
         END IF
      END IF

      END SUBROUTINE IP_INTERPOLATION
!------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE Ip_r

      USE libmtx
      IMPLICIT NONE
      double precision,dimension(NRMAX,NSAMAX)::RJ_G

      CALL FPCURRENT(RJ_G)
      IF(NRANK.eq.0)THEN
         CALL FPRIPP(RJ_G)
      END IF
      CALL mtx_broadcast_real8(RIPP,NSAMAX*NRMAX)

      END SUBROUTINE Ip_r
!----------------------------------
      SUBROUTINE FPCURRENT(RJ_G)
!
      USE fpmpi
      USE libmtx
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NSW, N
      integer:: IERR
      real(8):: RSUM2, FACT, PV
      double precision,dimension(NRSTART:NREND,NSAMAX)::RJ_L
      double precision,dimension(NRMAX,NSAMAX),INTENT(OUT)::RJ_G

      CALL mtx_set_communicator(comm_np) 
      DO NR=NRSTART,NRENDX
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)

            RSUM2=0.D0
            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
!                  DO NP=1,NPMAX
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)
                     END DO
                  ENDDO
               ELSE
!                  DO NP=1,NPMAX
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)/PV
                     END DO
                  END DO
               ENDIF
            ELSE
               IF(MODELR.EQ.0) THEN
!                  DO NP=1,NPMAX
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)*RLAMDA(NTH,NR)
                     END DO
                  ENDDO
               ELSE
!                  DO NP=1,NPMAX
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)/PV*RLAMDA(NTH,NR)
                     END DO
                  END DO
               ENDIF               
            END IF

            CALL p_theta_integration(RSUM2) 

            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)*RCOEFNG(NR)
            RJ_L(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)!*1.D-6

         ENDDO ! NSA
      ENDDO ! NR

      CALL mtx_set_communicator(comm_nsanr)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RJ_L,SAVLEN(NRANK+1),RJ_G,N,NSA)
      END DO
      CALL mtx_reset_communicator 

      RETURN
      END SUBROUTINE FPCURRENT
!----------------------------------

      SUBROUTINE FPRIPP(RJ_G)
!
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH 
      real(8):: EAVE, EAVE2, rtemp, rtemp2, THETAL, THETAL2
      double precision,dimension(NRMAX,NSAMAX),INTENT(IN)::RJ_G

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RIPP(NR,NSA)=0.D0
         END DO
      ENDDO

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            IF(NR.eq.1)THEN ! current with in 0<rho<RM(NR)
               RIPP(NR,NSA)  =RJ_G(NR,NSA)*VOLR(NR)/(2.D0*PI*RR)
            ELSE
               RIPP(NR,NSA)  =RIPP(NR-1,NSA)+RJ_G(NR,NSA)*VOLR(NR)/(2.D0*PI*RR)
            END IF
         ENDDO
      ENDDO
!      IF(NRANK.eq.0) WRITE(*,*) "B_POL", RIPP(NRMAX,1)*RMU0/(2.D0*PI*RA)
         
      RETURN
      END SUBROUTINE FPRIPP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

    end MODULE fpcaleind




