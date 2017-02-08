!     $Id: fpsave.f90,v 1.41 2013/02/08 07:36:24 nuga Exp $
!
! *************************
!     SAVE DATA ROUTINE
! *************************
!
      MODULE fpsave

      USE fpcomm
      USE plprof
      USE fpexec
      USE libbes,ONLY: besekn,besjnv

      contains

!----------------------------------

      SUBROUTINE FPSSUB
!
      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NPS, NPE
      integer:: IERR
      real(8):: RSUM1, RSUM2, RSUM3, RSUM4, RSUM5, RSUM6, RSUM7, RSUM_FS2
      real(8):: RSUM8, RSUM9, RSUM123, RSUM11B,RSUM11F,RSUM11S,RSUM11L
      real(8):: RSUM12, RSUM_IC, rsum_test, RSUM1R,RSUM2R, RSUM_synch, RSUM_loss
      real(8):: PV, WPL, WPM, WPP
      real(8):: DFP, DFT, FFP, testa, testb, FACT, DFDP, WRL, WRH, DINT_DFDT_R1, DINT_DFDT_R2
      real(8):: DFDR_R1, DFDR_R2, DFDT_R1, DFDT_R2, DINT_DR, RSUM_DR,RGAMA,F_R1,F_R2, RSUMN_DR
      real(8),dimension(NSBMAX):: RSUM10
      real(8),dimension(NPMAX):: RSUMNP_E
      real(8),dimension(NTHMAX,NPMAX):: T_BULK
      real(8),dimension(NPSTART:NPEND):: RPL_BULK_send
      real(8),dimension(NPMAX):: RPL_BULK_recv
!      integer,dimension(NRSTART:NREND,NSBMAX):: NP_BULK
      real(8):: ratio, RSUM_T, RSUM_V, P_BULK_R, RATIO_0_S
      real(8):: SRHOR1, SRHOR2, RSUM_DRS, RSUMN_DRS
      real(8):: DSDR, SRHOP, SRHOM, RSUM63
      real:: gut1,gut2
      real(8)::FLUXS_PMAX,FLUXS_PC, rtemp, DRRM, DRRp, RL, PMAX_BULK

      IF(ISAVE.NE.0) RETURN
      CALL GUTIME(gut1)

!      DO NSA=NSASTART,NSAEND
      DO NSA=1,NSAMAX
         RNDRS(NSA) =0.D0
         RPDRS(NSA) =0.D0
         DO NR=NRSTART, NREND
            RPDRL(NR,NSA)=0.D0
            RNDRL(NR,NSA)=0.D0
         END DO
      END DO

      CALL mtx_set_communicator(comm_np) 
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)

            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            RSUM5=0.D0
            RSUM6=0.D0
            RSUM63=0.D0
            RSUM7=0.D0
            RSUM8=0.D0
            RSUM9=0.D0
            RSUM_IC=0.D0
            RSUM_synch=0.D0
            RSUM_loss=0.D0
            DO NSB=1,NSBMAX
               RSUM10(NSB)=0.D0
            ENDDO
            RSUM11B=0.D0
            RSUM11F=0.D0
            RSUM11S=0.D0
            RSUM11L=0.D0

            RSUM12=0.D0

            IF(MODELA.eq.0)THEN
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)
                  END DO
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                          *RLAMDA(NTH,NR)*RFSADG(NR)
                  END DO
               ENDDO
            END IF

!           DEFINE BULK MOMENTUM RANGE: 0 < NP < NP_BULK(NR,NSA) 
            IF(NT_init.eq.0)THEN
               PMAX_BULK = FACT_BULK
            ELSE
               PMAX_BULK = SQRT( RT_BULK(NR,NSA)/RTFP(NR,NSA) )*FACT_BULK
            END IF
            RATIO_0_S=SQRT(RTFPS(NSA)/RTFP0(NSA)) ! p_th(r)/p_th(0)
            P_BULK_R = PMAX_BULK* ( (1.D0-RATIO_0_S)*(1.D0-RM(NR)**2)+RATIO_0_S )
            DO NP=NPMAX, 1, -1
               IF(P_BULK_R.lt.PG(NP,NSA)) NP_BULK(NR,NSA)=NP
            END DO
            IF(PMAX(NSA).lt.P_BULK_R)THEN
               NP_BULK(NR,NSA) = NPMAX
            END IF
!            IF(NRANK.eq.0) WRITE(*,'(A,3I5,3E14.6)') "TEST NP_BULK =",NSA, NR, NP_BULK(NR,NSA),FACT_BULK,PMAX_BULK,P_BULK_R
               
!
            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *0.5D0*PM(NP,NSBA)**2
                     END DO
                     RSUMNP_E(NP)=RSUM3
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)/PV
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *(PV-1.D0)/THETA0(NSA)
                     END DO
                     RSUMNP_E(NP)=RSUM3
                  END DO
               ENDIF
            ELSE ! MODELA=1
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)*RLAMDA(NTH,NR)*RFSADG(NR)
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *0.5D0*PM(NP,NSBA)**2*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                     RSUMNP_E(NP)=RSUM3
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1, ITL(NR)-1
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)/PV &
                             *(1.D0+EPSRM2(NR))
                     END DO
                     DO NTH=ITU(NR)+1, NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)/PV &
                             *(1.D0+EPSRM2(NR))
                     END DO
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *(PV-1.D0)/THETA0(NSA)*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                     RSUMNP_E(NP)=RSUM3
                  END DO
               ENDIF               
            END IF
! --------- current density without runaway region 
!            NPS=NPSTART
!            NPE=NPEND
!            IF(NPC_runaway-1.le.NPEND) NPE=NPC_runaway-1
!            IF(NPC_runaway-1.le.NPSTART) NPS=NPSTART
!            RSUM1R=0.D0
!            RSUM2R=0.D0
!            IF(MODELR.eq.0)THEN
!               DO NP=NPS, NPE
!                  DO NTH=1,NTHMAX
!                     RSUM1R = RSUM1R+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)                     
!                     RSUM2R = RSUM2R                       &
!                          +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
!                          *PM(NP,NSBA)*COSM(NTH)
!                  END DO
!               END DO
!            ELSEIF(MODELR.eq.1)THEN
!               DO NP=NPS, NPE
!                  PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
!                  DO NTH=1,NTHMAX
!                     RSUM1R = RSUM1R+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)                     
!                     RSUM2R = RSUM2R                       &
!                          +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
!                          *PM(NP,NSBA)*COSM(NTH)/PV
!                  END DO
!               END DO
!            END IF
! --------- 
            IF(NPSTART.eq.1)THEN
               NPS=2
            ELSE
               NPS=NPSTART
            END IF
            DO NP=NPS,NPEND
               PV=SQRT(1.D0+THETA0(NSA)*PG(NP,NSBA)**2)
               DO NTH=1,NTHMAX
                  WPL=WEIGHP(NTH  ,NP,NR,NSA)
                  IF(NTH.EQ.1) THEN
                     WPM=0.D0
                  ELSE
                     WPM=WEIGHP(NTH-1,NP,NR,NSA)
                  ENDIF
                  IF(NTH.EQ.NTHMAX) THEN
                     WPP=0.D0
                  ELSE
                     WPP=WEIGHP(NTH+1,NP,NR,NSA)
                  ENDIF
                  DFP=    PG(NP,NSBA) &
                       /DELP(NSBA)*(FNSP(NTH,NP,NR,NSBA)-FNSP(NTH,NP-1,NR,NSBA))
                  IF(NTH.EQ.1) THEN
                     DFT=1.D0/DELTH                             &
                          *(                                     &
                          ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSBA)   &
                          +WPP *FNSP(NTH+1,NP-1,NR,NSBA)) &
                          -                                    &
                          ((1.D0-WPM)*FNSP(NTH,NP  ,NR,NSBA)     &
                          +WPM *FNSP(NTH,NP-1,NR,NSBA)) &
                          )
                     
                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/DELTH                         & 
                          *(-                                    &
                          ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSBA)   &
                          +WPM *FNSP(NTH-1,NP-1,NR,NSBA)) &
                          +                                     &
                          ((1.D0-WPP)*FNSP(NTH,NP  ,NR,NSBA)     &
                          +WPP *FNSP(NTH,NP-1,NR,NSBA)) &
                          )
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)                  &
                          *(                                     &
                          ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSBA)   &
                          +WPP *FNSP(NTH+1,NP-1,NR,NSBA)) &
                          -                                    &
                          ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSBA)   &
                          +WPM *FNSP(NTH-1,NP-1,NR,NSBA)) &
                          )
                  ENDIF
                  FFP=    PG(NP,NSBA)                           &
                       *((1.D0-WPL)*FNSP(NTH  ,NP  ,NR,NSBA)  &
                       +WPL *FNSP(NTH  ,NP-1,NR,NSBA))
!
                  RSUM4 = RSUM4+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                       *(DCPP(NTH,NP,NR,NSA)*DFP           &
                        +DCPT(NTH,NP,NR,NSA)*DFT           &
                        -FCPP(NTH,NP,NR,NSA)*FFP           &
                          )
                  RSUM5 = RSUM5+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                          *(DWPP(NTH,NP,NR,NSA)*DFP           &
                           +DWPT(NTH,NP,NR,NSA)*DFT)
                  RSUM6 = RSUM6-PG(NP,NSBA)**2*SINM(NTH)/PV   &
                          *(FEPP(NTH,NP,NR,NSA)*FFP)
                  RSUM7 = RSUM7+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                          *(DWLHPP(NTH,NP,NR,NSA)*DFP         &
                           +DWLHPT(NTH,NP,NR,NSA)*DFT)
                  RSUM8 = RSUM8+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                          *(DWFWPP(NTH,NP,NR,NSA)*DFP         &
                       +DWFWPT(NTH,NP,NR,NSA)*DFT)
                  RSUM9 = RSUM9+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                          *(DWECPP(NTH,NP,NR,NSA)*DFP         &
                       +DWECPT(NTH,NP,NR,NSA)*DFT)
                  RSUM_IC = RSUM_IC+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                              *(DWICPP(NTH,NP,NR,NSA)*DFP         &
                               +DWICPT(NTH,NP,NR,NSA)*DFT)
                  RSUM_synch = RSUM_synch + &
                       PG(NP,NSBA)**2*SINM(NTH)/PV   &
                       *( -FSPP(NTH,NP,NR,NSA)*FFP )
                  RSUM_loss = RSUM_loss + &
                       PG(NP,NSBA)**2*SINM(NTH)/PV   &
                       *(DLPP(NTH,NP,NR,NSA)*DFP           &
                        -FLPP(NTH,NP,NR,NSA)*FFP )

                  DO NSB=1,NSBMAX
                     RSUM10(NSB)=RSUM10(NSB)+PG(NP,NSBA)**2*SINM(NTH)/PV &
                             *(DCPP2(NTH,NP,NR,NSB,NSA)*DFP              &
                              +DCPT2(NTH,NP,NR,NSB,NSA)*DFT              & 
                              -FCPP2(NTH,NP,NR,NSB,NSA)*FFP)
                  END DO
                  RSUM12 = RSUM12+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                             *(DPP(NTH,NP,NR,NSA)*DFP           &
                              +DPT(NTH,NP,NR,NSA)*DFT           &
                              -FPP(NTH,NP,NR,NSA)*FFP           &
                              )
               ENDDO ! NTH
            ENDDO

!      SOURCE POWER
            IF(MODELR.eq.1)THEN
               DO NP=NPSTART,NPEND
                  PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                  DO NTH=1,NTHMAX
                     RSUM11B = RSUM11B + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)*SPPB(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
                     RSUM11F = RSUM11F + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)*SPPF(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
!                     IF(MODEL_SINK.eq.0)THEN
!                        RSUM11S = RSUM11S + PM(NP,NSBA)**2*SINM(NTH) &
!                             *(PV-1.D0)/THETA0(NSA)*SPPS(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
!                     ELSE
                        RSUM11S = RSUM11S + PM(NP,NSBA)**2*SINM(NTH) &
                             *(PV-1.D0)/THETA0(NSA)*SPPL(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
!                     END IF
                     RSUM11L = RSUM11L + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)* PPL(NTH,NP,NR,NSA)&
                          *FNSP(NTH,NP,NR,NSBA)!*RLAMDA(NTH,NR)!*RFSADG(NR) PPL includes RLAMDA
                  END DO
               END DO
            ELSE ! MODELR=0
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM11B = RSUM11B + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*SPPB(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
                     RSUM11F = RSUM11F + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*SPPF(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
!                     IF(MODEL_SINK.eq.0)THEN
!                        RSUM11S = RSUM11S + PM(NP,NSBA)**2*SINM(NTH) &
!                             *0.5D0*PM(NP,NSBA)**2*SPPS(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
!                     ELSE
                        RSUM11S = RSUM11S + PM(NP,NSBA)**2*SINM(NTH) &
                             *0.5D0*PM(NP,NSBA)**2*SPPL(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
!                     END IF
                     RSUM11L = RSUM11L + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*PPL(NTH,NP,NR,NSA) &
                          *FNSP(NTH,NP,NR,NSBA)!*RLAMDA(NTH,NR)!*RFSADG(NR)
                  END DO
               END DO
            END IF


!-------    Calculation of bulk temperature
            CALL BULK_TEMPERATURE(NP_BULK(NR,NSA),NR,NSA) 
            CALL mtx_set_communicator(comm_np) 

!-------    Calculation of bulk temperature

! --------- radial transport 
            RSUM_DR=0.D0
            RSUM_DRS=0.D0
            RSUMN_DR=0.D0
            RSUMN_DRS=0.D0
            DINT_DFDT_R1=0.D0
            DINT_DFDT_R2=0.D0
            SRHOR1 = 0.D0
            SRHOR2 = 0.D0
            IF(MODELA.eq.0)THEN
               DRRM=RG(NR)
               DRRP=RG(NR+1)
               RL=RM(NR)
            ELSE
               DRRM=1.D0
               DRRP=1.D0
               RL=1.D0
            END IF
            IF(MODELD.ne.0)THEN
               DO NP=NPSTART,NPEND
                  RGAMA=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                  DO NTH=1,NTHMAX
                     WRL=WEIGHR(NTH,NP,NR,NSA)
                     WRH=WEIGHR(NTH,NP,NR+1,NSA)
                     IF(NR.ne.1.and.NR.ne.NRMAX)THEN
                        DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR
                        F_R1 = ( (1.D0-WRL)*FNSP(NTH,NP,NR,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) )
                        
                        DFDR_R2 = ( FNSP(NTH,NP,NR+1,NSA)-FNSP(NTH,NP,NR,NSA) ) / DELR
                        F_R2 = ( (1.D0-WRH)*FNSP(NTH,NP,NR+1,NSA) + WRH*FNSP(NTH,NP,NR,NSA) )

                        DFDT_R1 = ( DRR(NTH,NP,NR,NSA)*DFDR_R1 - FRR(NTH,NP,NR,NSA)*F_R1 )*DRRM
                        DFDT_R2 = ( DRR(NTH,NP,NR+1,NSA)*DFDR_R2 - FRR(NTH,NP,NR+1,NSA)*F_R2)*DRRP
                     ELSEIF(NR.eq.1)THEN
                        F_R2 = ( (1.D0-WRH)*FNSP(NTH,NP,NR+1,NSA) + WRH*FNSP(NTH,NP,NR,NSA) )
                        DFDR_R2 = ( FNSP(NTH,NP,NR+1,NSA)-FNSP(NTH,NP,NR,NSA) ) / DELR

                        DFDT_R1 = 0.D0
                        DFDT_R2 = ( DRR(NTH,NP,NR+1,NSA)*DFDR_R2 - FRR(NTH,NP,NR+1,NSA)*F_R2)*DRRP
                     ELSEIF(NR.eq.NRMAX)THEN
                        F_R1 = ( (1.D0-WRL)*FNSP(NTH,NP,NR,NSA) + WRL*FNSP(NTH,NP,NR-1,NSA) )
                        DFDR_R1 = ( FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP,NR-1,NSA) ) / DELR

                        F_R2 = ( (1.D0-WRH)*FS2(NTH,NP,NSA) + WRH*FNSP(NTH,NP,NR,NSA) )
                        DFDR_R2 = ( FS2(NTH,NP,NSA)-FNSP(NTH,NP,NR,NSA) ) / DELR
                        
                        DFDT_R1 = ( DRR(NTH,NP,NR,NSA)*DFDR_R1 - FRR(NTH,NP,NR,NSA)*F_R1 )*DRRM
                        DFDT_R2 = ( DRR(NTH,NP,NR+1,NSA)*DFDR_R2 - FRR(NTH,NP,NR+1,NSA)*F_R2)*DRRP
                     END IF
                     
                     DINT_DR = ( DFDT_R2 - DFDT_R1 )/DELR/RL*RFSADG(NR)

                     RSUMN_DR=RSUMN_DR +             DINT_DR*VOLP(NTH,NP,NSBA)
                     IF(MODELR.eq.1)THEN
                        RSUM_DR=RSUM_DR+(RGAMA-1.D0)/THETA0(NSA)*DINT_DR*VOLP(NTH,NP,NSBA)
                     ELSEIF(MODELR.eq.0)THEN
                        RSUM_DR=RSUM_DR + 0.5D0*PG(NP,NSBA)**2*DINT_DR*VOLP(NTH,NP,NSBA)
                     END IF
!
                     IF(NR.eq.NRMAX)THEN
                        RSUMN_DRS=RSUMN_DRS +      DFDT_R2*VOLP(NTH,NP,NSBA)/DELR *RFSADG(NR)/RL
                        IF(MODELR.eq.1)THEN
                           RSUM_DRS=RSUM_DRS+(RGAMA-1.D0)/THETA0(NSA)*DFDT_R2*VOLP(NTH,NP,NSBA)/DELR *RFSADG(NR)
                        ELSE
                           RSUM_DRS=RSUM_DRS+0.5D0*PG(NP,NSBA)**2*DFDT_R2*VOLP(NTH,NP,NSBA)/DELR *RFSADG(NR)
                        END IF
                     END IF
                  END DO

               END DO
!     REDUCE RSUM
               CALL p_theta_integration(RSUMN_DR)
               CALL p_theta_integration(RSUM_DR)

               RNDRL(NR,NSA)=RNFP0(NSA)*RSUMN_DR/(DELR*RM(NR))!*RFSADG(NR)
               RPDRL(NR,NSA)=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*1.D-6 &
                    *RSUM_DR/(DELR*RM(NR))!*RFSADG(NR)
               IF(NR.eq.NRMAX)THEN
                  CALL p_theta_integration(RSUMN_DRS)
                  CALL p_theta_integration(RSUM_DRS)
                  RNDRS(NSA) = RNFP0(NSA)*RSUMN_DRS/(DELR*RM(NR))!* 2.D0
                  RPDRS(NSA) = RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*1.D-6 &
                       *RSUM_DRS/(DELR*RM(NR))
               END IF               
            END IF ! MODELD

! --------- end of radial transport
            CALL p_theta_integration(RSUM1)
            CALL p_theta_integration(RSUM2)
            CALL p_theta_integration(RSUM3)
            CALL p_theta_integration(RSUM4)
            CALL p_theta_integration(RSUM5)
            CALL p_theta_integration(RSUM6)
            CALL p_theta_integration(RSUM7)
            CALL p_theta_integration(RSUM8)
            CALL p_theta_integration(RSUM9)
            CALL p_theta_integration(RSUM_IC)
            CALL p_theta_integration(RSUM_synch)
            CALL p_theta_integration(RSUM_loss)
            CALL p_theta_integration(RSUM11B)
            CALL p_theta_integration(RSUM11F)
            CALL p_theta_integration(RSUM11S)
            CALL p_theta_integration(RSUM11L)
            DO NSB=1,NSBMAX
               CALL p_theta_integration(RSUM10(NSB))
            END DO
!            CALL p_theta_integration(RSUM1R)

 888        FORMAT(2I4,12E14.6)
            FACT=RNFP0(NSA)*1.D20
            RNSL(NR,NSA) = RSUM1*FACT*1.D-20
            RJSL(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)*1.D-6

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RWSL(NR,NSA) = RSUM3*FACT               *1.D-6
            RWS123L(NR,NSA) =-RSUM12*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*RFSADG(NR)
            RPCSL(NR,NSA)=-RSUM4*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPWSL(NR,NSA)=-RSUM5*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPESL(NR,NSA)=-RSUM6*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RLHSL(NR,NSA)=-RSUM7*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RFWSL(NR,NSA)=-RSUM8*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RECSL(NR,NSA)=-RSUM9*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RICSL(NR,NSA)=-RSUM_IC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPSSL(NR,NSA)=-RSUM_synch*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPLSL(NR,NSA)=-RSUM_loss*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            DO NSB=1,NSBMAX
               RPCS2L(NR,NSB,NSA)=-RSUM10(NSB) &
                    *FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            END DO

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*RFSADG(NR)
            RSPBL(NR,NSA)= RSUM11B*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPFL(NR,NSA)= RSUM11F*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPSL(NR,NSA)= RSUM11S*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPLL(NR,NSA)= RSUM11L*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
         ENDDO ! NSA
      ENDDO ! NR
      CALL mtx_reset_communicator 

      CALL FPSAVECOMM

      CALL GUTIME(gut2)

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            IF(RNS(NR,NSA).NE.0.D0) THEN
               IF(MODELR.eq.0)THEN
                  RT_T(NR,NSA) = RWS(NR,NSA)*1.D6 &
                       /(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
               ELSEIF(MODELR.eq.1)THEN
                  CALL FPNEWTON(NR,NSA,RNS(NR,NSA),RWS(NR,NSA),rtemp)
                  RT_T(NR,NSA) = rtemp
               END IF
            ELSE
               RT_T(NR,NSA) = 0.D0
            ENDIF
         END DO
      END DO

!      DO NSA=1,NSAMAX
!        DO NR=NRSTART,NREND
!            RNFP(NR,NSA)=RNS(NR,NSA)
!            RTFP(NR,NSA)=RT_T(NR,NSA)
!            NSB=NSB_NSA(NSA)
!            IF(NSB.NE.0) THEN
!               RNFD(NR,NSB)=RNFP(NR,NSA)
!               RTFD(NR,NSB)=RTFP(NR,NSA)
!            END IF
!         END DO
!      END DO

!      IF(NRANK.eq.0) WRITE(6,'(A,2E16.8)') "FNSP ",FNSP(1,2,1,1), FNSP(NTHMAX,2,1,1)
!      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "-----FPSSUB=", gut2-gut1

      ISAVE=1
      RETURN
      END SUBROUTINE FPSSUB
!
! *************************
!     SAVE GLOBAL DATA
! *************************
!
      SUBROUTINE FPSGLB
!
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH 
      real(8):: EAVE, EAVE2, rtemp, rtemp2, THETAL, THETAL2

      NTG1=NTG1+1
      call fp_adjust_ntg1

      PTG(NTG1)=TIMEFP

      DO NSA=1,NSAMAX
         PNT(NSA,NTG1)=0.D0
         PIT(NSA,NTG1)=0.D0
         PWT(NSA,NTG1)=0.D0
         PPCT(NSA,NTG1)=0.D0
         PPWT(NSA,NTG1)=0.D0
         PPET(NSA,NTG1)=0.D0
         PLHT(NSA,NTG1)=0.D0
         PFWT(NSA,NTG1)=0.D0
         PECT(NSA,NTG1)=0.D0
         PICT(NSA,NTG1)=0.D0
         PWT2(NSA,NTG1)=0.D0
         PSPBT(NSA,NTG1)=0.D0
         PSPFT(NSA,NTG1)=0.D0
         PSPST(NSA,NTG1)=0.D0
         PSPLT(NSA,NTG1)=0.D0
         PWTD(NSA,NTG1)=0.D0
         DO NSB=1,NSBMAX
            PPCT2(NSB,NSA,NTG1)= 0.D0
         END DO
         PDR(NSA,NTG1)=0.D0
         PNDR(NSA,NTG1)=0.D0
         PTT_BULK(NSA,NTG1)=0.D0
      ENDDO

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            PNT(NSA,NTG1) =PNT(NSA,NTG1) +RNS(NR,NSA)*VOLR(NR)
            PIT(NSA,NTG1) =PIT(NSA,NTG1) +RJS(NR,NSA)*VOLR(NR)
            PWT(NSA,NTG1) =PWT(NSA,NTG1) +RWS(NR,NSA)*VOLR(NR)
            PWTD(NSA,NTG1)=PWTD(NSA,NTG1)+RWS123(NR,NSA)*VOLR(NR)
            PPCT(NSA,NTG1)=PPCT(NSA,NTG1)+RPCS(NR,NSA)*VOLR(NR)
            PPWT(NSA,NTG1)=PPWT(NSA,NTG1)+RPWS(NR,NSA)*VOLR(NR)
            PPET(NSA,NTG1)=PPET(NSA,NTG1)+RPES(NR,NSA)*VOLR(NR)
            PLHT(NSA,NTG1)=PLHT(NSA,NTG1)+RLHS(NR,NSA)*VOLR(NR)
            PFWT(NSA,NTG1)=PFWT(NSA,NTG1)+RFWS(NR,NSA)*VOLR(NR)
            PECT(NSA,NTG1)=PECT(NSA,NTG1)+RECS(NR,NSA)*VOLR(NR)
            PICT(NSA,NTG1)=PICT(NSA,NTG1)+RICS(NR,NSA)*VOLR(NR)
            PSPBT(NSA,NTG1)=PSPBT(NSA,NTG1)+RSPB(NR,NSA)*VOLR(NR)
            PSPFT(NSA,NTG1)=PSPFT(NSA,NTG1)+RSPF(NR,NSA)*VOLR(NR)
            PSPST(NSA,NTG1)=PSPST(NSA,NTG1)+RSPS(NR,NSA)*VOLR(NR)
            PSPLT(NSA,NTG1)=PSPLT(NSA,NTG1)+RSPL(NR,NSA)*VOLR(NR)
            PDR(NSA,NTG1) = PDR(NSA,NTG1) + RPDR(NR,NSA)*VOLR(NR)
            PNDR(NSA,NTG1) = PNDR(NSA,NTG1) + RNDR(NR,NSA)*VOLR(NR)
            PTT_BULK(NSA,NTG1)=PTT_BULK(NSA,NTG1) + RT_BULK(NR,NSA)*VOLR(NR)*RNS(NR,NSA)
            PPST(NSA,NTG1)=PPST(NSA,NTG1)+RPSS(NR,NSA)*VOLR(NR)
            PPLT(NSA,NTG1)=PPLT(NSA,NTG1)+RPLS(NR,NSA)*VOLR(NR)

            IF(MODELR.eq.1) then
               CALL FPNEWTON(NR,NSA,RNS(NR,NSA),RWS(NR,NSA),rtemp)
            else
               rtemp=0.D0
               rtemp2=0.D0
            endif
            PWT2(NSA,NTG1) =PWT2(NSA,NTG1) +rtemp*VOLR(NR)/1.D6 &
                                  *(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
            DO NSB=1,NSBMAX
               PPCT2(NSB,NSA,NTG1)=PPCT2(NSB,NSA,NTG1)       &
                                  +RPCS2(NR,NSB,NSA)*VOLR(NR)
            END DO
         ENDDO

         PIT(NSA,NTG1) =PIT(NSA,NTG1)/(2.D0*PI*RR)
         PSPT(NSA,NTG1)=PSPBT(NSA,NTG1)+PSPFT(NSA,NTG1) &
                       +PSPST(NSA,NTG1)+PSPLT(NSA,NTG1)

         IF(PNT(NSA,NTG1).NE.0.d0) THEN
            PTT(NSA,NTG1) =PWT(NSA,NTG1)*1.D6   &
                 /(1.5D0*PNT(NSA,NTG1)*1.D20*AEE*1.D3)
            PTT2(NSA,NTG1) =PWT2(NSA,NTG1)*1.D6 &
                 /(1.5D0*PNT(NSA,NTG1)*1.D20*AEE*1.D3)
!            IF(MODELR.eq.1)THEN
!               PTT_BULK(NSA,NTG1) =PWT_BULK2(NSA,NTG1)*1.D6 &
!                    /(1.5D0*PNT_BULK(NSA,NTG1)*1.D20*AEE*1.D3)
!            ELSE
!               PTT_BULK(NSA,NTG1) =PWT_BULK(NSA,NTG1)*1.D6 &
!                    /(1.5D0*PNT_BULK(NSA,NTG1)*1.D20*AEE*1.D3)               
!            END IF
         ELSE
            PTT(NSA,NTG1)=0.D0
            PTT2(NSA,NTG1)=0.D0
!            PTT_BULK(NSA,NTG1)=0.D0
         ENDIF
         PTT_BULK(NSA,NTG1) = PTT_BULK(NSA,NTG1)/PNT(NSA,NTG1)
         PNT(NSA,NTG1) =PNT(NSA,NTG1)/TVOLR
!         PNT_BULK(NSA,NTG1) =PNT_BULK(NSA,NTG1)/TVOLR
!         PNDR(NSA,NTG1)=PNDR(NSA,NTG1)/TVOLR
      ENDDO

!      IF(MODEL_NBI.NE.0) THEN
!        Q_ENG(NTG1)=5*PSPFT(4,NTG1)/PSPBT(2,NTG1)
!         WRITE(6,*) 'Q_ENG', Q_ENG(NTG1)
!      ENDIF         

      RETURN
      END SUBROUTINE FPSGLB
!
! *************************
!     SAVE PROFILE DATA
! *************************
!
      SUBROUTINE FPSPRF
!
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NS
      real(8):: RS, rtemp, rtemp2, EAVE2, THETAL2

      NTG2=NTG2+1
      call fp_adjust_ntg2

      RTG(NTG2)=TIMEFP

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RNT(NR,NSA,NTG2) = RNS(NR,NSA)
            rate_runaway(NR,NTG2)=RFP(NR)
            RJT(NR,NSA,NTG2) = RJS(NR,NSA)
            RWT(NR,NSA,NTG2) = RWS(NR,NSA)
            RPCT(NR,NSA,NTG2)= RPCS(NR,NSA)
            RPWT(NR,NSA,NTG2)= RPWS(NR,NSA)
            RPET(NR,NSA,NTG2)= RPES(NR,NSA)
            RLHT(NR,NSA,NTG2)= RLHS(NR,NSA)
            RFWT(NR,NSA,NTG2)= RFWS(NR,NSA)
            RECT(NR,NSA,NTG2)= RECS(NR,NSA)
            RICT(NR,NSA,NTG2)= RICS(NR,NSA)
            DO NSB=1,NSBMAX
               RPCT2(NR,NSB,NSA,NTG2)= RPCS2(NR,NSB,NSA)
            END DO
            RSPBT(NR,NSA,NTG2)= RSPB(NR,NSA)
            RSPFT(NR,NSA,NTG2)= RSPF(NR,NSA)
            RSPST(NR,NSA,NTG2)= RSPS(NR,NSA)
            RSPLT(NR,NSA,NTG2)= RSPL(NR,NSA)
            RPDRT(NR,NSA,NTG2)= RPDR(NR,NSA)
            RNDRT(NR,NSA,NTG2)= RNDR(NR,NSA)
            RTT_BULK(NR,NSA,NTG2) = RT_BULK(NR,NSA)
!
            IF(RNS(NR,NSA).NE.0.D0) THEN
               IF(MODELR.eq.0)THEN
                  RTT(NR,NSA,NTG2) = RWS(NR,NSA)*1.D6 &
                       /(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
               ELSEIF(MODELR.eq.1)THEN
                  CALL FPNEWTON(NR,NSA,RNS(NR,NSA),RWS(NR,NSA),rtemp)
                  RTT(NR,NSA,NTG2) = rtemp
               END IF
            ELSE
               RTT(NR,NSA,NTG2) = 0.D0
            ENDIF
            RET(NR,NTG2) = E1(NR)
            RS=RSRHON(RM(NR))
            RQT(NR,NTG2) = RS*BB*2.D0/(RR*(BP(NR)+BP(NR+1)))
            RT_T(NR,NSA) = RTT(NR,NSA,NTG2)
         ENDDO
      ENDDO
      DO NR=1,NRMAX
         EPTR(NR,NTG2)=E1(NR)
      END DO

      RETURN
      END SUBROUTINE FPSPRF
!
! ***********************************************************
!
!                         RESULT
!
! ***********************************************************
!
      SUBROUTINE FPWRTGLB
!
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH
      real(8):: rtotalPW, rtotalPC,rtotalSP,rtotalPC2
      real(8):: rtotalDR,rtotalEC,rtotalIC,rtotalIP
      character:: fmt0*50
!
      WRITE(6,*)"--------------------------------------------"
      WRITE(6,*)"-----Global data"
      WRITE(6,101) PTG(NTG1)*1000

      DO NSA=1,NSAMAX
         IF(MODELR.eq.0)THEN
            WRITE(6,112) NSA,NS_NSA(NSA), &
              PNT(NSA,NTG1),PTT(NSA,NTG1),PWT(NSA,NTG1),PIT(NSA,NTG1),RNDRS(NSA),PTT_BULK(NSA,NTG1)
         ELSE
            WRITE(6,112) NSA,NS_NSA(NSA), &
              PNT(NSA,NTG1),PTT2(NSA,NTG1),PWT(NSA,NTG1),PIT(NSA,NTG1),RNDRS(NSA),PTT_BULK(NSA,NTG1)
         END IF
      ENDDO

      rtotalPW=0.D0
      rtotalPC=0.D0
      rtotalSP=0.D0
      rtotalPC2=0.D0
      rtotalEC=0.D0
      rtotalIC=0.D0
      rtotalIP=0.D0

      WRITE(fmt0,'(a48)') &
           '(8X,2I2," PC,PW,PE,PDR,PPLT=",6X,1P5E12.4)'
      DO NSA=1,NSAMAX
         WRITE(6,fmt0) NSA,NS_NSA(NSA), &
              PPCT(NSA,NTG1),PPWT(NSA,NTG1),PPET(NSA,NTG1),RPDRS(NSA),PPLT(NSA,NTG1)
         rtotalPW=rtotalPW + PPWT(NSA,NTG1)
         rtotalPC=rtotalPC + PPCT(NSA,NTG1)
         rtotalSP=rtotalSP + PSPT(NSA,NTG1)
         rtotalPC2 = rtotalPC2 +PPCT(NSA,NTG1)-PPCT2(NSA,NSA,NTG1)
         rtotalEC=rtotalEC + PECT(NSA,NTG1)
         rtotalIC=rtotalIC + PICT(NSA,NTG1)
         rtotalIP=rtotalIP + PIT(NSA,NTG1)
      END DO
      DO NSA=1,NSAMAX
         IF(NSBMAX.GT.1) THEN
            write(6,104) NSA,NS_NSA(NSA), &
                 (PPCT2(NSB,NSA,NTG1),NSB=1,NSBMAX)
         ENDIF
      END DO

      DO NSA=1,NSAMAX
         write(6,108) NSA,NS_NSA(NSA),PSPBT(NSA,NTG1),PSPFT(NSA,NTG1), &
                                      PSPST(NSA,NTG1),PSPLT(NSA,NTG1), &
                                      PSPST(NSA,NTG1)+PSPLT(NSA,NTG1)
      END DO
      write(6,105) rtotalpw,rtotalEC,rtotalIC
      write(6,107) rtotalPC
      write(6,109) rtotalSP
      write(6,110) rtotalPC2
      WRITE(6,'("total plasma current   [MA]",1PE12.4)') rtotalIP

      RETURN

  101 FORMAT(' TIME=',F12.3,' ms')
  112 FORMAT(' NSA,NS=',2I2,' n,T,W,I,dn,T2=',1PE11.4,1P6E12.4)
  104 FORMAT('        ',2I2,' PCAB    =',10X,1P14E12.4)

 105  FORMAT('Total absorption power [MW]', 1PE12.4,'    EC:',1PE12.4,'    IC:',1PE12.4)
 106  FORMAT(F12.4, 8E12.4)
 107  FORMAT('total collision power  [MW]', 1PE12.4)
 108  FORMAT('        ',2I2,' PSPB/F/S/L/S+L=',4X,1P5E12.4) 
 109  FORMAT('total source power     [MW]', 1PE12.4)
 110  FORMAT('collision balance      [MW]', 1PE12.4)

      END SUBROUTINE FPWRTGLB

! ***********************************************************

      SUBROUTINE FPWRTPRF
!
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH
      real(8):: RTFDL, RTFD0L, THETAL, rtemp, rtemp2, temp
      character:: fmt0*50
!
!      WRITE(fmt0,'(a15)') '(2I3,1P20E13.4)'
!      WRITE(fmt0,'(a44)') '(2I3,1P8E16.8,1P5E12.3e3,1PE12.4,1P5E12.3e3)'

      WRITE(6,*)"-----Radial profile data"
      WRITE(6,'(A,F12.3)') " TIME=", TIMEFP*1000

      IF(MODEL_DISRUPT.eq.0)THEN
         WRITE(fmt0,'(a44)') '(2I3,1P8E12.4,1P5E12.3e3,1PE12.4,1P5E12.3e3)'
         WRITE(6,106) 
      ELSE
         WRITE(fmt0,'(a44)') '(2I3,1P8E12.4,1P5E12.3e3,1PE12.4,1P5E12.3e3)'
         WRITE(6,108) 
      END IF

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            IF(MODELR.eq.1)THEN
               CALL FPNEWTON(NR,NSA,RNS(NR,NSA),RWS(NR,NSA),rtemp)
               RTT(NR,NSA,NTG2)=rtemp
            END IF

            IF(MODEL_DISRUPT.ne.0)THEN
               WRITE(6,fmt0) NSA,NS_NSA(NSA),&
                    RM(NR),RNT(NR,NSA,NTG2),RTT(NR,NSA,NTG2), &
                    POST_tau_ta(NR,NSA),     &
                    RSPS(NR,NSA), &
                    E1(NR), conduct_sp(NR), RN_disrupt(NR),   &
                    RN_runaway(NR), RN_runaway(NR)-RN_drei(NR), RJ_ohm(NR), RJ_runaway(NR), &
                    RJ_bs(NR), &
                    RT_quench(NR), &
                    RFP(NR), Rconnor(NR), RFP_ava(NR), &
                    ER_crit(NR)
            ELSE
!               IF(E1(NR).NE.0.D0) THEN
!                  TEMP=RJS(NR,NSA)/E1(NR)
!               ELSE
!                  TEMP=0.D0
!               ENDIF
               WRITE(6,fmt0) NSA,NS_NSA(NSA),&
                    RM(NR),RNT(NR,NSA,NTG2),RTT(NR,NSA,NTG2), &
                    RJT(NR,NSA,NTG2),RPCT(NR,NSA,NTG2),       &
                    RPET(NR,NSA,NTG2), &
                    RPWT(NR,NSA,NTG2),&
                    RSPBT(NR,NSA,NTG2), &
!                    RDIDT(NR,NSA), &
!                    RPCT2(NR,NSAMAX-NSA+1,NSA,NTG2), &
!                    RPCT2(NR,NSA,NSA,NTG2), &
!                    RECT(NR,NSA,NTG2),    &
!RSPFT(NR,NSA,NTG2),RPDRT(NR,NSA,NTG2), &
                    RPDR(NR,NSA), &
                    RT_BULK(NR,NSA)!, &
!                    TEMP, conduct_sp(NR)
            END IF
         ENDDO
      ENDDO
      RETURN
  106 FORMAT( &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j    ',5X,'PC     ',5X,'PE     ',5X,   &
           'PW ',9X,'PNBI',8X,'PDRP',8X,'TBULK' )
  107 FORMAT( &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j     ',5X,'PC     ',5X,'PC12     ',5X,  &
           'PC11   ',5X,'PE     ',5X,'E_ind    ',5X,  &
           'PE_IND ',4X,'PSIP   ',5X,'DPSIP ',5X,'SIGMA' )
  108 FORMAT( &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' NBI   ',5X,'j_fp   ',5X,'E1    ',5X,  &
           'sigma  ',5X,'n_bulk ',5X,'n_run  ',5X,'n_second',4X,  &
           'j_ohm  ',5X,'j_run  ',5X,'j_bs   ',5X,'T      ',5X, &
           'r_rate ',5X,'Rconnor',5X,'rate_a ',5X,'E_C')

      END SUBROUTINE FPWRTPRF
! ***********************************************************
!
      SUBROUTINE FPWRTSNAP
!
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NP, NTH
      real(8):: rnute, resist
      real(8):: FACT1, FACT2, rntv
      real(8):: taue_col, sigma_sp, FACT
      real(8):: taue_col2, sigma_sp2

!-----check of conductivity--------
      IF(NTG1.ne.1.and.NRANK.eq.0)then
         NR=1
         NSA=1
         NSB=2
!         Do NSA=1,NSAMAX
!         Do NSB=1,NSBMAX
            rnute=RNUD(NR,NSA,NSA)*SQRT(2.D0)*RNFP(NR,NSA)/RNFP0(NSA)     &
                 *(PTFP0(NSA)/PTFD(1,NSA))**2
            resist=RNUTE*AMFP(NSA)/RNFP(1,NSA)/AEFP(NSA)**2/1.D20
!----------
            FACT=AEFP(NSA)**2*AEFD(NSB)**2*LNLAM(NR,NSB,NSA)/(EPS0**2)
            taue_col=3.D0*SQRT((2.D0*PI)**3)/FACT*SQRT(AMFP(1)*(AEE*RTT(NR,NSA,NTG2)*1.D3)**3)/RNS(NR,NSB)*1.D-20
            taue_col2=3.D0*SQRT((2.D0*PI)**3)/FACT*SQRT(AMFP(1)*(AEE*PTT_BULK(NSA,NTG1)*1.D3)**3)/RNS(NR,NSB)*1.D-20
            sigma_sp=1.96D0*RNS(NR,NSA)*1.D20*AEFP(NSA)**2*taue_col/AMFP(NSA) ! P. 174
            sigma_sp2=1.96D0*RNS(NR,NSA)*1.D20*AEFP(NSA)**2*taue_col2/AMFP(NSA) ! P. 174
!            IF(MODELA.eq.0)THEN
!               sigma_sp=sigma_sp*ZEFF*(1.D0+0.27D0*ZEFF-0.27D0)/(1.D0+0.47D0*ZEFF-0.47D0)
!            END IF

            if(E0.ne.0)THEN 
               write(6,'(A,1PE16.8)') &
                    " THETA0= ", THETA0(NSA)
               write(6,'(A,1PE16.8,A,1PE16.8,A,1PE16.8,A,1PE16.8)') &
                    " whole space sigma_norm= ",(RJS(1,1)+RJS(1,2))/E1(1)*1.D6*resist, &
                    " sigma=j/E*1.D6=", (RJS(1,NSA)+RJS(1,2))/E1(1)*1.D6
               write(6,'(A,1PE16.8,A,1PE16.8)') &
                    " sigma_Spitzer= ",sigma_sp, &
                    " sigma_Spitzer_bulk=", sigma_sp2
               WRITE(6,'(A,1PE14.6,A,1PE14.6,A,1PE14.6)') &
                    " j_norm=j/q_e/n(t)/v_ta0=", RJS(1,NSA)/(ABS(AEFP(1))*RNS(1,NSA)*1.D20*VTFP0(1))*1.D6, &
                    " j_norm=sigma_norm*E0=",(RJS(1,1))/E1(1)*1.D6*resist*E0, &
                    " j_norm=(j_e+j_i)=", (RJS(1,1)+RJS(1,2))/(ABS(AEFP(1))*RNS(1,NSA)*1.D20*VTFP0(1))*1.D6
               WRITE(6,'(A,1PE14.6,A,1PE14.6)') &
                    " sigma_bulk_karney =    ", &
                    ( RJS(1,NSA)/(ABS(AEFP(1))*RNS(1,NSA)*1.D20*VTFP0(1))*1.D6 &
                    -0.5D0*(rate_runaway(1,NTG2)/E0 )*PG(NPMAX+1,1)**2 )/ &
                    E0/resist, &
                    " J_bulk_Karney=", &
                    RJS(1,NSA)/(ABS(AEFP(1))*RNS(1,NSA)*1.D20*VTFP0(1))*1.D6 &
                    -0.5D0*(rate_runaway(1,NTG2)/E0 )*PG(NPMAX+1,1)**2

            END if
!         END DO
!         END DO
      END IF
!----end of conductivity check---------

      RETURN

      END SUBROUTINE FPWRTSNAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine FPNEWTON(NR,NSA,RNSL_,RWSL_,rtemp)

      IMPLICIT NONE
      INTEGER,INTENT(IN)::NR,NSA
      double precision,intent(in):: RNSL_, RWSL_ 
      double precision,intent(out)::rtemp
      integer:: NSB, NP, NTH, ncount
      real(8):: rtemp2, xeave, xeave2
      real(8):: xtemp, xtemp2, thetal,thetal2, EAVE, EAVE2

!-----Average kinetic energy
!      EAVE=RWS(NR,NSA)*AMFP(NSA)*THETA0(NSA) &
!           /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
      EAVE=RWSL_*AMFP(NSA)*THETA0(NSA) & 
           /(RNSL_*1.D20*PTFP0(NSA)**2*1.D-6)

!-----initial value of THETAL
      THETAL=2.d0*EAVE/3.d0
      xtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)

      CALL XNEWTON(EAVE,THETAL,ncount)

      rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
      xeave=AMFP(NSA)*VC**2*EAVE/(AEE*1.D3)

      RETURN

      CONTAINS

      SUBROUTINE xnewton(eave,thetal,ncount)
      IMPLICIT NONE
      INTEGER NR,NSA
      REAL(8),intent(in):: eave
      REAL(8),intent(inout):: thetal
      INTEGer,intent(out):: ncount
      REAL(8),parameter:: eps=1.d-10
      REAL(8):: delthetal,thetalnew,epsthetal

!--------iteration start
      ncount=0
      DO while(ncount.le.100)
         ncount=ncount+1
         delthetal=-(rfunc(thetal)-eave)/dfunc(thetal)
         thetalnew=thetal+delthetal
         epsthetal=ABS(delthetal/thetal)

         thetal=thetalnew
         IF(epsthetal.le.eps) EXIT
      END DO
      RETURN
      END SUBROUTINE xnewton

      FUNCTION rfunc(thetal)
      IMPLICIT NONE
      INTEGER NR,NSA
      REAL(8):: thetal,rfunc
      REAL(8):: z,dkbsl1,dkbsl2
      z=1.D0/thetal
      dkbsl1=BESEKN(1,Z)
      dkbsl2=BESEKN(2,Z)
      rfunc= (dkbsl1 /dkbsl2 -1.D0+3.D0/Z)
      RETURN
      END FUNCTION rfunc

      FUNCTION rfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,rfuncp
      REAL(8):: z,dkbsl1,dkbsl2
      z=1.D0/thetal
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      rfuncp= dkbsl1 /dkbsl2 -1.D0+3.D0/Z
      RETURN
      END FUNCTION rfuncp
      
      FUNCTION dfunc(thetal)
      IMPLICIT NONE
      INTEGER NR,NSA
      REAL(8):: thetal,dfunc
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3
      z=1.D0/thetal
      dkbsl0=BESEKN(0,z)
      dkbsl1=BESEKN(1,z)
      dkbsl2=BESEKN(2,z)
      dkbsl3=BESEKN(3,z)
      dfunc =( (dkbsl0 +dkbsl2 )/dkbsl2                           &
                -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2 &
            +3.d0
      RETURN
      END FUNCTION dfunc

      FUNCTION dfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,dfuncp
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3
      z=1.D0/thetal
      dkbsl0=1.D0 -  1.D0/8.D0/z +   9.D0/128.D0/z**2
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      dkbsl3=1.D0 + 35.D0/8.D0/z + 945.D0/128.D0/z**2
      dfuncp =( (dkbsl0 +dkbsl2 )/dkbsl2                          &
               -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2  & 
            +3.d0
      RETURN
      END FUNCTION dfuncp
      
      end Subroutine FPNEWTON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      SUBROUTINE FPSAVECOMM

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NSW, N
      double precision,dimension(:),allocatable::DDATA
      double precision,dimension(:),allocatable::DDATA_R
      INTEGER:: NOFFSET, NR2, NLENGTH, NALL
      INTEGER,dimension(nsize)::COMPOS, COMLEN
      double precision,dimension(NSAMAX*2)::tempr_in, tempr_out
      double precision,dimension(NRSTART:NREND,NSAMAX):: work
      double precision,dimension(NRMAX,NSAMAX):: workg
      double precision,dimension(NSAMAX):: temp_nsanr
      double precision,dimension(NPMAX,NRSTART:NREND):: temp_npnr1
      double precision,dimension(NPMAX,NRMAX):: temp_npnr2
      double precision,dimension(NPMAX,NRMAX,NSASTART:NSAEND):: temp_npnr3
      integer,dimension(NSAMAX):: vloc
      INTEGER:: isave_sw

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL,SAVLEN(NRANK+1),RNS,N,NSA)
!         CALL fp_gatherv_real8_sav(RFPL,SAVLEN(NRANK+1),RFP,N,NSA)
         CALL fp_gatherv_real8_sav(RJSL,SAVLEN(NRANK+1),RJS,N,NSA)
         CALL fp_gatherv_real8_sav(RWSL,SAVLEN(NRANK+1),RWS,N,NSA)
         CALL fp_gatherv_real8_sav(RWS123L,SAVLEN(NRANK+1),RWS123,N,NSA)
         CALL fp_gatherv_real8_sav(RPCSL,SAVLEN(NRANK+1),RPCS,N,NSA)
         CALL fp_gatherv_real8_sav(RPWSL,SAVLEN(NRANK+1),RPWS,N,NSA)
         CALL fp_gatherv_real8_sav(RPESL,SAVLEN(NRANK+1),RPES,N,NSA)
         CALL fp_gatherv_real8_sav(RLHSL,SAVLEN(NRANK+1),RLHS,N,NSA)
         CALL fp_gatherv_real8_sav(RFWSL,SAVLEN(NRANK+1),RFWS,N,NSA)
         CALL fp_gatherv_real8_sav(RECSL,SAVLEN(NRANK+1),RECS,N,NSA)
         CALL fp_gatherv_real8_sav(RICSL,SAVLEN(NRANK+1),RICS,N,NSA)
         CALL fp_gatherv_real8_sav(RSPBL,SAVLEN(NRANK+1),RSPB,N,NSA)
         CALL fp_gatherv_real8_sav(RSPFL,SAVLEN(NRANK+1),RSPF,N,NSA)
         CALL fp_gatherv_real8_sav(RSPSL,SAVLEN(NRANK+1),RSPS,N,NSA)
         CALL fp_gatherv_real8_sav(RSPLL,SAVLEN(NRANK+1),RSPL,N,NSA)
         CALL fp_gatherv_real8_sav(RPDRL,SAVLEN(NRANK+1),RPDR,N,NSA)
         CALL fp_gatherv_real8_sav(RNDRL,SAVLEN(NRANK+1),RNDR,N,NSA)
         CALL fp_gatherv_real8_sav(RTL_BULK,SAVLEN(NRANK+1),RT_BULK,N,NSA) 
         CALL fp_gatherv_real8_sav(RDIDTL,SAVLEN(NRANK+1),RDIDT,N,NSA)
         CALL fp_gatherv_real8_sav(RPSSL,SAVLEN(NRANK+1),RPSS,N,NSA)
         CALL fp_gatherv_real8_sav(RPLSL,SAVLEN(NRANK+1),RPLS,N,NSA)
      END DO
      DO NSB=1,NSBMAX
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               work(NR,NSA) = RPCS2L(NR,NSB,NSA)
            END DO
         END DO
         DO N=1,NSW
            NSA=N+NSASTART-1
            CALL fp_gatherv_real8_sav(work,SAVLEN(NRANK+1),workg,N,NSA)
         ENDDO

         DO NSA=1,NSAMAX
            DO NR=1,NRMAX
               RPCS2(NR,NSB,NSA) = workg(NR,NSA)
            END DO
         END DO
      ENDDO
      CALL mtx_reset_communicator 

      CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RWS,NRMAX*NSAMAX)
!      CALL mtx_broadcast_real8(RFP,NRMAX)

      IF(MODELD.ne.0)THEN
         temp_nsanr(:)=0.D0
         CALL mtx_set_communicator(comm_nsanr) 
         CALL mtx_reduce_real8(RPDRS,NSAMAX,3,temp_nsanr,vloc)
         RPDRS(:)=temp_nsanr
         CALL mtx_reduce_real8(RNDRS,NSAMAX,3,temp_nsanr,vloc)
         RNDRS(:)=temp_nsanr
         CALL mtx_reset_communicator 
      END IF

      CALL mtx_set_communicator(comm_nr) 
      DO NSA=NSASTART, NSAEND
         DO NR=NRSTART, NREND
            DO NP=1,NPMAX
               temp_npnr1(NP,NR) = RPL_BULK(NP,NR,NSA)
            END DO
         END DO
         CALL mtx_allgather_real8(temp_npnr1,NPMAX*(NREND-NRSTART+1),temp_npnr2)
         DO NR=1, NRMAX
            DO NP=1, NPMAX
               temp_npnr3(np,nr,nsa) = temp_npnr2(np,nr)
            END DO
         END DO
      END DO
      CALL mtx_set_communicator(comm_nsa) 
      CALL mtx_allgather_real8(temp_npnr3,NPMAX*NRMAX*(NSAEND-NSASTART+1),RP_BULK)
      CALL mtx_reset_communicator 

!      IF(NRANK.eq.0)THEN
!         NSA=2
!         NR=1
!         WRITE(19,'(A,E14.6)') "# RT_BULK", TIMEFP
!         WRITE(19,*) " "
!         WRITE(19,*) " "
!         DO NP=NPSTART, NPEND
!            IF(RP_BULK(NP,NR,NSA).ne.0.D0)THEN
!               WRITE(19,'(I4,2E14.6)') NP, RP_BULK(NP,1,2), RP_BULK(NP,2,2)
!            END IF
!         END DO
!      END IF

      END SUBROUTINE FPSAVECOMM
!^------------------------------ 
!==============================================================
      SUBROUTINE FPSSUB2
!
      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NPS
      integer:: IERR
      real(8):: RSUM1, RSUM3, fact
      real(8):: PV, WPL, WPM, WPP
      real:: gut1,gut2

      CALL mtx_set_communicator(comm_np) 
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)

            RSUM1=0.D0
            RSUM3=0.D0

            IF(MODELA.eq.0)THEN
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)
                  END DO
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                          *RLAMDA(NTH,NR)*RFSADG(NR) 
                  END DO
               ENDDO
            END IF
!
            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *0.5D0*PM(NP,NSBA)**2
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *(PV-1.D0)/THETA0(NSA)
                     END DO
                  END DO
               ENDIF
            ELSE
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *0.5D0*PM(NP,NSBA)**2*RLAMDA(NTH,NR)*RFSADG(NR) 
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *(PV-1.D0)/THETA0(NSA)*RLAMDA(NTH,NR)*RFSADG(NR) 
                     END DO
                  END DO
               ENDIF 
            END IF
!     REDUCE RSUM
            CALL p_theta_integration(RSUM1)
            CALL p_theta_integration(RSUM3)
! --------- end of radial transport
               
            FACT=RNFP0(NSA)*1.D20
            RNSL(NR,NSA) = RSUM1*FACT*1.D-20

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RWSL(NR,NSA) = RSUM3*FACT               *1.D-6
         ENDDO ! NSA
      ENDDO ! NR
      CALL mtx_reset_communicator 

      CALL FPSAVECOMM2

      RETURN
      END SUBROUTINE FPSSUB2
!
! *************************
!     SAVE PROFILE DATA
! *************************
!
      SUBROUTINE FPSPRF2
!
      USE plprof
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NS
      real(8):: rtemp, rhon
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF

      DO NS=1,NSMAX
         DO NR=1,NRMAX
            RHON=RM(NR) 
            CALL PL_PROF(RHON,PLF) 
            RN_TEMP(NR,NS)=PLF(NS)%RN
            RT_TEMP(NR,NS)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
         ENDDO
      END DO
      Do NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         IF(NS.NE.0)THEN
            DO NR=1,NRMAX
               RN_TEMP(NR,NS) = RNS(NR,NSA)
               RT_TEMP(NR,NS) = RT_BULK(NR,NSA)
!               RT_TEMP(NR,NS) = RT_T(NR,NSA)
            END DO
         END IF
      ENDDO

      RETURN
      END SUBROUTINE FPSPRF2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FPSAVECOMM2

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NSW, N

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL,SAVLEN(NRANK+1),RNS,N,NSA)
         CALL fp_gatherv_real8_sav(RWSL,SAVLEN(NRANK+1),RWS,N,NSA)
      END DO
      CALL mtx_reset_communicator 

      CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RWS,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RT_BULK,NRMAX*NSAMAX)

      END SUBROUTINE FPSAVECOMM2
!==============================================================
      SUBROUTINE update_RN_RT

      IMPLICIT NONE

      CALL FPSSUB2
      CALL FPSPRF2

      END SUBROUTINE update_RN_RT
!==============================================================
      SUBROUTINE BULK_TEMPERATURE(NPB,NR,NSA)

      USE fpmpi
      IMPLICIT NONE
      INTEGER,intent(in):: NPB, NR, NSA
!      REAL(8),INTENT(OUT):: RTL
      integer:: ISW_BULK, NP, NTH, NS, NSBA
      real(8):: RSUM_T, RSUM_V, PV, DFDP, WPL, FFP, RSUMN, RSUMW, FACT
      real(8),dimension(NTHMAX,NPMAX):: T_BULK
      real(8),dimension(NPSTART:NPEND):: RPL_BULK_send 
      real(8),dimension(NPMAX):: RPL_BULK_recv
      real(8):: RNL_BULK, RWL_BULK, rtemp

      ISW_BULK=0 ! sometimes DFDP becomes 0 and then it makes density NaN
!      ISW_BULK=1 ! requires higher FACT_BULK (FACT_BULK > 5) to obtain accurate T_bulk

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

      CALL mtx_set_communicator(comm_np)
      IF(ISW_BULK.eq.0)THEN ! BULK T calculation using dfdp
         RSUM_T=0.D0
         RSUM_V=0.D0
         DO NP=NPSTART,NPEND
            IF(NP.ge.2.and.NP.le.NPB)THEN
               PV=SQRT(1.D0+THETA0(NSA)*PG(NP,NSBA)**2)
               DO NTH=1,NTHMAX
                  IF(FNSP(NTH,NP,NR,NSA).gt.0.D0.and.FNSP(NTH,NP-1,NR,NSA).gt.0.D0)THEN
                     DFDP=DELP(NSA)/ &
                          ( log(FNSP(NTH,NP,NR,NSA))-log(FNSP(NTH,NP-1,NR,NSA)) )
                  ELSE
                     WPL=WEIGHP(NTH  ,NP,NR,NSA)
                     FFP=   ( (1.D0-WPL)*FNSP(NTH  ,NP  ,NR,NSBA)  &
                          +WPL *FNSP(NTH  ,NP-1,NR,NSBA) )
                     DFDP=DELP(NSA)*FFP/(                         &
                          FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP-1,NR,NSA) )
                  END IF
                  IF(MODEL_DISRUPT.ne.1)THEN
                     IF(DFDP.ne.DFDP)THEN ! NaN never equals to any other variables.
!                        WRITE(16,'(A,E14.6,3I4,3E14.6)') "DFDP is NaN in fpsave. TIMEFP= ", TIMEFP, NP, NTH, NSA, &
!                             FNSP(NTH,NP,NR,NSA), FNSP(NTH,NP-1,NR,NSA), FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP-1,NR,NSA)
                        DFDP=0.D0
                     END IF
                     IF(DFDP.gt.0.D0)THEN
!                        WRITE(16,'(A,E14.6,3I4,3E14.6)') "DFDP is positive in fpsave. TIMEFP= ", TIMEFP, NP, NTH, NSA, &
!                             FNSP(NTH,NP,NR,NSA), FNSP(NTH,NP-1,NR,NSA), FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP-1,NR,NSA)
                        DFDP=0.D0
                     END IF
                  END IF
                  
                  T_BULK(NTH,NP)=-PG(NP,NSA)*PTFP0(NSA)*DFDP &
                       /AEE/1.D3*VTFP0(NSA)/PV 
                  RSUM_T = RSUM_T + T_BULK(NTH,NP)*VOLP(NTH,NP,NSBA)
                  RSUM_V = RSUM_V + VOLP(NTH,NP,NSBA)
               END DO
               RPL_BULK_send(NP) = RSUM_T/RSUM_V
            ELSE
               RPL_BULK_send(NP) = 0.D0
            END IF
         END DO
         CALL mtx_allgather_real8(RPL_BULK_send,NPEND-NPSTART+1,RPL_BULK_recv)
         DO NP=1, NPMAX
            RPL_BULK(NP,NR,NSA) = RPL_BULK_recv(NP)
         END DO
         CALL p_theta_integration(RSUM_T)
         CALL p_theta_integration(RSUM_V)
         
         RTL_BULK(NR,NSA)=RSUM_T/RSUM_V
      ELSE ! ISW_BULK=1  BULK T using T=W/n 
         RSUMN=0.D0
         RSUMW=0.D0

         IF(MODELA.eq.0)THEN
            DO NP=NPSTART,NPEND
               IF(NP.le.NPB)THEN
                  DO NTH=1,NTHMAX
                     RSUMN = RSUMN+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)
                  END DO
               END IF
            ENDDO
         ELSE
            DO NP=NPSTART,NPEND
               IF(NP.le.NPB)THEN
                  DO NTH=1,NTHMAX
                     RSUMN = RSUMN+VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                          *RLAMDA(NTH,NR)*RFSADG(NR)
                  END DO
               END IF
            ENDDO
         END IF
         IF(MODELA.eq.0) THEN
            IF(MODELR.EQ.0) THEN
               DO NP=NPSTART,NPEND
                  IF(NP.le.NPB)THEN
                     DO NTH=1,NTHMAX
                        RSUMW = RSUMW                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *0.5D0*PM(NP,NSBA)**2
                     END DO
                  END IF
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  IF(NP.le.NPB)THEN
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUMW = RSUMW                       &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA) &
                             *(PV-1.D0)/THETA0(NSA)
                     END DO
                  END IF
               END DO
            ENDIF
         ELSE ! MODELA=1
            IF(MODELR.EQ.0) THEN
               DO NP=NPSTART,NPEND
                  IF(NP.le.NPB)THEN
                     DO NTH=1,NTHMAX
                        RSUMW = RSUMW                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *0.5D0*PM(NP,NSBA)**2*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                  END IF
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  IF(NP.le.NPB)THEN
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUMW = RSUMW                        &
                             +VOLP(NTH,NP,NSBA)*FNSP(NTH,NP,NR,NSBA)  &
                             *(PV-1.D0)/THETA0(NSA)*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                  END IF
               END DO
            ENDIF
         END IF
         CALL p_theta_integration(RSUMN)
         CALL p_theta_integration(RSUMW)

         FACT=RNFP0(NSA)*1.D20 
         RNL_BULK = RSUMN*FACT*1.D-20
         FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
         RWL_BULK = RSUMW*FACT*1.D-6

         IF(RNL_BULK.NE.0.D0) THEN
            IF(MODELR.eq.0)THEN
               RTL_BULK(NR,NSA) = RWL_BULK*1.D6 &
                    /(1.5D0*RNL_BULK*1.D20*AEE*1.D3)
            ELSEIF(MODELR.eq.1)THEN
               CALL FPNEWTON(NR,NSA,RNL_BULK,RWL_BULK,rtemp)
               RTL_BULK(NR,NSA) = rtemp
            END IF
         END IF

!         WRITE(*,'(A,I3,3E14.6)') "TEST", NRANK, RWL_BULK, RNL_BULK, RTL_BULK(NR,NSA)
         DO NP=1, NPMAX
            RPL_BULK(NP,NR,NSA) = 0.D0
         END DO
      END IF
      CALL mtx_reset_communicator 
      

      END SUBROUTINE BULK_TEMPERATURE
!==============================================================


      end MODULE fpsave




