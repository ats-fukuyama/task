!     $Id$
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
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS
      integer:: IERR
      real(8):: RSUM1, RSUM2, RSUM3, RSUM4, RSUM5, RSUM6, RSUM7, RSUM_FS2
      real(8):: RSUM8, RSUM9, RSUM123, RSUM11B,RSUM11F,RSUM11S,RSUM11L
      real(8):: RSUM12, RSUM_IC, rsum_test
      real(8):: PV, WPL, WPM, WPP
      real(8):: DFP, DFT, FFP, testa, testb, FACT, DFDP, WRL, WRH, DINT_DFDT_R1, DINT_DFDT_R2
      real(8):: DFDR_R1, DFDR_R2, DFDT_R1, DFDT_R2, DINT_DR, RSUM_DR,RGAMA,F_R1,F_R2, RSUMN_DR
      real(8),dimension(NSBMAX):: RSUM10
      real(8),dimension(NTHMAX+1,NPMAX+1):: R_FLUX
      real(8),dimension(NPMAX):: RSUMNP_E
      real(8),dimension(NTHMAX,NPMAX):: T_BULK
      integer,dimension(NRSTART:NRENDX,NSBMAX):: NP_BULK
      real(8):: ratio, RSUM_T, RSUM_V, P_BULK_R, FACT_BULK, RATIO_0_S
      real(8):: SRHOR1, SRHOR2, RSUM_DRS, RSUMN_DRS
!      real(8),dimension(NRMAX, NSAMAX):: RFWS, RECS
!      real(8),dimension(NRMAX, NSBMAX, NSAMAX):: RPCS2
!      real(8),dimension(NSAMAX, 0:NTMAX):: PWT2, PTT2
      INTEGER,dimension(NPROCS)::NMTXPOS, NMTXLEN
      real(8):: DSDR, SRHOP, SRHOM

      IF(ISAVE.NE.0) RETURN

      DO NSA=1,NSAMAX
         RNDRS(NSA) =0.D0
         RPDRS(NSA) =0.D0
      END DO
      DO NR=NRSTART,NRENDX
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)

            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            RSUM5=0.D0
            RSUM6=0.D0
            RSUM7=0.D0
            RSUM8=0.D0
            RSUM9=0.D0
            RSUM_IC=0.D0
            DO NSB=1,NSBMAX
               RSUM10(NSB)=0.D0
            ENDDO
            RSUM11B=0.D0
            RSUM11F=0.D0
            RSUM11S=0.D0
            RSUM11L=0.D0

            RSUM12=0.D0
!            RSUM_test=0.D0

            IF(MODELA.eq.0)THEN
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA)
                  END DO
               ENDDO
            ELSE
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA) &
                          *RLAMDA(NTH,NR)
                  END DO
               ENDDO
            END IF

!           DEFINE BULK MOMENTUM RANGE
            FACT_BULK=3.D0
            RATIO_0_S=SQRT(RTFPS(NSA)/RTFP0(NSA)) ! p_th(r)/p_th(0)
            DO NP=NPMAX, 1, -1
               P_BULK_R = FACT_BULK* ( (1.D0-RATIO_0_S)*(1.D0-RM(NR)**2)+RATIO_0_S )
               IF(PG(NP,NSA).gt.P_BULK_R) NP_BULK(NR,NSA)=NP
            END DO
!            WRITE(*,*) NSA, NR, NP_BULK(NR,NSA)
               
!
            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=1,NPMAX
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA) &
                             *0.5D0*PM(NP,NSBA)**2
                     END DO
                  RSUMNP_E(NP)=RSUM3
                  ENDDO
               ELSE
                  DO NP=1,NPMAX
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)/PV
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA) &
                             *(PV-1.D0)/THETA0(NSA)
                     END DO
                  RSUMNP_E(NP)=RSUM3
                  END DO
               ENDIF
            ELSE
               IF(MODELR.EQ.0) THEN
                  DO NP=1,NPMAX
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)*RLAMDA(NTH,NR)
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA)  &
                             *0.5D0*PM(NP,NSBA)**2*RLAMDA(NTH,NR)
                     END DO
                  RSUMNP_E(NP)=RSUM3
                  ENDDO
               ELSE
                  DO NP=1,NPMAX
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)/PV*RLAMDA(NTH,NR)
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSBA)  &
                             *(PV-1.D0)/THETA0(NSA)*RLAMDA(NTH,NR)
                     END DO
                  RSUMNP_E(NP)=RSUM3
                  END DO
               ENDIF               
            END IF

!            open(8,file='F_DFDP_r0c4_w.dat')
!            open(8,file='nth1_r1c4_x.dat')
!            open(8,file='T_BULK.dat')
            DO NP=2,NPMAX
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
                       /DELP(NSBA)*(FNS1(NTH,NP,NR,NSBA)-FNS1(NTH,NP-1,NR,NSBA))
                  IF(NTH.EQ.1) THEN
                     DFT=1.D0/DELTH                             &
                         *(                                     &
                            ((1.D0-WPP)*FNS1(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNS1(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNS1(NTH,NP  ,NR,NSBA)     &
                                  +WPM *FNS1(NTH,NP-1,NR,NSBA))&
                          )

                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/DELTH                         & 
                         *(-                                    &
                            ((1.D0-WPM)*FNS1(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNS1(NTH-1,NP-1,NR,NSBA))&
                          +                                     &
                            ((1.D0-WPP)*FNS1(NTH,NP  ,NR,NSBA)     &
                                  +WPP *FNS1(NTH,NP-1,NR,NSBA))&
                          )
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)                  &
                         *(                                     &
                            ((1.D0-WPP)*FNS1(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNS1(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNS1(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNS1(NTH-1,NP-1,NR,NSBA))&
                                  )
                  ENDIF
                  FFP=    PG(NP,NSBA)                           &
                         *((1.D0-WPL)*FNS1(NTH  ,NP  ,NR,NSBA)  &
                                +WPL *FNS1(NTH  ,NP-1,NR,NSBA))
                  RSUM4 = RSUM4+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                       *(DCPP(NTH,NP,NR,NSA)*DFP           &
                        +DCPT(NTH,NP,NR,NSA)*DFT           &
                        -FCPP(NTH,NP,NR,NSA)*FFP           &
                          )
                  R_FLUX(NTH,NP) = -  &
                       (DPP(NTH,NP,NR,NSA)*DFP           &
                       +DPT(NTH,NP,NR,NSA)*DFT           &
                       -FPP(NTH,NP,NR,NSA)*FFP           &
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
               ENDDO
            ENDDO

!      SOURCE POWER
            IF(MODELR.eq.1)THEN
               DO NP=1,NPMAX
                  PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                  DO NTH=1,NTHMAX
                     RSUM11B = RSUM11B + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)*SPPB(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)
                     RSUM11F = RSUM11F + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)*SPPF(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)
                     RSUM11S = RSUM11S + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)*SPPS(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)
                     RSUM11L = RSUM11L + PM(NP,NSBA)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NSA)* PPL(NTH,NP,NR,NSA)&!*RLAMDA(NTH,NR) &
                          *FNS1(NTH,NP,NR,NSBA)
                  END DO
               END DO
            ELSE ! MODELR=0
               DO NP=1,NPMAX
                  PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                  DO NTH=1,NTHMAX
                     RSUM11B = RSUM11B + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*RLAMDAG(NTH,NR)*SPPB(NTH,NP,NR,NSA)
                     RSUM11F = RSUM11F + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*RLAMDAG(NTH,NR)*SPPF(NTH,NP,NR,NSA)
                     RSUM11S = RSUM11S + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*RLAMDAG(NTH,NR)*SPPS(NTH,NP,NR,NSA)
                     RSUM11L = RSUM11L + PM(NP,NSBA)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NSBA)**2*RLAMDAG(NTH,NR)*PPL(NTH,NP,NR,NSA) &
                          *FNS1(NTH,NP,NR,NSBA)
                  END DO
               END DO
            END IF


!-------    Calculation of bulk temperature
            RSUM_T=0.D0
            RSUM_V=0.D0
            DO NP=2,NP_BULK(NR,NSA)
               PV=SQRT(1.D0+THETA0(NSA)*PG(NP,NSBA)**2)
               DO NTH=1,NTHMAX
                  IF(FNS1(NTH,NP,NR,NSA).gt.0.D0.and.FNS1(NTH,NP-1,NR,NSA).gt.0.D0)THEN
                     DFDP=DELP(NSA)/ &
                       ( log(FNS1(NTH,NP,NR,NSA))-log(FNS1(NTH,NP-1,NR,NSA)) )
                  ELSE
                     WPL=WEIGHP(NTH  ,NP,NR,NSA)
                     FFP=   ( (1.D0-WPL)*FNS1(NTH  ,NP  ,NR,NSBA)  &
                          +WPL *FNS1(NTH  ,NP-1,NR,NSBA) )
                     DFDP=DELP(NSA)*FFP/(                         &
                          FNS1(NTH,NP,NR,NSA)-FNS1(NTH,NP-1,NR,NSA) )
                  END IF
                  T_BULK(NTH,NP)=-PG(NP,NSA)*PTFP0(NSA)*DFDP &
                          /AEE/1.D3*VTFP0(NSA)/PV 
                  RSUM_T = RSUM_T + T_BULK(NTH,NP)*VOLP(NTH,NP,NSBA)
                  RSUM_V = RSUM_V + VOLP(NTH,NP,NSBA)
               END DO
            END DO

            RTL_BULK(NR,NSA)=RSUM_T/RSUM_V
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

            DO NP=1,NPMAX
               RGAMA=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
               DO NTH=1,NTHMAX
                  WRL=WEIGHR(NTH,NP,NR,NSA)
                  WRH=WEIGHR(NTH,NP,NR+1,NSA)
                  IF(NR.ne.1.and.NR.ne.NRMAX)THEN
                     DFDR_R1 = ( FNS(NTH,NP,NR,NSA)&!*RLAMDAG(NTH,NR)/RFSADG(NR) &
                          -FNS(NTH,NP,NR-1,NSA)&!*RLAMDAG(NTH,NR-1)/RFSADG(NR-1) &
                          ) / DELR
                     F_R1 = ( (1.D0-WRL)*FNS(NTH,NP,NR,NSA) + WRL*FNS(NTH,NP,NR-1,NSA) ) !&
!                          *RLAMDA_GG(NTH,NR)/RFSAD_GG(NR)

                     DFDR_R2 = ( FNS(NTH,NP,NR+1,NSA)&!*RLAMDAG(NTH,NR+1)/RFSADG(NR+1) &
                          -FNS(NTH,NP,NR,NSA)&!*RLAMDAG(NTH,NR)/RFSADG(NR) &
                          ) / DELR
                     F_R2 = ( (1.D0-WRH)*FNS(NTH,NP,NR+1,NSA) + WRH*FNS(NTH,NP,NR,NSA) )!&
!                          *RLAMDA_GG(NTH,NR+1)/RFSAD_GG(NR+1)

                     DFDT_R1 = ( DRR(NTH,NP,NR,NSA)*DFDR_R1 - FRR(NTH,NP,NR,NSA)*F_R1 )*RG(NR)
                     DFDT_R2 = ( DRR(NTH,NP,NR+1,NSA)*DFDR_R2 - FRR(NTH,NP,NR+1,NSA)*F_R2)*RG(NR+1)
                  ELSEIF(NR.eq.1)THEN
                     F_R2 = ( (1.D0-WRH)*FNS(NTH,NP,NR+1,NSA) + WRH*FNS(NTH,NP,NR,NSA) )!&
!                          *RLAMDA_GG(NTH,NR+1)/RFSAD_GG(NR+1)
                     DFDR_R2 = ( FNS(NTH,NP,NR+1,NSA)&!*RLAMDAG(NTH,NR+1)/RFSADG(NR+1) &
                          -FNS(NTH,NP,NR,NSA)&!*RLAMDAG(NTH,NR)/RFSADG(NR) &
                          ) / DELR

                     DFDT_R1 = 0.D0
                     DFDT_R2 = ( DRR(NTH,NP,NR+1,NSA)*DFDR_R2 - FRR(NTH,NP,NR+1,NSA)*F_R2)*RG(NR+1)
                  ELSEIF(NR.eq.NRMAX)THEN
                     F_R1 = ( (1.D0-WRL)*FNS(NTH,NP,NR,NSA) + WRL*FNS(NTH,NP,NR-1,NSA) )!&
!                          *RLAMDA_GG(NTH,NR)/RFSAD_GG(NR)
                     DFDR_R1 = ( FNS(NTH,NP,NR,NSA)&!*RLAMDAG(NTH,NR)/RFSADG(NR) &
                          -FNS(NTH,NP,NR-1,NSA)&!*RLAMDAG(NTH,NR-1)/RFSADG(NR-1) &
                          ) / DELR

                     F_R2 = ( (1.D0-WRH)*FS2(NTH,NP,NSA) + WRH*FNS(NTH,NP,NR,NSA) )!&
!                          *RLAMDA_GG(NTH,NR+1)/RFSAD_GG(NR+1)
                     DFDR_R2 = ( FS2(NTH,NP,NSA)&!*RLAMDAG(NTH,NR+1)/RFSADG(NR+1) &
                          -FNS(NTH,NP,NR,NSA)&!*RLAMDAG(NTH,NR)/RFSADG(NR) &
                          ) / DELR

                     DFDT_R1 = ( DRR(NTH,NP,NR,NSA)*DFDR_R1 - FRR(NTH,NP,NR,NSA)*F_R1 )*RG(NR)
                     DFDT_R2 = ( DRR(NTH,NP,NR+1,NSA)*DFDR_R2 - FRR(NTH,NP,NR+1,NSA)*F_R2)*RG(NR+1)
                  END IF

                  DINT_DR = ( DFDT_R2 - DFDT_R1 )
!                  SRHOR1 = SRHOR1 + DFDT_R1*VOLP(NTH,NP,NSBA)
!                  SRHOR2 = SRHOR2 + DFDT_R2*VOLP(NTH,NP,NSBA)

                  RSUMN_DR=RSUMN_DR +             DINT_DR*VOLP(NTH,NP,NSBA)
                  IF(MODELR.eq.1)THEN
                     RSUM_DR=RSUM_DR+(RGAMA-1.D0)/THETA0(NSA)*DINT_DR*VOLP(NTH,NP,NSBA)
                  ELSEIF(MODELR.eq.0)THEN
                     RSUM_DR=RSUM_DR + 0.5D0*PG(NP,NSBA)**2*DINT_DR*VOLP(NTH,NP,NSBA)
                  END IF
!
                  IF(NR.eq.NRMAX)THEN
                     RSUMN_DRS=RSUMN_DRS +      DFDT_R2*VOLP(NTH,NP,NSBA)
                     IF(MODELR.eq.1)THEN
                        RSUM_DRS=RSUM_DRS+(RGAMA-1.D0)/THETA0(NSA)*DFDT_R2*VOLP(NTH,NP,NSBA)
                     ELSE
                        RSUM_DRS=RSUM_DRS+0.5D0*PG(NP,NSBA)**2*DFDT_R2*VOLP(NTH,NP,NSBA)
                     END IF
                  END IF

               END DO
            END DO
!            IF(NRSTART.eq.28.or.NRSTART.eq.29.or.NRSTART.eq.30.or.NRSTART.eq.31)THEN
!               IF(NSA.eq.1) &
!               WRITE(*,*) NRSTART, NR, SRHOR1, SRHOR2
!            END IF

            RNDRL(NR,NSA)=RNFP0(NSA)*RSUMN_DR/(DELR*RM(NR))!*RFSADG(NR)
            RPDRL(NR,NSA)=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*1.D-6 &
                 *RSUM_DR/(DELR*RM(NR))!*RFSADG(NR)
            IF(NR.eq.NRMAX)THEN
               RNDRS(NSA) = RNFP0(NSA)*RSUMN_DRS/(DELR*RM(NR))!* 2.D0
               RPDRS(NSA) = RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*1.D-6 &
                 *RSUM_DRS/(DELR*RM(NR))
            END IF


! --------- end of radial transport
               
 888        FORMAT(2I4,12E14.6)
            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)*RCOEFNG(NR)
            RNSL(NR,NSA) = RSUM1*FACT                   *1.D-20
            RJSL(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)*1.D-6

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)/RFSADG(NR) *RCOEFNG(NR)
            RWSL(NR,NSA) = RSUM3*FACT               *1.D-6
            RWS123L(NR,NSA) =-RSUM12*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPCSL(NR,NSA)=-RSUM4*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPWSL(NR,NSA)=-RSUM5*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPESL(NR,NSA)=-RSUM6*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RLHSL(NR,NSA)=-RSUM7*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RFWSL(NR,NSA)=-RSUM8*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RECSL(NR,NSA)=-RSUM9*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RICSL(NR,NSA)=-RSUM_IC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            DO NSB=1,NSBMAX
               RPCS2L(NR,NSB,NSA)=-RSUM10(NSB) &
                    *FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            END DO

!            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*0.5D0!*RCOEFNG(NR)
            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)!/RFSADG(NR)*RCOEFNG(NR)
            RSPBL(NR,NSA)= RSUM11B*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPFL(NR,NSA)= RSUM11F*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPSL(NR,NSA)= RSUM11S*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPLL(NR,NSA)= RSUM11L*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
         ENDDO ! NSA
      ENDDO ! NR

      CALL FPSAVECOMM

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

            IF(MODELR.eq.1) then
               CALL FPNEWTON(NR,NSA,rtemp)
            else
!               EAVE=RWS(NR,NSA)*AMFP(NSA)*THETA0(NSA)     &
!                    /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
!               EAVE2=RWS_BULK(NR,NSA)*AMFP(NSA)*THETA0(NSA)     &
!                    /(RNS_BULK(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
!               THETAL=2.d0*EAVE/3.d0
!               THETAL2=2.d0*EAVE2/3.d0
!               rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
!               rtemp2=AMFP(NSA)*VC**2*THETAL2/(AEE*1.D3)
               rtemp=0.D0
               rtemp2=0.D0
            endif
            PWT2(NSA,NTG1) =PWT2(NSA,NTG1) +rtemp*VOLR(NR)/1.D6 &
                                  *(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
!            PWT_BULK2(NSA,NTG1) =PWT_BULK2(NSA,NTG1) +rtemp2*VOLR(NR)/1.D6 &
!                                  *(1.5D0*RNS_BULK(NR,NSA)*1.D20*AEE*1.D3)
            DO NSB=1,NSBMAX
               PPCT2(NSB,NSA,NTG1)=PPCT2(NSB,NSA,NTG1)       &
                                  +RPCS2(NR,NSB,NSA)*VOLR(NR)
            END DO
!            PDR(NSA,NTG1) = RPDR(NRMAX,NSA)*VOLR(NRMAX)
!            PDR(NSA,NTG1) = RPDRS(NSA)*TVOLR
!            PNDR(NSA,NTG1) = RNDRS(NSA)*TVOLR
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
!            RNT_test(NR,NSA,NTG2) = RNS_test(NR,NSA)
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
!                  RTT_BULK(NR,NSA,NTG2) = RWS_BULK(NR,NSA)*1.D6 &
!                       /(1.5D0*RNS_BULK(NR,NSA)*1.D20*AEE*1.D3)
               ELSEIF(MODELR.eq.1)THEN
!                  CALL FPNEWTON(NR,NSA,rtemp,rtemp2)
                  CALL FPNEWTON(NR,NSA,rtemp)
                  RTT(NR,NSA,NTG2) = rtemp
!                  RTT_BULK(NR,NSA,NTG2) = rtemp2
!                  RWS2(NR,NSA)*1.D6 &
!                       /(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
               END IF
            ELSE
               RTT(NR,NSA,NTG2) = 0.D0
            ENDIF
            RET(NR,NTG2) = E1(NR)
            RS=RSRHON(RM(NR))
            RQT(NR,NTG2) = RS*BB*2.D0/(RR*(BP(NR)+BP(NR+1)))
         ENDDO
      ENDDO

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
      real(8):: rtotalPW, rtotalPC,rtotalSP,rtotalPC2, rtotalDR,rtotalEC,rtotalIC
!      INCLUDE '../wr/wrcom1.inc'
!
      WRITE(6,*)"--------------------------------------------"
      WRITE(6,*)"-----Global data"
      WRITE(6,101) PTG(NTG1)*1000

      DO NSA=1,NSAMAX
         IF(MODELR.eq.0)THEN
            WRITE(6,112) NSA,NS_NSA(NSA), &
              PNT(NSA,NTG1),PTT(NSA,NTG1),PWT(NSA,NTG1),PIT(NSA,NTG1),PNDR(NSA,NTG1),PTT_BULK(NSA,NTG1)
         ELSE
            WRITE(6,112) NSA,NS_NSA(NSA), &
              PNT(NSA,NTG1),PTT2(NSA,NTG1),PWT(NSA,NTG1),PIT(NSA,NTG1),PNDR(NSA,NTG1),PTT_BULK(NSA,NTG1)
         END IF
!         WRITE(6,103) PPCT(NSA,NTG1),PPWT(NSA,NTG1),PPET(NSA,NTG1)
         IF(NSBMAX.GT.1) THEN
!            write(6,104)(PPCT2(NSB,NSA,NTG1),NSB=1,NSBMAX)
         ENDIF
      ENDDO

      rtotalPW=0.D0
      rtotalPC=0.D0
      rtotalSP=0.D0
      rtotalPC2=0.D0
      rtotalEC=0.D0
      rtotalIC=0.D0
      DO NSA=1,NSAMAX
!         WRITE(6,103) NSA,NS_NSA(NSA), &
!              PPCT(NSA,NTG1),PPWT(NSA,NTG1),PPET(NSA,NTG1)
         WRITE(6,113) NSA,NS_NSA(NSA), &
              PPCT(NSA,NTG1),PPWT(NSA,NTG1),PPET(NSA,NTG1),PDR(NSA,NTG1)
         rtotalPW=rtotalPW + PPWT(NSA,NTG1)
         rtotalPC=rtotalPC + PPCT(NSA,NTG1)
         rtotalSP=rtotalSP + PSPT(NSA,NTG1)
         rtotalPC2 = rtotalPC2 +PPCT(NSA,NTG1)-PPCT2(NSA,NSA,NTG1)
         rtotalEC=rtotalEC + PECT(NSA,NTG1)
         rtotalIC=rtotalIC + PICT(NSA,NTG1)
      END DO
      DO NSA=1,NSAMAX
         IF(NSBMAX.GT.1) THEN
            write(6,104) NSA,NS_NSA(NSA), &
                 (PPCT2(NSB,NSA,NTG1),NSB=1,NSBMAX)
         ENDIF
      END DO

!      DO NSA=1, NSAMAX
!         IF(NSBMAX.GT.1) THEN
!            write(6,104) NSA,NS_NSA(NSA), &
!                 (PPCT2(NSB,NSA,NTG1)+PPCT2(NSA,NSB,NTG1),NSB=1,NSBMAX)
!         ENDIF
!      END DO

      DO NSA=1,NSAMAX
         write(6,108) NSA,NS_NSA(NSA),PSPBT(NSA,NTG1),PSPFT(NSA,NTG1), &
                                      PSPST(NSA,NTG1),PSPLT(NSA,NTG1), &
                                      PSPST(NSA,NTG1)+PSPLT(NSA,NTG1)
      END DO
      write(6,105) rtotalpw,rtotalEC,rtotalIC
      write(6,107) rtotalPC
      write(6,109) rtotalSP
      write(6,110) rtotalPC2
      IF(NTG1.gt.1)THEN
         write(6,'("Steady State Criterion ",6E14.6)')  (DEPS_SS(NSA),NSA=1,NSAMAX)
      END IF

!      write(6,1111) PECT(1,NTG1)

 1001 FORMAT(I4,9E14.6)
      RETURN
  101 FORMAT(' TIME=',F12.3,' ms')
  102 FORMAT(' NSA,NS=',2I2,' n,T,W,I=',1PE11.4,1P3E12.4)
  112 FORMAT(' NSA,NS=',2I2,' n,T,W,I,dn,T2=',1PE11.4,1P5E12.4)
  103 FORMAT('        ',2I2,' PC,PW,PE=',10X,1P4E12.4)
  104 FORMAT('        ',2I2,' PCAB    =',10X,1P14E12.4)
  113 FORMAT('        ',2I2,' PC,PW,PE,PDR=',6X,1P4E12.4)

 105  FORMAT('Total absorption power [MW]', 1PE12.4,'    EC:',1PE12.4,'    IC:',1PE12.4)
 106  FORMAT(F12.4, 8E12.4)
 107  FORMAT('total collision power  [MW]', 1PE12.4)
 108  FORMAT('        ',2I2,' PSPB/F/S/L/S+L=',4X,1P5E12.4) 
 109  FORMAT('total source power     [MW]', 1PE12.4)
 110  FORMAT('collision balance      [MW]', 1PE12.4)
! 1111 FORMAT('Absorption by EC       [MW]', 1PE12.4)
      END SUBROUTINE FPWRTGLB

! ***********************************************************

      SUBROUTINE FPWRTPRF
!
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH
      real(8):: RTFDL, RTFD0L, THETAL, rtemp, rtemp2
!      INCLUDE '../wr/wrcom1.inc'
!
      WRITE(6,*)"-----Radial profile data"
      WRITE(6,106) TIMEFP*1000

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
!            RTFDL=RTFP(NR,NSA)
!            RTFD0L=(PTPR(NSA)+2.D0*PTPP(NSA))/3.D0
!            THETAL=THETA0(NSA)*RTFDL/RTFD0L

            IF(MODELR.eq.1)THEN
!               CALL FPNEWTON(NR,NSA,rtemp,rtemp2)
               CALL FPNEWTON(NR,NSA,rtemp)
               RTT(NR,NSA,NTG2)=rtemp
            END IF


!!            IF(RSPBT(NR,NSA,NTG2).GT.0.D0)THEN
!!               WRITE(6,104) NSA,NS_NSA(NSA),                  &
!!                    RM(NR),RNT(NR,NSA,NTG2),RTT(NR,NSA,NTG2), &
!!                    RJT(NR,NSA,NTG2),RPCT(NR,NSA,NTG2),       &
!!                    RPET(NR,NSA,NTG2),RPWT(NR,NSA,NTG2),      &
!!                    RSPBT(NR,NSA,NTG2),RSPFT(NR,NSA,NTG2)
!            IF(RPWT(NR,NSA,NTG2).GT.0.D0) THEN
               WRITE(6,104) NSA,NS_NSA(NSA),                  &
                    RM(NR),RNT(NR,NSA,NTG2), &
                    !RNT_test(NR,NSA,NTG2), &
                    RTT(NR,NSA,NTG2), &
                    RJT(NR,NSA,NTG2),RPCT(NR,NSA,NTG2),       &
                    RPET(NR,NSA,NTG2),RPWT(NR,NSA,NTG2),RECT(NR,NSA,NTG2),    &
                    RSPBT(NR,NSA,NTG2),RSPFT(NR,NSA,NTG2),RPDRT(NR,NSA,NTG2), &
                    RTT_BULK(NR,NSA,NTG2),RNDRT(NR,NSA,NTG2)
!                    ( RWT(NR,NSA,NTG2)-RWT(NR,NSA,NTG2-1) )/DELT, &
!                    RNDRT(NR,NSA,NTG2), &
!                    ( RNT(NR,NSA,NTG2)-RNT(NR,NSA,NTG2-1) )/DELT

!                    RPCT2(NR,1,NSA,NTG2),RPCT2(NR,2,NSA,NTG2)
!!                    RLHT(NR,NSA,NTG2),                        &
!!                    RFWT(NR,NSA,NTG2),RECT(NR,NSA,NTG2)
!            ELSE
!               WRITE(6,102) NSA,NS_NSA(NSA),                  &
!                    RM(NR),RNT(NR,NSA,NTG2),RTT(NR,NSA,NTG2), &
!                    RJT(NR,NSA,NTG2),RPCT(NR,NSA,NTG2),       &
!                    RPET(NR,NSA,NTG2)
!            ENDIF
!!               RTT(NR,NSA,NTG2)=rtemp
         ENDDO
      ENDDO
      RETURN
  101 FORMAT(' TIME=',F12.3,' ms'/                   &
           'NSA/NS',5X,'RM',10X,' n',8X,' T//PW',6X, &
           ' j//PLH',5X,'PC//PIC',5X,'PE//PEC')
  102 FORMAT(2I3,1P6E12.4)
  103 FORMAT(30X,1P4E14.6)
  104 FORMAT(2I3,1P20E12.4) 
  105 FORMAT(' TIME=',F12.3,' ms'/                   &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j     ',5X,'PC     ',5X,'PE     ')
!  106 FORMAT(' TIME=',F12.3,' ms'/                   &
!           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
!           ' j     ',5X,'PC     ',5X,'PE     ',5X,   &
!           'PW     ',5X,'PNB//PLH',5X,'PNF//PIC',5X,'PEC   ')
  106 FORMAT(' TIME=',F12.3,' ms'/                   &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j     ',5X,'PC     ',5X,'PE     ',5X,   &
           'PW     ',5X,'PEC    ',5X,'PNB//PLH',4X,  &
           ' PNF/PIC',4X,'RPDR   ',5X,'T_BULK ',5X,'RNDR' )
      END SUBROUTINE FPWRTPRF
! ***********************************************************
!
      SUBROUTINE FPWRTSNAP
!
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NP, NTH
      real(8):: rnute, resist2, resist, dtaue, dtaui 
      real(8):: FACT1, FACT2, rntv
!      INCLUDE '../wr/wrcom1.inc'
!
      WRITE(6,101) TIMEFP*1000

!$$$      if(NTG1.ne.1)then
!$$$c-----check of conductivity--------
!      IF(NTG1.ne.1)then
!!!     Spitzer
!         FACT1 =  &
!              2.D0*PI*RSRHON(RM(NRMAX))*(RSRHON(RG(NRMAX+1))-RSRHON(RG(NRMAX)))
!!              PI*RSRHON(RG(NRMAX+1))**2
!         FACT2 = 2.D0*PI*RR

!      Do NSA=1,NSAMAX
!      Do NSB=1,NSBMAX
!         rnute=RNUD(1,1,NSA)*SQRT(2.D0)*RNFP(1,NSA)/RNFP0(NSA)     &
!              *(PTFP0(NSA)/PTFD(1,NSA))**2
!         resist2=RNUTE*AMFP(NSA)/RNFP(1,NSA)/AEFP(NSA)**2/1.D20
!      resist = 1.D0/(RNFP0(NSA)/RNFP(1,NSA)*1.D20*AEFP(NSA)**2/AMFP(NSA)&
!          /RNUD(1,1,NSA)*RNFD(1,NSA)/RNFP0(NSA) ) *SQRT(2.D0)/rnfp0(NSA)
!      if(E0.ne.0.d0) &
!         write(6,*) "J/E*eta*1.D6", PIT(NSA,NTG1)/E0*1.D6/FACT1*resist2 &
!     ,"THETA0", THETA0(NSA)
!      END DO
!      END DO

!      dtaue=1.09D16*(PTT(1,NTG1))**(1.5D0)/PNT(1,NTG1)/PZ(2)**2 &
!           /15.D0*1.D-20
!      if(nsamax.gt.1) then
!      dtaui=6.6D17*(PTT(2,NTG1))**(1.5D0)/PNT(2,NTG1)/PZ(2)**4  &
!           /15.D0*1.D-20*(AMFP(2)/AMFP(2))**0.5D0
!      else
!         dtaui=0.d0
!      endif
!!      write(6,*)"tau_e tau_i[ms]",dtaue*1.D3,dtaui*1.D3

!      end if
!----end of conductivity check---------

      RETURN
!
  101 FORMAT(1H ,' TIME=',F12.3,' ms')

      END SUBROUTINE FPWRTSNAP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine FPNEWTON(NR,NSA,rtemp)

      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH, ncount
      real(8):: rtemp, rtemp2, xeave, xeave2
      real(8):: xtemp, xtemp2, thetal,thetal2, EAVE, EAVE2

!-----Average kinetic energy
      EAVE=RWS(NR,NSA)*AMFP(NSA)*THETA0(NSA) &
           /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)

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
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS
      REAL(8),dimension(:),allocatable::DDATA
      REAL(8),dimension(:),allocatable::DDATA_R
      INTEGER:: NOFFSET, NR2, NLENGTH, NALL
      INTEGER,dimension(NPROCS)::NMTXPOS, NMTXLEN
      real(8),dimension(NSAMAX*2)::tempr_in, tempr_out
      INTEGER:: isave_sw

      isave_sw=0
      IF(isave_sw.eq.1)THEN

!         NLENGTH=18
!         NOFFSET=NREND-NRSTART+1
!         NALL = NRMAX*NLENGTH
!         allocate(DDATA(NOFFSET*NLENGTH))
!         allocate(DDATA_R(NALL))
!         
!         DO NR=0,NPROCS-1
!            NMTXLEN(NR+1)=NOFFSET*NLENGTH
!            NMTXPOS(NR+1)=NOFFSET*NLENGTH*NR
!         END DO
!         DO NSA=1,NSAMAX
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2) = RNSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET) = RJSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*2) = RWSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*3) = RWS123L(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*4) = RPCSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*5) = RPWSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*6) = RPESL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*7) = RLHSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*8) = RFWSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*9) = RECSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*10) = RICSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*11) = RSPBL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*12) = RSPFL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*13) = RSPSL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*14) = RSPLL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*15) = RPDRL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*16) = RNDRL(NR,NSA)
!            END DO
!            DO NR = NRSTART, NREND
!               NR2= NR - NRSTART + 1
!               DDATA(NR2+NOFFSET*17) = RTL_BULK(NR,NSA)
!            END DO
!
!            CALL mtx_gatherv_real8(DDATA,NOFFSET*NLENGTH, &
!                 DDATA_R,NALL,NMTXLEN,NMTXPOS)
!            IF(NRANK.eq.0)THEN
!!               DO NR = 1, NRMAX
!               DO NR = 1, NPROCS
!                  RNS(NR,NSA) = DDATA_R( 1 + (NR-1)*NOFFSET*NLENGTH )
!                  RJS(NR,NSA) = DDATA_R( 2 + (NR-1)*NOFFSET*NLENGTH )
!                  RWS(NR,NSA) = DDATA_R( 3 + (NR-1)*NOFFSET*NLENGTH )
!                  RWS123(NR,NSA) = DDATA_R( 4 + (NR-1)*NOFFSET*NLENGTH )
!                  RPCS(NR,NSA) = DDATA_R( 5 + (NR-1)*NOFFSET*NLENGTH )
!                  RPWS(NR,NSA) = DDATA_R( 6 + (NR-1)*NOFFSET*NLENGTH )
!                  RPES(NR,NSA) = DDATA_R( 7 + (NR-1)*NOFFSET*NLENGTH )
!                  RLHS(NR,NSA) = DDATA_R( 8 + (NR-1)*NOFFSET*NLENGTH )
!                  RFWS(NR,NSA) = DDATA_R( 9 + (NR-1)*NOFFSET*NLENGTH )
!                  RECS(NR,NSA) = DDATA_R(10 + (NR-1)*NOFFSET*NLENGTH )
!                  RICS(NR,NSA) = DDATA_R(11 + (NR-1)*NOFFSET*NLENGTH )
!                  RSPB(NR,NSA) = DDATA_R(12 + (NR-1)*NOFFSET*NLENGTH )
!                  RSPF(NR,NSA) = DDATA_R(13 + (NR-1)*NOFFSET*NLENGTH )
!                  RSPS(NR,NSA) = DDATA_R(14 + (NR-1)*NOFFSET*NLENGTH )
!                  RSPL(NR,NSA) = DDATA_R(15 + (NR-1)*NOFFSET*NLENGTH )
!                  RPDR(NR,NSA) = DDATA_R(16 + (NR-1)*NOFFSET*NLENGTH )
!                  RNDR(NR,NSA) = DDATA_R(17 + (NR-1)*NOFFSET*NLENGTH )
!                  RT_BULK(NR,NSA) = DDATA_R(18 + (NR-1)*NOFFSET*NLENGTH )
!               END DO
!            END IF
!         END DO ! NSA
!         DO NSA=1,NSAMAX
!            tempr_in(NSA) = RNDRS(NSA)
!            tempr_in(NSA+NSAMAX) = RPDRS(NSA)
!         END DO
!         CALL mtx_reduce_v_real8(tempr_in,NSAMAX*2,3,tempr_out)
!         DO NSA=1,NSAMAX
!            RNDRS(NSA) = tempr_out(NSA)
!            RPDRS(NSA) = tempr_out(NSA+NSAMAX)
!         END DO
      ELSE
         DO NSA=1,NSAMAX
            CALL mtx_gatherv_real8(RNSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RNS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RJSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RJS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RWSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RWS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RWS123L(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RWS123(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RPCSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RPCS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RPWSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RPWS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RPESL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RPES(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RLHSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RLHS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RFWSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                                RFWS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RECSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RECS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RICSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RICS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RSPBL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RSPB(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RSPFL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RSPF(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RSPSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RSPS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RSPLL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RSPL(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RPDRL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RPDR(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RNDRL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RNDR(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
            CALL mtx_gatherv_real8(RTL_BULK(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                 RT_BULK(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
         END DO
      END IF

      DO NSA=1,NSAMAX
         DO NSB=1,NSBMAX
            CALL mtx_gatherv_real8(RPCS2L(NRSTART:NRENDX,NSB,NSA), &
                                   MTXLEN(NRANK+1), &
                                   RPCS2(1:NRMAX,NSB,NSA),NRMAX,MTXLEN,MTXPOS)
         ENDDO
      ENDDO

      END SUBROUTINE FPSAVECOMM

!^-------------------------------
    end MODULE fpsave




