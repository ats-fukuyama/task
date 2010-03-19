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

      interface
         DOUBLE PRECISION FUNCTION BESEKN(N,X)
           real(8) :: X
           integer :: N
         end function BESEKN
      end interface 

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
      real(8):: PV, WPL, WPM, WPP
      real(8):: DFP, DFT, FFP, testa, testb, FACT
      real(8),dimension(NSBMAX):: RSUM10
!      real(8),dimension(NRMAX, NSAMAX):: RWS123, RPCS, RPWS, RPES, RLHS
!      real(8),dimension(NRMAX, NSAMAX):: RFWS, RECS
!      real(8),dimension(NRMAX, NSBMAX, NSAMAX):: RPCS2
!      real(8),dimension(NSAMAX, 0:NTMAX):: PWT2, PTT2

      IF(ISAVE.NE.0) RETURN

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
            RSUM123=0.D0
            DO NSB=1,NSBMAX
               RSUM10(NSB)=0.D0
            ENDDO
            RSUM11B=0.D0
            RSUM11F=0.D0
            RSUM11S=0.D0
            RSUM11L=0.D0

            NSBA=NSB_NSA(NSA)
            IF(MODELA.eq.0)THEN
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA)
               END DO
            ENDDO
            ELSE
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA) &
                       *RLAMDAG(NTH,NR)
               END DO
            ENDDO
            END IF
!
            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=1,NPMAX
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA) &
                             *0.5D0*PM(NP,NSBA)**2
                     END DO
                  ENDDO
               ELSE
                  DO NP=1,NPMAX
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)/PV
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA) &
                             *(PV-1.D0)/THETA0(NSA)
                        RSUM123 = RSUM123                   &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA) &
                             *(PV-1.D0)/THETA0(NSA)
                     END DO
                  END DO
               ENDIF
            ELSE
               IF(MODELR.EQ.0) THEN
                  DO NP=1,NPMAX
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)*RLAMDAG(NTH,NR)
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA)  &
                             *0.5D0*PM(NP,NSBA)**2*RLAMDAG(NTH,NR)
                     END DO
                  ENDDO
               ELSE
                  DO NP=1,NPMAX
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)/PV*RLAMDAG(NTH,NR)
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA)  &
                             *(PV-1.D0)/THETA0(NSA)*RLAMDAG(NTH,NR)
                        RSUM123 = RSUM123                    &
                             +VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSBA)  &
                             *(PV-1.D0)/THETA0(NSA)*RLAMDAG(NTH,NR)
                     END DO
                  END DO
               ENDIF               
            END IF

!            open(8,file='dfp_dft_r6_a1_50_200.dat')
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
                       /DELP(NSBA)*(FNS(NTH,NP,NR,NSBA)-FNS(NTH,NP-1,NR,NSBA))
                  IF(NTH.EQ.1) THEN
                     DFT=1.D0/DELTH                             &
                         *(                                     &
                            ((1.D0-WPP)*FNS(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNS(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNS(NTH,NP  ,NR,NSBA)     &
                                  +WPM *FNS(NTH,NP-1,NR,NSBA))&
                          )

                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/DELTH                         & 
                         *(-                                    &
                            ((1.D0-WPM)*FNS(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNS(NTH-1,NP-1,NR,NSBA))&
                          +                                     &
                            ((1.D0-WPP)*FNS(NTH,NP  ,NR,NSBA)     &
                                  +WPP *FNS(NTH,NP-1,NR,NSBA))&
                          )
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)                  &
                         *(                                     &
                            ((1.D0-WPP)*FNS(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNS(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNS(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNS(NTH-1,NP-1,NR,NSBA))&
                                  )
                  ENDIF
                  FFP=    PG(NP,NSBA)                         &
                         *((1.D0-WPL)*FNS(NTH  ,NP  ,NR,NSBA)  &
                                +WPL *FNS(NTH  ,NP-1,NR,NSBA))
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
                  DO NSB=1,NSBMAX
                     RSUM10(NSB)=RSUM10(NSB)+PG(NP,NSBA)**2*SINM(NTH)/PV &
                         *(DCPP2(NTH,NP,NR,NSB,NSA)*DFP              &
                          +DCPT2(NTH,NP,NR,NSB,NSA)*DFT              & 
                          -FCPP2(NTH,NP,NR,NSB,NSA)*FFP)
                  END DO

                  RSUM11B=RSUM11B+PG(NP,NSBA)**4*SINM(NTH)/PV    &
                                  *SPPB(NTH,NP,NR,NSA)
                  RSUM11F=RSUM11F+PG(NP,NSBA)**4*SINM(NTH)/PV    &
                                  *SPPF(NTH,NP,NR,NSA)
                  RSUM11S=RSUM11S+PG(NP,NSBA)**4*SINM(NTH)/PV    &
                                  *SPPS(NTH,NP,NR,NSA)
                  RSUM11L=RSUM11L+PG(NP,NSBA)**4*SINM(NTH)/PV    &
                                  *PPL(NTH,NP,NR,NSA)*FNS(NTH,NP,NR,NSBA)

                  IF(NSA.eq.3)THEN
                  testa=PG(NP,NSBA)**2*SINM(NTH)/PV  &
                         *(DCPP(NTH,NP,NR,NSA)*DFP  &
                          +DCPT(NTH,NP,NR,NSA)*DFT  &
                          -FCPP(NTH,NP,NR,NSA)*FFP  &
                       )

                  testb=-PG(NP,NSBA)**2*SINM(NTH)/PV     &
                         *(FEPP(NTH,NP,NR,NSA)*FFP)

                  END IF
!                  IF(NR.eq.6.and.NSA.eq.1)THEN
!                     WRITE(8,'(3I4,1P12E12.4)') NR,NP,NTH             &
!                          ,PM(NP,NSA)*COSM(NTH),PM(NP,NSA)*SINM(NTH) &
!                          ,DFP,DFT,FFP,FNS(NTH,NP,NR,NSA)            &
!                          ,DCPP(NTH,NP,NR,NSA),DCPT(NTH,NP,NR,NSA)   &
!                          ,FCPP(NTH,NP,NR,NSA)
!                  END IF
               ENDDO
!            IF(NR.eq.6.and.NSA.eq.1)then
!               write(8,*)" "
!               write(8,*)" "
!            END IF

            ENDDO
!            close(8)

!            IF(NTG1.eq.NTMAX.and.NSA.eq.1.and.NR.eq.40) THEN
!               NP=20
!               DO NTH=1,NTHMAX
!                  WRITE(9,888) NR,NTH, PG(np,NSA), PM(NP,NSA) &
!                       ,DCPP(NTH,NP,NR,NSA),DCPT(NTH,NP,NR,NSA) &
!                       ,DCTT(NTH,NP,NR,NSA),FCPP(NTH,NP,NR,NSA)
!               END DO
!            END IF
               
 888        FORMAT(2I4,12E14.6)
            FACT=RNFP0(NSA)*1.D20
            RNSL(NR,NSA) = RSUM1*FACT                   *1.D-20
            RJSL(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)*1.D-6

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RWSL(NR,NSA) = RSUM3*FACT               *1.D-6
            RWS123L(NR,NSA) = RSUM123*FACT                   *1.D-6
            RPCSL(NR,NSA)=-RSUM4*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPWSL(NR,NSA)=-RSUM5*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6 
            RPESL(NR,NSA)=-RSUM6*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RLHSL(NR,NSA)=-RSUM7*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RFWSL(NR,NSA)=-RSUM8*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RECSL(NR,NSA)=-RSUM9*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6

            DO NSB=1,NSBMAX
               RPCS2L(NR,NSB,NSA)=-RSUM10(NSB) &
                                *FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            END DO
            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*0.5D0
            RSPBL(NR,NSA)= RSUM11B*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPFL(NR,NSA)= RSUM11F*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPSL(NR,NSA)= RSUM11S*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
            RSPLL(NR,NSA)= RSUM11L*FACT*2.D0*PI*DELP(NSBA)*DELTH*1.D-6
         ENDDO
      ENDDO

!      IF(NREND.eq.NRMAX)THEN
!      DO NSA=1,NSAMAX
!         NS=NS_NSA(NSA)
!         NSBA=NSB_NSA(NSA)
!         RSUM_FS2=0.D0
!         DO NP=1,NPMAX
!            DO NTH=1,NTHMAX
!               RSUM_FS2 = RSUM_FS2+VOLP(NTH,NP,NSBA)*FS2(NTH,NP,NSA) &
!                    *RLAMDAG(NTH,NRMAX+1)
!            END DO
!         ENDDO
!         FACT=RNFP0(NSA)*1.D20
!         RNS_S2(NSA) = RSUM_FS2*FACT                   *1.D-20
!         WRITE(*,'(2I2,1P3E14.6)') NRANK, NSA, RNS_S2(NSA)
!      END DO
!      END IF
      

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
         DO NSB=1,NSBMAX
            CALL mtx_gatherv_real8(RPCS2L(NRSTART:NRENDX,NSB,NSA), &
                                   MTXLEN(NRANK+1), &
                                   RPCS2(1:NRMAX,NSB,NSA),NRMAX,MTXLEN,MTXPOS)
         ENDDO
         CALL mtx_gatherv_real8(RSPBL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                                RSPB(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
         CALL mtx_gatherv_real8(RSPFL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                                RSPF(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
         CALL mtx_gatherv_real8(RSPSL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                                RSPS(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
         CALL mtx_gatherv_real8(RSPLL(NRSTART:NRENDX,NSA),MTXLEN(NRANK+1), &
                                RSPL(1:NRMAX,NSA),NRMAX,MTXLEN,MTXPOS)
      ENDDO

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
      real(8):: EAVE, rtemp, THETAL

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
         PWT2(NSA,NTG1)=0.D0
         PSPBT(NSA,NTG1)=0.D0
         PSPFT(NSA,NTG1)=0.D0
         PSPST(NSA,NTG1)=0.D0
         PSPLT(NSA,NTG1)=0.D0
         DO NSB=1,NSBMAX
            PPCT2(NSB,NSA,NTG1)= 0.D0
         END DO
      ENDDO

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            PNT(NSA,NTG1) =PNT(NSA,NTG1) +RNS(NR,NSA)*VOLR(NR)
            PIT(NSA,NTG1) =PIT(NSA,NTG1) +RJS(NR,NSA)*VOLR(NR)
            PWT(NSA,NTG1) =PWT(NSA,NTG1) +RWS(NR,NSA)*VOLR(NR)
            PPCT(NSA,NTG1)=PPCT(NSA,NTG1)+RPCS(NR,NSA)*VOLR(NR)
            PPWT(NSA,NTG1)=PPWT(NSA,NTG1)+RPWS(NR,NSA)*VOLR(NR)
            PPET(NSA,NTG1)=PPET(NSA,NTG1)+RPES(NR,NSA)*VOLR(NR)
            PLHT(NSA,NTG1)=PLHT(NSA,NTG1)+RLHS(NR,NSA)*VOLR(NR)
            PFWT(NSA,NTG1)=PFWT(NSA,NTG1)+RFWS(NR,NSA)*VOLR(NR)
            PECT(NSA,NTG1)=PECT(NSA,NTG1)+RECS(NR,NSA)*VOLR(NR)
            PSPBT(NSA,NTG1)=PSPBT(NSA,NTG1)+RSPB(NR,NSA)*VOLR(NR)
            PSPFT(NSA,NTG1)=PSPFT(NSA,NTG1)+RSPF(NR,NSA)*VOLR(NR)
            PSPST(NSA,NTG1)=PSPST(NSA,NTG1)+RSPS(NR,NSA)*VOLR(NR)
            PSPLT(NSA,NTG1)=PSPLT(NSA,NTG1)+RSPL(NR,NSA)*VOLR(NR)
            IF(MODELR.eq.1) then
               CALL FPNEWTON(NR,NSA,rtemp)
            else
               EAVE=RWS123(NR,NSA)*AMFP(NSA)*THETA0(NSA)     &
                    /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
               THETAL=2.d0*EAVE/3.d0
               rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
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
         ELSE
            PTT(NSA,NTG1)=0.D0
            PTT2(NSA,NTG1)=0.D0
         ENDIF
         PNT(NSA,NTG1) =PNT(NSA,NTG1)/TVOLR
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
      real(8):: RS

      NTG2=NTG2+1
      call fp_adjust_ntg2

      RTG(NTG2)=TIMEFP

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RNT(NR,NSA,NTG2) = RNS(NR,NSA)
            RJT(NR,NSA,NTG2) = RJS(NR,NSA)
            RWT(NR,NSA,NTG2) = RWS(NR,NSA)
            RPCT(NR,NSA,NTG2)= RPCS(NR,NSA)
            RPWT(NR,NSA,NTG2)= RPWS(NR,NSA)
            RPET(NR,NSA,NTG2)= RPES(NR,NSA)
            RLHT(NR,NSA,NTG2)= RLHS(NR,NSA)
            RFWT(NR,NSA,NTG2)= RFWS(NR,NSA)
            RECT(NR,NSA,NTG2)= RECS(NR,NSA)
            DO NSB=1,NSBMAX
               RPCT2(NR,NSB,NSA,NTG2)= RPCS2(NR,NSB,NSA)
            END DO
            RSPBT(NR,NSA,NTG2)= RSPB(NR,NSA)
            RSPFT(NR,NSA,NTG2)= RSPF(NR,NSA)
            RSPST(NR,NSA,NTG2)= RSPS(NR,NSA)
            RSPLT(NR,NSA,NTG2)= RSPL(NR,NSA)
!
            IF(RNS(NR,NSA).NE.0.D0) THEN
               RTT(NR,NSA,NTG2) = RWS(NR,NSA)*1.D6 &
                                  /(1.5D0*RNS(NR,NSA)*1.D20*AEE*1.D3)
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
      real(8):: rtotalPW, rtotalPC,rtotalSP,rtotalPC2
!      INCLUDE '../wr/wrcom1.inc'
!
      WRITE(6,*)"--------------------------------------------"
      WRITE(6,*)"-----Global data"
      WRITE(6,101) PTG(NTG1)*1000

      DO NSA=1,NSAMAX
         IF(MODELR.eq.0)THEN
            WRITE(6,102) NSA,NS_NSA(NSA), &
              PNT(NSA,NTG1),PTT(NSA,NTG1),PWT(NSA,NTG1),PIT(NSA,NTG1)
         ELSE
            WRITE(6,102) NSA,NS_NSA(NSA), &
              PNT(NSA,NTG1),PTT2(NSA,NTG1),PWT(NSA,NTG1),PIT(NSA,NTG1)
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
      DO NSA=1,NSAMAX
         WRITE(6,103) NSA,NS_NSA(NSA), &
              PPCT(NSA,NTG1),PPWT(NSA,NTG1),PPET(NSA,NTG1)
         rtotalPW=rtotalPW + PPWT(NSA,NTG1)
         rtotalPC=rtotalPC + PPCT(NSA,NTG1)
         rtotalSP=rtotalSP + PSPT(NSA,NTG1)
         rtotalPC2 = rtotalPC2 +PPCT(NSA,NTG1)-PPCT2(NSA,NSA,NTG1)
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
                                      PSPST(NSA,NTG1),PSPLT(NSA,NTG1)
      END DO
      write(*,105) rtotalpw
      write(*,107) rtotalPC
      write(*,109) rtotalSP
      write(*,110) rtotalPC2


 1001 FORMAT(I4,9E14.6)
      RETURN
  101 FORMAT(' TIME=',F12.3,' ms')
  102 FORMAT(' NSA,NS=',2I2,' n,T,W,I=',1PE11.4,1P3E12.4)
  103 FORMAT('        ',2I2,' PC,PW,PE=',10X,1P4E12.4)
  104 FORMAT('        ',2I2,' PCAB    =',10X,1P4E12.4)

 105  FORMAT('total absorption power [MW]', 1PE12.4)
 106  FORMAT(F12.4, 8E12.4)
 107  FORMAT('total collision power  [MW]', 1PE12.4)
 108  FORMAT('        ',2I2,' PSPB/F/S/L=',8X,1P4E12.4) 
 109  FORMAT('total source power     [MW]', 1PE12.4)
 110  FORMAT('collision balance      [MW]', 1PE12.4)
      END SUBROUTINE FPWRTGLB

! ***********************************************************

      SUBROUTINE FPWRTPRF
!
      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH
      real(8):: RTFDL, RTFD0L, THETAL, rtemp 
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
                    RM(NR),RNT(NR,NSA,NTG2),RTT(NR,NSA,NTG2), &
                    RJT(NR,NSA,NTG2),RPCT(NR,NSA,NTG2),       &
                    RPET(NR,NSA,NTG2),RPWT(NR,NSA,NTG2),      &
                    RSPBT(NR,NSA,NTG2),RSPFT(NR,NSA,NTG2)
!                    RPCT2(NR,NSA,NSA,NTG2),RPCT2(NR,2,NSA,NTG2)
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
  104 FORMAT(2I3,1P10E12.4) 
  105 FORMAT(' TIME=',F12.3,' ms'/                   &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j     ',5X,'PC     ',5X,'PE     ')
  106 FORMAT(' TIME=',F12.3,' ms'/                   &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j     ',5X,'PC     ',5X,'PE     ',5X,   &
           'PW     ',5X,'PNB//PLH',5X,'PNF//PIC',5X,'PEC   ')
      END SUBROUTINE FPWRTPRF
! ***********************************************************
!
      SUBROUTINE FPWRTSNAP
!
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NP, NTH
      real(8):: rnute, resist2, resist, dtaue, dtaui 
!      INCLUDE '../wr/wrcom1.inc'
!
      WRITE(6,101) TIMEFP*1000

!$$$      if(NTG1.ne.1)then
!$$$c-----check of conductivity--------
      IF(NTG1.ne.1)then
!         FACT1 =  &
!      2.D0*PI*RSRHON(RM(NRMAX))*(RSRHON(RG(2))-RSRHON(RG(NRMAX)))
!         FACT2 = 2.D0*PI*RR
!         rntv = 1.6D0/1.D19*1.5D0*
!     &FACT1*FACT2*(PTT(NTG1)-PTT(NTG1-1))*PNT(NTG1)*1.D3*1.D20/DELT

!      write(6,*) "Pw+Pc",(PPWT(NTG1)+PPCT(NTG1))*1.D6,"n*delta T*V",rntv
!      write(6,*) "c*epsilon*RabsE*S",RABSE**2*FACT1*8.8D0*VC/1.D12
!     & ,"Pw[W]",PPWT(NTG1)*1.D6
!      write(6,*) "Pc+PE+Pw",PPCT(NTG1)+PPET(NTG1)+PPWT(NTG1),
!     &     "(Pc+Pe+Pw)/Pc",(PPCT(NTG1)+PPET(NTG1)+PPWT(NTG1))/PPCT(NTG1)
!      write(6,*) "Pc+-Pc-/Pc",(PPCT(NTG1)-PPCT(NTG1-1))/PPCT(NTG1)
!     Spitzer

!      Do NSA=1,NSAMAX
!      Do NSB=1,NSBMAX
!         rnute=RNUD(1,1,NSA)*SQRT(2.D0)*RNFP(1,NSA)/RNFP0(NSA)     &
!              *(PTFP0(NSA)/PTFD(1,NSA))**2
!         resist2=RNUTE*AMFP(NSA)/RNFP(1,NSA)/AEFP(NSA)**2/1.D20
!      resist = 1.D0/(RNFP0(NSA)/RNFP(1,NSA)*1.D20*AEFP(NSA)**2/AMFP(NSA)&
!          /RNUD(1,1,NSA)*RNFD(1,NSA)/RNFP0(NSA) ) *SQRT(2.D0)/rnfp0(NSA)
!      if(E0.ne.0.d0) 
!     &    write(6,*) "J/E*eta*1.D6", PIT(NSA,NTG1)/E0*1.D6/FACT1*resist2
!     &     ,resist,resist2
!     &     ,RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)
!     & ,"THETA0", (PTFP0(NSA)/(AMFP(NSA)*VC))**2
!     &     ,DCTT2(2,2,1,NSB,NSA),NSA,NSB
!      END DO
!      END DO

      dtaue=1.09D16*(PTT(1,NTG1))**(1.5D0)/PNT(1,NTG1)/PZ(2)**2 &
           /15.D0*1.D-20
      if(nsamax.gt.1) then
      dtaui=6.6D17*(PTT(2,NTG1))**(1.5D0)/PNT(2,NTG1)/PZ(2)**4  &
           /15.D0*1.D-20*(AMFP(2)/AMFP(2))**0.5D0
      else
         dtaui=0.d0
      endif
!      write(6,*)"tau_e tau_i[ms]",dtaue*1.D3,dtaui*1.D3
      end if
!----end of conductivity check---------

      RETURN
!
  101 FORMAT(1H ,' TIME=',F12.3,' ms')

      END SUBROUTINE FPWRTSNAP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine FPNEWTON(NR,NSA,rtemp)

      IMPLICIT NONE
      integer:: NSA, NSB, NR, NP, NTH, ncount
      real(8):: rtemp, xeave, xtemp, thetal, EAVE

!-----Average kinetic energy
      EAVE=RWS123(NR,NSA)*AMFP(NSA)*THETA0(NSA) &
           /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)

!-----initial value of THETAL
      THETAL=2.d0*EAVE/3.d0
      xtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)

      CALL XNEWTON(EAVE,THETAL,ncount)

      rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
      xeave=AMFP(NSA)*VC**2*EAVE/(AEE*1.D3)
!      write(6,'(3I5,1P3E12.4)') NSA,NR,ncount,xeave,xtemp,rtemp

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

!^-------------------------------
    end MODULE fpsave




