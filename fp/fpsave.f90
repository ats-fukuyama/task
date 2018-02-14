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
      integer:: NR, NSA, NS

      IF(ISAVE.NE.0) RETURN

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

      CALL Define_Bulk_NP

      CALL MOMENT_0TH_ORDER(FNSP,RNSL)
      CALL MOMENT_1ST_ORDER(FNSP,RJSL)
      CALL MOMENT_2ND_ORDER(FNSP,RWSL)
      CALL POWER_FROM_DIFFUSION_COEF
      CALL POWER_FROM_SOURCE_TERM
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
!-------    Calculation of bulk temperature
            CALL BULK_TEMPERATURE(NP_BULK(NR,NSA),NR,NSA) 
            CALL mtx_set_communicator(comm_np) 
!-------    Calculation of bulk temperature
         ENDDO ! NSA
      ENDDO ! NR
      CALL POWER_FROM_RADIAL_TRANSPORT

      CALL mtx_reset_communicator 

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
      integer:: NSA, NSB, NR
      real(8):: rtemp, rtemp2

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
         PWRT(NSA,NTG1)=0.D0
         PWMT(NSA,NTG1)=0.D0
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
         PPST(NSA,NTG1)=0.D0
         PPLT(NSA,NTG1)=0.D0
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
            PWRT(NSA,NTG1)=PWRT(NSA,NTG1)+RWRS(NR,NSA)*VOLR(NR)
            PWMT(NSA,NTG1)=PWMT(NSA,NTG1)+RWMS(NR,NSA)*VOLR(NR)
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
         ELSE
            PTT(NSA,NTG1)=0.D0
            PTT2(NSA,NTG1)=0.D0
         ENDIF
         PTT_BULK(NSA,NTG1) = PTT_BULK(NSA,NTG1)/PNT(NSA,NTG1)
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
      integer:: NR, NSA, NSB
      real(8):: RS, rtemp

      NTG2=NTG2+1
      call fp_adjust_ntg2

      RTG(NTG2)=TIMEFP

      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RNT(NR,NSA,NTG2) = RNS(NR,NSA)
            IF(MODEL_DISRUPT.eq.1) rate_runaway(NR,NTG2)=RFP(NR)
            RJT(NR,NSA,NTG2) = RJS(NR,NSA)
            RWT(NR,NSA,NTG2) = RWS(NR,NSA)
            RPCT(NR,NSA,NTG2)= RPCS(NR,NSA)
            RPWT(NR,NSA,NTG2)= RPWS(NR,NSA)
            RPET(NR,NSA,NTG2)= RPES(NR,NSA)
            RLHT(NR,NSA,NTG2)= RLHS(NR,NSA)
            RFWT(NR,NSA,NTG2)= RFWS(NR,NSA)
            RECT(NR,NSA,NTG2)= RECS(NR,NSA)
            RWRT(NR,NSA,NTG2)= RWRS(NR,NSA)
            RWMT(NR,NSA,NTG2)= RWMS(NR,NSA)
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
!            RT_T(NR,NSA) = RTT(NR,NSA,NTG2)
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
      integer:: NSA, NSB
      real(8):: rtotalPW, rtotalPC,rtotalSP,rtotalPC2
      real(8):: rtotalDR,rtotalLH,rtotalFW,rtotalEC,rtotalWR,rtotalWM,rtotalIP
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
      rtotalLH=0.D0
      rtotalFW=0.D0
      rtotalEC=0.D0
      rtotalWR=0.D0
      rtotalWM=0.D0
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
         rtotalLH=rtotalLH + PLHT(NSA,NTG1)
         rtotalFW=rtotalFW + PFWT(NSA,NTG1)
         rtotalEC=rtotalEC + PECT(NSA,NTG1)
         rtotalWR=rtotalWR + PWRT(NSA,NTG1)
         rtotalWM=rtotalWM + PWMT(NSA,NTG1)
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
      write(6,105) rtotalpw,rtotalWR,rtotalWM
      write(6,115) rtotalLH,rtotalFW,rtotalEC
      write(6,107) rtotalPC
      write(6,109) rtotalSP
      write(6,110) rtotalPC2
      WRITE(6,'("total plasma current   [MA]",1PE12.4)') rtotalIP

      RETURN

  101 FORMAT(' TIME=',F12.3,' ms')
  112 FORMAT(' NSA,NS=',2I2,' n,T,W,I,dn,T2=',1PE11.4,1P6E12.4)
  104 FORMAT('        ',2I2,' PCAB    =',10X,1P14E12.4)
  105 FORMAT('Total absorption power [MW]', 1PE12.4,'    WR:',1PE12.4,'    WM:',1PE12.4)
  115 FORMAT('   absorption power [MW] LH', 1PE12.4,'    FW:',1PE12.4,'    EC:',1PE12.4)
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
      integer:: NSA, NR, NS
      real(8):: rtemp
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
         NS=NS_NSA(NSA)
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
               WRITE(6,fmt0) NSA,NS_NSA(NSA),&
                    RM(NR),RNT(NR,NSA,NTG2),RTT(NR,NSA,NTG2), &
                    RJT(NR,NSA,NTG2),RPCT(NR,NSA,NTG2),       &
                    RPET(NR,NSA,NTG2), &
!                    RPWT(NR,NSA,NTG2),&
                    RNS_DELF(NR,NS), &
                    RSPBT(NR,NSA,NTG2), &
                    RSPS_CX(NR,NSA), &
                    RT_BULK(NR,NSA)!, &
            END IF
         ENDDO
      ENDDO
      RETURN
  106 FORMAT( &
           'NSA/NS',5X,'RM',10X,' n',8X,' T    ',6X, &
           ' j    ',5X,'PC     ',5X,'PE     ',5X,   &
           'n_b',9X,'PNBI',8X,'PDRP',8X,'TBULK' )
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
      integer:: NR, NSA, NSB, NS
      real(8):: rnute, resist
      real(8):: taue_col, sigma_sp, FACT
      real(8):: taue_col2, sigma_sp2

!-----check of conductivity--------
      IF(NTG1.ne.1.and.NRANK.eq.0)then
         NR=1
         NSA=1
         NSB=2
         NS=NS_NSA(NSA)
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
                    " THETA0= ", THETA0(NS)
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
      integer:: ncount, NS
      real(8):: xeave
      real(8):: xtemp, thetal, EAVE

      NS=NS_NSA(NSA)
!-----Average kinetic energy
!      EAVE=RWS(NR,NSA)*AMFP(NSA)*THETA0(NS) &
!           /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
      EAVE=RWSL_*AMFP(NSA)*THETA0(NS) & 
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
      integer:: NR, NSA, NSB, NP, NSW, N
      double precision,dimension(NRSTART:NREND,NSAMAX):: work
      double precision,dimension(NRMAX,NSAMAX):: workg
      double precision,dimension(NSAMAX):: temp_nsanr
      double precision,dimension(NPMAX,NRSTART:NREND):: temp_npnr1
      double precision,dimension(NPMAX,NRMAX):: temp_npnr2
      double precision,dimension(NPMAX,NRMAX,NSASTART:NSAEND):: temp_npnr3
      integer,dimension(NSAMAX):: vloc

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL,SAVLEN(NRANK+1),RNS,N,NSA)
         CALL fp_gatherv_real8_sav(RJSL,SAVLEN(NRANK+1),RJS,N,NSA)
         CALL fp_gatherv_real8_sav(RWSL,SAVLEN(NRANK+1),RWS,N,NSA)
         CALL fp_gatherv_real8_sav(RWS123L,SAVLEN(NRANK+1),RWS123,N,NSA)
         CALL fp_gatherv_real8_sav(RPCSL,SAVLEN(NRANK+1),RPCS,N,NSA)
         CALL fp_gatherv_real8_sav(RPWSL,SAVLEN(NRANK+1),RPWS,N,NSA)
         CALL fp_gatherv_real8_sav(RPESL,SAVLEN(NRANK+1),RPES,N,NSA)
         CALL fp_gatherv_real8_sav(RLHSL,SAVLEN(NRANK+1),RLHS,N,NSA)
         CALL fp_gatherv_real8_sav(RFWSL,SAVLEN(NRANK+1),RFWS,N,NSA)
         CALL fp_gatherv_real8_sav(RECSL,SAVLEN(NRANK+1),RECS,N,NSA)
         CALL fp_gatherv_real8_sav(RWRSL,SAVLEN(NRANK+1),RWRS,N,NSA)
         CALL fp_gatherv_real8_sav(RWMSL,SAVLEN(NRANK+1),RWMS,N,NSA)
         CALL fp_gatherv_real8_sav(RSPBL,SAVLEN(NRANK+1),RSPB,N,NSA)
         CALL fp_gatherv_real8_sav(RSPFL,SAVLEN(NRANK+1),RSPF,N,NSA)
         CALL fp_gatherv_real8_sav(RSPSL,SAVLEN(NRANK+1),RSPS,N,NSA)
         CALL fp_gatherv_real8_sav(RSPLL,SAVLEN(NRANK+1),RSPL,N,NSA)
         CALL fp_gatherv_real8_sav(RPDRL,SAVLEN(NRANK+1),RPDR,N,NSA)
         CALL fp_gatherv_real8_sav(RNDRL,SAVLEN(NRANK+1),RNDR,N,NSA)
         CALL fp_gatherv_real8_sav(RTL_BULK,SAVLEN(NRANK+1),RT_BULK,N,NSA) 
!         CALL fp_gatherv_real8_sav(RDIDTL,SAVLEN(NRANK+1),RDIDT,N,NSA)
         CALL fp_gatherv_real8_sav(RPSSL,SAVLEN(NRANK+1),RPSS,N,NSA)
         CALL fp_gatherv_real8_sav(RPLSL,SAVLEN(NRANK+1),RPLS,N,NSA)
         CALL fp_gatherv_real8_sav(RSPSL_CX,SAVLEN(NRANK+1),RSPS_CX,N,NSA)
      END DO
      work(:,:)=0.D0
      workg(:,:)=0.D0
      DO NSB=1,NSBMAX
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               work(NR,NSA) = RPCS2L(NR,NSB,NSA)
!               WRITE(*,'(A,3I5,2E14.6)') "TEST WORK ", NPSTART, NSB, NSA, work(NR,NSA), RPCS2L(NR,NSB,NSA)
            END DO
         END DO
         DO N=1,NSW
            NSA=N+NSASTART-1
            CALL fp_gatherv_real8_sav(work,SAVLEN(NRANK+1),workg,N,NSA)
         ENDDO

         DO NSA=1,NSAMAX
            DO NR=1,NRMAX
               RPCS2(NR,NSB,NSA) = workg(NR,NSA)
!               WRITE(*,'(A,3I5,2E14.6)') "TEST WORKG ", NPSTART, NSB, NSA, workg(NR,NSA), RPCS2(NR,NSB,NSA)
            END DO
         END DO
      ENDDO

      work(:,:)=0.D0
      workg(:,:)=0.D0
      DO NSB=1,NSBMAX
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               work(NR,NSA) = RPCS2L_DEL(NR,NSB,NSA)
!               WRITE(*,'(A,3I5,2E14.6)') "TEST WORK ", NPSTART, NSB, NSA, work(NR,NSA), RPCS2L(NR,NSB,NSA)
            END DO
         END DO
         DO N=1,NSW
            NSA=N+NSASTART-1
            CALL fp_gatherv_real8_sav(work,SAVLEN(NRANK+1),workg,N,NSA)
         ENDDO

         DO NSA=1,NSAMAX
            DO NR=1,NRMAX
               RPCS2_DEL(NR,NSB,NSA) = workg(NR,NSA)
!               WRITE(*,'(A,3I5,2E14.6)') "TEST WORKG ", NPSTART, NSB, NSA, workg(NR,NSA), RPCS2(NR,NSB,NSA)
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

      END SUBROUTINE FPSAVECOMM
!^------------------------------ 
!==============================================================
!     update RNS_DELF, RWS_PARA, RWS_PERP
      SUBROUTINE COUNT_BEAM_DENSITY

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NTH, NP, NR, NSA, NS, N, NSW
      double precision:: FACT, RSUM3_PARA, RSUM3_PERP

      RNS_DELF_NSA(:,:)=0.D0
      RWS_DELF_PARA(:,:)=0.D0
      RWS_DELF_PERP(:,:)=0.D0
      RNSL_DELF(:,:)=0.D0
      RWSL_PARA(:,:)=0.D0
      RWSL_PERP(:,:)=0.D0

      CALL MOMENT_0TH_ORDER(FNSP_DEL,RNSL_DELF)
      CALL mtx_set_communicator(comm_np) 

      DO NSA=NSASTART, NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTART, NREND
            RSUM3_PARA=0.D0
            RSUM3_PERP=0.D0
!Pressure
            IF(MODELA.eq.0) THEN
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM3_PARA = RSUM3_PARA                   &
                          +VOLP(NTH,NP,NS)*FNSP_DEL(NTH,NP,NR,NSA) &
                          *(PM(NP,NS)*COSM(NTH))**2
                     RSUM3_PERP = RSUM3_PERP                   &
                          +VOLP(NTH,NP,NS)*FNSP_DEL(NTH,NP,NR,NSA) &
                          *0.5D0*(PM(NP,NS)*SINM(NTH))**2
                  END DO
               ENDDO
            ELSE ! MODELA=1
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM3_PARA = RSUM3_PARA                   &
                          +VOLP(NTH,NP,NS)*FNSP_DEL(NTH,NP,NR,NSA) &
                          *(PM(NP,NS)*COSM(NTH))**2*RLAMDA(NTH,NR)*RFSADG(NR)
                     RSUM3_PERP = RSUM3_PERP                   &
                          +VOLP(NTH,NP,NS)*FNSP_DEL(NTH,NP,NR,NSA) &
                          *0.5D0*(PM(NP,NS)*SINM(NTH))**2*RLAMDA(NTH,NR)*RFSADG(NR)
                  END DO
               ENDDO
            END IF

            CALL p_theta_integration(RSUM3_PARA)
            CALL p_theta_integration(RSUM3_PERP)
            
            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RWSL_PARA(NR,NSA) = RSUM3_PARA*FACT
            RWSL_PERP(NR,NSA) = RSUM3_PERP*FACT

         END DO
      END DO
      CALL mtx_reset_communicator 

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL_DELF,SAVLEN(NRANK+1),RNS_DELF_NSA,N,NSA)
         CALL fp_gatherv_real8_sav(RWSL_PARA,SAVLEN(NRANK+1),RWS_DELF_PARA,N,NSA)
         CALL fp_gatherv_real8_sav(RWSL_PERP,SAVLEN(NRANK+1),RWS_DELF_PERP,N,NSA)
      END DO
      CALL mtx_reset_communicator 

      RNS_DELF(:,:)=0.D0
      DO NSA=1, NSAMAX
         NS=NS_NSA(NSA)
         DO NR=1, NRMAX
            RNS_DELF(NR,NS)=RNS_DELF_NSA(NR,NSA)
            RWS_PARA(NR,NS)=RWS_DELF_PARA(NR,NSA)
            RWS_PERP(NR,NS)=RWS_DELF_PERP(NR,NSA)
         END DO
      END DO

      CALL mtx_broadcast_real8(RNS_DELF,NRMAX*NSMAX)
      CALL mtx_broadcast_real8(RWS_PARA,NRMAX*NSMAX)
      CALL mtx_broadcast_real8(RWS_PERP,NRMAX*NSMAX)

      END SUBROUTINE COUNT_BEAM_DENSITY
!==============================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FPSAVECOMM2

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NSA, NSW, N

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
!      DO N=1,NRMAX
!         WRITE(6,'(A,I5,1pE12.4)') 'NR,RTL_BULK=', &
!              N,RTL_BULK(N,1)
!      END DO
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RNSL,SAVLEN(NRANK+1),RNS,N,NSA)
         CALL fp_gatherv_real8_sav(RWSL,SAVLEN(NRANK+1),RWS,N,NSA)
         CALL fp_gatherv_real8_sav(RTL_BULK,SAVLEN(NRANK+1),RT_BULK,N,NSA) 
      END DO
!      DO N=1,NRMAX
!         WRITE(6,'(A,I5,1p2E12.4)') 'NR,RTL_BULK,RT_BULK=', &
!              N,RTL_BULK(N,1),RT_BULK(N,1)
!      END DO
      CALL mtx_reset_communicator 


      CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RWS,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RT_BULK,NRMAX*NSAMAX)

      END SUBROUTINE FPSAVECOMM2
!==============================================================
      SUBROUTINE update_RN_RT(recv)

      USE libmtx
      USE fpmpi
      USE EG_READ
      USE plprof
      IMPLICIT NONE
      integer:: NR, NSA, NS
      real(8):: rhon
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND),intent(in):: recv

      CALL Define_Bulk_NP
      CALL MOMENT_0TH_ORDER(recv,RNSL)
      CALL MOMENT_2ND_ORDER(recv,RWSL)
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
!-------    Calculation of bulk temperature
            CALL BULK_TEMPERATURE(NP_BULK(NR,NSA),NR,NSA) 
            CALL mtx_set_communicator(comm_np) 
!-------    Calculation of bulk temperature
         ENDDO ! NSA
      ENDDO ! NR
      CALL mtx_reset_communicator
      CALL FPSAVECOMM2

      IF(MODEL_EX_READ_Tn.eq.0)THEN
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               RHON=RM(NR) 
               CALL PL_PROF(RHON,PLF) 
               RN_TEMP(NR,NS)=PLF(NS)%RN
               RT_TEMP(NR,NS)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
!               IF(NRANK.EQ.0) THEN
!                  IF(RT_TEMP(NR,NS).LE.0.D0) &
!                       WRITE(6,'(A,2I5,1P3E12.4)') &
!                       'NS,NR,RTPR,RTPP=',NS,NR,PLF(NS)%RTPR,PLF(NS)%RTPP
!               END IF
            ENDDO
         END DO
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1,NRMAX
               RN_TEMP(NR,NS) = RNS(NR,NSA)
               RT_TEMP(NR,NS) = RT_BULK(NR,NSA)
!               IF(NRANK.EQ.0) THEN
!                  IF(RT_TEMP(NR,NS).LE.0.D0) &
!                       WRITE(6,'(A,2I5,1P4E12.4)') &
!                       'NSA,NR,RT_BULK=',NSA,NR,RNS(NR,NSA),RT_BULK(NR,NSA)
!               END IF
            END DO
         ENDDO
      ELSEIF(MODEL_EX_READ_Tn.ne.0)THEN
         CALL MAKE_EXP_PROF(timefp+delt)
      END IF

      END SUBROUTINE update_RN_RT
!==============================================================
      SUBROUTINE BULK_TEMPERATURE(NPB,NR,NSA)

      USE fpmpi
      IMPLICIT NONE
      INTEGER,intent(in):: NPB, NR, NSA
!      REAL(8),INTENT(OUT):: RTL
      integer:: ISW_BULK, NP, NTH, NS
      real(8):: RSUM_T, RSUM_V, PV, DFDP, WPL, FFP, RSUMN, RSUMW, FACT
      real(8),dimension(NTHMAX,NPMAX):: T_BULK
      real(8),dimension(NPSTART:NPEND):: RPL_BULK_send 
      real(8),dimension(NPMAX):: RPL_BULK_recv
      real(8):: RNL_BULK, RWL_BULK, rtemp

!      ISW_BULK=0 ! sometimes DFDP becomes 0 and then it makes density NaN. (for FACT_BULK < 4)
      ISW_BULK=1 ! requires higher FACT_BULK (FACT_BULK >= 4) to obtain accurate RT_bulk

      NS=NS_NSA(NSA)
!      RTL_BULK(:,:)=0.D0
!      RPL_BULK(:,:,:)=0.D0

      CALL mtx_set_communicator(comm_np)
      IF(MODEL_EX_READ_Tn.eq.0)THEN
      IF(ISW_BULK.eq.0)THEN ! BULK T calculation using dfdp
         RSUM_T=0.D0
         RSUM_V=0.D0
         DO NP=NPSTART,NPEND
            IF(NP.ge.2.and.NP.le.NPB)THEN
               PV=SQRT(1.D0+THETA0(NS)*PG(NP,NS)**2)
               DO NTH=1,NTHMAX
                  IF(FNSP(NTH,NP,NR,NSA).gt.0.D0.and.FNSP(NTH,NP-1,NR,NSA).gt.0.D0)THEN
                     DFDP=DELP(NS)/ &
                          ( log(FNSP(NTH,NP,NR,NSA))-log(FNSP(NTH,NP-1,NR,NSA)) )
                  ELSE
                     WPL=WEIGHP(NTH  ,NP,NR,NSA)
                     FFP=   ( (1.D0-WPL)*FNSP(NTH  ,NP  ,NR,NSA)  &
                          +WPL *FNSP(NTH  ,NP-1,NR,NSA) )
                     DFDP=DELP(NS)*FFP/(                         &
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
                  
                  T_BULK(NTH,NP)=-PG(NP,NS)*PTFP0(NSA)*DFDP &
                       /AEE/1.D3*VTFP0(NSA)/PV 
                  RSUM_T = RSUM_T + T_BULK(NTH,NP)*VOLP(NTH,NP,NS)
                  RSUM_V = RSUM_V + VOLP(NTH,NP,NS)
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
                     RSUMN = RSUMN+VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA)
                  END DO
               END IF
            ENDDO
         ELSE
            DO NP=NPSTART,NPEND
               IF(NP.le.NPB)THEN
                  DO NTH=1,NTHMAX
                     RSUMN = RSUMN+VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA) &
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
                             +VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA) &
                             *0.5D0*PM(NP,NS)**2
                     END DO
                  END IF
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  IF(NP.le.NPB)THEN
                     PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                     DO NTH=1,NTHMAX
                        RSUMW = RSUMW                       &
                             +VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA) &
                             *(PV-1.D0)/THETA0(NS)
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
                             +VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA)  &
                             *0.5D0*PM(NP,NS)**2*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                  END IF
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  IF(NP.le.NPB)THEN
                     PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                     DO NTH=1,NTHMAX
                        RSUMW = RSUMW                        &
                             +VOLP(NTH,NP,NS)*FNSP(NTH,NP,NR,NSA)  &
                             *(PV-1.D0)/THETA0(NS)*RLAMDA(NTH,NR)*RFSADG(NR)
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
!         write(6,'(A,4I5,1P3E12.4)') &
!              'NR,NSA,A,R,RNL,RWL,RTL_BULK', &
!              NR,NSA,MODELA,MODELR,RNL_BULK,RWL_BULK,RTL_BULK(NR,NSA)
      END IF
      ELSE ! MODEL_EX_READ_Tn!=0
         RTL_BULK(NR,NSA)=RT_READ(NR,NS)
      END IF
      CALL mtx_reset_communicator       

      END SUBROUTINE BULK_TEMPERATURE
!==============================================================
      SUBROUTINE Define_Bulk_NP

      IMPLICIT NONE
      integer:: NP, NR, NSA, NS
      double precision:: pmax_bulk, p_bulk_r, rhon, RTFPL
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF

!     DEFINE BULK MOMENTUM RANGE: 0 < NP < NP_BULK(NR,NSA) 
      IF(MODEL_BULK_CONST.eq.0)THEN
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1, NRMAX
               IF(NT_init.eq.0)THEN
                  PMAX_BULK = FACT_BULK
                  P_BULK_R = PMAX_BULK* RT_TEMP(NR,NS)/RTFP0(NSA)
               ELSE
                  RHON=RM(NR)
                  CALL PL_PROF(RHON,PLF)
                  RTFPL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
                  PMAX_BULK = SQRT( RT_BULK(NR,NSA)/RTFPL )*FACT_BULK
                  P_BULK_R = PMAX_BULK* RT_BULK(NR,NSA)/RTFP0(NSA)
               END IF
               DO NP=NPMAX, 1, -1
                  IF(P_BULK_R.lt.PG(NP,NS)) NP_BULK(NR,NSA)=NP
               END DO
               IF(PMAX(NS).lt.P_BULK_R)THEN
                  NP_BULK(NR,NSA) = NPMAX
               END IF
            END DO
         END DO
      ELSE
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1, NRMAX
               PMAX_BULK = FACT_BULK
               IF(MODEL_EX_READ_Tn.eq.0)THEN
                  P_BULK_R = PMAX_BULK* RTFP(NR,NSA)/RTFP0(NSA)
               ELSE
                  P_BULK_R = PMAX_BULK* RT_TEMP(NR,NS)/RTFP0(NSA)
               END IF
               DO NP=NPMAX, 1, -1
                  IF(P_BULK_R.lt.PG(NP,NS)) NP_BULK(NR,NSA)=NP
               END DO
               IF(PMAX(NS).lt.P_BULK_R)THEN
                  NP_BULK(NR,NSA) = NPMAX
               END IF
            END DO
         END DO
      END IF
!      IF(NRANK.eq.0) WRITE(*,'(A,E14.6,2I5,5E14.6)') &
!           "NR, time, NP_BULK ", timefp, (NP_BULK(1,i),i=1,NSAMAX)

      END SUBROUTINE Define_Bulk_NP
!==============================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MOMENT_0TH_ORDER(SEND,RECV)
!     SEND=FNSP, FNSP_DEL,    RECV=RNSL, RNSL_DEL

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      INTEGER:: NP, NTH, NR, NSA, NS
      double precision:: RSUM1, FACT
      double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND),intent(in)::SEND
      double precision,dimension(NRSTART:NREND,NSAMAX),intent(out):: RECV

      RSUM1=0.D0
      CALL mtx_set_communicator(comm_np)

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTART,NREND
            IF(MODELA.eq.0)THEN
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA)
                  END DO
               ENDDO
            ELSE
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA) &
                          *RLAMDA(NTH,NR)*RFSADG(NR)
                  END DO
               ENDDO
            END IF

            CALL p_theta_integration(RSUM1)

            FACT=RNFP0(NSA)*1.D20
            RECV(NR,NSA) = RSUM1*FACT*1.D-20
!            RNSL(NR,NSA) = RSUM1*FACT*1.D-20
         END DO
      END DO

      CALL mtx_reset_communicator

      END SUBROUTINE MOMENT_0TH_ORDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MOMENT_1ST_ORDER(SEND,RECV)

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      INTEGER:: NP, NTH, NR, NSA, NS
      double precision:: FACT, RSUM2, PV
      double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND),intent(in)::SEND
      double precision,dimension(NRSTART:NREND,NSAMAX),intent(out):: RECV

      CALL mtx_set_communicator(comm_np) 

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTART,NREND
            RSUM2=0.D0

            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA) &
                             *PM(NP,NS)*COSM(NTH)
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA) &
                             *PM(NP,NS)*COSM(NTH)/PV
                     END DO
                  END DO
               ENDIF
            ELSE ! MODELA=1
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA)  &
                             *PM(NP,NS)*COSM(NTH)*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                     DO NTH=1, ITL(NR)-1
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA)  &
                             *PM(NP,NS)*COSM(NTH)/PV &
                             *(1.D0+EPSRM2(NR))
                     END DO
                     DO NTH=ITU(NR)+1, NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA)  &
                             *PM(NP,NS)*COSM(NTH)/PV &
                             *(1.D0+EPSRM2(NR))
                     END DO
                  END DO
               ENDIF           
            END IF

            CALL p_theta_integration(RSUM2)

            FACT=RNFP0(NSA)*1.D20
            RECV(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)*1.D-6
         END DO
      END DO

      CALL mtx_reset_communicator 
      END SUBROUTINE MOMENT_1ST_ORDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE MOMENT_2ND_ORDER(SEND,RECV)

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      INTEGER:: NP, NTH, NR, NSA, NS
      double precision:: FACT, RSUM3, PV
      double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND),intent(in)::SEND
      double precision,dimension(NRSTART:NREND,NSAMAX),intent(out):: RECV

      CALL mtx_set_communicator(comm_np) 

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTART,NREND
            RSUM3=0.D0

            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA) &
                             *0.5D0*PM(NP,NS)**2
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                       &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA) &
                             *(PV-1.D0)/THETA0(NS)
                     END DO
                  END DO
               ENDIF
            ELSE ! MODELA=1
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA)  &
                             *0.5D0*PM(NP,NS)**2*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                     DO NTH=1,NTHMAX
                        RSUM3 = RSUM3                        &
                             +VOLP(NTH,NP,NS)*SEND(NTH,NP,NR,NSA)  &
                             *(PV-1.D0)/THETA0(NS)*RLAMDA(NTH,NR)*RFSADG(NR)
                     END DO
                  END DO
               ENDIF               
            END IF

            CALL p_theta_integration(RSUM3)

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RECV(NR,NSA) = RSUM3*FACT               *1.D-6
         END DO
      END DO

      CALL mtx_reset_communicator 
      END SUBROUTINE MOMENT_2ND_ORDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE POWER_FROM_DIFFUSION_COEF

      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      INTEGER:: NP, NTH, NR, NSA, NS, NPS, NSB
      double precision:: FACT, DFP, DFT, FFP
      double precision,dimension(NTHMAX,NPSTART:NPENDWG,NRSTART:NRENDWM,NSAMAX):: no_coef
      double precision:: RSUM1, RSUM2, RSUM3, RSUM4, RSUM5, RSUM6, RSUM8, RSUM9, RSUM10
      double precision,dimension(NSBMAX):: RSUM11, RSUM12
      double precision:: RSUM_WR, RSUM_WM

      no_coef(:,:,:,:)=0.D0
      CALL mtx_set_communicator(comm_np)

      IF(NPSTART.eq.1)THEN
         NPS=2
      ELSE
         NPS=NPSTART
      END IF
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTART,NREND
            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            RSUM5=0.D0
            RSUM6=0.D0
            RSUM8=0.D0
            RSUM9=0.D0
            RSUM10=0.D0
            RSUM11(:)=0.D0
            RSUM12(:)=0.D0
            RSUM_WR=0.D0
            RSUM_WM=0.D0
            DO NP=NPS,NPEND
               DO NTH=1,NTHMAX
                  CALL PHASE_SPACE_DERIVATIVE_F(NTH,NP,NR,NSA,FNSP,DFP,DFT,FFP)

                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DCPP,DCPT,FCPP,RSUM1)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DWPP,DWPT,no_coef,RSUM2)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,no_coef,no_coef,FEPP,RSUM3)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DWLHPP,DWLHPT,no_coef,RSUM4)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DWFWPP,DWFWPT,no_coef,RSUM5)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DWECPP,DWECPT,no_coef,RSUM6)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,no_coef,no_coef,FSPP,RSUM8)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DLPP,no_coef,FLPP,RSUM9)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DPP,DPT,FPP,RSUM10)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DWWRPP,DWWRPT,no_coef,RSUM_WR)
                  CALL INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DWWMPP,DWWMPT,no_coef,RSUM_WM)
                  DO NSB=1,NSBMAX
                     CALL INTEGRAND_IN_M2OD_NSB(NTH,NP,NR,NSA,NSB,DFP,DFT,FFP,DCPP2,DCPT2,FCPP2,RSUM11)
                  END DO
               END DO
            END DO
            CALL p_theta_integration(RSUM1)
            CALL p_theta_integration(RSUM2)
            CALL p_theta_integration(RSUM3)
            CALL p_theta_integration(RSUM4)
            CALL p_theta_integration(RSUM5)
            CALL p_theta_integration(RSUM6)
            CALL p_theta_integration(RSUM8)
            CALL p_theta_integration(RSUM9)
            CALL p_theta_integration(RSUM10)
            CALL p_theta_integration(RSUM_WR)
            CALL p_theta_integration(RSUM_WM)
            DO NSB=1,NSBMAX
               CALL p_theta_integration(RSUM11(NSB))
            END DO

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*RFSADG(NR)

            RPCSL(NR,NSA)=-RSUM1*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RPWSL(NR,NSA)=-RSUM2*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RPESL(NR,NSA)=-RSUM3*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RLHSL(NR,NSA)=-RSUM4*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RFWSL(NR,NSA)=-RSUM5*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RECSL(NR,NSA)=-RSUM6*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RPSSL(NR,NSA)=-RSUM8*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RPLSL(NR,NSA)=-RSUM9*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RWS123L(NR,NSA) =-RSUM10*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RWRSL(NR,NSA)=-RSUM_WR*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            RWMSL(NR,NSA)=-RSUM_WM*FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            DO NSB=1,NSBMAX
               RPCS2L(NR,NSB,NSA)=-RSUM11(NSB) &
                    *FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
            END DO

!      for delta f
            IF(MODEL_DELTA_F(NS).eq.1)THEN
!      Collisional power transfer delta f
               DO NP=NPS,NPEND
                  DO NTH=1,NTHMAX
                     CALL PHASE_SPACE_DERIVATIVE_F(NTH,NP,NR,NSA,FNSP_DEL,DFP,DFT,FFP)
                     DO NSB=1,NSBMAX
                        CALL INTEGRAND_IN_M2OD_NSB(NTH,NP,NR,NSA,NSB,DFP,DFT,FFP,DCPP2,DCPT2,FCPP2,RSUM12)
                     END DO
                  END DO
               END DO
               DO NSB=1,NSBMAX
                  CALL p_theta_integration(RSUM12(NSB))
               END DO
               DO NSB=1,NSBMAX
                  RPCS2L_DEL(NR,NSB,NSA)=-RSUM12(NSB) &
                       *FACT*2.D0*PI*DELP(NS)*DELTH *1.D-6
               END DO
            END IF
            
         END DO ! NR
      END DO ! NSA

      CALL mtx_reset_communicator 

      END SUBROUTINE POWER_FROM_DIFFUSION_COEF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PHASE_SPACE_DERIVATIVE_F(NTH,NP,NR,NSA,SEND,DFP,DFT,FFP)

      IMPLICIT NONE
      integer,intent(in):: NTH,NP,NR,NSA
      integer:: NS
      double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND),intent(in)::SEND
      double precision,intent(out):: DFP, DFT, FFP
      double precision:: WPP, WPM, WPL

      NS=NS_NSA(NSA)

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
      DFP=    PG(NP,NS) &
           /DELP(NS)*(SEND(NTH,NP,NR,NSA)-SEND(NTH,NP-1,NR,NSA))
      IF(NTH.EQ.1) THEN

         DFT=1.D0/DELTH                             &
              *(                                     &
              ((1.D0-WPP)*SEND(NTH+1,NP  ,NR,NSA)   &
              +WPP *SEND(NTH+1,NP-1,NR,NSA)) &
              -                                    &
              ((1.D0-WPM)*SEND(NTH,NP  ,NR,NSA)     &
              +WPM *SEND(NTH,NP-1,NR,NSA)) &
              )
         
      ELSE IF(NTH.EQ.NTHMAX) THEN
         DFT=    1.D0/DELTH                         & 
              *(-                                    &
              ((1.D0-WPM)*SEND(NTH-1,NP  ,NR,NSA)   &
              +WPM *SEND(NTH-1,NP-1,NR,NSA)) &
              +                                     &
              ((1.D0-WPP)*SEND(NTH,NP  ,NR,NSA)     &
              +WPP *SEND(NTH,NP-1,NR,NSA)) &
              )
      ELSE
         DFT=    1.D0/(2.D0*DELTH)                  &
              *(                                     &
              ((1.D0-WPP)*SEND(NTH+1,NP  ,NR,NSA)   &
              +WPP *SEND(NTH+1,NP-1,NR,NSA)) &
              -                                    &
              ((1.D0-WPM)*SEND(NTH-1,NP  ,NR,NSA)   &
              +WPM *SEND(NTH-1,NP-1,NR,NSA)) &
              )
      ENDIF
      FFP=    PG(NP,NS)                           &
           *((1.D0-WPL)*SEND(NTH  ,NP  ,NR,NSA)  &
           +WPL *SEND(NTH  ,NP-1,NR,NSA))
      
      END SUBROUTINE PHASE_SPACE_DERIVATIVE_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INTEGRAND_IN_M2OD(NTH,NP,NR,NSA,DFP,DFT,FFP,DIFPP,DIFPT,FRICP,RSUM)

      IMPLICIT NONE
      integer,intent(in):: NTH,NP,NR,NSA 
      integer:: NS
      double precision,intent(in):: DFP,DFT,FFP      
      double precision,dimension(NTHMAX,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX),intent(in)::DIFPP, DIFPT, FRICP
      double precision,intent(inout):: RSUM
      double precision:: PV

      NS=NS_NSA(NSA)
      PV=SQRT(1.D0+THETA0(NS)*PG(NP,NS)**2)
      RSUM = RSUM+PG(NP,NS)**2*SINM(NTH)/PV   &
           *(DIFPP(NTH,NP,NR,NSA)*DFP           &
           +DIFPT(NTH,NP,NR,NSA)*DFT           &
           -FRICP(NTH,NP,NR,NSA)*FFP           &
           )

      END SUBROUTINE INTEGRAND_IN_M2OD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INTEGRAND_IN_M2OD_NSB(NTH,NP,NR,NSA,NSB,DFP,DFT,FFP,DIFPP,DIFPT,FRICP,RSUM)

      IMPLICIT NONE
      integer,intent(in):: NTH,NP,NR,NSA,NSB
      integer:: NS
      double precision,intent(in):: DFP,DFT,FFP      
      double precision,dimension(NTHMAX,NPSTART :NPENDWG,NRSTART:NRENDWM,NSBMAX,NSAMAX),intent(in)::DIFPP, DIFPT, FRICP
      double precision,dimension(NSBMAX),intent(inout):: RSUM
      double precision:: PV

      NS=NS_NSA(NSA)
      PV=SQRT(1.D0+THETA0(NS)*PG(NP,NS)**2)
      RSUM(NSB) = RSUM(NSB)+PG(NP,NS)**2*SINM(NTH)/PV   &
           *(DIFPP(NTH,NP,NR,NSB,NSA)*DFP           &
           +DIFPT(NTH,NP,NR,NSB,NSA)*DFT           &
           -FRICP(NTH,NP,NR,NSB,NSA)*FFP           &
           )

      END SUBROUTINE INTEGRAND_IN_M2OD_NSB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE POWER_FROM_SOURCE_TERM

      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      INTEGER:: NTH, NP, NR, NSA, NS
      double precision:: FACT, RSUM11B, RSUM11F, RSUM11S, RSUM11L, RSUM11S_CX, PV

      CALL mtx_set_communicator(comm_np)
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            RSUM11B=0.D0
            RSUM11F=0.D0
            RSUM11S=0.D0
            RSUM11L=0.D0
            RSUM11S_CX=0.D0

            IF(MODELR.eq.1)THEN
               DO NP=NPSTART,NPEND
                  PV=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
                  DO NTH=1,NTHMAX
                     RSUM11B = RSUM11B + PM(NP,NS)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NS)*SPPB(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR) 
                     RSUM11F = RSUM11F + PM(NP,NS)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NS)*SPPF(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR) 
                     RSUM11S = RSUM11S + PM(NP,NS)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NS)*SPPL(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR) 
                     IF(MODEL_DELTA_F(NS).eq.0)THEN
                        RSUM11L = RSUM11L + PM(NP,NS)**2*SINM(NTH) &
                             *(PV-1.D0)/THETA0(NS)* PPL(NTH,NP,NR,NSA)&
                             *FNSP(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)!*RFSADG(NR) PPL includes RLAMDA 
                     ELSEIF(MODEL_DELTA_F(NS).eq.1)THEN
                        RSUM11L = RSUM11L + PM(NP,NS)**2*SINM(NTH) &
                             *(PV-1.D0)/THETA0(NS)* PPL(NTH,NP,NR,NSA)&
                             *FNSP_DEL(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)!*RFSADG(NR) PPL includes RLAMDA 
                     END IF
                     RSUM11S_CX = RSUM11S_CX + PM(NP,NS)**2*SINM(NTH) &
                          *(PV-1.D0)/THETA0(NS)*SPPL_CX(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR) 
                  END DO
               END DO
            ELSE ! MODELR=0
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     RSUM11B = RSUM11B + PM(NP,NS)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NS)**2*SPPB(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
                     RSUM11F = RSUM11F + PM(NP,NS)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NS)**2*SPPF(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
                     RSUM11S = RSUM11S + PM(NP,NS)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NS)**2*SPPL(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
                     IF(MODEL_DELTA_F(NS).eq.0)THEN
                        RSUM11L = RSUM11L + PM(NP,NS)**2*SINM(NTH) &
                             *0.5D0*PM(NP,NS)**2*PPL(NTH,NP,NR,NSA) &
                             *FNSP(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)!*RFSADG(NR)
                     ELSEIF(MODEL_DELTA_F(NS).eq.1)THEN
                        RSUM11L = RSUM11L + PM(NP,NS)**2*SINM(NTH) &
                             *0.5D0*PM(NP,NS)**2*PPL(NTH,NP,NR,NSA) &
                             *FNSP_DEL(NTH,NP,NR,NSA)!*RLAMDA(NTH,NR)!*RFSADG(NR)
                     END IF
                     RSUM11S_CX = RSUM11S_CX + PM(NP,NS)**2*SINM(NTH) &
                          *0.5D0*PM(NP,NS)**2*SPPL_CX(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)!*RFSADG(NR)
                  END DO
               END DO
            END IF
            CALL p_theta_integration(RSUM11B) ! beam power
            CALL p_theta_integration(RSUM11F) ! fusion power
            CALL p_theta_integration(RSUM11S) ! tloss power
            CALL p_theta_integration(RSUM11L) ! PPL power including TLOSS and thermalization
            CALL p_theta_integration(RSUM11S_CX) ! CX loss power

            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)*RFSADG(NR)
            RSPBL(NR,NSA)= RSUM11B*FACT*2.D0*PI*DELP(NS)*DELTH*1.D-6
            RSPFL(NR,NSA)= RSUM11F*FACT*2.D0*PI*DELP(NS)*DELTH*1.D-6
            RSPSL(NR,NSA)= RSUM11S*FACT*2.D0*PI*DELP(NS)*DELTH*1.D-6
            RSPLL(NR,NSA)= RSUM11L*FACT*2.D0*PI*DELP(NS)*DELTH*1.D-6
            RSPSL_CX(NR,NSA)= RSUM11S_CX*FACT*2.D0*PI*DELP(NS)*DELTH*1.D-6

         END DO
      END DO

      CALL mtx_reset_communicator 

      END SUBROUTINE POWER_FROM_SOURCE_TERM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE POWER_FROM_RADIAL_TRANSPORT
      
      USE fpmpi
      USE libmpi
      IMPLICIT NONE
      INTEGER:: NTH, NP, NR, NSA, NS
      real(8):: WRL, WRH, DINT_DFDT_R1, DINT_DFDT_R2
      real(8):: DFDR_R1, DFDR_R2, DFDT_R1, DFDT_R2, DINT_DR, RSUM_DR,RGAMA,F_R1,F_R2, RSUMN_DR
      real(8):: SRHOR1, SRHOR2, RSUM_DRS, RSUMN_DRS
      real(8)::DRRM, DRRP, RL

      IF(MODELD.ne.0)THEN
         CALL mtx_set_communicator(comm_np)
         DO NR=NRSTART,NREND
            DO NSA=NSASTART,NSAEND
               NS=NS_NSA(NSA)

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
               DO NP=NPSTART,NPEND
                  RGAMA=SQRT(1.D0+THETA0(NS)*PM(NP,NS)**2)
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

                     RSUMN_DR=RSUMN_DR +             DINT_DR*VOLP(NTH,NP,NS)
                     IF(MODELR.eq.1)THEN
                        RSUM_DR=RSUM_DR+(RGAMA-1.D0)/THETA0(NS)*DINT_DR*VOLP(NTH,NP,NS)
                     ELSEIF(MODELR.eq.0)THEN
                        RSUM_DR=RSUM_DR + 0.5D0*PG(NP,NS)**2*DINT_DR*VOLP(NTH,NP,NS)
                     END IF
!
                     IF(NR.eq.NRMAX)THEN
                        RSUMN_DRS=RSUMN_DRS +      DFDT_R2*VOLP(NTH,NP,NS)/DELR *RFSADG(NR)/RL
                        IF(MODELR.eq.1)THEN
                           RSUM_DRS=RSUM_DRS+(RGAMA-1.D0)/THETA0(NS)*DFDT_R2*VOLP(NTH,NP,NS)/DELR *RFSADG(NR)
                        ELSE
                           RSUM_DRS=RSUM_DRS+0.5D0*PG(NP,NS)**2*DFDT_R2*VOLP(NTH,NP,NS)/DELR *RFSADG(NR)
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
! --------- end of radial transport
            ENDDO ! NSA
         ENDDO ! NR
         CALL mtx_reset_communicator 
      END IF ! MODELD

      END SUBROUTINE POWER_FROM_RADIAL_TRANSPORT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================================================


      end MODULE fpsave




