!     ***********************************************************
!
!           TO OUTPUT SEVERAL DATA TO .dat FILE
!
!     ***********************************************************
!
      MODULE fpoutdata

      USE fpcomm
      use libmpi
      use fpmpi

      contains
!----------------------------------------------------------------
      SUBROUTINE OUT_TXT_FNS_DEL

      USE fpfunc
      USE eg_read
      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NS
      double precision:: FL

      
      IF(NT_init.eq.0)THEN
         WRITE(35,'(A)') "# Number of mesh, theta, p, rho, s, total= "
         WRITE(35,'(4I5,I8)') NTHMAX, NPMAX, NRMAX, NSAMAX, NTHMAX*NPMAX*NRMAX*NSAMAX
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            WRITE(35,'(A)') "# species="
            WRITE(35,'(I4)') NSA
            WRITE(35,'(A)') "# norm cont v, p, n= "
            WRITE(35,'(3E14.6)') VTFP0(NSA), PTFP0(NSA), RNFP0(NSA)*1.D20
            WRITE(35,'(A)') "# normalized p values"
            WRITE(35,'(200E14.6)') (PM(NP,NS), NP=1, NPMAX)
         END DO
         WRITE(35,'(A)') "# theta values"
         WRITE(35,'(200E14.6)') (THM(NTH), NTH=1, NTHMAX)
         WRITE(35,'(A)') "# rho values"
         WRITE(35,'(200E14.6)') (RM(NR), NR=1, NRMAX)
      END IF
      
      IF(NT_init.ne.0)THEN
         WRITE(33,'(A)') " "
         WRITE(33,'(A)') " "
         WRITE(33,'(A,E14.6)') "# TIME= ", TIMEFP+time_exp_offset
         WRITE(33,'(A)') "# NSA, NR, NP, NTHETA, f, f_maxwell"
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1, NRMAX
               DO NP=1, NPMAX
                  FL=FPMXWL_EXP(PM(NP,NS),NR,NS)
                  DO NTH=1, NTHMAX
                     WRITE(33,'(4I5,3E16.8)') NSA, NR, NP, NTH, &
                          FNS(NTH,NP,NR,NSA), FL
                  END DO
               END DO
            END DO
         END DO
      END IF

      END SUBROUTINE OUT_TXT_FNS_DEL
!------------------------------------------      
      SUBROUTINE OUT_TXT_F1

      USE fpfunc
      USE EG_READ
      IMPLICIT NONE
      INTEGER:: NP, NTH, NS_F1, DNTH
      double precision:: pitch_angle_av, FL, PML, av, sum, av2

      NS_F1=NS_NSA(NSA_F1) 
      WRITE(9,'(A,E14.6,I4,E14.6)') "# TIME, NP_BULK, E_bulk ", TIMEFP+time_exp_offset, NP_BULK(NR_F1,NSA_F1), &
           0.5D0*PTFP0(NSA_F1)**2*PM(NP_BULK(NR_F1,NSA_F1),NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3)
      DO NP=1,NPMAX
!pitch angle average
         pitch_angle_av = 0.D0
         av2=0.D0
         DO NTH=1,NTHMAX
            pitch_angle_av = pitch_angle_av + FNS(NTH,NP,NR_F1,NSA_F1)*SINM(NTH)*RLAMDAG(NTH,NR_F1)*RFSADG(NR_F1)
            av2 = av2 + SINM(NTH)**RLAMDAG(NTH,NR_F1)*RFSADG(NR_F1)
         END DO
         pitch_angle_av = pitch_angle_av/av2
         PML=PM(NP,NS_F1)
         FL=FPMXWL_EXP(PML,NR_F1,NS_F1)
!         sum=0.D0
!         av=0.D0
!         DO NTH=1, NTHMAX
!            DNTH=ABS(NTH_F1-NTH)
!            IF(DNTH.lt.NTHMAX/2)THEN
!               sum = sum + FNS(NTH,NP,NR_F1,NSA_F1)*COSG(DNTH)
!               av = av + COSG(DNTH)
!            END IF
!         END DO
         WRITE(9,'(I6,1P30E17.8e3)') NP, PM(NP,NS_F1), &
              PM(NP,NS_F1)*PTFP0(NSA_F1)/AMFP(NSA_F1)/VC/SQRT(1.D0+PM(NP,NS_F1)**2*THETA0(NS_F1)), &
              0.5D0*PTFP0(NSA_F1)**2*PM(NP,NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3), &
              FNS(NTHMAX+1-NTH_F1,NP,NR_F1,NSA_F1), FNS(NTH_F1,NP,NR_F1,NSA_F1), &
              FNS(NTHMAX/2,NP,NR_F1,NSA_F1), pitch_angle_av, FL!, sum/av
      END DO
      WRITE(9,*) " "
      WRITE(9,*) " "

      END SUBROUTINE OUT_TXT_F1
!------------------------------------------      
      SUBROUTINE OUT_TXT_NF_RADIAL_PROF(ID)

      USE fpnfrr
      IMPLICIT NONE
      integer,intent(in):: ID
      integer:: NR, ndata, NS
      double precision,dimension(NRSTART:NREND):: send1, send2, send3
      double precision,dimension(NRMAX):: recv1, recv2, recv3
      double precision:: tot1, tot2, tot3

      CALL mtx_set_communicator(comm_nr)
      ndata=NREND-NRSTART+1
      DO NR=NRSTART, NREND
         send1(NR)=RATE_NF(NR,ID)
         send2(NR)=RATE_NF_BT(NR,ID)
         send3(NR)=RATE_NF_TT(NR,ID)
      END DO
      CALL mtx_allgather_real8(send1,ndata,recv1)
      CALL mtx_allgather_real8(send2,ndata,recv2)
      CALL mtx_allgather_real8(send3,ndata,recv3)

      CALL mtx_reset_communicator

      IF(NRANK.eq.0)THEN
         WRITE(27,*) " "
         WRITE(27,*) " "
         WRITE(27,'(A,2E14.6)') "# NF_PROF ", time_exp_offset, TIMEFP
         DO NR=1, NRMAX
            WRITE(27,'(99E14.5e3)') RM(NR), recv1(NR)*VOLR(NR), recv2(NR)*VOLR(NR), recv3(NR)*VOLR(NR), &
                 recv1(NR), recv2(NR), recv3(NR)
         END DO
      END IF

      END SUBROUTINE OUT_TXT_NF_RADIAL_PROF
!------------------------------------------      
      SUBROUTINE OUT_TXT_F1_PITCH

      USE fpfunc
      IMPLICIT NONE
      INTEGER:: NTH, NP, NS_F1
      double precision:: RSUM_deno, RSUM_nume, PML, F_MX, RSUM_deno_high, RSUM_nume_high, Energy

      NS_F1=NS_NSA(NSA_F1) 
      WRITE(36,'(A,E14.6,I4,E14.6)') "# TIME, NP_BULK, E_bulk ", TIMEFP+time_exp_offset, NP_BULK(NR_F1,NSA_F1), &
           0.5D0*PTFP0(NSA_F1)**2*PM(NP_BULK(NR_F1,NSA_F1),NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3)
      DO NTH=1, NTHMAX
         RSUM_deno=0.D0
         RSUM_nume=0.D0
         RSUM_deno_high=0.D0
         RSUM_nume_high=0.D0
         DO NP=1, NPMAX
            PML=PM(NP,NS_F1)
            F_MX=FPMXWL_EXP(PML,NR_F1,NS_F1)
            RSUM_deno = RSUM_deno + PM(NP,NS_F1)**2*DELP(NS_F1)
            RSUM_nume = RSUM_nume + PM(NP,NS_F1)**2*DELP(NS_F1)*(FNS(NTH,NP,NR_F1,NSA_F1)-F_MX)
            Energy=0.5D0*PTFP0(NSA_F1)**2*PM(NP,NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3)
            IF(Energy.ge.80.D0)THEN 
               RSUM_deno_high = RSUM_deno_high + PM(NP,NS_F1)**2*DELP(NS_F1)
               RSUM_nume_high = RSUM_nume_high + PM(NP,NS_F1)**2*DELP(NS_F1)*(FNS(NTH,NP,NR_F1,NSA_F1)-F_MX)
            END IF
         END DO
         WRITE(36,'(3E14.5e3)') THM(NTH), RSUM_nume/RSUM_deno, RSUM_nume_high/RSUM_deno_high
      END DO
      WRITE(36,*) " "
      WRITE(36,*) " "

      END SUBROUTINE OUT_TXT_F1_PITCH
!------------------------------------------      
      SUBROUTINE OUT_TXT_BEAM_WIDTH

      IMPLICIT NONE
      INTEGER:: NP
      double precision:: beam_peak_value
      integer:: NP_2e_h, NP_2e_l, NP_1e_h, NP_1e_l, NP_half_l, NP_half_h, NP_B_peak, NS_F1, i, j
      double precision:: log10_neu0, log10_neus, alpha, beta, k_energy, log_energy
      double precision:: N_NEUT, sigma_cx

      NS_F1=NS_NSA(NSA_F1)
      NP_B_peak=1
      NP_2e_l=1
      NP_2e_h=1
      NP_1e_l=1
      NP_1e_h=1
      NP_half_l=1
      NP_half_h=1
      beam_peak_value=0.D0
      DO NP=1, NPMAX-1
         IF(FNS(NTH_F1,NP,NR_F1,NSA_F1).le.2.D-5.and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).ge.2.D-5)THEN
            NP_2e_l=NP
         ELSEIF(FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.2.D-5.and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).le.2.D-5)THEN
            NP_2e_h=NP
         END IF
         IF(FNS(NTH_F1,NP,NR_F1,NSA_F1).le.5.D-6 &
              .and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).ge.5.D-6)THEN
            NP_1e_l=NP
         ELSEIF(FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.5.D-6 &
              .and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).le.5.D-6)THEN
            NP_1e_h=NP
         END IF
      END DO

      DO NP=NP_BULK(NR_F1,NSA_F1), NPMAX-1 ! half value width
         IF(FNS(NTH_F1,NP-1,NR_F1,NSA_F1).le.FNS(NTH_F1,NP,NR_F1,NSA_F1) &
         .and.FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.FNS(NTH_F1,NP+1,NR_F1,NSA_F1))THEN
            beam_peak_value=FNS(NTH_F1,NP,NR_F1,NSA_F1)*0.5D0
            NP_B_peak=NP
         END IF
      END DO
      DO NP=1, NPMAX-1
         IF(FNS(NTH_F1,NP,NR_F1,NSA_F1).le.beam_peak_value &
              .and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).ge.beam_peak_value)THEN
            NP_half_l=NP
         ELSEIF(FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.beam_peak_value &
              .and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).le.beam_peak_value)THEN
            NP_half_h=NP
         END IF
      END DO

      WRITE(24,'(E14.6,3I5,20E14.6)') TIMEFP, &
!              NP_2e_l, NP_2e_h, &
           NP_half_l, NP_half_h, &
           NP_B_peak, &
!              PTFP0(2)**2*PM(NP_2e_l,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              PTFP0(2)**2*PM(NP_2e_h,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              FNS(NTH_F1,NP_2e_l,NR_F1,2), FNS(NTH_F1,NP_2e_h,NR_F1,2), &
!              PTFP0(2)**2*PM(NP_1e_l,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              PTFP0(2)**2*PM(NP_1e_h,2)**2/(AEE*AMFP(2)*1.D3)*0.5D0, &
!              FNS(NTH_F1,NP_1e_l,NR_F1,2), FNS(NTH_F1,NP_1e_h,NR_F1,2), &
           PTFP0(NSA_F1)**2*PM(NP_half_l,NS_F1)**2 &
           /(AEE*AMFP(NSA_F1)*1.D3)*0.5D0, &
           PTFP0(NSA_F1)**2*PM(NP_half_h,NS_F1)**2 &
           /(AEE*AMFP(NSA_F1)*1.D3)*0.5D0, &
           PTFP0(NSA_F1)**2*PM(NP_B_peak,NS_F1)**2 &
           /(AEE*AMFP(NSA_F1)*1.D3)*0.5D0, &
           FNS(NTH_F1,NP_half_l,NR_F1,NSA_F1), &
           FNS(NTH_F1,NP_half_h,NR_F1,NSA_F1), beam_peak_value*2.D0
      WRITE(*,'(A,E14.6,6I5,E14.6,I5)') "BEAM_HALF ",TIMEFP, &
           NP_2e_l, NP_2e_h, NP_1e_l, NP_1e_h, NP_half_l, NP_half_h, &
           beam_peak_value, NP_BULK(NR_F1,NSA_F1)


!      log10_neu0=LOG10(RN_NEU0)
!      log10_neus=LOG10(RN_NEUS)
!!     neutral gas profile
!      alpha=2.5D0
!      beta=0.8D0
!      N_NEUT=10**( (log10_neu0-log10_neus)*(1.D0-RM(NR_F1)**alpha)**beta+log10_neus )
!
!      k_energy = 0.5D0*(PTFP0(NSA_F1)*PM(NP_B_peak,NS_F1))**2/(AMFP(NSA_F1)*AEE)
!      log_energy = dlog10(k_energy)
!      
!      sigma_cx = 0.6937D-14*(1.D0-0.155D0*log_energy)**2/ &
!           (1.D0+0.1112D-14*k_energy**3.3D0)*1.D-4
!      WRITE(34,'(1P5E14.6)') TIMEFP, k_energy*1.D-3, sigma_cx, &
!           1.D0/(sigma_cx*VTFP0(NSA_F1)*PM(NP_B_peak,NS_F1)*N_NEUT*1.D20), &
!           1.D0/(sigma_cx*VTFP0(NSA_F1)*PM(NP_B_peak,NS_F1)*N_NEUT*1.D20*EXP(1.D0))

      END SUBROUTINE OUT_TXT_BEAM_WIDTH
!------------------------------------------      
      SUBROUTINE NUMBER_OF_NONTHERMAL_IONS ! beam direction

      IMPLICIT NONE
      integer:: NP, NS_F1, NR, NS
      double precision:: RSUM1, RSUM2, RSUM3, RSUM4, RSUM5, RSUM_ZEFF, RSUM_deno
!      double precision,dimension(NRSTART:NREND):: RNSL_BEAM
!      double precision,dimension(NRMAX):: RNS_BEAM

      RSUM1=0.D0
      RSUM2=0.D0
      RSUM3=0.D0
      RSUM4=0.D0
      RSUM5=0.D0
      RSUM_ZEFF=0.D0
      RSUM_deno=0.D0
      NS_F1=NS_NSA(NSA_F1)
      DO NR=1, NRMAX
         RSUM1=RSUM1+VOLR(NR)*RNS_DELF(NR,NS_F1)
         RSUM2=RSUM2+VOLR(NR)*RWS_PARA(NR,NS_F1)
         RSUM3=RSUM3+VOLR(NR)*RWS_PERP(NR,NS_F1)
         RSUM4=RSUM4+VOLR(NR)*RNS(NR,NSA_F1)
         RSUM5=RSUM5+VOLR(NR)*RN_TEMP(NR,1)
      END DO
      IF(MODEL_DELTA_F_CN.eq.0)THEN
         DO NS=2, NSBMAX
            RSUM_ZEFF = RSUM_ZEFF + (RN_BULK(1,NS))*PZ(NS)**2
            RSUM_deno = RSUM_deno + (RN_BULK(1,NS))*PZ(NS)
         END DO
      ELSEIF(MODEL_DELTA_F_CN.eq.1)THEN
         DO NS=2, NSBMAX
            RSUM_ZEFF = RSUM_ZEFF + (RN_BULK(1,NS)+RNS_DELF(1,NS))*PZ(NS)**2
            RSUM_deno = RSUM_deno + (RN_BULK(1,NS)+RNS_DELF(1,NS))*PZ(NS)
         END DO
      END IF
!      RSUM_ZEFF = RSUM_ZEFF/((RN_TEMP(1,1))*AEFD(1)**2)
      RSUM_ZEFF = RSUM_ZEFF/RSUM_deno

      WRITE(31,'(210E14.6)') TIMEFP, (RNS_DELF(NR,NS_F1),NR=1,NRMAX), RSUM1,&
           (RWS_PARA(NR,NS_F1),NR=1,NRMAX), RSUM2, &
           (RWS_PERP(NR,NS_F1),NR=1,NRMAX), RSUM3, &
           (RNS(NR,NSA_F1)/RN_TEMP(NR,1),NR=1,NRMAX), RSUM4/RSUM5, RSUM4, RSUM_ZEFF
!      WRITE(31,'(52E14.6,16I5)') TIMEFP, (RNS_DELF(NR,NS_F1),NR=1,NRMAX), RSUM1,&
!           (RWS_PARA(NR,NS_F1),NR=1,NRMAX), RSUM2,&
!           (RWS_PERP(NR,NS_F1),NR=1,NRMAX), RSUM3, (NP_BULK(NR,1),NR=1,NRMAX)

      END SUBROUTINE NUMBER_OF_NONTHERMAL_IONS
!------------------------------------------      
      SUBROUTINE PROF_OF_NF_REACTION_RATE(ID)

      USE fpnfrr
      IMPLICIT NONE
      integer,intent(in):: ID
      integer:: NR, ndata, NS
      double precision,dimension(NRSTART:NREND):: send1, send2, send3
      double precision,dimension(NRMAX):: recv1, recv2, recv3
      double precision:: tot1, tot2, tot3

      CALL mtx_set_communicator(comm_nr)
      ndata=NREND-NRSTART+1
      DO NR=NRSTART, NREND
         send1(NR)=RATE_NF(NR,ID)
         send2(NR)=RATE_NF_BT(NR,ID)
         send3(NR)=RATE_NF_TT(NR,ID)
      END DO
      CALL mtx_allgather_real8(send1,ndata,recv1)
      CALL mtx_allgather_real8(send2,ndata,recv2)
      CALL mtx_allgather_real8(send3,ndata,recv3)

      CALL mtx_reset_communicator

      tot1=0.D0
      tot2=0.D0
      tot3=0.D0
      NS=NS_NSA(NSA_F1)
      DO NR=1, NRMAX
         tot1 = tot1 + recv1(NR)*VOLR(NR)
         tot2 = tot2 + recv2(NR)*VOLR(NR)
         tot3 = tot3 + recv3(NR)*VOLR(NR)
      END DO

      IF(NRANK.eq.0)THEN
         WRITE(25,'(99E14.6)') TIMEFP, (recv1(NR), NR=1, NRMAX), (recv2(NR), NR=1, NRMAX), (recv3(NR), NR=1, NRMAX)
         WRITE(26,'(99E14.5e3)') TIMEFP, tot1, tot2, tot3, (recv1(NR)*VOLR(NR), NR=1, NRMAX)
      END IF

      IF(NRANK.eq.0.and.TIMEFP+time_exp_offset.eq.4.8D0)THEN
         DO NR=1, NRMAX
            WRITE(27,'(99E14.6)') RM(NR), recv1(NR)*VOLR(NR), recv2(NR)*VOLR(NR), recv3(NR)*VOLR(NR), &
                 recv1(NR), recv2(NR), recv3(NR)
         END DO
      END IF

      END SUBROUTINE PROF_OF_NF_REACTION_RATE
!------------------------------------------      
      SUBROUTINE OUT_TXT_HEAT_PROF

      IMPLICIT NONE
      INTEGER:: NSA, NSB, NR
      double precision,dimension(NSBMAX):: RSUM_PC, RSUM_PC_DEL
      double precision:: RSUM_RSPL, RSUM_RSPB

!      IF(NRANK.eq.0.and.MODEL_NBI.eq.2)THEN
      NSA=NSA_F1
      RSUM_PC(:)=0.D0
      RSUM_PC_DEL(:)=0.D0
      RSUM_RSPL=0.D0
      RSUM_RSPB=0.D0
      DO NSB=1, NSBMAX
         DO NR=1, NRMAX
            RSUM_PC(NSB)=RSUM_PC(NSB)+RPCS2(NR,NSB,NSA)*VOLR(NR)
            RSUM_PC_DEL(NSB)=RSUM_PC_DEL(NSB)+RPCS2_DEL(NR,NSB,NSA)*VOLR(NR)
         END DO
      END DO
      DO NR=1, NRMAX
         RSUM_RSPL=RSUM_RSPL+RSPL(NR,NSA)*VOLR(NR)
         RSUM_RSPB=RSUM_RSPB+RSPB(NR,NSA)*VOLR(NR)
      END DO

      WRITE(29,'(99E14.6)') TIMEFP, ((RPCS2(NR,NSB,NSA),NR=1,NRMAX),NSB=1,NSBMAX), (RSUM_PC(NSB),NSB=1,NSBMAX)
      WRITE(30,'(99E14.6)') TIMEFP, ((RPCS2_DEL(NR,NSB,NSA),NR=1,NRMAX),NSB=1,NSBMAX), (RSUM_PC_DEL(NSB),NSB=1,NSBMAX)
      WRITE(32,'(99E14.6)') TIMEFP, (RSPL(NR,NSA),NR=1,NRMAX), RSUM_RSPL
      WRITE(34,'(99E14.6)') TIMEFP, (RSPB(NR,NSA),NR=1,NRMAX), RSUM_RSPB
!      END IF

      END SUBROUTINE OUT_TXT_HEAT_PROF
!------------------------------------------      

      END MODULE fpoutdata
