!
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

      USE eg_read
      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NS
      double precision:: FL

      
      IF(NT_init.eq.0)THEN
         WRITE(33,'(A,4I5,I8)') "# Number of mesh, theta, p, rho, s, total ", &
              NTHMAX, NPMAX, NRMAX, NSAMAX, NTHMAX*NPMAX*NRMAX*NSAMAX
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            WRITE(33,'(A,I4)') "# species", NSA
            WRITE(33,'(A,3E14.6)') "# norm cont v, p, n= ", VTFP0(NSA), PTFP0(NSA), RNFP0(NSA)*1.D20
            WRITE(33,'(A,I5)') "# p grid", NPMAX
            WRITE(33,'(200E14.6)') (PM(NP,NS), NP=1, NPMAX)
         END DO
         WRITE(33,'(A,I5)') "# theta grid", NTHMAX
         WRITE(33,'(200E14.6)') (THM(NTH), NTH=1, NTHMAX)
         WRITE(33,'(A,I5)') "# rho grid", NRMAX
         WRITE(33,'(200E14.6)') (RM(NR), NR=1, NRMAX)
         
         WRITE(33,'(A)') "# NSA, NR, NP, NTHETA, del_f"
      END IF
      
      IF(NT_init.ne.0)THEN
         WRITE(33,'(A,E14.6)') "# TIME= ", TIMEFP
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1, NRMAX
               DO NP=1, NPMAX
                  FL=FPMXWL_EXP(PM(NP,NS),NR,NS)
                  DO NTH=1, NTHMAX
                     WRITE(33,'(4I5,3E16.8)') NSA, NR, NP, NTH, &
                          FNS(NTH,NP,NR,NSA)-FL
                  END DO
               END DO
            END DO
         END DO
      END IF

      END SUBROUTINE OUT_TXT_FNS_DEL
!------------------------------------------      
      SUBROUTINE OUT_TXT_F1
      
      IMPLICIT NONE
      INTEGER:: NP, NTH, NS_F1
      double precision:: pitch_angle_av

      NS_F1=NS_NSA(NSA_F1) 
      DO NP=1,NPMAX
!pitch angle average
         pitch_angle_av = 0.D0
         DO NTH=1,NTHMAX
            pitch_angle_av = pitch_angle_av + FNS(NTH,NP,NR_F1,NSA_F1)*SINM(NTH)*0.5D0*RLAMDAG(NTH,NR_F1)*RFSADG(NR_F1)
         END DO
         pitch_angle_av = pitch_angle_av/NTHMAX
         WRITE(9,'(1PE12.4,I6,1P30E17.8e3)') PTG(NTG1)*1000, NP, PM(NP,NS_F1), &
              PM(NP,NS_F1)*PTFP0(NSA_F1)/AMFP(NSA_F1)/VC/SQRT(1.D0+PM(NP,NS_F1)**2*THETA0(NS_F1)), &
              PM(NP,NS_F1)**2, &
              0.5D0*PTFP0(NSA_F1)**2*PM(NP,NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3), &
              FNS(1,NP,NR_F1,NSA_F1), FNS(NTH_F1,NP,NR_F1,NSA_F1), &
              FNS(NTHMAX/2,NP,NR_F1,NSA_F1), pitch_angle_av
      END DO
      WRITE(9,*) " "
      WRITE(9,*) " "

      END SUBROUTINE OUT_TXT_F1
!------------------------------------------      
      SUBROUTINE OUT_TXT_BEAM_WIDTH

      IMPLICIT NONE
      INTEGER:: NP
      double precision:: beam_peak_value
      integer:: NP_2e_h, NP_2e_l, NP_1e_h, NP_1e_l, NP_half_l, NP_half_h, NP_B_peak, NS_F1

      NS_F1=NS_NSA(NSA_F1)
      NP_B_peak=0
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
         IF(FNS(NTH_F1,NP,NR_F1,NSA_F1).le.5.D-6.and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).ge.5.D-6)THEN
            NP_1e_l=NP
         ELSEIF(FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.5.D-6.and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).le.5.D-6)THEN
            NP_1e_h=NP
         END IF
      END DO

      DO NP=NP_BULK(NR_F1,NSA_F1), NPMAX-1 ! half value width
         IF(FNS(NTH_F1,NP-1,NR_F1,NSA_F1).le.FNS(NTH_F1,NP,NR_F1,NSA_F1).and.FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.FNS(NTH_F1,NP+1,NR_F1,NSA_F1))THEN
            beam_peak_value=FNS(NTH_F1,NP,NR_F1,NSA_F1)*0.5D0
            NP_B_peak=NP
         END IF
      END DO
      DO NP=1, NPMAX-1
         IF(FNS(NTH_F1,NP,NR_F1,NSA_F1).le.beam_peak_value.and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).ge.beam_peak_value)THEN
            NP_half_l=NP
         ELSEIF(FNS(NTH_F1,NP,NR_F1,NSA_F1).ge.beam_peak_value.and.FNS(NTH_F1,NP+1,NR_F1,NSA_F1).le.beam_peak_value)THEN
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
           PTFP0(NSA_F1)**2*PM(NP_half_l,NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3)*0.5D0, &
           PTFP0(NSA_F1)**2*PM(NP_half_h,NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3)*0.5D0, &
           PTFP0(NSA_F1)**2*PM(NP_B_peak,NS_F1)**2/(AEE*AMFP(NSA_F1)*1.D3)*0.5D0, &
           FNS(NTH_F1,NP_half_l,NR_F1,NSA_F1), FNS(NTH_F1,NP_half_h,NR_F1,NSA_F1), beam_peak_value*2.D0
      WRITE(*,'(A,E14.6,6I5,E14.6,I5)') "BEAM_HALF ",TIMEFP, NP_2e_l, NP_2e_h, NP_1e_l, NP_1e_h, NP_half_l, NP_half_h, beam_peak_value, NP_BULK(NR_F1,NSA_F1)

      
      END SUBROUTINE OUT_TXT_BEAM_WIDTH
!------------------------------------------      
      SUBROUTINE NUMBER_OF_NONTHERMAL_IONS

      IMPLICIT NONE
      integer:: NP, NS_F1
      double precision:: RSUM_BEAM

      NS_F1=NS_NSA(NSA_F1)
      CALL mtx_set_communicator(comm_np)
      RSUM_BEAM=0.D0

      IF(NSASTART.eq.NSA_F1.and.NRSTART.eq.NR_F1)THEN
         DO NP=NPSTART, NPEND
            RSUM_BEAM = RSUM_BEAM + FNSP_DEL(NTH_F1,NP,NR_F1,NSA_F1)*VOLP(NTH_F1,NP,NS_F1)
         END DO
      END IF
      CALL p_theta_integration(RSUM_BEAM)
      IF(NSASTART.eq.NSA_F1.and.NRSTART.eq.NR_F1.and.OUTPUT_TXT_BEAM_WIDTH.eq.1)THEN
         WRITE(31,'(2E14.6)') TIMEFP, RSUM_BEAM
      END IF

      CALL mtx_reset_communicator

      END SUBROUTINE NUMBER_OF_NONTHERMAL_IONS
!------------------------------------------      

      END MODULE fpoutdata
