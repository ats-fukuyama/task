!
! ************************************************************
!
!            CALCULATION OF DW
!
! ************************************************************

MODULE fpcalw

  USE fpcomm
  USE fpcalwm
  USE fpcalwr

contains

!------------------------------------------

  SUBROUTINE FP_CALW

    IMPLICIT NONE
    integer:: NSA, NR, NTH, NP, NS
    integer:: NCONST_RF
    real(kind8):: DWTTEC, DWTTIC, DWTPEC, DWTPIC

!     ----- Initialize ------------------------------------- 

    DO NSA=NSASTART,NSAEND
       DO NR=NRSTART,NREND
          DO NP=NPSTART,NPENDWG
             DO NTH=1,NTHMAX
                DWECPP(NTH,NP,NR,NSA)=0.D0
                DWECPT(NTH,NP,NR,NSA)=0.D0
             END DO
          END DO
          DO NP=NPSTARTW,NPENDWM
             DO NTH=1,NTHMAX+1
                DWECTP(NTH,NP,NR,NSA)=0.D0
                DWECTT(NTH,NP,NR,NSA)=0.D0
             END DO
          END DO
       END DO
    END DO

!     ECRF 
      IF(DEC.ne.0.and.NSA.eq.1) THEN
         CALL FP_CALW0
         DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  DWECPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA)
                  DWECPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA)
                  DWPP(NTH,NP,NR,NSA)=0.D0
                  DWPT(NTH,NP,NR,NSA)=0.D0
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  DWECTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA)
                  DWECTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA)
                  DWTP(NTH,NP,NR,NSA)=0.D0
                  DWTT(NTH,NP,NR,NSA)=0.D0
               END DO
            END DO
         END DO
         END DO
      END IF

!     ICRF
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         IF(MODELW(NS).EQ.1) THEN
            CALL FP_CALWR
         ELSEIF(MODELW(NS).EQ.2) THEN
            CALL FP_CALWR
         ELSEIF(MODELW(NS).EQ.3) THEN
            CALL FP_CALWM
         ELSEIF(MODELW(NS).EQ.4) THEN
            CALL FP_CALWM
         ELSEIF(MODELW(NS).ne.0) THEN
            IF(nrank.eq.0) WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW(NS)
         ENDIF
      END DO

      DO NSA=NSASTART,NSAEND
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DWICPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA)
               DWICPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA)
               DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA) + DWECPP(NTH,NP,NR,NSA)
               DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA) + DWECPT(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) + DWECTP(NTH,NP,NR,NSA)
               DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) + DWECTT(NTH,NP,NR,NSA)
            END DO
         END DO
      END DO
      END DO

!     POOL coef DW in order to reduce the number of call fp_calwm
      IF(N_IMPL.EQ.0)THEN ! N_IMPL=0
      DO NSA=NSASTART, NSAEND
      DO NR=NRSTART, NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DWPP_P(NTH,NP,NR,NSA) = DWPP(NTH,NP,NR,NSA)
               DWPT_P(NTH,NP,NR,NSA) = DWPT(NTH,NP,NR,NSA)
               DWICPP_P(NTH,NP,NR,NSA) = DWICPP(NTH,NP,NR,NSA)
               DWICPT_P(NTH,NP,NR,NSA) = DWICPT(NTH,NP,NR,NSA)
               DWECPP_P(NTH,NP,NR,NSA) = DWECPP(NTH,NP,NR,NSA)
               DWECPT_P(NTH,NP,NR,NSA) = DWECPT(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWTP_P(NTH,NP,NR,NSA) = DWTP(NTH,NP,NR,NSA)
               DWTT_P(NTH,NP,NR,NSA) = DWTT(NTH,NP,NR,NSA)
!              IC, EC is not updated
            END DO
         END DO
      END DO
      END DO
      END IF ! N_IMPL=0

      IF(N_IMPL.ne.0)THEN ! N_IMPL!=0
      DO NSA=NSASTART, NSAEND
      DO NR=NRSTART, NREND
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DWPP(NTH,NP,NR,NSA) = DWPP_P(NTH,NP,NR,NSA)
               DWPT(NTH,NP,NR,NSA) = DWPT_P(NTH,NP,NR,NSA)
               DWICPP(NTH,NP,NR,NSA) = DWICPP_P(NTH,NP,NR,NSA)
               DWICPT(NTH,NP,NR,NSA) = DWICPT_P(NTH,NP,NR,NSA)
               DWECPP(NTH,NP,NR,NSA) = DWECPP_P(NTH,NP,NR,NSA)
               DWECPT(NTH,NP,NR,NSA) = DWECPT_P(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DWTP(NTH,NP,NR,NSA) = DWTP_P(NTH,NP,NR,NSA)
               DWTT(NTH,NP,NR,NSA) = DWTT_P(NTH,NP,NR,NSA)
!              IC, EC is not updated
            END DO
         END DO
      END DO
      END DO
      END IF ! N_IMPL!=0
!!!   end of the reduction of the number of calling fp_calwm and fp_calw

      IF(N_IMPL.ne.0) CALL FPWAVE_CONST
!     ----- Constant Dw
      NCONST_RF=3

      DO NSA=NSASTART,NSAEND

! TOTAL Pabs(r) invariant
      IF(MODELW(NSA).eq.4.and.NCONST_RF.eq.2.and.N_IMPL.ne.0)THEN 
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  IF(RPW_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWPP(NTH,NP,NR,NSA)=DWPP(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA) &
                                        /RPW_IMPL(NR,NSA,N_IMPL)
                     DWPT(NTH,NP,NR,NSA)=DWPT(NTH,NP,NR,NSA)*RPW_INIT(NR,NSA) &
                                        /RPW_IMPL(NR,NSA,N_IMPL)
                     DWECPP(NTH,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA) &
                                          *RPW_INIT(NR,NSA) &
                                          /RPW_IMPL(NR,NSA,N_IMPL)
                     DWECPT(NTH,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA) &
                                          *RPW_INIT(NR,NSA) &
                                          /RPW_IMPL(NR,NSA,N_IMPL)
                  END IF
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  IF(RPW_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWTP(NTH,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA) &
                                        *RPW_INIT(NR,NSA) &
                                        /RPW_IMPL(NR,NSA,N_IMPL)
                     DWTT(NTH,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA) &
                                        *RPW_INIT(NR,NSA) &
                                        /RPW_IMPL(NR,NSA,N_IMPL)
                  END IF
               END DO
            END DO
         END DO
! Pabs_EC(r), Pabs_IC(r) invariant
      ELSEIF(MODELW(NSA).eq.4.and.NCONST_RF.eq.3.and.N_IMPL.ne.0)THEN 
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  IF(RPWEC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWECPP(NTH,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA) &
                                  *RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL)
                     DWECPT(NTH,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA) &
                                  *RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL)
                  END IF
                  IF(RPWIC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWICPP(NTH,NP,NR,NSA)=DWICPP(NTH,NP,NR,NSA) &
                                 *RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL)
                     DWICPT(NTH,NP,NR,NSA)=DWICPT(NTH,NP,NR,NSA) &
                                 *RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL)
                  END IF
                  DWPP(NTH,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA) &
                                     +DWICPP(NTH,NP,NR,NSA)
                  DWPT(NTH,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA) &
                                     +DWICPT(NTH,NP,NR,NSA)
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  DWTPEC=0.D0
                  DWTPIC=0.D0
                  DWTTEC=0.D0
                  DWTTIC=0.D0
                  IF(RPWEC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWTPEC = DWECTP(NTH,NP,NR,NSA) &
                                *RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL) 
                     DWTTEC = DWECTT(NTH,NP,NR,NSA) &
                                *RPWEC_INIT(NR,NSA)/RPWEC_IMPL(NR,NSA,N_IMPL) 
                  END IF
                  IF(RPWIC_IMPL(NR,NSA,N_IMPL).gt.0.D0)THEN
                     DWTPIC =( DWTP(NTH,NP,NR,NSA)-DWECTP(NTH,NP,NR,NSA) ) &
                          *RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL) 
                     DWTTIC =( DWTT(NTH,NP,NR,NSA)-DWECTT(NTH,NP,NR,NSA) ) &
                          *RPWIC_INIT(NR,NSA)/RPWIC_IMPL(NR,NSA,N_IMPL) 
                  END IF
                  DWTP(NTH,NP,NR,NSA)=DWTPEC+DWTPIC
                  DWTT(NTH,NP,NR,NSA)=DWTTEC+DWTTIC
               END DO
            END DO
         END DO
      END IF
   END DO
 END SUBROUTINE FP_CALW

      SUBROUTINE FP_CALW0
!
      USE plprof, only: rsrhon
      IMPLICIT NONE
      integer:: NSA, NR, NP, NTH
      real(8):: DLHA, DFWA, DECA, DECB, DECC, DLHL, DFWL, DECL
      real(8):: FACT
!
! =============  CALCULATION OF DWPP AND DWPT  ===============
!
      FACT=0.5D0
      DO NSA=NSASTART,NSAEND
      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
!            DO NP=1,NPMAX+1
            DO NP=NPSTART,NPENDWG
               CALL FPSUMW(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP,NSA),NR, &
                           DLHA,DFWA,DECA,DECB,DECC,NSA)
               DWPP(NTH,NP,NR,NSA)  =ABS(COSM(NTH))  &
                                     *(DLHA+DFWA+SINM(NTH)**2*DECA)
               DWLHPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*DLHA
               DWFWPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*DFWA
               DWECPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*SINM(NTH)**2*DECA
               IF(NTH.LE.NTHMAX/2) THEN
                  DWPT(NTH,NP,NR,NSA)  =-SINM(NTH)*(DLHA+DFWA-DECB)
                  DWLHPT(NTH,NP,NR,NSA)=-SINM(NTH)*DLHA
                  DWFWPT(NTH,NP,NR,NSA)=-SINM(NTH)*DFWA
                  DWECPT(NTH,NP,NR,NSA)= SINM(NTH)*DECB
               ELSE
                  DWPT(NTH,NP,NR,NSA)  = SINM(NTH)*(DLHA+DFWA-DECB)
                  DWLHPT(NTH,NP,NR,NSA)= SINM(NTH)*DLHA
                  DWFWPT(NTH,NP,NR,NSA)= SINM(NTH)*DFWA
                  DWECPT(NTH,NP,NR,NSA)=-SINM(NTH)*DECB
               ENDIF
            ENDDO
  101       CONTINUE
         ENDDO
!
         IF(MODELA.EQ.1) THEN
!            DO NP=1,NPMAX+1
            DO NP=NPSTART,NPENDWG
               DO NTH=ITL(NR)+1,NTHMAX/2
                  DWPP(NTH,NP,NR,NSA)  =(DWPP(NTH,NP,NR,NSA) &
                                  +DWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWPT(NTH,NP,NR,NSA)  =(DWPT(NTH,NP,NR,NSA) &
                                  +DWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWLHPP(NTH,NP,NR,NSA)=(DWLHPP(NTH,NP,NR,NSA) &
                                  +DWLHPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWLHPT(NTH,NP,NR,NSA)=(DWLHPT(NTH,NP,NR,NSA) &
                                  +DWLHPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWFWPP(NTH,NP,NR,NSA)=(DWFWPP(NTH,NP,NR,NSA) &
                                  +DWFWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWFWPT(NTH,NP,NR,NSA)=(DWFWPT(NTH,NP,NR,NSA) &
                                  +DWFWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWECPP(NTH,NP,NR,NSA)=(DWECPP(NTH,NP,NR,NSA) &
                                  +DWECPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWECPT(NTH,NP,NR,NSA)=(DWECPT(NTH,NP,NR,NSA) &
                                  +DWECPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWPP(NTH,NP,NR,NSA)
                  DWPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWPT(NTH,NP,NR,NSA)
                  DWLHPP(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPP(NTH,NP,NR,NSA)
                  DWLHPT(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPT(NTH,NP,NR,NSA)
                  DWFWPP(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPP(NTH,NP,NR,NSA)
                  DWFWPT(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPT(NTH,NP,NR,NSA)
                  DWECPP(NTHMAX-NTH+1,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)
                  DWECPT(NTHMAX-NTH+1,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)
               ENDDO
               DWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0     &
                        *( DWPP(ITL(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWPP(ITL(NR)+1,NP,NR,NSA)               &
                                            /RLAMDA(ITL(NR)+1,NR)  &
                          +DWPP(ITU(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWPP(ITU(NR)+1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)+1,NR)) 
               DWLHPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWLHPP(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWLHPP(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWLHPP(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWLHPP(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWFWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWFWPP(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWFWPP(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWFWPP(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWFWPP(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWECPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWECPP(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWECPP(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWECPP(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWECPP(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
!
               DWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0     &
                        *( DWPT(ITL(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWPT(ITL(NR)+1,NP,NR,NSA)               &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWPT(ITU(NR)-1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWPT(ITU(NR)+1,NP,NR,NSA)               &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWLHPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWLHPT(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWLHPT(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWLHPT(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWLHPT(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWFWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWFWPT(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWFWPT(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWFWPT(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWFWPT(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWECPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0   &
                        *( DWECPT(ITL(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)-1,NR) &
                          +DWECPT(ITL(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITL(NR)+1,NR) &
                          +DWECPT(ITU(NR)-1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)-1,NR) &
                          +DWECPT(ITU(NR)+1,NP,NR,NSA)             &
                                             /RLAMDA(ITU(NR)+1,NR))
               DWPP(ITU(NR),NP,NR,NSA)  =DWPP(ITL(NR),NP,NR,NSA)
               DWPT(ITU(NR),NP,NR,NSA)  =DWPT(ITL(NR),NP,NR,NSA)
               DWLHPP(ITU(NR),NP,NR,NSA)=DWLHPP(ITL(NR),NP,NR,NSA)
               DWLHPT(ITU(NR),NP,NR,NSA)=DWLHPT(ITL(NR),NP,NR,NSA)
               DWFWPP(ITU(NR),NP,NR,NSA)=DWFWPP(ITL(NR),NP,NR,NSA)
               DWFWPT(ITU(NR),NP,NR,NSA)=DWFWPT(ITL(NR),NP,NR,NSA)
               DWECPP(ITU(NR),NP,NR,NSA)=DWECPP(ITL(NR),NP,NR,NSA)
               DWECPT(ITU(NR),NP,NR,NSA)=DWECPT(ITL(NR),NP,NR,NSA)

            END DO
         ENDIF
      END DO
      END DO
!
! =============  CALCULATION OF DWTP AND DWTT  ===============
!
      DO NSA=NSASTART,NSAEND
      DO NR=NRSTART,NREND
!
         DO NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  CALL FPSUMW(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP,NSA), &
                              NR,DLHA,DFWA,DECA,DECB,DECC,NSA)
!
                  IF(NTH.LE.NTHMAX/2) THEN
                     DWTP(NTH,NP,NR,NSA)=-SING(NTH)*(DLHA+DFWA-DECB)
                  ELSE
                     DWTP(NTH,NP,NR,NSA)= SING(NTH)*(DLHA+DFWA-DECB)
                  ENDIF
                  DWTT(NTH,NP,NR,NSA)=(SING(NTH)**2*(DLHA+DFWA)+DECC) &
                                     /ABS(COSG(NTH))
               ENDDO
            ENDIF
         ENDDO
!
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            CALL FPWAVE(PM(NP,NSA),0.D0,NR,RSRHON(RM(NR)),0.0D0, &
                        DLHL,DFWL,DECL,NSA)
            DWTP(NTHMAX/2+1,NP,NR,NSA)=0.D0
            DWTT(NTHMAX/2+1,NP,NR,NSA)=DLHL+DFWL
         ENDDO
!
         IF(MODELA.EQ.1) THEN
            DO NTH=ITL(NR)+1,NTHMAX/2
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  DWTP(NTH,NP,NR,NSA)=(DWTP(NTH,NP,NR,NSA) &
                                 -DWTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWTT(NTH,NP,NR,NSA)=(DWTT(NTH,NP,NR,NSA) &
                                 +DWTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWTP(NTHMAX-NTH+2,NP,NR,NSA)=-DWTP(NTH,NP,NR,NSA)
                  DWTT(NTHMAX-NTH+2,NP,NR,NSA)= DWTT(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE FP_CALW0

!-------------------------------------------------------------

      SUBROUTINE FPWAVE_CONST
!
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NPS
      integer:: IERR
      real(kind8):: RSUM_W,RSUM_EC,RSUM_IC
      real(kind8):: PV, WPL, WPM, WPP
      real(kind8):: DFP, DFT, FFP, FACT

      CALL mtx_set_communicator(comm_np) 
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)

            RSUM_W=0.D0
            RSUM_EC=0.D0
            RSUM_IC=0.D0

            IF(NPSTART.eq.1)THEN
               NPS=2
            ELSE
               NPS=NPSTART
            END IF
!            DO NP=2,NPMAX
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
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNSP(NTH,NP  ,NR,NSBA)     &
                                  +WPM *FNSP(NTH,NP-1,NR,NSBA))&
                          )

                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT=    1.D0/DELTH                         & 
                         *(-                                    &
                            ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSBA))&
                          +                                     &
                            ((1.D0-WPP)*FNSP(NTH,NP  ,NR,NSBA)     &
                                  +WPP *FNSP(NTH,NP-1,NR,NSBA))&
                          )
                  ELSE
                     DFT=    1.D0/(2.D0*DELTH)                  &
                         *(                                     &
                            ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSBA)   &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSBA))&
                           -                                    &
                            ((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSBA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSBA))&
                                  )
                  ENDIF

                  RSUM_W = RSUM_W+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                         *(DWPP(NTH,NP,NR,NSA)*DFP           &
                          +DWPT(NTH,NP,NR,NSA)*DFT)
                  RSUM_IC = RSUM_IC+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                         *(DWICPP(NTH,NP,NR,NSA)*DFP         &
                          +DWICPT(NTH,NP,NR,NSA)*DFT)
                  RSUM_EC = RSUM_EC+PG(NP,NSBA)**2*SINM(NTH)/PV   &
                         *(DWECPP(NTH,NP,NR,NSA)*DFP         &
                          +DWECPT(NTH,NP,NR,NSA)*DFT)
               ENDDO
            ENDDO
            CALL p_theta_integration(RSUM_W)
            CALL p_theta_integration(RSUM_IC)
            CALL p_theta_integration(RSUM_EC)
               
            FACT=RNFP0(NSA)*1.D20*PTFP0(NSA)**2/AMFP(NSA)
            RPW_IMPL(NR,NSA,N_IMPL)=-RSUM_W*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6 
            RPWIC_IMPL(NR,NSA,N_IMPL)=-RSUM_IC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            RPWEC_IMPL(NR,NSA,N_IMPL)=-RSUM_EC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            IF(N_IMPL.eq.0)THEN
!               WRITE(6,'("ALERT ", 3I4)') NR, NSA, N_IMPL
               RPW_INIT(NR,NSA)=-RSUM_W*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6 
               RPWIC_INIT(NR,NSA)=-RSUM_IC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
               RPWEC_INIT(NR,NSA)=-RSUM_EC*FACT*2.D0*PI*DELP(NSBA)*DELTH *1.D-6
            END IF
         ENDDO
      ENDDO

      CALL mtx_reset_communicator
      RETURN
      END SUBROUTINE FPWAVE_CONST

!
! =======================================================
!
      SUBROUTINE FPSUMW(ETA,RSIN,RCOS,P,NR,SUM1,SUM2,SUM3,SUM4,SUM5,NSA)
!
      USE plprof, only: rsrhon
      IMPLICIT NONE
      integer:: NR, NSA, N
      real(8):: ETA, RSIN, RCOS, P, SUM1, SUM2, SUM3, SUM4, SUM5
      real(8):: DELH, ETAL, X, Y, PSI, PSIN, PCOS, PPERP, PPARA
      real(8):: DLHL, DFWL, DECL, XM
!
      DELH=2.D0*ETA/NAVMAX
!
      SUM1=0.D0
      SUM2=0.D0
      SUM3=0.D0
      SUM4=0.D0
      SUM5=0.D0
!
      DO N=1,NAVMAX
         ETAL=DELH*(N-0.5D0)
         X=RSRHON(RM(NR))*COS(ETAL)
         Y=RSRHON(RM(NR))*SIN(ETAL)
         XM=EPSRM(NR)*COS(ETAL)*RR
!
         IF(MODELA.EQ.0) THEN
            PSI=1.D0
            PSIN=RSIN
!            PCOS=ABS(RCOS)
            PCOS=RCOS
         ELSE
            PSI=(1.D0+EPSRM(NR))/(1.D0+XM/RR)
            PSIN=SQRT(PSI)*RSIN
!            PCOS=SQRT(1.D0-PSI*RSIN**2)
!            write(6,'(A,I5,1P5E12.4)') 'NR:',NR,&
!                 EPSRM(NR),X/RR,PSI,RSIN,1.D0-PSI*RSIN**2
            IF (RCOS.GT.0.0D0) THEN
               PCOS= SQRT(1.D0-PSI*RSIN**2)
            ELSE
               PCOS=-SQRT(1.D0-PSI*RSIN**2)
            END IF
         ENDIF
!
         PPERP=P*PSIN
!         IF (RCOS.GT.0.0) THEN
            PPARA =  P*PCOS
!         ELSE
!            PPARA = -P*PCOS
!         END IF
!
         CALL FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL,NSA)
!
         SUM1=SUM1+PCOS       *DLHL
         SUM2=SUM2+PCOS       *DFWL
!         SUM3=SUM3+PSI/PCOS   *DECL
!         SUM4=SUM4+PCOS       *DECL
!         SUM5=SUM5+PCOS**3/PSI*DECL
         SUM3=SUM3+DECL*RCOS/PCOS
         SUM4=SUM4+DECL/SQRT(PSI)
         SUM5=SUM5+DECL*PCOS/RCOS/PSI
!!! unverified for LH, FW
      END DO
      SUM1=SUM1*DELH
      SUM2=SUM2*DELH
      SUM3=SUM3*DELH
      SUM4=SUM4*DELH
      SUM5=SUM5*DELH

      RETURN
      END SUBROUTINE FPSUMW
!
!     ******************************************
!        QUASI-LINEAR WAVE DIFFUSION COEF.
!     ******************************************
!
      SUBROUTINE FPWAVE(PPARA,PPERP,NR,X,Y,DLHL,DFWL,DECL,NSA)
!
      IMPLICIT NONE
      integer:: NR, NSA, NSB
      real(8):: PPARA, PPERP, X, Y, DLHL, DFWL, DECL
      real(8):: P2, PVPARA, RNUDL, RNUFL, AMI, AEI, WPI2, FACT, FACT2
      real(8):: DFWL1, DFWL2, ARG, ARG1, FACT1, W, PARAN, FN, DELF, ARG2
      real(8):: WFW2, ARG3, FACT3
!
      P2=PPARA**2+PPERP**2
      PVPARA=PPARA/SQRT(1.D0+P2*THETA0(NSA))
      RNUDL=0.D0 
      DO NSB=1,NSBMAX
         IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
            RNUDL=RNUD(NR,NSB,NSA)
            RNUFL=RNUF(NR,NSB,NSA)
         ENDIF
      ENDDO
!
!        ***** LOWER HTYBRID WAVE *****
!
!        IF (PVPARA.GT.PLH1.AND.PVPARA.LT.PLH2.AND.
!    &       RM(NR).GE.RLH)THEN
!           DLHL=DLH*RNUD(NR,NSA)
!        ELSE
!           DLHL=0.D0
!        END IF
!
         IF (DLH.GT.0.D0.AND.RM(NR).GE.RLH)THEN
            IF(ABS(PVPARA).GT.1.D-70) THEN
               ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-PLH1)/PLH2)**2
               IF(ARG.LT.20.D0) THEN
                  DLHL=DLH*RNUDL+EXP(-ARG)
               ELSE
                  DLHL=0.D0
               ENDIF
            ELSE
               DLHL=0.D0
            ENDIF
         ELSE
            DLHL=0.D0
         END IF
!
!        ***** FAST WAVE *****
!
!        IF (PVPARA.GT.PFW1.AND.PVPARA.LT.PFW2.AND.
!    &       RM(NR).GE.RFW)THEN
!           DFWL=DFW*RNUD(NR,NSA)/(PPARA**2)
!    &         *(PPERP**2-2.D0*(TE(NR)+0.025D0*551.D0/TE0))**2
!        ELSE
!           DFWL=0.D0
!        END IF
!
         IF (DFW.GT.0.D0.AND.RM(NR).GE.RFW)THEN
            IF(ABS(PVPARA).GT.1.D-70) THEN
               AMI=2*AMP
               AEI=1*AEE
               WPI2=RNFP0(NSA)*1.D20*RNFP(NR,NSA)*AEI*AEE/(AMI*EPS0)
               WFW2=(RFW*1.D6*2*PI)**2
               FACT=RTFP(NR,NSA)+WFW2*AMFP(NSA)*VC*VC &
                                   /(WPI2*RTFP0(NSA)*1.D3)
               IF(ABS(PFW1-2.2D0).LT.1.D-30) THEN
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=DFW*RNUDL*EXP(-ARG)/(PVPARA**2) &
                         *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)+2.2D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=DFW*RNUFL*EXP(-ARG)/(PVPARA**2) &
                         *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSEIF(ABS(PFW1-1.6D0).LT.1.D-30) THEN
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-1.6D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL1=1.5D0*DFW*RNUDL*EXP(-ARG)   &
                          /(PVPARA**2)                 &
                         *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL1=0.D0
                  ENDIF
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)+2.8D0)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL2=0.5D0*DFW*RNUDL*EXP(-ARG)   &
                          /(PVPARA**2)                 &
                          *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL2=0.D0
                  ENDIF
                  DFWL=DFWL1+DFWL2
               ELSE
                  ARG=((AMFP(NSA)*VC/(PTFP0(NSA)*PVPARA)-PFW1)/PFW2)**2
                  IF(ARG.LT.20.D0) THEN
                     DFWL=DFW*RNUDL*EXP(-ARG)          &
                          /(PVPARA**2)                 &
                          *(PPERP**2-2.D0*FACT)**2
                  ELSE
                     DFWL=0.D0
                  ENDIF
               ENDIF
            ELSE
               DFWL=0.D0
            ENDIF
         ELSE
            DFWL=0.D0
         END IF
!
!        ***** ELECTRON CYCLOTRON WAVE *****
!
         IF(DEC.GT.0.D0.AND.ABS(PPARA).GT.1.D-70) THEN
            ARG1=(Y/DELYEC)**2
            IF(ARG1.LT.20.0) THEN
               FACT1=EXP(-ARG1)
            ELSE
               FACT1=0.D0
            END IF

            W=RFEC/(1.D0+X/RR)
            PARAN=PPARA*SQRT(RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC))
            FN=PEC1*PARAN-SQRT(1.D0+THETA0(NSA)*P2)+W
            DELF=PEC2*PARAN
            ARG2=(FN/DELF)**2
            IF(ARG2.LT.20.0) THEN
               FACT2=EXP(-ARG2)
            ELSE
               FACT2=0.D0
            END IF

            ARG3=(X-PEC3)/PEC4
            IF(ARG3.GT.20.D0) THEN
               FACT3=1.D0
            ELSE IF(ARG3.LT.-20.D0) THEN
               FACT3=0.D0
            ELSE
               FACT3=0.5D0*(1.D0+TANH(ARG3))
            ENDIF
              
            DECL=DEC*RNUDL*FACT1*FACT2*FACT3
         ELSE
            DECL=0.D0
         ENDIF
!
      RETURN
      END SUBROUTINE FPWAVE
      END MODULE fpcalw
      
