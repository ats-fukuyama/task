! fpcalw.f90

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
    real(rkind):: FACTOR,PABSX_LH,PABSX_FW,PABSX_EC,PABSX_WR,PABSX_WM

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         IF(MODELW(NS).EQ.0) THEN
            CALL FP_CALW0(NSA)
         ELSEIF(MODELW(NS).EQ.1) THEN
            CALL FP_CALWR(NSA)
         ELSEIF(MODELW(NS).EQ.2) THEN
            CALL FP_CALWR(NSA)
         ELSEIF(MODELW(NS).EQ.3) THEN
            CALL FP_CALWM(NSA)
         ELSEIF(MODELW(NS).EQ.4) THEN
            CALL FP_CALWM(NSA)
         ELSE
            IF(nrank.eq.0) WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW(NS)
         ENDIF
      END DO

!     ----- Adjust DW

      IF(PABS_LH.NE.0.D0) THEN
         CALL FP_CALPABS(DWLHPP,DWLHPT,PABSX_LH)
         IF(PABSX_LH.NE.0.D0) THEN
            FACTOR=PABS_LH/PABSX_LH
            DWLHPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWLHPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWLHPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWLHPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
         END IF
      END IF

      IF(PABS_FW.NE.0.D0) THEN
         CALL FP_CALPABS(DWFWPP,DWFWPT,PABSX_FW)
         IF(PABSX_FW.NE.0.D0) THEN
            FACTOR=PABS_FW/PABSX_FW
            DWFWPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWFWPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWFWPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWFWPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
         END IF
      END IF

      IF(PABS_EC.NE.0.D0) THEN
         CALL FP_CALPABS(DWECPP,DWECPT,PABSX_EC)
         IF(PABSX_EC.NE.0.D0) THEN
            FACTOR=PABS_EC/PABSX_EC
            DWECPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWECPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWECPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWECPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWECTP(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWECTP(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND)
            DWECTT(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWECTT(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND)
         END IF
      END IF

      IF(PABS_WR.NE.0.D0) THEN
         CALL FP_CALPABS(DWWRPP,DWWRPT,PABSX_WR)
         IF(PABSX_WR.NE.0.D0) THEN
            FACTOR=PABS_WR/PABSX_WR
            DWWRPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWRPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWWRPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWRPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWWRTP(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWRTP(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND)
            DWWRTT(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWRTT(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND)
         END IF
      END IF

      IF(PABS_WM.NE.0.D0) THEN
         CALL FP_CALPABS(DWWMPP,DWWMPT,PABSX_WM)
         IF(PABSX_WM.NE.0.D0) THEN
            FACTOR=PABS_WM/PABSX_WM
            DWWMPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWMPP(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWWMPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWMPT(1:NTHMAX,NPSTART:NPENDWG,NRSTART:NREND,NSASTART:NSAEND)
            DWWMTP(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWMTP(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND)
            DWWMTT(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND) &
            =FACTOR* &
            DWWMTT(1:NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NREND,NSASTART:NSAEND)
         END IF
      END IF

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  DWPP(NTH,NP,NR,NSA)=DWLHPP(NTH,NP,NR,NSA) &
                                     +DWFWPP(NTH,NP,NR,NSA) &
                                     +DWECPP(NTH,NP,NR,NSA) &
                                     +DWWRPP(NTH,NP,NR,NSA) &
                                     +DWWMPP(NTH,NP,NR,NSA)
                  DWPT(NTH,NP,NR,NSA)=DWLHPT(NTH,NP,NR,NSA) &
                                     +DWFWPT(NTH,NP,NR,NSA) &
                                     +DWECPT(NTH,NP,NR,NSA) &
                                     +DWWRPT(NTH,NP,NR,NSA) &
                                     +DWWMPT(NTH,NP,NR,NSA)
               END DO
            END DO
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  DWTP(NTH,NP,NR,NSA)=DWECTP(NTH,NP,NR,NSA) &
                                     +DWWRTP(NTH,NP,NR,NSA) &
                                     +DWWMTP(NTH,NP,NR,NSA)
                  DWTT(NTH,NP,NR,NSA)=DWECTT(NTH,NP,NR,NSA) &
                                     +DWWRTT(NTH,NP,NR,NSA) &
                                     +DWWMTT(NTH,NP,NR,NSA)
            END DO
         END DO
      END DO
   END DO

   RETURN

 END SUBROUTINE FP_CALW

 SUBROUTINE FP_CALW0(NSA)
!
      USE plprof, only: rsrhon
      IMPLICIT NONE
      integer:: NSA, NR, NP, NTH, NS
      REAL(rkind):: DLHA, DFWA, DECA, DECB, DECC, DLHL, DFWL, DECL
      REAL(rkind):: FACT
!
! =============  CALCULATION OF DWPP AND DWPT  ===============
!
      FACT=0.5D0
      NS=NS_NSA(NSA)

      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
!            DO NP=1,NPMAX+1
            DO NP=NPSTART,NPENDWG
               CALL FPSUMW(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP,NS),NR, &
                           DLHA,DFWA,DECA,DECB,DECC,NSA)
               DWLHPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*DLHA
               DWFWPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*DFWA
               DWECPP(NTH,NP,NR,NSA)=ABS(COSM(NTH))*SINM(NTH)**2*DECA
               IF(NTH.LE.NTHMAX/2) THEN
                  DWLHPT(NTH,NP,NR,NSA)=-SINM(NTH)*DLHA
                  DWFWPT(NTH,NP,NR,NSA)=-SINM(NTH)*DFWA
                  DWECPT(NTH,NP,NR,NSA)= SINM(NTH)*DECB
               ELSE
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
                  DWLHPP(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPP(NTH,NP,NR,NSA)
                  DWLHPT(NTHMAX-NTH+1,NP,NR,NSA)=DWLHPT(NTH,NP,NR,NSA)
                  DWFWPP(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPP(NTH,NP,NR,NSA)
                  DWFWPT(NTHMAX-NTH+1,NP,NR,NSA)=DWFWPT(NTH,NP,NR,NSA)
                  DWECPP(NTHMAX-NTH+1,NP,NR,NSA)=DWECPP(NTH,NP,NR,NSA)
                  DWECPT(NTHMAX-NTH+1,NP,NR,NSA)=DWECPT(NTH,NP,NR,NSA)
               ENDDO
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
               DWLHPP(ITU(NR),NP,NR,NSA)=DWLHPP(ITL(NR),NP,NR,NSA)
               DWLHPT(ITU(NR),NP,NR,NSA)=DWLHPT(ITL(NR),NP,NR,NSA)
               DWFWPP(ITU(NR),NP,NR,NSA)=DWFWPP(ITL(NR),NP,NR,NSA)
               DWFWPT(ITU(NR),NP,NR,NSA)=DWFWPT(ITL(NR),NP,NR,NSA)
               DWECPP(ITU(NR),NP,NR,NSA)=DWECPP(ITL(NR),NP,NR,NSA)
               DWECPT(ITU(NR),NP,NR,NSA)=DWECPT(ITL(NR),NP,NR,NSA)

            END DO
         ENDIF
      END DO
!
! =============  CALCULATION OF DWTP AND DWTT  ===============
!
      DO NR=NRSTART,NREND
!
         DO NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  CALL FPSUMW(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP,NS), &
                              NR,DLHA,DFWA,DECA,DECB,DECC,NSA)
!
                  IF(NTH.LE.NTHMAX/2) THEN
                     DWECTP(NTH,NP,NR,NSA)=+SING(NTH)*DECB
                  ELSE
                     DWECTP(NTH,NP,NR,NSA)=-SING(NTH)*DECB
                  ENDIF
                  DWECTT(NTH,NP,NR,NSA)=             DECC/ABS(COSG(NTH))
               ENDDO
            ENDIF
         ENDDO
!
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            CALL FPWAVE(PM(NP,NSA),0.D0,NR,RSRHON(RM(NR)),0.0D0, &
                        DLHL,DFWL,DECL,NSA)
            DWECTP(NTHMAX/2+1,NP,NR,NSA)=0.D0
            DWECTT(NTHMAX/2+1,NP,NR,NSA)=0.D0
         ENDDO
!
         IF(MODELA.EQ.1) THEN
            DO NTH=ITL(NR)+1,NTHMAX/2
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  DWECTP(NTH,NP,NR,NSA)=(DWECTP(NTH,NP,NR,NSA) &
                                        -DWECTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWECTT(NTH,NP,NR,NSA)=(DWECTT(NTH,NP,NR,NSA) &
                                        +DWECTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWECTP(NTHMAX-NTH+2,NP,NR,NSA)=-DWECTP(NTH,NP,NR,NSA)
                  DWECTT(NTHMAX-NTH+2,NP,NR,NSA)= DWECTT(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE FP_CALW0

!-------------------------------------------------------------

      SUBROUTINE FP_CALPABS(DQPP,DQPT,PABSX)
!
      USE fpmpi
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: &
           DQPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX), &
           DQPT(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX)
      REAL(rkind),INTENT(OUT):: PABSX
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NPS, N, NSW
      integer:: IERR
      real(rkind):: RSUML,RSUM(NRSTART:NREND,NSAMAX),RSUMG(NRMAX,NSAMAX)
      real(rkind):: PV, WPL, WPM, WPP
      real(rkind):: DFP, DFT, FFP, FACTOR

!----- sum up over NTH and NP

      CALL mtx_set_communicator(comm_np) 
      DO NR=NRSTART,NREND
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)

            RSUML=0.D0

            IF(NPSTART.eq.1)THEN
               NPS=2
            ELSE
               NPS=NPSTART
            END IF
            DO NP=NPS,NPEND
               PV=SQRT(1.D0+THETA0(NS)*PG(NP,NS)**2)
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
                  DFP=(FNSP(NTH,NP,NR,NSA)-FNSP(NTH,NP-1,NR,NSA)) &
                      *PG(NP,NS)/DELP(NS)
                  IF(NTH.EQ.1) THEN
                     DFT= 1.D0/DELTH                              &
                         *( ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSA)  &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSA)) &
                           -((1.D0-WPM)*FNSP(NTH  ,NP  ,NR,NSA)    &
                                  +WPM *FNSP(NTH  ,NP-1,NR,NSA)))  
                  ELSE IF(NTH.EQ.NTHMAX) THEN
                     DFT= 1.D0/DELTH                               & 
                         *(-((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSA))  &
                          + ((1.D0-WPP)*FNSP(NTH  ,NP  ,NR,NSA)   &
                                  +WPP *FNSP(NTH  ,NP-1,NR,NSA)))
                  ELSE
                     DFT= 1.D0/(2.D0*DELTH)                        &
                         *( ((1.D0-WPP)*FNSP(NTH+1,NP  ,NR,NSA)   &
                                  +WPP *FNSP(NTH+1,NP-1,NR,NSA))  &
                           -((1.D0-WPM)*FNSP(NTH-1,NP  ,NR,NSA)   &
                                  +WPM *FNSP(NTH-1,NP-1,NR,NSA)))
                  ENDIF

                  RSUML = RSUML+PG(NP,NS)**2*SINM(NTH)/PV &
                                *(DQPP(NTH,NP,NR,NSA)*DFP    &
                                 +DQPT(NTH,NP,NR,NSA)*DFT)
               ENDDO
            ENDDO
            CALL p_theta_integration(RSUML)
               
            FACTOR=-RNFP0(NSA)*1.D20*PTFP0(NSA)**2*RFSADG(NR) &
                   *2.D0*PI*DELP(NS)*DELTH*1.D-6/AMFP(NSA)
            RSUM(NR,NSA)=RSUML*FACTOR
         ENDDO
      ENDDO
      CALL mtx_reset_communicator

!----- sum up over NR and NSA

      CALL mtx_set_communicator(comm_nsanr) 
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RSUM,SAVLEN(NRANK+1),RSUMG,N,NSA)
      END DO
      IF(nrank.EQ.0) THEN
         PABSX=0.D0
         DO NSA=1,NSAMAX
            DO NR=1,NRMAX
               PABSX=PABSX+RSUMG(NR,NSA)*VOLR(NR)
            END DO
         END DO
      END IF
      CALL mtx_reset_communicator 

!----- broadcast over world

      CALL mtx_broadcast1_real8(PABSX)

      RETURN
      END SUBROUTINE FP_CALPABS

!
! =======================================================
!
      SUBROUTINE FPSUMW(ETA,RSIN,RCOS,P,NR,SUM1,SUM2,SUM3,SUM4,SUM5,NSA)
!
      USE plprof, only: rsrhon
      IMPLICIT NONE
      integer:: NR, NSA, N
      REAL(rkind):: ETA, RSIN, RCOS, P, SUM1, SUM2, SUM3, SUM4, SUM5
      REAL(rkind):: DELH, ETAL, X, Y, PSI, PSIN, PCOS, PPERP, PPARA
      REAL(rkind):: DLHL, DFWL, DECL, XM
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
      integer:: NR, NSA, NSB, NS
      REAL(rkind):: PPARA, PPERP, X, Y, DLHL, DFWL, DECL
      REAL(rkind):: P2, PVPARA, RNUDL, RNUFL, AMI, AEI, WPI2, FACT, FACT2
      REAL(rkind):: DFWL1, DFWL2, ARG, ARG1, FACT1, W, PARAN, FN, DELF, ARG2
      REAL(rkind):: WFW2, ARG3, FACT3
!
      NS=NS_NSA(NSA)

      P2=PPARA**2+PPERP**2
      PVPARA=PPARA/SQRT(1.D0+P2*THETA0(NS))
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
            FN=PEC1*PARAN-SQRT(1.D0+THETA0(NS)*P2)+W
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
      
