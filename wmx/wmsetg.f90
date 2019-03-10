! wmsetg.f90

MODULE wmsetg

  PRIVATE
  PUBLIC wm_setg

CONTAINS

!     ****** RADIAL MESH AND METRIC TENSOR ******

  SUBROUTINE wm_setg(ierr)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    IERR=0

    IF(NTHMAX.EQ.1) THEN
       NTHMAX_F=1
    ELSE
       NTHMAX?F=NINT(NTHMAX*factor_nth)
    END IF
    IF(NHHMAX.EQ.1) THEN
       NHHMAX_F=1
    ELSE
       NHHMAX_F=NINT(NHHMAX*factor_nhh)
    END IF

    SELECT CASE(MODELG)

!     ****** Rectangular COORDINATES ******

    CASE(0)
!       CALL wmsetg_rect(ierr)
       WRITE(6,*) 'XX wmsetg: MODELG=0 is not supported'

!     ****** CYLINDRICAL COORDINATES ******

    CASE(1)
       CALL wmsetg_cyl(IERR)

!     ****** TOROIDAL COORDINATES ******

    CASE(2)
       CALL wmsetg_tor(IERR)

!     ****** EQUILIBRIUM (TASK/EQ) ******

    CASE(3)
       CALL wmsetg_eq(IERR)

!     ****** EQUILIBRIUM (VMEC) ******

    CASE(4)
       CALL wmsetg_vmec(IERR)

!     ****** EQUILIBRIUM (EQDSK) ******

    CASE(5)
       CALL wmsetg_eqdsk(IERR)

!     ****** EQUILIBRIUM (BOOZER) ******

    CASE(6)
       CALL wmsetg_boozer(IERR)
    END SELECT

    IF(IERR.NE.0) RETURN

! --- calculate inverse metric tensor RGI ---

    CALL wm_invert_rg

! --- load plasma profile data ---

    IF(MODELN.EQ.7) CALL WMDPRF(IERR)
    IF(MODELN.EQ.8) CALL WMXPRF(IERR)
    IF(MODELN.EQ.9) CALL WMTRLOAD(IERR)

    IF(IERR.NE.0) RETURN

! --- calculate mode range

    IF(NHHMAX.EQ.1) THEN
       NDSIZ  = 1
       NDMIN  = 0
       NDMAX  = 0
       KDSIZ  = 1
       KDMIN  = 0
       KDMAX  = 0
       NDSIZ_F= 1
       NDMIN_F= 0
       NDMAX_F= 0
       KDSIZ_F= 1
       KDMIN_F= 0
       KDMAX_F= 0
    ELSE
       NDSIZ  = NHHMAX
       NDMIN  =-NHHMAX/2+1
       NDMAX  = NHHMAX/2
       KDSIZ  = NHHMAX
       KDMIN  =-NHHMAX/2+1
       KDMAX  = NHHMAX/2
       NDSIZ_F= NHHMAX_F
       NDMIN_F=-NHHMAX_F/2+1
       NDMAX_F= NHHMAX_F/2
       KDSIZ_F= NHHMAX_F
       KDMIN_F=-NHHMAX_F/2+1
       KDMAX_F= NHHMAX_F/2
    ENDIF

    IF(NTHMAX.EQ.1) THEN
       MDSIZ  = 1
       MDMIN  = 0
       MDMAX  = 0
       LDSIZ  = 1
       LDMIN  = 0
       LDMAX  = 0
       MDSIZ_F= 1
       MDMIN_F= 0
       MDMAX_F= 0
       LDSIZ_F= 1
       LDMIN_F= 0
       LDMAX_F= 0
    ELSE
       MDSIZ  = NTHMAX
       MDMIN  =-NTHMAX/2+1
       MDMAX  = NTHMAX/2
       LDSIZ  = NTHMAX
       LDMIN  =-NTHMAX/2+1
       LDMAX  = NTHMAX/2
       MDSIZ_F= NTHMAX_F
       MDMIN_F=-NTHMAX_F/2+1
       MDMAX_F= NTHMAX_F/2
       LDSIZ_F= NTHMAX_F
       LDMIN_F=-NTHMAX_F/2+1
       LDMAX_F= NTHMAX_F/2
    ENDIF

    IF(MODELG.NE.6) CALL WMICRS
    CALL BTHOB
    IERR=0
    RETURN
  END SUBROUTINE wm_setg

!     ****** RADIAL MESH (CYLINDRICAL COORDINATES) ******

  SUBROUTINE wmsetg_cyl(IERR)

    USE wmcomm
    IMPLICIT NONE

    IERR=0

    NSUMAX=31
    NSWMAX=31
    NHHMAX=1
    NHHMAX_F=1

    PSIPA=RA*RA*BB/(Q0+QA)
    PSIPB=SQRT(RB**2/RA**2+(RB**2/RA**2-1.D0)*Q0/QA)*PSIPA
    DRHO=SQRT(PSIPB/PSIPA)/NRMAX

    DO NR=1,NRMAX+1
       RHOL=DRHO*(NR-1)
       XRHO(NR)=RHOL

       IF(NR.EQ.1) THEN
          XR(NR)=0.D0
       ELSEIF(RHOL.LT.1.D0) THEN
          XR(NR)=RA*SQRT((2.D0*Q0*RHOL**2+(QA-Q0)*RHOL**4)/(Q0+QA))
       ELSE
          XR(NR)=RA*SQRT((Q0+QA*RHOL**4)/(Q0+QA))
       ENDIF

!            WRITE(6,*) 'NR,RHO,XR=',NR,XRHO(NR),XR(NR)

       IF(RHOL.LT.1.D0) THEN
          QPS(NR)=Q0+(QA-Q0)*RHOL**2
       ELSE
          QPS(NR)=QA*RHOL**2
       ENDIF
    ENDDO

    DTH=2.D0*PI/NTHMAX_F
    DTHG=2.D0*PI/nthgmax
    DO NR=1,NRMAX+1
       IF(NR.EQ.1) THEN
          RS=XR(2)/XRHO(2)
       ELSE
          RS=XR(NR)/XRHO(NR)
       ENDIF
       RSD=QPS(NR)/(BB*RS)
       DO NTH=1,NTHMAX_F
          RCOS=COS(DTH*(NTH-1))
          RSIN=SIN(DTH*(NTH-1))
          RPS(NTH,NR)    = RR + XR(NR)*RCOS
          ZPS(NTH,NR)    =      XR(NR)*RSIN
          DRPSI(NTH,NR)  =      RSD   *RCOS
          DZPSI(NTH,NR)  =      RSD   *RSIN
          DRCHI(NTH,NR)  =     -RS    *RSIN
          DZCHI(NTH,NR)  =      RS    *RCOS
       ENDDO
       DO NTH=1,nthgmax
          RCOS=COS(DTHG*(NTH-1))
          RSIN=SIN(DTHG*(NTH-1))
          IF(MODELG.EQ.0) THEN
             RPSG(NTH,NR)  =      XR(NR)*RCOS
          ELSE IF(MODELG.EQ.1) THEN
             RPSG(NTH,NR)  = RR + XR(NR)*RCOS
          ENDIF
          ZPSG(NTH,NR)  =         XR(NR)*RSIN
       ENDDO
    ENDDO

    P0=0.D0
    DO NS=1,NSMAX
       P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
    ENDDO
    P0=P0*1.D20*AEE*1.D3/1.D6

    DO NR=1,NRMAX+1
       RHOL=XRHO(NR)
       IF(RHOL.LE.1.D0) THEN
          IF(PN(1).LE.0.D0) THEN
             FEDGE=0.D0
          ELSE
             FEDGE=PNS(1)/PN(1)
          END IF
          FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)+FEDGE
          PT=(PTPR(1)+2*PTPP(1))/3.D0
          FEDGE=PTS(1)/PT
          FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)+FEDGE
          PPS(NR)=P0*FACTN*FACTT
          PPS(NR)=0.D0
       ENDIF
       RBPS(NR)=BB*RR
       VPS(NR)=2*PI*RR*PI*XR(NR)**2
       RLEN(NR)=2*PI*XR(NR)
    ENDDO

    DTHU=2.D0*PI/(NSUMAX-1)
    DO NSU=1,NSUMAX
       RCOS=COS(DTHU*(NSU-1))
       RSIN=SIN(DTHU*(NSU-1))
       RSU(NSU,1)=RR+RA*RCOS
       RSW(NSU,1)=RR+RB*RCOS
       ZSU(NSU,1)=RA*RSIN
       ZSW(NSU,1)=RB*RSIN
    ENDDO

    RGMIN= RR-RB*1.01D0
    RGMAX= RR+RB*1.01D0
    ZGMIN=-RB*1.01D0
    ZGMAX= RB*1.01D0

    DO NR=1,NRMAX+1
       DO NTH=1,NTHMAX_F
          RRG=RPS(NTH,NR)
          DO NHH=1,NHHMAX_F
             RPST(NTH,NHH,NR)=RPS(NTH,NR)
             ZPST(NTH,NHH,NR)=ZPS(NTH,NR)

             RG11(NTH,NHH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
             RG12(NTH,NHH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR) &
                              +DZPSI(NTH,NR)*DZCHI(NTH,NR)
             RG13(NTH,NHH,NR)= 0.D0
             RG22(NTH,NHH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
             RG23(NTH,NHH,NR)= 0.D0
             RG33(NTH,NHH,NR)= RRG**2
             RJ  (NTH,NHH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR) &
                                    -DRCHI(NTH,NR)*DZPSI(NTH,NR))
             BFLD(2,NTH,NHH,NR)=1.D0/RJ(NTH,NHH,NR)
             BFLD(3,NTH,NHH,NR)=RBPS(NR)/RRG**2
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE wmsetg_cyl

!     ****** RADIAL MESH (TOROIDAL COORDINATES) ******

  SUBROUTINE wmsetg_tor(IERR)

!      USE plfile_prof_mod
      USE plprof,ONLY: pl_qprf
      INCLUDE 'wmcomm.inc'

!
         IERR=0

         NSUMAX=31
         NSWMAX=31
         NHHMAX=1
         NHHMAX_F=1
         NHHMAX_IPS_F=1

!         CALL plfile_prof_read(modeln,modelq,ierr)

         PSIPA=RA*RA*BB/(Q0+QA)
         DRHO=(RB/RA)/NRMAX

         DO NR=1,NRMAX+1
            RHOL=DRHO*(NR-1)
            XRHO(NR)=RHOL
            XR(NR)=RA*RHOL
            CALL PL_QPRF(RHOL,QPS(NR))
         ENDDO

         DTH=2.D0*PI/NTHMAX_F
         DO NR=1,NRMAX+1
            RSD=RA/(2.D0*PSIPA)
            DO NTH=1,NTHMAX_F
               RCOS=COS(DTH*(NTH-1))
               RSIN=SIN(DTH*(NTH-1))
               RPS(NTH,NR) = RR + XR(NR)*RCOS
               ZPS(NTH,NR)    =      XR(NR)*RSIN
               DRPSI(NTH,NR)  =      RSD   *RCOS
               DZPSI(NTH,NR)  =      RSD   *RSIN
               DRCHI(NTH,NR)  =     -RA    *RSIN
               DZCHI(NTH,NR)  =      RA    *RCOS
            ENDDO
         END DO

         DTHG=2.D0*PI/NTHGMAX
         DO NR=1,NRMAX+1
            DO NTH=1,NTHGMAX
               RCOS=COS(DTHG*(NTH-1))
               RSIN=SIN(DTHG*(NTH-1))
               IF(MODELG.EQ.0) THEN
                  RPSG(NTH,NR)  =      XR(NR)*RCOS
               ELSE
                  RPSG(NTH,NR)  = RR + XR(NR)*RCOS
               ENDIF
               ZPSG(NTH,NR)  =         XR(NR)*RSIN
            ENDDO
         ENDDO

         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6

         DO NR=1,NRMAX+1
            RHOL=XRHO(NR)
            IF(RHOL.LE.1.D0) THEN
               IF(PN(1).LE.0.D0) THEN
                  FEDGE=0.D0
               ELSE
                  FEDGE=PNS(1)/PN(1)
               END IF
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)+FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)+FEDGE
               PPS(NR)=P0*FACTN*FACTT
            ELSE
               PPS(NR)=0.D0
            ENDIF
            RBPS(NR)=BB*RR
            VPS(NR)=2*PI*RR*PI*XR(NR)**2
            RLEN(NR)=2*PI*XR(NR)
         ENDDO

         DTHU=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            RCOS=COS(DTHU*(NSU-1))
            RSIN=SIN(DTHU*(NSU-1))
            RSU(NSU,1)=RR+RA*RCOS
            RSW(NSU,1)=RR+RB*RCOS
            ZSU(NSU,1)=RA*RSIN
            ZSW(NSU,1)=RB*RSIN
         ENDDO

         RGMIN=RR-RB*1.01D0
         RGMAX=RR+RB*1.01D0
         ZGMIN=-RB*1.01D0
         ZGMAX= RB*1.01D0
         RAXIS=RR

         DO NR=1,NRMAX+1
            DO NTH=1,NTHMAX_F
               RRG=RPS(NTH,NR)
               DO NHH=1,NHHMAX_F
                  RPST(NTH,NHH,NR)=RPS(NTH,NR)
                  ZPST(NTH,NHH,NR)=ZPS(NTH,NR)

                  RG11(NTH,NHH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
                  RG12(NTH,NHH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR) &
                                   +DZPSI(NTH,NR)*DZCHI(NTH,NR)
                  RG13(NTH,NHH,NR)= 0.D0
                  RG22(NTH,NHH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
                  RG23(NTH,NHH,NR)= 0.D0
                  RG33(NTH,NHH,NR)= RRG**2
                  RJ  (NTH,NHH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR) &
                                         -DRCHI(NTH,NR)*DZPSI(NTH,NR))

                  BFLD(2,NTH,NHH,NR)=BB*RR/RRG**2/QPS(NR)
                  BFLD(3,NTH,NHH,NR)=BB*RR/RRG**2
               ENDDO
            ENDDO
         ENDDO

         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6

         DO NR=1,NRMAX+1
            RHOL=XRHO(NR)
            IF(RHOL.LE.1.D0) THEN
               IF(PN(1).LE.0.D0) THEN
                  FEDGE=0.D0
               ELSE
                  FEDGE=PNS(1)/PN(1)
               END IF
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)+FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)+FEDGE
               PPS(NR)=P0*FACTN*FACTT
            ELSE
               PPS(NR)=0.D0
            ENDIF
            RBPS(NR)=BB*RR
            VPS(NR)=2*PI*RR*PI*XR(NR)**2
            RLEN(NR)=2*PI*XR(NR)
         ENDDO
      RETURN
  END SUBROUTINE wmsetg_tor

!     ****** TASK/EQ wequilibrium ******

  SUBROUTINE wmsetg_eq(IERR)

    USE wmcomm
    IMPLICIT NONE
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: NRMAXSV=0,NTHMAXSV=0,NSUMAXSV=0

    IERR=0

    IF(NRMAXSV.EQ.NRMAX.AND. &
       NTHMAXSV.EQ.NTHMAX_F.AND. &
       NSUMAXSV.EQ.NSUMAX.AND. &
       KNAMEQ_SAVE.EQ.KNAMEQ) RETURN

      NSUMAX=41
      NSWMAX=NSUMAX
      NHHMAX=1
      NHHMAX_F=1

      IF(NTHMAX_F.LT.4) THEN
         IF(NRANK.EQ.0) &
              WRITE(6,*) 'XX WMXRZF: NTHMAX MUST BE GREATER THAN 4'
         IERR=1
         RETURN
      ENDIF

      IF(NRANK.EQ.0) THEN
         CALL EQLOAD(MODELG,KNAMEQ,IERR)
      ENDIF
      CALL mtx_broadcast1_integer(IERR)
      IF(IERR.EQ.1) RETURN

      NRMAXSV=NRMAX
      NTHMAXSV=NTHMAX_F
      NSUMAXSV=NSUMAX
      KNAMEQ_SAVE=KNAMEQ

      IF(NRANK.EQ.0) THEN
         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHMAX_F
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
         CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
         CALL EQGETP(RHOT,PSIP,NRMAX+1)
         CALL EQGETR(RPS,DRPSI,DRCHI,NTHMAX_F,NTHMAX_F,NRMAX+1)
         CALL EQGETZ(ZPS,DZPSI,DZCHI,NTHMAX_F,NTHMAX_F,NRMAX+1)
         CALL EQGETBB(BPR,BPZ,BPT,BTP,NTHMAX_F,NTHMAX_F,NRMAX+1)
         CALL EQGETQ(PPS,QPS,RBPS,VPS,RLEN,NRMAX+1)
         CALL EQGETU(RSU,ZSU,RSW,ZSW,NSUMAX)
         CALL EQGETF(RGMIN,RGMAX,ZGMIN,ZGMAX)
         CALL EQGETA(RAXIS,ZAXIS,PSIPA,PSITA,Q0,QA)

         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHGM
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
         CALL EQGETG(RPSG,ZPSG,NTHGM,NTHGM,NRMAX+1)
      ENDIF
      CALL mtx_broadcast1_real8(BB)
      CALL mtx_broadcast1_real8(RR)
      CALL mtx_broadcast1_real8(RIP)
      CALL mtx_broadcast1_real8(RA)
      CALL mtx_broadcast1_real8(RKAP)
      CALL mtx_broadcast1_real8(RDLT)
      CALL mtx_broadcast1_real8(RB)
      CALL mtx_broadcast_real8(RHOT,NRMAX+1)
      CALL mtx_broadcast_real8(PSIP,NRMAX+1)

      CALL mtx_broadcast_real8(RPS,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(DRPSI,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(DRCHI,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(ZPS,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(DZPSI,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(DZCHI,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(BPR,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(BPZ,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(BPT,NTHMF*(NRMAX+1))
      CALL mtx_broadcast_real8(BTP,NTHMF*(NRMAX+1))

      CALL mtx_broadcast_real8(PPS,NRMAX+1)
      CALL mtx_broadcast_real8(QPS,NRMAX+1)
      CALL mtx_broadcast_real8(RBPS,NRMAX+1)
      CALL mtx_broadcast_real8(VPS,NRMAX+1)
      CALL mtx_broadcast_real8(RLEN,NRMAX+1)
      CALL mtx_broadcast_real8(RSU,NSUMAX)
      CALL mtx_broadcast_real8(ZSU,NSUMAX)
      CALL mtx_broadcast_real8(RSW,NSUMAX)
      CALL mtx_broadcast_real8(ZSW,NSUMAX)
      CALL mtx_broadcast1_real8(RGMIN)
      CALL mtx_broadcast1_real8(RGMAX)
      CALL mtx_broadcast1_real8(ZGMIN)
      CALL mtx_broadcast1_real8(ZGMAX)
      CALL mtx_broadcast1_real8(RAXIS)
      CALL mtx_broadcast1_real8(ZAXIS)
      CALL mtx_broadcast1_real8(PSIPA)
      CALL mtx_broadcast1_real8(PSITA)
      CALL mtx_broadcast1_real8(Q0)
      CALL mtx_broadcast1_real8(QA)

         RMAX=RSU(1,1)
         DO NSU=2,NSUMAX
            RMAX=MAX(RMAX,RSU(NSU,1))
         ENDDO

         DO NR=1,NRMAX+1
            ARG=RHOT(NR)
            IF(ARG.LE.0.D0) THEN
               XRHO(NR)=0.D0
            ELSE
               XRHO(NR)=ARG
            ENDIF
            XR(NR)=RA*XRHO(NR)
         ENDDO

         DO NR=1,NRMAX+1
         DO NHH=1,NHHMAX_F
         DO NTH=1,NTHMAX_F
            RPST(NTH,NHH,NR)=RPS(NTH,NR)
            ZPST(NTH,NHH,NR)=ZPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO

         DO NR=2,NRMAX+1
         DO NTH=1,NTHMAX_F
         DO NHH=1,NHHMAX_F
            RG11(NTH,NHH,NR)= (DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2)
            RG12(NTH,NHH,NR)= (DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                        +DZPSI(NTH,NR)*DZCHI(NTH,NR))/XRHO(NR)
            RG13(NTH,NHH,NR)=0.D0
            RG22(NTH,NHH,NR)= (DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2)
     &                        /(XRHO(NR)*XRHO(NR))
            RG23(NTH,NHH,NR)=0.D0
            RG33(NTH,NHH,NR)= RPS(NTH,NR)**2
            RJ  (NTH,NHH,NR)= RPS(NTH,NR)
     &                      *( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                        -DRCHI(NTH,NR)*DZPSI(NTH,NR))/XRHO(NR)
            BPTL=(BPR(NTH,NR)*DZPSI(NTH,NR)
     &           -BPZ(NTH,NR)*DRPSI(NTH,NR))/XRHO(NR)
     &           /SQRT(RG11(NTH,NHH,NR))
     &           /SQRT(RG22(NTH,NHH,NR))

            BFLD(2,NTH,NHH,NR)=BPTL
            BFLD(3,NTH,NHH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO

         NR=1
         DO NHH=1,NHHMAX_F
         DO NTH=1,NTHMAX_F
            RG11(NTH,NHH,NR)= RG11(NTH,NHH,2)
            RG12(NTH,NHH,NR)= RG12(NTH,NHH,2)
            RG13(NTH,NHH,NR)= 0.D0
            RG22(NTH,NHH,NR)= RG22(NTH,NHH,2)
            RG23(NTH,NHH,NR)= 0.D0
            RG33(NTH,NHH,NR)= RPST(NTH,NHH,NR)**2
            RJ  (NTH,NHH,NR)= RJ(NTH,NHH,2)
            BPTL=0.D0
            BFLD(2,NTH,NHH,NR)=BFLD(2,NTH,NHH,2)
            BFLD(3,NTH,NHH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
         ENDDO
         ENDDO

!   *** recalculate for IPS ***

      IF(NRANK.EQ.0) THEN
         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHMAX_IPS_F
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
         CALL EQGETR(RPS_IPS,DRPSI_IPS,DRCHI_IPS,
     &               NTHMIPSF,NTHMAX_IPS_F,NRMAX+1)
         CALL EQGETZ(ZPS_IPS,DZPSI_IPS,DZCHI_IPS,
     &               NTHMIPSF,NTHMAX_IPS_F,NRMAX+1)
         CALL EQGETBB(BPR_IPS,BPZ_IPS,BPT_IPS,BTP_IPS,
     &                NTHMIPSF,NTHMAX_IPS_F,NRMAX+1)
      ENDIF

      CALL mtx_broadcast_real8(RPS_ips,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(DRPSI_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(DRCHI_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(ZPS_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(DZPSI_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(DZCHI_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(BPR_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(BPZ_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(BPT_IPS,NTHMIPSF*(NRMAX+1))
      CALL mtx_broadcast_real8(BTP_IPS,NTHMIPSF*(NRMAX+1))

         DO NR=1,NRMAX+1
         DO NHH=1,NHHMAX_IPS_F
         DO NTH=1,NTHMAX_IPS_F
            RPST_IPS(NTH,NHH,NR)=RPS_IPS(NTH,NR)
            ZPST_IPS(NTH,NHH,NR)=ZPS_IPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO

         DO NR=2,NRMAX+1
         DO NTH=1,NTHMAX_IPS_F
         DO NHH=1,NHHMAX_IPS_F
            RG11_IPS(NTH,NHH,NR)
     &           = (DRPSI_IPS(NTH,NR)**2+DZPSI_IPS(NTH,NR)**2)
            RG12_IPS(NTH,NHH,NR)
     &           = (DRPSI_IPS(NTH,NR)*DRCHI_IPS(NTH,NR)
     &             +DZPSI_IPS(NTH,NR)*DZCHI_IPS(NTH,NR))/XRHO(NR)
            RG13_IPS(NTH,NHH,NR)=0.D0
            RG22_IPS(NTH,NHH,NR)
     &           = (DRCHI_IPS(NTH,NR)**2+DZCHI_IPS(NTH,NR)**2)
     &                        /(XRHO(NR)*XRHO(NR))
            RG23_IPS(NTH,NHH,NR)=0.D0
            RG33_IPS(NTH,NHH,NR)= RPS_IPS(NTH,NR)**2
!            RJ_IPS(NTH,NHH,NR)  = RPS_IPS(NTH,NR)
!     &                      *( DRPSI_IPS(NTH,NR)*DZCHI_IPS(NTH,NR)
!     &                        -DRCHI_IPS(NTH,NR)*DZPSI_IPS(NTH,NR))/XRHO(NR)
            BPTL_IPS=(BPR_IPS(NTH,NR)*DZPSI_IPS(NTH,NR)
     &               -BPZ_IPS(NTH,NR)*DRPSI_IPS(NTH,NR))/XRHO(NR)
     &               /SQRT(RG11_IPS(NTH,NHH,NR))
     &               /SQRT(RG22_IPS(NTH,NHH,NR))

            BFLD_IPS(2,NTH,NHH,NR)=BPTL_IPS
            BFLD_IPS(3,NTH,NHH,NR)=BTP_IPS(NTH,NR)/RPS_IPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO

         NR=1
         DO NHH=1,NHHMAX_IPS_F
         DO NTH=1,NTHMAX_IPS_F
            RG11_IPS(NTH,NHH,NR)= RG11_IPS(NTH,NHH,2)
            RG12_IPS(NTH,NHH,NR)= RG12_IPS(NTH,NHH,2)
            RG13_IPS(NTH,NHH,NR)= 0.D0
            RG22_IPS(NTH,NHH,NR)= RG22_IPS(NTH,NHH,2)
            RG23_IPS(NTH,NHH,NR)= 0.D0
            RG33_IPS(NTH,NHH,NR)= RPST_IPS(NTH,NHH,NR)**2
!            RJ_IPS(NTH,NHH,NR)= RJ_IPS(NTH,NHH,2)
            BFLD_IPS(2,NTH,NHH,NR)=BFLD_IPS(2,NTH,NHH,2)
            BFLD_IPS(3,NTH,NHH,NR)=BTP_IPS(NTH,NR)/RPS_IPS(NTH,NR)
         ENDDO
         ENDDO
      RETURN
  END SUBROUTINE wmsetg_eq

!     ****** CALCULATE ION CYCROTRON RESONANCE SURFACE ******

  SUBROUTINE WMICRS

      INCLUDE 'wmcomm.inc'

      DO NR=1,NRMAX+1
!      DO NHH=1,NHHMAX
!      DO NTH=1,NTHMAX
!  seki
      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
!  seki
         BSUPTH=BFLD(2,NTH,NHH,NR)
         BSUPPH=BFLD(3,NTH,NHH,NR)
         BABS=SQRT(     RG22(NTH,NHH,NR)*BSUPTH*BSUPTH*XRHO(NR)**2
     &            +2.D0*RG23(NTH,NHH,NR)*BSUPTH*BSUPPH*XRHO(NR)
     &            +     RG33(NTH,NHH,NR)*BSUPPH*BSUPPH)
!         BABS=SQRT(BPT(NTH,NR)**2+BTP(NTH,NR)**2)
         BPST(NTH,NHH,NR)=BABS
!         IF(NTH.EQ.1) THEN
!            WRITE(6,'(I5,1P4E12.4)') NR,BABS,BABS1,BSUPTH,BSUPPH
!         ENDIF
      ENDDO
      ENDDO
      DO NHH=1,NHHMAX_IPS_F
      DO NTH=1,NTHMAX_IPS_F
         BSUPTH=BFLD_IPS(2,NTH,NHH,NR)
         BSUPPH=BFLD_IPS(3,NTH,NHH,NR)
         BABS=SQRT(     RG22_IPS(NTH,NHH,NR)*BSUPTH*BSUPTH*XRHO(NR)**2
     &            +2.D0*RG23_IPS(NTH,NHH,NR)*BSUPTH*BSUPPH*XRHO(NR)
     &            +     RG33_IPS(NTH,NHH,NR)*BSUPPH*BSUPPH)
         BPST_IPS(NTH,NHH,NR)=BABS

      ENDDO
      ENDDO
      ENDDO

      RETURN
  END SUBROUTINE WMICRS

!     ****** CALCULATE INVERT MATRIX ******

  SUBROUTINE WMRGINVERT

      INCLUDE 'wmcomm.inc'
      DIMENSION RGA(3,3),RGB(3,3)
!        ----- Invert matrix to obtain mu^(-1)=RMB and g^(-1)=RGB ----
      DO NR  =1,NRMAX+1
      DO NHHF=1,NHHMAX_F
      DO NTHF=1,NTHMAX_F
         RGA(1,1)=RG11(NTHF,NHHF,NR)
         RGA(1,2)=RG12(NTHF,NHHF,NR)
         RGA(1,3)=RG13(NTHF,NHHF,NR)
         RGA(2,1)=RG12(NTHF,NHHF,NR)
         RGA(2,2)=RG22(NTHF,NHHF,NR)
         RGA(2,3)=RG23(NTHF,NHHF,NR)
         RGA(3,1)=RG13(NTHF,NHHF,NR)
         RGA(3,2)=RG23(NTHF,NHHF,NR)
         RGA(3,3)=RG33(NTHF,NHHF,NR)
         DO J=1,3
         DO I=1,3
            RGB(I,J)=RGA(I,J)
         ENDDO
         ENDDO
         CALL INVMRD(RGB,3,3,ILL)
         RGI11(NTHF,NHHF,NR)=RGB(1,1)
         RGI12(NTHF,NHHF,NR)=RGB(1,2)
         RGI13(NTHF,NHHF,NR)=RGB(1,3)
         RGI22(NTHF,NHHF,NR)=RGB(2,2)
         RGI23(NTHF,NHHF,NR)=RGB(2,3)
         RGI33(NTHF,NHHF,NR)=RGB(3,3)
      ENDDO
      ENDDO
      ENDDO

      DO NR  =1,NRMAX+1
      DO NHHF=1,NHHMAX_F
      DO NTHF=1,NTHMAX_F
         RGA(1,1)=RGH11(NTHF,NHHF,NR)
         RGA(1,2)=RGH12(NTHF,NHHF,NR)
         RGA(1,3)=RGH13(NTHF,NHHF,NR)
         RGA(2,1)=RGH12(NTHF,NHHF,NR)
         RGA(2,2)=RGH22(NTHF,NHHF,NR)
         RGA(2,3)=RGH23(NTHF,NHHF,NR)
         RGA(3,1)=RGH13(NTHF,NHHF,NR)
         RGA(3,2)=RGH23(NTHF,NHHF,NR)
         RGA(3,3)=RGH33(NTHF,NHHF,NR)
         DO J=1,3
         DO I=1,3
            RGB(I,J)=RGA(I,J)
         ENDDO
         ENDDO
         CALL INVMRD(RGB,3,3,ILL)
         RGIH11(NTHF,NHHF,NR)=RGB(1,1)
         RGIH12(NTHF,NHHF,NR)=RGB(1,2)
         RGIH13(NTHF,NHHF,NR)=RGB(1,3)
         RGIH22(NTHF,NHHF,NR)=RGB(2,2)
         RGIH23(NTHF,NHHF,NR)=RGB(2,3)
         RGIH33(NTHF,NHHF,NR)=RGB(3,3)
      ENDDO
      ENDDO
      ENDDO

      DO NR  =1,NRMAX+1
      DO NHHF=1,NHHMAX_IPS_F
      DO NTHF=1,NTHMAX_IPS_F
         RGA(1,1)=RG11_IPS(NTHF,NHHF,NR)
         RGA(1,2)=RG12_IPS(NTHF,NHHF,NR)
         RGA(1,3)=RG13_IPS(NTHF,NHHF,NR)
         RGA(2,1)=RG12_IPS(NTHF,NHHF,NR)
         RGA(2,2)=RG22_IPS(NTHF,NHHF,NR)
         RGA(2,3)=RG23_IPS(NTHF,NHHF,NR)
         RGA(3,1)=RG13_IPS(NTHF,NHHF,NR)
         RGA(3,2)=RG23_IPS(NTHF,NHHF,NR)
         RGA(3,3)=RG33_IPS(NTHF,NHHF,NR)
         DO J=1,3
         DO I=1,3
            RGB(I,J)=RGA(I,J)
         ENDDO
         ENDDO
         CALL INVMRD(RGB,3,3,ILL)
         RGI11_IPS(NTHF,NHHF,NR)=RGB(1,1)
         RGI12_IPS(NTHF,NHHF,NR)=RGB(1,2)
         RGI13_IPS(NTHF,NHHF,NR)=RGB(1,3)
         RGI22_IPS(NTHF,NHHF,NR)=RGB(2,2)
         RGI23_IPS(NTHF,NHHF,NR)=RGB(2,3)
         RGI33_IPS(NTHF,NHHF,NR)=RGB(3,3)
      ENDDO
      ENDDO
      ENDDO
  END SUBROUTINE WMRGINVERT

!

  SUBROUTINE BTHOB
      INCLUDE 'wmcomm.inc'
      DIMENSION RGA(3,3),RGB(3,3)
!        ----- Invert matrix to obtain mu^(-1)=RMB and g^(-1)=RGB ----
      DO NR  =1,NRMAX+1
         BTHOBN(NR)=0d0
         RHON=XRHO(NR)
      DO NHHF=1,NHHMAX_F
      DO NTHF=1,NTHMAX_F
         BSUPTH=BFLD(2,NTHF,NHHF,NR)
         BSUPPH=BFLD(3,NTHF,NHHF,NR)
         BABS  =BPST(NTHF,NHHF,NR)
         BTH=ABS(BSUPTH*RA*RHON)
         BTHOBN(NR)=BTHOBN(NR) + BTH/BABS
      ENDDO
      ENDDO
         BTHOBN(NR)=ABS(BTHOBN(NR))/DBLE(NHHMAX_F*NTHMAX_F)
         print *,"BTHOBN",BTHOBN(NR)
      ENDDO
  END SUBROUTINE BTHOB
END MODULE wmsetg
