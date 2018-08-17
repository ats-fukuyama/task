C     $Id: wmsetg.f,v 1.43 2014/10/11 05:33:14 fukuyama Exp $
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE WMSETG(IERR)
C
      INCLUDE 'wmcomm.inc'
C
C
      IERR=0
C
      IF(NHHMAX.EQ.1) THEN
         NHHMAX_F = 1
         NHHMAX_IPS_F = 1
      ELSE
         NHHMAX_F  = NHHMAX*NFACT
         NHHMAX_IPS_F  = NDMAX_IPS
      ENDIF
C
      IF(NTHMAX.EQ.1) THEN
         NTHMAX_F  = 1
         NTHMAX_IPS_F  = 1
      ELSE
         NTHMAX_F  = NTHMAX*MFACT
         NTHMAX_IPS_F  = MDMAX_IPS
      ENDIF
C
C     ****** Cylindrical COORDINATES ******
C
      IF(MODELG.EQ.0) THEN
         CALL WMXRZC(ierr)
C
C     ****** CYLINDRICAL COORDINATES ******
C
      ELSEIF(MODELG.EQ.1) THEN
         CALL WMXRZT(IERR)
C
C     ****** TOROIDAL COORDINATES ******
C
      ELSEIF(MODELG.EQ.2) THEN
         CALL WMXRZT(IERR)
C
C     ****** EQUILIBRIUM (TASK/EQ) ******
C
      ELSEIF(MODELG.EQ.3) THEN
         CALL WMXRZF(IERR)
C         IF(MODEL_UHR.EQ.0) THEN
C            CALL WMXRZF(IERR)
C         ELSEIF(MODEL_UHR.EQ.1) THEN
C            CALL WMXRZF_UHR(IERR)
C         ENDIF
C
C     ****** EQUILIBRIUM (VMEC) ******
C
C     atode
      ELSEIF(MODELG.EQ.4) THEN
         CALL WMXRZV(IERR)
C
C     ****** EQUILIBRIUM (EQDSK) ******
C
      ELSEIF(MODELG.EQ.5) THEN
         CALL WMXRZF(IERR)
C
C     ****** EQUILIBRIUM (BOOZER) ******
C
      ELSEIF(MODELG.EQ.6) THEN
         CALL WMXRZB(IERR)
      ENDIF
C
      
      IF(IERR.NE.0) RETURN
C
      IF(MODELN.EQ.7) CALL WMDPRF(IERR)
      IF(MODELN.EQ.8) CALL WMXPRF(IERR)
      IF(MODELN.EQ.9) THEN
         IF(KNAMTR.NE.KNAMTR_SAVE) THEN
C            CALL WMTRLOAD(KNAMTR,IERR)
            KNAMTR_SAVE=KNAMTR
         ENDIF
      ENDIF
      IF(IERR.NE.0) RETURN
C
      IF(NHHMAX.EQ.1) THEN
         NHHMAX_F = 1
         NHHMAX_IPS_F = 1
         NDSIZ  = 1
         NDMSIZ  = 1   !
         NDSIZ_F  = 1
         NDMIN  = 0
         NDMMIN  = 0   !
         NDMIN_F  = 0
         NDMAX  = 0
         NDMMAX  = 0   !
         NDMAX_F  = 0
         KDSIZ  = 1
         KDMSIZ  = 1   !
         KDSIZ_F  = 1
         KDSIZ_IPS_F  = 1
         KDMIN  = 0
         KDMMIN  = 0   !
         KDMIN_F  = 0
         KDMIN_IPS_F  = 0
         KDMAX  = 0
         KDMMAX  = 0   !
         KDMAX_F  = 0
         KDMAX_IPS_F  = 0
         NDSIZX = 1
         NDSIZX_F = 1
      ELSE
         NDSIZ  = NHHMAX
         NDMSIZ  = NHHMAX_SV    !
         NHHMAX_F  = NHHMAX*NFACT
!         NHHMAX_IPS_F  = NHHMAX_F*NIPSFACT
         NHHMAX_IPS_F  = NDMIPSF
         NDSIZ_F  = NHHMAX_F
         NDMIN  =-NHHMAX/2+1
         NDMMIN  =-NHHMAX_SV/2+1  !
         NDMIN_F  =-NHHMAX_F/2+1
         NDMAX  = NHHMAX/2
         NDMMAX  = NHHMAX_SV/2    !
         NDMAX_F  = NHHMAX_F/2
         KDSIZ  = NHHMAX
         KDMSIZ  = NHHMAX_SV      !
         KDSIZ_F  = NHHMAX_F
         KDSIZ_IPS_F  = NHHMAX_IPS_F
         KDMIN  =-NHHMAX/2+1
         KDMMIN  =-NHHMAX_SV/2+1  !
         KDMIN_F  =-NHHMAX_F/2+1
         KDMIN_IPS_F  =-NHHMAX_IPS_F/2+1
         KDMAX  = NHHMAX/2
         KDMMAX  = NHHMAX_SV/2    !
         KDMAX_F  = NHHMAX_F/2
         KDMAX_IPS_F  = NHHMAX_IPS_F/2
C         NDSIZX = 3*NHHMAX/2
C         NDSIZX_F = 3*NHHMAX_F/2
         NDSIZX = (NHHMAX + NHHMAX_F)/2
!         NDSIZX = NHHMAX + (NHHMAX_F)/2
         NDSIZX_F = NDSIZX
C         NDSIZX_F = 3*NHHMAX_F
      ENDIF
C
      IF(NTHMAX.EQ.1) THEN
         MDSIZ  = 1
         MDMSIZ  = 1    !
         MDSIZ_F  = 1
         MDMIN  = 0
         MDMMIN  = 0    !
         MDMIN_F  = 0
         MDMAX  = 0
         MDMMAX  = 0    !
         MDMAX_F  = 0
         LDSIZ  = 1
         LDMSIZ  = 1    !
         LDSIZ_F  = 1
         LDSIZ_IPS_F  = 1
         LDMIN  = 0
         LDMMIN  = 0    !
         LDMIN_F  = 0
         LDMIN_IPS_F  = 0
         LDMAX  = 0
         LDMMAX  = 0    !
         LDMAX_F  = 0
         LDMAX_IPS_F  = 0
         MDSIZX = 1
         MDSIZX_F = 1
      ELSE
         MDSIZ  = NTHMAX
         NTHMAX_F  = NTHMAX*MFACT
!         NTHMAX_IPS_F  = NTHMAX_F*MIPSFACT
         NTHMAX_IPS_F  = MDMIPSF
         MDMSIZ  = NTHMAX        !
         MDSIZ_F  = NTHMAX_F
         MDMIN  =-NTHMAX/2+1
         MDMMIN  =-NTHMAX/2+1    !
         MDMIN_F  =-NTHMAX_F/2+1
         MDMAX  = NTHMAX/2
         MDMMAX  = NTHMAX/2      !
         MDMAX_F  = NTHMAX_F/2
         LDSIZ  = NTHMAX
         LDMSIZ  = NTHMAX        !
         LDSIZ_F  = NTHMAX_F
         LDSIZ_IPS_F  = NTHMAX_IPS_F
         LDMIN  =-NTHMAX/2+1
         LDMMIN  =-NTHMAX/2+1    !
         LDMIN_F  =-NTHMAX_F/2+1
         LDMIN_IPS_F  =-NTHMAX_IPS_F/2+1
         LDMAX  = NTHMAX/2
         LDMMAX  = NTHMAX/2      !
         LDMAX_F  = NTHMAX_F/2
         LDMAX_IPS_F  = NTHMAX_IPS_F/2
C         LDMAX  = NTHMAX/2-1
C         LDMAX_F  = NTHMAX_F/2-1
C         MDSIZX = 3*NTHMAX/2
C         MDSIZX_F = 3*NTHMAX_F/2
         MDSIZX = (NTHMAX+NTHMAX_F)/2
!         MDSIZX = NTHMAX+ (NTHMAX_F)/2
         MDSIZX_F = MDSIZX
C         MDSIZX_F = 3*NTHMAX_F
      ENDIF

C
      IF(MODELG.NE.6) CALL WMICRS
      IERR=0
      RETURN
      END
C
C     ****** RADIAL MESH (CYLINDRICAL COORDINATES) ******
C
      SUBROUTINE WMXRZC(IERR)
C
      INCLUDE 'wmcomm.inc'
C
         IERR=0
C
         NSUMAX=31
         NSWMAX=31
         NHHMAX=1
         NHHMAX_F=1
C
         PSIPA=RA*RA*BB/(Q0+QA)
         PSIPB=SQRT(RB**2/RA**2+(RB**2/RA**2-1.D0)*Q0/QA)*PSIPA
         DRHO=SQRT(PSIPB/PSIPA)/NRMAX
C
         DO NR=1,NRMAX+1
            RHOL=DRHO*(NR-1)
            XRHO(NR)=RHOL
C
            IF(NR.EQ.1) THEN
               XR(NR)=0.D0
            ELSEIF(RHOL.LT.1.D0) THEN
               XR(NR)=RA*SQRT((2.D0*Q0*RHOL**2+(QA-Q0)*RHOL**4)
     &                        /(Q0+QA))
            ELSE
               XR(NR)=RA*SQRT((Q0+QA*RHOL**4)/(Q0+QA))
            ENDIF
C
C            WRITE(6,*) 'NR,RHO,XR=',NR,XRHO(NR),XR(NR)
C
            IF(RHOL.LT.1.D0) THEN
               QPS(NR)=Q0+(QA-Q0)*RHOL**2
            ELSE
               QPS(NR)=QA*RHOL**2
            ENDIF
         ENDDO
C
!         DTH=2.D0*PI/NTHMAX
!  seki
         DTH=2.D0*PI/NTHMAX_F
!  seki
         DTHG=2.D0*PI/NTHGM
         DO NR=1,NRMAX+1
            IF(NR.EQ.1) THEN
               RS=XR(2)/XRHO(2)
            ELSE
               RS=XR(NR)/XRHO(NR)
            ENDIF
            RSD=QPS(NR)/(BB*RS)
C            WRITE(6,*) 'NR,RS,RSD=',NR,RS,RSD
!  seki
            DO NTH=1,NTHMAX_F
!  seki
C            DO NTH=1,NTHMAX
               RCOS=COS(DTH*(NTH-1))
               RSIN=SIN(DTH*(NTH-1))
               RPS(NTH,NR)    = RR + XR(NR)*RCOS
               ZPS(NTH,NR)    =      XR(NR)*RSIN
               DRPSI(NTH,NR)  =      RSD   *RCOS
               DZPSI(NTH,NR)  =      RSD   *RSIN
               DRCHI(NTH,NR)  =     -RS    *RSIN
               DZCHI(NTH,NR)  =      RS    *RCOS
            ENDDO
            DO NTH=1,NTHGM
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
C
         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6
C
         DO NR=1,NRMAX+1
            RHOL=XRHO(NR)
            IF(RHOL.LE.1.D0) THEN
               IF(PN(1).LE.0.D0) THEN
                  FEDGE=0.D0
               ELSE
                  FEDGE=PNS(1)/PN(1)
               END IF
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)
     &              +FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)
     &              +FEDGE
               PPS(NR)=P0*FACTN*FACTT
               PPS(NR)=0.D0
            ENDIF
            RBPS(NR)=BB*RR
            VPS(NR)=2*PI*RR*PI*XR(NR)**2
            RLEN(NR)=2*PI*XR(NR)
         ENDDO
C
         DTHU=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            RCOS=COS(DTHU*(NSU-1))
            RSIN=SIN(DTHU*(NSU-1))
            RSU(NSU,1)=RR+RA*RCOS
            RSW(NSU,1)=RR+RB*RCOS
            ZSU(NSU,1)=RA*RSIN
            ZSW(NSU,1)=RB*RSIN
         ENDDO
C
         RGMIN= RR-RB*1.01D0
         RGMAX= RR+RB*1.01D0
         ZGMIN=-RB*1.01D0
         ZGMAX= RB*1.01D0

C
         DO NR=1,NRMAX+1
C         DO NTH=1,NTHMAX
!  seki
         DO NTH=1,NTHMAX_F
!  seki
            RRG=RPS(NTH,NR)
C         DO NHH=1,NHHMAX
!  seki
         DO NHH=1,NHHMAX_F
!  seki
            RPST(NTH,NHH,NR)=RPS(NTH,NR)
            ZPST(NTH,NHH,NR)=ZPS(NTH,NR)
C
            RG11(NTH,NHH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
            RG12(NTH,NHH,NR)= 0d0
C            RG12(NTH,NHH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR)
C     &                       +DZPSI(NTH,NR)*DZCHI(NTH,NR)
            RG13(NTH,NHH,NR)= 0.D0
            RG22(NTH,NHH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
            RG23(NTH,NHH,NR)= 0.D0
            RG33(NTH,NHH,NR)= RRG**2
            RJ  (NTH,NHH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                             -DRCHI(NTH,NR)*DZPSI(NTH,NR))
C
!            BFLD(2,NTH,NHH,NR)=1.D0/RJ(NTH,NHH,NR)
            BFLD(2,NTH,NHH,NR)=0d0
            BFLD(3,NTH,NHH,NR)=RBPS(NR)/RRG**2
C            WRITE(6,*) 'NR,RJ,BFLD2,BFLD3=',NR,RJ(NTH,NHH,NR),
C     &                 BFLD(2,NTH,NHH,NR),BFLD(3,NTH,NHH,NR)
         ENDDO
         ENDDO
         ENDDO
      RETURN
      END
C
C     ****** RADIAL MESH (TOROIDAL COORDINATES) ******
C
      SUBROUTINE WMXRZT(IERR)
C
C      USE plfile_prof_mod
      USE plprof,ONLY: pl_qprf
      INCLUDE 'wmcomm.inc'
C
C
         IERR=0
C
         NSUMAX=31
         NSWMAX=31
         NHHMAX=1
         NHHMAX_F=1
         NHHMAX_IPS_F=1

C         CALL plfile_prof_read(modeln,modelq,ierr)
C
         PSIPA=RA*RA*BB/(Q0+QA)
         DRHO=(RB/RA)/NRMAX
C
         DO NR=1,NRMAX+1
            RHOL=DRHO*(NR-1)
            XRHO(NR)=RHOL
            XR(NR)=RA*RHOL
            CALL PL_QPRF(RHOL,QPS(NR))
         ENDDO
C
!  seki
         DTH=2.D0*PI/NTHMAX_F
!  seki
         DTHG=2.D0*PI/NTHGM
         DO NR=1,NRMAX+1
            RSD=RA/(2.D0*PSIPA)
C            DO NTH=1,NTHMAX
!  seki
            DO NTH=1,NTHMAX_F
!  seki
               RCOS=COS(DTH*(NTH-1))
               RSIN=SIN(DTH*(NTH-1))
               RPS(NTH,NR) = RR + XR(NR)*RCOS
               ZPS(NTH,NR)    =      XR(NR)*RSIN
               DRPSI(NTH,NR)  =      RSD   *RCOS
               DZPSI(NTH,NR)  =      RSD   *RSIN
               DRCHI(NTH,NR)  =     -RA    *RSIN
               DZCHI(NTH,NR)  =      RA    *RCOS
            ENDDO
            DO NTH=1,NTHGM
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
C
         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6
C
         DO NR=1,NRMAX+1
            RHOL=XRHO(NR)
            IF(RHOL.LE.1.D0) THEN
               IF(PN(1).LE.0.D0) THEN
                  FEDGE=0.D0
               ELSE
                  FEDGE=PNS(1)/PN(1)
               END IF
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)
     &              +FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)
     &              +FEDGE
               PPS(NR)=P0*FACTN*FACTT
            ELSE
               PPS(NR)=0.D0
            ENDIF
            RBPS(NR)=BB*RR
            VPS(NR)=2*PI*RR*PI*XR(NR)**2
            RLEN(NR)=2*PI*XR(NR)
         ENDDO
C
         DTHU=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            RCOS=COS(DTHU*(NSU-1))
            RSIN=SIN(DTHU*(NSU-1))
            RSU(NSU,1)=RR+RA*RCOS
            RSW(NSU,1)=RR+RB*RCOS
            ZSU(NSU,1)=RA*RSIN
            ZSW(NSU,1)=RB*RSIN
         ENDDO
C
         RGMIN=RR-RB*1.01D0
         RGMAX=RR+RB*1.01D0
         ZGMIN=-RB*1.01D0
         ZGMAX= RB*1.01D0
         RAXIS=RR

C
         DO NR=1,NRMAX+1
!  seki
C         DO NTH=1,NTHMAX
         DO NTH=1,NTHMAX_F
!  seki
            RRG=RPS(NTH,NR)
!  seki
C         DO NHH=1,NHHMAX
         DO NHH=1,NHHMAX_F
!  seki
            RPST(NTH,NHH,NR)=RPS(NTH,NR)
            ZPST(NTH,NHH,NR)=ZPS(NTH,NR)
C
            RG11(NTH,NHH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
            RG12(NTH,NHH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                       +DZPSI(NTH,NR)*DZCHI(NTH,NR)
            RG13(NTH,NHH,NR)= 0.D0
            RG22(NTH,NHH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
C            if(nth.eq.1) write(6,'(A,I5,1P3E12.4)') 
C     &           '-- ',NR,DRCHI(NTH,NR),DZCHI(NTH,NR),RG22(NTH,NHH,NR)
            RG23(NTH,NHH,NR)= 0.D0
            RG33(NTH,NHH,NR)= RRG**2
            RJ  (NTH,NHH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                             -DRCHI(NTH,NR)*DZPSI(NTH,NR))
C
            BFLD(2,NTH,NHH,NR)=BB*RR/RRG**2/QPS(NR)
            BFLD(3,NTH,NHH,NR)=BB*RR/RRG**2
C            WRITE(6,*) 'NR,RJ,BFLD2,BFLD3=',NR,RJ(NTH,NHH,NR),
C     &                 BFLD(2,NTH,NHH,NR),BFLD(3,NTH,NHH,NR)
         ENDDO
         ENDDO
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  seki
         DTH=2.D0*PI/NTHMAX_IPS_F
         DO NR=1,NRMAX+1
            RSD=RA/(2.D0*PSIPA)
C            DO NTH=1,NTHMAX
            DO NTH=1,NTHMAX_IPS_F
               RCOS=COS(DTH*(NTH-1))
               RSIN=SIN(DTH*(NTH-1))
               RPS_IPS(NTH,NR) = RR + XR(NR)*RCOS
               ZPS_IPS(NTH,NR)    =      XR(NR)*RSIN
               DRPSI_IPS(NTH,NR)  =      RSD   *RCOS
               DZPSI_IPS(NTH,NR)  =      RSD   *RSIN
               DRCHI_IPS(NTH,NR)  =     -RA    *RSIN
               DZCHI_IPS(NTH,NR)  =      RA    *RCOS
            ENDDO
         ENDDO
C
         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6
C
         DO NR=1,NRMAX+1
            RHOL=XRHO(NR)
            IF(RHOL.LE.1.D0) THEN
               IF(PN(1).LE.0.D0) THEN
                  FEDGE=0.D0
               ELSE
                  FEDGE=PNS(1)/PN(1)
               END IF
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)
     &              +FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)
     &              +FEDGE
               PPS(NR)=P0*FACTN*FACTT
            ELSE
               PPS(NR)=0.D0
            ENDIF
            RBPS(NR)=BB*RR
            VPS(NR)=2*PI*RR*PI*XR(NR)**2
            RLEN(NR)=2*PI*XR(NR)
         ENDDO
C
C
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX_IPS_F
!  seki
            RRG=RPS_IPS(NTH,NR)
         DO NHH=1,NHHMAX_IPS_F
!  seki
            RPST_IPS(NTH,NHH,NR)=RPS_IPS(NTH,NR)
            ZPST_IPS(NTH,NHH,NR)=ZPS_IPS(NTH,NR)
C
            RG11_IPS(NTH,NHH,NR)= DRPSI_IPS(NTH,NR)**2
     &                       +DZPSI_IPS(NTH,NR)**2
            RG12_IPS(NTH,NHH,NR)= DRPSI_IPS(NTH,NR)*DRCHI_IPS(NTH,NR)
     &                       +DZPSI_IPS(NTH,NR)*DZCHI_IPS(NTH,NR)
            RG13_IPS(NTH,NHH,NR)= 0.D0
            RG22_IPS(NTH,NHH,NR)= DRCHI_IPS(NTH,NR)**2
     &                       +DZCHI_IPS(NTH,NR)**2
C            if(nth.eq.1) write(6,'(A,I5,1P3E12.4)') 
C     &           '-- ',NR,DRCHI(NTH,NR),DZCHI(NTH,NR),RG22(NTH,NHH,NR)
            RG23_IPS(NTH,NHH,NR)= 0.D0
            RG33_IPS(NTH,NHH,NR)= RRG**2
C
            BFLD_IPS(2,NTH,NHH,NR)=BB*RR/RRG**2/QPS(NR)
            BFLD_IPS(3,NTH,NHH,NR)=BB*RR/RRG**2
C            WRITE(6,*) 'NR,RJ,BFLD2,BFLD3=',NR,RJ(NTH,NHH,NR),
C     &                 BFLD(2,NTH,NHH,NR),BFLD(3,NTH,NHH,NR)
         ENDDO
         ENDDO
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RETURN
      END
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE WMXRZF(IERR)
C
C      USE plfile_prof_mod
      INCLUDE 'wmcomm.inc'
C
      CHARACTER*(80) LINE
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV
      DATA NRMAXSV,NTHMAXSV,NSUMAXSV/0,0,0/
C
      IERR=0
C
      IF(NRMAXSV.EQ.NRMAX.AND.
C     &   NTHMAXSV.EQ.NTHMAX.AND.
!  seki
     &   NTHMAXSV.EQ.NTHMAX_F.AND.
!  seki
     &   NSUMAXSV.EQ.NSUMAX.AND.
     &   KNAMEQ_SAVE.EQ.KNAMEQ) RETURN
C
      NSUMAX=41
      NSWMAX=NSUMAX
      NHHMAX=1
      NHHMAX_F=1
C
      IF(NTHMAX_F.LT.4) THEN
         IF(NRANK.EQ.0) 
     &        WRITE(6,*) 'XX WMXRZF: NTHMAX MUST BE GREATER THAN 4'
         IERR=1
         RETURN
      ENDIF
C
      IF(NRANK.EQ.0) THEN
         CALL EQLOAD(MODELG,KNAMEQ,IERR)
      ENDIF
      CALL mtx_broadcast1_integer(IERR)
      IF(IERR.EQ.1) RETURN
C
      NRMAXSV=NRMAX
!  seki
      NTHMAXSV=NTHMAX_F
!  seki
      NSUMAXSV=NSUMAX
      KNAMEQ_SAVE=KNAMEQ
C
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
         CALL EQGETR(RPS,DRPSI,DRCHI,NTHMF,NTHMAX_F,NRMAX+1)
         CALL EQGETZ(ZPS,DZPSI,DZCHI,NTHMF,NTHMAX_F,NRMAX+1)
         CALL EQGETBB(BPR,BPZ,BPT,BTP,NTHMF,NTHMAX_F,NRMAX+1)
         CALL EQGETQ(PPS,QPS,RBPS,VPS,RLEN,NRMAX+1)
         CALL EQGETU(RSU,ZSU,RSW,ZSW,NSUMAX)
         CALL EQGETF(RGMIN,RGMAX,ZGMIN,ZGMAX)
         CALL EQGETA(RAXIS,ZAXIS,PSIPA,PSITA,Q0,QA)
C
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
C
         RMAX=RSU(1,1)
         DO NSU=2,NSUMAX
            RMAX=MAX(RMAX,RSU(NSU,1))
         ENDDO
C
         DO NR=1,NRMAX+1
            ARG=RHOT(NR)
            IF(ARG.LE.0.D0) THEN
               XRHO(NR)=0.D0
            ELSE
               XRHO(NR)=ARG
            ENDIF
            XR(NR)=RA*XRHO(NR)
         ENDDO
C
         DO NR=1,NRMAX+1
         DO NHH=1,NHHMAX_F
         DO NTH=1,NTHMAX_F
            RPST(NTH,NHH,NR)=RPS(NTH,NR)
            ZPST(NTH,NHH,NR)=ZPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO
C
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
C
            BFLD(2,NTH,NHH,NR)=BPTL
            BFLD(3,NTH,NHH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO
C
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

C   *** recalculate for IPS ***

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
C
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
C            RJ_IPS(NTH,NHH,NR)  = RPS_IPS(NTH,NR)
C     &                      *( DRPSI_IPS(NTH,NR)*DZCHI_IPS(NTH,NR)
C     &                        -DRCHI_IPS(NTH,NR)*DZPSI_IPS(NTH,NR))/XRHO(NR)
            BPTL_IPS=(BPR_IPS(NTH,NR)*DZPSI_IPS(NTH,NR)
     &               -BPZ_IPS(NTH,NR)*DRPSI_IPS(NTH,NR))/XRHO(NR)
     &               /SQRT(RG11_IPS(NTH,NHH,NR))
     &               /SQRT(RG22_IPS(NTH,NHH,NR))
C
            BFLD_IPS(2,NTH,NHH,NR)=BPTL_IPS
            BFLD_IPS(3,NTH,NHH,NR)=BTP_IPS(NTH,NR)/RPS_IPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO
C
         NR=1
         DO NHH=1,NHHMAX_IPS_F
         DO NTH=1,NTHMAX_IPS_F
            RG11_IPS(NTH,NHH,NR)= RG11_IPS(NTH,NHH,2)
            RG12_IPS(NTH,NHH,NR)= RG12_IPS(NTH,NHH,2)
            RG13_IPS(NTH,NHH,NR)= 0.D0
            RG22_IPS(NTH,NHH,NR)= RG22_IPS(NTH,NHH,2)
            RG23_IPS(NTH,NHH,NR)= 0.D0
            RG33_IPS(NTH,NHH,NR)= RPST_IPS(NTH,NHH,NR)**2
C            RJ_IPS(NTH,NHH,NR)= RJ_IPS(NTH,NHH,2)
            BFLD_IPS(2,NTH,NHH,NR)=BFLD_IPS(2,NTH,NHH,2)
            BFLD_IPS(3,NTH,NHH,NR)=BTP_IPS(NTH,NR)/RPS_IPS(NTH,NR)
         ENDDO
         ENDDO
      RETURN
      END
C
C     ****** CALCULATE ION CYCROTRON RESONANCE SURFACE ******
C
      SUBROUTINE WMICRS
C
      INCLUDE 'wmcomm.inc'
C
      DO NR=1,NRMAX+1
C      DO NHH=1,NHHMAX
C      DO NTH=1,NTHMAX
!  seki
      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
!  seki
         BSUPTH=BFLD(2,NTH,NHH,NR)
         BSUPPH=BFLD(3,NTH,NHH,NR)
         BABS=SQRT(     RG22(NTH,NHH,NR)*BSUPTH*BSUPTH*XRHO(NR)**2
     &            +2.D0*RG23(NTH,NHH,NR)*BSUPTH*BSUPPH*XRHO(NR)
     &            +     RG33(NTH,NHH,NR)*BSUPPH*BSUPPH)
C         BABS=SQRT(BPT(NTH,NR)**2+BTP(NTH,NR)**2)
         BPST(NTH,NHH,NR)=BABS
C         IF(NTH.EQ.1) THEN
C            WRITE(6,'(I5,1P4E12.4)') NR,BABS,BABS1,BSUPTH,BSUPPH
C         ENDIF
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
C
      RETURN
      END
