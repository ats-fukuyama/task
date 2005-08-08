C     $Id$
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE WMSETG(IERR)
C
      INCLUDE 'wmcomm.inc'
C
C     ****** RECTANGULAR COORDINATES ******
C
      IF(MODELG.EQ.0) THEN
         CALL WMXRZC(IERR)
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
C
C     ****** EQUILIBRIUM (VMEC) ******
C
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
      IF(MODELN.EQ.8) CALL WMDPRF(IERR)
      IF(MODELN.EQ.9) CALL WMXPRF(IERR)
      IF(IERR.NE.0) RETURN
C
      IF(NPHMAX.EQ.1) THEN
         NDSIZ  = 1
         NDMIN  = 0
         NDMAX  = 0
         KDSIZ  = 1
         KDMIN  = 0
         KDMAX  = 0
         NDSIZX = 1
      ELSE
         NDSIZ  = NPHMAX
         NDMIN  =-NPHMAX/2+1
         NDMAX  = NPHMAX/2
         KDSIZ  = NPHMAX
         KDMIN  =-NPHMAX/2+1
         KDMAX  = NPHMAX/2
         NDSIZX = 3*NPHMAX/2
      ENDIF
C
      IF(NTHMAX.EQ.1) THEN
         MDSIZ  = 1
         MDMIN  = 0
         MDMAX  = 0
         LDSIZ  = 1
         LDMIN  = 0
         LDMAX  = 0
         MDSIZX = 1
      ELSE
         MDSIZ  = NTHMAX
         MDMIN  =-NTHMAX/2+1
         MDMAX  = NTHMAX/2
         LDSIZ  = NTHMAX
         LDMIN  =-NTHMAX/2+1
         LDMAX  = NTHMAX/2
         MDSIZX = 3*NTHMAX/2
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
         NPHMAX=1
C
         IF(RD.LE.RA.OR.RD.GE.RB) THEN
            IF(MYRANK.EQ.0) WRITE(6,*) 'XX WMXRZC: RD = (RA+RB)/2'
            RD=0.5D0*(RA+RB)
            IF(MYRANK.EQ.0) 
     &           WRITE(6,'(A,1P3E12.4)') 'RA,RB,RD=',RA,RB,RD
         ENDIF
C
         PSIA=RA*RA*BB/(Q0+QA)
         PSIB=SQRT(RB**2/RA**2+(RB**2/RA**2-1.D0)*Q0/QA)*PSIA
         DRHO=SQRT(PSIB/PSIA)/NRMAX
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
         DTH=2.D0*PI/NTHMAX
         DTHG=2.D0*PI/NTHGM
         DO NR=1,NRMAX+1
            IF(NR.EQ.1) THEN
               RS=XR(2)/XRHO(2)
            ELSE
               RS=XR(NR)/XRHO(NR)
            ENDIF
            RSD=QPS(NR)/(BB*RS)
C            WRITE(6,*) 'NR,RS,RSD=',NR,RS,RSD
            DO NTH=1,NTHMAX
               RCOS=COS(DTH*(NTH-1))
               RSIN=SIN(DTH*(NTH-1))
               IF(MODELG.EQ.0) THEN
                  RPS(NTH,NR)  =XR(NR)*RCOS
               ELSE IF(MODELG.EQ.1) THEN
                  RPS(NTH,NR) = RR + XR(NR)*RCOS
               ENDIF
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
               FEDGE=PNS(1)/PN(1)
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1)**PROFN2+FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1)**PROFT2+FEDGE
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
            IF(MODELG.EQ.0) THEN
               RSU(NSU,1)=RA*RCOS
               RSW(NSU,1)=RB*RCOS
            ELSE IF(MODELG.EQ.1) THEN
               RSU(NSU,1)=RR+RA*RCOS
               RSW(NSU,1)=RR+RB*RCOS
            ENDIF
            ZSU(NSU,1)=RA*RSIN
            ZSW(NSU,1)=RB*RSIN
         ENDDO
C
         IF(MODELG.EQ.0) THEN
            RGMIN=-RB*1.01D0
            RGMAX=+RB*1.01D0
         ELSE IF(MODELG.EQ.1) THEN
            RGMIN=RR-RB*1.01D0
            RGMAX=RR+RB*1.01D0
         ENDIF
         ZGMIN=-RB*1.01D0
         ZGMAX= RB*1.01D0

C
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX
            IF(MODELG.EQ.0) THEN
               RRG=RR
            ELSE
               RRG=RPS(NTH,NR)
            ENDIF
         DO NPH=1,NPHMAX
            RPST(NTH,NPH,NR)=RPS(NTH,NR)
            ZPST(NTH,NPH,NR)=ZPS(NTH,NR)
C
            RG11(NTH,NPH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
            RG12(NTH,NPH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                       +DZPSI(NTH,NR)*DZCHI(NTH,NR)
            RG13(NTH,NPH,NR)= 0.D0
            RG22(NTH,NPH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
            RG23(NTH,NPH,NR)= 0.D0
            RG33(NTH,NPH,NR)= RRG**2
            RJ  (NTH,NPH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                             -DRCHI(NTH,NR)*DZPSI(NTH,NR))
C
            BFLD(2,NTH,NPH,NR)=1.D0/RJ(NTH,NPH,NR)
            BFLD(3,NTH,NPH,NR)=RBPS(NR)/RRG**2
C            WRITE(6,*) 'NR,RJ,BFLD2,BFLD3=',NR,RJ(NTH,NPH,NR),
C     &                 BFLD(2,NTH,NPH,NR),BFLD(3,NTH,NPH,NR)
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
      INCLUDE 'wmcomm.inc'
C
C
         IERR=0
C
         NSUMAX=31
         NSWMAX=31
         NPHMAX=1
C
         IF(RD.LT.RA.OR.RD.GT.RB) THEN
            IF(MYRANK.EQ.0) WRITE(6,*) 'XX WMXRZT: RD = (RA+RB)/2'
            RD=0.5D0*(RA+RB)
            IF(MYRANK.EQ.0) 
     &           WRITE(6,'(A,1P3E12.4)') 'RA,RB,RD=',RA,RB,RD
         ENDIF
C
         PSIA=RA*RA*BB/(Q0+QA)
         DRHO=(RB/RA)/NRMAX
C
         DO NR=1,NRMAX+1
            RHOL=DRHO*(NR-1)
            XRHO(NR)=RHOL
            XR(NR)=RA*RHOL
C
C            WRITE(6,*) 'NR,RHO,XR=',NR,XRHO(NR),XR(NR)
C
            IF(RHOL.GT.1.D0) THEN
               QPS(NR)  = QA*RHOL**2
            ELSEIF(RHOMIN.LE.0.D0)THEN
               QPS(NR)  =(Q0-QA)*(1-RHOL**2)+QA
            ELSE
               QSA0    =1/Q0-1/QMIN
               QSAA    =1/QA-1/QMIN
               IF(RHOL.LE.RHOMIN)THEN
                  QPS(NR) =1/(1/Q0-QSA0*(3*RHOL**2/RHOMIN**2
     &                            -2*RHOL**3/RHOMIN**3))
               ELSE
                  QPS(NR) =1/(1/QMIN+3*QSA0*(RHOL-RHOMIN)**2/RHOMIN**2
     &                  +(QSAA-3*QSA0*(1-RHOMIN)**2/RHOMIN**2)
     &                  *(RHOL-RHOMIN)**3/(1-RHOMIN)**3)
               ENDIF
            ENDIF
         ENDDO
C
         DTH=2.D0*PI/NTHMAX
         DTHG=2.D0*PI/NTHGM
         DO NR=1,NRMAX+1
            RSD=RA/(2.D0*PSIA)
            DO NTH=1,NTHMAX
               RCOS=COS(DTH*(NTH-1))
               RSIN=SIN(DTH*(NTH-1))
               IF(MODELG.EQ.0) THEN
                  RPS(NTH,NR)  =XR(NR)*RCOS
               ELSE
                  RPS(NTH,NR) = RR + XR(NR)*RCOS
               ENDIF
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
               FEDGE=PNS(1)/PN(1)
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1)**PROFN2+FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1)**PROFT2+FEDGE
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
            IF(MODELG.EQ.0) THEN
               RSU(NSU,1)=RA*RCOS
               RSW(NSU,1)=RB*RCOS
            ELSE
               RSU(NSU,1)=RR+RA*RCOS
               RSW(NSU,1)=RR+RB*RCOS
            ENDIF
            ZSU(NSU,1)=RA*RSIN
            ZSW(NSU,1)=RB*RSIN
         ENDDO
C
         IF(MODELG.EQ.0) THEN
            RGMIN=-RB*1.01D0
            RGMAX=+RB*1.01D0
         ELSE
            RGMIN=RR-RB*1.01D0
            RGMAX=RR+RB*1.01D0
         ENDIF
         ZGMIN=-RB*1.01D0
         ZGMAX= RB*1.01D0

C
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX
            IF(MODELG.EQ.0) THEN
               RRG=RR
            ELSE
               RRG=RPS(NTH,NR)
            ENDIF
         DO NPH=1,NPHMAX
            RPST(NTH,NPH,NR)=RPS(NTH,NR)
            ZPST(NTH,NPH,NR)=ZPS(NTH,NR)
C
            RG11(NTH,NPH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
            RG12(NTH,NPH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                       +DZPSI(NTH,NR)*DZCHI(NTH,NR)
            RG13(NTH,NPH,NR)= 0.D0
            RG22(NTH,NPH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
            RG23(NTH,NPH,NR)= 0.D0
            RG33(NTH,NPH,NR)= RRG**2
            RJ  (NTH,NPH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                             -DRCHI(NTH,NR)*DZPSI(NTH,NR))
C
            BFLD(2,NTH,NPH,NR)=BB*RR/RRG**2/QPS(NR)
            BFLD(3,NTH,NPH,NR)=BB*RR/RRG**2
C            WRITE(6,*) 'NR,RJ,BFLD2,BFLD3=',NR,RJ(NTH,NPH,NR),
C     &                 BFLD(2,NTH,NPH,NR),BFLD(3,NTH,NPH,NR)
         ENDDO
         ENDDO
         ENDDO
      RETURN
      END
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE WMXRZF(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      CHARACTER KNAMEQSV*80
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV,KNAMEQSV
      DATA NRMAXSV,NTHMAXSV,NSUMAXSV/0,0,0/
      DATA KNAMEQSV/'                                '/
C
      IERR=0
C
      IF(NRMAXSV.EQ.NRMAX.AND.
     &   NTHMAXSV.EQ.NTHMAX.AND.
     &   NSUMAXSV.EQ.NSUMAX.AND.
     &   KNAMEQSV.EQ.KNAMEQ) RETURN
C
      NSUMAX=41
      NSWMAX=NSUMAX
      NPHMAX=1
C
      IF(NTHMAX.LT.4) THEN
         IF(MYRANK.EQ.0) 
     &        WRITE(6,*) 'XX WMXRZF: NTHMAX MUST BE GREATER THAN 4'
         IERR=1
         RETURN
      ENDIF
C
      IF(MYRANK.EQ.0) THEN
         CALL EQLOAD(MODELG,KNAMEQ,IERR)
      ENDIF
      CALL MPBCIA(IERR)
      IF(IERR.EQ.1) RETURN
C
      NRMAXSV=NRMAX
      NTHMAXSV=NTHMAX
      NSUMAXSV=NSUMAX
      KNAMEQSV=KNAMEQ
C
      IF(MYRANK.EQ.0) THEN
         CALL EQSETP
         CALL EQCALQ(NRMAX+1,NTHMAX,NSUMAX,IERR)
C
         CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
C
C         WRITE(6,'(1P6E12.4)') BB,RR,RIP,RA,RKAP,RB
C
         CALL EQGETP(PSIPS,NRMAX+1)
         CALL EQGETR(RPS,DRPSI,DRCHI,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETZ(ZPS,DZPSI,DZCHI,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETBB(BPR,BPZ,BPT,BTP,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETQ(PPS,QPS,RBPS,VPS,RLEN,NRMAX+1)
         CALL EQGETU(RSU,ZSU,RSW,ZSW,NSUMAX)
         CALL EQGETF(RGMIN,RGMAX,ZGMIN,ZGMAX)
         CALL EQGETA(RAXIS,ZAXIS,SAXIS,Q0,QA)
C
C         WRITE(6,'(1P5E12.4)') (RPS(NTH,NRMAX+1),NTH=1,NTHMAX)
C
         CALL EQCALQ(NRMAX+1,NTHGM,NSUMAX,IERR)
         CALL EQGETG(RPSG,ZPSG,NTHGM,NTHGM,NRMAX+1)
C
C         WRITE(6,'(1P5E12.4)') (RPSG(NTH,NRMAX+1),NTH=1,NTHGM)
      ENDIF
      CALL MPBCDA(BB)
      CALL MPBCDA(RR)
      CALL MPBCDA(RIP)
      CALL MPBCDA(RA)
      CALL MPBCDA(RKAP)
      CALL MPBCDA(RDLT)
      CALL MPBCDA(RB)
      CALL MPBCDN(PSIPS,NRMAX+1)
      CALL MPBCDN(RPS,NTHM*(NRMAX+1))
      CALL MPBCDN(DRPSI,NTHM*(NRMAX+1))
      CALL MPBCDN(DRCHI,NTHM*(NRMAX+1))
      CALL MPBCDN(ZPS,NTHM*(NRMAX+1))
      CALL MPBCDN(DZPSI,NTHM*(NRMAX+1))
      CALL MPBCDN(DZCHI,NTHM*(NRMAX+1))
      CALL MPBCDN(BPT,NTHM*(NRMAX+1))
      CALL MPBCDN(BTP,NTHM*(NRMAX+1))
      CALL MPBCDN(PPS,NRMAX+1)
      CALL MPBCDN(QPS,NRMAX+1)
      CALL MPBCDN(RBPS,NRMAX+1)
      CALL MPBCDN(VPS,NRMAX+1)
      CALL MPBCDN(RLEN,NRMAX+1)
      CALL MPBCDN(RSU,NSUMAX)
      CALL MPBCDN(ZSU,NSUMAX)
      CALL MPBCDN(RSW,NSUMAX)
      CALL MPBCDN(ZSW,NSUMAX)
      CALL MPBCDA(RGMIN)
      CALL MPBCDA(RGMAX)
      CALL MPBCDA(ZGMIN)
      CALL MPBCDA(ZGMAX)
      CALL MPBCDA(RAXIS)
      CALL MPBCDA(ZAXIS)
      CALL MPBCDA(SAXIS)
      CALL MPBCDA(Q0)
      CALL MPBCDA(QA)
C
C      WRITE(6,'(1P5E12.4)') RAXIS,ZAXIS,SAXIS,Q0,QA
C
         IF(RD.LT.RA.OR.RD.GT.RB) THEN
            IF(MYRANK.EQ.0) WRITE(6,*) 'XX WMXRZF: RD = (RA+RB)/2'
            RD=0.5D0*(RA+RB)
            IF(MYRANK.EQ.0) 
     &           WRITE(6,'(A,1P3E12.4)') 'RA,RB,RD=',RA,RB,RD
         ENDIF
C
         RMAX=RSU(1,1)
         DO NSU=2,NSUMAX
            RMAX=MAX(RMAX,RSU(NSU,1))
         ENDDO
C
         PSIA=-SAXIS
C
         DO NR=1,NRMAX+1
            ARG=(PSIPS(NR)-SAXIS)/PSIA
            IF(ARG.LE.0.D0) THEN
               XRHO(NR)=0.D0
            ELSE
               XRHO(NR)=SQRT(ARG)
            ENDIF
            XR(NR)=RA*XRHO(NR)
         ENDDO
C
         DO NR=1,NRMAX+1
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            RPST(NTH,NPH,NR)=RPS(NTH,NR)
            ZPST(NTH,NPH,NR)=ZPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO
C
         DO NR=2,NRMAX+1
         DO NTH=1,NTHMAX
         DO NPH=1,NPHMAX
            RG11(NTH,NPH,NR)= (DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2)
            RG12(NTH,NPH,NR)= (DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                        +DZPSI(NTH,NR)*DZCHI(NTH,NR))/XRHO(NR)
            RG13(NTH,NPH,NR)=0.D0
            RG22(NTH,NPH,NR)= (DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2)
     &                        /(XRHO(NR)*XRHO(NR))
            RG23(NTH,NPH,NR)=0.D0
            RG33(NTH,NPH,NR)= RPS(NTH,NR)**2
            RJ  (NTH,NPH,NR)= RPS(NTH,NR)
     &                      *( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                        -DRCHI(NTH,NR)*DZPSI(NTH,NR))/XRHO(NR)
C
C            BFLD(2,NTH,NPH,NR)=1.D0/RJ(NTH,NPH,NR)
C            BFLD(3,NTH,NPH,NR)=RBPS(NR)/RPS(NTH,NR)**2
C
            BPTL=(BPR(NTH,NR)*DZPSI(NTH,NR)
     &           -BPZ(NTH,NR)*DRPSI(NTH,NR))/XRHO(NR)
     &           /SQRT(RG11(NTH,NPH,NR))
     &           /SQRT(RG22(NTH,NPH,NR))
C     &           /RPS(NTH,NR)
            BFLD(2,NTH,NPH,NR)=BPTL
            BFLD(3,NTH,NPH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
C            WRITE(6,'(2I3,1P4E12.4)') NTH,NR,1.D0/RJ(NTH,NPH,NR),
C     &           RBPS(NR)/RPS(NTH,NR)**2,
C     &           BPTL,BPTL*RJ(NTH,NPH,NR)
C     &           BPT,BTP(NTH,NR)/RPS(NTH,NR)
C
C            IF((NR.EQ.2).OR.(NR.EQ.3)) THEN
C            WRITE(6,*) 'NR,NTH,NPH=',NR,NTH,NPH
C            WRITE(6,'(1P3E21.4)') 
C     &           RG11(NTH,NPH,NR),RG12(NTH,NPH,NR),RG13(NTH,NPH,NR)
C            WRITE(6,'(1P3E21.4)') 
C     &           RG22(NTH,NPH,NR),RG23(NTH,NPH,NR),RG33(NTH,NPH,NR)
C            WRITE(6,'(1P3E21.4)') 
C     &           RJ(NTH,NPH,NR),BFLD(2,NTH,NPH,NR),BFLD(3,NTH,NPH,NR)
C            ENDIF
         ENDDO
         ENDDO
         ENDDO
C
         NR=1
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            RG11(NTH,NPH,NR)= RG11(NTH,NPH,2)
            RG12(NTH,NPH,NR)= RG12(NTH,NPH,2)
            RG13(NTH,NPH,NR)= 0.D0
            RG22(NTH,NPH,NR)= RG22(NTH,NPH,2)
            RG23(NTH,NPH,NR)= 0.D0
            RG33(NTH,NPH,NR)= RPST(NTH,NPH,NR)**2
            RJ  (NTH,NPH,NR)= RJ(NTH,NPH,2)
            BPTL=0.D0
C            BFLD(2,NTH,NPH,NR)=1.D0/RJ(NTH,NPH,NR)
C            BFLD(3,NTH,NPH,NR)=RBPS(NR)/RPS(NTH,NR)**2
C            BPT=(BPR(NTH,NR)*DRPSI(NTH,NR)
C     &          +BPZ(NTH,NR)*DZPSI(NTH,NR))/XRHO(NR)
            BFLD(2,NTH,NPH,NR)=BFLD(2,NTH,NPH,2)
            BFLD(3,NTH,NPH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
C            WRITE(6,'(2I3,1P4E12.4)') NTH,NR,1.D0/RJ(NTH,NPH,NR),
C     &           RBPS(NR)/RPS(NTH,NR)**2,
C     &           BPTL,BPTL*RJ(NTH,NPH,NR)
C     &           BPT,BTP(NTH,NR)/RPS(NTH,NR)
C
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
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         BSUPTH=BFLD(2,NTH,NPH,NR)
         BSUPPH=BFLD(3,NTH,NPH,NR)
         BABS=SQRT(     RG22(NTH,NPH,NR)*BSUPTH*BSUPTH*XRHO(NR)**2
     &            +2.D0*RG23(NTH,NPH,NR)*BSUPTH*BSUPPH*XRHO(NR)
     &            +     RG33(NTH,NPH,NR)*BSUPPH*BSUPPH)
         BPST(NTH,NPH,NR)=BABS
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
