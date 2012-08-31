C     $Id$
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE WMSETG(IERR)
C
      INCLUDE 'wmcomm.inc'
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
         NDSIZ  = 1
         NDMIN  = 0
         NDMAX  = 0
         KDSIZ  = 1
         KDMIN  = 0
         KDMAX  = 0
         NDSIZX = 1
      ELSE
         NDSIZ  = NHHMAX
         NDMIN  =-NHHMAX/2+1
         NDMAX  = NHHMAX/2
         KDSIZ  = NHHMAX
         KDMIN  =-NHHMAX/2+1
         KDMAX  = NHHMAX/2
         NDSIZX = 3*NHHMAX/2
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
         NHHMAX=1
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
         DO NTH=1,NTHMAX
            RRG=RPS(NTH,NR)
         DO NHH=1,NHHMAX
            RPST(NTH,NHH,NR)=RPS(NTH,NR)
            ZPST(NTH,NHH,NR)=ZPS(NTH,NR)
C
            RG11(NTH,NHH,NR)= DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2
            RG12(NTH,NHH,NR)= DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                       +DZPSI(NTH,NR)*DZCHI(NTH,NR)
            RG13(NTH,NHH,NR)= 0.D0
            RG22(NTH,NHH,NR)= DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2
            RG23(NTH,NHH,NR)= 0.D0
            RG33(NTH,NHH,NR)= RRG**2
            RJ  (NTH,NHH,NR)= RRG*( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                             -DRCHI(NTH,NR)*DZPSI(NTH,NR))
C
            BFLD(2,NTH,NHH,NR)=1.D0/RJ(NTH,NHH,NR)
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
      USE plfile_prof_mod
      INCLUDE 'wmcomm.inc'
C
C
         IERR=0
C
         NSUMAX=31
         NSWMAX=31
         NHHMAX=1

C         CALL plfile_prof_read(modeln,modelq,ierr)
C
         PSIPA=RA*RA*BB/(Q0+QA)
         DRHO=(RB/RA)/NRMAX
C
         DO NR=1,NRMAX+1
            RHOL=DRHO*(NR-1)
            XRHO(NR)=RHOL
            XR(NR)=RA*RHOL
            CALL PLQPRF(RHOL,QPS(NR))
         ENDDO
C
         DTH=2.D0*PI/NTHMAX
         DTHG=2.D0*PI/NTHGM
         DO NR=1,NRMAX+1
            RSD=RA/(2.D0*PSIPA)
            DO NTH=1,NTHMAX
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
         DO NTH=1,NTHMAX
            RRG=RPS(NTH,NR)
         DO NHH=1,NHHMAX
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
      RETURN
      END
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE WMXRZF(IERR)
C
      USE plfile_prof_mod
      INCLUDE 'wmcomm.inc'
C
      CHARACTER*(80) LINE
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV
      DATA NRMAXSV,NTHMAXSV,NSUMAXSV/0,0,0/
C
      IERR=0
C
      IF(NRMAXSV.EQ.NRMAX.AND.
     &   NTHMAXSV.EQ.NTHMAX.AND.
     &   NSUMAXSV.EQ.NSUMAX.AND.
     &   KNAMEQ_SAVE.EQ.KNAMEQ) RETURN
C
      NSUMAX=41
      NSWMAX=NSUMAX
      NHHMAX=1
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
      KNAMEQ_SAVE=KNAMEQ
C
      IF(MYRANK.EQ.0) THEN
         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHMAX
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
C
         CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
C
C         WRITE(6,'(1P6E12.4)') BB,RR,RIP,RA,RKAP,RB
C
         CALL EQGETP(RHOT,PSIP,NRMAX+1)
         CALL EQGETR(RPS,DRPSI,DRCHI,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETZ(ZPS,DZPSI,DZCHI,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETBB(BPR,BPZ,BPT,BTP,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETQ(PPS,QPS,RBPS,VPS,RLEN,NRMAX+1)
         CALL EQGETU(RSU,ZSU,RSW,ZSW,NSUMAX)
         CALL EQGETF(RGMIN,RGMAX,ZGMIN,ZGMAX)
         CALL EQGETA(RAXIS,ZAXIS,PSIPA,PSITA,Q0,QA)
C
C         WRITE(6,'(A,1P2E12.4)') 'PSIPA,PSITA=',PSIPA,PSITA
C         WRITE(6,'(1P5E12.4)') (RPS(NTH,NRMAX+1),NTH=1,NTHMAX)
C
         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHGM
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
         CALL EQGETG(RPSG,ZPSG,NTHGM,NTHGM,NRMAX+1)
C
C         WRITE(6,'(1P5E12.4)') (RPSG(NTH,NRMAX+1),NTH=1,NTHGM)

C         CALL plfile_prof_read(modeln,modelq,ierr)

      ENDIF
      CALL MPBCDA(BB)
      CALL MPBCDA(RR)
      CALL MPBCDA(RIP)
      CALL MPBCDA(RA)
      CALL MPBCDA(RKAP)
      CALL MPBCDA(RDLT)
      CALL MPBCDA(RB)
      CALL MPBCDN(RHOT,NRMAX+1)
      CALL MPBCDN(PSIP,NRMAX+1)
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
      CALL MPBCDA(PSIPA)
      CALL MPBCDA(PSITA)
      CALL MPBCDA(Q0)
      CALL MPBCDA(QA)
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
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            RPST(NTH,NHH,NR)=RPS(NTH,NR)
            ZPST(NTH,NHH,NR)=ZPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO
C
         DO NR=2,NRMAX+1
         DO NTH=1,NTHMAX
         DO NHH=1,NHHMAX
            RG11(NTH,NHH,NR)= (DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2)
C     &                       *(2.D0*PI)**2
            RG12(NTH,NHH,NR)= (DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                        +DZPSI(NTH,NR)*DZCHI(NTH,NR))/XRHO(NR)
C     &                       * 2.D0*PI
            RG13(NTH,NHH,NR)=0.D0
            RG22(NTH,NHH,NR)= (DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2)
     &                        /(XRHO(NR)*XRHO(NR))
            RG23(NTH,NHH,NR)=0.D0
            RG33(NTH,NHH,NR)= RPS(NTH,NR)**2
            RJ  (NTH,NHH,NR)= RPS(NTH,NR)
     &                      *( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                        -DRCHI(NTH,NR)*DZPSI(NTH,NR))/XRHO(NR)
C     &                      / 2.D0*PI
C
C            BFLD(2,NTH,NHH,NR)=1.D0/(2.D0*PI*RJ(NTH,NHH,NR))
C            BFLD(3,NTH,NHH,NR)=RBPS(NR)/RPS(NTH,NR)**2
C     &                        /(2.D0*PI)
C
            BPTL=(BPR(NTH,NR)*DZPSI(NTH,NR)
     &           -BPZ(NTH,NR)*DRPSI(NTH,NR))/XRHO(NR)
     &           /SQRT(RG11(NTH,NHH,NR))
     &           /SQRT(RG22(NTH,NHH,NR))
C     &                       * 2.D0*PI
C
            BFLD(2,NTH,NHH,NR)=BPTL
            BFLD(3,NTH,NHH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
C
C            IF(NTH.EQ.1) WRITE(6,'(2I3,1P6E12.4)') 
C     &           NR,NTH,1.D0/(2.D0*PI*RJ(NTH,NHH,NR)),
C     &           BPTL,
C     &           RBPS(NR)/RPS(NTH,NR)**2/(2.D0*PI),
C     &           BTP(NTH,NR)/RPS(NTH,NR),
C     &           BPZ(NTH,NR),XRHO(NR)
C
C            IF((NR.EQ.2).OR.(NR.EQ.3)) THEN
C            WRITE(6,*) 'NR,NTH,NHH=',NR,NTH,NHH
C            WRITE(6,'(1P3E21.4)') 
C     &           RG11(NTH,NHH,NR),RG12(NTH,NHH,NR),RG13(NTH,NHH,NR)
C            WRITE(6,'(1P3E21.4)') 
C     &           RG22(NTH,NHH,NR),RG23(NTH,NHH,NR),RG33(NTH,NHH,NR)
C            WRITE(6,'(1P3E21.4)') 
C     &           RJ(NTH,NHH,NR),BFLD(2,NTH,NHH,NR),BFLD(3,NTH,NHH,NR)
C            ENDIF
         ENDDO
         ENDDO
         ENDDO
C
         NR=1
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            RG11(NTH,NHH,NR)= RG11(NTH,NHH,2)
            RG12(NTH,NHH,NR)= RG12(NTH,NHH,2)
            RG13(NTH,NHH,NR)= 0.D0
            RG22(NTH,NHH,NR)= RG22(NTH,NHH,2)
            RG23(NTH,NHH,NR)= 0.D0
            RG33(NTH,NHH,NR)= RPST(NTH,NHH,NR)**2
            RJ  (NTH,NHH,NR)= RJ(NTH,NHH,2)
            BPTL=0.D0
C
C            BFLD(2,NTH,NHH,NR)=1.D0/(2.D0*PI*RJ(NTH,NHH,NR))
C            BFLD(3,NTH,NHH,NR)=RBPS(NR)/RPS(NTH,NR)**2
C     &                        /(2.D0*PI)
C
            BFLD(2,NTH,NHH,NR)=BFLD(2,NTH,NHH,2)
            BFLD(3,NTH,NHH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
C
C            WRITE(6,'(2I3,1P4E12.4)') NTH,NR,1.D0/RJ(NTH,NHH,NR),
C     &           RBPS(NR)/RPS(NTH,NR)**2,
C     &           BPTL,BPTL*RJ(NTH,NHH,NR)
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
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
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
      ENDDO
C
      RETURN
      END
