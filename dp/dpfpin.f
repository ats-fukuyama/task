C     $Id$
C
C****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPFMFL(NS)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      RHON_MIN=0.D0
      RHON_MAX=1.D0
C
      PN0 = PN(NS)
      PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
      PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))
C
      RN0 = RN(NS)
      TPR = RTPR(NS)
      TPP = RTPP(NS)
      RT0 = (TPR+2.D0*TPP)/3.D0
C
      IF(ID.EQ.0) THEN
         PTH0W=0.D0
      ELSE
         PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      ENDIF
C
      DELP=PMAX/NPMAX
      DELTH=PI/NTHMAX
C
      DO NP=1,NPMAX
         PM(NP)=DELP*(NP-0.5D0)
         PG(NP)=DELP* NP
      ENDDO
C
      DO NTH=1,NTHMAX
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO
C
      IF(ID.EQ.0) THEN
         TNPR=TPR/RT0
         TNPP=TPP/RT0
         SUM=0.D0
         DO NP=1,NPMAX
            PML=PM(NP)
            DO NTH=1,NTHMAX
               PPP=PML*TSNM(NTH)
               PPR=PML*TCSM(NTH)
               EX=-(PPR**2/TNPR+PPP**2/TNPP)
               FM(NP,NTH) = EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP)*PM(NP)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP*DELTH
      ELSE
         TN00=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
         TNPR=TPR*1.D3*AEE/(AMP*PA(NS)*VC**2)
         TNPP=TPP*1.D3*AEE/(AMP*PA(NS)*VC**2)
         SUM=0.D0
         DO NP=1,NPMAX
            PML=PM(NP)
            DO NTH=1,NTHMAX
               PPP=PML*TSNM(NTH)
               PPR=PML*TCSM(NTH)
               EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP)
     &           -SQRT(1.D0/TN00**2)
               FM(NP,NTH) = EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP)*PM(NP)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP*DELTH
      ENDIF

      FACTOR=1.D0/SUM
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            FM(NP,NTH) = FACTOR*FM(NP,NTH)
         ENDDO
      ENDDO
      RETURN
      END                     
C
C     ****** SET LOCAL VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPFPFL(NS)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION THT(NTHM+2),FPT(NTHM+2),FPTX(NTHM+2)
      DIMENSION U2(4,NTHM+2)
      DIMENSION FPR(NRM+2),FPRX(NRM+2)
C
      DO NSA1=1,NSAMAX
         IF(NS.EQ.NS_NSA(NSA1)) THEN
            NSA=NSA1
            GO TO 1
         ENDIF
      ENDDO
      RETURN

    1 CONTINUE

      PN0   =RNFP0(NSA)
      PT0   =RTFP0(NSA)
      PTH0  =SQRT(PT0*1.D3*AEE*AMFP(NSA))
      DELR  =FPDATA(1)
      DELP  =FPDATA(2)
      DELTH =FPDATA(3)
      NRMAX =NFPDAT(1)
      NPMAX =NFPDAT(2)
      NTHMAX=NFPDAT(3)
C
      PTH0W=(PTH0/(AMFP(NSA)*VC))**2
C
      DO NP=1,NPMAX
         PM(NP)=DELP*(NP-0.5D0)
         PG(NP)=DELP* NP
      ENDDO
C
      DO NTH=1,NTHMAX
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO
C
      DO NP=1,NPMAX
        DO NTH=1,NTHMAX
            DO NR=1,NRMAX
               FP(NTH,NP,NR)=FNS(NTH,NP,NR,NSA)
            ENDDO
            IF(RHON_MIN.EQ.0.D0)THEN
               FPR(1)=(9.D0*FP(NTH,NP,1)-FP(NTH,NP,2))/8.D0
               FPRX(1)=0.D0
               ID=1
            ELSE
               FPR(1)=(4.D0*FP(NTH,NP,1)-FP(NTH,NP,2))/3.D0
               ID=0
            ENDIF
            DO NR=1,NRMAX
               FPR(NR+1)=FP(NTH,NP,NR)
            ENDDO
            FPR(NRMAX+2)=(4.D0*FP(NTH,NP,NRMAX+1)-FP(NTH,NP,NRMAX))/3.D0
            CALL SPL1D(RFP,FPR,FPRX,URFP(1,1,NP,NTH),NRMAX+2,ID,IERR)
         ENDDO
      ENDDO
C
      RHON=RHON_LOC
      CALL PLBMIN(RHON,BMINL)
      PSIS=BABS/BMINL

      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            CALL SPL1DF(RHON,Y,RFP,URFP(1,1,NP,NTH),NRMAX+2,IERR)
            FPT(NTH+1)=Y
         ENDDO
         THT(1)=0.D0
         DO NTH=1,NTHMAX
            THT(NTH+1)=DELTH*(NTH-0.5D0)
         ENDDO
         THT(NTHMAX+2)=PI
         FPT(1)=(9.D0*FPT(2)-FPT(3))/8.D0
         FPT(NTHMAX+2)=(9.D0*FPT(NTHMAX+1)-FPT(NTHMAX))/8.D0         
         FPTX(1)=0.D0
         FPTX(NTHMAX+2)=0.D0
         CALL SPL1D(THT,FPT,FPTX,U2,NTHMAX+2,3,IERR)            
         DO NTH2=1,NTHMAX
            TH=DELTH*(NTH2-0.5D0)
            XX=SQRT(1/PSIS)*SIN(TH)
            IF(XX.GE.1.D0) THEN
               X=0.5D0*PI
            ELSE
               IF(COS(TH).GT.0.D0) THEN
                  X=   ASIN(XX)
               ELSE
                  X=PI-ASIN(XX)
               ENDIF
            ENDIF
            CALL SPL1DF(X,Y,THT,U2,NTHMAX+2,IERR)
            FM(NP,NTH2)=Y
         ENDDO
      ENDDO
      RETURN
      END
C
C****** LOAD VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPLDFP
C
      INCLUDE 'dpcomm.inc'
C      
      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF
      REWIND(21)

      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELP,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX),
     &                NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'NRMAX,NPMAX,NTHMAX,NSAMAX =',NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(6,'(A,1P5E12.4)') 'DELR/P/TH,RMIN/RMAX =',
     &                         DELR,DELP,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         WRITE(6,*) 'NSA,NS =',NSA,NS_NSA(NSA)
         WRITE(6,'(A,1P4E12.4)') 'AEFP,AMFP,RNFP0,RTFP0=',
     &        AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
      ENDDO
C
      RHON_MIN=RMIN
      RHON_MAX=RMAX

      IF(NRMAX.EQ.1) THEN
         DELR=1.D0
      ENDIF
C
      RFP(1)=RMIN
      DO NR=1,NRMAX
         RFP(NR+1)=RMIN+DELR*(NR-0.5D0)
      ENDDO
      RFP(NRMAX+2)=RMIN+DELR*NRMAX
C
      FPDATA(1)=DELR
      FPDATA(2)=DELP
      FPDATA(3)=DELTH
      NFPDAT(1)=NRMAX
      NFPDAT(2)=NPMAX
      NFPDAT(3)=NTHMAX
C
  900 RETURN
      END
C
C****** LOAD Maxwellian VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPLDFM(ID,NCHMAX)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      RHON_MIN=0.D0
      RHON_MAX=1.D0
C
      IF(NRMAX.EQ.1) THEN
         DELR=0.1D0
      ELSE
         DELR=(RMAX-RMIN)/(NRMAX-1)
      ENDIF

      DELCH=2.D0*PI/NCHMAX
      DO NCH=1,NCHMAX
         CHI(NCH)=DELCH*(NCH-1)
         CCHI(NCH)=COS(CHI(NCH))
         SCHI(NCH)=SIN(CHI(NCH))
      ENDDO

      DELP=PMAX/NPMAX
      DELTH=PI/NTHMAX
C
      DO NR=1,NRMAX
         RM(NR)=RMIN+DELR*(NR-1)
      ENDDO
C
      DO NP=1,NPMAX
         PM(NP)=DELP*(NP-0.5D0)
         PG(NP)=DELP* NP
      ENDDO
C
      DO NTH=1,NTHMAX
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO
C
      DO NS=1,NSAMAX
         PN0 = PN(NS)
         PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
         PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))

         NS_NSA(NS)=NS
         RNFP0(NS)=PN0
         RTFP0(NS)=PT0
         AMFP(NS)=PA(NS)
C
         IF(ID.EQ.0) THEN
            PTH0W=0.D0
         ELSE
            PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
         ENDIF
C
         DO NR=1,NRMAX
            RHON=RM(NR)
            CALL PLPROF(RHON)

            RN0 = RN(NS)
            TPR = RTPR(NS)
            TPP = RTPP(NS)
C
            IF(ID.EQ.0) THEN
               TNPR=TPR/PT0
               TNPP=TPP/PT0
               SUM=0.D0
               DO NP=1,NPMAX
                  PML=PM(NP)
                  DO NTH=1,NTHMAX
                     PPP=PML*TSNM(NTH)
                     PPR=PML*TCSM(NTH)
                     EX=0.5D0*(PPR**2/TNPR+PPP**2/TNPP)
                     FM(NP,NTH) = EXP(-EX)
                     SUM=SUM+FM(NP,NTH)*PM(NP)*PM(NP)*TSNM(NTH)
                  ENDDO
               ENDDO
               SUM=SUM*2.D0*PI*DELP*DELTH
            ELSE
               TN00=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
               TNPR=TPR*1.D3*AEE/(AMP*PA(NS)*VC**2)
               TNPP=TPP*1.D3*AEE/(AMP*PA(NS)*VC**2)
               SUM=0.D0
               DO NP=1,NPMAX
                  PML=PM(NP)
                  DO NTH=1,NTHMAX
                     PPP=PML*TSNM(NTH)
                     PPR=PML*TCSM(NTH)
                     EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP)
     &                 -SQRT(1.D0/TN00**2)
                     FM(NP,NTH) = EXP(-EX)
                     SUM=SUM+FM(NP,NTH)*PM(NP)*PM(NP)*TSNM(NTH)
                  ENDDO
               ENDDO
               SUM=SUM*2.D0*PI*DELP*DELTH
            ENDIF

            FACTOR=RN0/(PN0*SUM)
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNS(NTH,NP,NR,NS) = FACTOR*FM(NP,NTH)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      FPDATA(1)=DELR
      FPDATA(2)=DELP
      FPDATA(3)=DELTH
      NFPDAT(1)=NRMAX
      NFPDAT(2)=NPMAX
      NFPDAT(3)=NTHMAX
C
  900 RETURN
      END
C
C****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA with radial dep. ******
C
      SUBROUTINE DPFMFLR(NS)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      NRMAX=3
      NCHMAX=16
      RHON_MIN=0.D0
      RHON_MAX=1.D0

      DELR=DR
      DELCH=2.D0*PI/NCHMAX

      DO NR=1,3
         SELECT CASE(NR)
         CASE(1)
            RS=R1-DR
         CASE(2)
            RS=R1
         CASE(3)
            RS=R1+DR
         END SELECT
         RM(NR)=RS
      ENDDO

      DO NCH=1,NCHMAX
         CHI(NCH)=DELCH*(NCH-1)
         CCHI(NCH)=COS(CHI(NCH))
         SCHI(NCH)=SIN(CHI(NCH))
      ENDDO

      NPMAX=50
      NTHMAX=50
      PMAX=10.D0

      DELP=PMAX/NPMAX
      DELTH=PI/NTHMAX

      IF(ID.EQ.0) THEN
         PTH0W=0.D0
      ELSE
         PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      ENDIF

      DO NP=1,NPMAX
         PM(NP)=DELP*(NP-0.5D0)
         PG(NP)=DELP* NP
      ENDDO
C
      DO NTH=1,NTHMAX
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO

      RHON=0.D0
      CALL PLPROF(RHON)

      PN0 = RN(NS)
      PT0 = (RTPR(NS)+2*RTPP(NS))/3.D0
      PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))

      DO NR=1,NRMAX
         RHON=RM(NR)/RA

         CALL PLPROF(RHON)

         RN0 = RN(NS)
         TPR = RTPR(NS)
         TPP = RTPP(NS)
         RT0 = (TPR+2.D0*TPP)/3.D0

         IF(ID.EQ.0) THEN
            FACT = 1.D0/SQRT(2.D0*PI)**3
            SUM=0.D0
            DO NP=1,NPMAX
               PML=PM(NP)
               EX =PML*PML/2.D0
               DO NTH=1,NTHMAX
                  FP(NTH,NP,NR) = FACT * EXP(-EX)
                  SUM=SUM+FP(NTH,NP,NR)*PM(NP)*PM(NP)*TSNM(NTH)
               ENDDO
            ENDDO
            SUM=SUM*2.D0*PI*DELP*DELTH
         ELSE
            TN00=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
            TNPR=TPR*1.D3*AEE/(AMP*PA(NS)*VC**2)
            TNPP=TPP*1.D3*AEE/(AMP*PA(NS)*VC**2)
            SUM=0.D0
            DO NP=1,NPMAX
               PML=PM(NP)
               DO NTH=1,NTHMAX
                  PPP=PML*TSNM(NTH)
                  PPR=PML*TCSM(NTH)
                  EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP)
     &                 -SQRT(1.D0/TN00**2)
                  FP(NTH,NP,NR) = EXP(-EX)
                  SUM=SUM+FP(NTH,NP,NR)*PM(NP)*PM(NP)*TSNM(NTH)
               ENDDO
            ENDDO
            SUM=SUM*2.D0*PI*DELP*DELTH
         ENDIF
         FACTOR=1.D0/SUM
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               FP(NTH,NP,NR) = FACTOR*FP(NTH,NP,NR)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END                     
