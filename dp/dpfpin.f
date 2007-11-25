C     $Id$
C
C****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPFMFL(NS,ID)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      NPMAX=100
      NTHMAX=100
      PMAX=10.D0
C
      PN0 = RN(NS)
      PT0 = (RTPR(NS)+2*RTPP(NS))/3.D0
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
         RGMM(NP)=SQRT(1+PTH0W*PM(NP)**2) 
         RGMG(NP)=SQRT(1+PTH0W*PG(NP)**2)
      ENDDO
C
      DO NTH=1,NTHMAX
         THM=DELTH*(NTH-0.5D0)
         THG=DELTH* NTH
         TSNM(NTH) = SIN(THM)
         TSNG(NTH) = SIN(THG)
         TCSM(NTH) = COS(THM)
         TCSG(NTH) = COS(THG)
         TTNM(NTH) = TAN(THM)
         TTNG(NTH) = TAN(THG)
      ENDDO
C
      IF(ID.EQ.0) THEN
         FACT = 1.D0/SQRT(2.D0*PI)**3
         SUM=0.D0
         DO NP=1,NPMAX
            PML=PM(NP)
            EX =PML*PML/2.D0
            DO NTH=1,NTHMAX
               FM(NP,NTH) = FACT * EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP)*PM(NP)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP*DELTH
C         WRITE(6,*) 'SUM=',SUM
C         WRITE(6,*) 'FM(10,10)=',FM(10,10)
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
C         WRITE(6,*) 'SUM=',SUM
         FACTOR=1.D0/SUM
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               FM(NP,NTH) = FACTOR*FM(NP,NTH)
            ENDDO
         ENDDO
C         WRITE(6,*) 'FM(10,10)=',FM(10,10)
      ENDIF
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
      DIMENSION THT(NTHM+2),FPT(NTHM+2),FPX(NTHM+2)
      DIMENSION U2(4,NTHM+2)
C
      PN0   =FPDATA(1)
      PT0   =FPDATA(2)
      PTH0  =FPDATA(3)
      DELR  =FPDATA(4)
      DELP  =FPDATA(5)
      DELTH =FPDATA(6)
      NRMAX =NFPDAT(1)
      NPMAX =NFPDAT(2)
      NTHMAX=NFPDAT(3)
C
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
C
      DO 110 NP=1,NPMAX
         PM(NP)=DELP*(NP-0.5D0)
         PG(NP)=DELP* NP
         RGMM(NP)=SQRT(1+PTH0W*PM(NP)**2) 
         RGMG(NP)=SQRT(1+PTH0W*PG(NP)**2)
  110 CONTINUE
C
      DO 120 NTH=1,NTHMAX
         THM=DELTH*(NTH-0.5D0)
         THG=DELTH* NTH
         TSNM(NTH) = SIN(THM)
         TSNG(NTH) = SIN(THG)
         TCSM(NTH) = COS(THM)
         TCSG(NTH) = COS(THG)
         TTNM(NTH) = TAN(THM)
         TTNG(NTH) = TAN(THG)
  120 CONTINUE
C
      RHON=RHON_LOC
      CALL PLBMIN(RHON,BMINL)
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
C            WRITE(6,'(6I8)') NRMAX,NP,NTH,NRM,NPM,NTHM
            CALL SPL1DF(RHON,Y,RFP,URFP(1,1,NP,NTH),NRMAX+2,IERR)
            FPT(NTH+1)=Y
         ENDDO
C
         PSIS=BABS/BMINL
C
         THT(1)=0.D0
         DO NTH=1,NTHMAX
            THT(NTH+1)=DELTH*(NTH-0.5D0)
         ENDDO
         THT(NTHMAX+2)=PI
         FPT(1)=(9.D0*FPT(2)-FPT(3))/8.D0
         FPT(NTHMAX+2)=(9.D0*FPT(NTHMAX+1)-FPT(NTHMAX))/8.D0         
         FPX(1)=0.D0
         FPX(NTHMAX+2)=0.D0
         CALL SPL1D(THT,FPT,FPX,U2,NTHMAX+2,3,IERR)            
         DO NTH2=1,NTHMAX
            TH=DELTH*(NTH2-0.5D0)
            XX=SQRT(1/PSIS)*SIN(TH)
            IF(XX.GE.1.D0) THEN
               X=0.5D0*PI
C               write(6,'(A,1P4E12.4)') 
C     &              'X,PSIS,TH,SIN=',X,PSIS,TH,SIN(TH)
            ELSE
               IF(COS(TH).GT.0.D0) THEN
                  X=   ASIN(SQRT(1/PSIS)*SIN(TH))
               ELSE
                  X=PI-ASIN(SQRT(1/PSIS)*SIN(TH))
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
      DIMENSION FPR(NRM+2),FPX(NRM+2)
C
      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF
      REWIND(21)

      READ(21) PN0,PT0,PTH0
      READ(21) DELR,DELP,DELTH,RMIN,RMAX
      READ(21) NRMAX,NPMAX,NTHMAX
      READ(21) (((FP(NTH,NP,NR),NTH=1,NTHMAX),NP=1,NPMAX),NR=1,NRMAX)
      CLOSE(21)
C
      PN0=0.2D0*PN0
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'PN0,PT0,PTH0 =',PN0,PT0,PTH0
      WRITE(6,*) 'DELR,DELP,DELTH =',DELR,DELP,DELTH
      WRITE(6,*) 'NRMAX,NPMAX,NTHMAX =',NRMAX,NPMAX,NTHMAX
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
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            IF(NRMAX.EQ.1)THEN
               FPR(1)=4.D0*FP(NTH,NP,1)/3.D0
            ELSE
               FPR(1)=(9.D0*FP(NTH,NP,1)-FP(NTH,NP,2))/8.D0
            ENDIF
            DO NR=1,NRMAX
               FPR(NR+1)=FP(NTH,NP,NR)
            ENDDO
            FPR(NRMAX+2)=0.D0
            FPX(1)=0.D0
            CALL SPL1D(RFP,FPR,FPX,URFP(1,1,NP,NTH),NRMAX+2,1,IERR)
         ENDDO
      ENDDO
C
      FPDATA(1)=PN0
      FPDATA(2)=PT0
      FPDATA(3)=PTH0
      FPDATA(4)=DELR
      FPDATA(5)=DELP
      FPDATA(6)=DELTH
      NFPDAT(1)=NRMAX
      NFPDAT(2)=NPMAX
      NFPDAT(3)=NTHMAX
C
  900 RETURN
      END
