C     $Id$
C
C****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPFMFL(NS)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      IF(MODELV.EQ.1) RETURN
C
      NPMAX=100
      NTHMAX=100
      PMAX=10.D0
C
      RNE0 = RN(NS)
      TE0  = RTPR(NS)
      PTH0 = SQRT(TE0*1.D3*AEE*AMP*PA(NS))
C
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
C
      DELP=PMAX/NPMAX
      DELTH=PI/NTHMAX
C
      DO 10 NP=1,NPMAX
         PM(NP)=DELP*(NP-0.5D0)
         PG(NP)=DELP* NP
         RGMM(NP)=SQRT(1+PTH0W*PM(NP)**2) 
         RGMG(NP)=SQRT(1+PTH0W*PG(NP)**2)
   10 CONTINUE
C
      DO 20 NTH=1,NTHMAX
         THM=DELTH*(NTH-0.5D0)
         THG=DELTH* NTH
         TSNM(NTH) = SIN(THM)
         TSNG(NTH) = SIN(THG)
         TCSM(NTH) = COS(THM)
         TCSG(NTH) = COS(THG)
         TTNM(NTH) = TAN(THM)
         TTNG(NTH) = TAN(THG)
   20 CONTINUE
C
      FACT = 1.D0/SQRT(2.D0*PI)**3
      DO 100 NP=1,NPMAX
         PML=PM(NP)
         EX =PML*PML/2.D0
      DO 100 NTH=1,NTHMAX
         FM(NP,NTH) = FACT * EXP(-EX)
  100 CONTINUE
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
      RNE0  =FPDATA(1)
      TE0   =FPDATA(2)
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
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            X=RHOFP*RA
            CALL SPL1DF(X,Y,RFP,URFP(1,1,NP,NTH),NRMAX+2,IERR)
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
            IF(COS(TH).GT.0.D0) THEN
               X=   ASIN(SQRT(1/PSIS)*SIN(TH))
            ELSE
               X=PI-ASIN(SQRT(1/PSIS)*SIN(TH))
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
      CHARACTER KNAM*72
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : FPDATA FILE NAME : ',KNAMFP
      READ(5,'(A72)',ERR=1,END=900) KNAM
      IF(KNAM(1:2).NE.'/ ') KNAMFP=KNAM
C
      INQUIRE(FILE=KNAMFP,EXIST=LEX)
      IF(LEX) THEN
         OPEN(21,FILE=KNAMFP,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMFP,') IS ASSIGNED FOR INPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         STOP
C        GOTO 1
      ELSE
         WRITE(6,*) 'XX FILE (',KNAMFP,') NOT FOUND' 
         STOP
C        GOTO 1
      ENDIF
C
   30 READ(21) RNE0,TE0,PTH0
      READ(21) DELR,DELP,DELTH
      READ(21) NRMAX,NPMAX,NTHMAX
      READ(21) (((FP(NTH,NP,NR),NTH=1,NTHMAX),NP=1,NPMAX),NR=1,NRMAX)
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'RNE0,TE0,PTH0 =',RNE0,TE0,PTH0
      WRITE(6,*) 'DELR,DELP,DELTH =',DELR,DELP,DELTH
      WRITE(6,*) 'NRMAX,NPMAX,NTHMAX =',NRMAX,NPMAX,NTHMAX
C
      IF(NRMAX.EQ.1) THEN
         DELR=RA
      ENDIF
C
      RFP(1)=0.D0
      DO NR=1,NRMAX
         RFP(NR+1)=DELR*(NR-0.5D0)
      ENDDO
      RFP(NRMAX+2)=DELR*NRMAX
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
      FPDATA(1)=RNE0
      FPDATA(2)=TE0
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
