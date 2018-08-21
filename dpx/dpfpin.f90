MODULE dpfpin

CONTAINS

!****** LOAD VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPLDFP(IERR)

      USE dpcomm
      USE commpi
      USE libmpi
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: IERR
      CHARACTER(LEN=80),SAVE::  KNAMFP_SAVE=' '
      REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: FNSDP
      INTEGER:: NSA,NTH,NP,NR,NS
      
      
      IERR=0
      IF(KNAMFP.EQ.KNAMFP_SAVE) RETURN

      IF(nrank.eq.0) THEN
         CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
            RETURN
         ENDIF
         REWIND(21)

         READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
         CALL DPFP_ALLOCATE
         READ(21) DELR,DELTH,RMIN,RMAX
         DO NSA=1,NSAMAX
            READ(21) NS_NSA(NSA)
            READ(21) DELP(NSA)
            READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
            READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
         ENDDO
         CLOSE(21)
      ENDIF
      
      CALL mtx_broadcast1_integer(NRMAX)
      CALL mtx_broadcast1_integer(NPMAX)
      CALL mtx_broadcast1_integer(NTHMAX)
      CALL mtx_broadcast1_integer(NSAMAX)
      CALL mtx_broadcast1_real8(DELR)
      CALL mtx_broadcast1_real8(DELTH)
      CALL mtx_broadcast1_real8(RMIN)
      CALL mtx_broadcast1_real8(RMAX)
      CALL mtx_broadcast_integer(NS_NSA,NSAMAX)
      CALL mtx_broadcast_real8(DELP,NSAMAX)
      CALL mtx_broadcast_real8(AEFP,NSAMAX)
      CALL mtx_broadcast_real8(AMFP,NSAMAX)
      CALL mtx_broadcast_real8(RNFP0,NSAMAX)
      CALL mtx_broadcast_real8(RTFP0,NSAMAX)

      DO NSA=1,NSAMAX
         CALL mtx_broadcast_real8(FNS(1:NTHMAX,1:NPMAX,1:NRMAX,NSA), &
                                  NTHMAX*NPMAX*NRMAX)
      ENDDO

      IF(nrank.eq.0) THEN
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'NRMAX,NPMAX,NTHMAX,NSAMAX =',NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(6,'(A,1P4E12.4)') 'DELR/TH,RMIN/MAX =', &
                               DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         WRITE(6,*) 'NSA,NS =',NSA,NS_NSA(NSA)
         WRITE(6,'(A,1P5E12.4)') 'AE,AM,RN0,RT0,DELP=', &
               AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA),DELP(NSA)
      ENDDO
      ENDIF
!     
      DELTH=PI/NTHMAX

      DO NSA=1,NSAMAX
      DO NP=1,NPMAX
         PM(NP,NSA)=DELP(NSA)*(NP-0.5D0)
         PG(NP,NSA)=DELP(NSA)* NP
      ENDDO
      ENDDO

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

      RHON_MIN=RMIN
      RHON_MAX=RMAX

      IF(NRMAX.EQ.1) THEN
         DELR=1.D0
      ENDIF

      RFP(1)=RMIN
      DO NR=1,NRMAX
         RFP(NR+1)=RMIN+DELR*(NR-0.5D0)
      ENDDO
      RFP(NRMAX+2)=RMIN+DELR*NRMAX

      FPDATA(1)=DELR
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         FPDATAS(NS)=DELP(NSA)
      ENDDO
      FPDATA(2)=DELTH
      NFPDATA(1)=NRMAX
      NFPDATA(2)=NPMAX
      NFPDATA(3)=NTHMAX

  900 RETURN
  END SUBROUTINE DPLDFP

!****** LOAD Maxwellian VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPLDFM(ID,IERR)

      USE dpcomm
      USE pllocal
      USE plprof,ONLY: pl_prof_old
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ID
      INTEGER,INTENT(OUT):: IERR
      INTEGER NR,NTH,NS,NSA,NP
      REAL(rkind):: PTH0W,RHON,RN0,TPR,TPP,RT0
      REAL(rkind):: TNPR,TNPP,SUM,PML,PPP,PPR,EX,TN00,FACTOR
      
      CALL dpfp_allocate

      IF(NRMAX.EQ.1) THEN
         DELR=0.1D0
      ELSE
         DELR=(RMAX-RMIN)/(NRMAX-1)
      ENDIF
      DO NR=1,NRMAX
         RM(NR)=RMIN+DELR*(NR-1)
      ENDDO

      DELTH=PI/NTHMAX
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

      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         DELP(NSA)=PMAX(NSA)/NPMAX
         DO NP=1,NPMAX
            PM(NP,NSA)=DELP(NSA)*(NP-0.5D0)
            PG(NP,NSA)=DELP(NSA)* NP
         ENDDO

         PN0 = PN(NSA)
         PT0 = (PTPR(NSA)+2*PTPP(NSA))/3.D0
         PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NSA))

         RNFP0(NSA)=PN0
         RTFP0(NSA)=PT0
         AMFP(NSA)=PA(NS)

         IF(ID.EQ.0) THEN
            PTH0W=0.D0
         ELSE
            PTH0W=(PTH0/(AMP*PA(NSA)*VC))**2
         ENDIF

         DO NR=1,NRMAX
            RHON=RM(NR)
            CALL PL_PROF_OLD(RHON)

            RN0 = RN(NSA)
            TPR = RTPR(NSA)
            TPP = RTPP(NSA)
            RT0 = (TPR+2.D0*TPP)/3.D0

            IF(ID.EQ.0) THEN
               TNPR=TPR/PT0
               TNPP=TPP/PT0
               SUM=0.D0
               DO NP=1,NPMAX
                  PML=PM(NP,NSA)
                  DO NTH=1,NTHMAX
                     PPP=PML*TSNM(NTH)
                     PPR=PML*TCSM(NTH)
                     EX=0.5D0*(PPR**2/TNPR+PPP**2/TNPP)
                     FM(NP,NTH) = EXP(-EX)
                     SUM=SUM+FM(NP,NTH)*PM(NP,NSA)*PM(NP,NSA)*TSNM(NTH)
                  ENDDO
               ENDDO
               SUM=SUM*2.D0*PI*DELP(NSA)*DELTH
            ELSE
               TN00=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
               TNPR=TPR*1.D3*AEE/(AMP*PA(NS)*VC**2)
               TNPP=TPP*1.D3*AEE/(AMP*PA(NS)*VC**2)
               SUM=0.D0
               DO NP=1,NPMAX
                  PML=PM(NP,NSA)
                  DO NTH=1,NTHMAX
                     PPP=PML*TSNM(NTH)
                     PPR=PML*TCSM(NTH)
                     EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP) &
                       -SQRT(1.D0/TN00**2)
                     FM(NP,NTH) = EXP(-EX)
                     SUM=SUM+FM(NP,NTH)*PM(NP,NSA)*PM(NP,NSA)*TSNM(NTH)
                  ENDDO
               ENDDO
               SUM=SUM*2.D0*PI*DELP(NSA)*DELTH
            ENDIF

            FACTOR=RN0/(PN0*SUM)
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNS(NTH,NP,NR,NSA) = FACTOR*FM(NP,NTH)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      FPDATA(1)=DELR
      FPDATA(2)=DELTH
      DO NSA=1,NSAMAX
         FPDATA(NSA)=DELP(NSA)
      END DO
      NFPDATA(1)=NRMAX
      NFPDATA(2)=NPMAX
      NFPDATA(3)=NTHMAX

      IERR=0
      RETURN
  END SUBROUTINE DPLDFM

!     ****** SET LOCAL VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPFPFL(NS1,IERR)

      USE dpcomm
      USE pllocal
      USE plprof,ONLY: pl_bminmax
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NS1
      INTEGER,INTENT(OUT):: IERR
      REAL(rkind),DimensION(:),ALLOCATABLE:: THT,FPT,FPTX,FPR,FPRX
      REAL(rkind),DimensION(:,:),ALLOCATABLE:: U2
      INTEGER:: NSA,NSA1,NS,NP,NTH,NR,ID,NTH2
      REAL(rkind):: PTH0W,RHON,BMINL,BMAXL,PSIS,TH,XX,X,Y

      ALLOCATE(THT(NTHMAX+2),FPT(NTHMAX+2),FPTX(NTHMAX+2))
      ALLOCATE(FPR(NRMAX+2),FPRX(NRMAX+2))
      ALLOCATE(U2(4,NTHMAX+2))

      DO NSA1=1,NSAMAX
         IF(NS1.EQ.NS_NSA(NSA1)) THEN
            NSA=NSA1
            GO TO 1
         ENDIF
      ENDDO
      WRITE(6,*) 'XX DPFPFL: unknown NS1: NS1=',NS1
      IERR=6001
      RETURN

    1 CONTINUE

      PN0   =RNFP0(NSA)
      PT0   =RTFP0(NSA)
      PTH0  =SQRT(PT0*1.D3*AEE*AMFP(NSA))
      DELR  =FPDATA(1)
      DO NS=1,NSMAX
         DELP(NS) =FPDATAS(NS)
      ENDDO
      DELTH =FPDATA(2)
      NRMAX =NFPDATA(1)
      NPMAX =NFPDATA(2)
      NTHMAX=NFPDATA(3)

      PTH0W=(PTH0/(AMFP(NSA)*VC))**2

      DO NP=1,NPMAX
         PM(NP,NS)=DELP(NS)*(NP-0.5D0)
         PG(NP,NS)=DELP(NS)* NP
         RGMM(NP)=SQRT(1.D0+PTH0W*PM(NP,NS)**2)
         RGMG(NP)=SQRT(1.D0+PTH0W*PG(NP,NS)**2)
      ENDDO

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
            IF(IERR.NE.0) RETURN
         ENDDO
      ENDDO

      RHON=RHON_LOC
      CALL PL_BMINMAX(RHON,BMINL,BMAXL)
      PSIS=BABS/BMINL

      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            CALL SPL1DF(RHON,Y,RFP,URFP(1,1,NP,NTH),NRMAX+2,IERR)
            IF(IERR.NE.0) RETURN
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
         IF(IERR.NE.0) RETURN
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
            IF(IERR.NE.0) RETURN
            FM(NP,NTH2)=Y
         ENDDO
      ENDDO
      IERR=0
      RETURN
  END SUBROUTINE DPFPFL

!****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPFMFL(NS,ID)

      USE dpcomm
      USE pllocal
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NS,ID
      REAL(RKIND):: RN0,TPR,TPP,RT0,PTH0W,TNPR,TNPP,SUM,PPP,PPR,EX
      REAL(RKIND):: TN00,PML,FACTOR
      INTEGER:: NP,NTH

      RHON_MIN=0.D0
      RHON_MAX=1.D0

      PN0 = PN(NS)
      PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
      PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))

      RN0 = RN(NS)
      TPR = RTPR(NS)
      TPP = RTPP(NS)
      RT0 = (TPR+2.D0*TPP)/3.D0

      IF(ID.EQ.0) THEN
         PTH0W=0.D0
      ELSE
         PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      ENDIF

      DELP(NS)=PMAX(NS)/NPMAX
      DELTH=PI/NTHMAX

      DO NP=1,NPMAX
         PM(NP,NS)=DELP(NS)*(NP-0.5D0)
         PG(NP,NS)=DELP(NS)* NP
         RGMM(NP)=SQRT(1.D0+PTH0W*PM(NP,NS)**2)
         RGMG(NP)=SQRT(1.D0+PTH0W*PG(NP,NS)**2)
      ENDDO

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

      IF(ID.EQ.0) THEN
         TNPR=TPR/RT0
         TNPP=TPP/RT0
         SUM=0.D0
         DO NP=1,NPMAX
            PML=PM(NP,NS)
            DO NTH=1,NTHMAX
               PPP=PML*TSNM(NTH)
               PPR=PML*TCSM(NTH)
               EX=(PPR**2/TNPR+PPP**2/TNPP)*0.5D0
               FM(NP,NTH)=EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP,NS)*PM(NP,NS)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP(NS)*DELTH
      ELSE
         TN00=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
         TNPR=TPR*1.D3*AEE/(AMP*PA(NS)*VC**2)
         TNPP=TPP*1.D3*AEE/(AMP*PA(NS)*VC**2)
         SUM=0.D0
         DO NP=1,NPMAX
            PML=PM(NP,NS)
            DO NTH=1,NTHMAX
               PPP=PML*TSNM(NTH)
               PPR=PML*TCSM(NTH)
               EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP) &
                 -SQRT(1.D0/TN00**2)
               FM(NP,NTH) = EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP,NS)*PM(NP,NS)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP(NS)*DELTH
      ENDIF

      FACTOR=1.D0/SUM
      DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            FM(NP,NTH) = FACTOR*FM(NP,NTH)
         ENDDO
      ENDDO
      RETURN
  END SUBROUTINE DPFMFL

END MODULE dpfpin
