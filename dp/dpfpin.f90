MODULE dpfpin

CONTAINS

!****** LOAD VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPLDFP(IERR)

    USE dpcomm
    USE plcomm_type
    USE libfio
    USE commpi
    USE libmpi
    USE libfio
    IMPLICIT NONE
      INTEGER,INTENT(OUT):: IERR
      CHARACTER(LEN=80),SAVE::  KNAMFP_SAVE=' '
      INTEGER:: NSA,NTH,NP,NR,NS,NSBMAX_DP
      REAL(rkind):: PT0,PTH0,PTH0W,RMIN,RMAX
      
      IERR=0
      IF(KNAMFP.EQ.KNAMFP_SAVE) RETURN

      IF(nrank.EQ.0) THEN

         CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
            RETURN
         ENDIF
         REWIND(21)

         READ(21) NRMAX_DP,NPMAX_DP,NTHMAX_DP,NSAMAX_DP,NSBMAX_DP,NSMAX

         CALL DPFP_ALLOCATE

         READ(21) DELR,DELTH,RMIN,RMAX
         DO NSA=1,NSAMAX_DP
            RHON_MIN(NSA)=RMIN
            RHON_MAX(NSA)=RMAX
            READ(21) NS_NSA_DP(NSA)
            READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         END DO

         NSA_NS_DP(1:NSMAX)=0
         DO NSA=1,NSAMAX_DP
            NS=NS_NSA_DP(NSA)
            NSA_NS_DP(NS)=NSA
         END DO

         DO NSA=1,NSAMAX_DP
            READ(21) DELP(NSA),PMAX_dp(NSA),EMAX_dp(NSA)
         END DO

         DO NSA=1, NSAMAX_DP
            DO NR=1, NRMAX_DP
               DO NP=1, NPMAX_DP
                  READ(21) (FNS(NTH,NP,NR,NSA), NTH=1, NTHMAX_DP)
               END DO
            END DO
         END DO

         CLOSE(21)
      END IF
      
      CALL mtx_broadcast1_integer(NRMAX_DP)
      CALL mtx_broadcast1_integer(NPMAX_DP)
      CALL mtx_broadcast1_integer(NTHMAX_DP)
      CALL mtx_broadcast1_integer(NSAMAX_DP)
      CALL mtx_broadcast1_integer(NSBMAX_DP)
      CALL mtx_broadcast1_integer(NSMAX)
      CALL mtx_broadcast1_real8(DELR)
      CALL mtx_broadcast1_real8(DELTH)
      CALL mtx_broadcast1_real8(RMIN)
      CALL mtx_broadcast1_real8(RMAX)
      CALL mtx_broadcast_integer(NS_NSA_DP,NSAMAX_DP)
      CALL mtx_broadcast_integer(NSA_NS_DP,NSMAX)
      CALL mtx_broadcast_real8(AEFP,NSAMAX_DP)
      CALL mtx_broadcast_real8(AMFP,NSAMAX_DP)
      CALL mtx_broadcast_real8(RNFP0,NSAMAX_DP)
      CALL mtx_broadcast_real8(RTFP0,NSAMAX_DP)
      CALL mtx_broadcast_real8(DELP,NSAMAX_DP)
      CALL mtx_broadcast_real8(pmax_dp,NSAMAX_DP)
      CALL mtx_broadcast_real8(emax_dp,NSAMAX_DP)

      DO NSA=1,NSAMAX_DP
         CALL mtx_broadcast_real8(FNS(1:NTHMAX_DP,1:NPMAX_DP,1:NRMAX_DP,NSA), &
                                  NTHMAX_DP*NPMAX_DP*NRMAX_DP)
      ENDDO

      IF(nrank.eq.0) THEN
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'NRMAX_DP,NPMAX_DP,NTHMAX_DP,NSAMAX_DP =', &
           NRMAX_DP,NPMAX_DP,NTHMAX_DP,NSAMAX_DP
      WRITE(6,'(A,1P4E12.4)') 'DELR/TH,RMIN/MAX =', &
                               DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX_DP
         WRITE(6,*) 'NSA,NS =',NSA,NS_NSA_DP(NSA)
         WRITE(6,'(A,1P5E12.4)') 'AE,AM,RN0,RT0,DELP=', &
               AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA),DELP(NSA)
      ENDDO
      ENDIF

      DELTH=PI/NTHMAX_DP

      DO NSA=1,NSAMAX_DP
         PT0=RTFP0(NSA)
         PTH0=SQRT(PT0*1.D3*AEE*AMFP(NSA))
         PTH0W=(PTH0/(AMFP(NSA)*VC))**2
      DO NP=1,NPMAX_DP
         PM(NP,NSA)=DELP(NSA)*(NP-0.5D0)
         PG(NP,NSA)=DELP(NSA)* NP
         RGMM(NP,NSA)=SQRT(1.D0+PTH0W*PM(NP,NSA)**2)
         RGMG(NP,NSA)=SQRT(1.D0+PTH0W*PG(NP,NSA)**2)
      ENDDO
      ENDDO

      DO NTH=1,NTHMAX_DP
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO

      IF(NRMAX_DP.EQ.1) THEN
         DELR=1.D0
      ENDIF

      RFP(1)=RMIN
      DO NR=1,NRMAX_DP
         RFP(NR+1)=RMIN+DELR*(NR-0.5D0)
      ENDDO
      RFP(NRMAX_DP+2)=RMIN+DELR*NRMAX_DP

      FPDATA(1)=DELR
      FPDATA(2)=DELTH
      DO NSA=1,NSAMAX_DP
         NS=NS_NSA_DP(NSA)
         FPDATAS(NS)=DELP(NSA)
      ENDDO
      NFPDATA(1)=NRMAX_DP
      NFPDATA(2)=NPMAX_DP
      NFPDATA(3)=NTHMAX_DP

      RETURN
  END SUBROUTINE DPLDFP

!****** LOAD Maxwellian VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPLDFM(NSA,ID,IERR)

    USE dpcomm
    USE plprofw
    USE libbes
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NSA,ID  ! ID=0 : non-relativistic, ID=1: relativistic
    INTEGER,INTENT(OUT):: IERR
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    INTEGER NR,NTH,NS,NP,NSA1
    REAL(rkind):: PTH0W,RHON,RN0,TPR,TPP,RT0,PN0,PT0,PTH0
    REAL(rkind):: TNPR,TNPP,SUM,PML,PPP,PPR,EX,TN00,FACTOR,TNL

    NS=NS_NSA_DP(NSA)

    CALL dpfp_allocate

    IF(NRMAX_DP.EQ.1) THEN
       DELR=0.1D0
    ELSE
       DELR=(RHON_MAX(NSA)-RHON_MIN(NSA))/(NRMAX_DP-1)
    ENDIF
    DO NR=1,NRMAX_DP
       RM(NR)=RHON_MIN(NSA)+DELR*(NR-1)
    ENDDO

    DELTH=PI/NTHMAX_DP
    DO NTH=1,NTHMAX_DP
       THM(NTH)=DELTH*(NTH-0.5D0)
       THG(NTH)=DELTH* NTH
       TSNM(NTH) = SIN(THM(NTH))
       TSNG(NTH) = SIN(THG(NTH))
       TCSM(NTH) = COS(THM(NTH))
       TCSG(NTH) = COS(THG(NTH))
       TTNM(NTH) = TAN(THM(NTH))
       TTNG(NTH) = TAN(THG(NTH))
    ENDDO

    DELP(NSA)=PMAX_dp(NSA)/NPMAX_DP
         
    DO NP=1,NPMAX_DP
       PM(NP,NSA)=DELP(NS)*(NP-0.5D0)
       PG(NP,NSA)=DELP(NS)* NP
    ENDDO

    PN0 = PN(NS)
    PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
    PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))

    RNFP0(NSA)=PN0
    RTFP0(NSA)=PT0
    AMFP(NSA)=AMP*PA(NS)

    IF(ID.EQ.0) THEN
       PTH0W=0.D0
    ELSE
       PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
    ENDIF

    DO NR=1,NRMAX_DP
       IF(NRMAX_DP.EQ.1) THEN
          RHON=0.D0
       ELSE
          RHON=RM(NR)
       END IF
       CALL PL_PROFW(RHON,plfw)

       RN0 = plfw(NSA)%RN
       TPR = plfw(NSA)%RTPR
       TPP = plfw(NSA)%RTPP
       RT0 = (TPR+2.D0*TPP)/3.D0

       IF(ID.EQ.0) THEN
          TNPR=TPR/PT0
          TNPP=TPP/PT0
          SUM=0.D0
          DO NP=1,NPMAX_DP
             PML=PM(NP,NSA)
             DO NTH=1,NTHMAX_DP
                PPP=PML*TSNM(NTH)
                PPR=PML*TCSM(NTH)
                EX=0.5D0*(PPR**2/TNPR+PPP**2/TNPP)
                FM(NP,NTH) = EXP(-EX)
                SUM=SUM+FM(NP,NTH)*PM(NP,NSA)*PM(NP,NSA)*TSNM(NTH)
             ENDDO
          ENDDO
          SUM=SUM*2.D0*PI*DELP(NS)*DELTH
          SUM=1.D0/(SQRT(2.D0*PI)**3*SQRT(TNPR)*(TNPP))
       ELSE
          TN00=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
          TNPR=TPR*1.D3*AEE/(AMP*PA(NS)*VC**2)
          TNPP=TPP*1.D3*AEE/(AMP*PA(NS)*VC**2)
          SUM=0.D0
          DO NP=1,NPMAX_DP
             PML=PM(NP,NSA)
             DO NTH=1,NTHMAX_DP
                PPP=PML*TSNM(NTH)
                PPR=PML*TCSM(NTH)
                EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP) &
                  -SQRT(1.D0/TN00**2)
                FM(NP,NTH) = EXP(-EX)
                SUM=SUM+FM(NP,NTH)*PM(NP,NSA)*PM(NP,NSA)*TSNM(NTH)
             ENDDO
          ENDDO
          SUM=SUM*2.D0*PI*DELP(NS)*DELTH
          TNL=(TNPR+2.D0*TNPP)/3.D0
          SUM=SUM*SQRT(TNL)/(4.D0*PI*BESEKNX(2,1.D0/TNL) &
                           *SQRT(TPR/PT0)*(TPP/PT0))
       ENDIF

       FACTOR=RN0/(PN0*SUM)
       DO NP=1,NPMAX_DP
          DO NTH=1,NTHMAX_DP
             FNS(NTH,NP,NR,NSA) = FACTOR*FM(NP,NTH)
          ENDDO
       ENDDO
    ENDDO

    FPDATA(1)=DELR
    FPDATA(2)=DELTH
    DO NSA1=1,NSAMAX_DP
       FPDATAS(NSA1)=DELP(NSA1)
    END DO
    NFPDATA(1)=NRMAX_DP
    NFPDATA(2)=NPMAX_DP
    NFPDATA(3)=NTHMAX_DP
 
    IERR=0
    RETURN
  END SUBROUTINE DPLDFM

!     ****** SET LOCAL VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPFPFL(NSA,mag,ierr)

    USE dpcomm
    USE plcomm_type
      USE plprof
      USE libspl1d
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSA
      TYPE(pl_mag_type),INTENT(IN):: mag
      INTEGER,INTENT(OUT):: IERR
      REAL(rkind),DimensION(:),ALLOCATABLE:: THT,FPT,FPTX,FPR,FPRX
      REAL(rkind),DimensION(:,:),ALLOCATABLE:: U2
      INTEGER:: NS,NP,NTH,NR,ID,NTH2,NSA1
      REAL(rkind):: PN0,PT0,PTH0,PTH0W,RHON,BMINL,BMAXL,PSIS,TH,XX,X,Y

      ALLOCATE(THT(NTHMAX_DP+2),FPT(NTHMAX_DP+2),FPTX(NTHMAX_DP+2))
      ALLOCATE(FPR(NRMAX_DP+2),FPRX(NRMAX_DP+2))
      ALLOCATE(U2(4,NTHMAX_DP+2))

      NS=NS_NSA_DP(NSA)

      PN0   =RNFP0(NSA)
      PT0   =RTFP0(NSA)
      PTH0  =SQRT(PT0*1.D3*AEE*AMFP(NSA))
      DELR  =FPDATA(1)
      DELTH =FPDATA(2)
      DELP(NS) =FPDATAS(NSA)
      NRMAX_DP =NFPDATA(1)
      NPMAX_DP =NFPDATA(2)
      NTHMAX_DP=NFPDATA(3)

      PTH0W=(PTH0/(AMFP(NSA)*VC))**2

      DO NP=1,NPMAX_DP
         PM(NP,NS)=DELP(NSA)*(NP-0.5D0)
         PG(NP,NS)=DELP(NSA)* NP
         RGMM(NP,NSA)=SQRT(1.D0+PTH0W*PM(NP,NSA)**2)
         RGMG(NP,NSA)=SQRT(1.D0+PTH0W*PG(NP,NSA)**2)
      ENDDO

      DO NTH=1,NTHMAX_DP
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO

      DO NP=1,NPMAX_DP
        DO NTH=1,NTHMAX_DP
            DO NR=1,NRMAX_DP
               FP(NTH,NP,NR)=FNS(NTH,NP,NR,NSA)
            ENDDO
            IF(RHON_MIN(NSA).EQ.0.D0)THEN
               FPR(1)=(9.D0*FP(NTH,NP,1)-FP(NTH,NP,2))/8.D0
               FPRX(1)=0.D0
               ID=1
            ELSE
               FPR(1)=(4.D0*FP(NTH,NP,1)-FP(NTH,NP,2))/3.D0
               ID=0
            ENDIF
            DO NR=1,NRMAX_DP
               FPR(NR+1)=FP(NTH,NP,NR)
            ENDDO
            FPR(NRMAX_DP+2) &
                 =(4.D0*FP(NTH,NP,NRMAX_DP+1)-FP(NTH,NP,NRMAX_DP))/3.D0
            CALL SPL1D(RFP,FPR,FPRX,URFP(1,1,NP,NTH),NRMAX_DP+2,ID,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX DPFPFL: SPL1D: RFP-FPR: IERR=',IERR
               STOP
            END IF
         ENDDO
      ENDDO

      RHON=mag%rhon
      CALL PL_BMINMAX(RHON,BMINL,BMAXL)
      PSIS=mag%BABS/BMINL

      DO NP=1,NPMAX_DP
         DO NTH=1,NTHMAX_DP
            CALL SPL1DF(RHON,Y,RFP,URFP(1,1,NP,NTH),NRMAX_DP+2,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX DPFPFL: SPL1DF: RHON-RFP: IERR=',IERR
               STOP
            END IF
            FPT(NTH+1)=Y
         ENDDO
         THT(1)=0.D0
         DO NTH=1,NTHMAX_DP
            THT(NTH+1)=DELTH*(NTH-0.5D0)
         ENDDO
         THT(NTHMAX_DP+2)=PI
         FPT(1)=(9.D0*FPT(2)-FPT(3))/8.D0
         FPT(NTHMAX_DP+2)=(9.D0*FPT(NTHMAX_DP+1)-FPT(NTHMAX_DP))/8.D0         
         FPTX(1)=0.D0
         FPTX(NTHMAX_DP+2)=0.D0
         CALL SPL1D(THT,FPT,FPTX,U2,NTHMAX_DP+2,3,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX DPFPFL: SPL1D: THT-FPT: IERR=',IERR
            STOP
         END IF
         DO NTH2=1,NTHMAX_DP
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
            CALL SPL1DF(X,Y,THT,U2,NTHMAX_DP+2,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX DPFPFL: SPL1DF: THT-FPT: IERR=',IERR
               STOP
            END IF
            IF(IERR.NE.0) RETURN
            FM(NP,NTH2)=Y
         ENDDO
      ENDDO
      IERR=0
      RETURN
  END SUBROUTINE DPFPFL

!****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA ******

  SUBROUTINE DPFMFL(NS,plfw,ID)

      USE dpcomm
      USE plprofw
      USE libbes
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NS,ID ! ID=0 : non-relativistic, ID=1: relativistic
      TYPE(pl_prfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      REAL(RKIND):: RN0,TPR,TPP,RT0,PN0,PT0,PTH0,PTH0W,TNPR,TNPP,SUM
      REAL(RKIND):: PPP,PPR,EX,EXX,TN0,TNL,PML,FACTOR
      INTEGER:: NP,NTH,NSA

      CALL dpfp_allocate

      RHON_MIN=0.D0
      RHON_MAX=1.D0

      PN0 = PN(NS)
      PT0 = (PTPR(NS)+2*PTPP(NS))/3.D0
      PTH0 = SQRT(PT0*1.D3*AEE*AMP*PA(NS))

      NSA=NSA_NS_DP(NS)
      RNFP0(NSA)=PN0
      RTFP0(NSA)=PT0
      AMFP(NSA)=AMP*PA(NS)

      RN0 = plfw(NS)%RN
      TPR = plfw(NS)%RTPR
      TPP = plfw(NS)%RTPP
      RT0 = (TPR+2.D0*TPP)/3.D0

      IF(ID.EQ.0) THEN
         PTH0W=0.D0
      ELSE
         PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      ENDIF

      DELP(NSA)=PMAX_dp(NSA)/NPMAX_DP
      DELTH=PI/NTHMAX_DP

      DO NP=1,NPMAX_DP
         PM(NP,NSA)=DELP(NSA)*(NP-0.5D0)
         PG(NP,NSA)=DELP(NSA)* NP
         RGMM(NP,NSA)=SQRT(1.D0+PTH0W*PM(NP,NSA)**2)
         RGMG(NP,NSA)=SQRT(1.D0+PTH0W*PG(NP,NSA)**2)
      ENDDO

      DO NTH=1,NTHMAX_DP
         THM(NTH)=DELTH*(NTH-0.5D0)
         THG(NTH)=DELTH* NTH
         TSNM(NTH) = SIN(THM(NTH))
         TSNG(NTH) = SIN(THG(NTH))
         TCSM(NTH) = COS(THM(NTH))
         TCSG(NTH) = COS(THG(NTH))
         TTNM(NTH) = TAN(THM(NTH))
         TTNG(NTH) = TAN(THG(NTH))
      ENDDO

      TNPR=TPR/PT0
      TNPP=TPP/PT0
      IF(ID.EQ.0) THEN
         FACTOR=1.D0/(SQRT(2.D0*PI)**3*SQRT(TNPR)*(TNPP))
         SUM=0.D0
         DO NP=1,NPMAX_DP
            PML=PM(NP,NSA)
            DO NTH=1,NTHMAX_DP
               PPR=PML*TCSM(NTH)
               PPP=PML*TSNM(NTH)
               EX=0.5D0*(PPR**2/TNPR+PPP**2/TNPP)
               FM(NP,NTH)=FACTOR*EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP,NSA)*PM(NP,NSA)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP(NSA)*DELTH
      ELSE
         TN0=PT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
         TNL=RT0*1.D3*AEE/(AMP*PA(NS)*VC**2)
         FACTOR=SQRT(TNL)/(4.D0*PI*BESEKNX(2,1.D0/TNL) &
               *SQRT(TNPR)*(TNPP))
         SUM=0.D0
         DO NP=1,NPMAX_DP
            PML=PM(NP,NSA)
            DO NTH=1,NTHMAX_DP
               PPP=PML*TSNM(NTH)
               PPR=PML*TCSM(NTH)
               EXX=PPR**2/TNPR+PPP**2/TNPP
               EX=(SQRT(1.D0+TN0*EXX)-1.D0)/TN0
               FM(NP,NTH) = FACTOR*EXP(-EX)
               SUM=SUM+FM(NP,NTH)*PM(NP,NSA)*PM(NP,NSA)*TSNM(NTH)
            ENDDO
         ENDDO
         SUM=SUM*2.D0*PI*DELP(NSA)*DELTH
      ENDIF

      FACTOR=RN0/PN0
      DO NP=1,NPMAX_DP
         DO NTH=1,NTHMAX_DP
            FM(NP,NTH) = FACTOR*FM(NP,NTH)
         ENDDO
      ENDDO

      RETURN
  END SUBROUTINE DPFMFL

END MODULE dpfpin
