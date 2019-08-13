C     $Id$
C
C****** SET MAXWELLIAN VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPFMFL(NS,ID)
C
      USE plcomm
      USE pllocal
      INCLUDE 'dpcomm.inc'
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
      DELP(NS)=PMAX(NS)/NPMAX
      DELTH=PI/NTHMAX
C
      DO NP=1,NPMAX
         PM(NP,NS)=DELP(NS)*(NP-0.5D0)
         PG(NP,NS)=DELP(NS)* NP
         RGMM(NP)=SQRT(1.D0+PTH0W*PM(NP,NS)**2)
         RGMG(NP)=SQRT(1.D0+PTH0W*PG(NP,NS)**2)
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
            PML=PM(NP,NS)
            DO NTH=1,NTHMAX
               PPP=PML*TSNM(NTH)
               PPR=PML*TCSM(NTH)
               EX=-(PPR**2/TNPR+PPP**2/TNPP)*0.5D0 !+ *0.5D0 2015/07/28
!               FM(NP,NTH) = EXP(-EX) !original
               FM(NP,NTH)=EXP(EX)    !2015/07/28
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
               EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP)
     &           -SQRT(1.D0/TN00**2)
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
      END                     
C
C     ****** SET LOCAL VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPFPFL(NS)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: pl_bminmax
      INCLUDE 'dpcomm.inc'
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
      DO NS=1,NSMAX
         DELP(NS) =FPDATAS(NS)
      ENDDO
      DELTH =FPDATA(2)
      NRMAX =NFPDAT(1)
      NPMAX =NFPDAT(2)
      NTHMAX=NFPDAT(3)
C
      PTH0W=(PTH0/(AMFP(NSA)*VC))**2
C
      DO NP=1,NPMAX
         PM(NP,NS)=DELP(NS)*(NP-0.5D0)
         PG(NP,NS)=DELP(NS)* NP
         RGMM(NP)=SQRT(1.D0+PTH0W*PM(NP,NS)**2)
         RGMG(NP)=SQRT(1.D0+PTH0W*PG(NP,NS)**2)
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
      CALL PL_BMINMAX(RHON,BMINL,BMAXL)
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
      USE plcomm
      USE commpi
      USE libmpi
      USE libfio
      INCLUDE 'dpcomm.inc'
      CHARACTER(LEN=80),SAVE::  KNAMFP_SAVE=' '
      REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: FNSDP
C     
      IF(KNAMFP.EQ.KNAMFP_SAVE) RETURN
C
      IF(nrank.eq.0) THEN
      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
       IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
       ENDIF
      REWIND(21)

      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) DELP(NSA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX),
     &                NP=1,NPMAX),NR=1,NRMAX)
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

      ALLOCATE(FNSDP(NTHMAX,NPMAX,NRMAX,NSAMAX))
      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNSDP(NTH,NP,NR,NSA)=
     &                 FNS(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      DO NSA=1,NSAMAX
      CALL mtx_broadcast_real8(FNSDP(1:NTHMAX,1:NPMAX,1:NRMAX,NSA),
     &                         NTHMAX*NPMAX*NRMAX)
      ENDDO

!      DO NSA=1,NSAMAX
!         DO NR=1,NRMAX
!            DO NP=1,NPMAX
!               DO NTH=1,NTHMAX
!                  FNS(NTH,NP,NR,NSA)=
!     &                 FNSDP(NTH,NP,NR,NSA)
!               ENDDO
!            ENDDO
!         ENDDO
!      ENDDO

C     
      IF(nrank.eq.0) THEN
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'NRMAX,NPMAX,NTHMAX,NSAMAX =',NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(6,'(A,1P4E12.4)') 'DELR/TH,RMIN/RMAX =',
     &                         DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         WRITE(6,*) 'NSA,NS =',NSA,NS_NSA(NSA)
         WRITE(6,'(A,1P5E12.4)') 'AE,AM,RN0,RT0,DELP=',
     &        AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA),DELP(NSA)
      ENDDO
      ENDIF
C     
      DELTH=PI/NTHMAX
C
      DO NS=1,NSMAX
      DO NP=1,NPMAX
         PM(NP,NS)=DELP(NS)*(NP-0.5D0)
         PG(NP,NS)=DELP(NS)* NP
      ENDDO
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
!         write(6,'(A,I5,1P3E12.4)') 
!     &        'NTH=',NTH,THM(NTH),TSNM(NTH),TCSM(NTH)
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
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         FPDATAS(NS)=DELP(NSA)
      ENDDO
      FPDATA(2)=DELTH
      NFPDAT(1)=NRMAX
      NFPDAT(2)=NPMAX
      NFPDAT(3)=NTHMAX
C
!!!--fa
      IF(MODEFA.EQ.4) THEN
       NSA=3
       NRMAXFP=NRMAX
       NTHMAXFP=NTHMAX
        PMa0(1)=0.D0
       DO NP=1,NPMAX
        PMa0(NP+1)=DELP(NSA)*(REAL(NP)-0.5D0)
       ENDDO
        PMa0(NPMAX+2)=DELP(NSA)*NPMAX

        RHOa0(1)=RMIN/RMAX
       DO NR=1,NRMAXFP
        RHOa0(NR+1)=RMIN/RMAX+DELR*(REAL(NR)-0.5D0)
       ENDDO
        RHOa0(NRMAXFP+2)=RMAX/RMAX

       DO NTH=1,NTHMAXFP
       DO NR=1,NRMAXFP
       DO NP=1,NPMAX
          fa0(NP+1,NR+1,NTH)=FNSDP(NTH,NP,NR,NSA)
          fa0(NP+1,1,NTH) =FNSDP(NTH,NP,1,NSA)
          fa0(NP+1,NRMAXFP+2,NTH)=FNSDP(NTH,NP,NRMAXFP,NSA)*1.D-1!0.D0
       ENDDO
          fa0(1,NR+1,NTH) =FNSDP(NTH,1,NR,NSA)
          fa0(NPMAX+2,NR+1,NTH)=FNSDP(NTH,NPMAX,NR,NSA)*1.D-1 !2.5D-1 ! f(NP+2)=f(NP+1)*1/4
       ENDDO
          fa0(1,1,NTH)=FNSDP(NTH,1,1,NSA)
          fa0(1,NRMAXFP+2,NTH)=FNSDP(NTH,1,NRMAX,NSA)*1.D-1 !0.D0
          fa0(NPMAX+2,1,NTH)=FNSDP(NTH,NPMAX,1,NSA)*1.D-1 !2.5D-1          
          fa0(NPMAX+2,NRMAXFP+2,NTH)=FNSDP(NTH,NPMAX,NRMAX,NSA)*1.D-2!0.D0
       ENDDO

       DO NTH=1,NTHMAXFP
       DO NR=1,NRMAXFP+2
       DO NP=1,NPMAX+2
          dfpa0(NP,NR,NTH)=0.D0
          dfra0(NP,NR,NTH)=0.D0 !0.D0
       ENDDO
       ENDDO
       ENDDO

       DO NTH=1,NTHMAXFP
        CALL SPL2D(PMa0,RHOa0,fa0(1:NPM,1:NRM,NTH),dfpa0,dfra0,
     &  dfpra0,US(1:4,1:4,1:NPM,1:NRM,NTH),
     &  NPM,NPMAX+2,NRMAX+2,1,1,IERR)  !,1,1,IERR)
       ENDDO   
      ELSE
       RETURN
      ENDIF
      DEALLOCATE(FNSDP)
!!!   --fa
  900 RETURN
      END
C
C****** LOAD Maxwellian VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPLDFM(ID,NCHMAX,NRMAX_1,RMIN_1,RMAX_1)
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: pl_prof_old
      INCLUDE 'dpcomm.inc'
C
      NRMAX=NRMAX_1
      RMIN=RMIN_1
      RMAX=RMAX_1
      RHON_MIN=0.D0
      RHON_MAX=1.D0
C
      IF(NRMAX.EQ.1) THEN
         DELR=0.1D0
      ELSE
         DELR=(RMAX-RMIN)/(NRMAX-1)
      ENDIF
      DO NR=1,NRMAX
         RM(NR)=RMIN+DELR*(NR-1)
      ENDDO

      DELCH=2.D0*PI/NCHMAX
      DO NCH=1,NCHMAX
         CHI(NCH)=DELCH*(NCH-1)
         CCHI(NCH)=COS(CHI(NCH))
         SCHI(NCH)=SIN(CHI(NCH))
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

      DO NS=1,NSAMAX

         DELP(NS)=PMAX(NS)/NPMAX
         DO NP=1,NPMAX
            PM(NP,NS)=DELP(NS)*(NP-0.5D0)
            PG(NP,NS)=DELP(NS)* NP
         ENDDO

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
            CALL PL_PROF_OLD(RHON)

            RN0 = RN(NS)
            TPR = RTPR(NS)
            TPP = RTPP(NS)
C
            IF(ID.EQ.0) THEN
               TNPR=TPR/PT0
               TNPP=TPP/PT0
               SUM=0.D0
               DO NP=1,NPMAX
                  PML=PM(NP,NS)
                  DO NTH=1,NTHMAX
                     PPP=PML*TSNM(NTH)
                     PPR=PML*TCSM(NTH)
                     EX=0.5D0*(PPR**2/TNPR+PPP**2/TNPP)
                     FM(NP,NTH) = EXP(-EX)
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
                     EX=SQRT(1.D0/TN00**2+PPR**2/TNPR+PPP**2/TNPP)
     &                 -SQRT(1.D0/TN00**2)
                     FM(NP,NTH) = EXP(-EX)
                     SUM=SUM+FM(NP,NTH)*PM(NP,NS)*PM(NP,NS)*TSNM(NTH)
                  ENDDO
               ENDDO
               SUM=SUM*2.D0*PI*DELP(NS)*DELTH
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
      
      FPDATA(2)=DELP(NS)
      FPDATA(3)=DELTH
      NFPDAT(1)=NRMAX
      NFPDAT(2)=NPMAX
      NFPDAT(3)=NTHMAX
C
  900 RETURN
      END
