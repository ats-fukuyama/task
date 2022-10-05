!     $Id: fpgout.f90,v 1.13 2013/01/20 23:24:03 fukuyama Exp $
!
! ***********************************************************
!
!                    SELECTION OF GRAPH
!
! ***********************************************************
!
      MODULE FPGOUT

      USE fpcomm
      USE fpcont
      USE fpfout
      USE fpgsub
      USE libgrf
      USE libmpi

      contains
!---------------------------------------
      SUBROUTINE FP_GOUT
!
!     F1    FPGRAPA   FNS 
!     FR1   FPGRAPRA  FNS 
!     F2    FPGRACB   FNS 
!     FX2   FPGRAC    TEMP
!     FY2   FPGRAC    TEMP
!     FS11  FPGRAPA   FS1   ?  FS21
!     FS12  FPGRACXA  FS1   4  FS22
!     FP    FPGRACA   FPP   1  FPP FCP FEP
!     FT    FPGRACA   FTH   2  FTH FCT FET
!     FR    FPGRACA   FRR   0  FRR
!     FCP   FPGRACAB  FCPP2 1
!     FCT   FPGRACAB  FCTH2 2
!     RN    FPGRARA   RNT      RI,RW,RPC,RPW,RPE,RLH,RFW,REC,RT
!     RPC   FPGRARAB  RPCT2    RPC
!     RQ    FPGRAR    RQT      RE
!     TN    FPGRATA   PNT      TI,TW,TPC,TPW,TPE,TLH,TFW,TEC,TT
!     TPC   FPGRATAB  PPCT2    RPC
!     TQ    FPGRAT    PQT      PET
!     DPP   FPGRACA   DPP   1  DPT,DCPP,DCPT,DWPP,DWPT
!     DTT   FPGRACA   DTT   2  DTP,DCTT,DCTP,DWTT,DWTP
!     DRR   FPGRACA   DRR   3  
!     DCPP  FPGRACAB  DCPP2 1  DCPP,DCPT
!     DCTT  FPGRACAB  DCTT2 2  DCTT,DCTP
!     PFP   FPGRACPA  FPP   1  PFPP PFCP PFEP
!     PFT   FPGRACPA  FTH   2  PFTH PFCT PFET
!     PFR   FPGRACPA  FRR   0  PFRR
!     PFCP  FPGRACPAB FCPP2 1
!     PFCT  FPGRACPAB FCTH2 2
!     PDPP  FPGRACPA  DPP   1  PDPT,PDCPP,PDCPT,PDWPP,PDWPT
!     PDTT  FPGRACPA  DTT   2  PDTP,PDCTT,PDCTP,PDWTT,PDWTP
!     PDRR  FPGRACPA  DRR   3  PDRR
!     PDCPP FPGRACPAB DCPP2 1  PDCPP,PDCPT
!     PDCTT FPGRACPAB DCTT2 2  PDCTT,PDCTP
!
!     FPGRARA,FPGRARAB   -> FPGRAR
!     FPGRATA,FPGRATAB   -> FPGRAT
!     FPGRAPA,FPGRAPAB   -> FPGRAP
!     FPGRAPRA,FPGRAPRAB -> FPGRAPR -> FPGR3D
!     FPGR3D
!     FPGRACA,FPGRACAB   -> FPGRAC
!     FPGRACB            -> FPGRAC
!     FPGRACAX           -> FPGRACX
!     FPGRACPA,FPGRACPAB -> FPGRACP

      USE libchar
      IMPLICIT NONE
      integer,SAVE:: NSA=1,NSB=0
      integer:: NR, NP, NTH, NSA1, NS
!
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX)::TEMP
      CHARACTER KID*5,KID1*1,KID2*3
      CHARACTER(LEN=80):: STRING1
!
1     CONTINUE
      IF(nrank.EQ.0) THEN
         WRITE(6,*)'INPUT GRAPH TYPE : NSA,NSB=',NSA,NSB
         WRITE(6,*)'    : F/FX/FR/FS1/FS2 1/2, Nn:NSA=n, Nnm:NSA=n,NSB=m,'
         WRITE(6,*)'    : D/DC/DW PP/PT/TP/TT/RR, F/FC/FE P/T/R'
         WRITE(6,*)'    : PD/PDC/PDW PP/PT/TP/TT/RR, PF/PFC/FE P/T/R'
         WRITE(6,*)'    : R/T N/I/W/PC/PW/PE/T/Q/E,  Gn,  X:exit'
         READ(5,'(A5)',ERR=1,END=9000) KID
      END IF
      CALL mtx_broadcast_character(KID,4)

      KID1=KID(1:1)
      CALL toupper(KID1)
         
      IF(KID1.NE.'P') THEN
         KID2=KID(2:4)
         CALL toupper(KID2(1:1))
         CALL toupper(KID2(2:2))
         CALL toupper(KID2(3:3))
!
      NS=NS_NSA(NSA)
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'1  ') THEN
            IF(nrank.EQ.0) CALL FPGRAPA('F1',FNS,NSA)
         ELSE IF(KID2.EQ.'R1 ') THEN
            IF(nrank.EQ.0) CALL FPGRAPRA('FR1',FNS,NSA)
         ELSE IF(KID2.EQ.'2  ') THEN
            CALL fp_grac(4,FNS,'F2  ')
         ELSE IF(KID2.EQ.'X2  ') THEN
            WRITE(STRING1,'(A,A1,I2,A1)') 'FX2','(',NSA,')'
            IF(nrank.EQ.0) CALL FPGRACD(TRIM(STRING1),FNS,4,NSA)
         ELSE IF(KID2.EQ.'Y2  ') THEN
            WRITE(STRING1,'(A,A1,I2,A1)') 'FY2','(',NSA,')'
            IF(nrank.EQ.0) CALL FPGRACD(TRIM(STRING1),FNS,0,NSA)
         ELSE IF(KID2.EQ.'S11') THEN
            IF(nrank.EQ.0) CALL FPGRAPA('FS11',FS1,NSA)
         ELSE IF(KID2.EQ.'S12') THEN
            IF(nrank.EQ.0) CALL FPGRACXA('FS12',FS1,4,NSA)
         ELSE IF(KID2.EQ.'S21') THEN
            IF(nrank.EQ.0) CALL FPGRAPA('FS21',FS2,NSA)
         ELSE IF(KID2.EQ.'S22') THEN
            IF(nrank.EQ.0) CALL FPGRACXA('FS22',FS2,4,NSA)
         ELSE IF(KID2.EQ.'P  ') THEN
            CALL fp_grac(1,FPP,'FPP ')
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL fp_grac(2,FTH,'FTH ')
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL fp_grac(0,FRR,'FRR ')
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL fp_grac(1,FPP,'FPP ')
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL fp_grac(2,FTH,'FTH ')
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL fp_grac(0,FRR,'FRR ')
         ELSE IF(KID2.EQ.'CP ') THEN
            IF(NSB.EQ.0) THEN
               CALL fp_grac(1,FCPP,'FCP ')
            ELSE
               IF(nrank.EQ.0) CALL FPGRACAB('FCP ',FCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CT ') THEN
            IF(NSB.eq.0) THEN
               CALL fp_grac(2,FCTH,'FCT ')
            ELSE
               IF(nrank.EQ.0) CALL FPGRACAB('FCT ',FCTH2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL fp_grac(1,FEPP,'FEP ')
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL fp_grac(2,FETH,'FET ')
         ELSE
            IF(nrank.EQ.0)WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'R') THEN
         IF(KID2.EQ.'N  ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RN  ',RNT,NSA)
         ELSE IF(KID2.EQ.'I  ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RI  ',RJT,NSA)
         ELSE IF(KID2.EQ.'W  ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RW  ',RWT,NSA)
         ELSE IF(KID2.EQ.'PC ') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRARA('RPC ',RPCT,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRARAB('RPC ',RPCT2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'PW ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RPW ',RPWT,NSA)
         ELSE IF(KID2.EQ.'PE ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RPE ',RPET,NSA)
         ELSE IF(KID2.EQ.'LH ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RLH ',RLHT,NSA)
         ELSE IF(KID2.EQ.'FW ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RFW ',RFWT,NSA)
         ELSE IF(KID2.EQ.'EC ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('REC ',RECT,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RT  ',RTT,NSA)
         ELSE IF(KID2.EQ.'TB ') THEN
            IF(nrank.EQ.0) CALL FPGRARA('RTB ',RTT_BULK,NSA)
         ELSE IF(KID2.EQ.'Q  ') THEN
            IF(nrank.EQ.0) CALL FPGRAR('RQ  ',RQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            IF(nrank.EQ.0) CALL FPGRAR('RE  ',RET)
         ELSE
            IF(nrank.EQ.0)WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'T') THEN
         IF(KID2.EQ.'N  ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TN  ',PNT,NSA)
         ELSE IF(KID2.EQ.'I  ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TI  ',PIT,NSA)
         ELSE IF(KID2.EQ.'W  ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TW  ',PWT,NSA)
         ELSE IF(KID2.EQ.'PC ') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRATA('TPC ',PPCT,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRATAB('TPC ',PPCT2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'PW ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TPW ',PPWT,NSA)
         ELSE IF(KID2.EQ.'PE ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TPE ',PPET,NSA)
         ELSE IF(KID2.EQ.'LH ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TLH ',PLHT,NSA)
         ELSE IF(KID2.EQ.'FW ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TFW ',PFWT,NSA)
         ELSE IF(KID2.EQ.'EC ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TEC ',PECT,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TT  ',PTT,NSA)
         ELSE IF(KID2.EQ.'TB ') THEN
            IF(nrank.EQ.0) CALL FPGRATA('TTB ',PTT_BULK,NSA)
         ELSE IF(KID2.EQ.'Q  ') THEN
            IF(nrank.EQ.0) CALL FPGRAT('TQ  ',PQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            IF(nrank.EQ.0) CALL FPGRAT('TE  ',PET)
         ELSE IF(KID2.EQ.'QE ') THEN
            IF(nrank.EQ.0) CALL FPGRAT('TQE ',Q_ENG)
         ELSE
            IF(nrank.EQ.0)WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL fp_grac(1,DPP,'DPP ')
!            IF(nrank.EQ.0) CALL FPGRACA('DPP ',DPP ,1,NSA)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL fp_grac(1,DPT,'DPT ')
!            IF(nrank.EQ.0) CALL FPGRACA('DPT ',DPT ,1,NSA)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL fp_grac(2,DTP,'DTP ')
!            IF(nrank.EQ.0) CALL FPGRACA('DTP ',DTP ,2,NSA)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL fp_grac(2,DTT,'DTT ')
!            IF(nrank.EQ.0) CALL FPGRACA('DTT ',DTT ,2,NSA)
         ELSE IF(KID2.EQ.'CPP') THEN
            IF(NSB.EQ.0) THEN
               CALL fp_grac(1,DCPP,'DCPP')
!               IF(nrank.EQ.0) CALL FPGRACA('DCPP',DCPP,1,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACAB('DCPP',DCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CPT') THEN
            IF(NSB.EQ.0) THEN
               CALL fp_grac(1,DCPT,'DCPT')
!              IF(nrank.EQ.0) CALL FPGRACA('DCPT',DCPT,1,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACAB('DCPT',DCPT2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTP') THEN
            IF(NSB.EQ.0) THEN
               CALL fp_grac(2,DCTP,'DCTP')
!              IF(nrank.EQ.0) CALL FPGRACA('DCTP',DCTP,2,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACAB('DCTP',DCTP2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTT') THEN
            IF(NSB.EQ.0) THEN
               CALL fp_grac(2,DCTT,'DCTT')
!               IF(nrank.EQ.0) CALL FPGRACA('DCTT',DCTT,2,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACAB('DCTT',DCTT2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL fp_grac(1,DWPP,'DWPP')
!            IF(nrank.EQ.0) CALL FPGRACA('DWPP',DWPP,1,NSA)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL fp_grac(1,DWPT,'DWPT')
!            IF(nrank.EQ.0) CALL FPGRACA('DWPT',DWPT,1,NSA)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL fp_grac(2,DWTP,'DWTP')
!            IF(nrank.EQ.0) CALL FPGRACA('DWTP',DWTP,2,NSA)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL fp_grac(2,DWTT,'DWTT')
!            IF(nrank.EQ.0) CALL FPGRACA('DWTT',DWTT,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL fp_grac(3,DRR ,'DRR ')
!            IF(nrank.EQ.0) CALL FPGRACA('DRR ',DRR ,3,NSA)
         ELSE
            IF(nrank.EQ.0) WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'N') THEN
         IF(nrank.EQ.0) THEN
            READ(KID2,*,ERR=1000,END=1000) NSA1
            IF(NSA1.GT.0) THEN
               IF(NSA1.GE.10) THEN
                  NSA=NSA1/10
                  NSB=MOD(NSA1,10)
               ELSE
                  NSA=NSA1
                  NSB=0
               ENDIF
               WRITE(6,'(A,2I3)') '# NSA and NSB are changed to',NSA,NSB
            ELSE
               WRITE(6,'(A)') 'XX Negative number is not allowed.'
            ENDIF
1000        CONTINUE
         END IF
         CALL mtx_broadcast1_integer(NSA)
         CALL mtx_broadcast1_integer(NSB)
      ELSE IF (KID1.EQ.'W') THEN
         IF(KID2.EQ.'P  ') THEN
            CALL fp_grac(1,WEIGHP,'WP  ')
!            IF(nrank.EQ.0) CALL FPGRACA('WP  ',WEIGHP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL fp_grac(2,WEIGHT,'WT  ')
!            IF(nrank.EQ.0) CALL FPGRACA('WT  ',WEIGHT,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL fp_grac(3,WEIGHR,'WR  ')
!            IF(nrank.EQ.0) CALL FPGRACA('WR  ',WEIGHR,3,NSA)
         ENDIF
      ELSE IF (KID1.EQ.'G') THEN
         IF(nrank.EQ.0) THEN
            READ(KID2,*,ERR=2000,END=2000) NGRAPH
            WRITE(6,'(A,I3)') '# NGRAPH is changed to',NGRAPH
2000        CONTINUE
         END IF
         CALL mtx_broadcast1_integer(NGRAPH)
      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000
      ELSE IF (KID1.EQ.'Q') THEN
         GO TO 9000
      ELSE
         IF(nrank.EQ.0) WRITE(6,*) 'XX UNKNOWN KID1'
      END IF
!
      ELSE
         KID1=KID(2:2)
         CALL toupper(KID1)
         KID2=KID(3:5)
         CALL toupper(KID2(1:1))
         CALL toupper(KID2(2:2))
         CALL toupper(KID2(3:3))
!
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'P  ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FP  ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FT  ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FR  ',FRR,3,NSA)
         ELSE IF(KID2.EQ.'PP ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FPP ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'TH ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FTH ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FRR ',FRR,3,NSA)
         ELSE IF(KID2.EQ.'CP ') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRACPA('FCP ',FCPP,1,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACPAB('FCP ',FCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CT ') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRACPA('FCT ',FCTH,2,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACPAB('FCT ',FCTH2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'EP ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FEP ',FEPP,1,NSA)
         ELSE IF(KID2.EQ.'ET ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('FET ',FETH,2,NSA)
         ELSE
            IF(nrank.EQ.0) WRITE(6,*) 'XX UNKNOWN KID3'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DPP ',DPP ,1,NSA)
         ELSE IF(KID2.EQ.'PT ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DPT ',DPT ,1,NSA)
         ELSE IF(KID2.EQ.'TP ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DTP ',DTP ,2,NSA)
         ELSE IF(KID2.EQ.'TT ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DTT ',DTT ,2,NSA)
         ELSE IF(KID2.EQ.'CPP') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRACPA('DCPP',DCPP,1,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACPAB('DCPP',DCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CPT') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRACPA('DCPT',DCPT,1,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACPAB('DCPT',DCPT2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTP') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRACPA('DCTP',DCTP,2,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACPAB('DCTP',DCTP2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTT') THEN
            IF(NSB.EQ.0) THEN
               IF(nrank.EQ.0) CALL FPGRACPA('DCTT',DCTT,2,NSA)
            ELSE
               IF(nrank.EQ.0) CALL FPGRACPAB('DCTT',DCTT2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'WPP') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DWPP',DWPP,1,NSA)
         ELSE IF(KID2.EQ.'WPT') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DWPT',DWPT,1,NSA)
         ELSE IF(KID2.EQ.'WTP') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DWTP',DWTP,2,NSA)
         ELSE IF(KID2.EQ.'WTT') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DWTT',DWTT,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            IF(nrank.EQ.0) CALL FPGRACPA('DRR ',DRR ,3,NSA)
         ELSE
            IF(nrank.EQ.0) WRITE(6,*) 'XX UNKNOWN KID3'
         ENDIF
      ELSE
         IF(nrank.EQ.0) WRITE(6,*) 'XX UNKNOWN KID2'
      END IF
      END IF
      CALL mtx_barrier
      GO TO 1

 9000 RETURN
    END SUBROUTINE FP_GOUT
!
! ***********************************************************
!
!                        R-GRAPHIC
!
! ***********************************************************
!
      SUBROUTINE FPGRARA(STRING,FRA,NSA)

      IMPLICIT NONE
      real(rkind),DIMENSION(NRMAX,NSAMAX,NTG2MAX):: FRA
      real(rkind),dimension(NRMAX,NTG2MAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NT2, NR, NSA

      DO NT2=1,NTG2
         DO NR=1,NRMAX
            TEMP(NR,NT2)=FRA(NR,NSA,NT2)
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAR(TRIM(STRING1),TEMP)
      RETURN
      END SUBROUTINE FPGRARA
!--------------------------------------------------------
      SUBROUTINE FPGRARAB(STRING,FRAB,NSA,NSB)

      IMPLICIT NONE
      real(rkind),DIMENSION(NRMAX,NSBMAX,NSAMAX,NTG2MAX):: FRAB
      real(rkind),dimension(NRMAX,NTG2MAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NT2, NR, NSA, NSB

      DO NT2=1,NTG2
         DO NR=1,NRMAX
            TEMP(NR,NT2)=FRAB(NR,NSB,NSA,NT2)
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAR(TRIM(STRING1),TEMP)
      RETURN
      END SUBROUTINE FPGRARAB
!--------------------------------------------------
      SUBROUTINE FPGRAR(STRING,FR)
!
      IMPLICIT NONE
      real(rkind),dimension(NRMAX,NTG2MAX)::FR
      real,dimension(NRMAX):: GX,GY
      CHARACTER(LEN=*),INTENT(IN):: STRING
      integer:: NT2, NR
      real:: GXMIN, GXMAX, GYMAX0, GYMIN0, GYMIN, GYMAX
      real:: GYMIN1, GYMAX1, GYSTEP, GXMIN1, GXMAX1, GXSTEP
      real:: GXORG
!
      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTR(STRING,FR)
         RETURN
      ENDIF

      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(0.3,0.)
      CALL SETFNT(32)

      DO NR=1,NRMAX
          GX(NR)=gdclip(RM(NR))
      ENDDO
      CALL GMNMX1(GX,1,NRMAX,1,GXMIN,GXMAX)
      IF(RGMIN.NE.0.0) GXMIN=gdclip(RGMIN)
      IF(RGMAX.NE.1.0) GXMIN=gdclip(RGMAX)

      GYMAX0=-1.E30
      GYMIN0= 1.E30
      DO NT2=1,NTG2
         DO NR=1,NRMAX
            GY(NR)=gdclip(FR(NR,NT2))
         END DO
         CALL GMNMX1(GY,1,NRMAX,1,GYMIN,GYMAX)
         GYMAX0=MAX(GYMAX0,GYMAX)
         GYMIN0=MIN(GYMIN0,GYMIN)
      END DO
      IF(GYMIN0 .GE. 0.) GYMIN0=0.
      IF(GYMAX0 .LE. 0.) GYMAX0=0.

      CALL GQSCAL(GYMIN0,GYMAX0,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(GXMIN,GXMAX,GXMIN1,GXMAX1,GXSTEP)
      CALL GDEFIN(3.0,23.0,2.0,17.0,GXMIN1,GXMAX1,GYMIN1,GYMAX1)
      CALL GFRAME
      GXORG=(INT(GXMIN1/(2*GXSTEP)+1))*2*GXSTEP
      CALL GSCALE(GXORG,GXSTEP,0.0,GYSTEP,0.1,9)    
      CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,ngslen(2*GXSTEP))
      CALL GVALUE(0.,0.0,0.0,2*GYSTEP,ngslen(2*GYSTEP))
      DO NT2=1,NTG2
         DO NR=1,NRMAX
            GY(NR)=gdclip(FR(NR,NT2))
            GX(NR)=gdclip(RM(NR))
         END DO
         CALL SETLIN(0,2,7-MOD(NT2-1,5))
         CALL GPLOTP(GX,GY,1,NRMAX,1,0,0,0)
      END DO
      CALL SETLIN(0,2,7)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(24.0,1.0)
      CALL TEXT('r',1)
      CALL PAGEE
!
      RETURN
      END SUBROUTINE FPGRAR
!
! ***********************************************************
!
!                        T-GRAPHIC
!
! ***********************************************************
!
      SUBROUTINE FPGRATA(STRING,FTA,NSA)

      IMPLICIT NONE
      real(rkind),DIMENSION(NSAMAX,NTG1MAX):: FTA
      real(rkind),dimension(NTG1MAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NT1, NSA

      DO NT1=1,NTG1
         TEMP(NT1)=FTA(NSA,NT1)
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAT(TRIM(STRING1),TEMP)
      RETURN
      END SUBROUTINE FPGRATA
!---------------------------------------------------
      SUBROUTINE FPGRATAB(STRING,FTAB,NSA,NSB)

      IMPLICIT NONE
      real(rkind),DIMENSION(NSBMAX,NSAMAX,NTG1MAX):: FTAB
      real(rkind),dimension(NTG1MAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NT1, NSA, NSB

      DO NT1=1,NTG1
         TEMP(NT1)=FTAB(NSB,NSA,NT1)
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAT(TRIM(STRING1),TEMP)
      RETURN
      END SUBROUTINE FPGRATAB
!---------------------------------------------------
      SUBROUTINE FPGRAT(STRING,FT)
!
      IMPLICIT NONE
      real(rkind),DIMENSION(NTG1MAX):: FT
      real,DIMENSION(NTG1MAX):: GX,GY
      CHARACTER(LEN=*),INTENT(IN):: STRING
      integer:: NT1
      real:: GYMIN, GYMAX, GYMIN1, GYMAX1, GYSTEP
      real:: GXMIN, GXMAX, GXSTEP

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTT(STRING,FT)
         RETURN
      ENDIF

      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
!
      DO  NT1=1,NTG1
         GX(NT1)=gdclip(PTG(NT1))
         GY(NT1)=gdclip(FT(NT1))
      END DO
      CALL GMNMX1(GY,1,NTG1,1,GYMIN,GYMAX)
      IF(GYMIN.GT.0.) GYMIN=0.0
      IF(GYMAX.LT.0.) GYMAX=0.0
      CALL GQSCAL(GYMIN,GYMAX,GYMIN1,GYMAX1,GYSTEP)
      CALL GQSCAL(0.,GX(NTG1),GXMIN,GXMAX,GXSTEP)
      CALL GDEFIN(3.,23.,2.,17.,0.,GX(NTG1),GYMIN1,GYMAX1)
      CALL SETLNW(0.035)
      CALL GSCALE(0.,GXSTEP,0.,GYSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GFRAME
      CALL GVALUE(0.,GXSTEP*2,0.,0.,ngslen(2*GXSTEP))
      CALL GVALUE(0.,0.,0.,GYSTEP*2,ngslen(2*GYSTEP))
      CALL GPLOTP(GX,GY,1,NTG1,1,0,0,0)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(24.0,1.0)
      CALL TEXT('t',1)
      CALL PAGEE
      RETURN
      END SUBROUTINE FPGRAT
!
! ***********************************************************
!
!                        P-GRAPHIC
!
! ***********************************************************
!
      SUBROUTINE FPGRAPA(STRING,FGA,NSA)

      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX,NSAMAX):: FGA
      real(rkind),dimension(NTHMAX,NPMAX,NRMAX)::TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSA)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAP(TRIM(STRING1),TEMP,NSA)
      RETURN
      END SUBROUTINE FPGRAPA
!--------------------------------------------------------------
      SUBROUTINE FPGRAPAB(STRING,FGAB,NSA,NSB)

      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX,NSBMAX,NSAMAX):: FGAB
      real(rkind),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA, NSB

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGAB(NTH,NP,NR,NSB,NSA)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAP(TRIM(STRING1),TEMP,NSA)
      RETURN
      END SUBROUTINE FPGRAPAB
!---------------------------------------------------------------
      SUBROUTINE FPGRAP(STRING,FG,NSA)
!
      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      real,dimension(NPMAX):: GX
      real,dimension(NPMAX,NRMAX):: GY
      CHARACTER(LEN=*),INTENT(IN):: STRING
      INTEGER,INTENT(IN):: NSA
      CHARACTER(LEN=80):: STRING1
!      real(rkind):: PGMAX, RGMAX, RGMIN
      real:: GXMAX, GXMINP, GYMIN, GYMAX, GYMINP, GYSTEP, GXSTEP, GXMAXP, GYMAXP
      integer:: NR, NP, NTH, NPM, NS

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTP(STRING,FG)
         RETURN
      ENDIF

    1 WRITE(6,'(A,I4,A)') '# INPUT NTH (1..',NTHMAX,' or 0) :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
      WRITE(STRING1,'(A,A,I4)') STRING,' : NTH=',NTH
      NS=NS_NSA(NSA)

      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(0.3,0.)
      CALL SETFNT(32)
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            IF(FG(NTH,NP,NR).LT.1.D-14) THEN
               GY(NP,NR)=-14.0
            ELSE
               GY(NP,NR)=gdclip(LOG10(FG(NTH,NP,NR)))
            ENDIF
         END DO
      END DO
      DO NP=1,NPMAX
         GX(NP)=gdclip(PM(NP,NS)**2)
      ENDDO
      GXMAX=gdclip(PMAX(NS)**2)
      IF(PGMAX.NE.0.0) GXMAX=gdclip(PGMAX**2)

      CALL GMNMX2(GY,NPM,1,NPMAX,1,1,NRMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GYMIN,GYMAX,GYMINP,GYMAXP,GYSTEP)
      CALL GQSCAL(0.0,GXMAX,GXMINP,GXMAXP,GXSTEP)
!
      CALL GDEFIN(3.0,23.0,2.0,17.0,0.0,GXMAX,-8.0,0.0)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,GXSTEP*2,0.0,0.0,1.0,0)
      CALL SETLNW(0.07)
      CALL GSCALE(0.,GXSTEP,0.0,0.0,0.2,9)
      CALL GSCALL(0.,0,0.0,1,1.0,0)
      CALL GSCALL(0.,0,0.0,2,0.2,9)
      CALL GFRAME
      CALL GVALUE(0.,GXSTEP*2,0.0,0.0,ngslen(2*GXSTEP))
      CALL GVALUL(0.,0.0,0.0,1,0)
      DO NR=1,NRMAX
         CALL SETLIN(0,0,7-MOD(NR-1,5))
         CALL GPLOTP(GX,GY(1,NR),1,NPMAX,1,0,0,0)
      END DO
      CALL SETLIN(0,0,7)
      CALL MOVE(1.0,17.5)
      CALL TEXT(STRING1,LEN(STRING1))
      CALL MOVE(24.0,1.0)
      CALL TEXT('p**2',4)
      CALL PAGEE
!
 9000 RETURN
      END SUBROUTINE FPGRAP
!
! ***********************************************************
!
!                        PR-GRAPHIC
!
! ***********************************************************
!
      SUBROUTINE FPGRAPRA(STRING,FGA,NSA)

      IMPLICIT NONE
      REAL(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX,NSAMAX):: FGA
      REAL(rkind),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSA)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAPR(TRIM(STRING1),TEMP,NSA)
      RETURN
      END SUBROUTINE FPGRAPRA
!-------------------------------------------
      SUBROUTINE FPGRAPRAB(STRING,FGAB,NSA,NSB)

      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX,NSBMAX,NSAMAX):: FGAB
      real(rkind),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA, NSB

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FGAB(NTH,NP,NR,NSB,NSA)
            ENDDO
         ENDDO
      ENDDO
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAPR(TRIM(STRING1),TEMP,NSA)
      RETURN
      END SUBROUTINE FPGRAPRAB
!-------------------------------------------------------
      SUBROUTINE FPGRAPR(STRING,FG,NSA)
!       
      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      real,dimension(NPMAX):: GX
      real,dimension(NRMAX):: GY
      real,dimension(NPMAX,NRMAX):: GZ
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA, NS, NPGMAX, NPM
      real:: GX1, GX2, GY1, GY2

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTP(STRING,FG)
         RETURN
      ENDIF
      NS=NS_NSA(NSA)

    1 WRITE(6,'(A,I4,A)') '# INPUT NTH (1..',NTHMAX,' or 0) :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
      WRITE(STRING1,'(A,A,I4)') STRING,' : NTH=',NTH

      NPGMAX=1
      DO NP=1,NPMAX
         IF(PGMAX.NE.0.0) THEN
            IF(PM(NP,NS).GT.PGMAX) GOTO 10
         ENDIF
         NPGMAX=NP
         GX(NP)=gdclip(PM(NP,NS)**2)
      ENDDO
   10 CONTINUE
      DO NR=1,NRMAX
         GY(NR)=gdclip(RM(NR))
      ENDDO
      DO NR=1,NRMAX
      DO NP=1,NPMAX
         IF(FG(NTH,NP,NR).LT.1.D-14) THEN
            GZ(NP,NR)=-14.0
         ELSE
            GZ(NP,NR)=gdclip(LOG10(FG(NTH,NP,NR)))
         ENDIF
      ENDDO
      ENDDO

      GX1=3.0
      GX2=20.0
      GY1=2.0
      GY2=17.0

      CALL PAGES
      CALL FPGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NPM,NPGMAX,NRMAX,  &
           gdclip(PGMAX**2),gdclip(RGMIN),gdclip(RGMAX),TRIM(STRING1))
      CALL PAGEE

 9000 RETURN
      END SUBROUTINE FPGRAPR

!     ***********************************************************

!           SUBPROGRAM FOR 2D PROFILE

!     ***********************************************************

      SUBROUTINE FPGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX, &
                        GXMAX1,GYMIN1,GYMAX1,STR)

      IMPLICIT NONE
      REAL,        INTENT(IN):: GX1,GX2,GY1,GY2,GXMAX1,GYMIN1,GYMAX1
      INTEGER,     INTENT(IN):: NXM,NXMAX,NYMAX
      REAL,DIMENSION(NXMAX),      INTENT(IN):: GX
      REAL,DIMENSION(NYMAX),      INTENT(IN):: GY
      REAL,DIMENSION(NXM,NYMAX),  INTENT(IN):: GZ
      CHARACTER(LEN=*),             INTENT(IN):: STR
      INTEGER :: I
      REAL    :: GOX, GOY, GOZ, GPHI, GRADIUS, GSTEPX, GSTEPY,   &
                    GSTEPZ, GSXMAX, GSXMIN, GSYMAX, GSYMIN, GSZMAX, &
                    GSZMIN, GTHETA, GXL, GXMAX, GXMIN, GXORG, GYL,  &
                    GYMAX, GYMIN, GYORG, GZL, GZMAX, GZMIN, GZVAL
      EXTERNAL R2W2B,R2Y2W,W2G2B,WHITE

      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.07)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLIN(-1,-1,7)

      CALL MOVE(GX1,GY2+0.2)
      CALL TEXT(STR,LEN(STR))

      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
      GZVAL=0.5*(ABS(GZMIN)+ABS(GZMAX))
      IF((GZVAL.LE.1.D-6).OR.(ABS(GZMAX-GZMIN)/GZVAL.LT.1.D-6)) THEN
         WRITE(6,*) GZMIN,GZMAX
         CALL GUFLSH
         RETURN
      ENDIF

      IF(GZMIN.GE.0.0) THEN
         GZMIN=0.0
      ELSEIF(GZMAX.LE.0.0) THEN
         GZMAX=0.0
      ENDIF

      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)

      IF(GXMAX1.NE.0.0) GXMAX=GXMAX1
      IF(GYMIN1.NE.0.0) GYMIN=GYMIN1
      IF(GYMAX1.NE.1.0) GYMAX=GYMAX1

      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
!      WRITE(6,*) GSZMIN,GSZMAX,GSTEPZ
!      CALL GUFLSH

      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GSXMIN
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GSYMIN
      ENDIF
      IF(GZMIN*GZMAX.LT.0.0) THEN
         GSZMAX=MAX(ABS(GSZMIN),ABS(GSZMAX))
         GSZMIN=-GSZMAX
      ENDIF

      GXL=10.0*1.5
      GYL=20.0*1.5
      GZL=10.0*1.5
      CALL GDEFIN3D(GX1,GX2,GY1,GY2,GXL,GYL,GZL)
      GPHI=-60.0
      GTHETA=65.0
      GRADIUS=100.0
      GOX=0.5*(GSXMIN+GSXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GSZMIN+GSZMAX)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,1.0,1,GOX,GOY,GOZ)
      CALL GDATA3D3(GZ,GX,GY,NXM,NXMAX,NYMAX, &
                    GSXMIN,GSXMAX,GSYMIN,GSYMAX,GSZMIN,GSZMAX)
      CALL SETCHS(0.2,0.0)
      CALL SETLIN(0,0,7)

      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
      CALL GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
      CALL GSCALE3DZ(GSZMIN,GSTEPZ,0.3,0)
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,ngslen(GSTEPX))
      CALL GVALUE3DY(GSYMIN,GSTEPY,1,ngslen(GSTEPY))
      CALL GVALUE3DZ(GSZMIN,GSTEPZ,2,ngslen(GSTEPZ))

      IF(GZMIN*GZMAX.LT.0.0) THEN
         CALL CPLOT3D1(7,R2W2B)
      ELSEIF(GZMIN.LT.0.0) THEN
         CALL CPLOT3D1(7,W2G2B)
!         CALL CPLOT3D1(7,WHITE)
      ELSE
         CALL CPLOT3D1(7,R2Y2W)
      ENDIF

      CALL GAXIS3D(0)

      CALL SETLIN(0,0,7)
      RETURN
      END SUBROUTINE FPGR3D
!
! ***********************************************************
!
!                        C-GRAPHIC
!
! ***********************************************************

      SUBROUTINE FPGRACA(STRING,FGA,MODE,NSA)
       
      IMPLICIT NONE
      REAL(rkind),DIMENSION(:,:,:,:):: FGA
      REAL(rkind),dimension(NTHMAX+1,NPMAX+1,NRMAX+1):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA, NS
      integer:: MODE
      NS=NS_NSA(NSA)
!
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.3) THEN
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRAC(TRIM(STRING1),TEMP,MODE,NSA)
      RETURN
      END SUBROUTINE FPGRACA
!---------------------------------------------------
      SUBROUTINE FPGRACXA(STRING,FGA,MODE,NSA)
       
      IMPLICIT NONE
      REAL(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NSAMAX):: FGA
      REAL,dimension(NTHMAX+1,NPMAX+1):: GF
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NPM, NTHM, NRM, NR, NP, NTH, NSA, NS
      integer:: NM, NM1, MODE
      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1
      NS=NS_NSA(NSA)
!
      IF(MODE.EQ.0) THEN
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            GF(NP,NTH)=gdclip(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            GF(NP,NTH)=gdclip(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            GF(NP,NTH)=gdclip(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.3) THEN
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            GF(NP,NTH)=gdclip(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               IF(FGA(NTH,NP,NSA).LT.1.D-70) THEN
                  GF(NP,NTH)=-70.0
               ELSE
                  GF(NP,NTH)=gdclip(LOG10(ABS(FGA(NTH,NP,NSA))))
               ENDIF
            END DO
         END DO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRACX(TRIM(STRING1),GF,MODE,NSA)
      RETURN
      END SUBROUTINE FPGRACXA
!---------------------------------------------------
      SUBROUTINE FPGRACB(STRING,FGA,MODE,NSB)
       
      IMPLICIT NONE
      REAL(rkind),DIMENSION(:,:,:,:):: FGA
!      REAL(rkind),dimension(NTHMAX+1,NPMAX+1,NRMAX+1):: TEMP
      REAL(rkind),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSB
      integer:: NM, NM1, MODE
!
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSB)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSB)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSB)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.3) THEN
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            TEMP(NTH,NP,NR)=FGA(NTH,NP,NR,NSB)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSB,')'
      CALL FPGRAC_2(TRIM(STRING1),TEMP,MODE,NSB)
      RETURN
      END SUBROUTINE FPGRACB
!---------------------------------------------------
      SUBROUTINE FPGRACAB(STRING,FGAB,MODE,NSA,NSB)
       
      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1,NSBMAX,NSAMAX):: FGAB
      real(rkind),dimension((NRMAX+1)*(NTHMAX+1)*(NPMAX+1)):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA,NSB,NS,NM, NM1, NPM, NTHM, NRM, NSBM, MODE
      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1
      NSBM=NSBMAX
      NS=NS_NSA(NSA)
!
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.3) THEN
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1,I2,A1)') STRING,'(',NSA,',',NSB,')'
      CALL FPGRAC(TRIM(STRING1),TEMP,MODE,NSA)
      RETURN
      END SUBROUTINE FPGRACAB
!--------------------------------------------------
      SUBROUTINE FPGRAC(STRING,FG,MODE,NSB)
!       
      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1):: FG
      real,DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real,dimension(NPMAX+1):: GP
      real,dimension(NTHMAX+1):: GTH
      real(rkind),dimension(8,NPMAX+1,NTHMAX+1)::KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL:: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL:: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real:: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real:: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
      integer:: NR, NP, NTH, NSB, NRG, MODE
      integer:: NPG, NTHG, NPM, NTHM, NRM, NM, NGLMAX, NGL
      integer:: I

      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTC(STRING,FG,MODE)
         RETURN
      ENDIF

    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.3) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
      WRITE(STRING1,'(A,A,I4)') STRING,' : NR=',NR
!
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.3) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               IF(FG(NTH,NP,NR).LT.1.D-70) THEN
                  GF(NP,NTH)=-70.0
               ELSE
                  GF(NP,NTH)=gdclip(LOG10(ABS(FG(NTH,NP,NR))))
               ENDIF
            END DO
         END DO
      ENDIF
!
      CALL FPGRACX(STRING1,GF,MODE,NSB)

      IF(NRMAX.GT.1) GOTO 1
!
 9000 RETURN
      END SUBROUTINE FPGRAC
!--------------------------------------------------
      SUBROUTINE FPGRAC_2(STRING,FG,MODE,NSB)
!       
      IMPLICIT NONE
!      real(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1):: FG
!      real,DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      real,DIMENSION(NPMAX,NTHMAX):: GF
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL:: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL:: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real:: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real:: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
      integer:: NR, NP, NTH, NRG, MODE, NSB
      integer:: NPG, NTHG, NPM, NTHM, NRM, NM, NGLMAX, NGL
      integer:: I

      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1

      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTC_2(STRING,FG,MODE)
         RETURN
      ENDIF

    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.3) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
      WRITE(STRING1,'(A,A,I4)') STRING,' : NR=',NR
!
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.3) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               GF(NP,NTH)=gdclip(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               IF(FG(NTH,NP,NR).LT.1.D-70) THEN ! TEST
                  GF(NP,NTH)=-70.0
!               IF(FG(NTH,NP,NR).LT.1.D-50) THEN ! TEST
!                  GF(NP,NTH)=-50.0
!               IF(FG(NTH,NP,NR).LT.1.D-12) THEN 
!                  GF(NP,NTH)=-12.0
               ELSE
                  GF(NP,NTH)=gdclip(LOG10(ABS(FG(NTH,NP,NR))))
               ENDIF
            END DO
         END DO
      ENDIF
!
      CALL FPGRACX_2(STRING,GF,MODE,NSB)

      IF(NRMAX.GT.1) GOTO 1
!
 9000 RETURN
      END SUBROUTINE FPGRAC_2
!--------------------------------------------------
      SUBROUTINE FPGRACX_2(STRING,GF,MODE,NSB)
!       
      IMPLICIT NONE
      real,DIMENSION(NPMAX,NTHMAX):: GF
      real,dimension(NPMAX):: GP
      real,dimension(NTHMAX):: GTH
      real(rkind),dimension(8,NPMAX,NTHMAX)::KA
!      real,DIMENSION(NPMAX+1,NTHMAX+1):: GF
!      real,dimension(NPMAX+1):: GP
!      real,dimension(NTHMAX+1):: GTH
!      real(rkind),dimension(8,NPMAX+1,NTHMAX+1)::KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      INTEGER,PARAMETER:: NGLM=30
      REAL:: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL:: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real:: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real:: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
      integer:: NR, NP, NTH, NSB, NRG, MODE, LMODE
      integer:: NPG, NTHG, NPM, NTHM, NGLMAX, NGL
      integer:: I, NS

      NPM=NPMAX+1
      NTHM=NTHMAX+1

      LMODE=MODE/4
      NS=NS_NSA(NSB)
!
      IF(MODE.EQ.1) THEN
         DO NP=1,NPMAX+1
            GP(NP)=gdclip(PG(NP,NS))
         END DO
         NPG=NPMAX+1
      ELSE
         DO NP=1,NPMAX
            GP(NP)=gdclip(PM(NP,NS))
         END DO
         NPG=NPMAX
      ENDIF

      IF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            GTH(NTH)=gdclip(THG(NTH))
         END DO
         NTHG=NTHMAX+1
      ELSE
         DO NTH=1,NTHMAX
            GTH(NTH)=gdclip(THM(NTH))
         END DO
         NTHG=NTHMAX
      ENDIF
!
      GPMAX=gdclip(PMAX(NS))
!
      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
!
      CALL GMNMX2(GF,NPM-1,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
!      CALL GMNMX2(GF,NPM,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
!
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,ngslen(2*GPSTEP))
!
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=GFMIN1+0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,     &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
!c$$$               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                     GFMIN1,0.5*GFSTEP,30,0,KA)
            ELSE
               NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=GFMAX1-0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=1.D0
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,     &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
!c$$$               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                     GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
               GFFMAX=MAX(ABS(GFMAX1),ABS(GFMIN1))
               NGLMAX=INT(GFFMAX/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=-0.25*GFSTEP-0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=1.D0
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,   &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
               DO NGL=1,NGLMAX
                  ZL(NGL)= 0.25*GFSTEP+0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,   &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
!c$$$            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                   0.25*GFSTEP, 0.5*GFSTEP,15,0,KA)
!c$$$            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                  -0.25*GFSTEP,-0.5*GFSTEP,15,2,KA)
         ENDIF
      ELSE
         DO NGL=1,NGLINE
            ZL(NGL)=GFMAX-0.020*(NGL-1)**2
            CALL SETRGBFP(1.0-FLOAT(NGL-1)/FLOAT(NGLINE-1),RGB(1,NGL))
!            WRITE(6,'(I5,1P5E12.4)') &
!               NGL,ZL(NGL),1.0-FLOAT(NGL-1)/FLOAT(NGLINE-1), &
!               RGB(1,NGL),RGB(2,NGL),RGB(3,NGL)
            ILN(NGL)=0
            WLN(NGL)=0.07
         ENDDO
!         CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,ZL,RGB,ILN,WLN,NGLINE,0)
         CALL CONTG4X(GF,GP,GTH,NPG,NPG,NTHG,ZL,RGB,ILN,WLN,NGLINE,0)

!            CALL SETLIN(0,0,7-MOD(I-1,5))
!            CALL SETLNW(0.07)
!            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,        &
!                 GLIN,GFSTEP*100,1,0,KA)

!            CALL CONTQ4(GF,GP,GTH,NPM-1,NPG,NTHG,        &
!                 GLIN,GFSTEP*100,1,0,KA)
      ENDIF
!
      CALL SETLIN(0,2,7)
      CALL MOVE(24.0,1.0)
      CALL TEXT('PPARA',5)
      CALL MOVE(1.0,13.5)
      CALL TEXT('PPERP',5)
!
      CALL MOVE(3.0,12.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(8.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(13.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(18.0,12.5)
         CALL TEXT('STEP =',6)
         CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
      ENDIF
      CALL PAGEE
!
      RETURN
      END SUBROUTINE FPGRACX_2
!--------------------------------------------------
      SUBROUTINE FPGRACX(STRING,GF,MODE,NSB)
!       
      IMPLICIT NONE
      real,DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real,dimension(NPMAX+1):: GP
      real,dimension(NTHMAX+1):: GTH
      real(rkind),dimension(8,NPMAX+1,NTHMAX+1)::KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      INTEGER,PARAMETER:: NGLM=30
      REAL:: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL:: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real:: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real:: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
      integer:: NR, NP, NTH, NSB, NRG, MODE, LMODE
      integer:: NPG, NTHG, NPM, NTHM, NRM, NM, NGLMAX, NGL
      integer:: I, NS

      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1

      LMODE=MODE/4
      NS=NS_NSB(NSB)
!
      IF(MODE.EQ.1) THEN
         DO NP=1,NPMAX+1
            GP(NP)=gdclip(PG(NP,NS))
         END DO
         NPG=NPMAX+1
      ELSE
         DO NP=1,NPMAX
            GP(NP)=gdclip(PM(NP,NS))
         END DO
         NPG=NPMAX
      ENDIF

      IF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            GTH(NTH)=gdclip(THG(NTH))
         END DO
         NTHG=NTHMAX+1
      ELSE
         DO NTH=1,NTHMAX
            GTH(NTH)=gdclip(THM(NTH))
         END DO
         NTHG=NTHMAX
      ENDIF
!
      GPMAX=gdclip(PMAX(NS))
!
      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
!
      CALL GMNMX2(GF,NPM,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
!
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,ngslen(2*GPSTEP))
!
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=GFMIN1+0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,     &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
!c$$$               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                     GFMIN1,0.5*GFSTEP,30,0,KA)
            ELSE
               NGLMAX=INT((GFMAX1-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=GFMAX1-0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=1.D0
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,     &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
!c$$$               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                     GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
               GFFMAX=MAX(ABS(GFMAX1),ABS(GFMIN1))
               NGLMAX=INT(GFFMAX/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  ZL(NGL)=-0.25*GFSTEP-0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=1.D0
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,   &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
               DO NGL=1,NGLMAX
                  ZL(NGL)= 0.25*GFSTEP+0.5*GFSTEP*(NGL-1)
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,   &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
!c$$$            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                   0.25*GFSTEP, 0.5*GFSTEP,15,0,KA)
!c$$$            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,
!c$$$     &                  -0.25*GFSTEP,-0.5*GFSTEP,15,2,KA)
         ENDIF
      ELSE
         DO NGL=1,NGLINE
            ZL(NGL)=GFMAX-0.020*(NGL-1)**2
            CALL SETRGBFP(1.0-FLOAT(NGL-1)/FLOAT(NGLINE-1),RGB(1,NGL))
!            WRITE(6,'(I5,1P5E12.4)') &
!               NGL,ZL(NGL),1.0-FLOAT(NGL-1)/FLOAT(NGLINE-1), &
!               RGB(1,NGL),RGB(2,NGL),RGB(3,NGL)
            ILN(NGL)=0
            WLN(NGL)=0.07
         ENDDO
         CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG,ZL,RGB,ILN,WLN,NGLINE,0)
!            CALL SETLIN(0,0,7-MOD(I-1,5))
!            CALL SETLNW(0.07)
!            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,        &
!                 GLIN,GFSTEP*100,1,0,KA)
      ENDIF
!
      CALL SETLIN(0,2,7)
      CALL MOVE(24.0,1.0)
      CALL TEXT('PPARA',5)
      CALL MOVE(1.0,13.5)
      CALL TEXT('PPERP',5)
!
      CALL MOVE(3.0,12.5)
      CALL TEXT(STRING,LEN(STRING))
      CALL MOVE(8.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(13.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(18.0,12.5)
         CALL TEXT('STEP =',6)
         CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
      ENDIF
      CALL PAGEE
!
      RETURN
      END SUBROUTINE FPGRACX

      SUBROUTINE SETRGBFP(F,RGB)
        IMPLICIT NONE
        REAL,INTENT(IN):: F
        REAL,DIMENSION(3),INTENT(OUT):: RGB
        INTEGER,PARAMETER:: NFMAX=8
        REAL,DIMENSION(3,NFMAX):: RGBC
        DATA RGBC/ 0.0,0.0,0.0, &
                   0.0,0.0,1.0, &
                   0.0,0.8,1.0, &
                   0.0,0.8,0.0, &
                   1.0,0.8,0.0, &
                   1.0,0.4,0.0, &
                   1.0,0.0,0.0, &
                   1.0,1.0,1.0/
        REAL(rkind):: GF,DF
        INTEGER:: IM
!
        GF=F*DBLE(NFMAX-1)+1
        IM=MIN(INT(GF),NFMAX-1)
        DF=GF-IM
        RGB(1)=RGBC(1,IM)*(1.D0-DF)+RGBC(1,IM+1)*DF
        RGB(2)=RGBC(2,IM)*(1.D0-DF)+RGBC(2,IM+1)*DF
        RGB(3)=RGBC(3,IM)*(1.D0-DF)+RGBC(3,IM+1)*DF
        RETURN
        END SUBROUTINE SETRGBFP
          
!
! ***********************************************************
!
!                        C-GRAPHIC muliplied p
!
! ***********************************************************
!
      SUBROUTINE FPGRACPA(STRING,FGA,MODE,NSA)
       
      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1,NSAMAX):: FGA
      real(rkind),dimension((NRMAX+1)*(NPMAX+1)*(NTHMAX+1)):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA, MODE, NPM, NTHM, NRM, NM, NM1
      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1
!
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*(NSA-1)+NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGA(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRACP(TRIM(STRING1),TEMP,MODE,NSA)
      RETURN
      END SUBROUTINE FPGRACPA
!------------------------------------------------------------
      SUBROUTINE FPGRACPAB(STRING,FGAB,MODE,NSA,NSB)
       
      IMPLICIT NONE
      real(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1,NSBMAX,NSAMAX):: FGAB
      real(rkind),dimension((NRMAX+1)*(NPMAX+1)*(NTHMAX+1)):: TEMP
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NTH, NP, NSA, NSB, MODE, NPM, NTHM, NRM, NM, NM1, NSBM
      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1
      NSBM=NSBMAX
!
      IF(MODE.EQ.0.OR.MODE.EQ.4) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1) &
               +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX
         DO NP=1,NPMAX+1
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1) &
               +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
         DO NTH=1,NTHMAX+1
         DO NP=1,NPMAX
            NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            NM1=NPM*NTHM*NRM*NSBM*(NSA-1)+NPM*NTHM*NRM*(NSB-1) &
               +NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
            TEMP(NM)=FGAB(NTH,NP,NR,NSB,NSA)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      WRITE(STRING1,'(A,A1,I2,A1)') STRING,'(',NSA,')'
      CALL FPGRACP(TRIM(STRING1),TEMP,MODE,NSA)
      RETURN
      END SUBROUTINE FPGRACPAB
!-------------------------------------------------
      SUBROUTINE FPGRACP(STRING,FG,MODE,NSA)
!
      IMPLICIT NONE
      real(rkind),DIMENSION((NRMAX+1)*(NPMAX+1)*(NTHMAX+1)):: FG
      real,DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real,dimension(NPMAX+1):: GP
      real,dimension(NTHMAX+1):: GTH
      integer,dimension(8,NPMAX+1,NTHMAX+1):: KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL:: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL:: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real:: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real:: GPMIN1, GPMAX1, GPSTEP, GLIN
      integer:: NR, NTH, NP, NSA, MODE, NRG, LMODE, NTHG, NS
      integer:: NPM, NTHM, NRM, NM, NGL, NGLMAX, I, NPG
      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1

      NS=NS_NSA(NSA)
      IF(NGRAPH.EQ.0) THEN
         CALL FPFOTC(STRING,FG,MODE)
         RETURN
      ENDIF

    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
      WRITE(STRING1,'(A,A,I4)') STRING,' : NR=',NR
!
      LMODE=MODE/4
!
      IF(MOD(MODE,2).EQ.0) THEN
         DO NP=1,NPMAX
            GP(NP)=gdclip(PM(NP,NS))
         END DO
         NPG=NPMAX
      ELSE
         DO NP=1,NPMAX+1
            GP(NP)=gdclip(PG(NP,NS))
         END DO
         NPG=NPMAX+1
      ENDIF
!
      IF(MOD(MODE/2,2).EQ.0) THEN
         DO NTH=1,NTHMAX
            GTH(NTH)=gdclip(THM(NTH))
         END DO
         NTHG=NTHMAX
      ELSE
         DO NTH=1,NTHMAX+1
            GTH(NTH)=gdclip(THG(NTH))
         END DO
         NTHG=NTHMAX+1
      ENDIF
!
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               GF(NP,NTH)=gdclip(FG(NM)*PM(NP,NS))
            END DO
         END DO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               GF(NP,NTH)=gdclip(FG(NM)*PG(NP,NS))
            END DO
         END DO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               GF(NP,NTH)=gdclip(FG(NM)*PM(NP,NS))
            END DO
         END DO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               IF(FG(NM)*PM(NP,NS).LT.1.D-70) THEN
                  GF(NP,NTH)=-70.0
               ELSE
                  GF(NP,NTH)=gdclip(LOG10(ABS(FG(NM)*PM(NP,NS))))
               ENDIF
            END DO
         END DO
      ENDIF
      GPMAX=gdclip(PMAX(NS))

      CALL PAGES
      CALL SETLNW(0.07)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
!
      CALL GMNMX2(GF,NPM,1,NPG,1,1,NTHG,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFMIN1,GFMAX1,GFSTEP)
      CALL GQSCAL(0.0,GPMAX,GPMIN1,GPMAX1,GPSTEP)
!
      CALL GDEFIN(3.,23.,2.,12.,-GPMAX,GPMAX,0.,GPMAX)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,GPSTEP,0.,GPSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,ngslen(2*GPSTEP))
!
      IF(LMODE.EQ.0) THEN
         IF(GFMIN*GFMAX.GE.0.0) THEN
            IF(GFMIN.GE.0.0) THEN
               NGLMAX=INT((GFMAX-GFMIN1)/(0.5*GFSTEP))
               DO NGL=1,NGLMAX
                  RGB(1,NGL)=1.D0
                  RGB(2,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  RGB(3,NGL)=0.9*FLOAT(NGLMAX-NGL)/FLOAT(NGLMAX-1)
                  ILN(NGL)=0
                  WLN(NGL)=0.07
               ENDDO
               CALL CONTG4X(GF,GP,GTH,NPM,NPG,NTHG, &
                           ZL,RGB,ILN,WLN,NGLMAX,0)
            ELSE
               CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,  &
                           GFMIN1,0.5*GFSTEP,30,2,KA)
            ENDIF
         ELSE
            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,     &
                         0.25*GFSTEP, 0.5*GFSTEP,15,0,KA)
            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,     &
                        -0.25*GFSTEP,-0.5*GFSTEP,15,2,KA)
         ENDIF
      ELSE
         DO I=1,NGLINE
            GLIN=GFMAX-0.020*(I-1)**2
            CALL SETLIN(0,0,7-MOD(I-1,5))
            CALL SETLNW(0.07)
            CALL CONTQ4(GF,GP,GTH,NPM,NPG,NTHG,  &
                 GLIN,GFSTEP*100,1,0,KA)
         END DO
      ENDIF
!
      CALL SETLIN(0,2,7)
      CALL MOVE(24.0,1.0)
      CALL TEXT('PPARA',5)
      CALL MOVE(1.0,13.5)
      CALL TEXT('PPERP',5)
!
      CALL MOVE(3.0,12.5)
      CALL TEXT(STRING1,LEN(STRING1))
      CALL MOVE(8.0,12.5)
      CALL TEXT('FMIN =',6)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(13.0,12.5)
      CALL TEXT('FMAX =',6)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      IF(LMODE.EQ.0) THEN
         CALL MOVE(18.0,12.5)
         CALL TEXT('STEP =',6)
         CALL NUMBR(0.5*GFSTEP,'(1PE12.4)',12)
      ENDIF
      CALL PAGEE
      IF(NRMAX.GT.1) GOTO 1
!
 9000 RETURN
      END SUBROUTINE FPGRACP

      SUBROUTINE FPGRACD(STRING,FG,MODE,NSA)

        IMPLICIT NONE
        CHARACTER(LEN=*),INTENT(IN):: STRING
        REAL(rkind),INTENT(IN),DIMENSION(NTHMAX,NPMAX,NRMAX,NSAMAX):: FG
        INTEGER,INTENT(IN):: MODE,NSA
        REAL(rkind),ALLOCATABLE:: TEMP(:,:,:)
        INTEGER:: NR,NP,NTH

        ALLOCATE(TEMP(NTHMAX+1,NPMAX+1,NRMAX+1))
        
        DO NR=1,NRMAX
           DO NP=1,NPMAX
              DO NTH=1,NTHMAX
                 TEMP(NTH,NP,NR)=FG(NTH,NP,NR,NSA)         &
                                -FG(NTHMAX+1-NTH,NP,NR,NSA)
              ENDDO
           ENDDO
        ENDDO

        CALL FPGRAC(STRING,TEMP,MODE,NSA)

        DEALLOCATE(TEMP)

        RETURN
      END SUBROUTINE FPGRACD

    END MODULE FPGOUT
