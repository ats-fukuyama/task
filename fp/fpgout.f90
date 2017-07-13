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
      interface
         real(4) function GUCLIP(X)
           real(8):: X
         end function GUCLIP
         integer(4) function NGULEN(Y)
           real(4):: Y
         end function NGULEN
      end interface

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

      IMPLICIT NONE
      integer,SAVE:: NSA=1,NSB=0
      DATA NSA,NSB/1,0/
      integer:: NR, NP, NTH, NSA1, NS
!
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX)::TEMP
      CHARACTER KID*5,KID1*1,KID2*3
      CHARACTER(LEN=80):: STRING1
!
    1 WRITE(6,*)'INPUT GRAPH TYPE : NSA,NSB=',NSA,NSB
      WRITE(6,*)'    : F/FX/FR/FS1/FS2 1/2, Nn:NSA=n, Nnm:NSA=n,NSB=m,'
      WRITE(6,*)'    : D/DC/DW PP/PT/TP/TT/RR, F/FC/FE P/T/R'
      WRITE(6,*)'    : PD/PDC/PDW PP/PT/TP/TT/RR, PF/PFC/FE P/T/R'
      WRITE(6,*)'    : R/T N/I/W/PC/PW/PE/T/Q/E,  Gn,  X:exit'
      READ(5,'(A5)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      CALL GUCPTL(KID1)

      IF(KID1.NE.'P') THEN
         KID2=KID(2:4)
         CALL GUCPTL(KID2(1:1))
         CALL GUCPTL(KID2(2:2))
         CALL GUCPTL(KID2(3:3))
!
      NS=NS_NSA(NSA)
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'1  ') THEN
            CALL FPGRAPA('F1',FNS,NSA)
         ELSE IF(KID2.EQ.'R1 ') THEN
            CALL FPGRAPRA('FR1',FNS,NSA)
         ELSE IF(KID2.EQ.'2  ') THEN
            CALL FPGRACB('F2',FNS,4,NSA)
         ELSE IF(KID2.EQ.'X2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FNS(NTH,NP,NR,NSA)         &
                              -FNS(NTHMAX+1-NTH,NP,NR,NSA)
            ENDDO
            ENDDO
            ENDDO
            WRITE(STRING1,'(A,A1,I2,A1)') 'FX2','(',NSA,')'
            CALL FPGRAC(TRIM(STRING1),TEMP,4,NSA)
         ELSE IF(KID2.EQ.'Y2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=FNS(NTH,NP,NR,NSA)         &
                              -FNS(NTHMAX+1-NTH,NP,NR,NSA)
            ENDDO
            ENDDO
            ENDDO
            WRITE(STRING1,'(A,A1,I2,A1)') 'FY2','(',NSA,')'
            CALL FPGRAC(TRIM(STRING1),TEMP,0,NSA)
         ELSE IF(KID2.EQ.'S11') THEN
            CALL FPGRAPA('FS11',FS1,NSA)
         ELSE IF(KID2.EQ.'S12') THEN
            CALL FPGRACXA('FS12',FS1,4,NSA)
         ELSE IF(KID2.EQ.'S21') THEN
            CALL FPGRAPA('FS21',FS2,NSA)
         ELSE IF(KID2.EQ.'S22') THEN
            CALL FPGRACXA('FS22',FS2,4,NSA)
         ELSE IF(KID2.EQ.'P  ') THEN
            CALL FPGRACA('FP  ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRACA('FT  ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRACA('FR  ',FRR,0,NSA)
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL FPGRACA('FPP ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL FPGRACA('FTH ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACA('FRR ',FRR,0,NSA)
         ELSE IF(KID2.EQ.'CP ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('FCP ',FCPP,1,NSA)
            ELSE
               CALL FPGRACAB('FCP ',FCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CT ') THEN
            IF(NSB.eq.0) THEN
               CALL FPGRACA('FCT ',FCTH,2,NSA)
            ELSE
               CALL FPGRACAB('FCT ',FCTH2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL FPGRACA('FEP ',FEPP,1,NSA)
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL FPGRACA('FET ',FETH,2,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'R') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPGRARA('RN  ',RNT,NSA)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPGRARA('RI  ',RJT,NSA)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRARA('RW  ',RWT,NSA)
         ELSE IF(KID2.EQ.'PC ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRARA('RPC ',RPCT,NSA)
            ELSE
               CALL FPGRARAB('RPC ',RPCT2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPGRARA('RPW ',RPWT,NSA)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPGRARA('RPE ',RPET,NSA)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPGRARA('RLH ',RLHT,NSA)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPGRARA('RFW ',RFWT,NSA)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPGRARA('REC ',RECT,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRARA('RT  ',RTT,NSA)
         ELSE IF(KID2.EQ.'TB ') THEN
            CALL FPGRARA('RTB ',RTT_BULK,NSA)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPGRAR('RQ  ',RQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPGRAR('RE  ',RET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'T') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPGRATA('TN  ',PNT,NSA)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPGRATA('TI  ',PIT,NSA)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPGRATA('TW  ',PWT,NSA)
         ELSE IF(KID2.EQ.'PC ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRATA('TPC ',PPCT,NSA)
            ELSE
               CALL FPGRATAB('TPC ',PPCT2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPGRATA('TPW ',PPWT,NSA)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPGRATA('TPE ',PPET,NSA)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPGRATA('TLH ',PLHT,NSA)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPGRATA('TFW ',PFWT,NSA)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPGRATA('TEC ',PECT,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRATA('TT  ',PTT,NSA)
         ELSE IF(KID2.EQ.'TB ') THEN
            CALL FPGRATA('TTB ',PTT_BULK,NSA)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPGRAT('TQ  ',PQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPGRAT('TE  ',PET)
         ELSE IF(KID2.EQ.'QE ') THEN
            CALL FPGRAT('TQE ',Q_ENG)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL FPGRACA('DPP ',DPP ,1,NSA)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL FPGRACA('DPT ',DPT ,1,NSA)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL FPGRACA('DTP ',DTP ,2,NSA)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL FPGRACA('DTT ',DTT ,2,NSA)
         ELSE IF(KID2.EQ.'CPP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCPP',DCPP,1,NSA)
            ELSE
               CALL FPGRACAB('DCPP',DCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CPT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCPT',DCPT,1,NSA)
            ELSE
               CALL FPGRACAB('DCPT',DCPT2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCTP',DCTP,2,NSA)
            ELSE
               CALL FPGRACAB('DCTP',DCTP2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACA('DCTT',DCTT,2,NSA)
            ELSE
               CALL FPGRACAB('DCTT',DCTT2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL FPGRACA('DWPP',DWPP,1,NSA)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL FPGRACA('DWPT',DWPT,1,NSA)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL FPGRACA('DWTP',DWTP,2,NSA)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL FPGRACA('DWTT',DWTT,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACA('DRR ',DRR ,3,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'N') THEN
         READ(KID2,*,ERR=1,END=1) NSA1
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
      ELSE IF (KID1.EQ.'W') THEN
         IF(KID2.EQ.'P  ') THEN
            CALL FPGRACA('WP  ',WEIGHP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRACA('WT  ',WEIGHT,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRACA('WR  ',WEIGHR,3,NSA)
         ENDIF
      ELSE IF (KID1.EQ.'G') THEN
         READ(KID2,*,ERR=1,END=1) NGRAPH
         WRITE(6,'(A,I3)') '# NGRAPH is changed to',NGRAPH
      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000
      ELSE IF (KID1.EQ.'Q') THEN
         GO TO 9000
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
      END IF
!
      ELSE
         KID1=KID(2:2)
         CALL GUCPTL(KID1)
         KID2=KID(3:5)
         CALL GUCPTL(KID2(1:1))
         CALL GUCPTL(KID2(2:2))
         CALL GUCPTL(KID2(3:3))
!
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'P  ') THEN
            CALL FPGRACPA('FP  ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPGRACPA('FT  ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPGRACPA('FR  ',FRR,3,NSA)
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL FPGRACPA('FPP ',FPP,1,NSA)
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL FPGRACPA('FTH ',FTH,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACPA('FRR ',FRR,3,NSA)
         ELSE IF(KID2.EQ.'CP ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('FCP ',FCPP,1,NSA)
            ELSE
               CALL FPGRACPAB('FCP ',FCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CT ') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('FCT ',FCTH,2,NSA)
            ELSE
               CALL FPGRACPAB('FCT ',FCTH2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL FPGRACPA('FEP ',FEPP,1,NSA)
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL FPGRACPA('FET ',FETH,2,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID3'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL FPGRACPA('DPP ',DPP ,1,NSA)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL FPGRACPA('DPT ',DPT ,1,NSA)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL FPGRACPA('DTP ',DTP ,2,NSA)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL FPGRACPA('DTT ',DTT ,2,NSA)
         ELSE IF(KID2.EQ.'CPP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCPP',DCPP,1,NSA)
            ELSE
               CALL FPGRACPAB('DCPP',DCPP2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CPT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCPT',DCPT,1,NSA)
            ELSE
               CALL FPGRACPAB('DCPT',DCPT2,1,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTP') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCTP',DCTP,2,NSA)
            ELSE
               CALL FPGRACPAB('DCTP',DCTP2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'CTT') THEN
            IF(NSB.EQ.0) THEN
               CALL FPGRACPA('DCTT',DCTT,2,NSA)
            ELSE
               CALL FPGRACPAB('DCTT',DCTT2,2,NSA,NSB)
            ENDIF
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL FPGRACPA('DWPP',DWPP,1,NSA)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL FPGRACPA('DWPT',DWPT,1,NSA)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL FPGRACPA('DWTP',DWTP,2,NSA)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL FPGRACPA('DWTT',DWTT,2,NSA)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPGRACPA('DRR ',DRR ,3,NSA)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID3'
         ENDIF
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
      END IF
      END IF

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
      real(8),DIMENSION(NRMAX,NSAMAX,NTG2MAX):: FRA
      real(8),dimension(NRMAX,NTG2MAX):: TEMP
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
      real(8),DIMENSION(NRMAX,NSBMAX,NSAMAX,NTG2MAX):: FRAB
      real(8),dimension(NRMAX,NTG2MAX):: TEMP
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
      real(8),dimension(NRMAX,NTG2MAX)::FR
      real(4),dimension(NRMAX):: GX,GY
      CHARACTER(LEN=*),INTENT(IN):: STRING
      integer:: NT2, NR
      real(4):: GXMIN, GXMAX, GYMAX0, GYMIN0, GYMIN, GYMAX
      real(4):: GYMIN1, GYMAX1, GYSTEP, GXMIN1, GXMAX1, GXSTEP
      real(4):: GXORG
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
          GX(NR)=GUCLIP(RM(NR))
      ENDDO
      CALL GMNMX1(GX,1,NRMAX,1,GXMIN,GXMAX)
      IF(RGMIN.NE.0.0) GXMIN=GUCLIP(RGMIN)
      IF(RGMAX.NE.1.0) GXMIN=GUCLIP(RGMAX)

      GYMAX0=-1.E30
      GYMIN0= 1.E30
      DO NT2=1,NTG2
         DO NR=1,NRMAX
            GY(NR)=GUCLIP(FR(NR,NT2))
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
      CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUE(0.,0.0,0.0,2*GYSTEP,NGULEN(2*GYSTEP))
      DO NT2=1,NTG2
         DO NR=1,NRMAX
            GY(NR)=GUCLIP(FR(NR,NT2))
            GX(NR)=GUCLIP(RM(NR))
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
      real(8),DIMENSION(NSAMAX,NTG1MAX):: FTA
      real(8),dimension(NTG1MAX):: TEMP
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
      real(8),DIMENSION(NSBMAX,NSAMAX,NTG1MAX):: FTAB
      real(8),dimension(NTG1MAX):: TEMP
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
      real(8),DIMENSION(NTG1MAX):: FT
      real(4),DIMENSION(NTG1MAX):: GX,GY
      CHARACTER(LEN=*),INTENT(IN):: STRING
      integer:: NT1
      real(4):: GYMIN, GYMAX, GYMIN1, GYMAX1, GYSTEP
      real(4):: GXMIN, GXMAX, GXSTEP

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
         GX(NT1)=GUCLIP(PTG(NT1))
         GY(NT1)=GUCLIP(FT(NT1))
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
      CALL GVALUE(0.,GXSTEP*2,0.,0.,NGULEN(2*GXSTEP))
      CALL GVALUE(0.,0.,0.,GYSTEP*2,NGULEN(2*GYSTEP))
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
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX,NSAMAX):: FGA
      real(8),dimension(NTHMAX,NPMAX,NRMAX)::TEMP
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
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX,NSBMAX,NSAMAX):: FGAB
      real(8),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
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
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      real(4),dimension(NPMAX):: GX
      real(4),dimension(NPMAX,NRMAX):: GY
      CHARACTER(LEN=*),INTENT(IN):: STRING
      INTEGER,INTENT(IN):: NSA
      CHARACTER(LEN=80):: STRING1
!      real(8):: PGMAX, RGMAX, RGMIN
      real(4):: GXMAX, GXMINP, GYMIN, GYMAX, GYMINP, GYSTEP, GXSTEP, GXMAXP, GYMAXP
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
               GY(NP,NR)=GUCLIP(LOG10(FG(NTH,NP,NR)))
            ENDIF
         END DO
      END DO
      DO NP=1,NPMAX
         GX(NP)=GUCLIP(PM(NP,NS)**2)
      ENDDO
      GXMAX=GUCLIP(PMAX(NS)**2)
      IF(PGMAX.NE.0.0) GXMAX=GUCLIP(PGMAX**2)

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
      CALL GVALUE(0.,GXSTEP*2,0.0,0.0,NGULEN(2*GXSTEP))
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
      REAL(8),DIMENSION(NTHMAX,NPMAX,NRMAX,NSAMAX):: FGA
      REAL(8),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
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
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX,NSBMAX,NSAMAX):: FGAB
      real(8),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
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
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      real(4),dimension(NPMAX):: GX
      real(4),dimension(NRMAX):: GY
      real(4),dimension(NPMAX,NRMAX):: GZ
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      integer:: NR, NP, NTH, NSA, NS, NPGMAX, NPM
      real(4):: GX1, GX2, GY1, GY2

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
         GX(NP)=GUCLIP(PM(NP,NS)**2)
      ENDDO
   10 CONTINUE
      DO NR=1,NRMAX
         GY(NR)=GUCLIP(RM(NR))
      ENDDO
      DO NR=1,NRMAX
      DO NP=1,NPMAX
         IF(FG(NTH,NP,NR).LT.1.D-14) THEN
            GZ(NP,NR)=-14.0
         ELSE
            GZ(NP,NR)=GUCLIP(LOG10(FG(NTH,NP,NR)))
         ENDIF
      ENDDO
      ENDDO

      GX1=3.0
      GX2=20.0
      GY1=2.0
      GY2=17.0

      CALL PAGES
      CALL FPGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NPM,NPGMAX,NRMAX,  &
           GUCLIP(PGMAX**2),GUCLIP(RGMIN),GUCLIP(RGMAX),TRIM(STRING1))
      CALL PAGEE

 9000 RETURN
      END SUBROUTINE FPGRAPR

!     ***********************************************************

!           SUBPROGRAM FOR 2D PROFILE

!     ***********************************************************

      SUBROUTINE FPGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX, &
                        GXMAX1,GYMIN1,GYMAX1,STR)

      IMPLICIT NONE
      REAL(4),        INTENT(IN):: GX1,GX2,GY1,GY2,GXMAX1,GYMIN1,GYMAX1
      INTEGER(4),     INTENT(IN):: NXM,NXMAX,NYMAX
      REAL(4),DIMENSION(NXMAX),      INTENT(IN):: GX
      REAL(4),DIMENSION(NYMAX),      INTENT(IN):: GY
      REAL(4),DIMENSION(NXM,NYMAX),  INTENT(IN):: GZ
      CHARACTER(LEN=*),             INTENT(IN):: STR
      INTEGER(4) :: I, NGULEN
      REAL(4)    :: GOX, GOY, GOZ, GPHI, GRADIUS, GSTEPX, GSTEPY,   &
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
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,NGULEN(GSTEPX))
      CALL GVALUE3DY(GSYMIN,GSTEPY,1,NGULEN(GSTEPY))
      CALL GVALUE3DZ(GSZMIN,GSTEPZ,2,NGULEN(GSTEPZ))

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
      REAL(8),DIMENSION(:,:,:,:):: FGA
      REAL(8),dimension(NTHMAX+1,NPMAX+1,NRMAX+1):: TEMP
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
      REAL(8),DIMENSION(NTHMAX+1,NPMAX+1,NSAMAX):: FGA
      REAL(4),dimension(NTHMAX+1,NPMAX+1):: GF
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
            GF(NP,NTH)=GUCLIP(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            GF(NP,NTH)=GUCLIP(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            GF(NP,NTH)=GUCLIP(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.3) THEN
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            GF(NP,NTH)=GUCLIP(FGA(NTH,NP,NSA))
         ENDDO
         ENDDO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               IF(FGA(NTH,NP,NSA).LT.1.D-70) THEN
                  GF(NP,NTH)=-70.0
               ELSE
                  GF(NP,NTH)=GUCLIP(LOG10(ABS(FGA(NTH,NP,NSA))))
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
      REAL(8),DIMENSION(:,:,:,:):: FGA
!      REAL(8),dimension(NTHMAX+1,NPMAX+1,NRMAX+1):: TEMP
      REAL(8),dimension(NTHMAX,NPMAX,NRMAX):: TEMP
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
      real(8),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1,NSBMAX,NSAMAX):: FGAB
      real(8),dimension((NRMAX+1)*(NTHMAX+1)*(NPMAX+1)):: TEMP
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
      real(8),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1):: FG
      real(4),DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real(4),dimension(NPMAX+1):: GP
      real(4),dimension(NTHMAX+1):: GTH
      real(8),dimension(8,NPMAX+1,NTHMAX+1)::KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real(4):: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real(4):: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
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
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.3) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               IF(FG(NTH,NP,NR).LT.1.D-70) THEN
                  GF(NP,NTH)=-70.0
               ELSE
                  GF(NP,NTH)=GUCLIP(LOG10(ABS(FG(NTH,NP,NR))))
               ENDIF
            END DO
         END DO
      ENDIF
!
      CALL FPGRACX(STRING,GF,MODE,NSB)

      IF(NRMAX.GT.1) GOTO 1
!
 9000 RETURN
      END SUBROUTINE FPGRAC
!--------------------------------------------------
      SUBROUTINE FPGRAC_2(STRING,FG,MODE,NSB)
!       
      IMPLICIT NONE
!      real(8),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1):: FG
!      real(4),DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real(8),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      real(4),DIMENSION(NPMAX,NTHMAX):: GF
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real(4):: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real(4):: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
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
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
            END DO
         END DO
      ELSEIF(MODE.EQ.3) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               GF(NP,NTH)=GUCLIP(FG(NTH,NP,NR))
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
                  GF(NP,NTH)=GUCLIP(LOG10(ABS(FG(NTH,NP,NR))))
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
      real(4),DIMENSION(NPMAX,NTHMAX):: GF
      real(4),dimension(NPMAX):: GP
      real(4),dimension(NTHMAX):: GTH
      real(8),dimension(8,NPMAX,NTHMAX)::KA
!      real(4),DIMENSION(NPMAX+1,NTHMAX+1):: GF
!      real(4),dimension(NPMAX+1):: GP
!      real(4),dimension(NTHMAX+1):: GTH
!      real(8),dimension(8,NPMAX+1,NTHMAX+1)::KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real(4):: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real(4):: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
      integer:: NR, NP, NTH, NSB, NRG, MODE, LMODE
      integer:: NPG, NTHG, NPM, NTHM, NGLMAX, NGL
      integer:: I

      NPM=NPMAX+1
      NTHM=NTHMAX+1

      LMODE=MODE/4
!
      IF(MODE.EQ.1) THEN
         DO NP=1,NPMAX+1
            GP(NP)=GUCLIP(PG(NP,NSB))
         END DO
         NPG=NPMAX+1
      ELSE
         DO NP=1,NPMAX
            GP(NP)=GUCLIP(PM(NP,NSB))
         END DO
         NPG=NPMAX
      ENDIF

      IF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            GTH(NTH)=GUCLIP(THG(NTH))
         END DO
         NTHG=NTHMAX+1
      ELSE
         DO NTH=1,NTHMAX
            GTH(NTH)=GUCLIP(THM(NTH))
         END DO
         NTHG=NTHMAX
      ENDIF
!
      GPMAX=GUCLIP(PMAX(NSB))
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
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
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
!
      RETURN
      END SUBROUTINE FPGRACX_2
!--------------------------------------------------
      SUBROUTINE FPGRACX(STRING,GF,MODE,NSB)
!       
      IMPLICIT NONE
      real(4),DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real(4),dimension(NPMAX+1):: GP
      real(4),dimension(NTHMAX+1):: GTH
      real(8),dimension(8,NPMAX+1,NTHMAX+1)::KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER:: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real(4):: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real(4):: GPMIN1, GPMAX1, GPSTEP, GLIN, GFFMAX
      integer:: NR, NP, NTH, NSB, NRG, MODE, LMODE
      integer:: NPG, NTHG, NPM, NTHM, NRM, NM, NGLMAX, NGL
      integer:: I

      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1

      LMODE=MODE/4
!
      IF(MODE.EQ.1) THEN
         DO NP=1,NPMAX+1
            GP(NP)=GUCLIP(PG(NP,NSB))
         END DO
         NPG=NPMAX+1
      ELSE
         DO NP=1,NPMAX
            GP(NP)=GUCLIP(PM(NP,NSB))
         END DO
         NPG=NPMAX
      ENDIF

      IF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            GTH(NTH)=GUCLIP(THG(NTH))
         END DO
         NTHG=NTHMAX+1
      ELSE
         DO NTH=1,NTHMAX
            GTH(NTH)=GUCLIP(THM(NTH))
         END DO
         NTHG=NTHMAX
      ENDIF
!
      GPMAX=GUCLIP(PMAX(NSB))
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
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
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
!
      RETURN
      END SUBROUTINE FPGRACX

      SUBROUTINE SETRGBFP(F,RGB)
        IMPLICIT NONE
        REAL(4),INTENT(IN):: F
        REAL(4),DIMENSION(3),INTENT(OUT):: RGB
        INTEGER,PARAMETER:: NFMAX=8
        REAL(4),DIMENSION(3,NFMAX):: RGBC
        DATA RGBC/ 0.0,0.0,0.0, &
                   0.0,0.0,1.0, &
                   0.0,0.8,1.0, &
                   0.0,0.8,0.0, &
                   1.0,0.8,0.0, &
                   1.0,0.4,0.0, &
                   1.0,0.0,0.0, &
                   1.0,1.0,1.0/
        REAL(8):: GF,DF
        INTEGER(4):: IM
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
      real(8),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1,NSAMAX):: FGA
      real(8),dimension((NRMAX+1)*(NPMAX+1)*(NTHMAX+1)):: TEMP
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
      real(8),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1,NSBMAX,NSAMAX):: FGAB
      real(8),dimension((NRMAX+1)*(NPMAX+1)*(NTHMAX+1)):: TEMP
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
      real(8),DIMENSION((NRMAX+1)*(NPMAX+1)*(NTHMAX+1)):: FG
      real(4),DIMENSION(NPMAX+1,NTHMAX+1):: GF
      real(4),dimension(NPMAX+1):: GP
      real(4),dimension(NTHMAX+1):: GTH
      integer(4),dimension(8,NPMAX+1,NTHMAX+1):: KA
      CHARACTER(LEN=*),INTENT(IN):: STRING
      CHARACTER(LEN=80):: STRING1
      INTEGER,PARAMETER:: NGLM=30
      REAL(4):: ZL(NGLM),RGB(3,NGLM),WLN(NGLM)
      INTEGER(4):: ILN(NGLM)
      REAL(4):: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
      real(4):: GPMAX, GFMIN, GFMAX, GFMIN1, GFMAX1, GFSTEP
      real(4):: GPMIN1, GPMAX1, GPSTEP, GLIN
      integer:: NR, NTH, NP, NSA, MODE, NRG, LMODE, NTHG
      integer:: NPM, NTHM, NRM, NM, NGL, NGLMAX, I, NPG
      NPM=NPMAX+1
      NTHM=NTHMAX+1
      NRM=NRMAX+1


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
            GP(NP)=GUCLIP(PM(NP,NSA))
         END DO
         NPG=NPMAX
      ELSE
         DO NP=1,NPMAX+1
            GP(NP)=GUCLIP(PG(NP,NSA))
         END DO
         NPG=NPMAX+1
      ENDIF
!
      IF(MOD(MODE/2,2).EQ.0) THEN
         DO NTH=1,NTHMAX
            GTH(NTH)=GUCLIP(THM(NTH))
         END DO
         NTHG=NTHMAX
      ELSE
         DO NTH=1,NTHMAX+1
            GTH(NTH)=GUCLIP(THG(NTH))
         END DO
         NTHG=NTHMAX+1
      ENDIF
!
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               GF(NP,NTH)=GUCLIP(FG(NM)*PM(NP,NSA))
            END DO
         END DO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               GF(NP,NTH)=GUCLIP(FG(NM)*PG(NP,NSA))
            END DO
         END DO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               GF(NP,NTH)=GUCLIP(FG(NM)*PM(NP,NSA))
            END DO
         END DO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               IF(FG(NM)*PM(NP,NSA).LT.1.D-70) THEN
                  GF(NP,NTH)=-70.0
               ELSE
                  GF(NP,NTH)=GUCLIP(LOG10(ABS(FG(NM)*PM(NP,NSA))))
               ENDIF
            END DO
         END DO
      ENDIF
      GPMAX=GUCLIP(PMAX(NSA))

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
      CALL GVALUE(0.,GPSTEP*2,0.,GPSTEP*2,NGULEN(2*GPSTEP))
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

      END MODULE FPGOUT
