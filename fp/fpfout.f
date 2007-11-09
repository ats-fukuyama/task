C     $Id$
C
C ***********************************************************
C
C                    SELECTION OF FILE OUTPUT
C
C ***********************************************************
C
      SUBROUTINE FPFOUT
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION TEMP(NTHM,NPM,NRM)
      CHARACTER KID*4,KID1*1,KID2*3,KNAM*72
      LOGICAL LEX,LOP
C
      INQUIRE(22,OPENED=LOP)
      IF(.NOT.LOP) THEN
C 1001    WRITE(6,*) '# INPUT : FPFOUT FILE NAME : ',KNAMFO
C         READ(5,'(A72)',ERR=1001,END=9000) KNAM
C         IF(KNAM(1:2).NE.'/ ') KNAMFO=KNAM
C
         INQUIRE(FILE=KNAMFO,EXIST=LEX)
         IF(LEX) THEN
C            WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
C     &                 'ARE YOU SURE {Y/N} ?'
C            READ(5,'(A1)') KID
C            CALL GUCPTL(KID)
C            IF(KID.NE.'Y') GOTO 1001
            OPEN(22,FILE=KNAMFO,IOSTAT=IST,STATUS='OLD',ERR=1002,
     &           FORM='FORMATTED')
            WRITE(6,*) '# OLD FILE (',KNAMFO,') IS ASSIGNED FOR OUTPUT.'
            GOTO 1005
 1002       WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
            GOTO 9000
         ELSE
            OPEN(22,FILE=KNAMFO,IOSTAT=IST,STATUS='NEW',ERR=1003,
     &           FORM='FORMATTED')
            WRITE(6,*) '# NEW FILE (',KNAMFO,') IS CREATED FOR OUTPUT.'
            GOTO 1005
 1003       WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
            GOTO 9000
         ENDIF
      ENDIF
C
 1005 CONTINUE
C
    1 WRITE(6,*)'INPUT DATA TYPE : F/FX/FS1/FS2 1/2, N:newfile, X:exit,'
      WRITE(6,*)'             : D/DC/DW PP/PT/TP/TT/RR, F/FC/FE P/T/R'
      WRITE(6,*)'             : R/T N/I/W/PC/PW/PE/T/Q/E'
      READ(5,'(A4)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      KID2=KID(2:4)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2(1:1))
      CALL GUCPTL(KID2(2:2))
      CALL GUCPTL(KID2(3:3))
C
      IF (KID1.EQ.'F') THEN
         IF(KID2.EQ.'1  ') THEN
            CALL FPFOTP('F1  ',F)
         ELSE IF(KID2.EQ.'2  ') THEN
            CALL FPFOTC('F2  ',F,4)
         ELSE IF(KID2.EQ.'X2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=F(NTH,NP,NR)-F(NTHMAX+1-NTH,NP,NR)
            ENDDO
            ENDDO
            ENDDO
            CALL FPFOTC('FX2 ',TEMP,4)
         ELSE IF(KID2.EQ.'Y2  ') THEN
            DO NR=1,NRMAX
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               TEMP(NTH,NP,NR)=F(NTH,NP,NR)-F(NTHMAX+1-NTH,NP,NR)
            ENDDO
            ENDDO
            ENDDO
            CALL FPFOTC('FY2 ',TEMP,0)
         ELSE IF(KID2.EQ.'S11') THEN
            CALL FPFOTP('FS11',FS1)
         ELSE IF(KID2.EQ.'S12') THEN
            CALL FPFOTC('FS12',FS1,4)
         ELSE IF(KID2.EQ.'S21') THEN
            CALL FPFOTP('FS21',FS2)
         ELSE IF(KID2.EQ.'S22') THEN
            CALL FPFOTC('FS22',FS2,4)
         ELSE IF(KID2.EQ.'P  ') THEN
            CALL FPFOTC('FP  ',FPP,1)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPFOTC('FT  ',FTH,2)
         ELSE IF(KID2.EQ.'R  ') THEN
            CALL FPFOTC('FR  ',FRR,0)
         ELSE IF(KID2.EQ.'PP ') THEN
            CALL FPFOTC('FPP ',FPP,1)
         ELSE IF(KID2.EQ.'TH ') THEN
            CALL FPFOTC('FTH ',FTH,2)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPFOTC('FRR ',FRR,0)
         ELSE IF(KID2.EQ.'CP ') THEN
            CALL FPFOTC('FCP ',FCPP,1)
         ELSE IF(KID2.EQ.'CT ') THEN
            CALL FPFOTC('FCT ',FCTH,2)
         ELSE IF(KID2.EQ.'EP ') THEN
            CALL FPFOTC('FEP ',FEPP,1)
         ELSE IF(KID2.EQ.'ET ') THEN
            CALL FPFOTC('FET ',FETH,2)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'R') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPFOTR('RN  ',RNT)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPFOTR('RI  ',RJT)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPFOTR('RW  ',RWT)
         ELSE IF(KID2.EQ.'PC ') THEN
            CALL FPFOTR('RPC ',RPCT)
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPFOTR('RPW ',RPWT)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPFOTR('RPE ',RPET)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPFOTR('RLH ',RLHT)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPFOTR('RFW ',RFWT)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPFOTR('REC ',RECT)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPFOTR('RT  ',RTT)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPFOTR('RQ  ',RQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPFOTR('RE  ',RET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'T') THEN
         IF(KID2.EQ.'N  ') THEN
            CALL FPFOTT('TN  ',PNT)
         ELSE IF(KID2.EQ.'I  ') THEN
            CALL FPFOTT('TI  ',PIT)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPFOTT('TW  ',PWT)
         ELSE IF(KID2.EQ.'PC ') THEN
            CALL FPFOTT('TPC ',PPCT)
         ELSE IF(KID2.EQ.'PW ') THEN
            CALL FPFOTT('TPW ',PPWT)
         ELSE IF(KID2.EQ.'PE ') THEN
            CALL FPFOTT('TPE ',PPET)
         ELSE IF(KID2.EQ.'LH ') THEN
            CALL FPFOTT('TLH ',PLHT)
         ELSE IF(KID2.EQ.'FW ') THEN
            CALL FPFOTT('TFW ',PFWT)
         ELSE IF(KID2.EQ.'EC ') THEN
            CALL FPFOTT('TEC ',PECT)
         ELSE IF(KID2.EQ.'T  ') THEN
            CALL FPFOTT('TT  ',PTT)
         ELSE IF(KID2.EQ.'Q  ') THEN
            CALL FPFOTT('TQ  ',PQT)
         ELSE IF(KID2.EQ.'E  ') THEN
            CALL FPFOTT('TE  ',PET)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'D') THEN
         IF(KID2.EQ.'PP ') THEN
            CALL FPFOTC('DPP ',DPP ,1)
         ELSE IF(KID2.EQ.'PT ') THEN
            CALL FPFOTC('DPT ',DPT ,1)
         ELSE IF(KID2.EQ.'TP ') THEN
            CALL FPFOTC('DTP ',DTP ,2)
         ELSE IF(KID2.EQ.'TT ') THEN
            CALL FPFOTC('DTT ',DTT ,2)
         ELSE IF(KID2.EQ.'CPP') THEN
            CALL FPFOTC('DCPP',DCPP,1)
         ELSE IF(KID2.EQ.'CPT') THEN
            CALL FPFOTC('DCPT',DCPT,1)
         ELSE IF(KID2.EQ.'CTP') THEN
            CALL FPFOTC('DCTP',DCTP,2)
         ELSE IF(KID2.EQ.'CTT') THEN
            CALL FPFOTC('DCTT',DCTT,2)
         ELSE IF(KID2.EQ.'W  ') THEN
            CALL FPFOTC2('DW  ',DWPP,DWTT,0)
         ELSE IF(KID2.EQ.'WPP') THEN
            CALL FPFOTC('DWPP',DWPP,1)
         ELSE IF(KID2.EQ.'WPT') THEN
            CALL FPFOTC('DWPT',DWPT,1)
         ELSE IF(KID2.EQ.'WTP') THEN
            CALL FPFOTC('DWTP',DWTP,2)
         ELSE IF(KID2.EQ.'WTT') THEN
            CALL FPFOTC('DWTT',DWTT,2)
         ELSE IF(KID2.EQ.'RR ') THEN
            CALL FPFOTC('DRR ',DRR ,0)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2'
         ENDIF
      ELSE IF (KID1.EQ.'N') THEN
         CLOSE(22)
 2001    WRITE(6,*) '# INPUT : FPFOUT FILE NAME : ',KNAMFO
         READ(5,'(A72)',ERR=2001,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMFO=KNAM
C
         INQUIRE(FILE=KNAMFO,EXIST=LEX)
         IF(LEX) THEN
C            WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
C     &                 'ARE YOU SURE {Y/N} ?'
C            READ(5,'(A1)') KID
C            CALL GUCPTL(KID)
C            IF(KID.NE.'Y') GOTO 2001
            OPEN(22,FILE=KNAMFO,IOSTAT=IST,STATUS='OLD',ERR=2002,
     &           FORM='FORMATTED')
            WRITE(6,*) '# OLD FILE (',KNAMFO,') IS ASSIGNED FOR OUTPUT.'
            GOTO 2005
 2002       WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
            GOTO 2001
         ELSE
            OPEN(22,FILE=KNAMFO,IOSTAT=IST,STATUS='NEW',ERR=2003,
     &           FORM='FORMATTED')
            WRITE(6,*) '# NEW FILE (',KNAMFO,') IS CREATED FOR OUTPUT.'
            GOTO 2005
 2003       WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
            GOTO 2001
         ENDIF
C
 2005    CONTINUE
         GO TO 1
      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000
      ELSE IF (KID1.EQ.'Q') THEN
         GO TO 9000
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
      END IF
C
      GO TO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        RADIAL PROFILE
C
C ***********************************************************
C
      SUBROUTINE FPFOTR(STRING,FR)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FR(NRM,NTG1M)
      CHARACTER STRING*4
C
      CALL FPFOUX(STRING,NRMAX,NTG1,FR,NRM)
C
      RETURN
      END
C
C ***********************************************************
C
C                        TIME EVOLUTION
C
C ***********************************************************
C
      SUBROUTINE FPFOTT(STRING,FT)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION  FT(NTG2M)
      CHARACTER STRING*4
C
      CALL FPFOUX(STRING,NTG2,1,FT,NTG2M)
      RETURN
      END
C
C ***********************************************************
C
C                        Momentum Dependence
C
C ***********************************************************
C
      SUBROUTINE FPFOTP(STRING,FG)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION FG(NTHM,NPM,NRM),FL(NPM,NRM)
      CHARACTER STRING*4
C
    1 WRITE(6,*) '# INPUT NTH (1..',NTHMAX,') :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1
C
      DO NR=1,NRMAX
         DO NP=1,NPMAX
            FL(NP,NR)=FG(NTH,NP,NR)
         ENDDO
      ENDDO
C
      CALL FPFOUX(STRING,NPMAX,NRMAX,FL,NPM)
      GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPFOTC(STRING,FG,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*)
      DIMENSION FL(NPM,NTHM)
      CHARACTER STRING*4
C
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,*) '# INPUT NR (1..',NRG,') :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
C
      IF(MOD(MODE,2).EQ.0) THEN
         NPG=NPMAX
      ELSE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         NTHG=NTHMAX
      ELSE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=FG(NM)
            ENDDO
         ENDDO
      ENDIF         
C
      CALL FPFOUX(STRING,NPG,NTHG,FL,NPM)
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        C-GRAPHIC
C
C ***********************************************************
C
      SUBROUTINE FPFOTC2(STRING,FG,FH,MODE)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FG(*),FH(*)
      DIMENSION FL(NPM,NTHM)
      CHARACTER STRING*4
C
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,*) '# INPUT NR (1..',NRG,') :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF
C
      IF(MOD(MODE,2).EQ.0) THEN
         NPG=NPMAX
      ELSE
         NPG=NPMAX+1
      ENDIF
C
      IF(MOD(MODE/2,2).EQ.0) THEN
         NTHG=NTHMAX
      ELSE
         NTHG=NTHMAX+1
      ENDIF
C
      IF(MODE.EQ.0) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NMP1=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               NMP2=NPM*NTHM*(NR-1)+NTHM* NP   +NTH
               FGA=0.5D0*(FG(NMP1)+FG(NMP2))
               NMT1=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               NMT2=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH+1
               FHA=0.5D0*(FH(NMT1)+FH(NMT2))
               FL(NP,NTH)=SQRT(FGA**2+FHA**2)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX+1
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=SQRT(FG(NM)**2+FH(NM)**2)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NTH=1,NTHMAX+1
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=SQRT(FG(NM)**2+FH(NM)**2)
            ENDDO
         ENDDO
      ELSEIF(MODE.EQ.4) THEN
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
               NM=NPM*NTHM*(NR-1)+NTHM*(NP-1)+NTH
               FL(NP,NTH)=LOG10(ABS(SQRT(FG(NM)**2+FH(NM)**2)))
            ENDDO
         ENDDO
      ENDIF         
C
      CALL FPFOUX(STRING,NPG,NTHG,FL,NPM)
      IF(NRMAX.GT.1) GOTO 1
C
 9000 RETURN
      END
C
C ***********************************************************
C
C                        FILE OUTPUT EXEC
C
C ***********************************************************
C
      SUBROUTINE FPFOUX(STRING,N1MAX,N2MAX,FL,N1M)
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION FL(N1M,*)
      CHARACTER STRING*4
C
      WRITE(22,'(A4)') STRING
      WRITE(22,'(2I10)') N1MAX,N2MAX
      DO N2=1,N2MAX
         WRITE(22,'(1P5E15.7)') (FL(N1,N2),N1=1,N1MAX)
      ENDDO
      RETURN
      END
