C     $Id$
C
C     *************************** TASK.TR ***************************
C     *                                                             *
C     *  RADIAL TRANSPORT IN A TOKAMAK                              *
C     *                                                             *
C     *  PROGRAMMED BY                                              *
C     *                                                             *
C     *        A. FUKUYAMA                                          *
C     *        DEPARTMENT OF NUCLEAR ENGINEERING                    *
C     *        GRADUATE SCHOOL OF ENGINEERING                       *
C     *        KYOTO UNIVERSITY                                     *
C     *        KYOTO 606-8501, JAPAN                                *
C     *                                                             *
C     *  ACKNOWLEDGEMENTS FOR                                       *
C     *                                                             *
C     *        T. KASAI                                             *
C     *        M. ISHIZU                                            *
C     *        S. TAKATSUKA                                         *
C     *                                                             *
C     *  PROGRAM VERSION                                            *
C     *                                                             *
C     *        V1.00 : 89/04/26                                     *
C     *        V1.10 : 89/05/20                                     *
C     *        V1.20 : 90/06/12                                     *
C     *        V1.30 : 90/08/18                                     *
C     *        V1.40 : 91/07/24                                     *
C     *        V2.00 : 93/06/18                                     *
C     *        V2.01 : 93/07/28                                     *
C     *        V2.02 : 93/08/05                                     *
C     *        V2.03 : 93/08/17                                     *
C     *        V2.10 : 97/08/23                                     *
C     *        V2.20 : 98/09/21                                     *
C     *        V2.21 : 00/04/19                                     *
C     *                                                             *
C     ***************************************************************
C
      INCLUDE 'trcomm.h'
      PARAMETER(NURM=51,NUTM=61)
      COMMON /TRBCH1/ RAD(NURM),FRFHE(NURM,NUTM),FRFHI(NURM,NUTM)
      COMMON /TRBCH3/ FUT(NUTM),NUFMAX,NTXMAX
      CHARACTER KFILE*10
C
      CHARACTER KID*1,LINE*80
      CHARACTER KIG*5,K1*1,K2*1,K3*1,K4*1,K5*1,KK*3
      DATA INIT/0/,INQG/0/
C
C     ------ INITIALIZATION ------
C
      WRITE(6,600)
  600 FORMAT(' ')
      WRITE(6,601) NRM,NSM,NFM,NGM,NTM
  601 FORMAT(' ','# TASK/TR --- V2.30 : 04/03/21 ',
     &       '*** NRM, NSM, NFM, NGM, NTM  ***'/
     &       ' ',33X,I5,I5,I5,I5,I5)
      OPEN(33,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL GSOPEN
      CALL GUTIME(GTCPU1)
      CALL TRINIT
      CALL TRPARF
C
C     ------ SELECTION OF TASK TYPE ------
C
      IERR=0
      NTSMAX=NTMAX
    1 IF(INIT.EQ.0) THEN
         WRITE(6,*) '# INPUT : P,V,U/PARM  R/RUN  L/LOAD  ',
     &                        'D/DATA  H/HELP  Q/QUIT'
      ELSEIF(INIT.EQ.1) THEN
         WRITE(6,*) '# INPUT : G/GRAPH  ' 
         WRITE(6,*) '          P,V,U/PARM  R/RUN  L/LOAD  ',
     &                        'D/DATA  H/HELP  Q/QUIT'
      ELSE
         WRITE(6,*) '# INPUT : C/CONT  E/EQ  G/GRAPH  W/WRITE  S/SAVE ',
     &                        ' O/UFILEOUT '
         WRITE(6,*) '          P,V,U/PARM  R/RUN  L/LOAD  ',
     &                        'D/DATA  H/HELP  Q/QUIT'
      ENDIF
C
      CALL TRKLIN(LINE,KID,MODE,NTMAX,NTSMAX,IERR)
      IF(MODE.EQ.1.AND.IERR.EQ.2) GOTO 100
      IF(MODE.EQ.1.AND.IERR.EQ.3) GOTO 9000
      IF(MODE.NE.1) GOTO 1
C
  100 CONTINUE
      IF(KID.EQ.'P') THEN
         CALL TRPARM(KID)
         IF(KID.EQ.'Q') GOTO 9000
C
      ELSE IF(KID.EQ.'V') THEN
         CALL TRVIEW(0)
      ELSE IF(KID.EQ.'U') THEN
         CALL TRVIEW(1)
      ELSE IF(KID.EQ.'L') THEN
         CALL TRLOAD
         INIT=2
      ELSE IF(KID.EQ.'S'.AND.INIT.EQ.2) THEN
         CALL TRSAVE
      ELSE IF(KID.EQ.'R') THEN
         CALL TR_EQS_SELECT
         IF(MDLUF.EQ.1) THEN
            IF(INIT.EQ.2.AND.NT.NE.0) THEN
               NT=0
               NTMAX=NTSMAX
            ENDIF
            CALL TR_TIME_UFILE
         ELSEIF(MDLUF.EQ.2) THEN
            CALL TR_STEADY_UFILE
         ELSEIF(MDLUF.EQ.3) THEN
            KFILE='RFHE'
            CALL UFREAD2(KFILE,RAD,FUT,FRFHE,NUFMAX,NTXMAX,MDRFHE,IERR)
            KFILE='RFHI'
            CALL UFREAD2(KFILE,RAD,FUT,FRFHI,NUFMAX,NTXMAX,MDRFHI,IERR)
         ENDIF
         CALL TRPROF
         CALL TRLOOP
C
         INIT=2
         NTMOLD=NTMAX
      ELSE IF(KID.EQ.'E'.AND.INIT.EQ.2) THEN
         CALL TRCONV(L,IERR)
C   
      ELSE IF(KID.EQ.'C'.AND.INIT.EQ.2) THEN
         IF(MDLUF.EQ.1) THEN
            NT=NTMOLD
            NTMAX=NTSMAX+NTMOLD
         ENDIF
         CALL TRLOOP
         NTMOLD=NTMAX
C
      ELSE IF(KID.EQ.'G'.AND.INIT.GE.1) THEN
  101    WRITE(6,*) '# SELECT : R1-R9, T1-T9, G1-G5, P1-P5, Z1, Y1,',
     &                        ' A1-A2, E1-E9, D1-D32, M1-M3'
         WRITE(6,*) '           S/SAVE  L/LOAD  H/HELP  C/CLEAR  ',
     &              'I/INQ  X/EXIT'
         READ(5,'(A5)',ERR=101,END=9000) KIG
         K1=KIG(1:1)
         K2=KIG(2:2)
         K3=KIG(3:3)
         K4=KIG(4:4)
         K5=KIG(5:5)
         CALL GUCPTL(K1)
         CALL GUCPTL(K2)
         CALL GUCPTL(K3)
         CALL GUCPTL(K4)
         CALL GUCPTL(K5)
         KK=K3//K4//K5
C
         IF(K1.EQ.'C') THEN
            NGR=0
            NGT=0
         ELSEIF(K1.EQ.'I') THEN
            IF(INQG.EQ.0) THEN
               INQG=4      
               WRITE(6,*) '## GRAPHIC SCALE INQUIRE MODE : ON'
            ELSE
               INQG=0      
               WRITE(6,*) '## GRAPHIC SCALE INQUIRE MODE : OFF'
            ENDIF
         ELSEIF(K1.EQ.'S') THEN
            CALL TRGRSV
         ELSEIF(K1.EQ.'L') THEN
            CALL TRGRLD
         ELSEIF(K1.EQ.'R') THEN
            CALL TRGRR0(K2,INQG)
         ELSEIF(K1.EQ.'Y') THEN
            CALL TRGRY0(K2,INQG)
         ELSEIF(K1.EQ.'T') THEN
            CALL TRGRT0(K2,INQG)
         ELSEIF(K1.EQ.'Z') THEN
            CALL TRGRX0(K2,INQG)
         ELSEIF(K1.EQ.'G') THEN
            CALL TRGRG0(K2,INQG)
         ELSEIF(K1.EQ.'P') THEN
            CALL TRGRP0(K2,INQG)
         ELSEIF(K1.EQ.'A') THEN
            CALL TRGRA0(K2,INQG)
         ELSEIF(K1.EQ.'E') THEN
            CALL TRGRE0(K2,INQG)
         ELSEIF(K1.EQ.'D') THEN
            CALL TRGRD0(K2,KK,INQG)
         ELSEIF(K1.EQ.'H') THEN
            CALL TRHELP('G')
         ELSEIF(K1.EQ.'M') THEN
            CALL TRCOMP(K2,INQG)
         ELSEIF(K1.EQ.'X') THEN
            GOTO 1
         ELSE
            WRITE(6,*) 'UNSUPPORTED GRAPH ID'
         ENDIF
         GOTO 101
C
      ELSE IF(KID.EQ.'W'.AND.INIT.EQ.2) THEN
  102    WRITE(6,*) '# SELECT ',
     &              ': PRINT TYPE (1..9)  H/HELP  X/EXIT'
         READ(5,'(A1)',ERR=102,END=1) KID
         CALL GUCPTL(KID)
         IF(KID.EQ.'H') THEN
            CALL TRHELP('W')
         ELSEIF(KID.EQ.'X') THEN
            GOTO 1
         ELSE
            CALL TRPRNT(KID)
         ENDIF
         GOTO 102
C
      ELSE IF(KID.EQ.'D') THEN
         NFLMAX=0
    4    WRITE(6,*) '## HOW MANY DATA FILES ?'
         READ(5,*,END=1,ERR=4) NFLMAX
         DO 1000 NFL=1,NFLMAX
            CALL TRLOAD
            NGR=NFL-1
            NGT=NFL-1
            CALL TRCALC(IERR)
            CALL TRGLOB
            CALL TRATOT
            CALL TRATOG
 1000    CONTINUE
         INIT=1
C
      ELSE IF(KID.EQ.'O') THEN
         CALL TRXOUT
      ELSE IF(KID.EQ.'H') THEN
         CALL TRHELP('M')
      ELSE IF(KID.EQ.'Q') THEN
         GOTO 9000
C
      ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
         KID=' '
      ELSE
         WRITE(6,*) 'XX TRMAIN: UNKNOWN KID'
         KID=' '
      END IF
C
      GOTO 1
C
C     ------ END OF RUN ------
C
 9000 CALL GSCLOS
      CLOSE(33)
      CALL GUTIME(GTCPU2)
      WRITE(6,666) GTCPU2-GTCPU1
  666 FORMAT(' ','#      CPU TIME :   ',F8.3,' SEC ')
      STOP
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE TRKLIN(LINE,KID,MODE,NTMAX,NTSMAX,IERR)
C
      CHARACTER LINE*80,KID*1
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
      IF(IERR.NE.0) GOTO 4
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL TRPARL(LINE)
         IF(NTSMAX.NE.NTMAX) NTSMAX=NTMAX
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'P'.OR.
     &   KID.EQ.'V'.OR.
     &   KID.EQ.'U'.OR.
     &   KID.EQ.'L'.OR.
     &   KID.EQ.'S'.OR.
     &   KID.EQ.'R'.OR.
     &   KID.EQ.'E'.OR.
     &   KID.EQ.'C'.OR.
     &   KID.EQ.'G'.OR.
     &   KID.EQ.'W'.OR.
     &   KID.EQ.'D'.OR.
     &   KID.EQ.'O'.OR.
     &   KID.EQ.'H'.OR.
     &   KID.EQ.'Q') THEN
         MODE=1
         RETURN
      ENDIF
C
      KID=' '
      MODE=0
      RETURN
C
    2 WRITE(6,*) 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
    3 KID='Q'
      MODE=1
      RETURN
C
    4 KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'G') THEN
         MODE=1
         IERR=2
      ELSEIF (KID.EQ.'Q') THEN
         MODE=1
         IERR=3
      ENDIF
      RETURN
      END
