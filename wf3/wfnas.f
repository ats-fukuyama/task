C     $Id$
C
C     ********** NAS FILE PROCESSOR **********
C
      SUBROUTINE WFNAS
C
      USE libchar
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
      NSMAX=0
      NAMAX=0
C
    1 WRITE(6,*) '## INPUT: L/LOAD  G,N/DRAW  P,V/PARM  X/EXIT'
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL toupper(KID)
C
      IF(KID.EQ.'L') THEN
         CALL WFRNAS(IERR)
         IF(IERR.NE.0) GOTO 1
         WRITE(6,*) '--- WFMNAS start ---'
         CALL WFMNAS(IERR)
         IF(IERR.NE.0) GOTO 1
C
      ELSEIF(KID.EQ.'G') THEN
         CALL WFGNAS(0)
         CALL WFGNAS(1)
         CALL WFGNAS(2)
         CALL WFGNAS(3)
         CALL WFGNAS(4)
         CALL WFGNAS(5)
      ELSEIF(KID.EQ.'N') THEN
         CALL WFGNAS(1)
      ELSEIF(KID.EQ.'P') THEN
         CALL WFPARM(KID)
      ELSEIF(KID.EQ.'V') THEN
         CALL WFVIEW
      ELSEIF(KID.EQ.'S') THEN
         CALL WFWELM(1)
      ELSEIF(KID.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ******* INPUT NAS ELEMENT DATA *******
C
      SUBROUTINE WFRNAS(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      LOGICAL LEX
      CHARACTER KDF*25,KL*36
      CHARACTER K0*8,K1*8,K2*8,K3*8,K4*8,K5*8,K6*8,K7*8,K8*8,K9*8
      DIMENSION TANDLT(NMM),MODE_TANDL(NMM)
C
      IERR=0
C
      DO NB=1,NBM
         DO NBP=1,NBPM
            NENBP(NBP,NB)=0
            NDNBP(NBP,NB)=0
         ENDDO
      ENDDO
      DO NM=1,NMM
         MODE_TANDL(NM)=0
      ENDDO
C
      INQUIRE(FILE=KFNAMN,EXIST=LEX)
      IF(LEX) THEN
         OPEN(26,FILE=KFNAMN,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='FORMATTED')
         WRITE(6,*) '## FILE (',KFNAMN,') IS ASSIGNED FOR ELM INPUT'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX FILE (',KFNAMN,') IS NOT FOUND'
         GOTO 9000
      ENDIF
C
   30 CONTINUE
C
      READ(26,'(10A8)',ERR=9100,END=9000) 
     &     K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
      IF(K0(1:2).EQ.'ID') GOTO 101
      GOTO 201
C
  100 READ(26,'(10A8)',ERR=9100,END=9000) 
     &     K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
  101 IF(K0.NE.'+FEMAPC2') GOTO 100
      NNMAX=0
      NEMAX=0
      NMMAX=0
      NBMAX=0
      NKMAX=0
C
  200 READ(26,'(10A8)',ERR=9100,END=9000) 
     &     K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
  201 CONTINUE
C
      IF(K0.EQ.'FREQENCY') THEN
         WRITE(6,*) 'RF section'
         READ(K1,*) RF
         WRITE(6,'(A,1PE12.4)') ' --- RF=',RF
         RF=RF*1.D-6
C
      ELSEIF(K0.EQ.'UNIT_LEN') THEN
         WRITE(6,*) 'UNIT_LEN section'
         FACT_LEN=1.D0
         IF(K1(5:8).EQ.'km  ') FACT_LEN=1.D3
         IF(K1(5:8).EQ.'m   ') FACT_LEN=1.D0
         IF(K1(5:8).EQ.'cm  ') FACT_LEN=1.D-2
         IF(K1(5:8).EQ.'mm  ') FACT_LEN=1.D-3
         IF(K1(5:8).EQ.'mcrm') FACT_LEN=1.D-6
         IF(K1(5:8).EQ.'nm  ') FACT_LEN=1.D-9
         IF(K1(5:8).EQ.'A   ') FACT_LEN=1.D-10
         WRITE(6,'(A,1PE12.4)') ' --- FACT_LEN=',FACT_LEN
C
      ELSEIF(K0.EQ.'GRID    ') THEN
         WRITE(6,*) 'GRID section'
  301    IF(NNMAX.GE.NNM) THEN
            WRITE(6,*) 'XX WFRNAX: NNMAX .GT. NNM'
            IERR=1
            RETURN
         ENDIF
         NNMAX=NNMAX+1
         READ(K1,*) IDND(NNMAX)
         READ(K3,*) XND(NNMAX)
         READ(K4,*) YND(NNMAX)
         READ(K5,*) ZND(NNMAX)
         READ(26,'(10A8)',ERR=9100,END=9000) 
     &        K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
C         WRITE(6,'(A,2I8,1P3E12.4)') 
C     &   ' --- NODE ',NNMAX,IDND(NNMAX),
C     &                XND(NNMAX),YND(NNMAX),ZND(NNMAX)
         IF(K0.EQ.'GRID    ') GOTO 301
         WRITE(6,*) '--- NNMAX=',NNMAX
         GOTO 201
C
      ELSEIF(K0.EQ.'CTETRA  ') THEN
         WRITE(6,*) 'CTETRA section'
  302    IF(NEMAX.GE.NEM) THEN
            WRITE(6,*) 'XX WFRNAX: NEMAX .GT. NEM'
            IERR=1
            RETURN
         ENDIF
         NEMAX=NEMAX+1
         READ(K1,*) IDELM(NEMAX)
         READ(K2,*) KAELM(NEMAX)
         READ(K3,*) NDELM(1,NEMAX)
         READ(K4,*) NDELM(2,NEMAX)
         READ(K5,*) NDELM(3,NEMAX)
         READ(K6,*) NDELM(4,NEMAX)

C         WRITE(6,'(A,3I8,4I8)') 
C     &   ' --- ELM  ',NEMAX,IDELM(NEMAX),KAELM(NEMAX),
C     &        (NDELM(I,NEMAX),I=1,4)
         READ(26,'(10A8)',ERR=9100,END=9000) 
     &        K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
         IF(K0.EQ.'CTETRA  ') GOTO 302
         WRITE(6,*) '--- NEMAX=',NEMAX
         GOTO 201
C
      ELSEIF(K0.EQ.'$ FEMAP ') THEN
         KL=K2(2:8)//K3//K4//K5//K6(1:5)
         DO I=1,9
            IF(KL(I:I).EQ.' ') THEN
               ISP=I
               GOTO 401
            ENDIF
         ENDDO
  401    CONTINUE
         READ(KL(1:ISP-1),*) IDF
         KDF=KL(ISP+3:ISP+27)
C
         IF(K1.EQ.'Property') THEN
            WRITE(6,*) 'Property section: ',KDF
            IF(NKMAX.GE.NKM) THEN
               WRITE(6,*) 'XX WFRNAX: NKMAX .GT. NKM'
               IERR=1
               RETURN
            ENDIF
            NKMAX=NKMAX+1
            IDKA(NKMAX)=IDF
            KDKA(NKMAX)=KDF
            READ(26,'(10A8)',ERR=9100,END=9000) 
     &           K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
            IF(K0.EQ.'PSOLID  ') THEN
               READ(K1,*) IDTEMP
               IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
               READ(K2,*) NMKA(NKMAX)
               WRITE(6,*) '--- PSOLID section: ',NKMAX,NMKA(NKMAX)
               GOTO 200
            ENDIF
            GOTO 201
         ELSEIF(K1.EQ.'Material') THEN
            WRITE(6,*) 'Material section: ',KDF
            IF(NMMAX.GE.NMM) THEN
               WRITE(6,*) 'XX WFRNAX: NMMAX .GT. NMM'
               IERR=1
               RETURN
            ENDIF
            NMMAX=NMMAX+1
            IDMAT(NMMAX)=IDF
            KDMAT(NMMAX)=KDF
C
            IF(KDF(1:6).EQ.'plasma') THEN
  501          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'MAT1    ') THEN
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
                  READ(K2,*) PNL
                  READ(K3,*) PZCLL
                  WRITE(6,'(A,I5,1P2E12.4)') 
     &                 ' --- NM,PN,PZCL=',NMMAX,PNL,PZCLL
                  NSMAX=2
                  GOTO 501
               ENDIF
            ELSE
  502          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'MAT1    ') THEN
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
                  READ(K6,*) EPSDM(NMMAX)
                  GOTO 502
               ELSEIF(K0.EQ.'MAT4    ') THEN
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
                  IF(K2.NE.'        ') THEN
                     READ(K2,*) TANDLT(NMMAX)
                     MODE_TANDL(NMMAX)=1
                  ENDIF
                  READ(K3,*) AMUDM(NMMAX)
                  IF(K7.NE.'        ') THEN
                     READ(K7,*) SIGDM(NMMAX)
                     MODE_TANDL(NMMAX)=2
                  ENDIF
                  GOTO 502
               ENDIF
               IF(MODE_TANDL(NMMAX).EQ.1) THEN
                  WRITE(6,'(A,I5,1P3E12.4)') 
     &                 ' --- NM,EPS,AMU,SIG=',
     &                 NMMAX,EPSDM(NMMAX),AMUDM(NMMAX),SIGDM(NMMAX)
               ELSE
                  WRITE(6,'(A,I5,1P3E12.4)') 
     &                 ' --- NM,EPS,AMU,TDL=',
     &                 NMMAX,EPSDM(NMMAX),AMUDM(NMMAX),TANDLT(NMMAX)
               ENDIF
               GOTO 201
            ENDIF
         ELSEIF(K1.EQ.'Load Set') THEN
            WRITE(6,*) 'Load Set section: ',KDF
            IF(NBMAX.GE.NBM) THEN
               WRITE(6,*) 'XX WFRNAX: NBMAX .GT. NBM'
               IERR=1
               RETURN
            ENDIF
            NBMAX=NBMAX+1
            IDBDY(NBMAX)=IDF
            KDBDY(NBMAX)=KDF
C
            IF(KDF(1:8).EQ.'DotaiMen') THEN
               KABDY(NBMAX)=1
               NBP=0
  601          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'PLOAD4  ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NENBP(NBP,NBMAX)
                  READ(K8,*) NDNBP(NBP,NBMAX)
                  IF(NBMAX.EQ.1) WRITE(6,*) NBMAX,NBP,
     &                 NENBP(NBP,NBMAX),NDNBP(NBP,NBMAX)
                  GOTO 601
               ENDIF
               WRITE(6,'(A,2I5)') ' --- NB,NBPMAX=',NBMAX,NBP
               NBPMAX(NBMAX)=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'KyushuMen') THEN
               KABDY(NBMAX)=4
               NBP=0
  602          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'PLOAD4  ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NENBP(NBP,NBMAX)
                  READ(K8,*) NDNBP(NBP,NBMAX)
                  GOTO 602
               ENDIF
               WRITE(6,'(A,2I5)') ' --- NB,NBPMAX=',NBMAX,NBP
               NBPMAX(NBMAX)=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Kuke_TE10') THEN
               KABDY(NBMAX)=31
               READ(KDF(11:18),*) PWRWG(NBMAX)
               NBP=0
  603          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP(NBP,NBMAX)
                  READ(K5,*) EX1WG(NBMAX)
                  READ(K6,*) EY1WG(NBMAX)
                  READ(K7,*) EZ1WG(NBMAX)
                  GOTO 603
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  GOTO 603
               ENDIF
               WRITE(6,'(A,2I5,1PE12.4)') 
     &              ' --- NB,NBPMAX,PWRWG=',NBMAX,NBP,PWRWG(NBMAX)
               NBPMAX(NBMAX)=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Enke_TE11') THEN
               KABDY(NBMAX)=21
               READ(KDF(11:18),*) PWRWG(NBMAX)
               READ(KDF(20:23),*) PHAWG(NBMAX)
               NBP=0
  604          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP(NBP,NBMAX)
                  READ(K5,*) EX1WG(NBMAX)
                  READ(K6,*) EY1WG(NBMAX)
                  READ(K7,*) EZ1WG(NBMAX)
                  GOTO 604
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  READ(K5,*) EX2WG(NBMAX)
                  READ(K6,*) EY2WG(NBMAX)
                  READ(K7,*) EZ2WG(NBMAX)
                  GOTO 604
               ENDIF
               WRITE(6,'(A,2I5,1P2E12.4)') 
     &              ' --- NB,NBPMAX,PWRWG,PHAWG=',NBMAX,NBP,
     &              PWRWG(NBMAX),PHAWG(NBMAX)
               NBPMAX(NBMAX)=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Enke_TM01') THEN
               KABDY(NBMAX)=22
               READ(KDF(11:18),*) PWRWG(NBMAX)
               NBP=0
  605          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP(NBP,NBMAX)
                  READ(K5,*) EX1WG(NBMAX)
                  READ(K6,*) EY1WG(NBMAX)
                  READ(K7,*) EZ1WG(NBMAX)
                  GOTO 605
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  GOTO 605
               ENDIF
               WRITE(6,'(A,2I5,1PE12.4)') 
     &              ' --- NB,NBPMAX,PWRWG=',NBMAX,NBP,PWRWG(NBMAX)
               NBPMAX(NBMAX)=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Dojik_TEM') THEN
               KABDY(NBMAX)=11
               READ(KDF(11:18),*) PWRWG(NBMAX)
               NBP=0
  606          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP(NBP,NBMAX)
                  READ(K5,*) EX1WG(NBMAX)
                  READ(K6,*) EY1WG(NBMAX)
                  READ(K7,*) EZ1WG(NBMAX)
                  GOTO 606
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  GOTO 606
               ENDIF
               WRITE(6,'(A,2I5,1PE12.4)') 
     &              ' --- NB,NBPMAX,PWRWG=',NBMAX,NBP,PWRWG(NBMAX)
               NBPMAX(NBMAX)=NBP
               GOTO 201
C
            ELSEIF(KDF(1:7).EQ.'Denastu') THEN
               KABDY(NBMAX)=2
               NBP=0
  607          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'TEMP    ') THEN
                  IF(NBP.GE.NBPM) THEN
                     WRITE(6,*) 'XX WFRNAX: NBP .GT. NBPM AT NB=',NB
                     IERR=1
                     RETURN
                  ENDIF
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP(NBP,NBMAX)
                  READ(K3,*) PHIBDY(NBMAX)
                  GOTO 607
               ENDIF
               WRITE(6,'(A,2I5,1PE12.4)') 
     &              ' --- NB,NDNBP,PWRWG=',
     &              NBMAX,NDNBP(NBP,NBMAX),PHIBDY(NBMAX)
               NBPMAX(NBMAX)=NBP
               GOTO 201
            ENDIF
C
         ENDIF
C
      ELSEIF(K0.EQ.'ENDDATA ') THEN
         GOTO 9000
      ELSE
         WRITE(6,'(A)') K0
      ENDIF
      GOTO 200
C
 9000 CLOSE(25)
      DO NB=1,NBMAX
         WRITE(6,'(A,2I8,4X,A,2I8)') 
     &        'NB=',NB,IDBDY(NB),KDBDY(NB),KABDY(NB),NBPMAX(NB)
      ENDDO
      DO NK=1,NKMAX
         WRITE(6,'(A,2I8,4X,A)') 
     &        'NK=',NK,IDKA(NK),KDKA(NK)
      ENDDO
      DO NM=1,NMMAX
         WRITE(6,'(A,2I8,4X,A)') 
     &        'NM=',NM,IDMAT(NM),KDMAT(NM)
      ENDDO
C
      IF(NSMAX.GT.0) THEN
         PN(1)=PNL*1.D-14
         PN(2)=PNL*1.D-14*ABS(PZ(1))/PZ(2)
         PNS(1)=PNL*1.D-14
         PNS(2)=PNL*1.D-14*ABS(PZ(1))/PZ(2)
         PZCL(1)=PZCLL/(2.D0*PI*RF*1.D6)
         PZCL(2)=PZCLL/(2.D0*PI*RF*1.D6)
         WRITE(6,'(A,I8,1P2E12.4)') 'NSMAX,PNL,PZCLL=',NSMAX,PNL,PZCLL
         WRITE(6,'(A,1P4E12.4)') 'PN,PZCL=',PN(1),PN(2),PZCL(1),PZCL(2)
      ENDIF
      DO NM=1,NMMAX
         IF(MODE_TANDL(NM).EQ.1) THEN
            SIGDM(NM)=2.D0*PI*RF*1.D6*EPS0
     &               *EPSDM(NM)*TANDLT(NM)
         ENDIF
      ENDDO
C
      RETURN
C
 9100 WRITE(6,*) 'XX ERROR IN READING FILE: ',KFNAMN
      GOTO 9000
C
      END
C
C     ******* MODIFY NAS ELEMENT DATA *******
C
      SUBROUTINE WFMNAS(IERR)
C
      INCLUDE 'wfcomm.inc'
      EXTERNAL WFSRTS,WFSRTX
C
      IERR=0
C
C     ----- MULTIPLY UNIT LENGTH -----
C
      DO NN=1,NNMAX
         XND(NN)=FACT_LEN*XND(NN)
         YND(NN)=FACT_LEN*YND(NN)
         ZND(NN)=FACT_LEN*ZND(NN)
      ENDDO
C
C     ----- PREPARE TO REPLACE NODE INDEX BY NODE NUMBER -----
C
      IDNMAX=0
      DO NN=1,NNMAX
         IF(IDND(NN).GT.IDNMAX) IDNMAX=IDND(NN)
      ENDDO
      IDNMIN=IDNMAX
      DO NN=1,NNMAX
         IF(IDND(NN).LT.IDNMIN) IDNMIN=IDND(NN)
      ENDDO
      WRITE(6,*) 'IDNMIN,IDNMAX=',IDNMIN,IDNMAX
C
C     ----- PREPARE TO REPLACE ELM INDEX BY ELM NUMBER -----
C
      IDEMAX=0
      DO NE=1,NEMAX
         IF(IDELM(NE).GT.IDEMAX) IDEMAX=IDELM(NE)
      ENDDO
      IDEMIN=IDEMAX
      DO NE=1,NEMAX
         IF(IDELM(NE).LT.IDEMIN) IDEMIN=IDELM(NE)
      ENDDO
      WRITE(6,*) 'IDEMIN,IDEMAX=',IDEMIN,IDEMAX
C
C     ----- REPLACE ELM INDEX BY ELM NUMBER -----
C
      DO NE=1,NEMAX
         DO I=1,4
            ID=NDELM(I,NE)
            CALL WFIDNN(ID,NN)
            IF(IDND(NN).NE.ID) THEN
               WRITE(6,*) 'XX NN,IDND,ID=',NN,IDND(NN),ID
            ENDIF
            NDELM(I,NE)=NN
         ENDDO
      ENDDO
C
C     ----- Modify element numbers in boundary block -----
C
      DO NB=1,NBMAX
         DO NBP=1,NBPMAX(NB)
            ID=NENBP(NBP,NB)
            IF(ID.NE.0) THEN
               CALL WFIDNE(ID,NE)
               IF(IDELM(NE).NE.ID) THEN
                  WRITE(6,*) 'XX NE,IDELM,ID=',NE,IDELM(NE),ID
               ENDIF
               NENBP(NBP,NB)=NE
            ENDIF
C
            ID=NDNBP(NBP,NB)
            IF(ID.NE.0) THEN
               CALL WFIDNN(ID,NN)
               IF(IDND(NN).NE.ID) THEN
                  WRITE(6,*) 'XX NN,IDND,ID=',NN,IDND(NN),ID
               ENDIF
               NDNBP(NBP,NB)=NN
            ENDIF
         ENDDO
      ENDDO
C
      DO NE=1,NEMAX
         CALL WFCHKV(NE,V)
         IF(V.LT.0.D0) THEN
            ND=NDELM(3,NE)
            NDELM(3,NE)=NDELM(2,NE)
            NDELM(2,NE)=ND
         ENDIF
      ENDDO
C
C     ----- Sort elements by SINDEX -----
C
      WRITE(6,*) '--- WFINDX start ---'
      CALL WFINDX
      WRITE(6,*) '--- WFFEPI start ---'
      CALL WFFEPI
C
C     ----- Modify element numbers in boundary block -----
C
      DO NB=1,NBMAX
         DO NBP=1,NBPMAX(NB)
            ID=NENBP(NBP,NB)
            IF(ID.NE.0) THEN
               NENBP(NBP,NB)=IWELM(ID)
            ENDIF
         ENDDO
      ENDDO
C
C     ----- Adjust material number -----
C
      DO NK1=1,NKMAX
         ID=IDKA(NK1)
         NM=NMKA(NK1)
         DO N=1,NMMAX
            IF(IDMAT(N).EQ.NM) THEN
               NM=N
               GOTO 1000
            ENDIF
         ENDDO
         WRITE(6,*) 'XX UNDEFINED NMKA: NK1,NMKA=',NK1,NM
         STOP
C
 1000    IF(KDMAT(NM)(1:6).EQ.'plasma') THEN
            NM=0
         ENDIF
         NMKA(NK1)=NM
         WRITE(6,*) 'NK,ID,NM=',NK1,ID,NM
      ENDDO
C
C     ----- Modify material number of element -----
C
      DO NE=1,NEMAX
         KA=KAELM(NE)
         DO NK1=1,NKMAX
            IF(IDKA(NK1).EQ.KA) GOTO 1100
         ENDDO
         WRITE(6,*) 'XX UNDEFINED KAELM: NE,KAELM=',NE,KA
         STOP
C
 1100    CONTINUE
      ENDDO
C
C     ----- Set up WG data -----
C
      DO NB=1,NBMAX
         IF(KABDY(NB).GE.10) THEN
            XC=0.D0
            YC=0.D0
            ZC=0.D0
            DO NBP=1,NBPMAX(NB)
               NN=NDNBP(NBP,NB)
               XC=XC+XND(NN)
               YC=YC+YND(NN)
               ZC=ZC+ZND(NN)
            ENDDO
            XC=XC/NBPMAX(NB)
            YC=YC/NBPMAX(NB)
            ZC=ZC/NBPMAX(NB)
            IF(KABDY(NB).EQ.21) THEN
CCCC
CCCC Assuming second vector is always given
CCCC
               ANX1=EX1WG(NB)
               ANY1=EY1WG(NB)
               ANZ1=EZ1WG(NB)
               AN1=SQRT(ANX1*ANX1+ANY1*ANY1+ANZ1*ANZ1)
               ANX1=ANX1/AN1
               ANY1=ANY1/AN1
               ANZ1=ANZ1/AN1
               ANX2=EX2WG(NB)
               ANY2=EY2WG(NB)
               ANZ2=EZ2WG(NB)
               AN2=SQRT(ANX2*ANX2+ANY2*ANY2+ANZ2*ANZ2)
               ANX2=ANX2/AN2
               ANY2=ANY2/AN2
               ANZ2=ANZ2/AN2
               ANX3=ANY1*ANZ2-ANZ1*ANY2
               ANY3=ANZ1*ANX2-ANX1*ANZ2
               ANZ3=ANX1*ANY2-ANY1*ANX2
               SIZE1=0.D0
               SIZE2=0.D0
               DO NBP=1,NBPMAX(NB)
                  NN=NDNBP(NBP,NB)
                  SIZE1L=(XND(NN)-XC)*ANX1
     &                  +(YND(NN)-YC)*ANY1
     &                  +(ZND(NN)-ZC)*ANZ1
                  IF(SIZE1L.GT.SIZE1) SIZE1=SIZE1L
                  SIZE2L=(XND(NN)-XC)*ANX2
     &                  +(YND(NN)-YC)*ANY2
     &                  +(ZND(NN)-ZC)*ANZ2
                  IF(SIZE2L.GT.SIZE2) SIZE2=SIZE2L
               ENDDO
            ELSE
CCCC
CCCC Assuming first vector is always given
CCCC
               ANX1=EX1WG(NB)
               ANY1=EY1WG(NB)
               ANZ1=EZ1WG(NB)
               AN1=SQRT(ANX1*ANX1+ANY1*ANY1+ANZ1*ANZ1)
               ANX1=ANX1/AN1
               ANY1=ANY1/AN1
               ANZ1=ANZ1/AN1
               NN1=NDNBP(1,NB)
               NN2=NDNBP(NBPMAX(NB)/3,NB)
               NN3=NDNBP(NBPMAX(NB),NB)
               ANX3=(YND(NN1)-YND(NN2))*(ZND(NN3)-ZND(NN2))
     &             -(ZND(NN1)-ZND(NN2))*(YND(NN3)-YND(NN2))
               ANY3=(ZND(NN1)-ZND(NN2))*(XND(NN3)-XND(NN2))
     &             -(XND(NN1)-XND(NN2))*(ZND(NN3)-ZND(NN2))
               ANZ3=(XND(NN1)-XND(NN2))*(YND(NN3)-YND(NN2))
     &             -(YND(NN1)-YND(NN2))*(XND(NN3)-XND(NN2))
               AN3=SQRT(ANX3*ANX3+ANY3*ANY3+ANZ3*ANZ3)
               ANX3=ANX3/AN3
               ANY3=ANY3/AN3
               ANZ3=ANZ3/AN3
CCCCCC
               IF(ANZ3.GT.0.D0) THEN
                  ANX3=-ANX3
                  ANY3=-ANY3
                  ANZ3=-ANZ3
               ENDIF
CCCCCC
               ANX2=ANY3*ANZ1-ANZ3*ANY1
               ANY2=ANZ3*ANX1-ANX3*ANZ1
               ANZ2=ANX3*ANY1-ANY3*ANX1
               SIZE1=0.D0
               SIZE2=1.D16
               DO NBP=1,NBPMAX(NB)
                  NN=NDNBP(NBP,NB)
                  SIZEL=(XND(NN)-XC)**2
     &                 +(YND(NN)-YC)**2
     &                 +(ZND(NN)-ZC)**2
                  IF(SIZEL.GT.SIZE1) SIZE1=SIZEL
                  IF(SIZEL.LT.SIZE2) SIZE2=SIZEL
               ENDDO
               SIZE1=SQRT(SIZE1)
               SIZE2=SQRT(SIZE2)
            ENDIF
            WRITE(6,'(A,1P3E12.4)') ' --- XC,YC,ZC=',XC,YC,ZC
            WRITE(6,'(A,1P3E12.4)') ' --- ANX/Y/Z1=',ANX1,ANY1,ANZ1
            WRITE(6,'(A,1P3E12.4)') ' --- ANX/Y/Z2=',ANX2,ANY2,ANZ2
            WRITE(6,'(A,1P3E12.4)') ' --- ANX/Y/Z3=',ANX3,ANY3,ANZ3
            WRITE(6,'(A,1P2E12.4)') ' --- SIZE1/2 =',SIZE1,SIZE2
C
            PWRBDY(NB)=PWRWG(NB)
            PHABDY(NB)=PHAWG(NB)
            XPBDY(NB)=XC
            YPBDY(NB)=YC
            ZPBDY(NB)=ZC
            XNBDY(1,NB)=ANX1
            YNBDY(1,NB)=ANY1
            ZNBDY(1,NB)=ANZ1
            XNBDY(2,NB)=ANX2
            YNBDY(2,NB)=ANY2
            ZNBDY(2,NB)=ANZ2
            XNBDY(3,NB)=ANX3
            YNBDY(3,NB)=ANY3
            ZNBDY(3,NB)=ANZ3
            SZBDY(1,NB)=SIZE1
            SZBDY(2,NB)=SIZE2
            NMBDY(NB)=1
         ENDIF
      ENDDO
C      
      RETURN
      END
C
C     ******* RETURN NN FOR IDND *******
C
      SUBROUTINE WFIDNN(ID,NN)
C
      INCLUDE 'wfcomm.inc'
C
      NN=INT(NNMAX*DBLE(ID-IDNMIN)/DBLE(IDNMAX-IDNMIN+1))+1
    1 ID0=IDND(NN)
      IF(ID0.GT.ID) THEN
         NN=NN-1
         IF(NN.LE.0) GOTO 8000
         GOTO 1
      ELSEIF(ID0.LT.ID) THEN
         NN=NN+1
         IF(NN.GT.NNMAX) GOTO 8000
         GOTO 1
      ENDIF
      GOTO 9001
C
 8000 DO NNDO=1,NNMAX
         IF(IDND(NNDO).EQ.ID) THEN
            NN=NNDO
            GOTO 9002
         ENDIF
      ENDDO
      WRITE(6,*) 'XX WFIDNN ERROR: ID=',ID
      STOP
C
 9001 IF(IDND(NN).NE.ID) WRITE(6,*) '9001: ID0,ID=',ID0,ID
      RETURN
 9002 IF(IDND(NN).NE.ID) WRITE(6,*) '9002'
      RETURN
      END
C
C     ******* RETURN NE FOR IDELM *******
C
      SUBROUTINE WFIDNE(ID,NE)
C
      INCLUDE 'wfcomm.inc'
C
      NE=ID-IDEMIN+1
      IF(IDELM(NE).EQ.ID) RETURN
C
    1 ID0=IDELM(NE)
      IF(ID0.GT.ID) THEN
         NE=NE-1
         IF(NE.LE.0) GOTO 8000
         GOTO 1
      ELSEIF(ID0.LT.ID) THEN
         NE=NE+1
         IF(NE.GT.NEMAX) GOTO 8000
      ENDIF
      RETURN
C
 8000 DO NEDO=1,NEMAX
         IF(IDELM(NEDO).EQ.ID) THEN
            NE=NEDO
            RETURN
         ENDIF
      ENDDO
      WRITE(6,*) 'XX WFIDNE ERROR: ID=',ID
      STOP
      END
C
C     ******* A,B,C,D,V  CALCULATION *******
C
      SUBROUTINE WFCHKV(NE,V)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION XE(4),YE(4),ZE(4)
C
      CALL WFNODE(NE,XE,YE,ZE)
C
      V=(XE(2)-XE(1))*(YE(3)-YE(1))*(ZE(4)-ZE(1))
     & -(XE(2)-XE(1))*(YE(4)-YE(1))*(ZE(3)-ZE(1))
     & +(XE(3)-XE(1))*(YE(4)-YE(1))*(ZE(2)-ZE(1))
     & -(XE(3)-XE(1))*(YE(2)-YE(1))*(ZE(4)-ZE(1))
     & +(XE(4)-XE(1))*(YE(2)-YE(1))*(ZE(3)-ZE(1))
     & -(XE(4)-XE(1))*(YE(3)-YE(1))*(ZE(2)-ZE(1))
      V=V/6.D0
      RETURN
      END
