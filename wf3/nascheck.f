C     $Id$
C
      CALL WFRNAS
      STOP
      END
C
C     ******* INPUT NAS ELEMENT DATA *******
C
      SUBROUTINE WFRNAS
C
      LOGICAL LEX
      CHARACTER KDF*25,KL*36
      CHARACTER K0*8,K1*8,K2*8,K3*8,K4*8,K5*8,K6*8,K7*8,K8*8,K9*8
      CHARACTER KFNAMN*32
C
C      KFNAMN='Sample-01.NAS.txt'
C      KFNAMN='../data/Sample-en.NAS.txt'
C
      XNDMIN= 1.D5
      XNDMAX=-1.D5
      YNDMIN= 1.D5
      YNDMAX=-1.D5
      ZNDMIN= 1.D5
      ZNDMAX=-1.D5
      NDCOUNT=0
C
      WRITE(6,*) '# Input file name ?'
      READ(5,'(A32)') KFNAMN
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
  301    NNMAX=NNMAX+1
         READ(K1,*) IDND
         READ(K3,*) XND
         READ(K4,*) YND
         READ(K5,*) ZND
         XNDMIN=MIN(XNDMIN,XND)
         XNDMAX=MAX(XNDMAX,XND)
         YNDMIN=MIN(YNDMIN,YND)
         YNDMAX=MAX(YNDMAX,YND)
         ZNDMIN=MIN(ZNDMIN,ZND)
         ZNDMAX=MAX(ZNDMAX,ZND)
         IF(ZND.LE.1.D-6) NDCOUNT=NDCOUNT+1
         READ(26,'(10A8)',ERR=9100,END=9000) 
     &        K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
C         WRITE(6,'(A,2I8,1P3E12.4)') 
C     &   ' --- NODE ',NNMAX,IDND(NNMAX),
C     &                XND(NNMAX),YND(NNMAX),ZND(NNMAX)
         IF(K0.EQ.'GRID    ') GOTO 301
         WRITE(6,*) '--- NNMAX=',NNMAX
         WRITE(6,'(A,1P2E12.4)') '---    XNDMIN/MAX=',XNDMIN,XNDMAX
         WRITE(6,'(A,1P2E12.4)') '---    YNDMIN/MAX=',YNDMIN,YNDMAX
         WRITE(6,'(A,1P2E12.4)') '---    ZNDMIN/MAX=',ZNDMIN,ZNDMAX
         WRITE(6,'(A,I8)')       '---    NDCOUNT   =',NDCOUNT
         GOTO 201
C
      ELSEIF(K0.EQ.'CTETRA  ') THEN
         WRITE(6,*) 'CTETRA section'
  302    NEMAX=NEMAX+1
         READ(K1,*) IDELM
         READ(K2,*) KAELM
         READ(K3,*) NDELM
         READ(K4,*) NDELM
         READ(K5,*) NDELM
         READ(K6,*) NDELM
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
            NKMAX=NKMAX+1
            READ(26,'(10A8)',ERR=9100,END=9000) 
     &           K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
            IF(K0.EQ.'PSOLID  ') THEN
               READ(K1,*) IDTEMP
               IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
               READ(K2,*) NMKA
               WRITE(6,*) '--- PSOLID section: ',NKMAX,NMKA
               GOTO 200
            ENDIF
            GOTO 201
         ELSEIF(K1.EQ.'Material') THEN
            WRITE(6,*) 'Material section: ',KDF
            NMMAX=NMMAX+1
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
                  GOTO 501
               ENDIF
            ELSE
  502          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'MAT1    ') THEN
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
                  READ(K6,*) EPSDM
                  GOTO 502
               ELSEIF(K0.EQ.'MAT4    ') THEN
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDKA conflict'
                  IF(K2.NE.'        ') READ(K2,*) TANDLT
                  READ(K3,*) AMUDM
                  IF(K7.NE.'        ') READ(K7,*) SIGDM
                  IF(SIGDM.EQ.0.D0) THEN
                     SIGDM=2.D0*PI*RF*1.D6*EPS0
     &                            *EPSDM*TANDLT
                  ENDIF
                  GOTO 502
               ENDIF
               WRITE(6,'(A,I5,1P3E12.4)') 
     &              ' --- NM,EPS,AMU,SIG=',
     &              NMMAX,EPSDM,AMUDM,SIGDM
               GOTO 201
            ENDIF
         ELSEIF(K1.EQ.'Load Set') THEN
            WRITE(6,*) 'Load Set section: ',KDF
            NBMAX=NBMAX+1
C
            IF(KDF(1:8).EQ.'DotaiMen') THEN
               KABDY=1
               NBP=0
  601          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'PLOAD4  ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NENBP
                  READ(K8,*) NDNBP
                  GOTO 601
               ENDIF
               WRITE(6,'(A,2I8)') ' --- NB,NBPMAX=',NBMAX,NBP
               NBPMAX=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'KyushuMen') THEN
               KABDY=4
               NBP=0
  602          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'PLOAD4  ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NENBP
                  READ(K8,*) NDNBP
                  GOTO 602
               ENDIF
               WRITE(6,'(A,2I8)') ' --- NB,NBPMAX=',NBMAX,NBP
               NBPMAX=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Kuke_TE10') THEN
               KABDY=8
               READ(KDF(11:18),*) PWRWG
               NBP=0
  603          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP
                  READ(K5,*) EX1WG
                  READ(K6,*) EY1WG
                  READ(K7,*) EZ1WG
                  GOTO 603
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  GOTO 603
               ENDIF
               WRITE(6,'(A,2I8,1PE12.4)') 
     &              ' --- NB,NBPMAX,PWRWG=',NBMAX,NBP,PWRWG
               NBPMAX=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Enke_TE11') THEN
               READ(KDF(11:18),*) PWRWG
               READ(KDF(20:23),*) PHAWG
               NBP=0
  604          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP
                  READ(K5,*) EX1WG
                  READ(K6,*) EY1WG
                  READ(K7,*) EZ1WG
                  GOTO 604
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  READ(K5,*) EX2WG
                  READ(K6,*) EY2WG
                  READ(K7,*) EZ2WG
                  GOTO 604
               ENDIF
               WRITE(6,'(A,2I8,1P2E12.4)') 
     &              ' --- NB,NBPMAX,PWRWG,PHAWG=',NBMAX,NBP,
     &              PWRWG,PHAWG
               NBPMAX=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Enke_TM01') THEN
               KABDY=11
               READ(KDF(11:18),*) PWRWG
               NBP=0
  605          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP
                  READ(K5,*) EX1WG
                  READ(K6,*) EY1WG
                  READ(K7,*) EZ1WG
                  GOTO 605
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  GOTO 605
               ENDIF
               WRITE(6,'(A,2I8,1PE12.4)') 
     &              ' --- NB,NBPMAX,PWRWG=',NBMAX,NBP,PWRWG
               NBPMAX=NBP
               GOTO 201
C
            ELSEIF(KDF(1:9).EQ.'Dojik_TEM') THEN
               KABDY=12
               READ(KDF(11:18),*) PWRWG
               NBP=0
  606          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'FORCE   ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP
                  READ(K5,*) EX1WG
                  READ(K6,*) EY1WG
                  READ(K7,*) EZ1WG
                  GOTO 606
               ELSEIF(K0.EQ.'MOMENT  ') THEN
                  GOTO 606
               ENDIF
               WRITE(6,'(A,2I8,1PE12.4)') 
     &              ' --- NB,NBPMAX,PWRWG=',NBMAX,NBP,PWRWG
               NBPMAX=NBP
               GOTO 201
C
            ELSEIF(KDF(1:7).EQ.'Denastu') THEN
               KABDY=13
               NBP=0
  607          READ(26,'(10A8)',ERR=9100,END=9000) 
     &              K0,K1,K2,K3,K4,K5,K6,K7,K8,K9
               IF(K0.EQ.'TEMP    ') THEN
                  NBP=NBP+1
                  READ(K1,*) IDTEMP
                  IF(IDTEMP.NE.IDF) WRITE(6,*) '!! IDBDY conflict'
                  READ(K2,*) NDNBP
                  READ(K3,*) PHIBDY
                  GOTO 607
               ENDIF
               WRITE(6,'(A,2I8,1PE12.4)') 
     &              ' --- NB,NDNBP,PWRWG=',
     &              NBMAX,NDNBP,PHIBDY
               NBPMAX=NBP
               GOTO 201
            ENDIF
C
         ENDIF
C
      ELSEIF(K0.EQ.'ENDDATA ') THEN
         GOTO 8000
      ELSE
         WRITE(6,'(A)') K0
      ENDIF
      GOTO 200
C
 8000 CLOSE(25)
C
 9000 RETURN
C
 9100 WRITE(6,*) 'XX ERROR IN READING FILE: ',KFNAMN
      GOTO 9000
C
      END
