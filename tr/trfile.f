C     $Id$
C
C     ***********************************************************
C
C           SAVE TRANSPORT DATA
C
C     ***********************************************************
C
      SUBROUTINE TRSAVE
C
      INCLUDE 'trcomm.h'
C
      CHARACTER*32 TRFNAM
C      CHARACTER*1 KID
      CHARACTER*3 K1,K2,K3,K4,K5,K6
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : SAVE FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFNAM
      IF(TRFNAM.EQ.'                                ') GOTO 900
      INQUIRE(FILE=TRFNAM,EXIST=LEX)
      IF(LEX) THEN
C         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
C     &              'ARE YOU SURE {Y/N} ?'
C         READ(5,'(A1)',ERR=1,END=900) KID
C         CALL GUCPTL(KID)
C         IF(KID.NE.'Y') GOTO 1
         OPEN(21,FILE=TRFNAM,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',TRFNAM,') IS ASSIGNED FOR OUTPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ELSE
         OPEN(21,FILE=TRFNAM,IOSTAT=IST,STATUS='NEW',ERR=20,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (',TRFNAM,') IS CREATED FOR OUTPUT.'
         GOTO 30
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ENDIF
C
   30 WRITE(21) NRMAX,DT,NGPST,TSST
      WRITE(21) RR,RA,RKAP,BB,RIPS,RIPE
      WRITE(21) PA,PZ,PN,PNS,PT,PTS
      WRITE(21) PNC,PNFE,PNNU,PNNUS
      WRITE(21) PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,PROFJ1,PROFJ2,
     &          ALP,AD0,CNC,CDW,MDLKAI,MDLETA,MDLAD,MDLAVK
      WRITE(21) TPRST,MDLST,MDLNF,IZERO
      WRITE(21) PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB
      WRITE(21) PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC
      WRITE(21) PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH
      WRITE(21) PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC
      WRITE(21) PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD
      WRITE(21) PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,PELPAT,MDLPEL
      WRITE(21) DR,FKAP,PNSS,T,TST,VSEC,WPPRE,TPRE
      WRITE(21) RG,RM,RN,RT,RW,BP,ANC,ANFE,ANNU
      WRITE(21) CHP,CK0,CKALFA,CKBETA,CKGUMA
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

      OPEN(16,ACCESS='APPEND',FILE=KFNLOG)
C
         CALL GUDATE(NDY,NDM,NDD,NTH1,NTM1,NTS1)
         WRITE(K1,'(I3)') 100+NDY
         WRITE(K2,'(I3)') 100+NDM
         WRITE(K3,'(I3)') 100+NDD
         WRITE(K4,'(I3)') 100+NTH1
         WRITE(K5,'(I3)') 100+NTM1
         WRITE(K6,'(I3)') 100+NTS1
         WRITE(16,1670) K1(2:3),K2(2:3),K3(2:3),K4(2:3),K5(2:3),K6(2:3),
     &                  TRFNAM,
     &                  RIPS,RIPE,PN(1),PN(2),BB,PICTOT,PLHTOT,PLHNPR
 1670    FORMAT(1H /
     &          1H ,'## DATE: ',
     &              A2,'-',A2,'-',A2,'  ',A2,':',A2,':',A2,' : ',
     &              '  FILE: ',A40/
     &          1H ,3X,'RIPS  =',1PD10.3,'  RIPE  =',1PD10.3,
     &               '  PNE   =',1PD10.3,'  PNI   =',1PD10.3/
     &          1H ,3X,'BB    =',1PD10.3,'  PICTOT=',1PD10.3,
     &               '  PLHTOT=',1PD10.3,'  PLHNPR=',1PD10.3)
         WRITE(16,1671) T,
     &                WPT,TAUE1,TAUE2,TAUEP,
     &                BETAP0,BETAPA,BETA0,BETAA
 1671    FORMAT(1H ,'# TIME : ',F7.3,' SEC'/
     &          1H ,3X,'WPT   =',1PD10.3,'  TAUE  =',1PD10.3,
     &               '  TAUED =',1PD10.3,'  TAUEP =',1PD10.3/
     &          1H ,3X,'BETAP0=',1PD10.3,'  BETAPA=',1PD10.3,
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)
C
         WRITE(16,1672) WST(1),TS0(1),TSAV(1),ANSAV(1),
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
 1672    FORMAT(1H ,3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3,
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/
     &          1H ,3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3,
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)
C
         WRITE(16,1673) AJT,VLOOP,ALI,Q0,
     &                AJOHT,AJNBT,AJRFT,AJBST
 1673    FORMAT(1H ,3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3,
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/
     &          1H ,3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3,
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)
C
         WRITE(16,1674) PINT,POHT,PNBT,
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4),
     &                POUT,PRLT,PCXT,PIET
 1674    FORMAT(1H ,3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3,
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/
     &          1H ,3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3,
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3)
C
      CLOSE(16)
C
  900 RETURN
      END
C
C     ***********************************************************
C
C           LOAD TRANSPORT DATA
C
C     ***********************************************************
C
      SUBROUTINE TRLOAD
C
      INCLUDE 'trcomm.h'
C
      CHARACTER*32 TRFNAM
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : LOAD FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFNAM
      IF(TRFNAM.EQ.'                                ') GOTO 900
      INQUIRE(FILE=TRFNAM,EXIST=LEX)
      IF(LEX) THEN
         OPEN(21,FILE=TRFNAM,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',TRFNAM,') IS ASSIGNED FOR INPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ELSE
         WRITE(6,*) 'XX FILE (',TRFNAM,') NOT FOUND'
         GOTO 1
      ENDIF
C
   30 READ(21) NRMAX,DT,NGPST,TSST
      READ(21) RR,RA,RKAP,BB,RIPS,RIPE
      READ(21) PA,PZ,PN,PNS,PT,PTS
      READ(21) PNC,PNFE,PNNU,PNNUS
      READ(21) PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,PROFJ1,PROFJ2,
     &         ALP,AD0,CNC,CDW,MDLKAI,MDLETA,MDLAD,MDLAVK
      READ(21) TPRST,MDLST,MDLNF,IZERO
      READ(21) PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB
      READ(21) PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC
      READ(21) PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH
      READ(21) PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC
      READ(21) PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD
      READ(21) PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,PELPAT,MDLPEL
      READ(21) DR,FKAP,PNSS,T,TST,VSEC,WPPRE,TPRE
      READ(21) RG,RM,RN,RT,RW,BP,ANC,ANFE,ANNU
      READ(21) CHP,CK0,CKALFA,CKBETA,CKGUMA
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
C
      NGR=0
      NGT=0
      NGST=0
      RIPS=RIPE
      GRG(1)=0.0
      DO 800 NR=1,NRMAX
         GRM(NR)  =GCLIP(RM(NR))
         GRG(NR+1)=GCLIP(RG(NR))
         QP(NR)  = FKAP*RG(NR)*BB/(RR*BP(NR))
  800 CONTINUE
      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
C
  900 RETURN
      END
C
C     ***********************************************************
C
C           SAVE GRAPHIC DATA
C
C     ***********************************************************
C
      SUBROUTINE TRGRSV
C
      INCLUDE 'trcomm.h'
C
      CHARACTER*32 TRFLNM
C      CHARACTER*1  KID
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : GRAPHIC SAVE FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFLNM
      IF(TRFLNM.EQ.'                                ') GOTO 900
      INQUIRE (FILE=TRFLNM,EXIST=LEX)
      IF(LEX) THEN
C         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
C     &              'ARE YOU SURE {Y/N}?'
C         READ(5,'(A1)',ERR=1,END=900) KID
C         CALL GUCPTL(KID)
C         IF(KID.NE.'Y') GOTO 1
         OPEN(22,FILE=TRFLNM,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',TRFLNM,') IS ASSIGNED FOR OUTPUT.'
         GOTO 30
   10    WRITE(6,*) '# XX OLD FILE OPEN ERROR : IOSTAT=',IST
         GOTO 1
      ELSE
         OPEN(22,FILE=TRFLNM,IOSTAT=IST,STATUS='NEW',ERR=20,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (',TRFLNM,') IS CREATED FOR OUTPUT.'
         GOTO 30
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT=',IST
         GOTO 1
      ENDIF
C
   30 WRITE(22) GVR,GRM,GRG,GTR,NGR
      WRITE(22) GVT,GT,NGT
      CLOSE(22)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
  900 RETURN
      END
C
C     ***********************************************************
C
C           LOAD GRAPHIC DATA
C
C     ***********************************************************
C
      SUBROUTINE TRGRLD
C
      INCLUDE 'trcomm.h'
C
      CHARACTER*32 TRFLNM
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : GRAPHIC LOAD FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFLNM
      IF(TRFLNM.EQ.'                                ') GOTO 900
      INQUIRE (FILE=TRFLNM,EXIST=LEX)
      IF(LEX) THEN
         OPEN(22,FILE=TRFLNM,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',TRFLNM,') IS ASSIGNED FOR INPUT.'
         GOTO 30
   10    WRITE(6,*) '# XX OLD FILE OPEN ERROR : IOSTAT=',IST
         GOTO 1
      ELSE
         WRITE(6,*) 'XX FILE (',TRFLNM,') NOT FOUND'
         GOTO 1
      ENDIF
   30 READ(22) GVR,GRM,GRG,GTR,NGR
      READ(22) GVT,GT,NGT
      CLOSE(22)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
C
  900 RETURN
      END
