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
      INCLUDE 'trcomm.inc'
C
      CHARACTER TRFNAM*32
C      CHARACTER KID*1
      CHARACTER K1*3,K2*3,K3*3,K4*3,K5*3,K6*3
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
     &          ALP,AD0,CNP,CNH,CDP,CDH,CDW,MDLKAI,MDLETA,MDLAD,MDLAVK
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

C      OPEN(16,POSITION='APPEND',FILE=KFNLOG)
C      OPEN(16,ACCESS='APPEND',FILE=KFNLOG)
      OPEN(16,ACCESS='SEQUENTIAL',FILE=KFNLOG)
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
 1670    FORMAT(' '/
     &          ' ','## DATE: ',
     &              A2,'-',A2,'-',A2,'  ',A2,':',A2,':',A2,' : ',
     &              '  FILE: ',A40/
     &          ' ',3X,'RIPS  =',1PD10.3,'  RIPE  =',1PD10.3,
     &               '  PNE   =',1PD10.3,'  PNI   =',1PD10.3/
     &          ' ',3X,'BB    =',1PD10.3,'  PICTOT=',1PD10.3,
     &               '  PLHTOT=',1PD10.3,'  PLHNPR=',1PD10.3)
         WRITE(16,1671) T,
     &                WPT,TAUE1,TAUE2,TAUEP,
     &                BETAP0,BETAPA,BETA0,BETAA
 1671    FORMAT(' ','# TIME : ',F7.3,' SEC'/
     &          ' ',3X,'WPT   =',1PD10.3,'  TAUE  =',1PD10.3,
     &               '  TAUED =',1PD10.3,'  TAUEP =',1PD10.3/
     &          ' ',3X,'BETAP0=',1PD10.3,'  BETAPA=',1PD10.3,
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)
C
         WRITE(16,1672) WST(1),TS0(1),TSAV(1),ANSAV(1),
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
 1672    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3,
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3,
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)
C
         WRITE(16,1673) AJT,VLOOP,ALI,Q0,
     &                AJOHT,AJNBT,AJRFT,AJBST
 1673    FORMAT(' ',3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3,
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3,
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)
C
         WRITE(16,1674) PINT,POHT,PNBT,
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4),
     &                POUT,PRLT,PCXT,PIET
 1674    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3,
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3,
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
      INCLUDE 'trcomm.inc'
C
      CHARACTER TRFNAM*32
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
     &         ALP,AD0,CNP,CNH,CDP,CDH,CDW,MDLKAI,MDLETA,MDLAD,MDLAVK
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
      DO NR=1,NRMAX
         GRM(NR)  =GUCLIP(RM(NR))
         GRG(NR+1)=GUCLIP(RG(NR))
         QP(NR)   =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
      ENDDO
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
      INCLUDE 'trcomm.inc'
C
      CHARACTER TRFLNM*32
C      CHARACTER KID*1
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
      INCLUDE 'trcomm.inc'
C
      CHARACTER TRFLNM*32
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
C
C     ***********************************************************
C
C           CHECKING WHETHER IMPURITY EXISTS
C
C     ***********************************************************
C
      SUBROUTINE CHECK_IMPURITY(MDSLCT)
C
      INCLUDE 'trcomm.inc'
      CHARACTER KDIRX*80
      CHARACTER KDIRR2*80
      CHARACTER KFILE*80,KFID*20
      COMMON /TRUFC1/ KDIRX
      COMMON /TRUFC3/ NM2CHK
      LOGICAL LEX
C
      CALL KTRIM(KUFDEV,IKNDEV)
      CALL KTRIM(KUFDCG,IKNDCG)
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'2d'//KUFDCG(1:IKNDCG)//'.'
      CALL KTRIM(KDIRR2,KL2)
C
      IF(MDNI.LT.0.OR.MDNI.GT.3) MDNI=0
      MDSLCT=0
C
      KFID='NM2'
      KFILE=KDIRR2(1:KL2)//KFID
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF((LEX.EQV..TRUE.).AND.
     &   (KUFDEV.EQ.'tftr'.OR.KUFDEV.EQ.'d3d')) THEN
         MDSLCT=MDSLCT+1
         NM2CHK=1
      ELSE
         NM2CHK=0
         KFID='NM1'
         KFILE=KDIRR2(1:KL2)//KFID
         INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
         IF(LEX) MDSLCT=MDSLCT+1
      ENDIF
      KFID='ZEFFR'
      KFILE=KDIRR2(1:KL2)//KFID
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) MDSLCT=MDSLCT+2
      KFID='NIMP'
      KFILE=KDIRR2(1:KL2)//KFID
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) MDSLCT=MDSLCT+4
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           CONTROL ROUTINE FOR UFILE READING
C
C     ***********************************************************
C
      SUBROUTINE TR_UFILE_CONTROL(NSW)
C
      INCLUDE 'trcomm.inc'
C
      NRAMAX=INT(RHOA*NRMAX)
      DR = 1.D0/DBLE(NRMAX)
C
      IF(NSW.NE.0) THEN
         CALL CHECK_IMPURITY(MDSLCT)
C
         IF(MDSLCT.EQ.0) THEN
            MDNI=0
            NEQL=0
            DO NEQ=1,NEQMAX
               IF(NST(NEQ).NE.0) NEQL=NEQL+1
            ENDDO
            IF(NEQL.EQ.1) THEN
               NSMAX=1
            ELSE
               NSMAX=2
            ENDIF
         ELSEIF(MDSLCT.EQ.1) THEN
            IF(MDNI.EQ.2.OR.MDNI.EQ.3) MDNI=1
         ELSEIF(MDSLCT.EQ.2) THEN
            IF(MDNI.EQ.1.OR.MDNI.EQ.3) MDNI=2
         ELSEIF(MDSLCT.EQ.3) THEN
            IF(MDNI.EQ.3) MDNI=1
         ELSEIF(MDSLCT.EQ.4) THEN
            IF(MDNI.EQ.1.OR.MDNI.EQ.2) MDNI=3
         ELSEIF(MDSLCT.EQ.5) THEN
            IF(MDNI.EQ.2) MDNI=1
         ELSEIF(MDSLCT.EQ.6) THEN
            IF(MDNI.EQ.1) MDNI=3
         ENDIF
C
         IF(MDLFLX.NE.0) THEN
            CNP=0.D0
            CDP=0.D0
            AD0=0.D0
         ENDIF
C
         WRITE(6,'(A7,A8,A6,A15)') 'DEVICE=',KUFDEV,'SHOT#=',KUFDCG
C
         IF(NSW.EQ.1) THEN
            CALL TR_TIME_UFILE
         ELSEIF(NSW.EQ.2) THEN
            CALL TR_STEADY_UFILE
         ELSEIF(NSW.EQ.3) THEN
            CALL TR_TIME_UFILE_TOPICS
         ENDIF
      ENDIF
C
      NROMAX = NRMAX
      NRMAX  = NRAMAX
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           STEADY STATE UFILE DATA INPUT ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TR_STEADY_UFILE
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
      COMMON /TRUFC3/ NM2CHK
      COMMON /TRUFC4/ NREMAX(2),GRE(NRM,2)
      COMMON /TRERU1/ RTEXU(NTUM,NRMU),   RTIXU(NTUM,NRMU),
     &                RNEXU(NTUM,NRMU)
      COMMON /TRERU2/ RTEXEU(NTUM,NRMU),  RTIXEU(NTUM,NRMU),
     &                RNEXEU(NTUM,NRMU)
      COMMON /TRUFC1/ KDIRX
      DIMENSION RUF(NRMU),F1(NTUM),FAS(NRMP)
      CHARACTER KFID*10,KDIRX*80
C
      MDRGEO=0
      MDAMIN=0
      MDIP=0
      MDBT=0
      MDKAPPA=0
      MDPHIA=0
C
      NTSLS=0
C
C     *** 1D VALUE ***
C
      KFID='RGEO'
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,
     &                 NTXMAX1,NTUM,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDRGEO=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            RRU(NTX)=F1(NTX)
         ENDDO
      ENDIF
C
      KFID='AMIN'
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,
     &                 NTXMAX1,NTUM,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDAMIN=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            RAU(NTX)=F1(NTX)
         ENDDO
      ENDIF
C
      KFID='IP'
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,
     &                 NTXMAX1,NTUM,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDIP=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            RIPU(NTX)=ABS(F1(NTX)*1.D-6)
         ENDDO
      ENDIF
C
      KFID='BT'
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,
     &                 NTXMAX1,NTUM,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDBT=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            BBU(NTX)=ABS(F1(NTX))
         ENDDO
      ENDIF
C
      IF(KUFDEV.EQ.'tftr'.AND.(KUFDCG.EQ.'50862'
     &     .OR.KUFDCG.EQ.'50921'.OR.KUFDCG.EQ.'52527')) THEN
      DO NTX=1,NTXMAX1
         RKAPU(NTX)=1.D0
      ENDDO
      ELSE
      KFID='KAPPA'
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,
     &                 NTXMAX1,NTUM,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDKAPPA=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            RKAPU(NTX)=F1(NTX)
         ENDDO
      ENDIF
      ENDIF
      IF(IERR.NE.0) STOP 'SOME 1D UFILES DO NOT EXIST.'
C
      KFID='PHIA'
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,
     &                 NTXMAX1,NTUM,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDPHIA=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            PHIAU(NTX)=F1(NTX)
         ENDDO
      ENDIF
C
C     *** 2D VALUE ***
C
      AMP=1.D-3
      KFID='TE'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,
     &            NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RTU(1,NR,1)=FAS(NR)
         RT(NR,1)=FAS(NR)
      ENDDO
      PT(1)=RT(1,1)
      PTS(1)=PTMP1
      IF(RHOA.NE.1.D0) PTSA(1)=PTMP2
C
      KFID='TEXP'
      CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTEXU,
     &                   NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      IF(IERR.EQ.0) THEN
         NREMAX(1)=NRFMAX
         DO NR=1,NRFMAX
            GRE(NR,1)=GUCLIP(RUF(NR))
         ENDDO
      ENDIF
      KFID='TEXPEB'
      CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTEXEU,
     &                   NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
C     
      KFID='TI'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,
     &            NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RTU(1,NR,2)=FAS(NR)
         RTU(1,NR,3)=FAS(NR)
         RT(NR,2)=FAS(NR)
         RT(NR,3)=FAS(NR)
      ENDDO
      PT  (2)=RT(1,2)
      PT  (3)=RT(1,3)
      PT  (4)=RT(1,2)
      PTS (2)=PTMP1
      PTS (3)=PTMP1
      PTS (4)=PTMP1
      IF(RHOA.NE.1.D0) THEN
         PTSA(2)=PTMP2
         PTSA(3)=PTMP2
         PTSA(4)=PTMP2
      ENDIF
C
      KFID='TIXP'
      CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTIXU,
     &                   NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      KFID='TIXPEB'
      CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTIXEU,
     &                   NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
C
      AMP=1.D-20
      KFID='NE'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,
     &            NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RNU(1,NR,1)=FAS(NR)
         RNU(1,NR,2)=FAS(NR)
      ENDDO
      PN(1)  =RNU(1,1,1)
      PN(2)  =RNU(1,1,2)-2.D-7
      PN(3)  =1.D-7
      PN(4)  =1.D-7
      PNS(1) =PTMP1
      PNS(2) =PNS(1)-2.D-8
      PNS(3) =1.D-8
      PNS(4) =1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(1)=PTMP2
         PNSA(2)=PNSA(1)-2.D-8
         PNSA(3)=1.D-8
         PNSA(4)=1.D-8
      ENDIF
C
      KFID='NEXP'
      CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RNEXU,
     &                   NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      IF(IERR.EQ.0) THEN
         NREMAX(2)=NRFMAX
         DO NR=1,NRFMAX
            GRE(NR,2)=GUCLIP(RUF(NR))
         ENDDO
      ENDIF
      KFID='NEXPEB'
      CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RNEXEU,
     &                   NRFMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
C
      IF(MDNI.EQ.1) THEN !!!
      IF(NM2CHK.EQ.1) THEN
         KFID='NM2'
      ELSE
         KFID='NM1'
      ENDIF
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,
     &            NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RNU(1,NR,2)=FAS(NR)
         RNU(1,NR,3)=(RNU(1,NR,1)-PZ(2)*RNU(1,NR,2))/PZ(3)-1.D-7
         ZEFFU(1,NR)=(PZ(2)**2*RNU(1,NR,2)+PZ(3)**2*RNU(1,NR,3))
     &              / RNU(1,NR,1)
      ENDDO
      PN(2) =RNU(1,1,2)
      PN(3) =RNU(1,1,3)-1.D-7
      PN(4) =1.D-7
      IF(PTMP1.LE.1.D-3) THEN
         PNS(1)=1.D-2
         PTMP1=1.D-2-1.D-3
      ENDIF
      PNS(2)=PTMP1
      PNS(3)=(PNS(1)-PZ(2)*PNS(2))/PZ(3)-1.D-8
      PNS(4)=1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(2)=PTMP2
         PNSA(3)=(PNSA(1)-PZ(2)*PNSA(2))/PZ(3)-1.D-8
         PNSA(4)=1.D-8
      ENDIF
C
      ELSEIF(MDNI.EQ.2) THEN !!!
C
      AMP=1.D0
      KFID='ZEFFR'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,
     &            NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         ZEFFU(1,NR)=FAS(NR)
         RNU(1,NR,2)=(PZ(3)-FAS(NR))/(PZ(2)*(PZ(3)-PZ(2)))*RNU(1,NR,1)
         RNU(1,NR,3)=(PZ(2)-FAS(NR))/(PZ(3)*(PZ(2)-PZ(3)))*RNU(1,NR,1)
      ENDDO
      PN(2)=RNU(1,1,2)
      PN(3)=RNU(1,1,3)-1.D-7
      PN(4)=1.D-7
      PNS(2)=(PZ(3)-PTMP1)/(PZ(2)*(PZ(3)-PZ(2)))*PNS(1)
      PNS(3)=(PZ(2)-PTMP1)/(PZ(3)*(PZ(2)-PZ(3)))*PNS(1)-1.D-8
      PNS(4)=1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(2)=(PZ(3)-PTMP2)/(PZ(2)*(PZ(3)-PZ(2)))*PNSA(1)
         PNSA(3)=(PZ(2)-PTMP2)/(PZ(3)*(PZ(2)-PZ(3)))*PNSA(1)-1.D-8
         PNSA(4)=1.D-8
      ENDIF
C
      ELSEIF(MDNI.EQ.3) THEN !!!
C
      KFID='NIMP'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,
     &            NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RNU(1,NR,2)=(RNU(1,NR,1)-PZ(3)*FAS(NR))/PZ(2)
         RNU(1,NR,3)=FAS(NR)
         ZEFFU(1,NR)=PZ(2)+PZ(3)*(PZ(3)-PZ(2))*(FAS(NR)/RNU(1,NR,1))
      ENDDO
      PNS(2)=(PNS(1)-PZ(3)*PTMP1)/PZ(2)
      PNS(3)=PTMP1-1.D-8
      PNS(4)=1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(2)=(PNSA(1)-PZ(3)*PTMP2)/PZ(2)
         PNSA(3)=PTMP2-1.D-8
         PNSA(4)=1.D-8
      ENDIF
      ENDIF !!!
C
      DO NR=1,NRMAX
         RNFU(1,NR)=0.D0
         PBMU(1,NR)=0.D0
      ENDDO
      AMP=1.D-20
      KFID='NFAST'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(IERR.EQ.0.AND.FAS(NR).LE.0.D0) FAS(NR)=1.D-10
         RNFU(1,NR)=FAS(NR)
      ENDDO
C
      AMP=1.D0
      KFID='PBEAM'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         PBMU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='Q'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,1,IERR)
      IF(IERR.NE.0) MDLJQ=0
      DO NR=1,NRMAX
         QPU(1,NR)=FAS(NR)
      ENDDO
c$$$      IF(KUFDEV.EQ.'jt60u') THEN
c$$$         DO NR=NRMAX-1,NRMAX
c$$$            QPU(1,NR)=FEDG(DBLE(NR)*DR,DBLE(NR-1)*DR,DBLE(NR-2)*DR,
c$$$     &                     QPU(1,NR-1),QPU(1,NR-2))
c$$$         ENDDO
c$$$         WRITE(6,*) KFID,'IS MODIFIED.'
c$$$      ENDIF
C
      IF(MDNI.NE.2) THEN
      KFID='ZEFFR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,1,IERR)
      DO NR=1,NRMAX
         ZEFFU(1,NR)=FAS(NR)
      ENDDO
      ENDIF
C
      KFID='CURTOT'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      IF(IERR.NE.0) MDLJQ=1
      DO NR=1,NRMAX
         AJU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='CURNBI'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         AJNBU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='BPOL'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,0,IERR)
      DO NR=1,NRMAX
         BPU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='QNBIE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PNBU(1,NR,1)=0.D0
         ELSE
            PNBU(1,NR,1)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QNBII'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PNBU(1,NR,2)=0.D0
         ELSE
            PNBU(1,NR,2)=FAS(NR)
         ENDIF
      ENDDO
C
      DO NR=1,NRMAX
         PICU(1,NR,1)=0.D0
         PICU(1,NR,2)=0.D0
         PECU(1,NR  )=0.D0
         POHU(1,NR  )=0.D0
         WROTU(1,NR )=0.D0
         VTORU(1,NR )=0.D0
      ENDDO
      KFID='QICRHE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PICU(1,NR,1)=0.D0
         ELSE
            PICU(1,NR,1)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QICRHI'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PICU(1,NR,2)=0.D0
         ELSE
            PICU(1,NR,2)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QECH'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PECU(1,NR)=0.D0
         ELSE
            PECU(1,NR)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QECHE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PECU(1,NR)=0.D0
         ELSE
            PECU(1,NR)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QECHI'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PICU(1,NR,2)=PICU(1,NR,2)
         ELSE
            PICU(1,NR,2)=PICU(1,NR,2)+FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QRAD'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         PRLU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='QOHM'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         POHU(1,NR)=FAS(NR)
      ENDDO
C
      AMP=1.D-20
      KFID='SNBIE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         SNBU(1,NR,1)=FAS(NR)
      ENDDO

      KFID='SNBII'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         SNBU(1,NR,2)=FAS(NR)
      ENDDO
C
      KFID='SWALL'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         SWLU(1,NR)=FAS(NR)
      ENDDO
C
      AMP=1.D0
      KFID='VROT'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         WROTU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='VTOR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDVTOR)
      DO NR=1,NRMAX
         VTORU(1,NR)=FAS(NR)
      ENDDO
C
C     *****
C
      MDSUM=MDRGEO+MDAMIN+MDIP+MDBT+MDKAPPA+MDPHIA
      IF(MDSUM.NE.0) THEN
         DO NTX=1,NTXMAX
           IF(MDRGEO .NE.0) RRU  (NTX)=RR
           IF(MDAMIN .NE.0) RAU  (NTX)=RA
           IF(MDIP   .NE.0) RIPU (NTX)=RIPS
           IF(MDBT   .NE.0) BBU  (NTX)=BB
           IF(MDKAPPA.NE.0) RKAPU(NTX)=RKAP
           IF(MDPHIA .NE.0) PHIAU(NTX)=0.D0
         ENDDO
      ENDIF
C
      IF(TMU(2).EQ.0.D0) THEN
         DO NTX1=1,NTXMAX1
            IF(TMU1(NTX1)-TMU(1).LE.1.D-5) NTSLS=NTX1
         ENDDO
         IF(NTSLS.EQ.0) NTSLS=1
      ENDIF
C
C     *** GEOMETORY FACTORS ***
C
      KFID='RMAJOR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      IF(TMU(2).EQ.0.D0) NTSL=NTSLS
      DO NR=1,NRMAX
         RMJRHOU(1,NR)=FAS(NR)
C         ARRHOU(1,NR)=1.D0/RMJRHOU(1,NR)**2
         ARRHOU(1,NR)=1.D0/RRU(NTSL)**2
         TTRHOU(1,NR)=BBU(NTSL)*RRU(NTSL)
      ENDDO
C
      KFID='RMINOR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,0,IERR)
      DO NR=1,NRMAX
         RMNRHOU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='KAPPAR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         EKAPU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='GRHO1'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         AR1RHOU(1,NR)=FAS(NR)
      ENDDO
C
      KFID='GRHO2'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,IERR)
      DO NR=1,NRMAX
         AR2RHOU(1,NR)=FAS(NR)
         ABRHOU(1,NR)=AR2RHOU(1,NR)*ARRHOU(1,NR)
      ENDDO
c$$$      IF(KUFDEV.EQ.'jt60u') THEN
c$$$         DO NR=NRMAX-2,NRMAX
c$$$            AR2RHOU(1,NR)=FEDG(DBLE(NR)*DR,DBLE(NR-1)*DR,DBLE(NR-2)*DR,
c$$$     &                         AR2RHOU(1,NR-1),AR2RHOU(1,NR-2))
c$$$            ABRHOU(1,NR)=AR2RHOU(1,NR)*ARRHOU(1,NR)
c$$$         ENDDO
c$$$         WRITE(6,*) KFID,'IS MODIFIED.'
c$$$      ENDIF
C
      KFID='SURF'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,0,IERR)
      IF(TMU(2).EQ.0.D0) NTSL=NTSLS
      DO NR=1,NRMAX
         DVRHOU(1,NR)=FAS(NR)/AR1RHOU(1,NR)
C         DSRHOU(1,NR)=DVRHOU(1,NR)/(2.D0*PI*RMJRHOU(1,NR))
         DSRHOU(1,NR)=DVRHOU(1,NR)/(2.D0*PI*RRU(NTSL))
      ENDDO
C
      KFID='VOLUME'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,0,IERR)
      DO NR=1,NRMAX
         RJCBU(1,NR)=SQRT(2.D0*PI*PI*RRU(NTSL)/FAS(NRMAX))
      ENDDO
C
C     *****
C
      RR   = RRU(NTSL)
      RA   = RAU(NTSL)
      RIPS = RIPU(NTSL)
      BB   = BBU(NTSL)
      RKAP = RKAPU(NTSL)
      PHIA = PHIAU(NTSL)
C
c$$$      IF(NTXMAX1.EQ.0) THEN
c$$$         NTXMAX1=NTXMAX
c$$$         DO NTX=1,NTXMAX
c$$$            TMU1(NTX)=TMU(NTX)
c$$$         ENDDO
c$$$      ENDIF
c$$$C
c$$$      DO NTX=1,NTXMAX1
c$$$         IF(ABS(TMU1(NTX)-TSLC).LE.1.D-5) THEN
c$$$            RR   = RRU(NTX)
c$$$            RA   = RAU(NTX)
c$$$            RIPS = RIPU(NTX)
c$$$            BB   = BBU(NTX)
c$$$            RKAP = RKAPU(NTX)
c$$$            PHIA = PHIAU(NTX)
c$$$            GOTO 2000
c$$$         ENDIF
c$$$      ENDDO
c$$$      WRITE(6,*) 'XX 1D UFILES NO SLICE TIME THE SAME AS ',
c$$$     &           'THAT OF 2D UFILES!'
c$$$      WRITE(6,*) '## THE PROCESS CONTINUES THROUGH USING THE LAGLANGE ',
c$$$     &           'INTERPOLATION...'
c$$$      CALL LAGLANGE(TSLC,RR  ,TMU1,RRU  ,NTXMAX1,NTUM,IERR)
c$$$      CALL LAGLANGE(TSLC,RA  ,TMU1,RAU  ,NTXMAX1,NTUM,IERR)
c$$$      CALL LAGLANGE(TSLC,RIPS,TMU1,RIPU ,NTXMAX1,NTUM,IERR)
c$$$      CALL LAGLANGE(TSLC,BB  ,TMU1,BBU  ,NTXMAX1,NTUM,IERR)
c$$$      CALL LAGLANGE(TSLC,RKAP,TMU1,RKAPU,NTXMAX1,NTUM,IERR)
c$$$      CALL LAGLANGE(TSLC,PHIA,TMU1,PHIAU,NTXMAX1,NTUM,IERR)
c$$$ 2000 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
C           TIME EVOLUTIONAL UFILE DATA INPUT ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TR_TIME_UFILE
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC2/ NTAMAX
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
      COMMON /TRUFC3/ NM2CHK
      DIMENSION FAT(NTUM,NRMP),PV(NTUM),PVA(NTUM)
      CHARACTER KFID*10
C
      ICK=0
      TMUMAX=0.D0
C
      MDRGEO=0
      MDAMIN=0
      MDIP=0
      MDBT=0
      MDKAPPA=0
      MDPHIA=0
      MDPNBI=0
C
C     *** 1D VALUE ***
C
      KFID='RGEO'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RRU  ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDRGEO=IERR
         IERR=0
      ENDIF
C
      KFID='AMIN'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RAU  ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDAMIN=IERR
         IERR=0
      ENDIF
C
      KFID='IP'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RIPU ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDIP=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            RIPU(NTX)=ABS(RIPU(NTX)*1.D-6)
         ENDDO
      ENDIF
C
      KFID='BT'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,BBU  ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDBT=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX1
            BBU(NTX)=ABS(BBU(NTX))
         ENDDO
      ENDIF
C
      KFID='KAPPA'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RKAPU,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDKAPPA=IERR
         IERR=0
      ENDIF
c$$$      IF(IERR.EQ.1) THEN
c$$$         DO NTX=1,NTXMAX1
c$$$            RKAPU(NTX)=1.D0
c$$$         ENDDO
c$$$      ENDIF
C
      KFID='PHIA'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,PHIAU,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDPHIA=IERR
         IERR=0
      ENDIF
C
      KFID='PNBI'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,PNBIU,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDPNBI=IERR
         IERR=0
      ENDIF
C
C     *** 2D VALUE ***
C
      AMP=1.D-3
      KFID='TE'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NR=1,NRMAX
         DO NTX=1,NTXMAX
            RTU(NTX,NR,1)=FAT(NTX,NR)
         ENDDO
         RT(NR,1)=RTU(1,NR,1)
      ENDDO
      PT(1)=RT(1,1)
      DO NTX=1,NTXMAX
         PTSU(NTX,1)=PV(NTX)
      ENDDO
      PTS(1)=PTSU(1,1)
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PTSUA(NTX,1)=PVA(NTX)
         ENDDO
         PTSA(1)=PTSUA(1,1)
      ENDIF
C
      KFID='TI'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NR=1,NRMAX
         DO NS=2,3
            DO NTX=1,NTXMAX
               RTU(NTX,NR,NS)=FAT(NTX,NR)
            ENDDO
            RT(NR,NS)=RTU(1,NR,NS)
         ENDDO
      ENDDO
      PT(2)=RT(1,2)
      PT(3)=RT(1,3)
      PT(4)=RT(1,2)
      DO NS=2,4
         DO NTX=1,NTXMAX
            PTSU (NTX,NS)=PV(NTX)
         ENDDO
         PTS (NS)=PTSU(1,NS)
      ENDDO
      IF(RHOA.NE.1.D0) THEN
         DO NS=2,4
            DO NTX=1,NTXMAX
               PTSUA(NTX,NS)=PVA(NTX)
            ENDDO
            PTSA(NS)=PTSUA(1,NS)
         ENDDO
      ENDIF
C
      AMP=1.D-20
      KFID='NE'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RNU(NTX,NR,1)=FAT(NTX,NR)
            RNU(NTX,NR,2)=FAT(NTX,NR)-2.D-7
            RNU(NTX,NR,3)=1.D-7
            RNU(NTX,NR,4)=1.D-7
         ENDDO
      ENDDO
      PN(1)=RNU(1,1,1)
      PN(2)=RNU(1,1,2)-2.D-7
      PN(3)=1.D-7
      PN(4)=1.D-7
      DO NTX=1,NTXMAX
         PNSU(NTX,1)=PV(NTX)
         PNSU(NTX,2)=PV(NTX)-2.D-8
         PNSU(NTX,3)=1.D-8
         PNSU(NTX,4)=1.D-8
      ENDDO
      PNS(1)=PNSU(1,1)
      PNS(2)=PNSU(1,2)
      PNS(3)=PNSU(1,3)
      PNS(4)=PNSU(1,4)
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PNSUA(NTX,1)=PVA(NTX)
            PNSUA(NTX,2)=PVA(NTX)-2.D-8
            PNSUA(NTX,3)=1.D-8
            PNSUA(NTX,4)=1.D-8
         ENDDO
         PNSA(1)=PNSUA(1,1)
         PNSA(2)=PNSUA(1,2)
         PNSA(3)=PNSUA(1,3)
         PNSA(4)=PNSUA(1,4)
      ENDIF
C
      IF(MDNI.EQ.1) THEN !!!
      IF(NM2CHK.EQ.1) THEN
         KFID='NM2'
      ELSE
         KFID='NM1'
      ENDIF
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RNU(NTX,NR,2)=FAT(NTX,NR)
            RNU(NTX,NR,3)=(RNU(NTX,NR,1)-PZ(2)*RNU(NTX,NR,2))/PZ(3) 
            ZEFFU(NTX,NR)=(PZ(2)**2*RNU(NTX,NR,2)
     &                    +PZ(3)**2*RNU(NTX,NR,3))/RNU(NTX,NR,1)
         ENDDO
      ENDDO
      PN(2)=RNU(1,1,2)
      PN(3)=RNU(1,1,3)-1.D-7
      PN(4)=1.D-7
      DO NTX=1,NTXMAX
         PNSU(NTX,2)=PV(NTX)
         PNSU(NTX,3)=(PNSU(NTX,1)-PZ(2)*PNSU(NTX,2))/PZ(3)-1.D-8
         PNSU(NTX,4)=1.D-8
      ENDDO
      PNS(2)=PNSU(1,2)
      PNS(3)=PNSU(1,3)
      PNS(4)=PNSU(1,4)
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PNSUA(NTX,2)=PVA(NTX)
            PNSUA(NTX,3)=(PNSUA(NTX,1)-PZ(2)*PNSUA(NTX,2))/PZ(3)-1.D-8
            PNSUA(NTX,4)=1.D-8
         ENDDO
         PNSA(2)=PNSUA(1,2)
         PNSA(3)=PNSUA(1,3)
         PNSA(4)=PNSUA(1,4)
      ENDIF
C
      ELSEIF(MDNI.EQ.2) THEN !!!
C
      AMP=1.D0
      KFID='ZEFFR'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            ZEFFU(NTX,NR)=FAT(NTX,NR)
            RNU(NTX,NR,2)=(PZ(3)-FAT(NTX,NR))/(PZ(2)*(PZ(3)-PZ(2)))
     &                   *RNU(NTX,NR,1)
            RNU(NTX,NR,3)=(PZ(2)-FAT(NTX,NR))/(PZ(3)*(PZ(2)-PZ(3)))
     &                   *RNU(NTX,NR,1)
            IF(RNU(NTX,NR,3).LT.0.D0) STOP 'NIMP IS BELOW ZERO.'
         ENDDO
      ENDDO
      PN(2)=RNU(1,1,2)
      PN(3)=RNU(1,1,3)-1.D-7
      PN(4)=1.D-7
      DO NTX=1,NTXMAX
         PNSU(NTX,2)=(PZ(3)-PV(NTX))/(PZ(2)*(PZ(3)-PZ(2)))*PNSU(NTX,1)
         PNSU(NTX,3)=(PZ(2)-PV(NTX))/(PZ(3)*(PZ(2)-PZ(3)))*PNSU(NTX,1)
     &              -1.D-8
         PNSU(NTX,4)=1.D-8
      ENDDO
      PNS(2)=PNSU(1,2)
      PNS(3)=PNSU(1,3)
      PNS(4)=PNSU(1,4)
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PNSUA(NTX,2)=(PZ(3)-PVA(NTX))/(PZ(2)*(PZ(3)-PZ(2)))
     &                  *PNSUA(NTX,1)
            PNSUA(NTX,3)=(PZ(2)-PVA(NTX))/(PZ(3)*(PZ(2)-PZ(3)))
     &                  *PNSUA(NTX,1)-1.D-8
            PNSUA(NTX,4)=1.D-8
         ENDDO
         PNSA(2)=PNSUA(1,2)
         PNSA(3)=PNSUA(1,3)
         PNSA(4)=PNSUA(1,4)
      ENDIF
C
      ELSEIF(MDNI.EQ.3) THEN !!!
C
      KFID='NIMP'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RNU(NTX,NR,2)=(RNU(NTX,NR,1)-PZ(3)*FAT(NTX,NR))/PZ(2)
            RNU(NTX,NR,3)=FAT(NTX,NR)-1.D-7
            ZEFFU(NTX,NR)=PZ(2)+PZ(3)
     &                  *(PZ(3)-PZ(2))*(FAT(NTX,NR)/RNU(NTX,NR,1))
         ENDDO
      ENDDO
      PN(2)=RNU(1,1,2)
      PN(3)=RNU(1,1,3)-1.D-7
      PN(4)=1.D-7
      DO NTX=1,NTXMAX
         PNSU(NTX,2)=(PNSU(NTX,1)-PZ(3)*PV(NTX))/PZ(2)
         PNSU(NTX,3)=PV(NTX)-1.D-8
         PNSU(NTX,4)=1.D-8
      ENDDO
      PNS(2)=PNSU(1,2)
      PNS(3)=PNSU(1,3)
      PNS(4)=PNSU(1,4)
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PNSUA(NTX,2)=(PNSUA(NTX,1)-PZ(3)*PVA(NTX))/PZ(2)
            PNSUA(NTX,3)=PVA(NTX)-1.D-8
            PNSUA(NTX,4)=1.D-8
         ENDDO
         PNSA(2)=PNSUA(1,2)
         PNSA(3)=PNSUA(1,3)
         PNSA(4)=PNSUA(1,4)
      ENDIF
      ENDIF !!!
C
      AMP=1.D-20
      KFID='NFAST'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RNFU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      AMP=1.D0
      KFID='PBEAM'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,PBMU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      KFID='Q'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,QPU ,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,IERR)
      IF(IERR.NE.0) MDLJQ=0
      IF(MDNI.NE.2) THEN
      KFID='ZEFFR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,ZEFFU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,IERR)
      ENDIF
      KFID='CURTOT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AJU ,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      IF(IERR.NE.0) MDLJQ=1
      KFID='CURNBI'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AJNBU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      KFID='BPOL'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,BPU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,1,0,ICK,IERR)
      KFID='QNBIE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PNBU(NTX,NR,1)=0.D0
            ELSE
               PNBU(NTX,NR,1)=FAT(NTX,NR)*CNB
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QNBII'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT ,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PNBU(NTX,NR,2)=0.D0
            ELSE
               PNBU(NTX,NR,2)=FAT(NTX,NR)*CNB
            ENDIF
         ENDDO
      ENDDO
C
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            PICU(NTX,NR,1)=0.D0
            PICU(NTX,NR,2)=0.D0
            PECU(NTX,NR  )=0.D0
            POHU(NTX,NR  )=0.D0
         ENDDO
      ENDDO
      KFID='QICRHE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT ,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PICU(NTX,NR,1)=0.D0
            ELSE
               PICU(NTX,NR,1)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QICRHI'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT ,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PICU(NTX,NR,2)=0.D0
            ELSE
               PICU(NTX,NR,2)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QECH'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT ,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PECU(NTX,NR)=0.D0
            ELSE
               PECU(NTX,NR)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QRAD'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,PRLU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
      KFID='QOHM'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,POHU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
C     *****
C
      MDSUM=MDRGEO+MDAMIN+MDIP+MDBT+MDKAPPA+MDPHIA+MDPNBI
      IF(MDSUM.NE.0) THEN
         NTXMAX1=NTXMAX
         DO NTX=1,NTXMAX
            TMU1(NTX)=TMU(NTX)
            IF(MDRGEO .NE.0) RRU  (NTX)=RR
            IF(MDAMIN .NE.0) RAU  (NTX)=RA
            IF(MDIP   .NE.0) RIPU (NTX)=RIPS
            IF(MDBT   .NE.0) BBU  (NTX)=BB
            IF(MDKAPPA.NE.0) RKAPU(NTX)=RKAP
            IF(MDPHIA .NE.0) PHIAU(NTX)=0.D0
            IF(MDPNBI .NE.0) PNBIU(NTX)=0.D0
         ENDDO
      ENDIF
C
C     *** GEOMETORY FACTORS ***
C
      KFID='RMAJOR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RMJRHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
C            ARRHOU(NTX,NR)=1.D0/RMJRHOU(NTX,NR)**2
C            EPS=RMN*RAU(NTX)/RRU(NTX)
C            ARRHOU(NTX,NR)=1.D0/(RRU(NTX)**2*(1.D0-EPS**2)**1.5)
C            ARRHOU(NTX,NR)=(1.D0+1.5D0*EPS**2)/RRU(NTX)**2
            ARRHOU(NTX,NR)=1.D0/RRU(NTX)**2
            TTRHOU(NTX,NR)=RRU(NTX)*BBU(NTX)
         ENDDO
         IF(MDRGEO.NE.0.AND.KUFDEV.EQ.'jt60u') RRU(NTX)=RMJRHOU(NTX,1)
      ENDDO
C
      KFID='RMINOR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RMNRHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,0,ICK,IERR)
      IF(MDAMIN.NE.0.AND.KUFDEV.EQ.'jt60u') THEN
         DO NTX=1,NTXMAX
            RAU(NTX)=RMNRHOU(NTX,NRMAX)
         ENDDO
      ENDIF
C
      KFID='KAPPAR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,EKAPU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      IF(MDKAPPA.NE.0.AND.KUFDEV.EQ.'jt60u') THEN
         DO NTX=1,NTXMAX
            RKAPU(NTX)=EKAPU(NTX,NRMAX)
         ENDDO
      ENDIF
C
      KFID='GRHO1'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AR1RHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      IF(IERR.NE.0.AND.KUFDEV.EQ.'jt60u') THEN
         DO NTX=1,NTXMAX
            AR1RHOU(NTX,NR)=1.D0/RMNRHOU(NTX,NRMAX)
         ENDDO
      ENDIF
C
      KFID='GRHO2'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AR2RHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            ABRHOU(NTX,NR)=AR2RHOU(NTX,NR)*ARRHOU(NTX,NR)
            IF(IERR.NE.0.AND.KUFDEV.EQ.'jt60u') THEN
               AR2RHOU(NTX,NR)=1.D0/RMNRHOU(NTX,NRMAX)**2
               ABRHOU(NTX,NR)=AR2RHOU(NTX,NR)*ARRHOU(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='SURF'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            DVRHOU(NTX,NR)=FAT(NTX,NR)/AR1RHOU(NTX,NR)
C            DSRHOU(NTX,NR)=DVRHOU(NTX,NR)/(2.D0*PI*RMJRHOU(NTX,NR))
            DSRHOU(NTX,NR)=DVRHOU(NTX,NR)/(2.D0*PI*RRU(NTX))
         ENDDO
      ENDDO
C
      KFID='VOLUME'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RJCBU(NTX,NR)=SQRT(2.D0*PI*PI*RRU(NTX)/FAT(NTX,NRMAX))
         ENDDO
      ENDDO
C
      IF(NTAMAX.LT.NTAMAX1) NTAMAX=NTAMAX1
      RETURN
      END
C
C
C     *****************************************************************
C
C           TIME EVOLUTIONAL UFILE DATA INPUT ROUTINE FOR TOPICS
C
C     *****************************************************************
C
      SUBROUTINE TR_TIME_UFILE_TOPICS
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC2/ NTAMAX
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
      DIMENSION FAT(NTUM,NRMP),TMP(NTUM),PV(NTUM),PVA(NTUM)
      CHARACTER KFID*10
C
      ICK=0
      TMUMAX=0.D0
C
C     *** 1D VALUE ***
C
      KFID='RGEO'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RRU  ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      MDRGEO=IERR
C
      KFID='AMIN'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RAU  ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      MDAMIN=IERR
C
      KFID='IP'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RIPU ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      DO NTX=1,NTXMAX1
         RIPU(NTX)=ABS(RIPU(NTX)*1.D-6)
      ENDDO
      MDIP=IERR
C
      KFID='BT'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,BBU  ,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      MDBT=IERR
C
      KFID='KAPPA'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RKAPU,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      MDKAPPA=IERR
C
      KFID='VOL'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,TMP,
     &          NTAMAX1,NTXMAX1,TMUMAX,ICK,IERR)
      MDKAPPA=IERR
C
C     *** 2D VALUE ***
C
      AMP=1.D-3
      KFID='TE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RTU(NTX,NR,1)=FAT(NTX,NR)
         ENDDO
      ENDDO
C
      KFID='TI'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RTU(NTX,NR,2)=FAT(NTX,NR)
            RTU(NTX,NR,3)=FAT(NTX,NR)
         ENDDO
      ENDDO
C
      AMP=1.D-20
      KFID='NE'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RNU(NTX,NR,1)=FAT(NTX,NR)
C            RNU(NTX,NR,2)=FAT(NTX,NR)
         ENDDO
         PNSU (NTX,1)=PV (NTX)
C         PNSU (NTX,2)=PV (NTX)
      ENDDO
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PNSUA(NTX,1)=PVA(NTX)
C            PNSUA(NTX,2)=PVA(NTX)
         ENDDO
      ENDIF
C
      KFID='NIMP'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,
     &            NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RNU(NTX,NR,2)=(RNU(NTX,NR,1)-PZ(3)*FAT(NTX,NR))/PZ(2)
            RNU(NTX,NR,3)=FAT(NTX,NR)
         ENDDO
         PNSU(NTX,2)=(PNSU(NTX,1)-PZ(3)*PV(NTX))/PZ(2)
         PNSU(NTX,3)=PV(NTX)
      ENDDO
      IF(RHOA.NE.1.D0) THEN
         DO NTX=1,NTXMAX
            PNSUA(NTX,2)=(PNSUA(NTX,1)-PZ(3)*PVA(NTX))/PZ(2)
            PNSUA(NTX,3)=PVA(NTX)
         ENDDO
      ENDIF
C
      KFID='NM1'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            RNU(NTX,NR,4)=FAT(NTX,NR)
         ENDDO
      ENDDO
C
      AMP=1.D0
      KFID='ZEFFR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            ZEFFU(NTX,NR)=FAT(NTX,NR)
         ENDDO
      ENDDO
C
      KFID='Q'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,QPU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,IERR)
      KFID='CURTOT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AJU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      IF(KUFDEV.EQ.'X'.AND.KUFDCG.EQ.'14') THEN
         KFID='CURBS'
         CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AJBSU,AMP,
     &              NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      ENDIF
C
      KFID='QNBIE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PNBU(NTX,NR,1)=0.D0
            ELSE
               PNBU(NTX,NR,1)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
      KFID='QNBII'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PNBU(NTX,NR,2)=0.D0
            ELSE
               PNBU(NTX,NR,2)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            PICU(NTX,NR,1)=0.D0
            PICU(NTX,NR,2)=0.D0
            PECU(NTX,NR  )=0.D0
            POHU(NTX,NR  )=0.D0
         ENDDO
      ENDDO
C      KFID='QICRHE'
      KFID='QLHE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PICU(NTX,NR,1)=0.D0
            ELSE
               PICU(NTX,NR,1)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C      KFID='QICRHI'
      KFID='QLHI'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PICU(NTX,NR,2)=0.D0
            ELSE
               PICU(NTX,NR,2)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QECH'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            IF(FAT(NTX,NR).LT.0.D0) THEN
               PECU(NTX,NR)=0.D0
            ELSE
               PECU(NTX,NR)=FAT(NTX,NR)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QRAD'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,PRLU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
      KFID='QOHM'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,POHU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
C     *** GEOMETORY FACTORS ***
C
      KFID='RMAJOR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RMJRHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
      KFID='RMINOR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RMNRHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,0,ICK,IERR)
C
      KFID='GRHO1'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AR1RHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
      KFID='GRHO2'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AR2RHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
      KFID='VRO'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,DVRHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
C            DSRHOU(NTX,NR)=DVRHOU(NTX,NR)/(2.D0*PI*RMJRHOU(NTX,NR))
            DSRHOU(NTX,NR)=DVRHOU(NTX,NR)/(2.D0*PI*RR)
         ENDDO
      ENDDO
C
      KFID='AAT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,ARRHOU,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
C
      KFID='HDT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            TTRHOU(NTX,NR)=FAT(NTX,NR)/ARRHOU(NTX,NR)
         ENDDO
      ENDDO
C
      KFID='CKT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,
     &           NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,IERR)
      DO NTX=1,NTXMAX
         DO NR=1,NRMAX
            ABRHOU(NTX,NR)=FAT(NTX,NR)/DVRHOU(NTX,NR)**2
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         RG(NR)    =DBLE(NR)*DR
         EPSRHO(NR)=RA*RG(NR)/RR
         EKAP(NR)=RKAP
      ENDDO
C
C     *** 1D VALUE ***
C
      MDALL=MDRGEO+MDAMIN+MDIP+MDBT+MDKAPPA
      IF(MDALL.NE.0) THEN
         NTXMAX1=NTXMAX
         DO NTX=1,NTXMAX1
            IF(MDRGEO.NE.0)  RRU(NTX)=RR
            IF(MDAMIN.NE.0)  RAU(NTX)=RA
            IF(MDIP.NE.0)    RIPU(NTX)=RIPS
            IF(MDBT.NE.0)    BBU(NTX)=BB
            IF(MDKAPPA.NE.0) RKAPU(NTX)=RKAP
            TMU1(NTX)=TMU(NTX)
         ENDDO
      ENDIF
C
      DO NTX=1,NTXMAX1
         DO NR=1,NRMAX
            RJCBU(NTX,NR)=SQRT(2.D0*PI*PI*RRU(NTX)/TMP(NTX))
         ENDDO
      ENDDO

C     ***
C
      RETURN
      END
C
C
C     ***********************************************************
C
C           1D UFILE READER
C
C     ***********************************************************
C
C     input:
C
C     KFID        : Variable Name
C     DT          : Time Step Width
C     NTLMAX      : Maximum Time Step for TASK/TR
C     TLMAX       : Maximum Time for TASK/TR
C     ICK         : Check Indicator of UFILE Consistency
C
C     output:
C
C     TL(NTUM)    : Total Time Data (The Number of DT)
C     F1(NTUM)    : Functional Values
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF1D(KFID,KUFDEV,KUFDCG,DT,TL,F1,
     &                NTLMAX,NTXMAX,TLMAX,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      COMMON /TRUFC1/ KDIRX
      DIMENSION TL(NTUM),F1(NTUM)
      CHARACTER KFID*10,KDIRX*80
C
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TL,F1,
     &                 NTXMAX,NTUM,MDCHK,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX UF1D:',KFID,'UFILE HAS AN ERROR!'
         RETURN
      ENDIF
C
      DO NTX=2,NTXMAX
         TL(NTX)=(TL(NTX)-TL(1))
      ENDDO
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(ICK.EQ.2.AND.TLMAX.EQ.0.D0) GOTO 1000
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX TRFILE:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
            STOP
         ENDIF
      ENDIF
 1000 CONTINUE
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
      IF(ICK.NE.2) ICK=1
C
      RETURN
      END
C
C     *** FOR GRAPHIC ROUTINE ***
C
      SUBROUTINE UF1DG(KFID,KUFDEV,KUFDCG,GTL,TL,FOUT,AMP,NINMAX,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION TL(NTUM),F1(NTUM),FOUT(NTUM),GTL(NTM)
      COMMON /TRUFC1/ KDIRX
      CHARACTER KFID*10,KDIRX*80
C
      CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TL,F1,
     &                 NTXMAX,NTUM,MDCHK,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX UF1DG:',KFID,'UFILE HAS AN ERROR!'
         RETURN
      ENDIF
C
      DO NTX=2,NTXMAX
         TL(NTX)=(TL(NTX)-TL(1))
      ENDDO
      TL(1)=0.D0
      DO NIN=1,NINMAX
         TLN=DBLE(GTL(NIN))
         CALL LAGLANGE(TLN,F0,TL,F1,NTXMAX,NTUM,IERRL)
C         IF(IERRL.NE.0) 
C     &        WRITE(6,600) "XX TRFILE: LAGLAN ",KFID,": IERR=",IERR
         FOUT(NIN)=F0
      ENDDO
C     
      IF(IERR.EQ.1) THEN
         DO NIN=1,NINMAX
            FOUT(NIN)=0.D0
         ENDDO
      ELSE
         DO NIN=1,NINMAX
            FOUT(NIN)=FOUT(NIN)*AMP
         ENDDO
      ENDIF
C
 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END
C
C     ***********************************************************
C
C           2D UFILE READER
C
C     ***********************************************************
C
C     *** FOR STEADY STATE SIMULATION ***
C
C     input:
C
C     KFID        : Variable Name
C     DR          : Radial Step Width
C     AMP         : Amplitude Factor of FOUT
C     NRMAX       : Radial Node Number
C     NSW         : Mesh Selector (0:RM, 1:RG)
C     ID          : Boundary Condition in the center of edge for Spline
C
C     output:
C
C     TL(NTUM)    : Total Time Data (The Number of DT)
C     FOUT(NRMP)  : Functional Values
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF2DS(KFID,KUFDEV,KUFDCG,DR,TL,FOUT,AMP,
     &                 NRMAX,NSW,ID,IERR)
C
      INCLUDE 'trcom0.inc'
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
      COMMON /TRUFC1/ KDIRX
      DIMENSION RL(NRMU),TL(NTUM),F2(NTUM,NRMU),FOUT(NRMP)
      DIMENSION U(4,NRMU)
      CHARACTER KFID*10,KDIRX*80
C
      CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,
     &                  NRLMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      IF(KUFDEV.EQ.'jet'.AND.(KUFDCG.EQ.'35156'.OR.KUFDCG.EQ.'35171')
     &     .AND.KFID.EQ.'GRHO1') THEN
         DO NTX=1,NTXMAX
            F2(NTX,2)=FCTR(RL(3),RL(4),F2(NTX,3),F2(NTX,4))
            F2(NTX,1)=FCTR(RL(2),RL(3),F2(NTX,2),F2(NTX,3))
         ENDDO
      ENDIF
      IERRP=IERR
      CALL PRETREAT0(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,ID,IERRP)
      IF(NRLMAX.LE.1) RETURN
      DO NRL=1,NRMAX
         IF(NSW.EQ.0) THEN
            RSL=(DBLE(NRL)-0.5D0)*DR
         ELSEIF(NSW.EQ.1) THEN
            RSL= DBLE(NRL)       *DR
         ENDIF
         CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERRS)
         IF(IERRS.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
         FOUT(NRL)=F0*AMP
      ENDDO
C
 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END
C
C     *****************************************************
C
C     input:
C
C     KFID        : Variable name
C     DR          : Radial Step Width
C     AMP         : Amplitude Factor of FOUT
C     RHOA        : Normalized Radius at Arbitrary Surface
C     NRAMAX      : Radial Node Number Corresponding to RHOA
C     NRMAX       : Radial Node Number
C
C     output:
C
C     PV          : Surface(Peripheral) Functional Value
C     PVS         : Surface(Peripheral) Functional Value 
C                   Corresponding to RHOA
C     FOUT(NRMP)  : Functional Values
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF2DSP(KFID,KUFDEV,KUFDCG,DR,PV,PVA,FOUT,AMP,RHOA,
     &                  NRAMAX,NRMAX,IERR)
C
      INCLUDE 'trcom0.inc'
      COMMON /TRUFC1/ KDIRX
      DIMENSION RL(NRMU),TL(NTUM),F2(NTUM,NRMU),FOUT(NRMP)
      DIMENSION U(4,NRMU)
      CHARACTER KFID*10,KDIRX*80
C
      CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,
     &                  NRLMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      CALL PRETREAT0(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,1,IERR)
      IF(NRLMAX.LE.1) RETURN
      DO NRL=1,NRMAX
         RMN=(DBLE(NRL)-0.5D0)*DR
         CALL SPL1DF(RMN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
         FOUT(NRL)=F0*AMP
      ENDDO
C
      RGN=DBLE(NRMAX)*DR
      CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
      IF(IERR.NE.0)
     &     WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
      PV=F0*AMP
      IF(RHOA.NE.1.D0) THEN
         RGN=DBLE(NRAMAX)*DR
         CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
         PVA=F0*AMP
      ENDIF
C
 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END
C
C     *** FOR TIME EVOLUTION SIMULATION ***
C
C     input:
C
C     KFID        : Variable name
C     DR          : Radial Step Width
C     DT          : Time Step Width
C     AMP         : Amplitude Factor of FOUT
C     NRMAX       : Radial Node Number
C     TLMAX       : Maximum Time
C     NSW         : Mesh Selector (0:RM, 1:RG)
C     ID          : Boundary Condition in the center of edge for Spline
C     ICK         : Check Indicator
C
C     output:
C
C     TL(NTUM)    : Total Time Data (The Number of DT)
C     FOUT(NRMP)  : Functional Values
C     NTAMAX      : Maximum Time Step Number Corresponding to RHOA
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TL,FOUT,AMP,
     &                 NTLMAX,NTXMAX,NRMAX,TLMAX,NSW,ID,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION RL(NRMU),TL(NTUM),F2(NTUM,NRMU),FOUT(NTUM,NRMP)
      DIMENSION U(4,NRMU),DERIV(NRMU)
      DIMENSION TMP0(NRMU)
      COMMON /TRUFC1/ KDIRX
      CHARACTER KFID*10,KDIRX*80
C
      CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,
     &                  NRLMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      IF(NRLMAX.LE.1) RETURN
      IF(IERR.NE.0) THEN
         DO NTX=1,NTUM
            DO NRX=1,NRMP
               FOUT(NTX,NRX)=0.D0
            ENDDO
         ENDDO
         RETURN
      ENDIF
      DERIV(1)=0.D0
      DERIV(NRLMAX)=0.D0
C
      DO NTX=2,NTXMAX
         TL(NTX)=(TL(NTX)-TL(1))
      ENDDO
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(ICK.EQ.2.AND.TLMAX.EQ.0.D0) GOTO 1000
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX TRFILE:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
            STOP
         ENDIF
      ENDIF
 1000 CONTINUE
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
      IF(ICK.NE.2) ICK=1
C
      DO NTX=1,NTXMAX
         DO NRL=1,NRLMAX
            TMP0(NRL)=F2(NTX,NRL)
         ENDDO
         CALL SPL1D(RL,TMP0,DERIV,U,NRLMAX,ID,IERR)
         DO NRL=1,NRMAX
            IF(NSW.EQ.0) THEN
               RSL=(DBLE(NRL)-0.5D0)*DR
            ELSEIF(NSW.EQ.1) THEN
               RSL= DBLE(NRL)       *DR
            ENDIF
            CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
            FOUT(NTX,NRL)=F0*AMP
         ENDDO
      ENDDO
C
 600  FORMAT(' ',A18,A10,A7,I2)      
      RETURN
      END
C
C     *** FOR THE VARIABLE WE'D LIKE TO USE PERIPHERAL VALUES ***
C
      SUBROUTINE UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TL,FOUT,AMP,
     &                  NTLMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TLMAX,
     &                  NSW,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION RL(NRMU),TL(NTUM),F2(NTUM,NRMU),FOUT(NTUM,NRMP)
      DIMENSION U(4,NRMU),DERIV(NRMU)
      DIMENSION TMP0(NRMU)
      DIMENSION PV(NTUM),PVA(NTUM)
      COMMON /TRUFC1/ KDIRX
      CHARACTER KFID*10,KDIRX*80
C
      CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,
     &                  NRLMAX,NTXMAX,NRMU,NTUM,MDCHK,IERR)
      IF(NRLMAX.LE.1) RETURN
      IF(IERR.NE.0) THEN
         DO NTX=1,NTUM
            DO NRX=1,NRMP
               FOUT(NTX,NRX)=0.D0
            ENDDO
         ENDDO
         RETURN
      ENDIF
C
      DO NTX=2,NTXMAX
         TL(NTX)=(TL(NTX)-TL(1))
      ENDDO
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(ICK.EQ.2.AND.TLMAX.EQ.0.D0) GOTO 1000
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX TRFILE:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
C            STOP
         ENDIF
      ENDIF
 1000 CONTINUE
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
      IF(ICK.NE.2) ICK=1
C
      DO NTX=1,NTXMAX
C         if(kfid.eq.'TI') write(6,*) ntx,f2(ntx,nrlmax)
         DO NRL=1,NRLMAX
            TMP0(NRL)=F2(NTX,NRL)
         ENDDO
         CALL SPL1D(RL,TMP0,DERIV,U,NRLMAX,1,IERR)
         DO NRL=1,NRMAX
            IF(NSW.EQ.0) THEN
               RSL=(DBLE(NRL)-0.5D0)*DR
            ELSEIF(NSW.EQ.1) THEN
               RSL= DBLE(NRL)       *DR
            ENDIF
            CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
            FOUT(NTX,NRL)=F0*AMP
         ENDDO
C
         RGN=DBLE(NRMAX)*DR
         CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
         PV(NTX)=F0*AMP
C         if(kfid.eq.'TI') write(6,*) ntx,pv(ntx)
         IF(RHOA.NE.1.D0) THEN
            RGN=DBLE(NRAMAX)*DR
            CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
            PVA(NTX)=F0*AMP
         ENDIF
      ENDDO
C
 600  FORMAT(' ',A18,A10,A7,I2)      
      RETURN
      END
C
C     ***********************************************************
C
C           PRETREATMENT SUBROUTINE FOR UFILE INTERFACE
C
C     ***********************************************************
C
C     This subroutine is only used if MDLUF=2 and any subroutines
C     do not call this except UF2DS and UF2DSP.
C     In this subroutine, one can choose arbitrary slice time if one
C     did not choose it ahead of time. If one choose the time,
C     this checks its value consistent with the range of time.
C     One of the main function in this section is to make "Spline Array"
C     in order to interpolate various profiles radially.
C
      SUBROUTINE PRETREAT0(KFID,RL,TL,F2,U,NRFMAX,NTXMAX,ID,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION RL(NRMU),TL(NTUM),F2(NTUM,NRMU)
      DIMENSION DERIV(NRMU),U(4,NRMU),TMP(NRMU)
      CHARACTER KFID*10
C
      DERIV(1)=0.D0
      DERIV(NRFMAX)=0.D0
C
      IF(IERR.EQ.1) THEN
         DO NRF=1,NRMU
         DO NA=1,4
            U(NA,NRF)=0.D0
         ENDDO
         ENDDO
         IF(KFID.EQ.'CURTOT'.AND.IERR.NE.0) MDLJQ=1
         IERR=0
         RETURN
      ENDIF
C
      NTSL=1
      IF(NTXMAX.NE.1) THEN
         IF(TIME_INT.LE.0.D0) THEN
 100        WRITE(6,500) 'INPUT ARBITRARY TIME:',TL(1),' -',TL(NTXMAX)
            READ(5,*,ERR=100) TIME_INT
            IF(TIME_INT.LT.TL(1).OR.TIME_INT.GT.TL(NTXMAX)) GOTO 100
            DO NTX=1,NTXMAX
               IF(ABS(TL(NTX)-TIME_INT).LE.1.D-5) THEN
                  NTSL=NTX
                  GOTO 1000
               ENDIF
            ENDDO
C
            TL_MIN=TL(NTXMAX)
            DO NTX=1,NTXMAX
               TL_MIN_OLD=TL_MIN
               TL_MIN=MIN(ABS(TL(NTX)-TIME_INT),TL_MIN)
               IF(TL_MIN_OLD.EQ.TL_MIN) THEN
                  NTX_MIN=NTX-1
                  GOTO 200
               ENDIF
            ENDDO
         ELSE
            IF(TIME_INT.GE.TL(1).AND.TIME_INT.LE.TL(NTXMAX)) THEN
               DO NTX=1,NTXMAX
                  IF(ABS(TL(NTX)-TIME_INT).LE.1.D-5) THEN
                     NTSL=NTX
                     GOTO 1000
                  ENDIF
               ENDDO
C
               TL_MIN=TL(NTXMAX)
               DO NTX=1,NTXMAX
                  TL_MIN_OLD=TL_MIN
                  TL_MIN=MIN(ABS(TL(NTX)-TIME_INT),TL_MIN)
                  IF(TL_MIN_OLD.EQ.TL_MIN) THEN
                     NTX_MIN=NTX-1
                     GOTO 200
                  ENDIF
               ENDDO
            ENDIF
            WRITE(6,*)
     &           '## DESIGNATED TIME_INT IS NOT IN THE RANGE.'
            WRITE(6,*) 'TIME_INT=',SNGL(TIME_INT),
     &           'RANGE=',SNGL(TL(1)),'-',SNGL(TL(NTXMAX))
            STOP
         ENDIF
C
 200     WRITE(6,500) 'TIME_INT=',TIME_INT,' HAS BEEN REPLACED BY',
     &              TL(NTX_MIN)
         TIME_INT=TL(NTX_MIN)
         NTSL=NTX_MIN
 1000    CONTINUE
      ENDIF
C
      DO NRL=1,NRFMAX
         TMP(NRL)=F2(NTSL,NRL)
      ENDDO
      CALL SPL1D(RL,TMP,DERIV,U,NRFMAX,ID,IERR)
      IF(IERR.NE.0)
     &     WRITE(6,*) 'XX TRFILE: SPL1D',KFID,': IERR=',IERR
C
      RETURN
 500  FORMAT(' ',A,F9.5,A,F9.5)
      END
C
C     ***********************************************************
C
C           READING DATA FROM UFILES WITH NT INCREMENT
C
C     ***********************************************************
C
      SUBROUTINE TR_UFREAD
C
      INCLUDE 'trcomm.inc'
      COMMON /TRINS1/ INS
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC2/ NTAMAX
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
      COMMON /TMSLC4/ PNBI
C
      TSL=DT*DBLE(NT+1)
      CALL LAGLANGE(TSL,RKAP,TMU1,RKAPU,NTXMAX1,NTUM,IERR)
      CALL LAGLANGE(TSL,RR  ,TMU1,RRU  ,NTXMAX1,NTUM,IERR)
      CALL LAGLANGE(TSL,RA  ,TMU1,RAU  ,NTXMAX1,NTUM,IERR)
      CALL LAGLANGE(TSL,BB  ,TMU1,BBU  ,NTXMAX1,NTUM,IERR)
      CALL LAGLANGE(TSL,PHIA,TMU1,PHIAU,NTXMAX1,NTUM,IERR)
      CALL LAGLANGE(TSL,PNBI,TMU1,PNBIU,NTXMAX1,NTUM,IERR)
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IF(MDLUF.EQ.1) THEN
      DO NR=1,NRMAX
         IF(MDLEQN.EQ.0) THEN
            CALL LAGLANGE(TSL,RNEL ,TMU,RNU(1,NR,1),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,RNDL ,TMU,RNU(1,NR,2),NTXMAX,NTUM,IERR)
            RN(NR,1)=RNEL
            RN(NR,2)=RNDL
            IF(MDNI.NE.0) THEN
               CALL LAGLANGE(TSL,RNIL ,TMU,RNU(1,NR,3),NTXMAX,NTUM,IERR)
               RN(NR,3)=RNIL
            ENDIF
         ENDIF
         IF(INS.NE.0) THEN
            IF(MDLEOI.EQ.1) THEN
               CALL LAGLANGE(TSL,RTDL ,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
               RT(NR,2)=RTDL
            ELSEIF(MDLEOI.EQ.2) THEN
               CALL LAGLANGE(TSL,RTEL ,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
               RT(NR,1)=RTEL
            ENDIF
            CALL LAGLANGE(TSL,RTIL ,TMU,RTU(1,NR,3),NTXMAX,NTUM,IERR)
            IF(INS.EQ.2) RT(NR,3)=RTIL
         ENDIF
         CALL LAGLANGE(TSL,QPL ,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
         QP(NR)=QPL
C
         IF(KUFDEV.EQ.'X') THEN
            IF(PNBI.LT.12.D6) THEN
               PEX(NR,1)=PNBU(1,NR,1)
               PEX(NR,2)=PNBU(1,NR,2)
            ELSE
               PEX(NR,1)=PNBU(2,NR,1)
               PEX(NR,2)=PNBU(2,NR,2)
               IF(NT.EQ.NTAMAX) THEN
                  PEX(NR,1)=PNBU(3,NR,1)
                  PEX(NR,2)=PNBU(3,NR,2)
               ENDIF
            ENDIF
         ELSE
            CALL LAGLANGE(TSL,PNBEL,TMU,PNBU(1,NR,1),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PNBDL,TMU,PNBU(1,NR,2),NTXMAX,NTUM,IERR)
            PEX(NR,1)=PNBEL
            PEX(NR,2)=PNBDL
         ENDIF
C
         CALL LAGLANGE(TSL,PICEL,TMU,PICU(1,NR,1),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,PICDL,TMU,PICU(1,NR,2),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,PECL ,TMU,PECU(1,NR  ),NTXMAX,NTUM,IERR)
         PRF(NR,1)=PICEL+PECL
         PRF(NR,2)=PICDL
         CALL LAGLANGE(TSL,TTRL ,TMU,TTRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,DVRL ,TMU,DVRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,DSRL ,TMU,DSRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,ABRL ,TMU,ABRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,ARRL ,TMU,ARRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,AR1RL,TMU,AR1RHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,AR2RL,TMU,AR2RHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RMJRL,TMU,RMJRHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RMNRL,TMU,RMNRHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RJCBL,TMU,RJCBU  (1,NR),NTXMAX,NTUM,IERR)
         TTRHO(NR)=TTRL
         DVRHO(NR)=DVRL
         DSRHO(NR)=DSRL
         ABRHO(NR)=ABRL
         ARRHO(NR)=ARRL
         AR1RHO(NR)=AR1RL
         AR2RHO(NR)=AR2RL
         RJCB(NR)=RJCBL
         RMJRHO(NR)=RMJRL
         RMNRHO(NR)=RMNRL
         EKAP(NR)=RKAP
         CALL LAGLANGE(TSL,PBML,TMU,PBMU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RNFL,TMU,RNFU(1,NR),NTXMAX,NTUM,IERR)
         PBM(NR)=PBML
         RNF(NR,1)=RNFL
      ENDDO
      CALL TRGFRG
      IF(MDLJQ.EQ.0) THEN
         NR=1
            CALL LAGLANGE(TSL,AJL  ,TMU,AJU  (1,NR),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,AJNBL,TMU,AJNBU(1,NR),NTXMAX,NTUM,IERR)
            AJ(NR)=AJL
            AJNB(NR)=AJNBL
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         DO NR=2,NRMAX
            CALL LAGLANGE(TSL,AJL  ,TMU,AJU  (1,NR),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,AJNBL,TMU,AJNBU(1,NR),NTXMAX,NTUM,IERR)
            AJ(NR)=AJL
            AJNB(NR)=AJNBL
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
      ELSE
         DO NR=1,NRMAX
            RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
      ENDIF
      ELSEIF(MDLUF.EQ.3) THEN
      DO NR=1,NRMAX
         CALL LAGLANGE(TSL,RNEL ,TMU,RNU(1,NR,1),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RNDL ,TMU,RNU(1,NR,2),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RNIL ,TMU,RNU(1,NR,3),NTXMAX,NTUM,IERR)
         RN(NR,1)=RNEL
         RN(NR,2)=RNDL
         RN(NR,3)=RNIL
C         CALL LAGLANGE(TSL,QPL ,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
C         QP(NR)=QPL
         CALL LAGLANGE(TSL,PNBEL,TMU,PNBU(1,NR,1),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,PNBDL,TMU,PNBU(1,NR,2),NTXMAX,NTUM,IERR)
         PEX(NR,1)=PNBEL
         PEX(NR,2)=PNBDL
         CALL LAGLANGE(TSL,PICEL,TMU,PICU(1,NR,1),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,PICDL,TMU,PICU(1,NR,2),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,PECL ,TMU,PECU(1,NR  ),NTXMAX,NTUM,IERR)
         IF(TSL.LE.1.D0) THEN
            PRF(NR,1)=0.D0
            PRF(NR,2)=0.D0
         ELSE
            DO NTL=1,NTXMAX
               IF(TMU(NTL).GT.TSL) THEN
                  NTLL=NTL
                  GOTO 1000
               ENDIF
            ENDDO
            NTLL=NTXMAX
 1000       PRF(NR,1)=PICU(NTLL,NR,1)+PECU(NTLL,NR)
            PRF(NR,2)=PICU(NTLL,NR,2)
         ENDIF
         CALL LAGLANGE(TSL,TTRL ,TMU,TTRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,DVRL ,TMU,DVRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,DSRL ,TMU,DSRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,ABRL ,TMU,ABRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,ARRL ,TMU,ARRHOU (1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,AR1RL,TMU,AR1RHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,AR2RL,TMU,AR2RHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RMJRL,TMU,RMJRHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RMNRL,TMU,RMNRHOU(1,NR),NTXMAX,NTUM,IERR)
         CALL LAGLANGE(TSL,RJCBL,TMU,RJCBU  (1,NR),NTXMAX,NTUM,IERR)
         TTRHO(NR)=TTRL
         DVRHO(NR)=DVRL
         DSRHO(NR)=DSRL
         ABRHO(NR)=ABRL
         ARRHO(NR)=ARRL
         AR1RHO(NR)=AR1RL
         AR2RHO(NR)=AR2RL
         RJCB(NR)=RJCBL
         RMJRHO(NR)=RMJRL
         RMNRHO(NR)=RMNRL
         EKAP(NR)=RKAP
      ENDDO
      ENDIF
C      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
      CALL TRGFRG
C
      IF(MDLUF.EQ.1) THEN
         DO NS=1,2
            CALL LAGLANGE(TSL,PNSL ,TMU,PNSU (1,NS),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PTSL ,TMU,PTSU (1,NS),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PNSAL,TMU,PNSUA(1,NS),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PTSAL,TMU,PTSUA(1,NS),NTXMAX,NTUM,IERR)
            PNS (NS)=PNSL
            PTS (NS)=PTSL
            PNSA(NS)=PNSAL
            PTSA(NS)=PTSAL
         ENDDO
         IF(MDNI.NE.0) THEN 
            NS=3
            CALL LAGLANGE(TSL,PNSL ,TMU,PNSU (1,NS),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PTSL ,TMU,PTSU (1,NS),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PNSAL,TMU,PNSUA(1,NS),NTXMAX,NTUM,IERR)
            CALL LAGLANGE(TSL,PTSAL,TMU,PTSUA(1,NS),NTXMAX,NTUM,IERR)
            PNS (NS)=PNSL
            PTS (NS)=PTSL
            PNSA(NS)=PNSAL
            PTSA(NS)=PTSAL
         ENDIF
         IF(RHOA.NE.1.D0) THEN
            CALL LAGLANGE(TSL,RTDM,TMU,RTU(1,NRMAX,2),NTXMAX,NTUM,IERR)
            DO NR=NRAMAX+1,NRMAX
               PROF    =(1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               CALL LAGLANGE(TSL,RTEL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
               CALL LAGLANGE(TSL,RTDL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
               RT(NR,1)= RTEL
               RT(NR,2)= RTDL
               RT(NR,3)= RTDL
               RT(NR,4)=(RTDL-RTDM)*PROF+RTDM
            ENDDO
            IF(MDNI.EQ.0) THEN
               DO NR=NRAMAX+1,NRMAX
                  CALL LAGLANGE(TSL,RTDL,TMU,RTU(1,NR,2),NTXMAX,NTUM,
     &                          IERR)
                  PROF    =(1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
                  RT(NR,3)=(RTDL-RTDM)*PROF+RTDM
               ENDDO
            ENDIF
         ENDIF
      ENDIF
C
C     *** CALCULATE PZC,PZFE ***
C
      CALL TRZEFF
C
C     *** CALCULATE ANEAVE ***
C
      ANESUM=0.D0
      DO NR=1,NRMAX
         ANESUM=ANESUM+RN(NR,1)*RM(NR)
      ENDDO 
      ANEAVE=ANESUM*2.D0*DR
C
C     *** CALCULATE IMPURITY DENSITY
C                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***
C
      IF(MDLUF.NE.3) THEN
         DO NR=1,NRMAX
            ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC
     &               *1.D-2*RN(NR,1)
            ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE
     &               *1.D-2*RN(NR,1)
            ANI = 0.D0
            DO NS=2,NSM
               ANI=ANI+PZ(NS)*RN(NR,NS)
            ENDDO
            ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
            DILUTE = 1.D0-ANZ/ANI
            DO NS=2,NSM
               RN(NR,NS) = RN(NR,NS)*DILUTE
            ENDDO
         ENDDO
         PNSS(1)=PNS(1)
         DO NS=2,NSM
            PNSS(NS)=PNS(NS)*DILUTE
         ENDDO
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1)=PNSA(1)
            DO NS=2,NSM
               PNSSA(NS)=PNSA(NS)*DILUTE
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ELSE
c$$$         DO NR=1,NRMAX
c$$$            ANC(NR)=1.D0/3.D1*RN(NR,1)
c$$$         ENDDO
         DO NS=1,NSM
            PNSS(NS)=PNS(NS)
         ENDDO
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            DO NS=1,NSM
               PNSSA(NS)=PNSA(NS)*DILUTE
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ENDIF
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
C
      RETURN
      END
C
C     ******
C
      SUBROUTINE TR_UFREAD_S
C
      INCLUDE 'trcomm.inc'
      COMMON /TRINS1/ INS
C
      RKAP=RKAPU(NTSL)
      RR=RRU(NTSL)
      RA=RAU(NTSL)
      BB=BBU(NTSL)
      PHIA=PHIAU(NTSL)
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      DO NR=1,NRMAX
         IF(MDLEQN.EQ.0) THEN
            RN(NR,1)=RNU(1,NR,1)
            RN(NR,2)=RNU(1,NR,2)
            IF(MDNI.NE.0) RN(NR,3)=RNU(1,NR,3)
         ENDIF
         IF(INS.NE.0) THEN
            IF(MDLEOI.EQ.1) THEN
               RT(NR,2)=RTU(1,NR,2)
            ELSEIF(MDLEOI.EQ.2) THEN
               RT(NR,1)=RTU(1,NR,1)
            ENDIF
            IF(INS.EQ.2) RT(NR,3)=RTU(1,NR,3)
         ENDIF
         PEX(NR,1)=PNBU(1,NR,1)
         PEX(NR,2)=PNBU(1,NR,2)
         PRF(NR,1)=PICU(1,NR,1)+PECU(1,NR)
         PRF(NR,2)=PICU(1,NR,2)
         WROT(NR) =WROTU(1,NR)
         IF(MDVTOR.EQ.0) VTOR(NR) =VTORU(1,NR)
         TTRHO(NR)=TTRHOU(1,NR)
         DVRHO(NR)=DVRHOU(1,NR)
         DSRHO(NR)=DSRHOU(1,NR)
         ABRHO(NR)=ABRHOU(1,NR)
         ARRHO(NR)=ARRHOU(1,NR)
         AR1RHO(NR)=AR1RHOU(1,NR)
         AR2RHO(NR)=AR2RHOU(1,NR)
         RJCB(NR)=RJCBU(1,NR)
         RMJRHO(NR)=RMJRHOU(1,NR)
         RMNRHO(NR)=RMNRHOU(1,NR)
         EKAP(NR)=RKAP
         RNF(NR,1)=RNFU(1,NR)
      ENDDO
C      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
      CALL TRGFRG
      DO NR=1,NRMAX
         RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
         BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
      ENDDO
C
C     *** CALCULATE PZC,PZFE ***
C
      CALL TRZEFF
C
C     *** CALCULATE ANEAVE ***
C
      ANESUM=0.D0
      DO NR=1,NRMAX
         ANESUM=ANESUM+RN(NR,1)*RM(NR)
      ENDDO 
      ANEAVE=ANESUM*2.D0*DR
C
C     *** CALCULATE IMPURITY DENSITY
C                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***
C
      DO NR=1,NRMAX
         ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC
     &            *1.D-2*RN(NR,1)
         ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE
     &            *1.D-2*RN(NR,1)
         ANI = 0.D0
         DO NS=2,NSM
            ANI=ANI+PZ(NS)*RN(NR,NS)
         ENDDO
         ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DILUTE = 1.D0-ANZ/ANI
         DO NS=2,NSM
            RN(NR,NS) = RN(NR,NS)*DILUTE
         ENDDO
      ENDDO
      PNSS(1)=PNS(1)
      DO NS=2,NSM
         PNSS(NS)=PNS(NS)*DILUTE
      ENDDO
      PNSS(7)=PNS(7)
      PNSS(8)=PNS(8)
      IF(RHOA.NE.1.D0) THEN
         PNSSA(1)=PNSA(1)
         DO NS=2,NSM
            PNSSA(NS)=PNSA(NS)*DILUTE
         ENDDO
         PNSSA(7)=PNSA(7)
         PNSSA(8)=PNSA(8)
      ENDIF
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
C
      RETURN
      END
C
C     ***********************************************************
C
C           CALCULATING FLUX FROM SOURCE TERM FILES
C
C     ***********************************************************
C
      SUBROUTINE FLUX
C
      INCLUDE 'trcomm.inc'
      DIMENSION SALIL(NRM),SALEL(NRM)
C
      IF(MDLFLX.EQ.0) THEN
         DO NR=1,NRMAX
            DO NS=1,NSM
               RGFLX(NR,NS)=0.D0
            ENDDO
         ENDDO
      ELSE
         DO NR=1,NRMAX
            SALEL(NR)=SNBU(1,NR,1)+SWLU(1,NR)/PZ(2)
            CALL TRSUMD(SALEL,DVRHO,NR,RGESUM)
            RGFLX(NR,1)=RGESUM*DR
            SALIL(NR)=SNBU(1,NR,2)+SWLU(1,NR)
            CALL TRSUMD(SALIL,DVRHO,NR,RGISUM)
            RGFLX(NR,2)=RGISUM*DR
            DO NS=3,NSM
               RGFLX(NR,NS)=0.D0
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C
C   *******************************************
C   **    LAGRANGEAN INTERPOLATION METHOD    **
C   *******************************************
C
C     input:
C
C     SLT     : Interesting Time
C     T(NTM)  : Total Time Data
C     F1(NTM) : Functional Values
C     NTMAX   : Maximum Number of the Time Mesh from UFILE
C     NTM     : Maximum Array Width for Time
C
C     output:
C
C     FOUT    : Interpolated Value
C     IERR    : Error Indicator
C     
C     ***********************************************************
C
      SUBROUTINE LAGLANGE(SLT,FOUT,T,F1,NTMAX,NTM,IERR)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      COMMON /COMEPS/ EPS,IERRL
      DIMENSION T(NTM),F1(NTM)
C
C      M=5
      M=1
      EPS=1.D-5
C
      FOUT=FITLAG(SLT,T,F1,M,NTMAX,NTM)
      IERR=IERRL
C
      RETURN
      END
C
      FUNCTION FITLAG(X,A,F,M,N,NTM)
****************************************************
*                                                  *
*     Compute iterated Lagrange interpolation      *
*     at X based on data given in A(I) and F(I)    *
*                                                  *
*     ==== input data for arguments  ====          *
*                                                  *
*     N...total number of data (A(I),F(I)) given   *
*         in the table                             *
*     A(I)...I-th sampling point                   *
*            A(I) must be given in ascending       *
*            order with respect to I               *
*     F(I)...data at I-th point, i.e. FUNC(A(I))   *
*     M...number of sampling points used in both   *
*         sides of X for interpolation             *
*         M must be less than or equal to 5        *
*                                                  *
*     ==== input data for common block ====        *
*                                                  *
*     EPS...absolute error tolerance               *
*                                                  *
*     ==== output data for common block ====       *
*                                                  *
*     IERR...if IERR = 0 then converged            *
*            if IERR = 1 then not converged        *
*                                                  *
*     ==== work arrays ====                        *
*                                                  *
*     B(J)...A(I) - X                              *
*     V(I,J)...triangular table for interpolation  *
*                                                  *
****************************************************
      IMPLICIT REAL*8 (A-F,H,O-Z)
      DIMENSION V(10,10),B(10)
C
      COMMON /COMEPS/ EPS,IERR
      DIMENSION A(NTM),F(NTM)
C
      IF(M.GT.5) STOP 'M MUST BE LESS THAN 6.'
C
      M2=2*M
C
C     ---- find the nearest sampling point to x ----
C                 by bisection method
C
      IS=1
      IB=N
C
 1    CONTINUE
C
      IM=(IS+IB)/2
      IF(X.LT.A(IM)) THEN
         IB=IM
      ELSEIF(X.GT.A(IM)) THEN
         IS=IM
      ELSE
         FITLAG=F(IM)
         IERR=0
         RETURN
      ENDIF
C
      IF(IS.LT.IB-1) GOTO 1
C
C     ---- set sampling points to use ----
C
      IL=IS-M+1
      IR=IB+M-1
      IF(IL.LT.1) THEN
         IL=1
      ELSEIF(IR.GT.N) THEN
         IL=N-M2+1
      ENDIF
C
      DO I=1,M2
         B(I)=A(I+IL-1)-X
         V(1,I)=F(I+IL-1)
      ENDDO
C
C     ---- sort ABS(B(I)) in ascending order ----
C
      DO I=2,M2
         K=I-1
         DO J=1,M2
            IF(ABS(B(J)).LT.ABS(B(K))) K=J
         ENDDO
         IF(K.NE.I-1) THEN
            BV=B(K)
            B(K)=B(I-1)
            B(I-1)=BV
            VV=V(1,K)
            V(1,K)=V(1,I-1)
            V(1,I-1)=VV
         ENDIF
      ENDDO
C
C     ---- compute iterated interpolation ----
C
      DO I=2,M2
C
         DO J=2,I
C            write(6,*) I,J-1,B(I),B(J-1)
            V(J,I)=(V(J-1,J-1)*B(I)
     &            - V(J-1,I  )*B(J-1))/(B(I)-B(J-1))
         ENDDO
C
         IF(ABS(V(I,I)-V(I-1,I-1)).LE.EPS) THEN
C
C     ---- converged ----
C
            FITLAG=V(I,I)
            IERR=0
            RETURN
         ENDIF
C
      ENDDO
C
C     ---- not converged ----
C
      FITLAG=V(M2,M2)
      IERR=1
      RETURN
C
      END
