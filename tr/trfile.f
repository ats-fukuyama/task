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
         GRM(NR)  =GUCLIP(RM(NR))
         GRG(NR+1)=GUCLIP(RG(NR))
         QP(NR)   =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
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
C           CONTROL ROUTINE FOR UFILE READING
C
C     ***********************************************************
C
      SUBROUTINE TR_UFILE_CONTROL(NSW)
C
      INCLUDE 'trcomm.inc'
      CHARACTER KDIRX*80
      CHARACTER KDIRR2*80
      CHARACTER KFILE*80,KFID*20
      COMMON /TRUFC1/ KDIRX
      COMMON /TRUFC2/ IKNDEV,IKNDCG
      COMMON /TRUFC3/ NM2CHK
      LOGICAL DIR,LEX
C
      NRAMAX=INT(RHOA*NRMAX)
      DR = 1.D0/DBLE(NRMAX)
C
      IF(NSW.NE.0) THEN
      CALL KTRIM(KUFDEV,IKNDEV)
      CALL KTRIM(KUFDCG,IKNDCG)
C
      KDIRX='../../tr.new/data/'//KUFDEV(1:IKNDEV)//'/'
     &                          //KUFDCG(1:IKNDCG)//'/in/'
C      KDIRX='../../../profiledb/profile_data/'//KUFDEV(1:IKNDEV)//'/'
C     &                          //KUFDCG(1:IKNDCG)//'/in/'
      INQUIRE(FILE=KDIRX,EXIST=DIR,ERR=9000)
      IF(DIR.NEQV..TRUE.) THEN
         WRITE(6,'(A25,A34,A17)') 
     &        '## DESIGNATED DIRECTORY( ',KDIRX,' ) DOES NOT EXIST!'
         STOP
      ENDIF
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
      IF(LEX.EQV..TRUE.) THEN
         IF(KUFDEV.EQ.'tftr'.OR.KUFDEV.EQ.'d3d') THEN
            MDSLCT=MDSLCT+1
            NM2CHK=1
         ELSE
            GOTO 1000
         ENDIF
      ELSE
 1000    CONTINUE
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
      IF(MDSLCT.EQ.0) THEN
         MDNI=0
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
      COMMON /PRETREAT1/ RUF(NRMU),TMU(NTURM),F1(NTURM),F2(NRMU,NTURM)
      COMMON /TRUFC3/ NM2CHK
      COMMON /TRUFC4/ NREMAX(2),GRE(NRM,2)
      COMMON /TRSVUC/ RTEXU(NRMP,NTUM),   RTIXU(NRMP,NTUM),
     &                RNEXU(NRMP,NTUM)
      COMMON /TRSVUD/ RTEXEU(NRMP,NTUM),  RTIXEU(NRMP,NTUM),
     &                RNEXEU(NRMP,NTUM)
      DIMENSION FAS(NRMP),TMUS(NTURM)
      CHARACTER KFID*10
C
      MDRGEO=0
      MDAMIN=0
      MDIP=0
      MDBT=0
      MDKAPPA=0
C
C     *** 1D VALUE ***
C
      KFID='RGEO'
      CALL UFREAD_TIME(KFID,TMU,F1,NTXMAX,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDRGEO=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX
            RRU(NTX)=F1(NTX)
         ENDDO
      ENDIF
C
      KFID='AMIN'
      CALL UFREAD_TIME(KFID,TMU,F1,NTXMAX,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDAMIN=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX
            RAU(NTX)=F1(NTX)
         ENDDO
      ENDIF
C
      KFID='IP'
      CALL UFREAD_TIME(KFID,TMU,F1,NTXMAX,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDIP=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX
            RIPU(NTX)=ABS(F1(NTX)*1.D-6)
         ENDDO
      ENDIF
C
      KFID='BT'
      CALL UFREAD_TIME(KFID,TMU,F1,NTXMAX,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDBT=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX
            BBU(NTX)=ABS(F1(NTX))
         ENDDO
      ENDIF
C
      KFID='KAPPA'
      CALL UFREAD_TIME(KFID,TMU,F1,NTXMAX,MDCHK,IERR)
      IF(IERR.EQ.1.AND.KUFDEV.EQ.'jt60u') THEN
         MDKAPPA=IERR
         IERR=0
      ELSE
         DO NTX=1,NTXMAX
            RKAPU(NTX)=F1(NTX)
         ENDDO
      ENDIF
C
      IF(IERR.NE.0) STOP 'SOME 1D UFILES DO NOT EXIST.'
      DO NTX=1,NTXMAX
         TMUS(NTX)=TMU(NTX)
      ENDDO
      NTXSMAX=NTXMAX
C
C     *** 2D VALUE ***
C
      AMP=1.D-3
      KFID='TE'
      CALL UF2DSP(KFID,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RTU(NR,1,1)=FAS(NR)
         RT(NR,1)=FAS(NR)
      ENDDO
      PT(1)=RT(1,1)
      PTS(1)=PTMP1
      IF(RHOA.NE.1.D0) PTSA(1)=PTMP2
C
      KFID='TEXP'
      CALL UFREAD2_ERROR(KFID,RUF,TMU,RTEXU,NRFMAX,NTXMAX,MDCHK,IERR)
      NREMAX(1)=NRFMAX
      DO NR=1,NRFMAX
         GRE(NR,1)=GUCLIP(RUF(NR))
      ENDDO
      KFID='TEXPEB'
      CALL UFREAD2_ERROR(KFID,RUF,TMU,RTEXEU,NRFMAX,NTXMAX,MDCHK,IERR)
C     
      KFID='TI'
      CALL UF2DSP(KFID,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RTU(NR,2,1)=FAS(NR)
         RTU(NR,3,1)=FAS(NR)
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
      CALL UFREAD2_ERROR(KFID,RUF,TMU,RTIXU,NRFMAX,NTXMAX,MDCHK,IERR)
      KFID='TIXPEB'
      CALL UFREAD2_ERROR(KFID,RUF,TMU,RTIXEU,NRFMAX,NTXMAX,MDCHK,IERR)
C
      AMP=1.D-20
      KFID='NE'
      CALL UF2DSP(KFID,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RNU(NR,1,1)=FAS(NR)
         RNU(NR,2,1)=FAS(NR)
      ENDDO
      PN(1)  =RNU(1,1,1)
      PN(2)  =RNU(1,2,1)-2.D-7
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
      CALL UFREAD2_ERROR(KFID,RUF,TMU,RNEXU,NRFMAX,NTXMAX,MDCHK,IERR)
      NREMAX(2)=NRFMAX
      DO NR=1,NRFMAX
         GRE(NR,2)=GUCLIP(RUF(NR))
      ENDDO
      KFID='NEXPEB'
      CALL UFREAD2_ERROR(KFID,RUF,TMU,RNEXEU,NRFMAX,NTXMAX,MDCHK,IERR)
C
      IF(MDNI.EQ.1) THEN !!!
      IF(NM2CHK.EQ.1) THEN
         KFID='NM2'
      ELSE
         KFID='NM1'
      ENDIF
      CALL UF2DSP(KFID,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RNU(NR,2,1)=FAS(NR)
         RNU(NR,3,1)=(RNU(NR,1,1)-PZ(2)*RNU(NR,2,1))/PZ(3)-1.D-7
         ZEFFU(NR,1)=(PZ(2)**2*RNU(NR,2,1)+PZ(3)**2*RNU(NR,3,1))
     &              / RNU(NR,1,1)
      ENDDO
      PN(2) =RNU(1,2,1)
      PN(3) =RNU(1,3,1)-1.D-7
      PN(4) =1.D-7
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
      CALL UF2DSP(KFID,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         ZEFFU(NR,1)=FAS(NR)
         RNU(NR,2,1)=(PZ(3)-FAS(NR))/(PZ(2)*(PZ(3)-PZ(2)))*RNU(NR,1,1)
         RNU(NR,3,1)=(PZ(2)-FAS(NR))/(PZ(3)*(PZ(2)-PZ(3)))*RNU(NR,1,1)
      ENDDO
      PN(2)=RNU(1,2,1)
      PN(3)=RNU(1,3,1)-1.D-7
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
      CALL UF2DSP(KFID,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,IERR)
      DO NR=1,NRMAX
         RNU(NR,2,1)=(RNU(NR,1,1)-PZ(3)*FAS(NR))/PZ(2)
         RNU(NR,3,1)=FAS(NR)
         ZEFFU(NR,1)=PZ(2)+PZ(3)*(PZ(3)-PZ(2))*(FAS(NR)/RNU(NR,1,1))
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
      AMP=1.D0
      KFID='Q'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,1,IERR)
      DO NR=1,NRMAX
         QPU(NR,1)=FAS(NR)
      ENDDO
c$$$      IF(KUFDEV.EQ.'jt60u') THEN
c$$$         DO NR=NRMAX-1,NRMAX
c$$$            QPU(NR,1)=FEDG(DBLE(NR)*DR,DBLE(NR-1)*DR,DBLE(NR-2)*DR,
c$$$     &                     QPU(NR-1,1),QPU(NR-2,1))
c$$$         ENDDO
c$$$         WRITE(6,*) KFID,'IS MODIFIED.'
c$$$      ENDIF
C
      KFID='CURTOT'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         AJU(NR,1)=FAS(NR)
      ENDDO
C
      KFID='BPOL'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,1,IERR)
      DO NR=1,NRMAX
         BPU(NR,1)=FAS(NR)
      ENDDO
C
      KFID='QNBIE'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PNBU(NR,1,1)=0.D0
         ELSE
            PNBU(NR,1,1)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QNBII'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PNBU(NR,2,1)=0.D0
         ELSE
            PNBU(NR,2,1)=FAS(NR)
         ENDIF
      ENDDO
C
      DO NR=1,NRMAX
         PICU(NR,1,1)=0.D0
         PICU(NR,2,1)=0.D0
         PECU(NR,  1)=0.D0
      ENDDO
      KFID='QICRHE'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PICU(NR,1,1)=0.D0
         ELSE
            PICU(NR,1,1)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QICRHI'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         IF(FAS(NR).LT.0.D0) THEN
            PICU(NR,2,1)=0.D0
         ELSE
            PICU(NR,2,1)=FAS(NR)
         ENDIF
      ENDDO
C
      KFID='QRAD'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         PRLU(NR,1)=FAS(NR)
      ENDDO
C
      KFID='SNBIE'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         SNBU(NR,1,1)=FAS(NR)
      ENDDO
C
      KFID='SNBII'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         SNBU(NR,2,1)=FAS(NR)
      ENDDO
C
C     *** GEOMETORY FACTORS ***
C
      KFID='RMAJOR'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         RMJRHOU(NR,1)=FAS(NR)
C         ARRHOU(NR,1)=1.D0/RMJRHOU(NR,1)**2
         ARRHOU(NR,1)=1.D0/RRU(1)**2
         TTRHOU(NR,1)=BBU(1)*RRU(1)
      ENDDO
C
      KFID='RMINOR'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         RMNRHOU(NR,1)=FAS(NR)
      ENDDO
C
      KFID='GRHO1'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         AR1RHOU(NR,1)=FAS(NR)
      ENDDO
C
      KFID='GRHO2'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         AR2RHOU(NR,1)=FAS(NR)
         ABRHOU(NR,1)=AR2RHOU(NR,1)*ARRHOU(NR,1)
      ENDDO
c$$$      IF(KUFDEV.EQ.'jt60u') THEN
c$$$         DO NR=NRMAX-2,NRMAX
c$$$            AR2RHOU(NR,1)=FEDG(DBLE(NR)*DR,DBLE(NR-1)*DR,DBLE(NR-2)*DR,
c$$$     &                         AR2RHOU(NR-1,1),AR2RHOU(NR-2,1))
c$$$            ABRHOU(NR,1)=AR2RHOU(NR,1)*ARRHOU(NR,1)
c$$$         ENDDO
c$$$         WRITE(6,*) KFID,'IS MODIFIED.'
c$$$      ENDIF
C
      KFID='SURF'
      CALL UF2DS(KFID,DR,TMU,FAS,AMP,NRMAX,0,IERR)
      DO NR=1,NRMAX
         DVRHOU(NR,1)=FAS(NR)/AR1RHOU(NR,1)
C         DSRHOU(NR,1)=DVRHOU(NR,1)/(2.D0*PI*RMJRHOU(NR,1))
         DSRHOU(NR,1)=DVRHOU(NR,1)/(2.D0*PI*RRU(1))
      ENDDO
C
C     *****
C
      IF(NTXSMAX.EQ.0) THEN
         NTXSMAX=NTXMAX
         DO NTX=1,NTXMAX
            TMUS(NTX)=TMU(NTX)
         ENDDO
      ENDIF
C
      MDSUM=MDRGEO+MDAMIN+MDIP+MDBT+MDKAPPA
      IF(MDSUM.NE.0) THEN
         DO NTX=1,NTXMAX
           IF(MDRGEO .NE.0) RRU  (NTX)=RR
           IF(MDAMIN .NE.0) RAU  (NTX)=RA
           IF(MDIP   .NE.0) RIPU (NTX)=RIPS
           IF(MDBT   .NE.0) BBU  (NTX)=BB
           IF(MDKAPPA.NE.0) RKAPU(NTX)=RKAP
         ENDDO
      ENDIF
C
      DO NTX=1,NTXSMAX
         IF(ABS(TMUS(NTX)-TMU(1)).LE.1.D-6) THEN
            RR   = RRU(NTX)
            RA   = RAU(NTX)
            RIPS = RIPU(NTX)
            BB   = BBU(NTX)
            RKAP = RKAPU(NTX)
            GOTO 2000
         ENDIF
      ENDDO
      STOP 'UFILE HAS AN ERROR!'
 2000 CONTINUE
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
      COMMON /PRETREAT1/ RUF(NRMU),TMU(NTURM),F1(NTURM),F2(NRMU,NTURM)
      COMMON /PRETREAT2/ NTAMAX
      COMMON /TRUFC3/ NM2CHK
      DIMENSION FAT(NRMP,NTUM),PV(NTUM),PVA(NTUM)
      CHARACTER KFID*10
C
      ICK=0
      TMUMAX=0.D0
C
C     *** 1D VALUE ***
C
      KFID='RGEO'
      CALL UF1D(KFID,DT,TMU,RRU  ,NTAMAX,TMUMAX,ICK,IERR)
      KFID='AMIN'
      CALL UF1D(KFID,DT,TMU,RAU  ,NTAMAX,TMUMAX,ICK,IERR)
      KFID='IP'
      CALL UF1D(KFID,DT,TMU,RIPU ,NTAMAX,TMUMAX,ICK,IERR)
      KFID='BT'
      CALL UF1D(KFID,DT,TMU,BBU  ,NTAMAX,TMUMAX,ICK,IERR)
      KFID='KAPPA'
      CALL UF1D(KFID,DT,TMU,RKAPU,NTAMAX,TMUMAX,ICK,IERR)
      IF(IERR.EQ.1) THEN
         DO NTA=1,NTAMAX
            RKAPU(NTA)=1.D0
         ENDDO
      ENDIF
C
C     *** 2D VALUE ***
C
      AMP=1.D-3
      KFID='TE'
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,NRAMAX,
     &            NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NR=1,NRMAX
         DO NTA=1,NTAMAX
            RTU(NR,1,NTA)=FAT(NR,NTA)
         ENDDO
         RT(NR,1)=RTU(NR,1,1)
      ENDDO
      PT(1)=RT(1,1)
      DO NTA=1,NTAMAX
         PTSU(1,NTA)=PV(NTA)
      ENDDO
      PTS(1)=PTSU(1,1)
      IF(RHOA.NE.1.D0) THEN
         DO NTA=1,NTAMAX
            PTSUA(1,NTA)=PVA(NTA)
         ENDDO
         PTSA(1)=PTSUA(1,1)
      ENDIF
C
      KFID='TI'
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,NRAMAX,
     &            NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NR=1,NRMAX
         DO NS=2,3
            DO NTA=1,NTAMAX
               RTU(NR,NS,NTA)=FAT(NR,NTA)
            ENDDO
            RT(NR,NS)=RTU(NR,NS,1)
         ENDDO
      ENDDO
      PT(2)=RT(1,2)
      PT(3)=RT(1,3)
      PT(4)=RT(1,2)
      DO NS=2,4
         DO NTA=1,NTAMAX
            PTSU (NS,NTA)=PV(NTA)
         ENDDO
         PTS (NS)=PTSU(NS,1)
      ENDDO
      IF(RHOA.NE.1.D0) THEN
         DO NS=2,4
            DO NTA=1,NTAMAX
               PTSUA(NS,NTA)=PVA(NTA)
            ENDDO
            PTSA(NS)=PTSUA(NS,1)
         ENDDO
      ENDIF
C
      AMP=1.D-20
      KFID='NE'
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,NRAMAX,
     &            NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            RNU(NR,1,NTA)=FAT(NR,NTA)
            RNU(NR,2,NTA)=FAT(NR,NTA)-2.D-7
            RNU(NR,3,NTA)=1.D-7
            RNU(NR,4,NTA)=1.D-7
         ENDDO
      ENDDO
      PN(1)=RNU(1,1,1)
      PN(2)=RNU(1,2,1)-2.D-7
      PN(3)=1.D-7
      PN(4)=1.D-7
      DO NTA=1,NTAMAX
         PNSU(1,NTA)=PV(NTA)
         PNSU(2,NTA)=PV(NTA)-2.D-8
         PNSU(3,NTA)=1.D-8
         PNSU(4,NTA)=1.D-8
      ENDDO
      PNS(1)=PNSU(1,1)
      PNS(2)=PNSU(2,1)
      PNS(3)=PNSU(3,1)
      PNS(4)=PNSU(4,1)
      IF(RHOA.NE.1.D0) THEN
         DO NTA=1,NTAMAX
            PNSUA(1,NTA)=PVA(NTA)
            PNSUA(2,NTA)=PVA(NTA)-2.D-8
            PNSUA(3,NTA)=1.D-8
            PNSUA(4,NTA)=1.D-8
         ENDDO
         PNSA(1)=PNSUA(1,1)
         PNSA(2)=PNSUA(2,1)
         PNSA(3)=PNSUA(3,1)
         PNSA(4)=PNSUA(4,1)
      ENDIF
C
      IF(MDNI.EQ.1) THEN !!!
      IF(NM2CHK.EQ.1) THEN
         KFID='NM2'
      ELSE
         KFID='NM1'
      ENDIF
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,NRAMAX,
     &            NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            RNU(NR,2,NTA)=FAT(NR,NTA)
            RNU(NR,3,NTA)=(RNU(NR,1,NTA)-PZ(2)*RNU(NR,2,NTA))/PZ(3) 
            ZEFFU(NR,NTA)=(PZ(2)**2*RNU(NR,2,NTA)
     &                    +PZ(3)**2*RNU(NR,3,NTA))/RNU(NR,1,NTA)
         ENDDO
      ENDDO
      PN(2)=RNU(1,2,1)
      PN(3)=RNU(1,3,1)-1.D-7
      PN(4)=1.D-7
      DO NTA=1,NTAMAX
         PNSU(2,NTA)=PV(NTA)
         PNSU(3,NTA)=(PNSU(1,NTA)-PZ(2)*PNSU(2,NTA))/PZ(3)-1.D-8
         PNSU(4,NTA)=1.D-8
      ENDDO
      PNS(2)=PNSU(2,1)
      PNS(3)=PNSU(3,1)
      PNS(4)=PNSU(4,1)
      IF(RHOA.NE.1.D0) THEN
         DO NTA=1,NTAMAX
            PNSUA(2,NTA)=PVA(NTA)
            PNSUA(3,NTA)=(PNSUA(1,NTA)-PZ(2)*PNSUA(2,NTA))/PZ(3)-1.D-8
            PNSUA(4,NTA)=1.D-8
         ENDDO
         PNSA(2)=PNSUA(2,1)
         PNSA(3)=PNSUA(3,1)
         PNSA(4)=PNSUA(4,1)
      ENDIF
C
      ELSEIF(MDNI.EQ.2) THEN !!!
C
      AMP=1.D0
      KFID='ZEFFR'
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,NRAMAX,
     &            NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            ZEFFU(NR,NTA)=FAT(NR,NTA)
            RNU(NR,2,NTA)=(PZ(3)-FAT(NR,NTA))/(PZ(2)*(PZ(3)-PZ(2)))
     &                   *RNU(NR,1,NTA)
            RNU(NR,3,NTA)=(PZ(2)-FAT(NR,NTA))/(PZ(3)*(PZ(2)-PZ(3)))
     &                   *RNU(NR,1,NTA)
            IF(RNU(NR,3,NTA).LT.0.D0) STOP 'NIMP IS BELOW ZERO.'
         ENDDO
      ENDDO
      PN(2)=RNU(1,2,1)
      PN(3)=RNU(1,3,1)-1.D-7
      PN(4)=1.D-7
      DO NTA=1,NTAMAX
         PNSU(2,NTA)=(PZ(3)-PV(NTA))/(PZ(2)*(PZ(3)-PZ(2)))*PNSU(1,NTA)
         PNSU(3,NTA)=(PZ(2)-PV(NTA))/(PZ(3)*(PZ(2)-PZ(3)))*PNSU(1,NTA)
     &              -1.D-8
         PNSU(4,NTA)=1.D-8
      ENDDO
      PNS(2)=PNSU(2,1)
      PNS(3)=PNSU(3,1)
      PNS(4)=PNSU(4,1)
      IF(RHOA.NE.1.D0) THEN
         DO NTA=1,NTAMAX
            PNSUA(2,NTA)=(PZ(3)-PVA(NTA))/(PZ(2)*(PZ(3)-PZ(2)))
     &                  *PNSUA(1,NTA)
            PNSUA(3,NTA)=(PZ(2)-PVA(NTA))/(PZ(3)*(PZ(2)-PZ(3)))
     &                  *PNSUA(1,NTA)-1.D-8
            PNSUA(4,NTA)=1.D-8
         ENDDO
         PNSA(2)=PNSUA(2,1)
         PNSA(3)=PNSUA(3,1)
         PNSA(4)=PNSUA(4,1)
      ENDIF
C
      ELSEIF(MDNI.EQ.3) THEN !!!
C
      KFID='NIMP'
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,NRAMAX,
     &            NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            RNU(NR,2,NTA)=(RNU(NR,1,NTA)-PZ(3)*FAT(NR,NTA))/PZ(2)
            RNU(NR,3,NTA)=FAT(NR,NTA)-1.D-7
            ZEFFU(NR,NTA)=PZ(2)+PZ(3)
     &                  *(PZ(3)-PZ(2))*(FAT(NR,NTA)/RNU(NR,1,NTA))
         ENDDO
      ENDDO
      PN(2)=RNU(1,2,1)
      PN(3)=RNU(1,3,1)-1.D-7
      PN(4)=1.D-7
      DO NTA=1,NTAMAX
         PNSU(2,NTA)=(PNSU(1,NTA)-PZ(3)*PV(NTA))/PZ(2)
         PNSU(3,NTA)=PV(NTA)-1.D-8
         PNSU(4,NTA)=1.D-8
      ENDDO
      PNS(2)=PNSU(2,1)
      PNS(3)=PNSU(3,1)
      PNS(4)=PNSU(4,1)
      IF(RHOA.NE.1.D0) THEN
         DO NTA=1,NTAMAX
            PNSUA(2,NTA)=(PNSUA(1,NTA)-PZ(3)*PVA(NTA))/PZ(2)
            PNSUA(3,NTA)=PVA(NTA)-1.D-8
            PNSUA(4,NTA)=1.D-8
         ENDDO
         PNSA(2)=PNSUA(2,1)
         PNSA(3)=PNSUA(3,1)
         PNSA(4)=PNSUA(4,1)
      ENDIF
      ENDIF !!!
C
      AMP=1.D0
      KFID='Q'
      CALL UF2DT(KFID,DR,DT,TMU,QPU ,AMP,NTAMAX,NRMAX,TMUMAX,1,0,
     &           ICK,IERR)
      KFID='CURTOT'
      CALL UF2DT(KFID,DR,DT,TMU,AJU ,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &           ICK,IERR)
      KFID='BPOL'
      CALL UF2DT(KFID,DR,DT,TMU,BPU ,AMP,NTAMAX,NRMAX,TMUMAX,1,1,
     &     ICK,IERR)
      KFID='QNBIE'
      CALL UF2DT(KFID,DR,DT,TMU,FAT ,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PNBU(NR,1,NTA)=0.D0
            ELSE
               PNBU(NR,1,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QNBII'
      CALL UF2DT(KFID,DR,DT,TMU,FAT ,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PNBU(NR,2,NTA)=0.D0
            ELSE
               PNBU(NR,2,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            PICU(NR,1,NTA)=0.D0
            PICU(NR,2,NTA)=0.D0
            PECU(NR,  NTA)=0.D0
         ENDDO
      ENDDO
      KFID='QICRHE'
      CALL UF2DT(KFID,DR,DT,TMU,FAT ,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PICU(NR,1,NTA)=0.D0
            ELSE
               PICU(NR,1,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QICRHI'
      CALL UF2DT(KFID,DR,DT,TMU,FAT ,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PICU(NR,2,NTA)=0.D0
            ELSE
               PICU(NR,2,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QECH'
      CALL UF2DT(KFID,DR,DT,TMU,FAT ,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PECU(NR,NTA)=0.D0
            ELSE
               PECU(NR,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QRAD'
      CALL UF2DT(KFID,DR,DT,TMU,PRLU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
C     *** GEOMETORY FACTORS ***
C
      KFID='RMAJOR'
      CALL UF2DT(KFID,DR,DT,TMU,RMJRHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
C            ARRHOU(NR,NTA)=1.D0/RMJRHOU(NR,NTA)**2
C            EPS=RMN*RAU(NTA)/RRU(NTA)
C            ARRHOU(NR,NTA)=1.D0/(RRU(NTA)**2*(1.D0-EPS**2)**1.5)
C            ARRHOU(NR,NTA)=(1.D0+1.5D0*EPS**2)/RRU(NTA)**2
            ARRHOU(NR,NTA)=1.D0/RRU(NTA)**2
            TTRHOU(NR,NTA)=RRU(NTA)*BBU(NTA)
         ENDDO
         IF(KUFDEV.EQ.'jt60u') RRU(NTA)=RMJRHOU(1,NTA)
      ENDDO
C
      KFID='RMINOR'
      CALL UF2DT(KFID,DR,DT,TMU,RMNRHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,1,
     &     ICK,IERR)
      IF(KUFDEV.EQ.'jt60u') THEN
         DO NTA=1,NTAMAX
            RAU(NTA)=RMNRHOU(NRMAX,NTA)
         ENDDO
      ENDIF
C
      KFID='GRHO1'
      CALL UF2DT(KFID,DR,DT,TMU,AR1RHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      IF(KUFDEV.EQ.'jt60u') THEN
         DO NTA=1,NTAMAX
            AR1RHOU(NR,NTA)=1.D0/RMNRHOU(NRMAX,NTA)
         ENDDO
      ENDIF
C
      KFID='GRHO2'
      CALL UF2DT(KFID,DR,DT,TMU,AR2RHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            ABRHOU(NR,NTA)=AR2RHOU(NR,NTA)*ARRHOU(NR,NTA)
            IF(KUFDEV.EQ.'jt60u') THEN
               AR2RHOU(NR,NTA)=1.D0/RMNRHOU(NRMAX,NTA)**2
               ABRHOU(NR,NTA)=AR2RHOU(NR,NTA)*ARRHOU(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='SURF'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,1,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            DVRHOU(NR,NTA)=FAT(NR,NTA)/AR1RHOU(NR,NTA)
C            DSRHOU(NR,NTA)=DVRHOU(NR,NTA)/(2.D0*PI*RMJRHOU(NR,NTA))
            DSRHOU(NR,NTA)=DVRHOU(NR,NTA)/(2.D0*PI*RRU(NTA))
         ENDDO
      ENDDO
C
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
      COMMON /PRETREAT1/ RUF(NRMU),TMU(NTURM),F1(NTURM),F2(NRMU,NTURM)
      COMMON /PRETREAT2/ NTAMAX
      DIMENSION FAT(NRMP,NTUM),PV(NTUM),PVA(NTUM)
      CHARACTER KFID*10
C
      ICK=0
      TMUMAX=0.D0
C
C     *** 2D VALUE ***
C
      AMP=1.D-3
      KFID='TE'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &           ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            RTU(NR,1,NTA)=FAT(NR,NTA)
         ENDDO
      ENDDO
C
      KFID='TI'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &           ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            RTU(NR,2,NTA)=FAT(NR,NTA)
         ENDDO
      ENDDO
C
      AMP=1.D-20
      KFID='NE'
      CALL UF2DTP(KFID,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,RHOA,
     &            NRAMAX,NRMAX,TMUMAX,0,0,ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            RNU(NR,1,NTA)=FAT(NR,NTA)
            RNU(NR,2,NTA)=FAT(NR,NTA)
         ENDDO
         PNSU (1,NTA)=PV (NTA)
         PNSU (2,NTA)=PV (NTA)
      ENDDO
      IF(RHOA.NE.1.D0) THEN
         DO NTA=1,NTAMAX
            PNSUA(1,NTA)=PVA(NTA)
            PNSUA(2,NTA)=PVA(NTA)
         ENDDO
      ENDIF
C
      AMP=1.D0
      KFID='Q'
      CALL UF2DT(KFID,DR,DT,TMU,QPU,AMP,NTAMAX,NRMAX,TMUMAX,1,0,
     &           ICK,IERR)
      KFID='CURTOT'
      CALL UF2DT(KFID,DR,DT,TMU,AJU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &           ICK,IERR)
      IF(KUFDEV.EQ.'X'.AND.KUFDCG.EQ.'14') THEN
         KFID='CURBS'
         CALL UF2DT(KFID,DR,DT,TMU,AJBSU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &              ICK,IERR)
      ENDIF
C
      KFID='QNBIE'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &           ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PNBU(NR,1,NTA)=0.D0
            ELSE
               PNBU(NR,1,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
      KFID='QNBII'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PNBU(NR,2,NTA)=0.D0
            ELSE
               PNBU(NR,2,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            PICU(NR,1,NTA)=0.D0
            PICU(NR,2,NTA)=0.D0
            PECU(NR,  NTA)=0.D0
         ENDDO
      ENDDO
      KFID='QICRHE'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PICU(NR,1,NTA)=0.D0
            ELSE
               PICU(NR,1,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
      KFID='QICRHI'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PICU(NR,2,NTA)=0.D0
            ELSE
               PICU(NR,2,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QECH'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            IF(FAT(NR,NTA).LT.0.D0) THEN
               PECU(NR,NTA)=0.D0
            ELSE
               PECU(NR,NTA)=FAT(NR,NTA)
            ENDIF
         ENDDO
      ENDDO
C
      KFID='QRAD'
      CALL UF2DT(KFID,DR,DT,TMU,PRLU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
C     *** GEOMETORY FACTORS ***
C
      KFID='RMAJOR'
      CALL UF2DT(KFID,DR,DT,TMU,RMJRHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
      KFID='RMINOR'
      CALL UF2DT(KFID,DR,DT,TMU,RMNRHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
      KFID='GRHO1'
      CALL UF2DT(KFID,DR,DT,TMU,AR1RHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
      KFID='GRHO2'
      CALL UF2DT(KFID,DR,DT,TMU,AR2RHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
      KFID='VRO'
      CALL UF2DT(KFID,DR,DT,TMU,DVRHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
C            DSRHOU(NR,NTA)=DVRHOU(NR,NTA)/(2.D0*PI*RMJRHOU(NR,NTA))
            DSRHOU(NR,NTA)=DVRHOU(NR,NTA)/(2.D0*PI*RR)
         ENDDO
      ENDDO
C
      KFID='AAT'
      CALL UF2DT(KFID,DR,DT,TMU,ARRHOU,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
C
      KFID='HDT'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            TTRHOU(NR,NTA)=FAT(NR,NTA)/ARRHOU(NR,NTA)
         ENDDO
      ENDDO
C
      KFID='CKT'
      CALL UF2DT(KFID,DR,DT,TMU,FAT,AMP,NTAMAX,NRMAX,TMUMAX,0,0,
     &     ICK,IERR)
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
            ABRHOU(NR,NTA)=FAT(NR,NTA)/DVRHOU(NR,NTA)**2
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         RG(NR)    =DBLE(NR)*DR
         EPSRHO(NR)=RA*RG(NR)/RR
         EKAPPA(NR)=RKAP
      ENDDO
C
C     *** 1D VALUE ***
C
      DO NTA=1,NTAMAX
         RRU(NTA)=RR
         RAU(NTA)=RA
         RIPU(NTA)=RIPS
         BBU(NTA)=BB
         RKAPU(NTA)=RKAP
      ENDDO
C
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
      SUBROUTINE UF1D(KFID,DT,TL,FOUT,NTLMAX,TLMAX,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION TL(NTURM),F1(NTURM),FOUT(NTUM),U(4,NTURM)
      CHARACTER KFID*10
C
      CALL UFREAD_TIME(KFID,TL,F1,NTXMAX,MDCHK,IERR)
      CALL PRETREATMENT1(KFID,TL,F1,U,NTXMAX,NTLMAX,TLMAX,ICK,IERR)
      DO NTL=1,NTLMAX
         TLN=DT*DBLE(NTL)
         CALL SPL1DF(TLN,F0,TL,U,NTXMAX,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,600) "XX TRFILE: SPL1DF",KFID,": IERR=",IERR
         FOUT(NTL)=F0
      ENDDO
C
 600  FORMAT(' ',A17,A10,A7,I2)
      RETURN
      END
C
      SUBROUTINE UF1DG(KFID,GTL,TL,FOUT,AMP,NINMAX,TLMAX,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION TL(NTURM),F1(NTURM),FOUT(NTUM),U(4,NTURM),GTL(NTM)
      CHARACTER KFID*10
C
      CALL UFREAD_TIME(KFID,TL,F1,NTXMAX,MDCHK,IERR)
      CALL PRETREATMENT1(KFID,TL,F1,U,NTXMAX,NTLMAX,TLMAX,ICK,IERR)
      DO NIN=1,NINMAX
         TLN=DBLE(GTL(NIN))
         CALL SPL1DF(TLN,F0,TL,U,NTXMAX,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,600) "XX TRFILE: SPL1DF",KFID,": IERR=",IERR
         FOUT(NIN)=F0
      ENDDO
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
 600  FORMAT(' ',A17,A10,A7,I2)
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
C
C     output:
C
C     TL(NTURM)   : Total Time Data (The Number of DT)
C     FOUT(NRMP)  : Function Values
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF2DS(KFID,DR,TL,FOUT,AMP,NRMAX,NSW,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION RL(NRMU),TL(NTURM),F2(NRMU,NTURM),FOUT(NRMP)
      DIMENSION U(4,NRMU)
      CHARACTER KFID*10
C
      CALL UFREAD2_TIME(KFID,RL,TL,F2,NRLMAX,NTXMAX,MDCHK,IERR)
      CALL PRETREATMENT0(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,IERR)
      DO NRL=1,NRMAX
         IF(NSW.EQ.0) THEN
            RSL=(DBLE(NRL)-0.5D0)*DR
         ELSEIF(NSW.EQ.1) THEN
            RSL= DBLE(NRL)       *DR
         ENDIF
         CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF",KFID,": IERR=",IERR
         FOUT(NRL)=F0*AMP
      ENDDO
C
 600  FORMAT(' ',A17,A10,A7,I2)
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
C     PV          : Surface(Peripheral) Function Value
C     PVS         : Surface(Peripheral) Function Value 
C                   Corresponding to RHOA
C     FOUT(NRMP)  : Function Values
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF2DSP(KFID,DR,PV,PVA,FOUT,AMP,RHOA,NRAMAX,NRMAX,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION RL(NRMU),TL(NTURM),F2(NRMU,NTURM),FOUT(NRMP)
      DIMENSION U(4,NRMU)
      CHARACTER KFID*10
C
      CALL UFREAD2_TIME(KFID,RL,TL,F2,NRLMAX,NTXMAX,MDCHK,IERR)
      CALL PRETREATMENT0(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,IERR)
      DO NRL=1,NRMAX
         RMN=(DBLE(NRL)-0.5D0)*DR
         CALL SPL1DF(RMN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF",KFID,": IERR=",IERR
         FOUT(NRL)=F0*AMP
      ENDDO
C
      RGN=DBLE(NRMAX)*DR
      CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
      IF(IERR.NE.0) WRITE(6,600) "XX TRFILE: SPL1DF",KFID,": IERR=",IERR
      PV=F0*AMP
      IF(RHOA.NE.1.D0) THEN
         RGN=DBLE(NRAMAX)*DR
         CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,600) "XX TRFILE: SPL1DF",KFID,": IERR=",IERR
         PVA=F0*AMP
      ENDIF
C
 600  FORMAT(' ',A17,A10,A7,I2)
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
C     MODE        : Boundary Condition in the center for Spline
C     ICK         : Check Indicator
C
C     output:
C
C     TL(NTURM)   : Total Time Data (The Number of DT)
C     FOUT(NRMP)  : Function Values
C     NTAMAX      : Maximum Time Step Number Corresponding to RHOA
C     IERR        : Error Indicator
C
C     *****************************************************
C
      SUBROUTINE UF2DT(KFID,DR,DT,TL,FOUT,AMP,NTAMAX,NRMAX,TLMAX,NSW,
     &                 MODE,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION RL(NRMU),TL(NTURM),F2(NRMU,NTURM),FOUT(NRMP,NTUM)
      DIMENSION U(4,4,NRMU,NTURM)
      CHARACTER KFID*10
C
      CALL UFREAD2_TIME(KFID,RL,TL,F2,NRLMAX,NTXMAX,MDCHK,IERR)
      CALL PRETREATMENT2(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,NTAMAX,TLMAX,
     &                   ICK,MODE,IERR)
      DO NRL=1,NRMAX
         IF(NSW.EQ.0) THEN
            RSL=(DBLE(NRL)-0.5D0)*DR
         ELSEIF(NSW.EQ.1) THEN
            RSL= DBLE(NRL)       *DR
         ENDIF
         DO NTA=1,NTAMAX
            TSL=DT*DBLE(NTA)
            CALL SPL2DF(RSL,TSL,F0,RL,TL,U,NRMU,NRLMAX,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,600) "XX TRFILE: SPL2DF",KFID,": IERR=",IERR
            FOUT(NRL,NTA)=F0*AMP
         ENDDO
      ENDDO
C
 600  FORMAT(' ',A17,A10,A7,I2)      
      RETURN
      END
C
      SUBROUTINE UF2DTP(KFID,DR,DT,PV,PVA,TL,FOUT,AMP,NTAMAX,RHOA,
     &                  NRAMAX,NRMAX,TLMAX,NSW,MODE,ICK,IERR)
C
      INCLUDE 'trcom0.inc'
      DIMENSION RL(NRMU),TL(NTURM),F2(NRMU,NTURM),FOUT(NRMP,NTUM)
      DIMENSION U(4,4,NRMU,NTURM)
      DIMENSION PV(NTUM),PVA(NTUM)
      CHARACTER KFID*10
C
      CALL UFREAD2_TIME(KFID,RL,TL,F2,NRLMAX,NTXMAX,MDCHK,IERR)
      CALL PRETREATMENT2(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,NTAMAX,TLMAX,
     &                   ICK,MODE,IERR)
      DO NRL=1,NRMAX
         IF(NSW.EQ.0) THEN
            RSL=(DBLE(NRL)-0.5D0)*DR
         ELSEIF(NSW.EQ.1) THEN
            RSL= DBLE(NRL)       *DR
         ENDIF
         DO NTA=1,NTAMAX
            TSL=DT*DBLE(NTA)
            CALL SPL2DF(RSL,TSL,F0,RL,TL,U,NRMU,NRLMAX,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,600) "XX TRFILE: SPL2DF",KFID,": IERR=",IERR
            FOUT(NRL,NTA)=F0*AMP
         ENDDO
      ENDDO
      DO NTA=1,NTAMAX
         DO NR=1,NRMAX
C            if(kfid.eq.'NE'.and.nta.le.3) write(6,*) NTA,NR,FOUT(NR,NTA)
         ENDDO
      ENDDO
C
      RGN=DBLE(NRMAX)*DR
      DO NTA=1,NTAMAX
         TSL=DT*DBLE(NTA)
         CALL SPL2DF(RGN,TSL,F0,RL,TL,U,NRMU,NRLMAX,NTXMAX,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) "XX TRFILE: SPL2DF",KFID,": IERR=",IERR
         PV(NTA)=F0*AMP
      ENDDO
      IF(RHOA.NE.1.D0) THEN
         RGN=DBLE(NRAMAX)*DR
         DO NTA=1,NTAMAX
            TSL=DT*DBLE(NTA)
            CALL SPL2DF(RGN,TSL,F0,RL,TL,U,NRMU,NRLMAX,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,600) "XX TRFILE: SPL2DF",KFID,": IERR=",IERR
            PVA(NTA)=F0*AMP
         ENDDO
      ENDIF
C
 600  FORMAT(' ',A17,A10,A7,I2)      
      RETURN
      END
C
C     ***********************************************************
C
C           PRETREATMENT SUBROUTINE FOR TR_TIME_UFILE
C
C     ***********************************************************
C
      SUBROUTINE PRETREATMENT0(KFID,RL,TL,F2,U,NRFMAX,NTXMAX,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION RL(NRMU),TL(NTURM),F2(NRMU,NTURM)
      DIMENSION DERIV(NRMU),U(4,NRMU)
      CHARACTER KFID*10
C
      IF(IERR.EQ.1) THEN
         DO NRF=1,NRFMAX
         DO NA=1,4
            U(NA,NRF)=0.D0
         ENDDO
         ENDDO
         IF(KFID.EQ.'CURTOT') MDCURT=1
         IERR=0
         RETURN
      ENDIF
C
      IF(KFID.EQ.'NM1') THEN
         MDNM1=1
C         IF(NSMAX.EQ.2) THEN
C            WRITE(6,*) "XX NSMAX=3 SHOULD BE RECOMMENDED."
C            PAUSE
C         ENDIF
      ENDIF
C
      DERIV(1)=0.D0
      DERIV(NRFMAX)=0.D0
      ID=0
C
      IF(NTXMAX.NE.1) THEN
         IF(TIME_INT.LE.0.D0) THEN
 100        WRITE(6,500) 'INPUT ARBITRARY TIME:',TL(1),' -',TL(NTXMAX)
            READ(5,*,ERR=100) TIME_INT
            IF(TIME_INT.LT.TL(1).OR.TIME_INT.GT.TL(NTXMAX)) GOTO 100
            DO NTX=1,NTXMAX
               IF(ABS(TL(NTX)-TIME_INT).LE.1.D-6) THEN
                  NISTP=NTX
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
                  IF(ABS(TL(NTX)-TIME_INT).LE.1.D-6) THEN
                     NISTP=NTX
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
 200     WRITE(6,510) 'TIME_INT=',TIME_INT,' HAS BEEN REPLACED BY',
     &              TL(NTX_MIN)
         TIME_INT=TL(NTX_MIN)
         NISTP=NTX_MIN
 1000    CONTINUE
C
         CALL SPL1D(RL,F2(1,NISTP),DERIV,U,NRFMAX,ID,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TRFILE: SPL1D',KFID,': IERR=',IERR
         
      ELSE
         CALL SPL1D(RL,F2(1,1),DERIV,U,NRFMAX,ID,IERR)
         IF(IERR.NE.0)
     &        WRITE(6,*) 'XX TRFILE: SPL1D',KFID,': IERR=',IERR
      ENDIF
C
      RETURN
 500  FORMAT(' ',A,F8.4,A,F9.5)
 510  FORMAT(' ',A,F8.4,A,F8.4)
      END
C
C     *****
C
      SUBROUTINE PRETREATMENT1(KFID,TL,F1,U,NTXMAX,NTLMAX,TLMAX,
     &                         ICK,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION TL(NTURM),DERIV(NTURM),U(4,NTURM)
      CHARACTER KFID*10
C
      IF(IERR.EQ.1) THEN
         IERR=0
         DO NTX=1,NTXMAX
         DO NA=1,4
            U(NA,NTX)=0.D0
         ENDDO
         ENDDO
         RETURN
      ENDIF
C
      DERIV(1)=0.D0
      DERIV(NTXMAX)=0.D0
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
      CALL SPL1D(TL,F1,DERIV,U,NTXMAX,3,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TRFILE: SPL1D',KFID,': IERR=',IERR
C
      RETURN
      END
C
C     *****
C
      SUBROUTINE PRETREATMENT2(KFID,RL,TL,F2,U,NRFMAX,NTXMAX,NTLMAX,
     &                         TLMAX,ICK,MODE,IERR)
C
      INCLUDE 'trcomm.inc'
C      COMMON /PRETREAT1/ RUF(NRMU),TMU(NTURM),F1(NTURM),F2(NRMU,NTURM)
C      COMMON /PRETREAT2/ NTAMAX
      DIMENSION RL(NRMU),TL(NTURM),F2(NRMU,NTURM)
      DIMENSION DERIVX(NRMU,NTURM),DERIVY(NRMU,NTURM)
      DIMENSION DERIVXY(NRMU,NTURM)
      DIMENSION U(4,4,NRMU,NTURM)
      CHARACTER KFID*10
      DIMENSION F3(NRMU),DERIV(NRMU),UT(4,NRMU)
C
      IF(IERR.EQ.1) THEN
         IERR=0
         DO NTX=1,NTXMAX
         DO NRF=1,NRFMAX
         DO NB=1,4
         DO NA=1,4
            U(NA,NB,NRF,NTX)=0.D0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         RETURN
      ENDIF
C
      IF(KFID.EQ.'NM1') THEN
         MDNM1=1
C         IF(NSMAX.EQ.2) THEN
C            WRITE(6,*) "XX NSMAX=3 SHOULD BE RECOMMENDED."
C            PAUSE
C         ENDIF
      ENDIF
C
      DO NTX=1,NTXMAX
         DERIVX(1,NTX)=0.D0
      ENDDO
      DERIVXY(1,     1)=0.D0
      DERIVXY(1,NTXMAX)=0.D0
C
      DO NTX=2,NTXMAX
         TL(NTX)=(TL(NTX)-TL(1))
      ENDDO
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX TRFILE:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
            STOP
         ENDIF
      ENDIF
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
C
      IF(MODE.EQ.0) THEN
         CALL SPL2D(RL,TL,F2,DERIVX,DERIVY,DERIVXY,U,
     &              NRMU,NRFMAX,NTXMAX,1,0,IERR)
      ELSEIF(MODE.EQ.1) THEN
         CALL SPL2D(RL,TL,F2,DERIVX,DERIVY,DERIVXY,U,
     &              NRMU,NRFMAX,NTXMAX,0,0,IERR)
      ELSE
         DERIV(1)=0.D0
         DERIV(NTXMAX)=0.D0
C
         NTN=MODE-1
         IF(NTN.LE.0) STOP
         DO NRF=1,NRFMAX
            F3(NRF)=F2(NRF,NTN)
         ENDDO
         CALL SPL1D(RL,F3,DERIV,UT,NRFMAX,3,IERR)
         DO NR=1,NRMAX
            RMN=(DBLE(NR)-0.5D0)*DR
            CALL SPL1DF(RMN,F0,RL,UT,NRFMAX,IERR)
            write(6,*) RMN,F0
         ENDDO
      ENDIF
      IF(IERR.NE.0) WRITE(6,*) 'XX TRFILE: SPL2D',KFID,': IERR=',IERR
      ICK=1
c$$$C
c$$$      TMLCL=DT*DBLE(4)
c$$$      DO NR=1,NRMAX
c$$$         RMN=(DBLE(NR)-0.5D0)*DR
c$$$         CALL SPL2DF(RMN,TMLCL,F0,RL,TL,U,
c$$$     &                  NRMU,NRFMAX,NTXMAX,IERR)
c$$$         write(6,*) NR,F0*1.D-3
c$$$      ENDDO
C
      RETURN
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
C
      RKAP=RKAPU(NT)
      RKAPS=SQRT(RKAPU(NT))
C
      RR=RRU(NT)
      RA=RAU(NT)
      BB=BBU(NT)
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IF(MDLUF.EQ.1) THEN
      DO NR=1,NRMAX
         RN(NR,1)=RNU(NR,1,NT)
         RN(NR,2)=RNU(NR,2,NT)
         IF(MDNI.NE.0) RN(NR,3)=RNU(NR,3,NT)
         IF(INS.NE.0) THEN
            RT(NR,2)=RTU(NR,2,NT)
            RT(NR,3)=RTU(NR,3,NT)
         ENDIF
         QP(NR)=QPU(NR,NT)
         PEX(NR,1)=PNBU(NR,1,NT)
         PEX(NR,2)=PNBU(NR,2,NT)
         PRF(NR,1)=PICU(NR,1,NT)+PECU(NR,NT)
         PRF(NR,2)=PICU(NR,2,NT)
         TTRHO(NR)=TTRHOU(NR,NT)
         DVRHO(NR)=DVRHOU(NR,NT)
         DSRHO(NR)=DSRHOU(NR,NT)
         ABRHO(NR)=ABRHOU(NR,NT)
         ARRHO(NR)=ARRHOU(NR,NT)
         AR1RHO(NR)=AR1RHOU(NR,NT)
         AR2RHO(NR)=AR2RHOU(NR,NT)
C         RJCB(NR)=1.D0/(RKAPS*RA)
         RJCB(NR)=AR1RHOU(NR,NT)
         RMJRHO(NR)=RMJRHOU(NR,NT)
         RMNRHO(NR)=RMNRHOU(NR,NT)
         EKAPPA(NR)=RKAP
C         BP(NR)=BPU(NR,NT)
C         ZEFF(NR)=ZEFFU(NR,NT)
      ENDDO
      CALL TRGFRG
      IF(MDLJQ.EQ.0) THEN
         NR=1
            AJ(NR)=AJU(NR,NT)
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         DO NR=2,NRMAX
            AJ(NR)=AJU(NR,NT)
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         ENDDO
         NR=1
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(AMYU0*DVRHO(NR))
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
C         RN(NR,1)=RNU(NR,1,NT)
C         RN(NR,2)=RNU(NR,2,NT)
C         QP(NR)=QPU(NR,NT)
         PEX(NR,1)=PNBU(NR,1,NT)
         PEX(NR,2)=PNBU(NR,2,NT)
         PRF(NR,1)=PICU(NR,1,NT)+PECU(NR,NT)
         PRF(NR,2)=PICU(NR,2,NT)
         TTRHO(NR)=TTRHOU(NR,NT)
         DVRHO(NR)=DVRHOU(NR,NT)
         DSRHO(NR)=DSRHOU(NR,NT)
         ABRHO(NR)=ABRHOU(NR,NT)
         ARRHO(NR)=ARRHOU(NR,NT)
         AR1RHO(NR)=AR1RHOU(NR,NT)
         AR2RHO(NR)=AR2RHOU(NR,NT)
C         RJCB(NR)=1.D0/(RKAPS*RA)
         RJCB(NR)=AR1RHOU(NR,NT)
         RMJRHO(NR)=RMJRHOU(NR,NT)
         RMNRHO(NR)=RMNRHOU(NR,NT)
         EKAPPA(NR)=RKAP
      ENDDO
      ENDIF
C      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
      CALL TRGFRG
C
      IF(MDLUF.EQ.1) THEN
         DO NS=1,2
            PNS (NS)=PNSU (NS,NT)
            PTS (NS)=PTSU (NS,NT)
            PNSA(NS)=PNSUA(NS,NT)
            PTSA(NS)=PTSUA(NS,NT)
         ENDDO
         IF(MDNI.NE.0) THEN 
            PNS (3)=PNSU (3,NT)
            PTS (3)=PTSU (3,NT)
            PNSA(3)=PNSUA(3,NT)
            PTSA(3)=PTSUA(3,NT)
         ENDIF
         IF(RHOA.NE.1.D0) THEN
            DO NR=NRAMAX+1,NRMAX
               PROF    =(1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,1)= RTU(NR,1,NT)
               RT(NR,2)= RTU(NR,2,NT)
               RT(NR,3)= RTU(NR,2,NT)
               RT(NR,4)=(RTU(NR,2,NT)-RTU(NRMAX,2,NT))*PROF
     &                  +RTU(NRMAX,2,1)
            ENDDO
            IF(MDNI.EQ.0) THEN
               DO NR=NRAMAX+1,NRMAX
                  PROF    =(1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
                  RT(NR,3)=(RTU(NR,2,NT)-RTU(NRMAX,2,NT))*PROF
     &                     +RTU(NRMAX,2,1)
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
      RKAP=RKAPU(1)
      RKAPS=SQRT(RKAP)
C
      RR=RRU(1)
      RA=RAU(1)
      BB=BBU(1)
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      DO NR=1,NRMAX
         RN(NR,1)=RNU(NR,1,1)
         RN(NR,2)=RNU(NR,2,1)
         IF(MDNI.NE.0) RN(NR,3)=RNU(NR,3,1)
         IF(INS.NE.0) THEN
            RT(NR,2)=RTU(NR,2,1)
            RT(NR,3)=RTU(NR,3,1)
         ENDIF
         QP(NR)=QPU(NR,1)
         PEX(NR,1)=PNBU(NR,1,1)
         PEX(NR,2)=PNBU(NR,2,1)
         PRF(NR,1)=PICU(NR,1,1)+PECU(NR,1)
         PRF(NR,2)=PICU(NR,2,1)
         TTRHO(NR)=TTRHOU(NR,1)
         DVRHO(NR)=DVRHOU(NR,1)
         DSRHO(NR)=DSRHOU(NR,1)
         ABRHO(NR)=ABRHOU(NR,1)
         ARRHO(NR)=ARRHOU(NR,1)
         AR1RHO(NR)=AR1RHOU(NR,1)
         AR2RHO(NR)=AR2RHOU(NR,1)
C         RJCB(NR)=1.D0/(RKAPS*RA)
         RJCB(NR)=AR1RHOU(NR,1)
         RMJRHO(NR)=RMJRHOU(NR,1)
         RMNRHO(NR)=RMNRHOU(NR,1)
         EKAPPA(NR)=RKAP
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
C   **************************************************
C   **    UFILE read for TR (Time Evolution UFILE)  **
C   **************************************************
C
C     input:
C
C     KFID     : Variable name in UFILEs
C
C     output:
C
C     RUF(NRMU)      : Equally Spaced Normalized Radial Data
C     TMU(NTURM)     : Total Time Data (The Number of DT)
C     F1(NTURM)      : Functional Values
C     F2(NRMU,NTURM) : Functional Values
C     NRFMAX         : Maximum Number of the Radial Mesh
C     NTXMAX         : Maximum Number of the Time Mesh
C     MDCHK          : Loop Check Value
C     IERR           : Error Indicator
C
C   ***************************************************************
C
      SUBROUTINE UFREAD_TIME(KFID,TT,F1,NTXMAX,MDCHK,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION TT(NTURM),F1(NTURM)
      CHARACTER KDIRX*80
      CHARACTER KDIRR1*80
      CHARACTER KFILE*80,KFID*10
      COMMON /TRUFC1/ KDIRX
      COMMON /TRUFC2/ IKNDEV,IKNDCG
      LOGICAL LEX
C
C      IF(MDCHK.NE.0) GOTO 9000
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR1=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'1d'//KUFDCG(1:IKNDCG)//'.'
C
      CALL KTRIM(KDIRR1,KL1)
      KFILE=KDIRR1(1:KL1)//KFID
C
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000) 
      IF(LEX) THEN
         CALL TRXR1D(KDIRR1,KFID,TT,F1,NTURM,NTXMAX,0)
         MDCHK=1
         IERR=0
      ELSE
         DO NTA=1,NTURM
            F1(NTA)=0.D0
         ENDDO
         IERR=1
      ENDIF
C
 9000 RETURN
      END
C
C     *****
C
      SUBROUTINE UFREAD2_TIME(KFID,RUF,TMU,F2,NRFMAX,NTXMAX,MDCHK,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION RUF(NRMU),TMU(NTURM),F2(NRMU,NTURM)
      DIMENSION F2CTR(NTURM),F2EDG(NTURM)
      CHARACTER KDIRX*80
      CHARACTER KDIRR2*80
      CHARACTER KFILE*80,KFID*10
      COMMON /TRUFC1/ KDIRX
      COMMON /TRUFC2/ IKNDEV,IKNDCG
      LOGICAL LEX
C
C      IF(MDCHK.NE.0) GOTO 9000
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'2d'//KUFDCG(1:IKNDCG)//'.'
C
      CALL KTRIM(KDIRR2,KL2)
      KFILE=KDIRR2(1:KL2)//KFID
C
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         NRFMAX=52              ! equal to NRMU
         CALL TRXR2D(KDIRR2,KFID,TMU,RUF,F2,NRMU,NTURM,NRFMAX,NTXMAX,0)
         MDCHK=1
         IERR=0
      ELSE
         DO NTA=1,NTURM
            DO NRF=1,NRMU
               F2(NRF,NTA)=0.D0
            ENDDO
         ENDDO
         IERR=1
         GOTO 9000
      ENDIF
C     
      MD=0
      IF(RUF(1).NE.0.D0.AND.RUF(NRFMAX).NE.1.D0) THEN
         DO NTX=1,NTXMAX
            F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(1,NTX),F2(2,NTX))
            F2EDG(NTX)=FEDG(1.D0,RUF(NRFMAX-1),RUF(NRFMAX),
     &                      F2(NRFMAX-1,NTX),F2(NRFMAX,NTX))
         ENDDO
         MD=1
      ELSEIF(RUF(1).NE.0.D0.AND.RUF(NRFMAX).EQ.1.D0) THEN
         DO NTX=1,NTXMAX
            F2CTR(NTX)=FCTR(RUF(1),RUF(2),F2(1,NTX),F2(2,NTX))
         ENDDO
         MD=2
      ELSEIF(RUF(1).EQ.0.D0.AND.RUF(NRFMAX).NE.1.D0) THEN
         DO NTX=1,NTXMAX
            F2EDG(NTX)=FEDG(1.D0,RUF(NRFMAX-1),RUF(NRFMAX),
     &                      F2(NRFMAX-1,NTX),F2(NRFMAX,NTX))
         ENDDO
         MD=3
      ELSEIF(RUF(1).EQ.0.D0.AND.RUF(NRFMAX).EQ.1.D0) THEN
         MD=4
      ELSE
         STOP 'XX TRFILE: TRXR2D: ERROR'
      ENDIF
C     
      CALL DATA_ERROR_CORRECT(KFID,RUF,F2,NTXMAX)
C
      IF(MD.EQ.1) THEN
         NRFMAX=NRFMAX+2
         DO NRF=NRFMAX-1,2,-1
            RUF(NRF)=RUF(NRF-1)
         ENDDO
         RUF(1)=0.D0
         RUF(NRFMAX)=1.D0
         DO NTX=1,NTXMAX
            DO NRF=NRFMAX-1,2,-1
               F2(NRF,NTX)=F2(NRF-1,NTX)
            ENDDO
            F2(1,NTX)=F2CTR(NTX)
            F2(NRFMAX,NTX)=F2EDG(NTX)
         ENDDO
      ELSEIF(MD.EQ.2) THEN
         NRFMAX=NRFMAX+1
         DO NRF=NRFMAX,2,-1
            RUF(NRF)=RUF(NRF-1)
         ENDDO
         RUF(1)=0.D0
         DO NTX=1,NTXMAX
            DO NRF=NRFMAX,2,-1
               F2(NRF,NTX)=F2(NRF-1,NTX)
            ENDDO
            F2(1,NTX)=F2CTR(NTX)
         ENDDO
      ELSEIF(MD.EQ.3) THEN
         NRFMAX=NRFMAX+1
         DO NTX=1,NTXMAX
            F2(NRFMAX,NTX)=F2EDG(NTX)
         ENDDO
      ENDIF
C
 9000 RETURN
      END
C
C     *****
C
      FUNCTION FCTR(R1,R2,F1,F2)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      FCTR = (R2**2*F1-R1**2*F2)/(R2**2-R1**2)
C
      RETURN
      END
C
      FUNCTION FEDG(R0,R1,R2,F1,F2)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      FEDG = ((R2-R0)*F1-(R1-R0)*F2)/(R2-R1)
C
      RETURN
      END
C
C
C     *****
C
      SUBROUTINE UFREAD2_ERROR(KFID,RUF,TMU,F2,NRFMAX,NTXMAX,MDCHK,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION RUF(NRMU),TMU(NTURM),F2(NRMU,NTURM)
      CHARACTER KDIRX*80
      CHARACTER KDIRR2*80
      CHARACTER KFILE*80,KFID*10
      COMMON /TRUFC1/ KDIRX
      COMMON /TRUFC2/ IKNDEV,IKNDCG
      LOGICAL LEX
C
C      IF(MDCHK.NE.0) GOTO 9000
C
      CALL KTRIM(KDIRX,IKDIRX)
      KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)
     &       //'2d'//KUFDCG(1:IKNDCG)//'.'
C
      CALL KTRIM(KDIRR2,KL2)
      KFILE=KDIRR2(1:KL2)//KFID
C
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         NRFMAX=NRMU           ! equal to NRMU
         CALL TRXR2D(KDIRR2,KFID,TMU,RUF,F2,NRMU,NTURM,NRFMAX,NTXMAX,0)
         MDCHK=1
         IERR=0
      ELSE
         DO NTA=1,NTURM
            DO NRF=1,NRMU
               F2(NRF,NTA)=0.D0
            ENDDO
         ENDDO
         IERR=1
         GOTO 9000
      ENDIF
C
      IF(RUF(NRFMAX).GT.1.D0) THEN
         DO NRF=1,NRFMAX
            IF(RUF(NRF).GT.1.D0) THEN
               NRFMAX=NRF-1
               GOTO 1000
            ENDIF
         ENDDO
      ENDIF
 1000 CONTINUE
      IF(KFID.EQ.'NEXP'.OR.KFID.EQ.'NEXPEB') THEN
         DO NTX=1,NTXMAX
            DO NRF=1,NRMU
               F2(NRF,NTX)=F2(NRF,NTX)*1.D-20
            ENDDO
         ENDDO
      ELSE
         DO NTX=1,NTXMAX
            DO NRF=1,NRMU
               F2(NRF,NTX)=F2(NRF,NTX)*1.D-3
            ENDDO
         ENDDO
      ENDIF
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           WRONG DATA CORRECT
C
C     ***********************************************************
C
      SUBROUTINE DATA_ERROR_CORRECT(KFID,RUF,F2,NTXMAX)
C
      INCLUDE 'trcomm.inc'
      DIMENSION RUF(NRMU),F2(NRMU,NTURM)
      CHARACTER KFID*10
C
      IF(KUFDEV.EQ.'jet') THEN
         IF(KFID.EQ.'GRHO1') THEN
            IF(KUFDCG.EQ.'57987'.OR.KUFDCG.EQ.'58159'.OR.
     &         KUFDCG.EQ.'58323') THEN
               DO NTX=1,NTXMAX
                  F2(1,NTX)=FCTR(RUF(2),RUF(3),F2(2,NTX),F2(3,NTX))
               ENDDO
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
