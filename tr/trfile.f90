! trfile.f90

MODULE trfile

  PRIVATE
  PUBLIC tr_save,tr_load,tr_grsv,tr_grld

CONTAINS

!     ***********************************************************

!           SAVE TRANSPORT DATA

!     ***********************************************************

      SUBROUTINE tr_save

      USE TRCOMM, ONLY : ABB2RHOG, ABRHO, ABRHOG, AD0, AIB2RHOG, AJBST, AJNBT, AJOH, AJOHT, AJRFT, AJT, ALI, ALP, AMZ, ANC,       &
     &                   ANFE, ANNU, ANSAV, AR1RHO, AR1RHOG, AR2RHO, AR2RHOG, ARHBRHOG, ARRHO, ARRHOG, AV0, BB, BETA0, BETAA,     &
     &                   BETAP0, BETAPA, BP, CALF, CDH, CDP, CDW, CHP, CK0, CK1, CKALFA, CKBETA, CKGUMA, CNB, CNH, CNP, CWEB,     &
     &                   DR, DT, DVRHO, DVRHOG, EPSLTR, EPSRHO, ETA, EZOH, IZERO, KFNLOG, KNAMTR, KUFDCG, KUFDEV, LMAXTR, MDCD05, &
     &                   MDDIAG, MDDW, MDLAD, MDLAVK, MDLCD, MDLEC, MDLEOI, MDLEQ0, MDLEQB, MDLEQE, MDLEQN, MDLEQT, MDLEQU,       &
     &                   MDLEQZ, MDLER, MDLETA, MDLFLX, MDLIC, MDLJBS, MDLJQ, MDLKAI, MDLKNC, MDLLH, MDLNB, MDLNF, MDLPCK,        &
     &                   MDLPEL, MDLST, MDLTPF, MDLUF, MDLWLD, MDLXP, MDNCLS, MDNI, MDTC, MODELG, MODEP, NEA, NEQMAX, NGPST,      &
     &                   NGRSTP, NGTSTP, NNS, NRAMAX, NRMAX, NROMAX, NSMAX, NSNMAX, NSS, NST, NSV, NSZMAX, NTEQIT, NTMAX,         &
     &                   NTMAX_SAVE, NTSTEP, PA, PBSCD, PCXT, PECCD, PECNPR, PECR0, PECRW, PECTOE, PECTOT, PELPAT, PELR0, PELRAD, &
     &                   PELRW, PELTIM, PELTOT, PELVEL, PHIA, PICCD, PICNPR, PICR0, PICRW, PICTOE, PICTOT, PIET, PINT, PLHCD,     &
     &                   PLHNPR, PLHR0, PLHRW, PLHTOE, PLHTOT, PN, PNBCD, PNBENG, PNBR0, PNBRTG, PNBRW, PNBT, PNBTOT, PNC, PNFE,  &
     &                   PNNU, PNNUS, PNS, PNSS, POHT, POUT, PRFT, PRLT, PROFJ1, PROFJ2, PROFN1, PROFN2, PROFT1, PROFT2, PROFU1,  &
     &                   PROFU2, PT, PTS, PZ, Q0, RA, RDLT, RDP, RG, RHOA, RHOG, RHOM, RIPE, RIPS, RJCB, RKAP, RKPRHO, RKPRHOG,   &
     &                   RM, RMJRHO, RMNRHO, RN, RNF, RPSI, RR, RT, RTF, RU, RW, T, TAUE1, TAUE2, TAUE89, TIME_INT, TPRE, TPRST,  &
     &                   TS0, TSAV, TSST, TST, TTRHO, TTRHOG, VLOOP, VPOL, VSEC, VTOR, WPPRE, WPT, WST, &
     &                   ABVRHO, ABVRHOG, RDPVRHOG, AJTTOR
      USE libfio
      IMPLICIT NONE
      INTEGER(4):: IERR, NDD, NDM, NDY, NTH1, NTM1, NTS1
      CHARACTER(LEN=3):: K1, K2, K3, K4, K5, K6


      CALL FWOPEN(21,KNAMTR,0,4,'TR',IERR)
      IF(IERR.NE.0) RETURN

      WRITE(21) MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI
      WRITE(21) NRMAX,NRAMAX,NROMAX,DT,NGPST,TSST
      WRITE(21) NTMAX,NTSTEP,NGTSTP,NGRSTP,NSMAX,NSZMAX,NSNMAX
      WRITE(21) NSS,NSV,NNS,NST,NEA,NEQMAX
      WRITE(21) RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA,PHIA
      WRITE(21) PA,PZ,PN,PNS,PT,PTS
      WRITE(21) PNC,PNFE,PNNU,PNNUS
      WRITE(21) PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,PROFJ1,PROFJ2, &
     &          ALP,AD0,AV0,CNP,CNH,CDP,CDH,CDW,CWEB,CALF,CNB
      WRITE(21) MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF,MDCD05
      WRITE(21) TPRST,MDLST,MDLNF,IZERO
      WRITE(21) MODELG,NTEQIT
      WRITE(21) MDLXP,MDLUF,MDNCLS,MDLWLD,MDDIAG,MDDW,MDLFLX,MDLER
      WRITE(21) KUFDEV,KUFDCG,TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK
      WRITE(21) PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB
      WRITE(21) PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC
      WRITE(21) PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH
      WRITE(21) PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC
      WRITE(21) PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD
      WRITE(21) PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,PELPAT,MDLPEL
      WRITE(21) DR,PNSS,T,TST,VSEC,WPPRE,TPRE,KFNLOG,NTMAX_SAVE
      WRITE(21) RG,RM,RN,RT,RU,RW,BP,RDP,RPSI,RNF,RTF,ANC,ANFE,ANNU,RDPVRHOG
      WRITE(21) VTOR,VPOL,AJOH,EZOH,ETA,AMZ
      WRITE(21) EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA,CNB,CALF
      WRITE(21) DVRHO,TTRHO,ABRHO,ABVRHO,ARRHO,AR1RHO,AR2RHO,RMJRHO,RMNRHO,RKPRHO,RJCB,EPSRHO
      WRITE(21) DVRHOG,TTRHOG,ABRHOG,ABVRHOG,ARRHOG,AR1RHOG,AR2RHOG,ABB2RHOG,AIB2RHOG,ARHBRHOG,RKPRHOG
      WRITE(21) RHOM,RHOG
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'


!      OPEN(16,POSITION='APPEND',FILE=KFNLOG)
!      OPEN(16,ACCESS='APPEND',FILE=KFNLOG)
      OPEN(16,ACCESS='SEQUENTIAL',FILE=KFNLOG)

         CALL GUDATE(NDY,NDM,NDD,NTH1,NTM1,NTS1)
         WRITE(K1,'(I3)') 100+NDY
         WRITE(K2,'(I3)') 100+NDM
         WRITE(K3,'(I3)') 100+NDD
         WRITE(K4,'(I3)') 100+NTH1
         WRITE(K5,'(I3)') 100+NTM1
         WRITE(K6,'(I3)') 100+NTS1
         WRITE(16,1670) K1(2:3),K2(2:3),K3(2:3),K4(2:3),K5(2:3),K6(2:3), &
     &                  RIPS,RIPE,PN(1),PN(2),BB,PICTOT,PLHTOT,PLHNPR
 1670    FORMAT(' '/ &
     &          ' ','## DATE: ', &
     &              A2,'-',A2,'-',A2,'  ',A2,':',A2,':',A2,' : ', &
     &              '  FILE: ',A40/ &
     &          ' ',3X,'RIPS  =',1PD10.3,'  RIPE  =',1PD10.3, &
     &               '  PNE   =',1PD10.3,'  PNI   =',1PD10.3/ &
     &          ' ',3X,'BB    =',1PD10.3,'  PICTOT=',1PD10.3, &
     &               '  PLHTOT=',1PD10.3,'  PLHNPR=',1PD10.3)
         WRITE(16,1671) T, &
     &                WPT,TAUE1,TAUE2,TAUE89, &
     &                BETAP0,BETAPA,BETA0,BETAA
 1671    FORMAT(' ','# TIME : ',F7.3,' SEC'/ &
     &          ' ',3X,'WPT   =',1PD10.3,'  TAUE  =',1PD10.3, &
     &               '  TAUED =',1PD10.3,'  TAUE89=',1PD10.3/ &
     &          ' ',3X,'BETAP0=',1PD10.3,'  BETAPA=',1PD10.3, &
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)

         WRITE(16,1672) WST(1),TS0(1),TSAV(1),ANSAV(1), &
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
 1672    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3, &
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/ &
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3, &
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)

!!$         WRITE(16,1673) AJT,VLOOP,ALI,Q0, &
!!$     &                AJOHT,AJNBT,AJRFT,AJBST
!!$ 1673    FORMAT(' ',3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3, &
!!$     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/ &
!!$     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3, &
!!$     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)
         WRITE(16,1673) AJTTOR,VLOOP,ALI,Q0, &
     &                AJT,AJOHT,AJNBT,AJBST
 1673    FORMAT(' ',3X,'AJTTOR=',1PD10.3,'  VLOOP =',1PD10.3, &
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/ &
     &          ' ',3X,'AJT   =',1PD10.3,'  AJOHT =',1PD10.3, &
     &               '  AJNBT =',1PD10.3,'  AJBST =',1PD10.3)

         WRITE(16,1674) PINT,POHT,PNBT, &
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4), &
     &                POUT,PRLT,PCXT,PIET
 1674    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3, &
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/ &
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3, &
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3)

      CLOSE(16)

  900 RETURN
      END SUBROUTINE tr_save

!     ***********************************************************

!           LOAD TRANSPORT DATA

!     ***********************************************************

      SUBROUTINE tr_load

      USE TRCOMM, ONLY : ABB2RHOG, ABRHO, ABRHOG, AD0, AIB2RHOG, AJOH, ALP, AMZ, ANC, ANFE, ANNU, AR1RHO, AR1RHOG, AR2RHO,      &
     &                   AR2RHOG, ARHBRHOG, ARRHO, ARRHOG, AV0, BB, BP, CALF, CDH, CDP, CDW, CHP, CK0, CK1, CKALFA, CKBETA,     &
     &                   CKGUMA, CNB, CNH, CNP, CWEB, DR, DT, DVRHO, DVRHOG, EPSLTR, EPSRHO, ETA, EZOH, GRG, GRM, IREAD, IZERO, &
     &                   KFNLOG, KNAMTR, KUFDCG, KUFDEV, LMAXTR, MDCD05, MDDIAG, MDDW, MDLAD, MDLAVK, MDLCD, MDLEC, MDLEOI,     &
     &                   MDLEQ0, MDLEQB, MDLEQE, MDLEQN, MDLEQT, MDLEQU, MDLEQZ, MDLER, MDLETA, MDLFLX, MDLIC, MDLJBS, MDLJQ,   &
     &                   MDLKAI, MDLKNC, MDLLH, MDLNB, MDLNF, MDLPCK, MDLPEL, MDLST, MDLTPF, MDLUF, MDLWLD, MDLXP, MDNCLS,      &
     &                   MDNI, MDTC, MODELG, MODEP, NEA, NEQMAX, NGPST, NGR, NGRSTP, NGST, NGT, NGTSTP, NNS, NRAMAX, NRMAX,     &
     &                   NROMAX, NSMAX, NSNMAX, NSS, NST, NSV, NSZMAX, NTEQIT, NTMAX, NTMAX_SAVE, NTSTEP, PA, PBSCD, PECCD,     &
     &                   PECNPR, PECR0, PECRW, PECTOE, PECTOT, PELPAT, PELR0, PELRAD, PELRW, PELTIM, PELTOT, PELVEL, PHIA, PI,  &
     &                   PICCD, PICNPR, PICR0, PICRW, PICTOE, PICTOT, PLHCD, PLHNPR, PLHR0, PLHRW, PLHTOE, PLHTOT, PN, PNBCD,   &
     &                   PNBENG, PNBR0, PNBRTG, PNBRW, PNBTOT, PNC, PNFE, PNNU, PNNUS, PNS, PNSS, PROFJ1, PROFJ2, PROFN1,       &
     &                   PROFN2, PROFT1, PROFT2, PROFU1, PROFU2, PT, PTS, PZ, Q0, QP, RA, RDLT, RDP, RG, RHOA, RHOG, RHOM, RIPE,&
     &                   RIPS, RJCB, RKAP, RKPRHO, RKPRHOG, RM, RMJRHO, RMNRHO, RN, RNF, RPSI, RR, RT, RTF, RU, RW, T, TIME_INT,&
     &                   TPRE, TPRST, TSST, TST, TTRHO, TTRHOG, VPOL, VSEC, VTOR, WPPRE, &
     &                   ALLOCATE_TRCOMM, ABVRHO, ABVRHOG, RDPVRHOG
      USE libfio
      IMPLICIT NONE
      INTEGER(4):: IERR
      REAL(8)   :: FCTR


      CALL FROPEN(21,KNAMTR,0,1,'TR',IERR)
      IF(IERR.NE.0) RETURN


      READ(21) MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI
      READ(21) NRMAX,NRAMAX,NROMAX,DT,NGPST,TSST
      READ(21) NTMAX,NTSTEP,NGTSTP,NGRSTP,NSMAX,NSZMAX,NSNMAX
        CALL ALLOCATE_TRCOMM(IERR)
        IF(IERR.NE.0) RETURN
      READ(21) NSS,NSV,NNS,NST,NEA,NEQMAX
      READ(21) RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA,PHIA
      READ(21) PA,PZ,PN,PNS,PT,PTS
      READ(21) PNC,PNFE,PNNU,PNNUS
      READ(21) PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,PROFJ1,PROFJ2, &
     &         ALP,AD0,AV0,CNP,CNH,CDP,CDH,CDW,CWEB,CALF,CNB
      READ(21) MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF,MDCD05
      READ(21) TPRST,MDLST,MDLNF,IZERO
      READ(21) MODELG,NTEQIT
      READ(21) MDLXP,MDLUF,MDNCLS,MDLWLD,MDDIAG,MDDW,MDLFLX,MDLER
      READ(21) KUFDEV,KUFDCG,TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK
      READ(21) PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB
      READ(21) PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC
      READ(21) PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH
      READ(21) PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC
      READ(21) PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD
      READ(21) PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,PELPAT,MDLPEL
      READ(21) DR,PNSS,T,TST,VSEC,WPPRE,TPRE,KFNLOG,NTMAX_SAVE
      READ(21) RG,RM,RN,RT,RU,RW,BP,RDP,RPSI,RNF,RTF,ANC,ANFE,ANNU,RDPVRHOG
      READ(21) VTOR,VPOL,AJOH,EZOH,ETA,AMZ
      READ(21) EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA,CNB,CALF
      READ(21) DVRHO,TTRHO,ABRHO,ABVRHO,ARRHO,AR1RHO,AR2RHO,RMJRHO,RMNRHO,RKPRHO,RJCB,EPSRHO
      READ(21) DVRHOG,TTRHOG,ABRHOG,ABVRHOG,ARRHOG,AR1RHOG,AR2RHOG,ABB2RHOG,AIB2RHOG,ARHBRHOG,RKPRHOG
      READ(21) RHOM,RHOG
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      IREAD=1

      NGR=0
      NGT=0
      NGST=0
      RIPS=RIPE
      GRG(1)=0.0
      GRM(1:NRMAX)  =SNGL(RM(1:NRMAX))
      GRG(2:NRMAX+1)=SNGL(RG(1:NRMAX))
      QP(1:NRMAX)   =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
!     *** calculate q_axis ***
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))

  900 RETURN
      END SUBROUTINE tr_load

!     ***********************************************************

!           SAVE GRAPHIC DATA

!     ***********************************************************

      SUBROUTINE tr_grsv

      USE TRCOMM, ONLY : GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT
      IMPLICIT NONE
      INTEGER(4):: IST
      CHARACTER(LEN=32):: TRFLNM
      LOGICAL   :: LEX


    1 WRITE(6,*) '# INPUT : GRAPHIC SAVE FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFLNM
      IF(TRFLNM.EQ.'                                ') GOTO 900
      INQUIRE (FILE=TRFLNM,EXIST=LEX)
      IF(LEX) THEN
!         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
!     &              'ARE YOU SURE {Y/N}?'
!         READ(5,'(A1)',ERR=1,END=900) KID
!         CALL GUCPTL(KID)
!         IF(KID.NE.'Y') GOTO 1
         OPEN(22,FILE=TRFLNM,IOSTAT=IST,STATUS='OLD',ERR=10,FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',TRFLNM,') IS ASSIGNED FOR OUTPUT.'
         GOTO 30
   10    WRITE(6,*) '# XX OLD FILE OPEN ERROR : IOSTAT=',IST
         GOTO 1
      ELSE
         OPEN(22,FILE=TRFLNM,IOSTAT=IST,STATUS='NEW',ERR=20,FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (',TRFLNM,') IS CREATED FOR OUTPUT.'
         GOTO 30
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT=',IST
         GOTO 1
      ENDIF

   30 WRITE(22) GVR,GRM,GRG,GTR,NGR
      WRITE(22) GVT,GT,NGT
      CLOSE(22)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

  900 RETURN
      END SUBROUTINE tr_grsv

!     ***********************************************************

!           LOAD GRAPHIC DATA

!     ***********************************************************

      SUBROUTINE tr_grld

      USE TRCOMM, ONLY : GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT
      IMPLICIT NONE
      INTEGER(4):: IST
      CHARACTER(LEN=32):: TRFLNM
      LOGICAL   :: LEX


    1 WRITE(6,*) '# INPUT : GRAPHIC LOAD FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFLNM
      IF(TRFLNM.EQ.'                                ') GOTO 900
      INQUIRE (FILE=TRFLNM,EXIST=LEX)
      IF(LEX) THEN
         OPEN(22,FILE=TRFLNM,IOSTAT=IST,STATUS='OLD',ERR=10,FORM='UNFORMATTED')
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

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  900 RETURN
      END SUBROUTINE tr_grld
    END MODULE trfile
