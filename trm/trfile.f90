! trfile.f90

MODULE trfile

  PRIVATE
  PUBLIC tr_save,tr_load,tr_grsv,tr_grld

CONTAINS

!     ***********************************************************

!           SAVE TRANSPORT DATA

!     ***********************************************************

      SUBROUTINE tr_save

      USE trcomm
      USE libfio
      IMPLICIT NONE
      INTEGER:: nnb,nec,nlh,nic,npel,npsc,nnf,ns
      INTEGER:: IERR, NDD, NDM, NDY, NTH1, NTM1, NTS1
      CHARACTER(LEN=3):: K1, K2, K3, K4, K5, K6


      CALL FWOPEN(21,KNAMTR,0,0,'TR',IERR)
      IF(IERR.NE.0) RETURN

      WRITE(21) MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI
      WRITE(21) NRMAX,NRAMAX,NROMAX,DT,NGPST,TSST
      WRITE(21) NTMAX,NTSTEP,NGTSTP,NGRSTP,NSMAX,NSZMAX,NSNMAX
      WRITE(21) NSS,NSV,NNS,NST,NEA,NEQMAX
      WRITE(21) RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA,PHIA
      WRITE(21) PA,PZ,PN,PNS,PT,PTS
      WRITE(21) PNC,PNFE,PNNU,PNNUS
      WRITE(21) PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,PROFJ1,PROFJ2, &
     &          ALP,AD0,AV0,CNP,CNH,CDP,CDH,CDW,CWEB,CALF
      WRITE(21) MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF,MDLCD05
      WRITE(21) TPRST,MDLST,model_nnf,IZERO
      WRITE(21) MODELG,NTEQIT
      WRITE(21) MDLXP,MDLUF,MDLNCL,MDLWLD,MDDIAG,MDLDW,MDLFLX,MDLER
      WRITE(21) KUFDEV,KUFDCG,TIME_INT,MODEP,MDLNI,MDLJQ,MDLTC,MDLPCK
      WRITE(21) nnbmax,necmax,nlhmax,nicmax,npelmax,npscmax,nnfmax
      WRITE(21) (PNBIN(nnb),PNBR0(nnb),PNBRW(nnb),PNBCD(nnb),PNBVY(nnb), &
           PNBVW(nnb),PNBENG(nnb),PNBRTG(nnb), &
           model_nnb(nnb),ns_nnb(nnb),nrmax_nnb(nnb),nnb=1,nnbmax)
      WRITE(21) (PECIN(nec),PECR0(nec),PECRW(nec),PECCD(nec),PECTOE(nec), &
           PECNPR(nec), &
           MDLEC(nec),nec=1,necmax)
      WRITE(21) (PLHIN(nlh),PLHR0(nlh),PLHRW(nlh),PLHCD(nlh),PLHTOE(nlh), &
           PLHNPR(nlh), &
           MDLLH(nlh),nlh=1,nlhmax)
      WRITE(21) (PICIN(nic),PICR0(nic),PICRW(nic),PICCD(nic),PICTOE(nic), &
           PICNPR(nic), &
           MDLIC(nic),nic=1,nicmax)
      WRITE(21) (PELIN(npel),PELR0(npel),PELRW(npel),PELRAD(npel), &
           PELVEL(npel),PELTIM(npel), &
           pellet_time_start(npel),pellet_time_interval(npel), &
           (PELPAT(ns,npel),ns=1,nsm), &
           MDLPEL(npel),npel=1,npelmax)
      WRITE(21) (PSCIN(npsc),PSCR0(npsc),PSCRW(npsc), &
           MDLPSC(npsc),NSPSC(npsc),npsc=1,npscmax)
      WRITE(21) (model_nnf(nnf),ns_nnf(nnf),nnf=1,nnfmax)
      WRITE(21) PBSCD,MDLCD
      WRITE(21) DR,PNSS,T,TST,VSEC,WPPRE,TPRE,KFNLOG,NTMAX_SAVE
      WRITE(21) RG,RM,RN,RT,RU,RW,BP,RDP,RPSI,RNF,RTF,ANC,ANFE,ANNU,RDPVRHOG
      WRITE(21) VTOR,VPOL,AJOH,EZOH,ETA,AMZ
      WRITE(21) EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA,CALF
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

      RETURN
      END SUBROUTINE tr_save

!     ***********************************************************

!           LOAD TRANSPORT DATA

!     ***********************************************************

      SUBROUTINE tr_load

        USE trcomm
      USE libitp
      USE libfio
      IMPLICIT NONE
      INTEGER:: nnb,nec,nlh,nic,npel,npsc,nnf,ns
      INTEGER:: IERR


      CALL FROPEN(21,KNAMTR,0,0,'TR',IERR)
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
     &         ALP,AD0,AV0,CNP,CNH,CDP,CDH,CDW,CWEB,CALF
      READ(21) MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF,MDLCD05
      READ(21) TPRST,MDLST,IZERO
      READ(21) MODELG,NTEQIT
      READ(21) MDLXP,MDLUF,MDLNCL,MDLWLD,MDDIAG,MDLDW,MDLFLX,MDLER
      READ(21) KUFDEV,KUFDCG,TIME_INT,MODEP,MDLNI,MDLJQ,MDLTC,MDLPCK
      READ(21) nnbmax,necmax,nlhmax,nicmax,npelmax,npscmax,nnfmax
      READ(21) (PNBIN(nnb),PNBR0(nnb),PNBRW(nnb),PNBCD(nnb),PNBVY(nnb), &
           PNBVW(nnb),PNBENG(nnb),PNBRTG(nnb), &
           model_nnb(nnb),ns_nnb(nnb),nrmax_nnb(nnb),nnb=1,nnbmax)
      READ(21) (PECIN(nec),PECR0(nec),PECRW(nec),PECCD(nec),PECTOE(nec), &
           PECNPR(nec), &
           MDLEC(nec),nec=1,necmax)
      READ(21) (PLHIN(nlh),PLHR0(nlh),PLHRW(nlh),PLHCD(nlh),PLHTOE(nlh), &
           PLHNPR(nlh), &
           MDLLH(nlh),nlh=1,nlhmax)
      READ(21) (PICIN(nic),PICR0(nic),PICRW(nic),PICCD(nic),PICTOE(nic), &
           PICNPR(nic), &
           MDLIC(nic),nic=1,nicmax)
      READ(21) (PELIN(npel),PELR0(npel),PELRW(npel),PELRAD(npel), &
           PELVEL(npel),PELTIM(npel), &
           pellet_time_start(npel),pellet_time_interval(npel), &
           (PELPAT(ns,npel),ns=1,nsm), &
           MDLPEL(npel),npel=1,npelmax)
      READ(21) (PSCIN(npsc),PSCR0(npsc),PSCRW(npsc), &
           MDLPSC(npsc),NSPSC(npsc),npsc=1,npscmax)
      READ(21) (model_nnf(nnf),ns_nnf(nnf),nnf=1,nnfmax)
      READ(21) PBSCD,MDLCD
      READ(21) DR,PNSS,T,TST,VSEC,WPPRE,TPRE,KFNLOG,NTMAX_SAVE
      READ(21) RG,RM,RN,RT,RU,RW,BP,RDP,RPSI,RNF,RTF,ANC,ANFE,ANNU,RDPVRHOG
      READ(21) VTOR,VPOL,AJOH,EZOH,ETA,AMZ
      READ(21) EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA,CALF
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

      RETURN
      END SUBROUTINE tr_load

!     ***********************************************************

!           SAVE GRAPHIC DATA

!     ***********************************************************

      SUBROUTINE tr_grsv

      USE TRCOMM, ONLY : GRG, GRM, GT, GTR, GVR, GVT, NGR, NGT
      IMPLICIT NONE
      INTEGER:: IST
      CHARACTER(LEN=32):: TRFLNM
      LOGICAL   :: LEX


    1 WRITE(6,*) '# INPUT : GRAPHIC SAVE FILE NAME (CR TO CANCEL)'
      READ(5,'(A32)',ERR=1,END=900) TRFLNM
      IF(TRFLNM.EQ.'                                ') GOTO 900
      INQUIRE (FILE=TRFLNM,EXIST=LEX)
      IF(LEX) THEN
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
      INTEGER:: IST
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
