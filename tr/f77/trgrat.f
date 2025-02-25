C     $Id$
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRT0(K2,INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER K2*1
C
      IF(NT.LT.NGTSTP) RETURN
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
C
      IF(K2.EQ.'1') CALL TRGRT1(INQ)
      IF(K2.EQ.'2') CALL TRGRT2(INQ)
      IF(K2.EQ.'5') CALL TRGRT5(INQ)
      IF(K2.EQ.'6') CALL TRGRT6(INQ)
      IF(K2.EQ.'7') CALL TRGRT7(INQ)
      IF(K2.EQ.'8') CALL TRGRT8(INQ)
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
C
      RETURN
      END
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRX0(K2,INQ)
C
      CHARACTER K2*1
C
      IF(K2.EQ.'1') CALL TRGRX1(INQ)
      RETURN
      END
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRA0(K2,INQ)
C
      CHARACTER K2*1
C
      IF(K2.EQ.'1') THEN
         CALL TRGRT6(INQ)
         CALL TRGRT7(INQ)
         CALL TRGRX1(INQ)
         CALL TRGRR1(INQ)
         CALL TRGRR2(INQ)
         CALL TRGRR4(INQ)
         CALL TRGRR7(INQ)
      ELSEIF(K2.EQ.'2') THEN
         CALL TRGRT6(INQ)
         CALL TRGRT7(INQ)
         CALL TRGRX1(INQ)
         CALL TRGRG1(INQ)
         CALL TRGRG3(INQ)
         CALL TRGRR4(INQ)
         CALL TRGRG5(INQ)
      ENDIF
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TEMPORAL CONTOUR
C                     TE, TD, NE, AJ, TE(XP), TD(XP), AJBS, AJCD
C
C     ***********************************************************
C
      SUBROUTINE TRGRT1(INQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION GYL(NTM,NRM),KATR(8,NTM,NRM)
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      CALL PAGES
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I, 1)
         ENDDO
      ENDDO
      CALL TRGR1DC( 3.0,12.0,14.0,17.0,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@TE [keV]  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I, 2)
         ENDDO
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@TD [keV]  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I, 5)
         ENDDO
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 5.4, 8.4,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@NE [10$+20$=/m$+3$=]  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I, 9) * 1.E-6
         ENDDO
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 1.1, 4.1,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@AJ [MA]  vs t@',KATR)
C
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         DO I=1,NGT
            TSL = DBLE(GT(I))
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,RTEL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
               GYL(I,NR) = GUCLIP(RTEL)
            ENDDO
         ENDDO
         CALL TRGR1DC(15.0,24.0,14.0,17.0,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &                '@TE [keV] (XP) vs t@',KATR)
C
         DO I=1,NGT
            TSL = DBLE(GT(I))
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,RTDL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
               GYL(I,NR) = GUCLIP(RTDL)
            ENDDO
         ENDDO
         CALL TRGR1DC(15.0,24.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &                '@TD [keV] (XP) vs t@',KATR)
      ENDIF
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I,13) * 1.E-6
         ENDDO
      ENDDO
      CALL TRGR1DC(15.0,24.0, 5.4, 8.4,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@AJBS [MA]  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = (G3D(NR,I,11) + G3D(NR,I,12)) * 1.E-6
         ENDDO
      ENDDO
      CALL TRGR1DC(15.0,24.0, 1.1, 4.1,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@AJCD [MA]  vs t@',KATR)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TEMPORAL CONTOUR
C                     QP, VTOR, VLOOP, ER
C
C     ***********************************************************
C
      SUBROUTINE TRGRT2(INQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION GYL(NTM,NRM),KATR(8,NTM,NRM)
C
      CALL PAGES
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I,27)
         ENDDO
      ENDDO
      CALL TRGR1DC( 3.0,12.0,14.0,17.0,GT,GRG,GYL,NTM,NGT,NRM,NRMAX,
     &            '@QP  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I,31)
         ENDDO
      ENDDO
      CALL TRGR1DC(15.0,24.0,14.0,17.0,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@VLOOP [V]  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I,60)
         ENDDO
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@VTOR [m/s]  vs t@',KATR)
C
      DO I=1,NGT
         DO NR=1,NRMAX
            GYL(I,NR) = G3D(NR,I,63)
         ENDDO
      ENDDO
      CALL TRGR1DC(15.0,24.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRM,NRMAX,
     &            '@ER [V/m] vs t@',KATR)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME DEPENDENCE : RR,RA,RKAP,IP,BB,H98Y2
C                     TIME DEPENDENCE : TRANSIENT
C
C     ***********************************************************
C
      SUBROUTINE TRGRT5(INQ)
C
      INCLUDE 'trcomm.inc'
C
      IF(NGST.EQ.0) THEN
C
      CALL PAGES
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,97)
         GYT(I,2)=GVT(I,98)
      ENDDO
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,2,
     &            '@RR, RA [M]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,99)
         GYT(I,2)=GVT(I,34)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,2,
     &            '@BT, IP  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,100)
      ENDDO
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@RKAP  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,81)
         GYT(I,2)=GVT(I,102)
         GYT(I,3)=GVT(I,103)
      ENDDO
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,3,
     &            '@TAUE89,TAUE98 [s], H98Y2  vs t@',2+INQ)
C
      CALL PAGEE
C
      ELSE
C
    1 WRITE(6,*) ' CHOSE ONE (TE=1,TD=2,TT=3,TA=4) '
      READ(5,*,ERR=1,END=900) M
C      
      CALL PAGES
C
      IX = INT((NRMAX+1-IZERO)/(NGPST-1))
      GW = 150.0/(15.0*NGPST/2.0-5.0)
      GD = 0.5*GW
C      IR = INT((TSST/DT)/NGTSTP)
C
      DO N=0,NGPST/2-1
         K=4*N
C
         DO I=1,NGST
            GYT(I,1) = GVT(I,K+89)
            GYT(I,2) = GVT(I,K+90)
            GYT(I,3) = GVT(I,K+91)
            GYT(I,4) = GVT(I,K+92)
C
         ENDDO
         CALL TRGR1D(3.0,12.0,17-(N+1)*GW-N*GD,17-N*(GW+GD),
     &            GTS,GYT,NTM,
     &            NGST,M,'@TE,TD,TT,TA [keV]@',1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(10.0,17.1-N*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)
C        
      ENDDO
C
      DO N=NGPST/2,NGPST-1
         K=4*N
C
         DO I=1,NGST
            GYT(I,1) = GVT(I,K+89)
            GYT(I,2) = GVT(I,K+90)
            GYT(I,3) = GVT(I,K+91)
            GYT(I,4) = GVT(I,K+92)
C
         ENDDO
         NGPSTH=N-NGPST/2
         CALL TRGR1D(15.5,24.5,17-(NGPSTH+1)*GW-NGPSTH*GD,
     &               17-NGPSTH*(GW+GD),
     &               GTS,GYT,NTM,NGST,M,'@TE,TD,TT,TA [keV]  vs t@',
     &               1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(22.5,17.1-NGPSTH*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)
C        
      ENDDO
C
      CALL PAGEE
C
      ENDIF
C
  900 RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME DEPENDENCE : TEMPERATURE,CURRENT,
C                                       INPUT POWER,OUTPUT POWER
C
C     ***********************************************************
C
      SUBROUTINE TRGRT6(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO I=1,NGT
         GYT(I,1)=GVT(I, 9)
         GYT(I,2)=GVT(I,10)
         GYT(I,3)=GVT(I,13)
         GYT(I,4)=GVT(I,14)
      ENDDO
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@TE,TD,<TE>,<TD> [keV]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,34)
         GYT(I,2)=GVT(I,35)
         GYT(I,3)=GVT(I,36)
         GYT(I,4)=GVT(I,37)
         GYT(I,5)=GVT(I,38)
         GYT(I,6)=GVT(I,101)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,6,
     &            '@IPPRL,IOH,INB,IRF,IBS,IP [MA]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,39)
         GYT(I,2)=GVT(I,40)
         GYT(I,3)=GVT(I,41)
         GYT(I,4)=GVT(I,42)+GVT(I,43)+GVT(I,44)+GVT(I,45)
         GYT(I,5)=GVT(I,46)
         GYT(I,6)=GVT(I,89)+GVT(I,90)+GVT(I,91)+GVT(I,92)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,6,
     &            '@PIN,POH,PNB,PRF,PNF,PEX [MW]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,57)
         GYT(I,2)=GVT(I,58)
         GYT(I,3)=GVT(I,59)
         GYT(I,4)=GVT(I,60)
         GYT(I,5)=GVT(I,61)+GVT(I,62)+GVT(I,63)+GVT(I,64)
C     &           +GVT(I,89)+GVT(I,90)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,5,
     &            '@POUT,PCX,PIE,PRL,PCON [MW]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,87)
      ENDDO
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@QF  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,13)
         GYT(I,2)=GVT(I,14)
         GYT(I,3)=GVT(I,15)
         GYT(I,4)=GVT(I,16)
      ENDDO
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,4,
     &            '@TEAV,TDAV,TTAV,TAAV [keV]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I, 9)
         GYT(I,2)=GVT(I,10)
         GYT(I,3)=GVT(I,11)
         GYT(I,4)=GVT(I,12)
      ENDDO
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,4,
     &            '@TE0,TD0,TT0,TA0 [keV]  vs t@',2+INQ)
C
      IF(MDLNF.EQ.0) THEN
      DO I=1,NGT
         GYT(I,1)=GVT(I,1)
         GYT(I,2)=GVT(I,5)
         GYT(I,3)=GVT(I,104)
      ENDDO
      CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,3,
     &            '@NE0,<NE>,<NEL> [10$+20$=/m$+3$=]  vs t@',2+INQ)
      ELSE
      DO I=1,NGT
         GYT(I,1)=GVT(I,1)
         GYT(I,2)=GVT(I,2)
         GYT(I,3)=GVT(I,3)
         GYT(I,4)=GVT(I,4)
      ENDDO
      CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,4,
     &            '@NE0,ND0,NT0,NA0 [10$+20$=/m$+3$=]  vs t@',2+INQ)
      ENDIF
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME DEPENDENCE : TAUE,QF,BETAP,BETA
C
C     ***********************************************************
C
      SUBROUTINE TRGRT7(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,33)
         GYT(I,2)=GVT(I,31)+GVT(I,29)
         GYT(I,3)=GVT(I,31)
         GYT(I,4)=GVT(I,17)
      ENDDO
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@WF,WB,WI,WE [MJ]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,85)*100.0
         GYT(I,2)=GVT(I,84)*100.0
         GYT(I,3)=GVT(I,85)*100.0/(GVT(I,34)/GUCLIP(RA*BB))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,3,
     &            '@BETAa,BETA0,[%],BETAN  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,79)
         GYT(I,2)=GVT(I,80)
         GYT(I,3)=GVT(I,81)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,3,
     &            '@TAUE1,TAUE2,TAUE89  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,83)
         GYT(I,2)=GVT(I,82)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,2,
     &            '@BETAPa,BETAP0  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,74)
      ENDDO
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@VLOOP [V]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,77)
      ENDDO
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,1,
     &            '@Q(0)  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,75)
      ENDDO
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,1,
     &            '@ALI  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,76)
      ENDDO
      CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,1,
     &            '@RQ1 [m]  vs t@',2+INQ)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME DEPENDENCE : TRANSIENT
C
C     ***********************************************************
C
      SUBROUTINE TRGRT8(INQ)
C
      INCLUDE 'trcomm.inc'
C
      IF(NGST.EQ.0) RETURN
C
    1 WRITE(6,*) ' CHOSE ONE (TE=1,TD=2,TT=3,TA=4) '
      READ(5,*,ERR=1,END=900) M
C      
      CALL PAGES
C
      IX = INT((NRMAX+1-IZERO)/(NGPST-1))
      GW = 150.0/(15.0*NGPST/2.0-5.0)
      GD = 0.5*GW
C      IR = INT((TSST/DT)/NGTSTP)
C
      DO N=0,NGPST/2-1
         K=4*N
C
         DO I=1,NGST
            GYT(I,1) = GVT(I,K+89)
            GYT(I,2) = GVT(I,K+90)
            GYT(I,3) = GVT(I,K+91)
            GYT(I,4) = GVT(I,K+92)
C
         ENDDO
         CALL TRGR1D(3.0,12.0,17-(N+1)*GW-N*GD,17-N*(GW+GD),
     &            GTS,GYT,NTM,
     &            NGST,M,'@TE,TD,TT,TA [keV]@',1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(10.0,17.1-N*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)
C        
      ENDDO
C
      DO N=NGPST/2,NGPST-1
         K=4*N
C
         DO I=1,NGST
            GYT(I,1) = GVT(I,K+89)
            GYT(I,2) = GVT(I,K+90)
            GYT(I,3) = GVT(I,K+91)
            GYT(I,4) = GVT(I,K+92)
C
         ENDDO
         NGPSTH=N-NGPST/2
         CALL TRGR1D(15.5,24.5,17-(NGPSTH+1)*GW-NGPSTH*GD,
     &               17-NGPSTH*(GW+GD),
     &               GTS,GYT,NTM,NGST,M,'@TE,TD,TT,TA [keV]  vs t@',
     &               1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(22.5,17.1-NGPSTH*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)
C        
      ENDDO
C
      CALL PAGEE
C
  900 RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME LISSAGE : tau vs eps*betap,li,
C
C     ***********************************************************
C
      SUBROUTINE TRGRX1(INQ)
C
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO I=1,NGT
         GYT(I,1)=GUCLIP(RA/RR)*GVT(I,83)
         GYT(I,2)=GVT(I,80)/GVT(I,81)
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@tauE/tauE89 vs eps*betap@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,75)
         GYT(I,2)=GVT(I,80)/GVT(I,81)
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@tauE/tauE89 vs li@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GUCLIP(RA/RR)*GVT(I,83)
         GYT(I,2)=GVT(I,38)/GVT(I,34)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@Ibs/Ip vs eps*betap@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,77)
         GYT(I,2)=GVT(I,80)/GVT(I,81)
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@tauE/tauE89 vs q0@',2+INQ)
C
      CALL PAGEE
C
      RETURN
      END
