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
      CHARACTER K2*1
C
      IF(K2.EQ.'1') CALL TRGRT1(INQ)
      IF(K2.EQ.'2') CALL TRGRT2(INQ)
      IF(K2.EQ.'3') CALL TRGRT3(INQ)
      IF(K2.EQ.'4') CALL TRGRT4(INQ)
      IF(K2.EQ.'5') CALL TRGRT5(INQ)
      IF(K2.EQ.'6') CALL TRGRT6(INQ)
      IF(K2.EQ.'7') CALL TRGRT7(INQ)
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
         CALL TRGRT1(INQ)
         CALL TRGRT3(INQ)
         CALL TRGRT4(INQ)
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
C           GRAPHIC : TIME DEPENDENCE : TEMPERATURE,CURRENT,
C                                       INPUT POWER,OUTPUT POWER
C
C     ***********************************************************
C
      SUBROUTINE TRGRT1(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO I=1,NGT
         GYT(I,1)=GVT(I, 9)
         GYT(I,2)=GVT(I,10)
         GYT(I,3)=GVT(I,13)
         GYT(I,4)=GVT(I,14)
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@TE,TD,<TE>,<TD> [keV]  vs t@',2+INQ)
C
      DO I=1,NGT
         GYT(I,1)=GVT(I,34)
         GYT(I,2)=GVT(I,35)
         GYT(I,3)=GVT(I,36)
         GYT(I,4)=GVT(I,37)
         GYT(I,5)=GVT(I,38)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GT,GYT,NTM,NGT,5,
     &            '@IP,IOH,INB,IRF,IBS [MA]  vs t@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GVT(I,39)
         GYT(I,2)=GVT(I,40)
         GYT(I,3)=GVT(I,41)
         GYT(I,4)=GVT(I,42)+GVT(I,43)+GVT(I,44)+GVT(I,45)+GVT(I,89)
         GYT(I,5)=GVT(I,46)
  300 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GT,GYT,NTM,NGT,5,
     &            '@PIN,POH,PNB,PRF,PNF [MW]  vs t@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,57)
         GYT(I,2)=GVT(I,58)
         GYT(I,3)=GVT(I,59)
         GYT(I,4)=GVT(I,60)
         GYT(I,5)=GVT(I,61)+GVT(I,62)+GVT(I,63)+GVT(I,64)
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GT,GYT,NTM,NGT,5,
     &            '@POUT,PCX,PIE,PRL,PCON [MW]  vs t@',2+INQ)
C
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME DEPENDENCE : STORED ENERGY,TE,TD,NE
C
C     ***********************************************************
C
      SUBROUTINE TRGRT2(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 I=1,NGT
         GYT(I,1)=GVT(I,87)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@QF  vs t@',2+INQ)
C
C      DO 500 I=1,NGT
C         GYT(I,1)=GVT(I,33)
C         GYT(I,2)=GVT(I,31)+GVT(I,29)
C         GYT(I,3)=GVT(I,31)
C         GYT(I,4)=GVT(I,17)
C  500 CONTINUE
C      CALL TRGR1D( 3.0,12.0,11.0,17.0,GT,GYT,NTM,NGT,4,
C     &            '@WF,WB,WI,WE [MJ]  vs t@',2+INQ)
C
      DO 200 I=1,NGT
         GYT(I,1)=GVT(I, 9)
         GYT(I,2)=GVT(I,10)
         GYT(I,3)=GVT(I,11)
         GYT(I,4)=GVT(I,12)
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@TE0,TD0,TT0,TA0 [keV]  vs t@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GVT(I,13)
         GYT(I,2)=GVT(I,14)
         GYT(I,3)=GVT(I,15)
         GYT(I,4)=GVT(I,16)
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GT,GYT,NTM,NGT,4,
     &            '@TEAV,TDAV,TTAV,TAAV [keV]  vs t@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,1)
         GYT(I,2)=GVT(I,5)
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GT,GYT,NTM,NGT,2,
     &            '@NE0,<NE> [10^20/m^3]  vs t@',2+INQ)
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
      SUBROUTINE TRGRT3(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 I=1,NGT
         GYT(I,1)=GVT(I,33)
         GYT(I,2)=GVT(I,31)+GVT(I,29)
         GYT(I,3)=GVT(I,31)
         GYT(I,4)=GVT(I,17)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@WF,WB,WI,WE [MJ]  vs t@',2+INQ)
C
      DO 200 I=1,NGT
         GYT(I,1)=GVT(I,79)
         GYT(I,2)=GVT(I,80)
         GYT(I,3)=GVT(I,81)
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GT,GYT,NTM,NGT,3,
     &            '@TAUE1,TAUE2,TAUEP [s]  vs t@',2+INQ)
C
C      DO 200 I=1,NGT
C         GYT(I,1)=GVT(I,87)
C  200 CONTINUE
C      CALL TRGR1D(15.5,24.5,11.0,17.0,GT,GYT,NTM,NGT,1,
C     &            '@QF  vs t@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GVT(I,85)
         GYT(I,2)=GVT(I,84)
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GT,GYT,NTM,NGT,2,
     &            '@BETAa,BETA0  vs t@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,83)
         GYT(I,2)=GVT(I,82)
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GT,GYT,NTM,NGT,2,
     &            '@BETAPa,BETAP0  vs t@',2+INQ)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : TIME DEPENDENCE : VLOOP,ALI,Q0,RQ1
C
C     ***********************************************************
C
      SUBROUTINE TRGRT4(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 I=1,NGT
         GYT(I,1)=GVT(I,74)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@VLOOP [V]  vs t@',2+INQ)
C
      DO 200 I=1,NGT
         GYT(I,1)=GVT(I,75)
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@ALI  vs t@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GVT(I,77)
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GT,GYT,NTM,NGT,1,
     &            '@Q(0)  vs t@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,76)
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GT,GYT,NTM,NGT,1,
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
      SUBROUTINE TRGRT5(INQ)
C
      INCLUDE 'trcomm.h'
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
      DO 10 N=0,NGPST/2-1
         K=4*N
C
         DO 100 I=1,NGST
            GYT(I,1) = GVT(I,K+89)
            GYT(I,2) = GVT(I,K+90)
            GYT(I,3) = GVT(I,K+91)
            GYT(I,4) = GVT(I,K+92)
C
  100    CONTINUE
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
   10 CONTINUE
C
      DO 20 N=NGPST/2,NGPST-1
         K=4*N
C
         DO 200 I=1,NGST
            GYT(I,1) = GVT(I,K+89)
            GYT(I,2) = GVT(I,K+90)
            GYT(I,3) = GVT(I,K+91)
            GYT(I,4) = GVT(I,K+92)
C
  200    CONTINUE
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
   20 CONTINUE
C
      CALL PAGEE
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
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 I=1,NGT
         GYT(I,1)=GVT(I, 9)
         GYT(I,2)=GVT(I,10)
         GYT(I,3)=GVT(I,13)
         GYT(I,4)=GVT(I,14)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@TE,TD,<TE>,<TD> [keV]  vs t@',2+INQ)
C
      DO 200 I=1,NGT
         GYT(I,1)=GVT(I,34)
         GYT(I,2)=GVT(I,35)
         GYT(I,3)=GVT(I,36)
         GYT(I,4)=GVT(I,37)
         GYT(I,5)=GVT(I,38)
  200 CONTINUE
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,5,
     &            '@IP,IOH,INB,IRF,IBS [MA]  vs t@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GVT(I,39)
         GYT(I,2)=GVT(I,40)
         GYT(I,3)=GVT(I,41)
         GYT(I,4)=GVT(I,42)+GVT(I,43)+GVT(I,44)+GVT(I,45)
         GYT(I,5)=GVT(I,46)
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,5,
     &            '@PIN,POH,PNB,PRF,PNF [MW]  vs t@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,57)
         GYT(I,2)=GVT(I,58)
         GYT(I,3)=GVT(I,59)
         GYT(I,4)=GVT(I,60)
         GYT(I,5)=GVT(I,61)+GVT(I,62)+GVT(I,63)+GVT(I,64)
  400 CONTINUE
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,5,
     &            '@POUT,PCX,PIE,PRL,PCON [MW]  vs t@',2+INQ)
C
      DO 1100 I=1,NGT
         GYT(I,1)=GVT(I,87)
 1100 CONTINUE
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@QF  vs t@',2+INQ)
C
      DO 1300 I=1,NGT
         GYT(I,1)=GVT(I,13)
         GYT(I,2)=GVT(I,14)
         GYT(I,3)=GVT(I,15)
         GYT(I,4)=GVT(I,16)
 1300 CONTINUE
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,4,
     &            '@TEAV,TDAV,TTAV,TAAV [keV]  vs t@',2+INQ)
C
      DO 1200 I=1,NGT
         GYT(I,1)=GVT(I, 9)
         GYT(I,2)=GVT(I,10)
         GYT(I,3)=GVT(I,11)
         GYT(I,4)=GVT(I,12)
 1200 CONTINUE
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,4,
     &            '@TE0,TD0,TT0,TA0 [keV]  vs t@',2+INQ)
C
      IF(MDLNF.EQ.0) THEN
      DO 1400 I=1,NGT
         GYT(I,1)=GVT(I,1)
         GYT(I,2)=GVT(I,5)
 1400 CONTINUE
      CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,2,
     &            '@NE0,<NE> [10^20/m^3]  vs t@',2+INQ)
      ELSE
      DO 1410 I=1,NGT
         GYT(I,1)=GVT(I,1)
         GYT(I,2)=GVT(I,2)
         GYT(I,3)=GVT(I,3)
         GYT(I,4)=GVT(I,4)
 1410 CONTINUE
      CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,4,
     &            '@NE0,ND0,NT0,NA0 [10^20/m^3]  vs t@',2+INQ)
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
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 I=1,NGT
         GYT(I,1)=GVT(I,33)
         GYT(I,2)=GVT(I,31)+GVT(I,29)
         GYT(I,3)=GVT(I,31)
         GYT(I,4)=GVT(I,17)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@WF,WB,WI,WE [MJ]  vs t@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GVT(I,85)*100.0
         GYT(I,2)=GVT(I,84)*100.0
         GYT(I,3)=GVT(I,85)*100.0/(GVT(I,34)/GCLIP(RA*BB))
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,3,
     &            '@BETAa,BETA0,[%],BETAN  vs t@',2+INQ)
C
      DO 200 I=1,NGT
         GYT(I,1)=GVT(I,79)
         GYT(I,2)=GVT(I,80)
         GYT(I,3)=GVT(I,81)
  200 CONTINUE
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,3,
     &            '@TAUE1,TAUE2,TAUEP [s]  vs t@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,83)
         GYT(I,2)=GVT(I,82)
  400 CONTINUE
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,2,
     &            '@BETAPa,BETAP0  vs t@',2+INQ)
C
      DO 1100 I=1,NGT
         GYT(I,1)=GVT(I,74)
 1100 CONTINUE
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,
     &            '@VLOOP [V]  vs t@',2+INQ)
C
      DO 1300 I=1,NGT
         GYT(I,1)=GVT(I,77)
 1300 CONTINUE
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,1,
     &            '@Q(0)  vs t@',2+INQ)
C
      DO 1200 I=1,NGT
         GYT(I,1)=GVT(I,75)
 1200 CONTINUE
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,1,
     &            '@ALI  vs t@',2+INQ)
C
      DO 1400 I=1,NGT
         GYT(I,1)=GVT(I,76)
 1400 CONTINUE
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
C           GRAPHIC : TIME LISSAGE : tau vs eps*betap,li,
C
C     ***********************************************************
C
      SUBROUTINE TRGRX1(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 I=1,NGT
         GYT(I,1)=GCLIP(RA/RR)*GVT(I,83)
         GYT(I,2)=GVT(I,80)/GVT(I,81)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@tauE/tauE89 vs eps*betap@',2+INQ)
C
      DO 200 I=1,NGT
         GYT(I,1)=GVT(I,75)
         GYT(I,2)=GVT(I,80)/GVT(I,81)
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@tauE/tauE89 vs li@',2+INQ)
C
      DO 300 I=1,NGT
         GYT(I,1)=GCLIP(RA/RR)*GVT(I,83)
         GYT(I,2)=GVT(I,38)/GVT(I,34)
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@Ibs/Ip vs eps*betap@',2+INQ)
C
      DO 400 I=1,NGT
         GYT(I,1)=GVT(I,77)
         GYT(I,2)=GVT(I,80)/GVT(I,81)
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GYT(1,1),GYT(1,2),NTM,NGT,1,
     &            '@tauE/tauE89 vs q0@',2+INQ)
C
      CALL PAGEE
C
      RETURN
      END
