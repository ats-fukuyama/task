C     $Id$
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRR0(K2,INQ)
C
      CHARACTER K2*1
C
      IF(K2.EQ.'1') CALL TRGRR1(INQ)
      IF(K2.EQ.'2') CALL TRGRR2(INQ)
      IF(K2.EQ.'3') CALL TRGRR3(INQ)
      IF(K2.EQ.'4') CALL TRGRR4(INQ)
      IF(K2.EQ.'5') CALL TRGRR5(INQ)
      IF(K2.EQ.'6') CALL TRGRR6(INQ)
      IF(K2.EQ.'7') CALL TRGRR7(INQ)
      IF(K2.EQ.'8') CALL TRGRR8(INQ)
      IF(K2.EQ.'9') CALL TRGRR9(INQ)
      RETURN
      END
C  
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRY0(K2,INQ)
C
      CHARACTER K2*1
C
      IF(K2.EQ.'1') CALL TRGRY1(INQ)
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : RNS,RNF,RTS,RTF
C
C     ***********************************************************
C
      SUBROUTINE TRGRR1(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NS=1,NSM
      DO 100 NR=1,NRMAX
         GYR(NR,NS) = GCLIP(RN(NR,NS))
  100 CONTINUE
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@NE,ND [10^20/m^3]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
     &               '@NE,ND,NT,NA [10^20/m^3]  vs r@',2+INQ)
      ENDIF
      DO 200 NF=1,NFM
      DO 200 NR=1,NRMAX
         GYR(NR,NF) = GCLIP(RNF(NR,NF))
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,
     &            '@NB,NF [10^20 /m^3]  vs r@',2+INQ)
C
      DO 300 NS=1,NSM
      DO 300 NR=1,NRMAX
         GYR(NR,NS) = GCLIP(RT(NR,NS))
  300 CONTINUE
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TE,TD, [keV]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,
     &               '@TE,TD,TT,TA [keV]  vs r@',2+INQ)
      ENDIF
C
      DO 400 NS=1,NSM
      DO 400 NR=1,NRMAX
         GYR(NR,NS) =  GCLIP(RN(NR,NS)*RT(NR,NS)*RKEV*1.D14)
  400 CONTINUE
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@PE,PD [MPa]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,
     &               '@PE,PD,PT,PA [MPa]  vs r@',2+INQ)
      ENDIF
      CALL TRGRTM
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : WS,WF,BETA,BETAP
C
C     ***********************************************************
C
      SUBROUTINE TRGRR2(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NR=1,NRMAX
         GYR(NR,1) = GCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GCLIP(PNB(NR) * 1.D-6)
         GYR(NR,3) = GCLIP(PNF(NR) * 1.D-6)
  100 CONTINUE
      DO 110 NS=1,NSM
      DO 110 NR=1,NRMAX
         GYR(NR,NS+3) = GCLIP(PRF(NR,NS) * 1.D-6)
  110 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM+3,
     &            '@POH,PNB,PNF,PRF [MW/m^3]  vs r@',2+INQ)
C
      DO 200 NR=1,NRMAX
         GYR(NR,1) = GCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GCLIP(PRL(NR) * 1.D-6)
         GYR(NR,3) = GCLIP(PCX(NR) * 1.D-6)
         GYR(NR,4) = GCLIP(PIE(NR) * 1.D-6)
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,
     &            '@POH,PRL,PCX,PIE [MW/m^3]  vs r@',2+INQ)
C
      DO 300 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(BETAL(NR))
         GYR(NR+1,2) = GCLIP(BETA(NR))
  300 CONTINUE
      GYR(1,1)=GCLIP(BETA0)
      GYR(1,2)=GCLIP(BETA0)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &            '@BETA,<BETA>  vs r@',2+INQ)
C
      DO 400 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(BETAPL(NR))
         GYR(NR+1,2) = GCLIP(BETAP(NR))
         GYR(NR+1,3) = GCLIP(BETAQ(NR))
  400 CONTINUE
      GYR(1,1)=GCLIP(BETAP0)
      GYR(1,2)=GCLIP(BETAP0)
      GYR(1,3)=GCLIP(BETAQ0)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,
     &            '@BETAP,<BETAP>,BETAP1  vs r@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : POWER
C
C     ***********************************************************
C
      SUBROUTINE TRGRR3(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NF=1,NFM
      DO 100 NR=1,NRMAX
         GYR(NR,NF) = GCLIP(RTF(NR,NF))
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,
     &            '@TB,TF [keV]  vs r@',2+INQ)
C
      DO 200 NF=1,NFM
      DO 200 NR=1,NRMAX
         GYR(NR,NF) = GCLIP(1.5D0*RW(NR,NF)*RKEV*1.D14)
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,
     &            '@WB,WF [MJ]  vs r@',2+INQ)
C
      DO 300 NR=1,NRMAX
         GYR(NR,1) = GCLIP(PBIN(NR)   * 1.D-6)
  300 CONTINUE
      DO 310 NS=1,NSM
      DO 310 NR=1,NRMAX
         GYR(NR,NS+1) = GCLIP(PBCL(NR,NS)   * 1.D-6)
  310 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM+1,
     &            '@PBIN,PBCL [MW/m^3]  vs r@',2+INQ)
C
      DO 400 NR=1,NRMAX
         GYR(NR,1) = GCLIP(PFIN(NR) * 1.D-6)
  400 CONTINUE
      DO 410 NS=1,NSM
      DO 410 NR=1,NRMAX
         GYR(NR,NS+1) = GCLIP(PFCL(NR,NS) * 1.D-6)
  410 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM+1,
     &            '@PFIN,PFCL [MW/m^3]  vs r@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : Q,EZ,AJ,S
C
C     ***********************************************************
C
      SUBROUTINE TRGRR4(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(QP(NR))
  100 CONTINUE
      GYR(1,1) = GCLIP(Q0)
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,
     &            '@QP  vs r@',2+INQ)
C
      DO 200 NR=1,NRMAX
         GYR(NR,1) = GCLIP(EZOH(NR))
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@EZOH [V/m]  vs r@',2+INQ)
C
      DO 300 NR=1,NRMAX
         GYR(NR,1) = GCLIP(AJ(NR)   * 1.D-6)
         GYR(NR,2) = GCLIP(AJOH(NR) * 1.D-6)
         GYR(NR,3) = GCLIP(AJNB(NR) * 1.D-6)
         GYR(NR,4) = GCLIP(AJRF(NR) * 1.D-6)
         GYR(NR,5) = GCLIP(AJBS(NR) * 1.D-6)
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,5,
     &            '@JTOT,JOH,JNB,JRF,JBS [MA/m^2]  vs r@',2+INQ)
C
      DO 400 NR=1,NRMAX
         GYR(NR,1) = GLOG(ETA(NR),1.D-10,1.D0)
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@LOG:ETA  vs r @',11+INQ)
C
C      DO 400 NR=1,NRMAX
C         GYR(NR,1) = GCLIP(SIE(NR))
C         GYR(NR,2) = GCLIP(SNB(NR))
C         GYR(NR,3) = GCLIP(SNF(NR))
C  400 CONTINUE
C      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,3,
C     &            '@SIE,SNB,SNF [/sm^3]  vs r@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : ZEFF,IMPURITY
C
C     ***********************************************************
C
      SUBROUTINE TRGRR5(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NR=1,NRMAX
         GYR(NR,1) = GCLIP(ZEFF(NR))
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ZEFF  vs r@',2+INQ)
C
      DO 200 NR=1,NRMAX
         GYR(NR,1) = GCLIP(PZC(NR))
         GYR(NR,2) = GCLIP(PZFE(NR))
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@PZC,PZFE  vs r@',2+INQ)
C
      DO 300 NR=1,NRMAX
         GYR(NR,1) = GCLIP(ANC(NR))
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ANC [10^20/m^3]  vs r@',2+INQ)
C
      DO 400 NR=1,NRMAX
         GYR(NR,1) = GCLIP(ANFE(NR))
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ANFE [10^20/m^3]  vs r@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : PIN,SSIN,PELLET
C
C     ***********************************************************
C
      SUBROUTINE TRGRR6(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NS=1,NSM
      DO 100 NR=1,NRMAX
         GYR(NR,NS) = GCLIP(PIN(NR,NS)*1.D-6)
  100 CONTINUE
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
     &            '@PIN [MW/m^3]  vs r@',2+INQ)
C
      DO 200 NS=1,NSM
      DO 200 NR=1,NRMAX
         GYR(NR,NS) = GCLIP(SSIN(NR,NS))
  200 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,
     &            '@SSIN [/sm^3]  vs r@',2+INQ)
C
      DO 300 NS=1,NSM
      DO 300 NR=1,NRMAX
         GYR(NR,NS) = GCLIP(SPE(NR,NS))
  300 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
     &            '@SPE [/m^3]  vs r@',2+INQ)
C
      CALL TRGRTM
C
      CALL MOVE(17.5,4.0)
      CALL TEXT('PELVEL=',7)
      CALL NUMBD(PELVEL,'(1PE10.3)',10)
      CALL TEXT('[m/s]',5)
      CALL MOVE(17.5,3.0)
      CALL TEXT('PELRAD=',7)
      CALL NUMBD(PELRAD,'(1PE10.3)',10)
      CALL TEXT('[m]',3)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : ETA,KAI
C
C     ***********************************************************
C
      SUBROUTINE TRGRR7(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(VGR1(NR,1))
         GYR(NR+1,2) = GCLIP(VGR1(NR,2))
         GYR(NR+1,3) = GCLIP(VGR1(NR,3))
  100 CONTINUE
      GYR(1,1) = GCLIP((4.D0*VGR1(1,1)-VGR1(2,1))/3.D0)
      GYR(1,2) = 0.0
      GYR(1,3) = 0.0
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,3,
     &           '@G,s,alpha vs r@',2+INQ)
C
      DO 200 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(AK(NR,1))
         GYR(NR+1,2) = GCLIP(AK(NR,2))
  200 CONTINUE
      GYR(1,1) = GCLIP(AK(1,1))
      GYR(1,2) = GCLIP(AK(1,2))
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,2,
     &           '@AKE,AKI vs r @',2+INQ)
C
      DO 300 NR=1,NRMAX-1
C         GYR(NR+1,1) = GLOG(AK  (NR,1),1.D-2,1.D2)
C         GYR(NR+1,2) = GLOG(AKNC(NR,1),1.D-2,1.D2)
C         GYR(NR+1,3) = GLOG(AKDW(NR,1),1.D-2,1.D2)
         GYR(NR+1,1) = GCLIP(AK(NR,1))
         GYR(NR+1,2) = GCLIP(AKNC(NR,1))
         GYR(NR+1,3) = GCLIP(AKDW(NR,1))
  300 CONTINUE
         GYR(1,1) = GCLIP(AK(1,1))
         GYR(1,2) = GCLIP(AKNC(1,1))
         GYR(1,3) = GCLIP(AKDW(1,1))
C         GYR(NRMAX+1,1) = GCLIP(AK(NRMAX,1)/2)
C         GYR(NRMAX+1,2) = GCLIP(AKNC(NRMAX,1)/2)
C         GYR(NRMAX+1,3) = GCLIP(AKDW(NRMAX,1)/2)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,
     &            '@AKE,AKNCE,AKDWE [m^2/s]  vs r@',2+INQ)
C
      DO 400 NR=1,NRMAX-1
C         GYR(NR+1,1) = GLOG(AK  (NR,2),1.D-2,1.D2)
C         GYR(NR+1,2) = GLOG(AKNC(NR,2),1.D-2,1.D2)
C         GYR(NR+1,3) = GLOG(AKDW(NR,2),1.D-2,1.D2)
         GYR(NR+1,1) = GCLIP(AK(NR,2))
         GYR(NR+1,2) = GCLIP(AKNC(NR,2))
         GYR(NR+1,3) = GCLIP(AKDW(NR,2))
  400 CONTINUE
         GYR(1,1) = GCLIP(AK(1,2))
         GYR(1,2) = GCLIP(AKNC(1,2))
         GYR(1,3) = GCLIP(AKDW(1,2))
C         GYR(NRMAX+1,1) = GCLIP(AK(NRMAX,2)/2)
C         GYR(NRMAX+1,2) = GCLIP(AKNC(NRMAX,2)/2)
C         GYR(NRMAX+1,3) = GCLIP(AKDW(NRMAX,2)/2)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMp,NRMAX,3,
     &            '@AKD,AKNCD,AKDWD [m^2/s]  vs r@',2+INQ)
C
C      DO 400 NR=1,NRMAX-1
C         GYR(NR+1,1) = GLOG(AK  (NR,3),1.D-2,1.D2)
C         GYR(NR+1,2) = GLOG(AKNC(NR,3),1.D-2,1.D2)
C         GYR(NR+1,3) = GLOG(AKDW(NR,3),1.D-2,1.D2)
C         GYR(NR+1,1) = GCLIP(AK(NR,3)))
C         GYR(NR+1,2) = GCLIP(AKNC(NR,3))
C         GYR(NR+1,3) = GCLIP(AKDW(NR,3))
C  400 CONTINUE
C         GYR(1,1) = GCLIP(AK(1,3))
C         GYR(1,2) = GCLIP(AKNC(1,3))
C         GYR(1,3) = GCLIP(AKDW(1,3))
C         GYR(NRMAX+1,1) = GCLIP(AK(NRMAX,3)/2)
C         GYR(NRMAX+1,2) = GCLIP(AKNC(NRMAX,3)/2)
C         GYR(NRMAX+1,3) = GCLIP(AKDW(NRMAX,3)/2)
C      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,
C     &            '@AKT,AKNCT,AKDWT [m^2/s]  vs r@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : VGR
C
C     ***********************************************************
C
      SUBROUTINE TRGRR8(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(VGR2(NR,1))
  100 CONTINUE
      GYR(1,1) = GCLIP(VGR2(1,1))
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,
     &            KGR1,2+INQ)
C
      DO 200 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(VGR2(NR,2))
  200 CONTINUE
      GYR(1,1) = GCLIP(VGR2(1,2))
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,
     &            KGR2,2+INQ)
C
      DO 300 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(VGR3(NR,1))
         GYR(NR+1,2) = GCLIP(VGR3(NR,2))
         GYR(NR+1,3) = GCLIP(VGR3(NR,3))
  300 CONTINUE
      GYR(1,1) = GCLIP(VGR3(1,1))
      GYR(1,2) = GCLIP(VGR3(1,2))
      GYR(1,3) = GCLIP(VGR3(1,3))
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,
     &            KGR3,2+INQ)
C
      DO 400 NR=1,NRMAX
         GYR(NR+1,1) = GCLIP(VGR4(NR,1))
         GYR(NR+1,2) = GCLIP(VGR4(NR,2))
         GYR(NR+1,3) = GCLIP(VGR4(NR,3))
  400 CONTINUE
      GYR(1,1) = GCLIP(VGR4(1,1))
      GYR(1,2) = GCLIP(VGR4(1,2))
      GYR(1,3) = GCLIP(VGR4(1,3))
C
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,
     &            KGR4,2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL PROFILE : AD,AV,AVK,TAUB,TAUF
C
C     ***********************************************************
C
      SUBROUTINE TRGRR9(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NS=1,NSM
      DO 100 NR=1,NRMAX+1
         GYR(NR+1,NS) = GCLIP(AD(NR,NS))
  100 CONTINUE
         GYR(1,NS) = GCLIP(AD(2,NS))
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,NSM,
     &            '@AD [m^2/s]  vs r@',2+INQ)
C
      DO 200 NS=1,NSM
      DO 200 NR=1,NRMAX
         GYR(NR+1,NS) = GCLIP(AV(NR,NS))
  200 CONTINUE
         GYR(1,NS) = 0.0
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,NSM,
     &            '@AV [m/s]  vs r@',2+INQ)
C
      DO 300 NS=1,NSM
      DO 300 NR=1,NRMAX
         GYR(NR+1,NS) = GCLIP(AVK(NR,NS))
  300 CONTINUE
         GYR(1,NS) = 0.0
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,NSM,
     &            '@AVK [m/s]  vs r@',2+INQ)
C
      DO 400 NR=1,NRMAX
         GYR(NR,1) = GCLIP(TAUB(NR))
         GYR(NR,2) = GCLIP(TAUF(NR))
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@TAUB,TAUF [s]  vs r@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           GRAPHIC : RADIAL LISSAGE : S ALPHA
C
C     ***********************************************************
C
      SUBROUTINE TRGRY1(INQ)
C
      INCLUDE 'trcomm.h'
C
      CALL PAGES
C
      DO 100 NR=1,NRMAX
         GYR(NR,1) = GCLIP(VGR1(NR,3))
         GYR(NR,2) = GCLIP(VGR1(NR,2))
  100 CONTINUE
C      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
C     &            '@s  vs alpha@',2+INQ)
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
     &            '@s  vs alpha@',2+4)
C
      DO 200 NR=1,NRMAX
         GYR(NR,1) = GCLIP(AK(NR,1))
         GYR(NR,2) = GCLIP(AK(NR,2))
  200 CONTINUE
      CALL TRGR1D(15.5,24.5,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
     &           '@AKD vs AKE@',2+INQ)
C
      DO 300 NR=1,NRMAX
         GYR(NR,1) = GCLIP(QP(NR))
         GYR(NR,2) = GCLIP(AK(NR,1))
  300 CONTINUE
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
     &            '@AKD vs q@',2+INQ)
C
      DO 400 NR=1,NRMAX-1
         GYR(NR,1) = GCLIP(QP(NR))
         GYR(NR,2) = GCLIP(RT(NR,1))
  400 CONTINUE
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GYR(1,1),GYR(1,2),NRMP,NRMAX-1,1,
     &            '@Te vs q@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
