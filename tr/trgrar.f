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
      INCLUDE 'trcomm.inc'
      CHARACTER K2*1
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(RN(NR,NS))
      ENDDO
      ENDDO
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@NE,ND [10$+20$=/m$+3$=]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
     &               '@NE,ND,NT,NA [10$+20$=/m$+3$=]  vs r@',2+INQ)
      ENDIF
c$$$      IF(MDLEQ0.EQ.0) THEN
      DO NF=1,NFM
      DO NR=1,NRMAX
         GYR(NR,NF) = GUCLIP(RNF(NR,NF))
      ENDDO
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,
     &            '@NB,NF [10$+20$=/m$+3$=]  vs r@',2+INQ)
c$$$      ELSEIF(MDLEQ0.EQ.1) THEN
c$$$      DO NR=1,NRMAX
c$$$         GYR(NR,1) = GUCLIP(RN(NR,7)*1.D5)
c$$$         GYR(NR,2) = GUCLIP(RN(NR,8)*1.D5)
c$$$      ENDDO
c$$$      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
c$$$     &            '@NNC,NNH [10$+15$=/m$+3$=]  vs r@',2+INQ)
c$$$      ENDIF
C
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(RT(NR,NS))
      ENDDO
      ENDDO
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TE,TD [keV]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,
     &               '@TE,TD,TT,TA [keV]  vs r@',2+INQ)
      ENDIF
C
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) =  GUCLIP(RN(NR,NS)*RT(NR,NS)*RKEV*1.D14)
      ENDDO
      ENDDO
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GUCLIP(PNB(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(PNF(NR) * 1.D-6)
      ENDDO
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS+3) = GUCLIP((PRF(NR,NS)+PEX(NR,NS)) * 1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM+3,
     &            '@POH,PNB,PNF,PRF [MW/m$+3$=]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GUCLIP(PRL(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(PCX(NR) * 1.D-6)
         GYR(NR,4) = GUCLIP(PIE(NR) * 1.D-6)
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,
     &            '@POH,PRL,PCX,PIE [MW/m$+3$=]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(BETAL(NR))
         GYR(NR+1,2) = GUCLIP(BETA(NR))
         GYR(NR+1,3) = GUCLIP(BETAQ(NR))
      ENDDO
      GYR(1,1)=GUCLIP(BETA0)
      GYR(1,2)=GUCLIP(BETA0)
      GYR(1,3)=GUCLIP(BETAQ0)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,
     &            '@BETA,<BETA>,BETAQ  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(BETAPL(NR))
         GYR(NR+1,2) = GUCLIP(BETAP(NR))
      ENDDO
      GYR(1,1)=GUCLIP(BETAP0)
      GYR(1,2)=GUCLIP(BETAP0)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &            '@BETAP,<BETAP>  vs r@',2+INQ)
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NF=1,NFM
      DO NR=1,NRMAX
         GYR(NR,NF) = GUCLIP(RTF(NR,NF))
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,
     &            '@TB,TF [keV]  vs r@',2+INQ)
C
      DO NF=1,NFM
      DO NR=1,NRMAX
         GYR(NR,NF) = GUCLIP(1.5D0*RW(NR,NF)*RKEV*1.D14)
      ENDDO
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,
     &            '@WB,WF [MJ]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PBIN(NR)   * 1.D-6)
      ENDDO
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS+1) = GUCLIP(PBCL(NR,NS)   * 1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM+1,
     &            '@PBIN,PBCL [MW/m$+3$=]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PFIN(NR) * 1.D-6)
      ENDDO
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS+1) = GUCLIP(PFCL(NR,NS) * 1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM+1,
     &            '@PFIN,PFCL [MW/m$+3$=]  vs r@',2+INQ)
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(QP(NR))
      ENDDO
      GYR(1,1) = GUCLIP(Q0)
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,
     &            '@QP  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(EZOH(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@EZOH [V/m]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AJ(NR)   * 1.D-6)
         GYR(NR,2) = GUCLIP(AJOH(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(AJNB(NR) * 1.D-6)
         GYR(NR,4) = GUCLIP(AJRF(NR) * 1.D-6)
         GYR(NR,5) = GUCLIP(AJBS(NR) * 1.D-6)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,5,
     &            '@JTOT,JOH,JNB,JRF,JBS [MA/m$+2$=]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GLOG(ETA(NR),1.D-10,1.D0)
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@LOG:ETA  vs r @',11+INQ)
C
C      DO NR=1,NRMAX
C         GYR(NR,1) = GUCLIP(SIE(NR))
C         GYR(NR,2) = GUCLIP(SNB(NR))
C         GYR(NR,3) = GUCLIP(SNF(NR))
C      ENDDO
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ZEFF(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ZEFF  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PZC(NR))
         GYR(NR,2) = GUCLIP(PZFE(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@PZC,PZFE  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ANC(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ANC [10$+20$=/m$+3$=]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ANFE(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@ANFE [10$+20$=/m$+3$=]  vs r@',2+INQ)
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(PIN(NR,NS)*1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
     &            '@PIN [MW/m$+3$=]  vs r@',2+INQ)
C
c$$$      DO NS=1,NSM
c$$$      DO NR=1,NRMAX
c$$$C         GYR(NR,NS) = GUCLIP(SSIN(NR,NS))
c$$$      ENDDO
c$$$      ENDDO
c$$$      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,
c$$$     &            '@SSIN [/sm$+3$=]  vs r@',2+INQ)
C
      NS=1
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(SIE(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@SIE [10$+20$=/sm$+3$=]  vs r@',2+INQ)
C
c$$$      DO NS=1,NSMAX
c$$$      DO NR=1,NRMAX
c$$$         GYR(NR,NS) = GUCLIP(SEX(NR,NS))
c$$$      ENDDO
c$$$      ENDDO
c$$$      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,
c$$$     &            '@SEX [/sm$+3$=]  vs r@',2+INQ)
C
c$$$      DO NS=1,NSM
c$$$      DO NR=1,NRMAX
c$$$         GYR(NR,NS) = GUCLIP(SPE(NR,NS))
c$$$      ENDDO
c$$$      ENDDO
c$$$      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
c$$$     &            '@SPE [/m$+3$=]  vs r@',2+INQ)
C
      NS=1
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(SCX(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@SCX [10$+20$=/sm$+3$=]  vs r@',2+INQ)
C
      NS=1
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(RPSI(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
     &            '@PSI [Wb]  vs r@',2+INQ)
c$$$      CALL TRGRTM
c$$$C
c$$$      CALL MOVE(17.5,4.0)
c$$$      CALL TEXT('PELVEL=',7)
c$$$      CALL NUMBD(PELVEL,'(1PE10.3)',10)
c$$$      CALL TEXT('[m/s]',5)
c$$$      CALL MOVE(17.5,3.0)
c$$$      CALL TEXT('PELRAD=',7)
c$$$      CALL NUMBD(PELRAD,'(1PE10.3)',10)
c$$$      CALL TEXT('[m]',3)
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
      INCLUDE 'trcomm.inc'
      DIMENSION RQFLSUM(NRMP,NSM)
      DIMENSION RNN(NRM,NSM),DTN(NRM,NSM)
      DIMENSION AKNCG(NRM,NSM)
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR1(NR,1))
         GYR(NR+1,2) = GUCLIP(VGR1(NR,2))
         GYR(NR+1,3) = GUCLIP(VGR1(NR,3))
      ENDDO
      GYR(1,1) = GUCLIP((4.D0*VGR1(1,1)-VGR1(2,1))/3.D0)
      GYR(1,2) = 0.0
      GYR(1,3) = 0.0
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,3,
     &           '@G,s,alpha vs r@',2+INQ)
C
      IF(MDNCLS.EQ.0) THEN
C
      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            GYR(NR+1,1) = GLOG(AK(NR,1),1.D-3,1.D2)
            GYR(NR+1,2) = GLOG(AK(NR,2),1.D-3,1.D2)
         ENDDO
         GYR(1,1) = GLOG(AK(1,1),1.D-3,1.D2)
         GYR(1,2) = GLOG(AK(1,2),1.D-3,1.D2)
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,2,
     &               '@LOG:AKE,AKI vs r @',11+INQ)
      ELSE
         DO NR=1,NRMAX
            GYR(NR+1,1) = GUCLIP(AK(NR,3))
            GYR(NR+1,2) = GUCLIP(AKNC(NR,3))
            GYR(NR+1,3) = GUCLIP(AKDW(NR,3))
         ENDDO
         GYR(1,1) = GUCLIP(AK(1,3))
         GYR(1,2) = GUCLIP(AKNC(1,3))
         GYR(1,3) = GUCLIP(AKDW(1,3))
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,3,
     &               '@AKQ,AKNCQ,AKDWQ [m$+2$=/s]  vs r@',2+INQ)
      ENDIF
C
      DO NR=1,NRMAX
C         GYR(NR+1,1) = GLOG(AK  (NR,1),1.D-2,1.D2)
C         GYR(NR+1,2) = GLOG(AKNC(NR,1),1.D-2,1.D2)
C         GYR(NR+1,3) = GLOG(AKDW(NR,1),1.D-2,1.D2)
         GYR(NR+1,1) = GUCLIP(AK(NR,1))
         GYR(NR+1,2) = GUCLIP(AKNC(NR,1))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,1))
      ENDDO
         GYR(1,1) = GUCLIP(AK(1,1))
         GYR(1,2) = GUCLIP(AKNC(1,1))
         GYR(1,3) = GUCLIP(AKDW(1,1))
C         GYR(NRMAX+1,1) = GUCLIP(AK(NRMAX,1)/2)
C         GYR(NRMAX+1,2) = GUCLIP(AKNC(NRMAX,1)/2)
C         GYR(NRMAX+1,3) = GUCLIP(AKDW(NRMAX,1)/2)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,
     &            '@AKE,AKNCE,AKDWE [m$+2$=/s]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
C         GYR(NR+1,1) = GLOG(AK  (NR,2),1.D-2,1.D2)
C         GYR(NR+1,2) = GLOG(AKNC(NR,2),1.D-2,1.D2)
C         GYR(NR+1,3) = GLOG(AKDW(NR,2),1.D-2,1.D2)
         GYR(NR+1,1) = GUCLIP(AK(NR,2))
         GYR(NR+1,2) = GUCLIP(AKNC(NR,2))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,2))
      ENDDO
         GYR(1,1) = GUCLIP(AK(1,2))
         GYR(1,2) = GUCLIP(AKNC(1,2))
         GYR(1,3) = GUCLIP(AKDW(1,2))
C         GYR(NRMAX+1,1) = GUCLIP(AK(NRMAX,2)/2)
C         GYR(NRMAX+1,2) = GUCLIP(AKNC(NRMAX,2)/2)
C         GYR(NRMAX+1,3) = GUCLIP(AKDW(NRMAX,2)/2)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,
     &            '@AKD,AKNCD,AKDWD [m$+2$=/s]  vs r@',2+INQ)
C
      ELSE
C
      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            GYR(NR+1,1) = GLOG(AKLP(NR,1,1),1.D-3,1.D2)
            GYR(NR+1,2) = GLOG(AKLP(NR,2,2),1.D-3,1.D2)
         ENDDO
         GYR(1,1) = GLOG(AKLP(1,1,1),1.D-3,1.D2)
         GYR(1,2) = GLOG(AKLP(1,2,2),1.D-3,1.D2)
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,2,
     &               '@LOG:AKE,AKI vs r @',11+INQ)
      ELSE
         DO NR=1,NRMAX
            GYR(NR+1,1) = GUCLIP(AK(NR,3))
            GYR(NR+1,2) = GUCLIP(AKNC(NR,3))
            GYR(NR+1,3) = GUCLIP(AKDW(NR,3))
         ENDDO
         GYR(1,1) = GUCLIP(AK(1,3))
         GYR(1,2) = GUCLIP(AKNC(1,3))
         GYR(1,3) = GUCLIP(AKDW(1,3))
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,3,
     &               '@AKQ,AKNCQ,AKDWQ [m$+2$=/s]  vs r@',2+INQ)
      ENDIF
C
      DO NR=1,NRMAX
         DO NS=1,NSMAX
            RQFLSUM(NR,NS)=0.D0
            DO NA=1,5
               RQFLSUM(NR,NS)=RQFLSUM(NR,NS)+RQFLS(NR,NA,NS)
            ENDDO
         ENDDO
      ENDDO
      DO NR=1,NRMAX-1
         DO NS=1,NSMAX
            RNN(NR,NS)=(RN(NR+1,NS)+RN(NR,NS))*0.5D0
            DTN(NR,NS)=(RT(NR+1,NS)-RT(NR,NS))*RKEV*RJCB(NR)/DR
         ENDDO
      ENDDO
      NR=NRMAX
      DO NS=1,NSMAX
         RNN(NR,NS)=PNSS(NS)
         DTN(NR,NS)=2.D0*(PTS (NS)-RT(NR,NS))*RKEV*RJCB(NR)/DR
      ENDDO
      DO NR=1,NRMAX
         DO NS=1,NSMAX
            AKNCG(NR,NS)=-RQFLSUM(NR,NS)/(RNN(NR,NS)*DTN(NR,NS))
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AKDW(NR,1)+AKNCG(NR,1))
         GYR(NR+1,2) = GUCLIP(AKNCG(NR,1))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,1))
         GYR(NR+1,4) = GUCLIP(AKLP(NR,1,1))
      ENDDO
         GYR(1,1) = GUCLIP(AKDW(1,1)+AKNCG(1,1))
         GYR(1,2) = GUCLIP(AKNCG(1,1))
         GYR(1,3) = GUCLIP(AKDW(1,1))
         GYR(1,4) = GUCLIP(AKLP(1,1,1))
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,4,
     &            '@AKE,AKNCE,AKDWE,AKEDIAG [m$+2$=/s]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AKDW(NR,2)+AKNCG(NR,2))
         GYR(NR+1,2) = GUCLIP(AKNCG(NR,2))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,2))
         GYR(NR+1,4) = GUCLIP(AKLP(NR,2,2))
      ENDDO
         GYR(1,1) = GUCLIP(AKDW(1,2)+AKNCG(1,2))
         GYR(1,2) = GUCLIP(AKNCG(1,2))
         GYR(1,3) = GUCLIP(AKDW(1,2))
         GYR(1,4) = GUCLIP(AKLP(1,2,2))
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,4,
     &            '@AKD,AKNCD,AKDWD,AKDDIAG [m$+2$=/s]  vs r@',2+INQ)
C
      ENDIF
C
C      DO NR=1,NRMAX-1
C         GYR(NR+1,1) = GLOG(AK  (NR,3),1.D-2,1.D2)
C         GYR(NR+1,2) = GLOG(AKNC(NR,3),1.D-2,1.D2)
C         GYR(NR+1,3) = GLOG(AKDW(NR,3),1.D-2,1.D2)
C         GYR(NR+1,1) = GUCLIP(AK(NR,3)))
C         GYR(NR+1,2) = GUCLIP(AKNC(NR,3))
C         GYR(NR+1,3) = GUCLIP(AKDW(NR,3))
C      ENDDO
C         GYR(1,1) = GUCLIP(AK(1,3))
C         GYR(1,2) = GUCLIP(AKNC(1,3))
C         GYR(1,3) = GUCLIP(AKDW(1,3))
C         GYR(NRMAX+1,1) = GUCLIP(AK(NRMAX,3)/2)
C         GYR(NRMAX+1,2) = GUCLIP(AKNC(NRMAX,3)/2)
C         GYR(NRMAX+1,3) = GUCLIP(AKDW(NRMAX,3)/2)
C      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,
C     &            '@AKT,AKNCT,AKDWT [m$+2$=/s]  vs r@',2+INQ)
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR2(NR,1))
      ENDDO
      GYR(1,1) = GUCLIP(VGR2(1,1))
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,
     &            KGR1,2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR2(NR,2))
      ENDDO
      GYR(1,1) = GUCLIP(VGR2(1,2))
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,
     &            KGR2,2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR3(NR,1))
         GYR(NR+1,2) = GUCLIP(VGR3(NR,2))
         GYR(NR+1,3) = GUCLIP(VGR3(NR,3))
      ENDDO
      GYR(1,1) = GUCLIP(VGR3(1,1))
      GYR(1,2) = GUCLIP(VGR3(1,2))
      GYR(1,3) = GUCLIP(VGR3(1,3))
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,
     &            KGR3,2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR4(NR,1))
         GYR(NR+1,2) = GUCLIP(VGR4(NR,2))
         GYR(NR+1,3) = GUCLIP(VGR4(NR,3))
      ENDDO
      GYR(1,1) = GUCLIP(VGR4(1,1))
      GYR(1,2) = GUCLIP(VGR4(1,2))
      GYR(1,3) = GUCLIP(VGR4(1,3))
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
      INCLUDE 'trcomm.inc'
      DIMENSION RGFLSUM(NRMP,NSM)
      DIMENSION DNN(NRM,NSM)
      DIMENSION ADNCG(NRMP,NSM)
C
      CALL PAGES
C
      IF(MDNCLS.EQ.0) THEN
         DO NS=1,NSMAX
         DO NR=1,NRMAX
            GYR(NR+1,NS) = GUCLIP(AD(NR,NS))
         ENDDO
            GYR(1,NS) = GUCLIP(AD(2,NS))
         ENDDO
      ELSE
         MODE=0
         IF(MODE.EQ.0) THEN
            DO NR=1,NRMAX
               DO NS=1,NSMAX
                  ADNCG(NR,NS)=ADNCS(NR,NS)
               ENDDO
            ENDDO
         ELSE
            DO NR=1,NRMAX
               DO NS=1,NSMAX
                  RGFLSUM(NR,NS)=0.D0
                  DO NA=1,5
                     RGFLSUM(NR,NS)=RGFLSUM(NR,NS)+RGFLS(NR,NA,NS)
                  ENDDO
               ENDDO
            ENDDO
            DO NR=1,NRMAX-1
               DO NS=1,NSMAX
                  DNN(NR,NS)=(RN(NR+1,NS)-RN(NR,NS))*RJCB(NR)/DR
               ENDDO
            ENDDO
            NR=NRMAX
            DO NS=1,NSMAX
               DNN(NR,NS)=2.D0*(PNSS(NS)-RN(NR,NS))*RJCB(NR)/DR
            ENDDO
            DO NR=1,NRMAX
               DO NS=1,NSMAX
                  ADNCG(NR,NS)=-RGFLSUM(NR,NS)/DNN(NR,NS)
               ENDDO
            ENDDO
         ENDIF
C     
         DO NS=1,NSMAX
         DO NR=1,NRMAX
            GYR(NR+1,NS) = GUCLIP(ADNCG(NR,NS)+ADDW(NR,NS))
         ENDDO
            GYR(1,NS) = GUCLIP(ADNCG(2,NS)+ADDW(2,NS))
         ENDDO
      ENDIF
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,NSMAX,
     &            '@AD [m$+2$=/s]  vs r@',2+INQ)
C
      IF(MDNCLS.EQ.0) THEN
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               GYR(NR+1,NS) = GUCLIP(AV(NR,NS))
            ENDDO
            GYR(1,NS) = 0.0
         ENDDO
      ELSE
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               GYR(NR+1,NS) = GUCLIP(AV(NR,NS)+AVNCS(NR,NS))
            ENDDO
            GYR(1,NS) = 0.0
         ENDDO
      ENDIF
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,NSMAX,
     &            '@AV [m/s]  vs r@',2+INQ)
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         GYR(NR+1,NS) = GUCLIP(AVK(NR,NS))
      ENDDO
         GYR(1,NS) = 0.0
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,NSMAX,
     &            '@AVK [m/s]  vs r@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(TAUB(NR))
         GYR(NR,2) = GUCLIP(TAUF(NR))
      ENDDO
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
      INCLUDE 'trcomm.inc'
C
      CALL PAGES
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(VGR1(NR,3))
         GYR(NR,2) = GUCLIP(VGR1(NR,2))
      ENDDO
C      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
C     &            '@s  vs alpha@',2+INQ)
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
     &            '@s  vs alpha@',2+4)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AK(NR,1))
         GYR(NR,2) = GUCLIP(AK(NR,2))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
     &           '@AKD vs AKE@',2+INQ)
C
      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(QP(NR))
         GYR(NR,2) = GUCLIP(AK(NR,1))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
     &            '@AKD vs q@',2+INQ)
C
      DO NR=1,NRMAX-1
         GYR(NR,1) = GUCLIP(QP(NR))
         GYR(NR,2) = GUCLIP(RT(NR,1))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GYR(1,1),GYR(1,2),NRMP,NRMAX-1,1,
     &            '@Te vs q@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
