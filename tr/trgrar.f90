!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRR0(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN) :: K2
      INTEGER(4),       INTENT(IN) :: INQ

      IF(RHOA.NE.1.D0) NRMAX=NROMAX

      IF(K2.EQ.'1') CALL TRGRR1(INQ)
      IF(K2.EQ.'2') CALL TRGRR2(INQ)
      IF(K2.EQ.'3') CALL TRGRR3(INQ)
      IF(K2.EQ.'4') CALL TRGRR4(INQ)
      IF(K2.EQ.'5') CALL TRGRR5(INQ)
      IF(K2.EQ.'6') CALL TRGRR6(INQ)
      IF(K2.EQ.'7') CALL TRGRR7(INQ)
      IF(K2.EQ.'8') CALL TRGRR8(INQ)
      IF(K2.EQ.'9') CALL TRGRR9(INQ)
      IF(K2.EQ.'A') CALL TRGRRA(INQ)

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRGRR0
!
!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRY0(K2,INQ)

      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN) :: K2
      INTEGER(4),       INTENT(IN) :: INQ

      IF(K2.EQ.'1') CALL TRGRY1(INQ)
      RETURN
      END SUBROUTINE TRGRY0

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : RNS,RNF,RTS,RTF

!     ***********************************************************

      SUBROUTINE TRGRR1(INQ)

      USE TRCOMM, ONLY : GRM, GYR, MDLEQ0, MDLNF, NFM, NRMAX, NRMP, NSM, RKEV, RN, RNF, RT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NF, NR, NS
      REAL(4)    :: GUCLIP

      CALL PAGES

      IF(MDLNF.EQ.0) THEN
         DO NS=1,2
            DO NR=1,NRMAX
               GYR(NR,NS) = GUCLIP(RN(NR,NS))
            ENDDO
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@NE,ND [10$+20$=/m$+3$=]  vs r@',2+INQ)
      ELSE
         DO NS=1,NSM
            DO NR=1,NRMAX
               GYR(NR,NS) = GUCLIP(RN(NR,NS))
            ENDDO
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM, '@NE,ND,NT,NA [10$+20$=/m$+3$=]  vs r@',2+INQ)
      ENDIF

      IF(MDLEQ0.EQ.0) THEN
         DO NF=1,NFM
            DO NR=1,NRMAX
               GYR(NR,NF) = GUCLIP(RNF(NR,NF))
            ENDDO
         ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM, '@NB,NF [10$+20$=/m$+3$=]  vs r@',2+INQ)
      ELSEIF(MDLEQ0.EQ.1) THEN
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RN(NR,7)*1.D5)
            GYR(NR,2) = GUCLIP(RN(NR,8)*1.D5)
         ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2, &
     &        '@NNC [10$+15$=/m$+3$=], NNH [10$+15$=/m$+3$=]  vs r@',2+INQ)
      ENDIF

      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(RT(NR,NS))
      ENDDO
      ENDDO
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@TE,TD [keV]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,'@TE,TD,TT,TA [keV]  vs r@',2+INQ)
      ENDIF

      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) =  GUCLIP(RN(NR,NS)*RT(NR,NS)*RKEV*1.D14)
      ENDDO
      ENDDO
      IF(MDLNF.EQ.0) THEN
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@PE,PD [MPa]  vs r@',2+INQ)
      ELSE
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM, '@PE,PD,PT,PA [MPa]  vs r@',2+INQ)
      ENDIF
      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRR1

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : WS,WF,BETA,BETAP

!     ***********************************************************

      SUBROUTINE TRGRR2(INQ)

      USE TRCOMM, ONLY : &
           AK, AKNC,AKDW, GRG, GRM, GYR, NRMAX, NRMP, NSM, PNB, PNF, POH, &
           PRF, PRB, PRC, PRL, PRSUM, PCX, PIE, PEX 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NS, NR
      REAL(4)    :: GUCLIP


      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GUCLIP(PNB(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(PNF(NR) * 1.D-6)
         GYR(NR,4) = GUCLIP(PRSUM(NR) * 1.D-6)
      ENDDO
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS+4) = GUCLIP((PRF(NR,NS)+PEX(NR,NS)) * 1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM+4, &
     &            '@POH,PNB,PNF,PR,PRF [MW/m$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GUCLIP(PRB(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(PRC(NR) * 1.D-6)
         GYR(NR,4) = GUCLIP(PRL(NR) * 1.D-6)
         GYR(NR,5) = GUCLIP(PCX(NR) * 1.D-6)
         GYR(NR,6) = GUCLIP(PIE(NR) * 1.D-6)
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,6, &
     &            '@POH,PRB,PRC,PRL,PCX,PIE [MW/m$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AK(NR,1))
         GYR(NR+1,2) = GUCLIP(AKNC(NR,1))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,1))
      ENDDO
         GYR(1,1) = GUCLIP(AK(1,1))
         GYR(1,2) = GUCLIP(AKNC(1,1))
         GYR(1,3) = GUCLIP(AKDW(1,1))
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3, &
      &           '@AKE,AKNCE,AKDWE [m$+2$=/s]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(AK(NR,2))
         GYR(NR+1,2) = GUCLIP(AKNC(NR,2))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,2))
      ENDDO
         GYR(1,1) = GUCLIP(AK(1,2))
         GYR(1,2) = GUCLIP(AKNC(1,2))
         GYR(1,3) = GUCLIP(AKDW(1,2))
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3, &
      &           '@AKD,AKNCD,AKDWD [m$+2$=/s]  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRR2

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : POWER

!     ***********************************************************

      SUBROUTINE TRGRR3(INQ)

      USE TRCOMM, ONLY : GRM, GYR, NFM, NRMAX, NRMP, NSM, PBCL, PBIN, PFCL, PFIN, RKEV, RTF, RW
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NF, NR, NS
      REAL(4)    :: GUCLIP

      CALL PAGES

      DO NF=1,NFM
      DO NR=1,NRMAX
         GYR(NR,NF) = GUCLIP(RTF(NR,NF))
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM,'@TB,TF [keV]  vs r@',2+INQ)

      DO NF=1,NFM
      DO NR=1,NRMAX
         GYR(NR,NF) = GUCLIP(1.5D0*RW(NR,NF)*RKEV*1.D14)
      ENDDO
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NFM, '@WB,WF [MJ]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PBIN(NR)   * 1.D-6)
      ENDDO
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS+1) = GUCLIP(PBCL(NR,NS)   * 1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM+1,'@PBIN,PBCL [MW/m$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PFIN(NR) * 1.D-6)
      ENDDO
      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS+1) = GUCLIP(PFCL(NR,NS) * 1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM+1,'@PFIN,PFCL [MW/m$+3$=]  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR3

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : Q,EZ,AJ,S

!     ***********************************************************

      SUBROUTINE TRGRR4(INQ)

      USE TRCOMM, ONLY : AJ, AJBS, AJNB, AJOH, AJRF, ETA, EZOH, GRG, GRM, GYR, NRMAX, NRMP, Q0, QP
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR
      REAL(4)    :: GUCLIP, GLOG

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(QP(NR))
      ENDDO
      GYR(1,1) = GUCLIP(Q0)

      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,'@QP  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(EZOH(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,1, '@EZOH [V/m]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AJ(NR)   * 1.D-6)
         GYR(NR,2) = GUCLIP(AJOH(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(AJNB(NR) * 1.D-6)
         GYR(NR,4) = GUCLIP(AJRF(NR) * 1.D-6)
         GYR(NR,5) = GUCLIP(AJBS(NR) * 1.D-6)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,5, &
     &            '@JTOT,JOH,JNB,JRF,JBS [MA/m$+2$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GLOG(ETA(NR),1.D-10,1.D0)
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1, '@LOG:ETA  vs r @',11+INQ)

!      DO NR=1,NRMAX
!         GYR(NR,1) = GUCLIP(SIE(NR))
!         GYR(NR,2) = GUCLIP(SNB(NR))
!         GYR(NR,3) = GUCLIP(SNF(NR))
!      ENDDO
!      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,3,
!     &            '@SIE,SNB,SNF [/sm^3]  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR4

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : ZEFF,IMPURITY

!     ***********************************************************

      SUBROUTINE TRGRR5(INQ)

      USE TRCOMM, ONLY : GRM, GYR, NRMAX, NRMP, PZC, PZFE, VPOL, VTOR, ZEFF
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR
      REAL(4)    :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(ZEFF(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,1,'@ZEFF  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(PZC(NR))
         GYR(NR,2) = GUCLIP(PZFE(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@PZC,PZFE  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(VTOR(NR))
!CC         GYR(NR,2) = GUCLIP(VPAR(NR))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,'@V$-tor$=, V$-para$= [m/s]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(VPOL(NR))
!CC         GYR(NR,2) = GUCLIP(VPRP(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,'@V$-pol$=, V$-perp$= [m/s]  vs r@',2+INQ)

!$$$      DO NR=1,NRMAX
!$$$         GYR(NR,1) = GUCLIP(ANC(NR))
!$$$      ENDDO
!$$$      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
!$$$     &            '@ANC [10$+20$=/m$+3$=]  vs r@',2+INQ)
!$$$C
!$$$      DO NR=1,NRMAX
!$$$         GYR(NR,1) = GUCLIP(ANFE(NR))
!$$$      ENDDO
!$$$      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,
!$$$     &            '@ANFE [10$+20$=/m$+3$=]  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR5

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : PIN,SSIN,PELLET

!     ***********************************************************

      SUBROUTINE TRGRR6(INQ)

      USE TRCOMM, ONLY : GRM, GYR, MDLUF, NRMAX, NRMP, NSM, PIN, RPSI, SCX, SEX, SIE, SNBU, SSIN, SWLU
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR, NS
      REAL(4)    :: GUCLIP

      CALL PAGES

      DO NS=1,NSM
      DO NR=1,NRMAX
         GYR(NR,NS) = GUCLIP(PIN(NR,NS)*1.D-6)
      ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,'@PIN [MW/m$+3$=]  vs r@',2+INQ)

      IF(MDLUF.NE.0) THEN
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(SEX(NR,1))
            GYR(NR,2) = GUCLIP(SEX(NR,2))
            GYR(NR,3) = GUCLIP(SNBU(1,NR,1))
            GYR(NR,4) = GUCLIP(SNBU(1,NR,2))
            GYR(NR,5) = GUCLIP(SWLU(1,NR))
         ENDDO
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,5,'@SEX [10$+20$=/sm$+3$=]  vs r@',2+INQ)
      ELSE
         DO NS=1,NSM
            DO NR=1,NRMAX
               GYR(NR,NS) = GUCLIP(SSIN(NR,NS))
            ENDDO
         ENDDO
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,NSM,'@SSIN [/sm$+3$=]  vs r@',2+INQ)
      ENDIF

!$$$      DO NS=1,NSM
!$$$      DO NR=1,NRMAX
!$$$         GYR(NR,NS) = GUCLIP(SPE(NR,NS))
!$$$      ENDDO
!$$$      ENDDO
!$$$      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM,
!$$$     &            '@SPE [/m$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(SCX(NR))
         GYR(NR,2) = GUCLIP(SIE(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@SCX, SIE [10$+20$=/sm$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(RPSI(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,1,'@PSI [Wb]  vs r@',2+INQ)
!$$$      CALL TRGRTM
!$$$C
!$$$      CALL MOVE(17.5,4.0)
!$$$      CALL TEXT('PELVEL=',7)
!$$$      CALL NUMBD(PELVEL,'(1PE10.3)',10)
!$$$      CALL TEXT('[m/s]',5)
!$$$      CALL MOVE(17.5,3.0)
!$$$      CALL TEXT('PELRAD=',7)
!$$$      CALL NUMBD(PELRAD,'(1PE10.3)',10)
!$$$      CALL TEXT('[m]',3)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR6

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : ETA,KAI

!     ***********************************************************

      SUBROUTINE TRGRR7(INQ)

      USE TRCOMM, ONLY : AK, AKDW, AKLP, AKNC, DR, GRG, GYR, MDLUF, MDNCLS, NRMAX, NRMP, NSM, NSMAX, PNSS, &
     &                   PTS, RHOG, RHOM, RJCB, RKEV, RN, RQFLS, RT, VGR1
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR, NS
      REAL(8), DIMENSION(NRMP,NSMAX) :: RQFLSUM
      REAL(8), DIMENSION(NRMAX,NSMAX)  :: AKNCG, DTN, RNN
      REAL(4)    :: GUCLIP, GLOG
      REAL(8)  :: DERIV3P

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR1(NR,1))
         GYR(NR+1,2) = GUCLIP(VGR1(NR,2))
         GYR(NR+1,3) = GUCLIP(VGR1(NR,3))
      ENDDO
      GYR(1,1) = GUCLIP((4.D0*VGR1(1,1)-VGR1(2,1))/3.D0)
      GYR(1,2) = 0.0
      GYR(1,3) = 0.0
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,3,'@G,s,alpha vs r@',2+INQ)

      IF(MDNCLS.EQ.0) THEN

      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            GYR(NR+1,1) = GLOG(AK(NR,1),1.D-3,1.D2)
            GYR(NR+1,2) = GLOG(AK(NR,2),1.D-3,1.D2)
         ENDDO
         GYR(1,1) = GLOG(AK(1,1),1.D-3,1.D2)
         GYR(1,2) = GLOG(AK(1,2),1.D-3,1.D2)
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,2,'@LOG:AKE,AKI vs r @',11+INQ)
      ELSE
         DO NR=1,NRMAX
            GYR(NR+1,1) = GUCLIP(AK(NR,3))
            GYR(NR+1,2) = GUCLIP(AKNC(NR,3))
            GYR(NR+1,3) = GUCLIP(AKDW(NR,3))
         ENDDO
         GYR(1,1) = GUCLIP(AK(1,3))
         GYR(1,2) = GUCLIP(AKNC(1,3))
         GYR(1,3) = GUCLIP(AKDW(1,3))
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,3,'@AKQ,AKNCQ,AKDWQ [m$+2$=/s]  vs r@',2+INQ)
      ENDIF

      DO NR=1,NRMAX
!         GYR(NR+1,1) = GLOG(AK  (NR,1),1.D-2,1.D2)
!         GYR(NR+1,2) = GLOG(AKNC(NR,1),1.D-2,1.D2)
!         GYR(NR+1,3) = GLOG(AKDW(NR,1),1.D-2,1.D2)
         GYR(NR+1,1) = GUCLIP(AK(NR,1))
         GYR(NR+1,2) = GUCLIP(AKNC(NR,1))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,1))
      ENDDO
         GYR(1,1) = GUCLIP(AK(1,1))
         GYR(1,2) = GUCLIP(AKNC(1,1))
         GYR(1,3) = GUCLIP(AKDW(1,1))
!         GYR(NRMAX+1,1) = GUCLIP(AK(NRMAX,1)/2)
!         GYR(NRMAX+1,2) = GUCLIP(AKNC(NRMAX,1)/2)
!         GYR(NRMAX+1,3) = GUCLIP(AKDW(NRMAX,1)/2)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,'@AKE,AKNCE,AKDWE [m$+2$=/s]  vs r@',2+INQ)

      DO NR=1,NRMAX
!         GYR(NR+1,1) = GLOG(AK  (NR,2),1.D-2,1.D2)
!         GYR(NR+1,2) = GLOG(AKNC(NR,2),1.D-2,1.D2)
!         GYR(NR+1,3) = GLOG(AKDW(NR,2),1.D-2,1.D2)
         GYR(NR+1,1) = GUCLIP(AK(NR,2))
         GYR(NR+1,2) = GUCLIP(AKNC(NR,2))
         GYR(NR+1,3) = GUCLIP(AKDW(NR,2))
      ENDDO
         GYR(1,1) = GUCLIP(AK(1,2))
         GYR(1,2) = GUCLIP(AKNC(1,2))
         GYR(1,3) = GUCLIP(AKDW(1,2))
!         GYR(NRMAX+1,1) = GUCLIP(AK(NRMAX,2)/2)
!         GYR(NRMAX+1,2) = GUCLIP(AKNC(NRMAX,2)/2)
!         GYR(NRMAX+1,3) = GUCLIP(AKDW(NRMAX,2)/2)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,'@AKD,AKNCD,AKDWD [m$+2$=/s]  vs r@',2+INQ)

      ELSE

      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            GYR(NR+1,1) = GLOG(AKLP(NR,1,1),1.D-3,1.D2)
            GYR(NR+1,2) = GLOG(AKLP(NR,2,2),1.D-3,1.D2)
         ENDDO
         GYR(1,1) = GLOG(AKLP(1,1,1),1.D-3,1.D2)
         GYR(1,2) = GLOG(AKLP(1,2,2),1.D-3,1.D2)
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,2,'@LOG:AKE,AKI vs r @',11+INQ)
      ELSE
         DO NR=1,NRMAX
            GYR(NR+1,1) = GUCLIP(AK(NR,3))
            GYR(NR+1,2) = GUCLIP(AKNC(NR,3))
            GYR(NR+1,3) = GUCLIP(AKDW(NR,3))
         ENDDO
         GYR(1,1) = GUCLIP(AK(1,3))
         GYR(1,2) = GUCLIP(AKNC(1,3))
         GYR(1,3) = GUCLIP(AKDW(1,3))
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX,3,'@AKQ,AKNCQ,AKDWQ [m$+2$=/s]  vs r@',2+INQ)
      ENDIF

      RQFLSUM(1:NRMAX,1:NSMAX)=0.D0
      DO NR=1,NRMAX
         DO NS=1,NSMAX
            RQFLSUM(NR,NS)=RQFLSUM(NR,NS)+SUM(RQFLS(NR,1:5,NS))
         ENDDO
      ENDDO
      DO NS=1,NSMAX
         RNN(1:NRMAX-1,NS)=(RN(2:NRMAX,NS)+RN(1:NRMAX-1,NS))*0.5D0
         DTN(1:NRMAX-1,NS)=(RT(2:NRMAX,NS)-RT(1:NRMAX-1,NS))*RKEV*RJCB(1:NRMAX-1)/DR
      ENDDO
      NR=NRMAX
      DO NS=1,NSMAX
         RNN(NR,NS)=PNSS(NS)
         DTN(NR,NS)=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS),RHOG(NR),RHOM(NR),RHOM(NR-1))*RKEV
      ENDDO
      AKNCG(1:NRMAX,1:NSMAX)=-RQFLSUM(1:NRMAX,1:NSMAX)/(RNN(1:NRMAX,1:NSMAX)*DTN(1:NRMAX,1:NSMAX))

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
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,4,'@AKE,AKNCE,AKDWE,AKEDIAG [m$+2$=/s]  vs r@',2+INQ)

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
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,4,'@AKD,AKNCD,AKDWD,AKDDIAG [m$+2$=/s]  vs r@',2+INQ)

      ENDIF

!      DO NR=1,NRMAX-1
!         GYR(NR+1,1) = GLOG(AK  (NR,3),1.D-2,1.D2)
!         GYR(NR+1,2) = GLOG(AKNC(NR,3),1.D-2,1.D2)
!         GYR(NR+1,3) = GLOG(AKDW(NR,3),1.D-2,1.D2)
!         GYR(NR+1,1) = GUCLIP(AK(NR,3)))
!         GYR(NR+1,2) = GUCLIP(AKNC(NR,3))
!         GYR(NR+1,3) = GUCLIP(AKDW(NR,3))
!      ENDDO
!         GYR(1,1) = GUCLIP(AK(1,3))
!         GYR(1,2) = GUCLIP(AKNC(1,3))
!         GYR(1,3) = GUCLIP(AKDW(1,3))
!         GYR(NRMAX+1,1) = GUCLIP(AK(NRMAX,3)/2)
!         GYR(NRMAX+1,2) = GUCLIP(AKNC(NRMAX,3)/2)
!         GYR(NRMAX+1,3) = GUCLIP(AKDW(NRMAX,3)/2)
!      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX,3,
!     &            '@AKT,AKNCT,AKDWT [m$+2$=/s]  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR7

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : VGR

!     ***********************************************************

      SUBROUTINE TRGRR8(INQ)

      USE TRCOMM, ONLY : GRG, GYR, KGR1, KGR2, KGR3, KGR4, NRMAX, NRMP, &
           VGR1, VGR2, VGR3, VGR4
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR
      REAL(4)    :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR1(NR,1))
      ENDDO
      GYR(1,1) = GUCLIP(VGR1(1,1))
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,KGR1,2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR2(NR,1))
      ENDDO
      GYR(1,1) = GUCLIP(VGR2(1,1))
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,1,KGR2,2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR3(NR,1))
         GYR(NR+1,2) = GUCLIP(VGR3(NR,2))
         GYR(NR+1,3) = GUCLIP(VGR3(NR,3))
      ENDDO
      GYR(1,1) = GUCLIP(VGR3(1,1))
      GYR(1,2) = GUCLIP(VGR3(1,2))
      GYR(1,3) = GUCLIP(VGR3(1,3))
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,KGR3,2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(VGR4(NR,1))
         GYR(NR+1,2) = GUCLIP(VGR4(NR,2))
         GYR(NR+1,3) = GUCLIP(VGR4(NR,3))
      ENDDO
!$$$      GYR(1,1) = GUCLIP(VGR4(1,1))
!$$$      GYR(1,2) = GUCLIP(VGR4(1,2))
!$$$      GYR(1,3) = GUCLIP(VGR4(1,3))
      GYR(1,1) = 0.0
      GYR(1,2) = 0.0
      GYR(1,3) = 0.0

      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,KGR4,2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR8

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : AD,AV,AVK,TAUB,TAUF

!     ***********************************************************

      SUBROUTINE TRGRR9(INQ)

      USE TRCOMM, ONLY : AD, ADDW, ADNCG, AV, AVK, AVNCG, CDP, CNP, DR, GRG, GRM, GYR, MDNCLS, NRMAX, NRMP, NSM, &
     &                   NSMAX, PNSS, RGFLS, RHOG, RHOM, RJCB, RN, TAUB, TAUF
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: MODE, NR, NS
      REAL(8), DIMENSION(NRMP,NSMAX):: RGFLSUM
      REAL(8), DIMENSION(NRMAX,NSMAX) :: DNN
      REAL(4)  :: GUCLIP
      REAL(8)  :: DERIV3P

      CALL PAGES

      IF(MDNCLS.EQ.0) THEN
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               GYR(NR+1,NS) = GUCLIP(AD(NR,NS))
            ENDDO
            GYR(1,NS) = GUCLIP(AD(2,NS))
         ENDDO
      ELSE
         MODE=0
         IF(MODE.NE.0) THEN
            RGFLSUM(1:NRMAX,1:NSMAX)=0.D0
            DO NR=1,NRMAX
               DO NS=1,NSMAX
                  RGFLSUM(NR,NS)=RGFLSUM(NR,NS)+SUM(RGFLS(NR,1:5,NS))
               ENDDO
            ENDDO
            DO NS=1,NSMAX
               DNN(1:NRMAX-1,NS)=(RN(1:NRMAX-1+1,NS)-RN(1:NRMAX-1,NS))*RJCB(1:NRMAX-1)/DR
            ENDDO
            NR=NRMAX
            DO NS=1,NSMAX
               DNN(NR,NS)=DERIV3P(PNSS(NS),RN(NR,NS),RN(NR-1,NS),RHOG(NR),RHOM(NR),RHOM(NR-1))
            ENDDO
            ADNCG(1:NRMAX,1:NSMAX)=-RGFLSUM(1:NRMAX,1:NSMAX)/DNN(1:NRMAX,1:NSMAX)
         ENDIF

         DO NS=1,NSMAX
            DO NR=1,NRMAX
               GYR(NR+1,NS) = GUCLIP(CNP*ADNCG(NR,NS)+CDP*ADDW(NR,NS))
            ENDDO
            GYR(1,NS) = GUCLIP(CNP*ADNCG(2,NS)+CDP*ADDW(2,NS))
         ENDDO
      ENDIF
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,NSMAX,'@AD [m$+2$=/s]  vs r@',2+INQ)

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
               GYR(NR+1,NS) = GUCLIP(AV(NR,NS)+AVNCG(NR,NS))
            ENDDO
            GYR(1,NS) = 0.0
         ENDDO
      ENDIF
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRG,GYR,NRMP,NRMAX+1,NSMAX,'@AV [m/s]  vs r@',2+INQ)

      DO NS=1,NSMAX
      DO NR=1,NRMAX
         GYR(NR+1,NS) = GUCLIP(AVK(NR,NS))
      ENDDO
         GYR(1,NS) = 0.0
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,NSMAX,'@AVK [m/s]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(TAUB(NR))
         GYR(NR,2) = GUCLIP(TAUF(NR))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@TAUB,TAUF [s]  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRR9

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : WS,WF,BETA,BETAP

!     ***********************************************************

      SUBROUTINE TRGRRA(INQ)

      USE TRCOMM, ONLY : BETA, BETA0, BETAL, BETAP, BETAP0, BETAPL, BETAQ, BETAQ0, GRG, GRM, GYR, NRMAX, NRMP, NSM, PCX, &
     &                   PEX, PIE, PNB, PNF, POH, PRF, PRL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NS, NR
      REAL(4)    :: GUCLIP


      CALL PAGES

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
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,NSM+3, &
     &            '@POH,PNB,PNF,PRF [MW/m$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(POH(NR) * 1.D-6)
         GYR(NR,2) = GUCLIP(PRL(NR) * 1.D-6)
         GYR(NR,3) = GUCLIP(PCX(NR) * 1.D-6)
         GYR(NR,4) = GUCLIP(PIE(NR) * 1.D-6)
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,4, &
     &            '@POH,PRL,PCX,PIE [MW/m$+3$=]  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(BETAL(NR))
         GYR(NR+1,2) = GUCLIP(BETA(NR))
         GYR(NR+1,3) = GUCLIP(BETAQ(NR))
      ENDDO
      GYR(1,1)=GUCLIP(BETA0)
      GYR(1,2)=GUCLIP(BETA0)
      GYR(1,3)=GUCLIP(BETAQ0)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,3,'@BETA,<BETA>,BETAQ  vs r@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR+1,1) = GUCLIP(BETAPL(NR))
         GYR(NR+1,2) = GUCLIP(BETAP(NR))
      ENDDO
      GYR(1,1)=GUCLIP(BETAP0)
      GYR(1,2)=GUCLIP(BETAP0)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,'@BETAP,<BETAP>  vs r@',2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRRA

!     ***********************************************************

!           GRAPHIC : RADIAL LISSAGE : S ALPHA

!     ***********************************************************

      SUBROUTINE TRGRY1(INQ)

      USE TRCOMM, ONLY : AK, GYR, NRMAX, NRMP, QP, RT, VGR1
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR
      REAL(4)    :: GUCLIP

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(VGR1(NR,3))
         GYR(NR,2) = GUCLIP(VGR1(NR,2))
      ENDDO
!      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,
!     &            '@s  vs alpha@',2+INQ)
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,'@s  vs alpha@',2+4)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(AK(NR,1))
         GYR(NR,2) = GUCLIP(AK(NR,2))
      ENDDO
      CALL TRGR1D(15.5,24.5,11.0,17.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,'@AKD vs AKE@',2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1) = GUCLIP(QP(NR))
         GYR(NR,2) = GUCLIP(AK(NR,1))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GYR(1,1),GYR(1,2),NRMP,NRMAX,1,'@AKD vs q@',2+INQ)

      DO NR=1,NRMAX-1
         GYR(NR,1) = GUCLIP(QP(NR))
         GYR(NR,2) = GUCLIP(RT(NR,1))
      ENDDO
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GYR(1,1),GYR(1,2),NRMP,NRMAX-1,1,'@Te vs q@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRY1

!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRN0(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA
      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(IN) :: K2
      INTEGER(4),      INTENT(IN) :: INQ

      IF(RHOA.NE.1.D0) NRMAX=NROMAX

      IF(K2.EQ.'1') THEN
         CALL TRGRN1(INQ)
         CALL TRGRN2(INQ)
         CALL TRGRN3(INQ)
      ELSE IF(K2.EQ.'2') THEN
         CALL TRGRN4(INQ)
         CALL TRGRN5(INQ)
         CALL TRGRN6(INQ)
      END IF

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRGRN0

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : NBI chords part 1

!     ***********************************************************

      SUBROUTINE TRGRN1(INQ)

      USE TRCOMM, ONLY : GPNB, GRM, GYR, NRMAX, NRMP, RTG
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR

      CHARACTER(LEN=40):: KFID
      CHARACTER(LEN=5) :: KRTG

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,1)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,1)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,1)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,1)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(1)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,2)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,2)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,2)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,2)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(2)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,3)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,3)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,3)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,3)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(3)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,4)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,4)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,4)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,4)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(4)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRN1

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : NBI chords part 2

!     ***********************************************************

      SUBROUTINE TRGRN2(INQ)

      USE TRCOMM, ONLY : GPNB, GRM, GYR, NRMAX, NRMP, RTG
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR
      CHARACTER(LEN=40) :: KFID
      CHARACTER(LEN=5)  :: KRTG

      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,5)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,5)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,5)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,5)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(5)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,6)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,6)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,6)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,6)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(6)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,7)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,7)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,7)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,7)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(7)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,8)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,8)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,8)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,8)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(8)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRN2

!     ***********************************************************

!           GRAPHIC : RADIAL PROFILE : NBI chords part 3

!     ***********************************************************

      SUBROUTINE TRGRN3(INQ)

      USE TRCOMM, ONLY : GPNB, GRM, GYR, NRMAX, NRMP, RTG
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NR
      CHARACTER(LEN=40) :: KFID
      CHARACTER(LEN=5)  :: KRTG


      CALL PAGES

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,9)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,9)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,9)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,9)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(9)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      DO NR=1,NRMAX
         GYR(NR,1)=GPNB(NR        ,10)*1.E-6
         GYR(NR,2)=GPNB(NR+  NRMAX,10)*1.E-6
         GYR(NR,3)=GPNB(NR+2*NRMAX,10)*1.E-6
         GYR(NR,4)=GPNB(NR+3*NRMAX,10)*1.E-6
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(10)
      KFID='@PNB [MW/m$+3$=] vs r, RTG='//KRTG//' m@'
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,4,KFID,2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRN3

!     ***********************************************************

!           GRAPHIC : PROFILE ALONG THE LINE OF SIGHT: NBI chords part 1

!     ***********************************************************

      SUBROUTINE TRGRN4(INQ)

      USE TRCOMM, ONLY : GBAN, GBL, GBP1, GBR, GBRH, NLM, NLMAX, RTG
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NA, NL
      REAL(4),DIMENSION(NLM)  :: GYBLA
      REAL(4),DIMENSION(NLM,4):: GYBLB
      CHARACTER(LEN=50) :: KFID
      CHARACTER(LEN=5)  :: KRTG

      CALL PAGES

      NA=1
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD( 3.0,12.0,11.0,17.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=2
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD(15.5,24.5,11.0,17.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=3
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD( 3.0,12.0, 2.0, 8.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=4
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD(15.5,24.5, 2.0, 8.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRN4

!     ***********************************************************

!           GRAPHIC : PROFILE ALONG THE LINE OF SIGHT: NBI chords part 2

!     ***********************************************************

      SUBROUTINE TRGRN5(INQ)

      USE TRCOMM, ONLY : GBAN, GBL, GBP1, GBR, GBRH, NLM, NLMAX, RTG
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NA, NL
      REAL(4),DIMENSION(NLM)  :: GYBLA
      REAL(4),DIMENSION(NLM,3):: GYBLB
      CHARACTER(LEN=50) :: KFID
      CHARACTER(LEN=5)  :: KRTG


      CALL PAGES

      NA=5
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD( 3.0,12.0,11.0,17.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=6
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD(15.5,24.5,11.0,17.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=7
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD( 3.0,12.0, 2.0, 8.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=8
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD(15.5,24.5, 2.0, 8.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRN5

!     ***********************************************************

!           GRAPHIC : PROFILE ALONG THE LINE OF SIGHT: NBI chords part 3

!     ***********************************************************

      SUBROUTINE TRGRN6(INQ)

      USE TRCOMM, ONLY : GBAN, GBL, GBP1, GBR, GBRH, NLM, NLMAX, RTG
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INQ
      INTEGER(4) :: NA, NL
      REAL(4),DIMENSION(NLM)  :: GYBLA
      REAL(4),DIMENSION(NLM,3):: GYBLB
      CHARACTER(LEN=50) :: KFID
      CHARACTER(LEN=5)  :: KRTG

      CALL PAGES

      NA=9
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD( 3.0,12.0,11.0,17.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      NA=10
      DO NL=1,NLMAX(NA)
         GYBLA(NL)  =GBR (NL,NA)
         GYBLB(NL,1)=GBRH(NL,NA)
         GYBLB(NL,2)=GBP1(NL,NA)*1.E2
         GYBLB(NL,3)=GBAN(NL,NA)
      ENDDO
      WRITE(KRTG,'(F5.3)') RTG(NA)
      KFID='@(L)R, (R)RHO, P1*10$+2$=, ANL vs r, RTG='//KRTG//' m@'
      CALL TRGR1DD(15.5,24.5,11.0,17.0,GBL,GYBLA,GYBLB,NLM,NLMAX(NA),1,3,KFID,3+INQ,2+INQ)

      CALL TRGRTM
      CALL PAGEE
      RETURN
      END SUBROUTINE TRGRN6
