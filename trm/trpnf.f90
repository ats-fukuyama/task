!     ***********************************************************

!           Nuclear reaction (DT)

!     ***********************************************************

      SUBROUTINE TRNFDT

      USE TRCOMM, ONLY : &
           AME, AMP, ANC, ANFE, MDLNF, NRMAX, PA, PBIN, PFCL, &
           PFIN, PI, PNBENG, PNF, PZ, PZC, PZFE, RKEV, &
           RN, RNF, RT, RTF, RW, SNF, TAUF, rkind, NNBMAX, PNB_NNB
      IMPLICIT NONE
!      INCLUDE 'trcomm.inc'
      REAL(rkind)   :: &
           AMA, AMD, AMT, ANE, EC, HYF, P1, PTNT, SS, SSB, TAUS, &
           TD, TE, TT, VC3, VCA3, VCD3, VCR, VCT3, VF, WF, ZEFFM, PB
      INTEGER:: NR,NNB
      REAL(rkind)   :: SIGMAM, COULOG, SIGMAB, HY   !FUNCTION


      AMD=PA(2)*AMP
      AMT=PA(3)*AMP
      AMA=PA(4)*AMP
      VF =SQRT(2.D0*3.5D3 *RKEV/AMA)

      DO NR=1,NRMAX
         ANE= RN(NR,1)
         TE = RT(NR,1)
         TD = RT(NR,2)
         TT = RT(NR,3)
!         WRITE(6,*) NR,TD,TT
         SS = SIGMAM(TD,TT)
         IF(MDLNF.GE.3) THEN
            ZEFFM = (PZ(2)*PZ(2)*RN(NR,2)/PA(2) &
                    +PZ(3)*PZ(3)*RN(NR,3)/PA(3) &
                    +PZ(4)*PZ(4)*RN(NR,4)/PA(4) &
                    +PZC(NR)*PZC(NR) *ANC(NR) /12.D0 &
                    +PZFE(NR)*PZFE(NR)*ANFE(NR)/52.D0)/ANE
            EC  = 14.8D0*TE*PA(2)*ZEFFM**(2.D0/3.D0)
            TAUS= 0.2D0*PA(2)*ABS(TE)**1.5D0 /(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))
            SSB=0.D0
            PB=0.D0
            DO NNB=1,NNBMAX
               PTNT= PBIN(NR)*TAUS/(RN(NR,2)*1.D20*PNBENG(NNB)*RKEV)
               SSB = SSB+PNB_NNB(NNB,NR)*SIGMAB(PNBENG(NNB),EC,TT,PTNT)
               PB  = PB +PNB_NNB(NNB,NR)
            END DO
            IF(PB.NE.0.D0) THEN
               SSB=SSB/PB
            ELSE
               SSB=0.D0
            END IF
         ELSE
            SSB=0.D0
         ENDIF
         SNF(NR) = (SS+SSB)*RN(NR,2)*RN(NR,3)*1.D20
         PNF(NR) = SNF(NR)*3.5D3*RKEV*1.D20
         IF(MOD(MDLNF,2).EQ.1) SNF(NR) = 0.D0
      ENDDO

      DO NR=1,NRMAX
         ANE= RN(NR,1)
         TE = RT(NR,1)
         WF = RW(NR,2)
         P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
         VCD3 = P1*RN(NR,2)*PZ(2)**2/AMD
         VCT3 = P1*RN(NR,3)*PZ(3)**2/AMT
         VCA3 = P1*RN(NR,4)*PZ(4)**2/AMA
         VC3  = VCD3+VCT3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYF=HY(VF/VCR)
         TAUS = 0.2D0*PA(4)*ABS(TE)**1.5D0 /(PZ(4)**2*ANE*COULOG(1,2,ANE,TE))
         TAUF(NR)= 0.5D0*TAUS*(1.D0-HYF)
         RNF(NR,2)= 2.D0*LOG(1.D0+(VF/VCR)**3)*WF /(3.D0*(1.D0-HYF)*3.5D3)
         IF(RNF(NR,2).GT.0.D0) THEN
            RTF(NR,2)= WF/RNF(NR,2)
         ELSE
            RTF(NR,2)= 0.D0
         ENDIF
         PFIN(NR) = WF*RKEV*1.D20/TAUF(NR)
         PFCL(NR,1)=    (1.D0-HYF)*PFIN(NR)
         PFCL(NR,2)=(VCD3/VC3)*HYF*PFIN(NR)
         PFCL(NR,3)=(VCT3/VC3)*HYF*PFIN(NR)
         PFCL(NR,4)=(VCA3/VC3)*HYF*PFIN(NR)
      ENDDO

      RETURN
      END SUBROUTINE TRNFDT

!     ***********************************************************

!           REACTION CROSS SECTION (MAXELLIAN) DT

!     ***********************************************************

      FUNCTION SIGMAM(TD,TT)

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) TD,TT,SIGMAM
      REAL(rkind) TI,H,ARG

      TI = (3.D0*ABS(TD)+2.D0*ABS(TT))/5.D0
      H  = TI/37.D0 + 5.45D0/(3.D0+TI*(1.D0+(TI/37.5D0)**2.8D0))
      ARG= -20.D0/TI**(1.D0/3.D0)
      IF(ARG.GE.-100.D0)  THEN
         SIGMAM = 3.7D-18*TI**(-2.D0/3.D0)*EXP(ARG)/H
      ELSE
         SIGMAM = 0.D0
      ENDIF

      RETURN
      END FUNCTION SIGMAM

!     ***********************************************************

!           REACTION CROSS SECTION (MAXELLIAN) DD

!     ***********************************************************

      FUNCTION SIGMAMDD(TD)

        ! not completed

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) TD,SIGMAMDD
      REAL(rkind) TI,H,ARG

      H  = TI/37.D0 + 5.45D0/(3.D0+TI*(1.D0+(TI/37.5D0)**2.8D0))
      ARG= -20.D0/TI**(1.D0/3.D0)
      IF(ARG.GE.-100.D0)  THEN
         SIGMAMDD = 3.7D-18*TI**(-2.D0/3.D0)*EXP(ARG)/H
      ELSE
         SIGMAMDD = 0.D0
      ENDIF

      RETURN
    END FUNCTION SIGMAMDD

!     ***********************************************************

!      REACTION RATE : TAIL

!     ***********************************************************

      FUNCTION SIGMAB(EB,EC,TI,PTNT)

!      APPROXIMATE FORMULA OF FUSION REACTION RATE
!         FOR SLOWING DOWN ION DISTRIBUTION
!      REF. TAKIZUKA AND YAMAGIWA, JAERI-M 87-066

!      EB   : BEAM ENERGY (KEV)
!      EC   : CRITICAL ENERGY (KEV)
!      TI   : TRITIUM TEMPERATURE (KEV)
!      PTNT : PB * TAUS / (ND * EB)

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) EB,EC,TI,PTNT,SIGMAB
      REAL(rkind) SIGMBS
      REAL(rkind) XB,XC,AG1,AG2,AG3,AL1,AL2,AL3,X1,X2,X3,X4,SA

      XB=SQRT(EB/127.D0)
      XC=SQRT(EC/127.D0)

      AG1= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG2= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG3= 0.33D0
      AL1= 1.D0/(0.40D0+0.032D0*TI)
      AL2=-1.D0/(0.91D0+0.016D0*SQRT(ABS(TI)))
      AL3=-0.11D0
      X1=0.97D0-AG1/AL1
      X2=0.97D0
      X3=0.97D0+(AG2-AG3)/(AL3-AL2)
      X4=0.97D0+3.D0

      IF(XB.LT.X1) THEN
         SA=0.D0
      ELSE
         SA=-SIGMBS(X1,AG1,AL1,XC)
         IF(XB.LT.X2) THEN
            SA=SA+SIGMBS(XB,AG1,AL1,XC)
         ELSE
            SA=SA+SIGMBS(X2,AG1,AL1,XC)-SIGMBS(X2,AG2,AL2,XC)
            IF(XB.LT.X3) THEN
               SA=SA+SIGMBS(XB,AG2,AL2,XC)
            ELSE
               SA=SA+SIGMBS(X3,AG2,AL2,XC)-SIGMBS(X3,AG3,AL3,XC)
               IF(XB.LT.X4) THEN
                  SA=SA+SIGMBS(XB,AG3,AL3,XC)
               ELSE
                  SA=SA+SIGMBS(X4,AG3,AL3,XC)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      SIGMAB=PTNT*1.67D-21*SA
      RETURN
      END FUNCTION SIGMAB

      FUNCTION SIGMBS(XX,RGG,RGL,XC)

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) :: XX,RGG,RGL,XC,SIGMBS
      REAL(rkind) :: X

      X=XX/XC
      SIGMBS=((RGG-0.97D0*RGL)/3.D0+XC*RGL/6.D0)*LOG(X*X*X+1.D0) &
     &      +XC*RGL*(X-LOG(X+1.D0)/2.D0 -ATAN((2.D0*X-1.D0)/SQRT(3.D0))/SQRT(3.D0))
      RETURN
      END FUNCTION SIGMBS

!     ***********************************************************

!           Nuclear reaction (D-He3)

!     ***********************************************************

      SUBROUTINE TRNFDHe3

      USE TRCOMM, ONLY : &
           AME, AMP, MDLNF, NRMAX, PA, PFCL, &
           PFIN, PI, PNF, PZ, RKEV, &
           RN, RNF, RT, RTF, RW, SNF, TAUF, rkind
      IMPLICIT NONE
      REAL(rkind)   :: &
           AMA, AMD, AMHe3, ANE, HYF, P1, SS, TAUS, &
           TD, TE, THe3, VC3, VCA3, VCD3, VCR, VCHe3, VF, WF
      INTEGER:: NR
      REAL(rkind)   :: SIGMADHe3, COULOG, HY   !FUNCTION

      AMD=  PA(2)*AMP
      AMHe3=PA(3)*AMP
      AMA=  PA(4)*AMP
      VF =SQRT(2.D0*3.6D3 *RKEV/AMA)

      DO NR=1,NRMAX
         ANE  = RN(NR,1)
         TE   = RT(NR,1)
         TD   = RT(NR,2)
         THe3 = RT(NR,3)
         SS = SIGMADHe3(TD,THe3)
         SNF(NR) = SS*RN(NR,2)*RN(NR,3)*1.D20
         PNF(NR) = SNF(NR)*(3.5D3+14.7D3)*RKEV*1.D20  ! proton energy added
                                                      ! for simplicity
         IF(MOD(MDLNF,2).EQ.1) SNF(NR) = 0.D0
      ENDDO

      DO NR=1,NRMAX
         ANE= RN(NR,1)
         TE = RT(NR,1)
         WF = RW(NR,2)
         P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
         VCD3  = P1*RN(NR,2)*PZ(2)**2/AMD
         VCHe3 = P1*RN(NR,3)*PZ(3)**2/AMHe3
         VCA3  = P1*RN(NR,4)*PZ(4)**2/AMA
         VC3  = VCD3+VCHe3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYF=HY(VF/VCR)
         TAUS = 0.2D0*PA(4)*ABS(TE)**1.5D0 /(PZ(4)**2*ANE*COULOG(1,2,ANE,TE))
         TAUF(NR)= 0.5D0*TAUS*(1.D0-HYF)
         RNF(NR,2)= 2.D0*LOG(1.D0+(VF/VCR)**3)*WF /(3.D0*(1.D0-HYF)*3.6D3)
         IF(RNF(NR,2).GT.0.D0) THEN
            RTF(NR,2)= WF/RNF(NR,2)
         ELSE
            RTF(NR,2)= 0.D0
         ENDIF
         PFIN(NR) = WF*RKEV*1.D20/TAUF(NR)
         PFCL(NR,1)=    (1.D0-HYF)*PFIN(NR)
         PFCL(NR,2)=(VCD3 /VC3)*HYF*PFIN(NR)
         PFCL(NR,3)=(VCHe3/VC3)*HYF*PFIN(NR)
         PFCL(NR,4)=(VCA3 /VC3)*HYF*PFIN(NR)
      ENDDO

      RETURN
      END SUBROUTINE TRNFDHe3

!     ***********************************************************

!           REACTION CROSS SECTION (MAXELLIAN)

!     ***********************************************************

      FUNCTION SIGMADHe3(TD,THe3)

      USE trcomm,ONLY: rkind
      USE libspl1d
      IMPLICIT NONE
      REAL(rkind) TD,THe3,TI,TIL,XRATEL,SIGMADHe3
      REAL(rkind),DIMENSION(10),save:: RENG,RRATE
      REAL(rkind),DIMENSION(10),save:: RENGL,RRATEL,DIFF
      REAL(rkind),DIMENSION(4,10),save:: URRATE
      INTEGER:: NX,IERR
      INTEGER,save:: INIT=0
      DATA RENG/1.D0, 2.D0, 5.D0, 10.D0, 20.D0, &
                50.D0,100.D0,200.D0,500.D0,1000.D0/
      DATA RRATE/1.0D-32, 1.4D-29, 6.7D-27, 2.3D-25, 3.8D-24, &
                 5.4D-23, 1.6D-22, 2.4D-22, 2.3D-22, 1.8D-22/

      IF(INIT.EQ.0) THEN
         
         DO NX=1,10
            RENGL(NX)=LOG(RENG(NX))
            RRATE(NX)=LOG(RRATEL(NX))
         ENDDO
         CALL SPL1D(RENGL,RRATEL,DIFF,URRATE,10,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SIGMADHe3: SPL1D: IERR=',IERR
         INIT=1
      ENDIF

      TI = (3.D0*ABS(TD)+2.D0*ABS(THe3))/5.D0
      IF(TI.GT.0.D0) THEN
         TIL=LOG(MAX(TI,1.D0))
         CALL SPL1DF(TIL,XRATEL,RENGL,URRATE,10,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SIGMADHe3: SPL1DF: IERR=',IERR
         IF(IERR.NE.0) WRITE(6,*) TIL,RENGL(1),RENGL(10)
         IF(IERR.NE.0) WRITE(6,*) TI,TD,THe3
         SIGMADHe3=EXP(XRATEL)
      ELSE
         SIGMADHe3=0.D0
      ENDIF
      RETURN
      END FUNCTION SIGMADHe3
