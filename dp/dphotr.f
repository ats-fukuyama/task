C     $Id$
C
C ******************************************************
C                       DPHOTR
C      dielectric tensor with relativistic effects
C                 
C                      94/05/14             
C                           programed by K.TANAKA
C
C      revised version
C
C                      96/01/17
C                           programed by N.KAIHARA
C ******************************************************
C
      SUBROUTINE DPHOTR(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION CLDISP(6),CLDISP1(6),CLDISP2(6)
C
      CALL DPHOTRR(CW,CKPR,CKPP,NS,CLDISP1)
      CALL DPHOTRI(CW,CKPR,CKPP,NS,CLDISP2)
      DO 100 I=1,6
         CLDISP(I)=CLDISP1(I)+CLDISP2(I)
  100 CONTINUE
      RETURN
      END
C
C ******************************************************
C                       DPHOTRR
C ******************************************************
C
      SUBROUTINE DPHOTRR(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION CLDISP(6)
C
      NCMIN = NDISP1(NS)
      NCMAX = NDISP2(NS)
      NHMAX=MAX(ABS(NCMIN),ABS(NCMAX),2)+5
      DELPL  = 0.5D0
C
      CWP=RNE0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      DKPP=DBLE(CKPP)
C
      DO 50 NTH=1,NTHMAX
      DO 50 NP=1,NPMAX-1
         DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP
   50 CONTINUE 
      DO 60 NP=1,NPMAX
      DO 60 NTH=1,NTHMAX-1
         DFT(NP,NTH) = (FM(NP,NTH+1) - FM(NP,NTH))/DELTH
   60 CONTINUE
C
C***********DGP1,DGP2,DGT1,DGT2************
C
      DO 65 NP=1,NPMAX
         DO 65 NTH=1,NTHMAX
         DGP1(NP,NTH)=PTH0W*PG(NP)/SQRT(1+PTH0W*PG(NP)**2)
     &               -DKPRW*TCSM(NTH)
         DGP2(NP,NTH)=PTH0W*PM(NP)/SQRT(1+PTH0W*PM(NP)**2)
     &               -DKPRW*TCSG(NTH)
         DGT1(NP,NTH)=DKPRW*PG(NP)*TSNM(NTH)
         DGT2(NP,NTH)=DKPRW*PM(NP)*TSNG(NTH)
   65 CONTINUE 
C
C*****************PRINCIPAL VALUE***********************
C
C*************SUM1********************
C
      CINTG111 = (0.D0,0.D0)
      CINTG112 = (0.D0,0.D0)
      CINTG113 = (0.D0,0.D0)
      CINTG121 = (0.D0,0.D0)
      CINTG122 = (0.D0,0.D0)
      CINTG123 = (0.D0,0.D0)
      CINTG131 = (0.D0,0.D0)
      CINTG132 = (0.D0,0.D0)
      CINTG133 = (0.D0,0.D0)
C
      DO 120 NP=1,NPMAX-1
      DO 120 NTH=1,NTHMAX
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
C
         X = DKPP*PTH0*PG(NP)*TSNM(NTH)/WCM
         CALL BESSJN(X,NHMAX,ADJ,ADJD)
C
         DO 110 NC=NCMIN,NCMAX
            NCD = ABS(NC)
            CDENX= RGMG(NP)-CKPRW*PG(NP)*TCSM(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP1(NP,NTH)*DELP)**2
     &                             +(DELPL*DGT1(NP,NTH)*DELTH)**2)
C
            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            PAI1  = NC*INC*ADJ(NCD)/X
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)
C
            CSM11 = CSM11 + PAI1          *PAI1*CDEN
            CSM12 = CSM12 + PAI1         *CPAI2*CDEN
            CSM13 = CSM13 + PAI1          *PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CDEN
            CSM33 = CSM33 + PAI3          *PAI3*CDEN
  110    CONTINUE
C
         PART1= DFP(NP,NTH)*PG(NP)*PG(NP)*PG(NP)
     &                     *TSNM(NTH)*TSNM(NTH)*TSNM(NTH)
     &         *DELTH*DELP
C     
         CINTG111= CINTG111 + CSM11*PART1
         CINTG112= CINTG112 + CSM12*PART1
         CINTG113= CINTG113 + CSM13*PART1
         CINTG121= CINTG121 - CSM12*PART1
         CINTG122= CINTG122 + CSM22*PART1
         CINTG123= CINTG123 + CSM23*PART1
         CINTG131= CINTG131 + CSM13*PART1
         CINTG132= CINTG132 - CSM23*PART1
         CINTG133= CINTG133 + CSM33*PART1
  120 CONTINUE 
C
C*************SUM2********************
C
      CINTG211 = (0.D0,0.D0)
      CINTG212 = (0.D0,0.D0)
      CINTG213 = (0.D0,0.D0)
      CINTG221 = (0.D0,0.D0)
      CINTG222 = (0.D0,0.D0)
      CINTG223 = (0.D0,0.D0)
      CINTG231 = (0.D0,0.D0)
      CINTG232 = (0.D0,0.D0)
      CINTG233 = (0.D0,0.D0)
C
      DO 140 NP=1,NPMAX
      DO 140 NTH=1,NTHMAX-1
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
C
         X = DKPP*PTH0*PM(NP)*TSNG(NTH)/WCM
         CALL BESSJN(X,NHMAX,ADJ,ADJD)
C
         DO 130 NC=NCMIN,NCMAX
            NCD = ABS(NC)
            CDENX = RGMM(NP)-CKPRW*PM(NP)*TCSG(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP2(NP,NTH)*DELP)**2
     &                             +(DELPL*DGT2(NP,NTH)*DELTH)**2)
C
            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            PAI1  = NC*INC*ADJ(NCD)/X
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)
C
            CSM11 = CSM11 + PAI1          *PAI1*CDEN
            CSM12 = CSM12 + PAI1         *CPAI2*CDEN
            CSM13 = CSM13 + PAI1          *PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CDEN
            CSM33 = CSM33 + PAI3          *PAI3*CDEN
  130       CONTINUE
C 
         CPART2= DFT(NP,NTH)*PM(NP)*PM(NP)
     &                      *TSNG(NTH)*TSNG(NTH)
     &          *(TCSG(NTH)-CKPRW*PM(NP)/RGMM(NP))
     &          *DELTH*DELP
C 
         CINTG211= CINTG211 + CSM11*CPART2
         CINTG212= CINTG212 + CSM12*CPART2
         CINTG213= CINTG213 + CSM13*CPART2
         CINTG221= CINTG221 - CSM12*CPART2
         CINTG222= CINTG222 + CSM22*CPART2
         CINTG223= CINTG223 + CSM23*CPART2
         CINTG231= CINTG231 + CSM13*CPART2
         CINTG232= CINTG232 - CSM23*CPART2
         CINTG233= CINTG233 + CSM33*CPART2
     &                      - PM(NP)*PM(NP)*TCSG(NTH)
     &                        *DFT(NP,NTH)/RGMM(NP)
     &                        *DELTH*DELP
  140 CONTINUE 
C
         FACT=2.D0*PI*DBLE(CWP)
C
         CLDISP(1)=FACT*(CINTG111+CINTG211)
         CLDISP(2)=FACT*(CINTG133+CINTG233)-CLDISP(1)
         CLDISP(3)=FACT*(CINTG122+CINTG222)-CLDISP(1)
         CLDISP(4)=FACT*(CINTG113+CINTG213)
         CLDISP(5)=FACT*(CINTG112+CINTG212)
         CLDISP(6)=FACT*(CINTG123+CINTG223)
C 
      RETURN
      END                     
C
C ******************************************************
C                       DPHOTRI
C ******************************************************
C
      SUBROUTINE DPHOTRI(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
C
      NCMIN = NDISP1(NS)
      NCMAX = NDISP2(NS)
C
      CWP=RNE0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      DCWC=DBLE(CWC)
      DNPR=DBLE(CKPR*VC/CW)
      DKPP=DBLE(CKPP)
C
      DO 50 NTH=1,NTHMAX
      DO 50 NP=1,NPMAX-1
         DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP
   50 CONTINUE 
      DO 60 NP=1,NPMAX
      DO 60 NTH=1,NTHMAX-1
         DFT(NP,NTH) = (FM(NP,NTH+1) - FM(NP,NTH))/DELTH
   60 CONTINUE
C
C***************SINGULAR POINT***************************
C
C
C*****************SUM3************************
C
      CINTG311 = (0.D0,0.D0)
      CINTG312 = (0.D0,0.D0)
      CINTG313 = (0.D0,0.D0)
      CINTG321 = (0.D0,0.D0)
      CINTG322 = (0.D0,0.D0)
      CINTG323 = (0.D0,0.D0)
      CINTG331 = (0.D0,0.D0)
      CINTG332 = (0.D0,0.D0)
      CINTG333 = (0.D0,0.D0)
C
      DO 300 NTH=1,NTHMAX
C
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
C
         DO 310 NC=NCMIN,NCMAX
C
            D = (NC*DCWC)**2+(DNPR*TCSM(NTH))**2
            IF(D.LE.1.D0) GOTO 310
C
C   PNEAR1
C
            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2)*
     &               (DNPR*NC*DCWC*TCSM(NTH)+SQRT(D-1))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP*NPMAX) THEN
               GOTO 302
            ELSEIF (DKPRW*TCSM(NTH)*PNEAR1+NC*DCWC.LT.0.D0) THEN
               GOTO 302
            END IF
            NP1 = INT(PNEAR1/DELP)
            IF (NP1.LT.0.OR.NP1.GE.NPMAX) THEN
               GOTO 302
            ELSEIF (NP1.EQ.0) THEN
               DIF = PNEAR1/DELP
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP1.EQ.NPMAX-1) THEN
               DIF = (PNEAR1 - PG(NP1))/DELP
               DFP3 = (1.D0-DIF)*DFP(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PG(NP1))/DELP
               DFP3 = DIF*DFP(NP1+1,NTH)+(1.D0-DIF)*DFP(NP1,NTH)
            ENDIF
C
            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNM(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)
C
            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART31= DFP3*PNEAR1**3
     &               /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSM(NTH))
C
            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            PAI1  = NC*INC*ADJ(NCD)/X
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)
C
            CSM11 = CSM11 + PAI1          *PAI1*CPART31
            CSM12 = CSM12 + PAI1         *CPAI2*CPART31
            CSM13 = CSM13 + PAI1          *PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART31
            CSM33 = CSM33 + PAI3          *PAI3*CPART31
C
C     PNEAR2
C
  302       PNEAR2 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2)*
     &               (DNPR*NC*DCWC*TCSM(NTH)-SQRT(D-1))/PTH0
            IF (PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP*NPMAX) THEN 
               GOTO 310
            ELSEIF (DKPRW*TCSM(NTH)*PNEAR2+NC*DCWC.LT.0.D0) THEN
               GOTO 310
            END IF
            NP2 = INT(PNEAR2/DELP)
            IF (NP2.LT.0.OR.NP2.GE.NPMAX) THEN
               GOTO 310
            ELSE IF(NP2.EQ.0) THEN
               DIF = PNEAR2/DELP
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP2.EQ.NPMAX-1) THEN
               DIF = (PNEAR2 - PG(NP2))/DELP
               DFP3 = (1.D0-DIF)*DFP(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PG(NP2))/DELP
               DFP3 = DIF*DFP(NP2+1,NTH)+(1.D0-DIF)*DFP(NP2,NTH)
            ENDIF
C
            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNM(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)
C
            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART31= DFP3*PNEAR2**3
     &               /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSM(NTH))
C
            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            PAI1  = NC*INC*ADJ(NCD)/X
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)
C
            CSM11 = CSM11 + PAI1          *PAI1*CPART31
            CSM12 = CSM12 + PAI1         *CPAI2*CPART31
            CSM13 = CSM13 + PAI1          *PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART31
            CSM33 = CSM33 + PAI3          *PAI3*CPART31
  310    CONTINUE 
C
         CPART32= -CI*PI*TSNM(NTH)**3*DELTH
C              
         CINTG311 = CINTG311 + CSM11*CPART32
         CINTG312 = CINTG312 + CSM12*CPART32
         CINTG313 = CINTG313 + CSM13*CPART32
         CINTG321 = CINTG321 - CSM12*CPART32
         CINTG322 = CINTG322 + CSM22*CPART32
         CINTG323 = CINTG323 + CSM23*CPART32
         CINTG331 = CINTG331 + CSM13*CPART32
         CINTG332 = CINTG332 - CSM23*CPART32
         CINTG333 = CINTG333 + CSM33*CPART32
  300 CONTINUE
C 
C*****************SUM4************************
C 
      CINTG411 = (0.D0,0.D0)
      CINTG412 = (0.D0,0.D0)
      CINTG413 = (0.D0,0.D0)
      CINTG421 = (0.D0,0.D0)
      CINTG422 = (0.D0,0.D0)
      CINTG423 = (0.D0,0.D0)
      CINTG431 = (0.D0,0.D0)
      CINTG432 = (0.D0,0.D0)
      CINTG433 = (0.D0,0.D0)
C
      DO 400 NTH=1,NTHMAX-1
C               
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
C               
         DO 410 NC=NCMIN,NCMAX
C
            D = (NC*DCWC)**2+(DNPR*TCSG(NTH))**2
            IF(D.LE.1.D0) GOTO 410
C
C   PNEAR1
C
            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2)*
     &               (DNPR*NC*DCWC*TCSG(NTH)+SQRT(D-1))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP*NPMAX) THEN
               GOTO 402
            ELSEIF (DKPRW*TCSG(NTH)*PNEAR1+NC*DCWC.LT.0.D0) THEN
               GOTO 402
            END IF
            NP1 = INT(PNEAR1/DELP+0.5D0)
            IF (NP1.LT.0.OR.NP1.GE.NPMAX) THEN
               GOTO 402
            ELSE IF(NP1.EQ.0) THEN
               DIF = (PNEAR1 - PM(1))/DELP
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP1.EQ.NPMAX-1) THEN
               DIF = (PNEAR1 - PM(NP1))/DELP
               DFT4  = (1.D0-DIF)*DFT(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PM(NP1))/DELP
               DFT4  = DIF*DFT(NP1+1,NTH)+(1.D0-DIF)*DFT(NP1,NTH)
            ENDIF
C
            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNG(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)
C
            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART41=DFT4*PNEAR1**2*(TCSG(NTH)-PNEAR1*CKPRW/RGM)
     &              /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSG(NTH))
C
            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            PAI1  = NC*INC*ADJ(NCD)/X
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)
C
            CSM11 = CSM11 + PAI1          *PAI1*CPART41
            CSM12 = CSM12 + PAI1         *CPAI2*CPART41
            CSM13 = CSM13 + PAI1          *PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART41
            CSM33 = CSM33 + PAI3          *PAI3*CPART41
C
C   PNEAR2
C
 402        PNEAR2 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2)*
     &               (DNPR*NC*DCWC*TCSG(NTH)-SQRT(D-1))/PTH0
            IF(PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP*NPMAX) THEN
               GOTO 410
            ELSEIF (DKPRW*TCSG(NTH)*PNEAR2+NC*DCWC.LT.0.D0) THEN
               GOTO 410
            END IF
            NP2 = INT(PNEAR2/DELP+0.5D0)
            IF (NP2.LT.0.OR.NP2.GE.NPMAX) THEN
               GOTO 410
            ELSE IF(NP2.EQ.0) THEN
               DIF = (PNEAR2 - PM(1))/DELP
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP2.EQ.NPMAX-1) THEN
               DIF = (PNEAR2 - PM(NP2))/DELP
               DFT4  = (1.D0-DIF)*DFT(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PM(NP2))/DELP
               DFT4  = DIF*DFT(NP2+1,NTH)+(1.D0-DIF)*DFT(NP2,NTH)
            ENDIF
C
            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNG(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)
C
            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART41=DFT4*PNEAR2**2*(TCSG(NTH)-PNEAR2*CKPRW/RGM)
     &              /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSG(NTH))
C
            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            PAI1  = NC*INC*ADJ(NCD)/X
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)
C
            CSM11 = CSM11 + PAI1          *PAI1*CPART41
            CSM12 = CSM12 + PAI1         *CPAI2*CPART41
            CSM13 = CSM13 + PAI1          *PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART41
            CSM33 = CSM33 + PAI3          *PAI3*CPART41

  410    CONTINUE 
C
         CPART42= -CI*PI*TSNG(NTH)**2*DELTH
C
         CINTG411 = CINTG411 + CSM11*CPART42
         CINTG412 = CINTG412 + CSM12*CPART42
         CINTG413 = CINTG413 + CSM13*CPART42
         CINTG421 = CINTG421 - CSM12*CPART42
         CINTG422 = CINTG422 + CSM22*CPART42
         CINTG423 = CINTG423 + CSM23*CPART42
         CINTG431 = CINTG431 + CSM13*CPART42
         CINTG432 = CINTG432 - CSM23*CPART42
         CINTG433 = CINTG433 + CSM33*CPART42
  400 CONTINUE
C
      FACT=2.D0*PI*DBLE(CWP)
C
      CLDISP(1)=FACT*(CINTG311+CINTG411)
      CLDISP(2)=FACT*(CINTG333+CINTG433)-CLDISP(1)
      CLDISP(3)=FACT*(CINTG322+CINTG422)-CLDISP(1)
      CLDISP(4)=FACT*(CINTG313+CINTG413)
      CLDISP(5)=FACT*(CINTG312+CINTG412)
      CLDISP(6)=FACT*(CINTG323+CINTG423)
C
C      WRITE(6,'(1P6E12.4)') CLDISP(1),CLDISP(2),CLDISP(3)
C      WRITE(6,'(1P6E12.4)') CLDISP(4),CLDISP(5),CLDISP(6)
C 
      RETURN
      END                     
