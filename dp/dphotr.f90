!
! ******************************************************
!                       DPHOTR
!      dielectric tensor with relativistic effects
!                 
!                      94/05/14             
!                           programed by K.TANAKA
!
!      revised version
!
!                      96/01/17
!                           programed by N.KAIHARA
! ******************************************************

MODULE dphotr

  PRIVATE
  PUBLIC DP_HOTR,DP_HOTRR,DP_HOTRI

CONTAINS

  SUBROUTINE DP_HOTR(CW,CKPR,CKPP,NSA,mag,CLDISP)

    USE dpcomm
    USE plprof
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NSA
    TYPE(pl_mag_type),INTENT(IN):: mag
!    TYPE(pl_prf_type),DIMENSION(nsmax),INTENT(IN):: plf
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    COMPLEX(rkind):: CLDISP1(6),CLDISP2(6)
    INTEGER:: I
      
    CALL DP_HOTRR(CW,CKPR,CKPP,NSA,mag,CLDISP1)
    CALL DP_HOTRI(CW,CKPR,CKPP,NSA,mag,CLDISP2)
    DO I=1,6
       CLDISP(I)=CLDISP1(I)+CLDISP2(I)
    ENDDO
    RETURN
  END SUBROUTINE DP_HOTR

! ******************************************************
!                       DPHOTRR
! ******************************************************

  SUBROUTINE DP_HOTRR(CW,CKPR,CKPP,NSA,mag,CLDISP)

    USE dpcomm
    USE plprof
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NSA
    TYPE(pl_mag_type),INTENT(IN):: mag
!    TYPE(pl_prf_type),DIMENSION(nsmax),INTENT(IN):: plf
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NHMAX,NTH,NP,NC,NCD,INC,NS
    REAL(rkind):: DELPL,WCM,DKPRW,PTH0W,DKPP,X,PAI1,PAI3,PART1,FACT
    REAL(rkind):: PN0,PT0,PTH0
    COMPLEX(rkind):: CWP,CWC,CKPRW,CDENX,CDEN,CPAI2,CPART2
    COMPLEX(rkind):: CINTG111,CINTG112,CINTG113
    COMPLEX(rkind):: CINTG121,CINTG122,CINTG123
    COMPLEX(rkind):: CINTG131,CINTG132,CINTG133
    COMPLEX(rkind):: CINTG211,CINTG212,CINTG213
    COMPLEX(rkind):: CINTG221,CINTG222,CINTG223
    COMPLEX(rkind):: CINTG231,CINTG232,CINTG233
    COMPLEX(rkind):: CSM11,CSM12,CSM13,CSM22,CSM23,CSM33

    NS=NS_NSA_DP(NS)
      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)),2)+2
      ALLOCATE(ADJ(0:NHMAX),ADJD(0:NHMAX))
      DELPL  = 0.5D0

      PN0=RNFP0(NSA)
      PT0=RTFP0(NSA)
      PTH0=SQRT(PT0*1.D3*AEE*AMFP(NSA))

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=mag%BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      DKPP=DBLE(CKPP)

      DO NTH=1,NTHMAX_DP
      DO NP=1,NPMAX_DP-1
         DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP(NS)
      ENDDO
      ENDDO
      DO NP=1,NPMAX_DP
      DO NTH=1,NTHMAX_DP-1
         DFT(NP,NTH) = (FM(NP,NTH+1) - FM(NP,NTH))/DELTH
      ENDDO
      ENDDO

!***********DGP1,DGP2,DGT1,DGT2************

      DO NP=1,NPMAX_DP
      DO NTH=1,NTHMAX_DP
         DGP1(NP,NTH)=PTH0W*PG(NP,NSA)/SQRT(1+PTH0W*PG(NP,NSA)**2) &
                     -DKPRW*TCSM(NTH)
         DGP2(NP,NTH)=PTH0W*PM(NP,NSA)/SQRT(1+PTH0W*PM(NP,NSA)**2) &
                     -DKPRW*TCSG(NTH)
         DGT1(NP,NTH)=DKPRW*PG(NP,NSA)*TSNM(NTH)
         DGT2(NP,NTH)=DKPRW*PM(NP,NSA)*TSNG(NTH)
      ENDDO
      ENDDO

!*****************PRINCIPAL VALUE***********************

!*************SUM1********************

      CINTG111 = (0.D0,0.D0)
      CINTG112 = (0.D0,0.D0)
      CINTG113 = (0.D0,0.D0)
      CINTG121 = (0.D0,0.D0)
      CINTG122 = (0.D0,0.D0)
      CINTG123 = (0.D0,0.D0)
      CINTG131 = (0.D0,0.D0)
      CINTG132 = (0.D0,0.D0)
      CINTG133 = (0.D0,0.D0)

      DO NP=1,NPMAX_DP-1
      DO NTH=1,NTHMAX_DP
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         X = DKPP*PTH0*PG(NP,NSA)*TSNM(NTH)/WCM
         CALL BESSJN(X,NHMAX,ADJ,ADJD)

         DO NC=NCMIN(NS),NCMAX(NS)
            NCD = ABS(NC)
            CDENX= RGMG(NP,NSA)-CKPRW*PG(NP,NSA)*TCSM(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP1(NP,NTH)*DELP(NSA))**2 &
                                   +(DELPL*DGT1(NP,NTH)*DELTH)**2)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CDEN
            CSM12 = CSM12 + PAI1         *CPAI2*CDEN
            CSM13 = CSM13 + PAI1          *PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CDEN
            CSM33 = CSM33 + PAI3          *PAI3*CDEN
         ENDDO

         PART1= DFP(NP,NTH)*PG(NP,NSA)*PG(NP,NSA)*PG(NP,NSA) &
                           *TSNM(NTH)*TSNM(NTH)*TSNM(NTH) &
               *DELTH*DELP(NSA)

         CINTG111= CINTG111 + CSM11*PART1
         CINTG112= CINTG112 + CSM12*PART1
         CINTG113= CINTG113 + CSM13*PART1
         CINTG121= CINTG121 - CSM12*PART1
         CINTG122= CINTG122 + CSM22*PART1
         CINTG123= CINTG123 + CSM23*PART1
         CINTG131= CINTG131 + CSM13*PART1
         CINTG132= CINTG132 - CSM23*PART1
         CINTG133= CINTG133 + CSM33*PART1
      ENDDO
      ENDDO

!*************SUM2********************

      CINTG211 = (0.D0,0.D0)
      CINTG212 = (0.D0,0.D0)
      CINTG213 = (0.D0,0.D0)
      CINTG221 = (0.D0,0.D0)
      CINTG222 = (0.D0,0.D0)
      CINTG223 = (0.D0,0.D0)
      CINTG231 = (0.D0,0.D0)
      CINTG232 = (0.D0,0.D0)
      CINTG233 = (0.D0,0.D0)

      DO NP=1,NPMAX_DP
      DO NTH=1,NTHMAX_DP-1
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         X = DKPP*PTH0*PM(NP,NSA)*TSNG(NTH)/WCM
         CALL BESSJN(X,NHMAX,ADJ,ADJD)

         DO NC=NCMIN(NS),NCMAX(NS)
            NCD = ABS(NC)
            CDENX = RGMM(NP,NSA)-CKPRW*PM(NP,NSA)*TCSG(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP2(NP,NTH)*DELP(NSA))**2 &
                                   +(DELPL*DGT2(NP,NTH)*DELTH)**2)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CDEN
            CSM12 = CSM12 + PAI1         *CPAI2*CDEN
            CSM13 = CSM13 + PAI1          *PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CDEN
            CSM33 = CSM33 + PAI3          *PAI3*CDEN
         ENDDO

         CPART2= DFT(NP,NTH)*PM(NP,NSA)*PM(NP,NSA) &
                            *TSNG(NTH)*TSNG(NTH) &
                *(TCSG(NTH)-CKPRW*PM(NP,NS)/RGMM(NP,NSA)) &
                *DELTH*DELP(NSA)

         CINTG211= CINTG211 + CSM11*CPART2
         CINTG212= CINTG212 + CSM12*CPART2
         CINTG213= CINTG213 + CSM13*CPART2
         CINTG221= CINTG221 - CSM12*CPART2
         CINTG222= CINTG222 + CSM22*CPART2
         CINTG223= CINTG223 + CSM23*CPART2
         CINTG231= CINTG231 + CSM13*CPART2
         CINTG232= CINTG232 - CSM23*CPART2
         CINTG233= CINTG233 + CSM33*CPART2 &
                            - PM(NP,NSA)*PM(NP,NSA)*TCSG(NTH) &
                              *DFT(NP,NTH)/RGMM(NP,NSA) &
                              *DELTH*DELP(NSA)
      ENDDO
      ENDDO

         FACT=2.D0*PI*DBLE(CWP)

         CLDISP(1)=FACT*(CINTG111+CINTG211)
         CLDISP(2)=FACT*(CINTG133+CINTG233)-CLDISP(1)
         CLDISP(3)=FACT*(CINTG122+CINTG222)-CLDISP(1)
         CLDISP(4)=FACT*(CINTG113+CINTG213)
         CLDISP(5)=FACT*(CINTG112+CINTG212)
         CLDISP(6)=FACT*(CINTG123+CINTG223)
! 
         DEALLOCATE(ADJ,ADJD)
      RETURN
  END SUBROUTINE DP_HOTRR

! ******************************************************
!                       DPHOTRI
! ******************************************************

  SUBROUTINE DP_HOTRIX(CW,CKPR,CKPP,NSA,mag,CLDISP)

    USE dpcomm
    USE plprof
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NSA
    TYPE(pl_mag_type),INTENT(IN):: mag
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NHMAX,NTH,NP,NC,NCD,NP1,INC,NP2,NS
    REAL(rkind):: WCM,DKPRW,PTH0C,PTH0W,DCWC,DNPR,DKPP,D
    REAL(rkind):: PNEAR1,DIF,DFP3,X,RGM,PAI1,PAI3,PNEAR2,DFT4,FACT
    REAL(rkind):: PN0,PT0,PTH0
    COMPLEX(rkind):: CWP,CWC,CKPRW,CPAI2
    COMPLEX(rkind):: CPART31,CPART32,CPART41,CPART42
    COMPLEX(rkind):: CINTG311,CINTG312,CINTG313
    COMPLEX(rkind):: CINTG321,CINTG322,CINTG323
    COMPLEX(rkind):: CINTG331,CINTG332,CINTG333
    COMPLEX(rkind):: CINTG411,CINTG412,CINTG413
    COMPLEX(rkind):: CINTG421,CINTG422,CINTG423
    COMPLEX(rkind):: CINTG431,CINTG432,CINTG433
    COMPLEX(rkind):: CSM11,CSM12,CSM13,CSM22,CSM23,CSM33
    COMPLEX(rkind):: cdelta

!    cdelta=CI*0.01D0
    cdelta=CI*0.003D0
!    cdelta=CI*0.0D0

    NS=NS_NSA_DP(NSA)
      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)),2)+2
      ALLOCATE(ADJ(0:NHMAX),ADJD(0:NHMAX))

      PN0=RNFP0(NSA)
      PT0=RTFP0(NSA)
      PTH0=SQRT(PT0*1.D3*AEE*AMFP(NSA))

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=mag%BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0C=PTH0/(AMP*PA(NS)*VC)
      PTH0W=PTH0C**2
      DCWC=DBLE(CWC)
      DNPR=DBLE(CKPR*VC/CW)
      DKPP=DBLE(CKPP)

      DO NTH=1,NTHMAX_DP
      DO NP=1,NPMAX_DP-1
         DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP(NSA)
      ENDDO
      ENDDO
      DO NP=1,NPMAX_DP
      DO NTH=1,NTHMAX_DP-1
         DFT(NP,NTH) = (FM(NP,NTH+1) - FM(NP,NTH))/DELTH
      ENDDO
      ENDDO

!***************SINGULAR POINT***************************

!
!*****************SUM3************************

      CINTG311 = (0.D0,0.D0)
      CINTG312 = (0.D0,0.D0)
      CINTG313 = (0.D0,0.D0)
      CINTG321 = (0.D0,0.D0)
      CINTG322 = (0.D0,0.D0)
      CINTG323 = (0.D0,0.D0)
      CINTG331 = (0.D0,0.D0)
      CINTG332 = (0.D0,0.D0)
      CINTG333 = (0.D0,0.D0)

      DO NTH=1,NTHMAX_DP

         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         DO NC=NCMIN(NS),NCMAX(NS)

            D = (NC*DCWC)**2+(DNPR*TCSM(NTH))**2-1.D0
            IF(D.LT.0.D0) GOTO 310

!   PNEAR1

            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2) &
                    *(DNPR*NC*DCWC*TCSM(NTH)+SQRT(D))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP(NSA)*NPMAX_DP) GOTO 302

            NP1 = INT(PNEAR1/DELP(NSA))
            IF (NP1.LT.0.OR.NP1.GE.NPMAX_DP) GOTO 302
            IF (NP1.EQ.0) THEN             
               DIF = PNEAR1/DELP(NSA)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP1.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR1 - PG(NP1,NSA))/DELP(NSA)
               DFP3 = (1.D0-DIF)*DFP(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PG(NP1,NSA))/DELP(NSA)
               DFP3 = DIF*DFP(NP1+1,NTH)+(1.D0-DIF)*DFP(NP1,NTH)
            ENDIF

            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNM(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART31= DFP3*PNEAR1**3 &
                     /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSM(NTH)+cdelta)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART31
            CSM12 = CSM12 + PAI1         *CPAI2*CPART31
            CSM13 = CSM13 + PAI1          *PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART31
            CSM33 = CSM33 + PAI3          *PAI3*CPART31

!     PNEAR2

  302       PNEAR2 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2) &
                    *(DNPR*NC*DCWC*TCSM(NTH)-SQRT(D))/PTH0
            IF (PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP(NSA)*NPMAX_DP) GOTO 310

            NP2 = INT(PNEAR2/DELP(NSA))
            IF (NP2.LT.0.OR.NP2.GE.NPMAX_DP) GOTO 310
            IF(NP2.EQ.0) THEN
               DIF = PNEAR2/DELP(NSA)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP2.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR2 - PG(NP2,NSA))/DELP(NSA)
               DFP3 = (1.D0-DIF)*DFP(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PG(NP2,NSA))/DELP(NSA)
               DFP3 = DIF*DFP(NP2+1,NTH)+(1.D0-DIF)*DFP(NP2,NTH)
            ENDIF

            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNM(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART31= DFP3*PNEAR2**3 &
                     /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSM(NTH)+cdelta)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART31
            CSM12 = CSM12 + PAI1         *CPAI2*CPART31
            CSM13 = CSM13 + PAI1          *PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART31
            CSM33 = CSM33 + PAI3          *PAI3*CPART31

  310       CONTINUE
         ENDDO

         CPART32= -CI*PI*TSNM(NTH)**3*DELTH
!              
         CINTG311 = CINTG311 + CSM11*CPART32
         CINTG312 = CINTG312 + CSM12*CPART32
         CINTG313 = CINTG313 + CSM13*CPART32
         CINTG321 = CINTG321 - CSM12*CPART32
         CINTG322 = CINTG322 + CSM22*CPART32
         CINTG323 = CINTG323 + CSM23*CPART32
         CINTG331 = CINTG331 + CSM13*CPART32
         CINTG332 = CINTG332 - CSM23*CPART32
         CINTG333 = CINTG333 + CSM33*CPART32
      ENDDO
! 
!*****************SUM4************************
! 
      CINTG411 = (0.D0,0.D0)
      CINTG412 = (0.D0,0.D0)
      CINTG413 = (0.D0,0.D0)
      CINTG421 = (0.D0,0.D0)
      CINTG422 = (0.D0,0.D0)
      CINTG423 = (0.D0,0.D0)
      CINTG431 = (0.D0,0.D0)
      CINTG432 = (0.D0,0.D0)
      CINTG433 = (0.D0,0.D0)

      DO NTH=1,NTHMAX_DP-1
!               
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
!               
         DO NC=NCMIN(NS),NCMAX(NS)

            D = (NC*DCWC)**2+(DNPR*TCSG(NTH))**2-1.D0
            IF(D.LT.0.D0) GOTO 410

!   PNEAR1

            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2) &
                    *(DNPR*NC*DCWC*TCSG(NTH)+SQRT(D))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP(NSA)*NPMAX_DP) GOTO 402

            NP1 = INT(PNEAR1/DELP(NSA)+0.5D0)
            IF (NP1.LT.0.OR.NP1.GE.NPMAX_DP) GOTO 402
            IF(NP1.EQ.0) THEN
               DIF = (PNEAR1 - PM(1,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP1.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR1 - PM(NP1,NSA))/DELP(NSA)
               DFT4  = (1.D0-DIF)*DFT(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PM(NP1,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(NP1+1,NTH)+(1.D0-DIF)*DFT(NP1,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNG(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART41=DFT4*PNEAR1**2*(TCSG(NTH)-PNEAR1*CKPRW/RGM) &
                    /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSG(NTH)+cdelta)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART41
            CSM12 = CSM12 + PAI1         *CPAI2*CPART41
            CSM13 = CSM13 + PAI1          *PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART41
            CSM33 = CSM33 + PAI3          *PAI3*CPART41

!   PNEAR2

 402        PNEAR2 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2) &
                    *(DNPR*NC*DCWC*TCSG(NTH)-SQRT(D))/PTH0
            IF(PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP(NSA)*NPMAX_DP) GOTO 410
!            IF (DKPRW*TCSG(NTH)*PNEAR2+NC*DCWC.LT.0.D0) GOTO 410

            NP2 = INT(PNEAR2/DELP(NSA)+0.5D0)
            IF (NP2.LT.0.OR.NP2.GE.NPMAX_DP) GOTO 410
            IF(NP2.EQ.0) THEN
               DIF = (PNEAR2 - PM(1,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP2.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR2 - PM(NP2,NSA))/DELP(NSA)
               DFT4  = (1.D0-DIF)*DFT(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PM(NP2,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(NP2+1,NTH)+(1.D0-DIF)*DFT(NP2,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNG(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART41=DFT4*PNEAR2**2*(TCSG(NTH)-PNEAR2*CKPRW/RGM) &
                    /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSG(NTH)+cdelta)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART41
            CSM12 = CSM12 + PAI1         *CPAI2*CPART41
            CSM13 = CSM13 + PAI1          *PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART41
            CSM33 = CSM33 + PAI3          *PAI3*CPART41

  410       CONTINUE 
         ENDDO

         CPART42= -CI*PI*TSNG(NTH)**2*DELTH

         CINTG411 = CINTG411 + CSM11*CPART42
         CINTG412 = CINTG412 + CSM12*CPART42
         CINTG413 = CINTG413 + CSM13*CPART42
         CINTG421 = CINTG421 - CSM12*CPART42
         CINTG422 = CINTG422 + CSM22*CPART42
         CINTG423 = CINTG423 + CSM23*CPART42
         CINTG431 = CINTG431 + CSM13*CPART42
         CINTG432 = CINTG432 - CSM23*CPART42
         CINTG433 = CINTG433 + CSM33*CPART42
      ENDDO

      FACT=2.D0*PI*DBLE(CWP)

      CLDISP(1)=FACT*(CINTG311+CINTG411)
      CLDISP(2)=FACT*(CINTG333+CINTG433)-CLDISP(1)
      CLDISP(3)=FACT*(CINTG322+CINTG422)-CLDISP(1)
      CLDISP(4)=FACT*(CINTG313+CINTG413)
      CLDISP(5)=FACT*(CINTG312+CINTG412)
      CLDISP(6)=FACT*(CINTG323+CINTG423)

      DEALLOCATE(ADJ,ADJD)
      RETURN
  END SUBROUTINE DP_HOTRIX


! ******************************************************
!                       DPHOTRI
! ******************************************************

  SUBROUTINE DP_HOTRI(CW,CKPR,CKPP,NSA,mag,CLDISP)

    USE dpcomm
    USE plprof
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NSA
    TYPE(pl_mag_type),INTENT(IN):: mag
!    TYPE(pl_prf_type),DIMENSION(nsmax),INTENT(IN):: plf
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NHMAX,NTH,NP,NC,NCD,NP1,INC,NP2,NS
    REAL(rkind):: WCM,DKPRW,PTH0W,DCWC,DNPR,DKPP,D
    REAL(rkind):: PNEAR1,DIF,DFP3,X,RGM,PAI1,PAI3,PNEAR2,DFT4,FACT
    REAL(rkind):: PN0,PT0,PTH0
    COMPLEX(rkind):: CWP,CWC,CKPRW,CPAI2
    COMPLEX(rkind):: CPART31,CPART32,CPART41,CPART42
    COMPLEX(rkind):: CINTG311,CINTG312,CINTG313
    COMPLEX(rkind):: CINTG321,CINTG322,CINTG323
    COMPLEX(rkind):: CINTG331,CINTG332,CINTG333
    COMPLEX(rkind):: CINTG411,CINTG412,CINTG413
    COMPLEX(rkind):: CINTG421,CINTG422,CINTG423
    COMPLEX(rkind):: CINTG431,CINTG432,CINTG433
    COMPLEX(rkind):: CSM11,CSM12,CSM13,CSM22,CSM23,CSM33

    NS=NS_NSA_DP(NSA)
    NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)),2)+2
      ALLOCATE(ADJ(0:NHMAX),ADJD(0:NHMAX))

      PN0=RNFP0(NSA)
      PT0=RTFP0(NSA)
      PTH0=SQRT(PT0*1.D3*AEE*AMFP(NSA))

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=mag%BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      DCWC=DBLE(CWC)
      DNPR=DBLE(CKPR*VC/CW)
      DKPP=DBLE(CKPP)

      DO NTH=1,NTHMAX_DP
         DO NP=1,NPMAX_DP-1
            DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP(NSA)
         END DO
      END DO
      DO NP=1,NPMAX_DP
         DO NTH=1,NTHMAX_DP-1
            DFT(NP,NTH) = (FM(NP,NTH+1) - FM(NP,NTH))/DELTH
         END DO
      END DO

!***************SINGULAR POINT***************************

!*****************SUM3************************

      CINTG311 = (0.D0,0.D0)
      CINTG312 = (0.D0,0.D0)
      CINTG313 = (0.D0,0.D0)
      CINTG321 = (0.D0,0.D0)
      CINTG322 = (0.D0,0.D0)
      CINTG323 = (0.D0,0.D0)
      CINTG331 = (0.D0,0.D0)
      CINTG332 = (0.D0,0.D0)
      CINTG333 = (0.D0,0.D0)

      DO 300 NTH=1,NTHMAX_DP

         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         DO 310 NC=NCMIN(NS),NCMAX(NS)

            D = (NC*DCWC)**2+(DNPR*TCSM(NTH))**2
            IF(D.LE.1.D0) GOTO 310

!   PNEAR1

            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2)* &
                    (DNPR*NC*DCWC*TCSM(NTH)+SQRT(D-1))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP(NSA)*NPMAX_DP) THEN
               GOTO 302
            ELSEIF (DKPRW*TCSM(NTH)*PNEAR1+NC*DCWC.LT.0.D0) THEN
               GOTO 302
            END IF
            NP1 = INT(PNEAR1/DELP(NSA))
            IF (NP1.LT.0.OR.NP1.GE.NPMAX_DP) THEN
               GOTO 302
            ELSEIF (NP1.EQ.0) THEN
               DIF = PNEAR1/DELP(NSA)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP1.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR1 - PG(NP1,NSA))/DELP(NSA)
               DFP3 = (1.D0-DIF)*DFP(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PG(NP1,NSA))/DELP(NSA)
               DFP3 = DIF*DFP(NP1+1,NTH)+(1.D0-DIF)*DFP(NP1,NTH)
            ENDIF

            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNM(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART31= DFP3*PNEAR1**3 &
                    /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSM(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART31
            CSM12 = CSM12 + PAI1         *CPAI2*CPART31
            CSM13 = CSM13 + PAI1          *PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART31
            CSM33 = CSM33 + PAI3          *PAI3*CPART31

!     PNEAR2

  302       PNEAR2 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2)* &
                    (DNPR*NC*DCWC*TCSM(NTH)-SQRT(D-1))/PTH0
            IF (PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP(NSA)*NPMAX_DP) THEN 
               GOTO 310
            ELSEIF (DKPRW*TCSM(NTH)*PNEAR2+NC*DCWC.LT.0.D0) THEN
               GOTO 310
            END IF
            NP2 = INT(PNEAR2/DELP(NSA))
            IF (NP2.LT.0.OR.NP2.GE.NPMAX_DP) THEN
               GOTO 310
            ELSE IF(NP2.EQ.0) THEN
               DIF = PNEAR2/DELP(NSA)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP2.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR2 - PG(NP2,NSA))/DELP(NSA)
               DFP3 = (1.D0-DIF)*DFP(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PG(NP2,NSA))/DELP(NSA)
               DFP3 = DIF*DFP(NP2+1,NTH)+(1.D0-DIF)*DFP(NP2,NTH)
            ENDIF

            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNM(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART31= DFP3*PNEAR2**3 &
                    /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSM(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART31
            CSM12 = CSM12 + PAI1         *CPAI2*CPART31
            CSM13 = CSM13 + PAI1          *PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART31
            CSM33 = CSM33 + PAI3          *PAI3*CPART31
  310    CONTINUE 

         CPART32= -CI*PI*TSNM(NTH)**3*DELTH

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
 
!*****************SUM4************************
 
      CINTG411 = (0.D0,0.D0)
      CINTG412 = (0.D0,0.D0)
      CINTG413 = (0.D0,0.D0)
      CINTG421 = (0.D0,0.D0)
      CINTG422 = (0.D0,0.D0)
      CINTG423 = (0.D0,0.D0)
      CINTG431 = (0.D0,0.D0)
      CINTG432 = (0.D0,0.D0)
      CINTG433 = (0.D0,0.D0)

      DO 400 NTH=1,NTHMAX_DP-1

         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         DO 410 NC=NCMIN(NS),NCMAX(NS)

            D = (NC*DCWC)**2+(DNPR*TCSG(NTH))**2
            IF(D.LE.1.D0) GOTO 410

!   PNEAR1

            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2)* &
                     (DNPR*NC*DCWC*TCSG(NTH)+SQRT(D-1))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP(NSA)*NPMAX_DP) THEN
               GOTO 402
            ELSEIF (DKPRW*TCSG(NTH)*PNEAR1+NC*DCWC.LT.0.D0) THEN
               GOTO 402
            END IF
            NP1 = INT(PNEAR1/DELP(NSA)+0.5D0)
            IF (NP1.LT.0.OR.NP1.GE.NPMAX_DP) THEN
               GOTO 402
            ELSE IF(NP1.EQ.0) THEN
               DIF = (PNEAR1 - PM(1,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP1.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR1 - PM(NP1,NSA))/DELP(NSA)
               DFT4  = (1.D0-DIF)*DFT(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PM(NP1,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(NP1+1,NTH)+(1.D0-DIF)*DFT(NP1,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNG(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART41=DFT4*PNEAR1**2*(TCSG(NTH)-PNEAR1*CKPRW/RGM) &
                   /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSG(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART41
            CSM12 = CSM12 + PAI1         *CPAI2*CPART41
            CSM13 = CSM13 + PAI1          *PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART41
            CSM33 = CSM33 + PAI3          *PAI3*CPART41

!   PNEAR2

 402        PNEAR2 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2)* &
                    (DNPR*NC*DCWC*TCSG(NTH)-SQRT(D-1))/PTH0
            IF(PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP(NSA)*NPMAX_DP) THEN
               GOTO 410
            ELSEIF (DKPRW*TCSG(NTH)*PNEAR2+NC*DCWC.LT.0.D0) THEN
               GOTO 410
            END IF
            NP2 = INT(PNEAR2/DELP(NSA)+0.5D0)
            IF (NP2.LT.0.OR.NP2.GE.NPMAX_DP) THEN
               GOTO 410
            ELSE IF(NP2.EQ.0) THEN
               DIF = (PNEAR2 - PM(1,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP2.EQ.NPMAX_DP-1) THEN
               DIF = (PNEAR2 - PM(NP2,NSA))/DELP(NSA)
               DFT4  = (1.D0-DIF)*DFT(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PM(NP2,NSA))/DELP(NSA)
               DFT4  = DIF*DFT(NP2+1,NTH)+(1.D0-DIF)*DFT(NP2,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNG(NTH)/WCM
            CALL BESSJN(X,NHMAX,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART41=DFT4*PNEAR2**2*(TCSG(NTH)-PNEAR2*CKPRW/RGM) &
                   /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSG(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               IF(NCD.EQ.1) THEN
                  PAI1=NC*INC*0.5D0
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*INC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CPART41
            CSM12 = CSM12 + PAI1         *CPAI2*CPART41
            CSM13 = CSM13 + PAI1          *PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CPART41
            CSM33 = CSM33 + PAI3          *PAI3*CPART41

  410    CONTINUE 

         CPART42= -CI*PI*TSNG(NTH)**2*DELTH

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

      FACT=2.D0*PI*DBLE(CWP)

      CLDISP(1)=FACT*(CINTG311+CINTG411)
      CLDISP(2)=FACT*(CINTG333+CINTG433)-CLDISP(1)
      CLDISP(3)=FACT*(CINTG322+CINTG422)-CLDISP(1)
      CLDISP(4)=FACT*(CINTG313+CINTG413)
      CLDISP(5)=FACT*(CINTG312+CINTG412)
      CLDISP(6)=FACT*(CINTG323+CINTG423)

      RETURN
    END  SUBROUTINE DP_HOTRI
END MODULE dphotr
