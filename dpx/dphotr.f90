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

  SUBROUTINE DP_HOTR(CW,CKPR,CKPP,NS,CLDISP)

    USE dpcomm
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    COMPLEX(rkind):: CLDISP1(6),CLDISP2(6)
    INTEGER:: I
      
    CALL DP_HOTRR(CW,CKPR,CKPP,NS,CLDISP1)
    CALL DP_HOTRI(CW,CKPR,CKPP,NS,CLDISP2)
    DO I=1,6
       CLDISP(I)=CLDISP1(I)+CLDISP2(I)
    ENDDO
    RETURN
  END SUBROUTINE DP_HOTR

! ******************************************************
!                       DPHOTRR
! ******************************************************

  SUBROUTINE DP_HOTRR(CW,CKPR,CKPP,NS,CLDISP)

    USE dpcomm
    USE pllocal
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NCMIN,NCMAX,NHMAX,NTH,NP,NC,NCD,INC
    REAL(rkind):: DELPL,WCM,DKPRW,PTH0W,DKPP,X,PAI1,PAI3,PART1,FACT
    COMPLEX(rkind):: CWP,CWC,CKPRW,CDENX,CDEN,CPAI2,CPART2
    COMPLEX(rkind):: CINTG111,CINTG112,CINTG113
    COMPLEX(rkind):: CINTG121,CINTG122,CINTG123
    COMPLEX(rkind):: CINTG131,CINTG132,CINTG133
    COMPLEX(rkind):: CINTG211,CINTG212,CINTG213
    COMPLEX(rkind):: CINTG221,CINTG222,CINTG223
    COMPLEX(rkind):: CINTG231,CINTG232,CINTG233
    COMPLEX(rkind):: CSM11,CSM12,CSM13,CSM22,CSM23,CSM33

      NCMIN = NDISP1(NS)
      NCMAX = NDISP2(NS)
      NHMAX=MAX(ABS(NCMIN),ABS(NCMAX),2)+5
      ALLOCATE(ADJ(0:NHMAX),ADJD(0:NHMAX))
      DELPL  = 0.5D0

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0W=(PTH0/(AMP*PA(NS)*VC))**2
      DKPP=DBLE(CKPP)

      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP(NS)
      ENDDO
      ENDDO
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX-1
         DFT(NP,NTH) = (FM(NP,NTH+1) - FM(NP,NTH))/DELTH
      ENDDO
      ENDDO

!***********DGP1,DGP2,DGT1,DGT2************

      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         DGP1(NP,NTH)=PTH0W*PG(NP,NS)/SQRT(1+PTH0W*PG(NP,NS)**2) &
                     -DKPRW*TCSM(NTH)
         DGP2(NP,NTH)=PTH0W*PM(NP,NS)/SQRT(1+PTH0W*PM(NP,NS)**2) &
                     -DKPRW*TCSG(NTH)
         DGT1(NP,NTH)=DKPRW*PG(NP,NS)*TSNM(NTH)
         DGT2(NP,NTH)=DKPRW*PM(NP,NS)*TSNG(NTH)
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

      DO NP=1,NPMAX-1
      DO NTH=1,NTHMAX
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         X = DKPP*PTH0*PG(NP,NS)*TSNM(NTH)/WCM
         CALL BESSJN(X,NHMAX,ADJ,ADJD)

         DO NC=NCMIN,NCMAX
            NCD = ABS(NC)
            CDENX= RGMG(NP)-CKPRW*PG(NP,NS)*TCSM(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP1(NP,NTH)*DELP(NS))**2 &
                                   +(DELPL*DGT1(NP,NTH)*DELTH)**2)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               if(NCD.eq.1) THEN
                  PAI1=NC*INC*0.5D0
               else
                  PAI1=0.D0
               endif
            else
               PAI1  = NC*INC*ADJ(NCD)/X
            endif
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CDEN
            CSM12 = CSM12 + PAI1         *CPAI2*CDEN
            CSM13 = CSM13 + PAI1          *PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CDEN
            CSM33 = CSM33 + PAI3          *PAI3*CDEN
         ENDDO

         PART1= DFP(NP,NTH)*PG(NP,NS)*PG(NP,NS)*PG(NP,NS) &
                           *TSNM(NTH)*TSNM(NTH)*TSNM(NTH) &
               *DELTH*DELP(NS)

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

      DO NP=1,NPMAX
      DO NTH=1,NTHMAX-1
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         X = DKPP*PTH0*PM(NP,NS)*TSNG(NTH)/WCM
         CALL BESSJN(X,NHMAX,ADJ,ADJD)

         DO NC=NCMIN,NCMAX
            NCD = ABS(NC)
            CDENX = RGMM(NP)-CKPRW*PM(NP,NS)*TCSG(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP2(NP,NTH)*DELP(NS))**2 &
                                   +(DELPL*DGT2(NP,NTH)*DELTH)**2)

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               if(NCD.eq.1) THEN
                  PAI1=NC*INC*0.5D0
               else
                  PAI1=0.D0
               endif
            else
               PAI1  = NC*INC*ADJ(NCD)/X
            endif
            CPAI2 = CI*INC*ADJD(NCD)
            PAI3  =    INC*ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1          *PAI1*CDEN
            CSM12 = CSM12 + PAI1         *CPAI2*CDEN
            CSM13 = CSM13 + PAI1          *PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2) *PAI3*CDEN
            CSM33 = CSM33 + PAI3          *PAI3*CDEN
         ENDDO

         CPART2= DFT(NP,NTH)*PM(NP,NS)*PM(NP,NS) &
                            *TSNG(NTH)*TSNG(NTH) &
                *(TCSG(NTH)-CKPRW*PM(NP,NS)/RGMM(NP)) &
                *DELTH*DELP(NS)

         CINTG211= CINTG211 + CSM11*CPART2
         CINTG212= CINTG212 + CSM12*CPART2
         CINTG213= CINTG213 + CSM13*CPART2
         CINTG221= CINTG221 - CSM12*CPART2
         CINTG222= CINTG222 + CSM22*CPART2
         CINTG223= CINTG223 + CSM23*CPART2
         CINTG231= CINTG231 + CSM13*CPART2
         CINTG232= CINTG232 - CSM23*CPART2
         CINTG233= CINTG233 + CSM33*CPART2 &
                            - PM(NP,NS)*PM(NP,NS)*TCSG(NTH) &
                              *DFT(NP,NTH)/RGMM(NP) &
                              *DELTH*DELP(NS)
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

  SUBROUTINE DP_HOTRI(CW,CKPR,CKPP,NS,CLDISP)

    USE dpcomm
    USE pllocal
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NCMIN,NCMAX,NHMAX,NTH,NP,NC,NCD,NP1,INC,NP2
    REAL(rkind):: WCM,DKPRW,PTH0C,PTH0W,DCWC,DNPR,DKPP,D
    REAL(rkind):: PNEAR1,DIF,DFP3,X,RGM,PAI1,PAI3,PNEAR2,DFT4,FACT
    COMPLEX(rkind):: CWP,CWC,CKPRW,CPAI2
    COMPLEX(rkind):: CPART31,CPART32,CPART41,CPART42
    COMPLEX(rkind):: CINTG311,CINTG312,CINTG313
    COMPLEX(rkind):: CINTG321,CINTG322,CINTG323
    COMPLEX(rkind):: CINTG331,CINTG332,CINTG333
    COMPLEX(rkind):: CINTG411,CINTG412,CINTG413
    COMPLEX(rkind):: CINTG421,CINTG422,CINTG423
    COMPLEX(rkind):: CINTG431,CINTG432,CINTG433
    COMPLEX(rkind):: CSM11,CSM12,CSM13,CSM22,CSM23,CSM33

      NCMIN = NDISP1(NS)
      NCMAX = NDISP2(NS)
      NHMAX=MAX(ABS(NCMIN),ABS(NCMAX),2)+5
      ALLOCATE(ADJ(0:NHMAX),ADJD(0:NHMAX))

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
      PTH0C=PTH0/(AMP*PA(NS)*VC)
      PTH0W=PTH0C**2
      DCWC=DBLE(CWC)
      DNPR=DBLE(CKPR*VC/CW)
      DKPP=DBLE(CKPP)

      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DFP(NP,NTH) = (FM(NP+1,NTH) - FM(NP,NTH))/DELP(NS)
      ENDDO
      ENDDO
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX-1
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

      DO NTH=1,NTHMAX

         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)

         DO NC=NCMIN,NCMAX

            D = (NC*DCWC)**2+(DNPR*TCSM(NTH))**2-1
            IF(D.LT.0.D0) GOTO 310

!   PNEAR1

            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSM(NTH))**2) &
                    *(DNPR*NC*DCWC*TCSM(NTH)+SQRT(D))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP(NS)*NPMAX) GOTO 302

!            IF(DNPR**2.LE.1.D0) THEN
!               FL = NC*DCWC/SQRT(1-DNPR**2)
!               IF(FL.LT.1.D0) THEN
!                  WRITE(6,'(A,I5,1p3E12.4)') 
!     &                 'FL:NC,DCWC,DNPR,FL=',NC,DCWC,DNPR,FL
!               ENDIF
!            ENDIF

            NP1 = INT(PNEAR1/DELP(NS))
            IF (NP1.LT.0.OR.NP1.GE.NPMAX) GOTO 302
            IF (NP1.EQ.0) THEN
               DIF = PNEAR1/DELP(NS)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP1.EQ.NPMAX-1) THEN
               DIF = (PNEAR1 - PG(NP1,NS))/DELP(NS)
               DFP3 = (1.D0-DIF)*DFP(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PG(NP1,NS))/DELP(NS)
               DFP3 = DIF*DFP(NP1+1,NTH)+(1.D0-DIF)*DFP(NP1,NTH)
            ENDIF

            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNM(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART31= DFP3*PNEAR1**3 &
                     /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSM(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               if(NCD.eq.1) THEN
                  PAI1=NC*INC*0.5D0
               else
                  PAI1=0.D0
               endif
            else
               PAI1  = NC*INC*ADJ(NCD)/X
            endif
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
            IF (PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP(NS)*NPMAX) GOTO 310

!            IF(DNPR**2.LE.1.D0) THEN
!               FL = NC*DCWC/SQRT(1-DNPR**2)
!               IF(FL.LT.1.D0) THEN
!                  WRITE(6,'(A,I5,1p3E12.4)') 
!     &                 'FL:NC,DCWC,DNPR,FL=',NC,DCWC,DNPR,FL
!               ENDIF
!            ENDIF

            NP2 = INT(PNEAR2/DELP(NS))
            IF (NP2.LT.0.OR.NP2.GE.NPMAX) GOTO 310
            IF(NP2.EQ.0) THEN
               DIF = PNEAR2/DELP(NS)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP2.EQ.NPMAX-1) THEN
               DIF = (PNEAR2 - PG(NP2,NS))/DELP(NS)
               DFP3 = (1.D0-DIF)*DFP(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PG(NP2,NS))/DELP(NS)
               DFP3 = DIF*DFP(NP2+1,NTH)+(1.D0-DIF)*DFP(NP2,NTH)
            ENDIF

            NCD = ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNM(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART31= DFP3*PNEAR2**3 &
                     /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSM(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               if(NCD.eq.1) THEN
                  PAI1=NC*INC*0.5D0
               else
                  PAI1=0.D0
               endif
            else
               PAI1  = NC*INC*ADJ(NCD)/X
            endif
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

      DO NTH=1,NTHMAX-1
!               
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
!               
         DO NC=NCMIN,NCMAX

            D = (NC*DCWC)**2+(DNPR*TCSG(NTH))**2-1.D0
            IF(D.LT.0.D0) GOTO 410

!   PNEAR1

            PNEAR1 = (AMP*PA(NS)*VC)/(1-(DNPR*TCSG(NTH))**2) &
                    *(DNPR*NC*DCWC*TCSG(NTH)+SQRT(D))/PTH0
            IF (PNEAR1.LT.0.D0.OR.PNEAR1.GT.DELP(NS)*NPMAX) GOTO 402

!            IF(DNPR**2.LE.1.D0) THEN
!               FL = NC*DCWC/SQRT(1-DNPR**2)
!               IF(FL.LT.1.D0) THEN
!                  WRITE(6,'(A,I5,1p3E12.4)') 
!     &                 'FL:NC,DCWC,DNPR,FL=',NC,DCWC,DNPR,FL
!               ENDIF
!            ENDIF

            NP1 = INT(PNEAR1/DELP(NS)+0.5D0)
            IF (NP1.LT.0.OR.NP1.GE.NPMAX) GOTO 402
            IF(NP1.EQ.0) THEN
               DIF = (PNEAR1 - PM(1,NS))/DELP(NS)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP1.EQ.NPMAX-1) THEN
               DIF = (PNEAR1 - PM(NP1,NS))/DELP(NS)
               DFT4  = (1.D0-DIF)*DFT(NP1,NTH)
            ELSE
               DIF = (PNEAR1 - PM(NP1,NS))/DELP(NS)
               DFT4  = DIF*DFT(NP1+1,NTH)+(1.D0-DIF)*DFT(NP1,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR1*TSNG(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR1**2) 
            CPART41=DFT4*PNEAR1**2*(TCSG(NTH)-PNEAR1*CKPRW/RGM) &
                    /ABS(PTH0W*PNEAR1/RGM-CKPRW*TCSG(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               if(NCD.eq.1) THEN
                  PAI1=NC*INC*0.5D0
               else
                  PAI1=0.D0
               endif
            else
               PAI1  = NC*INC*ADJ(NCD)/X
            endif
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
            IF(PNEAR2.LT.0.D0.OR.PNEAR2.GT.DELP(NS)*NPMAX) GOTO 410

!            IF(DNPR**2.LE.1.D0) THEN
!               FL = NC*DCWC/SQRT(1-DNPR**2)
!               IF(FL.LT.1.D0) THEN
!                  WRITE(6,'(A,I5,1p3E12.4)') 
!     &                 'FL:NC,DCWC,DNPR,FL=',NC,DCWC,DNPR,FL
!               ENDIF
!            ENDIF

            NP2 = INT(PNEAR2/DELP(NS)+0.5D0)
            IF (NP2.LT.0.OR.NP2.GE.NPMAX) GOTO 410
            IF(NP2.EQ.0) THEN
               DIF = (PNEAR2 - PM(1,NS))/DELP(NS)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSE IF(NP2.EQ.NPMAX-1) THEN
               DIF = (PNEAR2 - PM(NP2,NS))/DELP(NS)
               DFT4  = (1.D0-DIF)*DFT(NP2,NTH)
            ELSE
               DIF = (PNEAR2 - PM(NP2,NS))/DELP(NS)
               DFT4  = DIF*DFT(NP2+1,NTH)+(1.D0-DIF)*DFT(NP2,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR2*TSNG(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)

            RGM=SQRT(1+PTH0W*PNEAR2**2) 
            CPART41=DFT4*PNEAR2**2*(TCSG(NTH)-PNEAR2*CKPRW/RGM) &
                    /ABS(PTH0W*PNEAR2/RGM-CKPRW*TCSG(NTH))

            IF(NC.LT.0.AND.MOD(-NC,2).EQ.1) THEN
               INC=-1
            ELSE
               INC=1
            ENDIF
            IF(X.EQ.0.d0) THEN
               if(NCD.eq.1) THEN
                  PAI1=NC*INC*0.5D0
               else
                  PAI1=0.D0
               endif
            else
               PAI1  = NC*INC*ADJ(NCD)/X
            endif
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

!      WRITE(6,'(1P6E12.4)') CLDISP(1),CLDISP(2),CLDISP(3)
!      WRITE(6,'(1P6E12.4)') CLDISP(4),CLDISP(5),CLDISP(6)

      DEALLOCATE(ADJ,ADJD)
      RETURN
  END SUBROUTINE DP_HOTRI
END MODULE dphotr
