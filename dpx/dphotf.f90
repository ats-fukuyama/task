!
! ******************************************************
!                       DPHOTF              
!      dielectric tensor without relativistic effects
!                 
!                      94/05/14             
!                           programed by K.TANAKA
! ******************************************************

MODULE dphotf

  PRIVATE
  PUBLIC DP_HOTF,DP_HOTFR,DP_HOTFI

CONTAINS

  SUBROUTINE DP_HOTF(CW,CKPR,CKPP,NS,mag,CLDISP)

    USE dpcomm
    USE plprof
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    TYPE(pl_mag_type),INTENT(IN):: mag
!    TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    COMPLEX(rkind):: CLDISP1(6),CLDISP2(6)
    INTEGER:: I
      
    CALL DP_HOTFR(CW,CKPR,CKPP,NS,mag,CLDISP1)
    CALL DP_HOTFI(CW,CKPR,CKPP,NS,mag,CLDISP2)
    DO I=1,6
       CLDISP(I)=CLDISP1(I)+CLDISP2(I)
    ENDDO
    RETURN
  END SUBROUTINE DP_HOTF

! ******************************************************
!                       DPHOTFR
! ******************************************************

  SUBROUTINE DP_HOTFR(CW,CKPR,CKPP,NS,mag,CLDISP)

    USE dpcomm
    USE plprof
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    TYPE(pl_mag_type),INTENT(IN):: mag
!    TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NCMIN,NCMAX,NHMAX,NTH,NP,NC,NCD
    REAL(rkind):: RGM,WCM,DKPRW,DKPP,X,PAI1,PAI3,PART1,FACT,DELPL
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
      RGM   = 1.D0
      DELPL = 0.5D0

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=mag%BABS*PZ(NS)*AEE
      CKPRW= CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
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
         DGP1(NP,NTH)=-DKPRW*TCSM(NTH)
         DGP2(NP,NTH)=-DKPRW*TCSG(NTH)
         DGT1(NP,NTH)= DKPRW*PG(NP,NS)*TSNM(NTH)
         DGT2(NP,NTH)= DKPRW*PM(NP,NS)*TSNG(NTH)
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
            CDENX= RGM-CKPRW*PG(NP,NS)*TCSM(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP1(NP,NTH)*DELP(NS))**2 &
                                   +(DELPL*DGT1(NP,NTH)*DELTH)**2)
            IF(X.EQ.0.D0) THEN
               IF(NCD.EQ.0) THEN
                  PAI1=0.D0
               ELSEIF(NCD.EQ.1) THEN
                  PAI1=0.5D0*NC
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*ADJD(NCD)
            PAI3  = ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1* PAI1*CDEN
            CSM12 = CSM12 + PAI1*CPAI2*CDEN
            CSM13 = CSM13 + PAI1* PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2)* PAI3*CDEN
            CSM33 = CSM33 + PAI3* PAI3*CDEN
         ENDDO

         PART1= DFP(NP,NTH)*PG(NP,NS)*PG(NP,NS)*PG(NP,NS) &
                           *TSNM(NTH)*TSNM(NTH)*TSNM(NTH) &
               *DELTH*DELP(NS)
!     
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
            CDENX = RGM-CKPRW*PM(NP,NS)*TCSG(NTH)-NC*CWC
            CDEN  = CDENX/(CDENX**2+(DELPL*DGP2(NP,NTH)*DELP(NS))**2 &
                                   +(DELPL*DGT2(NP,NTH)*DELTH)**2)
            IF(X.EQ.0.D0) THEN
               IF(NCD.EQ.0) THEN
                  PAI1=0.D0
               ELSEIF(NCD.EQ.1) THEN
                  PAI1=0.5D0*NC
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*ADJD(NCD)
            PAI3  = ADJ(NCD)/TTNG(NTH)
            CSM11 = CSM11 + PAI1* PAI1*CDEN
            CSM12 = CSM12 + PAI1*CPAI2*CDEN
            CSM13 = CSM13 + PAI1* PAI3*CDEN
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CDEN
            CSM23 = CSM23 + DCONJG(CPAI2)* PAI3*CDEN
            CSM33 = CSM33 + PAI3* PAI3*CDEN
         ENDDO
! 
         CPART2= DFT(NP,NTH)*PM(NP,NS)*PM(NP,NS) &
                            *TSNG(NTH)*TSNG(NTH) &
                *(TCSG(NTH)-CKPRW*PM(NP,NS)/RGM) &
               *DELTH*DELP(NS)
! 
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
                              *DFT(NP,NTH)/RGM &
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

         DEALLOCATE(ADJ,ADJD)
      RETURN
  END SUBROUTINE DP_HOTFR

! ******************************************************
!                       DPHOTFI
! ******************************************************

  SUBROUTINE DP_HOTFI(CW,CKPR,CKPP,NS,mag,CLDISP)

    USE dpcomm
    USE plprof
    USE libbes,ONLY: bessjn
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    TYPE(pl_mag_type),INTENT(IN):: mag
!    TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: ADJ,ADJD
    INTEGER:: NCMIN,NCMAX,NHMAX,NTH,NP,NC,NCD
    REAL(rkind):: RGM,WCM,DKPRW,DKPP,PNEAR,DIF,DFP3,X,PAI1,PAI3,DFT4,FACT
    COMPLEX(rkind):: CWP,CWC,CKPRW,CPAI2
    COMPLEX(rkind):: CPART31,CPART32,CPART41,CPART411,CPART42
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
      RGM   = 1.D0

      CWP=PN0*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      WCM=mag%BABS*PZ(NS)*AEE
      CKPRW=CKPR*PTH0/(AMP*PA(NS)*CW)
      DKPRW=DBLE(CKPRW)
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
         DGP1(NP,NTH)=-DKPRW*TCSM(NTH)
         DGP2(NP,NTH)=-DKPRW*TCSG(NTH)
         DGT1(NP,NTH)= DKPRW*PG(NP,NS)*TSNM(NTH)
         DGT2(NP,NTH)= DKPRW*PM(NP,NS)*TSNG(NTH)
      ENDDO
      ENDDO

!
!***************SINGULAR POINT***************************
!  
! 
!*****************SUM3************************
! 
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
!               
         CSM11 = (0.D0,0.D0)
         CSM12 = (0.D0,0.D0)
         CSM13 = (0.D0,0.D0)
         CSM22 = (0.D0,0.D0)
         CSM23 = (0.D0,0.D0)
         CSM33 = (0.D0,0.D0)
!               
         DO NC=NCMIN,NCMAX

            PNEAR = DBLE((RGM-CWC*NC)/(CKPRW*TCSM(NTH)))
            IF(PNEAR.LT.0.D0.OR.PNEAR.GT.DELP(NS)*NPMAX) GOTO 310
            NP = INT(PNEAR/DELP(NS))
            IF (NP.LT.0.OR.NP.GE.NPMAX) GOTO 310
            IF (NP.EQ.0) THEN
               DIF = PNEAR/DELP(NS)
               DFP3 = DIF*DFP(1,NTH)
            ELSE IF(NP.EQ.NPMAX-1) THEN
               DIF = (PNEAR - PG(NP,NS))/DELP(NS)
               DFP3 = (1.D0-DIF)*DFP(NP,NTH)
            ELSE
               DIF = (PNEAR - PG(NP,NS))/DELP(NS)
               DFP3 = DIF*DFP(NP+1,NTH)+(1.D0-DIF)*DFP(NP,NTH)
            ENDIF

            NCD = ABS(NC)

            X = DKPP*PTH0*PNEAR*TSNM(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)

            CPART31= DFP3*(RGM-NC*CWC)**3

            IF(X.EQ.0.D0) THEN
               IF(NCD.EQ.0) THEN
                  PAI1=0.D0
               ELSEIF(NCD.EQ.1) THEN
                  PAI1=0.5D0*NC
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*ADJD(NCD)
            PAI3  = ADJ(NCD)/TTNM(NTH)

            CSM11 = CSM11 + PAI1* PAI1*CPART31
            CSM12 = CSM12 + PAI1*CPAI2*CPART31
            CSM13 = CSM13 + PAI1* PAI3*CPART31
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART31
            CSM23 = CSM23 + DCONJG(CPAI2)* PAI3*CPART31
            CSM33 = CSM33 + PAI3* PAI3*CPART31
  310    CONTINUE 
         ENDDO
         CPART32= -CI*PI*ABS(1.D0/(CKPRW*TCSM(NTH))) &
                 *(TTNM(NTH)/CKPRW)**3*DELTH
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
!               
            IF(NTH*2.EQ.NTHMAX) GOTO 410
            PNEAR = DBLE((RGM-CWC*NC)/(CKPRW*TCSG(NTH)))
            IF(PNEAR.LT.0.D0.OR.PNEAR.GT.DELP(NS)*NPMAX) GOTO 410
            NP = INT(PNEAR/DELP(NS)+0.5D0)
            IF(NP.LT.0.OR.NP.GE.NPMAX) GOTO 410
            IF (NP.EQ.0) THEN
               DIF = (PNEAR - PM(1,NS))/DELP(NS)
               DFT4  = DIF*DFT(1,NTH)-(1.D0-DIF)*DFT(1,NTH)
            ELSEIF(NP.EQ.NPMAX-1) THEN
               DIF = (PNEAR - PM(NP,NS))/DELP(NS)
               DFT4  = (1.D0-DIF)*DFT(NP,NTH)
            ELSE
               DIF = (PNEAR - PM(NP,NS))/DELP(NS)
               DFT4  = DIF*DFT(NP+1,NTH)+(1.D0-DIF)*DFT(NP,NTH)
            ENDIF

            NCD=ABS(NC)
            X = DKPP*PTH0*PNEAR*TSNG(NTH)/WCM
            CALL BESSJN(X,MAX(NCD,2)+5,ADJ,ADJD)

            CPART411=RGM-CWC*NC
            CPART41= DFT4*CPART411*CPART411 &
                    *(TCSG(NTH)-CPART411/(RGM*TCSG(NTH)))
            IF(X.EQ.0.D0) THEN
               IF(NCD.EQ.0) THEN
                  PAI1=0.D0
               ELSEIF(NCD.EQ.1) THEN
                  PAI1=0.5D0*NC
               ELSE
                  PAI1=0.D0
               ENDIF
            ELSE
               PAI1  = NC*ADJ(NCD)/X
            ENDIF
            CPAI2 = CI*ADJD(NCD)
            PAI3  = ADJ(NCD)/TTNG(NTH)

            CSM11 = CSM11 + PAI1*PAI1*CPART41
            CSM12 = CSM12 + PAI1*CPAI2*CPART41
            CSM13 = CSM13 + PAI1*PAI3*CPART41
            CSM22 = CSM22 + DCONJG(CPAI2)*CPAI2*CPART41
            CSM23 = CSM23 + DCONJG(CPAI2)*PAI3*CPART41
            CSM33 = CSM33 + PAI3*PAI3*CPART41
  410    CONTINUE 
         ENDDO

         CPART42= -CI*PI*ABS(1.D0/(CKPRW*TCSG(NTH))) &
                 *(TTNG(NTH)/CKPRW)**2*DELTH

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
  END SUBROUTINE DP_HOTFI
END MODULE dphotf
