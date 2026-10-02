C     $Id: wmpout.f,v 1.33 2013/01/22 16:21:46 fukuyama Exp $
C
C     ****** CALCULATE ELECTRIC FIELD ******
C
      SUBROUTINE WMEFLD
C
      USE libinv
      INCLUDE 'wmcomm.inc'
      DIMENSION CEF1(MDM,NDM),CEF2(MDM,NDM),RMA(3,3)
      DIMENSION CEFM1(MDMM,NDMM),CEFM2(MDMM,NDMM)
      DIMENSION CEFLDR(MDM,NDM), CRMARHF(MDM,NDM)
      DIMENSION CEFLDM(3,MDM,NDMM,NRM)
C
      DRHO1=(XRHO(2)-XRHO(1))**2
      DRHO2=(XRHO(3)-XRHO(1))**2
      A1= DRHO2/(DRHO2-DRHO1)
      A2=-DRHO1/(DRHO2-DRHO1)
C
      DO NR=1,NRMAX+2
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
         IG=3*MDSIZ*NDSIZ*(NR-1)+3*MDSIZ*(NDX-1)+3*(MDX-1)
         CEFLDK(1,MDX,NDX,NR+1)=CFVG(IG+1)
         CEFLDK(2,MDX,NDX,NR+1)=CFVG(IG+2)
         CEFLDK(3,MDX,NDX,NR+1)=CFVG(IG+3)
      ENDDO
      ENDDO
      ENDDO

      NR=1
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
         MM=NTH0+MD
         IF(ABS(MM).EQ.1) THEN
            CEFLDK(1,MDX,NDX,NR)=CEFLDK(1,MDX,NDX,NR+1)
            CEFLDK(2,MDX,NDX,NR)=A1*CEFLDK(2,MDX,NDX,NR+1)
     &                          +A2*CEFLDK(2,MDX,NDX,NR+2)
         ELSE
            CEFLDK(1,MDX,NDX,NR)=0.D0
            CEFLDK(2,MDX,NDX,NR)=0.D0
         ENDIF
         IF(MM.EQ.0) THEN
            CEFLDK(3,MDX,NDX,NR)=A1*CEFLDK(3,MDX,NDX,NR+1)
     &                          +A2*CEFLDK(3,MDX,NDX,NR+2)
         ELSE
            CEFLDK(3,MDX,NDX,NR)=0.D0
         ENDIF
      ENDDO
      ENDDO
C
      CEFLDR=0d0
      DO NR=1,NRMAX+3
         NRP=MIN(NR+1,NRMAX+3)
         IF (NR < NR_S)THEN 
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               CEFLD(1,MDX,NDX,NR)=0.5D0*(CEFLDK(1,MDX,NDX,NR )
     &                                   +CEFLDK(1,MDX,NDX,NRP))
               CEFLD(2,MDX,NDX,NR)=CEFLDK(2,MDX,NDX,NR)
               CEFLD(3,MDX,NDX,NR)=CEFLDK(3,MDX,NDX,NR)
            ENDDO
            ENDDO
         ELSEIF (NR > NR_S + 1)THEN
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               CEFLD(1,MDX,NDX,NR-2)=0.5D0*(CEFLDK(1,MDX,NDX,NR )
     &                                     +CEFLDK(1,MDX,NDX,NRP))
               CEFLD(2,MDX,NDX,NR-2)=CEFLDK(2,MDX,NDX,NR)
               CEFLD(3,MDX,NDX,NR-2)=CEFLDK(3,MDX,NDX,NR)
            ENDDO
            ENDDO
         ENDIF
      ENDDO
C
      NR=1
      DO ND=NDMIN,NDMAX
        NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
        MDX=MD-MDMIN+1
         MM=NTH0+MD
         IF(ABS(MM).EQ.1) THEN
            CEFLD(1,MDX,NDX,NR)=A1*CEFLD(1,MDX,NDX,NR+1)
     &                         +A2*CEFLD(1,MDX,NDX,NR+2)
            CEFLD(2,MDX,NDX,NR)=A1*CEFLD(2,MDX,NDX,NR+1)
     &                         +A2*CEFLD(2,MDX,NDX,NR+2)
         ELSE
            CEFLD(1,MDX,NDX,NR)=0.D0
            CEFLD(2,MDX,NDX,NR)=0.D0
         ENDIF
         IF(MM.EQ.0) THEN
            CEFLD(3,MDX,NDX,NR)=A1*CEFLD(3,MDX,NDX,NR+1)
     &                         +A2*CEFLD(3,MDX,NDX,NR+2)
            CEFLD(3,MDX,NDX,NR)=CEFLD(3,MDX,NDX,NR+1)
        ELSE
            CEFLD(3,MDX,NDX,NR)=0.D0
         ENDIF
      ENDDO
      ENDDO
C
C     ----- Inverse Fourier transform ----
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO NR=1,NRMAX+1
      DO I=1,3
         CEFM1=0d0
         DO NDX=1,NDSIZ
         NDXM=NHC*(NDX-1)+1
         IF (NHC /= 1)THEN
             NPH  = NPH0 + NPHMAX/2 + 1 - NPH0_SV
             NDXM = NDXM + NPH + (NHCF-NHC)*NDSIZ/2
         ENDIF
         DO MDX=1,MDMSIZ
            CEFM1(MDX,NDXM)=CEFLD(I,MDX,NDX,NR)
         ENDDO
         ENDDO
         CALL WMSUBEM(CEFM1,CEFM2)
         DO NHH=1,NDMSIZ
         DO NTH=1,NTHMAX
            CEFLDM(I,NTH,NHH,NR)=CEFM2(NTH,NHH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO NR=1,NRMAX+1
      DO I=1,3
         DO NDX=1,NDSIZ
         DO MDX=1,MDSIZ
            CEF1(MDX,NDX)=CEFLD(I,MDX,NDX,NR)
         ENDDO
         ENDDO
         CALL WMSUBE(CEF1,CEF2)
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CEFLD(I,NTH,NHH,NR)=CEF2(NTH,NHH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- Calculate CEN from CEsup -----
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CEN(1,NTH,NHH,NR)=0.0d0
         CEN(2,NTH,NHH,NR)=0.0d0
         CEN(3,NTH,NHH,NR)=0.0d0

         CEP(1,NTH,NHH,NR)=0.0d0
         CEP(2,NTH,NHH,NR)=0.0d0
         CEP(3,NTH,NHH,NR)=0.0d0
      ENDDO
      ENDDO
      ENDDO

      NR0=1

      DO NR=NR0,NRMAX+1
C
      IF(NR.EQ.1) THEN
         XRI=1.D6/XRHO(2)
         XRL=XRHO(2)/1.D6
      ELSE
         XRI=1.D0/XRHO(NR)
         XRL=XRHO(NR)
      ENDIF
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         NHHF=(NHH-1)*NFACT +1
         NTHF=(NTH-1)*MFACT +1
         NHHD=(NHH-1)*NHCF +1
C
C        ----- Calculate rotation matrix mu=RMA -----
C
         CALL WMCMAG_F(NR,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
         TC2=BSUPTH/BABS
         TC3=BSUPPH/BABS
C
C        ***** RF11=RJ*SQRT(G^11)/XR *****
C
         RF11=SQRT(RG22(NTHF,NHHF,NR)*RG33(NTHF,NHHF,NR)
     &            -RG23(NTHF,NHHF,NR)*RG23(NTHF,NHHF,NR))
         RMA(1,1)= RJ(NTHF,NHHF,NR)/RF11*XRI
         RMA(2,1)= 0.D0
         RMA(3,1)= 0.D0
         RMA(1,2)= (TC2*(RG23(NTHF,NHHF,NR)*RG12(NTHF,NHHF,NR)
     &                  -RG22(NTHF,NHHF,NR)*RG13(NTHF,NHHF,NR))
     &             +TC3*(RG33(NTHF,NHHF,NR)*RG12(NTHF,NHHF,NR)
     &                  -RG23(NTHF,NHHF,NR)*RG13(NTHF,NHHF,NR))*XRI)
     &             /RF11
         RMA(2,2)= TC3*RF11*XRL
         RMA(3,2)=-TC2*RF11*XRL
         RMA(1,3)=TC2*RG12(NTHF,NHHF,NR)
     &           +TC3*RG13(NTHF,NHHF,NR)*XRI
         RMA(2,3)=TC2*RG22(NTHF,NHHF,NR)*XRL*XRL
     &           +TC3*RG23(NTHF,NHHF,NR)*XRL
         RMA(3,3)=TC2*RG23(NTHF,NHHF,NR)*XRL
     &           +TC3*RG33(NTHF,NHHF,NR)
C
         CALL INVMRD(RMA,3,3,ILL)
         IF(ILL.NE.0) THEN
            WRITE(6,*) 'XX WEFLD: INVMRD(RMA) : SINGULAR MATRIX'
            WRITE(6,'(3I5,1P2E12.4)') NR,NTH,NHH,TC3,TC2
         ENDIF
C
C         IF(NR.EQ.2.OR.NR.EQ.3) THEN
C            WRITE(6,*) 'NR,NTH,NHH=',NR,NTH,NHH
C            WRITE(6,'(1P3E12.4)') RMA(1,1),RMA(1,2),RMA(1,3)
C            WRITE(6,'(1P3E12.4)') RMA(2,1),RMA(2,2),RMA(2,3)
C            WRITE(6,'(1P3E12.4)') RMA(3,1),RMA(3,2),RMA(3,3)
C            WRITE(6,'(1P3E12.4)') RJ(NTH,NHH,NR),RF11,XRI
C            WRITE(6,'(1P3E12.4)') TC2,TC3,XRL
C         ENDIF
C
C         IF (XRHO(NR) <=20d0 +1D-20 )THEN
C         IF (NR  <= NR_S )THEN

         IF (NR  < NR_S )THEN
C         IF (NR  <= -1 )THEN
           CEN(1,NTH,NHH,NR) =CEFLD(1,NTH,NHH,NR)
           CEN(2,NTH,NHH,NR) =CEFLD(2,NTH,NHH,NR)
           CEN(3,NTH,NHH,NR) =CEFLD(3,NTH,NHH,NR)
           CEND(1,NTH,NHH,NR)=CEND(1,NTH,NHH,NR)
     &                       +CEFLDM(1,NTH,NHHD,NR)
           CEND(2,NTH,NHH,NR)=CEND(2,NTH,NHH,NR)
     &                       +CEFLDM(2,NTH,NHHD,NR)
           CEND(3,NTH,NHH,NR)=CEND(3,NTH,NHH,NR)
     &                       +CEFLDM(3,NTH,NHHD,NR)
         ELSE
         CEN(1,NTH,NHH,NR) =RMA(1,1)*CEFLD(1,NTH,NHH,NR)
     &                     +RMA(1,2)*CEFLD(2,NTH,NHH,NR)
     &                     +RMA(1,3)*CEFLD(3,NTH,NHH,NR)
         CEN(2,NTH,NHH,NR) =RMA(2,1)*CEFLD(1,NTH,NHH,NR)
     &                     +RMA(2,2)*CEFLD(2,NTH,NHH,NR)
     &                     +RMA(2,3)*CEFLD(3,NTH,NHH,NR)
         CEN(3,NTH,NHH,NR) =RMA(3,1)*CEFLD(1,NTH,NHH,NR)
     &                     +RMA(3,2)*CEFLD(2,NTH,NHH,NR)
     &                     +RMA(3,3)*CEFLD(3,NTH,NHH,NR)

         CEND(1,NTH,NHH,NR)=CEND(1,NTH,NHH,NR)
     &                      +RMA(1,1)*CEFLDM(1,NTH,NHHD,NR)
     &                      +RMA(1,2)*CEFLDM(2,NTH,NHHD,NR)
     &                      +RMA(1,3)*CEFLDM(3,NTH,NHHD,NR)
         CEND(2,NTH,NHH,NR)=CEND(2,NTH,NHH,NR)
     &                      +RMA(2,1)*CEFLDM(1,NTH,NHHD,NR)
     &                      +RMA(2,2)*CEFLDM(2,NTH,NHHD,NR)
     &                      +RMA(2,3)*CEFLDM(3,NTH,NHHD,NR)
         CEND(3,NTH,NHH,NR)=CEND(3,NTH,NHH,NR)
     &                      +RMA(3,1)*CEFLDM(1,NTH,NHHD,NR)
     &                      +RMA(3,2)*CEFLDM(2,NTH,NHHD,NR)
     &                      +RMA(3,3)*CEFLDM(3,NTH,NHHD,NR)
         ENDIF

C         CEN(1,NTH,NHH,NR)=RMA(1,1)*CEFLD(1,NTH,NHH,NR)*XRI
C     &                    +RMA(1,2)*CEFLD(2,NTH,NHH,NR)*XRL
C     &                    +RMA(1,3)*CEFLD(3,NTH,NHH,NR)
C         CEN(2,NTH,NHH,NR)=RMA(2,1)*CEFLD(1,NTH,NHH,NR)*XRI
C     &                    +RMA(2,2)*CEFLD(2,NTH,NHH,NR)*XRL
C     &                    +RMA(2,3)*CEFLD(3,NTH,NHH,NR)
C         CEN(3,NTH,NHH,NR)=RMA(3,1)*CEFLD(1,NTH,NHH,NR)*XRI
C     &                    +RMA(3,2)*CEFLD(2,NTH,NHH,NR)*XRL
C     &                    +RMA(3,3)*CEFLD(3,NTH,NHH,NR)

         CEP(1,NTH,NHH,NR) =(   CEN(1,NTH,NHH,NR)
     &                     + CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
         CEP(2,NTH,NHH,NR) =(   CEN(1,NTH,NHH,NR)
     &                     - CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
         CEP(3,NTH,NHH,NR) =    CEN(3,NTH,NHH,NR)

C         CEPD(1,NTH,NHH,NR)=
C     &                        (  CEND(1,NTH,NHH,NR)
C     &                      + CI*CEND(2,NTH,NHH,NR))/SQRT(2.D0)
C         CEPD(2,NTH,NHH,NR)=
C     &                        (  CEND(1,NTH,NHH,NR)
C     &                      - CI*CEND(2,NTH,NHH,NR))/SQRT(2.D0)
C         CEPD(3,NTH,NHH,NR)=
C     &                        CEND(3,NTH,NHH,NR)
C
C         IF(NR.EQ.1) THEN
C            WRITE(6,'(1P4E12.4)') CEN(1,NTH,NHH,NR),CEP(1,NTH,NHH,NR)
C            WRITE(6,'(1P4E12.4)') CEN(2,NTH,NHH,NR),CEP(2,NTH,NHH,NR)
C            WRITE(6,'(1P4E12.4)') CEN(3,NTH,NHH,NR),CEP(3,NTH,NHH,NR)
C         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
C     ------- Set E+, E- and Ez0 to those outside the plasma ---------
      DO K=1,3
         DO NR=1,NRMAX+1
            IF(XRHO(NR).GT.1.0D0) THEN
               DO NTH=1,NTHMAX
                  DO NHH=1,NHHMAX
                     CEPD(K,NTH,NHH,NR)=(0.D0,0.D0)
                  ENDDO
               ENDDO
            ELSE
            ENDIF
         ENDDO
      ENDDO
C

      DO NHH=1,NHHMAX
         CEN1 =0.D0
         CEN2 =0.D0
         CEN3 =0.D0
         CEP1 =0.D0
         CEP2 =0.D0
         CEP3 =0.D0
         DO NTH=1,NTHMAX
            CEN1 =CEN1 +CEN(1,NTH,NHH,2)
            CEN2 =CEN2 +CEN(2,NTH,NHH,2)
            CEN3 =CEN3 +CEN(3,NTH,NHH,2)
            CEP1 =CEP1 +CEP(1,NTH,NHH,2)
            CEP2 =CEP2 +CEP(2,NTH,NHH,2)
            CEP3 =CEP3 +CEP(3,NTH,NHH,2)
         ENDDO
         CEN1 =CEN1 /NTHMAX
         CEN2 =CEN2 /NTHMAX
         CEN3 =CEN3 /NTHMAX
         CEP1 =CEP1 /NTHMAX
         CEP2 =CEP2 /NTHMAX
         CEP3 =CEP3 /NTHMAX
         DO NTH=1,NTHMAX
            CEN(1,NTH,NHH,1) =CEN1
            CEN(2,NTH,NHH,1) =CEN2
            CEN(3,NTH,NHH,1) =CEN3
            CEP(1,NTH,NHH,1) =CEP1
            CEP(2,NTH,NHH,1) =CEP2
            CEP(3,NTH,NHH,1) =CEP3
         ENDDO
      ENDDO
C
C     ----- Normalize CEFLD -----
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         NHHF=(NHH-1)*NFACT +1
         NTHF=(NTH-1)*MFACT +1
         NHHD=(NHH-1)*NHCF +1
C 
         RF11=(RG22(NTHF,NHHF,NR)*RG33(NTHF,NHHF,NR)
     &        -RG23(NTHF,NHHF,NR)**2)/RJ(NTHF,NHHF,NR)**2
         RF22=(RG33(NTHF,NHHF,NR)*RG11(NTHF,NHHF,NR)
     &        -RG13(NTHF,NHHF,NR)**2)/RJ(NTHF,NHHF,NR)**2
         RF33=(RG11(NTHF,NHHF,NR)*RG22(NTHF,NHHF,NR)
     &         -RG12(NTHF,NHHF,NR)**2)/RJ(NTHF,NHHF,NR)**2
         RG011=SQRT(RF11)
         RG022=SQRT(RF22)
         RG033=SQRT(RF33)
C
C         IF (XRHO(NR) <=20d0 +1D-20 )THEN
         IF (NR  < NR_S )THEN
C         IF (NR  <= -1 )THEN
           CEFLD(1,NTH,NHH,NR) =CEFLD(1,NTH,NHH,NR)
           CEFLD(2,NTH,NHH,NR) =CEFLD(2,NTH,NHH,NR)
           CEFLD(3,NTH,NHH,NR) =CEFLD(3,NTH,NHH,NR)
           CEFLDD(1,NTH,NHH,NR)=CEFLDD(1,NTH,NHH,NR)
     &                         +CEFLDM(1,NTH,NHHD,NR)
           CEFLDD(2,NTH,NHH,NR)=CEFLDD(2,NTH,NHH,NR)
     &                         +CEFLDM(2,NTH,NHHD,NR)
           CEFLDD(3,NTH,NHH,NR)=CEFLDD(3,NTH,NHH,NR)
     +                         +CEFLDM(3,NTH,NHHD,NR)
         ELSE
           CEFLD(1,NTH,NHH,NR)=CEFLD(1,NTH,NHH,NR)*RG011
           CEFLD(2,NTH,NHH,NR)=CEFLD(2,NTH,NHH,NR)*RG022
           CEFLD(3,NTH,NHH,NR)=CEFLD(3,NTH,NHH,NR)*RG033
           CEFLDD(1,NTH,NHH,NR)=CEFLDD(1,NTH,NHH,NR)
     &                         +CEFLDM(1,NTH,NHHD,NR)*RG011
           CEFLDD(2,NTH,NHH,NR)=CEFLDD(1,NTH,NHH,NR)
     &                         +CEFLDM(2,NTH,NHHD,NR)*RG022
           CEFLDD(3,NTH,NHH,NR)=CEFLDD(3,NTH,NHH,NR)
     &                         +CEFLDM(3,NTH,NHHD,NR)*RG033
         ENDIF
C
C
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE ELECTRIC FIELD (ONLY loop)******
C
      SUBROUTINE WMEFLD_POST
C
      INCLUDE 'wmcomm.inc'

      NR0=1
      DO NR=NR0,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CEPD(1,NTH,NHH,NR)=
     &                        (  CEND(1,NTH,NHH,NR)
     &                      + CI*CEND(2,NTH,NHH,NR))/SQRT(2.D0)
         CEPD(2,NTH,NHH,NR)=
     &                        (  CEND(1,NTH,NHH,NR)
     &                      - CI*CEND(2,NTH,NHH,NR))/SQRT(2.D0)
         CEPD(3,NTH,NHH,NR)=
     &                        CEND(3,NTH,NHH,NR)
      ENDDO
      ENDDO
      ENDDO
      DO K=1,3
         DO NR=1,NRMAX+1
            IF(XRHO(NR).GT.1.0D0) THEN
               DO NTH=1,NTHMAX
                  DO NHH=1,NHHMAX
                     CEPD(K,NTH,NHH,NR)=(0.D0,0.D0)
                  ENDDO
               ENDDO
            ELSE
            ENDIF
         ENDDO
      ENDDO
      DO NHH=1,NHHMAX
         CEND1 =0.D0
         CEND2 =0.D0
         CEND3 =0.D0
         CEPD1 =0.D0
         CEPD2 =0.D0
         CEPD3 =0.D0
         DO NTH=1,NTHMAX
            CEND1 =CEND1 +CEND(1,NTH,NHH,2)
            CEND2 =CEND2 +CEND(2,NTH,NHH,2)
            CEND3 =CEND3 +CEND(3,NTH,NHH,2)
            CEPD1 =CEPD1 +CEPD(1,NTH,NHH,2)
            CEPD2 =CEPD2 +CEPD(2,NTH,NHH,2)
            CEPD3 =CEPD3 +CEPD(3,NTH,NHH,2)
         ENDDO
         CEND1 =CEND1 /NTHMAX
         CEND2 =CEND2 /NTHMAX
         CEND3 =CEND3 /NTHMAX
         CEPD1 =CEPD1 /NTHMAX
         CEPD2 =CEPD2 /NTHMAX
         CEPD3 =CEPD3 /NTHMAX
         DO NTH=1,NTHMAX
            CEND(1,NTH,NHH,1) =CEND1
            CEND(2,NTH,NHH,1) =CEND2
            CEND(3,NTH,NHH,1) =CEND3
            CEPD(1,NTH,NHH,1) =CEPD1
            CEPD(2,NTH,NHH,1) =CEPD2
            CEPD(3,NTH,NHH,1) =CEPD3
         ENDDO
      ENDDO
      END
C
C     ****** CALCULATE MAGNETIC FIELD ******
C
      SUBROUTINE WMBFLD
C
      INCLUDE 'wmcomm.inc'
      DIMENSION CBF1(MDM,NDM),CBF2(MDM,NDM)
      DIMENSION CBFM1(MDMM,NDMM),CBFM2(MDMM,NDMM)
      DIMENSION CBFLDM(3,MDM,NDMM,NRM)
C
      CW=2.D0*PI*CRF*1.D6
C
      DRHO1=(XRHO(2)-XRHO(1))**2
      DRHO2=(XRHO(3)-XRHO(1))**2
      A1= DRHO2/(DRHO2-DRHO1)
      A2=-DRHO1/(DRHO2-DRHO1)
C
      DO NR=2,NRMAX+1
C
         XRHOM=       XRHO(NR-1)
         XRHOMH=0.5D0*(XRHO(NR-1)+XRHO(NR))
         XRHOC=       XRHO(NR)
         IF(MODELG.EQ.3) THEN
            QPMH=0.5D0*(QPS(NR-1)+QPS(NR))
            DPSIPDRHOMH=2.D0*PSITA*XRHOMH/QPMH
         ELSE
            DPSIPDRHOMH=2.D0*PSIPA*XRHOMH
         ENDIF
C
         DRHOM=XRHO(NR)-XRHO(NR-1)
C
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            MM=NTH0+MD
C
            CBFLD1=-(CI/CW)
     &            *(CI*MM*CEFLDK(3,MDX,NDX,NR)
     &             -CI*NN*CEFLDK(2,MDX,NDX,NR))/XRHOC
            CBFLD2=-(CI/CW)
     &            *(CI*NN*CEFLDK(1,MDX,NDX,NR)
     &             -(CEFLDK(3,MDX,NDX,NR  )
     &              -CEFLDK(3,MDX,NDX,NR-1))/(DPSIPDRHOMH*DRHOM)
     &             )*XRHOM
            CBFLD3=-(CI/CW)
     &            *((CEFLDK(2,MDX,NDX,NR  )
     &              -CEFLDK(2,MDX,NDX,NR-1))/(DPSIPDRHOMH*DRHOM)
     &            -CI*MM*CEFLDK(1,MDX,NDX,NR))
C
            CBFLDK(1,MDX,NDX,NR)=VC*CBFLD1
            CBFLDK(2,MDX,NDX,NR)=VC*CBFLD2
            CBFLDK(3,MDX,NDX,NR)=VC*CBFLD3
C
         ENDDO
         ENDDO
C
      ENDDO
C
      NR=1
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
         MM=NTH0+MD
C
         IF(ABS(MM).EQ.1) THEN
            CBFLDK(1,MDX,NDX,NR)=A1*CBFLDK(1,MDX,NDX,NR+1)
     &                          +A2*CBFLDK(1,MDX,NDX,NR+2)
         ELSE
            CBFLDK(1,MDX,NDX,NR)=0.D0
         ENDIF
         CBFLDK(2,MDX,NDX,NR)=CBFLDK(2,MDX,NDX,2)
         CBFLDK(3,MDX,NDX,NR)=CBFLDK(3,MDX,NDX,2)
C
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX+1
      DO I=1,3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO NDX=1,NDSIZ
         NDXM=(NDX-1)*NHC + 1
         IF (NHC /= 1)THEN
             NPH  = NPH0 + NPHMAX/2 + 1 - NPH0_SV
             NDXM = NDXM + NPH + (NHCF-NHC)*NDSIZ/2
         ENDIF
         DO MDX=1,MDMSIZ
            CBFM1(MDX,NDXM)=CBFLDK(I,MDX,NDX,NR)
         ENDDO
         ENDDO
         CALL WMSUBEM(CBFM1,CBFM2)
         DO NHH=1,NDMSIZ
         DO NTH=1,NTHMAX
            CBFLDM(I,NTH,NHH,NR)=CBFM2(NTH,NHH)
         ENDDO
         ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO NDX=1,NDSIZ
         DO MDX=1,MDSIZ
            CBF1(MDX,NDX)=CBFLDK(I,MDX,NDX,NR)
         ENDDO
         ENDDO
         CALL WMSUBE(CBF1,CBF2)
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CBFLD(I,NTH,NHH,NR)=CBF2(NTH,NHH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         NHHF=(NHH-1)*NFACT +1
         NTHF=(NTH-1)*MFACT +1
         NHHD=(NHH-1)*NHCF +1
C
         CBFLD(1,NTH,NHH,NR)=CBFLD(1,NTH,NHH,NR)/RJ(NTHF,NHHF,NR)
     &                      *SQRT(RG11(NTHF,NHHF,NR))
         CBFLD(2,NTH,NHH,NR)=CBFLD(2,NTH,NHH,NR)/RJ(NTHF,NHHF,NR)
     &                      *SQRT(RG22(NTHF,NHHF,NR))
         CBFLD(3,NTH,NHH,NR)=CBFLD(3,NTH,NHH,NR)/RJ(NTHF,NHHF,NR)
     &                      *SQRT(RG33(NTHF,NHHF,NR))
         CBFLDD(1,NTH,NHH,NR)=CBFLDD(1,NTH,NHH,NR)
     &                       +CBFLDM(1,NTH,NHHD,NR)/RJ(NTHF,NHHF,NR)
     &                      *SQRT(RG11(NTHF,NHHF,NR))
         CBFLDD(2,NTH,NHH,NR)=CBFLDD(2,NTH,NHH,NR)
     &                       +CBFLDM(2,NTH,NHHD,NR)/RJ(NTHF,NHHF,NR)
     &                      *SQRT(RG22(NTHF,NHHF,NR))
         CBFLDD(3,NTH,NHH,NR)=CBFLDD(3,NTH,NHH,NR)
     &                       +CBFLDM(3,NTH,NHHD,NR)/RJ(NTHF,NHHF,NR)
     &                      *SQRT(RG33(NTHF,NHHF,NR))
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** HERMITIAN ******
C
      COMPLEX*16 FUNCTION CHERMIT(CX,CY)
C
      COMPLEX*16 CX,CY
C
C      CHERMIT=CX
      CHERMIT=0.5D0*(CX-DCONJG(CY))
      RETURN
      END
C
C     ****** CALCULATE ENERGY FLUX ******
C
      SUBROUTINE WMPFLX
C
      INCLUDE 'wmcomm.inc'
C
      RETURN
      END
C
C     ****** CALCULATE ANTENNA IMPEDANCE ******
C
      SUBROUTINE WMPANT
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CFVP(NDM,MDM,3)
C
      IF(MODEEG.EQ.0.AND.MODEWG.EQ.0) THEN
         NRANT=0
         DO NRI=1,NRMAX 
            IF(XR(NRI)/RD.LT.1.D0) NRANT=NRI
         ENDDO
C
         DTH=2.D0*PI/DBLE(NTHMAX)
         DPH=2.D0*PI/DBLE(NHHMAX)/DBLE(NHC)
         CPRAD=(0.D0,0.D0)
C
         DO NR=NRANT-1,NRMAX
C
            CALL WMSETM_V(NR,CFVP)
C
            DO MDX=1,MDSIZ
               DO NDX=1,NDSIZ
                  CCE1=CEFLDK(1,MDX,NDX,NR+1)
                  CCE2=CEFLDK(2,MDX,NDX,NR+1)
                  CCE3=CEFLDK(3,MDX,NDX,NR+1)
                  CPRAD=CPRAD+DCONJG(CCE1)*CFVP(NDX,MDX,1)
     &                       +DCONJG(CCE2)*CFVP(NDX,MDX,2)
     &                       +DCONJG(CCE3)*CFVP(NDX,MDX,3)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         CPRAD=0.D0
      ENDIF
      CRADTT=CPRAD
      CRADDTT=CRADDTT
     &       + CRADTT
      IF (NPHMAX > 1)THEN
        NPH = NPH0+NPHMAX/2+1 - NPH0_SV
        CRADFTT(NPH)=CRADDTT
      ENDIF

C      WRITE(6,'(A,1P2E12.4)') '# CPRAD=',CPRAD

      RETURN
      END
C
C     ****** CALCULATE CURRENT DRIVE EFFICIENCY ******
C
C      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
C      Z  = ZEFF   : EFFECTIVE Z
C      XL = X / RR : NORMALIZED X
C      YL = Y / RR : NORMALIZED Y
C      ID : 0 : LANDAU DAMPING
C           1 : TTMP
C
      FUNCTION W1CDEF(WT,Z,XL,YL,ID)
C
      INCLUDE 'wmcomm.inc'
C
      R=SQRT(XL*XL+YL*YL)
      IF(ID.EQ.0) THEN
         D=3.D0/Z
         QC=3.83D0
         A=0.D0
         RM=1.38D0
         RC=0.389D0
      ELSE
         D=11.91D0/(0.678D0+Z)
         QC=4.13D0
         A=12.3D0
         RM=2.48D0
         RC=0.0987D0
      ENDIF
      IF(WT.LE.1.D-20) THEN
         W=1.D-20
      ELSE
         W=WT
      ENDIF
      EFF0=D/W+QC/Z**0.707D0+4.D0*W*W/(5.D0+Z)
      EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)
C
      Y2=(R+XL)/(1.D0+R)
      IF(Y2.LT.0.D0) Y2=0.D0
      Y1=SQRT(Y2)
      EFF2=1.D0+A*(Y1/W)**3
C
      IF(Y2.LE.1.D-20) THEN
         YT=(1.D0-Y2)*WT*WT/1.D-60
      ELSE
         YT=(1.D0-Y2)*WT*WT/Y2
      ENDIF
      IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
         ARG=(RC*YT)**RM
         IF(ARG.LE.100.D0) THEN
            EFF3=1.D0-MIN(EXP(-ARG),1.D0)
         ELSE
            EFF3=1.D0
         ENDIF
      ELSE
         EFF3=1.D0
      ENDIF
C
      W1CDEF=EFF0*EFF1*EFF2*EFF3
C
      RETURN
      END
C
C     ****** DISPLAY OUTPUT DATA ******
C
      SUBROUTINE WMPOUT
C
      INCLUDE 'wmcomm.inc'
C
      IF(NPRINT.LT.1) RETURN
C
      IF(NRANK.EQ.0) THEN
         WRITE(6,601) CRADTT,PABSTT,PCURT
         WRITE(6,602) (PABST(NS),NS=1,NSMAX)
      ENDIF
C
      IF(NPRINT.LT.2) RETURN
C
      IF(NRANK.EQ.0) THEN
C         WRITE(6,*) '   MD   IMP                        PABSKT'
C         DO ND=NDMIN,NDMAX
C            NDX=ND-NDMIN+1
C            NN=NPH0+NHC*ND
C         DO MD=MDMIN,MDMAX
C            MDX=MD-MDMIN+1
C            MM=NTH0+MD
C            WRITE(6,603) NN,MM,CRADKT(MDX,NDX),
C     &                   (PABSKT(MDX,NDX,NS),NS=1,NSMAX)
C         ENDDO
C         ENDDO
      ENDIF
C
      IF(NPRINT.LT.4) RETURN
C
      DO NR=1,NRMAX+1
         IF(NRANK.EQ.0) 
     &        WRITE(6,'(A,I3,1P6E10.2)') 'NR,E1,E2,E3=',NR,
     &               CEFLD(1,1,1,NR),CEFLD(2,1,1,NR),CEFLD(3,1,1,NR)
      ENDDO
C
      RETURN
C
  601 FORMAT(' ',5X,3X,'RANT=',1PE12.4,10X,'LANT=',1PE12.4/
     &       ' ',5X,3X,'PABS=',1PE12.4,10X,'IDRV=',1PE12.4)
  602 FORMAT(' ',5X,3X,27X,'PABS=',1P3E12.4)
C  603 FORMAT(' ',I5,I5,3X,1P2E12.4,3X,1P3E12.4)
C  604 FORMAT(' ',8X,3F8.4,24X,2F8.4)
C  605 FORMAT(' ',9F8.4)
C
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBH(RF1,CF2)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RF1(MDM,NDM),CF2(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=RF1(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF2(LDX,NHH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         DO KDX=1,KDSIZ
            CF2(LDX,KDX)=CFN(KDX)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBC(CF1)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CF1(MDM,NDM),CF2(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF2(LDX,NHH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         DO KDX=1,KDSIZ
            CF1(LDX,KDX)=CFN(KDX)
         ENDDO
      ENDDO
      RETURN
      END
C     
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBE(CF1,CF2)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CF1(MDM,NDM),CF2(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)
C
      DO NDX=1,NDSIZ
         DO MDX=1,MDSIZ
            CFM(MDX)=CF1(MDX,NDX)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,1)
         DO NTH=1,NTHMAX
            CF2(NTH,NDX)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         DO NDX=1,NDSIZ
            CFN(NDX)=CF2(NTH,NDX)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,1)
         DO NHH=1,NHHMAX
            CF2(NTH,NHH)=CFN(NHH)
         ENDDO
      ENDDO
      RETURN
      END
C     
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBE_F(CF1,CF2)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CF1(MDMF,NDMF),CF2(MDMF,NDMF)
      DIMENSION CFM(MDMF),CFN(NDMF)
C
      DO NDX=1,NDSIZ_F
         DO MDX=1,MDSIZ_F
            CFM(MDX)=CF1(MDX,NDX)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX_F,1)
         DO NTH=1,NTHMAX_F
            CF2(NTH,NDX)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX_F
         DO NDX=1,NDSIZ_F
            CFN(NDX)=CF2(NTH,NDX)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX_F,1)
         DO NHH=1,NHHMAX_F
            CF2(NTH,NHH)=CFN(NHH)
         ENDDO
      ENDDO
      RETURN
      END
C     
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBEM(CF1,CF2)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CF1(MDMM,NDMM),CF2(MDMM,NDMM)
      DIMENSION CFM(MDMM),CFN(NDMM)
C
      DO NDX=1,NDMSIZ
         DO MDX=1,MDMSIZ
            CFM(MDX)=CF1(MDX,NDX)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,1)
         DO NTH=1,NTHMAX
            CF2(NTH,NDX)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         DO NDX=1,NDMSIZ
            CFN(NDX)=CF2(NTH,NDX)
         ENDDO
         CALL WMXFFT(CFN,NDMSIZ,1)
         DO NHH=1,NDMSIZ
            CF2(NTH,NHH)=CFN(NHH)
         ENDDO
      ENDDO
      RETURN
      END

C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBD(RF1,CF2)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RF1(MDM,NDM),CF2(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=1.D0/RF1(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF2(LDX,NHH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         DO KDX=1,KDSIZ
            CF2(LDX,KDX)=CFN(KDX)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** CALCULATE ABSORBED POWER ******
C
      SUBROUTINE WMPABS
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION RNC(NSM),RTPRC(NSM),RTPPC(NSM),RUC(NSM)
      DIMENSION DS(NRM),DSS(NTHM,NHHM,NRM)
      DIMENSION CPF1(MDM,NDM),CPF2(MDM,NDM)
      DIMENSION CPMF1(MDMM,NDMM),CPMF2(MDMM,NDMM)
      DIMENSION CPF1_B(MDM,NDM)
      DIMENSION CDV(3,3,3),CDW(3,3,3)
      DIMENSION iposa(nsize),ilena(nsize)
      COMPLEX*16, ALLOCATABLE,
     &           DIMENSION(:,:,:,:) ::CPABSKM
      COMPLEX*16, ALLOCATABLE,
     &           DIMENSION(:,:,:,:) ::CPABSKC
      COMPLEX*16, ALLOCATABLE,
     &           DIMENSION(:) ::CFA,CFB
C
      CW=2.D0*PI*CRF*1.D6
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NHHMAX)/DBLE(NHC)
      DPHM=2.D0*PI/DBLE(NHHMAX_SV)
C
      ALLOCATE(CPABSKM(MDM,MDM,NDM,NDM))
      ALLOCATE(CPABSKC(MDM,MDM,NDM,NDM))
      NM=NRM*NSM*MDM*NDM

      IF(NRANK.EQ.0) THEN
         NBSTX=NBST
      ELSE
         NBSTX=NBST-1
      ENDIF
      IF(NRANK.EQ.NSIZE-1) THEN
         NBEDX=NBED+1
      ELSE
         NBEDX=NBED
      ENDIF
C         
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         PCUR(NTH,NHH,NR)=0.D0
      ENDDO
      ENDDO
      ENDDO

      NR=NBSTX
      DO NS=1,NSMAX
C
         IF(NR.LE.NR_S-1) THEN
         CALL WMSETF(NR,NS)
         DO NDX=1,NDSIZ
         DO KDX=1,KDSIZ_F
         DO MDX=1,MDSIZ
         DO LDX=1,LDSIZ_F
         DO J=1,3
         DO I=1,3
            CPSF(I,J,LDX,MDX,KDX,NDX,2)=CPSF(I,J,LDX,MDX,KDX,NDX,3)
            CPSHF(I,J,LDX,MDX,KDX,NDX,2)=CPSHF(I,J,LDX,MDX,KDX,NDX,3)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ELSEIF(NR.GE.NR_S+2) THEN
         CALL WMSETF_OUT(NR,NS)
         DO NDX=1,NDSIZ
         DO KDX=1,KDSIZ_F
         DO MDX=1,MDSIZ
         DO LDX=1,LDSIZ_F
         DO J=1,3
         DO I=1,3
            CGDD(I,J,LDX,MDX,KDX,NDX,2)=CGDD(I,J,LDX,MDX,KDX,NDX,3)
            CGDDH(I,J,LDX,MDX,KDX,NDX,2)=CGDDH(I,J,LDX,MDX,KDX,NDX,3)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDIF
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CPABS(NTH,NHH,NR,NS)=(0.D0,0.D0)
      ENDDO
      ENDDO
      ENDDO

      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CPABSK(NTH,NHH,NR,NS)=(0.D0,0.D0)
      ENDDO
      ENDDO
      ENDDO

      DO NB=NBSTX,NBED
         PRINT *,NB
C
         IF(NB.LE.NR_S-1) THEN
            NR=NB
         ELSEIF(NB.GE.NR_S+2) THEN
            NR=NB-2
            IF (NB .EQ. NR_S +2)THEN
            CALL WMSETF_OUT(NR,NS)
            DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ_F
            DO MDX=1,MDSIZ
            DO LDX=1,LDSIZ_F
            DO J=1,3
            DO I=1,3
               CGDD(I,J,LDX,MDX,KDX,NDX,2)=CGDD(I,J,LDX,MDX,KDX,NDX,3)
               CGDDH(I,J,LDX,MDX,KDX,NDX,2)=CGDDH(I,J,LDX,MDX,KDX,NDX,3)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

            ENDIF
         ELSE
            CYCLE
         END IF
         DO NKX=1,NDSIZ
         DO KDX=1,KDSIZ
         DO MLX=1,MDSIZ
         DO LDX=1,LDSIZ
            CPABSKM(LDX,MLX,KDX,NKX)=(0.D0,0.D0)
            CPABSKC(LDX,MLX,KDX,NKX)=(0.D0,0.D0)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

         IF (NB <= NR_S - 1)THEN
            CALL WMSETF(NR+1,NS)
            DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ_F
            DO MDX=1,MDSIZ
            DO LDX=1,LDSIZ_F
            DO J=1,3
            DO I=1,3
               CPSF(I,J,LDX,MDX,KDX,NDX,1)=CPSF(I,J,LDX,MDX,KDX,NDX,2)
               CPSF(I,J,LDX,MDX,KDX,NDX,2)=CPSF(I,J,LDX,MDX,KDX,NDX,3)
               CPSHF(I,J,LDX,MDX,KDX,NDX,1)=CPSHF(I,J,LDX,MDX,KDX,NDX,2)
               CPSHF(I,J,LDX,MDX,KDX,NDX,2)=CPSHF(I,J,LDX,MDX,KDX,NDX,3)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO


         ELSE
C            CALL WMSETF_OUT(NR,NS)
C            DO NDX=1,NDSIZ
C            DO KDX=1,KDSIZ_F
C            DO MDX=1,MDSIZ
C            DO LDX=1,LDSIZ_F
C            DO J=1,3
C            DO I=1,3
C               CGDD(I,J,LDX,MDX,KDX,NDX,1)=CGDD(I,J,LDX,MDX,KDX,NDX,2)
C               CGDD(I,J,LDX,MDX,KDX,NDX,2)=CGDD(I,J,LDX,MDX,KDX,NDX,3)
C               CGDDH(I,J,LDX,MDX,KDX,NDX,1)=CGDDH(I,J,LDX,MDX,KDX,NDX,2)
C               CGDDH(I,J,LDX,MDX,KDX,NDX,2)=CGDDH(I,J,LDX,MDX,KDX,NDX,3)
C            ENDDO
C            ENDDO
C            ENDDO
C            ENDDO
C            ENDDO
C            ENDDO

            CALL WMSETF_OUT(NR+1,NS)
            DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ_F
            DO MDX=1,MDSIZ
            DO LDX=1,LDSIZ_F
            DO J=1,3
            DO I=1,3
               CGDD(I,J,LDX,MDX,KDX,NDX,1)=CGDD(I,J,LDX,MDX,KDX,NDX,2)
               CGDD(I,J,LDX,MDX,KDX,NDX,2)=CGDD(I,J,LDX,MDX,KDX,NDX,3)
               CGDDH(I,J,LDX,MDX,KDX,NDX,1)=CGDDH(I,J,LDX,MDX,KDX,NDX,2)
               CGDDH(I,J,LDX,MDX,KDX,NDX,2)=CGDDH(I,J,LDX,MDX,KDX,NDX,3)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
C
            IF(NR.EQ.1) THEN
               XRHOM=XRHO(2)/1.D6
            ELSE
               XRHOM =XRHO(NR)
            ENDIF
            XRHOC =       XRHO(NR+1)
            IF(NR.EQ.NRMAX) THEN
               XRHOP=XRHO(NR+1)
            ELSE
               XRHOP=XRHO(NR+2)
            ENDIF
            XRHOMH=0.5D0*(XRHOC+XRHOM)
            XRHOPH=0.5D0*(XRHOC+XRHOP)
C
            DRHOM =XRHO(NR+1)-XRHO(NR)
            IF(NR.EQ.NRMAX) THEN
               NRPP=NR+1
               DRHOP=XRHO(NR+1)-XRHO(NR)
            ELSE
               NRPP=NR+2
               DRHOP=XRHO(NR+2)-XRHO(NR+1)
            ENDIF
C
            IF(MODELG.EQ.3) THEN
               QPMH=0.5D0*(QPS(NR)+QPS(NR+1))
               QPC=QPS(NR+1)
               DPSIPDRHOMH=2.D0*PSITA*XRHOMH/QPMH
               DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
            ELSE
               DPSIPDRHOMH=2.D0*PSIPA*XRHOMH
               DPSIPDRHOC =2.D0*PSIPA*XRHOC
            ENDIF
C
            DRHOPM= 0.5D0*(DRHOM+DRHOP)
C
            FMHM=0.5D0
            FMHC=0.5D0
C
            FCMH  = DRHOM/(2.D0*DRHOPM)
            FCPH  = DRHOP/(2.D0*DRHOPM)
C
C        ND : (n - n0) / Np
C        KD : n'              = (k - n) / Np
C        NN : n               = n0 +  ND       * Np
C        NK : k = n + n' * Np = n0 + (ND + KD) * Np
C
C        MD : m - m0
C        LD : m'              = l - m
C        MM : m               = m0 +  MD
C        ML : l = m + m'      = m0 + (MD + LD)
C
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
            DO NK=NDMIN,NDMAX
               NKX=NK-NDMIN+1
            KDMAX_TMP=KDMAX-1
            IF(KDMAX==0)KDMAX_TMP=0
            DO KK=KDMIN,KDMAX_TMP
               KKX=KK-KDMIN+1
               KD=NK-ND
               IF (KD+KK .GE. KDMIN_F .AND.
     &            (KD+KK .LE. KDMAX_F.OR.KDMAX_F==0))THEN
                 KX1=KD+KK-KDMIN_F+1
                 NX1=MOD( ND   -NDMIN+4*NDSIZ,NDSIZ)+1
C
           DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
            DO ML=MDMIN,MDMAX
               MLX=ML-MDMIN+1
            LDMAX_TMP=LDMAX-1
            IF(LDMAX==0)LDMAX_TMP=0
            DO LL=LDMIN,LDMAX_TMP
               LLX=LL-LDMIN+1
               LD= ML-MD
               IF (LD+LL .GE. LDMIN_F .AND.
     &            (LD+LL .LE. LDMAX_F .OR. LDMAX_F==0))THEN
                  LX1=LD+LL-LDMIN_F+1

                 MX1=MOD( MD   -MDMIN+4*MDSIZ,MDSIZ)+1
               IF (NR <= NR_S - 1)THEN
                 DO K=1,2
                 DO J=1,3
                 DO I=1,3
                   CDV(I,J,K)=CPSF(I,J,LX1,MX1,KX1,NX1,K)
                   CDW(I,J,K)=CPSF(I,J,LX1,MX1,KX1,NX1,K)
                 ENDDO
                 ENDDO
                 ENDDO
C
                FACT1M=XRHOMH/XRHOM
                FACT1C=XRHOMH/XRHOC
                FACT2C=XRHOMH/XRHOC
                FACT2P=XRHOPH/XRHOC
                FACT3C=XRHOMH/XRHOC
                FACT3P=XRHOPH/XRHOC
C
                CDV11M=CDV(1,1,1)
                CDV11C=CDV(1,1,2)
                CDV12M=CDV(1,2,1)
                CDV12C=CDV(1,2,2)
                CDV13M=CDV(1,3,1)
                CDV13C=CDV(1,3,2)
                CDV21C=CDV(2,1,2)
                CDV22C=CDV(2,2,2)
                CDV23C=CDV(2,3,2)
                CDV31C=CDV(3,1,2)
                CDV32C=CDV(3,2,2)
                CDV33C=CDV(3,3,2)
C
C     --- R COMPONENT OF MAXWELL EQUATION ---
C

                CEMM11=0.5D0*CDV11M*FACT1M ! /XRHOMH/XRHOMH 
                CEMC11=0.5D0*CDV11C*FACT1C ! /XRHOMH/XRHOMH
                CEMM12= FMHM*CDV12M        !       /XRHOMH
                CEMC12= FMHC*CDV12C        !       /XRHOMH
                CEMM13= FMHM*CDV13M        !       /XRHOMH
                CEMC13= FMHC*CDV13C        !       /XRHOMH

C     --- THETA COMPONENT OF MAXWELL EQUATION ---
C

                CEMC21= FCMH*CDV21C*FACT2C  !      /XRHOMH
                CEMP21= FCPH*CDV21C*FACT2P  !      /XRHOPH
                CEMC22=      CDV22C       
                CEMC23=      CDV23C        

C     --- PHI COMPONENT OF MAXWELL EQUATION ---
C

                CEMC31= FCMH*CDV31C*FACT3C !       /XRHOMH
                CEMP31= FCPH*CDV31C*FACT3P !       /XRHOPH
                CEMC32=      CDV32C       
                CEMC33=      CDV33C

               ELSE
                 DO K=1,2
                 DO J=1,3
                 DO I=1,3
                   CDV(I,J,K)=CGDD(I,J,LX1,MX1,KX1,NX1,K)
                   CDW(I,J,K)=CGDD(I,J,LX1,MX1,KX1,NX1,K)
                 ENDDO
                 ENDDO
                 ENDDO
C
                FACT1M=XRHOMH/XRHOM
                FACT1C=XRHOMH/XRHOC
                FACT2C=XRHOMH/XRHOC
                FACT2P=XRHOPH/XRHOC
                FACT3C=XRHOMH/XRHOC
                FACT3P=XRHOPH/XRHOC
C
                CDV11M=CDV(1,1,1)
                CDV11C=CDV(1,1,2)
                CDV12M=CDV(1,2,1)
                CDV12C=CDV(1,2,2)
                CDV13M=CDV(1,3,1)
                CDV13C=CDV(1,3,2)
                CDV21C=CDV(2,1,2)
                CDV22C=CDV(2,2,2)
                CDV23C=CDV(2,3,2)
                CDV31C=CDV(3,1,2)
                CDV32C=CDV(3,2,2)
                CDV33C=CDV(3,3,2)
C
C     --- R COMPONENT OF MAXWELL EQUATION ---
C
                CEMM11=0.5D0*CDV11M*FACT1M /XRHOMH/XRHOMH
                CEMC11=0.5D0*CDV11C*FACT1C /XRHOMH/XRHOMH
                CEMM12= FMHM*CDV12M        *XRHOM /XRHOMH
                CEMC12= FMHC*CDV12C        *XRHOC /XRHOMH
                CEMM13= FMHM*CDV13M               /XRHOMH
                CEMC13= FMHC*CDV13C               /XRHOMH
C
C     --- THETA COMPONENT OF MAXWELL EQUATION ---
C
                CEMC21= FCMH*CDV21C*FACT2C /XRHOMH*XRHOC
                CEMP21= FCPH*CDV21C*FACT2P /XRHOPH*XRHOC
                CEMC22=      CDV22C        *XRHOC *XRHOC
                CEMC23=      CDV23C               *XRHOC
C
C     --- PHI COMPONENT OF MAXWELL EQUATION ---
C
                CEMC31= FCMH*CDV31C*FACT3C /XRHOMH
                CEMP31= FCPH*CDV31C*FACT3P /XRHOPH
                CEMC32=      CDV32C        *XRHOC
                CEMC33=      CDV33C
               ENDIF


               CCE1=DCONJG(CEFLDK(1,MLX,NKX,NR+1))
               CCE2=DCONJG(CEFLDK(2,MLX,NKX,NR+1))
               CCE3=DCONJG(CEFLDK(3,MLX,NKX,NR+1))
               CPM11=CCE1*CEMM11*CEFLDK(1,MDX,NDX,NR+1)
               CPC11=CCE1*CEMC11*CEFLDK(1,MDX,NDX,NR+1)
               CPM12=CCE1*CEMM12*CEFLDK(2,MDX,NDX,NR)
               CPC12=CCE1*CEMC12*CEFLDK(2,MDX,NDX,NR+1)
               CPM13=CCE1*CEMM13*CEFLDK(3,MDX,NDX,NR)
               CPC13=CCE1*CEMC13*CEFLDK(3,MDX,NDX,NR+1)
               CPC21=CCE2*CEMC21*CEFLDK(1,MDX,NDX,NR+1)
               CPP21=CCE2*CEMP21*CEFLDK(1,MDX,NDX,NRPP)
               CPC22=CCE2*CEMC22*CEFLDK(2,MDX,NDX,NR+1)
               CPC23=CCE2*CEMC23*CEFLDK(3,MDX,NDX,NR+1)
               CPC31=CCE3*CEMC31*CEFLDK(1,MDX,NDX,NR+1)
               CPP31=CCE3*CEMP31*CEFLDK(1,MDX,NDX,NRPP)
               CPC32=CCE3*CEMC32*CEFLDK(2,MDX,NDX,NR+1)
               CPC33=CCE3*CEMC33*CEFLDK(3,MDX,NDX,NR+1)

C     
               CPABSM=CI*EPS0*(VC*VC/CW)
     &            *(0d0
     &             +CPM11
     &             +CPM12
     &             +CPM13
     &             )*DPSIPDRHOMH*DRHOM

               CPABSC=CI*EPS0*(VC*VC/CW)
     &            *(0d0
     &             +CPC11
     &             +CPC12
     &             +CPC13
     &             )*DPSIPDRHOMH*DRHOM
     &            +CI*EPS0*(VC*VC/CW)
     &            *(0d0
     &             +CPC22*2d0
     &             +CPC33*2d0
     &             +CPC23*2d0
     &             +CPC32*2d0
     &             +CPC21
     &             +CPC31
     &             +CPP21
     &             +CPP31
     &             )*DPSIPDRHOC*DRHOPM
               WRITE(6,'(A,2I6,4ES12.4)')
     &              'CPABSMC:',NS,NR,CPABSM,CPABSC
               CPABSKM(LLX,MDX,KKX,NDX)
     &        =CPABSKM(LLX,MDX,KKX,NDX)
     &        +0.5D0*CPABSM
C
               CPABSKC(LLX,MDX,KKX,NDX)
     &        =CPABSKC(LLX,MDX,KKX,NDX)
     &        +0.5D0*CPABSC
            ENDIF

            ENDDO
            ENDDO
            ENDDO
C
            ENDIF

            ENDDO
            ENDDO
            ENDDO
C
C     +++++ CALCULATE PABS IN MODE NUMBER SPACE +++++
C

       IF (NR < NR_S)THEN
            DO NDX=1,NDSIZ
               KK=0
               KKX=KK-KDMIN+1
            DO MDX=1,MDSIZ
              LL=0
              LLX=LL-LDMIN+1
              CPABSK(MDX,NDX,NR,NS)  = CPABSK(MDX,NDX,NR,NS)
     &                                 + CPABSKM(LLX,MDX,KKX,NDX)
              CPABSK(MDX,NDX,NR+1,NS)= CPABSK(MDX,NDX,NR+1,NS)
     &                                 + CPABSKC(LLX,MDX,KKX,NDX)
            ENDDO
            ENDDO
       ELSEIF (NR == NR_S)THEN
            DO NDX=1,NDSIZ
               KK=0
               KKX=KK-KDMIN+1
            DO MDX=1,MDSIZ
              LL=0
              LLX=LL-LDMIN+1
              CPABSK(MDX,NDX,NR,NS)  = CPABSK(MDX,NDX,NR,NS)
     &                                 + CPABSKM(LLX,MDX,KKX,NDX)
            ENDDO
            ENDDO
       ENDIF
C
C     +++++ CALCULATE PABS IN REAL SPACE +++++
C
            DO NDX=1,NDSIZ
            DO MDX=1,MDSIZ
               CPF1=0d0
               CPF2=0d0
               DO KK=KDMIN,KDMAX
                  IF (KDMAX /=0 .and. kk==KDMAX)CYCLE
                  KKX =KK-KDMIN+1
               DO LL=LDMIN,LDMAX
                  IF (LDMAX /= 0 .and. LL==LDMAX)CYCLE
                  LLX=LL-LDMIN+1
                  CPF1(LLX,KKX)= CPABSKM (LLX,MDX,KKX,NDX)
               ENDDO
               ENDDO
               CALL WMSUBE(CPF1,CPF2)
               DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  NHHF=NHH
                  NTHF=NTH
                  CPABS(NTH,NHH,NR,NS)=CPABS(NTH,NHH,NR,NS)
     &                               +CPF2(NTH,NHH)
                  CPABS3D(NTH,NHH,NR,NS)=CPABS3D(NTH,NHH,NR,NS)
     &                                +CPF2(NTH,NHH)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
       IF (NR < NR_S)THEN
            DO NDX=1,NDSIZ
            DO MDX=1,MDSIZ
               CPF1=0d0
               CPF2=0d0
               DO KK=KDMIN,KDMAX
                  IF (KDMAX /=0 .and. kk==KDMAX)CYCLE
                  KKX =KK-KDMIN+1
                  DO LL=LDMIN,LDMAX
                  IF (LDMAX /= 0 .and. LL==LDMAX)CYCLE
                  LLX=LL-LDMIN+1
                  CPF1(LLX,KKX)=CPABSKC(LLX,MDX,KKX,NDX)
               ENDDO
               ENDDO
               CALL WMSUBE(CPF1,CPF2)
               DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  NHHF=NHH
                  NTHF=NTH
                  CPABS(NTH,NHH,NR+1,NS)=CPABS(NTH,NHH,NR+1,NS)
     &                               +CPF2(NTH,NHH)
                  CPABS3D(NTH,NHH,NR+1,NS)=CPABS3D(NTH,NHH,NR+1,NS)
     &                               +CPF2(NTH,NHH)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
       ENDIF
      
C
C     +++++ CALCULATE DRIVEN CURRENT IN REAL SPACE +++++
C
            IF (NS .EQ. 1) THEN
              CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
              VTE=SQRT(RTPR(1)*AEE*1.D3/(PA(1)*AMP))
              WW=DBLE(CW)
              IF(RN(1).LE.0.D0) THEN
                 RLNLMD=15.D0
              ELSE
                 RT=(RTPR(1)+2*RTPP(1))/3.D0
                 RLNLMD=16.1D0 - 1.15D0*LOG10(RN(1))
     &               + 2.30D0*LOG10(RT)
              ENDIF
              CALL WMCDEN(NR+1,RNC,RTPRC,RTPPC,RUC)
              VTEC=SQRT(RTPRC(1)*AEE*1.D3/(PA(1)*AMP))
              WW=DBLE(CW)
              IF(RNC(1).LE.0.D0) THEN
                 RLNLMDC=15.D0
              ELSE
                 RTC=(RTPRC(1)+2*RTPPC(1))/3.D0
                 RLNLMDC=16.1D0 - 1.15D0*LOG10(RNC(1))
     &               + 2.30D0*LOG10(RTC)
              ENDIF
              DO ND=NDMIN,NDMAX
                 NDX=ND-NDMIN+1
                 NN=NPH0+NHC*ND
              DO MD=MDMIN,MDMAX
                 MDX=MD-MDMIN+1
                 MM=NTH0+MD
                 DO KKX=1,KDSIZ
                   DO LLX=1,LDSIZ
                     CPF1(LLX,KKX)=CPABSKM(LLX,MDX,KKX,NDX)
                   ENDDO
                   ENDDO
                   CALL WMSUBE(CPF1,CPF2)
                   DO NHH=1,NHHMAX
                   DO NTH=1,NTHMAX
                     NHHF=(NHH-1)*NFACT +1
                     NTHF=(NTH-1)*MFACT +1
                     CALL WMCMAG_F(NR,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
                     RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
                    IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
                    IF(ABS(WW/RKPR).LT.VC) THEN
                     W=WW/(RKPR*VTE)
                     XL=(RPST(NTHF,NHHF,NR)-RR  )/RR
                     YL=(ZPST(NTHF,NHHF,NR)-0.D0)/RR
                     EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                     IF(W.LT.0.D0) EFCD=-EFCD
                     IF (RN(1).GT.0.D0) THEN
                        PCUR(NTH,NHH,NR)=PCUR(NTH,NHH,NR)
     &                       +0.384D0*RTPR(1)/(AEE*1.D3)*EFCD
     &                       /((RN(1)/1.D20)*RLNLMD)*DBLE(CPF2(NTH,NHH))
     &                       /(2.D0*PI*RPST(NTHF,NHHF,NR))
                     END IF
                    ENDIF
            ENDDO
            ENDDO
               DO KKX=1,KDSIZ
               DO LLX=1,LDSIZ
                  CPF1(LLX,KKX)=CPABSKC(LLX,MDX,KKX,NDX)
               ENDDO
               ENDDO
               CALL WMSUBE(CPF1,CPF2)
               DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  NHHF=(NHH-1)*NFACT +1
                  NTHF=(NTH-1)*MFACT +1
                  CALL WMCMAG_F(NR+1,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
                  RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
                 IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
                 IF(ABS(WW/RKPR).LT.VC) THEN
                 W=WW/(RKPR*VTEC)
                 XL=(RPST(NTHF,NHHF,NR+1)-RR  )/RR
                 YL=(ZPST(NTHF,NHHF,NR+1)-0.D0)/RR
                 EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                 IF(W.LT.0.D0) EFCD=-EFCD
                 IF (RNC(1).GT.0.D0) THEN
                    PCUR(NTH,NHH,NR+1)=PCUR(NTH,NHH,NR+1)
     &                   +0.384D0*RTPRC(1)/(AEE*1.D3)*EFCD
     &                   /((RNC(1)/1.D20)*RLNLMDC)*DBLE(CPF2(NTH,NHH))
     &                   /(2.D0*PI*RPST(NTHF,NHHF,NR+1))
                 END IF
                 ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDIF
         ENDDO
         ENDDO

      DEALLOCATE(CPABSKM)
      DEALLOCATE(CPABSKC)

      NRS=NBST
      IF(NRANK.EQ.NSIZE-1) THEN
         NRE=NBED+1
      ELSE
         NRE=NBED
      ENDIF

      NM=NRM*NSM*MDM*NDM
      ALLOCATE(CFA(NM))
      ALLOCATE(CFB(NM))

      MN=0
      DO NS=1,NSMAX
      DO NR=NRS,NRE
      DO NDX=1,NDSIZ
      DO MDX=1,MDSIZ
         MN=MN+1
         CFB(MN)=CPABSK(MDX,NDX,NR,NS)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL mtx_allgather1_integer(MN,ilena)
      ntot=0
      DO i=1,nsize
         iposa(i)=ntot
         ntot=ntot+ilena(i)
      END DO
      CFA=CFB
      NM=ntot

      IF(NRANK.EQ.0) THEN
         MN=0
         DO NS=1,NSMAX
         DO NR=1,NRMAX+1
         DO NDX=1,NDSIZ
         DO MDX=1,MDSIZ
            MN=MN+1
            CPABSK(MDX,NDX,NR,NS)=CFA(MN)
             PABSK(MDX,NDX,NR,NS)=DBLE(CFA(MN))
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      DEALLOCATE(CFA)
      DEALLOCATE(CFB)


      NM=NRM*NSM*MDM*NDM
      ALLOCATE(CFA(NM))
      ALLOCATE(CFB(NM))

      MN=0
      DO NS=1,NSMAX
      DO NR=NRS,NRE
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         MN=MN+1
         CFB(MN)=CPABS(NTH,NHH,NR,NS)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL mtx_allgather1_integer(MN,ilena)
      ntot=0
      DO i=1,nsize
         iposa(i)=ntot
         ntot=ntot+ilena(i)
      END DO
      CFA=CFB
      NM=ntot

C
      IF(NRANK.EQ.0) THEN
         MN=0
         DO NS=1,NSMAX
         DO NR=1,NRMAX+1
            write(1212,*)"#",NR,NS
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            MN=MN+1
            CPABS(NTH,NHH,NR,NS)=CFA(MN)
            PABS(NTH,NHH,NR,NS)=DBLE(CFA(MN))
            write(1212,*)NTH,NHH,PABS(NTH,NHH,NR,NS) !,CPABS(NTH,NHH,NR,NS)
         ENDDO
            write(1212,*)
         ENDDO
            write(1212,*)
            write(1212,*)
         ENDDO
         ENDDO
      ENDIF
C
      DEALLOCATE(CFA)
      DEALLOCATE(CFB)

      NM=NRM*NSM*MDM*NDM*NHC
      ALLOCATE(CFA(NM))
      ALLOCATE(CFB(NM))

      MN=0
      DO NS=1,NSMAX
      DO NR=NRS,NRE
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         MN=MN+1
         CFB(MN)=CPABS3D(NTH,NHH,NR,NS)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL mtx_allgather1_integer(MN,ilena)
      ntot=0
      DO i=1,nsize
         iposa(i)=ntot
         ntot=ntot+ilena(i)
      END DO
      CFA=CFB
      NM=ntot

C
      IF(NRANK.EQ.0) THEN
         MN=0
         DO NS=1,NSMAX
         DO NR=1,NRMAX+1
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            MN=MN+1
            CPABS3D(NTH,NHH,NR,NS)=CFA(MN)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF
C
      DEALLOCATE(CFA)
      DEALLOCATE(CFB)

      NM=NRM*MDM*NDM
      ALLOCATE(CFA(NM))
      ALLOCATE(CFB(NM))

      MN=0
      DO NR=NRS,NRE
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         MN=MN+1
         CFB(MN)=PCUR(NTH,NHH,NR)
      ENDDO
      ENDDO
      ENDDO

      CALL mtx_allgather1_integer(MN,ilena)
      ntot=0
      DO i=1,nsize
         iposa(i)=ntot
         ntot=ntot+ilena(i)
      END DO
      CALL mtx_gatherv_complex8(CFB,MN,CFA,ntot,ilena,iposa)
      NM=ntot

C
      IF(NRANK.EQ.0) THEN
         MN=0
         DO NR=1,NRMAX+1
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            MN=MN+1
            PCUR(NTH,NHH,NR)=CFA(MN)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      DEALLOCATE(CFA)
      DEALLOCATE(CFB)

      IF(NRANK.EQ.0) THEN
C
      DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH*DPH*NHC
         ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         PABSR(NR,NS)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DPH*NHC
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      PCURT=0.D0
      DO NR=1,NRMAX
         PCURT=PCURT+PCURR(NR)
      ENDDO
C
      DO NS=1,NSMAX
         PABST(NS)=0.D0
         DO NR=1,NRMAX
            PABST(NS)=PABST(NS)+PABSR(NR,NS)
         ENDDO
      ENDDO
C
      PABSTT=0.D0
      DO NS=1,NSMAX
         PABSTT=PABSTT+PABST(NS)
      ENDDO
C
      IF (NPHMAX > 1)THEN
        NPH = NPH0+NPHMAX/2+1 - NPH0_SV
        PABSFTT(NPH)=PABSTT
      ENDIF
      DO NS=1,NSMAX
         IF (NPHMAX > 1)THEN
           PABSFT(NS,NPH)=PABST(NS)
         ENDIF
         DO NR=1,NRMAX
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PABSD(NTH,NHH,NR,NS) =dble(CPABS3D(NTH,NHH,NR,NS))
               CPABSD(NTH,NHH,NR,NS)=     CPABS3D(NTH,NHH,NR,NS)
            ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO NR=1,NRMAX
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCURD(NTH,NHH,NR)=PCURD(NTH,NHH,NR)+PCUR(NTH,NHH,NR)
         ENDDO
         ENDDO
      ENDDO
      FACT=1.D0
C
      IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
         FACT=PRFIN/PABSTT
          print *,"FACT",FACT,PRFIN,PABSTT
      ENDIF
      FACTSQ=SQRT(FACT)
C
      NR=1
         DS(NR)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            DSS(NTH,NHH,NR)=0.D0
         ENDDO
         ENDDO
      DO NR=2,NRMAX
         DS(NR)=0.D0
         DRHO=0.5D0*(XRHO(NR+1)-XRHO(NR-1))
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
         NHHF=(NHH-1)*NFACT +1
         NTHF=(NTH-1)*MFACT +1
            IF(MODELG.EQ.3) THEN
               DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
            ELSE
               DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
            ENDIF
            DSSS=DPSIPDRHO*DRHO*RJ(NTHF,NHHF,NR)
            DSS(NTH,NHH,NR)=1.D0/DSSS
C            DS(NR)=DS(NR)+DSSS*DTH*DPH
            DS(NR)=DS(NR)+DSSS*DTH*DPH*NHC
         ENDDO
         ENDDO
         DS(NR)=1.D0/DS(NR)
      ENDDO
C
      PABSTT=FACT*PABSTT
      DO NS=1,NSMAX
         PABST(NS)=FACT*PABST(NS)
         DO NR=1,NRMAX
            PABSR(NR,NS)=FACT*PABSR(NR,NS)*DS(NR)
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NHH,NR,NS)=FACT*PABS(NTH,NHH,NR,NS)
     &                            *DSS(NTH,NHH,NR)
               CPABS(NTH,NHH,NR,NS)=FACT*CPABS(NTH,NHH,NR,NS)
     &                            *DSS(NTH,NHH,NR)
            ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO NR=1,NRMAX
          WRITE(313,'(10(es15.6,1x))')XRHO(NR),(PABSR(NR,NS),NS=1,NSMAX)
          DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
          WRITE(314,'(10(es15.6,1x))')DTH*(NTH-1),DPH*(NHH-1),XRHO(NR), 
     &             (PABS(NTH,NHH,NR,NS),NS=1,NSMAX)
          ENDDO
          WRITE(314,'(10(es15.6,1x))')
          ENDDO
          WRITE(314,'(10(es15.6))')
       ENDDO
          WRITE(313,'(10(es15.6,1x))')
          WRITE(313,'(10(es15.6,1x))')

C
      ENDIF
C
      DO NR=1,NRMAX+1
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NDXH=(NDX-1)*NHC+1
            IF (NPHMAX > 1)THEN
              NPH = NPH0+NPHMAX/2+1 - NPH0_SV
              NDXH=NDXH+NPH
            ENDIF
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDKF(1,MDX,NDXH,NR)=CEFLDK(1,MDX,NDX,NR)
            CEFLDKF(2,MDX,NDXH,NR)=CEFLDK(2,MDX,NDX,NR)
            CEFLDKF(3,MDX,NDXH,NR)=CEFLDK(3,MDX,NDX,NR)
            CBFLDKF(1,MDX,NDXH,NR)=CBFLDK(1,MDX,NDX,NR)
            CBFLDKF(2,MDX,NDXH,NR)=CBFLDK(2,MDX,NDX,NR)
            CBFLDKF(3,MDX,NDXH,NR)=CBFLDK(3,MDX,NDX,NR)
            DO NS=1,NSMAX
              CPABSKF(MDX,NDXH,NR,NS)=CPABSK(MDX,NDX,NR,NS)
               PABSKF(MDX,NDXH,NR,NS)= PABSK(MDX,NDX,NR,NS)
            ENDDO
         ENDDO
         ENDDO
      ENDDO

      DO NR=1,NRMAX+1
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDK(1,MDX,NDX,NR)=FACTSQ*CEFLDK(1,MDX,NDX,NR)
            CEFLDK(2,MDX,NDX,NR)=FACTSQ*CEFLDK(2,MDX,NDX,NR)
            CEFLDK(3,MDX,NDX,NR)=FACTSQ*CEFLDK(3,MDX,NDX,NR)
            CBFLDK(1,MDX,NDX,NR)=FACTSQ*CBFLDK(1,MDX,NDX,NR)
            CBFLDK(2,MDX,NDX,NR)=FACTSQ*CBFLDK(2,MDX,NDX,NR)
            CBFLDK(3,MDX,NDX,NR)=FACTSQ*CBFLDK(3,MDX,NDX,NR)
         ENDDO
         ENDDO
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CEFLD(1,NTH,NHH,NR)=FACTSQ*CEFLD(1,NTH,NHH,NR)
            CEFLD(2,NTH,NHH,NR)=FACTSQ*CEFLD(2,NTH,NHH,NR)
            CEFLD(3,NTH,NHH,NR)=FACTSQ*CEFLD(3,NTH,NHH,NR)
            CBFLD(1,NTH,NHH,NR)=FACTSQ*CBFLD(1,NTH,NHH,NR)
            CBFLD(2,NTH,NHH,NR)=FACTSQ*CBFLD(2,NTH,NHH,NR)
            CBFLD(3,NTH,NHH,NR)=FACTSQ*CBFLD(3,NTH,NHH,NR)
            CEN(1,NTH,NHH,NR)  =FACTSQ*CEN(1,NTH,NHH,NR)
            CEN(2,NTH,NHH,NR)  =FACTSQ*CEN(2,NTH,NHH,NR)
            CEN(3,NTH,NHH,NR)  =FACTSQ*CEN(3,NTH,NHH,NR)
            CEP(1,NTH,NHH,NR)  =FACTSQ*CEP(1,NTH,NHH,NR)
            CEP(2,NTH,NHH,NR)  =FACTSQ*CEP(2,NTH,NHH,NR)
            CEP(3,NTH,NHH,NR)  =FACTSQ*CEP(3,NTH,NHH,NR)
         ENDDO
         ENDDO
      ENDDO
      CALL MTX_BARRIER
C
      RETURN
      END
      SUBROUTINE WMPABS_POST
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION DS(NRM),DSS(NTHM,NHHM,NRM)
      CW=2.D0*PI*CRF*1.D6
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NHHMAX)/DBLE(NHC)
      DPHM=2.D0*PI/DBLE(NHHMAX_SV)
C
      IF(NRANK.EQ.0) THEN
      DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NHH,NR)=PCURD(NTH,NHH,NR)
            PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH*DPH*NHC
         ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         PABSR(NR,NS)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PABS(NTH,NHH,NR,NS) =PABSD(NTH,NHH,NR,NS)
            CPABS(NTH,NHH,NR,NS)=CPABSD(NTH,NHH,NR,NS)
            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DPH*NHC
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      PCURT=0.D0
      DO NR=1,NRMAX
         PCURT=PCURT+PCURR(NR)
      ENDDO
C
      DO NS=1,NSMAX
         PABST(NS)=0.D0
         DO NR=1,NRMAX
            PABST(NS)=PABST(NS)+PABSR(NR,NS)
         ENDDO
      ENDDO

      PABSTT=0.D0
      DO NS=1,NSMAX
         PABSTT=PABSTT+PABST(NS)
      ENDDO
      FACT=1.D0
C
      IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
         FACT=PRFIN/PABSTT
      ENDIF
      FACTSQ=SQRT(FACT)

      NR=1
         DS(NR)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            DSS(NTH,NHH,NR)=0.D0
         ENDDO
         ENDDO
      DO NR=2,NRMAX
         DS(NR)=0.D0
         DRHO=0.5D0*(XRHO(NR+1)-XRHO(NR-1))
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
         NHHF=(NHH-1)*NFACT +1
         NTHF=(NTH-1)*MFACT +1
            IF(MODELG.EQ.3) THEN
               DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
            ELSE
               DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
            ENDIF
            DSSS=DPSIPDRHO*DRHO*RJ(NTHF,NHHF,NR)
            DSS(NTH,NHH,NR)=1.D0/DSSS
C            DS(NR)=DS(NR)+DSSS*DTH*DPH
            DS(NR)=DS(NR)+DSSS*DTH*DPH*NHC
         ENDDO
         ENDDO
         DS(NR)=1.D0/DS(NR)
      ENDDO
C
      PABSTT=FACT*PABSTT
      DO NS=1,NSMAX
         PABST(NS)=FACT*PABST(NS)
         DO NR=1,NRMAX
            PABSR(NR,NS)=FACT*PABSR(NR,NS)*DS(NR)
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NHH,NR,NS)=FACT*PABS(NTH,NHH,NR,NS)
     &                            *DSS(NTH,NHH,NR)
               CPABS(NTH,NHH,NR,NS)=FACT*CPABS(NTH,NHH,NR,NS)
     &                            *DSS(NTH,NHH,NR)
            ENDDO
            ENDDO
         ENDDO
      ENDDO

      

      DO NR=1,NRMAX
          WRITE(313,'(10(es15.6,1x))')XRHO(NR),(PABSR(NR,NS),NS=1,NSMAX)
          DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
          WRITE(314,'(10(es15.6,1x))')DTH*(NTH-1),DPH*(NHH-1),XRHO(NR), 
     &             (PABS(NTH,NHH,NR,NS),NS=1,NSMAX)
          ENDDO
          WRITE(314,'(10(es15.6,1x))')
          ENDDO
          WRITE(314,'(10(es15.6))')
       ENDDO
          WRITE(313,'(10(es15.6,1x))')
          WRITE(313,'(10(es15.6,1x))')

C
      ENDIF
C
      DO NR=1,NRMAX+1
C         DO ND=1,NPHMAX_SV
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDKF(1,MDX,NDX,NR)=FACTSQ*CEFLDKF(1,MDX,NDX,NR)
            CEFLDKF(2,MDX,NDX,NR)=FACTSQ*CEFLDKF(2,MDX,NDX,NR)
            CEFLDKF(3,MDX,NDX,NR)=FACTSQ*CEFLDKF(3,MDX,NDX,NR)
            CBFLDKF(1,MDX,NDX,NR)=FACTSQ*CBFLDKF(1,MDX,NDX,NR)
            CBFLDKF(2,MDX,NDX,NR)=FACTSQ*CBFLDKF(2,MDX,NDX,NR)
            CBFLDKF(3,MDX,NDX,NR)=FACTSQ*CBFLDKF(3,MDX,NDX,NR)
         ENDDO
         ENDDO
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CEFLD(1,NTH,NHH,NR)=FACTSQ*CEFLDD(1,NTH,NHH,NR)
            CEFLD(2,NTH,NHH,NR)=FACTSQ*CEFLDD(2,NTH,NHH,NR)
            CEFLD(3,NTH,NHH,NR)=FACTSQ*CEFLDD(3,NTH,NHH,NR)
            CBFLD(1,NTH,NHH,NR)=FACTSQ*CBFLDD(1,NTH,NHH,NR)
            CBFLD(2,NTH,NHH,NR)=FACTSQ*CBFLDD(2,NTH,NHH,NR)
            CBFLD(3,NTH,NHH,NR)=FACTSQ*CBFLDD(3,NTH,NHH,NR)
            CEN(1,NTH,NHH,NR)  =FACTSQ*CEND(1,NTH,NHH,NR)
            CEN(2,NTH,NHH,NR)  =FACTSQ*CEND(2,NTH,NHH,NR)
            CEN(3,NTH,NHH,NR)  =FACTSQ*CEND(3,NTH,NHH,NR)
            CEP(1,NTH,NHH,NR)  =FACTSQ*CEPD(1,NTH,NHH,NR)
            CEP(2,NTH,NHH,NR)  =FACTSQ*CEPD(2,NTH,NHH,NR)
            CEP(3,NTH,NHH,NR)  =FACTSQ*CEPD(3,NTH,NHH,NR)
         ENDDO
         ENDDO
      ENDDO
      CALL mtx_barrier
C
      RETURN
      END


      SUBROUTINE WMPOUT_INIT
C
      INCLUDE 'wmcomm.inc'
C
      CRADDTT=0d0
C
      DO NS=1,NSMAX
         DO NR=1,NRMAX
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PABSD(NTH,NHH,NR,NS)=0d0
               CPABSD(NTH,NHH,NR,NS)=0d0
            ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO NR=1,NRMAX
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCURD(NTH,NHH,NR)=0d0
         ENDDO
         ENDDO
      ENDDO
      DO NR=1,NRMAX+1
C         DO ND=1,NPHMAX_SV
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDKF(1,MDX,NDX,NR)=0d0
            CEFLDKF(2,MDX,NDX,NR)=0d0
            CEFLDKF(3,MDX,NDX,NR)=0d0
            CBFLDKF(1,MDX,NDX,NR)=0d0
            CBFLDKF(2,MDX,NDX,NR)=0d0
            CBFLDKF(3,MDX,NDX,NR)=0d0
         ENDDO
         ENDDO
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CEFLDD(1,NTH,NHH,NR)=0d0
            CEFLDD(2,NTH,NHH,NR)=0d0
            CEFLDD(3,NTH,NHH,NR)=0d0
            CBFLDD(1,NTH,NHH,NR)=0d0
            CBFLDD(2,NTH,NHH,NR)=0d0
            CBFLDD(3,NTH,NHH,NR)=0d0
            CEND(1,NTH,NHH,NR)  =0d0
            CEND(2,NTH,NHH,NR)  =0d0
            CEND(3,NTH,NHH,NR)  =0d0
            CEPD(1,NTH,NHH,NR)  =0d0
            CEPD(2,NTH,NHH,NR)  =0d0
            CEPD(3,NTH,NHH,NR)  =0d0
         ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
