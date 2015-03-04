C     $Id$
C
C     ****** CALCULATE ELECTRIC FIELD ******
C
      SUBROUTINE WMEFLD
C
      INCLUDE 'wmcomm.inc'
      DIMENSION CEF1(MDM,NDM),CEF2(MDM,NDM),RMA(3,3)
C
      DRHO1=(XRHO(2)-XRHO(1))**2
      DRHO2=(XRHO(3)-XRHO(1))**2
      A1= DRHO2/(DRHO2-DRHO1)
      A2=-DRHO1/(DRHO2-DRHO1)
C     
      DO NR=1,NRMAX
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
C
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

      DO NR=2,NRMAX+1
         NRP=MIN(NR+1,NRMAX+1)
C
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLD(1,MDX,NDX,NR)=0.5D0*(CEFLDK(1,MDX,NDX,NR )
     &                                +CEFLDK(1,MDX,NDX,NRP))
            CEFLD(2,MDX,NDX,NR)=CEFLDK(2,MDX,NDX,NR)
            CEFLD(3,MDX,NDX,NR)=CEFLDK(3,MDX,NDX,NR)
         ENDDO
         ENDDO
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
         ELSE
            CEFLD(3,MDX,NDX,NR)=0.D0
         ENDIF
      ENDDO
      ENDDO
C
C     ----- Inverse Fourier transform ----
C
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
      DO NR=1,2
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CEN(1,NTH,NHH,NR)=0.D0
         CEN(2,NTH,NHH,NR)=0.D0
         CEN(3,NTH,NHH,NR)=0.D0
         CEP(1,NTH,NHH,NR)=0.D0
         CEP(2,NTH,NHH,NR)=0.D0
         CEP(3,NTH,NHH,NR)=0.D0
      END DO
      END DO
      END DO


      DO NR=3,NRMAX+1
C
         XRI=1.D0/XRHO(NR)
         XRL=XRHO(NR)
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
C
C        ----- Calculate rotation matrix mu=RMA -----
C
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         TC2=BSUPTH/BABS
         TC3=BSUPPH/BABS
C
C        ***** RF11=RJ*SQRT(G^11)/XR *****
C
         RF11=SQRT(RG22(NTH,NHH,NR)*RG33(NTH,NHH,NR)
     &            -RG23(NTH,NHH,NR)*RG23(NTH,NHH,NR))
         RMA(1,1)= RJ(NTH,NHH,NR)/RF11*XRI
         RMA(2,1)= 0.D0
         RMA(3,1)= 0.D0
         RMA(1,2)= (TC2*(RG23(NTH,NHH,NR)*RG12(NTH,NHH,NR)
     &                  -RG22(NTH,NHH,NR)*RG13(NTH,NHH,NR))
     &             +TC3*(RG33(NTH,NHH,NR)*RG12(NTH,NHH,NR)
     &                  -RG23(NTH,NHH,NR)*RG13(NTH,NHH,NR))*XRI)
     &             /RF11
         RMA(2,2)= TC3*RF11*XRL
         RMA(3,2)=-TC2*RF11*XRL
         RMA(1,3)=TC2*RG12(NTH,NHH,NR)
     &           +TC3*RG13(NTH,NHH,NR)*XRI
         RMA(2,3)=TC2*RG22(NTH,NHH,NR)*XRL*XRL
     &           +TC3*RG23(NTH,NHH,NR)*XRL
         RMA(3,3)=TC2*RG23(NTH,NHH,NR)*XRL
     &           +TC3*RG33(NTH,NHH,NR)
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
         CEN(1,NTH,NHH,NR)=RMA(1,1)*CEFLD(1,NTH,NHH,NR)*XRI
     &                    +RMA(1,2)*CEFLD(2,NTH,NHH,NR)*XRL
     &                    +RMA(1,3)*CEFLD(3,NTH,NHH,NR)
         CEN(2,NTH,NHH,NR)=RMA(2,1)*CEFLD(1,NTH,NHH,NR)*XRI
     &                    +RMA(2,2)*CEFLD(2,NTH,NHH,NR)*XRL
     &                    +RMA(2,3)*CEFLD(3,NTH,NHH,NR)
         CEN(3,NTH,NHH,NR)=RMA(3,1)*CEFLD(1,NTH,NHH,NR)*XRI
     &                    +RMA(3,2)*CEFLD(2,NTH,NHH,NR)*XRL
     &                    +RMA(3,3)*CEFLD(3,NTH,NHH,NR)
         CEP(1,NTH,NHH,NR)=(   CEN(1,NTH,NHH,NR)
     &                     +CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
         CEP(2,NTH,NHH,NR)=(   CEN(1,NTH,NHH,NR)
     &                     -CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
         CEP(3,NTH,NHH,NR)=    CEN(3,NTH,NHH,NR)
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
                     CEP(K,NTH,NHH,NR)=(0.D0,0.D0)
                  ENDDO
               ENDDO
            ELSE
            ENDIF
         ENDDO
      ENDDO
C
      DO NHH=1,NHHMAX
         CEN1=0.D0
         CEN2=0.D0
         CEN3=0.D0
         CEP1=0.D0
         CEP2=0.D0
         CEP3=0.D0
         DO NTH=1,NTHMAX
            CEN1=CEN1+CEN(1,NTH,NHH,2)
            CEN2=CEN1+CEN(2,NTH,NHH,2)
            CEN3=CEN1+CEN(3,NTH,NHH,2)
            CEP1=CEP1+CEP(1,NTH,NHH,2)
            CEP2=CEP1+CEP(2,NTH,NHH,2)
            CEP3=CEP1+CEP(3,NTH,NHH,2)
         ENDDO
         CEN1=CEN1/NTHMAX
         CEN2=CEN2/NTHMAX
         CEN3=CEN3/NTHMAX
         CEP1=CEP1/NTHMAX
         CEP2=CEP2/NTHMAX
         CEP3=CEP3/NTHMAX
         DO NTH=1,NTHMAX
            CEN(1,NTH,NHH,1)=CEN1
            CEN(2,NTH,NHH,1)=CEN2
            CEN(3,NTH,NHH,1)=CEN3
            CEP(1,NTH,NHH,1)=CEP1
            CEP(2,NTH,NHH,1)=CEP2
            CEP(3,NTH,NHH,1)=CEP3
         ENDDO
      ENDDO
C      NR=1
C      NTH=1
C      NHH=1
C      WRITE(6,'(1P4E12.4)') CEN(1,NTH,NHH,NR),CEP(1,NTH,NHH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(2,NTH,NHH,NR),CEP(2,NTH,NHH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(3,NTH,NHH,NR),CEP(3,NTH,NHH,NR)
C      NR=2
C      NTH=1
C      NHH=1
C      WRITE(6,'(1P4E12.4)') CEN(1,NTH,NHH,NR),CEP(1,NTH,NHH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(2,NTH,NHH,NR),CEP(2,NTH,NHH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(3,NTH,NHH,NR),CEP(3,NTH,NHH,NR)
C
C     ----- Normalize CEFLD -----
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
C
         RF11=(RG22(NTH,NHH,NR)*RG33(NTH,NHH,NR)
     &        -RG23(NTH,NHH,NR)**2)/RJ(NTH,NHH,NR)**2
         RF22=(RG33(NTH,NHH,NR)*RG11(NTH,NHH,NR)
     &        -RG13(NTH,NHH,NR)**2)/RJ(NTH,NHH,NR)**2
         RF33=(RG11(NTH,NHH,NR)*RG22(NTH,NHH,NR)
     &         -RG12(NTH,NHH,NR)**2)/RJ(NTH,NHH,NR)**2
         RG011=SQRT(RF11)
         RG022=SQRT(RF22)
         RG033=SQRT(RF33)
C
         CEFLD(1,NTH,NHH,NR)=CEFLD(1,NTH,NHH,NR)*RG011
         CEFLD(2,NTH,NHH,NR)=CEFLD(2,NTH,NHH,NR)*RG022
         CEFLD(3,NTH,NHH,NR)=CEFLD(3,NTH,NHH,NR)*RG033
C
C         WRITE(21,*) nth,nph,nr,CEFLD(1,NTH,NHH,NR)
C         WRITE(21,*) nth,nph,nr,CEFLD(2,NTH,NHH,NR)
C         WRITE(21,*) nth,nph,nr,CEFLD(3,NTH,NHH,NR)
C
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE MAGNETIC FIELD ******
C
      SUBROUTINE WMBFLD
C
      INCLUDE 'wmcomm.inc'
      DIMENSION CBF1(MDM,NDM),CBF2(MDM,NDM)
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
C
         CBFLD(1,NTH,NHH,NR)=CBFLD(1,NTH,NHH,NR)/RJ(NTH,NHH,NR)
     &                      *SQRT(RG11(NTH,NHH,NR))
         CBFLD(2,NTH,NHH,NR)=CBFLD(2,NTH,NHH,NR)/RJ(NTH,NHH,NR)
     &                      *SQRT(RG22(NTH,NHH,NR))
         CBFLD(3,NTH,NHH,NR)=CBFLD(3,NTH,NHH,NR)/RJ(NTH,NHH,NR)
     &                      *SQRT(RG33(NTH,NHH,NR))
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
         DPH=2.D0*PI/DBLE(NHHMAX)
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
C
C      WRITE(6,'(A,1P2E12.4)') '# CPRAD=',CPRAD
C
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
      DIMENSION DS(NRM),DSS(NTHM,NHHM,NRM)
      DIMENSION CPF1(MDM,NDM),CPF2(MDM,NDM)
      DIMENSION CDV(3,3,3),CDW(3,3,3)
      DIMENSION CFA(NRM*NSM*MDM*MDM*NDM*NDM)
      DIMENSION CFB(NRM*NSM*MDM*MDM*NDM*NDM)
      DIMENSION iposa(nsize),ilena(nsize)
C
      NM=NRM*NSM*MDM*MDM*NDM*NDM
      CW=2.D0*PI*CRF*1.D6
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NHHMAX)
C
      IF(NRANK.EQ.0) THEN
         NRS=NBST
      ELSE
         NRS=NBST-1
      ENDIF
      IF(NRANK.EQ.NSIZE-1) THEN
         NRE=NBED+1
      ELSE
         NRE=NBED
      ENDIF
C         
      DO NR=1,NRMAX+1
      DO NS=1,NSMAX
      DO NKX=1,NDSIZ
      DO KDX=1,KDSIZ
      DO MLX=1,MDSIZ
      DO LDX=1,LDSIZ
         CPABSK(LDX,MLX,KDX,NKX,NS,NR)=(0.D0,0.D0)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSMAX
C
         CALL WMSETF(NRS,NS)
C
         DO NDX=1,NDSIZ
         DO KDX=1,KDSIZ
         DO MDX=1,MDSIZ
         DO LDX=1,LDSIZ
         DO J=1,3
         DO I=1,3
            CGD(I,J,LDX,MDX,KDX,NDX,2)=CGD(I,J,LDX,MDX,KDX,NDX,3)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
         DO NR=NRS,NBED
C
            CALL WMSETF(NR+1,NS)
C
            DO NDX=1,NDSIZ
            DO KDX=1,KDSIZ
            DO MDX=1,MDSIZ
            DO LDX=1,LDSIZ
            DO J=1,3
            DO I=1,3
               CGD(I,J,LDX,MDX,KDX,NDX,1)=CGD(I,J,LDX,MDX,KDX,NDX,2)
               CGD(I,J,LDX,MDX,KDX,NDX,2)=CGD(I,J,LDX,MDX,KDX,NDX,3)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
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
            DO KK=KDMIN,KDMAX
               KKX=KK-KDMIN+1
               KD=NK-ND
               KX1=MOD( KD+KK-KDMIN+4*KDSIZ,KDSIZ)+1
               NX1=MOD( NK   -NDMIN+4*NDSIZ,NDSIZ)+1
C
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
            DO ML=MDMIN,MDMAX
               MLX=ML-MDMIN+1
            DO LL=LDMIN,LDMAX
               LLX=LL-LDMIN+1
               LD= ML-MD
               LX1=MOD( LD+LL-LDMIN+4*LDSIZ,LDSIZ)+1
               MX1=MOD( ML   -MDMIN+4*MDSIZ,MDSIZ)+1
C
               DO K=1,2
               DO J=1,3
               DO I=1,3
                  CDV(I,J,K)=CGD(I,J,LX1,MX1,KX1,NX1,K)
                  CDW(I,J,K)=CGD(I,J,LX1,MX1,KX1,NX1,K)
               ENDDO
               ENDDO
               ENDDO
C
               FACT1M=XRHOMH/XRHOM
               FACT1C=XRHOMH/XRHOC
C               FACT1P=XRHOPH/XRHOC
C               FACT2M=XRHOMH/XRHOM
               FACT2C=XRHOMH/XRHOC
               FACT2P=XRHOPH/XRHOC
C               FACT3M=XRHOMH/XRHOM
               FACT3C=XRHOMH/XRHOC
               FACT3P=XRHOPH/XRHOC
C
C               CDV11M=CHERMIT(CDV(1,1,1),CDW(1,1,1))
C               CDV11C=CHERMIT(CDV(1,1,2),CDW(1,1,2))
C               CDV12M=CHERMIT(CDV(1,2,1),CDW(2,1,1))
C               CDV12C=CHERMIT(CDV(1,2,2),CDW(2,1,2))
C               CDV13M=CHERMIT(CDV(1,3,1),CDW(3,1,1))
C               CDV13C=CHERMIT(CDV(1,3,2),CDW(3,1,2))
C               CDV21C=CHERMIT(CDV(2,1,2),CDW(1,2,2))
C               CDV22C=CHERMIT(CDV(2,2,2),CDW(2,2,2))
C               CDV23C=CHERMIT(CDV(2,3,2),CDW(3,2,2))
C               CDV31C=CHERMIT(CDV(3,1,2),CDW(1,3,2))
C               CDV32C=CHERMIT(CDV(3,2,2),CDW(2,3,2))
C               CDV33C=CHERMIT(CDV(3,3,2),CDW(3,3,2))
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
C
               CCE1=DCONJG(CEFLDK(1,MLX,NKX,NR+1))
               CCE2=DCONJG(CEFLDK(2,MLX,NKX,NR+1))
               CCE3=DCONJG(CEFLDK(3,MLX,NKX,NR+1))
C
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
     &            *(0
     &             +CPM11
     &             +CPM12
     &             +CPM13
     &             )*DPSIPDRHOMH*DRHOM
C
               CPABSC=CI*EPS0*(VC*VC/CW)
     &            *(0
     &             +CPC11
     &             +CPC12
     &             +CPC13
     &             )*DPSIPDRHOMH*DRHOM
     &            +CI*EPS0*(VC*VC/CW)
     &            *(
     &             +CPC22
     &             +CPC33
     &             +CPC23
     &             +CPC32
     &             +CPC21
     &             +CPC31
     &             +CPP21
     &             +CPP31
     &             )*DPSIPDRHOC*DRHOPM
C
               CPABSK(LLX,MDX,KKX,NDX,NS,NR)
     &        =CPABSK(LLX,MDX,KKX,NDX,NS,NR)
     &        +0.5D0*CPABSM
C               CPABSK(LLX,MLX,KKX,NKX,NS,NR)
C     &        =CPABSK(LLX,MLX,KKX,NKX,NS,NR)
C     &        +0.5D0*CPABSM
C
               CPABSK(LLX,MDX,KKX,NDX,NS,NR+1)
     &        =CPABSK(LLX,MDX,KKX,NDX,NS,NR+1)
     &        +0.5D0*CPABSC
C               CPABSK(LLX,MLX,KKX,NKX,NS,NR+1)
C     &        =CPABSK(LLX,MLX,KKX,NKX,NS,NR+1)
C     &        +0.5D0*CPABSC
C
            ENDDO
            ENDDO
            ENDDO
C
            ENDDO
            ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      NRS=NBST
      IF(NRANK.EQ.NSIZE-1) THEN
         NRE=NBED+1
      ELSE
         NRE=NBED
      ENDIF
C
      MN=0
      DO NR=NRS,NRE
      DO NS=1,NSMAX
      DO NKX=1,NDSIZ
      DO KDX=1,KDSIZ
      DO MLX=1,MDSIZ
      DO LDX=1,LDSIZ
         MN=MN+1
         CFB(MN)=CPABSK(LDX,MLX,KDX,NKX,NS,NR)
C         if(nrank.eq.0) write(23,*) nr,ns,CPABSK(MLX,LDX,NKX,KDX,NS,NR)
C         if(nrank.eq.1) write(24,*) nr,ns,CPABSK(MLX,LDX,NKX,KDX,NS,NR)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      CALL mtx_allgather1_integer(MN,ilena)
      ntot=0
      DO i=1,nsize
         iposa(i)=ntot
         ntot=ntot+ilena(i)
      END DO
      CALL mtx_gatherv_complex8(CFB,MN,CFA,ntot,ilena,iposa)
      NM=ntot
C      CALL MPGTCV(CFB,MN,CFA,NVTOT,NM)
C
      IF(NRANK.EQ.0) THEN
         MN=0
         DO NR=1,NRMAX
         DO NS=1,NSMAX
         DO NKX=1,NDSIZ
         DO KDX=1,KDSIZ
         DO MLX=1,MDSIZ
         DO LDX=1,LDSIZ
            MN=MN+1
            CPABSK(LDX,MLX,KDX,NKX,NS,NR)=CFA(MN)
C            if(nrank.eq.0) write(21,*) nr,ns,
C     &           CPABSK(MLX,LDX,NKX,KDX,NS,NR)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE PABS IN MODE NUMBER SPACE +++++
C
         DO NR=1,NRMAX+1
         DO NS=1,NSMAX
         DO NDX=1,NDSIZ
            KK=0
            KKX=KK-KDMIN+1
         DO MDX=1,MDSIZ
            LL=0
            LLX=LL-LDMIN+1
            PABSK(MDX,NDX,NR,NS)=DBLE(CPABSK(LLX,MDX,KKX,NDX,NS,NR))
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE PABS IN REAL SPACE +++++
C
         DO NS=1,NSMAX
         DO NR=1,NRMAX+1
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NHH,NR,NS)=0.D0
               CPABS(NTH,NHH,NR,NS)=0.D0
            ENDDO
            ENDDO
            DO NDX=1,NDSIZ
            DO MDX=1,MDSIZ
               DO KK=KDMIN,KDMAX
                  KKX=KK-KDMIN+1
               DO LL=LDMIN,LDMAX
                  LLX=LL-LDMIN+1
                  CPF1(LLX,KKX)=CPABSK(LLX,MDX,KKX,NDX,NS,NR)
               ENDDO
               ENDDO
               CALL WMSUBE(CPF1,CPF2)
               DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PABS(NTH,NHH,NR,NS)=PABS(NTH,NHH,NR,NS)
     &                               +DBLE(CPF2(NTH,NHH))
                  CPABS(NTH,NHH,NR,NS)=CPABS(NTH,NHH,NR,NS)
     &                               +CPF2(NTH,NHH)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE DRIVEN CURRENT IN REAL SPACE +++++
C
         DO NR=1,NRMAX+1
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NHH,NR)=0.D0
         ENDDO
         ENDDO
         ENDDO

      IF(NSMAX.GE.1) THEN
      NS=1
      DO NR=1,NRMAX+1
         CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
         VTE=SQRT(RTPR(1)*AEE*1.D3/(PA(1)*AMP))
         WW=DBLE(CW)
         IF(RN(1).LE.0.D0) THEN
            RLNLMD=15.D0
         ELSE
            RT=(RTPR(1)+2*RTPP(1))/3.D0
            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(1))
     &           + 2.30D0*LOG10(RT)
         ENDIF
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            MM=NTH0+MD
            DO KKX=1,KDSIZ
            DO LLX=1,LDSIZ
               CPF1(LLX,KKX)=CPABSK(LLX,MDX,KKX,NDX,NS,NR)
            ENDDO
            ENDDO
            CALL WMSUBE(CPF1,CPF2)
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
               RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
              IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
              IF(ABS(WW/RKPR).LT.VC) THEN
                 W=WW/(RKPR*VTE)
                 XL=(RPST(NTH,NHH,NR)-RR  )/RR
                 YL=(ZPST(NTH,NHH,NR)-0.D0)/RR
                 EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                 IF(W.LT.0.D0) EFCD=-EFCD
                 IF (RN(1).GT.0.D0) THEN
                    PCUR(NTH,NHH,NR)=PCUR(NTH,NHH,NR)
     &                   +0.384D0*RTPR(1)/(AEE*1.D3)*EFCD
     &                   /((RN(1)/1.D20)*RLNLMD)*DBLE(CPF2(NTH,NHH))
     &                   /(2.D0*PI*RPST(NTH,NHH,NR))
                 END IF
              ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDIF
C
      DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH
         ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         PABSR(NR,NS)=0.D0
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DPH
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
      FACT=1.D0
C
      IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
         FACT=PRFIN/PABSTT
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
            IF(MODELG.EQ.3) THEN
               DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
            ELSE
               DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
            ENDIF
            DSSS=DPSIPDRHO*DRHO*RJ(NTH,NHH,NR)
            DSS(NTH,NHH,NR)=1.D0/DSSS
            DS(NR)=DS(NR)+DSSS*DTH*DPH
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
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            PABSK(MDX,NDX,NR,NS)=FACT*PABSK(MDX,NDX,NR,NS)*DS(NR)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      PCURT=FACT*PCURT
      DO NR=1,NRMAX
         PCURR(NR)=FACT*PCURR(NR)*DS(NR)
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NHH,NR)=FACT*PCUR(NTH,NHH,NR)*DSS(NTH,NHH,NR)
         ENDDO
         ENDDO
      ENDDO
      ENDIF
C
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
      CALL MPSYNC
C
      RETURN
      END
