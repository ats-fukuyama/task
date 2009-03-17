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
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CEFLD(I,NTH,NPH,NR)=CEF2(NTH,NPH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- Calculate CEN from CEsup -----
C
      DO NR=2,NRMAX+1
C
         XRI=1.D0/XRHO(NR)
         XRL=XRHO(NR)
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C
C        ----- Calculate rotation matrix mu=RMA -----
C
         CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
         TC2=BSUPTH/BABS
         TC3=BSUPPH/BABS
C
C        ***** RF11=RJ*SQRT(G^11)/XR *****
C
         RF11=SQRT(RG22(NTH,NPH,NR)*RG33(NTH,NPH,NR)
     &            -RG23(NTH,NPH,NR)*RG23(NTH,NPH,NR))
         RMA(1,1)= RJ(NTH,NPH,NR)/RF11*XRI
         RMA(2,1)= 0.D0
         RMA(3,1)= 0.D0
         RMA(1,2)= (TC2*(RG23(NTH,NPH,NR)*RG12(NTH,NPH,NR)
     &                  -RG22(NTH,NPH,NR)*RG13(NTH,NPH,NR))
     &             +TC3*(RG33(NTH,NPH,NR)*RG12(NTH,NPH,NR)
     &                  -RG23(NTH,NPH,NR)*RG13(NTH,NPH,NR))*XRI)
     &             /RF11
         RMA(2,2)= TC3*RF11*XRL
         RMA(3,2)=-TC2*RF11*XRL
         RMA(1,3)=TC2*RG12(NTH,NPH,NR)
     &           +TC3*RG13(NTH,NPH,NR)*XRI
         RMA(2,3)=TC2*RG22(NTH,NPH,NR)*XRL*XRL
     &           +TC3*RG23(NTH,NPH,NR)*XRL
         RMA(3,3)=TC2*RG23(NTH,NPH,NR)*XRL
     &           +TC3*RG33(NTH,NPH,NR)
C
         CALL INVMRD(RMA,3,3,ILL)
         IF(ILL.NE.0) THEN
            WRITE(6,*) 'XX WEFLD: INVMRD(RMA) : SINGULAR MATRIX'
            WRITE(6,'(3I5,1P2E12.4)') NR,NTH,NPH,TC3,TC2
         ENDIF
C
C         IF(NR.EQ.2.OR.NR.EQ.3) THEN
C            WRITE(6,*) 'NR,NTH,NPH=',NR,NTH,NPH
C            WRITE(6,'(1P3E12.4)') RMA(1,1),RMA(1,2),RMA(1,3)
C            WRITE(6,'(1P3E12.4)') RMA(2,1),RMA(2,2),RMA(2,3)
C            WRITE(6,'(1P3E12.4)') RMA(3,1),RMA(3,2),RMA(3,3)
C            WRITE(6,'(1P3E12.4)') RJ(NTH,NPH,NR),RF11,XRI
C            WRITE(6,'(1P3E12.4)') TC2,TC3,XRL
C         ENDIF
C
         CEN(1,NTH,NPH,NR)=RMA(1,1)*CEFLD(1,NTH,NPH,NR)*XRI
     &                    +RMA(1,2)*CEFLD(2,NTH,NPH,NR)*XRL
     &                    +RMA(1,3)*CEFLD(3,NTH,NPH,NR)
         CEN(2,NTH,NPH,NR)=RMA(2,1)*CEFLD(1,NTH,NPH,NR)*XRI
     &                    +RMA(2,2)*CEFLD(2,NTH,NPH,NR)*XRL
     &                    +RMA(2,3)*CEFLD(3,NTH,NPH,NR)
         CEN(3,NTH,NPH,NR)=RMA(3,1)*CEFLD(1,NTH,NPH,NR)*XRI
     &                    +RMA(3,2)*CEFLD(2,NTH,NPH,NR)*XRL
     &                    +RMA(3,3)*CEFLD(3,NTH,NPH,NR)
         CEP(1,NTH,NPH,NR)=(   CEN(1,NTH,NPH,NR)
     &                     +CI*CEN(2,NTH,NPH,NR))/SQRT(2.D0)
         CEP(2,NTH,NPH,NR)=(   CEN(1,NTH,NPH,NR)
     &                     -CI*CEN(2,NTH,NPH,NR))/SQRT(2.D0)
         CEP(3,NTH,NPH,NR)=    CEN(3,NTH,NPH,NR)
C
C         IF(NR.EQ.1) THEN
C            WRITE(6,'(1P4E12.4)') CEN(1,NTH,NPH,NR),CEP(1,NTH,NPH,NR)
C            WRITE(6,'(1P4E12.4)') CEN(2,NTH,NPH,NR),CEP(2,NTH,NPH,NR)
C            WRITE(6,'(1P4E12.4)') CEN(3,NTH,NPH,NR),CEP(3,NTH,NPH,NR)
C         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
C     ------- プラズマの外側のE+,E-,Ez0にする ---------
      DO K=1,3
         DO NR=1,NRMAX+1
            IF(XRHO(NR).GT.1.0D0) THEN
               DO NTH=1,NTHMAX
                  DO NPH=1,NPHMAX
                     CEP(K,NTH,NPH,NR)=(0.D0,0.D0)
                  ENDDO
               ENDDO
            ELSE
            ENDIF
         ENDDO
      ENDDO
C
      DO NPH=1,NPHMAX
         CEN1=0.D0
         CEN2=0.D0
         CEN3=0.D0
         CEP1=0.D0
         CEP2=0.D0
         CEP3=0.D0
         DO NTH=1,NTHMAX
            CEN1=CEN1+CEN(1,NTH,NPH,2)
            CEN2=CEN1+CEN(2,NTH,NPH,2)
            CEN3=CEN1+CEN(3,NTH,NPH,2)
            CEP1=CEP1+CEP(1,NTH,NPH,2)
            CEP2=CEP1+CEP(2,NTH,NPH,2)
            CEP3=CEP1+CEP(3,NTH,NPH,2)
         ENDDO
         CEN1=CEN1/NTHMAX
         CEN2=CEN2/NTHMAX
         CEN3=CEN3/NTHMAX
         CEP1=CEP1/NTHMAX
         CEP2=CEP2/NTHMAX
         CEP3=CEP3/NTHMAX
         DO NTH=1,NTHMAX
            CEN(1,NTH,NPH,1)=CEN1
            CEN(2,NTH,NPH,1)=CEN2
            CEN(3,NTH,NPH,1)=CEN3
            CEP(1,NTH,NPH,1)=CEP1
            CEP(2,NTH,NPH,1)=CEP2
            CEP(3,NTH,NPH,1)=CEP3
         ENDDO
      ENDDO
C      NR=1
C      NTH=1
C      NPH=1
C      WRITE(6,'(1P4E12.4)') CEN(1,NTH,NPH,NR),CEP(1,NTH,NPH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(2,NTH,NPH,NR),CEP(2,NTH,NPH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(3,NTH,NPH,NR),CEP(3,NTH,NPH,NR)
C      NR=2
C      NTH=1
C      NPH=1
C      WRITE(6,'(1P4E12.4)') CEN(1,NTH,NPH,NR),CEP(1,NTH,NPH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(2,NTH,NPH,NR),CEP(2,NTH,NPH,NR)
C      WRITE(6,'(1P4E12.4)') CEN(3,NTH,NPH,NR),CEP(3,NTH,NPH,NR)
C
C     ----- Normalize CEFLD -----
C
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C
         RF11=(RG22(NTH,NPH,NR)*RG33(NTH,NPH,NR)
     &        -RG23(NTH,NPH,NR)**2)/RJ(NTH,NPH,NR)**2
         RF22=(RG33(NTH,NPH,NR)*RG11(NTH,NPH,NR)
     &        -RG13(NTH,NPH,NR)**2)/RJ(NTH,NPH,NR)**2
         RF33=(RG11(NTH,NPH,NR)*RG22(NTH,NPH,NR)
     &         -RG12(NTH,NPH,NR)**2)/RJ(NTH,NPH,NR)**2
         RG011=SQRT(RF11)
         RG022=SQRT(RF22)
         RG033=SQRT(RF33)
C
         CEFLD(1,NTH,NPH,NR)=CEFLD(1,NTH,NPH,NR)*RG011
         CEFLD(2,NTH,NPH,NR)=CEFLD(2,NTH,NPH,NR)*RG022
         CEFLD(3,NTH,NPH,NR)=CEFLD(3,NTH,NPH,NR)*RG033
C
C         WRITE(21,*) nth,nph,nr,CEFLD(1,NTH,NPH,NR)
C         WRITE(21,*) nth,nph,nr,CEFLD(2,NTH,NPH,NR)
C         WRITE(21,*) nth,nph,nr,CEFLD(3,NTH,NPH,NR)
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
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CBFLD(I,NTH,NPH,NR)=CBF2(NTH,NPH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C
         CBFLD(1,NTH,NPH,NR)=CBFLD(1,NTH,NPH,NR)/RJ(NTH,NPH,NR)
     &                      *SQRT(RG11(NTH,NPH,NR))
         CBFLD(2,NTH,NPH,NR)=CBFLD(2,NTH,NPH,NR)/RJ(NTH,NPH,NR)
     &                      *SQRT(RG22(NTH,NPH,NR))
         CBFLD(3,NTH,NPH,NR)=CBFLD(3,NTH,NPH,NR)/RJ(NTH,NPH,NR)
     &                      *SQRT(RG33(NTH,NPH,NR))
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE ABSORBED POWER ******
C
      SUBROUTINE WMPABS
C
      INCLUDE 'wmcomm.inc'
C
      IF(MDLWMX.EQ.0) THEN
         CALL WMPABS0
      ELSE
         CALL WMPABS1
      ENDIF
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
         DPH=2.D0*PI/DBLE(NPHMAX)
         CPRAD=(0.D0,0.D0)
C
         DO NR=NRANT-1,NRMAX
C
            IF(MDLWMX.EQ.0) THEN
               CALL WMSETM0_V(NR,CFVP)
            ELSE
               CALL WMSETM1_V(NR,CFVP)
            ENDIF
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
      IF(MYRANK.EQ.0) THEN
         WRITE(6,601) CRADTT,PABSTT,PCURT
         WRITE(6,602) (PABST(NS),NS=1,NSMAX)
      ENDIF
C
      IF(NPRINT.LT.2) RETURN
C
      IF(MYRANK.EQ.0) THEN
C         WRITE(6,*) '   MD   IMP                        PABSKT'
C         DO ND=NDMIN,NDMAX
C            NDX=ND-NDMIN+1
C            NN=NPH0+ND
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
         IF(MYRANK.EQ.0) 
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
      DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=RF1(NTH,NPH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF2(LDX,NPH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NPH=1,NPHMAX
            CFN(NPH)=CF2(LDX,NPH)
         ENDDO
         CALL WMXFFT(CFN,NPHMAX,0)
         DO KDX=1,KDSIZ
            CF2(LDX,KDX)=CFN(KDX)
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
         CALL WMXFFT(CFN,NPHMAX,1)
         DO NPH=1,NPHMAX
            CF2(NTH,NPH)=CFN(NPH)
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
      DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=1.D0/RF1(NTH,NPH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF2(LDX,NPH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NPH=1,NPHMAX
            CFN(NPH)=CF2(LDX,NPH)
         ENDDO
         CALL WMXFFT(CFN,NPHMAX,0)
         DO KDX=1,KDSIZ
            CF2(LDX,KDX)=CFN(KDX)
         ENDDO
      ENDDO
      RETURN
      END
