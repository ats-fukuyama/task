! wmemfp.f90

MODULE wmemfp

  PRIVATE
  PUBLIC wm_efield,wm_bfield,wm_pwrflux,wm_pwrant,wm_pabs

CONTAINS
  
!     ****** CALCULATE ELECTRIC FIELD ******

  SUBROUTINE wm_efield

    USE wmcomm
    USE wmprof
    USE wmsub
    IMPLICIT NONE
    COMPLEX(rkind),ALLOCATABLE:: CEF1(:,:),CEF2(:,:),CEFM1(:,:),CEFM2(:,:)
    COMPLEX(rkind),ALLOCATABLE:: CEFLDR(:,:),CRMARHF(:,:),CEFLDM(:,:,:,:)
    COMPLEX(rkind):: CEN1,CEN2,CEN3,CEP1,CEP2,CEP3
    REAL(rkind):: RMA(3,3)
    REAL(rkind):: DRHO1,DRHO2,A1,A2,XRI,XRL,BABS,BSUPTH,BSUPPH
    REAL(rkind):: TC2,TC3,RF11,RF22,RF33,RG011,RG022,RG033
    INTEGER:: NR,ND,MD,NDX,MDX,IG,NRP,I,NDXM,NPH,NPP,NTH,NHH,NTHF,NHHF
    INTEGER:: ILL,K,NHHD,MM

    ALLOCATE(CEF1(nthmax,nhhmax),CEF2(nthmax,nhhmax))
    ALLOCATE(CEFM1(nthmax,nhhmax),CEFM2(nthmax,nhhmax))
    ALLOCATE(CEFLDR(nthmax,nhhmax), CRMARHF(nthmax,nhhmax))
    ALLOCATE(CEFLDM(3,nthmax,nhhmax,nrmax+1))

    DRHO1=(XRHO(2)-XRHO(1))**2
    DRHO2=(XRHO(3)-XRHO(1))**2
    A1= DRHO2/(DRHO2-DRHO1)
    A2=-DRHO1/(DRHO2-DRHO1)

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

    NR=1
    DO ND=NDMIN,NDMAX
       NDX=ND-NDMIN+1
       DO MD=MDMIN,MDMAX
          MDX=MD-MDMIN+1
          MM=NTH0+MD
          IF(ABS(MM).EQ.1) THEN
             CEFLDK(1,MDX,NDX,NR)=CEFLDK(1,MDX,NDX,NR+1)
             CEFLDK(2,MDX,NDX,NR)=A1*CEFLDK(2,MDX,NDX,NR+1) &
                                 +A2*CEFLDK(2,MDX,NDX,NR+2)
          ELSE
             CEFLDK(1,MDX,NDX,NR)=0.D0
             CEFLDK(2,MDX,NDX,NR)=0.D0
          ENDIF
          IF(MM.EQ.0) THEN
             CEFLDK(3,MDX,NDX,NR)=A1*CEFLDK(3,MDX,NDX,NR+1) &
                                 +A2*CEFLDK(3,MDX,NDX,NR+2)
          ELSE
             CEFLDK(3,MDX,NDX,NR)=0.D0
          ENDIF
       ENDDO
    ENDDO

    DO NR=1,NRMAX
       NRP=MIN(NR+1,NRMAX)
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             CEFLD(1,MDX,NDX,NR)=0.5D0*(CEFLDK(1,MDX,NDX,NR ) &
                                       +CEFLDK(1,MDX,NDX,NRP))
             CEFLD(2,MDX,NDX,NR)=CEFLDK(2,MDX,NDX,NR)
             CEFLD(3,MDX,NDX,NR)=CEFLDK(3,MDX,NDX,NR)
          ENDDO
       ENDDO
    END DO

    NR=1
    DO ND=NDMIN,NDMAX
       NDX=ND-NDMIN+1
       DO MD=MDMIN,MDMAX
          MDX=MD-MDMIN+1
          MM=NTH0+MD
          IF(ABS(MM).EQ.1) THEN
             CEFLD(1,MDX,NDX,NR)=A1*CEFLD(1,MDX,NDX,NR+1) &
                                 +A2*CEFLD(1,MDX,NDX,NR+2)
             CEFLD(2,MDX,NDX,NR)=A1*CEFLD(2,MDX,NDX,NR+1) &
                                 +A2*CEFLD(2,MDX,NDX,NR+2)
          ELSE
             CEFLD(1,MDX,NDX,NR)=0.D0
             CEFLD(2,MDX,NDX,NR)=0.D0
          ENDIF
          IF(MM.EQ.0) THEN
             CEFLD(3,MDX,NDX,NR)=A1*CEFLD(3,MDX,NDX,NR+1) &
                                +A2*CEFLD(3,MDX,NDX,NR+2)
             CEFLD(3,MDX,NDX,NR)=CEFLD(3,MDX,NDX,NR+1)
          ELSE
             CEFLD(3,MDX,NDX,NR)=0.D0
          ENDIF
       ENDDO
    ENDDO

!    ----- Inverse Fourier transform ----

    DO NR=1,NRMAX+1
       DO I=1,3
          CEFM1(1:nthmax,1:nhhmax)=(0.D0,0.D0)
          DO NDX=1,NDSIZ
             DO MDX=1,MDSIZ
                CEFM1(MDX,NDX)=CEFLD(I,MDX,NDX,NR)
             ENDDO
          ENDDO
          CALL WMSUBE(CEFM1,CEFM2)
          DO NHH=1,NDSIZ
             DO NTH=1,NTHMAX
                CEFLDM(I,NTH,NHH,NR)=CEFM2(NTH,NHH)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

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

!     ----- Calculate CEN from CEsup -----

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

    DO NR=1,NRMAX+1

       IF(NR.EQ.1) THEN
          XRI=1.D6/XRHO(2)
          XRL=XRHO(2)/1.D6
       ELSE
          XRI=1.D0/XRHO(NR)
          XRL=XRHO(NR)
       ENDIF

       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             NHHF=(NHH-1)*FACTOR_NHH +1
             NTHF=(NTH-1)*FACTOR_NTH +1

!        ----- Calculate rotation matrix mu=RMA -----

             CALL WMCMAG(NR,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
             TC2=BSUPTH/BABS
             TC3=BSUPPH/BABS

!        ***** RF11=RJ*SQRT(G^11)/XR *****

             RF11=SQRT(RG22(NTHF,NHHF,NR)*RG33(NTHF,NHHF,NR) &
                      -RG23(NTHF,NHHF,NR)*RG23(NTHF,NHHF,NR))
             RMA(1,1)= RJ(NTHF,NHHF,NR)/RF11*XRI
             RMA(2,1)= 0.D0
             RMA(3,1)= 0.D0
             RMA(1,2)= (TC2*(RG23(NTHF,NHHF,NR)*RG12(NTHF,NHHF,NR) &
                            -RG22(NTHF,NHHF,NR)*RG13(NTHF,NHHF,NR)) &
                       +TC3*(RG33(NTHF,NHHF,NR)*RG12(NTHF,NHHF,NR) &
                            -RG23(NTHF,NHHF,NR)*RG13(NTHF,NHHF,NR))*XRI) &
                       /RF11
             RMA(2,2)= TC3*RF11*XRL
             RMA(3,2)=-TC2*RF11*XRL
             RMA(1,3)=TC2*RG12(NTHF,NHHF,NR) &
                     +TC3*RG13(NTHF,NHHF,NR)*XRI
             RMA(2,3)=TC2*RG22(NTHF,NHHF,NR)*XRL*XRL &
                     +TC3*RG23(NTHF,NHHF,NR)*XRL
             RMA(3,3)=TC2*RG23(NTHF,NHHF,NR)*XRL &
                     +TC3*RG33(NTHF,NHHF,NR)

             CALL INVMRD(RMA,3,3,ILL)
             IF(ILL.NE.0) THEN
                WRITE(6,*) 'XX WEFLD: INVMRD(RMA) : SINGULAR MATRIX'
                WRITE(6,'(3I5,1P2E12.4)') NR,NTH,NHH,TC3,TC2
                WRITE(6,'(15X,3ES12.4)') (RMA(1,I),I=1,3)
                WRITE(6,'(15X,3ES12.4)') (RMA(2,I),I=1,3)
                WRITE(6,'(15X,3ES12.4)') (RMA(3,I),I=1,3)
                STOP
             ENDIF

             CEN(1,NTH,NHH,NR) =CEFLD(1,NTH,NHH,NR)
             CEN(2,NTH,NHH,NR) =CEFLD(2,NTH,NHH,NR)
             CEN(3,NTH,NHH,NR) =CEFLD(3,NTH,NHH,NR)

             CEP(1,NTH,NHH,NR) =(   CEN(1,NTH,NHH,NR) &
                               + CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
             CEP(2,NTH,NHH,NR) =(   CEN(1,NTH,NHH,NR) &
                               - CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
             CEP(3,NTH,NHH,NR) =    CEN(3,NTH,NHH,NR)
          ENDDO
       ENDDO
    ENDDO

!     ------- Set E+, E- and Ez0 to those outside the plasma ---------

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

!     ----- Normalize CEFLD -----

    DO NR=1,NRMAX+1
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             NHHF=(NHH-1)*FACTOR_NHH+1
             NTHF=(NTH-1)*FACTOR_NTH +1
 
             RF11=(RG22(NTHF,NHHF,NR)*RG33(NTHF,NHHF,NR) &
                  -RG23(NTHF,NHHF,NR)**2)/RJ(NTHF,NHHF,NR)**2
             RF22=(RG33(NTHF,NHHF,NR)*RG11(NTHF,NHHF,NR) &
                  -RG13(NTHF,NHHF,NR)**2)/RJ(NTHF,NHHF,NR)**2
             RF33=(RG11(NTHF,NHHF,NR)*RG22(NTHF,NHHF,NR) &
                  -RG12(NTHF,NHHF,NR)**2)/RJ(NTHF,NHHF,NR)**2
             RG011=SQRT(RF11)
             RG022=SQRT(RF22)
             RG033=SQRT(RF33)

             CEFLD(1,NTH,NHH,NR) =CEFLD(1,NTH,NHH,NR)
             CEFLD(2,NTH,NHH,NR) =CEFLD(2,NTH,NHH,NR)
             CEFLD(3,NTH,NHH,NR) =CEFLD(3,NTH,NHH,NR)
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE wm_efield

!     ****** CALCULATE MAGNETIC FIELD ******

  SUBROUTINE wm_bfield
    USE wmcomm
    USE wmsub
    IMPLICIT NONE
    COMPLEX(rkind),ALLOCATABLE:: CBF1(:,:),CBF2(:,:)
    COMPLEX(rkind):: CW,CBFLD1,CBFLD2,CBFLD3
    REAL(rkind):: DRHO1,DRHO2,A1,A2,XRHOM,XRHOMH,XRHOC,QPMH,DPSIPDRHOMH
    REAL(rkind):: DRHOM
    INTEGER:: NR,ND,NDX,NN,MD,MDX,MM,NDXM,I,NTH,NHH,NTHF,NHHF

    ALLOCATE(CBF1(nthmax,nhhmax),CBF2(nthmax,nhhmax))

    CW=2.D0*PI*CMPLX(RF,RFI)*1.D6

    DRHO1=(XRHO(2)-XRHO(1))**2
    DRHO2=(XRHO(3)-XRHO(1))**2
    A1= DRHO2/(DRHO2-DRHO1)
    A2=-DRHO1/(DRHO2-DRHO1)

    DO NR=2,NRMAX+1
       XRHOM=       XRHO(NR-1)
       XRHOMH=0.5D0*(XRHO(NR-1)+XRHO(NR))
       XRHOC=       XRHO(NR)
       IF(MODELG.EQ.3) THEN
          QPMH=0.5D0*(QPS(NR-1)+QPS(NR))
          DPSIPDRHOMH=2.D0*PSITA*XRHOMH/QPMH
       ELSE
          DPSIPDRHOMH=2.D0*PSIPA*XRHOMH
       ENDIF

       DRHOM=XRHO(NR)-XRHO(NR-1)

       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          NN=NPH0+NHC*ND
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             MM=NTH0+MD
             CBFLD1=-(CI/CW) &
                   *(CI*MM*CEFLDK(3,MDX,NDX,NR) &
                    -CI*NN*CEFLDK(2,MDX,NDX,NR))/XRHOC
             CBFLD2=-(CI/CW) &
                   *(CI*NN*CEFLDK(1,MDX,NDX,NR) &
                    -(CEFLDK(3,MDX,NDX,NR  ) &
                    -CEFLDK(3,MDX,NDX,NR-1))/(DPSIPDRHOMH*DRHOM) &
                    )*XRHOM
             CBFLD3=-(CI/CW) &
                   *((CEFLDK(2,MDX,NDX,NR  ) &
                     -CEFLDK(2,MDX,NDX,NR-1))/(DPSIPDRHOMH*DRHOM) &
                    -CI*MM*CEFLDK(1,MDX,NDX,NR))

             CBFLDK(1,MDX,NDX,NR)=VC*CBFLD1
             CBFLDK(2,MDX,NDX,NR)=VC*CBFLD2
             CBFLDK(3,MDX,NDX,NR)=VC*CBFLD3

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
             CBFLDK(1,MDX,NDX,NR)=A1*CBFLDK(1,MDX,NDX,NR+1) &
                                 +A2*CBFLDK(1,MDX,NDX,NR+2)
          ELSE
             CBFLDK(1,MDX,NDX,NR)=0.D0
          ENDIF
          CBFLDK(2,MDX,NDX,NR)=CBFLDK(2,MDX,NDX,2)
          CBFLDK(3,MDX,NDX,NR)=CBFLDK(3,MDX,NDX,2)
       ENDDO
    ENDDO

    DO NR=1,NRMAX+1
       DO I=1,3
          DO NDX=1,NDSIZ
             DO MDX=1,MDSIZ
                CBF1(MDX,NDX)=CBFLDK(I,MDX,NDX,NR)
             ENDDO
          ENDDO
          CALL WMSUBE(CBF1,CBF2)
          DO NHH=1,NDSIZ
             DO NTH=1,NTHMAX
                CBFLD(I,NTH,NHH,NR)=CBF2(NTH,NHH)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO NR=1,NRMAX+1
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             NHHF=(NHH-1)*FACTOR_NHH +1
             NTHF=(NTH-1)*FACTOR_NTH +1
             CBFLD(1,NTH,NHH,NR)=CBFLD(1,NTH,NHH,NR)/RJ(NTHF,NHHF,NR) &
                                *SQRT(RG11(NTHF,NHHF,NR))
             CBFLD(2,NTH,NHH,NR)=CBFLD(2,NTH,NHH,NR)/RJ(NTHF,NHHF,NR) &
                                *SQRT(RG22(NTHF,NHHF,NR))
             CBFLD(3,NTH,NHH,NR)=CBFLD(3,NTH,NHH,NR)/RJ(NTHF,NHHF,NR) &
                                *SQRT(RG33(NTHF,NHHF,NR))
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE wm_bfield

!     ****** HERMITIAN ******

  FUNCTION CHERMIT(CX,CY)
    USE wmcomm,ONLY: rkind
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CX,CY
    COMPLEX(rkind):: CHERMIT

    CHERMIT=0.5D0*(CX-DCONJG(CY))
    RETURN
  END FUNCTION CHERMIT

!     ****** CALCULATE ENERGY FLUX ******

  SUBROUTINE wm_pwrflux

    USE wmcomm
    IMPLICIT NONE

    RETURN
  END SUBROUTINE wm_pwrflux

!     ****** CALCULATE ANTENNA IMPEDANCE ******

  SUBROUTINE wm_pwrant

    USE wmcomm
    USE wmsetm
    IMPLICIT NONE
    COMPLEX(rkind),ALLOCATABLE:: CFVP(:,:,:)
    COMPLEX(rkind):: CCE1,CCE2,CCE3,CPRAD_L
    INTEGER:: NRANT,NR,MDX,NDX
    

    CPRAD=0.D0
    DO NDX=1,ndsiz
       DO MDX=1,mdsiz
          CPRADK(MDX,NDX)=(0.D0,0.D0)
       END DO
    END DO

    IF(MODEEG.EQ.0.AND.MODEWG.EQ.0) THEN
       ALLOCATE(CFVP(nthmax,nhhmax,3))
       NRANT=0
       DO NR=1,NRMAX 
          IF(XR(NR)/RD.LT.1.D0) NRANT=NR
       ENDDO

       DO NR=NRANT-1,NRMAX
          CALL wm_setv(NR,CFVP)
          DO MDX=1,MDSIZ
             DO NDX=1,NDSIZ
                CCE1=CEFLDK(1,MDX,NDX,NR+1)
                CCE2=CEFLDK(2,MDX,NDX,NR+1)
                CCE3=CEFLDK(3,MDX,NDX,NR+1)
                CPRAD_L=DCONJG(CCE1)*CFVP(MDX,NDX,1) &
                       +DCONJG(CCE2)*CFVP(MDX,NDX,2) &
                       +DCONJG(CCE3)*CFVP(MDX,NDX,3)
                CPRADK(MDX,NDX)=CPRADK(MDX,NDX)+CPRAD_L
             ENDDO
          ENDDO
       ENDDO
       DO NDX=1,ndsiz
          DO MDX=1,mdsiz
             CPRAD=CPRAD+CPRADK(MDX,NDX)
          END DO
       END DO
    ENDIF

    RETURN
  END SUBROUTINE wm_pwrant

!     ****** CALCULATE CURRENT DRIVE EFFICIENCY ******

!      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
!      Z  = ZEFF   : EFFECTIVE Z
!      XL = X / RR : NORMALIZED X
!      YL = Y / RR : NORMALIZED Y
!      ID : 0 : LANDAU DAMPING
!           1 : TTMP

  FUNCTION W1CDEF(WT,Z,XL,YL,ID)
    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: WT,Z,XL,YL
    INTEGER,INTENT(IN):: ID
    REAL(rkind):: W1CDEF
    REAL(rkind):: R,D,QC,A,RM,RC,W,EFF0,EFF1,Y2,Y1,EFF2,YT,ARG,EFF3

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

    Y2=(R+XL)/(1.D0+R)
    IF(Y2.LT.0.D0) Y2=0.D0
    Y1=SQRT(Y2)
    EFF2=1.D0+A*(Y1/W)**3

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

    W1CDEF=EFF0*EFF1*EFF2*EFF3

    RETURN
  END FUNCTION W1CDEF

!     ****** CALCULATE ABSORBED POWER ******

  SUBROUTINE wm_pabs
    USE wmcomm
    USE wmsetf
    USE wmprof
    USE wmsub
    IMPLICIT NONE
    INTEGER,ALLOCATABLE:: iposa(:),ilena(:)
    REAL(rkind),ALLOCATABLE:: DS(:),DSS(:,:,:)
    COMPLEX(rkind),ALLOCATABLE:: CPF1(:,:),CPF2(:,:)
    COMPLEX(rkind),ALLOCATABLE:: CPF1_B(:,:)
    COMPLEX(rkind),ALLOCATABLE,DIMENSION(:,:,:,:):: CPABSKM
    COMPLEX(rkind),ALLOCATABLE,DIMENSION(:,:,:,:):: CPABSKC
    COMPLEX(rkind),ALLOCATABLE,DIMENSION(:):: CFA,CFB
    REAL(rkind):: RN(NSMAX),RTPR(NSMAX),RTPP(NSMAX),RU(NSMAX)
    REAL(rkind):: RNC(NSMAX),RTPRC(NSMAX),RTPPC(NSMAX),RUC(NSMAX)
    COMPLEX(rkind):: CDV(3,3,3),CDW(3,3,3)
    COMPLEX(rkind):: CW,CDV11M,CDV11C,CDV12M,CDV12C,CDV13M,CDV13C
    COMPLEX(rkind):: CDV21C,CDV22C,CDV23C,CDV31C,CDV32C,CDV33C
    COMPLEX(rkind):: CEMM11,CEMC11,CEMM12,CEMC12,CEMM13,CEMC13
    COMPLEX(rkind):: CEMC21,CEMP21,CEMC22,CEMC23
    COMPLEX(rkind):: CEMP31,CEMC31,CEMC32,CEMC33
    COMPLEX(rkind):: CCE1,CCE2,CCE3
    COMPLEX(rkind):: CPM11,CPC11,CPM12,CPC12,CPM13,CPC13
    COMPLEX(rkind):: CPC21,CPP21,CPC22,CPP22,CPC23
    COMPLEX(rkind):: CPC31,CPP31,CPC32,CPP32,CPC33
    COMPLEX(rkind):: CPABSM,CPABSC
    REAL(rkind):: XRHOC,XRHOM,XRHOP,XRHOMH,XRHOPH,DRHOM,DRHOP
    REAL(rkind):: QPMH,QPC,FMHM,FMHC,FCMH,FCPH
    REAL(rkind):: DPSIPDRHOMH,DPSIPDRHOC,DRHOPM
    REAL(rkind):: FACT1M,FACT1C,FACT2C,FACT2P,FACT3C,FACT3P
    REAL(rkind):: VTE,WW,RLNLMD,RT,VTEC,RLNLMDC,RTC
    REAL(rkind):: BABS,BSUPTH,BSUPPH,RKPR,W,XL,YL,EFCD
    INTEGER:: NM,NBSTX,NBEDX,NR,NS,NHH,NTH,NDX,KDX,MDX,LDX,I,J,MLX,NKX
    INTEGER:: ND,NK,KDMAX_TMP,LL,KD,KK,KKX,KX1,NX1,MD,NRPP,ML,LDMAX_TMP,LLX
    INTEGER:: LD,LX1,MX1,K,NHHF,NTHF,NN,MM

    ALLOCATE(DS(nrmax),DSS(nthmax,nhhmax,nrmax))
    ALLOCATE(CPF1(nthmax,nhhmax),CPF2(nthmax,nhhmax))
    ALLOCATE(CPABSKM(nthmax,nthmax,nhhmax,nhhmax))
    ALLOCATE(CPABSKC(nthmax,nthmax,nhhmax,nhhmax))

    CW=2.D0*PI*CMPLX(RF,RFI)*1.D6

    NM=nrmax*nsmax*nthmax*nhhmax

    DO NS=1,NSMAX
       DO NR=1,NRMAX
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                CPABS(NTH,NHH,NR,NS)=(0.D0,0.D0)
             ENDDO !nth
          ENDDO ! nhh
       ENDDO !nr
    END DO !ns
    
    DO NR=1,NRMAX
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             PCUR(NTH,NHH,NR)=0.D0
          ENDDO !nth
       ENDDO !nhh
    ENDDO !nr

    NR=nr_start
    DO NS=1,NSMAX
       CALL wm_setf(NR,NS)
       DO NDX=1,NDSIZ
          DO KDX=1,KDSIZ_F
             DO MDX=1,MDSIZ
                DO LDX=1,LDSIZ_F
                   DO J=1,3
                      DO I=1,3
                         CPSF(I,J,LDX,MDX,KDX,NDX,2) &
                              =CPSF(I,J,LDX,MDX,KDX,NDX,3)
                      ENDDO !i
                   ENDDO !j
                ENDDO !ldx
             ENDDO !mdx
          ENDDO !kdx
       ENDDO !ndx

       DO NR=nr_start,MIN(nr_end,nrmax-1)
          DO NKX=1,NDSIZ
             DO KDX=1,KDSIZ
                DO MLX=1,MDSIZ
                   DO LDX=1,LDSIZ
                      CPABSKM(LDX,MLX,KDX,NKX)=(0.D0,0.D0)
                      CPABSKC(LDX,MLX,KDX,NKX)=(0.D0,0.D0)
                   ENDDO !ldx
                ENDDO !mlx
             ENDDO !kdx
          ENDDO !nkx

          CALL wm_setf(NR+1,NS)
          DO NDX=1,NDSIZ
             DO KDX=1,KDSIZ_F
                DO MDX=1,MDSIZ
                   DO LDX=1,LDSIZ_F
                      DO J=1,3
                         DO I=1,3
                            CPSF(I,J,LDX,MDX,KDX,NDX,1) &
                                 =CPSF(I,J,LDX,MDX,KDX,NDX,2)
                            CPSF(I,J,LDX,MDX,KDX,NDX,2) &
                                 =CPSF(I,J,LDX,MDX,KDX,NDX,3)
                         ENDDO !i
                      ENDDO !j
                   ENDDO !ldx
                ENDDO !mdx
             ENDDO !kdx
          ENDDO !ndx
         


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

          DRHOM =XRHO(NR+1)-XRHO(NR)
          IF(NR.EQ.NRMAX) THEN
             NRPP=NR+1
             DRHOP=XRHO(NR+1)-XRHO(NR)
          ELSE
             NRPP=NR+2
             DRHOP=XRHO(NR+2)-XRHO(NR+1)
          ENDIF

          IF(MODELG.EQ.3) THEN
             QPMH=0.5D0*(QPS(NR)+QPS(NR+1))
             QPC=QPS(NR+1)
             DPSIPDRHOMH=2.D0*PSITA*XRHOMH/QPMH
             DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
          ELSE
             DPSIPDRHOMH=2.D0*PSIPA*XRHOMH
             DPSIPDRHOC =2.D0*PSIPA*XRHOC
          ENDIF

          DRHOPM= 0.5D0*(DRHOM+DRHOP)

          FMHM=0.5D0
          FMHC=0.5D0
          FCMH=DRHOM/(2.D0*DRHOPM)
          FCPH=DRHOP/(2.D0*DRHOPM)

!        ND : (n - n0) / Np
!        KD : n'              = (k - n) / Np
!        NN : n               = n0 +  ND       * Np
!        NK : k = n + n' * Np = n0 + (ND + KD) * Np

!        MD : m - m0
!        LD : m'              = l - m
!        MM : m               = m0 +  MD
!        ML : l = m + m'      = m0 + (MD + LD)

          DO ND=NDMIN,NDMAX
             NDX=ND-NDMIN+1
             DO NK=NDMIN,NDMAX
                NKX=NK-NDMIN+1
                KDMAX_TMP=KDMAX-1
                IF(KDMAX==0) KDMAX_TMP=0
                DO KK=KDMIN,KDMAX_TMP
                   KKX=KK-KDMIN+1
                   KD=NK-ND
                   IF (KD+KK .GE. KDMIN_F .AND. &
                      (KD+KK .LE. KDMAX_F .OR. KDMAX_F==0))THEN
                      KX1=KD+KK-KDMIN_F+1
                      NX1=MOD( ND   -NDMIN+4*NDSIZ,NDSIZ)+1

                      DO MD=MDMIN,MDMAX
                         MDX=MD-MDMIN+1
                         DO ML=MDMIN,MDMAX
                            MLX=ML-MDMIN+1
                            LDMAX_TMP=LDMAX-1
                            IF(LDMAX==0) LDMAX_TMP=0
                            DO LL=LDMIN,LDMAX_TMP
                               LLX=LL-LDMIN+1
                               LD= ML-MD
                               IF (LD+LL .GE. LDMIN_F .AND. &
                                  (LD+LL .LE. LDMAX_F .OR. LDMAX_F==0))THEN
                                  LX1=LD+LL-LDMIN_F+1

                                  MX1=MOD( MD   -MDMIN+4*MDSIZ,MDSIZ)+1
                                  DO K=1,2
                                     DO J=1,3
                                        DO I=1,3
                                           CDV(I,J,K) &
                                                =CPSF(I,J,LX1,MX1,KX1,NX1,K)
                                           CDW(I,J,K) &
                                                =CPSF(I,J,LX1,MX1,KX1,NX1,K)
                                        ENDDO
                                     ENDDO
                                  ENDDO

                                  FACT1M=XRHOMH/XRHOM
                                  FACT1C=XRHOMH/XRHOC
                                  FACT2C=XRHOMH/XRHOC
                                  FACT2P=XRHOPH/XRHOC
                                  FACT3C=XRHOMH/XRHOC
                                  FACT3P=XRHOPH/XRHOC

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

!     --- R COMPONENT OF MAXWELL EQUATION ---

                                  CEMM11=0.5D0*CDV11M*FACT1M ! /XRHOMH/XRHOMH 
                                  CEMC11=0.5D0*CDV11C*FACT1C ! /XRHOMH/XRHOMH
                                  CEMM12= FMHM*CDV12M        !       /XRHOMH
                                  CEMC12= FMHC*CDV12C        !       /XRHOMH
                                  CEMM13= FMHM*CDV13M        !       /XRHOMH
                                  CEMC13= FMHC*CDV13C        !       /XRHOMH

!     --- THETA COMPONENT OF MAXWELL EQUATION ---

                                  CEMC21= FCMH*CDV21C*FACT2C  !      /XRHOMH
                                  CEMP21= FCPH*CDV21C*FACT2P  !      /XRHOPH
                                  CEMC22=      CDV22C       
                                  CEMC23=      CDV23C        

!     --- PHI COMPONENT OF MAXWELL EQUATION ---

                                  CEMC31= FCMH*CDV31C*FACT3C !       /XRHOMH
                                  CEMP31= FCPH*CDV31C*FACT3P !       /XRHOPH
                                  CEMC32=      CDV32C       
                                  CEMC33=      CDV33C

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

                                  CPABSM=CI*EPS0*(VC*VC/CW) &
                                       *(CPM11 &
                                        +CPM12 &
                                        +CPM13 &
                                        )*DPSIPDRHOMH*DRHOM

                                  CPABSC=CI*EPS0*(VC*VC/CW) &
                                       *(CPC11 &
                                        +CPC12 &
                                        +CPC13 &
                                        )*DPSIPDRHOMH*DRHOM &
                                       +CI*EPS0*(VC*VC/CW) &
                                       *(CPC22*2d0 &
                                        +CPC33*2d0 &
                                        +CPC23*2d0 &
                                        +CPC32*2d0 &
                                        +CPC21 &
                                        +CPC31 &
                                        +CPP21 &
                                        +CPP31 &
                                        )*DPSIPDRHOC*DRHOPM

                                  CPABSKM(LLX,MDX,KKX,NDX) &
                                       =CPABSKM(LLX,MDX,KKX,NDX) &
                                       +0.5D0*CPABSM

                                  CPABSKC(LLX,MDX,KKX,NDX) &
                                       =CPABSKC(LLX,MDX,KKX,NDX) &
                                       +0.5D0*CPABSC
                               ENDIF

                            ENDDO !ll
                         ENDDO !ml
                      ENDDO !md
                   ENDIF
                ENDDO !kk
             ENDDO !nk
          ENDDO !nd

!     +++++ CALCULATE PABS IN MODE NUMBER SPACE +++++


          DO NDX=1,NDSIZ
             KK=0
             KKX=KK-KDMIN+1
             DO MDX=1,MDSIZ
                LL=0
                LLX=LL-LDMIN+1
                CPABSK(MDX,NDX,NR,NS)  = CPABSK(MDX,NDX,NR,NS) &
                                       + CPABSKM(LLX,MDX,KKX,NDX)
                CPABSK(MDX,NDX,NR+1,NS)= CPABSK(MDX,NDX,NR+1,NS) &
                                       + CPABSKC(LLX,MDX,KKX,NDX)
             ENDDO !mdx
          ENDDO !ndx


          !     +++++ CALCULATE PABS IN REAL SPACE +++++

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
                   ENDDO !ll
                ENDDO !kk
                CALL WMSUBE(CPF1,CPF2)
                DO NHH=1,NHHMAX
                   DO NTH=1,NTHMAX
                      CPABS(NTH,NHH,NR,NS)=CPABS(NTH,NHH,NR,NS) &
                           +CPF2(NTH,NHH)
                      CPABS3D(NTH,NHH,NR,NS)=CPABS3D(NTH,NHH,NR,NS) &
                           +CPF2(NTH,NHH)
                   ENDDO !nth
                ENDDO !nhh
             ENDDO !mdx
          ENDDO !ndx


!     +++++ CALCULATE DRIVEN CURRENT IN REAL SPACE +++++

          IF (NS .EQ. 1) THEN
             CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
             VTE=SQRT(RTPR(1)*AEE*1.D3/(PA(1)*AMP))
             WW=RF*2.D0*PI
             IF(RN(1).LE.0.D0) THEN
                RLNLMD=15.D0
             ELSE
                RT=(RTPR(1)+2*RTPP(1))/3.D0
                RLNLMD=16.1D0 - 1.15D0*LOG10(RN(1)) &
                              + 2.30D0*LOG10(RT)
             ENDIF
             CALL WMCDEN(NR+1,RNC,RTPRC,RTPPC,RUC)
             VTEC=SQRT(RTPRC(1)*AEE*1.D3/(PA(1)*AMP))
             WW=RF*2.D0*PI
             IF(RNC(1).LE.0.D0) THEN
                RLNLMDC=15.D0
             ELSE
                RTC=(RTPRC(1)+2*RTPPC(1))/3.D0
                RLNLMDC=16.1D0 - 1.15D0*LOG10(RNC(1)) &
                       + 2.30D0*LOG10(RTC)
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
                         NHHF=(NHH-1)*FACTOR_NHH+1
                         NTHF=(NTH-1)*FACTOR_NTH+1
                         CALL WMCMAG(NR,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
                         RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
                         IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
                         IF(ABS(WW/RKPR).LT.VC) THEN
                            W=WW/(RKPR*VTE)
                            XL=(RPST(NTHF,NHHF,NR)-RR  )/RR
                            YL=(ZPST(NTHF,NHHF,NR)-0.D0)/RR
                            EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                            IF(W.LT.0.D0) EFCD=-EFCD
                            IF (RN(1).GT.0.D0) THEN
                               PCUR(NTH,NHH,NR)=PCUR(NTH,NHH,NR) &
                                    +0.384D0*RTPR(1)*EFCD &
                                    /(RN(1)*RLNLMD) &
                                    *DBLE(CPF2(NTH,NHH)) &
                                    /(2.D0*PI*RPST(NTHF,NHHF,NR))
                            END IF
                         ENDIF
                      ENDDO !nth
                   ENDDO !nhh
                   DO KKX=1,KDSIZ
                      DO LLX=1,LDSIZ
                         CPF1(LLX,KKX)=CPABSKC(LLX,MDX,KKX,NDX)
                      ENDDO !llx
                   ENDDO !kkx
                   CALL WMSUBE(CPF1,CPF2)
                   DO NHH=1,NHHMAX
                      DO NTH=1,NTHMAX
                         NHHF=(NHH-1)*FACTOR_NHH +1
                         NTHF=(NTH-1)*FACTOR_NTH +1
                         CALL WMCMAG(NR+1,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
                         RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
                         IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
                         IF(ABS(WW/RKPR).LT.VC) THEN
                            W=WW/(RKPR*VTEC)
                            XL=(RPST(NTHF,NHHF,NR+1)-RR  )/RR
                            YL=(ZPST(NTHF,NHHF,NR+1)-0.D0)/RR
                            EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                            IF(W.LT.0.D0) EFCD=-EFCD
                            IF (RNC(1).GT.0.D0) THEN
                               PCUR(NTH,NHH,NR+1)=PCUR(NTH,NHH,NR+1) &
                                    +0.384D0*RTPRC(1)*EFCD &
                                    /(RNC(1)*RLNLMDC) &
                                    *DBLE(CPF2(NTH,NHH)) &
                                    /(2.D0*PI*RPST(NTHF,NHHF,NR+1))
                            END IF
                         ENDIF
                      ENDDO !nth
                   ENDDO !nhh
                ENDDO !md
             ENDDO !nd
          ENDIF !ns=1
       ENDDO !nr
    ENDDO !ns

    DO NS=1,NSMAX
       DO NR=1,NRMAX
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                PABS(NTH,NHH,NR,NS)=DBLE(CPABS(NTH,NHH,NR,NS))
             END DO
          END DO
          DO NDX=1,NDSIZ
             DO MDX=1,MDSIZ
                PABSK(MDX,NDX,NR,NS)=DBLE(CPABSK(MDX,NDX,NR,NS))
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE(CPABSKM)
    DEALLOCATE(CPABSKC)

    RETURN
  END SUBROUTINE wm_pabs
END MODULE wmemfp
