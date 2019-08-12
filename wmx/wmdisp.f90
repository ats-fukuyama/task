! wmdisp.f90

MODULE wmdisp

  PRIVATE
  PUBLIC wmtnsr

CONTAINS

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE WMTNSR(NR,NS)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS
    INTEGER:: NHH,NTH,ND,MD
    INTEGER:: MODELP_save,MODELV_save
    
    DO NHH=1,NHHMAX_F
    DO NTH=1,NTHMAX_F
       DO ND=NDMIN_F,NDMAX_F
       DO MD=MDMIN_F,MDMAX_F
          CTNSR(1,1,MD,ND,NTH,NHH)=0.D0
          CTNSR(1,2,MD,ND,NTH,NHH)=0.D0
          CTNSR(1,3,MD,ND,NTH,NHH)=0.D0
          CTNSR(2,1,MD,ND,NTH,NHH)=0.D0
          CTNSR(2,2,MD,ND,NTH,NHH)=0.D0
          CTNSR(2,3,MD,ND,NTH,NHH)=0.D0
          CTNSR(3,1,MD,ND,NTH,NHH)=0.D0
          CTNSR(3,2,MD,ND,NTH,NHH)=0.D0
          CTNSR(3,3,MD,ND,NTH,NHH)=0.D0
       ENDDO
       ENDDO
    ENDDO
    ENDDO


    IF((MOD(MODELA,2).EQ.1).AND.(NS.EQ.3)) THEN
       CALL WMTNAX(NR)
    ELSEIF((MOD(MODELA/2,2).EQ.1).AND.(NS.EQ.1)) THEN
       CALL WMTNEX(NR)
    ELSEIF(NS.EQ.5.OR.NS.EQ.6) THEN
       CALL WMTNDK(NR,NS)
    ELSE
       IF(MODELP(NS).LT.0) THEN
          CALL WMTNSX(NR,NS)
       ELSE
          IF(MODELV(NS).EQ.5) THEN
             CALL WMDPDK(NR,NS)
          ELSE
             CALL WMDPIN(NR,NS)
          ENDIF
       ENDIF
    ENDIF

!      IF(NR.EQ.1) THEN
!      WRITE(6,*) 'WMDISP: NR,NS=',NR,NS
!      WRITE(6,'(1P6E12.4)') 
!     &     CTNSR(1,1,0,0,1,1),CTNSR(1,2,0,0,1,1),CTNSR(1,3,0,0,1,1),
!     &     CTNSR(2,1,0,0,1,1),CTNSR(2,2,0,0,1,1),CTNSR(2,3,0,0,1,1),
!     &     CTNSR(3,1,0,0,1,1),CTNSR(3,2,0,0,1,1),CTNSR(3,3,0,0,1,1)
!      ENDIF
!      IF(NR.EQ.2) STOP

    RETURN
  END SUBROUTINE WMTNSR

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE WMTNSX(NR,NS)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    USE libdsp,ONLY: DSPFN
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS
    REAL(rkind):: RN(NSMAX),RTPR(NSMAX),RTPP(NSMAX),RU(NSMAX)
    INTEGER:: NHH,NTH,ND,NN,MD,MM,NSS,NC
    REAL(rkind):: WW,BABS,BSUPTH,BSUPH,RKTH,RKPH,RKPR
    REAL(rkind):: RNPR,DTT,DTX,AM,AE,WP,WC,RNPP2,RKPP2,RKX2,RKPP,RKT2,RKR2
    REAL(rkind):: UXX2,UYY2,T,RWC,RKTPR,RKTPP
    COMPLEX(rkind):: CW,CWN,CWU,CWP,CWC,CPERP,CPARA,CCROS,CPERM
    COMPLEX(rkind):: CGZ0,CGZ,CADD,CPERP1,CPERP2,CPERM2,CPARA1,CPARA2
    COMPLEX(rkind):: CCROS1,CCROS2



      CW=2.D0*PI*CRF*1.D6
      WW=DBLE(CW)

      CALL WMCDEN(NR,RN,RTPR,RTPP,RU)

!      IF(NS.EQ.1.AND.NR.EQ.1) THEN
!         WRITE(6,'(A,1P6E12.4)') 'RN  :',(RN(NS1),NS1=1,NSMAX)
!         WRITE(6,'(A,1P6E12.4)') 'RTPR:',(RTPR(NS1),NS1=1,NSMAX)
!         WRITE(6,'(A,1P6E12.4)') 'RTPP:',(RTPP(NS1),NS1=1,NSMAX)
!      ENDIF

      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
!         IF(NR.EQ.1.AND.NTH.EQ.1.AND.NHH.EQ.1) THEN
!            WRITE(6,*) 'BABS:',BABS,BSUPTH,BSUPPH
!         ENDIF

         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD

            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH

!            WRITE(11,'(5I5,1P2E12.4)') NR,NN,MM,NHH,NTH,RKPR,BABS

            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            RNPR=VC*RKPR/WW

            IF(MODELP(NS).EQ.-4) THEN
               DTT=1.D0
               DTX=0.D0
               DO NSS=1,NSMAX
                  AM=PA(NSS)*AMP
                  AE=PZ(NSS)*AEE
                  WP=AE*AE*RN(NSS)*1.D20/(AM*EPS0)
                  WC=AE*BABS/AM
                  DTT=DTT-WP/(WW*WW-WC*WC)
                  DTX=DTX-WP*WC/((WW*WW-WC*WC)*WW)
               ENDDO
               RNPP2=((DTT-RNPR**2)**2-DTX**2)/(DTT-RNPR**2)
               RKPP2=RNPP2*WW*WW/(VC*VC)
               RKX2=     MM*MM*RG22(NTH,NHH,NR)*XRHO(NR)**2 &
                   +2.D0*MM*NN*RG23(NTH,NHH,NR)*XRHO(NR) &
                   +     NN*NN*RG33(NTH,NHH,NR)
               IF(RKPP2.GT.0.D0) THEN
!                  RKPP=SQRT(RKPP2)
                  RKT2=RKX2-RKPR**2
                  IF(RKT2.GT.0.D0) THEN
                     IF(RKT2.LE.RKPP2) THEN
                        RKR2=RKPP2-RKT2
                     ELSE
                        RKR2=0.D0
                     ENDIF
                  ELSE
                     RKT2=0.D0
                     RKR2=RKPP2
                  ENDIF
                  UXX2= RKX2/RKPP2
                  UYY2= RKR2/RKPP2
               ELSE
                  RKPP2=0.D0
                  RKPP=0.D0
                  RKT2=0.D0
                  RKR2=0.D0
                  UXX2=0.D0
                  UYY2=0.D0
               ENDIF
            ELSE
               RKPP2=0.D0
               RKPP=0.D0
               RKT2=0.D0
               RKR2=0.D0
               UXX2=0.D0
               UYY2=0.D0
            ENDIF

            AM=PA(NS)*AMP
            AE=PZ(NS)*AEE
            IF(MODELP(NS).EQ.-1) THEN
               CWN=CW-RU(NS)* RKPR+CI*PZCL(NS)*WW
               CWU=CW-RU(NS)* RKPR
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0)
               CWC=AE*BABS/AM
               IF(NS.EQ.1) THEN
                  CPERP=(0.D0,0.D0)
                  CPARA=CI*CWP/(CW*CW*PZCL(NS))
               ELSE
                  CPERP= CWP*CWU*CWN/(CW*CW*CWC**2)
                  CPARA=CI*CWP/(CW*CW*PZCL(NS))
               END IF
               CCROS=(0.D0,0.D0)
               CPERM=(0.D0,0.D0)
            ELSE IF(MODELP(NS).EQ.-2) THEN
               CWN=CW-RU(NS)* RKPR+CI*PZCL(NS)*WW
               CWU=CW-RU(NS)* RKPR
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0)
               CWC=AE*BABS/AM
               CPERP=-   CWP*CWU*CWN/(CW*CW*(CWN**2-CWC**2))
               CCROS= CI*CWP*CWC*CWU/(CW*CW*(CWN**2-CWC**2))
               CPARA=-   CWP/(CW*CW) &
                         *(CW*CW/(CWU*CWN) &
                          +2*RKPP2*RU(NS)*RU(NS)/(CWN*CWN-CWC*CWC) &
                           *CWN/CWU)
               CPERM=(0.D0,0.D0)
!               WRITE(6,*) NR,NS,CPERP,CPARA
            ELSE IF(MODELP(NS).EQ.-3) THEN
               RT=(RTPR(NS)+2*RTPP(NS))/3.D0
               RKTPR=ABS(RKPR)*SQRT(2.D0*RT*AEE*1.D3/AM)
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0*CW*CW)
               CWC=AE*BABS/AM
               CGZ0=(CW       -RKPR*RU(NS))/RKTPR
               CPERP=0.D0
               CCROS=0.D0
               CPARA=0.D0
               DO NC=-1,1
                  CGZ=(CW-NC*CWC-RKPR*RU(NS))/RKTPR
                  CALL DSPFN(CGZ,CZ,CDZ)
                  IF(NC.EQ.-1) THEN
                     CPERP=CPERP+   CWP*CGZ0*CZ/2
                     CCROS=CCROS-CI*CWP*CGZ0*CZ/2
                  ELSEIF(NC.EQ.0) THEN
!                     CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS))
                     cadd=1+RKPR*RU(NS)*CW/(CW-RKPR*RU(NS))**2 !
!                     CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD*CADD
                     CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP=CPERP+   CWP*CGZ0*CZ/2
                     CCROS=CCROS+CI*CWP*CGZ0*CZ/2
                  ENDIF
               ENDDO
               CPERM=(0.D0,0.D0)
            ELSE IF(MODELP(NS).EQ.-4) THEN
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0*CW*CW)
               RWC=AE*BABS/AM
               RKTPR=ABS(RKPR)*SQRT(2.D0*RTPR(NS)*AEE*1.D3/AM)
               RKTPP=RTPP(NS)*AEE*1.D3/(AM*RWC*RWC)
               CWC=DCMPLX(RWC,0.D0)
               CGZ0=(CW       -RKPR*RU(NS))/RKTPR
               CPERP1= 0.D0
               CPERP2= 0.D0
               CPERM2= 0.D0
               CPARA1= 0.D0
               CPARA2= 0.D0
               CCROS1= 0.D0
               CCROS2= 0.D0
               DO NC=-2,2
                  CGZ =(CW-NC*CWC-RKPR*RU(NS))/RKTPR
                  CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS)-NC*CWC)
                  CALL DSPFN(CGZ,CZ,CDZ)
                  IF(NC.EQ.-2) THEN 
                     CPERP2=CPERP2+CWP*RKTPP*CGZ0*CZ
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*CZ
                  ELSEIF(NC.EQ.-1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
                     CPARA2=CPARA2 &
                           -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
                     CCROS1=CCROS1-CI*CWP*CGZ0*CZ/2
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*2*CZ
                  ELSEIF(NC.EQ.0) THEN
                     CPERM2=CPERM2+CWP*RKTPP*CGZ0*4*CZ
                     CPARA1=CPARA1-CWP*CDZ*CGZ*CGZ0*CADD*CADD
                     CPARA2=CPARA2 &
                           +CWP*RKTPP*CGZ0*2*CGZ*CDZ*CADD*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
                     CPARA2=CPARA2 &
                           -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
                     CCROS1=CCROS1+CI*CWP*CGZ0*CZ/2
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*2*CZ
                  ELSEIF(NC.EQ.2) THEN
                     CPERP2=CPERM2+CWP*RKTPP*CGZ0*CZ
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*CZ
                  ENDIF
               ENDDO

!               CPERP1= CWP*CGZ(0)*(CZ(1)+CZ(-1))/2
!               CPERP2=-CWP*RKTPP*CGZ(0)
!     &                    *(CZ(1)+CZ(-1)-CZ(2)-CZ(-2))    
!               CPERM2=-CWP*RKTPP*CGZ(0)
!     &                    *(2*CZ(1)+2*CZ(-1)-4*CZ(0))
!               CPARA1=-CWP*CDZ(0)*CGZ(0)*CGZ(0)*CADD(0)*CADD(0)
!               CPARA2= CWP*RKTPP*CGZ(0)
!     &                 *(2*CGZ( 0)*CDZ( 0)*CADD( 0)*CADD( 0)
!     &                    -CGZ( 1)*CDZ( 1)*CADD( 1)*CADD( 1)
!     &                    -CGZ(-1)*CDZ(-1)*CADD(-1)*CADD(-1))
!               CCROS1= CI*CWP*CGZ(0)*(CZ(1)-CZ(-1))/2
!               CCROS2=-CI*CWP*RKTPP*CGZ(0)
!     &                   *(2*CZ(1)-2*CZ(-1)-CZ(2)+CZ(-2))

               CPERP=CPERP1+0.5D0*CPERP2*RKPP2
               CPERM=       0.5D0*CPERM2*RKPP2
               CPARA=CPARA1+0.5D0*CPARA2*RKPP2
               CCROS=CCROS1+0.5D0*CCROS2*RKPP2
            ELSE
               CPERP=(0.D0,0.D0)
               CPARA=(0.D0,0.D0)
               CCROS=(0.D0,0.D0)
               CPERM=(0.D0,0.D0)
            ENDIF

!      IF(NR.EQ.1) THEN
!         WRITE(6,*) 'CPERP,CPERM,CPARA:',CPERP,CPERM,CPARA
!         WRITE(6,*) 'UXX2,UYY2:',UXX2,UYY2
!      ENDIF

!            CTNSR(1,1,MD,ND,NTH,NHH)= CPERP + UXX2*CPERM
!            CTNSR(1,2,MD,ND,NTH,NHH)=-CCROS
!            CTNSR(1,3,MD,ND,NTH,NHH)= 0.D0
!            CTNSR(2,1,MD,ND,NTH,NHH)= CCROS
!            CTNSR(2,2,MD,ND,NTH,NHH)= CPERP + UYY2*CPERM
!            CTNSR(2,3,MD,ND,NTH,NHH)= 0.D0
!            CTNSR(3,1,MD,ND,NTH,NHH)= 0.D0
!            CTNSR(3,2,MD,ND,NTH,NHH)= 0.D0
!            CTNSR(3,3,MD,ND,NTH,NHH)= CPARA

            CTNSR(1,1,MD,ND,NTH,NHH) &
           =CTNSR(1,1,MD,ND,NTH,NHH) + CPERP + UXX2*CPERM
            CTNSR(1,2,MD,ND,NTH,NHH) &
           =CTNSR(1,2,MD,ND,NTH,NHH) - CCROS
            CTNSR(2,1,MD,ND,NTH,NHH) &
           =CTNSR(2,1,MD,ND,NTH,NHH) + CCROS
            CTNSR(2,2,MD,ND,NTH,NHH) &
           =CTNSR(2,2,MD,ND,NTH,NHH) + CPERP + UYY2*CPERM
            CTNSR(3,3,MD,ND,NTH,NHH) &
           =CTNSR(3,3,MD,ND,NTH,NHH) + CPARA

!            IF(NR.EQ.30.AND.NTH.EQ.1.AND.MDX.EQ.1) THEN
!               WRITE(6,*) 'RN,RTPR,RTPP=',RN(1)/1.D20,
!     &                    RTPR(1)/(AEE*1.D3),RTPP(1)/(AEE*1.D3)
!               WRITE(6,*) 'BABS,QR=',BABS,QR
!               WRITE(6,*) 'RKPR,RKPP2=',RKPR,RKPP2
!               WRITE(6,*) 'CGZ(0)=',CGZ(0)
!               WRITE(6,*) 'CDZ(0)=',CDZ(0)
!               WRITE(6,*) 'CTNSR(3,3)=',CTNSR(3,3,NTH,MDX)
!            ENDIF
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      RETURN
  END SUBROUTINE WMTNSX

!     ****** IMPORT FROM TASK/DP ******

  SUBROUTINE WMDPIN(NR,NS)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS
    COMPLEX(rkind):: CDTNS(3,3)
    INTEGER NDSIZX_TMP, MDSIZX_TMP

      CW=2.D0*PI*CRF*1.D6
      WW=DBLE(CW)

      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON)
      NDSIZX_TMP=NDSIZX_F
      MDSIZX_TMP=MDSIZX_F
      IF (NDSIZX_TMP==1)NDSIZX_TMP=0
      IF (MDSIZX_TMP==1)MDSIZX_TMP=0

      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
         CALL WMCMAG_F(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         BTH=ABS(BSUPTH*RA*RHON)
!         BTH=ABS(BSUPTH*RA)
         WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
!         RKPR_EFF=SQRT(CW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
         RKPR_EFF=SQRT(CW*BTHOBN(NR)/SQRT(8.D0*WTPR)/RR)
         IF (RKPR_EFF < 1d-5)then
           print *,RKPR_EFF,BTH,WKPT
           RKPR_EFF=1.D-5
           STOP
         endif
         DO ND=-NDSIZX_TMP,NDSIZX_TMP
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX_TMP,MDSIZX_TMP
            MM=NTH0+MD

            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH

            IF(ABS(RKPR).LE.RKPR_EFF) RKPR=SIGN(RKPR_EFF,RKPR)
            RNPR=VC*RKPR/WW

            IF(MODELP(NS).EQ.5.OR.MODELP(NS).EQ.15) THEN
               DTT=1.D0
               DTX=0.D0
               DO NSS=1,NSMAX
                  AM=PA(NSS)*AMP
                  AE=PZ(NSS)*AEE
                  WP=AE*AE*RN(NSS)*1.D20/(AM*EPS0)
                  WC=AE*BABS/AM
                  DTT=DTT-WP/(WW*WW-WC*WC)
                  DTX=DTX-WP*WC/((WW*WW-WC*WC)*WW)
               ENDDO
               RNPP2=((DTT-RNPR**2)**2-DTX**2)/(DTT-RNPR**2)
               RKPP2=RNPP2*WW*WW/(VC*VC)
               IF (NR .eq. 1)then
                  RKX2=     MM*MM*RGI22(NTH,NHH,NR)*(1d6/XRHO(2))**2 &
                      +2.D0*MM*NN*RGI23(NTH,NHH,NR)*(1d6/XRHO(2)) &
                      +     NN*NN*RGI33(NTH,NHH,NR)
               ELSE
                  RKX2=     MM*MM*RGI22(NTH,NHH,NR)/XRHO(NR)**2 &
                      +2.D0*MM*NN*RGI23(NTH,NHH,NR)/XRHO(NR) &
                      +     NN*NN*RGI33(NTH,NHH,NR)
               ENDIF

               IF(RKPP2.GT.0.D0) THEN
                  RKPP=SQRT(RKPP2)
                  RKT2=RKX2-RKPR**2
                  IF(RKT2.GE.0.D0) THEN
                     IF(RKT2.LE.RKPP2) THEN
                        RKR2=RKPP2-RKT2
                     ELSE
                        RKR2=0.D0
                     ENDIF
                  ELSE
                     RKT2=0.D0
                     RKR2=RKPP2
                  ENDIF
                  UXX2= RKT2/RKPP2
                  UYY2= RKR2/RKPP2
               ELSE
                  RKPP=0.D0
                  RKT2=0.D0
                  RKR2=0.D0
                  UXX2=0.D0
                  UYY2=0.D0
               ENDIF
            ELSE
               RKPP2=0.D0
               RKPP=0.D0
               RKT2=0.D0
               RKR2=0.D0
               UXX2=0.D0
               UYY2=0.D0
            ENDIF

            CKPR=RKPR
            CKPP=RKPP

            CALL DPCALC_2(CW,CKPR,CKPP,RHON,BABS,BTH,NS,CDTNS)

            CPERM=CDTNS(2,2)-CDTNS(1,1)

            CTNSR(1,1,MD,ND,NTH,NHH) &
           =CTNSR(1,1,MD,ND,NTH,NHH) + CDTNS(1,1) + UXX2*CPERM
            CTNSR(1,2,MD,ND,NTH,NHH) &
           =CTNSR(1,2,MD,ND,NTH,NHH) + CDTNS(1,2)
            CTNSR(2,1,MD,ND,NTH,NHH) &
           =CTNSR(2,1,MD,ND,NTH,NHH) + CDTNS(2,1)
            CTNSR(2,2,MD,ND,NTH,NHH) &
           =CTNSR(2,2,MD,ND,NTH,NHH) + CDTNS(1,1) + UYY2*CPERM
            CTNSR(3,3,MD,ND,NTH,NHH) &
           =CTNSR(3,3,MD,ND,NTH,NHH) + CDTNS(3,3)

         ENDDO
         ENDDO
      ENDDO
      ENDDO

      RETURN
  END SUBROUTINE WMDPIN

!     ****** IMPORT FROM TASK/DP ******

  SUBROUTINE WMDPDK(NR,NS)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS
      DIMENSION CDTNS(3,3)

      CW=2.D0*PI*CRF*1.D6

      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON)
      IF(RN(NS).EQ.0.D0) RETURN

      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)

         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD

            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CKPR=RKPR

            CALL DPDKDT(CW,CKPR,NS,NR,NTH,NTH,MM,CDTNS)

            DO j=1,3
            DO i=1,3
               CTNSR(i,j,MD,ND,NTH,NHH) &
                    =CTNSR(i,j,MD,ND,NTH,NHH)+CDTNS(i,j)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      RETURN
  END SUBROUTINE WMDPDK

!     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******

  SUBROUTINE WMTNAX(NR)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR

      CW=2.D0*PI*CRF*1.D6

      CALL WMCPOS(NR,XL)
      IF(XL.LT.RA) THEN
         RNA=PNA*EXP(-(XL/PNAL)**2)*1.D20
         DRN=-2.D0*XL/(PNAL**2)
      ELSE
         RETURN
      ENDIF
      RTA=PTA*AEE*1.D3
      AM=PA(3)*AMP
      AE=PZ(3)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)

      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
!         ANGPH=(NHH-1)*2.D0*PI/NHHMAX
         WC=AE*BABS/AM
         RHOA=VTA/WC
         RHOR=RHOA/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTA*DRN/(BABS*CW*AE*XL)
         ENDIF

         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CX=CW/(ABS(RKPR)*VTA)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

            CTNSR(1,1,MD,ND,NTH,NHH) &
           =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH) &
           =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH) &
           =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH) &
           =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH) &
           =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH) &
           =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH) &
           =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH) &
           =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NHH) &
           =CTNSR(3,3,MD,ND,NTH,NHH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE WMTNAX

!     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******

  SUBROUTINE WMTNEX(NR)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)

      CW=2.D0*PI*CRF*1.D6

      CALL WMCPOS(NR,XL)
      CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
      RNX=RN(1)
      RTX=RTPP(1)
      IF(NR.LT.NRMAX) THEN
         CALL WMCPOS(NR+1,XLP)
         CALL WMCDEN(NR+1,RN,RTPR,RTPP,RU)
         RNXP=RN(1)
      ELSE
         XLP=XL
         RNXP=RNX
      ENDIF
      IF(NR.GT.1) THEN
         CALL WMCPOS(NR-1,XLM)
         CALL WMCDEN(NR-1,RN,RTPR,RTPP,RU)
         RNXM=RN(1)
      ELSE
         XLM=XL
         RNXM=RNX
      ENDIF
      IF(XL.LT.RA) THEN
         RNE= RNX
         DRN=(RNXP-RNXM)/(XLP-XLM)
      ELSE
         RETURN
      ENDIF
      RTE=RTX*AEE*1.D3
      AM=AME
      AE=-AEE
      VTE=SQRT(2.D0*RTE/AM)
      WP2=AE*AE*RNE*1.D20/(AM*EPS0)

      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
!         ANGPH=(NHH-1)*2.D0*PI/NHHMAX
         WC=AE*BABS/AM
         RHOE=VTE/WC
         RHOR=RHOE/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTE*DRN/(BABS*CW*AE*XL)
         ENDIF

         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CX=CW/(ABS(RKPR)*VTE)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

            CTNSR(1,1,MD,ND,NTH,NHH) &
           =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH) &
           =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH) &
           =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH) &
           =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH) &
           =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH) &
           =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH) &
           =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH) &
           =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NHH) &
           =CTNSR(3,3,MD,ND,NTH,NHH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE WMTNEX

!     ****** CALCULATE MAGNETIC DRIFT DIELECTRIC TENSOR ******

  SUBROUTINE WMTNDK(NR,NS)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS

      CW=2.D0*PI*CRF*1.D6

      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON+0.01D0)
      RNAP=RN(NS)*1.D20
      RTAP=(RTPR(NS)+RTPP(NS))/3.D0*1.D3*AEE
      CALL PL_PROF_OLD(RHON)
      RNA=RN(NS)*1.D20
      RTA=(RTPR(NS)+RTPP(NS))/3.D0*1.D3*AEE
      DRN=(RNAP-RNA)/(0.01D0*RA*RNA)
!      WRITE(6,'(I5,1P3E12.4)') NR,RNA,RTA,DRN

      XL=RHON*RA
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)

      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
!         ANGPH=(NHH-1)*2.D0*PI/NHHMAX
         WC=AE*BABS/AM
         RHOA=VTA/WC
         RHOR=RHOA/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTA*DRN/(BABS*CW*AE*XL)
         ENDIF

         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CX=CW/(ABS(RKPR)*VTA)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

            CTNSR(1,1,MD,ND,NTH,NHH) &
           =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH) &
           =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH) &
           =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH) &
           =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH) &
           =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH) &
           =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH) &
           =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH) &
           =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NHH) &
           =CTNSR(3,3,MD,ND,NTH,NHH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE WMTNDK

END MODULE wmdisp
