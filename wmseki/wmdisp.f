C     $Id: wmdisp.f,v 1.26 2014/10/11 05:33:13 fukuyama Exp $
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE WMTNSR(NR,NS)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcomm.inc'
C
C      DO NHH=1,NHHMAX
C      DO NTH=1,NTHMAX
C         DO ND=-NDSIZX,NDSIZX
C         DO MD=-MDSIZX,MDSIZX
C
!!!seki
      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
         DO ND=-NDSIZX_F,NDSIZX_F
         DO MD=-MDSIZX_F,MDSIZX_F
!!!seki
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
C

!      print *,"seki"
      IF((MOD(MODELA,2).EQ.1).AND.(NS.EQ.3)) THEN
         CALL WMTNAX(NR)
      ELSEIF((MOD(MODELA/2,2).EQ.1).AND.(NS.EQ.1)) THEN
         CALL WMTNEX(NR)
      ELSEIF(NS.EQ.5.OR.NS.EQ.6) THEN
         CALL WMTNDK(NR,NS)
      ELSE
         IF(MODELP(NS).LT.0) THEN
            CALL WMTNSX(NR,NS)
         ELSEIF(MODELP(NS).EQ.9) THEN
            MODELPS=MODELP(NS)
            MODELP(NS)=MODELPR(NR,NS)
            IF(MODELP(NS).EQ.8) THEN
               IF(MODELV(NS).EQ.9) THEN
                  MODELVS=MODELV(NS)
                  MODELV(NS)=MODELVR(NR,NS)
                  CALL WMDPIN(NR,NS)
                  MODELV(NS)=MODELVS
               ELSE
                  CALL WMDPIN(NR,NS)
               ENDIF
            ELSE
               CALL WMDPIN(NR,NS)
            ENDIF
            MODELP(NS)=MODELPS
         ELSE
            IF(MODELV(NS).EQ.5) THEN
               WRITE(6,*) NR,NS
               CALL WMDPDK(NR,NS)
            ELSE
               CALL WMDPIN(NR,NS)
            ENDIF
         ENDIF
      ENDIF
C
C      IF(NR.EQ.1) THEN
C      WRITE(6,*) 'WMDISP: NR,NS=',NR,NS
C      WRITE(6,'(1P6E12.4)') 
C     &     CTNSR(1,1,0,0,1,1),CTNSR(1,2,0,0,1,1),CTNSR(1,3,0,0,1,1),
C     &     CTNSR(2,1,0,0,1,1),CTNSR(2,2,0,0,1,1),CTNSR(2,3,0,0,1,1),
C     &     CTNSR(3,1,0,0,1,1),CTNSR(3,2,0,0,1,1),CTNSR(3,3,0,0,1,1)
C      ENDIF
C      IF(NR.EQ.2) STOP
C
      RETURN
      END
      SUBROUTINE WMTNSR_IPS(NR,NS,MD,ND)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcomm.inc'
C
C      DO NHH=1,NHHMAX
C      DO NTH=1,NTHMAX
C         DO ND=-NDSIZX,NDSIZX
C         DO MD=-MDSIZX,MDSIZX
!!!seki
      DO NHH=1,NHHMAX_IPS_F
      DO NTH=1,NTHMAX_IPS_F
!!!seki
            CTNSR_1(1,1,NTH,NHH)=0.D0
            CTNSR_1(1,2,NTH,NHH)=0.D0
            CTNSR_1(1,3,NTH,NHH)=0.D0
            CTNSR_1(2,1,NTH,NHH)=0.D0
            CTNSR_1(2,2,NTH,NHH)=0.D0
            CTNSR_1(2,3,NTH,NHH)=0.D0
            CTNSR_1(3,1,NTH,NHH)=0.D0
            CTNSR_1(3,2,NTH,NHH)=0.D0
            CTNSR_1(3,3,NTH,NHH)=0.D0
      ENDDO
      ENDDO
C
!      print *,"seki"
      IF((MOD(MODELA,2).EQ.1).AND.(NS.EQ.3)) THEN
         CALL WMTNAX(NR)
      ELSEIF((MOD(MODELA/2,2).EQ.1).AND.(NS.EQ.1)) THEN
         CALL WMTNEX(NR)
      ELSEIF(NS.EQ.5.OR.NS.EQ.6) THEN
         CALL WMTNDK(NR,NS)
      ELSE
         IF(MODELP(NS).LT.0) THEN
            CALL WMTNSX(NR,NS)
         ELSEIF(MODELP(NS).EQ.9) THEN
            MODELPS=MODELP(NS)
            MODELP(NS)=MODELPR(NR,NS)
            IF(MODELP(NS).EQ.8) THEN
               IF(MODELV(NS).EQ.9) THEN
                  MODELVS=MODELV(NS)
                  MODELV(NS)=MODELVR(NR,NS)
                  CALL WMDPIN_IPS(NR,NS,MD,ND)
                  MODELV(NS)=MODELVS
               ELSE
                  CALL WMDPIN_IPS(NR,NS,MD,ND)
               ENDIF
            ELSE
               CALL WMDPIN_IPS(NR,NS,MD,ND)
            ENDIF
            MODELP(NS)=MODELPS
         ELSE
            IF(MODELV(NS).EQ.5) THEN
               WRITE(6,*) NR,NS
               CALL WMDPDK(NR,NS)
            ELSE
               CALL WMDPIN_IPS(NR,NS,MD,ND)
            ENDIF
         ENDIF
      ENDIF
C
C      IF(NR.EQ.1) THEN
C      WRITE(6,*) 'WMDISP: NR,NS=',NR,NS
C      WRITE(6,'(1P6E12.4)') 
C     &     CTNSR(1,1,0,0,1,1),CTNSR(1,2,0,0,1,1),CTNSR(1,3,0,0,1,1),
C     &     CTNSR(2,1,0,0,1,1),CTNSR(2,2,0,0,1,1),CTNSR(2,3,0,0,1,1),
C     &     CTNSR(3,1,0,0,1,1),CTNSR(3,2,0,0,1,1),CTNSR(3,3,0,0,1,1)
C      ENDIF
C      IF(NR.EQ.2) STOP
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE WMTNSX(NR,NS)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      USE libdsp,ONLY: DSPFN
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
C
      CW=2.D0*PI*CRF*1.D6
      WW=DBLE(CW)
C
      CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
C
C      IF(NS.EQ.1.AND.NR.EQ.1) THEN
C         WRITE(6,'(A,1P6E12.4)') 'RN  :',(RN(NS1),NS1=1,NSMAX)
C         WRITE(6,'(A,1P6E12.4)') 'RTPR:',(RTPR(NS1),NS1=1,NSMAX)
C         WRITE(6,'(A,1P6E12.4)') 'RTPP:',(RTPP(NS1),NS1=1,NSMAX)
C      ENDIF
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
C         IF(NR.EQ.1.AND.NTH.EQ.1.AND.NHH.EQ.1) THEN
C            WRITE(6,*) 'BABS:',BABS,BSUPTH,BSUPPH
C         ENDIF
C
         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
C
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
C
C            WRITE(11,'(5I5,1P2E12.4)') NR,NN,MM,NHH,NTH,RKPR,BABS
C
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            RNPR=VC*RKPR/WW
C
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
               RKX2=     MM*MM*RG22(NTH,NHH,NR)*XRHO(NR)**2
     &             +2.D0*MM*NN*RG23(NTH,NHH,NR)*XRHO(NR)
     &             +     NN*NN*RG33(NTH,NHH,NR)
               IF(RKPP2.GT.0.D0) THEN
C                  RKPP=SQRT(RKPP2)
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
C
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
               CPARA=-   CWP/(CW*CW)
     &                   *(CW*CW/(CWU*CWN)
     &                    +2*RKPP2*RU(NS)*RU(NS)/(CWN*CWN-CWC*CWC)
     &                     *CWN/CWU)
               CPERM=(0.D0,0.D0)
C               WRITE(6,*) NR,NS,CPERP,CPARA
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
ccc                     CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS))
                     cadd=1+RKPR*RU(NS)*CW/(CW-RKPR*RU(NS))**2 !
ccc                     CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD*CADD
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
                     CPARA2=CPARA2
     &                     -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
                     CCROS1=CCROS1-CI*CWP*CGZ0*CZ/2
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*2*CZ
                  ELSEIF(NC.EQ.0) THEN
                     CPERM2=CPERM2+CWP*RKTPP*CGZ0*4*CZ
                     CPARA1=CPARA1-CWP*CDZ*CGZ*CGZ0*CADD*CADD
                     CPARA2=CPARA2
     &                     +CWP*RKTPP*CGZ0*2*CGZ*CDZ*CADD*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
                     CPARA2=CPARA2
     &                     -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
                     CCROS1=CCROS1+CI*CWP*CGZ0*CZ/2
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*2*CZ
                  ELSEIF(NC.EQ.2) THEN
                     CPERP2=CPERM2+CWP*RKTPP*CGZ0*CZ
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*CZ
                  ENDIF
               ENDDO
C
C               CPERP1= CWP*CGZ(0)*(CZ(1)+CZ(-1))/2
C               CPERP2=-CWP*RKTPP*CGZ(0)
C     &                    *(CZ(1)+CZ(-1)-CZ(2)-CZ(-2))    
C               CPERM2=-CWP*RKTPP*CGZ(0)
C     &                    *(2*CZ(1)+2*CZ(-1)-4*CZ(0))
C               CPARA1=-CWP*CDZ(0)*CGZ(0)*CGZ(0)*CADD(0)*CADD(0)
C               CPARA2= CWP*RKTPP*CGZ(0)
C     &                 *(2*CGZ( 0)*CDZ( 0)*CADD( 0)*CADD( 0)
C     &                    -CGZ( 1)*CDZ( 1)*CADD( 1)*CADD( 1)
C     &                    -CGZ(-1)*CDZ(-1)*CADD(-1)*CADD(-1))
C               CCROS1= CI*CWP*CGZ(0)*(CZ(1)-CZ(-1))/2
C               CCROS2=-CI*CWP*RKTPP*CGZ(0)
C     &                   *(2*CZ(1)-2*CZ(-1)-CZ(2)+CZ(-2))
C
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
C
C      IF(NR.EQ.1) THEN
C         WRITE(6,*) 'CPERP,CPERM,CPARA:',CPERP,CPERM,CPARA
C         WRITE(6,*) 'UXX2,UYY2:',UXX2,UYY2
C      ENDIF
C
C            CTNSR(1,1,MD,ND,NTH,NHH)= CPERP + UXX2*CPERM
C            CTNSR(1,2,MD,ND,NTH,NHH)=-CCROS
C            CTNSR(1,3,MD,ND,NTH,NHH)= 0.D0
C            CTNSR(2,1,MD,ND,NTH,NHH)= CCROS
C            CTNSR(2,2,MD,ND,NTH,NHH)= CPERP + UYY2*CPERM
C            CTNSR(2,3,MD,ND,NTH,NHH)= 0.D0
C            CTNSR(3,1,MD,ND,NTH,NHH)= 0.D0
C            CTNSR(3,2,MD,ND,NTH,NHH)= 0.D0
C            CTNSR(3,3,MD,ND,NTH,NHH)= CPARA
C
            CTNSR(1,1,MD,ND,NTH,NHH)
     &     =CTNSR(1,1,MD,ND,NTH,NHH) + CPERP + UXX2*CPERM
!     &     =CTNSR(1,1,MD,ND,NTH,NHH) + CPERP + UYY2*CPERM
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH) - CCROS
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH) + CCROS
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CPERP + UYY2*CPERM
!     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CPERP + UXX2*CPERM
            CTNSR(3,3,MD,ND,NTH,NHH)
     &     =CTNSR(3,3,MD,ND,NTH,NHH) + CPARA
C
C            IF(NR.EQ.30.AND.NTH.EQ.1.AND.MDX.EQ.1) THEN
C               WRITE(6,*) 'RN,RTPR,RTPP=',RN(1)/1.D20,
C     &                    RTPR(1)/(AEE*1.D3),RTPP(1)/(AEE*1.D3)
C               WRITE(6,*) 'BABS,QR=',BABS,QR
C               WRITE(6,*) 'RKPR,RKPP2=',RKPR,RKPP2
C               WRITE(6,*) 'CGZ(0)=',CGZ(0)
C               WRITE(6,*) 'CDZ(0)=',CDZ(0)
C               WRITE(6,*) 'CTNSR(3,3)=',CTNSR(3,3,NTH,MDX)
C            ENDIF
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** IMPORT FROM TASK/DP ******
C
CCCCCCCCCCCCCcseki取りあえずここのみ
      SUBROUTINE WMDPIN(NR,NS)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcoml.inc'
      DIMENSION CDTNS(3,3)
      INTEGER NDSIZX_TMP, MDSIZX_TMP
C
      CW=2.D0*PI*CRF*1.D6
      WW=DBLE(CW)
C
      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON)
!      IF (RHON > 1.0)then
!         CALL PL_PROF_OLD(1.0D0)
!C         RN(NS)=1d-6
!         RTPR(NS)=1d-6
!         RTPP(NS)=1d-6
!      ENDIF
C
C      IF(NS.EQ.1.AND.NR.EQ.1) THEN
C      IF(NS.EQ.1) THEN
C         WRITE(6,'(A,I5)')       'NR  :',NR
C         WRITE(6,'(A,1P6E12.4)') 'RN  :',(RN(NS1),NS1=1,NSMAX)
C         WRITE(6,'(A,1P6E12.4)') 'RTPR:',(RTPR(NS1),NS1=1,NSMAX)
C         WRITE(6,'(A,1P6E12.4)') 'RTPP:',(RTPP(NS1),NS1=1,NSMAX)
C      ENDIF
C
C      IF(NR.EQ.10) THEN
Ca         WRITE(6,'(A,1P3E12.4)') 'RN  :',RN(1),RN(2),RN(3)
C         WRITE(6,'(A,1P3E12.4)') 'RTPR:',RTPR(1),RTPR(2),RTPR(3)
C         WRITE(6,'(A,1P3E12.4)') 'RTPP:',RTPP(1),RTPP(2),RTPP(3)
C         WRITE(6,'(A,1P3E12.4)') 'RU  :',RU(1)  ,RU(2)  ,RU(3)
C      ENDIF
C
!      print *,"NR",NR,NS
!      if (ns == 3 )then
!           WRITE(112,*)"# NR",NR
!           WRITE(113,*)"# NR",NR
!      endif
      NDSIZX_TMP=NDSIZX_F
      MDSIZX_TMP=MDSIZX_F
      IF (NDSIZX_TMP==1)NDSIZX_TMP=0
      IF (MDSIZX_TMP==1)MDSIZX_TMP=0
!!!!!!!! seki
      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
!!!!!!!! seki
         CALL WMCMAG_F(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
C         CALL WMCMAG_IPS_F(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
C         IF(NR.EQ.10.AND.NTH.EQ.1.AND.NHH.EQ.1) THEN
C            WRITE(6,'(A,1P3E12.4)') 'BABS:',BABS,BSUPTH,BSUPPH
C         ENDIF
C
!!!!!!!! seki
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
!!!!!!!! seki
C         DO ND=0,0
            NN=NPH0+NHC*ND
!!!!!!!! seki
         DO MD=-MDSIZX_TMP,MDSIZX_TMP
!!!!!!!! seki
            MM=NTH0+MD
C
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
C            RKPR=PARAK(NTH,NHH,NR)
C            RKPR=NN/(2D0*3.14D0*(RR+RA*XRHO(NR)))
C
Cseki            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            IF(ABS(RKPR).LE.RKPR_EFF) RKPR=SIGN(RKPR_EFF,RKPR)
            RNPR=VC*RKPR/WW
C
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
C               RKX2=     MM*MM*RG22(NTH,NHH,NR)*XRHO(NR)**2
C     &             +2.D0*MM*NN*RG23(NTH,NHH,NR)*XRHO(NR)
C     &             +     NN*NN*RG33(NTH,NHH,NR)
C               RKX2=     MM*MM*RGI22(NTH,NHH,NR)/XRHO(NR)**2
C     &             +2.D0*MM*NN*RGI23(NTH,NHH,NR)/XRHO(NR)
C     &             +     NN*NN*RGI33(NTH,NHH,NR)
              IF (NR .eq. 1)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               RKX2=     MM*MM*RGI22(NTH,NHH,NR)*(1d6/XRHO(2))**2
     &             +2.D0*MM*NN*RGI23(NTH,NHH,NR)*(1d6/XRHO(2))
     &             +     NN*NN*RGI33(NTH,NHH,NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                 RKX2=  NN*NN*(RGI33(NTH,NHH,NR))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ELSE
               RKX2=     MM*MM*RGI22(NTH,NHH,NR)/XRHO(NR)**2
     &             +2.D0*MM*NN*RGI23(NTH,NHH,NR)/XRHO(NR)
     &             +     NN*NN*RGI33(NTH,NHH,NR)
              ENDIF
!              IF (NR .eq. 1)then
!                 RKX2= MM*MM/(RG22(NTH,NHH,2)*XRHO(2)**2)
!     &             +     NN*NN/(RG33(NTH,NHH,NR))
!                 IF (ABS(RG23(NTH,NHH,2)) .gt. 1d-15)THEN
!                    RKX2= RKX2
!     &                  +2.D0*MM*NN/(RG23(NTH,NHH,2)*XRHO(2))
!                 ENDIF
!               ELSE
!                 RKX2= MM*MM/(RG22(NTH,NHH,NR)*XRHO(NR)**2)
!     &             +     NN*NN/(RG33(NTH,NHH,NR))
!                 IF (ABS(RG23(NTH,NHH,NR)) .gt. 1d-15)THEN
!                    RKX2= RKX2
!     &                  +2.D0*MM*NN/(RG23(NTH,NHH,NR)*XRHO(NR))
!                 ENDIF
!               ENDIF
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
!                     RKPP=0d0
!                     RKR2=0d0
                  ENDIF
!                  UXX2= RKX2/RKPP2
                  UXX2= RKT2/RKPP2
                  UYY2= RKR2/RKPP2
               ELSE
!                  print *,"seki",RKPP2
!                  RKPP2=0.D0
                  RKPP=0.D0
                  RKT2=0.D0
                  RKR2=0.D0
                  UXX2=0.D0
                  UYY2=0.D0
!                  RKPP2=0.D0
!                  RKPP=0.D0
!                  RKT2=0.D0
!                  RKR2=0.D0
!                  UXX2=0.D0
!                  UYY2=0.D0
               ENDIF
            ELSE
               RKPP2=0.D0
               RKPP=0.D0
               RKT2=0.D0
               RKR2=0.D0
               UXX2=0.D0
               UYY2=0.D0
            ENDIF
C
C            RKPP=0.D0
C            UXX2=0.D0
C            UYY2=0.D0
C
            CKPR=RKPR
            CKPP=RKPP

            CALL DPCALC_2(CW,CKPR,CKPP,RHON,BABS,BTH
     &                     ,NS,CDTNS)
C            print *,"1"
C            print *,CDTNS
C            print *,"2"
C            print *, UXX2,UYY2
C            print *,"3"
C            print *,CKPR,CKPP

C
C      IF(NR.EQ.1.AND.
C     &   MD.EQ.0.AND.
C     &   ND.EQ.0.AND.
C     &   NTH.EQ.1.AND.
C     &   NHH.EQ.1) THEN
C         WRITE(6,'(A,2I5,1PE12.4)') 
C     &        'NS,MODELP,RHON',NS,MODELP(NS),RHON
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CW,R,P=',CW,CKPR,CKRR
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CDTNS1=',CDTNS(1,1),CDTNS(1,2),CDTNS(1,3)
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CDTNS2=',CDTNS(2,1),CDTNS(2,2),CDTNS(2,3)
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CDTNS3=',CDTNS(3,1),CDTNS(3,2),CDTNS(3,3)
C      ENDIF
C
            CPERM=CDTNS(2,2)-CDTNS(1,1)
C
            CTNSR(1,1,MD,ND,NTH,NHH)
     &     =CTNSR(1,1,MD,ND,NTH,NHH) + CDTNS(1,1) + UXX2*CPERM
!     &     =CTNSR(1,1,MD,ND,NTH,NHH) + CDTNS(1,1) + UYY2*CPERM
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH) + CDTNS(1,2)
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH) + CDTNS(2,1)
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CDTNS(1,1) + UYY2*CPERM
!     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CDTNS(1,1) + UXX2*CPERM
            CTNSR(3,3,MD,ND,NTH,NHH)
     &     =CTNSR(3,3,MD,ND,NTH,NHH) + CDTNS(3,3)

C
C            CTNSR(1,1,MD,ND,NTH,NHH)
C     &     =CTNSR(1,1,MD,ND,NTH,NHH) + CDTNS(1,1) + UXX2*CPERM
C            CTNSR(1,2,MD,ND,NTH,NHH)
C     &     =CTNSR(1,2,MD,ND,NTH,NHH) + CDTNS(1,2)
C            CTNSR(1,3,MD,ND,NTH,NHH)
C     &     =CTNSR(1,3,MD,ND,NTH,NHH) + CDTNS(1,3)
C            CTNSR(2,1,MD,ND,NTH,NHH)
C     &     =CTNSR(2,1,MD,ND,NTH,NHH) + CDTNS(2,1)
C            CTNSR(2,2,MD,ND,NTH,NHH)
C     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CDTNS(2,2) + UYY2*CPERM
C            CTNSR(2,3,MD,ND,NTH,NHH)
C     &     =CTNSR(2,3,MD,ND,NTH,NHH) + CDTNS(2,3)
C            CTNSR(3,1,MD,ND,NTH,NHH)
C     &     =CTNSR(3,1,MD,ND,NTH,NHH) + CDTNS(3,1)
C            CTNSR(3,2,MD,ND,NTH,NHH)
C     &     =CTNSR(3,2,MD,ND,NTH,NHH) + CDTNS(3,2)
C            CTNSR(3,3,MD,ND,NTH,NHH)
C     &     =CTNSR(3,3,MD,ND,NTH,NHH) + CDTNS(3,3)
C
C            IF(NR.EQ.30.AND.NTH.EQ.1.AND.MDX.EQ.1) THEN
C               WRITE(6,*) 'RN,RTPR,RTPP=',RN(1)/1.D20,
C     &                    RTPR(1)/(AEE*1.D3),RTPP(1)/(AEE*1.D3)
C               WRITE(6,*) 'BABS,QR=',BABS,QR
C               WRITE(6,*) 'RKPR,RKPP2=',RKPR,RKPP2
C               WRITE(6,*) 'CGZ(0)=',CGZ(0)
C               WRITE(6,*) 'CDZ(0)=',CDZ(0)
C               WRITE(6,*) 'CTNSR(3,3)=',CTNSR(3,3,NTH,MDX)
C            ENDIF
         ENDDO
         ENDDO
!         if (ns == 3 )then
!         DO MD=-MDSIZX,MDSIZX
!            MM=NTH0+MD
!           WRITE(112,"(2(i4,1x),10(e15.6,1x))")MM,NTH,
!     &                    real(CTNSR(1,1,MD,0,NTH,NHH)),
!     &                    real(CTNSR(1,2,MD,0,NTH,NHH)),
!     &                    real(CTNSR(2,1,MD,0,NTH,NHH)),
!     &                    real(CTNSR(2,2,MD,0,NTH,NHH)),
!     &                    real(CTNSR(3,3,MD,0,NTH,NHH))
!           WRITE(113,"(2(i4,1x),10(e15.6,1x))")MM,NTH,
!     &                    imag(CTNSR(1,1,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(1,2,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(2,1,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(2,2,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(3,3,MD,0,NTH,NHH))
!         ENDDO
!           WRITE(112,"(i4,1x,(e15.6,1x))")
!           WRITE(113,"(i4,1x,(e15.6,1x))")
!          endif
      ENDDO
      ENDDO
!!         if (ns == 3 .and. nr==20)then
!         if (ns == 3 )then
!           WRITE(112,"(i4,1x,(e15.6,1x))")
!           WRITE(112,"(i4,1x,(e15.6,1x))")
!           WRITE(113,"(i4,1x,(e15.6,1x))")
!           WRITE(113,"(i4,1x,(e15.6,1x))")
!          endif
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc

CCCCCCCCCCCCCcseki取りあえずここのみ
      SUBROUTINE WMDPIN_IPS(NR,NS,MD,ND)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcoml.inc'
      DIMENSION CDTNS(3,3)
      INTEGER NDSIZX_TMP, MDSIZX_TMP
C
      CW=2.D0*PI*CRF*1.D6
      WW=DBLE(CW)
C
      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON)
!      IF (RHON > 1.0)then
!         CALL PL_PROF_OLD(1.0D0)
!C         RN(NS)=1d-6
!         RTPR(NS)=1d-6
!         RTPP(NS)=1d-6
!      ENDIF
C
C      IF(NS.EQ.1.AND.NR.EQ.1) THEN
C      IF(NS.EQ.1) THEN
C         WRITE(6,'(A,I5)')       'NR  :',NR
C         WRITE(6,'(A,1P6E12.4)') 'RN  :',(RN(NS1),NS1=1,NSMAX)
C         WRITE(6,'(A,1P6E12.4)') 'RTPR:',(RTPR(NS1),NS1=1,NSMAX)
C         WRITE(6,'(A,1P6E12.4)') 'RTPP:',(RTPP(NS1),NS1=1,NSMAX)
C      ENDIF
C
C      IF(NR.EQ.10) THEN
Ca         WRITE(6,'(A,1P3E12.4)') 'RN  :',RN(1),RN(2),RN(3)
C         WRITE(6,'(A,1P3E12.4)') 'RTPR:',RTPR(1),RTPR(2),RTPR(3)
C         WRITE(6,'(A,1P3E12.4)') 'RTPP:',RTPP(1),RTPP(2),RTPP(3)
C         WRITE(6,'(A,1P3E12.4)') 'RU  :',RU(1)  ,RU(2)  ,RU(3)
C      ENDIF
C
!      print *,"NR",NR,NS
!      if (ns == 3 )then
!           WRITE(112,*)"# NR",NR
!           WRITE(113,*)"# NR",NR
!      endif
!      NDSIZX_TMP=NDSIZX_F
!      MDSIZX_TMP=MDSIZX_F
!      IF (NDSIZX_TMP==1)NDSIZX_TMP=0
!      IF (MDSIZX_TMP==1)MDSIZX_TMP=0
!!!!!!!! seki
      DO NHH=1,NHHMAX_IPS_F
!$OMP PARALLEL DO DEFAULT(shared)
!$OMP+private(BABS,BSUPTH,BSUPPH)
!$OMP+private(BTH,WTPR,RKPR_EFF,RKTH,RKPH,RKPR)
!$OMP+private(RNPR,NSS,AM,AE,WP,WC,DTT,DTX)
!$OMP+private(RNPP2,RKPP2,RKX2,CDTNS)
!$OMP+private(RKT2,RKR2,UXX2,UYY2)
!$OMP+private(CKPR,CKPP,CPERM)
!$OMP+private(CKPPF,CKPPS)
!$OMP+private(MM,NN)
      DO NTH=1,NTHMAX_IPS_F
!!!!!!!! seki
C         CALL WMCMAG_F(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         CALL WMCMAG_IPS_F(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
C         IF(NR.EQ.10.AND.NTH.EQ.1.AND.NHH.EQ.1) THEN
C            WRITE(6,'(A,1P3E12.4)') 'BABS:',BABS,BSUPTH,BSUPPH
C         ENDIF
C
!!!!!!!! seki
!         BTH=ABS(BSUPTH*RA)
         BTH=ABS(BSUPTH*RA*RHON)
         WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
         RKPR_EFF=SQRT(CW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
!         RKPR_EFF=SQRT(CW*BTHOBN(NR)/SQRT(8.D0*WTPR)/RR)
         IF (RKPR_EFF < 1d-5)then
           RKPR_EFF=1.D-5
!         IF (RKPR_EFF < 1d1*RHON)then
!           RKPR_EFF=1.D1*RHON
!           if (NR .NE. 1)then 
!             print *,RKPR_EFF,BTH,WKPT
!             STOP
!           ENDIF
         endif
C         DO ND=-NDSIZX_TMP,NDSIZX_TMP
!!!!!!!! seki
C         DO ND=0,0
            NN=NPH0+NHC*ND
!!!!!!!! seki
C         DO MD=-MDSIZX_TMP,MDSIZX_TMP
!!!!!!!! seki
            MM=NTH0+MD
C
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
C            RKPR=PARAK(NTH,NHH,NR)
C            RKPR=NN/(2D0*3.14D0*(RR+RA*XRHO(NR)))
C
Cseki            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            IF(ABS(RKPR).LE.RKPR_EFF)RKPR=SIGN(RKPR_EFF,RKPR)
            RNPR=VC*RKPR/WW
C
            CKPPF=0d0
C            IF(MM .eq. 0 .AND. NN .eq. 0)cycle
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
               CKPR=RKPR
               CKPPF=RKPP

               call DPCOLD_RKPERP_2(CW,CKPR,CKPPF,CKPPS,BABS,BTH)
               RKPP2=REAL(CKPPF**2)
    

C               RKX2=     MM*MM*RG22_IPS(NTH,NHH,NR)*XRHO(NR)**2
C     &             +2.D0*MM*NN*RG23_IPS(NTH,NHH,NR)*XRHO(NR)
C     &             +     NN*NN*RG33_IPS(NTH,NHH,NR)
C               RKX2=     MM*MM*RGI22_IPS(NTH,NHH,NR)/XRHO(NR)**2
C     &             +2.D0*MM*NN*RGI23_IPS(NTH,NHH,NR)/XRHO(NR)
C     &             +     NN*NN*RGI33_IPS(NTH,NHH,NR)
              IF (NR .eq. 1)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               RKX2=     MM*MM*RGI22_IPS(NTH,NHH,NR)*(1d6/XRHO(2))**2
     &             +2.D0*MM*NN*RGI23_IPS(NTH,NHH,NR)*(1d6/XRHO(2))
     &             +     NN*NN*RGI33_IPS(NTH,NHH,NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                 RKX2=  NN*NN*(RGI33_IPS(NTH,NHH,NR))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ELSE
               RKX2=     MM*MM*RGI22_IPS(NTH,NHH,NR)/XRHO(NR)**2
     &             +2.D0*MM*NN*RGI23_IPS(NTH,NHH,NR)/XRHO(NR)
     &             +     NN*NN*RGI33_IPS(NTH,NHH,NR)
              ENDIF
!              IF (NR .eq. 1)then
!                 RKX2= MM*MM/(RG22_IPS(NTH,NHH,2)*XRHO(2)**2)
!     &             +     NN*NN/(RG33_IPS(NTH,NHH,NR))
!                 IF (ABS(RG23_IPS(NTH,NHH,2)) .gt. 1d-15)THEN
!                    RKX2= RKX2
!     &                  +2.D0*MM*NN/(RG23_IPS(NTH,NHH,2)*XRHO(2))
!                 ENDIF
!                 RKX2=  NN*NN/(RG33_IPS(NTH,NHH,NR))
!               ELSE
!                 RKX2= MM*MM/(RG22_IPS(NTH,NHH,NR)*XRHO(NR)**2)
!     &             +     NN*NN/(RG33_IPS(NTH,NHH,NR))
!                 IF (ABS(RG23_IPS(NTH,NHH,NR)) .gt. 1d-15)THEN
!                    RKX2= RKX2
!     &                  +2.D0*MM*NN/(RG23_IPS(NTH,NHH,NR)*XRHO(NR))
!                 ENDIF
!               ENDIF

               IF(RKPP2.GT.0.D0) THEN
                  RKPP=SQRT(RKPP2)
                  RKT2=RKX2-RKPR**2
!                  print *,RKX2,RKPR**2,RKT2
                  IF(RKT2.GT.0.D0) THEN
                     IF(RKT2.LE.RKPP2) THEN
                        RKR2=RKPP2-RKT2
                     ELSE
                        RKR2=0.D0
!                        RKT2=RKPP2
                     ENDIF
                  ELSE
                     RKT2=0.D0
                     RKR2=RKPP2
!                     RKPP=0d0
!                     RKR2=0d0
                  ENDIF
!                  UXX2= RKX2/RKPP2
                  UXX2= RKT2/RKPP2
                  UYY2= RKR2/RKPP2
               ELSE
!                  print *,"seki2",XRHO(NR),RKPP2
!                  RKPP2=0.D0
                  RKPP=0.D0
                  RKT2=0.D0
                  RKR2=0.D0
                  UXX2=0.D0
                  UYY2=0.D0
!                  RKPP2=0.D0
!                  RKPP=0.D0
!                  RKT2=0.D0
!                  RKR2=0.D0
!                  UXX2=0.D0
!                  UYY2=0.D0
                   CKPPF=RKPP
               ENDIF
            ELSE
               RKPP2=0.D0
               RKPP=0.D0
               RKT2=0.D0
               RKR2=0.D0
               UXX2=0.D0
               UYY2=0.D0
               CKPPF=RKPP
            ENDIF
C
C            RKPP=0.D0
C            UXX2=0.D0
C            UYY2=0.D0
C
            CKPR=RKPR
!            CKPP=(0d0,0d0) !RKPP
!            CKPR=(0d0,0d0) !RKPP

!            CKPP=RKPP2
!            CKPP=sqrt(CKPP)

!            CKPP=RKPP
            CKPP=CKPPF

!            CALL DPCALC(CW,CKPR,CKPP,RHON,BABS,NS,CDTNS)
!            print *,CW,CKPR,CKPP,RHON,BABS,BSUPTH*RA,NS
!            print "(10(es10.3,1x))",XRHO(NR), CKPP,CKPR,UXX2,UYY2
            CALL DPCALC_2(CW,CKPR,CKPP,RHON,BABS,BTH
     &                     ,NS,CDTNS)
C            print *,"1"
C            print *,CDTNS
C            print *,"2"
C            print *, UXX2,UYY2
C            print *,"3"
C            print *,CKPR,CKPP

C
C      IF(NR.EQ.1.AND.
C     &   MD.EQ.0.AND.
C     &   ND.EQ.0.AND.
C     &   NTH.EQ.1.AND.
C     &   NHH.EQ.1) THEN
C         WRITE(6,'(A,2I5,1PE12.4)') 
C     &        'NS,MODELP,RHON',NS,MODELP(NS),RHON
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CW,R,P=',CW,CKPR,CKRR
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CDTNS1=',CDTNS(1,1),CDTNS(1,2),CDTNS(1,3)
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CDTNS2=',CDTNS(2,1),CDTNS(2,2),CDTNS(2,3)
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CDTNS3=',CDTNS(3,1),CDTNS(3,2),CDTNS(3,3)
C      ENDIF
C
            CPERM=CDTNS(2,2)-CDTNS(1,1)
C
            CTNSR_1(1,1,NTH,NHH)
     &     =CTNSR_1(1,1,NTH,NHH) + CDTNS(1,1)  + UXX2*CPERM
!     &     =CTNSR_1(1,1,NTH,NHH) + CDTNS(1,1)  + UYY2*CPERM
            CTNSR_1(1,2,NTH,NHH)
     &     =CTNSR_1(1,2,NTH,NHH) + CDTNS(1,2)
            CTNSR_1(2,1,NTH,NHH)
     &     =CTNSR_1(2,1,NTH,NHH) + CDTNS(2,1)
            CTNSR_1(2,2,NTH,NHH)
     &     =CTNSR_1(2,2,NTH,NHH) + CDTNS(1,1)  + UYY2*CPERM
!     &     =CTNSR_1(2,2,NTH,NHH) + CDTNS(1,1)  + UXX2*CPERM
            CTNSR_1(3,3,NTH,NHH)
     &     =CTNSR_1(3,3,NTH,NHH) + CDTNS(3,3)

C
C            CTNSR(1,1,MD,ND,NTH,NHH)
C     &     =CTNSR(1,1,MD,ND,NTH,NHH) + CDTNS(1,1) + UXX2*CPERM
C            CTNSR(1,2,MD,ND,NTH,NHH)
C     &     =CTNSR(1,2,MD,ND,NTH,NHH) + CDTNS(1,2)
C            CTNSR(1,3,MD,ND,NTH,NHH)
C     &     =CTNSR(1,3,MD,ND,NTH,NHH) + CDTNS(1,3)
C            CTNSR(2,1,MD,ND,NTH,NHH)
C     &     =CTNSR(2,1,MD,ND,NTH,NHH) + CDTNS(2,1)
C            CTNSR(2,2,MD,ND,NTH,NHH)
C     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CDTNS(2,2) + UYY2*CPERM
C            CTNSR(2,3,MD,ND,NTH,NHH)
C     &     =CTNSR(2,3,MD,ND,NTH,NHH) + CDTNS(2,3)
C            CTNSR(3,1,MD,ND,NTH,NHH)
C     &     =CTNSR(3,1,MD,ND,NTH,NHH) + CDTNS(3,1)
C            CTNSR(3,2,MD,ND,NTH,NHH)
C     &     =CTNSR(3,2,MD,ND,NTH,NHH) + CDTNS(3,2)
C            CTNSR(3,3,MD,ND,NTH,NHH)
C     &     =CTNSR(3,3,MD,ND,NTH,NHH) + CDTNS(3,3)
C
C            IF(NR.EQ.30.AND.NTH.EQ.1.AND.MDX.EQ.1) THEN
C               WRITE(6,*) 'RN,RTPR,RTPP=',RN(1)/1.D20,
C     &                    RTPR(1)/(AEE*1.D3),RTPP(1)/(AEE*1.D3)
C               WRITE(6,*) 'BABS,QR=',BABS,QR
C               WRITE(6,*) 'RKPR,RKPP2=',RKPR,RKPP2
C               WRITE(6,*) 'CGZ(0)=',CGZ(0)
C               WRITE(6,*) 'CDZ(0)=',CDZ(0)
C               WRITE(6,*) 'CTNSR(3,3)=',CTNSR(3,3,NTH,MDX)
C            ENDIF
C         ENDDO
C         ENDDO
!         if (ns == 3 )then
!         DO MD=-MDSIZX,MDSIZX
!            MM=NTH0+MD
!           WRITE(112,"(2(i4,1x),10(e15.6,1x))")MM,NTH,
!     &                    real(CTNSR(1,1,MD,0,NTH,NHH)),
!     &                    real(CTNSR(1,2,MD,0,NTH,NHH)),
!     &                    real(CTNSR(2,1,MD,0,NTH,NHH)),
!     &                    real(CTNSR(2,2,MD,0,NTH,NHH)),
!     &                    real(CTNSR(3,3,MD,0,NTH,NHH))
!           WRITE(113,"(2(i4,1x),10(e15.6,1x))")MM,NTH,
!     &                    imag(CTNSR(1,1,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(1,2,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(2,1,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(2,2,MD,0,NTH,NHH)),
!     &                    imag(CTNSR(3,3,MD,0,NTH,NHH))
!         ENDDO
!           WRITE(112,"(i4,1x,(e15.6,1x))")
!           WRITE(113,"(i4,1x,(e15.6,1x))")
!          endif
      ENDDO
!$OMP END PARALLEL DO
      ENDDO
!!         if (ns == 3 .and. nr==20)then
!         if (ns == 3 )then
!           WRITE(112,"(i4,1x,(e15.6,1x))")
!           WRITE(112,"(i4,1x,(e15.6,1x))")
!           WRITE(113,"(i4,1x,(e15.6,1x))")
!           WRITE(113,"(i4,1x,(e15.6,1x))")
!          endif
C
      RETURN
      END
C
C     ****** IMPORT FROM TASK/DP ******
C
      SUBROUTINE WMDPDK(NR,NS)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcoml.inc'
      DIMENSION CDTNS(3,3)
C
      CW=2.D0*PI*CRF*1.D6
C
      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON)
      IF(RN(NS).EQ.0.D0) RETURN
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)

         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
C
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
C
            CKPR=RKPR

            CALL DPDKDT(CW,CKPR,NS,NR,NTH,NTH,MM,CDTNS)

            DO j=1,3
            DO i=1,3
               CTNSR(i,j,MD,ND,NTH,NHH)
     &              =CTNSR(i,j,MD,ND,NTH,NHH)+CDTNS(i,j)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******
C
      SUBROUTINE WMTNAX(NR)
C
      INCLUDE 'wmcomm.inc'
C
      CW=2.D0*PI*CRF*1.D6
C
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
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
C         ANGPH=(NHH-1)*2.D0*PI/NHHMAX
         WC=AE*BABS/AM
         RHOA=VTA/WC
         RHOR=RHOA/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTA*DRN/(BABS*CW*AE*XL)
         ENDIF
C
         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
C
            CX=CW/(ABS(RKPR)*VTA)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0
C
C            CTNSR(1,1,MD,ND,NTH,NHH)=CPM*COS(ANGTH)**2
C            CTNSR(1,2,MD,ND,NTH,NHH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(1,3,MD,ND,NTH,NHH)=CQM*COS(ANGTH)
C            CTNSR(2,1,MD,ND,NTH,NHH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(2,2,MD,ND,NTH,NHH)=CPM*SIN(ANGTH)**2
C            CTNSR(2,3,MD,ND,NTH,NHH)=CQM*SIN(ANGTH)
C            CTNSR(3,1,MD,ND,NTH,NHH)=CQM*COS(ANGTH)
C            CTNSR(3,2,MD,ND,NTH,NHH)=CQM*SIN(ANGTH)
C            CTNSR(3,3,MD,ND,NTH,NHH)=CRM
C
            CTNSR(1,1,MD,ND,NTH,NHH)
     &     =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH)
     &     =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH)
     &     =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH)
     &     =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH)
     &     =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NHH)
     &     =CTNSR(3,3,MD,ND,NTH,NHH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
      END
C
C     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******
C
      SUBROUTINE WMTNEX(NR)
C
      INCLUDE 'wmcomm.inc'
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
C
      CW=2.D0*PI*CRF*1.D6
C
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
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
C         ANGPH=(NHH-1)*2.D0*PI/NHHMAX
         WC=AE*BABS/AM
         RHOE=VTE/WC
         RHOR=RHOE/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTE*DRN/(BABS*CW*AE*XL)
         ENDIF
C
         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
C
            CX=CW/(ABS(RKPR)*VTE)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0
C
C            CTNSR(1,1,MD,ND,NTH,NHH)=CPM*COS(ANGTH)**2
C            CTNSR(1,2,MD,ND,NTH,NHH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(1,3,MD,ND,NTH,NHH)=CQM*COS(ANGTH)
C            CTNSR(2,1,MD,ND,NTH,NHH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(2,2,MD,ND,NTH,NHH)=CPM*SIN(ANGTH)**2
C            CTNSR(2,3,MD,ND,NTH,NHH)=CQM*SIN(ANGTH)
C            CTNSR(3,1,MD,ND,NTH,NHH)=CQM*COS(ANGTH)
C            CTNSR(3,2,MD,ND,NTH,NHH)=CQM*SIN(ANGTH)
C            CTNSR(3,3,MD,ND,NTH,NHH)=CRM
C
            CTNSR(1,1,MD,ND,NTH,NHH)
     &     =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH)
     &     =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH)
     &     =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH)
     &     =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH)
     &     =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NHH)
     &     =CTNSR(3,3,MD,ND,NTH,NHH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
      END
C
C     ****** CALCULATE MAGNETIC DRIFT DIELECTRIC TENSOR ******
C
      SUBROUTINE WMTNDK(NR,NS)
C
      INCLUDE 'wmcoml.inc'
C
      CW=2.D0*PI*CRF*1.D6
C
      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON+0.01D0)
      RNAP=RN(NS)*1.D20
      RTAP=(RTPR(NS)+RTPP(NS))/3.D0*1.D3*AEE
      CALL PL_PROF_OLD(RHON)
      RNA=RN(NS)*1.D20
      RTA=(RTPR(NS)+RTPP(NS))/3.D0*1.D3*AEE
      DRN=(RNAP-RNA)/(0.01D0*RA*RNA)
C      WRITE(6,'(I5,1P3E12.4)') NR,RNA,RTA,DRN
C
      XL=RHON*RA
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
C         ANGPH=(NHH-1)*2.D0*PI/NHHMAX
         WC=AE*BABS/AM
         RHOA=VTA/WC
         RHOR=RHOA/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTA*DRN/(BABS*CW*AE*XL)
         ENDIF
C
         DO ND=-NDSIZX,NDSIZX
            NN=NPH0+NHC*ND
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
C
            CX=CW/(ABS(RKPR)*VTA)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0
C
C            CTNSR(1,1,MD,ND,NTH,NHH)=CPM*COS(ANGTH)**2
C            CTNSR(1,2,MD,ND,NTH,NHH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(1,3,MD,ND,NTH,NHH)=CQM*COS(ANGTH)
C            CTNSR(2,1,MD,ND,NTH,NHH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(2,2,MD,ND,NTH,NHH)=CPM*SIN(ANGTH)**2
C            CTNSR(2,3,MD,ND,NTH,NHH)=CQM*SIN(ANGTH)
C            CTNSR(3,1,MD,ND,NTH,NHH)=CQM*COS(ANGTH)
C            CTNSR(3,2,MD,ND,NTH,NHH)=CQM*SIN(ANGTH)
C            CTNSR(3,3,MD,ND,NTH,NHH)=CRM
C
            CTNSR(1,1,MD,ND,NTH,NHH)
     &     =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH)
     &     =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH)
     &     =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH)
     &     =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH)
     &     =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NHH)
     &     =CTNSR(3,3,MD,ND,NTH,NHH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE WMINIKPARA
C
      INCLUDE 'wmcoml.inc'
      complex*8 CEF1,CEF2,CEF3
      complex*8 CPHI,CTHETA
      READ(1122)CEFLDK
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/NHHMAX

      DO NR =1,NRMAX
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
C         IF(NR.EQ.10.AND.NTH.EQ.1.AND.NHH.EQ.1) THEN
C            WRITE(6,'(A,1P3E12.4)') 'BABS:',BABS,BSUPTH,BSUPPH
C         ENDIF
C
         THETA=DTH*NTH
         PHI  =DPHI*NHH
         ARKPR=0d0
         AREF=0d0
         ARKPR1=0d0
         AREF1=0d0
         ARKPR2=0d0
         AREF2=0d0
         ARKPR3=0d0
         AREF3=0d0
         PARAK(NTH,NHH,NR)=0d0
C         DO ND=-NDSIZX,NDSIZX
         DO ND=0,0
            NN=NPH0+NHC*ND
            NDV=ND + NDSIZX+1
            CPHI=cos(NN*phi)+ci*sin(NN*phi)
         DO MD=-MDSIZX,MDSIZX
            MM=NTH0+MD
            MDV=MD + MDSIZX+1
            CTHETA=cos(mm*theta)+ci*sin(mm*theta)
            CEF1=CEFLDK(1,MDV,NDV,NR)*CPHI*CTHETA/(XRHO(NR)+ 0.0001d0)
            CEF2=CEFLDK(2,MDV,NDV,NR)*CPHI*CTHETA
            CEF3=CEFLDK(3,MDV,NDV,NR)*CPHI*CTHETA
!            REF=sqrt((real(CEF1)**2 + real(CEF2)**2+real(CEF3)**2))
!            REF=(real(CEF1) + real(CEF2)+real(CEF3))
            REF=sqrt((CEF1)*conjg(CEF1)
     &              +(CEF2)*conjg(CEF2)
     &              +(CEF3)*conjg(CEF3))
C
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            ARKPR=ARKPR + RKPR*REF
            AREF=AREF + REF
            ARKPR1=ARKPR1 + RKPR*REAL(CEF1)
            AREF1=AREF1 +  REAL(CEF1)
            ARKPR2=ARKPR2 + RKPR*REAL(CEF2)
            AREF2=AREF2 +  REAL(CEF2)
            ARKPR3=ARKPR3 + RKPR*REAL(CEF3)
            AREF3=AREF3 +  REAL(CEF3)
          ENDDO
          ENDDO
          print *,NR,NTH,NHH
!          print *,ARKPR1/(AREF1+1d-36)
!          print *,ARKPR2/(AREF2+1d-36)
!          print *,ARKPR3/(AREF3+1d-36)
          PARAK(NTH,NHH,NR)=ARKPR/(AREF+1d-36)
!          PARAK(NTH,NHH,NR)=max(ARKPR1/(AREF1+1d-36),
!     &                          ARKPR2/(AREF2+1d-36),
!     &                          ARKPR3/(AREF3+1d-36))
          print *,PARAK(NTH,NHH,NR)
      ENDDO
      ENDDO
      ENDDO

      END
