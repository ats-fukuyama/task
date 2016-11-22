C     $Id$
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
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         DO ND=NDMIN,NDMAX
         DO MD=MDMIN,MDMAX
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
      IF((MOD(MODELA,2).EQ.1).AND.(NS.EQ.3)) THEN
         CALL WMTNAX(NR)
      ELSEIF((MOD(MODELA/2,2).EQ.1).AND.(NS.EQ.1)) THEN
         CALL WMTNEX(NR)
!      ELSEIF(NS.EQ.5.OR.NS.EQ.6) THEN
!         CALL WMTNDK(NR,NS)
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
         DO ND=NDMIN,NDMAX
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
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
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH) - CCROS
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH) + CCROS
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CPERP + UYY2*CPERM
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
      SUBROUTINE WMDPIN(NR,NS)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcoml.inc'
      DIMENSION CDTNS(3,3)
C
      CW=2.D0*PI*CRF*1.D6
      WW=DBLE(CW)
C
      RHON=XRHO(NR)
      CALL PL_PROF_OLD(RHON)
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
C         WRITE(6,'(A,1P3E12.4)') 'RN  :',RN(1),RN(2),RN(3)
C         WRITE(6,'(A,1P3E12.4)') 'RTPR:',RTPR(1),RTPR(2),RTPR(3)
C         WRITE(6,'(A,1P3E12.4)') 'RTPP:',RTPP(1),RTPP(2),RTPP(3)
C         WRITE(6,'(A,1P3E12.4)') 'RU  :',RU(1)  ,RU(2)  ,RU(3)
C      ENDIF
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
C         IF(NR.EQ.10.AND.NTH.EQ.1.AND.NHH.EQ.1) THEN
C            WRITE(6,'(A,1P3E12.4)') 'BABS:',BABS,BSUPTH,BSUPPH
C         ENDIF
C
         DO ND=NDMIN,NDMAX
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
            MM=NTH0+MD
C
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
C
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
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
               RKX2=     MM*MM*RG22(NTH,NHH,NR)*XRHO(NR)**2
     &             +2.D0*MM*NN*RG23(NTH,NHH,NR)*XRHO(NR)
     &             +     NN*NN*RG33(NTH,NHH,NR)
               IF(RKPP2.GT.0.D0) THEN
                  RKPP=SQRT(RKPP2)
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
C            RKPP=0.D0
C            UXX2=0.D0
C            UYY2=0.D0
C
            CKPR=RKPR
            CKPP=RKPP
            CALL DPCALC(CW,CKPR,CKPP,RHON,BABS,NS,CDTNS)
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
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH) + CDTNS(1,2)
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH) + CDTNS(2,1)
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH) + CDTNS(1,1) + UYY2*CPERM
            CTNSR(3,3,MD,ND,NTH,NHH)
     &     =CTNSR(3,3,MD,ND,NTH,NHH) + CDTNS(3,3)
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

         DO ND=NDMIN,NDMAX
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
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
!      IMPLICIT NONE
      COMPLEX*16 :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2 
!      COMPLEX*16 :: CXM(1:6)
      DOUBLE PRECISION :: RHOL,WP02,AE2N0
      INTEGER :: NS,NRWM
      DOUBLE PRECISION :: DELRWM,RL
C
!      WRITE(6,'(A,I5)') 'NR=',NR
      NS=3
      CW=2.D0*PI*CRF*1.D6
C
      CALL WMCPOS(NR,XL)
      IF(XL.LT.RA) THEN
!         RNA=PNA*EXP(-(XL/PNAL)**2)*1.D20    !gaussian
!         DRN=-2.D0*XL/(PNAL**2)              !gaussian
!          RNA=(PNA-1.D-5)*1.D20*(1.D0-XL**2)+1.D-5*1.D20 ! PNA*1.D20
!          DRN=-2.D0*XL*(1.D0-1.D-5*1.D20/RNA)/(1.D0-XL**2) !0.D0
          RNA=(PNA-1.D-5)*1.D20*((1.00001D0-XL**2)**0.5D0)+ 1.D-5*1.D20    !parabola ()^1 or ()^1/2
          DRN=-XL*(1.D0-1.D-5*1.D20/RNA)/(1.00001D0-XL**2) !-XL/((1.00001D0-XL**2)**0.5D0) !parabola
      ELSE
         RETURN
      ENDIF
      IF(XL.EQ.0.D0) THEN
         RHOL=RB/RA*(NR-0.5D0)/NRMAX
      ELSE
         RHOL=XL                ! RB/RA*(NR-0.5D0)/NRMAX
      ENDIF
!      RTA=PTA*AEE*1.D3          ! original
      RTA=(PTA-5.D-1)*(1.D0-XL**2)*AEE*1.D3 +5.D0*AEE*1.D3  ! parabola
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
         DO ND=NDMIN,NDMAX
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
C
            CX=CW/(ABS(RKPR)*VTA)

            IF(MODEFA.EQ.0) THEN
               IF(ABS(CX).GT.5.D0) THEN
                  CEX=(0.D0,0.D0)
               ELSE
                  CEX=CX*EXP(-CX*CX)
               ENDIF
               CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+CX*CX*(2.D0+CX*CX))
               CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
               CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0
            ELSEIF((MODEFA.EQ.1).OR.(MODEFA.EQ.2)) THEN 

      CALL WMDPFA(CX,CFN,COEF,RHOR,CPM,CQM,CRM,MODEFA)

            ELSEIF(MODEFA.EQ.3) THEN 

               RHOL=1.2D0*(NR-0.5D0)/NRMAX  
!      CALL WMDPFA2(NS,CW,RHOL,RKPR,VTA,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)  !normarized p,theta
      CALL WMDPFA2(CW,AM,RHOL,RKPR,VTA,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2) !not normarized p,theta

!             WRITE(6,'(i5,1P2E12.4)') NS,RHOL,RKPR
!             WRITE(6,'(1P6E12.4)') CPM1,CPM2,CQM1
!             WRITE(6,'(1P6E12.4)') CQM2,CRM1,CRM2
             CPM=5.D-1*PI*WP2/(CW*CW*WC*WC*RR*RR)*
     &           (-CPM1+2.D0*CPM2*MM/(AM*CW*WC*PNAL*PNAL))/4.D0
             CQM=1.D0 *PI*WP2/(CW*CW*WC*RR)*
     &           (-CQM1+2.D0*CQM2*MM/(AM*CW*WC*PNAL*PNAL))/2.D0
             CRM=2.D0 *PI*WP2/(CW*CW)*
     &           (-CRM1+2.D0*CRM2*MM/(AM*CW*WC*PNAL*PNAL))

            ELSEIF(MODEFA.EQ.4) THEN 

!              WRITE(6,'(A,5I5)',ADVANCE='NO') 
!              WRITE(6,'(A,5I5)')
!     &              'MD,ND,NTH,NHH,NR=',MD,ND,NTH,NHH,NR
               
               IF(RHOL.LT.1.D0) THEN
          CALL WMDPFAA(CW,RHOL,RKPR,AE2N0,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)
          WP02=AE2N0*1.D20/(AM*EPS0)
             CPM=0.5D0*PI*WP02/(CW*CW*WC*WC*RR*RR) *
     &           (CPM1+CPM2*MM/(RHOL*AM*CW*WC)) !0.5D0*CPM1
             CQM=PI*WP02/(CW*CW*WC*RR)*
     &           (CQM1+CQM2*MM/(RHOL*AM*CW*WC))  !0.5D0*CQM1
             CRM=2.D0 *PI*WP02/(CW*CW)*
     &           (CRM1+CRM2*MM/(RHOL*AM*CW*WC))  !0.5D0*CQR1
             ELSE
               CPM=(0.D0,0.D0)
               CQM=(0.D0,0.D0)
               CRM=(0.D0,0.D0)
             ENDIF

          ELSEIF(MODEFA.EQ.5) THEN
                IF((NR.EQ.2).OR.(NR.EQ.20).OR.(NR.EQ.40)) THEN
                  CALL WMDPFA(CX,CFN,COEF,RHOR,CPM,CQM,CRM,1)
                  WRITE(6,'(3I5)') NR,NN,MM
                  WRITE(6,'(1P6E12.4)') CPM,CQM,CRM
                  
                  IF(RHOL.LT.1.D0) THEN
          CALL WMDPFAA(CW,RHOL,RKPR,AE2N0,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)
                 WP02=AE2N0*1.D20/(AM*EPS0)
                 CPM=(PI/2.D0)*WP02/(CW*CW*WC*WC*RR*RR)*
     &            (CPM1+CPM2*MM/(RHOL*AM*CW*WC))  
                 CQM=PI*WP02/(CW*CW*WC*RR)*
     &            (CQM1+CQM2*MM/(RHOL*AM*CW*WC)) 
                 CRM=2.D0 *PI*WP02/(CW*CW)*
     &            (CRM1+CRM2*MM/(RHOL*AM*CW*WC))

                 WRITE(6,'(1P6E12.4)') CPM,CQM,CRM
                   ELSE
                  CPM=(0.D0,0.D0)
                  CQM=(0.D0,0.D0)
                  CRM=(0.D0,0.D0)

                   ENDIF
                ELSE
                  CPM=(0.D0,0.D0)
                  CQM=(0.D0,0.D0)
                  CRM=(0.D0,0.D0)
                ENDIF 
              ELSE
            RETURN
!              IF(RL.LT.RA) THEN
!       CALL WMDPFA3(CW,RHOL,RKPR,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)
!            CPM=5.D-1*PI*WP2/(CW*CW*WC*WC*RR*RR)*                         !FAR(NR)
!     &          (-CPM1+CPM2*MM/(AM*CW*WC*RL))/4.D0
!            CQM=1.D0 *PI*WP2/(CW*CW*WC*RR)*
!     &          (-CQM1+CQM2*MM/(AM*CW*WC*RL))/2.D0
!            CRM=2.D0 *PI*WP2/(CW*CW)*                                     !(RNA/(PNA*1.D20))*
!     &          (-CRM1+CRM2*MM/(AM*CW*WC*RL))
!               ELSE
!                RETURN
!              ENDIF
             ENDIF
!             WRITE(6,'(1P5E12.4)') RHOL,RKPR,AE2N0,WP02
!             WRITE(6,'(1P6E12.4)') CPM,CQM,CRM

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
     &     =CTNSR(1,1,MD,ND,NTH,NHH)+CPM*SIN(ANGTH)**2!COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NHH)
     &     =CTNSR(1,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NHH)
     &     =CTNSR(1,3,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)!COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NHH)
     &     =CTNSR(2,1,MD,ND,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NHH)
     &     =CTNSR(2,2,MD,ND,NTH,NHH)+CPM*COS(ANGTH)**2!SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NHH)
     &     =CTNSR(2,3,MD,ND,NTH,NHH)+CQM*COS(ANGTH)!SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NHH)
     &     =CTNSR(3,1,MD,ND,NTH,NHH)+CQM*SIN(ANGTH)!COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NHH)
     &     =CTNSR(3,2,MD,ND,NTH,NHH)+CQM*COS(ANGTH)!SIN(ANGTH)
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
         DO ND=NDMIN,NDMAX
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
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
         DO ND=NDMIN,NDMAX
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MDMAX
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
