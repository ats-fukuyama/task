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
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         DO ND=-NDSIZX,NDSIZX
         DO MD=-MDSIZX,MDSIZX
            CTNSR(1,1,MD,ND,NTH,NPH)=0.D0
            CTNSR(1,2,MD,ND,NTH,NPH)=0.D0
            CTNSR(1,3,MD,ND,NTH,NPH)=0.D0
            CTNSR(2,1,MD,ND,NTH,NPH)=0.D0
            CTNSR(2,2,MD,ND,NTH,NPH)=0.D0
            CTNSR(2,3,MD,ND,NTH,NPH)=0.D0
            CTNSR(3,1,MD,ND,NTH,NPH)=0.D0
            CTNSR(3,2,MD,ND,NTH,NPH)=0.D0
            CTNSR(3,3,MD,ND,NTH,NPH)=0.D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      IF((MOD(MODELA,2).EQ.1).AND.(NS.EQ.3)) THEN
         CALL WMTNAX(NR)
      ELSEIF((MOD(MODELA/2,2).EQ.1).AND.(NS.EQ.1)) THEN
         CALL WMTNEX(NR)
      ELSE
         CALL WMTNSX(NR,NS)
      ENDIF
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
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION CGZ(1),CZ(1),CDZ(1)
C
      CW=2*PI*CRF*1.D6
      WW=DBLE(CW)
C
      CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
C
C      IF(NR.EQ.1) THEN
C         WRITE(6,*) 'RN  :',RN(1),RN(2),RN(3)
C         WRITE(6,*) 'RTPR:',RTPR(1),RTPR(2),RTPR(3)
C         WRITE(6,*) 'RTPP:',RTPP(1),RTPP(2),RTPP(3)
C         WRITE(6,*) 'RU  :',RU(1)  ,RU(2)  ,RU(3)
C      ENDIF
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
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
C            WRITE(11,'(5I5,1P2E12.4)') NR,NN,MM,NPH,NTH,RKPR,BABS
C
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            RNPR=VC*RKPR/WW
C
            IF(MODELP.EQ.4) THEN
               DTT=1.D0
               DTX=0.D0
               DO NSS=1,NSMAX
                  AM=PA(NSS)*AMP
                  AE=PZ(NSS)*AEE
                  WP=AE*AE*RN(NSS)/(AM*EPS0)
                  WC=AE*BABS/AM
                  DTT=DTT-WP/(WW*WW-WC*WC)
                  DTX=DTX-WP*WC/((WW*WW-WC*WC)*WW)
               ENDDO
               RNPP2=((DTT-RNPR**2)**2-DTX**2)/(DTT-RNPR**2)
               RKPP2=RNPP2*WW*WW/(VC*VC)
               RKX2=     MM*MM*RG22(NTH,NPH,NR)*XRHO(NR)**2
     &             +2.D0*MM*NN*RG23(NTH,NPH,NR)*XRHO(NR)
     &             +     NN*NN*RG33(NTH,NPH,NR)
               IF(RKPP2.GT.0.D0) THEN
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
                  RKT2=0.D0
                  RKR2=0.D0
                  UXX2=0.D0
                  UYY2=0.D0
               ENDIF
            ELSE
               RKPP2=0.D0
               RKT2=0.D0
               RKR2=0.D0
               UXX2=0.D0
               UYY2=0.D0
            ENDIF
C
            AM=PA(NS)*AMP
            AE=PZ(NS)*AEE
            IF(MODELP.EQ.1) THEN
               CWN=CW-RU(NS)* RKPR+CI*PZCL(NS)*WW
               CWU=CW-RU(NS)* RKPR
               CWP=AE*AE*RN(NS)/(AM*EPS0)
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
            ELSE IF(MODELP.EQ.2) THEN
               CWN=CW-RU(NS)* RKPR+CI*PZCL(NS)*WW
               CWU=CW-RU(NS)* RKPR
               CWP=AE*AE*RN(NS)/(AM*EPS0)
               CWC=AE*BABS/AM
               CPERP=-   CWP*CWU*CWN/(CW*CW*(CWN**2-CWC**2))
               CCROS= CI*CWP*CWC*CWU/(CW*CW*(CWN**2-CWC**2))
               CPARA=-   CWP/(CW*CW)
     &                   *(CW*CW/(CWU*CWN)
     &                    +2*RKPP2*RU(NS)*RU(NS)/(CWN*CWN-CWC*CWC)
     &                     *CWN/CWU)
               CPERM=(0.D0,0.D0)
C               WRITE(6,*) NR,NS,CPERP,CPARA
            ELSE IF(MODELP.EQ.3) THEN
               RT=(RTPR(NS)+2*RTPP(NS))/3.D0
               RKTPR=ABS(RKPR)*SQRT(2.D0*RT/AM)
               CWP=AE*AE*RN(NS)/(AM*EPS0*CW*CW)
               CWC=AE*BABS/AM
               CGZ0=(CW       -RKPR*RU(NS))/RKTPR
               CPERP=0.D0
               CCROS=0.D0
               CPARA=0.D0
               DO NC=-1,1
                  CGZ(1)=(CW-NC*CWC-RKPR*RU(NS))/RKTPR
                  CALL DSPFNA(1,CGZ,CZ,CDZ)
                  IF(NC.EQ.-1) THEN
                     CPERP=CPERP+   CWP*CGZ0*CZ(1)/2
                     CCROS=CCROS-CI*CWP*CGZ0*CZ(1)/2
                  ELSEIF(NC.EQ.0) THEN
                     CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS))
                     CPARA=CPARA-CWP*CDZ(1)*CGZ(1)*CGZ0*CADD*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP=CPERP+   CWP*CGZ0*CZ(1)/2
                     CCROS=CCROS+CI*CWP*CGZ0*CZ(1)/2
                  ENDIF
               ENDDO
               CPERM=(0.D0,0.D0)
            ELSE IF(MODELP.EQ.4) THEN
               RT=(RTPR(NS)+2*RTPP(NS))/3.D0
               RKTPR=ABS(RKPR)*SQRT(2.D0*RT/AM)
               CWP=AE*AE*RN(NS)/(AM*EPS0*CW*CW)
               RWC=AE*BABS/AM
               CWC=DCMPLX(RWC,0.D0)
               RKTPP=RT/(AM*RWC*RWC)
               CGZ0=(CW       -RKPR*RU(NS))/RKTPR
               CPERP1= 0.D0
               CPERP2= 0.D0
               CPERM2= 0.D0
               CPARA1= 0.D0
               CPARA2= 0.D0
               CCROS1= 0.D0
               CCROS2= 0.D0
               DO NC=-2,2
                  CGZ(1) =(CW-NC*CWC-RKPR*RU(NS))/RKTPR
                  CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS)-NC*CWC)
                  CALL DSPFNA(1,CGZ,CZ,CDZ)
                  IF(NC.EQ.-2) THEN 
                     CPERP2=CPERP2+CWP*RKTPP*CGZ0*CZ(1)
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*CZ(1)
                  ELSEIF(NC.EQ.-1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ(1)/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ(1)
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ(1)
                     CPARA2=CPARA2
     &                     -CWP*RKTPP*CGZ0*CGZ(1)*CDZ(1)*CADD*CADD
                     CCROS1=CCROS1-CI*CWP*CGZ0*CZ(1)/2
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*2*CZ(1)
                  ELSEIF(NC.EQ.0) THEN
                     CPERM2=CPERM2+CWP*RKTPP*CGZ0*4*CZ(1)
                     CPARA1=CPARA1-CWP*CDZ(1)*CGZ(1)*CGZ0*CADD*CADD
                     CPARA2=CPARA2
     &                     +CWP*RKTPP*CGZ0*2*CGZ(1)*CDZ(1)*CADD*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ(1)/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ(1)
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ(1)
                     CPARA2=CPARA2
     &                     -CWP*RKTPP*CGZ0*CGZ(1)*CDZ(1)*CADD*CADD
                     CCROS1=CCROS1+CI*CWP*CGZ0*CZ(1)/2
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*2*CZ(1)
                  ELSEIF(NC.EQ.2) THEN
                     CPERP2=CPERM2+CWP*RKTPP*CGZ0*CZ(1)
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*CZ(1)
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
C            CTNSR(1,1,MD,ND,NTH,NPH)= CPERP + UXX2*CPERM
C            CTNSR(1,2,MD,ND,NTH,NPH)=-CCROS
C            CTNSR(1,3,MD,ND,NTH,NPH)= 0.D0
C            CTNSR(2,1,MD,ND,NTH,NPH)= CCROS
C            CTNSR(2,2,MD,ND,NTH,NPH)= CPERP + UYY2*CPERM
C            CTNSR(2,3,MD,ND,NTH,NPH)= 0.D0
C            CTNSR(3,1,MD,ND,NTH,NPH)= 0.D0
C            CTNSR(3,2,MD,ND,NTH,NPH)= 0.D0
C            CTNSR(3,3,MD,ND,NTH,NPH)= CPARA
C
            CTNSR(1,1,MD,ND,NTH,NPH)
     &     =CTNSR(1,1,MD,ND,NTH,NPH) + CPERP + UXX2*CPERM
            CTNSR(1,2,MD,ND,NTH,NPH)
     &     =CTNSR(1,2,MD,ND,NTH,NPH) - CCROS
            CTNSR(2,1,MD,ND,NTH,NPH)
     &     =CTNSR(2,1,MD,ND,NTH,NPH) + CCROS
            CTNSR(2,2,MD,ND,NTH,NPH)
     &     =CTNSR(2,2,MD,ND,NTH,NPH) + CPERP + UYY2*CPERM
            CTNSR(3,3,MD,ND,NTH,NPH)
     &     =CTNSR(3,3,MD,ND,NTH,NPH) + CPARA
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
C     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******
C
      SUBROUTINE WMTNAX(NR)
C
      INCLUDE 'wmcomm.inc'
C
      CW=2*PI*CRF*1.D6
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
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
C         ANGPH=(NPH-1)*2.D0*PI/NPHMAX
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
C            CTNSR(1,1,MD,ND,NTH,NPH)=CPM*COS(ANGTH)**2
C            CTNSR(1,2,MD,ND,NTH,NPH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(1,3,MD,ND,NTH,NPH)=CQM*COS(ANGTH)
C            CTNSR(2,1,MD,ND,NTH,NPH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(2,2,MD,ND,NTH,NPH)=CPM*SIN(ANGTH)**2
C            CTNSR(2,3,MD,ND,NTH,NPH)=CQM*SIN(ANGTH)
C            CTNSR(3,1,MD,ND,NTH,NPH)=CQM*COS(ANGTH)
C            CTNSR(3,2,MD,ND,NTH,NPH)=CQM*SIN(ANGTH)
C            CTNSR(3,3,MD,ND,NTH,NPH)=CRM
C
            CTNSR(1,1,MD,ND,NTH,NPH)
     &     =CTNSR(1,1,MD,ND,NTH,NPH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NPH)
     &     =CTNSR(1,2,MD,ND,NTH,NPH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NPH)
     &     =CTNSR(1,3,MD,ND,NTH,NPH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NPH)
     &     =CTNSR(2,1,MD,ND,NTH,NPH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NPH)
     &     =CTNSR(2,2,MD,ND,NTH,NPH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NPH)
     &     =CTNSR(2,3,MD,ND,NTH,NPH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NPH)
     &     =CTNSR(3,1,MD,ND,NTH,NPH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NPH)
     &     =CTNSR(3,2,MD,ND,NTH,NPH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NPH)
     &     =CTNSR(3,3,MD,ND,NTH,NPH)+CRM
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
      CW=2*PI*CRF*1.D6
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
      WP2=AE*AE*RNE/(AM*EPS0)
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
         ANGTH=(NTH-1)*2.D0*PI/NTHMAX
         ANGPH=(NPH-1)*2.D0*PI/NPHMAX
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
C            CTNSR(1,1,MD,ND,NTH,NPH)=CPM*COS(ANGTH)**2
C            CTNSR(1,2,MD,ND,NTH,NPH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(1,3,MD,ND,NTH,NPH)=CQM*COS(ANGTH)
C            CTNSR(2,1,MD,ND,NTH,NPH)=CPM*COS(ANGTH)*SIN(ANGTH)
C            CTNSR(2,2,MD,ND,NTH,NPH)=CPM*SIN(ANGTH)**2
C            CTNSR(2,3,MD,ND,NTH,NPH)=CQM*SIN(ANGTH)
C            CTNSR(3,1,MD,ND,NTH,NPH)=CQM*COS(ANGTH)
C            CTNSR(3,2,MD,ND,NTH,NPH)=CQM*SIN(ANGTH)
C            CTNSR(3,3,MD,ND,NTH,NPH)=CRM
C
            CTNSR(1,1,MD,ND,NTH,NPH)
     &     =CTNSR(1,1,MD,ND,NTH,NPH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,MD,ND,NTH,NPH)
     &     =CTNSR(1,2,MD,ND,NTH,NPH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,MD,ND,NTH,NPH)
     &     =CTNSR(1,3,MD,ND,NTH,NPH)+CQM*COS(ANGTH)
            CTNSR(2,1,MD,ND,NTH,NPH)
     &     =CTNSR(2,1,MD,ND,NTH,NPH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,MD,ND,NTH,NPH)
     &     =CTNSR(2,2,MD,ND,NTH,NPH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,MD,ND,NTH,NPH)
     &     =CTNSR(2,3,MD,ND,NTH,NPH)+CQM*SIN(ANGTH)
            CTNSR(3,1,MD,ND,NTH,NPH)
     &     =CTNSR(3,1,MD,ND,NTH,NPH)+CQM*COS(ANGTH)
            CTNSR(3,2,MD,ND,NTH,NPH)
     &     =CTNSR(3,2,MD,ND,NTH,NPH)+CQM*SIN(ANGTH)
            CTNSR(3,3,MD,ND,NTH,NPH)
     &     =CTNSR(3,3,MD,ND,NTH,NPH)+CRM
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      RETURN
      END
