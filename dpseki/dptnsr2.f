C     $Id: dptnsr2.f,v 1.7 2013/01/20 23:24:02 fukuyama Exp $
C
C     ***************************************
C         COMPONENTS OF DIELECTRIC TENSOR
C             MAGNETIC FIELD     (0,   0,   B)
C             WAVE NUMBER VECTOR (k_x, 0, k_z)
C
C             CLDISP(1)=EPS_XX
C             CLDISP(2)=EPS_ZZ - EPS_XX
C             CLDISP(3)=EPS_YY - EPS_XX
C             CLDISP(4)=EPS_ZX
C             CLDISP(5)=EPS_XY
C             CLDISP(6)=EPS_YZ

C     ****** CALCULATE HOT DIELECTRIC TENSOR ******
C
      SUBROUTINE DPTNFK2(CW,CKPR,CKPP,NS,CLDISP)
C
      USE libdsp,ONLY: dspfn
      USE plcomm
      USE pllocal
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CLDISP(6)
C
      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      CKPP2=CKPP**2
C
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      IF(MODELP(NS).EQ.11) THEN
         CWN=CW-RU(NS)* CKPR+CI*PZCL(NS)*WW
         CWU=CW-RU(NS)* CKPR
         CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0)
         CWC=AE*BABS/AM
         IF(NS.EQ.1) THEN
            CPERP=(0.D0,0.D0)
            CPARA=CI*CWP/(CW*CW*PZCL(NS))
         ELSE
            CPERP=   CWP*CWU*CWN/(CW*CW*CWC**2)
            CPARA=CI*CWP/(CW*CW*PZCL(NS))
         END IF
         CCROS=(0.D0,0.D0)
         CPERM=(0.D0,0.D0)
      ELSE IF(MODELP(NS).EQ.12) THEN
         CWN=CW-RU(NS)* CKPR+CI*PZCL(NS)*WW
         CWU=CW-RU(NS)* CKPR
         CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0)
         CWC=AE*BABS/AM
         CPERP=-   CWP*CWU*CWN/(CW*CW*(CWN**2-CWC**2))
         CCROS= CI*CWP*CWC*CWU/(CW*CW*(CWN**2-CWC**2))
         CPARA=-   CWP/(CW*CW)*(CW*CW/(CWU*CWN)
     &         +2*CKPP2*RU(NS)*RU(NS)/(CWN*CWN-CWC*CWC)*CWN/CWU)
         CPERM=(0.D0,0.D0)
      ELSE IF(MODELP(NS).EQ.13) THEN
         RT=(RTPR(NS)+2*RTPP(NS))/3.D0
         RKTPR=ABS(CKPR)*SQRT(2.D0*RT*AEE*1.D3/AM)
         CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0*CW*CW)
         CWC=AE*BABS/AM
         CGZ0=(CW       -CKPR*RU(NS))/RKTPR
         CPERP=0.D0
         CCROS=0.D0
         CPARA=0.D0
         DO NC=-1,1
            CGZ=(CW-NC*CWC-CKPR*RU(NS))/RKTPR
            CALL DSPFN(CGZ,CZ,CDZ)
            IF(NC.EQ.-1) THEN
               CPERP=CPERP+   CWP*CGZ0*CZ/2
               CCROS=CCROS-CI*CWP*CGZ0*CZ/2
            ELSEIF(NC.EQ.0) THEN
               CADD=1+CKPR*RU(NS)*CW/(CW-CKPR*RU(NS))**2 !
               CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD
            ELSEIF(NC.EQ.1) THEN
               CPERP=CPERP+   CWP*CGZ0*CZ/2
               CCROS=CCROS+CI*CWP*CGZ0*CZ/2
            ENDIF
         ENDDO
         CPERM=(0.D0,0.D0)
      ELSE IF(MODELP(NS).EQ.14.OR.MODELP(NS).EQ.15) THEN
         CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0*CW*CW)
         RWC=AE*BABS/AM
         RKTPR=ABS(CKPR)*SQRT(2.D0*RTPR(NS)*AEE*1.D3/AM)
         RKTPP=RTPP(NS)*AEE*1.D3/(AM*RWC*RWC)
         CWC=DCMPLX(RWC,0.D0)
         CGZ0=(CW       -CKPR*RU(NS))/RKTPR
          
         CPERP1= 0.D0
         CPERP2= 0.D0
         CPERM2= 0.D0
         CPARA1= 0.D0
         CPARA2= 0.D0
         CCROS1= 0.D0
         CCROS2= 0.D0
         DO NC=-2,2
            CGZ =(CW-NC*CWC-CKPR*RU(NS))/RKTPR
            CADD=1+CKPR*RU(NS)/(CW-CKPR*RU(NS)-NC*CWC)
            CALL DSPFN(CGZ,CZ,CDZ)
            IF(NC.EQ.-2) THEN 
               CPERP2=CPERP2+CWP*RKTPP*CGZ0*CZ
               CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*CZ
            ELSEIF(NC.EQ.-1) THEN
               CPERP1=CPERP1+CWP*CGZ0*CZ/2
               CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
               CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
               CPARA2=CPARA2
     &               -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
               CCROS1=CCROS1-CI*CWP*CGZ0*CZ/2
               CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*2*CZ
            ELSEIF(NC.EQ.0) THEN
               CPERM2=CPERM2+CWP*RKTPP*CGZ0*4*CZ
               CPARA1=CPARA1-CWP*CDZ*CGZ*CGZ0*CADD*CADD
               CPARA2=CPARA2
     &               +CWP*RKTPP*CGZ0*2*CGZ*CDZ*CADD*CADD
            ELSEIF(NC.EQ.1) THEN
               CPERP1=CPERP1+CWP*CGZ0*CZ/2
               CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
               CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
               CPARA2=CPARA2
     &               -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
               CCROS1=CCROS1+CI*CWP*CGZ0*CZ/2
               CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*2*CZ
            ELSEIF(NC.EQ.2) THEN
               CPERP2=CPERM2+CWP*RKTPP*CGZ0*CZ
               CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*CZ
            ENDIF
         ENDDO
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
         CPERP=CPERP1+0.5D0*CPERP2*CKPP2
         CPERM=       0.5D0*CPERM2*CKPP2
         CPARA=CPARA1+0.5D0*CPARA2*CKPP2
         CCROS=CCROS1+0.5D0*CCROS2*CKPP2
      ELSE
         CPERP=(0.D0,0.D0)
         CPARA=(0.D0,0.D0)
         CCROS=(0.D0,0.D0)
         CPERM=(0.D0,0.D0)
      ENDIF
C
      CLDISP(1)=CLDISP(1)+CPERP
      CLDISP(2)=CLDISP(2)+CPARA-CPERP
      CLDISP(3)=CLDISP(3)+CPERM
      CLDISP(5)=CLDISP(5)-CCROS
      RETURN
      END
C
C     ****** CALCULATE DRIFT KINETIC DIELECTRIC TENSOR ******
C
      SUBROUTINE DPTNDK0(CW,CKPR,CKPP,NS,CLDISP)
C
      USE plcomm
      USE pllocal
      INCLUDE '../dp/dpcomm.inc'
      DIMENSION CLDISP(6)
C
      RLA=RLN(NS) ! inverse of density scale lenth

      RNA=RN(NS)*1.D20
      RTA=RTPP(NS)*AEE*1.D3
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)
C
      WC=AE*BABS/AM
      RHOA=VTA/WC
      RHOR=RHOA/RR
      COEF=CI*WP2*SQRT(PI)/(CW*CW)
      IF(NR.EQ.1) THEN
         CWAST=0.D0
      ELSE
         CWAST=RTA*RLA/(BABS*CW*AE)
      ENDIF
C
      CFN=1.D0-CKPP*CWAST
      IF(ABS(CKPR).LT.1.D-5) CKPR=1.D-5
C
      CX=CW/(ABS(CKPR)*VTA)
      IF(ABS(CX).GT.5.D0) THEN
         CEX=(0.D0,0.D0)
      ELSE
         CEX=CX*EXP(-CX*CX)
      ENDIF
      CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
      CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
      CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0
C
      CLDISP(1)=CLDISP(1)+CPM
      CLDISP(2)=CLDISP(2)+CRM-CPM
      CLDISP(4)=CLDISP(4)+CQM
      RETURN
      END
