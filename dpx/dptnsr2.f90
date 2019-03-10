MODULE dptnsr2
!
!     ***************************************
!         COMPONENTS OF DIELECTRIC TENSOR
!             MAGNETIC FIELD     (0,   0,   B)
!             WAVE NUMBER VECTOR (k_x, 0, k_z)
!
!             CLDISP(1)=EPS_XX
!             CLDISP(2)=EPS_ZZ - EPS_XX
!             CLDISP(3)=EPS_YY - EPS_XX
!             CLDISP(4)=EPS_ZX
!             CLDISP(5)=EPS_XY
!             CLDISP(6)=EPS_YZ

  PRIVATE
  PUBLIC DPTNFK2,DPTNDK0

CONTAINS

!     ****** CALCULATE HOT DIELECTRIC TENSOR ******

  SUBROUTINE DPTNFK2(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE dpcomm
      USE plprof
      USE libdsp,ONLY: dspfn
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CKPP2,CWN,CWU,CPERP,CPARA,CCROS,CPERM,CWNU
      COMPLEX(rkind):: CGZ0,CGZ,CZ,CDZ,CADD
      COMPLEX(rkind):: CPERP1,CCROS1,CPARA1
      COMPLEX(rkind):: CPERP2,CCROS2,CPERM2,CPARA2
      REAL(rkind):: AM,AE,RT,RKTPR,RKTPP,RWC
      INTEGER:: I,NC

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF

      CKPP2=CKPP**2

      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      IF(MODELP(NS).EQ.11) THEN
         CWN=CW-plf(NS)%RU*CKPR+CI*plf(NS)%RZCL*CW
         CWU=CW-plf(NS)%RU*CKPR
         CWP=AE*AE*plf(NS)%RN*1.D20/(AM*EPS0)
         CWC=AE*mag%BABS/AM
         IF(NS.EQ.1) THEN
            CPERP=(0.D0,0.D0)
            CPARA=CI*CWP/(CW*CW*plf(NS)%RZCL)
         ELSE
            CPERP=   CWP*CWU*CWN/(CW*CW*CWC**2)
            CPARA=CI*CWP/(CW*CW*plf(NS)%RZCL)
         END IF
         CCROS=(0.D0,0.D0)
         CPERM=(0.D0,0.D0)
      ELSE IF(MODELP(NS).EQ.12) THEN
         CWN=CW-plf(ns)%RU*CKPR+CI*plf(NS)%RZCL*CW
         CWU=CW-plf(NS)%RU*CKPR
         CWP=AE*AE*plf(NS)%RN*1.D20/(AM*EPS0)
         CWC=AE*mag%BABS/AM
         CPERP=-   CWP*CWU*CWN/(CW*CW*(CWN**2-CWC**2))
         CCROS= CI*CWP*CWC*CWU/(CW*CW*(CWN**2-CWC**2))
         CPARA=-   CWP/(CW*CW)*(CW*CW/(CWU*CWN) &
               +2*CKPP2*plf(NS)%RU*plf(NS)%RU/(CWN*CWN-CWC*CWC)*CWN/CWU)
         CPERM=(0.D0,0.D0)
      ELSE IF(MODELP(NS).EQ.13) THEN
         RT=(plf(NS)%RTPR+2*plf(NS)%RTPP)/3.D0
         RKTPR=ABS(CKPR)*SQRT(2.D0*RT*AEE*1.D3/AM)
         CWP=AE*AE*plf(NS)%RN*1.D20/(AM*EPS0*CW*CW)
         CWC=AE*mag%BABS/AM
         CGZ0=(CW       -CKPR*plf(NS)%RU)/RKTPR
         CPERP=0.D0
         CCROS=0.D0
         CPARA=0.D0
         DO NC=-1,1
            CGZ=(CW*CWNU-NC*CWC-CKPR*plf(NS)%RU)/RKTPR
            CALL DSPFN(CGZ,CZ,CDZ)
            IF(NC.EQ.-1) THEN
               CPERP=CPERP+   CWP*CGZ0*CZ/2
               CCROS=CCROS-CI*CWP*CGZ0*CZ/2
            ELSEIF(NC.EQ.0) THEN
               CADD=1+CKPR*plf(NS)%RU*CW/(CW-CKPR*plf(NS)%RU)**2 !
               CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD
            ELSEIF(NC.EQ.1) THEN
               CPERP=CPERP+   CWP*CGZ0*CZ/2
               CCROS=CCROS+CI*CWP*CGZ0*CZ/2
            ENDIF
         ENDDO
         CPERM=(0.D0,0.D0)
      ELSE IF(MODELP(NS).EQ.14.OR.MODELP(NS).EQ.15) THEN
         CWP=AE*AE*plf(NS)%RN*1.D20/(AM*EPS0*CW*CW)
         RWC=AE*mag%BABS/AM
         RKTPR=ABS(CKPR)*SQRT(2.D0*plf(NS)%RTPR*AEE*1.D3/AM)
         RKTPP=plf(NS)%RTPP*AEE*1.D3/(AM*RWC*RWC)
         CWC=DCMPLX(RWC,0.D0)
         CGZ0=(CW*CWNU  -CKPR*plf(NS)%RU)/RKTPR
          
         CPERP1= 0.D0
         CPERP2= 0.D0
         CPERM2= 0.D0
         CPARA1= 0.D0
         CPARA2= 0.D0
         CCROS1= 0.D0
         CCROS2= 0.D0
         DO NC=-2,2
            CGZ =(CW*CWNU-NC*CWC-CKPR*plf(NS)%RU)/RKTPR
            CADD=1+CKPR*plf(NS)%RU/(CW-CKPR*plf(NS)%RU-NC*CWC)
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

      CLDISP(1)=CLDISP(1)+CPERP
      CLDISP(2)=CLDISP(2)+CPARA-CPERP
      CLDISP(3)=CLDISP(3)+CPERM
      CLDISP(5)=CLDISP(5)-CCROS
      RETURN
  END SUBROUTINE DPTNFK2

!     ****** CALCULATE DRIFT KINETIC DIELECTRIC TENSOR ******

  SUBROUTINE DPTNDK0(CW,CKPR,CKPP,NS,mag,plf,grd,CLDISP)

      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      TYPE(pl_grd_type),DIMENSION(nsmax),INTENT(IN):: grd
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: COEF,CWAST,CFN,CX,CEX,CPM,CQM,CRM,CKPR1
      REAL(rkind):: RLA,RNA,RTA,AM,AE,VTA,WP2,WC,RHOA,RHOR

      RLA=grd(NS)%grdn ! inverse of density scale lenth

      RNA=plf(NS)%RN*1.D20
      RTA=plf(NS)%RTPP*AEE*1.D3
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)

      WC=AE*mag%BABS/AM
      RHOA=VTA/WC
      RHOR=RHOA/RR
      COEF=CI*WP2*SQRT(PI)/(CW*CW)
      CWAST=RTA*RLA/(mag%BABS*CW*AE)

      CFN=1.D0-CKPP*CWAST
      CKPR1=CKPR
      IF(ABS(CKPR1).LT.1.D-5) CKPR1=1.D-5

      CX=CW/(ABS(CKPR1)*VTA)
      IF(ABS(CX).GT.5.D0) THEN
         CEX=(0.D0,0.D0)
      ELSE
         CEX=CX*EXP(-CX*CX)
      ENDIF
      CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
      CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
      CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

      CLDISP(1)=CLDISP(1)+CPM
      CLDISP(2)=CLDISP(2)+CRM-CPM
      CLDISP(4)=CLDISP(4)+CQM
      RETURN
  END SUBROUTINE DPTNDK0
END MODULE dptnsr2
