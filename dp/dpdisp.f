C     $Id$
C
C     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******
C
      COMPLEX*16 FUNCTION CFDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)
C
      INCLUDE 'dpcomm.h'
      DIMENSION CDET(3,3)
C
      CALL DPDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)
C
      CDET11=CDET(1,1)
      CDET12=CDET(1,2)
      CDET13=CDET(1,3)
      CDET21=CDET(2,1)
      CDET22=CDET(2,2)
      CDET23=CDET(2,3)
      CDET31=CDET(3,1)
      CDET32=CDET(3,2)
      CDET33=CDET(3,3)
C
      CFDISP =CDET33*(CDET11*CDET22-CDET12*CDET21)
     &       +CDET12*CDET23*CDET31-CDET23*CDET32*CDET11
     &       +CDET13*CDET32*CDET21-CDET13*CDET31*CDET22
C
      RETURN
      END
C
C     ****** CALCULATE REAL PART OF DETERMINANT OF DISPERSION TENSOR ******
C
      COMPLEX*16 FUNCTION CFDISPR(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS)
C
      INCLUDE 'dpcomm.h'
      DIMENSION CDET(3,3)
C
      CALL DPDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)
C
      CDET11=0.5D0*(CDET(1,1)+DCONJG(CDET(1,1)))
      CDET12=0.5D0*(CDET(1,2)+DCONJG(CDET(2,1)))
      CDET13=0.5D0*(CDET(1,3)+DCONJG(CDET(3,1)))
      CDET21=0.5D0*(CDET(2,1)+DCONJG(CDET(1,2)))
      CDET22=0.5D0*(CDET(2,2)+DCONJG(CDET(2,2)))
      CDET23=0.5D0*(CDET(2,3)+DCONJG(CDET(3,2)))
      CDET31=0.5D0*(CDET(3,1)+DCONJG(CDET(1,3)))
      CDET32=0.5D0*(CDET(3,2)+DCONJG(CDET(2,3)))
      CDET33=0.5D0*(CDET(3,3)+DCONJG(CDET(3,3)))
C
      CFDISPR=CDET33*(CDET11*CDET22-CDET12*CDET21)
     &       +CDET12*CDET23*CDET31-CDET23*CDET32*CDET11
     &       +CDET13*CDET32*CDET21-CDET13*CDET31*CDET22
C
      RETURN
      END
C
C     ****** CALCULATE DISPERSION TENSOR ******
C
      SUBROUTINE DPDISP(CRF,CKX,CKY,CKZ,XPOS,YPOS,ZPOS,CDET)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CDISP(6),CLDISP(6),CDET(3,3)
C
      CALL PLMAG(XPOS,YPOS,ZPOS,PSIN)
      CALL PLPROF(PSIN)
C
      CW=2.D0*PI*1.D6*CRF
      CKPR=BNX*CKX+BNY*CKY+BNZ*CKZ
      IF(ABS(CKPR).LE.1.D-8) CKPR=1.D-8
      CKPP=SQRT(CKX**2+CKY**2+CKZ**2-CKPR**2)
      IF(ABS(CKPP).LE.1.D-8) CKPP=1.D-8
C
      CL2=VC*VC/(CW*CW)
      CDV11 =1.D0-CL2*(CKY*CKY+CKZ*CKZ)
      CDV12 =     CL2* CKX*CKY
      CDV13 =     CL2* CKX*CKZ
      CDV21 =     CL2* CKY*CKX
      CDV22 =1.D0-CL2*(CKZ*CKZ+CKX*CKX)
      CDV23 =     CL2* CKY*CKZ
      CDV31 =     CL2* CKZ*CKX
      CDV32 =     CL2* CKZ*CKY
      CDV33 =1.D0-CL2*(CKX*CKX+CKY*CKY)
C
      DO I=1,6
         CDISP(I)=0.D0
      ENDDO
C
      DO NS=1,NSMAX
         CALL DPTENS(NS,CLDISP)
         DO I=1,6
            CDISP(I)=CDISP(I)+CLDISP(I)
         ENDDO
      ENDDO
C
      CU11=( (1.D0-BNX**2)*CKX-BNX*BNY*CKY-BNX*BNZ*CKZ)/CKPP
      CU12=(-BNY*BNX*CKX+(1.D0-BNY**2)*CKY-BNY*BNZ*CKZ)/CKPP
      CU13=(-BNZ*BNX*CKX-BNZ*BNY*CKY+(1.D0-BNZ**2)*CKZ)/CKPP
      CU21=(BNY*CKZ-BNZ*CKY)/CKPP
      CU22=(BNZ*CKX-BNX*CKZ)/CKPP
      CU23=(BNX*CKY-BNY*CKX)/CKPP
      CU31=BNX
      CU32=BNY
      CU33=BNZ
C
      CDET11=CDV11+CDISP(1)
     &            +     CU31*CU31*CDISP(2)
     &            +     CU21*CU21*CDISP(3)
     &            +2.D0*CU31*CU11*CDISP(4)
      CDET22=CDV22+CDISP(1)
     &            +     CU32*CU32*CDISP(2)
     &            +     CU22*CU22*CDISP(3)
     &            +2.D0*CU32*CU12*CDISP(4)
      CDET33=CDV33+CDISP(1)
     &            +     CU33*CU33*CDISP(2)
     &            +     CU23*CU23*CDISP(3)
     &            +2.D0*CU33*CU13*CDISP(4)
      CDET1P=           CU31*CU32*CDISP(2)
     &            +     CU21*CU22*CDISP(3)
     &            +     CU31*CU12*CDISP(4)
     &            +     CU32*CU11*CDISP(4)
      CDET2P=           CU32*CU33*CDISP(2)
     &            +     CU22*CU23*CDISP(3)
     &            +     CU32*CU13*CDISP(4)
     &            +     CU33*CU12*CDISP(4)
      CDET3P=           CU33*CU31*CDISP(2)
     &            +     CU23*CU21*CDISP(3)
     &            +     CU33*CU11*CDISP(4)
     &            +     CU31*CU13*CDISP(4)
      CDET1M=           CU33*CDISP(5)
     &            +     CU13*CDISP(6)
      CDET2M=           CU31*CDISP(5)
     &            +     CU11*CDISP(6)
      CDET3M=           CU32*CDISP(5)
     &            +     CU12*CDISP(6)
      CDET(1,1)=CDET11
      CDET(2,2)=CDET22
      CDET(3,3)=CDET33
      CDET(1,2)=CDV12+CDET1P+CDET1M
      CDET(2,1)=CDV21+CDET1P-CDET1M
      CDET(2,3)=CDV23+CDET2P+CDET2M
      CDET(3,2)=CDV32+CDET2P-CDET2M
      CDET(3,1)=CDV31+CDET3P+CDET3M
      CDET(1,3)=CDV13+CDET3P-CDET3M
C
      RETURN
      END
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE DPTENS(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
C
      DIMENSION CLDISP(6),CLDISP1(6),CLDISP2(6)
      DATA CI/(0.D0,1.D0)/
C
      IF(RN(NS).LE.0.D0) THEN
         CALL DPCOLD0(NS,CLDISP)
      ELSE
         ID1=MOD(MODELP(NS),10)
         ID2=MODELP(NS)/10
C
C        0-- 9 : BOTH HERMITE AND ANTI-HERMITE : ID1 MODEL
C       10--19 : BOTH HERMITE AND ANTI-HERMITE : ID1 MODEL
C       30--39 : BOTH HERMITE AND ANTI-HERMITE : ID1 MODEL
C       50--59 : BOTH HERMITE AND ANTI-HERMITE : ID1 MODEL
C
         IF(ID2.EQ.0.OR.ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5) THEN
            IF(ID1.EQ.0) THEN
               CALL DPCOLD(NS,CLDISP)
            ELSE IF(ID1.EQ.1) THEN
               CALL DPWARM(NS,CLDISP)
            ELSE IF(ID1.EQ.2) THEN
               CALL DPHOT1(NS,CLDISP)
            ELSE IF(ID1.EQ.3) THEN
               CALL DPHOT2(NS,CLDISP)
            ELSE IF(ID1.EQ.4) THEN
               CALL DPHOTN(NS,CLDISP)
            ELSE IF(ID1.EQ.5) THEN
               CALL DPHOTNR(NS,CLDISP)
            ELSE IF(ID1.EQ.6) THEN
               CALL DPFMFL(NS)
               CALL DPHOTF(NS,CLDISP)
            ELSE IF(ID1.EQ.7) THEN
               CALL DPFMFL(NS)
               CALL DPHOTR(NS,CLDISP)
            ELSE IF(ID1.EQ.8) THEN
               CALL DPFPFL(NS)
               CALL DPHOTF(NS,CLDISP)
            ELSE IF(ID1.EQ.9) THEN
               CALL DPFPFL(NS)
               CALL DPHOTR(NS,CLDISP)  
            ENDIF
C
C        20--29 : HERMITE      : COLD0
C        40--49 : HERMITE      : HOTN
C        60--69 : HERMITE      : HOTNR
C                 ANTI-HERMITE : ID1 MODEL
C
         ELSE IF(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6) THEN
            IF(ID2.EQ.2) THEN
               CALL DPCOLD0(NS,CLDISP1)
            ELSEIF(ID2.EQ.4) THEN
               CALL DPHOTN(NS,CLDISP1)
            ELSEIF(ID2.EQ.6) THEN
               CALL DPHOTNR(NS,CLDISP1)
            ENDIF
C
            IF(ID1.EQ.0) THEN
               CALL DPCOLD(NS,CLDISP2)
            ELSE IF(ID1.EQ.1) THEN
               CALL DPWARM(NS,CLDISP2)
            ELSE IF(ID1.EQ.2) THEN
               CALL DPHOT1(NS,CLDISP2)
            ELSE IF(ID1.EQ.3) THEN
               CALL DPHOT2(NS,CLDISP2)
            ELSE IF(ID1.EQ.4) THEN
               CALL DPHOTN(NS,CLDISP2)
            ELSE IF(ID1.EQ.5) THEN
               CALL DPHOTNR(NS,CLDISP2)
            ELSE IF(ID1.EQ.6) THEN
               CALL DPFMFL(NS)
               CALL DPHOTFI(NS,CLDISP2)
            ELSE IF(ID1.EQ.7) THEN
               CALL DPFMFL(NS)
               CALL DPHOTRI(NS,CLDISP2)
            ELSE IF(ID1.EQ.8) THEN
               CALL DPFPFL(NS)
               CALL DPHOTFI(NS,CLDISP2)
            ELSE IF(ID1.EQ.9) THEN
               CALL DPFPFL(NS)
               CALL DPHOTRI(NS,CLDISP2)  
            ENDIF
            CLDISP(1)=DBLE(CLDISP1(1))+CI*DIMAG(CLDISP2(1))
            CLDISP(2)=DBLE(CLDISP1(2))+CI*DIMAG(CLDISP2(2))
            CLDISP(3)=DBLE(CLDISP1(3))+CI*DIMAG(CLDISP2(3))
            CLDISP(4)=DBLE(CLDISP1(4))+CI*DIMAG(CLDISP2(4))
            CLDISP(5)=CI*DIMAG(CLDISP1(5))+DBLE(CLDISP2(5))
            CLDISP(6)=CI*DIMAG(CLDISP1(6))+DBLE(CLDISP2(6))
         ELSE
            CALL DPCOLD0(NS,CLDISP)
         ENDIF
      ENDIF
      RETURN
      END
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
C           
C     ****** CALCULATE COLD DISPERSION ******
C
      SUBROUTINE DPCOLD0(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)=0.D0
      CLDISP(4)=0.D0
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)=0.D0
      RETURN
      END
C
C     ****** CALCULATE COLD DISPERSION ******
C
      SUBROUTINE DPCOLD(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      RNYU=PZCL(NS)
      CWP=CWP/DCMPLX(1.D0,RNYU)
      CWC=CWC/DCMPLX(1.D0,RNYU)
C
      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)=0.D0
      CLDISP(4)=0.D0
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)=0.D0
      RETURN
      END
C
C     ****** CALCULATE WARM DISPERSION ******
C
      SUBROUTINE DPWARM(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPP=RTPP(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
C
      RNYU=RZCL(NS)
      CWP=CWP/DCMPLX(1.D0,RNYU)
      CWC=CWC/DCMPLX(1.D0,RNYU)
      CTPR=1.D0/(1.D0-CKPR**2*WTPR/(CW*CW*DCMPLX(1.D0,RNYU)))
      CTPP=WTPP/((1.D0-CWC*CWC)*CW*CW*DCMPLX(1.D0,RNYU))
C
      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP*CTPR*(1.D0-CKPP**2*CTPP)
     &          -CLDISP(1)
      CLDISP(3)= CWP*CKPP**2*CTPP*CTPR
      CLDISP(4)=-CKPP*CKPR*CTPP*CTPR*WTPX
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)=(0.D0, 1.D0)*CWC*CLDISP(4)
      RETURN
      END
C
C     ****** CALCULATE UNMAGNETIZE HOT DISPERSION ******
C
      SUBROUTINE DPUMAG(NS ,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
C
      RT=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
      CK2=SQRT(CKPR**2+CKPP**2)
      VT2=RT*AEE*1.D3/(AMP*PA(NS))
      CGZ=CW/SQRT(2.D0*CK2*VT2)
      CALL DSPFNA(1,CGZ,CZ,CDZ)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CFX=CKPP/CK2
      CFZ=CKPR/CK2
      CLP=-CWP*CGZ**2*CDZ
      CLT= CWP*CGZ*CZ
C
      CLDISP(1)=  CFZ**2*CLT+CFX**2*CLP
      CLDISP(2)= (CFZ**2-CFX**2)*(CLP-CLT)
      CLDISP(3)= -CFX**2*(CLP-CLT)
      CLDISP(4)= CFX*CFZ*(CLP-CLT)
      CLDISP(5)= 0.D0
      CLDISP(6)= 0.D0
      RETURN
      END
C
C     ****** CALCULATE HOT DISPERSION UPTO 1 HARMONICS ******
C
      SUBROUTINE DPHOT1(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
      DATA CI/(0.D0,1.D0)/
C
      DO 10 I=1,6
         CLDISP(I)=0.D0
   10 CONTINUE
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPP=RTPP(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      CALAM(0)=1-CBETA
      CALAM(1)=0.5D0*CBETA
      CALAM(2)=0.D0
C
      DO 1000 IC=-1,1
C
         CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
         CGZ= (1.D0-IC*CWC)*CPR
         CALL DSPFNA(1,CGZ,CZ,CDZ)
C
         CK=CKPP/CKPR
C
         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-IC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-IC*(1.D0-1.D0/WTPX))*CGZ*CDZ
C
         CLAM=CALAM(ABS(IC))
         CLAMM=CALAM(ABS(IC-1))
         CLAMP=CALAM(ABS(IC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)
C
         CLDISP(1)=CLDISP(1)+CWP*   IC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-IC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*IC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
 1000 CONTINUE
      RETURN
      END
C
C     ****** CALCULATE HOT DISPERSION UPTO 2 HARMONICS ******
C
      SUBROUTINE DPHOT2(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
      DATA CI/(0.D0,1.D0)/
C
      DO 10 I=1,6
         CLDISP(I)=0.D0
   10 CONTINUE
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPP=RTPP(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      CALAM(0)=1.D0-      CBETA+0.750D0*CBETA*CBETA
      CALAM(1)=     0.5D0*CBETA-0.500D0*CBETA*CBETA
      CALAM(2)=                 0.125D0*CBETA*CBETA
      CALAM(3)=0.D0
C
      DO 1000 IC=-2,2
C
         CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
         CGZ= (1.D0-IC*CWC)*CPR
         CALL DSPFNA(1,CGZ,CZ,CDZ)
C
         CK=CKPP/CKPR
C
         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-IC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-IC*(1.D0-1.D0/WTPX))*CGZ*CDZ
C
         CLAM=CALAM(ABS(IC))
         CLAMM=CALAM(ABS(IC-1))
         CLAMP=CALAM(ABS(IC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)
C
         CLDISP(1)=CLDISP(1)+CWP*   IC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-IC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*IC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
 1000 CONTINUE
      RETURN
      END
C
C     ****** CALCULATE HOT DISPERSION UPTO N HARMONICS ******
C
      SUBROUTINE DPHOTN(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
      DATA CI/(0.D0,1.D0)/
C
      DO 10 I=1,6
         CLDISP(I)=0.D0
   10 CONTINUE
C
      NMIN=NDISP1(NS)
      NMAX=NDISP2(NS)
      IF(MAX(ABS(NMIN),ABS(NMAX))+1.GT.NHM) THEN
        WRITE(6,*) 'XX NDISP+1 EXCEEDS NHM: NHM = ',NHM
        NMIN=-NHM+1
        NMAX= NHM-1
      ENDIF
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPP=RTPP(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      CALL LAMBDA(MAX(ABS(NMIN),ABS(NMAX))+1,CBETA,CALAM)
C
      DO 1000 IC=NMIN,NMAX
C
         IF(ABS(CKPR).LE.0.D0) THEN
            CPR=CW/SQRT(2.D0*1.D-4**2*WTPR)
            CK=CKPP/1.D-4
         ELSE
            CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
            CK=CKPP/CKPR
        ENDIF
         CGZ= (1.D0-IC*CWC)*CPR
         CALL DSPFNA(1,CGZ,CZ,CDZ)
C
         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-IC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-IC*(1.D0-1.D0/WTPX))*CGZ*CDZ
C
         CLAM=CALAM(ABS(IC))
         CLAMM=CALAM(ABS(IC-1))
         CLAMP=CALAM(ABS(IC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)
C
         CLDISP(1)=CLDISP(1)+CWP*   IC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-IC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*IC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
 1000 CONTINUE
      RETURN
      END
C
C     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******
C     ******                    IN                      ******
C     ******     WEAKLY RELATIVISTIC THERMAL PLASMA     ******
C
      SUBROUTINE DPHOTNR(NS,CLDISP)
C
      INCLUDE 'dpcomm.h'
      INCLUDE '../pl/plcom2.h'
      DIMENSION CLDISP(6)
C      
      DATA CI/(0.D0,1.D0)/
C
      DELZ=1.D-6
C
      NMIN=NDISP1(NS)
      NMAX=NDISP2(NS)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=ABS(BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW))
C
      WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPP=RTPP(NS)*1.D3*AEE/(AMP*PA(NS))
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      IF(IDEBUG.EQ.1) THEN
      WRITE(6,*) 'XXX CWC,CKPP,CBETA=',
     &            DBLE(CWC),DBLE(CKPP),DBLE(CBETA)
      ENDIF
C
      RMU=VC**2/WTPR
      CNPP=CKPP*VC/CW
      CNPR=CKPR*VC/CW
      IF(IDEBUG.EQ.1) THEN
      WRITE(6,*) 'XXX RMU,CNPP,CNPR=',
     &            DBLE(RMU),DBLE(CNPP),DBLE(CNPR)
      ENDIF
C
      CSUM1=(0.D0,0.D0)
      CSUM2=(0.D0,0.D0)
      CSUM3=(0.D0,0.D0)
      CSUM4=(0.D0,0.D0)
      CSUM5=(0.D0,0.D0)
C
      DO IC=NMIN,NMAX
         IF(IC.NE.0)THEN
            N=ABS(IC)
            NSIG=IC/N
            CZ=IC*RMU*CWC
            DKAI=DKAIJO(N)
            CN1=N**2*CBETA**(N-1)/(2**N*DKAI)
            CN2=N   *CBETA**(N-1)/(2**N*DKAI)
            CN3=     CBETA**N    /(2**N*DKAI)
C
            CF1=CFQ(N+3.D0/2.D0,CZ     ,CNPR     ,RMU)
            CF2=CFQ(N+5.D0/2.D0,CZ+DELZ,CNPR     ,RMU)
            CF3=CFQ(N+5.D0/2.D0,CZ-DELZ,CNPR     ,RMU)
            CF4=CFQ(N+5.D0/2.D0,CZ     ,CNPR+DELZ,RMU)
            CF5=CFQ(N+5.D0/2.D0,CZ     ,CNPR-DELZ,RMU)
C
            CSUM1=CSUM1+     CN1*CF1
            CSUM2=CSUM2+NSIG*CN1*CF1
            CSUM3=CSUM3+     CN2*(CF2-CF3)/(DELZ*2.D0)
            CSUM4=CSUM4+NSIG*CN2*(CF2-CF3)/(DELZ*2.D0)
            CSUM5=CSUM5+     CN3*(CF4*(CNPR+DELZ)-CF5*(CNPR-DELZ))
     &                       /(DELZ*2.D0)
         ENDIF
         IF(IDEBUG.EQ.1) THEN
            WRITE(6,*) 'IC,CZ,DELZ=',IC,CZ,DELZ
            WRITE(6,*) 'CF2,CF3=',CF2,CF3
            WRITE(6,*) 'CSUM3,CSUM4=',CSUM3,CSUM4
         ENDIF
      ENDDO
      CF6=CFQ(  5.D0/2.D0,(0.D0,0.D0),CNPR+DELZ,RMU)
      CF7=CFQ(  5.D0/2.D0,(0.D0,0.D0),CNPR-DELZ,RMU)
      CPART=(CF6*(CNPR+DELZ)-CF7*(CNPR-DELZ))/(DELZ*2.D0)
C
      CE11=   -CWP*RMU          *CSUM1
      CE12= CI*CWP*RMU          *CSUM2
      CE22=   -CWP*RMU          *CSUM1
      CE23=-CI*CWP/CWC*RMU*CNPP*CNPR*CSUM4
      CE31=   -CWP/CWC*RMU*CNPP*CNPR*CSUM3
      CE33=   -CWP*RMU*(CPART+CSUM5)   
C
      CLDISP(1)=CE11
      CLDISP(2)=CE33-CE11
      CLDISP(3)=CE22-CE11
      CLDISP(4)=CE31
      CLDISP(5)=CE12
      CLDISP(6)=CE23
      IF(IDEBUG.EQ.1) THEN
         WRITE(6,*) 'CLDISP=',(CLDISP(I),I=1,6)
      ENDIF
C
      RETURN
      END
C
C     ****** CALCULATE ! ******
C
      FUNCTION DKAIJO(N)
C
      REAL*8 DKAIJO,D
      DIMENSION K(0:10)
      DATA K/1,1,2,6,24,120,720,5040,40320,362880,3628800/
C
      IF(N.LT.0) THEN
         WRITE(6,*) 'XX DKAIJO: WRONG ARGUMENT: ',N
      ELSEIF(N.LE.10) THEN
         DKAIJO=DBLE(K(N))
      ELSE
         D=DBLE(K(10))
         DO I=11,N
            D=D*DBLE(I)
         ENDDO
         DKAIJO=D
      ENDIF
C
      RETURN
      END
C
C     ****** CALCULATE F ******
C
      FUNCTION CFQ(Q,CZ,CNPR,RMU)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C
      CFQ0=CFQZ(Q,RMU-CZ)
      CFQ=CFQ0
     &   +RMU*CNPR**2/2.D0*(   CFQZ(Q-1,RMU-CZ)
     &                      -2*CFQ0
     &                      +  CFQZ(Q+1,RMU-CZ))
      RETURN
      END
C     
C     ****** CALCULATE SHKAROFSKY ******
C           *** WITH Z-FUNCTION ***
C
      FUNCTION CFQZ(Q,CZ)
C
      INCLUDE 'dpcomm.h'
C
      DATA CI/(0.D0,1.D0)/
C
      IF(ABS(CZ).GT.15.D0) THEN
         NUMAX=20
         CSUM=0.D0
         DO NU=0,NUMAX
            CTERM=-DGAMM(Q+NU)/(-CZ)**(NU+1)
            CSUM=CSUM+CTERM
            IF(ABS(CTERM).LE.1.D-12) GOTO 100
         ENDDO
  100    CONTINUE
         TEMP=DBLE(SQRT(CZ))
         IF(ABS(TEMP).LT.1.D-12) THEN
            CSUM=CSUM-  CI*PI*(-CZ)**(Q-1)*EXP(CZ)
         ELSEIF(TEMP.LT.0.D0) THEN
            CSUM=CSUM-2*CI*PI*(-CZ)**(Q-1)*EXP(CZ)
         ENDIF
         CFQZ=CSUM/DGAMM(Q)
      ELSE
         NUMAX=NINT(Q-3.D0/2.D0)
         CSUM=(0.D0,0.D0)
         DO NU=0,NUMAX 
            CTERM=(-CZ)**NU*DGAMM(Q-1-NU)
            CSUM=CSUM+CTERM
         ENDDO
         CGZ=CI*SQRT(CZ)
         CALL DSPFNA(1,CGZ,CZ2,CDZ2)
         CSUM=CSUM+SQRT(PI)*(-CZ)**NUMAX*(CI*SQRT(CZ)*CZ2)
         CFQZ=CSUM/DGAMM(Q)
      ENDIF
      RETURN
      END
C     
C     ****** CALCULATE SHKAROFSKY ******
C           *** WITH Z-FUNCTION ***
C
      FUNCTION CFQZ_Z(Q,CZ)
C
      INCLUDE 'dpcomm.h'
C
      DATA CI/(0.D0,1.D0)/
C
      NUMAX=NINT(Q-3.D0/2.D0)
      CSUM=(0.D0,0.D0)
      DO NU=0,NUMAX 
         CTERM=(-CZ)**NU*DGAMM(Q-1-NU)
         CSUM=CSUM+CTERM
      ENDDO
      CGZ=CI*SQRT(CZ)
      CALL DSPFNA(1,CGZ,CZ2,CDZ2)
      CSUM=CSUM+SQRT(PI)*(-CZ)**NUMAX*(CI*SQRT(CZ)*CZ2)
      CFQZ_Z=CSUM/DGAMM(Q)
      RETURN
      END
C     
C     ****** CALCULATE SHKAROFSKY ******
C     *** WITH ASYMPTOTIC EXPANSION ***
C
      FUNCTION CFQZ_EXP(Q,CZ)
C
      INCLUDE 'dpcomm.h'
C
      CFQZ=(0.D0,0.D0)
C
         RREZ=DBLE(CZ)
         IF(RREZ.GE.-(11.D0+(Q-5.D0/2.D0)*1.35D0).AND.RREZ.LE.11.D0)THEN
            DO NU=0,50
               CSUM=(-CZ)**NU*DGAMM(Q-1-NU)
               CFQZ=CFQZ+CSUM
            ENDDO
            CFQZ=(CFQZ-PI*(-CZ)**(Q-3.D0/2.D0)*SQRT(CZ)*EXP(CZ))
     &           /DGAMM(Q)
         ELSE
            RIMZ=DIMAG(CZ)
            IF(RIMZ.EQ.0.D0.AND.RREZ.LT.0.D0)THEN
               SIG=1.D0
            ELSE
               SIG=0.D0
            ENDIF
            DO NU=0,10
               CSUM=DGAMM(Q+NU)*((-CZ)**(-1-NU))
               CFQZ=CFQZ+CSUM
            ENDDO
            CFQZ=(CFQZ+(CFQZ+DGAMM(Q+(NU+1))*(-CZ)**(-1-NU-1)))/2.D0
            CFQZ_EXP
     &          =-(CFQZ-PI*SIG*(-CZ)**(Q-3.D0/2.D0)*SQRT(CZ)*EXP(CZ))
     &           /DGAMM(Q)
         ENDIF
      RETURN
      END
C
C****************************************************
C                   gamma function                  *
C                in double precision                *
C      COPYRIGHT : M.Mori  JUNE 30 1989  V.1        *
C****************************************************
C
      FUNCTION DGAMM(X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION C(0:19)
      DATA IN / 19 /
C
C     ---- for single precision ----
C     DATA IN / 10 /
C
      DATA C /  1.0                   D0,
     &         -0.42278 43350 98467 1 D0,
     &         -0.23309 37364 21786 7 D0,
     &          0.19109 11013 87691 5 D0,
     &         -0.24552 49000 54000 2 D-1,
     &         -0.17645 24455 01443 2 D-1,
     &          0.80232 73022 26734 7 D-2,
     &         -0.80432 97756 04247 0 D-3,
     &         -0.36083 78162 548     D-3,
     &          0.14559 61421 399     D-3,
     &         -0.17545 85975 17      D-4,
     &         -0.25889 95022 4       D-5,
     &          0.13385 01546 6       D-5,
     &         -0.20547 43152         D-6,
     &         -0.15952 68            D-9,
     &          0.62756 218           D-8,
     &         -0.12736 143           D-8,
     &          0.92339 7             D-10,
     &          0.12002 8             D-10,
     &         -0.42202               D-11 /
C
      IF (X .GT. 57.0D0) GO TO 901
C
      XX = X
      IF (XX .LE. 1.5D0) THEN
        IF (XX .GE. 0.5D0) THEN
          A = XX - 1.0D0
          FCTR = 1.0D0
        ELSE
          M = INT(XX)
          A = XX - M
          IF (A .EQ. 0.0D0) THEN
            GO TO 902
          ELSE IF (A .GE. -0.5D0) THEN
            MG = IABS(M) + 1
          ELSE
            MG = IABS(M) + 2
            A = A + 1.0D0
          END IF
          Z = 1.0D0
          DO 20 I = 1, MG
            Z = Z * XX
            XX = XX + 1.0D0
   20     CONTINUE
          FCTR = 1.0D0 / Z
        END IF
C
      ELSE
        M = INT (XX)
        A = XX - M
        IF (A .LE. 0.5D0) THEN
          MG = M - 1
        ELSE
          MG = M
          A = A - 1.0D0
        END IF
        Z = 1.0D0
        DO 30 I = 1, MG
          Z = Z * (XX - 1.0D0)
          XX = XX - 1.0D0
   30   CONTINUE
        FCTR = Z
      END IF
C
      Y = C(IN)
      DO 10 I = IN - 1, 0, -1
        Y = C(I) + A * Y
   10 CONTINUE
C
      DGAMM = FCTR / ((1.0D0 + A) * Y)
      RETURN
C
  901 CONTINUE
      WRITE (6,2001) X
 2001 FORMAT (' (FUNC.DGAMM) X(=',D23.16,')',
     &        ' must be smaller than 57.0')
      DGAMM = 1.0D75
      RETURN
C
  902 CONTINUE
      WRITE (6,2002) X
 2002 FORMAT (' (FUNC.DGAMM) invalid argument',
     &        ' X =',D23.16)
      DGAMM = 1.0D75
      RETURN
C
      END
