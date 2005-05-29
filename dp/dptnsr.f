C     $Id$
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE DPTENS(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION CLDISP(6),CLDISP1(6)
C
      IF(RN(NS).LE.0.D0) THEN
         CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP)
      ELSE
         ID1=MOD(MODELP(NS),10)
         ID2=MODELP(NS)/10
C
         IF(ID2.EQ.2) THEN
            CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP1)
         ELSEIF(ID2.EQ.3) THEN
            CALL DPTNKP(CW,CKPR,CKPP,NS,CLDISP1)
         ENDIF
C
         IF(ID1.EQ.0) THEN
            CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.1) THEN
            CALL DPTNCC(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.2) THEN
            CALL DPTNIM(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.3) THEN
            CALL DPTNRM(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.4) THEN
            CALL DPTNHP(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.5) THEN
            CALL DPTNKP(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.6) THEN
            CALL DPTNKR(CW,CKPR,CKPP,NS,CLDISP)
         ELSE IF(ID1.EQ.7) THEN
            IF(MODELV.EQ.0) THEN
               CALL DPFMFL(NS)
               IF(ID2.EQ.2) THEN
                  CALL DPHOTFI(CW,CKPR,CKPP,NS,CLDISP)
               ELSE
                  CALL DPHOTF(CW,CKPR,CKPP,NS,CLDISP)
               ENDIF
            ELSE IF(MODELV.EQ.1) THEN
               CALL DPFPFL(NS)
               IF(ID2.EQ.2) THEN
                  CALL DPHOTFI(CW,CKPR,CKPP,NS,CLDISP)
               ELSE
                  CALL DPHOTF(CW,CKPR,CKPP,NS,CLDISP)
               ENDIF
            ELSE IF(MODELV.EQ.2) THEN
               CALL DPFMFL(NS)
               IF(ID2.EQ.2) THEN
                  CALL DPHOTRI(CW,CKPR,CKPP,NS,CLDISP)
               ELSE
                  CALL DPHOTR(CW,CKPR,CKPP,NS,CLDISP)
               ENDIF
            ELSE IF(MODELV.EQ.3) THEN
               CALL DPFPFL(NS)
               IF(ID2.EQ.2) THEN
                  CALL DPHOTRI(CW,CKPR,CKPP,NS,CLDISP)
               ELSE
                  CALL DPHOTR(CW,CKPR,CKPP,NS,CLDISP)
               ENDIF
            ENDIF
         ENDIF
C
         IF(ID2.EQ.2) THEN
            CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP1)
         ELSEIF(ID2.EQ.3) THEN
            CALL DPTNKP(CW,CKPR,CKPP,NS,CLDISP1)
         ENDIF
C
         IF(ID2.EQ.2.OR.ID2.eq.3) THEN
            CLDISP(1)=DBLE(CLDISP1(1))+CI*DIMAG(CLDISP(1))
            CLDISP(2)=DBLE(CLDISP1(2))+CI*DIMAG(CLDISP(2))
            CLDISP(3)=DBLE(CLDISP1(3))+CI*DIMAG(CLDISP(3))
            CLDISP(4)=DBLE(CLDISP1(4))+CI*DIMAG(CLDISP(4))
            CLDISP(5)=CI*DIMAG(CLDISP1(5))+DBLE(CLDISP(5))
            CLDISP(6)=CI*DIMAG(CLDISP1(6))+DBLE(CLDISP(6))
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
C     ****** COLLISIONLESS COLD MODEL ******
C
      SUBROUTINE DPTNCL(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)= 0.D0
      RETURN
      END
C
C     ****** COLLISIONAL COLD MODEL ******
C
      SUBROUTINE DPTNCC(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
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
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)= 0.D0
      RETURN
      END
C
C     ****** IDEAL MHD MODEL ******
C
      SUBROUTINE DPTNIM(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS))
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
C
      CLDISP(1)= CWP/(CWC*CWC)
      CLDISP(2)=-1.D16-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)= 0.D0
      CLDISP(6)= 0.D0
      RETURN
      END
C
C     ****** RESISTIVE MHD MODEL ******
C
      SUBROUTINE DPTNRM(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      RNYU=PZCL(NS)
      CWP=CWP/DCMPLX(1.D0,RNYU)
      CWC=CWC/DCMPLX(1.D0,RNYU)
C
      CLDISP(1)= CWP/(CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)= 0.D0
      CLDISP(6)= 0.D0
      RETURN
      END
C
C     ****** WARM PLAMSA MODEL ******
C
      SUBROUTINE DPTNWP(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
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
C     ******  UNMAGNETIZE KINETIC DISPERSION ******
C
      SUBROUTINE DPTNUP(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
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
C     ****** KINETIC MODEL WITHOUT FLR ******
C
      SUBROUTINE DPTNHP(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
      DIMENSION CGZ(-1:1),CZ(-1:1),CDZ(-1:1)
C
      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO
C
      CWP=RN(NS)*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
C
      WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPP=RTPP(NS)*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
C
      CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
      CGZ(-1)= (1.D0+CWC)*CPR
      CGZ( 0)=            CPR
      CGZ(+1)= (1.D0-CWC)*CPR
      CALL DSPFNA(3,CGZ,CZ,CDZ)
C
      CK=CKPP/CKPR
C
      CFWM=CPR*CZ(-1)+0.5D0*(1.D0-WTPX)*CDZ(-1)
      CFW0=CPR*CZ( 0)+0.5D0*(1.D0-WTPX)*CDZ( 0)
      CFWP=CPR*CZ(+1)+0.5D0*(1.D0-WTPX)*CDZ(+1)
      CFFM=0.5D0*(WTPX/CWC+(WTPX-1.D0))*CDZ(-1)
      CFF0=0.5D0* WTPX/CWC             *CDZ( 0)
      CFFP=0.5D0*(WTPX/CWC-(WTPX-1.D0))*CDZ(+1)
      CFGM=+CPR/WTPX*CGZ(-1)*CDZ(-1)
      CFG0=-CPR     *CGZ( 0)*CDZ( 0)
      CFGP=-CPR/WTPX*CGZ(+1)*CDZ(+1)
C
      CLAM=CALAM(ABS(0))
      CLAMM=CALAM(ABS(-1))
      CLAMP=CALAM(ABS(+1))
      CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
      CBLAM=0.5D0*(CLAMM-CLAMP)
C
      CLDISP(1)=CLDISP(1)-CWP*      CBLAM*CFWM
      CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFGM+CBLAM*CFWM)
      CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFWM
      CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFFM
      CLDISP(5)=CLDISP(5)-CWP*CI*   CDLAM*CFWM
      CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFFM
C
      CLDISP(2)=CLDISP(2)+CWP*      CLAM*CFG0
      CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW0
      CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF0
      CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF0
C
      CLDISP(1)=CLDISP(1)+CWP*      CBLAM*CFWP
      CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFGP-CBLAM*CFWP)
      CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFWP
      CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFFP
      CLDISP(5)=CLDISP(5)+CWP*CI*   CDLAM*CFWP
      CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFFP
C
      RETURN
      END
C
C     ****** KINETIC MODEL WITH LOWEST FLR ******
C
      SUBROUTINE DPTNKL(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
C
      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO
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
      DO NC=-2,2
C
         CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
         CGZ= (1.D0-NC*CWC)*CPR
         CALL DSPFNA(1,CGZ,CZ,CDZ)
C
         CK=CKPP/CKPR
C
         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-NC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-NC*(1.D0-1.D0/WTPX))*CGZ*CDZ
C
         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)
C
         CLDISP(1)=CLDISP(1)+CWP*   NC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-NC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
      ENDDO
      RETURN
      END
C
C     ****** KINETIC MODEL WITH FLR ******
C
      SUBROUTINE DPTNKP(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
C
      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO
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
      DO NC=NMIN,NMAX
C
         IF(ABS(CKPR).LE.0.D0) THEN
            CPR=CW/SQRT(2.D0*1.D-4**2*WTPR)
            CK=CKPP/1.D-4
         ELSE
            CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
            CK=CKPP/CKPR
        ENDIF
         CGZ= (1.D0-NC*CWC)*CPR
         CALL DSPFNA(1,CGZ,CZ,CDZ)
C
         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-NC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-NC*(1.D0-1.D0/WTPX))*CGZ*CDZ
C
         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)
C
         CLDISP(1)=CLDISP(1)+CWP*   NC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-NC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
      ENDDO
      RETURN
      END
C
C     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******
C     ******                    IN                      ******
C     ******     WEAKLY RELATIVISTIC THERMAL PLASMA     ******
C
      SUBROUTINE DPTNKR(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION CLDISP(6)
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
      INCLUDE '../dp/dpcomm.inc'
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
      INCLUDE '../dp/dpcomm.inc'
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
      INCLUDE '../dp/dpcomm.inc'
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
          DO I = 1, MG
            Z = Z * XX
            XX = XX + 1.0D0
          ENDDO
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
        DO I = 1, MG
          Z = Z * (XX - 1.0D0)
          XX = XX - 1.0D0
       ENDDO
        FCTR = Z
      END IF
C
      Y = C(IN)
      DO I = IN - 1, 0, -1
        Y = C(I) + A * Y
      ENDDO
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
