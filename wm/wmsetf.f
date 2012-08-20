C     $Id$
C
C     ****** CALCULATE LOCAL MATRIX ******
C
      SUBROUTINE WMSETF(NR,NS0)
C
      INCLUDE 'wmcomm.inc'
      DIMENSION RMA(3,3),RMB(3,3),RGA(3,3),RGB(3,3)
      DIMENSION CEP0(3,3,-MDMX:MDMX,-NDMX:NDMX,MDM,NDM)
      DIMENSION CRA(MDM,NDM,3,3),CFA(MDM,NDM,3,3)
      DIMENSION CRB(MDM,NDM,3,3,-MDMX:MDMX,-NDMX:NDMX)
      DIMENSION CFB(MDM,NDM,3,3,-MDMX:MDMX,-NDMX:NDMX)
      DIMENSION CRC(MDM,NDM,3,3),CFC(MDM,NDM,3,3)
      DIMENSION CSUMA(3,3,MDM,MDM,NDM,NDM),CS(3,3,MDM,NDM)
C
C      DATA IDEBUG/0/
C
      CW=2.D0*PI*CRF*1.D6
      CWC2=CW**2/VC**2
C
      IF(NR.EQ.1) THEN
         XRI=1.D6/XRHO(2)
         XRL=XRHO(2)/1.D6
      ELSE
         XRI=1.D0/XRHO(NR)
         XRL=XRHO(NR)
      ENDIF
C
C     ----- Fourier decompose metric tensor g -----
C
      CALL WMSUBG(RG11(1,1,NR),RJ(1,1,NR),CGF11(1,1,3))
      CALL WMSUBG(RG12(1,1,NR),RJ(1,1,NR),CGF12(1,1,3))
      CALL WMSUBG(RG13(1,1,NR),RJ(1,1,NR),CGF13(1,1,3))
      CALL WMSUBG(RG22(1,1,NR),RJ(1,1,NR),CGF22(1,1,3))
      CALL WMSUBG(RG23(1,1,NR),RJ(1,1,NR),CGF23(1,1,3))
      CALL WMSUBG(RG33(1,1,NR),RJ(1,1,NR),CGF33(1,1,3))
C
C     ----- Calculate dielectric tensor -----
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
      DO ND=-NDSIZX,NDSIZX
      DO MD=-MDSIZX,MDSIZX
         DO J=1,3
         DO I=1,3
            CEP0(I,J,MD,ND,NTH,NHH)=0.D0
         ENDDO
         ENDDO
         IF(NS0.EQ.0) THEN
            CEP0(1,1,MD,ND,NTH,NHH)=1.D0
            CEP0(2,2,MD,ND,NTH,NHH)=1.D0
            CEP0(3,3,MD,ND,NTH,NHH)=1.D0
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      IF(NS0.EQ.0) THEN
         NS1=1
         NS2=NSMAX
      ELSE
         NS1=NS0
         NS2=NS0
      ENDIF
C
      DO NS=NS1,NS2
         CALL WMTNSR(NR,NS)
C
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            DO ND=-NDSIZX,NDSIZX
            DO MD=-MDSIZX,MDSIZX
               DO J=1,3
               DO I=1,3
                  CEP0(I,J,MD,ND,NTH,NHH)
     &           =CEP0(I,J,MD,ND,NTH,NHH)
     &           +CTNSR(I,J,MD,ND,NTH,NHH)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
C     ----- Calculate A, B, C -----
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
C
C        ----- Calculate rotation matrix mu=RMA -----
C
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         TC2=BSUPTH/BABS
         TC3=BSUPPH/BABS
C
C        ***** RF11=RJ*SQRT(G^11)/XR *****
C
         RF11=SQRT(RG22(NTH,NHH,NR)*RG33(NTH,NHH,NR)
     &            -RG23(NTH,NHH,NR)*RG23(NTH,NHH,NR))
         RMA(1,1)= RJ(NTH,NHH,NR)/RF11*XRI
         RMA(2,1)= 0.D0
         RMA(3,1)= 0.D0
         RMA(1,2)= (TC2*(RG23(NTH,NHH,NR)*RG12(NTH,NHH,NR)
     &                  -RG22(NTH,NHH,NR)*RG13(NTH,NHH,NR))
     &             +TC3*(RG33(NTH,NHH,NR)*RG12(NTH,NHH,NR)
     &                  -RG23(NTH,NHH,NR)*RG13(NTH,NHH,NR))*XRI)
     &             /RF11
         RMA(2,2)= TC3*RF11*XRL
         RMA(3,2)=-TC2*RF11*XRL
         RMA(1,3)=TC2*RG12(NTH,NHH,NR)
     &           +TC3*RG13(NTH,NHH,NR)*XRI
         RMA(2,3)=TC2*RG22(NTH,NHH,NR)*XRL*XRL
     &           +TC3*RG23(NTH,NHH,NR)*XRL
         RMA(3,3)=TC2*RG23(NTH,NHH,NR)*XRL
     &           +TC3*RG33(NTH,NHH,NR)
C
C         IF(NR.EQ.2.AND.NS0.EQ.0) THEN
C            WRITE(6,*) 'NR,NTH,NHH=',NR,NTH,NHH
C            WRITE(6,'(1P3E12.4)') RMA(1,1),RMA(1,2),RMA(1,3)
C            WRITE(6,'(1P3E12.4)') RMA(2,1),RMA(2,2),RMA(2,3)
C            WRITE(6,'(1P3E12.4)') RMA(3,1),RMA(3,2),RMA(3,3)
C            WRITE(6,'(1P3E12.4)') RJ(NTH,NHH,NR),RF11,XRI
C         ENDIF
C
C        ----- Set metric matrix g=RGA -----
C
         RGA(1,1)=RG11(NTH,NHH,NR)
         RGA(1,2)=RG12(NTH,NHH,NR)
         RGA(1,3)=RG13(NTH,NHH,NR)
         RGA(2,1)=RG12(NTH,NHH,NR)
         RGA(2,2)=RG22(NTH,NHH,NR)
         RGA(2,3)=RG23(NTH,NHH,NR)
         RGA(3,1)=RG13(NTH,NHH,NR)
         RGA(3,2)=RG23(NTH,NHH,NR)
         RGA(3,3)=RG33(NTH,NHH,NR)
C
C        ----- Invert matrix to obtain mu^(-1)=RMB and g^(-1)=RGB -----
C
         DO J=1,3
         DO I=1,3
            RMB(I,J)=RMA(I,J)
            RGB(I,J)=RGA(I,J)
         ENDDO
         ENDDO
C
         CALL INVMRD(RMB,3,3,ILL)
         IF(ILL.NE.0) THEN
            WRITE(6,*) 'XX WMSETF: INVMRD(RMB) : SINGULAR MATRIX'
            WRITE(6,'(3I5,1P2E12.4)') NR,NTH,NHH,TC3,TC2
            GOTO 9000
         ENDIF
C
         CALL INVMRD(RGB,3,3,ILL)
         IF(ILL.NE.0) THEN
            WRITE(6,*) 'XX WMSETF: INVMRD(RGB) : SINGULAR MATRIX'
            WRITE(6,*) '   NTH,NHH,NR = ',NTH,NHH,NR
            GOTO 9000
         ENDIF
C
C        ----- RGB = g^(-1)
C
         RGB(1,1)=RGB(1,1)*XRL**2
         RGB(1,2)=RGB(1,2)
         RGB(1,3)=RGB(1,3)*XRL   
         RGB(2,1)=RGB(2,1)
         RGB(2,2)=RGB(2,2)*XRI**2
         RGB(2,3)=RGB(2,3)*XRI   
         RGB(3,1)=RGB(3,1)*XRL   
         RGB(3,2)=RGB(3,2)*XRI   
         RGB(3,3)=RGB(3,3)
C
C        ----- Setup Matrix A=CRA, B=CRB, C=CRC -----
C
         DO J=1,3
         DO I=1,3
            CSUM=0.D0
            DO K=1,3
               CSUM=CSUM+RGB(I,K)*RMA(K,J)
            ENDDO
            CRA(NTH,NHH,I,J)=CSUM*RJ(NTH,NHH,NR)
         ENDDO
         ENDDO
C
         DO ND=-NDSIZX,NDSIZX
         DO MD=-MDSIZX,MDSIZX
            DO J=1,3
            DO I=1,3
               CRB(NTH,NHH,I,J,MD,ND)=CEP0(I,J,MD,ND,NTH,NHH)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
         DO J=1,3
         DO I=1,3
            CRC(NTH,NHH,I,J)=RMB(I,J)
         ENDDO
         ENDDO
C
      ENDDO
      ENDDO
C
C     ----- Fourier decompose A, B, C -----
C
      DO J=1,3
      DO I=1,3
         CALL WMSUBF(CRA(1,1,I,J),CFA(1,1,I,J))
      ENDDO
      ENDDO
C
      DO ND=-NDSIZX,NDSIZX
      DO MD=-MDSIZX,MDSIZX
      DO J=1,3
      DO I=1,3
         CALL WMSUBF(CRB(1,1,I,J,MD,ND),CFB(1,1,I,J,MD,ND))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO J=1,3
      DO I=1,3
         CALL WMSUBF(CRC(1,1,I,J),CFC(1,1,I,J))
      ENDDO
      ENDDO
C
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO KD=KDMIN,KDMAX
         KDX=KD-KDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO LD=LDMIN,LDMAX
         LDX=LD-LDMIN+1
      DO J=1,3
      DO I=1,3
         CSUMA(I,J,LDX,MDX,KDX,NDX)=0.D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
C
         DO KAB=KDMIN,KDMAX
            KABX=KAB-KDMIN+1
         DO LAB=LDMIN,LDMAX
            LABX=LAB-LDMIN+1
         DO K=1,3
         DO I=1,3
            CS(I,K,LABX,KABX)=0.D0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
         DO KA=KDMIN,KDMAX
            KAX=KA-KDMIN+1
         DO KB=KDMIN,KDMAX
            KBX=KB-KDMIN+1
C
            KAB=KA+KB
            IF(MODELK.EQ.0.OR.
     &         (KAB.GE.KDMIN.AND.KAB.LE.KDMAX)) THEN
               KABX=MOD(KAB-KDMIN+2*KDSIZ,KDSIZ)+1
C
               IF(MOD(KB,2).EQ.0) THEN
                  NB1=ND-(KA+KB/2)
                  NB2=ND-(KA+KB/2)
               ELSE
                  NB1=ND-(KA+(KB-1)/2)
                  NB2=ND-(KA+(KB+1)/2)
               ENDIF
C
         DO LA=LDMIN,LDMAX
            LAX=LA-LDMIN+1
         DO LB=LDMIN,LDMAX
            LBX=LB-LDMIN+1
C
            LAB=LA+LB
            IF(MODELK.EQ.0.OR.
     &         (LAB.GE.LDMIN.AND.LAB.LE.LDMAX)) THEN
               LABX=MOD(LAB-LDMIN+2*LDSIZ,LDSIZ)+1
C
               IF(MOD(LB,2).EQ.0) THEN
                  MB1=MD-(LA+LB/2)
                  MB2=MD-(LA+LB/2)
               ELSE
                  MB1=MD-(LA+(LB-1)/2)
                  MB2=MD-(LA+(LB+1)/2)
               ENDIF
C
                  DO K=1,3
                  DO I=1,3
                  DO L=1,3
                     CS(I,K,LABX,KABX)=CS(I,K,LABX,KABX)
     &                      +CFA(LAX,KAX,I,L)
     &                      *0.25D0*(CFB(LBX,KBX,L,K,MB1,NB1)
     &                              +CFB(LBX,KBX,L,K,MB1,NB2)
     &                              +CFB(LBX,KBX,L,K,MB2,NB1)
     &                              +CFB(LBX,KBX,L,K,MB2,NB2))
                  ENDDO
                  ENDDO
                  ENDDO
C
            ENDIF
         ENDDO
         ENDDO
C
            ENDIF
         ENDDO
         ENDDO
C
         DO KD=KDMIN,KDMAX
            KDX=KD-KDMIN+1
         DO KAB=KDMIN,KDMAX
            KABX=KAB-KDMIN+1
            KC=KD-KAB
            IF(MODELK.EQ.0.OR.
     &         (KC.GE.KDMIN.AND.KC.LE.KDMAX)) THEN
               KCX=MOD(KC-KDMIN+2*KDSIZ,KDSIZ)+1
C
         DO LD=LDMIN,LDMAX
            LDX=LD-LDMIN+1
         DO LAB=LDMIN,LDMAX
            LABX=LAB-LDMIN+1
            LC=LD-LAB
            IF(MODELK.EQ.0.OR.
     &         (LC.GE.LDMIN.AND.LC.LE.LDMAX)) THEN
               LCX=MOD(LC-LDMIN+2*LDSIZ,LDSIZ)+1
C
               DO J=1,3
               DO I=1,3
               DO K=1,3   
                  CSUMA(I,J,LDX,MDX,KDX,NDX)
     &           =CSUMA(I,J,LDX,MDX,KDX,NDX)
     &           +CS(I,K,LABX,KABX)*CFC(LCX,KCX,K,J)
               ENDDO
               ENDDO
               ENDDO
C
            ENDIF
         ENDDO
         ENDDO
C
            ENDIF
         ENDDO
         ENDDO
C
      ENDDO
      ENDDO
C
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO KD=KDMIN,KDMAX
         KDX=KD-KDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO LD=LDMIN,LDMAX
         LDX=LD-LDMIN+1
      DO J=1,3
      DO I=1,3
         CGD(I,J,LDX,MDX,KDX,NDX,3)=-CWC2*CSUMA(I,J,LDX,MDX,KDX,NDX)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
C      IF(NR.EQ.2.AND.NS0.EQ.0) THEN
C         WRITE(6,*) 'CRC  =',CRC(1,1,1,1)
C         WRITE(6,*) 'CFA  =',CFA(1,1,1,1)
C         WRITE(6,*) 'CFB  =',CFB(1,1,1,1,1,1)
C         WRITE(6,*) 'CFC  =',CFC(1,1,1,1)
C         WRITE(6,*) 'CS   =',CS(1,1,1,1)
C         WRITE(6,*) 'CSUMA=',CSUMA(1,1,1,1,1,1)
C         WRITE(6,*) 'CGD  =',CGD(1,1,1,1,1,1,3)
C      ENDIF
C
      RETURN
C
 9000 CONTINUE
      RETURN
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBG(RF1,RF2,CF)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RF1(MDM,NDM),RF2(MDM,NDM),CF(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=RF1(NTH,NHH)/RF2(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF(LDX,NHH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         DO KDX=1,KDSIZ
            CF(LDX,KDX)=CFN(KDX)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBF(CF1,CF2)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CF1(MDM,NDM),CF2(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         DO LDX=1,LDSIZ
            CF2(LDX,NHH)=CFM(LDX)
         ENDDO
      ENDDO
C
      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         DO KDX=1,KDSIZ
            CF2(LDX,KDX)=CFN(KDX)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** INTERFACE FOR FFT ******
C
      SUBROUTINE WMXFFT(CA,N,KEY)
C
      USE libfft,ONLY: FFT2L
      INCLUDE 'wmcomm.inc'
C
      COMPLEX*16 CA(N)
      DATA NS/0/
C
      IF(N.NE.1) THEN
C         LP=NINT(LOG(DBLE(N))/LOG(2.D0))
         IF(N.EQ.NS) THEN
            IND=0
         ELSE
            IND=1
            NS=N
         ENDIF
         IF(KEY.EQ.0) THEN
            CALL FFT2L(CA,CT,RFFT,LFFT,N,IND,KEY)
            DO I=1,N
               IX=I+N/2-1
               IF(IX.GT.N) IX=IX-N
               CA(IX)=CT(I)
            ENDDO
         ELSE
            DO I=1,N
               IX=I+N/2-1
               IF(IX.GT.N) IX=IX-N
               CT(I)=CA(IX)
            ENDDO
            CALL FFT2L(CT,CA,RFFT,LFFT,N,IND,KEY)
         ENDIF
      ENDIF
C
      RETURN
      END
