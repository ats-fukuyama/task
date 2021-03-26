! wmsetf0.f90

MODULE wmsetf0

  PRIVATE
  PUBLIC wm_setf0

CONTAINS

!     ****** CALCULATE LOCAL MATRIX ******

  SUBROUTINE wm_setf0(NR,NS0)

    USE wmcomm
    USE wmprof
    USE wmdisp
    USE wmsub
    USE libinv
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS0
    COMPLEX(rkind),ALLOCATABLE:: &
         CEP0(:,:,:,:,:,:), &
         CRA(:,:,:,:),CFA(:,:,:,:), &
         CRB(:,:,:,:,:,:),CFB(:,:,:,:,:,:), &
         CRC(:,:,:,:),CFC(:,:,:,:), &
         CS(:,:,:,:),CSUMA(:,:,:,:,:,:)
    REAL(rkind):: RMA(3,3),RMB(3,3),RGA(3,3),RGB(3,3)
    COMPLEX(rkind):: CW,CWC2,CSUM
    REAL(rkind):: XRI,XRL,BABS,BSUPTH,BSUPPH,TC2,TC3,RF11
    INTEGER:: ND,MD,NS1,NS2,NS,NHH,NTH,I,J,K,L,KD,LD,KDX,LDX,NDX,MDX
    INTEGER:: KA,KAX,KB,KBX,NB1,NB2,LA,LAX,LB,LBX,MB1,MB2,LAB,LABX,KAB,KABX
    INTEGER:: KC,KCX,LC,LCX
    INTEGER:: ILL

    ALLOCATE(CEP0(3,3,-nthmax_f:nthmax_f,-nhhmax_f:nhhmax_f,nthmax_f,nhhmax_f))
    ALLOCATE(CRA(nthmax_f,nhhmax_f,3,3),CFA(nthmax_f,nhhmax_f,3,3))
    ALLOCATE(CRB(nthmax_f,nhhmax_f,3,3,-nthmax_f:nthmax_f,-nhhmax_f:nhhmax_f))
    ALLOCATE(CFB(nthmax_f,nhhmax_f,3,3,-nthmax_f:nthmax_f,-nhhmax_f:nhhmax_f))
    ALLOCATE(CRC(nthmax_f,nhhmax_f,3,3),CFC(nthmax_f,nhhmax_f,3,3))
    ALLOCATE(CSUMA(3,3,nthmax_f,nthmax_f,nhhmax_f,nhhmax_f))
    ALLOCATE(CS(3,3,nthmax_f,nhhmax_f))

    CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
    CWC2=CW**2/VC**2

    IF(NR.EQ.1) THEN
       XRI=1.D6/XRHO(2)
       XRL=XRHO(2)/1.D6
    ELSE
       XRI=1.D0/XRHO(NR)
       XRL=XRHO(NR)
    ENDIF

!     ----- Fourier decompose metric tensor g -----

    CALL WMSUBG_F(RG11(:,:,NR),RJ(:,:,NR),CGF11(:,:,3))
    CALL WMSUBG_F(RG12(:,:,NR),RJ(:,:,NR),CGF12(:,:,3))
    CALL WMSUBG_F(RG13(:,:,NR),RJ(:,:,NR),CGF13(:,:,3))
    CALL WMSUBG_F(RG22(:,:,NR),RJ(:,:,NR),CGF22(:,:,3))
    CALL WMSUBG_F(RG23(:,:,NR),RJ(:,:,NR),CGF23(:,:,3))
    CALL WMSUBG_F(RG33(:,:,NR),RJ(:,:,NR),CGF33(:,:,3))

!     ----- Calculate dielectric tensor -----

    DO NHH=1,NHHMAX_F
       DO NTH=1,NTHMAX_F
          DO ND=NDMIN_F,NDMAX_F
             DO MD=MDMIN_F,MDMAX_F
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
       END DO
    END DO

    IF(NS0.EQ.0) THEN
       NS1=1
       NS2=NSMAX
    ELSE
       NS1=NS0
       NS2=NS0
    ENDIF

    DO NS=NS1,NS2
       DO ND=NDMIN_F,NDMAX_F
          DO MD=MDMIN_F,MDMAX_F
             CALL wm_tnsr(NR,NS,MD,ND)
             DO NHH=1,NHHMAX_F
                DO NTH=1,NTHMAX_F
                   DO J=1,3
                      DO I=1,3
                         CEP0(I,J,MD,ND,NTH,NHH) &
                              =CEP0(I,J,MD,ND,NTH,NHH) &
                              +CTNSR(I,J,NTH,NHH)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          END DO
       END DO
    END DO

!        ----- Calculate rotation matrix mu=RMA -----

    DO NHH=1,NHHMAX_F
       DO NTH=1,NTHMAX_F

          CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
          TC2=BSUPTH/BABS
          TC3=BSUPPH/BABS

!        ***** RF11=RJ*SQRT(G^11)/XR *****

          RF11=SQRT(RG22(NTH,NHH,NR)*RG33(NTH,NHH,NR) &
                   -RG23(NTH,NHH,NR)*RG23(NTH,NHH,NR))
          RMA(1,1)= RJ(NTH,NHH,NR)/RF11*XRI
          RMA(2,1)= 0.D0
          RMA(3,1)= 0.D0
          RMA(1,2)= (TC2*(RG23(NTH,NHH,NR)*RG12(NTH,NHH,NR) &
                         -RG22(NTH,NHH,NR)*RG13(NTH,NHH,NR)) &
                    +TC3*(RG33(NTH,NHH,NR)*RG12(NTH,NHH,NR) &
                         -RG23(NTH,NHH,NR)*RG13(NTH,NHH,NR))*XRI) &
                    /RF11
          RMA(2,2)= TC3*RF11*XRL
          RMA(3,2)=-TC2*RF11*XRL
          RMA(1,3)=TC2*RG12(NTH,NHH,NR) &
                  +TC3*RG13(NTH,NHH,NR)*XRI
          RMA(2,3)=TC2*RG22(NTH,NHH,NR)*XRL*XRL &
                  +TC3*RG23(NTH,NHH,NR)*XRL
          RMA(3,3)=TC2*RG23(NTH,NHH,NR)*XRL &
                  +TC3*RG33(NTH,NHH,NR)

!        ----- Set metric matrix g=RGA -----

          RGA(1,1)=RG11(NTH,NHH,NR)
          RGA(1,2)=RG12(NTH,NHH,NR)
          RGA(1,3)=RG13(NTH,NHH,NR)
          RGA(2,1)=RG12(NTH,NHH,NR)
          RGA(2,2)=RG22(NTH,NHH,NR)
          RGA(2,3)=RG23(NTH,NHH,NR)
          RGA(3,1)=RG13(NTH,NHH,NR)
          RGA(3,2)=RG23(NTH,NHH,NR)
          RGA(3,3)=RG33(NTH,NHH,NR)

!        ----- Invert matrix to obtain mu^(-1)=RMB and g^(-1)=RGB -----

          DO J=1,3
             DO I=1,3
                RMB(I,J)=RMA(I,J)
                RGB(I,J)=RGA(I,J)
             ENDDO
          ENDDO

          CALL INVMRD(RMB,3,3,ILL)
          IF(ILL.NE.0) THEN
             WRITE(6,*) 'XX WMSETF: INVMRD(RMB) : SINGULAR MATRIX'
             WRITE(6,'(3I5,1P2E12.4)') NR,NTH,NHH,TC3,TC2
             GOTO 9000
          ENDIF

!        ----- RGB = g^(-1) calculated in wmsetg

          CALL INVMRD(RGB,3,3,ILL)
          IF(ILL.NE.0) THEN
             WRITE(6,*) 'XX WMSETF: INVMRD(RGB) : SINGULAR MATRIX'
             WRITE(6,'(3I5)') NR,NTH,NHH
             GOTO 9000
          ENDIF

          RGB(1,1)=RGB(1,1)*XRL**2
          RGB(1,2)=RGB(1,2)
          RGB(1,3)=RGB(1,3)*XRL   
          RGB(2,1)=RGB(2,1)
          RGB(2,2)=RGB(2,2)*XRI**2
          RGB(2,3)=RGB(2,3)*XRI   
          RGB(3,1)=RGB(3,1)*XRL   
          RGB(3,2)=RGB(3,2)*XRI   
          RGB(3,3)=RGB(3,3)

!        ----- Setup Matrix A=CRA, B=CRB, C=CRC -----

          DO J=1,3
             DO I=1,3
                CSUM=0.D0
                DO K=1,3
                   CSUM=CSUM+RGB(I,K)*RMA(K,J)
                ENDDO
                CRA(NTH,NHH,I,J)=CSUM*RJ(NTH,NHH,NR)
             ENDDO
          ENDDO

          DO ND=NDMIN_F,NDMAX_F
             DO MD=MDMIN_F,MDMAX_F
                DO J=1,3
                   DO I=1,3
                      CRB(NTH,NHH,I,J,MD,ND)=CEP0(I,J,MD,ND,NTH,NHH)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          DO J=1,3
             DO I=1,3
                CRC(NTH,NHH,I,J)=RMB(I,J)
             ENDDO
          ENDDO

       END DO ! NTH
    END DO ! NHH
          
!     ----- Fourier decompose A, B, C -----

    DO J=1,3
       DO I=1,3
          CALL WMSUBF_F(CRA(:,:,I,J),CFA(:,:,I,J))
          CALL WMSUBF_F(CRC(:,:,I,J),CFC(:,:,I,J))
       ENDDO
    ENDDO

    DO ND=NDMIN_F,NDMAX_F
       DO MD=MDMIN_F,MDMAX_F
          DO J=1,3
             DO I=1,3
                CALL WMSUBF_F(CRB(:,:,I,J,MD,ND),CFB(:,:,I,J,MD,ND))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

!     ----- Fourier decompose A, B, C -----

    DO ND=NDMIN,NDMAX
       NDX=ND-NDMIN+1
       DO MD=MDMIN,MDMAX
          MDX=MD-MDMIN+1
          DO KDX=1,nhhmax_f
             DO LDX=1,nthmax_f
                DO J=1,3
                   DO I=1,3
                      CSUMA(I,J,LDX,MDX,KDX,NDX)=0d0
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO ND=NDMIN,NDMAX
       NDX=ND-NDMIN+1
       DO MD=MDMIN,MDMAX
          MDX=MD-MDMIN+1

          DO KAB=KDMIN_F,KDMAX_F
             KABX=KAB-KDMIN_F+1
             DO LAB=LDMIN_F,LDMAX_F
                LABX=LAB-LDMIN_F+1
                DO K=1,3
                   DO I=1,3
                      CS(I,K,LABX,KABX)=0.D0
                   ENDDO
                ENDDO
             ENDDO
          END DO
          
          DO KA=KDMIN_F,KDMAX_F
             KAX=KA-KDMIN_F+1
             DO KB=KDMIN_F,KDMAX_F
                KBX=KB-KDMIN_F+1

                KAB=KA+KB
                IF(KAB.GE.KDMIN_F.AND.(KAB.LE.KDMAX_F.OR.KDMAX_F.EQ.0))THEN
                   KABX=KAB-KDMIN_F+1

                   IF(MOD(KB,2).EQ.0) THEN
                      NB1=ND-(KA+KB/2)
                      NB2=ND-(KA+KB/2)
                   ELSE
                      NB1=ND-(KA+(KB-1)/2)
                      NB2=ND-(KA+(KB+1)/2)
                   ENDIF

                   DO LA=LDMIN_F,LDMAX_F
                      LAX=LA-LDMIN_F+1
                      DO LB=LDMIN_F,LDMAX_F
                         LBX=LB-LDMIN_F+1

                         LAB=LA+LB
                         IF(LAB.GE.LDMIN_F.AND. &
                           (LAB.LE.LDMAX_F.OR.LDMAX_F==0))THEN
                            LABX=LAB-LDMIN_F+1

                            IF(MOD(LB,2).EQ.0) THEN
                               MB1=MD-(LA+LB/2)
                               MB2=MD-(LA+LB/2)
                            ELSE
                               MB1=MD-(LA+(LB-1)/2)
                               MB2=MD-(LA+(LB+1)/2)
                            ENDIF

                            DO K=1,3
                               DO I=1,3
                                  DO L=1,3
                                     CS(I,K,LABX,KABX) &
                                          =CS(I,K,LABX,KABX) &
                                          +CFA(LAX,KAX,I,L) &
                                          *0.25D0*( &
                                              CFB(LBX,KBX,L,K,MB1,NB1) &
                                             +CFB(LBX,KBX,L,K,MB1,NB2) &
                                             +CFB(LBX,KBX,L,K,MB2,NB1) &
                                             +CFB(LBX,KBX,L,K,MB2,NB2))
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

          DO KD=KDMIN_F,KDMAX_F
             KDX=KD-KDMIN+1
             DO KAB=KDMIN_F,KDMAX_F
                KABX=KAB-KDMIN_F+1
                KC=KD-KAB
                IF(KC.GE.KDMIN_F.AND.KC.LE.KDMAX_F) THEN
                   KCX=MOD(KC-KDMIN_F+2*KDSIZ_F,KDSIZ_F)+1

                   DO LD=LDMIN_F,LDMAX_F
                      LDX=LD-LDMIN_F+1
                      DO LAB=LDMIN_F,LDMAX_F
                         LABX=LAB-LDMIN_F+1
                         LC=LD-LAB
                         IF(LC.GE.LDMIN_F.AND.LC.LE.LDMAX_F) THEN
                            LCX=MOD(LC-LDMIN_F+2*LDSIZ_F,LDSIZ_F)+1
                            
                            DO J=1,3
                               DO I=1,3
                                  DO K=1,3
                                     CSUMA(I,J,LDX,MDX,KDX,NDX) &
                                    =CSUMA(I,J,LDX,MDX,KDX,NDX) &
                                    +CS(I,K,LABX,KABX)*CFC(LCX,KCX,K,J)
                                  ENDDO
                               ENDDO
                            ENDDO
                            
                         ENDIF
                      ENDDO
                   ENDDO
                   
                ENDIF
             ENDDO
          ENDDO
          
          DO KD=KDMIN_F,KDMAX_F
             KDX=KD-KDMIN_F+1
             DO LD=LDMIN_F,LDMAX_F
                LDX=LD-LDMIN_F+1
                DO J=1,3
                   DO I=1,3
                      CGD(I,J,LDX,MDX,KDX,NDX,3) &
                           =-CWC2*CSUMA(I,J,LDX,MDX,KDX,NDX)
                      CPSF(I,J,LDX,MDX,KDX,NDX,3) &
                           =-CWC2*CSUMA(I,J,LDX,MDX,KDX,NDX)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       ENDDO
    ENDDO

    RETURN

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_setf0

END MODULE wmsetf0
