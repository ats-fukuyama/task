! wmsetf.f90

MODULE wmsetf

  PRIVATE
  PUBLIC wm_setf

CONTAINS

!     ****** CALCULATE LOCAL MATRIX ******

  SUBROUTINE wm_setf(NR1,NS0)

    USE wmcomm
    IMPLICIT NONE

    REAL(rkind):: RMA(3,3),RMB(3,3),RGA(3,3),RGB(3,3)

    DIMENSION CEP0(3,3,MDMIPSF,NDMIPSF)
      DIMENSION CRA(MDMF,NDMF,3,3),CFA(MDMF,NDMF,3,3)
      DIMENSION CRB(MDMIPSF,NDMIPSF,3,3)
      DIMENSION CFB(MDMF,NDMF,3,3,-MDMXF:MDMXF,-NDMXF:NDMXF)
      DIMENSION CRC(MDMF,NDMF,3,3),CFC(MDMF,NDMF,3,3)
      DIMENSION CSUMA(3,3,MDMF,MDM,NDMF,NDM),CS(3,3,MDMF,NDMF)

      CW=2.D0*PI*CRF*1.D6
      CWC2=CW**2/VC**2

      NR=NR1 !+50
      XRHI=1d0
      XRHL=1d0
      IF(NR.EQ.1) THEN
         XRI=1.D6/XRHO(2)
         XRL=XRHO(2)/1.D6
          XRHI=2d0/(2D0*XRHO(NR)+XRHO(NR)-XRHO(NR+1))
          XRHL=(2D0*XRHO(NR)+XRHO(NR)-XRHO(NR+1))/2d0
      ELSE
         XRI=1.D0/XRHO(NR)
         XRL=XRHO(NR)
         IF(NR .GT. 1 )THEN
           XRHI=2d0/(XRHO(NR)+XRHO(NR-1))
           XRHL=(XRHO(NR)+XRHO(NR-1))/2d0
         ENDIF
      ENDIF
      RGH11=0d0
      RGH12=0d0
      RGH13=0d0
      RGH22=0d0
      RGH23=0d0
      RGH33=0d0
      RJH =0d0

!     ----- Fourier decompose metric tensor g -----

      IF (NR .GT. 1)THEN
      DO ND=1,NDMF
      DO MD=1,MDMF
      RGH11(MD,ND,NR)=(RG11(MD,ND,NR) + RG11(MD,ND,NR-1))/2d0
      RGH12(MD,ND,NR)=(RG12(MD,ND,NR) + RG12(MD,ND,NR-1))/2d0
      RGH13(MD,ND,NR)=(RG13(MD,ND,NR) + RG13(MD,ND,NR-1))/2d0
      RGH22(MD,ND,NR)=(RG22(MD,ND,NR) + RG22(MD,ND,NR-1))/2d0
      RGH23(MD,ND,NR)=(RG23(MD,ND,NR) + RG23(MD,ND,NR-1))/2d0
      RGH33(MD,ND,NR)=(RG33(MD,ND,NR) + RG33(MD,ND,NR-1))/2d0
      RJH  (MD,ND,NR)=  (RJ(MD,ND,NR) + RJ(MD,ND,NR-1))/2d0
      ENDDO
      ENDDO
      ELSE
      DO ND=1,NDMF
      DO MD=1,MDMF
        RGH11(MD,ND,NR)=(3D0*RG11(MD,ND,NR) - RG11(MD,ND,NR+1))/2d0
        RGH12(MD,ND,NR)=(3D0*RG12(MD,ND,NR) - RG12(MD,ND,NR+1))/2d0
        RGH13(MD,ND,NR)=(3D0*RG13(MD,ND,NR) - RG13(MD,ND,NR+1))/2d0
        RGH22(MD,ND,NR)=(3D0*RG22(MD,ND,NR) - RG22(MD,ND,NR+1))/2d0
        RGH23(MD,ND,NR)=(3D0*RG23(MD,ND,NR) - RG23(MD,ND,NR+1))/2d0
        RGH33(MD,ND,NR)=(3D0*RG33(MD,ND,NR) - RG33(MD,ND,NR+1))/2d0
        RJH  (MD,ND,NR)=  (3D0*RJ(MD,ND,NR) -   RJ(MD,ND,NR+1))/2d0
      ENDDO
      ENDDO
      ENDIF

      CALL WMSUBG_F(RG11(1,1,NR),RJ(1,1,NR),CGF11(1,1,3))
      CALL WMSUBG_F(RG12(1,1,NR),RJ(1,1,NR),CGF12(1,1,3))
      CALL WMSUBG_F(RG13(1,1,NR),RJ(1,1,NR),CGF13(1,1,3))
      CALL WMSUBG_F(RG22(1,1,NR),RJ(1,1,NR),CGF22(1,1,3))
      CALL WMSUBG_F(RG23(1,1,NR),RJ(1,1,NR),CGF23(1,1,3))
      CALL WMSUBG_F(RG33(1,1,NR),RJ(1,1,NR),CGF33(1,1,3))

      CALL WMSUBG_F(RGH11(1,1,NR),RJH(1,1,NR),CGHF11(1,1,3))
      CALL WMSUBG_F(RGH12(1,1,NR),RJH(1,1,NR),CGHF12(1,1,3))
      CALL WMSUBG_F(RGH13(1,1,NR),RJH(1,1,NR),CGHF13(1,1,3))
      CALL WMSUBG_F(RGH22(1,1,NR),RJH(1,1,NR),CGHF22(1,1,3))
      CALL WMSUBG_F(RGH23(1,1,NR),RJH(1,1,NR),CGHF23(1,1,3))
      CALL WMSUBG_F(RGH33(1,1,NR),RJH(1,1,NR),CGHF33(1,1,3))

!     ----- Calculate dielectric tensor -----

      IF(NS0.EQ.0) THEN
         NS1=1
         NS2=NSMAX
      ELSE
         NS1=NS0
         NS2=NS0
      ENDIF

      DO ND=NDMIN_F,NDMAX_F
      DO MD=MDMIN_F,MDMAX_F
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            DO J=1,3
            DO I=1,3
              CEP0_1 (I,J,NTH,NHH)=0.D0
              CEPH0_1(I,J,NTH,NHH)=0.D0
            ENDDO
            ENDDO
            IF(NS0.EQ.0) THEN
              CEP0_1 (1,1,NTH,NHH)=1.D0
              CEP0_1 (2,2,NTH,NHH)=1.D0
              CEP0_1 (3,3,NTH,NHH)=1.D0
              CEPH0_1(1,1,NTH,NHH)=1.D0
              CEPH0_1(2,2,NTH,NHH)=1.D0
              CEPH0_1(3,3,NTH,NHH)=1.D0
            ENDIF
         ENDDO
         ENDDO
         DO NS=NS1,NS2
            CALL WMTNSR_IPS(NR,NS,MD,ND)

            DO NHH=1,NHHMAX_IPS_F
            DO NTH=1,NTHMAX_IPS_F
               DO J=1,3
               DO I=1,3
                    CEP0_1(I,J,NTH,NHH)
     &             =CEP0_1(I,J,NTH,NHH)
     &             +CTNSR_1(I,J,NTH,NHH)

                    CEPH0_1(I,J,NTH,NHH)
     &             =CEPH0_1(I,J,NTH,NHH)
     &             +CTNSR_1(I,J,NTH,NHH)/2d0
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         IF (NR .EQ. 1 )THEN
            DO NHH=1,NHHMAX_IPS_F
            DO NTH=1,NTHMAX_IPS_F
               DO J=1,3
               DO I=1,3
                    CEPH0_1(I,J,NTH,NHH)
     &             =CEPH0_1(I,J,NTH,NHH)
     &             +CTNSR_1(I,J,NTH,NHH)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDIF

         ENDDO

         IF (NR .GT. 1)THEN
         DO NS=NS1,NS2
            CALL WMTNSR_IPS(NR-1,NS,MD,ND)

             DO NHH=1,NHHMAX_IPS_F
             DO NTH=1,NTHMAX_IPS_F
               DO J=1,3
               DO I=1,3
                    CEPH0_1(I,J,NTH,NHH)
     &             =CEPH0_1(I,J,NTH,NHH)
     &             +CTNSR_1(I,J,NTH,NHH)/2d0
               ENDDO
               ENDDO
             ENDDO
             ENDDO
         ENDDO
         ELSE
         DO NS=NS1,NS2
            CALL WMTNSR_IPS(NR+1,NS,MD,ND)

            DO NHH=1,NHHMAX_IPS_F
            DO NTH=1,NTHMAX_IPS_F
               DO J=1,3
               DO I=1,3
                    CEPH0_1(I,J,NTH,NHH)
     &             =CEPH0_1(I,J,NTH,NHH)
     &             -CTNSR_1(I,J,NTH,NHH)/2d0
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDIF

         IF(NR .EQ. NR_S)THEN
         DO NHH=1,NHHMAX_IPS_F
         DO NTH=1,NTHMAX_IPS_F
           DO J=1,3
           DO I=1,3
             CEP0_1 (I,J,NTH,NHH)=0.D0
           ENDDO
           ENDDO
           IF(NS0.EQ.0) THEN
              CEP0_1 (1,1,NTH,NHH)=1.D0
              CEP0_1 (2,2,NTH,NHH)=1.D0
              CEP0_1 (3,3,NTH,NHH)=1.D0
           ENDIF
        ENDDO
        ENDDO
        ENDIF

        DO NHHF=1,NHHMAX_IPS_F
        DO NTHF=1,NTHMAX_IPS_F
           DO J=1,3
           DO I=1,3
              CRB_1(NTHF,NHHF,I,J)=CEP0_1(I,J,NTHF,NHHF)
              CRHB_1(NTHF,NHHF,I,J)=CEPH0_1(I,J,NTHF,NHHF)
           ENDDO
           ENDDO
        ENDDO
        ENDDO
        DO J=1,3
        DO I=1,3
            CALL WMSUBF_IPS_F(CRB_1(1,1,I,J),CFB_1(1,1,I,J))
            CALL WMSUBF_IPS_F(CRHB_1(1,1,I,J),CFHB_1(1,1,I,J))
        ENDDO
        ENDDO
          KDMAX_TMP=KDMAX_F
          LDMAX_TMP=LDMAX_F
          IF(KDMAX_F==0)KDMAX_TMP=0
          IF(LDMAX_F==0)LDMAX_TMP=0
        DO J=1,3
        DO I=1,3
          DO KD=KDMIN_F,KDMAX_TMP
            KDX_1=KD-KDMIN_IPS_F+1
            KDX=KD-KDMIN_F+1
          DO LD=LDMIN_F,LDMAX_TMP
            LDX_1=LD-LDMIN_IPS_F+1
            LDX=LD-LDMIN_F+1
            CFB(LDX,KDX,I,J,MD,ND) =CFB_1(LDX_1,KDX_1,I,J)
            CFHB(LDX,KDX,I,J,MD,ND)=CFHB_1(LDX_1,KDX_1,I,J)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO

      TCH2=0d0
      TCH3=0d0

      DO NHHF=1,NHHMAX_F
      DO NTHF=1,NTHMAX_F

!        ----- Calculate rotation matrix mu=RMA -----

         CALL WMCMAG_F(NR,NTHF,NHHF,BABS,BSUPTH,BSUPPH)
         TC2=BSUPTH/BABS
         TC3=BSUPPH/BABS

         IF(NR .GT. 1)THEN
         CALL WMCMAG_F(NR-1,NTHF,NHHF,BABSH,BSUPTHH,BSUPPHH)
         TCH2=(BSUPTH+BSUPTHH)/(BABS+BABSH)
         TCH3=(BSUPPH+BSUPPHH)/(BABS+BABSH)
         ELSE
         CALL WMCMAG_F(NR+1,NTHF,NHHF,BABSH,BSUPTHH,BSUPPHH)
         TCH2=(3D0*BSUPTH-BSUPTHH)/(3D0*BABS-BABSH)
         TCH3=(3D0*BSUPPH-BSUPPHH)/(3D0*BABS-BABSH)
         ENDIF

!        ***** RF11=RJ*SQRT(G^11)/XR *****

         RF11=SQRT(RG22(NTHF,NHHF,NR)*RG33(NTHF,NHHF,NR)
     &            -RG23(NTHF,NHHF,NR)*RG23(NTHF,NHHF,NR))
         RMA(1,1)= RJ(NTHF,NHHF,NR)/RF11*XRI
         RMA(2,1)= 0.D0
         RMA(3,1)= 0.D0
         RMA(1,2)= (TC2*(RG23(NTHF,NHHF,NR)*RG12(NTHF,NHHF,NR)
     &                  -RG22(NTHF,NHHF,NR)*RG13(NTHF,NHHF,NR))
     &             +TC3*(RG33(NTHF,NHHF,NR)*RG12(NTHF,NHHF,NR)
     &                  -RG23(NTHF,NHHF,NR)*RG13(NTHF,NHHF,NR))*XRI)
     &             /RF11
         RMA(2,2)= TC3*RF11*XRL
         RMA(3,2)=-TC2*RF11*XRL
         RMA(1,3)=TC2*RG12(NTHF,NHHF,NR)
     &           +TC3*RG13(NTHF,NHHF,NR)*XRI
         RMA(2,3)=TC2*RG22(NTHF,NHHF,NR)*XRL*XRL
     &           +TC3*RG23(NTHF,NHHF,NR)*XRL
         RMA(3,3)=TC2*RG23(NTHF,NHHF,NR)*XRL
     &           +TC3*RG33(NTHF,NHHF,NR)
        
         RFH11=SQRT(RGH22(NTHF,NHHF,NR)*RGH33(NTHF,NHHF,NR)
     &             -RGH23(NTHF,NHHF,NR)*RGH23(NTHF,NHHF,NR))
         RMHA(1,1)= RJH(NTHF,NHHF,NR)/RFH11*XRHI
         RMHA(2,1)= 0.D0
         RMHA(3,1)= 0.D0
         RMHA(1,2)=(TCH2*(RGH23(NTHF,NHHF,NR)*RGH12(NTHF,NHHF,NR)
     &                   -RGH22(NTHF,NHHF,NR)*RGH13(NTHF,NHHF,NR))
     &             +TCH3*(RGH33(NTHF,NHHF,NR)*RGH12(NTHF,NHHF,NR)
     &                   -RGH23(NTHF,NHHF,NR)*RGH13(NTHF,NHHF,NR))*XRHI)
     &              /RFH11
         RMHA(2,2)= TCH3*RFH11*XRHL
         RMHA(3,2)=-TCH2*RFH11*XRHL
         RMHA(1,3)=TCH2*RGH12(NTHF,NHHF,NR)
     &            +TCH3*RGH13(NTHF,NHHF,NR)*XRHI
         RMHA(2,3)=TCH2*RGH22(NTHF,NHHF,NR)*XRHL*XRHL
     &            +TCH3*RGH23(NTHF,NHHF,NR)*XRHL
         RMHA(3,3)=TCH2*RGH23(NTHF,NHHF,NR)*XRHL
     &            +TCH3*RGH33(NTHF,NHHF,NR)

!        ----- Set metric matrix g=RGA -----

         RGA(1,1)=RG11(NTHF,NHHF,NR)
         RGA(1,2)=RG12(NTHF,NHHF,NR)
         RGA(1,3)=RG13(NTHF,NHHF,NR)
         RGA(2,1)=RG12(NTHF,NHHF,NR)
         RGA(2,2)=RG22(NTHF,NHHF,NR)
         RGA(2,3)=RG23(NTHF,NHHF,NR)
         RGA(3,1)=RG13(NTHF,NHHF,NR)
         RGA(3,2)=RG23(NTHF,NHHF,NR)
         RGA(3,3)=RG33(NTHF,NHHF,NR)

         RGHA(1,1)=RGH11(NTHF,NHHF,NR)
         RGHA(1,2)=RGH12(NTHF,NHHF,NR)
         RGHA(1,3)=RGH13(NTHF,NHHF,NR)
         RGHA(2,1)=RGH12(NTHF,NHHF,NR)
         RGHA(2,2)=RGH22(NTHF,NHHF,NR)
         RGHA(2,3)=RGH23(NTHF,NHHF,NR)
         RGHA(3,1)=RGH13(NTHF,NHHF,NR)
         RGHA(3,2)=RGH23(NTHF,NHHF,NR)
         RGHA(3,3)=RGH33(NTHF,NHHF,NR)

!        ----- Invert matrix to obtain mu^(-1)=RMB and g^(-1)=RGB -----

         DO J=1,3
         DO I=1,3
            RMB(I,J)=RMA(I,J)
!            RGB(I,J)=RGA(I,J)

            RMHB(I,J)=RMHA(I,J)
!            RGHB(I,J)=RGHA(I,J)
         ENDDO
         ENDDO

         CALL INVMRD(RMB,3,3,ILL)
         IF(ILL.NE.0) THEN
            WRITE(6,*) 'XX WMSETF: INVMRD(RMB) : SINGULAR MATRIX'
            WRITE(6,'(3I5,1P2E12.4)') NR,NTH,NHH,TC3,TC2
            print *,BSUPTH,BABS
            GOTO 9000
         ENDIF

         CALL INVMRD(RMHB,3,3,ILL)

!         CALL INVMRD(RGB,3,3,ILL)
!         IF(ILL.NE.0) THEN
!            WRITE(6,*) 'XX WMSETF: INVMRD(RGB) : SINGULAR MATRIX'
!            WRITE(6,*) '   NTH,NHH,NR = ',NTH,NHH,NR
!            GOTO 9000
!         ENDIF
!
!         CALL INVMRD(RGHB,3,3,ILL)

!        ----- RGB = g^(-1)
         RGB(1,1)=RGI11(NTHF,NHHF,NR)
         RGB(1,2)=RGI12(NTHF,NHHF,NR)
         RGB(1,3)=RGI13(NTHF,NHHF,NR)
         RGB(2,1)=RGI12(NTHF,NHHF,NR)
         RGB(2,2)=RGI22(NTHF,NHHF,NR)
         RGB(2,3)=RGI23(NTHF,NHHF,NR)
         RGB(3,1)=RGI13(NTHF,NHHF,NR)
         RGB(3,2)=RGI23(NTHF,NHHF,NR)
         RGB(3,3)=RGI33(NTHF,NHHF,NR)

         RGHB(1,1)=RGIH11(NTHF,NHHF,NR)
         RGHB(1,2)=RGIH12(NTHF,NHHF,NR)
         RGHB(1,3)=RGIH13(NTHF,NHHF,NR)
         RGHB(2,1)=RGIH12(NTHF,NHHF,NR)
         RGHB(2,2)=RGIH22(NTHF,NHHF,NR)
         RGHB(2,3)=RGIH23(NTHF,NHHF,NR)
         RGHB(3,1)=RGIH13(NTHF,NHHF,NR)
         RGHB(3,2)=RGIH23(NTHF,NHHF,NR)
         RGHB(3,3)=RGIH33(NTHF,NHHF,NR)


         RGB(1,1)=RGB(1,1)*XRL**2
         RGB(1,2)=RGB(1,2)
         RGB(1,3)=RGB(1,3)*XRL   
         RGB(2,1)=RGB(2,1)
         RGB(2,2)=RGB(2,2)*XRI**2
         RGB(2,3)=RGB(2,3)*XRI   
         RGB(3,1)=RGB(3,1)*XRL   
         RGB(3,2)=RGB(3,2)*XRI   
         RGB(3,3)=RGB(3,3)

         RGHB(1,1)=RGHB(1,1)*XRHL**2
         RGHB(1,2)=RGHB(1,2)
         RGHB(1,3)=RGHB(1,3)*XRHL   
         RGHB(2,1)=RGHB(2,1)
         RGHB(2,2)=RGHB(2,2)*XRHI**2
         RGHB(2,3)=RGHB(2,3)*XRHI   
         RGHB(3,1)=RGHB(3,1)*XRHL   
         RGHB(3,2)=RGHB(3,2)*XRHI   
         RGHB(3,3)=RGHB(3,3)

!        ----- Setup Matrix A=CRA, B=CRB, C=CRC -----

           DO J=1,3
           DO I=1,3
              CSUM=0.D0
              DO K=1,3
                 CSUM=CSUM+RGB(I,K)*RMA(K,J)
              ENDDO
              CRA(NTHF,NHHF,I,J)=CSUM*RJ(NTHF,NHHF,NR)
           ENDDO
           ENDDO

           DO J=1,3
           DO I=1,3
              CSUM=0.D0
              DO K=1,3
                 CSUM=CSUM+RGHB(I,K)*RMHA(K,J)
              ENDDO
              CRHA(NTHF,NHHF,I,J)=CSUM*RJH(NTHF,NHHF,NR)
           ENDDO
           ENDDO

           DO J=1,3
           DO I=1,3
              CRC(NTHF,NHHF,I,J)=RMB(I,J)
              CMC(NTHF,NHHF,I,J)=RMA(I,J)
              CRHC(NTHF,NHHF,I,J)=RMHB(I,J)
              CMHC(NTHF,NHHF,I,J)=RMHA(I,J)
           ENDDO
           ENDDO

      ENDDO
      ENDDO

!     ----- Fourier decompose A, B, C -----

      DO J=1,3
      DO I=1,3
         CALL WMSUBF_F(CRA(1,1,I,J),CFA(1,1,I,J))
         CALL WMSUBF_F(CRHA(1,1,I,J),CFHA(1,1,I,J))
      ENDDO
      ENDDO

      DO J=1,3
      DO I=1,3
         CALL WMSUBF_F(CRC(1,1,I,J),CFC(1,1,I,J))
         CALL WMSUBF_F(CMC(1,1,I,J),CMF(1,1,I,J))
         CALL WMSUBF_F(CRHC(1,1,I,J),CFHC(1,1,I,J))
         CALL WMSUBF_F(CMHC(1,1,I,J),CMHF(1,1,I,J))
         DO NDX=1,NDMF
         DO MDX=1,MDMF
         CMAF(I,J,MDX,NDX,3)=CMF(MDX,NDX,I,J)
         CRMAF(I,J,MDX,NDX,3)=CFC(MDX,NDX,I,J)
         CMAHF(I,J,MDX,NDX,3)=CMHF(MDX,NDX,I,J)
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO KD=KDMIN_F,KDMAX_F
         KDX=KD-KDMIN_F+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO LD=LDMIN_F,LDMAX_F
         LDX=LD-LDMIN_F+1
      DO J=1,3
      DO I=1,3
         CSUMA(I,J,LDX,MDX,KDX,NDX)=0.D0
         CSUMAH(I,J,LDX,MDX,KDX,NDX)=0.D0
         CSUMPF(I,J,LDX,MDX,KDX,NDX)=0.d0
         CSUMPHF(I,J,LDX,MDX,KDX,NDX)=0.d0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=NDMIN,NDMAX
         MDX=MD-MDMIN+1
      DO NDN=1,NDMF
      DO MDN=1,MDMF
            DO J=1,3
            DO I=1,3
        CSUMA(I,J,MDN,MDX,NDN,NDX)=0d0
        CSUMAH(I,J,MDN,MDX,NDN,NDX)=0d0
        CSUMPF(I,J,MDN,MDX,NDN,NDX)=0d0
        CSUMPHF(I,J,MDN,MDX,NDN,NDX)=0d0
            ENDDO
            ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
!$OMP PARALLEL DO default (shared)
!$OMP+private(MDX)
!$OMP+private(KA, KB, KAB,NB1,NB2)
!$OMP+private(KAX,KBX,KABX)
!$OMP+private(LA, LB, LAB,MB1,MB2)
!$OMP+private(LAX,LBX,LABX,K,I,L,J)
!$OMP+private(LBXK,KBXK)
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1

!         DO KAB=KDMIN_F,KDMAX_F
!            KABX=KAB-KDMIN_F+1
!         DO LAB=LDMIN_F,LDMAX_F
!            LABX=LAB-LDMIN_F+1
!         DO K=1,3
!        DO I=1,3
!           CS(I,K,LABX,KABX)=0.D0
!            CSH(I,K,LABX,KABX)=0.D0
!            CPF(I,K,LABX,KABX)=0.D0
!            CPHF(I,K,LABX,KABX)=0.D0
!         ENDDO
!         ENDDO
!         ENDDO
!         ENDDO

         DO KA=KDMIN_F,KDMAX_F
            KAX=KA-KDMIN_F+1
         DO KB=KDMIN_F,KDMAX_F
            KBX=KB-KDMIN_F+1
            KAB=KA+KB
           IF(KAB.GE.KDMIN_F.AND.(KAB.LE.KDMAX_F.OR.KDMAX_F==0))THEN
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
          IF(LAB.GE.LDMIN_F.AND.(LAB.LE.LDMAX_F.OR.LDMAX_F==0))THEN
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
!                     CS(I,K,LABX,KABX)=CS(I,K,LABX,KABX)
        CSUMA(I,K,LABX,MDX,KABX,NDX)=CSUMA(I,K,LABX,MDX,KABX,NDX)
     &                      +CFA(LAX,KAX,I,L)
     &                      *0.25D0*(CFB(LBX,KBX,L,K,MB1,NB1)
     &                              +CFB(LBX,KBX,L,K,MB1,NB2)
     &                              +CFB(LBX,KBX,L,K,MB2,NB1)
     &                              +CFB(LBX,KBX,L,K,MB2,NB2))
!                     CSH(I,K,LABX,KABX)=CSH(I,K,LABX,KABX)
        CSUMAH(I,K,LABX,MDX,KABX,NDX)=CSUMAH(I,K,LABX,MDX,KABX,NDX)
     &                      +CFHA(LAX,KAX,I,L)
     &                      *0.25D0*(CFHB(LBX,KBX,L,K,MB1,NB1)
     &                              +CFHB(LBX,KBX,L,K,MB1,NB2)
     &                              +CFHB(LBX,KBX,L,K,MB2,NB1)
     &                              +CFHB(LBX,KBX,L,K,MB2,NB2))
                  ENDDO
                  ENDDO
                  ENDDO
            ENDIF
         ENDDO
         ENDDO
            ENDIF
         ENDDO
         ENDDO

         KA= 0 
         KAX=KA-KDMIN_F+1
         DO KB=KDMIN_F,KDMAX_F
            KBX=KB-KDMIN_F+1
            KAB=KA+KB

              IF(MOD(KB,2).EQ.0) THEN
                  NB1=ND-(KA+KB/2)
                  NB2=ND-(KA+KB/2)
               ELSE
                  NB1=ND-(KA+(KB-1)/2)
                  NB2=ND-(KA+(KB+1)/2)
               ENDIF
               LA=0
               LAX=LA-LDMIN_F+1
               DO LB=LDMIN_F,LDMAX_F
                  LBX=LB-LDMIN_F+1
                  LAB=LA+LB

                     IF(MOD(LB,2).EQ.0) THEN
                        MB1=MD-(LA+LB/2)
                        MB2=MD-(LA+LB/2)
                     ELSE
                        MB1=MD-(LA+(LB-1)/2)
                        MB2=MD-(LA+(LB+1)/2)
                     ENDIF
                     LBXK=LBX
                     KBXK=KBX
                     DO K=1,3
                     DO L=1,3
!                     CPF(L,K,LBXK,KBXK)=CPF(L,K,LBXK,KBXK)
      CSUMPF(L,K,LBXK,MDX,KBXK,NDX)=CSUMPF(L,K,LBXK,MDX,KBXK,NDX)
     &                       +0.25D0*(CFB(LBX,KBX,L,K,MB1,NB1)
     &                               +CFB(LBX,KBX,L,K,MB1,NB2)
     &                               +CFB(LBX,KBX,L,K,MB2,NB1)
     &                               +CFB(LBX,KBX,L,K,MB2,NB2))
!                     CPHF(L,K,LBXK,KBXK)=CPHF(L,K,LBXK,KBXK)
      CSUMPHF(L,K,LBXK,MDX,KBXK,NDX)=CSUMPHF(L,K,LBXK,MDX,KBXK,NDX)
     &                       +0.25D0*(CFHB(LBX,KBX,L,K,MB1,NB1)
     &                               +CFHB(LBX,KBX,L,K,MB1,NB2)
     &                               +CFHB(LBX,KBX,L,K,MB2,NB1)
     &                               +CFHB(LBX,KBX,L,K,MB2,NB2))
                     ENDDO
                     ENDDO
               ENDDO
         ENDDO

!            DO J=1,3
!            DO I=1,3
!            CSUMA(I,J,1:MDMF,MDX,1:NDMF,NDX)=CS(I,J,1:MDMF,1:NDMF)
!            CSUMAH(I,J,1:MDMF,MDX,1:NDMF,NDX)=CSH(I,J,1:MDMF,1:NDMF)
!            CSUMPF(I,J,1:MDMF,MDX,1:NDMF,NDX)=CPF(I,J,1:MDMF,1:NDMF)
!            CSUMPHF(I,J,1:MDMF,MDX,1:NDMF,NDX)=CPHF(I,J,1:MDMF,1:NDMF)
!            ENDDO
!            ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENDDO


      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO KD=KDMIN_F,KDMAX_F
         KDX=KD-KDMIN_F+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO LD=LDMIN_F,LDMAX_F
         LDX=LD-LDMIN_F+1
      DO J=1,3
      DO I=1,3
         CGD(I,J,LDX,MDX,KDX,NDX,3)=-CWC2*CSUMA(I,J,LDX,MDX,KDX,NDX)
         CGDH(I,J,LDX,MDX,KDX,NDX,3)=-CWC2*CSUMAH(I,J,LDX,MDX,KDX,NDX)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO KD=KDMIN_F,KDMAX_F
         KDX=KD-KDMIN_F+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO LD=LDMIN_F,LDMAX_F
         LDX=LD-LDMIN_F+1
      DO J=1,3
      DO I=1,3
       CPSF(I,J,LDX,MDX,KDX,NDX,3)=-CWC2*CSUMPF(I,J,LDX,MDX,KDX,NDX)
       CPSHF(I,J,LDX,MDX,KDX,NDX,3)=-CWC2*CSUMPHF(I,J,LDX,MDX,KDX,NDX)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN

 9000 CONTINUE
      RETURN
      END

!     ****** 2D FOURIER TRANSFORM ******

      SUBROUTINE WMSUBG(RF1,RF2,CF)

      INCLUDE 'wmcomm.inc'

      DIMENSION RF1(MDM,NDM),RF2(MDM,NDM),CF(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)

      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=RF1(NTH,NHH)/RF2(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         IF (LDSIZ == 1)THEN
            LDX=LDSIZ
            CF(LDX,NHH)=CFM(LDX)
         ELSE
            DO LDX=1,LDSIZ
               CF(LDX,NHH)=CFM(LDX)
            ENDDO
         ENDIF
      ENDDO

      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         IF (KDSIZ == 1)THEN
            KDX=KDSIZ
            CF(LDX,KDX)=CFN(KDX)
         ELSE
            DO KDX=1,KDSIZ
               CF(LDX,KDX)=CFN(KDX)
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END

!     ****** 2D FOURIER TRANSFORM ******

      SUBROUTINE WMSUBG_F(RF1,RF2,CF)

      INCLUDE 'wmcomm.inc'

      DIMENSION RF1(MDMF,NDMF),RF2(MDMF,NDMF),CF(MDMF,NDMF)
      DIMENSION CFM(MDMF),CFN(NDMF)

      DO NHH=1,NHHMAX_F
         DO NTH=1,NTHMAX_F
            CFM(NTH)=RF1(NTH,NHH)/RF2(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX_F,0)
         IF (LDSIZ_F == 1)THEN
            LDX=LDSIZ_F
            CF(LDX,NHH)=CFM(LDX)
         ELSE
            DO LDX=1,LDSIZ_F
               CF(LDX,NHH)=CFM(LDX)
            ENDDO
         ENDIF
      ENDDO

      DO LDX=1,LDSIZ_F
         DO NHH=1,NHHMAX_F
            CFN(NHH)=CF(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX_F,0)
         IF (KDSIZ_F == 1)THEN
            KDX=KDSIZ_F
            CF(LDX,KDX)=CFN(KDX)
         ELSE
            DO KDX=1,KDSIZ_F
               CF(LDX,KDX)=CFN(KDX)
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END

!     ****** 2D FOURIER TRANSFORM ******

      SUBROUTINE WMSUBF(CF1,CF2)

      INCLUDE 'wmcomm.inc'

      DIMENSION CF1(MDM,NDM),CF2(MDM,NDM)
      DIMENSION CFM(MDM),CFN(NDM)

      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX,0)
         IF (LDSIZ == 1)THEN
            LDX=LDSIZ
            CF2(LDX,NHH)=CFM(LDX)
         ELSE
            DO LDX=1,LDSIZ
               CF2(LDX,NHH)=CFM(LDX)
            ENDDO
         ENDIF
      ENDDO

      DO LDX=1,LDSIZ
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX,0)
         IF (KDSIZ == 1)THEN
            KDX=KDSIZ
            CF2(LDX,KDX)=CFN(KDX)
         ELSE
            DO KDX=1,KDSIZ
               CF2(LDX,KDX)=CFN(KDX)
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END

!     ****** 2D FOURIER TRANSFORM 2 ******

      SUBROUTINE WMSUBF_F(CF1,CF2)

      INCLUDE 'wmcomm.inc'

      DIMENSION CF1(MDMF,NDMF),CF2(MDMF,NDMF)
      DIMENSION CFM(MDMF),CFN(NDMF)

      DO NHH=1,NHHMAX_F
         DO NTH=1,NTHMAX_F
            CFM(NTH)=CF1(NTH,NHH)
         ENDDO
         CALL WMXFFT(CFM,NTHMAX_F,0)
         IF (LDSIZ_F == 1)THEN
            LDX=LDSIZ_F
            CF2(LDX,NHH)=CFM(LDX)
         ELSE
            DO LDX=1,LDSIZ_F
               CF2(LDX,NHH)=CFM(LDX)
            ENDDO
         ENDIF
      ENDDO

      DO LDX=1,LDSIZ_F
         DO NHH=1,NHHMAX_F
            CFN(NHH)=CF2(LDX,NHH)
         ENDDO
         CALL WMXFFT(CFN,NHHMAX_F,0)
         IF (KDSIZ_F == 1)THEN
            KDX=KDSIZ_F
            CF2(LDX,KDX)=CFN(KDX)
         ELSE
            DO KDX=1,KDSIZ_F
               CF2(LDX,KDX)=CFN(KDX)
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END

!     ****** INTERFACE FOR FFT ******

      SUBROUTINE WMXFFT(CA,N,KEY)

      USE libfft,ONLY: FFT2L
      INCLUDE 'wmcomm.inc'

      COMPLEX*16 CA(N)
      DATA NS/0/
!
      IF(N.NE.1) THEN
         IF(N.EQ.NS) THEN
            IND=0
         ELSE
            IND=1
            NS=N
         END IF
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

      RETURN
      END
