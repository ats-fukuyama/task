C     $Id$

!---- interface for wm parameter

      subroutine get_wmparm(crf_,nth0_,nph0_,idbgwm_)
      
      include '../wm/wmcomm.inc'
      complex(8),intent(out):: crf_
      integer,intent(out):: nth0_,nph0_
      integer,intent(out):: idbgwm_
      crf_=crf
      nth0_=nth0
      nph0_=nph0
      idbgwm_=idbgwm
      return
      end subroutine get_wmparm

!---- interface for wm parameter

      subroutine get_wmparm1(rr_,ra_,rb_)
      
      include '../wm/wmcomm.inc'
      real(8),intent(out):: rr_,ra_,rb_
      rr_=rr
      ra_=ra
      rb_=rb
      return
      end subroutine get_wmparm1
C
C     ****** CALCULATE METRIC AND CONVERSION TENSOR ******
C
      SUBROUTINE wmfem_metric(gma,mma,gj)
C
      INCLUDE 'wmcomm.inc'
      real(8),intent(out):: gma(3,3,nthmax,nphmax,nrmax+1)
      real(8),intent(out):: mma(3,3,nthmax,nphmax,nrmax+1)
      real(8),intent(out):: gj(nthmax,nphmax,nrmax+1)
      real(8),dimension(3,3):: RMA
C
      DO NR=1,NRMAX+1
C
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
C
            gma(1,1,nth,nph,nr)=RG11(nth,nph,nr)
            gma(1,2,nth,nph,nr)=RG12(nth,nph,nr)
            gma(1,3,nth,nph,nr)=RG13(nth,nph,nr)
            gma(2,1,nth,nph,nr)=RG12(nth,nph,nr)
            gma(2,2,nth,nph,nr)=RG22(nth,nph,nr)
            gma(2,3,nth,nph,nr)=RG23(nth,nph,nr)
            gma(3,1,nth,nph,nr)=RG13(nth,nph,nr)
            gma(3,2,nth,nph,nr)=RG23(nth,nph,nr)
            gma(3,3,nth,nph,nr)=RG33(nth,nph,nr)
            gj(nth,nph,nr)=RJ(nth,nph,nr)
C
C        ----- Calculate rotation matrix mu=RMA -----
C
            BSUPTH=BFLD(2,NTH,NPH,NR)
            BSUPPH=BFLD(3,NTH,NPH,NR)
            BABS=SQRT(     RG22(NTH,NPH,NR)*BSUPTH*BSUPTH
     &               +2.D0*RG23(NTH,NPH,NR)*BSUPTH*BSUPPH
     &               +     RG33(NTH,NPH,NR)*BSUPPH*BSUPPH)
            TC2=BSUPTH/BABS
            TC3=BSUPPH/BABS
C
C        ***** RF11=RJ*SQRT(G^11) *****
C
            if(nr.eq.1) then
               NRL=3
               RJb  =RJ(NTH,NPH,NRL)
               RF11b=SQRT(RG22(NTH,NPH,NRL)*RG33(NTH,NPH,NRL)
     &                   -RG23(NTH,NPH,NRL)*RG23(NTH,NPH,NRL))
               RMAb = (TC2*(RG23(NTH,NPH,NRL)*RG12(NTH,NPH,NRL)
     &                     -RG22(NTH,NPH,NRL)*RG13(NTH,NPH,NRL))
     &                +TC3*(RG33(NTH,NPH,NRL)*RG12(NTH,NPH,NRL)
     &                     -RG23(NTH,NPH,NRL)*RG13(NTH,NPH,NRL)))
               NRL=2
               RJa  =RJ(NTH,NPH,NRL)
               RF11a=SQRT(RG22(NTH,NPH,NRL)*RG33(NTH,NPH,NRL)
     &                   -RG23(NTH,NPH,NRL)*RG23(NTH,NPH,NRL))
               RMAa = (TC2*(RG23(NTH,NPH,NRL)*RG12(NTH,NPH,NRL)
     &                     -RG22(NTH,NPH,NRL)*RG13(NTH,NPH,NRL))
     &                +TC3*(RG33(NTH,NPH,NRL)*RG12(NTH,NPH,NRL)
     &                     -RG23(NTH,NPH,NRL)*RG13(NTH,NPH,NRL)))

               RJL  =(RJa*xrho(3)**2  -RJb*xrho(2)**2)
     &               /(xrho(3)*xrho(2)*(xrho(3)-xrho(2)))
               RF11L=(RF11a*xrho(3)**2-RF11b*xrho(2)**2)
     &               /(xrho(3)*xrho(2)*(xrho(3)-xrho(2)))
               RMAL =(RMAa*xrho(3)**2 -RMAb*xrho(2)**2)
     &               /(xrho(3)*xrho(2)*(xrho(3)-xrho(2)))
            else
               RJL  =RJ(NTH,NPH,NR)
               RF11L=SQRT(RG22(NTH,NPH,NR)*RG33(NTH,NPH,NR)
     &                  -RG23(NTH,NPH,NR)*RG23(NTH,NPH,NR))
               RMAL = (TC2*(RG23(NTH,NPH,NR)*RG12(NTH,NPH,NR)
     &                     -RG22(NTH,NPH,NR)*RG13(NTH,NPH,NR))
     &                +TC3*(RG33(NTH,NPH,NR)*RG12(NTH,NPH,NR)
     &                     -RG23(NTH,NPH,NR)*RG13(NTH,NPH,NR)))
            endif
            
            RMA(1,1)= RJL/RF11L
            RMA(2,1)= 0.D0
            RMA(3,1)= 0.D0
            RMA(1,2)= RMAL/RF11L
            RMA(2,2)= TC3*RF11
            RMA(3,2)=-TC2*RF11
            RMA(1,3)=TC2*RG12(NTH,NPH,NR)
     &              +TC3*RG13(NTH,NPH,NR)
            RMA(2,3)=TC2*RG22(NTH,NPH,NR)
     &              +TC3*RG23(NTH,NPH,NR)
            RMA(3,3)=TC2*RG23(NTH,NPH,NR)
     &              +TC3*RG33(NTH,NPH,NR)

            do j=1,3
               do i=1,3
                  mma(i,j,nth,nph,nr)=RMA(i,j)
               enddo
            enddo

         enddo
         enddo
            write(6,*) 'gj(1,1,',NR,'1)=',gj(1,1,NR)
            write(6,*) 'mma(1,1,1,1,',NR,')=',mma(1,1,1,1,NR)
         enddo
         return
         end
      
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE wmfem_disp(nr,ns,fms)
C
C           NR : NODE NUMBER (RADIAL POSITION)
C           NS : PARTICLE SPECIES 
C
      INCLUDE 'wmcomm.inc'
      complex(8):: fms(3,3,nthmax*nphmax,nthmax*nphmax)
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
      ELSEIF(NS.EQ.5.OR.NS.EQ.6) THEN
         CALL WMTNDK(NR,NS)
      ELSE
         IF(MODELP(NS).LT.0) THEN
            CALL WMTNSX(NR,NS)
         ELSEIF(MODELP(NS).EQ.8) THEN
            IF(MODELV(NS).EQ.9) THEN
               MODELVS=MODELV(NS)
               MODELV(NS)=MODELVR(NR,NS)
               CALL WMDPIN(NR,NS)
               MODELV(NS)=MODELVS
            ELSE
               CALL WMDPIN(NR,NS)
            ENDIF
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
            CALL WMDPIN(NR,NS)
         ENDIF
      ENDIF
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         nfc1=nthmax*(nph-1)+nth
         DO ND=NDMIN,NDMAX
            ndx=nd-ndmin+1
         DO MD=MDMIN,MDMAX
            mdx=md-mdmin+1
            nfc2=nthmax*(ndx-1)+mdx
            do j=1,3
            do i=1,3
               fms(i,j,nfc1,nfc2)=CTNSR(i,j,MD,ND,NTH,NPH)
            enddo
            enddo
         ENDDO
         ENDDO
      ENDDO
      ENDDO
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
C     ****** ASSEMBLE TOTAL ELEMENT FREE VECTOR ******
C
      SUBROUTINE get_wmfvb(NR,CFVP)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CFVP(nphmax,nthmax,3)
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
C
         DO MDX=1,MDSIZ
            DO NDX=1,NDSIZ
               CFVP(NDX,MDX,1)=0.D0
               CFVP(NDX,MDX,2)=0.D0
               CFVP(NDX,MDX,3)=0.D0
            ENDDO
         ENDDO
C
         IF(MODEEG.EQ.0) THEN
            DO NRI=1,NRMAX 
               IF(XR(NRI)/RD.LT.1.D0) NRANT=NRI
            ENDDO
C
            IF(NR+1.GE.NRANT) THEN
C
C               write(6,*) 'nr,nrant=',nr,nrant
C
               CW=2*PI*CRF*1.D6
               CC=CI*CW*RMU0
               DPH=2.D0*PI/NPHMAX
               DTH=2.D0*PI/NTHMAX
C
               IF(NR+1.EQ.NRANT.OR.NR+1.EQ.NRANT+1) THEN
                  XRHO1=XRHO(NRANT)
                  XRHO2=XRHO(NRANT+1)
                  DRHO=XRHO2-XRHO1
                  XRHOC=0.5D0*(XRHO2+XRHO1)
C
                  IF(MODELG.EQ.3) THEN
                     QPC=0.5D0*(QPS(NRANT)+QPS(NRANT+1))
                     DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
                  ELSE
                     DPSIPDRHOC =2.D0*PSIPA*XRHOC
                  ENDIF
C
                  FACTM=(XRHO2-RD/RA)/DRHO
                  FACTP=(RD/RA-XRHO1)/DRHO
C
                  DO ND=NDMIN,NDMAX
                     NDX=ND-NDMIN+1
                     NN=NPH0+NHC*ND
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     MM=NTH0+MD
                     CJTHM=CC*CJANT(2,MDX,NDX)*FACTM*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*FACTM*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)
     &                                         *DPSIPDRHOC*DRHO/XRHOC
                     CJTHP=CC*CJANT(2,MDX,NDX)*FACTP*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DPH)
                     CJPHP=CC*CJANT(3,MDX,NDX)*FACTP*XRHOC
     &                                        /(DPSIPDRHOC*DRHO*DTH)
                     CFVP(NDX,MDX,1)=0.D0
                     CFVP(NDX,MDX,2)=0.D0
                     CFVP(NDX,MDX,3)=0.D0
                     IF(NR+1.EQ.NRANT) THEN
                        CFVP(NDX,MDX,2)=CJTHM
                        CFVP(NDX,MDX,3)=CJPHM
                     ELSE IF(NR+1.EQ.NRANT+1) THEN
                        CFVP(NDX,MDX,1)=CJR
                        CFVP(NDX,MDX,2)=CJTHP
                        CFVP(NDX,MDX,3)=CJPHP
                     ENDIF
c$$$                     WRITE(6,'(A,3I5,A,1P2E12.4)') 
c$$$     &                    'CFVP(',NDX,MDX,1,')',CFVP(NDX,MDX,1)
c$$$                     WRITE(6,'(A,3I5,A,1P2E12.4)')
c$$$     &                    'CFVP(',NDX,MDX,2,')',CFVP(NDX,MDX,2)
c$$$                     WRITE(6,'(A,3I5,A,1P2E12.4)')
c$$$     &                    'CFVP(',NDX,MDX,3,')',CFVP(NDX,MDX,3)
                  ENDDO
                  ENDDO
               ELSE
                  XRHO1=XRHO(NR+1)
                  IF(NR.LT.NRMAX) THEN
                     XRHO2=XRHO(NR+2)
                  ELSE
                     XRHO2=2*XRHO(NR+1)-XRHO(NR)
                  ENDIF
                  DRHO=XRHO2-XRHO1
                  XRHOC=0.5D0*(XRHO2+XRHO1)
C
                  IF(MODELG.EQ.3) THEN
                     IF(NR.LT.NRMAX) THEN
                        QPC=0.5D0*(QPS(NR+1)+QPS(NR+2))
                     ELSE
                        QPC=0.5D0*(3*QPS(NR+1)-QPS(NR))
                     ENDIF
                     DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
                  ELSE
                     DPSIPDRHOC =2.D0*PSIPA*XRHOC
                  ENDIF
C
                  DO ND=NDMIN,NDMAX
                     NDX=ND-NDMIN+1
                     NN=NPH0+NHC*ND
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     MM=NTH0+MD
                     CJTHM=CC*CJANT(2,MDX,NDX)*XRHOC
     &                             /(DPSIPDRHOC*DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*XRHOC
     &                             /(DPSIPDRHOC*DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DPSIPDRHOC*DRHO
     &                             /XRHOC
                     CFVP(NDX,MDX,1)=CJR
c$$$                     WRITE(6,'(A,3I5,A,1P2E12.4)') 
c$$$     &                    'CFVP(',NDX,MDX,1,')',CFVP(NDX,MDX,1)
c$$$                     WRITE(6,'(A,3I5,A,1P2E12.4)')
c$$$     &                    'CFVP(',NDX,MDX,2,')',CFVP(NDX,MDX,2)
c$$$                     WRITE(6,'(A,3I5,A,1P2E12.4)')
c$$$     &                    'CFVP(',NDX,MDX,3,')',CFVP(NDX,MDX,3)
                  ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ELSE
            DPH=2.D0*PI/NPHMAX
            DTH=2.D0*PI/NTHMAX
            XRHO1=XRHO(NR+1)
            IF(NR.LT.NRMAX) THEN
               XRHO2=XRHO(NR+2)
            ELSE
               XRHO2=2*XRHO(NR+1)-XRHO(NR)
            ENDIF
            DRHO=XRHO2-XRHO1
            XRHOC=0.5D0*(XRHO2+XRHO1)
C
            IF(MODELG.EQ.3) THEN
               IF(NR.LT.NRMAX) THEN
                  QPC=0.5D0*(QPS(NR+1)+QPS(NR+2))
               ELSE
                  QPC=0.5D0*(3*QPS(NR+1)-QPS(NR))
               ENDIF
               DPSIPDRHOC =2.D0*PSITA*XRHOC /QPC
            ELSE
               DPSIPDRHOC =2.D0*PSIPA*XRHOC
            ENDIF
C
            CALL WMCDEN(NR+1,RN,RTPR,RTPP,RU)
            RT=(RTPR(1)+2*RTPP(1))
            RJFACT=RN(1)*RT
            RJFACT=RJFACT*XRHO(NR+1)
            DO ND=NDMIN,NDMAX
               NDX=ND-NDMIN+1
               NN=NPH0+NHC*ND
            DO MD=MDMIN,MDMAX
               MDX=MD-MDMIN+1
               MM=NTH0+MD
               CJTHM=RJFACT*XRHOC/(DPSIPDRHOC*DRHO*DPH)
               CJPHM=RJFACT*XRHOC/(DPSIPDRHOC*DRHO*DTH)
               CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DPSIPDRHOC*DRHO/XRHOC
C               CFVP(NDX,MDX,1)=CJR
               CFVP(NDX,MDX,1)=0.D0
               CFVP(NDX,MDX,2)=CJTHM
               CFVP(NDX,MDX,3)=CJPHM
            ENDDO
            ENDDO
         ENDIF
C
      RETURN
      END
C
C     ****** CALCULATE ELECTRIC FIELD ******
C
      SUBROUTINE WMFEM_EFLD(cef)
C
      INCLUDE 'wmcomm.inc'
      dimension cef(3,nthmax,nphmax,nrmax+1)
C
      DIMENSION CEF1(MDM,NDM),CEF2(MDM,NDM),RMA(3,3)
C
      do nr=1,nrmax+1
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDK(1,MDX,NDX,NR )=cef(1,mdx,ndx,nr)
            CEFLDK(2,MDX,NDX,NR )=cef(2,mdx,ndx,nr)
            CEFLDK(3,MDX,NDX,NR )=cef(3,mdx,ndx,nr)
            CEFLD(1,MDX,NDX,NR )=cef(1,mdx,ndx,nr)
            CEFLD(2,MDX,NDX,NR )=cef(2,mdx,ndx,nr)
            CEFLD(3,MDX,NDX,NR )=cef(3,mdx,ndx,nr)
         enddo
         enddo
      enddo
C
C     ----- Inverse Fourier transform ----
C
      DO NR=1,NRMAX+1
      DO I=1,3
         DO NDX=1,NDSIZ
         DO MDX=1,MDSIZ
            CEF1(MDX,NDX)=CEFLD(I,MDX,NDX,NR)
         ENDDO
         ENDDO
         CALL WMSUBEX(CEF1,CEF2,NTHMAX,NPHMAX)
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CEFLD(I,NTH,NPH,NR)=CEF2(NTH,NPH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- Calculate CEN from CEsup -----
C
C
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         CEN(1,NTH,NPH,NR)=CEFLD(1,NTH,NPH,NR)
         CEN(2,NTH,NPH,NR)=CEFLD(1,NTH,NPH,NR)
         CEN(3,NTH,NPH,NR)=CEFLD(1,NTH,NPH,NR)
         CEP(1,NTH,NPH,NR)=(   CEN(1,NTH,NPH,NR)
     &                     +CI*CEN(2,NTH,NPH,NR))/SQRT(2.D0)
         CEP(2,NTH,NPH,NR)=(   CEN(1,NTH,NPH,NR)
     &                     -CI*CEN(2,NTH,NPH,NR))/SQRT(2.D0)
         CEP(3,NTH,NPH,NR)=    CEN(3,NTH,NPH,NR)
      ENDDO
      ENDDO
      ENDDO
c$$$C
c$$$      DO K=1,3
c$$$         DO NR=1,NRMAX+1
c$$$            IF(XRHO(NR).GT.1.0D0) THEN
c$$$               DO NTH=1,NTHMAX
c$$$                  DO NPH=1,NPHMAX
c$$$                     CEP(K,NTH,NPH,NR)=(0.D0,0.D0)
c$$$                  ENDDO
c$$$               ENDDO
c$$$            ELSE
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDDO
C
C
      RETURN
      END
C
C     ****** CALCULATE ABSORBED POWER ******
C
      SUBROUTINE WMFEM_PABS(cpp)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION DS(NRM),DSS(NTHM,NPHM,NRM)
      DIMENSION CPF1(MDM,NDM),CPF2(MDM,NDM)
      dimension cpp(nthmax,nphmax,nthmax,nphmax,nrmax+1,nsmax)
C
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NPHMAX)
C
      DO NR=1,nrmax
      DO NS=1,NSMAX
      DO NKX=1,NDSIZ
      DO KDX=1,KDSIZ
      DO MLX=1,MDSIZ
      DO LDX=1,LDSIZ
         CPABS(LDX,MLX,KDX,NKX,NS,NR)=cpp(MLX,NKX,LDX,KDX,NR,NS)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
C
C     +++++ CALCULATE PABS IN MODE NUMBER SPACE +++++
C
         DO NR=1,NRMAX
         DO NS=1,NSMAX
         DO NDX=1,NDSIZ
            KK=0
            KKX=KK-KDMIN+1
         DO MDX=1,MDSIZ
            LL=0
            LLX=LL-LDMIN+1
            PABSK(MDX,NDX,NR,NS)=DBLE(CPABS(LLX,MDX,KKX,NDX,NS,NR))
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE PABS IN REAL SPACE +++++
C
         DO NS=1,NSMAX
         DO NR=1,NRMAX
            DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NPH,NR,NS)=0.D0
            ENDDO
            ENDDO
            DO NDX=1,NDSIZ
            DO MDX=1,MDSIZ
               DO KK=KDMIN,KDMAX
                  KKX=KK-KDMIN+1
               DO LL=LDMIN,LDMAX
                  LLX=LL-LDMIN+1
                  CPF1(LLX,KKX)=CPABS(LLX,MDX,KKX,NDX,NS,NR)
               ENDDO
               ENDDO
               CALL WMSUBEX(CPF1,CPF2,NTHMAX,NPHMAX)
               DO NPH=1,NPHMAX
               DO NTH=1,NTHMAX
                  PABS(NTH,NPH,NR,NS)=PABS(NTH,NPH,NR,NS)
     &                               +DBLE(CPF2(NTH,NPH))
C                  write(6,'(3I8,1P3E12.4)') NR,NTH,NPH,CPF2(NTH,NPH),
C     &                 PABS(NTH,NPH,NR,NS)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
C     +++++ CALCULATE DRIVEN CURRENT IN REAL SPACE +++++
C
      NS=1
      DO NR=1,NRMAX
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NPH,NR)=0.D0
         ENDDO
         ENDDO
         CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
         VTE=SQRT(RTPR(1)*AEE*1.D3/(PA(1)*AMP))
         WW=DBLE(CW)
         IF(RN(1).LE.0.D0) THEN
            RLNLMD=15.D0
         ELSE
            RT=(RTPR(1)+2*RTPP(1))/3.D0
            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(1))
     &           + 2.30D0*LOG10(RT)
         ENDIF
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NN=NPH0+ND
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            MM=NTH0+MD
            DO KKX=1,KDSIZ
            DO LLX=1,LDSIZ
               CPF1(LLX,KKX)=CPABS(LLX,MDX,KKX,NDX,NS,NR)
            ENDDO
            ENDDO
            CALL WMSUBEX(CPF1,CPF2,NTHMAX,NPHMAX)
            DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
               RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
              IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
              IF(ABS(WW/RKPR).LT.VC) THEN
                 W=WW/(RKPR*VTE)
                 XL=(RPST(NTH,NPH,NR)-RR  )/RR
                 YL=(ZPST(NTH,NPH,NR)-0.D0)/RR
                 EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
                 IF(W.LT.0.D0) EFCD=-EFCD
                 IF (RN(1).GT.0.D0) THEN
                    PCUR(NTH,NPH,NR)=PCUR(NTH,NPH,NR)
     &                   +0.384D0*RTPR(1)/(AEE*1.D3)*EFCD
     &                   /((RN(1)/1.D20)*RLNLMD)*DBLE(CPF2(NTH,NPH))
     &                   /(2.D0*PI*RPST(NTH,NPH,NR))
                 END IF
              ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PCURR(NR)=PCURR(NR)+PCUR(NTH,NPH,NR)*DTH
         ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         PABSR(NR,NS)=0.D0
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NPH,NR,NS)*DTH*DPH
         ENDDO
         ENDDO
C         write(6,'(2I8,1P3E12.4)') NR,NS,PABSR(NR,NS),DTH,DPH
      ENDDO
      ENDDO
C
      PCURT=0.D0
      DO NR=1,NRMAX
         PCURT=PCURT+PCURR(NR)
      ENDDO
C
      DO NS=1,NSMAX
         PABST(NS)=0.D0
         DO NR=1,NRMAX
            PABST(NS)=PABST(NS)+PABSR(NR,NS)
         ENDDO
      ENDDO
C
      PABSTT=0.D0
      DO NS=1,NSMAX
         PABSTT=PABSTT+PABST(NS)
      ENDDO
C
      FACT=1.D0
C
      IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
         FACT=PRFIN/PABSTT
      ENDIF
      FACTSQ=SQRT(FACT)
C
      NR=1
         DS(NR)=0.D0
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            DSS(NTH,NPH,NR)=0.D0
         ENDDO
         ENDDO
      DO NR=2,NRMAX
         DS(NR)=0.D0
         DRHO=0.5D0*(XRHO(NR+1)-XRHO(NR-1))
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            IF(MODELG.EQ.3) THEN
               DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
            ELSE
               DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
            ENDIF
            DSSS=DPSIPDRHO*DRHO*RJ(NTH,NPH,NR)
            DSS(NTH,NPH,NR)=1.D0/DSSS
            DS(NR)=DS(NR)+DSSS*DTH*DPH
         ENDDO
         ENDDO
         DS(NR)=1.D0/DS(NR)
      ENDDO
C
      PABSTT=FACT*PABSTT
      DO NS=1,NSMAX
         PABST(NS)=FACT*PABST(NS)
         DO NR=1,NRMAX
            PABSR(NR,NS)=FACT*PABSR(NR,NS)*DS(NR)
            DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NPH,NR,NS)=FACT*PABS(NTH,NPH,NR,NS)
     &                            *DSS(NTH,NPH,NR)
            ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO NS=1,NSMAX
      DO NR=1,NRMAX
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            PABSK(MDX,NDX,NR,NS)=FACT*PABSK(MDX,NDX,NR,NS)*DS(NR)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      PCURT=FACT*PCURT
      DO NR=1,NRMAX
         PCURR(NR)=FACT*PCURR(NR)*DS(NR)
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            PCUR(NTH,NPH,NR)=FACT*PCUR(NTH,NPH,NR)*DSS(NTH,NPH,NR)
         ENDDO
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX+1
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
         DO MD=MDMIN,MDMAX
            MDX=MD-MDMIN+1
            CEFLDK(1,MDX,NDX,NR)=FACTSQ*CEFLDK(1,MDX,NDX,NR)
            CEFLDK(2,MDX,NDX,NR)=FACTSQ*CEFLDK(2,MDX,NDX,NR)
            CEFLDK(3,MDX,NDX,NR)=FACTSQ*CEFLDK(3,MDX,NDX,NR)
            CBFLDK(1,MDX,NDX,NR)=FACTSQ*CBFLDK(1,MDX,NDX,NR)
            CBFLDK(2,MDX,NDX,NR)=FACTSQ*CBFLDK(2,MDX,NDX,NR)
            CBFLDK(3,MDX,NDX,NR)=FACTSQ*CBFLDK(3,MDX,NDX,NR)
         ENDDO
         ENDDO
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CEFLD(1,NTH,NPH,NR)=FACTSQ*CEFLD(1,NTH,NPH,NR)
            CEFLD(2,NTH,NPH,NR)=FACTSQ*CEFLD(2,NTH,NPH,NR)
            CEFLD(3,NTH,NPH,NR)=FACTSQ*CEFLD(3,NTH,NPH,NR)
            CBFLD(1,NTH,NPH,NR)=FACTSQ*CBFLD(1,NTH,NPH,NR)
            CBFLD(2,NTH,NPH,NR)=FACTSQ*CBFLD(2,NTH,NPH,NR)
            CBFLD(3,NTH,NPH,NR)=FACTSQ*CBFLD(3,NTH,NPH,NR)
            CEN(1,NTH,NPH,NR)  =FACTSQ*CEN(1,NTH,NPH,NR)
            CEN(2,NTH,NPH,NR)  =FACTSQ*CEN(2,NTH,NPH,NR)
            CEN(3,NTH,NPH,NR)  =FACTSQ*CEN(3,NTH,NPH,NR)
            CEP(1,NTH,NPH,NR)  =FACTSQ*CEP(1,NTH,NPH,NR)
            CEP(2,NTH,NPH,NR)  =FACTSQ*CEP(2,NTH,NPH,NR)
            CEP(3,NTH,NPH,NR)  =FACTSQ*CEP(3,NTH,NPH,NR)
         ENDDO
         ENDDO
      ENDDO
      CALL MPSYNC
C
      RETURN
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBFX(CF1,CF2,NTHMAX,NPHMAX)
C
      implicit none
      complex(8),dimension(nthmax,nphmax),intent(in):: CF1
      complex(8),dimension(nthmax,nphmax),intent(out):: CF2
      integer,intent(in):: nthmax,nphmax
      complex(8),dimension(nthmax):: CFM
      complex(8),dimension(nphmax):: CFN
      integer:: NPH,NTH
C
      DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NPH)
         ENDDO
         CALL WMFFFT(CFM,NTHMAX,0)
         DO NTH=1,NTHMAX
            CF2(NTH,NPH)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         DO NPH=1,NPHMAX
            CFN(NPH)=CF2(NTH,NPH)
         ENDDO
         CALL WMFFFT(CFN,NPHMAX,0)
         DO NPH=1,NPHMAX
            CF2(NTH,NPH)=CFN(NPH)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBEX(CF1,CF2,NTHMAX,NPHMAX)
C
      implicit none
      complex(8),dimension(nthmax,nphmax),intent(in):: CF1
      complex(8),dimension(nthmax,nphmax),intent(out):: CF2
      integer,intent(in):: nthmax,nphmax
      complex(8),dimension(nthmax):: CFM
      complex(8),dimension(nphmax):: CFN
      integer:: NPH,NTH
C
      DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NPH)
         ENDDO
         CALL WMFFFT(CFM,NTHMAX,1)
         DO NTH=1,NTHMAX
            CF2(NTH,NPH)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         DO NPH=1,NPHMAX
            CFN(NPH)=CF2(NTH,NPH)
         ENDDO
         CALL WMFFFT(CFN,NPHMAX,1)
         DO NPH=1,NPHMAX
            CF2(NTH,NPH)=CFN(NPH)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** INTERFACE FOR FFT ******
C
      SUBROUTINE WMFFFT(CA,N,KEY)
C
      INCLUDE 'wmcomm.inc'
C
      COMPLEX*16 CA(N)
      DATA NS/0/
C
      IF(N.NE.1) THEN
         IF(N.EQ.NS) THEN
            IND=0
         ELSE
            IND=1
            NS=N
         ENDIF
         IF(KEY.EQ.0) THEN
            CALL FFT2L(CA,CT,RFFT,LFFT,N,IND,KEY)
            DO I=1,N
               CA(I)=CT(I)
            ENDDO
         ELSE
            DO I=1,N
               CT(I)=CA(I)
            ENDDO
            CALL FFT2L(CT,CA,RFFT,LFFT,N,IND,KEY)
         ENDIF
      ENDIF
C
      RETURN
      END
