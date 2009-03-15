C     $Id$

!---- interface for wm parameter

      subroutine get_wmfem_parm(crf_,nth0_,nph0_,mdlwmf_,mdlwmd_)
      
      include '../wm/wmcomm.inc'
      complex(8),intent(out):: crf_
      integer,intent(out):: nth0_,nph0_
      integer,intent(out):: mdlwmf_,mdlwmd_
      crf_=crf
      nth0_=nth0
      nph0_=nph0
      mdlwmf_=mdlwmf
      mdlwmd_=mdlwmd
      return
      end subroutine get_wmfem_parm

!---- interface for wm parameter

      subroutine get_wmfem_parm1(rr_,ra_,rb_)
      
      include '../wm/wmcomm.inc'
      real(8),intent(out):: rr_,ra_,rb_
      rr_=rr
      ra_=ra
      rb_=rb
      return
      end subroutine get_wmfem_parm1

!---- interface for wm parameter

      subroutine get_wmfem_size(nrmax_,nthmax_,nphmax_,nsmax_,ierr)
      
      use wmfem_com, only: rhoa
      include '../wm/wmcomm.inc'
      integer,intent(out):: nrmax_,nthmax_,nphmax_,nsmax_
      integer,intent(out):: ierr
      integer,save::  nrmax_save=0
      integer:: nr

      nrmax_=nrmax
      nthmax_=nthmax
      nphmax_=nphmax
      nsmax_=nsmax

      if(nrmax.ne.nrmax_save) then
         if(associated(rhoa)) deallocate(rhoa)
         allocate(rhoa(nrmax_))
      endif

      do nr=1,nrmax_
         rhoa(nr)=xrho(nr)
      end do
      ierr=0

      return
      end subroutine get_wmfem_size

!     ****** CALCULATE METRIC TENSOR ******

      SUBROUTINE wmfem_metrics(rho,th,ph,gm,gj)

      INCLUDE 'wmcomm.inc'
      real(8),intent(in):: rho,th,ph
      real(8),intent(out),dimension(3,3):: gm
      real(8),intent(out):: gj
      real(8):: rrl,zzl,drrrho,dzzrho,drrchi,dzzchi

      select case(modelg)
      case(0)
         gm(1,1)= ra**2
         gm(1,2)= 0.d0
         gm(1,3)= 0.d0
         gm(2,1)= 0.d0
         gm(2,2)= (ra*rho)**2
         gm(2,3)= 0.d0
         gm(3,1)= 0.d0
         gm(3,2)= 0.d0
         gm(3,3)= rr**2
         gj=rr*ra**2*rho
      case(1,2)
         rcos=cos(th)
         rsin=sin(th)
         rrl   = rr + ra*rho*rcos
         zzl   =      ra*rho*rsin
         drrrho=      ra*rcos
         dzzrho=      ra*rsin
         drrchi=     -ra*rho*rsin
         dzzchi=      ra*rho*rcos
         gm(1,1)= drrrho**2+dzzrho**2
         gm(1,2)= drrrho*drrchi+dzzrho*dzzchi
         gm(1,3)= 0.d0
         gm(2,1)= gm(1,2)
         gm(2,2)= drrchi**2+dzzchi**2
         gm(2,3)= 0.d0
         gm(3,1)= gm(1,3)
         gm(3,2)= gm(2,3)
         gm(3,3)= rrl**2
         gj     = rrl*(drrrho*dzzchi-drrchi*dzzrho)
      case(3,5)
         call wmeq_get_posrz(rho,th,rrl,zzl,
     &                        drrrho,dzzrho,drrchi,dzzchi)
         gm(1,1)= drrrho**2+dzzrho**2
         gm(1,2)= drrrho*drrchi+dzzrho*dzzchi
         gm(1,3)= 0.d0
         gm(2,1)= gm(1,2)
         gm(2,2)= drrchi**2+dzzchi**2
         gm(2,3)= 0.d0
         gm(3,1)= gm(1,3)
         gm(3,2)= gm(2,3)
         gm(3,3)= rrl**2
         gj     = rrl*(drrrho*dzzchi-drrchi*dzzrho)
      case default
         stop 'wmfem_metrics: undefined modelg is input'
      end select
      return
      end subroutine wmfem_metrics
         
!     ****** CALCULATE MAGNETIC FIELD ******

      SUBROUTINE wmfem_magnetic(rho,th,ph,babs,bsupth,bsupph)

      INCLUDE 'wmcomm.inc'
      real(8),intent(in):: rho,th,ph
      real(8),intent(out):: babs,bsupth,bsupph
      real(8):: rrl,qinv

      select case(modelg)
      case(0)
         call wmfem_qprofile(rho,qinv)
         bsupth=(bb*qinv)/rr
         bsupph=bb/rr
         babs=bb*sqrt(1.d0+(ra*rho*qinv/rr)**2)
      case(1,2)
         call wmfem_qprofile(rho,qinv)
         rrl   = rr + ra*rho*cos(th)
         bsupth=(bb*qinv)/rrl
         bsupph=bb/rrl
         babs=bb*sqrt(1.d0+(ra*rho*qinv/rr)**2)*rr/rrl
      case(3,5)
         call wmeq_get_magnetic(rho,th,babs,bsupth,bsupph)
      case default
         stop 'wmfem_magnetic: undefined modelg is input'
      end select
      return
      end subroutine wmfem_magnetic

!     ****** CALCULATE PLASMA DENSITY AND TEMPERATURE ******

      SUBROUTINE wmfem_plasma(rho,nsmax_,rn_,rtpr_,rtpp_,ru_)

      INCLUDE 'wmcomm.inc'
      real(8),intent(in):: rho
      integer,intent(in):: nsmax_
      real(8),dimension(nsmax),intent(out):: rn_,rtpr_,rtpp_,ru_
      INCLUDE '../pl/plcom2.inc'
      
      CALL PLPROF(rho)
      DO ns=1,nsmax_
         rn_(ns)=rn(ns)
         rtpr_(ns)=ptpr(ns)
         rtpp_(ns)=rtpp(ns)
         ru_(ns)=ru(ns)
      ENDDO
      return
      end subroutine wmfem_plasma

!     ***** calculate position R and Z from eqdata ****

      subroutine wmeq_get_posrz(rho,th,rrl,zzl,
     &                          drrrho,dzzrho,drrchi,dzzchi)

      INCLUDE '../eq/eqcomq.inc'
      real(8),intent(in):: rho,th
      real(8),intent(out):: rrl,zzl,drrrho,dzzrho,drrchi,dzzchi

      psipl=fnpsip(rho)
      dpsipl=fndpsip(rho)
      CALL spl2dd(th,psipl,rrl,drrchi,drrpsi,
     &                  THIT,PSIP,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      CALL spl2dd(th,psipl,zzl,dzzchi,dzzpsi,
     &                  THIT,PSIP,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      drrrho=drrpsi*dpsipl
      dzzrho=dzzpsi*dpsipl
      return
      end subroutine wmeq_get_posrz

!     ***** calculated magnetic field from eqdata ****

      subroutine wmeq_get_magnetic(rho,th,babs,bsupth,bsupph)

      INCLUDE '../eq/eqcomq.inc'
      real(8),intent(in):: rho,th
      real(8),intent(out):: babs,bsupth,bsupph
      real(8):: rrl,zzl,drrpsi,dzzpsi,drrchi,dzzchi,rhol
      real(8):: bprr,bpzz,bthl,bphl,ttl
      real(8),dimension(3,3):: gm

      IF(rho.LE.1.D-8) THEN
         rhol=1.D-8
      ELSE
         rhol=rho
      ENDIF

      psipl=fnpsip(rhol)
      CALL spl2dd(th,psipl,rrl,drrchi,drrpsi,
     &                  THIT,PSIP,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      CALL spl2dd(th,psipl,zzl,dzzchi,dzzpsi,
     &                  THIT,PSIP,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      gm(2,2)= drrchi**2+dzzchi**2
      gm(2,3)= 0.d0
      gm(3,3)= rrl**2

      CALL EQPSID(rrl,zzl,DPSIDR,DPSIDZ)

      bprr= DPSIDZ/(2.D0*PI*rrl)
      bpzz=-DPSIDR/(2.D0*PI*rrl)
      ttl=FNTTS(rho)
      bthl=SQRT(bprr**2+bpzz**2)
      bphl= ttl/(2.d0*PI*rrl)
      bsupth=bthl/sqrt(gm(2,2))
      bsupph=bphl/sqrt(gm(3,3))
      babs=sqrt(     gm(2,2)*bsupth*bsupth
     &         +2.d0*gm(2,3)*bsupth*bsupph
     &         +     gm(3,3)*bsupph*bsupph)
!      write(6,'(A,1P4E12.4)') '%%%--',bthl,bphl,bsupth,bsupph
      return

      end subroutine wmeq_get_magnetic

!     ****** CALCULATE Q PROFILE ******

      SUBROUTINE wmfem_qprofile(rho,qinv)

      INCLUDE 'wmcomm.inc'
      real(8),intent(in):: rho
      real(8),intent(out):: qinv
      real(8):: ql,qsa0,qsaa

      select case(modelg)
      case(0,1,2)
         select case(modelq)
         case(0)
            if(rho.gt.1.d0) then
               ql=qa*rho**2
            else if(rhomin.le.0.d0) then
               ql=(q0-qa)*(1.d0-rho**2)+qa
            else
               qsa0    =1/q0-1/qmin
               qsaa    =1/qa-1/qmin
               if(rho.le.rhomin)then
                  ql=1/(1/q0-qsa0*(3*rho**2/rhomin**2
     &                  -2*rho**3/rhomin**3))
               else
                  ql=1/(1/qmin+3*qsa0*(rho-rhomin)**2/rhomin**2
     &               +(qsaa-3*qsa0*(1-rhomin)**2/rhomin**2)
     &                *(rho-rhomin)**3/(1-rhomin)**3)
               end if
            end if
            qinv=1.d0/ql
         case(1)
            if(rip.eq.0.d0) then
               qinv=0.d0
            else
               qa=2.d0*pi*ra*ra*bb/(rmu0*rip*1.d6*rr)
               q0=qa/(1.d0+profj)
               if(rho.ge.1.d0) then
                  ql= qa*rho**2
               else
                  ql= qa*rho**2/(1.d0-(1.d0-rho**2)**(profj+1.d0))
               end if
               qinv=1.d0/ql
            end if
         case(2)
            if(rho.gt.1.d0) then
               qinv=1.d0/(qa*rho**2)
            else
               qinv=(1.d0/q0-1.d0/qa)*(1.d0-rho**2)+1.d0/qa
            end if
         end select
      case(3,5)
         call getqp(rho,ql)
         qinv=1.d0/ql
      case default
         stop 'wmfem_qprofile: undefined modelg is input'
      end select
      return
      end subroutine wmfem_qprofile
C
C     ****** CALCULATE ANTENNA CURRENT ******
C     
      SUBROUTINE WMFEM_SETJ(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RHO(0:3,1:3,2:3)
C
      DATA RHO/ 3.832, 1.841, 3.054, 4.201,
     &          7.016, 5.332, 6.706, 8.015,
     &         10.173, 8.536, 9.969,11.346,
     &          2.405, 3.832, 5.136, 6.380,
     &          5.520, 7.016, 8.417, 9.761,
     &          8.654,10.173,11.620,13.015/
C
         IF(RD.LE.RA.OR.RD.GE.RB) THEN
            IF(MYRANK.EQ.0) WRITE(6,*) '!! WMSETJ: RD = (RA+RB)/2'
            RD=0.5D0*(RA+RB)
            IF(MYRANK.EQ.0) 
     &           WRITE(6,'(A,1P3E12.4)') 'RA,RB,RD=',RA,RB,RD
         ENDIF
C
      DO NDX=1,NDSIZ
      DO MDX=1,MDSIZ
         CJANT(1,MDX,NDX)=(0.D0,0.D0)
         CJANT(2,MDX,NDX)=(0.D0,0.D0)
         CJANT(3,MDX,NDX)=(0.D0,0.D0)
      ENDDO
      ENDDO
C
      MODELWG=0
      IF(MODELJ.EQ.0) THEN
         CALL WMFEM_CANT
      ELSEIF(MODELJ.EQ.1) THEN
         MODELWG=1
      ELSEIF(MODELJ.EQ.2) THEN
         CJANT(2,1,1)= 1.D0
      ELSEIF(MODELJ.EQ.3) THEN
         CJANT(3,1,1)= 1.D0
      ELSE
         NMODE=MOD(MODELJ,10)
         IMODE=MODELJ/10
         IF(NTH0.LT.0.OR.NTH0.GT.3) GOTO 9000
         IF(NMODE.LT.1.OR.NMODE.GT.3) GOTO 9000
         IF(IMODE.LT.2.OR.IMODE.GT.3) GOTO 9000
         RF =0.5D0*VC/PI
     &            *SQRT((NPH0/RR)**2+(RHO(NTH0,NMODE,IMODE)/RB)**2)
     &            *1.D-6
         RFI=0.D0
         CRF=DCMPLX(RF,RFI)
         CJANT(IMODE,1,1)=(1.D0,0.D0)
      ENDIF
C
      IERR=0
      RETURN
C
 9000 WRITE(6,*) 'XX WMCALJ ERROR'
      IERR=1
      RETURN
      END
C
C     ****** CALCULATE ANTENNA CURRENT ******
C
      SUBROUTINE WMFEM_CANT
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CJT(MDM,NDM,NAM)
      DIMENSION CJZ(MDM,NDM,NAM)
C
      DO NA=1,NAMAX
         TH1=THJ1(NA)*PI/180.D0
         TH2=THJ2(NA)*PI/180.D0
         PH1=PHJ1(NA)*PI/180.D0
         PH2=PHJ2(NA)*PI/180.D0
C
         CAJ=AJ(NA)*EXP(DCMPLX(0.D0,APH(NA)*PI/180.D0))
C   
      DO NDX=1,NPHMAX
         ND=NDX-1
         IF(NPHMAX.GT.1.AND.NDX.GT.NPHMAX/2) ND=ND-NPHMAX
         NN=NPH0+NHC*ND
         IF(NN.EQ.0.OR.ABS(PH2-PH1).LE.1.D-15) THEN
            CJN=-CI
         ELSE
            CJN=(EXP(-CI*NN*PH2)-EXP(-CI*NN*PH1))/(NN*(PH2-PH1))
         ENDIF
      DO MDX=1,NTHMAX
         MD=MDX-1
         IF(NTHMAX.GT.1.AND.MDX.GT.NTHMAX/2) MD=MD-NTHMAX
         MM=NTH0+MD
         IF(ABS(MM+BETAJ).LE.0.D0) THEN
            CJMP=-CI*(TH2-TH1)
         ELSE
            CJMP=(EXP(-CI*(MM+BETAJ)*TH2)-EXP(-CI*(MM+BETAJ)*TH1))
     &           /(MM+BETAJ)
         ENDIF
         IF(ABS(MM-BETAJ).LE.0.D0) THEN
            CJMM=-CI*(TH2-TH1)
         ELSE
            CJMM=(EXP(-CI*(MM-BETAJ)*TH2)-EXP(-CI*(MM-BETAJ)*TH1))
     &           /(MM-BETAJ)
         ENDIF
         CJTEMP=CAJ/(8*PI**2)*CJN*(CJMP+CJMM)
         CJT(MDX,NDX,NA)=CJTEMP*COS(2*PI*ANTANG/360.D0)
         CJZ(MDX,NDX,NA)=CJTEMP*SIN(2*PI*ANTANG/360.D0)
      ENDDO
      ENDDO
      ENDDO
C
      DO NDX=1,NPHMAX
      DO MDX=1,NTHMAX
      DO NA=1,NAMAX
         CJANT(2,MDX,NDX)=CJANT(2,MDX,NDX)+CJT(MDX,NDX,NA)
         CJANT(3,MDX,NDX)=CJANT(3,MDX,NDX)+CJZ(MDX,NDX,NA)
         IF(MYRANK.EQ.0) THEN
            IF(NPRINT.GE.3) WRITE(6,'(A,2I4,1P2E15.7)') 
     &                   'NN,MM,CJANT=',
     &                   NPH0+NHC*ND,NTH0+MD,CJANT(2,MDX,NDX)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
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
                  FACTM=(XRHO2-RD/RA)/DRHO
                  FACTP=(RD/RA-XRHO1)/DRHO
C
                  DO NDX=1,NPHMAX
                     ND=NDX-1
                     IF(NPHMAX.GT.1.AND.NDX.GT.NPHMAX/2) ND=ND-NPHMAX
                     NN=NPH0+NHC*ND
                  DO MDX=1,NTHMAX
                     MD=MDX-1
                     IF(NTHMAX.GT.1.AND.MDX.GT.NTHMAX/2) MD=MD-NTHMAX
                     MM=NTH0+MD

                     CJTHM=CC*CJANT(2,MDX,NDX)*FACTM*XRHOC
     &                                        /(DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*FACTM*XRHOC
     &                                        /(DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)
     &                                         *DRHO/XRHOC
                     CJTHP=CC*CJANT(2,MDX,NDX)*FACTP*XRHOC
     &                                        /(DRHO*DPH)
                     CJPHP=CC*CJANT(3,MDX,NDX)*FACTP*XRHOC
     &                                        /(DRHO*DTH)
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
                  DO NDX=1,NPHMAX
                     ND=NDX-1
                     IF(NPHMAX.GT.1.AND.NDX.GT.NPHMAX/2) ND=ND-NPHMAX
                     NN=NPH0+NHC*ND
                  DO MDX=1,NTHMAX
                     MD=MDX-1
                     IF(NTHMAX.GT.1.AND.MDX.GT.NTHMAX/2) MD=MD-NTHMAX
                     MM=NTH0+MD

                     CJTHM=CC*CJANT(2,MDX,NDX)*XRHOC
     &                             /(DRHO*DPH)
                     CJPHM=CC*CJANT(3,MDX,NDX)*XRHOC
     &                             /(DRHO*DTH)
                     CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DRHO
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
               CJTHM=RJFACT*XRHOC/(DRHO*DPH)
               CJPHM=RJFACT*XRHOC/(DRHO*DTH)
               CJR  =-(CI*MM*CJTHM+CI*NN*CJPHM)*DRHO/XRHOC
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
C     ****** CALCULATE WAVE ELECTRIC FIELD ******
C
      SUBROUTINE WMFEM_EFLD
C
      use wmfem_com, only: cef
      INCLUDE 'wmcomm.inc'

      DIMENSION CEF1(MDM,NDM),CEF2(MDM,NDM),RMA(3,3)
C
      do nr=1,nrmax
         DO NDX=1,nphmax
            if(nphmax.eq.1) then
               NDX1=NDX
            else
               NDX1=NDX+nphmax/2-1
               IF(NDX1.gt.nphmax) NDX1=NDX1-nphmax
            endif
         DO MDX=1,nthmax
            if(nthmax.eq.1) then
               MDX1=MDX
            else
               MDX1=MDX+nthmax/2-1
               IF(MDX1.gt.nthmax) MDX1=MDX1-nthmax
            endif
            CEFLDK(1,MDX1,NDX1,NR )=cef(1,mdx,ndx,nr)
            CEFLDK(2,MDX1,NDX1,NR )=cef(2,mdx,ndx,nr)
            CEFLDK(3,MDX1,NDX1,NR )=cef(3,mdx,ndx,nr)
            CEFLD(1,MDX,NDX,NR )=cef(1,mdx,ndx,nr)
            CEFLD(2,MDX,NDX,NR )=cef(2,mdx,ndx,nr)
            CEFLD(3,MDX,NDX,NR )=cef(3,mdx,ndx,nr)
         enddo
         enddo
      enddo
C
C     ----- Inverse Fourier transform ----
C
      DO NR=1,NRMAX
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
      DO NR=1,NRMAX
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
c$$$         DO NR=1,NRMAX
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
C     ****** CALCULATE WAVE MAGNETIC FIELD ******
C
      SUBROUTINE WMFEM_BFLD
C
      use wmfem_com, only: cbf
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CBF1(MDM,NDM),CBF2(MDM,NDM),RMA(3,3)
C
      do nr=1,nrmax
         DO NDX=1,nphmax
            if(nphmax.eq.1) then
               NDX1=NDX
            else
               NDX1=NDX+nphmax/2-1
               IF(NDX1.gt.nphmax) NDX1=NDX1-nphmax
            endif
         DO MDX=1,nthmax
            if(nthmax.eq.1) then
               MDX1=MDX
            else
               MDX1=MDX+nthmax/2-1
               IF(MDX1.gt.nthmax) MDX1=MDX1-nthmax
            endif
            CBFLDK(1,MDX1,NDX1,NR )=cbf(1,mdx,ndx,nr)
            CBFLDK(2,MDX1,NDX1,NR )=cbf(2,mdx,ndx,nr)
            CBFLDK(3,MDX1,NDX1,NR )=cbf(3,mdx,ndx,nr)
            CBFLD(1,MDX,NDX,NR )=cbf(1,mdx,ndx,nr)
            CBFLD(2,MDX,NDX,NR )=cbf(2,mdx,ndx,nr)
            CBFLD(3,MDX,NDX,NR )=cbf(3,mdx,ndx,nr)
         enddo
         enddo
      enddo
C
C     ----- Inverse Fourier transform ----
C
      DO NR=1,NRMAX
      DO I=1,3
         DO NDX=1,NDSIZ
         DO MDX=1,MDSIZ
            CBF1(MDX,NDX)=CBFLD(I,MDX,NDX,NR)
         ENDDO
         ENDDO
         CALL WMSUBEX(CBF1,CBF2,NTHMAX,NPHMAX)
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            CBFLD(I,NTH,NPH,NR)=CBF2(NTH,NPH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- Calculate CBN from CBsup -----
C
C
      DO NR=1,NRMAX
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         CBN(1,NTH,NPH,NR)=CBFLD(1,NTH,NPH,NR)
         CBN(2,NTH,NPH,NR)=CBFLD(2,NTH,NPH,NR)
         CBN(3,NTH,NPH,NR)=CBFLD(3,NTH,NPH,NR)
         CBP(1,NTH,NPH,NR)=(   CBN(1,NTH,NPH,NR)
     &                     +CI*CBN(2,NTH,NPH,NR))/SQRT(2.D0)
         CBP(2,NTH,NPH,NR)=(   CBN(1,NTH,NPH,NR)
     &                     -CI*CBN(2,NTH,NPH,NR))/SQRT(2.D0)
         CBP(3,NTH,NPH,NR)=    CBN(3,NTH,NPH,NR)
      ENDDO
      ENDDO
      ENDDO
c$$$C
c$$$      DO K=1,3
c$$$         DO NR=1,NRMAX
c$$$            IF(XRHO(NR).GT.1.0D0) THEN
c$$$               DO NTH=1,NTHMAX
c$$$                  DO NPH=1,NPHMAX
c$$$                     CBP(K,NTH,NPH,NR)=(0.D0,0.D0)
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
      SUBROUTINE WMFEM_PABS
C
      use wmfem_com, only: cpp,cpa,nthmax2,nphmax2
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION DS(NRM),DSS(NTHM,NPHM,NRM)
      DIMENSION CPF1(nthmax2,nphmax2),CPF2(nthmax2,nphmax2)
      real(8),dimension(3,3)::  gm
      real(8):: gj1,gj2
C
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NPHMAX)
C
      DO NR=1,nrmax
      DO NS=1,NSMAX
      DO NDX=1,NDSIZ
      DO KKX=1,KDSIZ
      DO MDX=1,MDSIZ
      DO LLX=1,LDSIZ
         CPABS(LLX,MDX,KKX,NDX,NS,NR)=cpp(MDX,NDX,LLX,KKX,NR,NS)
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
            KKX=1
            LLX=1
         DO NDX=1,nphmax
            if(nphmax.eq.1) then
               NDX1=NDX
            else
               NDX1=NDX+nphmax/2-1
               IF(NDX1.gt.nphmax) NDX1=NDX1-nphmax
            endif
         DO MDX=1,nthmax
            if(nthmax.eq.1) then
               MDX1=MDX
            else
               MDX1=MDX+nthmax/2-1
               IF(MDX1.gt.nthmax) MDX1=MDX1-nthmax
            endif
            PABSK(MDX1,NDX1,NR,NS)=DBLE(cpp(MDX,NDX,LLX,KKX,NR,NS))
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
               DO KKX=1,nphmax2
               DO LLX=1,nthmax2
                  CPF1(LLX,KKX)=cpp(MDX,NDX,LLX,KKX,NR,NS)
               ENDDO
               ENDDO
               CALL WMSUBEX(CPF1,CPF2,NTHMAX2,NPHMAX2)
               DO NPH=1,NPHMAX
               DO NTH=1,NTHMAX
                  PABS(NTH,NPH,NR,NS)=PABS(NTH,NPH,NR,NS)
     &                               +DBLE(CPF2(2*NTH-1,2*NPH-1))
               ENDDO
               ENDDO
            ENDDO
            ENDDO
c$$$            if(nr.eq.1.or.nr.eq.2) then
c$$$            DO NPH=1,NPHMAX
c$$$            DO NTH=1,NTHMAX
c$$$               write(6,'(3I8,1P2E12.4)') NR,NTH,NPH,
c$$$     &                                   PABS(NTH,NPH,NR,NS)
c$$$            ENDDO
c$$$            ENDDO
c$$$            endif
         ENDDO
         ENDDO

C     +++++ Antenna impedance +++++

      cradtt=(0.d0,0.d0)
      do nph=1,nphmax
         do nth=1,nthmax
            cradtt=cradtt+cpa(nth,nph)
         end do
      end do

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
c$$$         CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
c$$$         VTE=SQRT(RTPR(1)*AEE*1.D3/(PA(1)*AMP))
c$$$         WW=DBLE(CW)
c$$$         IF(RN(1).LE.0.D0) THEN
c$$$            RLNLMD=15.D0
c$$$         ELSE
c$$$            RT=(RTPR(1)+2*RTPP(1))/3.D0
c$$$            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(1))
c$$$     &           + 2.30D0*LOG10(RT)
c$$$         ENDIF
c$$$         DO ND=NDMIN,NDMAX
c$$$            NDX=ND-NDMIN+1
c$$$            NN=NPH0+ND
c$$$         DO MD=MDMIN,MDMAX
c$$$            MDX=MD-MDMIN+1
c$$$            MM=NTH0+MD
c$$$            DO KKX=1,KDSIZ
c$$$            DO LLX=1,LDSIZ
c$$$               CPF1(LLX,KKX)=CPABS(LLX,MDX,KKX,NDX,NS,NR)
c$$$            ENDDO
c$$$            ENDDO
c$$$            CALL WMSUBEX(CPF1,CPF2,NTHMAX,NPHMAX)
c$$$            DO NPH=1,NPHMAX
c$$$            DO NTH=1,NTHMAX
c$$$               CALL WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
c$$$               RKPR=MM*BSUPTH/BABS+NN*BSUPPH/BABS
c$$$              IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
c$$$              IF(ABS(WW/RKPR).LT.VC) THEN
c$$$                 W=WW/(RKPR*VTE)
c$$$                 XL=(RPST(NTH,NPH,NR)-RR  )/RR
c$$$                 YL=(ZPST(NTH,NPH,NR)-0.D0)/RR
c$$$                 EFCD=W1CDEF(ABS(W),ZEFF,XL,YL,1)
c$$$                 IF(W.LT.0.D0) EFCD=-EFCD
c$$$                 IF (RN(1).GT.0.D0) THEN
c$$$                    PCUR(NTH,NPH,NR)=PCUR(NTH,NPH,NR)
c$$$     &                   +0.384D0*RTPR(1)/(AEE*1.D3)*EFCD
c$$$     &                   /((RN(1)/1.D20)*RLNLMD)*DBLE(CPF2(NTH,NPH))
c$$$     &                   /(2.D0*PI*RPST(NTH,NPH,NR))
c$$$                 END IF
c$$$              ENDIF
c$$$            ENDDO
c$$$            ENDDO
c$$$         ENDDO
c$$$         ENDDO
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
      DO NR=1,NRMAX-1
         DS(NR)=0.D0
         DRHO=XRHO(NR+1)-XRHO(NR)
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            th=dth*(nth-1)
            ph=dph*(nph-1)
            call wmfem_metrics(xrho(nr),th,ph,gm,gj1)
            call wmfem_metrics(xrho(nr+1),th,ph,gm,gj2)
!            write(6,'(I5,1P6E12.4)') 
!     &           nr,xrho(nr),xrho(nr+1),th,ph,gj1,gj2
            DSSS=DRHO*0.5d0*(gj1+gj2)
            DSS(NTH,NPH,NR)=1.D0/DSSS
            DS(NR)=DS(NR)+DSSS*DTH*DPH
         ENDDO
         ENDDO
         DS(NR)=1.D0/DS(NR)
      ENDDO
      DS(NRMAX)=0.D0
      DO NTH=1,NTHMAX
         DO NPH=1,NPHMAX
            DSS(NTH,NPH,NRMAX)=0.d0
         ENDDO
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
!      write(6,'(2I5,1P2E12.4)') 
!     &     (NR,NS,DS(NR),PABSR(NR,NS),NR=1,NRMAX)
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
      DO NR=1,NRMAX
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
