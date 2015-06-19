C     $Id$

!---- interface for wm parameter

      subroutine get_wmfem_parm(crf_,nth0_,nph0_,mdlwmf_,mdlwmd_)
      
      include 'wmcomm.inc'
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
      
      include 'wmcomm.inc'
      real(8),intent(out):: rr_,ra_,rb_
      rr_=rr
      ra_=ra
      rb_=rb
      return
      end subroutine get_wmfem_parm1

!---- interface for wm parameter

      subroutine get_wmfem_size(nrmax_,nthmax_,nhhmax_,nsmax_)
      
      use wmfem_comm, only: rhoa
      include 'wmcomm.inc'
      integer,intent(out):: nrmax_,nthmax_,nhhmax_,nsmax_
      integer,save::  nrmax_save=0
      integer:: nr

      nrmax_=nrmax
      nthmax_=nthmax
      nhhmax_=nhhmax
      nsmax_=nsmax

      if(nrmax.ne.nrmax_save) then
         if(allocated(rhoa)) deallocate(rhoa)
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
!      SUBROUTINE wmfem_magnetic(rho,th,ph,babs,bsupr,bsupth,bsupph)

      INCLUDE 'wmcomm.inc'
      real(8),intent(in):: rho,th,ph
      real(8),intent(out):: babs,bsupth,bsupph
!      real(8),intent(out):: babs,bsupr,bsupth,bsupph
      real(8):: rrl,qinv

      select case(modelg)
      case(0)
         call wmfem_qprofile(rho,qinv)
!         bsupr=0.d0
!         bsupth=(bb*qinv)/rr
         bsupth=0d0
         bsupph=bb/rr
         babs=bb*sqrt(1.d0+(ra*rho*qinv/rr)**2)
      case(1,2)
         call wmfem_qprofile(rho,qinv)
         rrl   = rr + ra*rho*cos(th)
!         bsupr=0.d0
         bsupth=(bb*qinv)/rrl
         bsupph=bb/rrl
         babs=bb*sqrt(1.d0+(ra*rho*qinv/rr)**2)*rr/rrl
      case(3,5)
!         bsupr=0.d0
         call wmeq_get_magnetic(rho,th,babs,bsupth,bsupph)
      case default
         stop 'wmfem_magnetic: undefined modelg is input'
      end select
      return
      end subroutine wmfem_magnetic

c$$$!     ****** CALCULATE PLASMA DENSITY AND TEMPERATURE ******
c$$$
c$$$      SUBROUTINE wmfem_plasma(rho,nsmax_,rn_,rtpr_,rtpp_,ru_)
c$$$
c$$$      INCLUDE 'wmcomm.inc'
c$$$      integer,intent(in):: nsmax_
c$$$      real(8),intent(in):: rho
c$$$      real(8),dimension(nsmax),intent(out):: rn_,rtpr_,rtpp_,ru_
c$$$      INCLUDE '../pl/plcom2.inc'
c$$$      
c$$$      CALL PLPROF(rho)
c$$$      DO ns=1,nsmax_
c$$$         rn_(ns)=rn(ns)
c$$$         rtpr_(ns)=ptpr(ns)
c$$$         rtpp_(ns)=rtpp(ns)
c$$$         ru_(ns)=ru(ns)
c$$$      ENDDO
c$$$      return
c$$$      end subroutine wmfem_plasma

!     ***** calculate position R and Z from eqdata ****

      subroutine wmeq_get_posrz(rho,th,rrl,zzl,
     &                          drrrho,dzzrho,drrchi,dzzchi)

      INCLUDE '../eq/eqcomq.inc'
      real(8),intent(in):: rho,th
      real(8),intent(out):: rrl,zzl,drrrho,dzzrho,drrchi,dzzchi

!!!!seki       CALL spl2dd(th,rho,rrl,drrchi,drrrho,
!!!!seki      &                  THIT,RHOT,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
!!!!seki       CALL spl2dd(th,rho,zzl,dzzchi,dzzrho,
!!!!seki      &                  THIT,RHOT,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      pause
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
!!!!seki       CALL spl2dd(th,rhol,rrl,drrchi,drrrho,
!!!!seki      &                  THIT,RHOT,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
!!!!seki       CALL spl2dd(th,rhol,zzl,dzzchi,dzzrho,
!!!!seki      &                  THIT,RHOT,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      pause
      gm(2,2)= drrchi**2+dzzchi**2
      gm(2,3)= 0.d0
      gm(3,3)= rrl**2

      CALL psigd(rrl,zzl,DPSIDR,DPSIDZ)

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
c$$$      if(abs(th).le.1.d-6.and.rhol.le.0.01) then
c$$$         write(6,'(A,1P4E12.4)') 'rhol,th,ph,babs:',rhol,th,ph,babs
c$$$         write(6,'(A,1P3E12.4)') 'gm:22,33:ttl:   ',gm(2,2),gm(3,3),ttl
c$$$         write(6,'(A,1P4E12.4)') 'rr,zz,bprr,bpzz:',rrl,zzl,bprr,bpzz
c$$$         write(6,'(A,1P4E12.4)') 'bthl,bphl,bsup; ',bthl,bphl,
c$$$     &        bsupth,bsupph
c$$$         write(6,*)
c$$$      endif
      return

      end subroutine wmeq_get_magnetic

!     ***** calculated magnetic field from eqdata ****

      subroutine wmeq_get_mtxCL(nthmax2,nhhmax2,mtxcl)

      INCLUDE '../eq/eqcomq.inc'
      INTEGER,INTENT(IN):: nthmax2,nhhmax2
      COMPLEX(8),DIMENSION(3,3,nthmax2,nhhmax2),INTENT(OUT):: mtxcl
      real(8):: rrl,zzl,drrpsi,dzzpsi,drrchi,dzzchi,rhol
      real(8):: bprr,bpzz,bthl,bphl,ttl,absdrho,bbl
      real(8):: dth2,dph2,rho
      INTEGER:: nth2,nhh2,i,j
      real(8),dimension(3,3):: em
      COMPLEX(8),dimension(nthmax2,nhhmax2):: cf1,cf2

      rho=0.d0
      dth2=2.d0*pi/nthmax2
      dph2=2.d0*pi/nhhmax2
      do nhh2=1,nhhmax2
         do nth2=1,nthmax2
            th=dth2*(nth2-1)
            ph=dph2*(nhh2-1)
         
            select case(modelg)
            case(0,1,2)
               drrrho=ra*cos(th)
               dzzrho=ra*sin(th)
            case(3,5)
               psipl=fnpsip(rho)
!!!!seki                CALL spl2dd(th,rho,rrl,drrchi,drrrho,
!!!!seki      &                     THIT,RHOT,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
!!!!seki                CALL spl2dd(th,rho,zzl,dzzchi,dzzrho,
!!!!seki      &                     THIT,RHOT,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
              pause
            end select
            absdrho=sqrt(drrrho**2+dzzrho**2)
            em(1,1)=drrrho/absdrho
            em(1,2)=dzzrho/absdrho
            em(1,3)=0.D0

            SELECT CASE(modelg)
            CASE(0,1,2)
               call wmfem_qprofile(rho,qinv)
               rrl=rr+ra*rho*cos(th)
               zzl=   ra*rho*sin(th)
               bphl=bb*rr/rrl
               bprr=-bb*qinv*ra*rho/rrl*sin(th)
               bpzz= bb*qinv*ra*rho/rrl*cos(th)
            CASE(3,5)
               CALL psigd(rrl,zzl,DPSIDR,DPSIDZ)
               bprr= DPSIDZ/(2.D0*PI*rrl)
               bpzz=-DPSIDR/(2.D0*PI*rrl)
               ttl=FNTTS(rhol)
               bphl= ttl/(2.d0*PI*rrl)
               bprr= DPSIDZ/(2.D0*PI*rrl)
               bpzz=-DPSIDR/(2.D0*PI*rrl)
            END SELECT
            bbl=sqrt(bphl**2+bprr**2+bpzz**2)
            em(3,1)=bprr/bbl
            em(3,2)=bpzz/bbl
            em(3,3)=bphl/bbl
            em(2,1)=em(3,2)*em(1,3)-em(3,3)*em(1,2)
            em(2,2)=em(3,3)*em(1,1)-em(3,1)*em(1,3)
            em(2,3)=em(3,1)*em(1,2)-em(3,2)*em(1,1)
            DO j=1,3
               DO i=1,3
                  mtxCL(i,j,nth2,nhh2)=em(i,j)
               END DO
            END DO
         END DO
      ENDDO

!      DO nhh2=1,nhhmax2
!         DO nth2=1,nthmax2
!            th=dth2*(nth2-1)
!            ph=dph2*(nhh2-1)
!            write(21,'(A,1P2E12.4)') 'th,ph=',th,ph
!            write(21,'(A,1P6E12.4)') 
!     &           '  e1=',(mtxcl(1,j,nth2,nhh2),j=1,3)
!            write(21,'(A,1P6E12.4)') 
!     &           '  e2=',(mtxcl(2,j,nth2,nhh2),j=1,3)
!            write(21,'(A,1P6E12.4)') 
!     &           '  e3=',(mtxcl(3,j,nth2,nhh2),j=1,3)
!         END DO
!      END DO

!     ----- Fourier transform -----
      DO j=1,3
         DO i=1,3
            DO nhh2=1,nhhmax2
               DO nth2=1,nthmax2
                  cf1(nth2,nhh2)=mtxCL(i,j,nth2,nhh2)
               END DO
            END DO
            CALL WMSUBFX(cf1,cf2,nthmax2,nhhmax2)
            DO nhh2=1,nhhmax2
               DO nth2=1,nthmax2
                  mtxCL(i,j,nth2,nhh2)=cf2(nth2,nhh2)
               END DO
            END DO
         END DO
      ENDDO

!      DO nhh2=1,nhhmax2
!         DO nth2=1,nthmax2
!            write(21,'(A,2I6)') 'mm2,nn2=',nth2-1,nhh2-1
!            write(21,'(A,1P6E12.4)') 
!     &           '  e1=',(mtxcl(1,j,nth2,nhh2),j=1,3)
!            write(21,'(A,1P6E12.4)') 
!     &           '  e2=',(mtxcl(2,j,nth2,nhh2),j=1,3)
!            write(21,'(A,1P6E12.4)') 
!     &           '  e3=',(mtxcl(3,j,nth2,nhh2),j=1,3)
!         END DO
!      END DO

      RETURN
      END SUBROUTINE wmeq_get_mtxCL

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
            IF(NRANK.EQ.0) WRITE(6,*) '!! WMSETJ: RD = (RA+RB)/2'
            RD=0.5D0*(RA+RB)
            IF(NRANK.EQ.0) 
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
      DO NDX=1,NHHMAX
         ND=NDX-1
         IF(NHHMAX.GT.1.AND.NDX.GT.NHHMAX/2) ND=ND-NHHMAX
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
      DO NDX=1,NHHMAX
      DO MDX=1,NTHMAX
      DO NA=1,NAMAX
         CJANT(2,MDX,NDX)=CJANT(2,MDX,NDX)+CJT(MDX,NDX,NA)
         CJANT(3,MDX,NDX)=CJANT(3,MDX,NDX)+CJZ(MDX,NDX,NA)
         IF(NRANK.EQ.0) THEN
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
      USE plprof,ONLY: pl_prof2
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CFVP(nhhmax,nthmax,3)
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
               CW=2.D0*PI*CRF*1.D6
               CC=CI*CW*RMU0
               DPH=2.D0*PI/NHHMAX
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
                  DO NDX=1,NHHMAX
                     ND=NDX-1
                     IF(NHHMAX.GT.1.AND.NDX.GT.NHHMAX/2) ND=ND-NHHMAX
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
                  DO NDX=1,NHHMAX
                     ND=NDX-1
                     IF(NHHMAX.GT.1.AND.NDX.GT.NHHMAX/2) ND=ND-NHHMAX
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
            DPH=2.D0*PI/NHHMAX
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
            CALL pl_prof2(XRHOC,RN,RTPR,RTPP,RU)
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
      use wmfem_comm, only: cef
      INCLUDE 'wmcomm.inc'

      DIMENSION CEF1(MDM,NDM),CEF2(MDM,NDM),RMA(3,3)
C
      do nr=1,nrmax
         DO NDX=1,nhhmax
            if(nhhmax.eq.1) then
               NDX1=NDX
            else
               NDX1=NDX+nhhmax/2-1
               IF(NDX1.gt.nhhmax) NDX1=NDX1-nhhmax
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
         CALL WMSUBEX(CEF1,CEF2,NTHMAX,NHHMAX)
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CEFLD(I,NTH,NHH,NR)=CEF2(NTH,NHH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- Calculate CEN from CEsup -----
C
C
      DO NR=1,NRMAX
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CEN(1,NTH,NHH,NR)=CEFLD(1,NTH,NHH,NR)
         CEN(2,NTH,NHH,NR)=CEFLD(1,NTH,NHH,NR)
         CEN(3,NTH,NHH,NR)=CEFLD(1,NTH,NHH,NR)
         CEP(1,NTH,NHH,NR)=(   CEN(1,NTH,NHH,NR)
     &                     +CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
         CEP(2,NTH,NHH,NR)=(   CEN(1,NTH,NHH,NR)
     &                     -CI*CEN(2,NTH,NHH,NR))/SQRT(2.D0)
         CEP(3,NTH,NHH,NR)=    CEN(3,NTH,NHH,NR)
      ENDDO
      ENDDO
      ENDDO
c$$$C
c$$$      DO K=1,3
c$$$         DO NR=1,NRMAX
c$$$            IF(XRHO(NR).GT.1.0D0) THEN
c$$$               DO NTH=1,NTHMAX
c$$$                  DO NHH=1,NHHMAX
c$$$                     CEP(K,NTH,NHH,NR)=(0.D0,0.D0)
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
      use wmfem_comm, only: cbf
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CBF1(MDM,NDM),CBF2(MDM,NDM),RMA(3,3)
C
      do nr=1,nrmax
         DO NDX=1,nhhmax
            if(nhhmax.eq.1) then
               NDX1=NDX
            else
               NDX1=NDX+nhhmax/2-1
               IF(NDX1.gt.nhhmax) NDX1=NDX1-nhhmax
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
         CALL WMSUBEX(CBF1,CBF2,NTHMAX,NHHMAX)
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CBFLD(I,NTH,NHH,NR)=CBF2(NTH,NHH)
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- Calculate CBN from CBsup -----
C
C
      DO NR=1,NRMAX
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         CBN(1,NTH,NHH,NR)=CBFLD(1,NTH,NHH,NR)
         CBN(2,NTH,NHH,NR)=CBFLD(2,NTH,NHH,NR)
         CBN(3,NTH,NHH,NR)=CBFLD(3,NTH,NHH,NR)
         CBP(1,NTH,NHH,NR)=(   CBN(1,NTH,NHH,NR)
     &                     +CI*CBN(2,NTH,NHH,NR))/SQRT(2.D0)
         CBP(2,NTH,NHH,NR)=(   CBN(1,NTH,NHH,NR)
     &                     -CI*CBN(2,NTH,NHH,NR))/SQRT(2.D0)
         CBP(3,NTH,NHH,NR)=    CBN(3,NTH,NHH,NR)
      ENDDO
      ENDDO
      ENDDO
c$$$C
c$$$      DO K=1,3
c$$$         DO NR=1,NRMAX
c$$$            IF(XRHO(NR).GT.1.0D0) THEN
c$$$               DO NTH=1,NTHMAX
c$$$                  DO NHH=1,NHHMAX
c$$$                     CBP(K,NTH,NHH,NR)=(0.D0,0.D0)
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
      USE wmfem_comm, only: cpp,cpa,nthmax2,nhhmax2
      USE plprof,ONLY: pl_prof2
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION CPF1(nthmax2,nhhmax2),CPF2(nthmax2,nhhmax2)
C
      CW=2.D0*PI*CRF*1.D6
      DTH=2.D0*PI/DBLE(NTHMAX)
      DPH=2.D0*PI/DBLE(NHHMAX)
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
         DO NDX=1,nhhmax
            if(nhhmax.eq.1) then
               NDX1=NDX
            else
               NDX1=NDX+nhhmax/2-1
               IF(NDX1.gt.nhhmax) NDX1=NDX1-nhhmax
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
            DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PABS(NTH,NHH,NR,NS)=0.D0
            ENDDO
            ENDDO
            DO NDX=1,NDSIZ
            DO MDX=1,MDSIZ
               DO KKX=1,nhhmax2
               DO LLX=1,nthmax2
                  CPF1(LLX,KKX)=cpp(MDX,NDX,LLX,KKX,NR,NS)
               ENDDO
               ENDDO
               CALL WMSUBEX(CPF1,CPF2,NTHMAX2,NHHMAX2)
               DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PABS(NTH,NHH,NR,NS)=PABS(NTH,NHH,NR,NS)
     &                               +DBLE(CPF2(2*NTH-1,2*NHH-1))
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO

C     +++++ Antenna impedance +++++

      cradtt=(0.d0,0.d0)
      do nhh=1,nhhmax
         do nth=1,nthmax
            cradtt=cradtt+cpa(nth,nhh)
         end do
      end do

C
C     +++++ CALCULATE DRIVEN CURRENT IN REAL SPACE +++++
C
      ns=1
      DO nr=1,nrmax
         rho=xrho(nr)
         DO nhh=1,nhhmax
         DO nth=1,nthmax
            pcur(nth,nhh,nr)=0.D0
         ENDDO
         ENDDO

         CALL pl_prof2(rho,rn,rtpr,rtpp,ru)
         vte=SQRT(rtpr(1)*aee*1.D3/(pa(1)*amp))
         ww=DBLE(cw)
         IF(rn(1).LE.0.D0) THEN
            rlnlmd=15.D0
         ELSE
            rt=(rtpr(1)+2*rtpp(1))/3.D0
            rlnlmd=16.1D0 - 1.15D0*LOG10(rn(1)) + 2.30D0*LOG10(rt)
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
            CALL WMSUBEX(CPF1,CPF2,NTHMAX,NHHMAX)

            DO nhh=1,nhhmax
               ph=dph*(nhh-1)
            DO nth=1,nthmax
               th=dth*(nth-1)
               CALL wmfem_magnetic(rho,th,ph,babs,bsupth,bsupph)
               rkpara=mm*bsupth/babs+nn*bsupph/babs
               IF(ABS(rkpara).LT.1.D-8) rkpara=1.D-8

               IF(ABS(ww/rkpara).LT.VC) THEN
                  w=ww/(rkpara*vte)
                  xl=(rpst(nth,nhh,nr)-raxis)/rr
                  yl=(ZPST(NTH,NHH,NR)-zaxis)/RR
                  efcd=w1cdef(ABS(w),zeff,xl,yl,1)
                  IF(w.LT.0.D0) efcd=-efcd
                  IF (rn(1).GT.0.D0) THEN
                     pcur(nth,nhh,nr)=pcur(nth,nhh,nr)
     &                    +0.384D0*rtpr(1)*efcd
     &                    /(rn(1)*rlnlmd)*DBLE(cpf2(nth,nhh))
     &                    /(2.D0*PI*rpst(nth,nhh,nr))
                  END IF
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE CURRENT DRIVE EFFICIENCY ******
C
C      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
C      Z  = ZEFF   : EFFECTIVE Z
C      XL = X / RR : NORMALIZED X
C      YL = Y / RR : NORMALIZED Y
C      ID : 0 : LANDAU DAMPING
C           1 : TTMP
C
C      Reference: D.A. Ehst and C.F.F. Karney, 
C                 Nucl. Fusion, 31 (10) 1933-1938 (1991)
C
      FUNCTION W1CDEF(WT,Z,XL,YL,ID)
C
      INCLUDE 'wmcomm.inc'
C
      R=SQRT(XL*XL+YL*YL)
      IF(ID.EQ.0) THEN
         D=3.D0/Z
         QC=3.83D0
         A=0.D0
         RM=1.38D0
         RC=0.389D0
      ELSE
         D=11.91D0/(0.678D0+Z)
         QC=4.13D0
         A=12.3D0
         RM=2.48D0
         RC=0.0987D0
      ENDIF
      IF(WT.LE.1.D-20) THEN
         W=1.D-20
      ELSE
         W=WT
      ENDIF
      EFF0=D/W+QC/Z**0.707D0+4.D0*W*W/(5.D0+Z)
      EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)
C
      Y2=(R+XL)/(1.D0+R)
      IF(Y2.LT.0.D0) Y2=0.D0
      Y1=SQRT(Y2)
      EFF2=1.D0+A*(Y1/W)**3
C
      IF(Y2.LE.1.D-20) THEN
         YT=(1.D0-Y2)*WT*WT/1.D-60
      ELSE
         YT=(1.D0-Y2)*WT*WT/Y2
      ENDIF
      IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
         ARG=(RC*YT)**RM
         IF(ARG.LE.100.D0) THEN
            EFF3=1.D0-MIN(EXP(-ARG),1.D0)
         ELSE
            EFF3=1.D0
         ENDIF
      ELSE
         EFF3=1.D0
      ENDIF
C
      W1CDEF=EFF0*EFF1*EFF2*EFF3
C
      RETURN
      END
C
C     ****** 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBFX(CF1,CF2,NTHMAX,NHHMAX)
C
      implicit none
      integer,intent(in):: nthmax,nhhmax
      complex(8),dimension(nthmax,nhhmax),intent(in):: CF1
      complex(8),dimension(nthmax,nhhmax),intent(out):: CF2
      complex(8),dimension(nthmax):: CFM
      complex(8),dimension(nhhmax):: CFN
      integer:: NHH,NTH
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NHH)
         ENDDO
         CALL WMFFFT(CFM,NTHMAX,0)
         DO NTH=1,NTHMAX
            CF2(NTH,NHH)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(NTH,NHH)
         ENDDO
         CALL WMFFFT(CFN,NHHMAX,0)
         DO NHH=1,NHHMAX
            CF2(NTH,NHH)=CFN(NHH)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** INVERSE 2D FOURIER TRANSFORM ******
C
      SUBROUTINE WMSUBEX(CF1,CF2,NTHMAX,NHHMAX)
C
      implicit none
      integer,intent(in):: nthmax,nhhmax
      complex(8),dimension(nthmax,nhhmax),intent(in):: CF1
      complex(8),dimension(nthmax,nhhmax),intent(out):: CF2
      complex(8),dimension(nthmax):: CFM
      complex(8),dimension(nhhmax):: CFN
      integer:: NHH,NTH
C
      DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            CFM(NTH)=CF1(NTH,NHH)
         ENDDO
         CALL WMFFFT(CFM,NTHMAX,1)
         DO NTH=1,NTHMAX
            CF2(NTH,NHH)=CFM(NTH)
         ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         DO NHH=1,NHHMAX
            CFN(NHH)=CF2(NTH,NHH)
         ENDDO
         CALL WMFFFT(CFN,NHHMAX,1)
         DO NHH=1,NHHMAX
            CF2(NTH,NHH)=CFN(NHH)
         ENDDO
      ENDDO
      RETURN
      END
C
C     ****** INTERFACE FOR FFT ******
C
      SUBROUTINE WMFFFT(CA,N,KEY)
C
      USE libfft,ONLY: fft2l
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
