module fprad

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE fp_radiation
!           use fpsave
!           use fpcomm
!           use libgrf
!           USE libbes
!
!           implicit none
!
!
! !for bermsstralung
!           double precision :: h,g,alpha1,beta,pv,sumr,rsum1,rsum2!,eps00
!           double precision,dimension(nrmax) :: pbre,pave,pave_perp
!           double precision,dimension(nrmax,nsamax) :: pbr,rhot,temps
!           integer :: nsa,nr,ns,np,nth
!
! !for cyclotron
!           double precision,dimension(nrmax,nsamax)::pcyc
!           double precision,dimension(nrmax)::pcyce
!           double precision,dimension(nthmax,npmax,nrmax,nsamax)::plcy
!           double precision,dimension(:,:,:),allocatable::sump
!           double precision,dimension(:,:,:,:),allocatable::pcycle
!           double precision :: b, bpara,bperp,omega,omega0,di,x,jint,jdeff,alpha2
!           integer::l,i,imax,j,j2
!
!
!
!           !eps00=eps0*aee*1.d-3 !in keV/m
!           !h=4.135667662d-18 !in keV*sec
!           h=6.626070040d-34 !in J*sec
!           g=2.d0*sqrt(3.d0)/pi
!
!             do nr=1,nrmax
!               rsum1=0.d0
!               rsum2=0.d0
!
!               do np=1,npmax
!                 beta=SQRT(1.D0+THETA0(1)*PM(NP,1)**2)
!                 do nth=1,nthmax
!                   RSUM1 = RSUM1+VOLP(NTH,NP,1)*fnsp(NTH,NP,NR,1)*PM(NP,1)*sinm(nth)/beta &
!                        *(1.D0+EPSRM2(NR))
!
!                   RSUM2 = RSUM2+VOLP(NTH,NP,1)*fnsp(NTH,NP,NR,1)*PM(NP,1)/beta &
!                        *(1.D0+EPSRM2(NR))
!                  end do
!                end do
!                pave_perp(nr)=rsum1*ptfp0(1)/rnsl(nr,1)
!                pave(nr)=rsum2*ptfp0(1)/rnsl(nr,1)
!              end do
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           ! bremsstralung
!
!             alpha1=g*aee**6.d0/(6.d0*sqrt(1.5d0)*pi**1.5d0*eps0**3.d0*vc**3.d0*h*ame**1.5d0)
!
!             ! do nsa=1,nsamax
!             !   ns=ns_nsa(nsa)
!             !   do nr=1,nrmax
!             !     sumr=0.d0
!             !     do np=1,npmax
!             !       pv=(1.d0+theta0(ns)*pm(np,ns)**2.d0)
!             !       do nth=1,nthmax
!             !         sumr=sumr+volp(nth,np,nsa)*fnsp(nth,np,nr,nsa) &
!             !               *pa(ns)*amp*(pv-1.d0)*vc**2.d0*RLAMDA(NTH,NR)*RFSADG(NR)
!             !       end do
!             !     end do
!             !     rhot(nr,nsa)=sumr/aee*1.d-3! *1d20 keV/m^3
!             !     temps(nr,nsa)=2.d0*rhot(nr,nsa)/(3.d0*rnsl(nr,nsa))!keV
!             !   end do
!             ! end do
!
!
!           do nsa=1,nsamax
!             ns=ns_nsa(nsa)
!             do nr=1,nrmax
!               if(ns.eq.1) then
!                 pbre(nr)=0.d0
!               else
!                 pbre(nr)=pbre(nr) &
!                         -alpha1*pz(ns)**2.d0*rnsl(nr,1)*rnsl(nr,nsa)*1.d40*pave(nr)/ame*(ame/3.d0)**0.5
!                 ! pbre(nr)=pbre(nr)-5.34d-37*pz(ns)**2.d0*rnsl(nr,nsa)*rnsl(nr,1)*1.d40*temps(nr,1)**0.5d0
!               end if
!             end do
!           end do
!           do nsa=1,nsamax
!             ns=ns_nsa(nsa)
!             do nr=1,nrmax
!               if(ns.eq.1) then
!                 pbr(nr,nsa)=pbre(nr)*1.d-6
!               else
!                 pbr(nr,nsa)=0.d0
!               end if
!             end do
!           end do
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! cyclotron radiation
!           if(modelr.eq.0)then
!             alpha2=aee**2.d0/(6.d0*pi*eps0*vc**3)
!             do nr=1,nrmax
!               pcyce(nr)=alpha2*rnsl(nr,1)*1.d20*(aee*bb/ame**2.d0*pave_perp(nr))**2.d0*1.d-3 !kW/m^3
!             end do
!             do nsa=1,nsamax
!               ns=ns_nsa(nsa)
!               do nr=1,nrmax
!                 if(ns.eq.1)then
!                   pcyc(nr,nsa)=-pcyce(nr)
!                 else
!                   pcyc(nr,nsa)=0.d0
!                 end if
!               end do
!             end do
!
!           else !relativestic
!             omega0=aee*bb/ame
!             alpha2=aee**2.d0/(2.d0*pi*eps0*vc)
!             l=4
!             allocate(pcycle(nthmax,npmax,nrmax,l))
!             allocate(sump(nthmax,npmax,nrmax))
!             pcyce=0.d0
!
!             do nr=1,nrmax
!               do np=1,npmax
!                 do nth=1,nthmax
!                   b=(1.d0-1.d0/(1+theta0(1)*pm(np,1)**2.d0))**0.5d0
!                   bpara=(1.d0-1.d0/(1+theta0(1)*pm(np,1)**2.d0*cosm(nth)**2.d0))**0.5d0
!                   bperp=(1.d0-1.d0/(1+theta0(1)*pm(np,1)**2.d0*sinm(nth)**2.d0))**0.5d0
!                   imax=100
!                   sump(nth,np,nr)=0.d0
!                   do j=1,l
!                     j2=2*j
!                     x=2.d0*j*bperp/(1-bpara**2.d0)**0.5d0
!                     di=x/imax
!                     jdeff=(BESJN(j2-1,x)+BESJN(j2+1,x))*0.5d0
!                     jint=0.d0
!                     do i=1,imax-1
!                       jint=jint+di*BESJN(j2,(i-0.5d0)*di)/j2
!                     end do
!                     pcycle(nth,np,nr,j)=alpha2*omega0**2.d0 &
!                                         *(1-b**2.d0)/bperp/(1.d0-bpara**2)**1.5d0 &
!                                         *(j*bperp**2.d0*jdeff &
!                                         -j**2.d0*(1-b**2.d0)*jint)
!                     sump(nth,np,nr)=sump(nth,np,nr)+pcycle(nth,np,nr,j)
!                   end do
!                   ! if(nr.eq.1)then
!                   !   pcyce(nr)=volp(nth,np,1)*sump(nth,np,nr)*fnsp(nth,np,nr,1)*1.d20*RLAMDA(NTH,NR)*RFSADG(NR)
!                   ! else
!                     pcyce(nr)=pcyce(nr)+volp(nth,np,1)*sump(nth,np,nr)*fnsp(nth,np,nr,1)*1.d20*RLAMDA(NTH,NR)*RFSADG(NR)
!                   ! end if
!                 end do
!               end do
!             end do
!
!             do nsa=1,nsamax
!               ns=ns_nsa(nsa)
!                 if(ns.eq.1)then
!                   do nr=1,nrmax
!                     pcyc(nr,nsa)=-pcyce(nr)*1.d-6
!                     do np=1,npmax
!                       do nth=1,nthmax
!                         plcy(nth,np,nr,nsa)=-sump(nth,np,nr)*1.d-3/aee/rtfd0(nsa)
!                       end do
!                     end do
!                   end do
!                 else
!                     do nr=1,nrmax
!                       pcyc(nr,nsa)=0.d0
!                       do np=1,npmax
!                         do nth=1,nthmax
!                           plcy(nth,np,nr,nsa)=0.d0
!                         end do
!                       end do
!                     end do
!                 end if
!               end do
!             end if
!
!
!           call pages
!             CALL GRD1D(1,RM,pbr,  NRMAX,  NRMAX,  NSAMAX,'@brem MW/m-3@')
!             CALL GRD1D(2,RM,pcyc,  NRMAX,  NRMAX,  NSAMAX,'@cyclotron MW/m-3@')
!           call pagee
!
!           deallocate(pcycle,sump)
!
!           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           !       SUBROUTINE synchrotron
!           !
!           !       IMPLICIT NONE
!           !       real(8):: alpha, rgama_para, u, rgama
!           !       integer:: NSA, NTH, NP, NR, NS
!           !
!           !       DO NSA=NSASTART,NSAEND
!           !          NS=NS_NSA(NSA)
!           !       DO NR=NRSTART, NREND
!           !
!           !          DO NP=NPSTART, NPENDWG
!           !             rgama=SQRT(1.D0+THETA0(NS)*PG(NP,NSA)**2)
!           !             alpha=(2.D0*AEE**4)/(3*AMFP(NSA)**3*VC**5*RGAMA)*1.e30
!           !             u=PTFP0(NSA)*PG(NP,NSA)/AMFP(NSA)
!           !             DO NTH=1, NTHMAX
!           !                rgama_para=SQRT(1.D0+THETA0(NS)*PG(NP,NSA)**2*COSM(NTH)**2)
!           !
!           !                FSPP(NTH,NP,NR,NSA)= -&
!           !                     alpha*BB**2*rgama_para**2*U*SINM(NTH)**2* &
!           !                     (1.D0+ (U/(rgama_para*VC))**2*COSM(NTH)**2 ) &
!           !                     *AMFP(NSA)/PTFP0(NSA)
!           !             END DO
!           !          END DO
!           !
!           !          DO NP=NPSTARTW, NPENDWM
!           !             rgama=SQRT(1.D0+THETA0(NS)*PM(NP,NSA)**2)
!           !             alpha=(2.D0*AEE**4)/(3*AMFP(NSA)**3*VC**5*RGAMA)*1.e30
!           !             u=PTFP0(NSA)*PM(NP,NSA)/AMFP(NSA)
!           !             DO NTH=1, NTHMAX+1
!           !                rgama_para=SQRT(1.D0+THETA0(NS)*PM(NP,NSA)**2*COSG(NTH)**2)
!           !
!           !                FSTH(NTH,NP,NR,NSA)= -&
!           !                     alpha*BB**2*rgama_para**2*U*SING(NTH)*COSG(NTH)* &
!           !                     (1.D0-U**2/(rgama_para*VC)**2*SING(NTH)**2 ) &
!           !                     *AMFP(NSA)/PTFP0(NSA)
!           !             END DO
!           !          END DO
!           !
!           !       END DO
!           !       END DO
!           !
!           !       END SUBROUTINE synchrotron
!           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        END SUBROUTINE fp_radiation
end module fprad
