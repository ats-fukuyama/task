c
      module eqfct_mod
      use eqsub_mod
      use tpxssl_mod
      public
      contains
c=======================================================================
      subroutine eqfct
c=======================================================================
c     MAIN PROGRAM                                                 JAERI
c              OF FREE BOUNDARY FCT EQUILIBRIUM SOLVER
c=======================================================================
      use aaa_mod
      use equ_mod
      use geo_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, ia, ie
      real*8    b2, b3, dv1, dv2, psistr(irzdm), v1, v2, v3
c=======================================================================
c step > 1 :  go
c-----------------------------------------------------------------------
      if(ieqmax.le.0)return
      call cpusec(cputim)
c-----------------------------------------------------------------------
      call eqode
c-----------------------------------------------------------------------
      v1=0.
      v2=0.
      v3=0.
      b2=0.
      b3=0.
      if(bav.gt.bavmax)bav=bavmax
c-----
      ie=0
      ieq=1
      do while(ie.lt.ieqmax.and.ieq.gt.0)
      ie=ie+1
      ieq=ie
      call eqrcu
      call eqbnd
      call eqpde
      do i=1,nrz
      psistr(i)=psi(i)
      enddo
c-----
      ia=0
      iad=1
      do while(ia.lt.iadmax.and.iad.gt.0)
      ia=ia+1
      iad=ia
      do i=1,nrz
      psi(i)=psistr(i)
      enddo
      call eqsep
      call eqadj
      if(iad.gt.0)then 
      call eqaxi
      call eqsep
      call eqrbp
      call eqlin
      call eqode
      endif
      enddo
c-----
      call eqchk
c-----
      v1=v2
      v2=v3
      v3=vlv(nv)
      b2=b3
      b3=bav
      if(ieq.gt.2)then
c-----
      if((v1-v2)*(v2-v3).gt.0.)then
      if(b2.eq.b3)bav=2.*b3
      if(bav.gt.bavmax)bav=bavmax
      else
      dv1=dabs(v1-v2)
      dv2=dabs(v2-v3)
      if(dv2.gt..5*dv1)bav=.5*bav
      if(bav.lt.bavmin)bav=bavmin
      endif
      endif
      enddo
c-----------------------------------------------------------------------
      call cpusec(cpuequ)
c=======================================================================
      return
      end subroutine eqfct
c    
c=======================================================================
      subroutine eqode
c=======================================================================
c     solve ode
c=======================================================================
      use aaa_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      implicit none
! local variables
      integer   n
      real*8    aag(ivdm), akn(ivdm), amg(ivdm), d(ivdm), dakn(ivdm)
     >        , damg(ivdm), e, er, f(ivdm), r(ivdm), rbg(ivdm), rbt
     >        , x, y
c-----------------------------------------------------------------------
      do n=1,nv
      aag(n)=aav(n)**gam
      enddo
      iod=0
      errod=1.d0
      do while(iod.lt.iodmax.and.errod.gt.eodmax)
      iod=iod+1

c%%%%%     call spln(muv,hiv,nv,xxhit,zzhit,nnhit,wwmu,1)
c%%%%%     call spln(nuv,hiv,nv,yyhit,zzhit,nnhit,wwnu,1)

      do n=1,nv
      akn(n)=nuv(n)*aav(n)*ckv(n)
      amg(n)=muv(n)*aag(n)
      rbg(n)=(hdv(n)/aav(n))**gam2
      enddo
      call deriv(dakn,akn,siw,nv)
      call deriv(damg,amg,siw,nv)
      d(1)=(nuv(1)*aav(1)*dakn(1)+rbg(1)*damg(1))
     >    /(nuv(1)*aav(1)*akn(1)+gam*amg(1)*rbg(1)+aav(1))
      r(1)=aav(1)
      f(1)=0.
      e=0.
      do n=2,nv
      d(n)=(nuv(n)*aav(n)*dakn(n)+rbg(n)*damg(n))
     >    /(nuv(n)*aav(n)*akn(n)+gam*amg(n)*rbg(n)+aav(n))
      e=e-0.5*(d(n)+d(n-1))*(siw(n)-siw(n-1))
      r(n)=aav(n)*dexp(e)
      f(n)=f(n-1)+0.5*(r(n)+r(n-1))*(vlv(n)-vlv(n-1))
      enddo
      errod=0.
      x=hiv(nv)/f(nv)
      do n=1,nv
      y=x*r(n)
      er=dabs(1.-y/hdv(n))
      if(er.gt.errod)errod=er
      hdv(n)=y
      hiv(n)=x*f(n)
      enddo
      enddo
      bts=hdv(nv)/aav(nv)
      do n=nv,1,-1
      sdv(n)=nuv(n)*hdv(n)
      rbt=hdv(n)/aav(n)
      pds(n)=(-gam*amg(n)*d(n)+damg(n))*rbt**gam
      fds(n)=-d(n)*rbt*rbt
      if(n.ne.nv)siv(n)=siv(n+1)+0.5*(nuv(n)+nuv(n+1))*(hiv(n)-hiv(n+1))
      enddo
c..
      cpl=ckv(nv)*sdv(nv)/cnpi2
c..
      return
      end subroutine eqode
c 
c=======================================================================
      subroutine eqchk
c=======================================================================
c     convergence check
c=======================================================================
      use aaa_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, j, k, n
      real*8    e, sm
C      real*8    sav1(ivdm,4), sav2(ivdm,2), sav3(ivdm,5), sav4(ivdm,2)
c=======================================================================
C      equivalence(hiv(1),sav1(1,1))
C      equivalence(siw(1),sav2(1,1))
C      equivalence(vlv(1),sav3(1,1))
C      equivalence(pds(1),sav4(1,1))
c=======================================================================
      erreq=0.
      do j=1,13
      erch(j)=0.
      enddo
c-----------------------------------------------------------------------
      call eqchk_sub(hiv,nv,1)
      call eqchk_sub(hdv,nv,2)
      call eqchk_sub(siv,nv,3)
      call eqchk_sub(sdv,nv,4)
C
      call eqchk_sub(siw,nv,5)
      call eqchk_sub(sdw,nv,6)
C
      call eqchk_sub(vlv,nv,7)
      call eqchk_sub(ckv,nv,8)
      call eqchk_sub(ssv,nv,9)
      call eqchk_sub(aav,nv,10)
      call eqchk_sub(rrv,nv,11)
C
      call eqchk_sub(pds,nv,12)
      call eqchk_sub(fds,nv,13)
c-----------------------------------------------------------------------
c     i=-dlog10(dabs(errcu))+0.9999
c     e=0.
c     if(i.gt.0)e=i+0.1*dabs(errcu)*(10.**i)
c     if(errcu.lt.0.)e=-e
c     errcu=e
      do j=1,13
      i=-dlog10(dabs(erch(j)))+0.9999
      e=0.
      if(i.gt.0)e=i+0.1*dabs(erch(j))*(10.**i)
      if(erch(j).lt.0.)e=-e
      erch(j)=e
      enddo
c-----------------------------------------------------------------------
      ccpl=cpl/cnmu0
      ttcu=tcu/cnmu0
      if(ieqout(2).ge.2)
     >write(ft06,'(i3,2f6.3,4f8.3,14f5.1)')
     >        ieq,bav,raxis,saxis,vlv(nv),ccpl,ttcu
     >       ,(erch(k),k=1,13)
c-----------------------------------------------------------------------
      if(erreq.lt.eeqmax)ieq=-ieq
      if(ieq.gt.ieqmax)ieq=0
      return
      end subroutine eqchk
c 
c=======================================================================
      subroutine eqchk_sub(data,nv,j)
c=======================================================================
c     convergence check
c=======================================================================
      use aaa_mod
      use equ_mod
      implicit none
! local variables
      integer   j, k, n, nv
      real*8    e, sm
      real*8    data(ivdm)
c=======================================================================
      sm=0.
      do n=1,nv
      if(dabs(data(n)).gt.sm)sm=dabs(data(n))
      enddo
      if(sm.le.0.)sm=1.
      do n=1,nv
      e=(data(n)-save(n,j))/sm/bav
      save(n,j)=data(n)
      if(dabs(e).gt.dabs(erch(j)))erch(j)=e
      if(dabs(e).gt.erreq)erreq=dabs(e)
      enddo
c-----
      return
      end subroutine eqchk_sub
c
      end module eqfct_mod
