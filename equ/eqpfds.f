c
      module eqpfds_mod
      use tpxssl_mod
      public
      contains
c
c=======================================================================
      subroutine eqpfds(pds0,fds0)
c=======================================================================
c     load pds0 and fds0
c=======================================================================
      use aaa_mod
      use par_mod
      use eqv_mod
      use com_mod
      implicit none
! argument
      real*8    pds0(ivdm), fds0(ivdm)
!local variables
      integer   n 
      real*8    frs0(ivdm), prs(ivdm), prs0(ivdm), rho0(ivdm), x
c=======================================================================
      do n=1,nv
      x=dfloat(n-1)/dfloat(nvm)
      rho(n)=dsqrt(vlv(n)/vlv(nv))
      rho0(n)=x
      prs0(n)=(1.-cp(4))*(1.-x**cp(2))**cp(3)+cp(4)*(1.-x**2)**cp(5)
      frs0(n)=-((1.-cp(8))*(1.-x**cp(6))**cp(7)+cp(8))
      enddo
c-----
      call spln(prs,rho,nv,prs0,rho0,nv,ww1,0)
      call spln(fds0,rho,nv,frs0,rho0,nv,ww1,0)
c-----
      do n=2,nvm
      pds0(n)=(prs(n+1)-prs(n-1))/(siw(n+1)-siw(n-1))
      enddo
      pds0(1)=2.*pds0(2)-pds0(3)
      pds0(nv)=2.*pds0(nv-1)-pds0(nv-2)
c-----
      return
      end subroutine eqpfds
c
      end module eqpfds_mod
