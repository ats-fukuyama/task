c
      module eqset_mod
      use tpxssl_mod
      use eqsub_mod
      use eqfct_mod
      public
      contains
c=======================================================================
      subroutine eqset
c=======================================================================
c     set up free bondary equilibriium
c             for fiven functional forms of dp/ds and TdT/ds
c=======================================================================
      use aaa_mod
      use geo_mod
      use par_mod
      use equ_mod
      use eqv_mod
      use cnt_mod
      use com_mod
      implicit none
! local variables
      real*8   bav0
      integer keqmax,i
c=======================================================================
c<<restart>>
      if(iread.gt.0)then
      do i=1,nsr*nsz
      rbp(i)=psi(i)
      enddo
      bav0=bav
      bav=1.
      call eqvac
      do i=1,nsr*nsz
      psi(i)=rbp(i)
      enddo
      call eqrcu
      bav=bav0
c--confirm equilibrium accuracy for metrics
      write(6,*)
      write(6,*)'== check of equilibrium accuracy through ode =='
      call eqeqv
      keqmax=ieqmax
      ieqmax=10
      call eqfct
      ieqmax=keqmax
      return
      endif
c=======================================================================
      call cpusec(cputim)
c-----------------------------------------------------------------------
c     load Solovev equilibrium as initial guess
c-----------------------------------------------------------------------
      call eqvac
      call eqsol
      call eqrbp
      call eqlin
      call eqequ
      call eqeqv
      call eqout
c--confirm equilibrium accuracy for metrics
      write(6,*)
      write(6,*)'== check of equilibrium accuracy through ode =='
      keqmax=ieqmax
      ieqmax=10
      call eqfct
      ieqmax=keqmax
      call cpusec(cpuequ)
c=======================================================================
      return
      end subroutine eqset
c
      end module eqset_mod
