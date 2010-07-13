c
c=======================================================================
      subroutine nbout
c=======================================================================
      implicit none
      include 'AAA'
      include 'CNT'
      include 'TRN'
      include 'NBI'
! local variables
      integer n,is
      real*8 xnbie
c=======================================================================
      xnbie=0.d0
      do is=0,mion
      xnbie=xnbie+tnbie(is)
      enddo
c-----
      write(6,'(5(a7,1pd10.3))')'  time=',time,' tnbip=',tpnbi
     >    ,' tnbie=',xnbie,' tnbic=',tnbic
      write(66,'(5(a7,1pd10.3))')'  time=',time,' tnbip=',tpnbi
     >    ,' tnbie=',xnbie,' tnbic=',tnbic
      do n=1,nro
      write(66,'(1p10d10.3)')den(n,0),tem(n,0)
     >    ,prnbi(n),snbic(n),(snbie(n,is),is=0,nion)
      enddo
      return
      end
