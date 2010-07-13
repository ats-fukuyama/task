c
c=======================================================================
      subroutine trsvq
c=======================================================================
c     TEMPORAL EVOLUTION OF CURRENT PROFILE
c
c     SOLVE FARADAYS LAW                                           JAERI
c=======================================================================
      implicit none
      include 'AAA'
      include 'PAR'
      include 'EQT'
      include 'TRN'
      include 'CNT'
! local variables
      integer   m, n
      real*8    a(itdm), a0, a1, b(itdm), b0, b1, c(itdm), c0, c1, cext
     >        , cfac, cpal, ctor, d(itdm), dkah(0:itdm), dro2(0:itdm)
     >        , facx, fkah(0:itdm), gcur(0:itdm), rez, u, x
     >        , xqi0(0:itdm), y, z, zqimin
     >        ,eta(itdm)
c=======================================================================
      cfac=1.0
c=======================================================================
c-----BOUNDARY CONDITION
c-----------------------------------------------------------------------
c::ICQI(1)>0   FIXED Q-PROFILE : FCT SCHEME
c::ICQI(1)=0   TOTAL CURRENT IS GIVEN : IP=TCUR(MA)
c-----------------------------------------------------------------------
      qi(nt)=2.*cnpi*(cnmu0*tcur)/(ckt(nt)*hdt(nt))
c=======================================================================
c::ICQI(2)=1  STEADY STATE SOLUTION OF Q-PROFILE
c-----------------------------------------------------------------------
      zqimin=0.
 101  continue
      x=0.
      y=0.
      do 100 n=1,ntm
      eta(n)=1./(0.01+(1.-roh(n)**2)**2)**1.5
      u=(ro(n+1)-ro(n))/(eta(n)/cnmu*roh(n)/(aah(n)*vrh(n))**2)
      x=x+u
 100  continue
      facx=(aat(nt)*ckt(nt)*qi(nt)-y)/x
      x=0.
      y=0.
      do 110 n=ntm,2,-1
      u=(ro(n+1)-ro(n))/(eta(n)/cnmu*roh(n)/(aah(n)*vrh(n))**2)
      x=x+u
      z=facx*x+y
      qi(n)=(aat(nt)*ckt(nt)*qi(nt)-z)/(aat(n)*ckt(n))
      if(qi(n).lt.zqimin)zqimin=qi(n)
 110  continue
      qi(1)=(4.*qi(2)-qi(3))/3.
      if(qi(1).lt.zqimin)zqimin=qi(1)
c-----
      if(zqimin.lt.0.)then
      write(6,fmt='(20X,A35)')'>>>>> CUATION : QAXIS < 0. AT TRSVQ'
      if(cfac.eq.1.)then
      cfac=0.
      write(6,fmt='(20X,A35)')'                CFAC = 0.  AT TRSVQ'
      goto 101
      endif
      endif
c--
      do 120 n=1,nro
      qic(n)=qi(n)
 120  qi0(n)=qi(n)
c-----------------------------------------------------------------------
      do 300 n=1,nt
 300  if(qi(n).le.0.)goto 310
      return
 310  continue
      do 380 n=1,nt
      write(6,*) 'N=',n,'QI=',qi(n)
 380  continue
      stop'=====  ERROR STOP AT TRSVQ  ====='
      end
