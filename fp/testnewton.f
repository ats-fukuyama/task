C
      IMPLICIT NONE
      INTEGER,PARAMETER:: nthlm=101
      INTEGER:: nthlmax,nthl,nthlmaxl
      REAL(8):: thlmin,thlmax,dthl,thl,thetal,eps
      REAL(8):: rfunc,rfuncp,dfunc,dfuncp
      REAL(8),DIMENSION(nthlm):: gx
      REAL(8),DIMENSION(nthlm,4):: gy0,gy1,gy2

      CALL GSOPEN

      nthlmax=81
      nthlmaxl=61
      thlmin=-6
      thlmax=+2
      dthl=(thlmax-thlmin)/(nthlmax-1)
      eps=1.d-8
      
      DO nthl=1,nthlmax
         thl=thlmin+dthl*(nthl-1)
         thetal=10.d0**thl

         gx(nthl)=thl
         gy0(nthl,1)=LOG10(rfunc(thetal))
         gy0(nthl,2)=LOG10(rfuncp(thetal))
         gy1(nthl,1)=rfunc(thetal)/thetal
         gy1(nthl,2)=rfuncp(thetal)/thetal
         gy2(nthl,1)=dfunc(thetal)
         gy2(nthl,2)=dfuncp(thetal)
         gy2(nthl,3)=(rfunc(thetal+eps)-rfunc(thetal-eps))/(2.d0*eps)
         gy2(nthl,4)=(rfuncp(thetal+eps)-rfuncp(thetal-eps))/(2.d0*eps)
      ENDDO

      CALL PAGES
      CALL GRD1D(1,gx,gy0,nthlm,nthlmax,2,'@LOG(rfunc)@',3)
      CALL PAGEE

      CALL PAGES
      CALL GRD1D(1,gx,gy1,nthlm,nthlmax,2,'@rfunc*z@',1)
      CALL PAGEE

      CALL PAGES
      CALL GRD1D(1,gx,gy2,nthlm,nthlmax,4,'@dfunc@',1)
      CALL PAGEE

      CALL GSCLOS
      STOP

      CONTAINS

      FUNCTION rfunc(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,rfunc
      REAL(8):: z,dkbsl1,dkbsl2,BESEKN
      z=1.D0/thetal
      dkbsl1=BESEKN(1,Z)
      dkbsl2=BESEKN(2,Z)
      rfunc= dkbsl1 /dkbsl2 -1.D0+3.D0/Z
      RETURN
      END FUNCTION rfunc

      FUNCTION rfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,rfuncp
      REAL(8):: z,dkbsl1,dkbsl2
      z=1.D0/thetal
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      rfuncp= dkbsl1 /dkbsl2 -1.D0+3.D0/Z
      RETURN
      END FUNCTION rfuncp
      
      FUNCTION dfunc(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,dfunc
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3,BESEKN
      z=1.D0/thetal
      dkbsl0=BESEKN(0,z)
      dkbsl1=BESEKN(1,z)
      dkbsl2=BESEKN(2,z)
      dkbsl3=BESEKN(3,z)
      dfunc =( (dkbsl0 +dkbsl2 )/dkbsl2
     &        -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2
     &      +3.d0
      RETURN
      END FUNCTION dfunc

      FUNCTION dfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,dfuncp
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3
      z=1.D0/thetal
      dkbsl0=1.D0 -  1.D0/8.D0/z +   9.D0/128.D0/z**2
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      dkbsl3=1.D0 + 35.D0/8.D0/z + 945.D0/128.D0/z**2
      dfuncp =( (dkbsl0 +dkbsl2 )/dkbsl2
     &         -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2
     &      +3.d0
      RETURN
      END FUNCTION dfuncp

      END






