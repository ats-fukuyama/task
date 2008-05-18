!======================================================================c

      subroutine  rkhn(y,x,dt)

!======================================================================c

  USE mcnmod, ONLY : ntest

  IMPLICIT NONE
  INTEGER(4),PARAMETER:: ndim=4
  REAL(8),INTENT(INOUT):: y(ntest,ndim), x
  REAL(8),INTENT(IN)   :: dt
  INTEGER(4):: i, ipt
  REAL(8):: c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, &
    c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27
  REAL(8):: f1(ntest,ndim), f2(ntest,ndim), f3(ntest,ndim), f4(ntest,ndim),       &
    f5(ntest,ndim), f6(ntest,ndim), f7(ntest,ndim), f8(ntest,ndim),  f9(ntest,ndim)


      c01 = dt/9.d0
      c02 = dt*.4166666666666667d-1
      c03 = dt*.125d0
      c04 = dt/6.d0
      c05 = dt*0.5d0
      c06 = dt*.6666666666666667d0
      c07 = dt/3.d0
      c08 = dt*.375d0
      c09 = dt*.1333333333333333d1
      c10 = dt*.3333333333333333d1
      c11 = dt*.7d1
      c12 = dt*.9666666666666667d1
      c13 = dt*.1533333333333333d2
      c14 = dt*.6111111111111111d0
      c15 = dt*.1166666666666667d1
      c16 = dt*.1375d1
      c17 = dt*.8333333333333333d0
      c18 = dt*.4390243902439024d0
      c19 = dt*.8780487804878049d0
      c20 = dt*.1304878048780488d1
      c21 = dt*.2097560975609756d1
      c22 = dt*.2963414634146341d1
      c23 = dt*.4317073170731707d1
      c24 = dt*.3214285714285714d-1
      c25 = dt*.4880952380952381d-1
      c26 = dt*.2571428571428571d0
      c27 = dt*.3238095238095238d0

!         ***** main *****

      DO i = 1,ndim
        DO ipt = 1,ntest
          f1(ipt,i) = y(ipt,i)
        END DO
      END DO

!-----------------------------c
!         start rkhn
!-----------------------------c
!          1st step
!-----------------------------c

      CALL subfn( x,f1,f3)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i) = c01*f3(ipt,i) +f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!          2nd step
!-----------------------------c

      CALL subfn(x+c01,f2,f4)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i) = c02*f3(ipt,i) +c03*f4(ipt,i) +f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!          3rd step
!-----------------------------c

      CALL subfn(x+c04,f2,f5)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i) = c04*f3(ipt,i) -c05*f4(ipt,i) +c06*f5(ipt,i) +  f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!          4th step
!-----------------------------c

      CALL subfn(x+c07,f2,f6)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i) = c03*f3(ipt,i) +c08*f6(ipt,i) +f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!          5th step
!-----------------------------c

      CALL subfn(x+c05,f2,f7)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i) = -c09*f3(ipt,i) +c10*f7(ipt,i) -c11*f4(ipt,i) -c12*f6(ipt,i) &
                      +c13*f5(ipt,i) +    f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!          6th step
!-----------------------------c

      CALL subfn(x+c06,f2,f8)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i) = -c01*f3(ipt,i) +c03*f8(ipt,i) +c14*f7(ipt,i) -c15*f5(ipt,i) &
                      +c16*f4(ipt,i) +    f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!          7th step
!-----------------------------c

      CALL subfn(x+c17,f2,f9)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f2(ipt,i)= -c18*f8(ipt,i) +c19*f9(ipt,i) +c20*f3(ipt,i) -c21*f7(ipt,i) &
                     -c22*f4(ipt,i) +c23*f6(ipt,i) +f1(ipt,i)
        END DO
      END DO

!-----------------------------c
!         last step
!-----------------------------c

      CALL subfn(x+dt,f2,f4)

      DO i = 1,ndim
        DO ipt = 1,ntest
          f1(ipt,i) = (f6(ipt,i) +f8(ipt,i))*c24 +(f3(ipt,i) +f4(ipt,i))*c25   &
                     +(f5(ipt,i) +f9(ipt,i))*c26 +f7(ipt,i)*c27 +f1(ipt,i)
          y(ipt,i) = f1(ipt,i)
        END DO
      END DO


!-----------------------------c
!         end rkhn
!-----------------------------c


      x = x + dt

      RETURN
      END

!======================================================================c
      SUBROUTINE subfn(x,y,yp)
!======================================================================c

  USE mcnmod, ONLY : adinv, bco, c1bf, c1et, c1g, c1i, c2bf, c2et, c2g, c2i, c3bf,&
    c3et, c3g, c3i, cm, cn, cpi, cug, cui, eot, esp0, isw, kmsh, mmx, ntest, psi, &
    psia, tchg, tmass
  IMPLICIT NONE
  REAL(8),INTENT(INOUT):: x, y(ntest,4),yp(ntest,4)
  INTEGER(4):: i, ip, jpip
  INTEGER(4):: ipsw(ntest), jp(ntest)
  REAL(8):: cxip, sxip, th1ip
  REAL(8):: b(ntest), b1(ntest), c1(ntest), c21(ntest), c22(ntest), c3(ntest),     &
    cug1(ntest), cui1(ntest), dcug1(ntest), dcui1(ntest),db(ntest), db1(ntest),    &
    dbfi(ntest), dbth(ntest), def(ntest), delt(ntest), ds(ntest), eot1(ntest),     &
    g1(ntest), g2(ntest), gamm(ntest), sj(ntest), sp(ntest)

!CDIR NODEP
      DO ip=1, ntest
        IF (y(ip,1).lt.0.d0) THEN
          y(ip,1) = -y(ip,1)
          y(ip,2) =  y(ip,2) +2.d0*cpi
          jp(ip) = 1
          ds(ip) = 0.d0
        ELSE IF(y(ip,1).ge.psia) THEN
          sp(ip) = sqrt(psia)
          jp(ip) = kmsh +1
          ds(ip) = 0.d0
        ELSE
          sp(ip) = sqrt(y(ip,1))
          jp(ip) = int( y(ip,1)*float(kmsh)/psia ) +1
          sj(ip) = sqrt(psi(jp(ip)))
          ds(ip) = sp(ip) -sj(ip)
        END IF
      END DO
!
!CDIR NODEP
      DO ip=1, ntest
        jpip = jp(ip)
!
         eot1(ip)  = eot(jpip) +(c1et(jpip) +(c2et(jpip) &
                    +c3et(jpip)*ds(ip))*ds(ip))*ds(ip)
         cui1(ip)  = cui(jpip) +( c1i(jpip) +( c2i(jpip) &
                    +c3i(jpip)*ds(ip))*ds(ip))*ds(ip)
         cug1(ip)  = cug(jpip) +( c1g(jpip) +( c2g(jpip) &
                    +c3g(jpip)*ds(ip))*ds(ip))*ds(ip)
!
         dcui1(ip) = ( c1i(jpip) +( 2.d0*c2i(jpip) &
                    + 3.d0*c3i(jpip)*ds(ip))*ds(ip) )/(2.d0*sp(ip))
         dcug1(ip) = ( c1g(jpip) +( 2.d0*c2g(jpip) &
                    + 3.d0*c3g(jpip)*ds(ip))*ds(ip) )/(2.d0*sp(ip))
      END DO
!
!CDIR NODEP
      DO ip=1, ntest
        ipsw(ip) = int((isw(ip)+1.d0)/2.d0)
        b(ip) = 0.d0
        db(ip) = 0.d0
        dbth(ip) = 0.d0
        dbfi(ip) = 0.d0
      END DO

      DO i =1, mmx+1
!CDIR NODEP
        DO ip=1, ntest
          th1ip  = cm(i)*y(ip,2) -cn(i)*y(ip,3)
          cxip = cos(th1ip)
          sxip = sin(th1ip)

          jpip   = jp(ip)

          b1(ip)  = bco(jpip,i) +(c1bf(jpip,i) +(c2bf(jpip,i)      &
                   +c3bf(jpip,i)*ds(ip))*ds(ip))*ds(ip)
          db1(ip) = ( c1bf(jpip,i) +(2.d0*c2bf(jpip,i)             &
                   +3.d0*c3bf(jpip,i)*ds(ip))*ds(ip) )/(2.d0*sp(ip))
!
          b(ip) =    b(ip) + b1(ip)*cxip
          db(ip) =   db(ip) +db1(ip)*cxip

          dbth(ip) = dbth(ip) -cm(i)*b1(ip)*sxip
          dbfi(ip) = dbfi(ip) +cn(i)*b1(ip)*sxip
!
        END DO
      END DO

!CDIR NODEP
      DO ip=1, ntest
!
        g1(ip) = y(ip,4)*dcui1(ip) +1.d0
        g2(ip) = y(ip,4)*dcug1(ip) -eot1(ip)
!
        delt(ip) =  tchg*tchg*y(ip,4)*y(ip,4)*b(ip)/tmass +adinv(ip)
        gamm(ip) =  tchg*(cug1(ip)*g1(ip) -cui1(ip)*g2(ip))

!--------------------------------------------c
!       radial electric field
!--------------------------------------------c

!       esp(ip) =       esp0*(1.d0 -y(ip,1)/psia)**2
        def(ip) = -2.d0*esp0*(1.d0 -y(ip,1)/psia)/psia
!
!--------------------------------------------c
!
        c1(ip)  = delt(ip)*db(ip) +tchg*def(ip)
        c21(ip) = c1(ip)*cug1(ip)/gamm(ip)
        c22(ip) = c1(ip)*cui1(ip)/gamm(ip)
        c3(ip)  = tchg*tchg*b(ip)*b(ip)/tmass*y(ip,4)/gamm(ip)
!
        yp(ip,1) = -delt(ip)*( dbth(ip)*cug1(ip) -dbfi(ip)*cui1(ip))/gamm(ip)
        yp(ip,2) =  c21(ip) -g2(ip)*c3(ip)
        yp(ip,3) = -c22(ip) +g1(ip)*c3(ip)
        yp(ip,4) =  delt(ip)*( g2(ip)*dbth(ip)  -g1(ip)*dbfi(ip))/gamm(ip)
!
!-------------------------------------c
!     isw(ip) = -1     yp = 0
!     isw(ip) =  0     yp = 0
!-------------------------------------c
!
        yp(ip,1) = ipsw(ip)*yp(ip,1)
        yp(ip,2) = ipsw(ip)*yp(ip,2)
        yp(ip,3) = ipsw(ip)*yp(ip,3)
        yp(ip,4) = ipsw(ip)*yp(ip,4)
!
!
      END DO
      RETURN
      END
