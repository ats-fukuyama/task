c*uptodate c05pcftext
      subroutine c05pcf(fcn, n, x, fvec, fjac, ldfjac, xtol,maxfev,
     * diag, mode, factor, nprint, nfev, njev, r, lr, qtf,w, ifail)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     subroutine c05pcf
      implicit real*8(a-h,o-z)
c
c     the purpose of c05pcf is to interface to c05pcz.
c     the latter is based upon minpack routine hybrj.
c     .. scalar arguments ..
      double precision factor, xtol
      integer ifail, ldfjac, lr, maxfev, mode, n, nfev, njev, nprint
c     .. array arguments ..
      double precision diag(n), fjac(ldfjac,n), fvec(n), qtf(n), r(lr),
     * w(n,4),x(n)
c     .. subroutine arguments ..
c     fcn
c     ..
c     .. local scalars ..
c$p 1
      double precision srname
      integer if
c     .. function references ..
      integer p01aaf
c     .. subroutine references ..
c     c05pcz
c     ..
      external fcn
      data srname /8h c05pcf /
      if = 1
      if (n.le.0 .or. ldfjac.lt.n .or. xtol.lt.0.0d0 .or. maxfev.
     *le.0.or. factor.le.0.0d0 .or. lr.lt.(n*(n+1))/2) go to 20
      call c05pcz(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev,diag,
     * mode, factor, nprint, if, nfev, njev, r, lr, qtf,w(1,1), w(1,2),
     * w(1,3), w(1,4))
   20 ifail = p01aaf(ifail,if,srname)
      return
      end
c*uptodate c05pcztext
      subroutine c05pcz(fcn, n, x, fvec, fjac, ldfjac, xtol,maxfev,
     * diag, mode, factor, nprint, info, nfev, njev, r, lr,qtf, wa1,
     * wa2, wa3, wa4)
      implicit real*8(a-h,o-z)
c     mark 9 release.  nag copyright 1981
c     mark 11 revised. ier-439 (feb 1984).
c     **********
c
c     subroutine c05pcz(based on minpack routine hybrj )
c
c     the purpose of c05pcz is to find a zero of a system of
c     n nonlinear functions in n variables by a modification
c     of the powell hybrid method. the user must provide a
c     subroutine which calculates the functions and the jacobian.
c
c     **********
c     .. scalar arguments ..
      double precision factor, xtol
      integer info, ldfjac, lr, maxfev, mode, n, nfev, njev, nprint
c     .. array arguments ..
      double precision diag(n), fjac(ldfjac,n), fvec(n), qtf(n), r(lr),
     * wa1(n),wa2(n), wa3(n), wa4(n), x(n)
c     .. subroutine arguments ..
c     fcn
c     ..
c     .. local scalars ..
      double precision actred, delta, epsmch, fnorm, fnorm1, one, pnorm,
     *prered, p0001, p001, p1, p5, ratio, srname, sum, temp,xnorm, zero
      integer i, iflag, iter, j, jm1, l, ncfail, ncsuc, nslow1,nslow2
      logical jeval, sing
c     .. local arrays ..
      integer iwa(1)
c     .. function references ..
      double precision c05nct, x02aaf
c     .. subroutine references ..
c     c05ncu, c05ncw, c05ncx, c05ncy, c05ncz
c     ..
      external fcn
      data one, p1, p5, p001, p0001, zero /1.0d0,1.0d-1,5.0d-1,1.0d-3,
     *1.0d-4,0.0d0/
      epsmch = x02aaf(0.0d0)
c
      info = 1
      iflag = 0
      nfev = 0
      njev = 0
c
c     check the input parameters for errors.
c
      if (n.le.0 .or. ldfjac.lt.n .or. xtol.lt.zero .or.maxfev.le.0 .or.
     * factor.le.zero .or. lr.lt.(n*(n+1))/2) goto 600
      if (mode.ne.2) go to 40
      do 20 j=1,n
         if (diag(j).le.zero) go to 600
   20 continue
   40 continue
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
      call fcn(n, x, fvec, fjac, ldfjac, iflag)
      nfev = 1
      if (iflag.lt.0) go to 600
      fnorm = c05nct(n,fvec)
c
c     initialize iteration counter and monitors.
c
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
c
c     beginning of the outer loop.
c
   60 continue
      jeval = .true.
c
c     calculate the jacobian matrix.
c
      iflag = 2
      call fcn(n, x, fvec, fjac, ldfjac, iflag)
      njev = njev + 1
      if (iflag.lt.0) go to 600
c
c     compute the qr factorization of the jacobian.
c
      call c05ncx(n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2,wa3)
c
c     on the first iteration and if mode is 1, scale according
c     to the norms of the columns of the initial jacobian.
c
      if (iter.ne.1) go to 140
      if (mode.eq.2) go to 100
      do 80 j=1,n
         diag(j) = wa2(j)
         if (wa2(j).eq.zero) diag(j) = one
   80 continue
  100 continue
c
c     on the first iteration, calculate the norm of the scaled x
c     and initialize the step bound delta.
c
      do 120 j=1,n
         wa3(j) = diag(j)*x(j)
  120 continue
      xnorm = c05nct(n,wa3)
      delta = factor*xnorm
      if (delta.eq.zero) delta = factor
  140 continue
c
c     form (q transpose)*fvec and store in qtf.
c
      do 160 i=1,n
         qtf(i) = fvec(i)
  160 continue
      do 240 j=1,n
         if (fjac(j,j).eq.zero) go to 220
         sum = zero
         do 180 i=j,n
            sum = sum + fjac(i,j)*qtf(i)
  180    continue
         temp = -sum/fjac(j,j)
         do 200 i=j,n
            qtf(i) = qtf(i) + fjac(i,j)*temp
  200    continue
  220    continue
  240 continue
c
c     copy the triangular factor of the qr factorization into r.
c
      sing = .false.
      do 300 j=1,n
         l = j
         jm1 = j - 1
         if (jm1.lt.1) go to 280
         do 260 i=1,jm1
            r(l) = fjac(i,j)
            l = l + n - i
  260    continue
  280    continue
         r(l) = wa1(j)
         if (wa1(j).eq.zero) sing = .true.
  300 continue
c
c     accumulate the orthogonal factor in fjac.
c
      call c05ncw(n, n, fjac, ldfjac, wa1)
c
c     rescale if necessary.
c
      if (mode.eq.2) go to 340
      do 320 j=1,n
         diag(j) = dmax1(diag(j),wa2(j))
  320 continue
  340 continue
c
c     beginning of the inner loop.
c
  360 continue
c
c     if requested, call fcn to enable printing of iterates.
c
      if (nprint.le.0) go to 380
      iflag = 0
      if (mod(iter-1,nprint).eq.0) call fcn(n, x, fvec, fjac,ldfjac,
     * iflag)
      if (iflag.lt.0) go to 600
  380 continue
c
c     determine the direction p.
c
      call c05ncu(n, r, lr, diag, qtf, delta, wa1, wa2, wa3)
c
c     store the direction p and x + p. calculate the norm of p.
c
      do 400 j=1,n
         wa1(j) = -wa1(j)
         wa2(j) = x(j) + wa1(j)
         wa3(j) = diag(j)*wa1(j)
  400 continue
      pnorm = c05nct(n,wa3)
c
c     on the first iteration, adjust the initial step bound.
c
      if (iter.eq.1) delta = dmin1(delta,pnorm)
c
c     evaluate the function at x + p and calculate its norm.
c
      iflag = 1
      call fcn(n, wa2, wa4, fjac, ldfjac, iflag)
      nfev = nfev + 1
      if (iflag.lt.0) go to 600
      fnorm1 = c05nct(n,wa4)
c
c     compute the scaled actual reduction.
c
      actred = -one
      if (fnorm1.lt.fnorm) actred = one - (fnorm1/fnorm)**2
c
c     compute the scaled predicted reduction.
c
      l = 1
      do 440 i=1,n
         sum = zero
         do 420 j=i,n
            sum = sum + r(l)*wa1(j)
            l = l + 1
  420    continue
         wa3(i) = qtf(i) + sum
  440 continue
      temp = c05nct(n,wa3)
      prered = 1.0d0
      if (temp.lt.fnorm) prered = one - (temp/fnorm)**2
c
c     compute the ratio of the actual to the predicted
c     reduction.
c
      ratio = zero
      if (prered.gt.zero) ratio = actred/prered
c
c     update the step bound.
c
      if (ratio.ge.p1) go to 460
      ncsuc = 0
      ncfail = ncfail + 1
      delta = p5*delta
      go to 480
  460 continue
      ncfail = 0
      ncsuc = ncsuc + 1
      if (ratio.ge.p5 .or. ncsuc.gt.1) delta = dmax1(delta,pnorm/p5)
      if (dabs(ratio-one).le.p1) delta = pnorm/p5
  480 continue
c
c     test for successful iteration.
c
      if (ratio.lt.p0001) go to 520
c
c     successful iteration. update x, fvec, and their norms.
c
      do 500 j=1,n
         x(j) = wa2(j)
         wa2(j) = diag(j)*x(j)
         fvec(j) = wa4(j)
  500 continue
      xnorm = c05nct(n,wa2)
      fnorm = fnorm1
      iter = iter + 1
  520 continue
c
c     determine the progress of the iteration.
c
      nslow1 = nslow1 + 1
      if (actred.ge.p001) nslow1 = 0
      if (jeval) nslow2 = nslow2 + 1
      if (actred.ge.p1) nslow2 = 0
c
c     test for convergence.
c
      if (delta.le.xtol*xnorm .or. fnorm.eq.zero) info = 0
      if (info.ne.1) go to 600
c
c     tests for termination and stringent tolerances.
c
      if (nfev.ge.maxfev) info = 2
      if (delta.le.epsmch*xnorm) info = 3
      if (nslow2.eq.5) info = 4
      if (nslow1.eq.10) info = 5
      if (info.ne.1) go to 600
c
c     criterion for recalculating jacobian.
c
      if (ncfail.eq.2) go to 580
c
c     calculate the rank one modification to the jacobian
c     and update qtf if necessary.
c
      do 560 j=1,n
         sum = zero
         do 540 i=1,n
            sum = sum + fjac(i,j)*wa4(i)
  540    continue
         wa2(j) = (sum-wa3(j))/pnorm
         wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
         if (ratio.ge.p0001) qtf(j) = sum
  560 continue
c
c     compute the qr factorization of the updated jacobian.
c
      call c05ncz(n, n, r, lr, wa1, wa2, wa3, sing)
      call c05ncy(n, n, fjac, ldfjac, wa2, wa3)
      call c05ncy(1, n, qtf, 1, wa2, wa3)
c
c     end of the inner loop.
c
      jeval = .false.
      go to 360
  580 continue
c
c     end of the outer loop.
c
      go to 60
  600 continue
c
c     termination, either normal or user imposed.
c
      if (iflag.lt.0) info = iflag
      iflag = 0
      if (nprint.gt.0) call fcn(n, x, fvec, fjac, ldfjac, iflag)
      return
      end
c *uptodate p01aaftext
      integer function p01aaf(ifail, error, srname)
      common /cipacc/ ipacc
c     mark 1 release.  nag copyright 1971
c     mark 3 revised
c     mark 4a revised, ier-45
c     mark 4.5 revised
c     mark 7 revised (dec 1978)
c     mark 11 revised (feb 1984)
c     returns the value of error or terminates the program.
      integer error, ifail, nout
c$p 1
      double precision srname
c     test if no error detected
      if (error.eq.0) go to 20
c     determine output unit for message
      call x04aaf (0,nout)
c     test for soft failure
      if (mod(ifail,10).eq.1) go to 10
c     hard failure
cx    write (nout,99999) srname, error
cx    write (7,99999) srname, error ,ipacc
c     ******************** implementation note ********************
c     the following stop statement may be replaced by a call to an
c     implementation-dependent routine to display a message and/or
c     abort the program.
c     *************************************************************
      go to 20
cx    return
c     stop
c     soft fail
c     test if error messages suppressed
   10 if (mod(ifail/10,10).eq.0) go to 20
cx    write (nout,99999) srname, error
cx    write (7,99999) srname, error,ipacc
   20 p01aaf = error
      return
99999 format (1h , 38herror detected by nag library routine , a8,
     & 11h - ifail = , i5,'  ipacc =',i5)
      end
cc uptodate c05ncttext
      double precision function c05nct(n, x)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     function c05nct
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c     real function c05nct(n,x)
c
c     where
c
c     n is a positive integer input variable.
c
c     x is an input array of length n.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c     .. scalar arguments ..
      integer n
c     .. array arguments ..
      double precision x(n)
c     ..
c     .. local scalars ..
      double precision agiant, floatn, one, rdwarf, rgiant, s1, s2, s3,
     * xabs,x1max, x3max, zero
      integer i
c     .. function references ..
      double precision dsqrt
c     ..
      data one, zero, rdwarf, rgiant /1.0d0,0.0d0,1.0d-19,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 180 i=1,n
         xabs = dabs(x(i))
         if (xabs.gt.rdwarf .and. xabs.lt.agiant) go to 140
         if (xabs.le.rdwarf) go to 60
c
c     sum for large components.
c
         if (xabs.le.x1max) go to 20
         s1 = one + s1*(x1max/xabs)**2
         x1max = xabs
         go to 40
   20    continue
         s1 = s1 + (xabs/x1max)**2
   40    continue
         go to 120
   60    continue
c
c     sum for small components.
c
         if (xabs.le.x3max) go to 80
         s3 = one + s3*(x3max/xabs)**2
         x3max = xabs
         go to 100
   80    continue
         if (xabs.ne.zero) s3 = s3 + (xabs/x3max)**2
  100    continue
  120    continue
         go to 160
  140    continue
c
c     sum for intermediate components.
c
         s2 = s2 + xabs**2
  160    continue
  180 continue
c
c     calculation of norm.
c
      if (s1.eq.zero) go to 200
      c05nct = x1max*dsqrt(s1+(s2/x1max)/x1max)
      go to 260
  200 continue
      if (s2.eq.zero) go to 220
      if (s2.ge.x3max) c05nct = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
      if (s2.lt.x3max) c05nct = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
      go to 240
  220 continue
      c05nct = x3max*dsqrt(s3)
  240 continue
  260 continue
      return
      end
c  *uptodate c05ncutext
      subroutine c05ncu(n, r, lr, diag, qtb, delta, x, wa1, wa2)
      implicit real*8(a-h,o-z)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     subroutine c05ncu(based on minpack routine dogleg)
c
c     given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta, the
c     problem is to determine the convex combination x of the
c     gauss-newton and scaled gradient directions that minimizes
c     (a*x - b) in the least squares sense, subject to the
c     restriction that the euclidean norm of d*x be at most delta.
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization of a. that is, if a = q*r, where q has
c     orthogonal columns and r is an upper triangular matrix,
c     then c05ncu expects the full upper triangle of r and
c     the first n components of (q transpose)*b.
c
c     the subroutine statement is
c
c     subroutine c05ncu(n,r,lr,diag,qtb,delta,x,wa1,wa2)
c
c     where
c
c     n is a positive integer input variable set to the order of r.
c
c     r is an input array of length lr which must contain the upper
c     triangular matrix r stored by rows.
c
c     lr is a positive integer input variable not less than
c     (n*(n+1))/2.
c
c     diag is an input array of length n which must contain the
c     diagonal elements of the matrix d.
c
c     qtb is an input array of length n which must contain the first
c     n elements of the vector (q transpose)*b.
c
c     delta is a positive input variable which specifies an upper
c     bound on the euclidean norm of d*x.
c
c     x is an output array of length n which contains the desired
c     convex combination of the gauss-newton direction and the
c     scaled gradient direction.
c
c     wa1 and wa2 are work arrays of length n.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c     .. scalar arguments ..
      double precision delta
      integer lr, n
c     .. array arguments ..
      double precision diag(n), qtb(n), r(lr), wa1(n), wa2(n), x(n)
c     ..
c     .. local scalars ..
      double precision alpha, bnorm, epsmch, gnorm, one, qnorm, sgnorm,
     * sum,temp, zero
      integer i, j, jj, jp1, k, l
c     .. function references ..
      double precision c05nct, dsqrt, x02aaf
c     ..
      data one, zero /1.0d0,0.0d0/
      epsmch = x02aaf(0.0d0)
c
c     first, calculate the gauss-newton direction.
c
      jj = (n*(n+1))/2 + 1
      do 100 k=1,n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if (n.lt.jp1) go to 40
         do 20 i=jp1,n
            sum = sum + r(l)*x(i)
            l = l + 1
   20    continue
   40    continue
         temp = r(jj)
         if (temp.ne.zero) go to 80
         l = j
         do 60 i=1,j
            temp = dmax1(temp,dabs(r(l)))
            l = l + n - i
   60    continue
         temp = epsmch*temp
         if (temp.eq.zero) temp = epsmch
   80    continue
         x(j) = (qtb(j)-sum)/temp
  100 continue
c
c     test whether the gauss-newton direction is acceptable.
c
      do 120 j=1,n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
  120 continue
      qnorm = c05nct(n,wa2)
      if (qnorm.le.delta) go to 280
c
c     the gauss-newton direction is not acceptable.
c     next, calculate the scaled gradient direction.
c
      l = 1
      do 160 j=1,n
         temp = qtb(j)
         do 140 i=j,n
            wa1(i) = wa1(i) + r(l)*temp
            l = l + 1
  140    continue
         wa1(j) = wa1(j)/diag(j)
  160 continue
c
c     calculate the norm of the scaled gradient and test for
c     the special case in which the scaled gradient is zero.
c
      gnorm = c05nct(n,wa1)
      sgnorm = zero
      alpha = delta/qnorm
      if (gnorm.eq.zero) go to 240
c
c     calculate the point along the scaled gradient
c     at which the quadratic is minimized.
c
      do 180 j=1,n
         wa1(j) = (wa1(j)/gnorm)/diag(j)
  180 continue
      l = 1
      do 220 j=1,n
         sum = zero
         do 200 i=j,n
            sum = sum + r(l)*wa1(i)
            l = l + 1
  200    continue
         wa2(j) = sum
  220 continue
      temp = c05nct(n,wa2)
      sgnorm = (gnorm/temp)/temp
c
c     test whether the scaled gradient direction is acceptable.
c
      alpha = zero
      if (sgnorm.ge.delta) go to 240
c
c     the scaled gradient direction is not acceptable.
c     finally, calculate the point along the c05ncu
c     at which the quadratic is minimized.
c
      bnorm = c05nct(n,qtb)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2 +dsqrt((temp-(delta/
     *qnorm))**2+(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(one-(sgnorm/delta)**2))/temp
  240 continue
c
c     form appropriate convex combination of the gauss-newton
c     direction and the scaled gradient direction.
c
      temp = (one-alpha)*dmin1(sgnorm,delta)
      do 260 j=1,n
         x(j) = temp*wa1(j) + alpha*x(j)
  260 continue
  280 continue
      return
      end
c*uptodate c05ncwtext
      subroutine c05ncw(m, n, q, ldq, wa)
      implicit real*8(a-h,o-z)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     subroutine c05ncw(based on minpack routine qform )
c
c     this subroutine proceeds from the computed qr factorization of
c     an m by n matrix a to accumulate the m by m orthogonal matrix
c     q from its factored form.
c
c     the subroutine statement is
c
c     subroutine c05ncw(m,n,q,ldq,wa)
c
c     where
c
c     m is a positive integer input variable set to the number
c     of rows of a and the order of q.
c
c     n is a positive integer input variable set to the number
c     of columns of a.
c
c     q is an m by m array. on input the full lower trapezoid in
c     the first min(m,n) columns of q contains the factored form.
c     on output q has been accumulated into a square matrix.
c
c     ldq is a positive integer input variable not less than m
c     which specifies the leading dimension of the array q.
c
c     wa is a work array of length m.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c     .. scalar arguments ..
      integer ldq, m, n
c     .. array arguments ..
      double precision q(ldq,m), wa(m)
c     ..
c     .. local scalars ..
      double precision one, sum, temp, zero
      integer i, j, jm1, k, l, minmn, np1
c     ..
      data one, zero /1.0d0,0.0d0/
      minmn = min0(m,n)
      if (minmn.lt.2) go to 60
      do 40 j=2,minmn
         jm1 = j - 1
         do 20 i=1,jm1
            q(i,j) = zero
   20    continue
   40 continue
   60 continue
c
c     initialize remaining columns to those of the identity matrix.
c
      np1 = n + 1
      if (m.lt.np1) go to 120
      do 100 j=np1,m
         do 80 i=1,m
            q(i,j) = zero
   80    continue
         q(j,j) = one
  100 continue
  120 continue
c
c     accumulate q from its factored form.
c
      do 240 l=1,minmn
         k = minmn - l + 1
         do 140 i=k,m
            wa(i) = q(i,k)
            q(i,k) = zero
  140    continue
         q(k,k) = one
         if (wa(k).eq.zero) go to 220
         do 200 j=k,m
            sum = zero
            do 160 i=k,m
               sum = sum + q(i,j)*wa(i)
  160       continue
            temp = sum/wa(k)
            do 180 i=k,m
               q(i,j) = q(i,j) - temp*wa(i)
  180       continue
  200    continue
  220    continue
  240 continue
      return
      end
c*uptodate c05ncxtext
      subroutine c05ncx(m, n, a, lda, pivot, ipvt, lipvt, rdiag,acnorm,
     * wa)
      implicit real*8(a-h,o-z)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     subroutine c05ncx(based on minpack routine qrfac )
c
c     this subroutine uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, c05ncx determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c     t
c     i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack subroutine.
c
c     the subroutine statement is
c
c     subroutine c05ncx(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     where
c
c     m is a positive integer input variable set to the number
c     of rows of a.
c
c     n is a positive integer input variable set to the number
c     of columns of a.
c
c     a is an m by n array. on input a contains the matrix for
c     which the qr factorization is to be computed. on output
c     the strict upper trapezoidal part of a contains the strict
c     upper trapezoidal part of r, and the lower trapezoidal
c     part of a contains a factored form of q (the non-trivial
c     elements of the u vectors described above).
c
c     lda is a positive integer input variable not less than m
c     which specifies the leading dimension of the array a.
c
c     pivot is a logical input variable. if pivot is set true,
c     then column pivoting is enforced. if pivot is set false,
c     then no column pivoting is done.
c
c     ipvt is an integer output array of length lipvt. ipvt
c     defines the permutation matrix p such that a*p = q*r.
c     column j of p is column ipvt(j) of the identity matrix.
c     if pivot is false, ipvt is not referenced.
c
c     lipvt is a positive integer input variable. if pivot is false,
c     then lipvt may be as small as 1. if pivot is true, then
c     lipvt must be at least n.
c
c     rdiag is an output array of length n which contains the
c     diagonal elements of r.
c
c     acnorm is an output array of length n which contains the
c     norms of the corresponding columns of the input matrix a.
c     if this information is not needed, then acnorm can coincide
c     with rdiag.
c
c     wa is a work array of length n. if pivot is false, then wa
c     can coincide with rdiag.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c     .. scalar arguments ..
      integer lda, lipvt, m, n
      logical pivot
c     .. array arguments ..
      double precision a(lda,n), acnorm(n), rdiag(n), wa(n)
      integer ipvt(lipvt)
c     ..
c     .. local scalars ..
      double precision ajnorm, epsmch, one, p05, sum, temp, zero
      integer i, j, jp1, k, kmax, minmn
c     .. function references ..
      double precision c05nct, dsqrt, x02aaf
c     ..
      data one, p05, zero /1.0d0,5.0d-2,0.0d0/
      epsmch = x02aaf(0.0d0)
c
c     compute the initial column norms and initialize several
c     arrays.
c
      do 20 j=1,n
         acnorm(j) = c05nct(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   20 continue
c
c     reduce a to r with householder transformations.
c
      minmn = min0(m,n)
      do 220 j=1,minmn
         if (.not.pivot) go to 80
c
c     bring the column of largest norm into the pivot position.
c
         kmax = j
         do 40 k=j,n
            if (rdiag(k).gt.rdiag(kmax)) kmax = k
   40    continue
         if (kmax.eq.j) go to 80
         do 60 i=1,m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   60    continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   80    continue
c
c     compute the householder transformation to reduce the
c     j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = c05nct(m-j+1,a(j,j))
         if (ajnorm.eq.zero) go to 200
         if (a(j,j).lt.zero) ajnorm = -ajnorm
         do 100 i=j,m
            a(i,j) = a(i,j)/ajnorm
  100    continue
         a(j,j) = a(j,j) + one
c
c     apply the transformation to the remaining columns
c     and update the norms.
c
         jp1 = j + 1
         if (n.lt.jp1) go to 200
         do 180 k=jp1,n
            sum = zero
            do 120 i=j,m
               sum = sum + a(i,j)*a(i,k)
  120       continue
            temp = sum/a(j,j)
            do 140 i=j,m
               a(i,k) = a(i,k) - temp*a(i,j)
  140       continue
            if (.not.pivot .or. rdiag(k).eq.zero) go to 160
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2.gt.epsmch) go to 160
            rdiag(k) = c05nct(m-j,a(jp1,k))
            wa(k) = rdiag(k)
  160       continue
  180    continue
  200    continue
         rdiag(j) = -ajnorm
  220 continue
      return
      end
c*uptodate c05ncytext
      subroutine c05ncy(m, n, a, lda, v, w)
      implicit real*8(a-h,o-z)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     subroutine c05ncy(based on minpack routine r1mpyq)
c
c     given an m by n matrix a, this subroutine computes a*q where
c     q is the product of 2*(n - 1) transformations
c
c     gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
c     and gv(i), gw(i) are givens rotations in the (i,n) plane which
c     eliminate elements in the i-th and n-th planes, respectively.
c     q itself is not given, rather the information to recover the
c     gv, gw rotations is supplied.
c
c     the subroutine statement is
c
c     subroutine c05ncy(m,n,a,lda,v,w)
c
c     where
c
c     m is a positive integer input variable set to the number
c     of rows of a.
c
c     n is a positive integer input variable set to the number
c     of columns of a.
c
c     a is an m by n array. on input a must contain the matrix
c     to be postmultiplied by the orthogonal matrix q
c     described above. on output a*q has replaced a.
c
c     lda is a positive integer input variable not less than m
c     which specifies the leading dimension of the array a.
c
c     v is an input array of length n. v(i) must contain the
c     information necessary to recover the givens rotation gv(i)
c     described above.
c
c     w is an input array of length n. w(i) must contain the
c     information necessary to recover the givens rotation gw(i)
c     described above.
c
c     subroutines called
c
c     fortran-supplied ... abs,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
c     .. scalar arguments ..
      integer lda, m, n
c     .. array arguments ..
      double precision a(lda,n), v(n), w(n)
c     ..
c     .. local scalars ..
      double precision aacos, one, aasin, temp
      integer i, j, nmj, nm1
c     .. function references ..
      double precision dsqrt
c     ..
      data one /1.0d0/
c
c     apply the first set of givens rotations to a.
c
      nm1 = n - 1
      if (nm1.lt.1) go to 100
      do 40 nmj=1,nm1
         j = n - nmj
         if (dabs(v(j)).gt.one) aacos = one/v(j)
         if (dabs(v(j)).gt.one) aasin = dsqrt(one-aacos**2)
         if (dabs(v(j)).le.one) aasin = v(j)
         if (dabs(v(j)).le.one) aacos = dsqrt(one-aasin**2)
         do 20 i=1,m
            temp = aacos*a(i,j) - aasin*a(i,n)
            a(i,n) = aasin*a(i,j) + aacos*a(i,n)
            a(i,j) = temp
   20    continue
   40 continue
c
c     apply the second set of givens rotations to a.
c
      do 80 j=1,nm1
         if (dabs(w(j)).gt.one) aacos = one/w(j)
         if (dabs(w(j)).gt.one) aasin = dsqrt(one-aacos**2)
         if (dabs(w(j)).le.one) aasin = w(j)
         if (dabs(w(j)).le.one) aacos = dsqrt(one-aasin**2)
         do 60 i=1,m
            temp = aacos*a(i,j) + aasin*a(i,n)
            a(i,n) = -aasin*a(i,j) + aacos*a(i,n)
            a(i,j) = temp
   60    continue
   80 continue
  100 continue
      return
      end
c*uptodate c05ncztext
      subroutine c05ncz(m, n, s, ls, u, v, w, sing)
      implicit real*8(a-h,o-z)
c     mark 9 release.  nag copyright 1981
c     **********
c
c     subroutine c05ncz(based on minpack routine r1updt)
c
c     given an m by n lower trapezoidal matrix s, an m-vector u,
c     and an n-vector v, the problem is to determine an
c     orthogonal matrix q such that
c
c     t
c     (s + u*v )*q
c
c     is again lower trapezoidal.
c
c     this subroutine determines q as the product of 2*(n - 1)
c     transformations
c
c     gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
c     where gv(i), gw(i) are givens rotations in the (i,n) plane
c     which eliminate elements in the i-th and n-th planes,
c     respectively. q itself is not accumulated, rather the
c     information to recover the gv, gw rotations is returned.
c
c     the subroutine statement is
c
c     subroutine c05ncz(m,n,s,ls,u,v,w,sing)
c
c     where
c
c     m is a positive integer input variable set to the number
c     of rows of s.
c
c     n is a positive integer input variable set to the number
c     of columns of s. n must not exceed m.
c
c     s is an array of length ls. on input s must contain the lower
c     trapezoidal matrix s stored by columns. on output s contains
c     the lower trapezoidal matrix produced as described above.
c
c     ls is a positive integer input variable not less than
c     (n*(2*m-n+1))/2.
c
c     u is an input array of length m which must contain the
c     vector u.
c
c     v is an array of length n. on input v must contain the vector
c     v. on output v(i) contains the information necessary to
c     recover the givens rotation gv(i) described above.
c
c     w is an output array of length m. w(i) contains information
c     necessary to recover the givens rotation gw(i) described
c     above.
c
c     sing is a logical output variable. sing is set true if any
c     of the diagonal elements of the output s are zero. otherwise
c     sing is set false.
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more,
c     john l. nazareth
c
c     **********
c     .. scalar arguments ..
      integer ls, m, n
      logical sing
c     .. array arguments ..
      double precision s(ls), u(m), v(n), w(m)
c     ..
c     .. local scalars ..
      double precision aacos, cotan, giant, one, p25, p5, aasin, tan, tau,
     * temp,zero
      integer i, j, jj, l, nmj, nm1
c     .. function references ..
      double precision dsqrt, x02acf
c     ..
      data one, p5, p25, zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
      giant = x02acf(0.0d0)
c
c     initialize the diagonal element pointer.
c
      jj = (n*(2*m-n+1))/2 - (m-n)
c
c     move the nontrivial part of the last column of s into w.
c
      l = jj
      do 20 i=n,m
         w(i) = s(l)
         l = l + 1
   20 continue
c
c     rotate the vector v into a multiple of the n-th unit vector
c     in such a way that a spike is introduced into w.
c
      nm1 = n - 1
      if (nm1.lt.1) go to 140
      do 120 nmj=1,nm1
         j = n - nmj
         jj = jj - (m-j+1)
         w(j) = zero
         if (v(j).eq.zero) go to 100
c
c     determine a givens rotation which eliminates the
c     j-th element of v.
c
         if (dabs(v(n)).ge.dabs(v(j))) go to 40
         cotan = v(n)/v(j)
         aasin = p5/dsqrt(p25+p25*cotan**2)
         aacos = aasin*cotan
         tau = one
         if (dabs(aacos)*giant.gt.one) tau = one/aacos
         go to 60
   40    continue
         tan = v(j)/v(n)
         aacos = p5/dsqrt(p25+p25*tan**2)
         aasin = aacos*tan
         tau = aasin
   60    continue
c
c     apply the transformation to v and store the information
c     necessary to recover the givens rotation.
c
         v(n) = aasin*v(j) + aacos*v(n)
         v(j) = tau
c
c     apply the transformation to s and extend the spike in w.
c
         l = jj
         do 80 i=j,m
            temp = aacos*s(l) - aasin*w(i)
            w(i) = aasin*s(l) + aacos*w(i)
            s(l) = temp
            l = l + 1
   80    continue
  100    continue
  120 continue
  140 continue
c
c     add the spike from the rank 1 update to w.
c
      do 160 i=1,m
         w(i) = w(i) + v(n)*u(i)
  160 continue
c
c     eliminate the spike.
c
      sing = .false.
      if (nm1.lt.1) go to 280
      do 260 j=1,nm1
         if (w(j).eq.zero) go to 240
c
c     determine a givens rotation which eliminates the
c     j-th element of the spike.
c
         if (dabs(s(jj)).ge.dabs(w(j))) go to 180
         cotan = s(jj)/w(j)
         aasin = p5/dsqrt(p25+p25*cotan**2)
         aacos = aasin*cotan
         tau = one
         if (dabs(aacos)*giant.gt.one) tau = one/aacos
         go to 200
  180    continue
         tan = w(j)/s(jj)
         aacos = p5/dsqrt(p25+p25*tan**2)
         aasin = aacos*tan
         tau = aasin
  200    continue
c
c     apply the transformation to s and reduce the spike in w.
c
         l = jj
         do 220 i=j,m
            temp = aacos*s(l) + aasin*w(i)
            w(i) = -aasin*s(l) + aacos*w(i)
            s(l) = temp
            l = l + 1
  220    continue
c
c     store the information necessary to recover the
c     givens rotation.
c
         w(j) = tau
  240    continue
c
c     test for zero diagonal elements in the output s.
c
         if (s(jj).eq.zero) sing = .true.
         jj = jj + (m-j+1)
  260 continue
  280 continue
c
c     move w back into the last column of the output s.
c
      l = jj
      do 300 i=n,m
         s(l) = w(i)
         l = l + 1
  300 continue
      if (s(jj).eq.zero) sing = .true.
      return
      end
c*uptodate x02aaftext
      double precision function x02aaf(x)
c     nag copyright 1975
c     mark 4.5 release
      double precision x
c     * eps *
c
c     ibm double precision version
c
c     returns the value eps where eps is the smallest
c     positive
c     number such that 1.0 + eps > 1.0
c     the x parameter is not used
c     for icl 1900
c     x02aaf = 2.0**(-37.0)
c     for ibm 360/370
c     x02aaf = 2.0d0**(-52.0d0)
      double precision z
c      data z/z3410000000000000/
      x02aaf = 2.2204460492503131d-16
      return
      end
c*uptodate x02acftext
      double precision function x02acf(x)
c     nag copyright 1975
c     mark 4.5 release
      double precision x
c     * rmax *
c
c     ibm double precision version
c
c     returns the value of the largest positive real  floating-
c     point number representable on the computer
c     for icl 1900
c     x02acf = (2.0 - 2.0**(-36.0))*2.0**254.0
c     for ibm 360/370
c     x02acf = (1.0d0-16.0d0**(-14.0d0))*16.0d0**63.0d0
      double precision z
c      data z/z7fffffffffffffff/
c NEC
c     x02acf = 7.2370055773322621d+75
      x02acf = 1.797693134862315D+308
      return
      end
c*uptodate x04aaftext
      subroutine x04aaf(i,nerr)
      implicit real*8(a-h,o-z)
c     mark 7 release. nag copyright 1978
c     mark 7c revised ier-190 (may 1979)
c     if i = 0, sets nerr to current error message unit number
c     (stored in nerr1).
c     if i = 1, changes current error message unit number to
c     value specified by nerr.
c
c     *** note ***
c     this routine assumes that the value of nerr1 is saved
c     between calls.  in some implementations it may be
c     necessary to store nerr1 in a labelled common
c     block /ax04aa/ to achieve this.
c
c     .. scalar arguments ..
      integer i, nerr
c     ..
c     .. local scalars ..
      integer nerr1
c     ..
      data nerr1 /6/
      if (i.eq.0) nerr = nerr1
      if (i.eq.1) nerr1 = nerr
      return
      end
