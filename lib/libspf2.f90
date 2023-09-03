! libspf2.f90

MODULE libspf2

  PRIVATE

  PUBLIC airya  ! Airy functions and their derivatives.
  PUBLIC airyb  ! Airy functions and their derivatives.
  PUBLIC airyzo ! the first NT zeros of Ai(x) and Ai'(x).
  PUBLIC ajyik  ! Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
  PUBLIC aswfa  ! prolate and oblate spheroidal angular functions
                !   of the first kind.
  PUBLIC aswfb  ! prolate and oblate spheroidal angular functions
                !   of the first kind.
  PUBLIC bernoa ! the Bernoulli number Bn.
  PUBLIC bernob ! the Bernoulli number Bn.
  PUBLIC beta   ! the Beta function B(p,q).
  PUBLIC bjndd  ! Bessel functions Jn(x) and first and second derivatives.
  PUBLIC cbk    ! coefficients for oblate radial functions with small argument.
  PUBLIC cchg   ! the confluent hypergeometric function.
  PUBLIC cerf   ! the error function and derivative for a complex argument.
  PUBLIC cerror ! the error function for a complex argument.
  PUBLIC cerzo  ! the complex zeros of the error function.
  PUBLIC cfc    ! the complex Fresnel integral C(z) and C'(z).
  PUBLIC cfs    ! the complex Fresnel integral S(z) and S'(z).
  PUBLIC cgama  ! the Gamma function for complex argument.
  PUBLIC ch12n  ! Hankel functions of first and second kinds, complex argument.
  PUBLIC chgm   ! the confluent hypergeometric function M(a,b,x).
  PUBLIC chgu   ! the confluent hypergeometric function U(a,b,x).
  PUBLIC chgubi ! confluent hypergeometric function with integer argument B.
  PUBLIC chguit ! the hypergeometric function using Gauss-Legendre integration.
  PUBLIC chgul  ! confluent hypergeometric function U(a,b,x)
                !   for large argument X.
  PUBLIC chgus  ! confluent hypergeometric function U(a,b,x)
                !   for small argument X.
  PUBLIC cik01  ! modified Bessel I0(z), I1(z), K0(z) and K1(z)
                !   for complex argument.
  PUBLIC ciklv  ! modified Bessel functions Iv(z), Kv(z),
                !   for complex argument, large order.
  PUBLIC cikna  ! modified Bessel functions In(z), Kn(z), derivatives,
                !   for complex argument.
  PUBLIC ciknb  ! complex modified Bessel functions In(z) and Kn(z).
  PUBLIC cikva  ! modified Bessel functions Iv(z), Kv(z), arbitrary order
                !   for complex argument.
  PUBLIC cikvb  ! modified Bessel functions,Iv(z), Kv(z), arbitrary order
                !   for complex argument.
  PUBLIC cisia  ! cosine Ci(x) and sine integrals Si(x).
  PUBLIC cisib  ! cosine and sine integrals.
  PUBLIC cjk    ! asymptotic expansion coefficients for Bessel functions
                !   of large order.
  PUBLIC cjy01  ! complexBessel functions, derivatives,
                !   J0(z), J1(z), Y0(z), Y1(z).
  PUBLIC cjylv  ! Bessel functions Jv(z), Yv(z) of complex argument
                !    and large order v.
  PUBLIC cjyna  ! Bessel functions and derivatives, Jn(z) and Yn(z)
                !   of complex argument.
  PUBLIC cjynb  ! Bessel functions, derivatives, Jn(z) and Yn(z)
                !   of complex argument.
  PUBLIC cjyva  ! Bessel functions and derivatives, Jv(z) and Yv(z)
                !   of complex argument.
  PUBLIC cjyvb  ! Bessel functions and derivatives, Jv(z) and Yv(z)
                !   of complex argument.
  PUBLIC clpmn  ! associated Legendre functions and derivatives
                !   for complex argument.
  PUBLIC clpn   ! Legendre functions and derivatives for complex argument.
  PUBLIC clqmn  ! associated Legendre functions and derivatives
                !   for complex argument.
  PUBLIC clqn   ! Legendre function Qn(z) and derivative Wn'(z)
                !   for complex argument.
  PUBLIC comelp ! complete elliptic integrals K(k) and E(k).
  PUBLIC cpbdn  ! parabolic cylinder function Dn(z) and Dn'(z)
                !   for complex argument.
  PUBLIC cpdla  ! complex parabolic cylinder function Dn(z) for large argument.
  PUBLIC cpdsa  ! complex parabolic cylinder function Dn(z) for small argument.
  PUBLIC cpsi   ! the psi function for a complex argument.
  PUBLIC csphik ! complex modified spherical Bessel functions and derivatives.
  PUBLIC csphjy ! spherical Bessel functions jn(z) and yn(z)
                !   for complex argument.
  PUBLIC cv0    ! the initial characteristic value of Mathieu functions.
  PUBLIC cva1   ! a sequence of characteristic values of Mathieu functions.
  PUBLIC cva2   ! a specific characteristic value of Mathieu functions.
  PUBLIC cvf    ! F for the characteristic equation of Mathieu functions.
  PUBLIC cvql   ! the characteristic value of Mathieu functions for q <= 3*m.
  PUBLIC cvqm   ! the characteristic value of Mathieu functions for q <= m*m.
  PUBLIC cy01   ! complex Bessel functions Y0(z) and Y1(z) and derivatives.
  PUBLIC cyzo   ! zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
  PUBLIC dvla   ! parabolic cylinder functions Dv(x) for large argument.
  PUBLIC dvsa   ! parabolic cylinder functions Dv(x) for small argument.
  PUBLIC e1xa   ! the exponential integral E1(x).
  PUBLIC e1xb   ! the exponential integral E1(x).
  PUBLIC e1z    ! the complex exponential integral E1(z).
  PUBLIC eix    ! the exponential integral Ei(x).
  PUBLIC elit   ! complete and incomplete elliptic integrals
                !    F(k,phi) and E(k,phi).
  PUBLIC elit3  ! the elliptic integral of the third kind.
!  PUBLIC envj   ! a utility function used by MSTA1 and MSTA2.
  PUBLIC enxa   ! the exponential integral En(x).
  PUBLIC enxb   ! the exponential integral En(x).
  PUBLIC error  ! the error function.
  PUBLIC eulera ! the Euler number En.
  PUBLIC eulerb ! the Euler number En.
  PUBLIC fcoef  ! expansion coefficients for Mathieu and modified
                !    Mathieu functions.
  PUBLIC fcs    ! Fresnel integrals C(x) and S(x).
  PUBLIC fcszo  ! complex zeros of Fresnel integrals C(x) or S(x).
  PUBLIC ffk    ! modified Fresnel integrals F+/-(x) and K+/-(x).
  PUBLIC gaih   ! the GammaH function.
  PUBLIC gam0   ! the Gamma function for the LAMV function.
  PUBLIC gamma  ! the Gamma function.
  PUBLIC gmn    ! quantities for oblate radial functions with small argument.
  PUBLIC herzo  ! the zeros the Hermite polynomial Hn(x).
  PUBLIC hygfx  ! the hypergeometric function F(A,B,C,X).
  PUBLIC hygfz  ! the hypergeometric function F(a,b,c,x) for complex argument.
  PUBLIC ik01a  ! modified Bessel function I0(x), I1(x), K0(x), and K1(x).
  PUBLIC ik01b  ! modified Bessel function I0(x), I1(x), K0(x), and K1(x)
                !    and derivatives.
  PUBLIC ikna   ! modified Bessel function In(x) and Kn(x), and derivatives.
  PUBLIC iknb   ! modified Bessel function In(x) and Kn(x).
  PUBLIC ikv    ! modified Bessel function Iv(x) and Kv(x), and derivatives.
  PUBLIC incob  ! the incomplete beta function Ix(a,b).
  PUBLIC incog  ! the incomplete gamma function r(a,x), ,(a,x), P(a,x).
  PUBLIC itairy ! the integrals of Airy functions.
  PUBLIC itika  ! the integral of the modified Bessel functions
                !    I0(t) and K0(t).
  PUBLIC itikb  ! the integral of the modified Bessel functions
                !    I0(t) and K0(t).
  PUBLIC itjya  ! integrals of Bessel functions J0(t) and Y0(t).
  PUBLIC itjyb  ! integrals of Bessel functions J0(t) and Y0(t).
  PUBLIC itsh0  ! the Struve function H0(t) from 0 to x.
  PUBLIC itsl0  ! the Struve function L0(t) from 0 to x.
  PUBLIC itth0  ! integral of H0(t)/t from x to oo.
  PUBLIC ittika ! integral of (I0(t)-1)/t from 0 to x, K0(t)/t from x to oo.
  PUBLIC ittikb ! integral of (I0(t)-1)/t from 0 to x, K0(t)/t from x to oo.
  PUBLIC ittjya ! integral of (1-J0(t))/t from 0 to x, Y0(t)/t from x to oo.
  PUBLIC ittjyb ! integral of (1-J0(t))/t from 0 to x, Y0(t)/t from x to oo.
  PUBLIC jdzo   ! the zeros of Bessel functions Jn(x) and Jn'(x).
  PUBLIC jelp   ! Jacobian elliptic functions SN(u), CN(u), DN(u).
  PUBLIC jy01a  ! Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  PUBLIC jy01b  ! Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  PUBLIC jyna   ! Bessel functions Jn(x) and Yn(x) and derivatives.
  PUBLIC jynb   ! Bessel functions Jn(x) and Yn(x) and derivatives.
  PUBLIC jyndd  ! Bessel functions Jn(x) and Yn(x),
                !    first and second derivatives.
  PUBLIC jyv    ! Bessel functions Jv(x) and Yv(x) and their derivatives.
  PUBLIC jyzo   ! the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
  PUBLIC klvna  ! Kelvin functions ber(x), bei(x), ker(x), and kei(x),
                !    and derivatives.
  PUBLIC klvnb  ! Kelvin functions ber(x), bei(x), ker(x), and kei(x),
                !    and derivatives.
  PUBLIC klvnzo ! zeros of the Kelvin functions.
  PUBLIC kmn    ! expansion coefficients of prolate or oblate spheroidal
                !    functions.
  PUBLIC lagzo  ! zeros of the Laguerre polynomial, and integration weights.
  PUBLIC lamn   ! lambda functions and derivatives.
  PUBLIC lamv   ! lambda functions and derivatives of arbitrary order.
  PUBLIC legzo  ! the zeros of Legendre polynomials, and integration weights.
  PUBLIC lgama  ! the gamma function or its logarithm.
  PUBLIC lpmn   ! associated Legendre functions Pmn(X) and derivatives P'mn(x).
  PUBLIC lpmns  ! associated Legendre functions Pmn(X) and derivatives P'mn(x).
  PUBLIC lpmv   ! associated Legendre functions Pmv(X) with arbitrary degree.
  PUBLIC lpn    ! Legendre polynomials Pn(x) and derivatives Pn'(x).
  PUBLIC lpni   ! Legendre polynomials Pn(x), derivatives, and integrals.
  PUBLIC lqmn   ! associated Legendre functions Qmn(x) and derivatives.
  PUBLIC lqmns  ! associated Legendre functions Qmn(x) and derivatives Qmn'(x).
  PUBLIC lqna   ! Legendre function Qn(x) and derivatives Qn'(x).
  PUBLIC lqnb   ! Legendre function Qn(x) and derivatives Qn'(x).
!  PUBLIC msta1  ! a backward recurrence starting point for Jn(x).
!  PUBLIC msta2  ! a backward recurrence starting point for Jn(x).
  PUBLIC mtu0   ! Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
  PUBLIC mtu12  ! modified Mathieu functions of the first and second kind.
  PUBLIC othpl  ! orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
  PUBLIC pbdv   ! parabolic cylinder functions Dv(x) and derivatives.
  PUBLIC pbvv   ! parabolic cylinder functions Vv(x) and their derivatives.
  PUBLIC pbwa   ! parabolic cylinder functions W(a,x) and derivatives.
  PUBLIC psi    ! the PSI function.
  PUBLIC qstar  ! Q*mn(-ic) for oblate radial functions with a small argument.
  PUBLIC r8_gamma_log ! the logarithm of the gamma function.
  PUBLIC rctj   ! Riccati-Bessel function of the first kind, and derivatives.
  PUBLIC rcty   ! Riccati-Bessel function of the second kind, and derivatives.
  PUBLIC refine ! an estimate of the characteristic value of Mathieu functions.
  PUBLIC rmn1   ! prolate and oblate spheroidal functions of the first kind.
  PUBLIC rmn2l  ! prolate and oblate spheroidal functions, second kind,
                !    with large CX.
  PUBLIC rmn2so ! oblate radial functions of the second kind
                !    with small argument.
  PUBLIC rmn2sp ! prolate, oblate spheroidal radial functions
                !    of the second, kind 2, with small argument.
  PUBLIC rswfo  ! prolate spheroidal radial function of first and second kinds.
  PUBLIC rswfp  ! prolate spheroidal radial function of first and second kinds.
  PUBLIC scka   ! expansion coefficients for prolate and oblate
                !    spheroidal functions.
  PUBLIC sckb   ! expansion coefficients for prolate and oblate
                !    spheroidal functions.
  PUBLIC sdmn   ! expansion coefficients for prolate and oblate
                !    spheroidal functions.
  PUBLIC segv   ! the characteristic values of spheroidal wave functions.
  PUBLIC sphi   ! mpdified spherical Bessel functions in(x) and derivatives.
  PUBLIC sphj   ! spherical Bessel functions jn(x) and their derivatives.
  PUBLIC sphk   ! modified spherical Bessel functions kn(x) and derivatives.
  PUBLIC sphy   ! spherical Bessel functions yn(x) and their derivatives.
  PUBLIC stvh0  ! the Struve function H0(x).
  PUBLIC stvh1  ! the Struve function H1(x).
  PUBLIC stvhv  ! the Struve function Hv(x) with arbitrary order v.
  PUBLIC stvl0  ! the modified Struve function L0(x).
  PUBLIC stvl1  ! the modified Struve function L1(x).
  PUBLIC stvlv  ! the modified Struve function Lv(x) with arbitary order v.
  PUBLIC timestamp ! prints the current YMDHMS date as a time stamp.
  PUBLIC vvla   ! parabolic cylinder function Vv(x) for large arguments.
  PUBLIC vvsa   ! parabolic cylinder function V(nu,x) for small arguments.

CONTAINS
  
  INCLUDE 'special_functions.f90'

END MODULE libspf2
