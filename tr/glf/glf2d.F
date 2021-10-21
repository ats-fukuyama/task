! 12mar2003 jek added rms-theta for retuned version
! 02feb2001 jek fixed section w/ xparam(24) for power law ExB
! 08may2000 jek merged ITG and ETG loops over ky for better MPI optimization
! 23nov1999 jek added zgeev eigenvalue solver (eigen_gf=2)
! 29mar1999 fgtok -s cgg.table "dmc:  rename eispack routines"
! 29mar1999 fgtok -s rr.table "dmc:  rename intrinsic REAL -> AREAL"
!
! dmc -- Cray/workstation portable real*8<-->complex*16 conversion routines
 
#include "f77_dcomplx.h"
 
c
cglf2d.f 12-mar-03 Kinsey
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
 
           subroutine glf2d(iglf)
 
c
c***********************************************************************
c questions  should be addressed to
c  Ron Waltz 619-455-4584  or email: waltz@gav.gat.com
c***********************************************************************
 
 
c 2D GLF equations with massless isothermal passing electrons from
c  Waltz et al, Phys. of Plasmas 6(1995)2408
 
c In Eq. 24  p_u_par is replaced by n_u/(1-reps) +t_u and
c the isothermal conditions t_u= (betae/2)*w_s/k_par*(rlt)*a_par is
c used. Thus n_u = (1-reps)*ph (adiabatic passing electrons) for betae=0
c
c In Eq. 23 (p_u_par+p_u_per) is replaced by n_u+t_u
c using the isothermal condition the mHD beta limit is too low by
c 1/(1+ reps(-0.25+0.75rlte)/(taui*(rln+rlti)+(rln+rlte))
c It is possible to patch this up by replacing with x*n_u+y*t_u
c then solving for x and y to obtain
c universal MHD betae_crit=2*k_par**2/(w_d*(taui*(rln+rlti)+w_d*(rln+rlti)))
c beta_crit=(1+taui)*betae_crit=(1/2)(1/q**2)*(L_p/rmajor)
c 1/2 is replaced by s_hat in models with shear
 
c EVERYTHING else including normalizing units follows the paper.
 
c  unit of microlength is rho_s, unit of macrolength is a
c   a is typically rho value at separatrix
c  unit of time is a/c_s; unit of diffusion is (c_s/a)*rho_s**2
c  c_s=sqrt(Te/M_i), omega=eB/cM_i, rho_s=c_s/omega
 
c example balance equations to clarify meaning  of diffusivities
c
c       chie_hat is effective energy diffusivity
c
c  (3/2) d ne Te/ dt =
c -1/V(rho)_prime d /d rho V(rho)_prime
c             |grad rho|**2 ne chie_hat (c_s rho_s**2/a)(-d Te/ d rho)
c      -exch_hat (ne Te) c_s/a (rho_s/a)**2 +heating density
c
c and similarly for ion equation
c note that no convective part is added, ie "convection" is included
c inside chie_hat
c note no impurity energy flow is computed although it could be easily done
 
c        d_hat is effective plasma diffusivity for ions
c
c        d ni / dt =
c -1/V(rho)_prime d /d rho V(rho)_prime
c               |grad rho|**2 ne d_hat (c_s rho_s**2/a) (-d ni/ d rho)
c        + plasma source density
 
c        d_im_hat is the effective plasma diffusivity for impurities
 
c        eta_phi is the toroidal viscosity or momentum diffusivity
c
c       M ni d v_phi/ dt =
c error found 5/12/98  should be d (M ni v_phi)/ dt =
c -1/V(rho)_prime d /d rho V(rho)_prime
c      |grad rho|**2 ( M ni eta_phi_hat (c_s rho_s**2/a) (-d v_phi/ d rho)
c                    +    M v_phi ne d_hat (c_s rho_s**2/a) (-d ne/ d rho))
c        + toroidal momentum source density
c
c note that a convective part was added
c
c  eta_par_hat and eta_per_hat are diagnostic. See CAUTION on eta_phi_hat
c  at large gamma_p=  (-d v_phi/ d rho) /(c_s/a)
c
c  chie_e_gf is the eta_e mode electron transport which is te <-> ti
c  and mi <-> me isomorphic to eta_i (ITG) ion transport
c  with adiabatic electrons.
c  these mode obtain at high-k where the ions are adiabatic from
c  the gyro cut-off.
c  their wave numbers are sqrt(mi/me) larger than ITG modes and
c  since their frequencies are sqrt(mi/me) larger, they are not
c  rotationally shaer satbilized.
c  when xparam_gf(10).eq.0 xparam_gf(10)*chie_e_gf is added to
c  chie_gf and chie_e_gf is a diagnostic.
 
c input

c  eigen_gf = 0 use cgg eigenvalue solver (default)
c           = 1 use generalized tomsqz eigenvalue solver
c           = 2 use zgeev eigenvalue solver
c  nroot number of equations
c  iflagin(1:20) control flags
c   iflagin(1) 0 use ky=ky0; 1 use landau damping point
c   iflagin(2) 0. local w_d and k_par "2d"; 1 fit to trial function "3d"
c   iflagin(3) 0,1,and 2 fix up park low high beta and beta lim elong factor
c   iflagin(4) 0 trapped electron Waltz EoS 1 weiland EoS
c   iflagin(5) rms_theta 0:fixed; 1 inverse to q/2 ; 2 inverse to root q/2
c                        3: inverse to xparam(13)*(q/2-1)+1.
c              5 for retuned rms-theta
c  xparam(1:20) control parameters
c   xparam(1:2): idelta=xi*xparam(1)+xparam(2) nonadiabatic electron response
c   xparam(3) multiplier park_gf(high betae)/ park_gf(low betae) -1
c   xparam(6)+1. is enhancement of xnueff
c   xparam(7) coef of resistivity
c   xparam(8) cut off on rotational stabilization
c   xparam(9)+1. is shape (triangularity) enhancement to beta_crit
c   xparam(10) is high k electron mode enhancement
c   xparam(11:12) lamda parameters
c   xparam(13) rms_theta q-dependence
c   xparam(14)  adjustment to gamma_p avoiding negative viscosity
c   xparam(15)   (1+xparam(15)*reps trapped electron fraction
c   xparam(16) rms_theta shat dependence
c   xparam(17) ""
c   xparam(18) rms_theta betae dependence
c   xparam(19:20)  extra
c   xparam(21) 1 add impurity energy diffusivity to ion energy diffusivity
c   xparam(22) >0 keeps gamma_e from changeing spectrum
c   xparam(23) 1. kills kx**2 in k_m**2
c   xparam(24) exb damping model
c  ky0=k_theta*rho_s; k_theta= nq/r; normally 0.3
c  rms_theta width of phi**2 mode function for best fit near pi/3
c  rlti=a/L_Ti   a/L_f= sqrt(kappa) a d ln f / d rho
c  rlte=a/L_Te
c  rlne= a/L_ne
c  rlni= a/L_ni
c  rlnimp= a/L_nim
c  dil=1.-ni_0/ne_0  dilution
c  apwt = ni_0/ne_0
c  aiwt = nim_0/ne_0
c  taui=Ti/Te
c  rmin=r/a
c  rmaj=Rmaj/a
c  xnu=nu_ei/(c_s/a)
c  betae=neTe/(B**2/(8pi))  0 is electrostatic
c  shat= dlnr/drho used only for parallel dynamics part
c  alpha local shear parameter or MHD pressure grad (s-alpha diagram)
c  elong= local elongation or kappa
c  xwell amount of magnetic well xwell*min(alpha,1)
c  park=1  (0) is a control parameter to turn on (off) parallel motion
c       0.405 best at zero beta and 2.5x larger at high beta..see iflagin(3)
c  ghat=1  (0) is a control parameter to turn on (off) curvature drift
c  gchat=1 (0) is a control parameter to turn on (off) div EXB motion
c  adamp= radial mode damping exponent  1/4 < adamp < 3/4
c       0.25 from direct fit of simulations varying radial mode damping
c   but 0.75 is better fit to rlti dependence
c  alpha_star O(1-3)  gyyrobohm breaking coef for diamg. rot. shear
c  gamma_star ion diamagnetic rot shear rate in units of c_s/a
c  alpha_e O(1-3)   doppler rot shear coef
c  gamma_e    doppler rot shear rate in units of c_s/a
c  alpha_p 1.5  fit for parallel velocity shear effect at rmaj=3 and q=2
c  gamma_p    parallel velocity shear rate (-d v_phi/ drho) in units of c_s/a
c  kdamp model damping normally 0.
 
c output
 
c  yparam(20) output diagnostics
c kyf  value of ky used
c gamma   leading mode growth rate in c_s/a
c freq    leading mode freq rate in c_s/a
c ph_m    (e phi /T_e)/(rho_s/a)  saturation value
c d_hat    plasma diffusivity for ions
c d_im_hat    plasma diffusivity for impurities
c chii_hat ion energy diffusivity
c chie_hat electron energy diffusivity
c exch_hat anomalous e to i energy exchange
c eta_par_hat parallel component of toroidal momentum diffusivity
c eta_per_hat perpendicular    ""
c eta_phi_hat toroidal momentun diffusivity
 
c internal definitions
c nroot = number of equations,
c   nroot=12 full impurity dynamics
c   nroot=9 exb convective impurity dynamics
c   nroot=8 full pure plasma, nrout=6 (betae=0), nrout=5 (betae=0 and park=0)
c v(i)  12 solution vector
c   v(1)=n_i,v(2)=p_par,v(3)=p_per,v(4)=n_t,v(5)=p_t
c   v(6)=u_par, v(7)=n_u, v(8)=a_par
c   v(9)=n_im, v(10)=p_im_par,v(11)=p_im_per,v(12)=u_im_par
c -i*omega v(i)= sum_j amat(i,j) v(j) where omega=freq+xi*gamma
c quasineitrality is
c  (-idelta+(1-dil)*(1/taui)*(1-g1))*ph=f0*ph=(1-dil)*n_i-n_t-n_u
c  or (1-idelta-reps+(1-dil)*(1/taui)*(1-g1))*ph=f1*ph=(1-dil)*n_i-n_t
c  for betae=0
c k_m  inverse mixing length
 
c numerical controls and flags
c
c
c...Dimensions
c
c neq maximum number of equations
c ieq actual number used
 
c***********************************************************************
c***********************************************************************
c
      implicit none
c
      include 'glf.m'
c
c Glf is common block, which must contain all the _gf inputs and outputs
c
      character*1 jobvr, jobvl
      integer neq, iflagin(30), ilhmax, ilh, ikymaxtot,
     >  lprint, ieq, j1, j2, j, i, jmax,
     >  iar, ifail, jroot(4), itheta, iky, iky0, iroot, iglf
      double precision epsilon
      parameter ( neq = 12, epsilon = 1.D-34 )
c
      double precision pi, xparam(30),yparam(2*nmode),
     >  nroot,ky0,rms_theta,rlti,rlte,rlne,rlni,dil,taui,
     >  rmin,rmaj,q,rlnimp,amassimp,zimp,mimp,
     >  aikymax, aiky, apwt, aiwt,
     >  alpha_mode, gamma_mode, alpha_p, gamma_p,
     >  amassgas, chi_par_1, chi_per_1, x43,
     >  anorm, ave_g, ave_g0,
     >  ave_cos, ave_theta2, ave_kxdky2,
     >  dtheta, theta, phi2, ave_k_par, chk,
     >  alpha_n, yk, byk, fnn, fnp, fnf, fpn, fpp, fpf,
     >  xnu_col, amass_e, chkf, xadiabat,
     >  eta_par_hat, eta_per_hat, anorm_k, del_k,
     >  gamma_k_max, gamma_gross_net,
     >  xnu,betae,shat,alpha,elong,xwell,
     >  park,ghat,gchat,kdamp,
     >  adamp,alpha_star,gamma_star,alpha_e,gamma_e,
     >  kyf,gamma,freq,ph_m,d_hat,d_im_hat,
     >  chii_hat,chie_hat,exch_hat
      double COMPLEX xi, idelta,
     >  v(1:12), amat(1:12,1:12),
     >  n_i,p_par,p_per,n_t,p_t,u_par,n_u,a_par,ph,t_u,n_e,
     >  n_im,p_im_par,p_im_per
c     complex u_im_par
      double precision b0,g0,g1,g2,g3,g12,g23,
     >  b0i,g0i,g1i,g2i,g3i,g12i,g23i
      double COMPLEX f0,f1
      double precision k_par,ky,kx,k_per,k_m,
     >  w_s, w_d, w_d0, w_cd,
     >  reps,xnueff,betae0,k_par0
      double COMPLEX xmu,lamda_d,
     >  xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
      double precision gam_par,gam_per,x_par,x_per,xt_mhd,yt_mhd,
     >  th,tc,fh,fc,
     >  phi_norm,gamma_r
      double COMPLEX chknu,chknt,chknt2
      double precision phi_renorm,gamma_net
c
c...Declarations for eigenvaluesolver
c
      double precision zgamax
c
c... solver varaibles
c
      parameter ( iar=neq )
c     integer iai, ivr, ivi, intger(neq) ! if NAG solver f02ake used
c     parameter ( iai=neq, ivr=neq, ivi=neq )
      double precision ar(iar,neq), ai(iar,neq), rr(neq), ri(neq)
     &  , vr(iar,neq), vi(iar,neq)
      double precision br(iar,neq), bi(iar,neq), beta_tom(neq), ztemp1
 
      integer matz
      double precision fv1(neq),fv2(neq),fv3(neq)
c
c amat(i,j) = complex matrix A
c zevec(j) = complex eigenvector
c
      integer lwork
      parameter ( lwork=198 )
      double complex mata(iar,neq),cvr(iar,neq),cvl(iar,neq),w(neq)
      double complex work(lwork)
      double precision rwork(2*neq)
      double COMPLEX zevec(neq,neq), zomega(neq)
      double precision gammaroot(4),freqroot(4),phi_normroot(4)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c
c   return zero if no unstable modes
c
      if(ipert_gf.eq.1.and.ngrow_k_gf(0).eq.0)go to 888   
c
c...initialize variables
c
c inputs.........................................................

      do i=1,30
       iflagin(i)=iflagin_gf(i)
       xparam(i)=xparam_gf(i)
      enddo
 
      ilhmax=1
      ikymaxtot=ikymax_gf
c     if (xparam_gf(10).gt.0.) ilhmax=2
c
c If ETG modes included, then double ky spectrum
c Inside ky loop, ilh=1 low k ion modes and ilh=2 high k electron modes
c For ETG mode, use complete te <-> ti and mi <-> me isomorphism
c with adiabatic electron ITG mode then chii_hat will be electron 
c transport in gyrobohm electron units with T_i.
c chie_e_gf converted back to c_s*rho_s**2/a units and added 
c to chie_gf after ky loop
c
      if (xparam_gf(10).gt.0.) then
        ilhmax=2
        ikymaxtot=2*ikymax_gf
      endif
c
      nroot=nroot_gf
      ky0=xky0_gf
      rms_theta=rms_theta_gf
      rlti=rlti_gf
      rlte=rlte_gf
      rlne=rlne_gf
      rlni=rlni_gf
      rlnimp=rlnimp_gf
      dil=dil_gf
      apwt=apwt_gf
      aiwt=aiwt_gf
      taui=taui_gf
      rmin=rmin_gf
      rmaj=rmaj_gf
      q=q_gf
      xnu=xnu_gf
      betae=betae_gf
      shat=shat_gf
      alpha=alpha_gf
      elong=elong_gf
      xwell=xwell_gf
      park=park_gf
      ghat=ghat_gf
      gchat=gchat_gf
      adamp=adamp_gf
      alpha_star=alpha_star_gf
      gamma_star=gamma_star_gf
      alpha_e=alpha_e_gf
      gamma_e=gamma_e_gf
      alpha_mode=alpha_mode_gf
      gamma_mode=gamma_mode_gf
      alpha_p=alpha_p_gf
      gamma_p=gamma_p_gf
      kdamp=xkdamp_gf
      lprint=lprint_gf
      amassgas=amassgas_gf
      amassimp=amassimp_gf
      zimp=zimp_gf
      if(ipert_gf.eq.0)then
       do j=0,nmode
        ngrow_k_gf(j) = 0
       enddo
      endif
c
      idelta=0.D0
c     if(ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
 
c.................................................................
c
      if (lprint.gt.0) open(1)
      ieq  = nroot
c
      if (lprint.eq.99) then
      write(1,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ',
     >    ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
      write(1,*)'xnu,beta,shat,alpha,elong,xwell: ',
     >    xnu,betae,shat,alpha,elong,xwell
      write(1,*)'park, ghat, gchat: ',
     >    park, ghat, gchat
      write(1,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ',
     >    adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
      if (lprint.eq.98) then
      write(2,*) 'ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q: ',
     >    ky0,rms_theta,rlti,rlte,rlne,rlni,taui,rmin,rmaj,q
        write(2,*)'xnu,betae,shat,alpha,elong,xwell: ',
     >    xnu,betae,shat,alpha,elong,xwell
        write(2,*)'park, ghat, gchat: ',
     >    park, ghat, gchat
        write(2,*)'adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp: ',
     >    adamp,alpha_star,gamma_star,alpha_e,gamma_e,kdamp
 
      endif
c
      xi=(0.D0,1.D0)
      pi=atan2 ( 0.0D0, -1.0D0 )
 
c GLF model coefficients
 
      chi_par_1=2.D0*sqrt(2.D0/pi)
      chi_per_1=sqrt(2.D0/pi)
 
      gam_par=3.D0
      gam_per=1.D0
      x_par=2.D0
      x_per=3.D0/2.D0
 
      xmu=(0.80D0+.57D0*xi)
 
      xnu_par_par=(1.D0+xi)
      xnu_par_per=0.D0
      xnu_per_par=0.D0
      xnu_per_per=(1.D0+xi)
 
      lamda_d=(-0.7D0-0.80D0*xi)
      if(xparam(11).ne.0..or.xparam(12).ne.0.)
     >  lamda_d=xparam(11)-xi*xparam(12)
      x43=1.D0
      if (iflagin(4).eq.1) then
       lamda_d=5.D0/3.D0
       x43=4.D0/3.D0
      endif
c
c 3d trial wave function analysis
c
      if(iflagin(2).ge.1) then
       if(rms_theta.eq.0.) rms_theta=pi/3.D0
       if(iflagin(5).eq.1) rms_theta=rms_theta_gf*(2.D0/q_gf)
       if(iflagin(5).eq.2) rms_theta=rms_theta_gf*(2.D0/q_gf)**0.5D0
       if(iflagin(5).eq.3) rms_theta=rms_theta_gf/
     >      (xparam(13)*(q_gf/2.D0-1.D0)+1.D0)
     >  /sqrt(1.D0+xparam(16)*(shat_gf**2-1.D0)+xparam(17)*
     >    (shat_gf-1.D0)**2)
     >  /(1.D0+xparam(18)*sqrt(betae/.006D0))
       if(iflagin(5).eq.4) rms_theta=rms_theta_gf/
     >      (xparam(13)*(q_gf/2.D0-1.D0)+1.D0)
     >  /sqrt(1.D0+xparam(16)*(shat_gf**2-1.0D0)+xparam(17)*
     >    (shat_gf-1.D0)**2+xparam(19)*(alpha-0.5D0)**2.D0)
     >  /(1.D0+xparam(18)*sqrt(betae/.006D0))
       if(iflagin(5).eq.5) rms_theta=rms_theta_gf/
     >      (xparam(13)*((q_gf/2.D0)-1.D0)+1.D0)
     >  /sqrt(1.D0+xparam(16)*((shat_gf-
     >  xparam(19)*alpha)**2-0.5D0)+xparam(17)*
     >  (shat_gf-xparam(19)*alpha-0.5D0)**2)/
     >  taui_gf**0.25D0
c along the filed line physics with wave function
c phi=exp(-theta**2/(4.*rms_theta**2))=W_even
c ave_F= [int(0 to inf) F phi**2 d_theta]/[int(0 to inf)  phi**2 d_theta]
c ave_theta**2=(rms_theta)**2
 
c phi, densities and pressures are even functions of theta
c the odd functions like u_par can be represented by
c W_odd= W*i*theta/rms_theta*W_even
c then for W=-1, the k_par operator = i*k_par=1/(rmaj*q) d/dtheta
c becomes 1/(rmaj*q)/(2.*rms_theta)*park  (park=1) in every equation
c ie ave_k_par=1/(2.*rms_theta)
c park is tuned to best fit.
c
c parallel velocity shear gamma_p breaks parity so wave functions
c become mixed W_even->(1-xi*gamma_p*alpha_n*theta/rms_theta)*W_even
c
c gamma_p*alpha_n mustbe dimensionless and independent of norm a
c hence alpha_n=(rmaj/3)*alpha_p since gamma_p is in units
c of c_s/a  and rmaj=rmajor/a. Since in the slab limit where
c where we can have the parallel shear drive, rmaj enters only with
c parameter rmaj*q, we further assume
c alpha_n=(rmaj/3.)*(q/2)*alpha_p as the appropriate scaling
c q=2 and rmaj=3 are the norm points for alpha_p=1.5
c For the extreme toroidal limit q-> infinity where rmaj and
c q are not product associated, we will lose the instability.
c
c to first order in gamma_p*alpha_n
c this leads to a weighting  factor [gamma_p*alpha_n] in the
c xi*ky*ph1*gamma_p linear drive and in the
c eta_phi_hat=conjg(u_par)*(-xi*ky*ph)/gamma_p toroidal vocosity.
c
c the correct dependence gamma-gamma0 going like gamma_p**2 is found
c but QLT eta_phi_hat goes like gamma_p**2 also
c CAUTION: this is worked out only to first order in gamma_p*alpha_n
c small.  It seems likely that there is a higher order saturation factor
c something like 1/(1+(gamma_p*alpha_n)**2) in  eta_phi_hat
c
c doppler (EXB) rotational shear also breaks parity
c thus there should a term egamma*alpha_n_e added to gamma_p*alpha_n
c but it is unclear how to weight alpha_n_e compared to alpha_n
c
c see R.R.Dominguez and G.M Staebler Phys. Fluids B5 (1993) 3876
c for a discussion of QLT theory of anomalous momentum transport in
c slab geometry
 
c compute weight factors
c fix later so these are computed only once per j grid point
 
c      ave_theta2=rms_theta**2
 
      anorm=0.D0
      ave_g=0.D0
      ave_g0=0.D0
      ave_cos=0.D0
      ave_theta2=0.D0
      ave_kxdky2=0.D0
 
      dtheta=4.D0*rms_theta/100.D0
      theta=0.D0
c
      do itheta=1,100
       theta=theta+dtheta
       phi2=exp(-theta**2/(2.D0*rms_theta**2))
       anorm=anorm+phi2*dtheta
       ave_theta2=ave_theta2+
     >  theta**2*phi2*dtheta
       ave_g=ave_g +
     > (-xwell*min(1.D0,alpha)+cos(theta)+
     >    (shat*theta-alpha*sin(theta))*sin(theta))*phi2*dtheta
       ave_g0=ave_g0 + phi2*dtheta
       ave_kxdky2=ave_kxdky2+
     >  (abs(shat*theta-alpha*sin(theta)))**2*phi2*dtheta
       ave_cos=ave_cos +
     >  cos(theta)*phi2*dtheta
      enddo
c
      ave_theta2=ave_theta2/anorm
      ave_g=ave_g/anorm
      ave_g0=ave_g0/anorm
      ave_kxdky2=ave_kxdky2/anorm
      ave_cos=ave_cos/anorm
c 
      ave_k_par=1/(2.D0*rms_theta)
 
      chk=abs(ave_theta2-rms_theta**2)/rms_theta**2
      if (chk.gt..02) write (6,*) 'chk:', chk
 
      alpha_n=(rmaj/3.D0)*(q/2.D0)*alpha_p
 
      if(lprint.eq.2) then
       write(6,*) 'rms_theta,chk :', rms_theta, chk
       write(6,*) 'ave_theta2,ave_g,ave_k_par,ave_cos:',
     >   ave_theta2,ave_g,ave_k_par,ave_cos
      endif
      endif
      if(iflagin(2).eq.0) then
       shat=1.D0
       ave_theta2=1.D0
       ave_g=1.D0
       ave_g0=1.D0
       ave_kxdky2=1.D0
       ave_k_par=1.D0
       ave_cos=1.D0
      endif
c
c start ky loop 
c first half ITG, second half high-k ETG ... each with ikymax_gf modes
c ilh=1 low k ion modes  ilh=2 high k electron modes

      do iky0=1,ikymaxtot
c
cgms      iky=iky0
      iky = ikymax_gf+1-iky0
      ilh=1
c
c offset iky if in high-k range and set ilh=2
c
      if (iky0.gt.ikymax_gf) then
cgms         iky=iky0-ikymax_gf
         iky = ikymaxtot+1-iky0
         ilh=2
      endif
c
      if (ilh.eq.2) then
       nroot=6
       ieq=nroot
       xnu=0.D0
       betae=1.D-6
       rlte=rlti_gf
       rlti=rlte_gf
       rlne=rlni_gf
       rlni=rlne_gf
       rlnimp=epsilon
       dil=1.D0-1.D0/(1.D0-dil_gf)
       apwt=1.D0
       aiwt=0.D0
       taui=1.D0/taui_gf
       rmin=epsilon
       xparam(7)=0.D0
       xparam(6)=-1.D0
       alpha_star=0.D0
       alpha_e=0.D0
       alpha_p=0.D0
       alpha_n=0.D0
       alpha_mode=0.D0
c check this for current driven mode
      endif
c
      idelta=0.D0
      if (ilh.eq.1) idelta=xi*xparam(1)+xparam(2)
c
c logarithmic ky grid
c
      if(ikymax_gf.gt.1) then
       aikymax=ikymax_gf
       aiky=iky
       yk=aiky/aikymax
 
       byk=log(xkymax_gf/xkymin_gf)/(1.D0-1.D0/aikymax)
       ky=xkymax_gf*exp(byk*(yk-1.D0))
      endif
      if(ikymax_gf.eq.1) then
c     ky=sqrt(2.*taui)/rlti/(rmaj*q)
c     from w_star_ti=v_i_th*k_par
c  possible physics basis of q (ie current) scaling ..to be determined
       if(iflagin(1).eq.0) ky=ky0
       if(iflagin(1).eq.1) ky=ky0*sqrt(taui)*(3.D0/rlti)*
     >    (3.D0*2.D0/rmaj/q)
       if(iflagin(1).eq.2) ky=ky0*(2.D0/q)
       if(ky0.eq.0.) ky=0.3D0
      endif
 
 
      kyf=ky
c
 
      kx=ky*sqrt(ave_kxdky2)
      k_per=sqrt(ky**2+kx**2)
      k_m=sqrt(ky**2+(1.D0-xparam(23))*kx**2) ! inverse mixing length model
 
       do iroot=1,4
        gammaroot(iroot)=0.D0
        freqroot(iroot)=0.D0
        phi_normroot(iroot)=0.D0
       enddo
        d_hat=0.D0
        d_im_hat=0.D0
        chie_hat=0.D0
        chii_hat=0.D0
        exch_hat=0.D0
        eta_par_hat=0.D0
        eta_per_hat=0.D0
        jroot(1)=0
        jroot(2)=0
        jroot(3)=0
c
c skip this k for perturbation if last call was stable
c
      if(ipert_gf.eq.1.and.ngrow_k_gf(iky0).eq.0)go to 777 
c skip this k if the previous k was stable and 4 k's have been done      
      if(iky.lt.ikymax_gf-4.and.ngrow_k_gf(iky0-1).eq.0)go to 777
 
c primary ions
 
      b0=taui*k_per**2
 
c     Pade aproximates...may use gamma functions later
      g0=1.D0
      g1=1.D0/(1+b0)
      g2=1.D0/(1+b0)*g1
      g3=1.D0/(1+b0)*g2
 
      g12=(g1+g2)/2.D0
      g23=(g2+g3)/2.D0
 
c impurity ions
 
      b0i=taui*k_per**2*amassimp/amassgas/zimp**2
 
c     Pade aproximates...may use gamma functions later
      g0i=1.D0
      g1i=1.D0/(1+b0i)
      g2i=1.D0/(1+b0i)*g1i
      g3i=1.D0/(1+b0i)*g2i
 
      g12i=(g1i+g2i)/2.D0
      g23i=(g2i+g3i)/2.D0
 
      mimp=amassimp/amassgas
 
 
      w_s=ky
      w_d=(ghat*2.D0/rmaj)*ky*ave_g
      w_d0=(ghat*2.D0/rmaj)*ky*ave_g0
      w_cd=(gchat*2.D0/rmaj)*ky*ave_g
 
      k_par=park/(rmaj*q)*ave_k_par*sqrt((1.D0+elong**2)/2.D0)
 
c     sqrt((1.+elong**2)/2.) to get higher beta_crit prop to k_par**2
c     roughly same as betae-> betae/((1.+elong**2)/2.)
c     physically like shortening the connection length to good curv.
 
      if (iflagin(3).eq.2) then
       betae=betae_gf/(1.D0+xparam(3))**2/(1.D0+xparam(9))
      endif
 
 
      if (iflagin(3).eq.1) then
       k_par=park/(rmaj*q)*ave_k_par
       betae=betae_gf/(1.D0+xparam(3))**2/(1.D0+xparam(9))
     >      /((1.D0+elong**2)/2.D0)
      endif
 
c     we put the park enhancement directy into betae
c     ie park=.5 best at low beta and 2.5x.5=1.25 at high beta
 
c     option iglagin(3)=1 puts beta_crit elongation enhancement
c     directly into betae
 
c     option iflagin(3)=2 puts beta_crit elongation factor into
c     the connection length
 
c     an extra shape  factor 2 (triangularity) enhancement
c     is optained by (1.+xparam(9))=2.
c     if w_d negative flip sign of dissipative parts
      if(w_d.lt.0.) then
c error 12/21       lamda_d=-conjg(lamda_d)
       lamda_d=conjg(lamda_d)
 
       xmu=-conjg(xmu)
 
       xnu_par_par=-conjg(xnu_par_par)
       xnu_par_per=-conjg(xnu_par_per)
       xnu_per_par=-conjg(xnu_per_par)
       xnu_per_per=-conjg(xnu_par_per)
      endif
 
      reps=(1.D0+xparam(15))*
     >   sqrt((rmin/rmaj)*(1.D0+ave_cos)/(1.D0+(rmin/rmaj)*ave_cos))
 
      if(nroot.le.3) reps=0.D0
 
c fix trapped eletron MHD limit
c 3/4*reps*(1+rlte) + xt_mhd*(1-reps)+yt_mhd*rlte=1+rlte
c solve for xt_mhd and yt_mhd
c 3/4*reps+yt_mhd=1; 3/4*reps+xt_mhd*(1-reps)=1
 
      yt_mhd=(1-x43*(3.D0/4.D0)*reps)
      xt_mhd=(1.D0-x43*(3.D0/4.D0)*reps)/(1.D0-reps)
 
c collision detrapping retrapping model
 
      xnueff=(1.D0+xparam(6))*xnu/(reps**2+1.D-6)
 
c very difficult get xnueff correct hince add enhancement factor
c and fit to Beer or GKS
 
      th=4.08D0
      tc=0.918D0
      fh=0.184D0
      fc=0.816D0
 
      fnn=xnueff*((th/tc**(3.D0/2.D0))-(tc/th**(3.D0/2.D0)))/(th-tc)
      fnp=xnueff*(3.D0/2.D0)*
     >   ((1.D0/th)**(3.D0/2.D0)-(1.D0/tc)**(3.D0/2.D0))/(th-tc)
      fnf=xnueff*((fh/th**(3.D0/2.D0))+(fc/tc**(3.D0/2.D0)))
 
      fpn=xnueff*(2.D0/3.D0)*
     >   ((th/tc**(1.D0/2.D0))-(tc/th**(1.D0/2.D0)))/(th-tc)
      fpp=xnueff*((1.D0/th)**(1.D0/2.D0)-(1.D0/tc)**(1.D0/2.D0))/(th-tc)
      fpf=xnueff*(2.D0/3.D0)*((fh/th**(1.D0/2.D0))+(fc/tc**(1.D0/2.D0)))
 
c  collisional modes added with xnu_col
c  must fix for atomic mass dependence other than deuterium
      xnu_col=xparam(7)
      amass_e=2.7D-4*(2.D0/amassgas)
 
c check adiabatic property that chkf should be 1.0  (finite xnu)
 
      chkf=(fnn*fpp-fpn*fnp)/((fnf*fpp-fpf*fnp)+epsilon)
      if (lprint.eq.2) write(6,*) 'chkf:', chkf
 
      if(neq.le.3) reps=0.D0
 
      f0=-idelta+(1.D0-dil)*apwt*(1/taui)*(1.D0-g1)
     >           +zimp**2*aiwt*(1/taui)*(1.D0-g1i)
      f1=1.D0-reps + f0
 
      xadiabat=0.D0
      if(nroot.le.6) then
        betae=0.D0
        f0=f1
        xadiabat=1.D0
      endif
      if(nroot.le.5) k_par=0.D0
 
      betae0=betae+epsilon
      k_par0=k_par+epsilon
 
      if (lprint.eq.98) then
        write(2,*) 'ky,g1,g2,g3,g12,g23,w_s,w_d,w_cd: ',
     >    ky,g1,g2,g12,g23,w_s,w_d,w_cd
        write(2,*) 'f0,taui,k_par,reps: ',
     >    f0,taui,k_par,k_per,reps
        write(2,*) 'chi_par_1,chi_per_1,gam_par,gam_per:',
     >    chi_par_1,chi_per_1,gam_par,gam_per
        write(2,*) 'x_par,x_per,xmu:',
     >    x_par,x_per,xmu
        write(2,*) 'xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per:',
     >    xnu_par_par,xnu_par_per,xnu_per_par,xnu_per_per
        write(2,*) 'lamda_d,betae,xadiabat:',
     >    lamda_d,betae,xadiabat
        write(2,*) 'yt_mhd,xt_mhd:',
     >    yt_mhd,xt_mhd
 
      endif
c 
c matrix in order
c note ph=(n_i-n_t-n_u)/f0 results in (i,1)-(i,4)-(i,7) parts
c
c n_i equ #1
 
      amat(1,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,2)=
     >  +xi*w_d*taui*0.5D0
 
      amat(1,3)=
     >  +xi*w_d*taui*0.5D0
 
      amat(1,4)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,5)= 0.D0
 
      amat(1,6)=
     >  -xi*k_par
 
      amat(1,7)= -(-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,8)= 0.D0
 
      amat(1,9)= aiwt*zimp*
     >  (-xi*w_s*((rlni-rlti)*g1+rlti*g2)+xi*w_cd*g12)/f0
 
      amat(1,10)=0.D0
 
      amat(1,11)=0.D0
 
      amat(1,12)=0.D0
 
c p_par equ #2
 
      amat(2,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
     >  +k_par*chi_par_1
     >  -(xi*w_d*taui*3.D0/2.D0-w_d*taui*xnu_par_par)
     >  -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)
 
      amat(2,2)=
     >  -k_par*chi_par_1
     >  +xi*w_d*taui*x_par +
     >  (xi*w_d*taui*3.D0/2.D0-w_d*taui*xnu_par_par)
 
      amat(2,3)=
     >  (xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)
 
      amat(2,4)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,5)=0.D0
 
      amat(2,6)=
     >  -xi*gam_par*k_par
 
      amat(2,7)= -(-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,8)=0.D0
 
      amat(2,9)= aiwt*zimp*
     >  (-xi*w_s*(rlni*g1+rlti*g2)+xi*x_par*w_cd*g12)/f0
 
      amat(2,10)=0.D0
 
      amat(2,11)=0.D0
 
      amat(2,12)=0.D0
 
c p_per equ #3
 
      amat(3,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
     >  +k_par*chi_per_1
     >  -(xi*w_d*taui-w_d*taui*xnu_per_per)
     >  -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_per_par)
 
 
      amat(3,2)=
     >  +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)
 
      amat(3,3)=
     >  -k_par*chi_per_1
     >  +xi*w_d*taui*x_per   +(xi*w_d*taui-w_d*taui*xnu_per_per)
 
      amat(3,4)=
     > -(-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,5)=0.D0
 
      amat(3,6)=
     >  -xi*gam_per*k_par
 
      amat(3,7)=
     >  -(-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,8)=0.D0
 
      amat(3,9)= aiwt*zimp*
     >   (-xi*w_s*((rlni-rlti)*g2+2.D0*rlti*g3)+xi*x_per*w_cd*g23)/f0
 
      amat(3,10)=0.D0
 
      amat(3,11)=0.D0
 
      amat(3,12)=0.D0
 
c n_t equ #4
 
      amat(4,1)=(1.D0-dil)*apwt*
     >  (-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >  -(1.D0-dil)*apwt*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(4,2)=0.D0
 
      amat(4,3)=0.D0
 
      amat(4,4)=-(-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >  -((1.D0-reps)*fnn)
     >  -(-(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat))
 
      amat(4,5)=
     >  -xi*w_d*x43*3.D0/4.D0
     >  -((1.D0-reps)*fnp)
 
      amat(4,6)=0.D0
 
      amat(4,7)=-(-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >   -(-reps*fnf)
 
      amat(4,8)=0.D0
 
      amat(4,9)=aiwt*zimp*
     >   (-xi*w_s*rlne*reps*g0+xi*x43*3.D0/4*w_cd*reps*g0)/f0
     >   -aiwt*zimp*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(4,10)=0.D0
 
      amat(4,11)=0.D0
 
      amat(4,12)=0.D0
 
c p_t equ #5
 
      amat(5,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >  -(1.D0-dil)*apwt*(-reps*fpf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(5,2)=0.D0
 
      amat(5,3)=0.D0
 
      amat(5,4)=
     >  -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >            +xi*w_d*lamda_d
     >  -((1.D0-reps)*fpn)
     >  -(-(-reps*fpf*(1.D0-reps)*g0/f0*xadiabat))
 
      amat(5,5)=
     >  -xi*w_d*x43*5.D0/4.D0-xi*w_d*lamda_d
     >  -((1.D0-reps)*fpp)
 
      amat(5,6)=0.D0
 
      amat(5,7)=
     >  -(-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >  -(-reps*fpf)
 
      amat(5,8)=0.D0
 
      amat(5,9)= aiwt*zimp*
     >  (-xi*w_s*(rlni+rlte)*reps*g0+xi*x43*5.D0/4*w_cd*reps*g0)/f0
     >  -aiwt*zimp*(-reps*fpf*(1.D0-reps)*g0/f0*xadiabat)

      amat(5,10)=0.D0

      amat(5,11)=0.D0

      amat(5,12)=0.D0

c u_par equ #6
 
      amat(6,1)=(1.D0-dil)*apwt*
     >   (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0
     >   -(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
 
      amat(6,2)=-xi*k_par*taui
 
      amat(6,3)=0.D0
 
      amat(6,4)=
     >  -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)
     >  -(-(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
 
      amat(6,5)=0.D0
 
      amat(6,6)=
     > +xi*w_d*(gam_par+gam_per)/2.D0 -w_d*xmu
 
      amat(6,7)=
     >  -(-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0)
     >  -(-(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
     >  -(betae/2.D0)*xi*k_par*(2.D0/betae0)/(1.D0-reps)
 
      amat(6,8)=
     >  -(betae/2.D0)*(-xi*w_s*(rlni*g1+rlti*g2))
     >  -(betae/2.D0)*(-xi*w_s*rlne)
     >  +amass_e*xnu*xnu_col*k_per**2
     >  -(betae/2.D0)*(-2.D0/betae0*amass_e*xnu*xnu_col*k_per**2)
c note there is no double counting in last two terms
 
      amat(6,9)=aiwt*zimp*
     >  (-xi*k_par*g1/f0-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1/f0
     >  -(betae/2.D0)*(-xi*k_par*(2.D0/betae0)*g0)/f0)
 
      amat(6,10)=0.D0
 
      amat(6,11)=0.D0
 
      amat(6,12)=0.D0
 
c n_u equ #7
 
      amat(7,1)=(1.D0-dil)*apwt*
     >  (-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  +(1.D0-dil)*apwt*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(7,2)=0.D0
 
      amat(7,3)=0.D0
 
      amat(7,4)=
     >  -(-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  +((1.D0-reps)*fnn)
     >  +(-(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat))
 
      amat(7,5)=0.D0
     >  +((1.D0-reps)*fnp)
 
      amat(7,6)=-xi*k_par
 
      amat(7,7)=
     >  -(-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  -xi*w_d*xt_mhd
     >  +(-reps*fnf)
 
      amat(7,8)=
     >  -xi*k_par*(-k_per**2)
     >  -xi*w_d*yt_mhd*(w_s*(betae/2.D0)/k_par0*rlte)
 
      amat(7,9)=aiwt*zimp*
     >  (-xi*w_s*rlne*(1.D0-reps)*g0+xi*w_cd*
     >  (1.D0-x43*(3.D0/4.D0)*reps)*g0)/f0
     >  +aiwt*zimp*(-reps*fnf*(1.D0-reps)*g0/f0*xadiabat)
 
      amat(7,10)=0.D0
 
      amat(7,11)=0.D0
 
      amat(7,12)=0.D0
 
c a_par equ #8
 
      amat(8,1)=(1.D0-dil)*apwt*(-xi*k_par*(2.D0/betae0)*g0/f0)
 
      amat(8,2)=0.D0
 
      amat(8,3)=0.D0
 
      amat(8,4)=-(-xi*k_par*(2.D0/betae0)*g0/f0)
 
      amat(8,5)=0.D0
 
      amat(8,6)=0.D0
 
      amat(8,7)=-(-xi*k_par*(2.D0/betae0)*g0/f0)
     >  +xi*k_par*(2.D0/betae0)/(1.D0-reps)
 
      amat(8,8)=-xi*w_s*rlne
     >  -(2.D0/betae0)*amass_e*xnu*xnu_col*(k_per**2)
 
      amat(8,9)=aiwt*zimp*(-xi*k_par*(2.D0/betae0)*g0/f0)
 
      amat(8,10)=0.D0
 
      amat(8,11)=0.D0
 
      amat(8,12)=0.D0
 
c n_im equ #9
 
      amat(9,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,2)=0.D0
 
      amat(9,3)=0.D0
 
      amat(9,4)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,5)= 0.D0
 
      amat(9,6)= 0.D0
 
      amat(9,7)= -(-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,8)= 0.D0
 
      amat(9,9) = aiwt*zimp*
     >  (-xi*w_s*((rlnimp-rlti)*g1i+rlti*g2i)+xi*w_cd*g12i)/f0
 
      amat(9,10)=
     >  +xi*w_d*taui*0.5D0/zimp
 
      amat(9,11)=
     >  +xi*w_d*taui*0.5D0/zimp
 
      amat(9,12)=
     >  -xi*k_par
 
c pim_par equ #10
 
      amat(10,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,2)=0.D0
 
      amat(10,3)=0.D0
 
      amat(10,4)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,5)=0.D0
 
      amat(10,6)=0.D0
 
      amat(10,7)= -(-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
 
      amat(10,8)=0.D0
 
      amat(10,9)= aiwt*zimp*
     >   (-xi*w_s*(rlnimp*g1i+rlti*g2i)+xi*x_par*w_cd*g12i)/f0
     >   +k_par*chi_par_1/sqrt(mimp)
     >   -(xi*w_d*taui*3.D0/2.D0-w_d*taui*xnu_par_par)/zimp
     >   -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)/zimp
 
      amat(10,10)=
     >  -k_par*chi_par_1/sqrt(mimp)
     >  +xi*w_d*taui*x_par/zimp +(xi*w_d*taui*3.D0/2.D0
     >  -w_d*taui*xnu_par_par)/zimp
 
      amat(10,11)=
     >  (xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_par_per)/zimp
 
      amat(10,12)=
     >  -xi*gam_par*k_par
 
c pim_per equ #11
 
      amat(11,1)= (1.D0-dil)*apwt*
     >  (-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
 
      amat(11,2)= 0.D0
 
      amat(11,3)= 0.D0
 
      amat(11,4)=
     >  -(-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
 
      amat(11,5)=0.D0
 
      amat(11,6)=0.D0
 
      amat(11,7)=
     >  -(-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
 
      amat(11,8)=0.D0
 
      amat(11,9)= aiwt*zimp*
     >  (-xi*w_s*((rlnimp-rlti)*g2i+2.D0*rlti*g3i)+
     >  xi*x_per*w_cd*g23i)/f0
     >  +k_par*chi_per_1/sqrt(mimp)
     >  -(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp
     >  -(xi*w_d*taui*1.D0/2.D0-w_d*taui*xnu_per_par)/zimp
 
      amat(11,10)=
     >  +(xi*w_d*taui*1/2-w_d*taui*xnu_per_par)/zimp
 
      amat(11,11)=
     >  -k_par*chi_per_1/sqrt(mimp)
     >  +xi*w_d*taui*x_per/zimp 
     >  +(xi*w_d*taui-w_d*taui*xnu_per_per)/zimp
 
      amat(11,12)=
     >  -xi*gam_per*k_par
 
c uim_par equ #12
cgms 5/21/99 added gamma_p to amat(12,1),amat(12,4),amat(12,7),amat(12,9)
c    added xnu_col term to amat(12,8)
c    fixed mimp factor in amat(12,12)
 
      amat(12,1)=(1.D0/mimp)*(1.D0-dil)*apwt*
     >  ((-xi*k_par*g1i/f0)*zimp
     >  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0
     >  -(betae/2.D0)*zimp*(-xi*k_par*(2.D0/betae0)*g0i)/f0)

      amat(12,2)=0.D0
 
      amat(12,3)=0.D0
 
      amat(12,4)=
     > -(1.D0/mimp)*(-xi*k_par*g1i/f0)*zimp
     > -(1.D0/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)
     > -(1.D0/mimp)*(-(betae/2.D0)*(-xi*k_par
     > *(2.D0/betae0)*g0i)/f0)*zimp

      amat(12,5)=0.D0
 
      amat(12,6)=0.D0
 
       amat(12,7)=
     > -(1.D0/mimp)*(-xi*k_par*g1i/f0)*zimp
     > -(1.D0/mimp)*(-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0)
     > -(1.D0/mimp)*(-(betae/2.D0)*
     >  (-xi*k_par*(2.D0/betae0)*g0i)/f0)*zimp
     > -(1.D0/mimp)*(betae/2.D0)*xi*
     >  k_par*(2.D0/betae0)/(1.D0-reps)*zimp
 
      amat(12,8)=
     >  -(1.D0/mimp)*(betae/2.D0)*(-xi*w_s*(rlnimp*g1i+rlti*g2i))
     >  -(1.D0/mimp)*(betae/2.D0)*(-xi*w_s*rlne)*zimp
     >  +(1.D0/mimp)*zimp*amass_e*xnu*xnu_col*k_per**2
     >  -(1.D0/mimp)*(betae/2.D0)*
     >   (-2.D0/betae0*amass_e*xnu*xnu_col*k_per**2)*zimp

      amat(12,9)=(1.D0/mimp)*aiwt*zimp*
     >   ((-xi*k_par*g1i/f0)*zimp
     >  -xi*ky*gamma_p*(-gamma_p*alpha_n)*g1i/f0
     >  -(betae/2.D0)*zimp*(-xi*k_par*(2.D0/betae0)*g0i)/f0)

      amat(12,10)=-(1.D0/mimp)*xi*k_par*taui
 
      amat(12,11)=0.D0
 
      amat(12,12)= (1.D0/mimp)*(xi*w_d*(gam_par+gam_per)
     >  /2.D0/zimp -w_d*xmu/zimp)
c
c put in rot shear stabilization and possible source of gyrobohm breaking
c and model damping kdamp
c
c***********************************************************************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c solve 12x12 complex
c -xi*omega*v(i)=sum_j amat(i,j)*v(j)  omega=freq+xi*gamma
c upto nroot
c order with max gamma and find eigenvector v(i) with ant fixed norm.
c
c...Fill matricies for eigenvalue equation
c
      do j1=1,neq
        rr(j1) = 0.0D0
        ri(j1) = 0.0D0
        do j2=1,neq
          ar(j1,j2) = AREAL(  amat(j1,j2) )
          ai(j1,j2) = aimag( amat(j1,j2) )
c...test tmp
c         ai(j1,j2) = 0.0
c         ar(j1,j2) = 0.0
c         if (j1.eq.j2) ar(j1,j2)=j1
c
          vr(j1,j2) = 0.0D0
          vi(j1,j2) = 0.0D0
        enddo
      enddo
c
c...diagnostic output
c
      if ( lprint .gt. 6 ) then
        write (1,*)
        write (1,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (1,192) (ar(j1,j2),j2=1,neq)
        enddo
c
        write (1,*)
        write (1,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (1,192) (ai(j1,j2),j2=1,neq)
        enddo
 192    format (1p8e10.2)
 193    format (1p8e12.4)
      endif
c
c..find the eigenvalues and eigenvectors 
c
c.. eigen_gf = 0 use cgg solver (default)
c..          = 1 use tomsqz solver
c..          = 2 use zgeev solver
c.. not longer used:
c        call f02ake( ar,iar,ai,iai,ieq,rr,ri,vr,ivr,vi,ivi,
c     >               intger, ifail )
c
        ifail = 0
c
        if (eigen_gf .eq. 2 ) then
c
        jobvl = 'N'
        jobvr = 'V'
        do j1=1,neq
         do j2=1,ieq
           mata(j1,j2) = cmplx(ar(j1,j2),ai(j1,j2))
         enddo
        enddo
c
        call zgeev(jobvl,jobvr,ieq,mata,neq,w,cvl,neq,cvr,
     &             neq,work,lwork,rwork,ifail)
        do j1=1,neq
         rr(j1) = areal(w(j1))
         ri(j1) = aimag(w(j1))
         do j2=1,ieq
           vr(j1,j2) = areal(cvr(j1,j2))
           vi(j1,j2) = aimag(cvr(j1,j2))
         enddo
        enddo
c
        elseif (eigen_gf .eq. 1 ) then
c
        do j2=1,neq
           do j1=1,neq
              bi(j1,j2)=0.0D0
              if(j1.eq.j2) then
                 br(j1,j2)=1.0D0
              else
                 br(j1,j2)=0.0D0
              endif
           enddo
        enddo
c
        call r8tomsqz(neq,ieq,ar,ai,br,bi, rr,ri,beta_tom, vr,vi, ifail)
c
        do j1=1,ieq
           ztemp1 = beta_tom(j1)
           if ( abs(beta_tom(j1)) .lt. epsilon ) ztemp1 = epsilon
           rr(j1)=rr(j1) / ztemp1
           ri(j1)=ri(j1) / ztemp1
        enddo
c
        else
c
        matz=1
c       write(*,*) 'neq = ',neq
c       write(*,*) 'ieq = ',ieq
c       write(*,*) 'matz = ',matz
c       write (*,*) ' ar(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (ar(j1,j2),j2=1,neq)
c       enddo
c       write (*,*) ' ai(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (ai(j1,j2),j2=1,neq)
c       enddo
c
        call cgg_glf(neq,ieq,ar,ai,rr,ri,matz,vr,vi,fv1,fv2,fv3,ifail)
c
c       write (*,*) ' wr(j1) and wi(j1)'
c       do j1=1,neq
c         write (*,193) rr(j1), ri(j1)
c       enddo
c       write (*,*) ' zr(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (vr(j1,j2),j2=1,neq)
c       enddo
c       write (*,*) ' zi(j1,j2)  j2 ->'
c       do j1=1,neq
c         write (*,193) (vi(j1,j2),j2=1,neq)
c       enddo

        endif
c
        if ( lprint .gt. 1 ) then
          write (1,*) ifail,' = ifail routine '
        endif
c
c..print eigenvalues
c
        if ( lprint .gt. 6 ) then
          write (1,121)
          do j=1,ieq
            write (1,122) rr(j), ri(j)
          enddo
 121      format (/' Solution of the eigenvalue equations'
     &     /t4,'real   ',t18,'imag   ')
 122      format (1p2e14.5)
        endif
c
c...Store the complex eigenvectors and eigenvalues
c...Note the routines here solve A.v = lambda v
c...but that the equation solved is A.v = -i omega v
c...The i-th column of v is the i-th eigenvector
c
        do j1=1,ieq
           zomega(j1) = xi*(rr(j1)+xi*ri(j1))
          do j2=1,ieq
            zevec(j2,j1) = vr(j2,j1) + xi*vi(j2,j1)
          enddo
        enddo
c
        if ( lprint .gt. 6 ) then
          write (6,123)
          do j=1,ieq
            write (6,122) AREAL(zomega(j)), aimag(zomega(j))
          enddo
 123      format (/' Multiplied by i: '
     &     /t4,'zomegar',t18,'zomegai')
        endif
        do iroot=1,4
c
c
c..save growth rates and frequencies in real variables
c
        zgamax = 0.0D0
ctemp
        zgamax = -1.D10
        jmax=0
        gamma=0.D0
        do j=1,ieq
         if(j.ne.jroot(1).and.j.ne.jroot(2).and.j.ne.jroot(3)) then
          if (aimag(zomega(j)).gt. zgamax) then
            zgamax = aimag(zomega(j))
            jmax=j
          endif
         endif
        enddo
c
        IF(jmax.ne.0) THEN
           freq = AREAL( zomega(jmax) )
           gamma = aimag( zomega(jmax) )
        END IF
c
c skip stable modes
c        if(gamma.lt.0.D0)go to 775

         jroot(iroot)=jmax
 
ctemp        if(zgamax.lt.zepsqrt) gamma=0.
 
        if (jmax.ne.0) then
         gammaroot(iroot)=gamma
         freqroot(iroot)=freq
 
 
        do j=1,12
         v(j)=0.D0
        enddo
        do j=1,ieq
          v(j) = zevec(j,jmax)
        enddo
c
c***********************************************************************
 
      n_i=0.D0
      p_par=0.D0
      p_per=0.D0
      n_t=0.D0
      p_t=0.D0
      u_par=0.D0
      n_u=0.D0
      a_par=0.D0
      n_im=0.D0
      p_im_par=0.D0
      p_im_per=0.D0
c     u_im_par=0.
 
      t_u=0.D0
 
 
      n_i=v(1)
      p_par=v(2)
      p_per=v(3)
      if(ieq.ge.5) n_t=v(4)
      if(ieq.ge.5) p_t=v(5)
      if(ieq.ge.6) u_par=v(6)
      if(ieq.ge.8) n_u=v(7)
      if(ieq.ge.8) a_par=v(8)
      if(ieq.ge.9) n_im=v(9)
      if(ieq.eq.12) p_im_par=v(10)
      if(ieq.eq.12) p_im_per=v(11)
c     if(ieq.eq.12) u_im_par=v(12)
 
      if (ieq.ge.8) ph=((1.D0-dil)*apwt*n_i+aiwt*zimp*n_im-n_t-n_u)/f0
      if (ieq.lt.8) ph= ((1.D0-dil)*apwt*n_i-n_t)/f0
      if (ieq.le.3) ph= ((1.D0-dil)*apwt*n_i)/f0
      t_u=(betae/2.D0)*w_s/k_par0*rlte*a_par
 
      n_e=(1.D0-dil)*apwt*(n_i-(g0-g1)/taui*ph)
     >     +aiwt*zimp*(n_im-zimp*(g0i-g1i)/taui*ph)
 
c impurity trace convective limit
      if (aiwt.lt.-epsilon) then
       n_im=0.D0
       do j=1,8
        n_im=n_im+amat(9,j)*v(j)/(-xi*freq+gamma)
       enddo
      endif
 
c idelta=xi*yparam(1)+yparam(2)   for trapped electrons
 
      yparam(1)=aimag(-(n_t-reps*ph)/ph)
      yparam(2)=AREAL(-(n_t-reps*ph)/ph)
 
 
      chknu=n_u/(1.D0-reps)/ph
      chknt=n_t/reps/ph
      chknt2=n_t*(f0+reps)/(reps*((1.D0-dil)*apwt*n_i+aiwt*zimp*n_im))
      if (lprint.eq.2) write (6,*) 'chknu,chknt,chknt2:',
     >    chknu,chknt,chknt2
 
c non linear saturation rule
 
      gamma_r= 0.2D0*3.D0/2.D0*abs(w_d)*taui   !only scaling important
      if(iglf.eq.1) gamma_r= 0.2D0*3.D0/2.D0*abs(w_d0)*taui
c
      gamma_net=gamma-abs(alpha_star*gamma_star
     >         +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp

      gamma_net=max(gamma_net,xparam(8)*gamma)
      if( gamma_net.gt.0.D0)then
       ph_m=gamma_net**(1.D0-adamp)*gamma_r**adamp/(k_m*ky)
c set flag ngrow_k_gf: found at least one unstable mode for this k
       if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
      endif
 
      if( gamma_net.le.0.) ph_m=0.D0
c
      if(xparam(24).gt.0.D0) then
        if(gamma.gt.0.D0)then
          ph_m=dabs(gamma)**(1.D0-adamp)*gamma_r**adamp/(k_m*ky)/
     >    dsqrt(1.D0+(dabs(alpha_star*gamma_star+
     >    alpha_e*gamma_e+alpha_mode*gamma_mode)/
     >    (dabs(gamma)+.00001D0))**xparam(24))
          if(ipert_gf.eq.0)ngrow_k_gf(iky0)=1
        else
           ph_m=0.D0 
        endif
      endif
 
c 7.17.96
      if(xparam(22).gt.0) then
        if(gamma.gt.0.) ph_m=gamma**(1.D0-adamp-xparam(22))
     >   *gamma_r**adamp/(k_m*ky)
        if(gamma.le.0.) ph_m=0.D0
      endif
 
         phi_norm=0
         phi_normroot(iroot)=0.D0
 
       if( ph_m.gt.0.) then
 
       phi_norm=ph_m*ph_m/ABS((conjg(ph)*ph))
 
       phi_normroot(iroot)=phi_norm
 
c note only real part survives in diffusivities
c    ...units are c_s*rho_s**2/a
c magnetic futter component is too small to worry about
 
      d_hat    = phi_norm*AREAL(conjg(n_i)*(-xi*ky*ph))/rlni
     >+d_hat
 
      d_im_hat    = phi_norm*AREAL(conjg(n_im)*(-xi*ky*ph))/(rlnimp+
     >epsilon)
     >+d_im_hat
 
      chii_hat = phi_norm*3.D0/2.D0*
     >AREAL(conjg((1.D0/3.D0)*p_par+(2.D0/3.D0)*p_per)*(-xi*ky*ph))/rlti
     >+chii_hat
      chii_hat=chii_hat + aiwt/apwt*xparam(21)*phi_norm*3.D0/2.D0*
     >AREAL(conjg((1.D0/3.D0)*p_im_par+(2.D0/3.D0)*p_im_per)*
     >   (-xi*ky*ph))/rlti
 
      chie_hat = phi_norm*3.D0/2.D0*
     >AREAL(conjg(p_t+n_u+t_u)*(-xi*ky*ph))/rlte
     >+chie_hat
 
c electron to ion energy exchange in units n0*t0*c_s/a*(rho_a/a)**2
c note here we interpret QLT d/dt=-xi*freq dropping gamma part to
c avoid getting nonzero result for n_e->ph adiabatic limit
c ie <(d n_e/dt)*conjg(ph)>_time ave -> 0 adiabatic limit
c note  (-1) means exch_hat is electron to ion rate or
c ion heating rate, ie positive exch_hat cools electrons and heats ions
 
      exch_hat = phi_norm*(-1.D0)*
     >AREAL(conjg(-xi*freq*n_e)*ph)
     >+exch_hat
 
      eta_par_hat=(1.D0-xparam(14))*phi_norm*
     >AREAL(conjg(u_par)
     >*(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)
     >+xparam(14)*phi_norm*AREAL(conjg(
     >-xi*ky*gamma_p*(-gamma_p*alpha_n)*g1*ph/(-xi*freq+gamma))
     >*(-xi*ky*ph))/(gamma_p+epsilon)*(-gamma_p*alpha_n)
     >+eta_par_hat
 
      eta_per_hat = phi_norm*
     >AREAL(conjg(-ky*(ky*shat*rms_theta)*ph)*
     >   (ph+taui*((1.D0/3.D0)*p_par+(2.D0/3.D0)*p_per)))*
     >   (-gamma_p*alpha_n)
     >/(gamma_p+epsilon)
     >+eta_per_hat
 
       endif
       endif
      enddo
 777  continue
c
      if(ilh.eq.1) then
 
       do i=1,nmode
        yparam_k_gf(i,iky)=yparam(i)
       enddo
       xkyf_k_gf(iky)=kyf
 
       do iroot=1,4
        gamma_k_gf(iroot,iky)=gammaroot(iroot)
        freq_k_gf(iroot,iky)=freqroot(iroot)
        phi_norm_k_gf(iroot,iky)=phi_normroot(iroot)
       enddo
 
       diff_k_gf(iky)=d_hat
       diff_im_k_gf(iky)=d_im_hat
       chii_k_gf(iky)=chii_hat
       chie_k_gf(iky)=chie_hat
       exch_k_gf(iky)=exch_hat
 
 
ctest       exch_gf=-freq_gf(1)/xkyf_gf*diff_gf*rlni
 
       eta_par_k_gf(iky)=eta_par_hat
       eta_per_k_gf(iky)=eta_per_hat
 
c b_pol/b_phi=rmin/(rmaj*q)
       eta_phi_k_gf(iky)=eta_par_k_gf(iky)+
     >  rmin/(rmaj*q)*eta_per_k_gf(iky)
 
      endif
c
c computed high-k ETG electron transport at each ky
c Note: not added to ITG chi-e here ... done after ky loop
      chie_e_k_gf(iky)=0.D0
      if(ilh.eq.2) then
        chie_e_k_gf(iky)=xparam(10)*chii_hat*
     >                   taui_gf**(3.D0/2.D0)/
     >                   (1836.D0*amassgas_gf)**.5D0
      endif
c end ky loop
      enddo
 
c end big loop on ilh ... no longer used
c     enddo
c
c check to see if any unstable modes were found
c
      if(ipert_gf.eq.0)then
        do j=1,nmode
         if(ngrow_k_gf(j).ne.0)ngrow_k_gf(0)=1
        enddo
        if(ngrow_k_gf(0).eq.0)go to 888
      endif
c
c
c***********************************************************************
c initializations for summations over ky
c
      anorm_k=0.D0
      diff_gf=0.D0
      diff_im_gf=0.D0
      chii_gf=0.D0
      chie_gf=0.D0
      exch_gf=0.D0
      eta_par_gf=0.D0
      eta_per_gf=0.D0
      eta_phi_gf=0.D0
      chie_e_gf=0.D0
c
c Sum ITG and ETG transport
c over logarithmic ky grid (d ky=ky*d yk)
c
      do iky=1,ikymax_gf
       del_k=xkyf_k_gf(iky)
       anorm_k=anorm_k+del_k
       diff_gf=diff_gf+diff_k_gf(iky)*del_k
       diff_im_gf=diff_im_gf+diff_im_k_gf(iky)*del_k
       chii_gf=chii_gf+chii_k_gf(iky)*del_k
       chie_gf=chie_gf+chie_k_gf(iky)*del_k
       exch_gf=exch_gf+exch_k_gf(iky)*del_k
       eta_par_gf=eta_par_gf+eta_par_k_gf(iky)*del_k
       eta_per_gf=eta_per_gf+eta_per_k_gf(iky)*del_k
       eta_phi_gf=eta_phi_gf+eta_phi_k_gf(iky)*del_k
       chie_e_gf=chie_e_gf+chie_e_k_gf(iky)*del_k
      enddo
c
c Add ITG and ETG electron transport
c
      chie_gf=chie_gf + 1.D0*chie_e_gf
c
      diff_gf=diff_gf/anorm_k
      diff_im_gf=diff_im_gf/anorm_k
      chii_gf=chii_gf/anorm_k
      chie_gf=chie_gf/anorm_k
      exch_gf=exch_gf/anorm_k
      eta_par_gf=eta_par_gf/anorm_k
      eta_per_gf=eta_per_gf/anorm_k
      eta_phi_gf=eta_phi_gf/anorm_k
      chie_e_gf=chie_e_gf/anorm_k
c
c
c pick off maximum gamma
c 
      do iroot=1,4
       gamma_k_max=-1.D6
       do iky=1,ikymax_gf
        if(gamma_k_gf(iroot,iky).gt.gamma_k_max) then
         gamma_k_max=gamma_k_gf(iroot,iky)
         gamma_gf(iroot)=gamma_k_gf(iroot,iky)
         freq_gf(iroot)=freq_k_gf(iroot,iky)
         xky_gf(iroot)=xkyf_k_gf(iky)
        endif
       enddo
      enddo
c
c pick off 2nd maximum gamma
c
       gamma_k_max=-1.D6
       do iky=1,ikymax_gf
        if( (gamma_k_gf(1,iky).gt.gamma_k_max) .and.
     >      (gamma_k_gf(1,iky).lt.gamma_gf(1)) ) then
         gamma_k_max=gamma_k_gf(1,iky)
         gamma_gf(2)=gamma_k_gf(1,iky)
         freq_gf(2)=freq_k_gf(1,iky)
         xky_gf(2)=xkyf_k_gf(iky)
        endif
       enddo
c
c       write(6,*) gamma_gf(1), gamma_gf(2), xky_gf(1), xky_gf(2)
c
c print to file log
c      write(*,66)chii_gf,(gamma_gf(j),j=1,4)
 66    format(f14.9,4f14.9)
 67    format(2i2,f14.9)

 
      if(xparam(22).gt.0.) then
       phi_renorm=1.D0
       gamma_gross_net=gamma_gf(1)-abs(alpha_star*gamma_star
     >         +alpha_e*gamma_e+alpha_mode*gamma_mode)-kdamp
       if(gamma_gross_net.gt.0.)
     >  phi_renorm=gamma_gross_net**(xparam(22)*2.D0)
       if(gamma_gross_net.le.0.) phi_renorm=0.D0
 
      diff_gf=diff_gf*phi_renorm
      diff_im_gf=diff_im_gf*phi_renorm
      chii_gf=chii_gf*phi_renorm
      chie_gf=chie_gf*phi_renorm
      exch_gf=exch_gf*phi_renorm
      eta_par_gf=eta_par_gf*phi_renorm
      eta_per_gf=eta_per_gf*phi_renorm
      eta_phi_gf=eta_phi_gf*phi_renorm
      chie_e_gf=chie_e_gf*1.0D0
 
 
      endif
 
c put in cnorm_gf 12/22/96
 
 
      diff_gf=cnorm_gf*diff_gf
      diff_im_gf=cnorm_gf*diff_im_gf
      chii_gf=cnorm_gf*chii_gf
      chie_gf=cnorm_gf*chie_gf
      exch_gf=cnorm_gf*exch_gf
      eta_par_gf=cnorm_gf*eta_par_gf
      eta_per_gf=cnorm_gf*eta_per_gf
      eta_phi_gf=cnorm_gf*eta_phi_gf
      chie_e_gf=cnorm_gf*chie_e_gf
 
 
      if(lprint.gt.0) then
          write(1,*) 'gamma_gf=',  gamma_gf
          write(1,*) 'freq_gf=',  freq_gf
          write(1,*) 'ph_m=',  ph_m
          write(1,*) 'diff_gf=', diff_gf
          write(1,*) 'diff_im_gf=', diff_im_gf
          write(1,*) 'chii_gf=', chii_gf
          write(1,*) 'chie_gf=', chie_gf
          write(1,*) 'exch_gf=', exch_gf
      endif
 
 
      if (lprint.eq.98) then
        write(6,*) 'rlti,rlte,rlne,rlni,rlnimp: ',
     >    rlti,rlte,rlne,rlni,rlnimp
       write(6,*) 'chii,chie,diff,diff_im: ',
     >    chii_gf,chie_gf,diff_gf,diff_im_gf
       write(6,*) 'gamma_gf,freq_gf,ph_m: ',
     >    gamma_gf,freq_gf,ph_m
       write(6,*) 'jmax: ',
     >    jmax
        write(2,*) 'rlti,rlte,rlne,rlni,rlnimp: ',
     >    rlti,rlte,rlne,rlni,rlnimp
       write(2,*) 'chii,chie,diff,diff_im: ',
     >    chii_gf,chie_gf,diff_gf,diff_im_gf
       write(2,*) 'gamma_gf,freq_gf,ph_m: ',
     >    gamma_gf,freq_gf,ph_m
       write(2,*) 'jmax: ',
     >    jmax
 
        write (2,*) ' ar(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (ar(j1,j2),j2=1,neq)
        enddo
c
        write (2,*)
        write (2,*) ' ai(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (ai(j1,j2),j2=1,neq)
        enddo
c
        write (2,*) ' vr(j1,j2)  j2 ->'
        do j1=1,neq
          write (2,*) (vr(j1,j2),j2=1,neq)
        enddo
c
        write (2,*)
        write (2,*) ' vi(j1,j2)  j2->'
        do j1=1,neq
          write (2,*) (vi(j1,j2),j2=1,neq)
        enddo
 
      endif
 
 999  continue
      if (lprint.gt.0) close(1)

      return
c
c return for case with  no unstable modes
c
 888  continue  
      diff_gf=0.D0
      diff_im_gf=0.D0
      chii_gf=0.D0
      chie_gf=0.D0
      exch_gf=0.D0
      eta_par_gf=0.D0
      eta_per_gf=0.D0
      eta_phi_gf=0.D0
      chie_e_gf=0.D0
      do j1=1,4
        gamma_gf(j1)=0.D0
        freq_gf(j1)=0.D0
        xky_gf(j1)=0.D0
      enddo
      return

      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccgg
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c 
      subroutine cgg_glf(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)

      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)

c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.

c     on input

c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.

c        n  is the order of the matrix  a=(ar,ai).

c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.

c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.

c     on output

c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.

c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.

c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.

c        fv1, fv2, and  fv3  are temporary storage arrays.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50

   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)

      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s

c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by  cbal.

c        scale contains information determining the permutations
c          and scaling factors used by  cbal.

c        m is the number of eigenvectors to be back transformed.

c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.

c     on output

c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120

      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.000/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue

  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140

         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue

  140 continue

  200 return
      end
      subroutine cbal(nm,n,ar,ai,low,igh,scale)

      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv

c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.

c     on output

c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.

c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.

c        scale contains information determining the
c           permutations and scaling factors used.

c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.

c     note that 1 is returned for igh if igh is zero formally.

c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)

c     arithmetic is real throughout.

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      radix = 16.000

      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50

      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue

      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue

   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj

         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.000 .or. ai(j,i) .ne. 0.000) go to 120
  110    continue

         m = l
         iexc = 1
         go to 20
  120 continue

      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1

  140 do 170 j = k, l

         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.000 .or. ai(i,j) .ne. 0.000) go to 170
  150    continue

         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.000
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.

      do 270 i = k, l
         c = 0.000
         r = 0.000

         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.000 .or. r .eq. 0.000) go to 270
         g = r / radix
         f = 1.000
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.000 / f
         scale(i) = scale(i) * f
         noconv = .true.

         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue

         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue

  270 continue

      if (noconv) go to 190

  280 low = k
      igh = l
      return
      end
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci

c     complex division, (cr,ci) = (ar,ai)/(br,bi)

      double precision s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)

      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag,dlapy3gf

c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.

c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.

c     on output

c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.

c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.

c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.

c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  sqrt(a*a + b*b) .

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1

      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.000) go to 170
         norm = dlapy3gf(hr(i,i-1),hi(i,i-1))
crew inserted norm+1.d-100
         yr = hr(i,i-1) / (norm+1.d-100)
         yi = hi(i,i-1) / (norm+1.d-100)
         hr(i,i-1) = norm
         hi(i,i-1) = 0.000

         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue

         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue

  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue

      en = igh
      tr = 0.000
      ti = 0.000
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.000 .and. xi .eq. 0.000) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.000
      yi = (hi(enm1,enm1) - si) / 2.000
      call csroot(yr**2-yi**2+xr,2.000*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.000) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.000

  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue

      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1

      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.000
         norm = dlapy3gf(dlapy3gf(hr(i-1,i-1),hi(i-1,i-1)),sr)
crew inserted norm+1.d-100
         xr = hr(i-1,i-1) / (norm+1.d-100)
         wr(i-1) = xr
         xi = hi(i-1,i-1) / (norm+1.d-100)
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.000
         hi(i,i-1) = sr / (norm+1.d-100)

         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue

  500 continue

      si = hi(en,en)
      if (si .eq. 0.000) go to 540
      norm = dlapy3gf(hr(en,en),si)
crew inserted norm+1.d-100
      sr = hr(en,en) / (norm+1.d-100)
      si = si / (norm+1.d-100)
      hr(en,en) = norm
      hi(en,en) = 0.000
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)

         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.000
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue

  600 continue

      if (si .eq. 0.000) go to 240

      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue

      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)

      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag, dlapy3gf

c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.

c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.000 for these elements.

c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.

c     on output

c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.

c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.

c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.

c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.

c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  sqrt(a*a + b*b) .

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated october 1989.

c     ------------------------------------------------------------------

      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n

         do 100 i = 1, n
            zr(i,j) = 0.000
            zi(i,j) = 0.000
  100    continue
         zr(j,j) = 1.000
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.000 .and. orti(i) .eq. 0.000) go to 140
         if (hr(i,i-1) .eq. 0.000 .and. hi(i,i-1) .eq. 0.000) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1

         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue

         do 130 j = i, igh
            sr = 0.000
            si = 0.000
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
crew inserted norm+1.d-100
            sr = sr / (norm+1.d-100)
            si = si / (norm+1.d-100)

            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue

  130    continue

  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1

      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.000) go to 170
         norm = dlapy3gf(hr(i,i-1),hi(i,i-1))
crew     inserted norm+1.d-100
         yr = hr(i,i-1) / (norm+1.d-100)
         yi = hi(i,i-1) / (norm+1.d-100)
         hr(i,i-1) = norm
         hi(i,i-1) = 0.000

         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue

         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue

         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue

  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue

      en = igh
      tr = 0.000
      ti = 0.000
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.000 .and. xi .eq. 0.000) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.000
      yi = (hi(enm1,enm1) - si) / 2.000
      call csroot(yr**2-yi**2+xr,2.000*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.000) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.000

  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue

      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1

      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.000
         norm = dlapy3gf(dlapy3gf(hr(i-1,i-1),hi(i-1,i-1)),sr)
crew inserted norm+1.d-100
         xr = hr(i-1,i-1) / (norm+1.d-100)
         wr(i-1) = xr
         xi = hi(i-1,i-1) / (norm+1.d-100)
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.000
         hi(i,i-1) = sr / (norm+1.d-100)

         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue

  500 continue

      si = hi(en,en)
      if (si .eq. 0.000) go to 540
      norm = dlapy3gf(hr(en,en),si)
crew inserted norm+1.d-100
      sr = hr(en,en) / (norm+1.d-100)
      si = si / (norm+1.d-100)
      hr(en,en) = norm
      hi(en,en) = 0.000
      if (en .eq. n) go to 540
      ip1 = en + 1

      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)

         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.000
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue

         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue

  600 continue

      if (si .eq. 0.000) go to 240

      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue

      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue

      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.000

      do 720 i = 1, n

         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue

      if (n .eq. 1 .or. norm .eq. 0.000) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.000
         hi(en,en) = 0.000
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.000
            zzi = 0.000
            ip1 = i + 1

            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue

            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.000 .or. yi .ne. 0.000) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.000) go to 780
            tst1 = tr
            tst2 = tst1 + 1.000/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue

  780    continue

  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840

         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue

  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)

         do 880 i = low, igh
            zzr = 0.000
            zzi = 0.000

            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue

            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue

      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)

      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag,dlapy3gf

c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.

c     on input

c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.

c        n is the order of the matrix.

c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.

c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.

c     on output

c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.

c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.

c     calls pythag for  sqrt(a*a + b*b) .

c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory

c     this version dated august 1983.

c     ------------------------------------------------------------------

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200

      do 180 m = kp1, la
         h = 0.000
         ortr(m) = 0.000
         orti(m) = 0.000
         scale = 0.000
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))

         if (scale .eq. 0.000) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue

         g = sqrt(h)
         f = dlapy3gf(ortr(m),orti(m))
         if (f .eq. 0.000) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.000 + g) * ortr(m)
         orti(m) = (1.000 + g) * orti(m)
         go to 105

  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.000
            fi = 0.000
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue

            fr = fr / h
            fi = fi / h

            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue

  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.000
            fi = 0.000
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue

            fr = fr / h
            fi = fi / h

            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue

  160    continue

         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue

  200 return
      end
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi

c     (yr,yi) = complex sqrt(xr,xi)
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)

      double precision s,tr,ti,pythag,dlapy3gf
      tr = xr
      ti = xi
      s = sqrt(0.5d0*(dlapy3gf(tr,ti) + abs(tr)))
      if (tr .ge. 0.000) yr = s
      if (ti .lt. 0.000) s = -s
      if (tr .le. 0.000) yi = s
      if (tr .lt. 0.000) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.000) yi = 0.5d0*(ti/yr)
      return
      end
      double precision function pythag(a,b)
      double precision a,b

c     finds sqrt(a**2+b**2) without overflow or destructive underflow

      double precision p,r,s,t,u
crew changed dmax1 to max
      p = max(abs(a),abs(b))
      if (p .eq. 0.000) go to 20
crew changed dmin1 to min
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.000 + r
c        write(*,*) 't = ',t
         if (abs(t-4.000) .lt. 1.e-5) go to 20
         s = r / t
         u = 1.000 + 2.000*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
      DOUBLE PRECISION FUNCTION DLAPY3GF( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
*     ..
*
*  Purpose
*  =======
*
*  DLAPY3GF returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*  Z       (input) DOUBLE PRECISION
*          X, Y and Z specify the values x, y and z.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      Z = 0
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
         DLAPY3GF = ZERO
      ELSE
         DLAPY3GF = W*SQRT( ( XABS / W )**2+( YABS / W )**2+
     $            ( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY3GF
c
      end
