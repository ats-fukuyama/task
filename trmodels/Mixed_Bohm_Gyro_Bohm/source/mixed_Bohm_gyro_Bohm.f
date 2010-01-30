!| %
!| %  This is a LaTeX ASCII file.  To typeset document type:
!| % latex theory
!| %  To extract the Fortran source code, type:
!| % python s2tex.py mixed_model.tex
!| %
!| % The following lines control how the LaTeX document is typeset
!| 
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| \begin{document}
!| 
!| \begin{center}
!| \Large {\bf Mixed Bohm/gyro-Bohm Transport Model} \\
!| \vspace{1pc} \normalsize
!| Glenn Bateman, Alexei Pankin, and Arnold Kritz \\
!| Lehigh University Physics Department \\
!| 16 Memorial Drive East, Bethlehem PA 18015 \\
!| bateman@fusion.physics.lehigh.edu \\
!| pankin@fusion.physics.lehigh.edu \\
!| kritz@fusion.physics.lehigh.edu
!| \end{center}
!| 
!| The \verb'mixed_Bohm_gyro_Bohm' module is used to compute
!| anomalous plasma transport coefficients using the Mixed-Bohm/gyro-Bohm
!| transport model.
!| The Mixed Bohm/gyro-Bohm anomalous transport model (also called the JET
!| or JETTO transport model) contains Bohm and gyro-Bohm contributions.
!| The Bohm contribution, which is linear in the gyro-radius, is a
!| non-local transport model, in which the transport throughout the plasma
!| depends on a finite difference approximation to the electron temperature
!| gradient at the edge of the plasma.  The gyro-Bohm contribution is a
!| local transport model, which is proportional to the square of the
!| gyro-radius.  The version of the Mixed Bohm/gyro-Bohm model in this
!| module is described in detail in Nuclear Fusion 38 (1998) 1013.
!| The magnetic and flow shear stabilization are described in Plasma
!| Physics and Controlled Fusion 44 (2002) A495.
!| 
!| There are two subroutines in this module:
!| 
!| The \verb'mixed_model' subroutine computes arrays of effective
!| diffusivities using arrays representing the profiles from the
!| center to the edge of the plasma.
!| 
!| The \verb'mixed' subroutine computes scalar effective diffusivities
!| given the plasma properties at a given point, and also given the
!| finite difference approximation to the electron temperature gradient
!| at the edge of the plasma.
!| 
!| The version of the Mixed-Bohm/gyro-Bohm
!| model that is used in this subroutine is described
!| in Nuclear Fusion {\bf 38} (1998) 1013 by
!| M.~Erba, T.~Aniel, V.~Basiuk, A.~Becoulet, and X.~Litaudon.
!| Both the electron and ion thermal diffusivities consist
!| of two terms.  One term has Bohm scaling
!| \begin{equation}
!| \chi^{\rm Bohm} \equiv \rho_s c_s  q^2
!| \frac{ a ( d p_e / d r ) }{ p_e } \Delta_{T_e},
!| \label{eq-Bohm}
!| \end{equation}
!| while the other term has gyro-Bohm scaling
!| \begin{equation}
!| \chi^{\rm gyro-Bohm} \equiv \frac{ \rho_s^2 c_s }{ a }
!|   \frac{ a ( d T_e / d r ) }{ T_e }.
!| \label{eq-gyro-Bohm}
!| \end{equation}
!| The notation is described in Table \ref{table-notation}.
!| In the Bohm diffusivity expression, $ \Delta_{T_e} $ is a finite
!| difference approximation to the normalized temperature electron
!| temperature difference at the plasma edge
!| \begin{equation}
!| \Delta_{T_e} \equiv \frac{T_e(r/a=0.8) - T_e(r/a=1)}{T_e(r/a=1)}.
!| \label{eq-Delta}
!| \end{equation}
!| The resulting anomalous
!| ion and electron thermal diffusivities are constructed from the
!| sum of these Bohm and gyro-Bohm terms, with empirically determined
!| coefficients \cite{erb98}
!| \begin{equation}
!| \chi_i^{\rm JET} = 1.6 \times 10^{-4} \chi^{\rm Bohm} +
!|  1.75 \times 10^{-2} \chi^{\rm gyro-Bohm}
!| \label{eq-chi-i}
!| \end{equation}
!| \begin{equation}
!| \chi_e^{\rm JET} = 8 \times 10^{-5} \chi^{\rm Bohm} +
!|  3.5 \times 10^{-2} \chi^{\rm gyro-Bohm}
!| \label{eq-chi-e}
!| \end{equation}
!| and the hydrogenic charged particle diffusivity is given by
!| \begin{equation}
!| D^{\rm JET} \propto  \frac{ \chi_i \chi_e }{ \chi_i + \chi_e }
!| \label{eq-D}
!| \end{equation}
!| 
!| It should be noted that the Mixed-Bohm/gyro-Bohm model is not a
!| local transport model --- it cannot be evaluated with just the
!| local plasma parameters.
!| The Bohm contribution to the core transport $ \chi^{\rm Bohm} $
!| is proportional to a finite difference approximation to the edge
!| electron temperature gradient through the term $ \Delta_{T_e} $.
!| 
!| 
!| \begin{table}
!| \caption{Notation}
!| 
!| % \renewcommand{\arraystretch}{0.7}
!| 
!| \label{table-notation}
!| \begin{center}
!| \begin{tabular}{lll}
!| \hline \hline
!| Variable \hspace{3pt} & Units \hspace{3pt} & Meaning \\
!| \hline % \\[-2mm]
!| $a$& m & minor radius (half-width) of plasma \\
!| $B_T$ & Tesla & vacuum toroidal magnetic field at \\ &
!|    & major radius $R$ along flux surface \\
!| $c_s$ & m/s & $ [k_b T_e / m_i]^{1/2} $ speed of sound \\
!| $D$ & m$^2$/s & effective charged particle diffusivity \\
!|  & & charged particle flux divided by density gradient \\
!| $e$ & C & electron charge \\
!| $I_p$ & MA & toroidal plasma current \\
!| $k_b$ &  & conversion from keV to Joules \\
!| $m_i$ & kg  & average ion mass \\
!| $n_e$ & m$^{-3}$ & electron density \\
!| $q$  & & magnetic q-value \\
!| $r$& m & minor radius (half-width) of each flux surface \\
!| $R$ & m & major radius to geometric \\
!|  & & center of flux surface \\
!| $T_e$ & keV & electron temperature \\
!| $T_i$ & keV & ion temperature \\
!| $Z_{\rm eff}$ & & $\sum_s n_s Z_s^2 / n_e $ summed over each species \\
!| $\chi$  & m$^2$/s & effective thermal diffusivity \\
!|  & & heat flux divided by density time temperature gradient \\
!| $\delta$ &  & plasma triangularity \\
!| $\kappa$ &  & plasma elongation \\
!| $\nu_{*}$ & & collision frequency divided by bounce frequency \\
!| $\rho_s$ & m &  gyroradius [$ c_s m_i / ( e B_T ) $] \\
!| $\rho_*$ & & normalized gyroradius ($ \rho_s / a $) \\
!| \hline \hline
!| \end{tabular}
!| \end{center}
!| \end{table}
!| 
c@mixed_model.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
      module mixed_Bohm_gyro_Bohm
      private
      public mixed_model
c
c..physical constants
c
      real, parameter  :: zckb = 1.60210e-16 ! energy conversion factor [Joule/keV]
      real, parameter  :: zcme = 9.1091e-31  ! electron mass [kg]
      real, parameter  :: zcmp = 1.67252e-27 ! proton mass [kg]
      real, parameter  :: zce  = 1.60210e-19 ! electron charge [Coulomb]
c
      contains
c
      subroutine mixed_model (
     &   rminor,  rmajor,  tekev,   tikev,   q
     & , btor,    aimass,  charge,  wexbs
     & , grdte,   grdne,   shear
     & , t_e_kev_edge, rminor_edge, npoints
     & , chi_i_mix,  themix,   thdmix
     & , thigb,   thegb,    thibohm, thebohm
     & , ierr
     & , lflowshear)
c
c
c    All the following 1-D arrays are assumed to be defined on flux
c    surfaces called zone boundaries where the transport fluxes are
c    to be computed.  The number of flux surfaces is given by npoints
c    (see below).  For example, if you want to compute the transport
c    on only one flux surface, set npoints = 1.
c
c  Input arrays:
c  -------------
c
c  rminor(jz)   = minor radius (half-width) of zone boundary [m]
c  rmajor(jz)   = major radius to geometric center of zone bndry [m]
c
c  tekev(jz)    = T_e (electron temperature) [keV]
c  tikev(jz)    = T_i (temperature of thermal ions) [keV]
c  q(jz)        = magnetic q-value
c  btor(jz)     = ( R B_tor ) / rmajor(jz)  [tesla]
c
c  aimass(jz)   = mean atomic mass of main thermal ions [AMU]
c               = ( sum_i n_i M_i ) / ( sum_i n_i ) where
c                 sum_i = sum over all ions, each with mass M_i
c
c  charge(jz)   = charge number of the main thermal ions
c                 = 1.0 for hydrogenic ions
c
c  wexbs(jz)    = ExB shearing rate in [rad/s]
c
c    All of the following normalized gradients are at zone boundaries.
c    r = half-width, R = major radius to center of flux surface
c
c  grdte(jz) = -R ( d T_e / d r ) / T_e
c  grdne(jz) = -R ( d n_e / d r ) / n_e
c  shear(jz) =  r ( d q   / d r ) / q    magnetic shear
c
c    The following variables are scalars:
c
c  t_e_kev_edge = electron temperature at the edge of the plasma [keV]
c  rminor_edge  = minor radius at the edge of the plasma [m]
c
c  npoints = number of values of jz in all of the above arrays [integer]
c
c  Control variable (input):
c  -------------------------
c
c  lflowshear = 0 for no magnetic and flow shear stabilization
c             = 1 to use magnetic and flow shear stabilzation by
c                 [T.J. Tala et al Plasma Phys. Controlled Fusion 44
c                 (2002) A495]
c
c  Output:
c  -------
c
c    The following effective diffusivities are given in MKS units m^2/sec
c
c  chi_i_mix(jz) = total ion thermal diffusivity from the MIXED model
c  themix(jz) = total electron thermal diffusivity from the MIXED model
c  thdmix(jz) = total hydrogenic ion diffusivity from the MIXED model
c
c    The following contributions to the effective diffusivities are
c  for diagnostic purposes:
c
c  thigb(jz) = gyro-Bohm contribution to the ion thermal diffusivity
c  thegb(jz) = gyro-Bohm contribution to the electron thermal diffusivity
c
c  thibohm(jz) = Bohm contribution to the ion thermal diffusivity
c  thebohm(jz) = Bohm contribution to the electron thermal diffusivity
c
c  ierr    = returning with value .ne. 0 indicates error
c
c***********************************************************************
c
c-----------------------------------------------------------------------
c
c  Compile this routine and routines that it calls with a compiler
c  option, such as -r8, to convert real to double precision when used on
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: mixed_model calls the following routines
c
c  mixed             - Computes diffusivity from the MIXED model
c
c-----------------------------------------------------------------------

      implicit none
c
c-----------------------------------------------------------------------
c..input variables
c
      integer, intent(in)  :: npoints ! number of radial points
c
      real, intent(in) ::
     &   rminor(*),  rmajor(*)
     & , tekev(*),   tikev(*),    q(*),       btor(*)
     & , wexbs(*),   aimass(*),   charge(*)
c
      real, intent(in) ::                    ! gradients
     &   grdne(*)
     & , grdte(*), shear(*)
c
      real, intent(in) ::  t_e_kev_edge, rminor_edge
c
c  optional input switch for flow-shear stabilization option
c
      integer, intent(in), optional :: lflowshear
c
c-----------------------------------------------------------------------
c..output variables
c
      real, intent(out) ::
     &   thigb(*),   thegb(*),    thibohm(*), thebohm(*)
c
      real, intent(out) ::
     &   chi_i_mix(*),   themix(*),  thdmix(*)
c
      integer, intent(out) :: ierr
c
c-----------------------------------------------------------------------
c..local variables
c
      integer  :: jz, j1, j2, jm
c
c  npoints1 = value of jz corresponding to 80 0f the normalized radius
c
      integer  :: npoints1, llflow

      integer  :: lswitch5
c
      real     :: zte_p8, zte_edge, zi
c
c..local variables connected to the mixed module
c

      real :: zq,    zsound,     zgyrfi,    zrhos, zra
     & , zchii, zchie, zdhyd
     & , zrmaj, zaimass, zcharge,   zbtor, zrmin
     & , zgte,  zti,  zte, zgne
     & , zrlpe, zshear, zgradte
     & , zchbe, zchbi, zchgbe, zchgbi
c
c.. variables for exb model
c
      real  zwexb
c
c  zwexb    = local copy of ExB shearing rate
c
c..initialize arrays
c
c
      thigb(1:npoints)  = 0.
      thegb(1:npoints)  = 0.
      thibohm(1:npoints)= 0.
      thebohm(1:npoints)= 0.
      chi_i_mix(1:npoints) = 0.
      themix(1:npoints) = 0.
      thdmix(1:npoints) = 0.
c
c..initialize switches
c
      if (present(lflowshear)) then ! flow shear correction
        llflow = lflowshear
      else
        llflow = 0
      endif
c
c-----------------------------------------------------------------------
c
c..physical constants
c
c     define the jz value at 0.8 0f normalized radius
c
      npoints1 = int(npoints*0.8)
c
c.. start the main do-loop over the radial index "jz"..........
c
c
      do 300 jz = 1, npoints
c
c  compute scalar quantities necessary for mixed module
c
        zshear   = shear(jz)
        zte      = tekev(jz)
        zti      = tikev(jz)
        zaimass  = aimass(jz)
        zcharge  = charge(jz)
        zrmin    = rminor(jz)
        zrmaj    = rmajor(jz)
        zgte     = grdte(jz)
        zgne     = grdne(jz)
        zgradte   = abs( ( zgte / zrmaj ) * zte)
        zbtor    = btor(jz)
        zq       = q(jz)
        zra      = rminor(jz) / rminor(npoints)
        zrlpe    = (zgte + zgne)
        zte_p8   = tekev(npoints1)
        zte_edge = tekev(npoints)
c
        zwexb = wexbs(jz)
c
c
      call mixed(
     &  zbtor,          zgradte,        zq,
     &  zra,            zrlpe,
     &  zrmaj,          zshear,         zte,            zte_p8,
     &  zte_edge,       zti,            zwexb,          zaimass,
     &  zcharge,        llflow,
     &  zchie,          zchii,          zdhyd,
     &  zchbe,          zchbi,          zchgbe,         zchgbi,
     &  ierr)
c
c
c  If ierr not equal to 0 an error has occured
c
        if (ierr .ne. 0) return
c
c  compute effective diffusivites for diagnostic purposes only
c
         thdmix(jz)  = zdhyd
         themix(jz)  = zchie
         chi_i_mix(jz)  = zchii
c
c  put value of gbohm and bohm terms into kb and rb terms for
c  diagnostic purposes
c
         thebohm(jz) = zchbe
         thibohm(jz) = zchbi
c
         thegb(jz)   = zchgbe
         thigb(jz)   = zchgbi
c
c..end of mixed model
c
c
 300  continue
c
c
c   end of the main do-loop over the radial index, "jz"----------
c
      return
      end subroutine mixed_model
c--------1---------2---------3---------4---------5---------6---------7-c
c
!| \newpage
!| \begin{center}
!| {\LARGE Subroutine for Computing Particle and Energy Fluxes\\ \vskip8pt
!| Using the JET Mixed Bohm/gyro-Bohm\\ \vskip8pt
!| Transport Model
!| }\vskip1.0cm
!| Version 1.2: 22 August 2003 \\
!| Implemented by M. Erba, G. Bateman, A. H. Kritz, T. Onjun, and A. Pankin\\
!|  Lehigh University
!| \end{center}
!| For questions about this routine, please contact: \\
!| Arnold Kritz, Lehigh: {\tt kritz@plasma.physics.lehigh.edu}\\
!| Glenn Bateman, Lehigh: {\tt glenn@plasma.physics.lehigh.edu}\\
!| Thawatchai Onjun, Lehigh: {\tt onjun@fusion.physics.lehigh.edu}\\
!| Alexei Pankin, Lehigh: {\tt pankin@fusion.physics.lehigh.edu}\vskip8pt
!| 
!| 
c@mixed.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine mixed(
     &  btor,           gradte,         q_safety,
     &  ra,             rlpe,
     &  rmaj,           shear,          tekev,          te_p8,
     &  te_edge,        tikev,          wexb,           aimass,
     &  zi,             lflow,
     &  chi_e,          chi_i,          d_hyd,
     &  chi_e_bohm,     chi_i_bohm,
     &  chi_e_gyro_bohm, chi_i_gyro_bohm,
     &  ierr)
c
c Revision History
c ----------------
c      date          description
c
c   25-Aug-2003      New version of flow shear sabilization correction
c                    Converted to Fortran-90 module
c   27-Jul-2000      Major rewrite by Bateman
c   12-Jun-2000      gyro-Bohm term restored
c   06-May-1999      Revised as module by Thawatchai Onjun
c   24-Feb-1999      First Version by Matteo Erba
c
c Mixed Transport Model (by M.Erba, V.V.Parail, A.Taroni)
c
c Inputs:
c    btor:      Toroidal magnetic field strength in Tesla
c                 at geometric center along magnetic flux surface
c    gradte:    Local electron temperature gradient [keV/m]
c    q_safety:  Local value of q, the safety factor
c    ra:        Normalized minor radius r/a
c    rlpe:      a/L_pe, where a is the minor radius of the plasma,
c                    and L_pe is the local electron pressure scale length
c    rmaj:      Major radius [m]
c    shear:     Local magnetic Shear
c    tekev:     Local electron temperature [keV]
c    tikev:     Ion temperature [keV]
c    te_p8:     Electron temperature at r/a = 0.8
c    te_edge:   Electron temperature at the edge [r/a = 1]
c    wexb:      EXB Rotation
c                    "Effects of {ExB} velocity shear and magnetic shear
c                    on turbulence and transport in magnetic confinement
c                    devices", Phys. of Plasmas, 4, 1499 (1997).
c    zi:        Ion charge
c
c Outputs:
c    chi_e:      Total electron thermal diffusivity [m^2/sec]
c    chi_i:      Total ion thermal diffusivity [m^2/sec]
c    d_hyd:      Hydrogenic ion particle diffusivity [m^2/sec]
c
c    chi_e_bohm:      Bohm contribution to electron thermal diffusivity
c                         [m^2/sec]
c    chi_i_bohm:      Bohm contribution to ion thermal diffusivity [m^2/sec]
c    chi_e_gyro_bohm: gyro-Bohm contribution to electron thermal diffusivity
c                 [m^2/sec]
c    chi_i_gyro_bohm:    gyro-Bohm contribution to ion thermal diffusivity
c                 [m^2/sec]
c    ierr:      Status code returned; 0 = OK, .ne.0 indicates error
c
c Coefficients set internally:
c
c    alfa_be:   Bohm contribution to electron thermal diffusivity
c    alfa_bi:   Bohm contribution to ion thermal diffusivity
c    alfa_gbe:  gyro-Bohm contribution to electron thermal diffusivity
c    alfa_gbi:  gyro-Bohm contribution to ion thermal diffusivity
c
c    coef1:     First coefficient for empirical hydrogen diffusivity
c    coef2:     Second coefficient for empirical hydrogen diffusivity
c
c Other internal variables:
c
c    gamma:     The characteristic growthrate for the ITG type of
c                    electrostatic turbulence
c    func:      Function for EXB and magnetic shear stabilization
c
      IMPLICIT NONE
c
c Declare variables
c
c..Input variables
c
      real, intent(in) ::
     &  btor,           ra,             rmaj,           rlpe,
     &  shear,          tekev,          te_p8,          te_edge,
     &  tikev,          gradte,         aimass,
     &  wexb,           zi,             q_safety
c
      integer, intent(in)  :: lflow
c
c..Output variables
c
      real, intent(out)    ::
     &  chi_e,               chi_i,               d_hyd,
     &  chi_e_gyro_bohm,     chi_i_gyro_bohm,
     &  chi_e_bohm,          chi_i_bohm
c
c
      integer, intent(out) :: ierr
c
c..Local variables
c
      REAL
     &  alfa_be,        alfa_bi,        alfa_gbe,       alfa_gbi,
     &  alfa_d,         chi0,           coef1,          coef2,
     &  delta_edge,
     &  em_i,           func,           gamma,
     &  omega_ce,       omega_ci,       rho,
     &  v_sound,        vte_sq,         vti,
     &  zepsilon
c
c
c..initialize diffusivities
c
      chi_e = 0.0
      chi_i = 0.0
      d_hyd = 0.0
c
      chi_e_bohm = 0.0
      chi_i_bohm = 0.0
      chi_e_gyro_bohm = 0.0
      chi_i_gyro_bohm = 0.0
c
c check input for validity
c
      zepsilon = 1.e-10
c
      ierr = 0
      if ( tekev .lt. zepsilon ) then
         ierr=1
         return
      elseif ( tikev .lt. zepsilon ) then
         ierr=2
         return
      elseif ( te_p8 .lt. zepsilon ) then
         ierr=3
         return
      elseif ( te_edge .lt. zepsilon ) then
         ierr=4
         return
      elseif ( rmaj  .lt. zepsilon ) then
         ierr=5
         return
      endif
!| 
!| \section{The Mixed Bohm/gyro-Bohm model}
!| 
!| The Mixed Bohm/gyro-Bohm transport model derives from an
!| originally purely Bohm-like model for electron transport
!| developed for the JET Tokamak\cite{tar94}. This preliminary model has
!| subsequently been extended to describe ion transport\cite{erb95},
!| and a gyro-Bohm term has been added in order to simulate data from
!| different machines\cite{erb98}.
!| 
!| \subsection{Bohm term}
!| 
!| The mixed model is derived using the dimensional analysis approach,
!| whereby the diffusivity in a Tokamak plasma can be written as:\hfill
!| \[ \chi = \chi_0 F(x_1, x_2, x_3, ...)\]
!| where $\chi_0$ is some basic transport coefficient and F is a function
!| of the plasma dimensionless parameters $(x_1,\ x_2,\ x_3,\ ...)$. We choose
!| for $\chi_0$ the Bohm diffusivity:\hfill
!| \[ \chi_0 = \frac{cT_e}{eB}\]
!| 
!| The expression of the dimensionless function F is chosen according to
!| the following criteria:\
!| \begin{itemize}
!| \item{The diffusivity must be bowl-shaped, increasing towards the plasma
!| boundary}
!| \item{The functional dependencies of F must be in agreement with
!| scaling relationships of the global confinement time, reflecting
!| trends such as power degradation and linear dependence on plasma
!| current}
!| \item{The diffusivity must provide the right degree of resilience
!| of the temperature profile}
!| \end{itemize}
!| 
!| It easily shown that a very simple expression of F that satisfies
!| the above requirements is:\hfill
!| \[ F = q^2/|L_{pe}^*|\]
!| 
!| where q is the safety factor and $L_{pe}^*=p_e(dp_e/dr)^{-1}/a$, being a the
!| plasma minor radius. The resulting expression of the diffusivity
!| can be written as:\hfill
!| \[ \chi \propto |v_d| \Delta G\]
!| where $v_d$ is the plasma diamagnetic velocity, $\Delta=a$ and $G=q^2$,
!| so that it is clear that this model represents transport due to
!| long-wavelength turbulence.\\
!| The evidence coming up from the simulation of non-stationary
!| JET experiments \cite{erb97}(such as ELMs, cold pulses, sawteeth, {\sl etc}.)
!| suggested that the above Bohm term should depend non-locally
!| on the plasma edge conditions through the temperature
!| gradient averaged over a region near the edge:\hfill
!| 
!| \[ <L_{T_e}^*>_{\Delta V}^{-1} = \frac{T_e(x=0.8) - T_e(x=1)}{T_e(x=1)}\]
!| where x is the normalized toroidal flux coordinate. The final
!| expression of the Bohm-like model is:\hfill
!| 
!| \[ \chi_{e,i}^B = \alpha_{Be,i} \frac{cT_e}{eB} L_{pe}^{*-1} q^2 <L_{T_e}^*>_{\Delta V}^{-1}\]
!| where $\alpha_{e,i}^B$ is a parameter to be determined empirically,
!| both for ions and electrons.\hfill
!| 
!| \subsection{gyro-Bohm term}
!| 
!| The Bohm-like expression so derived proved to be very successful
!| in simulating JET discharges, but failed badly in smaller Tokamaks
!| such as START\cite{roa96}.\\
!| For this reason a simple gyro-Bohm-like term, also based on
!| dimensional analysis, was added:
!| 
!| \[ \chi_{e,i}^{gB} = \alpha_{e,i}^gB \frac{cT_e}{eB} L_{Te}^{*-1} \rho^*\]
!| where $\rho^*$ is the normalized larmor radius:
!| 
!| \[ \rho^* =  \frac {M^{1/2}cT_e^{1/2}}{aZ_ieB_t}\]
!| 
!| This expression is what can be expected from small scale
!| drift-wave turbulence. It is important to note that in large
!| Tokamaks such as JET and TFTR the gyro-Bohm term is negligible,
!| while in smaller machines, with larger values of $\rho^*$, the
!| gyro-Bohm term can play a role especially near the plasma centre.\hfill
!| 
!| \subsection{Final Model}
!| The resulting expressions of the diffusivities are:
!| 
!| \[ \chi_{e,i}=\chi_{Be,i}+\chi_{gBe,i}\]
!| where the Bohm and gyro-Bohm terms are defined above and the adopted values
!| of the empirical parameters are:\hfill
!| 
!| \[ \alpha_{Be} = 8\times10^{-5} ,\alpha_{Bi} = 2\times\alpha_{Be}\]
!| \[ \alpha_{gBe} = 3.5\times10^{-2} , \alpha_{gBi} = \alpha_{gBe}/2\]\vskip8pt
!| 
!| 
!| The following definitions are used:\vskip8pt
!| 
!| \begin{tabular}{lll}
!| el\_mass &electron mass    &$m_e$ \\
!| c       &velocity of light &$c$ \\
!| e       &electron charge  &$e$\\
!| vte\_sq &eletron thermal velocity squared &$v_{\rm te}^2$\\
!| omega\_ce &electron cyclotron frequency &$\omega_{\rm ce}$\\
!| chi0 &Bohm diffusitivity & $\chi_0$\\
!| em\_i & ion atomic mass [kg] & $M_{\rm i}$\\
!| %v_sound & ion
!| \end{tabular}\vskip8pt
!| 
c *
c * Definition of the mixed model
c *
c
c Coefficients
c
      alfa_be  =  8.00000000000000E-05
      alfa_bi  =  1.60000000000000E-04
      alfa_gbe =  3.50000000000000E-02
      alfa_gbi =  1.75000000000000E-02
c
      coef1    =  1.00000000000000E+00
      coef2    =  3.00000000000000E-01
c
c Calculate chi0
c
      vte_sq    = tekev * zckb / zcme
      omega_ce  = zce * btor / zcme
      chi0      = vte_sq / omega_ce
c
c Calculate chibohm
c
      delta_edge = abs((te_p8 - te_edge) / te_edge)
      chi_e_bohm
     &   = abs(alfa_be * chi0 * rlpe * (q_safety**2) * delta_edge)
      chi_i_bohm
     &   = abs(alfa_bi * chi0 * rlpe * (q_safety**2) * delta_edge)
c
c Calculate chi_gyrobohm
c
      em_i      = aimass * zcmp
      v_sound   = sqrt (tekev * zckb /em_i) !                   [m/sec]
      omega_ci  = zi * zce * btor / em_i
      rho       = v_sound / omega_ci
      chi_e_gyro_bohm = abs(alfa_gbe * chi0 * gradte * rho / tekev)
      chi_i_gyro_bohm = abs(alfa_gbi * chi0 * gradte * rho / tekev)
c
c Calculate function for EXB and magnetic shear stabilization
c   use lflow = 1 for the stabilization term described by
c   T.J. Tala et al Plasma Phys. Controlled Fusion 44 (2002) A495
c
c
      vti       = sqrt (2.0 * tikev * zckb / em_i)
      gamma     = vti / (q_safety * rmaj)
c
      func = 1.0
      if ( lflow .eq. 1 ) then
        func = -0.14 + shear - 1.47 * abs(wexb*rmaj/vti)
      elseif ( lflow .eq. 2 ) then
        func = 0.1 + shear - abs (wexb / gamma)
      endif
c
      if ( func .lt. 0.0 ) then
          chi_e_bohm = 0.0
          chi_i_bohm = 0.0
      endif
c
c Now determine the actual electron and ion thermal and particle
c diffusivities.
c
      chi_e      = chi_e_bohm + chi_e_gyro_bohm
      chi_i      = chi_i_bohm + chi_i_gyro_bohm
      alfa_d    = coef1 + (coef2 - coef1) * ra
      if (abs (chi_e + chi_i) .lt. 1e-10) then
        d_hyd = 0.0
      else
        d_hyd  = abs(alfa_d * chi_e * chi_i / (chi_e + chi_i))
      endif
c
c The impurity diffusivity is not included in the mixed
c model described in Ref. [1].
c
        return
        end subroutine mixed
        end module mixed_Bohm_gyro_Bohm
!| 
!| %**********************************************************************c
!| 
!| \begin{thebibliography}{99}
!| \bibitem{tar94}
!| A. Taroni, M. Erba, E. Springmann and Tibone F.,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 36} (1994) 1629.
!| \bibitem{erb95}
!| M. Erba, V. Parail, E. Springmann and A. Taroni,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 37} (1995) 1249.
!| \bibitem{erb98}
!| M. Erba, et al.,
!| {\em Nuclear Fusion,} {\bf 38} (1998) 1013.
!| \bibitem{erb97}
!| M. Erba, et al.,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 39} (1997) 261.
!| \bibitem{roa96}
!| C.M. Roach,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 38} (1996) 2187.
!| \end{thebibliography}
!| 
!| %**********************************************************************c
!| \end{document}
