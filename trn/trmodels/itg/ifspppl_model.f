c
c     The 1995 and 1994 IFS/PPPL transport models
c     as implemented by Aaron J Redd, February 1998
c
c     For problems with this routine, please contact:
c        Aaron Redd, Lehigh U:  aredd@plasma.physics.lehigh.edu
c        Arnold Kritz, Lehigh U:  kritz@plasma.physics.lehigh.edu
c     or Glenn Bateman, Lehigh U:  bateman@plasma.physics.lehigh.edu
c
c     Last revision:  July 7, 1998
c
      subroutine IFSPPPL( znine, zncne, znbne, zrlt, zrln,
     &                    zq, zshat, zeps, zkappa, omegaexb,
     &                    ne19, tekev, tikev, rmajor, btesla,
     &                    switches, grhoi, gvti, gnu, gtau,
     &                    chii, chie,
     &                    zchiicyc, zchii1, zchii2, zchie1, zchie2,
     &                    zrlt1, zrlt2, ierr )
c
      IMPLICIT NONE
c
c
c     switches(j) settings:
c
c        switches(1) = 0  No output
c                    = 1  Print page of inputs, computational steps and output
c        switches(2) = 0  Use first set of dimensioned input
c                         (ne19, tekev, tikev, btesla)
c                    = 1  Use second set of dimensioned/dimensionless input
c                         (grhoi, gvti, gnu, gtau)
c                         Note that, with either approach, rmajor is required.
c        switches(3) = 0  Use 1995 IFS/PPPL model
c                    = 1  Use 1994 IFS/PPPL model (IAEA paper)
c        switches(4) = 0  Use the value of gnu as given (if applicable)
c                    = 1  Multiply gnu by 0.84; this allows direct comparison
c                         with the fortran code provided by B. Dorland.
c                         There is a factor of 0.84 difference between the
c                         definition of collisionality in the 1995 published
c                         paper and the coded implementation provided by
c                         B. Dorland
c        switches(5) = 0  Do not relax restrictions on the value of znu
c                    = 1  Allow the value of znu to be larger than 10.0
c
c        switches(30), switches(31), switches(32)
c                       Reserved as output channel numbers
c
c        switches(j) should be defaulted to zero.
c
c     Dimensionless input parameters:
c        znine:     Ratio of thermal deuterium density to electron density
c        zncne:     Ratio of thermal carbon density to electron density
c        znbne:     Ratio of beam deuterium density to electron density
c        zrlt:      Ratio of major radius to ion temperature scale length
c        zrln:      Ratio of major radius to electron density scale length
c        zq:        Plasma safety factor
c        zshat:     Magnetic shear (defined as r/q dq/dr)
c        zeps:      Local inverse aspect ratio
c        zkappa:    Local plasma elongation
c        gtau:      Ratio of Ti/Te (used only when switches(2)=1)
c        gnu:       Local collisionality (used only when switches(2)=1)
c
c     Dimensioned input parameters:
c        ne19:      Electron density, in 10^19 meters^-3
c        tekev:     Electron temperature, in keV
c        tikev:     Deuterium ion temperature, in keV
c        rmajor:    Major radius, in meters
c        btesla:    Toroidal magnetic field, in Tesla
c        grhoi:     Local ion gyroradius, in meters (used when switches(2)=1)
c                   Input to replace calculated value for vthermi/omegaci
c        gvti:      Ion thermal speed, in m/s (used only when switches(2)=1),
c                      defined as sqrt(Ti/mi) -- nonstandard definition
c        omegaexb:  The E x B shearing rate, in inverse seconds
c
c     Dimensioned outputs:
c        chii:      Ion thermal diffusivity, in m^2/sec
c        chie:      Electron thermal diffusivity, in m^2/sec
c
c     Dimensionless outputs:
c        zchiicyc:  \chi_i, normalized by L_n/(\rho_i^2 v_{ti})
c        zchii1:    Hydrogenic \chi_i, normalized by R/(\rho_i^2 v_{ti})
c        zchii2:    Carbon-driven \chi_i, normalized by R/(\rho_i^2 v_{ti})
c        zchie1:    Hydrogenic \chi_e, normalized by R/(\rho_i^2 v_{ti})
c        zchie2:    Carbon-driven \chi_e, normalized by R/(\rho_i^2 v_{ti})
c        zrlt1:     Normalized threshold gradient for hydrogenic mode
c        zrlt2:     Normalized threshold gradient for carbon-driven mode
c        ierr:      Error code = 0 if all is OK
c                              = 1 if inputs are outside the ranges prescribed
c                                    in Kotschenreuther, et al,
c                                    Physics of Plasmas, v2, p2381 (1995)
c
c
c     Declare variables
c     NAMING CONVENTION: Dimensionless variables begin with a 'z'
c                        Dorland-style input variables begin with a 'g'
c
      INTEGER switches(32), ierr
      REAL znine, zncne, znbne, zrlt, zrln, zq, zshat, zeps,
     &       ne19, tekev, tikev, rmajor, rhoi, vthermi, omegaexb,
     &       ztau, znu, zzeffstr, zrlnstr, ztaub, zkappa,
     &       zchii, zchii1, zchii2, zscriptg1, zscriptg2, zx,
     &       zscriptz, zf, zg, zh, zscriptd, zscripte,
     &       zrlt1, zrlt2, zchie, zchie1, zchie2, zscriptj,
     &       btesla, chie, chii, omegaci, zcorrect, zkapfact,
     &       grhoi, gvti, gtau, gnu, zscripty, gamma,
     &       zc1, zc2, zc3, zscriptg, zchiicyc, zscripts
      CHARACTER*35 calcflag
c
c
c     Set some initial values
c
      ierr = 0
      omegaexb = abs(omegaexb)
      calcflag = '                                   '
      if ( switches(2) .eq. 1)
     &     calcflag = '* Calcd from gnu,gvti,gtau,grhoi * '
c
c
c     Check some of the inputs for validity
c
      if (zq .lt. 0.7) then
         zq = 0.7
         ierr = 1
      end if
      if (zq .gt. 8.0) then
         zq = 8.0
         ierr = 1
      end if
      if (zshat .lt. 0.5) then
         zshat = 0.5
         ierr = 1
      end if
      if (zshat .gt. 2.0) then
         zshat = 2.0
         ierr = 1
      end if
      if (zeps .lt. 0.1) then
         zeps = 0.1
         ierr = 1
      end if
      if (zeps .gt. 0.3) then
         zeps = 0.3
         ierr = 1
      end if
      if (zrln .lt. 1.0e-6) then
         zrln = 1.0e-6
         ierr = 1
      end if
      if (zrln .gt. 6.0) then
         zrln = 6.0
         ierr = 1
      end if
c
c
c     Initial calculations (check tau, nu, zeff for validity)
c
      if ( switches(2) .eq. 0 ) then
         vthermi = sqrt(tikev*1.6022e-16/3.3436e-27)
         omegaci = 1.6022e-19 * btesla / 3.3436e-27
         rhoi = vthermi/omegaci
         if ( (tikev .lt. 1.0e-4) .or. (tekev .lt. 1.0e-4) ) then
            znu = 10.0
            ierr = 1
         else
            znu = 2.1*rmajor*ne19/(tekev**1.5 * tikev**0.5)
            if ((znu .lt. 0.5) .and. (switches(3) .eq. 0)) then
               znu = 0.5
               ierr = 1
            end if
            if ( (znu .gt. 10.0) .and. (switches(5) .eq. 0) ) then
               znu = 10.0
               ierr = 1
            end if
         end if
         if ( tekev .lt. 1.0e-4 ) then
            ztau = 4.0
            ierr = 1
         else
            ztau = tikev/tekev
            if (ztau .lt. 0.5) then
               ztau = 0.5
               ierr = 1
            end if
            if (ztau .gt. 4.0) then
               ztau = 4.0
               ierr = 1
            end if
         end if
         gvti = vthermi
         grhoi = rhoi
         gtau = ztau
         gnu = znu
      end if
c
      if ( switches(2) .eq. 1 ) then
         vthermi = gvti
         rhoi = grhoi
         if (rhoi .lt. 1.0e-10) then
            rhoi = 1.0e-10
            ierr = 1
         end if
         znu = gnu
         if (switches(4) .eq. 1) znu = 0.84*gnu
         if ((znu .lt. 0.5) .and. (switches(3) .eq. 0)) then
            znu = 0.5
            ierr = 1
         end if
         if ( (znu .gt. 10.0) .and. (switches(5) .eq. 0) ) then
            znu = 10.0
            ierr = 1
         end if
         ztau = gtau
         if (ztau .lt. 0.5) then
            ztau = 0.5
            ierr = 1
         end if
         if (ztau .gt. 4.0) then
            ztau = 4.0
            ierr = 1
         end if
         tikev = 3.3436e-27*gvti*gvti/(1.6022e-16)
         tekev = tikev/ztau
         omegaci = gvti/grhoi
         btesla = omegaci*3.3436e-27/1.6022e-19
         ne19 = znu*tekev**1.5*tikev**0.5/(2.1*rmajor)
      end if
c
      ztaub = ztau/(1.0 - znbne)
      zzeffstr = 1.0 + 30.0*zncne/(znine + 6.0*zncne)
      if (zzeffstr .lt. 1.0) then
         zzeffstr = 1.0
         ierr = 1
      end if
      if (zzeffstr .gt. 4.0) then
         zzeffstr = 4.0
         ierr = 1
      end if
      zrlnstr = min(6.0,zrln)
c
c
c     Are we using the 1994 IFS/PPPL model?
c     If so, go to the start of that section
c
      if ( switches(3) .eq. 1 ) go to 2000
c
c     The 1995 IFS/PPPL model
c
c     The threshold temperature gradient for the hydrogenic mode
c
      zf = 1. - ( 0.942*zzeffstr**0.516/zshat**0.671 )*
     &          ( 2.95*zeps**1.257/znu**0.235 - 0.2126 )
      zg = ( 0.671 + 0.570*zshat - 0.189*zrlnstr )**2 + 0.392 +
     &          0.335*zrlnstr - 0.779*zshat + 0.210*zshat*zshat
      zh = 2.46 * (zzeffstr/2.0)**0.7 * ztaub**0.52 *
     &          (1.0 + 2.78/(zq*zq))**0.26
c
      zrlt1 = zf*zg*zh
c
 
c
c     The threshold temperature gradient for the carbon mode
c
      zscriptd = max( 1.0, 3.0-0.6667*zrlnstr )
      zscripte = 1.0 + 6.0*max( 0.0, 2.9-zzeffstr )
c
      zrlt2 = 0.75*(1.0+ztaub)*(1.0+zshat)*zscriptd*zscripte
c
 
c
c     Calculate zkapfact (=1.0 when kappa=1.0 -- circular)
c
      zkapfact = 1.0  / ( 1. + ((zkappa-1.)*zq/3.6)**2 )
c
c     Calculate the shear flow stabilization factor
c
      gamma = (vthermi/rmajor)*(0.25/(ztau*(1.+0.5*zshat*zshat)))*
     &     ((1.+3.*max(0.0,zeps-0.16667)) / (1.+max(0.,zq-3.)/15.))*
     &     (zrlt-zrlt1)
      zscripts = 0.0
      if (gamma .gt. omegaexb) zscripts = 1.0 - omegaexb/gamma
c
c     Calculate zchii1
c
      zx = zrlt - zrlt1
      zscriptg1 = zx
      if ( zx .lt. 0.0 ) zscriptg1 = 0.0
      if ( zx .gt. 1.0 ) zscriptg1 = sqrt(zx)
c
      zscriptz = min( 1.0, (3.0/zzeffstr)**1.8 )
c
      zchii1 = ( 11.8*zq**1.13/ztaub**1.07) / (1.0 + zshat**0.84) *
     &        (1.0 + 6.72*zeps/(zq**0.96*znu**0.26)) *
     &        zscriptz * zscriptg1 * zkapfact * zscripts
c
c     Calculate zchii2
c
      zx = zrlt - zrlt2
      zscriptg2 = zx
      if ( zx .lt. 0.0 ) zscriptg2 = 0.0
      if ( zx .gt. 1.0 ) zscriptg2 = sqrt(zx)
c
      zchii2 = (7.88)/( ztaub**0.8 * (1.0 + zshat) )*
     &        max( 0.25, zzeffstr-3.0 ) * zscriptg2 *
     &        zkapfact * zscripts
c
c     Now calculate chii
c
      zchii = max(zchii1,zchii2)
      chii = (rhoi*rhoi*vthermi/rmajor)*zchii
c
c     Get chii normalized in Cyclone units
c
      zchiicyc = zchii/zrln
c
c
c     Calculate chie correction factor
c
      zcorrect = (7.0 - zzeffstr)/6.0
c
c     Calculate zscriptj
c
      zscriptj = max ( 2.0, (1.0 + 0.3333*zrlnstr) )
c
c     Calculate zchie1 and zchie2
c
      zchie1 = zchii1 * 0.72 * max( 0.1667, zeps ) * znu**0.14 *
     &        (zq/zshat)**0.3 * ztau**0.4 *
     &        ( 1. + 0.333 * zrlnstr)
c
      zchie2 = zchii2 * 0.263 * ztau * znu**0.22 * zscriptj
c
c     Finally, calculate chie
c
      zchie = max(zchie1,zchie2) * zcorrect
      chie = (rhoi*rhoi*vthermi/rmajor)*zchie
c
c
c     (Optional) Print out page of calculation results
c
      if ( switches(1) .eq. 1 ) then
 
         if ( switches(30) .gt. 0) then
            write (switches(30),1000) chii, ierr, chie,
     &           rhoi, vthermi,
     &           omegaci, calcflag,
     &           tikev, calcflag,
     &           tekev, calcflag,
     &           ne19, calcflag,
     &           btesla, calcflag,
     &           omegaexb,
     &           zchiicyc,
     &           zchii1, zchii2,
     &           zchie1, zchie2, zrlt, zrlt1, zrlt2,
     &           zkapfact, zcorrect,
     &           gamma, zscripts,
     &           zscriptj, zscriptz,
     &           zscriptg1, zscriptg2,
     &           zscriptd, zscripte, zf, zg, zh,
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(31) .gt. 0) then
            write (switches(31),1000) chii, ierr, chie,
     &           rhoi, vthermi,
     &           omegaci, calcflag,
     &           tikev, calcflag,
     &           tekev, calcflag,
     &           ne19, calcflag,
     &           btesla, calcflag,
     &           omegaexb,
     &           zchiicyc,
     &           zchii1, zchii2,
     &           zchie1, zchie2, zrlt, zrlt1, zrlt2,
     &           zkapfact, zcorrect,
     &           gamma, zscripts,
     &           zscriptj, zscriptz,
     &           zscriptg1, zscriptg2,
     &           zscriptd, zscripte, zf, zg, zh,
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(32) .gt. 0) then
            write (switches(32),1000) chii, ierr, chie,
     &           rhoi, vthermi,
     &           omegaci, calcflag,
     &           tikev, calcflag,
     &           tekev, calcflag,
     &           ne19, calcflag,
     &           btesla, calcflag,
     &           omegaexb,
     &           zchiicyc,
     &           zchii1, zchii2,
     &           zchie1, zchie2, zrlt, zrlt1, zrlt2,
     &           zkapfact, zcorrect,
     &           gamma, zscripts,
     &           zscriptj, zscriptz,
     &           zscriptg1, zscriptg2,
     &           zscriptd, zscripte, zf, zg, zh,
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
 1000    format (/,' Results from the 1995 IFS/PPPL Model',/,
     &       /,' Dimensional:',
     &       /,' chi-i =',1pe21.14,' (m^2/s)','   error code =',i2,
     &       /,' chi-e =',1pe21.14,' (m^2/s)',/,
     &       /,' rho-i =',1pe21.14,' meters',
     &       /,' v-therm-i =',1pe21.14,' m/s',
     &       /,' omega-ci =',1pe21.14,' rad/s  ',a35,
     &       /,' T-i =',1pe21.14,' keV         ',a35,
     &       /,' T-e =',1pe21.14,' keV         ',a35,
     &       /,' n-e =',1pe21.14,' 10^19 m^-3  ',a35,
     &       /,' B-tor =',1pe21.14,' Tesla     ',a35,
     &       /,' omega-ExB =',1pe21.14,' rad/s  ',/,
     &       /,' Dimensionless:',
     &       /,' chi-i =',1pe21.14,' in Cyclone units',
     &       /,' zchi-i1 =',1pe21.14,' zchi-i2 =',1pe21.14,
     &       /,' zchi-e1 =',1pe21.14,' zchi-e2 =',1pe21.14,/,
     &       /,' R/L_T =',1pe21.14,
     &       /,' R/L_Tcrit1 =',1pe21.14,' R/L_Tcrit2 =',1pe21.14,/,
     &       /,'  zkapfact =',1pe21.14,'  zcorrect =',1pe21.14,
     &       /,' gamma s^-1=',1pe21.14,'  zscripts =',1pe21.14,
     &       /,'  zscriptj =',1pe21.14,'  zscriptz =',1pe21.14,
     &       /,' zscriptg1 =',1pe21.14,' zscriptg2 =',1pe21.14,
     &       /,'  zscriptd =',1pe21.14,'  zscripte =',1pe21.14,/,
     &       /,'       zf =',1pe21.14,'      zg =',1pe21.14,
     &       /,'       zh =',1pe21.14,'    ztau =',1pe21.14,
     &       /,'    ztaub =',1pe21.14,'     znu =',1pe21.14,
     &       /,' zzeffstr =',1pe21.14,' zrlnstr =',1pe21.14,
     &       /)
 
      end if
c
      go to 3000
c
c
c     The 1994 IFS/PPPL model
c
 2000 continue
c
      zc1 = 1.26 + abs(zzeffstr-3.)*
     &             (-0.27+0.075*zzeffstr-0.044*zzeffstr*zzeffstr)
      zc1 = max(0.57,zc1)
c
      zc2 = 2.0*(3.5-zzeffstr) + 0.1/ztau
      if ( zzeffstr .lt. 3.0) zc2 = 1.0 + 0.1/ztau
      if ( zzeffstr .gt. 3.5) zc2 = (0.1 + 0.2*(zzeffstr-3.5))/ztau
c
      zc3 = 0.61 - 0.27/( 1.0 + exp( 8.0*(3.3-zzeffstr) ) )
c
c
c     Calculate zscripty
c
      zscripty = abs(0.1976-0.4550*zshat+0.1616*zrln)**0.769
     &           + 0.7813 + 0.2762*zshat + 0.3967*zshat*zshat
c
c     Calculate the critical temperature gradient
c
      zrlt1 = 2.778 * zc1 * ztaub**zc3 * sqrt(0.5+1./zq) *
     &         abs(1. - 0.85*zeps/zshat**0.25) * zscripty
c
c
c     Calculate zscriptg and chii
c
      zx = zrlt - zrlt1
      zscriptg = zx
      if ( zx .lt. 0.0 ) zscriptg = 0.0
      if ( zx .gt. 1.0 ) zscriptg = sqrt(zx)
c
      zchii1 = 14.*zq*zc2/(ztaub*(2.+zshat))
     &         * (1. + zeps/3.) * zscriptg
      chii = zchii1 * rhoi*rhoi*vthermi/rmajor
c
c     Get chii normalized in Cyclone units
c
      zchiicyc = zchii1/zrln
c
c     Zero out zchii2, zchie2 and zrlt2 (they are not used in the 1994 model)
c
      zchii2 = 0.0
      zchie2 = 0.0
      zrlt2 = 0.0
c
c
c     Calculate chie and zchie
c
      zchie1 = zchii1 * 0.27 * ztau**0.7
      chie = chii * 0.27 * ztau**0.7
c
c
c     (Optional) Print out page of calculation results
c
      if ( switches(1) .eq. 1 ) then
 
         if ( switches(30) .gt. 0) then
            write (switches(30),5000) chii, ierr, chie,
     &           rhoi, vthermi,
     &           omegaci, calcflag,
     &           tikev, calcflag,
     &           tekev, calcflag,
     &           ne19, calcflag,
     &           btesla, calcflag,
     &           zchiicyc, zrlt, zrlt1,
     &           zscriptg, zc1, zc2, zc3,
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(31) .gt. 0) then
            write (switches(31),5000) chii, ierr, chie,
     &           rhoi, vthermi,
     &           omegaci, calcflag,
     &           tikev, calcflag,
     &           tekev, calcflag,
     &           ne19, calcflag,
     &           btesla, calcflag,
     &           zchiicyc, zrlt, zrlt1,
     &           zscriptg, zc1, zc2, zc3,
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
         if ( switches(32) .gt. 0) then
            write (switches(32),5000) chii, ierr, chie,
     &           rhoi, vthermi,
     &           omegaci, calcflag,
     &           tikev, calcflag,
     &           tekev, calcflag,
     &           ne19, calcflag,
     &           btesla, calcflag,
     &           zchiicyc, zrlt, zrlt1,
     &           zscriptg, zc1, zc2, zc3,
     &           ztau, ztaub, znu, zzeffstr, zrlnstr
         end if
 
 5000    format (/,' Results from the 1994 IFS/PPPL Model',/,
     &       /,' Dimensional:',
     &       /,' chi-i =',1pe21.14,' (m^2/s)','   error code=',i2,
     &       /,' chi-e =',1pe21.14,' (m^2/s)',/,
     &       /,' rho-i =',1pe21.14,' meters',
     &       /,' v-therm-i =',1pe21.14,' m/s',
     &       /,' omega-ci =',1pe21.14,' rad/s  ',a35,
     &       /,' T-i =',1pe21.14,' keV         ',a35,
     &       /,' T-e =',1pe21.14,' keV         ',a35,
     &       /,' n-e =',1pe21.14,' 10^19 m^-3  ',a35,
     &       /,' B-tor =',1pe21.14,' Tesla     ',a35,/,
     &       /,' Dimensionless:',
     &       /,' chi-i =',1pe21.14,' in Cyclone units',
     &       /,' R/L_T =',1pe21.14,'  R/L_Tcrit1 =',1pe21.14,/,
     &       /,'  zscriptg =',1pe21.14,/,
     &       /,' zc1 =',1pe21.14,' zc2 =',1pe21.14,
     &       /,' zc3 =',1pe21.14,/,
     &       /,'    ztau =',1pe21.14,'    ztaub =',1pe21.14,
     &       /,'     znu =',1pe21.14,' zzeffstr =',1pe21.14,
     &       /,' zrlnstr =',1pe21.14,
     &       /)
 
      end if
c
c
 3000 continue
      return
      end
c
c--------1---------2---------3---------4---------5---------6---------7-
c
c Modification History
c
c ajr 08-Jul-98 Fixed reference rewoldt98.  Changed an if condition on
c         switches(2) from .ne.1 to .eq.0.  More wordsmithing.
c ajr 07-Jul-98 Added table of contents, more section headings.
c         Changed variable declarations so that compiler controls precision.
c         Added more digits to diagnostic output (1pe12.4 -> 1pe21.14).
c         Added some opening remarks re: the 1994 model, including citation.
c         Added better descriptions of I/O parameters.  Some wordsmithing.
c
c--------1---------2---------3---------4---------5---------6---------7-
