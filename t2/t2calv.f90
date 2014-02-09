!C-----------------------------------------------------------------------
!C
!C MODULE T2CALV
!C 
!C
!C      CALCULATION OF COEFFICIENTS 
!C      OF ADVECTION-DIFFUSION EQUATIONS 
!C      FOR MULTI-FLUID TRANSPORT MODEL
!C 
!C 
!C                       2014-02-04
!C
!C----------------------------------------------------------------------
MODULE T2CALV
  
  USE T2CNST, ONLY:i0ikind,i0rkind
  USE T2COMM, ONLY:i0smax
  
  IMPLICIT NONE
  
  INTEGER(i0ikind)::&
       i0midi
  REAL(   i0rkind)::&
       d0cogrr,d0cogrp,d0cogpp,d0cogtt,d0sqrtg,&
       d0ctgrr,d0ctgrp,d0ctgpp,d0ctgtt,d0ugr,d0ugp,&
       d0ctbp,d0ctbpi,d0cobt,d0ctbt,d0bp2i,d0bp2,d0bt2,d0bb,d0bb2
  
  PRIVATE
  
  PUBLIC T2_CALV
  
CONTAINS
  
  SUBROUTINE T2_CALV
    
    USE T2COMM, ONLY: i0mmax
    
    DO i0midi = 1, i0mmax
       
       CALL T2CALV_PQ
       CALL T2CALV_MS
       CALL T2CALV_AV
       CALL T2CALV_AT
       CALL T2CALV_DT
       CALL T2CALV_GV
       CALL T2CALV_GT
       CALL T2CALV_ES
       CALL T2CALV_EV
       CALL T2CALV_ET
       CALL T2CALV_SS
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2_CALV
  
  !C---------------------------------------------------------
  !C
  !C CALCULATION OF FUNDAMENTAL PHYSICAL QUANTITIES 
  !C 
  !C     2014-01-26 H.SETO
  !C                                        
  !C---------------------------------------------------------
  
  SUBROUTINE T2CALV_PQ
    
    USE T2CNST
    
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         i0xa,i0vmax,d0rmjr,d0qc,d0qs,&
         i2crt,d2mtrc,&
         d2xvec_befor,d2mfc1,d2rzm,d2jm1,d2ws,&
         d1ee,d1mm,d1nn,d1ni,d1pp,d1pi,d1tt,d1ti,d1pa,d1pz,&
         d1ur,d1up,d1ut,d1ub,d1u2,&
         d1qr,d1qp,d1qt,d1qb,d1vb,d1vt,d1hex,&
         d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4,&
         d2x,d2y,d2z,d2bcf
    
    USE LIBT2, ONLY:integ_f
    
    INTEGER(i0ikind)::&
         i0xid1d,i0xid2d,i0widi,i0vidi,i0sidi,i0sidj,i0sidk,i0nflag
    
    REAL(   i0rkind)::&
         d0mm_a,d0ee_a,d0ee_b,d0vti_a,&
         d0nn_a,d0ni_a,d0nn_b,d0ni_b,&
         d0pp_a,d0pi_a,       d0pi_b,&
         d0tt_a,d0ti_a,&
         d0clog_ab,d0cfreq_ab,d0cfreq_ac,&
         d0x_ab,d0y_ab,d0z_ab,d0zi_ab,&
         d0bmem00_ab,d0bmem01_ab,            d0bmem11_ab, &
         d0bmem00_ac,d0bmem01_ac,d0bmem10_ac,d0bmem11_ac, &
         d0bmen00_ab,d0bmen01_ab,d0bmen10_ab,d0bmen11_ab, &
         d0nfcl11_ab,d0nfcl12_ab,d0nfcl21_ab,d0nfcl22_ab, &
         d0nfcf1_ab, d0nfcf2_ab, d0nfcf3_ab, d0nfcf4_ab,  &
         d0nvcm1_a,  d0nvcm2_a,  d0nvcm3_a,               &
         d0nvcc1_a,  d0nvcc2_a,  d0nvcc3_a,  d0nvcc4_a
    
    REAL(i0rkind)::&
         d0rzcr,d0mfcr,d0mfcp,d0psip,d0g1,d0g2,d0q0,d0liar,&
         d0cps,d0cbn,d0tcr,d0wv1,d0wv2,d0wv3,d0wv4,&
         d0k11,  d0k12,  d0k22,&
         d0k11ps,d0k12ps,d0k22ps,&
         d0k11bn,d0k12bn,d0k22bn
    
    REAL(i0rkind),DIMENSION(1:i0smax)::&
         d1fr,d1fb,d1ft,d1vti,&
         d1nvcm1,d1nvcm2,d1nvcm3
    REAL(i0rkind),DIMENSION(1:i0smax,1:i0smax)::&
         d2clog,d2cfreq,&
         d2bmem00,d2bmem01,d2bmem10,d2bmem11,&
         d2bmen00,d2bmen01,d2bmen10,d2bmen11,&
         d2nfcl11,d2nfcl12,d2nfcl21,d2nfcl22
    REAL(i0rkind)::d0err
    
    !C
    !C d0psip: DERIVATIVE POLOIDAL FLUX FUNCTION
    !C         WITH RESPECT TO RHO               : \psi'
    !C d0pcf : POLOIDAL CURRENT FUNCTION         : I
    !C d1nn  : PARTICLE DENSITY                  : n_{a}
    !C d1pp  : PRESSURE                          : p_{a}   
    !C
    !C d1ur  : CT RADIAL   FLOW                  : u_{a}^{\rho} 
    !C d1up  : CT POLOIDAL FLOW                  : u_{a}^{\chi}
    !C d1ut  : CO TOROIDAL FLOW                  : u_{a\zeta}
    !C d1ub  :    PARALLEL FLOW                  : u_{a\para}
    !C
    !C d1qr  : CT RADIAL   TOTAL HEAT FLUX       : Q_{a}^{\rho}
    !C d1qp  : CT POLOIDAL TOTAL HEAT FLUX       : Q_{a}^{\chi}
    !C d1qt  : CO TOROIDAL TOTAL HEAT FLUX       : Q_{a\zeta}
    !C d1qb  :    PARALLEL TOTAL HEAT FLUX       : Q_{a\para}
    !C 
    !C d1fr  : CT RADIAL   FLOW                  : n_{a}u_{a}^{\rho} 
    !C d1fp  : CT POLOIDAL FLOW                  : n_{a}u_{a}^{\chi}
    !C d1ft  : CO TOROIDAL FLOW                  : n_{a}u_{a\zeta}
    !C d1fb  :    PARALLEL PARTICLE FLUX         : n_{a}n_{u\para}
    !C
    !C d01vb : PARALLEL TOTAL HEAT FLOW FUNCTION : Q_{a\para}/p_{a}
    !C 
    !C * CO = COVARIANT, CT = CONTRAVARIANT
    
    !C
    !C MASS AND ELECTRIC CHARGE
    !C

    DO i0sidi = 1,i0smax
       d1mm(i0sidi) = d1pa(i0sidi)*d0amp
       d1ee(i0sidi) = d1pz(i0sidi)*d0aee
    ENDDO
    
    !C
    !C INITIALIZE     
    !C
    
    i0xid1d = i2crt( 3,i0midi)
    i0xid2d = i2crt( 2,i0midi)
    d0rzcr  = d2rzm( 1,i0midi)
    d0mfcr  = d2mfc1(1,i0midi)
    d0mfcp  = d2mfc1(2,i0midi)

    !C
    !C CONVERT VARIABLES TO SI-UNIT
    !C
       
    DO i0sidi = 1, i0smax
       i0nflag = 0
       i0vidi =  8*i0sidi - 3
       
       d0nn_a = d2xvec_befor(i0vidi+1,i0xid2d)*1.D-20
       d0pp_a = d2xvec_befor(i0vidi+5,i0xid2d)*1.D-23/d0aee
       
       d1nn(i0sidi) = d0nncst*d2xvec_befor(i0vidi+1,i0xid2d)
       d1fr(i0sidi) = d0frcst*d2xvec_befor(i0vidi+2,i0xid2d)
       d1fb(i0sidi) = d0fbcst*d2xvec_befor(i0vidi+3,i0xid2d)
       d1ft(i0sidi) = d0ftcst*d2xvec_befor(i0vidi+4,i0xid2d)
       d1pp(i0sidi) = d0ppcst*d2xvec_befor(i0vidi+5,i0xid2d)
       d1qr(i0sidi) = d0qrcst*d2xvec_befor(i0vidi+6,i0xid2d)
       d1qb(i0sidi) = d0qbcst*d2xvec_befor(i0vidi+7,i0xid2d)
       d1qt(i0sidi) = d0qtcst*d2xvec_befor(i0vidi+8,i0xid2d)
       
       IF(        d0nn_a .GT.0.D0 )THEN
          d1ni(i0sidi) = 1.D0/d1nn(i0sidi)
       ELSEIF(ABS(d0nn_a).LE.1.D-4)THEN 
          d1ni(i0sidi) = 0.D0
          d1nn(i0sidi) = 0.D0
       ELSE
          i0nflag = 1
       ENDIF
       
       IF(        d0pp_a .GT.0.D0 )THEN
          d1pi(i0sidi) = 1.D0/d1pp(i0sidi)
       ELSEIF(ABS(d0pp_a).LE.1.D-4)THEN
          d1pi(i0sidi) = 0.D0
          d1pp(i0sidi) = 0.D0
       ELSE
          i0nflag = 1
       ENDIF
       
       IF(i0nflag.EQ.1)THEN
          WRITE(6,*)'NEGATIVE  DENSITY or PRESSURE'
          WRITE(6,*)'SPECIS=',i0sidi,'NODE=',i0xid2d,&
               'N',d0nn_a,'10^{20} /m3','P=',d0pp_a/d0nn_a,'keV'
          STOP       
       ENDIF
    ENDDO
    


    !C
    !C GEOMETRICAL COEFFICIENTS
    !C
    
    !C 
    !C d0sqrtg : \sqrt{g}
    !C
    !C d0cogrr : g_{\rho\rho}
    !C d0cogrp : g_{\rho\chi}
    !C d0cogpp : g_{\chi\chi}
    !C d0cogtt : g_{\zeta\zeta}
    !C
    !C d0ctgrr : g^{\rho\rho}
    !C d0ctgrp : g^{\rho\chi}
    !C d0ctgpp : g^{\chi\chi}
    !C d0ctgtt : g^{\zeta\zeta} 
    !C
    !C d0ugr   : CT RADIAL   FLAME MOVING VELOCITY      : u_{g}^{\rho}
    !C d0ugr   : CT POLOIDAL FLAME MOVING VELOCITY      : u_{g}^{\rho}
    !C
    !C d0ctbp  : CT POLOIDAL MAGNETIC FIELD             : B^{\chi}
    !C d0ctbpi : INVERSE CONTRA POLOIDAL MAGNETIC FIELD : 1/B^{\chi}  
    !C d0bp2   : SQUARED INTENSITY 
    !C                   of POLOIDAL MAGENTIC FIELD     : Bp^{2}
    !C d0bt2   : SQUARED INTENSITY 
    !C                   of TOROIDAL MAGENTIC FIELD     : Bt^{2}
    !C d0bp2i  : INVERSE SQUARED INTENSITY 
    !C                   of POLOIDAL MAGNETIC FIELD     : 1/Bp^{2}
    !C d0bb    : INTENSITY OF MAGNETIC FIELD            : B 
    !C d0bb2   : SQUARED INTENSITY of MAGNETIC FIELD    : B^{2} 
    !C
    
    d0psip = d0mfcst*d2xvec_befor(1,i0xid1d)
    d0cobt = d0btcst*d2xvec_befor(2,i0xid1d)
    d0sqrtg = d2jm1(1,i0midi)
    d0ctgrr = d2jm1(2,i0midi)
    d0ctgrp = d2jm1(3,i0midi)
    d0ctgpp = d2jm1(4,i0midi)
    d0ctgtt = d2jm1(5,i0midi)
    
    d0cogtt =  1.D0/d0ctgtt
    
    ! modified by H.SETO 2014-02-03
    IF(d0sqrtg.GT.0.D0)THEN
       !C d0wv1  : WORKING VARIABLE: \sqrt{g}^{2}/R^{2}
       d0wv1   =  (d0sqrtg**2)*d0ctgtt
       d0cogrr =  d0ctgpp*d0wv1
       d0cogrp = -d0ctgrp*d0wv1
       d0cogpp =  d0ctgrr*d0wv1
       d0ctbp  =  d0psip/d0sqrtg
       d0ctbpi = 1.D0/d0ctbp
    ELSE
       d0cogrr = 0.D0
       d0cogrp = 0.D0
       d0cogpp = 0.D0
       d0ctbp  = 0.D0
       d0ctbpi = 0.D0
    ENDIF
    
    d2mtrc(1,i0midi) = d0cogrr
    d2mtrc(2,i0midi) = d0cogpp

    
    d0ugr = 0.D0
    d0ugp = 0.D0
    
    d0ctbt = d0cobt*d0ctgtt
    
    d0bp2 = (d0ctbp**2)*d0cogpp
    d0bt2 = (d0cobt**2)*d0ctgtt
    d0bb2 = d0bp2 + d0bt2
    d0bb  = SQRT(d0bb2)

    IF(d0bp2.GT.0.D0)THEN
       d0bp2i= 1.D0/d0bp2
    ELSE
       d0bp2i = 0.D0
    ENDIF
    
    d0g1 =  d0ctbp*d0bp2i*d0bb
    d0g2 = -d0ctbp*d0bp2i*d0ctbt
    
    DO i0sidi = 1, i0smax 
       
       d1ur(i0sidi) = d1fr(i0sidi)*d1ni(i0sidi)
       d1ub(i0sidi) = d1fb(i0sidi)*d1ni(i0sidi)
       d1ut(i0sidi) = d1ft(i0sidi)*d1ni(i0sidi)
       
       d1up(i0sidi) = d0g1*d1ub(i0sidi) + d0g2*d1ut(i0sidi)
       d1u2(i0sidi) = d1ub(i0sidi)*d1ub(i0sidi)
       
       d1qp(i0sidi) = d0g1*d1qb(i0sidi) + d0g2*d1qt(i0sidi)
       d1vb(i0sidi) = d1qb(i0sidi)*d1pi(i0sidi)
       
       d1tt(i0sidi) = d1pp(i0sidi)*d1ni(i0sidi)
       d1ti(i0sidi) = d1nn(i0sidi)*d1pi(i0sidi)
       
    ENDDO
    
    !C SET WORKING SCALAR FOR DIFFERENTIAL (GROBAL)
    !C
    !C D2WS(1   ,:) : MAGNETIC FIELD INTENSITY    : B   
    !C D2WS(2   ,:) : MAJOR RADIUS                : R   
    !C D2WS(2N+1,:) : PARALLEL FLOW               : u_{a\para}  
    !C D2WS(2N+2,:) : PARALLEL HEAT FLOW FUNCTION : Q_{a\para}/p_{a}
    !C
    
    d2ws(         1,i0midi) = d0bb
    d2ws(         2,i0midi) = d0rzcr
    
    DO i0sidi = 1, i0smax
       i0widi = 2*i0sidi
       d2ws(i0widi+1,i0midi) = d1ub(i0sidi)
       d2ws(i0widi+2,i0midi) = d1vb(i0sidi)
    ENDDO
    
    !C
    !C FOR COLLISION AND VISCOUS TERMS
    !C
    
    !C THERMAL VELOCITY [m/s]
    !C v_{ta} = \sqrt{2T_{a}/M_{a}}
    !C
    !C      checked 2014-02-03 by H.Seto
    !C
    DO i0sidi = 1, i0smax
       d0mm_a = d1mm(i0sidi)
       d0tt_a = d1tt(i0sidi)
       d0ti_a = d1ti(i0sidi)
       d1vt( i0sidi)= SQRT(2.0D0/d0mm_a*d0tt_a)
       d1vti(i0sidi)= SQRT(0.5D0*d0mm_a*d0ti_a)
    ENDDO
    
    !C COULOMB LOGARITHM
    !C ref: NRL PLASMA FORMULARY 2011
    !C FOR DEBUG lnA = 17
    !C d2clog: COULOMB LOGARITHM
    
    DO i0sidi = 1, i0smax
    DO i0sidj = 1, i0smax
       d2clog(i0sidi,i0sidj) = 17.D0
    ENDDO
    ENDDO
    
    !C
    !C BASIC COLLISION FREQUENCY  [1/s]
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002)  P.277
    !C
    !C      checked 2014-02-03 by H.Seto
    !C
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       d0clog_ab = d2clog(i0sidi,i0sidj)
       d0nn_b    = d1nn(i0sidj)
       d0ee_b    = d1ee(i0sidj)
       d0mm_a    = d1mm( i0sidi)
       d0ee_a    = d1ee( i0sidi)       
       d0vti_a   = d1vti(i0sidi)
       d2bcf(i0sidi,i0sidj)&
            = (d0nn_b*(d0ee_a**2)*(d0ee_b**2)*d0clog_ab*(d0vti_a**3))&
            / (4.D0*d0pi*(d0eps0**2)*(d0mm_a**2))
    ENDDO
    ENDDO

    !C
    !C COLLISION FREQUENCY [1/s]
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002)  P.277
    !C
    !C      checked 2014-02-03 by H.Seto
    !C
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       d2cfreq(i0sidi,i0sidj) &
            = d2bcf(i0sidi,i0sidj)/(0.75D0*SQRT(d0pi))
    ENDDO
    ENDDO
    !C 
    !C FOR FRICTION COEFFICIENTS
    !C 

    !C
    !C BRAGINSKII'S MATRIX ELEMENT OF COLLISION OPERATOR
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002)
    !C      P.278
    !C
    !C DEFINITION OF VARIABLES FOR DIMENSIONLESS PARAMETERS
    !C
    !C d2x   = x_ab = v_{tb}/v_{ta}
    !C d2y   = y_ab = m_{a} /m_{b}
    !C d2z   = t_ab = t_{a} /t_{b}
    !C d0wv2 = 1/SQRT(1 + x_{ab}^{2})
    !C
    !C DEFINITION OF VARIABLRS FOR BRAGINSKII MATRIX ELEMENTS
    !C 
    !C d2bmem00: M_{ab}^{00}
    !C d2bmem01: M_{ab}^{01}
    !C d2bmem10: M_{ab}^{10}
    !C d2bmem11: M_{ab}^{11}
    !C d2bmen00: N_{ab}^{00}
    !C d2bmen01: N_{ab}^{01}
    !C d2bmen10: N_{ab}^{10}
    !C d2bmen11: N_{ab}^{11}    
    !C

    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax 
       d2x(i0sidi,i0sidj) = d1vti(i0sidi)*d1vt(i0sidj)
       d2y(i0sidi,i0sidj) = d1mm( i0sidi)/d1mm(i0sidj)
       d2z(i0sidi,i0sidj) = d1tt( i0sidi)*d1ti(i0sidj)
    ENDDO
    ENDDO
    
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       
       d0x_ab  = d2x(i0sidi,i0sidj) 
       d0y_ab  = d2y(i0sidi,i0sidj)
       d0z_ab  = d2z(i0sidi,i0sidj)
       d0zi_ab = d2z(i0sidj,i0sidi)
       
       d0wv2   = 1.D0/SQRT(1.D0 + d0x_ab**2)
       
       d0bmem00_ab = -        (1.D0+d0y_ab)*(d0wv2**3)
       d0bmem01_ab = -  1.5D0*(1.D0+d0y_ab)*(d0wv2**5)
       d0bmem11_ab = - (3.25D0+4.D0*(d0x_ab**2)&
                     +  7.50D0*(d0x_ab**4))*(d0wv2**5)
       
       d2bmem00(i0sidi,i0sidj) = d0bmem00_ab
       d2bmem01(i0sidi,i0sidj) = d0bmem01_ab
       d2bmem10(i0sidi,i0sidj) = d0bmem01_ab
       d2bmem11(i0sidi,i0sidj) = d0bmem11_ab
       
       d2bmen00(i0sidi,i0sidj) = - d0bmem00_ab
       d2bmen01(i0sidj,i0sidi) = - d0bmem01_ab*d0x_ab*d0zi_ab
       d2bmen10(i0sidi,i0sidj) = - d0bmem01_ab
       d2bmen11(i0sidi,i0sidj) = 6.75D0*d0z_ab*(d0x_ab**2)*(d0wv2**5)
       
    ENDDO
    ENDDO
    
    !C
    !C PARALLEL FRICTION COEFFICIENTS
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002) P.239 
    !C 
    !C DEFINITION OF VARIABLRS FOR NEOCLASSICAL FRICTION COEFFICIENTS
    !C
    !C d2nfcl11: l^{ab}_{11}
    !C d2nfcl12: l^{ab}_{12}
    !C d2nfcl21: l^{ab}_{21}
    !C d2nfcl22: l^{ab}_{22}

    d0nfcl11_ab = 0.D0
    d0nfcl12_ab = 0.D0
    d0nfcl21_ab = 0.D0
    d0nfcl22_ab = 0.D0
    
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       
       d0cfreq_ab  = d2cfreq( i0sidi,i0sidj)
       
       d0bmen00_ab = d2bmen00(i0sidi,i0sidj)
       d0bmen01_ab = d2bmen01(i0sidi,i0sidj)
       d0bmen10_ab = d2bmen10(i0sidi,i0sidj)
       d0bmen11_ab = d2bmen11(i0sidi,i0sidj)
       
       d0nfcl11_ab = d0bmen00_ab*d0cfreq_ab
       d0nfcl12_ab = d0bmen01_ab*d0cfreq_ab
       d0nfcl21_ab = d0bmen10_ab*d0cfreq_ab
       d0nfcl22_ab = d0bmen11_ab*d0cfreq_ab
       
       IF(i0sidi.EQ.i0sidj)THEN
          
          DO i0sidk = 1, i0smax
             
             d0cfreq_ac  = d2cfreq( i0sidi,i0sidk)
             
             d0bmem00_ac = d2bmem00(i0sidi,i0sidk)
             d0bmem01_ac = d2bmem01(i0sidi,i0sidk)
             d0bmem10_ac = d2bmem10(i0sidi,i0sidk)
             d0bmem11_ac = d2bmem11(i0sidi,i0sidk)
             
             d0nfcl11_ab = d0nfcl11_ab + d0bmem00_ac*d0cfreq_ac
             d0nfcl12_ab = d0nfcl12_ab + d0bmem01_ac*d0cfreq_ac
             d0nfcl21_ab = d0nfcl21_ab + d0bmem10_ac*d0cfreq_ac
             d0nfcl22_ab = d0nfcl22_ab + d0bmem11_ac*d0cfreq_ac
             
          ENDDO
          
       ENDIF
       
       d0wv3 = d1mm(i0sidi)*d1nn(i0sidi)
       
       d2nfcl11(i0sidi,i0sidj) = d0nfcl11_ab*d0wv3 ! l_{11}^{ab}
       d2nfcl12(i0sidi,i0sidj) = d0nfcl12_ab*d0wv3 ! l_{12}^{ab}
       d2nfcl21(i0sidi,i0sidj) = d0nfcl21_ab*d0wv3 ! l_{21}^{ab}
       d2nfcl22(i0sidi,i0sidj) = d0nfcl22_ab*d0wv3 ! l_{22}^{ab}
       
    ENDDO
    ENDDO

!    DO i0sidi =1,i0smax
!    DO i0sidj =1,i0smax
!       print*,i0sidi,i0sidj
!       print*,d2nfcl11(i0sidi,i0sidj),d2nfcl12(i0sidi,i0sidj)
!       print*,d2nfcl21(i0sidi,i0sidj),d2nfcl22(i0sidi,i0sidj)
!    ENDDO
!    ENDDO


    !C PARALLEL FRICTION COEFFICIENTS 
    !C          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
    !C
    !C DEFINITION OF VARIABLRS 
    !C            FOR NEOCLASSICAL FRICTION COEFFICIENTS
    !C
    !C d2nfcf1: f^{ab}_{1}
    !C d2nfcf2: f^{ab}_{2}
    !C d2nfcf3: f^{ab}_{3}
    !C d2nfcf4: f^{ab}_{4}
    
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       
       d0nfcl11_ab = d2nfcl11(i0sidi,i0sidj)
       d0nfcl12_ab = d2nfcl12(i0sidi,i0sidj)
       d0nfcl21_ab = d2nfcl21(i0sidi,i0sidj)
       d0nfcl22_ab = d2nfcl22(i0sidi,i0sidj)
       
       d0ni_b = d1ni(i0sidj)
       d0pi_b = d1pi(i0sidj)
       
       d0wv4  = d1tt(i0sidi)/d1mm(i0sidi)
       
       d0nfcf1_ab =   (d0nfcl11_ab+d0nfcl12_ab)*d0ni_b 
       d0nfcf2_ab = - 0.4D0*d0nfcl12_ab*d0pi_b 
       d0nfcf3_ab = - (d0nfcl21_ab+d0nfcl22_ab)*d0ni_b
       d0nfcf4_ab =   0.4D0*d0nfcl22_ab*d0pi_b
       
       d0nfcf3_ab = d0wv4*(2.5D0*d0nfcf1_ab+d0nfcf3_ab)
       d0nfcf4_ab = d0wv4*(2.5D0*d0nfcf2_ab+d0nfcf4_ab)
       
       d2nfcf1(i0sidi,i0sidj) = d0nfcf1_ab
       d2nfcf2(i0sidi,i0sidj) = d0nfcf2_ab
       d2nfcf3(i0sidi,i0sidj) = d0nfcf3_ab
       d2nfcf4(i0sidi,i0sidj) = d0nfcf4_ab
       
    ENDDO
    ENDDO
    
    !C
    !C HEAT EXCHANGE COEFFICIENT
    !C CHECKED 2013-11-14
    !C
    
    DO i0sidi = 1, i0smax
       d1hex(i0sidi) = 0.D0
    ENDDO
    
    DO i0sidi = 1, i0smax
    DO i0sidj = 1, i0smax
       d0zi_ab       = d2z(i0sidj,i0sidi)
       d0cfreq_ab    = d2cfreq(i0sidi,i0sidj)
       d1hex(i0sidi) = d1hex(i0sidi)+1.5D0*(1.D0 - d0zi_ab)*d0cfreq_ab
    ENDDO
    ENDDO
    
    !C
    !C NEOCLASSICAL PARALLEL COEFFICIENTS
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002)  CHAP.12 AND P.279
    !C
    !C d0iar: INVERSE ASPECT RATIO: r/R
    !C d0tcr: TRAPPED/CIRCULATING RATE: f_t/f_c
    !C d0cbn: COEFFICIENT FOR BANANA REGIME
    !C d0cps: COEFFICIENT FOR PS     REGIME
    !C
    !C d0k11ps: MOD VISCOSITY COEFFICIENT IN PS REGIME : K^{11}_{aPS}
    !C d0k12ps: MOD VISCOSITY COEFFICIENT IN PS REGIME : K^{12}_{aPS}
    !C d0k22ps: MOD VISCOSITY COEFFICIENT IN PS REGIME : K^{22}_{aPS}
    !C
    !C d0k11bn: MOD VISCOSITY COEFFICIENT IN BN REGIME : K^{11}_{aBN}
    !C d0k12bn: MOD VISCOSITY COEFFICIENT IN BN REGIME : K^{12}_{aBN}
    !C d0k22bn: MOD VISCOSITY COEFFICIENT IN BN REGIME : K^{22}_{aBN}
    !C
    !C d0k11  : MOD VISCOSITY COEFFICIENT              : K^{11}_{a}
    !C d0k12  : MOD VISCOSITY COEFFICIENT              : K^{12}_{a}
    !C d0k22  : MOD VISCOSITY COEFFICIENT              : K^{22}_{a}
    !C
    !C d1nvc1 : NEOCLASSICAL VISCOSITY COEFFICIENT     : \mu_{1a}
    !C d1nvc2 : NEOCLASSICAL VISCOSITY COEFFICIENT     : \mu_{2a} 
    !C d1nvc3 : NEOCLASSICAL VISCOSITY COEFFICIENT     : \mu_{3a} 
    !C
    
    IF(d0mfcr.GT.0.D0)THEN
       
       d0q0 = d0ctbt*d0ctbpi

       d0liar = d0mfcr/d0rmjr     !C LOCAL INVERSE ASPECT RATIO (r/R0)
       d0tcr = 1.46D0*SQRT(d0liar)/(1.D0-1.46D0*SQRT(d0liar)) 
       
       DO i0xa = 1, i0smax
          
          d0mm_a = d1mm(i0xa)
          d0nn_a = d1nn(i0xa)
          d0pp_a = d1pp(i0xa)
          
          !C
          !C MODIFIED VISCOSITY COEFFICIENT IN PS REGIME
          !C
          
          CALL INTEG_F(fd0k11ps,d0k11ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k12ps,d0k12ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k22ps,d0k22ps,d0err,EPS=1.D-8,ILST=0)
          
          d0cps   = 0.4D0*d0pp_a*d0de
          d0k11ps = d0cps*d0k11ps
          d0k12ps = d0cps*d0k12ps
          d0k22ps = d0cps*d0k22ps

          !C
          !C MODIFIED VISCOSITY COEFFICIENT IN BN REGIME
          !C 
          
          CALL INTEG_F(fd0k11bn,d0k11bn,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k12bn,d0k12bn,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k22bn,d0k22bn,d0err,EPS=1.D-8,ILST=0)
       
          d0cbn = 2.D0*d0mm_a*d0nn_a*d0de*d0tcr*(d0q0**2)*(d0rzcr**2)
          d0cbn = d0cbn / (3.D0*(d0liar**2))
          
          d0k11bn = d0cbn*d0k11bn
          d0k12bn = d0cbn*d0k12bn
          d0k22bn = d0cbn*d0k22bn
          
          !C
          !C MODIFIED VISCOSITY COEFFICIENT IN INTER-REGIMES 
          !C
          
          d0k11 = d0k11ps*d0k11bn/(d0k11ps+d0k11bn)
          d0k12 = d0k12ps*d0k12bn/(d0k12ps+d0k12bn)
          d0k22 = d0k22ps*d0k22bn/(d0k22ps+d0k22bn)

          !C
          !C NEOCLASSICAL VISCOSITY COEFFICIENTS
          !C
          d1nvcm1(i0xa) = d0k11 
          d1nvcm2(i0xa) = d0k12 - 2.5D0*d0k11 
          d1nvcm3(i0xa) = d0k22 - 5.0D0*d0k12 + 6.25D0*d0k11
          
       ENDDO
       
       !C NEOCLASSICAL VISCOSITY COEFFICIENTS 
       !C          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
       !C
       !C DEFINITION OF VARIABLRS 
       !C            FOR NEOCLASSICAL VISCOSITY COEFFICIENTS
       !C
       !C d2nfcf1: f^{ab}_{1}
       !C d2nfcf2: f^{ab}_{2}
       !C d2nfcf3: f^{ab}_{3}
       !C d2nfcf4: f^{ab}_{4}
       
       DO i0sidi = 1, i0smax
          
          d0nvcm1_a = d1nvcm1(i0sidi)
          d0nvcm2_a = d1nvcm2(i0sidi)
          d0nvcm3_a = d1nvcm3(i0sidi)
          
          d0ni_a = d1ni(i0sidi)
          d0pi_a = d1pi(i0sidi)
          d0tt_a = d1tt(i0sidi)
          d0mm_a = d1mm(i0sidi)
          
          d0nvcc1_a =      (d0nvcm1_a-d0nvcm2_a)*d0ni_a
          d0nvcc2_a = 0.4D0*d0nvcm2_a*d0pi_a
          d0nvcc3_a =      (d0nvcm2_a-d0nvcm3_a)*d0ni_a
          d0nvcc4_a = 0.4D0*d0nvcm3_a*d0pi_a
          
          d0nvcc3_a = (d0tt_a/d0mm_a)*(2.5D0*d0nvcc1_a + d0nvcc3_a)
          d0nvcc4_a = (d0tt_a/d0mm_a)*(2.5D0*d0nvcc2_a + d0nvcc4_a)
          
          d1nvcc1(i0sidi) = d0nvcc1_a
          d1nvcc2(i0sidi) = d0nvcc2_a
          d1nvcc3(i0sidi) = d0nvcc3_a
          d1nvcc4(i0sidi) = d0nvcc4_a
      
       ENDDO

    ELSE
       
       DO i0sidi = 1, i0smax 
          d1nvcc1(i0sidi) = 0.D0
          d1nvcc2(i0sidi) = 0.D0
          d1nvcc3(i0sidi) = 0.D0
          d1nvcc4(i0sidi) = 0.D0
       ENDDO
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2CALV_PQ


  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF MASS SCALAR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  
  SUBROUTINE T2CALV_MS
    
    USE T2CNST, ONLY:&
         d0vci2

    USE T2COMM, ONLY:&
         i0vmax,i0smax,d3ms
    INTEGER(i0ikind)::&
         i0vidi,i0vidj,i0sidi
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       d3ms(i0vidi,i0vidj,i0midi) = 0.D0
    ENDDO
    ENDDO

    !C
    !C EQUATION FOR PSI'
    !C

    !C PSI'
    i0vidi = 1
    i0vidj = 1
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg
    
    !C
    !C EQUATION FOR I
    !C
    
    !C I
    i0vidi = 2
    i0vidj = 2
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg
    

    !C
    !C EQUATION FOR Et
    !C
    
    !C Et
    i0vidi = 3
    i0vidj = 3
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0vci2
    
    
    !C
    !C EQUATION FOR Ep
    !C
    
    !C Ep
    i0vidi = 4
    i0vidj = 4
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0vci2
        
    DO i0sidi = 1, i0smax
       
       !C
       !C EQUATION FOR N
       !C
       i0vidi = 8*i0sidi - 2
       
       !C N 
       i0vidj = i0vidi
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg
              
       !C
       !C EQUATION FOR Fb
       !C 
       i0vidi = 8*i0sidi

       !C Fb
       i0vidj = i0vidi
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0bb
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vidi = 8*i0sidi + 1
       
       !C Ft
       i0vidj = i0vidi
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg
       
       !C
       !C EQUATION FOR P
       !C
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = i0vidi
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*1.5D0
       
       !C
       !C EQUATION FOR Qb
       !C
       i0vidi = 8*i0sidi + 4
       
       !C Qb       
       i0vidj = i0vidi
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0bb
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vidi = 8*i0sidi + 5
       
       !C Qt       
       i0vidj = i0vidi
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_MS
  
  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ADVECTION VECTOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_AV
    
    USE T2CNST, ONLY:&
         d0vci2
    
    USE T2COMM, ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1mm,d1tt,d1nn,d1pp,d1ni,d1pi,&
         d1ub,d1ur,d1ut,d1up,d1u2,d1qr,d1qb,d1qt,d1qp,&
         i0vmax,i0dmax,d4av
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0vidi,i0vidj
    
    REAL(   i0rkind)::&
         d0mm_a,d0tt_a,d0nn_a,d0pp_a,d0ni_a,d0pi_a,&
         d0ur_a,d0up_a,d0ub_a,d0u2_a,&
         d0qr_a,d0qp_a,d0qt_a,d0qb_a

    REAL(   i0rkind)::&
         d0x1,d0x2
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0didi = 1, i0dmax
          d4av(i0didi,i0vidi,i0vidj,i0midi) = 0.D0
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C EQUATION FOR PSI'
    !C
    
    i0vidi = 1
    
    !C PSI'
    i0vidj = 1
    d4av(1,i0vidi,i0vidj,i0midi) = d0ugr
    d4av(2,i0vidi,i0vidj,i0midi) = d0ugp
    
    !C
    !C EQUATION FOR I
    !C

    i0vidi = 2
    
    !C I
    i0vidj = 2
    d4av(1,i0vidi,i0vidj,i0midi) = d0ugr
    d4av(2,i0vidi,i0vidj,i0midi) = d0ugp
    
    !C
    !C EQUATION FOR Et
    !C

    i0vidi = 3
    
    !C PSI'
    i0vidj = 1
    d4av(1,i0vidi,i0vidj,i0midi) = -d0mfcst*d0sqrtg*d0ctgrr/d0etcst
    d4av(2,i0vidi,i0vidj,i0midi) = -d0mfcst*d0sqrtg*d0ctgrp/d0etcst
    
    !C Et
    i0vidj = 3
    d4av(1,i0vidi,i0vidj,i0midi) =  d0ugr*d0vci2
    d4av(2,i0vidi,i0vidj,i0midi) =  d0ugp*d0vci2
    
    !C
    !C EQUATION FOR Ep
    !C

    i0vidi = 4

    !C Ep
    i0vidj = 4
    d4av(1,i0vidi,i0vidj,i0midi) =  d0ugr*d0vci2
    d4av(2,i0vidi,i0vidj,i0midi) =  d0ugp*d0vci2
    
    !C
    !C EQUATION FOR Er
    !C
    
    i0vidi = 5
    
    !C Ep
    i0vidj = 4
    d4av(1,i0vidi,i0vidj,i0midi) = d0epcst*d0sqrtg*d0ctgrp/d0ercst
    d4av(2,i0vidi,i0vidj,i0midi) = d0epcst*d0sqrtg*d0ctgpp/d0ercst
    
    !C Er
    i0vidj = 5
    d4av(1,i0vidi,i0vidj,i0midi) = d0sqrtg*d0ctgrr
    d4av(2,i0vidi,i0vidj,i0midi) = d0sqrtg*d0ctgrp
    
    DO i0sidi = 1, i0smax
       
       d0mm_a = d1mm(i0sidi)
       d0tt_a = d1tt(i0sidi)
       d0nn_a = d1nn(i0sidi)
       d0ni_a = d1ni(i0sidi)
       d0pi_a = d1pi(i0sidi)
       d0ur_a = d1ur(i0sidi)
       d0up_a = d1up(i0sidi)
       d0ub_a = d1ub(i0sidi)
       d0pp_a = d1pp(i0sidi)
       d0u2_a = d1u2(i0sidi)
       d0qr_a = d1qr(i0sidi)
       d0qp_a = d1qp(i0sidi)
       d0qb_a = d1qb(i0sidi)
       d0qt_a = d1qt(i0sidi)

       !C
       !C EQUATION FOR N
       !C
       
       i0vidi = 8*i0sidi - 2
       
       !C N
       i0vidj = i0vidi 
       d4av(1,i0vidi,i0vidj,i0midi) = d0ugr+d0sqrtg*d0ur_a
       d4av(2,i0vidi,i0vidj,i0midi) = d0ugp+d0sqrtg*d0up_a
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb
       i0vidj = i0vidi
       d4av(1,i0vidi,i0vidj,i0midi) = d0ugr*d0bb
       d4av(2,i0vidi,i0vidj,i0midi) = d0ugr*d0bb + d0sqrtg*d0ctbp*d0ub_a
       
       !C P
       i0vidj = i0vidi + 2
       d4av(2,i0vidi,i0vidj,i0midi) = d0ppcst*d0sqrtg*d0ctbp/d0mm_a/d0fbcst
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vidi = 8*i0sidi + 1
       
       !C Ft
       i0vidj = i0vidi 
       d4av(1,i0vidi,i0vidj,i0midi) = d0ugr + d0sqrtg*d0ur_a
       d4av(2,i0vidi,i0vidj,i0midi) = d0ugp + d0sqrtg*d0up_a
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       d0x1 = -0.5D0*d0sqrtg*d0mm_a*d0nn_a*d0u2_a*d0pi_a
       d0x2 =  d0sqrtg*d0pi_a
       
       !C P
       i0vidj = i0vidi
       d4av(1,i0vidi,i0vidj,i0midi)&
            = 1.5D0*d0ugr + d0x1*d0ur_a + d0x2*d0qr_a
       d4av(2,i0vidi,i0vidj,i0midi)&
            = 1.5D0*d0ugp + d0x1*d0up_a + d0x2*d0qp_a
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C Fb
       i0vidj = i0vidi - 4
       d4av(2,i0vidi,i0vidj,i0midi)&
            = d0fbcst*(d0qb_a - 1.5D0*d0pp_a*d0ub_a)&
            * d0sqrtg*d0ctbp*d0ni_a/d0qbcst
       
       !C P
       i0vidj = i0vidi - 2
       d4av(2,i0vidi,i0vidj,i0midi)&
            = d0ppcst*2.5D0*d0tt_a*d0sqrtg*d0ctbp/d0mm_a/d0qbcst
       
       !C Qb
       i0vidj = i0vidi
       d4av(1,i0vidi,i0vidj,i0midi) = d0ugr*d0bb
       d4av(2,i0vidi,i0vidj,i0midi) = d0ugp*d0bb+d0sqrtg*d0ctbp*d0ub_a
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5
       
       !C Ft
       i0vidj = i0vidi - 4
       d4av(1,i0vidi,i0vidj,i0midi)&
            = d0ftcst*d0sqrtg*d0ni_a*(d0qr_a - 1.5D0*d0pp_a*d0ur_a)/d0qtcst
       d4av(2,i0vidi,i0vidj,i0midi)&
            = d0ftcst*d0sqrtg*d0ni_a*(d0qp_a - 1.5D0*d0pp_a*d0up_a)/d0qtcst
       
       !C Qt
       i0vidj = i0vidj
       d4av(1,i0vidi,i0vidj,i0midi) = d0ugr + d0sqrtg*d0ur_a
       d4av(2,i0vidi,i0vidj,i0midi) = d0ugp + d0sqrtg*d0up_a
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_AV

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ADVECTION TENSOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C
  !C---------------------------------------------------------  
  SUBROUTINE T2CALV_AT
    
    USE T2CNST
    
    USE T2COMM, ONLY:&
         i0vmax,i0dmax,i0wmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1mm,d1tt,d1nn,d1pp,d1ni,d1pi,&
         d1ub,d1ur,d1ut,d1up,d1u2,d1qr,d1qb,d1qt,d1qp,&
         d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,d6at    
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0widi,i0vidi,i0vidj
    
    REAL(   i0rkind)::&
         d0mm_a,d0mi_a,d0ur_a,d0up_a,d0ub_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a

    REAL(   i0rkind)::&
         d0x1,d0x2,d0x3,d0x4,d0x5,&
         d0c1r_a,d0c1p_a,d0c2r_a,d0c2p_a

    !C
    !C INITIALIZATION
    !C

    DO i0didi = 1, i0dmax
    DO i0didj = 1, i0dmax
       DO i0widi = 1, i0wmax
          DO i0vidi = 1, i0vmax
          DO i0vidj = 1, i0vmax
             d6at(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
           
    d0x1 = d0sqrtg*d0bt2*d0bp2i*d0ctbp/d0bb
    d0x2 = 2.D0*d0ctbp/d0bb
    d0x3 = 2.D0*d0ctbp/d0cobt
    d0x4 = 3.D0*d0cobt*d0ctbp/(d0bb**3)
    d0x5 = 3.D0*d0ctbp/(d0bb**2)
    
    DO i0sidi = 1, i0smax
       
       d0mm_a = d1mm(i0sidi)
       d0ur_a = d1ur(i0sidi)
       d0up_a = d1up(i0sidi)
       d0ub_a = d1ub(i0sidi)
       
       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       d0nvcc3_a = d1nvcc3(i0sidi)
       d0nvcc4_a = d1nvcc4(i0sidi) 
      
       d0mi_a = 1.D0/d0mm_a
       
       !C
       !C EQUATION FOR Fb
       !C
       
       i0vidi = 8*i0sidi 
       
       !C Fb (B)
       i0vidj = i0vidi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0x2*d0nvcc1_a*d0mi_a
       
       !C Ft (B)
       i0vidj = i0vidi + 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x1*d0x3*d0nvcc1_a*d0mi_a/d0fbcst
       
       !C Qb (B)
       i0vidj = i0vidi + 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0x2*d0nvcc2_a*d0mi_a/d0fbcst
       
       !C Qt (B)
       i0vidj = i0vidi + 5
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0qtcst*d0x1*d0x3*d0nvcc2_a*d0mi_a/d0fbcst
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vidi = 8*i0sidi + 1
       
       !C Fb (B)
       i0vidj = i0vidi - 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x1*d0x4*d0nvcc1_a*d0mi_a/d0ftcst

       !C Ft (B)
       i0vidj = i0vidi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0x5*d0nvcc1_a*d0mi_a
       
       !C Qb (B)
       i0vidj = i0vidi + 3
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0x4*d0nvcc2_a*d0mi_a/d0ftcst

       !C Qt (B)
       i0vidj = i0vidi + 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0qtcst*d0x1*d0x5*d0nvcc2_a*d0mi_a/d0ftcst
       
       !C
       !C EQUATION FOR P
       !C

       i0vidi = 8*i0sidi + 2
       
       d0c1r_a =  d0ur_a/d0bb
       d0c1p_a = (d0up_a - 3.D0*d0ctbp*d0ub_a/d0bb)/d0bb
       d0c2r_a =  d0ur_a/d0cobt
       d0c2p_a = (d0up_a - 3.D0*d0ctbp*d0ub_a/d0bb)/d0cobt

       !C Fb (B)
       i0vidj = i0vidi - 2
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x1*d0c1r_a*d0nvcc1_a/d0ppcst
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x1*d0c1p_a*d0nvcc1_a/d0ppcst


       !C Ft
       i0vidj = i0vidi - 1
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x1*d0c2r_a*d0nvcc1_a/d0ppcst
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x2*d0c2p_a*d0nvcc1_a /d0ppcst
       
       !C Qb
       i0vidj = i0vidi + 2
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0c1r_a*d0nvcc2_a/d0ppcst
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0c1p_a*d0nvcc2_a/d0ppcst
       
       !C Qt
       i0vidj = i0vidi + 3
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            = d0qtcst*d0x1*d0c2r_a*d0nvcc2_a/d0ppcst
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = d0qtcst*d0x1*d0c2p_a*d0nvcc2_a/d0ppcst
       
       !C
       !C EQUATION FOR Qb
       !C

       i0vidi = 8*i0sidi + 4
       
       !C Fb
       i0vidj = i0vidi - 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x1*d0x2*d0nvcc3_a/d0qbcst
       
       !C Ft
       i0vidj = i0vidi - 3
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x1*d0x3*d0nvcc3_a/d0qbcst
       
       !C Qb
       i0vidj = i0vidi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0x2*d0nvcc4_a/d0qbcst
       
       !C Qt
       i0vidj = i0vidi + 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0qtcst*d0x1*d0x3*d0nvcc4_a/d0qbcst
       
       !C
       !C EQUATION FOR Qt
       !C

       i0vidi = 8*i0sidi + 5
       
       !C Fb
       i0vidj = i0vidi - 5
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x1*d0x4*d0nvcc3_a/d0qtcst
       
       !C Ft
       i0vidj = i0vidi - 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x1*d0x5*d0nvcc3_a/d0qtcst
       
       !C Qb
       i0vidj = i0vidi - 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0x4*d0nvcc4_a/d0qtcst
       
       !C Qt
       i0vidj = i0vidi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qtcst*d0x1*d0x5*d0nvcc4_a/d0qtcst
       
    ENDDO
    
    RETURN

  END SUBROUTINE T2CALV_AT

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF DIFFUSION TENSOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------  
  SUBROUTINE T2CALV_DT
    
    USE T2CNST
    USE T2COMM,ONLY:&
         i0vmax,i0dmax,i0vmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1mm,d1pp,d1ur,d1up,d1ub,d1qb,d1vb,&
         d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         d5dt
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0vidi,i0vidj
    
    REAL(   i0rkind)::&
         d0mm_a,d0mi_a,d0pp_a,d0ur_a,d0up_a,d0ub_a,d0qb_a,d0vb_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    REAL(   i0rkind)::&
         d0x1,d0x2,d0x3,d0c1r_a,d0c1p_a

    !C
    !C INITIALIZATION
    !C 
    DO i0didi = 1, i0dmax
    DO i0didj = 1, i0dmax
       DO i0vidi = 1, i0vmax
       DO i0vidj = 1, i0vmax
          d5dt(i0didi,i0didj,i0vidi,i0vidj,i0midi) = 0.D0 
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    d0x1 = 2.D0*d0sqrtg*(d0ctbp**2)/d0bb
    d0x2 = 3.D0*d0sqrtg*(d0ctbp**2)*d0cobt/(d0bb**3)
    d0x3 = d0sqrtg*d0ctbp/d0bb
    
    DO i0sidi = 1, i0smax
       
       d0mm_a     = d1mm(i0sidi)
       d0mi_a     = 1.D0/d0mm_a
       d0pp_a     = d1pp(i0sidi)
       d0ur_a     = d1ur(i0sidi)
       d0up_a     = d1up(i0sidi)
       d0ub_a     = d1ub(i0sidi)
       d0qb_a     = d1qb(i0sidi)
       d0vb_a     = d1vb(i0sidi)
       d0nvcc1_a  = d1nvcc1(i0sidi)
       d0nvcc2_a  = d1nvcc2(i0sidi)
       d0nvcc3_a  = d1nvcc3(i0sidi)
       d0nvcc4_a  = d1nvcc4(i0sidi)
       
       !C   
       !C EQUATION FOR Fb
       !C
       
       i0vidi = 8*i0sidi 
       
       !C N
       i0vidj = i0vidi - 2
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0nncst*d0x1*d0nvcc1_a*d0mi_a*d0ub_a/d0fbcst
       !C Fb
       i0vidj = i0vidi
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0fbcst*d0x1*d0nvcc1_a*d0mi_a/d0fbcst
       !C P
       i0vidj = i0vidi + 2
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0ppcst*d0x1*d0nvcc2_a*d0mi_a*d0vb_a/d0fbcst
       !C Qb
       i0vidj = i0vidi + 4
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0qbcst*d0x1*d0nvcc2_a*d0mi_a/d0fbcst
       
       !C
       !C  EQUATION FOR Ft
       !C 

       i0vidi = 8*i0sidi + 1
       
       !C N
       i0vidj = i0vidi - 3
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0nncst*d0x2*d0nvcc1_a*d0mi_a*d0ub_a/d0ftcst
       !C Fb
       i0vidj = i0vidi - 1
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0fbcst*d0x2*d0nvcc1_a*d0mi_a/d0ftcst
       !C P
       i0vidj = i0vidi + 1
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0ppcst*d0x2*d0nvcc2_a*d0mi_a*d0vb_a/d0ftcst
       !C Qb
       i0vidj = i0vidi + 3
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0qbcst*d0x2*d0nvcc2_a*d0mi_a/d0ftcst

       !C
       !C EQUATION FOR P
       !C

       i0vidi = 8*i0sidi + 2
       
       d0c1r_a = d0ur_a
       d0c1p_a = d0up_a - 3.D0*d0ctbp*d0ub_a/d0bb
       
       !C N
       i0vidj = i0vidi - 4
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            = -d0nncst*d0x3*d0c1r_a*d0nvcc1_a*d0ub_a/d0ppcst
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0nncst*d0x3*d0c1p_a*d0nvcc1_a*d0ub_a/d0ppcst
       !C Fb
       i0vidj = i0vidi - 2
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            =  d0fbcst*d0x3*d0c1r_a*d0nvcc1_a/d0ppcst
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0fbcst*d0x3*d0c1p_a*d0nvcc1_a/d0ppcst
       !C P
       i0vidj = i0vidi
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            = -d0ppcst*d0x3*d0c1r_a*d0nvcc2_a*d0vb_a/d0ppcst
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0ppcst*d0x3*d0c1p_a*d0nvcc2_a*d0vb_a/d0ppcst
       !C Qb
       i0vidj = i0vidi + 2
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            =  d0qbcst*d0x3*d0c1r_a*d0nvcc2_a/d0ppcst
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0qbcst*d0x3*d0c1p_a*d0nvcc2_a/d0ppcst

       !C   
       !C EQUATION FOR Qb
       !C

       i0vidi = 8*i0sidi + 4
       
       !C N
       i0vidj = i0vidi - 6
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0nncst*d0x1*d0nvcc3_a*d0ub_a/d0qbcst
       !C Fb
       i0vidj = i0vidi - 4
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0fbcst*d0x1*d0nvcc3_a/d0qbcst
       !C P
       i0vidj = i0vidi - 2
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0ppcst*d0x1*d0nvcc4_a*d0vb_a/d0qbcst
       !C Qb
       i0vidj = i0vidi
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0qbcst*d0x1*d0nvcc4_a/d0qbcst
       
       !C
       !C  EQUATION FOR Qt
       !C 

       i0vidi = 8*i0sidi + 5
       
       !C N
       i0vidj = i0vidi - 7
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0nncst*d0x2*d0nvcc3_a*d0ub_a/d0qtcst
       !C Fb
       i0vidj = i0vidi - 5
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0fbcst*d0x2*d0nvcc3_a/d0qtcst
       !C P
       i0vidj = i0vidi - 3
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0ppcst*d0x2*d0nvcc4_a*d0vb_a/d0qtcst
       !C Qb
       i0vidj = i0vidi - 1
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0qbcst*d0x2*d0nvcc4_a/d0qtcst
       
    ENDDO

    RETURN
    
  END SUBROUTINE T2CALV_DT

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF GRADIENT VECTOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_GV

    USE T2CNST
    USE T2COMM, ONLY:&
         i0vmax,i0dmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1pp,d1tt,d1ur,d1up,d4gv
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0vidi,i0vidj
    REAL(   i0rkind)::&
         d0tt_a,d0ur_a,d0up_a,d0c1_a
    
    !C
    !C INITIALIZATION
    !C
    DO i0didi = 1, i0dmax
       DO i0vidi = 1, i0vmax
       DO i0vidj = 1, i0vmax
          d4gv(i0didi,i0vidi,i0vidj,i0midi) = 0.D0
       ENDDO
       ENDDO
    ENDDO

    !C
    !C EQUATION FOR PSI
    !C

    i0vidi = 1
    
    !C Et
    i0vidj = 3
    d4gv(1,i0vidi,i0vidj,i0midi) = -d0etcst*d0sqrtg/d0mfcst
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vidi = 2
    
    !C Ep
    i0vidj = 4 
    d4gv(1,i0vidi,i0vidj,i0midi) =  d0epcst*d0cogtt/d0btcst
    
    !C Er
    i0vidj = 5
    d4gv(2,i0vidi,i0vidj,i0midi) = -d0ercst*d0cogtt/d0btcst
    
    
    !C 
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    d4gv(1,i0vidi,i0vidj,i0midi) =  d0btcst*d0cogpp/d0epcst
    
    DO i0sidi = 1, i0smax
       
       d0tt_a = d1tt(i0sidi)
       d0ur_a = d1ur(i0sidi)
       d0up_a = d1up(i0sidi)
       
       !C
       !C EQUATION FOR Fr
       !C

       i0vidi = 8*i0sidi - 1
       
       !C P
       i0vidj = i0vidi + 3
       d4gv(1,i0vidi,i0vidj,i0midi) =  d0ppcst*d0sqrtg*d0ctgrr/d0frcst
       d4gv(2,i0vidi,i0vidj,i0midi) =  d0ppcst*d0sqrtg*d0ctgrp/d0frcst
       
       
       !C 
       !C EQUATION FOR P
       !C

       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = i0vidi
       d4gv(1,i0vidi,i0vidj,i0midi) = -d0ppcst*d0sqrtg*d0ur_a/d0ppcst
       d4gv(2,i0vidi,i0vidj,i0midi) = -d0ppcst*d0sqrtg*d0up_a/d0ppcst
       
       !C
       !C EQUATION FOR Qr
       !C

       i0vidi = 8*i0sidi + 3
       
       d0c1_a = 2.5D0*d0sqrtg*d0tt_a
       
       !C N
       i0vidj = i0vidi - 5
       d4gv(1,i0vidi,i0vidj,i0midi) = -d0nncst*d0c1_a*d0ctgrr*d0tt_a/d0qrcst
       d4gv(2,i0vidi,i0vidj,i0midi) = -d0nncst*d0c1_a*d0ctgrp*d0tt_a/d0qrcst
       
       !C P
       i0vidj = i0vidi - 1
       d4gv(1,i0vidi,i0vidj,i0midi) =  d0ppcst*d0c1_a*d0ctgrr*2.D0/d0qrcst
       d4gv(2,i0vidi,i0vidj,i0midi) =  d0ppcst*d0c1_a*d0ctgrp*2.D0/d0qrcst
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_GV

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF GRADIENT TENSOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_GT
    
    USE T2CNST
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i0wmax,i0dmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1mm, d1ut,d1ub,d1vb,&
         d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         d6gt 
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0widi,i0vidi,i0vidj
    
    REAL(   i0rkind)::&
         d0mm_a,d0mi_a,d0ub_a,d0ut_a,d0vb_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    REAL(   i0rkind)::&
         d0x1,d0x2,d0c1_a
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0didi = 1, i0dmax
    DO i0didj = 1, i0dmax
       DO i0widi = 1, i0wmax
          DO i0vidi = 1, i0vmax
          DO i0vidj = 1, i0vmax
             d6gt(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    d0x1 = 3.D0*d0sqrtg*(d0ctbp**2)/(d0bb**2)
    d0x2 = d0x1*d0bt2*d0bp2i
    
    DO i0sidi = 1, i0smax
       
       d0mm_a    = d1mm(   i0sidi)
       d0mi_a    = d0mm_a
       d0ut_a    = d1ut(   i0sidi)
       d0ub_a    = d1ub(   i0sidi)
       d0vb_a    = d1vb(   i0sidi)
       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       d0nvcc3_a = d1nvcc3(i0sidi)
       d0nvcc4_a = d1nvcc4(i0sidi)
       
       d0c1_a = d0ub_a/d0bb - d0ut_a/d0cobt
       
       !C
       !C EQUATION FOR Fb
       !C
       
       i0vidi = 8*i0sidi
       
       !C N (B)
       i0vidj = i0vidi - 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0nncst*d0x1*d0nvcc1_a*d0mi_a*d0ub_a/d0fbcst
       
       !C Fb (B)
       i0vidj = i0vidi
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0fbcst*d0x1*d0nvcc1_a*d0mi_a/d0fbcst
       !C P (B)
       i0vidj = i0vidi + 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0ppcst*d0x1*d0nvcc2_a*d0mi_a*d0vb_a/d0fbcst
       !C Qb
       i0vidj = i0vidi + 4
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0qbcst*d0x1*d0nvcc2_a*d0mi_a/d0fbcst
       
       !C
       !C EQUATION FOR P
       !C

       i0vidi = 8*i0sidi + 2

       !C N (B)
       i0vidj = i0vidi - 4
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0nncst*d0x2*d0nvcc1_a*d0ub_a*d0c1_a/d0ppcst

       !C N (Ub)
       i0vidj = i0vidi - 4
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0nncst*d0x1*d0nvcc1_a*d0ub_a/d0ppcst
       
       !C Fb (B)
       i0vidj = i0vidi - 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0fbcst*d0x2*d0nvcc1_a*d0c1_a/d0ppcst

       !C Fb (Ub)
       i0vidj = i0vidi - 2
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x1*d0nvcc1_a/d0ppcst
       

       !C P (B)
       i0vidj = i0vidi 
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0ppcst*d0x2*d0nvcc2_a*d0c1_a*d0vb_a/d0ppcst

       !C P (Ub)
       i0vidj = i0vidi
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0ppcst*d0x1*d0nvcc2_a*d0vb_a/d0ppcst
       
       !C Qb (B)
       i0vidj = i0vidi + 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0qbcst*d0x2*d0nvcc2_a*d0c1_a/d0ppcst

       !C Qb (Ub)
       i0vidj = i0vidi + 2
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x1*d0nvcc2_a/d0ppcst
       
       !C
       !C EQUATION FOR Qb
       !C
     
       i0vidi = 8*i0sidi + 4
       
       !C N  (B)
       i0vidj = i0vidi - 6
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0nncst*d0x1*d0nvcc3_a*d0ub_a/d0qbcst
       
       !C Fb (B)
       i0vidj = i0vidi - 4
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0fbcst*d0x1*d0nvcc3_a/d0qbcst
       
       !C P  (B)
       i0vidj = i0vidi - 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0ppcst*d0x1*d0nvcc4_a*d0vb_a/d0qbcst

       !C Qb (B)
       i0vidj = i0vidi
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0qbcst*d0x1*d0nvcc4_a/d0qbcst
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_GT

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ECITATION SCALAR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_ES
    
    USE T2CNST,ONLY:&
         d0eps0,d0rmu0
    USE T2COMM,ONLY:&
         i0vmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1ee, d1mm, d1nn, d1pp,d1tt,&
         d1hex,d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4,&
         d3es
    
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj
    
    REAL(   i0rkind)::&
         d0mm_a,d0mi_a,d0ee_a,d0ee_b,d0nn_a,d0pp_a,d0tt_a,&
         d0nfcf1_ab,d0nfcf2_ab,d0nfcf3_ab,d0nfcf4_ab,d0hex_a
    REAL(   i0rkind)::&
         d0x1,&
         d0c07a_1,d0c07a_2,d0c07a_3,d0c07a_4,&
         d0c08a_1,d0c08a_2,d0c08a_3,&
         d0c09a_1,d0c09a_2,d0c1_b

    !C
    !C INITIALIZATION
    !C

    DO i0vidi = 1, i0vmax
    DO i0vidj = 1, i0vmax
       d3es(i0vidi,i0vidj,i0midi) = 0.D0
    ENDDO
    ENDDO
    
    d0x1 = (d0sqrtg**2)*d0bb*d0ctbpi
    
    !C
    !C EQUATION FOR Et
    !C
    
    i0vidi = 3
    
    DO i0sidj = 1, i0smax
       
       d0ee_b  = d1ee(i0sidj)
       
       !C Ft
       i0vidj = 8*i0sidj + 1
       d3es(i0vidi,i0vidj,i0midi) = d0ftcst*d0sqrtg*d0rmu0*d0ee_b/d0etcst
       
    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    DO i0sidj = 1, i0smax
       
       d0ee_b = d1ee(i0sidj)
       d0c1_b = d0sqrtg*d0rmu0*d0ctbpi*d0ee_b

       !C Fb
       i0vidj = 8*i0sidj 
       d3es(i0vidi,i0vidj,i0midi) =  d0fbcst*d0c1_b*d0bb/d0epcst

       !C Ft
       i0vidj = 8*i0sidj + 1
       d3es(i0vidi,i0vidj,i0midi) = -d0ftcst*d0c1_b*d0ctbt/d0epcst
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C
    
    i0vidi = 5
    
    DO i0sidj = 1, i0smax
       
       d0ee_b  = d1ee(i0sidj)
       
       !C N
       i0vidj = 8*i0sidj - 2
       d3es(i0vidi,i0vidj,i0midi) = -d0nncst*d0sqrtg*d0ee_b/d0eps0/d0ercst
       
    ENDDO
    
    DO i0sidi = 1, i0smax
       
       d0nn_a  = d1nn( i0sidi)
       d0pp_a  = d1pp( i0sidi)
       d0tt_a  = d1tt( i0sidi)
       d0mm_a  = d1mm( i0sidi)
       d0ee_a  = d1ee( i0sidi)
       d0hex_a = d1hex(i0sidi)
       d0mi_a  = d0mm_a

       d0c07a_1 = d0sqrtg*d0ee_a*d0nn_a*d0ctgrp
       d0c07a_2 = d0sqrtg*d0ee_a*d0nn_a*d0ctgrr
       d0c07a_3 = d0x1*d0ee_a*d0cobt
       d0c07a_4 = d0x1*d0ee_a*d0bb
       
       d0c08a_1 = d0sqrtg*d0ee_a*d0nn_a*d0mi_a*d0ctbt
       d0c08a_2 = d0sqrtg*d0ee_a*d0nn_a*d0mi_a*d0ctbp
       d0c08a_3 = d0sqrtg*d0bb
       
       d0c09a_1 = d0sqrtg*d0ee_a*d0nn_a*d0mi_a
       d0c09a_2 = (d0sqrtg**2)*d0ee_a*d0ctbp*d0mi_a

       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi = 8*i0sidi - 1
       
       !C Ep C_07a.04
       i0vidj = 4
       d3es(i0vidi,i0vidj,i0midi) = -d0epcst*d0c07a_1/d0frcst

       
       !C Er C_07a.05
       i0vidj = 5
       d3es(i0vidi,i0vidj,i0midi) = -d0ercst*d0c07a_2/d0frcst


       !C Fb C_07a.08a
       i0vidj = i0vidi + 1
       d3es(i0vidi,i0vidj,i0midi) = -d0fbcst*d0c07a_3/d0frcst
       
       !C Ft C_07a.09a
       i0vidj = i0vidi + 2
       d3es(i0vidi,i0vidj,i0midi) =  d0ftcst*d0c07a_4/d0frcst
       
       !C
       !C EQUATION FOR Fb
       !C
       
       i0vidi = 8*i0sidi
       
       DO i0sidj = 1, i0smax
          
          d0nfcf1_ab = d2nfcf1(i0sidi,i0sidj)
          d0nfcf2_ab = d2nfcf2(i0sidi,i0sidj)
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0etcst*d0c08a_1/d0fbcst
             
             !C Ep
             i0vidj = 4
             d3es(i0vidi,i0vidj,i0midi) = -d0epcst*d0c08a_2/d0fbcst
             
          ENDIF
          
          !C Fb
          i0vidj = 8*i0sidj
          d3es(i0vidi,i0vidj,i0midi) &
               = -d0fbcst*d0c08a_3*d0nfcf1_ab*d0mi_a/d0fbcst
          
          !C Qb
          i0vidj = 8*i0sidj + 4
          d3es(i0vidi,i0vidj,i0midi) &
               = -d0qbcst*d0c08a_3*d0nfcf2_ab*d0mi_a/d0fbcst
          
       ENDDO
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vidi = 8*i0sidi + 1
       
       DO i0sidj = 1, i0smax
          
          d0nfcf1_ab = d2nfcf1(i0sidi,i0sidj)
          d0nfcf2_ab = d2nfcf2(i0sidi,i0sidj)
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0etcst*d0c09a_1/d0ftcst
             
             !C Fr
             i0vidj = 8*i0sidj - 1
             d3es(i0vidi,i0vidj,i0midi) = -d0frcst*d0c09a_2/d0ftcst
             
          ENDIF
          
          !C Ft
          i0vidj = 8*i0sidj + 1
          d3es(i0vidi,i0vidj,i0midi) &
               = -d0ftcst*d0sqrtg*d0nfcf1_ab*d0mi_a/d0ftcst
          
          !C Qt
          i0vidj = 8*i0sidj + 5
          d3es(i0vidi,i0vidj,i0midi) &
               = -d0qtcst*d0sqrtg*d0nfcf2_ab*d0mi_a/d0ftcst
          
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = i0vidi 
       d3es(i0vidi,i0vidj,i0midi) = d0ppcst*d0sqrtg*d0hex_a/d0ppcst
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 8*i0sidi + 3
       
       !C Ep
       i0vidj = 4
       d3es(i0vidi,i0vidj,i0midi) = -d0epcst*2.5D0*d0tt_a*d0c07a_1/d0qrcst
       
       !C Er
       i0vidj = 5
       d3es(i0vidi,i0vidj,i0midi) = -d0ercst*2.5D0*d0tt_a*d0c07a_2/d0qrcst
       
       !C Qb
       i0vidj = i0vidi + 1
       d3es(i0vidi,i0vidj,i0midi) = -d0qbcst*d0c07a_3/d0qrcst
       
       !C Qt
       i0vidj = i0vidi + 2
       d3es(i0vidi,i0vidj,i0midi) =  d0qtcst*d0c07a_4/d0qrcst
       
       !C
       !C EQUATION FOR Qb
       !C

       i0vidi = 8*i0sidi + 4 
       
       DO i0sidj = 1, i0smax
          
          d0nfcf3_ab = d2nfcf3(i0sidi,i0sidj)
          d0nfcf4_ab = d2nfcf4(i0sidi,i0sidj)
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0etcst*2.5D0*d0tt_a*d0c08a_1/d0qbcst
             
             !C Ep
             i0vidj = 4
             d3es(i0vidi,i0vidj,i0midi) = -d0epcst*2.5D0*d0tt_a*d0c08a_2/d0qbcst

          ENDIF
          
          !C Fb  
          i0vidj = 8*i0sidj 
          d3es(i0vidi,i0vidj,i0midi) = -d0fbcst*d0c08a_3*d0nfcf3_ab/d0qbcst
          
          !C Qb
          i0vidj = 8*i0sidj + 4
          d3es(i0vidi,i0vidj,i0midi) = -d0qbcst*d0c08a_3*d0nfcf4_ab/d0qbcst
          
       ENDDO
       
       !C
       !C EQUATION FOR Qt
       !C

       i0vidi = 8*i0sidi + 5
       
       DO i0sidj = 1, i0smax
          
          d0nfcf3_ab = d2nfcf3(i0sidi,i0sidj)
          d0nfcf4_ab = d2nfcf4(i0sidi,i0sidj)
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi)&
                  = -d0etcst*2.5D0*d0tt_a*d0c09a_1/d0qtcst
          ENDIF
          
          !C Ft
          i0vidj = 8*i0sidj + 1
          d3es(i0vidi,i0vidj,i0midi) = -d0ftcst*d0sqrtg*d0nfcf3_ab/d0qtcst
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Qr
             i0vidj = 8*i0sidj + 3
             d3es(i0vidi,i0vidj,i0midi) = -d0qrcst*d0c09a_2/d0qtcst

          ENDIF

          !C Qt
          i0vidj = 8*i0sidj + 5
          d3es(i0vidi,i0vidj,i0midi) = -d0qtcst*d0sqrtg*d0nfcf4_ab/d0qtcst
          
       ENDDO
       
    ENDDO

    RETURN
    
  END SUBROUTINE T2CALV_ES

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ECITATION VECTOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_EV
    
    USE T2CNST
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i0dmax,i0wmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1ee,d1mm,d1nn,d1ut,d1ub,d1pp,d1qb,d1qt,d1ni,&
         d1nvcc1,d1nvcc2,&
         d5ev
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0vidj,i0widi,i0didi
    REAL(   i0rkind)::&
         d0ee_a,d0mm_a,d0mi_a,d0nn_a,d0ni_a,d0ut_a,d0ub_a,d0pp_a,&
         d0qb_a,d0qt_a,&
         d0nvcc1_a,d0nvcc2_a
    REAL(   i0rkind)::&
         d0x1,d0x2,d0x3,d0x4,d0x5,d0x6,d0x7
    REAL(   i0rkind)::&
         d0c1_a,d0c2_a,d0c3_a
    
    !C 
    !C INITIALIZATION
    !C
    DO i0didi = 1, i0dmax   
       DO i0widi = 1, i0wmax 
          DO i0vidi = 1, i0vmax
          DO i0vidj = 1, i0vmax
             d5ev(i0didi,i0widi,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    d0x1 = 2.D0*d0sqrtg*sqrt(d0ctgtt)
    d0x2 = d0sqrtg*d0ctbp/d0bb
    d0x3 = d0bt2*d0bp2i
    d0x4 = 3.D0*d0sqrtg*(d0ctbp/d0bb)*d0ctbt
    d0x5 = 3.D0*d0sqrtg*(d0ctbp/d0bb)*d0ctbp   
    d0x6 = d0sqrtg*(3.D0*(d0bt2/d0bb2) - 1.D0)*(d0ctbp/d0bb)
    d0x7 = 3.D0*d0sqrtg*d0cobt*(d0ctbp**2)/(d0bb**3)
    
    !C
    !C EQUATION FOR Et
    !C

    i0vidi = 3
 
    !C PSI' (R)
    i0vidj = 1
    i0widi = 2
    d5ev(1,i0widi,i0vidi,i0vidj,i0midi) = d0mfcst*d0x1*d0ctgrr/d0etcst
    d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = d0mfcst*d0x1*d0ctgrp/d0etcst
    
    DO i0sidi = 1, i0smax
       
       d0ee_a    = d1ee(   i0sidi)
       d0mm_a    = d1mm(   i0sidi)
       d0nn_a    = d1nn(   i0sidi)
       d0ni_a    = d1ni(   i0sidi)
       d0ut_a    = d1ut(   i0sidi)
       d0ub_a    = d1ub(   i0sidi)
       d0pp_a    = d1pp(   i0sidi)
       d0qb_a    = d1qb(   i0sidi)
       d0qt_a    = d1qt(   i0sidi)
       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)

       d0mi_a = 1.D0/d0mm_a       
       
       d0c1_a = d0ee_a*d0mi_a&
              *(d0nvcc1_a*d0nn_a*(d0ub_a/d0bb - d0ut_a/d0cobt)&
              + d0nvcc2_a*       (d0qb_a/d0bb - d0qt_a/d0cobt))
       d0c2_a = d0ee_a*d0mi_a*d0nvcc1_a*d0nn_a
       d0c3_a = d0ee_a*d0mi_a*d0nvcc2_a*d0pp_a
       
       !C
       !C EQUATION FOR Fb
       !C
       
       i0vidi = 8*i0sidi
       
       !C Fb (B)
       
       i0vidj = i0vidi
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = -d0fbcst*d0x2*d0ub_a/d0fbcst
       
       !C
       !C EQUATION FOR Qb
       !C

       i0vidi = 8*i0sidi + 4
       
       !C Et (B)
       i0vidj = 3
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = -d0etcst*d0x3*d0x4*d0c1_a/d0qbcst
       
       !C Et (Ub)
       i0vidj = 3
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0etcst     *d0x4*d0c2_a/d0qbcst
              
       !C Et (Qb/P)
       i0vidj = 3
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0etcst     *d0x4*d0c3_a/d0qbcst
       
       !C Ep (B)
       i0vidj = 4
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = -d0epcst*d0x3*d0x5*d0c1_a/d0qbcst
       
       !C Ep (Ub)
       i0vidj = 4
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0epcst     *d0x5*d0c2_a/d0qbcst
       
       !C Ep (Qb/P)
       i0vidj = 4
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0epcst     *d0x5*d0c3_a/d0qbcst
       
       !C Fb (B)
       i0vidj = 8*i0sidi
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x2*d0ni_a*(d0qb_a - 1.5D0*d0pp_a*d0ub_a)/d0qbcst
       
       !C Qb (B) 
       i0vidj = 8*i0sidi + 4
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = -d0qbcst*d0x2*d0ub_a/d0qbcst
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5
       
       !C Et (B)
       i0vidj = 3
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = -d0etcst*d0x3*d0x6*d0c1_a/d0qtcst
       
       !C Et (Ub)
       i0vidj = 3
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0etcst     *d0x6*d0c2_a/d0qtcst
       
       !C Et (Qb/P)
       i0vidj = 3
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0etcst     *d0x6*d0c3_a/d0qtcst

       !C Ep (B)
       i0vidj = 4
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) = -d0epcst*d0x3*d0x7*d0c1_a/d0qtcst
              
       !C Ep (Ub)
       i0vidj = 4
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0epcst     *d0x7*d0c2_a/d0qtcst
       
       !C Ep (Qb)
       i0vidj = 4
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi) =  d0epcst     *d0x7*d0c3_a/d0qtcst
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_EV
  
  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ECITATION VECTOR COEFFICIENTS
  !C
  !C             2014-01-26 H.SETO
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_ET
    
    USE T2CNST
    USE T2COMM,ONLY:&
         i0smax,i0dmax,i0vmax,i0wmax,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1mm,d1ut,d1ub,d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         d7et
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0widi,i0widj,i0vidi,i0vidj
    
    REAL(   i0rkind)::&
         d0mm_a,d0mi_a,d0ut_a,d0ub_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    REAL(   i0rkind)::&
         d0x1,d0x2,d0x3,d0x4,d0x5,d0c1_a
  
    !C
    !C
    !C
    DO i0didi = 1, i0dmax
    DO i0didj = 1, i0dmax
       DO i0widi = 1, i0wmax
       DO i0widj = 1, i0wmax
          DO i0vidi = 1, i0vmax
          DO i0vidj = 1, i0vmax
             d7et(i0didi,i0didj,i0widi,i0didj,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    d0x1 = 3.D0*d0sqrtg*(d0ctbp**2)*d0bt2*d0bp2i/d0bb2
    d0x2 = d0x1/d0bb
    d0x3 = d0x1/d0cobt
    d0x4 = d0x2*d0bt2*d0bp2i
    d0x5 = d0x3*d0bt2*d0bp2i

    DO i0sidi = 1, i0smax
       
       d0mm_a    = d1mm(   i0sidi)
       d0ub_a    = d1ub(   i0sidi)
       d0ut_a    = d1ut(   i0sidi)
       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       d0nvcc3_a = d1nvcc3(i0sidi)
       d0nvcc4_a = d1nvcc4(i0sidi)
       
       d0mi_a = 1.D0/d0mm_a
       d0c1_a = d0ub_a/d0bb - d0ut_a/d0cobt

       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x2*d0nvcc1_a*d0mi_a/d0fbcst
       
       !C Ft (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi + 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x3*d0nvcc1_a*d0mi_a/d0fbcst

       !C Qb (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi + 4
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x2*d0nvcc2_a*d0mi_a/d0fbcst

       !C Qt (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi + 5
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0qtcst*d0x3*d0nvcc2_a*d0mi_a/d0fbcst

       !C
       !C EQUATION FOR P
       !C

       i0vidi = 8*i0sidi + 2

       !C Fb (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi - 2
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x4*d0nvcc1_a*d0c1_a/d0ppcst

       !C Fb (Ub,B)
       i0widi = 2*i0sidi + 1
       i0widj = 1
       i0vidj = i0vidi - 2
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0fbcst*d0x2*d0nvcc1_a/d0ppcst

       !C Ft (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi - 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x5*d0nvcc1_a*d0c1_a/d0ppcst

       !C Ft (Ub,B)
       i0widi = 2*i0sidi + 1
       i0widj = 1
       i0vidj = i0vidi - 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0ftcst*d0x3*d0nvcc1_a/d0ppcst

       !C Qb (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi + 2
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x4*d0nvcc2_a*d0c1_a/d0ppcst

       !C Qb (ub,B)
       i0widi = 2*i0sidi + 1
       i0widj = 1
       i0vidj = i0vidi + 2
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0qbcst*d0x2*d0nvcc2_a/d0ppcst

       !C Qt (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi + 3
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0qtcst*d0x5*d0nvcc2_a*d0c1_a/d0ppcst

       !C Qt (ub,B)
       i0widi = 2*i0sidi + 1
       i0widj = 1
       i0vidj = i0vidi + 3
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0qtcst*d0x3*d0nvcc2_a/d0ppcst

       !C
       !C EQUATION FOR Qb
       !C

       i0vidi = 8*i0sidi + 4

       !C Fb (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi - 4
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0fbcst*d0x2*d0nvcc3_a/d0qbcst

       !C Ft (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi - 3
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0ftcst*d0x3*d0nvcc3_a/d0qbcst

       !C Qb (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0qbcst*d0x2*d0nvcc4_a/d0qbcst

       !C Qt (B,B)
       i0widi = 1
       i0widj = 1
       i0vidj = i0vidi + 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0qtcst*d0x3*d0nvcc4_a/d0qbcst
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_ET

  SUBROUTINE T2CALV_SS
    
    USE T2COMM, ONLY: i0vmax,d3ss
    
!    USE T2COMM, ONLY:&
!         i0vmax,i0smax,d3ss,&
!         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
!         d0nncst,d0frcst,d0fbcst,d0ftcst,&
!         d0ppcst,d0qrcst,d0qbcst,d0qtcst
    
    INTEGER(i0ikind)::&
         i0vidi,i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       d3ss(i0vidi,i0vidj,i0midi) = 0.D0
    ENDDO
    ENDDO

    RETURN
  END SUBROUTINE T2CALV_SS
  
  FUNCTION fd0k11ps(d0xa2)
 
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:i0xa,d2x,d2y,d2z,d2bcf    
    INTEGER(i0ikind)::i0xb
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k11ps
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2) ! V/Vt_a
    
    d0nut = 0.D0
    
    DO i0xb  = 1, i0smax

       d0x   = d2x(  i0xb,i0xa)  ! Vt_a/Vt_b
       d0y   = d2y(  i0xb,i0xa)  ! m_b/m_a
       d0z   = d2z(  i0xa,i0xb)  ! T_a/T_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x          ! V/Vt_b
       d0nud =      d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nus = 2.D0*d0bcf*d0z*(1.D0+d0y)*func_G(d0xb)/d0xa 
       d0nub = 2.D0*d0bcf*func_G(d0xb)/(d0xa**3)
       
       d0nut = d0nut + 2.D0*d0nus - d0nub + d0nud
       
    ENDDO
    
    IF(    d0nut.GT.0.D0)THEN 
       d0nut = 1.D0/d0nut
    ELSEIF(d0nut.EQ.0.D0)THEN
       d0nut = 0.D0
    ELSE
       WRITE(6,*)'ERROR IN FD0K11PS'
       STOP
    END IF
    
    fd0k11ps = EXP(-d0xa**2)*(d0xa**5)*d0nut
    
    RETURN
    
  END FUNCTION fd0k11ps
 
  FUNCTION fd0k12ps(d0xa2)
    
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:i0xa,d2x,d2y,d2z,d2bcf    
    
    INTEGER(i0ikind)::i0xb
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k12ps
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2) ! V/Vt_a

    d0nut = 0.D0

    DO i0xb  = 1, i0smax

       d0x   = d2x(  i0xb,i0xa)  ! Vt_a/Vt_b
       d0y   = d2y(  i0xb,i0xa)  ! m_b/m_a
       d0z   = d2z(  i0xa,i0xb)  ! T_a/T_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x          ! V/Vt_b
       d0nud =      d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nus = 2.D0*d0bcf*d0z*(1.D0+d0y)*func_G(d0xb)/d0xa 
       d0nub = 2.D0*d0bcf*func_G(d0xb)/(d0xa**3)
       
       d0nut = d0nut + 2.D0*d0nus - d0nub + d0nud
    ENDDO

    IF(    d0nut.GT.0.D0)THEN 
       d0nut = 1.D0/d0nut
    ELSEIF(d0nut.EQ.0.D0)THEN
       d0nut = 0.D0
    ELSE
       WRITE(6,*)'ERROR IN FD0K12PS'
       STOP
    END IF
    
    fd0k12ps = EXP(-d0xa**2)*(d0xa**7)*d0nut
    
    RETURN
    
  END FUNCTION fd0k12ps

  FUNCTION fd0k22ps(d0xa2)
    
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:i0xa,d2x,d2y,d2z,d2bcf    
    
    INTEGER(i0ikind)::i0xb
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k22ps
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2) ! V/Vt_a

    d0nut = 0.D0

    DO i0xb  = 1, i0smax

       d0x   = d2x(  i0xb,i0xa)  ! Vt_a/Vt_b
       d0y   = d2y(  i0xb,i0xa)  ! m_b/m_a
       d0z   = d2z(  i0xa,i0xb)  ! T_a/T_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x          ! V/Vt_b
       d0nud =      d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nus = 2.D0*d0bcf*d0z*(1.D0+d0y)*func_G(d0xb)/d0xa 
       d0nub = 2.D0*d0bcf*func_G(d0xb)/(d0xa**3)
       
       d0nut = d0nut + 2.D0*d0nus - d0nub + d0nud

    ENDDO
    
    IF(    d0nut.GT.0.D0)THEN 
       d0nut = 1.D0/d0nut
    ELSEIF(d0nut.EQ.0.D0)THEN
       d0nut = 0.D0
    ELSE
       WRITE(6,*)'ERROR IN FD0K22PS'
       STOP
    END IF

    fd0k22ps = EXP(-d0xa**2)*(d0xa**9)*d0nut
    
    RETURN
    
  END FUNCTION fd0k22ps

  FUNCTION fd0k11bn(d0xa2)
        
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:i0xa,d2x,d2bcf
    INTEGER(i0ikind)::i0xb
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k11bn
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nud,d0nut
    
    d0xa = SQRT(d0xa2) !C V/Vt_a

    d0nut = 0.D0

    DO i0xb  = 1, i0smax

       d0x   = d2x(  i0xb,i0xa)  ! Vt_a/Vt_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x          ! V/Vt_b
       d0nud = d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
      
       d0nut = d0nut + d0nud
    ENDDO
    
    fd0k11bn = EXP(-d0xa**2)*(d0xa**3)*d0nut
    
    RETURN
    
  END FUNCTION fd0k11bn

  FUNCTION fd0k12bn(d0xa2)
    
    USE LIBT2,ONLY:func_phi,func_G
    USE T2COMM,ONLY:i0xa,d2x,d2bcf

    INTEGER(i0ikind)::i0xb
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k12bn
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0nud,d0nut
    
    d0xa = SQRT(d0xa2) ! V/Vt_a

    d0nut = 0.D0

    DO i0xb  = 1, i0smax
       d0x   = d2x(  i0xb,i0xa)  ! Vt_a/Vt_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! V/Vt_b
       d0nud = d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nut = d0nut + d0nud
    ENDDO
    
    fd0k12bn = EXP(-d0xa**2)*(d0xa**5)*d0nut
    
    RETURN
    
  END FUNCTION fd0k12bn

  FUNCTION fd0k22bn(d0xa2)
    
    USE LIBT2,ONLY:func_phi,func_G
    USE T2COMM,ONLY:i0xa,d2x,d2bcf

    INTEGER(i0ikind)::i0xb
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k22bn
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0nud,d0nut
    
    d0xa = SQRT(d0xa2) ! V/Vt_a
    
    d0nut = 0.D0
    
    DO i0xb  = 1, i0smax
       d0x   = d2x(  i0xb,i0xa)  ! Vt_a/Vt_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! V/Vt_b
       d0nud = d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nut = d0nut + d0nud
    ENDDO
    
    fd0k22bn = EXP(-d0xa**2)*(d0xa**7)*d0nut
    
    RETURN
    
  END FUNCTION fd0k22bn
END MODULE T2CALV
