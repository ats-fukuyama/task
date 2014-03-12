!C--------------------------------------------------------------------
!C
!C MODULE T2CALV
!C 
!C
!C      CALCULATION OF COEFFICIENTS 
!C      OF ADVECTION-DIFFUSION EQUATIONS 
!C      FOR MULTI-FLUID TRANSPORT MODEL 
!C 
!C                       2014-02-23
!C
!C-------------------------------------------------------------------
MODULE T2CALV
  
  USE T2CNST, ONLY:i0ikind,i0rkind
  USE T2COMM, ONLY:i0smax
  
  IMPLICIT NONE
  
  INTEGER(i0ikind)::&
       i0midi
  REAL(   i0rkind)::&
       d0cogrr,d0cogrp,d0cogpp,d0cogtt,&
       d0ctgrr,d0ctgrp,d0ctgpp,d0ctgtt,&
       d0sqrtg,d0mfcr, d0rzcr, d0ugr,  d0ugp,  &
       d0ctbp, d0cobt, d0ctbt,d0bp2,d0bt2,d0bb,d0bb2
  
  PRIVATE
  
  PUBLIC T2_CALV
  
CONTAINS
  
  SUBROUTINE T2_CALV
    
    USE T2COMM, ONLY: i0mmax,i0cchk
    USE T2CCHK, ONLY: T2_CCHK

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
    
    IF(i0cchk.EQ.1) CALL T2_CCHK
    
    RETURN
    
  END SUBROUTINE T2_CALV
  
  !C---------------------------------------------------------
  !C
  !C CALCULATION OF FUNDAMENTAL PHYSICAL QUANTITIES 
  !C 
  !C     2014-03-12 H.SETO
  !C            
  !C---------------------------------------------------------
  
  SUBROUTINE T2CALV_PQ
    
    USE T2CNST
    
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,d0fpcst,d0ubcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,d0qpcst,d0wbcst,&
         i0xa,i0vmax,d0rmjr,d0iar,&
         i2crt,d2mtrc,&
         d2xvec_befor,d2mfc1,d2rzm,d2jm1,d2ws,&
         d1ee,d1mm,d1nn,d1ni,d1pp,d1pi,d1tt,d1ti,d1pa,d1pz,&
         d1ur,d1up,d1ut,d1ub,d1u2,&
         d1qr,d1qp,d1qt,d1qb,&
         d1wb,d1wt,d1wp,d1vt,d1hex,&
         d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4,&
         d2x,d2y,d2z,d2bcf
    
    USE LIBT2, ONLY:integ_f
    
    INTEGER(i0ikind)::&
         i0xid1d,i0xid2d,i0widi,i0vidi,i0sidi,i0sidj,i0sidk,i0nflag
    
    REAL(   i0rkind)::&
         d0mm_a,d0ee_a,d0ee_b,d0vti_a,&
         d0nn_a,d0nn_b,d0pp_a,d0tt_a,d0ti_a,&
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
         d0psip,&
         d0cps,d0cbn,d0wv2,d0wv3,d0temp,d0temp2,d0temp3,&
         d0k11,d0k11ps,d0k11bn,&
         d0k12,d0k12ps,d0k12bn,&
         d0k22,d0k22ps,d0k22bn
    
    REAL(i0rkind),DIMENSION(1:i0smax)::&
         d1fr,d1fb,d1ft,d1fp,d1vti,&
         d1nvcm1,d1nvcm2,d1nvcm3
    REAL(i0rkind),DIMENSION(1:i0smax,1:i0smax)::&
         d2clog,d2cfreq,&
         d2bmem00,d2bmem01,d2bmem10,d2bmem11,&
         d2bmen00,d2bmen01,d2bmen10,d2bmen11,&
         d2nfcl11,d2nfcl12,d2nfcl21,d2nfcl22
    REAL(i0rkind)::d0err

    
    !C
    !C d0psip : DERIVATIVE POLOIDAL FLUX FUNCTION
    !C          WITH RESPECT TO RHO               : \psi'
    !C d0pcf  : POLOIDAL CURRENT FUNCTION         : I
    !C
    !C d1nn   : PARTICLE DENSITY                  : n_{a}
    !C d1ur   : CT RADIAL   FLOW                  : u_{a}^{\rho} 
    !C d1ub   :    PARALLEL FLOW                  : u_{a\para}
    !C d1ut   : CO TOROIDAL FLOW                  : u_{a\zeta}
    !C d1up   : CT POLOIDAL FLOW                  : u_{a}^{\chi}
    !C
    !C d1pp   : PRESSURE                          : p_{a}   
    !C d1qr   : CT RADIAL   TOTAL HEAT FLUX       : Q_{a}^{\rho}
    !C d1qb   :    PARALLEL TOTAL HEAT FLUX       : Q_{a\para}
    !C d1qt   : CO TOROIDAL TOTAL HEAT FLUX       : Q_{a\zeta}
    !C d1qp   : CT POLOIDAL TOTAL HEAT FLUX       : Q_{a}^{\chi}
    !C 
    !C d1fr   : CT RADIAL   FLOW                  : n_{a}u_{a}^{\rho} 
    !C d1fb   :    PARALLEL PARTICLE FLUX         : n_{a}n_{u\para}
    !C d1ft   : CO TOROIDAL FLOW                  : n_{a}u_{a\zeta}
    !C d1fp   : CT POLOIDAL FLOW                  : n_{a}u_{a}^{\chi}
    !C
    !C d1wb  : PARALLEL TOTAL HEAT FLOW FUNCTION  : Q_{a\para}/p_{a}
    !C 
    !C * CO = COVARIANT, CT = CONTRAVARIANT
    !C
    !C      checked 2014-02-20 by H.Seto
    !C
   
    !C
    !C MASS AND ELECTRIC CHARGE [SI]
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


    !C
    !C CONVERT VARIABLES TO SI-UNIT
    !C
    
    DO i0sidi = 1, i0smax
       
       i0nflag = 0
       i0vidi =  10*i0sidi - 5
       
       !C d0nn_a: 10^{20}m^{-3}
       !C d0pp_a: keV*10^{20}m^{-3}
       d0nn_a = d2xvec_befor(i0vidi+1,i0xid2d)*d0nncst*1.D-20
       d0pp_a = d2xvec_befor(i0vidi+6,i0xid2d)*d0ppcst*1.D-23/d0aee
       
       d1nn(i0sidi) = d0nncst*d2xvec_befor(i0vidi+ 1,i0xid2d)
       d1fr(i0sidi) = d0frcst*d2xvec_befor(i0vidi+ 2,i0xid2d)*d0mfcr
       d1fb(i0sidi) = d0fbcst*d2xvec_befor(i0vidi+ 3,i0xid2d)
       d1ft(i0sidi) = d0ftcst*d2xvec_befor(i0vidi+ 4,i0xid2d)
       d1fp(i0sidi) = d0fpcst*d2xvec_befor(i0vidi+ 5,i0xid2d)
       
       d1pp(i0sidi) = d0ppcst*d2xvec_befor(i0vidi+ 6,i0xid2d)
       d1qr(i0sidi) = d0qrcst*d2xvec_befor(i0vidi+ 7,i0xid2d)*d0mfcr
       d1qb(i0sidi) = d0qbcst*d2xvec_befor(i0vidi+ 8,i0xid2d)
       d1qt(i0sidi) = d0qtcst*d2xvec_befor(i0vidi+ 9,i0xid2d)
       d1qp(i0sidi) = d0qpcst*d2xvec_befor(i0vidi+10,i0xid2d)
       
       !write(10,'(2(A4,I5),8(A4,D15.8))')'ND=',i0midi,'SP=',i0sidi,&
       !     'NN=',d1nn(i0sidi),'FR=',d1fr(i0sidi),&
       !     'FB=',d1fb(i0sidi),'FT=',d1ft(i0sidi),&
       !     'PP=',d1pp(i0sidi),'QR=',d1qr(i0sidi),&
       !     'QB=',d1qb(i0sidi),'QT=',d1qt(i0sidi)
       
       IF(        d0nn_a .GT.0.D0 )THEN
          d1ni(i0sidi) = 1.D0/d1nn(i0sidi)
       ELSE
          i0nflag = 1
       ENDIF
       
       IF(        d0pp_a .GT.0.D0 )THEN
          d1pi(i0sidi) = 1.D0/d1pp(i0sidi)
       ELSE
          i0nflag = 1
       ENDIF
       
       IF(i0nflag.EQ.1)THEN
          WRITE(6,*)'NEGATIVE  DENSITY or PRESSURE'
          WRITE(6,*)'SPECIS=',i0sidi,'NODE=',i0xid2d,&
               'N',d0nn_a,'10^{20} /m3','P=',d0pp_a,'keV*10^{20}/m^{3}'
          STOP       
       ENDIF
    ENDDO
    
    !C
    !C GEOMETRICAL COEFFICIENTS
    !C
    !C      checked 2014-02-20 by H.Seto
    !C

    !C 
    !C d0sqrtg : \sqrt{g}
    !C
    !C d0cogrr : g_{\rho\rho}*\rho
    !C d0cogrp : g_{\rho\chi}
    !C d0cogpp : g_{\chi\chi}
    !C d0cogtt : g_{\zeta\zeta}
    !C
    !C d0ctgrr : g^{\rho\rho}
    !C d0ctgrp : g^{\rho\chi}
    !C d0ctgpp : g^{\chi\chi}*rho
    !C d0ctgtt : g^{\zeta\zeta} 
    !C
    !C d0ugr   : CT RADIAL   FLAME MOVING VELOCITY      : u_{g}^{\rho}
    !C d0ugr   : CT POLOIDAL FLAME MOVING VELOCITY      : u_{g}^{\rho}
    !C
    !C d0ctbp  : CT POLOIDAL MAGNETIC FIELD             : B^{\chi}
    !C d0bp2   : SQUARED INTENSITY 
    !C                   of POLOIDAL MAGENTIC FIELD     : Bp^{2}
    !C d0bt2   : SQUARED INTENSITY 
    !C                   of TOROIDAL MAGENTIC FIELD     : Bt^{2}
    !C d0bb    : INTENSITY OF MAGNETIC FIELD            : B 
    !C d0bb2   : SQUARED INTENSITY of MAGNETIC FIELD    : B^{2} 
    !C
    
    d0psip  = d0mfcst*d2xvec_befor(1,i0xid1d)
    d0cobt  = d0btcst*d2xvec_befor(2,i0xid1d)

    d0sqrtg = d2jm1(1,i0midi)
    d0cogrr = d2jm1(2,i0midi)
    d0cogrp = d2jm1(3,i0midi)
    d0cogpp = d2jm1(4,i0midi)
    d0cogtt = d2jm1(5,i0midi)
    d0ctgrr = d2jm1(6,i0midi)
    d0ctgrp = d2jm1(7,i0midi)
    d0ctgpp = d2jm1(8,i0midi)
    d0ctgtt = d2jm1(9,i0midi)

    d0ctbp = d0psip/d0sqrtg
    d0ctbt = d0cobt*d0ctgtt
    
    d0ugr = 0.D0
    d0ugp = 0.D0
    
    d0bp2 = (d0ctbp**2)*d0cogpp
    d0bt2 = (d0cobt**2)*d0ctgtt

    d0bb2 = d0bp2 + d0bt2
    d0bb  = SQRT(d0bb2)
    
    DO i0sidi = 1, i0smax 
       
       d1ur(i0sidi) = d1fr(i0sidi)*d1ni(i0sidi)
       d1ub(i0sidi) = d1fb(i0sidi)*d1ni(i0sidi)
       d1ut(i0sidi) = d1ft(i0sidi)*d1ni(i0sidi)
       d1up(i0sidi) = d1fp(i0sidi)*d1ni(i0sidi)
    
       d1u2(i0sidi) = d1ub(i0sidi)*d1ub(i0sidi)
       
       d1wb(i0sidi) = d1qb(i0sidi)*d1pi(i0sidi)
       d1wt(i0sidi) = d1qt(i0sidi)*d1pi(i0sidi)
       d1wp(i0sidi) = d1qp(i0sidi)*d1pi(i0sidi)
       
       d1tt(i0sidi) = d1pp(i0sidi)*d1ni(i0sidi)
       d1ti(i0sidi) = d1nn(i0sidi)*d1pi(i0sidi)
       
    ENDDO
    
    !C
    !C SET WORKING SCALAR FOR DIFFERENTIAL (GROBAL)
    !C
    !C D2WS(1   ,:) : MAGNETIC FIELD INTENSITY    : B   
    !C D2WS(2   ,:) : MAJOR RADIUS                : R   
    !C D2WS(2N+1,:) : PARALLEL FLOW               : u_{a\para}  
    !C D2WS(2N+2,:) : PARALLEL HEAT FLOW FUNCTION : Q_{a\para}/p_{a}
    !C
    !C      checked 2014-02-20 by H.Seto
    !C
    
    d2ws(         1,i0midi) = d0bb
    d2ws(         2,i0midi) = d0rzcr
    
    DO i0sidi = 1, i0smax
       i0widi = 2*i0sidi
       d2ws(i0widi+1,i0midi) = d1ub(i0sidi)/d0ubcst
       d2ws(i0widi+2,i0midi) = d1wb(i0sidi)/d0wbcst
    ENDDO
    
    !C
    !C FOR COLLISION AND VISCOUS TERMS
    !C
    !C
    !C THERMAL VELOCITY [m/s]
    !C v_{ta} = \sqrt{2T_{a}/M_{a}}
    !C
    !C      checked 2014-02-20 by H.Seto
    !C
    
    DO i0sidi = 1, i0smax
       d0mm_a = d1mm(i0sidi)
       d0tt_a = d1tt(i0sidi)
       d0ti_a = d1ti(i0sidi)
       d1vt( i0sidi) = SQRT(2.0D0/d0mm_a*d0tt_a)
       d1vti(i0sidi) = SQRT(0.5D0*d0mm_a*d0ti_a)
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
    !C      checked 2014-02-12 by H.Seto
    !C
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       d0clog_ab = d2clog(i0sidi,i0sidj)
       d0nn_b    = d1nn( i0sidj)
       d0ee_b    = d1ee( i0sidj)
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
    !C      checked 2014-02-12 by H.Seto
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
       

    !C PARALLEL FRICTION COEFFICIENTS 
    !C          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
    !C
    !C DEFINITION OF VARIABLRS 
    !C            FOR NEOCLASSICAL FRICTION COEFFICIENTS
    !C
    !C d2nfcf1: \bar{l}^{ab}_{1}
    !C d2nfcf2: \bar{l}^{ab}_{2}
    !C d2nfcf3: \bar{l}^{ab}_{3}
    !C d2nfcf4: \bar{l}^{ab}_{4}
    
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       
       d0nfcl11_ab = d2nfcl11(i0sidi,i0sidj)
       d0nfcl12_ab = d2nfcl12(i0sidi,i0sidj)
       d0nfcl21_ab = d2nfcl21(i0sidi,i0sidj)
       d0nfcl22_ab = d2nfcl22(i0sidi,i0sidj)
       
       d0nfcf1_ab =   d0nfcl11_ab + d0nfcl12_ab
       d0nfcf2_ab = - d0nfcl12_ab * 0.4D0 
       d0nfcf3_ab = - d0nfcl21_ab + d0nfcl22_ab
       d0nfcf4_ab =   d0nfcl22_ab * 0.4D0
       
       d0nfcf3_ab = 2.5D0*d0nfcf1_ab+d0nfcf3_ab
       d0nfcf4_ab = 2.5D0*d0nfcf2_ab+d0nfcf4_ab
       
       d2nfcf1(i0sidi,i0sidj) = d0nfcf1_ab
       d2nfcf2(i0sidi,i0sidj) = d0nfcf2_ab
       d2nfcf3(i0sidi,i0sidj) = d0nfcf3_ab
       d2nfcf4(i0sidi,i0sidj) = d0nfcf4_ab
       
    ENDDO
    ENDDO
    
    !C
    !C HEAT EXCHANGE COEFFICIENT
    !C CHECKED 2013-02-20
    !C
    
    DO i0sidi = 1, i0smax
       d1hex(i0sidi) = 0.D0
       DO i0sidj = 1, i0smax
          d0zi_ab       = d2z(    i0sidj,i0sidi)
          d0cfreq_ab    = d2cfreq(i0sidi,i0sidj)
          d1hex(       i0sidi) &
               = d1hex(i0sidi) + 1.5D0*(1.D0 - d0zi_ab)*d0cfreq_ab
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
    !C d1nvcm1 : NEOCLASSICAL VISCOSITY COEFFICIENT     : \mu_{1a}
    !C d1nvcm2 : NEOCLASSICAL VISCOSITY COEFFICIENT     : \mu_{2a} 
    !C d1nvcm3 : NEOCLASSICAL VISCOSITY COEFFICIENT     : \mu_{3a} 
    !C

    d0temp  = SQRT(d0mfcr)
    d0temp  = d0iar*d0temp
    d0temp  = SQRT(d0temp)
    d0temp2 = (d0temp**3)*(d0ctbp**2)
    d0temp3 = 3.D0*(1.D0-1.46D0*d0temp)
   
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
       
       d0cbn = 2.92D0*d0mm_a*d0nn_a*d0bt2*d0de/d0temp3
       
       d0k11bn = d0cbn*d0k11bn
       d0k12bn = d0cbn*d0k12bn
       d0k22bn = d0cbn*d0k22bn
       
       !C
       !C MODIFIED VISCOSITY COEFFICIENT IN INTER-REGIMES 
       !C
       
       d0k11 = d0k11bn*d0k11ps/(d0k11bn + d0k11ps*d0temp2)
       d0k12 = d0k12bn*d0k12ps/(d0k12bn + d0k12ps*d0temp2)
       d0k22 = d0k22bn*d0k22ps/(d0k22bn + d0k22ps*d0temp2)
       
       !C
       !C NEOCLASSICAL VISCOSITY COEFFICIENTS: \mu_{ai}
       !C
       
       d1nvcm1(i0xa) = d0k11 
       d1nvcm2(i0xa) = d0k12 - 2.5D0*d0k11 
       d1nvcm3(i0xa) = d0k22 - 5.0D0*d0k12 + 6.25D0*d0k11
       
    ENDDO
    

    d0temp  = SQRT(d0mfcr)
    d0temp  = d0iar*d0temp
    d0temp  = SQRT(d0temp)
    d0temp2 = (d0temp**3)*(d0ctbp**2)
    d0temp3 = 3.D0*(1.D0-1.46D0*d0temp)
   
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
       
       d0cbn = 2.92D0*d0mm_a*d0nn_a*d0bt2*d0de/d0temp3
       
       d0k11bn = d0cbn*d0k11bn
       d0k12bn = d0cbn*d0k12bn
       d0k22bn = d0cbn*d0k22bn
       
       !C
       !C MODIFIED VISCOSITY COEFFICIENT IN INTER-REGIMES 
       !C
       
       d0k11 = d0k11bn*d0k11ps/(d0k11bn + d0k11ps*d0temp2)
       d0k12 = d0k12bn*d0k12ps/(d0k12bn + d0k12ps*d0temp2)
       d0k22 = d0k22bn*d0k22ps/(d0k22bn + d0k22ps*d0temp2)
       
       !C
       !C NEOCLASSICAL VISCOSITY COEFFICIENTS: \mu_{ai}
       !C
       
       d1nvcm1(i0xa) = d0k11 
       d1nvcm2(i0xa) = d0k12 - 2.5D0*d0k11 
       d1nvcm3(i0xa) = d0k22 - 5.0D0*d0k12 + 6.25D0*d0k11
       
    ENDDO
    
    !C d1nvcc1 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{1a}
    !C d1nvcc2 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{2a} 
    !C d1nvcc3 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{3a} 
    !C d1nvcc3 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{3a} 
    
    DO i0sidi = 1, i0smax
       
       d0nvcm1_a = d1nvcm1(i0sidi)
       d0nvcm2_a = d1nvcm2(i0sidi)
       d0nvcm3_a = d1nvcm3(i0sidi)
       
       d0nvcc1_a = d0nvcm1_a - d0nvcm2_a
       d0nvcc2_a = d0nvcm2_a * 0.4D0
       d0nvcc3_a = d0nvcm2_a - d0nvcm3_a
       d0nvcc4_a = d0nvcm3_a * 0.4D0
       
       d0nvcc3_a = 2.5D0*d0nvcc1_a + d0nvcc3_a
       d0nvcc4_a = 2.5D0*d0nvcc2_a + d0nvcc4_a
       
       d1nvcc1(i0sidi) = d0nvcc1_a
       d1nvcc2(i0sidi) = d0nvcc2_a
       d1nvcc3(i0sidi) = d0nvcc3_a
       d1nvcc4(i0sidi) = d0nvcc4_a

    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_PQ
  
  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF MASS SCALAR COEFFICIENTS
  !C
  !C             2014-03-12 H.SETO
  !C
  !C---------------------------------------------------------
  
  SUBROUTINE T2CALV_MS
    
    USE T2CNST, ONLY:&
         d0vci2
    
    USE T2COMM, ONLY:&
         & d0mfcst,d0btcst,d0etcst,d0epcst,        &
         & d0nncst,        d0fbcst,d0ftcst,        &
         & d0ppcst,        d0qbcst,d0qtcst,        &
         & d0mffct,d0btfct,d0etfct,d0epfct,        &
         & d0nnfct,        d0fbfct,d0ftfct,        &
         & d0ppfct,        d0qbfct,d0qtfct,        &
         & i0vmax,i0smax,d1mm,d3ms

    REAL(   i0rkind)::&
         d0mm_a
    INTEGER(i0ikind)::&
         i0vidi,i0vidj,i0sidi,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       d3ms(i0vidi,i0vidj,i0midi) = 0.D0
    ENDDO
    ENDDO

    !C
    !C EQ_001
    !C
    
    !C PSI'
    i0vidi = 1
    i0vidj = 1
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg                * d0mfcst/d0mffct
    
    !C
    !C EQ_002
    !C
    
    !C I
    
    i0vidi = 2
    i0vidj = 2
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg                * d0btcst/d0btfct
    

    !C
    !C EQ_003
    !C
    
    !C Et
    i0vidi = 3
    i0vidj = 3
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0vci2         * d0etcst/d0etfct
    
    
    !C
    !C EQ_004
    !C
    
    !C Ep
    i0vidi = 4
    i0vidj = 4
    d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0mfcr*d0vci2  * d0epcst/d0epfct
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       d0mm_a = d1mm(i0sidi)
       
       !C
       !C EQ_006
       !C
       
       i0vidi = i0vofi - 4
       
       !C N 
       i0vidj = i0vofi - 4
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg             * d0nncst/d0nnfct
       
       !C
       !C EQ_008
       !C 
       
       i0vidi = i0vofi - 2
       
       !C Fb
       i0vidj = i0vofi - 2
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0mm_a*d0bb * d0fbcst/d0fbfct
       
       !C
       !C EQ_009
       !C

       i0vidi = i0vofi - 1
       
       !C Ft
       i0vidj = i0vofi - 1
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0mm_a      * d0ftcst/d0ftfct
       
       !C
       !C EQ_011
       !C

       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*1.5D0       * d0ppcst/d0ppfct 
       
       !C
       !C EQ_013
       !C

       i0vidi = i0vofi + 3
       
       !C Qb       
       i0vidj = i0vofi + 3
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0mm_a*d0bb * d0qbcst/d0qbfct
       
       !C
       !C EQ_014
       !C

       i0vidi = i0vofi + 4
       
       !C Qt       
       i0vidj = i0vofi + 4
       d3ms(i0vidi,i0vidj,i0midi) = d0sqrtg*d0mm_a      * d0qtcst/d0qtfct 
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_MS
  
  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ADVECTION VECTOR COEFFICIENTS
  !C
  !C             2014-03-12 H.SETO
  !C
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_AV
    
    USE T2CNST, ONLY:&
         d0vci2
    
    USE T2COMM, ONLY:&
         & d0mfcst,d0btcst,d0etcst,d0epcst,d0ercst,&
         & d0nncst,        d0fbcst,d0ftcst,        &
         & d0ppcst,        d0qbcst,d0qtcst,        &
         & d0mffct,d0btfct,d0etfct,d0epfct,d0erfct,&
         & d0nnfct,        d0fbfct,d0ftfct,        &
         & d0ppfct,        d0qbfct,d0qtfct,        &
         & d1mm,d1tt,d1nn,d1pp,d1ni,d1pi,&
         & d1ub,d1ur,d1ut,d1up,d1u2,d1qr,d1qb,d1qt,d1qp,&
         & i0vmax,i0dmax,d4av
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0vidi,i0vidj,i0vofi
    
    REAL(   i0rkind)::&
         d0mm_a,d0tt_a,d0nn_a,d0pp_a,&
         d0ni_a,d0pi_a,&
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
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C PSI'
    i0vidj = 1
    d4av(1,i0vidi,i0vidj,i0midi) = d0ugr * d0mfcst/d0mffct
    d4av(2,i0vidi,i0vidj,i0midi) = d0ugp * d0mfcst/d0mffct
    
    !C
    !C EQ_002
    !C

    i0vidi = 2
    
    !C I
    i0vidj = 2
    d4av(1,i0vidi,i0vidj,i0midi) = d0ugr  * d0btcst/d0btfct
    d4av(2,i0vidi,i0vidj,i0midi) = d0ugp  * d0btcst/d0btfct
    
    !C
    !C EQ_003
    !C

    i0vidi = 3
    
    !C PSI'
    i0vidj = 1
    d4av(1,i0vidi,i0vidj,i0midi) = -d0sqrtg*d0ctgrr * d0mfcst/d0etfct
    d4av(2,i0vidi,i0vidj,i0midi) = -d0sqrtg*d0ctgrp * d0mfcst/d0etfct
    
    !C Et
    i0vidj = 3
    d4av(1,i0vidi,i0vidj,i0midi) =  d0ugr*d0vci2    * d0etcst/d0etfct
    d4av(2,i0vidi,i0vidj,i0midi) =  d0ugp*d0vci2    * d0etcst/d0etfct
    
    !C
    !C EQ_004
    !C
    
    i0vidi = 4
    
    !C Ep
    i0vidj = 4
    d4av(1,i0vidi,i0vidj,i0midi)&
         =  d0ugr*d0vci2*d0mfcr * d0epcst/d0epfct
    d4av(2,i0vidi,i0vidj,i0midi)&
         =  d0ugp*d0vci2*d0mfcr * d0epcst/d0epfct
    
    !C
    !C EQ_005
    !C
    
    i0vidi = 5
    
    !C Ep
    i0vidj = 4
    d4av(1,i0vidi,i0vidj,i0midi)&
         =  d0sqrtg*d0ctgrp*d0mfcr * d0epcst/d0erfct
    d4av(2,i0vidi,i0vidj,i0midi)&
         =  d0sqrtg*d0ctgpp        * d0epcst/d0erfct
    
    !C Er
    i0vidj = 5
    d4av(1,i0vidi,i0vidj,i0midi)&
         =  d0sqrtg*d0ctgrr * d0ercst/d0erfct
    d4av(2,i0vidi,i0vidj,i0midi)&
         =  d0sqrtg*d0ctgrp * d0ercst/d0erfct
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
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
       !C EQ_006
       !C
       
       i0vidi = i0vofi - 4
       
       !C N
       i0vidj = i0vofi - 4 
       d4av(1,i0vidi,i0vidj,i0midi)&
            = (d0ugr + d0sqrtg*d0ur_a) * d0nncst/d0nnfct
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (d0ugp + d0sqrtg*d0up_a) * d0nncst/d0nnfct
       
       !C
       !C EQ_008
       !C

       i0vidi = i0vofi - 2
       
       !C Fb
       i0vidj = i0vofi - 2
       d4av(1,i0vidi,i0vidj,i0midi)&
            = d0ugr*d0bb*d0mm_a  * d0fbcst/d0fbfct
       d4av(2,i0vidi,i0vidj,i0midi) &
            = (d0ugr*d0bb + d0sqrtg*d0ctbp*d0ub_a)*d0mm_a&
            * d0fbcst/d0fbfct
       
       !C P
       i0vidj = i0vofi + 1
       d4av(2,i0vidi,i0vidj,i0midi)&
            = d0sqrtg*d0ctbp * d0ppcst/d0fbfct
       
       !C
       !C EQ_009
       !C
       
       i0vidi = i0vofi - 1
       
       !C Ft
       i0vidj = i0vofi - 1 
       d4av(1,i0vidi,i0vidj,i0midi)&
            = (d0ugr + d0sqrtg*d0ur_a)*d0mm_a * d0ftcst/d0ftfct
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (d0ugp + d0sqrtg*d0up_a)*d0mm_a * d0ftcst/d0ftfct
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       d0x1 =  0.5D0*d0sqrtg*d0mm_a*d0nn_a*d0u2_a*d0pi_a
       d0x2 =  d0sqrtg*d0pi_a
       
       !C P
       i0vidj = i0vofi + 1
       d4av(1,i0vidi,i0vidj,i0midi)&
            = (1.5D0*d0ugr - d0x1*d0ur_a + d0x2*d0qr_a)&
            * d0ppcst/d0ppfct
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (1.5D0*d0ugp - d0x1*d0up_a + d0x2*d0qp_a)&
            * d0ppcst/d0ppfct
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C Fb
       i0vidj = i0vofi - 2
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (d0qb_a - 1.5D0*d0pp_a*d0ub_a)&
            * d0sqrtg*d0ctbp*d0ni_a*d0mm_a * d0fbcst/d0qbfct
       
       !C P
       i0vidj = i0vofi + 1
       d4av(2,i0vidi,i0vidj,i0midi)&
            = 2.5D0*d0sqrtg*d0tt_a*d0ctbp  * d0ppcst/d0qbfct
       
       !C Qb
       i0vidj = i0vofi + 3
       d4av(1,i0vidi,i0vidj,i0midi)&
            =  d0ugr*d0bb*d0mm_a           * d0qbcst/d0qbfct
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (d0ugp*d0bb + d0sqrtg*d0ctbp*d0ub_a)&
            * d0mm_a                       * d0qbcst/d0qbfct
       
       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4
       
       !C Ft
       i0vidj = i0vofi - 1
       d4av(1,i0vidi,i0vidj,i0midi)&
            = (d0qr_a - 1.5D0*d0pp_a*d0ur_a)&
            * d0sqrtg*d0ni_a*d0mm_a           * d0ftcst/d0qtfct
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (d0qp_a - 1.5D0*d0pp_a*d0up_a)&
            * d0sqrtg*d0ni_a*d0mm_a           * d0ftcst/d0qtfct
       
       !C Qt
       i0vidj = i0vofi + 4
       d4av(1,i0vidi,i0vidj,i0midi)&
            = (d0ugr + d0sqrtg*d0ur_a)*d0mm_a * d0qtcst/d0qtfct
       d4av(2,i0vidi,i0vidj,i0midi)&
            = (d0ugp + d0sqrtg*d0up_a)*d0mm_a * d0qtcst/d0qtfct
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_AV

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ADVECTION TENSOR COEFFICIENTS
  !C
  !C             2014-03-05 H.SETO
  !C
  !C---------------------------------------------------------  
  SUBROUTINE T2CALV_AT
    
    USE T2CNST
    
    USE T2COMM, ONLY:&
         & i0vmax,i0dmax,i0wmax,&
         &                         d0ftcst,d0fpcst,&
         &                         d0qtcst,d0qpcst,&
         &                 d0fbfct,d0ftfct,        &
         & d0ppfct,        d0qbfct,d0qtfct,        &
         & d1tt,d1ni,d1pi,d1ub,d1ur,d1up,&
         & d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,d6at    
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0widi,i0vidi,i0vidj,i0vofi
    
    REAL(   i0rkind)::&
         d0mm_a,d0ni_a,d0pi_a,d0tt_a,d0ur_a,d0up_a,d0ub_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    
    REAL(   i0rkind)::&
         d0y1,  d0y2,  d0y3,&
         d0x1,  d0x2,  d0x3,  d0x4,  d0x5,  d0x6,&
         d0c1_a,d0c2_a,d0c3_a,d0c4_a,&
         d0usr_a,d0usp_a
    
    !C
    !C INITIALIZATION
    !C

    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0widi = 1, i0wmax
          DO i0didj = 1, i0dmax
          DO i0didi = 1, i0dmax
             d6at(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
           
    d0y1 = (d0sqrtg*d0bt2)/(d0bb2*d0bb)
    d0y2 = d0ctbp/d0cobt
    d0y3 = d0ctbp/d0bb2
    
    d0x1 =  d0y1*d0y2     *d0ctbp       *2.D0
    d0x2 =  d0y1          *d0ctbp       *2.D0
    d0x3 =  d0y1     *d0y3*d0ctbp       *3.D0
    d0x4 =  d0y1     *d0y3       *d0cobt*3.D0
    d0x5 =  d0y1*d0y2
    d0x6 =  d0y1
    
    DO i0sidi = 1, i0smax
    
       i0vofi = 10*i0sidi

       d0ur_a = d1ur(i0sidi)
       d0up_a = d1up(i0sidi)
       d0ub_a = d1ub(i0sidi)
       d0ni_a = d1ni(i0sidi)
       d0pi_a = d1pi(i0sidi)
       d0tt_a = d1tt(i0sidi)


       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       d0nvcc3_a = d1nvcc3(i0sidi)
       d0nvcc4_a = d1nvcc4(i0sidi) 
       
       d0usr_a = d0ur_a
       d0usp_a = d0up_a -3.D0*d0ub_a*d0ctbp/d0bb 
       
       d0c1_a = d0nvcc1_a*d0ni_a
       d0c2_a = d0nvcc2_a*d0pi_a
       d0c3_a = d0nvcc3_a*d0ni_a*d0tt_a
       d0c4_a = d0nvcc4_a*d0pi_a*d0tt_a
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C Ft: (B) 
       i0vidj = i0vofi - 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c1_a         * d0ftcst/d0fbfct
       
       !C Fp: (B)
       i0vidj = i0vofi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c1_a         * d0fpcst/d0fbfct
       
       !C Qt: (B)
       i0vidj = i0vofi + 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c2_a         * d0qtcst/d0fbfct
       
       !C Qp: (B)
       i0vidj = i0vofi + 5
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c2_a         * d0qpcst/d0fbfct
       
       !C
       !C EQ_009
       !C
       
       i0vidi = i0vofi - 1
       
       !C Ft: (B)
       i0vidj = i0vofi - 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x3*d0c1_a         * d0ftcst/d0ftfct

       !C Fp (B)
       i0vidj = i0vofi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x4*d0c1_a         * d0fpcst/d0ftfct
       
       !C Qt (B)
       i0vidj = i0vofi + 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x3*d0c2_a         * d0qtcst/d0ftfct

       !C Qp (B)
       i0vidj = i0vofi + 5
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x4*d0c2_a         * d0qpcst/d0ftfct
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C Ft: (B)
       i0vidj = i0vofi - 1
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x5*d0c1_a*d0usr_a * d0ftcst/d0ppfct
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x5*d0c1_a*d0usp_a * d0ftcst/d0ppfct
       
       !C Fp: (B)
       i0vidj = i0vofi
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x6*d0c1_a*d0usr_a * d0fpcst/d0ppfct
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x6*d0c1_a*d0usp_a * d0fpcst/d0ppfct
       
       !C Qt: (B)
       i0vidj = i0vofi + 4
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x5*d0c2_a*d0usr_a * d0qtcst/d0ppfct
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x5*d0c2_a*d0usp_a * d0qtcst/d0ppfct

       !C Qp: (B)
       i0vidj = i0vofi + 5
       i0widi = 1
       d6at(1,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x6*d0c2_a*d0usr_a * d0qpcst/d0ppfct
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x6*d0c2_a*d0usp_a * d0qpcst/d0ppfct
       
       !C
       !C EQ_013
       !C

       i0vidi = i0vofi + 3
       
       !C Ft: (B)
       i0vidj = i0vofi - 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c3_a         * d0ftcst/d0qbfct
       
       !C Fp: (B)
       i0vidj = i0vofi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c3_a         * d0fpcst/d0qbfct
       
       !C Qt: (B)
       i0vidj = i0vofi + 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c4_a         * d0qtcst/d0qbfct
       
       !C Qp: (B)
       i0vidj = i0vofi + 5
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c4_a         * d0qpcst/d0qbfct
       
       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4
       
       !C Ft: (B)
       i0vidj = i0vofi - 1
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x3*d0c3_a         * d0ftcst/d0qtfct

       !C Fp: (B)
       i0vidj = i0vofi
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x4*d0c3_a         * d0fpcst/d0qtfct
       
       !C Qt: (B)
       i0vidj = i0vofi + 4
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x3*d0c4_a         * d0qtcst/d0qtfct
       
       !C Qp: (B)
       i0vidj = i0vofi + 5
       i0widi = 1
       d6at(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x4*d0c4_a         * d0qpcst/d0qtfct
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_AT
  
  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF DIFFUSION TENSOR COEFFICIENTS
  !C
  !C             2014-03-05 H.SETO
  !C
  !C---------------------------------------------------------  
  SUBROUTINE T2CALV_DT
    
    USE T2CNST
    USE T2COMM,ONLY:&
         & i0vmax,i0dmax,i0vmax,&
         & d0nncst,        d0fbcst,                &
         & d0ppcst,        d0qbcst,                &
         &                 d0fbfct,d0ftfct,        &
         & d0ppfct,        d0qbfct,d0qtfct,        &
         & d1ni,d1pi,d1tt,d1ur,d1up,d1ub,d1wb,&
         & d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         & d5dt
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0vidi,i0vidj,i0vofi
    
    REAL(   i0rkind)::&
         d0ni_a,d0pi_a,d0tt_a,&
         d0ur_a,d0up_a,d0ub_a,d0wb_a,&
         d0usr_a,  d0usp_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    REAL(   i0rkind)::&
         d0y1,&
         d0x1,  d0x2,  d0x3,&
         d0c1_a,d0c2_a,d0c3_a,d0c4_a

    !C
    !C INITIALIZATION
    !C
    
    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d5dt(i0didi,i0didj,i0vidi,i0vidj,i0midi) = 0.D0 
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    d0y1 = d0sqrtg*d0ctbp/d0bb
    d0x1 = d0y1*2.D0*d0ctbp
    d0x2 = d0y1*3.D0*d0ctbp*d0cobt/d0bb2
    d0x3 = d0y1
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       d0tt_a     = d1tt(i0sidi)
       d0ni_a     = d1ni(i0sidi)
       d0pi_a     = d1pi(i0sidi)
       d0ur_a     = d1ur(i0sidi)
       d0up_a     = d1up(i0sidi)
       d0ub_a     = d1ub(i0sidi)
       d0wb_a     = d1wb(i0sidi)

       d0nvcc1_a  = d1nvcc1(i0sidi)
       d0nvcc2_a  = d1nvcc2(i0sidi)
       d0nvcc3_a  = d1nvcc3(i0sidi)
       d0nvcc4_a  = d1nvcc4(i0sidi)
       
       d0usr_a = d0ur_a
       d0usp_a = d0up_a - d0ub_a*3.D0*d0ctbp/d0bb
       
       d0c1_a = d0nvcc1_a*d0ni_a
       d0c2_a = d0nvcc2_a*d0pi_a
       d0c3_a = d0nvcc3_a*d0ni_a*d0tt_a
       d0c4_a = d0nvcc4_a*d0pi_a*d0tt_a
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C N
       i0vidj = i0vofi - 4
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c1_a        *d0ub_a * d0nncst/d0fbfct
       
       !C Fb
       i0vidj = i0vofi - 2
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c1_a                * d0fbcst/d0fbfct
       
       !C P
       i0vidj = i0vofi + 1
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c2_a        *d0wb_a * d0ppcst/d0fbfct
       
       !C Qb
       i0vidj = i0vofi + 3
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c2_a                * d0qbcst/d0fbfct
       
       !C
       !C  EQ_009
       !C 
       
       i0vidi = i0vofi - 1
       
       !C N
       i0vidj = i0vofi - 4
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c1_a        *d0ub_a * d0nncst/d0ftfct
       !C Fb
       i0vidj = i0vofi - 2
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c1_a                * d0fbcst/d0ftfct
       !C P
       i0vidj = i0vofi + 1
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c2_a        *d0wb_a * d0ppcst/d0ftfct
       !C Qb
       i0vidj = i0vofi + 3
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c2_a                * d0qbcst/d0ftfct
       
       !C
       !C EQ_011
       !C

       i0vidi = i0vofi + 1
       
       !C N
       i0vidj = i0vofi - 4
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            = -d0x3*d0c1_a*d0usr_a*d0ub_a * d0nncst/d0ppfct
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0x3*d0c1_a*d0usp_a*d0ub_a * d0nncst/d0ppfct
       !C Fb
       i0vidj = i0vofi - 2
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            =  d0x3*d0c1_a*d0usr_a        * d0fbcst/d0ppfct
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0x3*d0c1_a*d0usp_a        * d0fbcst/d0ppfct
       !C P
       i0vidj = i0vofi + 1
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            = -d0x3*d0c2_a*d0usr_a*d0wb_a * d0ppcst/d0ppfct
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            = -d0x3*d0c2_a*d0usp_a*d0wb_a * d0ppcst/d0ppfct
       !C Qb
       i0vidj = i0vofi + 3
       d5dt(1,2,i0vidi,i0vidj,i0midi) &
            =  d0x3*d0c2_a*d0usr_a        * d0qbcst/d0ppfct
       d5dt(2,2,i0vidi,i0vidj,i0midi) &
            =  d0x3*d0c2_a*d0usp_a        * d0qbcst/d0ppfct
       
       !C   
       !C  EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C N 
       i0vidj = i0vofi - 4
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c3_a*d0ub_a         * d0nncst/d0qbfct
       !C Fb
       i0vidj = i0vofi - 2
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c3_a                * d0fbcst/d0qbfct
       !C P 
       i0vidj = i0vofi + 1
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c4_a*d0wb_a         * d0ppcst/d0qbfct
       !C Qb
       i0vidj = i0vofi + 3
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c4_a                * d0qbcst/d0qbfct
       
       !C
       !C  EQ_014
       !C 
       
       i0vidi = i0vofi + 4
       
       !C N
       i0vidj = i0vofi - 4
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c3_a*d0ub_a         * d0nncst/d0qtfct
       !C Fb
       i0vidj = i0vofi - 2
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c3_a                * d0fbcst/d0qtfct
       !C P
       i0vidj = i0vofi + 1
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c4_a*d0wb_a         * d0ppcst/d0qtfct
       !C Qb
       i0vidj = i0vofi + 3
       d5dt(2,2,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c4_a                * d0qbcst/d0qtfct
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_DT

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF GRADIENT VECTOR COEFFICIENTS
  !C
  !C             2014-03-05 H.SETO
  !C
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_GV

    USE T2CNST
    USE T2COMM, ONLY:&
         & i0vmax,i0dmax,&
         &         d0btcst,d0etcst,d0epcst,d0ercst,&
         & d0nncst,                                &
         & d0ppcst,                                &
         & d0mffct,d0btfct,        d0epfct,        &
         &         d0frfct,                        &
         & d0ppfct,d0qrfct,                        &
         & d1pp,d1tt,d1ur,d1up,d4gv
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0vidi,i0vidj,i0vofi
    REAL(   i0rkind)::&
         d0tt_a,d0ur_a,d0up_a,&
         d0c1_a,d0c2_a,d0x1,d0x2
    
    !C
    !C INITIALIZATION
    !C

    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0didi = 1, i0dmax
          d4gv(i0didi,i0vidi,i0vidj,i0midi) = 0.D0
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C Et
    i0vidj = 3
    d4gv(1,i0vidi,i0vidj,i0midi) = -d0sqrtg        * d0etcst/d0mffct

    !C
    !C EQ_002
    !C
    
    i0vidi = 2
    
    !C Ep
    i0vidj = 4 
    d4gv(1,i0vidi,i0vidj,i0midi) =  d0cogtt*d0mfcr * d0epcst/d0btfct
    
    !C Er
    i0vidj = 5
    d4gv(2,i0vidi,i0vidj,i0midi) = -d0cogtt        * d0ercst/d0btfct
    
    !C 
    !C EQ_004
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    d4gv(1,i0vidi,i0vidj,i0midi) =  d0cogpp        * d0btcst/d0epfct
    
    d0x1 = d0sqrtg*d0ctgrr*d0ctbp
    d0x2 = d0sqrtg*d0ctgrp*d0ctbp
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       d0tt_a = d1tt(i0sidi)
       d0ur_a = d1ur(i0sidi)
       d0up_a = d1up(i0sidi)

       d0c1_a = 2.5D0*(d0tt_a**2)
       d0c2_a = 5.0D0* d0tt_a
       
       !C
       !C EQ_007
       !C
       
       i0vidi = i0vofi - 3
       
       !C P
       i0vidj = i0vofi + 1
       d4gv(1,i0vidi,i0vidj,i0midi) =  d0x1           * d0ppcst/d0frfct
       d4gv(2,i0vidi,i0vidj,i0midi) =  d0x2           * d0ppcst/d0frfct
       
       !C 
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       d4gv(1,i0vidi,i0vidj,i0midi) = -d0sqrtg*d0ur_a * d0ppcst/d0ppfct
       d4gv(2,i0vidi,i0vidj,i0midi) = -d0sqrtg*d0up_a * d0ppcst/d0ppfct
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi + 2
       
       !C N
       i0vidj = i0vofi - 4
       d4gv(1,i0vidi,i0vidj,i0midi) = -d0x1*d0c1_a    * d0nncst/d0qrfct
       d4gv(2,i0vidi,i0vidj,i0midi) = -d0x2*d0c1_a    * d0nncst/d0qrfct
       
       !C P
       i0vidj = i0vofi + 1
       d4gv(1,i0vidi,i0vidj,i0midi) =  d0x1*d0c2_a    * d0ppcst/d0qrfct
       d4gv(2,i0vidi,i0vidj,i0midi) =  d0x2*d0c2_a    * d0ppcst/d0qrfct
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_GV

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF GRADIENT TENSOR COEFFICIENTS
  !C
  !C             2014-03-06 H.SETO
  !C
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_GT
    
    USE T2CNST
    USE T2COMM,ONLY:&
         & i0smax,i0vmax,i0wmax,i0dmax,            &
         & d0nncst,        d0fbcst,                &
         & d0ubcst,                                &
         & d0ppcst,        d0qbcst,                & 
         &                 d0fbfct,                &
         & d0ppfct,        d0qbfct,                &
         & d1ni,d1pi,d1tt,d1ub,d1ut,d1up,d1wb,     &
         & d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,d6gt 
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0widi,i0vidi,i0vidj,i0vofi
    
    REAL(   i0rkind)::&
         d0ni_a,d0pi_a,d0tt_a,&
         d0up_a,d0ub_a,d0ut_a,d0wb_a,d0udp_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    REAL(   i0rkind)::&
         d0x1,d0x2,&
         d0c1_a,d0c2_a,d0c3_a,d0c4_a,d0c5_a,d0c6_a
    
    !C
    !C INITIALIZATION
    !C

    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0widi = 1, i0wmax
          DO i0didj = 1, i0dmax
          DO i0didi = 1, i0dmax
             d6gt(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    d0x1 = 3.D0*d0sqrtg*(d0ctbp**2)/d0bb2
    d0x2 = 3.D0*d0sqrtg*d0bt2*d0ctbp/(d0bb2**2)
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       d0ni_a = d1ni(i0sidi)
       d0pi_a = d1pi(i0sidi)
       d0tt_a = d1tt(i0sidi)
       d0up_a = d1up(i0sidi)
       d0ut_a = d1ut(i0sidi)
       d0ub_a = d1ub(i0sidi)
       d0wb_a = d1wb(i0sidi)

       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       d0nvcc3_a = d1nvcc3(i0sidi)
       d0nvcc4_a = d1nvcc4(i0sidi)

       d0c1_a = d0nvcc1_a*d0ni_a
       d0c2_a = d0nvcc2_a*d0pi_a
       d0c3_a = d0nvcc3_a*d0ni_a*d0tt_a
       d0c4_a = d0nvcc4_a*d0pi_a*d0tt_a
       
       d0udp_a = d0ut_a*d0ctbp/d0cobt - d0up_a       

       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C N : (B)
       i0vidj = i0vofi - 4
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c1_a        *d0ub_a * d0nncst/d0fbfct
       
       !C Fb: (B)
       i0vidj = i0vofi - 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c1_a                * d0fbcst/d0fbfct
       
       !C P : (B)
       i0vidj = i0vofi + 1
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c2_a        *d0wb_a * d0ppcst/d0fbfct
       
       !C Qb: (B)
       i0vidj = i0vofi + 3
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c2_a                * d0qbcst/d0fbfct
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C N: (B)
       i0vidj = i0vofi - 4
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c1_a*d0udp_a*d0ub_a * d0nncst/d0ppfct
       
       !C Fb (B)
       i0vidj = i0vofi - 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c1_a*d0udp_a        * d0fbcst/d0ppfct

       !C P (B)
       i0vidj = i0vofi + 1 
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c2_a*d0udp_a*d0wb_a * d0ppcst/d0ppfct

       !C Qb (B)
       i0vidj = i0vofi + 3
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c2_a*d0udp_a        * d0qbcst/d0ppfct
       
       !C N: (Ub)
       i0vidj = i0vofi - 4
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c1_a        *d0ub_a * d0ubcst*d0nncst/d0ppfct
              
       !C Fb (Ub)
       i0vidj = i0vofi - 2
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c1_a                * d0ubcst*d0fbcst/d0ppfct
       
       !C P (Ub)
       i0vidj = i0vofi + 1
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c2_a        *d0wb_a * d0ubcst*d0ppcst/d0ppfct
       
       !C Qb (Ub)
       i0vidj = i0vofi + 3
       i0widi = 2*i0sidi + 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c2_a                * d0ubcst*d0qbcst/d0ppfct
       
       !C
       !C EQ_013
       !C
     
       i0vidi = i0vofi + 3
       
       !C N  (B)
       i0vidj = i0vofi - 4
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c3_a        *d0ub_a * d0nncst/d0qbfct
       
       !C Fb (B)
       i0vidj = i0vofi - 2
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c3_a                * d0fbcst/d0qbfct
       
       !C P  (B)
       i0vidj = i0vofi + 1
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c4_a        *d0wb_a * d0ppcst/d0qbfct

       !C Qb (B)
       i0vidj = i0vofi + 3
       i0widi = 1
       d6gt(2,2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c4_a                * d0qbcst/d0qbfct
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_GT

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF EXCITATION SCALAR COEFFICIENTS
  !C
  !C             2014-03-05 H.SETO
  !C
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_ES
    
    USE T2CNST,ONLY:&
         d0eps0,d0rmu0
    USE T2COMM,ONLY:&
         & i0vmax,&
         &                 d0ercst,d0epcst,d0etcst,&
         & d0nncst,d0frcst,d0fbcst,d0ftcst,d0fpcst,&
         & d0ppcst,d0qrcst,d0qbcst,d0qtcst,d0qpcst,&
         &         d0btfct,d0erfct,d0epfct,d0etfct,&
         &         d0frfct,d0fbfct,d0ftfct,d0fpfct,&
         & d0ppfct,d0qrfct,d0qbfct,d0qtfct,d0qpfct,&
         & d1ee,d1nn,d1pp,d1tt,d1ni,d1pi,d1hex,&
         & d2nfcf1,d2nfcf2,d2nfcf3,d2nfcf4,d3es
    
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj,i0vofi,i0vofj
    
    REAL(   i0rkind)::&
         d0ee_a,d0ee_b,d0nn_a,d0pp_a,d0tt_a,&
         d0ni_b,d0pi_b,&
         d0nfcf1_ab,d0nfcf2_ab,d0nfcf3_ab,d0nfcf4_ab,d0hex_a

    REAL(   i0rkind)::&
         d0y1,d0x1,d0x2,d0x3,d0x4,d0x5,d0x6,d0x7,d0x8,d0x9,d0x10
    
    REAL(   i0rkind)::&
         d0c1_a, d0c2_a, d0c3_a,&
         d0c1_ab,d0c2_ab,d0c3_ab,d0c4_ab

    !C
    !C INITIALIZATION
    !C

    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       d3es(i0vidi,i0vidj,i0midi) = 0.D0
    ENDDO
    ENDDO
    
    d0x1 = d0sqrtg*d0rmu0
    d0x2 = d0sqrtg/d0eps0

    !C
    !C EQ_002
    !C

    i0vidi = 2
    
    !C Ep
    i0vidj = 4
    d3es(i0vidi,i0vidj,i0midi) =  d0cogtt * d0epcst/d0btfct
    
    DO i0sidj = 1, i0smax
       
       i0vofj = 10*i0sidj
       
       d0ee_b  = d1ee(i0sidj)
       
       !C
       !C EQ_003
       !C
       
       i0vidi = 3
       
       !C Ft
       i0vidj = i0vofj - 1
       d3es(i0vidi,i0vidj,i0midi)&
            =  d0x1*d0ee_b                * d0ftcst/d0etfct
       
       !C
       !C EQ_004
       !C
       
       i0vidi = 4
       
       !C Fr
       i0vidj = i0vofj - 3
       d3es(i0vidi,i0vidj,i0midi)&
            =  d0x1*d0ee_b*d0cogrp*d0mfcr * d0frcst/d0epfct
       
       !C Fp
       i0vidj =i0vofj
       d3es(i0vidi,i0vidj,i0midi)&
            =  d0x1*d0ee_b*d0cogpp        * d0fpcst/d0epfct
       
       !C
       !C EQ_005
       !C
       
       i0vidi = 5
       
       !C N
       i0vidj = i0vofj - 4
       d3es(i0vidi,i0vidj,i0midi)&
            = -d0x2*d0ee_b                * d0nncst/d0erfct
       
    ENDDO
    d0y1 = d0sqrtg*d0ctbp
    d0x3  = d0y1*d0ctgrp*d0mfcr
    d0x4  = d0y1*d0ctgrr
    d0x5  = d0cobt*d0bb
    d0x6  = d0sqrtg*d0ctbt
    d0x7  = d0y1*d0mfcr
    d0x8  = d0sqrtg*d0bb
    d0x9  = d0sqrtg*d0x7
    d0x10 = d0y1*d0cogpp
 
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi

       d0nn_a  = d1nn( i0sidi)
       d0pp_a  = d1pp( i0sidi)
       d0tt_a  = d1tt( i0sidi)
       d0ee_a  = d1ee( i0sidi)
       d0hex_a = d1hex(i0sidi)


       d0c1_a = d0ee_a*d0nn_a
       d0c2_a = d0ee_a*d0pp_a*2.5D0

       !C
       !C EQ_007
       !C
       
       i0vidi = i0vofi - 3
       
       !C Ep
       i0vidj = 4
       d3es(i0vidi,i0vidj,i0midi) = -d0c1_a*d0x3  * d0epcst/d0frfct
       
       !C Er
       i0vidj = 5
       d3es(i0vidi,i0vidj,i0midi) = -d0c1_a*d0x4  * d0ercst/d0frfct
       
       !C Fb
       i0vidj = i0vofi - 2
       d3es(i0vidi,i0vidj,i0midi) = -d0ee_a*d0x5  * d0fbcst/d0frfct
       
       !C Ft
       i0vidj = i0vofi - 1
       d3es(i0vidi,i0vidj,i0midi) =  d0ee_a*d0bb2 * d0ftcst/d0frfct
       
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          d0ni_b = d1ni(i0sidj)
          d0pi_b = d1pi(i0sidj)
          
          d0nfcf1_ab = d2nfcf1(i0sidi,i0sidj)
          d0nfcf2_ab = d2nfcf2(i0sidi,i0sidj)
          
          d0c1_ab = d0nfcf1_ab*d0ni_b
          d0c2_ab = d0nfcf2_ab*d0pi_b

          !C
          !C EQ_008
          !C
          
          i0vidi = i0vofi - 2
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0x6*d0c1_a * d0etcst/d0fbfct
             
             !C Ep
             i0vidj = 4
             d3es(i0vidi,i0vidj,i0midi) = -d0x7*d0c1_a * d0epcst/d0fbfct
             
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          d3es(i0vidi,i0vidj,i0midi) = -d0x8*d0c1_ab   * d0fbcst/d0fbfct
          
          !C Qb
          i0vidj = i0vofj + 3
          d3es(i0vidi,i0vidj,i0midi) = -d0x8*d0c2_ab   * d0qbcst/d0fbfct
          
          !C
          !C EQ_009
          !C
          
          i0vidi = i0vofi - 1
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0sqrtg*d0c1_a * d0etcst/d0ftfct
             
             !C Fr
             i0vidj = i0vofj - 3
             d3es(i0vidi,i0vidj,i0midi) = -d0x9*d0ee_a    * d0frcst/d0ftfct
             
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          d3es(i0vidi,i0vidj,i0midi) = -d0sqrtg*d0c1_ab   * d0ftcst/d0ftfct
          
          !C Qt
          i0vidj = i0vofj + 4
          d3es(i0vidi,i0vidj,i0midi) = -d0sqrtg*d0c2_ab   * d0qtcst/d0ftfct
          
       ENDDO
       
       !C 
       !C EQ_010
       !C
       
       i0vidi = i0vofi 
       
       !C Fb
       i0vidj = i0vofi - 2
       d3es(i0vidi,i0vidj,i0midi) = -d0x8                 * d0fbcst/d0fpfct

       !C Ft
       i0vidj = i0vofi - 1
       d3es(i0vidi,i0vidj,i0midi) =  d0x6                 * d0ftcst/d0fpfct

       !C Fp
       i0vidj = i0vofi
       d3es(i0vidi,i0vidj,i0midi) =  d0x10                * d0fpcst/d0fpfct 
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1 
       d3es(i0vidi,i0vidj,i0midi) = d0sqrtg*d0hex_a       * d0ppcst/d0ppfct
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi  + 2
       
       !C Ep
       i0vidj = 4
       d3es(i0vidi,i0vidj,i0midi) = -d0c2_a*d0x3          * d0epcst/d0qrfct
       
       !C Er
       i0vidj = 5
       d3es(i0vidi,i0vidj,i0midi) = -d0c2_a*d0x4          * d0ercst/d0qrfct
       
       !C Qb
       i0vidj = i0vofi + 3
       d3es(i0vidi,i0vidj,i0midi) = -d0ee_a*d0x5          * d0qbcst/d0qrfct
       
       !C Qt
       i0vidj = i0vofi + 4
       d3es(i0vidi,i0vidj,i0midi) =  d0ee_a*d0bb2         * d0qtcst/d0qrfct
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          d0ni_b = d1ni(i0sidj)
          d0pi_b = d1pi(i0sidj)
          
          d0nfcf3_ab = d2nfcf3(i0sidi,i0sidj)
          d0nfcf4_ab = d2nfcf4(i0sidi,i0sidj)
          
          d0c3_ab = d0nfcf3_ab*d0tt_a*d0ni_b
          d0c4_ab = d0nfcf4_ab*d0tt_a*d0pi_b
          
          !C
          !C EQ_013
          !C
          
          i0vidi = i0vofi + 3 
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0c2_a*d0x6    * d0etcst/d0qbfct
             
             !C Ep
             i0vidj = 4
             d3es(i0vidi,i0vidj,i0midi) = -d0c2_a*d0x7    * d0epcst/d0qbfct
             
          ENDIF
          
          !C Fb  
          i0vidj = i0vofj - 2
          d3es(i0vidi,i0vidj,i0midi) = -d0c3_ab*d0x8      * d0fbcst/d0qbfct
          
          !C Qb
          i0vidj = i0vofj + 3
          d3es(i0vidi,i0vidj,i0midi) = -d0c4_ab*d0x8      * d0qbcst/d0qbfct
          
       
          !C
          !C EQ_014
          !C
          
          i0vidi = i0vofi + 4
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             d3es(i0vidi,i0vidj,i0midi) = -d0c2_a*d0sqrtg * d0etcst/d0qtfct
             
             !C Qr
             i0vidj = i0vofj + 2
             d3es(i0vidi,i0vidj,i0midi) = -d0ee_a*d0x9    * d0qrcst/d0qtfct
             
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          d3es(i0vidi,i0vidj,i0midi) = -d0c3_ab*d0sqrtg   * d0ftcst/d0qtfct
          
          !C Qt
          i0vidj = i0vofj + 4
          d3es(i0vidi,i0vidj,i0midi) = -d0c4_ab*d0sqrtg   * d0qtcst/d0qtfct
          
       ENDDO
       
       !C 
       !C EQUATION FOR Qp
       !C
       
       i0vidi = i0vofi + 5
       
       !C Qb
       i0vidj = i0vofi + 3
       d3es(i0vidi,i0vidj,i0midi) = -d0x8                 * d0qbcst/d0qpfct

       !C Qt
       i0vidj = i0vofi + 4
       d3es(i0vidi,i0vidj,i0midi) =  d0x6                 * d0qtcst/d0qpfct

       !C Qp
       i0vidj = i0vofi + 5
       d3es(i0vidi,i0vidj,i0midi) =  d0x10                * d0qpcst/d0qpfct 
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_ES

  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ECITATION VECTOR COEFFICIENTS
  !C
  !C             2014-03-05 H.SETO
  !C
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_EV
    
    USE T2CNST
    
    USE T2COMM,ONLY:&
         & i0smax,i0vmax,i0dmax,i0wmax,&
         & d0mfcst,        d0etcst,d0epcst,        &
         &                 d0fbcst,                &
         & d0ubcst,                                &
         &                 d0qbcst,                &
         & d0wbcst,                                &
         &                 d0etfct,                &
         &                 d0fbfct,                &
         &                 d0qbfct,d0qtfct,        &
         & d1ee,d1mm,d1nn,d1ub,d1ut,d1up,d1pp,d1qb,d1qt,d1ni,&
         & d1wb,d1wt,d1wp,&
         & d1nvcc1,d1nvcc2,d5ev
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0vidj,i0widi,i0didi,i0vofi
    REAL(   i0rkind)::&
         d0ee_a,d0mm_a,d0ni_a,d0pp_a,&
         d0ub_a,d0ut_a,d0up_a,&
         d0qb_a,&
         d0wb_a,d0wt_a,d0wp_a,&
         d0udp_a,d0wdp_a,&
         d0nvcc1_a,d0nvcc2_a
    REAL(   i0rkind)::&
         d0x1,d0x2,d0x3,d0x4,d0x5,d0x6,d0x7,d0x8,d0x9,&
         d0y1,d0y2,d0y3,d0y4,d0y5
    REAL(   i0rkind)::&
         d0c1_a,d0c2_a,d0c3_a,d0c4_a,d0c5_a
    
    !C 
    !C INITIALIZATION
    !C

    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0widi = 1, i0wmax 
          DO i0didi = 1, i0dmax   
             d5ev(i0didi,i0widi,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C EQ_003
    !C

    i0vidi = 3
 
    !C PSI' (R)
    i0vidj = 1
    i0widi = 2
    d5ev(1,i0widi,i0vidi,i0vidj,i0midi)&
         = d0sqrtg*2.D0*SQRT(d0ctgtt)*d0ctgrr * d0mfcst/d0etcst
    d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
         = d0sqrtg*2.D0*SQRT(d0ctgtt)*d0ctgrp * d0mfcst/d0etcst

    d0y1 = d0sqrtg*d0ctbt*3.D0              /d0bb
    d0y2 = d0sqrtg*d0ctbp*3.D0*d0mfcr       /d0bb
    d0y3 = d0sqrtg*(3.D0*d0bt2/d0bb2-1.D0)  /d0bb
    d0y4 = d0sqrtg*d0cobt*d0ctbp*d0mfcr*3.D0/(d0bb2*d0bb)
    d0y5 = d0ctbp/d0cobt
    
    d0x1 = d0sqrtg*d0ctbp/d0bb
    d0x2 = d0y1*d0bt2/d0bb2
    d0x3 = d0y1*d0ctbp
    d0x4 = d0y2*d0bt2/d0bb2
    d0x5 = d0y2*d0ctbp
    d0x6 = d0y3*d0bt2/d0bb2
    d0x7 = d0y3*d0ctbp
    d0x8 = d0y4*d0bt2/d0bb2
    d0x9 = d0y4*d0ctbp
  
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi

       d0ee_a    = d1ee(i0sidi)
       d0mm_a    = d1mm(i0sidi)
       d0pp_a    = d1pp(i0sidi)
       d0ni_a    = d1ni(i0sidi)
       d0ub_a    = d1ub(i0sidi)
       d0ut_a    = d1ut(i0sidi)
       d0up_a    = d1up(i0sidi)
       d0qb_a    = d1qb(i0sidi)
       d0wb_a    = d1wb(i0sidi)
       d0wt_a    = d1wt(i0sidi)
       d0wp_a    = d1wp(i0sidi)
       
       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       
       d0udp_a = d0y5*d0ut_a - d0up_a
       d0wdp_a = d0y5*d0wt_a - d0wp_a
       
       d0c1_a = d0ee_a*(d0nvcc1_a*d0udp_a+d0nvcc2_a*d0wdp_a) 
       d0c2_a = d0ee_a*d0nvcc1_a
       d0c3_a = d0ee_a*d0nvcc2_a
       d0c4_a = (d0qb_a-1.5*d0pp_a*d0ub_a)*d0mm_a*d0ni_a
       d0c5_a = d0ub_a*d0mm_a
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C Fb: (B)
       
       i0vidj = i0vofi - 2
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0ub_a*d0mm_a * d0fbcst/d0fbfct
       
       !C
       !C EQ_013
       !C

       i0vidi = i0vofi + 3
       
       !C Et (B)
       i0vidj = 3
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x2*d0c1_a * d0etcst/d0qbfct
       
       !C Et (Ub)
       i0vidj = 3
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x3*d0c2_a * d0ubcst*d0etcst/d0qbfct
              
       !C Et (Qb/P)
       i0vidj = 3
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x3*d0c3_a * d0wbcst*d0etcst/d0qbfct
       
       !C Ep (B)
       i0vidj = 4
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x4*d0c1_a * d0epcst/d0qbfct
       
       !C Ep (Ub)
       i0vidj = 4
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x5*d0c2_a * d0ubcst*d0epcst/d0qbfct
       
       !C Ep (Qb/P)
       i0vidj = 4
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x5*d0c3_a * d0wbcst*d0epcst/d0qbfct
       
       !C Fb (B)
       i0vidj = i0vofi - 2
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c4_a * d0fbcst/d0qbfct
       
       !C Qb (B) 
       i0vidj = i0vofi + 3
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            = -d0x1*d0c5_a * d0qbcst/d0qbfct
       
       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4
       
       !C Et (B)
       i0vidj = 3
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x6*d0c1_a * d0etcst/d0qtfct
       
       !C Et (Ub)
       i0vidj = 3
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x7*d0c2_a * d0ubcst*d0etcst/d0qtfct
       
       !C Et (Qb/P)
       i0vidj = 3
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x7*d0c3_a * d0wbcst*d0etcst/d0qtfct
       
       !C Ep (B)
       i0vidj = 4
       i0widi = 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x8*d0c1_a * d0epcst/d0qtfct
       
       !C Ep (Ub)
       i0vidj = 4
       i0widi = 2*i0sidi + 1
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x9*d0c2_a * d0ubcst*d0epcst/d0qtfct
       
       !C Ep (Qb)
       i0vidj = 4
       i0widi = 2*i0sidi + 2
       d5ev(2,i0widi,i0vidi,i0vidj,i0midi)&
            =  d0x9*d0c3_a * d0wbcst*d0epcst/d0qtfct
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_EV
  
  !C---------------------------------------------------------
  !C 
  !C CALCULATION OF ECITATION VECTOR COEFFICIENTS
  !C
  !C             2014-03-06 H.SETO
  !C
  !C---------------------------------------------------------
  SUBROUTINE T2CALV_ET
    
    USE T2CNST
    USE T2COMM,ONLY:&
         & i0smax,i0dmax,i0vmax,i0wmax,            &
         &                         d0ftcst,d0fpcst,&
         & d0ubcst,                                &
         &                         d0qtcst,d0qpcst,&
         & d0wbcst,                                &
         &                 d0fbfct,                &
         & d0ppfct,        d0qbfct,                &
         & d1ni,d1pi,d1tt,d1ut,d1ub,d1up,&
         & d1nvcc1,d1nvcc2,d1nvcc3,d1nvcc4,&
         & d7et
    
    INTEGER(i0ikind)::&
         i0sidi,i0didi,i0didj,i0widi,i0widj,i0vidi,i0vidj,i0vofi
    
    REAL(   i0rkind)::&
         d0tt_a,d0pi_a,d0ni_a,d0ut_a,d0ub_a,d0up_a,d0udp_a,&
         d0nvcc1_a,d0nvcc2_a,d0nvcc3_a,d0nvcc4_a
    REAL(   i0rkind)::&
         d0x1,  d0x2,  d0x3,  d0x4,  d0x5,  d0x6,&
         d0y1,  d0y2,  d0y3,  d0y4,&
         d0c1_a,d0c2_a,d0c3_a,d0c4_a,d0c5_a,d0c6_a
  
    !C
    !C
    !C
    DO i0vidj = 1, i0vmax
    DO i0vidi = 1, i0vmax
       DO i0widj = 1, i0wmax
       DO i0widi = 1, i0wmax
          DO i0didj = 1, i0dmax
          DO i0didi = 1, i0dmax
             d7et(i0didi,i0didj,i0widi,i0didj,i0vidi,i0vidj,i0midi) = 0.D0
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    d0y1 = 3.D0*d0sqrtg*d0bt2*d0ctbp/(d0bb2**2)
    d0y2 = 3.D0*d0sqrtg*(d0bt2**2)  /(d0bb2**3)
    d0y3 = 3.D0*d0sqrtg*d0bt2*d0ctbp/(d0bb2**2)
    d0y4 = d0ctbp/d0cobt

    d0x1 = d0y1*d0y4
    d0x2 = d0y1
    d0x3 = d0y2*d0y4
    d0x4 = d0y2
    d0x5 = d0y3*d0y4
    d0x6 = d0y3

    DO i0sidi = 1, i0smax

       i0vofi = 10*i0sidi

       d0ni_a    = d1ni(   i0sidi)
       d0pi_a    = d1pi(   i0sidi)
       d0tt_a    = d1tt(   i0sidi)
       d0ub_a    = d1ub(   i0sidi)
       d0ut_a    = d1ut(   i0sidi)
       d0up_a    = d1up(   i0sidi)
       d0nvcc1_a = d1nvcc1(i0sidi)
       d0nvcc2_a = d1nvcc2(i0sidi)
       d0nvcc3_a = d1nvcc3(i0sidi)
       d0nvcc4_a = d1nvcc4(i0sidi)
       

       d0udp_a = d0y4*d0ut_a - d0up_a 

       d0c1_a = d0nvcc1_a*d0ni_a
       d0c2_a = d0nvcc2_a*d0pi_a
       d0c3_a = d0nvcc3_a*d0ni_a*d0tt_a
       d0c4_a = d0nvcc4_a*d0pi_a*d0tt_a

       
       
       !C
       !C EQ_008
       !C

       i0vidi = i0vofi - 2
       
       !C Ft: (B,B)
       i0vidj = i0vofi - 1
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c1_a * d0ftcst/d0fbfct
       
       !C Fp: (B,B)
       i0vidj = i0vofi
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c1_a * d0fpcst/d0fbfct
       
       !C Qt: (B,B)
       i0vidj = i0vofi + 4
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c2_a * d0qtcst/d0fbfct

       !C Qp: (B,B)
       i0vidj = i0vofi + 5
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0x2*d0c2_a * d0qtcst/d0fbfct

       !C
       !C EQ_011
       !C

       i0vidi = i0vofi + 1

       !C Ft: (B,B)
       i0vidj = i0vofi - 1
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0x3*d0c1_a*d0udp_a * d0ftcst/d0ppfct 


       !C Fp: (B,B)
       i0vidj = i0vofi
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x4*d0c1_a*d0udp_a * d0fpcst/d0ppfct 

       !C Qt: (B,B)
       i0vidj = i0vofi + 4
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0x3*d0c2_a*d0udp_a * d0qtcst/d0ppfct 

       !C Qp (B,B)
       i0vidj = i0vofi + 5
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x4*d0c2_a*d0udp_a * d0qpcst/d0ppfct 

       !C Ft: (Ub,B)
       i0vidj = i0vofi - 1
       i0widi = 2*i0sidi + 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0x5*d0c1_a         * d0ubcst*d0ftcst/d0ppfct 

       !C Fp: (Ub,B)
       i0vidj = i0vofi
       i0widi = 2*i0sidi + 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x6*d0c1_a        * d0ubcst*d0fpcst/d0ppfct 


       !C Qt: (ub,B)
       i0vidj = i0vofi + 4
       i0widi = 2*i0sidi + 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            = -d0x5*d0c2_a        * d0ubcst*d0qtcst/d0ppfct 

       !C Qp (ub,B)
       i0vidj = i0vofi + 5
       i0widi = 2*i0sidi + 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x6*d0c2_a        * d0ubcst*d0qpcst/d0ppfct 

       !C
       !C EQUATION FOR Qb
       !C

       i0vidi = i0vofi + 3

       !C Ft (B,B)
       i0vidj = i0vofi - 1
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c3_a        * d0ftcst/d0qbfct

       !C Fp (B,B)
       i0vidj = i0vofi
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  -d0x2*d0c3_a       * d0fpcst/d0qbfct


       !C Qt (B,B)
       i0vidj = i0vofi + 4
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  d0x1*d0c4_a        * d0qtcst/d0qbfct

       !C Qp (B,B)
       i0vidj = i0vofi + 5
       i0widi = 1
       i0widj = 1
       d7et(2,2,i0widi,i0widj,i0vidi,i0vidj,i0midi)&
            =  -d0x2*d0c4_a       * d0qpcst/d0qbfct
       
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
