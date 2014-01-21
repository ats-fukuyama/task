!C 
!C
!C MODULE T2CALV
!C 
!C
!C      CALCULATION OF COEFFICIENTS 
!C      OF ADVECTION-DIFFUSION EQUATIONS 
!C      FOR MULTI-FLUID TRANSPORT MODEL
!C 
!C 
MODULE T2CALV
  
  USE T2CNST, ONLY:i0ikind,i0rkind
  USE T2COMM, ONLY:i0spcs
  
  IMPLICIT NONE
  
  INTEGER(i0ikind)::&
       i0nid
  !C
  !C FUNDAMENTAL VARIABLES
  !C

  REAL(   i0rkind)::&
       d0mf,&            !C DERIVATIVE OF POROIDAL FLUX IN RADIAL
       d0bp,&            !C CONTRAVARIANT POROIDAL MAGNETIC FIELD        
       d0bt,&            !C COVARIANT     TOROIDAL MAGNETIC FIELD 
       d0jm1,&           !C JACOBIAN OF MSC                       
       d0jm2,&           !C CONTRAVARIANT METRIC RHO-RHO   OF MSC 
       d0jm3,&           !C CONTRAVARIANT METRIC RHO-CHI   OF MSC 
       d0jm4,&           !C CONTRAVARIANT METRIC CHI-CHI   OF MSC 
       d0jm5             !C CONTRAVARIANT METRIC ZETA-ZETA OF MSC 
  REAL(   i0rkind)::&
       !C 
       !C GEOMETRICAL AND MAGNETIC VARIABLES
       !C
       d0ug1,d0ug2,&
       d0ugr,d0ugp,&
       d0bb,&
       d0bp2,&
       d0bt2,&
       d0bpb,&
       d0btb,&
       d0bpi,&
       d0q0,&
       d0g1,d0g2,d0g3,d0g4,d0g5
  
  PRIVATE
  
  PUBLIC T2_CALV
  
CONTAINS
  
  SUBROUTINE T2_CALV
    
    USE T2COMM, ONLY: i0nmax1
    REAL(4)::e0time1,e0time2
    
    DO i0nid = 1,i0nmax1
       
       CALL T2CALV_WV       
       CALL T2CALV_MS
       CALL T2CALV_AV
       CALL T2CALV_AT
       CALL T2CALV_DT
       CALL T2CALV_GV
       CALL T2CALV_GT
       CALL T2CALV_ES
       CALL T2CALV_EV
       CALL T2CALV_ET
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2_CALV
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF BASIC PHYSICAL QUANTITIES 
  !C 
  !C                                                2013-11-03 CHECKED
  !C------------------------------------------------------------------
  
  SUBROUTINE T2CALV_WV
    
    USE T2CNST
    USE T2COMM
    USE LIBT2, ONLY:integ_f
    
    INTEGER(i0ikind)::&
         i1,j1,k1,&
         i0nid1d,i0nid2d,i0vid1d,i0vid2d,i0wid,i0sid
    
    REAL(   i0rkind)::&
         d0w1,d0w2,d0w3,d0w4,d0w5,d0w6,d0w7,d0w8,d0w9,d0w3i,&
         d0n, d0p, d0fb,d0ub,d0qb,&
         d0ur,d0up,d0ut,&
         d0qr,d0qp,d0qt,&
         d0ni,d0pib,d0mfcr,d0rzcr,&
         d0k11ps,d0k11bn,d0k11,&
         d0k12ps,d0k12bn,d0k12,&
         d0k22ps,d0k22bn,d0k22,&
         d0err,d0nu,d0x,d0y,d0z,d0zi,d0nt,d0pt,&
         d0m00,d0m01,d0m11,&
         d0l11, d0l12, d0l21, d0l22,&
         d0lx11,d0lx12,d0lx21,d0lx22
    !C
    !C d0pf1d: DERIVATIVE POLOIDAL FLUX FUNCTION
    !C        WITH RESPECT TO RHO            : \psi'
    !C d0pc1d: POLOIDAL CURRENT FUNCTION     : I
    !C d0n2d : PARTICLE DENSITY              : n_{a}
    !C d0p2d : PRESSURE                      : p_{a}
    !C
    !C d0ur2d: CNT RADIAL   FLOW             : u_{a}^{\rho} 
    !C d0up2d: CNT POLOIDAL FLOW             : u_{a}^{\chi}
    !C d0ut2d: CO  TOROIDAL FLOW             : u_{a\zeta}
    !C d0ub2d:     PARALLEL FLOW             : u_{a\para}
    !C
    !C d0qr2d: CNT RADIAL   TOTAL HEAT FLUX  : Q_{a}^{\rho}
    !C d0qp2d: CNT POLOIDAL TOTAL HEAT FLUX  : Q_{a}^{\chi}
    !C d0qt2d: CO  TOROIDAL TOTAL HEAT FLUX  : Q_{a\zeta}
    !C d0qb2d:     PARALLEL TOTAL HEAT FLUX  : Q_{a\para}
    !C 
    !C d0fr2d: CNT RADIAL   FLOW             : n_{a}u_{a}^{\rho} 
    !C d0fp2d: CNT POLOIDAL FLOW             : n_{a}u_{a}^{\chi}
    !C d0ft2d: CO  TOROIDAL FLOW             : n_{a}u_{a\zeta}
    !C d0fb2d:     PARALLEL PARTICLE FLUX    : n_{a}n_{u\para}
    !C
    !C d0vb2d: PARALLEL TOTAL HEAT FLOW      : 2Q_{a\para}/5p_{a}
    !C 
    !C
    !C
    !C 
    !C
    !C 

    !C
    !C
    !C
    !C
    !C
    REAL(   i0rkind),DIMENSION(1:i0spcs)::&
         d1mu1,d1mu2,d1mu3,d1fr,d1fb,d1ft,d1ti,d1vti
    
    REAL(   i0rkind),DIMENSION(1:i0spcs,1:i0spcs)::&
         d2m00,d2m01,d2m10,d2m11,&
         d2n00,d2n01,d2n10,d2n11,&
         d2cl,d2nu
    
    !C MASS AND ELECTRIC CHARGE
    
    DO i0sida = 1,i0spcs
       d1m(i0sida) = d1pa(i0sida)*d0amp
       d1e(i0sida) = d1pz(i0sida)*d0aee
    ENDDO
    
    !C INITIALIZE     
    !i0nid2 = i2crt( i0nid,2) - 1
    !i0nid3 = i1mfc1(i0nid  ) - 1
    !d0rzcr = d2rzc1(i0nid,1)
    
    !C INITIALIZE 
    
    i0nid1d = i2crt( i0nid,3) - 1
    i0nid2d = i2crt( i0nid,2) - 1
    d0rzcr  = d2rzc1(i0nid,1)
    d0mfcr  = d2mfc1(i0nid,1)    
    
    !C CONVERT VARIABLES TO SI-UNIT
    !i0vid2 = i0vmax*i0nid2
    !i0vid3 = i0vmax*i0nid3
    
    !C CONVERT VARIABLES TO SI-UNIT
    i0vid1d = i0vmax*i0nid1d
    i0vid2d = i0vmax*i0nid2d
    
    d0psi = d0mfcst*d1guv_befor(i0vid1d+1)
    d0pcf = d0btcst*d1guv_befor(i0vid1d+2)
    
    DO i0sida = 1, i0spcs
       
       i0vid2d = i0vmax*i0nid2d + 8*i0sida - 3
       
       d0nn = d1guv_befor(i0vid2d+1)
       d0pp = d1guv_befor(i0vid2d+5)
       
       d1nn(i0sid) = d0nncst*d1guv_befor(i0vid2d+1)
       d1fr(i0sid) = d0frcst*d1guv_befor(i0vid2d+2)
       d1fb(i0sid) = d0fbcst*d1guv_befor(i0vid2d+3)
       d1ft(i0sid) = d0ftcst*d1guv_befor(i0vid2d+4)
       d1pp(i0sid) = d0ppcst*d1guv_befor(i0vid2d+5)
       d1qr(i0sid) = d0qrcst*d1guv_befor(i0vid2d+6)
       d1qb(i0sid) = d0qbcst*d1guv_befor(i0vid2d+7)
       d1qt(i0sid) = d0qtcst*d1guv_befor(i0vid2d+8)
       
       IF((d0nt.GE.0.D0).AND.(d0pt.GE.0.D0))THEN
          d1ni(i0sid) = 1.D0/d1nn(i0sid)
          d1pi(i0sid) = 1.D0/d1pp(i0sid)
       ELSEIF((ABS(d0nt).LE.1.D-5).AND.(ABS(d0pt).LE.1.D-5))THEN
          d1ni(i0sid) = 0.D0
          d1pi(i0sid) = 0.D0
       ELSE
          WRITE(6,*)'NEGATIVE TEMP or DENS'
          WRITE(6,*)'i1=',i1,'N',d0nn,'/m3','P=',d0pp,'keV/m3'
          STOP
       END IF
    ENDDO
    
    !C
    !C GEOMETRICAL COEFFICIENTS
    !C

    !C 
    !C d0sqrtg: \sqrt{g}
    !C
    !C d0cvgrr: g_{\rho\rho}
    !C d0cvgrc: g_{\rho\chi}
    !C d0cvgcc: g_{\chi\chi}
    !C d0cvgtt: g_{\zeta\zeta}
    !C
    !C d0ctgrr: g^{\rho\rho}
    !C d0ctgrc: g^{\rho\chi}
    !C d0ctgcc: g^{\chi\chi}
    !C d0ctgtt: g^{\zeta\zeta} 
    !C
    !C d0ugr  : CONTRA RADIAL FLAME MOVING VELOCITY: u_{g}^{\rho}
    !C d0ugr  : CONTRA PO FLAME MOVING VELOCITY: u_{g}^{\rho}
    !C
    !C d0bp  : CONTRA POLOIDAL MAGNETIC FIELD        : B^{\chi}
    !C d0bpi : INVERSE CONTRA POLOIDAL MAGNETIC FIELD: 1/B^{\chi}  
    !C d0bp2 : SQUARED INTENSITY 
    !C        of POLOIDAL MAGENTIC FIELD            : Bp*Bp
    !C d0bt2 : SQUARED INTENSITY 
    !C        of TOROIDAL MAGENTIC FIELD            : Bt*Bt
    !C d0bp2i: INVERSE SQUARED INTENSITY 
    !C        of POLOIDAL MAGNETIC FIELD            : 1/Bp^{2}
    !C d0bb  : INTENSITY OF MAGNETIC FIELD           : B 
    
    d0sqrtg = d2mtrcs(1,i0nid)
    d0cvgrr = d2mtrcs(2,i0nid)
    d0cvgrp = d2mtrcs(3,i0nid)
    d0cvgpp = d2mtrcs(4,i0nid)
    d0cvgtt = d2mtrcs(5,i0nid)
    
    d0ctgtt =  1.D0/d0cvgtt
    
    IF(d0sqrtg.GT.0.D0)THEN
       
       !C d0wv01  : WORKING VARIABLE: R^{2}/\sqrt{g}^{2}   
       d0wv01    =  d0cvgtt/(d0sqrtg*d0sqrtg)
       
       d0ctgrr =  d0cvgpp*d0wv01
       d0ctgrp = -d0cvgrp*d0wv01
       d0ctgpp =  d0cvgrr*d0wv01
       
       d0bp    = d0psi/d0sqrtg
       d0bpi   = 1.D0/d0bp
       
    ELSE
       
       d0ctgrr = 0.D0
       d0ctgrp = 0.D0
       d0ctgpp = 0.D0
       
       d0bp    = 0.D0
       d0bpi   = 0.D0
       
    ENDIF
    
    d0ugr = 0.D0
    d0ugp = 0.D0
    
    d0bp2 = (d0bp**2)*d0cvgcc
    d0bt2 = (d0bt**2)*d0ctgtt
    d0bb  = d0bp2 + d0bt2
    d0bb  = SQRT(d0bb)
    
    IF(d0bp2.GT.0.D0)THEN
       d0bp2i = 1.D0/d0bp2
    ELSE
       d0bp2i = 0.D0
    ENDIF
    
    d0g1 =  d0bp*d0bp2i*d0bb
    d0g2 = -d0bp*d0bp2i*d0bt*d0ctgtt
    
    DO i0sid = 1, i0spcs 
       
       d0ur = d1fr(i0sid)*d1ni(i0sid)
       d0ub = d1fb(i0sid)*d1ni(i0sid)
       d0ut = d1ft(i0sid)*d1ni(i0sid)
       
       d1ur(i0sid) = d0ur
       d1ub(i0sid) = d0ub
       d1ut(i0sid) = d0ut
       d1up(i0sid) = d0g1*d0ub + d0g2*d0ut
       d1u2(i0sid) = d0ub*d0ub !C FOR DEBUG 
       
       d1qp(i0sid) = d0g1*d1qb(i0sid) + d0g2*d1qt(i0sid)
       d1vb(i0sid) = 0.4D0*d1qb(i0sid)*d1pi(i0sid)
       
       d1tt(i0sid) = d1pp(i0sid)*d1ni(i0sid)
       d1ti(i0sid) = d1nn(i0sid)*d1pi(i0sid)
       
    ENDDO
    
    !C SET WORKING ARRAY FOR DIFFERENTIAL (GROBAL)
    !C
    !C D2WA(1   ,:) : MAGNETIC FIELD INTENSITY         : B   
    !C D2WA(2   ,:) : MAJOR RADIUS OF TOKAMKAK DEVICE  : R   
    !C D2WA(4N-1,:) : PARTICLE DENSITY                 : n_{a}   
    !C D2WA(4N-2,:) : PARALLEL FLOW                    : u_{a\para}  
    !C D2WA(4N+1,:) : PRESSURE                         : p_{a}   
    !C D2WA(4N+2,:) : PARALLEL TOTAL HEAT FLOW VELOCITY: 2Q_{a\para}/5p_{a}  
    !C
    
    d2wa(         1,i0nid) = d0bb
    d2wa(         2,i0nid) = d0rzcr
    
    DO i0sida = 1, i0spcs
       
       i0wid = 4*i0sida - 2
       
       d2wa(i0wid+1,i0nid) = d1nn(i0sid)
       d2wa(i0wid+2,i0nid) = d1ub(i0sid)
       d2wa(i0wid+3,i0nid) = d1pp(i0sid)
       d2wa(i0wid+4,i0nid) = d1vb(i0sid)
       
    ENDDO
    
    !C
    !C FOR COLLISION AND VISCOUS TERMS
    !C
    
    !C THERMAL VELOCITY [m/s]
    !C VT*VT = 2T/M
    DO i0sida = 1, i0spcs
       
       d0m_a  = d1m( i0sida)
       d0tt_a = d1tt(i0sida)
       d0ti_a = d1ti(i0sida)
       
       d1vt( i0sida)= SQRT(2.0D0/d0m_a*d0tt_a)
       d1vti(i0sida)= SQRT(0.5D0*d0m_a*d0ti_a)
    ENDDO
    
    
    !C COULOMB LOGARITHM
    !C ref: NRL PLASMA FORMULARY 2011
    !C FOR DEBUG lnA = 17
    !C d2clog: COULOMB LOGARITHM
    
    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs
       d2clog(i1,j1) = 17.D0
    ENDDO
    ENDDO
    
    DO i0sidb = 1, i0spcs
    DO i0sida = 1, i0spcs
       
       d0clog_ab = d2clog(i0sida,i0sidb)

       d0nn_b    = d1nn(i0sidb)
       d0e_b     = d1e(i0sidb)

       d0m_a     = d1m(i0sida)
       d0e_a     = d1e(i0sida)       
       d0vti_a   = d1vti(i0sida)
       
       d2bcf(i0sida,i0sidb)&
            = (d0nn_b*(d0e_a**2)*(d0e_b**2)*d0clog_ab*(d0vti_a**3))&
            / (4.D0*d0pi*(d0eps0**2)*(d0m_a**2))
    ENDDO
    ENDDO

    !DO i0sidb = 1, i0spcs
    !DO i0sida = 1, i0spcs
    !   
    !   d2bcf(i0sida,i0sidb)&
    !        = (d1nn(i0sidb)*(d1e(i0sida)**2)*(d1e(i0sidb)**2)&
    !        * d2cl(i0sida,i0sidb)*(d1vti(i0sida)**3))&
    !        / (4.D0*d0pi*(d0eps0**2)*(d1m(i0sida)**2))
    !ENDDO
    !ENDDO
    
    !C COLLISION FREQUENCY [1/s]
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002)  P.277
    !C CHECKED 2013-11-13
    
    DO i0sidb = 1, i0spcs
    DO i0sida = 1, i0spcs
       d2cfreq(i0sida,i0sidb) &
            = d2bcf(i0sida,i0sidb)/(0.75D0*SQRT(d0pi))
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
    !C d2x  = x_ab = v_{tb}/v_{ta}
    !C d2y  = y_ab = ma /mb
    !C d2z  = t_ab = ta /tb
    !C d0wv2   = 1/SQRT(1 + x_{ab}^{2})
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
    DO i0sidb = 1, i0spcs
    DO i0sida = 1, i0spcs 
       d2x(i0sida,i0sidb) = d1vti(i0sida)*d1vt(i0sidb)
       d2y(i0sida,i0sidb) = d1m(  i0sida)/d1m( i0sidb)
       d2z(i0sida,i0sidb) = d1tt( i0sida)*d1ti(i0sidb)
    ENDDO
    ENDDO


    DO i0sidb = 1, i0spcs
    DO i0sida = 1, i0spcs
       
       d0x_ab  = d2x(i0sida,i0sidb) 
       d0y_ab  = d2y(i0sida,i0sidb)
       d0z_ab  = d2z(i0sida,i0sidb)
       d0zi_ab = d2z(i0sidb,i0sida)
       
       d0wv2   = 1.D0/SQRT(1.D0 + d0x_ab**2)
       
       d0bmem00_ab = -        (1.D0+d0y_ab)*(d0wv2**3)
       d0bmem01_ab = -  1.5D0*(1.D0+d0y_ab)*(d0wv2**5)
       d0bmem11_ab = - (3.25D0+4.D0*(d0x_ab**2)&
                     +  7.50D0*(d0x_ab**4))*(d0wv2**5)
       
       d2bmem00(i0sida,i0sidb) = d0bmem00_ab
       d2bmem01(i0sida,i0sidb) = d0bmem01_ab
       d2bmem10(i0sida,i0sidb) = d0bmem01_ab
       d2bmem11(i0sida,i0sidb) = d0bmem11_ab
       
       d2bmen00(i0sida,i0sidb) = - d0bmem00_ab
       d2bmen01(i0sidb,i0sida) = - d0bmem01_ab*d0x_ab*d0zi_ab
       d2bmen10(i0sida,i0sidb) = - d0bmem01_ab
       d2bmen11(i0sida,i0sidb) =   6.75D0*d0z_ab*(d0x_ab**2)*(d0wv2**5)
       
    ENDDO
    ENDDO
    
    !C
    !C PARALLEL FRICTION COEFFICIENTS
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002) P.239 
    !C 
    !C DEFINITION OF VARIABLRS FOR PARALLEL FRICTION COEFFICIENTS
    !C
    !C d2pfcl11: l^{ab}_{11}
    !C d2pfcl12: l^{ab}_{12}
    !C d2pfcl21: l^{ab}_{21}
    !C d2pfcl22: l^{ab}_{22}

    
    d0w1 = 0.D0; d0w2 = 0.D0
    
    d0lx11 = 0.D0; d0lx12 = 0.D0; d0lx21 = 0.D0; d0lx22 = 0.D0

    d0pfcl11_ab = 0.D0
    d0pfcl12_ab = 0.D0
    d0pfcl21_ab = 0.D0
    d0pfcl22_ab = 0.D0
    
    DO i0sidb = 1, i0spcs
    DO i0sida = 1, i0spcs
          
       d0cfreq_ab = d2cfreq(i0sida,i0sidb)
       d0bmen00_ab = d2bmen00(i0sida,i0sidb)
       d0bmen01_ab = d2bmen01(i0sida,i0sidb)
       d0bmen10_ab = d2bmen10(i0sida,i0sidb)
       d0bmen11_ab = d2bmen11(i0sida,i0sidb)
       
       d0pfcl11_ab = d0bmen00_ab*d0cfreq_ab
       d0pfcl12_ab = d0bmen01_ab*d0cfreq_ab
       d0pfcl21_ab = d0bmen10_ab*d0cfreq_ab
       d0pfcl22_ab = d0bmen11_ab*d0cfreq_ab
       
       IF(i0sida.EQ.i0sidb)THEN
          
          DO i0sidc = 1, i0spcs
             
             d0cfreq_ac  = d2cfreq( i0sida,i0sidc)
             d0bmem00_ac = d2bmem00(i0sida,i0sidc)
             d0bmem01_ac = d2bmem01(i0sida,i0sidc)
             d0bmem10_ac = d2bmem10(i0sida,i0sidc)
             d0bmem11_ac = d2bmem11(i0sida,i0sidc)

             d0pfcl11_ab = d0pfcl11_ab + d0bmem00_ac*d0cfreq_ac
             d0pfcl12_ab = d0pfcl12_ab + d0bmem01_ac*d0cfreq_ac
             d0pfcl21_ab = d0pfcl21_ab + d0bmem10_ac*d0cfreq_ac
             d0pfcl22_ab = d0pfcl22_ab + d0bmem11_ac*d0cfreq_ac
             
          ENDDO
          
       ENDIF
       
       d0wv3 = d1m(isida)*d1nn(i0sida)
       
       d0pfcl11_ab = d0pfcl11_ab*d0wv3 ! l_{11}^{ab}
       d0pfcl12_ab = d0pfcl12_ab*d0wv3 ! l_{12}^{ab}
       d0pfcl21_ab = d0pfcl21_ab*d0wv3 ! l_{21}^{ab}
       d0pfcl22_ab = d0pfcl22_ab*d0wv3 ! l_{22}^{ab}
       
       d2pfcl11(i0sida,i0sidb) = d0pfcl11_ab
       d2pfcl12(i0sida,i0sidb) = d0pfcl12_ab
       d2pfcl21(i0sida,i0sidb) = d0pfcl21_ab
       d2pfcl22(i0sida,i0sidb) = d0pfcl22_ab
              
    ENDDO
    ENDDO


    !C PARALLEL FRICTION COEFFICIENTS 
    !C          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
    !C
    !C DEFINITION OF VARIABLRS FOR PARALLEL FRICTION COEFFICIENTS
    !C
    !C d2pfcf1: f^{ab}_{1}
    !C d2pfcf2: f^{ab}_{2}
    !C d2pfcf3: f^{ab}_{3}
    !C d2pfcf4: f^{ab}_{4}
    
    DO i0sidb = 1, i0spcs
    DO i0sida = 1, i0spcs
       
       d0pfcl11_ab = d2pfcl11(i0sida,i0sidb)
       d0pfcl12_ab = d2pfcl12(i0sida,i0sidb)
       d0pfcl21_ab = d2pfcl21(i0sida,i0sidb)
       d0pfcl22_ab = d2pfcl22(i0sida,i0sidb)
       
       d0ni_b = d1ni(i0sidb)
       d0pi_b = d1pi(i0sidb)
       
       d0wv4  = d1t(i0sida)/d1m(i0sdia)
       
       d0pfcf1_ab =   (d0pfcl11_ab+d0pfcl12_ab)*d0ni_b 
       d0pfcf2_ab = - 0.4D0*d0pfcl12_ab*d0pi_b 
       d0pfcf3_ab = - (d0pfcl21_ab+d0pfcl22_ab)*d0ni_b
       d0pfcf4_ab =   0.4D0*d0pfcl22_ab*d0pi_b
       
       d0pfcf3_ab = d0wv4*(2.5D0*pfcf1_ab+d0pfcf3_ab)
       d0pfcf4_ab = d0wv4*(2.5D0*pfcf2_ab+d0pfcf4_ab)
       
       d2pfcf1(i0sida,i0sidb) = d0pfcf1_ab
       d2pfcf2(i0sida,i0sidb) = d0pfcf2_ab
       d2pfcf3(i0sida,i0sidb) = d0pfcf3_ab
       d2pfcf4(i0sida,i0sidb) = d0pfcf4_ab
       
    ENDDO
    ENDDO
    
    !C
    !C HEAT EXCHANGE COEFFICIENT
    !C CHECKED 2013-11-14
    !C
    
    DO i0sida = 1, i0spcs
       d1hex(i0sida) = 0.D0
    ENDDO
    
    DO i0sida = 1, i0spcs
    DO i0sidb = 1, i0spcs
       d0zi_ab    = d2z(i0sidb,i0sidb)
       d0cfreq_ab = d2cfreq(i0sida,i0sidb)
       d1hex(i0sida) = d1hex(i0sida)+1.5D0*(1.D0 - d0zi_ab)*d0cfreq_ab
    ENDDO
    ENDDO
    
    !C
    !C NEOCLASSICAL PARALLEL COEFFICIENTS
    !C REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !C      P. HELANDER AND D.J. SIGMAR (2002)  CHAP.12 AND P.279
    !C
    !C ft/fc
    !C
    !C
    !C
    IF(d0mfcr.GT.0.D0)THEN
       IF(   (d0mfcr.GT.0.D0).AND.(d0mfcr.LE.1.D0))THEN
          !d0q0 = (d0qc-d0qs)*(1.D0 - d0mfcr**2)+d0qs
          d0q0 = (d0qc-d0qs)*((1.D0 - d0mfcr**2)**2)+d0qs
       ELSEIF(d0mfcr.GT.1.D0)THEN
          !d0q0 = (d0qs-d0qc)*(       d0mfcr**2)+d0qc
          d0q0 = d0qs
       ELSE
          WRITE(6,*)'WRONG RHO INPUT'
          print*,d0mfcr
          STOP
       ENDIF
       
       d0iar = d0mfcr/d0rmjr     !C INVERSE ASPECT RATIO (r/R0)
       d0w3 = 1.46D0*SQRT(d0iar)/(1.D0-1.46D0*SQRT(d0iar)) 
       d0w4 = (d0de*d0w3*2.D0*(d0q0**2)*(d0rzcr**2))/(3.D0*(d0w1**2))

       !d0w1 = d0mfcr/d0rmjr     !C INVERSE ASPECT RATIO (r/R0)
       !d0w2 = 1.46D0*SQRT(d0w1) !C ft
       !d0w3 = d0w2/(1.D0-d0w2)  !C ft/fc, fc = 1 - ft
       !d0w4 = (d0de*d0w3*2.D0*(d0q0**2)*(d0rzcr**2))/(3.D0*(d0w1**2))
       
       DO i0xa=1,i0spcs
          !C
          !C NC VISCOSITY COEFFICIENTS IN PFIRSCH-SCHLUTER REGIME 
          !C 
          CALL INTEG_F(fd0k11ps,d0k11ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k12ps,d0k12ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k22ps,d0k22ps,d0err,EPS=1.D-8,ILST=0)
          
          d0w5    = 0.4D0*d1p(i0xa)*d0de
          d0k11ps = d0w5*d0k11ps
          d0k12ps = d0w5*d0k12ps
          d0k22ps = d0w5*d0k22ps

          !C
          !C NC VISCOSITY COEFFICIENTS IN BANANA REGIME 
          !C 
          CALL INTEG_F(fd0k11bn,d0k11bn,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k12bn,d0k12bn,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k22bn,d0k22bn,d0err,EPS=1.D-8,ILST=0)
          
          d0w6 = d0w4*d1m(i0xa)*d1n(i0xa)
          
          d0k11bn = d0w6*d0k11bn
          d0k12bn = d0w6*d0k12bn
          d0k22bn = d0w6*d0k22bn
          
          !C
          !C INTER-REGIME NC VISCOSITY COEFFICIENTS 
          !C
          d0k11 = d0k11ps*d0k11bn/(d0k11ps+d0k11bn)
          d0k12 = d0k12ps*d0k12bn/(d0k12ps+d0k12bn)
          d0k22 = d0k22ps*d0k22bn/(d0k22ps+d0k22bn)
          
          d1mu1(i0xa) = d0k11 
          d1mu2(i0xa) = d0k12 - 2.5D0*d0k11 
          d1mu3(i0xa) = d0k22 - 5.0D0*d0k12 + 6.25D0*d0k11
       ENDDO

       DO i1=1,i0spcs
          
          d0w1 = 3.0D0*(d1mu1(i1)-d1mu2(i1))*d1ni(i1)
          d0w2 = 1.2D0*d1mu2(i1)*d1pi(i1)
          d0w3 = 3.0D0*(d1mu2(i1)-d1mu3(i1))*d1ni(i1)
          d0w4 = 1.2D0*d1mu3(i1)*d1pi(i1)
          d0w5 = d1t(i1)/d1m(i1)
       
          d1c1(i1) = d0w1
          d1c2(i1) = d0w2
          d1c3(i1) = d0w5*(2.5D0*d0w1+d0w3)
          d1c4(i1) = d0w5*(2.5D0*d0w2+d0w4)
      
       ENDDO

    ELSE
       
       DO i1 = 1, i0spcs 
          d1c1(i1) = 0.D0
          d1c2(i1) = 0.D0
          d1c3(i1) = 0.D0
          d1c4(i1) = 0.D0
       ENDDO
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2CALV_WV

  !C------------------------------------------------------------------
  !C 
  !C CALCULATION OF MASS SCALAR COEFFICIENTS
  !C
  !C                                                  CHECK 2013-11-03
  !C------------------------------------------------------------------
  
  SUBROUTINE T2CALV_MS
    
    USE T2CNST, ONLY:&
         d0vci2

    USE T2COMM, ONLY:&
         i0mscmx,d2ms,i0nmax1,i0dbg,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst
    INTEGER(i0ikind)::&
         i0cid,i1,i0vc

    !INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1msidv

    !ALLOCATE(i1msidv(1:i0mscmx))
    !i1msidv(1:i0mscmx)=0
    
    !C
    !C EQUATION OF PSI
    !C
    
    !C PSI'
    d2ms(1,i0nid) = d0mfcst*d0jm1
    !i1msidv(1)=1

    !C
    !C EQUATION OF I
    !C
    
    !C I
    d2ms(2,i0nid) = d0btcst*d0jm1
    !i1msidv(2)=2

    !C
    !C EQUATION OF Et
    !C
    
    !C Et
    d2ms(3,i0nid) = d0etcst*d0jm1*d0vci2
    !i1msidv(3)=3
    
    !C
    !C EQUATION OF Ep
    !C
    
    !C Ep
    d2ms(4,i0nid) = d0epcst*d0jm1*d0vci2
    !i1msidv(4)=4
    
    !C
    !C EQUATION OF Er
    !C

    !C Er
    d2ms(5,i0nid) = 0.D0
    !i1msidv(5)=5
    
    i0cid = 5
    
    DO i1 = 1, i0spcs
       
       i0vc = 8*i1-3
       !C
       !C EQUATION FOR N
       !C
       
       !C N
       d2ms(i0cid+1,i0nid) = d0nncst*d0jm1
       !i1msidv(i0cid+1) = i0vc+1
       
       !C
       !C EQUATION FOR Fr
       !C 
       
       !C Fr
       d2ms(i0cid+2,i0nid) = 0.D0
       !i1msidv(i0cid+2) = i0vc+2
       
       !C
       !C EQUATION FOR Fb
       !C 

       !C Fb
       d2ms(i0cid+3,i0nid) = d0fbcst*d0jm1*d0bb
       !i1msidv(i0cid+3) = i0vc+3
       
       !C
       !C EQUATION FOR Ft
       !C

       !C Ft
       d2ms(i0cid+4,i0nid) = d0ftcst*d0jm1
       !i1msidv(i0cid+4) = i0vc+4
       
       !C
       !C EQUATION FOR P
       !C

       !C P
       d2ms(i0cid+5,i0nid) = d0ppcst*d0jm1*1.5D0
       !i1msidv(i0cid+5) = i0vc+5
       
       !C
       !C EQUATION FOR Qr
       !C

       !C Qr
       d2ms(i0cid+6,i0nid) = 0.D0
       !i1msidv(i0cid+6) = i0vc+6       
       
       !C
       !C EQUATION FOR Qb
       !C

       !C Qb
       d2ms(i0cid+7,i0nid) = d0qbcst*d0jm1*d0bb
       !i1msidv(i0cid+7) = i0vc+7
       
       !C
       !C EQUATION FOR Qt
       !C

       !C Qt
       d2ms(i0cid+8,i0nid) = d0qtcst*d0jm1
       !i1msidv(i0cid+8) = i0vc+8
       
       i0cid = i0cid + 8
       
    ENDDO

    !OPEN(10,FILE='TEST_MSV.dat')
    !DO i1=1,i0mscmx
    !   WRITE(10,*)'i1=',i1,'I1MSIDC=',i1msidv(i1)
    !ENDDO
    !CLOSE(10)
    !DEALLOCATE(i1msidv)
    
    IF(i0cid.EQ.i0mscmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_MS_COEFF'
    STOP
    
  END SUBROUTINE T2CALV_MS

  !C------------------------------------------------------------------
  !C 
  !C CALCULATION OF ADVECTION VECTOR COEFFICIENTS
  !C
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2CALV_AV
    
    USE T2CNST, ONLY:&
         d0vci2
    
    USE T2COMM, ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1m, d1t, d1n, d1p, d1ni,d1pi,&
         d1ub,d1ur,d1ut,d1up,d1u2,d1qr,d1qb,d1qt,d1qp,&
         i0avcmx,d3av
    
    
    INTEGER(i0ikind)::&
         i0cid,i1,i0vc
    REAL(   i0rkind)::&
         d0m,d0t,d0n,d0p,d0ni,d0pi,&
         d0ur,d0up,d0ub,d0u2,&
         d0qr,d0qp,d0qt,d0qb,&
         d0x1,d0x2
    
    !INTEGER,DIMENSION(:),ALLOCATABLE::i1avidv

    !ALLOCATE(i1avidv(1:i0avcmx))
    !i1avidv(1:i0avcmx) = 0


    !C
    !C EQUATION FOR PSI
    !C
    
    !C PSI'
    d3av(1,1,i0nid) = d0mfcst*d0ugr
    d3av(2,1,i0nid) = d0mfcst*d0ugp
    !i1avidv(1) = 1

    !C
    !C EQUATION FOR I
    !C
    
    !C I
    d3av(1,2,i0nid) = d0btcst*d0ugr
    d3av(2,2,i0nid) = d0btcst*d0ugp
    !i1avidv(2) = 2

    !C
    !C EQUATION FOR Et
    !C
    
    !C PSI'
    d3av(1,3,i0nid) = -d0mfcst*d0jm1*d0jm2
    d3av(2,3,i0nid) = -d0mfcst*d0jm1*d0jm3
    !i1avidv(3) = 1
    
    !C Et
    d3av(1,4,i0nid) =  d0etcst*d0ugr*d0vci2
    d3av(2,4,i0nid) =  d0etcst*d0ugp*d0vci2
    !i1avidv(4) = 3

    !C
    !C EQUATION FOR Ep
    !C

    !C Ep
    d3av(1,5,i0nid) = d0epcst*d0ugr*d0vci2
    d3av(2,5,i0nid) = d0epcst*d0ugp*d0vci2
    !i1avidv(5) = 4

    !C
    !C EQUATION FOR Er
    !C

    !C Ep
    d3av(1,6,i0nid) = d0epcst*d0jm1*d0jm3
    d3av(2,6,i0nid) = d0epcst*d0jm1*d0jm4
    !i1avidv(6) = 4
    
    !C Er
    d3av(1,7,i0nid) = d0ercst*d0jm1*d0jm2
    d3av(2,7,i0nid) = d0ercst*d0jm1*d0jm3
    !i1avidv(7) = 5
    
    i0cid = 7
    
    DO i1 = 1, i0spcs
       
       d0m  = d1m( i1)
       d0t  = d1t( i1)
       d0n  = d1n( i1)
       d0ni = d1ni(i1)
       d0pi = d1pi(i1)
       d0ur = d1ur(i1)
       d0up = d1up(i1)
       d0ub = d1ub(i1)
       d0p  = d1p( i1)
       d0u2 = d1u2(i1)
       d0qr = d1qr(i1)
       d0qp = d1qp(i1)
       d0qb = d1qb(i1)
       d0qt = d1qt(i1)

       i0vc = 8*i1 - 3

       !C
       !C EQUATION FOR N
       !C
       
       !C N
       d3av(1,i0cid+1,i0nid) = d0nncst*(d0ugr+d0jm1*d0ur)
       d3av(2,i0cid+1,i0nid) = d0nncst*(d0ugp+d0jm1*d0up)
       !i1avidv(i0cid+1) = i0vc+1
       i0cid = i0cid + 1
       
       !C
       !C EQUATION FOR Fb
       !C
       
       !C Fb
       d3av(1,i0cid+1,i0nid) = d0fbcst* d0ug1*d0bb
       d3av(2,i0cid+1,i0nid) = d0fbcst*(d0ug2*d0bb + d0jm1*d0bp*d0ub)
       !i1avidv(i0cid+1) = i0vc+3
       
       !C P
       d3av(2,i0cid+2,i0nid) = d0ppcst*d0jm1*d0bp/d0m
       !i1avidv(i0cid+2) = i0vc+5
       
       i0cid = i0cid + 2
       
       !C
       !C EQUATION FOR Ft
       !C
       
       !C Ft
       d3av(1,i0cid+1,i0nid) = d0ftcst*(d0ugr + d0jm1*d0ur)
       d3av(2,i0cid+1,i0nid) = d0ftcst*(d0ugp + d0jm1*d0up)
       !i1avidv(i0cid+1) = i0vc + 4
       
       i0cid = i0cid + 1
       
       !C
       !C EQUATION FOR P
       !C
       
       d0x1 = 0.5D0*d0jm1*d0m*d0n*d0u2*d0pi
       d0x2 = d0jm1*d0pi
       
       !C P
       d3av(1,i0cid+1,i0nid)&
            = d0ppcst*(1.5D0*d0ugr - d0x1*d0ur + d0x2*d0qr)
       d3av(2,i0cid+1,i0nid)&
            = d0ppcst*(1.5D0*d0ugp - d0x1*d0up + d0x2*d0qp)
       !i1avidv(i0cid+1) = i0vc+5
       i0cid = i0cid + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       !C Fb
       d3av(2,i0cid+1,i0nid)&
            = d0fbcst*d0jm1*d0ni*(d0qb - 1.5D0*d0p*d0ub)*d0bp
       !i1avidv(i0cid+1) = i0vc+3

       !C P
       d3av(2,i0cid+2,i0nid) = d0ppcst*2.5D0*d0jm1*d0t*d0bp/d0m
       !i1avidv(i0cid+2) = i0vc+5
       
       !C Qb
       d3av(1,i0cid+3,i0nid) = d0qbcst* d0bb*d0ugr
       d3av(2,i0cid+3,i0nid) = d0qbcst*(d0bb*d0ugp+d0jm1*d0ub*d0bp)
       !i1avidv(i0cid+3) = i0vc+7
       
       i0cid = i0cid + 3
       
       !C
       !C EQUATION FOR Qt
       !C
       
       !C Ft
       d3av(1,i0cid+1,i0nid) = d0ftcst*d0jm1*d0ni*(d0qr-1.5D0*d0p*d0ur)
       d3av(2,i0cid+1,i0nid) = d0ftcst*d0jm1*d0ni*(d0qp-1.5D0*d0p*d0up)
       !i1avidv(i0cid+1) = i0vc+4

       !C Qt
       d3av(1,i0cid+2,i0nid) = d0qtcst*(d0ugr + d0jm1*d0ur)
       d3av(2,i0cid+2,i0nid) = d0qtcst*(d0ugp + d0jm1*d0up)
       !i1avidv(i0cid+2) = i0vc+8
       
       i0cid = i0cid + 2
       
    ENDDO
    
    !OPEN(10,FILE='TEST_AVV.dat')
    !DO i1=1,i0avcmx
    !   WRITE(10,*)'i1=',i1,'I1AVIDC=',i1avidv(i1)
    !ENDDO
    !CLOSE(10)
    !DEALLOCATE(i1avidv)

    IF(i0cid.EQ.i0avcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_AV_COEFF'
    STOP
    
  END SUBROUTINE T2CALV_AV
  
  !C
  !C
  !C  
  SUBROUTINE T2CALV_AT
    
    USE T2CNST
    
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1e,d1m,d1ur,d1up,d1ub,d2ug,&
         d1c1,d1c2,d1c3,d1c4,&
         d4at,i0atcmx
    
    INTEGER(i0ikind)::i1,i0cid,i0vc
    REAL(   i0rkind)::&
         d0m,d0ur,d0up,d0ub,d0c1,d0c2,d0c3,d0c4,&
         d0x1,d0x2
    !INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1atidv
    
    !ALLOCATE(i1atidv(1:i0atcmx))
    !i1atidv(1:i0atcmx) = 0

    i0cid = 0
    
    DO i1 = 1, i0spcs
       
       d0m  = d1m( i1)
       d0ur = d1ur(i1)
       d0up = d1up(i1)
       d0ub = d1ub(i1)
       d0c1 = d1c1(i1)
       d0c2 = d1c2(i1)
       d0c3 = d1c3(i1)
       d0c4 = d1c4(i1)
       
       i0vc = 8*i1 - 3
       
       !C
       !C EQUATION FOR Fb
       !C
       
       d0x1 = (2.D0*d0jm1*d0bp)/(3.D0*d0m)
       
       !C Fb (B)
       d4at(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c1*d0g4
       !i1atidv(i0cid+1) = i0vc+3

       !C Ft (B)
       d4at(2,2,i0cid+2,i0nid) = -d0ftcst*d0x1*d0c1*d0g5
       !i1atidv(i0cid+2) = i0vc+4

       !C Qb (B)
       d4at(2,2,i0cid+3,i0nid) = -d0qbcst*d0x1*d0c2*d0g4
       !i1atidv(i0cid+3) = i0vc+7

       !C Qt (B)
       d4at(2,2,i0cid+4,i0nid) = -d0qtcst*d0x1*d0c2*d0g5
       !i1atidv(i0cid+4) = i0vc+8

       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Ft
       !C
       
       d0x1 = d0jm1*d0bpb*d0btb/d0m
       
       !C Fb (B)
       d4at(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c1*d0g4
       !i1atidv(i0cid+1) = i0vc+3
       !C Ft (B)
       d4at(2,2,i0cid+2,i0nid) = -d0ftcst*d0x1*d0c1*d0g5
       !i1atidv(i0cid+2) = i0vc+4

       !C Qb (B)
       d4at(2,2,i0cid+3,i0nid) = -d0qbcst*d0x1*d0c2*d0g4
       !i1atidv(i0cid+3) = i0vc+7

       !C Qt (B)
       d4at(2,2,i0cid+4,i0nid) = -d0qtcst*d0x1*d0c2*d0g5
       !i1atidv(i0cid+4) = i0vc+8       

       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR P
       !C

       d0x1 = -d0jm1*d0ur/3.D0 
       d0x2 = -d0jm1*d0up/3.D0 + d0jm1*d0ub*d0bpb
       
       !C Fb
       d4at(1,2,i0cid+1,i0nid) = d0fbcst*d0x1*d0c1*d0g4
       d4at(2,2,i0cid+1,i0nid) = d0fbcst*d0x2*d0c1*d0g4
       !i1atidv(i0cid+1) = i0vc+3

       !C Ft
       d4at(1,2,i0cid+2,i0nid) = d0ftcst*d0x1*d0c1*d0g5
       d4at(2,2,i0cid+2,i0nid) = d0ftcst*d0x2*d0c1*d0g5 
       !i1atidv(i0cid+2) = i0vc+4
       
       !C Qb
       d4at(1,2,i0cid+3,i0nid) = d0qbcst*d0x1*d0c2*d0g4
       d4at(2,2,i0cid+3,i0nid) = d0qbcst*d0x2*d0c2*d0g4
       !i1atidv(i0cid+3) = i0vc+7

       !C Qt
       d4at(1,2,i0cid+4,i0nid) = d0qtcst*d0x1*d0c2*d0g5
       d4at(2,2,i0cid+4,i0nid) = d0qtcst*d0x2*d0c2*d0g5
       !i1atidv(i0cid+4) = i0vc+8
       
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Qb
       !C 
       
       d0x1 = 2.D0*d0jm1*d0bp/3.D0
       
       !C Fb
       d4at(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c3*d0g4
       !i1atidv(i0cid+1) = i0vc+3

       !C Ft
       d4at(2,2,i0cid+2,i0nid) = -d0ftcst*d0x1*d0c3*d0g5
       !i1atidv(i0cid+2) = i0vc+4

       !C Qb
       d4at(2,2,i0cid+3,i0nid) = -d0qbcst*d0x1*d0c4*d0g4
       !i1atidv(i0cid+3) = i0vc+7

       !C Qt
       d4at(2,2,i0cid+4,i0nid) = -d0qtcst*d0x1*d0c4*d0g5
       !i1atidv(i0cid+4) = i0vc+8
       
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Qt
       !C

       d0x1 = d0jm1*d0bpb*d0btb

       !C Fb
       d4at(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c3*d0g4
       !i1atidv(i0cid+1) = i0vc+3

       !C Ft
       d4at(2,2,i0cid+2,i0nid) = -d0ftcst*d0x1*d0c3*d0g5
       !i1atidv(i0cid+2) = i0vc+4

       !C Qb
       d4at(2,2,i0cid+3,i0nid) = -d0qbcst*d0x1*d0c4*d0g4
       !i1atidv(i0cid+3) = i0vc+7

       !C Qt
       d4at(2,2,i0cid+4,i0nid) = -d0qtcst*d0x1*d0c4*d0g5
       !i1atidv(i0cid+4) = i0vc+8
       
       i0cid = i0cid + 4
       
    ENDDO
    
    !OPEN(10,FILE='TEST_ATV.dat')
    !DO i1=1,i0atcmx
    !   WRITE(10,*)'i1=',i1,'I1ATIDC=',i1atidv(i1)
    !ENDDO
    !CLOSE(10)
    !DEALLOCATE(i1atidv)
    
    IF(i0cid.EQ.i0atcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_AT_COEFF'
    
    STOP
    
  END SUBROUTINE T2CALV_AT

  !C
  !C
  !C  
  SUBROUTINE T2CALV_DT
    
    USE T2CNST
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1m,d1p,d1ur,d1up,d1ub,d1qb,d1qbp,&
         d1c1,d1c2,d1c3,d1c4,&
         i0dtcmx,d4dt
    
    INTEGER(i0ikind)::i1,i0cid,i0vc
    
    REAL(   i0rkind)::&
         d0m,d0p,d0ur,d0up,d0ub,d0qb,d0qbp,&
         d0c1,d0c2,d0c3,d0c4,d0x1,d0x2
    
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1dtidv

    !ALLOCATE(i1dtidv(1:i0dtcmx))
    !i1dtidv(1:i0dtcmx) = 0
    
    i0cid = 0
    
    DO i1 = 1, i0spcs
       
       d0m   = d1m(  i1)
       d0p   = d1p(  i1)
       d0ur  = d1ur( i1)
       d0up  = d1up( i1)
       d0ub  = d1ub( i1)
       d0qb  = d1qb( i1)
       d0qbp = d1qbp(i1)
       d0c1  = d1c1( i1)
       d0c2  = d1c2( i1)
       d0c3  = d1c3( i1)
       d0c4  = d1c4( i1)
       
       i0vc = 8*i1 - 3
       
       !C   
       !C EQUATION FOR Fb
       !C
       d0x1 =  (2.D0*d0jm1*d0bp*d0bpb)/(3.D0*d0m)
       
       !C N
       d4dt(2,2,i0cid+1,i0nid) = - d0nncst*d0x1*d0c1*d0ub
       !C Fb
       d4dt(2,2,i0cid+2,i0nid) =   d0fbcst*d0x1*d0c1      
       !C P
       d4dt(2,2,i0cid+3,i0nid) = - d0ppcst*d0x1*d0c2*d0qbp
       !C Qb
       d4dt(2,2,i0cid+4,i0nid) =   d0qbcst*d0x1*d0c2      
       
       !i1dtidv(i0cid+1) = i0vc+1
       !i1dtidv(i0cid+2) = i0vc+3
       !i1dtidv(i0cid+3) = i0vc+5
       !i1dtidv(i0cid+4) = i0vc+7

       i0cid = i0cid + 4

       !C
       !C  EQUATION FOR Ft
       !C 
       d0x1 = d0jm1*d0btb*d0bpb*d0bpb/d0m
       
       !C N
       d4dt(2,2,i0cid+1,i0nid) = -d0nncst*d0x1*d0c1*d0ub
       !C Fb
       d4dt(2,2,i0cid+2,i0nid) =  d0fbcst*d0x1*d0c1
       !C P
       d4dt(2,2,i0cid+3,i0nid) = -d0ppcst*d0x1*d0c2*d0qbp
       !C Qb
       d4dt(2,2,i0cid+4,i0nid) =  d0qbcst*d0x1*d0c2

       !i1dtidv(i0cid+1) = i0vc+1
       !i1dtidv(i0cid+2) = i0vc+3
       !i1dtidv(i0cid+3) = i0vc+5
       !i1dtidv(i0cid+4) = i0vc+7
       
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR P
       !C
       
       d0x1 = d0jm1*d0bpb*d0ur/3.D0 
       d0x2 = d0jm1*d0bpb*d0up/3.D0 - d0jm1*d0bpb*d0bpb*d0ub       
       !C N
       d4dt(1,2,i0cid+1,i0nid) = -d0nncst*d0x1*d0c1*d0ub
       d4dt(2,2,i0cid+1,i0nid) = -d0nncst*d0x2*d0c1*d0ub
       !C Fb
       d4dt(1,2,i0cid+2,i0nid) =  d0fbcst*d0x1*d0c1
       d4dt(2,2,i0cid+2,i0nid) =  d0fbcst*d0x2*d0c1
       !C P
       d4dt(1,2,i0cid+3,i0nid) = -d0ppcst*d0x1*d0c2*d0qbp
       d4dt(2,2,i0cid+3,i0nid) = -d0ppcst*d0x2*d0c2*d0qbp
       !C Qb
       d4dt(1,2,i0cid+4,i0nid) =  d0qbcst*d0x1*d0c2
       d4dt(2,2,i0cid+4,i0nid) =  d0qbcst*d0x2*d0c2

       !i1dtidv(i0cid+1) = i0vc+1
       !i1dtidv(i0cid+2) = i0vc+3
       !i1dtidv(i0cid+3) = i0vc+5
       !i1dtidv(i0cid+4) = i0vc+7

       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Qb
       !C
       
       d0x1 = 2.D0*d0jm1*d0bp*d0bpb/3.D0
       
       !C N
       d4dt(2,2,i0cid+1,i0nid) = - d0nncst*d0x1*d0c3*d0ub
       !C Fb
       d4dt(2,2,i0cid+2,i0nid) =   d0fbcst*d0x1*d0c3
       !C P
       d4dt(2,2,i0cid+3,i0nid) = - d0ppcst*d0x1*d0c4*d0qbp
       !C Qb
       d4dt(2,2,i0cid+4,i0nid) =   d0qbcst*d0x1*d0c4

       !i1dtidv(i0cid+1) = i0vc+1
       !i1dtidv(i0cid+2) = i0vc+3
       !i1dtidv(i0cid+3) = i0vc+5
       !i1dtidv(i0cid+4) = i0vc+7
       
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Qt
       !C

       d0x1 = d0jm1*d0btb*d0bpb*d0bpb

       !C  
       d4dt(2,2,i0cid+1,i0nid) = - d0nncst*d0x1*d0c3*d0ub
       !C
       d4dt(2,2,i0cid+2,i0nid) =   d0fbcst*d0x1*d0c3
       !C
       d4dt(2,2,i0cid+3,i0nid) = - d0ppcst*d0x1*d0c4*d0qbp
       !C
       d4dt(2,2,i0cid+4,i0nid) =   d0qbcst*d0x1*d0c4

       !i1dtidv(i0cid+1) = i0vc+1
       !i1dtidv(i0cid+2) = i0vc+3
       !i1dtidv(i0cid+3) = i0vc+5
       !i1dtidv(i0cid+4) = i0vc+7
       
       i0cid = i0cid + 4

    ENDDO

    !OPEN(10,FILE='TEST_DTV.dat')
    !DO i1=1,i0dtcmx
    !   WRITE(10,*)'i1=',i1,'I1DTIDC=',i1dtidv(i1)
    !ENDDO
    !CLOSE(10)
    !DEALLOCATE(i1dtidv)
    IF(i0cid.EQ.i0dtcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_DT_COEFF'
    STOP
    
  END SUBROUTINE T2CALV_DT

  !C
  !C
  !C
  SUBROUTINE T2CALV_GV

    USE T2CNST
    USE T2COMM, ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1m,d1p,d1t,d1ur,d1up,d2f1,d2f2,d2f3,d2f4,i1pdn1,&
         i0gvcmx,d3gv
    
    INTEGER(i0ikind)::&
         i0cid,i1,i0vc
    REAL(   i0rkind)::&
         d0m,d0t,d0ur,d0up,d0x1    
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1gvidv
    
    !ALLOCATE(i1gvidv(1:i0gvcmx))
    !i1gvidv(1:i0gvcmx)         = 0    
    
    !C
    !C EQUATION FOR PSI
    !C
    
    !C Et
    d3gv(1,1,i0nid) = - d0etcst*d0jm1
    !i1gvidv(1) = 3
    
    !C
    !C EQUATION FOR I
    !C
    
    d0x1 = 1.D0/d0jm5
    
    !C Ep
    d3gv(1,2,i0nid) =  d0epcst*d0x1
    !i1gvidv(2) = 4
    
    !C Er
    d3gv(2,3,i0nid) = -d0ercst*d0x1
    !i1gvidv(3) = 5
    
    !C 
    !C EQUATION FOR Ep
    !C
    
    !C I
    d3gv(1,4,i0nid) =  d0btcst*(d0jm1**2)*d0jm2*d0jm5
    !i1gvidv(4) = 2
    
    i0cid = 4
    
    DO i1 = 1,i0spcs
       
       d0m  = d1m( i1)
       d0t  = d1t( i1)
       d0ur = d1ur(i1)
       d0up = d1up(i1)
       
       i0vc = 8*i1 - 3
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C P
       d3gv(1,i0cid+1,i0nid) =  d0ppcst*d0jm1*d0jm2/d0m
       d3gv(2,i0cid+1,i0nid) =  d0ppcst*d0jm1*d0jm3/d0m
       !i1gvidv(i0cid+1) = i0vc + 5
       i0cid = i0cid + 1
       
       !C 
       !C EQUATION FOR P
       !C

       !C P
       d3gv(1,i0cid+1,i0nid) = -d0ppcst*d0jm1*d0ur
       d3gv(2,i0cid+1,i0nid) = -d0ppcst*d0jm1*d0up
       !i1gvidv(i0cid+1) = i0vc + 5
       
       i0cid = i0cid + 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       d0x1 = 2.5D0*d0jm1*d0t/d0m
       
       !C N
       d3gv(1,i0cid+1,i0nid) = -d0nncst     *d0x1*d0jm2*d0t
       d3gv(2,i0cid+1,i0nid) = -d0nncst     *d0x1*d0jm3*d0t
       !i1gvidv(i0cid+1) = i0vc + 1
       
       !C P
       d3gv(1,i0cid+2,i0nid) =  d0ppcst*2.D0*d0x1*d0jm2
       d3gv(2,i0cid+2,i0nid) =  d0ppcst*2.D0*d0x1*d0jm3
       !i1gvidv(i0cid+2) = i0vc + 5
       
       i0cid = i0cid + 2
       
    ENDDO
    
    !OPEN(10,FILE='TEST_GVV.dat')
    !DO i1=1,i0gvcmx
    !   WRITE(10,*)'i1=',i1,'I1GVIDC=',i1gvidv(i1)
    !ENDDO
    !CLOSE(10)
    !DEALLOCATE(i1gvidv)
    
    IF(i0cid.EQ.i0gvcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_GV_COEFF'
    STOP
    
    RETURN
  
  END SUBROUTINE T2CALV_GV
  

  !C
  !C
  !C
  SUBROUTINE T2CALV_GT
    
    USE T2CNST
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
       d1m, d1ut,d1ub,d1qbp,&
       d1c1,d1c2,d1c3,d1c4,&
       i0gtcmx,d4gt 
    
    INTEGER(i0ikind)::i1,i0cid,i0vc
    REAL(   i0rkind)::&
         d0m,d0ub,d0ut,d0qbp,&
         d0c1,d0c2,d0c3,d0c4,d0x1,d0x2
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1gtidv

    !ALLOCATE(i1gtidv(1:i0gtcmx))
    !i1gtidv(1:i0gtcmx) = 0
    
    i0cid = 0
    
    DO i1 = 1, i0spcs

       d0m   = d1m(  i1)
       d0ut  = d1ut( i1)
       d0ub  = d1ub( i1)
       d0qbp = d1qbp(i1)
       d0c1  = d1c1( i1)
       d0c2  = d1c2( i1)
       d0c3  = d1c3( i1)
       d0c4  = d1c4( i1)
       
       i0vc = 8*i1 - 3
       
       !C
       !C EQUATION FOR Fb
       !C
       
       d0x1 = d0jm1*d0bpb*d0bpb/d0m

       !C N
       d4gt(2,2,i0cid+1,i0nid) = - d0nncst*d0x1*d0c1*d0ub
       !C Fb
       d4gt(2,2,i0cid+2,i0nid) =   d0fbcst*d0x1*d0c1     
       !C P
       d4gt(2,2,i0cid+3,i0nid) = - d0ppcst*d0x1*d0c2*d0qbp
       !C Qb
       d4gt(2,2,i0cid+4,i0nid) =   d0qbcst*d0x1*d0c2

       !i1gtidv(i0cid+1) = i0vc+1
       !i1gtidv(i0cid+2) = i0vc+3
       !i1gtidv(i0cid+3) = i0vc+5
       !i1gtidv(i0cid+4) = i0vc+7
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR P
       !C
       
       d0x1 = d0jm1*d0bpb*(d0g4*d0ub+d0g5*d0ut)
       d0x2 = d0jm1*d0bpb*d0bpb
       
       !C N (B)
       d4gt(2,2,i0cid+1,i0nid) = -d0nncst*d0x1*d0c1*d0ub
       !C N (Ub)
       d4gt(2,2,i0cid+2,i0nid) =  d0nncst*d0x2*d0c1*d0ub

       !C Fb (B)
       d4gt(2,2,i0cid+3,i0nid) =  d0fbcst*d0x1*d0c1
       !C Fb (Ub)
       d4gt(2,2,i0cid+4,i0nid) = -d0fbcst*d0x2*d0c1

       !C P (B)
       d4gt(2,2,i0cid+5,i0nid) = -d0ppcst*d0x1*d0c2*d0qbp
       !C P (Ub)
       d4gt(2,2,i0cid+6,i0nid) =  d0ppcst*d0x2*d0c2*d0qbp

       !C Qb (B)
       d4gt(2,2,i0cid+7,i0nid) =  d0qbcst*d0x1*d0c2
       !C Qb (Ub)
       d4gt(2,2,i0cid+8,i0nid) = -d0qbcst*d0x2*d0c2

       !i1gtidv(i0cid+1) = i0vc+1
       !i1gtidv(i0cid+2) = i0vc+1
       !i1gtidv(i0cid+3) = i0vc+3
       !i1gtidv(i0cid+4) = i0vc+3
       !i1gtidv(i0cid+5) = i0vc+5
       !i1gtidv(i0cid+6) = i0vc+5
       !i1gtidv(i0cid+7) = i0vc+7
       !i1gtidv(i0cid+8) = i0vc+7
       i0cid = i0cid + 8

       !C
       !C EQUATION FOR Qb
       !C
       
       d0x1 = d0jm1*d0bpb*d0bpb

       !C N  (B)
       d4gt(2,2,i0cid+1,i0nid) = -d0nncst*d0x1*d0c3*d0ub
       !C Fb (B)
       d4gt(2,2,i0cid+2,i0nid) =  d0fbcst*d0x1*d0c3
       !C P  (B)
       d4gt(2,2,i0cid+3,i0nid) = -d0ppcst*d0x1*d0c4*d0qbp
       !C Qb (B)
       d4gt(2,2,i0cid+4,i0nid) =  d0qbcst*d0x1*d0c4
       
       !i1gtidv(i0cid+1) = i0vc+1
       !i1gtidv(i0cid+2) = i0vc+3
       !i1gtidv(i0cid+3) = i0vc+5
       !i1gtidv(i0cid+4) = i0vc+7
       i0cid = i0cid + 4
       
    ENDDO
    
    !OPEN(10,FILE='TEST_GTV.dat')
    !DO i1=1,i0gtcmx
    !   WRITE(10,*)'i1=',i1,'I1GTIDC=',i1gtidv(i1)
    !ENDDO
    !CLOSE(10)
    !DEALLOCATE(i1gtidv)
    
    IF(i0cid.EQ.i0gtcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_GT_COEFF'
    STOP
    
    RETURN

  END SUBROUTINE T2CALV_GT
  
  !C
  !C
  !C
  SUBROUTINE T2CALV_ES
    
    USE T2CNST,ONLY:&
         d0eps0,d0rmu0
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1e, d1m, d1n, d1p,&
         d1hx,d2f1,d2f2,d2f3,d2f4,&
         i0escmx,d2es
    
    INTEGER(i0ikind)::i1,j1,i0cid,i0vc
    
    REAL(   i0rkind)::&
         d0m,d0e,d0n,d0p,&
         d0f1,d0f2,d0f3,d0f4,d0hx,&
         d0x1,d0x2
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1esidv
    
    !ALLOCATE(i1esidv(1:i0escmx))
    !i1esidv(1:i0escmx) = 0
    
    i0cid = 0
    
    !C
    !C EQUATION FOR Et
    !C 
    DO j1 = 1, i0spcs
       
       d0e  = d1e(j1)
       i0vc = 8*j1-3
       
       !C Ft
       d2es(i0cid+1,i0nid) = d0ftcst*d0jm1*d0rmu0*d0e
       !i1esidv(i0cid+1) = i0vc+4
       i0cid = i0cid + 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    DO j1 =1, i0spcs
       
       d0e  = d1e(j1)
       d0x1 = d0jm1*d0rmu0*d0e*d0bpi
       i0vc = 8*j1-3
       
       !C Fb
       d2es(i0cid+1,i0nid) =   d0fbcst*d0x1*d0bb
       !i1esidv(i0cid+1) = i0vc+3 
       
       !C Ft
       d2es(i0cid+2,i0nid) = - d0ftcst*d0x1*d0bt*d0jm5
       !i1esidv(i0cid+2) = i0vc+4 
       
       i0cid = i0cid + 2
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C
    DO j1 = 1, i0spcs
       
       d0e  = d1e(j1)
       i0vc = 8*j1-3
       
       !C N
       d2es(i0cid+1,i0nid) = - d0nncst*d0jm1*d0e/d0eps0
       !i1esidv(i0cid+1) = i0vc+1 
       i0cid = i0cid + 1
       
    ENDDO
    
    DO i1 = 1, i0spcs
       
       d0n  = d1n(i1)
       d0p  = d1p(i1)
       d0m  = d1m(i1)
       d0e  = d1e(i1)
       
       !C
       !C EQUATION FOR Fr
       !C
       
       d0x1 = d0e/d0m
       i0vc = 8*i1-3
       
       !C Ep
       d2es(i0cid+1,i0nid) = -d0epcst*d0x1*d0jm1*d0n*d0jm3
       !i1esidv(i0cid+1) = 4
       
       !C Er
       d2es(i0cid+2,i0nid) = -d0ercst*d0x1*d0jm1*d0n*d0jm2
       !i1esidv(i0cid+2) = 5
       
       !C Fb
       d2es(i0cid+3,i0nid) = -d0fbcst*d0x1*d0bt*d0bb*d0bpi
       !i1esidv(i0cid+3) = i0vc+3 
       
       !C Ft
       d2es(i0cid+4,i0nid) =  d0ftcst*d0x1*d0bb*d0bb*d0bpi
       !i1esidv(i0cid+4) = i0vc+4 
       
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Fb
       !C
       
       DO j1 = 1, i0spcs
          
          d0x1 = d0jm1*d0e*d0n/d0m 
          d0x2 = d0jm1*d0bb/d0m
          i0vc = 8*j1 - 3
          d0f1 = d2f1(i1,j1)
          d0f2 = d2f2(i1,j1)
          
          IF(i1.EQ.j1)THEN
             
             !C Et
             d2es(i0cid+1,i0nid) = - d0etcst*d0x1*d0bt*d0jm5
             !i1esidv(i0cid+1) = 3
             
             !C Ep
             d2es(i0cid+2,i0nid) = - d0epcst*d0x1*d0bp
             !i1esidv(i0cid+2) = 4
             
             !C Fb
             d2es(i0cid+3,i0nid) = - d0fbcst*d0x2*d0f1
             !i1esidv(i0cid+3) = i0vc+3 
             
             !C Qb
             d2es(i0cid+4,i0nid) = - d0qbcst*d0x2*d0f2
             !i1esidv(i0cid+4) = i0vc+7 
             
             i0cid = i0cid + 4
             
          ELSE
             
             !C Fb
             d2es(i0cid+1,i0nid) = - d0fbcst*d0x2*d0f1
             !i1esidv(i0cid+1) = i0vc+3 
             
             !C Qb
             d2es(i0cid+2,i0nid) = - d0qbcst*d0x2*d0f2            
             !i1esidv(i0cid+2) = i0vc+7 
             
             i0cid = i0cid + 2
             
          ENDIF
       ENDDO

       !C
       !C EQUATION FOR Ft
       !C
       
       DO j1 = 1, i0spcs

          d0x1 = d0jm1*d0e/d0m
          d0x2 = d0jm1/d0m
          i0vc = 8*j1 - 3
          d0f1 = d2f1(i1,j1)
          d0f2 = d2f2(i1,j1)
          
          IF(i1.EQ.j1)THEN
             
             !C LORENTZ FORCE
             
             !C Et
             d2es(i0cid+1,i0nid) = -d0etcst*d0x1*d0n
             !i1esidv(i0cid+1) = 3

             !C Fr
             d2es(i0cid+2,i0nid) = -d0frcst*d0x1*d0jm1*d0bp
             !i1esidv(i0cid+2) = i0vc+2

             !C FRICTION FORCE

             !C Ft
             d2es(i0cid+3,i0nid) = -d0ftcst*d0x2*d0f1
             !i1esidv(i0cid+3) = i0vc+4 
             
             !C Qt
             d2es(i0cid+4,i0nid) = -d0qtcst*d0x2*d0f2
             !i1esidv(i0cid+4) = i0vc+8 
             
             i0cid = i0cid + 4
             
          ELSE
             !C FRICTION FORCE

             !C Ft
             d2es(i0cid+1,i0nid) = -d0ftcst*d0x2*d0f1
             !i1esidv(i0cid+1) = i0vc+4
             
             !C Qt
             d2es(i0cid+2,i0nid) = -d0qtcst*d0x2*d0f2
             !i1esidv(i0cid+2) = i0vc+8 
             
             i0cid = i0cid + 2
             
          ENDIF
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       
       !C P
       i0vc = 8*i1 - 3
       
       d0hx = d1hx(i1)
       d2es(i0cid+1,i0nid) = d0ppcst*d0jm1*d0hx
       !i1esidv(i0cid+1) = i0vc+5 
       i0cid = i0cid + 1

       !C
       !C EQUATION FOR Qr
       !C
       
       i0vc = 8*i1 - 3
       d0x1 = d0e/d0m
       d0x2 = 2.5D0*d0jm1*d0x1*d0p
       
       !C Ep
       d2es(i0cid+1,i0nid) = -d0epcst*d0x2*d0jm3
       !i1esidv(i0cid+1) = 4

       !C Er
       d2es(i0cid+2,i0nid) = -d0ercst*d0x2*d0jm2
       !i1esidv(i0cid+2) = 5 

       !C Qb
       d2es(i0cid+3,i0nid) = -d0qbcst*d0x1*d0bt*d0bb*d0bpi
       !i1esidv(i0cid+3) = i0vc+7 
       
       !C Qt
       d2es(i0cid+4,i0nid) =  d0qtcst*d0x1*d0bb*d0bb*d0bpi
       !i1esidv(i0cid+4) = i0vc+8               
       
       i0cid = i0cid + 4
       
       !C
       !C EQUATION FOR Qb
       !C
       
       DO j1 = 1, i0spcs
          
          d0x1 = 2.5D0*d0jm1*d0e*d0p/d0m
          d0x2 = d0jm1*d0bb
          d0f3 = d2f3(i1,j1)
          d0f4 = d2f4(i1,j1)
          i0vc = 8*j1 - 3

          IF(i1.EQ.j1)THEN
             
             !C LORENTZ FORCE
             
             !C Et
             d2es(i0cid+1,i0nid) = -d0etcst*d0x1*d0bt*d0jm5
             !i1esidv(i0cid+1) = 3

             !C Ep
             d2es(i0cid+2,i0nid) = -d0epcst*d0x1*d0bp
             !i1esidv(i0cid+2) = 4

             !C FRICTION FORCE
                      
             !C Fb  
             d2es(i0cid+3,i0nid) = -d0fbcst*d0x2*d0f3
             !i1esidv(i0cid+3) = i0vc+3

             !C Qb
             d2es(i0cid+4,i0nid) = -d0qbcst*d0x2*d0f4
             !i1esidv(i0cid+4) = i0vc+7
             
             i0cid = i0cid + 4
             
          ELSE
             
             !C FRICTION FORCE
             
             !C Fb
             d2es(i0cid+1,i0nid) = -d0fbcst*d0x2*d0f3
             !i1esidv(i0cid+1) = i0vc+3 
             
             !C Qb
             d2es(i0cid+2,i0nid) = -d0qbcst*d0x2*d0f4            
             !i1esidv(i0cid+2) = i0vc+7
             
             i0cid = i0cid + 2
             
          END IF
       ENDDO
       
       !C
       !C EQUATION FOR Qt
       !C
       
       DO j1 = 1, i0spcs

          d0x1 = d0jm1*d0e/d0m
          d0f3 = d2f3(i1,j1)
          d0f4 = d2f4(i1,j1)
          i0vc = 8*j1 - 3

          IF(i1.EQ.j1)THEN
                          
             !C Et
             d2es(i0cid+1,i0nid) = -d0etcst*d0x1*2.5D0*d0p
             !i1esidv(i0cid+1) = 3

             !C Ft
             d2es(i0cid+2,i0nid) = -d0ftcst*d0jm1*d0f3
             !i1esidv(i0cid+2) = i0vc+4

             !C Qr
             d2es(i0cid+3,i0nid) = -d0qrcst*d0x1*d0jm1*d0bp
             !i1esidv(i0cid+3) = i0vc+6

             !C Qt
             d2es(i0cid+4,i0nid) = -d0qtcst*d0jm1*d0f4
             !i1esidv(i0cid+4) = i0vc+8
             
             i0cid = i0cid + 4
             
          ELSE
             
             !C Ft
             d2es(i0cid+1,i0nid) = -d0ftcst*d0jm1*d0f3
             !i1esidv(i0cid+1) = i0vc+4 
             
             !C Qt
             d2es(i0cid+2,i0nid) = -d0qtcst*d0jm1*d0f4
             !i1esidv(i0cid+2) = i0vc+8              
             
             i0cid = i0cid + 2
             
          ENDIF
       ENDDO
    ENDDO
    
    !OPEN(10,FILE='TEST_ESV.dat')
    !DO i1=1,i0escmx
    !   WRITE(10,*)'i1=',i1,'I1ESIDC=',i1esidv(i1)
    !ENDDO
    !CLOSE(10)
    
    IF(i0cid.EQ.i0escmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_ES_COEFF'
    STOP 
    
  END SUBROUTINE T2CALV_ES
  
  
  SUBROUTINE T2CALV_EV
    
    USE T2CNST
    USE T2COMM,ONLY:&
         i0spcs,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1e, d1m, d1n, d1ut,d1ub,d1p, d1qb,d1qt,d1qbp,d1ni,&
         d1c1,d1c2,&
         i0evcmx,d3ev
    
    INTEGER(i0ikind)::&
         i1,i0cid,i0vc
    REAL(   i0rkind)::&
         d0e, d0m, d0n, d0ut,d0ub,d0p ,d0qb,d0qt,d0qbp,d0ni,&
         d0c1,d0c2,d0x1,d0x2,d0x3,d0x4
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1evidv

    
    !ALLOCATE(i1evidv(1:i0evcmx))
    !i1evidv(1:i0evcmx) = 0

    !C
    !C EQUATION FOR Et
    !C
    
    !C PSI' (ln R)
    d3ev(1,1,i0nid) = d0mfcst*2.D0*d0jm1*d0jm2
    d3ev(2,1,i0nid) = d0mfcst*2.D0*d0jm1*d0jm3
    !i1evidv(1) = 1
    
    i0cid = 1
    
    DO i1 = 1,i0spcs
       
       d0e   = d1e(  i1)
       d0m   = d1m(  i1)
       d0n   = d1n(  i1)
       d0ni  = d1ni( i1)
       d0ut  = d1ut( i1)
       d0ub  = d1ub( i1)
       d0p   = d1p(  i1)
       d0qb  = d1qb( i1)
       d0qbp = d1qbp(i1)
       d0qt  = d1qt( i1)
       d0c1  = d1c1( i1)
       d0c2  = d1c2( i1)
       
       i0vc = 8*i1 - 3
       
       !C
       !C EQUATION FOR Fb
       !C
       
       !C Fb (B)
       d3ev(2,i0cid+1,i0nid) = -d0fbcst*d0jm1*d0ub*d0bpb
       
       !i1evidv(i0cid+1) = i0vc+3
       i0cid = i0cid + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       d0x1 = d0jm1*d0e*d0bp/d0m
       d0x2 = d0jm1*d0e*d0bt*d0jm5/d0m
       d0x3 = d0c1*d0n*(d0ub*d0g4+d0ut*d0g5)&
            + d0c2*(d0qb*d0g4+d0qt*d0g5)
       d0x4 = d0qb-1.5D0*d0p*d0ub
       
       !C Et (B)
       d3ev(2,i0cid+ 1,i0nid) = -d0etcst*d0x2*d0x3
       !i1evidv(i0cid+ 1) = 3
       
       !C Et (N)
       d3ev(2,i0cid+ 2,i0nid) = -d0etcst*d0x2*d0c1*d0bpb*d0ub
       !i1evidv(i0cid+ 2) = 3
       
       !C Et (Fb)
       d3ev(2,i0cid+ 3,i0nid) =  d0etcst*d0x2*d0c1*d0bpb
       !i1evidv(i0cid+ 3) = 3
       
       !C Et (P)
       d3ev(2,i0cid+ 4,i0nid) = -d0etcst*d0x2*d0c2*d0bpb*d0qbp
       !i1evidv(i0cid+ 4) = 3
       
       !C Et (Qb)
       d3ev(2,i0cid+ 5,i0nid) =  d0etcst*d0x2*d0c2*d0bpb
       !i1evidv(i0cid+ 5) = 3
       
       !C Ep (B)
       d3ev(2,i0cid+ 6,i0nid) = -d0epcst*d0x1*d0x3
       !i1evidv(i0cid+ 6) = 4
       
       !C Ep (N)
       d3ev(2,i0cid+ 7,i0nid) = -d0epcst*d0x1*d0c1*d0bpb*d0ub
       !i1evidv(i0cid+ 7) = 4
       
       !C Ep (Fb)
       d3ev(2,i0cid+ 8,i0nid) =  d0epcst*d0x1*d0c1*d0bpb
       !i1evidv(i0cid+ 8) = 4
       
       !C Ep (P)
       d3ev(2,i0cid+ 9,i0nid) = -d0epcst*d0x1*d0c2*d0bpb*d0qbp
       !i1evidv(i0cid+ 9) = 4
       
       !C Ep (Qb)
       d3ev(2,i0cid+10,i0nid) =  d0epcst*d0x1*d0c2*d0bpb
       !i1evidv(i0cid+10) = 4
       
       !C Fb (B)
       d3ev(2,i0cid+11,i0nid) = -d0fbcst*d0jm1*d0ni*d0x4*d0bpb
       !i1evidv(i0cid+11) = i0vc+3
       
       !C Qb (B) 
       d3ev(2,i0cid+12,i0nid) = -d0qbcst*d0jm1*d0ub*d0bpb
       !i1evidv(i0cid+12) = i0vc+7
       
       i0cid = i0cid + 12

       !C
       !C EQUATION FOR Qt
       !C
       
       d0x1 = d0jm1*d0e*d0bpb*d0btb/d0m
       d0x2 = d0jm1*d0e*((d0btb**2)*d0jm5 -1.D0/3.D0)/d0m
       d0x3 = d0c1*d0n*(d0ub*d0g4+d0ut*d0g5)&
            + d0c2*(d0qb*d0g4+d0qt*d0g5)

       !C Et (B)
       d3ev(2,i0cid+ 1,i0nid) = -d0etcst*d0x2*d0x3
       !i1evidv(i0cid+ 1) = 3
       
       !C Et (N)
       d3ev(2,i0cid+ 2,i0nid) = -d0etcst*d0x2*d0c1*d0bpb*d0ub
       !i1evidv(i0cid+ 2) = 3
       
       !C Et (Fb)
       d3ev(2,i0cid+ 3,i0nid) =  d0etcst*d0x2*d0c1*d0bpb
       !i1evidv(i0cid+ 3) = 3

       !C Et (P)
       d3ev(2,i0cid+ 4,i0nid) = -d0etcst*d0x2*d0c2*d0bpb*d0qbp
       !i1evidv(i0cid+ 4) = 3

       !C Et (Qb)
       d3ev(2,i0cid+ 5,i0nid) =  d0etcst*d0x2*d0c2*d0bpb
       !i1evidv(i0cid+ 5) = 3

       !C Ep (B)
       d3ev(2,i0cid+ 6,i0nid) = -d0epcst*d0x1*d0x3
       !i1evidv(i0cid+ 6) = 4

       !C Ep (N)
       d3ev(2,i0cid+ 7,i0nid) = -d0epcst*d0x1*d0c1*d0bpb*d0ub
       !i1evidv(i0cid+ 7) = 4
       
       !C Ep (Fb)
       d3ev(2,i0cid+ 8,i0nid) =  d0epcst*d0x1*d0c1*d0bpb
       !i1evidv(i0cid+ 8) = 4
       
       !C Ep (P)
       d3ev(2,i0cid+ 9,i0nid) = -d0epcst*d0x1*d0c2*d0bpb*d0qbp
       !i1evidv(i0cid+ 9) = 4
       
       !C Ep (Qb)
       d3ev(2,i0cid+10,i0nid) =  d0epcst*d0x1*d0c2*d0bpb
       !i1evidv(i0cid+10) = 4
       
       i0cid = i0cid + 10
       
    ENDDO
    
    ! OPEN(10,FILE='TEST_EVV.dat')
    ! DO i1=1,i0evcmx
    !    WRITE(10,*)'i1=',i1,'I1EVIDC=',i1evidv(i1)
    ! ENDDO
    ! CLOSE(10)
    ! DEALLOCATE(i1evidv)
    IF(i0cid.EQ.i0evcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_EV_COEFF'
    STOP
    
  END SUBROUTINE T2CALV_EV
  

  SUBROUTINE T2CALV_ET
    
    USE T2CNST
    USE T2COMM,ONLY:&
         i0spcs,&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         d1m, d1ut,d1ub,d1c1,d1c2,d1c3,d1c4,&
         i0etcmx,d4et
    
    INTEGER(i0ikind)::&
         i1,i0cid,i0vc
    
    REAL(   i0rkind)::&
         d0m, d0ut,d0ub,&
         d0c1,d0c2,d0c3,d0c4,d0x1,d0x2
    
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1etidv
    ALLOCATE(i1etidv(1:i0etcmx))
    i1etidv(1:i0etcmx) = 0

    i0cid = 0
    
    DO i1 = 1,i0spcs
       
       d0m  = d1m( i1)
       d0ub = d1ub(i1)
       d0ut = d1ut(i1)
       d0c1 = d1c1(i1)
       d0c2 = d1c2(i1)
       d0c3 = d1c3(i1)
       d0c4 = d1c4(i1)

       i0vc = 8*i1 - 3

       !C
       !C EQUATION FOR Fb
       !C
       
       d0x1 = d0jm1*d0bpb/d0m
       
       !C Fb (B,B)
       d4et(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c1*d0g4
       !C Ft (B,B)
       d4et(2,2,i0cid+2,i0nid) = -d0ftcst*d0x1*d0c1*d0g5
       !C Qb (B,B)
       d4et(2,2,i0cid+3,i0nid) = -d0qbcst*d0x1*d0c2*d0g4
       !C Qt (B,B)
       d4et(2,2,i0cid+4,i0nid) = -d0qtcst*d0x1*d0c2*d0g5

  !     i1etidv(i0cid+1) = i0vc+3
  !     i1etidv(i0cid+2) = i0vc+4
  !     i1etidv(i0cid+3) = i0vc+7
  !     i1etidv(i0cid+4) = i0vc+8
       i0cid = i0cid + 4

       !C
       !C EQUATION FOR P
       !C
       
       d0x1 = d0jm1*(d0ub*d0g4+d0ut*d0g5)
       d0x2 = d0jm1*d0bpb

       !C Fb (B,B)
       d4et(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c1*d0g4
       !C Fb (ub,B)
       d4et(2,2,i0cid+2,i0nid) =  d0fbcst*d0x2*d0c1*d0g4
       !C Ft (B,B)
       d4et(2,2,i0cid+3,i0nid) = -d0ftcst*d0x1*d0c1*d0g5
       !C Ft (ub,B)
       d4et(2,2,i0cid+4,i0nid) =  d0ftcst*d0x2*d0c1*d0g5
       !C Qb (B,B)
       d4et(2,2,i0cid+5,i0nid) = -d0qbcst*d0x1*d0c2*d0g4
       !C Qb (ub,B)
       d4et(2,2,i0cid+6,i0nid) =  d0qbcst*d0x2*d0c2*d0g4
       !C Qt (B,B)
       d4et(2,2,i0cid+7,i0nid) = -d0qtcst*d0x1*d0c2*d0g5
       !C Qt (ub,B)
       d4et(2,2,i0cid+8,i0nid) =  d0qtcst*d0x2*d0c2*d0g5
       
  !     i1etidv(i0cid+1) = i0vc+3
  !     i1etidv(i0cid+2) = i0vc+3
  !     i1etidv(i0cid+3) = i0vc+4
  !     i1etidv(i0cid+4) = i0vc+4
  !     i1etidv(i0cid+5) = i0vc+7
  !     i1etidv(i0cid+6) = i0vc+7
  !     i1etidv(i0cid+7) = i0vc+8
  !     i1etidv(i0cid+8) = i0vc+8
       i0cid = i0cid + 8

       !C
       !C EQUATION FOR Qb
       !C
       
       d0x1 = d0jm1*d0bpb
       
       !C Fb (B,B)
       d4et(2,2,i0cid+1,i0nid) = -d0fbcst*d0x1*d0c3*d0g4
       !C Ft (B,B)
       d4et(2,2,i0cid+2,i0nid) = -d0ftcst*d0x1*d0c3*d0g5
       !C Qb (B,B)
       d4et(2,2,i0cid+3,i0nid) = -d0qbcst*d0x1*d0c4*d0g4
       !C Qt (B,B)
       d4et(2,2,i0cid+4,i0nid) = -d0qtcst*d0x1*d0c4*d0g5

 !      i1etidv(i0cid+1) = i0vc+3
 !      i1etidv(i0cid+2) = i0vc+4
 !      i1etidv(i0cid+3) = i0vc+7
 !      i1etidv(i0cid+4) = i0vc+8
       i0cid = i0cid + 4
       
    ENDDO

!    OPEN(10,FILE='TEST_ETV.dat')
!    DO i1=1,i0etcmx
!       WRITE(10,*)'i1=',i1,'I1ETIDC=',i1etidv(i1)
!    ENDDO
!    CLOSE(10)
    
    If(i0cid.EQ.i0etcmx) RETURN
    
    WRITE(6,*)'ERROR IN T2CALV_ET_COEFF'
    STOP 

    RETURN
  
  END SUBROUTINE T2CALV_ET

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
    
    DO i0xb  = 1, i0spcs

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

    DO i0xb  = 1, i0spcs

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

    DO i0xb  = 1, i0spcs

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
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2) !C V/Vt_a

    d0nut = 0.D0

    DO i0xb  = 1, i0spcs

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

    DO i0xb  = 1, i0spcs
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
    
    DO i0xb  = 1, i0spcs
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
