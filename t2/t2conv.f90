MODULE T2CONV
  
  USE T2CNST

  
  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE T2_CONV
    
    USE T2COMM, ONLY:&
         i0smax,i0vmax,i0xmax,i0mmax,i2crt,&
         d0mfcst,d0btcst,d0etcst,d0epcst,d0ercst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,d0fpcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,d0qpcst,&
         d2xvec,d2xout,d0rmnr,d2jm1,d2mfc1

    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0midi,i0xid1d,i0xid2d
    
    REAL(   i0rkind)::&
         d0bp_pu,d0bt_pu,d0et_pu,d0ep_pu,d0er_pu
    
    REAL(   i0rkind),DIMENSION(1:i0smax)::&
         d1nn_pu,d1ur_pu,d1ub_pu,d1ut_pu,d1up_pu,&
         d1tt_pu,d1qr_pu,d1qb_pu,d1qt_pu,d1qp_pu
   
    REAL(   i0rkind)::&
         d0sqrtg,d0mfcr,d0sqrtr,d0sqrtri,&
         d0ctgrr,d0ctgrp,d0ctgpp,d0ctgtt,&
         d0cogrr,d0cogrp,d0cogpp,d0cogtt,&
         d0nn,d0ni,d0tt,d0qr,d0qb,d0qt,d0qp,&
         d0psip,d0cobt,d0coet,d0coep,d0coer 
   
    REAL(   i0rkind),DIMENSION(1:i0smax)::&
         d1nn,d1fr,d1fb,d1ft,d1fp,&
         d1pp,d1qr,d1qb,d1qt,d1qp
    
    d2xout(1:i0vmax,1:i0xmax) = 0.D0
    
    DO i0midi = 1, i0mmax
       
       i0xid2d = i2crt( 2,i0midi)
       i0xid1d = i2crt( 3,i0midi)
       d0mfcr  = d2mfc1(1,i0midi)
       d0sqrtr = SQRT(d0mfcr)
       d0sqrtg = d2jm1( 1,i0midi)
       d0cogrr = d2jm1( 2,i0midi)
       d0cogrp = d2jm1( 3,i0midi)
       d0cogpp = d2jm1( 4,i0midi)
       d0cogtt = d2jm1( 5,i0midi)
       d0ctgrr = d2jm1( 6,i0midi)
       d0ctgrp = d2jm1( 7,i0midi)
       d0ctgpp = d2jm1( 8,i0midi)
       d0ctgtt = d2jm1( 9,i0midi)
       
       IF(d0sqrtr.GT.0.D0)THEN
          d0sqrtri = 1.D0/d0sqrtr
       ELSE
          d0sqrtri = 0.D0
       ENDIF
       
       !C
       !C INITTIALIZATION
       !C
       
       d0psip = d2xvec(1,i0xid1d)*d0mfcst
       d0cobt = d2xvec(2,i0xid1d)*d0btcst
       d0coet = d2xvec(3,i0xid1d)*d0etcst
       d0coep = d2xvec(4,i0xid2d)*d0epcst*d0sqrtr
       d0coer = d2xvec(5,i0xid2d)*d0ercst
       
       DO i0sidi = 1, i0smax
          i0vidi = 10*i0sidi - 5
          d1nn(i0sidi) = d2xvec(i0vidi+ 1,i0xid2d)*d0nncst
          d1fr(i0sidi) = d2xvec(i0vidi+ 2,i0xid2d)*d0frcst
          d1fb(i0sidi) = d2xvec(i0vidi+ 3,i0xid2d)*d0fbcst
          d1ft(i0sidi) = d2xvec(i0vidi+ 4,i0xid2d)*d0ftcst
          d1fp(i0sidi) = d2xvec(i0vidi+ 5,i0xid2d)*d0fpcst
          d1pp(i0sidi) = d2xvec(i0vidi+ 6,i0xid2d)*d0ppcst
          d1qr(i0sidi) = d2xvec(i0vidi+ 7,i0xid2d)*d0qrcst
          d1qb(i0sidi) = d2xvec(i0vidi+ 8,i0xid2d)*d0qbcst
          d1qt(i0sidi) = d2xvec(i0vidi+ 9,i0xid2d)*d0qtcst
          d1qp(i0sidi) = d2xvec(i0vidi+10,i0xid2d)*d0qpcst
       ENDDO
       
       !C
       !C CONVERSION TO PHYSICAL UNIT 
       !C
       
       !C
       !C d0bp_pu: Poroidal Magnetic Field [T]
       !C
       
       !d0bp_pu = d0psip*SQRT(d0cogpp)/d0sqrtg
       d0bp_pu = d0psip*SQRT(d0cogpp)/d0sqrtg
       
       !C
       !C d0bt_pu: Toroidal Magnetic Field [T]
       !C
       
       !d0bt_pu = d0cobt*SQRT(d0ctgtt)
       d0bt_pu = d0cobt*SQRT(d0ctgtt)
       
       !C
       !C d0et_pu: Toroidal Electric Field [V/m]
       !C
       
       !d0et_pu = d0coet*SQRT(d0ctgtt)
       d0et_pu = d0coet*SQRT(d0ctgtt)
       
       !C
       !C d0ep_pu: Poroidal Electric Field [mV/m]
       !C
       
       d0ep_pu = d0coep*SQRT(d0ctgpp)*1.D3
       
       !C
       !C d0er_pu: Radial   Electric Field [kV/m]
       !C
       
       d0er_pu = d0coer*SQRT(d0ctgrr)*1.D-3
       
       DO i0sidi = 1, i0smax
          
          d0nn = d1nn(i0sidi)
          
          IF(d0nn.GT.0.D0)THEN
             d0ni = 1.D0/d0nn
          ELSE
             d0ni = 0.D0
          ENDIF
          
          !C
          !C d1nn_pu: Particle Density [10^20 m^-3]
          !C
          
          d1nn_pu(i0sidi) = d0nn*1.D-20
          
          !C
          !C d1ur_pu: Radial velocity  [m/s]
          !C
          
          d1ur_pu(i0sidi) = d1fr(i0sidi)*d0ni*SQRT(d0cogrr)*d0sqrtri
          
          !C
          !C d1ub_pu: Parallel velocity  [km/s]
          !C
          
          d1ub_pu(i0sidi) = d1fb(i0sidi)*d0ni*1.D-3
          
          !C
          !C d1ut_pu: Toroidal velocity  [km/s]
          !C
          
          d1ut_pu(i0sidi) = d1ft(i0sidi)*d0ni*SQRT(d0ctgtt)*1.D-3
          
          !C
          !C d1up_pu: Poloidal velocity  [km/s]
          !C
          
          d1up_pu(i0sidi) = d1fp(i0sidi)*d0ni*SQRT(d0cogpp)*1.D-3
          
          !C
          !C d1tt_pu: Temperature [keV]
          !C
          
          d0tt =  d1pp(i0sidi)*d0ni ! Joule
          d1tt_pu(i0sidi) = d0tt/d0aee*1.D-3 ! keV
          
          !C 
          !C d1qr_pu: Radial heat Flux [kJ*m/s]
          !C
          
          d0qr = d1qr(i0sidi) - 2.5D0*d0tt*d1fr(i0sidi)
          d1qr_pu(i0sidi) = d0qr*SQRT(d0cogrr)*d0sqrtri*1.D-3 
          
          
          !C 
          !C d1qb_pu: Parallel heat Flux [MJ*m/s]
          !C
          
          d0qb = d1qb(i0sidi) - 2.5D0*d0tt*d1fb(i0sidi)
          d1qb_pu(i0sidi) = d0qb*1.D-6
          
          !C 
          !C d1qb_pu: Toroidal heat Flux [MJ*m/s]
          !C
          
          d0qt = d1qt(i0sidi) - 2.5D0*d0tt*d1ft(i0sidi)
          d1qt_pu(i0sidi) = d0qt*SQRT(d0ctgtt)*1.D-6
          
          !C 
          !C d1qp_pu: Toroidal heat Flux [MJ*m/s]
          !C
          
          d0qp = d1qp(i0sidi) - 2.5D0*d0tt*d1fp(i0sidi)
          d1qp_pu(i0sidi) = d0qp*SQRT(d0cogpp)*1.D-6
          
       ENDDO
       
       d2xout(1,i0xid2d) = d0bp_pu
       d2xout(2,i0xid2d) = d0bt_pu
       d2xout(3,i0xid2d) = d0et_pu
       d2xout(4,i0xid2d) = d0ep_pu
       d2xout(5,i0xid2d) = d0er_pu
       
       DO i0sidi = 1, i0smax
          i0vidi = 10*i0sidi - 5
          d2xout(i0vidi+ 1,i0xid2d) = d1nn_pu(i0sidi)
          d2xout(i0vidi+ 2,i0xid2d) = d1ur_pu(i0sidi)
          d2xout(i0vidi+ 3,i0xid2d) = d1ub_pu(i0sidi)
          d2xout(i0vidi+ 4,i0xid2d) = d1ut_pu(i0sidi)
          d2xout(i0vidi+ 5,i0xid2d) = d1up_pu(i0sidi)
          d2xout(i0vidi+ 6,i0xid2d) = d1tt_pu(i0sidi)
          d2xout(i0vidi+ 7,i0xid2d) = d1qr_pu(i0sidi)
          d2xout(i0vidi+ 8,i0xid2d) = d1qb_pu(i0sidi)
          d2xout(i0vidi+ 9,i0xid2d) = d1qt_pu(i0sidi)
          d2xout(i0vidi+10,i0xid2d) = d1qp_pu(i0sidi)
       ENDDO
       
    END DO
    
    RETURN
    
  END SUBROUTINE T2_CONV
  
END MODULE T2CONV
