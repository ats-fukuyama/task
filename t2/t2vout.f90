MODULE T2VOUT
  
  USE T2CNST
  
  
  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE T2VOUT_EXECUTE
    
    USE T2COMM, ONLY:&
         NSMAX,NVMAX,NXMAX,NMMAX,NFMAX,i2crt,&
         BpNF,BtNF,EtNF,EpNF,ErNF,&
         NnNF,FrNF,FbNF,FtNF,FpNF,&
         PpNF,QrNF,QbNF,QtNF,QpNF,&
         Xvec,d2xout,Metric,GlobalCrd,EqSet

    INTEGER(ikind)::&
         i0sidi,i0vidi,i0midi,i0xid1d,i0xid2d
    
    REAL(   rkind)::&
         d0bp_pu,d0bt_pu,d0et_pu,d0ep_pu,d0er_pu
    
    REAL(   rkind),DIMENSION(1:NSMAX)::&
         d1nn_pu,d1ur_pu,d1ub_pu,d1ut_pu,d1up_pu,&
         d1tt_pu,d1qr_pu,d1qb_pu,d1qt_pu,d1qp_pu
   
    REAL(   rkind)::&
         d0sqrtg,d0mfcr,d0sqrtr,d0sqrtri,&
         d0ctgrr,d0ctgrp,d0ctgpp,d0ctgtt,&
         d0cogrr,d0cogrp,d0cogpp,d0cogtt,&
         d0nn,d0ni,d0tt,d0qr,d0qb,d0qt,d0qp,&
         d0psip,d0cobt,d0coet,d0coep,d0coer 
   
    REAL(   rkind),DIMENSION(1:NSMAX)::&
         d1nn,d1fr,d1fb,d1ft,d1fp,&
         d1pp,d1qr,d1qb,d1qt,d1qp
    
    d2xout(1:NVMAX,1:NXMAX) = 0.D0
    
    DO i0midi = 1, NMMAX
       
       i0xid2d = i2crt( 2,i0midi)
       i0xid1d = i2crt( 3,i0midi)
       d0mfcr  = GlobalCrd(1,i0midi)
       d0sqrtr = SQRT(d0mfcr)
       d0sqrtg = Metric( 1,i0midi)
       d0cogrr = Metric( 2,i0midi)
       d0cogrp = Metric( 3,i0midi)
       d0cogpp = Metric( 4,i0midi)
       d0cogtt = Metric( 5,i0midi)
       d0ctgrr = Metric( 6,i0midi)
       d0ctgrp = Metric( 7,i0midi)
       d0ctgpp = Metric( 8,i0midi)
       d0ctgtt = Metric( 9,i0midi)
       
       IF(d0sqrtr.GT.0.D0)THEN
          d0sqrtri = 1.D0/d0sqrtr
       ELSE
          d0sqrtri = 0.D0
       ENDIF
       
       !C
       !C INITTIALIZATION
       !C
       SELECT CASE (EqSet)
       
       CASE (1)
          d0psip = Xvec(1,i0xid1d)*BpNF
          d0cobt = Xvec(2,i0xid1d)*BtNF
          d0coet = Xvec(3,i0xid1d)*EtNF
          d0coep = Xvec(4,i0xid2d)*EpNF*d0mfcr
          d0coer = Xvec(5,i0xid2d)*ErNF
          
          DO i0sidi = 1, NSMAX
             i0vidi = 10*(i0sidi-1) + NFMAX
             d1nn(i0sidi) = Xvec(i0vidi+ 1,i0xid2d)*NnNF
             d1fr(i0sidi) = Xvec(i0vidi+ 2,i0xid2d)*FrNF*d0mfcr
             d1fb(i0sidi) = Xvec(i0vidi+ 3,i0xid2d)*FbNF
             d1ft(i0sidi) = Xvec(i0vidi+ 4,i0xid2d)*FtNF
             d1fp(i0sidi) = Xvec(i0vidi+ 5,i0xid2d)*FpNF
             d1pp(i0sidi) = Xvec(i0vidi+ 6,i0xid2d)*PpNF
             d1qr(i0sidi) = Xvec(i0vidi+ 7,i0xid2d)*QrNF*d0mfcr
             d1qb(i0sidi) = Xvec(i0vidi+ 8,i0xid2d)*QbNF
             d1qt(i0sidi) = Xvec(i0vidi+ 9,i0xid2d)*QtNF
             d1qp(i0sidi) = Xvec(i0vidi+10,i0xid2d)*QpNF
          ENDDO
       CASE (2)
          
          d0psip = Xvec(1,i0xid1d)*BpNF
          d0cobt = Xvec(2,i0xid1d)*BtNF
          d0coet = Xvec(3,i0xid1d)*EtNF
          d0coep = Xvec(4,i0xid2d)*EpNF
          d0coer = Xvec(5,i0xid2d)*ErNF
          
          DO i0sidi = 1, NSMAX
             i0vidi = 10*(i0sidi-1) + NFMAX
             d1nn(i0sidi) = Xvec(i0vidi+ 1,i0xid2d)*NnNF
             d1fr(i0sidi) = Xvec(i0vidi+ 2,i0xid2d)*FrNF
             d1fb(i0sidi) = Xvec(i0vidi+ 3,i0xid2d)*FbNF
             d1ft(i0sidi) = Xvec(i0vidi+ 4,i0xid2d)*FtNF
             d1fp(i0sidi) = Xvec(i0vidi+ 5,i0xid2d)*FpNF
             d1pp(i0sidi) = Xvec(i0vidi+ 6,i0xid2d)*PpNF
             d1qr(i0sidi) = Xvec(i0vidi+ 7,i0xid2d)*QrNF
             d1qb(i0sidi) = Xvec(i0vidi+ 8,i0xid2d)*QbNF
             d1qt(i0sidi) = Xvec(i0vidi+ 9,i0xid2d)*QtNF
             d1qp(i0sidi) = Xvec(i0vidi+10,i0xid2d)*QpNF
          ENDDO
          
       END SELECT
       !
       ! CONVERSION TO PHYSICAL UNIT 
       !
       
       ! d0bp_pu: Poroidal Magnetic Field [T]      
       d0bp_pu = d0psip*SQRT(d0cogpp)/d0sqrtg
       
       ! d0bt_pu: Toroidal Magnetic Field [T]
       d0bt_pu = d0cobt*SQRT(d0ctgtt)
       
       ! d0et_pu: Toroidal Electric Field [V/m]
       d0et_pu = d0coet*SQRT(d0ctgtt)
       
       ! d0ep_pu: Poroidal Electric Field [mV/m]
       IF(d0mfcr.EQ.0.D0)THEN
          d0ep_pu = 0.D0
       ELSE
          d0ep_pu = d0coep*SQRT(d0ctgpp/d0mfcr)*1.D3
       ENDIF
       
       ! d0er_pu: Radial   Electric Field [kV/m]
       d0er_pu = d0coer*SQRT(d0ctgrr)*1.D-3
       
       DO i0sidi = 1, NSMAX
          
          d0nn = d1nn(i0sidi)
          
          IF(d0nn.GT.0.D0)THEN
             d0ni = 1.D0/d0nn
          ELSE
             d0ni = 0.D0
          ENDIF
          
          ! d1nn_pu: Particle Density [10^20 m^-3]         
          d1nn_pu(i0sidi) = d0nn*1.D-20
          
          ! d1ur_pu: Radial velocity  [m/s]
          IF(d0mfcr.EQ.0)THEN
             d1ur_pu(i0sidi) = 0.D0
          ELSE
             d1ur_pu(i0sidi) = d1fr(i0sidi)*d0ni*SQRT(d0cogrr/d0mfcr)
          END IF
          
          ! d1ub_pu: Parallel velocity  [km/s]          
          d1ub_pu(i0sidi) = d1fb(i0sidi)*d0ni*1.D-3
          
          ! d1ut_pu: Toroidal velocity  [km/s]
          d1ut_pu(i0sidi) = d1ft(i0sidi)*d0ni*SQRT(d0ctgtt)*1.D-3

          ! d1up_pu: Poloidal velocity  [km/s]
          d1up_pu(i0sidi) = d1fp(i0sidi)*d0ni*SQRT(d0cogpp)*1.D-3
          
          ! d1tt_pu: Temperature [keV]
          d0tt =  d1pp(i0sidi)*d0ni ! Joule
          d1tt_pu(i0sidi) = d0tt/d0aee*1.D-3 ! keV
          
          ! d1qr_pu: Radial heat Flux [kJ*m/s]
          d0qr = d1qr(i0sidi) - 2.5D0*d0tt*d1fr(i0sidi)
          IF(d0mfcr.EQ.0.D0)THEN
             d1qr_pu(i0sidi) = 0.D0
          ELSE
             d1qr_pu(i0sidi) = d0qr*SQRT(d0cogrr/d0mfcr)*1.D-3 
          END IF
          
          ! d1qb_pu: Parallel heat Flux [MJ*m/s]
          d0qb = d1qb(i0sidi) - 2.5D0*d0tt*d1fb(i0sidi)
          d1qb_pu(i0sidi) = d0qb*1.D-6
          
          ! d1qt_pu: Toroidal heat Flux [MJ*m/s]
          d0qt = d1qt(i0sidi) - 2.5D0*d0tt*d1ft(i0sidi)
          d1qt_pu(i0sidi) = d0qt*SQRT(d0ctgtt)*1.D-6
          
          !C d1qp_pu: Toroidal heat Flux [MJ*m/s]    
          d0qp = d1qp(i0sidi) - 2.5D0*d0tt*d1fp(i0sidi)
          d1qp_pu(i0sidi) = d0qp*SQRT(d0cogpp)*1.D-6
          
       ENDDO
       
       d2xout(1,i0xid2d) = d0bp_pu
       d2xout(2,i0xid2d) = d0bt_pu
       d2xout(3,i0xid2d) = d0et_pu
       d2xout(4,i0xid2d) = d0ep_pu
       d2xout(5,i0xid2d) = d0er_pu
       
       DO i0sidi = 1, NSMAX
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
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VOUT_EXECUTE
END MODULE T2VOUT
