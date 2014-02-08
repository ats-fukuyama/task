!C--------------------------------------------------------------------
!C 
!C  T2MFCS
!C
!C
!C
!C
!C
!C
!C--------------------------------------------------------------------
MODULE T2PROF
  
  USE T2CNST,ONLY:&
       i0ikind,i0rkind

  IMPLICIT NONE

CONTAINS

  SUBROUTINE T2_PROF


    USE T2COMM,ONLY:&
         i0mfcs
    SELECT CASE(i0mfcs)
       ! SELECT COORDINATE TYPE
       ! 1:     TOROIDAL COORDINATE W/O EQUILIBRIUM (\rho = r  )
       
    CASE(1)
       CALL T2PROF_TOROIDAL
    CASE DEFAULT
       WRITE(6,*)'IMPROPER INPUT >> I0MFCS'
       STOP
    END SELECT


    RETURN
    
  END SUBROUTINE T2_PROF
  
  SUBROUTINE T2PROF_TOROIDAL
    
    USE T2COMM,ONLY:&
         d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
         d0nncst,d0frcst,d0fbcst,d0ftcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         i0smax,i0xmax,i0vmax,i0mmax,i1pdn1,i2crt,&
         d0rmnr,d0rmjr,d2mfc1,d2rzm, d2rzx, d2jm1,d2xvec

  USE T2GOUT, ONLY: T2_GOUT
  USE T2CONV, ONLY: T2_CONV
    INTEGER(i0ikind)::&
         i0vidi,i0sidi,i0midi,i0xidi,i0xid1d,i0xid2d
    
    REAL(   i0rkind)::d0mfcr,d0mfcp,d0jm1,d0cobt
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1n0,d1p0
    REAL(   i0rkind),DIMENSION(1:6,1:i0smax)::d2f0
    
100 FORMAT( 6E15.8)
110 FORMAT(10E15.8)
120 FORMAT( 5E15.8)
    
    d2xvec(1:i0vmax,1:i0xmax) = 0.D0
    
    DO i0midi = 1, i0mmax
       
       !C
       !C INITIIALIZATION
       !C
       
       d0mfcr = 0.D0
       d0mfcp = 0.D0
       
       !C
       !C CYLINDRICAL COORDINATES
       !C
       
       d0mfcr = d2mfc1(1,i0midi)
       d0mfcp = d2mfc1(2,i0midi)
       
       d2rzm(1,i0midi)  = fd0rzcr(d0mfcr,d0mfcp)
       d2rzm(2,i0midi)  = fd0rzcz(d0mfcr,d0mfcp)
       
       !C
       !C CONTRAVARIANT GEOMETRIC TENSOR
       !C
       
       d2jm1(1:5,i0midi)= fd1mc(d0mfcr,d0mfcp)
       
       i0xid2d = i2crt(2,i0midi)
       i0xid1d = i2crt(3,i0midi)
       
       !C PSI'
       !d2xvec(1,i0xid1d) = fd0psip(d0mfcr,d0mfcp)/d0mfcst 
       d2xvec(1,i0xid2d) = fd0psip(d0mfcr,d0mfcp)/d0mfcst 
       !C I
       !d2xvec(2,i0xid1d) = fd0cobt(d0mfcr,d0mfcp)/d0btcst 
       d2xvec(2,i0xid2d) = fd0cobt(d0mfcr,d0mfcp)/d0btcst 

       !C E_{\zeta}
       !d2xvec(3,i0xid1d) = fd0coet(d0mfcr,d0mfcp)/d0etcst 
       d2xvec(3,i0xid2d) = fd0coet(d0mfcr,d0mfcp)/d0etcst 
       !C E_{\chi }
       d2xvec(4,i0xid2d) = fd0coep(d0mfcr,d0mfcp)/d0epcst
       !C E_{\rho }
       d2xvec(5,i0xid2d) = fd0coer(d0mfcr,d0mfcp)/d0ercst
       
       !C INITIAL PROFILE: Fr Fb Fb Qr Qb Qt (DIMENSIONLESS)
       
       d1n0 = fd1n0(d0mfcr,d0mfcp)
       d1p0 = fd1p0(d0mfcr,d0mfcp)
       d2f0 = fd2f0(d0mfcr,d0mfcp)
       
       DO i0sidi = 1,i0smax
          i0vidi = 8*i0sidi - 3
          d2xvec(i0vidi+1,i0xid2d) = d1n0(  i0sidi)/d0nncst
          d2xvec(i0vidi+2,i0xid2d) = d2f0(1,i0sidi)/d0frcst
          d2xvec(i0vidi+3,i0xid2d) = d2f0(2,i0sidi)/d0fbcst
          d2xvec(i0vidi+4,i0xid2d) = d2f0(3,i0sidi)/d0ftcst
          d2xvec(i0vidi+5,i0xid2d) = d1p0(  i0sidi)/d0ppcst
          d2xvec(i0vidi+6,i0xid2d) = d2f0(4,i0sidi)/d0qrcst
          d2xvec(i0vidi+7,i0xid2d) = d2f0(5,i0sidi)/d0qbcst
          d2xvec(i0vidi+8,i0xid2d) = d2f0(6,i0sidi)/d0qtcst
       ENDDO
    ENDDO
    
    DO i0midi = 1, i0mmax
       i0xidi = i2crt(2,i0midi)
       d2rzx(1,i0xidi) = d2rzm(1,i0midi)
       d2rzx(2,i0xidi) = d2rzm(2,i0midi)
    ENDDO
!    CALL T2_CONV
!    CALL T2_GOUT
    
    RETURN
    
  END SUBROUTINE T2PROF_TOROIDAL

  !C------------------------------------------------------------------
  !C
  !C RADIAL PROFILE CALCURATION 
  !C 
  !C for 0 < r < 1
  !C   f(r) = (fc-fs)*(1-r^n)^m + fs 
  !C
  !C for 1 < r < rw
  !C   f(r) = fw + a1*(r-rw)^2 + a2*(r-rw)^3 + a1*(r-rw)^4 
  !C
  !C
  !C------------------------------------------------------------------
  FUNCTION fd1rf(i0m0,i0n0,d0fc,d0fs,d0fw,d0rw,d0rr)
    
    INTEGER(i0ikind),INTENT(IN )::i0m0,i0n0
    REAL(   i0rkind),INTENT(IN )::d0fc,d0fs,d0fw,d0rw,d0rr
    REAL(   i0rkind),DIMENSION(1:3)::fd1rf
    REAL(   i0rkind)::&
         d0a1,d0a2,d0a3,d0b1,d0b2,d0b3,&
         d0m0,d0n0,d0rx,&
         d0fcs,d0fsw,d0rsw
    INTEGER(i0ikind)::&
         i0m1,i0n1,i0n2
    
    IF(i0n0.LT.2)THEN 
       WRITE(6,*)'INAPPROPRIATE PROFILE PARAMETER: n',i0n0
       STOP
    ENDIF

    d0n0  = DBLE(i0n0)
    d0m0  = DBLE(i0m0)
       
    d0fcs = d0fc - d0fs
    i0n1  = i0n0 - 1
    i0n2  = i0n0 - 2
    i0m1  = i0m0 - 1
    
    !C
    !C  0 <= r <= 1
    !C
    
    IF((d0rr.GE.0.D0).AND.(d0rr.LE.1.D0))THEN
       
       d0rx  = 1.D0-(d0rr**i0n0)
       
       SELECT CASE(i0m0)
       CASE(1)
          fd1rf(1)    =  d0fcs*d0rx + d0fs
          fd1rf(2)    = -d0fcs      *(d0rr**i0n1)*d0n0
          IF(    i0n2.NE.0)THEN
             fd1rf(3) = -d0fcs      *(d0rr**i0n2)*d0n0
          ELSEIF(i0n2.EQ.0)THEN
             fd1rf(3) = -d0fcs                   *d0n0
          ENDIF
       CASE(2:)
          fd1rf(1)    =  d0fcs*(d0rx**i0m0) + d0fs
          fd1rf(2)    = -d0fcs*(d0rx**i0m1)*(d0rr**i0n1)*d0n0*d0m0
          IF(    i0n2.NE.0)THEN
             fd1rf(3) = -d0fcs*(d0rx**i0m1)*(d0rr**i0n2)*d0n0*d0m0
          ELSEIF(i0n2.EQ.0)THEN
             fd1rf(3) = -d0fcs*(d0rx**i0m1)             *d0n0*d0m0
          ENDIF
       END SELECT
    ELSEIF(d0rr.GT.1.D0)THEN
       
       d0rsw = 1.D0/(1.D0 - d0rw)
       d0fsw = d0fs - d0fw
       d0rx = d0rr - d0rw

       SELECT CASE(i0m0)
       CASE (1)
          d0b1 =  d0fsw*(d0rsw**2)
          d0b2 = -d0fcs* d0rsw    *d0n0
          d0b3 = -d0fcs           *d0n0*(d0n0-1.D0)
          
       CASE (2)
          d0b1 =  d0fsw*(d0rsw**2)
          d0b2 =  0.D0
          d0b3 =  d0fcs           *d0n0*d0n0*2.D0
       CASE (3:)
          d0b1 = d0fsw*(d0rsw**2)
          d0b2 = 0.D0
          d0b3 = 0.D0
       CASE DEFAULT
          WRITE(6,*)'INAPPROPRIATE PROFILE PARAMETER: m',i0m0
          STOP
       END SELECT
       
       d0a1 =  6.0D0*d0b1 - 3.0D0*d0b2 + 0.5D0*d0b3
       d0a2 = -8.0D0*d0b1 + 5.0D0*d0b2 - 1.0D0*d0b3
       d0a3 =  3.0D0*d0b1 - 2.0D0*d0b2 + 0.5D0*d0b3
       
       d0a2 = d0a2*d0rsw
       d0a3 = d0a3*d0rsw*d0rsw
       
       fd1rf(1) =      d0a1*(d0rx**2) +      d0a2*(d0rx**3)&
            &   +      d0a3*(d0rx**4) + d0fw  
       fd1rf(2) = 2.D0*d0a1* d0rx     + 3.D0*d0a2*(d0rx**2) &
            &   + 4.D0*d0a3*(d0rx**3) 
       fd1rf(3) = fd1rf(2)/d0rr
    
    ENDIF
    
    RETURN
    
  END FUNCTION fd1rf

  FUNCTION fd0rzcr(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY: d0rmjr,d0rmnr
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0rzcr
    
    fd0rzcr = d0rmjr + d0rmnr*d0mfcr * COS(d0mfcp)
    
    RETURN
    
  END FUNCTION fd0rzcr

  FUNCTION fd0rzcz(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY: d0rmjr,d0rmnr
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0rzcz
    
    fd0rzcz =        - d0rmnr*d0mfcr * SIN(d0mfcp)
    
    RETURN

  END FUNCTION fd0rzcz
  
  FUNCTION fd1mc(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY: d0rmjr,d0rmnr
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind),DIMENSION(1:5)::fd1mc
    REAL(i0rkind)::d0rzcr

    !C CYLINDRICAL COORDINATES AND MSCS GEOMETRIC TENSOR
    !C
    !C    MFC: MAGNETIC SURFACE COORDINATE SYSTEM [rho,chi,zeta]
    !C    RZC: CYLINDRICAL COORDINATE SYSTEM      [  R,phi,Z]
    !C
    !C    R   = R0+a*rho*COS(chi)
    !C    Z   =   -a*rho*SIN(chi)
    !C    phi = zeta
    !C
    !C    FD1MC(1) : SQRT{g} 
    !C    FD1MC(2) : g^{rho  rho } 
    !C    FD1MC(3) : g^{rho  chi } 
    !C    FD1MC(4) : g^{chi  chi } 
    !C    FD1MC(5) : g^{zeta zeta} 
    !C
    
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    
    IF(d0mfcr.GT.0.D0)THEN
       fd1mc(1)= (d0rmnr**2)*d0mfcr*d0rzcr
       fd1mc(2)= 1.d0/d0rmnr**2
       fd1mc(3)= 0.d0 
       fd1mc(4)= 1.d0/((d0rmnr**2)*(d0mfcr**2))
       fd1mc(5)= 1.d0/(d0rzcr**2)
    ELSE
       fd1mc(1)= (d0rmnr**2)*d0mfcr*d0rzcr
       fd1mc(2)= 1.d0/d0rmnr**2
       fd1mc(3)= 0.d0 
       fd1mc(4)= 0.D0
       fd1mc(5)= 1.d0/(d0rzcr**2)
    ENDIF
    
    RETURN
    
  END FUNCTION fd1mc
  
  FUNCTION fd0q0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0rmnr,d0qc,d0qs
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0q0
    
    IF((d0mfcr.GE.0.D0).AND.(d0mfcr.LE.1.D0))THEN
       fd0q0 = (d0qc-d0qs)*(1.D0 - d0mfcr**2)+d0qs
    ELSEIF(d0mfcr.GT.1.D0)THEN
       fd0q0 = (d0qs-d0qc)*(       d0mfcr**2)+d0qc
    ELSE
       WRITE(6,*)'WRONG RHO INPUT'
       STOP
    ENDIF
    
    RETURN
    
  END FUNCTION fd0q0

  FUNCTION fd0q1(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0rmnr,d0qc,d0qs
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0q1
    
    fd0q1 = 2.D0*(d0qs-d0qc)*d0mfcr
    
    RETURN
    
  END FUNCTION fd0q1

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF d\psi/d\rho
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0psip(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmnr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0psip
    REAL(i0rkind)::d0q0
    
    d0q0    = fd0q0(d0mfcr,d0mfcp)
    fd0psip = d0bc*(d0rmnr**2)*d0mfcr/d0q0
    
    RETURN
    
  END FUNCTION fd0psip
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF I
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0cobt(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0cobt
    
    fd0cobt = d0rmjr*d0bc

    RETURN
    
  END FUNCTION fd0cobt

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF B
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0bb(d0mfcr,d0mfcp)
        
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0bb
    REAL(i0rkind),DIMENSION(1:5)::d1mc
    REAL(i0rkind)::d0psip,d0cobt,d0bb
    
    d0psip = fd0psip(d0mfcr,d0mfcp)
    d0cobt = fd0cobt(d0mfcr,d0mfcp)
    d1mc   = fd1mc(  d0mfcr,d0mfcp)
    
    d0bb = (d0psip**2)*d1mc(2)+(d0cobt**2)
    d0bb = d0bb*d1mc(5)
    fd0bb = SQRT(d0bb)
    
    RETURN

  END FUNCTION fd0bb

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF E_{\zeta}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0coet(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d0iar
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0coet
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1t0
    REAL(   i0rkind)::d0t0,d0clog,d0eta,d0jt
    
    d1t0   = fd1t0(d0mfcr,d0mfcp)
    d0t0   = d1t0(1)/(1.D3*d0aee) !C  Electron temperature in keV
    d0clog = 17.D0 ! Coulomb logarithm for debug 
    d0jt   = fd0jt(d0mfcr,d0mfcp)  
    
    d0eta = (1.65D-9)*d0clog/SQRT(d0t0)**3
    d0eta = d0eta/((1.D0-SQRT(d0iar*d0mfcr))**2)
    
    fd0coet = d0eta*d0jt
    
    RETURN
    
  END FUNCTION fd0coet

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF E_{\chi}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0coep(d0mfcr,d0mfcp)
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0coep
    
    fd0coep = 0.D0
    
    RETURN
    
  END FUNCTION fd0coep
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF E_{\rho}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0coer(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0coer
    INTEGER(i0ikind)::i0sidi
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1p1,d1n0
    REAL(   i0rkind)::d0n0,d0p1
    
    fd0coer = 0.D0
    d1p1 = fd1p1(d0mfcr,d0mfcp)
    d1n0 = fd1n0(d0mfcr,d0mfcp)
    
    DO i0sidi = 2, i0smax
       d0n0 = d1n0(1) 
       d0p1 = d1p1(i0sidi) 
       fd0coer = fd0coer + d0p1/(d0aee*d0n0)
    ENDDO
    
    RETURN
    
  END FUNCTION fd0coer
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF n_{a}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1n0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1n0
    REAL(   i0rkind),DIMENSION(1:3)::d1rn
    REAL(   i0rkind)::d0nc,d0ns,d0nw
    INTEGER(i0ikind)::i0sidi
        
    DO i0sidi = 1,i0smax
       
       d0nc = d1nc(i0sidi)
       d0ns = d1ns(i0sidi)
       d0nw = d1nw(i0sidi)
       
       !C
       !C density in 10^{20} m^{-3}
       !C
       
       d1rn = fd1rf(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr)
       
       !C
       !C density in m^{-3}
       !C
       
       fd1n0(i0sidi) = d1rn(1)*1.D20
       
    ENDDO
    
    RETURN

  END FUNCTION fd1n0

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF dn_{a}/d\rho
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1n1(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1n1
    REAL(   i0rkind),DIMENSION(1:3)::d1rn
    REAL(   i0rkind)::d0nc,d0ns,d0nw
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1, i0smax
       
       d0nc = d1nc(i0sidi)
       d0ns = d1ns(i0sidi)
       d0nw = d1nw(i0sidi)
       
       !C
       !C d n_{a} /d \rho in 10^{20} m^{-3}
       !C
       
       d1rn = fd1rf(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr)
       
       !C
       !C d n_{a} /d \rho in  m^{-3}
       !C
       
       fd1n1(i0sidi) = d1rn(2)*1.D20
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1n1

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF (1/\rho)(dn_{a}/d\rho)
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1n2(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1n2
    REAL(   i0rkind),DIMENSION(1:3)::d1rn
    REAL(   i0rkind)::d0nc,d0ns,d0nw
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1, i0smax
       
       d0nc = d1nc(i0sidi)
       d0ns = d1ns(i0sidi)
       d0nw = d1nw(i0sidi)
       
       !C
       !C (1/\rho)(d n_{a} /d \rho) in 10^{20} m^{-3}
       !C
    
       d1rn = fd1rf(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr)
       
       !C
       !C (1/\rho)(d n_{a} /d \rho) in m^{-3}
       !C
       
       fd1n2(i0sidi) = d1rn(3)*1.D20
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1n2
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF T_{a}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1t0(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1t0
    REAL(   i0rkind),DIMENSION(1:3)::d1rt
    REAL(   i0rkind)::d0tc,d0ts,d0tw
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1, i0smax
       
       d0tc = d1tc(i0sidi)
       d0ts = d1ts(i0sidi)
       d0tw = d1tw(i0sidi)
       
       !C
       !C T_{a} in keV
       !C
       
       d1rt = fd1rf(2,2,d0tc,d0ts,d0tw,d0rw,d0mfcr)
       
       !C
       !C T_{a} in J
       !C
       
       fd1t0(i0sidi) = d1rt(1)*1.D3*d0aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1t0

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF dT_{a}/d\rho
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------  
  FUNCTION fd1t1(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1t1
    REAL(   i0rkind),DIMENSION(1:3)::d1rt
    REAL(   i0rkind)::d0tc,d0ts,d0tw
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1,i0smax
       
       d0tc = d1tc(i0sidi)
       d0ts = d1ts(i0sidi)
       d0tw = d1tw(i0sidi)

       !C
       !C dT_{a}/d\rho in keV
       !C
       
       d1rt = fd1rf(2,2,d0tc,d0ts,d0tw,d0rw,d0mfcr)

       !C
       !C dT_{a}/d\rho in J
       !C
       
       fd1t1(i0sidi) = d1rt(2)*1.D3*d0aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1t1

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF (1/\rho)(dT_{a}/d\rho)
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------  
  FUNCTION fd1t2(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1t2
    REAL(   i0rkind),DIMENSION(1:3)::d1rt
    REAL(   i0rkind)::d0tc,d0ts,d0tw
    INTEGER(i0ikind)::i0sidi
        
    DO i0sidi = 1,i0smax
       
       d0tc = d1tc(i0sidi)
       d0ts = d1ts(i0sidi)
       d0tw = d1tw(i0sidi)

       !C
       !C (1/\rho) dT_{a}/d\rho in keV
       !C

       d1rt = fd1rf(2,2,d0tc,d0ts,d0tw,d0rw,d0mfcr)

       !C
       !C dT_{a}/d\rho in J
       !C
       
       fd1t2(i0sidi) = d1rt(3)*1.D3*d0aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1t2

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF p_{a}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1p0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1p0,d1n0,d1t0
    INTEGER(i0ikind)::i0sidi
    
    d1n0( 1:i0smax) = fd1n0(d0mfcr,d0mfcp)
    d1t0( 1:i0smax) = fd1t0(d0mfcr,d0mfcp)
    
    DO i0sidi = 1, i0smax

       !C
       !C p_{a} in J/m^{3}
       !C

       fd1p0(i0sidi) = d1n0(i0sidi)*d1t0(i0sidi)

    ENDDO

    RETURN
    
  END FUNCTION fd1p0

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF dp_{a}/d\rho
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1p1(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1p1
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1n0,d1n1,d1t0,d1t1
    INTEGER(i0ikind)::i0sidi
    
    d1n0( 1:i0smax) = fd1n0(d0mfcr,d0mfcp)    
    d1n1( 1:i0smax) = fd1n1(d0mfcr,d0mfcp)

    d1t0( 1:i0smax) = fd1t0(d0mfcr,d0mfcp)
    d1t1( 1:i0smax) = fd1t1(d0mfcr,d0mfcp)
    

    DO i0sidi = 1, i0smax

       !C
       !C dp_{a}/d\rho in J/m^{3}
       !C

       fd1p1(i0sidi) &
            = d1n0(i0sidi)*d1t1(i0sidi)&
            + d1n1(i0sidi)*d1t0(i0sidi)
    ENDDO
    
    RETURN
    
  END FUNCTION fd1p1

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF (1/\rho)(dp_{a}/d\rho)
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1p2(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1p2
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1n0,d1n2,d1t0,d1t2
    INTEGER(i0ikind)::i0sidi
    
    d1n0( 1:i0smax) = fd1n0(d0mfcr,d0mfcp)    
    d1n2( 1:i0smax) = fd1n2(d0mfcr,d0mfcp)

    d1t0( 1:i0smax) = fd1t0(d0mfcr,d0mfcp)
    d1t2( 1:i0smax) = fd1t2(d0mfcr,d0mfcp)
    

    DO i0sidi = 1, i0smax

       !C
       !C (1/\rho)(dp_{a}/d\rho) in J/m^{3}
       !C
       
       fd1p2(i0sidi) &
            = d1n0(i0sidi)*d1t2(i0sidi)&
            + d1n2(i0sidi)*d1t0(i0sidi)
    ENDDO
    
    RETURN
    
  END FUNCTION fd1p2


  FUNCTION fd0jt(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d0iar,d0bc
    
    INTEGER(i0ikind)::i0sidi
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jt
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1p2
    REAL(   i0rkind)::d0q0

    fd0jt = 0.D0
    d1p2  = fd1p2(d0mfcr,d0mfcp)
    d0q0  = fd0q0(d0mfcr,d0mfcp)

    DO i0sidi = 1, i0smax
       fd0jt = fd0jt + d1p2(i0sidi)
    ENDDO
    
    fd0jt = -fd0jt*(1.D0 + 1.5D0*(d0iar**2)*(d0mfcr**2))*d0q0&
         &        /(d0bc*(d0iar**2))
    
    RETURN
    
  END FUNCTION fd0jt
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF j_{\para}
  !C
  !C                     2014-02-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0jb(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d0rmnr,d0bc
    
    INTEGER(i0ikind)::i0sidi
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jb
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1p2
    REAL(   i0rkind)::d0q0,d0cobt,d0bb
    
    fd0jb  = 0.D0
    d0q0   = fd0q0(  d0mfcr,d0mfcp)
    d0cobt = fd0cobt(d0mfcr,d0mfcp)
    d0bb   = fd0bb(  d0mfcr,d0mfcp)
    d1p2   = fd1p2(  d0mfcr,d0mfcp)

    
    DO i0sidi = 1, i0smax
       fd0jb = fd0jb + d1p2(i0sidi)
    ENDDO
    
    fd0jb = -fd0jb*d0cobt*d0q0/(d0bb*d0bc*(d0rmnr**2))
    
    RETURN
    
    
    RETURN

  END FUNCTION fd0jb
  
  
  FUNCTION fd2f0(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:6,1:i0smax)::fd2f0
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1t0
    REAL(   i0rkind)::d0fb,d0ft,d0t0
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1,i0smax
       IF(i0sidi.EQ.1)THEN
          d1t0 =   fd1t0(d0mfcr,d0mfcp)
          d0fb = - fd0jb(d0mfcr,d0mfcp)/d0aee
          d0ft = - fd0jt(d0mfcr,d0mfcp)/d0aee
          d0t0 = d1t0(1)
          
          fd2f0(1,i0sidi) = 0.D0
          fd2f0(2,i0sidi) = d0fb
          fd2f0(3,i0sidi) = d0ft
          fd2f0(4,i0sidi) = 0.D0
          fd2f0(5,i0sidi) = 2.5D0*d0t0*d0fb
          fd2f0(6,i0sidi) = 2.5D0*d0t0*d0ft
       ELSE
          fd2f0(1,i0sidi) = 0.D0
          fd2f0(2,i0sidi) = 0.D0
          fd2f0(3,i0sidi) = 0.D0
          fd2f0(4,i0sidi) = 0.D0
          fd2f0(5,i0sidi) = 0.D0
          fd2f0(6,i0sidi) = 0.D0
       ENDIF
    ENDDO
    
    RETURN
  
  END FUNCTION fd2f0
END MODULE T2PROF
