!C--------------------------------------------------------------------
!C 
!C  T2PROF
!C
!C
!C
!C
!C                       2014-02-22 H.SETO
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
         d0nncst,d0frcst,d0fbcst,d0ftcst,d0fpcst,&
         d0ppcst,d0qrcst,d0qbcst,d0qtcst,d0qpcst,&
         i0smax,i0xmax,i0vmax,i0mmax,i1pdn1,i2crt,&
         d0rmnr,d0rmjr,d2mfc1,d2rzm, d2rzx, d2jm1,d2xvec

  USE T2GOUT, ONLY: T2_GOUT
  USE T2CONV, ONLY: T2_CONV
    INTEGER(i0ikind)::&
         i0vidi,i0sidi,i0midi,i0xidi,i0xid1d,i0xid2d
    
    REAL(   i0rkind)::d0mfcr,d0mfcp,d0jm1,d0cobt
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::d2rn,d2rp
    REAL(   i0rkind),DIMENSION(1:8,1:i0smax)::d2f0
    
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
       
       d2jm1(1:9,i0midi)= fd1mc(d0mfcr,d0mfcp)
       i0xid2d = i2crt(2,i0midi)
       i0xid1d = i2crt(3,i0midi)
       
       !C PSI'
       d2xvec(1,i0xid1d) = fd0psip(d0mfcr,d0mfcp)/d0mfcst 
       !C I
       d2xvec(2,i0xid1d) = fd0cobt(d0mfcr,d0mfcp)/d0btcst 
       !C E_{\zeta}
       d2xvec(3,i0xid1d) = fd0coet(d0mfcr,d0mfcp)/d0etcst 
       !C \bar{E}_{\chi }
       d2xvec(4,i0xid2d) = fd0coep(d0mfcr,d0mfcp)/d0epcst
       !C E_{\rho }
       d2xvec(5,i0xid2d) = fd0coer(d0mfcr,d0mfcp)/d0ercst
       
       !C INITIAL PROFILE: Fr Fb Ft Qr Qb Qt 
       
       d2rn = fd2rn(d0mfcr,d0mfcp)
       d2rp = fd2rp(d0mfcr,d0mfcp)
       d2f0 = fd2f0(d0mfcr,d0mfcp)
       
       DO i0sidi = 1,i0smax
          i0vidi = 10*i0sidi - 5
          d2xvec(i0vidi+ 1,i0xid2d) = d2rn(1,i0sidi)/d0nncst
          d2xvec(i0vidi+ 2,i0xid2d) = d2f0(1,i0sidi)/d0frcst
          d2xvec(i0vidi+ 3,i0xid2d) = d2f0(2,i0sidi)/d0fbcst
          d2xvec(i0vidi+ 4,i0xid2d) = d2f0(3,i0sidi)/d0ftcst
          d2xvec(i0vidi+ 5,i0xid2d) = d2f0(4,i0sidi)/d0fpcst
          d2xvec(i0vidi+ 6,i0xid2d) = d2rp(1,i0sidi)/d0ppcst
          d2xvec(i0vidi+ 7,i0xid2d) = d2f0(5,i0sidi)/d0qrcst
          d2xvec(i0vidi+ 8,i0xid2d) = d2f0(6,i0sidi)/d0qbcst
          d2xvec(i0vidi+ 9,i0xid2d) = d2f0(7,i0sidi)/d0qtcst
          d2xvec(i0vidi+10,i0xid2d) = d2f0(8,i0sidi)/d0qpcst
       ENDDO
    ENDDO
    
    DO i0midi = 1, i0mmax
       i0xidi = i2crt(2,i0midi)
       d2rzx(1,i0xidi) = d2rzm(1,i0midi)
       d2rzx(2,i0xidi) = d2rzm(2,i0midi)
    ENDDO
    CALL T2_CONV
    CALL T2_GOUT
    
    RETURN
    
  END SUBROUTINE T2PROF_TOROIDAL

  !C------------------------------------------------------------------
  !C
  !C RADIAL PROFILE CALCURATION 
  !C 
  !C for 0 < \rho < 1
  !C   f(r) = (fc-fs)*(1-\sqrt{\rho}^n)^m + fs 
  !C
  !C for 1 < \rho < rw^2
  !C   f(r) = fs
  !C                     2014-03-15 H.SETO
  !C
  !C------------------------------------------------------------------
  FUNCTION fd1rf(i0m0,i0n0,d0fc,d0fs,d0rw,d0rho)
    
    INTEGER(i0ikind),INTENT(IN )::i0m0,i0n0
    REAL(   i0rkind),INTENT(IN )::d0fc,d0fs,d0rw,d0rho
    REAL(   i0rkind),DIMENSION(1:2)::fd1rf
    REAL(   i0rkind)::&
         d0m0,d0n0,d0rr,d0rx,d0fcs
    INTEGER(i0ikind)::&
         i0m1,i0n2
    
    d0rr = SQRT(d0rho)
    
    IF((i0n0.LT.2).OR.(i0m0.LT.2))THEN 
       WRITE(6,*)'INAPPROPRIATE PROFILE PARAMETER:'
       WRITE(6,*)'n=',i0n0,'m=',i0m0
       STOP
    ENDIF
    
    i0m1  = i0m0 - 1
    i0n2  = i0n0 - 2
    d0m0  = DBLE(i0m0)
    d0n0  = DBLE(i0n0)
    d0fcs = d0fc-d0fs
    d0rx = 1-d0rr**i0n0
    

    IF(d0rr.LE.1)THEN
       fd1rf(1) = d0fcs*(d0rx**i0m0) + d0fs
       IF(i0n2.EQ.0)THEN
          fd1rf(2) = -0.5*d0m0*d0n0*d0fcs*(d0rx**i0m1)
       ELSE
          fd1rf(2) = -0.5*d0m0*d0n0*d0fcs*(d0rx**i0m1)*(d0rr**i0n2)
       ENDIF
    ELSE
       fd1rf(1) = d0fs
       fd1rf(2) = 0.D0
    END IF

    RETURN
    
  END FUNCTION fd1rf

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF R (R,\phi,Z)
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0rzcr(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY: d0rmjr,d0rmnr
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0rzcr
    
    fd0rzcr = d0rmjr + d0rmnr*SQRT(d0mfcr) * COS(d0mfcp)
    
    RETURN
    
  END FUNCTION fd0rzcr
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF Z (R,\phi, Z)
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0rzcz(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY: d0rmjr,d0rmnr
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0rzcz
    
    fd0rzcz =        - d0rmnr*SQRT(d0mfcr)* SIN(d0mfcp)
    
    RETURN

  END FUNCTION fd0rzcz
  

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF METRIC COEFFICIENTS
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd1mc(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY: d0rmjr,d0rmnr
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind),DIMENSION(1:9)::fd1mc
    REAL(i0rkind)::d0rzcr

    !C CYLINDRICAL COORDINATES AND MSCS GEOMETRIC TENSOR
    !C
    !C    MFC: MAGNETIC SURFACE COORDINATE SYSTEM [rho,chi,zeta]
    !C    RZC: CYLINDRICAL COORDINATE SYSTEM      [  R,phi,Z]
    !C
    !C    R   = R0+a*\sqrt{\rho}*COS(\chi)
    !C    phi = \zeta
    !C    Z   =   -a*\sqrt{\rho}*SIN(\chi)
    !C
    !C    FD1MC(1) : SQRT{g} 
    !C    FD1MC(2) : g_{\rho  \rho }*\rho 
    !C    FD1MC(3) : g_{\rho  \chi } 
    !C    FD1MC(4) : g_{\chi  \chi }
    !C    FD1MC(5) : g_{\zeta \zeta} 
    !C    FD1MC(6) : g^{\rho  \rho } 
    !C    FD1MC(7) : g^{\rho  \chi } 
    !C    FD1MC(8) : g^{\chi  \chi }*\rho 
    !C    FD1MC(9) : g^{\zeta \zeta} 
    !C
    
    d0rzcr   = fd0rzcr(d0mfcr,d0mfcp)
    fd1mc(1) = (d0rmnr**2)*d0rzcr*0.5D0
    fd1mc(2) = (d0rmnr**2)/4.D0
    fd1mc(3) = 0.D0 
    fd1mc(4) = (d0rmnr**2)*d0mfcr
    fd1mc(5) = d0rzcr**2
    fd1mc(6) = 4.D0*d0mfcr/(d0rmnr**2)
    fd1mc(7) = 0.D0
    fd1mc(8) = 1.D0/(d0rmnr**2)
    fd1mc(9) = 1.D0/(d0rzcr**2)
    
    RETURN
    
  END FUNCTION fd1mc
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF d\psi/d\rho
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0psip(d0mfcr,d0mfcp)
    USE T2CNST, ONLY:d0rmu0,d0aee
    USE T2COMM, ONLY:&
         i0smax,d0rmnr,d0rmjr,&
         i0nm,i0nn,i0tm,i0tn,d0rw,&
         d1nc,d1ns,d1tc,d1ts
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0psip
    REAL(   i0rkind)::&
         d0a,d0b,d0c,d0d,d0e,d0f,&
         d0coef0,d0coef1,d0coef2,d0coef3,d0coef4,d0coef5,d0coef6
    INTEGER(i0ikind)::i0sidi
    !IF(  (i0nm.NE.3).OR.(i0nn.NE.2).OR.&
    !     (i0tm.NE.3).OR.(i0tn.NE.2))THEN
    !   WRITE(6,*)'UNDERCONSTRUCTION: fd0psip'
    !   STOP
    !ENDIF
    d0a = 0.D0
    d0b = 0.D0
    d0c = 0.D0
    d0d = 0.D0
    d0e = 0.D0
    d0f = 0.D0
    
    DO i0sidi = 1,i0smax
       d0a = d0a + d1nc(i0sidi)
       d0b = d0b + d1ns(i0sidi)
       d0c = d0c + d1tc(i0sidi)
       d0d = d0d + d1ts(i0sidi)
    ENDDO
    
    d0a = d0a*1.D20
    d0b = d0b*1.D20
    d0c = d0c*1.D3*d0aee
    d0d = d0d*1.D3*d0aee
    d0e = d0rmu0*(d0rmnr**2)*(d0rmjr**2)*0.25D0
    d0f = 1.50D0*(d0rmnr**2)/(d0rmjr**2)
    
    
    d0coef0 =-d0e*3.0D0*(-2.D0*d0a*d0c+d0a*d0d+d0b*d0c)
    
    d0coef1 =-d0e*2.0D0*(&
         &    d0a*( 1.0D1*d0c-8.0D0*d0d-2.0D0*d0c*d0f+1.0D0*d0d*d0f)&
         &   +d0b*(-8.0D0*d0c+6.0D0*d0d+1.0D0*d0c*d0f              ))
    
    d0coef2 =-d0e*1.5D0*(&
         &    d0a*(-2.0D1*d0c+1.9D1*d0d+1.0D1*d0c*d0f-8.0D0*d0d*d0f)&
         &   +d0b*( 1.9D1*d0c-1.8D1*d0d-8.0D0*d0c*d0f+6.0D0*d0d*d0f))
    
    d0coef3 =-d0e*1.2D0*(&
         &    d0a*( 2.0D1*d0c-2.0D1*d0d-2.0D1*d0c*d0f+1.9D1*d0d*d0f)&
         &   +d0b*(-2.0D1*d0c+2.0D1*d0d+1.9D1*d0c*d0f-1.8D1*d0d*d0f))
    
    d0coef4 =-d0e*1.0D1       *(d0a-d0b)*(d0c-d0d)*( 2.D0*d0f-1.D0)

    d0coef5 =-d0e*(1.2D1/7.D0)*(d0a-d0b)*(d0c-d0d)*(-5.D0*d0f+1.D0)
         
    d0coef6 =-d0e*1.5D0       *(d0a-d0b)*(d0c-d0d)*( 1.D0*d0f     )
    
    IF(d0mfcr.LE.1.D0)THEN
       fd0psip = d0coef0             + d0coef1* d0mfcr     &
            &  + d0coef2*(d0mfcr**2) + d0coef3*(d0mfcr**3) &
            &  + d0coef4*(d0mfcr**4) + d0coef5*(d0mfcr**5) &
            &  + d0coef6*(d0mfcr**6)
    ELSE
       fd0psip = (d0coef0+d0coef1+d0coef2+d0coef3&
            &    +d0coef4+d0coef5+d0coef6)/(d0mfcr**2)
    ENDIF
    
    RETURN
    
  END FUNCTION fd0psip
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF I
  !C
  !C                     2014-03-15 H.SETO
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
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0bb(d0mfcr,d0mfcp)
        
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0bb
    REAL(i0rkind),DIMENSION(1:9)::d1mc
    REAL(i0rkind)::d0psip,d0cobt
    
    d0psip = fd0psip(d0mfcr,d0mfcp)
    d0cobt = fd0cobt(d0mfcr,d0mfcp)
    d1mc   = fd1mc(  d0mfcr,d0mfcp)
    
    fd0bb = (d0psip**2)*d1mc(6)+(d0cobt**2)
    fd0bb = fd0bb*d1mc(9)
    fd0bb = SQRT(fd0bb)
    
    RETURN

  END FUNCTION fd0bb
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF E_{\zeta}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0coet(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d0iar,d0rmjr
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0coet
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::d2rt,d2rp
    REAL(   i0rkind)::d0jt1d,d0psip,d0t0,d0clog,d0eta,d0temp
    INTEGER(i0ikind)::i0sidi
    d0psip = fd0psip(d0mfcr,d0mfcp)
    d2rp   = fd2rp(  d0mfcr,d0mfcp)
    d2rt   = fd2rt(  d0mfcr,d0mfcp)
    d0t0   = d2rt(1,1)/(1.D3*d0aee) !C  Electron temperature in keV
    
    d0jt1d = 0.D0
    
    DO i0sidi = 1, i0smax
       d0jt1d = d0jt1d + d2rp(2,i0sidi)
    ENDDO
    
    d0jt1d = -d0jt1d*(1.D0 + 1.5D0*(d0iar**2)*d0mfcr) &
         &  * (d0rmjr**2)/d0psip
    
    

    d0clog = 17.D0 ! Coulomb logarithm for debug 
    d0temp = (1.D0-SQRT(d0iar*SQRT(d0mfcr)))**2
    d0eta = (1.65D-9)*d0clog/(SQRT(d0t0)**3)
    d0eta = d0eta/d0temp
    
    fd0coet = d0eta*d0jt1d
    
    RETURN
    
  END FUNCTION fd0coet
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF \bar{E}_{\chi}
  !C
  !C                     2014-03-15 H.SETO
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
  !C                     2014-02-14 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0coer(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0coer
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::d2rp,d2rn
    REAL(   i0rkind)::d0n0,d0p1
    INTEGER(i0ikind)::i0sidi
    
    fd0coer = 0.D0
    
    d2rp = fd2rp(d0mfcr,d0mfcp)
    d2rn = fd2rn(d0mfcr,d0mfcp)
    
    DO i0sidi = 2, i0smax
       fd0coer = fd0coer + d2rp(2,i0sidi)
    ENDDO
    
    fd0coer = fd0coer/(d0aee*d2rn(1,1))
    
    RETURN
    
  END FUNCTION fd0coer
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF n_{a}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd2rn(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,i0nm,i0nn,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::fd2rn
    REAL(   i0rkind),DIMENSION(1:2)::d1rn
    REAL(   i0rkind)::d0nc,d0ns
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1,i0smax
       
       d0nc = d1nc(i0sidi)
       d0ns = d1ns(i0sidi)
       
       !C
       !C density in 10^{20} m^{-3}
       !C
       
       d1rn = fd1rf(i0nm,i0nn,d0nc,d0ns,d0rw,d0mfcr)
       
       !C
       !C density in m^{-3}
       !C
       
       fd2rn(1,i0sidi) = d1rn(1)*1.D20
       fd2rn(2,i0sidi) = d1rn(2)*1.D20
       
    ENDDO
    
    RETURN

  END FUNCTION fd2rn
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF T_{a}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd2rt(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,i0tm,i0tn,d1tc,d1ts,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::fd2rt
    REAL(   i0rkind),DIMENSION(1:2)::d1rt
    REAL(   i0rkind)::d0tc,d0ts,d0tw
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1, i0smax
       
       d0tc = d1tc(i0sidi)
       d0ts = d1ts(i0sidi)
       
       !C
       !C T and dT/d\rho in keV 
       !C
       
       d1rt = fd1rf(i0tm,i0tn,d0tc,d0ts,d0rw,d0mfcr)
       
       !C
       !C T and dT/d\rho in J
       !C
       
       fd2rt(1,i0sidi) = d1rt(1)*1.D3*d0aee
       fd2rt(2,i0sidi) = d1rt(2)*1.D3*d0aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd2rt

  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF p_{a}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd2rp(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::fd2rp,d2rn,d2rt
    INTEGER(i0ikind)::i0sidi
    
    d2rn = fd2rn(d0mfcr,d0mfcp)
    d2rt = fd2rt(d0mfcr,d0mfcp)
    
    DO i0sidi = 1, i0smax

       !C
       !C p_{a} in J/m^{3}
       !C
       
       fd2rp(1,i0sidi) = d2rn(1,i0sidi)*d2rt(1,i0sidi)
       fd2rp(2,i0sidi) = d2rn(1,i0sidi)*d2rt(2,i0sidi)&
            &           +d2rn(2,i0sidi)*d2rt(1,i0sidi)
    ENDDO
    
    RETURN
    
  END FUNCTION fd2rp
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF j_{\zeta}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0jt(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    INTEGER(i0ikind)::i0sidi
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jt
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::d2rp
    REAL(   i0rkind)::d0psip,d0rzcr
    
    fd0jt  = 0.D0
    
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    d0psip = fd0psip(d0mfcr,d0mfcp)
    d2rp   = fd2rp(  d0mfcr,d0mfcp)
    
    DO i0sidi = 1, i0smax
       fd0jt = fd0jt + d2rp(2,i0sidi)
    ENDDO
    
    fd0jt = -fd0jt*(d0rzcr**2)/d0psip
    
    RETURN
    
  END FUNCTION fd0jt
  
  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF j_{\para}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd0jb(d0mfcr,d0mfcp)
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jb
    REAL(   i0rkind)::d0cobt,d0bb,d0rzcr,d0jt
    
    fd0jb  = 0.D0
    
    d0cobt = fd0cobt(d0mfcr,d0mfcp)
    d0bb   = fd0bb(  d0mfcr,d0mfcp)
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    d0jt   = fd0jt(  d0mfcr,d0mfcp)
    
    fd0jb  = d0jt*d0cobt/(d0bb*(d0rzcr**2))
    
    RETURN
    
    
    RETURN

  END FUNCTION fd0jb


  !C-------------------------------------------------------------------
  !C
  !C INITITIAL PROFILE OF n_{a}\bar{u}_{a}^{\rho}, 
  !C                      n_{a}u_{a\para},
  !C                      n_{a}u_{a\zeta},
  !C                      n_{a}u_{a}^{\chi},
  !C                      \bar{Q}_{a}^{\rho}
  !C                      Q_{a\para}
  !C                      Q_{a\zeta}
  !C                      Q_{a}^{\chi}
  !C
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  FUNCTION fd2f0(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:8,1:i0smax)::fd2f0
    REAL(   i0rkind),DIMENSION(1:2,1:i0smax)::d2rt
    REAL(   i0rkind)::d0fb,d0ft,d0t0
    INTEGER(i0ikind)::i0sidi
    
    DO i0sidi = 1,i0smax
       IF(i0sidi.EQ.1)THEN
          d2rt =   fd2rt(d0mfcr,d0mfcp)
          d0fb = - fd0jb(d0mfcr,d0mfcp)/d0aee
          d0ft = - fd0jt(d0mfcr,d0mfcp)/d0aee
          d0t0 = d2rt(1,1)
          
          fd2f0(1,i0sidi) = 0.D0
          fd2f0(2,i0sidi) = d0fb
          fd2f0(3,i0sidi) = d0ft
          fd2f0(4,i0sidi) = 0.D0

          fd2f0(5,i0sidi) = 0.D0
          fd2f0(6,i0sidi) = 2.5D0*d0t0*d0fb
          fd2f0(7,i0sidi) = 2.5D0*d0t0*d0ft
          fd2f0(8,i0sidi) = 0.D0
       ELSE
          fd2f0(1:8,i0sidi) = 0.D0
       ENDIF
    ENDDO
    
    RETURN
  
  END FUNCTION fd2f0
END MODULE T2PROF
