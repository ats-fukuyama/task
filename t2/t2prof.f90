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
         i0smax,i0vmax,i0mmax,i1pdn1,i2crt,&
         d0rmnr,d0rmjr,d2mfc1,d2rzm, d2rzx, d2jm1,d2xvec
    
    INTEGER(i0ikind)::&
         i0vidi,i0sidi,i0midi,i0xidi,i0xid1d,i0xid2d
    
    REAL(   i0rkind)::d0mfcr,d0mfcp,d0jm1
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1n0,d1p0
    REAL(   i0rkind),DIMENSION(1:6,1:i0smax)::d2f0
    
100 FORMAT( 6E15.8)
110 FORMAT(10E15.8)
120 FORMAT( 5E15.8)
    
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
       d0jm1   = d2jm1(1,i0midi)
       
       i0xid2d = i2crt(2,i0midi)
       i0xid1d = i2crt(3,i0midi)
       
       !C PSI'
       d2xvec(1,i0xid1d) = fd0ctbp(d0mfcr,d0mfcp)*d0jm1/d0mfcst 
       !C I
       d2xvec(2,i0xid1d) = fd0cobt(d0mfcr,d0mfcp)/d0btcst 
       !C E_{\zeta}
       d2xvec(3,i0xid1d) = fd0coet(d0mfcr,d0mfcp)/d0etcst 
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
    
    RETURN
    
  END SUBROUTINE T2PROF_TOROIDAL
  
  SUBROUTINE T2RPROF(i0m0,i0n0,d0fc,d0fs,d0fw,d0rw,d0r,d0f0,d0f1)
    
    INTEGER(i0ikind),INTENT(IN )::i0m0,i0n0
    REAL(   i0rkind),INTENT(IN )::d0fc,d0fs,d0fw,d0rw,d0r
    REAL(   i0rkind),INTENT(OUT)::d0f0,d0f1
    REAL(   i0rkind)::&
         d0a1,d0a2,d0a3,d0a4,&
         d0b1,d0b2,d0b3,&
         d0m0,d0m1,d0m2,d0n0,d0n1,d0n2,d0rx,d0xx
    INTEGER(i0ikind)::&
         i0m1,i0m2,i0n1,i0n2
    
    i0m1 = i0m0-1
    i0m2 = i0m0-2
    i0n1 = i0n0-1
    i0n2 = i0n0-2
    
    d0m0 = DBLE(i0m0)
    d0m1 = DBLE(i0m1)
    d0m2 = DBLE(i0m2)
    d0n0 = DBLE(i0n0)
    d0n1 = DBLE(i0n1)
    d0n2 = DBLE(i0n2)
    
    d0xx = d0n0*d0m0*(d0fc-d0fs)
    d0rx = 1.D0 - d0rw
    
    d0a1 = 0.D0
    d0a2 = 0.D0
    d0a3 = 0.D0
    d0a4 = 0.D0
       
    IF(d0rx.GT.0.D0)THEN
       IF(i0n0.GE.2)THEN
          IF(    i0m0.GE.3)THEN
             d0b1 =   0.D0
             d0b2 =   0.D0
             d0b3 =  (d0fs-d0fw)/(d0rx**2)
          ELSEIF(i0m0.EQ.2)THEN
             d0b1 =   0.D0
             d0b2 =   d0n0*d0m1*d0xx
             d0b3 =  (d0fs-d0fw)/(d0rx**2)
          ELSEIF(i0m0.EQ.1)THEN
             d0b1 = - d0xx/d0rx
             d0b2 = - d0n1*d0xx
             d0b3 =  (d0fs-d0fw)/(d0rx**2)
          ELSE
             WRITE(6,*)'WRONG I0M'
             STOP
          ENDIF
       ELSE
          WRITE(6,*)'WRONG I0N'
          STOP
       ENDIF
       
       d0a1 = 0.D0
       d0a2 =  -3.0D0*d0b1 + 0.5D0*d0b2 + 6.0D0*d0b3
       d0a3 = ( 5.0D0*d0b1 - 1.0D0*d0b2 - 8.0D0*d0b3)/ d0rx
       d0a4 = (-2.0D0*d0b1 + 0.5D0*d0b2 + 3.0D0*d0b3)/(d0rx**2)
    ENDIF
    
    IF(    (d0r.GE.0.D0).AND.(d0r.LE.1.D0))THEN
       d0f0 = (d0fc-d0fs)*((1.D0-(d0r**i0n0))**i0m0) + d0fs
    ELSEIF(d0r.GT.1.D0)THEN
       d0rx = d0r - d0rw
       d0f0 = d0fw + d0a1*d0rx      + d0a2*(d0rx**2)&
            + d0a3*(d0rx**3) + d0a4*(d0rx**4) 
    ENDIF
    
    
    d0rx = d0r - d0rw
    
    IF(     d0r.EQ.0.D0)THEN
       d0f1 = 0.D0
    ELSEIF((d0r.GT.0.D0).AND.(d0r.LT.1.D0))THEN
       d0f1 = - d0xx*((1.D0-(d0r**i0n0))**i0m1)*(d0r**i0n1)
    ELSEIF((d0r.GE.1.D0).AND.(d0r.LT.d0rw))THEN
       d0f1 = d0a1 + 2.D0*d0a2*d0rx + 3.D0*d0a3*(d0rx**2)&
            + 4.D0*d0a4*(d0rx**3)
    ELSEIF(d0r.GE.d0rw)THEN
       d0f1 = 0.D0
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2RPROF

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
    
    IF((d0mfcr.GE.0.D0).AND.(d0mfcr.LE.1.D0))THEN
       fd0q1 = 2.D0*(d0qs-d0qc)*d0mfcr
    ELSEIF(d0mfcr.GT.1.D0)THEN
       fd0q1 = 2.D0*(d0qs-d0qc)*d0mfcr
    ELSE
       WRITE(6,*)'WRONG RHO INPUT'
       STOP
    ENDIF
    
    RETURN
    
  END FUNCTION fd0q1
  
  FUNCTION fd0cobt(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0cobt
    
    fd0cobt = d0rmjr*d0bc
    
    RETURN
    
  END FUNCTION fd0cobt
  
  FUNCTION fd0ctbp(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0ctbp
    REAL(i0rkind)::d0rzcr,d0cobt,d0q0
    
    d0rzcr  = fd0rzcr(d0mfcr,d0mfcp)
    d0cobt  = fd0cobt( d0mfcr,d0mfcp)
    d0q0    = fd0q0( d0mfcr,d0mfcp)
    fd0ctbp = d0cobt/((d0rzcr**2)*d0q0)

    RETURN
    
  END FUNCTION fd0ctbp

  FUNCTION fd0bb(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0bb
    REAL(i0rkind),DIMENSION(1:5)::d1mc
    REAL(i0rkind)::d0ctbp,d0cobt,d0bb
    
    d0ctbp = fd0ctbp(d0mfcr,d0mfcp)
    d0cobt = fd0cobt(d0mfcr,d0mfcp)
    d1mc   = fd1mc(  d0mfcr,d0mfcp)
    
    d0bb = (d0ctbp**2)*(d1mc(1)**2)*d1mc(5)*d1mc(2)+(d0cobt**2)*d1mc(5)
    fd0bb = SQRT(d0bb)
    
    RETURN

  END FUNCTION fd0bb
  
  FUNCTION fd0coer(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind),DIMENSION(1:i0smax)::d1p1,d1n0
    REAL(i0rkind)::d0n0,d0p1,fd0coer
    
    d1p1 = fd1p1(d0mfcr,d0mfcp)
    d1n0 = fd1n0(d0mfcr,d0mfcp)
    d0n0 = d1n0(2) 
    d0p1 = d1p1(2) 
    
    fd0coer = d0p1/(d0aee*d0n0)
    
    RETURN
    
  END FUNCTION fd0coer
  
  FUNCTION fd0coep(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0coep
    REAL(i0rkind),DIMENSION(1:i0smax)::d1t0
    REAL(i0rkind)::&
         d0jb,d0jt,d0ctbp,d0cobt,d0bb,d0n0,d0t0,&
         d0rzcr,d0ctbpi
    
    d0jb   = fd0jb(  d0mfcr,d0mfcp)
    d0jt   = fd0jt(  d0mfcr,d0mfcp)
    
    d0cobt   = fd0cobt(d0mfcr,d0mfcp)
    d0ctbp   = fd0ctbp(d0mfcr,d0mfcp)
    d0bb     = fd0bb(  d0mfcr,d0mfcp)
    
    d1t0   = fd1t0(  d0mfcr,d0mfcp)
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    
    d0t0   = d1t0(1)/(1.D3*d0aee) !C keV
    
    IF(ABS(d0ctbp).GT.0.D0)THEN
       d0ctbpi = 1.D0/d0ctbp
    ELSE
       d0ctbpi = 0.D0
    ENDIF
    
    fd0coep = (1.65D-9*15.D0/(SQRT(d0t0)**3))&
            * (d0jb*d0bb - (d0cobt/(d0rzcr**2))*d0jt)*d0ctbpi
    
    RETURN
    
  END FUNCTION fd0coep
  
  FUNCTION fd0coet(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d0rmjr,d0bc
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0coet
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1t0
    REAL(   i0rkind)::d0jb,d0jt,d0ctbp,d0cobt,d0t0,d0n0,d0ctbpi,d0rzcr
    INTEGER(i0ikind)::i1
    
    d0jb   = fd0jb(  d0mfcr,d0mfcp)
    d0jt   = fd0jt(  d0mfcr,d0mfcp)
    d0cobt = fd0cobt(d0mfcr,d0mfcp)
    d0ctbp = fd0ctbp(d0mfcr,d0mfcp)
    d1t0   = fd1t0(  d0mfcr,d0mfcp)
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    
    d0t0   = d1t0(1)/(1.D3*d0aee) !C keV
    
    IF(ABS(d0ctbp).GT.0.D0)THEN
       d0ctbpi = 1.D0/d0ctbp
    ELSE
       d0ctbpi = 0.D0
    ENDIF
    
    fd0coet = (1.65D-9*15.D0/(SQRT(d0t0)**3))*d0jt
    
    RETURN
    
  END FUNCTION fd0coet
  
  FUNCTION fd1n0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1n0
    REAL(   i0rkind)::d0nc,d0ns,d0nw,d0n0,d0n1
    INTEGER(i0ikind)::i1
    fd1n0(1:i0smax) = 0.D0
    
    DO i1 = 1,i0smax
       
       d0nc = d1nc(i1); d0ns = d1ns(i1); d0nw = d1nw(i1)
       CALL T2RPROF(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr,d0n0,d0n1)
       fd1n0(i1) = d0n0*1.D20
       
    ENDDO
    
    RETURN

  END FUNCTION fd1n0

  FUNCTION fd1n1(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1n1
    REAL(   i0rkind)::d0nc,d0ns,d0nw,d0n0,d0n1
    INTEGER(i0ikind)::i1

    fd1n1(1:i0smax) = 0.D0
    
    DO i1 = 1,i0smax
       
       d0nc = d1nc(i1); d0ns = d1ns(i1); d0nw = d1nw(i1)
       CALL T2RPROF(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr,d0n0,d0n1)
       fd1n1(i1) = d0n1*1.D20
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1n1

  FUNCTION fd1t0(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1t0
    REAL(   i0rkind)::d0tc,d0ts,d0tw,d0t0,d0t1
    INTEGER(i0ikind)::i1

    fd1t0(1:i0smax) = 0.D0
    
    DO i1 = 1,i0smax
       
       d0tc = d1tc(i1); d0ts = d1ts(i1); d0tw = d1tw(i1)
       CALL T2RPROF(1,3,d0tc,d0ts,d0tw,d0rw,d0mfcr,d0t0,d0t1)
       fd1t0(i1) = d0t0*1.D3*d0aee
       
    ENDDO
    
    RETURN

  END FUNCTION fd1t0

  FUNCTION fd1t1(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0smax,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1t1
    REAL(   i0rkind)::d0tc,d0ts,d0tw,d0t0,d0t1
    INTEGER(i0ikind)::i1

    fd1t1(1:i0smax) = 0.D0
    
    DO i1 = 1,i0smax
       
       d0tc = d1tc(i1); d0ts = d1ts(i1); d0tw = d1tw(i1)
       CALL T2RPROF(1,3,d0tc,d0ts,d0tw,d0rw,d0mfcr,d0t0,d0t1)
       fd1t1(i1) = d0t1*1.D3*d0aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1t1

  FUNCTION fd1p0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1p0,d1n0,d1t0
    INTEGER(i0ikind)::i1
    
    fd1p0(1:i0smax) = 0.D0
    d1n0( 1:i0smax) = 0.D0
    d1t0( 1:i0smax) = 0.D0
    
    d1n0( 1:i0smax) = fd1n0(d0mfcr,d0mfcp)
    d1t0( 1:i0smax) = fd1t0(d0mfcr,d0mfcp)
    
    DO i1 = 1, i0smax
       fd1p0(i1) = d1n0(i1)*d1t0(i1)
    ENDDO

    RETURN
    
  END FUNCTION fd1p0

  FUNCTION fd1p1(d0mfcr,d0mfcp)

    USE T2COMM, ONLY:i0smax
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0smax)::fd1p1
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1n0,d1n1,d1t0,d1t1
    INTEGER(i0ikind)::i1
    
    fd1p1(1:i0smax) = 0.D0

    d1n0( 1:i0smax) = 0.D0
    d1n1( 1:i0smax) = 0.D0
    d1t0( 1:i0smax) = 0.D0
    d1t1( 1:i0smax) = 0.D0

    d1n0( 1:i0smax) = fd1n0(d0mfcr,d0mfcp)    
    d1n1( 1:i0smax) = fd1n1(d0mfcr,d0mfcp)

    d1t0( 1:i0smax) = fd1t0(d0mfcr,d0mfcp)
    d1t1( 1:i0smax) = fd1t1(d0mfcr,d0mfcp)

    
    DO i1 = 1, i0smax
       fd1p1(i1) = d1n0(i1)*d1t1(i1) +  d1n1(i1)*d1t0(i1)
    ENDDO
    
    RETURN
    
  END FUNCTION fd1p1
  
  FUNCTION fd0jb(d0mfcr,d0mfcp)

    USE T2COMM, ONLY:i0smax,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jb
    REAL(   i0rkind),DIMENSION(1:i0smax)::d1p1
    REAL(   i0rkind)::&
         d0ctbp,d0cobt,d0bb,d0jt,d0jbs,d0rzcr
    INTEGER(i0ikind)::i2
    
    
    IF((d0mfcr.GE.0.D0).AND.(d0mfcr.LE.1.D0))THEN
       
       d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
       d0ctbp = fd0ctbp(d0mfcr,d0mfcp)
       d0cobt = fd0cobt(d0mfcr,d0mfcp)
       d1p1   = fd1p1(  d0mfcr,d0mfcp)
       
       d0bb   = fd0bb(  d0mfcr,d0mfcp)
       d0jt   = fd0jt(  d0mfcr,d0mfcp)
       
       fd0jb = 0.D0
       
       DO i2 = 1,i0smax
          fd0jb = fd0jb + d1p1(i2)
       ENDDO
       
       fd0jb = (fd0jb*d0mfcr*d0rzcr*d0ctbp)/(d0cobt*d0bb)&
            + d0bb*d0jt/d0cobt
       
    ELSEIF(d0mfcr.GT.1.D0)THEN
       
       d0rzcr = fd0rzcr(1.D0,d0mfcp)
       d0ctbp = fd0ctbp(1.D0,d0mfcp)
       d0cobt = fd0cobt(1.D0,d0mfcp)
       d1p1   = fd1p1(  1.D0,d0mfcp)
       
       d0bb   = fd0bb(  1.D0,d0mfcp)
       d0jt   = fd0jt(  1.D0,d0mfcp)
       
       d0jbs = 0.D0
       
       DO i2 = 1,i0smax
          d0jbs = d0jbs + d1p1(i2)
       ENDDO
       
       d0jbs = (d0jbs*1.D0*d0rzcr*d0ctbp)/(d0cobt*d0bb)&
            + d0bb*d0jt/d0cobt
       
       fd0jb = d0jbs*((d0mfcr-d0rw)**2)/((1.D0-d0rw)**2) 

    ENDIF
    
    RETURN

  END FUNCTION fd0jb
  
  FUNCTION fd0jt(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0rmu0
    USE T2COMM, ONLY:i0smax,d0rw,d0rmnr
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jt
    REAL(   i0rkind)::&
         d0ctbp,d0cobt,d0bb,d0q0,d0q1,d0jts,d0r0
    
    IF((d0mfcr.GE.0.D0).AND.(d0mfcr.LE.1.D0))THEN

       d0r0   = fd0rzcr(d0mfcr,d0mfcp)
       d0ctbp = fd0ctbp(d0mfcr,d0mfcp)
       d0cobt = fd0cobt(d0mfcr,d0mfcp)
       d0q0   = fd0q0(  d0mfcr,d0mfcp)
       d0q1   = fd0q1(  d0mfcr,d0mfcp)
       
       fd0jt = d0ctbp*d0r0*(2.D0 - d0mfcr*d0q1/d0q0&
            - d0mfcr*d0rmnr*COS(d0mfcp)/d0r0)/d0rmu0
       
    ELSEIF((d0mfcr.GE.1.D0).AND.(d0mfcr.LE.d0rw))THEN

       d0r0   = fd0rzcr(1.D0,d0mfcp)
       d0ctbp = fd0ctbp(1.D0,d0mfcp)
       d0cobt = fd0cobt(1.D0,d0mfcp)
       d0q0   = fd0q0(  1.D0,d0mfcp)
       d0q1   = fd0q1(  1.D0,d0mfcp)
       
       d0jts = d0ctbp*d0r0*(2.D0 - 1.D0*d0q1/d0q0&
            - 1.D0*d0rmnr*COS(d0mfcp)/d0r0)/d0rmu0
       
       fd0jt = d0jts*((d0mfcr-d0rw)**2)/((1.D0-d0rw)**2)
    ENDIF
    
    RETURN
    
  END FUNCTION fd0jt
  
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
