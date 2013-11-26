MODULE TEST

  USE T2CNST
  IMPLICIT NONE
  INTEGER(i0ikind)::&
       i0vmax = 21,i0spcs=2,i0nid3=0,i0nid=1

  INTEGER(i0ikind)::&
       i0xa,i0xb
  
  REAL(   i0rkind),DIMENSION(:,:),ALLOCATABLE::&
       d2x,d2y,d2z,d2bcf
  
CONTAINS  

  SUBROUTINE ALLOCATE_TEST
    ALLOCATE(d2x(  1:2,1:2))
    ALLOCATE(d2y(  1:2,1:2))
    ALLOCATE(d2z(  1:2,1:2))
    ALLOCATE(d2bcf(1:2,1:2))
    RETURN
  END SUBROUTINE ALLOCATE_TEST
  SUBROUTINE T2CALV_WORKING_VARIABLES
    
    USE PLCOMM
    USE PLPROF
    USE T2CNST
    USE LIBT2,ONLY:integ_f
    implicit none
      REAL(   i0rkind)::&
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
       d0qf,&
       d0g1,d0g2,d0g3,d0g4,d0g5,d0err,&
       d0mfcr,d0mfcp,d0rzcr,d0rmjr=3.D0

    REAL(   i0rkind)::&
         d1guv_befor(1:21),d2jm1(1:5,1),d2ws(1:11,1),&
         d1e(1:2),d1m(1:2),d1n(1:2),d1p(1:2),d1t(1:2),&
         d2ug(1:2,1),d1vt(1:2),d2cl(1:2,1:2)
    REAL(   i0rkind)::&
         d1ur(1:2),d1up(1:2),d1ut(1:2),d1ub(1:2),d1u2(1:2),&
         d1qr(1:2),d1qp(1:2),d1qt(1:2),d1qb(1:2),d1qbp(1:2),d1ni(1:2),&
         d1c1(1:2),d1c2(1:2),d1c3(1:2),d1c4(1:2),&
         d2f1(1:2,1:2),d2f2(1:2,1:2),d2f3(1:2,1:2),d2f4(1:2,1:2),&
         d1k11ps(1:2),d1k12ps(1:2),d1k22ps(1:2),&
         d1k11bn(1:2),d1k12bn(1:2),d1k22bn(1:2)
    
    INTEGER(i0ikind)::&
         i1,j1,k1,i0vid3,i0wid
    REAL(   i0rkind)::&
         d0w1,d0w2,d0w3,d0w4,d0w5,d0w6,d0w7,d0w3i,&
         d0n, d0p, d0fb, d0qb, d0up, d0ur,d0ub,d0ut,&
         d0ni,d0pib,&
         d0k11ps,d0k12ps,d0k22ps,&
         d0k11bn,d0k12bn,d0k22bn
    REAL(   i0rkind),DIMENSION(1:i0spcs)::&
         d1mu1,d1mu2,d1mu3,d1fr,d1fb,d1ft,d1ti,d1pi,d1vti
    REAL(   i0rkind),DIMENSION(1:i0spcs,1:i0spcs)::&
         d2m00,d2m01,d2m10,d2m11,&
         d2n00,d2n01,d2n10,d2n11,&
         d2nu

    !C
    !C INITIALIZE
    !C
    CALL ALLOCATE_TEST

    d2jm1(1,1)= 3.1D-2
    d2jm1(2,1)= 1.D0
    d2jm1(3,1)= 0.D0
    d2jm1(4,1)= 1.D2
    d2jm1(5,1)= 1.D0/(3.1D0*3.1D0)

    i0vid3 = i0vmax*i0nid3
  
    d1guv_befor( 1) = 0.1
    d1guv_befor( 2) = 9.D0 
    d1guv_befor( 3) = 0.D0
    d1guv_befor( 4) = 0.D0
    d1guv_befor( 5) = 0.D0

    d1guv_befor( 6) = 1.D0
    d1guv_befor( 7) = 1.D0
    d1guv_befor( 8) = 2.D0
    d1guv_befor( 9) = 3.D0
    d1guv_befor(10) = 3.D0
    d1guv_befor(11) = 1.D0
    d1guv_befor(12) = 2.D0
    d1guv_befor(13) = 3.D0

    d1guv_befor(14) = 1.D0
    d1guv_befor(15) = 1.D0
    d1guv_befor(16) = 2.D0
    d1guv_befor(17) = 3.D0
    d1guv_befor(18) = 3.D0
    d1guv_befor(19) = 1.D0
    d1guv_befor(20) = 2.D0
    d1guv_befor(21) = 3.D0
    


    d0jm1          = d2jm1(1,i0nid) ! sqrt{g}
    d0jm2          = d2jm1(2,i0nid) ! g^{rr}
    d0jm3          = d2jm1(3,i0nid) ! g^{rp}
    d0jm4          = d2jm1(4,i0nid) ! g^{pp}
    d0jm5          = d2jm1(5,i0nid) ! g^{tt}

    !C <<<

    print*,'JACS=',d0jm1,d0jm2,d0jm3,d0jm4,d0jm5

    d0bp        = d0bpcst*d1guv_befor(i0vid3+1)
    d0bt        = d0btcst*d1guv_befor(i0vid3+2)

    !C
    !C GEOMETRICAL COEFFICIENTS
    !C


!    d0ugr    = - d0jm1*d2ug(1,i0nid)
!    d0ugp    = - d0jm1*d2ug(2,i0nid)

    d0ugr = 0.D0
    d0ugp = 0.D0

    d0bp2 = (d0bp**2)*d0jm5*(d0jm1**2)*d0jm2
    d0bt2 = (d0bt**2)*d0jm5
    d0bb  = d0bp2 + d0bt2
    d0bb  = SQRT(d0bb)
    d0bpb = d0bp/d0bb
    d0btb = d0bt/d0bb
    
    print*,'B   =',d0bp,d0bt,d0bb

    IF(i0nid.GT.7)THEN
       d0g1  =   d0jm3/d0jm2
       d0g2  =   d0bb*d0bp/d0bp2
       d0g3  = - d0bt*d0bp*d0jm5/d0bp2
       d0g4  =  (d0bt2*d0bpb)/(d0bp2*d0bb)
       d0g5  = -(d0bt2*d0bpb)/(d0bp2*d0bt)
       d0bpi = 1.D0/d0bp
    ELSE
       d0g1  = 0.D0
       d0g2  = 0.D0
       d0g3  = 0.D0
       d0g4  = 0.D0
       d0g5  = 0.D0
       d0bpi = 0.D0
    END IF
    
    DO i1 = 1, i0spcs

       i0vid3 = i0vmax*i0nid3 + 8*i1

       d1n( i1) = d0nncst*d1guv_befor(i0vid3-2)
       d1fr(i1) = d0frcst*d1guv_befor(i0vid3-1)
       d1fb(i1) = d0fbcst*d1guv_befor(i0vid3  )
       d1ft(i1) = d0ftcst*d1guv_befor(i0vid3+1)
       d1p( i1) = d0ppcst*d1guv_befor(i0vid3+2)
       d1qr(i1) = d0qrcst*d1guv_befor(i0vid3+3)
       d1qb(i1) = d0qbcst*d1guv_befor(i0vid3+4)
       d1qt(i1) = d0qtcst*d1guv_befor(i0vid3+5)

       IF((d1n(i1).GT.1.D-8).AND.(d1p(i1).GT.1.D-8))THEN
          d1ur(i1) = d1fr(i1)/d1n(i1)
          d1ub(i1) = d1fb(i1)/d1n(i1)
          d1ut(i1) = d1ft(i1)/d1n(i1)
          d1t( i1) = d1p( i1)/d1n(i1)
          d1ti(i1) = 1.D0    /d1t(i1)
          d1ni(i1) = 1.D0    /d1n(i1)
          d1pi(i1) = 1.D0    /d1p(i1)
       ELSEIF(ABS((d1n(i1)).LE.1.D-8).AND.(ABS(d1p(i1)).LE.1.D-8))THEN
          d1ur(i1) = 0.D0
          d1ub(i1) = 0.D0
          d1ut(i1) = 0.D0
          d1t( i1) = 0.D0
          d1ti(i1) = 0.D0
          d1ni(i1) = 0.D0
          d1pi(i1) = 0.D0
       ELSE
          WRITE(6,*)'NEGATIVE TEMP or DENS'
          WRITE(6,*)'i1=',i1,'N',d1n(i1),'P=',d1p(i1)
          STOP
       END IF

       IF(d1p(i1).GT.1.D-15)THEN
          d1qbp(i1) = d1qb(i1)/d1p(i1)
       ELSE
          d1qbp(i1) = 0.D0
       END IF
    ENDDO
    
    print*,'DENS=',d1n(1),d1n(2)
    print*,'TEMP=',d1t(1),d1t(2)
    print*,'PRES=',d1p(1),d1p(2)
    print*,'Vr  =',d1ur(1),d1ur(2)
    print*,'Fr  =',d1fr(1),d1fr(2)
    print*,'Vp  =',d1ub(1),d1ub(2)
    print*,'Fp  =',d1fb(1),d1fb(2)
    print*,'Vt  =',d1ut(1),d1ut(2)
    print*,'Ft  =',d1ft(1),d1ft(2)

    !C >>> cheked

    d1e( 1) = -d0aee
    d1e( 2) =  d0aee
    d1m( 1) =  d0ame
    d1m( 2) =  d0amp

    print*,'CHARGE=',d1e(1),d1e(2)
    print*,'MASS  =',d1m(1),d1m(2)

  
    d2ws(         1,i0nid) = d0bb

    DO i1 = 1, i0spcs
       i0wid = 5*i0spcs - 4
       d2ws(i0wid+1,i0nid) = d1ub(i1)
       d2ws(i0wid+2,i0nid) = d1n( i1)
       d2ws(i0wid+3,i0nid) = d1fb(i1)
       d2ws(i0wid+4,i0nid) = d1p( i1)
       d2ws(i0wid+5,i0nid) = d1qb(i1)
    ENDDO


    DO i1 = 1, i0spcs
       d0ur = d1ur(i1)
       d0ub = d1ub(i1)
       d0ut = d1ut(i1)
       d0up = d0g1*d0ur + d0g2*d0ub + d0g3*d0ut
       d1up(i1) = d0up
       d1u2(i1) = (d0ur**2)*d0jm2 + (d0up**2)*d0jm4 + (d0ut**2)/d0jm5

    ENDDO

    !C
    !C
    !C FOR COLLISION AND VISCOUS TERMS
    !C
    !C

    !C THERMAL VELOCITY [m/s]
    !C CHECKED
    DO i1 = 1, i0spcs
       d1vt(i1)= SQRT(2.D0*d1t(i1)/d1m(i1))
       d1vti(i1)= SQRT(d1m(i1)*d1ti(i1)/2.D0)
    ENDDO

    print*,'Vt =',d1vt( 1),d1vt( 2)
    print*,'Vti=',d1vti(1),d1vti(2)
    !C COULOMB LOGARITHM
    !C ref: NRL PLASMA FORMULARY 2011
    !C FOR DEBUG lnA = 15

    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs
       d2cl(i1,j1) = 15.D0
    ENDDO
    ENDDO

    !C BASIC COLLISION FREQUENCY
    !C CHECKED
    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs
       IF(d1vt(i1).GT.0.D0)THEN
          d2bcf(i1,j1)&
               = (d1n(j1)*(d1e(i1)**2)*(d1e(j1)**2)*d2cl(i1,j1))&
               / (4.D0*d0pi*(d0eps0**2)*(d1m(i1)**2)*(d1vt(i1)**3))
       ELSE
          d2bcf(i1,j1) = 0.D0
       ENDIF
    ENDDO
    ENDDO
    
    print*,'d2nu1 =',d2bcf(1,1),d2bcf(1,2),d2bcf(2,1),d2bcf(2,2)
    
    !C COLLISION TIME
    !C CHECKED

    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs
       d2nu(i1,j1) = d2bcf(i1,j1)/(0.75D0*SQRT(d0pi))
    ENDDO
    ENDDO
    
    print*,'d2nu =',d2nu(1,1),d2nu(1,2),d2nu(2,1),d2nu(2,2)
    
    !C
    !C FOR FRICTION COEFFICIENTS
    !C

    !C
    !C BRAGINSKII'S MATRIX ELEMENT OF COLLISION OPERATOR
    !C Collisional Transport in Tokamak Plasmas pp278
    !C
    !C d0w1  = x_ab = vtb/vta
    !C d0w2  = y_ab = ma /mb
    !C d0w3  = t_ab = ta /tb
    !C d0w3i = t_ba = tb /ta (for avoiding zero divide) 
    !C d0w4  = 1.D0/SQRT(1 + x_ab**2)
    
    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs
       d0w1  = d1vt(j1)*d1vti(i1)
       d0w2  = d1m(i1)/d1m(j1)
       d0w3  = d1t(i1)*d1ti(j1)
       d0w3i = d1ti(i1)*d1t(j1)
       d0w4  = 1.D0/SQRT(1.D0 + d0w1**2)

       d0w5 = -      (1.D0+d0w2)*(d0w4**3)
       d0w6 = -1.5D0*(1.D0+d0w2)*(d0w4**5)

       d2m00(i1,j1) = d0w5
       d2m01(i1,j1) = d0w6
       d2m10(i1,j1) = d0w6
       d2m11(i1,j1) = -(3.25D0+4.D0*(d0w1**2)+7.5D0*(d0w1**4))*(d0w4**5)

       d2n00(i1,j1) = - d0w5
       d2n01(j1,i1) = - d0w6*d0w1*d0w3i
       d2n10(i1,j1) = - d0w6
       d2n11(i1,j1) =   6.75D0*d0w3*(d0w1**2)*(d0w4**5)
    ENDDO
    ENDDO
    
    print*,'m00=',d2m00(1,1),d2m00(1,2),d2m00(2,1),d2m00(2,2)
    print*,'m10=',d2m10(1,1),d2m10(1,2),d2m10(2,1),d2m10(2,2)
    print*,'m01=',d2m01(1,1),d2m01(1,2),d2m01(2,1),d2m01(2,2)
    print*,'m11=',d2m11(1,1),d2m11(1,2),d2m11(2,1),d2m11(2,2)

    print*,'n00=',d2n00(1,1),d2n00(1,2),d2n00(2,1),d2n00(2,2)
    print*,'n10=',d2n10(1,1),d2n10(1,2),d2n10(2,1),d2n10(2,2)
    print*,'n01=',d2n01(1,1),d2n01(1,2),d2n01(2,1),d2n01(2,2)
    print*,'n11=',d2n11(1,1),d2n11(1,2),d2n11(2,1),d2n11(2,2)

    !C FRICTION COEFFICIENTS
    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs

       d0w1 = d2n00(i1,j1)*d2nu(i1,j1)
       d0w2 = d2n01(i1,j1)*d2nu(i1,j1)
       d0w3 = d2n10(i1,j1)*d2nu(i1,j1)
       d0w4 = d2n11(i1,j1)*d2nu(i1,j1)

       IF(j1.EQ.i1)THEN
          DO k1 = 1, i0spcs
             d0w1 = d0w1 + d2m00(i1,k1)*d2nu(i1,k1)
             d0w2 = d0w2 + d2m01(i1,k1)*d2nu(i1,k1)
             d0w3 = d0w3 + d2m10(i1,k1)*d2nu(i1,k1)
             d0w4 = d0w4 + d2m11(i1,k1)*d2nu(i1,k1)
          ENDDO
       ENDIF
       
       d0w5 = d1m(i1)*d1n(i1)
       
       d0w1 = d0w1*d0w5 ! l^{11}_{ab}
       d0w2 = d0w2*d0w5 ! l^{12}_{ab}
       d0w3 = d0w3*d0w5 ! l^{21}_{ab}
       d0w4 = d0w4*d0w5 ! l^{22}_{ab}
       
       print*,'l^ab_ij=',i1,j1,d0w1,d0w2,d0w3,d0w4
       
       d0ni  = d1ni(j1)
       d0pib = d1pi(j1)
       d0w5  = d1t(i1)/d1m(i1)
       d0w6  = (d0w1+d0w2)*d0ni
       d0w7  = -0.4D0*d0w2*d0pib
       d2f1(i1,j1) = d0w6
       d2f2(i1,j1) = d0w7
       d2f3(i1,j1) = d0w5*(2.5D0*d0w6 - (d0w3+d0w4)*d0ni )
       d2f4(i1,j1) = d0w5*(2.5D0*d0w7 +  0.4D0*d0w4*d0pib)
       
       print*,'f^ab_ij=',i1,j1,&
            d2f1(i1,j1),d2f2(i1,j1),d2f3(i1,j1),d2f4(i1,j1) 
       
    ENDDO
    ENDDO
    
    !C
    !C NEOCLASSICAL COEFFICIENTS
    !C
    
    DO j1 = 1, i0spcs
    DO i1 = 1, i0spcs 
       d2x(i1,j1) = d1vti(i1)*d1vt(j1)
       d2y(i1,j1) = d1m(  i1)/d1m( j1)
       d2z(i1,j1) = d1t(  i1)*d1ti(j1)
    ENDDO
    ENDDO

    MODELG=2
    RR =3.D0
    RA =1.D0
    !C FOR BANANA
    d0mfcr = 0.1
    IF(d0mfcr.NE.0.D0)THEN
       CALL PL_QPRF(d0mfcr,d0qf)
       d0rzcr = 3.1
       
       d0w1 = d0mfcr/d0rmjr
       d0w2 = 1.46D0*SQRT(d0w1)
       d0w3 = d0w2/(1.D0-d0w2)
       d0w4 = (d0de*d0w3*2.D0*(d0qf**2)*(d0rzcr**2))/(3.D0*(d0w1**2))
       
       
       DO i0xa=1,i0spcs
          
          CALL INTEG_F(fd0k11ps,d0k11ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k12ps,d0k12ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k22ps,d0k22ps,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k11bn,d0k11bn,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k12bn,d0k12bn,d0err,EPS=1.D-8,ILST=0)
          CALL INTEG_F(fd0k22bn,d0k22bn,d0err,EPS=1.D-8,ILST=0)
          
!          print*,'i0xa=',i0xa,'k11ps=',d0k11ps,'ERR=',d0err
!          print*,'i0xa=',i0xa,'k12ps=',d0k12ps,'ERR=',d0err
!          print*,'i0xa=',i0xa,'k22ps=',d0k22ps,'ERR=',d0err
!          print*,'i0xa=',i0xa,'k11bn=',d0k11bn,'ERR=',d0err
!          print*,'i0xa=',i0xa,'k12bn=',d0k12bn,'ERR=',d0err
!          print*,'i0xa=',i0xa,'k22bn=',d0k22bn,'ERR=',d0err          
          d0w5 = 0.4D0*d1p(i0xa)*d0de
          
          d1k11ps(i0xa) = d0w5*d0k11ps
          d1k12ps(i0xa) = d0w5*d0k12ps
          d1k22ps(i0xa) = d0w5*d0k22ps
          
          d0w6 = d0w4*d1m(i0xa)*d1n(i0xa)
          
          d1k11bn(i0xa) = d0w6*d0k11bn
          d1k12bn(i0xa) = d0w6*d0k12bn
          d1k22bn(i0xa) = d0w6*d0k22bn
          
       ENDDO

       print*,d1k11ps(2)/d1p(2)*d2nu(2,2)
       print*,d1k12ps(2)/d1p(2)*d2nu(2,2)
       print*,d1k22ps(2)/d1p(2)*d2nu(2,2)
       print*,d1k11bn(2)/d1p(2)*d2nu(2,2)
       print*,d1k12bn(2)/d1p(2)*d2nu(2,2)
       print*,d1k22bn(2)/d1p(2)*d2nu(2,2)
       print*,d0w3,d0w4,d0qf
       d1mu1(1:i0spcs)= 1.D0
       d1mu2(1:i0spcs)= 1.D0
       d1mu3(1:i0spcs)= 1.D0
       
       DO i1=1,i0spcs
          d0w1 = 3.D0*(d1mu1(i1)-d1mu2(i1))*d1ni(i1)
          d0w2 = 1.2D0*d1mu2(i1)*d1pi(i1)
          d0w3 = 3.D0*(d1mu2(i1)-d1mu3(i1))*d1ni(i1)
          d0w4 = 1.2D0*d1mu3(i1)*d1pi(i1)
          d0w5 = d1t(i1)/d1m(i1)
          d1c1(i1) = d0w1
          d1c2(i1) = d0w2
          d1c3(i1) = d0w5*(2.5D0*d0w1+d0w3)
          d1c4(i1) = d0w5*(2.5D0*d0w2+d0w4)
          d1c1(i1) = 0.D0
          d1c2(i1) = 0.D0
          d1c3(i1) = 0.D0
          d1c4(i1) = 0.D0
          !       print*,d1c1(i1),d1c2(i1),d1c3(i1),d1c4(i1),d1pi(i1)
       ENDDO
    ELSE
       
    END IF
    !C
    !C
    !C
    
    
    RETURN

  END SUBROUTINE T2CALV_WORKING_VARIABLES
  
  FUNCTION fd0k11ps(d0xa2)
        
 !   USE T2COMM,ONLY:d0xa
    USE LIBT2,ONLY:func_phi,func_G
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k11ps
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2)

    d0nut = 0.D0

    DO i0xb  = 1, i0spcs

       d0x   = d2x(i0xb,i0xa)  ! x_a/x_b
       d0y   = d2y(i0xb,i0xa)  ! m_b/m_a
       d0z   = d2z(i0xa,i0xb)  ! T_a/T_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! x_b
       d0nud =      d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nus = 2.D0*d0bcf*d0z*(1.D0+d0y)*func_G(d0xb)/d0xa 
       d0nub = 2.D0*d0bcf*func_G(d0xb)/(d0xa**3)
       
       d0nut = d0nut + 2.D0*d0nus - d0nub + d0nud
    ENDDO
    
    fd0k11ps = EXP(-d0xa**2)*(d0xa**5)/d0nut
    
    RETURN
    
  END FUNCTION fd0k11ps
 
  FUNCTION fd0k12ps(d0xa2)
        
 !   USE T2COMM,ONLY:d0xa
    USE LIBT2,ONLY:func_phi,func_G
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k12ps
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2)

    d0nut = 0.D0

    DO i0xb  = 1, i0spcs

       d0x   = d2x(i0xb,i0xa)  ! x_a/x_b
       d0y   = d2y(i0xb,i0xa)  ! m_b/m_a
       d0z   = d2z(i0xa,i0xb)  ! T_a/T_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! x_b
       d0nud =      d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nus = 2.D0*d0bcf*d0z*(1.D0+d0y)*func_G(d0xb)/d0xa 
       d0nub = 2.D0*d0bcf*func_G(d0xb)/(d0xa**3)
       
       d0nut = d0nut + 2.D0*d0nus - d0nub + d0nud
    ENDDO
    
    fd0k12ps = EXP(-d0xa**2)*(d0xa**7)/d0nut
    
    RETURN
    
  END FUNCTION fd0k12ps

  FUNCTION fd0k22ps(d0xa2)
        
 !   USE T2COMM,ONLY:d0xa
    USE LIBT2,ONLY:func_phi,func_G
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k22ps
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2)

    d0nut = 0.D0

    DO i0xb  = 1, i0spcs

       d0x   = d2x(i0xb,i0xa)  ! x_a/x_b
       d0y   = d2y(i0xb,i0xa)  ! m_b/m_a
       d0z   = d2z(i0xa,i0xb)  ! T_a/T_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! x_b
       d0nud =      d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nus = 2.D0*d0bcf*d0z*(1.D0+d0y)*func_G(d0xb)/d0xa 
       d0nub = 2.D0*d0bcf*func_G(d0xb)/(d0xa**3)
       
       d0nut = d0nut + 2.D0*d0nus - d0nub + d0nud
    ENDDO
    
    fd0k22ps = EXP(-d0xa**2)*(d0xa**9)/d0nut
    
    RETURN
    
  END FUNCTION fd0k22ps

  FUNCTION fd0k11bn(d0xa2)
        
 !   USE T2COMM,ONLY:d0xa
    USE LIBT2,ONLY:func_phi,func_G
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k11bn
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0y,d0z,&
         d0nus,d0nub,d0nud,d0nut
    
    d0xa = SQRT(d0xa2)

    d0nut = 0.D0

    DO i0xb  = 1, i0spcs

       d0x   = d2x(i0xb,i0xa)  ! x_a/x_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! x_b
       d0nud = d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
      
       d0nut = d0nut + d0nud
    ENDDO
    
    fd0k11bn = EXP(-d0xa**2)*(d0xa**3)*d0nut
    
    RETURN
    
  END FUNCTION fd0k11bn

  FUNCTION fd0k12bn(d0xa2)
        
 !   USE T2COMM,ONLY:d0xa
    USE LIBT2,ONLY:func_phi,func_G
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k12bn
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0nud,d0nut
    
    d0xa = SQRT(d0xa2)

    d0nut = 0.D0

    DO i0xb  = 1, i0spcs
       d0x   = d2x(i0xb,i0xa)  ! x_a/x_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! x_b
       d0nud = d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nut = d0nut + d0nud
    ENDDO
    
    fd0k12bn = EXP(-d0xa**2)*(d0xa**5)*d0nut
    
    RETURN
    
  END FUNCTION fd0k12bn

  FUNCTION fd0k22bn(d0xa2)
        
    !   USE T2COMM,ONLY:d0xa
    USE LIBT2,ONLY:func_phi,func_G
    REAL(   i0rkind),INTENT(IN)::d0xa2
    REAL(   i0rkind)::fd0k22bn
    REAL(   i0rkind)::&
         d0xa,d0xb,d0bcf,d0x,d0nud,d0nut
    
    d0xa = SQRT(d0xa2)
    
    d0nut = 0.D0
    
    DO i0xb  = 1, i0spcs
       d0x   = d2x(i0xb,i0xa)  ! x_a/x_b
       d0bcf = d2bcf(i0xa,i0xb)
       d0xb  = d0xa*d0x ! x_b
       d0nud = d0bcf*(func_phi(d0xb)-func_G(d0xb))/(d0xa**3)
       d0nut = d0nut + d0nud
    ENDDO
    
    fd0k22bn = EXP(-d0xa**2)*(d0xa**7)*d0nut
    
    RETURN
    
  END FUNCTION fd0k22bn
END MODULE TEST

PROGRAM TEST_CALV
  USE TEST
  CALL T2CALV_WORKING_VARIABLES
  STOP
END PROGRAM TEST_CALV

