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
         i0spcs, i0vmax,&
         i0nmax1,i2crt,d0rmnr,d0rmjr,&
         d2mfc1,d2rzc1,d2rzc3,d2jm1,d1guv,i1pdn1,&
         i1mfc1,i0nmax4
    
    
    INTEGER(i0ikind)::i1,i2,i0nid3,i0vid3,i0nid4,i0vid4
    REAL(   i0rkind)::d0mfcr,d0mfcp,d0jm1
    REAL(   i0rkind),DIMENSION(1:i0spcs)::d1n0,d1p0
    REAL(   i0rkind),DIMENSION(1:6,1:i0spcs)::d2f0
    
100 FORMAT( 6E15.8)
110 FORMAT(10E15.8)
120 FORMAT( 5E15.8)
    
    
    DO i1=1,i0nmax1
       
       !C INITIIALIZATION
       
       d0mfcr = 0.D0
       d0mfcp = 0.D0
           
       !C CYLINDRICAL COORDINATES
       
       d0mfcr = d2mfc1(i1,1)
       d0mfcp = d2mfc1(i1,2)
       
       d2rzc1(i1,1)  = fd0rzcr(d0mfcr,d0mfcp)
       d2rzc1(i1,2)  = fd0rzcz(d0mfcr,d0mfcp)
       
       !C CONTRAVARIANT GEOMETRIC TENSOR
       
       d2jm1(1:5,i1)= fd1mc(d0mfcr,d0mfcp)

       d0jm1 = d2jm1(1,i1)
       
       i0nid3 = i2crt(i1,2) -1
       i0vid3 = i0vmax*i0nid3
       i0nid4 = i1mfc1(i1) - 1

       i0vid4 = i0vmax*i0nid4

       d1guv(i0vid4+1) = fd0bp(d0mfcr,d0mfcp)*d0jm1/d0mfcst
       d1guv(i0vid4+2) = fd0bt(d0mfcr,d0mfcp)/d0btcst
       d1guv(i0vid4+3) = fd0et(d0mfcr,d0mfcp)/d0etcst
       
       d1guv(i0vid3+4) = fd0ep(d0mfcr,d0mfcp)/d0epcst
       d1guv(i0vid3+5) = fd0er(d0mfcr,d0mfcp)/d0ercst
       

       
       
       !C INITIAL PROFILE: Fr Fb Fb Qr Qb Qt (DIMENSIONLESS)

       d1n0 = fd1n0(d0mfcr,d0mfcp)
       d1p0 = fd1p0(d0mfcr,d0mfcp)
       d2f0 = fd2f0(d0mfcr,d0mfcp)
       
       DO i2 = 1,i0spcs

          d1guv(i0vid3+8*i2-2) = d1n0(  i2)/d0nncst
          d1guv(i0vid3+8*i2-1) = d2f0(1,i2)/d0frcst
          d1guv(i0vid3+8*i2  ) = d2f0(2,i2)/d0fbcst
          d1guv(i0vid3+8*i2+1) = d2f0(3,i2)/d0ftcst
          d1guv(i0vid3+8*i2+2) = d1p0(  i2)/d0ppcst
          d1guv(i0vid3+8*i2+3) = d2f0(4,i2)/d0qrcst
          d1guv(i0vid3+8*i2+4) = d2f0(5,i2)/d0qbcst
          d1guv(i0vid3+8*i2+5) = d2f0(6,i2)/d0qtcst
          
       ENDDO
       
    ENDDO
    
    DO i1=1,i0nmax1
       d2rzc3(i2crt(i1,2),1) = d2rzc1(i2crt(i1,1),1)
       d2rzc3(i2crt(i1,2),2) = d2rzc1(i2crt(i1,1),2)
    ENDDO
    
    !CALL T2READ
    
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

  SUBROUTINE T2READ
    USE T2COMM,ONLY:&
         i0fnum,i0spcs,i0nmax0,i0nmax3,i0xmax,i0vmax,&
         i1pdn2,i1rdn2,i0emax,d2rzc3,i3enr,d1guv,&
         d0rmjr,i0lmax,i0stm2,i1mlvl,&
         !C
         d1bp3,d1bt3,d1er3,d1ep3,d1et3,&
         d2n3, d2fr3,d2fb3,d2ft3,&
         d2p3, d2qr3,d2qb3,d2qt3,&
         d1jm1,d1jm2,d1jm3,d1jm4,d1jm5
    
    INTEGER(i0ikind)::&
         i1,i2,j2,i0nnc, i0csz, i0ctype,&
         i0rdn2,i0pdn2,i0ofst,i0mpt,i0ecnt,&
         i0mlva,i0mlvb,i0mlvc,i0vnumb,i0fnumb,i0vid
    REAL(   i0rkind)::d0qf

    CHARACTER(40)::c40fname,c40tname,c40lname
    CHARACTER(10)::c10fnumb
    CHARACTER( 2)::c2vnumb

100 FORMAT(A26)
110 FORMAT(A40)
120 FORMAT(A5)
130 FORMAT(A25)
140 FORMAT(A6,I7,1X,A5)
150 FORMAT(3E15.6)
160 FORMAT(A5,I7,I7)
170 FORMAT(4I7)
175 FORMAT(5I7)
176 FORMAT(6I7)
180 FORMAT(A10,I7)
190 FORMAT(I7)
200 FORMAT(A20)
210 FORMAT(E15.6)
    PRINT*,'T2READ'
    OPEN(i0fnum,file='T2PROF.inc')
    
    !C READ INITIAL PROFILE
    
    !C READ HEADER
    READ(i0fnum,*)
    READ(i0fnum,*)
    READ(i0fnum,*)
    READ(i0fnum,*)

    !C WRITE NUMBER OF POINTS AND DATATYPE
    READ(i0fnum,*)
    
    !C WRITE COODINATES OF EACH POINTS 
    !C 0    , R-R_0      , PHI      ,Z
    !C ................
    !C ................
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO

    !C WRITE ELEMENT - NODE RELATIONS 
   
    READ(i0fnum,*)
    
    i0ecnt=0
    
    DO i1=1,i0lmax
       i0mlva=i1mlvl(i1-1)
       i0mlvb=i1mlvl(i1)
       i0mlvc=i1mlvl(i1+1)
       i0rdn2=i1rdn2(i1)
       i0pdn2=i1pdn2(i1)
       DO i2=1,i0rdn2
       DO j2=1,i0pdn2
          i0ecnt= i0ecnt+1
          IF(    (i0mlva.EQ.0).AND.(i2.EQ.1))THEN
             READ(i0fnum,*)
          ELSEIF((i0mlvb.NE.i0mlvc).AND.&
                 (i0mlvc.NE.0).AND.&
                 (i2.EQ.i0rdn2))THEN
             READ(i0fnum,*)
          ELSE
             READ(i0fnum,*)
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    
    !C WRITE CELL TYPE: i0ctype
    !C  5: VTK TRIANGLE
    !C  9: VTK QUAD
    READ(i0fnum,*)
    DO i1=1,i0lmax
       i0rdn2=i1rdn2(i1)
       i0pdn2=i1pdn2(i1)
       DO i2=1,i0rdn2
       DO j2=1,i0pdn2
          READ(i0fnum,*)
       ENDDO
       ENDDO
    ENDDO
    
    !C WRITE SCALAR
    READ(i0fnum,*)

    !C gm1

    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO
    !C gm2
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO

    !C gm3
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO

    !C gm4
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO

    !C gm5
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO

    !C qprof
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       READ(i0fnum,*)
    ENDDO

    !C Bp
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       i0vid = i0vmax*(i1 - 1)
       READ(i0fnum,210)d1guv(i0vid+1)
    ENDDO


    !C Bt
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       i0vid = i0vmax*(i1 - 1)
       READ(i0fnum,210)d1guv(i0vid+2)
    ENDDO
    !C Er
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       i0vid = i0vmax*(i1 - 1)
       READ(i0fnum,210)d1guv(i0vid+3)
    ENDDO

    !C Ep
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       i0vid = i0vmax*(i1 - 1)
       READ(i0fnum,210)d1guv(i0vid+4)
    ENDDO

    !C Et
    READ(i0fnum,*)
    READ(i0fnum,*)
    DO i1=1,i0nmax3
       i0vid = i0vmax*(i1 - 1)
       READ(i0fnum,210)d1guv(i0vid+5)
    ENDDO

    !C N
    DO i2 = 1, i0spcs
       
       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+1)
       ENDDO
       
    ENDDO

    !C Fr
    DO i2 = 1, i0spcs

       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+2)
       ENDDO
       
    ENDDO
    
    !C Fb
    DO i2 = 1, i0spcs
       
       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+3)
       ENDDO
    ENDDO

    !C Ft
    DO i2 = 1, i0spcs

       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+4)
       ENDDO

    ENDDO

    !C P
    DO i2 = 1, i0spcs

       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+5)
       ENDDO

    ENDDO

    !C Qr
    DO i2 = 1, i0spcs

       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+6)
       ENDDO

    ENDDO

    !C Qb
    DO i2 = 1, i0spcs

       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+7)
       ENDDO

    ENDDO

    !C Qt
    DO i2 = 1, i0spcs

       READ(i0fnum,*)
       READ(i0fnum,*)
       DO i1=1,i0nmax3
          i0vid = i0vmax*(i1 - 1) + 5 + 8*(i2-1)
          READ(i0fnum,210)d1guv(i0vid+8)
       ENDDO

    ENDDO

    CLOSE(i0fnum)

    RETURN

  END SUBROUTINE T2READ
  
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
    !C    FD2MC(1) : SQRT{g} 
    !C    FD2MC(2) : g^{rho  rho } 
    !C    FD2MC(3) : g^{rho  chi } 
    !C    FD2MC(4) : g^{chi  chi } 
    !C    FD2MC(5) : g^{zeta zeta} 
    !C
  
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    
    IF(d0mfcr.GT.1.D-15)THEN
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
       !fd0q0 = (d0qc-d0qs)*((1.D0 - d0mfcr**2)**2)+d0qs
    ELSEIF(d0mfcr.GT.1.D0)THEN
       fd0q0 = (d0qs-d0qc)*(       d0mfcr**2)+d0qc
       !fd0q0 = d0qs
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
       !fd0q1 = 4.D0*(d0qs-d0qc)*(1.D0 - d0mfcr**2)*d0mfcr
    ELSEIF(d0mfcr.GT.1.D0)THEN
       fd0q1 = 2.D0*(d0qs-d0qc)*d0mfcr
       !d0q1 = 0.D0
    ELSE
       WRITE(6,*)'WRONG RHO INPUT'
       STOP
    ENDIF
    
    RETURN
    
  END FUNCTION fd0q1
  
  FUNCTION fd0bt(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0bt
    
    fd0bt = d0rmjr*d0bc
    
    RETURN
    
  END FUNCTION fd0bt
  
  FUNCTION fd0bp(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0bp
    REAL(i0rkind)::d0rzcr,d0bt,d0q0
    
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    d0bt   = fd0bt( d0mfcr,d0mfcp)
    d0q0   = fd0q0( d0mfcr,d0mfcp)
    fd0bp  = d0bt/((d0rzcr**2)*d0q0)
    
    RETURN
    
  END FUNCTION fd0bp

  FUNCTION fd0bb(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0spcs
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0bb
    REAL(i0rkind),DIMENSION(1:5)::d1mc
    REAL(i0rkind)::d0bp,d0bt,d0bb
    
    d0bp   = fd0bp(  d0mfcr,d0mfcp)
    d0bt   = fd0bt(  d0mfcr,d0mfcp)
    d1mc   = fd1mc(  d0mfcr,d0mfcp)
    
    d0bb = (d0bp**2)*(d1mc(1)**2)*d1mc(5)*d1mc(2)+(d0bt**2)*d1mc(5)
    fd0bb = SQRT(d0bb)

    RETURN

  END FUNCTION fd0bb
  
  FUNCTION fd0er(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0spcs,d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind),DIMENSION(1:i0spcs)::d1p1,d1n0
    REAL(i0rkind)::d0n0,d0p1,fd0er
    
    d1p1 = fd1p1(d0mfcr,d0mfcp)
    d1n0 = fd1n0(d0mfcr,d0mfcp)
    d0n0 = d1n0(2) 
    d0p1 = d1p1(2) 
    
    fd0er = d0p1/(d0aee*d0n0)
    
    RETURN
    
  END FUNCTION fd0er
  
  FUNCTION fd0ep(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0spcs,d0rmjr,d0bc
    
    REAL(i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(i0rkind)::fd0ep
    REAL(i0rkind),DIMENSION(1:i0spcs)::d1t0
    REAL(i0rkind)::&
         d0jb,d0jt,d0bp,d0bt,d0bb,d0n0,d0t0,&
         d0rzcr,d0bpi
    
    d0jb   = fd0jb(  d0mfcr,d0mfcp)
    d0jt   = fd0jt(  d0mfcr,d0mfcp)

    d0bt   = fd0bt(  d0mfcr,d0mfcp)
    d0bp   = fd0bp(  d0mfcr,d0mfcp)
    d0bb   = fd0bb(  d0mfcr,d0mfcp)
    
    d1t0   = fd1t0(  d0mfcr,d0mfcp)
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    
    d0t0   = d1t0(1)/(1.D3*d0aee) !C keV
    
    IF(ABS(d0bp).GT.0.D0)THEN
       d0bpi = 1.D0/d0bp
    ELSE
       d0bpi = 0.D0
    ENDIF
    
    fd0ep = (1.65D-9*15.D0/(SQRT(d0t0)**3))&
         * (d0jb*d0bb - (d0bt/(d0rzcr**2))*d0jt)*d0bpi
    
    RETURN
    
  END FUNCTION fd0ep
  
  FUNCTION fd0et(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0spcs,d0rmjr,d0bc
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0et
    REAL(   i0rkind),DIMENSION(1:i0spcs)::d1t0
    REAL(   i0rkind)::d0jb,d0jt,d0bp,d0bt,d0t0,d0n0,d0p1,d0bpi,d0rzcr
    INTEGER(i0ikind)::i1
    
    d0jb   = fd0jb(  d0mfcr,d0mfcp)
    d0jt   = fd0jt(  d0mfcr,d0mfcp)
    d0bt   = fd0bt(  d0mfcr,d0mfcp)
    d0bp   = fd0bp(  d0mfcr,d0mfcp)
    d1t0   = fd1t0(  d0mfcr,d0mfcp)
    d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
    
    d0t0   = d1t0(1)/(1.D3*d0aee) !C keV
    
    IF(ABS(d0bp).GT.0.D0)THEN
       d0bpi = 1.D0/d0bp
    ELSE
       d0bpi = 0.D0
    ENDIF
    
    fd0et = (1.65D-9*15.D0/(SQRT(d0t0)**3))*d0jt
    
    RETURN
    
  END FUNCTION fd0et
  
  FUNCTION fd1n0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0spcs,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0spcs)::fd1n0
    REAL(   i0rkind)::d0nc,d0ns,d0nw,d0n0,d0n1
    INTEGER(i0ikind)::i1
    fd1n0(1:i0spcs) = 0.D0
    
    DO i1 = 1,i0spcs
       
       d0nc = d1nc(i1); d0ns = d1ns(i1); d0nw = d1nw(i1)
       CALL T2RPROF(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr,d0n0,d0n1)
       fd1n0(i1) = d0n0*1.D20
       
    ENDDO
    
    RETURN

  END FUNCTION fd1n0

  FUNCTION fd1n1(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0spcs,d1nc,d1ns,d1nw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0spcs)::fd1n1
    REAL(   i0rkind)::d0nc,d0ns,d0nw,d0n0,d0n1
    INTEGER(i0ikind)::i1

    fd1n1(1:i0spcs) = 0.D0
    
    DO i1 = 1,i0spcs
       
       d0nc = d1nc(i1); d0ns = d1ns(i1); d0nw = d1nw(i1)
       CALL T2RPROF(1,3,d0nc,d0ns,d0nw,d0rw,d0mfcr,d0n0,d0n1)
       fd1n1(i1) = d0n1*1.D20
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1n1

  FUNCTION fd1t0(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0spcs,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0spcs)::fd1t0
    REAL(   i0rkind)::d0tc,d0ts,d0tw,d0t0,d0t1
    INTEGER(i0ikind)::i1

    fd1t0(1:i0spcs) = 0.D0
    
    DO i1 = 1,i0spcs
       
       d0tc = d1tc(i1); d0ts = d1ts(i1); d0tw = d1tw(i1)
       CALL T2RPROF(1,3,d0tc,d0ts,d0tw,d0rw,d0mfcr,d0t0,d0t1)
       fd1t0(i1) = d0t0*1.D3*d0aee
       
    ENDDO
    
    RETURN

  END FUNCTION fd1t0

  FUNCTION fd1t1(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0spcs,d1tc,d1ts,d1tw,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0spcs)::fd1t1
    REAL(   i0rkind)::d0tc,d0ts,d0tw,d0t0,d0t1
    INTEGER(i0ikind)::i1

    fd1t1(1:i0spcs) = 0.D0
    
    DO i1 = 1,i0spcs
       
       d0tc = d1tc(i1); d0ts = d1ts(i1); d0tw = d1tw(i1)
       CALL T2RPROF(1,3,d0tc,d0ts,d0tw,d0rw,d0mfcr,d0t0,d0t1)
       fd1t1(i1) = d0t1*1.D3*d0aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION fd1t1

  FUNCTION fd1p0(d0mfcr,d0mfcp)
    
    USE T2COMM, ONLY:i0spcs
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0spcs)::fd1p0,d1n0,d1t0
    INTEGER(i0ikind)::i1
    
    fd1p0(1:i0spcs) = 0.D0
    d1n0( 1:i0spcs) = 0.D0
    d1t0( 1:i0spcs) = 0.D0
    
    d1n0( 1:i0spcs) = fd1n0(d0mfcr,d0mfcp)
    d1t0( 1:i0spcs) = fd1t0(d0mfcr,d0mfcp)
    
    DO i1 = 1, i0spcs
       fd1p0(i1) = d1n0(i1)*d1t0(i1)
    ENDDO

    RETURN
    
  END FUNCTION fd1p0

  FUNCTION fd1p1(d0mfcr,d0mfcp)

    USE T2COMM, ONLY:i0spcs
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:i0spcs)::fd1p1
    REAL(   i0rkind),DIMENSION(1:i0spcs)::d1n0,d1n1,d1t0,d1t1
    INTEGER(i0ikind)::i1
    
    fd1p1(1:i0spcs) = 0.D0

    d1n0( 1:i0spcs) = 0.D0
    d1n1( 1:i0spcs) = 0.D0
    d1t0( 1:i0spcs) = 0.D0
    d1t1( 1:i0spcs) = 0.D0

    d1n0( 1:i0spcs) = fd1n0(d0mfcr,d0mfcp)    
    d1n1( 1:i0spcs) = fd1n1(d0mfcr,d0mfcp)

    d1t0( 1:i0spcs) = fd1t0(d0mfcr,d0mfcp)
    d1t1( 1:i0spcs) = fd1t1(d0mfcr,d0mfcp)

    
    DO i1 = 1, i0spcs
       fd1p1(i1) = d1n0(i1)*d1t1(i1) +  d1n1(i1)*d1t0(i1)
    ENDDO
    
    RETURN
    
  END FUNCTION fd1p1
  
  FUNCTION fd0jb(d0mfcr,d0mfcp)

    USE T2COMM, ONLY:i0spcs,d0rw
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jb
    REAL(   i0rkind),DIMENSION(1:i0spcs)::d1p1
    REAL(   i0rkind)::&
         d0bp,d0bt,d0bb,d0jt,d0jbs,d0rzcr
    INTEGER(i0ikind)::i2
    
    
    IF((d0mfcr.GE.0.D0).AND.(d0mfcr.LE.1.D0))THEN
       
       d0rzcr = fd0rzcr(d0mfcr,d0mfcp)
       d0bp   = fd0bp(  d0mfcr,d0mfcp)
       d0bt   = fd0bt(  d0mfcr,d0mfcp)
       d1p1   = fd1p1(  d0mfcr,d0mfcp)
       
       d0bb   = fd0bb(  d0mfcr,d0mfcp)
       d0jt   = fd0jt(  d0mfcr,d0mfcp)
       
       fd0jb = 0.D0
       
       DO i2 = 1,i0spcs
          fd0jb = fd0jb + d1p1(i2)
       ENDDO
       
       fd0jb = (fd0jb*d0mfcr*d0rzcr*d0bp)/(d0bt*d0bb)&
            + d0bb*d0jt/d0bt
       
    ELSEIF(d0mfcr.GT.1.D0)THEN
       
       d0rzcr = fd0rzcr(1.D0,d0mfcp)
       d0bp   = fd0bp(  1.D0,d0mfcp)
       d0bt   = fd0bt(  1.D0,d0mfcp)
       d1p1   = fd1p1(  1.D0,d0mfcp)
       
       d0bb   = fd0bb(  1.D0,d0mfcp)
       d0jt   = fd0jt(  1.D0,d0mfcp)
       
       d0jbs = 0.D0
       
       DO i2 = 1,i0spcs
          d0jbs = d0jbs + d1p1(i2)
       ENDDO
       
       d0jbs = (d0jbs*1.D0*d0rzcr*d0bp)/(d0bt*d0bb)&
            + d0bb*d0jt/d0bt
       
       fd0jb = d0jbs*((d0mfcr-d0rw)**2)/((1.D0-d0rw)**2) 

    ENDIF
    
    RETURN

  END FUNCTION fd0jb
  
  FUNCTION fd0jt(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0rmu0
    USE T2COMM, ONLY:i0spcs,d0rw,d0rmnr
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind)::fd0jt
    REAL(   i0rkind)::&
         d0bp,d0bt,d0bb,d0q0,d0q1,d0jts,d0r0
    
    IF((d0mfcr.GE.0.D0).AND.(d0mfcr.LE.1.D0))THEN

       d0r0 = fd0rzcr(d0mfcr,d0mfcp)
       d0bp = fd0bp(  d0mfcr,d0mfcp)
       d0bt = fd0bt(  d0mfcr,d0mfcp)
       d0q0 = fd0q0(  d0mfcr,d0mfcp)
       d0q1 = fd0q1(  d0mfcr,d0mfcp)
       
       fd0jt = d0bp*d0r0*(2.D0 - d0mfcr*d0q1/d0q0&
            - d0mfcr*d0rmnr*COS(d0mfcp)/d0r0)/d0rmu0
       
    ELSEIF((d0mfcr.GE.1.D0).AND.(d0mfcr.LE.d0rw))THEN

       d0r0 = fd0rzcr(1.D0,d0mfcp)
       d0bp = fd0bp(  1.D0,d0mfcp)
       d0bt = fd0bt(  1.D0,d0mfcp)
       d0q0 = fd0q0(  1.D0,d0mfcp)
       d0q1 = fd0q1(  1.D0,d0mfcp)
       
       d0jts = d0bp*d0r0*(2.D0 - 1.D0*d0q1/d0q0&
            - 1.D0*d0rmnr*COS(d0mfcp)/d0r0)/d0rmu0
       
       fd0jt = d0jts*((d0mfcr-d0rw)**2)/((1.D0-d0rw)**2)
    ENDIF
    
    RETURN
    
  END FUNCTION fd0jt
  
  FUNCTION fd2f0(d0mfcr,d0mfcp)
    
    USE T2CNST, ONLY:d0aee
    USE T2COMM, ONLY:i0spcs
    
    REAL(   i0rkind),INTENT(IN)::d0mfcr,d0mfcp
    REAL(   i0rkind),DIMENSION(1:6,1:i0spcs)::fd2f0
    REAL(   i0rkind),DIMENSION(1:i0spcs)::d1t0
    REAL(   i0rkind)::d0fb,d0ft,d0t0
    INTEGER(i0ikind)::i1
    DO i1 = 1,i0spcs
       IF(i1.EQ.1)THEN
          d1t0 =   fd1t0(d0mfcr,d0mfcp)
          d0fb = - fd0jb(d0mfcr,d0mfcp)/d0aee
          d0ft = - fd0jt(d0mfcr,d0mfcp)/d0aee
          d0t0 = d1t0(1)
          
          fd2f0(1,i1) = 0.D0
          fd2f0(2,i1) = d0fb
          fd2f0(3,i1) = d0ft
          fd2f0(4,i1) = 0.D0
          fd2f0(5,i1) = 2.5D0*d0t0*d0fb
          fd2f0(6,i1) = 2.5D0*d0t0*d0ft
       ELSE
          fd2f0(1,i1) = 0.D0
          fd2f0(2,i1) = 0.D0
          fd2f0(3,i1) = 0.D0
          fd2f0(4,i1) = 0.D0
          fd2f0(5,i1) = 0.D0
          fd2f0(6,i1) = 0.D0
       ENDIF
    ENDDO
  
    RETURN
  
  END FUNCTION fd2f0
END MODULE T2PROF
