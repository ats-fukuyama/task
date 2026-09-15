! wiwkxm.f90

MODULE w1wkxm

  CONTAINS

    SUBROUTINE w1_wkxm
      USE w1comm
      IMPLICIT NONE
      COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CSOL
      REAL(rkind),PARAMETER:: EXPARG=80.D0
      INTEGER:: NX,N,NK,KK
      REAL(rkind):: RKV,DX
      COMPLEX(rkind):: CDSP0,CDSP1,CDSP2,CDET,CK
      COMPLEX(rkind):: CDTXX,CDTXY,CDTXZ,CDTYY,CDTYZ,CDTZZ,CDETIP

    MWID = 3*4 - 1
    MLEN = 6*NXPMAX+10

    IF(ALLOCATED(CF)) DEALLOCATE(CF)
    IF(ALLOCATED(CSOL)) DEALLOCATE(CSOL)
    ALLOCATE(CF(MWID,MLEN))
    ALLOCATE(CSOL(2,NXPMAX))

    DO NX = 1 , NXPMAX
       CDSP2 = CD0(1,NX)+CD1(1,NX)**2
       CDSP1 =-CD0(1,NX)*(CD0(3,NX)+CD0(4,NX))-CD0(2,NX)**2 &
              +2.D0*CD0(2,NX)*CD1(1,NX)*CD1(2,NX) &
              +CD0(1,NX)*CD1(2,NX)**2 &
              -CD0(3,NX)*CD1(1,NX)**2
       CDSP0 = CD0(4,NX)*(CD0(1,NX)*CD0(3,NX)+CD0(2,NX)**2)

       CDET  = SQRT(CDSP1**2-4.D0*CDSP2*CDSP0)
       CSOL(1,NX) = (-CDSP1+CDET)/(2.D0*CDSP2)
       CSOL(2,NX) = (-CDSP1-CDET)/(2.D0*CDSP2)
    END DO

    DO N = 1 , 2
       DO NX = 1 , NXPMAX
          NK = (NX-1)*4 + (N-1)*2 + 1
          KK = (NX-1)*2 + N
          DX = RKV*(XA(NX+1)-XA(NX))
          CK = SQRT(CSOL(N,NX))
          IF(ABS(AIMAG(CK)*DX).GT.EXPARG) THEN
             IF(AIMAG(CK).GE.0.D0) THEN
                CK = DCMPLX(DBLE(CK), EXPARG/DX)
             ELSE
                CK = DCMPLX(DBLE(CK),-EXPARG/DX)
             ENDIF
          ENDIF
          CSKX ( KK ) =   CK

          CDTXX = CD0(1,NX) + CSOL(N,NX)*CD2(1,NX)
          CDTXY = CD0(2,NX) + CSOL(N,NX)*CD2(2,NX)
          CDTXZ = CD1(1,NX) * CSKX( KK )
          CDTYY = CD0(3,NX) + CSOL(N,NX)*CD2(3,NX)
          CDTYZ = CD1(2,NX) * CSKX( KK )
          CDTZZ = CD0(4,NX) + CSOL(N,NX)*CD2(4,NX)

          CDETIP = 1.D0 / ( CDTXX*CDTYZ + CDTXZ*CDTXY )
          CSPX ( KK ) = (  CDTYY*CDTXZ - CDTXY*CDTYZ ) * CDETIP
          CSPZ ( KK ) =-(  CDTXY*CDTXY + CDTYY*CDTXX ) * CDETIP

          CF( 1 , NK+10 ) =  ( 1.D0 , 0.D0 )
          CF( 1 , NK+11 ) =  ( 1.D0 , 0.D0 )
          CF( 2 , NK+10 ) =  (  CD2(2,NX)*CSPX( KK ) &
                              - CD2(3,NX)          )*CSKX( KK ) &
                              - CD1(2,NX)*CSPZ( KK )
          CF( 2 , NK+11 ) = -(  CD2(2,NX)*CSPX( KK ) &
                              - CD2(3,NX)          )*CSKX( KK ) &
                              + CD1(2,NX)*CSPZ( KK )
          CF( 3 , NK+10 ) =   CSPZ( KK )
          CF( 3 , NK+11 ) = - CSPZ( KK )
          CF( 4 , NK+10 ) =   CD1(1,NX)*CSPX( KK ) &
                            + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
          CF( 4 , NK+11 ) =   CD1(1,NX)*CSPX( KK ) &
                            + CD2(4,NX)*CSPZ( KK )*CSKX( KK )
       END DO
    END DO
    DEALLOCATE(CSOL)
    RETURN
  END SUBROUTINE w1_wkxm
END MODULE w1wkxm
