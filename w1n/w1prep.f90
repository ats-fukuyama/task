MODULE w1prep

CONTAINS

!     ******* DIVISION IN X DIRECTION *******

  SUBROUTINE W1SETX(IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NXVH,NHM,NHM1,NXM,NXQ,NS,NXPMAX1,NSACM,NX,NX1,NX2,IERR,NZH,I
    REAL(rkind):: RW,FWC,WCH,WCL,XACM,DX,FVT,DXACM,WT,X,FACT,X1,DX1,DX2

    NXTMAX=NXPMAX+2*NXVMAX
    NXVH=NXVMAX/2
    NXM=NXPMAX*6+10
    NXQ=NXPMAX+1

    IF(MOD(NMODEL,2).EQ.1) THEN
       NXPMAX1=NXPMAX
    ELSE
       NXPMAX1=NXPMAX-1
    ENDIF

    RW=2.D6*PI*RF
    NSACM=0
    DO NS=1,NSMAX
       FWC=AEE*PZ(nS)*BB/(AMP*PA(nS))
       WCH=FWC/(1.D0+RA/RR)
       WCL=FWC/(1.D0-RA/RR)
       IF(ABS(WCH).LT.RW.AND.ABS(WCL).GT.RW) THEN
          XACM=RR*(FWC/RW-1.D0)
          NSACM=NS
       ENDIF
    END DO

    IF(NSACM.EQ.0.OR.DXFACT.LT.1.D-6) THEN
       DX=2.D0*RA/DBLE(NXPMAX1)
       DO NX=1,NXPMAX1+1
          XA(NX)=DBLE(NX-1)*DX-RA
       END DO
    ELSE
       FVT=AEE*1.D3/(AMP*PA(NSACM))
       FWC=AEE*PZ(NSACM)*BB/(AMP*PA(NSACM))
       DXACM=DXWDTH*SQRT(FVT*PTPP(NSACM))/FWC
       WT=2.D0*RA+DXFACT*DXACM &
                 *(ATAN((RA-XACM)/DXACM)-ATAN((-RA-XACM)/DXACM))
       X=-RA
       DO NX=1,NXPMAX1
          XA(NX)=X
          FACT=1.D0+DXFACT/(1.D0+((X-XACM)/DXACM)**2)
          X1=X+0.5D0*WT/(FACT*(NXPMAX1))
          FACT=1.D0+DXFACT/(1.D0+((X1-XACM)/DXACM)**2)
          X=X+WT/(FACT*(NXPMAX1))
       END DO
       WRITE(6,601) NSACM,XACM,DXACM,RA-X
       XA(NXPMAX1+1)=RA
    END IF

    IF(MOD(NMODEL,2).EQ.1) THEN
       DO NX=1,NXPMAX
          XAM(NX)=0.5D0*(XA(NX)+XA(NX+1))
       END DO
    ELSE
       DO NX=1,NXPMAX
          XAM(NX)=XA(NX)
       END DO
    ENDIF

    DX1=(RB-RD)/DBLE(NXVH-1)
    DX2=(RD-RA)/DBLE(NXVH-1)
    DO I=1,NXVMAX
       IF(I.LE.NXVH) THEN
          X=DX1*DBLE(I-1)-RB
       ELSE
          X=DX2*DBLE(I-NXVH-1)-RD
       ENDIF
       NX1=NXPMAX+NXVMAX+I
       XAM(NX1)=X
       NX2=NXPMAX+NXVMAX+1-I
       XAM(NX2)=-X
    END DO
    IERR=0
    RETURN

  601 FORMAT(1H ,'** DX-ACCUM. : NS,X,XWIDTH,XERR = ', &
                 I3,1P2E12.4,1PE12.4)
  END SUBROUTINE W1SETX

!     ******* DIVISION IN Z DIRECTION *******

  SUBROUTINE W1SETZ

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZH,NZ
    REAL(rkind):: AKZ0

    NZH=NZPMAX/2
    DZ=RZ/DBLE(NZPMAX)

    IF(NZPMAX.EQ.1) THEN
       ZA(1)=0.D0
       AKZ(1)=RKZ
       IF(ABS(RKZ).LT.1.D-32) RKZ=0.001D0
    ELSE
       DO NZ=1,NZPMAX
          ZA(NZ)=DZ*DBLE(NZ-1) - 0.5D0*RZ
       END DO
       AKZ0=2.D0*PI/RZ
       AKZ(1)   =0.001D0
       AKZ(NZH+1)=DBLE(NZH)*AKZ0
       DO NZ=2,NZH
          AKZ(NZ)      =DBLE(NZ-1)*AKZ0
          AKZ(NZPMAX-NZ+2)=-AKZ(NZ)
       END DO
    END IF

    RETURN

  END SUBROUTINE W1SETZ

!     ******* SETING ANTENNA CURRENT DENSITY *******

  SUBROUTINE W1ANTS

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZ, NA
    REAL(rkind):: PHASE,ANTPOS
    COMPLEX(rkind):: CPHASE

    DO NZ = 1 , NZPMAX
       CJ1( NZ ) = ( 0.D0 , 0.D0 )
       CJ2( NZ ) = ( 0.D0 , 0.D0 )
       CJ3( NZ ) = ( 0.D0 , 0.D0 )
       CJ4( NZ ) = ( 0.D0 , 0.D0 )
    END DO

    IF(NZPMAX.EQ.1) THEN
       NZANT1(1)=1
       CJ1(1)   =AJYH(1)/DZ
       NZANT2(1)=1
       CJ2(1)   =AJYL(1)/DZ
       CJ3(1)   =AJZH(1)
       CJ4(1)   =AJZL(1)
    ELSE
       DO NA = 1 , NAMAX
          ANTPOS            = RZ*ALYH(NA)/360.D0 + (RZ + DZ)*.5D0
          NZANT1(NA)        = INT( ANTPOS/DZ ) + 1
          PHASE             = APYH(NA)*PI/180.D0
          CPHASE            = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
          CJ1( NZANT1(NA) ) = AJYH( NA ) / DZ * CPHASE &
                            + CJ1( NZANT1(NA) )
          CJ3( NZANT1(NA) ) = AJZH( NA ) * CPHASE &
                            + CJ3( NZANT1(NA) )

          ANTPOS            = RZ*ALYL(NA)/360.D0 + (RZ + DZ)*.5D0
          NZANT2(NA)        = INT( ANTPOS/DZ ) + 1
          PHASE             = APYL(NA)*PI/180.D0
          CPHASE            = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
          CJ2( NZANT2(NA) ) = AJYL( NA ) / DZ * CPHASE &
                            + CJ2( NZANT2(NA) )
          CJ4( NZANT2(NA) ) = AJZL( NA ) * CPHASE &
                            + CJ4( NZANT2(NA) )
       END DO
    END IF

    RETURN
  END SUBROUTINE W1ANTS

!     ****** INTERFACE FOR FFT ******

  SUBROUTINE W1FFTL(CAF,N,KEY)
    USE libfft
    USE w1comm,ONLY: rkind
    IMPLICIT NONE

    COMPLEX(rkind),INTENT(INOUT):: CAF(N)
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CBF
    INTEGER,INTENT(IN):: N,KEY

    ALLOCATE(CBF(N))
    CALL CFFT1D(CAF,CBF,N,KEY)
    CAF(1:N)=CBF(1:N)
    DEALLOCATE(CBF)

    RETURN
  END SUBROUTINE W1FFTL

!     ****** SUBSTITUTION BY VIRUE OF SYMMETRY IN Z-DIRECTION ******

  SUBROUTINE W1SYMS(NZ)

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZ 
    INTEGER:: NZ0,NX
    INTEGER:: NS

    NZ0 = NZPMAX +2 -NZ

    IF(NSYM.EQ.1) THEN
       DO NX = 1, NXTMAX
          CE2DA(NZ,NX,1) =  CE2DA(NZ0,NX,1)
          CE2DA(NZ,NX,2) =  CE2DA(NZ0,NX,2)
          CE2DA(NZ,NX,3) = -CE2DA(NZ0,NX,3)
       END DO
    ELSE IF(NSYM.EQ.-1) THEN
       DO NX = 1, NXTMAX
          CE2DA(NZ,NX,1) = -CE2DA(NZ0,NX,1)
          CE2DA(NZ,NX,2) = -CE2DA(NZ0,NX,2)
          CE2DA(NZ,NX,3) =  CE2DA(NZ0,NX,3)
       END DO
    END IF
    PAKT(NZ,1) =  PAKT(NZ0,1)
    PAKT(NZ,2) =  PAKT(NZ0,2)
    PAKT(NZ,3) =  PAKT(NZ0,3)
    PANTK(NZ)  =  PANTK(NZ0)
    DO NS = 1, NSMAX
       PAK(NZ,NS) =  PAK(NZ0,NS)
    END DO

    RETURN
  END SUBROUTINE W1SYMS
END MODULE w1prep
