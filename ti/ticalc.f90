! MODULE ticalc

MODULE ticalc

  PRIVATE
  PUBLIC ti_calc

CONTAINS

!   --- Assemble matrix equation ---

  SUBROUTINE ti_calc(NR)
    USE ticomm
    USE ticoef,ONLY: ti_coef
    USE tisource,ONLY: ti_source
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind):: RHON
    REAL(rkind):: DV11,DV11M,DV11P,DDM,DDP,VVM,VVP
    REAL(rkind):: DV23,DV53,CC25,CC83
    INTEGER:: NSA,NV,NEQ,I

    RHON=RM(NR)/RA

    DV11=DVRHO(NR)
    DV23=DVRHO(NR)**(2.D0/3.D0)
    DV53=DVRHO(NR)**(5.D0/3.D0)
    CC25=2.5D0
    CC83=8.D0/3.D0

    IF(NR.EQ.1) THEN
       DV11M=0.D0
       DV11P=0.5D0*(DVRHO(NR)+DVRHO(NR+1))
    ELSEIF(NR.EQ.NRMAX) THEN
       DV11M=0.5D0*(DVRHO(NR-1)+DVRHO(NR))
       DV11P=0.D0
    ELSE
       DV11M=0.5D0*(DVRHO(NR-1)+DVRHO(NR))
       DV11P=0.5D0*(DVRHO(NR)+DVRHO(NR+1))
    ENDIF

    MAT_LOCAL(1:NEQMAX,1:3*NEQMAX)=0.D0
    VEC_LOCAL(1:NEQMAX)=0.D0

    DO NEQ=1,NEQMAX
       NSA=NSA_NEQ(NEQ)
       NV=NV_NEQ(NEQ)
       SELECT CASE(NV)
       CASE(1)
          IF(NR.EQ.1) THEN
             DDM=0.D0
             DDP=0.5D0*(DD(NEQ,NR)+DD(NEQ,NR+1))
             VVM=0.D0
             VVP=0.5D0*(VV(NEQ,NR)+VV(NEQ,NR+1))
          ELSEIF(NR.EQ.NRMAX) THEN
             DDM=0.5D0*(DD(NEQ,NR-1)+DD(NEQ,NR))
             DDP=0.D0
             VVM=0.5D0*(VV(NEQ,NR-1)+VV(NEQ,NR))
             VVP=0.D0
          ELSE
             DDM=0.5D0*(DD(NEQ,NR-1)+DD(NEQ,NR))
             DDP=0.5D0*(DD(NEQ,NR)+DD(NEQ,NR+1))
             VVM=0.5D0*(VV(NEQ,NR-1)+VV(NEQ,NR))
             VVP=0.5D0*(VV(NEQ,NR)+VV(NEQ,NR+1))
          ENDIF

!  --- time derivative ---

          MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               =MAT_LOCAL(NEQ,  NEQMAX+NEQ)+DV11
          VEC_LOCAL(NEQ)=VEC_LOCAL(NEQ)+DV11*SOL_ORG(NEQ)

!  --- radial partilcle transport ---

          MAT_LOCAL(NEQ,         NEQ) = MAT_LOCAL(NEQ,         NEQ) &
               +       DV11M*DDM/(DR)**2*DT &
               - 0.5D0*DV11M*VVM/DR*DT
          MAT_LOCAL(NEQ,  NEQMAX+NEQ) = MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               -       DV11M*DDM/(DR)**2*DT &
               -       DV11P*DDP/(DR)**2*DT &
               - 0.5D0*DV11M*VVM/DR*DT &
               + 0.5D0*DV11P*VVP/DR*DT
          MAT_LOCAL(NEQ,2*NEQMAX+NEQ) = MAT_LOCAL(NEQ,2*NEQMAX+NEQ) &
               +       DV11P*DDP/(DR)**2*DT &
               + 0.5D0*DV11P*VVP/DR*DT

!  --- particle source and sink ---

          VEC_LOCAL(NEQ) = VEC_LOCAL(NEQ) + DV11*SSIN(NSA,NR)

       END SELECT
    END DO

    WRITE(6,'(A,I5,1P2E12.4)') 'NR,VEC = ',NR,(VEC_LOCAL(NEQ),NEQ=1,NEQMAX)
    DO NEQ=1,NEQMAX
       WRITE(6,'(1P6E12.4)') (MAT_LOCAL(NEQ,I),I=1,3*NEQMAX)
    END DO

    RETURN
  END SUBROUTINE ti_calc
END MODULE ticalc
