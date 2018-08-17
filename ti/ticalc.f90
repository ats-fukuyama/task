! MODULE ticalc

MODULE ticalc

  PUBLIC ti_calc

CONTAINS

  SUBROUTINE ti_calc(NR)
    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind):: RHON
    REAL(rkind):: DV11,DV11M,DV11P,DDM,DDP,VVM,VVP
    REAL(rkind):: DV23,DV53,CC25,CC83
    INTEGER:: NSA1,NV1,NEQ1

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

    RHON=RM(NR)/RA

    DD_LOCAL


    MAT_LOCAL(1:3*NEQMAX,1:NEQMAX)=0.D0
    VEC_LOCAL(1:NEQMAX)=0.D0

    DO NEQ1=1,NEQMAX
       NSA1=NSA_NEQ(NEQ1)
       NV1=NV_NEQ(NEQ1)
       SELECT CASE(NV1)
       CASE(1)
          IF(NR.EQ.1) THEN
             DDM=0.D0
             DDP=0.5D0*(DD(NEQ1,NR)+DD(NEQ1,NR+1))
             VVM=0.D0
             VVP=0.5D0*(VV(NEQ1,NR)+VV(NEQ1,NR+1))
          ELSEIF(NR.EQ.NRMAX) THEN
             DDM=0.5D0*(DD(NEQ1,NR-1)+DD(NEQ1,NR))
             DDP=0.D0
             VVM=0.5D0*(VV(NEQ1,NR-1)+VV(NEQ1,NR))
             VVP=0.D0
          ELSE
             DDM=0.5D0*(DD(NEQ1,NR-1)+DD(NEQ1,NR))
             DDP=0.5D0*(DD(NEQ1,NR)+DD(NEQ1,NR+1))
             VVM=0.5D0*(VV(NEQ1,NR-1)+VV(NEQ1,NR))
             VVP=0.5D0*(VV(NEQ1,NR)+VV(NEQ1,NR+1))
          ENDIF

          




          MAT_LOCAL(NEQ1,  NEQMAX+NEQ1) &
               =MAT_LOCAL(NEQ1,  NEQMAX+NEQ1)+DV11
          VEC_LOCAL(NEQ1)=VEC_LOCAL(NEQ1)+DV11*SOL_PREV(NEQ1)
          MAT_LOCAL(NEQ1,         NEQ1) &
               =MAT_LOCAL(NEQ1,         NEQ1)+DV11M*DDM/(DR)**2*DT &
                                             -0.5D0*DV11M*VVM/DR*DT
          MAT_LOCAL(NEQ1,  NEQMAX+NEQ1) &
               =MAT_LOCAL(NEQ1,  NEQMAX+NEQ1)-DV11M*DDM/(DR)**2*DT &
                                             -DV11P*DDP/(DR)**2*DT &
                                             -0.5D0*DV11M*VVM/DR*DT &
                                             +0.5D0*DV11P*VVP/DR*DT
          MAT_LOCAL(NEQ1,2*NEQMAX+NEQ1) &
               =MAT_LOCAL(NEQ1,2*NEQMAX+NEQ1)+DV11P*DDP/(DR)**2*DT &
                                             +0.5D0*DV11P*VVP/DR*DT
          VEC_LOCAL(NEQ1)=VEC_LOCAL(NEQ1)+DV11*SSIN(NSA1,NR)
       END SELECT
    END DO
    RETURN
  END SUBROUTINE ti_calc
END MODULE ticalc
