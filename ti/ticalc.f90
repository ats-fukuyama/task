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
    INTEGER:: NSA,NV,NEQ,i,j,NEQ1,NSA1,NV1

    RHON=RM(NR)/RA

    DV11=DVRHO(NR)
    DV23=DVRHO(NR)**(2.D0/3.D0)
    DV53=DVRHO(NR)**(5.D0/3.D0)
    CC25=2.5D0
    CC83=8.D0/3.D0

    IF(NR.EQ.1) THEN
       DV11M=0.D0    ! gradient=0 at rho=0
       DV11P=0.5D0*(DVRHO(NR)+DVRHO(NR+1))
    ELSE IF(NR.EQ.NRMAX) THEN
       DV11M=0.5D0*(DVRHO(NR-1)+DVRHO(NR))
       DV11P=DVRHOS  ! evaluated at RG(NRMAX)=RA
    ELSE
       DV11M=0.5D0*(DVRHO(NR-1)+DVRHO(NR))
       DV11P=0.5D0*(DVRHO(NR)+DVRHO(NR+1))
    ENDIF

    MAT_LOCAL(1:NEQMAX,1:3*NEQMAX)=0.D0
    VEC_LOCAL(1:NEQMAX)=0.D0

    DO NEQ=1,NEQMAX
       NSA=NSA_NEQ(NEQ)
       NV=NV_NEQ(NEQ)
       i=(NR-1)*NEQMAX+NEQ

       SELECT CASE(NV)
       CASE(1)

!  --- time derivative ---

          MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               =MAT_LOCAL(NEQ,  NEQMAX+NEQ)+DV11
          VEC_LOCAL(NEQ)=VEC_LOCAL(NEQ)+DV11*SOL_BASE(i)


          IF(NR.EQ.1) THEN
             DDM=       DDN(NEQ,NR)                      ! symmetric around 0
             DDP=0.5D0*(DDN(NEQ,NR)+DDN(NEQ,NR+1))
             VVM=       0.D0                            ! anti-symmetric 
             VVP=0.5D0*(VVN(NEQ,NR)+VVN(NEQ,NR+1))
          ELSEIF(NR.EQ.NRMAX) THEN
             DDM=0.5D0*( DDN(NEQ,NR-1)+DDN(NEQ,NR))
             DDP=0.5D0*(-DDN(NEQ,NR-1)+3.D0*DDN(NEQ,NR))  ! extrapolate to NR+1
             VVM=0.5D0*( VVN(NEQ,NR-1)+VVN(NEQ,NR))
             VVP=0.5D0*(-VVN(NEQ,NR-1)+3.D0*DDN(NEQ,NR))  ! extrapolate to NR+1
          ELSE
             DDM=0.5D0*(DDN(NEQ,NR-1)+DDN(NEQ,NR  ))
             DDP=0.5D0*(DDN(NEQ,NR  )+DDN(NEQ,NR+1))
             VVM=0.5D0*(VVN(NEQ,NR-1)+VVN(NEQ,NR  ))
             VVP=0.5D0*(VVN(NEQ,NR  )+VVN(NEQ,NR+1))
          ENDIF

          IF(NR.EQ.NRMAX) THEN
             NSA=NSA_NEQ(NEQ)
             NV=NV_NEQ(NEQ)
             SELECT CASE(NV)
             CASE(1)
                VAL_EDGE(NEQ)=PNS(NS_NSA(NSA))
             CASE(2)
                VAL_EDGE(NEQ)=PUS(NS_NSA(NSA))
             CASE(3)
                VAL_EDGE(NEQ)=PTS(NS_NSA(NSA))
             END SELECT
          END IF

!  --- radial partilcle transport ---

          MAT_LOCAL(NEQ,         NEQ) = MAT_LOCAL(NEQ,         NEQ) &
               -       DV11M*DDM/(DR)**2*DT &
               + 0.5D0*DV11M*VVM/DR*DT
          MAT_LOCAL(NEQ,  NEQMAX+NEQ) = MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               +       DV11M*DDM/(DR)**2*DT &
               +       DV11P*DDP/(DR)**2*DT &
               + 0.5D0*DV11M*VVM/DR*DT &
               - 0.5D0*DV11P*VVP/DR*DT
          MAT_LOCAL(NEQ,2*NEQMAX+NEQ) = MAT_LOCAL(NEQ,2*NEQMAX+NEQ) &
               -       DV11P*DDP/(DR)**2*DT &
               - 0.5D0*DV11P*VVP/DR*DT

!  --- particle radial partilcle transport ---

          DO NEQ1=1,NEQMAX
             NSA1=NSA_NEQ(NEQ1)
             NV1=NV_NEQ(NEQ1)
             SELECT CASE(NV1)
             CASE(1) 
                MAT_LOCAL(NEQ,NEQMAX+NEQ1)=MAT_LOCAL(NEQ,NEQMAX+NEQ1) &
                     -DV11*CCN(NSA,NSA1,NR)*DT
!                IF(ABS(CCN(NSA,NSA1,NR)).GT.0.D0) THEN
!                   WRITE(6,'(A,3I5,1PE12.4)') &
!                        'NR,NSA,NSA1,CCN=',NR,NSA,NSA1,CCN(NSA,NSA1,NR)
!                END IF
             END SELECT
          END DO

!  --- particle source and sink ---

          VEC_LOCAL(NEQ) = VEC_LOCAL(NEQ) + DV11*SSIN(NSA,NR)*DT

       END SELECT

       IF(NR.EQ.NRMAX) THEN
          MAT_LOCAL(NEQ,  NEQMAX+NEQ)=MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               - MAT_LOCAL(NEQ,2*NEQMAX+NEQ)
          VEC_LOCAL(NEQ)=VEC_LOCAL(NEQ) &
               - 2.D0*MAT_LOCAL(NEQ,2*NEQMAX+NEQ)*VAL_EDGE(NEQ)
          MAT_LOCAL(NEQ,2*NEQMAX+NEQ)=0.D0
       END IF
    END DO

!    IF(NRANK.EQ.0.AND.IDEBUG.EQ.2) THEN
!       WRITE(6,'(A,I5,1P2E12.4)') 'NR,VEC = ',NR,(VEC_LOCAL(NEQ),NEQ=1,NEQMAX)
!       DO NEQ=1,NEQMAX
!          WRITE(6,'(1P6E12.4)') (MAT_LOCAL(NEQ,J),J=1,3*NEQMAX)
!       END DO
!    END IF

    RETURN
  END SUBROUTINE ti_calc
END MODULE ticalc
