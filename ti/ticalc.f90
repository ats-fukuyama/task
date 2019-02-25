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
    REAL(rkind):: RHON,temp
    REAL(rkind):: DV11,DV11M,DV11P,DDM,DDP,VVM,VVP
    REAL(rkind):: DV23,DV53,CC25,CC83
    INTEGER:: NSA,NV,NEQ,i,j,NEQ1,NSA1,NV1

    RHON=RM(NR)/RA

! --- evalutaion of metric ---

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

! --- Initialize coefficient matrix and right-hand-side vector ---

    MAT_LOCAL(1:NEQMAX,1:3*NEQMAX)=0.D0
    VEC_LOCAL(1:NEQMAX)=0.D0

! --- Calculate for each line ---

    DO NEQ=1,NEQMAX
       NSA=NSA_NEQ(NEQ)
       NV=NV_NEQ(NEQ)
       i=(NR-1)*NEQMAX+NEQ

       SELECT CASE(NV)
       CASE(1) ! density equation

!  --- time derivative ---

          MAT_LOCAL(NEQ,  NEQMAX+NEQ)=MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               +DV11
          VEC_LOCAL(NEQ)=VEC_LOCAL(NEQ) &
               +DV11*SOL_BASE(i)

!  --- radial transport (diffusion and advection) ---

          IF(NR.EQ.1) THEN
             DDM=       DDN(NEQ,NR)                      ! symmetric around 0
             DDP=0.5D0*(DDN(NEQ,NR)+DDN(NEQ,NR+1))
             VVM=       0.D0                             ! anti-symmetric 
             VVP=0.5D0*(VVN(NEQ,NR)+VVN(NEQ,NR+1))
          ELSEIF(NR.EQ.NRMAX) THEN
             DDM=0.5D0*( DDN(NEQ,NR-1)+DDN(NEQ,NR))
             DDP=0.5D0*(-DDN(NEQ,NR-1)+3.D0*DDN(NEQ,NR)) ! extrapolation
             VVM=0.5D0*( VVN(NEQ,NR-1)+VVN(NEQ,NR))
             VVP=0.5D0*(-VVN(NEQ,NR-1)+3.D0*VVN(NEQ,NR)) ! extrapolation 
          ELSE
             DDM=0.5D0*(DDN(NEQ,NR-1)+DDN(NEQ,NR  ))
             DDP=0.5D0*(DDN(NEQ,NR  )+DDN(NEQ,NR+1))
             VVM=0.5D0*(VVN(NEQ,NR-1)+VVN(NEQ,NR  ))
             VVP=0.5D0*(VVN(NEQ,NR  )+VVN(NEQ,NR+1))
          ENDIF

!  --- radial partilcle transport ---

          MAT_LOCAL(NEQ,         NEQ) = MAT_LOCAL(NEQ,         NEQ) &
               -       DV11M*DDM/(DR)**2*DT &
               + 0.5D0*DV11M*VVM/DR*DT
          MAT_LOCAL(NEQ,  NEQMAX+NEQ) = MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
               +       DV11M*DDM/(DR)**2*DT &
               + 0.5D0*DV11M*VVM/DR*DT
          IF(NR.NE.NRMAX) THEN
             MAT_LOCAL(NEQ,  NEQMAX+NEQ) = MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
                  +       DV11P*DDP/(DR)**2*DT &
                  - 0.5D0*DV11P*VVP/DR*DT
             MAT_LOCAL(NEQ,2*NEQMAX+NEQ) = MAT_LOCAL(NEQ,2*NEQMAX+NEQ) &
                  -       DV11P*DDP/(DR)**2*DT &
                  - 0.5D0*DV11P*VVP/DR*DT
          ELSE
             SELECT CASE(MODEL_BNDA(NV,NSA))
             CASE(0)      ! Refelection (no gradient)
                CONTINUE
             CASE(1)      ! Fixed valude
                MAT_LOCAL(NEQ,  NEQMAX+NEQ) = MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
                     +       DV11P*DDP/(DR)**2*DT &
                     - 0.5D0*DV11P*VVP/DR*DT
                TEMP=-       DV11P*DDP/(DR)**2*DT &
                     - 0.5D0*DV11P*VVP/DR*DT
                MAT_LOCAL(NEQ,  NEQMAX+NEQ)=MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
                     - temp
                VEC_LOCAL(NEQ)=VEC_LOCAL(NEQ) &
                     - 2.D0*temp*BND_VALUEA(NV,NSA)
!                WRITE(6,'(A,3I5,1P3E12.4)') &
!                     'bnd:',NV,NSA,MODEL_BNDA(NV,NSA),BND_VALUEA(NV,NSA), &
!                     temp,VEC_LOCAL(NEQ)
             CASE(2)      ! Fixed influx
                VEC_LOCAL(NEQ)=VEC_LOCAL(NEQ) &
                      + DV11P*BND_VALUEA(NV,NSA)/DR*DT
!                WRITE(6,'(A,3I5,1P2E12.4)') &
!                     'bnd:',NV,NSA,MODEL_BNDA(NV,NSA),BND_VALUEA(NV,NSA), &
!                     VEC_LOCAL(NEQ)
             CASE(3)      ! Fixed decay length
                MAT_LOCAL(NEQ,  NEQMAX+NEQ) = MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
                     +       DV11P*DDP/(DR)**2*DT &
                     - 0.5D0*DV11P*VVP/DR*DT
                TEMP=-       DV11P*DDP/(DR)**2*DT &
                     - 0.5D0*DV11P*VVP/DR*DT
                MAT_LOCAL(NEQ,  NEQMAX+NEQ)=MAT_LOCAL(NEQ,  NEQMAX+NEQ) &
                     + temp &
                     - temp*EXP(-DR/BND_VALUEA(NV,NEQ))
             END SELECT
          END IF

!  --- particle generation and destruction ---

          DO NEQ1=1,NEQMAX
             NSA1=NSA_NEQ(NEQ1)
             NV1=NV_NEQ(NEQ1)
             MAT_LOCAL(NEQ,NEQMAX+NEQ1)=MAT_LOCAL(NEQ,NEQMAX+NEQ1) &
                  -DV11*CCN(NSA,NSA1,NR)*DT
          END DO

!  --- particle source and sink ---

          VEC_LOCAL(NEQ) = VEC_LOCAL(NEQ) &
               + DV11*SSIN(NSA,NR)*DT

       END SELECT

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
