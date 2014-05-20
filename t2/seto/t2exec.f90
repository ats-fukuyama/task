!C-------------------------------------------------------------------- 
!C MODULE FOR SOLVING ADVECTION-DIFFUSION PROBLEM BY FEM
!C
!C                  LAST UPDATE 2014-05-20 H.SETO
!C
!C T2INTG requires following variables:
!C
!C        rkind,ikind [from T2CNST]
!C        NNMAX,NDMAX,NVMAX                           [from T2COMM] 
!C
!C        MassScaIntgPG, AdveVecIntgPG, AdveTenIntgPG
!C        DiffTenIntgPG, GradVecIntgPG, GradTenIntgPG
!C        ExciScaIntgPG, ExciVecIntgPG, ExciTenIntgPG
!C        SourScaIntgPG                               [from T2INTG]
!C
!C        MassScaCoef, AdveVecCoef, AdveTenCoef
!C        DiffTenCoef, GradVecCoef, GradTenCoef
!C        ExciScaCoef, ExciVecCoef, ExciTenCoef
!C        SourScaCoef                                 [from T2CALV]
!C
!C
!C T2INTG provides following variables:
!C
!C        xvec: 
!C
!C  and subroutines:
!C
!C        T2_EXEC
!C        T2_EXEC_TERMINATE
!C
!C -------------------------------------------------------------------

MODULE T2EXEC
  
  USE T2CNST, ONLY:&
       i0ikind,i0rkind,ikind,rkind
  
  IMPLICIT NONE
  
  PRIVATE

  INTEGER(ikind)::&
       i0vidi,i0vidj,i0widi,i0widj,i0eidi

  INTEGER(ikind),SAVE::&
       I_v,J_v,&! variable       index
       I_k,J_k,&! known variable index
       I_e      ! element        index
  
  INTEGER(ikind),SAVE,ALLOCATABLE::& 
       NEGraphElem(:,:)   ! node-element graph table in en element       
  REAL(   rkind),SAVE::&
       JacDetLocCrd       ! jacobian of local coordinates
  REAL(   rkind),SAVE,ALLOCATABLE::&
       JacInvLocCrd(:,:),&! inverse jacobi matrix of local coordinates
       ValKnownElem(:,:),&! known variables at nodes in an element 
       BvecElem(:,:),&    ! element right hand side vector  (Ax=b)
       AmatElem(:,:,:,:)  ! element stiffness matrix        (Ax=b)

  PUBLIC T2_EXEC
  
CONTAINS
  
  !C-------------------------------------------------------------------
  !C
  !C T2EXEC (PUBLIC)
  !C FEM SOLVER FOR SIMULTANEOUS ADVECTION-DIFFUSION EQUATIONS
  !C 
  !C                2014-05-20 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2_EXEC
    
    USE T2COMM,ONLY:&
         nnmax,nmmax,nemax,nkmax,nvmax,&
         i2msvt,i2avvt,i2atvt,i2dtvt,i2gvvt,i2gtvt,&
         i2esvt,i2evvt,i2etvt,i2ssvt,&
         i3atwt,i3gtwt,i3evwt,i4etwt
    
    
    INTEGER(ikind)::&
         i0msvt,i0avvt,i0atvt,i0dtvt,i0gvvt,i0gtvt,&
         i0esvt,i0evvt,i0etvt,i0ssvt,&
         i0atwt,i0gtwt,i0evwt,i0etwt

!    LOGICAL::&
!         hasMassSca, hasAdveVec, hasAdveTen, hasDiffTen,&
!         hasGradVec, hasGradTen, hasExctSca, hasExctVec,&
!         hasExctTen, hasSourSca
    
    REAL(4)::e0time_0,e0time_1
    

    CALL CPU_TIME(e0time_0)

    DO I_e = 1, NEMAX
       
       CALL T2EXEC_SETUP_ELEMENT_VARIABLES
       
       DO J_v = 1, NVMAX
       DO I_v = 1, NVMAX
          
          i0msvt = i2msvt(I_v,J_v)
          i0avvt = i2avvt(I_v,J_v)
          i0atvt = i2atvt(I_v,J_v)
          i0dtvt = i2dtvt(I_v,J_v)
          i0gvvt = i2gvvt(I_v,J_v)
          i0gtvt = i2gtvt(I_v,J_v)
          i0esvt = i2esvt(I_v,J_v)
          i0evvt = i2evvt(I_v,J_v)
          i0etvt = i2etvt(I_v,J_v)
          i0ssvt = i2ssvt(I_v,J_v)

          IF(i0msvt.EQ.1)       CALL T2EXEC_MS
          
          IF(i0avvt.EQ.1)       CALL T2EXEC_AV
          
          IF(i0atvt.EQ.1)THEN
             DO i_k = 1, NKMAX
                i0atwt = i3atwt(I_k,I_v,J_v)
                IF(i0atwt.EQ.1) CALL T2EXEC_AT
             ENDDO
          ENDIF
          
          IF(i0dtvt.EQ.1)       CALL T2EXEC_DT
          
          IF(i0gvvt.EQ.1)       CALL T2EXEC_GV
          
          IF(i0gtvt.EQ.1)THEN
             DO I_k = 1, NKMAX
                i0gtwt = i3gtwt(I_k,I_v,J_v)
                IF(i0gtwt.EQ.1) CALL T2EXEC_GT
             ENDDO
          ENDIF
          
          IF(i0esvt.EQ.1)       CALL T2EXEC_ES
          
          IF(i0evvt.EQ.1)THEN
             DO I_k = 1, NKMAX
                i0evwt = i3evwt(I_k,I_v,J_v)
                IF(i0evwt.EQ.1) CALL T2EXEC_EV
             ENDDO
          ENDIF
          
          IF(i0etvt.EQ.1)THEN          
             DO J_k = 1, NKMAX
             DO I_k = 1, NKMAX
                i0etwt = i4etwt(I_k,J_k,I_v,J_v)
                IF(i0etwt.EQ.1) CALL T2EXEC_ET
             ENDDO
             ENDDO
          ENDIF
          
          IF(i0ssvt.EQ.1)       CALL T2EXEC_SS
          
       ENDDO
       ENDDO

       CALL T2EXEC_STORE
       
    ENDDO
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC: completed:          cpu=', &
         e0time_1-e0time_0,' [s]'

    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_BCOND
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_BCOND completed:     cpu=', &
         e0time_1-e0time_0,' [s]'

    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_SOLVE
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_SOLVE completed:     cpu=', &
         e0time_1-e0time_0,' [s]'
    
    RETURN
    
  END SUBROUTINE T2_EXEC
  
  !C------------------------------------------------------------------ 
  !C
  !C T2EXEC_SETUP_ELEMENT_VARIABLES
  !C  　
  !C    * SET LOCAL NODE-ELEMENT GRAPH
  !C    * CALCULATE JACOBIAN OF LOCAL COORDINATES
  !C    * SET KNOWN VARIABLES FOR DIFFERENTIAL
  !C    * SET STABILIZATION FACTORS FOR SUPG (killed)
  !C
  !C                2014-05-20 H.SETO
  !C
  !C------------------------------------------------------------------  
  SUBROUTINE T2EXEC_SETUP_ELEMENT_VARIABLES
    
    USE T2COMM, ONLY:&
         NNMAX,NDMAX,NVMAX,NSMAX,&
         d2ws, d2xvec,d2wrks,d2kwns,d2mtrc,&
         i3enr,d2mfc1

         
    INTEGER(ikind)::i_r,i_m,i_s,i_n,j_n,k_n
    REAL(   rkind)::gridSizeRad,gridSizePol
    
    CALL T2EXEC_PRIVATE_ALLOCATE
    
    !
    ! SET LOCAL NODE-ELEMENT GRAPH
    !
    !   NEGraphElem(1:NNMAX,1) : FOR COEFFICIENT CALCULATION
    !   NEGraphElem(1:NNMAX,2) : FOR 2D-2D 
    !   NEGraphElem(1:NNMAX,3) : FOR 1D-2D,1D-1D 
    !   NEGraphElem(1:NNMAX,4) : FOR 2D-1D
    !
    
    DO i_r = 1, 4
       DO i_n = 1, nnmax
          !NEgraphElem(i_n,i_r)=NEgraph(i_n,i_r,i_e)
          NEGraphElem(i_n,i_r)=i3enr(i_n,i_r,i_e)
       ENDDO
    ENDDO
    
    
    !C
    !C CALCULATE JACOBIAN OF LOCAL COORDINATE 
    !C 
    !C JacDetLocCrd: JACOBIAN OF LOCAL COORDINATES
    !C JacInvLocCrd: INVERSE JACOBI MATRIX OF LOCAL COORDINATES
    !C 
    
    i_n = NEGraphElem(1,1)
    j_n = NEGraphElem(2,1)
    k_n = NEGraphElem(4,1)

    gridSizeRad = d2mfc1(1,j_n)-d2mfc1(1,i_n)
    gridSizeRad  = ABS(gridSizeRad)

    gridSizePol = d2mfc1(2,k_n)-d2mfc1(2,i_n)
    gridSizePol  = ABS(gridSizePol)

    JacDetLocCrd = gridSizeRad*gridSizePol*0.25D0

    IF(JacDetLocCrd.LE.0.D0)THEN
       WRITE(6,'("ERROR:: Jacobian of Local Coords is singular")')
       STOP
    ENDIF
    
    JacInvLocCrd(1,1) = 2.D0/gridSizeRad
    JacInvLocCrd(1,2) = 0.D0
    JacInvLocCrd(2,1) = 0.D0
    JacInvLocCrd(2,2) = 2.D0/gridSizePol
    
    !
    ! STORE WORKING ARRAY FOR DIFFERENTIAL
    !
    ! ValKnownElem(:,1   ) : B    AT L-TH PICARD ITERATION
    ! ValKnownElem(:,2   ) : R    AT L-TH PICARD ITERATION
    ! ValKnownElem(:,2N+1) : Ub   AT L-TH PICARD ITERATION 
    ! ValKnownElem(:,2N+2) : Qb/P AT L-TH PICARD ITERATION 
    !
    
    DO i_n = 1, NNMAX
       i_m = NEGraphElem(i_n,1)
       !ValKnownElem(i_n,1) = ValKnown(1,i_m)
       !ValKnownElem(i_n,2) = ValKnown(2,i_m)
       ValKnownElem(i_n,1) = d2ws(1,i_m)
       ValKnownElem(i_n,2) = d2ws(2,i_m)

       DO i_s = 1, NSMAX
          I_k = 2*i_s
          !ValKnownElem(i_n,I_k+1) = ValKnown(I_k+1,i_m)
          !ValKnownElem(i_n,I_k+2) = ValKnown(I_k+2,i_m)
          ValKnownElem(i_n,I_k+1) = d2ws(I_k+1,i_m)
          ValKnownElem(i_n,I_k+2) = d2ws(I_k+2,i_m)
       ENDDO
    ENDDO
    
    !C
    
    AmatElem(1:NNMAX,1:NNMAX,1:NVMAX,1:NVMAX) = 0.D0 
    BvecElem(1:NNMAX,        1:NVMAX        ) = 0.D0
    
    RETURN
    
  END SUBROUTINE T2EXEC_SETUP_ELEMENT_VARIABLES

  !-------------------------------------------------------------------
  !
  ! CALCULATION OF SUBMATRCES FROM MASS SCALAR COEFFICIENTS
  !
  !                     2014-05-20 H.Seto
  !
  !-------------------------------------------------------------------

  SUBROUTINE T2EXEC_MS
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,Dt,Xvec,MassScaCoef
    USE T2INTG,ONLY: MassScaIntgPG
    !USE T2CALV,ONLY: MassScaCoef
    
    INTEGER(ikind)::i_n,j_n,k_n,i_m,i_b,i_x
    REAL(   rkind)::&
         massScaMatElem( 1:NNMAX,1:NNMAX),&
         massScaVecElem( 1:NNMAX        ),&
         massScaCoefElem(1:NNMAX        ),&
         valPriorElem(   1:NNMAX        )
    
    ! INITIALIZATION
    
    DO i_n = 1, NNMAX
       i_m = NEgraphElem( i_n,1)
       massScaCoefElem(   i_n            ) &
            = MassScaCoef(    i_m,I_v,J_v) &
            * JacDetLocCrd/Dt
    ENDDO
    
    SELECT CASE(J_v)
       
    CASE (1:3)   ! for FSA variables 
       DO i_n = 1, NNMAX
          i_b = NEGraphElem(i_n,4)
          valPriorElem(i_n) = Xvec(i_b,J_v)
       ENDDO
    CASE DEFAULT ! for 2D dependent variables
       DO i_n = 1, NNMAX
          i_x = NEgraphElem(i_n,2)
          valPriorElem(i_n) = Xvec(i_x,J_v)          
       ENDDO
    END SELECT
    
    ! MAIN LOOP 
    
    DO i_n = 1, NNMAX
    DO j_n = 1, NNMAX
       massScaMatElem(i_n,j_n) = 0.D0
       DO k_n = 1, NNMAX
          massScaMatElem(            i_n,j_n) &
               = massScaMatElem(     i_n,j_n) &
               + MassScaIntgPG(  k_n,i_n,j_n) &
               * massScaCoefElem(k_n        ) 
       ENDDO
    ENDDO
    ENDDO

    DO i_n = 1, NNMAX
       massScaVecElem(i_n) = 0.D0
       DO j_n = 1, NNMAX  
          massScaVecElem(       i_n    ) &
               = massScaVecElem(i_n    ) &
               + massScaMatElem(i_n,j_n) &
               * valPriorElem(      j_n)
       ENDDO
    ENDDO
    
    ! store submatrix
   
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + massScaMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    ! store subvector
    
    DO i_n = 1, NNMAX
       BvecElem(             i_n,I_v) &
            = BvecElem(      i_n,I_v) & 
            + massScaVecElem(i_n    )
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_MS
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF ADVECTION SUBMATRCES: V1
  !C
  !C                     2014-02-12 H.SETO
  !C
  !C-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_AV
    
    USE T2COMM,ONLY:&
         NNMAX,NDMAX,NVMAX,AdveVecCoef
    USE T2INTG,ONLY:&
         AdveVecIntgPG
    INTEGER(i0ikind)::&
         i_d,j_d,&
         i_n,j_n,k_n,&
         i_m
    
    REAL(   i0rkind)::&
         adveVecMatElem( 1:NNMAX,1:NNMAX),&
         adveVecCoefElem(1:NDMAX,1:NNMAX),&
         adveVecCoefLoc( 1:NDMAX,1:NNMAX)
    
    
    !C
    !C INTITIALIZATION
    !C
    
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO i_d = 1, NDMAX
          adveVecCoefElem(   i_d,i_n            ) &
               = AdveVecCoef(i_d,    I_v,J_v,i_m) &
               * JacDetLocCrd
          
       ENDDO
    ENDDO
    
    !C
    !C MAIN LOOP
    !C

    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX    
          adveVecCoefLoc(j_d,k_n) = 0.D0
          DO i_d = 1, NDMAX
             adveVecCoefLoc(            j_d,k_n) &
                  = adveVecCoefLoc(     j_d,k_n) &
                  + adveVecCoefElem(i_d,    k_n) &
                  * JacInvLocCrd(   i_d,j_d    ) 
          ENDDO
       ENDDO
    ENDDO

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       adveVecMatElem(i_n,j_n) = 0.D0
       DO k_n = 1, NNMAX
          DO j_d = 1, NDMAX
             adveVecMatElem(               i_n,j_n) &
                  = adveVecMatElem(        i_n,j_n) &
                  + AdveVecIntgPG( j_d,k_n,i_n,j_n) & 
                  * adveVecCoefLoc(j_d,k_n        )            
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C
    !C STORE SUBMATRIX
    !C
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + adveVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AV
  
  !-------------------------------------------------------------------
  !
  ! CALCULATION OF ADVECTION TENSOR SUBMATRCES: V2
  !
  !                     2014-05-20 H.SETO
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_AT
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX ,AdveTenCoef
    USE T2INTG,ONLY: AdveTenIntgPG
    !USE T2COEF,ONLY: AdveTenCoef
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         adveTenCoefElem( 1:NDMAX,1:NDMAX,1:NNMAX),&
         adveTenCoefLoc(  1:NDMAX,1:NDMAX,1:NNMAX),&
         adveTenMatElem(  1:NNMAX,1:NNMAX),&
         adveTenKnownElem(1:NNMAX)

    ! initialization
    
    DO i_n = 1, NNMAX
       adveTenKnownElem(       i_n)&
            = ValKnownElem(I_k,i_n)
    ENDDO
    
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          adveTenCoefElem(   i_d,j_d,i_n                ) &
               = AdveTenCoef(i_d,j_d,    I_k,I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
       ENDDO
    ENDDO
         
    ! main loop
    
    DO l_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          adveTenCoefLoc(k_d,l_d,l_n) = 0.D0          
          DO i_d = 1, NDMAX
          DO j_d = 1, NDMAX
             adveTenCoefLoc(                k_d,l_d,l_n) &
                  = adveTenCoefLoc(         k_d,l_d,l_n) &
                  + adveTenCoefElem(i_d,j_d,        l_n) &
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
       
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       adveTenMatElem(i_n,j_n) = 0.D0
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             adveTenMatElem(                         i_n,j_n) &
                  = adveTenMatElem(                  i_n,j_n) &
                  + AdveTenIntgPG(   k_d,l_d,k_n,l_n,i_n,j_n) &
                  * adveTenCoefLoc(  k_d,l_d,    l_n        ) &
                  * adveTenKnownElem(        k_n            )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    ! store submatrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + adveTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AT

  
  !-------------------------------------------------------------------
  !
  ! CALCULATION OF DIFFUSION SUBMATRCES: V1
  !
  !                     2014-05-20 H.Seto
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_DT
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,DiffTenCoef
    USE T2INTG,ONLY: DiffTenIntgPG
    !USE T2COEF,ONLY: DiffTenCoef
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         diffTenMatElem( 1:NNMAX,1:NNMAX),&
         diffTenCoefElem(1:NDMAX,1:NDMAX,1:NNMAX), &
         diffTenCoefLoc( 1:NDMAX,1:NDMAX,1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          diffTenCoefElem(   i_d,j_d,i_n            ) &
               = diffTenCoef(i_d,j_d,    I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
       ENDDO
    ENDDO

    !C
    !C MAIN LOOP
    !C

    DO k_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          diffTenCoefLoc(k_d,l_d,k_d) = 0.D0
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             diffTenCoefLoc(                k_d,l_d,k_n) &
                  = diffTenCoefLoc(         k_d,l_d,k_n) &
                  + diffTenCoefElem(i_d,j_d,        k_n) &
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
          END DO
          END DO     
       END DO
       END DO
    END DO

    DO i_n = 1, NNMAX
    DO j_n = 1, NNMAX
       diffTenMatElem(i_n,j_n) = 0.D0
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             diffTenMatElem(                   i_n,j_n) &
                  = diffTenMatElem(            i_n,j_n) &
                  + DiffTenIntgPG( k_d,l_d,k_n,i_n,j_n) &
                  * diffTenCoefLoc(k_d,l_d,k_n        )
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C store submatrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + diffTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_DT

  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF GRADIENT SUBMATRCES: A1
  !C
  !C                     2014-05-20 H.Seto
  !C
  !C-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GV
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,GradVecCoef
    USE T2INTG,ONLY: GradVecIntgPG
    !USE T2CALV,ONLY: GradVecCoef
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,&
         i_m
    
    REAL(   i0rkind)::&
         gradVecMatElem( 1:NNMAX,1:NNMAX),&
         gradVecCoefElem(1:NDMAX,1:NNMAX),&
         gradVecCoefLoc( 1:NDMAX,1:NNMAX)

    ! initialization
          
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO i_d = 1, NDMAX
          gradVecCoefElem(   i_d,i_n            ) &
               = GradVecCoef(i_d,    I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
    ENDDO
    
    ! main loop
    
    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX
          gradVecCoefLoc(j_d,k_n) = 0.D0
          DO i_d = 1, NDMAX
             gradVecCoefLoc(            j_d,k_n) &
                  = gradVecCoefLoc(     j_d,k_n) &
                  + gradVecCoefElem(i_d,    k_n) &
                  * JacInvLocCrd(   i_d,j_d    ) 
          ENDDO
       ENDDO
    ENDDO

    DO i_n = 1, NNMAX
    DO j_n = 1, NNMAX
       gradVecMatElem(i_n,j_n) = 0.D0
       DO k_n = 1, NNMAX
          DO j_d = 1, NDMAX
             gradVecMatElem(               i_n,j_n) &
                  = gradVecMatElem(        i_n,j_n) &
                  + GradVecIntgPG( j_d,k_n,i_n,j_n) &
                  * gradVecCoefLoc(j_d,k_n        )
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    ! store submatrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + gradVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GV
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF GRADIENT SUBMATRCES: A2
  !C
  !C                     2014-05-20 H.SETO
  !C
  !C-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GT
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,GradTenCoef
    USE T2INTG,ONLY: GradTenIntgPG
    !USE T2CALV,ONLY: GradTenCoef
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         gradTenMatElem(  1:NNMAX,1:NNMAX),&
         gradTenCoefElem( 1:NDMAX,1:NDMAX,1:NNMAX),&
         gradTenCoefLoc(  1:NDMAX,1:NDMAX,1:NNMAX),&
         gradTenKnownElem(1:NNMAX)
    
    ! INITIALIZATION
       
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          gradTenCoefElem(   i_d,j_d,i_n                ) &
               = GradTenCoef(i_d,j_d,    I_k,I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
       ENDDO
    ENDDO
    
    DO i_n = 1, NNMAX
       gradTenKnownElem(i_n) = ValKnownElem(i_n,I_k)
    ENDDO
    
    !C
    !C MAIN LOOP
    !C
    
    DO l_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          gradTenCoefLoc(k_d,l_d,l_n) = 0.D0
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             gradTenCoefLoc(                k_d,l_d,l_n) &
                  = gradTenCoefLoc(         k_d,l_d,l_n) &
                  + gradTenCoefElem(i_d,j_d,        l_n) &
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       gradTenMatElem(i_n,j_n) = 0.D0
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             gradTenMatElem(i_n,j_n)&
                  = gradTenMatElem(                  i_n,j_n) &
                  + GradTenIntgPG(   k_d,l_d,k_n,l_n,i_n,j_n) &
                  * gradTenKnownElem(        k_n            ) &
                  * gradTenCoefLoc(  k_d,l_d,    l_n        ) 
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO 
    
    ! store matrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + gradTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GT
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF EXCITATION SUBMATRCES: C1
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_ES

    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,ExciScaCoef
    USE T2INTG,ONLY: ExciScaIntgPG
    !USE T2CALV,ONLY: ExciScaCoef
         
    INTEGER(ikind)::&
         i_n,j_n,k_n,&
         i_d,&
         i_m


    REAL(   rkind)::&
         exciScaMatElem( 1:NNMAX,1:NNMAX),&
         exciScaCoefElem(1:NNMAX        )

    ! initialization
        
    DO i_n = 1,NNMAX
       i_m = NEgraphElem(i_n,1)
       exciScaCoefElem(   i_n            ) &
            = ExciScaCoef(    I_v,J_v,i_m) &
            * JacDetLocCrd
    ENDDO

    ! main loop
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       exciScaMatElem(i_n,j_n) = 0.D0
       DO k_n = 1, NNMAX
          exciScaMatElem(i_n,j_n)&
               = exciScaMatElem(     i_n,j_n) &
               + ExciScaIntgPG(  k_n,i_n,j_n) &
               * exciScaCoefElem(k_n        )
       ENDDO
    ENDDO
    ENDDO
    
    ! store matrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + exciScaMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_ES

  !------------------------------------------------------------------
  !
  !  CALCULATION OF EXCITATION SUBMATRCES
  !
  !                     2014-05-20 H.Seto
  !
  !-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_EV
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,ExciVecCoef
    USE T2INTG,ONLY: ExciVecIntgPG
    !USE T2CALV,ONLY: ExciVecCoef
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,m_n,&
         i_d,j_d,&
         i_m

    REAL(   rkind)::&
         exciVecMatElem(  1:NNMAX,1:NNMAX),&
         exciVecCoefElem( 1:NDMAX,1:NNMAX),&
         exciVecCoefLoc(  1:NDMAX,1:NNMAX),&    
         exciVecKnownElem(1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO i_d = 1, NDMAX
          exciVecCoefElem(   i_d,i_n                ) &
               = ExciVecCoef(i_d,    I_k,I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
    ENDDO
    
    DO i_n = 1, NNMAX
       exciVecKnownElem(i_n) = ValKnownElem(i_n,I_k)
    ENDDO
    
    ! main loop
    
    DO l_n = 1, NNMAX   
       DO j_d = 1, NDMAX
          exciVecCoefLoc(j_d,l_n) = 0.D0     
          DO i_d = 1, NDMAX
             exciVecCoefLoc(            j_d,l_n) &
                  = exciVecCoefLoc(     j_d,l_n) &
                  + exciVecCoefElem(i_d,    l_n) &
                  * JacInvLocCrd(   i_d,j_d    )
          ENDDO
       ENDDO
    ENDDO
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       exciVecMatElem(i_n,j_n) = 0.D0
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO j_d = 1, NDMAX
             exciVecMatElem(                   i_n,j_n  ) &
                  = exciVecMatElem(            i_n,j_n  ) &
                  + ExciVecIntgPG(   j_d,k_n,l_n,i_n,j_n) &
                  * exciVecKnownElem(    k_n            ) &
                  * exciVecCoefLoc(  j_d,    l_n        )
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    
    ! store submatrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + exciVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_EV

  !-------------------------------------------------------------------
  !
  ! CALCULATION OF EXCITATION SUBMATRCES
  !
  !                     2015-05-20 H.Seto
  !
  !-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_ET
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,ExciTenCoef
    USE T2INTG,ONLY: ExciTenIntgPG
    !USE T2CALV,ONLY: ExciTenCoef
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,m_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         exciTenMatElem(   1:NNMAX,1:NNMAX),&
         exciTenCoefElem(  1:NDMAX,1:NDMAX,1:NNMAX),&
         exciTenCoefLoc(   1:NDMAX,1:NDMAX,1:NNMAX),&    
         exciTenKnown1Elem(1:NNMAX),&
         exciTenKnown2Elem(1:NNMAX)

    ! initialization
  
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          exciTenCoefElem(   i_d,j_d,i_n)&
               = ExciTenCoef(i_d,j_d,    I_k,J_k,I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
       ENDDO
    ENDDO
    
    DO i_n = 1, NNMAX
       exciTenKnown1Elem(i_n) = ValKnownElem(i_n,I_k)
       exciTenKnown2Elem(i_n) = ValKnownElem(i_n,J_k)
    ENDDO
    
    ! main loop
    
    DO m_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          exciTenCoefLoc(k_d,l_d,m_n) = 0.D0
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             exciTenCoefLoc(                k_d,l_d,m_n) &
                  = exciTenCoefLoc(         k_d,l_d,m_n) &
                  + exciTenCoefElem(i_d,j_d,        m_n) &
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       exciTenMatElem(i_n,j_n) = 0.D0
       DO m_n = 1, NNMAX
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             exciTenMatElem(                              i_n,j_n) &
                  = exciTenMatElem(                       i_n,j_n) &
                  + ExciTenIntgPG(    k_d,l_d,k_n,l_n,m_n,i_n,j_n) &
                  * exciTenKnown1Elem(        k_n                ) &
                  * exciTenKnown2Elem(            l_n            ) &
                  * exciTenCoefLoc(   k_d,l_d,        m_n        )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    ! store submatrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + exciTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_ET
  
  SUBROUTINE T2EXEC_SS
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,SourScaCoef
    USE T2INTG,ONLY: SourScaIntgPG
    !USE T2CALV,ONLY: SourScaCoef

    INTEGER(ikind)::&
         i_n,j_n,k_n,&
         i_d,&
         i_m
    
    REAL(   rkind)::&
         sourScaVecElem( 1:NNMAX),&
         sourScaCoefElem(1:NNMAX)

    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = NEgraphElem(i_n,1)
       sourScaCoefElem(   i_n            ) &
            = SourScaCoef(    I_v,J_v,i_m) &
            * JacDetLocCrd
    ENDDO

    ! main loop
    
    DO i_n = 1, NNMAX
       sourScaVecElem(i_n) = 0.D0
       DO j_n = 1, NNMAX
          sourScaVecElem(            i_n) &
               = sourScaVecElem(     i_n) &
               + SourScaIntgPG(  j_n,i_n) &
               * sourScaCoefElem(j_n    )
       ENDDO
    ENDDO
    
    ! store subvector

    DO i_n = 1, NNMAX
       BvecElem(             i_n,I_v) &
            = BvecElem(      i_n,I_v) & 
            + sourScaVecElem(i_n    )
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_SS
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR STORE SUBMATRIX 
  !C FOR BI-LINEAR RECTANGULAR ELEMENT
  !C
  !C
  !C
  !C-------------------------------------------------------------------
  
  SUBROUTINE T2EXEC_STORE
    
    USE T2COMM,ONLY:&
         i0nmax,i0vmax,i0bmax,&
         i1nidr,i1nidc,i2hbc,i2enr0,&
         d4smat,d3amat,d2svec,d2bvec
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,&
         i0aidi,&
         i0bidi,i0bidj,i0bidk,&
         i0hidi,i0hidj,i0hidk,i0hidl,&
         i0xidi,i0xidj


    !C
    !C
    !C MATRIX
    !C
    !C
    
    !C
    !C 1Dx1D
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax         
       i0bidi  = i2enr0(i0nidi,4)
       i0bidj  = i2enr0(i0nidj,4)
       DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
          i0bidk = i1nidc(i0aidi)
          IF(i0bidk.EQ.i0bidj)THEN
             DO i0vidj = 1,3
             DO i0vidi = 1,3
                d3amat(                     i0vidi,i0vidj,i0aidi) &
                     = d3amat(              i0vidi,i0vidj,i0aidi) &
                     + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
             ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C 1Dx2D
    !C
   
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
              
       i0bidi = i2enr0(i0nidi,3)
       i0xidj = i2enr0(i0nidj,2)
       
       IF(    i0xidj.LE.i0bmax)THEN
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 1, 3
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF(i0xidj.GT.i0bmax)THEN
          i0xidj = i0xidj - i0bmax
          i0hidi = i2hbc(1,i0xidj)
          i0hidj = i2hbc(2,i0xidj)
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidi).OR.(i0bidk.EQ.i0hidj))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 1, 3
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    
    !C
    !C 2Dx1D
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       i0xidi = i2enr0(i0nidi,2)
       i0bidj = i2enr0(i0nidj,4)
       IF(    i0xidi.LE.i0bmax)THEN
          DO i0aidi = i1nidr(i0xidi), i1nidr(i0xidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF(i0bidk.EQ.i0bidj)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       ELSEIF(i0xidi.GT.i0bmax)THEN
       
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          
          DO i0aidi = i1nidr(i0hidi), i1nidr(i0hidi+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0bidj)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i0aidi = i1nidr(i0hidj), i1nidr(i0hidj+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0bidj)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       
       ENDIF
    ENDDO
    ENDDO
      
    !C 
    !C 2Dx2D
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       
       i0xidi = i2enr0(i0nidi,2)
       i0xidj = i2enr0(i0nidj,2)
       
       IF((i0xidi.LE.i0bmax).AND.(i0xidj.LE.i0bmax))THEN
          
          DO i0aidi = i1nidr(i0xidi), i1nidr(i0xidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
       ELSEIF((i0xidi.LE.i0bmax).AND.(i0xidj.GT.i0bmax))THEN
          
          i0xidj = i0xidj - i0bmax
          i0hidi = i2hbc(1,i0xidj)
          i0hidj = i2hbc(2,i0xidj)
          
          DO i0aidi = i1nidr(i0xidi), i1nidr(i0xidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidi).OR.(i0bidk.EQ.i0hidj))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i0xidi.GT.i0bmax).AND.(i0xidj.LE.i0bmax))THEN
          
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          
          DO i0aidi = i1nidr(i0hidi), i1nidr(i0hidi+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
          DO i0aidi = i1nidr(i0hidj), i1nidr(i0hidj+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i0xidi.GT.i0bmax).AND.(i0xidj.GT.i0bmax))THEN
          
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          
          i0xidj = i0xidj - i0bmax
          i0hidk = i2hbc(1,i0xidj)
          i0hidl = i2hbc(2,i0xidj)
          
          DO i0aidi = i1nidr(i0hidi), i1nidr(i0hidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidk).OR.(i0bidk.EQ.i0hidl))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i0aidi = i1nidr(i0hidj), i1nidr(i0hidj+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidk).OR.(i0bidk.EQ.i0hidl))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ENDIF
       
    ENDDO
    ENDDO
    
    !C
    !C
    !C RIGHT HANDSIDE VECTOR
    !C
    !C
    
    !C
    !C 1D
    !C
    DO i0nidi = 1, i0nmax
       i0bidi = i2enr0(i0nidi,4)
       DO i0vidi = 1, 3
          d2bvec(              i0vidi,i0bidi) &
               = d2bvec(       i0vidi,i0bidi) &
               + d2svec(i0nidi,i0vidi       )
       ENDDO
    ENDDO
    
    !C
    !C 2D
    !C
    
    DO i0nidi = 1, i0nmax
       i0xidi = i2enr0(i0nidi,2)
       IF(    i0xidi.LE.i0bmax)THEN
          DO i0vidi = 4, i0vmax
             d2bvec(              i0vidi,i0xidi) &
                  = d2bvec(       i0vidi,i0xidi) &
                  + d2svec(i0nidi,i0vidi       )
          ENDDO
       ELSEIF(i0xidi.GT.i0bmax)THEN
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          DO i0vidi = 4, i0vmax
             d2bvec(              i0vidi,i0hidi) &
                  = d2bvec(       i0vidi,i0hidi) &
                  + d2svec(i0nidi,i0vidi       ) &
                  * 0.50D0
             d2bvec(              i0vidi,i0hidj) &
                  = d2bvec(       i0vidi,i0hidj) &
                  + d2svec(i0nidi,i0vidi       ) &
                  * 0.50D0
          ENDDO
       ENDIF
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_STORE
 
  !C-------------------------------------------------------------------
  !C 
  !C BOUNDARY CONDITIONS
  !C
  !C                     2014-03-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2EXEC_BCOND
    
    USE T2COMM, ONLY:&
         i0solv,i0lmax,i0bmax,i0vmax,&
         i1pdn2,i1nidr,i1nidc,&
         d3amat,d2bvec,d2xvec
    
    INTEGER::&
         i0aidi,&
         i0bidi,i0bidj,&
         i0vidi,i0vidj,&
         i0bsta,i0bend,i0pdn2
    
    i0pdn2 = i1pdn2(i0lmax)
    i0bsta = i0bmax - i0pdn2 + 1
    i0bend = i0bmax

    !C
    !C
    !C SET FIXED VARIALES 
    !C 
    !C 
    
    SELECT CASE(i0solv)
       
       !C i0solv = 1
       !C SOLVE ONLY ELECTRON DENSITY 
       !C
       
    CASE(1)
       
       !C
       !C
       !C SET FIXED VALUES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,7:)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,7:)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,7:)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,7:)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
       ENDDO
       
       !C
       !C i0solv = 2
       !C SOLVE ONLY ELECTRON AND
       !C            ION DENSITIES AND MOMENTUMS
       !C
       
    CASE(2)
       
       !C
       !C
       !C SET FIXED VALUES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,11:15,21:25)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,11:15,21:25)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
          
       ENDDO

       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
       
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:16,20:25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:16,20:)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
          
       ENDDO
       
       !C
       !C i0slov = 3
       !C SOLVE ELECTRON DENSITY, MOMENTUMS,
       !C            AND PRESSURE
       !C
       
    CASE(3)

       !C
       !C
       !C SET FIXED VALUES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,11,16,21)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,11,16,21)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:11,15:16,20:21,25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                END SELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:11,15:16,20:21,25)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
       ENDDO
       
       !C
       !C i0solv = 4: 
       !C SOLVE ELECTRON 
       !C
       
    CASE(4)
       
       !C
       !C
       !C SET FIXED VARIABLES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,16:)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,16:)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C SET DIRICHLET CONDITION 
       !C
    
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:11,15:)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
        
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:11,15:)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
       ENDDO

       !C
       !C i0solv = 5: 
       !C SOLVE ELECTRON AND ION DENSITIES, 
       !C                MOMENTUMS AND HEATFLUXES
       !C

    CASE(5)

       !C
       !C
       !C SET FIXED VARIABLES
       !C
       !C

       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       

       !C
       !C
       !C SET DIRICHLET CONDITION (MAGNETIC AXIS)
       !C
       !C
       
       i0bidi = 1
       
       !C
       !C STIFFNESS MATRIX
       !C
       
       DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
          i0bidj = i1nidc(i0aidi)
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,8:11,13:16,18:21,23:25)
                CYCLE
             CASE DEFAULT
                DO i0vidj = 1, i0vmax
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                ENDDO
             ENDSELECT
          ENDDO
       ENDDO
       
       !C
       !C RHS VECTOR 
       !C
       
       DO i0vidi = 1, i0vmax
          SELECT CASE(i0vidi)
          CASE(1:6,8:11,13:16,18:21,23:25)
             CYCLE
          CASE DEFAULT
             d2bvec(i0vidi,i0bidi) = 0.D0
          ENDSELECT
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION (FIRST WALL)
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                !CASE(1:6,10:11,15:16,20:21,25)
                CASE(1:5,7:8,10,12:13,15,17:18,20,22:23,25)
                !CASE(1:6,8,10:11,13,15:16,18,20:21,23,25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             !CASE(1:6,10:11,15:16,20:21,25)
             !CASE(1:7,10:12,15:17,20:22,25)
             CASE(1:5,7:8,10,12:13,15,17:18,20,22:23,25)
             !CASE(1:6,8,10:11,13,15:16,18,20:21,23,25)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
       ENDDO

    CASE(10)

       !C
       !C
       !C SET FIXED VARIABLES
       !C
       !C

       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(6:)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(6:)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(6:)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
        
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(6:)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
       ENDDO
    CASE DEFAULT
       WRITE(6,*)'INCORRECT I0SOLV'
    ENDSELECT
    
    RETURN
    
  END SUBROUTINE T2EXEC_BCOND

  SUBROUTINE T2EXEC_PRIVATE_ALLOCATE
    RETURN
  END SUBROUTINE T2EXEC_PRIVATE_ALLOCATE

  !
  !
  ! MODIFIED 2014-03-09
  !
  !
  SUBROUTINE T2EXEC_SOLVE
    
    USE T2COMM, ONLY:&
         i0vmax,i0bmax,i0xmax,i0amax,i0lmax,i0bvmax,i0avmax,&
         i0dbg,i2vvvt,i2hbc,idebug,&
         i1nidr,i1nidc,i1pdn2,i1rdn2,&
         d3amat,d2bvec,d2xvec,d2xvec_befor,d2xvec_after
    
    USE LIBMPI
    USE COMMPI
    USE LIBMTX
    

    INTEGER(i0ikind)::istart,iend,its
    INTEGER(i0ikind)::itype, m1, m2
    REAL(   i0rkind)::tolerance,d0val
    REAL(   i0rkind),POINTER,SAVE::x(:)
    INTEGER(i0ikind)::&
         i0nr,i0nc,i0pdn2,i0rdn2
    INTEGER(i0ikind)::&
         i0lidi,i0ridi,i0pidi,i0aidi,i0bidi,&
         i0xidi,i0xidc,i0xidd,i0xidu,i0xrd,i0xru,&
         i0vvvt,i0offset,&
         i0br,i0xr,i0ar,i0ac,&
         i0arc,i0arc1,i0arc2,i0arc3,&
         i0acc,i0acc1,i0acc2,i0acc3,i0acl,i0acl1,i0acl2,i0acl3

100 FORMAT(A5,I3,A5,I3,A5,I3,A5,D15.6,A5,D15.6)

    itype = 0
    m1    = 4
    
    IF(nsize.EQ.1)THEN
       m2 = 5
    ELSE
       m2 = 0
    ENDIF
    
    tolerance=1.D-7
    idebug = 0

    i0bvmax = i0bmax*i0vmax
    i0avmax = i0amax*i0vmax*i0vmax
    
       
    ALLOCATE(x(i0bvmax))
    
    CALL MTX_SETUP(i0bvmax,istart,iend,nzmax=i0avmax,idebug=0)
    
    !C 
    !C STORE GLOBAL STIFFNESS MATRIX  
    !C 
    DO i0nr = 1, i0bmax   
       DO i0aidi = i1nidr(i0nr), i1nidr(i0nr+1)-1
          i0nc = i1nidc(i0aidi) 
          DO i0vidj = 1, i0vmax
          DO i0vidi = 1, i0vmax
             i0vvvt = i2vvvt(i0vidi,i0vidj)
             IF(i0vvvt.EQ.1) THEN
                d0val = d3amat(i0vidi,i0vidj,i0aidi)
                i0ar  = i0vmax*(i0nr-1) + i0vidi
                i0ac  = i0vmax*(i0nc-1) + i0vidj
                CALL MTX_SET_MATRIX(i0ar,i0ac,d0val)
                IF(IDEBUG.EQ.1) THEN
                   WRITE(18,'(2I5,I10,2I3,2I7,1PE12.4)') &
                        i0nr,i0nc,i0aidi,i0vidi,i0vidj,i0ar,i0ac,d0val
                END IF
                
             END IF
             
          ENDDO
          ENDDO
       ENDDO
    ENDDO

    !C
    !C ADDITIONAL COMPONENTS
    !C FOR FLUX SURFACE AVERAGING
    !C

    GOTO 2000

    i0offset  = 1

    DO i0lidi = 1, i0lmax
       
       i0pdn2 = i1pdn2(i0lidi)
       i0rdn2 = i1rdn2(i0lidi)
       
       DO i0ridi = 1, i0rdn2
          DO i0pidi = 1, i0pdn2
                
             i0arc  = i0vmax*i0offset
             i0arc1 = i0arc + 1
             i0arc2 = i0arc + 2
             i0arc3 = i0arc + 3
                
             i0acc  = i0arc
             i0acc1 = i0arc1
             i0acc2 = i0arc2
             i0acc3 = i0arc3
                
             i0acl  = i0vmax*(i0offset-1)
             i0acl1 = i0acl + 1
             i0acl2 = i0acl + 2
             i0acl3 = i0acl + 3
                
             IF((i0pidi.GE.1).AND.(i0pidi.LT.i0pdn2))THEN
                CALL MTX_SET_MATRIX(i0arc1,i0acc1,-1.D0)
                CALL MTX_SET_MATRIX(i0arc2,i0acc2,-1.D0)
                CALL MTX_SET_MATRIX(i0arc3,i0acc3,-1.D0)
             ENDIF
             
             IF((i0pidi.GT.1).AND.(i0pidi.LE.i0pdn2))THEN!
                CALL MTX_SET_MATRIX(i0arc1,i0acl1, 1.D0)
                CALL MTX_SET_MATRIX(i0arc2,i0acl2, 1.D0)
                CALL MTX_SET_MATRIX(i0arc3,i0acl3, 1.D0)
            ENDIF
             
             i0offset = i0offset + 1
             
          ENDDO
       ENDDO
    ENDDO
    
2000 CONTINUE

    !C
    !C SET GLOBAL RIGHT HAND SIDE VECTOR
    !C
    
    DO i0bidi = 1, i0bmax
       DO i0vidi = 1, i0vmax
          i0br  = i0vmax*(i0bidi-1) + i0vidi
          d0val = d2bvec(i0vidi,i0bidi)
          CALL MTX_SET_SOURCE(i0br,d0val)
       ENDDO
    ENDDO
    
    !C
    !C SET INITIAL VALUE IN X
    !C
    
    DO i0xidi = 1, i0bmax
       DO i0vidi = 1, i0vmax
          i0xr  = i0vmax*(i0xidi-1) + i0vidi
          d0val = d2xvec_befor(i0vidi,i0xidi)
          CALL MTX_SET_VECTOR(i0xr,d0val)
       ENDDO
    ENDDO
    
    CALL MTX_SOLVE(itype,tolerance,its,methodKSP=m1,methodPC=m2)
    
    CALL MTX_GATHER_VECTOR(x)
    
    DO i0xidi = 1, i0xmax
       IF(    i0xidi.LE.i0bmax)THEN
          DO i0vidi = 1, i0vmax
             i0xr = i0vmax*(i0xidi - 1) + i0vidi
             d2xvec_after(i0vidi,i0xidi) = x(i0xr)
          ENDDO
       ELSEIF(i0xidi.GT.i0bmax)THEN
          i0xidc = i0xidi - i0bmax
          i0xidd = i2hbc(1,i0xidc)
          i0xidu = i2hbc(2,i0xidc)
          DO i0vidi = 1, i0vmax 
             i0xrd = i0vmax*(i0xidd-1) + i0vidi
             i0xru = i0vmax*(i0xidu-1) + i0vidi
             d2xvec_after(i0vidi,i0xidi) = 0.5D0*(x(i0xrd)+x(i0xru))
          ENDDO
       ENDIF
    ENDDO
    
    CALL MTX_CLEANUP
    
    DEALLOCATE(x)
    
    RETURN
    
  END SUBROUTINE T2EXEC_SOLVE
END MODULE T2EXEC
