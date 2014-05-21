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
!C        RowCRS, ColCRS, ENGraph, HangedNodeTable    [from T2NGRA]
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
       ikind,rkind,ikind,rkind
  
  IMPLICIT NONE
  
  PRIVATE


  INTEGER(ikind),SAVE::&
       I_v,J_v,&! variable       index
       I_k,J_k,&! known variable index
       I_e      ! element        index
  
  INTEGER(ikind),SAVE,ALLOCATABLE::& 
       ENGraphElem(:,:)   ! node-element graph table in en element       
  REAL(   rkind),SAVE::&
       JacDetLocCrd       ! jacobian of local coordinates
  REAL(   rkind),SAVE,ALLOCATABLE::&
       JacInvLocCrd(:,:),&! inverse jacobi matrix of local coordinates
       ValKnownElem(:,:),&! known variables at nodes in an element 
       BvecElem(:,:),&    ! element right hand side vector  (Ax=b)
       AmatElem(:,:,:,:),&! element stiffness matrix        (Ax=b)
       Amat(:,:,:),&
       Bvec(:,:),&
       Xvec(:)

  INTEGER(ikind),SAVE,ALLOCATABLE::&
       RowCRS(:),ColCRS(:),HangedNodeTable(:,:),ENGraph(:,:,:)
 
  REAL(   rkind),SAVE,ALLOCATABLE::&
       MassScaCoef(:,:,:    ),AdveVecCoef(:,:,:,:  ),&
       AdveTenCoef(:,:,:,:,:),DiffTenCoef(:,:,:,:,:),&
       GradVecCoef(:,:,:,:  ),GradTenCoef(:,:,:,:,:),&
       ExciScaCoef(:,:,:    ),ExciVecCoef(:,:,:,:  ),&
       ExciTenCoef(:,:,:,:,:),SourScaCoef(:,:,:    )
  
  LOGICAL,ALLOCATABLE::&
       VariableGraphMat(:,:),&
       VariableGraphVec(:,:),&
       EquationTable(:,:),&
       DirichletAxisTable(:,:),&
       DirichletWallTable(:,:)
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
    
    CALL T2EXEC_ADHOC!this will be removed soon

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

          IF(i0msvt.EQ.1)       CALL T2EXEC_MS_SUBMATRIX
          
          IF(i0avvt.EQ.1)       CALL T2EXEC_AV_SUBMATRIX
          
          IF(i0atvt.EQ.1)THEN
             DO i_k = 1, NKMAX
                i0atwt = i3atwt(I_k,I_v,J_v)
                IF(i0atwt.EQ.1) CALL T2EXEC_AT_SUBMATRIX
             ENDDO
          ENDIF
          
          IF(i0dtvt.EQ.1)       CALL T2EXEC_DT_SUBMATRIX
          
          IF(i0gvvt.EQ.1)       CALL T2EXEC_GV_SUBMATRIX
          
          IF(i0gtvt.EQ.1)THEN
             DO I_k = 1, NKMAX
                i0gtwt = i3gtwt(I_k,I_v,J_v)
                IF(i0gtwt.EQ.1) CALL T2EXEC_GT_SUBMATRIX
             ENDDO
          ENDIF
          
          IF(i0esvt.EQ.1)       CALL T2EXEC_ES_SUBMATRIX
          
          IF(i0evvt.EQ.1)THEN
             DO I_k = 1, NKMAX
                i0evwt = i3evwt(I_k,I_v,J_v)
                IF(i0evwt.EQ.1) CALL T2EXEC_EV_SUBMATRIX
             ENDDO
          ENDIF
          
          IF(i0etvt.EQ.1)THEN          
             DO J_k = 1, NKMAX
             DO I_k = 1, NKMAX
                i0etwt = i4etwt(I_k,J_k,I_v,J_v)
                IF(i0etwt.EQ.1) CALL T2EXEC_ET_SUBMATRIX
             ENDDO
             ENDDO
          ENDIF
          
          IF(i0ssvt.EQ.1)       CALL T2EXEC_SS_SUBMATRIX
          
       ENDDO
       ENDDO

       CALL T2EXEC_STORE_SUBMATRIX
       
    ENDDO
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)')&
         '-- T2EXEC: matrix construction was completed: cpu=', &
         e0time_1-e0time_0,' [s]'
    
    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_SET_EQUATIONS_TO_BE_SOLVED
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)')&
         '-- T2EXEC: equation setting was completed:    cpu=', &
         e0time_1-e0time_0,' [s]'
    
    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_SET_DIRICHLET_BOUNDARY
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)')&
         '-- T2EXEC: boundary setting was completed:    cpu=', &
         e0time_1-e0time_0,' [s]'
    
    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_SOLVE_MATRIX_EQUATION
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_SOLVE completed:  cpu=', &
         e0time_1-e0time_0,' [s]'
    
    RETURN
    
  END SUBROUTINE T2_EXEC
  
  !------------------------------------------------------------------ 
  !
  !  T2EXEC_SETUP_ELEMENT_VARIABLES (PRIVATE)
  !   　
  !     * SET LOCAL NODE-ELEMENT GRAPH
  !     * CALCULATE JACOBIAN OF LOCAL COORDINATES
  !     * SET KNOWN VARIABLES FOR DIFFERENTIAL
  !     * SET STABILIZATION FACTORS FOR SUPG (killed)
  ! 
  !                 2014-05-21 H.SETO 
  !
  !------------------------------------------------------------------  
  SUBROUTINE T2EXEC_SETUP_ELEMENT_VARIABLES
    
    USE T2COMM, ONLY:&
         NNMAX,NDMAX,NVMAX,NSMAX,&
         d2ws, d2xvec,d2wrks,&
         i3enr,d2mfc1
    !USE T2CALV, ONLY: ValKnown
    !USE T2NGRA, ONLY: ENGraph
    INTEGER(ikind)::i_r,i_m,i_s,i_n,j_n,k_n
    REAL(   rkind)::gridSizeRad,gridSizePol
    
    CALL T2EXEC_PRIVATE_ALLOCATE
    
    !
    ! SET LOCAL NODE-ELEMENT GRAPH
    !
    !   ENGraphElem(1:NNMAX,1) : FOR COEFFICIENT CALCULATION
    !   ENGraphElem(1:NNMAX,2) : FOR 2D-2D 
    !   ENGraphElem(1:NNMAX,3) : FOR 1D-2D,1D-1D 
    !   ENGraphElem(1:NNMAX,4) : FOR 2D-1D
    !
    
    DO i_r = 1, 4
       DO i_n = 1, NNMAX
          ENGraphElem(i_n,i_r)=ENGraph(i_n,i_r,i_e)
       ENDDO
    ENDDO
    
    
    !
    ! CALCULATE JACOBIAN OF LOCAL COORDINATE 
    ! 
    ! JacDetLocCrd: JACOBIAN OF LOCAL COORDINATES
    ! JacInvLocCrd: INVERSE JACOBI MATRIX OF LOCAL COORDINATES
    ! 
    
    i_n = ENGraphElem(1,1)
    j_n = ENGraphElem(2,1)
    k_n = ENGraphElem(4,1)

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
       i_m = ENGraphElem(i_n,1)
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
  !    MASS SCALAR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_MS_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,Dt
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
       i_m = ENGraphElem( i_n,1)
       massScaCoefElem(   i_n            ) &
            = MassScaCoef(    i_m,I_v,J_v) &
            * JacDetLocCrd/Dt
    ENDDO
    
    SELECT CASE(J_v)
       
    CASE (1:3)   ! for FSA variables 
       DO i_n = 1, NNMAX
          i_b = ENGraphElem(i_n,4)
          valPriorElem(i_n) = Xvec(i_b,J_v)
       ENDDO
    CASE DEFAULT ! for 2D dependent variables
       DO i_n = 1, NNMAX
          i_x = ENGraphElem(i_n,2)
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
    
  END SUBROUTINE T2EXEC_MS_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !    ADVECTION VECTOR SUBMATRCES (PRIVATE)
  !
  !                 2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_AV_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: AdveVecIntgPG
    !USE T2CALV,ONLY: AdveVecCoef

    INTEGER(ikind)::&
         i_d,j_d,&
         i_n,j_n,k_n,&
         i_m
    
    REAL(   rkind)::&
         adveVecMatElem( 1:NNMAX,1:NNMAX),&
         adveVecCoefElem(1:NDMAX,1:NNMAX),&
         adveVecCoefLoc( 1:NDMAX,1:NNMAX)
    
    ! INTITIALIZATION
        
    DO i_n = 1, NNMAX
       i_m = ENGraphElem(i_n,1)
       DO i_d = 1, NDMAX
          adveVecCoefElem(   i_d,i_n            ) &
               = AdveVecCoef(i_d,    I_v,J_v,i_m) &
               * JacDetLocCrd
          
       ENDDO
    ENDDO
    
    ! MAIN LOOP
    
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

    ! STORE SUBMATRIX
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + adveVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AV_SUBMATRIX
  
  !-------------------------------------------------------------------
  ! 
  !    ADVECTION TENSOR SUBMATRCES (PRIVATE)
  ! 
  !                 2014-05-21 H.SETO
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_AT_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
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
       adveTenKnownElem(i_n) = ValKnownElem(i_n,I_k)
    ENDDO
    
    DO i_n = 1, NNMAX
       i_m = ENGraphElem(i_n,1)
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
          DO j_d = 1, NDMAX        
          DO i_d = 1, NDMAX
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
    
  END SUBROUTINE T2EXEC_AT_SUBMATRIX

  !-------------------------------------------------------------------
  !
  !  DIFFUSION TENSOR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-20 H.Seto
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_DT_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: DiffTenIntgPG
    !USE T2COEF,ONLY: DiffTenCoef
    INTEGER(ikind)::&
         i_n,j_n,k_n,&!l_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         diffTenMatElem( 1:NNMAX,1:NNMAX),&
         diffTenCoefElem(1:NDMAX,1:NDMAX,1:NNMAX), &
         diffTenCoefLoc( 1:NDMAX,1:NDMAX,1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = ENGraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          diffTenCoefElem(   i_d,j_d,i_n            ) &
               = diffTenCoef(i_d,j_d,    I_v,J_v,i_m) &
               * JacDetLocCrd
       ENDDO
       ENDDO
    ENDDO

    ! MAIN LOOP

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

    ! store submatrix
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
            + diffTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_DT_SUBMATRIX

  !-------------------------------------------------------------------
  !
  !   GRADIENT VECTOR SUBMATRICES (PRIVATE)
  !
  !                 2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GV_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: GradVecIntgPG
    !USE T2CALV,ONLY: GradVecCoef
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,&!l_n,&
         i_d,j_d,&
         i_m
    
    REAL(   rkind)::&
         gradVecMatElem( 1:NNMAX,1:NNMAX),&
         gradVecCoefElem(1:NDMAX,1:NNMAX),&
         gradVecCoefLoc( 1:NDMAX,1:NNMAX)

    ! initialization
          
    DO i_n = 1, NNMAX
       i_m = ENGraphElem(i_n,1)
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
    
  END SUBROUTINE T2EXEC_GV_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !   GRADIENT TENSOR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GT_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
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
    
    ! initialization
       
    DO i_n = 1, NNMAX
       i_m = ENGraphElem(i_n,1)
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
    
    ! main loop
    
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
    
  END SUBROUTINE T2EXEC_GT_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !  EXCITATION SCALAR SUBMATRCES (PRIVATE)
  !
  !                 2014-05-21 H.SETO
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_ES_SUBMATRIX

    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
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
       i_m = ENGraphElem(i_n,1)
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
    
  END SUBROUTINE T2EXEC_ES_SUBMATRIX

  !------------------------------------------------------------------
  !
  !    EXCITATION SUBMATRCES (PRIVATE)
  !
  !                 2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_EV_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: ExciVecIntgPG
    !USE T2CALV,ONLY: ExciVecCoef
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&!m_n,&
         i_d,j_d,&
         i_m

    REAL(   rkind)::&
         exciVecMatElem(  1:NNMAX,1:NNMAX),&
         exciVecCoefElem( 1:NDMAX,1:NNMAX),&
         exciVecCoefLoc(  1:NDMAX,1:NNMAX),&    
         exciVecKnownElem(1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = ENGraphElem(i_n,1)
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
             exciVecMatElem(                     i_n,j_n) &
                  = exciVecMatElem(              i_n,j_n) &
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
    
  END SUBROUTINE T2EXEC_EV_SUBMATRIX

  !-------------------------------------------------------------------
  !
  !     EXCITATION TENESOR SUBMATRCES (PRIVATE)
  !
  !                 2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_ET_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
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
       i_m = ENGraphElem(i_n,1)
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
    
  END SUBROUTINE T2EXEC_ET_SUBMATRIX
  
  SUBROUTINE T2EXEC_SS_SUBMATRIX
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
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
       i_m = ENGraphElem(i_n,1)
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
    
  END SUBROUTINE T2EXEC_SS_SUBMATRIX
  
  !-------------------------------------------------------------------
  !
  ! SUBROUTINE FOR STORE SUBMATRIX 
  ! FOR BI-LINEAR RECTANGULAR ELEMENT
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_STORE_SUBMATRIX
    
    USE T2COMM,ONLY:NNMAX,NVMAX,NBMAX
    !USE T2NGRA,ONLY: RowCRS,ColCRS,HangedNodeTable
    
    INTEGER(ikind)::&
         i_n,j_n,&
         i_a,&
         i_b,j_b,k_b,&
         i_h,j_h,k_h,l_h,&
         i_x,j_x
    REAL(rkind)::&
         amatElemTF(1:NVMAX,1:NVMAX,1:NNMAX,1:NNMAX),&
         bvecElemTF(1:NVMAX,1:NNMAX)
    
    LOGICAL::isThisValSolved
    ! transform and clear AmatElem 
    DO J_v = 1, NVMAX
    DO I_v = 1, NVMAX
       isValueStored = VariableGraphMat(I_v,J_v)
       IF(isValueStored)THEN
          DO j_n = 1, NNMAX
          DO i_n = 1, NNMAX
             amatElemTF(I_v,J_v,i_n,j_n) = AmatElem(i_n,j_n,I_v,J_v)
             AmatElem(  i_n,j_n,I_v,J_v) = 0.D0
          ENDDO
          ENDDO
       ELSE
          amatElemTF(I_v,J_v,1:NNMAX,1:NNMAX) = 0.D0
          AmatElem(  1:NNMAX,1:NNMAX,I_v,J_v) = 0.D0
       ENDIF
    ENDDO
    ENDDO
    
    ! transform and clear BvecElem 
    
    DO I_v = 1, NVMAX
       isValueStored = VariableGraphVec(I_v)
       IF(isValueStored)THEN! this will be modified soon
          DO i_n = 1, NNMAX
             bvecElemTF(I_v,i_n) = BvecElem(i_n,I_v)
             BvecElem(  i_n,I_v) = 0.D0
          ENDDO
       ELSE
          bvecElemTF(I_v,1:NNMAX) = 0.D0
          BvecElem(  1:NNMAX,I_v) = 0.D0
       ENDIF
    ENDDO
    
    !
    ! for stiffness matrix
    !
    
    ! I_v: 1D value (FSA)
    ! J_v: 1D value (FSA)
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX         
       i_b  = ENGraphElem(i_n,4)
       j_b  = ENGraphElem(j_n,4)
       DO i_a = RowCRS(i_b), RowCRS(i_b+1)-1
          k_b = ColCRS(i_a)
          IF(k_b.EQ.j_b)THEN
             DO J_v = 1,3
             DO I_v = 1,3
                Amat(             I_v,J_v,        i_a) &
                     = Amat(      I_v,J_v,        i_a) &
                     + amatElemTF(I_v,J_v,i_n,j_n    )
             ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    
    
    ! I_v: 1D value (FSA)
    ! J_v: 2D value
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       
       i_b = ENGraphElem(i_n,3)
       j_x = ENGraphElem(j_n,2)
       
       IF(    j_x.LE.NBMAX)THEN
          DO i_a = RowCRS(i_b), RowCRS(i_b+1)-1
             k_b = ColCRS(i_a)
             IF(k_b.EQ.j_x)THEN
                DO J_v = 4, NVMAX
                DO I_v = 1, 3
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF(j_x.GT.NBMAX)THEN
          j_x = j_x - NBMAX
          i_h = HangedNodeTable(1,j_x)
          j_h = HangedNodeTable(2,j_x)
          DO i_a = RowCRS(i_b), RowCRS(i_b+1)-1
             k_b = ColCRS(i_a)
             IF((k_b.EQ.i_h).OR.(k_b.EQ.j_h))THEN
                DO J_v = 4, NVMAX
                DO I_v = 1, 3
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    

    ! I_v: 2D value 
    ! J_v: 1D value (FSA) 
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       i_x = ENGraphElem(i_n,2)
       j_b = ENGraphElem(j_n,4)
       IF(    i_x.LE.NBMAX)THEN
          DO i_a = RowCRS(i_x), RowCRS(i_x+1)-1
             k_b = ColCRS(i_a)
             IF(k_b.EQ.j_b)THEN
                DO J_v = 1, 3
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       ELSEIF(i_x.GT.NBMAX)THEN
       
          i_x = i_x - NBMAX
          i_h = HangedNodeTable(1,i_x)
          j_h = HangedNodeTable(2,i_x)
          
          DO i_a = RowCRS(i_h), RowCRS(i_h+1)-1
             k_b = ColCRS(i_a) 
             IF(k_b.EQ.j_b)THEN
                DO J_v = 1, 3
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i_a = RowCRS(j_h), RowCRS(j_h+1)-1
             k_b = ColCRS(i_a) 
             IF(k_b.EQ.j_b)THEN
                DO J_v = 1, 3
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO

       ENDIF
    ENDDO
    ENDDO
      

    ! I_v: 2D value
    ! J_v: 2D value
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       
       i_x = ENGraphElem(i_n,2)
       j_x = ENGraphElem(j_n,2)
       
       IF((i_x.LE.NBMAX).AND.(j_x.LE.NBMAX))THEN
          
          DO i_a = RowCRS(i_x), RowCRS(i_x+1)-1
             k_b = ColCRS(i_a)
             IF(k_b.EQ.j_x)THEN
                DO J_v = 4, NVMAX
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
       ELSEIF((i_x.LE.NBMAX).AND.(j_x.GT.NBMAX))THEN
          
          j_x = j_x - NBMAX
          i_h = HangedNodeTable(1,j_x)
          j_h = HangedNodeTable(2,j_x)
          
          DO i_a = RowCRS(i_x), RowCRS(i_x+1)-1
             k_b = ColCRS(i_a)
             IF((k_b.EQ.i_h).OR.(k_b.EQ.j_h))THEN
                DO J_v = 4, NVMAX
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i_x.GT.NBMAX).AND.(j_x.LE.NBMAX))THEN
          
          i_x = i_x - NBMAX
          i_h = HangedNodeTable(1,i_x)
          j_h = HangedNodeTable(2,i_x)
          
          DO i_a = RowCRS(i_h), RowCRS(i_h+1)-1
             k_b = ColCRS(i_a) 
             IF(k_b.EQ.j_x)THEN
                DO J_v = 4, NVMAX
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
          DO i_a = RowCRS(j_h), RowCRS(j_h+1)-1
             k_b = ColCRS(i_a) 
             IF(k_b.EQ.j_x)THEN
                DO J_v = 4, NVMAX
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i_x.GT.NBMAX).AND.(j_x.GT.NBMAX))THEN
          
          i_x = i_x - NBMAX
          i_h = HangedNodeTable(1,i_x)
          j_h = HangedNodeTable(2,i_x)
          
          j_x = j_x - NBMAX
          k_h = HangedNodeTable(1,j_x)
          l_h = HangedNodeTable(2,j_x)
          
          DO i_a = RowCRS(i_h), RowCRS(i_h+1)-1
             k_b = ColCRS(i_a)
             IF((k_b.EQ.k_h).OR.(k_b.EQ.l_h))THEN
                DO J_v = 4, NVMAX
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i_a = RowCRS(j_h), RowCRS(j_h+1)-1
             k_b = ColCRS(i_a)
             IF((k_b.EQ.k_h).OR.(k_b.EQ.l_h))THEN
                DO J_v = 4, NVMAX
                DO I_v = 4, NVMAX
                   Amat(             I_v,J_v,        i_a) &
                        = Amat(      I_v,J_v,        i_a) &
                        + amatElemTF(I_v,J_v,i_n,j_n    ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       
       ENDIF
       
    ENDDO
    ENDDO
    
    ! RIGHT HANDSIDE VECTOR
    
    ! I_v: 1D value (FSA)
    
    DO i_n = 1, NNMAX
       i_b = ENGraphElem(i_n,4)
       DO I_v = 1, 3
          Bvec(             I_v,    i_b) &
               = Bvec(      I_v,    i_b) &
               + bvecElemTF(I_v,i_n    )
       ENDDO
    ENDDO
    
    ! I_v 2D value
    
    DO i_n = 1, NNMAX
       i_x = ENGraphElem(i_n,2)
       IF(    i_x.LE.NBMAX)THEN
          DO I_v = 4, NVMAX
             Bvec(              I_v,    i_x) &
                  = Bvec(       I_v,    i_x) &
                  + bvecElemTF( I_v,i_n    )
          ENDDO
       ELSEIF(i_x.GT.NBMAX)THEN
          i_x = i_x - NBMAX
          i_h = HangedNodeTable(1,i_x)
          j_h = HangedNodeTable(2,i_x)
          DO I_v = 4, NVMAX
             Bvec(             I_v,    i_h) &
                  = Bvec(      I_v,    i_h) &
                  + bvecElemTF(I_v,i_n    ) &
                  * 0.50D0
             Bvec(             I_v,    j_h) &
                  = Bvec(      I_v,    j_h) &
                  + bvecElemTF(I_v,i_n   ) &
                  * 0.50D0
          ENDDO
       ENDIF
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_STORE_SUBMATRIX
 
  !-------------------------------------------------------------------
  ! 
  !  BOUNDARY CONDITIONS
  !
  !                     2014-05-21 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_BCOND
    
    USE T2COMM, ONLY:NLMAX,NBMAX,NVMAX,i1pdn2
    !USE T2NGRA, ONLY:RowCRS,ColCRS
    
    INTEGER(ikind)::&
         i_a,&
         i_b,j_b,&
         I_v,J_v,&
         i0bsta,i0bend,i0pdn2
    
    i0pdn2 = i1pdn2(NLMAX)
    i0bsta = NBMAX - i0pdn2 + 1
    i0bend = NBMAX

    
    !
    ! SET UNSOLVED VARIALES 
    !
    
    DO i_b= 1, NBMAX
       
       ! STIFFNESS MATRIX
       
       DO i_a = RowCRS(i_b), RowCRS(i_b+1)-1
          j_b = ColCRS(i_a)
          DO I_v = 1, NVMAX
             
             IF(.TRUE.) CYCLE  ! solve    equation I_v 
             
             DO J_v = 1, NVMAX ! unsolve  equation I_v 
                IF((i_b.EQ.j_b).AND.(I_v.EQ.J_v))THEN
                   Amat(I_v,J_v,i_a) = 1.D0
                ELSE
                   Amat(I_v,J_v,i_a) = 0.D0
                ENDIF
             ENDDO
             
          ENDDO
       ENDDO
       
       ! RHS VECTOR 
       
       DO I_v = 1, NVMAX
          IF(.TRUE.) CYCLE             ! solve   equation I_v
          Bvec(I_v,i_b) = Xvec(I_v,i_b)! unsolve equation I_v
       ENDDO
       
    ENDDO
    
    !
    ! set dirichlet condition on magnetic axis
    !
    
    !CALL T2EXEC_DIRICHLET(start,end,flag)
    
    i_b = 1
    
    DO i_a = RowCRS(i_b), RowCRS(i_b+1)-1
       j_b = ColCRS(i_a)
       DO I_v = 1, NVMAX
          
          IF(.TRUE.)THEN ! modify stiffness matrix
             DO J_v = 1, NVMAX
                IF((i_b.EQ.j_b).AND.(I_v.EQ.J_v))THEN
                   Amat(I_v,J_v,i_a) = 1.D0
                ELSE
                   Amat(I_v,J_v,i_a) = 0.D0
                ENDIF
             ENDDO
          ENDIF
          
       ENDDO
    ENDDO
    
    DO I_v = 1, NVMAX
       IF(.TRUE.)THEN ! modify RHS vector
          Bvec(I_v,i_b) = Xvec(I_v,i_b)
       ENDIF
    ENDDO
    
    !
    ! set dirichlet boundary condition on first wall
    !
    
    DO i_b = i0bsta, i0bend
       
       DO i_a = RowCRS(i_b), RowCRS(i_b+1)-1
          j_b = ColCRS(i_a)
          DO I_v = 1, NVMAX
             
             IF(.TRUE.)THEN ! modify
             DO J_v = 1, i0vmax
                   IF((i_b.EQ.j_b).AND.(I_v.EQ.J_v))THEN
                      d3amat(I_v,J_v,i_a) = 1.D0
                   ELSE
                      d3amat(I_v,J_v,i_a) = 0.D0
                   ENDIF
                ENDDO
             ENDSELECT
          ENDDO
       ENDDO
       
       ! RHS VECTOR 
      
       DO I_v = 1, i0vmax
          SELECT CASE(I_v)
          CASE(1:5,7:)
             CYCLE
          CASE DEFAULT
             d2bvec(I_v,i_b) = d2xvec(I_v,i_b)
          END SELECT
       ENDDO

    ENDDO
       
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
         i0dbg,i2vvvt,idebug,i1pdn2,i1rdn2,&
         d3amat,d2bvec,d2xvec,d2xvec_befor,d2xvec_after
    
    USE LIBMPI
    USE COMMPI
    USE LIBMTX
    

    INTEGER(ikind)::istart,iend,its
    INTEGER(ikind)::itype, m1, m2
    REAL(   rkind)::tolerance,d0val
    REAL(   rkind),POINTER,SAVE::x(:)
    INTEGER(ikind)::&
         i0nr,i0nc,i0pdn2,i0rdn2
    INTEGER(ikind)::&
         i0lidi,i0ridi,i0pidi,i_a,i_b,&
         i_x,i0xidc,i0xidd,i0xidu,i0xrd,i0xru,&
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
       DO i_a = RowCRS(i0nr), RowCRS(i0nr+1)-1
          i0nc = ColCRS(i_a) 
          DO J_v = 1, i0vmax
          DO I_v = 1, i0vmax
             i0vvvt = i2vvvt(I_v,J_v)
             IF(i0vvvt.EQ.1) THEN
                d0val = d3amat(I_v,J_v,i_a)
                i0ar  = i0vmax*(i0nr-1) + I_v
                i0ac  = i0vmax*(i0nc-1) + J_v
                CALL MTX_SET_MATRIX(i0ar,i0ac,d0val)
                IF(IDEBUG.EQ.1) THEN
                   WRITE(18,'(2I5,I10,2I3,2I7,1PE12.4)') &
                        i0nr,i0nc,i_a,I_v,J_v,i0ar,i0ac,d0val
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
    
    DO i_b = 1, i0bmax
       DO I_v = 1, i0vmax
          i0br  = i0vmax*(i_b-1) + I_v
          d0val = d2bvec(I_v,i_b)
          CALL MTX_SET_SOURCE(i0br,d0val)
       ENDDO
    ENDDO
    
    !C
    !C SET INITIAL VALUE IN X
    !C
    
    DO i_x = 1, i0bmax
       DO I_v = 1, i0vmax
          i0xr  = i0vmax*(i_x-1) + I_v
          d0val = d2xvec_befor(I_v,i_x)
          CALL MTX_SET_VECTOR(i0xr,d0val)
       ENDDO
    ENDDO
    
    CALL MTX_SOLVE(itype,tolerance,its,methodKSP=m1,methodPC=m2)
    
    CALL MTX_GATHER_VECTOR(x)
    
    DO i_x = 1, i0xmax
       IF(    i_x.LE.i0bmax)THEN
          DO I_v = 1, i0vmax
             i0xr = i0vmax*(i_x - 1) + I_v
             d2xvec_after(I_v,i_x) = x(i0xr)
          ENDDO
       ELSEIF(i_x.GT.i0bmax)THEN
          i0xidc = i_x - i0bmax
          i0xidd = HangedNodeTable(1,i0xidc)
          i0xidu = HangedNodeTable(2,i0xidc)
          DO I_v = 1, i0vmax 
             i0xrd = i0vmax*(i0xidd-1) + I_v
             i0xru = i0vmax*(i0xidu-1) + I_v
             d2xvec_after(I_v,i_x) = 0.5D0*(x(i0xrd)+x(i0xru))
          ENDDO
       ENDIF
    ENDDO
    
    CALL MTX_CLEANUP
    
    DEALLOCATE(x)
    
    RETURN
    
  END SUBROUTINE T2EXEC_SOLVE
END MODULE T2EXEC
