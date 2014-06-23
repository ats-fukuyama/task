!--------------------------------------------------------------------
!  MODULE FOR SOLVING ADVECTION-DIFFUSION PROBLEM BY FEM
! 
!                   LAST UPDATE 2014-06-20 H.SETO
!
!  T2EXEC provides following variables:
!
!
!    and subroutines:
!
!         T2EXEC_EXECUTE
! 
!  T2EXEC requires following variables:
!  
!  [from T2CNST]
!
!        rkind,ikind
!
!  [from T2COMM] 
!
!        NEMAX,NBMAX,NVMAX,NKMAX,NNMAX,NDMAX,Dt  
!
!  [from T2NGRA]
!        NodeRowCRS, NodeColCRS, ElementNodeGraph, HangedNodeTable    
!
!  [from T2PROF] 
!        GlobalCrd
!
!  [from T2INTG]
!        MassScaIntgPG AdveVecIntgPG AdveTenIntgPG DiffTenIntgPG
!        GradVecIntgPG GradTenIntgPG ExciScaIntgPG ExciVecIntgPG
!        ExciTenIntgPG SourScaIntgPG
!
!  [from T2CALV]
!        MassScaCoef AdveVecCoef AdveTenCoef DiffTenCoef GradVecCoef
!        GradTenCoef ExciScaCoef ExciVecCoef ExciTenCoef SourScaCoef
!
!  [from T2STEP]
!        XvecIn,XvecOut
!
!  [from T2VGRA]
!        HaveMassScaCoef HaveAdveVecCoef HaveAdveTenCoef
!        HaveDiffTenCoef HaveGradVecCoef HaveGradTenCoef
!        HaveExciScaCoef HaveExciVecCoef HaveExciTenCoef
!        HaveSourScaCoef
!        HaveAdveTenKval HaveGradTenKval HaveExciVecKval
!        HaveExciTenKval
!        HaveMat         HaveVec
!        LockEqs         LockAxi         LockWal
!
! -------------------------------------------------------------------

MODULE T2EXEC
  
  USE T2CNST, ONLY: ikind,rkind
  
  IMPLICIT NONE
  
  PRIVATE
  
  INTEGER(ikind),SAVE,ALLOCATABLE::& 
       nodeTableA(:),&    ! node-element graph table in element
       nodeTableB(:),&    ! node-element graph table in element
       nodeTableC(:),&    ! node-element graph table in element
       nodeTableD(:)      ! node-element graph table in element
  REAL(   rkind),SAVE::&
       jacDetLocCrd       ! jacobian of local coordinates

  REAL(   rkind),SAVE,ALLOCATABLE::&
       jacInvLocCrd(:,:),  &! inverse jacobi matrix of local coordinates
       knownVarElem(:,:),  &! known variables at nodes in an element 
       bvecElem(:,:),      &! element right hand side vector  (Ax=b)
       amatElem(:,:,:,:),  &! element stiffness matrix        (Ax=b)
       bvecElemTF(:,:),    &
       amatElemTF(:,:,:,:),&
       amat(:,:,:),        &
       bvec(:,:)
  
  ! for debug
  INTEGER(ikind),SAVE::i_elm,nvmx
  
  PUBLIC T2EXEC_EXECUTE,&
       & T2EXEC_DEALLOCATE
  
CONTAINS
  
  !------------------------------------------------------------------
  !
  !  T2EXEC_EXECUTE
  !  FEM SOLVER FOR SIMULTANEOUS ADVECTION-DIFFUSION EQUATIONS
  ! 
  !                2014-06-23 H.SETO
  !
  !------------------------------------------------------------------
  SUBROUTINE T2EXEC_EXECUTE(nvmx_in)
    
    USE T2COMM,ONLY: NNMAX,NEMAX,NKMAX,NVMAX,NBMAX,NAMAX,NDMAX,NVFMX,&
         &           CoordinateSwitch,&
         &           HaveMassScaCoef,HaveAdveVecCoef,HaveAdveTenCoef,&
         &           HaveDiffTenCoef,HaveGradVecCoef,HaveGradTenCoef,&
         &           HaveExciScaCoef,HaveExciVecCoef,HaveExciTenCoef,&
         &           HaveSourScaCoef,&
         &           HaveAdveTenKval,HaveGradTenKval,HaveExciVecKval,&
         &           HaveExciTenKval,&
         &           HaveMat,&
         &           StartEqs,EndEqs,LockEqs,StartAxi,EndAxi,LockAxi,&
         &           StartWal,EndWal,LockWal
    
    INTEGER(ikind),INTENT(IN)::nvmx_in
    INTEGER(ikind)::&
         i_v,j_v,i_k,j_k!,i_e
    REAL(4)::e0time_0,e0time_1

    CALL CPU_TIME(e0time_0)

    nvmx = nvmx_in

    CALL T2EXEC_ALLOCATE
    
    ! INITIALIZATION

    CALL T2EXEC_INITIALIZE_MATRIX

    DO i_elm = 1, NEMAX
       
       amatElem(1:NNMAX,1:NNMAX,1:NVMAX,1:NVMAX) = 0.D0
       bvecElem(1:NNMAX,1:NVMAX)                 = 0.D0
       
       CALL T2EXEC_SETUP_ELEMENT_VARIABLES(i_elm)
       
       DO j_v = 1, NVMAX
       DO i_v = 1, NVMAX
          
          IF(HaveMassScaCoef(i_v,j_v))&
               &     CALL T2EXEC_MS_SUBMATRIX(i_v,j_v        )
          
          IF(HaveAdveVecCoef(i_v,j_v))&
               &     CALL T2EXEC_AV_SUBMATRIX(i_v,j_v        )
           
          IF(HaveAdveTenCoef(i_v,j_v))THEN
             DO i_k = 1, NKMAX
                IF(HaveAdveTenKval(i_k,i_v,j_v))&
                     CALL T2EXEC_AT_SUBMATRIX(i_v,j_v,i_k    )
             ENDDO
          ENDIF
          
          IF(HaveDiffTenCoef(i_v,j_v))&
               &     CALL T2EXEC_DT_SUBMATRIX(i_v,j_v        )
           
          IF(HaveGradVecCoef(i_v,j_v))&
               &     CALL T2EXEC_GV_SUBMATRIX(i_v,j_v        )
          
          IF(HaveGradTenCoef(i_v,j_v))THEN
             DO i_k = 1, NKMAX
                IF(HaveGradTenKval(i_k,i_v,j_v))&
                     CALL T2EXEC_GT_SUBMATRIX(i_v,j_v,i_k    )
             ENDDO
          ENDIF
          
          IF(HaveExciScaCoef(i_v,j_v))&
               &     CALL T2EXEC_ES_SUBMATRIX(i_v,j_v        )
          
          IF(HaveExciVecCoef(i_v,j_v))THEN
             DO i_k = 1, NKMAX
                IF(HaveExciVecKval(i_k,i_v,j_v))&
                     CALL T2EXEC_EV_SUBMATRIX(i_v,j_v,i_k    )
             ENDDO
          ENDIF
          
          IF(HaveExciTenCoef(i_v,j_v))THEN          
             DO j_k = 1, NKMAX
             DO i_k = 1, NKMAX
                IF(HaveExciTenKval(i_k,j_k,i_v,j_v))&
                     CALL T2EXEC_ET_SUBMATRIX(i_v,j_v,i_k,j_k)
             ENDDO
             ENDDO
          ENDIF
          
          IF(HaveSourScaCoef(i_v,j_v))&
               &     CALL T2EXEC_SS_SUBVECTOR(i_v,j_v        )
          
       ENDDO
       ENDDO

       SELECT CASE(CoordinateSwitch)
          
       CASE (1)
          CALL T2EXEC_STORE
       CASE (2)
          CALL T2EXEC_STORE_TEST
       END SELECT
    
    ENDDO
    
    !CALL T2EXEC_CHECK
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)')&
         '-- T2EXEC: matrix construction was completed: cpu=', &
         e0time_1-e0time_0,' [s]'
    

    CALL CPU_TIME(e0time_0)
    ! set equations to be locked

    CALL T2EXEC_LOCK_VALUES(StartEqs,EndEqs,LockEqs)
    ! set dirichlet boundary condition on magnetic axis
    CALL T2EXEC_LOCK_VALUES(StartAxi,EndAxi,LockAxi)
    ! set dirichlet boundary condition on first wall
    CALL T2EXEC_LOCK_VALUES(StartWal,EndWal,LockWal)
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)')&
         '-- T2EXEC: boundary setting was completed:    cpu=', &
         e0time_1-e0time_0,' [s]'
    
    CALL CPU_TIME(e0time_0)
    
    CALL T2EXEC_SOLVE
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_SOLVE completed:  cpu=', &
         e0time_1-e0time_0,' [s]'
    
    CALL T2EXEC_DEALLOCATE

    RETURN
    
  END SUBROUTINE T2EXEC_EXECUTE
  
  SUBROUTINE T2EXEC_INITIALIZE_MATRIX
    
    USE T2COMM,ONLY:&
         NVMAX,NAMAX,NBMAX,NPMIN,NVFMX,&
         NodeDiaCRS,CoordinateSwitch
    
    INTEGER(ikind)::i_b,i_v,i_a
    amat(1:NVMAX,1:NVMAX,1:NAMAX) = 0.D0
    bvec(1:NVMAX,1:NBMAX)         = 0.D0
    
    SELECT CASE (CoordinateSwitch)
    CASE (1)
       DO i_b = 1,NBMAX
          IF((i_b.NE.1).AND.(MOD(i_b-1,NPMIN).NE.0))THEN
             i_a = NodeDiaCRS(i_b)
             DO i_v =1,NVFMX
                amat(i_v,i_v,i_a) = 1.D0
                !print*,i_v,i_v,i_a,i_b,NVFMX
             ENDDO
          ENDIF
       ENDDO
    CASE(2)
       DO i_b = 1,NBMAX
          IF(MOD(i_b,NPMIN).NE.0)THEN
             i_a = NodeDiaCRS(i_b)
             DO i_v =1,NVFMX
                amat(i_v,i_v,i_a) = 1.D0
             ENDDO
          END IF
       ENDDO
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2EXEC_INITIALIZE_MATRIX
  
  SUBROUTINE T2EXEC_CHECK
    
    USE T2COMM, ONLY:NAMAX,NBMAX,NVMAX,NXMAX,Xvec
    
    INTEGER(ikind)::i_v,j_v,i_a,i_b,i_x
    
110 FORMAT('AMAT',1X,'iv=',I2,1X,'jv=',I2,1X,'ia=',I8,1X,'val=',D15.6)
120 FORMAT('BVEC',1X,'iv=',I2,1X,'jv=',I2,1X,'ia=',I8,1X,'val=',D15.6)
130 FORMAT('XVEC',1X,'iv=',I2,1X,'jv=',I2,1X,'ix=',I8,1X,'val=',D15.6)

    OPEN(32,FILE='TEST_AMAT.txt')

    DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(32,110)i_v,j_v,12,amat(i_v,j_v,12)
       ENDDO
    ENDDO

    CLOSE(32)
    !OPEN(32,FILE='TEST_BVEC.txt')
   ! 
   ! DO i_v = 1, NVMAX
   !    !DO j_v = 1, NVMAX
   !    WRITE(32,120)i_v,0,7,bvec(i_v,7)
   !    !ENDDO
   ! ENDDO
   
    !CLOSE(32)
    !OPEN(32,FILE='TEST_XVEC.txt')
   ! 
   ! DO i_v = 1, NVMAX
   !    !DO j_v = 1, NVMAX
   !    WRITE(32,130)i_v,0,7,Xvec(i_v,7)
   !    !ENDDO
   ! ENDDO
   ! CLOSE(32)
    !OPEN(32,FILE='TEST_XVEC.txt')
    !DO i_x = 1,NXMAX
    !   DO i_v = 1, NVMAX
    !   !DO j_v = 1, NVMAX
    !      WRITE(32,130)i_v,0,i_a,Xvec(i_v,i_x)
    !   !ENDDO
    !   ENDDO
    !ENDDO
    !CLOSE(32)
    !STOP
    ! <<<<< ***** for debug  ***** <<<<<
    
    RETURN

  END SUBROUTINE T2EXEC_CHECK
  !-------------------------------------------------------------------
  !
  !       ALLOCATOR OF GLOBAL VARIABLES  FOR T2EXEC
  !
  !        *  Their scopes are inside of T2EXEC 
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_ALLOCATE
    
    USE T2COMM,ONLY: NDMAX,NNMAX,NVMAX,NKMAX,NBMAX,NXMAX,NAMAX

    INTEGER(ikind):: ierr
    
    DO 
       ALLOCATE(nodeTableA(  1:NNMAX),STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(nodeTableB(  1:NNMAX),STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(nodeTableC(  1:NNMAX),STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(nodeTableD(  1:NNMAX),STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(jacInvLocCrd(1:NDMAX,1:NDMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(knownVarElem(1:NNMAX,1:NKMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(bvecElem(    1:NNMAX,1:NVMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(amatElem(    1:NNMAX,1:NNMAX,1:NVMAX,1:NVMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(bvecElemTF(  1:NVMAX,1:NNMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(amatElemTF(  1:NVMAX,1:NVMAX,1:NNMAX,1:NNMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(amat(        1:NVMAX,1:NVMAX,1:NAMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       ALLOCATE(bvec(        1:NVMAX,1:NBMAX),&
            &                         STAT=ierr); IF(ierr.NE.0) EXIT
       
       nodeTableA(1:NNMAX) = 0
       nodeTableB(1:NNMAX) = 0
       nodeTableC(1:NNMAX) = 0
       nodeTableD(1:NNMAX) = 0

       jacInvLocCrd(1:NDMAX,1:NDMAX)     = 0.D0
       knownVarElem(1:NNMAX,1:NKMAX)     = 0.D0
       bvecElem(    1:NNMAX,1:NVMAX)     = 0.D0
       amatElem(    1:NNMAX,1:NNMAX,1:NVMAX,1:NVMAX)  = 0.D0
       bvecElemTF(  1:NVMAX,1:NNMAX)                  = 0.D0
       amatElemTF(  1:NVMAX,1:NVMAX,1:NNMAX,1:NNMAX ) = 0.D0
       amat(        1:NVMAX,1:NVMAX,1:NAMAX)          = 0.D0
       bvec(        1:NVMAX,1:NBMAX)                  = 0.D0

       RETURN
       
    ENDDO
    
    WRITE(6,'(A)') 'T2EXEC_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
    
    STOP
    
  END SUBROUTINE T2EXEC_ALLOCATE

  !-------------------------------------------------------------------
  !
  !       DEALLOCATOR OF GLOBAL VARIABLES FOR T2EXEC
  !
  !        *  Their scopes are limited to the inside of T2EXEC 
  !
  !                                     LAST UPDATE 2014-05-27
  ! 
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_DEALLOCATE
    
    IF(ALLOCATED(nodeTableA))   DEALLOCATE(nodeTableA)
    IF(ALLOCATED(nodeTableB))   DEALLOCATE(nodeTableB)
    IF(ALLOCATED(nodeTableC))   DEALLOCATE(nodeTableC)
    IF(ALLOCATED(nodeTableD))   DEALLOCATE(nodeTableD)
    IF(ALLOCATED(jacInvLocCrd)) DEALLOCATE(jacInvLocCrd)
    IF(ALLOCATED(knownVarElem)) DEALLOCATE(knownVarElem)
    IF(ALLOCATED(bvecElem))     DEALLOCATE(bvecElem)
    IF(ALLOCATED(amatElem))     DEALLOCATE(amatElem)
    IF(ALLOCATED(bvecElemTF))   DEALLOCATE(bvecElemTF)
    IF(ALLOCATED(amatElemTF))   DEALLOCATE(amatElemTF)
    IF(ALLOCATED(amat))         DEALLOCATE(amat)
    IF(ALLOCATED(bvec))         DEALLOCATE(bvec)
    
    RETURN
    
  END SUBROUTINE T2EXEC_DEALLOCATE
  
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
  SUBROUTINE T2EXEC_SETUP_ELEMENT_VARIABLES(i_e)
    
    USE T2COMM,ONLY:NNMAX,NDMAX,NVMAX,NKMAX,&
         &         KnownVar,ElementNodeGraph,GlobalCrd
    
    INTEGER(ikind),INTENT(IN)::i_e
    
    INTEGER(ikind)::i_r,i_m,i_k,i_n,j_n,k_n
    REAL(   rkind)::gridSizeRad,gridSizePol
        
    !
    ! SET LOCAL NODE-ELEMENT GRAPH
    !
    !   ElementNodeGraphElem(1:NNMAX,1): FOR COEFFICIENT CALCULATION
    !   ElementNodeGraphElem(1:NNMAX,2): FOR 2D-2D 
    !   ElementNodeGraphElem(1:NNMAX,3): FOR 1D-2D,1D-1D 
    !   ElementNodeGraphElem(1:NNMAX,4): FOR 2D-1D
    !
    
    DO i_n = 1, NNMAX
       nodeTableA(i_n) = ElementNodeGraph(i_n,1,i_e)
       nodeTableB(i_n) = ElementNodeGraph(i_n,2,i_e)
       nodeTableC(i_n) = ElementNodeGraph(i_n,3,i_e)
       nodeTableD(i_n) = ElementNodeGraph(i_n,4,i_e)
    ENDDO
    
    !
    ! CALCULATE JACOBIAN OF LOCAL COORDINATE 
    ! 
    ! JacDetLocCrd: JACOBIAN OF LOCAL COORDINATES
    ! JacInvLocCrd: INVERSE JACOBI MATRIX OF LOCAL COORDINATES
    ! 
    
    i_n = nodeTableA(1)
    j_n = nodeTableA(2)
    k_n = nodeTableA(4)
    
    gridSizeRad = GlobalCrd(1,j_n)-GlobalCrd(1,i_n)
    gridSizePol = GlobalCrd(2,k_n)-GlobalCrd(2,i_n)

    gridSizePol  = ABS(gridSizePol)
    gridSizeRad  = ABS(gridSizeRad)
    
    jacDetLocCrd = gridSizeRad*gridSizePol*0.25D0

    IF(JacDetLocCrd.LE.0.D0)THEN
       WRITE(6,'("ERROR:: Jacobian of Local Coords is singular")')
       STOP
    ENDIF
    
    jacInvLocCrd(1,1) = 2.D0/gridSizeRad
    jacInvLocCrd(1,2) = 0.D0
    jacInvLocCrd(2,1) = 0.D0
    jacInvLocCrd(2,2) = 2.D0/gridSizePol
    
    !
    ! STORE WORKING ARRAY FOR DIFFERENTIAL
    !
    ! KnownVarElem(:,1   ) : B  AT L-TH PICARD ITERATION
    ! KnownVarElem(:,2   ) : R  AT L-TH PICARD ITERATION
    ! KnownVarElem(:,2N+1) : Ub AT L-TH PICARD ITERATION 
    ! KnownVarElem(:,2N+2) : Wb AT L-TH PICARD ITERATION 
    !
    
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO i_k = 1, NKMAX
          knownVarElem(i_n,i_k) = KnownVar(i_k,i_m)
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_SETUP_ELEMENT_VARIABLES
  
  !-------------------------------------------------------------------
  !
  !    MASS SCALAR SUBMATRCES (PRIVATE)
  !
  !                     2014-06-09 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_MS_SUBMATRIX(i_v,j_v)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,NVFMX,CoordinateSwitch,&
         &           MassScaIntgPG,MassScaCoef,Xvec,TestCase
    
    INTEGER(ikind),INTENT(IN):: i_v,j_v

    INTEGER(ikind)::&
         i_n,j_n,k_n,i_m,i_b,i_x

    REAL(   rkind)::&
         massScaMatElem( 1:NNMAX,1:NNMAX),&
         massScaVecElem( 1:NNMAX        ),&
         massScaCoefElem(1:NNMAX        ),&
         xvecElem(       1:NNMAX        )
    
    ! initialization   
    DO i_n = 1, NNMAX
       i_m = nodeTableA(  i_n)
       massScaCoefElem(   i_n            ) &
            = MassScaCoef(    i_v,j_v,i_m) &
            * JacDetLocCrd
    ENDDO
    
    SELECT CASE(CoordinateSwitch)
    CASE (1)
       IF(    (j_v.GE.1    ).AND.(j_v.LE.NVFMX))THEN
          DO i_n = 1, NNMAX
             i_b = nodeTableD(i_n)
             xvecElem(i_n) = Xvec(j_v,i_b)
             !IF((i_elm.EQ.4).AND.(i_v.EQ.4).AND.(j_v.EQ.2))THEN
             !   print*,i_n,i_b,i_v,j_v,Xvec(j_v,i_b)
             !END IF
          ENDDO
       ELSEIF((j_v.GT.NVFMX).AND.(j_v.LE.NVMAX))THEN
          DO i_n = 1, NNMAX
             i_x = nodeTableB(i_n)
             xvecElem(i_n) = Xvec(j_v,i_x)
          ENDDO
       ENDIF
    CASE (2)
       SELECT CASE (TestCase)
       CASE(1:2)          
          DO i_n = 1, NNMAX
             i_x = nodeTableB(i_n)
             xvecElem(i_n) = Xvec(j_v,i_x)          
          ENDDO
       CASE (3)
          DO i_n = 1, NNMAX
             i_b = nodeTableD(i_n)
             xvecElem(i_n) = Xvec(j_v,i_b)
          ENDDO
       END SELECT
    END SELECT

    ! main loop    
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
               * xvecElem(          j_n)
       ENDDO
    ENDDO
    
    ! store submatrix
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
            + massScaMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    ! store subvector
    DO i_n = 1, NNMAX
       BvecElem(             i_n,i_v) &
            = BvecElem(      i_n,i_v) & 
            + massScaVecElem(i_n    )
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_MS_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !    ADVECTION VECTOR SUBMATRCES (PRIVATE)
  !
  !                 2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_AV_SUBMATRIX(i_v,j_v)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           AdveVecIntgPG,AdveVecCoef

    INTEGER(ikind),INTENT(IN)::i_v,j_v
    INTEGER(ikind)::&
         i_d,j_d,&
         i_n,j_n,k_n,&
         i_m
    
    REAL(   rkind)::&
         adveVecMatElem( 1:NNMAX,1:NNMAX),&
         adveVecCoefElem(1:NDMAX,1:NNMAX),&
         adveVecCoefLoc( 1:NDMAX,1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO i_d = 1, NDMAX
          adveVecCoefElem(   i_d,i_n            ) &
               = AdveVecCoef(i_d,    i_v,j_v,i_m) &
               * jacDetLocCrd
          
       ENDDO
    ENDDO

    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX    
          adveVecCoefLoc(j_d,k_n) = 0.D0
          DO i_d = 1, NDMAX
             adveVecCoefLoc(            j_d,k_n) &
                  = adveVecCoefLoc(     j_d,k_n) &
                  + adveVecCoefElem(i_d,    k_n) &
                  * jacInvLocCrd(   i_d,j_d    ) 
          ENDDO
       ENDDO
    ENDDO
    
    ! main loop
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

    ! store submatrix
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + adveVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    !IF((i_v.Eq.5).AND.(j_v.Eq.5))THEN
    !   DO i_n = 1, NNMAX
    !      DO j_n = 1, NNMAX
    !         print*,i_m,i_v,j_v,amatElem(             i_n,j_n,i_v,j_v)
    !      ENDDO
    !   ENDDO
    !END IF
    
    RETURN
    
  END SUBROUTINE T2EXEC_AV_SUBMATRIX
  
  !-------------------------------------------------------------------
  ! 
  !    ADVECTION TENSOR SUBMATRCES (PRIVATE)
  ! 
  !                 2014-05-22 H.SETO
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_AT_SUBMATRIX(i_v,j_v,i_k)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           AdveTenIntgPG,AdveTenCoef

    INTEGER(ikind),INTENT(IN)::i_v,j_v,i_k
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         adveTenCoefElem(1:NDMAX,1:NDMAX,1:NNMAX),&
         adveTenCoefLoc( 1:NDMAX,1:NDMAX,1:NNMAX),&
         adveTenMatElem( 1:NNMAX,1:NNMAX),&
         adveTenKvarElem(1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       adveTenKvarElem(i_n) = knownVarElem(i_n,i_k)
    ENDDO
    
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          adveTenCoefElem(   i_d,j_d,i_n                ) &
               = AdveTenCoef(i_d,j_d,    i_k,i_v,j_v,i_m) &
               * jacDetLocCrd
       ENDDO
       ENDDO
    ENDDO
    
    DO l_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          adveTenCoefLoc(k_d,l_d,l_n) = 0.D0
          DO j_d = 1, NDMAX        
          DO i_d = 1, NDMAX
             adveTenCoefLoc(                k_d,l_d,l_n) &
                  = adveTenCoefLoc(         k_d,l_d,l_n) &
                  + adveTenCoefElem(i_d,j_d,        l_n) &
                  * jacInvLocCrd(   i_d,    k_d        ) &
                  * jacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
         
    ! main loop
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       adveTenMatElem(i_n,j_n) = 0.D0
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             adveTenMatElem(                        i_n,j_n) &
                  = adveTenMatElem(                 i_n,j_n) &
                  + AdveTenIntgPG(  k_d,l_d,k_n,l_n,i_n,j_n) &
                  * adveTenKvarElem(        k_n            ) &
                  * adveTenCoefLoc( k_d,l_d,    l_n        ) 
             
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    ! store submatrix
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + adveTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AT_SUBMATRIX
  
  !-------------------------------------------------------------------
  !
  !  DIFFUSION TENSOR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_DT_SUBMATRIX(i_v,j_v)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           DiffTenIntgPG,DiffTenCoef
    
    INTEGER(ikind),INTENT(IN)::i_v,j_v
    INTEGER(ikind)::&
         i_n,j_n,k_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         diffTenMatElem( 1:NNMAX,1:NNMAX),&
         diffTenCoefElem(1:NDMAX,1:NDMAX,1:NNMAX), &
         diffTenCoefLoc( 1:NDMAX,1:NDMAX,1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          diffTenCoefElem(   i_d,j_d,i_n            ) &
               = diffTenCoef(i_d,j_d,    i_v,j_v,i_m) &
               * jacDetLocCrd
       ENDDO
       ENDDO
    ENDDO

    DO k_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          diffTenCoefLoc(k_d,l_d,k_n) = 0.D0
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             diffTenCoefLoc(                k_d,l_d,k_n) &
                  = diffTenCoefLoc(         k_d,l_d,k_n) &
                  + diffTenCoefElem(i_d,j_d,        k_n) &
                  * jacInvLocCrd(   i_d,    k_d        ) &
                  * jacInvLocCrd(       j_d,    l_d    )
          END DO
          END DO
       END DO
       END DO
    END DO

    ! main loop

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
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + diffTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_DT_SUBMATRIX
  
  !-------------------------------------------------------------------
  !
  !   GRADIENT VECTOR SUBMATRICES (PRIVATE)
  !
  !                 2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GV_SUBMATRIX(i_v,j_v)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           GradVecIntgPG,GradVecCoef
    
    INTEGER(ikind),INTENT(IN)::i_v,j_v
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
       i_m = nodeTableA(i_n)
       DO i_d = 1, NDMAX
          gradVecCoefElem(   i_d,i_n            ) &
               = GradVecCoef(i_d,    i_v,j_v,i_m) &
               * jacDetLocCrd
       ENDDO
    ENDDO
    
    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX
          gradVecCoefLoc(j_d,k_n) = 0.D0
          DO i_d = 1, NDMAX
             gradVecCoefLoc(            j_d,k_n) &
                  = gradVecCoefLoc(     j_d,k_n) &
                  + gradVecCoefElem(i_d,    k_n) &
                  * jacInvLocCrd(   i_d,j_d    ) 
          ENDDO
       ENDDO
    ENDDO

    ! main loop
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
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
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + gradVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GV_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !   GRADIENT TENSOR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-27 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GT_SUBMATRIX(i_v,j_v,i_k)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           GradTenIntgPG,GradTenCoef
    
    INTEGER(ikind),INTENT(IN)::i_v,j_v,i_k
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         gradTenMatElem( 1:NNMAX,1:NNMAX),&
         gradTenCoefElem(1:NDMAX,1:NDMAX,1:NNMAX),&
         gradTenCoefLoc( 1:NDMAX,1:NDMAX,1:NNMAX),&
         gradTenKvarElem(1:NNMAX)
    
    ! initialization  
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          gradTenCoefElem(   i_d,j_d,i_n                ) &
               = GradTenCoef(i_d,j_d,    i_k,i_v,j_v,i_m) &
               * jacDetLocCrd
       ENDDO
       ENDDO
    ENDDO
    
    DO l_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          gradTenCoefLoc(k_d,l_d,l_n) = 0.D0
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             gradTenCoefLoc(                k_d,l_d,l_n) &
                  = gradTenCoefLoc(         k_d,l_d,l_n) &
                  + gradTenCoefElem(i_d,j_d,        l_n) &
                  * jacInvLocCrd(   i_d,    k_d        ) &
                  * jacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    
    DO i_n = 1, NNMAX
       gradTenKvarElem(i_n) = knownVarElem(i_n,i_k)
    ENDDO
    
    ! main loop
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       gradTenMatElem(i_n,j_n) = 0.D0
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             gradTenMatElem(i_n,j_n)&
                  = gradTenMatElem(                 i_n,j_n) &
                  + GradTenIntgPG(  k_d,l_d,k_n,l_n,i_n,j_n) &
                  * gradTenKvarElem(        k_n            ) &
                  * gradTenCoefLoc( k_d,l_d,    l_n        ) 
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO 
    
    ! store matrix
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + gradTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GT_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !  EXCITATION SCALAR SUBMATRCES (PRIVATE)
  !
  !                 2014-05-22 H.SETO
  !
  !-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_ES_SUBMATRIX(i_v,j_v)

    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           ExciScaIntgPG,ExciScaCoef
         
    INTEGER(ikind),INTENT(IN)::i_v,j_v
    INTEGER(ikind)::&
         i_n,j_n,k_n,&
         i_m


    REAL(   rkind)::&
         exciScaMatElem( 1:NNMAX,1:NNMAX),&
         exciScaCoefElem(1:NNMAX        )

    ! initialization
        
    DO i_n = 1,NNMAX
       i_m = nodeTableA(i_n)
       exciScaCoefElem(   i_n            ) &
            = ExciScaCoef(    i_v,j_v,i_m) &
            * jacDetLocCrd
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
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + exciScaMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_ES_SUBMATRIX

  !------------------------------------------------------------------
  !
  !    EXCITATION SUBMATRCES (PRIVATE)
  !
  !                 2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_EV_SUBMATRIX(i_v,j_v,i_k)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           ExciVecIntgPG,ExciVecCoef
    
    INTEGER(ikind),INTENT(IN)::i_v,j_v,i_k
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,&
         i_d,j_d,&
         i_m

    REAL(   rkind)::&
         exciVecMatElem( 1:NNMAX,1:NNMAX),&
         exciVecCoefElem(1:NDMAX,1:NNMAX),&
         exciVecCoefLoc( 1:NDMAX,1:NNMAX),&    
         exciVecKvarElem(1:NNMAX)
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO i_d = 1, NDMAX
          exciVecCoefElem(   i_d,i_n                ) &
               = ExciVecCoef(i_d,    i_k,i_v,j_v,i_m) &
               * jacDetLocCrd
       ENDDO
    ENDDO
    
    DO l_n = 1, NNMAX   
       DO j_d = 1, NDMAX
          exciVecCoefLoc(j_d,l_n) = 0.D0     
          DO i_d = 1, NDMAX
             exciVecCoefLoc(            j_d,l_n) &
                  = exciVecCoefLoc(     j_d,l_n) &
                  + exciVecCoefElem(i_d,    l_n) &
                  * jacInvLocCrd(   i_d,j_d    )
          ENDDO
       ENDDO
    ENDDO
        
    DO i_n = 1, NNMAX
       exciVecKvarElem(i_n) = knownVarElem(i_n,i_k)
    ENDDO
    
    ! main loop
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       exciVecMatElem(i_n,j_n) = 0.D0
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO j_d = 1, NDMAX
             exciVecMatElem(                    i_n,j_n) &
                  = exciVecMatElem(             i_n,j_n) &
                  + ExciVecIntgPG(  j_d,k_n,l_n,i_n,j_n) &
                  * exciVecKvarElem(    k_n            ) &
                  * exciVecCoefLoc( j_d,    l_n        )
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
        
    ! store submatrix
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + exciVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_EV_SUBMATRIX
  
  !-------------------------------------------------------------------
  !
  !     EXCITATION TENESOR SUBMATRCES (PRIVATE)
  !
  !                 2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_ET_SUBMATRIX(i_v,j_v,i_k,j_k)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           ExciTenIntgPG,ExciTenCoef
    
    INTEGER(ikind),INTENT(IN)::i_v,j_v,i_k,j_k
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,m_n,&
         i_d,j_d,k_d,l_d,&
         i_m
    
    REAL(   rkind)::&
         exciTenMatElem(  1:NNMAX,1:NNMAX),&
         exciTenCoefElem( 1:NDMAX,1:NDMAX,1:NNMAX),&
         exciTenCoefLoc(  1:NDMAX,1:NDMAX,1:NNMAX),&    
         exciTenKvar1Elem(1:NNMAX),&
         exciTenKvar2Elem(1:NNMAX)

    ! initialization
  
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          exciTenCoefElem(   i_d,j_d,i_n)&
               = ExciTenCoef(i_d,j_d,    i_k,j_k,i_v,j_v,i_m) &
               * jacDetLocCrd
       ENDDO
       ENDDO
    ENDDO
    
    DO m_n = 1, NNMAX
       DO l_d = 1, NDMAX
       DO k_d = 1, NDMAX
          exciTenCoefLoc(k_d,l_d,m_n) = 0.D0
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             exciTenCoefLoc(                k_d,l_d,m_n) &
                  = exciTenCoefLoc(         k_d,l_d,m_n) &
                  + exciTenCoefElem(i_d,j_d,        m_n) &
                  * jacInvLocCrd(   i_d,    k_d        ) &
                  * jacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO

    DO i_n = 1, NNMAX
       exciTenKvar1Elem(i_n) = KnownVarElem(i_n,i_k)
       exciTenKvar2Elem(i_n) = KnownVarElem(i_n,j_k)
    ENDDO
    
    ! main loop
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       exciTenMatElem(i_n,j_n) = 0.D0
       DO m_n = 1, NNMAX
       DO l_n = 1, NNMAX
       DO k_n = 1, NNMAX
          DO l_d = 1, NDMAX
          DO k_d = 1, NDMAX
             exciTenMatElem(                             i_n,j_n) &
                  = exciTenMatElem(                      i_n,j_n) &
                  + ExciTenIntgPG(   k_d,l_d,k_n,l_n,m_n,i_n,j_n) &
                  * exciTenKvar1Elem(        k_n                ) &
                  * exciTenKvar2Elem(            l_n            ) &
                  * exciTenCoefLoc(  k_d,l_d,        m_n        )
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
       amatElem(             i_n,j_n,I_v,J_v) &
            = amatElem(      i_n,j_n,I_v,J_v) &
            + exciTenMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_ET_SUBMATRIX

  !
  !
  !
  !
  SUBROUTINE T2EXEC_SS_SUBVECTOR(i_v,j_v)
    
    USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,&
         &           SourScaIntgPG,SourScaCoef
    
    INTEGER(ikind),INTENT(IN)::i_v,j_v
    INTEGER(ikind)::&
         i_n,j_n,&
         i_m
    
    REAL(   rkind)::&
         sourScaVecElem( 1:NNMAX),&
         sourScaCoefElem(1:NNMAX)
    
    ! initialization
    DO i_n = 1, NNMAX
       i_m = nodeTableA(i_n)
       sourScaCoefElem(   i_n            ) &
            = SourScaCoef(    i_v,j_v,i_m) &
            * jacDetLocCrd
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
       bvecElem(             i_n,i_v) &
            = bvecElem(      i_n,i_v) & 
            + sourScaVecElem(i_n    )
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_SS_SUBVECTOR
  
  !-------------------------------------------------------------------
  !
  ! SUBROUTINE FOR STORE SUBMATRIX 
  ! FOR BI-LINEAR RECTANGULAR ELEMENT
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_STORE
    
    USE T2COMM,ONLY: NNMAX,NVMAX,NVFMX,HaveMat,&
         &           NodeRowCRS,NodeColCRS,HangedNodeTable
    INTEGER(ikind)::&
         i_n,j_n,&
         i_v,j_v
    
    ! transform AmatElem to AmatElemTF and clear AmatElem
    
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       IF(HaveMat(i_v,j_v))THEN
          DO j_n = 1, NNMAX
          DO i_n = 1, NNMAX
             AmatElemTF(i_v,j_v,i_n,j_n) = AmatElem(i_n,j_n,i_v,j_v)
             AmatElem(  i_n,j_n,i_v,j_v) = 0.D0
          ENDDO
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    
    ! transform BvecElem to BvecElemTF and clear BvecElem  
    
    DO I_v = 1, NVMAX
       DO i_n = 1, NNMAX
          BvecElemTF(i_v,i_n) = BvecElem(i_n,i_v)
          BvecElem(  i_n,i_v) = 0.D0
       ENDDO
    ENDDO
    
    !
    ! for stiffness matrix
    ! i_v: 1D value (FSA)
    ! j_v: 1D value (FSA)
    CALL T2EXEC_STORE_MATRIX(1      ,NVFMX,nodeTableD,& ! row
         &                   1      ,NVFMX,nodeTableD)  ! column
    ! i_v: 1D value (FSA)
    ! j_v: 2D value
    CALL T2EXEC_STORE_MATRIX(1      ,NVFMX,nodeTableD,& ! row
         &                   NVFMX+1,NVMAX,nodeTableB)  ! column
    ! i_v: 2D value 
    ! j_v: 1D value (FSA) 
    CALL T2EXEC_STORE_MATRIX(NVFMX+1,NVMAX,nodeTableB,& ! row
         &                   1      ,NVFMX,nodeTableD)  ! column
    !CALL T2EXEC_STORE_MATRIX(NVFMX+1,NVMAX,nodeTableB,& ! row
    !     &                   1      ,NVFMX,nodeTableB)  ! column
    ! i_v: 2D value
    ! j_v: 2D value
    CALL T2EXEC_STORE_MATRIX(NVFMX+1,NVMAX,nodeTableB,& ! row
         &                   NVFMX+1,NVMAX,nodeTableB)  ! column

    !
    ! for RHS vector
    ! i_v: 1D value (FSA)
    CALL T2EXEC_STORE_VECTOR(1      ,NVFMX,nodeTableD) ! row
    ! i_v 2D value
    CALL T2EXEC_STORE_VECTOR(NVFMX+1,NVMAX,nodeTableB) ! row
    RETURN
    
  END SUBROUTINE T2EXEC_STORE

  !-------------------------------------------------------------------
  !
  ! SUBROUTINE FOR STORE SUBMATRIX 
  ! FOR BI-LINEAR RECTANGULAR ELEMENT
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_STORE_TEST
    
    USE T2COMM,ONLY: NNMAX,NVMAX, HaveMat,TestCase,&
         &           NodeRowCRS,NodeColCRS,HangedNodeTable
    INTEGER(ikind)::&
         i_n,j_n,i_v,j_v
    
    ! transform AmatElem to AmatElemTF and clear AmatElem
    
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       IF(HaveMat(i_v,j_v))THEN
          DO j_n = 1, NNMAX
          DO i_n = 1, NNMAX
             AmatElemTF(i_v,j_v,i_n,j_n) = AmatElem(i_n,j_n,i_v,j_v)
             AmatElem(  i_n,j_n,i_v,j_v) = 0.D0
          ENDDO
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    
    ! transform BvecElem to BvecElemTF and clear BvecElem  
    
    DO I_v = 1, NVMAX
       DO i_n = 1, NNMAX
          BvecElemTF(i_v,i_n) = BvecElem(i_n,i_v)
          BvecElem(  i_n,i_v) = 0.D0
       ENDDO
    ENDDO
    
    ! set node tables 
        
    ! for stiffness matrix
    SELECT CASE (TestCase)
    CASE (1:2)
       CALL T2EXEC_STORE_MATRIX(1    ,1    ,nodeTableB,& ! row
            &                   1    ,1    ,nodeTableB)  ! column
       ! for RHS vector 
       CALL T2EXEC_STORE_VECTOR(1    ,1    ,nodeTableB) ! row
    CASE (3)
       CALL T2EXEC_STORE_MATRIX(1    ,1    ,nodeTableD,& ! row
            &                   1    ,1    ,nodeTableD)  ! column
       CALL T2EXEC_STORE_VECTOR(1    ,1    ,nodeTableD)  ! row
    END SELECT
    RETURN
    
  END SUBROUTINE T2EXEC_STORE_TEST

  !-------------------------------------------------------------------
  !
  !  T2EXEC_STORE_MATRIX
  !
  !         store submatrix into global stiffness matrix 
  !
  !                                node     data  CRS-format 
  !                                variable data  uncompressed-format
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_STORE_MATRIX(rowValStart,rowValEnd,rowNodeTable,&
       &                         colValStart,colValEnd,colNodeTable)

    USE T2COMM,ONLY: NNMAX,NBMAX,HaveMat,&
         &           NodeRowCRS, NodeColCRS,HangedNodeTable
    
    INTEGER,INTENT(IN)::&
         rowValStart,rowValEnd,rowNodeTable(1:NNMAX),&
         colValStart,colValEnd,colNodeTable(1:NNMAX)
    
    INTEGER(ikind)::&
         i_n,j_n,i_v,j_v,i_a,&
         !
         i_rowN,j_rowN,k_rowN,&
         i_colN,j_colN,k_colN,l_colN

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       
       i_rowN = rowNodeTable(i_n)
       i_colN = colNodeTable(j_n)
       
       IF((i_rowN.LE.NBMAX).AND.(i_colN.LE.NBMAX))THEN
          
          DO i_a = NodeRowCRS(i_rowN), NodeRowCRS(i_rowN+1)-1
             j_colN = NodeColCRS(i_a)
             IF(i_colN.EQ.j_colN)THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      !IF(i_elm.EQ.12)THEN
                         !print'(5(I3,1X),D15.8)',i_v,j_v,i_n,j_n,i_a,amatElemTF(i_v,j_v,i_n,j_n)
                      !END IF
                      amat(             i_v,j_v,        i_a) &
                           = amat(      i_v,j_v,        i_a) &
                           + amatElemTF(i_v,j_v,i_n,j_n    )
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
       ELSEIF((i_rowN.LE.NBMAX).AND.(i_colN.GT.NBMAX))THEN
          
          i_colN = i_colN - NBMAX
          j_colN = HangedNodeTable(1,i_colN)
          k_colN = HangedNodeTable(2,i_colN)
          
          DO i_a = NodeRowCRS(i_rowN), NodeRowCRS(i_rowN+1)-1
             l_colN = NodeColCRS(i_a)
             IF((l_colN.EQ.j_colN).OR.(l_colN.EQ.k_colN))THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      amat(             i_v,j_v,        i_a) &
                           = amat(      i_v,j_v,        i_a) &
                           + amatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.50D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i_rowN.GT.NBMAX).AND.(i_colN.LE.NBMAX))THEN
          
          i_rowN = i_rowN - NBMAX
          j_rowN = HangedNodeTable(1,i_rowN)
          k_rowN = HangedNodeTable(2,i_rowN)
          
          DO i_a = NodeRowCRS(j_rowN), NodeRowCRS(j_rowN+1)-1
             j_colN = NodeColCRS(i_a) 
             IF(j_colN.EQ.i_colN)THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      amat(             i_v,j_v,        i_a) &
                           = amat(      i_v,j_v,        i_a) &
                           + amatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.50D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
          DO i_a = NodeRowCRS(k_rowN), NodeRowCRS(k_rowN+1)-1
             j_colN = NodeColCRS(i_a) 
             IF(j_colN.EQ.i_colN)THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      amat(             i_v,j_v,        i_a) &
                           = amat(      i_v,j_v,        i_a) &
                           + amatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.50D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i_rowN.GT.NBMAX).AND.(i_colN.GT.NBMAX))THEN
          
          i_rowN = i_rowN - NBMAX
          j_rowN = HangedNodeTable(1,i_rowN)
          k_rowN = HangedNodeTable(2,i_rowN)
          
          i_colN = i_colN - NBMAX
          j_colN = HangedNodeTable(1,i_colN)
          k_colN = HangedNodeTable(2,i_colN)
          
          DO i_a = NodeRowCRS(j_rowN), NodeRowCRS(j_rowN+1)-1
             l_colN = NodeColCRS(i_a)
             IF((l_colN.EQ.j_colN).OR.(l_colN.EQ.k_colN))THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      amat(             i_v,j_v,        i_a) &
                           = amat(      i_v,j_v,        i_a) &
                           + amatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.25D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i_a = NodeRowCRS(k_rowN), NodeRowCRS(k_rowN+1)-1
             l_colN = NodeColCRS(i_a)
             IF((l_colN.EQ.j_colN).OR.(l_colN.EQ.k_colN))THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      amat(             i_v,j_v,        i_a) &
                           = amat(      i_v,j_v,        i_a) &
                           + amatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.25D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       
       ENDIF
       
    ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE T2EXEC_STORE_MATRIX

  !-------------------------------------------------------------------
  !
  !  T2EXEC_STORE_VECTOR
  !
  !         store subvector into global RHS vector
  !                                node     data  uncompressed-format 
  !                                variable data  uncompressed-format
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_STORE_VECTOR(rowValStart,rowValEnd,rowNodeTable)
    
    USE T2COMM,ONLY: NNMAX,NBMAX,HangedNodeTable
    INTEGER(ikind),INTENT(IN)::&
         rowValStart,rowValEnd,rowNodeTable(1:NNMAX)
    
    INTEGER(ikind)::&
         i_n,i_v,&
         i_row,j_row,k_row
    
    DO i_n = 1, NNMAX
       i_row = rowNodeTable(i_n)
       IF(    i_row.LE.NBMAX)THEN
          DO i_v = rowValStart,rowValEnd
             bvec(              i_v,    i_row) &
                  = bvec(       i_v,    i_row) &
                  + bvecElemTF( i_v,i_n      )
          ENDDO
       ELSEIF(i_row.GT.NBMAX)THEN
          i_row = i_row - NBMAX
          j_row = HangedNodeTable(1,i_row)
          k_row = HangedNodeTable(2,i_row)
          DO i_v = rowValStart, rowValEnd
             bvec(             i_v,    j_row) &
                  = bvec(      i_v,    j_row) &
                  + bvecElemTF(i_v,i_n    ) &
                  * 0.50D0
          ENDDO
          DO i_v = rowValStart, rowValEnd
             bvec(             i_v,    k_row) &
                  = bvec(      i_v,    k_row) &
                  + bvecElemTF(i_v,i_n      ) &
                  * 0.50D0
          ENDDO
       ENDIF
    ENDDO

    RETURN
    
  END SUBROUTINE T2EXEC_STORE_VECTOR
  
  !-------------------------------------------------------------------
  ! 
  !       T2EXEC_LOCK_VALUES
  !
  !           
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_LOCK_VALUES(nodeStart,nodeEnd,varLockTableIn)
    
    USE T2COMM,ONLY: NVMAX,NodeRowCRS, NodeColCRS,Xvec
    
    INTEGER(ikind),INTENT(IN)::nodeStart,nodeEnd
    LOGICAL       ,INTENT(IN)::varLockTableIn(1:NVMAX)
    
    LOGICAL::varLockTable(  1:NVMAX)
    INTEGER(ikind)::&
         ia,iv,jv,i_row,i_col
    
    varLockTable(1:NVMAX) = .TRUE.

    DO iv = 1,nvmx
       varLockTable(iv) = varLockTableIn(iv)
    ENDDO

    DO i_row = nodeStart, nodeEnd
       
       ! stiffness matrix
       DO ia = NodeRowCRS(i_row), NodeRowCRS(i_row+1)-1
          i_col = NodeColCRS(ia)
          DO jv = 1, NVMAX
          DO iv = 1, NVMAX
             IF(varLockTable(iv))THEN
                IF((i_row.EQ.i_col).AND.(iv.EQ.jv))THEN
                   amat(iv,jv,ia) = 1.D0
                ELSE
                   amat(iv,jv,ia) = 0.D0
                ENDIF
             ENDIF
          ENDDO
          ENDDO
       ENDDO
    
       ! RHS vector 
       DO iv = 1, NVMAX
          IF(varLockTable(iv))THEN
             bvec(iv,i_row) = Xvec(iv,i_row)
          ENDIF
       ENDDO
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_LOCK_VALUES
  
  !-------------------------------------------------------------------
  !
  !       SUBROUTINE T2EXEC_SOLVE
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_SOLVE
    
    USE T2COMM, ONLY:&
         NPMIN,NVMAX,NBMAX,NXMAX,NAMAX,NLMAX,NBVMX,NAVMX,&
         NodeRowCRS,NodeColCRS,HangedNodeTable,&
         i0dbg,idebug,i1pdn2,i1rdn2,HaveMat,&
         Xvec,XvecIn,XvecOut,HangedNodeTable,&
         !
         NVFMX,NFAMX,NodeFSA
    
    USE LIBMPI
    USE COMMPI
    USE LIBMTX

    INTEGER(ikind)::istart,iend,its
    INTEGER(ikind)::itype, m1, m2
    REAL(   rkind)::tolerance,val
    REAL(   rkind),POINTER,SAVE::x(:)
    INTEGER(ikind)::i0pdn2,i0rdn2
    INTEGER(ikind)::&
         i_rowN, j_rowN, k_rowN, &
         i_rowNV,j_rowNV,k_rowNV,&
         i_colN, i_colNV,i_rowV, i_colV,i_a,&
         i0offset,i0lidi,i0ridi,i0pidi
    
    INTEGER(ikind)::i_v,j_v,i_b,i_x,i1

    INTEGER(ikind)::im,iv,ifsa,jfsa

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
       
    ALLOCATE(x(NBVMX))
    x = 0.D0
    CALL MTX_SETUP(NBVMX,istart,iend,nzmax=NAVMX,idebug=0)
    !CALL MTX_SETUP(NBVMX,istart,iend,nzmax=NAVMX+2*NFAMX*NVFMX,idebug=0)
    ! 
    ! STORE GLOBAL STIFFNESS MATRIX  
    ! 
    !    OPEN(35,FILE="MATRIX.txt")
    DO i_rowN = 1, NBMAX   
       DO i_a = NodeRowCRS(i_rowN), NodeRowCRS(i_rowN+1)-1
          i_colN = NodeColCRS(i_a) 
          DO i_colV = 1, NVMAX
          DO i_rowV = 1, NVMAX
             IF(HaveMat(i_rowV,i_colV)) THEN
                val = amat(i_rowV,i_colV,i_a)
                i_rowNV  = NVMAX*(i_rowN-1) + i_rowV
                i_colNV  = NVMAX*(i_colN-1) + i_colV
                CALL MTX_SET_MATRIX(i_rowNV,i_colNV,val)
                !IF((i_rowV.EQ.1).AND.(i_colV.EQ.1))THEN
                !WRITE(35,'("in=",I5,1X,"jn=",I5,1X,"v=",D15.8)')&
                !                     i_rowN,i_colN,val
                !            ENDIF
                !IF(IDEBUG.EQ.1) &
                !                     & 
                !IF((i_rowV.EQ.7).OR.(i_rowV.Eq.17))&
                !WRITE(*,'(2I5,I10,2I3,2I7,1PE12.4)') &
                !     & i_rowN,i_colN,i_a,i_rowV,i_colV,i_rowNV,i_colNV,val
             END IF
          ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! >>>> for debug >>>>>
    !DO im = 1, NFAMX
    !   DO iv = 1, NVFMX
    !      ifsa = NodeFSA(1,im)*NVMAX
    !      jfsa = NodeFSA(2,im)*NVMAX
    !      CALL MTX_SET_MATRIX(ifsa+iv,ifsa+iv, 1.D0)
    !      CALL MTX_SET_MATRIX(ifsa+iv,jfsa+iv,-1.D0)
    !   ENDDO
    !ENDDO
    ! <<<< for debug <<<<<

    !CLOSE(35)
    ! SET GLOBAL RIGHT HAND SIDE VECTOR
    DO i_rowN = 1, NBMAX
       DO i_rowV = 1, NVMAX
          i_rowNV = NVMAX*(i_rowN-1) + i_rowV
          val  = bvec(i_rowV,i_rowN)
          CALL MTX_SET_SOURCE(i_rowNV,val)
       ENDDO
    ENDDO
    
    ! SET INITIAL GUESS
    DO i_rowN = 1, NBMAX
       DO i_rowV = 1, NVMAX
          i_rowNV  = NVMAX*(i_rowN-1) + i_rowV
          val = XvecIn(i_rowV,i_rowN)
          CALL MTX_SET_VECTOR(i_rowNV,val)
       ENDDO
    ENDDO

    CALL MTX_SOLVE(itype,tolerance,its,methodKSP=m1,methodPC=m2)
    
    CALL MTX_GATHER_VECTOR(x)
    
    DO i_rowN = 1, NXMAX
       IF(    i_rowN.LE.NBMAX)THEN
          DO i_rowV = 1, NVMAX
             i_rowNV = NVMAX*(i_rowN-1) + i_rowV
             XvecOut(i_rowV,i_rowN) = x(i_rowNV)
          ENDDO
       ELSEIF(i_rowN.GT.NBMAX)THEN
          k_rowN = i_rowN - NBMAX
          j_rowN = HangedNodeTable(1,k_rowN)
          k_rowN = HangedNodeTable(2,k_rowN)
          DO i_rowV = 1, NVMAX 
             j_rowNV = NVMAX*(j_rowN-1) + i_rowV
             k_rowNV = NVMAX*(k_rowN-1) + i_rowV
             XvecOut(i_rowV,i_rowN) = 0.5D0*(x(j_rowNV)+x(k_rowNV))
          ENDDO
       ENDIF
    ENDDO
    
    CALL MTX_CLEANUP
    
    DEALLOCATE(x)
    
    RETURN
    
  END SUBROUTINE T2EXEC_SOLVE
END MODULE T2EXEC
