!-------------------------------------------------------------------- 
!  MODULE FOR SOLVING ADVECTION-DIFFUSION PROBLEM BY FEM
! 
!                   LAST UPDATE 2014-05-20 H.SETO
!
!  T2INTG provides following variables:
!
!         xvec: 
!
!    and subroutines:
!
!         T2EXEC_EXECUTE
!         T2EXEC_TERMINATE
! 
!  T2INTG requires following variables:
!  
!  [from T2CNST]
!        rkind,ikind
!
!  [from T2COMM] 
!        NEMAX,NBMAX,NVMAX,NKMAX,NNMAX,NDMAX,Dt  
!
!  [from T2NGRA]
!        RowCRS, ColCRS, ElementNodeGraph, HangedNodeTable    
!
!  [from T2PROF] 
!        GlobalCrd
!
!  [from T2INTG]
!        MassScaIntgPG AdveVecIntgPG AdveTenIntgPG DiffTenIntgPG
!        GradVecIntgPG GradTenIntgPG ExciScaIntgPG ExciVecIntgPG
!        ExciTenIntgPG SourScaIntgPG

!  [from T2CALV]
!        MassScaCoef AdveVecCoef AdveTenCoef DiffTenCoef GradVecCoef
!        GradTenCoef ExciScaCoef ExciVecCoef ExciTenCoef SourScaCoef
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
       ElementNodeGraphElem(:,:)! node-element graph table in en element
  REAL(   rkind),SAVE::&
       JacDetLocCrd       ! jacobian of local coordinates
  REAL(   rkind),SAVE,ALLOCATABLE::&
       JacInvLocCrd(:,:),&! inverse jacobi matrix of local coordinates
       KnownVarElem(:,:),&! known variables at nodes in an element 
       BvecElem(:,:),&    ! element right hand side vector  (Ax=b)
       AmatElem(:,:,:,:),&! element stiffness matrix        (Ax=b)
       BvecElemTF(:,:),&
       AmatElemTF(:,:,:,:),&
       Amat(:,:,:),&
       Bvec(:,:),&
       Xvec(:,:)

  !
  ! for T2EXEC_ADHOC: they will be removed as soon as possible
  !
  
  ! from T2COMM
  INTEGER(ikind),SAVE::&
       NEMAX,NBMAX,NXMAX,NVMAX,NKMAX,NNMAX,NDMAX,NMMAX,Dt,&
       NHMAX,NAMAX,NNRMX,NLMAX

  ! from T2NGRA
  INTEGER(ikind),SAVE::&
       StartEqs,EndEqs,&
       StartAxi,EndAxi,&
       StartWal,EndWal

  ! form T2PROF
  REAL(   rkind),SAVE,ALLOCATABLE::&
       GlobalCrd(:,:)               ! d2mfc1

  ! from T2NGRA
  INTEGER(ikind),SAVE,ALLOCATABLE::&
       RowCRS(:),&                  ! i1nidr
       ColCRS(:),&                  ! i1nidc 
       HangedNodeTable(:,:),&       ! i2hbc
       ElementNodeGraph(:,:,:)      ! i3enr
  
  ! from T2VGRA
  LOGICAL,SAVE,ALLOCATABLE::&
       LockEqs(:),&                 ! new
       LockAxi(:),&                 ! new
       LockWal(:),&                 ! new
       HaveMassScaCoef(:,:),&       ! new     
       HaveAdveVecCoef(:,:),&       ! new     
       HaveAdveTenCoef(:,:),&       ! new     
       HaveDiffTenCoef(:,:),&       ! new     
       HaveGradVecCoef(:,:),&       ! new     
       HaveGradTenCoef(:,:),&       ! new     
       HaveExciScaCoef(:,:),&       ! new     
       HaveExciVecCoef(:,:),&       ! new     
       HaveExciTenCoef(:,:),&       ! new     
       HaveSourScaCoef(:,:),&       ! new     
       !
       HaveAdveTenKval(:,:,:),&     ! new     
       HaveGradTenKval(:,:,:),&     ! new     
       HaveExciVecKval(:,:,:),&     ! new     
       HaveExciTenKval(:,:,:,:),&   ! new     
       !
       HaveMat(:,:),&               ! new
       HaveVec(:)                   ! new


  ! from T2CALV
  REAL(   rkind),SAVE,ALLOCATABLE::&
       KnownVar(:,:),&              ! d2ws
       !
       MassScaCoef(:,:,:        ),& ! d3ms
       AdveVecCoef(:,:,:,:      ),& ! d4av
       AdveTenCoef(:,:,:,:,:,:  ),& ! d6at
       DiffTenCoef(:,:,:,:,:    ),& ! d5dt
       GradVecCoef(:,:,:,:      ),& ! d4gv
       GradTenCoef(:,:,:,:,:,:  ),& ! d6gt
       ExciScaCoef(:,:,:        ),& ! d3es
       ExciVecCoef(:,:,:,:,:    ),& ! d5ev
       ExciTenCoef(:,:,:,:,:,:,:),& ! d7et
       SourScaCoef(:,:,:        )   ! d3ss 
    
  PUBLIC T2EXEC_EXECUTE,T2EXEC_TERMINATE

CONTAINS
  
  !C------------------------------------------------------------------
  !C
  !C T2EXEC (PUBLIC)
  !C FEM SOLVER FOR SIMULTANEOUS ADVECTION-DIFFUSION EQUATIONS
  !C 
  !C                2014-05-22 H.SETO
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2EXEC_EXECUTE
    
    !USE T2COMM,ONLY: NNMAX,NEMAX,NKMAX,NVMAX
    !USE T2VGRA,ONLY: HaveMassScaCoef,HaveAdveVecCoef,HaveAdveTenCoef,&
    !     &           HaveDiffTenCoef,HaveGradVecCoef,HaveGradTenCoef,&
    !     &           HaveExciScaCoef,HaveExciVecCoef,HaveExciTenCoef,&
    !     &           HaveSourScaCoef,&
    !     &           HaveAdveTenKval,HaveGradTenKval,HaveExciVecKval,&
    !     &           HaveExciTenKval,&
    !     &           LockEqs,        LockAxi,        LockWal
    
    INTEGER(ikind)::&
         i_v,j_v,i_k,j_k,i_e
    
    REAL(4)::e0time_0,e0time_1
    
    CALL CPU_TIME(e0time_0)
    
    CALL T2EXEC_INITIALIZE

    DO i_e = 1, NEMAX
       
       CALL T2EXEC_SETUP_ELEMENT_VARIABLES(i_e)
       
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

       CALL T2EXEC_STORE
       
    ENDDO
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)')&
         '-- T2EXEC: matrix construction was completed: cpu=', &
         e0time_1-e0time_0,' [s]'
    
    !
    ! set boundary conditions 
    !
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
    
    RETURN
    
  END SUBROUTINE T2EXEC_EXECUTE
  
  !------------------------------------------------------------------ 
  !
  !
  !
  !
  !
  !
  !
  !------------------------------------------------------------------ 
  SUBROUTINE T2EXEC_INITIALIZE

    CALL T2EXEC_ADHOC
    CALL T2EXEC_PUBLIC_ALLOCATE
    CALL T2EXEC_
    RETURN
  END SUBROUTINE T2EXEC_INITIALIZE


  SUBROUTINE T2EXEC_ADHOC

    USE T2COMM,ONLY:i0vmax,i0wmax,i0qmax,i0dmax,i0mmax,i0nmax,&
         &          i0xmax,i0bmax,i0emax,i0hmax,i0lmax,i0amax,&
         &          i0nrmx,&
         ! T2NGRA
         &          i1nidr,i1nidc,i3enr,i2hbc,&
         ! 
         &          d2ws,d3ms,d4av,d6at,d5dt,d4gv,d6gt,d3es,&
         &          d5ev,d7et,d3ss,&
         ! T2PROF
         &          d2mfc1


    NMMAX = i0mmax
    NVMAX = i0vmax
    NKMAX = i0wmax
    NDMAX = i0dmax
    NNMAX = i0nmax
    NMMAX = i0mmax
    NXMAX = i0xmax
    NBMAX = i0bmax
    NEMAX = i0emax
    NHMAX = i0hmax
    NLMAX = i0lmax
    NAMAX = i0amax
    NNRMX = i0nrmx
    
    CALL T2EXEC_ADHOC_ALLOCATE
    
    ! for T2NGRA
    RowCRS               = i1nidr
    ColCRS               = i1nidc
    ElementNodeGraph     = i3enr
    IF(NHMAX.NE.0)&
         HangedNodeTable = i2hbc
    
    ! for T2CALV
    KnownVar    = d2ws
    MassScaCoef = d3ms
    AdveVecCoef = d4av
    AdveTenCoef = d6at
    DiffTenCoef = d5dt
    GradVecCoef = d4gv
    GradTenCoef = d6gt
    ExciScaCoef = d3es
    ExciVecCoef = d5ev
    ExciTenCoef = d7et
    SourScaCoef = d3ss 

    ! for T2PROF
    GlobalCrd = d2mfc1
    RETURN
  
  END SUBROUTINE T2EXEC_ADHOC

  SUBROUTINE T2EXEC_ADHOC_ALLOCATE
    
    INTEGER(ikind),SAVE::&
         NVMAX_save=0,NNMAX_save=0,NDMAX_save=0,NEMAX_save=0,&
         NHMAX_save=0,NMMAX_save=0,NAMAX_save=0,NNRMX_save=0,&
         NKMAX_save=0

    INTEGER(ikind):: ierr
    

    IF(  (NNMAX.NE.NNMAX_save ).OR.&
         (NDMAX.NE.NDMAX_save ).OR.&
         (NMMAX.NE.NMMAX_save ).OR.&
         (NKMAX.NE.NKMAX_save ).OR.&
         (NAMAX.NE.NAMAX_save ).OR.&
         (NEMAX.NE.NEMAX_save ).OR.&
         (NNRMX.NE.NNRMX_save ).OR.&
         (NHMAX.NE.NHMAX_save ).OR.&
         (NVMAX.NE.NVMAX_save ))THEN
       
       CALL T2EXEC_ADHOC_DEALLOCATE
       
       DO
          
          ! for T2NGRA
          ALLOCATE(RowCRS(1:NNRMX), STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ColCRS(1:NAMAX), STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ElementNodeGraph(1:NNMAX,1:4,1:NEMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          IF(NHMAX.NE.0)&
               ALLOCATE(HangedNodeTable(1:2,1:NHMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          !
          ALLOCATE(GlobalCrd(1:NDMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          
          ! T2VGRA
          ALLOCATE(HaveMassScaCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveAdveVecCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveAdveTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveDiffTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradVecCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciScaCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciVecCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciTenCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveSourScaCoef(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          !
          ALLOCATE(HaveAdveTenKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradTenKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciVecKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciTenKval(1:NKMAX,1:NKMAX,1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          !
          ALLOCATE(HaveMat(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveVec(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockEqs(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockAxi(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockWal(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT

          ! T2CALV
          ALLOCATE(KnownVar(   1:NKMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          !
          ALLOCATE(MassScaCoef(1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(AdveVecCoef(1:NDMAX,1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(AdveTenCoef(1:NDMAX,1:NDMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(DiffTenCoef(1:NDMAX,1:NDMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GradVecCoef(1:NDMAX,1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(GradTenCoef(1:NDMAX,1:NDMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ExciScaCoef(1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ExciVecCoef(1:NDMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(ExciTenCoef(1:NDMAX,1:NDMAX,1:NKMAX,1:NKMAX,&
               &               1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(SourScaCoef(1:NVMAX,1:NVMAX,1:NMMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          
          NNMAX_save = NNMAX
          NEMAX_save = NEMAX
          NKMAX_save = NKMAX
          NVMAX_save = NVMAX
          NDMAX_save = NDMAX
          NMMAX_save = NMMAX
          NAMAX_save = NAMAX
          NNRMX_save = NNRMX
          NHMAX_save = NHMAX 
          
          WRITE(6,'(A)') '-- T2EXEC_ADHOC_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')&
            '--T2EXEC_ADHOC_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2EXEC_ADHOC_ALLOCATE

  SUBROUTINE T2EXEC_ADHOC_DEALLOCATE

    ! T2NGRA
    IF(ALLOCATED(RowCRS))           DEALLOCATE(RowCRS)
    IF(ALLOCATED(ColCRS))           DEALLOCATE(ColCRS)
    IF(ALLOCATED(ElementNodeGraph)) DEALLOCATE(ElementNodeGraph)
    IF(ALLOCATED(HangedNodeTable))  DEALLOCATE(HangedNodeTable)

    ! T2PROF
    IF(ALLOCATED(GlobalCrd))        DEALLOCATE(GlobalCrd)
    
    ! T2VGRA
    IF(ALLOCATED(HaveMassScaCoef))  DEALLOCATE(HaveMassScaCoef)
    IF(ALLOCATED(HaveAdveVecCoef))  DEALLOCATE(HaveAdveVecCoef)
    IF(ALLOCATED(HaveAdveTenCoef))  DEALLOCATE(HaveAdveTenCoef)
    IF(ALLOCATED(HaveDiffTenCoef))  DEALLOCATE(HaveDiffTenCoef)
    IF(ALLOCATED(HaveGradVecCoef))  DEALLOCATE(HaveGradVecCoef)
    IF(ALLOCATED(HaveGradTenCoef))  DEALLOCATE(HaveGradTenCoef)
    IF(ALLOCATED(HaveExciScaCoef))  DEALLOCATE(HaveExciScaCoef)
    IF(ALLOCATED(HaveExciVecCoef))  DEALLOCATE(HaveExciVecCoef)
    IF(ALLOCATED(HaveExciTenCoef))  DEALLOCATE(HaveExciTenCoef)
    IF(ALLOCATED(HaveSourScaCoef))  DEALLOCATE(HaveSourScaCoef)
    !
    IF(ALLOCATED(HaveAdveTenKval))  DEALLOCATE(HaveAdveTenKval)
    IF(ALLOCATED(HaveGradTenKval))  DEALLOCATE(HaveGradTenKval)
    IF(ALLOCATED(HaveExciVecKval))  DEALLOCATE(HaveExciVecKval)
    IF(ALLOCATED(HaveExciTenKval))  DEALLOCATE(HaveExciTenKval)
    !
    IF(ALLOCATED(HaveMat))           DEALLOCATE(HaveMat)
    IF(ALLOCATED(HaveVec))           DEALLOCATE(HaveVec)
    IF(ALLOCATED(LockEqs))           DEALLOCATE(LockEqs)
    IF(ALLOCATED(LockAxi))           DEALLOCATE(LockAxi)
    IF(ALLOCATED(LockWal))           DEALLOCATE(LockWal)

    ! T2CALV
    IF(ALLOCATED(KnownVar))         DEALLOCATE(KnownVar)
    !
    IF(ALLOCATED(MassScaCoef))      DEALLOCATE(MassScaCoef)
    IF(ALLOCATED(AdveVecCoef))      DEALLOCATE(AdveVecCoef)
    IF(ALLOCATED(AdveTenCoef))      DEALLOCATE(AdveTenCoef)
    IF(ALLOCATED(DiffTenCoef))      DEALLOCATE(DiffTenCoef)
    IF(ALLOCATED(GradVecCoef))      DEALLOCATE(GradVecCoef)
    IF(ALLOCATED(GradTenCoef))      DEALLOCATE(GradTenCoef)
    IF(ALLOCATED(ExciScaCoef))      DEALLOCATE(ExciScaCoef)
    IF(ALLOCATED(ExciVecCoef))      DEALLOCATE(ExciVecCoef)
    IF(ALLOCATED(ExciTenCoef))      DEALLOCATE(ExciTenCoef)
    IF(ALLOCATED(SourScaCoef))      DEALLOCATE(SourScaCoef)

    RETURN
    
  END SUBROUTINE T2EXEC_ADHOC_DEALLOCATE

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
    
    !USE T2COMM, ONLY:NNMAX,NDMAX,NVMAX,NKMAX
    
    !USE T2CALV, ONLY: KnownVar
    !USE T2NGRA, ONLY: ElementNodeGraph
    
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
    
    DO i_r = 1, 4
       DO i_n = 1, NNMAX
          ElementNodeGraphElem(i_n,i_r)=ElementNodeGraph(i_n,i_r,i_e)
       ENDDO
    ENDDO
    
    
    !
    ! CALCULATE JACOBIAN OF LOCAL COORDINATE 
    ! 
    ! JacDetLocCrd: JACOBIAN OF LOCAL COORDINATES
    ! JacInvLocCrd: INVERSE JACOBI MATRIX OF LOCAL COORDINATES
    ! 
    
    i_n = ElementNodeGraphElem(1,1)
    j_n = ElementNodeGraphElem(2,1)
    k_n = ElementNodeGraphElem(4,1)

    gridSizeRad = GlobalCrd(1,j_n)-GlobalCrd(1,i_n)
    gridSizePol = GlobalCrd(2,k_n)-GlobalCrd(2,i_n)

    gridSizePol  = ABS(gridSizePol)
    gridSizeRad  = ABS(gridSizeRad)

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
    ! KnownVarElem(:,1   ) : B    AT L-TH PICARD ITERATION
    ! KnownVarElem(:,2   ) : R    AT L-TH PICARD ITERATION
    ! KnownVarElem(:,2N+1) : Ub   AT L-TH PICARD ITERATION 
    ! KnownVarElem(:,2N+2) : Qb/P AT L-TH PICARD ITERATION 
    !
    
 
    DO i_n = 1, NNMAX
       i_m = ElementNodeGraphElem(i_n,1)
       DO i_k = 1, NKMAX
          KnownVarElem(i_n,i_k) = KnownVar(i_k,i_m)
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_SETUP_ELEMENT_VARIABLES
  
  !-------------------------------------------------------------------
  !
  !    MASS SCALAR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2EXEC_MS_SUBMATRIX(i_v,j_v)
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX,Dt
    USE T2INTG,ONLY: MassScaIntgPG
    !USE T2CALV,ONLY: MassScaCoef
    
    INTEGER(ikind),INTENT(IN):: i_v,j_v

    INTEGER(ikind)::&
         i_n,j_n,k_n,i_m,i_b,i_x

    REAL(   rkind)::&
         massScaMatElem( 1:NNMAX,1:NNMAX),&
         massScaVecElem( 1:NNMAX        ),&
         massScaCoefElem(1:NNMAX        ),&
         valPriorElem(   1:NNMAX        )
    
    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = ElementNodeGraphElem(i_n,1)
       massScaCoefElem(   i_n            ) &
            = MassScaCoef(    i_v,j_v,i_m) &
            * JacDetLocCrd/Dt
    ENDDO
    
    SELECT CASE(j_v)
       
    CASE (1:3)   ! for FSA variables 
       DO i_n = 1, NNMAX
          i_b = ElementNodeGraphElem(i_n,4)
          valPriorElem(i_n) = Xvec(i_b,j_v)
       ENDDO
    CASE DEFAULT ! for 2D dependent variables
       DO i_n = 1, NNMAX
          i_x = ElementNodeGraphElem(i_n,2)
          valPriorElem(i_n) = Xvec(i_x,j_v)          
       ENDDO
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

    massScaVecElem(1:NNMAX) = 0.D0
    DO j_n = 1, NNMAX  
    DO i_n = 1, NNMAX
       massScaVecElem(       i_n    ) &
            = massScaVecElem(i_n    ) &
            + massScaMatElem(i_n,j_n) &
            * valPriorElem(      j_n)
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: AdveVecIntgPG
    !USE T2CALV,ONLY: AdveVecCoef

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
       i_m = ElementNodeGraphElem(i_n,1)
       DO i_d = 1, NDMAX
          adveVecCoefElem(   i_d,i_n            ) &
               = AdveVecCoef(i_d,    i_v,j_v,i_m) &
               * JacDetLocCrd
          
       ENDDO
    ENDDO

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

    ! STORE SUBMATRIX
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
            + adveVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: AdveTenIntgPG
    !USE T2COEF,ONLY: AdveTenCoef

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
       adveTenKvarElem(i_n) = KnownVarElem(i_n,i_k)
    ENDDO
    
    DO i_n = 1, NNMAX
       i_m = ElementNodeGraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          adveTenCoefElem(   i_d,j_d,i_n                ) &
               = AdveTenCoef(i_d,j_d,    i_k,i_v,j_v,i_m) &
               * JacDetLocCrd
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
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
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
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: DiffTenIntgPG
    !USE T2COEF,ONLY: DiffTenCoef
    
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
       i_m = ElementNodeGraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          diffTenCoefElem(   i_d,j_d,i_n            ) &
               = diffTenCoef(i_d,j_d,    i_v,j_v,i_m) &
               * JacDetLocCrd
       ENDDO
       ENDDO
    ENDDO

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
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: GradVecIntgPG
    !USE T2CALV,ONLY: GradVecCoef
    
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
       i_m = ElementNodeGraphElem(i_n,1)
       DO i_d = 1, NDMAX
          gradVecCoefElem(   i_d,i_n            ) &
               = GradVecCoef(i_d,    i_v,j_v,i_m) &
               * JacDetLocCrd
       ENDDO
    ENDDO
    
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

    ! main loop
    
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
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
            + gradVecMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GV_SUBMATRIX
  
  !------------------------------------------------------------------
  !
  !   GRADIENT TENSOR SUBMATRCES (PRIVATE)
  !
  !                     2014-05-22 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GT_SUBMATRIX(i_v,j_v,i_k)
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: GradTenIntgPG
    !USE T2CALV,ONLY: GradTenCoef
    
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
       i_m = ElementNodeGraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          gradTenCoefElem(   i_d,j_d,i_n                ) &
               = GradTenCoef(i_d,j_d,    i_k,i_v,j_v,i_m) &
               * JacDetLocCrd
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
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    
    DO i_n = 1, NNMAX
       gradTenKvarElem(i_n) = KnownVarElem(i_n,i_k)
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
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
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

    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: ExciScaIntgPG
    !USE T2CALV,ONLY: ExciScaCoef
         
    INTEGER(ikind),INTENT(IN)::i_v,j_v
    INTEGER(ikind)::&
         i_n,j_n,k_n,&
         i_m


    REAL(   rkind)::&
         exciScaMatElem( 1:NNMAX,1:NNMAX),&
         exciScaCoefElem(1:NNMAX        )

    ! initialization
        
    DO i_n = 1,NNMAX
       i_m = ElementNodeGraphElem(i_n,1)
       exciScaCoefElem(   i_n            ) &
            = ExciScaCoef(    i_v,j_v,i_m) &
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
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: ExciVecIntgPG
    !USE T2CALV,ONLY: ExciVecCoef
    
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
       i_m = ElementNodeGraphElem(i_n,1)
       DO i_d = 1, NDMAX
          exciVecCoefElem(   i_d,i_n                ) &
               = ExciVecCoef(i_d,    i_k,i_v,j_v,i_m) &
               * JacDetLocCrd
       ENDDO
    ENDDO
    
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
        
    DO i_n = 1, NNMAX
       exciVecKvarElem(i_n) = KnownVarElem(i_n,i_k)
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
       AmatElem(             i_n,j_n,i_v,j_v) &
            = AmatElem(      i_n,j_n,i_v,j_v) &
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: ExciTenIntgPG
    !USE T2CALV,ONLY: ExciTenCoef
    
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
       i_m = ElementNodeGraphElem(i_n,1)
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          exciTenCoefElem(   i_d,j_d,i_n)&
               = ExciTenCoef(i_d,j_d,    i_k,j_k,i_v,j_v,i_m) &
               * JacDetLocCrd
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
                  * JacInvLocCrd(   i_d,    k_d        ) &
                  * JacInvLocCrd(       j_d,    l_d    )
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
       AmatElem(             i_n,j_n,I_v,J_v) &
            = AmatElem(      i_n,j_n,I_v,J_v) &
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
    
    !USE T2COMM,ONLY: NNMAX,NDMAX,NVMAX
    USE T2INTG,ONLY: SourScaIntgPG
    !USE T2CALV,ONLY: SourScaCoef

    INTEGER(ikind),INTENT(IN)::i_v,j_v
    INTEGER(ikind)::&
         i_n,j_n,&
         i_m
    
    REAL(   rkind)::&
         sourScaVecElem( 1:NNMAX),&
         sourScaCoefElem(1:NNMAX)

    ! initialization
    
    DO i_n = 1, NNMAX
       i_m = ElementNodeGraphElem(i_n,1)
       sourScaCoefElem(   i_n            ) &
            = SourScaCoef(    i_v,j_v,i_m) &
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
       BvecElem(             i_n,i_v) &
            = BvecElem(      i_n,i_v) & 
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
    
    !USE T2COMM,ONLY:NNMAX,NVMAX
    !USE T2NGRA,ONLY: RowCRS,ColCRS,HangedNodeTable
    !USE T2VGRA,ONLY: HaveMat
    INTEGER(ikind)::&
         i_n,j_n,&
         i_v,j_v,&
         nodeTableA(1:NNMAX),& ! 1D node table
         nodeTableB(1:NNMAX),& ! 1D node table for FSA
         nodeTableC(1:NNMAX)   ! 2D node table
    
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
    
    DO i_n = 1, NNMAX
       nodeTableA(i_n) = ElementNodeGraphElem(i_n,4)
       nodeTableB(i_n) = ElementNodeGraphElem(i_n,3)
       nodeTableC(i_n) = ElementNodeGraphElem(i_n,2)
    ENDDO
    
    !
    ! for stiffness matrix
    !

    ! i_v: 1D value (FSA)
    ! j_v: 1D value (FSA)
    CALL T2EXEC_STORE_MATRIX(1    ,3    ,nodeTableA,& ! row
         &                   1    ,3    ,nodeTableA)  ! column
            
    ! i_v: 1D value (FSA)
    ! j_v: 2D value
    
    CALL T2EXEC_STORE_MATRIX(1    ,3    ,nodeTableB,& ! row
         &                   4    ,NVMAX,nodeTableC)  ! column

    ! i_v: 2D value 
    ! j_v: 1D value (FSA) 
    
    CALL T2EXEC_STORE_MATRIX(4    ,NVMAX,nodeTableC,& ! row
         &                   1    ,3    ,nodeTableA)  ! column

    ! i_v: 2D value
    ! j_v: 2D value
    
    CALL T2EXEC_STORE_MATRIX(4    ,NVMAX,nodeTableC,& ! row
         &                   4    ,NVMAX,nodeTableC)  ! column

    !
    ! for RHS vector
    !
    
    ! i_v: 1D value (FSA)
    CALL T2EXEC_STORE_VECTOR(1    ,3    ,nodeTableA) ! row
    
    ! i_v 2D value
    CALL T2EXEC_STORE_VECTOR(4    ,NVMAX,nodeTableC) ! row
    
    RETURN
    
  END SUBROUTINE T2EXEC_STORE

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
    
    INTEGER,INTENT(IN)::&
         rowValStart,rowValEnd,rowNodeTable(1:NNMAX),&
         colValStart,colValEnd,colNodeTable(1:NNMAX)

    !USE T2COMM,ONLY: NNMAX,NBMAX
    !USE T2NGRA,ONLY: RowCRS, ColCRS,HangedNodeTable
    !USE T2VGRA,ONLY:HaveMat
    
    INTEGER(ikind)::&
         i_rowN,j_rowN,k_rowN,&
         i_colN,j_colN,k_colN,l_colN

    INTEGER(ikind)::&
         i_n,j_n,i_v,j_v,i_a

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       
       i_rowN = rowNodeTable(i_n)
       i_colN = colNodeTable(j_n)
       
       IF((i_rowN.LE.NBMAX).AND.(i_colN.LE.NBMAX))THEN
          
          DO i_a = RowCRS(i_rowN), RowCRS(i_rowN+1)-1
             j_colN = ColCRS(i_a)
             IF(i_colN.EQ.j_colN)THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      Amat(             i_v,j_v,        i_a) &
                           = Amat(      i_v,j_v,        i_a) &
                           + AmatElemTF(i_v,j_v,i_n,j_n    )
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
       ELSEIF((i_rowN.LE.NBMAX).AND.(i_colN.GT.NBMAX))THEN
          
          i_colN = i_colN - NBMAX
          j_colN = HangedNodeTable(1,i_colN)
          k_colN = HangedNodeTable(2,i_colN)
          
          DO i_a = RowCRS(i_rowN), RowCRS(i_rowN+1)-1
             l_colN = ColCRS(i_a)
             IF((l_colN.EQ.j_colN).OR.(l_colN.EQ.k_colN))THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      Amat(             i_v,j_v,        i_a) &
                           = Amat(      i_v,j_v,        i_a) &
                           + AmatElemTF(i_v,j_v,i_n,j_n    ) &
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
          
          DO i_a = RowCRS(j_rowN), RowCRS(j_rowN+1)-1
             j_colN = ColCRS(i_a) 
             IF(j_colN.EQ.i_colN)THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      Amat(             i_v,j_v,        i_a) &
                           = Amat(      i_v,j_v,        i_a) &
                           + AmatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.50D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
          DO i_a = RowCRS(k_rowN), RowCRS(k_rowN+1)-1
             j_colN = ColCRS(i_a) 
             IF(j_colN.EQ.i_colN)THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      Amat(             i_v,j_v,        i_a) &
                           = Amat(      i_v,j_v,        i_a) &
                           + AmatElemTF(i_v,j_v,i_n,j_n    ) &
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
          
          DO i_a = RowCRS(j_rowN), RowCRS(j_rowN+1)-1
             l_colN = ColCRS(i_a)
             IF((l_colN.EQ.j_colN).OR.(l_colN.EQ.k_colN))THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      Amat(             i_v,j_v,        i_a) &
                           = Amat(      i_v,j_v,        i_a) &
                           + AmatElemTF(i_v,j_v,i_n,j_n    ) &
                           * 0.25D0
                   ENDIF
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i_a = RowCRS(k_rowN), RowCRS(k_rowN+1)-1
             l_colN = ColCRS(i_a)
             IF((l_colN.EQ.j_colN).OR.(l_colN.EQ.k_colN))THEN
                DO j_v = colValStart, colValEnd
                DO i_v = rowValStart, rowValEnd
                   IF(HaveMat(i_v,j_v))THEN
                      Amat(             i_v,j_v,        i_a) &
                           = Amat(      i_v,j_v,        i_a) &
                           + AmatElemTF(i_v,j_v,i_n,j_n    ) &
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
    
    INTEGER(ikind),INTENT(IN)::&
         rowValStart,rowValEnd,rowNodeTable(1:NNMAX)
    
    INTEGER(ikind)::&
         i_n,i_v,&
         i_row,j_row,k_row
    
    DO i_n = 1, NNMAX
       i_row = rowNodeTable(i_n)
       IF(    i_row.LE.NBMAX)THEN
          DO i_v = rowValStart,rowValEnd
             Bvec(              i_v,    i_row) &
                  = Bvec(       i_v,    i_row) &
                  + bvecElemTF( i_v,i_n      )
          ENDDO
       ELSEIF(i_row.GT.NBMAX)THEN
          i_row = i_row - NBMAX
          j_row = HangedNodeTable(1,i_row)
          k_row = HangedNodeTable(2,i_row)
          DO i_v = rowValStart, rowValEnd
             Bvec(             i_v,    j_row) &
                  = Bvec(      i_v,    j_row) &
                  + BvecElemTF(i_v,i_n    ) &
                  * 0.50D0
          ENDDO
          DO i_v = rowValStart, rowValEnd
             Bvec(             i_v,    k_row) &
                  = Bvec(      i_v,    k_row) &
                  + BvecElemTF(i_v,i_n      ) &
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
  SUBROUTINE T2EXEC_LOCK_VALUES(nodeStart,nodeEnd,varLockTable)
    
    !USE T2COMM,ONLY: NVMAX
    !USE T2NGRA,ONLY: RowCRS, colCRS

    INTEGER(ikind),INTENT(IN)::nodeStart,nodeEnd
    LOGICAL,INTENT(IN)::varLockTable(1:NVMAX)
    
   
    INTEGER(ikind)::&
         i_a,i_v,j_v,i_row,i_col
    
    DO i_row = nodeStart, nodeEnd
       
       ! stiffness matrix     
       DO i_a = RowCRS(i_row), RowCRS(i_row+1)-1
          i_col = ColCRS(i_a)
          DO j_v = 1, NVMAX
          DO i_v = 1, NVMAX
             IF(varLockTable(i_v))THEN
                IF((i_row.EQ.i_col).AND.(i_v.EQ.j_v))THEN
                   Amat(i_v,j_v,i_a) = 1.D0
                ELSE
                   Amat(i_v,j_v,i_a) = 0.D0
                ENDIF
             ENDIF
          ENDDO
          ENDDO
       ENDDO
       
       ! RHS vector 
       DO i_v = 1, NVMAX
          IF(varLockTable(i_v))THEN
             Bvec(i_v,i_row) = Xvec(i_v,i_row)
          ENDIF
       ENDDO
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_LOCK_VALUES
  
  !
  !
  !
  !
  !
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
    
    
    INTEGER(ikind)::i_v,j_v

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

  SUBROUTINE T2EXEC_TERMINATE
    RETURN
  END SUBROUTINE T2EXEC_TERMINATE
END MODULE T2EXEC
