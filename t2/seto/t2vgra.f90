!------------------------------------------------------------------
!
!   MODULE FOR GENERATING VARIABLE-VARIABLE GRAPH 
!       
!
!                       LAST UPDATE 2014-05-23 H.Seto
!
!   T2INTG requires following variables:
!
!   [from T2CNST]
!        rkind,ikind
!   [from T2COMM]
!        NSMAX,NVMAX,NKMAX
!
!   T2INTG provides following variables:
!
!        HaveMassScaCoef
!        HaveAdveVecCoef
!        HaveAdveTenCoef
!        HaveDiffTenCoef
!        HaveGradVecCoef
!        HaveGradTenCoef
!        HaveExciScaCoef
!        HaveExciVecCoef
!        HaveExciTenCoef
!        HaveSourScaCoef
!        HaveAdveTenKval
!        HaveGradTenKval
!        HaveExciVecKval
!        HaveExciTenKval
!        HaveMat
!        LockEqs
!        LockAxi
!        LockWal
!
!   and subroutines:
!
!        T2VGRA_EXCUTE
!        T2VGRA_TERMINATE
!
! -------------------------------------------------------------------
MODULE T2VGRA
  
  USE T2CNST,ONLY:&
       ikind,rkind
  
  IMPLICIT NONE
  
  PRIVATE 

  ! >>>> this section will be removed in final version  >>>>>
  INTEGER(ikind)::NVMAX,NKMAX,NSMAX
  LOGICAL::UsePotentialDescription = .FALSE.
  LOGICAL::&
       SolveField,&
       SolveElectron,&
       SolveIons,&
       SolveDensity,&
       SolveFlux,&
       SolvePressure,&
       SolveHeatFlux,&
       !
       LockPoloidalMageticFieldOnAxis,&
       LockToroidalMageticFieldOnAxis,&
       LockRadialElectricFieldOnAxis,&
       LockPoloidalElectricFieldOnAxis,&
       LockToroidalElectricFieldOnAxis,&
       LockDensityOnAxis,&
       LockRaidalFluxOnAxis,&
       LockParallelFluxOnAxis,&
       LockToroidalFluxOnAxis,&
       LockPoroidalFluxOnAxis,& 
       LockPressureOnAxis,& 
       LockRaidalHeatFluxOnAxis,& 
       LockParallelHeatFluxOnAxis,&
       LockToroidalHeatFluxOnAxis,&
       LockPoroidalHeatFluxOnAxis,&
       !
       LockPoloidalMageticFieldOnWall,&
       LockToroidalMageticFieldOnWall,&
       LockRadialElectricFieldOnWall,&
       LockPoloidalElectricFieldOnWall,&
       LockToroidalElectricFieldOnWall,&
       LockDensityOnWall,&
       LockRaidalFluxOnWall,&
       LockParallelFluxOnWall,&
       LockToroidalFluxOnWall,&
       LockPoroidalFluxOnWall,& 
       LockPressureOnWall,& 
       LockRaidalHeatFluxOnWall,& 
       LockParallelHeatFluxOnWall,&
       LockToroidalHeatFluxOnWall,&
       LockPoroidalHeatFluxOnWall

  LOGICAL,SAVE,ALLOCATABLE::&
       HaveMassScaCoef(:,:    ),HaveAdveVecCoef(:,:    ),&
       HaveAdveTenCoef(:,:    ),HaveDiffTenCoef(:,:    ),&
       HaveGradVecCoef(:,:    ),HaveGradTenCoef(:,:    ),&
       HaveExciScaCoef(:,:    ),HaveExciVecCoef(:,:    ),&
       HaveExciTenCoef(:,:    ),HaveSourScaCoef(:,:    ),&
       !
       HaveAdveTenKval(:,:,:  ),HaveGradTenKval(:,:,:  ),&
       HaveExciVecKval(:,:,:  ),HaveExciTenKval(:,:,:,:),&
       !
       HaveMat(:,:), LockEqs(:),  LockAxi(:),  LockWal(:)
  ! <<<< this section will be removed in final version  <<<<<

  PUBLIC T2VGRA_EXECUTE
  
CONTAINS 
  
  SUBROUTINE T2VGRA_EXECUTE
    
    USE T2COMM, ONLY: i0smax,i0vmax,i0wmax,idfile
    
    ! >>>> this section will be removed in final version  >>>>>
    NSMAX = i0smax
    NVMAX = i0vmax
    NKMAX = i0wmax
    ! <<<< this section will be removed in final version  <<<<<
    
    CALL T2VGRA_ALLOCATE
    
    IF(.NOT.UsePotentialDescription)THEN 
       
       CALL T2VGRA_MS_COEF_EB
       CALL T2VGRA_AV_COEF_EB
       CALL T2VGRA_AT_COEF_EB
       CALL T2VGRA_DT_COEF_EB
       CALL T2VGRA_GV_COEF_EB
       CALL T2VGRA_GT_COEF_EB
       CALL T2VGRA_ES_COEF_EB
       CALL T2VGRA_EV_COEF_EB
       CALL T2VGRA_ET_COEF_EB
       CALL T2VGRA_SS_COEF_EB

       CALL T2VGRA_LOCK_EQUA_EB
       CALL T2VGRA_LOCK_AXIS_EB
       CALL T2VGRA_LOCK_WALL_EB

    ELSE
       
       CALL T2VGRA_MS_COEF_PhiA
       CALL T2VGRA_AV_COEF_PhiA
       CALL T2VGRA_AT_COEF_PhiA
       CALL T2VGRA_DT_COEF_PhiA
       CALL T2VGRA_GV_COEF_PhiA
       CALL T2VGRA_GT_COEF_PhiA
       CALL T2VGRA_ES_COEF_PhiA
       CALL T2VGRA_EV_COEF_PhiA
       CALL T2VGRA_ET_COEF_PhiA
       CALL T2VGRA_SS_COEF_PhiA

       CALL T2VGRA_LOCK_EQUA_PhiA
       CALL T2VGRA_LOCK_AXIS_PhiA
       CALL T2VGRA_LOCK_WALL_PhiA
              
    ENDIF

    CALL T2VGRA_VAR_MAT
    
    IF(idfile.ge.5) CALL T2VGRA_OUTPUT
    
    RETURN
    
  END SUBROUTINE T2VGRA_EXECUTE
  
  SUBROUTINE T2VGRA_ALLOCATE
    
    INTEGER(ikind),SAVE::&
         NVMAX_save=0,NKMAX_save=0
    
    INTEGER(ikind):: ierr
    
    IF(  (NVMAX.NE.NVMAX_save).OR.&
         (NKMAX.NE.NKMAX_save))THEN
       
       CALL T2VGRA_DEALLOCATE
       
       DO
          
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
    
          ALLOCATE(HaveAdveTenKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveGradTenKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciVecKval(1:NKMAX,        1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(HaveExciTenKval(1:NKMAX,1:NKMAX,1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
    
          ALLOCATE(HaveMat(1:NVMAX,1:NVMAX),&
               &                    STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockEqs(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockAxi(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          ALLOCATE(LockWal(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
          
          NVMAX_save = NVMAX
          NKMAX_save = NKMAX 
                 
          WRITE(6,'(A)') '-- T2VGRA_ALLOCATE: completed'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)')&
            '--T2VGRA_ALLOCATE: ALLOCATION ERROR: ECODE=',ierr
       STOP
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2VGRA_ALLOCATE
  
  SUBROUTINE T2VGRA_DEALLOCATE
    
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
    
    IF(ALLOCATED(HaveAdveTenKval))  DEALLOCATE(HaveAdveTenKval)
    IF(ALLOCATED(HaveGradTenKval))  DEALLOCATE(HaveGradTenKval)
    IF(ALLOCATED(HaveExciVecKval))  DEALLOCATE(HaveExciVecKval)
    IF(ALLOCATED(HaveExciTenKval))  DEALLOCATE(HaveExciTenKval)
    
    IF(ALLOCATED(HaveMat))          DEALLOCATE(HaveMat)
    IF(ALLOCATED(LockEqs))          DEALLOCATE(LockEqs)
    IF(ALLOCATED(LockAxi))          DEALLOCATE(LockAxi)
    IF(ALLOCATED(LockWal))          DEALLOCATE(LockWal)
    
    RETURN
    
  END SUBROUTINE T2VGRA_DEALLOCATE
    
  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF MASS SCALAR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveMassScaCoef(i_v,j_v) = .TRUE.", 
  !                 the term "(d/dt)(M_{i_v,j_v}*f_{j_v})" 
  !                           is exist in equation system.
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_MS_COEF_EB
    
    !USE T2COMM,ONLY:NSMAX,NVMAX,HaveMassScaCoef
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
        
    !C initialization
    HaveMassScaCoef(1:NVMAX,1:NVMAX) = .FALSE.

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
 
    i_v = 1                   ! Equation for psi'
    j_v = 1;                  HaveMassScaCoef(i_v,j_v) = .TRUE.
    
    i_v = 2                   ! Equation for I
    j_v = 2;                  HaveMassScaCoef(i_v,j_v) = .TRUE.
 
    i_v = 3                   ! Equation for E_{\zeta}
    j_v = 3;                  HaveMassScaCoef(i_v,j_v) = .TRUE.
    
    i_v = 4                   ! Equation for E_{\chi}
    j_v = 4;                  HaveMassScaCoef(i_v,j_v) = .TRUE.
    
    i_v = 5                   ! Equation for E_{\rho}
    
    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !
    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s  
       vOffsetA = vOffsetB

       i_v =  6 + vOffsetA    ! Equation for n_{a}
       j_v =  6 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
       
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v = 13 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v = 14 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_MS_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF ADVECTION VECTOR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveAdveVecCoef(i_v,j_v) = .TRUE.", 
  !                 the term "div (\vec{V}_{i_v,j_v}*f_{j_v})" 
  !                           is exist in equation system.
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AV_COEF_EB
    
    !USE T2COMM,ONLY: NSMAX,NVMAX,HaveAdveVecCoef
    
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    HaveAdveVecCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    j_v = 1;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    
    i_v = 2                   ! Equation for I
    j_v = 2;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    
    i_v = 3                   ! Equation for E_{\zeta}
    j_v = 1;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    j_v = 3;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.

    i_v = 4                   ! Equation for E_{\chi}
    j_v = 4;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    
    i_v = 5                   ! Equation for E_{\rho}
    j_v = 4;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    j_v = 5;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.

    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !

    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       vOffsetB = vOffsetA
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
 
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.

       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       !C
       !C EQUATION FOR Qb
       !C
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  8 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       !C
       !C EQUATION FOR Qt
       !C
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  9 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_AV_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF ADVECTION TENSOR TERMS 
  ! 
  !                    EB VERSION
  !
  !       IF "HaveAdveTenCoef(iv,jv) = .TRUE.", 
  !       the term "div (vec{V}_{iv,jv}*f_{jv})"
  !                           is exist in equation system.
  !       where
  !       vec{V}_{iv,jv} = ten{V}_{iv,jv,ik} dot grad g_{ik} 
  !
  !       IF "HaveAdveTenKval(ik,iv,jv) = .TRUE.", 
  !       the physical quantitiy "g_{ik} " 
  !                           is exist in equation system.
  ! 
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AT_COEF_EB
    
    !USE T2COMM,ONLY: NSMAX, NVMAX, NKMAX,&
    !     &           HaveAdveTenCoef, HaveAdveTenKval
    
    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveAdveTenCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveAdveTenKval(1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I
    i_v = 3                   ! Equation for E_{\zeta}
    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}
    
    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       vOffsetB = vOffsetA
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.

       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.    
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.   
       
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_AT_COEF_EB

  
  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF DIFFUSION TENSOR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveDiffTenCoef(iv,jv) = .TRUE.", 
  !
  !       the term "div (ten{D}_{iv,jv}cdot grad f_{jv})"
  !
  !                           is exist in equation system.
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_DT_COEF_EB
    
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    
    HaveDiffTenCoef(1:NVMAX,1:NVMAX) = .FALSE.

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I
    i_v = 3                   ! Equation for E_{\zeta}
    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}

    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !
    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       vOffsetB = vOffsetA

       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
    
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_DT_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF GRADINET VECTOR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveGradVecCoef(iv,jv) = .TRUE.", 
  !
  !       the term "vec{A}_{iv,jv} dot grad f_{jv}"
  !
  !                           is exist in equation system.
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GV_COEF_EB
    
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    
    HaveGradVecCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !

    i_v =  1                  ! Equation for psi'
    j_v =  3;                 HaveGradVecCoef(i_v,j_v) = .TRUE.

    i_v =  2                  ! Equation for I
    j_v =  4;                 HaveGradVecCoef(i_v,j_v) = .TRUE.
    j_v =  5;                 HaveGradVecCoef(i_v,j_v) = .TRUE.

    i_v =  3                  ! Equation for E_{\zeta}

    i_v =  4                  ! Equation for E_{\chi}
    j_v =  2;                 HaveGradVecCoef(i_v,j_v) = .TRUE.
    
    i_v =  5                  ! Equation for E_{\rho}    

    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       vOffsetB = vOffsetA
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       j_v = 11 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
       j_v =  6;              HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11;              HaveGradVecCoef(i_v,j_v) = .TRUE.
       ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<

       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
       j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}       
       ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
       j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
       
       i_v = 15 + vOffsetA    ! Equation for Q_{a\zeta}
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GV_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF GRADINET TENSOR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveGradTenCoef(iv,jv) = .TRUE.", 
  !
  !       the term "vec{A}_{iv,jv} cdot grad f_{jv}"
  !
  !                           is exist in equation system.
  !       where 
  !            vec{A}_{iv,jv} = grad g_{ik} dot ten{A}_{iv,jv,ik}
  !
  !       IF "HaveGradTenKval(ik,iv,jv) = .TRUE.", 
  !       the physical quantitiy "g_{ik} " 
  !                           is exist in equation system.
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GT_COEF_EB
    
    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,kOffsetX,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveGradTenCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveGradTenKval(1:NVMAX,1:NVMAX,1:NVMAX) = .FALSE.

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I
    i_v = 3                   ! Equation for E_{\zeta}
    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}

    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       vOffsetB = vOffsetA 
       kOffsetX =  2*i_s

       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}

       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GT_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF EXCITATION SCALAR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveExciScaCoef(iv,jv) = .TRUE.", 
  !
  !       the term "{C}_{iv,jv} * f_{jv}"
  !
  !                           is exist in equation system.
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ES_COEF_EB
    
    INTEGER(ikind)::&
         i_s,i_v,vOffsetA,&
         j_s,j_v,vOffsetB
    
    !initialization
    
    HaveExciScaCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'

    i_v = 2                   ! Equation for I
    j_v = 4;                  HaveExciScaCoef(i_v,j_v) = .TRUE.
    
    i_v = 3                   ! Equation for E_{\zeta}
    DO j_s = 0, NSMAX-1
       vOffsetB = 10*j_s
       j_v =  9 + vOffsetB;   HaveExciScaCoef(i_v,j_v) = .TRUE.
    ENDDO
    
    i_v = 4                   ! Equation for E_{\chi}
    DO j_s = 0, NSMAX-1
       vOffsetB = 10*j_s
       j_v =  7 + vOffsetB;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetB;   HaveExciScaCoef(i_v,j_v) = .TRUE.
    ENDDO

    i_v = 5                   ! Equation for E_{\rho}
    DO j_s = 0, NSMAX-1
       vOffsetB = 10*j_s
       j_v =  6 + vOffsetB;   HaveExciScaCoef(i_v,j_v) = .TRUE.
    ENDDO
    
    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  5;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       DO j_s = 0, NSMAX-1
          vOffsetB = 10*j_s
          j_v =  8 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v = 13 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
       ENDDO

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  7 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
       j_v =  8;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  9;              HaveExciScaCoef(i_v,j_v) = .TRUE.       
       !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
       DO j_s = 0, NSMAX-1
          vOffsetB = 10*j_s
          j_v =  9 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v = 14 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
       ENDDO
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.

       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
       j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  5;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       DO j_s = 0, NSMAX-1
          vOffsetB = 10*j_s
          j_v =  8 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v = 13 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
       ENDDO

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 12 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
       j_v = 13;              HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 14;              HaveExciScaCoef(i_v,j_v) = .TRUE.       
       !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
       DO j_s = 0, NSMAX-1
          vOffsetB = 10*j_s
          j_v =  9 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v = 14 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
       ENDDO
       
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
       j_v = 13 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_ES_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF EXCITATION VECTOR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveExciVecCoef(iv,jv) = .TRUE.", 
  !
  !       the term "{C}_{iv,jv} * f_{jv}"
  !
  !                           is exist in equation system.
  !       where
  !        
  !       C_{iv,jv} = grad f_{ik} dot vec{C}_{iv,jv,ik} 
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EV_COEF_EB
    
    INTEGER(ikind)::&
         i_s,i_v,vOffsetA,&
         j_v
    
    ! initialization
    HaveExciVecCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I

    i_v = 3                   ! Equation for E_{\zeta}
    j_v = 1;                  HaveExciVecCoef(i_v,j_v) = .TRUE.

    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}
    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}

       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetA;   HaveExciVecCoef(i_v,j_v) = .TRUE.

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  3;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       j_v =  4;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveExciVecCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveExciVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  3;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       j_v =  4;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_EV_COEF_EB

  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF EXCITATION TENSOR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveExciTenCoef(iv,jv) = .TRUE.", 
  !
  !       the term "{C}_{iv,jv} * f_{jv}"
  !
  !                           is exist in equation system.
  !       where
  !        
  !       C_{iv,jv} = grad g_{ik} dot vec{C}_{iv,jv,ik,jk} 
  !                                                dot grad g_{jk} 
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ET_COEF_EB
    
    INTEGER(ikind)::&
         i_s,i_v,j_v,vOffsetA
   
    ! initialization
    HaveExciTenCoef(1:NVMAX,1:NVMAX) = .FALSE.

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I
    i_v = 3                   ! Equation for E_{\zeta}
    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}

    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !    
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}

       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  9 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  9 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.

       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  9 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
    
    ENDDO
    
    RETURN 

  END SUBROUTINE T2VGRA_ET_COEF_EB
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR SOURCE SCALAR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_SS_COEF_EB
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2ssvt
    
    INTEGER(ikind)::&
         i0sidi,i0vidi,i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2ssvt(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2ssvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2ssvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vidi = 3

    !C Et
    i0vidj = 3
    i2ssvt(i0vidi,i0vidj) = 1
        
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    !C Ep
    i0vidj = 4
    i2ssvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Er
    !C

    i0vidi = 5
    
    !C Ep
    i0vidj = 5
    i2ssvt(i0vidi,i0vidj) = 1
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       !C
       !C EQUATION FOR N
       !C

       i0vidi = 10*i0sidi - 4
              
       !C N
       i0vidj = 10*i0sidi - 4
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi = 10*i0sidi - 3
       
       !C Qb
       i0vidj = 10*i0sidi - 3
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 10*i0sidi - 2
       
       !C Fb
       i0vidj = 10*i0sidi - 2
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = 10*i0sidi - 1
       
       !C Ft
       i0vidj = 10*i0sidi - 1
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Fp
       !C

       i0vidi = 10*i0sidi
       
       !C Fp
       i0vidj = 10*i0sidi 
       i2ssvt(i0vidi,i0vidj) = 1

       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 10*i0sidi + 1
       
       !C P
       i0vidj = 10*i0sidi + 1
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 10*i0sidi + 2
       
       !C Qr
       i0vidj = 10*i0sidi + 2
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 10*i0sidi + 3
       
       !C Qb
       i0vidj = 10*i0sidi + 3
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 10*i0sidi + 4

       !C Qt
       i0vidj = 10*i0sidi + 4
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qp
       !C
       
       i0vidi = 10*i0sidi + 5

       !C Qp
       i0vidj = 10*i0sidi + 5
       i2ssvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN

  END SUBROUTINE T2VGRA_SS_COEF_EB

    
  !C------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR VALIABLE MATRIX ARRAY
  !C 
  !C          MODIFIED 2014-03-06
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2VGRA_VAR_MAT
    
    !USE T2COMM,ONLY:NVMAX
   
    INTEGER(ikind)::i_v
    
    HaveMat =  HaveMassScaCoef.OR. &
         &     HaveAdveVecCoef.OR. &
         &     HaveAdveTenCoef.OR. &
         &     HaveDiffTenCoef.OR. &
         &     HaveGradVecCoef.OR. &
         &     HaveGradTenCoef.OR. &
         &     HaveExciScaCoef.OR. &
         &     HaveExciVecCoef.OR. &
         &     HaveExciTenCoef
    
    ! for eliminate singularity in variable matrix 
    DO i_v = 1, NVMAX
       HaveMat(i_v,i_v) = .TRUE.
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_VAR_MAT

  !
  !
  !
  SUBROUTINE T2VGRA_LOCK_EQUA_EB
    
    INTEGER(ikind)::i_s,vOffsetA
    
    ! initialization
    LockEqs(1:NVMAX) = .TRUE.
    
    IF(SolveField) LockEqs(1:5) = .FALSE.
    
    IF(SolveElectron)THEN
       LockEqs( 6) = .NOT.SolveDensity
       LockEqs( 7) = .NOT.SolveFlux
       LockEqs( 8) = .NOT.SolveFlux
       LockEqs( 9) = .NOT.SolveFlux
       LockEqs(10) = .NOT.SolveFlux
       LockEqs(11) = .NOT.SolvePressure
       LockEqs(12) = .NOT.SolveHeatFlux
       LockEqs(13) = .NOT.SolveHeatFlux
       LockEqs(14) = .NOT.SolveHeatFlux
       LockEqs(15) = .NOT.SolveHeatFlux
    ENDIF
    
    IF(SolveIons)THEN
       DO i_s = 1, NSMAX-1
          vOffsetA = 10*i_s
          LockEqs( 6+vOffsetA) = .NOT.SolveDensity
          LockEqs( 7+vOffsetA) = .NOT.SolveFlux
          LockEqs( 8+vOffsetA) = .NOT.SolveFlux
          LockEqs( 9+vOffsetA) = .NOT.SolveFlux
          LockEqs(10+vOffsetA) = .NOT.SolveFlux
          LockEqs(11+vOffsetA) = .NOT.SolvePressure
          LockEqs(12+vOffsetA) = .NOT.SolveHeatFlux
          LockEqs(13+vOffsetA) = .NOT.SolveHeatFlux
          LockEqs(14+vOffsetA) = .NOT.SolveHeatFlux
          LockEqs(15+vOffsetA) = .NOT.SolveHeatFlux
       ENDDO
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2VGRA_LOCK_EQUA_EB
      
  SUBROUTINE T2VGRA_LOCK_AXIS_EB
    
    INTEGER(ikind)::i_s,vOffsetA

    ! initialization
    LockAxi(1:NVMAX) = .FALSE.    

    IF(SolveField)THEN
       LockAxi( 1) = LockPoloidalMageticFieldOnAxis
       LockAxi( 2) = LockToroidalMageticFieldOnAxis
       LockAxi( 3) = LockRadialElectricFieldOnAxis
       LockAxi( 4) = LockPoloidalElectricFieldOnAxis
       LockAxi( 5) = LockToroidalElectricFieldOnAxis
    ENDIF

    IF(SolveElectron)THEN
       LockAxi( 6) = LockDensityOnAxis
       LockAxi( 7) = LockRaidalFluxOnAxis
       LockAxi( 8) = LockParallelFluxOnAxis
       LockAxi( 9) = LockToroidalFluxOnAxis
       LockAxi(10) = LockPoroidalFluxOnAxis 
       LockAxi(11) = LockPressureOnAxis 
       LockAxi(12) = LockRaidalHeatFluxOnAxis 
       LockAxi(13) = LockParallelHeatFluxOnAxis
       LockAxi(14) = LockToroidalHeatFluxOnAxis
       LockAxi(15) = LockPoroidalHeatFluxOnAxis 
    ENDIF

    IF(SolveIons)THEN
       DO i_s = 1,NSMAX-1
          vOffsetA = 10*i_s
          LockAxi( 6+vOffsetA) = LockDensityOnAxis
          LockAxi( 7+vOffsetA) = LockRaidalFluxOnAxis
          LockAxi( 8+vOffsetA) = LockParallelFluxOnAxis
          LockAxi( 9+vOffsetA) = LockToroidalFluxOnAxis
          LockAxi(10+vOffsetA) = LockPoroidalFluxOnAxis 
          LockAxi(11+vOffsetA) = LockPressureOnAxis 
          LockAxi(12+vOffsetA) = LockRaidalHeatFluxOnAxis 
          LockAxi(13+vOffsetA) = LockParallelHeatFluxOnAxis
          LockAxi(14+vOffsetA) = LockToroidalHeatFluxOnAxis
          LockAxi(15+vOffsetA) = LockPoroidalHeatFluxOnAxis 
       ENDDO
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2VGRA_LOCK_AXIS_EB
  
  SUBROUTINE  T2VGRA_LOCK_WALL_EB
    
    INTEGER(ikind)::i_s,vOffsetA

    ! initialization
    LockWal(1:NVMAX) = .FALSE.    

    IF(SolveField)THEN
       LockWal( 1) = LockPoloidalMageticFieldOnWall
       LockWal( 2) = LockToroidalMageticFieldOnWall
       LockWal( 3) = LockRadialElectricFieldOnWall
       LockWal( 4) = LockPoloidalElectricFieldOnWall
       LockWal( 5) = LockToroidalElectricFieldOnWall
    ENDIF

    IF(SolveElectron)THEN
       LockWal( 6) = LockDensityOnWall
       LockWal( 7) = LockRaidalFluxOnWall
       LockWal( 8) = LockParallelFluxOnWall
       LockWal( 9) = LockToroidalFluxOnWall
       LockWal(10) = LockPoroidalFluxOnWall 
       LockWal(11) = LockPressureOnWall
       LockWal(12) = LockRaidalHeatFluxOnWall
       LockWal(13) = LockParallelHeatFluxOnWall
       LockWal(14) = LockToroidalHeatFluxOnWall
       LockWal(15) = LockPoroidalHeatFluxOnWall 
    ENDIF

    IF(SolveIons)THEN
       DO i_s = 1,NSMAX-1
          vOffsetA = 10*i_s
          LockWal( 6+vOffsetA) = LockDensityOnWall
          LockWal( 7+vOffsetA) = LockRaidalFluxOnWall
          LockWal( 8+vOffsetA) = LockParallelFluxOnWall
          LockWal( 9+vOffsetA) = LockToroidalFluxOnWall
          LockWal(10+vOffsetA) = LockPoroidalFluxOnWall
          LockWal(11+vOffsetA) = LockPressureOnWall 
          LockWal(12+vOffsetA) = LockRaidalHeatFluxOnWall 
          LockWal(13+vOffsetA) = LockParallelHeatFluxOnWall
          LockWal(14+vOffsetA) = LockToroidalHeatFluxOnWall
          LockWal(15+vOffsetA) = LockPoroidalHeatFluxOnWall
       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE T2VGRA_LOCK_WALL_EB
  
  SUBROUTINE T2VGRA_OUTPUT

    USE T2COMM
    INTEGER(ikind):: i0sidi,i0sidj,i0widi,i0widj

    OPEN(10,FILE='TEST_VV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'VV=',i2vvvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_MS.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'MS=',i2msvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'AV=',i2avvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'AT=',i2atvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_DT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'DT=',i2dtvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'GV=',i2gvvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'GT=',i2gtvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ES.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'ES=',i2esvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_EV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'EV=',i2evvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ET.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'ET=',i2etvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_SS.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'SS=',i2ssvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I3ATWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,&
               'ATWT=',i3atwt(i0widi,i0sidi,i0sidj)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I3GTWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,&
               'GTWT=',i3gtwt(i0widi,i0sidi,i0sidj)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I3EVWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,&
               'EVWT=',i3evwt(i0widi,i0sidi,i0sidj)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I4ETWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widj = 1, i0wmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,'wj=',i0widj,&
               'ETWT=',i4etwt(i0widi,i0widj,i0sidi,i0sidj)
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_SS.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'SS=',i2ssvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    RETURN

  END SUBROUTINE T2VGRA_OUTPUT

  !
  ! The following section will be implemented 
  ! if EM-type equations will not work.
  ! phi-A type is suitable for integrated simulation with BOUT++ 
  !
  SUBROUTINE T2VGRA_VV_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_VV_COEF_PhiA
  
  SUBROUTINE T2VGRA_MS_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_MS_COEF_PhiA
  
  SUBROUTINE T2VGRA_AV_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_AV_COEF_PhiA
  
  SUBROUTINE T2VGRA_AT_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_AT_COEF_PhiA
  
  SUBROUTINE  T2VGRA_DT_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_DT_COEF_PhiA
  
  SUBROUTINE T2VGRA_GV_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_GV_COEF_PhiA
  
  SUBROUTINE T2VGRA_GT_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_GT_COEF_PhiA
  
  SUBROUTINE T2VGRA_ES_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_ES_COEF_PhiA
  
  SUBROUTINE T2VGRA_EV_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_EV_COEF_PhiA
  
  SUBROUTINE T2VGRA_ET_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_ET_COEF_PhiA

  SUBROUTINE T2VGRA_SS_COEF_PhiA
    RETURN
  END SUBROUTINE T2VGRA_SS_COEF_PhiA

  SUBROUTINE T2VGRA_AT_KVAL_PhiA
    RETURN
  END SUBROUTINE T2VGRA_AT_KVAL_PhiA

  SUBROUTINE T2VGRA_GT_KVAL_PhiA
    RETURN
  END SUBROUTINE T2VGRA_GT_KVAL_PhiA

  SUBROUTINE T2VGRA_EV_KVAL_PhiA
    RETURN
  END SUBROUTINE T2VGRA_EV_KVAL_PhiA

  SUBROUTINE T2VGRA_ET_KVAL_PhiA
    RETURN
  END SUBROUTINE T2VGRA_ET_KVAL_PhiA

END MODULE T2VGRA
