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
!        HaveVec
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
       HaveMat(:,:),HaveVec(:),&
       LockEqs(:),  LockAxi(:),  LockWal(:)
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
       CALL T2VGRA_VV_COEF_EB

       CALL T2VGRA_AT_KVAL_EB
       CALL T2VGRA_GT_KVAL_EB
       CALL T2VGRA_EV_KVAL_EB
       CALL T2VGRA_ET_KVAL_EB
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
       CALL T2VGRA_VV_COEF_PhiA

       CALL T2VGRA_AT_KVAL_PhiA
       CALL T2VGRA_GT_KVAL_PhiA
       CALL T2VGRA_EV_KVAL_PhiA
       CALL T2VGRA_ET_KVAL_PhiA
   
    ENDIF

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
          ALLOCATE(HaveVec(1:NVMAX),STAT=ierr);IF(ierr.NE.0)EXIT
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
    IF(ALLOCATED(HaveVec))          DEALLOCATE(HaveVec)
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
    
    !USE T2COMM,ONLY:NSMAX,NVMAX
    INTEGER(ikind)::i_s,i_v,j_v,vOffsetA
        
    !C initialization
    DO j_v = 1,NVMAX
    DO i_v = 1,NVMAX
       HaveMassScaCoef(1:NVMAX,1:NVMAX) = .FALSE.
    ENDDO
    ENDDO

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
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       j_v =  6 + vOffsetA;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetA;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetA;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetA;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
       
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v = 13 + vOffsetA;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v = 14 + vOffsetA;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
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
    
    !USE T2COMM,ONLY: NSMAX, NVMAX
    
    INTEGER(ikind):: i_s,i_v,j_v,vOffsetA
    
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
       
       vOffsetA  = 10*i_s
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       j_v =  6 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
 
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.

       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       !C
       !C EQUATION FOR Qb
       !C
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  8 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       !C
       !C EQUATION FOR Qt
       !C
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  9 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
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
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AT_COEF_EB
  
    
    INTEGER(ikind)::&
         i_s,i_v,j_v,vOffsetA
    
    ! initialization
    HaveAdveTenCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
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
       j_v =  9 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.


       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  9 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.       
       
       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  9 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  9 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 10 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 14 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       j_v = 15 + vOffsetA;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
       
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
         i_s,i_v,j_v,vOffsetA
    
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
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  6 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  6 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  6 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  6 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  6 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

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
         i_s,i_v,j_v,vOffsetA
    
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
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       j_v = 11 + vOffsetA;   HaveGradVecCoef(i_v,j_v) = .TRUE.
      
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
       j_v =  6;              HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11;              HaveGradVecCoef(i_v,j_v) = .TRUE.
       ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<

       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetA;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
       j_v =  6 + vOffsetA;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveGradVecCoef(i_v,j_v) = .TRUE.

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}       
       ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
       j_v =  6 + vOffsetA;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveGradVecCoef(i_v,j_v) = .TRUE.
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
  !
  !       where 
  !            vec{A}_{iv,jv} = grad g_{ik} dot ten{A}_{iv,jv,ik}
  !
  !                     LAST UPDATE 2014-05-23 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GT_COEF_EB
    
    INTEGER(ikind)::&
         i_s,i_v,j_v,vOffsetA
    
    ! initialization
    HaveGradTenCoef(1:NVMAX,1:NVMAX) = .FALSE.

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I
    i_v = 3                   ! Equation for E_{\zeta}
    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}

    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  6 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}

       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  6 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  6 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetA;   HaveGradTenCoef(i_v,j_v) = .TRUE.

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}

       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GT_COEF_EB
  ! kokomade
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION SCALAR ARRAY
  !C 
  !C                     2014-03-16
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2VGRA_ES_COEF_EB
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2esvt
    
    INTEGER(ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj,&
         i0vofi,i0vofj
    
    !C
    !C INITIALIZATION
    !C
    
    i2esvt(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C EQ_002
    !C
    
    i0vidi = 2
 
    !C Ep
    
    i0vidj = 4
    i2esvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax
       
       i0vofj = 10*i0sidj
       
       !C
       !C EQ_003
       !C
       
       i0vidi = 3
       
       !C Ft
       i0vidj = i0vofj - 1
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_004
       !C
       
       i0vidi = 4
       
       !C Fr
       i0vidj = i0vofj - 3
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Fp
       i0vidj = i0vofj
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_005
       !C
       
       i0vidi = 5
       
       !C N
       i0vidj = i0vofj - 4
       i2esvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_007
       !C
       
       i0vidi= i0vofi - 3
       
       !C Ep
       i0vidj = 4
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Er
       i0vidj = 5
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2esvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2esvt(i0vidi,i0vidj) = 1
       
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj          
          
          !C
          !C EQ_008
          !C
          
          i0vidi = i0vofi - 2
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2esvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = i0vofj + 3
          i2esvt(i0vidi,i0vidj) = 1
          
          !C
          !C EQ_009
          !C
          
          i0vidi = i0vofi - 1
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Fr
             i0vidj = i0vofj - 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C >>>> ANOMALOUS TRANSPORT >>>>
             
             !C Fbe
             i0vidj = 8
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Fte
             i0vidj = 9
             i2esvt(i0vidi,i0vidj) = 1
             
             !C <<<< ANOMALOUS TRANSPORT <<<<
             
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2esvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQ_010
       !C

       i0vidi = i0vofi

       !C Fb
       i0vidj = i0vofi - 2
       i2esvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2esvt(i0vidi,i0vidj) = 1

       !C Fp
       i0vidj = i0vofi
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi + 2
       
       !C Ep
       i0vidj = 4
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Er
       i0vidj = 5
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2esvt(i0vidi,i0vidj) = 1

       !C Qt
       i0vidj = i0vofi + 4
       i2esvt(i0vidi,i0vidj) = 1

       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj

          !C
          !C EQ_013
          !C
          
          i0vidi = i0vofi + 3
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2esvt(i0vidi,i0vidj) = 1             
             
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = i0vofj + 3
          i2esvt(i0vidi,i0vidj) = 1
          
          !C
          !C EQ_014
          !C
          
          i0vidi = i0vofi + 4
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Qr
             i0vidj = i0vofj + 2
             i2esvt(i0vidi,i0vidj) = 1
             
!             !C >>>> ANOMALOUS TRANSPORT >>>>!
!
!             !C Qbe
!             i0vidj = 13
!             i2esvt(i0vidi,i0vidj) = 1
!             
!             !C Qte
!             i0vidj = 14
!             i2esvt(i0vidi,i0vidj) = 1!
!
!             !C <<<< ANOMALOUS TRANSPORT <<<<

          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2esvt(i0vidi,i0vidj) = 1              
          
       ENDDO
       
       !C
       !C EQ_015
       !C
       
       i0vidi = i0vofi + 5
       
       !C Qb
       i0vidj = i0vofi + 3
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qp
       i0vidj = i0vofi + 5
       i2esvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_ES_COEF_EB

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION VECTOR ARRAY
  !C 
  !C                     2014-03-06
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EV_COEF_EB
    
    USE T2COMM,ONLY:&
         i0smax,i0wmax,i0vmax,i2evvt,i3evwt
    
    INTEGER(ikind)::&
         i0sidi,i0widi,i0vidi,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2evvt(         1:i0vmax,1:i0vmax) = 0
    i3evwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
    !C
    !C EQ_003
    !C
    
    i0vidi = 3
    
    !C PSI': (R)
    i0vidj = 1
    i2evvt(i0vidi,i0vidj) = 1
    
    i0widi = 2
    i3evwt(i0widi,i0vidi,i0vidj) = 1
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_008
       !C

       i0vidi = i0vofi - 2
       
       !C Fb: (B)
       i0vidj = i0vofi - 2
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C Et: (B) (Ub) (Qb/P)
       i0vidj = 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Ep: (B) (Ub) (Qb/P)
       i0vidj = 4
       i2evvt(i0vidi,i0vidj) = 1             
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb: (B)
       i0vidj = i0vofi - 2
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb: (B)
       i0vidj = i0vofi + 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = i0vofi + 4
       
       !C Et: B, Ub, Wb
       i0vidj = 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Ep: B, Ub, Wb
       i0vidj = 4
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_EV_COEF_EB

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION TENSOR ARRAY
  !C 
  !C                     2014-02-20 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ET_COEF_EB

    USE T2COMM,ONLY:&
         i0smax,i0wmax,i0vmax,i2etvt,i4etwt
    
    INTEGER(ikind)::&
         i0sidi,i0widi,i0vidi,&
                i0widj,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2etvt(                  1:i0vmax,1:i0vmax) = 0
    i4etwt(1:i0wmax,1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C Ft: (B,B)
       i0vidj = i0vofi - 1
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Fp: (B,B)
       i0vidj = i0vofi
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qt (B,B)
       i0vidj = i0vofi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qp (B,B)
       i0vidj = i0vofi + 5
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C Ft: (B,B), (Ub,B)
       i0vidj = i0vofi - 1
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Fp: (B,B), (Ub,B)
       i0vidj = i0vofi
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qt: (B,B), (Ub,B)
       i0vidj = i0vofi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Qp: (B,B), (Ub,B)
       i0vidj = i0vofi + 5
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = i0vofi + 3
       
       !C Ft: (B,B)
       i0vidj = i0vofi - 1
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Fp: (B,B)
       i0vidj = i0vofi
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Qt: (B,B)
       i0vidj = i0vofi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qp: (B,B)
       i0vidj = i0vofi + 5
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
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
  SUBROUTINE T2VGRA_VV_COEF_EB
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2vvvt
    
    INTEGER(ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj,&
         i0vofi,i0vofj
    
    !C
    !C INITIALIZATION
    !C
    
    i2vvvt(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_002
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2vvvt(i0vidi,i0vidj) = 1

    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1

    !C Er 
    i0vidj = 5
    i2vvvt(i0vidi,i0vidj) = 1

    !C
    !C EQ_003
    !C
    
    i0vidi = 3

    !C PSI'
    i0vidj = 1
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2vvvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax

       i0vofj = 10*i0sidj

       !C Ft
       i0vidj = i0vofj - 1
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQ_004
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1


    DO i0sidj = 1, i0smax
    
       i0vofj = 10*i0sidj

       !C Fr
       i0vidj = i0vofj - 3
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Fp
       i0vidj = i0vofj 
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQ_005
    !C

    i0vidi = 5
    
    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1

    !C Er
    i0vidj = 5
    i2vvvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax
       
       i0vofj = 10*i0sidj
       
       !C N
       i0vidj = i0vofj - 4
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi

       !C
       !C EQ_006
       !C
       
       i0vidi = i0vofi - 4
       
       !C N
       i0vidj = i0vofi - 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_007
       !C
       
       i0vidi= i0vofi - 3
       
       !C Ep
       i0vidj = 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fr
       i0vidj = i0vofi - 3
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = i0vofi - 2 
       i2vvvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C P 
       i0vidj = i0vofi + 1
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
      
             !C Ft
             i0vidj = i0vofj - 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj 
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qt
             i0vidj = i0vofj + 4
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1
          
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2vvvt(i0vidi,i0vidj) = 1

          !C Qb
          i0vidj = i0vofj + 3
          i2vvvt(i0vidi,i0vidj) = 1

       ENDDO
       
       !C
       !C EQ_009
       !C

       i0vidi = i0vofi - 1
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fr
             i0vidj = i0vofj - 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fb
             i0vidj = i0vofj - 2
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fp
             i0vidj = i0vofj
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Qb
             i0vidj = i0vofj + 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C >>>> ANOMALOUS TRANSPORT >>>>
             
             !C Ne
             i0vidj = 6
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fbe
             i0vidj = 8
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fte 
             i0vidj = 9
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Pe
             i0vidj = 11
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C <<<< ANOMALOUS TRANSPORT <<<<
          
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2vvvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQ_010
       !C

       i0vidi = i0vofi

       !C Fb
       i0vidj = i0vofi - 2 
       i2vvvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fp
       i0vidj = i0vofi
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C N
       i0vidj = i0vofi - 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Ft
       i0vidj = i0vofi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fp
       i0vidj = i0vofi
       i2vvvt(i0vidi,i0vidj) = 1

       !C P
       i0vidj = i0vofi + 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = i0vofi + 3
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qp
       i0vidj = i0vofi + 5
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi + 2
       
       !C Ep
       i0vidj = 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C N
       i0vidj = i0vofi - 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C P
       i0vidj = i0vofi + 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qr
       i0vidj = i0vofi + 2
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = i0vofi + 3
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qt
       i0vidj = i0vofi + 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1             

             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Ft
             i0vidj = i0vofj - 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qt
             i0vidj = i0vofj + 4
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1

          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = i0vofj + 3
          i2vvvt(i0vidi,i0vidj) = 1

       ENDDO

       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4

       DO i0sidj = 1, i0smax
       
          i0vofj = 10*i0sidj

          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fb
             i0vidj = i0vofj - 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qr
             i0vidj = i0vofj + 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qb
             i0vidj = i0vofj + 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1

!             !C >>>>  ANOMALOUS TRANSPORT >>>>
!
!             !C Ne
!             i0vidj = 6
!             i2vvvt(i0vidi,i0vidj) = 1
!             
!             !C Pe
!             i0vidj = 11
!             i2vvvt(i0vidi,i0vidj) = 1
!             
!             !C Qbe
!             i0vidj = 13
!             i2vvvt(i0vidi,i0vidj) = 1
!             
!             !C Qte 
!             i0vidj = 14
!             i2vvvt(i0vidi,i0vidj) = 1
!             
!             !C <<<<  ANOMALOUS TRANSPORT <<<<

          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2vvvt(i0vidi,i0vidj) = 1              
                    
       ENDDO
       
       !C
       !C EQ_015
       !C

       i0vidi = i0vofi + 5
       
       !C Qb
       i0vidj = i0vofi + 3 
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qp
       i0vidj = i0vofi + 5
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_VV_COEF_EB

  SUBROUTINE T2VGRA_AT_KVAL_EB
    RETURN
  END SUBROUTINE T2VGRA_AT_KVAL_EB

  SUBROUTINE T2VGRA_GT_KVAL_EB
    RETURN
  END SUBROUTINE T2VGRA_GT_KVAL_EB

  SUBROUTINE T2VGRA_EV_KVAL_EB
    RETURN
  END SUBROUTINE T2VGRA_EV_KVAL_EB
    
  SUBROUTINE T2VGRA_ET_KVAL_EB
    RETURN
  END SUBROUTINE T2VGRA_ET_KVAL_EB
  
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
