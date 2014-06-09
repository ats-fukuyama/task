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
!        NSMAX
!        NVMAX
!        NKMAX
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
!   T2INTG sets up following variables:
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
!   through subroutine: T2VGRA_EXCUTE
!
! -------------------------------------------------------------------
MODULE T2VGRA
  
  USE T2CNST,ONLY: ikind,rkind
  
  IMPLICIT NONE
  
  PRIVATE 

  PUBLIC T2VGRA_EXECUTE
  
CONTAINS 
  
  !-------------------------------------------------------------------
  !
  !       T2VGRA_EXCUTE (PUBLIC)
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EXECUTE
        
    USE T2COMM,ONLY: UsePotentialDescription
    
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
    
    !IF(idfile.ge.5) &
    CALL T2VGRA_OUTPUT
    
    RETURN
    
  END SUBROUTINE T2VGRA_EXECUTE
  
      
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
  !                     
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_MS_COEF_EB
    
    USE T2COMM,ONLY:NSMAX,NVMAX,HaveMassScaCoef
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
       vOffsetB = vOffsetA

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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AV_COEF_EB
    
    USE T2COMM,ONLY: NSMAX,NVMAX,HaveAdveVecCoef
    
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

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  8 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 11 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.

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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AT_COEF_EB
    
    USE T2COMM,ONLY: NSMAX, NVMAX, NKMAX,&
         &           HaveAdveTenCoef, HaveAdveTenKval
    
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
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{rho}

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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_DT_COEF_EB
    
    USE T2COMM,ONLY: NSMAX, NVMAX, HaveDiffTenCoef
        
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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GV_COEF_EB

    USE T2COMM, ONLY: NSMAX,NVMAX,HaveGradVecCoef

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
       j_v =  6;              HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v = 11;              HaveGradVecCoef(i_v,j_v) = .TRUE.
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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GT_COEF_EB
    
    USE T2COMM, ONLY: NSMAX, NVMAX, NKMAX,&
         &            HaveGradTenCoef,HaveGradTenKval

    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,kOffsetX,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveGradTenCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveGradTenKval(1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.

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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ES_COEF_EB

    USE T2COMM, ONLY: NSMAX, NVMAX,&
         &            HaveExciScaCoef
    
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
  !                     LAST UPDATE     2014-05-23 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EV_COEF_EB

    USE T2COMM, ONLY: NSMAX, NVMAX, NKMAX,&
         &            HaveExciVecCoef,HaveExciVecKval

    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,kOffsetX,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveExciVecCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveExciVecKVal(1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I

    i_v = 3                   ! Equation for E_{\zeta}
    j_v = 1;                  HaveExciVecCoef(i_v,j_v) = .TRUE.
    i_k = 2;                  HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
   
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
       j_v =  8 + vOffsetB;   HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  3;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       j_v =  4;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       j_v =  8 + vOffsetB;   HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       j_v = 13 + vOffsetB;   HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  3;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       j_v =  4;              HaveExciVecCoef(i_v,j_v) = .TRUE.
       i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       
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
  !                     LAST UPDATE     2014-06-08 H.Seto
  !                     OPERATION CHECK 2014-06-08 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ET_COEF_EB

    USE T2COMM, ONLY: NSMAX, NVMAX, NKMAX,&
         &            HaveExciTenCoef,HaveExciTenKval
    
    INTEGER(ikind)::&
         & i_s,i_v,i_k,kOffsetX,vOffsetA,&
         &     j_v,j_k,kOffsetY,vOffsetB

    ! initialization
    HaveExciTenCoef(                1:NVMAX,1:NVMAX) = .FALSE.
    HaveExciTenKval(1:NKMAX,1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
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
       kOffsetY = kOffsetX
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}

       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  9 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k = 1
       j_k = 1;               HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 10 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k = 1
       j_k = 1;               HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 14 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k = 1
       j_k = 1;               HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 15 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k = 1
       j_k = 1;               HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  9 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k = 1
       j_k = 1;               HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       i_k = 3 + kOffsetX
       j_k = 1;               HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

       j_v = 10 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       i_k =  3 + kOffsetX
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 14 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       i_k =  3 + kOffsetX
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 15 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       i_k =  3 + kOffsetX
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  9 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 10 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 14 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
       j_v = 15 + vOffsetA;   HaveExciTenCoef(i_v,j_v) = .TRUE.
       i_k =  1
       j_k =  1;              HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
    
    ENDDO
    
    RETURN 

  END SUBROUTINE T2VGRA_ET_COEF_EB
  
  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF SOURCE SCALAR TERMS 
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveExciTenCoef(iv,jv) = .TRUE.", 
  !
  !       the term "{S}_{iv,jv}"
  !
  !                           is exist in equation system.
  !
  !                     LAST UPDATE     2014-05-26 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_SS_COEF_EB
    
    USE T2COMM,ONLY: NSMAX,NVMAX,&
         &           HaveSourScaCoef
    
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    HaveSourScaCoef(1:NVMAX,1:NVMAX) = .FALSE.

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    j_v = 1;                  HaveSourScaCoef(i_v,j_v) = .TRUE.

    i_v = 2                   ! Equation for I
    j_v = 2;                  HaveSourScaCoef(i_v,j_v) = .TRUE.

    i_v = 3                   ! Equation for E_{\zeta}
    j_v = 3;                  HaveSourScaCoef(i_v,j_v) = .TRUE.

    i_v = 4                   ! Equation for E_{\chi}
    j_v = 4;                  HaveSourScaCoef(i_v,j_v) = .TRUE.
    
    i_v = 5                   ! Equation for E_{\rho}    
    j_v = 5;                  HaveSourScaCoef(i_v,j_v) = .TRUE.

    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !       
    DO i_s = 0, NSMAX-1
       
       vOffsetA = 10*i_s
       vOffsetB = vOffsetA
 
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       j_v =  6 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.

       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       j_v =  7 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.

       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       j_v =  9 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       j_v = 10 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v = 11 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
       j_v = 12 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v = 13 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v = 14 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.
       
       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
       j_v = 15 + vOffsetB;   HaveSourScaCoef(i_v,j_v) = .TRUE.

    ENDDO
    
    RETURN

  END SUBROUTINE T2VGRA_SS_COEF_EB

    
  !-------------------------------------------------------------------
  !
  ! CREATE VARIABLE GRAPH OF STIFFNESS MATRIX
  ! 
  !                    E-B VERSION
  !
  !       IF "HaveMat(iv,jv) = .TRUE.", 
  !
  !       a component "mat{A}_{*,*,iv,jv}"
  !
  !                           is exist in stiffness matrix .
  !
  !                     LAST UPDATE     2014-05-26 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------  
  SUBROUTINE T2VGRA_VAR_MAT
    
    USE T2COMM,ONLY:NVMAX,HaveMat,&
         &          HaveMassScaCoef,HaveAdveVecCoef,HaveAdveTenCoef,&
         &          HaveDiffTenCoef,HaveGradVecCoef,HaveGradTenCoef,&
         &          HaveExciScaCoef,HaveExciVecCoef,HaveExciTenCoef,&
         &          HaveSourScaCoef
    
    INTEGER(ikind)::i_v
    
    HaveMat =  HaveMassScaCoef .OR. HaveAdveVecCoef .OR. &
         &     HaveAdveTenCoef .OR. HaveDiffTenCoef .OR. &
         &     HaveGradVecCoef .OR. HaveGradTenCoef .OR. &
         &     HaveExciScaCoef .OR. HaveExciVecCoef .OR. &
         &     HaveExciTenCoef
    
    ! eliminate singularity in variable matrix 
    DO i_v = 1, NVMAX
       HaveMat(i_v,i_v) = .TRUE.
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_VAR_MAT

  !-------------------------------------------------------------------
  ! 
  ! CREATE TABLE ABOUT WHICH EQUATIONS ARE TO BE LOCKED
  !
  !
  !                     LAST UPDATE     2014-05-24 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_LOCK_EQUA_EB
    
    USE T2COMM, ONLY: NSMAX,NVMAX, LockEqs,&
         &            SolveField,    SolveElectron, SolveIons,&
         &            SolveDensity,  SolveFlux,&
         &            SolvePressure, SolveHeatFlux

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
      
  !-------------------------------------------------------------------
  ! 
  ! CREATE TABLE ABOUT WHICH VARIABLES ARE TO BE LOCKED ON AXIS
  !
  !
  !                     LAST UPDATE     2014-05-24 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_LOCK_AXIS_EB
    
    USE T2COMM, ONLY: NSMAX,NVMAX,LockAxi,LockEqs,&
         !
         &            SolveField,&
         &            SolveElectron,&
         &            SolveIons,&
         !
         &            LockPoloidalMageticFieldOnAxis,&
         &            LockToroidalMageticFieldOnAxis,&
         &            LockRadialElectricFieldOnAxis,&
         &            LockPoloidalElectricFieldOnAxis,&
         &            LockToroidalElectricFieldOnAxis,&
         &            LockDensityOnAxis,&
         &            LockRadialFluxOnAxis,&
         &            LockParallelFluxOnAxis,&
         &            LockToroidalFluxOnAxis,&
         &            LockPoroidalFluxOnAxis,&
         &            LockPressureOnAxis,&
         &            LockRadialHeatFluxOnAxis,&
         &            LockParallelHeatFluxOnAxis,&
         &            LockToroidalHeatFluxOnAxis,&
         &            LockPoroidalHeatFluxOnAxis
    
    INTEGER(ikind)::i_s,vOffsetA

    ! initialization
    LockAxi(1:NVMAX) = .FALSE.    
    
    IF(SolveField)THEN
       IF(.NOT.LockEqs(1))&
            LockAxi(   1) = LockPoloidalMageticFieldOnAxis
       IF(.NOT.LockEqs(2))&
            LockAxi(   2) = LockToroidalMageticFieldOnAxis
       IF(.NOT.LockEqs(3))&
            LockAxi(   3) = LockRadialElectricFieldOnAxis
       IF(.NOT.LockEqs(4))&
            LockAxi(   4) = LockPoloidalElectricFieldOnAxis
       IF(.NOT.LockEqs(5))&
            LockAxi(   5) = LockToroidalElectricFieldOnAxis
    ENDIF

    IF(SolveElectron)THEN
       IF(.NOT.LockEqs( 6))&
            LockAxi(    6) = LockDensityOnAxis
       IF(.NOT.LockEqs( 7))&
            LockAxi(    7) = LockRadialFluxOnAxis
       IF(.NOT.LockEqs( 8))&
            LockAxi(    8) = LockParallelFluxOnAxis
       IF(.NOT.LockEqs( 9))&
            LockAxi(    9) = LockToroidalFluxOnAxis
       IF(.NOT.LockEqs(10))&
            LockAxi(   10) = LockPoroidalFluxOnAxis 
       IF(.NOT.LockEqs(11))&
            LockAxi(   11) = LockPressureOnAxis 
       IF(.NOT.LockEqs(12))&
            LockAxi(   12) = LockRadialHeatFluxOnAxis 
       IF(.NOT.LockEqs(13))&
            LockAxi(   13) = LockParallelHeatFluxOnAxis
       IF(.NOT.LockEqs(14))&
            LockAxi(   14) = LockToroidalHeatFluxOnAxis
       IF(.NOT.LockEqs(15))&
            LockAxi(   15) = LockPoroidalHeatFluxOnAxis 
    ENDIF

    IF(SolveIons)THEN
       DO i_s = 1,NSMAX-1
          vOffsetA = 10*i_s
          IF(.NOT.LockEqs( 6+vOffsetA))&
               &  LockAxi( 6+vOffsetA) = LockDensityOnAxis
          IF(.NOT.LockEqs( 7+vOffsetA))&
               &  LockAxi( 7+vOffsetA) = LockRadialFluxOnAxis
          IF(.NOT.LockEqs( 8+vOffsetA))&
               &  LockAxi( 8+vOffsetA) = LockParallelFluxOnAxis
          IF(.NOT.LockEqs( 9+vOffsetA))&
               &  LockAxi( 9+vOffsetA) = LockToroidalFluxOnAxis
          IF(.NOT.LockEqs(10+vOffsetA))&
               &  LockAxi(10+vOffsetA) = LockPoroidalFluxOnAxis 
          IF(.NOT.LockEqs(11+vOffsetA))&
               &  LockAxi(11+vOffsetA) = LockPressureOnAxis 
          IF(.NOT.LockEqs(12+vOffsetA))&
               &  LockAxi(12+vOffsetA) = LockRadialHeatFluxOnAxis 
          IF(.NOT.LockEqs(13+vOffsetA))&
               &  LockAxi(13+vOffsetA) = LockParallelHeatFluxOnAxis
          IF(.NOT.LockEqs(14+vOffsetA))&
               &  LockAxi(14+vOffsetA) = LockToroidalHeatFluxOnAxis
          IF(.NOT.LockEqs(15+vOffsetA))&
               &  LockAxi(15+vOffsetA) = LockPoroidalHeatFluxOnAxis 
       ENDDO
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2VGRA_LOCK_AXIS_EB
  
  !-------------------------------------------------------------------
  ! 
  ! CREATE TABLE ABOUT WHICH VARIABLES ARE TO BE LOCKED ON WALL
  !
  !
  !                     LAST UPDATE     2014-05-24 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE  T2VGRA_LOCK_WALL_EB
    
    USE T2COMM,ONLY: NSMAX,NVMAX,LockWal,LockEqs,&
         !
         &           SolveField,&
         &           SolveElectron,&
         &           SolveIons,&
         !
         &           LockPoloidalMageticFieldOnWall,&
         &           LockToroidalMageticFieldOnWall,&
         &           LockRadialElectricFieldOnWall,&
         &           LockPoloidalElectricFieldOnWall,&
         &           LockToroidalElectricFieldOnWall,&
         &           LockDensityOnWall,&
         &           LockRadialFluxOnWall,&
         &           LockParallelFluxOnWall,&
         &           LockToroidalFluxOnWall,&
         &           LockPoroidalFluxOnWall,&
         &           LockPressureOnWall,&
         &           LockRadialHeatFluxOnWall,&
         &           LockParallelHeatFluxOnWall,&
         &           LockToroidalHeatFluxOnWall,&
         &           LockPoroidalHeatFluxOnWall

    INTEGER(ikind)::i_s,vOffsetA

    ! initialization
    LockWal(1:NVMAX) = .FALSE.    

    IF(SolveField)THEN
       IF(.NOT.LockEqs( 1))&
            &  LockWal( 1) = LockPoloidalMageticFieldOnWall
       IF(.NOT.LockEqs( 2))&
            &  LockWal( 2) = LockToroidalMageticFieldOnWall
       IF(.NOT.LockEqs( 3))&
            &  LockWal( 3) = LockRadialElectricFieldOnWall
       IF(.NOT.LockEqs( 4))&
            &  LockWal( 4) = LockPoloidalElectricFieldOnWall
       IF(.NOT.LockEqs( 5))&
            &  LockWal( 5) = LockToroidalElectricFieldOnWall
    ENDIF

    IF(SolveElectron)THEN
       IF(.NOT.LockEqs( 6))&
            &  LockWal( 6) = LockDensityOnWall
       IF(.NOT.LockEqs( 7))&
            &  LockWal( 7) = LockRadialFluxOnWall
       IF(.NOT.LockEqs( 8))&
            &  LockWal( 8) = LockParallelFluxOnWall
       IF(.NOT.LockEqs( 9))&
            &  LockWal( 9) = LockToroidalFluxOnWall
       IF(.NOT.LockEqs(10))&
            &  LockWal(10) = LockPoroidalFluxOnWall 
       IF(.NOT.LockEqs(11))&
            &  LockWal(11) = LockPressureOnWall
       IF(.NOT.LockEqs(12))&
            &  LockWal(12) = LockRadialHeatFluxOnWall
       IF(.NOT.LockEqs(13))&
            &  LockWal(13) = LockParallelHeatFluxOnWall
       IF(.NOT.LockEqs(14))&
            &  LockWal(14) = LockToroidalHeatFluxOnWall
       IF(.NOT.LockEqs(15))&
            &  LockWal(15) = LockPoroidalHeatFluxOnWall 
    ENDIF

    IF(SolveIons)THEN
       DO i_s = 1,NSMAX-1
          vOffsetA = 10*i_s
          IF(.NOT.LockEqs( 6+vOffsetA))&
            &     LockWal( 6+vOffsetA) = LockDensityOnWall
          IF(.NOT.LockEqs( 7+vOffsetA))&
            &     LockWal( 7+vOffsetA) = LockRadialFluxOnWall
          IF(.NOT.LockEqs( 8+vOffsetA))&
            &     LockWal( 8+vOffsetA) = LockParallelFluxOnWall
          IF(.NOT.LockEqs( 9+vOffsetA))&
            &     LockWal( 9+vOffsetA) = LockToroidalFluxOnWall
          IF(.NOT.LockEqs(10+vOffsetA))&
            &     LockWal(10+vOffsetA) = LockPoroidalFluxOnWall
          IF(.NOT.LockEqs(11+vOffsetA))&
            &     LockWal(11+vOffsetA) = LockPressureOnWall 
          IF(.NOT.LockEqs(12+vOffsetA))&
            &     LockWal(12+vOffsetA) = LockRadialHeatFluxOnWall 
          IF(.NOT.LockEqs(13+vOffsetA))&
            &     LockWal(13+vOffsetA) = LockParallelHeatFluxOnWall
          IF(.NOT.LockEqs(14+vOffsetA))&
            &     LockWal(14+vOffsetA) = LockToroidalHeatFluxOnWall
          IF(.NOT.LockEqs(15+vOffsetA))&
            &     LockWal(15+vOffsetA) = LockPoroidalHeatFluxOnWall
       ENDDO
    ENDIF

    RETURN
  END SUBROUTINE T2VGRA_LOCK_WALL_EB
  
  SUBROUTINE T2VGRA_OUTPUT

    USE T2COMM,ONLY: NVMAX,NKMAX,&
         &           HaveMassScaCoef,HaveAdveVecCoef,HaveAdveTenCoef,&
         &           HaveDiffTenCoef,HaveGradVecCoef,HaveGradTenCoef,&
         &           HaveExciScaCoef,HaveExciVecCoef,HaveExciTenCoef,&
         &           HaveSourScaCoef,&
         !
         &           HaveAdveTenKval,HaveGradTenKval,HaveExciVecKval,&
         &           HaveExciTenKval,&
         !
         &           HaveMat, LockEqs, LockAxi,  LockWal

    INTEGER(ikind):: i_v,j_v,i_k,j_k

    OPEN(10,FILE='TEST_VMC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'HaveMat=',HaveMat(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_MSC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'MS=',HaveMassScaCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AVC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'AV=',HaveAdveVecCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ATC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'AT=',HaveAdveTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_DTC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'DT=',HaveDiffTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GVC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'GV=',HaveGradVecCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GTC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'GT=',HaveGradTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ESC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'ES=',HaveExciScaCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_EVC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'EV=',HaveExciVecCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ETC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'ET=',HaveExciTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_SSC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'SS=',HaveSourScaCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ATK.dat')
    DO i_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'iv=',i_v,'jv=',j_v,&
               'AT=',HaveAdveTenKval(i_k,i_v,j_v)
       ENDDO
       ENDDO
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_GTK.dat')
    DO i_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'iv=',i_v,'jv=',j_v,&
               'GT=',HaveGradTenKval(i_k,i_v,j_v)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_EVK.dat')
    DO i_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'iv=',i_v,'jv=',j_v,&
               'EV=',HaveExciVecKval(i_k,i_v,j_v)
       ENDDO
       ENDDO
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_ETK.dat')
    DO i_k = 1, NKMAX
    DO j_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'jk=',j_k,'iv=',i_v,'jv=',j_v,&
               'ET=',HaveExciTenKval(i_k,j_k,i_v,j_v)
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_LEQ.dat')
    DO i_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'LEQ=',LockEqs(i_v)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_LAX.dat')
    DO i_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'LAX=',LockAxi(i_v)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_LWL.dat')
    DO i_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'LEQ=',LockWal(i_v)
    ENDDO
    CLOSE(10)
    RETURN

  END SUBROUTINE T2VGRA_OUTPUT

  !
  ! The following section will be implemented 
  ! if EM-type equations will not work.
  ! phi-A type is also suitable for integrated simulation with BOUT++
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

  SUBROUTINE  T2VGRA_LOCK_EQUA_PhiA
    RETURN
  END SUBROUTINE T2VGRA_LOCK_EQUA_PhiA

  SUBROUTINE T2VGRA_LOCK_AXIS_PhiA
    RETURN
  END SUBROUTINE T2VGRA_LOCK_AXIS_PhiA
  
  SUBROUTINE T2VGRA_LOCK_WALL_PhiA
    RETURN
  END SUBROUTINE T2VGRA_LOCK_WALL_PhiA

END MODULE T2VGRA
