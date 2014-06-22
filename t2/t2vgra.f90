!------------------------------------------------------------------
!
!   MODULE FOR GENERATING VARIABLE-VARIABLE GRAPH 
!       
!                      
!                       LAST UPDATE 2014-06-22 H.Seto
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
  !       T2VGRA_EXCUTE 
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EXECUTE
    
    USE T2COMM,ONLY:&
         UsePotentialDescription,CoordinateSwitch,&
         HaveMassScaCoef,HaveAdveVecCoef,HaveAdveTenCoef,&
         HaveDiffTenCoef,HaveGradVecCoef,HaveGradTenCoef,&
         HaveExciScaCoef,HaveExciVecCoef,HaveExciTenCoef,&
         HaveSourScaCoef,&
         !
         HaveAdveTenKval,HaveGradTenKval,HaveExciVecKval,&
         HaveExciTenKval,&
         !
         HaveMat,LockEqs,LockAxi,LockWal,TestLEQ,TestLAX,TestLWL
    
    SELECT CASE (CoordinateSwitch)
    CASE(1)
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
       
    CASE (2)

       HaveMassScaCoef(1,1) =.TRUE.
       HaveAdveVecCoef(1,1) =.TRUE.
       HaveAdveTenCoef(1,1) =.TRUE.
       HaveDiffTenCoef(1,1) =.TRUE.
       HaveGradVecCoef(1,1) =.TRUE.
       HaveGradTenCoef(1,1) =.TRUE.
       HaveExciScaCoef(1,1) =.TRUE.
       HaveExciVecCoef(1,1) =.TRUE.
       HaveExciTenCoef(1,1) =.TRUE.
       HaveSourScaCoef(1,1) =.TRUE.
       !
       HaveAdveTenKval(1,1,1) =.TRUE.
       HaveGradTenKval(1,1,1) =.TRUE.
       HaveExciVecKval(1,1,1) =.TRUE.
       HaveExciTenKval(1,1,1,1) =.TRUE.
       !
       HaveMat(1,1) =.TRUE.
       LockEqs(1)   = TestLEQ
       LockAxi(1)   = TestLAX
       LockWal(1)   = TestLWL

    END SELECT
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_MS_COEF_EB
    
    USE T2COMM,ONLY:NSMAX,NVMAX,NFMAX,HaveMassScaCoef,EqSet
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    !C initialization
    HaveMassScaCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    SELECT CASE (EqSet)
    CASE (1)
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
          
          vOffsetA = 10*i_s   + NFMAX
          vOffsetB = vOffsetA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  8 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  9 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    CASE (2)
       i_v = 1                   ! Equation for psi'
       j_v = 3;                  HaveMassScaCoef(i_v,j_v) = .TRUE.
    
       i_v = 2                   ! Equation for I
       j_v = 4;                  HaveMassScaCoef(i_v,j_v) = .TRUE.

       i_v = 3                   ! Equation for E_{\tor}
       j_v = 1;                  HaveMassScaCoef(i_v,j_v) = .TRUE.

       i_v = 4                   ! Equation for {E}_{\pol}
       j_v = 2;                  HaveMassScaCoef(i_v,j_v) = .TRUE.

       i_v = 5                   ! Equation for E_{\rad}

       i_v = 6                   ! Equation for E^{\pol}

       !
       ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
       !
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s   + NFMAX
          vOffsetB = vOffsetA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
       
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
      
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  9 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  8 + vOffsetB;   HaveMassScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
       
    END SELECT
    

    
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AV_COEF_EB
    
    USE T2COMM,ONLY: NSMAX,NVMAX,NFMAX,HaveAdveVecCoef,EqSet
    
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    HaveAdveVecCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    SELECT CASE (EqSet)

    CASE (1)
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
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  4 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
       ENDDO
    CASE (2)
       i_v = 1                   ! Equation for psi'
       j_v = 1;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 3;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    
       i_v = 2                   ! Equation for I
       j_v = 4;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
    
       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 1;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.   
   
       i_v = 4                   ! Equation for E_{\chi}
       j_v = 2;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
       
       i_v = 5                   ! Equation for E_{\rho}
       j_v = 4;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 5;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.
       j_v = 6;                  HaveAdveVecCoef(i_v,j_v) = .TRUE.

       i_v = 6                   ! Equation for E^{\chi}
       
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
                 
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  4 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveAdveVecCoef(i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
          
       ENDDO

    END SELECT
    
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AT_COEF_EB
    
    USE T2COMM,ONLY: NSMAX,NVMAX,NFMAX,NKMAX,EqSet,&
         &           HaveAdveTenCoef, HaveAdveTenKval
    
    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveAdveTenCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveAdveTenKval(1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
    SELECT CASE (EqSet)

    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.    
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.   
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.

          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.    
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.   
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{rho}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  4 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetB;   HaveAdveTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveAdveTenKval(i_k,i_v,j_v) = .TRUE.
       
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
          
       ENDDO
    END SELECT
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_DT_COEF_EB
    
    USE T2COMM,ONLY: NSMAX, NVMAX, NFMAX, EqSet, HaveDiffTenCoef
    
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    
    HaveDiffTenCoef(1:NVMAX,1:NVMAX) = .FALSE.

    SELECT CASE (EqSet)

    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA

          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
       
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
    
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA

          i_v =  1 + vOffsetA    ! Equation for n_{a}

          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.

          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveDiffTenCoef(i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    END SELECT

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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GV_COEF_EB

    USE T2COMM, ONLY: NSMAX,NVMAX,NFMAX,HaveGradVecCoef,EqSet

    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         &     j_v,vOffsetB
    
    ! initialization
    
    HaveGradVecCoef(1:NVMAX,1:NVMAX) = .FALSE.
    
    !
    ! variables as field (from i_v= 1 to i_v = NFMAX)
    !
    SELECT CASE (EqSet)
    CASE (1)
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
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
       
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
       
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  1 + NFMAX;      HaveGradVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + NFMAX;      HaveGradVecCoef(i_v,j_v) = .TRUE.
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
       
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  1 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}

          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}       
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  1 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a\zeta}
          
       ENDDO
  
    CASE (2)
       i_v =  1                  ! Equation for psi'
       
       i_v =  2                  ! Equation for I
       j_v =  2;                 HaveGradVecCoef(i_v,j_v) = .TRUE.

       i_v =  3                  ! Equation for E_{\zeta}
       j_v =  3;                 HaveGradVecCoef(i_v,j_v) = .TRUE.
    
       i_v =  4                  ! Equation for I
       j_v =  4;                 HaveGradVecCoef(i_v,j_v) = .TRUE.
       j_v =  5;                 HaveGradVecCoef(i_v,j_v) = .TRUE.

       i_v =  5                  ! Equation for E_{\rho}
       i_v =  6                  ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
       
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  1 + NFMAX;      HaveGradVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + NFMAX;      HaveGradVecCoef(i_v,j_v) = .TRUE.
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
       
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
       
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  1 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}

          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}       
          j_v =  1 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradVecCoef(i_v,j_v) = .TRUE.
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a\zeta}
          
       ENDDO
    END SELECT
        
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GT_COEF_EB
    
    USE T2COMM, ONLY: NSMAX, NVMAX, NFMAX, NKMAX, EqSet,&
         &            HaveGradTenCoef,HaveGradTenKval

    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,kOffsetX,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveGradTenCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveGradTenKval(1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
    SELECT CASE (EqSet)

    CASE (1)
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
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA 
          kOffsetX =  2*i_s
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  1 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}      
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
       
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = 5)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
       !    
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA 
          kOffsetX =  2*i_s
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  1 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  6 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveGradTenCoef(i_v,j_v) = .TRUE.
          i_k =  1           ;   HaveGradTenKval(i_k,i_v,j_v) = .TRUE.
       
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}      
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
          
       ENDDO
    END SELECT

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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ES_COEF_EB
    
    USE T2COMM, ONLY: NSMAX, NVMAX,NFMAX,EqSet,&
         &            HaveExciScaCoef,EqSet
    
    INTEGER(ikind)::&
         i_s,i_v,vOffsetA,&
         j_s,j_v,vOffsetB
    
    !initialization
    
    HaveExciScaCoef(1:NVMAX,1:NVMAX) = .FALSE.

    SELECT CASE (EqSet)

    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !    
       i_v = 1                   ! Equation for psi'
       
       i_v = 2                   ! Equation for I
       j_v = 4;                  HaveExciScaCoef(i_v,j_v) = .TRUE.       
       
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}           
       i_v = 5                   ! Equation for E_{\rho}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  5;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  3 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  8 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO

          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  2 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  3 + NFMAX;      HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4 + NFMAX;      HaveExciScaCoef(i_v,j_v) = .TRUE.       
          !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  4 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  9 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          j_v =  3 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.

          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  5;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  3 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  8 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  7 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.       
          !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  4 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  9 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !    
       i_v = 1                   ! Equation for psi'   
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}

       i_v = 6                   ! Equation for E^{\chi}
       j_v = 4;                  HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 5;                  HaveExciScaCoef(i_v,j_v) = .TRUE.
       j_v = 6;                  HaveExciScaCoef(i_v,j_v) = .TRUE.

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 0, NSMAX-1
          
          vOffsetA = 10*i_s + NFMAX
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  2 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  3 + NFMAX;      HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4 + NFMAX;      HaveExciScaCoef(i_v,j_v) = .TRUE.       
          !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  4 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  9 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  3 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  8 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO

          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  5;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          j_v =  3 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.

          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  7 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.       
          !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  4 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  9 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO
       
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          DO j_s = 0, NSMAX-1
             vOffsetB = 10*j_s + NFMAX
             j_v =  3 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
             j_v =  8 + vOffsetB;HaveExciScaCoef(i_v,j_v) = .TRUE.
          ENDDO
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  4;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  5;              HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          j_v =  8 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v =  9 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
          j_v = 10 + vOffsetA;   HaveExciScaCoef(i_v,j_v) = .TRUE.
       ENDDO
    END SELECT
       
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EV_COEF_EB

    USE T2COMM, ONLY: NSMAX,NVMAX,NFMAX,NKMAX,&
         &            HaveExciVecCoef,HaveExciVecKval,EqSet
    
    INTEGER(ikind)::&
         & i_s,i_v,i_k,vOffsetA,kOffsetX,&
         &     j_v,    vOffsetB
    
    ! initialization
    HaveExciVecCoef(        1:NVMAX,1:NVMAX) = .FALSE.
    HaveExciVecKVal(1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
    
    SELECT CASE (EqSet)
    CASE (1)
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
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 0, NSMAX-1
       
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          kOffsetX =  2*i_s
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB;   HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  3;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
       
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO

    CASE (2)
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 0, NSMAX-1
       
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          kOffsetX =  2*i_s
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB;   HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  3;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4+kOffsetX;     HaveExciVecKval(i_k,i_v,j_v) = .TRUE.

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  4;              HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  3 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          i_k =  4 + kOffsetX;   HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  3 + vOffsetB;   HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          j_v =  8 + vOffsetB;   HaveExciVecCoef(i_v,j_v)     = .TRUE.
          i_k =  1;              HaveExciVecKval(i_k,i_v,j_v) = .TRUE.
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
       ENDDO
    END SELECT
    
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
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-06-22 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ET_COEF_EB
    
    USE T2COMM, ONLY: NSMAX,NVMAX,NKMAX,NFMAX,&
         &            HaveExciTenCoef,HaveExciTenKval,EqSet
    
    INTEGER(ikind)::&
         & i_s,i_v,i_k,kOffsetX,vOffsetA,&
         &     j_v,j_k,kOffsetY,vOffsetB

    ! initialization
    HaveExciTenCoef(                1:NVMAX,1:NVMAX) = .FALSE.
    HaveExciTenKval(1:NKMAX,1:NKMAX,1:NVMAX,1:NVMAX) = .FALSE.
    
    SELECT CASE (EqSet)

    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 0, NSMAX-1
       
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          kOffsetX =  2*i_s
          kOffsetY = kOffsetX
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA ! Equation for Gamma_{a\para}
          j_v =  4 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k = 1
          j_k = 1;            HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k = 1
          j_k = 1;            HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  9 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k = 1
          j_k = 1;            HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v = 10 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k = 1
          j_k = 1;            HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

          i_v =  4 + vOffsetA ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA ! Equation for p_{a}
          j_v =  4 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  5 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  9 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v = 10 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          
          i_v =  7 + vOffsetA ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA ! Equation for Q_{a\para}
          j_v =  4 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  5 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  9 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v = 10 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

          i_v =  9 + vOffsetA ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA ! Equation for Q_{a}^{\chi}
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 0, NSMAX-1
       
          vOffsetA = 10*i_s + NFMAX
          vOffsetB = vOffsetA
          kOffsetX =  2*i_s
          kOffsetY = kOffsetX
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA ! Equation for Gamma_{a\para}
          j_v =  4 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE.
          j_v =  5 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  9 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v = 10 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

          i_v =  4 + vOffsetA ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA ! Equation for p_{a}
          j_v =  4 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  5 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  9 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v = 10 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          i_k =  3 + kOffsetX
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          
          i_v =  7 + vOffsetA ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA ! Equation for Q_{a\para}
          j_v =  4 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  5 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v =  9 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 
          j_v = 10 + vOffsetA;HaveExciTenCoef(i_v,j_v) = .TRUE.
          i_k =  1
          j_k =  1;           HaveExciTenKval(i_k,j_k,i_v,j_v) = .TRUE. 

          i_v =  9 + vOffsetA ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA ! Equation for Q_{a}^{\chi}
       ENDDO
    END SELECT
    
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
  !                     LAST UPDATE     2014-06-12 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_LOCK_EQUA_EB
    
    USE T2COMM, ONLY:&
         & NSMAX,NVMAX,NFMAX,                      &
         & SolveElectron,SolveIons,                &
         & SolveBp,SolveBt,SolveEt,SolveEp,SolveEr,&
         & SolveNn,SolveFr,SolveFb,SolveFt,SolveFp,&
         & SolvePp,SolveQr,SolveQb,SolveQt,SolveQp,& 
         & LockEqs,EqSet
    
    INTEGER(ikind)::i_s,vOffsetA
    
    ! initialization
    LockEqs(1:NVMAX) = .TRUE.
    
    SELECT CASE (EqSet)
    CASE (1)
       ! field
       LockEqs( 1)    = .NOT.SolveBp
       LockEqs( 2)    = .NOT.SolveBt
       LockEqs( 3)    = .NOT.SolveEt
       LockEqs( 4)    = .NOT.SolveEp
       LockEqs( 5)    = .NOT.SolveEr
    
       
       IF(SolveElectron)THEN
          LockEqs( 1 + NFMAX) = .NOT.SolveNn
          LockEqs( 2 + NFMAX) = .NOT.SolveFr
          LockEqs( 3 + NFMAX) = .NOT.SolveFb
          LockEqs( 4 + NFMAX) = .NOT.SolveFt
          LockEqs( 5 + NFMAX) = .NOT.SolveFp
          LockEqs( 6 + NFMAX) = .NOT.SolvePp
          LockEqs( 7 + NFMAX) = .NOT.SolveQr
          LockEqs( 8 + NFMAX) = .NOT.SolveQb
          LockEqs( 9 + NFMAX) = .NOT.SolveQt
          LockEqs(10 + NFMAX) = .NOT.SolveQp
       ENDIF
    
       IF(SolveIons)THEN
          DO i_s = 1, NSMAX-1
             vOffsetA = 10*i_s + NFMAX
             LockEqs( 1 + vOffsetA) = .NOT.SolveNn
             LockEqs( 2 + vOffsetA) = .NOT.SolveFr
             LockEqs( 3 + vOffsetA) = .NOT.SolveFb
             LockEqs( 4 + vOffsetA) = .NOT.SolveFt
             LockEqs( 5 + vOffsetA) = .NOT.SolveFp
             LockEqs( 6 + vOffsetA) = .NOT.SolvePp
             LockEqs( 7 + vOffsetA) = .NOT.SolveQr
             LockEqs( 8 + vOffsetA) = .NOT.SolveQb
             LockEqs( 9 + vOffsetA) = .NOT.SolveQt
             LockEqs(10 + vOffsetA) = .NOT.SolveQp
          ENDDO
       ENDIF
    CASE (2)
       ! field
       LockEqs( 1)    = .NOT.SolveBp
       LockEqs( 2)    = .NOT.SolveBt
       LockEqs( 3)    = .NOT.SolveEt
       LockEqs( 4)    = .NOT.SolveEp
       LockEqs( 5)    = .NOT.SolveEr
       LockEqs( 6)    = .NOT.SolveEp
       
       IF(SolveElectron)THEN
          LockEqs( 1 + NFMAX) = .NOT.SolveNn
          LockEqs( 2 + NFMAX) = .NOT.SolveFr
          LockEqs( 3 + NFMAX) = .NOT.SolveFb
          LockEqs( 4 + NFMAX) = .NOT.SolveFt
          LockEqs( 5 + NFMAX) = .NOT.SolveFp
          LockEqs( 6 + NFMAX) = .NOT.SolvePp
          LockEqs( 7 + NFMAX) = .NOT.SolveQr
          LockEqs( 8 + NFMAX) = .NOT.SolveQb
          LockEqs( 9 + NFMAX) = .NOT.SolveQt
          LockEqs(10 + NFMAX) = .NOT.SolveQp
       ENDIF
    
       IF(SolveIons)THEN
          DO i_s = 1, NSMAX-1
             vOffsetA = 10*i_s + NFMAX
             LockEqs( 1 + vOffsetA) = .NOT.SolveNn
             LockEqs( 2 + vOffsetA) = .NOT.SolveFr
             LockEqs( 3 + vOffsetA) = .NOT.SolveFb
             LockEqs( 4 + vOffsetA) = .NOT.SolveFt
             LockEqs( 5 + vOffsetA) = .NOT.SolveFp
             LockEqs( 6 + vOffsetA) = .NOT.SolvePp
             LockEqs( 7 + vOffsetA) = .NOT.SolveQr
             LockEqs( 8 + vOffsetA) = .NOT.SolveQb
             LockEqs( 9 + vOffsetA) = .NOT.SolveQt
             LockEqs(10 + vOffsetA) = .NOT.SolveQp
          ENDDO
       ENDIF
    END SELECT

    RETURN
    
  END SUBROUTINE T2VGRA_LOCK_EQUA_EB
      
  !-------------------------------------------------------------------
  ! 
  ! CREATE TABLE ABOUT WHICH VARIABLES ARE TO BE LOCKED ON AXIS
  !
  !
  !                     LAST UPDATE     2014-06-22 H.Seto
  !                     OPERATION CHECK 2014-05-26 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2VGRA_LOCK_AXIS_EB
    
    USE T2COMM, ONLY:&
         & NSMAX,NVMAX,NFMAX,LockAxi,LockEqs,&
         !
         & LockBpOnAxis,LockBtOnAxis,LockEtOnAxis,LockEpOnAxis,&
         & LockErOnAxis,&
         & LockNnOnAxis,LockFrOnAxis,LockFbOnAxis,LockFtOnAxis,&
         & LockFpOnAxis,&
         & LockPpOnAxis,LockQrOnAxis,LockQbOnAxis,LockQtOnAxis,&
         & LockQpOnAxis
    
    INTEGER(ikind)::i_s,vOffsetA
    
    ! initialization
    LockAxi(1:NVMAX) = .FALSE.    
    
    IF(.NOT.LockEqs(1))&
         &  LockAxi(1) = LockBpOnAxis
    IF(.NOT.LockEqs(2))&
         &  LockAxi(2) = LockBtOnAxis
    IF(.NOT.LockEqs(3))&
         &  LockAxi(3) = LockEtOnAxis
    IF(.NOT.LockEqs(4))&
         &  LockAxi(4) = LockEpOnAxis
    IF(.NOT.LockEqs(5))&
         &  LockAxi(5) = LockErOnAxis
    
    DO i_s = 0, NSMAX-1
       vOffsetA = 10*i_s + NFMAX
       IF(.NOT.LockEqs( 1 + vOffsetA))&
            &  LockAxi( 1 + vOffsetA) = LockNnOnAxis
       IF(.NOT.LockEqs( 2 + vOffsetA))&
            &  LockAxi( 2 + vOffsetA) = LockFrOnAxis
       IF(.NOT.LockEqs( 3 + vOffsetA))&
            &  LockAxi( 3 + vOffsetA) = LockFbOnAxis
       IF(.NOT.LockEqs( 4 + vOffsetA))&
            &  LockAxi( 4 + vOffsetA) = LockFtOnAxis
       IF(.NOT.LockEqs( 5 + vOffsetA))&
            &  LockAxi( 5 + vOffsetA) = LockFpOnAxis 
       IF(.NOT.LockEqs( 6 + vOffsetA))&
            &  LockAxi( 6 + vOffsetA) = LockPpOnAxis 
       IF(.NOT.LockEqs( 7 + vOffsetA))&
            &  LockAxi( 7 + vOffsetA) = LockQrOnAxis 
       IF(.NOT.LockEqs( 8 + vOffsetA))&
            &  LockAxi( 8 + vOffsetA) = LockQbOnAxis
       IF(.NOT.LockEqs( 9 + vOffsetA))&
            &  LockAxi( 9 + vOffsetA) = LockQtOnAxis
       IF(.NOT.LockEqs(10 + vOffsetA))&
            &  LockAxi(10 + vOffsetA) = LockQpOnAxis
    ENDDO
    
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
    
    USE T2COMM,ONLY:&
         & NSMAX,NVMAX,NFMAX,LockWal,LockEqs,&
         !
         !
         & LockBpOnWall,LockBtOnWall,LockEtOnWall,LockEpOnWall,&
         & LockErOnWall,&
         & LockNnOnWall,LockFrOnWall,LockFbOnWall,LockFtOnWall,&
         & LockFpOnWall,&
         & LockPpOnWall,LockQrOnWall,LockQbOnWall,LockQtOnWall,&
         & LockQpOnWall

    INTEGER(ikind)::i_s,vOffsetA

    ! initialization
    LockWal(1:NVMAX) = .FALSE.    
    
    
    IF(.NOT.LockEqs( 1))&
         &  LockWal( 1) = LockBpOnWall
    IF(.NOT.LockEqs( 2))&
         &  LockWal( 2) = LockBtOnWall
    IF(.NOT.LockEqs( 3))&
         &  LockWal( 3) = LockEtOnWall
    IF(.NOT.LockEqs( 4))&
         &  LockWal( 4) = LockEpOnWall
    IF(.NOT.LockEqs( 5))&
         &  LockWal( 5) = LockEtOnWall
    

    DO i_s = 0,NSMAX-1
       vOffsetA = 10*i_s + NFMAX
       IF(.NOT.LockEqs( 1 + vOffsetA))&
            &  LockWal( 1 + vOffsetA) = LockNnOnWall
       IF(.NOT.LockEqs( 2 + vOffsetA))&
            &  LockWal( 2 + vOffsetA) = LockFrOnWall
       IF(.NOT.LockEqs( 3 + vOffsetA))&
            &  LockWal( 3 + vOffsetA) = LockFbOnWall
       IF(.NOT.LockEqs( 4 + vOffsetA))&
            &  LockWal( 4 + vOffsetA) = LockFtOnWall
       IF(.NOT.LockEqs( 5 + vOffsetA))&
            &  LockWal( 5 + vOffsetA) = LockFpOnWall
       IF(.NOT.LockEqs( 6 + vOffsetA))&
            &  LockWal( 6 + vOffsetA) = LockPpOnWall 
       IF(.NOT.LockEqs( 7 + vOffsetA))&
            &  LockWal( 7 + vOffsetA) = LockQrOnWall 
       IF(.NOT.LockEqs( 8 + vOffsetA))&
            &  LockWal( 8 + vOffsetA) = LockQbOnWall
       IF(.NOT.LockEqs( 9 + vOffsetA))&
            &  LockWal( 9 + vOffsetA) = LockQtOnWall
       IF(.NOT.LockEqs(10 + vOffsetA))&
            &  LockWal(10 + vOffsetA) = LockQpOnWall
    ENDDO
    
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

    OPEN(10,FILE='T2VGRA_HMMC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'MC=',HaveMat(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HMSC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'MS=',HaveMassScaCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HAVC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'AV=',HaveAdveVecCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HATC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'AT=',HaveAdveTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HDTC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'DT=',HaveDiffTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HGVC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'GV=',HaveGradVecCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HGTC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'GT=',HaveGradTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HESC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'ES=',HaveExciScaCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HEVC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'EV=',HaveExciVecCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HETC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'ET=',HaveExciTenCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HSSC.dat')
    DO i_v = 1, NVMAX
    DO j_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'jv=',j_v,'SS=',HaveSourScaCoef(i_v,j_v)
    ENDDO
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2VGRA_HATK.dat')
    DO i_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'iv=',i_v,'jv=',j_v,&
               'AT=',HaveAdveTenKval(i_k,i_v,j_v)
       ENDDO
       ENDDO
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2VGRA_HGTK.dat')
    DO i_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'iv=',i_v,'jv=',j_v,&
               'GT=',HaveGradTenKval(i_k,i_v,j_v)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_HEVK.dat')
    DO i_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          WRITE(10,*)'ik=',i_k,'iv=',i_v,'jv=',j_v,&
               'EV=',HaveExciVecKval(i_k,i_v,j_v)
       ENDDO
       ENDDO
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2VGRA_HETK.dat')
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

    OPEN(10,FILE='T2VGRA_LEQ.dat')
    DO i_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'LEQ=',LockEqs(i_v)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_LAX.dat')
    DO i_v = 1, NVMAX
       WRITE(10,*)'iv=',i_v,'LAX=',LockAxi(i_v)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2VGRA_LWL.dat')
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
