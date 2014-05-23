!-------------------------------------------------------------------- 
!
!   MODULE FOR INTERPORATION FUNCTION INTEGRATION 
!                          BY GAUSSIAN QUADRATURE
!
!                  LAST UPDATE 2014-05-20 H.SETO
!
!   T2INTG requires following variables:
!
!   [from T2CNST]
!        rkind,ikind,abcsArrayXX,wghtArrayXX 
!   [from T2COMM] 
!        NNMAX,NQMAX,NDMAX                   
!
!   T2INTG provides following variables:
!
!        massScaIntgPG: INTEGRATION ARRAY FOR MASS SCALAR SUBMATRIX
!        adveVecIntgPG: INTEGRATION ARRAY FOR ADVE VECTOR SUBMATRIX
!        adveTenIntgPG: INTEGRATION ARRAY FOR ADVE TENSOR SUBMATRIX
!        diffTenIntgPG: INTEGRATION ARRAY FOR DIFF TENSOR SUBMATRIX
!        gradVecIntgPG: INTEGRATION ARRAY FOR GRAD VECTOR SUBMATRIX
!        gradTenIntgPG: INTEGRATION ARRAY FOR GRAD TENSOR SUBMATRIX
!        exciScaIntgPG: INTEGRATION ARRAY FOR EXCI SCALAR SUBMATRIX
!        exciVecIntgPG: INTEGRATION ARRAY FOR EXCI VECTOR SUBMATRIX
!        exciTenIntgPG: INTEGRATION ARRAY FOR EXCI TENSOR SUBMATRIX
!        sourScaIntgPG: INTEGRATION ARRAY FOR SOUR SCALAR SUBMATRIX
!
!   and subroutines:
!
!        T2INTG_EXECUTE
!        T2INTG_TERMINATE
!
! -------------------------------------------------------------------
MODULE T2INTG
  
  USE T2CNST, ONLY: rkind,ikind
  !USE T2COMM, ONLY: NNMAX,NQMAX,NDMAX

  IMPLICIT NONE

  PRIVATE  
  
  REAL(   rkind),ALLOCATABLE,SAVE::&
       MassScaIntgPG(:,:,:        ), AdveVecIntgPG(:,:,:,:      ),&
       AdveTenIntgPG(:,:,:,:,:,:  ), DiffTenIntgPG(:,:,:,:,:    ),&
       GradVecIntgPG(:,:,:,:      ), GradTenIntgPG(:,:,:,:,:,:  ),&
       ExciScaIntgPG(:,:,:        ), ExciVecIntgPG(:,:,:,:,:    ),&
       ExciTenIntgPG(:,:,:,:,:,:,:), SourScaIntgPG(:,:          )


  ! WORKING ARRAY FOR GAUSSIAN INTEGRATION

  
  REAL(rkind),DIMENSION(:,:    ),ALLOCATABLE,SAVE:: wghtArray
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE,SAVE:: intgArray
  
  ! T2INTG_ADHOC
  INTEGER(ikind),SAVE::NNMAX,NQMAX,NDMAX
  
  PUBLIC &
       massScaIntgPG, adveVecIntgPG, adveTenIntgPG, diffTenIntgPG,&
       gradVecIntgPG, gradTenIntgPG, exciScaIntgPG, exciVecIntgPG,&
       exciTenIntgPG, sourScaIntgPG,&
       T2INTG_EXECUTE,&
       T2INTG_TERMINATE
  
CONTAINS

  !-------------------------------------------------------------
  !
  ! T2INTG_ADHOC (PRIVATE)
  !
  !                2014-05-22 H.Seto
  !
  !-------------------------------------------------------------
  SUBROUTINE T2INTG_ADHOC
    
    USE T2COMM,ONLY:i0nmax,i0qmax,i0dmax
    
    NNMAX = i0nmax
    NQMAX = i0qmax
    NDMAX = i0dmax
    
    RETURN

  END SUBROUTINE T2INTG_ADHOC
  
  !-------------------------------------------------------------
  !
  ! T2INTG_INITIALIZE (PUBLIC)
  !
  !                2014-05-22 H.Seto
  !
  !-------------------------------------------------------------
  SUBROUTINE T2INTG_EXECUTE
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,m_n,&
         i_d,j_d,&
         i_q,j_q
    REAL(   rkind)::&
         intgNi,intgNj,intgNk,intgNl,intgNm,&
         weight, sumGaussQuad
    !C------------------------------------------------------
    
    CALL T2INTG_ADHOC
    
    CALL T2INTG_PUBLIC_ALLOCATE
    
    CALL T2INTG_SETUP_WORKING_ARRAYS

    !
    ! INTEGRAL ARRAYS FOR PG-FEM
    !
    
    ! for mass scalar term
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX    
    DO k_n = 1, NNMAX
       sumGaussQuad = 0.D0
       DO i_q = 1, NQMAX
       DO j_q = 1, NQMAX
          intgNi = intgArray(i_q,j_q,0,i_n) 
          intgNj = intgArray(i_q,j_q,0,j_n) 
          intgNk = intgArray(i_q,j_q,0,k_n) 
          weight = wghtArray(i_q,j_q      )
          sumGaussQuad  = sumGaussQuad &
               &        + intgNi*intgNj*intgNk*weight
       ENDDO
       ENDDO
       massScaIntgPG(k_n,i_n,j_n) = sumGaussQuad
    ENDDO
    ENDDO
    ENDDO
    
    ! for advection vector term

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX    
    DO k_n = 1, NNMAX
       DO i_d = 1, NDMAX
          sumGaussQuad =  0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,0,  i_n) 
             intgNj = intgArray(i_q,j_q,i_d,j_n) 
             intgNk = intgArray(i_q,j_q,0,  k_n) 
             weight = wghtArray(i_q,j_q        ) 
             sumGaussQuad = sumGaussQuad + intgNi*intgNj*intgNk*weight
             intgNj = intgArray(i_q,j_q,0,  j_n) 
             intgNk = intgArray(i_q,j_q,i_d,k_n) 
             sumGaussQuad = sumGaussQuad + intgNi*intgNj*intgNk*weight
          ENDDO
          ENDDO
          adveVecIntgPG(i_d,k_n,i_n,j_n) = sumGaussQuad
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    ! for advection tensor term
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO l_n = 1, NNMAX
    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          sumGaussQuad = 0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,i_d,i_n)
             intgNj = intgArray(i_q,j_q,0,  j_n)
             intgNk = intgArray(i_q,j_q,j_d,k_n)
             intgNl = intgArray(i_q,j_q,0,  l_n)
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad  = sumGaussQuad &
                  &        + intgNi*intgNj*intgNk*intgNl*weight
          ENDDO
          ENDDO
          adveTenIntgPG(i_d,j_d,k_n,l_n,i_n,j_n) = sumGaussQuad
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    ! for diffusion tensor term
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          sumGaussQuad = 0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,i_d,i_n) 
             intgNj = intgArray(i_q,j_q,j_d,j_n) 
             intgNk = intgArray(i_q,j_q,0,  k_n) 
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad = sumGaussQuad &
                  &       + intgNi*intgNj*intgNk*weight
          ENDDO
          ENDDO
          diffTenIntgPG(i_d,j_d,k_n,i_n,j_n) = sumGaussQuad
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    ! for gradient vector term
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO k_n = 1, NNMAX
       DO i_d = 1, NDMAX
          sumGaussQuad = 0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,0,  i_n) 
             intgNj = intgArray(i_q,j_q,i_d,j_n) 
             intgNk = intgArray(i_q,j_q,0,  k_n) 
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad = sumGaussQuad &
                  &       + intgNi*intgNj*intgNk*weight
          ENDDO
          ENDDO
          gradVecIntgPG(i_d,k_n,i_n,j_n) = sumGaussQuad
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    ! for gradient tensor term
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO l_n = 1, NNMAX
    DO k_n = 1, NNMAX    
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          sumGaussQuad = 0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,0,  i_n)
             intgNj = intgArray(i_q,j_q,j_d,j_n) 
             intgNk = intgArray(i_q,j_q,i_d,k_n)
             intgNl = intgArray(i_q,j_q,0,  l_n)
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad = sumGaussQuad &
                  &       + intgNi*intgNj*intgNk*intgNl*weight
          ENDDO
          ENDDO
          gradTenIntgPG(i_d,j_d,k_n,l_n,i_n,j_n) = sumGaussQuad
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    ! for excitation scalar term

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO k_n = 1, NNMAX
       sumGaussQuad = 0.D0
       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX
          intgNi = intgArray(i_q,j_q,0,i_n)
          intgNj = intgArray(i_q,j_q,0,j_n)
          intgNk = intgArray(i_q,j_q,0,k_n)
          weight = wghtArray(i_q,j_q      )
          sumGaussQuad = sumGaussQuad &
               &       + intgNi*intgNj*intgNk*weight
       ENDDO
       ENDDO
       exciScaIntgPG(k_n,i_n,j_n) = sumGaussQuad
    ENDDO
    ENDDO
    ENDDO
    
    ! for excitation vecor term

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO l_n = 1, NNMAX
    DO k_n = 1, NNMAX
       DO i_d = 1, NDMAX
          sumGaussQuad = 0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,0,  i_n)
             intgNj = intgArray(i_q,j_q,0,  j_n)
             intgNk = intgArray(i_q,j_q,i_d,k_n)
             intgNl = intgArray(i_q,j_q,0,  l_n)
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad = sumGaussQuad &
                  &       + intgNi*intgNj*intgNk*intgNl*weight
          ENDDO
          ENDDO
          exciVecIntgPG(i_d,k_n,l_n,i_n,j_n) = sumGaussQuad
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    ! for excitation tensor term

    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
    DO m_n = 1, NNMAX
    DO l_n = 1, NNMAX
    DO k_n = 1, NNMAX
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          sumGaussQuad = 0.D0
          DO j_q = 1, NQMAX
          DO i_q = 1, NQMAX
             intgNi = intgArray(i_q,j_q,0,  i_n)
             intgNj = intgArray(i_q,j_q,0,  j_n)
             intgNk = intgArray(i_q,j_q,i_d,k_n)
             intgNl = intgArray(i_q,j_q,j_d,l_n)
             intgNm = intgArray(i_q,j_q,0,  m_n)
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad = sumGaussQuad &
                  &       + intgNi*intgNj*intgNk*intgNl*intgNm*weight
          ENDDO
          ENDDO
          exciTenIntgPG(i_d,j_d,k_n,l_n,m_n,i_n,j_n) = sumGaussQuad
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    ! for source scalar term
    
    DO j_n = 1, NNMAX
    DO i_n = 1, NNMAX
       sumGaussQuad = 0.D0
       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX
          intgNi = intgArray(i_q,j_q,0,i_n)
          intgNj = intgArray(i_q,j_q,0,j_n)
          weight = wghtArray(i_q,j_q      )
          sumGaussQuad = sumGaussQuad &
               &       + intgNi*intgNj*weight
       ENDDO
       ENDDO
       sourScaIntgPG(i_n,j_n) = sumGaussQuad
    ENDDO
    ENDDO
    
    CALL T2INTG_TERMINATE_WORKING_ARRAYS
    
    RETURN
    
  END SUBROUTINE T2INTG_EXECUTE
  
  !-------------------------------------------------------------
  !
  !  T2_INTG_TERMINATE (PUBLIC)
  !
  !                2014-05-22 H.Seto
  !
  !-------------------------------------------------------------
  SUBROUTINE T2INTG_TERMINATE
    CALL T2INTG_PUBLIC_DEALLOCATE
    CALL T2INTG_PRIVATE_DEALLOCATE
    RETURN
  END SUBROUTINE T2INTG_TERMINATE

  !-------------------------------------------------------------
  !
  !  T2INTG_PUBLIC_ALLOCATE (PRIVATE)
  !
  !                2014-05-22 H.Seto 
  !
  !-------------------------------------------------------------
  SUBROUTINE T2INTG_PUBLIC_ALLOCATE
    
    INTEGER(ikind),SAVE::&
         NNMAX_save=0,NDMAX_save=0,NQMAX_save=0
    
    INTEGER(ikind):: i0err
    
    IF(  (NNMAX .NE.NNMAX_save ).OR.&
         (NDMAX .NE.NDMAX_save ).OR.&
         (NQMAX .NE.NQMAX_save ))THEN
       
       CALL T2INTG_PUBLIC_DEALLOCATE
       
       DO
          ! for PG-FEM
          ALLOCATE(massScaIntgPG(1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(adveVecIntgPG(1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(adveTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(diffTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(gradVecIntgPG(1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(gradTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(exciScaIntgPG(1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(exciVecIntgPG(1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(exciTenIntgPG(1:NDMAX,1:NDMAX,&
               &                 1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(sourScaIntgPG(1:NNMAX,1:NNMAX),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          
          ! for SUPG-FEM

          !ALLOCATE(massScaIntgSUPG(1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=i0err); IF(i0err.NE.0) EXIT
          !ALLOCATE(adveVecIntgSUPG(1:NDMAX,1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=i0err); IF(i0err.NE.0) EXIT
          !ALLOCATE(sourScaIntgSUPG(1:NDMAX,&
          !     &                   1:NNMAX,1:NNMAX,1:NNMAX),&
          !     STAT=i0err); IF(i0err.NE.0) EXIT
          
          NNMAX_save = NNMAX
          NDMAX_save = NDMAX
          NQMAX_save = NQMAX
          
          WRITE(6,'(A)') '-- T2INTG_PUBLIC_ALLOCATE: SUCCESSED'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)') 'XX T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',i0err
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2INTG_PUBLIC_ALLOCATE
  
  !-------------------------------------------------------------
  !
  !  T2INTG_PUBLIC_DEALLOCATE (PRIVATE)
  !
  !                2014-05-20 H.Seto 
  !
  !-------------------------------------------------------------
  SUBROUTINE T2INTG_PUBLIC_DEALLOCATE

    ! for PG-FEM

    IF(ALLOCATED(massScaIntgPG)) DEALLOCATE(massScaIntgPG)
    IF(ALLOCATED(adveVecIntgPG)) DEALLOCATE(adveVecIntgPG)
    IF(ALLOCATED(adveTenIntgPG)) DEALLOCATE(adveTenIntgPG)
    IF(ALLOCATED(diffTenIntgPG)) DEALLOCATE(diffTenIntgPG)
    IF(ALLOCATED(gradVecIntgPG)) DEALLOCATE(gradVecIntgPG)
    IF(ALLOCATED(gradTenIntgPG)) DEALLOCATE(gradTenIntgPG)
    IF(ALLOCATED(exciScaIntgPG)) DEALLOCATE(exciScaIntgPG)
    IF(ALLOCATED(exciVecIntgPG)) DEALLOCATE(exciVecIntgPG)
    IF(ALLOCATED(exciTenIntgPG)) DEALLOCATE(exciTenIntgPG)
    IF(ALLOCATED(sourScaIntgPG)) DEALLOCATE(sourScaIntgPG)
    
    !C for SUPG-FEM
    
    !IF(ALLOCATED(massScaIntgSUPG)) DEALLOCATE(massScaIntgSUPG)
    !IF(ALLOCATED(adveVecIntgSUPG)) DEALLOCATE(adveVecIntgSUPG)
    !IF(ALLOCATED(sourScaIntgSUPG)) DEALLOCATE(sourScaIntgSUPG)
   
    RETURN

  END SUBROUTINE T2INTG_PUBLIC_DEALLOCATE

  !-------------------------------------------------------------------
  !
  ! T2INTG_SETUP_WORKING_ARRAYS (PRIVATE)
  ! 
  !                2014-05-20 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2INTG_SETUP_WORKING_ARRAYS
    
    USE T2CNST, ONLY: abscArray32,wghtArray32
    !USE T2COMM, ONLY: NNMAX,NQMAX,NDMAX
    
    INTEGER(ikind)::&
         i_q,i_q_odd,i_q_eve,&
         j_q,j_q_odd,j_q_eve
    
    REAL(   rkind):: abscissa,weight,x,y,abscArray(1:NQMAX)

    ! allocate and initialize working arrays
    
    ALLOCATE(wghtArray(1:NQMAX,1:NQMAX),&
         &   intgArray(1:NQMAX,1:NQMAX,1:NDMAX,1:NNMAX))
    
    abscArray(1:NQMAX) = 0.D0
    wghtArray(1:NQMAX,1:NQMAX) = 0.D0
    intgArray(1:NQMAX,1:NQMAX,1:NDMAX,1:NNMAX) = 0.D0
    
    
    !
    ! SET NUMBER OF ABSCISSAS FOR GAUSS INTEGARATION
    !
    
    SELECT CASE (NQMAX)
       
    CASE(32) !C 32*32 POINTS 2D GAUSSIAN INTEGRATION
       
       DO j_q = 1, 16
          
          j_q_odd = 2*j_q - 1
          j_q_eve = 2*j_q
          
          DO i_q = 1, 16
             
             i_q_odd = 2*i_q - 1
             i_q_eve = 2*i_q
             
             weight = wghtArray32(i_q)*wghtArray32(j_q)

             wghtArray(i_q_odd,j_q_odd) = weight
             wghtArray(i_q_eve,j_q_odd) = weight
             wghtArray(i_q_odd,j_q_eve) = weight
             wghtArray(i_q_eve,j_q_eve) = weight
             
          ENDDO
          
          abscissa = abscArray32(j_q)
          
          abscArray(j_q_odd) = - abscissa
          abscArray(j_q_eve) =   abscissa
          
       ENDDO
       
    CASE DEFAULT
       WRITE(6,*)'-------------------------------------------------'
       WRITE(6,*)'SUBROUTINE SET_INTEGRATED_INTERPOLATION FUNCTIONS'
       WRITE(6,*)'ERROR: ILLEGAL NQMAX                             '
       WRITE(6,*)'-------------------------------------------------'
       STOP
    END SELECT
    
    !
    ! SET INTERPOLATION FUNCTION 
    !

    SELECT CASE (NNMAX)
       
    CASE(4)
       !
       !   04-points Bilinear Lagrangian Rectangular Element
       !
       !   (-1, 1)            ( 1, 1)         
       !         04----------03            
       !          |           |            
       !          |           |            y
       !          |           |            ^ 
       !          |           |            | 
       !          |           |            | 
       !         01----------02          --+------> x
       !   (-1,-1)            ( 1,-1)      |
       !                                     

       ! PHI (x_i,y_i)

       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX

          x = abscArray(i_q)
          y = abscArray(j_q)
          
          intgArray(i_q,j_q,0,1) = (1.D0-x)*(1.D0-y)/4.D0
          intgArray(i_q,j_q,0,2) = (1.D0+x)*(1.D0-y)/4.D0
          intgArray(i_q,j_q,0,3) = (1.D0+x)*(1.D0+y)/4.D0
          intgArray(i_q,j_q,0,4) = (1.D0-x)*(1.D0+y)/4.D0
          
       ENDDO
       ENDDO
       
       ! dPHI/dx (x_i,y_i)

       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX
          
          y = abscArray(j_q)
          
          intgArray(i_q,j_q,1,1) = -(1.D0-y)/4.D0 
          intgArray(i_q,j_q,1,2) =  (1.D0-y)/4.D0 
          intgArray(i_q,j_q,1,3) =  (1.D0+y)/4.D0
          intgArray(i_q,j_q,1,4) = -(1.D0+y)/4.D0

       ENDDO
       ENDDO

       ! dPHI/dy (x_i,y_i)

       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX

          x = abscArray(i_q)

          intgArray(i_q,j_q,2,1) = -(1.D0-x)/4.D0 
          intgArray(i_q,j_q,2,2) = -(1.D0+x)/4.D0 
          intgArray(i_q,j_q,2,3) =  (1.D0+x)/4.D0
          intgArray(i_q,j_q,2,4) =  (1.D0-x)/4.D0

       ENDDO
       ENDDO

    CASE(8)! this part will be replaced 
           ! by quadradic lagrangian rectangular element 
           ! for the compatibility with FSA algorithm

       !C FOR 8-NODE QUADRADIC SERENDIPITY RECTANGULAR ELEMENT
       !C
       !C   y
       !C   ^
       !C   |   7---6---5
       !C   |   |       |
       !C   |   8       4
       !C   |   |       |
       !C   |   1---2---3
       !C --+---------------> x
       !C   |
       !C
       !C

       ! F_i(Q_j)  
  
       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX
          x = abscArray(i_q)
          y = abscArray(j_q)
          intgArray(i_q,j_q,0,1) = -(1.D0-x)*(1.D0-y)*(1.D0+x+y)/4.D0
          intgArray(i_q,j_q,0,2) =  (1.D0-x)*(1.D0+x)*(1.D0  -y)/2.D0
          intgArray(i_q,j_q,0,3) = -(1.D0+x)*(1.D0-y)*(1.D0-x+y)/4.D0
          intgArray(i_q,j_q,0,4) =  (1.D0-y)*(1.D0+y)*(1.D0+x  )/2.D0
          intgArray(i_q,j_q,0,5) = -(1.D0+x)*(1.D0+y)*(1.D0-x-y)/4.D0
          intgArray(i_q,j_q,0,6) =  (1.D0-x)*(1.D0+x)*(1.D0  +y)/2.D0
          intgArray(i_q,j_q,0,7) = -(1.D0-x)*(1.D0+y)*(1.D0+x-y)/4.D0
          intgArray(i_q,j_q,0,8) =  (1.D0-y)*(1.D0+y)*(1.D0-x  )/2.D0 
       ENDDO
       ENDDO
 
       
       ! d F_i/dx (Q_j)
       
       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX

          x = abscArray(i_q)
          y = abscArray(j_q)
          
          intgArray(i_q,j_q,1,1) =  (1.D0-y)*(2.D0*x+y)/4.D0
          intgArray(i_q,j_q,1,2) = -x*(1.D0-y)
          intgArray(i_q,j_q,1,3) =  (1.D0-y)*(2.D0*x-y)/4.D0 
          intgArray(i_q,j_q,1,4) =  (1.D0-y**2)/2.D0 
          intgArray(i_q,j_q,1,5) =  (1.D0+y)*(2.D0*x+y)/4.D0
          intgArray(i_q,j_q,1,6) = -x*(1.D0+y)
          intgArray(i_q,j_q,1,7) =  (1.D0+y)*(2.D0*x-y)/4.D0
          intgArray(i_q,j_q,1,8) = -(1.D0-y**2)/2.D0

       ENDDO
       ENDDO
       
       ! d F_i/dy (Q_j)

       DO j_q = 1, NQMAX
       DO i_q = 1, NQMAX
          
          x = abscArray(i_q)
          y = abscArray(j_q)
          
          intgArray(i_q,j_q,2,1) =  (1.D0-x)*(2.D0*y+x)/4.D0
          intgArray(i_q,j_q,2,2) = -(1.D0-x**2)/2.D0
          intgArray(i_q,j_q,2,3) =  (1.D0+x)*(2.D0*y-x)/4.D0
          intgArray(i_q,j_q,2,4) = -y*(1.D0+x)
          intgArray(i_q,j_q,2,5) =  (1.D0+x)*(2.D0*y+x)/4.D0 
          intgArray(i_q,j_q,2,6) =  (1.D0-x**2)/2.D0
          intgArray(i_q,j_q,2,7) =  (1.D0-x)*(2.D0*y-x)/4.D0
          intgArray(i_q,j_q,2,8) = -y*(1.D0-x)

       ENDDO
       ENDDO

       ! dd F_i/dxdx (Q_j)
       ! dd F_i/dxdy (Q_j)
       ! dd F_i/dydy (Q_j)

    !CASE ( 9)
       !
       !   09-points Biquadradic Lagrangian Rectangular Element
       !
       !   (-1, 1)            ( 1, 1)         
       !         04----07----03            
       !          |           |            
       !          |           |            y
       !         08    09    06            ^ 
       !          |           |            | 
       !          |           |            | 
       !         01----05----02          --+------> x
       !   (-1,-1)            ( 1,-1)      |
       !                                     

       !  d F_i/dx   (Q_j)
       !  d F_i/dy   (Q_j)
       ! dd F_i/dxdx (Q_j)
       ! dd F_i/dxdy (Q_j)
       ! dd F_i/dydy (Q_j)

    !CASE(16)
       !
       !   16-points Bicubic Lagrangian Rectangular Element
       !
       !   (-1, 1)            ( 1, 1)         
       !         04--11--07--03            
       !          |           |            
       !         08  16  15  10            y
       !          |           |            ^ 
       !         12  13  14  06            | 
       !          |           |            | 
       !         01--05--09--02          --+------> x
       !   (-1,-1)            ( 1,-1)      |
       !                                     

       !  d F_i/dx   (Q_j)
       !  d F_i/dy   (Q_j)
       ! dd F_i/dxdx (Q_j)
       ! dd F_i/dxdy (Q_j)
       ! dd F_i/dydy (Q_j)
       
    CASE DEFAULT
       
       WRITE(6,*)'-------------------------------------------------'
       WRITE(6,*)'SUBROUTINE SET_INTEGRATED_INTERPOLATION FUNCTIONS'
       WRITE(6,*)'ERROR: ILLEGAL I0NMAX0'
       WRITE(6,*)'-------------------------------------------------'
       STOP

    END SELECT
    
    RETURN
    
  END SUBROUTINE T2INTG_SETUP_WORKING_ARRAYS

  !C------------------------------------------------------------------
  !C
  !C T2INTG_TERMINATE_WORKING_ARRAYS (PRIVATE)
  !C 
  !C                2014-05-20 H.Seto
  !C
  !C------------------------------------------------------------------  
  SUBROUTINE T2INTG_TERMINATE_WORKING_ARRAYS
    
    IF(ALLOCATED(wghtArray)) DEALLOCATE(wghtArray)
    IF(ALLOCATED(intgArray)) DEALLOCATE(intgArray)
    
    RETURN
    
  END SUBROUTINE T2INTG_TERMINATE_WORKING_ARRAYS
  
END MODULE T2INTG
