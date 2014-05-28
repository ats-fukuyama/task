!-------------------------------------------------------------------- 
!
!    MODULE FOR INTERPORATION FUNCTION INTEGRATION 
!                             BY GAUSSIAN QUADRATURE
!
!                  LAST UPDATE 2014-05-27 H.Seto
!
!   T2INTG requires following variables:
!
!   [from T2CNST]
!
!        rkind,ikind,AbcsArrayXX,WghtArrayXX 
!
!   [from T2COMM] 
!
!        NNMAX: NUMBER OF NODES PAR ELEMENT
!        NQMAX: NUMBER OF FOR GAUSSIAN-QUAD PAR DIRECTION
!        NDMAX: NUMBER OF DIMENSIONS
!
!        MassScaIntgPG: INTEGRATION ARRAY FOR MASS SCALAR SUBMATRIX
!        AdveVecIntgPG: INTEGRATION ARRAY FOR ADVE VECTOR SUBMATRIX
!        AdveTenIntgPG: INTEGRATION ARRAY FOR ADVE TENSOR SUBMATRIX
!        DiffTenIntgPG: INTEGRATION ARRAY FOR DIFF TENSOR SUBMATRIX
!        GradVecIntgPG: INTEGRATION ARRAY FOR GRAD VECTOR SUBMATRIX
!        GradTenIntgPG: INTEGRATION ARRAY FOR GRAD TENSOR SUBMATRIX
!        ExciScaIntgPG: INTEGRATION ARRAY FOR EXCI SCALAR SUBMATRIX
!        ExciVecIntgPG: INTEGRATION ARRAY FOR EXCI VECTOR SUBMATRIX
!        ExciTenIntgPG: INTEGRATION ARRAY FOR EXCI TENSOR SUBMATRIX
!        SourScaIntgPG: INTEGRATION ARRAY FOR SOUR SCALAR SUBMATRIX
!
!   T2INTG sets up following variables 
!
!        MassScaIntgPG AdveVecIntgPG AdveTenIntgPG DiffTenIntgPG
!        GradVecIntgPG GradTenIntgPG ExciScaIntgPG ExciVecIntgPG
!        ExciTenIntgPG SourScaIntgPG
!
!   through  subroutine: T2INTG_EXECUTE
!
! -------------------------------------------------------------------
MODULE T2INTG
  
  USE T2CNST, ONLY: rkind,ikind
  
  IMPLICIT NONE
  
  PRIVATE  
  
  ! WORKING ARRAY FOR GAUSSIAN INTEGRATION
  
  REAL(rkind),DIMENSION(:,:    ),ALLOCATABLE,SAVE:: wghtArray
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE,SAVE:: intgArray
  
  PUBLIC T2INTG_EXECUTE
  
CONTAINS
  
  !-------------------------------------------------------------
  !
  ! T2INTG_INITIALIZE (PUBLIC)
  !
  !                2014-05-22 H.Seto
  !
  !-------------------------------------------------------------
  SUBROUTINE T2INTG_EXECUTE
    
    USE T2COMM, ONLY: &
         NNMAX,NQMAX,NDMAX,&
         !
         MassScaIntgPG, AdveVecIntgPG, AdveTenIntgPG, DiffTenIntgPG,&
         GradVecIntgPG, GradTenIntgPG, ExciScaIntgPG, ExciVecIntgPG,&
         ExciTenIntgPG, SourScaIntgPG

    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,m_n,&
         i_d,j_d,&
         i_q,j_q
    REAL(   rkind)::&
         intgNi,intgNj,intgNk,intgNl,intgNm,&
         weight, sumGaussQuad
    !C------------------------------------------------------
    
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
       MassScaIntgPG(k_n,i_n,j_n) = sumGaussQuad
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
          AdveVecIntgPG(i_d,k_n,i_n,j_n) = sumGaussQuad
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
          AdveTenIntgPG(i_d,j_d,k_n,l_n,i_n,j_n) = sumGaussQuad
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
          DiffTenIntgPG(i_d,j_d,k_n,i_n,j_n) = sumGaussQuad
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
          GradVecIntgPG(i_d,k_n,i_n,j_n) = sumGaussQuad
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
          GradTenIntgPG(i_d,j_d,k_n,l_n,i_n,j_n) = sumGaussQuad
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
       ExciScaIntgPG(k_n,i_n,j_n) = sumGaussQuad
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
          ExciVecIntgPG(i_d,k_n,l_n,i_n,j_n) = sumGaussQuad
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
          ExciTenIntgPG(i_d,j_d,k_n,l_n,m_n,i_n,j_n) = sumGaussQuad
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
       SourScaIntgPG(i_n,j_n) = sumGaussQuad
    ENDDO
    ENDDO
    
    IF(ALLOCATED(wghtArray)) DEALLOCATE(wghtArray)
    IF(ALLOCATED(intgArray)) DEALLOCATE(intgArray)
    
    RETURN
    
  END SUBROUTINE T2INTG_EXECUTE
  
  !-------------------------------------------------------------------
  !
  ! T2INTG_SETUP_WORKING_ARRAYS (PRIVATE)
  ! 
  !                2014-05-20 H.Seto
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2INTG_SETUP_WORKING_ARRAYS
    
    USE T2CNST, ONLY: abscArray32,wghtArray32
    USE T2COMM, ONLY: NNMAX,NQMAX,NDMAX
    
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
END MODULE T2INTG
