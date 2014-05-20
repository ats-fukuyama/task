!C-------------------------------------------------------------------- 
!C MODULE FOR INTERPORATION FUNCTION INTEGRATION 
!C                          BY GAUSSIAN QUADRATURE
!C
!C                  LAST UPDATE 2014-05-20 H.SETO
!C
!C T2INTG requires following variables:
!C
!C        rkind,ikind,abcsArrayXX,wghtArrayXX [from T2CNST]
!C        nnmax,nqmax,ndmax                   [from T2COMM] 
!C
!C T2INTG provides following variables:
!C
!C        massScaIntgPG: INTEGRATION ARRAY FOR MASS SCALAR SUBMATRIX
!C        adveVecIntgPG: INTEGRATION ARRAY FOR ADVE VECTOR SUBMATRIX
!C        adveTenIntgPG: INTEGRATION ARRAY FOR ADVE TENSOR SUBMATRIX
!C        diffTenIntgPG: INTEGRATION ARRAY FOR DIFF TENSOR SUBMATRIX
!C        gradVecIntgPG: INTEGRATION ARRAY FOR GRAD VECTOR SUBMATRIX
!C        gradTenIntgPG: INTEGRATION ARRAY FOR GRAD TENSOR SUBMATRIX
!C        exciScaIntgPG: INTEGRATION ARRAY FOR EXCI SCALAR SUBMATRIX
!C        exciVecIntgPG: INTEGRATION ARRAY FOR EXCI VECTOR SUBMATRIX
!C        exciTenIntgPG: INTEGRATION ARRAY FOR EXCI TENSOR SUBMATRIX
!C        sourScaIntgPG: INTEGRATION ARRAY FOR SOUR SCALAR SUBMATRIX
!C
!C  and subroutines:
!C
!C        T2_INTG
!C        T2_INTG_TERMINATE
!C
!C -------------------------------------------------------------------
MODULE T2INTG
  
  USE T2CNST, ONLY: rkind,ikind
  USE T2COMM, ONLY: nnmax,nqmax,ndmax

  IMPLICIT NONE

  PRIVATE  

  REAL(   rkind),ALLOCATABLE,SAVE::&
       MassScaIntgPG(:,:,:        ), AdveVecIntgPG(:,:,:,:      ),&
       AdveTenIntgPG(:,:,:,:,:,:  ), DiffTenIntgPG(:,:,:,:,:    ),&
       GradVecIntgPG(:,:,:,:      ), GradTenIntgPG(:,:,:,:,:,:  ),&
       ExciScaIntgPG(:,:,:        ), ExciVecIntgPG(:,:,:,:,:    ),&
       ExciTenIntgPG(:,:,:,:,:,:,:), SourScaIntgPG(:,:          )

  !C
  !C WORKING ARRAY FOR GAUSSIAN INTEGRATION
  !C
  
  REAL(rkind),DIMENSION(:,:    ),ALLOCATABLE,SAVE:: wghtArray
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE,SAVE:: intgArray

  PUBLIC &
       massScaIntgPG, adveVecIntgPG, adveTenIntgPG, diffTenIntgPG,&
       gradVecIntgPG, gradTenIntgPG, exciScaIntgPG, exciVecIntgPG,&
       exciTenIntgPG, sourScaIntgPG,&
       T2_INTG,&
       T2_INTG_TERMINATE

CONTAINS
  !C-------------------------------------------------------------
  !C
  !C T2INTG_INITIALIZE (PUBLIC)
  !C
  !C                2014-05-20 H.Seto
  !C
  !C-------------------------------------------------------------
  SUBROUTINE T2_INTG
    
    INTEGER(ikind)::&
         i_n,j_n,k_n,l_n,m_n,&
         i_d,j_d,&
         i_q,j_q
    REAL(   rkind)::&
         intgNi,intgNj,intgNk,intgNl,intgNm,&
         weight, sumGaussQuad
    !C------------------------------------------------------
    
    CALL T2INTG_PUBLIC_ALLOCATE

    CALL T2INTG_SETUP_WORKING_ARRAYS
    
    
    !C
    !C INTEGRAL ARRAYS FOR PG-FEM
    !C
    
    !C for mass scalar term
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax    
    DO k_n = 1, nnmax
       sumGaussQuad = 0.D0
       DO i_q = 1, nqmax
       DO j_q = 1, nqmax
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
    
    !C for advection vector term

    DO j_n = 1, nnmax
    DO i_n = 1, nnmax    
    DO k_n = 1, nnmax
       DO i_d = 1, ndmax
          sumGaussQuad =  0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
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

    !C for advection tensor term
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO l_n = 1, nnmax
    DO k_n = 1, nnmax
       DO j_d = 1, ndmax
       DO i_d = 1, ndmax
          sumGaussQuad = 0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
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
    
    !C for diffusion tensor term
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO k_n = 1, nnmax
       DO j_d = 1, ndmax
       DO i_d = 1, ndmax
          sumGaussQuad = 0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
             intgNi = intgArray(i_q,j_q,i_d,i_n) 
             intgNj = intgArray(i_q,j_q,j_d,j_n) 
             intgNk = intgArray(i_q,j_q,0,  k_n) 
             weight = wghtArray(i_q,j_q        )
             sumGaussQuad = sumGaussQuad + intgNi*intgNj*intgNk*weight
          ENDDO
          ENDDO
          diffTenIntgPG(i_d,j_d,k_n,i_n,j_n) = sumGaussQuad
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C for gradient vector term
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO k_n = 1, nnmax
       DO i_d = 1, ndmax
          sumGaussQuad = 0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
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

    !C for gradient tensor term
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO l_n = 1, nnmax
    DO k_n = 1, nnmax    
       DO j_d = 1, ndmax
       DO i_d = 1, ndmax
          sumGaussQuad = 0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
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

    !C for excitation scalar term

    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO k_n = 1, nnmax
       sumGaussQuad = 0.D0
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax
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
    
    !C for excitation vecor term

    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO l_n = 1, nnmax
    DO k_n = 1, nnmax
       DO i_d = 1, ndmax
          sumGaussQuad = 0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
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

    !C for excitation tensor term

    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
    DO m_n = 1, nnmax
    DO l_n = 1, nnmax
    DO k_n = 1, nnmax
       DO j_d = 1, ndmax
       DO i_d = 1, ndmax
          sumGaussQuad = 0.D0
          DO j_q = 1, nqmax
          DO i_q = 1, nqmax
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

    !C for source scalar term
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
       sumGaussQuad = 0.D0
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax
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
    
  END SUBROUTINE T2_INTG
  
  !C-------------------------------------------------------------
  !C
  !C T2_INTG_TERMINATE (PUBLIC)
  !C
  !C                2014-05-20 H.Seto
  !C
  !C-------------------------------------------------------------
  SUBROUTINE T2_INTG_TERMINATE
    CALL T2INTG_PUBLIC_DEALLOCATE
    RETURN
  END SUBROUTINE T2_INTG_TERMINATE

  !C-------------------------------------------------------------
  !C
  !C T2INTG_PUBLIC_ALLOCATE (PRIVATE)
  !C
  !C                2014-05-20 H.Seto 
  !C
  !C-------------------------------------------------------------
  SUBROUTINE T2INTG_PUBLIC_ALLOCATE
    
    INTEGER(ikind),SAVE::&
         nnmax_save=0,ndmax_save=0,nqmax_save=0
    
    INTEGER(ikind):: i0err
    
    IF(  (nnmax .NE.nnmax_save ).OR.&
         (ndmax .NE.ndmax_save ).OR.&
         (nqmax .NE.nqmax_save ))THEN
       
       CALL T2INTG_PUBLIC_DEALLOCATE
       
       DO
          ! for PG-FEM
          ALLOCATE(massScaIntgPG(1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(adveVecIntgPG(1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(adveTenIntgPG(1:ndmax,1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(diffTenIntgPG(1:ndmax,1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(gradVecIntgPG(1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(gradTenIntgPG(1:ndmax,1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(exciScaIntgPG(1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(exciVecIntgPG(1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(exciTenIntgPG(1:ndmax,1:ndmax,&
               &                 1:nnmax,1:nnmax,1:nnmax,1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          ALLOCATE(sourScaIntgPG(1:nnmax,1:nnmax),&
               STAT=i0err); IF(i0err.NE.0) EXIT
          
          !C
          !C for SUPG-FEM
          !C

          !ALLOCATE(massScaIntgSUPG(1:ndmax,&
          !     &                   1:nnmax,1:nnmax,1:nnmax,1:nnmax),&
          !     STAT=i0err); IF(i0err.NE.0) EXIT
          !ALLOCATE(adveVecIntgSUPG(1:ndmax,1:ndmax,&
          !     &                   1:nnmax,1:nnmax,1:nnmax,1:nnmax),&
          !     STAT=i0err); IF(i0err.NE.0) EXIT
          !ALLOCATE(sourScaIntgSUPG(1:ndmax,&
          !     &                   1:nnmax,1:nnmax,1:nnmax),&
          !     STAT=i0err); IF(i0err.NE.0) EXIT
          
          nnmax_save = nnmax
          ndmax_save = ndmax
          nqmax_save = nqmax
          
          WRITE(6,'(A)') '-- T2INTG_PUBLIC_ALLOCATE: SUCCESSED'
          
          RETURN
          
       ENDDO
       
       WRITE(6,'(A)') 'XX T2COMM_ALLOCATE: ALLOCATION ERROR: ECODE=',i0err
       STOP
       
    END IF
    
    RETURN
    
  END SUBROUTINE T2INTG_PUBLIC_ALLOCATE
  
  !C-------------------------------------------------------------
  !C
  !C T2INTG_PUBLIC_DEALLOCATE (PRIVATE)
  !C
  !C                2014-05-20 H.Seto 
  !C
  !C-------------------------------------------------------------
  SUBROUTINE T2INTG_PUBLIC_DEALLOCATE

    !C for PG-FEM

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

  !C------------------------------------------------------------------
  !C
  !C T2INTG_SETUP_WORKING_ARRAYS (PRIVATE)
  !C 
  !C                2014-05-20 H.Seto
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2INTG_SETUP_WORKING_ARRAYS
    
    USE T2CNST, ONLY: abscArray32,wghtArray32
    USE T2COMM, ONLY: nnmax,nqmax,ndmax
    
    INTEGER(ikind)::&
         i_q,i_q_odd,i_q_eve,&
         j_q,j_q_odd,j_q_eve
    
    REAL(   rkind):: abscissa,weight,x,y,abscArray(1:nqmax)

    ! allocate and initialize working arrays
    
    ALLOCATE(wghtArray(1:nqmax,1:nqmax),&
         &   intgArray(1:nqmax,1:nqmax,1:ndmax,1:nnmax))
    
    abscArray(1:nqmax) = 0.D0
    wghtArray(1:nqmax,1:nqmax) = 0.D0
    intgArray(1:nqmax,1:nqmax,1:ndmax,1:nnmax) = 0.D0
    
    
    !C
    !C SET NUMBER OF ABSCISSAS FOR GAUSS INTEGARATION
    !C
    
    SELECT CASE (nqmax)
       
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
    
    !C
    !C SET INTERPOLATION FUNCTION 
    !C

    SELECT CASE (nnmax)
       
    CASE(4)
       !C FOR 4-NODE LINEAR LAGANGIAN RECTANGULAR ELEMENT 
       !C
       !C     y
       !C     ^
       !C     |  4-----3
       !C     |  |     |
       !C     |  |     |
       !C     |  1-----2
       !C   --+------------> x
       !C     |
       !C
       !C

       !C PHI (x_i,y_i)
       !C
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax

          x = abscArray(i_q)
          y = abscArray(j_q)
          
          intgArray(i_q,j_q,0,1) = (1.D0-x)*(1.D0-y)/4.D0
          intgArray(i_q,j_q,0,2) = (1.D0+x)*(1.D0-y)/4.D0
          intgArray(i_q,j_q,0,3) = (1.D0+x)*(1.D0+y)/4.D0
          intgArray(i_q,j_q,0,4) = (1.D0-x)*(1.D0+y)/4.D0
          
       ENDDO
       ENDDO
       
       !C
       !C dPHI/dx (x_i,y_i)
       !
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax
          
          y = abscArray(j_q)
          
          intgArray(i_q,j_q,1,1) = -(1.D0-y)/4.D0 
          intgArray(i_q,j_q,1,2) =  (1.D0-y)/4.D0 
          intgArray(i_q,j_q,1,3) =  (1.D0+y)/4.D0
          intgArray(i_q,j_q,1,4) = -(1.D0+y)/4.D0

       ENDDO
       ENDDO

       !C
       !C dPHI/dy (x_i,y_i)
       !C

       DO j_q = 1, nqmax
       DO i_q = 1, nqmax

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

       !C PHI(x_i,y_i) 
       !C
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax
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
 
       !C
       !C dPHI/dx (x_i,y_i)
       !C
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax

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
       
       !C
       !C dPHI/dy (x_i,y_i)
       !C
       
       DO j_q = 1, nqmax
       DO i_q = 1, nqmax
          
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

    !CASE(12)
       !C FOR 12-NODE CUBIC SERENDIPITY RECTANGULAR ELEMENT
       !C
       !C    10--09--08--07
       !C     |           |
       !C    11          06
       !C     |           |
       !C    12          05
       !C     |           |
       !C    01--02--03--04
       !C
       
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
  SUBROUTINE T2INITG_TERMINATE_WORKING_ARRAYS
    
    IF(ALLOCATED(wghtArray)) DEALLOCATE(wghtArray)
    IF(ALLOCATED(intgArray)) DEALLOCATE(intgArray)
    
    RETURN

  END SUBROUTINE T2INITG_TERMINATE_WORKING_ARRAYS
  
END MODULE T2INTG
