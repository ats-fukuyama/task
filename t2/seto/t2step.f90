!--------------------------------------------------------------------
!
!    MODULE FOR NONLINEAR ITERATION 
!                             BY PICARD ITETRATION
!
!                   LAST UPDATE 2014-05-28 H.Seto
!
!--------------------------------------------------------------------
MODULE T2STEP
  
  USE T2CNST,ONLY:&
       ikind,rkind
  
  IMPLICIT NONE
  
  PUBLIC T2_STEP
  PRIVATE
  
CONTAINS
  
  SUBROUTINE T2_STEP(i_conv,residual_conv)
    
    USE T2COMM,ONLY: NXMAX, NVMAX,NRMAX,&
         &           XvecIn, XvecOut,&
         &           nconvmax,eps_conv
    
    USE T2CALV,ONLY: T2_CALV
    USE T2EXEC,ONLY: T2EXEC_EXECUTE
    USE T2CONV,ONLY: T2_CONV
    
    INTEGER(ikind),INTENT(OUT):: i_conv
    REAL(   rkind),INTENT(OUT):: residual_conv
    INTEGER(ikind):: nconv
    INTEGER(ikind):: i_x,i_v
    CHARACTER(10)::c10nl
    REAL(4):: e0time_0,e0time_1
    
    c10nl='NL'
    
    DO nconv=1,nconvmax

       !C
       !C CALCULATE PLASMA COEFFICIENTS
       !C
       CALL CPU_TIME(e0time_0)
       CALL T2_CALV
       CALL CPU_TIME(e0time_1)
       WRITE(6,'(A,F10.3,A)') '-- T2_CALV completed:          cpu=', &
                              e0time_1-e0time_0,' [s]'
       !C
       !C ADVECTION DIFFUSION EQUATION SOLVER (SUPG)
       !C
       
       CALL T2EXEC_EXECUTE

       !C 
       !C CONVERGENCE CHECK
       !C
       
       CALL T2STEP_CONV(residual_conv)
       
       ! UPDATRE VARIABLES FOR NEXT ITERATION 
       DO i_x = 1, NXMAX
          DO i_v = 1, NVMAX
             XvecIn(i_v,i_x) = XvecOut(i_v,i_x)
          ENDDO
       ENDDO
       
       i_conv=nconv
       IF(residual_conv.LE.eps_conv) EXIT
       
    ENDDO
    
    WRITE(6,*)'NLLOOP=',nconv,'EXIT'
    
    CALL T2_CONV
    
    RETURN
    
  END SUBROUTINE T2_STEP
  
  SUBROUTINE T2STEP_CONV(residualMax)
    
    USE T2COMM, ONLY:&
         NVMAX,NXMAX,NRMAX,&
         LockEqs,&
         XvecIn,XvecOut,&
         i1mc1d
    
    REAL(   rkind),INTENT(OUT)::residualMax

    INTEGER(ikind):: i_x, i_v, i_r
    REAL(   rkind)::&
         valIn,valOut,residual,resDenominator,resNumerator,&
         resNumeratorSquared(  1:NVMAX),&
         resDenominatorSquared(1:NVMAX)

    ! initialization
    residualMax = 0.D0
    resNumeratorSquared(  1:NVMAX) = 0.D0
    resDenominatorSquared(1:NVMAX) = 0.D0
    
    ! for 1D dependent variables (FSA)
    DO i_r = 1, NRMAX
       i_x = i1mc1d(i_r)
       DO i_v = 1, 3
          valOut = XvecOut(i_v,i_x)
          valIn  = XvecIn( i_v,i_x)
          resNumeratorSquared(  i_v)&
               = resNumeratorSquared(  i_v) + (valOut-valIn)**2
          resDenominatorSquared(i_v)&
               = resDenominatorSquared(i_v) + valOut**2
       ENDDO
    ENDDO
    
    ! for 2D dependent variables
    DO i_x = 1, NXMAX
       DO i_v = 4,NVMAX
          valOut = XvecOut(i_v,i_x)
          valIn  = XvecIn( i_v,i_x)
          resNumeratorSquared(  i_v)&
               = resNumeratorSquared(  i_v) + (valOut-valIn)**2
          resDenominatorSquared(i_v)&
               = resDenominatorSquared(i_v) + valOut**2
       ENDDO
    ENDDO
    
    ! CHECK CONVERGENCE
    DO i_v = 1,NVMAX
       IF(.NOT.LockEqs(i_v))THEN
          IF(resDenominatorSquared(i_v).LE.0)THEN
             WRITE(6,'("*********************************************")')
             WRITE(6,'("       ERROR IN T2STEP_CONVERGENCE           ")')
             WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
             WRITE(6,'("*********************************************")')
             WRITE(6,*)i_v
             STOP
          ENDIF

          resDenominator = SQRT(resDenominatorSquared(i_v))
          resNumerator   = SQRT(resNumeratorSquared(  i_v))
          residual       = resNumerator/resDenominator
          WRITE(6,*),'VARIABLES=',i_v,'RESIDUAL=',residual
          residualMax = MAX(residualMax,residual)
       ENDIF
    ENDDO
   
    RETURN
    
  END SUBROUTINE T2STEP_CONV
END MODULE T2STEP
