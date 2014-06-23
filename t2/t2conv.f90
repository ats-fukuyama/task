MODULE T2CONV
  
  USE T2CNST,ONLY: ikind,rkind
  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE T2CONV_EXECUTE(nvmx_in,residualMax)
    
    USE T2COMM, ONLY:&
         NVMAX,NXMAX,NRMAX,LockEqs,&
         XvecIn,XvecOut,i1mc1d
    
    INTEGER(ikind),INTENT( IN)::nvmx_in
    REAL(   rkind),INTENT(OUT)::residualMax
    
    INTEGER(ikind):: i_x, i_v, i_r
    REAL(   rkind)::&
         valIn,valOut,residual,resDenominator,resNumerator,&
         resNumeratorSquared(  1:NVMAX),&
         resDenominatorSquared(1:NVMAX)
    LOGICAL::&
         convTable(1:NVMAX)
    ! initialization
    residualMax = 0.D0
    resNumeratorSquared(  1:NVMAX) = 0.D0
    resDenominatorSquared(1:NVMAX) = 0.D0
    
    convTable(1:NVMAX) = .TRUE.
    
    DO i_v = 1, nvmx_in
       convTable(i_v) = LockEqs(i_v)
    ENDDO

    ! for 1D dependent variables (FSA)
    DO i_r = 1, NRMAX
       i_x = i1mc1d(i_r)
       DO i_v = 1, 3
          IF(convTable(i_v))CYCLE
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
          IF(convTable(i_v))CYCLE
          valOut = XvecOut(i_v,i_x)
          valIn  = XvecIn( i_v,i_x)
          resNumeratorSquared(         i_v)&
               = resNumeratorSquared(  i_v) + (valOut-valIn)**2
          resDenominatorSquared(i_v)&
               = resDenominatorSquared(i_v) + valOut**2
       ENDDO
    ENDDO
    
    ! CHECK CONVERGENCE
    DO i_v = 1,NVMAX
       IF(.NOT.convTable(i_v))THEN
          IF(resDenominatorSquared(i_v).EQ.0.D0)THEN
             IF(resNumeratorSquared(i_v).EQ.0.D0)THEN
                residual = 0.D0
             ELSE
                WRITE(6,'("*********************************************")')
                WRITE(6,'("       ERROR IN T2STEP_CONVERGENCE           ")')
                WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
                WRITE(6,'("*********************************************")')
                WRITE(6,*)i_v,resDenominatorSquared(i_v)
                STOP
             ENDIF
          END IF
          resDenominator = SQRT(resDenominatorSquared(i_v))
          resNumerator   = SQRT(resNumeratorSquared(  i_v))
          residual       = resNumerator/resDenominator
          WRITE(6,*),'VARIABLES=',i_v,'RESIDUAL=',residual
          residualMax = MAX(residualMax,residual)
       ENDIF
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CONV_EXECUTE
END MODULE T2CONV
