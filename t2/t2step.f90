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
    
    USE T2COMM,ONLY: NXMAX, NVMAX,NRMAX,CoordinateSwitch,&
         &           Xvec, XvecIn, XvecOut,d2xout,&
         &           nconvmax,eps_conv
    
    USE T2COEF,ONLY: T2COEF_EXECUTE
    USE T2EXEC,ONLY: T2EXEC_EXECUTE
    USE T2CONV,ONLY: T2CONV_EXECUTE
    USE T2VOUT,ONLY: T2VOUT_EXECUTE
    
    INTEGER(ikind),INTENT(OUT):: i_conv
    REAL(   rkind),INTENT(OUT):: residual_conv
    INTEGER(ikind):: nconv
    INTEGER(ikind):: i_x,i_v
    CHARACTER(10)::c10nl
    REAL(4):: e0time_0,e0time_1
    
    c10nl='NL'
    
    DO i_x = 1, NXMAX
       DO i_v = 1, NVMAX
          XvecIn(i_v,i_x) = Xvec(i_v,i_x)
       ENDDO
    ENDDO
    
    DO nconv=1,nconvmax
       
       ! CALCULATE PLASMA COEFFICIENTS
       CALL CPU_TIME(e0time_0)
       
       CALL T2COEF_EXECUTE
       
       CALL CPU_TIME(e0time_1)
       WRITE(6,'(A,F10.3,A)') '-- T2_CALV completed:          cpu=', &
            e0time_1-e0time_0,' [s]'
       
       ! ADVECTION DIFFUSION EQUATION SOLVER (SUPG)
       CALL T2EXEC_EXECUTE(NVMAX)
       
       ! CONVERGENCE CHECK
       SELECT CASE (CoordinateSwitch)
       CASE (1)
          CALL T2CONV_EXECUTE(NVMAX,residual_conv)
       CASE (2)
          CALL T2STEP_CONV_TEST(residual_conv)
       END SELECT

       i_conv=nconv
       IF(residual_conv.LE.eps_conv) EXIT
       
       ! UPDATRE VARIABLES FOR NEXT ITERATION 
       DO i_x = 1, NXMAX
          DO i_v = 1, NVMAX
             XvecIn(i_v,i_x) = XvecOut(i_v,i_x)
          ENDDO
       ENDDO
       
    ENDDO
    
    DO i_x = 1, NXMAX
       DO i_v = 1, NVMAX
          Xvec(i_v,i_x) = XvecOut(i_v,i_x)
       ENDDO
    ENDDO
    
    WRITE(6,*)'NLLOOP=',nconv,'EXIT'

    SELECT CASE (CoordinateSwitch)
    CASE (1)
       CALL T2VOUT_EXECUTE
       !d2xout = Xvec
    CASE (2)
       CALL T2STEP_GNUPLOT
    END SELECT

    RETURN
    
  END SUBROUTINE T2_STEP
  
  SUBROUTINE T2STEP_CONV_TEST(residualMax)
    
    USE T2COMM, ONLY:&
         NVMAX,NXMAX,NRMAX,&
         LockEqs,TestCase,&
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

    SELECT CASE (TestCase)

    CASE (1:2)
       ! for 2D dependent variables
       DO i_x = 1, NXMAX
          DO i_v = 1,NVMAX
             IF(LockEqs(i_v))CYCLE
             valOut = XvecOut(i_v,i_x)
             valIn  = XvecIn( i_v,i_x)
             resNumeratorSquared(         i_v)&
                  = resNumeratorSquared(  i_v) + (valOut-valIn)**2
             resDenominatorSquared(i_v)&
                  = resDenominatorSquared(i_v) + valOut**2
          ENDDO
       ENDDO
    CASE (3)
       ! for 1D dependent variables (FSA)
       DO i_r = 1, NRMAX
          i_x = i1mc1d(i_r)
          DO i_v = 1,NVMAX
             IF(LockEqs(i_v))CYCLE
             valOut = XvecOut(i_v,i_x)
             valIn  = XvecIn( i_v,i_x)
             resNumeratorSquared(  i_v)&
                  = resNumeratorSquared(  i_v) + (valOut-valIn)**2
             resDenominatorSquared(i_v)&
                  = resDenominatorSquared(i_v) + valOut**2
          ENDDO
       ENDDO
    END SELECT

    ! CHECK CONVERGENCE
    DO i_v = 1,NVMAX
       IF(.NOT.LockEqs(i_v))THEN
          IF(resDenominatorSquared(i_v).LE.0.D0)THEN

             WRITE(6,'("*********************************************")')
             WRITE(6,'("       ERROR IN T2STEP_CONVERGENCE           ")')
             WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
             WRITE(6,'("*********************************************")')
             WRITE(6,*)i_v,resDenominatorSquared(i_v)
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
    
  END SUBROUTINE T2STEP_CONV_TEST

  SUBROUTINE T2STEP_GNUPLOT
    
    USE T2COMM,ONLY:NMMAX,GlobalCrd,Xvec,i2crt

    INTEGER(ikind)::i_m,i_x
    REAL(   rkind)::x_crd,y_crd,val
    OPEN(40,FILE="gnuplot.dat")
    DO i_m = 1,NMMAX
       i_x   = i2crt(2,i_m)
       x_crd = GlobalCrd(1,i_m)
       y_crd = GlobalCrd(2,i_m)
       val   = Xvec(1,i_x)
       IF(ABS(val).LE.1.D-15) val = 0
       !print'(3(E15.8,1X))',x_crd,y_crd,val
       WRITE(40,'(3(E15.8,1X))')x_crd,y_crd,val
    ENDDO
    CLOSE(40)
    RETURN
  END SUBROUTINE T2STEP_GNUPLOT

END MODULE T2STEP
