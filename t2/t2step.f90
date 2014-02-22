!C--------------------------------------------------------------------
!C
!C T2STEP
!C
!C                       2014-02-22 H.SETO
!C
!C--------------------------------------------------------------------
MODULE T2STEP
  
  USE T2CNST,ONLY:&
       i0ikind,i0rkind
  
  IMPLICIT NONE
  
  PUBLIC T2_STEP
  PRIVATE
  
CONTAINS

  SUBROUTINE T2_STEP(i_conv,residual_conv)
    
    USE T2COMM, ONLY:&
         i0pmax, i0xmax, i0vmax, &
         d0eps,&
         d2xvec, d2xvec_befor, d2xvec_after, &
         nconvmax,eps_conv
    
    USE T2CALV, ONLY: T2_CALV
    USE T2EXEC, ONLY: T2_EXEC
    USE T2CONV, ONLY: T2_CONV

    INTEGER(i0ikind),INTENT(OUT):: i_conv
    REAL(i0rkind),INTENT(OUT):: residual_conv
    INTEGER(i0ikind):: nconv
    INTEGER(i0ikind):: i0xidi,i0vidi
    CHARACTER(10)::c10nl
    REAL(4):: e0time_0,e0time_1
    
    c10nl='NL'
    
    DO i0xidi = 1, i0xmax
       DO i0vidi = 1, i0vmax
          d2xvec_befor( i0vidi,i0xidi)&
               = d2xvec(i0vidi,i0xidi)
       ENDDO
    ENDDO
    
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
       
       CALL T2_EXEC

       !C 
       !C CONVERGENCE CHECK
       !C
       
       CALL T2STEP_CONV(residual_conv)
       
       i_conv=nconv
       IF(residual_conv.LE.eps_conv) EXIT

       !C
       !C UPDATRE VARIABLES FOR NEXT ITERATION 
       !C
       
       DO i0xidi = 1, i0xmax
          DO i0vidi = 1, i0vmax
             d2xvec_befor(       i0vidi,i0xidi)&
                  = d2xvec_after(i0vidi,i0xidi)
          ENDDO
       ENDDO
       

    ENDDO
    
    print*,'NLLOOP=',nconv,'EXIT'
    
    DO i0xidi = 1, i0xmax
       DO i0vidi = 1, i0vmax
          d2xvec(i0vidi,i0xidi) = d2xvec_after(i0vidi,i0xidi)
       ENDDO
    ENDDO
    
    CALL T2_CONV
    
    RETURN
    
  END SUBROUTINE T2_STEP
  
  SUBROUTINE T2STEP_CONV(d0dif)
    
    USE T2COMM, ONLY:&
         i0solv,i0vmax,i0xmax,i0rmax,&
         i1mc1d,d2xvec_after,d2xvec_befor
    
    REAL(   i0rkind),INTENT(OUT)::d0dif
    INTEGER(i0ikind):: i0xidi, i0vidi, i0ridi
    REAL(   i0rkind):: d0aft, d0bfr, d0ave, d0dif_tmp,d0dif_max
    REAL(   i0rkind):: d1dif(1:i0vmax),d1ave(1:i0vmax)
    
    d0dif_max = 0.D0    
    d1dif(1:i0vmax) = 0.D0
    d1ave(1:i0vmax) = 0.D0
    
    
    !C
    !C FOR 2D VARIABLES
    !C
    
    DO i0xidi = 1, i0xmax  
       DO i0vidi = 4,i0vmax
          d0aft = d2xvec_after(i0vidi,i0xidi)
          d0bfr = d2xvec_befor(i0vidi,i0xidi)
          d1dif(i0vidi) = d1dif(i0vidi) + (d0aft-d0bfr)**2
          d1ave(i0vidi) = d1ave(i0vidi) + d0aft**2
       ENDDO
    ENDDO
     
    !C
    !C FOR 1D VARIABLES
    !C
    
    DO i0ridi = 1, i0rmax  
       i0xidi = i1mc1d(i0ridi)
       DO i0vidi = 1, 3
          d0aft = d2xvec_after(i0vidi,i0xidi)
          d0bfr = d2xvec_befor(i0vidi,i0xidi)
          d1dif(i0vidi) = d1dif(i0vidi) + (d0aft-d0bfr)**2
          d1ave(i0vidi) = d1ave(i0vidi) + d0aft**2
       ENDDO
    ENDDO
  
    !C
    !C CHECK CONVERGENCE
    !C
    SELECT CASE(i0solv)
       
       !C
       !C ELECTRON
       !C

    CASE(1)
       DO i0vidi = 6, 10
          
          d0dif  = d1dif(i0vidi)
          d0ave  = d1ave(i0vidi)
          
          IF(d0ave.LE.0.D0)THEN
             WRITE(6,'("*********************************************")')
             WRITE(6,'("       ERROR IN T2_STEP_CONVERGENCE          ")')
             WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
             WRITE(6,'("*********************************************")')
             PRINT*,i0vidi
             STOP
          ENDIF
          
          d0dif_tmp = d0dif/d0ave
          d0dif_tmp = SQRT(d0dif_tmp)
          WRITE(6,*),'VARIABLES=',i0vidi,'RESIDUAL=',d0dif_tmp
          d0dif_max = MAX(d0dif_max,d0dif_tmp)
          
       ENDDO
       
       !C
       !C ELECTRON AND IONS
       !C

    CASE (2)
       
       DO i0vidi = 6,i0vmax
          
          d0dif  = d1dif(i0vidi)
          d0ave  = d1ave(i0vidi)
          
          IF(d0ave.LE.0.D0)THEN
             WRITE(6,'("*********************************************")')
             WRITE(6,'("       ERROR IN T2_STEP_CONVERGENCE          ")')
             WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
             WRITE(6,'("*********************************************")')
             PRINT*,i0vidi
             STOP
          ENDIF
          
          d0dif_tmp = d0dif/d0ave
          d0dif_tmp = SQRT(d0dif_tmp)
          WRITE(6,*),'VARIABLES=',i0vidi,'RESIDUAL=',d0dif_tmp
          d0dif_max = MAX(d0dif_max,d0dif_tmp)

       ENDDO

       !C
       !C ELECTRON, IONS AND FIELD
       !C

    CASE (3)
       
       DO i0vidi = 1,i0vmax
          
          d0dif  = d1dif(i0vidi)
          d0ave  = d1ave(i0vidi)
          
          IF(d0ave.LE.0.D0)THEN
             WRITE(6,'("*********************************************")')
             WRITE(6,'("       ERROR IN T2_STEP_CONVERGENCE          ")')
             WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
             WRITE(6,'("*********************************************")')
             PRINT*,i0vidi
             STOP
          ENDIF
          
          d0dif_tmp = d0dif/d0ave
          d0dif_tmp = SQRT(d0dif_tmp)
          WRITE(6,*),'VARIABLES=',i0vidi,'RESIDUAL=',d0dif_tmp
          d0dif_max = MAX(d0dif_max,d0dif_tmp)
          
       ENDDO
       
    ENDSELECT
    
    d0dif = d0dif_max
    
    RETURN
    
  END SUBROUTINE T2STEP_CONV
  
END MODULE T2STEP
