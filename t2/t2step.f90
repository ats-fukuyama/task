!C--------------------------------------------------------------------
!C
!C T2STEP
!C
!C
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
         i0pmax, i0cmax, i0vmax, i0xmax, i0bmax, i0nlct,&
         d0eps,&
         i1nlct, d1rsdl,&
         d1gsm,  d1grv,  d1guv,  d1guv_after, d1guv_befor, &
         nconvmax,eps_conv
    
    USE T2CALV, ONLY: T2_CALV
    USE T2EXEC, ONLY: T2_EXEC
    USE T2WRIT
   
    INTEGER(i0ikind),INTENT(OUT):: i_conv
    REAL(i0rkind),INTENT(OUT):: residual_conv
    INTEGER(i0ikind):: nconv
    INTEGER(i0ikind):: i0pflg, i1
    REAL(   i0rkind):: d0aft,  d0bfr, d0dif, d0ave, d0dif_max
    CHARACTER(10)::c10nl
    REAL(4):: e0time_0,e0time_1
    
    c10nl='NL'

    DO i1=1,i0xmax
       d1guv_befor(i1)= d1guv(i1)
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
       
       CALL T2_STEP_CONV(residual_conv)
       
       i_conv=nconv
       IF(residual_conv.LE.eps_conv) EXIT

       !C
       !C UPDATRE VARIABLES FOR NEXT ITERATION 
       !C
       
       DO i1=1,i0xmax
          d1guv_befor(i1) = d1guv_after(i1)
       ENDDO
       
    ENDDO
    
    DO i1=1,i0xmax
       d1guv(i1) = d1guv_after(i1)
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2_STEP
  
  SUBROUTINE T2_STEP_CONV(d0dif)
    
    USE T2COMM, ONLY:&
         i0spcs,i0dbg,i0vmax,i0nmax2,i0nmax4,&
         i1mfc4,d0eps,d1guv_after,d1guv_befor
    
    REAL(   i0rkind),INTENT(OUT)::d0dif
    INTEGER(i0ikind):: i0xnt, i1, i2
    REAL(   i0rkind):: d0aft, d0bfr, d0ave, d0dif_tmp,d0dif_max
    
    d0dif_max = 0.D0    

    DO i1=1,i0vmax
       !IF((i1.GE.6).AND.(MOD(i1-6,8).GE.i0dbg)) CYCLE
       IF(i1.NE.3) CYCLE
!       IF(    i1.GT.3)THEN

          d0dif     = 0.D0
          d0ave     = 0.D0
          
          DO i2=0,i0nmax2-1
             
             i0xnt = i0vmax*i2 + i1
             d0aft = d1guv_after(i0xnt)
             d0bfr = d1guv_befor(i0xnt)
             d0dif = d0dif + (d0aft-d0bfr)**2
             d0ave = d0ave + d0aft**2
             
          ENDDO
          
 !      ELSEIF(i1.LE.3)THEN
 !         
 !         d0dif     = 0.D0
 !         d0ave     = 0.D0
 !         
 !         DO i2=1,i0nmax4
 !            i0xnt = i1mfc4(i2) - 1
 !            i0xnt = i0vmax*i0xnt + i1
 !            d0aft = d1guv_after(i0xnt)
 !            d0bfr = d1guv_befor(i0xnt)
 !            d0dif = d0dif + (d0aft-d0bfr)**2
 !            d0ave = d0ave + d0aft**2
 !         ENDDO
 !         
 !      ENDIF
      
       d0dif  = SQRT(d0dif)
       d0ave  = SQRT(d0ave)
       
       IF(d0ave.LE.0.D0)THEN
          WRITE(6,'("*********************************************")')
          WRITE(6,'("       ERROR IN T2_STEP_CONVERGENCE          ")')
          WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
          WRITE(6,'("*********************************************")')
          PRINT*,i1
          RETURN
       ENDIF
       
       d0dif_tmp = d0dif/d0ave
       d0dif_max = MAX(d0dif_max,d0dif_tmp)
       WRITE(6,*),'VARIABLES=',i1,'RESIDUAL=',d0dif_tmp
    ENDDO
    
    d0dif = d0dif_max
    
    RETURN
    
  END SUBROUTINE T2_STEP_CONV
  
END MODULE T2STEP
