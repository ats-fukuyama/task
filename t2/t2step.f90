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

  PUBLIC T2STEP_MAIN
  PRIVATE
  
CONTAINS

  SUBROUTINE T2STEP_MAIN
    
    USE T2COMM, ONLY:&
         i0pmax, i0cmax, i0vmax, i0xmax, i0bmax, i0tstp,i0nlct,&
         d0eps,&
         i1nlct, d1rsdl,&
         d1gsm,  d1grv,  d1guv,  d1guv_after, d1guv_befor
    
!    USE T2CALV, ONLY: T2CALV_MAIN
    USE T2CALV, ONLY: T2_CALV
    USE T2EXEC, ONLY: T2_EXEC
    USE T2WRIT
    INTEGER(i0ikind):: i0pflg, i1
    REAL(   i0rkind):: d0aft,  d0bfr, d0dif, d0ave, d0dif_max
    CHARACTER(10)::c10nl
    i0nlct = 0
    c10nl='NL'

    DO i1=1,i0xmax
       d1guv_befor(i1)= d1guv(i1)
    ENDDO
    
    DO
       
       DO i1=1,i0cmax
          d1gsm(i1)=0.d0
       ENDDO
       
       DO i1=1,i0bmax
          d1grv(i1)=0.d0
       ENDDO
       
       i0nlct = i0nlct+1
       
       !C
       !C CALCULATE PLASMA COEFFICIENTS
       !C
       
       CALL T2_CALV
       
       !C
       !C ADVECTION DIFFUSION EQUATION SOLVER (SUPG)
       !C
       
       CALL T2_EXEC

       !C 
       !C CONVERGENCE CHECK
       !C
       print*,i0nlct

       CALL T2STEP_CONVERGENCE(i0pflg,d0dif)

       IF(i0pflg.EQ.1)THEN
          WRITE(6,'("         PICARD ITERATION LOOP CONVERGED          ")')
          WRITE(6,*)'RESIDUAL=',d0dif,'TOLERANCE=',d0eps
          WRITE(6,'("**************************************************")')
          EXIT
       ENDIF

       !C>>>>>> DEBUG
       print*,'DEBUG'
       
       CALL T2_WRIT_GPT(21,i0nlct,d1guv_after)
       CALL T2WRIT_MAIN(d1guv_after,i0nlct,c10nl)

       !<<<<<<
       
       !C
       !C UPDATRE VARIABLES FOR NEXT ITERATION 
       !C
    
       DO i1=1,i0xmax
          d1guv_befor(i1)= d1guv_after(i1)
          d1guv_after(i1)=0.d0
       ENDDO
              
       IF(i0nlct.EQ.i0pmax)THEN
          WRITE(6,'("PICARD ITERATION LOOP CANNOT CONVERGED")')
          WRITE(6,*)'ITERATION NUMBER=',i0nlct, 'RESIDUAL=',d0dif
          STOP
       ENDIF
    ENDDO
    
    i1nlct(i0tstp) = i0nlct
    d1rsdl(i0tstp) = d0dif
    
    DO i1=1,i0xmax
       d1guv(i1)=d1guv_after(i1)
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2STEP_MAIN
  
  SUBROUTINE T2STEP_CONVERGENCE(i0flg,d0dif)
    
    USE T2COMM, ONLY:&
         i0spcs,i0dbg,i0vmax,i0nmax2, d0eps, d1guv_after, d1guv_befor
    
    INTEGER(i0ikind),INTENT(OUT)::i0flg
    REAL(   i0rkind),INTENT(OUT)::d0dif
    INTEGER(i0ikind):: i0xnt, i1, i2
    REAL(   i0rkind):: d0aft, d0bfr, d0ave, d0dif_tmp,d0dif_max
    
    i0flg = 0   
    d0dif_max = 0.D0    

    DO i1=1,i0vmax
!        IF(i1.ne.i0dbg) cycle
        !IF((i1.GT.5).OR.(i1.EQ.3)) cycle
!       IF(i1.LE.5) cycle
       IF((i1.GE.6).AND.(MOD(i1-6,8).GE.i0dbg)) CYCLE
!       IF((i1.GE.6).AND.(MOD(i1-6,8).EQ.1)) CYCLE
       d0dif     = 0.D0
       d0ave     = 0.D0
       
       DO i2=0,i0nmax2-1
          i0xnt = i0vmax*i2 + i1
          d0aft = d1guv_after(i0xnt)
          d0bfr = d1guv_befor(i0xnt)
          d0dif = d0dif + (d0aft-d0bfr)**2
          d0ave = d0ave + d0aft**2
       ENDDO
       
       d0dif  = SQRT(d0dif)
       d0ave  = SQRT(d0ave)
       
       IF(d0ave.LE.0.D0)THEN
          WRITE(6,'("*********************************************")')
          WRITE(6,'("       ERROR IN T2STEP_CONVERGENCE           ")')
          WRITE(6,'("       INDETERMINATE PROBLEM                 ")')
          WRITE(6,'("*********************************************")')
          PRINT*,i1
          STOP
       ENDIF
       
       d0dif_tmp = d0dif/d0ave
       d0dif_max = MAX(d0dif_max,d0dif_tmp)
       WRITE(6,*),'VARIABLES=',i1,'RESIDUAL=',d0dif_tmp
    ENDDO

    d0dif = d0dif_max
    
    IF(d0dif_max.LT.d0eps)THEN
       i0flg = 1
    END IF
    
    RETURN
    
  END SUBROUTINE T2STEP_CONVERGENCE
  
END MODULE T2STEP
