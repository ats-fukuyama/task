!C--------------------------------------------------------------------
!C 
!C
!C
!C
!C
!C--------------------------------------------------------------------
MODULE T2LOOP
  
  USE T2CNST, ONLY:i0ikind,i0rkind
  IMPLICIT NONE
  
  PUBLIC T2_LOOP
  PRIVATE
  
CONTAINS
  
  SUBROUTINE T2_LOOP
    
    USE T2COMM, ONLY:&
         i0tmax,i0wstp,i0tstp,i0nlct,&
         d0time,d0tstp,d0tmax,d0eps,d1guv,c10rname
    USE T2PROF, ONLY:&
         T2_PROF
    USE T2STEP, ONLY:&
         T2_STEP
    USE T2WRIT, ONLY:&
         T2WRIT_MAIN,T2_WRIT_GPT,T2_WRIT_GP1
    
    INTEGER(i0ikind)     :: i0tflg,i0tsws
    REAL(4),SAVE         :: e0time1
    REAL(4)              :: e0time2, e0time3, e0time4, e0time5
    REAL(   i0rkind),SAVE:: d0tstp_save

101 FORMAT('TIMESTEP=',I10,1X,'LAP TIME=',F10.4,'[s/step]',1X,&
         'ELAPSED TIME=',F10.4,'[s]')
    
    i0tsws      = 0
    i0tstp      = 0
    d0time      = 0.D0
    d0tstp_save = d0tstp

    DO !C TIME EVOLUTION LOOP
       
       !C
       !C CALCULATE GEOMETRICAL COEFFICIENTS
       !C
       
       IF(i0tstp.EQ.0) CALL T2_PROF
       
       IF(i0tstp.EQ.0) CALL T2WRIT_MAIN(d1guv,i0tstp,c10rname)
       IF(i0tstp.EQ.0) CALL T2_WRIT_GPT(20,i0tstp,d1guv)
       IF(i0tstp.EQ.0) CALL T2_WRIT_GP1(22,i0tstp,d1guv)
       IF(i0tmax.EQ.0) EXIT
       i0tstp = i0tstp + 1
       i0tflg = 0
       IF(i0tstp.EQ.1) CALL CPU_TIME(e0time1)
       CALL CPU_TIME(e0time2)
       
       IF((d0time + d0tstp).GE.d0tmax)THEN
          d0tstp = d0tmax - d0time
          i0tflg = 1
       ELSEIF(i0tstp.EQ.(i0tmax-1))THEN
          i0tmax = i0tmax+1
       ENDIF

       d0time = d0time + d0tstp
       
       !C
       !C NON-LINEAR ITERATION BY PICARD METHOD
       !C
       
       WRITE(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(6,*)'     NON-LINEAR ITERATION BY PICARD METHOD        '
       WRITE(6,*)'TIMESTEPNUMBER=',i0tstp,'TIMESTEPWIDTH=',d0tstp
       WRITE(6,*)'TOLERANCE=',d0eps
       
       CALL T2_STEP
       
       CALL CPU_TIME(e0time3)
       
       e0time4 = e0time3 - e0time2
       e0time5 = e0time3 - e0time1
       
       WRITE(6,101),i0tstp,e0time4,e0time5
       
       !C
       !C WRITE CALCULATION RESULT IN VTK FORMAT
       !C
       
       IF(MOD(i0tstp,i0wstp).EQ.0) CALL T2WRIT_MAIN(d1guv,i0tstp,c10rname)
       IF(MOD(i0tstp,i0wstp).EQ.0) CALL T2_WRIT_GPT(20,i0tstp,d1guv)       
       WRITE(6,*),'IN T2LOOP','i0tstp=',i0tstp,'i0tmax=',i0tmax

       IF((i0tstp.EQ.i0tmax).OR.(i0tflg.EQ.1))THEN
          d0tstp= d0tstp_save
          EXIT
       END IF
    ENDDO

    RETURN

  END SUBROUTINE T2_LOOP

END MODULE T2LOOP
