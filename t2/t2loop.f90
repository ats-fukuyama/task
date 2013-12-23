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
         i0wstp,i1nlct,d1rsdl,&
         d1guv,c10rname, &
         time_t2,dt,ntmax,ntstep,nt0dstep,nt0dmax,nt2dstep,nt2dmax,idfile
    USE T2STEP, ONLY:&
         T2_STEP
    USE T2WRIT, ONLY:&
         T2_WRIT_MAIN,T2_WRIT_GPT,T2_WRIT_GP1
    USE T2SAVE, ONLY:&
         T2_SAVE_0D,T2_SAVE_2D
    USE LIBMTX, ONLY: &
         MTX_INITIALIZE,MTX_FINALIZE
    
    INTEGER(i0ikind)     :: i0tflg,i0tsws,nt
    REAL(4),SAVE         :: e0time1
    REAL(4)              :: e0time2, e0time3, e0time4, e0time5
    REAL(   i0rkind),SAVE:: d0tstp_save
    INTEGER(i0ikind)     :: i_conv
    REAL(i0rkind):: residual_conv

101 FORMAT('NT=',I6,2X,'TIME=',1P2E12.4,' [s]')
    
    CALL CPU_TIME(e0time1)

    CALL MTX_INITIALIZE

!   ----- TIME EVOLUTION LOOP -----

    DO nt=1,ntmax
       
       time_t2 = time_t2 + dt
       
       CALL T2_STEP(i_conv,residual_conv)
       i1nlct(nt) = i_conv
       d1rsdl(nt) = residual_conv
       
       CALL CPU_TIME(e0time2)

       IF(MOD(nt,ntstep).EQ.0) &
       WRITE(6,'(A,I6,2X,A,1P2E12.4,2X,A,I6,2X,A,1P2E12.4)') &
            'NT=',nt,'TIME=',time_t2,'Loop=',i_conv,'CPU=',e0time2-e0time1
       
       !C
       !C WRITE CALCULATION RESULT IN VTK FORMAT
       !C
       
       IF(idfile.GE.1) CALL T2_WRIT_MAIN(d1guv,nt,c10rname)
       IF(idfile.GE.2) CALL T2_WRIT_GPT(20,nt,d1guv)
       IF(idfile.GE.3) CALL T2_WRIT_GP1(22,nt,d1guv)

!      ----- save data -----

       IF(MOD(nt,nt0dstep).EQ.0) CALL T2_SAVE_0D
       IF(MOD(nt,nt2dstep).EQ.0) CALL T2_SAVE_2D

    ENDDO

    CALL MTX_FINALIZE

    RETURN

  END SUBROUTINE T2_LOOP

END MODULE T2LOOP
