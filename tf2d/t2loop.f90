!--------------------------------------------------------------------
! 
!         T2LOOP
!
!                    2014-06-05    H.Seto
!
!--------------------------------------------------------------------
MODULE T2LOOP
  
  USE T2CNST, ONLY:ikind,rkind
  IMPLICIT NONE
  
  PUBLIC T2_LOOP
  
  PRIVATE
  
CONTAINS
  
  SUBROUTINE T2_LOOP
    
    USE T2COMM,ONLY:&
         i0wstp,i1nlct,d1rsdl,LockEqs,&
         time_t2,dt,ntmax,ntstep,nt0dstep,nt0dmax,nt2dstep,nt2dmax,idfile
    USE T2PREC,ONLY: T2PREC_EXECUTE
    USE T2STEP,ONLY: T2_STEP
    USE T2SAVE,ONLY: T2_SAVE_0D,T2_SAVE_2D

    
    INTEGER(ikind):: nt
    REAL(4),SAVE  :: e0time1
    REAL(4)       :: e0time2!, e0time3, e0time4, e0time5
    INTEGER(ikind):: i_conv
    integer(ikind),save::i_init=0
    REAL(   rkind):: residual_conv

101 FORMAT('NT=',I6,2X,'TIME=',1P2E12.4,' [s]')
    
    CALL CPU_TIME(e0time1)
    
!    CALL T2PREC_EXECUTE
!    LockEqs(1) = .TRUE.
!    LockEqs(2) = .TRUE.
!    LockEqs(3) = .TRUE.
!    LockEqs(4) = .TRUE.
!    LockEqs(5) = .TRUE.
!    LockEqs(6) = .TRUE.
    ! time evolution loop

!    if(i_init.eq.0)THEN
!       LockEqs( 1)  = .FALSE.
!       LockEqs( 2)  = .FALSE.
!       LockEqs( 3)  = .FALSE.
!       LockEqs( 4)  = .TRUE.
!       LockEqs( 5)  = .TRUE.
!       LockEqs( 6)  = .TRUE.
!       
!       LockEqs( 7)  = .TRUE.
!       LockEqs( 8)  = .TRUE.
!       LockEqs( 9)  = .TRUE.
!       LockEqs(10)  = .TRUE.
!       LockEqs(11)  = .TRUE.
!       LockEqs(12)  = .TRUE.
!       LockEqs(13)  = .TRUE.
!       LockEqs(14)  = .TRUE.
!       LockEqs(15)  = .TRUE.
!       LockEqs(16)  = .TRUE.
       
!       LockEqs(17)  = .TRUE.
!       LockEqs(18)  = .TRUE.
!       LockEqs(19)  = .TRUE.
 !      LockEqs(20)  = .TRUE.
!       LockEqs(21)  = .TRUE.
!       LockEqs(22)  = .TRUE.
!       LockEqs(23)  = .TRUE.
!       LockEqs(24)  = .TRUE.
!       LockEqs(25)  = .TRUE.
!       LockEqs(26)  = .TRUE.
!       
!       CALL T2PREC_EXECUTE
!       
!       LockEqs( 1)  = .TRUE.
!       LockEqs( 2)  = .TRUE.
!       LockEqs( 3)  = .TRUE.
!       LockEqs( 4)  = .TRUE.
!       LockEqs( 5)  = .TRUE.
!       LockEqs( 6)  = .TRUE.
!       
!       LockEqs( 7)  = .TRUE.
!       LockEqs( 8)  = .FALSE.
!       LockEqs( 9)  = .TRUE.
!       LockEqs(10)  = .FALSE.
!       LockEqs(11)  = .FALSE.
!
!       LockEqs(12)  = .TRUE.
!       LockEqs(13)  = .TRUE.
!       LockEqs(14)  = .TRUE.
!       LockEqs(15)  = .TRUE.
!       LockEqs(16)  = .TRUE.
!
!       LockEqs(17)  = .TRUE.
!       LockEqs(18)  = .TRUE.
 !      LockEqs(19)  = .TRUE.
!       LockEqs(20)  = .TRUE.
!       LockEqs(21)  = .TRUE.!

!       LockEqs(22)  = .TRUE.
!       LockEqs(23)  = .TRUE.
!       LockEqs(24)  = .TRUE.
!       LockEqs(25)  = .TRUE.
!       LockEqs(26)  = .TRUE.
!       CALL T2PREC_EXECUTE
!       LockEqs(1) = .TRUE.
!       LockEqs(2) = .TRUE.
!       LockEqs(3) = .TRUE.
!       LockEqs( 4)  = .TRUE.
!       LockEqs( 5)  = .TRUE.
!       LockEqs( 6)  = .TRUE.

!       LockEqs( 7)  = .FALSE.
!       LockEqs( 8)  = .FALSE.
!       LockEqs( 9)  = .TRUE.
!       LockEqs(10)  = .FALSE.
!       LockEqs(11)  = .FALSE.
!       LockEqs(12)  = .TRUE.
!       LockEqs(13)  = .TRUE.
!       LockEqs(14)  = .TRUE.
!       LockEqs(15)  = .TRUE.
!       LockEqs(16)  = .TRUE.

!       LockEqs(17)  = .FALSE.
!       LockEqs(18)  = .FALSE.
!       LockEqs(19)  = .TRUE.
!       LockEqs(20)  = .FALSE.
!       LockEqs(21)  = .FALSE.
!       LockEqs(22)  = .TRUE.
!       LockEqs(23)  = .TRUE.
!       LockEqs(23)  = .TRUE.
!       LockEqs(25)  = .TRUE.
!       LockEqs(26)  = .TRUE.
!       i_init = 1
!    ENDIF
    
    DO nt=1,ntmax
       
       time_t2 = time_t2 + dt

       CALL T2_STEP(i_conv,residual_conv)
       i1nlct(nt) = i_conv
       d1rsdl(nt) = residual_conv
       

       CALL CPU_TIME(e0time2)
       
       IF(MOD(nt,ntstep).EQ.0) &
            WRITE(6,'(A,I6,2X,A,1PE12.4,2X,A,I6,2X,A,1PE12.4)') &
            'NT=',nt,'TIME=',time_t2,'Loop=',i_conv,'CPU=',e0time2-e0time1
       
!      ----- save data -----

       IF(MOD(nt,nt0dstep).EQ.0) CALL T2_SAVE_0D
       IF(MOD(nt,nt2dstep).EQ.0) CALL T2_SAVE_2D
       
    ENDDO

    RETURN
    
  END SUBROUTINE T2_LOOP
END MODULE T2LOOP
