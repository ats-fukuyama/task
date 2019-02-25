! MODULE tiexec

MODULE tiexec

  PRIVATE
  PUBLIC ti_exec
  
CONTAINS

  SUBROUTINE ti_exec(IERR)

    USE ticomm
    USE tirecord
    USE libmtx
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: v
    INTEGER:: NTL

    IERR=0

    ALLOCATE(v(imax))

    residual_loop_max=0.D0
    icount_loop_max=0
    icount_mat_max=0

    DO NTL=1,NTMAX
       NT=NT+1              ! total time step
       T=T+DT               ! total time 

       CALL ti_step(ierr)   ! one time-step caculcation

       residual_loop_max=MAX(residual_loop_max,residual_loop)
       icount_loop_max=MAX(icount_loop_max,icount_loop)
       icount_mat_max=MAX(icount_mat_max,icount_mat)

       IF(MOD(NT,NTSTEP ).EQ.0) CALL ti_snap         ! integrate and save
       IF(MOD(NT,NGTSTEP).EQ.0) CALL ti_record_ngt   ! save for time history
       IF(MOD(NT,NGRSTEP).EQ.0) CALL ti_record_ngr   ! save for radial profile
    END DO

    RETURN
  END SUBROUTINE ti_exec

  SUBROUTINE ti_step(ierr)

    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ICOUNT
    REAL(rkind):: residual
    INTEGER:: NR,NEQ,NSA,NV,i

    IERR=0

!   *** convert present value as a solution vector ***       

    DO NR=NR_START,NR_END
       DO NEQ=1,NEQMAX
          i=(NR-1)*NEQMAX+NEQ  ! position in solution vector
          IF(i.GE.istart.AND.i.LE.iend) THEN
             NSA=NSA_NEQ(NEQ)     ! species id
             NV=NV_NEQ(NEQ)       ! equation id
             SELECT CASE(NV)   ! setup solution vector at previous time
             CASE(1)
                SOL_BASE(i)=RNA(NSA,NR)
             CASE(2)
                SOL_BASE(i)=RTA(NSA,NR)
             CASE(3)
                SOL_BASE(i)=RUA(NSA,NR)
             END SELECT
          END IF
       END DO
    END DO
    CALL ti_distribute_sol(sol_base)

    DO i=istart,iend           ! setup sol_prev from sol_base
       SOL_PREV(i)=SOL_BASE(i)
    END DO
    CALL ti_distribute_sol(sol_prev)

!   *** start iteration using previous time step variables ***       

    ICOUNT=0

    DO WHILE(ICOUNT.LT.MAXLOOP)
       ICOUNT=ICOUNT+1

!   *** solve transport equation using previous iteration variables ***       

       CALL ti_solve           ! calculate sol_new from sol_prev

!   *** convergence check ***       

       CALL ti_convergence(residual)  ! check ABS(sol_new - sol_prev)
       IF(nrank.EQ.0.AND.idebug.EQ.1) THEN
          WRITE(6,'(A,I5,1PE12.4)') &
               'Convergence: ICOUNT,RESIDUAL=',ICOUNT,residual
       END IF

       IF(residual.LE.EPSLOOP) GO TO 8000

!   *** update sol_prev with new iteration vector sol_new ***       

       DO i=1,imax
          sol_prev(i)=sol_new(i)
       END DO

    END DO

    WRITE(6,*) 'XX ti_exec: No convergence'
    residual_loop=residual
    icount_loop=icount
    IERR=1
    RETURN

8000 CONTINUE

    residual_loop=residual
    icount_loop=icount

!   *** convert a new solution vector to present values *** 
!       since sol_new has all values, all RNA,RUA,RTA are calculated

    DO NR=1,NRMAX
       DO NEQ=1,NEQMAX
          i=(NR-1)*NEQMAX+NEQ
          NSA=NSA_NEQ(NEQ)
          NV=NV_NEQ(NEQ)
          SELECT CASE(NV)
          CASE(1)
             RNA(NSA,NR)=SOL_NEW(i)
          CASE(2)
             RTA(NSA,NR)=SOL_NEW(i)
          CASE(3)
             RUA(NSA,NR)=SOL_NEW(i)
          END SELECT
       END DO
    END DO
    IERR=0
    RETURN

  END SUBROUTINE ti_step

  SUBROUTINE ti_solve
    USE ticomm
    USE ticoef
    USE tisource
    USE ticalc
    USE libmtx
    USE libmpi
    IMPLICIT NONE
    INTEGER:: itype,i,j,NEQ,NEQ1,NR
    INTEGER:: istart_,iend_
    REAL(rkind):: tolerance

    DO NR=MAX(1,NR_START-1),MIN(NRMAX,NR_END+1)
       CALL ti_coef(NR)    ! calculate DD,VV,CC
       CALL ti_source(NR)  ! calculate SSIN,VSIN,PSIN,AJIN
    END DO

    CALL mtx_setup(imax,istart_,iend_,jwidth)
    IF(istart_.NE.istart.OR.iend_.NE.iend) THEN
       WRITE(6,*) 'XX ti_solve: mtx_setup inconsistency'
       STOP
    END IF

    DO NR=NR_START,NR_END

       CALL ti_calc(NR) ! assemble coeeficint matrix and RHS vector

       DO NEQ=1,NEQMAX
          i=(NR-1)*NEQMAX+NEQ
          IF(i.GE.istart.AND.i.LE.iend) THEN
             DO NEQ1=1,3*NEQMAX
                j=(NR-1)*NEQMAX+NEQ1-NEQMAX
                IF(j.ge.1.AND.j.LE.imax) THEN
                   IF(ABS(MAT_LOCAL(NEQ,NEQ1)).GT.0.D0) THEN
                      CALL mtx_set_matrix(i,j,MAT_LOCAL(NEQ,NEQ1))
!                      WRITE(6,'(A,2I5,1PE12.4)') &
!                           'mat_l:',i,j,MAT_LOCAL(NEQ,NEQ1)
                   END IF
                END IF
             END DO
             IF(ABS(VEC_LOCAL(NEQ)).GT.0.D0) THEN
                CALL mtx_set_source(i,VEC_LOCAL(NEQ))
!               WRITE(6,'(A,I5,5X,1PE12.4)') &
!                       'vec_l:',i,VEC_LOCAL(NEQ,NEQ1)
             END IF
             CALL mtx_set_vector(i,SOL_PREV(i))
          END IF
       END DO
    END DO

    itype=MATTYPE
    tolerance=EPSMAT
    CALL mtx_solve(itype,tolerance,icount_mat)
    
    CALL mtx_gather_vector(sol_new)

    CALL mtx_cleanup
  END SUBROUTINE ti_solve

  SUBROUTINE ti_convergence(RESIDUAL)
    USE ticomm
    USE libmpi
    IMPLICIT NONE
    REAL(rkind),INTENT(OUT):: RESIDUAL
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: bulk,diff,bulk_tot,diff_tot
    INTEGER,DIMENSION(:),ALLOCATABLE:: loc
    INTEGER:: NEQ,NR,i

    ALLOCATE(bulk(NEQMAX),diff(NEQMAX),bulk_tot(NEQMAX),diff_tot(NEQMAX))
    ALLOCATE(loc(NEQMAX))
    DO NEQ=1,NEQMAX
       bulk(NEQ)=0.D0
       diff(NEQ)=0.D0
    END DO

    DO NR=NR_START,NR_END
       DO NEQ=1,NEQMAX
          i=(NR-1)*NEQMAX+NEQ
          IF(i.GE.istart.AND.i.LE.iend) THEN
             bulk(NEQ)=bulk(NEQ)+SOL_NEW(i)**2
             diff(NEQ)=diff(NEQ)+(SOL_NEW(i)-SOL_PREV(i))**2
          END IF
       END DO
    END DO

    CALL mtx_reduce_real8(bulk,neqmax,3,bulk_tot,loc)
    CALL mtx_reduce_real8(diff,neqmax,3,diff_tot,loc)

    IF(nrank.eq.0) THEN
       residual=0.D0
       DO NEQ=1,NEQMAX
          residual=MAX(residual,SQRT(diff_tot(NEQ))/SQRT(1.D0+bulk_tot(NEQ)))
       END DO
    END IF
    CALL mtx_broadcast1_real8(residual)
    DEALLOCATE(bulk,diff,bulk_tot,diff_tot,loc)
    RETURN
  END SUBROUTINE ti_convergence

  SUBROUTINE ti_distribute_sol(sol)
    USE libmpi
    USE ticomm,ONLY: rkind,istart,iend,imax
    IMPLICIT NONE
    REAL(rkind),INTENT(INOUT):: sol(imax)
    REAL(rkind):: vsend(istart:iend)
    INTEGER:: i

    DO i=istart,iend
       vsend(i)=sol(i)
    END DO
    
    CALL mtx_allgather_real8(vsend,iend-istart+1,sol)
    RETURN
  END SUBROUTINE ti_distribute_sol
    
END MODULE tiexec
