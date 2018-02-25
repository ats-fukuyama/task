! MODULE tiloop

MODULE tiloop

  PRIVATE
  PUBLIC ti_loop

CONTAINS

  SUBROUTINE ti_loop(IERR)

    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NTL

    IERR=0

    DO NTL=1,NTMAX
       NT=NT+1
       T=T+DT

       CALL ti_exec(ierr)
       IF(ierr.ne.0) exit

       IF(MOD(NT,NTSTEP).EQ.0) CALL ti_snap
       IF(MOD(NT,NGTSTEP).EQ.0) CALL ti_save_ngt
       IF(MOD(NT,NGRSTEP).EQ.0) CALL ti_save_ngr
    END DO
    RETURN
  END SUBROUTINE ti_loop

  SUBROUTINE ti_exec(ierr)

    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ICOUNT
    REAL(rkind):: RESIDUAL
    INTEGER:: NR,NEQ,NSA,NV,i

    IERR=0

!   *** save present value as previous value ***       

    DO NR=1,NRMAX
       DO NEQ=1,NEQMAX
          NSA=NSA_NEQ(NEQ)
          NV=NV_NEQ(NEQ)
          i=(NR-1)*NEQMAX+NEQ
          SELECT CASE(NV)
          CASE(1)
             SOL_PREV(i)=RNA(NSA,NR)
          CASE(2)
             SOL_PREV(i)=RUA(NSA,NR)
          CASE(3)
             SOL_PREV(i)=RTA(NSA,NR)
          END SELECT
       END DO
    END DO

    ICOUNT=0
    DO
       ICOUNT=ICOUNT+1

       CALL ti_solve

!   *** update present value ***       

       DO NR=1,NRMAX
          DO NEQ=1,NEQMAX
             NSA=NSA_NEQ(NEQ)
             NV=NV_NEQ(NEQ)
             i=(NR-1)*NEQMAX+NEQ
             SELECT CASE(NV)
             CASE(1)
                RNA(NSA,NR)=SOL_NEW(i)
             CASE(2)
                RUA(NSA,NR)=SOL_NEW(i)
             CASE(3)
                RTA(NSA,NR)=SOL_NEW(i)
             END SELECT
          END DO
    END DO

       CALL ti_convergence(RESIDUAL)
       WRITE(6,'(A,I5,1PE12.4)') &
            'Convergence: ICOUNT,RESIDUAL=',ICOUNT,RESIDUAL
       IF(RESIDUAL.LE.EPSLTI) GOTO 8000
       IF(ICOUNT.GE.LMAXTI) GOTO 9000
    END DO

8000 RETURN

9000 WRITE(6,*) 'XX ti_exec: No convergence'
    IERR=1
    RETURN
  END SUBROUTINE ti_exec

  SUBROUTINE ti_solve
    USE ticomm
    USE ticoef
    USE tisource
    USE ticalc
    USE libmtx
    IMPLICIT NONE
    INTEGER:: imax,jwidth,istart,iend,NRSTART,NREND,isize
    INTEGER:: itype,iterations,i,j,IL,JL,NR
    REAL(rkind):: tolerance
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: v

    imax=NRMAX*NEQMAX
    ALLOCATE(v(imax))
    jwidth=4*NEQMAX-1
    CALL mtx_setup(imax,istart,iend,jwidth)
    
    NRSTART=(istart+NEQMAX-1)/NEQMAX
    NREND=(iend+NEQMAX-1)/NEQMAX

    DO NR=NRSTART,NREND
       CALL ti_coef(NR) 
       CALL ti_source(NR) 
       CALL ti_calc(NR)
       DO IL=1,NEQMAX
          i=(NR-1)*NEQMAX+IL
          IF(i.GE.istart.AND.i.LE.iend) THEN
             DO JL=1,3*NEQMAX
                j=(NR-1)*NEQMAX+JL-NEQMAX
                IF(j.ge.1.AND.j.LE.NRMAX*NEQMAX) THEN
                   IF(ABS(MAT_LOCAL(IL,JL)).GT.0.D0) THEN
                      CALL mtx_set_matrix(i,j,MAT_LOCAL(IL,JL))
                   END IF
                END IF
             END DO
             CALL mtx_set_source(i,VEC_LOCAL(IL))
             CALL mtx_set_vector(i,SOL_PREV(i))
          END IF
       END DO
    END DO

    itype=0
    tolerance=1.D-10
    CALL mtx_solve(itype,tolerance,iterations)
    WRITE(6,'(A,I5)') '## Iteration =',iterations

    CALL mtx_gather_vector(v)
    DO i=1,imax
       SOL_NEW(i)=v(i)
    END DO
    DEALLOCATE(v)
  END SUBROUTINE ti_solve

  SUBROUTINE ti_snap
    USE ticomm
    IMPLICIT NONE
    WRITE(6,'(A,I5,1PE12.4)') '#NT,T=',NT,T
    RETURN
  END SUBROUTINE ti_snap

  SUBROUTINE ti_save_ngt
    USE ticomm
    IMPLICIT NONE
    RETURN
  END SUBROUTINE ti_save_ngt

  SUBROUTINE ti_save_ngr
    USE ticomm
    IMPLICIT NONE
    RETURN
  END SUBROUTINE ti_save_ngr

  SUBROUTINE ti_convergence(RESIDUAL)
    USE ticomm
    IMPLICIT NONE
    REAL(rkind),INTENT(OUT):: RESIDUAL
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: bulk,diff
    INTEGER:: NEQ,NR,i

    ALLOCATE(bulk(NEQMAX),diff(NEQMAX))
    DO NEQ=1,NEQMAX
       bulk(NEQ)=0.D0
       diff(NEQ)=0.D0
    END DO
    DO NR=1,NRMAX
       DO NEQ=1,NEQMAX
          i=(NR-1)*NEQMAX+NEQ
          bulk(NEQ)=bulk(NEQ)+SOL_NEW(i)**2
          diff(NEQ)=diff(NEQ)+(SOL_NEW(i)-SOL_PREV(i))**2
       END DO
    END DO
    RESIDUAL=0.D0
    DO NEQ=1,NEQMAX
       RESIDUAL=MAX(RESIDUAL,SQRT(diff(NEQ))/SQRT(1.D0+bulk(NEQ)))
    END DO
    DEALLOCATE(bulk,diff)
    RETURN
  END SUBROUTINE ti_convergence
END MODULE tiloop
