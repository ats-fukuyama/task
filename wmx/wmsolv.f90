! wmsolv.f90

MODULE wmsolv

  PRIVATE
  PUBLIC wm_solv

CONTAINS
  
!     ****** SOLVE SIMULTANEOUS EQUATIONS ******

  SUBROUTINE wm_solv

    USE wmcomm
    IMPLICIT NONE
    INTEGER:: ierr

    mblock_size=3*MDSIZ*NDSIZ
    ! mlen=NRMAX*mblock_size+MWGMAX*NAMAX defined in wm_alloc in wmcomm.f90
    ! mbnd=4*mblock_size-1 defined in wm_alloc in wmcomm.f90
    
    CALL wm_solv_mtxp(CFVG,ierr)
    IF(IERR.NE.0) WRITE(6,*) 'XX wm_solv_mtxp error: ierr=',ierr

    RETURN
  END SUBROUTINE wm_solv

  ! --- solving complex matrix equation
  
  SUBROUTINE wm_solv_mtxp(svec,ierr)

    USE wmcomm
    USE libmtx
    USE wmsetm0
    USE wmsetm1
    USE wmsetm2
    USE libfio
    USE libmpi
    USE commpi
    IMPLICIT NONE
    COMPLEX(rkind),DIMENSION(MLEN),INTENT(OUT):: svec
    INTEGER,INTENT(OUT)::ierr 
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: A
    COMPLEX(rkind):: X
    INTEGER:: i,j,nr_previous,nfl,n
    INTEGER:: itype,its
    REAL(rkind):: tolerance

    nfl=30
    IF(idebuga(61).NE.0.AND.nrank.EQ.0) THEN
       CALL FWOPEN(nfl,knam_dump,1,0,'wm',ierr)
       IF(ierr.NE.0) STOP
    END IF

    ALLOCATE(A(MBND))        ! local coefficient matrix

    CALL mtxc_setup(MLEN,istart,iend,jwidth=MBND)

    nr_start=(istart-1)/mblock_size+1
    nr_end=(iend-1)/mblock_size+2
    IF(nrank.EQ.nsize-1) nr_end=NRMAX+1
    IF(ALLOCATED(nr_start_nrank)) DEALLOCATE(nr_start_nrank,nr_end_nrank)
    ALLOCATE(nr_start_nrank(0:nsize-1),nr_end_nrank(0:nsize-1))
    CALL mtx_allgather1_integer(nr_start,nr_start_nrank)
    CALL mtx_allgather1_integer(nr_end,nr_end_nrank)

    IF(idebuga(41).NE.0.AND.nrank.EQ.0) THEN
       DO n=0,nsize-1
          WRITE(6,'(A,4I8)') 'nrank,nr_start,nr_end,nr_len=', &
               n,nr_start_nrank(n),nr_end_nrank(n), &
               nr_end_nrank(n)-nr_start_nrank(n)+1
       END DO
    END IF

!   ***** CALCULATE MATRIX COEFFICIENTS *****

    nr_previous=0

    DO i=istart,iend
       X=(0.D0,0.D0)
       A(1:MBND)=(0.D0,0.D0)
       SELECT CASE(MDLWMX)
       CASE(0)
          CALL wm_setm0(A,X,i,MBND,nr_previous)
       CASE(1)
          CALL wm_setm1(A,X,i,MBND,nr_previous)
       CASE(2)
          CALL wm_setm2(A,X,i,MBND,nr_previous)
       END SELECT
       IF(idebuga(61).NE.0.AND.nrank.EQ.0) &
            WRITE(nfl,'(A,2I6,ES12.4)') &
            'wmsolv:',nr_previous,i,XRHO(nr_previous)
       IF(idebuga(61).NE.0.AND.nrank.EQ.1) &
            WRITE(nfl+1,'(A,2I6,ES12.4)') &
            'wmsolv:',nr_previous,i,XRHO(nr_previous)
       IF(idebuga(61).NE.0.AND.nrank.EQ.2) &
            WRITE(nfl+2,'(A,2I6,ES12.4)') &
            'wmsolv:',nr_previous,i,XRHO(nr_previous)
       IF(idebuga(61+3).NE.0.AND.nrank.EQ.3) &
            WRITE(nfl,'(A,2I6,ES12.4)') &
            'wmsolv:',nr_previous,i,XRHO(nr_previous)
       DO j=MAX(i-MCENT+1,1),MIN(MLEN,i+MCENT-1)
          IF(ABS(A(j-i+MCENT)).GT.0.D0) THEN
             CALL mtxc_set_matrix(i,j,A(j-i+MCENT))
             IF(idebuga(61).NE.0.AND.nrank.EQ.0) &
                  WRITE(nfl,'(A,2I6,2ES12.4)') 'A:',i,j,A(j-i+MCENT)
             IF(idebuga(61).NE.0.AND.nrank.EQ.1) &
                  WRITE(nfl+1,'(A,2I6,2ES12.4)') 'A:',i,j,A(j-i+MCENT)
             IF(idebuga(61).NE.0.AND.nrank.EQ.2) &
                  WRITE(nfl+2,'(A,2I6,2ES12.4)') 'A:',i,j,A(j-i+MCENT)
             IF(idebuga(61).NE.0.AND.nrank.EQ.3) &
                  WRITE(nfl+3,'(A,2I6,2ES12.4)') 'A:',i,j,A(j-i+MCENT)
          END IF
       END DO
       IF(ABS(X).GT.0.D0) THEN
          CALL mtxc_set_source(i,X)
          IF(idebuga(61).NE.0.AND.nrank.EQ.0) &
               WRITE(nfl,'(A,2I6,2ES12.4)') 'X:',i,0,X
          IF(idebuga(61).NE.0.AND.nrank.EQ.1) &
               WRITE(nfl+1,'(A,2I6,2ES12.4)') 'X:',i,0,X
          IF(idebuga(61).NE.0.AND.nrank.EQ.2) &
               WRITE(nfl+2,'(A,2I6,2ES12.4)') 'X:',i,0,X
          IF(idebuga(61).NE.0.AND.nrank.EQ.3) &
               WRITE(nfl+3,'(A,2I6,2ES12.4)') 'X:',i,0,X
       END IF
    END DO

    itype=1  ! infolevel for MUMPS
    tolerance=1.D-12
    CALL mtxc_solve(itype,tolerance,its)
    IF(nrank.EQ.0) WRITE(6,'(A,I8)') '## wm_solv: iteration=',its
      
    CALL mtxc_gather_vector(svec)

    IF(idebuga(61).NE.0.AND.nrank.EQ.0) THEN
       DO i=1,MLEN,3
          WRITE(61,'(I6,6ES12.4)') &
               i,svec(i),svec(i+1),svec(i+2)
       END DO
    END IF
    CALL mtxc_cleanup

    DEALLOCATE(A)
    IERR=0

    RETURN
  END SUBROUTINE wm_solv_mtxp
END MODULE wmsolv
