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

    nblock_size=3*MDSIZ*NDSIZ
    ! mlen=NRMAX*nblock_size+MWGMAX*NAMAX defined in wm_alloc in wmcomm.f90
    ! mbnd=4*nblock_size-1 defined in wm_alloc in wmcomm.f90
    
    CALL wm_solv_mtxp(CFVG,ierr)
    IF(IERR.NE.0) WRITE(6,*) 'XX wm_solv_mtxp error: ierr=',ierr

    RETURN
  END SUBROUTINE wm_solv

  ! --- solving complex matrix equation
  
  SUBROUTINE wm_solv_mtxp(svec,ierr)

    USE wmcomm
    USE libmtx
    USE wmsetm
    USE commpi
    IMPLICIT NONE
    COMPLEX(rkind),DIMENSION(MLEN),INTENT(OUT):: svec
    INTEGER,INTENT(OUT)::ierr 
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: A
    COMPLEX(rkind):: X
    INTEGER:: i,j,nr_previous
    INTEGER:: itype,its
    REAL(rkind):: tolerance

    ALLOCATE(A(MBND))        ! local coefficient matrix

    CALL mtxc_setup(MLEN,istart,iend,jwidth=MBND)

    nr_start=(istart-1)/nblock_size+1
    nr_end=(iend-1)/nblock_size+1
    IF(nrank.EQ.nsize-1) nr_end=NRMAX
    WRITE(6,'(A,4I8)') 'nrank,nr_start,nr_end,nblock_size=', &
                        nrank,nr_start,nr_end,nblock_size

!   ***** CALCULATE MATRIX COEFFICIENTS *****

    nr_previous=0

    DO i=istart,iend
       X=(0.D0,0.D0)
       A(1:MBND)=(0.D0,0.D0)
       CALL wm_setm(A,X,i,MBND,nr_previous)
       WRITE(21,'(A,2I6,ES12.4)') 'wmsolv:',nr_previous,i,XRHO(nr_previous)
       WRITE(21,'(6ES12.4)') (A(j),j=1,MBND)
       WRITE(21,'(2ES12.4)') X
       DO j=MAX(i-MCENT+1,1),MIN(MLEN,i+MCENT-1)
          IF(ABS(A(j-i+MCENT)).GT.0.D0) THEN
             CALL mtxc_set_matrix(i,j,A(j-i+MCENT))
          END IF
       END DO
       CALL mtxc_set_source(i,X)
    END DO

    itype=0
    tolerance=1.D-12
    CALL mtxc_solve(itype,tolerance,its)
    WRITE(6,'(A,I8)') '## wm_solv: iteration=',its
      
    CALL mtxc_gather_vector(svec)

    CALL mtxc_cleanup

    DEALLOCATE(A)
    IERR=0

    RETURN
  END SUBROUTINE wm_solv_mtxp
END MODULE wmsolv
