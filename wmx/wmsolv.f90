! wmsolv.f90

MODULE wmsolv

  PRIVATE
  PUBLIC wm_solv

CONTAINS
  
!     ****** SOLVE SIMULTANEOUS EQUATIONS ******

  SUBROUTINE wm_solv

    USE wmcomm
    IMPLICIT NONE
    INTEGER:: nblock_size
    INTEGER:: ierr

    nblock_size=3*MDSIZ*NDSIZ
    ! mlen=NRMAX*nblock_size+MWGMAX*NAMAX defined in wm_alloc in wmcomm.f90
    ! mbnd=4*nblock_size-1 defined in wm_alloc in wmcomm.f90
    
    CALL wm_solv_mtxp(CFVG,MLEN,MBND,nblock_size,ierr)
    IF(IERR.NE.0) WRITE(6,*) 'XX wm_solv_mtxp error: ierr=',ierr

    RETURN
  END SUBROUTINE wm_solv

  ! --- solving complex matrix equation
  
  SUBROUTINE wm_solve_mtxp(svec,MLEN,MBND,nblock_size,ierr)

    USE wmcomm,ONLY: rkind,NRMAX,MODEWG,MWGMAX,NAMAX
    USE libmtx
    USE wmsetm
    USE commpi
    IMPLICIT NONE
    COMPLEX(rkind),DIMENSION(MLEN),INTENT(OUT):: svec
    INTEGER,INTENT(IN):: MLEN,MBND,nblock_size
    INTEGER,INTENT(OUT)::ierr 
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: A
    COMPLEX(rkind):: X
    INTEGER:: nblock_start,nblock_end
    INTEGER:: istart,iend,i,j,nr_previous
    INTEGER:: itype,its
    REAL(rkind):: tolerance

    ALLOCATE(A(MBND))        ! local coefficient matrix

    CALL mtxc_setup(MLEN,istart,iend,jwidth=MBND)

    nblock_start=(istart-1)/nblock_size+1
    nblock_end=(iend-1)/nblock_size+1
    IF(nrank.EQ.nsize-1) nblock_end=NRMAX
    WRITE(6,'(A,3I8)') 'nrank,nblock_start,nblock_end=', &
                        nrank,nblock_start,nblock_end

!   ***** CALCULATE MATRIX COEFFICIENTS *****

    nr_previous=0

    DO i=istart,iend
       X=(0.D0,0.D0)
       A(1:MBND)=(0.D0,0.D0)
       CALL wm_setm(A,X,i,MBND,nr_previous)
       DO j=MAX(i-(MBND+1)/2+1,1),MIN(MLEN,i+(MBND+1)/2-1)
          IF(ABS(A(j-i+(MBND+1)/2)).GT.0.D0) THEN
             CALL mtxc_set_matrix(i,j,A(j-i+(MBND+1)/2))
          END IF
       END DO
       CALL mtxc_set_source(i,X)
    END DO

    itype=0
    tolerance=1.D-12
    CALL mtxc_solve(itype,tolerance,its)
    WRITE(6,'(A,I8)') '## wmsolv: iteration=',its
      
    CALL mtxc_gather_vector(svec)

    CALL mtxc_cleanup

    DEALLOCATE(A)
    IERR=0

    RETURN
  END SUBROUTINE wm_solve_mtxp
END MODULE wmsolv
