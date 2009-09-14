      PROGRAM testmtx

      USE libmtx
      IMPLICIT NONE
      INTEGER:: imax,jwidth,isource
      INTEGER:: nrank,nprocs,i,istart,iend,its
      INTEGER:: nblock,nblock_per_proc,nblock_addition
      REAL(8):: v
      REAL(8),DIMENSION(:),POINTER:: x
      INTEGER,DIMENSION(2):: idata

      CALL mtx_initialize(nrank,nprocs)
      imax=11
      jwidth=3
      isource=6

    1 IF(nrank.eq.0) then
         WRITE(6,'(A,3I5)') '# INPUT: imax,isource=', &
                                      imax,isource
         READ(5,*,END=9000,ERR=1) imax,isource
         idata(1)=imax
         idata(2)=isource
      ENDIF
      CALL mtx_broadcast_integer(idata,2)
      imax=idata(1)
      isource=idata(2)

      IF(imax.EQ.0) GO TO 9000

      ALLOCATE(x(imax))

      CALL mtx_setup(imax,istart,iend,jwidth)

      DO i=istart,iend
         if(i.ne.1) CALL mtx_set_matrix(i,i-1,1.d0)
         CALL mtx_set_matrix(i,i, -2.d0)
         if(i.ne.imax) CALL mtx_set_matrix(i,i+1,1.d0)
      ENDDO

      DO i=istart,iend
         CALL mtx_set_vector(i,0.d0)
      ENDDO
      CALL mtx_set_vector(isource,-1.d0)
      CALL mtx_barrier

      CALL mtx_solve(its)
      if(nrank.eq.0) then
         write(6,*) 'Iteration Number=',its
      endif

      CALL mtx_gather_vector(x)

      if(nrank.eq.0) then
         DO i=1,imax
            WRITE(6,'(A,I5,1PE12.4)') 'i,x=',i,x(i)
         ENDDO
      endif

      CALL mtx_cleanup
      IF(nrank.EQ.0) DEALLOCATE(x)

      GO TO 1

 9000 CALL mtx_finalize
      STOP
      END PROGRAM testmtx
