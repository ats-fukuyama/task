      PROGRAM testksp

      use libmtx
      implicit none
      integer:: i,n,m,istart,iend,irank,isize,its,id
      real(8):: v

      CALL mtx_initialize(irank,isize)

      n=11
      m=3
      id=0
    1 WRITE(6,'(A,2I5)') 'INPUT: n,id=',n,id
      READ(5,*,END=9000,ERR=1) n,id
      IF(n.EQ.0) GO TO 9000

      CALL mtx_setup(n,m,istart,iend)
      write(6,*) 'irank,isize=',irank,isize
      write(6,*) 'istart,iend=',istart,iend

      DO i=istart,iend
         if(i.ne.1) CALL mtx_set_matrix(i,i-1,1.d0)
         CALL mtx_set_matrix(i,i, -2.d0)
         if(i.ne.n) CALL mtx_set_matrix(i,i+1,1.d0)
      ENDDO

      DO i=istart,iend
         CALL mtx_set_vector(i,0.d0)
      ENDDO

      SELECT CASE(ID)
      CASE(0) 
         CALL mtx_set_vector((n+1)/2,1.d0)
      CASE(1) 
         CALL mtx_set_vector((n+1)/4,1.d0)
      CASE(2) 
         CALL mtx_set_vector(3*(n+1)/4,1.d0)
      CASE DEFAULT
         CALL mtx_set_vector((n+1)/2,1.d0)
      END SELECT
            
      CALL mtx_barrier

      CALL mtx_solve(its)
      if(irank.eq.0) then
         write(6,*) 'Iteration Number=',its
      endif

      DO i=istart,iend
         CALL mtx_get_vector(i,v)
         WRITE(6,'(A,I5,1PE12.4)') 'i,v=',i,v
      ENDDO
      
      CALL mtx_cleanup
      GO TO 1

 9000 CALL mtx_finalize
      STOP
      END PROGRAM testksp
