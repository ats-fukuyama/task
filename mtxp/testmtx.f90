      PROGRAM testmtx

      USE libmtx
      IMPLICIT NONE
      INTEGER:: idim,ilen,isource
      INTEGER:: nrank,nprocs,istart,iend,its
      INTEGER:: imax,jwidth,jsource
      INTEGER:: i,j,k,l,m,n,iskip
      INTEGER:: nblock,nblock_per_proc,nblock_addition
      REAL(8):: v
      REAL(8),DIMENSION(:),POINTER:: x
      INTEGER,DIMENSION(2):: idata

      CALL mtx_initialize(nrank,nprocs)
      idim=1
      ilen=11
      isource=6

    1 IF(nrank.eq.0) then
         WRITE(6,'(A,3I5)') '# INPUT: idim,ilen,isource=', &
                                      idim,ilen,isource
         READ(5,*,END=9000,ERR=1) idim,ilen,isource
         idata(1)=idim
         idata(2)=ilen
         idata(3)=isource
         IF(idim.LT.0.OR.idim.GT.3) THEN
            WRITE(6,*) 'XX idim: out of range'
            GO TO 1
         ENDIF
      ENDIF
      CALL mtx_broadcast_integer(idata,3)
      idim=idata(1)
      ilen=idata(2)
      isource=idata(3)

      IF(idim.EQ.0) GO TO 9000

      SELECT CASE(idim)
      CASE(1)
         imax=ilen
         jwidth=3
         jsource=isource
      CASE(2)
         imax=ilen*ilen
         jwidth=4*ilen-1
         jsource=ilen*(isource-1)+isource
      CASE(3)
         imax=ilen*ilen*ilen
         jwidth=4*ilen*ilen-1
         jsource=ilen*(ilen*(isource-1)+isource-1)+isource
      END SELECT
      ALLOCATE(x(imax))

      CALL mtx_setup(imax,istart,iend,jwidth)

      SELECT CASE(idim)
      CASE(1)
         DO i=istart,iend
            l=i
            if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
            CALL mtx_set_matrix(i,i,-2.d0)
            if(l.lt.ilen) CALL mtx_set_matrix(i,i+1,1.d0)
         ENDDO
      CASE(2)
         DO i=istart,iend
            l=mod(i-1,ilen)+1
            m=(i-1)/ilen+1
            if(m.gt.1) CALL mtx_set_matrix(i,i-ilen,1.d0)
            if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
            CALL mtx_set_matrix(i,i,-4.d0)
            if(l.lt.ilen) CALL mtx_set_matrix(i,i+1,1.d0)
            if(m.lt.ilen) CALL mtx_set_matrix(i,i+ilen,1.d0)
         ENDDO
      CASE(3)
         DO i=istart,iend
            l=mod(i-1,ilen)+1
            n=(i-1)/ilen+1
            m=mod(n-1,ilen)+1
            n=(n-1)/ilen+1
            if(n.gt.1) CALL mtx_set_matrix(i,i-ilen*ilen,1.d0)
            if(m.gt.1) CALL mtx_set_matrix(i,i-ilen,1.d0)
            if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
            CALL mtx_set_matrix(i,i,-6.d0)
            if(l.lt.ilen) CALL mtx_set_matrix(i,i+1,1.d0)
            if(m.lt.ilen) CALL mtx_set_matrix(i,i+ilen,1.d0)
            if(n.lt.ilen) CALL mtx_set_matrix(i,i+ilen*ilen,1.d0)
         ENDDO
      END SELECT

      DO i=istart,iend
         CALL mtx_set_vector(i,0.d0)
      ENDDO
      CALL mtx_set_vector(jsource,-1.d0)
      CALL mtx_barrier

      CALL mtx_solve(its)
      if(nrank.eq.0) then
         write(6,*) 'Iteration Number=',its
      endif

      CALL mtx_gather_vector(x)

      IF(nrank.eq.0) THEN
         SELECT CASE(idim)
         CASE(1)
            DO i=1,ilen
               WRITE(6,'(A,I5,1PE12.4)') 'i,x=',i,x(i)
            ENDDO
         CASE(2)
            iskip=MAX(ilen/10,1)
            DO i=1,ilen,iskip
               DO j=1,ilen,iskip
                  WRITE(6,'(A,2I5,1PE12.4)') 'i,j,x=',i,j,x(ilen*(i-1)+j)
               ENDDO
            ENDDO
         CASE(3)
            iskip=MAX(ilen/5,1)
            DO i=1,ilen,iskip
               DO j=1,ilen,iskip
                  DO k=1,ilen,iskip
                     WRITE(6,'(A,3I5,1PE12.4)') 'i,j,k,x=', &
                     i,j,k,x(ilen*(ilen*(i-1)+j-1)+k)
                  ENDDO
               ENDDO
            ENDDO
         END SELECT
      ENDIF

      CALL mtx_cleanup
      IF(nrank.EQ.0) DEALLOCATE(x)

      GO TO 1

 9000 CALL mtx_finalize
      STOP
      END PROGRAM testmtx
