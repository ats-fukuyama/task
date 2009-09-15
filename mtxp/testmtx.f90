!     $Id

!     Test program of libmtx for 1D/2D/3D Poisson equation

!     Input parameters
!        idim : number of dimension (i or 2 or 3),  0 for quit
!        isiz : number of mesh point in one dimension
!        isource : source position (isource, isource, isource)
!        itype: type of linear solver (0 for default)
!        tolerance : tolerance in iterative method

      PROGRAM testmtx

      USE libmtx
      IMPLICIT NONE
      INTEGER:: idim,isiz,isource,itype
      INTEGER:: nrank,nprocs,istart,iend,its
      INTEGER:: imax,jwidth,jsource
      INTEGER:: i,j,k,l,m,n,iskip
      INTEGER:: nblock,nblock_per_proc,nblock_addition
      REAL(8):: v,tolerance
      REAL(8),DIMENSION(:),POINTER:: x
      INTEGER,DIMENSION(2):: idata
      REAL(4):: cputime1,cputime2

      CALL mtx_initialize(nrank,nprocs)
      idim=1
      isiz=11
      isource=6
      itype=0
      tolerance=1.d-7

    1 CONTINUE
      IF(nrank.eq.0) then
    2    WRITE(6,'(A,4I5,1PE12.4)') &
              '# INPUT: idim,isiz,isource,itype,tolerance=', &
                        idim,isiz,isource,itype,tolerance
         READ(5,*,END=3,ERR=2) idim,isiz,isource,itype,tolerance
         idata(1)=idim
         idata(2)=isiz
         idata(3)=isource
         IF(idim.LT.0.OR.idim.GT.3) THEN
            WRITE(6,*) 'XX idim: out of range'
            GO TO 2
         ENDIF
         GO TO 4
    3    idata(1)=0
    4    CONTINUE
      ENDIF
      CALL mtx_broadcast_integer(idata,3)
      idim=idata(1)
      isiz=idata(2)
      isource=idata(3)

      IF(idim.EQ.0) GO TO 9000

      IF(nrank.EQ.0) CALL CPU_TIME(cputime1)

      SELECT CASE(idim)
      CASE(1)
         imax=isiz
         jwidth=3
         jsource=isource
      CASE(2)
         imax=isiz*isiz
         jwidth=4*isiz-1
         jsource=isiz*(isource-1)+isource
      CASE(3)
         imax=isiz*isiz*isiz
         jwidth=4*isiz*isiz-1
         jsource=isiz*(isiz*(isource-1)+isource-1)+isource
      END SELECT
      ALLOCATE(x(imax))

      CALL mtx_setup(imax,istart,iend,jwidth)

      SELECT CASE(idim)
      CASE(1)
         DO i=istart,iend
            l=i
            if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
            CALL mtx_set_matrix(i,i,-2.d0)
            if(l.lt.isiz) CALL mtx_set_matrix(i,i+1,1.d0)
         ENDDO
      CASE(2)
         DO i=istart,iend
            l=mod(i-1,isiz)+1
            m=(i-1)/isiz+1
            if(m.gt.1) CALL mtx_set_matrix(i,i-isiz,1.d0)
            if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
            CALL mtx_set_matrix(i,i,-4.d0)
            if(l.lt.isiz) CALL mtx_set_matrix(i,i+1,1.d0)
            if(m.lt.isiz) CALL mtx_set_matrix(i,i+isiz,1.d0)
         ENDDO
      CASE(3)
         DO i=istart,iend
            l=mod(i-1,isiz)+1
            n=(i-1)/isiz+1
            m=mod(n-1,isiz)+1
            n=(n-1)/isiz+1
            if(n.gt.1) CALL mtx_set_matrix(i,i-isiz*isiz,1.d0)
            if(m.gt.1) CALL mtx_set_matrix(i,i-isiz,1.d0)
            if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
            CALL mtx_set_matrix(i,i,-6.d0)
            if(l.lt.isiz) CALL mtx_set_matrix(i,i+1,1.d0)
            if(m.lt.isiz) CALL mtx_set_matrix(i,i+isiz,1.d0)
            if(n.lt.isiz) CALL mtx_set_matrix(i,i+isiz*isiz,1.d0)
         ENDDO
      END SELECT

      DO i=istart,iend
         CALL mtx_set_vector(i,0.d0)
      ENDDO
      CALL mtx_set_vector(jsource,-1.d0)
      CALL mtx_barrier

      CALL mtx_solve(itype,tolerance,its)
      if(nrank.eq.0) then
         write(6,*) 'Iteration Number=',its
      endif

      CALL mtx_gather_vector(x)

      IF(nrank.eq.0) THEN
         SELECT CASE(idim)
         CASE(1)
            DO i=1,isiz
               WRITE(6,'(A,I5,1PE12.4)') 'i,x=',i,x(i)
            ENDDO
         CASE(2)
            iskip=MAX(isiz/10,1)
            DO i=1,isiz,iskip
               DO j=1,isiz,iskip
                  WRITE(6,'(A,2I5,1PE12.4)') 'i,j,x=',i,j,x(isiz*(i-1)+j)
               ENDDO
            ENDDO
         CASE(3)
            iskip=MAX(isiz/5,1)
            DO i=1,isiz,iskip
               DO j=1,isiz,iskip
                  DO k=1,isiz,iskip
                     WRITE(6,'(A,3I5,1PE12.4)') 'i,j,k,x=', &
                     i,j,k,x(isiz*(isiz*(i-1)+j-1)+k)
                  ENDDO
               ENDDO
            ENDDO
         END SELECT
      ENDIF

      CALL mtx_cleanup
      IF(nrank.EQ.0) DEALLOCATE(x)

      IF(nrank.eq.0) THEN
         CALL CPU_TIME(cputime2)
         write(6,'(A,F12.3)') &
              '--cpu time =',cputime2-cputime1
      ENDIF
      GO TO 1

 9000 CALL mtx_finalize
      STOP
      END PROGRAM testmtx

