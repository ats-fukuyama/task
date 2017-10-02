!     $Id$

!     Test program of libmtx for 1D/2D/3D Poisson equation

!     Input parameters
!        idimen : number of dimension (i or 2 or 3),  0 for quit
!        isiz : number of mesh point in one dimension
!        isource : source position (isource, isource, isource)
!        itype: type of linear solver (0 for default)
!        m1: methodKSP 0..13 (default=4)
!        m2: methodPC  0..12 (default=5)
!        tolerance : tolerance in iterative method

      PROGRAM testmtx

      USE libmpi 
      USE commpi
      USE libmtx 
      IMPLICIT NONE 
      INTEGER:: idimen,isiz,isource,itype,m1,m2,idebug
      INTEGER:: istart,iend,its 
      INTEGER:: imax,jwidth,jsource 
      INTEGER:: i,j,k,l,m,n,iskip,ncom 
      REAL(8):: v,tolerance 
      REAL(8),DIMENSION(:),POINTER:: x
      INTEGER,DIMENSION(7):: idata 
      REAL(8),DIMENSION(1):: ddata
      REAL(4):: cputime1,cputime2

      CALL mtx_initialize
      idimen=1
      isiz=11
      isource=6
      itype=0
      m1=4
      IF(nsize == 1) THEN
         m2=5
      ELSE
         m2=0
      ENDIF
      tolerance=1.d-7
      idebug=0

    1 CONTINUE
      IF(nrank.eq.0) then
    2    WRITE(6,'(A/6I5,1PE12.4,I3)') &
              '# INPUT: idimen,isiz,isource,itype,m1,m2,tolerance,idebug=', &
                        idimen,isiz,isource,itype,m1,m2,tolerance,idebug
         READ(5,*,END=3,ERR=2) idimen,isiz,isource,itype,m1,m2,tolerance,idebug
         idata(1)=idimen
         idata(2)=isiz
         idata(3)=isource
         idata(4)=itype
         idata(5)=m1
         idata(6)=m2
         idata(7)=idebug
         ddata(1)=tolerance
         IF(idimen.LT.0.OR.idimen.GT.3) THEN
            WRITE(6,*) 'XX idimen: out of range'
            GO TO 2
         ENDIF
         GO TO 4
    3    idata(1)=0
    4    CONTINUE
      ENDIF
      CALL mtx_broadcast_integer(idata,7)
      CALL mtx_broadcast_real8(ddata,1)
      idimen=idata(1)
      isiz=idata(2)
      isource=idata(3)
      itype=idata(4)
      m1=idata(5)
      m2=idata(6)
      idebug=idata(7)
      tolerance=ddata(1)

      IF(idimen.EQ.0) GO TO 9000

      IF(nrank.EQ.0) CALL CPU_TIME(cputime1)

      SELECT CASE(idimen)
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

      CALL mtx_setup(imax,istart,iend,idebug=idebug)

      SELECT CASE(idimen)
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
         CALL mtx_set_source(i,0.d0)
      ENDDO
      IF(jsource.GE.istart.AND.jsource.LE.iend) &
           CALL mtx_set_source(jsource,-1.d0)

      CALL mtx_solve(itype,tolerance,its,methodKSP=m1,methodPC=m2)
      if(nrank.eq.0) then
         write(6,*) 'Iteration Number=',its
      endif

      CALL mtx_gather_vector(x)

      IF(nrank.eq.0) THEN
         SELECT CASE(idimen)
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
!      CALL mtx_barrier

      CALL mtx_cleanup
      DEALLOCATE(x)

      IF(nrank.eq.0) THEN
         CALL CPU_TIME(cputime2)
         write(6,'(A,F12.3)') &
              '--cpu time =',cputime2-cputime1
      ENDIF
      GO TO 1

 9000 CALL mtx_finalize
      STOP
      END PROGRAM testmtx

