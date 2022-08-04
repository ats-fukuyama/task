! wqsolv2.f90

MODULE wqsolv2

  PRIVATE
  PUBLIC wq_solv2

CONTAINS
  
  subroutine wq_solv2

    !  solve equation: A (X_n+1 - X_n)/dt = M ( f X_n+1 +(1-f) X_n)
    !                  (A - f dt M) X_n+1 = (A + (1-f) dt M) X_n
    !                        f: implicit parameter
    !                           0.0 : explicit
    !                           0.5 : Cranck-Nicholson
    !                           1.0 : full implicit
    !                        cvec_l=X_n
    !                        cvec_r=(A + (1-f) dt M) X_n
    !                        cvec_s=X_n+1

    USE wqcomm
    USE libmpi
    USE libmtx
    implicit none

    INTEGER :: nx,ny,i,j,i0,ii0,il,jl,mblock_size,mlen,mbnd
    REAL(rkind)    :: kz,diag_factor,dx,dy
    INTEGER:: istart,iend,nx_start,nx_end
    INTEGER:: n,nxl,nyl
    INTEGER,ALLOCATABLE:: istart_nrank(:),iend_nrank(:)
    INTEGER,ALLOCATABLE:: nx_start_nrank(:),nx_end_nrank(:)
    COMPLEX(rkind) :: P,Q,R1xy,R2xx,R2yy,R3x,R3y
    COMPLEX(rkind):: CM(3,3,-1:1,-1:1)
    COMPLEX(rkind),ALLOCATABLE:: cvec_l(:),cvec_r(:),cvec_s(:)
    COMPLEX(rkind):: cmat_r,cmat_l

    SELECT CASE(model_geometry)
    CASE(0)
       kz=rkz
    CASE(1)
       kz=nph/RR
    END SELECT


    P  = omega**2/VC**2
    Q  = CI*omega*RMU0

    mblock_size=3*nymax
    mlen=nxmax*mblock_size
    mbnd=4*mblock_size-1

    CALL mtxc_setup(mlen,istart,iend,jwidth=mbnd)

    nx_start=(istart-1)/mblock_size+1    ! nx_start is nx including istart
    nx_end=(iend-1)/mblock_size+1        ! nx_end is nx including iend
    IF(nrank.EQ.nsize-1) nx_end=nxmax

    IF(nrank.EQ.0) &
         WRITE(6,'(A,5I6)') 'mtx_setup: mle,is,ie,nxs,nxe=', &
         mlen,istart,iend,nx_start,nx_end

    IF(ALLOCATED(istart_nrank)) &
         DEALLOCATE(istart_nrank,iend_nrank,nx_start_nrank,nx_end_nrank)
    ALLOCATE(istart_nrank(0:nsize-1),iend_nrank(0:nsize-1))
    ALLOCATE(nx_start_nrank(0:nsize-1),nx_end_nrank(0:nsize-1))
    CALL mtx_allgather1_integer(istart,istart_nrank)
    CALL mtx_allgather1_integer(iend,iend_nrank)
    CALL mtx_allgather1_integer(nx_start,nx_start_nrank)
    CALL mtx_allgather1_integer(nx_end,nx_end_nrank)

    IF(idebuga(41).EQ.1) THEN
       DO n=0,nsize-1
          WRITE(6,'(A,6I8)') 'nrank,is,ie,nxs,nxe,nxl=', &
               n,istart_nrank(n),iend_nrank(n), &
               nx_start_nrank(n),nx_end_nrank(n), &
               nx_end_nrank(n)-nx_start_nrank(n)+1
       END DO
    END IF

    CM(1:3,1:3,-1:1,-1:1)=(0.D0,0.D0)
    ALLOCATE(cvec_r(mlen),cvec_l(mlen),cvec_s(mlen))
    cvec_r(1:mlen)=0.D0

    do nx=1,nxmax
       do ny=1,nymax
          i=3*(nymax*(nx-1)+ny-1)
          cvec_l(i+1)=EX(nx,ny)
          cvec_l(i+2)=EY(nx,ny)
          cvec_l(i+3)=EZ(nx,ny)
       end do
    end do

    DO i0=istart,iend
       nx=(i0-1)/(3*nymax)+1   ! 1..nxmax
       ii0=i0-3*nymax*(nx-1)    ! 1..3*nymax
       ny=(ii0-1)/3+1          ! 1..nymax
       il=ii0-3*(ny-1)

       IF(nx.GE.2.AND.nx.LE.nxmax-1.AND.ny.GE.2.AND.ny.LE.nymax-1) THEN
          dx=0.5D0*(xg_nx(nx+1)-xg_nx(nx-1))
          dy=0.5D0*(yg_ny(ny+1)-yg_ny(ny-1))

          R1xy = 0.25d0/(dx*dy)
          R2xx = 1.00d0/dx**2
          R2yy = 1.00d0/dy**2
          R3x = 0.50d0*CI*kz/dx
          R3y = 0.50d0*CI*kz/dy

          CM(1,2,+1,+1)= R1xy
          CM(1,2,+1,-1)=-R1xy
          CM(1,2,-1,+1)=-R1xy
          CM(1,2,-1,-1)= R1xy
          CM(1,3,+1, 0)= R3x
          CM(1,3,-1, 0)=-R3x
          CM(1,1, 0,+1)=-     R2yy
          CM(1,1, 0, 0)=+2.D0*R2yy+kz**2-P-Q*CD(1,1,nx,ny)
          CM(1,1, 0,+1)=-     R2yy
          CM(1,2, 0, 0)=                  -Q*CD(1,2,nx,ny)
          CM(1,3, 0, 0)=                  -Q*CD(1,3,nx,ny)
        
          CM(2,3, 0,+1)= R3y
          CM(2,3, 0,-1)=-R3y
          CM(2,1,+1,+1)= R1xy
          CM(2,1,+1,-1)=-R1xy
          CM(2,1,-1,+1)=-R1xy
          CM(2,1,-1,-1)= R1xy
          CM(2,2,+1, 0)=-     R2xx
          CM(2,2, 0, 0)=+2.D0*R2xx+kz**2-P-Q*CD(2,2,nx,ny)
          CM(2,2,-1, 0)=-     R2xx
          CM(2,1, 0, 0)=                  -Q*CD(2,1,nx,ny)
          CM(2,3, 0, 0)=                  -Q*CD(2,3,nx,ny)

          CM(3,1,+1, 0)= R3x
          CM(3,1,-1, 0)=-R3x
          CM(3,2, 0,+1)= R3y
          CM(3,2, 0,-1)=-R3y
          CM(3,3,+1, 0)=-     R2xx
          CM(3,3, 0,+1)=-               R2yy
          CM(3,3, 0, 0)=+2.D0*R2xx+2.D0*R2yy-P-Q*CD(3,3,nx,ny)
          CM(3,3,-1, 0)=-     R2xx
          CM(3,3, 0,-1)=-               R2yy
          CM(3,1, 0, 0)=                      -Q*CD(3,1,nx,ny)
          CM(3,2, 0, 0)=                      -Q*CD(3,2,nx,ny)

          DO nxl=-1,1
             DO nyl=-1,1
                IF(nxl.EQ.0.AND.nyl.EQ.0) THEN
                   diag_factor=1.D0
                ELSE
                   diag_factor=0.D0
                END IF
                i=3*(nymax*(nx-1    )+ny-1    )+il
                DO jl=1,3
                   j=3*(nymax*(nx-1+nxl)+ny-1+nyl)+jl
                   IF(j.GE.1.AND.j.LE.mlen) THEN
                      cmat_l=A(il,jl,nx,ny)*diag_factor &
                            -fimplicit*dt*CM(il,jl,nxl,nyl)
                      cmat_r=A(il,jl,nx,ny)*diag_factor &
                            +(1.D0-fimplicit)*dt*CM(il,jl,nxl,nyl)
                      IF(ABS(cmat_l).GT.0.D0) THEN
                         CALL mtxc_set_matrix(j,i,cmat_l)
                         IF(idebuga(61).EQ.1.AND.nrank.EQ.0) &
                              WRITE(61,'(A,2I6,2ES12.4)') 'cmat_l:',i,j,cmat_l
                      END IF
                      cvec_r(j)=cvec_r(j)+cmat_r*cvec_l(i)
                   END IF
                END DO
             END DO
          END DO
       ELSE  ! boundary (nx=1 or nxmax: Ex'=0,Ey=Ez=0)
             !          (ny=1 or nymax: Ey'=0,Ex=Ez=0)
          i=i0
          j=i0
          cmat_l=1.D0
          CALL mtxc_set_matrix(j,i,cmat_l)
          IF(idebuga(61).EQ.1.AND.nrank.EQ.0) &
               WRITE(61,'(A,2I6,2ES12.4)') 'cmat_l:',i,j,cmat_l
          jl=il
          nxl=0
          nyl=0
          cmat_r=0.D0
          SELECT CASE(il)
          CASE(1) ! Ex
             IF(nx.EQ.1) THEN
                nxl=1
                cmat_l=-1.D0
             ELSE IF(nx.EQ.nxmax) THEN
                nxl=-1
                cmat_l=-1.D0
             END IF
          CASE(2) ! Ey
             IF(ny.EQ.1) THEN
                nyl=1
                cmat_l=-1.D0
             ELSE IF(ny.EQ.nymax) THEN
                nyl=-1
                cmat_l=-1.D0
             END IF
          END SELECT
          IF(ABS(cmat_l).GT.0.D0) THEN
             j=3*(nymax*(nx-1+nxl)+ny-1+nyl)+jl
             CALL mtxc_set_matrix(j,i,cmat_l)
             IF(idebuga(61).EQ.1.AND.nrank.EQ.0) &
                  WRITE(61,'(A,2I6,2ES12.4)') 'cmat_l:',i,j,cmat_l
          END IF
       END IF
    END DO

!   --- setup right-hand-side vector ---

    DO i=istart,iend
       IF(ABS(cvec_r(i)).GT.0.D0) CALL mtxc_set_source(i,cvec_r(i))
       CALL mtxc_set_source(i,cvec_r(i))
    END DO

!   --- setup initial solution vector for iterative scheme of matrix solver ---

    DO i=istart,iend
       IF(ABS(cvec_l(i)).GT.0.D0) CALL mtxc_set_vector(i,cvec_l(i))
    END DO

    IF(idebuga(62).NE.0) THEN
       DO i=istart,iend
          IF(ABS(cvec_r(i)).GT.0.D0.OR.ABS(cvec_l(i)).GT.0.D0) &
               WRITE(62,'(A,I8,4ES12.4)') &
               'cvec_r,cvec_l:',i,cvec_r(i),cvec_l(i)
       END DO
    END IF
    
!   --- solve matrix equation  ---

   CALL mtxc_solve(ntype_mat,eps_mat,icount_mat)

!   --- obtain solution vector  ---

    CALL mtxc_gather_vector(cvec_s)

!   --- clear up matrix solver  ---

    CALL mtxc_cleanup

    IF(idebuga(63).NE.0) THEN
       WRITE(63,'(A)') '### cvec_l,cvec_s'
       DO i=istart,iend
          IF(ABS(cvec_l(i)).GT.0.D0.OR.ABS(cvec_s(i)).GT.0.D0) &
          WRITE(63,'(A,I8,4ES12.4)') 'cvec_l,cvec_s:',i,cvec_l(i),cvec_s(i)
       END DO
    END IF
    
    do ny=1,nymax
       do nx=1,nxmax
          i=3*(nymax*(nx-1)+ny-1)
          EX(nx,ny)=cvec_s(i+1)
          EY(nx,ny)=cvec_s(i+2)
          EZ(nx,ny)=cvec_s(i+3)
       end do
    end do

  return
end subroutine wq_solv2
END MODULE wqsolv2
