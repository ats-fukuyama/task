! plotpkf.f90

PROGRAM plotpkf
  USE task_kinds
  USE task_constants
  USE libpdkf
  USE libpgkf
  USE libgrf
  IMPLICIT NONE

  INTEGER:: nxmax
  REAL(rkind):: xmin,xmax
  COMPLEX(rkind):: cf_
  REAL(rkind),ALLOCATABLE:: xa(:),fa(:,:)
  REAL(rkind):: alpha,beta,rky,dx
  INTEGER:: mode,n1,nn,nf,nm,nx
  CHARACTER(LEN=80):: title
  
  CALL GSOPEN

  nxmax=101
  xmin=0.D0
  xmax=30.D0
  alpha=0.D0
  beta=1.D0
  rky=0.D0
  n1=1
  nf=0
  nn=0
  nm=0
  mode=0

1 CONTINUE
  
  WRITE(6,'(A)') '## mode: 0:pdkf 1:pgkf, 9:QUIT'
  READ(5,*,ERR=1,END=9000) mode

  SELECT CASE(mode)
  CASE(0)

     WRITE(6,'(A)') '## n1=1:2'
10   CONTINUE
     WRITE(6,'(A/A,3ES12.4,I4,2ES12.4,I6)') &
          '## INPUT : alpha,beta,rky,n1,xmin,xmax,nxmax:', &
          '##   ', alpha,beta,rky,n1,xmin,xmax,nxmax
     READ(5,*,END=1) alpha,beta,rky,n1,xmin,xmax,nxmax

     IF(nxmax.LT.1) GO TO 1
     
     IF(n1.LT.1.OR.N1.GT.2) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pkrf: n1 must be in [1:2]: n1=',n1
        GO TO 10
     END IF

     ALLOCATE(xa(nxmax),fa(nxmax,2))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
        cf_=pdkf(xa(nx),alpha,beta,rky,5,n1)
        fa(nx,1)=REAL(cf_)
        fa(nx,2)=AIMAG(cf_)
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF vs x: n1=',n1,'@'
     CALL grd1d(0,xa,fa,nxmax,nxmax,2,title)
     CALL PAGEE
     DEALLOCATE(xa,fa)
     GO TO 10

  CASE(1)
     WRITE(6,'(A)') '## nf=0:4, xmax<=10.D0'
20   CONTINUE
     WRITE(6,'(A/A,3I4,2ES12.4,I6)') &
          '## INPUT : nf,nm,nn,xmin,xmax,nxmax:', &
          '##   ', nf,nm,nn,xmin,xmax,nxmax
     READ(5,*,END=1) nf,nm,nn,xmin,xmax,nxmax
     IF(xmax.gt.10.D0) xmax=10.D0

     IF(nxmax.LT.1) GO TO 1
     
     IF(nf.LT.0.OR.nf.GT.4) THEN
        WRITE(6,'(A,I4)') 'XX plotpkf: pkrf: nf must be in [0:4]: n4=',nf
        GO TO 20
     END IF

     ALLOCATE(xa(nxmax),fa(nxmax,2))

     dx=(xmax-xmin)/(nxmax-1)
     DO nx=1,nxmax
        xa(nx)=xmin+dx*(nx-1)
        cf_=pgkf(xa(nx),nf,nn,nm)
        fa(nx,1)=REAL(cf_)
        fa(nx,2)=AIMAG(cf_)
     END DO
  
     CALL PAGES
     WRITE(title,'(A,i1,A)') '@PDKF vs x: n1=',n1,'@'
     CALL grd1d(0,xa,fa,nxmax,nxmax,2,title)
     CALL PAGEE
     DEALLOCATE(xa,fa)
     GO TO 20
  CASE(9)
     GO TO 9000
  END SELECT
  GO TO 1

9000 CALL GSCLOS
  STOP
END PROGRAM plotpkf
