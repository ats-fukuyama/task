! testlamb.f90

PROGRAM testlamb

  USE task_kinds,ONLY: dp
  USE libbes
  USE libgrf
  IMPLICIT NONE
  COMPLEX(rkind):: ca
  COMPLEX(rkind),ALLOCATABLE:: cl(:)
  REAL(rkind),ALLOCATABLE:: xdata(:),fdata(:,:),f2data(:,:)
  REAL(rkind):: xmin,xmax,dx,x
  INTEGER:: nmin,nmax,nxmax,n,nx,ierr

  CALL GSOPEN
  nmin=0
  nmax=9
  xmin=0.D0
  xmax=12.D0
  nxmax=201

1 CONTINUE
  WRITE(6,'(A,3I5,2ES12.4)') &
       '## nxmax,nmin,nmax,xmin,xmax=',nxmax,nmin,nmax,xmin,xmax
  READ(5,*,ERR=1,END=9000) nxmax,nmin,nmax,xmin,xmax
  IF(nxmax.LE.0) GOTO 9000
  IF(nxmax.EQ.1) GOTO 1

  ALLOCATE(xdata(nxmax),fdata(nxmax,nmin:nmax),cl(0:nmax))
  ALLOCATE(f2data(nxmax,nmin:nmax))
  dx=(xmax-xmin)/(nxmax-1)
  DO nx=1,nxmax
     xdata(nx)=xmin+dx*(nx-1)
  END DO
  DO n=nmin,nmax
     DO nx=1,nxmax
        ca=xdata(nx)
        CALL LAMBDA(nmax,ca,cl,ierr)
        IF(ierr.NE.0) THEN
           WRITE(6,'(A,I5,3ES12.4)') 'XX Error in LAMBDA: ierr,ca=',ierr,ca
           EXIT
        END IF
        fdata(nx,n)=REAL(cl(n))
     END DO
     DO nx=1,nxmax
        x=xdata(nx)
        f2data(nx,n)=BESEINX(n,x)
     END DO
  END DO

  CALL PAGES
  CALL GRD1D(1,xdata,fdata,nxmax,nxmax,nmax-nmin*1,'@LAMBDA(x)@',0)
  CALL GRD1D(2,xdata,f2data,nxmax,nxmax,nmax-nmin*1,'@BESEINX(x)@',0)
  CALL PAGEE
  DEALLOCATE(xdata,fdata,f2data,cl)
  GOTO 1

9000 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM testlamb
