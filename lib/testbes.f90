! testbes.f90

PROGRAM testmbes

  USE task_kinds,ONLY: dp
  USE libbes
  USE libgrf
  IMPLICIT NONE
  COMPLEX(dp):: ca
  COMPLEX(dp),ALLOCATABLE:: cl(:)
  REAL(dp),ALLOCATABLE:: xdata(:),fdata(:,:),vl(:)
  REAL(dp),ALLOCATABLE:: fjdata(:,:),fydata(:,:)
  REAL(dp),ALLOCATABLE:: fidata(:,:),fkdata(:,:)
  REAL(dp),ALLOCATABLE:: eidata(:,:),ekdata(:,:)
  REAL(dp):: xmin,xmax,dx,x
  INTEGER:: itype,nmin,nmax,nxmax,n,nx,ierr

  CALL GSOPEN
  itype=1
  nmin=0
  nmax=9
  xmin=0.D0
  xmax=12.D0
  nxmax=201
  WRITE(6,'(A)') '## itype: 0:end,1,3:J,Y,I,K, 2,4:J,Y,EI,EK, 5:LAMBDA,EI'

1 CONTINUE
  WRITE(6,'(A,4I4,2ES12.4)') &
       '## itype[0:3],nxmax,nmin,nmax,xmin,xmax=', &
           itype,nxmax,nmin,nmax,xmin,xmax
  READ(5,*,ERR=1,END=9000) itype,nxmax,nmin,nmax,xmin,xmax
  IF(itype.LE.0) GOTO 9000
  IF(nxmax.EQ.1) GOTO 1

  ALLOCATE(xdata(nxmax),cl(0:nmax),fdata(nxmax,nmin:nmax),vl(0:nmax))
  ALLOCATE(fjdata(nxmax,nmin:nmax),fydata(nxmax,nmin:nmax))
  ALLOCATE(fidata(nxmax,nmin:nmax),fkdata(nxmax,nmin:nmax))
  ALLOCATE(eidata(nxmax,nmin:nmax),ekdata(nxmax,nmin:nmax))

  dx=(xmax-xmin)/(nxmax-1)
  DO nx=1,nxmax
     xdata(nx)=xmin+dx*(nx-1)
  END DO

  SELECT CASE(itype)
  CASE(1)
     DO n=nmin,nmax
        DO nx=1,nxmax
           x=xdata(nx)
           fjdata(nx,n)=BESJNX(n,x)
           fydata(nx,n)=BESYNX(n,x)
           fidata(nx,n)=BESINX(n,x)
           fkdata(nx,n)=BESKNX(n,x)
        END DO
     END DO
     CALL PAGES
     CALL GRD1D(1,xdata,fjdata,nxmax,nxmax,nmax-nmin+1,'@BESJNX(x)@',0)
     CALL GRD1D(2,xdata,fydata,nxmax,nxmax,nmax-nmin+1,'@BESYNX(x)@',0)
     CALL GRD1D(3,xdata,fidata,nxmax,nxmax,nmax-nmin+1,'@BESINX(x)@',0)
     CALL GRD1D(4,xdata,fkdata,nxmax,nxmax,nmax-nmin+1,'@BESKNX(x)@',0)
     CALL PAGEE
  CASE(2)
     DO n=nmin,nmax
        DO nx=1,nxmax
           x=xdata(nx)
           fjdata(nx,n)=BESJNX(n,x)
           fydata(nx,n)=BESYNX(n,x)
           eidata(nx,n)=BESEINX(n,x)
           ekdata(nx,n)=BESEKNX(n,x)
        END DO
     END DO
     CALL PAGES
     CALL GRD1D(1,xdata,fjdata,nxmax,nxmax,nmax-nmin+1,'@BESJNX(x)@',0)
     CALL GRD1D(2,xdata,fydata,nxmax,nxmax,nmax-nmin+1,'@BESYNX(x)@',0)
     CALL GRD1D(3,xdata,eidata,nxmax,nxmax,nmax-nmin+1,'@BESEINX(x)@',0)
     CALL GRD1D(4,xdata,ekdata,nxmax,nxmax,nmax-nmin+1,'@BESEKNX(x)@',0)
     CALL PAGEE
  CASE(3)
     DO nx=1,nxmax
        x=xdata(nx)
        CALL BESJNV(nmax,x,vl,ierr)
        fjdata(nx,nmin:nmax)=vl(nmin:nmax)
        CALL BESYNV(nmax,x,vl,ierr)
        fydata(nx,nmin:nmax)=vl(nmin:nmax)
        CALL BESINV(nmax,x,vl,ierr)
        fidata(nx,nmin:nmax)=vl(nmin:nmax)
        CALL BESKNV(nmax,x,vl,ierr)
        fkdata(nx,nmin:nmax)=vl(nmin:nmax)
     END DO
     CALL PAGES
     CALL GRD1D(1,xdata,fjdata,nxmax,nxmax,nmax-nmin+1,'@BESJNX(x)@',0)
     CALL GRD1D(2,xdata,fydata,nxmax,nxmax,nmax-nmin+1,'@BESYNX(x)@',0)
     CALL GRD1D(3,xdata,fidata,nxmax,nxmax,nmax-nmin+1,'@BESINX(x)@',0)
     CALL GRD1D(4,xdata,fkdata,nxmax,nxmax,nmax-nmin+1,'@BESKNX(x)@',0)
     CALL PAGEE
  CASE(4)
     DO nx=1,nxmax
        x=xdata(nx)
        CALL BESJNV(nmax,x,vl,ierr)
        fjdata(nx,nmin:nmax)=vl(nmin:nmax)
        CALL BESYNV(nmax,x,vl,ierr)
        fydata(nx,nmin:nmax)=vl(nmin:nmax)
        CALL BESEINV(nmax,x,vl,ierr)
        eidata(nx,nmin:nmax)=vl(nmin:nmax)
        CALL BESEKNV(nmax,x,vl,ierr)
        ekdata(nx,nmin:nmax)=vl(nmin:nmax)
     END DO
     CALL PAGES
     CALL GRD1D(1,xdata,fjdata,nxmax,nxmax,nmax-nmin+1,'@BESJNX(x)@',0)
     CALL GRD1D(2,xdata,fydata,nxmax,nxmax,nmax-nmin+1,'@BESYNX(x)@',0)
     CALL GRD1D(3,xdata,eidata,nxmax,nxmax,nmax-nmin+1,'@BESEINX(x)@',0)
     CALL GRD1D(4,xdata,ekdata,nxmax,nxmax,nmax-nmin+1,'@BESEKNX(x)@',0)
     CALL PAGEE
  CASE(5)
     DO nx=1,nxmax
        ca=xdata(nx)
        CALL LAMBDA(nmax,ca,cl,ierr)
        IF(ierr.NE.0) THEN
           WRITE(6,'(A,I5,3ES12.4)') 'XX Error in LAMBDA: ierr,ca=',ierr,ca
           EXIT
        END IF
        fdata(nx,nmin:nmax)=REAL(cl(nmin:nmax))
        x=xdata(nx)
        CALL BESEINV(nmax,x,vl,ierr)
        eidata(nx,nmin:nmax)=vl(nmin:nmax)
        ekdata(nx,nmin:nmax)=fdata(nx,nmin:nmax)-eidata(nx,nmin:nmax)
     END DO
     CALL PAGES
     CALL GRD1D(1,xdata,fdata,nxmax,nxmax,nmax-nmin+1,'@LAMBDA(x)@',0)
     CALL GRD1D(2,xdata,eidata,nxmax,nxmax,nmax-nmin+1,'@BESEINX(x)@',0)
     CALL GRD1D(3,xdata,ekdata,nxmax,nxmax,nmax-nmin+1,'@difference(x)@',0)
     CALL PAGEE
  END SELECT
  DEALLOCATE(xdata,fdata,cl,vl)
  DEALLOCATE(fjdata,fydata,fidata,fkdata,eidata,ekdata)
  GOTO 1

9000 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM testmbes
