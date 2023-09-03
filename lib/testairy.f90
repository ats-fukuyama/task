! testairy.f90

PROGRAM testairy

  USE task_kinds,ONLY: dp
  USE libspf2,ONLY: airyb
  USE libgrf
  IMPLICIT NONE
  COMPLEX(dp):: ca
  REAL(dp),ALLOCATABLE:: xdata(:)
  REAL(dp),ALLOCATABLE:: Aidata(:),dAidata(:)
  REAL(dp),ALLOCATABLE:: Bidata(:),dBidata(:)
  REAL(dp):: xmin,xmax,dx,x
  INTEGER:: itype,nxmax,nx,ierr

  CALL GSOPEN
  itype=1
  xmin=-6.D0
  xmax=+6.D0
  nxmax=201

1 CONTINUE
  WRITE(6,'(A,2I4,2ES12.4)') &
       '## itype[0:1],nxmax,xmin,xmax=', &
           itype,nxmax,xmin,xmax
  READ(5,*,ERR=1,END=9000) itype,nxmax,xmin,xmax
  IF(itype.LE.0) GOTO 9000
  IF(nxmax.EQ.1) GOTO 1

  ALLOCATE(xdata(nxmax))
  ALLOCATE(Aidata(nxmax),dAidata(nxmax))
  ALLOCATE(Bidata(nxmax),dBidata(nxmax))

  dx=(xmax-xmin)/(nxmax-1)
  DO nx=1,nxmax
     xdata(nx)=xmin+dx*(nx-1)
  END DO

  DO nx=1,nxmax
     CALL airyb(xdata(nx),Aidata(nx),Bidata(nx), &
                          dAidata(nx),dBidata(nx))
  END DO
  CALL PAGES
  CALL GRD1D(1,xdata,Aidata,nxmax,nxmax,1,'@Ai(x)@',0)
  CALL GRD1D(2,xdata,dAidata,nxmax,nxmax,1,'@dAi/dx(x)@',0)
  CALL GRD1D(3,xdata,Bidata,nxmax,nxmax,1,'@Bi(x)@',0)
  CALL GRD1D(4,xdata,dBidata,nxmax,nxmax,1,'@dBi/dx(x)@',0)
  CALL PAGEE
  DEALLOCATE(xdata)
  DEALLOCATE(Aidata,dAidata,Bidata,dBidata)
  GOTO 1

9000 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM testairy
