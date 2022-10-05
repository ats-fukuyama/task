MODULE libdeab

! XN=XHI*(X=XM):   -1.D0<XN<1.D0
! X =XH*XN+XM:     A < X < B

  USE libde
  USE bpsd_kinds,ONLY: dp

  REAL(dp):: XM_  ! mid XM=0.5D0*(A+B)
  REAL(dp):: XH_  ! half len XH=0.5D0*(B-A)
  REAL(dp):: XHI_ ! inverse of half len XHI=1.D0/XH

CONTAINS

  SUBROUTINE DEFTAB(A,B,FUNC,KID,H0,EPS,ILST,CS,ES)
    IMPLICIT NONE
    REAL(dp),INTENT(IN)::  A           ! Lower boundary of integral
    REAL(dp),INTENT(IN)::  B           ! Upper boundary of integral
    CHARACTER(LEN=*),INTENT(IN):: KID ! function identifier string
    REAL(dp),INTENT(IN)::  H0          ! Initial step size
    REAL(dp),INTENT(IN)::  EPS         ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST        ! print out control: 0 for no print out
    REAL(dp),INTENT(OUT):: CS          ! Integral
    REAL(dp),INTENT(OUT):: ES          ! Estimated error 
    INTERFACE
       FUNCTION FUNC(X,XM,XP)
         USE bpsd_kinds,ONLY: dp
         REAL(dp):: FUNC
         REAL(dp),INTENT(IN):: X,XM,XP
       END FUNCTION FUNC
    END INTERFACE

    XM_=0.5D0*(A+B)
    XH_=0.5D0*(B-A)
    XHI_=1.D0/XH_

    CALL DEFT(CS,ES,H0,EPS,ILST,FUNC,KID)

    CS=XH_*CS

    RETURN
  END SUBROUTINE DEFTAB
END MODULE libdeab

MODULE libde_wi

  USE libdeab
  REAL(dp):: gamma_,b_

CONTAINS

  FUNCTION FUNCAB1(X,XA,BX)
    USE libdeab
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN)::  X  ! Integral variable
    REAL(dp),INTENT(IN)::  XA ! X-A  near X=A
    REAL(dp),INTENT(IN)::  BX ! B-X  near X=B
    REAL(dp):: FUNCAB1,DUMMY

    DUMMY=XA
    DUMMY=BX
    FUNCAB1=(X-X*X/GAMMA_)*COS(X)

    RETURN
  END FUNCTION FUNCAB1
    
  FUNCTION FUNC1(X,XM,XP)
    USE libdeab
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN)::  X  ! Integral variable
    REAL(dp),INTENT(IN)::  XM ! 1.D0-X  near x= 1.D0
    REAL(dp),INTENT(IN)::  XP ! 1.D0+X  near X=-1.D0
    REAL(dp):: FUNC1

    FUNC1=FUNCAB1(XH_*X+XM_,XH_*XP,XH_*XM)

    RETURN
  END FUNCTION FUNC1

  FUNCTION FUNCAB2(X,XA,BX)
    USE libdeab
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN)::  X  ! Integral variable
    REAL(dp),INTENT(IN)::  XA ! X-A  near X=A
    REAL(dp),INTENT(IN)::  BX ! B-X  near X=B
    REAL(dp):: FUNCAB2,DUMMY

    DUMMY=XA
    DUMMY=BX
    FUNCAB2=(X-X*X/GAMMA_)*SIN(X)

    RETURN
  END FUNCTION FUNCAB2
    
  FUNCTION FUNC2(X,XM,XP)
    USE libdeab
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN)::  X  ! Integral variable
    REAL(dp),INTENT(IN)::  XM ! 1.D0-X  near x= 1.D0
    REAL(dp),INTENT(IN)::  XP ! 1.D0+X  near X=-1.D0
    REAL(dp):: FUNC2

    FUNC2=FUNCAB2(XH_*X+XM_,XH_*XP,XH_*XM)

    RETURN
  END FUNCTION FUNC2

  FUNCTION FUNCA(X)
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: X
    REAL(dp):: FUNCA,H0,EPS,CS,ES
    INTEGER:: ILST

    H0=0.25D0
    EPS=1.D-8
    ILST=0
    GAMMA_=X
    CALL DEFTAB(0.D0,GAMMA_,FUNC1,'FUNC1',H0,EPS,ILST,CS,ES)
    FUNCA=CS
  END FUNCTION FUNCA

  FUNCTION FUNCB(X)
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: X
    REAL(dp):: FUNCB,H0,EPS,CS,ES
    INTEGER:: ILST

    H0=0.25D0
    EPS=1.D-8
    ILST=0
    GAMMA_=X
    CALL DEFTAB(0.D0,GAMMA_,FUNC2,'FUNC2',H0,EPS,ILST,CS,ES)
    FUNCB=CS
  END FUNCTION FUNCB

  FUNCTION FUNC0(X)
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN)::  X  ! Integral variable
    REAL(dp):: FUNC0,v

    v=SQRT(2.D0*X)
    FUNC0=(FUNCA(b_*v)**2+FUNCB(b_*v))*EXP(-x)

    RETURN
  END FUNCTION FUNC0

  FUNCTION FUNCF(X)
    USE bpsd_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: X
    REAL(dp):: FUNCF,H0,EPS,CS,ES
    INTEGER:: ILST

    H0=0.25D0
    EPS=1.D-8
    ILST=0
    b_=X
    CALL DEHIFE(CS,ES,H0,EPS,0,FUNC0,'FF')
    FUNCF=CS
  END FUNCTION FUNCF

END MODULE libde_wi
  
    
PROGRAM plotwi

  USE libgrf
  USE libde_wi
  IMPLICIT NONE
  REAL(dp):: xmin,xmax,x,dx,ymin,ymax,dy,y
  INTEGER:: nxmax,nx
  REAL(dp),DIMENSION(:),ALLOCATABLE:: xa
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: fa
  EXTERNAL GSOPEN,GSCLOS,PAGES,PAGEE

  CALL GSOPEN

  nxmax=101
  xmin=0.D0
  xmax=10.D0
  
1 WRITE(6,'(A)') '## Input nxmax,xmin,xmax (nxmax=0 for end)'
  WRITE(6,'(A,I6,1P2E12.4)') '## nxmax,xmin,xmax = ',nxmax,xmin,xmax
  READ(5,*,END=9000,ERR=1) nxmax,xmin,xmax
  IF(nxmax.LE.1) GOTO 9000

  ALLOCATE(xa(nxmax),fa(nxmax,4))
  DX=(xmax-xmin)/(nxmax-1)
  nx=1
  xa(nx)=0.D0
  fa(nx,1)=0.D0
  fa(nx,2)=0.D0
  fa(nx,3)=0.D0
  fa(nx,4)=0.D0
  DO nx=2,nxmax
     x=xmin+dx*(nx-1)
     xa(nx)=x
     fa(nx,1)=funca(x)
     fa(nx,2)=funcb(x)
     fa(nx,3)=fa(nx,1)**2+fa(nx,2)**2
     write(6,'(1P4E12.4)') xa(nx),fa(nx,1),fa(nx,2),fa(nx,3)
  END DO

  CALL PAGES
  CALL GRD1D(0,xa,fa,nxmax,nxmax,2,'@hc,hs@',0)
  CALL PAGEE
  CALL PAGES
  CALL GRD1D(0,xa,fa(1:nxmax,3),nxmax,nxmax,1,'@h2@',0)
  CALL PAGEE

  ymin=-2.D0
  ymax= 4.D0
  dy=(ymax-ymin)/(nxmax-1)
  DO nx=1,nxmax
     y=ymin+dy*(nx-1)
     xa(nx)=y
     x=10.D0**y
     fa(nx,4)=funcf(x)
     write(6,'(1P3E12.4)') xa(nx),x,fa(nx,4)
  END DO

  CALL PAGES
  CALL GRD1D(0,xa,fa(1:nxmax,4),nxmax,nxmax,1,'@ff@',1)
  CALL PAGEE
  

  DEALLOCATE(xa,fa)

  GOTO 1

9000 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM plotwi
