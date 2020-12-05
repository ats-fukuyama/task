MODULE libdeab

! XN=XHI*(X-XM):   -1.D0<XN<1.D0
! X =XH*XN+XM:     A < X < B

  USE libde

  REAL(8):: XM_  ! mid XM=0.5D0*(A+B)
  REAL(8):: XH_  ! half len XH=0.5D0*(B-A)
  REAL(8):: XHI_ ! inverse of half len XHI=1.D0/XH

CONTAINS

  SUBROUTINE DEFTAB(A,B,FUNC,KID,H0,EPS,ILST,CS,ES)
    IMPLICIT NONE
    REAL(8),INTENT(IN)::  A           ! Lower boundary of integral
    REAL(8),INTENT(IN)::  B           ! Upper boundary of integral
    CHARACTER(LEN=*),INTENT(IN):: KID ! function identifier string
    REAL(8),INTENT(IN)::  H0          ! Initial step size
    REAL(8),INTENT(IN)::  EPS         ! Convergence thrshold
    INTEGER,INTENT(IN)::  ILST        ! print out control: 0 for no print out
    REAL(8),INTENT(OUT):: CS          ! Integral
    REAL(8),INTENT(OUT):: ES          ! Estimated error 
    INTERFACE
       FUNCTION FUNC(X,XM,XP)
         REAL(8):: FUNC
         REAL(8),INTENT(IN):: X,XM,XP
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
  REAL(8):: gamma_,beta_,alpha_

CONTAINS

  FUNCTION FUNCAB1(X,XA,BX)
    USE libdeab
    IMPLICIT NONE
    REAL(8),INTENT(IN)::  X  ! Integral variable
    REAL(8),INTENT(IN)::  XA ! X-A  near X=A
    REAL(8),INTENT(IN)::  BX ! B-X  near X=B
    REAL(8):: FUNCAB1,DUMMY

    DUMMY=XA
    DUMMY=BX
    FUNCAB1=SIN(alpha_*gamma_*X*(1.D0-X/GAMMA_))*COS(X)

    RETURN
  END FUNCTION FUNCAB1
    
  FUNCTION FUNC1(X,XM,XP)
    USE libdeab
    IMPLICIT NONE
    REAL(8),INTENT(IN)::  X  ! Integral variable
    REAL(8),INTENT(IN)::  XM ! 1.D0-X  near x= 1.D0
    REAL(8),INTENT(IN)::  XP ! 1.D0+X  near X=-1.D0
    REAL(8):: FUNC1

    FUNC1=FUNCAB1(XH_*X+XM_,XH_*XP,XH_*XM)

    RETURN
  END FUNCTION FUNC1

  FUNCTION FUNCAB2(X,XA,BX)
    USE libdeab
    IMPLICIT NONE
    REAL(8),INTENT(IN)::  X  ! Integral variable
    REAL(8),INTENT(IN)::  XA ! X-A  near X=A
    REAL(8),INTENT(IN)::  BX ! B-X  near X=B
    REAL(8):: FUNCAB2,DUMMY

    DUMMY=XA
    DUMMY=BX
    FUNCAB2=SIN(alpha_*gamma_*X*(1.D0-X/GAMMA_))*SIN(X)

    RETURN
  END FUNCTION FUNCAB2
    
  FUNCTION FUNC2(X,XM,XP)
    USE libdeab
    IMPLICIT NONE
    REAL(8),INTENT(IN)::  X  ! Integral variable
    REAL(8),INTENT(IN)::  XM ! 1.D0-X  near x= 1.D0
    REAL(8),INTENT(IN)::  XP ! 1.D0+X  near X=-1.D0
    REAL(8):: FUNC2

    FUNC2=FUNCAB2(XH_*X+XM_,XH_*XP,XH_*XM)

    RETURN
  END FUNCTION FUNC2

  FUNCTION FUNCA(X)
    IMPLICIT NONE
    REAL(8),INTENT(IN):: X
    REAL(8):: FUNCA,H0,EPS,CS,ES
    INTEGER:: ILST

    H0=0.25D0
    EPS=1.D-8
    ILST=0
    gamma_=X
    CALL DEFTAB(0.D0,gamma_,FUNC1,'FUNC1',H0,EPS,ILST,CS,ES)
    FUNCA=CS
  END FUNCTION FUNCA

  FUNCTION FUNCB(X)
    IMPLICIT NONE
    REAL(8),INTENT(IN):: X
    REAL(8):: FUNCB,H0,EPS,CS,ES
    INTEGER:: ILST

    H0=0.25D0
    EPS=1.D-8
    ILST=0
    gamma_=X
    CALL DEFTAB(0.D0,gamma_,FUNC2,'FUNC2',H0,EPS,ILST,CS,ES)
    FUNCB=CS
  END FUNCTION FUNCB

  FUNCTION FUNC0(X)
    IMPLICIT NONE
    REAL(8),INTENT(IN)::  X  ! Integral variable
    REAL(8):: FUNC0,gamma

    gamma=(beta_/alpha_)*SQRT(2.D0*X)
    FUNC0=(FUNCA(gamma)**2+FUNCB(gamma))*EXP(-x)

    RETURN
  END FUNCTION FUNC0

  FUNCTION FUNCF(X)
    IMPLICIT NONE
    REAL(8),INTENT(IN):: X
    REAL(8):: FUNCF,H0,EPS,CS,ES
    INTEGER:: ILST

    H0=0.25D0
    EPS=1.D-8
    ILST=0
    alpha_=0.5D0*beta_**2/X
    CALL DEHIFE(CS,ES,H0,EPS,0,FUNC0,'FF')
    FUNCF=CS
  END FUNCTION FUNCF

END MODULE libde_wi
  
    
PROGRAM plotwi

  USE libgrf
  USE libde_wi
  IMPLICIT NONE
  REAL(8):: xmin,xmax,x,dx,ymin,ymax,dy,y,beta,alpha
  INTEGER:: nxmax,nx
  REAL(8),DIMENSION(:),ALLOCATABLE:: xa
  REAL(8),DIMENSION(:,:),ALLOCATABLE:: fa
  EXTERNAL GSOPEN,GSCLOS,PAGES,PAGEE

  CALL GSOPEN

  nxmax=101
  alpha=0.01D0
  beta=0.1D0
  xmin=0.D0
  xmax=10.D0
  
1 WRITE(6,'(A)') '## Input nxmax,alpha,beta,xmin,xmax (nxmax=0 for end)'
  WRITE(6,'(A,I6,1P4E12.4)') '## nxmax,alpha,beta,xmin,xmax = ', &
                                 nxmax,alpha,beta,xmin,xmax
  READ(5,*,END=9000,ERR=1) nxmax,alpha,beta,xmin,xmax
  IF(nxmax.LE.1) GOTO 9000

  ALLOCATE(xa(nxmax),fa(nxmax,4))
  DX=(xmax-xmin)/(nxmax-1)
  nx=1
  xa(nx)=0.D0
  fa(nx,1)=0.D0
  fa(nx,2)=0.D0
  fa(nx,3)=0.D0
  fa(nx,4)=0.D0
  alpha_=alpha
  beta_=beta
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
  ymax= 0.5D0
  alpha_=alpha
  dy=(ymax-ymin)/(nxmax-1)
  DO nx=1,nxmax
     y=ymin+dy*(nx-1)
     xa(nx)=y
     x=DEXP10(y)
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
