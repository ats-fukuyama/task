PROGRAM disp

! cyclotron resonance:     blue  
! 2nd cyclotron harmonics: blue dot-dash
! upper hybrid resonance:  red
! ordinary cutoff:         green  dash
! right hand cutoff:       yellow dash
! left hand cutoff:        light blue dash

  USE libgrf
  IMPLICIT NONE
  INTEGER:: nmax,n
  REAL(8):: Xmin,Xmax,Ymin,Ymax,Npara,dz,delX,delY,D
  REAL(8),DIMENSION(:),ALLOCATABLE:: X,Y,Z,WC,WCC,WUH,OC,RC1,RC2,LC1,LC2
  REAL(8),DIMENSION(:,:),ALLOCATABLE:: G
  REAL(8),DIMENSION(:,:),ALLOCATABLE:: LINE_RGB
  INTEGER,DIMENSION(:),ALLOCATABLE:: LINE_MARK,LINE_PAT

  CALL gsopen
  Xmin=1.0D0   ! X = omegap2^2/omegape0^2
  Xmax=1.0D0
  Ymin=1.0D0   ! Y = -omegace/omegape0
  Ymax=0.5D0
  Npara=0.5D0
  nmax=101

  ALLOCATE(LINE_MARK(8),LINE_PAT(8),LINE_RGB(3,8))
  DO n=1,8
     LINE_MARK(n)=0
     LINE_PAT(n)=0
  END DO
  LINE_PAT(2)=4
  LINE_PAT(4)=2
  LINE_PAT(5)=2
  LINE_PAT(6)=2
  LINE_PAT(7)=2
  LINE_PAT(8)=2
  LINE_RGB(1:3, 1)=(/0.0D0,0.0D0,1.0D0/)  ! WC
  LINE_RGB(1:3, 2)=(/0.0D0,0.0D0,1.0D0/)  ! 2WC
  LINE_RGB(1:3, 3)=(/1.0D0,0.0D0,0.0D0/)  ! WUH
  LINE_RGB(1:3, 4)=(/0.0D0,1.0D0,0.0D0/)  ! OC
  LINE_RGB(1:3, 5)=(/1.0D0,1.0D0,0.0D0/)  ! RC1
  LINE_RGB(1:3, 6)=(/1.0D0,1.0D0,0.0D0/)  ! RC2
  LINE_RGB(1:3, 7)=(/0.0D0,1.0D0,1.0D0/)  ! LC1
  LINE_RGB(1:3, 8)=(/0.0D0,1.0D0,1.0D0/)  ! LC2


1 WRITE(6,'(A          )') '## DISP: Xmin,Xmax,Ymin,Ymax,Npara,nmax'
  WRITE(6,'(1P5E12.4,I5)') Xmin,Xmax,Ymin,Ymax,Npara,nmax
  WRITE(6,'(A          )') 'Input: Xmin,Xmax,Ymin,Ymax,Npara,nmax'
  READ(5,*,END=9,ERR=1) Xmin,Xmax,Ymin,Ymax,Npara,nmax

  ALLOCATE(X(nmax),Y(nmax),Z(nmax),WC(nmax),WCC(nmax),WUH(nmax))
  ALLOCATE(OC(nmax),RC1(nmax),RC2(nmax),LC1(nmax),LC2(nmax))
  ALLOCATE(G(nmax,8))

  dz=1.d0/(nmax-1)
  delX=(Xmax-Xmin)/(nmax-1)
  delY=(Ymax-Ymin)/(nmax-1)
  DO n=1,nmax
     Z(n)=dz*(n-1)
     X(n)=Xmin+delX*(n-1)
     Y(n)=Ymin+delY*(n-1)
     WC(n)=Y(n)
     WCC(n)=2.D0*Y(n)
     WUH(n)=SQRT(X(n)**2+Y(n)**2)
     OC(n)=X(n)
     D=Y(n)**2+4.D0*X(n)/(1-Npara**2)

     IF(D.GE.0.D0) THEN
        RC1(n)=0.5D0*( Y(n)+DSQRT(D))
        RC2(n)=0.5D0*( Y(n)-DSQRT(D))
        LC1(n)=0.5D0*(-Y(n)+DSQRT(D))
        LC2(n)=0.5D0*(-Y(n)-DSQRT(D))
     ELSE
        RC1(n)=0.0D0
        RC2(n)=0.0D0
        LC1(n)=0.0D0
        LC2(n)=0.0D0
     END IF
  END DO

  DO n=1,nmax
     G(n,1)=WC(n)
     G(n,2)=WCC(n)
     G(n,3)=WUH(n)
     G(n,4)=OC(n)
     G(n,5)=RC1(n)
     G(n,6)=RC2(n)
     G(n,7)=LC1(n)
     G(n,8)=LC2(n)
  END DO

  CALL pages
  CALL GRD1D(0,Z,G,nmax,nmax,8,'@WC,2WC,WUH,OC,RC1/2,LC1/2@', &
             FMIN=0.D0,NLMAX=8, &
             LINE_MARK=LINE_MARK,LINE_RGB=LINE_RGB,LINE_PAT=LINE_PAT)
  CALL pagee
  DEALLOCATE(X,Y,Z,WC,WCC,WUH,OC,RC1,RC2,LC1,LC2,G)
  GO TO 1

9 CONTINUE
  CALL gsclos
  STOP
END PROGRAM disp

  
