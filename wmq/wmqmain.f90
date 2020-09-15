! 
!   wmq
!         ---ANALYSIS OF WAVE PROPAGATION USING QUASI-WAVE OPTICS
!
!                               coded by Y.Maruyama

program wmqmain

  use bpsd
  implicit none
  integer(4) :: W,INMODE,flag
  integer(4) :: i,j,k,l,m,nt,ntmax,N,NA,ILL,steps
  real(8)    :: dt,dx,nu,dtfactor,dxfactor,nufactor,B0,RR,RA,q0,qa,n0,TMN
  real(8)    :: t,FREQ,omega,period,wavelength,omegaplus,omegaminus,domega,APS
  real(8)   ,ALLOCATABLE :: ne(:,:),OCE(:,:),OPE(:,:),OUH(:,:),AP(:,:),OR(:,:),OL(:,:)
  complex(8),ALLOCATABLE :: EX(:,:),EY(:,:),EZ(:,:),A(:,:,:,:),AA(:,:),B(:,:)
  complex(8),ALLOCATABLE :: CD(:,:,:,:),CDplus(:,:,:,:),CDminus(:,:,:,:)
   
  ! Initialize

  OPEN(10,FILE='input.dat')
  NAMELIST /input/ FREQ,dtfactor,dxfactor,nufactor,B0,RR,RA,q0,qa,n0,W,INMODE,TMN
  READ(10,input)

  OPEN(5000,FILE='APS.d')

  omega      = 2.d0*PI*FREQ
  period     = 1.d0/FREQ
  wavelength = VC/FREQ

  dt = dtfactor*Period
  dx = dxfactor*WaveLength
  nu = nufactor*omega

  domega     = 1.0d-3*omega
  omegaplus  = omega+domega
  omegaminus = omega-domega

  ALLOCATE(EX(W,W),EY(W,W),EZ(W,W))
  ALLOCATE(A(W,W,3,3),AA(3,3),B(3,3))
  ALLOCATE(CD(W,W,3,3),CDplus(W,W,3,3),CDminus(W,W,3,3))
  ALLOCATE(ne(W,W),OCE(W,W),OPE(W,W),OUH(W,W),AP(W,W),OR(W,W),OL(W,W))

  CALL GSOPEN
  
  CALL wmqtens(CDplus ,ne,OCE,OPE,OUH,OR,OL,W,dx,omegaplus ,B0,nu,RR,RA,q0,qa,n0)
  CALL wmqtens(CDminus,ne,OCE,OPE,OUH,OR,OL,W,dx,omegaminus,B0,nu,RR,RA,q0,qa,n0)
  CALL wmqtens(CD     ,ne,OCE,OPE,OUH,OR,OL,W,dx,omega     ,B0,nu,RR,RA,q0,qa,n0)

  ! compute A

  DO j=1,W
     DO i=1,W
        DO m=1,3
           DO l=1,3
              IF(l.EQ.m) THEN
                 A(i,j,l,m)= (2.d0*CI*omega)/(VC**2)        &
                           &-RMU0*CD(i,j,l,m)               &
                           &-omega*RMU0*(CDplus(i,j,l,m)    &
                           &-CDminus(i,j,l,m))/(2.d0*domega)
              ELSE
                 A(i,j,l,m)=                                &
                           &-RMU0*CD(i,j,l,m)               &
                           &-omega*RMU0*(CDplus(i,j,l,m)    &
                           &-CDminus(i,j,l,m))/(2.d0*domega)
              END IF
           END DO
        END DO
     END DO
  END DO

! compute A inverse

  N  =3
  NA =3
  ILL=0

  DO j=1,W
     DO i=1,W

        DO m=1,3
           DO l=1,3
              AA(l,m)=A(i,j,l,m) ! AA is A(l,m)
           END DO
        END DO

        CALL INVMCD(AA,N,NA,ILL) ! compute AA inverse

        IF(ILL.EQ.900) THEN
           WRITE(*,*) "ERROR: Computing AA inverse is impossible"
           CALL GSCLOS
           STOP
        END IF

        DO m=1,3
           DO l=1,3
              A(i,j,l,m)=AA(l,m)
           END DO
        END DO

     END DO
  END DO

  t=0.d0

  DO j=1,W
     DO i=1,W
        EX(i,j)=0.D0
        EY(i,j)=0.D0
        EZ(i,j)=0.D0
     ENDDO
  ENDDO

  steps = 0
  flag = 1

  1 CONTINUE

  write(*,'(A10,1PE12.4,1X,A7)')  "delta_t is", dtfactor, "*Period"
  write(*,'(A10,1PE12.4,1X,A11)') "delta_x is", dxfactor, "*WaveLength"
  write(*,*)                      "----Please input ntmax----"
  read (*,*) ntmax
  if(ntmax.eq.0) goto 9000

  do nt=1,ntmax
     ! give Boundary Condition
     i=W
     do j=1,W
        if(INMODE.eq.1) then
           EY(i,j) = 1.0d0/sqrt(2.d0)*exp(-5.0d-4/wavelength*(((j-W/2.D0))**2))
           EZ(i,j) = 1.0d0/sqrt(2.d0)*exp(-5.0d-4/wavelength*(((j-W/2.D0))**2))
        else if(INMODE.eq.2) then
           EY(i,j) = 1.0d0*exp(-5.0d-4/wavelength*(((j-W/2.D0))**2))
        else if(INMODE.eq.3) then
           EZ(i,j) = 1.0d0*exp(-5.0d-4/wavelength*(((j-W/2.D0))**2))
        else
           write(*,*) "ERROR:INMODE needs to be 1or2or3"
        end if
     end do

     !compute next E
     call wmqsolv(EX,EY,EZ,W,A,CD,RR,dt,dx,omega,TMN)
     t=t+dt
     if(mod(nt,1000).eq.0)then
        write(6,'(I4)') (ntmax-nt)/1000

        ! compute Absorbed Power
        do j=2,W-1
           do i=2,W-1
              AP(i,j)=0.5d0*real( conjg(EX(i,j))*(CD(i,j,1,1)*EX(i,j) &
                                                 +CD(i,j,1,2)*EY(i,j) &
                                                 +CD(i,j,1,3)*EZ(i,j))&
                                 +conjg(EY(i,j))*(CD(i,j,2,1)*EX(i,j) &
                                                 +CD(i,j,2,2)*EY(i,j) &
                                                 +CD(i,j,2,3)*EZ(i,j))&
                                 +conjg(EZ(i,j))*(CD(i,j,3,1)*EX(i,j) &
                                                 +CD(i,j,3,2)*EY(i,j) &
                                                 +CD(i,j,3,3)*EZ(i,j)))
           end do
        end do
        APS=0.d0
        do j=1,W
           do i=1,W
              APS=APS+AP(i,j)
           end do
        end do
        write(5000,*) t, APS
     end if
  end do
  
  steps=steps+ntmax

  call wmqgout(W,dx,RR,omega,EX,EY,EZ,ne,OCE,OPE,OUH,AP,OR,OL,flag)  

  go to 1

9000  call GSCLOS

  write(*,*) "compute", steps, "step"

  close(5000)
  stop
end program wmqmain
