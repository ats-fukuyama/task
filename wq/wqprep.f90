! wqprep.f90

MODULE wqprep

  PRIVATE
  PUBLIC wq_prep

CONTAINS

  SUBROUTINE wq_prep
    USE wqcomm
    USE wqtens
    USE libinv
    IMPLICIT NONE
    INTEGER:: nx,ny,k,i,j,nt,N,NA,ILL
   
    omega      = 2.d0*PI*FREQ
    period     = 1.d0/FREQ
    wavelength = VC/FREQ

    dt = dtfactor*Period
    dx = dxfactor*WaveLength
    dy = dyfactor*WaveLength
    nu = nufactor*omega

    domega     = 1.0d-3*omega
    omegaplus  = omega+domega
    omegaminus = omega-domega

    CALL wq_tens(omegaplus, CDplus)
    CALL wq_tens(omegaminus,CDminus)
    CALL wq_tens(omega,     CD)

  ! compute A

    DO ny=1,nymax
       DO nx=1,nxmax
          DO j=1,3
             DO i=1,3
                IF(i.EQ.j) THEN
                   A(i,j,nx,ny)= (2.d0*CI*omega)/(VC**2)        &
                              -RMU0*CD(i,j,nx,ny)               &
                              -omega*RMU0*(CDplus(i,j,nx,ny)    &
                              -CDminus(i,j,nx,ny))/(2.d0*domega)
                ELSE
                   A(i,j,nx,ny)=                                &
                              -RMU0*CD(i,j,nx,ny)               &
                              -omega*RMU0*(CDplus(i,j,nx,ny)    &
                              -CDminus(i,j,nx,ny))/(2.d0*domega)
                END IF
             END DO
          END DO
       END DO
    END DO

! compute A inverse

    N  =3
    NA =3
    ILL=0

    DO ny=1,nymax
       DO nx=1,nxmax
          DO j=1,3
             DO i=1,3
                AA(i,j)=A(i,j,nx,ny) ! AA is A(l,m)
             END DO
          END DO

          CALL INVMCD(AA,N,NA,ILL) ! compute AA inverse

          IF(ILL.EQ.900) THEN
             WRITE(*,*) "ERROR: Computing AA inverse is impossible"
             CALL GSCLOS
             STOP
          END IF

          DO j=1,3
             DO i=1,3
                A(i,j,nx,ny)=AA(i,j)
             END DO
          END DO

       END DO
    END DO

    ttot=0.d0
    nttot=0
    ntplot=0

    DO ny=1,nymax
       DO nx=1,nxmax
          EX(nx,ny)=0.D0
          EY(nx,ny)=0.D0
          EZ(nx,ny)=0.D0
       ENDDO
    ENDDO

  END SUBROUTINE wq_prep
END MODULE wqprep
