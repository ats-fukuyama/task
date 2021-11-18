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
    INTEGER:: nx,ny,i,j,N,NA,ILL,ngr,nt
    REAL(rkind):: dx,dy,factor,cos_ang,sin_ang,phase,xs,ys,xn,yn,t
    COMPLEX(rkind):: AA(3,3)
   
    omega       = 2.d0*PI*freq
    period      = 1.d0/freq
    wave_length = VC/freq
    wave_number = 2.D0*PI/wave_length

    dt = dtfactor*period
    dx = dxfactor*wave_length
    dy = dyfactor*wave_length

    domega     = 1.0d-3*omega
    omegaplus  = omega+domega
    omegaminus = omega-domega

    nxmax=NINT((xnmax-xnmin)/dxfactor)+1
    nymax=NINT((ynmax-ynmin)/dyfactor)+1
    dx=(xnmax-xnmin)/(nxmax-1)*wave_length
    dy=(ynmax-ynmin)/(nymax-1)*wave_length
    WRITE(6,'(A,2ES12.4,I8,ES12.4)') &
         'xnmin,xnmax,nxmax,dx=',xnmin,xnmax,nxmax,dx
    WRITE(6,'(A,2ES12.4,I8,ES12.4)') &
         'ynmin,ynmax,nymax,dy=',ynmin,ynmax,nymax,dy

    CALL wq_allocate
    
    DO nx=1,nxmax
       xg(nx)=xnmin+dx*(nx-1)
       xgn(nx)=xg(nx)/wave_length
    END DO
    DO ny=1,nymax
       yg(ny)=ynmin+dy*(ny-1)
       ygn(ny)=yg(ny)/wave_length
    END DO

    DO ny=1,nymax
       DO nx=1,nxmax
          medium_nx_ny(nx,ny)=0
       END DO
    END DO

    ! calculate matrix CD

    CALL wq_tens(omegaplus, CDplus)
    CALL wq_tens(omegaminus,CDminus)
    CALL wq_tens(omega,     CD)

    ! calculate matrix A

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

    ! calulate A inverse

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
                Ainv(i,j,nx,ny)=AA(i,j)
             END DO
          END DO

       END DO
    END DO

    ! --- clear time variables ---

    t_tot=0.D0
    nt_tot=0
    ngt_max=0
    ngr_max=0

    ! --- clear wave electric field ---

    DO ny=1,nymax
       DO nx=1,nxmax
          EX(nx,ny)=0.D0
          EY(nx,ny)=0.D0
          EZ(nx,ny)=0.D0
       ENDDO
    ENDDO

    ! --- set initial wave pulse ---
    
    SELECT CASE(model_source)
    CASE(10,11)
       cos_ang=cos(PI*source_angle/180.D0)
       sin_ang=sin(PI*source_angle/180.D0)
       DO ny=2,nymax-1
          DO nx=2,nxmax-1
             xn=xgn(nx)
             yn=ygn(ny)
             xs= (xn-source_position_xn)*cos_ang &
                +(yn-source_position_yn)*sin_ang
             ys=-(xn-source_position_xn)*sin_ang &
                +(yn-source_position_yn)*cos_ang
             factor=EXP(-(xs/source_thickness)**2-(ys/source_width)**2)
             phase=2.D0*PI*xs
             SELECT CASE(model_source)
             CASE(10)
                EX(nx,ny)=-factor*EXP(CI*phase)*sin_ang
                EY(nx,ny)= factor*EXP(CI*phase)*cos_ang
             CASE(11)
                EZ(nx,ny)= factor*EXP(CI*phase)
             END SELECT
          END DO
       END DO

       nt=0
       ngr=1
       t_ngr(ngr)=0.D0
       DO ny=1,nymax
          DO nx=1,nxmax
             EX_save(nx,ny,ngr)=EX(nx,ny)
             EY_save(nx,ny,ngr)=EY(nx,ny)
             EZ_save(nx,ny,ngr)=EZ(nx,ny)
          END DO
       END DO
       ngr_max=ngr
       WRITE(6,'(A,I8,ES12.4,I8)') '## nt,t,ngr:',nt,t,ngr
    END SELECT
  END SUBROUTINE wq_prep
END MODULE wqprep
