! fowlib.f90

MODULE fowlib

  PRIVATE
  PUBLIC solve_quadratic_equation
  PUBLIC first_order_derivative
  PUBLIC second_order_derivative
  PUBLIC gauss_jordan
  PUBLIC fow_cal_spl
  PUBLIC fow_cal_spl2D

CONTAINS
  
  subroutine solve_quadratic_equation(z,C)
    ! solve C(3)*z**2+C(2)*z+C(1) = 0
    ! z(1) : -sqrt(D)
    ! z(2) : +sqrt(D)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(in) :: C(3)
    COMPLEX(rkind),intent(out) :: z(2)
    real(rkind) :: D
    COMPLEX(rkind) :: ei=(0.d0,1.d0)

    z(1) = (0.d0, 0.d0)
    z(2) = (0.d0, 0.d0)

    D = C(2)**2-4.d0*C(1)*C(3)

    if ( D >= 0.d0 ) then
      z(1) = (-C(2)-sqrt(D))/(2.d0*C(3))+0.d0*ei
      z(2) = (-C(2)+sqrt(D))/(2.d0*C(3))+0.d0*ei
    else
      z(1) = (-C(2)-ei*sqrt(-D))/(2.d0*C(3))
      z(2) = (-C(2)+ei*sqrt(-D))/(2.d0*C(3))
    end if

  end subroutine solve_quadratic_equation

  subroutine first_order_derivative(dfdx, f, x)
    ! calcurate dfdx, the first order derivative of f
    ! error order is O(dx**2)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: dfdx(:)
    real(rkind),intent(in)  :: f(:), x(:)
    real(rkind) :: s, t
    integer :: imax, i

    imax = size(f)
    if ( imax /= size(x) .or. imax /= size(dfdx) ) then
      write(*,*)"ERROR at subroutine first_order_derivative in fowprep.f90"
    end if

    do i = 1, imax
      if ( i /= imax .and. i /= 1) then
        t = x(i+1)-x(i)
        s = x(i)-x(i-1)
        dfdx(i) = (s**2*f(i+1)+(t**2-s**2)*f(i)-t**2*f(i-1))/(s*t*(s+t))
      else if ( i == 1 ) then
        s = x(2)-x(1)
        t = x(3)-x(2)
        dfdx(1) = ((s**2-(s+t)**2)*f(1)+(s+t)**2*f(2)-s**2*f(3))/(s*t*(s+t))
      else if ( i == imax ) then
        t = x(imax)-x(imax-1)
        s = x(imax-1)-x(imax-2)
        dfdx(imax) = (((s+t)**2-t**2)*f(imax)-(s+t)**2*f(imax-1)+t**2*f(imax-2))/(s*t*(s+t))
      end if
    end do

  end subroutine first_order_derivative

  subroutine second_order_derivative(d2fdx2, f, x)
    ! calcurate d2fdx2, the second order derivative of f
    ! error order is O(dx**2)

    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: d2fdx2(:)
    real(rkind),intent(in)  :: f(:), x(:)
    real(rkind) :: s, t, r, v, w, A(3,3), B(3)
    integer :: imax, i

    imax = size(f)

    do i = 1, imax
      if ( i /= imax .and. i /= 1) then
        t = x(i+1)-x(i)
        s = x(i)-x(i-1)
        d2fdx2(i) = (s*f(i+1)-(s+t)*f(i)+t*f(i-1))/(0.5d0*s*t*(s+t))
      else if ( i == 1 ) then
        r = x(2)-x(1)
        s = x(3)-x(2)
        t = x(4)-x(3)
        v = r+s
        w = r+s+t

        A(1,1) = r
        A(1,2) = r**2/2
        A(1,3) = r**3/6

        A(2,1) = v
        A(2,2) = v**2/2
        A(2,3) = v**3/6

        A(3,1) = w
        A(3,2) = w**2/2
        A(3,3) = w**3/6

        B(1) = f(2)-f(1)
        B(2) = f(3)-f(1)
        B(3) = f(4)-f(1)

        call gauss_jordan(A, B, 3)

        d2fdx2(1) = B(2)

      else if ( i == imax ) then
        r = x(imax-2)-x(imax-3)
        s = x(imax-1)-x(imax-2)
        t = x(imax)-x(imax-1)
        v = s+t
        w = r+s+t

        A(1,1) = -1.d0*t
        A(1,2) = t**2/2
        A(1,3) = -1.d0*t**3/6

        A(2,1) = -1.d0*v
        A(2,2) = v**2/2
        A(2,3) = -1.d0*v**3/6

        A(3,1) = -1.d0*w
        A(3,2) = w**2/2
        A(3,3) = -1.d0*w**3/6

        B(1) = f(imax-1)-f(imax)
        B(2) = f(imax-2)-f(imax)
        B(3) = f(imax-3)-f(imax)

        call gauss_jordan(A, B, 3)

        d2fdx2(imax) = B(2)
        
      end if
    end do

  end subroutine second_order_derivative

  subroutine gauss_jordan(A, B, n)
    use fpcomm,only:rkind

    implicit none
    real(rkind) :: A(n,n), B(n)
    integer n,i,j,k
  
    do k = 1, n
      do j = k + 1, n
        a(k,j) = a(k,j) / a(k,k)
      end do
      b(k) = b(k) / a(k,k)
  
      do i = 1, n
        if ( i .ne. k ) then
          do j = k + 1, n
            a(i,j) = a(i,j) - a(i,k) * a(k,j)
          end do
          b(i) = b(i) - a(i,k) * b(k)
        end if
      end do
  
    end do
  
  end subroutine gauss_jordan
  
  subroutine fow_cal_spl(f_out, x_in, f, x)
    ! Calculate spline coefficient y = f(x),
    ! then return f_out = f(x_in)
    use fpcomm,only:rkind
    USE libspl1d
    implicit none
    real(rkind),intent(out) :: f_out
    real(rkind),intent(in) :: x_in, f(:), x(:)
    integer :: imax,ierr
    real(rkind),allocatable :: U(:,:), fx(:)

    ierr=0
    imax = size(f,1)

    allocate(U(4,imax), fx(imax))

    call first_order_derivative(fx, f, x)

    call SPL1D(x,f,fx,U,imax,3,IERR)

    call SPL1DF(x_in,f_out,x,U,imax,IERR)

    return

  end subroutine fow_cal_spl

  subroutine fow_cal_spl2D(f_out, x_in, y_in, f, x, y)
    ! 2D-version of fow_cal_spl
    use fpcomm,only:rkind
    USE libspl2d
    implicit none
    real(rkind),intent(out) :: f_out
    real(rkind),intent(in) :: x_in, y_in, f(:,:), x(:), y(:)
    integer :: imax, jmax, ierr
    real(rkind),allocatable :: U(:,:,:,:), fx(:,:), fy(:,:) , fxy(:,:)

    ierr=0
    imax = size(x)
    jmax = size(y)

    allocate(U(4,4,imax,jmax),fx(imax,jmax),fy(imax,jmax),fxy(imax,jmax))

    if ( size(f,1) /= imax &
        .or. size(f,2) /= jmax ) then
      write(*,*)"error at fow_cal_spl2D"
      STOP
    end if

    call SPL2D(x,y,f,fx,fy,fxy,U,imax,imax,jmax,0,0,IERR)

    call SPL2DF(x_in,y_in,f_out,x,y,U,imax,imax,jmax,ierr)

  end subroutine fow_cal_spl2D

END MODULE fowlib
