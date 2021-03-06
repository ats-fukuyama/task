module tx_ntv
  implicit none

  integer(4), parameter :: NMNQM=446, M_POL_M=64
  ! Read Wnm table
  real(8), dimension(:),   allocatable :: Fmnq, Wnm, rNuNTV, UastNC
  real(8), dimension(:,:), allocatable :: Umnq, deltam

contains

  !***************************************************************
  !
  !   Allocate variables
  !
  !***************************************************************

  subroutine alloc_ntv
    use tx_commons, only : NRMAX
    implicit none
    integer(4) :: ier

    allocate(Fmnq(1:NMNQM), Wnm(1:NMNQM), Umnq(1:4,1:NMNQM), deltam(0:NRMAX,0:M_POL_M), &
         &   rNuNTV(0:NRMAX), UastNC(0:NRMAX), stat=ier)
    if(ier /= 0) stop 'alloc_ntv error.'

  end subroutine alloc_ntv

  !***************************************************************
  !
  !   Deallocate variables
  !
  !***************************************************************

  subroutine dealloc_ntv
    implicit none
    integer(4) :: ier

    if(allocated(Fmnq)) then
       deallocate(Fmnq, Wnm, Umnq, deltam, rNuNTV, UastNC, stat=ier)
       if(ier /= 0) stop 'dealloc_ntv error.'
    end if

  end subroutine dealloc_ntv

  !***************************************************************
  !
  !   Evaluate perturbed magntetic field
  !
  !      Output : deltam(0:nrmax,0:m_pol-1)
  !
  !***************************************************************

  subroutine perturb_mag
    use libfft, only : ddct
    use tx_commons, only : nrmax, m_pol, PI, rho
    use tx_ripple, only : ripple
    implicit none

    integer(4) :: m, m_max, max_ip, max_w, nr
    integer(4), dimension(:), allocatable :: ip
    real(8), dimension(:), allocatable :: theta, delta_l, w
    real(8) :: dtheta

    m_max = m_pol

    dtheta = 2.D0 * PI / m_max
    allocate(theta(0:m_max-1))
    do m = 0, m_max-1
       theta(m) = m * dtheta
    end do

    max_ip = 2+sqrt(0.5*m_max)
    max_w  = m_max*1.25-1
    allocate(delta_l(0:m_max),ip(0:max_ip),w(0:max_w))

    do nr = 0, nrmax
       ! ripple amplitude in the poloidal plane at a certain radial position
       do m = 0, m_max-1
          delta_l(m) = ripple(rho(nr),theta(m),1.d0)
       end do
       ! Discrete Cosine Transformation
       call ddct(m_max, -1, delta_l, ip, w)
       ! Scaling (summping up delta_l for all m makes original delta_l, i.e. ripple amp.)
       deltam(nr,0:m_max-1) = 2.d0/m_max * delta_l(0:m_max-1)
       deltam(nr,0) = 0.5d0 * delta_l(0)
    end do

    deallocate(theta,delta_l,ip,w)

  end subroutine perturb_mag

  !***************************************************************
  !
  !   Read the table of pitch angle integral from m-nq=0 to 300
  !   from "wnm_data.d" and make the spline matrix.
  !
  !      Output : fmnq(1:nmnqm), wnm(1:nmnqm), umnq(1:4,1:nmnqm)
  !
  !***************************************************************

  subroutine Wnm_spline
    USE libspl1d
    implicit none
    integer(4) :: ist, nmnq, ierr
    real(8), dimension(:), allocatable :: deriv

    ! allocate variables if not allocated
    if(.not. allocated(Fmnq)) call alloc_ntv

    if(maxval(fmnq) == 0.D0) then

       open(21,file='wnm_data.d',iostat=ist,status='old',form='formatted')
       if(ist /= 0) stop 'Wnm_spline: cannot read "wnm_data.d"'
       do nmnq = 1, nmnqm
          read(21,*) fmnq(nmnq), wnm(nmnq)
       end do
       close(21)

       allocate(deriv(1:nmnqm))
       call spl1d(fmnq,wnm,deriv,umnq,nmnqm,0,ierr)
       if(ierr /= 0) stop 'Wnm_spline: cannot make spline matrix'
       deallocate(deriv)
    end if

  end subroutine Wnm_spline

  !***************************************************************
  !
  !   Calculate neoclassical toroidal viscosity
  !
  !      Output : rNuNTV(0:NRMAX), UastNC(0:NRMAX)
  !
  !***************************************************************

  subroutine NTVcalc
    use tx_commons, only : m_pol, n_tor, q, nrmax, rho, r, &
         &                 epst, RR, rKeV, amas, amp, rNuii, achg, AEE, BB, Var
    USE libspl1d
    USE libitp
    implicit none

    integer(4) :: nr, m, ierr
    real(8) :: suml, f, wnm_plus, wnm_minus, Vti2, dPTiV, B_lambda

    do nr = 0, nrmax
       suml = 0.d0
       do m = 0, m_pol - 1
          f = m - n_tor * q(nr)
          if(abs(f) > 300.d0) stop 'm-nq goes outside the Wnm table given in "wnm_table.d"'
          call spl1df(abs(f),wnm_plus,fmnq,umnq,nmnqm,ierr)
          if(ierr /= 0) stop 'Wnm_spline_1: cannot make spline matrix'

          f = m + n_tor * q(nr)
          if(abs(f) > 300.d0) stop 'm-nq goes outside the Wnm table given in "wnm_table.d"'
          call spl1df(abs(f),wnm_minus,fmnq,umnq,nmnqm,ierr)
          if(ierr /= 0) stop 'Wnm_spline_2: cannot make spline matrix'

          suml = suml + deltam(nr,m)**2 * (wnm_plus + wnm_minus)
       end do
       !     write(6,*) rho(nr),sum(deltam(nr,0:m_pol-1)**2)
       B_lambda = 0.5d0 * n_tor**2 * suml
       !     write(6,*) rho(nr),B_lambda

       Vti2 = ABS(Var(nr,2)%T) * rKeV / (amas(2)*amp)
       rNuNTV(nr) = 1.74d0 * Vti2 / rNuii(nr) * epst(NR)**1.5d0 / RR**2 * B_lambda
       !     write(6,*) rho(nr),rnuntv(nr)

       dPTiV = DERIV4(NR,rho,Var(0:NRMAX,2)%T,NRMAX,0) * rKeV
       if(nr /= 0) UastNC(nr) = 3.5d0 * RR * Q(nr) / (achg(2) * AEE * R(nr) * BB) * dPTiV
    end do
    UastNC(0) = AITKEN2P(R(0),UastNC(1),UastNC(2),UastNC(3),R(1),R(2),R(3))

  end subroutine NTVcalc

end module tx_ntv
