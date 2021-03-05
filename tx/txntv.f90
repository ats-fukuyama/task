!***************************************************************
!
!   Evaluate perturbed magntetic field
!
!      Output : deltam(0:nrmax,0:m_pol-1)
!
!***************************************************************

subroutine perturb_mag
  use tx_commons, only : nrmax, m_pol, PI, n_tor, deltam, r
  use tx_ripple, only : ripple
  use libfft,only : ddct
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
        delta_l(m) = ripple(r(nr),theta(m),1.d0)
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

subroutine Wnm_spline(fmnq, wnm, umnq, nmnqm)
  USE libspl1d
  implicit none
  integer(4), intent(in) :: nmnqm
  real(8), dimension(1:nmnqm), intent(out) :: fmnq, wnm
  real(8), dimension(1:4,1:nmnqm), intent(out) :: umnq
  integer(4) :: ist, nmnq, ierr
  real(8), dimension(:), allocatable :: deriv

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

end subroutine Wnm_spline

!***************************************************************
!
!   Calculate neoclassical toroidal viscosity
!
!      Output : rNuNTV(0:NRMAX), UastNC(0:NRMAX)
!
!***************************************************************

subroutine NTVcalc
  use tx_commons, only : rNuNTV, UastNC, m_pol, n_tor, q, fmnq, umnq, nrmax, nmnqm, &
       &                 DltRP, R, RR, PTiV, rKeV, AMI, rNuii, PTiV, RA, PZ, AEE, BB, &
       &                 deltam, m_pol
  USE libspl1d
  implicit none

  integer(4) :: nr, m, ierr
  real(8) :: suml, f, wnm_plus, wnm_minus, EpsL, Vti2, dPTiV, B_lambda
  real(8) :: deriv4, aitken2p ! external function

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
!     write(6,*) r(nr)/ra,sum(deltam(nr,0:m_pol-1)**2)
     B_lambda = 0.5d0 * n_tor**2 * suml
!     write(6,*) r(nr)/ra,B_lambda

     EpsL = R(nr) / RR
     Vti2 = ABS(PTiV(nr)) * rKeV / AMI
     rNuNTV(nr) = 1.74d0 * Vti2 / rNuii(nr) * EpsL**1.5d0 / RR**2 * B_lambda
!     write(6,*) r(nr)/ra,rnuntv(nr)

     dPTiV = DERIV4(NR,R,PTiV,NRMAX,0) * RA * rKeV
     if(nr /= 0) UastNC(nr) = 3.5d0 * RR * Q(nr) / (PZ * AEE * R(nr) * BB) * dPTiV
  end do
  UastNC(0) = AITKEN2P(R(0),UastNC(1),UastNC(2),UastNC(3),R(1),R(2),R(3))

end subroutine NTVcalc
