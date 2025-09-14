module var_equil3D
  use bpsd, only: rkind,dp 
  implicit none

!... wout file name
  character(len=80)             :: wout_file_name ! name of wout file from vmec (input)

!... data from wout_file 
  character(len=4)              :: version        ! version of vmec2000 ; accept '6.90' or '8.40' at present
  integer           :: ns             ! the number of radial grids ; j=0:ns
  integer           :: nsp1           ! = ns + 1
  integer           :: nsm1           ! = ns - 1
  integer           :: nfp            ! the number of toroidal field periods
  integer           :: mpol           ! the number of poloidal modes ;  0 <= m <= mpol-1
  integer           :: ntor           ! the number of toroidal modes ;  -ntor <= n <= ntor
  integer           :: mnmax          ! the number of Fourier modes for R, Z, and lambda
  integer           :: mnmax2         ! the number of Fourier modes for B and gsqrt ( =mnmax for version 6.9)
  integer           :: isgn           ! sign of  gsqrth in vmec2000 (left-hand coordinates -> isgn=-1)
  real(kind=rkind), allocatable :: xm(:)          ! poloidal mode number for R, Z, and lambda ; j=1:mnmax
  real(kind=rkind), allocatable :: xn(:)          ! toroidal mode number for R, Z, and lambda ; j=1:mnmax
  real(kind=rkind), allocatable :: xm2(:)         ! poloidal mode number for B and gsqrt      ; j=1:mnmax2
  real(kind=rkind), allocatable :: xn2(:)         ! toroidal mode number for B and gsqrt      ; j=1:mnmax2
  real(kind=rkind), allocatable :: rmnc(:,:)      ! Rmn on full grids        ; cos,  (1:mnmax,0:ns)
  real(kind=rkind), allocatable :: zmns(:,:)      ! Zmn on full grids        ; sin,  (1:mnmax,0:ns)
  real(kind=rkind), allocatable :: lmnsh(:,:)     ! lambda_mn on half grids  ; sin,  (1:mnmax,0:ns)
  real(kind=rkind), allocatable :: bmnh(:,:)      ! Bmn on half grids        ; cos,  (1:mnmax2,0:ns)
  real(kind=rkind), allocatable :: gmnh(:,:)      ! gsqrt_mn on half grids   ; cos,  (1:mnmax2,0:ns)
  real(kind=rkind), allocatable :: aiotah(:)      ! rotational transform on half grids ; j=0,ns+1
  real(kind=rkind), allocatable :: presh(:)       ! pressure on half grids [Pa]        ; j=0,ns+1
  real(kind=rkind), allocatable :: Itorh(:)       ! twopi*isgn*Itorh/mu0 = net toroidal current[A] on half grids ; j=0,ns+1
  real(kind=rkind), allocatable :: Ipolh(:)       ! twopi*isgn*Ipolh/mu0 = net poloidal current[A] on half grids ; j=0,ns+1
  real(kind=rkind), allocatable :: vprimh(:)      ! (dV/ds)/twopi**2 on half grids, where V is volume [m3] ; j=0,ns+1
  real(kind=rkind)              :: phiedge        ! (total toroidal flux at the edge)/(isgn*twopi) 

!... equilibrium quantities from vmec
  real(kind=rkind), allocatable :: iota_vmec(:)   ! rotational transform from vmec2000 ; j=0,ns
  real(kind=rkind), allocatable :: pprim(:)       ! dp/ds from vmec2000, where p[Pa] ; j=0,ns
  real(kind=rkind), allocatable :: S11(:)         ! S11 on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: S12(:)         ! S12 on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: Bsqav(:)       ! <B^2> [T2] on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: vprim(:)       ! dV/ds [m3] on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: grdssq_av(:)   ! <|grad(s)|^2> on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: Itorf(:)       ! net toroidal current [A] on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: Ipolf(:)       ! net poloidal current [A] on a full mesh ; j=0,ns
  real(kind=rkind), allocatable :: s(:)           ! radial grid points (full mesh) = s ; s is normalized tor. flux ; j=0:ns

end module var_equil3D
