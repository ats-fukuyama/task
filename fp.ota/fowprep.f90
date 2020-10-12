module fowprep
  PRIVATE
  PUBLIC :: fow_prep
contains

  subroutine fow_prep

    use fowcomm
    use fpcomm,only:rkind,nrmax,nthmax,npmax

    implicit none

    integer :: nth,nr,nrg,nzg,ierr = 0,eps = 1.0d-8
    real(rkind),allocatable :: fx(:,:),fy(:,:),fxy(:,:),R(:),Z(:)

    do nth = 1,nthmax
      xi(nth) = (nth-0.5d0)/nthmax*(-2.d0)+1.d0
      xig(nth) = (nth-1.d0)/nthmax*(-2.d0)+1.d0
    end do
    xig(nthmax+1) = -1.d0

    ! do nr = 1,nrmax
    !   do nth = 1,nthmax
    !     if(xi(nth)<0.d0)then
    !       BBm(nth,nr) = Bin(nr)
    !     else
    !       BBm(nth,nr) = Bout(nr)
    !     end if
    !   end do
    ! end do

    ! allocate(fx(nrgmax,nzgmax),fy(nrgmax,nzgmax),fxy(nrgmax,nzgmax),R(nrgmax),Z(nzgmax))

    ! do nrg = 1,nrgmax
    !   R(nrg) = nrg*1.d0
    ! end do
    ! do nzg = 1,nzgmax
    !   Z(nzg) = nzg*1.d0
    ! end do

    ! call SPL2D(R,Z,Brz,fx,fy,fxy,UBspl,nrgmax,nrgmax,nzgmax,0,0,ierr)
    ! call SPL2D(R,Z,Frz,fx,fy,fxy,UFspl,nrgmax,nrgmax,nzgmax,0,0,ierr)

    call fow_eqload(ierr)

    ! deallocate(fx,fy,fxy,R,Z)

  end subroutine fow_prep

  subroutine fow_eqload(ierr)
    use fowcomm
    use fpcomm,only:rkind,nrmax,nthmax,npmax
    USE libgrf
    implicit none
    integer,intent(out):: ierr
    character(len = 80) :: line
    real(rkind) :: rr_axis,zz_axis,psit0,qaxis,qsurf
    real(rkind),allocatable,dimension(:) :: ppsi,qpsi,vpsi,rlen,ritpsi,rhot
    real(rkind),allocatable,dimension(:,:) :: Br,Bz,Bp,Bt
    integer :: nr,nthp

    allocate(ppsi(nrmax+1),qpsi(nrmax+1),vpsi(nrmax+1),rlen(nrmax+1),ritpsi(nrmax+1)&
            ,rhot(nrmax+1))
    allocate(Br(nthpmax,nrmax+1),Bz(nthpmax,nrmax+1),Bp(nthpmax,nrmax+1),Bt(nthpmax,nrmax+1))

    ierr = 0
    write(*,*)"FP------------------------------------"

    modelg=3
    call eqload(modelg,knameq,ierr)
    if(ierr.ne.0) return

    write(line,'(a)') 'mdleqc = 1'   ! set boozer poloidal angle
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nrmax = ',nrmax+1
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nthmax = ',nthpmax
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nsumax = ',nrmax+1
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    call eqcalq(ierr)
    if(ierr.ne.0) return

    call eqgetp(rhot,psimg,nrmax+1)                        ! normalized psit radius , use only psimg
    call eqgetqn(ppsi,qpsi,Fpsig,vpsi,rlen,ritpsi,nrmax+1) ! flux functions         , use only Fpsig
    call eqgetbb(Br,Bz,Bp,Bt,nthpmax,nthpmax,nrmax+1)      ! mag field              , use only Bt and Bp
    call eqgeta(rr_axis,zz_axis,psi0,psit0,qaxis,qsurf)     ! axis and mag parameters, use only psi0

    CALL pages
    CALL grd1d(1,rhot,psimg,nrmax+1,nrmax+1,1,'@psimg@')
    CALL grd1d(2,rhot,Fpsig,nrmax+1,nrmax+1,1,'@Fpsig@')
    CALL grd1d(3,rhot,ppsi, nrmax+1,nrmax+1,1,'@ppsi@')
    CALL grd1d(4,rhot,qpsi, nrmax+1,nrmax+1,1,'@qpsi@')
    CALL pagee
    psi0 = psi0/(2*pi)
    do nr = 1, nrmax+1
      Fpsig(nr) = Fpsig(nr)/(2*pi)
      psimg(nr) = psimg(nr)/(2*pi)
      do nthp = 1, nthpmax
        Babs(nthp,nr) = sqrt(Bt(nthp,nr)**2+Bp(nthp,nr)**2)
      end do
    end do

    write(*,*)"FP------------------------------------"

    call fow_calculate_equator_variable(ierr)

  end subroutine

  subroutine fow_calculate_equator_variable(ierr)
    use fowcomm
    use fpcomm,only:rkind,nrmax,nthmax,npmax
    implicit none
    integer,intent(inout) :: ierr
    real(rkind),allocatable,dimension(:) :: x,f,fx,g,gx,h1,h2,h1x,h2x
    real(rkind),allocatable,dimension(:,:) :: U,V,W1,W2
    real(rkind) :: dps0,dps
    integer :: nr,i,j

    ierr = 0

    allocate(x(nrmax+1),f(nrmax+1),fx(nrmax+1),U(4,nrmax+1))
    allocate(g(nrmax+1),gx(nrmax+1),V(4,nrmax+1))
    allocate(h1(nrmax+1),h1x(nrmax+1),W1(4,nrmax+1))
    allocate(h2(nrmax+1),h2x(nrmax+1),W2(4,nrmax+1))

    do nr = 1,nrmax+1
      x(nr) = (nr-1)*1.d0
      f(nr) = psimg(nr)
      g(nr) = Fpsig(nr)
      h1(nr) = Babs(1,nr)
      h2(nr) = Babs((nthpmax+1)/2,nr)
    end do

    call SPL1D(x,f,fx,U,nrmax+1,0,IERR)
    call SPL1D(x,g,gx,V,nrmax+1,0,IERR)
    call SPL1D(x,h1,h1x,W1,nrmax+1,0,IERR)
    call SPL1D(x,h2,h2x,W2,nrmax+1,0,IERR)

    Boutg(1) = h1(1)
    Bing(1)  = h2(1)

    do nr = 2,nrmax+1
      if ( psi0 <= cal_spl1D(U(:,nr),1.d0) ) then
        dps0 = 1.d0
        call newton_spl1D(U(:,nr),psi0,dps0)
        dps = (x(nr-1)+dps0)/nrmax

        ! define psimg and Fpsimg
        do i = 2, nrmax+1 ! i is label of psimg and Fpsig
          do j = 2, nrmax+1 
            if ( x(j-1) < (i-1)*dps .and. (i-1)*dps <= x(j) ) then
              psimg(i) = cal_spl1D(U(:,j), (i-1)*dps-x(j-1))
              Fpsig(i) = cal_spl1D(V(:,j), (i-1)*dps-x(j-1))
              Boutg(i) = cal_spl1D(W1(:,j), (i-1)*dps-x(j-1))
              Bing(i)  = cal_spl1D(W2(:,j), (i-1)*dps-x(j-1))
              exit
            end if
          end do
        end do

        ! define psim and Fpsim
        do i = 1, nrmax ! i is label of psim and Fpsi
          do j = 2, nrmax+1 
            if ( x(j-1) < (i-0.5d0)*dps .and. (i-0.5d0)*dps <= x(j) ) then
              psim(i) = cal_spl1D(U(:,j), (i-0.5d0)*dps-x(j-1))
              Fpsi(i) = cal_spl1D(V(:,j), (i-0.5d0)*dps-x(j-1))
              Bout(i) = cal_spl1D(W1(:,j), (i-0.5d0)*dps-x(j-1))
              Bin(i)  = cal_spl1D(W2(:,j), (i-0.5d0)*dps-x(j-1))
              exit
            end if
          end do
        end do

        exit
      end if
    end do

  end subroutine fow_calculate_equator_variable

  function cal_spl1D(U,dx) result(cal)
    use fpcomm,only:rkind
    implicit none
    real(rkind) :: cal
    real(rkind),intent(in) :: U(4),dx
    integer :: k
    cal = 0.d0
    do k = 1, 4
      cal = cal+U(k)*dx**(k-1)
    end do
    
  end function

  subroutine newton_spl1D(U,a,x)
    ! U4*x**3+U3*x**2+U2*x+U1 = a
    use fpcomm,only:rkind
    implicit none
    real(rkind),intent(in) :: U(4),a
    real(rkind),intent(inout):: x
    integer :: i = 0,j
    real(rkind) :: x0,e = 1.d10,eps = 1.d-18,V(4),d

    V = U
    V(1) = V(1)-a

    do
      d = (V(4)*x**3+V(3)*x**2+V(2)*x+V(1))/(3.d0*V(4)*x**2+2.d0*V(3)*x+V(2))
      if(abs(d)<eps)exit

      x0 = x
      x = x0-d

      if(i>e)then
        write(*,*)"Do not converge at subroutine solve_spl1D"
        stop
      end if
    end do

  end subroutine newton_spl1D

  subroutine mesh_to_grid1D(f,g)
    use fpcomm,only:rkind
    implicit none
    real(rkind) :: f(:),g(:)
    real(rkind),allocatable :: x(:),fx(:),U(:,:)
    integer :: i,j,imax,jmax,IERR = 0,k

    imax = size(f)
    jmax = size(g)

    if(imax/=jmax-1)then
      write(*,*)"imax/ = jmax-1 at subroutine mesh_to_grid1D"
      STOP
    end if

    allocate(x(imax),fx(imax),U(4,imax))

    do i = 1,imax
      x(i) = i*1.d0
    end do

    call SPL1D(x,f,fx,U,imax,0,IERR)

    do j = 2,jmax-1
      g(j) = 0.d0
      do k = 1,4
        g(j) = g(j)+U(k,j)*0.5d0**(k-1)
      end do
    end do

    g(1) = 0.d0
    g(jmax) = 0.d0
    do k = 1,4
      g(1) = g(1)+U(k,2)*(-0.5d0)**(k-1)
      g(jmax) = g(jmax)+U(k,imax)*(1.5d0)**(k-1)
    end do

    deallocate(x,fx,U)

  end subroutine mesh_to_grid1D

end module fowprep
