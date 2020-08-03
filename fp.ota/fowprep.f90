module fowprep
  PRIVATE
  PUBLIC :: fow_prep
contains

  subroutine fow_prep

    use fowcomm

    implicit none

    integer :: nze,nps,nrg,nzg,ierr = 0
    real(8),allocatable :: fx(:,:),fy(:,:),fxy(:,:),R(:),Z(:)

    do nze = 1,nzemax
      zeta(nze) = (nze-0.5d0)/nzemax*(-2.d0)+1.d0
      zetag(nze) = (nze-1.d0)/nzemax*(-2.d0)+1.d0
    end do
    zetag(nzemax+1) = -1.d0

    call mesh_to_grid1D(Bout,Boutg)
    call mesh_to_grid1D(Bin,Bing)

    do nps = 1,npsmax
      do nze = 1,nzemax
        if(zeta(nze)<0.d0)then
          BBm(nze,nps) = Bin(nps)
        else
          BBm(nze,nps) = Bout(nps)
        end if
      end do
    end do

    allocate(fx(nrgmax,nzgmax),fy(nrgmax,nzgmax),fxy(nrgmax,nzgmax),R(nrgmax),Z(nzgmax))

    do nrg = 1,nrgmax
      R(nrg) = nrg*1.d0
    end do
    do nzg = 1,nzgmax
      Z(nzg) = nzg*1.d0
    end do

    call SPL2D(R,Z,Brz,fx,fy,fxy,UBspl,nrgmax,nrgmax,nzgmax,0,0,ierr)
    call SPL2D(R,Z,Frz,fx,fy,fxy,UFspl,nrgmax,nrgmax,nzgmax,0,0,ierr)

    call fow_eqload(ierr)

    deallocate(fx,fy,fxy,R,Z)

  end subroutine fow_prep

  subroutine fow_eqload(ierr)
    use fowcomm
    use fpcomm,only:rkind
    implicit none
    integer,intent(out):: ierr
    character(len = 80) :: line
    real(rkind),allocatable,dimension(:) :: ppsi,qpsi,vpsi,rlen,ritpsi,rhot
    real(rkind),allocatable,dimension(:,:) :: Br,Bz,Bp,Bt
    integer :: nps,nthp

    allocate(ppsi(npsmax+1),qpsi(npsmax+1),vpsi(npsmax+1),rlen(npsmax+1),ritpsi(npsmax+1)&
            ,rhot(npsmax+1))
    allocate(Br(nthpmax,npsmax+1),Bz(nthpmax,npsmax+1),Bp(nthpmax,npsmax+1),Bt(nthpmax,npsmax+1))

    ierr = 0
    write(*,*)"FP------------------------------------"

    call eqload(3,knameq,ierr)
    if(ierr.ne.0) return

    write(line,'(a)') 'mdleqc = 1'   ! set boozer poloidal angle
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nrmax = ',npsmax+1
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nthmax = ',nthpmax
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nsumax = ',npsmax+1
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    call eqcalq(ierr)
    if(ierr.ne.0) return

    call eqgetp(rhot,psimg,npsmax+1)  ! normalized psit radius
    call eqgetqn(ppsi,qpsi,Fpsig,vpsi,rlen,ritpsi,npsmax+1) ! flux functions
    CALL eqgetbb(Br,Bz,Bp,Bt,nthpmax,nthpmax,npsmax+1) ! mag field

    do nps = 1, npsmax+1
      do nthp = 1, nthpmax
        Babs(nthp,nps) = sqrt(Bt(nthp,nps)**2+Bp(nthp,nps)**2)
      end do
    end do

    do nps=1,npsmax+1
      Fpsig(nps)=Fpsig(nps)/(2*pi)
    end do

    write(*,*)"FP------------------------------------"
    call fow_calcurate_equator_variable(ierr)

  end subroutine

  subroutine fow_calcurate_equator_variable(ierr)
    use fowcomm
    use fpcomm,only:rkind
    implicit none
    integer,intent(inout) :: ierr
    real(rkind),allocatable,dimension(:) :: x,f,fx,g,gx,h1,h2,h1x,h2x
    real(rkind),allocatable,dimension(:,:) :: U,V,W1,W2
    real(rkind) :: dps0,dps
    integer :: nps,i,j

    ierr = 0

    allocate(x(npsmax+1),f(npsmax+1),fx(npsmax+1),U(4,npsmax+1))
    allocate(g(npsmax+1),gx(npsmax+1),V(4,npsmax+1))
    allocate(h1(npsmax+1),h1x(npsmax+1),W1(4,npsmax+1))
    allocate(h2(npsmax+1),h2x(npsmax+1),W2(4,npsmax+1))

    do nps = 1,npsmax+1
      x(nps) = (nps-1)*1.d0
      f(nps) = psimg(nps)
      g(nps) = Fpsig(nps)
      h1(nps) = Babs(1,nps)
      h2(nps) = Babs((nthpmax+1)/2,nps)
    end do

    call SPL1D(x,f,fx,U,npsmax+1,0,IERR)
    call SPL1D(x,g,gx,V,npsmax+1,0,IERR)
    call SPL1D(x,h1,h1x,W1,npsmax+1,0,IERR)
    call SPL1D(x,h2,h2x,W2,npsmax+1,0,IERR)

    Boutg(1) = h1(1)
    Bing(1)  = h2(1)

    do nps = 2,npsmax+1
      if ( psi0 <= cal_spl1D(U(:,nps),1.d0) ) then
        dps0 = 1.d0
        call newton_spl1D(U(:,nps),psi0,dps0)
        dps = (x(nps-1)+dps0)/npsmax

        ! define psimg and Fpsimg
        do i = 2, npsmax+1 ! i is label of psimg and Fpsig
          do j = 2, npsmax+1 
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
        do i = 1, npsmax ! i is label of psim and Fpsi
          do j = 2, npsmax+1 
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

  end subroutine fow_calcurate_equator_variable

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
    implicit none
    real(8),intent(in) :: U(4),a
    real(8),intent(inout):: x
    integer :: i = 0,j
    real(8) :: x0,e = 1.d10,eps = 1.d-18,V(4),d

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
    implicit none
    real(8) :: f(:),g(:)
    real(8),allocatable :: x(:),fx(:),U(:,:)
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
