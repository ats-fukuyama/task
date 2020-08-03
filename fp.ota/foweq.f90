module foweq

contains

  subroutine interfaceEQandFP(psifp,psirzfp,Ffp,Frz,Bout,Bin,Brzfp,psi_max)

    use fpwrite

    include '../eq/eqcomc.inc'

    double precision,intent(out) :: psifp(:),Ffp(:),Brzfp(:,:),psirzfp(:,:),Bout(:),Bin(:),Frz(:,:),psi_max
    integer :: nrmfp,nrgmfp,nzgmfp,nr,nrg,nzg,nsg
    double precision,allocatable :: psitmp(:),Ftmp(:),Btmp(:,:),psirztmp(:,:)

    ierr=0

    call eqinit
    open(1996,file='../eq/eqparm', status='old')
    call eqnlin(1996,IST,IERR)
    close(1996)
    write(*,*)'read namelist /EQ/ IOSTAT=',IST
    ! READ(NID,EQ,IOSTAT=IST,ERR=9800,END=9900)

    !     ***** solve GS equation *****
    !
    call eqcalc(ierr)

    !     ***** calculate eqilibrium quantities *****
    !
    call eqcalq(ierr)

    allocate(Btmp(nrgmax,nzgmax),Ftmp(nsgmax),psitmp(nsgmax),psirztmp(nrgmax,nzgmax))
    do nrg=1,nrgmax
      do nzg=1,nzgmax
        call getrz(rg(nrg),zg(nzg),0,br,bz,bt,rhon)
        Btmp(nrg,nzg)=sqrt(br**2+bz**2+bt**2)
        psirztmp(nrg,nzg)=psirz(nrg,nzg)/(-psi0)
      end do
    end do

    do nsg=1,nsgmax
       psipnl=1.d0-psi(1,nsg)/psi0
       call eqfpsi(psipnl,fpsi,dfpsi)
       Ftmp(nsg)=FPSI/(2*pi)
       psitmp(nsg)=psi(1,nsg)-0.5d0*DTG&
             *(-3.d0*psi(1,nsg)+4.d0*psi(2,nsg)-psi(3,nsg))/(2.d0*DTG)
    end do
    psi=psi/psi0*(-1.d0)
    psitmp=psitmp/psi0*(-1.d0)

    nrmfp=size(psifp)
    nrgmfp=size(psirzfp,1)
    nzgmfp=size(psirzfp,2)

    call spline_interporation1D_to_fp(psitmp,psifp)
    call spline_interporation1D_to_fp(Ftmp,Ffp)
    call spline_interporation2D_to_fp(psirztmp,psirzfp)
    call spline_interporation2D_to_fp(Btmp,Brzfp)

    call calcurate_Bout_Bin(psirztmp,Btmp,Bout,Bin)

    call calrurate_Frz(Frz,Ffp,psirzfp,psifp)

    psi_max=psi0*(-1.d0)

    deallocate(Btmp,Ftmp,psitmp,psirztmp)

  end subroutine interfaceEQandFP

  subroutine calcurate_Bout_Bin(psirz,Brz,Bout,Bin)

    use fpwrite

    implicit none

    real(8),intent(in) :: psirz(:,:),Brz(:,:)
    real(8),intent(inout):: Bout(:),Bin(:)
    integer :: nrgmax,nzgmax,nrg,nzg,nrmax,nr

    integer :: naxis(2),nwallout,nwallin,min_index,ierr=0,k,l,nrin,nrout
    real(8) :: axis(2),dR,dZ,rwallout,rwallin,drin,drout,drfp
    real(8),allocatable :: U2D(:,:,:,:),dBdr(:,:),dBdz(:,:),dBdrdz(:,:),R(:),Z(:)&
                          ,dpsdr(:,:),dpsdz(:,:),dpsdrdz(:,:),psieq(:),U1D(:,:),dpsieqdr(:)

    nrmax=size(Bout)
    nrgmax=size(psirz,1)
    nzgmax=size(psirz,2)

    dR=1.d0/nrgmax
    dZ=1.d0/nrgmax

    allocate(R(nrgmax),Z(nzgmax))

    do nrg=1,nrgmax
      R(nrg)=(nrg*1.d0-0.5d0)*dR
    end do
    do nzg=1,nzgmax
      Z(nzg)=(nzg*1.d0-0.5d0)*dZ
    end do

    allocate(U2D(4,4,nrgmax,nzgmax),dpsdr(nrgmax,nzgmax)&
            ,dpsdz(nrgmax,nzgmax),dpsdrdz(nrgmax,nzgmax))
    allocate(psieq(nrgmax),U1D(4,nrgmax),dpsieqdr(nrgmax))

    call SPL2D(r,z,psirz,dpsdr,dpsdz,dpsdrdz,U2D,nrgmax,nrgmax,nzgmax,0,0,ierr)

    ! calcurate axis(R,Z)
    naxis=minloc(psirz)
    if(dpsdr(naxis(1),naxis(2))<0.d0)naxis(1)=naxis(1)+1
    if(dpsdz(naxis(1),naxis(2))<0.d0)naxis(2)=naxis(2)+1

    axis(1)=dr
    axis(2)=dz
    call get_minimal_location(U2D(:,:,naxis(1),naxis(2)),axis(:))

    psieq(1)=func_spl2D(U2D(:,:,2,naxis(2)),0.d0,axis(2))
    do nrg=2,nrgmax
      psieq(nrg)=func_spl2D(U2D(:,:,nrg,naxis(2)),dR,axis(2))
    end do

    call SPL1D(r,psieq,dpsieqdr,U1D,nrgmax,0,IERR)

    ! calcurate wall point
    do nrg=1,nrgmax
      if(psieq(nrg)<=0.d0)then
        nwallin=nrg
        rwallin=dr
        call solve_spl1D(U1D(:,nrg),0.d0,rwallin)
        rwallin=R(nwallin-1)+rwallin
        exit
      end if
    end do

    do nrg=nrgmax,1,-1
      if(psieq(nrg-1)<=0.d0)then
        nwallout=nrg
        rwallout=dr
        call solve_spl1D(U1D(:,nrg),0.d0,rwallout)
        rwallout=rwallout+R(nwallout-1)
        exit
      end if
    end do

    deallocate(U2D,dpsdr,dpsdz,dpsdrdz)
    deallocate(U1D,psieq,dpsieqdr)

    allocate(U2D(4,4,nrgmax,nzgmax),dBdr(nrgmax,nzgmax)&
            ,dBdz(nrgmax,nzgmax),dBdrdz(nrgmax,nzgmax))
    call SPL2D(r,z,Brz,dBdR,dBdZ,dBdRdZ,U2D,nrgmax,nrgmax,nzgmax,0,0,ierr)

    drfp=(R(naxis(1)-1)+axis(1)-rwallin)/nrmax
    do nr=1,nrmax
      drin=rwallin+drfp*nr-0.5d0*drfp
      nrin=int(drin/dr+0.5)+1
      Bin(nrmax-nr+1)=func_spl2D(U2D(:,:,nrin,naxis(2)),drin-R(nrin-1),axis(2))
    end do

    drfp=(rwallout-R(naxis(1)-1)-axis(1))/nrmax
    do nr=1,nrmax
      drout=R(naxis(1)-1)+axis(1)+drfp*nr-0.5d0*drfp
      nrout=int(drout/dr+0.5)+1
      Bout(nr)=func_spl2D(U2D(:,:,nrout,naxis(2)),drout-R(nrout-1),axis(2))
    end do

    deallocate(U2D,dBdr,dBdz,dBdrdz,R,Z)

  end subroutine calcurate_Bout_Bin

  subroutine calrurate_Frz(Frz,Fpsi,psirz,psi)
    implicit none
    real(8),intent(inout) :: Frz(:,:)
    real(8),intent(in) :: psirz(:,:),psi(:),Fpsi(:)
    integer :: nr,nrg,nzg,nrmax,nrgmax,nzgmax

    integer :: ierr=0,nd,add
    real(8) :: dR,d,F_add(-2:1000),psi_add(-2:1000)
    real(8),allocatable :: UF(:,:),Upsi(:,:),dFdr(:),dpsidr(:),R(:),psi_added(:),F_added(:)

    nrmax=size(Fpsi)
    nrgmax=size(Frz,1)
    nzgmax=size(Frz,2)
    dR=1.d0/nrmax

    do nr=1,1000
      psi_add(nr)=-1.d0
      F_add(nr)=-1.d0
    end do

    psi_add(-2)=psi(nrmax-2)
    F_add(-2)=Fpsi(nrmax-2)
    psi_add(-1)=psi(nrmax-1)
    F_add(-1)=Fpsi(nrmax-1)
    psi_add(0)=psi(nrmax)
    F_add(0)=Fpsi(nrmax)
    add=0

    do
      if(maxval(psirz)<=maxval(psi_add))exit
      psi_add(add+1)=psi_add(add)+dR*(3.d0*psi_add(add)-4.d0*psi_add(add-1)+psi_add(add-2))/(2.d0*dR)
      F_add(add+1)=F_add(add)+dR*(3.d0*F_add(add)-4.d0*F_add(add-1)+F_add(add-2))/(2.d0*dR)
      add=add+1
    end do

    allocate(psi_added(nrmax+add),F_added(nrmax+add))
    allocate(UF(4,nrmax+add),Upsi(4,nrmax+add),dFdr(nrmax+add),dpsidr(nrmax+add),R(nrmax+add))

    do nr=1,nrmax+add
      R(nr)=(nr*1.d0-0.5d0)*dR
    end do

    do nr=1,nrmax+add
      if(nr<=nrmax)then
        psi_added(nr)=psi(nr)
        F_added(nr)=Fpsi(nr)
      else
        psi_added(nr)=psi_add(nr-nrmax)
        F_added(nr)=F_add(nr-nrmax)
      end if
    end do

    call SPL1D(R,F_added,dFdr,UF,nrmax+add,0,IERR)
    call SPL1D(R,psi_added,dpsidr,Upsi,nrmax+add,0,IERR)

    do nzg=1,nzgmax
      do nrg=1,nrgmax
        do nr=1,nrmax+add
          if(psirz(nrg,nzg)<=psi_added(nr))then
            nd=nr
            exit
          end if
        end do
        d=dr
        call solve_spl1D(Upsi(:,nd),psirz(nrg,nzg),d)
        Frz(nrg,nzg)=UF(4,nd)*d**3+UF(3,nd)*d**2+UF(2,nd)*d+UF(1,nd)
      end do
    end do


  end subroutine calrurate_Frz

  subroutine spline_interporation1D_to_fp(f,g)

    implicit none

    double precision,intent(in) :: f(:)
    double precision,intent(out):: g(:)
    integer :: imax,jmax,IERR=0,i,j,k
    double precision,allocatable :: U1D(:,:),x(:),fx(:),y(:),f2(:)
    double precision :: r,dx,dy,fx2

    imax=size(f)
    jmax=size(g)
    allocate(x(imax),fx(imax),U1D(4,imax),y(jmax))

    do i=1,imax
      x(i)=(i*1.d0-0.5d0)/imax
    end do
    do j=1,jmax
      y(j)=(j*1.d0-0.5d0)/jmax
    end do

    dx=1.d0/imax
    dy=1.d0/jmax

    fx(1)   =(-3.d0*f(1)+4.d0*f(2)-f(3))/(2.d0*dx)
    fx(imax)=(3.d0*f(imax)-4.d0*f(imax-1)+f(imax-2))/(2.d0*dx)

    call SPL1D(x,f,fx,U1D,imax,3,IERR)

    allocate(f2(imax+1))
    fx2=(2.d0*f(1)-5.d0*f(2)+4.d0*f(3)-f(4))/(dx**2)
    f2(1)=f(1)+fx(1)*(-0.5d0)*dx+0.5d0*fx2*(0.5d0*dx)**2
    fx2=(2.d0*f(imax)-5.d0*f(imax-1)+4.d0*f(imax-2)-f(imax-3))/(dx**2)
    f2(imax+1)=f(imax)+fx(imax)*0.5d0*dx+0.5d0*fx2*(0.5d0*dx)**2
    do i=2,imax
      f2(i)=0.d0
      do k=1,4
        f2(i)=f2(i)+U1D(k,i)*(0.5d0*dx)**(k-1)
      end do
    end do

    deallocate(x,fx,U1D)
    allocate(x(imax+1),fx(imax+1),U1D(4,imax+1))

    do i=1,imax+1
      x(i)=(i*1.d0-1.d0)/imax
    end do
    fx(1)   =(-3.d0*f2(1)+4.d0*f2(2)-f2(3))/(2.d0*dx)
    fx(imax+1)=(3.d0*f2(imax+1)-4.d0*f2(imax)+f2(imax-2))/(2.d0*dx)

    call SPL1D(x,f2,fx,U1D,imax+1,3,IERR)

    do j=1,jmax
      do i=2,imax+1
        if(x(i-1)<y(j).and.y(j)<=x(i))then
          r=dx-(x(i)-y(j))
          g(j)=0.d0
          do k=1,4
            g(j)=g(j)+U1D(k,i)*r**(k-1)
          end do
          exit
        end if
      end do
    end do

    deallocate(x,fx,U1D,y,f2)

  end subroutine spline_interporation1D_to_fp

  subroutine spline_interporation2D_to_fp(f,g)

    implicit none

    double precision,intent(in) :: f(:,:)
    double precision,intent(out):: g(:,:)
    integer :: imax,jmax,nmax,mmax,IERR=0,i,j,n,m,k,l,ir,jz
    double precision,allocatable :: U2D(:,:,:,:),x(:),y(:),fx(:,:),fy(:,:),fxy(:,:)&
                                    ,f2(:,:),r(:),z(:),fx2(:,:),fy2(:,:)
    double precision :: dx,dy,dr,dz,pr,pz,di,dj

    imax=size(f,1)
    jmax=size(f,2)
    nmax=size(g,1)
    mmax=size(g,2)

    allocate(x(imax),y(jmax),r(nmax),z(mmax))

    dx=1.d0/imax
    dy=1.d0/jmax
    dr=1.d0/nmax
    dz=1.d0/mmax

    do i=1,imax
      x(i)=(i*1.d0-0.5d0)*dx
    end do
    do j=1,jmax
      y(j)=(j*1.d0-0.5d0)*dy
    end do
    do n=1,nmax
      r(n)=(n*1.d0-0.5d0)*dr
    end do
    do m=1,mmax
      z(m)=(m*1.d0-0.5d0)*dz
    end do

    allocate(fx(imax,jmax),fy(imax,jmax),fxy(imax,jmax),U2D(4,4,imax,jmax))

    do i=1,imax
      fy(i,1)    = (-3.d0*f(i,1)+4.d0*f(i,2)-f(i,3))/(2.d0*dy)
      fy(i,jmax) = (3.d0*f(i,jmax)-4.d0*f(i,jmax-1)+f(i,jmax-2))/(2.d0*dy)
    end do
    do j=1,jmax
      fx(1,j)    = (-3.d0*f(1,j)+4.d0*f(2,j)-f(3,j))/(2.d0*dx)
      fx(imax,j) = (3.d0*f(imax,j)-4.d0*f(imax-1,j)+f(imax-2,j))/(2.d0*dx)
    end do
    fxy(1,1)      = (-3.d0*fx(1,1)+4.d0*fx(1,2)-fx(1,3))/(2.d0*dy)
    fxy(1,jmax)   = (3.d0*fx(1,jmax)-4.d0*fx(1,jmax-1)+fx(1,jmax-2))/(2.d0*dy)
    fxy(imax,1)   = (-3.d0*fx(imax,1)+4.d0*fx(imax,2)-fx(imax,3))/(2.d0*dy)
    fxy(imax,jmax)= (3.d0*fx(imax,jmax)-4.d0*fx(imax,jmax-1)+fx(imax,jmax-2))/(2.d0*dy)

    call SPL2D(x,y,f,fx,fy,fxy,U2D,imax,imax,jmax,3,3,ierr)

    allocate(f2(imax+1,jmax+1),fx2(imax,jmax),fy2(imax,jmax))

    ! calcurate 2nd derivatives of f along the edge, fx2 and fy2
    do i=1,imax
      fy2(i,1)    = (2.d0*f(i,1)-5.d0*f(i,2)+4.d0*f(i,3)-f(i,4))/(dy**2)
      fy2(i,jmax) = (2.d0*f(i,jmax)-5.d0*f(i,jmax-1)+4.d0*f(i,jmax-2)-f(i,jmax-3))/(dy**2)
    end do
    do j=1,jmax
      fx2(1,j)    = (2.d0*f(1,j)-5.d0*f(2,j)+4.d0*f(3,j)-f(4,j))/(dx**2)
      fx2(imax,j) = (2.d0*f(imax,j)-5.d0*f(imax-1,j)+4.d0*f(imax-2,j)-f(imax-3,j))/(dx**2)
    end do

    ! extrapolate f to integer grid points, using the Taylor expansion
    di=0.5d0*dx
    dj=0.5d0*dy
    do j=1,jmax
      f2(1,j)    = f(1,j)-fx(1,j)*di-fy(1,j)*dj&
                  +0.5d0*fx2(1,j)*di**2+0.5d0*fy2(1,j)*dj**2&
                  +fxy(1,j)*di*dj
      f2(imax+1,j)= f(imax,j)+fx(imax,j)*di-fy(imax,j)*dj&
                  +0.5d0*fx2(imax,j)*di**2+0.5d0*fy2(imax,j)*dj**2&
                  -fxy(imax,j)*di*dj
    end do
    f2(1,jmax+1) = f(1,jmax)-fx(1,jmax)*di+fy(1,jmax)*dj&
                +0.5d0*fx2(1,jmax)*di**2+0.5d0*fy2(1,jmax)*dj**2&
                -fxy(1,jmax)*di*dj
    f2(imax+1,jmax+1)= f(imax,jmax)+fx(imax,jmax)*di+fy(imax,jmax)*dj&
                +0.5d0*fx2(imax,jmax)*di**2+0.5d0*fy2(imax,jmax)*dj**2&
                +fxy(imax,jmax)*di*dj

    do i=2,imax
      f2(i,1)    = f(i,1)-fx(i,1)*di-fy(i,1)*dj&
                  +0.5d0*fx2(i,1)*di**2+0.5d0*fy2(i,1)*dj**2&
                  +fxy(i,1)*di*dj
      f2(i,jmax+1)= f(i,jmax)-fx(i,jmax)*di+fy(i,jmax)*dj&
                  +0.5d0*fx2(i,jmax)*di**2+0.5d0*fy2(i,jmax)*dj**2&
                  -fxy(i,jmax)*di*dj
    end do

    ! interpolate f to integer grid points, using spline coeffcients, U2D
    do i=2,imax
      do j=2,jmax
        f2(i,j)=0.d0
        do k=1,4
          do l=1,4
            f2(i,j)=f2(i,j)+U2D(k,l,i,j)*(0.5d0*dx)**(k-1)*(0.5d0*dy)**(l-1)
          end do
        end do
      end do
    end do

    ! calcurate derivatives of f2 along the edge
    deallocate(x,y,fx,fy,fxy,U2D)
    allocate(x(imax+1),y(jmax+1),fx(imax+1,jmax+1),fy(imax+1,jmax+1),fxy(imax+1,jmax+1),U2D(4,4,imax+1,jmax+1))
    do i=1,imax+1
      fy(i,1)    = (-3.d0*f2(i,1)+4.d0*f2(i,2)-f2(i,3))/(2.d0*dy)
      fy(i,jmax+1) = (3.d0*f2(i,jmax+1)-4.d0*f2(i,jmax)+f2(i,jmax-1))/(2.d0*dy)
    end do
    do j=1,jmax+1
      fx(1,j)    = (-3.d0*f2(1,j)+4.d0*f2(2,j)-f2(3,j))/(2.d0*dx)
      fx(imax+1,j) = (3.d0*f2(imax+1,j)-4.d0*f2(imax,j)+f2(imax-1,j))/(2.d0*dx)
    end do
    fxy(1,1)      = (-3.d0*fx(1,1)+4.d0*fx(1,2)-fx(1,3))/(2.d0*dy)
    fxy(1,jmax+1)   = (3.d0*fx(1,jmax+1)-4.d0*fx(1,jmax)+fx(1,jmax-1))/(2.d0*dy)
    fxy(imax+1,1)   = (-3.d0*fx(imax+1,1)+4.d0*fx(imax+1,2)-fx(imax+1,3))/(2.d0*dy)
    fxy(imax+1,jmax+1)= (3.d0*fx(imax+1,jmax+1)-4.d0*fx(imax+1,jmax)+fx(imax+1,jmax-1))/(2.d0*dy)

    do i=1,imax+1
      x(i)=(i*1.d0-1.d0)*dx
    end do
    do j=1,jmax+1
      y(j)=(j*1.d0-1.d0)*dy
    end do

    call SPL2D(x,y,f2,fx,fy,fxy,U2D,imax+1,imax+1,jmax+1,3,3,ierr)

    do n=1,nmax
      do m=1,mmax
        g(n,m)=0.d0
        do i=2,imax+1
          if(x(i-1)<r(n).and.r(n)<=x(i))then
            pr=dx-(x(i)-r(n))
            ir=i
            exit
          endif
        end do
        do j=2,jmax+1
          if(y(j-1)<z(m).and.z(m)<=y(j))then
            pz=dy-(y(j)-z(m))
            jz=j
            exit
          end if
        end do
        do k=1,4
          do l=1,4
            g(n,m)=g(n,m)+U2D(k,l,ir,jz)*pr**(k-1)*pz**(l-1)
          end do
        end do
      end do
    end do

    deallocate(x,y,fx,fy,fxy,U2D,f2,r,z,fx2,fy2)

  end subroutine spline_interporation2D_to_fp

  subroutine get_minimal_location(U,v)
    ! Newton's method
    implicit none
    real(8),intent(in) :: U(4,4) ! spl2D coefficient
    ! func_slp2D(U,R,Z) is minimal value at v=(v_R,v_Z), input as initial value
    real(8),intent(inout) :: v(2)
    real(8) :: d(2),H_inv(2,2),H(2,2),grad(2),eps=1.d-12
    integer :: i,j,k

    k=0
    do
      grad=grad_func_spl2D(U,v(1),v(2))
      if(abs(grad(1))<=eps.and.abs(grad(2))<=eps)exit

      H=Hesse_func_spl2D(U,v(1),v(2))
      H_inv(1,1)=H(2,2)
      H_inv(2,2)=H(1,1)
      H_inv(2,1)=-1.d0*H(2,1)
      H_inv(1,2)=-1.d0*H(1,2)
      H_inv=H_inv/(H(1,1)*H(2,2)-H(1,2)*H(2,1))

      do i=1,2
        d(i)=0.d0
        do j=1,2
          d(i)=d(i)+H_inv(i,j)*grad(j)
        end do
      end do

      do i=1,2
        v(i)=v(i)-d(i)
      end do

      k=k+1
      if(k>=1.d8)then
        write(*,*)"Do not converge by Newton's method"
        stop
      end if
    end do

  end subroutine

  subroutine solve_spl1D(U,a,x)
    ! U4*x**3+U3*x**2+U2*x+U1=a
    implicit none
    real(8),intent(in) :: U(4),a
    real(8),intent(inout):: x
    integer :: i=0,j
    real(8) :: x0,e=1.d10,eps=1.d-12,V(4),d

    V=U
    V(1)=V(1)-a

    do
      d=(V(4)*x**3+V(3)*x**2+V(2)*x+V(1))/(3.d0*V(4)*x**2+2.d0*V(3)*x+V(2))
      if(abs(d)<eps)exit

      x0=x
      x=x0-d

      if(i>e)then
        write(*,*)"Do not converge at subroutine solve_spl1D"
        stop
      end if
    end do

  end subroutine solve_spl1D

  function func_spl2D(U,x,y)
    real(8),intent(in) :: U(4,4),x,y
    real(8) :: func_spl2D
    real(8) :: v
    integer :: i,j

    v=0.d0
    do i=1,4
      do j=1,4
        v=v+U(i,j)*x**(i-1)*y**(j-1)
      end do
    end do

    func_spl2D=v

  end function func_spl2D

  function grad_func_spl2D(U,x,y)
    implicit none
    real(8),intent(in) :: U(4,4),x,y
    real(8) :: grad_func_spl2D(2)
    real(8) :: v,w
    integer :: i,j

    v=0.d0
    w=0.d0
    do i=1,4
      do j=1,4
        v=v+(i-1)*U(i,j)*x**(i-2)*y**(j-1)
        w=w+(j-1)*U(i,j)*x**(i-1)*y**(j-2)
      end do
    end do

    grad_func_spl2D(1)=v
    grad_func_spl2D(2)=w

  end function grad_func_spl2D

  function Hesse_func_spl2D(U,x,y)
    implicit none
    real(8),intent(in) :: U(4,4),x,y
    real(8) :: Hesse_func_spl2D(2,2)
    real(8) :: v,w,z
    integer :: i,j

    v=0.d0
    w=0.d0
    z=0.d0
    do i=1,4
      do j=1,4
        v=v+(i-1)*(i-2)*U(i,j)*x**(i-3)*y**(j-1)
        w=w+(j-1)*(j-2)*U(i,j)*x**(i-1)*y**(j-3)
        z=z+(i-1)*(j-1)*U(i,j)*x**(i-2)*y**(j-2)
      end do
    end do
    Hesse_func_spl2D(1,1)=v
    Hesse_func_spl2D(2,2)=w
    Hesse_func_spl2D(1,2)=z
    Hesse_func_spl2D(2,1)=z

  end function Hesse_func_spl2D

end module
