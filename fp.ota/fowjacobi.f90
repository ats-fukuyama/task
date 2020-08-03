module fow_jacobi

contains

  subroutine fow_jacobian(J,dIdu,th,p,rg,zg,nze,nps,nsa)

    use fpcomm
    use fowcomm

    implicit none

    real(8),intent(out) :: dIdu(3,3),J
    integer,intent(in) :: nsa,nze,nps
    real(8),intent(in) :: th,p,rg,zg

    ! elements of Jacobi matrix, dIdu
    real(8) ::  dpdp,  dzedp,  dpsdp,&
                dpdth, dzedth, dpsdth,&
                dpdr,  dzedr,  dpsdr
    real(8) :: A,B,C,D,E
    real(8) :: dBdps,dFoBdps,FoB,dps,F,B !FoB=F(psi)/B(psi), F=RB_phi(rg,zg), B=B(rg,zg)
    real(8),allocatable :: Beq(:)
    integer :: ips,i,j

    REAL(8) :: AEFDL

    AEFDL=PZ(NS)*AEE

    F=0.d0
    B=0.d0
    do i=4
      do j=4
        B=B+UBspl(i,j,int(rg)+1,int(zg)+1)*(rg-int(rg)*1.d0)**i*(zg-int(zg)*1.d0)**j
        F=F+UFspl(i,j,int(rg)+1,int(zg)+1)*(rg-int(rg)*1.d0)**i*(zg-int(zg)*1.d0)**j
      end do
    end do

    allocate(Beq(npsmax))

    dpdp =1.d0
    dpdth=0.d0
    dpdr =0.d0
    dzedr=0.d0
    dpsdr=0.d0

    dps=psi0/npsmax

    do ips=1,npsmax
      if(zeta(nze)<=0.d0)then
        Beq(ips)=Bin(ips)
      else
        Beq(ips)=Bout(ips)
      end if
    end do

    FoB=Fpsi(nps)/Beq(nps)

    if(nps<=npsmax-2)then
      dBdps=(-3.d0*Beq(nps)+4.d0*Beq(nps+1)-Beq(nps+2))/(2.d0*dps)
      dFoBdps=(-3.d0*Fpsi(nps)/Beq(nps)+4.d0*Fpsi(nps+1)/Beq(nps+1)-Fpsi(nps+2)/Beq(nps+2))/(2.d0*dps)
    else
      dBdps=(3.d0*Beq(nps)-4.d0*Beq(nps-1)+Beq(nps-2))/(2.d0*dps)
      dFoBdps=(3.d0*Fpsi(nps)/Beq(nps-1)-4.d0*Fpsi(nps)/Beq(nps-1)+Fpsi(nps-2)/Beq(nps-2))/(2.d0*dps)
    end if

    A=F/B*cosm(nth)-FoB*zeta(nze)
    B=(1.d0-zeta(nze)**2)/(2.d0*zeta(nze))*dBdps/Bin(nps)
    C=dFoBdps*p*zeta(nze)-AEFDL
    D=p*sinm(nth)/B*(Fpsi(nps)*cosm(nth)/zeta(nze)-F)
    E=Bin(nps)/B*sinm(nth)*cosm(nth)/zeta(nze)

    dpsdp = A/(B*dFoBdps*p+C)
    dzedp = A/(B*dFoBdps*p+C)*B
    dpsdth= D/(B*dFoBdps*p+C)
    dzedth= D/(B*dFoBdps*p+C)*B-E

    dIdu(1,1)=dpdp
    dIdu(1,2)=dzedp
    dIdu(1,3)=dpsdp
    dIdu(2,1)=dpdth
    dIdu(2,2)=dzedth
    dIdu(2,3)=dpsdth
    dIdu(3,1)=dpdr
    dIdu(3,2)=dzedr
    dIdu(3,3)=dpsdr

    J= dIdu(1,1)*dIdu(2,2)*dIdu(3,3)&
      +dIdu(1,2)*dIdu(2,3)*dIdu(3,1)&
      +dIdu(1,3)*dIdu(2,1)*dIdu(3,2)&
      -dIdu(1,3)*dIdu(2,2)*dIdu(3,1)&
      -dIdu(1,1)*dIdu(2,3)*dIdu(3,2)&
      -dIdu(1,2)*dIdu(2,1)*dIdu(3,3)

    deallocate(Beq)

  end subroutine

  subroutine get_B_and_F(B,F,rg,zg)

  end subroutine

end module fow_jacobi
