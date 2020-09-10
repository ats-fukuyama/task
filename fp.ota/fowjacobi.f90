module fowjacobi
  implicit none
  private
  ! public ::  fow_transformation_matrix, fow_calcurate_jacobian

contains

  ! subroutine fow_transformation_matrix(J,dIdu,th,p,rg,zg,nth,nr,nsa)

  !   use fpcomm
  !   use fowcomm

  !   real(rkind),intent(out) :: dIdu(3,3),J
  !   integer,intent(in) :: nsa,nth,nr
  !   real(rkind),intent(in) :: th,p,rg,zg

  !   ! elements of Jacobi matrix, dIdu
  !   real(rkind) ::  dpdp,  dzedp,  dpsdp,&
  !               dpdth, dzedth, dpsdth,&
  !               dpdr,  dzedr,  dpsdr
  !   real(rkind) :: A,B,C,D,E
  !   real(rkind) :: dBdps,dFoBdps,FoB,dps,F,B !FoB=F(psi)/B(psi), F=RB_phi(rg,zg), B=B(rg,zg)
  !   real(rkind),allocatable :: Beq(:)
  !   integer :: ir,i,j

  !   REAL(8) :: AEFDL

  !   AEFDL=PZ(NS)*AEE

  !   F=0.d0
  !   B=0.d0
  !   do i=4
  !     do j=4
  !       B=B+UBspl(i,j,int(rg)+1,int(zg)+1)*(rg-int(rg)*1.d0)**i*(zg-int(zg)*1.d0)**j
  !       F=F+UFspl(i,j,int(rg)+1,int(zg)+1)*(rg-int(rg)*1.d0)**i*(zg-int(zg)*1.d0)**j
  !     end do
  !   end do

  !   allocate(Beq(nrmax))

  !   dpdp =1.d0
  !   dpdth=0.d0
  !   dpdr =0.d0
  !   dzedr=0.d0
  !   dpsdr=0.d0

  !   dps=psi0/nrmax

  !   do ir=1,nrmax
  !     if(xi(nth)<=0.d0)then
  !       Beq(ir)=Bin(ir)
  !     else
  !       Beq(ir)=Bout(ir)
  !     end if
  !   end do

  !   FoB=Fpsi(nr)/Beq(nr)

  !   if(nr<=nrmax-2)then
  !     dBdps=(-3.d0*Beq(nr)+4.d0*Beq(nr+1)-Beq(nr+2))/(2.d0*dps)
  !     dFoBdps=(-3.d0*Fpsi(nr)/Beq(nr)+4.d0*Fpsi(nr+1)/Beq(nr+1)-Fpsi(nr+2)/Beq(nr+2))/(2.d0*dps)
  !   else
  !     dBdps=(3.d0*Beq(nr)-4.d0*Beq(nr-1)+Beq(nr-2))/(2.d0*dps)
  !     dFoBdps=(3.d0*Fpsi(nr)/Beq(nr-1)-4.d0*Fpsi(nr)/Beq(nr-1)+Fpsi(nr-2)/Beq(nr-2))/(2.d0*dps)
  !   end if

  !   A=F/B*cosm(nth)-FoB*xi(nth)
  !   B=(1.d0-xi(nth)**2)/(2.d0*xi(nth))*dBdps/Bin(nr)
  !   C=dFoBdps*p*xi(nth)-AEFDL
  !   D=p*sinm(nth)/B*(Fpsi(nr)*cosm(nth)/xi(nth)-F)
  !   E=Bin(nr)/B*sinm(nth)*cosm(nth)/xi(nth)

  !   dpsdp = A/(B*dFoBdps*p+C)
  !   dzedp = A/(B*dFoBdps*p+C)*B
  !   dpsdth= D/(B*dFoBdps*p+C)
  !   dzedth= D/(B*dFoBdps*p+C)*B-E

  !   dIdu(1,1)=dpdp
  !   dIdu(1,2)=dzedp
  !   dIdu(1,3)=dpsdp
  !   dIdu(2,1)=dpdth
  !   dIdu(2,2)=dzedth
  !   dIdu(2,3)=dpsdth
  !   dIdu(3,1)=dpdr
  !   dIdu(3,2)=dzedr
  !   dIdu(3,3)=dpsdr

  !   J= dIdu(1,1)*dIdu(2,2)*dIdu(3,3)&
  !     +dIdu(1,2)*dIdu(2,3)*dIdu(3,1)&
  !     +dIdu(1,3)*dIdu(2,1)*dIdu(3,2)&
  !     -dIdu(1,3)*dIdu(2,2)*dIdu(3,1)&
  !     -dIdu(1,1)*dIdu(2,3)*dIdu(3,2)&
  !     -dIdu(1,2)*dIdu(2,1)*dIdu(3,3)

  !   deallocate(Beq)

  ! end subroutine fow_transformation_matrix

  subroutine get_B_and_F(B,F,rg,zg)

  end subroutine

end module fowjacobi
