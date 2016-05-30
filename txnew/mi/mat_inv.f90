module matrix_inversion
  implicit none

contains

!***************************************************************
!
!        Interface for Matrix Inversion
!
!***************************************************************

  subroutine tx_matrix_inversion(NR,ETAout,BJBSout, &
       &                         ChiNCpel,ChiNCtel,ChiNCpil,ChiNCtil, &
       &                         ddPhidpsi_in)
    use tx_commons
    use tx_interface, only : coll_freq, ftfunc
    integer(4), intent(in) :: NR
    real(8), intent(in) :: ddPhidpsi_in
    real(8), intent(out) :: ETAout, BJBSout, ChiNCpel, ChiNCtel, ChiNCpil, ChiNCtil
    integer(4) :: imodel(5), ibeam, NSMB, i, j, i1, i2
    real(8) :: fbeam, sqzfaccoef, coencc, bjsum, ftl, fac, facee, epsL, fac2
    real(8) :: vdiamg1, vohm1, vdiamg2, vohm2, renorm, rbanana
    real(8) :: smallvalue = 1.d-5
    real(8), dimension(:), allocatable :: PNsVL, PTsVL, amasL, achgL, sqzfac, xsta, coebdc
    real(8), dimension(:,:), allocatable :: ztau, coebsc, amat, bmat, cmat, alf, &
         & chipBP, chiTBP, chipPS, chiTPS

    imodel(1) = 0 ! Fast ion viscosity
    imodel(2) = 1 ! Required for NBCD
    imodel(3) = 2 ! Fast ion contribution to friction forces
    imodel(4) = 0 ! Unused
    imodel(5) = 0 ! PS contribution nil when 1

    fbeam = 0.d0 ! 0: fast ion viscosity computed, 1: nil

    if(PNBH == 0.D0 .and. PNbV(NR) < 1.D-8) then
       ibeam = 0
       NSMB  = NSM
    else
       ibeam = 1
       NSMB  = NSM + 1
    end if

    allocate(PNsVL(NSMB),PTsVL(NSMB),amasL(NSMB),achgL(NSMB),sqzfac(NSMB),xsta(NSMB))
    allocate(ztau(NSM,NSM))
    allocate(coebsc(NSMB,2),coebdc(NSM+3),amat(2*NSMB,2*NSMB),bmat(2*NSMB,2*NSMB),cmat(2*NSMB,2*NSMB),alf(2*NSMB,2*NSMB))
    allocate(chipBP(NSM,NSM),chiTBP(NSM,NSM),chipPS(NSM,NSM),chiTPS(NSM,NSM))

    amasL(1:NSM) = amas(1:NSM)
    achgL(1:NSM) = achg(1:NSM)
    PNsVL(1:NSM) = Var(NR,1:NSM)%n
    PTsVL(1:NSM) = Var(NR,1:NSM)%T

    if( NR /= 0 ) then
       epsL = epst(NR)
    else
       epsL = smallvalue
    end if
    ftl = ftfunc(epsL)

    !-- Orbit squeezing factor, whose definition is identical to that in NCLASS
    !     Note that sqzfac would become infinity at the magnetic axis due to d/dpsi(dPhi/dpsi).
    !               sqzfaccoef = 0 when NR = 0 because ddPhidpsi_in is set to zero.
    sqzfaccoef = (fipol(NR)/bb)**2*amqp*abs(ddPhidpsi_in)
    do i = 1, NSM
       fac = 1.d0
       if(NR /= 0) then
          rbanana = sqrt(epsL)*amas(i)*amqp/(abs(achg(i))*BthV(NR)) &
               &   *sqrt(2.d0*Var(NR,i)%T*rKeV/(amas(i)*amp)) ! banana widt
          if(r(NR) < rbanana) fac = 0.d0
       end if
       sqzfac(i)  = 1.d0 + fac * sqzfaccoef*amas(i)/abs(achg(i))
    end do

    !-- Effective self-collisionality; xsta, self-collision time; ztau
    do i = 1, NSM
       xsta(i) = coll_freq(NR,i,i,epsL) ! only for booth9
       do j = 1, NSM
          ztau(i,j) = 1.d0 / coll_freq(NR,i,j)
       end do
    end do

    if( ibeam == 1 ) then
       amasL(NSMB)  = amb
       achgL(NSMB)  = achgb
       PNsVL(NSMB)  = PNbV(NR)
       PTsVL(NSMB)  = PTbV(NR)
       fac = 1.d0
       if(NR /= 0) then
          rbanana = sqrt(epsL)*amb*amqp/(abs(achgb)*BthV(NR)) &
               &   *sqrt(2.d0*PTsVL(NSMB)*rKeV/(amb*amp)) ! banana width
          if(r(NR) < rbanana) fac = 0.d0
       end if
       sqzfac(NSMB) = 1.d0 + fac * sqzfaccoef*amb/abs(achgb)
       xsta(NSMB)   = 0.d0
    end if

    !-- Call Matrix Inversion (booth9)
    call booth9(coebsc,coencc,coebdc,amat,bmat,cmat,alf &
         &     ,PNsVL,PTsVL,amasL,achgL,ftl,sqzfac,xsta,ztau &
         &     ,NSM,ibeam,imodel,fbeam,epsL)

    !-- Neoclassical resistivity
    !     Note: coencc is normalized by n_e m_e/tau_ee
    i = 1
    ETAout = amas(i) * amp / (coencc * aee**2 * Var(NR,i)%n * 1.d20 * ztau(i,i))

    !-- Bootstrap current
    bjsum = 0.d0
    do i = 1, NSM
       bjsum = bjsum + Var(NR,i)%T * rKeV / abs(achg(i)) &
            & * (  coebsc(i,1) * dPsdpsi(NR,i) / Var(NR,i)%p &
            &    + coebsc(i,2) * dTsdpsi(NR,i) / Var(NR,i)%T)
    end do
    BJBSout = - fipol(NR) * Var(NR,1)%n * 1.d20 * bjsum

!    !     alternative way of computing bootstrap current
!    bjsum = 0.d0
!    do i1 = 1, NSM
!       do i2 = 1, NSM
!          bjsum = bjsum - achg(i1) * Var(NR,i1)%n &
!               &  * (  alf(i1     ,i2     ) * BVsdiag(NR,i2,1) &
!               &     - alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2))
!       end do
!    end do
!    BJBSout = bjsum * aee * 1.d20

    !-- Neoclassical viscosities
    !     Note: bmat is normalized by n_a m_a/tau_aa
    do i1 = 1, NSM
       fac = amas(i1) * amp * Var(NR,i1)%n * 1.d20 / ztau(i1,i1)
       xmu(NR,i1,1,1) = FSNC * bmat(i1     ,i1     ) * fac
       xmu(NR,i1,1,2) =-FSNC * bmat(i1     ,i1+NSMB) * fac
       xmu(NR,i1,2,1) =-FSNC * bmat(i1+NSMB,i1     ) * fac
       xmu(NR,i1,2,2) = FSNC * bmat(i1+NSMB,i1+NSMB) * fac
    enddo
    !-- Beam neoclassical viscosities
    if( ibeam == 0 ) then
       xmuf(NR,1:2) = 0.d0
    else
       facee = amas(1) * amp * Var(NR,1)%n * 1.d20 / ztau(1,1)
       xmuf(NR,1) = FSNCB * bmat(NSMB     ,NSMB     ) * facee ! fast ion viscosity
       xmuf(NR,2) = FSNCB * bmat(NSMB+NSMB,NSMB+NSMB) * facee ! fast ion heat viscosity
       if(NR == 0) xmuf(NR,1:2) = 0.d0
    end if

    !-- Fricrion coefficients normalized by m_a n_a : lab(NR,i,j,k,l)/(m_a n_a)
    !     Note: amat is normalized by n_a m_a/tau_aa
    do i1 = 1, NSM
       fac = 1.d0 / ztau(i1,i1)
       do i2 = 1, NSM
          lab(NR,i1,i2,1,1) = amat(i1     ,i2     ) * fac
          lab(NR,i1,i2,1,2) = amat(i1     ,i2+NSMB) * fac
          lab(NR,i1,i2,2,1) = amat(i1+NSMB,i2     ) * fac
          lab(NR,i1,i2,2,2) = amat(i1+NSMB,i2+NSMB) * fac
       end do
    end do
    !-- Beam friction coefficients
    !     laf : thermal and fast,   normalized by m_a n_a
    !     lfb : fast and thermal,   normalized by m_e n_e
    !     lff : fast and fast,      normalized by m_e n_e
    if( ibeam == 0 ) then
       laf(NR,1:NSM,1:2,1:2) = 0.d0
       lfb(NR,1:NSM,1:2,1:2) = 0.d0
       lff(NR,1:2,1:2) = 0.d0
    else
       facee = 1.d0 / ztau(1,1)
       do i1 = 1, NSM
          fac = 1.d0 / ztau(i1,i1)
          laf(NR,i1,1,1) = amat(i1     ,NSMB     ) * fac
          laf(NR,i1,1,2) = amat(i1     ,NSMB+NSMB) * fac
          laf(NR,i1,2,1) = amat(i1+NSMB,NSMB     ) * fac
          laf(NR,i1,2,2) = amat(i1+NSMB,NSMB+NSMB) * fac
          lfb(NR,i1,1,1) = amat(NSMB     ,i1     ) * facee
          lfb(NR,i1,1,2) = amat(NSMB     ,i1+NSMB) * facee
          lfb(NR,i1,2,1) = amat(NSMB+NSMB,i1     ) * facee
          lfb(NR,i1,2,2) = amat(NSMB+NSMB,i1+NSMB) * facee
       end do
       lff(NR,1,1) = amat(NSMB     ,NSMB     ) * facee
       lff(NR,1,2) = amat(NSMB     ,NSMB+NSMB) * facee
       lff(NR,2,1) = amat(NSMB+NSMB,NSMB     ) * facee
       lff(NR,2,2) = amat(NSMB+NSMB,NSMB+NSMB) * facee
    end if

    !-- Flows
    do i1 = 1, NSM
       vdiamg1 = 0.d0
       vohm1   = 0.d0
       vdiamg2 = 0.d0
       vohm2   = 0.d0
       do i2 = 1, NSM
          ! Renormalization factor
          renorm = ztau(i2,i2) / (amas(i2) * amp * Var(NR,i2)%n)
          ! Diamagnetic parallel particle flow
          vdiamg1 = vdiamg1 - (  alf(i1     ,i2     ) * BVsdiag(NR,i2,1) &
               &               - alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2))
!          if(i1==2.and.i2==2) write(6,'(I3,1P3E15.7)') nr,alf(i1     ,i2     ) * BVsdiag(NR,i2,1),- alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2),alf(i1     ,i2     ) * BVsdiag(NR,i2,1)- alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2)
          ! Ohmic parallel particle flow
          vohm1   = vohm1   -   cmat(i1     ,i2     ) * achg(i2) * Var(NR,i2)%n * renorm
          ! Diamagnetic parallel heat flow
          vdiamg2 = vdiamg2 - (- alf(i1+NSMB,i2     ) * BVsdiag(NR,i2,1) &
               &               + alf(i1+NSMB,i2+NSMB) * BVsdiag(NR,i2,2))
          ! Ohmic parallel heat flow
          vohm2   = vohm2   - (-cmat(i1+NSMB,i2     ))* achg(i2) * Var(NR,i2)%n * renorm
!          write(6,'(I3,2I2,1P3E15.7)') nr,i1,i2,amat(i1     ,i2     )-bmat(i1     ,i2     ),amat(i1     ,i2+NSMB)-bmat(i1     ,i2+NSMB),cmat(i1     ,i2     )
!          write(6,'(I3,2I2,1P3E15.7)') nr,i1,i2,amat(i1     ,i2     )-bmat(i1     ,i2     ),amat(i1     ,i2     ),-bmat(i1     ,i2     )
       end do
       vohm1 = vohm1 * aee * BEpara(NR)
       vohm2 = vohm2 * aee * BEpara(NR)

       ! Parallel particle and heat flows
       UsparNCL(NR,i1,1) = vdiamg1 + vohm1 ! particle
       UsparNCL(NR,i1,2) = vdiamg2 + vohm2 ! heat
!       if(NR<3) write(6,'(I3,I2,1P2E15.7)') nr,i1,vdiamg1,vohm1
       ! Poloidal particle flows
       UsthNCL(NR,i1,1) = ( vdiamg1 - BVsdiag(NR,i1,1) ) / bbt(NR) ! diamag component
       UsthNCL(NR,i1,2) = vohm1 / bbt(NR) ! <E.B>  component
    end do

    !-- <B.nabla.Pi> for viscous heating
    do i1 = 1, NSM
       BnablaPi(NR,i1) = sum(xmu(NR,i1,1,1:2) * ( UsparNCL(NR,i1,1:2) - BVsdiag(NR,i1,1:2) ))
    end do

    !-- Neoclassical heat diffusivities
    if( NR /= 0 ) then
       fac = fipol(NR)**2 / (bbt(NR) * (sst(NR) * sdt(NR)**2))
    else
       fac = 0.d0
    end if
    do i1 = 1, NSM
       do i2 = 1, NSM
          fac2 = fac * (    amas(i1) * amqp * Var(NR,i2)%T * rKilo &
               &        / ( achg(i1) * achg(i2) * ztau(i1,i1) ) )
          ! Classical + Pfirsch-Shluter
          chipPS(i1,i2) = (bri(NR) / fipol(NR)**2) * amat(i1+NSMB,i2     ) * fac2
          chiTPS(i1,i2) =-(bri(NR) / fipol(NR)**2) * amat(i1+NSMB,i2+NSMB) * fac2
          ! Banana-Plateau
          if( i1 /= i2 ) then
             chipBP(i1,i2) =-(  bmat(i1+NSMB,i1     ) *         alf(i1     ,i2     ) &
                  &           + bmat(i1+NSMB,i1+NSMB) *         alf(i1+NSMB,i2     ) ) * fac2
             chiTBP(i1,i2) = (  bmat(i1+NSMB,i1     ) *         alf(i1     ,i2+NSMB) &
                  &           + bmat(i1+NSMB,i1+NSMB) *         alf(i1+NSMB,i2+NSMB) ) * fac2
          else
             chipBP(i1,i2) =-(  bmat(i1+NSMB,i1     ) * (1.d0 + alf(i1     ,i2     )) &
                  &           + bmat(i1+NSMB,i1+NSMB) *         alf(i1+NSMB,i2     ) ) * fac2
             chiTBP(i1,i2) = (  bmat(i1+NSMB,i1     ) *         alf(i1     ,i2+NSMB) &
                  &           + bmat(i1+NSMB,i1+NSMB) * (1.d0 + alf(i1+NSMB,i2+NSMB))) * fac2
          end if
       end do
    end do
    ChiNCpel = ChiNC * (chipBP(1,1) + chipPS(1,1))
    ChiNCtel = ChiNC * (chiTBP(1,1) + chiTPS(1,1))
    ChiNCpil = ChiNC * (chipBP(2,2) + chipPS(2,2))
    ChiNCtil = ChiNC * (chiTBP(2,2) + chiTPS(2,2))
!    write(6,*) NR,ChiNCpil,ChiNCtil
!    write(6,*) NR,chiTBP(2,2),chiTPS(2,2)

    deallocate(PNsVL,PTsVL,amasL,achgL,sqzfac,xsta,ztau,coebsc,coebdc,amat,bmat,cmat,alf)
    deallocate(chipBP,chiTBP,chipPS,chiTPS)

  end subroutine tx_matrix_inversion

!=======================================================================
  subroutine booth9(coebsc,coencc,coebdc,amat,bmat,cmat,alf &
       &           ,den,tem,mas,zz,tpfneo,sqzfac,xsta,ztau &
       &           ,ispc,ibeam,imodel,fbeam,eps,omgtau)
!=======================================================================
!
!     neoclassical coefficients for bootstrap current etc.
!
!                        based on modified Hirshman & Sigmar formulation
!
!       *   solve matrix equations of
!                     parallel force balance
!
!       *   use the integral formula of viscosity tensor
!
!       *   friction coeffs and viscosities normalized by (na*ma/tauaa)
!
!   << OUTPUT >>
!     2000/07/07 cl31(i),cl32(i) -> coebsc(i,1),coebsc(i,2)
!    (coebsc(i,1&2) (i=1,ispc+ibeam) : bootstrap current coeff)
!    (coencc  : neoclassical enhancement factor for resistivity)
!    (coebdc  : shielding factor of beam driven current)
!     amat : friction coefficients
!     bmat : neoclassical viscosities
!     cmat : (L - M)^{-1}_ab
!     alf  : [(L - M)^{-1}M]_ab
!
!   << INPUT >>
!     den(i)   (i=1,ispc+ibeam)  :  density      den(1) ele. [arb.]
!     tem(i)   (i=1,ispc+ibeam)  :  temperature  tem(1) ele. [arb.]
!     mas(i)   (i=1,ispc+ibeam)  :  mass no.     mas(1) mele./mploton
!     zz(i)    (i=1,ispc+ibeam)  :  charge no.   zz(1)  -1.d0
!     tpfneo                     :  trapped particle fraction
!     sqzfac(i)(i=1,ispc+ibeam)  :  orbit squeezing factor
!     xsta(i)  (i=1,ispc+ibeam)  :  effective selfcollision frequency
!     ztau(i,i) (i=1,ispc)       :  collision time
!     ispc   : no. of bulk particle species with maxwellian distribution
!     ibeam  : no. of beam particle species
!     imodel(1:5) :  model for beam friction/viscosity
!     fbeam   :  factor for beam viscosity
!                       { fbeam=0 -->  full viscosity  }
!
!   << optional INPUT >>
!     eps                        :  inverse aspect ratio
!        for booth9 only
!     omgtau(i) (i=1,ispc+ibeam) :  transit frequency * tau_ii
!        for boothx only
!
!=======================================================================
  use tx_commons, only : cnpi => Pi
  integer(4), intent(in) :: ispc, ibeam, imodel(:)
  real(8), intent(in)  :: den(:), tem(:), mas(:), zz(:) &
       &                , tpfneo, xsta(:), ztau(:,:), fbeam &
       &                , sqzfac(:)
  real(8), intent(in), optional :: eps, omgtau(:)
  real(8), intent(out) :: coebsc(:,:), coencc, coebdc(:) &
       &                , amat(:,:), bmat(:,:) &
       &                , cmat(:,:), alf(:,:)
!:: local
  real(8), allocatable :: lab(:,:,:,:), mab(:,:,:,:), nab(:,:,:,:), vt(:), xmu(:,:,:)
  real(8) :: dx, eb, ec, err, fac, fdp, fp, gfun &
       &   , sint, sq3, tbdc, teff, teff0 &
       &   , vb, vb2, vb3, vbc3, vc, vc3, vxmax &
       &   , wc3, x, x0, x32, xa, xa2, xa3, xab2, xab4, xb, xc0 &
       &   , xfac, xtab, xgt, xke11, xke12, xke22, xl &
       &   , xp321, xx, yk0, yk1, yk2, z1, z2, z3 &
       &   , z30, zeff, zk, zmud, zmut, zp43, y, zkdenom
  real(8) :: xsign = -1.d0, void = 0.d0
  integer :: i, i1, i2, ib, ill, isp, isp2, ivmax &
       &   , j, j1, j2, k
!=======================================================================
    vxmax = 4.d0 ! org. 3.d0
    ivmax = 15   ! org. 100
!<<total no. of particle species and matrix size>>
    isp=ispc+ibeam
    allocate(lab(isp,isp,2,2), mab(isp,isp,2,2), nab(isp,isp,2,2), vt(isp), xmu(isp,2,2))
    isp2=2*isp
!<<ratio of trapped and passing particles>>
    xgt=tpfneo/(1.d0-tpfneo)
!<<thermal velocity of bulk particles>>
    do i=1,isp
       vt(i)=sqrt(2.0d0*tem(i)/mas(i))
    enddo
    zp43=0.75d0*sqrt(cnpi)
    sq3=sqrt(3.d0)
!=======================================================================
!  friction coefficients for bulk particles
!=======================================================================
    do i1=1,ispc
       do i2=1,ispc
          xab2=(vt(i2)/vt(i1))**2
          xab4=xab2**2
          mab(i1,i2,1,1)=-(1.d0+mas(i1)/mas(i2))/(1.d0+xab2)**1.5d0
          mab(i1,i2,1,2)=-1.5d0*(1.d0+mas(i1)/mas(i2)) &
               &        /(1.d0+xab2)**2.5d0
          mab(i1,i2,2,1)=mab(i1,i2,1,2)
          mab(i1,i2,2,2)=-(13.d0/4.d0+4.d0*xab2+15.d0/2.d0*xab4) &
               &        /(1.d0+xab2)**2.5d0
          nab(i1,i2,1,1)=-mab(i1,i2,1,1)
          nab(i1,i2,2,1)=-mab(i1,i2,1,2)
          nab(i1,i2,2,2)=27.d0/4.d0*tem(i1)/tem(i2)*xab2 &
               &        /(1.d0+xab2)**2.5d0
       enddo
    enddo
    do i1=1,ispc
       do i2=1,ispc
          nab(i1,i2,1,2)=tem(i1)*vt(i1)/tem(i2)/vt(i2)*nab(i2,i1,2,1)
       enddo
    enddo
!-----------------------------------------------------------------------
!     friction coef. lab normalized by (na*ma/tauaa)
!-----------------------------------------------------------------------
    lab(:,:,:,:) = 0.0_8
    do i1=1,ispc
       do i2=1,ispc
!          fac=(zz(i2)/zz(i1))**2*den(i2)/den(i1)
          fac=ztau(i1,i1)/ztau(i1,i2)
          do j1=1,2
             do j2=1,2
                xl=fac*nab(i1,i2,j1,j2)
                if(i1.eq.i2)then
                   do k=1,ispc
!                      xl=xl+(zz(k)/zz(i1))**2*den(k)/den(i1)*mab(i1,k,j1,j2)
                      xl=xl+ztau(i1,i1)/ztau(i1,k)*mab(i1,k,j1,j2)
                   enddo
                endif
                lab(i1,i2,j1,j2)=xl
!                if(i1==1.and.i2==1.and.j1==1.and.j2==1) print *,'xl = ',xl
             enddo
          enddo
       enddo
    enddo
!=======================================================================
!   viscosity coeff for bulk particles
!=======================================================================
    do i1=1,ispc
       xke11=0.d0
       xke12=0.d0
       xke22=0.d0
       dx=vxmax/ivmax
       do i=1,ivmax
          xa=dx*i-0.5d0*dx
          xa3=xa**3
          xa2=xa**2
          zmud=0.d0
          zmut=0.d0
          do i2=1,ispc
             xb=xa*vt(i1)/vt(i2)
             if(xb < 5.922d0) then ! above this value, fp is always 1 in double precision.
                fp=erf(xb)
                fdp=2.d0/sqrt(cnpi)*exp(-xb**2) ! derivative of erf
             else
                fp=1.d0
                fdp=0.d0
             endif
             gfun=(fp-xb*fdp)/(2.d0*xb**2) ! Chandrasekhar function
!             xtab=(zz(i2)/zz(i1))**2*den(i2)/den(i1)
             xtab=ztau(i1,i1)/ztau(i1,i2)
             zmud=zmud+xtab*(fp-gfun)
             zmut=zmut+xtab*((fp-3.d0*gfun)/xa3 &
                  &   +4.d0*(tem(i1)/tem(i2)+(vt(i1)/vt(i2))**2)*gfun/xa)
          enddo
          zmud=zmud/xa3 ! nu_D^a * tau_aa w/o 3*sqrt(pi)/4
          if(xsta(i1) /= 0.d0) then
             zkdenom=1.d0/(1.d0+2.48d0*zp43*xsta(i1)*zmud/xa)
          else
             zkdenom=0.d0
          end if
          if(present(eps)) then ! booth9
             zk=zmud*zkdenom
             zk=zk/(1.d0+15.d0*cnpi**1.5d0/32.d0*zmut*eps**1.5d0*xsta(i1)/xa)
             if(imodel(5).eq.1) zk=zmud*zkdenom
          else ! boothx
             zk=zmud/(1.d0+xsta(i1)*zp43*zmud/xa)
             zk=zk/(1.d0+15.d0*cnpi**1.5d0/32.d0*zmut/(xa*omgtau(i1)))
             if(imodel(5).eq.1) zk=zmud/(1.d0+xsta(i1)*zp43*zmud/xa)
          endif
          yk0=zk*exp(-xa2)*xa2*xa2
          yk1=yk0*xa2
          yk2=yk1*xa2
          xke11=xke11+yk0
          xke12=xke12+yk1
          xke22=xke22+yk2
       enddo
       xfac=xgt*sqzfac(i1)**(-1.5d0)
       xke11=2.d0*xke11*dx*xfac ! 2 = (8/3/sqrt(pi)) * (3*sqrt(pi)/4)
       xke12=2.d0*xke12*dx*xfac
       xke22=2.d0*xke22*dx*xfac
       xmu(i1,1,1)=xke11
       xmu(i1,1,2)=2.5d0*xke11-xke12
       xmu(i1,2,1)=xmu(i1,1,2)
       xmu(i1,2,2)=25.d0/4.d0*xke11-5.d0*xke12+xke22
    enddo
    do i1 = ispc+1, isp ! added 6 lines by kamata 2002/10/21
       xmu(i1,1,1) = 0.0d0
       xmu(i1,1,2) = 0.0d0
       xmu(i1,2,1) = 0.0d0
       xmu(i1,2,2) = 0.0d0
    enddo
!=======================================================================
!  fast ion contribution
!
!     The model of a set for beam friction forces is taken from
!        J.P. Wang, et al., JAERI-M 92-107
!        J.P. Wang, et al., Nucl. Fusion 34 (1994) 231.
!     The model of a minimum set for beam friction forces is taken from
!        S.P. Hirshman and D.J. Sigmar, Nucl. Fusion 21 (1981) p.1174
!=======================================================================
    if(ibeam.gt.0)then
       do ib=ispc+1,isp
          z1=0.d0
          z2=0.d0
          do i=2,ispc
             z1=z1+zz(i)**2*den(i)
             z2=z2+zz(i)**2*den(i)/mas(i)
          enddo
          z1=z1/den(1)
          z2=z2*mas(ib)/den(1)
          wc3=zp43*vt(1)**3*mas(1)/mas(ib)*z1
          vc3=wc3*z2/z1
          vc=vc3**(1.d0/3.d0) ! critical speed
          ec=mas(ib)*vc**2*0.5d0
!--effective beam energy eb <-- nb, tb=pb/nb
          sint=0.d0
          x=3.d0*tem(2)/ec
          x32=x**1.5d0
          xp321=x32+1.d0
          z3=x32/(1.d0+x32)
          dx=0.1d0
          teff=0.d0
          do while( teff < tem(ib) )
             x0=x
             teff0=teff
             x=x+dx
             eb=x*ec
             z30=z3
             x32=x**1.5d0
             z3=x32/(1.d0+x32)
             sint=sint+0.5d0*(z3+z30)*dx
             teff=ec*sint/log((x32+1.d0)/xp321)
          end do
          eb=ec*(x0+dx*(tem(ib)-teff0)/(teff-teff0))
          vb=sqrt(2.d0*eb/mas(ib))
          vb2=vb*vb
          vb3=vb2*vb
!-----
          vbc3=vb3+vc3
!          xvb3=vb3/vbc3/log(vbc3/vc3)
!<e-b,1,1>
!   Since the equation for beam heat flux is not solved, then
!          lab( *,ib,*,2) = void,
!          lab(ib, *,2,*) = void.
          lab( 1,ib,1,1)=zz(ib)**2*den(ib)/den(1)
          if(imodel(3).eq.1) then ! minimum set
             if(imodel(2).eq.1) lab( 1,ib,2,1)=1.5d0*lab(1,ib,1,1)
             lab(ib,ib,1,1)=1.d0
          else
             lab( 1, 1,1,1)=lab(1,1,1,1)-lab(1,ib,1,1)
             lab(ib, 1,1,1)=lab(1,ib,1,1)
             lab(ib,ib,1,1)=-lab(ib,1,1,1)
!<e-b,1,2>
             if(imodel(3).eq.2) then ! original
!_org                lab( 1,ib,1,2)=1.5d0*lab(1,ib,1,1)
!_org                lab( 1, 1,1,2)=lab( 1, 1,1,2)-lab(1,ib,1,2)
!_org                lab(ib, 1,1,2)=lab(ib, 1,1,2)-lab(1,ib,1,2)
!_org                lab(ib,ib,1,2)=lab(ib,ib,1,2)+lab(1,ib,1,2)
             else ! maximum set
                lab( 1,ib,1,2)=void
                lab( 1, 1,1,2)=lab( 1, 1,1,2)-1.5d0*lab( 1,ib,1,1)
                lab(ib, 1,1,2)=1.5d0*lab(1,ib,1,1)
                lab(ib,ib,1,2)=void
             endif
!<e-b,2,1> fast ion-electron thermal friction
             if(imodel(2).eq.1)then
                lab( 1,ib,2,1)=1.5d0*lab(1,ib,1,1)
                lab( 1, 1,2,1)=lab(1,1,2,1)-lab(1,ib,2,1)
                if(imodel(3).eq.2) then
!_org                   lab(ib, 1,2,1)=lab(ib, 1,2,1)-lab(1,ib,2,1)
!_org                   lab(ib,ib,2,1)=lab(ib,ib,2,1)+lab(1,ib,2,1)
                else
                   lab(ib, 1,2,1)=void
                   lab(ib,ib,2,1)=void
!<e-b,2,2>
                   lab( 1,ib,2,2)=void
                   lab( 1, 1,2,2)=lab(1,1,2,2)-3.25d0*lab( 1,ib,1,1)
                   lab(ib, 1,2,2)=void
                   lab(ib,ib,2,2)=void
                endif
             endif
!<i-b,1,1>
             x=(log(vbc3)-3.d0*log(vb+vc) &
                  &        +2.d0*sq3*atan((2.d0*vb-vc)/sq3/vc)+cnpi/sq3)/(6.d0*vc)
             do i=2,ispc
                y=(zz(ib)/zz(i))**2*den(ib)/den(i)*zp43/log(vbc3/vc3)
                lab( i,ib,1,1)=y*(1.d0+mas(i)/mas(ib))*vt(i)**3/vc3
                lab( i, i,1,1)=lab(i, i,1,1)-lab(i,ib,1,1)
                lab(ib, i,1,1)=lab(i,ib,1,1) & 
                     &        *(den(i)*zz(i)**2/den(1))**2 &
                     &        *sqrt(mas(i)/mas(1))*(tem(1)/tem(i))**1.5d0
                lab(ib,ib,1,1)=lab(ib,ib,1,1)-lab(ib,i,1,1)
!<i-b,1,2>
                if(imodel(3).eq.0) then
                   lab( i,ib,1,2)=void
                   lab( i, i,1,2)=lab( i, i,1,2)
                   lab(ib, i,1,2)=0.d0
                   lab(ib,ib,1,2)=void
!<i-b,2,1>
                   if(imodel(2).eq.1)then
                      lab( i,ib,2,1)=y*5.d0*vt(i)*x
!                      lab( i,ib,2,1)=y*vt(i)*x
                      lab( i, i,2,1)=lab( i, i,2,1)-lab( i,ib,2,1)
                      lab(ib, i,2,1)=void
                      lab(ib,ib,2,1)=void
!<i-b,2,2>
                      lab( i,ib,2,2)=void
                      lab( i, i,2,2)=lab( i, i,2,2)
                      lab(ib, i,2,2)=void
                      lab(ib,ib,2,2)=void
                   endif
                endif
             enddo
          endif
! Fast ion viscosity
          if(imodel(1).eq.1)then
             xmu(ib,1,1)=xgt*(zz(ib)/zz(1))**2*den(ib)/den(1) &
!             xmu(ib,1,1)=xgt*ztau(1,1)/ztau(1,ib) &
                  &     *(1.d0+(1.d0+z1/z2)*vc3/(vb2/x-2.d0*vc3))
          else
             xmu(ib,1,1)=xgt*(zz(ib)/zz(1))**2*den(ib)/den(1) &
!             xmu(ib,1,1)=xgt*ztau(1,1)/ztau(1,ib) &
                  &     *z1/z2*vc3/(vb2/x-2.d0*vc3)
          endif
!>>
          xmu(ib,1,1)=(1.d0-fbeam)*xmu(ib,1,1)*sqzfac(ib)**(-1.5d0)
!>>
!<<dummy set for parallel thermal friction>>
          lab(ib,ib,2,2)=1.d0
       enddo
    endif
!=======================================================================
!<amat> : friction coefficient matrix, L
    do i1=1,isp
       do i2=1,isp
          amat(i1    ,i2    )=lab(i1,i2,1,1)
          amat(i1    ,i2+isp)=lab(i1,i2,1,2)
          amat(i1+isp,i2    )=lab(i1,i2,2,1)
          amat(i1+isp,i2+isp)=lab(i1,i2,2,2)
       enddo
    enddo
!-----
    do i1=1,isp2
       do i2=1,isp2
          bmat(i1,i2)=amat(i1,i2)
       enddo
    enddo
    do i1=1,isp
       bmat(i1    ,i1    )=bmat(i1    ,i1    )-xmu(i1,1,1)
       bmat(i1    ,i1+isp)=bmat(i1    ,i1+isp)-xmu(i1,1,2)
       bmat(i1+isp,i1    )=bmat(i1+isp,i1    )-xmu(i1,2,1)
       bmat(i1+isp,i1+isp)=bmat(i1+isp,i1+isp)-xmu(i1,2,2)
    enddo
!=======================================================================
!<cmat> : (L - M)^{-1}
    call matslv(isp2,isp2,bmat,cmat,err,ill)
!=======================================================================
!<bmat> : viscosity matrix, M
    do i1=1,isp2
       do i2=1,isp2
          bmat(i1,i2)=0.d0
       enddo
    enddo
    do i1=1,isp
       bmat(i1    ,i1    )=xmu(i1,1,1)
       bmat(i1    ,i1+isp)=xmu(i1,1,2)
       bmat(i1+isp,i1    )=xmu(i1,2,1)
       bmat(i1+isp,i1+isp)=xmu(i1,2,2)
    enddo
!=======================================================================
!<alf> : (L - M)^{-1}M
    do i1=1,isp2
       do i2=1,isp2
          x=0.d0
          do k=1,isp2
             x=x+cmat(i1,k)*bmat(k,i2)
          enddo
          alf(i1,i2)=x
       enddo
    enddo
!=======================================================================
!<L_31,L_32>
    do i=1,isp
       j=i+isp
       coebsc(i,1)=-alf(1,i)
       coebsc(i,2)=-alf(1,j)
       do i1=2,isp
          coebsc(i,1)=coebsc(i,1)+zz(i1)*den(i1)*alf(i1,i)/den(1)
          coebsc(i,2)=coebsc(i,2)+zz(i1)*den(i1)*alf(i1,j)/den(1)
       enddo
    enddo
!-----
    coebsc(1,2)=-coebsc(1,2)
    do i=2,isp
       coebsc(i,1)=-coebsc(i,1)
    enddo
!=======================================================================
!     neoclassical conductivity
!=======================================================================
!--matrix method
    zeff=0.d0
    do i=2,isp
       zeff=zeff+zz(i)**2*den(i)
    enddo
    zeff=zeff/den(1)
    xc0=0.d0
    do i1=1,isp
       j1=i1
       ! beam ion is normalized by m_e n_e/tau_ee
       if(ibeam.eq.1.and.i1.eq.isp)j1=1
       xx=0.d0
       ! cmat(a,b) : (L - M)^{-1}_ab
       do i2=1,isp
          xx=xx+zz(i2)*den(i2)/den(1)*(-cmat(i2,i1))
       enddo
       ! renormalized by multiplying (m_e n_e/tau_ee) / (m_i n_i/tau_ii)
!       xc0=xc0+xx*sqrt(mas(1)/mas(j1))*(tem(j1)/tem(1))**1.5d0 &
!            &    *(den(1)/(zz(j1)**2*den(j1)))**2 &
       xc0=xc0+xx*(mas(1)*den(1)*ztau(j1,j1)/(mas(j1)*den(j1)*ztau(1,1))) &
            &    *zz(i1)*den(i1)/den(1)
    enddo
    coencc=xc0
!=======================================================================
!   coefficients for nbi current drive
!=======================================================================
    coebdc(:)=0.d0
    if(ibeam.gt.0)then
       tbdc=0.d0
       do i=1,isp
          coebdc(i)=zz(i)*den(i)/den(1)*cmat(i,isp)
          tbdc=tbdc+coebdc(i)
       enddo
       coebdc(isp+1)=tbdc
    endif
!
!     Analytical model, G
!
    xsign = -1.d0
    coebdc(isp+2)=((xmu(1,2,2) + sqrt(2.d0) &
         &     + 3.25d0 * zeff) * xmu(1,1,1) + ( -xsign &
         &     * xmu(1,1,2) + 1.5d0 * zeff) * xsign * xmu(1,1,2)) /  &
         &     ((xmu(1,2,2) + sqrt(2.d0) &
         &     + 3.25d0 * zeff) * (xmu(1,1,1) + zeff) &
         &     - (-xsign * xmu(1,1,2) + 1.5d0 * zeff)**2)
!=======================================================================

    deallocate(lab,mab,nab,vt,xmu)

  end subroutine booth9

!===================================================================
!
!  Inverting a matrix
!
!     solve b*x=c to get x or invert b
!
!===================================================================
  subroutine matslv(m,n,b,c,err,ill)
    integer(4), intent(in) :: m, n
    integer(4), intent(out) :: ill ! do nothing
    real(8), dimension(:,:), intent(inout) :: b, c
    real(8), intent(out) :: err ! do nothing
    integer(4) :: i, j, k, l, nn, nnn
    real(8) :: x
!===================================================================
!     lu decomposition
!===================================================================
    do j=1,m
!-------------------------------------------------------------------
!::b(j,k) : k <= j
!-------------------------------------------------------------------
       do k=1,j
          x=0.d0
          if(k > 1) then
             do l=1,k-1
                x=x+b(j,l)*b(l,k)
             end do
          endif
!-----
          b(j,k)=b(j,k)-x
       end do
!-------------------------------------------------------------------
!::b(j,k) : k > j
!-------------------------------------------------------------------
       if(j < m) then
          do k=j+1,m
             x=0.d0
!-----
             do l=1,j-1
                x=x+b(j,l)*b(l,k)
             end do
!-----
             if(b(j,j) == 0.d0) b(j,j)=epsilon(1.d0)
!-----
             b(j,k)=(b(j,k)-x)/b(j,j)
          enddo
       endif
!-------------------------------------------------------------------
    enddo
!===================================================================
    if(n > m) then
!===================================================================
!  solve simultaneous equation
!===================================================================
       nn=n-m
       do nnn=1,nn
          do j=1,m
             x=0.d0
!-----
             if(j > 1) then
                do k=1,j-1
                   x=x+b(j,k)*b(k,m+nnn)
                end do
             endif
!-----
             b(j,m+nnn)=(b(j,m+nnn)-x)/b(j,j)
          end do
!-------------------------------------------------------------------
          do j=m,1,-1
             x=0.d0
!-----
             if(j < m) then
                do k=j+1,m
                   x=x+b(j,k)*b(k,m+nnn)
                end do
             endif
!-----
             b(j,m+nnn)=b(j,m+nnn)-x
          end do
       end do
!-----
       ill=0
       err=0.d0
!===================================================================
!  inverse matrix : c(i,j)
!===================================================================
    else
       do i=1,m
          do j=1,m
             c(i,j)=0.d0
          end do
          c(i,i)=1.d0
       end do
!-----
       do nnn=1,m
          do j=1,m
             x=0.d0
!-----
             if(j > 1) then
                do k=1,j-1
                   x=x+b(j,k)*c(k,nnn)
                end do
             endif
!-----
             c(j,nnn)=(c(j,nnn)-x)/b(j,j)
          end do
!-------------------------------------------------------------------
          do j=m,1,-1
             x=0.d0
!-----
             if(j < m) then
                do k=j+1,m
                   x=x+b(j,k)*c(k,nnn)
                end do
             endif
!-----
             c(j,nnn)=c(j,nnn)-x
          end do
       end do
!===================================================================
       ill=0
       err=0.d0
    end if

  end subroutine matslv

end module matrix_inversion
