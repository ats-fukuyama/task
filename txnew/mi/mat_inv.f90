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
       &                         ddPhidpsi_in,MDLNEOL)
    use tx_commons
    use tx_interface, only : coll_freq, ftfunc, corr
    integer(4), intent(in) :: NR, MDLNEOL
    real(8), intent(in) :: ddPhidpsi_in
    real(8), intent(out) :: ETAout, BJBSout, ChiNCpel, ChiNCtel, ChiNCpil, ChiNCtil
    integer(4) :: imodel(10), ibeam, NSMB, i, j, i1, i2, idebug = 0, MDLNEOLgflux, icoebdc
    real(8) :: fbeam, sqzfaccoef, coencc, bjsum, ftl, fac, facee, epsL, &
         &     fac2, fac3, gflxbp, gflxps, gflxware, coefmneo
    real(8) :: vdiamg1, vdiamg2, renorm, rbanana
    real(8) :: smallvalue = 1.d-4
    real(8), dimension(:), allocatable :: PNsVL, PTsVL, amasL, achgL, sqzfac, xsta, coebdc
    real(8), dimension(:,:), allocatable :: ztau, coebsc, amat, bmat, cmat, dmat, alf, &
         & chipBP, chiTBP, chipPS, chiTPS, DpBP, DTBP, DpPS, DTPS, vohm

    imodel(1) = 0 ! Fast ion viscosity
    imodel(2) = 1 ! Required for NBCD
    imodel(3) = 2 ! Fast ion contribution to friction forces
    imodel(4) = 0 ! Unused
    imodel(5) = 0 ! PS contribution nil when 1
    imodel(6) = 3 ! Higher-order flow contribution, valid only for nccoe

    fbeam = 0.d0 ! 0: fast ion viscosity computed, 1: nil

    if(PNBH == 0.D0 .and. PNbV(NR) < 1.D-8) then
       ibeam = 0
       NSMB  = NSM
    else
       ibeam = 1
       NSMB  = NSM + 1
    end if

    icoebdc = NSM + 1 + 2 ! thermal species + beam species + 2 model

    allocate(PNsVL(NSMB),PTsVL(NSMB),amasL(NSMB),achgL(NSMB),sqzfac(NSMB),xsta(NSMB))
    allocate(ztau(NSM,NSM))
    allocate(coebsc(NSMB,2),coebdc(icoebdc),amat(2*NSMB,2*NSMB),bmat(2*NSMB,2*NSMB), &
         &   cmat(2*NSMB,2*NSMB),dmat(2*NSMB,2*NSMB),alf(2*NSMB,2*NSMB),vohm(NSM,2))
    allocate(chipBP(NSM,NSM),chiTBP(NSM,NSM),chipPS(NSM,NSM),chiTPS(NSM,NSM))
    allocate(DpBP(NSM,NSM),DTBP(NSM,NSM),DpPS(NSM,NSM),DTPS(NSM,NSM))

    amasL(1:NSM) = amas(1:NSM)
    achgL(1:NSM) = achg(1:NSM)
    PNsVL(1:NSM) = Var(NR,1:NSM)%n
    PTsVL(1:NSM) = Var(NR,1:NSM)%T * rKilo

    if( NR /= 0 ) then
       epsL = epst(NR)
       ftl = ft(NR)
    else
       ! --- Artificial potate orbit effect, i.e. finite f_t even at axis ---
       !     It yields the finite banana viscosity at axis.
       epsL = smallvalue
       ftl = ftfunc(epsL)
       ! --------------------------------------------------------------------
    end if

    !-- Orbit squeezing factor, whose definition is identical to that in NCLASS
    !     Note that sqzfac would become infinity at the magnetic axis due to d/dpsi(dPhi/dpsi).
    !               sqzfaccoef = 0 when NR = 0 because ddPhidpsi_in is set to zero.
    sqzfaccoef = (fipol(NR)/bb)**2*amqp*abs(ddPhidpsi_in)
    do i = 1, NSM
       fac = 1.d0
       if(NR /= 0) then
          rbanana = sqrt(epsL)*amas(i)*amqp/(abs(achg(i))*BthV(NR)) &
               &   *sqrt(2.d0*Var(NR,i)%T*rKeV/(amas(i)*amp)) ! banana width
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
       PTsVL(NSMB)  = PTbV(NR) * rKilo
       fac = 1.d0
       if(NR /= 0) then
          rbanana = sqrt(epsL)*amb*amqp/(abs(achgb)*BthV(NR)) &
               &   *sqrt(2.d0*PTsVL(NSMB)*aee/(amb*amp)) ! banana width
          if(r(NR) < rbanana) fac = 0.d0
       end if
       sqzfac(NSMB) = 1.d0 + fac * sqzfaccoef*amb/abs(achgb)
       xsta(NSMB)   = 0.d0
    end if

    if( MDLNEOL == 1 ) then
       !-- Call Matrix Inversion (booth9)
       call booth9(coebsc,coencc,coebdc,amat,bmat,cmat,alf &
            &     ,PNsVL,PTsVL,amasL,achgL,ftl,sqzfac,xsta,ztau &
            &     ,NSM,ibeam,imodel,fbeam,epsL)
       MDLNEOLgflux = 1
    else ! if( MDLNEOL == 2 ) then
       !-- Call Matrix Inversion (nccoe)
       call nccoe (coebsc,coencc,coebdc,amat,bmat,cmat,dmat,alf &
            &     ,PNsVL,PTsVL,amasL,achgL,mxneo(NR),fmneo(1,NR),gamneo(NR),ftl,sqzfac,ztau &
            &     ,NSM,ibeam,imodel,fbeam,idebug)
       MDLNEOLgflux = 2
    end if

    !-- Neoclassical resistivity
    !     Note: coencc is normalized by n_e m_e/tau_ee
    i = 1
    if( coencc /= 0.d0 ) then
       ETAout = amas(i) * amp / (coencc * aee**2 * Var(NR,i)%n * 1.d20 * ztau(i,i))
    else
       ! coencc sometimes becomes nought at axis.
       ! In such a case, eta is replaced by the Spitzer resistivity.
       ETAout = corr(1.D0) * amas(1) * amqp * rNuei(NR) / (Var(NR,1)%n * 1.d20 * AEE)
    end if

    !-- Bootstrap current
    bjsum = sum(Var(NR,1:NSM)%T * rKeV / abs(achg(1:NSM)) &
            & * (  coebsc(1:NSM,1) * dPsdpsi(NR,1:NSM) / Var(NR,1:NSM)%p &
            &    + coebsc(1:NSM,2) * dTsdpsi(NR,1:NSM) / Var(NR,1:NSM)%T))
    BJBSout = - fipol(NR) * Var(NR,1)%n * 1.d20 * bjsum
!    write(6,'(I3,1P6E15.7)') nr,coebsc(1:NSM,1),dPsdpsi(NR,1:NSM),Var(NR,1:NSM)%p

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
       vohm(i1,:) = 0.d0
       vdiamg2 = 0.d0
       do i2 = 1, NSM
          ! Renormalization factor
          renorm = ztau(i2,i2) / (amas(i2) * amp * Var(NR,i2)%n)
          ! Diamagnetic parallel particle flow
          vdiamg1 = vdiamg1 - (  alf(i1     ,i2     ) * BVsdiag(NR,i2,1) &
               &               - alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2))
!          if(i1==2.and.i2==2) write(6,'(I3,1P3E15.7)') nr,alf(i1     ,i2     ) * BVsdiag(NR,i2,1),- alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2),alf(i1     ,i2     ) * BVsdiag(NR,i2,1)- alf(i1     ,i2+NSMB) * BVsdiag(NR,i2,2)
          ! Ohmic parallel particle flow
          vohm(i1,1) = vohm(i1,1) -   cmat(i1     ,i2     ) * achg(i2) * Var(NR,i2)%n * renorm
          ! Diamagnetic parallel heat flow
          vdiamg2 = vdiamg2 - (- alf(i1+NSMB,i2     ) * BVsdiag(NR,i2,1) &
               &               + alf(i1+NSMB,i2+NSMB) * BVsdiag(NR,i2,2))
          ! Ohmic parallel heat flow
          vohm(i1,2) = vohm(i1,2) - (-cmat(i1+NSMB,i2     ))* achg(i2) * Var(NR,i2)%n * renorm
!          write(6,'(I3,2I2,1P3E15.7)') nr,i1,i2,amat(i1     ,i2     )-bmat(i1     ,i2     ),amat(i1     ,i2+NSMB)-bmat(i1     ,i2+NSMB),cmat(i1     ,i2     )
!          write(6,'(I3,2I2,1P3E15.7)') nr,i1,i2,amat(i1     ,i2     )-bmat(i1     ,i2     ),amat(i1     ,i2     ),-bmat(i1     ,i2     )
       end do
!!$       end do
       vohm(i1,:) = vohm(i1,:) * aee * BEpara(NR)

       ! Parallel particle and heat flows
       UsparNCL(NR,i1,1) = vdiamg1 + vohm(i1,1) ! particle
       UsparNCL(NR,i1,2) = vdiamg2 + vohm(i1,2) ! heat
       ! Poloidal particle flows
       UsthNCL(NR,i1,1) = ( vdiamg1 - BVsdiag(NR,i1,1) ) / bbt(NR) ! diamag component
       UsthNCL(NR,i1,2) = vohm(i1,1) / bbt(NR) ! <E.B>  component
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

    !-- Neoclassical particle diffusivities
    do i1 = 1, NSM
       do i2 = 1, NSM
          fac2 = fac * (    amas(i1) * amqp * Var(NR,i2)%T * rKilo &
               &        / ( achg(i1) * achg(i2) * ztau(i1,i1) ) )
          ! Classical + Pfirsch-Shluter
          DpPS(i1,i2) =-(bri(NR) / fipol(NR)**2) * amat(i1     ,i2     ) * fac2
          DTPS(i1,i2) = (bri(NR) / fipol(NR)**2) * amat(i1     ,i2+NSMB) * fac2
          ! Banana-Plateau
          if( i1 /= i2 ) then
             DpBP(i1,i2) = (  bmat(i1     ,i1     ) *         alf(i1     ,i2     ) &
                  &         + bmat(i1     ,i1+NSMB) *         alf(i1+NSMB,i2     ) ) * fac2
             DTBP(i1,i2) =-(  bmat(i1     ,i1     ) *         alf(i1     ,i2+NSMB) &
                  &         + bmat(i1     ,i1+NSMB) *         alf(i1+NSMB,i2+NSMB) ) * fac2
          else
             DpBP(i1,i2) = (  bmat(i1     ,i1     ) * (1.d0 + alf(i1     ,i2     )) &
                  &         + bmat(i1     ,i1+NSMB) *         alf(i1+NSMB,i2     ) ) * fac2
             DTBP(i1,i2) =-(  bmat(i1     ,i1     ) *         alf(i1     ,i2+NSMB) &
                  &         + bmat(i1     ,i1+NSMB) * (1.d0 + alf(i1+NSMB,i2+NSMB))) * fac2
          end if
       end do
    end do

    !-- Neoclassical particle flux: gflux = <Gamma . nabla psi>
    fac2 = fac / fipol(NR)
    fac3 = fac2 * BEpara(NR)
    do i1 = 1, NSM
       gflxbp   = 0.d0
       gflxps   = 0.d0
       gflxware = 0.d0
       do i2 = 1, NSM
          gflxbp   = gflxbp   + DpBP(i1,i2) * dPsdpsi(NR,i2) / Var(NR,i2)%p &
               &              + DTBP(i1,i2) * dTsdpsi(NR,i2) / Var(NR,i2)%T
          gflxps   = gflxps   + DpPS(i1,i2) * dPsdpsi(NR,i2) / Var(NR,i2)%p &
               &              + DTPS(i1,i2) * dTsdpsi(NR,i2) / Var(NR,i2)%T
          gflxware = gflxware + fac3 * (achg(i2) * Var(NR,i2)%n) / (achg(i1) * Var(NR,i1)%n) &
               &                     * (  bmat(i1     ,i1     ) * ( cmat(i2     ,i1     )) &
               &                        + bmat(i1     ,i1+NSMB) * ( cmat(i2     ,i1+NSMB)))! &
       end do
!!$       gflxware = - fac2 * (amas(i1) * amqp) / (ztau(i1,i1) * achg(i1)) &
!!$            &   * (  bmat(i1     ,i1     ) * vohm(i1,1) &
!!$            &      - bmat(i1     ,i1+NSMB) * vohm(i1,2) )

       gflux(NR,i1,MDLNEOLgflux) = Var(NR,i1)%n * 1.d20 * (sst(NR) * sdt(NR)**2) &
            &                    * ( -( gflxbp + gflxps ) + gflxware )
    end do

    deallocate(PNsVL,PTsVL,amasL,achgL,sqzfac,xsta,ztau,coebsc,coebdc,amat,bmat,cmat,dmat,alf,vohm)
    deallocate(chipBP,chiTBP,chipPS,chiTPS,DpBP,DTBP,DpPS,DTPS)

  end subroutine tx_matrix_inversion

!=======================================================================
  subroutine nccoe(coebsc,coencc,coebdc &
    &             ,amat,bmat,cmat,dmat,alf &
    &             ,den,tem,mas,zz &
    &             ,mxneo,fmneo,gmneo,tpfneo,sqzfac,ztau &
    &             ,ispc,ibeam,imodel,fbeam,idebug)
!=======================================================================
!
!     neoclassical coefficients for bootstrap current etc.
!
!                                based on Hirshman & Sigmar formulation
!                                with Shaing's viscosity formula
!
!       *   solve matrix equations of
!                     parallel force balance
!           S.P. Hirshman and D.J. Sigmar, Nucl. Fusion 21 (1981) 1079.
!
!       *   use the integral formula of viscosity tensor
!                     for banana-plateau-Pfirsch-Schluter transition
!       [a] K.C. Shaing, M. Yokoyama, M. Wakatani and C.T. Hsu
!              Phys. Plasmas 3 (1996) 965.
!
!       *   include the orbit squeezing effect in banana viscosity
!           K.C. Shaing and R.D. Hazeltine, Phys. Fluids B 4 (1992) 2547.
!           K.C. Shaing, C.T. Hsu and R.D. Hazeltine,
!              Phys. Plasmas 1 (1994) 3365.
!
!       *   include higher-order flow effects on friction forces
!           M. Honda, Phys. Plasmas 21 (2014) 092508.
!
!       *   friction coeffs and viscosities normalized by (na*ma/tauaa)
!
!   << OUTPUT >>
!     2000/07/07 cl31(i),cl32(i) -> coebsc(i,1),coebsc(i,2)
!    (coebsc(i,1&2) (i=1,ispc+ibeam) : bootstrap current coeff)
!    (coencc  : neoclassical resistivity enhancement factor; not functioning)
!    (coebdc  : stacking factor of beam driven current; not functioning)
!     amat : friction coefficients
!     bmat : neoclassical viscosities
!     cmat : (L - M)^{-1}_ab
!     dmat : L * [(L - M)^{-1}M]_ab
!     alf  : [(L - M)^{-1}M]_ab
!
!   << INPUT >>
!     den(i)   (i=1,ispc+ibeam)  :  density      den(1) ele. [arb.]
!     tem(i)   (i=1,ispc+ibeam)  :  temperature  tem(1) ele. [eV]
!     mas(i)   (i=1,ispc+ibeam)  :  mass no.     mas(1) mele./mploton
!     zz(i)    (i=1,ispc+ibeam)  :  charge no.   zz(1)  -1.d0
!     mxneo       : maximum poloidal fourier mode number for fmneo
!     fmneo(1:10) : fourier factors for PS diffusion, eq. (5) in [a]
!     gmneo       : < n grad theta >
!     tpfneo                     :  trapped particle fraction
!     sqzfac(i)(i=1,ispc+ibeam)  :  orbit squeezing factor
!     ztau(i,i) (i=1,ispc+ibeam) :  collision time
!     ispc  : no. of bulk particle species with maxwellian distribution
!     ibeam : no. of beam particle species
!     imodel(1:6)  :  model for beam friction/viscosity up to (5)
!                     model for friction coefs. for (6)
!     fbeam   :  factor for beam viscosity
!                       { fbeam=0 -->  full viscosity  }
!     idebug  :  debug option
!
!     imodel(6) = 0: approximate coll. operator + 2nd moment
!               = 1: full        coll. operator + 2nd moment
!               = 2: approximate coll. operator + 3rd moment
!               = 3: full        coll. operator + 3rd moment
!
!=======================================================================
    use tx_commons, only : cnpi => Pi, cnmp => AMP, cnec => AEE
    implicit none

!:: argument
    integer, intent(in) :: mxneo, ispc, ibeam, imodel(6), idebug
    real(8), intent(in) :: fmneo(10), gmneo, tpfneo, fbeam
    real(8), intent(in), dimension(:) :: den, tem, mas, zz, sqzfac
    real(8), intent(in), dimension(:,:) :: ztau
    real(8), intent(out) :: coebsc(:,:), coencc, coebdc(:) &
         &        , amat(:,:), bmat(:,:), cmat(:,:), dmat(:,:) &
         &        , alf(:,:)
!:: local
    real(8), allocatable :: lab(:,:,:,:), mab(:,:,:,:), nab(:,:,:,:), vt(:), vtn(:), xmu(:,:,:)
    real(8), allocatable :: alph(:,:,:), alphhat(:,:,:), emat(:,:), fmat(:,:) &
         &   , delt(:,:), sumbc2(:), zast(:), zmas(:)
    real(8), allocatable :: labtmp(:,:,:,:), labe(:,:,:,:), lab3(:,:,:)
    real(8) :: dx, eb, ec, err, fac, fdp, fm, fp, gfun &
         &   , sint, sq3, tbdc, teff &
         &   , teff0, vb, vb2, vb3, vbc3, vc, vc3  &
         &   , vxmax, wc3, x, x0, x32, xa, xa2, xa3, xab2 &
         &   , xab4, xb, xc0, xfac, xtab, xgt, xke11, xke12, xke22 &
         &   , xl, xp321, xx, yk0, yk1, yk2, z1, z2 &
         &   , z3, z30, zat, zeff, zk, zkb, zkps, zmud, zmui, zmut &
         &   , zomg, zp43, zu, zu2, zui, renorm, y &
         &   , xmu1, xmu2, xmu3, ell11, ell12, ell22, tmp &
         &   , yab, sqdnm, ytemab, denom &
         &   , xlm, xln, yt, xab, xb2
    real(8) :: xsign = -1.d0, void = 0.d0, wc13 = 1.d0 / 3.d0
    real(8) :: epsab, c2, f, alp, barK1, barK2, barK2mod &
         &   , c1, c1a, c3
    real(8) :: funp, fundp
    integer :: i, i1, i2, i3, ib, ill, is, isp, isp2, ivmax &
         &   , j, j1, j2, js, k, m

!=======================================================================
!!   Standard
!    vxmax = 4.d0 ! org. 3.d0
!    ivmax = 15   ! org. 100
!   High definition
    vxmax = 4.5d0 ! org. 3.d0
    ivmax = 20    ! org. 100
!<<total no. of particle species and matrix size>>
    isp=ispc+ibeam
    allocate(lab(isp,isp,2,2), mab(isp,isp,3,3), nab(isp,isp,3,3), vt(isp), vtn(isp), xmu(isp,2,2))
    allocate(alph(isp,isp,2), alphhat(isp,isp,2), emat(isp,isp), fmat(isp,isp), delt(isp,isp))
    isp2=2*isp
!<<ratio of trapped and passing particles>>
    xgt=tpfneo/(1.d0-tpfneo)
!<<thermal velocity of bulk particles>>
    allocate(zmas(isp), zast(isp))
    do i = 1, ispc
       zmas(i) = mas(i) * cnmp
       vt(i)   = sqrt( 2.0d0 * cnec * tem(i) / zmas(i) )
       vtn(i)  = sqrt( 2.0d0        * tem(i) /  mas(i) )
    enddo
    zp43=0.75d0*sqrt(cnpi)
    sq3=sqrt(3.d0)
!=======================================================================
!  friction coefficients for bulk particles
!=======================================================================
    if(mod(imodel(6),2)==0) then
       do i1=1,ispc
          do i2=1,ispc
             xab2=(vtn(i2)/vtn(i1))**2
             xab4=xab2**2
             yab=mas(i1)/mas(i2)
             ytemab=tem(i1)/tem(i2)
             sqdnm=1.d0/sqrt(1.d0+xab2)
             ! M_{ab}
             mab(i1,i2,1,1)=-           (1.d0+yab)/(1.d0+xab2)   *sqdnm ! ok
             mab(i1,i2,1,2)=-1.5d0     *(1.d0+yab)/(1.d0+xab2)**2*sqdnm ! sign differs from NCLASS
             mab(i1,i2,1,3)=-15.d0/8.d0*(1.d0+yab)/(1.d0+xab2)**3*sqdnm ! ok
             mab(i1,i2,2,1)=mab(i1,i2,1,2) ! ok
             mab(i1,i2,2,2)=-(13.d0/4.d0+4.d0*xab2+15.d0/2.d0*xab4) &
                  &        /(1.d0+xab2)**2*sqdnm ! ok
             mab(i1,i2,2,3)=-(69.d0/16.d0+6.d0*xab2+63.d0/4.d0*xab4) &
                  &        /(1.d0+xab2)**3*sqdnm ! sign differs
             mab(i1,i2,3,1)=mab(i1,i2,1,3) ! ok
             mab(i1,i2,3,2)=mab(i1,i2,2,3) ! ok
             mab(i1,i2,3,3)=-( 433.d0/64.d0+17.d0*xab2+459.d0/8.d0*xab4 &
                  &           +28.d0*xab2*xab4+175.d0/8.d0*xab4**2) &
                  &        /(1.d0+xab2)**4*sqdnm ! ok

!!$             ! N_{ab}: sqrt is introduced, like NCLASS
!!$             nab(i1,i2,1,1)=-mab(i1,i2,1,1) ! ok
!!$             nab(i1,i2,1,2)=-mab(i1,i2,1,2)*xab2 ! ok
!!$             nab(i1,i2,1,3)=-mab(i1,i2,1,3)*xab4 ! ok
!!$             nab(i1,i2,2,1)=-mab(i1,i2,2,1) ! ok
!!$             nab(i1,i2,2,2)=  27.d0/4.d0 *sqrt(ytemab)*xab2/(1.d0+xab2)**2 &
!!$                  &        *sqdnm ! sqrt
!!$             nab(i1,i2,2,3)= 225.d0/16.d0*ytemab*xab4/(1.d0+xab2)**3*sqdnm ! sign differs
!!$!!!             nab(i1,i2,2,3)= 225.d0/16.d0*sqrt(ytemab)*xab4/(1.d0+xab2)**3 &
!!$!!!                  &           *sqdnm ! sign differs
!!$             nab(i1,i2,3,1)=-mab(i1,i2,3,1) ! ok
!!$!             nab(i1,i2,3,2)= nab(i1,i2,2,3)/ytemab ! wrong (NCLASS)
!!$             nab(i1,i2,3,2)= nab(i1,i2,2,3)/ytemab/xab2
!!$!!!             nab(i1,i2,3,2)= nab(i1,i2,2,3)/xab2
!!$             nab(i1,i2,3,3)=2625.d0/64.d0*sqrt(ytemab)*xab4/(1.d0+xab2)**4 &
!!$                  &        *sqdnm ! sqrt (NCLASS paper ok, source wrong)

             ! N_{ab}
             nab(i1,i2,1,1)=-mab(i1,i2,1,1)
             nab(i1,i2,1,2)=-mab(i1,i2,1,2)*xab2
             nab(i1,i2,1,3)=-mab(i1,i2,1,3)*xab4
             nab(i1,i2,2,1)=-mab(i1,i2,2,1)
             nab(i1,i2,2,2)=  27.d0/4.d0 *ytemab*xab2/(1.d0+xab2)**2*sqdnm
             nab(i1,i2,2,3)= 225.d0/16.d0*ytemab*xab4/(1.d0+xab2)**3*sqdnm
             nab(i1,i2,3,1)=-mab(i1,i2,3,1)
             nab(i1,i2,3,3)=2625.d0/64.d0*ytemab*xab4/(1.d0+xab2)**4*sqdnm
          enddo
       enddo

       ! H&S(4.9) (Ta vTa)^-1N_{ab}^{ij}=(Tb vTb)^-1 N_{ba}^{ji}
       do i1=1,ispc
          do i2=1,ispc
             nab(i1,i2,3,2)=nab(i2,i1,2,3) &
                  &        *tem(i1)/tem(i2)*vtn(i1)/vtn(i2)
             if(i1>i2) then
                nab(i1,i2,2,2)=nab(i2,i1,2,2)&
                     &        *tem(i1)/tem(i2)*vtn(i1)/vtn(i2)
                nab(i1,i2,3,3)=nab(i2,i1,3,3)&
                     &        *tem(i1)/tem(i2)*vtn(i1)/vtn(i2)
             endif
          enddo
       enddo

!!$       ! Check M_{ab}^{0j}=-(Ta vTa)/(Tb vTb)N_{ba}^{0j}=-N_{ab}^{j0}
!!$       do i1=1,ispc
!!$          do i2=1,ispc
!!$             fac=tem(i1)/tem(i2)*vtn(i1)/vtn(i2)
!!$             do j2=1,3
!!$                write(6,'(3I2,1P3E15.7)') i1,i2,j2 &
!!$                     & ,mab(i1,i2,1,j2),-fac*nab(i2,i1,1,j2),-nab(i1,i2,j2,1)
!!$             enddo
!!$          enddo
!!$       enddo

    else

       do i1=1,ispc
          do i2=1,ispc
             xab2=(vtn(i2)/vtn(i1))**2
             xab4=xab2*xab2
             yab=mas(i1)/mas(i2)
             yt =tem(i1)/tem(i2)
             sqdnm=1.d0/sqrt(1.d0+xab2)
             ! M_{ab}
             mab(i1,i2,1,1)=-            (1.d0+yab)/(1.d0+xab2)   *sqdnm ! M^00
             mab(i1,i2,1,2)=- 1.5d0     *(1.d0+yab)/(1.d0+xab2)**2*sqdnm ! M^01
             mab(i1,i2,1,3)=- 15.d0/8.d0*(1.d0+yab)/(1.d0+xab2)**3*sqdnm ! M^02
             mab(i1,i2,2,1)=-(1.5d0     *(1.d0+yab) &
                  &          +2.d0*xab2*(1.d0-yt)*(1.d0+2.5d0*xab2)) &
                  &        /(1.d0+xab2)**2*sqdnm ! M^10
!!             mab(i1,i2,2,1)=-0.5d0*(3.d0+xab2*(4.d0-yt &
!!                  &                +10.d0*(1.d0-yt)*xab2))/(1.d0+xab2)**2*sqdnm ! M^10
             mab(i1,i2,2,2)=-0.25d0*(13.d0+xab2*(20.d0+9.d0*yt &
                  &        +xab2*(52.d0-6.d0*yt+xab2*30.d0*yt))) &
                  &        /(1.d0+xab2)**3*sqdnm ! M^11
             mab(i1,i2,2,3)=-3.d0/16.d0*(23.d0+xab2*(36.d0+19.d0*yt &
                  &                     +xab2*(2.d0*(59.d0-yt)+xab2*84.d0*yt))) &
                  &        /(1.d0+xab2)**4*sqdnm ! M^12
             mab(i1,i2,3,1)=-(15.d0/8.d0*(1.d0+yab) &
                  &          +xab2*3.d0*(1.d0-yt)*(1.d0+3.5d0*xab2)) &
                  &        /(1.d0+xab2)**3*sqdnm ! M^20
!!             mab(i1,i2,3,1)=-3.d0/8.d0*(5.d0+xab2*(8.d0-3.d0*yt &
!!                  &                    +28.d0*(1.d0-yt)*xab2))/(1.d0+xab2)**3*sqdnm ! M^20
             mab(i1,i2,3,2)=-1.d0/16.d0*(69.d0+xab2*(184.d0-19.d0*yt  &
                  &         +xab2*(348.d0*(2.d0-yt)+xab2*(84.d0*(4.d0-yt) &
                  &         +xab2* 280.d0*(1.d0-yt)))))      /(1.d0+xab2)**4*sqdnm ! M^21
             mab(i1,i2,3,3)=-1.d0/64.d0*(433.d0+xab2*(1288.d0+233.d0*yt &
                  &       +xab2*4.d0*(1329.d0-139.d0*yt+xab2* 2.d0*(380.d0+303.d0*yt &
                  &       +xab2*7.d0*(  59.d0-  2.d0*yt+xab2*25.d0*yt))))) &
                  &                                          /(1.d0+xab2)**5*sqdnm ! M^22
             ! N_{ab}
             nab(i1,i2,1,1)=-mab(i1,i2,1,1)      ! N^00, perfect sym.
             nab(i1,i2,1,2)=-mab(i1,i2,1,2)*xab2 ! N^01
             nab(i1,i2,1,3)=-mab(i1,i2,1,3)*xab4 ! N^02
             nab(i1,i2,2,1)= 1.5d0      *(1.d0+(3.d0*yt-2.d0)*xab2) &
                  &                     /(1.d0+xab2)**2*sqdnm ! N^10
             nab(i1,i2,2,2)= 2.25d0     *(3.d0+(5.d0*yt-2.d0)*xab2) &
                  &                      *xab2/(1.d0+xab2)**3*sqdnm ! N^11
             nab(i1,i2,2,3)= 45.d0/16.d0*(5.d0+(7.d0*yt-2.d0)*xab2) &
                  &                      *xab4/(1.d0+xab2)**4*sqdnm ! N^12
             nab(i1,i2,3,1)= 15.d0/8.d0 *(1.d0+(5.d0*yt-4.d0)*xab2) &
                  &                     /(1.d0+xab2)**3*sqdnm ! N^20
             nab(i1,i2,3,2)= 75.d0/16.d0*(3.d0+(7.d0*yt-4.d0)*xab2) &
                  &                      *xab2/(1.d0+xab2)**4*sqdnm ! N^21
             nab(i1,i2,3,3)=525.d0/64.d0*(5.d0+(9.d0*yt-4.d0)*xab2) &
                  &                      *xab4/(1.d0+xab2)**5*sqdnm ! N^22
          enddo
       enddo
!      write(6,*)

!!$       do i1=1,ispc
!!$          do i2=1,ispc
!!$             xab =vtn(i2)/vtn(i1)
!!$             xab2=xab*xab
!!$             yab=mas(i1)/mas(i2)
!!$             yt =tem(i1)/tem(i2)
!!$             if(i1>i2) then
!!$                do j2=1,3
!!$                   mab(i1,i2,2,j2)=mab(i1,i2,2,j2)*tem(i2)/tem(i1)
!!$                   mab(i1,i2,3,j2)=mab(i1,i2,3,j2)*tem(i2)/tem(i1)
!!$                enddo
!!$             endif
!!$             write(6,'(2I2,1P4E15.7)') i1,i2 &
!!$                  & ,nab(i1,i2,1,2),nab(i2,i1,2,1)*yt/xab &
!!$                  & ,nab(i1,i2,1,2)/(nab(i2,i1,2,1)*yt/xab),yt
!!$             write(6,'(2I2,1P4E15.7)') i1,i2 &
!!$                  & ,nab(i1,i2,1,3),nab(i2,i1,3,1)*yt/xab &
!!$                  & ,nab(i1,i2,1,3)/(nab(i2,i1,3,1)*yt/xab),yt
!!$             write(6,'(2I2,1P6E15.7)') i1,i2 &
!!$                  & ,nab(i1,i2,2,2),nab(i2,i1,2,2)*yt/xab &
!!$                  & ,nab(i1,i2,2,2)/(nab(i2,i1,2,2)*yt/xab),yt
!!$!                  & ,(3.d0+(5.d0*yt-2.d0)*xab2)/(3.d0*yt*xab2+(5.d0-2.d0*yt))
!!$             write(6,'(2I2,1P4E15.7)') i1,i2 &
!!$                  & ,nab(i1,i2,2,3),nab(i2,i1,3,2)*yt/xab &
!!$                  & ,nab(i1,i2,2,3)/(nab(i2,i1,3,2)*yt/xab),yt
!!$             write(6,'(2I2,1P4E15.7)') i1,i2 &
!!$                  & ,nab(i1,i2,3,3),nab(i2,i1,3,3)*yt/xab &
!!$                  & ,nab(i1,i2,3,3)/(nab(i2,i1,3,3)*yt/xab),yt
!!$          enddo
!!$       enddo

!!$       do i1=1,ispc
!!$          do i2=1,ispc
!!$             xab=vtn(i2)/vtn(i1) ! chi
!!$             yab=mas(i1)/mas(i2) ! 1/mu
!!$             yt =tem(i1)/tem(i2) ! 1/theta
!!$             fac=yab**2/xab**3 ! F
!!$          ! Check momentum conservation
!!$             write(6,'(2I3,1P7E15.7)') i1,i2 &
!!$                  & ,nab(i2,i1,1,1)*yt/xab+mab(i1,i2,1,1) &
!!$                  & ,nab(i2,i1,1,2)*yt/xab+mab(i1,i2,1,2) &
!!$                  & ,nab(i2,i1,1,3)*yt/xab+mab(i1,i2,1,3) &
!!$                  & ,nab(i1,i2,1,1)       +mab(i1,i2,1,1) & ! approx. sym. used
!!$                  & ,nab(i1,i2,2,1)       +mab(i1,i2,1,2) & ! approx. sym. used
!!$                  & ,nab(i1,i2,3,1)       +mab(i1,i2,1,3) & ! approx. sym. used
!!$                  enddo
!!$          enddo
!!$          stop
    endif
!-----------------------------------------------------------------------
!     friction coef. lab normalized by (na*ma/tauaa)
!-----------------------------------------------------------------------
    lab(1:isp,1:isp,1:2,1:2) = 0.0_8

    if( imodel(6) >= 2 ) then
       ! alpha^1, alpha^2, Delta
       do i1=1,ispc
          denom = nab(i1,i1,3,3)
          do k=1,ispc
!             denom=denom+(zz(k)/zz(i1))**2*den(k)/den(i1)*mab(i1,k,3,3)
             denom=denom+ztau(i1,i1)/ztau(i1,k)*mab(i1,k,3,3)
          enddo
          do i2=1,ispc
!             fac=(zz(i2)/zz(i1))**2*den(i2)/den(i1) ! tau_{i1i1}/tau_{i1i2}
             fac=ztau(i1,i1)/ztau(i1,i2)
             do j2=1,2
                alph(i1,i2,j2)=fac*nab(i1,i2,3,j2)/denom
             enddo
             if(i1 == i2)then
                do j2=1,2
                   alph(i1,i2,j2)=alph(i1,i2,j2) &
!                        &        +sum( (zz(1:ispc)/zz(i1))**2*den(1:ispc)/den(i1) &
                        &        +sum(ztau(i1,i1)/ztau(i1,1:ispc) &
                        &        *mab(i1,1:ispc,3,j2))/denom
                enddo
                emat(i1,i2)=1.d0
                delt(i1,i2)=0.d0 ! Delta_{ab}
             else
                emat(i1,i2)=fac*nab(i1,i2,3,3)/denom ! -Delta_{ab}
                delt(i1,i2)=-emat(i1,i2) ! Delta_{ab}
             endif
          enddo
       enddo
       ! compute inverse matrix for ub2 (higher order) correction
       call matslv(ispc,ispc,emat,fmat,err,ill)
       ! alpha_hat^1, alpha_hat^2
       do i1=1,ispc
          do i2=1,ispc
             do j2=1,2
                alphhat(i1,i2,j2)=sum(fmat(i1,1:ispc)*alph(1:ispc,i2,j2))
             enddo
          enddo
       enddo

!!$       ! check alpha_hat and Delta
!!$       i1=2;i2=3
!!$       write(200,'(1P6E15.7)') mas(i2)/mas(i1),alphhat(i1,i2,1) &
!!$            & ,alphhat(i1,i2,2),delt(i1,i2) &
!!$            & ,alph(i1,i2,1)+sum(delt(i1,1:ispc)*alph(1:ispc,i2,1)) &
!!$            & ,alph(i1,i2,2)+sum(delt(i1,1:ispc)*alph(1:ispc,i2,2))

!!$       allocate(lab3(isp,isp,2))
!!$       ! compute correction friction forces
!!$       do i1=1,ispc
!!$          do i2=1,ispc
!!$!             fac=(zz(i2)/zz(i1))**2*den(i2)/den(i1)
!!$             fac=ztau(i1,i1)/ztau(i1,i2)
!!$             do j1=1,2
!!$                xl=fac*nab(i1,i2,j1,3)
!!$                if(i1 == i2)then
!!$                   do k=1,ispc
!!$!                      xl=xl+(zz(k)/zz(i1))**2*den(k)/den(i1)*mab(i1,k,j1,3)
!!$                      xl=xl+ztau(i1,i1)/ztau(i1,k)*mab(i1,k,j1,3)
!!$                   enddo
!!$                endif
!!$                lab3(i1,i2,j1)=xl
!!$             enddo
!!$          enddo
!!$       enddo
!!$       do i1=1,ispc
!!$          do i2=1,ispc
!!$             do j1=1,2
!!$                do j2=1,2
!!$                   xl=0.d0
!!$                   do k=1,ispc
!!$                      xl=xl+lab3(i1,k,j1)*alphhat(k,i2,j2)
!!$                   enddo
!!$                   lab(i1,i2,j1,j2)=-xl
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$       deallocate(lab3)

       do i1=1,ispc
          do i2=1,ispc
             do j1=1,2
!                xlm=sum( (zz(1:ispc)/zz(i1))**2*den(1:ispc)/den(i1) &
                xlm=sum( ztau(i1,i1)/ztau(i1,1:ispc) &
                     &  *mab(i1,1:ispc,j1,3))
                do j2=1,2
!                   xln=sum( (zz(1:ispc)/zz(i1))**2*den(1:ispc)/den(i1) &
                   xln=sum( ztau(i1,i1)/ztau(i1,1:ispc) &
                        &  *nab(i1,1:ispc,j1,3)*alphhat(1:ispc,i2,j2))
                   lab(i1,i2,j1,j2)=-(xlm*alphhat(i1,i2,j2)+xln)
                enddo
             enddo
          enddo
       enddo
    endif

!!    Old way, when ztau was not used.
!    do i1=1,ispc
!       do i2=1,ispc
!          fac=(zz(i2)/zz(i1))**2*den(i2)/den(i1)
!          do j1=1,2
!             do j2=1,2
!                xl=fac*nab(i1,i2,j1,j2)
!                if(i1.eq.i2)then
!                   do k=1,ispc
!                      xl=xl+(zz(k)/zz(i1))**2*den(k)/den(i1)*mab(i1,k,j1,j2)
!                   enddo
!                endif
!                lab(i1,i2,j1,j2)=lab(i1,i2,j1,j2)+xl
!             enddo
!          enddo
!       enddo
!    enddo
    do i1=1,ispc
       do i2=1,ispc
          fac=ztau(i1,i1)/ztau(i1,i2)
          do j1=1,2
             do j2=1,2
                lab(i1,i2,j1,j2)=lab(i1,i2,j1,j2)+fac*nab(i1,i2,j1,j2)
             enddo
          enddo
       enddo
    enddo
    forall (i1=1:ispc,i2=1:ispc,j1=1:2,j2=1:2,i1==i2) &
         & lab(i1,i2,j1,j2)=lab(i1,i2,j1,j2)+sum(ztau(i1,i1)/ztau(i1,:)*mab(i1,:,j1,j2))

!!$    do i1=1,ispc
!!$       do i2=1,ispc
!!$          do j1=1,2
!!$             do j2=1,2
!!$                write(6,*) i1,i2,j1,j2,lab(i1,i2,j1,j2)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo

    ! Check dependence of alphhat on zast

    do i1=1,ispc
       zast(i1)=0.d0
       do i2=1,ispc
!          fac=(zz(i2)/zz(i1))**2*den(i2)/den(i1)
          fac=ztau(i1,i1)/ztau(i1,i2)
          if( mas(i2) > mas(i1) ) zast(i1)=zast(i1)+fac
       enddo
    enddo
!!$    i1=2
!!$    zeff=0.d0
!!$    do i=2,ispc
!!$       zeff=zeff+zz(i)**2*den(i)
!!$    enddo
!!$    zeff=zeff/den(1)
!!$!    write(6,'(A,1P4E15.7)') "== ",mas(2),mas(3),den(2),den(3)
!!$    write(6,'(4F9.5,1P8E13.5)') zast(i1),zeff &
!!$         & ,den(2)/den(1),den(3)/den(1) &!,den(4)/den(1)
!!$         & ,alphhat(i1,i1,1),0.28d0*zast(i1)/(0.59d0+zast(i1)) &! (32c)
!!$!         & ,alph(i1,i1,1)+sum(delt(i1,1:ispc)*alph(1:ispc,i1,1))
!!$         & ,alphhat(i1,i1,2),(0.16d0+0.64d0*zast(i1))/(0.59d0+zast(i1)) &! (32d)
!!$         & ,alph(i1,i1,2)+sum(delt(i1,1:ispc)*alph(1:ispc,i1,2))

!!$    ! Check Z^* dependence of chi, compared to 25/(4C_3)
!!$    alp =zz(3)**2*den(3)/den(2)
!!$!    write(996,'(F10.5,1P3E15.7)')alp,den(2),lab(2,2,2,1),-lab(2,2,2,2)
!!$    write(996,'(F10.5,1P4E15.7)')alp,lab(2,2,2,1),-lab(2,2,2,2) &
!!$         & ,sum(lab(2,1:3,2,1)*zz(2)/zz(1:3)) &
!!$         & ,-sum(lab(2,1:3,2,2)*zz(2)/zz(1:3))

!!$    ! Check Z^* dependence of D, compared to C_1+C_2^2/C_3
!!$    alp =zz(3)**2*den(3)/den(2)
!!$    write(996,'(F10.5,1P4E15.7)')alp,-lab(2,2,1,1) &
!!$         & ,-sum(lab(2,1:3,1,1)*zz(2)/zz(1:3))

!!$    ! Hirshman 1977, Eqs.(49)(50) for D-D friction coef.
!!$    xl=0.d0
!!$    do i1=1,ispc
!!$       if(mas(i1)<mas(2)) then
!!$          c1=1.d0-(15.d0/8.d0)**2*zast(i1)/(45.d0/16.d0*sqrt(2.d0) &
!!$               & +(433.d0/64.d0)*zast(i1))
!!$          xl=xl+sqrt(mas(i1)/mas(2)*tem(2)/tem(i1))*tem(2)/tem(i1) &
!!$               & *(c1+(1.d0-den(2)*zz(2)**2/(den(i1)*zz(i1)**2)/zast(i1)) &
!!$               & *(1.d0-c1))
!!$       endif
!!$    enddo
!!$    c1a=1.d0-(15.d0/8.d0)**2*zast(2)/(45.d0/16.d0*sqrt(2.d0) &
!!$         & +(433.d0/64.d0)*zast(2))
!!$    c2=1.5d0-(15.d0/8.d0)*(0.75d0*sqrt(2.d0)+(69.d0/16.d0)*zast(2)) &
!!$         & /(45.d0/16.d0*sqrt(2.d0)+433.d0/64.d0*zast(2))
!!$    c3=sqrt(2.d0)+13.d0/4.d0*zast(2)-(0.75d0*sqrt(2.d0) &
!!$         & +69.d0/16.d0*zast(2))**2 &
!!$         & /(45.d0/16.d0*sqrt(2.d0)+433.d0/64.d0*zast(2))
!!$    write(201,'(1P6E15.7)') tem(2)/tem(1),zast(2) &
!!$         & ,-zast(2)*c1a-xl,-zast(2)*c2 &
!!$         & ,-zast(2)*c2,-c3

!-----------------------------------------------------------------------      
    if( idebug == 1 ) then
       allocate(labe(isp,isp,3,3))
       ! check lab
       ! Normalized f_{ij}^{ab} by deuterium parameters
       do i1=1,ispc
          do i2=1,ispc
!             xl=den(i1)/den(i2)*(zz(i1)/zz(i2))**2
             xl=ztau(i1,i2)/ztau(i1,i1) 
             do j1=1,2
                do j2=1,2
                   ! renormalized by electron parameters
                   labe(i1,i2,j1,j2)=lab(i1,i2,j1,j2)*xl
!!!ele.renorm     &             *(den(1)*mas(1)*cnmp/ztau(1,1))
                enddo
             enddo
          enddo
       enddo

!!$       i1=2;i2=3
!!$       write(200,'(1P6E15.7)') mas(i2)/mas(i1) &
!!$            & ,labe(i1,i2,1,1)/zast(i1),labe(i1,i2,1,2)/zast(i1) &
!!$            & ,labe(i1,i2,2,1)/zast(i1),labe(i1,i2,2,2)/zast(i1)
!!$       write(200,'(1P6E15.7)') mas(i2)/mas(i1) &
!!$            & ,labe(i1,i2,1,1),labe(i1,i2,1,2) &
!!$            & ,labe(i1,i2,2,1),labe(i1,i2,2,2)

!!$       i1=2;i2=2
!!$       write(200,'(1P6E15.7)') tem(2)/tem(1),zast(i1) &
!!$            & ,labe(i1,i2,1,1),labe(i1,i2,1,2) &
!!$            & ,labe(i1,i2,2,1),labe(i1,i2,2,2)

       allocate(labtmp(isp,isp,3,3),sumbc2(isp))
       ! Recompute original l_{ij}^{ab} as labtmp(i1,i2,j1,j2)
       do i1=1,ispc
          zast(i1)=0.d0
          do i2=1,ispc
!             fac=(zz(i2)/zz(i1))**2*den(i2)/den(i1)
             fac=ztau(i1,i1)/ztau(i1,i2)
             if( mas(i2) > mas(i1) ) zast(i1)=zast(i1)+fac
             do j1=1,3
                do j2=1,3
                   xl=fac*nab(i1,i2,j1,j2)
                   if(i1.eq.i2)then
                      do k=1,ispc
!                         xl=xl+(zz(k)/zz(i1))**2*den(k)/den(i1)*mab(i1,k,j1,j2)
                         xl=xl+ztau(i1,i1)/ztau(i1,k)*mab(i1,k,j1,j2)
                      enddo
                   endif
                   labtmp(i1,i2,j1,j2)=xl
                enddo
             enddo
          enddo
       enddo

! To confirm l_{ij}^{ab} = l_{ji}^{ba}, we introduce the factor
! "xl" that normalizes l_{ij}^{ab} by electron parameters.
! Note that "xl" is originally normalized by parameters regarding
! species "a".
       do i1=1,ispc
!!$          xl=(den(i1)/den(1)*(zz(i1)/zz(1))**2)**2 &
!!$               & *sqrt(mas(i1)/mas(1))*(tem(1)/tem(i1))**1.5d0 ! tau_{11}/tau_{i1i1}
          do i2=1,ispc
!             xl=den(i1)/den(i2)*(zz(i1)/zz(i2))**2 &
!                  & *(tem(1)/tem(i1))*vtn(1)/vtn(i1) ! H&S
             xl=ztau(i1,i2)/ztau(i1,i1) &
                  & *(tem(1)/tem(i1))*vtn(1)/vtn(i1) ! H&S
             do j1=1,3
                do j2=1,3
                   labe(i1,i2,j1,j2)=labtmp(i1,i2,j1,j2)*xl
!!!ele.renorm     &             *(den(1)*mas(1)*cnmp/ztau(1,1))
                enddo
             enddo
          enddo
       enddo
       write(6,*) 

       ! Check l_{ij}^{ab} = l_{ji}^{ba}
       write(6,'(1X,4(A,4X),1X,A,2(9X,A))') "a b i j" &
            & ,"lorg_{ij}^{ab}","l_{ij}^{ab}","b a j i","l_{ji}^{ba}","ratio","difference"
       do i1=1,ispc
          do i2=1,ispc
             xab=vtn(i2)/vtn(i1)
             yab=mas(i1)/mas(i2)
             yt =tem(i1)/tem(i2)
             fac=xab/yt
             do j1=1,3
                do j2=1,3
                   write(6,'(4I2,2(2X,1PE15.7),2X,4I2,3(2X,E15.7))') &
                        &  i1,i2,j1,j2,labtmp(i1,i2,j1,j2),labe(i1,i2,j1,j2) &
                        & ,i2,i1,j2,j1,labe(i2,i1,j2,j1) &
                        & ,labe(i1,i2,j1,j2)/labe(i2,i1,j2,j1) &
                        & ,(labe(i1,i2,j1,j2)-labe(i2,i1,j2,j1))! &
!                        & /labe(i1,i2,j1,j2)
                enddo
             enddo
          enddo
       enddo
       write(6,*) 

       ! Normalized f_{ij}^{ab} by electron parameters
       do i1=1,ispc
!          xl=(den(i1)/den(1)*(zz(i1)/zz(1))**2)**2 &
!               & *sqrt(mas(i1)/mas(1))*(tem(1)/tem(i1))**1.5d0 ! tau_{11}/tau_{i1i1}
          xl=ztau(1,1)/ztau(i1,i1)
          do i2=1,ispc
             do j1=1,2
                do j2=1,2
                   ! renormalized by electron parameters
                   labe(i1,i2,j1,j2)=lab(i1,i2,j1,j2)*xl
!!!ele.renorm     &             *(den(1)*mas(1)*cnmp/ztau(1,1))
                enddo
             enddo
          enddo
       enddo

       ! Check f_{ij}^{ab} = f_{ji}^{ba}
       write(6,'(1X,4(A,4X),1X,A,9X,A)') "a b i j" &
            & ,"forg_{ij}^{ab}","f_{ij}^{ab}","b a j i","f_{ji}^{ba}","ratio"
       do i1=1,ispc
          do i2=1,ispc
             do j1=1,2
                do j2=1,2
                   write(6,'(4I2,2(2X,1PE15.7),2X,4I2,2(2X,E15.7))') &
                        &  i1,i2,j1,j2,lab(i1,i2,j1,j2),labe(i1,i2,j1,j2) &
                        & ,i2,i1,j2,j1,labe(i2,i1,j2,j1) &
                        & ,labe(i1,i2,j1,j2)/labe(i2,i1,j2,j1)
                enddo
             enddo
          enddo
       enddo
       write(6,*) 

       ! Compare f_{ij}^{ab} and l_{ij}^{ab}
       write(6,'(1X,A,3X,A,4X,A,7X,A)') "a b i j","f_{ij}^{ab}","l_{ij}^{ab}","ratio"
       do i1=1,ispc
          do i2=1,ispc
             do j1=1,2
                do j2=1,2
                   write(6,'(4I2,1P3E15.7)') i1,i2,j1,j2,labe(i1,i2,j1,j2) &
                        & ,labtmp(i1,i2,j1,j2),labe(i1,i2,j1,j2)/labtmp(i1,i2,j1,j2)
                enddo
             enddo
          enddo
       enddo
       write(6,*) 

       write(6,'(1X,A,6X,A)') "a","Z_a^*"
       do i1=1,ispc
          write(6,'(I2,1PE15.7)') i1,zast(i1)
       enddo
       write(6,*) 

       if( imodel(6) /= 0 ) then
          ! Check alpha_hat^1, alpha_hat^2 against Hirshman PoF (1978) 1295
          write(6,*) "Check alpha_hat and beta_hat"
          write(6,'(1X,A,4X,A,6X,A,5X,A)') "j a b","numerical","Hirshman","error%"
          do i1=1,ispc
             do i2=1,ispc
                epsab =(mas(i1)*zz(i2))/(mas(i2)*zz(i1)) &
                     & /( (tem(i1)*vt(i1))/(tem(i2)*vt(i2)) &
                     & *(den(i2)*zz(i2)**2/(den(i1)*zz(i1)**2))**2 )
                j2=1
                if(i1==i2) then
                   f=0.28d0*zast(i1)/(0.59d0+zast(i1)) ! (32c)
                   write(6,'(3I2,1P2E15.7,-1PF7.2)') j2,i1,i2 &
                        & ,alphhat(i1,i2,j2)*epsab &
                        & ,f,abs(f/(alphhat(i1,i2,j2)*epsab)-1.d0)*1.d3
                else
                   write(6,'(3I2,1P2E15.7,6X,A)') j2,i1,i2 &
                        & ,alphhat(i1,i2,j2)*epsab,0.d0,"-"
                endif
                j2=2
                if(i1==i2) then
                   f=(0.16d0+0.64d0*zast(i1))/(0.59d0+zast(i1)) ! (32d)
                   write(6,'(3I2,1P2E15.7,-1PF7.2)') j2,i1,i2 &
                        & ,alphhat(i1,i2,j2)*epsab,f,abs(f/(alphhat(i1,i2,j2)*epsab)-1.d0)*1.d3
                else
                   write(6,'(3I2,1P2E15.7,6X,A)') j2,i1,i2 &
                        & ,alphhat(i1,i2,j2)*epsab,0.d0,"-"
                endif
             enddo
          enddo
          write(6,*) 
       endif

       ! Check c_2^{ab} against Hirshman PoF (1978) 589, 1295
       do i1=1,ispc
          do i2=1,ispc
             emat(i1,i2)=labe(i1,i2,2,2)
          enddo
       enddo
       ! compute inverse matrix of l_{22}^{ab}, i.e., (l_{22}^{-1})^{ab}
       call matslv(ispc,ispc,emat,fmat,err,ill)

       write(6,*) "Check eps_{ab}*c_2^{ab}"
       write(6,'(1X,A,4X,A,6X,A,5X,A)') "a b","numerical","Hirshman","error%"
       do i1=1,ispc
          sumbc2(i1)=0.d0
          do i2=1,ispc
             epsab =(mas(i1)*zz(i2))/(mas(i2)*zz(i1)) &
                  & /( (tem(i1)*vt(i1))/(tem(i2)*vt(i2)) &
                  & *(den(i2)*zz(i2)**2/(den(i1)*zz(i1)**2))**2 )
             c2 = 0.d0
             do k=1,ispc
!                xl=(den(i1)/den(k)*(zz(i1)/zz(k))**2)**2 &
!                     & *sqrt(mas(i1)/mas(k))*(tem(k)/tem(i1))**1.5d0 ! tau_{kk}/tau_{i1i1}
!                c2 = c2 + fmat(i1,k)*(-lab(k,i2,2,1)*xl)
                c2 = c2 + fmat(i1,k)*(-labe(k,i2,2,1))
             enddo
             sumbc2(i1) = sumbc2(i1) + c2
             if(i1==i2) then
                f=-(0.884d0+0.458d0*zast(i1)) &
                     & /(1.d0+2.966d0*zast(i1)+0.753d0*zast(i1)**2)*zast(i1) ! (32b)
                write(6,'(2I2,1P2E15.7,-1PF7.2)') i1,i2,c2*epsab &
                     & ,f,abs(f/(c2*epsab)-1.d0)*1.d3
             else
                write(6,'(2I2,1P2E15.7,6X,A)') i1,i2,c2*epsab,0.d0,"-"
             endif
          enddo
       enddo

       write(6,'(1X,A,2X,A)') "a","sum_b c2^{ab}"
       do i1=1,ispc
          write(6,'(I2,1PE15.7)') i1,sumbc2(i1)
       enddo
!!$       write(6,*)
!!$
!!$       alp =zz(3)**2*den(3)/den(2)
!!$       barK1=-alp*(0.83d0+0.42d0*alp)/(0.58d0+alp)
!!$       barK2=1.58d0*(1.d0+1.33d0*alp*(1.d0+0.6d0*alp)/(1.d0+1.79d0*alp))
!!$       barK2mod=1.58d0*(1.13d0/alp+0.5d0+0.56d0/(0.56d0+alp))*alp
!!$       write(6,'(F10.5,1P12E15.7)') alp,lab(2,2,2,1),barK1,-lab(2,2,2,2),barK2,barK2mod

       deallocate(labtmp,sumbc2)
       deallocate(labe)
       stop
    endif
!-----------------------------------------------------------------------      

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
             xb2=xb**2
             if(xb < 5.922d0) then ! above this value, fp is always 1 in double precision.
                fp=erf(xb)
                fdp=2.d0/sqrt(cnpi)*exp(-xb**2) ! derivative of erf
             else
                fp=1.d0
                fdp=0.d0
             endif
             gfun=(fp-xb*fdp)/(2.d0*xb2) ! Chandrasekhar function
             xtab=1.d0/ztau(i1,i2)
             zmud=zmud+xtab*(fp-gfun)/xa3
             zmut=zmut+xtab*((fp-3.d0*gfun)/xa3 &
                  & +4.d0*(tem(i1)/tem(i2)+(vt(i1)/vt(i2))**2)*gfun/xa)
          enddo
          zmud=zmud*zp43
          zkb=xgt*zmud/sqzfac(i1)**1.5d0
          zmut=zmut*zp43
          zkps=0.d0
          do m=1,mxneo
             fm=dfloat(m)
             zomg=xa*vt(i1)*fm*gmneo
             zu=zmut/zomg
             zu2=zu*zu
             zui=1.d0/zu
             zat=datan(zui)
             zmui=-1.5d0*zu2-4.5d0*zu2*zu2 &
                  & +(0.25d0+(1.5d0+2.25d0*zu2)*zu2)*2.d0*zu*zat
             if(zu.gt.100.d0)zmui=0.4d0
             zkps=zkps+fmneo(m)*zmui
          enddo
          zkps=1.5d0*(vt(i1)*xa)**2*zkps/zmut
!-----
!         zkps becomes zero at axis because fmneo becomes zero despite finite zmui.
!         On the other hand, zkb is still finite at axis only if the potate orbit effect
!         is taken into account. In short, the trapped particle fraction f_t becomes zero 
!         unless the potate orbit effect is taken into account, leading to zkb being zero, 
!         but in fact f_t remains finite even at axis.
!
!         For such a case, the following form of zk unfavorably goes to zero.
!         To avoid this, finite fmneo is given even at axis.
!
          if( zkb+zkps /= 0.d0 ) then
             zk=zkb*zkps/(zkb+zkps)
          else
             zk=0.d0
          end if
!          if(i1==1) write(6,'(1P4E15.7)') xa,tpfneo,zkb,zkps
          if(imodel(5).eq. 1)zk=zkb
          if(imodel(5).eq.-1)zk=zkps
!-----
          yk0=zk*dexp(-xa2)*xa2*xa2
          yk1=yk0*xa2
          yk2=yk1*xa2
          xke11=xke11+yk0
          xke12=xke12+yk1
          xke22=xke22+yk2
       enddo
!-----------------------------------------------------------------------
!         renormalized by tau_aa / (ma * na) for consistency against booth9
       renorm = ztau(i1,i1) / (zmas(i1) * den(i1))
       xfac = dx * 2.0_8 / zp43 * zmas(i1) * den(i1) * renorm
       xke11 = xke11 * xfac
       xke12 = xke12 * xfac
       xke22 = xke22 * xfac
       xmu(i1,1,1)=xke11
       xmu(i1,1,2)=2.5d0*xke11-xke12
       xmu(i1,2,1)=xmu(i1,1,2)
       xmu(i1,2,2)=25.d0/4.d0*xke11-5.d0*xke12+xke22
    enddo
    do i1 = ispc+1, isp
       xmu(i1,1,1) = 0.0d0
       xmu(i1,1,2) = 0.0d0
       xmu(i1,2,1) = 0.0d0
       xmu(i1,2,2) = 0.0d0
    enddo
!=======================================================================
!  fast ion contribution
!
!     The model of a set for beam friction forces is partly taken from
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
          wc3=zp43*vtn(1)**3*mas(1)/mas(ib)*z1
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
             teff=ec*sint/dlog((x32+1.d0)/xp321)
          enddo
          eb=ec*(x0+dx*(tem(ib)-teff0)/(teff-teff0))
          vb=sqrt(2.d0*eb/mas(ib))
          vb2=vb*vb
          vb3=vb2*vb
!-----
          vbc3=vb3+vc3
!<e-b,1,1> fast ion-electron particle friction
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
!_org              lab( 1,ib,1,2)=1.5d0*lab(1,ib,1,1)
!_org              lab( 1, 1,1,2)=lab( 1, 1,1,2)-lab(1,ib,1,2)
!_org              lab(ib, 1,1,2)=lab(ib, 1,1,2)-lab(1,ib,1,2)
!_org              lab(ib,ib,1,2)=lab(ib,ib,1,2)+lab(1,ib,1,2)
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
!_org                lab(ib, 1,2,1)=lab(ib, 1,2,1)-lab(1,ib,2,1)
!_org                lab(ib,ib,2,1)=lab(ib,ib,2,1)+lab(1,ib,2,1)
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
             x=(dlog(vbc3)-3.d0*dlog(vb+vc) &
                  & +2.d0*sq3*datan((2.d0*vb-vc)/sq3/vc)+cnpi/sq3)/(6.d0*vc)
             do i=2,ispc
                y=(zz(ib)/zz(i))**2*den(ib)/den(i)*zp43/dlog(vbc3/vc3)
                lab( i,ib,1,1)=y*(1.d0+mas(i)/mas(ib))*vtn(i)**3/vc3
                lab( i, i,1,1)=lab(i, i,1,1)-lab(i,ib,1,1)
                lab(ib, i,1,1)=lab(i,ib,1,1) &
                     & *(den(i)*zz(i)**2/den(1))**2 &
                     & *sqrt(mas(i)/mas(1))*(tem(1)/tem(i))**1.5d0
                lab(ib,ib,1,1)=lab(ib,ib,1,1)-lab(ib,i,1,1)
!<i-b,1,2>
                if(imodel(3).eq.0) then
                   lab( i,ib,1,2)=void
                   lab( i, i,1,2)=lab( i, i,1,2)
                   lab(ib, i,1,2)=0.d0
                   lab(ib,ib,1,2)=void
!<i-b,2,1>
                   if(imodel(2).eq.1)then
                      lab( i,ib,2,1)=y*5.d0*vtn(i)*x
!                      lab( i,ib,2,1)=y*vtn(i)*x
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
                  & *(1.d0+(1.d0+z1/z2)*vc3/(vb2/x-2.d0*vc3))
          else
             xmu(ib,1,1)=xgt*(zz(ib)/zz(1))**2*den(ib)/den(1) &
                  & *z1/z2*vc3/(vb2/x-2.d0*vc3)
          endif
!>>
          xmu(ib,1,1)=(1.d0-fbeam)*xmu(ib,1,1)*sqzfac(ib)**(-1.5d0)
!>>
!<<dummy set for parallel thermal friction>>
          lab(ib,ib,2,2)=1.d0
       enddo
    endif
!=======================================================================
!<amat>
    do i1=1,isp
       do i2=1,isp
          amat(i1    ,i2    )=lab(i1,i2,1,1)
          amat(i1    ,i2+isp)=lab(i1,i2,1,2)
          amat(i1+isp,i2    )=lab(i1,i2,2,1)
          amat(i1+isp,i2+isp)=lab(i1,i2,2,2)
       enddo
    enddo
!-----
!!$    bmat=amat
!!$    do i1=1,isp
!!$       bmat(i1    ,i1    )=bmat(i1    ,i1    )-xmu(i1,1,1)
!!$       bmat(i1    ,i1+isp)=bmat(i1    ,i1+isp)-xmu(i1,1,2)
!!$       bmat(i1+isp,i1    )=bmat(i1+isp,i1    )-xmu(i1,2,1)
!!$       bmat(i1+isp,i1+isp)=bmat(i1+isp,i1+isp)-xmu(i1,2,2)
!!$    enddo
!!$!=======================================================================
!!$!<cmat>
!!$    call matslv(isp2,isp2,bmat,cmat,err,ill)
!-----
!!$    ! xmu borrowed from booth9
!!$    xmu=0.d0
!!$    do i1=1,isp
!!$       xmu(i1,1,1)=bmat(i1    ,i1    )
!!$       xmu(i1,1,2)=bmat(i1    ,i1+isp)
!!$       xmu(i1,2,1)=bmat(i1+isp,i1    )
!!$       xmu(i1,2,2)=bmat(i1+isp,i1+isp)
!!$    enddo

    cmat = amat
    do i1=1,isp
       cmat(i1    ,i1    )=amat(i1    ,i1    )-xmu(i1,1,1)
       cmat(i1    ,i1+isp)=amat(i1    ,i1+isp)-xmu(i1,1,2)
       cmat(i1+isp,i1    )=amat(i1+isp,i1    )-xmu(i1,2,1)
       cmat(i1+isp,i1+isp)=amat(i1+isp,i1+isp)-xmu(i1,2,2)
    enddo
!=======================================================================
!<cmat> : (L - M)^{-1}
    call invmrd(cmat,isp2,isp2,ill)
!    --- Replace invmrd by the following lines when using LAPACK ---
!    allocate(ipiv(isp2))
!    call getrf( cmat, ipiv, ill )
!    call getri( cmat, ipiv, ill )
!    deallocate(ipiv)
!=======================================================================
!<bmat>
    bmat=0.d0
    do i1=1,isp
       bmat(i1    ,i1    )=xmu(i1,1,1)
       bmat(i1    ,i1+isp)=xmu(i1,1,2)
       bmat(i1+isp,i1    )=xmu(i1,2,1)
       bmat(i1+isp,i1+isp)=xmu(i1,2,2)
    enddo
!=======================================================================
!<alf> : (L - M)^{-1}M
    forall (i1=1:isp2,i2=1:isp2) alf(i1,i2)=sum(cmat(i1,:)*bmat(:,i2))
!    alf=matmul(cmat,bmat) ! This is equivalent to the former line, 
!                          ! but a little bit slow.
!=======================================================================
! Neoclassical bootstrap coefficients : coebsc, i.e., <L_31,L_32>
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
!     neoclassical conductivity : coencc
!=======================================================================
!--matrix method
    zeff=sum(zz(2:isp)**2*den(2:isp))/den(1)

    xc0=0.0_8
    do i1=1,isp
       j1=i1
       ! beam ion is normalized by m_e n_e/tau_ee
       if(ibeam.eq.1.and.i1.eq.isp)j1=1
       ! cmat(a,b) : (L - M)^{-1}_ab
       ! renormalized by multiplying (m_e n_e/tau_ee) / (m_i n_i/tau_ii)
       xc0=xc0+sum(zz(1:isp)*den(1:isp)/den(1)*(-cmat(1:isp,i1))) &
            &    *(mas(1)*den(1)*ztau(j1,j1)/(mas(j1)*den(j1)*ztau(1,1))) &
            &    *zz(i1)*den(i1)/den(1)
    enddo
    coencc=xc0
!=======================================================================
!   coefficients for nbi current drive : coebdc
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
    xmu1  =         xmu(1,1,1)
    xmu2  = xsign * xmu(1,1,2)
    xmu3  =         xmu(1,2,2)
    ell11 = - zeff
    ell12 = 1.5d0 * ell11
    ell22 = ell11 / zeff * (sqrt(2.d0) + 3.25d0 * zeff)
    coebdc(isp+2)=((xmu3 + sqrt(2.d0) + 3.25d0 * zeff) &
         &     *  xmu1 + ( -xmu2 + 1.5d0 * zeff) * xmu2) &
         &     / ((xmu3 + sqrt(2.d0) + 3.25d0 * zeff) &
         &     * (xmu1 + zeff) - (-xmu2 + 1.5d0 * zeff)**2)
!!$    write(101,'(1P5E12.5)') coebdc(1)/coebdc(isp), &
!!$         &     coebdc(2)/coebdc(isp), &
!!$         &     coebdc(3)/coebdc(isp), &
!!$         &     coebdc(4)/coebdc(isp), &
!!$         &     coebdc(5)/coebdc(isp)
!!$    write(101,'(1P5E12.5)') cmat(1,isp),cmat(2,isp),cmat(3,isp),cmat(4,isp)
!!$    write(101,'(1P5E12.5)') xmu1/xgt,xmu2/xgt, &
!!$         &     xmu3/xgt,coebdc(isp+2),coebdc(isp+1)/coebdc(isp)

!=======================================================================
!     diffusion coefficient matrix : dmat
!=======================================================================
    do i1=1,isp2
       do i2=1,isp2
          x=0.d0
          do k=1,isp2
             x=x+amat(i1,k)*alf(k,i2)
          enddo
          dmat(i1,i2)=x
       enddo
    enddo
!-----
    do i1=1,ispc
       do i2=1,ispc
          i3=isp+i2
          dmat(i1,i3)=dmat(i1,i3)-dmat(i1,i2)
       enddo
    enddo
!-----
    do i1=isp+1,isp2
       do i2=1,ispc
          dmat(i1,i2)=-dmat(i1,i2)
       enddo
    enddo
!-----
    do i1=isp+1,isp2
       do i2=1,ispc
          i3=isp+i2
          dmat(i1,i3)=dmat(i1,i3)+dmat(i1,i2)
       enddo
    enddo
!=======================================================================

    deallocate(zmas,zast)
    deallocate(alph,alphhat,emat,fmat,delt)
    deallocate(lab,mab,nab,vt,vtn,xmu)

  end subroutine nccoe

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
!LAPACK  use lapack95, only : getrf, getri
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
!LAPACK  integer(4), allocatable :: ipiv(:)
  real(8) :: dx, eb, ec, err, fac, fdp, fp, gfun &
       &   , sint, sq3, tbdc, teff, teff0 &
       &   , vb, vb2, vb3, vbc3, vc, vc3, vxmax &
       &   , wc3, x, x0, x32, xa, xa2, xa3, xab2, xab4, xb, xc0 &
       &   , xfac, xtab, xgt, xke11, xke12, xke22, xl &
       &   , xp321, xx, yk0, yk1, yk2, z1, z2, z3 &
       &   , z30, zeff, zk, zmud, zmut, zp43, y, zkdenom, dsecnd, s_init
  real(8) :: xsign = -1.d0, void = 0.d0
  integer :: i, i1, i2, ib, ill, isp, isp2, ivmax &
       &   , j, j1, j2, k
!=======================================================================
!!   Standard
!    vxmax = 4.d0 ! org. 3.d0
!    ivmax = 15   ! org. 100
!   High definition
    vxmax = 4.5d0 ! org. 3.d0
    ivmax = 20    ! org. 100
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
          fac=ztau(i1,i1)/ztau(i1,i2)
          do j1=1,2
             do j2=1,2
                lab(i1,i2,j1,j2)=fac*nab(i1,i2,j1,j2)
             enddo
          enddo
       enddo
    enddo
    forall (i1=1:ispc,i2=1:ispc,j1=1:2,j2=1:2,i1==i2) &
         & lab(i1,i2,j1,j2)=lab(i1,i2,j1,j2)+sum(ztau(i1,i1)/ztau(i1,:)*mab(i1,:,j1,j2))

!!$    do i1=1,ispc
!!$       do i2=1,ispc
!!$          do j1=1,2
!!$             do j2=1,2
!!$                write(6,*) i1,i2,j1,j2,lab(i1,i2,j1,j2)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$    write(6,*) 

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
!!$    bmat=amat
!!$    do i1=1,isp
!!$       bmat(i1    ,i1    )=bmat(i1    ,i1    )-xmu(i1,1,1)
!!$       bmat(i1    ,i1+isp)=bmat(i1    ,i1+isp)-xmu(i1,1,2)
!!$       bmat(i1+isp,i1    )=bmat(i1+isp,i1    )-xmu(i1,2,1)
!!$       bmat(i1+isp,i1+isp)=bmat(i1+isp,i1+isp)-xmu(i1,2,2)
!!$    enddo
!!$!=======================================================================
!!$!<cmat> : (L - M)^{-1}
!!$    call matslv(isp2,isp2,bmat,cmat,err,ill)
!-----
    cmat = amat
    do i1=1,isp
       cmat(i1    ,i1    )=amat(i1    ,i1    )-xmu(i1,1,1)
       cmat(i1    ,i1+isp)=amat(i1    ,i1+isp)-xmu(i1,1,2)
       cmat(i1+isp,i1    )=amat(i1+isp,i1    )-xmu(i1,2,1)
       cmat(i1+isp,i1+isp)=amat(i1+isp,i1+isp)-xmu(i1,2,2)
    enddo
!=======================================================================
!<cmat> : (L - M)^{-1}
    call invmrd(cmat,isp2,isp2,ill)
!    --- Replace invmrd by the following lines when using LAPACK ---
!    allocate(ipiv(isp2))
!    call getrf( cmat, ipiv, ill )
!    call getri( cmat, ipiv, ill )
!    deallocate(ipiv)
!=======================================================================
!<bmat> : viscosity matrix, M
    bmat=0.d0
    do i1=1,isp
       bmat(i1    ,i1    )=xmu(i1,1,1)
       bmat(i1    ,i1+isp)=xmu(i1,1,2)
       bmat(i1+isp,i1    )=xmu(i1,2,1)
       bmat(i1+isp,i1+isp)=xmu(i1,2,2)
    enddo
!=======================================================================
!<alf> : (L - M)^{-1}M
    forall (i1=1:isp2,i2=1:isp2) alf(i1,i2)=sum(cmat(i1,:)*bmat(:,i2))
!    alf=matmul(cmat,bmat) ! This is equivalent to the former line, 
!                          ! but a little bit slow.
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
    zeff=sum(zz(2:isp)**2*den(2:isp))/den(1)

    xc0=0.d0
    do i1=1,isp
       j1=i1
       ! beam ion is normalized by m_e n_e/tau_ee
       if(ibeam.eq.1.and.i1.eq.isp)j1=1
       ! cmat(a,b) : (L - M)^{-1}_ab
       ! renormalized by multiplying (m_e n_e/tau_ee) / (m_i n_i/tau_ii)
       xc0=xc0+sum(zz(1:isp)*den(1:isp)/den(1)*(-cmat(1:isp,i1))) &
            &    *(mas(1)*den(1)*ztau(j1,j1)/(mas(j1)*den(j1)*ztau(1,1))) &
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
