!
!	t3d_nccoef_lhd.f90
!	task.20110214
!
!	Created by WAKASA Arimitsu on 11/03/02.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!
! based on ER module version beta-1.2(2009/07/24)
  module t3d_interface_tr2er
    implicit none
    public
    integer, save :: nrmax_tr
    real(8), save :: bb_tr, rr_tr, ra_tr
    real(8), dimension(:), allocatable, save:: rm_tr, epsrho_tr, eh_tr
    real(8), dimension(:), allocatable, save :: mm_tr, pz_tr
    real(8), dimension(:,:), allocatable, save :: rn_tr, rt_tr

  contains

    subroutine t3d_tr2er
      use TRCOMM, only : nrmax, rr, ra, bb, rm, rn, rt, pz, epsrho, pa
      implicit none
      integer :: nr, ns
      if (allocated(rm_tr)) then
         continue
      else
         allocate(rm_tr(nrmax), epsrho_tr(nrmax), eh_tr(nrmax))
         allocate(rn_tr(nrmax,2), rt_tr(nrmax,2), mm_tr(2), pz_tr(2))
      end if

      nrmax_tr = nrmax
      bb_tr = bb
      rr_tr = rr
      ra_tr = ra
      rm_tr(1:nrmax) = rm(1:nrmax)
      epsrho_tr(1:nrmax) = epsrho(1:nrmax)

      do ns = 1, 2
         rn_tr(1:nrmax,ns) = rn(1:nrmax,ns)
         rt_tr(1:nrmax,ns) = rt(1:nrmax,ns)
         mm_tr(ns) = PA(ns)
         pz_tr(ns) = pz(ns)
      end do

      do nr = 1, nrmax
         eh_tr(nr) =  0.2d0*rm(nr)**2
      end do
    end subroutine t3d_tr2er
  end module t3d_interface_tr2er
!---------------------------------------------------------------------------------------

 module erlib_mod
  contains
    subroutine gradient(df, f, x, imax)
      implicit none
      integer, intent(in) :: imax
      real(8), dimension(imax), intent(in) :: f, x
      real(8), dimension(imax), intent(out) :: df
      integer :: i
      real(8) :: dxp, dxm
      
      do i = 2, imax - 1
         dxp = x(i+1) - x(i)
         dxm = x(i) - x(i-1)
         df(i) = (f(i+1)*dxm*dxm - f(i-1)*dxp*dxp + f(i)*(dxp*dxp - dxm*dxm)) &
              / (dxp*dxm*(dxm+dxp))
      end do

      call t3d_interpol_p(df(1), x(1), x(2), x(3), x(4), df(2), df(3), df(4))
      call t3d_interpol_p(df(imax), x(imax), x(imax-3), x(imax-2), x(imax-1), &
           df(imax-3), df(imax-2), df(imax-1))

    end subroutine gradient
    
    subroutine collision_frequency(tau, z, an, tem, imax)
      implicit none
      integer, intent(in) :: imax
      real(8), intent(in) :: z
      real(8), dimension(imax), intent(in) :: tem, an
      real(8), dimension(imax), intent(out) :: tau

      if (z < 0) then
         tau(1:imax) = 2.d0*1.6d-10*(tem(1:imax)*1.d3)**1.5/an(1:imax)
      else
         tau(1:imax) =  7.d-9*(tem(1:imax)*1.d3)**1.5/an(1:imax)
      end if
    end subroutine collision_frequency
  end module erlib_mod

!-------------------------------------------------------------------

!---------------------------------------------------------------------------------------

  module param_mod
    implicit none
    integer, save :: imax
    integer, save :: num_ions
    integer, save :: sw_giota
    real(8), save :: bb
    real(8), save :: rr
    real(8), save :: ra
    integer, save :: sw_gmgs
    integer, save :: transp_model

  contains
    subroutine param_set
      use t3d_interface_tr2er, only : nrmax_tr, bb_tr, rr_tr, ra_tr
      use TRCOMM, only : MDLER
      implicit none
      sw_gmgs = 1
      sw_giota = 0
      num_ions = 1
      imax = nrmax_tr
      bb = bb_tr
      rr = rr_tr
      ra = ra_tr
!---------------------------------------------------------
!    Transport model label
!---------------------------------------------------------
!    Shaing single helisity :: MDLER=4 .or. transp_model=1
!    DCOM/NNW               :: MDLER=5 .or. transp_model=2
!---------------------------------------------------------

      if (MDLER == 4) then
         transp_model = 1
      else if (MDLER == 5) then
         transp_model = 2
      else
         print *, 'MDLER must be 4 or 5.'
         stop
      end if
    end subroutine param_set
  end module param_mod


  module prof_mod
    use param_mod, only : imax, num_ions
    implicit none
    real(8), dimension(:), allocatable :: z, mm
    real(8), dimension(:), allocatable :: rho, er
    real(8), dimension(:,:), allocatable :: an, tem
    
  contains
    subroutine prof_allocate
      implicit none
      if (allocated(z)) then
         continue
      else
         allocate(z(0:num_ions), mm(0:num_ions))
         allocate(an(imax, 0:num_ions), tem(imax, 0:num_ions))
         allocate(rho(imax), er(imax))
      end if
    end subroutine prof_allocate

    subroutine prof_set
      use t3d_interface_tr2er, only : rm_tr, rn_tr, rt_tr, pz_tr, mm_tr
      implicit none
      integer :: i, k
      
      call prof_allocate

      do i = 1, imax
         rho(i) = rm_tr(i)
         an(i,0) = rn_tr(i,1)
         an(i,1) = rn_tr(i,2)
         tem(i,0) = rt_tr(i,1)
         tem(i,1) = rt_tr(i,2)
      end do

      do i = 1, num_ions
         z(i) = pz_tr(i+1)
         mm(i) = mm_tr(i+1)
      end do
      z(0) = -1.0d0
      mm(0) = -1.0d0
      
!@2010.02.01 by wakasa
!      er = 0.0d0
    end subroutine prof_set
  end module prof_mod

  module prof_dn_mod
    use param_mod
    use prof_mod
    implicit none
    real(8) :: beta = 0.000d0
    
! beta interpolation flag 0/1 :: OFF/ON
    integer(4) :: betaip_flg =0
    real(8) :: beta_low, beta_high
  contains
    subroutine prof_dn_set
      implicit none
      integer :: i

      call prof_set
      
    end subroutine prof_dn_set
  end module prof_dn_mod
  
!---------------------------------------------------------------------------------------

!-------------------------------------------------------------------


!----------------------------------------------------------------------------------------

 module ncdn_mod
    use trcomm, only: nrmax,nsmax
    use param_mod
    use prof_dn_mod
    use erlib_mod
    implicit none
    private
!    real(8), dimension(:,:), allocatable, public,save :: dan, dtem
    public :: ncdn_allocate, ncdn_flux, ncdn_local_flux, ncdn_pre_proc, ncdn_analysis
    real(8), dimension(:,:), allocatable, public,save :: dan, dtem
  contains
!-------------------------------------------------------------------
!   subprogram 1.  ncsh_allocate.
!-------------------------------------------------------------------
    subroutine ncdn_allocate
      if (allocated(dan)) then
         continue
      else
         allocate(dan(imax, 0:num_ions), dtem(imax, 0:num_ions))
      end if
    end subroutine ncdn_allocate
!-------------------------------------------------------------------
!   subprogram 2.  ncdn_analysis.
!-------------------------------------------------------------------
    subroutine ncdn_analysis(pflx, hflx, d, v, chi0,chi, kv)
      implicit none
!      real(8), dimension(imax, 0:num_ions), intent(out) :: pflx, hflx, d, v, chi, kv
      real(8), dimension(1:300, 0:10), intent(out) :: pflx, hflx, d, v, chi0,chi, kv
      integer :: i, i0, k, nsw
      real(8) :: local_flx
!      real(8), dimension(imax, 0:num_ions) :: eflx
      real(8), dimension(1:3) :: dd
      real(8), dimension(1:imax, 0:num_ions) :: d2
      
      call ncdn_set_nt
      call prof_set
      call ncdn_pre_proc

      if (rho(1) == 0.0d0) then
         i0 = 2
         pflx(1,0:num_ions) = 0.0d0
         hflx(1,0:num_ions) = 0.0d0
!         eflx(1,0:num_ions) = 0.0d0
         d(1,0:num_ions) = 0.0d0
         d2(1,0:num_ions) = 0.0d0
         chi0(1,0:num_ions) = 0.0d0
         chi(1,0:num_ions) = 0.0d0
      else
         i0 = 1
      end if

      do k = 0, num_ions
         do i = 1, imax
!         do i = i0, imax
            call ncdn_coef(dd, er(i), i, k)
            
            d(i,k) = dd(1)
            d2(i,k)=dd(2)
            chi0(i,k)=dd(3)
            chi(i,k)=dd(3)-1.5*dd(2)
            nsw = 1
            call ncdn_local_flux(local_flx, nsw, er(i), i, k)
            pflx(i,k) = local_flx
            nsw = 2
            call ncdn_local_flux(local_flx, nsw, er(i), i, k)
            hflx(i,k) = local_flx

         end do
      end do
      
 !     do i=1,imax
 !       print'(i6,2e20.10)',i,pflx(i,0),pflx(i,1)
 !     enddo
      
      
      do k = 0, num_ions
         do i = i0, imax
             v(i,k) = (pflx(i,k) + d(i,k)*dan(i,k)) / an(i,k)
             kv(i,k)=(hflx(i,k) + (chi(i,k))*an(i,k)*dtem(i,k)) / (an(i,k)*tem(i,k)) - (1.5d0*pflx(i,k)/an(i,k))
         end do
      end do
!         print'(a80)','=============================================================================' 
!         print'(a4,8a20)','rho','d(i,0)','d(i,1)','v(i,0)','v(i,1)','chi(i,0)','chi(i,1)','kv(i,0)','kv(i,1)'  
!         do i = i0, imax
!             print'(i4,8e20.10)',i,d(i,0),d(i,1),v(i,0),v(i,1),chi(i,0),chi(i,1),kv(i,0),kv(i,1)   
!         enddo
    end subroutine ncdn_analysis
!-------------------------------------------------------------------
!   subprogram 3.  ncdn_flux.
!-------------------------------------------------------------------
    subroutine ncdn_flux(flx, nsw)
      implicit none
      integer, intent(in) :: nsw
      real(8), dimension(imax, 0:num_ions), intent(out) :: flx
      integer :: i, i0, k
      real(8) :: local_flx

      call ncdn_set_nt
      call ncdn_pre_proc

      if (rho(1) == 0.0d0) then
         i0 = 2
         flx(1,0:num_ions) = 0.0d0
      else
         i0 = 1
      end if

      do k = 0, num_ions
         do i = i0, imax
            call ncdn_local_flux(local_flx, nsw, er(i), i, k)
            flx(i,k) = local_flx
         end do
      end do
      
    end subroutine ncdn_flux
!-------------------------------------------------------------------
!   subprogram 4.  ncdn_local_flux.
!-------------------------------------------------------------------
    subroutine ncdn_local_flux(flx, nsw, er0, i, k)
      implicit none
      integer, intent(in) :: nsw
      integer, intent(in) :: i
      integer, intent(in) :: k
      real(8), intent(in) :: er0
      real(8), intent(out) :: flx
      real(8), dimension(3) :: dd
      real(8), dimension(3) :: dd_low,dd_high
      real(8) :: flx_low, flx_high
      integer(4) :: j

      if(betaip_flg .eq.0) then
   
         call ncdn_coef(dd, er0, i, k)
         
         ! nsw = 1 : particle flux
         ! nsw = 2 : energy flux
         select case(nsw)
         case(1)
            flx = - an(i,k)*dd(1)*( &
                 dan(i,k)/an(i,k) &
                 - z(k)*er0/(tem(i,k)*1.0d3) &
                 + (dd(2)/dd(1)-1.5d0)*dtem(i,k)/tem(i,k) &
                 )
!          print'(9e20.10)',flx,an(i,k),dd(1),dan(i,k),z(k),er0,tem(i,k),dd(2),dtem(i,k)
!          print'(9e20.10)',flx,an(i,k),dd(1),dan(i,k),z(k),er0,tem(i,k),dd(2),dtem(i,k)

!            print*,'ip=0',beta,flx
         case(2)
            flx = -an(i,k)*tem(i,k)*dd(2)*( &
                 dan(i,k)/an(i,k) &
                 - z(k)*er0/(tem(i,k)*1.0d3) &
                 +(dd(3)/dd(2)-1.5d0)*dtem(i,k)/tem(i,k) &
                 )
         case default
            print *, 'error in ncsh_flux : nsw'
         end select
      elseif(betaip_flg.eq.1)then
         if(0.00<=beta .and. beta<0.005)then
             beta_low=0.00d0
             beta_high=0.005d0
         elseif(0.005<=beta .and. beta<0.01)then
             beta_low=0.005d0
             beta_high=0.01d0
         elseif(0.01<=beta .and. beta<0.02)then
             beta_low=0.01d0
             beta_high=0.02d0
         elseif(0.02<=beta .and. beta <0.03)then
             beta_low=0.02d0
             beta_high=0.03d0
         else
             print*,'beta interpolation error.'
             stop
         endif
         

         call ncdn_coef_betaip(dd_low, er0, i, k,beta_low)
         call ncdn_coef_betaip(dd_high, er0, i, k,beta_high)
         ! nsw = 1 : particle flux
         ! nsw = 2 : energy flux
         select case(nsw)
         case(1)
            flx_low = - an(i,k)*dd_low(1)*( &
                 dan(i,k)/an(i,k) &
                 - z(k)*er0/(tem(i,k)*1.0d3) &
                 + (dd_low(2)/dd_low(1)-1.5d0)*dtem(i,k)/tem(i,k) &
                 )
            flx_high = - an(i,k)*dd_high(1)*( &
                 dan(i,k)/an(i,k) &
                 - z(k)*er0/(tem(i,k)*1.0d3) &
                 + (dd_high(2)/dd_high(1)-1.5d0)*dtem(i,k)/tem(i,k) &
                 )
            flx = (flx_high -flx_low)/(beta_high -beta_low)*(beta-beta_low) +flx_low
!           print*,'ip=1',beta,flx,flx_low,flx_high,beta_low,beta_high
         case(2)
            flx_low = -an(i,k)*tem(i,k)*dd_low(2)*( &
                 dan(i,k)/an(i,k) &
                 - z(k)*er0/(tem(i,k)*1.0d3) &
                 +(dd_low(3)/dd_low(2)-1.5d0)*dtem(i,k)/tem(i,k) &
                 )
            flx_high = -an(i,k)*tem(i,k)*dd_high(2)*( &
                 dan(i,k)/an(i,k) &
                 - z(k)*er0/(tem(i,k)*1.0d3) &
                 +(dd_high(3)/dd_high(2)-1.5d0)*dtem(i,k)/tem(i,k) &
                 )
            flx = (flx_high -flx_low)/(beta_high -beta_low)*(beta-beta_low) +flx_low
         case default
            print *, 'error in ncsh_flux : nsw'
         end select
      endif

    end subroutine ncdn_local_flux
!-------------------------------------------------------------------
!   subprogram 5.-1  ncdn_set_nt.
!-------------------------------------------------------------------
    subroutine ncdn_set_nt
      use TRCOMM, only : nrmax, rn, rt
      use t3d_interface_tr2er, only : rn_tr, rt_tr        
      implicit none
      integer :: ns

      do ns = 1, 2
         rn_tr(1:nrmax,ns) = rn(1:nrmax,ns)
         rt_tr(1:nrmax,ns) = rt(1:nrmax,ns)
      end do
      
    end subroutine ncdn_set_nt
!-------------------------------------------------------------------
!   subprogram 5-2.  ncdn_pre_proc.
!-------------------------------------------------------------------
    subroutine ncdn_pre_proc
      implicit none
      integer :: k
      real(8), dimension(imax) :: f, df

      do k = 0, num_ions
         f(1:imax) = an(1:imax,k)
         call gradient(df, f, rho, imax)
         dan(1:imax,k) = df(1:imax)/ra
         
         f(1:imax) = tem(1:imax,k)
         call gradient(df, f, rho, imax)
         dtem(1:imax,k) = df(1:imax)/ra
      end do
    end subroutine ncdn_pre_proc
!-------------------------------------------------------------------
!   subprogram 6.  ncsh_coef.

!-------------------------------------------------------------------
    subroutine ncdn_coef(d, er0, i, k)
      implicit none
      real(8), intent(in) :: er0
      integer, intent(in) :: i, k
      real(8), dimension(3), intent(out) :: d

      integer :: is, seed
      real(8), dimension(0:2) :: ZZ
      real(8), dimension(0:2) :: M
      real(8), dimension(0:2) :: TT
      real(8), dimension(0:2) :: NN
      real(8), dimension(0:3) :: DD
      real(8), dimension(0:3) :: DD_low, DD_high

      if (num_ions > 2) then
         print *, 'Error in ncdn_coef : num_ions'
         stop
      end if

      zz(0) = -1.0d0
      zz(1) = 1.0d0
      zz(2) = 2.0d0
      m(0) = -1.0d0
      m(1) = 1.0d0
      m(2) = 4.0d0

      tt(0:2) = 0.0d0
      nn(0:2) = 0.0d0
      do is = 0, num_ions
         tt(is) = tem(i,is)*1.0d3
         nn(is) = an(i,is)*1.0d20
         zz(is) = z(is)
         if (z(is) < 0) then
            m(is) = -1.0
         else
            m(is) = mm(is)
         end if
      end do

      if (zz(k) < 0) then
         seed = 0
      else
         seed = 1
      end if
      
      if(betaip_flg.eq.0) then
         call intrfcDCOMNNW(seed, num_ions, rho(i), er0, rr, &
              beta, bb, tt, nn, zz, m, dd)
         d(1:3) = dd(1:3)
      elseif(betaip_flg.eq.1)then
         if(0.00<=beta .and. beta<0.005)then
             beta_low=0.00d0
             beta_high=0.005d0
         elseif(0.005<=beta .and. beta<0.01)then
             beta_low=0.005d0
             beta_high=0.01d0
         elseif(0.01<=beta .and. beta<0.02)then
             beta_low=0.01d0
             beta_high=0.02d0
         elseif(0.02<=beta .and. beta <0.03)then
             beta_low=0.02d0
             beta_high=0.03d0
         else
             print*,'beta interpolation error.'
             stop
         endif

         call intrfcDCOMNNW(seed, num_ions, rho(i), er0, rr, &
              beta_low, bb, tt, nn, zz, m, dd_low)
         call intrfcDCOMNNW(seed, num_ions, rho(i), er0, rr, &
              beta_high, bb, tt, nn, zz, m, dd_high)
         do is=1,3
         d(is) = (dd_high(is)-dd_low(is))/(beta_high-beta_low)*(beta-beta_low) + dd_low(is)
        enddo
      endif
    end subroutine ncdn_coef

    subroutine ncdn_coef_betaip(d, er0, i, k,betatmp)
      implicit none
      real(8), intent(in) :: er0,betatmp
      integer, intent(in) :: i, k
      real(8), dimension(3), intent(out) :: d

      integer :: is, seed
      real(8), dimension(0:2) :: ZZ
      real(8), dimension(0:2) :: M
      real(8), dimension(0:2) :: TT
      real(8), dimension(0:2) :: NN
      real(8), dimension(0:3) :: DD
      real(8), dimension(0:3) :: DD_low, DD_high


      if (num_ions > 2) then
         print *, 'Error in ncdn_coef : num_ions'
         stop
      end if
      
      zz(0) = -1.0d0
      zz(1) = 1.0d0
      zz(2) = 2.0d0
      m(0) = -1.0d0
      m(1) = 1.0d0
      m(2) = 4.0d0
      
      tt(0:2) = 0.0d0
      nn(0:2) = 0.0d0
      do is = 0, num_ions
         tt(is) = tem(i,is)*1.0d3
         nn(is) = an(i,is)*1.0d20
         zz(is) = z(is)
         if (z(is) < 0) then
            m(is) = -1.0
         else
            m(is) = mm(is)
         end if
      end do

      if (zz(k) < 0) then
         seed = 0
      else
         seed = 1
      end if
       
         call intrfcDCOMNNW(seed, num_ions, rho(i), er0, rr, &
              betatmp, bb, tt, nn, zz, m, dd)
         d(1:3) = dd(1:3)
    end subroutine ncdn_coef_betaip
  end module ncdn_mod

