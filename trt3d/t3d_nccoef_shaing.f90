!
!	t3d_nccoef_shaing.f90
!	task.20110214
!
!	Created by WAKASA Arimitsu on 11/03/02.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!

! based on ER module version beta-1.2(2009/07/24)

!---------------------------------------------------------------------------------------

  module prof_sh_mod
    use param_mod
    use prof_mod

    implicit none
    real(8), dimension(:), allocatable :: eh, et, gmgs
  contains
    subroutine prof_sh_allocate
      implicit none
      if (allocated(eh)) then
         continue
      else
         allocate(eh(imax), et(imax), gmgs(imax))
      end if
    end subroutine prof_sh_allocate
    
    subroutine prof_sh_set
      use t3d_interface_tr2er, only : epsrho_tr, eh_tr
      implicit none
      integer :: i
      
      call prof_set
      
      call prof_sh_allocate
      
      do i = 1, imax
         et(i) = epsrho_tr(i)
         eh(i) = eh_tr(i)
      end do
      gmgs = 1.0d0
      
    end subroutine prof_sh_set
  end module prof_sh_mod

  


!-------------------------------------------------------------------
!  file ncsh.f90.
!  Neoclassical transport model ( Shaing's single helicity model )
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!  code organization
!-------------------------------------------------------------------
!   0. ncsh_mod.
!   1. ncsh_allocate.
!   2. ncsh_analysis.
!   3. ncsh_flux.
!   4. ncsh_local_flux.
!   5. ncsh_pre_proc.
!   6. func_A.
!   7. func_base.
!   8. func_hf.
!   9. func.pf.
!  10. integral_semi_infty.
!-------------------------------------------------------------------
!   subprogram 0.  ncsh_mod.
!-------------------------------------------------------------------
  module ncsh_mod
    use param_mod
    use prof_sh_mod
    use erlib_mod
    implicit none
    private
    integer :: ii, kk
    real(8) :: gEr0
    real(8), dimension(:), allocatable :: deh
    real(8), dimension(:, :), allocatable :: dan, dtem, tau

    public :: ncsh_flux, ncsh_local_flux, ncsh_pre_proc, ncsh_allocate, ncsh_analysis

  contains
!-------------------------------------------------------------------
!   subprogram 1.  ncsh_allocate.
!-------------------------------------------------------------------
    subroutine ncsh_allocate
      implicit none
      if (allocated(deh)) then
         continue
      else
         allocate(deh(imax))
         allocate(dan(imax, 0:num_ions), dtem(imax, 0:num_ions), tau(imax, 0:num_ions))
      end if
    end subroutine ncsh_allocate
!-------------------------------------------------------------------
!   subprogram 2.  ncsh_analysis.
!-------------------------------------------------------------------
    subroutine ncsh_analysis(pflx, hflx, d, v, chi)
      implicit none
      real(8), dimension(imax, 0:num_ions), intent(out) :: pflx, hflx, d, v, chi
      integer :: i, i0, k, nsw
      real(8) :: local_flx
      real(8), dimension(imax, 0:num_ions) :: pflx_n, eflx

      call ncsh_pre_proc

      if (rho(1) == 0.0d0) then
         i0 = 2
         pflx(1,0:num_ions) = 0.0d0
         pflx_n(1,0:num_ions) = 0.0d0
         hflx(1,0:num_ions) = 0.0d0
         eflx(1,0:num_ions) = 0.0d0
      else
         i0 = 1
      end if

      do k = 0, num_ions
         do i = i0, imax
            nsw = 3
            call ncsh_local_flux(local_flx, nsw, er(i), i, k)
            pflx_n(i,k) = local_flx
            nsw = 1
            call ncsh_local_flux(local_flx, nsw, er(i), i, k)
            pflx(i,k) = local_flx
            nsw = 2
            call ncsh_local_flux(local_flx, nsw, er(i), i, k)
            eflx(i,k) = local_flx
         end do
      end do
      
      do k = 0, num_ions
         do i = i0, imax
            d(i,k) = - pflx_n(i,k)/an(i,k)
            v(i,k) = (pflx(i,k) + d(i,k)*dan(i,k)) / an(i,k)
            hflx(i,k) = eflx(i,k) + 1.5d0*pflx(i,k)*tem(i,k)
            chi(i,k) = - hflx(i,k)/(an(i,k)*dtem(i,k))
         end do
      end do
    end subroutine ncsh_analysis
!-------------------------------------------------------------------
!   subprogram 3.  ncsh_flux.
!-------------------------------------------------------------------
    subroutine ncsh_flux(flx, nsw)
      implicit none
      integer, intent(in) :: nsw
      real(8), dimension(imax, 0:num_ions), intent(out) :: flx
      integer :: i, i0, k
      real(8) :: local_flx

      call ncsh_pre_proc

      if (rho(1) == 0.0d0) then
         i0 = 2
         flx(1,0:num_ions) = 0.0d0
      else
         i0 = 1
      end if

      do k = 0, num_ions
         do i = i0, imax
            call ncsh_local_flux(local_flx, nsw, er(i), i, k)
            flx(i,k) = local_flx
         end do
      end do
      
  end subroutine ncsh_flux
!-------------------------------------------------------------------
!   subprogram 4.  ncsh_local_flux.
!-------------------------------------------------------------------
    subroutine ncsh_local_flux(flx, nsw, er0, i, k)
      implicit none
      integer, intent(in) :: nsw
      integer, intent(in) :: i
      integer, intent(in) :: k
      real(8), intent(in) :: er0
      real(8), intent(out) :: flx

      ii = i
      kk = k
      gEr0 = er0

      ! nsw = 1 : particle flux
      ! nsw = 2 : energy flux
      ! nsw = 3 : 
      select case(nsw)
      case(1)
         call integral_semi_infty(flx, func_pf)
      case(2)
         call integral_semi_infty(flx, func_hf)
      case(3)
         call integral_semi_infty(flx, func_base)
      case default
         print *, 'error in ncsh_flux : nsw'
      end select

    end subroutine ncsh_local_flux
!-------------------------------------------------------------------
!   subprogram 5.  ncsh_pre_proc.
!-------------------------------------------------------------------
    subroutine ncsh_pre_proc
      implicit none
      integer :: k
      real(8), dimension(imax) :: f, g, df, dg, tau_tmp

      do k = 0, num_ions
         f(1:imax) = an(1:imax,k)
         call gradient(df, f, rho, imax)
         dan(1:imax,k) = df(1:imax)/ra
         
         f(1:imax) = tem(1:imax,k)
         call gradient(df, f, rho, imax)
         dtem(1:imax,k) = df(1:imax)/ra
         
         f(1:imax) = an(1:imax,k)
         g(1:imax) = tem(1:imax,k)
         call collision_frequency(tau_tmp, z(k), f, g, imax)
         tau(1:imax,k) = tau_tmp(1:imax)
      end do
      
      call gradient(deh, eh, rho, imax)
      deh(1:imax) = deh(1:imax)/ra
    end subroutine ncsh_pre_proc
!-------------------------------------------------------------------
!   subprogram 6.  func_A.
!-------------------------------------------------------------------
    real(8) function func_A(x)
      implicit none
      real(8), intent(in) :: x

      func_A = dan(ii,kk)/an(ii,kk) - z(kk)*gEr0/(tem(ii,kk)*1.d3) + (x-1.5d0)*dtem(ii,kk)/tem(ii,kk)

    end function func_A
!-------------------------------------------------------------------
!   subprogram 7.  func_base.
!-------------------------------------------------------------------
    real(8) function func_base(x)
      implicit none
      real(8), intent(in) :: x
      real(8) :: CNYU, EPS, OMEGAE, OMEGAB, VD, OMEGA2
      
      VD = -tem(ii,kk)*1.d3/et(ii)/bb/rho(ii)/rr/z(kk)
      OMEGAE = -gEr0/et(ii)/bb/rho(ii)/rr
      OMEGAB = - tem(ii,kk)*1.d3*deh(ii)/et(ii)/bb/rho(ii)/rr/z(kk)
      CNYU=1.d0/tau(ii,kk)/eh(ii)

      OMEGA2 = 3.d0*CNYU**2/(x**3)/gmgs(ii) &
           + 1.67d0*et(ii)*(OMEGAE+OMEGAB*x)**2/eh(ii) &
           + (et(ii)/eh(ii))**1.5d0*(OMEGAB*x)**2/4.d0 &
           + 0.6d0*dabs(OMEGAB*x)*CNYU/(x**1.5d0)

      func_base = x**2.5d0*dexp(-x)*CNYU/(x**1.5d0) / OMEGA2
      func_base = func_base*(-et(ii)**2*dsqrt(eh(ii))*vd**2*an(ii,kk))

    end function func_base
!-------------------------------------------------------------------
!   subprogram 8.  func_hf.
!-------------------------------------------------------------------
    real(8) function func_hf(x)
      implicit none
      real(8), intent(in) :: x

      func_hf = tem(ii,kk)*x*func_pf(x)

    end function func_hf
!-------------------------------------------------------------------
!   subprogram 9.  func_pf.
!-------------------------------------------------------------------
    real(8) function func_pf(x)
      implicit none
      real(8), intent(in) :: x

      func_pf = func_base(x)*func_A(x)

    end function func_pf
!-------------------------------------------------------------------
!   subprogram 10. integral_semi_infty
!-------------------------------------------------------------------
    subroutine integral_semi_infty(re1, func)
      implicit none
      real(8), intent(out) :: re1
      real(8), parameter :: errrel=1.0d-7
      integer :: ilst, ind, nm, nmax, nmd, nmmin, np, npd, npmin
      real(8) :: at, atm, atp, csi, csp, ct, es, eps1, h, h0, hn, hs, x1
      real(8) :: func

      H0=0.5D0
      ILST=0
      
      EPS1=errrel**0.75D0
      H=H0
      X1=EXP(-1.D0)
      CSI=2.D0*X1*func(x1)
      RE1=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1
51    IND=0
      ATP=ABS(CSI)
      ATM=ATP 
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

11    IF(IND.NE.1) THEN
         IF(NP.EQ.NPMIN+2) NPD=1
         NP=NP+NPD
         HN=DBLE(NP)*H
         HS=EXP(-HN)
         X1=EXP( HN-HS)
         CT=H*(1.D0+HS)*x1*func(x1)
         RE1=RE1+CT
         AT=ATP
         ATP=ABS(CT)/H
         IF(NP.GE.NPMIN) THEN
         IF(AT+ATP.LE.EPS1*MAX(errrel,ABS(RE1))) THEN
               IF(IND.EQ.-1) GO TO 101
               IND=1
            ENDIF
         ENDIF
      ENDIF

      IF(IND.NE.-1) THEN
         IF(NM.EQ.NMMIN+2) NMD=1
         NM=NM+NMD
         HN=DBLE(NM)*H
         HS=EXP( HN)
         X1=EXP(-HN-HS)
         CT=H*(1.D0+HS)*x1*func(x1)
         re1=re1+CT
         AT=ATM
         ATM=ABS(CT)/H
         IF(NM.GE.NMMIN) THEN
         IF(AT+ATM.LE.EPS1*MAX(errrel,ABS(re1))) THEN
               IF(IND.EQ.1) GO TO 101
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 11

101   ES=ABS(re1-CSP)

      CSP=re1
      IF(ES.LT.EPS1*MAX(errrel,ABS(re1))) GO TO 201
      NMAX=MAX0(NP,NM)
      IF(NMAX.GT.50000) go to 201

      H=0.5D0*H
      re1=0.5D0*re1
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 51

201   continue
    end subroutine integral_semi_infty
  end module ncsh_mod