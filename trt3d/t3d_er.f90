
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------
!  file er.f90.
!  calculate Er from ambipolar condition.
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!  code organization
!-------------------------------------------------------------------
!   0. er_mod.
!   1. er_ambplr.
!   2. local_ambplr
!   3. set_boundary
!   4. zbrent
!   5. func_ambplr
!-------------------------------------------------------------------
!   subprogram 0.  er_mod.
!-------------------------------------------------------------------
  module er_mod
    use param_mod
    use prof_sh_mod
    use ncsh_mod
    use ncdn_mod
    implicit none
    private
    integer :: ii
    public :: er_ambplr,er_dynamics
  contains
!-------------------------------------------------------------------
!   subprogram 1.  er_ambplr.
!-------------------------------------------------------------------
    subroutine er_ambplr(er_i,id)
      implicit none
      integer, intent(in) :: id
      real(8), dimension(imax), intent(in) :: er_i
      integer :: i, i0, ierr
      real(8) :: er0

      select case(transp_model)
      case(1)
         call ncsh_pre_proc
         
      case(2)
         call ncdn_pre_proc
      case default
         print *, 'error in er_ambplr : model'
      end select

      if (rho(1) == 0.0d0) then
         i0 = 2
         er(1) = 0.0d0
      else
         i0 = 1
      end if

      if (id == 0) then
         er0 = 0.0d0
         do i = i0, imax
            call local_ambplr(er0, i, ierr, id)
            if (ierr .ne. 0) then
               print *, 'error in set_boundary'
               stop
            end if
            er(i) = er0
         end do
      else
         do i = i0, imax
            er0 = er_i(i)
            call local_ambplr(er0, i, ierr, id)
            if (ierr == 0) then
               er(i) = er0
            else
               if (id == 1) then
                  call t3d_interpol_p(er(i), rho(i), &
                       rho(i-3), rho(i-2), rho(i-1), &
                       er(i-3), er(i-2), er(i-1))
                  print *, 'error in set_boudary: i = ', i
                  print *, 'er(',i,') was interpolated.'
               end if
            end if
         end do
      end if
!---------------------------------------------------------
!    Transport model label
!---------------------------------------------------------
!    Shaing single helisity :: MDLER=4 .or. transp_model=1
!    DCOM/NNW               :: MDLER=5 .or. transp_model=2
!---------------------------------------------------------
      if (transp_model == 2) then
         if (rho(imax) == 1.0d0) then
            call t3d_interpol_p(er(imax), rho(imax), &
                 rho(imax-3), rho(imax-2), rho(imax-1), &
                 er(imax-3), er(imax-2), er(imax-1))
         end if
      end if

    end subroutine er_ambplr
!-------------------------------------------------------------------
!   subprogram 2.  local_ambplr.
!-------------------------------------------------------------------
    subroutine local_ambplr(er0, i, ierr, id)
      implicit none
      integer, intent(in) :: i, id
      real(8), intent(inout) :: er0
      integer, intent(out) :: ierr
      real(8) :: er1, er2, f1, f2, der
      real(8)::pflxe,pflxi
      integer :: er_itmax
      integer ::Swch_method_AmbiER
      real(8)::initila_er0,initila_er1
      integer :: cntloop=0
#ifdef ErTEST_OP1
!      real(8) :: pflxe,pflxi
#endif

      ii = i
!      der = 500.0d0
      er_itmax=200
      initila_er0=-50000.0d0
      initila_er1=50000.0d0
      der=(initila_er1-initila_er0)/(1.d0*er_itmax)

      if (id == 0) then
         call set_boundary(er1, er2, f1, f2, initila_er0,der, er_itmax,func_ambplr, ierr)
!         print'(a18,i4,4e20.10)','@S.local_ambplr',i,er1,er2,f1,f2
      else
         call set_boundary_multi(er1, er2, f1, f2, er0, der, func_ambplr, ierr)
      end if

!//      Swich of calculation method for detemine the Ambipolar Er// 
!        Swch_method_AmbiER
!        Sato thpe::1
!        Wakasa type:2
        Swch_method_AmbiER=1
      if(Swch_method_AmbiER==1)then
      if (ierr == 0) then
         call zbrent(er0, er1, er2, f1, f2, func_ambplr)
      end if
      else if(Swch_method_AmbiER==2) then
         cntloop=0
         do while(dabs(f1)>1.d-6)
            der=(er2-er1)/(1.d0*er_itmax)
!            print'(a4,e20.10,a8,e20.10,a8,2e20.10)','er1',er1,'<=> er2',er2,'f1&f2',f1,f2
            call set_boundary(er1, er2, f1, f2, er1, der,er_itmax,func_ambplr, ierr)
            er0=er1
            cntloop=cntloop+1
            if(cntloop>10)then
                print*,'loop over',f1,f2
                exit
            endif
         enddo
      else
        print*,'Illegal number is set for =Swch_method_AmbiER=. '
        stop
      endif
     
!#ifdef ErTEST_OP1

      call ncdn_local_flux(pflxe,1,er0,ii,0)
      call ncdn_local_flux(pflxi,1,er0,ii,1)
!      print'(a10,i5,3e20.10)','TEST',ii,er0,pflxe,pflxi
!#endif
      
    end subroutine local_ambplr

!-------------------------------------------------------------------
!   subprogram 3.  set_boundary.
!-------------------------------------------------------------------      
    subroutine set_boundary(x1, x2, fa, fb, x0,dx0,itmax,Func, ierr)

      implicit none
      real(8), intent(out) :: x1
      real(8), intent(out) :: x2
      real(8), intent(out) :: fa
      real(8), intent(out) :: fb
      integer, intent(out) :: ierr
      integer, intent(in) :: itMax
      
      real(8), intent(in) :: x0
      real(8), intent(in) :: dx0
      real(8) :: Func
      
      integer :: nsw
      integer :: i
!      integer, parameter :: itMax = 200
      real(8) :: dx
      integer :: icntEr
      real(8),dimension(1:3) :: fa0,fb0,xx
      integer(4)::iflg_eroot

!!      iflg_eroot 0/1 :: i-root/e-root
      iflg_eroot=1
      dx=dx0
      ierr=1      
      icntEr=0
      
      if(iflg_eroot==1)then
         nsw=-1
         x1=x0+(1.d0*itmax)*dx
      else
         nsw=1
         x1=x0
      endif

      fa=Func(x1)

      do i=1,itMax
         x2=x1+nsw*dx
         fb=Func(x2)
!         print'(a20,4e20.10)','in S.set_boundar',x1,x2,fa,fb
!         print*,fb
         if(fa*fb<=0.0)then
            icntEr=icntEr+1
            if(icntEr>=4)then
              print*,'NUM Er root Error!'
              stop
            endif
            if(nsw<0)then
                xx(icntEr)=x2
                fa0(icntEr)=fb
                fb0(icntEr)=fa
            else
                xx(icntEr)=x1
                fa0(icntEr)=fa
                fb0(icntEr)=fb
            endif
            ierr=0
            exit
         endif
         x1=x2
         fa=fb
      enddo
      x1=xx(1)
      x2=xx(1)+dx
      fa=fa0(1)
      fb=fb0(1)

!2010.02.15      x=100000
!2010.02.15!       x=x0
!2010.02.15      fa = Func(x)
!2010.02.15      ierr = 1
!2010.02.15      if (fa <= 0.0d0) then
!2010.02.15!         x = x0
!2010.02.15         dx = -dx0
!2010.02.15         nsw = -1
!2010.02.15      else
!2010.02.15!         x = x0
!2010.02.15         dx = dx0
!2010.02.15         nsw = 1
!2010.02.15      end if
!2010.02.15      
!2010.02.15      do it = 1, itMax
!2010.02.15         x = x + dx
!2010.02.15         fb = Func(x)
!2010.02.15!         print('(i6,3e20.10)'),it,x,fa,fb
!2010.02.15         if (fa*fb <= 0.0d0) then
!2010.02.15            ierr = 0
!2010.02.15            exit
!2010.02.15         end if
!2010.02.15      end do
!2010.02.15
!2010.02.15      if ( nsw == 1 ) then
!2010.02.15         x2 = x
!2010.02.15         x1 = x2-dx0 
!2010.02.15!         x1 = x0
!2010.02.15      else
!2010.02.15!         x2 = x0
!2010.02.15         x1 = x
!2010.02.15         x2=  x1+dx0
!2010.02.15         x = fa
!2010.02.15         fa = fb
!2010.02.15         fb = x
!2010.02.15      end if
      
!      if (fa < 0.0d0) then
!         print *, x1, x2
!         print *, 'fa < 0'
!         stop
!      end if
    end subroutine set_boundary

    subroutine set_boundary_81079(x1, x2, fa, fb, x0, dx0, Func,inr, ierr)

      implicit none
      real(8), intent(out) :: x1
      real(8), intent(out) :: x2
      real(8), intent(out) :: fa
      real(8), intent(out) :: fb
      integer, intent(out) :: ierr
      
      real(8), intent(in) :: x0
      real(8), intent(in) :: dx0
      integer(4),intent(in)::inr
      real(8) :: Func
      
      integer :: nsw
      integer :: i
      integer, parameter :: itMax = 200
      real(8) :: dx, x
      integer :: icntEr
      real(8),dimension(1:3) :: fa0,fb0,xx
      integer(4)::iflg_eroot   

!!      iflg_eroot 0/1 :: i-root/e-root
      
      if(inr.gt.0)then
      iflg_eroot=1
      else
      iflg_eroot=0
      endif
      
      dx=dx0 
      ierr=1      
      icntEr=0
      
      if(iflg_eroot==1)then
         nsw=1
      else
         nsw=-1
      endif

      do i=0,itMax
         x1=nsw*50000.0d0-nsw*i*dx
         x2=x1-nsw*dx
!         print*,x1,x2 
         fa=Func(x1)
!         print*,fa
         fb=Func(x2)
!         print*,fb
!         print*,x1,x2
         if(fa*fb<=0.0)then
            icntEr=icntEr+1
            if(icntEr>=4)then
              print*,'NUM Er root Error!'
              stop
            endif
            if(nsw>0)then
            xx(icntEr)=x2
            fa0(icntEr)=fb
            fb0(icntEr)=fa
            else
            xx(icntEr)=x1
            fa0(icntEr)=fa
            fb0(icntEr)=fb
            endif
            ierr=0
            exit
         endif
      enddo
      x1=xx(1)
      x2=xx(1)+dx
      fa=fa0(1)
      fb=fb0(1)
    end subroutine set_boundary_81079


!-------------------------------------------------------------------
!   subprogram 3.  set_boundary.
!-------------------------------------------------------------------      
    subroutine set_boundary_multi(x1, x2, fa, fb, x0, dx0, Func, ierr)
      implicit none
      real(8), intent(out) :: x1
      real(8), intent(out) :: x2
      real(8), intent(out) :: fa
      real(8), intent(out) :: fb
      integer, intent(out) :: ierr
      
      real(8), intent(in) :: x0
      real(8), intent(in) :: dx0
      real(8) :: Func

      integer :: it
      integer, parameter :: itMax = 1000
      real(8) :: x

      fa = Func(x0)
      ierr = 1
      x = x0
      do it = 1, itMax
         x = x + dx0
         fb = Func(x)
         if (fa*fb <= 0.0d0) then
            ierr = 0
            x1 = x - dx0
            x2 = x
            exit
         end if
      end do

    end subroutine set_boundary_multi

!-------------------------------------------------------------------
!   subprogram 4.  zbrent.
!-------------------------------------------------------------------
    subroutine zbrent(x, x1, x2, fa, fb, Func)
      implicit none
      real(8), intent(in) :: x1, x2
      real(8), intent(inout) :: fa, fb
      real(8) :: Func
      real(8), intent(out) :: x

      integer :: itMax
      real(8) :: tol, eps
      
      parameter(itMax=100, eps = 3.0d-12, tol = 1.0d-12)
      integer :: iter
      real(8) :: a, b, c, d, e, fc, p, q, r, s, tol1, xm
      
!     er1 と er2 の間に解があるかどうか判定
!     解がある場合，fa(er1)*fb(er2) < 0  である.
      if ( (fa.gt.0.0d0 .and. fb.gt.0.0d0) .or. &
           (fa.lt.0.0d0 .and. fb.lt.0.0d0) ) then
         print *, 'ZBrent: error'
         stop
      end if
      
      a = x1
      b = x2
      c = b
      fc = fb
      
      do iter = 1, itMax
         if ( (fb.gt.0.0d0 .and. fc.gt.0.0d0) .or. &
              (fb.lt.0.0d0 .and. fc.lt.0.0d0) ) then
            c = a
            fc = fa
            d = b-a
            e = d
         end if
         if (abs(fc).lt.abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         end if
         tol1 = 2.0d0*eps*abs(b) + 0.5d0*tol
         xm = 0.5d0*(c-b)
         if (abs(xm).le.tol1 .or. fb.eq.0.0d0) then
            x = b
            return
         end if
         
         if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
            s = fb / fa
            if (a.eq.c) then
               p = 2.0d0*xm*s
               q = 1.0d0 - s
            else
               q = fa / fc
               r = fb / fc
               p = s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
               q = (q-1.0d0)*(r-1.0d0)*(s-1.0d0)
            end if
            if (p.gt.0.0d0) q = -q
            p = abs(p)
            if (2.0d0*p.lt.min(3.0d0*xm*q-abs(tol1*q),ABS(e*q))) then
               e = d
               d = p / q
            else
               d = xm
               e = d
            end if
         else
            d = xm
            e = d
         end if
         a = b
         fa = fb
         if (abs(d) .gt. tol1) then
            b = b + d
         else
            b = b + sign(tol1, xm)
         end if
         fb = Func(b)
      end do
      print *, 'convergence error'
      stop
      x = b
      return
      
    end subroutine zbrent

!-------------------------------------------------------------------
!   subprogram 5.  func_ambplr.
!-------------------------------------------------------------------
    real(8) function func_ambplr(er0)
      use TRCOMM, only:RG
      implicit none
      real(8), intent(in) :: er0

      integer :: k, nsw
      real(8) :: pflx

#ifdef ErTEST_OP1
!      integer ::i
!      real(8), dimension(0:2)::pflx_tmp 
#endif
      nsw = 1
      func_ambplr = 0.0d0
      select case(transp_model)
      case(1)
         do k = 0, num_ions
            call ncsh_local_flux(pflx, nsw, er0, ii, k)
            func_ambplr = func_ambplr + z(k)*pflx
         end do

         func_ambplr = - func_ambplr
      case(2)
         do k = 0, num_ions
            call ncdn_local_flux(pflx, nsw, er0, ii, k)
!            print*,'TESTTESTTEST',k,er0,z(k)*pflx
#ifdef ER_adjust_32940
            if(k==0 .and. rg(ii)>0.7)then
                print*, "set ER_adjust_32940 parameters."
                print*,"Are you OK?"
                stop
                print'(a13,i4,3e14.6)','FUNC_ER :: ',ii,rg(ii),pflx,(0.61391 + 3.4969*rg(ii)**2 -16.553*rg(ii)**4 +16.62*rg(ii)**6)*1.d-1
                print*,'NOW!'
                pflx=pflx+ (2.607 -8.2819*rg(ii)**2 +2.4616*rg(ii)**4 +7.3509*rg(ii)**6)*1.d-1
            endif
#endif
            
            func_ambplr = func_ambplr + z(k)*pflx
#ifdef ErTEST_OP1
!            pflx_tmp(k)=z(k)*pflx
#endif
         end do
         func_ambplr = - func_ambplr
#ifdef ErTEST_OP1
!        print*,'+'
!        print'(4e20.10)',er0,pflx_tmp(0),pflx_tmp(1),func_ambplr
!        print*,'='
#endif
      case default
         print *, 'error in func_ambplr : transp_model'
      end select
    end function func_ambplr

!--------------------------------------------------------------------------------------
! suboutine 6. er_dynamics
!--------------------------------------------------------------------------------------
!     subroutine er_dynamics(t,er)
     subroutine er_dynamics(t,er)
     use TRCOMM, only:DT,DR,PI,VC,AEE,QP,NRMAX,RG,RN
     use t3d_er_param
     implicit none
     integer, intent(in):: t
     real(8), dimension(nrmax), intent(inout):: er
     real*8,dimension(3,nrmax-2)::a
     real*8,dimension(nrmax)::bs,bc,pdea,g,pfe,pfi, per,pper,ep_perp
     real*8,dimension(nrmax)::coe_aa,coe_bb,coe_cc,coe_dd,er0
     real(8),dimension(1:nrmax,1:nrmax) :: Matrix_LHS, inv_MATRIX_LHS
     real*8,dimension(nrmax-2)::b1,er2
     real*8:: coe,dt_dr2
     real*8:: flx,flxe,flxi
     integer :: nr,k,kmax, ierr,i,j
     real(8):: eri1,eri2,eri3,eri4,eri5,eri6,ak1,ak2,ak3,ak4,ak5,ak6
     real(8):: pfe1, pfe2, pfe3, pfe4, pfe5, pfe6, pfi1, pfi2, pfi3, pfi4, pfi5, pfi6
     real(8),dimension(NRMAX) :: PP, QQ
     integer(4), save :: flg_log_TODA=0
#ifdef er_dynamics_log_20100716
     integer(4), save :: flg_log=0
     integer(4) :: NUM_OUT_LOG
     real*8,dimension(nrmax) :: bs_stock, bc_stock, er_stock,g_stock,pfe_stock, pfi_stock, c2va2, ge_minus_gi
     real*8,dimension(nrmax):: bc_1st, bc_2nd,bc_3rd, denominator, bc_2nd_1, bc_2nd_2, bc_2nd_3
     real*8,dimension(nrmax):: per_stock, pper_stock
#endif
#ifdef er_dynamics_log_QQ
     integer(4), save :: flg_log_QQ=0, icnt
#endif
     integer(4):: scheme_switch
! for case(4): Balance between Particle flux and Diffusion of Er.
     real(8)    :: delta_er, epsilon_er, sum_error_balance, er_change_factor, pre_sum_error_balance
	 real(8), dimension(1:nrmax)::error_balance, del_error_balance, delGAMMA, delDIFF, er1
	 real(8), dimension(1:nrmax)::del_pfe, del_pfi, del_coe_aa, del_coe_bb
	 integer(4) :: cnt_er_iteration,MAX_er_iteration
! for case(5): TEST Routine.
     real(8) :: dea_const
! fpr case(37): iteration
     integer(4) :: it_cnt,it_cnt_linesearch, it_cnt_max
     real(8):: gamma_flx_e_0, gamma_flx_e_1, gamma_flx_i_0, gamma_flx_i_1, der_flx_e, der_flx_i
     real(8) :: SUM_Error_Func, NEW_SUM_Error_Func, accelerator_factor,accelerator_factor_1,accelerator_factor_2,accelerator_factor_3
     real(8) :: SUM_Error_Func_1=0.d0,SUM_Error_Func_2=0.d0,SUM_Error_Func_3=0.d0
     real(8) :: threshold_relative_error_SUM=0.d0,threshold_NEW_error_SUM=0.d0, relative_error_SUM=0.d0
     real(8) :: MAX_grad_SUM_Error_Func=0.d0
     real(8) :: SUM_Error_Func_tmp=0.d0, accelerator_factor_tmp=0.d0
     
     real(8) ::Error_Func(nrmax), Error_Func_1st(nrmax), Error_Func_2nd(nrmax), grad_SUM_Error_Func(nrmax),er_dr2(nrmax)
     real(8) :: er_0(nrmax), er_1(nrmax), er_2(nrmax), er_3(nrmax), er_tmp(nrmax)
     real(8) :: pfe_1(nrmax),pfi_1(nrmax),er_1_dr2(nrmax),ERROR_FUNC_1(nrmax)
     real(8) :: pfe_2(nrmax),pfi_2(nrmax),er_2_dr2(nrmax),ERROR_FUNC_2(nrmax)
     real(8) :: pfe_3(nrmax),pfi_3(nrmax),er_3_dr2(nrmax),ERROR_FUNC_3(nrmax)
     real(8) :: pfe_tmp(nrmax), pfi_tmp(nrmax),er_tmp_dr2(nrmax),Error_Func_tmp(nrmax)
     real(8) :: dea_dr(nrmax), dea_drder(nrmax), dea_der(nrmax), er_dr(nrmax), er_1_dr(nrmax), er_2_dr(nrmax), er_3_dr(nrmax), er_tmp_dr(nrmax)
     real(8) :: Error_Func_der_1st(nrmax), Error_Func_der_2nd(nrmax)
     real(8) :: cpu_time_start, cpu_time_end
     real(8), dimension(1:nrmax) :: term_alpha, term_beta, term_gamma, term_a2b2g2, term_ab_bg,term_ag, func_gamma, func_diff, dfunc_gamma_dEr
     real(8), dimension(0:nrmax+1) :: er_for_tmp
     real(8) :: EF_normalized_factor=0.d0     
    
    do nr=1,nrmax
        pfe(nr)=0.d0
        pfi(nr)=0.d0
        g(nr)=0.0d0
        Error_Func(nr)=0.0d0
        Error_Func_1st(nr)=0.0d0
		Error_Func_2nd(nr)=0.0d0
		grad_SUM_Error_Func(nr)=0.0d0
		er_dr2(nr)=0.0d0
		er_0(nr)=0.0d0
		er_1(nr)=0.0d0
		er_2(nr)=0.0d0
		er_3(nr)=0.0d0
		er_tmp(nr)=0.0d0
		pfe_1(nr)=0.0d0
		pfi_1(nr)=0.0d0
		er_1_dr2(nr)=0.0d0
		ERROR_FUNC_1(nr)=0.0d0
		pfe_2(nr)=0.0d0
		pfi_2(nr)=0.0d0
		er_2_dr2(nr)=0.0d0
		ERROR_FUNC_2(nr)=0.0d0
		pfe_3(nr)=0.0d0
		pfi_3(nr)=0.0d0
		er_3_dr2(nr)=0.0d0
		ERROR_FUNC_3(nr)=0.0d0
		pfe_tmp(nr)=0.0d0
		pfi_tmp(nr)=0.0d0
		er_tmp_dr2(nr)=0.0d0
		Error_Func_tmp(nr)=0.0d0
        dea_drder(nr)=0.0d0
		dea_dr(nr)=0.0d0
        dea_der(nr)=0.d0
		er_dr(nr)=0.0d0
		er_1_dr(nr)=0.0d0
		er_2_dr(nr)=0.0d0
		er_3_dr(nr)=0.0d0
		er_tmp_dr(nr)=0.0d0 
        Error_Func_der_1st(nr)=0.d0
        Error_Func_der_2nd(nr)=0.d0
        term_alpha=0.d0
        term_beta=0.d0
        term_gamma=0.d0
        term_a2b2g2=0.d0
        term_ab_bg=0.d0
        term_ag=0.d0
        func_gamma=0.d0
        func_diff=0.d0
        dfunc_gamma_dEr=0.d0
        er_for_tmp=0.d0
      enddo
     er_for_tmp(0)=0.0d0
     er_for_tmp(nrmax+1)=0.0d0
     er0(1:nrmax)=er(1:nrmax)
     
     
#ifdef er_dynamics_log_QQ
     if(flG_log_QQ.eq.0)then
     open(812,file='er_dynamics_log_QQ.txt',status='replace')
     else
     open(812,file='er_dynamics_log_QQ.txt',position='append')
     endif     
#endif
#ifdef er_dynamics_log_20100716
     if(flG_log_QQ.eq.0)then
     open(820,file='er_dynamics_log.txt',status='replace')
     else
     open(820,file='er_dynamics_log.txt',position='append')
     endif     
#endif
!   DR=Ra/nrmax, QP: safety factor, BB: toroidal magnetic field, vc: speed of light,
!   r: radial coordinate, AEE: electron charge, kmax=tstep(dynamics of n,T)/dt
!   To avoid the Zero divide, multiply 1e10. 

    if(flg_log_TODA==0)then
    open(830,file='er_dynamics_log_TODA.txt',status='replace')
    else
    open(830,file='er_dynamics_log_TODA.txt',position='append')
    endif     
    

    call t3d_er_param_initialize
!    scheme_switch = 3 ! 1:Euler explicit / 2:Crank-Nicolson / 3:Euler implicit / 4: Balance (dEr/dt=0) / 5: TEST routine
    if(30 > ER_CALC_METHOD .or. ER_CALC_METHOD >= 40) then
        print*,'=============================================================================='
        print*,'ER_CALC_METHOD ERROR, you must select 30 <= ER_CALC_METHOD <40 in er_dynamics.'
        print*,'Now you selected ER =', ER_CALC_METHOD
        print*
        stop
    endif

    print*,'==================  er_dynamics parameter  ================'
    print*,'DEA(1)=',dea(1)
    print*,'GAMMA=',er_gamma_factor
    print*,'DT_ER=',DT_ER
    print*

    select case(ER_CALC_METHOD)
    
    case(31) ! Euler explicit 
    kmax=int(1.d0/dt_er)/int(1.d0/DT)
#ifdef er_dynamics_log_20100716    
    NUM_OUT_LOG=kmax
#endif    
!    kmax=int(1.d0/dt_er)/int(1.d0/DT)*500

    do k=1,kmax
    if(mod(k,100)==0 .or. k==1)then
!      print'(i5,a1,i5)',k,'/',kmax
!      do nr=1,nrmax
!      write(830,'(2i6,2e20.10)') k,nr,rg(nr),er(nr)
!      print '(2i6,2e20.10)' ,k,nr,rg(nr),er(nr)
!      enddo
      flg_log_TODA=1
    endif
      
  
    do nr=1,nrmax
!    er(1)=0d0
!    er(nrmax)=12.0d0
    enddo

    do nr=2,nrmax-1
    pdea(nr)=(dea(nr+1)-dea(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    pdea(1)=(dea(2)-dea(1))/(RA*(RG(2)-RG(1)))
    pdea(nrmax)=(dea(nrmax)-dea(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))
 
    do nr=3,nrmax-2
    per(nr)=(er(nr-2)-8.d0*er(nr-1)+8.d0*er(nr+1)-er(nr+2))/(RA*(RG(nr)-RG(nr-1)))/12.d0
    enddo
    per(1)=(-25.d0*er(1)+48.d0*er(2)-36.d0*er(3)+16.d0*er(4)-3.d0*er(5))/(RA*(RG(2)-RG(1)))/12.d0
    per(2)=(-25.d0*er(2)+48.d0*er(3)-36.d0*er(4)+16.d0*er(5)-3.d0*er(6))/(RA*(RG(3)-RG(2)))/12.d0
    per(nrmax-1)=(er(nrmax-3)-4.d0*er(nrmax-2)+3.d0*er(nrmax-1))/(RA*(RG(nrmax-1)-RG(nrmax-2)))/2.d0
    per(nrmax)=(er(nrmax-2)-4.d0*er(nrmax-1)+3.d0*er(nrmax))/(RA*(RG(nrmax)-RG(nrmax-1)))/2.d0

    do nr=3,nrmax-2
    pper(nr)=(-er(nr-2)+16.d0*er(nr-1)-30.d0*er(nr)+16.d0*er(nr+1)-er(nr+2))/((RA*(RG(nr)-RG(nr-1)))*(RA*(RG(nr)-RG(nr-1))))/12.d0
    enddo
    pper(1)=(2.d0*er(1)-5.d0*er(2)+4.d0*er(3)-er(4))/((RA*(RG(2)-RG(1)))*(RA*(RG(2)-RG(1))))
    pper(2)=(2.d0*er(2)-5.d0*er(3)+4.d0*er(4)-er(5))/((RA*(RG(3)-RG(2)))*(RA*(RG(3)-RG(2))))
    pper(nrmax-1)=(er(nrmax-2)-2.d0*er(nrmax-1)+er(nrmax))/((RA*(RG(nrmax-1)-RG(nrmax-2)))*(RA*(RG(nrmax-1)-RG(nrmax-2))))
    pper(nrmax)=(er(nrmax-2)-2.d0*er(nrmax-1)+er(nrmax))/((RA*(RG(nrmax)-RG(nrmax-1)))*(RA*(RG(nrmax)-RG(nrmax-1))))
  
    call ncdn_pre_proc
    do nr=1,nrmax
     g(nr)= 8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
    enddo
    do nr=2,nrmax-1
     call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
     pfe(nr)=flxe*1.0d+20
     call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
     pfi(nr)=flxi*1.0d+20
     bs(nr)=(1.d0/6.d0/dt_er-dea(nr)/2.d0/DR**2)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     bc(nr)=(er_gamma_factor*(pfe(nr)-pfi(nr))/g(nr)+dea(nr)*pper(nr)+dea(nr)*per(nr)/(rg(nr)*ra)+(er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)*(dt_er)
#ifdef er_dynamics_log_20100716
     if(mod(k,NUM_OUT_LOG).eq.0)then
     per_stock(nr)=per(nr)
     pper_stock(nr)=pper(nr)
     bs_stock(nr)=bs(nr)
     bc_stock(nr)=bc(nr)
     er_stock(nr)=er(nr)
     g_stock(nr)=g(nr)
     pfe_stock(nr)=pfe(nr)
     pfi_stock(nr)=pfi(nr)
     c2va2(nr)=vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2
     ge_minus_gi(nr)=(pfe(nr)-pfi(nr))
     bc_1st(nr)=(er_gamma_factor*(pfe(nr)-pfi(nr))/g(nr))*dt_er
     bc_2nd(nr)=(dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr))*dt_er
     bc_3rd(nr)=((er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)*dt_er
     denominator(nr)=(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     bc_2nd_1(nr)=dea(nr)*pper(nr)*dt_er
     bc_2nd_2(nr)=dea(nr)*per(nr)/(rg(nr)*ra)*dt_er
     bc_2nd_3(nr)=pdea(nr)*per(nr)*dt_er
     endif
#endif
    enddo
    
#ifdef er_dynamics_log_20100716
     er_stock(1)=er(1)
     er_stock(nrmax)=er(nrmax)
#endif
       do nr=2,nrmax-1
       er(nr)=bc(nr)
       enddo
       er(1)=0.d0
       er(nrmax)=per(nrmax-1)*(rg(nrmax)*ra-rg(nrmax-1)*ra)+er(nrmax-1)

    select case(ER_CALC_BOUNDARY_CONDITION)
    case(0)
    coe=g(nrmax)
!
	do nr=nrmax,nrmax
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 0)
	pfe1=flx*1.d20
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 1)
	pfi1=flx*1.d20
	ak1=dt_er*(pfe1-pfi1)/coe
	eri1=er(nr)+ak1/4.d0
	call ncdn_local_flux(flx, 1, eri1, nr, 0)
	pfe2=flx*1.d20
	call ncdn_local_flux(flx, 1, eri1, nr, 1)
	pfi2=flx*1.d20
	ak2=dt_er*(pfe2-pfi2)/coe
	eri2=er(nr)+3.d0*ak1/3.2d1+9.d0*ak2/3.2d1
	call ncdn_local_flux(flx, 1, eri2, nr, 0)
	pfe3=flx*1.d20
	call ncdn_local_flux(flx, 1, eri2, nr, 1)
	pfi3=flx*1.d20
	ak3=dt_er*(pfe3-pfi3)/coe
	eri3=er(nr)+1.932d3*ak1/2.197d3-7.2d3*ak2/2.197d3+7.296d3*ak3/2.197d3
	call ncdn_local_flux(flx, 1, eri3, nr, 0)
	pfe4=flx*1.d20
	call ncdn_local_flux(flx, 1, eri3, nr, 1)
	pfi4=flx*1.d20
	ak4=dt_er*(pfe4-pfi4)/coe
	eri4=er(nr)+4.39d2*ak1/2.16d2-8.d0*ak2+3.68d3*ak3/5.13d2-8.45d2*ak4/4.104d3
	call ncdn_local_flux(flx, 1, eri4, nr, 0)
	pfe5=flx*1.d20
	call ncdn_local_flux(flx, 1, eri4, nr, 1)
	pfi5=flx*1.d20
	ak5=dt_er*(pfe5-pfi5)/coe
	eri5=er(nr)-8.d0*ak1/2.7d1+2.d0*ak2-3.544d3*ak3/2.565d3+1.859d3*ak4/4.104d3-1.1d1*ak5/4.d1
	call ncdn_local_flux(flx, 1, eri5, nr, 0)
	pfe6=flx*1.d20
	call ncdn_local_flux(flx, 1, eri5, nr, 1)
	pfi6=flx*1.d20
	ak6=dt_er*(pfe6-pfi6)/coe
	er(nr)=er(nr)+1.6d1*ak1/1.35d2+6.656d3*ak3/1.2825d4+2.8561d4*ak4/5.643d4-9.d0*ak5/5.d1+2.d0*ak6/5.5d1
       enddo
    case(1)
      er(nrmax)=er_boundary_value
    case default
      print*,'ER BOUNDARY CONDITION ERROR. You must select ER_CALC_BOUNDARY_CONDITION = 0 .or. 1.'
      print*,'stopped @ ',__FILE__,__LINE__
      stop
    end select

#ifdef er_dynamics_log_20100716
!      if(k==1 .or. k==kmax) then
      write(820,'(a1,a19,22a20)')'#','1:time','2:rho','3:r','4:er(t)','5:g','6:bc','7:bs','8:c^2/vA^2','9:bc_1st','10:bc_2nd','11:bc_3rd','12:denominator','13:Ge-Gi','14:bc_2nd_1','15:bc_2nd_2','16:bc_2nd_3','17:pfe','18:pfi','19:er(t+1)','20:per','21:pper','22:pfe','23:pfi'
      do nr=1,nrmax
      write(820,'(23e20.10)') t*DT+dt_er*k,RG(nr),RG(nr)*RA,er_stock(nr),g_stock(nr),bc_stock(nr),bs_stock(nr),c2va2(nr),bc_1st(nr),bc_2nd(nr),bc_3rd(nr),denominator(nr),ge_minus_gi(nr),bc_2nd_1(nr),bc_2nd_2(nr),bc_2nd_3(nr),pfe_stock(nr),pfi_stock(nr),er(nr),per_stock(nr),pper_stock(nr),pfe_stock(nr),pfi_stock(nr)
      enddo
      write(820,*)
 !     endif
      flg_log=1   
#endif
#ifdef er_dynamics_log_QQ
       icnt=1
       if(flg_log_QQ==0)then
       write(812,'(a5,5a20)') '# i  ','RG','time','Er','Qe','Qi'
       flg_log_QQ=1
       endif
       if(k==kmax)then 
       do i=1,nrmax
       call ncdn_local_flux(flxe, 2, er(i), i, 0)
       call ncdn_local_flux(flxi, 2, er(i), i, 1)
       write(812,'(i5,5e20.10)') i,rg(i),t*DT+dt_er*k,er(i),flxe*1d20,flxi*1d20
       enddo
       write(812,*)
       endif
#endif
       enddo

    case(32) ! Euler implicit

!   TDMA: Tro-Diagonal Matrix Algorithm
!   | coe_bb  -coe_cc                           | |er|  |coe_dd|
!   |-coe_aa   coe_bb  -coe_cc                  | |  |  |coe_dd|
!   |         -coe_aa   coe_bb  -coe_cc         | |  |= |coe_dd|
!   |                   ......                  | |  |  |coe_dd|
!   |                           -coe_aa   coe_bb| |  |  |coe_dd|
!
!   PP[NRMAX], QQ[NEMAX] :: 係数格納用
    kmax=int(1.d0/dt_er)/int(1.d0/DT)
#ifdef er_dynamics_log_20100716    
    NUM_OUT_LOG=kmax
#endif    
!    kmax=int(1.d0/dt_er)/int(1.d0/DT)*500

    do k=1,kmax
    print*,'k/kmax',k,'/',kmax,':DT',DT
    if(mod(k,100)==0 .or. k==1)then
!      print'(i5,a1,i5)',k,'/',kmax
!      do nr=1,nrmax
!      write(830,'(2i6,2e20.10)') k,nr,rg(nr),er(nr)
!      print '(2i6,2e20.10)' ,k,nr,rg(nr),er(nr)
!      enddo
      flg_log_TODA=1
    endif

!    do nr=2,nrmax-1
!     pdea(nr)=(dea(nr+1)-dea(nr-1))/(2.d0*DR)
!    enddo


    call ncdn_pre_proc
    do nr=1,nrmax
!      g(nr)=(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
      g(nr)= 8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
!      coe_aa(nr)=(1.d0/dt_er)-2.d0*dea(nr)/DR**2
!      coe_bb(nr)=(1.d0/DR**2 + 1.d0/(2.d0*rg(nr)*ra*DR))*dea(nr)
!      coe_cc(nr)=(1.d0/DR**2 - 1.d0/(2.d0*rg(nr)*ra*DR))*dea(nr)
 
      dt_dr2=dt_er/(DR*DR)           
      coe_cc(nr)=dt_dr2*dea(nr)
      coe_aa(nr)=1.d0+2.d0*dt_dr2*dea(nr)
      coe_bb(nr)=dt_dr2*dea(nr)      
  
      

      call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
      pfe(nr)=flxe*1.0d+20
      call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
      pfi(nr)=flxi*1.0d+20
      coe_dd(nr)=er(nr) + dt_er/g(nr)*(pfe(nr)-pfi(nr))*er_gamma_factor
    enddo

    PP(1)=coe_bb(1)/coe_aa(1)
    QQ(1)=coe_dd(1)/coe_aa(1)
    do nr=2,nrmax
    PP(nr)=coe_bb(nr)/(coe_aa(nr)-coe_cc(nr)*PP(nr-1))
    QQ(nr)=(coe_dd(nr)+coe_cc(nr)*QQ(nr-1))/(coe_aa(nr)-coe_cc(nr)*PP(nr-1))
    enddo
    er(nrmax)=QQ(nrmax)
    do nr=nrmax-1,1,-1
    er(nr)=PP(nr)*er(nr+1)+QQ(nr)  
    enddo

  
    er(1)=0.e0
    
    select case(ER_CALC_BOUNDARY_CONDITION)
    case(0)
    coe=g(nrmax)
!
	do nr=nrmax,nrmax
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 0)
	pfe1=flx*1.d20
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 1)
	pfi1=flx*1.d20
	ak1=dt_er*(pfe1-pfi1)/coe
	eri1=er(nr)+ak1/4.d0
	call ncdn_local_flux(flx, 1, eri1, nr, 0)
	pfe2=flx*1.d20
	call ncdn_local_flux(flx, 1, eri1, nr, 1)
	pfi2=flx*1.d20
	ak2=dt_er*(pfe2-pfi2)/coe
	eri2=er(nr)+3.d0*ak1/3.2d1+9.d0*ak2/3.2d1
	call ncdn_local_flux(flx, 1, eri2, nr, 0)
	pfe3=flx*1.d20
	call ncdn_local_flux(flx, 1, eri2, nr, 1)
	pfi3=flx*1.d20
	ak3=dt_er*(pfe3-pfi3)/coe
	eri3=er(nr)+1.932d3*ak1/2.197d3-7.2d3*ak2/2.197d3+7.296d3*ak3/2.197d3
	call ncdn_local_flux(flx, 1, eri3, nr, 0)
	pfe4=flx*1.d20
	call ncdn_local_flux(flx, 1, eri3, nr, 1)
	pfi4=flx*1.d20
	ak4=dt_er*(pfe4-pfi4)/coe
	eri4=er(nr)+4.39d2*ak1/2.16d2-8.d0*ak2+3.68d3*ak3/5.13d2-8.45d2*ak4/4.104d3
	call ncdn_local_flux(flx, 1, eri4, nr, 0)
	pfe5=flx*1.d20
	call ncdn_local_flux(flx, 1, eri4, nr, 1)
	pfi5=flx*1.d20
	ak5=dt_er*(pfe5-pfi5)/coe
	eri5=er(nr)-8.d0*ak1/2.7d1+2.d0*ak2-3.544d3*ak3/2.565d3+1.859d3*ak4/4.104d3-1.1d1*ak5/4.d1
	call ncdn_local_flux(flx, 1, eri5, nr, 0)
	pfe6=flx*1.d20
	call ncdn_local_flux(flx, 1, eri5, nr, 1)
	pfi6=flx*1.d20
	ak6=dt_er*(pfe6-pfi6)/coe
	er(nr)=er(nr)+1.6d1*ak1/1.35d2+6.656d3*ak3/1.2825d4+2.8561d4*ak4/5.643d4-9.d0*ak5/5.d1+2.d0*ak6/5.5d1
    enddo
    case(1)
      er(nrmax)=er_boundary_value
    case default
      print*,'ER BOUNDARY CONDITION ERROR. You must select ER_CALC_BOUNDARY_CONDITION = 0 .or. 1.'
      print*,'stopped @ ',__FILE__,__LINE__
      stop
    end select


#ifdef er_dynamics_log_20100716
!      if(k==1 .or. k==kmax) then
      write(820,'(a1,a19,22a20)')'#','1:time','2:rho','3:r','4:er(t)','5:g','6:bc','7:bs','8:c^2/vA^2','9:bc_1st','10:bc_2nd','11:bc_3rd','12:denominator','13:Ge-Gi','14:bc_2nd_1','15:bc_2nd_2','16:bc_2nd_3','17:pfe','18:pfi','19:er(t+1)','20:per','21:pper','22:pfe','23:pfi'
      do nr=1,nrmax
      write(820,'(23e20.10)') t*DT+dt_er*k,RG(nr),RG(nr)*RA,er_stock(nr),g_stock(nr),bc_stock(nr),bs_stock(nr),c2va2(nr),bc_1st(nr),bc_2nd(nr),bc_3rd(nr),denominator(nr),ge_minus_gi(nr),bc_2nd_1(nr),bc_2nd_2(nr),bc_2nd_3(nr),pfe_stock(nr),pfi_stock(nr),er(nr),per_stock(nr),pper_stock(nr),pfe_stock(nr),pfi_stock(nr)
      enddo
      write(820,*)
 !     endif
      flg_log=1   
#endif
#ifdef er_dynamics_log_QQ
       icnt=1
       if(flg_log_QQ==0)then
       write(812,'(a5,5a20)') '# i  ','RG','time','Er','Qe','Qi'
       flg_log_QQ=1
       endif
       if(k==kmax)then 
       do i=1,nrmax
       call ncdn_local_flux(flxe, 2, er(i), i, 0)
       call ncdn_local_flux(flxi, 2, er(i), i, 1)
       write(812,'(i5,5e20.10)') i,rg(i),t*DT+dt_er*k,er(i),flxe*1d20,flxi*1d20
       enddo
       write(812,*)
       endif
#endif
       enddo


    case(33) ! Crank-Nicolson  
    kmax=int(1.d0/dt_er)/int(1.d0/DT)
#ifdef er_dynamics_log_20100716    
    NUM_OUT_LOG=kmax
#endif    
!    kmax=int(1.d0/dt_er)/int(1.d0/DT)*500

    do k=1,kmax
    if(mod(k,100)==0 .or. k==1)then
!      print'(i5,a1,i5)',k,'/',kmax
!      do nr=1,nrmax
!      write(830,'(2i6,2e20.10)') k,nr,rg(nr),er(nr)
!      print '(2i6,2e20.10)' ,k,nr,rg(nr),er(nr)
!      enddo
      flg_log_TODA=1
    endif
      
  
 !   do nr=1,nrmax
!    er(1)=0d0
!    er(nrmax)=12.0d0
!     dea(nr)=1.d-2/8e-14
!      dea(nr)=1.0d-3
!     dea(nr)=1.e-10
!     dea(nr)=0.d0
!    dea(nr)=(1.0d0/(RN(NR,1)*10.0d0))*(1.d0/10.d0)
!    enddo
!    do nr=2,nrmax-1
!     pdea(nr)=(dea(nr+1)-dea(nr-1))/(2.d0*DR)
!    enddo


    do nr=2,nrmax-1
    pdea(nr)=(dea(nr+1)-dea(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    pdea(1)=(dea(2)-dea(1))/(RA*(RG(2)-RG(1)))
    pdea(nrmax)=(dea(nrmax)-dea(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))
 
    do nr=2,nrmax-1
    per(nr)=(er(nr+1)-er(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    per(1)=(er(2)-er(1))/(RA*(RG(2)-RG(1)))
    per(nrmax)=(er(nrmax)-er(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))

    do nr=2,nrmax-1
    pper(nr)=(per(nr+1)-per(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    pper(1)=(per(2)-per(1))/(RA*(RG(2)-RG(1)))
    pper(nrmax)=(per(nrmax)-per(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))
   
    call ncdn_pre_proc
    do nr=1,nrmax
!    ep_perp(nr)=8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)*(1.d0+2.d0*QP(nr)**2)
     g(nr)= 8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
    enddo
    do nr=2,nrmax-1
     call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
     pfe(nr)=flxe*1.0d+20
     call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
     pfi(nr)=flxi*1.0d+20
     bs(nr)=(1.d0/6.d0/dt_er-dea(nr)/2.d0/DR**2)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     bc(nr)=(er_gamma_factor*(pfe(nr)-pfi(nr))/g(nr)+dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr)+(er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
!     bc(nr)=((pfe(nr)-pfi(nr))/g(nr)+dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr)+er(nr)/dt_er)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
!    bc(nr)=((pfe(nr)-pfi(nr))/g(nr)+(dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr))/AEE/g(nr)+(er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
!     bc(nr)=((pfe(nr)-pfi(nr))/g(nr)+dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr))/(2.d0/3.d0/dt_er+dea(nr)/DR**2)+(er(nr-1)/6.d0+4.d0*er(nr)/6.d0+er(nr+1)/6.d0)
!     bc(nr)=(AEE*(pfe(nr)-pfi(nr))/g(nr)+(dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr))+(er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
#ifdef er_dynamics_log_20100716
     per_stock(nr)=per(nr)
     pper_stock(nr)=pper(nr)
     bs_stock(nr)=bs(nr)
     bc_stock(nr)=bc(nr)
 !    bc_stock(nr)=((pfe(nr)-pfi(nr))/g(nr)+dea(nr)*pper(nr)+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr)+(er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)*(dt_er)
!     bc_stock(nr)=(-(pfe(nr)-pfi(nr))/g(nr)+dea(nr)*pper(nr)+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr)+(er(nr))/dt_er)*(dt_er)
     er_stock(nr)=er(nr)
     g_stock(nr)=g(nr)
     pfe_stock(nr)=pfe(nr)
     pfi_stock(nr)=pfi(nr)
     c2va2(nr)=vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2
     ge_minus_gi(nr)=(pfe(nr)-pfi(nr))
     bc_1st(nr)=(er_gamma_factor*(pfe(nr)-pfi(nr))/g(nr))/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     bc_2nd(nr)=(dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr))/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     bc_3rd(nr)=((er(nr+1)/6.d0+2.d0*er(nr)/3.d0+er(nr-1)/6.d0)/dt_er)/(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     denominator(nr)=(2.d0/3.d0/dt_er+dea(nr)/DR**2)
     bc_2nd_1(nr)=dea(nr)*pper(nr)/2.d0
     bc_2nd_2(nr)=dea(nr)*per(nr)/(rg(nr)*ra)
     bc_2nd_3(nr)=pdea(nr)*per(nr)
#endif
    enddo
    
#ifdef er_dynamics_log_20100716
     er_stock(1)=er(1)
     er_stock(nrmax)=er(nrmax)
#endif
    b1(1)=bc(2)-bs(2)*er(1)
  
    do nr=2,nrmax-2-1
    b1(nr)=bc(nr+1)
    enddo
  
    b1(nrmax-2)=bc(nrmax-2+1)-bs(nrmax-2+1)*er(nrmax)
!     b1(nrmax)=bc(nrmax-2+1)-bs(nrmax-2+1)*er(nrmax)
!
    do j=1,nrmax-2
    do i=1,3
    a(i,j)=0.d0
    enddo
    enddo
!
    a(2,1)=1.d0
    a(3,1)=bs(2)
!
    do i=2,nrmax-2-1
    a(1,i)=bs(i+1)
    a(2,i)=1.d0
    a(3,i)=bs(i+1)
    enddo
!
    a(2,nrmax-2)=1.d0
    a(1,nrmax-2)=bs(nrmax-2+1)

    call BANDRD(a,b1,nrmax-2,3,3,ierr)
!
    do nr=1,nrmax-2
     er2(nr)=b1(nr)
    enddo
    do nr=2,nrmax-1
     er(nr)=er2(nr-1)
    enddo
    er(1)=0.e0
!
    select case(ER_CALC_BOUNDARY_CONDITION)
    case(0)
    coe=g(nrmax)
!
	do nr=nrmax,nrmax
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 0)
	pfe1=flx*1.d20
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 1)
	pfi1=flx*1.d20
	ak1=dt_er*(pfe1-pfi1)/coe
	eri1=er(nr)+ak1/4.d0
	call ncdn_local_flux(flx, 1, eri1, nr, 0)
	pfe2=flx*1.d20
	call ncdn_local_flux(flx, 1, eri1, nr, 1)
	pfi2=flx*1.d20
	ak2=dt_er*(pfe2-pfi2)/coe
	eri2=er(nr)+3.d0*ak1/3.2d1+9.d0*ak2/3.2d1
	call ncdn_local_flux(flx, 1, eri2, nr, 0)
	pfe3=flx*1.d20
	call ncdn_local_flux(flx, 1, eri2, nr, 1)
	pfi3=flx*1.d20
	ak3=dt_er*(pfe3-pfi3)/coe
	eri3=er(nr)+1.932d3*ak1/2.197d3-7.2d3*ak2/2.197d3+7.296d3*ak3/2.197d3
	call ncdn_local_flux(flx, 1, eri3, nr, 0)
	pfe4=flx*1.d20
	call ncdn_local_flux(flx, 1, eri3, nr, 1)
	pfi4=flx*1.d20
	ak4=dt_er*(pfe4-pfi4)/coe
	eri4=er(nr)+4.39d2*ak1/2.16d2-8.d0*ak2+3.68d3*ak3/5.13d2-8.45d2*ak4/4.104d3
	call ncdn_local_flux(flx, 1, eri4, nr, 0)
	pfe5=flx*1.d20
	call ncdn_local_flux(flx, 1, eri4, nr, 1)
	pfi5=flx*1.d20
	ak5=dt_er*(pfe5-pfi5)/coe
	eri5=er(nr)-8.d0*ak1/2.7d1+2.d0*ak2-3.544d3*ak3/2.565d3+1.859d3*ak4/4.104d3-1.1d1*ak5/4.d1
	call ncdn_local_flux(flx, 1, eri5, nr, 0)
	pfe6=flx*1.d20
	call ncdn_local_flux(flx, 1, eri5, nr, 1)
	pfi6=flx*1.d20
	ak6=dt_er*(pfe6-pfi6)/coe
	er(nr)=er(nr)+1.6d1*ak1/1.35d2+6.656d3*ak3/1.2825d4+2.8561d4*ak4/5.643d4-9.d0*ak5/5.d1+2.d0*ak6/5.5d1
    enddo
    case(1)
      er(nrmax)=er_boundary_value
    case default
      print*,'ER BOUNDARY CONDITION ERROR. You must select ER_CALC_BOUNDARY_CONDITION = 0 .or. 1.'
      print*,'stopped @ ',__FILE__,__LINE__
      stop
    end select

#ifdef er_dynamics_log_20100716
!      if(k==1 .or. k==kmax) then
      write(820,'(a1,a19,22a20)')'#','1:time','2:rho','3:r','4:er(t)','5:g','6:bc','7:bs','8:c^2/vA^2','9:bc_1st','10:bc_2nd','11:bc_3rd','12:denominator','13:Ge-Gi','14:bc_2nd_1','15:bc_2nd_2','16:bc_2nd_3','17:pfe','18:pfi','19:er(t+1)','20:per','21:pper','22:pfe','23:pfi'
      do nr=1,nrmax
      write(820,'(23e20.10)') t*DT+dt_er*k,RG(nr),RG(nr)*RA,er_stock(nr),g_stock(nr),bc_stock(nr),bs_stock(nr),c2va2(nr),bc_1st(nr),bc_2nd(nr),bc_3rd(nr),denominator(nr),ge_minus_gi(nr),bc_2nd_1(nr),bc_2nd_2(nr),bc_2nd_3(nr),pfe_stock(nr),pfi_stock(nr),er(nr),per_stock(nr),pper_stock(nr),pfe_stock(nr),pfi_stock(nr)
      enddo
      write(820,*)
 !     endif
      flg_log=1   
#endif
#ifdef er_dynamics_log_QQ
       icnt=1
       if(flg_log_QQ==0)then
       write(812,'(a5,5a20)') '# i  ','RG','time','Er','Qe','Qi'
       flg_log_QQ=1
       endif
       if(k==kmax)then 
       do i=1,nrmax
       call ncdn_local_flux(flxe, 2, er(i), i, 0)
       call ncdn_local_flux(flxi, 2, er(i), i, 1)
       write(812,'(i5,5e20.10)') i,rg(i),t*DT+dt_er*k,er(i),flxe*1d20,flxi*1d20
       enddo
       write(812,*)
       endif
#endif
       enddo
       
    case(34) 
        ! balance between particle flux and diffusion of Er.
        ! CAUTION ! 
        ! Only assumed, Dea=const.
     open(781,file='Grad_LOG.txt')

      flg_log_TODA=1
      k=1
      kmax=1

!    do nr=1,nrmax
!      dea(nr)=0.d0
!    enddo 

     print*,'NOW CALCULATING ER'
     delta_er=1.d-2
     cnt_er_iteration=0
     epsilon_er=0.1d0
     error_balance=10000.0d0
     er_change_factor=5.d-17
     MAX_er_iteration=150
     pre_sum_error_balance=1.d40

    call ncdn_pre_proc

	do while(1)

    do nr=2,nrmax-1
    pdea(nr)=(dea(nr+1)-dea(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    pdea(1)=(dea(2)-dea(1))/(RA*(RG(2)-RG(1)))
    pdea(nrmax)=(dea(nrmax)-dea(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))
 
    do nr=2,nrmax-1
    per(nr)=(er(nr+1)-er(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    per(1)=(er(2)-er(1))/(RA*(RG(2)-RG(1)))
    per(nrmax)=(er(nrmax)-er(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))

    do nr=2,nrmax-1
    pper(nr)=(per(nr+1)-per(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    pper(1)=(per(2)-per(1))/(RA*(RG(2)-RG(1)))
    pper(nrmax)=(per(nrmax)-per(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))

	do nr=1,nrmax
      g(nr)= 8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
    enddo
	  
	sum_error_balance=0.0d0  
	do nr=1,nrmax
	  call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
	  pfe(nr)=flxe*1.0d+20
	  call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
	  pfi(nr)=flxi*1.0d+20
	  coe_aa(nr)= (pfe(nr)-pfi(nr))/g(nr)
      coe_bb(nr)= dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr)
	  error_balance(nr)=0.5d0*(coe_aa(nr)+coe_bb(nr))*(coe_aa(nr)+coe_bb(nr))
	  sum_error_balance=sum_error_balance+error_balance(nr)
!      print'(i4,3e20.10)',nr,coe_aa(nr),coe_bb(nr),error_balance(nr)
    enddo
	
!	print*,'sum_error_balance',sum_error_balance
	if(sum_error_balance < epsilon_er) exit

!   Calculation of gradiate of Particl flux and DIFF.
    sum_error_balance=0.0
    do nr=1,nrmax
	  er1(nr)=er(nr)+delta_er
	  call ncdn_local_flux(flxe, 1, er1(nr), nr, 0)
	  del_pfe(nr)=flxe*1.0d+20
	  call ncdn_local_flux(flxi, 1, er1(nr), nr, 1)
	  del_pfi(nr)=flxi*1.0d+20
	  del_coe_aa(nr)= ((del_pfe(nr)-pfe(nr))/delta_er-(del_pfi(nr)-pfi(nr)/delta_er))/g(nr)
      del_coe_bb(nr)=0.0d0
	  del_error_balance(nr)=(coe_aa(nr)+coe_bb(nr))*(del_coe_aa(nr))
!      print'(i4,7e20.10)', nr,pfe(nr),pfi(nr),pfe(nr)-pfi(nr),del_pfe(nr),del_pfi(nr),del_pfe(nr)-del_pfi(nr),del_error_balance(nr)
      error_balance(nr)=0.5d0*(del_coe_aa(nr)+del_coe_bb(nr))*(del_coe_aa(nr)+del_coe_bb(nr))
	  sum_error_balance=sum_error_balance+error_balance(nr)
	enddo
!	print*,'sum_error_balance',sum_error_balance
         
	sum_error_balance=0.0d0         
	do nr=1,nrmax
	  er(nr)=er(nr)+del_error_balance(nr)*er_change_factor
	enddo           
	er(1)=0.0d0

    select case(ER_CALC_BOUNDARY_CONDITION)
    case(0)
    coe=g(nrmax)
!
	do nr=nrmax,nrmax
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 0)
	pfe1=flx*1.d20
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 1)
	pfi1=flx*1.d20
	ak1=dt_er*(pfe1-pfi1)/coe
	eri1=er(nr)+ak1/4.d0
	call ncdn_local_flux(flx, 1, eri1, nr, 0)
	pfe2=flx*1.d20
	call ncdn_local_flux(flx, 1, eri1, nr, 1)
	pfi2=flx*1.d20
	ak2=dt_er*(pfe2-pfi2)/coe
	eri2=er(nr)+3.d0*ak1/3.2d1+9.d0*ak2/3.2d1
	call ncdn_local_flux(flx, 1, eri2, nr, 0)
	pfe3=flx*1.d20
	call ncdn_local_flux(flx, 1, eri2, nr, 1)
	pfi3=flx*1.d20
	ak3=dt_er*(pfe3-pfi3)/coe
	eri3=er(nr)+1.932d3*ak1/2.197d3-7.2d3*ak2/2.197d3+7.296d3*ak3/2.197d3
	call ncdn_local_flux(flx, 1, eri3, nr, 0)
	pfe4=flx*1.d20
	call ncdn_local_flux(flx, 1, eri3, nr, 1)
	pfi4=flx*1.d20
	ak4=dt_er*(pfe4-pfi4)/coe
	eri4=er(nr)+4.39d2*ak1/2.16d2-8.d0*ak2+3.68d3*ak3/5.13d2-8.45d2*ak4/4.104d3
	call ncdn_local_flux(flx, 1, eri4, nr, 0)
	pfe5=flx*1.d20
	call ncdn_local_flux(flx, 1, eri4, nr, 1)
	pfi5=flx*1.d20
	ak5=dt_er*(pfe5-pfi5)/coe
	eri5=er(nr)-8.d0*ak1/2.7d1+2.d0*ak2-3.544d3*ak3/2.565d3+1.859d3*ak4/4.104d3-1.1d1*ak5/4.d1
	call ncdn_local_flux(flx, 1, eri5, nr, 0)
	pfe6=flx*1.d20
	call ncdn_local_flux(flx, 1, eri5, nr, 1)
	pfi6=flx*1.d20
	ak6=dt_er*(pfe6-pfi6)/coe
	er(nr)=er(nr)+1.6d1*ak1/1.35d2+6.656d3*ak3/1.2825d4+2.8561d4*ak4/5.643d4-9.d0*ak5/5.d1+2.d0*ak6/5.5d1
    enddo
    case(1)
      er(nrmax)=er_boundary_value
    case default
      print*,'ER BOUNDARY CONDITION ERROR. You must select ER_CALC_BOUNDARY_CONDITION = 0 .or. 1.'
      print*,'stopped @ ',__FILE__,__LINE__
      stop
    end select 

    do nr=2,nrmax-1
    per(nr)=(er(nr+1)-er(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    per(1)=(er(2)-er(1))/(RA*(RG(2)-RG(1)))
    per(nrmax)=(er(nrmax)-er(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))

    do nr=2,nrmax-1
    pper(nr)=(per(nr+1)-per(nr-1))/(RA*(RG(nr+1)-RG(nr-1)))
    enddo
    pper(1)=(per(2)-per(1))/(RA*(RG(2)-RG(1)))
    pper(nrmax)=(per(nrmax)-per(nrmax-1))/(RA*(RG(nrmax)-RG(nrmax-1)))
	
    do nr=1,nrmax
	  call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
	  pfe(nr)=flxe*1.0d+20
	  call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
	  pfi(nr)=flxi*1.0d+20
	  coe_aa(nr)= (pfe(nr)-pfi(nr))/g(nr)
      coe_bb(nr)= dea(nr)*pper(nr)/2.d0+dea(nr)*per(nr)/(rg(nr)*ra)+pdea(nr)*per(nr)
	  error_balance(nr)=0.5d0*(coe_aa(nr)+coe_bb(nr))*(coe_aa(nr)+coe_bb(nr))
	  sum_error_balance=sum_error_balance+error_balance(nr)
!      print'(i4,3e20.10)',nr,coe_aa(nr),coe_bb(nr),error_balance(nr)
    enddo
!	print*,'sum_error_balance',sum_error_balance
!    write(781,'(i5,e20.10)')cnt_er_iteration,sum_error_balance
	cnt_er_iteration=cnt_er_iteration+1
    print'(a10,i4,3e20.10)','NOW LOOP', cnt_er_iteration,sum_error_balance,pre_sum_error_balance,abs((sum_error_balance-pre_sum_error_balance)/sum_error_balance)
	if(cnt_er_iteration > MAX_er_iteration .or. (abs((sum_error_balance-pre_sum_error_balance)/sum_error_balance)<0.005)) then
    print*,'ERROR eval',(sum_error_balance-pre_sum_error_balance)/sum_error_balance
!	  print*,'Can not balance between GAMMA and DIFF. '
!      close(781)
!	  stop
	  exit
	else
      pre_sum_error_balance=sum_error_balance
    endif
	
	enddo ! do while

      write(781,'(a1,a19,12a20)')&
      & '#','1:time','2:rho','3:r','4:er(t)','5:g','6:pfe','7:pfi','8:GAMMAterm','9:dea','10:per','11:pper','12:DIFFterm','13:error(nr)'
      do nr=1,nrmax
      write(781,'(13e20.10)') &
      & t*DT,RG(nr),RG(nr)*RA,er(nr),g(nr),pfe(nr),pfi(nr),coe_aa(nr),dea(nr),per(nr),pper(nr),coe_bb(nr),error_balance(nr)
      enddo
    write(781,*)


!      print'(2a4,3a20)','k','nr','rg(nr)','er0(nr)','er(nr)'
!      do nr=1,nrmax
!         print'(2i4,3e20.10)',k,nr,rg(nr),er0(nr),er(nr)
!      enddo
#ifdef er_dynamics_log_QQ
       icnt=1
       if(flg_log_QQ==0)then
       write(812,'(a5,5a20)') '# i  ','RG','time','Er','Qe','Qi'
       flg_log_QQ=1
       endif
       if(k==kmax)then 
       do i=1,nrmax
       call ncdn_local_flux(flxe, 2, er(i), i, 0)
       call ncdn_local_flux(flxi, 2, er(i), i, 1)
       write(812,'(i5,5e20.10)') i,rg(i),t*DT+dt_er*k,er(i),flxe*1d20,flxi*1d20
       enddo
       write(812,*)
       endif
#endif

    case(35) 
!   TEST routine 
!
!   TDMA: Tro-Diagonal Matrix Algorithm
!   | coe_bb  -coe_cc                           | |er|  |coe_dd|
!   |-coe_aa   coe_bb  -coe_cc                  | |  |  |coe_dd|
!   |         -coe_aa   coe_bb  -coe_cc         | |  |= |coe_dd|
!   |                   ......                  | |  |  |coe_dd|
!   |                           -coe_aa   coe_bb| |  |  |coe_dd|
!
!   PP[NRMAX], QQ[NEMAX] :: 係数格納用
  
!  SET PARAMETER
 
    open(755,file='er_dynamics_Erlog.txt',status='replace')    
    open(756,file='er_dynamics_Erlog_last.txt',status='replace')
    open(757,file='er_dynamics_Erlog_1st.txt',status='replace')

!    dt_er=1e-6
    kmax=1000
!  dea_const=1.d-2

    k=1  
    do nr=1,nrmax
        er0(NR)=er(nr)
        write(757,'(2i5,3e20.10)')k,nr,k*dt_er,rg(nr),Er(nr)
    enddo
  
    do k=1,kmax
        call ncdn_pre_proc
        do nr=1,nrmax
!            dea(nr)=dea_const
            g(nr)= 8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
            coe_aa(nr)=(1.d0/dt_er)-2.d0*dea(nr)/DR**2
            coe_bb(nr)=(1.d0/DR**2 + 1.d0/(2.d0*rg(nr)*ra*DR))*dea(nr)
            coe_cc(nr)=(1.d0/DR**2 - 1.d0/(2.d0*rg(nr)*ra*DR))*dea(nr)

            call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
            pfe(nr)=flxe*1.0d+20
            call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
            pfi(nr)=flxi*1.0d+20
            coe_dd(nr)=er(nr)/dt_er + (pfe(nr)-pfi(nr))/g(nr)*er_gamma_factor
        enddo

        PP(1)=coe_bb(1)/coe_aa(1)
        QQ(1)=coe_dd(1)/coe_aa(1)
        
        do nr=2,nrmax
            PP(nr)=coe_bb(nr)/(coe_aa(nr)-coe_cc(nr)*PP(nr-1))
            QQ(nr)=(coe_dd(nr)+coe_cc(nr)*QQ(nr-1))/(coe_aa(nr)-coe_cc(nr)*PP(nr-1))
        enddo
        
        er(nrmax)=QQ(nrmax)
        
        do nr=nrmax-1,1,-1
            er(nr)=PP(nr)*er(nr+1)+QQ(nr)  
        enddo
    
        er(1)=0.e0
    
    select case(ER_CALC_BOUNDARY_CONDITION)
    case(0)
        coe=g(nrmax)
        do nr=nrmax,nrmax
            call ncdn_local_flux(flx, 1, er(nrmax), nr, 0)
            pfe1=flx*1.d20
            call ncdn_local_flux(flx, 1, er(nrmax), nr, 1)
            pfi1=flx*1.d20
            ak1=dt_er*(pfe1-pfi1)/coe
            eri1=er(nr)+ak1/4.d0
            call ncdn_local_flux(flx, 1, eri1, nr, 0)
            pfe2=flx*1.d20
            call ncdn_local_flux(flx, 1, eri1, nr, 1)
            pfi2=flx*1.d20
            ak2=dt_er*(pfe2-pfi2)/coe
            eri2=er(nr)+3.d0*ak1/3.2d1+9.d0*ak2/3.2d1
            call ncdn_local_flux(flx, 1, eri2, nr, 0)
            pfe3=flx*1.d20
            call ncdn_local_flux(flx, 1, eri2, nr, 1)
            pfi3=flx*1.d20
            ak3=dt_er*(pfe3-pfi3)/coe
            eri3=er(nr)+1.932d3*ak1/2.197d3-7.2d3*ak2/2.197d3+7.296d3*ak3/2.197d3
            call ncdn_local_flux(flx, 1, eri3, nr, 0)
            pfe4=flx*1.d20
            call ncdn_local_flux(flx, 1, eri3, nr, 1)
            pfi4=flx*1.d20
            ak4=dt_er*(pfe4-pfi4)/coe
            eri4=er(nr)+4.39d2*ak1/2.16d2-8.d0*ak2+3.68d3*ak3/5.13d2-8.45d2*ak4/4.104d3
            call ncdn_local_flux(flx, 1, eri4, nr, 0)
            pfe5=flx*1.d20
            call ncdn_local_flux(flx, 1, eri4, nr, 1)
            pfi5=flx*1.d20
            ak5=dt_er*(pfe5-pfi5)/coe
            eri5=er(nr)-8.d0*ak1/2.7d1+2.d0*ak2-3.544d3*ak3/2.565d3+1.859d3*ak4/4.104d3-1.1d1*ak5/4.d1
            call ncdn_local_flux(flx, 1, eri5, nr, 0)
            pfe6=flx*1.d20
            call ncdn_local_flux(flx, 1, eri5, nr, 1)
            pfi6=flx*1.d20
            ak6=dt_er*(pfe6-pfi6)/coe
            er(nr)=er(nr)+1.6d1*ak1/1.35d2+6.656d3*ak3/1.2825d4+2.8561d4*ak4/5.643d4-9.d0*ak5/5.d1+2.d0*ak6/5.5d1
        enddo
    case(1)
      er(nrmax)=er_boundary_value
    case default
      print*,'ER BOUNDARY CONDITION ERROR. You must select ER_CALC_BOUNDARY_CONDITION = 0 .or. 1.'
      print*,'stopped @ ',__FILE__,__LINE__
      stop
    end select
  
        if(mod(k,kmax/50)==0)then
            print*,'k=',k
            do nr=1,nrmax
                write(755,'(2i5,3e20.10)')k,nr,k*dt_er,rg(nr),Er(nr)
            enddo
            write(755,*)
        endif
    
    enddo
    do nr=1,nrmax
        write(756,'(2i5,3e20.10)')k,nr,k*dt_er,rg(nr),Er(nr)
    enddo
    close(755)
    close(756)
    close(757)
    stop
    
    case(36) ! Euler implicit NEW @20110427

!   TDMA: Tro-Diagonal Matrix Algorithm
!   | coe_bb  -coe_cc                           | |er|  |coe_dd|
!   |-coe_aa   coe_bb  -coe_cc                  | |  |  |coe_dd|
!   |         -coe_aa   coe_bb  -coe_cc         | |  |= |coe_dd|
!   |                   ......                  | |  |  |coe_dd|
!   |                           -coe_aa   coe_bb| |  |  |coe_dd|
!
!   PP[NRMAX], QQ[NEMAX] :: 係数格納用
    kmax=int(1.d0/dt_er)/int(1.d0/DT)
#ifdef er_dynamics_log_20100716    
    NUM_OUT_LOG=kmax
#endif    
!    kmax=int(1.d0/dt_er)/int(1.d0/DT)*500

    do k=1,kmax
    print*,'k/kmax',k,'/',kmax,':DT',DT
    if(mod(k,100)==0 .or. k==1)then
!      print'(i5,a1,i5)',k,'/',kmax
!      do nr=1,nrmax
!      write(830,'(2i6,2e20.10)') k,nr,rg(nr),er(nr)
!      print '(2i6,2e20.10)' ,k,nr,rg(nr),er(nr)
!      enddo
      flg_log_TODA=1
    endif

    call ncdn_pre_proc
    do nr=1,nrmax
!      g(nr)=(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
      g(nr)= 8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2)/AEE*(1.d0+2.d0*QP(nr)**2)
      dt_dr2=dt_er/(DR*DR)
      coe_cc(nr)=-dt_dr2*dea(nr)
      coe_aa(nr)=1.d0+2.d0*dt_dr2*dea(nr)
      coe_bb(nr)=-dt_dr2*dea(nr)

      call ncdn_local_flux(flxe, 1, er(nr), nr, 0)
      pfe(nr)=flxe*1.0d+20
      call ncdn_local_flux(flxi, 1, er(nr), nr, 1)
      pfi(nr)=flxi*1.0d+20
      coe_dd(nr)=er(nr) 
      coe_dd(nr)=er(nr) + er_gamma_factor*(pfe(nr)-pfi(nr))/g(nr)*1.602d-19*dt_er

    enddo
    
!    do nr=1,nrmax
!    print'(i4,8e14.4)',nr,coe_aa(nr),coe_bb(nr),coe_cc(nr),coe_dd(nr), er(nr),er_gamma_factor*(pfe(nr)-pfi(nr))/g(nr)*1.602d-19*dt_er,(pfe(nr)-pfi(nr)),g(nr)
!    enddo

    PP(1)=coe_bb(1)/coe_aa(1)
    QQ(1)=coe_dd(1)/coe_aa(1)
    do nr=2,nrmax
    PP(nr)=coe_bb(nr)/(coe_aa(nr)-coe_cc(nr)*PP(nr-1))
    QQ(nr)=(coe_dd(nr)+coe_cc(nr)*QQ(nr-1))/(coe_aa(nr)-coe_cc(nr)*PP(nr-1))
    enddo
    er(nrmax)=QQ(nrmax)
    do nr=nrmax-1,1,-1
    er(nr)=PP(nr)*er(nr+1)+QQ(nr)  
    enddo

  
    er(1)=0.e0
    
    select case(ER_CALC_BOUNDARY_CONDITION)
    case(0)
    coe=g(nrmax)
!
	do nr=nrmax,nrmax
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 0)
	pfe1=flx*1.d20
	call ncdn_local_flux(flx, 1, er(nrmax), nr, 1)
	pfi1=flx*1.d20
	ak1=dt_er*(pfe1-pfi1)/coe
	eri1=er(nr)+ak1/4.d0
	call ncdn_local_flux(flx, 1, eri1, nr, 0)
	pfe2=flx*1.d20
	call ncdn_local_flux(flx, 1, eri1, nr, 1)
	pfi2=flx*1.d20
	ak2=dt_er*(pfe2-pfi2)/coe
	eri2=er(nr)+3.d0*ak1/3.2d1+9.d0*ak2/3.2d1
	call ncdn_local_flux(flx, 1, eri2, nr, 0)
	pfe3=flx*1.d20
	call ncdn_local_flux(flx, 1, eri2, nr, 1)
	pfi3=flx*1.d20
	ak3=dt_er*(pfe3-pfi3)/coe
	eri3=er(nr)+1.932d3*ak1/2.197d3-7.2d3*ak2/2.197d3+7.296d3*ak3/2.197d3
	call ncdn_local_flux(flx, 1, eri3, nr, 0)
	pfe4=flx*1.d20
	call ncdn_local_flux(flx, 1, eri3, nr, 1)
	pfi4=flx*1.d20
	ak4=dt_er*(pfe4-pfi4)/coe
	eri4=er(nr)+4.39d2*ak1/2.16d2-8.d0*ak2+3.68d3*ak3/5.13d2-8.45d2*ak4/4.104d3
	call ncdn_local_flux(flx, 1, eri4, nr, 0)
	pfe5=flx*1.d20
	call ncdn_local_flux(flx, 1, eri4, nr, 1)
	pfi5=flx*1.d20
	ak5=dt_er*(pfe5-pfi5)/coe
	eri5=er(nr)-8.d0*ak1/2.7d1+2.d0*ak2-3.544d3*ak3/2.565d3+1.859d3*ak4/4.104d3-1.1d1*ak5/4.d1
	call ncdn_local_flux(flx, 1, eri5, nr, 0)
	pfe6=flx*1.d20
	call ncdn_local_flux(flx, 1, eri5, nr, 1)
	pfi6=flx*1.d20
	ak6=dt_er*(pfe6-pfi6)/coe
	er(nr)=er(nr)+1.6d1*ak1/1.35d2+6.656d3*ak3/1.2825d4+2.8561d4*ak4/5.643d4-9.d0*ak5/5.d1+2.d0*ak6/5.5d1
    enddo
    case(1)
      er(nrmax)=er_boundary_value
    case default
      print*,'ER BOUNDARY CONDITION ERROR. You must select ER_CALC_BOUNDARY_CONDITION = 0 .or. 1.'
      print*,'stopped @ ',__FILE__,__LINE__
      stop
    end select

#ifdef er_dynamics_log_20100716
!      if(k==1 .or. k==kmax) then
      write(820,'(a1,a19,22a20)')'#','1:time','2:rho','3:r','4:er(t)','5:g','6:bc','7:bs','8:c^2/vA^2','9:bc_1st','10:bc_2nd','11:bc_3rd','12:denominator','13:Ge-Gi','14:bc_2nd_1','15:bc_2nd_2','16:bc_2nd_3','17:pfe','18:pfi','19:er(t+1)','20:per','21:pper','22:pfe','23:pfi'
      do nr=1,nrmax
      write(820,'(23e20.10)') t*DT+dt_er*k,RG(nr),RG(nr)*RA,er_stock(nr),g_stock(nr),bc_stock(nr),bs_stock(nr),c2va2(nr),bc_1st(nr),bc_2nd(nr),bc_3rd(nr),denominator(nr),ge_minus_gi(nr),bc_2nd_1(nr),bc_2nd_2(nr),bc_2nd_3(nr),pfe_stock(nr),pfi_stock(nr),er(nr),per_stock(nr),pper_stock(nr),pfe_stock(nr),pfi_stock(nr)
      enddo
      write(820,*)
 !     endif
      flg_log=1   
#endif
#ifdef er_dynamics_log_QQ
       icnt=1
       if(flg_log_QQ==0)then
       write(812,'(a5,5a20)') '# i  ','RG','time','Er','Qe','Qi'
       flg_log_QQ=1
       endif
!       if(k==kmax)then 
       if((mod(k,10))==0)then 
       do i=1,nrmax
       call ncdn_local_flux(flxe, 2, er(i), i, 0)
       call ncdn_local_flux(flxi, 2, er(i), i, 1)
       write(812,'(i5,5e20.10)') i,rg(i),t*DT+dt_er*k,er(i),flxe*1d20,flxi*1d20
       enddo
       write(812,*)
       endif
#endif
       enddo
        case(37)
            call cpu_time(cpu_time_start)
            open(839,file='er_dynamics_log_iteration_ALL.txt',status='replace')
            close(839)       
            open(842,file='er_dynamics_log_iteration.txt',status='replace')
            close(842) 
            open(841,file='er_dynamics_iteration_errorlog.txt',status='replace')
            close(841)
            open(448,file='er_dynamics_log_iteration_TEST.txt',status='replace')

            it_cnt=0
            delta_er=100.d0
            threshold_relative_error_SUM=1.0d-255
            threshold_NEW_error_SUM=1.0d-15
            relative_error_SUM=1.d100
            SUM_Error_Func=0.0d0
            it_cnt_max=10000
            EF_normalized_factor=(8.854d-12/AEE/(rn(1,1)*1.d20))
           
            call ncdn_pre_proc
            
            er(1)=0.0d0
            er(nrmax)=er_boundary_value
            er_0(1)=er(1)
            er_0(nrmax)=er(nrmax)            
            er_1(1)=er(1)
            er_1(nrmax)=er(nrmax)            
            er_2(1)=er(1)
            er_2(nrmax)=er(nrmax)            
            er_3(1)=er(1)
            er_3(nrmax)=er(nrmax)            

            do nr=2,nrmax-1
                er_0(nr)=er(nr)
                er_1(nr)=er(nr)
            enddo                                            

            if (DEA_SETTING_CONDITION == 1) then
                do nr=1,nrmax
                    dea(nr)=1d-10*exp(17.956*rg(nr))
                enddo
            endif



!!          ===== CAUTION !!   ======
!!          = 現在は dea(nr) のEr依存性は考慮していない．
!!          = そのため，dea_der(nr), dea_drder(nr) にはゼロを代入しています．
!!          = 将来的には，Er 依存性を考慮して計算する必要があります．

!            do nr=1,nrmax-1
!                dea_drder(nr)=0.d0
!            enddo
!            do nr=1,nrmax-1
!                dea_der(nr)=0.d0
!            enddo
            do nr=1,nrmax
                g(nr)= 1.d0/((8.854d-12*(1.d0+vc**2*4.d-7*PI*1.67d-27*rn(nr,1)*1.0d+20/BB**2))*(1.d0+2.d0*QP(nr)**2))
            enddo
            
            do nr=2,nrmax-1
                term_alpha(nr)=(dea(nr-1)+4.d0*dea(nr)-dea(nr+1))/(4.d0*dr*dr)
                term_beta(nr)=-2.d0*dea(nr)/(dr*dr)-(dea(nr+1)-dea(nr-1))/(rg(nr)*ra*2.d0*dr)
                term_gamma(nr)=(-dea(nr-1)+4.d0*dea(nr)+dea(nr+1))/(4.d0*dr*dr)
                func_diff(nr)=(term_alpha(nr)*er_1(nr-1)+term_beta(nr)*er_1(nr)+term_gamma(nr)*er_1(nr+1))
            enddo
            do nr=2,nrmax-1
                call ncdn_local_flux(flxe, 1, er_1(nr), nr, 0)
                pfe_1(nr)=flxe*1.0d+20
                call ncdn_local_flux(flxi, 1, er_1(nr), nr, 1)
                pfi_1(nr)=flxi*1.0d+20
                func_gamma(nr)=(pfe_1(nr)-pfi_1(nr))*g(nr)*AEE
            enddo

            SUM_Error_Func_1=0.0d0

            do nr=2,nrmax-1
                Error_Func_1(nr)= (EF_normalized_factor**2) *0.5d0*(func_gamma(nr)+func_diff(nr))**2
                SUM_Error_Func_1=SUM_Error_Func_1+Error_Func_1(nr)
            enddo

            SUM_Error_Func=SUM_Error_Func_1
      
!       Calculate Error_Func_1st
        do 
            it_cnt=it_cnt+1
            
            do nr=2,nrmax-1
                call ncdn_local_flux(gamma_flx_e_0, 1, er_1(nr)-delta_er, nr, 0)
                call ncdn_local_flux(gamma_flx_e_1, 1, er_1(nr)+delta_er, nr, 0)
                call ncdn_local_flux(gamma_flx_i_0, 1, er_1(nr)-delta_er, nr, 1)
                call ncdn_local_flux(gamma_flx_i_1, 1, er_1(nr)+delta_er, nr, 1)
                der_flx_e=(gamma_flx_e_1-gamma_flx_e_0)/(2.d0*delta_er)
                der_flx_i=(gamma_flx_i_1-gamma_flx_i_0)/(2.d0*delta_er)
                dfunc_gamma_dEr(nr)=(der_flx_e-der_flx_i)*AEE*1.0d20*g(nr)
            enddo
            
            do nr=1,nrmax
                er_for_tmp(nr)=er_1(nr)
            enddo
            er_for_tmp(0)=0.0d0
            er_for_tmp(nrmax+1)=0.0d0
            
            do nr=2,nrmax-1
            write(448,'(2i6,13e20.10)'),it_cnt, nr, rg(nr),er_1(nr),er_for_tmp(nr),pfe_1(nr),pfi_1(nr),g(nr)*AEE*EF_normalized_factor,&
            &  term_alpha(nr), term_beta(nr),term_gamma(nr),  &
            &  er_for_tmp(nr-1), &
            &  func_gamma(nr), dfunc_gamma_dEr(nr),(pfe_1(nr)-pfi_1(nr))
            enddo            

            do nr=2,nrmax-1
                grad_SUM_Error_Func(nr)= -1.d0 *( &
            &   (term_alpha(nr+1)*term_alpha(nr+1) + term_beta(nr)*term_beta(nr) + term_gamma(nr-1)*term_gamma(nr-1))*er_for_tmp(nr) &
            & + (term_alpha(nr)*term_beta(nr) + term_beta(nr-1)*term_gamma(nr-1)) *er_for_tmp(nr-1) &
            & + (term_alpha(nr+1)*term_beta(nr+1) + term_beta(nr)*term_gamma(nr)) *er_for_tmp(nr+1) &
            & + (term_alpha(nr-1)*term_gamma(nr-1))*er_for_tmp(nr-2) &
            & + (term_alpha(nr+1)*term_gamma(nr+1))*er_for_tmp(nr+2) &
            & + (term_alpha(nr+1)*func_gamma(nr+1) + term_beta(nr)*func_gamma(nr) + term_gamma(nr-1)*func_gamma(nr-1)) &
            & + (term_alpha(nr)*er_for_tmp(nr-1) + term_beta(nr)*er_for_tmp(nr) + term_gamma(nr)*er_for_tmp(nr+1) +func_gamma(nr))*dfunc_gamma_dEr(nr) &
            & ) * (EF_normalized_factor**2)
            enddo
            
            do nr=2,nrmax-1
            write(448,'(2i6,11e20.10)'),it_cnt, nr, er_1(nr),er_for_tmp(nr),pfe_1(nr),pfi_1(nr),g(nr)*AEE*EF_normalized_factor,&
            &  term_alpha(nr), term_beta(nr),term_gamma(nr),  &
            &  er_for_tmp(nr-1), &
            &  func_gamma(nr), dfunc_gamma_dEr(nr)
            enddo
            write(448,*)

!       直線探索
!            SUM_Error_Func_1=0.0d0
!            SUM_Error_Func_2=0.0d0
!            SUM_Error_Func_3=0.0d0
            it_cnt_linesearch=0

            accelerator_factor_1=0.0d0

            MAX_grad_SUM_Error_Func=-1.d0            
            do nr=2,nrmax-1
                if(MAX_grad_SUM_Error_Func<dabs(grad_SUM_Error_Func(nr))) MAX_grad_SUM_Error_Func=dabs(grad_SUM_Error_Func(nr))
            enddo            
            accelerator_factor_2=10.d0/MAX_grad_SUM_Error_Func

            do nr=2,nrmax-1
                er_2(nr)=er_1(nr)+grad_SUM_Error_Func(nr)*accelerator_factor_2
            enddo

            do nr=2,nrmax-1 
                call ncdn_local_flux(flxe, 1, er_2(nr), nr, 0)
                pfe_2(nr)=flxe*1.0d+20
                call ncdn_local_flux(flxi, 1, er_2(nr), nr, 1)
                pfi_2(nr)=flxi*1.0d+20
                func_gamma(nr)=(pfe_2(nr)-pfi_2(nr))*g(nr)*AEE
            enddo 

            do nr=2,nrmax-1
                term_alpha(nr)=(dea(nr-1)+4.d0*dea(nr)-dea(nr+1))/(4.d0*dr*dr)
                term_beta(nr)=-2.d0*dea(nr)/(dr*dr)-(dea(nr+1)-dea(nr-1))/(rg(nr)*ra*2.d0*dr)
                term_gamma(nr)=(-dea(nr-1)+4.d0*dea(nr)+dea(nr+1))/(4.d0*dr*dr)
                func_diff(nr)=(term_alpha(nr)*er_2(nr-1)+term_beta(nr)*er_2(nr)+term_gamma(nr)*er_2(nr+1))
            enddo
            
            SUM_Error_Func_2=0.0d0
            do nr=2,nrmax-1
                Error_Func_2(nr)= (EF_normalized_factor**2) * 0.5d0*(func_gamma(nr)+func_diff(nr))**2
                SUM_Error_Func_2=SUM_Error_Func_2+Error_Func_2(nr)
            enddo            

            it_cnt_linesearch=0          
            if(SUM_Error_Func_2>SUM_Error_Func_1)then
                do
                    it_cnt_linesearch=it_cnt_linesearch+1
                    accelerator_factor_3=accelerator_factor_2*0.5d0

                    do nr=2,nrmax-1
                        er_3(nr)=er_1(nr)+grad_SUM_Error_Func(nr)*accelerator_factor_3
                    enddo
                    do nr=2,nrmax-1
                        call ncdn_local_flux(flxe, 1, er_3(nr), nr, 0)
                        pfe_3(nr)=flxe*1.0d+20
                        call ncdn_local_flux(flxi, 1, er_3(nr), nr, 1)
                        pfi_3(nr)=flxi*1.0d+20
                        func_gamma(nr)=(pfe_3(nr)-pfi_3(nr))*g(nr)*AEE
                    enddo 

                    do nr=2,nrmax-1
                        term_alpha(nr)=(dea(nr-1)+4.d0*dea(nr)-dea(nr+1))/(4.d0*dr*dr)
                        term_beta(nr)=-2.d0*dea(nr)/(dr*dr)-(dea(nr+1)-dea(nr-1))/(rg(nr)*ra*2.d0*dr)
                        term_gamma(nr)=(-dea(nr-1)+4.d0*dea(nr)+dea(nr+1))/(4.d0*dr*dr)
                        func_diff(nr)=(term_alpha(nr)*er_3(nr-1)+term_beta(nr)*er_3(nr)+term_gamma(nr)*er_3(nr+1))
                    enddo
            
                    SUM_Error_Func_3=0.0d0
                    do nr=2,nrmax-1
                        Error_Func_3(nr)= (EF_normalized_factor**2) * 0.5d0*(func_gamma(nr)+func_diff(nr))**2
                        SUM_Error_Func_3=SUM_Error_Func_3+Error_Func_3(nr)
                    enddo  

                    accelerator_factor_tmp=accelerator_factor_3
                    accelerator_factor_3=accelerator_factor_2
                    accelerator_factor_2=accelerator_factor_tmp

                    SUM_Error_Func_tmp=SUM_Error_Func_3 
                    SUM_Error_Func_3=SUM_Error_Func_2  
                    SUM_Error_Func_2=SUM_Error_Func_tmp            
                    do nr=2,nrmax-1
                        er_tmp(nr)=er_3(nr)
                        Error_Func_tmp(nr)=Error_Func_3(nr)
                    enddo
                    do nr=2,nrmax-1
                        er_3(nr)=er_2(nr)
                        Error_Func_3(nr)=Error_Func_2(nr)
                    enddo
                    do nr=2,nrmax-1
                        er_2(nr)=er_tmp(nr)
                        Error_Func_2(nr)=Error_Func_tmp(nr)
                    enddo                                                          

!                    print'(a8,i5,4e20.10)','ls1: ',it_cnt_linesearch,SUM_Error_Func_1,SUM_Error_Func_2,SUM_Error_Func_3,accelerator_factor_3         
                
                    if(SUM_Error_Func_2<SUM_Error_Func_1)then
                        exit
                    endif
                enddo
            else
                do
                    it_cnt_linesearch=it_cnt_linesearch+1
                    accelerator_factor_3=accelerator_factor_2*2.0d0

                    do nr=2,nrmax-1
                        er_3(nr)=er_1(nr)+grad_SUM_Error_Func(nr)*accelerator_factor_3
                    enddo
                    do nr=2,nrmax-1
                        call ncdn_local_flux(flxe, 1, er_3(nr), nr, 0)
                        pfe_3(nr)=flxe*1.0d+20
                        call ncdn_local_flux(flxi, 1, er_3(nr), nr, 1)
                        pfi_3(nr)=flxi*1.0d+20
                        func_gamma(nr)=(pfe_3(nr)-pfi_3(nr))*g(nr)*AEE
                    enddo 

                    do nr=2,nrmax-1
                        term_alpha(nr)=(dea(nr-1)+4.d0*dea(nr)-dea(nr+1))/(4.d0*dr*dr)
                        term_beta(nr)=-2.d0*dea(nr)/(dr*dr)-(dea(nr+1)-dea(nr-1))/(rg(nr)*ra*2.d0*dr)
                        term_gamma(nr)=(-dea(nr-1)+4.d0*dea(nr)+dea(nr+1))/(4.d0*dr*dr)
                        func_diff(nr)=(term_alpha(nr)*er_3(nr-1)+term_beta(nr)*er_3(nr)+term_gamma(nr)*er_3(nr+1))
                    enddo
            
                    SUM_Error_Func_3=0.0d0
                    do nr=2,nrmax-1
                        Error_Func_3(nr)=(EF_normalized_factor**2) * 0.5d0*(func_gamma(nr)+func_diff(nr))**2
                        SUM_Error_Func_3=SUM_Error_Func_3+Error_Func_3(nr)
                    enddo  

!                    print'(a8,i5,4e20.10)','ls2: ',it_cnt_linesearch,SUM_Error_Func_1,SUM_Error_Func_2,SUM_Error_Func_3,accelerator_factor_3         

                    if(SUM_Error_Func_3>SUM_Error_Func_2)then
                        exit
                    endif
                
                    accelerator_factor_1=accelerator_factor_2
                    accelerator_factor_2=accelerator_factor_3

                    SUM_Error_Func_1=SUM_Error_Func_2
                    SUM_Error_Func_2=SUM_Error_Func_3            
                    do nr=2,nrmax-1
                        er_1(nr)=er_2(nr)
                        Error_Func_1(nr)=Error_Func_2(nr)
                        er_2(nr)=er_3(nr)
                        Error_Func_2(nr)=Error_Func_3(nr)
                    enddo                
                enddo
            endif
            accelerator_factor=(SUM_Error_Func_1*(accelerator_factor_3**2-accelerator_factor_2**2) &
           &                    +SUM_Error_Func_2*(accelerator_factor_1**2-accelerator_factor_3**2) &
           &                    +SUM_Error_Func_3*(accelerator_factor_2**2-accelerator_factor_1**2))/ &
           &             (SUM_Error_Func_1*(accelerator_factor_3-accelerator_factor_2) &
           &                    +SUM_Error_Func_2*(accelerator_factor_1-accelerator_factor_3) &
                               +SUM_Error_Func_3*(accelerator_factor_2-accelerator_factor_1))*0.5d0
!        next Er
            do nr=2,nrmax-1
                er_0(nr)=er_0(nr)+grad_SUM_Error_Func(nr)*accelerator_factor
            enddo

!       Calculate _SUM_Error_Func
            do nr=2,nrmax-1
                call ncdn_local_flux(flxe, 1, er_0(nr), nr, 0)
                pfe(nr)=flxe*1.0d+20
                call ncdn_local_flux(flxi, 1, er_0(nr), nr, 1)
                pfi(nr)=flxi*1.0d+20
                func_gamma(nr)=(pfe(nr)-pfi(nr))*g(nr)*AEE
            enddo 

            do nr=2,nrmax-1
                term_alpha(nr)=(dea(nr-1)+4.d0*dea(nr)-dea(nr+1))/(4.d0*dr*dr)
                term_beta(nr)=-2.d0*dea(nr)/(dr*dr)-(dea(nr+1)-dea(nr-1))/(rg(nr)*ra*2.d0*dr)
                term_gamma(nr)=(-dea(nr-1)+4.d0*dea(nr)+dea(nr+1))/(4.d0*dr*dr)
                func_diff(nr)=(term_alpha(nr)*er_0(nr-1)+term_beta(nr)*er_0(nr)+term_gamma(nr)*er_0(nr+1))
            enddo
            
        
            NEW_SUM_Error_Func=0.0d0
            do nr=2,nrmax-1
                Error_Func(nr)=(EF_normalized_factor**2) * 0.5d0*(func_gamma(nr)+func_diff(nr))**2
                NEW_SUM_Error_Func=NEW_SUM_Error_Func+Error_Func(nr)
            enddo
            
            if(NEW_SUM_Error_Func>SUM_Error_Func_3)then
                do nr=2,nrmax-1
                    er_0(nr)=er_3(nr)
                enddo
                do nr=2,nrmax-1
                    call ncdn_local_flux(flxe, 1, er_0(nr), nr, 0)
                    pfe(nr)=flxe*1.0d+20
                    call ncdn_local_flux(flxi, 1, er_0(nr), nr, 1)
                    pfi(nr)=flxi*1.0d+20
                    func_gamma(nr)=(pfe(nr)-pfi(nr))*g(nr)*AEE
                enddo 

                do nr=2,nrmax-1
                    term_alpha(nr)=(dea(nr-1)+4.d0*dea(nr)-dea(nr+1))/(4.d0*dr*dr)
                    term_beta(nr)=-2.d0*dea(nr)/(dr*dr)-(dea(nr+1)-dea(nr-1))/(rg(nr)*ra*2.d0*dr)
                    term_gamma(nr)=(-dea(nr-1)+4.d0*dea(nr)+dea(nr+1))/(4.d0*dr*dr)
                    func_diff(nr)=(term_alpha(nr)*er_0(nr-1)+term_beta(nr)*er_0(nr)+term_gamma(nr)*er_0(nr+1))
                enddo                
                NEW_SUM_Error_Func=SUM_Error_Func_3
            endif
            relative_error_SUM = (NEW_SUM_Error_Func-SUM_Error_Func)/SUM_Error_Func
            
            open(841,file='er_dynamics_iteration_errorlog.txt',position='APPEND')
                write(841,'(i8,3e20.10)') it_cnt,SUM_Error_Func, NEW_SUM_Error_Func, relative_error_SUM
            close(841)

            open(839,file='er_dynamics_log_iteration_ALL.txt',position='APPEND') 
            do nr=1,nrmax
                write(839,'(2i6,7e20.10)') it_cnt,nr,rg(nr),er_0(nr),Error_Func(nr),func_gamma(nr)*EF_normalized_factor,func_diff(nr)*EF_normalized_factor,grad_SUM_Error_Func(nr),SUM_Error_Func
            enddo
            write(839,*)
            close(839)

            if(mod(it_cnt-1, (it_cnt_max/1000))==0)then
                open(842,file='er_dynamics_log_iteration.txt',position='APPEND') 
                do nr=1,nrmax
                    write(842,'(2i6,7e20.10)') it_cnt,nr,rg(nr),er_0(nr),Error_Func(nr),func_gamma(nr)*EF_normalized_factor,func_diff(nr)*EF_normalized_factor,grad_SUM_Error_Func(nr),SUM_Error_Func
                enddo
                write(842,*)
                close(842)
            endif

           open(840,file='er_dynamics_iteration_lastest_Values.txt',status='replace')
           do nr=1,nrmax      
           write(840,'(i6,9e20.10)') nr,rg(nr),er(nr),er_0(nr),func_gamma(nr)*EF_normalized_factor,func_diff(nr)*EF_normalized_factor,grad_SUM_Error_Func(nr),g(nr),pfe(nr),pfi(nr)
           enddo
           close(840)           

            print'(i8,3e20.10)',it_cnt,SUM_Error_Func, NEW_SUM_Error_Func, relative_error_SUM

            if((dabs(relative_error_SUM)<threshold_relative_error_SUM) .or. (it_cnt>it_cnt_max) .or. (NEW_SUM_Error_Func<threshold_NEW_error_SUM))then
                print*
                print'(a70)','==================== In ER modules, iteration was finished. ===================='
                print'(l2,a25,e20.10,a4,e20.10,a)',(NEW_SUM_Error_Func<threshold_NEW_error_SUM),'NEW_SUM_Error_Func',NEW_SUM_Error_Func,'<',threshold_NEW_error_SUM,'   threshold_NEW_error_SUM'
                print'(l2,a25,e20.10,a4,e20.10,a)',(dabs(relative_error_SUM)<threshold_relative_error_SUM),'relative_error_SUM',relative_error_SUM,'<',threshold_relative_error_SUM,'   threshold_relative_error_SUM'
                print'(l2,a25,i20,a4,i20,a)',(it_cnt>it_cnt_max),'it_cnt',it_cnt,'>',it_cnt_max,'   it_cnt_max in ER module'
                print*
                print'(a43,i4)','DEA_SETTING_CONDITION (0:Cnst. 1:Variable)=',DEA_SETTING_CONDITION
                print'(a12,e20.10,a13,e20.10)','dea(nr=0) = ', dea(1), 'dea(nrmax) = ',dea(nrmax)
                print'(a67,i4)','ER_CALC_INITIAL_VALUE (0:Function, 1:Ambipolar, 2:Experimental ) : ', ER_CALC_INITIAL_VALUE
                
                exit
!            if((it_cnt>it_cnt_max) .or. (NEW_SUM_Error_Func<threshold_NEW_error_SUM))then
!                print*
!                print*,'==================== In ER modules, iteration was finished. ===================='
!                print'(l2,a25,e20.10,a4,e20.10,a)',(NEW_SUM_Error_Func<threshold_NEW_error_SUM),'NEW_SUM_Error_Func',NEW_SUM_Error_Func,'<',threshold_NEW_error_SUM,'   threshold_NEW_error_SUM'
!                print'(l2,a25,e20.10,a4,e20.10,a)',(dabs(relative_error_SUM)<threshold_relative_error_SUM),'relative_error_SUM',relative_error_SUM,'<',threshold_relative_error_SUM,'   threshold_relative_error_SUM'
!                print'(l2,a25,i20,a4,i20,a)',(it_cnt>it_cnt_max),'it_cnt',it_cnt,'>',it_cnt_max,'   it_cnt_max in ER module'
!                exit
            else
                SUM_Error_Func=NEW_SUM_Error_Func
                do nr=2,nrmax-1
                    SUM_Error_Func_1=SUM_Error_Func
                    er_1(nr)=er_0(nr)
                enddo
            endif
!            if(mod(it_cnt,1)==0)then
!            endif
        enddo
        
        do nr=2,nrmax-1
            er(nr)=er_0(nr)
        enddo
        
        call cpu_time(cpu_time_end)
        print*
        print*,'========================================================'
        print*,'CPU TIME @tr_loop', cpu_time_end - cpu_time_start, 'sec'
        print*
        print*,'TEST RUN is Stopped'
        stop
        end select

#ifdef er_dynamics_log_QQ
       close(812)
#endif 
    end subroutine er_dynamics
    
      SUBROUTINE BANDRD( A , X , N , L , LA , IERR )

!          INPUT : A(LA,N) : COEFFICIENT MATRIX
!                  X(N)    : RIGHT-HAND-SIDE VECTOR
!                  N       : MATRIX SIZE
!                  L       : BAND WIDTH (L.LE.LA)
!                  LA      : SIZE OF ARRAY A
!          OUTPUT: X(N)    : SOLUTION VECTOR
!                  IERR    : ERROR CODE : 0 : NORMAL END
!                                         10000 : L IS EVEN
!                                         30000 : SINGULAR MATRIX
!          NOTICE: ARRAY A AND X WILL BE DESTROYED.

      IMPLICIT NONE
      REAL(8), DIMENSION(LA, N),INTENT(INOUT) :: A
      REAL(8), DIMENSION(N),    INTENT(INOUT) :: X
      INTEGER(4),               INTENT(IN)    :: N, L, LA
      INTEGER(4),               INTENT(OUT)   :: IERR
      REAL(8)            :: ABS1, ABS2, TEMP
      REAL(8), PARAMETER :: EPS = 1.D-70
      INTEGER(4)         :: I, J, K, LH, LHM, NM, LHMK, NPMK, LPMI, IPIVOT, IP, JJ

      IF( MOD(L,2) .EQ. 0 ) GO TO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1

      DO K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO I = 1 , LHMK
            LPMI = L+1-I
            A( 1:L-1 , K ) = A( 2:L , K )
            A( L    , K    ) = 0.D0
            A( LPMI , NPMK ) = 0.D0
         ENDDO
      ENDDO

      DO I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = ABS( A(1,IPIVOT) )
         DO K = IP , LH
            ABS1 = ABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
         ENDDO

         IF( ABS2 .LT. EPS ) GO TO 9002
         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO J = 1 , L
               TEMP            = A( J , I      )
               A( J , I      ) = A( J , IPIVOT )
               A( J , IPIVOT ) = TEMP
            ENDDO
         END IF

         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP

         A( 2:L , I ) = A( 2:L , I ) * TEMP

         DO K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            A( 1:L-1 , K ) = A( 2:L , K ) - A( 2:L , I ) * TEMP

            A( L , K ) = 0.D0
         ENDDO
         IF( LH .LT. N ) LH = LH + 1
      ENDDO

      IF( ABS(A(1,N)) .LT. EPS ) GO TO 9002
      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO I = 1 , NM
         K = N-I
         TEMP = 0.D0
         DO J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
         ENDDO
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
      ENDDO

      IERR = 0
      RETURN

 9000 IERR = 10000
      RETURN
 9002 IERR = 30000+I
      RETURN
      END SUBROUTINE BANDRD
	
  end module er_mod

!---------------------------------------------------------------------------------------


  module t3d_er
    use t3d_er_param
    implicit none
    private
    real(8), dimension(1:300,0:10),public,save  :: d,v,chi,chi0,kv,pflx,hflx
!    public :: t3d_er_ambplr, t3d_er_nccoef, t3d_er_show_nccoef,t3d_write_file_evo,t3d_write_file
    public :: t3d_er_ambplr, t3d_er_nccoef, t3d_er_show_nccoef

  contains
    subroutine t3d_er_ambplr(er_tr, nrmax)
      use TRCOMM, only:nt
      use param_mod
      use prof_sh_mod
      use prof_dn_mod
      use er_mod
      use ncsh_mod
      use ncdn_mod

      implicit none
      integer, intent(in) :: nrmax
      real(8), dimension(nrmax), intent(out) :: er_tr

      integer :: nsw
!      real(8), dimension(:,:), allocatable :: pflx,hflx
      real(8), dimension(1:nrmax):: er0

      integer(4),save :: ntsave=-1
      
      !    model = 1 ! SHq
      !    model = 2 ! DCOM/NNW

      call param_set

      select case(transp_model)
      case(1)
         call prof_sh_set
         call ncsh_allocate
      case(2)
         call prof_dn_set
         call ncdn_allocate
      end select


!      if (allocated(pflx)) then
!1         continue
!      else
!         allocate(pflx(imax, 0:num_ions), hflx(imax, 0:num_ions))
!!         allocate(er0(imax))
!      end if
!      if (allocated(d)) then
!         continue
!      else
!         allocate(d(imax, 0:num_ions), v(imax, 0:num_ions), chi(imax, 0:num_ions), kv(imax, 0:num_ions))
!      end if


        select case(transp_model)
        case(1)
           nsw = 0
           er0(1:imax) = 0.0d0
           call er_ambplr(er0, nsw)
           er_tr(1:nrmax) = er(1:nrmax)
        case(2)
            select case(ER_CALC_METHOD)
            case(10)
                if(nt.eq.0 .and. ntsave.lt.0)then
                    nsw = 0
                    er(1:imax)=0.0d0
                    call er_initialize_exp_tmp(er)
                    er_tr(1:nrmax) = er(1:nrmax)
                    ntsave=ntsave+1
                else
                    nsw =  0
                    er_tr(1:nrmax)=er(1:nrmax)
                endif 
            case(20)
                if(nt.eq.0 .and. ntsave.lt.0)then
                    nsw = 0
                    er(1:imax) = 0.0d0
                    call er_ambplr(er, nsw)

                    er_tr(1:nrmax) = er(1:nrmax)
                    ntsave=ntsave+1                    
                else
                    nsw=0
                    call er_ambplr(er, nsw)
                    er_tr(1:nrmax)=er(1:nrmax)
                endif
                
            case(30:39)
                if(nt.eq.0 .and. ntsave.lt.0)then
                    nsw = 0
                    er(1:imax) = 0.0d3
                    select case(ER_CALC_INITIAL_VALUE)
                    case(0)
                        er(1:imax)=0.d3
                    case(1)
                        call er_ambplr(er, nsw)
                    case(2)
                        call er_initialize_exp_tmp(er)
                    case default
                        print*,'ER INITIALIZE ERROR. You must select ER_CALC_INITIAL_VALU = 0 or 1 or 2.'
                        print*,'stopped @',__FILE__,__LINE__
                        stop
                    end select    
                
                    call er_dynamics(nt,er)
                    er_tr(1:nrmax) = er(1:nrmax)
                    ntsave=ntsave+1
                else
                    nsw=0
                    call er_dynamics(nt,er)
                    er_tr(1:nrmax)=er(1:nrmax)
                endif
            case default
            print*,'[ER calculation mesthd] select ERROR [no.1]'
            stop
            end select
        end select
           
      if(ntsave .lt. 0) then 
         print*,'Er calculate Error, Not define method./EXP/AMBIPOLAR/DYNAMICS/'
         stop
      endif
      select case(transp_model)
      case(1)
         call ncsh_analysis(pflx, hflx, d, v, chi)
      case(2)
         call ncdn_analysis(pflx, hflx, d, v,chi0, chi, kv)
      case default
         print *, 'error in main program : model'
      end select
    end subroutine t3d_er_ambplr
      
    subroutine t3d_er_nccoef(nsw)
      use TRCOMM, only : AD, AV, AK, AVK, NRMAX,AKNC,AVNC,ADNC,AVKNC
      use ncdn_mod
      use T3D_FILE_IO      
      integer, intent(in) :: nsw
      integer :: nr, ns
       
    if(nsw==3)then
        call ncdn_analysis(pflx, hflx, d, v, chi0,chi, kv)
            do ns=0, 1 ! sssumed the number of ions =1 
            do nr=1,nrmax
                pflx_save(nr,ns)=pflx(nr,ns)
                hflx_save(nr,ns)=hflx(nr,ns)
                dd_save(nr,ns)=d(nr,ns)
                vv_save(nr,ns)=v(nr,ns)
                chi0_save(nr,ns)=chi0(nr,ns)
                chi_save(nr,ns)=chi(nr,ns)
                kv_save(nr,ns)=kv(nr,ns)
        enddo
        enddo        
    endif

!               ADNC(1:nrmax,1)=d(1:nrmax,0)
!               ADNC(1:nrmax,2)=d(1:nrmax,1)               
!               AVNC(1:nrmax,1)=v(1:nrmax,0)
!               AVNC(1:nrmax,2)=v(1:nrmax,1)               
!               AKNC(1:nrmax,1)=chi(1:nrmax,0)
!               AKNC(1:nrmax,2)=chi(1:nrmax,1) 
!               AVKNC(1:nrmax,1)=kv(1:nrmax,0)
!               AVKNC(1:nrmax,2)=kv(1:nrmax,1) 

      select case(nsw)
      case(1)
            do nr = 1, nrmax
 !              AD(nr,1) = AD(nr,1) + CNP_ANe*d(nr,0)
 !              AD(nr,2) = AD(nr,2) + CNP_ANi*d(nr,1)
                ADNC(nr,1)=d(nr,0)
                ADNC(nr,2)=d(nr,1)               
            end do
      case(2)
            do nr = 1, nrmax
!               AV(nr,1) = AV(nr,1) + CNP_ANe*v(nr,0)
!               AV(nr,2) = AV(nr,2) + CNP_ANi*v(nr,1)
                AVNC(nr,1)=v(nr,0)
                AVNC(nr,2)=v(nr,1)               
            end do
      case(3)
            do nr = 1, nrmax
!               AK(nr,1) = AK(nr,1) + CNH_ANe*chi(nr,0)
!               AK(nr,2) = AK(nr,2) + CNH_ANi*chi(nr,1)
!!!!!print'(a120)','NC_CHIに関して特殊処理を実行：r/a<0.15に関して常にr/a=0.15の値を使用中．注意せよ．'
!!!!!if(nr<10)then
!!!!!             AKNC(nr,1)=chi(10,0)
!!!!!             AKNC(nr,2)=chi(10,1)
!!!!!else
               AKNC(nr,1)=chi(nr,0)
               AKNC(nr,2)=chi(nr,1)
!!!!!endif               
!               print'(4e20.10)',aknc(nr,1),chi(nr,0),aknc(nr,2),chi(nr,1)               
            end do
      case(4)
            do nr = 1, nrmax
               AVKNC(nr,1)=kv(nr,0)
               AVKNC(nr,2)=kv(nr,1) 
!               AVK(nr,1) = AVK(nr,1) + CNH_ANe*kv(nr,0)
!               AVK(nr,2) = AVK(nr,2) + CNH_ANi*kv(nr,1)
            end do
      end select
    end subroutine t3d_er_nccoef

    subroutine t3d_er_show_nccoef(d_out, v_out, chi_out,kv_out,pflx_out,hflx_out)
      use TRCOMM, only : NRMAX,er
      use ncdn_mod
      implicit none
      real(8), dimension(nrmax,2), intent(out) :: d_out, v_out, chi_out,kv_out,pflx_out,hflx_out
      integer :: ns
      integer i

!      do k = 0,1 
!         do i = 1, nrmax
!            call ncdn_local_flux(local_flx, 2, er(i), i, k)
!            hflx_tmp(i,k) = local_flx
!         end do
!      end do
      do i=1,nrmax
      do ns = 1, 2
         d_out(i,ns) = d(i,ns-1)
         v_out(i,ns) = v(i,ns-1)
         chi_out(i,ns) = chi(i,ns-1)
         kv_out(i,ns) = kv(i,ns-1)
         pflx_out(i,ns) = pflx(i,ns-1)
         hflx_out(i,ns) = hflx(i,ns-1)
!         h_out(1:nrmax,ns)=hflx_tmp(1:nrmax,ns-1)
      end do
      enddo
!      print'(a60)','++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!      do i=1,nrmax
!         print'(i4,2e20.10)',i,d_out(i,1),d_out(i,2)
!      enddo 

!      do ns = 1, 2
!         d_out(1:nrmax,ns) = d(1:nrmax,ns-1)
!         v_out(1:nrmax,ns) = v(1:nrmax,ns-1)
!         chi_out(1:nrmax,ns) = chi(1:nrmax,ns-1)
!         kv_out(1:nrmax,ns) = kv(1:nrmax,ns-1)
!         pflx_out(1:nrmax,ns) = pflx(1:nrmax,ns-1)
!         hflx_out(1:nrmax,ns) = hflx(1:nrmax,ns-1)
!!         h_out(1:nrmax,ns)=hflx_tmp(1:nrmax,ns-1)
!      end do
    end subroutine t3d_er_show_nccoef

  subroutine er_initialize_exp_tmp(er)
    use TRCOMM, only:NRMAX,KUFDCG
    real(8) :: er(nrmax)
  
!    er(1)=	5.0293d1
!    er(2)=	5.4688d1
!    er(3)=	7.3730d1
!    er(4)=	1.2500d2
!    er(5)=	2.3311d2
!    er(6)=	4.2969d2
!    er(7)=	7.5342d2
!    er(8)=	1.2500d3
!    er(9)=	1.9948d3
!    er(10)=	2.9394d3
!    er(11)=	3.9822d3
!    er(12)=	5.0502d3
!    er(13)=	6.0602d3
!    er(14)=	6.9249d3
!    er(15)=	7.5603d3
!    er(16)=	7.8935d3
!    er(17)=	7.8717d3
!    er(18)=	7.4724d3
!    er(19)=	6.7145d3
!    er(20)=	5.6697d3
!    er(21)=	4.4760d3
!    er(22)=	3.3514d3
!    er(23)=	2.6085d3
!    er(24)=	2.3331d3
!    er(25)=	2.2290d3
!    er(26)=	2.2009d3
!    er(27)=	2.2335d3
!    er(28)=	2.3142d3
!    er(29)=	2.4323d3
!    er(30)=	2.5791d3
!    er(31)=	2.7476d3
!    er(32)=	2.9323d3
!    er(33)=	3.1290d3
!    er(34)=	3.3346d3
!    er(35)=	3.5468d3
!    er(36)=	3.7643d3
!    er(37)=	3.9862d3
!    er(38)=	4.2123d3
!    er(39)=	4.4426d3
!    er(40)=	4.6773d3
!    er(41)=	4.9170d3
!    er(42)=	5.1621d3
!    er(43)=	5.4133d3
!    er(44)=	5.6710d3
!    er(45)=	5.9356d3
!    er(46)=	6.2074d3
!    er(47)=	6.4864d3
!    er(48)=	6.7725d3
!    er(49)=	7.0657d3
!    er(50)=	7.3654d3
!    er(51)=	7.6713d3
!    er(52)=	7.9829d3
!    er(53)=	8.2995d3
!    er(54)=	8.6207d3
!    er(55)=	8.9463d3
!    er(56)=	9.2760d3
!    er(57)=	9.6102d3
!    er(58)=	9.9496d3
!    er(59)=	1.0296d4
!    er(60)=	1.0650d4



    if(KUFDCG=='032940')then
    if(nrmax == 60) then   
    er(1)=	1.41d3
    er(2)=	2.57d3
    er(3)=	3.57d3
    er(4)=	4.40d3
    er(5)=	5.09d3
    er(6)=	5.64d3
    er(7)=	6.06d3
    er(8)=	6.37d3
    er(9)=	6.58d3
    er(10)=	6.70d3
    er(11)=	6.74d3
    er(12)=	6.71d3
    er(13)=	6.61d3
    er(14)=	6.46d3
    er(15)=	6.27d3
    er(16)=	6.04d3
    er(17)=	5.78d3
    er(18)=	5.50d3
    er(19)=	5.21d3
    er(20)=	4.90d3
    er(21)=	4.60d3
    er(22)=	4.30d3
    er(23)=	4.01d3
    er(24)=	3.73d3
    er(25)=	3.47d3
    er(26)=	3.23d3
    er(27)=	3.02d3
    er(28)=	2.84d3
    er(29)=	2.70d3
    er(30)=	2.59d3
    er(31)=	2.51d3
    er(32)=	2.48d3
    er(33)=	2.48d3
    er(34)=	2.53d3
    er(35)=	2.62d3
    er(36)=	2.75d3
    er(37)=	2.92d3
    er(38)=	3.14d3
    er(39)=	3.39d3
    er(40)=	3.69d3
    er(41)=	4.02d3
    er(42)=	4.38d3
    er(43)=	4.78d3
    er(44)=	5.20d3
    er(45)=	5.66d3
    er(46)=	6.13d3
    er(47)=	6.62d3
    er(48)=	7.12d3
    er(49)=	7.63d3
    er(50)=	8.14d3
    er(51)=	8.65d3
    er(52)=	9.15d3
    er(53)=	9.64d3
    er(54)=	1.01d4
    er(55)=	1.05d4
    er(56)=	1.09d4
    er(57)=	1.13d4
    er(58)=	1.15d4
    er(59)=	1.18d4
    er(60)=	1.19d4    
    elseif (nrmax == 120) then
er(1)=	7.515d2
er(2)=	1.406d3
er(3)=	2.014d3
er(4)=	2.576d3
er(5)=	3.094d3
er(6)=	3.570d3
er(7)=	4.006d3
er(8)=	4.403d3
er(9)=	4.763d3
er(10)=	5.087d3
er(11)=	5.378d3
er(12)=	5.636d3
er(13)=	5.863d3
er(14)=	6.060d3
er(15)=	6.229d3
er(16)=	6.372d3
er(17)=	6.489d3
er(18)=	6.582d3
er(19)=	6.652d3
er(20)=	6.701d3
er(21)=	6.730d3
er(22)=	6.739d3
er(23)=	6.731d3
er(24)=	6.706d3
er(25)=	6.666d3
er(26)=	6.611d3
er(27)=	6.543d3
er(28)=	6.463d3
er(29)=	6.371d3
er(30)=	6.270d3
er(31)=	6.159d3
er(32)=	6.040d3
er(33)=	5.914d3
er(34)=	5.782d3
er(35)=	5.644d3
er(36)=	5.502d3
er(37)=	5.356d3
er(38)=	5.208d3
er(39)=	5.057d3
er(40)=	4.905d3
er(41)=	4.752d3
er(42)=	4.600d3
er(43)=	4.449d3
er(44)=	4.299d3
er(45)=	4.152d3
er(46)=	4.008d3
er(47)=	3.867d3
er(48)=	3.730d3
er(49)=	3.598d3
er(50)=	3.471d3
er(51)=	3.350d3
er(52)=	3.235d3
er(53)=	3.126d3
er(54)=	3.025d3
er(55)=	2.931d3
er(56)=	2.845d3
er(57)=	2.767d3
er(58)=	2.698d3
er(59)=	2.638d3
er(60)=	2.587d3
er(61)=	2.545d3
er(62)=	2.513d3
er(63)=	2.491d3
er(64)=	2.478d3
er(65)=	2.476d3
er(66)=	2.484d3
er(67)=	2.503d3
er(68)=	2.532d3
er(69)=	2.572d3
er(70)=	2.622d3
er(71)=	2.682d3
er(72)=	2.754d3
er(73)=	2.835d3
er(74)=	2.927d3
er(75)=	3.030d3
er(76)=	3.142d3
er(77)=	3.265d3
er(78)=	3.398d3
er(79)=	3.540d3
er(80)=	3.691d3
er(81)=	3.852d3
er(82)=	4.022d3
er(83)=	4.200d3
er(84)=	4.387d3
er(85)=	4.582d3
er(86)=	4.784d3
er(87)=	4.994d3
er(88)=	5.210d3
er(89)=	5.432d3
er(90)=	5.661d3
er(91)=	5.895d3
er(92)=	6.134d3
er(93)=	6.377d3
er(94)=	6.624d3
er(95)=	6.874d3
er(96)=	7.126d3
er(97)=	7.381d3
er(98)=	7.637d3
er(99)=	7.893d3
er(100)=	8.150d3
er(101)=	8.405d3
er(102)=	8.659d3
er(103)=	8.911d3
er(104)=	9.159d3
er(105)=	9.404d3
er(106)=	9.643d3
er(107)=	9.877d3
er(108)=	1.010d4
er(109)=	1.032d4
er(110)=	1.053d4
er(111)=	1.073d4
er(112)=	1.093d4
er(113)=	1.110d4
er(114)=	1.127d4
er(115)=	1.142d4
er(116)=	1.156d4
er(117)=	1.168d4
er(118)=	1.178d4
er(119)=	1.186d4
er(120)=	1.193d4
    else
    print*,'Er EXP setteing error'
    stop
    endif
    endif

    if(KUFDCG=='109695')then
    if (nrmax == 120) then !dummy nrmax=120 && EXP_ER
er(1)=	-660.5882089
er(2)=	-682.6613623
er(3)=	-719.34519
er(4)=	-770.4830731
er(5)=	-835.8567894
er(6)=	-915.1875505
er(7)=	-1008.137328
er(8)=	-1114.310464
er(9)=	-1233.255546
er(10)=	-1364.467565
er(11)=	-1507.390308
er(12)=	-1661.419012
er(13)=	-1825.903236
er(14)=	-2000.149964
er(15)=	-2183.42689
er(16)=	-2374.96592
er(17)=	-2573.966819
er(18)=	-2779.601016
er(19)=	-2991.015559
er(20)=	-3207.337174
er(21)=	-3427.676412
er(22)=	-3651.131897
er(23)=	-3876.794617
er(24)=	-4103.752235
er(25)=	-4331.093457
er(26)=	-4557.912356
er(27)=	-4783.31268
er(28)=	-5006.412132
er(29)=	-5226.34656
er(30)=	-5442.274064
er(31)=	-5653.379011
er(32)=	-5858.87591
er(33)=	-6058.01314
er(34)=	-6250.076533
er(35)=	-6434.392765
er(36)=	-6610.332548
er(37)=	-6777.313634
er(38)=	-6934.803573
er(39)=	-7082.322245
er(40)=	-7219.444152
er(41)=	-7345.800448
er(42)=	-7461.080692
er(43)=	-7565.03436
er(44)=	-7657.472051
er(45)=	-7738.266424
er(46)=	-7807.352862
er(47)=	-7864.729837
er(48)=	-7910.459001
er(49)=	-7944.665003
er(50)=	-7967.535026
er(51)=	-7979.318053
er(52)=	-7980.323879
er(53)=	-7970.921872
er(54)=	-7951.539483
er(55)=	-7922.660536
er(56)=	-7884.823301
er(57)=	-7838.618356
er(58)=	-7784.68627
er(59)=	-7723.715105
er(60)=	-7656.437773
er(61)=	-7583.629237
er(62)=	-7506.103597
er(63)=	-7424.71107
er(64)=	-7340.334866
er(65)=	-7253.887992
er(66)=	-7166.309991
er(67)=	-7078.563617
er(68)=	-6991.631476
er(69)=	-6906.512621
er(70)=	-6824.219105
er(71)=	-6745.772513
er(72)=	-6672.200441
er(73)=	-6604.53293
er(74)=	-6543.798832
er(75)=	-6491.022099
er(76)=	-6447.217947
er(77)=	-6413.38888
er(78)=	-6390.520523
er(79)=	-6379.577195
er(80)=	-6381.497194
er(81)=	-6397.187675
er(82)=	-6427.519076
er(83)=	-6473.318957
er(84)=	-6535.365158
er(85)=	-6614.378142
er(86)=	-6711.012363
er(87)=	-6825.846507
er(88)=	-6959.372421
er(89)=	-7111.982509
er(90)=	-7283.955372
er(91)=	-7475.439462
er(92)=	-7686.434427
er(93)=	-7916.76987
er(94)=	-8166.081197
er(95)=	-8433.78216
er(96)=	-8719.033687
er(97)=	-9020.708634
er(98)=	-9337.351899
er(99)=	-9667.13543
er(100)=	-10007.8076
er(101)=	-10356.63631
er(102)=	-10710.34519
er(103)=	-11065.04221
er(104)=	-11416.14002
er(105)=	-11758.26707
er(106)=	-12085.16887
er(107)=	-12389.59828
er(108)=	-12663.19399
er(109)=	-12896.34609
er(110)=	-13078.04767
er(111)=	-13195.73121
er(112)=	-13235.08855
er(113)=	-13179.87311
er(114)=	-13011.68295
er(115)=	-12709.72309
er(116)=	-12250.54565
er(117)=	-11607.766
er(118)=	-10751.7532
er(119)=	-9649.292764
er(120)=	-8263.22
  else  if (nrmax == 121) then !dummy nrmax=120 && AMBIER*2

er(1)=	-74.693464 
er(2)=	-146.247455 
er(3)=	-211.753902 
er(4)=	-272.406231 
er(5)=	-342.903274 
er(6)=	-424.754483 
er(7)=	-510.830921 
er(8)=	-595.669654 
er(9)=	-678.355963 
er(10)=	-761.651625 
er(11)=	-843.266311 
er(12)=	-922.392014 
er(13)=	-1000.890973 
er(14)=	-1077.787943 
er(15)=	-1151.869260 
er(16)=	-1223.819576 
er(17)=	-1293.530395 
er(18)=	-1359.960663 
er(19)=	-1423.313338 
er(20)=	-1483.407937 
er(21)=	-1540.938539 
er(22)=	-1597.078502 
er(23)=	-1650.081288 
er(24)=	-1698.849648 
er(25)=	-1745.590191 
er(26)=	-1789.178804 
er(27)=	-1829.331546 
er(28)=	-1868.235334 
er(29)=	-1903.852508 
er(30)=	-1936.067084 
er(31)=	-1966.725948 
er(32)=	-1994.887175 
er(33)=	-2019.655940 
er(34)=	-2041.686858 
er(35)=	-2061.756952 
er(36)=	-2078.806484 
er(37)=	-2092.733938 
er(38)=	-2104.169948 
er(39)=	-2112.821500 
er(40)=	-2118.496746 
er(41)=	-2121.807054 
er(42)=	-2122.942832 
er(43)=	-2120.135376 
er(44)=	-2114.760758 
er(45)=	-2108.408878 
er(46)=	-2099.579088 
er(47)=	-2088.525744 
er(48)=	-2075.821162 
er(49)=	-2062.396674 
er(50)=	-2048.467082 
er(51)=	-2034.022484 
er(52)=	-2020.290744 
er(53)=	-2006.874424 
er(54)=	-1993.679352 
er(55)=	-1983.140089 
er(56)=	-1975.954027 
er(57)=	-1970.104698 
er(58)=	-1967.357389 
er(59)=	-1969.253544 
er(60)=	-1975.899991 
er(61)=	-2014.512016 
er(62)=	-2089.415734 
er(63)=	-2174.951066 
er(64)=	-2264.134736 
er(65)=	-2355.003914 
er(66)=	-2448.919824 
er(67)=	-2544.604658 
er(68)=	-2642.766994 
er(69)=	-2744.927688 
er(70)=	-2849.585286 
er(71)=	-2957.833580 
er(72)=	-3071.042460 
er(73)=	-3187.375392 
er(74)=	-3306.578604 
er(75)=	-3430.145900 
er(76)=	-3557.912758 
er(77)=	-3689.241060 
er(78)=	-3824.587938 
er(79)=	-3964.587048 
er(80)=	-4108.188490 
er(81)=	-4255.014552 
er(82)=	-4407.089932 
er(83)=	-4563.231930 
er(84)=	-4722.161342 
er(85)=	-4885.193486 
er(86)=	-5051.517690 
er(87)=	-5220.738898 
er(88)=	-5393.303832 
er(89)=	-5569.289420 
er(90)=	-5749.002302 
er(91)=	-5932.336438 
er(92)=	-6119.022884 
er(93)=	-6308.399520 
er(94)=	-6500.452500 
er(95)=	-6696.719990 
er(96)=	-6896.749300 
er(97)=	-7098.804846 
er(98)=	-7304.131802 
er(99)=	-7512.946958 
er(100)=	-7723.818780 
er(101)=	-7937.183856 
er(102)=	-8153.321638 
er(103)=	-8372.754594 
er(104)=	-8595.888122 
er(105)=	-8823.064586 
er(106)=	-9055.500310 
er(107)=	-9294.183796 
er(108)=	-9540.378756 
er(109)=	-9796.010990 
er(110)=	-10064.502468 
er(111)=	-10351.924214 
er(112)=	-10666.337528 
er(113)=	-11019.494534 
er(114)=	-11429.028496 
er(115)=	-11921.518074 
er(116)=	-12543.169886 
er(117)=	-13379.443180 
er(118)=	-14600.279650 
er(119)=	-16605.242858 
er(120)=	-20221.223040 
    else
    print*,'Er EXP setteing error'
    stop
    endif
    endif

    
    end subroutine
  end module t3d_er
