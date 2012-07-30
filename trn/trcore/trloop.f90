MODULE trloop

! This module advances time upto ntmax*dt

  USE trcomm,ONLY: rkind, ikind

  PUBLIC tr_loop,tr_save_pvprev,tr_calc_dpdrho2j
  PRIVATE

CONTAINS

  SUBROUTINE tr_loop

    USE trcomm, ONLY: rkev,ntmax,nrmax,nsamax,t,dt,ntstep,ngtstp, &
         modelg,nteqit,rip,rips,ripe,rn,rt,rp,rp_tot
    USE trbpsd, ONLY: tr_bpsd_set,tr_bpsd_get
    USE trstep, ONLY: tr_step
    USE trresult, ONLY: tr_status,tr_calc_global,tr_save_ngt
    USE equnit_mod, ONLY: eq_calc,eq_load
!    USE equunit_mod, ONLY: equ_calc
    IMPLICIT NONE

    REAL(4)        :: t1,t2,drip
    INTEGER(ikind) :: nt,ierr

    CALL GUTIME(t1)

    DO nt = 1, ntmax ! main loop

       CALL tr_save_pvprev

       ! incremental addtions
       t=t+dt

       drip = (ripe - rips) / ntmax
       rip  = rip + drip
       IF(nt==ntmax)THEN
          rip  = ripe
          rips = ripe
       END IF
       ! set the boundary value of dpdrho in terms of plasma current value
       CALL tr_mag_boundary
!       write(*,*) 'rip = ',rip

       ! trcalc1.f90

       CALL tr_step(ierr); IF(ierr /= 0) GO TO 9000

       ! Interaction with EQ
       IF(nteqit /= 0 .AND. MOD(nt, nteqit) == 0)THEN
          SELECT CASE(modelg)
          CASE(8)
!!$          CALL tr_bpsd_set(ierr)
!!$          CALL equ_calc
!!$          CALL tr_bpsd_get(ierr)
          CASE(9)
             CALL tr_bpsd_set(ierr)
             CALL eq_calc
             CALL tr_bpsd_get(ierr)
          CASE DEFAULT
          END SELECT
       END IF
       ! cofirmation of the conservation of nV', nTV'^(5/3)

       CALL tr_calc_dpdrho2j

       IF(MOD(nt,ntstep) == 0 .OR. &
          MOD(nt,ngtstp) == 0) CALL tr_calc_global
       IF(MOD(nt,ntstep) == 0) CALL tr_status
       IF(MOD(nt,ngtstp) == 0) CALL tr_save_ngt

       CALL tr_bpsd_set(ierr)
       
    END DO

9000 CONTINUE

    CALL GUTIME(t2)

    write(6,'(A,F12.3,A)') '## cputime = ',T2 - T1,' (s)'

  END SUBROUTINE tr_loop

! **************************************************************************

  SUBROUTINE tr_save_pvprev
! ----------------------------------------------------
!      save plasma variables to previous values
! ----------------------------------------------------
    USE trcomm,ONLY: ikind,rkind,nrmax,nsamax, &
         dpdrho,rn,ru,rt,dpdrho_prev,rn_prev,ru_prev,rt_prev
    IMPLICIT NONE

    dpdrho_prev(0:nrmax)=dpdrho(0:nrmax)
    rn_prev(1:nsamax,0:nrmax)=rn(1:nsamax,0:nrmax)
    ru_prev(1:nsamax,0:nrmax)=ru(1:nsamax,0:nrmax)
    rt_prev(1:nsamax,0:nrmax)=rt(1:nsamax,0:nrmax)

    RETURN
  END SUBROUTINE tr_save_pvprev

! **************************************************************************

  SUBROUTINE tr_mag_boundary
    USE trcomm, ONLY: pi,rmu0,nrmax,dvrho,abrho,rip,dpdrho

    dpdrho(nrmax) = 2.d0*pi*rmu0*rip*1.d6 / (dvrho(nrmax)*abrho(nrmax))

    RETURN
  END SUBROUTINE tr_mag_boundary

! **************************************************************************

  SUBROUTINE tr_calc_dpdrho2j
! --------------------------------------------------------------------------
!  calculate following conversions:  d psi/d rho --> bp,qp,jtor,jtot
! --------------------------------------------------------------------------
    USE trcomm, ONLY: pi,rmu0,nrmax,RR,ar1rho,ttrho,rmjrho,arrho,dvrho,  &
         abb1rho,abrho,abvrho,rhog,dpdrho,rdpvrho,qp,q0,qa,bp,rip,       &
         jtot,joh,jtor,jbs_nc,jex_nc,jcd_nb,jcd_ec,jcd_lh,jcd_ic
!         ,nrd3,nrd4

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: FCTR    ! function defined in TASK/lib
    REAL(rkind) :: deriv3  ! function defined in TASK/lib
!    REAL(rkind) :: ipl,dr
    REAL(rkind),DIMENSION(0:nrmax) :: factor1,factor2

    ! dpdrho --> bp
    bp(0:nrmax) = ar1rho(0:nrmax)*dpdrho(0:nrmax)/rmjrho(0:nrmax)

    ! dpdrho --> qp
    qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                  /(4.d0*pi**2 * dpdrho(1:nrmax))
!    qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)    &
!                  /(4.d0*pi**2 * rdpvrho(1:nrmax))
    !   * FCTR4pt: func. in TASK/lib
!    qp(0)       = FCTR4pt(rhog(1),rhog(2),rhog(3),qp(1),qp(2),qp(3))
    qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))
    q0 = qp(0)
    qa = qp(nrmax)


    factor1(0:nrmax) = dvrho(0:nrmax)*abrho(0:nrmax)*dpdrho(0:nrmax)
    factor2(0:nrmax) = factor1(0:nrmax)/ttrho(0:nrmax)
    DO nr = 1, nrmax
       ! dpdrho --> jtor(j_toroidal)    
       jtor(nr) = rmjrho(nr)/(rmu0*dvrho(nr))               &
                  * deriv3(nr,rhog,factor1,nrmax,0)
       ! dpdrho --> jtot(j_para)
       jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                  * deriv3(nr,rhog,factor2,nrmax,0)
    END DO
!    jtor(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtor(1),jtor(2),jtor(3))
    jtor(0) = FCTR(rhog(1),rhog(2),jtor(1),jtor(2))
!    jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
    jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))


    joh(0:nrmax) = jtot(0:nrmax) - jbs_nc(0:nrmax) - jex_nc(0:nrmax)

    ! ***** inverse conversion for confirmation *****
!!$    ipl = 0.d0
!!$    DO nr = 1, nrmax
!!$       dr = rhog(nr) - rhog(nr-1)
!!$       ipl = ipl + 0.5d0*(dvrho(nr)+dvrho(nr-1))  &
!!$                  *0.5d0*(jtor(nr)+jtor(nr-1))*dr
!!$    END DO
!!$    ipl = ipl/(2.d0*pi*RR)
!!$    write(*,*) rip,ipl
    ! ***********************************************

!    nrd3(0:nrmax) = dpdrho(0:nrmax)
!    nrd4(0:nrmax) = jtot(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc_dpdrho2j


!!$  SUBROUTINE tr_eq_calc
!!$
!!$    DO
!!$    ! iteration for convergence of q(safety factor)
!!$       CALL eq_prof
!!$       CALL eq_calc
!!$       CALL tr_bpsd_get(ierr)
!!$
!!$       !convergence check of q
!!$    END DO
!!$
!!$    CALL tr_prof_q2dpdrho
!!$
!!$  END SUBROUTINE tr_eq_calc
!!$
!!$  SUBROUTINE tr_prof_q2dpdrho
!!$    USE trcomm, ONLY: rmu0,pi,nrmax,BB,RR,q0,qa,profj1,profj2,knameq, &
!!$         &         rhog,abb1rho,dvrho,ttrho,abrho,arrho,dpdrho,ar1rho,&
!!$         &         jtot,joh,bp,qp, nrd1,nrd2,nrd3,nrd4
!!$
!!$    IMPLICIT NONE
!!$    REAL(rkind) :: dr,prof,factor0,sumfact1
!!$    REAL(rkind) :: FCTR ! function defined in TASK/lib
!!$    INTEGER(ikind) :: nr
!!$
!!$! q --> dpdrho, j
!!$
!!$    dpdrho(
!!$
!!$  END SUBROUTINE tr_prof_q2dpdrho  

END MODULE trloop
