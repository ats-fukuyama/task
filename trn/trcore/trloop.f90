MODULE trloop

! This module advances time upto ntmax*dt

  USE trcomm,ONLY: rkind, ikind

  PUBLIC tr_loop,tr_save_pvprev,tr_calc_dpdrho2jtot
  PRIVATE

CONTAINS

  SUBROUTINE tr_loop

    USE trcomm, ONLY: rkev,ntmax,nrmax,nsamax,t,dt,ntstep,ngtstp, &
         mdluf,ntlmax
    USE trbpsd, ONLY: tr_bpsd_set,tr_bpsd_get
    USE trcalc1, ONLY: tr_calc1
    USE trstep, ONLY: tr_step
    USE trresult, ONLY: tr_status,tr_calc_global,tr_save_ngt,tr_exp_compare
    USE trcalv, ONLY: tr_calc_variables
    IMPLICIT NONE

    REAL(4)        :: t1,t2
    INTEGER(ikind) :: nt,ierr

    CALL GUTIME(t1)

    time_evolution: DO nt = 1, ntmax ! main loop

       CALL tr_save_pvprev

       ! incremental addtion
       t = t + dt
       
       CALL tr_calc1

       ! non-linear iteration
       CALL tr_step(ierr); IF(ierr /= 0) EXIT

       ! calculate associated variables and save values
       CALL tr_calc_bpqpj

       IF(MOD(nt,ntstep) == 0 .OR. &
          MOD(nt,ngtstp) == 0) CALL tr_calc_global
       IF(MOD(nt,ntstep) == 0) CALL tr_status
       IF(MOD(nt,ntstep) == 0 .AND. mdluf > 0) CALL tr_exp_compare
       IF(MOD(nt,ngtstp) == 0) CALL tr_save_ngt

       CALL tr_bpsd_set(ierr)
       IF(ierr /= 0)THEN
          WRITE(6,*) 'XX tr_loop: Error to set the variables to BPSD interface. NT= ', nt
       END IF

       IF(mdluf==2 .OR. mdluf==3)THEN
          IF(nt==ntlmax)THEN
             WRITE(6,*) ' - End of experimental data. Stop calculation. NT= ',nt
             EXIT
          END IF
       END IF
    END DO time_evolution

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

  SUBROUTINE tr_calc_bpqpj
! --------------------------------------------------------------------------
!  calculate following conversions:  d psi/d rho --> bp,qp,jtor,jtot
! --------------------------------------------------------------------------
    USE trcomm, ONLY: pi,rmu0,nrmax,RR,ar1rho,ttrho,rmjrho,arrho,dvrho,  &
         abb1rho,abrho,abvrho,rhog,dpdrho,rdpvrho,qp,q0,qa,bp,rip,       &
         jtot,joh,jtor,jbs_nc,jex_nc,jcd_nb,jcd_ec,jcd_lh,jcd_ic,        &
         id_neq,mdlijq
!         ,nrd3,nrd4

    IMPLICIT NONE
    INTEGER(ikind) :: nr, mdlid
    REAL(rkind) :: FCTR    ! function defined in TASK/lib
!    REAL(rkind) :: FCTR4pt ! function defined in TASK/lib
    REAL(rkind) :: deriv3  ! function defined in TASK/lib
!    REAL(rkind) :: ipl,dr
    REAL(rkind),DIMENSION(0:nrmax) :: factor1,factor2

    factor1(0:nrmax) = dvrho(0:nrmax)*abrho(0:nrmax)*dpdrho(0:nrmax)
    factor2(0:nrmax) = factor1(0:nrmax)/ttrho(0:nrmax)

    ! dpdrho --> bp
    bp(0:nrmax) = ar1rho(0:nrmax)*dpdrho(0:nrmax)/rmjrho(0:nrmax)

    IF(id_neq(1) == 0)THEN

       mdlid = MOD(mdlijq,2)
       SELECT CASE(mdlid)
       CASE(0) ! in the case (qp --> dpdrho)
          DO nr = 1, nrmax
             ! dpdrho --> jtot(j_para)
             jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                        * deriv3(nr,rhog,factor2,nrmax,0)
          END DO
!          jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
          jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))
          
       CASE(1) ! in the case (jtot --> dpdrho)
          ! dpdrho --> qp
          qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                        /(4.d0*pi**2 * dpdrho(1:nrmax))
!          qp(0)       = FCTR4pt(rhog(1),rhog(2),rhog(3),qp(1),qp(2),qp(3))
          qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))
!          qp(0) = qp(1)
          q0 = qp(0)
          qa = qp(nrmax)

       END SELECT

    ELSE IF(id_neq(1) == 2)THEN
       ! dpdrho --> qp
       qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                     /(4.d0*pi**2 * dpdrho(1:nrmax))
!       qp(0)       = FCTR4pt(rhog(1),rhog(2),rhog(3),qp(1),qp(2),qp(3))
       qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))
!       qp(0) = qp(1)
       q0 = qp(0)
       qa = qp(nrmax)

       DO nr = 1, nrmax
          ! dpdrho --> jtot(j_para)
          jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                     * deriv3(nr,rhog,factor2,nrmax,0)
       END DO
!       jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
       jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))

    END IF

    ! only for graphic output -----------------------------------------
    rdpvrho(0:nrmax) = ttrho(0:nrmax)*arrho(0:nrmax)          &
                      /(4.d0*pi**2*qp(0:nrmax)) ! d psi/d V 
    DO nr = 1, nrmax
          ! dpdrho --> jtor(j_toroidal)    
          jtor(nr) = rmjrho(nr)/(rmu0*dvrho(nr))               &
                     * deriv3(nr,rhog,factor1,nrmax,0)

    END DO
!    jtor(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtor(1),jtor(2),jtor(3))
    jtor(0) = FCTR(rhog(1),rhog(2),jtor(1),jtor(2))
    ! -----------------------------------------------------------------


    joh(0:nrmax) = jtot(0:nrmax)                &
         - (jbs_nc(0:nrmax) + jex_nc(0:nrmax) + &
            jcd_nb(0:nrmax) + jcd_ec(0:nrmax) + &
            jcd_ic(0:nrmax) + jcd_lh(0:nrmax))

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
  END SUBROUTINE tr_calc_bpqpj

END MODULE trloop
