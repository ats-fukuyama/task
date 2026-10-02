MODULE trcalc1
  USE trcomm, ONLY: ikind,rkind

  PRIVATE
  PUBLIC tr_calc1,tr_calc_geometry,tr_calc_dpdrho

CONTAINS

  SUBROUTINE tr_calc1
    USE trcomm,ONLY: t,dt,ntmax,nteqit,mdluf,mdlgmt,mdlglb,nrmax, &
         vtor,wrot,rmjrho
    USE trbpsd,ONLY: tr_bpsd_set,tr_bpsd_get
    USE trsource,ONLY: tr_source1
    USE trufin,ONLY: tr_ufin_density,tr_ufin_rotation,tr_ufin_temperature, &
                     tr_ufin_geometry,tr_ufin_global
    USE equnit_mod,ONLY: eq_calc,eq_load
!    USE equunit_mod,ONLY: equ_calc
    IMPLICIT NONE

    INTEGER(ikind) :: nt,ierr
    REAL(rkind)    :: time

    time = t
    nt   = INT(time/dt)

    ! read experimental profile data at every time step
    IF(mdluf == 2)THEN
       CALL tr_ufin_density(time,2,ierr)
       CALL tr_ufin_rotation(time,2,ierr)
       CALL tr_ufin_temperature(time,2,ierr)
    END IF


    SELECT CASE(mdlglb)
    CASE(6)
       CALL tr_ufin_global(time,0,ierr)
    CASE(7)
       CALL tr_ufin_global(time,1,ierr)
    CASE(9)
       ! read equilibrium code ouput
    END SELECT
    

    SELECT CASE(mdlgmt)
    CASE(6)
       CALL tr_ufin_geometry(time,0,ierr)
    CASE(7)
       CALL tr_ufin_geometry(time,1,ierr)
    CASE(9)
       ! Interaction with equilibrium codes
       IF(nteqit /= 0 .AND. MOD(nt, nteqit) == 0)THEN

          CALL tr_bpsd_set(ierr)
          IF(ierr /= 0)THEN
             WRITE(6,*) &
       'XX tr_calc1: Error to set the variables to BPSD interface. NT= ', nt
          END IF

          CALL eq_calc
             
          CALL tr_bpsd_get(ierr)
          IF(ierr /= 0)THEN
             WRITE(6,*) &
       'XX tr_calc1: Error to get the variables from BPSD interface. NT= ', nt
          END IF
       END IF
       ! cofirmation of the conservation of nV', nTV'^(5/3)

    CASE DEFAULT
       CONTINUE
    END SELECT


    CALL tr_calc_geometry

    CALL tr_calc_dpdrho

    CALL tr_source1 ! source calculation (unnecessary non-linear iteration)

    ! Interim way; 
    !  meaningful calculation only if 'wrot' is obtained from exp. data.
    vtor(0:nrmax) = wrot(0:nrmax) * rmjrho(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc1

! ************************************************************************

  SUBROUTINE tr_calc_geometry
! -------------------------------------------------------------------------
!
! -------------------------------------------------------------------------
    USE trcomm, ONLY: nrmax,BB,RR,dvrho,abrho,rmjrho,rmnrho, &
         epsrho,abb1rho,abb2rho,aib2rho,ttrho,arrho,abvrho,ar2rho,mdlgmt

    IMPLICIT NONE
    INTEGER(ikind) :: nr

    IF(mdlgmt==0) RETURN

    ! *** associated values ***
    DO nr = 0, nrmax
       epsrho(nr)  = rmnrho(nr)/rmjrho(nr)

       ! toroidal field for now
       !       abb1rho(nr) = BB*(1.d0 + 0.5d0*epsrho(nr)**2) ! <B>
       abb1rho(nr) = BB
       !       abb2rho(nr) = BB**2 *(1.d0+1.5d0*epsrho(nr)**2)
       abb2rho(nr) = BB**2
       !       aib2rho(nr) = (1.d0+1.5d0*epsrho(nr)**2)/BB**2
       aib2rho(nr) = 1/BB**2

       !       ttrho(nr)   = abb1rho(nr) * rr ??? what's the definition ???
!       ttrho(nr)   = BB*RR
       ttrho(nr)   = abb1rho(nr) * rr
       !       arrho(nr)   = 1.d0/rr**2 * (1+1.5d0*epsrho(nr)**2)
!       arrho(nr)   = 1.d0/RR*2
       arrho(nr)   = 1.d0/rmjrho(nr)**2                    ! const

       abrho(nr)   = ar2rho(nr) * arrho(nr)
       abvrho(nr)  = dvrho(nr)**2 * abrho(nr)
    END DO
    
    RETURN
  END SUBROUTINE tr_calc_geometry


  SUBROUTINE tr_calc_dpdrho
! -------------------------------------------------------------------------
!
! -------------------------------------------------------------------------
    USE trcomm,ONLY: rmu0,pi,id_neq,t,dt,ntmax,nrmax,rhog,rip,rips,ripe, &
         mdluf,mdlijq,abb1rho,dvrho,ttrho,abrho,arrho,ar1rho,abvrho,     &
         rmjrho,dpdrho,rdpvrho,jtot,qp,bp
    USE trufin,ONLY: tr_ufin_field
    IMPLICIT NONE
    
    INTEGER(ikind) :: nt,nr,id,ierr
    REAL(rkind)    :: FCTR,DERIV3 ! the function in TASK/lib
    REAL(rkind)    :: time,drip,dr,dpdrhos,factor0,factor0p,factor0m, &
                      factorp,factorm,fact,rdpvrhomax,dpdrhomax
    REAL(rkind),DIMENSION(0:nrmax) :: factor1, factor2

    id   = id_neq(1)
    time = t

    exp: IF(mdluf == 0)THEN ! --------------------------------------------

       equation_id1: SELECT CASE(id)
       CASE(0)
          CONTINUE ! constant RIP
       CASE(2)
          ! incremental addition of total current
          drip = (ripe - rips) / ntmax
          rip  = rip + drip
          nt   = INT(time/dt)
          IF(nt==ntmax)THEN
             rip  = ripe
             rips = ripe
          END IF

       END SELECT equation_id1

       dpdrho(nrmax) = 2.d0*pi*rmu0*rip*1.d6/(dvrho(nrmax)*abrho(nrmax))

    ELSE IF(mdluf > 0)THEN ! ---------------------------------------------
       equation_id2: SELECT CASE(id)
       CASE(0) ! not solve
          CALL tr_ufin_field(time,2,mdlijq,ierr)

          IF(MOD(mdlijq,2)==1)THEN ! jtot --> dpdrho
             dpdrho(0:nrmax)  = 0.d0
             rdpvrho(0:nrmax) = 0.d0 ! not defined at nr=0 in this case
             DO nr = 1, nrmax
                dr      = rhog(nr)-rhog(nr-1)
                factor0p=rmu0*abb1rho(nr)*dvrho(nr)*jtot(nr)/ttrho(nr)**2
                factor0m=rmu0*abb1rho(nr-1)*dvrho(nr-1)*jtot(nr-1)/ttrho(nr-1)**2
                factor0 = 0.5d0*(factor0p + factor0m)
                factorp = abvrho(nr  )/ttrho(nr  )
                factorm = abvrho(nr-1)/ttrho(nr-1)

                rdpvrho(nr) = (factorm*rdpvrho(nr-1) + factor0*dr)/factorp
                dpdrho(nr)  = rdpvrho(nr) * dvrho(nr)
             END DO

             IF(mdlijq==1)THEN
                ! set the boundary value of dpdrho in terms of RIP
                dpdrhos = 2.d0*pi*rmu0*rip*1.d6 / (dvrho(nrmax)*abrho(nrmax))
!                dpdrhos = 2.d0*pi*rmu0*rip*1.d6*dvrho(nrmax)/abvrho(nrmax)

                ! correction in terms of the boundary value of dpdrho
                fact = dpdrhos / dpdrho(nrmax)
!                write(*,*) fact

                dpdrho(0:nrmax)  = fact*dpdrho(0:nrmax)
                rdpvrho(0:nrmax) = fact*rdpvrho(0:nrmax) ! d psi/d V
                jtot(0:nrmax)    = fact*jtot(0:nrmax)
             ELSE IF(mdlijq==3)THEN
                rip = dpdrho(nrmax)*dvrho(nrmax)*abrho(nrmax) &
                     /(2.d0*pi*rmu0*1.d6)
             END IF

             ! dpdrho --> qp                                                
             qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                          /(4.d0*pi**2 * dpdrho(1:nrmax))
             qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))

             rdpvrho(0) = ttrho(0)*arrho(0)/(4.d0*pi**2*qp(0))


          ELSE IF(MOD(mdlijq,2)==0)THEN ! qp --> dpdrho
             dpdrho(0:nrmax) = ttrho(0:nrmax)*arrho(0:nrmax)*dvrho(0:nrmax) &
                             / (4.d0*pi**2 * qp(0:nrmax))

             IF(mdlijq==2)THEN
                ! set the boundary value of dpdrho in terms of RIP
                dpdrhos = 2.d0*pi*rmu0*rip*1.d6 / (dvrho(nrmax)*abrho(nrmax))
!                dpdrhos = 2.d0*pi*rmu0*rip*1.d6*dvrho(nrmax)/abvrho(nrmax)

                ! correction in terms of the boundary value of dpdrho
                fact = dpdrhos / dpdrho(nrmax)
!                write(*,*) fact

                dpdrho(0:nrmax) = fact*dpdrho(0:nrmax)
                qp(0:nrmax)     = qp(0:nrmax) / fact
             ELSE IF(mdlijq==4)THEN
                rip = dpdrho(nrmax)*dvrho(nrmax)*abrho(nrmax) &
                     /(2.d0*pi*rmu0*1.d6)
             END IF
             rdpvrho(0:nrmax) = ttrho(0:nrmax)*arrho(0:nrmax) &
                                /(4.d0*pi**2*qp(0:nrmax)) ! d psi/d V

             ! dpdrho --> jtot
             factor1(0:nrmax) = dvrho(0:nrmax)*abrho(0:nrmax)*dpdrho(0:nrmax)
             factor2(0:nrmax) = factor1(0:nrmax)/ttrho(0:nrmax)
             DO nr = 1, nrmax
                ! dpdrho --> jtot(j_para)
                jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                          * deriv3(nr,rhog,factor2,nrmax,0)
             END DO
!             jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
             jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))

          END IF

          ! dpdrho --> bp
          bp(0:nrmax) = ar1rho(0:nrmax)*dpdrho(0:nrmax)/rmjrho(0:nrmax)


       CASE(2) ! fixed to zero on axis and fixed at plasma surface
          CALL tr_ufin_field(time,2,mdlijq,ierr)

          IF(mdlijq <= 2)THEN ! rip --> dpdrho at the surface
             dpdrho(nrmax) = 2.d0*pi*rmu0*rip*1.d6/(dvrho(nrmax)*abrho(nrmax))

          ELSE IF(mdlijq == 3)THEN ! jtot --> dpdrho at the surface
             dpdrhomax  = 0.d0
             rdpvrhomax = 0.d0 ! not defined at nr=0 in this case
             DO nr = 1, nrmax
                dr      = rhog(nr)-rhog(nr-1)
                factor0p=rmu0*abb1rho(nr)*dvrho(nr)*jtot(nr)/ttrho(nr)**2
                factor0m=rmu0*abb1rho(nr-1)*dvrho(nr-1)*jtot(nr-1)/ttrho(nr-1)**2
                factor0 = 0.5d0*(factor0p + factor0m)
                factorp = abvrho(nr  )/ttrho(nr  )
                factorm = abvrho(nr-1)/ttrho(nr-1)

                rdpvrhomax = (factorm*rdpvrhomax + factor0*dr)/factorp
                dpdrhomax  = rdpvrhomax * dvrho(nr)
             END DO
             dpdrho(nrmax) = dpdrhomax
             rip = dpdrho(nrmax)*dvrho(nrmax)*abrho(nrmax) &
                   /(2.d0*pi*rmu0*1.d6)

          ELSE IF(mdlijq == 4)THEN ! qp --> dpdrho at the surfac
             dpdrho(nrmax) = ttrho(nrmax)*arrho(nrmax)*dvrho(nrmax) &
                            / (4.d0*pi**2 * qp(nrmax))
             rip = dpdrho(nrmax)*dvrho(nrmax)*abrho(nrmax) &
                   /(2.d0*pi*rmu0*1.d6)
          END IF

       END SELECT equation_id2
    END IF exp ! ---------------------------------------------------------

    RETURN
  END SUBROUTINE tr_calc_dpdrho

END MODULE trcalc1
