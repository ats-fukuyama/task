MODULE trsetup

! This module setup table for computation and 
! initializes the plasma profiles
  USE trcomm, ONLY: ikind, rkind

  PUBLIC tr_setup
  PRIVATE

CONTAINS

  SUBROUTINE tr_setup

    USE trcomm, ONLY: &
         tr_nit_allocate,tr_nsa_allocate,tr_nr_allocate,tr_ngt_allocate, &
         t,t_prev,ngt,kidnsa,ns_nsa,idnsa,nrmax,nsamax,pa,pz,pz0,        &
         nitmax,modelg,rip,rips,vtor,vpol,vpar,vprp

    USE trbpsd, ONLY: tr_bpsd_init
    USE eq_bpsd_mod, ONLY: eq_bpsd_init
    USE trloop, ONLY: tr_save_pvprev
    USE trresult, ONLY: tr_calc_global, tr_save_ngt
    USE trcalv, ONLY: tr_calc_zeff, tr_calc_clseta
    USE trcalc, ONLY: tr_calc_source
    IMPLICIT NONE
    INTEGER(ikind):: ierr

    CALL tr_bpsd_init(ierr)
    CALL eq_bpsd_init(ierr)

    t      = 0.d0           ! time is initialized
    t_prev = 0.d0
    ngt    = -1             ! save data count is initilaized
    nitmax = 0              ! iteration count is initialized
    rip    = rips

! +++ setup of equations and speices information +++
    CALL tr_set_idneq

    CALL tr_nsa_allocate

    CALL tr_set_species


! +++ Basic setup +++
    CALL tr_ngt_allocate    ! allocation for data save

    CALL tr_setup_table     ! calculate table for matrix generation

    CALL tr_nr_allocate     ! allocation for radial profile

! +++ Initialization for successive interactive calculation +++
    vtor(0:nrmax) = 0.d0
    vpol(0:nrmax) = 0.d0
    vpar(0:nrmax) = 0.d0
    vprp(0:nrmax) = 0.d0

! +++ setup of initial geometic factor and profiles +++
    ! Calculate initial geometric factor for cylindrical assumption
    CALL tr_setup_metric_init

    ! Calculate initial profile.( preparation for calling equilibrium code )
    CALL tr_setup_profile

    ! set geometric factor (EQ/EQDSK file, calc EQU, calc EQ, calc VMEC)
    IF(modelg >= 3) CALL tr_setup_geometric


! +++
    CALL tr_nit_allocate    ! allocation for diagnostics of iteration

    CALL tr_save_pvprev

    CALL tr_calc_zeff      ! *** dummy calculation for tr_calc_global 
    CALL tr_calc_clseta    !      and preparation of 'eta' for NCLASS calc
    CALL tr_calc_source    ! *** 

    CALL tr_calc_global

    CALL tr_save_ngt

    CALL tr_bpsd_init(ierr)

    RETURN
  END SUBROUTINE tr_setup

! ***************************************************************************
! ***************************************************************************

  SUBROUTINE tr_set_idneq
! --------------------------------------------------------------------------
!   This subroutine sets the table for matrix generation.
! --------------------------------------------------------------------------
    USE trcomm, ONLY: tr_neq_allocate, &
         nsamax,neqmax,nvmax,nrmax,nsa_neq,nva_neq,id_neq

    IMPLICIT NONE
    INTEGER(ikind) :: i,neq,nsa

    neqmax = 1 + 3*nsamax
    nvmax  = neqmax*(nrmax+1)

    CALL tr_neq_allocate

    !   nsa_neq = 0 : magnetic field
    !     otherwise : particle species

    !   nva_neq = 1 : density
    !             2 : toroidal velocity
    !             3 : temperature

    !   id_neq  = 0 : equation is not solved in any radius
    !             1 : flat on axis and fixed at plasma surface
    !             2 : fixed to zero on axis and fixed at plasma surface
    !             3 : flat on axis and fixed scale length at plasma surface
    !             4 : fixed to zero on axis and fixed scale length at surface  
    ! magnetic diffusion equation
    neq = 1
    nsa_neq(neq) = 0
    nva_neq(neq) = 0
!    id_neq(neq)  = 2
    id_neq(neq)  = 0

    DO nsa=1,nsamax
       DO i=1,3
          neq = neq+1
          nsa_neq(neq) = nsa
          nva_neq(neq) = i
          SELECT CASE(i)
          CASE(1) ! density
             id_neq(neq) = 0
          CASE(2) ! toroidal velocity
             id_neq(neq) = 0
          CASE(3) ! temperature
             id_neq(neq) = 1
          END SELECT
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE tr_set_idneq


  SUBROUTINE tr_set_species
! --------------------------------------------------------------------------
!   This subroutine sets the table for matrix generation.
! --------------------------------------------------------------------------
    USE trcomm, ONLY: nsamax,pz,pz0,pa,idnsa,kidnsa,ns_nsa
    IMPLICIT NONE

    INTEGER(ikind) :: nsa,ns

    DO nsa=1,nsamax
       ns=ns_nsa(nsa)
       IF(NINT(pz0(ns)) == -1) THEN
          idnsa(nsa) = -1
       ELSE
          IF(NINT(pz(ns)) == 0) THEN
             idnsa(nsa) = 0
          ELSE
             idnsa(nsa) = 1
          END IF
       END IF
    END DO
          
    DO nsa=1,nsamax
       ns=ns_nsa(nsa)
       SELECT CASE(NINT(pz0(ns)))
       CASE(-1)
          kidnsa(nsa)='E'
       CASE(1)
          SELECT CASE(NINT(pa(ns)))
          CASE(1)
             kidnsa(nsa)='H'
          CASE(2)
             kidnsa(nsa)='D'
          CASE(3)
             kidnsa(nsa)='T'
          CASE DEFAULT
             kidnsa(nsa)=' '
          END SELECT
       CASE(2)
          kidnsa(nsa)='A'
       CASE DEFAULT
          kidnsa(nsa)=' '
       END SELECT
    END DO
  
    RETURN
  END SUBROUTINE tr_set_species


  SUBROUTINE tr_setup_table
! --------------------------------------------------------------------------
!   This subroutine sets the table for matrix generation.
! --------------------------------------------------------------------------
    USE trcomm, ONLY: nrmax,nsamax,neqmax,neqrmax,nvrmax,nvmax, &
                      nsa_neq,nva_neq,id_neq,id_neqnr,neq_neqr,       &
                      rg,rg_fixed,neqr_neq,tr_neqr_allocate
    IMPLICIT NONE
    INTEGER(ikind):: neq,nsa,nva,neqr,nr

    neqr = 0
    DO neq = 1, neqmax
       SELECT CASE(id_neq(neq))
!   On EACH grid :
!         id_neqnr = 0 : fixed to zero
!                    1 : solved normally
!                    2 : fixed to a given value
!                    3 : fixed to a given scale length

       CASE(0) ! equation is not solved in any radius (fixed to zero)
          neqr_neq(neq) = 0
          id_neqnr(neq,0:nrmax) = 2
       CASE(1,11) ! flat on axis and fixed at plasma surface
          neqr = neqr+1
          neqr_neq(neq) = neqr
          id_neqnr(neq,0:nrmax-1) = 1
          id_neqnr(neq,    nrmax) = 2
       CASE(2,12) ! fixed to zero on axis and fixed at plasma surface
          neqr = neqr+1
          neqr_neq(neq) = neqr
          id_neqnr(neq,        0) = 0
          id_neqnr(neq,1:nrmax-1) = 1
          id_neqnr(neq,    nrmax) = 2
       CASE(3,13) ! flat on axis and fixed scale length at plasma surface
          neqr = neqr+1
          neqr_neq(neq) = neqr
          id_neqnr(neq,0:nrmax-1) = 1
          id_neqnr(neq,    nrmax) = 3
       CASE(4,14) ! fixed to zero on axis and fixed scale length at surface
          neqr = neqr+1
          neqr_neq(neq) = neqr
          id_neqnr(neq,        0) = 0
          id_neqnr(neq,1:nrmax-1) = 1
          id_neqnr(neq,    nrmax) = 3
       CASE DEFAULT
          WRITE(6,*) 'XX tr_setup_table: undefied id_neq:'
          WRITE(6,*) '   neq,id_neq(neq)=',neq,id_neq(neq)
          STOP
       END SELECT

       IF(id_neq(neq) >= 10) THEN  ! fixed for rg >= rg_fixed
          DO nr=1,nrmax
             nsa=nsa_neq(neq)
             nva=nva_neq(neq)
             IF(rg(nr) >= rg_fixed(nva,nsa)) THEN
                id_neqnr(neq,nr)=2
             END IF
          END DO
       END IF

    END DO
    neqrmax = neqr
    nvrmax  = neqrmax*(nrmax+1)
!    write(*,*) neqrmax

    CALL tr_neqr_allocate

    DO neq=1,neqmax
       IF(neqr_neq(neq) /= 0) neq_neqr(neqr_neq(neq))=neq
    ENDDO
          
    RETURN
  END SUBROUTINE tr_setup_table

! *************************************************************************
! *************************************************************************
  SUBROUTINE tr_setup_metric_init
! --------------------------------------------------------------------------
!   This subroutine initializes geometric factors.
!
!   In the case of calling equilibrium code, however, its substitution is 
!    carried out for making interim profiles for equilibrium code.
! --------------------------------------------------------------------------
    USE trcomm, ONLY: pi,rkap,nrmax,ra,rr,bb,rkap,rg,rm,rjcb,rhog,rhom, &
         rhoa,ttrho,dvrho,abrho,abvrho,arrho,ar1rho,ar2rho,rmjrho,      &
         rmnrho,rmnrhom,rkprho,rkprhom,rhog,rhom,epsrho,                &
         abb2rho,aib2rho,abb1rho,pvolrho,psurrho  ! ,nrd1,nrd2,nrd3

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: dr

    ! --- Interim substitution ---
    dr   = SQRT(rkap)*ra/dble(nrmax)
    rhoa = 1.D0
    DO nr = 0, nrmax
       rg(nr)=dble(nr)*dr
       IF(nr /= 0) rm(nr)=0.5d0*(rg(nr-1)+rg(nr))

       ! normalized variables
       rjcb(nr)    = 1.d0/(SQRT(rkap)*ra)
       rhog(nr)    = rg(nr)*rjcb(nr)
       IF(nr /= 0) rhom(nr) = 0.5d0*(rhog(nr-1)+rhog(nr))   
    END DO

    ! --- cylindrical assumption ---
    DO nr = 0, nrmax
       ar1rho(nr)  = 1.d0/(SQRT(rkap)*ra)          ! const
       ar2rho(nr)  = 1.d0/(SQRT(rkap)*ra)**2       ! const
       abrho(nr)   = 1.d0/(SQRT(rkap)*ra*rr)**2    ! const
       rmjrho(nr)  = rr                            ! const [m]
       rmnrho(nr)  = SQRT(rkap)*ra*rhog(nr) ! [m]
       IF(nr /= 0) rmnrhom(nr)=0.5d0*(rmnrho(nr-1)+rmnrho(nr))
       rkprho(nr)  = rkap
       IF(nr /= 0) rkprhom(nr)=0.5d0*(rkprho(nr-1)+rkprho(nr))
!
       pvolrho(nr) = pi*rkap*(ra*rhog(nr))**2*2.d0*pi*rr
       psurrho(nr) = pi*(rkap+1.d0)*ra*rhog(nr)*2.d0*pi*rr
       dvrho(nr)   = 2.d0*pi*rkap*ra**2*2.d0*pi*rr*rhog(nr)
!
       epsrho(nr)  = rmnrho(nr)/rmjrho(nr)

       ! toroidal field for now
!       abb1rho(nr) = BB*(1.d0 + 0.5d0*epsrho(nr)**2) ! <B>
       abb1rho(nr) = BB
!       abb2rho(nr) = BB**2 *(1.d0+1.5d0*epsrho(nr)**2)
       abb2rho(nr) = BB**2
!       aib2rho(nr) = (1.d0+1.5d0*epsrho(nr)**2)/BB**2
       aib2rho(nr) = 1/BB**2
!       ttrho(nr)   = abb1rho(nr) * rr
       ttrho(nr)   = BB * rr
!       arrho(nr)   = 1.d0/rr**2 * (1+1.5d0*epsrho(nr)**2)            
       arrho(nr)   = 1.d0/rr**2                    ! const
!       abb2rho(nr) = 
       abvrho(nr)  = dvrho(nr)**2 * abrho(nr)

    END DO

  END SUBROUTINE tr_setup_metric_init


  SUBROUTINE tr_setup_geometric
! --------------------------------------------------------------------------
!
! --------------------------------------------------------------------------
    USE trcomm, ONLY: pi,nrmax,ra,rr,rkap,bb,rg,rm,rhoa,                   &
       ttrho,dvrho,abrho,abvrho,arrho,ar1rho,ar2rho,rmjrho,rmnrho,rmnrhom, &
       rkprho,rkprhom,rjcb,rhog,rhom,epsrho,abb2rho,abb1rho,pvolrho,       &
       psurrho,modelg,knameq !, nrd1,nrd2

    USE trbpsd, ONLY: tr_bpsd_set,tr_bpsd_get
    USE equnit_mod, ONLY: eq_parm,eq_prof,eq_calc,eq_load
!    USE equunit_mod, ONLY: equ_prof,equ_calc
!    USE pl_vmec_mod, ONLY: pl_vmec

    IMPLICIT NONE
    INTEGER(ikind) :: ierr
    CHARACTER(len=80) :: line

    CALL tr_bpsd_set(ierr)

    SELECT CASE(modelg)
    CASE(3,5) ! TASK/EQ,EQDSK output geometry
       write(line,'(A,I5)') 'nrmax=',nrmax+1
       call eq_parm(2,line,ierr)
       write(line,'(A,I5)') 'nthmax=',64
       call eq_parm(2,line,ierr)
       write(line,'(A,I5)') 'nsumax=',0
       call eq_parm(2,line,ierr)

       call eq_load(modelg,knameq,ierr) ! load eq data and calculate eq
       IF(ierr.NE.0) THEN
          WRITE(6,*) 'XX eq_load: ierr=',ierr
          RETURN
       ENDIF
       call tr_bpsd_get(ierr)

!    CASE(7)
!       call pl_vmec(knameq,ierr) ! load vmec data
!       call tr_bpsd_get(ierr)
!       call trgout
!    CASE(8) ! CALL TASK/EQU
!       call equ_prof             ! initial calculation of eq
!       call equ_calc             ! recalculate eq
!       call tr_bpsd_get(ierr)
!       call trgout        
    CASE(9) ! CALL TASK/EQ
          CALL eq_calc              ! recalculate eq
          CALL tr_bpsd_get(ierr)
          ! --- here the convergence of q profile must be confirmed
          ! ---  and show the graph of Psi(R,Z)
!          write(*,*) rhog(0:nrmax)
          CALL tr_setup_profile
    END SELECT

  END SUBROUTINE tr_setup_geometric

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_setup_profile
! -------------------------------------------------------------------------
!   This subroutine calculates inital profiles.
! -------------------------------------------------------------------------
    USE trcomm, ONLY: pi,rkev,rkap,rdlt,nrmax,nsmax,nsamax,ns_nsa, &
         rg,rm,rhog,rhom,ra,rr,rn,ru,rt,rp,rp_tot,          &
         ttrho,dvrho,arrho,ar1rho,rdpvrho,dpdrho,           &
         mdluf,jtot,joh,jbs_nc,jex_nc
    USE trloop, ONLY: tr_calc_dpdrho2j
    USE plprof, ONLY: pl_prof2,pl_qprf
                      
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nsmax):: rn_ns,ru_ns,rtpr_ns,rtpp_ns
    INTEGER(ikind):: nr,nsa,ns

    IF(mdluf == 0)THEN
       rp_tot(0:nrmax) = 0.d0
       DO nr=0,nrmax
          CALL pl_prof2(rhog(nr),rn_ns,rtpr_ns,rtpp_ns,ru_ns)
          DO nsa=1,nsamax
             ns=ns_nsa(nsa)
             rn(nsa,nr)=rn_ns(ns)
             ru(nsa,nr)=ru_ns(ns)
             rt(nsa,nr)=(rtpr_ns(ns)+2.D0*rtpp_ns(ns))/3.D0

             ! the pressure of each species
             rp(nsa,nr)  = rn(nsa,nr)*1.d20 * rt(nsa,nr)*rkev
             rp_tot(nr)  = rp_tot(nr) + rp(nsa,nr)
             
          ! MDLUF = 0 : trprof.f90
!          pex(nsa,nr) = 0.d0
!          sex(nsa,nr) = 0.d0
!          prf(nsa,nr) = 0.d0

!          pbm(nr)     = 0.d0
!          wrot(nr)    = 0.d0
!          vtor(nr)    = 0.d0
          ENDDO
       ENDDO

       ! initialization of current profiles
       jtot(0:nrmax)   = 0.d0
       joh(0:nrmax)    = 0.d0
       jbs_nc(0:nrmax) = 0.d0
       jex_nc(0:nrmax) = 0.d0
       ! set profile of 'd psi/d rho' from given jtot profile
       CALL tr_prof_j2dpdrho
       CALL tr_calc_dpdrho2j

    END IF

    RETURN
  END SUBROUTINE tr_setup_profile


  SUBROUTINE tr_prof_j2dpdrho
! -------------------------------------------------------------------------
!   This subroutine gives initial profile of toroidal current density 
!    and d psi/d rho.
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rmu0,pi,nrmax,BB,RR,ra,rkap,q0,qa,         &
         rg,rhog,abb1rho,dvrho,ttrho,abrho,arrho,ar1rho,abvrho,  &
         profj1,profj2,rip,rips,                                 &
         dpdrho,rdpvrho,jtot,joh,eta,knameq! ,nrd1,nrd2

    IMPLICIT NONE
    REAL(rkind) :: dr,factor0,factorp,factorm,fact,dpdrhos
    INTEGER(ikind) :: nr

    DO nr = 0, nrmax
       IF(((SQRT(rkap)*ra)**ABS(profj1)-rg(nr)**ABS(profj1)).LE.0.d0) THEN
          jtot(nr) = 0.D0
       ELSE
          jtot(nr) = &
          ((SQRT(rkap)*ra)**ABS(profj1)-rg(nr)**ABS(profj1))**ABS(profj2)
       ENDIF
    ENDDO

    dpdrho(0:nrmax)  = 0.d0
    rdpvrho(0:nrmax) = 0.d0
    DO nr = 1, nrmax
       dr      = rhog(nr)-rhog(nr-1)
       factor0 = rmu0*0.5d0*(abb1rho(nr)+abb1rho(nr-1)) &
                     *0.5d0*(dvrho(nr)+dvrho(nr-1))     &
                     *0.5d0*(jtot(nr)+jtot(nr-1))       &
                  /(0.5d0*(ttrho(nr)+ttrho(nr-1)))**2
       factorp = abvrho(nr  )/ttrho(nr  )
       factorm = abvrho(nr-1)/ttrho(nr-1)

       rdpvrho(nr) = (factorm*rdpvrho(nr-1) + factor0*dr)/factorp
       dpdrho(nr)  = rdpvrho(nr)*dvrho(nr)
    END DO
!    rdpvrho(nr) = ttrho(nr)*arrho(nr)/(4.d0*pi**2*qp(nr))! d psi/d V

    ! set the boundary value of dpdrho in terms of plasma current value
    dpdrhos = 2.d0*pi*rmu0*rip*1.d6 / (dvrho(nrmax)*abrho(nrmax))
!    dpdrhos = 2.d0*pi*rmu0*rip*1.d6*dvrho(nrmax)/abvrho(nrmax)
    ! correction in terms of the boundary value of dpdrho
    fact = dpdrhos / dpdrho(nrmax)
!    write(*,*) fact

    dpdrho(0:nrmax)  = fact*dpdrho(0:nrmax) 
    rdpvrho(0:nrmax) = fact*rdpvrho(0:nrmax)

    jtot(0:nrmax) = fact*jtot(0:nrmax)
    joh(0:nrmax)  = jtot(0:nrmax)

    eta(0:nrmax) = 1.d-7
    
    ! diagnostic
!    nrd1(0:nrmax) = rdpvrho(0:nrmax)
!    nrd2(0:nrmax) = jtot(0:nrmax)

!    write(*,*) jtot(0:nrmax)

    RETURN
  END SUBROUTINE tr_prof_j2dpdrho

END MODULE trsetup
