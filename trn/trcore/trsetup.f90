MODULE trsetup

! This module setup table for computation and 
! initializes the plasma profiles
  USE trcomm, ONLY: ikind, rkind

  PUBLIC  tr_setup
  PRIVATE

  REAL(rkind) :: time

CONTAINS

  SUBROUTINE tr_setup

    USE trcomm, ONLY: &
         tr_nit_allocate,tr_nsa_allocate,tr_nr_allocate,tr_ngt_allocate, &
         t,t_prev,ngt,ns_nsa,idnsa,nrmax,nsamax,nsafmax,   &
         pa,pz,pz0,nitmax,rip,rips

         

    USE trbpsd, ONLY: tr_bpsd_init
    USE eq_bpsd_mod, ONLY: eq_bpsd_init
    USE trloop, ONLY: tr_save_pvprev
    USE trcalv, ONLY: tr_calc_variables
    USE trresult, ONLY: tr_calc_global, tr_save_ngt
    USE trsource, ONLY:tr_source1,tr_source2
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
    CALL tr_nsa_allocate

    CALL tr_set_species

    CALL tr_set_idneq

    CALL tr_setup_table     ! calculate table for matrix generation

    CALL tr_print_species

! +++ allocation of array and initilization +++
    CALL tr_ngt_allocate    ! allocation for data save
    CALL tr_nr_allocate     ! allocation for radial profile
    CALL tr_nit_allocate    ! allocation for diagnostics of iteration

    CALL tr_setup_varinit

! +++ setup of initial profiles and metric quantities +++
    CALL tr_setup_profile

    CALL tr_setup_geometry

    CALL tr_setup_field

! +++ calculate global quantities +++
    CALL tr_calc_variables
    CALL tr_source1
    CALL tr_source2

    CALL tr_save_pvprev
    CALL tr_calc_global
    CALL tr_save_ngt


    CALL tr_bpsd_init(ierr)

    RETURN
  END SUBROUTINE tr_setup

! ***************************************************************************
! ***************************************************************************

  SUBROUTINE tr_set_species
! --------------------------------------------------------------------------
!   Setup cenversion (fast ion --> bulk ion) table
! --------------------------------------------------------------------------
    USE trcomm, ONLY: tr_neq_allocate, &
         nsm,nsamax,nrmax,neqmax,nvmax,nsabmax,nsafmax,nsanmax,    &
         pz,pz0,pa,idion,idnsa,kidns,ns_nsa,idion,nsa_neq,nva_neq, &
         nsab_nsa,nsaf_nsa,nsan_nsa,nsab_nsaf
    IMPLICIT NONE

    CHARACTER(len=1) :: kidnsf
    INTEGER(ikind)   :: nsa,ns,nsa1,ns1,nva,neq

    !   nsa_neq = 0 : magnetic field
    !     otherwise : particle species

    !   nva_neq = 1 : density
    !             2 : toroidal velocity
    !             3 : energy

    neqmax = 1 + 3*nsamax
    nvmax  = neqmax*(nrmax+1)

    CALL tr_neq_allocate

    neq = 1
    nsa_neq(neq) = 0
    nva_neq(neq) = 0
    DO nsa = 1, nsamax
       DO nva = 1, 3
          neq = neq + 1
          nsa_neq(neq) = nsa
          nva_neq(neq) = nva
       ENDDO
    ENDDO

    ! ----------------------------------------------------------------

    ns_nsa(1:nsm)    = 0

    DO nsa = 1, nsamax
       ns_nsa(nsa) = nsa
    END DO

    ! ----------------------------------------------------------------

    nsabmax = 0
    nsafmax = 0
    nsanmax = 0

    nsab_nsa(1:nsm)  = 0
    nsaf_nsa(1:nsm)  = 0
    nsan_nsa(1:nsm)  = 0
    nsab_nsaf(1:nsm) = 0

    DO nsa = 1, nsamax
       ns=ns_nsa(nsa)
       IF(kidns(ns) == ' ') CYCLE ! exclude for dummy species (see plinit)

       IF(NINT(pz0(ns)) == -1) THEN
          idnsa(nsa) = -1 ! electron
          nsabmax = nsabmax + 1
          nsab_nsa(nsa) = nsabmax
       ELSE
          IF(NINT(pz(ns)) == 0) THEN
             idnsa(nsa) = 0 ! neutral particle
             nsanmax = nsanmax + 1
             nsan_nsa(nsa) = nsanmax
          ELSE
             IF(idion(ns) == 0)THEN
                idnsa(nsa) = 1 ! bulk ion particle
                nsabmax = nsabmax + 1
                nsab_nsa(nsa) = nsabmax
             ELSE IF(idion(ns) == 1)THEN
                idnsa(nsa) = 2 ! fast ion particle
                nsafmax = nsafmax + 1
                nsaf_nsa(nsa) = nsafmax

                kidnsf = kidns(ns)
                DO nsa1 = 1, nsamax
                   ns1 = ns_nsa(nsa1)
                   IF(ns /= ns1 .AND. kidnsf == kidns(ns1))THEN
                      nsab_nsaf(nsa) = nsa1
                   END IF
                END DO

             END IF

          END IF
       END IF
    END DO

!    write(*,*)'nsabmax',nsabmax, 'nsafmax',nsafmax,'nsanmax',nsanmax
          
    RETURN
  END SUBROUTINE tr_set_species


  SUBROUTINE tr_set_idneq
! --------------------------------------------------------------------------
!   This subroutine sets the table for matrix generation.
! --------------------------------------------------------------------------
    USE trcomm, ONLY: neqmax,idnsa,nsa_neq,nva_neq,id_neq

    IMPLICIT NONE
    INTEGER(ikind) :: i,neq,nsa,nva

    !   id_neq  = 0 : equation is not solved in any radius
    !             1 : flat on axis and fixed at plasma surface
    !             2 : fixed to zero on axis and fixed at plasma surface
    !             3 : flat on axis and fixed scale length at plasma surface
    !             4 : fixed to zero on axis and fixed scale length at surface

    ! magnetic diffusion equation
    neq = 1
    id_neq(neq) = 2
!    id_neq(neq) = 0

    DO neq = 2, neqmax
       nsa = nsa_neq(neq)
       nva = nva_neq(neq)
       
       SELECT CASE(nva)
       CASE(1) ! density
          id_neq(neq) = 0
       CASE(2) ! toroidal velocity
          id_neq(neq) = 0
       CASE(3) ! temperature
          IF(idnsa(nsa) == 0 .OR. idnsa(nsa) == 2)THEN ! neutral and fast ions
             id_neq(neq) = 0
          ELSE ! bulk species
             id_neq(neq) = 1
          END IF
       END SELECT
    ENDDO

    RETURN
  END SUBROUTINE tr_set_idneq


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


  SUBROUTINE tr_print_species
! ------------------------------------------------------------------------
!   print conversion table for equations, species
! ------------------------------------------------------------------------
    USE trcomm,ONLY: neqmax,kidns,ns_nsa,nsa_neq,nva_neq,neqr_neq, &
         nsab_nsaf

    IMPLICIT NONE
    CHARACTER(1)   :: kid
    CHARACTER(50)  :: fmt_table
    INTEGER(ikind) :: neq,neqr,nsa,ns,nva,nsab

    fmt_table = '(1X,I3,I4,I5,I5,A6,I4)'
    WRITE(6,*) ! spacing
    WRITE(6,*) '# Variables conversion table'
    WRITE(6,*) 'NEQ NEQR NSA NSAB KIDNS NVA'

    neqr = 0
    nva  = 0
    nsa  = 0
    nsab = 0
    kid  = ' '
    
    DO neq = 1, neqmax
       neqr = neqr_neq(neq)
       IF(neq /= 1)THEN
          nsa  = nsa_neq(neq)
          nva  = nva_neq(neq)
          ns   = ns_nsa(nsa)
          nsab = nsab_nsaf(nsa)
          kid  = kidns(ns)
       END IF

       WRITE(6,fmt_table) neq,neqr,nsa,nsab,kid,nva
    END DO
    WRITE(6,*) '---------------------------'
    WRITE(6,*) ! spacing
    
    RETURN
  END SUBROUTINE tr_print_species

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_setup_profile
! -------------------------------------------------------------------------
!   This subroutine calculates inital profiles.
! -------------------------------------------------------------------------
    USE trcomm, ONLY: rkev,rkap,nrmax,nsmax,nsamax,ns_nsa,t, &
         rg,rm,rhog,rhom,rjcb,ra,rr,rn,ru,rt,rp,rp_tot,           &
         mdluf,mdleqn,mdlequ,mdleqt,mdleqm,mdlgmt,mdlglb,mdlsrc,time_slc
    USE trufile, ONLY: tr_ufile
    USE trufin, ONLY: tr_uf_init,tr_ufin_global,tr_ufin_density, &
                      tr_ufin_rotation,tr_ufin_temperature
    USE plprof, ONLY: pl_prof2
                      
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nsmax):: rn_ns,ru_ns,rtpr_ns,rtpp_ns
    REAL(rkind)    :: dr
    INTEGER(ikind) :: nr,nsa,ns,ierr

    ! generate mesh
    dr = SQRT(rkap)*ra/dble(nrmax)
    DO nr = 0, nrmax
       rg(nr)=dble(nr)*dr
       IF(nr /= 0) rm(nr)=0.5d0*(rg(nr-1)+rg(nr))

       rjcb(nr)    = 1.d0/(SQRT(rkap)*ra)
       rhog(nr)    = rg(nr)*rjcb(nr)
       IF(nr /= 0) rhom(nr) = 0.5d0*(rhog(nr-1)+rhog(nr))   
    END DO

    ! setup to read experimental data
    DO
       IF(mdluf /= 0)THEN
          CALL tr_ufile(ierr)
          IF(ierr /= 0)THEN
             ! set to default values
             mdluf  = 0
             mdlgmt = 1
             mdlglb = 1
             mdlsrc = 1
             EXIT
          END IF

          CALL tr_uf_init(mdluf)
          IF(mdluf==1)THEN
             time = time_slc
             t    = time
          ELSE IF(mdluf==2)THEN
             time = t
          END IF
          EXIT
       END IF
       EXIT
    END DO

    ! -------------------------------------------------------------------

    IF(mdluf == 0)THEN ! simple profile
       DO nr=0,nrmax
          CALL pl_prof2(rhog(nr),rn_ns,rtpr_ns,rtpp_ns,ru_ns)
          DO nsa=1,nsamax
             ns=ns_nsa(nsa)
             rn(nsa,nr)=rn_ns(ns)
             ru(nsa,nr)=ru_ns(ns)
             rt(nsa,nr)=(rtpr_ns(ns)+2.D0*rtpp_ns(ns))/3.D0          
          ENDDO
       ENDDO

    ELSE IF(mdluf > 0)THEN! experimental data

       SELECT CASE(mdlglb)
       CASE(6:7)
          CALL tr_ufin_global(time,0,ierr)
       END SELECT

       CALL tr_ufin_density(time,1,ierr)
       CALL tr_ufin_rotation(time,1,ierr)
       CALL tr_ufin_temperature(time,1,ierr)
    END IF


    ! pressure ----------------------------------------------
    rp_tot(0:nrmax) = 0.d0
    DO nr = 0, nrmax
       DO nsa = 1, nsamax
          rp(nsa,nr)  = rn(nsa,nr)*1.d20 * rt(nsa,nr)*rkev
          rp_tot(nr)  = rp_tot(nr) + rp(nsa,nr)          
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE tr_setup_profile


  SUBROUTINE tr_setup_geometry
! --------------------------------------------------------------------------
!
! --------------------------------------------------------------------------
    USE trcomm, ONLY: pi,nrmax,ra,rr,rkap,bb,rg,rm,dpdrho,           &
       ttrho,dvrho,abrho,abvrho,arrho,ar1rho,ar2rho,rmjrho,rmnrho,   &
       rkprho,rjcb,rhog,rhom,epsrho,abb2rho,abb1rho,pvolrho,         &
       aib2rho,psurrho,time_slc,mdluf,mdlgmt,knameq, &
       profj1,profj2,qp,jtot!, nrd1,nrd2

    USE trbpsd, ONLY: tr_bpsd_set,tr_bpsd_get
    USE trufin, ONLY: tr_ufin_geometry
    USE equnit_mod, ONLY: eq_parm,eq_prof,eq_calc,eq_load
!    USE equunit_mod, ONLY: equ_prof,equ_calc
!    USE pl_vmec_mod, ONLY: pl_vmec

    IMPLICIT NONE
    INTEGER(ikind)    :: nr,ierr, modelg
    REAL(rkind)       :: FCTR
    CHARACTER(len=80) :: line

    ! modelg --> mdlgmt


    ! Cylindrical geometry
    ! MDLGMT = 1  or  interim calculation for equilibrium code
    DO nr = 0, nrmax
       ar1rho(nr)  = 1.d0/(SQRT(rkap)*ra)          ! const
       ar2rho(nr)  = 1.d0/(SQRT(rkap)*ra)**2       ! const
       abrho(nr)   = 1.d0/(SQRT(rkap)*ra*rr)**2    ! const
       rmjrho(nr)  = rr                            ! const [m]
       rmnrho(nr)  = SQRT(rkap)*ra*rhog(nr) ! [m]
       rkprho(nr)  = rkap
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

    SELECT CASE(MDLGMT)
    CASE(2) ! Toroidal geometry

    CASE(3) ! TASK/EQ,EQDSK output geometry
       write(line,'(A,I5)') 'nrmax=',nrmax+1
       call eq_parm(2,line,ierr)
       write(line,'(A,I5)') 'nthmax=',64
       call eq_parm(2,line,ierr)
       write(line,'(A,I5)') 'nsumax=',0
       call eq_parm(2,line,ierr)

       modelg = mdlgmt
       call eq_load(modelg,knameq,ierr) ! load eq data and calculate eq
       IF(ierr.NE.0) THEN
          WRITE(6,*) 'XX eq_load: ierr=',ierr
          RETURN
       ENDIF
       call tr_bpsd_get(ierr)

    CASE(6:7) ! from experimental dataset
       CALL tr_ufin_geometry(time,1,ierr)

    CASE(8:9) ! CALL TASK/EQ       
       
       ! create interim q profile 
       !  using metric quantities on cylindrical assumption
       CALL tr_setup_field

       CALL tr_bpsd_set(ierr)
       CALL eq_calc              ! recalculate eq
       CALL tr_bpsd_get(ierr)
       ! --- here the convergence of q profile must be confirmed
       ! ---  and show the graph of Psi(R,Z)
    END SELECT

    RETURN
  END SUBROUTINE tr_setup_geometry


  SUBROUTINE tr_setup_field
! -----------------------------------------------------------------------
!
! -----------------------------------------------------------------------
    USE trcomm, ONLY: pi,nrmax,mdluf,mdlijq,ra,rkap,qp,rg,rhog, &
         profj1,profj2,jtot,joh,jcd_nb,jcd_ec,jcd_ic,jcd_lh,jbs_nc
         ! ,nrd1,nrd2
    USE trufin,ONLY: tr_ufin_field
    USE plprof,ONLY: pl_qprf

    IMPLICIT NONE
    INTEGER(ikind) :: nr,ierr


    IF(mdluf == 0)THEN ! simple profile
       IF(MOD(mdlijq,2)==1)THEN ! jtot --> dpdrho

          ! RIP is prior to jtot when not using exp. data
          IF(mdlijq == 3) mdlijq = 1

          DO nr = 0, nrmax
!             IF(((SQRT(rkap)*ra)**ABS(profj1)-rg(nr)**ABS(profj1)).LE.0.d0)THEN
             IF((1.d0-rhog(nr)**ABS(profj1))**ABS(profj2) < 0)THEN
                jtot(nr) = 0.D0
             ELSE
                jtot(nr) = &
!               ((SQRT(rkap)*ra)**ABS(profj1)-rg(nr)**ABS(profj1))**ABS(profj2)
                (1.d0-rhog(nr)**ABS(profj1))**ABS(profj2)
             ENDIF
          ENDDO

       ELSE IF(MOD(mdlijq,2)==0)THEN ! qp --> dpdrho

          ! RIP is prior to qp when not using exp. data
          IF(mdlijq == 4) mdlijq = 2

          ! simple q profile from q0, qa
          DO nr = 0, nrmax
             CALL pl_qprf(rhog(nr),qp(nr))
          END DO

       END IF

    ELSE IF(mdluf > 0)THEN
       CALL tr_ufin_field(time,1,ierr)

    END IF

    ! create d psi/d rho profile
    CALL tr_prof_dpdrho

    joh(0:nrmax) = jtot(0:nrmax)                        &
                   - (jcd_nb(0:nrmax) + jcd_ec(0:nrmax) &
                     +jcd_ic(0:nrmax)+jcd_lh(0:nrmax))  &
                   -  jbs_nc(0:nrmax)

    RETURN
  END SUBROUTINE tr_setup_field

! *************************************************************************
! *************************************************************************

  SUBROUTINE tr_prof_dpdrho
! -------------------------------------------------------------------------
!   toroidal current density profile --> d psi/d rho
! -------------------------------------------------------------------------
    USE trcomm,ONLY: rmu0,pi,nrmax,mdlijq,rhog,rip,abb1rho,dvrho,ttrho, &
         rmjrho,abrho,arrho,ar1rho,abvrho,dpdrho,rdpvrho,qp,bp,jtot
         ! ,nrd1,nrd2
    IMPLICIT NONE

    INTEGER(ikind):: nr
    REAL(rkind)   :: FCTR,DERIV3 ! the functions in TASK/lib
    REAL(rkind)   :: dr,dpdrhos,factor0,factor0p,factor0m,factorp,factorm,fact
    REAL(rkind),DIMENSION(0:nrmax) :: factor1,factor2


    IF(MOD(mdlijq,2)==1)THEN ! jtot --> dpdrho                          
       dpdrho(0:nrmax)  = 0.d0
       rdpvrho(0:nrmax) = 0.d0
       DO nr = 1, nrmax
          dr      = rhog(nr)-rhog(nr-1)
!          factor0 = rmu0*0.5d0*(abb1rho(nr)+abb1rho(nr-1)) &
!                        *0.5d0*(dvrho(nr)  +  dvrho(nr-1)) &
!                        *0.5d0*(jtot(nr)   +   jtot(nr-1)) &
!                       /(0.5d0*(ttrho(nr)  +  ttrho(nr-1)))**2
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
!          dpdrhos = 2.d0*pi*rmu0*rip*1.d6*dvrho(nrmax)/abvrho(nrmax)   
          ! correction in terms of the boundary value of dpdrho         
          fact = dpdrhos / dpdrho(nrmax)
!       write(*,*) fact                                              

          dpdrho(0:nrmax)  = fact*dpdrho(0:nrmax)
          rdpvrho(0:nrmax) = fact*rdpvrho(0:nrmax)
          jtot(0:nrmax)    = fact*jtot(0:nrmax)
       END IF

       ! dpdrho --> qp                                                 
       qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax)*dvrho(1:nrmax)    &
                    /(4.d0*pi**2 * dpdrho(1:nrmax))
!       qp(1:nrmax) = ttrho(1:nrmax)*arrho(1:nrmax) &
!                    /(4.d0*pi**2 * rdpvrho(1:nrmax))
       qp(0)       = FCTR(rhog(1),rhog(2),qp(1),qp(2))


    ELSE IF(MOD(mdlijq,2)==0)THEN ! qp --> dpdrho                       
       dpdrho(0:nrmax) = ttrho(0:nrmax)*arrho(0:nrmax)*dvrho(0:nrmax) &
                        / (4.d0*pi**2 * qp(0:nrmax))
       IF(mdlijq==2)THEN
          ! set the boundary value of dpdrho in terms of RIP            
          dpdrhos = 2.d0*pi*rmu0*rip*1.d6 / (dvrho(nrmax)*abrho(nrmax))
!          dpdrhos = 2.d0*pi*rmu0*rip*1.d6*dvrho(nrmax)/abvrho(nrmax)   
          ! correction in terms of the boundary value of dpdrho         
          fact = dpdrhos / dpdrho(nrmax)
!          write(*,*) fact                 
                             
          dpdrho(0:nrmax) = fact*dpdrho(0:nrmax)
          qp(0:nrmax)     = qp(0:nrmax) / fact
       END IF

       rdpvrho(0:nrmax) = ttrho(0:nrmax)*arrho(0:nrmax) &
                         /(4.d0*pi**2*qp(0:nrmax))! d psi/d V 

       ! dpdrho --> jtot
       factor1(0:nrmax) = dvrho(0:nrmax)*abrho(0:nrmax)*dpdrho(0:nrmax)
       factor2(0:nrmax) = factor1(0:nrmax)/ttrho(0:nrmax)
       DO nr = 1, nrmax
          ! dpdrho --> jtot(j_para)
          jtot(nr) = ttrho(nr)**2/(rmu0*abb1rho(nr)*dvrho(nr)) &
                    * deriv3(nr,rhog,factor2,nrmax,0)
       END DO
!       jtot(0) = FCTR4pt(rhog(1),rhog(2),rhog(3),jtot(1),jtot(2),jtot(3))
       jtot(0) = FCTR(rhog(1),rhog(2),jtot(1),jtot(2))

    END IF

    ! dpdrho --> bp                                               
    bp(0:nrmax) = ar1rho(0:nrmax)*dpdrho(0:nrmax)/rmjrho(0:nrmax)

    RETURN
  END SUBROUTINE tr_prof_dpdrho


  SUBROUTINE tr_setup_varinit
! --------------------------------------------------------------------
!   Initialization for successive interactive calculation
! --------------------------------------------------------------------
    USE trcomm,ONLY: nrmax, &
         nrd1,nrd2,nrd3,nrd4,vtor,vpol,vpar,vprp,  &
         jbs_nc,jex_nc,jcd_nb,jcd_ec,jcd_ic,jcd_lh

    nrd1(0:nrmax) = 0.d0
    nrd2(0:nrmax) = 0.d0
    nrd3(0:nrmax) = 0.d0
    nrd4(0:nrmax) = 0.d0

    vtor(0:nrmax) = 0.d0
    vpol(0:nrmax) = 0.d0
    vpar(0:nrmax) = 0.d0
    vprp(0:nrmax) = 0.d0

    jbs_nc(0:nrmax) = 0.d0
    jex_nc(0:nrmax) = 0.d0
    jcd_nb(0:nrmax) = 0.d0
    jcd_ec(0:nrmax) = 0.d0
    jcd_ic(0:nrmax) = 0.d0
    jcd_lh(0:nrmax) = 0.d0

    RETURN
  END SUBROUTINE tr_setup_varinit

END MODULE trsetup
