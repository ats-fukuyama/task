MODULE trsetup

! This module setup table for computation and 
! initializes the plasma profiles

  PUBLIC tr_setup
  PRIVATE

CONTAINS

  SUBROUTINE tr_setup

    USE trcomm, ONLY: ikind,t,ngt,kidnsa,ns_nsa,idnsa,nsamax,pa,pz,pz0, &
         tr_nit_allocate,tr_nsa_allocate,tr_nr_allocate,tr_ngt_allocate,&
         nitmax
    USE trbpsd, ONLY: tr_bpsd_init
    USE trloop, ONLY: tr_save_pvprev
    USE trresult, ONLY: tr_calc_global,tr_save_ngt
    IMPLICIT NONE
    INTEGER(ikind):: nsa,ns,ierr

    CALL tr_nsa_allocate

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

    CALL tr_ngt_allocate

    CALL tr_setup_table     ! calculate table for matrix generation

    CALL tr_nr_allocate
    CALL tr_setup_profile   ! calculate initial profile

    CALL tr_nit_allocate

    t      = 0.D0           ! time is initialized
    ngt    = -1             ! save data count is initilaized
    nitmax = 0              ! iteration count is initialized

    CALL tr_calc_global
    CALL tr_save_ngt

    CALL tr_save_pvprev

    CALL tr_bpsd_init(ierr)

    ! read TASK/EQ output file at the initilization phase.
!    CALL TASK/EQ
    CALL tr_setup_geometric ! calculate geometric factor in TASK/TR

    RETURN
  END SUBROUTINE tr_setup

! ***** calculate table for matrix generation *****

  SUBROUTINE tr_setup_table

    USE trcomm, ONLY: ikind,nrmax,nsamax,neqmax,neqrmax,nvrmax,nvmax, &
                      nsa_neq,nva_neq,id_neq,id_neqnr,neq_neqr, &
                      rg,rg_fixed,neqr_neq,tr_neq_allocate,tr_neqr_allocate
    IMPLICIT NONE
    INTEGER(ikind):: neq,nsa,nva,i,neqr,nr

    neqmax=1+3*nsamax
    nvmax=neqmax*(nrmax+1)

    CALL tr_neq_allocate

!   nsa_neq = 0 : magnetic field
!             otherwise : particle species
!   nva_neq = 1 : density
!             2 : toroidal velocity
!             3 : temperature
!   id_neq  = 0 : not solved
!             1 : solved

    neq=1
    nsa_neq(neq)=0
    nva_neq(neq)=1
    id_neq(neq)=0

    DO nsa=1,nsamax
       DO i=1,3
          neq=neq+1
          nsa_neq(neq)=nsa
          nva_neq(neq)=i
          SELECT CASE(i)
          CASE(1)
             id_neq(neq)=0
          CASE(2)
             id_neq(neq)=0
          CASE(3)
             id_neq(neq)=1
          END SELECT
       ENDDO
    ENDDO

    neqr=0
    DO neq=1,neqmax
       SELECT CASE(id_neq(neq))
!         id_neqnr = 0 : fixed to zero
!                    1 : solved normally
!                    2 : fixed to a given value
!                    3 : fixed to a given scale length
       CASE(0) ! equation is not solved in any radius (fixed to zero)
          neqr_neq(neq)=0
          id_neqnr(neq,0:nrmax)=0
       CASE(1,11) ! flat on axis and fixed at plasma surface
          neqr=neqr+1
          neqr_neq(neq)=neqr
          id_neqnr(neq,0:nrmax-1)=1
          id_neqnr(neq,nrmax)=2
       CASE(2,12) ! fixed to zero on axis and fixed at plasma surface
          neqr=neqr+1
          neqr_neq(neq)=neqr
          id_neqnr(neq,0)=2
          id_neqnr(neq,1:nrmax-1)=1
          id_neqnr(neq,nrmax)=2
       CASE(3,13) ! flat on axis and fixed scale length at plasma surface
          neqr=neqr+1
          neqr_neq(neq)=neqr
          id_neqnr(neq,0:nrmax-1)=1
          id_neqnr(neq,nrmax)=3
       CASE(4,14) ! fixed to zero on axis and fixed scale length at surface
          neqr=neqr+1
          neqr_neq(neq)=neqr
          id_neqnr(neq,0)=2
          id_neqnr(neq,1:nrmax-1)=1
          id_neqnr(neq,nrmax)=3
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
    neqrmax=neqr
    nvrmax=neqrmax*(nrmax+1)

    CALL tr_neqr_allocate

    DO neq=1,neqmax
       IF(neqr_neq(neq) /= 0) neq_neqr(neqr_neq(neq))=neq
    ENDDO
          

    RETURN
  END SUBROUTINE tr_setup_table

! ***** calculate inital profile *****

  SUBROUTINE tr_setup_profile

    USE trcomm, ONLY: ikind,rkind,rkap,rdlt,nrmax,nsmax,nsamax, &
         rg,rm,ra,rn,ru,rt,ns_nsa,qp
    USE plprof, ONLY: pl_prof2,pl_qprf
                      
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nsmax):: rn_ns,ru_ns,rtpr_ns,rtpp_ns
    INTEGER(ikind):: nr,nsa,ns
    REAL(rkind):: dr,rhon

    dr=SQRT(rkap)*ra/dble(nrmax)
    DO nr=0,nrmax
       rg(nr)=dble(nr)*dr
       IF(nr /= 0) rm(nr)=0.5d0*(rg(nr-1)+rg(nr))
       rhon=rg(nr)/ra
       CALL pl_qprf(rhon,qp(nr))
       CALL pl_prof2(rhon,rn_ns,rtpr_ns,rtpp_ns,ru_ns)
       DO nsa=1,nsamax
          ns=ns_nsa(nsa)
          rn(nsa,nr)=rn_ns(ns)
          ru(nsa,nr)=ru_ns(ns)
          rt(nsa,nr)=(rtpr_ns(ns)+2.D0*rtpp_ns(ns))/3.D0

          ! MDLUF = 0 : trprof.f90
!          pex(nsa,nr) = 0.d0
!          sex(nsa,nr) = 0.d0
!          prf(nsa,nr) = 0.d0

!          pbm(nr)     = 0.d0
!          wrot(nr)    = 0.d0
!          vtor(nr)    = 0.d0
       ENDDO
    ENDDO
  END SUBROUTINE tr_setup_profile


  SUBROUTINE tr_setup_geometric

    USE trcomm, ONLY: ikind,rkind,pi,nrmax,ra,rr,rkap,bb,rg,qp,rhoa,  &
       ttrho,dvrho,abrho,abvrho,arrho,ar1rho,ar2rho,rmjrho,rmnrho,rkprho, &
       rjcb,rhog,epsrho,abb2rho,pvolrho,psurrho,rdpvrho,bp,rdp

    IMPLICIT NONE
    INTEGER(ikind) :: nr

    rhoa = 1.D0

    ! --- values on grid mesh ---
    DO nr = 0, nrmax

!       IF(MDLUF .ne. 0) THEN ! input UFILE
!          DO nr = 0, nrmax
!             ! normalized variables
!             rjcb(nr)    = 1.d0/(SQRT(rkap)*ra) ! const
!             rhog(nr)    = rg(nr)/rjcb(nr)
!             !       rhom(nr)    = rm(nr)/rjcb(nr)
!          END DO

!       ELSE IF ! not input UFILE
             ttrho(nr)   = bb*rr ! const

             pvolrho(nr) = pi*rkap*(ra*rg(nr))**2*2.d0*pi*rr
             psurrho(nr) = pi*(rkap+1.d0)*ra*rg(nr)*2.d0*pi*rr
             dvrho(nr)   = 2.d0*pi*rkap*ra**2*2.d0*pi*rr*rg(nr)
!
             arrho(nr)   = 1.d0/rr**2 ! const

             ar1rho(nr)  = 1.d0/(SQRT(rkap)*ra) ! const
             ar2rho(nr)  = 1.d0/(SQRT(rkap)*ra)**2 ! const
             abrho(nr)   = 1.d0/(SQRT(rkap)*ra*rr)**2   ! const
             rmjrho(nr)  = rr ! const
             rmnrho(nr)  = ra*rg(nr) 
             rkprho(nr)  = rkap
!
             epsrho(nr)  = rmnrho(nr)/rmjrho(nr)

             abb2rho(nr) = bb*(1.d0+0.25d0*epsrho(nr)**2)
             abvrho(nr)  = dvrho(nr)**2*abrho(nr)

             ! normalized variables
             rjcb(nr)    = 1.d0/(SQRT(rkap)*ra)
             rhog(nr)    = rg(nr)/rjcb(nr)
!             rhom(nr)    = 

!       END IF

! create BP from given Q profile
       ! d psi/d V
       rdpvrho(nr) = ttrho(nr)*arrho(nr)/(4.d0*pi**2*qp(nr))
       rdp(nr)     = dvrho(nr)*rdpvrho(nr) ! d psi/d rho

       ! This part should be calculated in another module ----------
       ! poloidal magnetic field ~ (kappa^-2 * r * BB)/(RR * q)
       bp(nr)      = ar1rho(nr)*rdp(nr)/rr ! poloidal magnetic field
       ! -----------------------------------------------------------

    END DO

  END SUBROUTINE tr_setup_geometric
END MODULE trsetup
