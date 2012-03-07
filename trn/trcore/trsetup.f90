MODULE trsetup

! This module setup table for computation and 
! initializes the plasma profiles
  USE trcomm, ONLY: ikind, rkind

  PUBLIC tr_setup
  PRIVATE

CONTAINS

  SUBROUTINE tr_setup

    USE trcomm, ONLY: t,ngt,kidnsa,ns_nsa,idnsa,nsamax,pa,pz,pz0,        &
         tr_nit_allocate,tr_nsa_allocate,tr_nr_allocate,tr_ngt_allocate, &
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

    CALL tr_ngt_allocate    ! allocation for data save

    CALL tr_setup_table     ! calculate table for matrix generation


    CALL tr_nr_allocate     ! allocation for radial profile

    ! read TASK/EQ output file at the initilization phase.
!    IF(
!    CALL TASK/EQ
    CALL tr_setup_geometric ! calculate geometric factor in TASK/TR
    CALL tr_setup_profile   ! calculate initial profile

    CALL tr_setup_htr_prof

    CALL tr_nit_allocate    ! allocation for diagnostics of iteration

    t      = 0.D0           ! time is initialized
    ngt    = -1             ! save data count is initilaized
    nitmax = 0              ! iteration count is initialized

    CALL tr_calc_global
    CALL tr_save_ngt

    CALL tr_save_pvprev

    CALL tr_bpsd_init(ierr)

    RETURN
  END SUBROUTINE tr_setup

! ***** calculate table for matrix generation *****

  SUBROUTINE tr_setup_table

    USE trcomm, ONLY: nrmax,nsamax,neqmax,neqrmax,nvrmax,nvmax, &
                      nsa_neq,nva_neq,id_neq,id_neqnr,neq_neqr,       &
                      rg,rg_fixed,neqr_neq,tr_neqr_allocate
    IMPLICIT NONE
    INTEGER(ikind):: neq,nsa,nva,neqr,nr

    neqr = 0
    DO neq = 1, neqmax
       SELECT CASE(id_neq(neq))
!         id_neqnr = 0 : fixed to zero
!                    1 : solved normally
!                    2 : fixed to a given value
!                    3 : fixed to a given scale length
       CASE(0) ! equation is not solved in any radius (fixed to zero)
          neqr_neq(neq) = 0
          id_neqnr(neq,0:nrmax) = 0
       CASE(1,11) ! flat on axis and fixed at plasma surface
          neqr = neqr+1
          neqr_neq(neq) = neqr
          id_neqnr(neq,0:nrmax-1) = 1
          id_neqnr(neq,    nrmax) = 2
       CASE(2,12) ! fixed to zero on axis and fixed at plasma surface
          neqr = neqr+1
          neqr_neq(neq) = neqr
          id_neqnr(neq,        0) = 2
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
          id_neqnr(neq,        0) = 2
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
    write(*,*) neqrmax

    CALL tr_neqr_allocate

    DO neq=1,neqmax
       IF(neqr_neq(neq) /= 0) neq_neqr(neqr_neq(neq))=neq
    ENDDO
          

    RETURN
  END SUBROUTINE tr_setup_table

! ***************************************************************************
  SUBROUTINE tr_setup_geometric

    USE trcomm, ONLY: pi,nrmax,ra,rr,rkap,bb,rg,rm,rhoa,                    &
       ttrho,dvrho,abrho,abvrho,arrho,ar1rho,ar2rho,rmjrho,rmnrho,rmnrhom,  &
       rkprho,rkprhom,rjcb,rhog,rhom,epsrho,abb2rho,abb1rho,pvolrho,psurrho,&
       modelg!, nrd1,nrd2

    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: dr

    rhoa = 1.D0
    dr   = SQRT(rkap)*ra/dble(nrmax)

    ! --- values on grid mesh ---
    DO nr = 0, nrmax
       rg(nr)=dble(nr)*dr
       IF(nr /= 0) rm(nr)=0.5d0*(rg(nr-1)+rg(nr))

       SELECT CASE(modelg)
          CASE(2) ! Toloidal geometry : cylindrical assumption
             ! normalized variables
             rjcb(nr)    = 1.d0/(SQRT(rkap)*ra)
             rhog(nr)    = rg(nr)*rjcb(nr)
             IF(nr /= 0) rhom(nr) = 0.5d0*(rhog(nr-1)+rhog(nr))
             
             pvolrho(nr) = pi*rkap*(ra*rhog(nr))**2*2.d0*pi*rr
             psurrho(nr) = pi*(rkap+1.d0)*ra*rhog(nr)*2.d0*pi*rr
             dvrho(nr)   = 2.d0*pi*rkap*ra**2*2.d0*pi*rr*rhog(nr)
!
             ar1rho(nr)  = 1.d0/(SQRT(rkap)*ra)          ! const
             ar2rho(nr)  = 1.d0/(SQRT(rkap)*ra)**2       ! const
             abrho(nr)   = 1.d0/(SQRT(rkap)*ra*rr)**2    ! const
             rmjrho(nr)  = rr                            ! const [m]
             rmnrho(nr)  = ra*rhog(nr) ! [m]
             IF(nr /= 0) rmnrhom(nr)=0.5d0*(rmnrho(nr-1)+rmnrho(nr))
             rkprho(nr)  = rkap
             IF(nr /= 0) rkprhom(nr)=0.5d0*(rkprho(nr-1)+rkprho(nr))
!
             epsrho(nr)  = rmnrho(nr)/rmjrho(nr)

             abb1rho(nr) = BB*(1.d0 + 0.5d0*epsrho(nr)**2) ! <B>
             ttrho(nr)   = abb1rho(nr) * rr
!             arrho(nr)   = 1.d0/rr**2                    ! const
             arrho(nr)   = 1.d0/rr**2 * (1+1.5d0*epsrho(nr)**2)
            
!             abb2rho(nr) = 
             abvrho(nr)  = dvrho(nr)**2*abrho(nr)

          CASE(3) ! TASK/EQ output geometry
             continue
          CASE(8) ! CALL TASK/EQ
             continue
       END SELECT
!       nrd1(0:nrmax) = abb1rho(0:nrmax)
!       nrd2(0:nrmax) = ttrho(0:nrmax)

    END DO

  END SUBROUTINE tr_setup_geometric


! ***************************************************************************
! ***** calculate inital profile *****
  SUBROUTINE tr_setup_profile

    USE trcomm, ONLY: ikind,rkind,pi,rkap,rdlt,nrmax,nsmax,nsamax, &
         rg,rm,rhog,rhom,ra,rr,rn,ru,rt,ns_nsa,qp, &
         ttrho,dvrho,arrho,ar1rho,rdpvrho,dpdrho,bp, &
         mdluf
    USE plprof, ONLY: pl_prof2,pl_qprf
                      
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nsmax):: rn_ns,ru_ns,rtpr_ns,rtpp_ns
    INTEGER(ikind):: nr,nsa,ns

    IF(mdluf == 0)THEN
       DO nr=0,nrmax
          CALL pl_qprf(rhog(nr),qp(nr))
          CALL pl_prof2(rhog(nr),rn_ns,rtpr_ns,rtpp_ns,ru_ns)
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
    END IF
  END SUBROUTINE tr_setup_profile

! ***************************************************************************
  SUBROUTINE tr_setup_htr_prof
    USE trcomm, ONLY: rmu0,nrmax,BB,RR,abb1rho,dvrho,ttrho,abrho,dpdrho, &
         ar1rho,rhog,jtot,joh,bp, nrd1,nrd2

    IMPLICIT NONE
    REAL(rkind) :: dr,prof,factor0,sumfact1,profj1,profj2
    INTEGER(ikind) :: nr


    ! create BP from given Q profile
!    rdpvrho(nr) = ttrho(nr)*arrho(nr)/(4.d0*pi**2*qp(nr))! d psi/d V
!    dpdrho(nr)  = dvrho(nr)*rdpvrho(nr)                  ! d psi/d rho
    
    profj1 = 2.d0
    profj2 = 1.d0

    DO nr = 0, nrmax
       IF((1.d0-rhog(nr)**ABS(profj1)).LE.0.d0) THEN
          prof = 0.D0
       ELSE
          prof = 1.d6 * (1.D0-rhog(nr)**ABS(profj1))**ABS(profj2)
       ENDIF
       joh(nr)  = prof
       jtot(nr) = prof
    ENDDO

    sumfact1 = 0.d0
    dpdrho   = 0.d0
    DO nr = 1, nrmax
       dr = rhog(nr)-rhog(nr-1)
       factor0 = ttrho(nr)/(dvrho(nr)*abrho(nr))
       sumfact1 = sumfact1 +                               &
                  rmu0*abb1rho(nr)*dvrho(nr)/ttrho(nr)**2  &
                   * 0.5d0*(jtot(nr)+jtot(nr-1)) * dr
       dpdrho(nr) = factor0*sumfact1

       ! poloidal magnetic field ~ (kappa^-2 * r * BB)/(RR * q)
       bp(nr)      = ar1rho(nr)*dpdrho(nr)/rr

    END DO

    nrd1(0:nrmax) = jtot(0:nrmax)
    nrd2(0:nrmax) = dpdrho(0:nrmax)


!!$    NR=1    
!!$    FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
!!$    RDPVRHOG(NR)=FACTOR0*DR/FACTORP
!!$    RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
!!$    BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
!!$    DO NR=2,NRMAX
!!$       FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
!!$       FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
!!$       FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
!!$       RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
!!$       RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
!!$       BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
!!$    ENDDO
!!$    NR=1
!!$    FACTOR0=RR/(RMU0*DVRHO(NR))
!!$    FACTORP=ABVRHOG(NR  )
!!$    AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
!!$    DO NR=2,NRMAX
!!$       FACTOR0=RR/(RMU0*DVRHO(NR))
!!$       FACTORM=ABVRHOG(NR-1)
!!$       FACTORP=ABVRHOG(NR  )
!!$       AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/D\
!!$       R
!!$    ENDDO
!!$    
!!$    RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
!!$    FACT=RDPS/RDP(NRMAX)
!!$    RDP(1:NRMAX)=FACT*RDP(1:NRMAX)
!!$    RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
!!$    AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
!!$    AJ(1:NRMAX)  =AJOH(1:NRMAX)
!!$    BP(1:NRMAX)  =FACT*BP(1:NRMAX)
!!$    QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:\
!!$    NRMAX))
    

    RETURN
  END SUBROUTINE tr_setup_htr_prof

END MODULE trsetup
