MODULE trcalc

  USE trcomm,ONLY: ikind,rkind

  PRIVATE
  PUBLIC tr_calc, tr_calc_dpdrho2j

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_calc
    USE trcomm, ONLY: nrmax,nsamax,neqmax,rkev,dtr,vtr,ctr,str,htr, &
         nsa_neq,nva_neq, &
         dtr_nc,vtr_nc,dtr_tb,vtr_tb,ctr_ex,str_simple,htr_simple,  &
         jbs_nc,jex_nc,eta,poh,nrd4
    USE trcalv, ONLY: tr_calv_nr_alloc,tr_calc_variables
    USE trcoeftb, ONLY: tr_coeftb
    USE trcoefnc, ONLY: tr_coefnc
    IMPLICIT NONE
    INTEGER(ikind):: nr,neq,nsa,nva

    CALL tr_calc_dpdrho2j

    CALL tr_calv_nr_alloc
    CALL tr_calc_variables

    CALL tr_coefnc         ! calculate neoclassical transport coefficients
    CALL tr_coeftb         ! calculate turbulent transport coefficients

    CALL tr_calc_exchange  ! calculate exchange rate
    CALL tr_calc_excurrent ! calculate external driven current term
    CALL tr_calc_source    ! calculate source and sink term

    dtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    vtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    ctr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    htr(1:neqmax,0:nrmax)=0.D0 ! external driven current
    str(1:neqmax,0:nrmax)=0.D0

    DO nr=1,nrmax

       call tr_calc_mag_diff
!       nrd4(0:nrmax) = dtr(1,1,0:nrmax)

       DO neq=2,neqmax
          dtr(2:neqmax,neq,nr) &
               =dtr_nc(2:neqmax,neq,nr) &
               +dtr_tb(2:neqmax,neq,nr)
          vtr(2:neqmax,neq,nr) &
               =vtr_nc(2:neqmax,neq,nr) &
               +vtr_tb(2:neqmax,neq,nr)
          ctr(2:neqmax,neq,nr) &
               =ctr_ex(2:neqmax,neq,nr)
       END DO
    END DO

!    str(2:neqmax,0:nrmax) = str_simple(2:neqmax,0:nrmax)
    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       nva = nva_neq(neq)
       IF(nva==3 .AND. nsa==1)THEN ! energy, electron
          str(neq,0:nrmax) = poh(0:nrmax)/(rkev*1.d20)
       END IF
    END DO


    ! bootstrap current
    htr(1,0:nrmax) = jbs_nc(0:nrmax) + jex_nc(0:nrmax)

!    htr(1,0:nrmax) = htr_simple(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc

! ----- calculate exchange rate -----

  SUBROUTINE tr_calc_exchange
    USE trcomm, ONLY: aee,ame,amp,pi,rkev,pa,pz,eps0,nrmax,neqmax, &
         ns_nsa,nva_neq,nsa_neq,rt,rn,ctr_ex
    USE trcoefnc, ONLY: coulog ! coulomb logarithm
    IMPLICIT NONE
    INTEGER(ikind) :: nr,neq,neq1,nsa,nsa1,ns,ns1
    REAL(rkind) :: coef1,coef2,ams,ams1

    ctr_ex(1:neqmax,1:neqmax,0:nrmax)=0.D0

    coef1 = aee**4.d0 / (3.d0*SQRT(2.d0)*pi**1.5d0 * eps0**2.d0)

    ! *** only heat exchange ***
    DO nr = 0, nrmax
       DO neq = 1, neqmax
          IF(nva_neq(neq) == 3)THEN ! only for energy equation
             DO neq1 = 1, neqmax
                IF(nva_neq(neq1) == 3 .AND. neq /= neq1)THEN
                   nsa  = nsa_neq(neq)
                   nsa1 = nsa_neq(neq1)
                   ns  = ns_nsa(nsa)
                   ns1 = ns_nsa(nsa1)

                   ! particle mass
                   IF(pz(ns) < 0.d0)THEN ! electron
                      ams = ame
                   ELSE IF(pz(ns) > 0.d0)THEN ! ions
                      ams = amp*pa(ns)
                   END IF
                   IF(pz(ns1) < 0.d0)THEN ! electron
                      ams1 = ame
                   ELSE IF(pz(ns1) > 0.d0)THEN ! ions
                      ams1 = amp*pa(ns1)
                   END IF

                   coef2 = (pz(ns)*pz(ns1))**2.d0 / (ams*ams1)             &
                   *(rt(nsa,nr)*rkev/ams + rt(nsa1,nr)*rkev/ams1)**(-1.5d0)&
                   *coulog(ns,ns1,rn(1,nr),rt(nsa,nr))                     &
                   *rn(nsa1,nr)*1.d20

                   ! coef1*coef2 : nu_ij (heat exchange frequency)
                   ctr_ex(neq,neq1,nr) =   1.5d0*rn(nsa,nr)*coef1*coef2 &
                                           /rn(nsa1,nr)
                   ctr_ex(neq,neq ,nr) = - 1.5d0*rn(nsa,nr)*coef1*coef2 &
                                           /rn(nsa ,nr)
                END IF
             END DO
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_calc_exchange

! ----- calculate external driven current -----

  SUBROUTINE tr_calc_excurrent
    USE trcomm, ONLY: nrmax,neqmax,nva_neq,rhog,htr_simple
    IMPLICIT NONE
    INTEGER(ikind) :: nr,neq

    htr_simple(0:nrmax)=0.D0
    DO nr = 0, nrmax
       DO neq = 1, neqmax
          IF(nva_neq(neq) == 0) THEN
             htr_simple(nr) &
                  = 5.d5 * (1.d0-rhog(nr)**2.d0)**1.5d0
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_calc_excurrent

! ----- calculate source -----

  SUBROUTINE tr_calc_source
    USE trcomm, ONLY: nrmax,nsamax,neqmax,nva_neq,ph0,phs,rhog,ra, &
         str_simple,joh,eta,ezoh,poh,nrd1,nrd2,nrd3,nrd4
    IMPLICIT NONE
    INTEGER(ikind) :: nr, neq

    str_simple = 0.d0
    DO nr = 0, nrmax
       DO neq=1,neqmax
          IF(nva_neq(neq) == 3) THEN
             str_simple(neq,nr) = phs+(ph0-phs)*(1.D0-(rhog(nr)/ra)**2)
!             str_simple(neq,nr) = phs+(ph0-phs)*(1.D0-(rg(nr)/ra)**2)
!             str_simple(neq,nr) = 0.d0
          END IF
       END DO
    END DO

    ! ohmic heating [W/m^3]
    DO nr = 0, nrmax
       ezoh(nr) = eta(nr)*joh(nr)
       poh(nr) = ezoh(nr)*joh(nr)    
    END DO

    nrd1(0:nrmax) = eta(0:nrmax)
    nrd2(0:nrmax) = ezoh(0:nrmax)
    nrd3(0:nrmax) = poh(0:nrmax)
    nrd4(0:nrmax) = joh(0:nrmax)

    RETURN
  END SUBROUTINE tr_calc_source

! --- calculate coefficients for poloidal magnetic diffusion equation ---

  SUBROUTINE tr_calc_mag_diff
    USE trcomm, ONLY: rmu0,nrmax,dvrho,ttrho,arrho,abb1rho,eta,dtr,htr
    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: etam,ttrhom,arrhom,dvrhom

    ! registivity term (half grid)
    DO nr = 1, nrmax
       etam   = 0.5d0*(eta(nr)+eta(nr-1))
       ttrhom = 0.5d0*(ttrho(nr)+ttrho(nr-1))
       arrhom = 0.5d0*(arrho(nr)+arrho(nr-1))
       dvrhom = 0.5d0*(dvrho(nr)+dvrho(nr-1))

       dtr(1,1,nr) = etam/rmu0 * ttrhom*arrhom/dvrhom
    END DO
!    dtr(1,1,1:nrmax) = 2.d0

  END SUBROUTINE tr_calc_mag_diff

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

    REAL(rkind) :: dr,factorp,factorm,factor0,fact,dpdrhos

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

! -----------------------------------------------------------------------
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
    jtor(0:nrmax) = fact*jtor(0:nrmax)
!-------------------------------------------------------------------------

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

END MODULE trcalc
