MODULE trcalc2

  USE trcomm,ONLY: ikind,rkind

  PRIVATE
  PUBLIC tr_calc2

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_calc2
    USE trcomm, ONLY: rkev,nrmax,nsamax,neqmax,nsa_neq,nva_neq,   &
         dtr,vtr,ctr,str,htr,dtr_nc,vtr_nc,dtr_tb,vtr_tb,ctr_ex,  &
         htr_simple,jtot,joh,jcd_nb,jcd_ec,jcd_ic,jcd_lh,jbs_nc
    USE trcalv, ONLY: tr_calc_variables
    USE trcoeftb, ONLY: tr_coeftb
    USE trcoefnc, ONLY: tr_coefnc
    USE trsource, ONLY: tr_source2
    IMPLICIT NONE

    INTEGER(ikind):: nr,neq,nsa,nva

    dtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    vtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    ctr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    htr(1:neqmax,0:nrmax)=0.D0 ! external driven current


    CALL tr_calc_variables

    ! *** dtr, vtr ***
    ! *** CAUTION: Physical quantities associated with neoclassical theory
    !               must be calculated before calculating 'htr' and 'str'.

    CALL tr_coefnc         ! calculate neoclassical transport coefficients
    CALL tr_coeftb         ! calculate turbulent transport coefficients
    CALL tr_coefmg


    ! *** ctr ***
    CALL tr_calc2_energy_ex


    ! *** htr *** this section should be in the trsource directory ??
    CALL tr_calc2_excurrent ! calculate external driven current term

    joh(0:nrmax) = jtot(0:nrmax) - jbs_nc(0:nrmax)     &
              -( jcd_nb(0:nrmax) + jcd_ec(0:nrmax)     &
               + jcd_ic(0:nrmax) + jcd_lh(0:nrmax))

    ! ***    

    ! *** str ***
    CALL tr_source2


    ! substitution
    DO nr=1,nrmax
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


    RETURN
  END SUBROUTINE tr_calc2

! ***********************************************************************
! ***********************************************************************

  SUBROUTINE tr_coefmg
! -----------------------------------------------------------------------
!   calculate coefficients for poloidal magnetic diffusion equation
! -----------------------------------------------------------------------
    USE trcomm, ONLY: rmu0,nrmax,mdltr_nc, &
                      dvrho,ttrho,arrho,abb1rho,eta,etam_nc,dtr,htr
    IMPLICIT NONE
    INTEGER(ikind) :: nr
    REAL(rkind) :: etam,ttrhom,arrhom,dvrhom

    ! registivity term (half grid)
    DO nr = 1, nrmax
       etam   = 0.5d0*(eta(nr)+eta(nr-1))
       ttrhom = 0.5d0*(ttrho(nr)+ttrho(nr-1))
       arrhom = 0.5d0*(arrho(nr)+arrho(nr-1))
       dvrhom = 0.5d0*(dvrho(nr)+dvrho(nr-1))
       
       dtr(1,1,nr) = etam*ttrhom/(rmu0*dvrhom*arrhom)
    END DO

  END SUBROUTINE tr_coefmg


  SUBROUTINE tr_calc2_excurrent
! ----------------------------------------------------------------------
!   calculate external driven current
! ----------------------------------------------------------------------
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
  END SUBROUTINE tr_calc2_excurrent


  SUBROUTINE tr_calc2_energy_ex
! ----------------------------------------------------------------------
!   calculate energy exchange rate
! ----------------------------------------------------------------------
    USE trcomm, ONLY: aee,ame,amp,pi,rkev,pa,pz,eps0,nrmax,neqmax,id_neq, &
         ns_nsa,nva_neq,nsa_neq,rt,rn,ctr_ex
    USE trcalv, ONLY: coulog ! coulomb logarithm
    IMPLICIT NONE
    INTEGER(ikind) :: nr,neq,neq1,nsa,nsa1,ns,ns1
    REAL(rkind) :: coef1,coef2,ams,ams1

    ctr_ex(1:neqmax,1:neqmax,0:nrmax) = 0.D0

    ! *** heat exchange                          ***
    ! *** relaxation process (Tokamaks 3rd p.68) ***
    coef1 = aee**4.d0 / (3.d0*SQRT(2.d0)*pi**1.5d0 * eps0**2.d0)

    DO nr = 0, nrmax
       DO neq = 1, neqmax
          IF(id_neq(neq)==0) CYCLE
          IF(nva_neq(neq) == 3)THEN ! only for energy equation
             DO neq1 = 1, neqmax
                IF(id_neq(neq1)==0) CYCLE
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
                   *(rt(nsa,nr)*rkev/ams + rt(nsa1,nr)*rkev/ams1)**(-1.5d0) &
                   *coulog(ns,ns1,rn(1,nr),rt(nsa,nr))

                   ! coef1*coef2*n_j : nu_ij (heat exchange frequency)
                   ctr_ex(neq,neq,nr)  = ctr_ex(neq,neq,nr) &
                                       - 1.5d0*coef1*coef2*rn(nsa1,nr)*1.d20
                   ctr_ex(neq,neq1,nr) = 1.5d0*coef1*coef2*rn(nsa, nr)*1.d20
!                   write(6,*) 1.d0/(coef1*coef2*rn(nsa1,nr))
                END IF
             END DO
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_calc2_energy_ex

END MODULE trcalc2
