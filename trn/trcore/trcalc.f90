MODULE trcalc

  USE trcomm,ONLY: ikind,rkind

  PRIVATE
  PUBLIC tr_calc,tr_calc_source

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_calc
    USE trcomm, ONLY: nrmax,nsamax,neqmax,rkev,dtr,vtr,ctr,str,htr, &
         nsa_neq,nva_neq, &
         dtr_nc,vtr_nc,dtr_tb,vtr_tb,ctr_ex,str_simple,htr_simple,  &
         jbs_nc,jex_nc,eta,poh,pnb,nrd4
    USE trcalv, ONLY: tr_calc_variables
    USE trcoeftb, ONLY: tr_coeftb
    USE trcoefnc, ONLY: tr_coefnc
    IMPLICIT NONE
    INTEGER(ikind):: nr,neq,nsa,nva

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

!    str(2:neqmax,0:nrmax) = str_simple(2:neqmax,0:nrmax)*1.d6
    DO neq = 1, neqmax
       nsa = nsa_neq(neq)
       nva = nva_neq(neq)
       !  (NBI) power is edeposited to electron and ion equally
       IF(nva==3 .AND. nsa==1)THEN ! energy, electron
          str(neq,0:nrmax) = (poh(0:nrmax)+0.5d0*pnb(0:nrmax))/(rkev*1.d20) 
       ELSE IF(nva==3 .AND. nsa==2)THEN ! energy, ion
          str(neq,0:nrmax) = 0.5d0*pnb(0:nrmax)/(rkev*1.d20)
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
    USE trcalv, ONLY: coulog ! coulomb logarithm
    IMPLICIT NONE
    INTEGER(ikind) :: nr,neq,neq1,nsa,nsa1,ns,ns1
    REAL(rkind) :: coef1,coef2,ams,ams1

    ctr_ex(1:neqmax,1:neqmax,0:nrmax)=0.D0

    coef1 = aee**4.d0 / (3.d0*SQRT(2.d0)*pi**1.5d0 * eps0**2.d0)

    ! *** only heat exchange ***
    ! *** relaxation process (Tokamaks 3rd p.68) ***
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
    USE trcomm, ONLY: rkev,nrmax,nsamax,neqmax,nva_neq,ph0,phs,rhog,ra, &
         dvrho,str_simple,joh,eta,ezoh,nrd1,nrd2,nrd3,nrd4, &
         poh,pnb,pnb_tot,pnb_r0,pnb_rw,pnb_eng,snb
    IMPLICIT NONE
    INTEGER(ikind) :: nr, neq

    REAL(rkind) :: sum_nb,pnb0,rhom,dvrhom,drhog

    str_simple = 0.d0
    DO nr = 0, nrmax
       DO neq=1,neqmax
          IF(nva_neq(neq) == 3) THEN
             str_simple(neq,nr) = phs+(ph0-phs)*(1.D0-(rhog(nr)/ra)**2)*1.d6
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


    ! --- for the time being ---
    ! simple RF : Gaussian profile : only energy
    sum_nb = 0.d0
    DO nr = 1, nrmax
       rhom   = 0.5d0*(rhog(nr)+rhog(nr-1))
       dvrhom = 0.5d0*(dvrho(nr)+dvrho(nr-1))
       drhog  = rhog(nr) - rhog(nr-1)

       sum_nb = sum_nb &
              + DEXP( -((ra*rhom-pnb_r0)/pnb_rw)**2 ) *dvrhom*drhog
    ENDDO
    pnb0 = pnb_tot*1.d6 / sum_nb
    DO nr = 0, nrmax
       pnb(nr) = pnb0*DEXP(-((ra*rhog(nr)-pnb_r0)/pnb_rw)**2)
       snb(nr) = pnb(nr)/(pnb_eng*rkev)*1.d-20
    ENDDO

!    nrd1(0:nrmax) = eta(0:nrmax)
!    nrd2(0:nrmax) = ezoh(0:nrmax)
!    nrd3(0:nrmax) = poh(0:nrmax)
!    nrd4(0:nrmax) = joh(0:nrmax)

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

       dtr(1,1,nr) = etam*ttrhom/(rmu0*arrhom*dvrhom)
    END DO

  END SUBROUTINE tr_calc_mag_diff

END MODULE trcalc
