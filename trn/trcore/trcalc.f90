MODULE trcalc

  USE trcomm,ONLY: ikind,rkind

  PRIVATE
  PUBLIC tr_calc

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_calc
    USE trcomm, ONLY: nrmax,nsamax,neqmax,dtr,vtr,ctr,str,htr, &
         dtr_nc,vtr_nc,dtr_tb,vtr_tb,ctr_ex,str_simple
    USE trcoef, ONLY: tr_coef
    IMPLICIT NONE
    INTEGER(ikind):: nr,neq

    CALL tr_coef          ! calculate transport coefficients
    CALL tr_calc_exchange ! calculate exchange rate
    CALL tr_calc_source   ! calculate source and sink term

    dtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    vtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    ctr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    htr(1:neqmax,0:nrmax)=0.D0
    str(1:neqmax,0:nrmax)=0.D0

    DO nr=1,nrmax

       call tr_calc_mag_diff

       DO neq=2,neqmax
          dtr(2:neqmax,neq,nr) &
               =dtr_nc(1:neqmax-1,neq-1,nr) &
               +dtr_tb(1:neqmax-1,neq-1,nr)
          vtr(2:neqmax,neq,nr) &
               =vtr_nc(1:neqmax-1,neq-1,nr) &
               +vtr_tb(1:neqmax-1,neq-1,nr)
          ctr(2:neqmax,neq,nr) &
               =ctr_ex(1:neqmax-1,neq-1,nr)
       END DO
    END DO
    str(2:neqmax,0:nrmax) = str_simple(2:neqmax,0:nrmax)

    RETURN
  END SUBROUTINE tr_calc

! ----- calculate exchange rate -----

  SUBROUTINE tr_calc_exchange
    USE trcomm, ONLY: nrmax,neqmax,ctr_ex
    IMPLICIT NONE

    ctr_ex(1:3*neqmax,1:3*neqmax,0:nrmax)=0.D0

    RETURN
  END SUBROUTINE tr_calc_exchange

! ----- calculate source -----

  SUBROUTINE tr_calc_source
    USE trcomm, ONLY: nrmax,nsamax,neqmax,nva_neq,ph0,phs,rhog,ra,str_simple
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
    RETURN
  END SUBROUTINE tr_calc_source

! --- calculate coefficients for poloidal magnetic diffusion equation ---

  SUBROUTINE tr_calc_mag_diff
    USE trcomm, ONLY: rmu0,nrmax,dvrho,ttrho,arrho,abb1rho,eta,dtr,htr
    IMPLICIT NONE

    eta = 1.d-4
    ! registivity term (half grid)
    dtr(1,1,1:nrmax) = eta(1:nrmax)/rmu0  &
                       * ttrho(1:nrmax)*arrho(1:nrmax)/dvrho(1:nrmax)

!    dtr(1,1,1:nrmax) = 1.d0

  END SUBROUTINE tr_calc_mag_diff

END MODULE trcalc
