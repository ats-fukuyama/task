MODULE trcalc

  USE bpsd_kinds

  PRIVATE
  PUBLIC tr_calc

CONTAINS

! ***** calculate transport coefficients and souce terms *****

  SUBROUTINE tr_calc
    USE trcomm, ONLY: nrmax,nsamax,neqmax,dtr,vtr,ctr,str,htr, &
         dtr_nc,vtr_nc,dtr_tb,vtr_tb,ctr_ex
    USE trcoef, ONLY: tr_coef
    IMPLICIT NONE
    INTEGER(ikind):: nr,neq

    CALL tr_coef          ! calculate transport coefficients
    CALL tr_calc_exchange ! calculate exchange rate
    CALL tr_calc_source   ! calculate source and sink term

    dtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    vtr(1:neqmax,1:neqmax,0:nrmax)=0.D0
    ctr(1:neqmax,1:neqmax,0:nrmax)=0.D0
!    htr(1:neqmax,0:nrmax)=0.D0
!    str(1:neqmax,0:nrmax)=0.D0


    DO nr=1,nrmax
       dtr(1,1,nr)=0.D0       ! registivity term
       htr(1,nr)=0.D0         ! driven current term

!       str(2:neqmax,nr)=0.D0  ! source and sink term
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

    RETURN
  END SUBROUTINE tr_calc

! ----- calculate exchange rate -----

  SUBROUTINE tr_calc_exchange
    USE trcomm, ONLY: nrmax,nsamax,ctr_ex
    IMPLICIT NONE

    ctr_ex(1:3*nsamax,1:3*nsamax,0:nrmax)=0.D0

    RETURN
  END SUBROUTINE tr_calc_exchange

! ----- calculate source -----

  SUBROUTINE tr_calc_source
    USE trcomm, ONLY: nrmax,nsamax,neqmax,ph0,phs,rhog,ra,str
    IMPLICIT NONE
    INTEGER(ikind) :: nr, neq

    DO nr = 0, nrmax
       DO neq=1,neqmax
          str(neq,nr) = phs+(ph0-phs)*(1.D0-(rhog(nr)/ra)**2)
!          str(neq,nr) = phs+(ph0-phs)*(1.D0-(rg(nr)/ra)**2)
       END DO
    END DO
    RETURN
  END SUBROUTINE tr_calc_source

END MODULE trcalc
