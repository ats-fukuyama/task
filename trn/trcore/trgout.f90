MODULE trgout

! This module produce graphic output

PUBLIC tr_gout
PRIVATE

CONTAINS

  SUBROUTINE tr_gout

    USE trcomm, ONLY: rkind,ikind,nrmax,nsamax,rg,rhog,rhom, &
         rn,ru,rt,qp,dtr,str,vtr, &
         nitmax,error_it,lmaxtr,epsltr,lt_save,neqrmax,neqmax,neq_neqr, &
         nsa_neq,gvt,gvts,gvrt,gvrts,gparts,ngt, &
         mdltr_prv,dtr_prv,vtr_prv,er
    USE libgrf,ONLY: grd1d,grf1d
    IMPLICIT NONE
    CHARACTER(LEN=20):: label
    INTEGER(ikind):: nr,neqr,neq,nsa,nit,nitmaxl,ngg,ngg_interval
    INTEGER(ikind),parameter:: nggmax=10
    REAL(rkind),DIMENSION(1:nrmax) :: rhomg
    REAL(rkind),DIMENSION(0:nrmax,nsamax):: rtg,phg,vg1,vg2,vg3,vg4
    REAL(rkind),DIMENSION(1:nrmax,nsamax):: dfg,vcg
    REAL(rkind),DIMENSION(0:nrmax,nsamax):: pdg,pvg
    REAL(rkind),DIMENSION(0:nrmax,0:nggmax):: gg1,gg2,gg3,gg4,gparg1,gparg2
    REAL(rkind),DIMENSION(0:ngt):: gt
    REAL(rkind),DIMENSION(0:ngt,nsamax):: gt1,gt2,gt3,gt4
    REAL(rkind),DIMENSION(nitmax):: ig,erg
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: temp

    REAL(rkind),DIMENSION(3*neqmax,1:nrmax) :: dtrg,vtrg

! ----- current radial profile -----

    DO nsa=1,nsamax
       vg1(0:nrmax,nsa)=rn(nsa,0:nrmax)
       vg2(0:nrmax,nsa)=ru(nsa,0:nrmax)
       vg3(0:nrmax,nsa)=rt(nsa,0:nrmax)
    END DO
!       vg4(0:nrmax,1)=qp(0:nrmax
       vg4(0:nrmax,1)=er(0:nrmax)

    CALL PAGES
       LABEL = '/n vs rho/'
       CALL GRD1D(1,rhog,vg1,nrmax+1,nrmax+1,nsamax,LABEL,0)
       LABEL = '/u vs rho/'
       CALL GRD1D(2,rhog,vg2,nrmax+1,nrmax+1,nsamax,LABEL,0)
       LABEL = '/T vs rho/'
       CALL GRD1D(3,rhog,vg3,nrmax+1,nrmax+1,nsamax,LABEL,0)
       LABEL = '/q vs rho/'
       CALL GRD1D(4,rhog,vg4,nrmax+1,nrmax+1,1,LABEL,0)
    CALL PAGEE

! ----- history of radial profile -----

    ngg_interval=ngt/(MOD(ngt-1,nggmax)+1)
    DO ngg=0,nggmax
       gg1(0:nrmax,ngg)=gvrts(0:nrmax,ngg*ngg_interval,1,3)
       gg2(0:nrmax,ngg)=gvrts(0:nrmax,ngg*ngg_interval,2,3)
       gg3(0:nrmax,ngg)=gvrts(0:nrmax,ngg*ngg_interval,1,1)
       gg4(0:nrmax,ngg)=gvrt(0:nrmax,ngg*ngg_interval,1)
    END DO

    CALL PAGES
       LABEL = '/T1(t) vs rho/'
       CALL GRD1D(1,rhog,gg1,nrmax+1,nrmax+1,nggmax+1,LABEL,0)
       LABEL = '/T2(t) vs rho/'
       CALL GRD1D(2,rhog,gg2,nrmax+1,nrmax+1,nggmax+1,LABEL,0)
       LABEL = '/n1(t) vs rho/'
       CALL GRD1D(3,rhog,gg3,nrmax+1,nrmax+1,nggmax+1,LABEL,0)
       LABEL = '/qp(t) vs rho/'
       CALL GRD1D(4,rhog,gg4,nrmax+1,nrmax+1,nggmax+1,LABEL,0)
    CALL PAGEE

! ----- time evolution -----

    gt(0:ngt)=gvt(0:ngt,0)
    DO nsa=1,nsamax
       gt1(0:ngt,nsa)=gvts(0:ngt,nsa,1)
       gt2(0:ngt,nsa)=gvts(0:ngt,nsa,2)
       gt3(0:ngt,nsa)=gvts(0:ngt,nsa,3)
    END DO
    gt4(0:ngt,1)=gvt(0:ngt,1)

    gt4(0:ngt,2)=gvt(0:ngt,2)
    
    CALL PAGES
       LABEL = '/n(0) vs t/'
       CALL GRD1D(1,gt,gt1,ngt+1,ngt+1,nsamax,LABEL,0)
       LABEL = '/u(0) vs t/'
       CALL GRD1D(2,gt,gt2,ngt+1,ngt+1,nsamax,LABEL,0)
       LABEL = '/T(0) vs t/'
       CALL GRD1D(3,gt,gt3,ngt+1,ngt+1,nsamax,LABEL,0)
       LABEL = '/q(0),q(a) vs t/'
       CALL GRD1D(4,gt,gt4,ngt+1,ngt+1,2,LABEL,0)
    CALL PAGEE

! ----- diffusion coefficients -----

    rhomg(1:nrmax) = rhom(1:nrmax)
    DO neqr=1,neqrmax
       neq=neq_neqr(neqr)
       IF(neq == 0) THEN
          rtg(0:nrmax,neqr)=0.D0
          dfg(1:nrmax,neqr)=0.D0
          phg(0:nrmax,neqr)=0.D0
          vcg(1:nrmax,neqr)=0.D0
       ELSE
          IF(mdltr_prv /= 0)THEN
          ! for Pereverzev method
             DO nr=1,nrmax
                dtrg(neq,nr) = dtr(neq,neq,nr) - dtr_prv(neq-1,nr)
                vtrg(neq,nr) = vtr(neq,neq,nr) - vtr_prv(neq-1,nr)
             END DO
             nsa=nsa_neq(neq)
             rtg(0:nrmax,neqr)=rt(nsa,0:nrmax)
             dfg(1:nrmax,neqr)=MIN(dtrg(neq,1:nrmax),20.D0)
             phg(0:nrmax,neqr)=str(neq,0:nrmax)
             vcg(1:nrmax,neqr)=vtrg(neq,1:nrmax)
          ELSE
             nsa=nsa_neq(neq)
             rtg(0:nrmax,neqr)=rt(nsa,0:nrmax)
             dfg(1:nrmax,neqr)=MIN(dtr(neq,neq,1:nrmax),20.D0)
!             dfg(1:nrmax,neqr)=dtr(neq,neq,1:nrmax)
             phg(0:nrmax,neqr)=str(neq,0:nrmax)
             vcg(1:nrmax,neqr)=vtr(neq,neq,1:nrmax)
          END IF
       ENDIF
    END DO

    CALL PAGES
    LABEL = '/T vs rho/'
    CALL GRD1D(1,rhog,rtg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/Diffusion vs rho/'
    CALL GRD1D(2,rhomg,dfg, NRMAX, NRMAX, neqrmax, LABEL, 0)
    LABEL = '/Heat_pw vs rho/'
    CALL GRD1D(3,rhog,phg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/Convection vs rho/'
    CALL GRD1D(4,rhomg,vcg, NRMAX, NRMAX, neqrmax, LABEL, 0)
    CALL PAGEE

!--- for Pereverzev method ---

    IF(mdltr_prv == 0)THEN
       pdg     = 0.D0
       pvg     = 0.D0
    ELSE
       DO neqr = 1, neqrmax
          neq=neq_neqr(neqr)       
          IF(neq == 0)THEN
             pdg(0:nrmax,neqr) = 0.D0
             pvg(0:nrmax,neqr) = 0.D0
          ELSE
             nsa=nsa_neq(neq)          
             pdg(0:nrmax,neqr) = dtr_prv(neq-1,0:nrmax)
             pvg(0:nrmax,neqr) = vtr_prv(neq-1,0:nrmax)
          END IF
       END DO
    END IF
    
    CALL PAGES
    LABEL = '/T vs rho/'
    CALL GRD1D(1,rhog,rtg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/add_Diff vs rho/'
    CALL GRD1D(2,rhog,pdg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/add_Conv vs rho/'
    CALL GRD1D(3,rhog,pvg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    CALL PAGEE

!--- history of additional diffusive profile

    ngg_interval = ngt/(MOD(ngt-1,nggmax)+1)
    DO ngg = 0, nggmax
       ! thermal diffusivity 
       gparg1(0:nrmax,ngg) = gparts(0:nrmax,ngg*ngg_interval,1,3)
       gparg2(0:nrmax,ngg) = gparts(0:nrmax,ngg*ngg_interval,2,3)
    END DO

    CALL PAGES
    LABEL = '/add_Net(T_e) vs rho'
    CALL GRD1D(1,rhog,gparg1,nrmax+1,nrmax+1,nggmax+1,LABEL,0)
    LABEL = '/add_Net(T_D) vs rho'
    CALL GRD1D(2,rhog,gparg2,nrmax+1,nrmax+1,nggmax+1,LABEL,0)
    CALL PAGEE

!--- convergence
    
    nitmaxl = MIN(nitmax,lmaxtr)
    DO nit=1,nitmaxl
       ig(nit) = dble(nit)
    END DO
    IF(nitmax > 1) THEN
       DO nit=1,nitmaxl
          ig(nit)  = dble(nit)
          erg(nit) = log10(error_it(nit)+epsltr*1.D-2)
       END DO

       CALL PAGES
       ! MODE = 2 ; X:LINEAR  Y:LOG
       LABEL = '/convergence vs NIT/'
       CALL GRD1D(1,ig,erg, NITMAXL, NITMAXL, 1, LABEL, 2)

       LABEL = '/LT vs rho/'
       allocate(temp(0:nrmax,1:nsamax))
!       write(*,*) 'XXX nsamax = ',nsamax
!       temp(0:nrmax,1:nsamax)=lt_save(1:nsamax,0:nrmax)
       temp = transpose(lt_save)

       CALL GRD1D(2,rhog,temp,nrmax+1,nrmax+1,nsamax,LABEL,1)
       deallocate(temp)
       
       CALL PAGEE

    END IF

  END SUBROUTINE tr_gout
END MODULE trgout
