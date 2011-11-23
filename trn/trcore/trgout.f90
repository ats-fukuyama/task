MODULE trgout

! This module produce graphic output

PUBLIC tr_gout
PRIVATE

CONTAINS

  SUBROUTINE tr_gout

    USE trcomm, ONLY: rkind,ikind,nrmax,nsamax,rg,rn,ru,rt,qp,dfa,pha,vca, &
         nitmax,error_it,lmaxtr,epsltr,lt_save,neqrmax,neqmax,neq_neqr, &
         nsa_neq,gvt,gvts,gvrt,gvrts,ngt
    USE libgrf,ONLY: grd1d
    IMPLICIT NONE
    CHARACTER(LEN=20):: label
    INTEGER(ikind):: neqr,neq,nsa,nit,nitmaxl,ngg,ngg_interval
    INTEGER(ikind),parameter:: nggmax=10
    REAL(rkind),DIMENSION(0:nrmax,nsamax):: rtg,dfg,phg,vcg,vg1,vg2,vg3,vg4
    REAL(rkind),DIMENSION(0:nrmax,0:nggmax):: gg1,gg2,gg3,gg4
    REAL(rkind),DIMENSION(0:ngt):: gt
    REAL(rkind),DIMENSION(0:ngt,nsamax):: gt1,gt2,gt3,gt4
    REAL(rkind),DIMENSION(nitmax):: ig,erg

! ----- current radial profile -----

    DO nsa=1,nsamax
       vg1(0:nrmax,nsa)=rn(nsa,0:nrmax)
       vg2(0:nrmax,nsa)=ru(nsa,0:nrmax)
       vg3(0:nrmax,nsa)=rt(nsa,0:nrmax)
    END DO
    vg4(0:nrmax,1)=qp(0:nrmax)

    CALL PAGES
       LABEL = '/n vs r/'
       CALL GRD1D(1,rg,vg1,nrmax+1,nrmax+1,nsamax,LABEL,0,FMIN=0.D0)
       LABEL = '/u vs r/'
       CALL GRD1D(2,rg,vg2,nrmax+1,nrmax+1,nsamax,LABEL,0)
       LABEL = '/T vs r/'
       CALL GRD1D(3,rg,vg3,nrmax+1,nrmax+1,nsamax,LABEL,0,FMIN=0.D0)
       LABEL = '/q vs r/'
       CALL GRD1D(4,rg,vg4,nrmax+1,nrmax+1,1,LABEL,0)
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
       LABEL = '/T1(t) vs r/'
       CALL GRD1D(1,rg,gg1,nrmax+1,nrmax+1,nggmax+1,LABEL,0,FMIN=0.D0)
       LABEL = '/T2(t) vs r/'
       CALL GRD1D(2,rg,gg2,nrmax+1,nrmax+1,nggmax+1,LABEL,0,FMIN=0.D0)
       LABEL = '/n1(t) vs r/'
       CALL GRD1D(3,rg,gg3,nrmax+1,nrmax+1,nggmax+1,LABEL,0,FMIN=0.D0)
       LABEL = '/qp(t) vs r/'
       CALL GRD1D(4,rg,gg4,nrmax+1,nrmax+1,nggmax+1,LABEL,0,FMIN=0.D0)
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
       CALL GRD1D(1,gt,gt1,ngt+1,ngt+1,nsamax,LABEL,0,FMIN=0.D0)
       LABEL = '/u(0) vs t/'
       CALL GRD1D(2,gt,gt2,ngt+1,ngt+1,nsamax,LABEL,0,FMIN=0.D0)
       LABEL = '/T(0) vs t/'
       CALL GRD1D(3,gt,gt3,ngt+1,ngt+1,nsamax,LABEL,0,FMIN=0.D0)
       LABEL = '/q(0),q(a) vs t/'
       CALL GRD1D(4,gt,gt4,ngt+1,ngt+1,2,LABEL,0)
    CALL PAGEE

! ----- diffusion coefficients -----

    DO neqr=1,neqrmax
       neq=neq_neqr(neqr)
       IF(neq == 0) THEN
          rtg(0:nrmax,neqr)=0.D0
          dfg(0:nrmax,neqr)=0.D0
          phg(0:nrmax,neqr)=0.D0
          vcg(0:nrmax,neqr)=0.D0
       ELSE
          nsa=nsa_neq(neq)
          rtg(0:nrmax,neqr)=rt(nsa,0:nrmax)
          dfg(0:nrmax,neqr)=dfa(neq,neq,0:nrmax)
          phg(0:nrmax,neqr)=pha(neq,0:nrmax)
          vcg(0:nrmax,neqr)=vca(neq,neq,0:nrmax)
       ENDIF
    END DO

    ! --- #1 ---
    CALL PAGES
    LABEL = '/T vs r/'
    CALL GRD1D(1,RG,rtg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/Diffusion vs r/'
    CALL GRD1D(2,RG,dfg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/Heat_pw vs r/'
    CALL GRD1D(3,RG,phg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    LABEL = '/Convection vs r/'
    CALL GRD1D(4,RG,vcg, NRMAX+1, NRMAX+1, neqrmax, LABEL, 0)
    CALL PAGEE

    
    ! --- #2 ---
    nitmaxl=MIN(nitmax,lmaxtr)
    DO nit=1,nitmaxl
       ig(nit)=dble(nit)
    END DO
    IF(nitmax > 1) THEN
       DO nit=1,nitmaxl
          ig(nit)=dble(nit)
          erg(nit)=log10(error_it(nit)+epsltr*1.D-2)
       END DO

       ! --- #2 ---
       CALL PAGES
       ! MODE = 2 ; X:LINEAR  Y:LOG
       LABEL = '/convergence vs NIT/'
       CALL GRD1D(1,ig,erg, NITMAXL, NITMAXL, 1, LABEL, 2)

       LABEL = '/LT vs r/'
       CALL GRD1D(2,rg,lt_save,nrmax+1,nrmax+1,1,LABEL,1)
       
       CALL PAGEE

    END IF



  END SUBROUTINE tr_gout
END MODULE trgout
