!     ********************************************

!           CDBM Transport model (2009/03/06)
!              Modified by M. Honda (2009/09/14)

!     ********************************************

MODULE trcdbm

  USE bpsd_kinds

  PRIVATE
  PUBLIC:: tr_cdbm

CONTAINS

  SUBROUTINE tr_cdbm

    USE trcomm, ONLY: &
         bb,rr,rg,rkap,qp,rn,rt,pa,amp,aee,cdtrn,cdtru,cdtrt,dtr_tb,vtr_tb, &
         nrmax,nsamax,ns_nsa,idnsa,mdltr_tb
    IMPLICIT NONE
    INTEGER(ikind):: nr,ns,nsa,model
    REAL(rkind):: rkev,calf,ckap,cexb,rsl,qpl,shearl,pnel,ppp,ppm, &
         dpdrl,rhoni,dvexbdrl,chi_cdbm

    dtr_tb(1:3*nsamax,1:3*nsamax,0:nrmax)=0.D0
    vtr_tb(1:3*nsamax,1:3*nsamax,0:nrmax)=0.D0

    model=mdltr_tb-130
    rkev=aee*1.D3
    calf=1.D0
    ckap=1.D0
    cexb=1.D0
    
    DO nr=1,nrmax
       rsl=0.5D0*(rg(nr-1)+rg(nr))
       qpl=0.5D0*(qp(nr-1)+qp(nr))
       shearl= (rg(nr)+rg(nr-1))*(qp(nr)-qp(nr-1)) &
                 /((rg(nr)-rg(nr-1))*(qp(nr)+qp(nr-1)))
       pnel=0.d0
       ppp=0.D0
       ppm=0.D0
       rhoni=0.D0
       DO nsa=1,nsamax
          SELECT CASE(idnsa(nsa))
          CASE(-1) ! electron
             pnel=pnel+rn(nsa,nr)
             ppp=ppp+rn(nsa,nr  )*rt(nsa,nr  )
             ppm=ppm+rn(nsa,nr-1)*rt(nsa,nr-1)
          CASE(1) ! ion
             ppp=ppp+rn(nsa,nr  )*rt(nsa,nr  )
             ppm=ppm+rn(nsa,nr-1)*rt(nsa,nr-1)
             rhoni=rhoni+pa(ns)*amp*rn(nsa,nr)*1.D20
          END SELECT
       END DO
       pnel=pnel*1.D20
       dpdrl=(ppp-ppm)*1.D20*rkev/(rg(nr)-rg(nr-1))
       dvexbdrl=0.d0

       CALL cdbm(bb,rr,rsl,rkap,qpl,shearl,pnel,rhoni,dpdrl, &
                 &    dvexbdrl,calf,ckap,cexb,model,chi_cdbm)

       DO nsa=1,nsamax
          IF(idnsa(nsa) /= 0) THEN
             dtr_tb(3*nsa-2,3*nsa-2,nr)=cdtrn*chi_cdbm
             dtr_tb(3*nsa-1,3*nsa-1,nr)=cdtru*chi_cdbm
             dtr_tb(3*nsa  ,3*nsa  ,nr)=cdtrT*chi_cdbm
          ENDIF
       END DO
    END DO
    RETURN

  END SUBROUTINE tr_cdbm
END MODULE trcdbm
