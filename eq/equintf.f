      SUBROUTINE equ_set_var1(nsr,nsz,nv,nsu,ilimt,btv,saxis,ell,trg)
      INCLUDE '../eq/eqcomq.inc'
      INTEGER  nsr,nsz,nv,nsu,ilimt
      REAL*8 btv,saxis,ell,trg

      NRGMAX=nsr
      NZGMAX=nsz
      NPSMAX=nv
      Bctr=btv
      PSI0=saxis
      PSIA=0.D0
      RKAP=ell
      RDLT=trg
      NSUMAX=nsu
      limitr=ilimt
      RETURN
      END

      SUBROUTINE equ_set_psi(rg_,zg_,psi_)

      INCLUDE '../eq/eqcomq.inc'
      REAL*8 rg_(NRGMAX),zg_(NZGMAX),psi_(NRGMAX*NZGMAX)
      INTEGER NZG,NRG

      DO NRG=1,NRGMAX
         RG(NRG)=rg_(NRG)
      END DO

      DO NZG=1,NZGMAX
         ZG(NZG)=zg_(NZG)
      END DO

      DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZ(NRG,NZG)=psi_((NZG-1)*NRGMAX+NRG)
         END DO
      END DO
      RETURN
      END

      SUBROUTINE equ_set_fluxfn(pds,fds,vlv,qqv,prv,xxx)

      INCLUDE '../eq/eqcomq.inc'
      REAL*8 pds(NPSMAX),fds(NPSMAX),vlv(NPSMAX)
      REAL*8 qqv(NPSMAX),prv(NPSMAX),xxx(NPSMAX)
      INTEGER NPS

      DO NPS=1,NPSMAX
         PSIPS(NPS)=xxx(NPS)
         PPPS(NPS)=prv(NPS)
      END DO
      RETURN
      END

