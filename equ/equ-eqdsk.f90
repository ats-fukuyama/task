! equ-eqdsk.f90

MODULE equ_eqdsk
  PRIVATE
  PUBLIC equ_eqdsk_write
CONTAINS
  
!     ***** WRITE EQDSK FORMAT FILE from EQU *****

  SUBROUTINE equ_eqdsk_write(IERR)

    USE plcomm,ONLY: rkind,KNAMEQ,MODEFW
    USE aaa_mod
    USE geo_mod
    USE equ_mod
    USE eqv_mod
    USE par_mod
    USE vac_mod
    USE libfio
    IMPLICIT NONE

    INTEGER:: neqdsk,ierr,idum,i,iv,ir,iz,nps,npsmax,nps_save,isu,ilim
    CHARACTER(LEN=10):: string(6)
    REAL(rkind):: xdum,rmin,rmax,zmin,zmax,rdim,zdim,rctr,rleft,zmid,bctr
    REAL(rkind):: dps,ps,pv,dpv,pf,dpf
    REAL(rkind),ALLOCATABLE:: fx(:),pfv(:),psirz(:,:)
    REAL(rkind),ALLOCATABLE:: uvlv(:,:),uprv(:,:),upfv(:,:),uqqv(:,:)
    REAL(rkind),ALLOCATABLE:: pres(:),pprime(:),fpol(:),ffprime(:),qpsi(:)

    neqdsk=21
    CALL FWOPEN(neqdsk,KNAMEQ,1,MODEFW,'EQ',ierr)
    IF(ierr.NE.0) RETURN

    idum=0
    xdum=0.D0
    DO i=1,6
       string(i)='          '
    END DO

    WRITE(neqdsk,2000) (string(i),i=1,6),idum,nr,nz

    rmin=rsu(1)
    rmax=rsu(1)
    zmin=zsu(1)
    zmax=zsu(1)
    DO isu=2,nsu
       rmin=MIN(rmin,rsu(isu))
       rmax=MAX(rmax,rsu(isu))
       zmin=MIN(zmin,zsu(isu))
       zmax=MAX(zmax,zsu(isu))
    END DO

    rdim=rwmx-rwmn
    zdim=2.D0*zwmx
    rctr=0.5D0*(rwmx+rwmn)
    rleft=rwmn
    zmid=0.D0
    WRITE(6,'(A,1P4E12.4)') 'rz/minmax:',rmin,rmax,zmin,zmax
    WRITE(6,'(A,1P4E12.4)') 'rz/dimctr:',rdim,zdim,rctr,zmid
    WRITE(6,'(A,1P3E12.4)') 'rzs/axis: ',raxis,zaxis,saxis
    WRITE(neqdsk,2020) rdim,zdim,rctr,rleft,zmid
    bctr=btv/rmaj
    WRITE(neqdsk,2020) raxis,zaxis,saxis,0.D0,bctr
    WRITE(neqdsk,2020) tcur*1.D6,saxis,xdum,raxis,xdum
    WRITE(neqdsk,2020) zaxis,xdum,0.D0,xdum,xdum

    ALLOCATE(fx(nv),pfv(nv))
    ALLOCATE(uvlv(4,nv),uprv(4,nv),upfv(4,nv),uqqv(4,nv))
    fx(1:nv)=0.D0

    CALL SPL1D(siv,vlv,fx,uvlv,nv,0,ierr)
    IF(ierr.NE.0) GOTO 9010
    CALL SPL1D(vlv,prv,fx,uprv,nv,0,ierr)
    IF(ierr.NE.0) GOTO 9020
    DO iv=1,nv
       pfv(iv)=hdv(iv)/aav(iv)
!       WRITE(6,'(A,I5,1P5E12.4)') &
!            'pfv:',iv,hdv(iv),aav(iv),pfv(iv),siv(iv),hiv(iv)
    END DO
    CALL SPL1D(vlv,pfv,fx,upfv,nv,0,ierr)
    IF(ierr.NE.0) GOTO 9030
    CALL SPL1D(vlv,qqv,fx,uqqv,nv,0,ierr)
    IF(ierr.NE.0) GOTO 9040

    npsmax=nr
    ALLOCATE(pres(npsmax),pprime(npsmax),fpol(npsmax),ffprime(npsmax))
    ALLOCATE(qpsi(npsmax))
    dps=(siv(nv)-siv(1))/(npsmax-1)
    DO nps=1,npsmax
       nps_save=nps
       ps=siv(1)+dps*(nps-1)
       CALL SPL1DD(ps,pv,dpv,siv,uvlv,nv,ierr)
       IF(ierr.NE.0) GOTO 9060
       CALL SPL1DD(pv,pres(nps),pprime(nps),vlv,uprv,nv,ierr)
       pprime(nps)=pprime(nps)*dpv
       IF(ierr.NE.0) GOTO 9070
       CALL SPL1DD(pv,pf,dpf,vlv,upfv,nv,ierr)
       IF(ierr.NE.0) GOTO 9080
       fpol(nps)=pf
       ffprime(nps)=pf*dpf*dpv
       CALL SPL1DF(pv,qpsi(nps),vlv,uqqv,nv,ierr)
       IF(ierr.NE.0) GOTO 9090
       WRITE(6,'(A,I5,1P5E12.4)') &
            'nps=',nps,ps,pres(nps),pprime(nps),fpol(nps),ffprime(nps)
    END DO
       
    WRITE(neqdsk,2020) (fpol(nps),nps=1,npsmax)
    WRITE(neqdsk,2020) (pres(nps),nps=1,npsmax)
    WRITE(neqdsk,2020) (ffprime(nps),nps=1,npsmax)
    WRITE(neqdsk,2020) (pprime(nps),nps=1,npsmax)

    ALLOCATE(psirz(nr,nz))
    DO ir=1,nr
       DO iz=1,nz
          psirz(ir,iz)=psi(ir+nr*(iz-1))
       END DO
    END DO
    WRITE(neqdsk,2020) ((psirz(ir,iz),ir=1,nr),iz=1,nz)

    WRITE(neqdsk,2020) (qpsi(nps),nps=1,npsmax)

    WRITE(neqdsk,2022) nsumax,ilimt
    WRITE(neqdsk,2020) (rsu(isu),zsu(isu),isu=1,nsu)
    WRITE(neqdsk,2020) (rlimt(ilim),zlimt(ilim),ilim=1,ilimt)

    WRITE(neqdsk,2024) 0,0.D0,0
      
    DEALLOCATE(fx,pfv,psirz)
    DEALLOCATE(uvlv,uprv,upfv,uqqv)
    DEALLOCATE(pres,pprime,fpol,ffprime,qpsi)
    REWIND(neqdsk)
    CLOSE(neqdsk)
    RETURN

9010 WRITE(6,*) 'XX equ_efit: SPL1D(siv,vlv...): ierr=',ierr
    RETURN
9020 WRITE(6,*) 'XX equ_efit: SPL1D(vlv,prv...): ierr=',ierr
    RETURN
9030 WRITE(6,*) 'XX equ_efit: SPL1D(vlv,pfv...): ierr=',ierr
    RETURN
9040 WRITE(6,*) 'XX equ_efit: SPL1D(vlv,qqv...): ierr=',ierr
    RETURN
9060 WRITE(6,*) 'XX equ_efit: SPL1DD(ps,pv...): nps,ierr=',nps,ierr
    RETURN
9070 WRITE(6,*) 'XX equ_efit: SPL1DD(pv,pres...): nps,ierr=',nps,ierr
    RETURN
9080 WRITE(6,*) 'XX equ_efit: SPL1DD(pv,pf...): nps,ierr=',nps,ierr
    RETURN
9090 WRITE(6,*) 'XX equ_efit: SPL1DF(pv,qpsi...): nps,ierr=',nps,ierr
    RETURN

 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)
  END SUBROUTINE equ_eqdsk_write
END MODULE equ_eqdsk
