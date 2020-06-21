! obprep.f90

MODULE obprep

  USE obcomm,ONLY: rkind
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       rhot,psip,psit,pps,qps,rbps,vps,rlen,ritps,rsu,zsu,rsw,zsw
  REAL(rkind),Allocatable,DIMENSION(:,:):: &
       rpsi,drpsi,drchi,zpsi,dzpsi,dzchi,bpr,bpz,bpt,btp
  REAL(rkind),Allocatable,DIMENSION(:,:):: &
       rg11,rg12,rg13,rg22,rg23,rg33,rj,bfld2,bfld3
  REAL(rkind),ALLOCATABLE,DIMENSION(:)::  &
       chi
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:,:,:):: &
       urpsi,uzpsi,ubpsi
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: &
       uqps,urbps,uritps
  REAL(rkind):: &
       raxis,zaxis,psipa,psita,q0,qa
  
  PRIVATE
  PUBLIC ob_prep,cal_r_pos,cal_z_pos,cal_b_pos, &
         cal_qps_pos,cal_rbps_pos,cal_ritps_pos

CONTAINS
  
!     ***** Load equilibrium data *****

  SUBROUTINE ob_eqload(ierr)
    USE obcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    CHARACTER(LEN=80) :: line
    INTEGER:: nr,nth

    ierr=0

    ! --- if no change eq parameter, skip loading eq data ---
    
    IF(nthmax_ob.LT.4) THEN
       WRITE(6,*) 'XX ob_eqload: nthmax must be greater than 4'
       ierr=1
       RETURN
    ENDIF

    ! --- load equilibrium data from 'knameq' ---
    
    CALL eqload(modelg,knameq,ierr)
    IF(ierr.NE.0) RETURN

    ! --- set eq paramters for boozer coordinates' ---
    
    WRITE(line,'(A)') 'mdleqc=1'   ! set boozer poloidal angle
    CALL eqparm(2,line,ierr)
    IF(ierr.NE.0) RETURN
    
    WRITE(line,'(A,I5)') 'nrmax=',nrmax_ob
    CALL eqparm(2,line,ierr)
    IF(ierr.NE.0) RETURN
    
    WRITE(line,'(A,I5)') 'nthmax=',nthmax_ob
    CALL eqparm(2,line,ierr)
    IF(ierr.NE.0) RETURN
    
    WRITE(line,'(A,I5)') 'nsumax=',nsumax_ob
    CALL eqparm(2,line,ierr)
    IF(ierr.NE.0) RETURN

    CALL eqcalq(ierr)
    IF(ierr.NE.0) RETURN

    ! --- allocate eq variables ---
    
    IF(ALLOCATED(rhot)) THEN
       DEALLOCATE(rhot,psip,psit,pps,qps,rbps,vps,rlen,ritps,rsu,zsu,rsw,zsw)
       DEALLOCATE(rpsi,drpsi,drchi,zpsi,dzpsi,dzchi,bpr,bpz,bpt,btp)
       DEALLOCATE(rg11,rg12,rg13,rg22,rg23,rg33,rj,bfld2,bfld3)
       DEALLOCATE(urpsi,uzpsi,ubpsi)
       DEALLOCATE(uqps,urbps,uritps)
       DEALLOCATE(chi)
    END IF
    ALLOCATE(rhot(nrmax_ob),psip(nrmax_ob),psit(nrmax_ob))
    ALLOCATE(rpsi(nthmax_ob,nrmax_ob))
    ALLOCATE(drpsi(nthmax_ob,nrmax_ob))
    ALLOCATE(drchi(nthmax_ob,nrmax_ob))
    ALLOCATE(zpsi(nthmax_ob,nrmax_ob))
    ALLOCATE(dzpsi(nthmax_ob,nrmax_ob))
    ALLOCATE(dzchi(nthmax_ob,nrmax_ob))
    ALLOCATE(bpr(nthmax_ob,nrmax_ob))
    ALLOCATE(bpz(nthmax_ob,nrmax_ob))
    ALLOCATE(bpt(nthmax_ob,nrmax_ob))
    ALLOCATE(btp(nthmax_ob,nrmax_ob))
    ALLOCATE(pps(nrmax_ob),qps(nrmax_ob),rbps(nrmax_ob))
    ALLOCATE(vps(nrmax_ob),rlen(nrmax_ob),ritps(nrmax_ob))
    ALLOCATE(rsu(nsumax_ob),zsu(nsumax_ob),rsw(nsumax_ob),zsw(nsumax_ob))
    ALLOCATE(rg11(nthmax_ob,nrmax_ob),rg12(nthmax_ob,nrmax_ob))
    ALLOCATE(rg13(nthmax_ob,nrmax_ob),rg22(nthmax_ob,nrmax_ob))
    ALLOCATE(rg23(nthmax_ob,nrmax_ob),rg33(nthmax_ob,nrmax_ob))
    ALLOCATE(rJ(nthmax_ob,nrmax_ob))
    ALLOCATE(bfld2(nthmax_ob,nrmax_ob),bfld3(nthmax_ob,nrmax_ob))

    ALLOCATE(urpsi(4,4,nthmax_ob+1,nrmax_ob))
    ALLOCATE(uzpsi(4,4,nthmax_ob+1,nrmax_ob))
    ALLOCATE(ubpsi(4,4,nthmax_ob+1,nrmax_ob))
    ALLOCATE(uqps(4,nrmax_ob),urbps(4,nrmax_ob),uritps(4,nrmax_ob))
    ALLOCATE(chi(nthmax_ob+1))

    ! --- read eq data from eq ---

    CALL eqgetb(bb,rr,rip,ra,rkap,rdlt,rb) ! get basic eq. parameter
    CALL eqgetp(rhot,psip,nrmax_ob)  ! normalized psit radius
    CALL eqgetr(rpsi,drpsi,drchi,nthmax_ob,nthmax_ob,nrmax_ob) ! R and deriv
    CALL eqgetz(zpsi,dzpsi,dzchi,nthmax_ob,nthmax_ob,nrmax_ob) ! Z and deriv
    CALL eqgetbb(bpr,bpz,bpt,btp,nthmax_ob,nthmax_ob,nrmax_ob) ! mag field
    CALL eqgetqn(pps,qps,rbps,vps,rlen,ritps,nrmax_ob) ! flux functions
    CALL eqgetu(rsu,zsu,rsw,zsw,nsumax_ob) ! plasma and wall surface
    CALL eqgeta(raxis,zaxis,psipa,psita,q0,qa) ! axis and mag parameters

    ! --- define toroidal flux with 0 on magnetci axis ---

    DO nr=1,nrmax_ob
       psit(nr)=psita*rhot(nr)**2
    END DO

    ! -- define metric tensor, Jacobian and magnetic field  ---

    DO nr=2,nrmax_ob
       DO nth=1,nthmax_ob
          rg11(nth,nr)=drpsi(nth,nr)**2+dzpsi(nth,nr)**2
          rg12(nth,nr)=drpsi(nth,nr)*drchi(nth,nr) &
                      +dzpsi(nth,nr)*dzchi(nth,nr)
          rg13(nth,nr)=0.D0
          rg22(nth,nr)=drchi(nth,nr)**2+dzchi(nth,nr)**2
          rg23(nth,nr)=0.D0
          rg33(nth,nr)=rpsi(nth,nr)**2
          rj  (nth,nr)=rpsi(nth,nr)*(drpsi(nth,nr)*dzchi(nth,nr) &
                                   -drchi(nth,nr)*dzpsi(nth,nr))

          bfld2(nth,nr)=(bpr(nth,nr)*dzpsi(nth,nr) &
                         -bpz(nth,nr)*drpsi(nth,nr)) &
                         /SQRT(rg11(nth,nr)) &
                         /SQRT(rg22(nth,nr))
          bfld3(nth,nr)=btp(nth,nr)/rpsi(nth,nr)
       ENDDO
    ENDDO

    nr=1
    DO nth=1,nthmax_ob
       rg11(nth,nr)= rg11(nth,2)
       rg12(nth,nr)= rg12(nth,2)
       rg13(nth,nr)= 0.d0
       rg22(nth,nr)= rg22(nth,2)
       rg23(nth,nr)= 0.D0
       rg33(nth,nr)= rpsi(nth,nr)**2
       rj  (nth,nr)= rj(nth,2)
       bfld2(nth,nr)=bfld2(nth,2)
       bfld3(nth,nr)=btp(nth,nr)/rpsi(nth,nr)
    ENDDO
    RETURN
  END SUBROUTINE ob_eqload

  ! --- calculate spline coefficiens ---

  SUBROUTINE ob_eqspline(ierr)
    USE obcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: f,fx,fy,fxy
    REAL(rkind),ALLOCATABLE,DIMENSION(:):: fd
    REAL(rkind):: dchi
    INTEGER:: nth,nr

    ierr=0
    
    DCHI=1.D0*PI/nthmax_ob
    DO nth=1,nthmax_ob+1     ! chi=0 for nchi=1, chi=2*pi for nci=nchimax_ob+1
       CHI(nth)=dchi*(nth-1)
    END DO

    ALLOCATE(f(nthmax_ob+1,nrmax_ob))
    ALLOCATE(fx(nthmax_ob+1,nrmax_ob))
    ALLOCATE(fy(nthmax_ob+1,nrmax_ob))
    ALLOCATE(fxy(nthmax_ob+1,nrmax_ob))
    ALLOCATE(fd(nrmax_ob))
    
    DO nr=1,nrmax_ob
       DO nth=1,nthmax_ob
          f(nth,nr)=rpsi(nth,nr)
       END DO
       f(nthmax_ob+1,nr)=f(1,nr)
    END DO
    CALL SPL2D(chi,psip,f,fx,fy,fxy,urpsi,nthmax_ob+1,nthmax_ob+1,nrmax_ob, &
         4,0,ierr) ! spline R(chi,psip)
    IF(ierr.NE.0) RETURN

    DO nr=1,nrmax_ob
       DO nth=1,nthmax_ob
          f(nth,nr)=zpsi(nth,nr)
       END DO
       f(nthmax_ob+1,nr)=f(1,nr)
    END DO
    CALL SPL2D(chi,psip,f,fx,fy,fxy,uzpsi,nthmax_ob+1,nthmax_ob+1,nrmax_ob, &
               4,0,ierr) ! spline Z(chi,psip)
    IF(ierr.NE.0) RETURN

    DO nr=1,nrmax_ob
       DO nth=1,nthmax_ob
          f(nth,nr)=SQRT(bpt(nth,nr)**2+btp(nth,nr)**2)
       END DO
       f(nthmax_ob+1,nr)=f(1,nr)
    END DO
    CALL SPL2D(chi,psip,f,fx,fy,fxy,ubpsi,nthmax_ob+1,nthmax_ob+1,nrmax_ob, &
               4,0,ierr) ! spline Z(chi,psip)
    IF(ierr.NE.0) RETURN

    CALL SPL1D(psip,qps,fd,uqps,nrmax_ob,0,ierr)
    IF(ierr.NE.0) RETURN
    CALL SPL1D(psip,rbps,fd,urbps,nrmax_ob,0,ierr)
    IF(ierr.NE.0) RETURN
    CALL SPL1D(psip,ritps,fd,uritps,nrmax_ob,0,ierr)
    IF(ierr.NE.0) RETURN

    DEALLOCATE(f,fx,fy,fxy,fd)
    RETURN
  END SUBROUTINE ob_eqspline
    
! --- main subroutine for preparating orbit calculation ---
  
  SUBROUTINE ob_prep(ierr)
    USE obcomm
    IMPLICIT NONE
    INTEGER, INTENT(OUT):: ierr
    CHARACTER(LEN=80) :: knameq_save=''
    INTEGER,SAVE:: nrmax_save=0,nthmax_save=0,nsumax_save=0

    IF(nrmax_ob.EQ.nrmax_save.AND. &
       nthmax_ob.EQ.nthmax_save.AND. &
       nsumax_ob.EQ.nsumax_save.AND. &
       knameq.EQ.knameq_save) RETURN

    CALL ob_eqload(ierr)
    IF(ierr.NE.0) RETURN

    CALL ob_eqspline(ierr)
    IF(ierr.NE.0) RETURN
    
    nrmax_save=nrmax_ob
    nthmax_save=nthmax_ob
    nsumax_save=nsumax_ob
    knameq_save=knameq

    RETURN
  END SUBROUTINE ob_prep

  ! R(chi,psip)
  
  SUBROUTINE cal_r_pos(chi_pos,psip_pos,r_pos,ierr)
    USE obcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: chi_pos,psip_pos
    REAL(rkind),INTENT(OUT):: r_pos
    INTEGER,INTENT(OUT):: ierr

    CALL SPL2DF(chi_pos,psip_pos,r_pos,chi,psip,urpsi, &
         nthmax_ob+1,nthmax_ob+1,nrmax_ob,ierr)
    RETURN
  END SUBROUTINE cal_r_pos
    
  ! Z(chi,psip)
  
  SUBROUTINE cal_z_pos(chi_pos,psip_pos,z_pos,ierr)
    USE obcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: chi_pos,psip_pos
    REAL(rkind),INTENT(OUT):: z_pos
    INTEGER,INTENT(OUT):: ierr

    CALL SPL2DF(chi_pos,psip_pos,z_pos,chi,psip,uzpsi, &
         nthmax_ob+1,nthmax_ob+1,nrmax_ob,ierr)
    RETURN
  END SUBROUTINE cal_z_pos
    
  ! B(chi,psip)
  
  SUBROUTINE cal_b_pos(chi_pos,psip_pos,b_pos,db_dchi,db_dpsip,ierr)
    USE obcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: chi_pos,psip_pos
    REAL(rkind),INTENT(OUT):: b_pos,db_dchi,db_dpsip
    INTEGER,INTENT(OUT):: ierr

    CALL SPL2DD(chi_pos,psip_pos,b_pos,db_dchi,db_dpsip,ubpsi, &
         nthmax_ob+1,nthmax_ob+1,nrmax_ob,ierr)
    RETURN
  END SUBROUTINE cal_b_pos
    
  ! qps(chi,psip)
  
  SUBROUTINE cal_qps_pos(psip_pos,qps_pos,dqps_dpsip,ierr)
    USE obcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: psip_pos
    REAL(rkind),INTENT(OUT):: qps_pos,dqps_dpsip
    INTEGER,INTENT(OUT):: ierr

    CALL SPL1DD(psip_pos,qps_pos,dqps_dpsip,psip,uqps,nrmax_ob,ierr)
    
    RETURN
  END SUBROUTINE cal_qps_pos
    
  ! rbps(chi,psip)
  
  SUBROUTINE cal_rbps_pos(psip_pos,rbps_pos,drbps_dpsip,ierr)
    USE obcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: psip_pos
    REAL(rkind),INTENT(OUT):: rbps_pos,drbps_dpsip
    INTEGER,INTENT(OUT):: ierr

    CALL SPL1DD(psip_pos,rbps_pos,drbps_dpsip,psip,urbps,nrmax_ob,ierr)
    
    RETURN
  END SUBROUTINE cal_rbps_pos
    
  ! ritps(chi,psip)
  
  SUBROUTINE cal_ritps_pos(psip_pos,ritps_pos,dritps_dpsip,ierr)
    USE obcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: psip_pos
    REAL(rkind),INTENT(OUT):: ritps_pos,dritps_dpsip
    INTEGER,INTENT(OUT):: ierr

    CALL SPL1DD(psip_pos,ritps_pos,dritps_dpsip,psip,uritps,nrmax_ob,ierr)
    
    RETURN
  END SUBROUTINE cal_ritps_pos
    
END MODULE obprep

  

