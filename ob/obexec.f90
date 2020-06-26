!   obexec.f90

MODULE oblocal
  USE obcomm,ONLY: rkind
  REAL(rkind):: peng,pmu,pcangl
END MODULE oblocal

MODULE obexec

  PRIVATE
  PUBLIC ob_exec

CONTAINS

!   ***** Orbit tracing module *****

  SUBROUTINE ob_exec(y_in,nobt,ierr)

    USE obcomm
    USE oblocal
    USE obprep
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: y_in(neq_max)
    INTEGER,INTENT(IN):: nobt
    INTEGER,INTENT(OUT):: ierr
    REAL(rkind),ALLOCATABLE:: y(:,:)
    REAL(rkind):: bb_pos,phi_pos
    INTEGER:: nstp,neq

    ierr=0
    ALLOCATE(y(0:neq_max,0:nstp_max))
    
    IF(mdlobq.EQ.0) THEN
       CALL ob_rkft(y_in,y,nstp)
    ELSEIF(mdlobq.EQ.1) THEN
!       CALL obrkft_ode(y_in,y,nstp)
    ELSEIF(mdlobq.EQ.2) THEN
!       CALL obsymp(y_in,y,nstp)
    ELSE
       WRITE(6,*) 'XX OBCALC: unknown mdlobq =', mdlobq
       ierr=1
       RETURN
    ENDIF

    nstp_max_nobt(nobt)=nstp
    obts(0,0,nobt)=0.D0
    DO neq=1,neq_max
       obts(neq,0,nobt)=y_in(neq)
    END DO
    DO nstp=1,nstp_max_nobt(nobt)
       DO neq=0,neq_max
          obts(neq,nstp,nobt)=y(neq,nstp)
       END DO
    END DO
    RETURN
  END SUBROUTINE ob_exec

!   ***** variable conversion *****

  SUBROUTINE ob_convert(nobt,ierr)

    USE obcomm
    USE oblocal
    USE obprep
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nobt
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nstp
    REAL(rkind):: zetab,thetab,psip,rhopara,pze,pzem,rr_pos,zz_pos,bb_pos, &
         qps_pos,psit_pos,fg_pos,fI_pos,phi_pos,ptheta_pos,pzeta_pos,dummy, &
         ppara_pos,pperp_pos

    DO nstp=1,nstp_max_nobt(nobt)
       zetab=obts(1,nstp,nobt)
       thetab=obts(2,nstp,nobt)
       psip=obts(3,nstp,nobt)
       rhopara=obts(4,nstp,nobt)

       CALL cal_rr_pos(thetab,psip,rr_pos,ierr)
       CALL cal_zz_pos(thetab,psip,zz_pos,ierr)
       CALL cal_bb_pos(thetab,psip,bb_pos,ierr)

       CALL cal_qps_pos(psip,qps_pos,dummy,ierr)
       CALL cal_psit_pos(psip,psit_pos,dummy,ierr)
       CALL cal_rbps_pos(psip,fg_pos,dummy,ierr)
       CALL cal_ritps_pos(psip,fI_pos,dummy,ierr)
    
       pze=PZ(ns_ob)*AEE
       pzem=pze/(PA(ns_ob)*AMP)

       phi_pos=0.d0
       ptheta_pos=pze*(psit_pos+rhopara*fI_pos)
       pzeta_pos= pze*(rhopara*fg_pos-psip)

       ppara_pos=SQRT(2.D0*PA(ns_ob)*AMP*(peng-pze*phi_pos))*pcangl
       pperp_pos=SQRT(2.D0*PA(ns_ob)*AMP*pmu*bb_pos)

       time_ob(nstp,nobt)   = obts(0,nstp,nobt)
       zetab_ob(nstp,nobt)  = obts(1,nstp,nobt)
       thetab_ob(nstp,nobt) = obts(2,nstp,nobt)
       psip_ob(nstp,nobt)   = obts(3,nstp,nobt)
       rhopara_ob(nstp,nobt)= obts(4,nstp,nobt)
       pzeta_ob(nstp,nobt)  = pzeta_pos
       ptheta_ob(nstp,nobt) = ptheta_pos
       babs_ob(nstp,nobt)   = bb_pos
       phi_ob(nstp,nobt)    = phi_pos
       ppara_ob(nstp,nobt)= ppara_pos
       pperp_ob(nstp,nobt)= pperp_pos
       psit_ob(nstp,nobt)   = psit_pos
       rr_ob(nstp,nobt)     = rr_pos
       zz_ob(nstp,nobt)     = zz_pos
       zeta_ob(nstp,nobt)   = zetab_ob(nstp,nobt)
       rs_ob(nstp,nobt)     = SQRT((rr_pos-rr_axis)**2+(zz_pos-zz_axis)**2)
       theta_ob(nstp,nobt)  = ATAN2(zz_pos-zz_axis,rr_pos-rr_axis)
    END DO
    RETURN
  END SUBROUTINE ob_convert

!  --- original Runge-Kutta method ---

  SUBROUTINE ob_rkft(y_in,ya,nstp_last)

    USE obcomm
    USE obprep
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: y_in(neq_max)
    REAL(rkind),INTENT(OUT):: ya(0:neq_max,0:nstp_max)
    INTEGER,INTENT(OUT):: nstp_last
    REAL(rkind):: ys(neq_max),ye(neq_max),work(neq_max,2)
    INTEGER:: nstp,neq,id,nstp_lim
    REAL(rkind):: xs,xe

    xs = 0.D0
    xe = dels
    nstp_lim=MIN(INT(smax/dels),nstp_max)

    nstp=0
    ya(0,nstp)=xs
    DO neq=1,neq_max
       ys(neq)=y_in(neq)
       ya(neq,nstp)=y_in(neq)
    ENDDO
    IF(mdlobw.GE.1) &
         WRITE(6,'(6X,A,A/I10,1P5E12.4)') 'nstp', &
         '  length [m]  zeta [deg]  theta [deg]     psipn  rhopara [m]', &
         0,xs,ys(1)*180.D0/Pi,ys(2)*180.D0/Pi,ys(3)/psipa,ys(4)

    DO nstp = 1,nstp_lim
       CALL ODERK(4,ob_fdrv,xs,xe,1,ys,ye,work)

       IF(mdlobw.GE.1) THEN
          id=0
          SELECT CASE(mdlobw)
          CASE(1)
             id=1
          CASE(2)
             IF(MOD(nstp,10).EQ.0) id=1
          CASE(3)
             IF(MOD(nstp,100).EQ.0) id=1
          CASE(4)
             IF(MOD(nstp,1000).EQ.0) id=1
          CASE(5)
             IF(MOD(nstp,10000).EQ.0) id=1
          END SELECT
          IF(id.EQ.1) &
               WRITE(6,'(I10,1P5E12.4)') &
               nstp,xe,ye(1)*180.D0/Pi,ye(2)*180.D0/Pi,ye(3)/psipa,ye(4)
       ENDIF

       ya(0,nstp)=xe
       DO neq=1,neq_max
          ya(neq,nstp)=ye(neq)
       ENDDO
       nstp_last=nstp

       !   --- bound check ---

       xs=xe
       xe=xs+dels
       DO neq=1,neq_max
          ys(neq)=ye(neq)
       ENDDO
       
    END DO
    nstp_last=nstp_lim
!     
    IF(mdlobw.GE.1) THEN
       IF(id.EQ.0) &
            WRITE(6,'(1P5E12.4)') xe,ye(1),ye(2),ye(3),ye(4)
    ENDIF

    RETURN
  END SUBROUTINE ob_rkft

!  --- slave routine for obt tracing ---

  SUBROUTINE ob_fdrv(x,y,f)

    USE obcomm
    USE oblocal
    USE obprep
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x,y(neq_max)
    REAL(rkind),INTENT(OUT):: F(neq_max)
    REAL(rkind):: pze,pzem,zetab,thetab,psip,rhopara,pepara,vpara
    REAL(rkind):: rr_pos,zz_pos,bb_pos,db_dthetab,db_dpsip
    REAL(rkind):: qps,dqps_dpsip,rbps,drbps_dpsip,ritps,dritps_dpsip
    REAL(rkind):: fg,fI,fq,dg,dI,fD
    REAL(rkind):: phi,dphi_dzetab,dphi_dthetab,dphi_dpsip,db_dzetab
    REAL(rkind):: coef1,coef2,coef3,coef4,coef5,coef6
    INTEGER:: ierr
    
    pze=PZ(ns_ob)*AEE
    pzem=pze/(PA(ns_ob)*AMP)

    zetab=y(1)
    thetab=y(2)
    psip=y(3)
    rhopara=y(4)

    CALL cal_rr_pos(thetab,psip,rr_pos,ierr)
    CALL cal_zz_pos(thetab,psip,zz_pos,ierr)
    CALL cal_bdb_pos(thetab,psip,bb_pos,db_dthetab,db_dpsip,ierr)

    CALL cal_qps_pos(psip,qps,dqps_dpsip,ierr)
    CALL cal_rbps_pos(psip,rbps,drbps_dpsip,ierr)
    CALL cal_ritps_pos(psip,ritps,dritps_dpsip,ierr)

    fg=rbps
    fI=ritps
    fq=qps
    dg=drbps_dpsip
    dI=dritps_dpsip
    
    fD=fg*fq+fI+rhopara*(fg*dI-fI*dg)
    
    phi=0.D0
    dphi_dzetab=0.D0
    dphi_dthetab=0.D0
    dphi_dpsip=0.D0
    db_dzetab=0.D0

    pepara=peng-pze*phi-pmu*bb_pos
    vpara=SQRT(2.D0*pepara/(PA(ns_ob)*AMP))

    coef1=pzem*rhopara*bb_pos**2
    coef2=(fq+rhopara*dI)/fD
    coef3=pmu/pze+pzem*rhopara**2*bb_pos
    coef4=fI/fD
    coef5=(1.D0-rhopara*dg)/fD
    coef6=fg/fD
      
    F(1)=(coef1*coef2-coef3*coef4*db_dpsip  -coef4*dphi_dpsip)/vpara
    F(2)=(coef1*coef5+coef3*coef6*db_dpsip  +coef6*dphi_dpsip)/vpara
    F(3)=(           -coef3*coef6*db_dthetab+coef4*dphi_dzetab &
                     +coef3*coef4*db_dzetab -coef6*dphi_dthetab)/vpara
    F(4)=(           -coef3*coef5*db_dthetab-coef5*dphi_dthetab &
                     -coef3*coef2*db_dzetab -coef2*dphi_dzetab)/vpara

    IF(idebug.EQ.9) THEN
       WRITE(6,'(A,1P6E12.4)') 'coef:',coef1,coef2,coef3,coef4,coef5,coef6
       WRITE(6,'(A,1P5E12.4)') 'Y:   ',Y(1),Y(2),Y(3),Y(4),X
       WRITE(6,'(A,1P4E12.4)') 'F:   ',F(1),F(2),F(3),F(4)
    END IF
    RETURN
  END SUBROUTINE ob_fdrv

END MODULE obexec
