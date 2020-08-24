!   obexec.f90

MODULE oblocal
  USE obcomm,ONLY: rkind
  REAL(rkind):: peng,pmu,pcangl,omega_bounce
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

    CALL ob_convert(nobt,ierr)

    RETURN
  END SUBROUTINE ob_exec

!  --- original Runge-Kutta method ---

  SUBROUTINE ob_rkft(y_in,ya,nstp_last)

    USE obcomm
    USE obprep
    USE oblocal
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: y_in(neq_max)
    REAL(rkind),INTENT(OUT):: ya(0:neq_max,0:nstp_max)
    INTEGER,INTENT(OUT):: nstp_last
    REAL(rkind):: ys(neq_max),ye(neq_max),work(neq_max,2)
    REAL(rkind):: fs(neq_max)
    INTEGER:: nstp,neq,id,nstp_lim,mode_theta
    REAL(rkind):: xs,xe,theta_init,factor
    LOGICAL:: l_end

    ! --- initialization ---
    
    xs = 0.D0
    xe = delt
    nstp_lim=MIN(INT(tmax/delt),nstp_max)

    nstp=0
    ya(0,nstp)=xs
    DO neq=1,neq_max
       ys(neq)=y_in(neq)
       ya(neq,nstp)=y_in(neq)
    ENDDO

    IF(omega_bounce.LE.1.D-80) THEN
       nstp_last=0
       RETURN
    END IF

    ! --- check initial theta-dot to identify orbit mode ---
    
    CALL ob_fdrv(xs,ys,fs)
    IF(fs(2).GT.0.D0) THEN
       mode_theta=0
    ELSE IF(fs(2).LT.0.D0) THEN
       mode_theta=1
    ELSE
       IF(mdlobw.GE.1) &
            WRITE(6,'(A)') '!! starting from bounce point in theta'
       nstp_last=0
       RETURN
    END IF

       
    theta_init=ys(2)

    IF(mdlobw.GE.1) &
         WRITE(6,'(A,A/I7,1P5E12.4)') '   nstp wb time [s]', &
         '  zeta [deg]  theta [deg]     psipn  rhopara [m]', &
         0,xs,ys(1)*180.D0/Pi,ys(2)*180.D0/Pi,ys(3)/psipa,ys(4)

    
    DO nstp = 1,nstp_lim
       CALL ODERKN(neq_max,ob_fdrv,xs,xe,1,ys,ye)

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
               WRITE(6,'(I7,1P5E12.4)') &
               nstp,xe,ye(1)*180.D0/Pi,ye(2)*180.D0/Pi,ye(3)/psipa,ye(4)
       ENDIF

       ya(0,nstp)=xe
       DO neq=1,neq_max
          ya(neq,nstp)=ye(neq)
       ENDDO
       nstp_last=nstp

       !   --- bound check ---

       l_end=.FALSE.
       SELECT CASE(mdlobc)
       CASE(0)
          IF(xe+delt.GT.tmax) l_end=.TRUE.
       CASE(1)
          SELECT CASE(mode_theta)
          CASE(0) ! increasing theta
             IF(ye(2).GT.theta_init+2.D0*Pi) THEN ! untrapped particle
                factor=((theta_init+2.D0*Pi)-ys(2))/(ye(2)-ys(2))
                xe=xs+(xe-xs)*factor
                DO neq=1,neq_max
                   ye(neq)=ys(neq)+(ye(neq)-ys(neq))*factor
                END DO
                l_end=.TRUE.
             END IF
             IF(ye(2).LT.theta_init) mode_theta=2
          CASE(1) ! decreasing theta
             IF(ye(2).LT.theta_init-2.D0*Pi) THEN ! untrapped particle
                factor=((theta_init-2.D0*Pi)-ys(2))/(ye(2)-ys(2))
                xe=xs+(xe-xs)*factor
                DO neq=1,neq_max
                   ye(neq)=ys(neq)+(ye(neq)-ys(neq))*factor
                END DO
                l_end=.TRUE.
             END IF
             IF(ye(2).GT.theta_init) mode_theta=3
          CASE(2) ! decreasing after increasing theta
             IF(ye(2).GT.theta_init) THEN ! trapped particle
                factor=(theta_init-ys(2))/(ye(2)-ys(2))
                xe=xs+(xe-xs)*factor
                DO neq=1,neq_max
                   ye(neq)=ys(neq)+(ye(neq)-ys(neq))*factor
                END DO
                l_end=.TRUE.
             END IF
          CASE(3) ! increasing after decreasing theta
             IF(ye(2).LT.theta_init) THEN ! trapped particle
                factor=(theta_init-ys(2))/(ye(2)-ys(2))
                xe=xs+(xe-xs)*factor
                DO neq=1,neq_max
                   ye(neq)=ys(neq)+(ye(neq)-ys(neq))*factor
                END DO
                l_end=.TRUE.
             END IF
          END SELECT
       END SELECT

       ya(0,nstp)=xe
       DO neq=1,neq_max
          ya(neq,nstp)=ye(neq)
       ENDDO
       nstp_last=nstp

       IF(l_end) EXIT

       ! --- prepare for next step ----

       xs=xe
       xe=xs+delt
       DO neq=1,neq_max
          ys(neq)=ye(neq)
       ENDDO
       
    END DO
!     
    IF(mdlobw.GE.1) THEN
       IF(id.EQ.0) &
               WRITE(6,'(I7,1P6E12.4)') &
               nstp,xe,ye(1)*180.D0/Pi,ye(2)*180.D0/Pi,ye(3)/psipa,ye(4)
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
    REAL(rkind),INTENT(OUT):: f(neq_max)
    REAL(rkind):: pze,pzem,zetab,thetab,psip,rhopara
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
    IF(ierr.NE.0) WRITE(6,'(A,1PE12.4,I5)') &
         'XX cal_rr_pos: x,y,ierr',thetab,psip,ierr
    CALL cal_zz_pos(thetab,psip,zz_pos,ierr)
    CALL cal_bdb_pos(thetab,psip,bb_pos,db_dthetab,db_dpsip,ierr)

    CALL cal_qps_pos(psip,qps,dqps_dpsip,ierr)
    IF(ierr.NE.0) WRITE(6,'(A,1PE12.4,I5)') &
         'XX cal_qps_pos: x,ierr',psip,ierr
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

    coef1=pzem*rhopara*bb_pos**2
    coef2=(fq+rhopara*dI)/fD
    coef3=pmu/pze+pzem*rhopara**2*bb_pos
    coef4=fI/fD
    coef5=(1.D0-rhopara*dg)/fD
    coef6=fg/fD
      
    F(1)=(coef1*coef2-coef3*coef4*db_dpsip  -coef4*dphi_dpsip)/omega_bounce
    F(2)=(coef1*coef5+coef3*coef6*db_dpsip  +coef6*dphi_dpsip)/omega_bounce
    F(3)=(           -coef3*coef6*db_dthetab+coef4*dphi_dzetab &
                     +coef3*coef4*db_dzetab -coef6*dphi_dthetab)/omega_bounce
    F(4)=(           -coef3*coef5*db_dthetab-coef5*dphi_dthetab &
                     -coef3*coef2*db_dzetab -coef2*dphi_dzetab)/omega_bounce

    IF(idebug.EQ.9) THEN
       WRITE(6,'(A,1P6E12.4)') 'coef:',coef1,coef2,coef3,coef4,coef5,coef6
       WRITE(6,'(A,1P5E12.4)') 'Y:   ',Y(1),Y(2),Y(3),Y(4),X
       WRITE(6,'(A,1P4E12.4)') 'F:   ',F(1),F(2),F(3),F(4)
    END IF
    RETURN
  END SUBROUTINE ob_fdrv

!   ***** variable conversion *****

  SUBROUTINE ob_convert(nobt,ierr)

    USE obcomm
    USE oblocal
    USE obprep
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nobt
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nstp
    REAL(rkind):: zetab,thetab,psip,rhopara,pze,pzem,rr_pos,zz_pos,bb_pos
    REAL(rkind):: qps_pos,psit_pos,fg_pos,fI_pos,phi_pos,omega_c
    REAL(rkind):: ptheta_pos,pzeta_pos,dummy,vpara_pos,vperp_pos

    DO nstp=0,nstp_max_nobt(nobt)
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
       omega_c=pzem*bb_pos

       phi_pos=0.d0
       ptheta_pos=rhopara*omega_c
       pzeta_pos= pze*(rhopara*fg_pos-psip)

       vpara_pos=SQRT(2.D0*(peng-pze*phi_pos)/(PA(ns_ob)*AMP))*pcangl
       vperp_pos=SQRT(2.D0*pmu*bb_pos/(PA(ns_ob)*AMP))

       time_ob(nstp,nobt)   = obts(0,nstp,nobt)
       zetab_ob(nstp,nobt)  = obts(1,nstp,nobt)
       thetab_ob(nstp,nobt) = obts(2,nstp,nobt)
       psip_ob(nstp,nobt)   = obts(3,nstp,nobt)
       rhopara_ob(nstp,nobt)= obts(4,nstp,nobt)
       pzeta_ob(nstp,nobt)  = pzeta_pos
       ptheta_ob(nstp,nobt) = ptheta_pos
       babs_ob(nstp,nobt)   = bb_pos
       phi_ob(nstp,nobt)    = phi_pos
       vpara_ob(nstp,nobt)  = vpara_pos
       vperp_ob(nstp,nobt)  = vperp_pos
       psit_ob(nstp,nobt)   = psit_pos
       zeta_ob(nstp,nobt)   = zetab_ob(nstp,nobt)
       rr_ob(nstp,nobt)     = rr_pos
       zz_ob(nstp,nobt)     = zz_pos
       rs_ob(nstp,nobt)     = SQRT((rr_pos-rr_axis)**2+(zz_pos-zz_axis)**2)
       theta_ob(nstp,nobt)  = ATAN2(zz_pos-zz_axis,rr_pos-rr_axis)
    END DO
    RETURN
  END SUBROUTINE ob_convert

END MODULE obexec
