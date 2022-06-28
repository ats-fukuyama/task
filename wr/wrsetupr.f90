!   wrsetupr.f90

MODULE wrsetupr

  PRIVATE
  PUBLIC wr_setup_ray

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE wr_setup_ray(ierr)

    USE wrcomm
    USE wrexecr
    USE wrexecb
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    DO NRAY=1,NRAYMAX
       CALL wr_exec_single_ray(NRAY,IERR)
       IF(IERR.NE.0) cycle

       WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
            '    RKRI=  ',RKRI,  '  PABS/PIN=', &
            1.D0-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
    ENDDO

    CALL wr_calc_pwr

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

    RETURN
  END SUBROUTINE wr_exec

  ! --- calculation of radial profile of power absorption ---

  SUBROUTINE wr_calc_pwr

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: pl_mag_old,pl_rzsu
    USE libgrf
    IMPLICIT NONE
    INTEGER,PARAMETER:: nsum=201
    REAL(rkind),ALLOCATABLE:: rsu(:),zsu(:)
    INTEGER:: nrs,nray,nstp,nrs1,nrs2,ndrs,locmax
    INTEGER:: nrl1,nrl2,ndrl,nsumax,nsu,nrl
    REAL(rkind):: drs,xl,yl,zl,rs1,rs2,sdrs,delpwr,pwrmax,dpwr,ddpwr
    REAL(rkind):: rlmin,rlmax,drl,rl1,rl2,sdrl
    INTEGER,SAVE:: nrsmax_save=0,nrlmax_save=0,nraymax_save=0,nstpmax_save=0

!   ----- evaluate plasma major radius range -----

    ALLOCATE(rsu(nsum),zsu(nsum))
    CALL pl_rzsu(rsu,zsu,nsum,nsumax)
    rlmin=rsu(1)
    rlmax=rsu(1)
    DO nsu=2,nsumax
       rlmin=MIN(rlmin,rsu(nsu))
       rlmax=MAX(rlmax,rsu(nsu))
    END DO
    DEALLOCATE(rsu,zsu)

    !   --- allocate variables for power deposition profile ---

    IF(nstpmax.NE.nstpmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(rs_nstp_nray)) DEALLOCATE(rs_nstp_nray)
       IF(ALLOCATED(rl_nstp_nray)) DEALLOCATE(rl_nstp_nray)
       ALLOCATE(rs_nstp_nray(nstpmax,nraymax),rl_nstp_nray(nstpmax,nraymax))
    END IF
    IF(nrsmax.NE.nrsmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_nrs)) DEALLOCATE(pos_nrs)
       IF(ALLOCATED(pwr_nrs)) DEALLOCATE(pwr_nrs)
       IF(ALLOCATED(pwr_nrs_nray)) DEALLOCATE(pwr_nrs_nray)
       ALLOCATE(pos_nrs(nrsmax),pwr_nrs(nrsmax),pwr_nrs_nray(nrsmax,nraymax))
    END IF
    IF(nrlmax.NE.nrlmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_nrl)) DEALLOCATE(pos_nrl)
       IF(ALLOCATED(pwr_nrl)) DEALLOCATE(pwr_nrl)
       IF(ALLOCATED(pwr_nrl_nray)) DEALLOCATE(pwr_nrl_nray)
       ALLOCATE(pos_nrl(nrlmax),pwr_nrl(nrlmax),pwr_nrl_nray(nrlmax,nraymax))
    END IF
    IF(nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_pwrmax_rs_nray)) DEALLOCATE(pos_pwrmax_rs_nray)
       IF(ALLOCATED(pwrmax_rs_nray)) DEALLOCATE(pwrmax_rs_nray)
       IF(ALLOCATED(pos_pwrmax_rl_nray)) DEALLOCATE(pos_pwrmax_rl_nray)
       IF(ALLOCATED(pwrmax_rl_nray)) DEALLOCATE(pwrmax_rl_nray)
       ALLOCATE(pos_pwrmax_rs_nray(nraymax),pwrmax_rs_nray(nraymax))
       ALLOCATE(pos_pwrmax_rl_nray(nraymax),pwrmax_rl_nray(nraymax))
    END IF
    nrsmax_save=nrsmax
    nrlmax_save=nrlmax
    nraymax_save=nraymax
    nstpmax_save=nstpmax

!     ----- Setup for RADIAL DEPOSITION PROFILE (Minor radius) -----

    drs=1.D0/nrsmax
    DO nrs=1,nrsmax
       pos_nrs(nrs)=(DBLE(nrs)-0.5D0)*drs
    ENDDO
    DO nray=1,nraymax
       DO nrs=1,nrsmax
          pwr_nrs_nray(nrs,nray)=0.D0
       ENDDO
    ENDDO
    DO nrs=1,nrsmax
       pwr_nrs(nrs)=0.D0
    ENDDO

!     ----- Setup for RADIAL DEPOSITION PROFILE (Major radius) -----

    drl=(rlmax-rlmin)/(nrlmax-1)
    DO nrl=1,nrlmax
       pos_nrl(nrl)=rlmin+(nrl-1)*drl
    ENDDO
    DO nray=1,nraymax
       DO nrl=1,nrlmax
          pwr_nrl_nray(nrl,nray)=0.D0
       ENDDO
    ENDDO
    DO nrl=1,nrlmax
       pwr_nrl(nrl)=0.D0
    ENDDO

!   --- calculate power deposition density ----

    DO nray=1,nraymax
       DO nstp=0,nstpmax_nray(nray)-1

          ! --- minor radius porfile ---
          
          xl=rays(1,nstp,nray)
          yl=rays(2,nstp,nray)
          zl=rays(3,nstp,nray)
          CALL pl_mag_old(xl,yl,zl,rs1)
          xl=rays(1,nstp+1,nray)
          yl=rays(2,nstp+1,nray)
          zl=rays(3,nstp+1,nray)
          CALL pl_mag_old(xl,yl,zl,rs2)
          IF(rs1.LE.1.D0.OR.rs1.LE.1.D0) THEN
             nrs1=INT(rs1/drs)+1
             nrs2=INT(rs2/drs)+1
             IF(nrs1.GT.nrsmax) THEN
                nrs1=nrsmax
                IF(nrs2.GT.nrsmax) EXIT
             ENDIF
             IF(nrs2.GT.nrsmax) nrs2=nrsmax
                   
             ndrs=ABS(nrs2-nrs1)
             IF(ndrs.EQ.0) THEN
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +rays(8,nstp+1,nray)
             ELSE IF(nrs1.lt.nrs2) THEN
                sdrs=(rs2-rs1)/drs
                delpwr=rays(8,nstp+1,nray)/sdrs
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +(DBLE(nrs1)-rs1/drs)*delpwr
                DO nrs=nrs1+1,nrs2-1
                   pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray)+delpwr
                ENDDO
                pwr_nrs_nray(nrs2,nray)=pwr_nrs_nray(nrs2,nray) &
                     +(rs2/drs-DBLE(nrs2-1))*delpwr
             ELSE
                sdrs=(rs1-rs2)/drs
                delpwr=rays(8,nstp+1,nray)/sdrs
                pwr_nrs_nray(nrs2,nray)=pwr_nrs_nray(nrs2,nray) &
                     +(DBLE(nrs2)-rs2/drs)*delpwr
                DO nrs=nrs2+1,nrs1-1
                   pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray)+delpwr
                ENDDO
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +(rs1/drs-DBLE(nrs1-1))*delpwr
             END IF
          ENDIF

          ! --- major radius profile ---

          xl=rays(1,nstp,nray)
          yl=rays(2,nstp,nray)
          rl1=SQRT(xl**2+yl**2)
          nrl1=INT((rl1-rlmin)/drl)+1
          xl=rays(1,nstp+1,nray)
          yl=rays(2,nstp+1,nray)
          rl2=SQRT(xl**2+yl**2)
          nrl2=INT((rl2-rlmin)/drl)+1

          IF((nrl1.GE.1.and.nrl1.LE.nrlmax).OR. &
             (nrl2.GE.1.and.nrl2.LE.nrlmax)) THEN
             IF(nrl1.LT.1) nrl1=1
             IF(nrl1.GT.nrlmax) nrl1=nrlmax
             IF(nrl2.LT.1) nrl2=1
             IF(nrl2.GT.nrlmax) nrl2=nrlmax
             ndrl=ABS(nrl2-nrl1)
             IF(ndrl.EQ.0) THEN
                pwr_nrl_nray(nrl1,nray)=pwr_nrl_nray(nrl1,nray) &
                     +rays(8,nstp+1,nray)
             ELSE IF(nrl1.LT.nrl2) THEN
                sdrl=(rl2-rl1)/drl
                delpwr=rays(8,nstp+1,nray)/sdrl
                pwr_nrl_nray(nrl1,nray)=pwr_nrl_nray(nrl1,nray) &
                     +(DBLE(nrl1)-(rl1-rlmin)/drl)*delpwr
                DO nrl=nrl1+1,nrl2-1
                   pwr_nrl_nray(nrl,nray) &
                        =pwr_nrl_nray(nrl,nray)+delpwr
                END DO
                pwr_nrl_nray(nrl2,nray)=pwr_nrl_nray(nrl2,nray) &
                     +((rl2-rlmin)/drl-dble(nrl2-1))*delpwr
             ELSE
                sdrl=(rl1-rl2)/drl
                delpwr=rays(8,nstp+1,nray)/sdrl
                pwr_nrl_nray(nrl2,nray)=pwr_nrl_nray(nrl2,nray) &
                     +(DBLE(nrl2)-(rl2-rlmin)/drl)*delpwr
                DO nrl=nrl2+1,nrl1-1
                   pwr_nrl_nray(nrl,nray) &
                        =pwr_nrl_nray(nrl,nray)+delpwr
                END DO
                pwr_nrl_nray(nrl1,nray)=pwr_nrl_nray(nrl1,nray) &
                     +((rl1-rlmin)/drl-DBLE(nrl1-1))*delpwr
             END IF
          END IF
       ENDDO

       ! --- power divided by division area ---

       DO nrs=1,nrsmax
          pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray) &
               /(2.D0*PI*(DBLE(nrs)-0.5D0)*drs*drs)
       ENDDO
       DO nrl=1,nrlmax
          pwr_nrl_nray(nrl,nray)=pwr_nrl_nray(nrl,nray) &
               /(2.D0*PI*(DBLE(nrl)-0.5D0)*drl*drl)
       ENDDO
    ENDDO

!     ----- sum of power for each ray-------

    DO nrs=1,nrsmax
       pwr_nrs(nrs)=0.D0
       DO nray=1,nraymax
          pwr_nrs(nrs)=pwr_nrs(nrs)+pwr_nrs_nray(nrs,nray)
       ENDDO
    ENDDO

    DO nrl=1,nrlmax
       pwr_nrl(nrl)=0.D0
       DO nray=1,nraymax
          pwr_nrl(nrl)=pwr_nrl(nrl)+pwr_nrl_nray(nrl,nray)
       ENDDO
    ENDDO

!     ----- find location of absorbed power peak -----

    DO nray=1,nraymax
       pwrmax=0.D0
       locmax=0
       DO nrs=1,nrsmax
          IF(pwr_nrs_nray(nrs,nray).GT.pwrmax) THEN
             pwrmax=pwr_nrs_nray(nrs,nray)
             locmax=nrs
          ENDIF
       END DO
       IF(locmax.LE.1) THEN
          pos_pwrmax_rs_nray(nray)=pos_nrs(1)
       ELSE IF(locmax.GE.nrsmax) THEN
          pos_pwrmax_rs_nray(nray)=pos_nrs(nrsmax)
       ELSE
          dpwr =(pwr_nrs_nray(locmax+1,nray) &
                -pwr_nrs_nray(locmax-1,nray))/(2.D0*drs)
          ddpwr=(pwr_nrs_nray(locmax+1,nray) &
              -2*pwr_nrs_nray(locmax  ,nray) &
                +pwr_nrs_nray(locmax-1,nray))/drs**2
          pos_pwrmax_rs_nray(nray)=(locmax-0.5D0)/(nrsmax-1.D0)-dpwr/ddpwr
          pwrmax_rs_nray(nray)=pwrmax-dpwr**2/(2.D0*ddpwr)
       ENDIF
    END DO
   
    pwrmax=0.D0
    locmax=0
    DO nrl=1,nrlmax
       IF(pwr_nrl(nrl).GT.pwrmax) THEN
          pwrmax=pwr_nrl(nrl)
          locmax=nrl
       ENDIF
    END DO
    IF(locmax.LE.1) THEN
       pos_pwrmax_rl=pos_nrl(1)
    ELSE IF(locmax.GE.nrlmax) THEN
       pos_pwrmax_rl=pos_nrl(nrlmax)
    ELSE
       dpwr =(pwr_nrl(locmax+1) &
             -pwr_nrl(locmax-1))/(2.D0*drl)
       ddpwr=(pwr_nrl(locmax+1) &
           -2*pwr_nrl(locmax  ) &
             +pwr_nrl(locmax-1))/drl**2
       pos_pwrmax_rl=(locmax-0.5D0)/(nrlmax-1.D0)-dpwr/ddpwr
       pwrmax_rl=pwrmax-dpwr**2/(2.D0*ddpwr)
    ENDIF
   
    WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
         '    PWRMAX=',pwrmax_rl,'  AT RL =',pos_pwrmax_rl

!    CALL PAGES
!    CALL grd1d(1,pos_nrs,pwr_nrs_nray,nrsmax,nrsmax,nraymax, &
!         '@pwr-nrs vs. pos-nrs@')
!    CALL grd1d(2,pos_nrl,pwr_nrl_nray,nrlmax,nrlmax,nraymax, &
!         '@pwr-nrl vs. pos-nrl@')
!    CALL PAGEE

    RETURN
  END SUBROUTINE wr_calc_pwr

END MODULE wrexec
