! wrcalpwr.f90

MODULE wrcalpwr

  PRIVATE
  PUBLIC wr_calc_pwr
  
CONTAINS

  ! --- calculation of radial profile of power absorption ---

  SUBROUTINE wr_calc_pwr

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: pl_mag_old,pl_rzsu
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nrs,nray,nstp,nrs1,nrs2,ndrs,locmax
    INTEGER:: nrl1,nrl2,ndrl,nsu,nrl,nsa
    REAL(rkind):: drs,xl,yl,zl,rs1,rs2,sdrs,delpwr,pwrmax,dpwr,ddpwr
    REAL(rkind):: rlmin,rlmax,drl,rl1,rl2,sdrl
    REAL(rkind):: VV,TT,RF,XP,YP,ZP,RKXP,RKYP,RKZP,UU,OMG,ROMG
    REAL(rkind):: RXP,RYP,RZP,RRKXP,RRKYP,RRKZP
    REAL(rkind):: DOMG,DXP,DYP,DZP,DKXP,DKYP,DKZP,DS

!   ----- evaluate plasma major radius range -----

    CALL pl_rzsu(rsu_wr,zsu_wr,nsumax)
    rlmin=rsu_wr(1)
    rlmax=rsu_wr(1)
    DO nsu=2,nsumax
       rlmin=MIN(rlmin,rsu_wr(nsu))
       rlmax=MAX(rlmax,rsu_wr(nsu))
    END DO


!     ----- Setup for RADIAL DEPOSITION PROFILE (Minor radius) -----

    drs=1.D0/nrsmax
    DO nrs=1,nrsmax
       pos_nrs(nrs)=(DBLE(nrs)-0.5D0)*drs
    ENDDO
    DO nray=1,nraymax
       DO nrs=1,nrsmax
          DO nsa=1,nsamax_wr
             pwr_nrs_nray(nsa,nrs,nray)=0.D0
          END DO
       ENDDO
    ENDDO
    DO nrs=1,nrsmax
       DO nsa=1,nsamax_wr
          pwr_nrs(nsa,nrs)=0.D0
       END DO
    ENDDO

    !     ----- Setup for RADIAL DEPOSITION PROFILE (Major radius) -----

    drl=(rlmax-rlmin)/(nrlmax-1)
    DO nrl=1,nrlmax
       pos_nrl(nrl)=rlmin+(nrl-1)*drl
    ENDDO
    DO nray=1,nraymax
       DO nrl=1,nrlmax
          DO nsa=1,nsamax_wr
             pwr_nrl_nray(nsa,nrl,nray)=0.D0
          END DO
       ENDDO
    ENDDO
    DO nrl=1,nrlmax
       DO nsa=1,nsamax_wr
          pwr_nrl(nsa,nrl)=0.D0
       END DO
    ENDDO

!   --- calculate power deposition density ----

    DO nray=1,nraymax
       DO nstp=0,nstpmax_nray(nray)-1

          ! --- minor radius profile ---
          
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
                DO nsa=1,nsamax_wr
                   pwr_nrs_nray(nsa,nrs1,nray)=pwr_nrs_nray(nsa,nrs1,nray) &
                        +pwr_nsa_nstp_nray(nsa,nstp+1,nray)
                END DO
             ELSE IF(nrs1.lt.nrs2) THEN
                sdrs=(rs2-rs1)/drs
                DO nsa=1,nsamax_wr
                   delpwr=pwr_nsa_nstp_nray(nsa,nstp+1,nray)/sdrs
                   pwr_nrs_nray(nsa,nrs1,nray)=pwr_nrs_nray(nsa,nrs1,nray) &
                        +(DBLE(nrs1)-rs1/drs)*delpwr
                   DO nrs=nrs1+1,nrs2-1
                      pwr_nrs_nray(nsa,nrs,nray) &
                           =pwr_nrs_nray(nsa,nrs,nray)+delpwr
                   ENDDO
                   pwr_nrs_nray(nsa,nrs2,nray)=pwr_nrs_nray(nsa,nrs2,nray) &
                        +(rs2/drs-DBLE(nrs2-1))*delpwr
                END DO
             ELSE
                sdrs=(rs1-rs2)/drs
                DO nsa=1,nsamax_wr
                   delpwr=pwr_nsa_nstp_nray(nsa,nstp+1,nray)/sdrs
                   pwr_nrs_nray(nsa,nrs2,nray)=pwr_nrs_nray(nsa,nrs2,nray) &
                        +(DBLE(nrs2)-rs2/drs)*delpwr
                   DO nrs=nrs2+1,nrs1-1
                      pwr_nrs_nray(nsa,nrs,nray) &
                           =pwr_nrs_nray(nsa,nrs,nray)+delpwr
                   ENDDO
                   pwr_nrs_nray(nsa,nrs1,nray)=pwr_nrs_nray(nsa,nrs1,nray) &
                        +(rs1/drs-DBLE(nrs1-1))*delpwr
                END DO
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
                DO nsa=1,nsamax_wr
                   pwr_nrl_nray(nsa,nrl1,nray)=pwr_nrl_nray(nsa,nrl1,nray) &
                        +pwr_nsa_nstp_nray(nsa,nstp+1,nray)
                END DO
             ELSE IF(nrl1.LT.nrl2) THEN
                sdrl=(rl2-rl1)/drl
                DO nsa=1,nsamax_wr
                   delpwr=pwr_nsa_nstp_nray(nsa,nstp+1,nray)/sdrl
                   pwr_nrl_nray(nsa,nrl1,nray)=pwr_nrl_nray(nsa,nrl1,nray) &
                        +(DBLE(nrl1)-(rl1-rlmin)/drl)*delpwr
                   DO nrl=nrl1+1,nrl2-1
                      pwr_nrl_nray(nsa,nrl,nray) &
                           =pwr_nrl_nray(nsa,nrl,nray)+delpwr
                   END DO
                   pwr_nrl_nray(nsa,nrl2,nray)=pwr_nrl_nray(nsa,nrl2,nray) &
                        +((rl2-rlmin)/drl-dble(nrl2-1))*delpwr
                END DO
             ELSE
                sdrl=(rl1-rl2)/drl
                DO nsa=1,nsamax_wr
                   delpwr=pwr_nsa_nstp_nray(nsa,nstp+1,nray)/sdrl
                   pwr_nrl_nray(nsa,nrl2,nray)=pwr_nrl_nray(nsa,nrl2,nray) &
                        +(DBLE(nrl2)-(rl2-rlmin)/drl)*delpwr
                   DO nrl=nrl2+1,nrl1-1
                      pwr_nrl_nray(nsa,nrl,nray) &
                           =pwr_nrl_nray(nsa,nrl,nray)+delpwr
                   END DO
                   pwr_nrl_nray(nsa,nrl1,nray)=pwr_nrl_nray(nsa,nrl1,nray) &
                        +((rl1-rlmin)/drl-DBLE(nrl1-1))*delpwr
                END DO
             END IF
          END IF
       ENDDO
    END DO

       ! --- total power  ---

    pwr_tot=0.D0
    DO nsa=1,nsamax_wr
       pwr_nsa(nsa)=0.D0
       DO nray=1,nraymax
          pwr_nsa_nray(nsa,nray)=0.D0
          DO nrs=1,nrsmax
             pwr_nsa_nray(nsa,nray) &
                  =pwr_nsa_nray(nsa,nray)+pwr_nrs_nray(nsa,nrs,nray)
          END DO
          pwr_nsa(nsa)=pwr_nsa(nsa)+pwr_nsa_nray(nsa,nray)
       END DO
       pwr_tot=pwr_tot+pwr_nsa(nsa)
    END DO

    DO nsa=1,nsamax_wr
       DO nrs=1,nrsmax
          pwr_nrs(nsa,nrs)=0.D0
          DO nray=1,nraymax
             pwr_nrs(nsa,nrs)=pwr_nrs(nsa,nrs)+pwr_nrs_nray(nsa,nrs,nray)
          END DO
       END DO
       DO nrl=1,nrlmax
          pwr_nrl(nsa,nrl)=0.D0
          DO nray=1,nraymax
             pwr_nrl(nsa,nrl)=pwr_nrl(nsa,nrl)+pwr_nrl_nray(nsa,nrl,nray)
          END DO
       END DO
    END DO

    DO nray=1,nraymax
       WRITE(6,'(A,I4,A,4ES12.4)') 'nray=',nray, &
            '    pwr_nsa_nray=',(pwr_nsa_nray(nsa,nray),nsa=1,MIN(nsamax_wr,4))
    END DO
    
       ! --- power divided by division area ---

    DO nray=1,nraymax
       DO nrs=1,nrsmax
          DO nsa=1,nsamax_wr
             pwr_nrs_nray(nsa,nrs,nray)=pwr_nrs_nray(nsa,nrs,nray) &
                  /(2.D0*PI*(DBLE(nrs)-0.5D0)*drs*drs)
          ENDDO
       END DO
       DO nrl=1,nrlmax
          DO nsa=1,nsamax_wr
             pwr_nrl_nray(nsa,nrl,nray)=pwr_nrl_nray(nsa,nrl,nray) &
                  /(2.D0*PI*(DBLE(nrl)-0.5D0)*drl*drl)
          END DO
       ENDDO
    ENDDO

!     ----- find location of absorbed power peak -----

    DO nray=1,nraymax
       DO nsa=1,nsamax_wr
          pwrmax=0.D0
          locmax=0
          DO nrs=1,nrsmax
             IF(pwr_nrs_nray(nsa,nrs,nray).GT.pwrmax) THEN
                pwrmax=pwr_nrs_nray(nsa,nrs,nray)
                locmax=nrs
             ENDIF
          END DO
          IF(locmax.LE.1) THEN
             pos_pwrmax_rs_nray(nsa,nray)=pos_nrs(1)
          ELSE IF(locmax.GE.nrsmax) THEN
             pos_pwrmax_rs_nray(nsa,nray)=pos_nrs(nrsmax)
          ELSE
             dpwr =(pwr_nrs_nray(nsa,locmax+1,nray) &
                   -pwr_nrs_nray(nsa,locmax-1,nray))/(2.D0*drs)
             ddpwr=(pwr_nrs_nray(nsa,locmax+1,nray) &
                   -2*pwr_nrs_nray(nsa,locmax  ,nray) &
                   +pwr_nrs_nray(nsa,locmax-1,nray))/drs**2
             pos_pwrmax_rs_nray(nsa,nray) &
                  =(locmax-0.5D0)/(nrsmax-1.D0)-dpwr/ddpwr
             pwrmax_rs_nray(nsa,nray)=pwrmax-dpwr**2/(2.D0*ddpwr)
          ENDIF
   
          pwrmax=0.D0
          locmax=0
          DO nrl=1,nrlmax
             IF(pwr_nrl(nsa,nrl).GT.pwrmax) THEN
                pwrmax=pwr_nrl(nsa,nrl)
                locmax=nrl
             ENDIF
          END DO
          IF(locmax.LE.1) THEN
             pos_pwrmax_rl_nray(nsa,nray)=pos_nrl(1)
          ELSE IF(locmax.GE.nrlmax) THEN
             pos_pwrmax_rl_nray(nsa,nray)=pos_nrl(nrlmax)
          ELSE
             dpwr =(pwr_nrl(nsa,locmax+1) &
                   -pwr_nrl(nsa,locmax-1))/(2.D0*drl)
             ddpwr=(pwr_nrl(nsa,locmax+1) &
                 -2*pwr_nrl(nsa,locmax  ) &
                 +pwr_nrl(nsa,locmax-1))/drl**2
             pos_pwrmax_rl_nray(nsa,nray) &
                  =(locmax-0.5D0)/(nrlmax-1.D0)-dpwr/ddpwr
             pwrmax_rl_nray(nsa,nray)=pwrmax-dpwr**2/(2.D0*ddpwr)
          ENDIF
       END DO
    ENDDO
   
!    CALL PAGES
!    CALL grd1d(1,pos_nrs,pwr_nrs_nray,nrsmax,nrsmax,nraymax, &
!         '@pwr-nrs vs. pos-nrs@')
!    CALL grd1d(2,pos_nrl,pwr_nrl_nray,nrlmax,nrlmax,nraymax, &
!         '@pwr-nrl vs. pos-nrl@')
!    CALL PAGEE

    RETURN
  END SUBROUTINE wr_calc_pwr
END MODULE wrcalpwr
