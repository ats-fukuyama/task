!   wrcalc.f90

MODULE wrcalc

  PRIVATE
  PUBLIC wr_calc

CONTAINS

!   ***** Ray tracing module *****

  SUBROUTINE WR_CALC(IERR)

    USE wrcomm
    USE wrexecr,ONLY: wr_execr
    USE wrsub,ONLY: wr_cold_rkr0,wrnwtn,wrcale_i
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    REAL(4):: TIME1,TIME2
    INTEGER:: NS,NSS,NRAY
    REAL(rkind):: ZA1,ZA2,SINP2,SINT2,ANGZ,ANGPH
    REAL(rkind):: RKR0_1,RKR0_2,EPARA,EPERP
    REAL(rkind):: RKR0_11,RKR0_12,RKR0_21,RKR0_22

    IERR=0

    CALL GUTIME(TIME1)

    CALL WRSETUP(IERR)

!     --- eliminate disp factor for same z/a species ---

    DO NS=1,NSMAX
       NSDP(NS)=1
       DO NSS=1,NS-1
          ZA1=PZ(NS)/PA(NS)
          ZA2=PZ(NSS)/PA(NSS)
          IF(ABS(ZA1-ZA2).LE.1.D-8) NSDP(NS)=0
       ENDDO
    ENDDO

    DO NRAY=1,NRAYMAX
       RF=RFIN(NRAY)
       RPI=RPIN(NRAY)
       ZPI=ZPIN(NRAY)
       PHII=PHIIN(NRAY)
       RKR0=RKRIN(NRAY)
       RNZI=RNZIN(NRAY)
       RNPHII=RNPHIIN(NRAY)
       ANGZ=ANGZIN(NRAY)
       ANGPH=ANGPHIN(NRAY)
       MODEW=MODEWIN(NRAY)
       UUI=UUIN(NRAY)

       IF(MDLWRI.EQ.100)THEN
          WRITE(6,*) &
               '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
       ELSE
          WRITE(6,*) &
               '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
          SINP2=SIN(ANGZ*PI/180.D0)**2
          SINT2=SIN(ANGPH*PI/180.D0)**2
          RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
          RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
          IF(ANGZ.LT.0.D0) RNZI=-RNZI
          IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
       ENDIF

       RAYIN(1,NRAY)=RF
       RAYIN(2,NRAY)=RPI
       RAYIN(3,NRAY)=ZPI
       RAYIN(4,NRAY)=PHII
       RAYIN(5,NRAY)=RKR0
       RAYIN(6,NRAY)=ANGZ
       RAYIN(7,NRAY)=ANGPH
       RAYIN(8,NRAY)=UUI

       RKZI  =2.D6*PI*RF*RNZI  /VC
       RKPHII=2.D6*PI*RF*RNPHII/VC

       IF(MODEW.NE.0) THEN
          CALL WR_COLD_RKR0(RKR0_1,RKR0_2)
          rkr0_11= rkr0_1
          rkr0_12=-rkr0_1
          rkr0_21= rkr0_2
          rkr0_22=-rkr0_2
          WRITE(6,'(A,1P2E12.4)') 'COLD: RKR0_1,RKR0_2',RKR0_1,RKR0_2
          RKR0=RKR0_11
          CALL WRCALE_I(EPARA,EPERP)
          WRITE(6,'(A,1P3E12.4)') &
               'RKR0_11,EPARA,EPERP=',RKR0_11,EPARA,EPERP
          RKR0=RKR0_12
          CALL WRCALE_I(EPARA,EPERP)
          WRITE(6,'(A,1P3E12.4)') &
               'RKR0_12,EPARA,EPERP=',RKR0_12,EPARA,EPERP
          RKR0=RKR0_21
          CALL WRCALE_I(EPARA,EPERP)
          WRITE(6,'(A,1P3E12.4)') &
               'RKR0_21,EPARA,EPERP=',RKR0_21,EPARA,EPERP
          RKR0=RKR0_22
          CALL WRCALE_I(EPARA,EPERP)
          WRITE(6,'(A,1P3E12.4)') &
               'RKR0_22,EPARA,EPERP=',RKR0_22,EPARA,EPERP

          IF(MODEW.EQ.1) THEN
             RKR0=RKR0_11
          ELSE IF(MODEW.EQ.-1) THEN
             RKR0=RKR0_12
          ELSE IF(MODEW.EQ.2) THEN
             RKR0=RKR0_21
          ELSE IF(MODEW.EQ.-2) THEN
             RKR0=RKR0_22
          END IF
       END IF

       RKRI  = RKR0
       RKZI  = 2.D6*PI*RF*RNZI  /VC
       RKPHII= 2.D6*PI*RF*RNPHII/VC
       CALL WRNWTN(IERR)  ! input RKR0,RKZI,RKPHII; output RKRI
       IF(IERR.NE.0) cycle

       CALL wr_execr(NRAY,IERR)
       IF(IERR.NE.0) cycle

       WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
            '    RKRI=  ',RKRI,  '  PABS/PIN=', &
            1.D0-RAYS(7,NSTPMAX_NRAY(NRAY),NRAY)
    ENDDO

    CALL wr_calc_pwr

    CALL GUTIME(TIME2)
    WRITE(6,*) '% CPU TIME = ',TIME2-TIME1,' sec'

    RETURN
  END SUBROUTINE WR_CALC

!     ***** ABSORBED POWER PROFILE*****

  SUBROUTINE WRSETUP(IERR)

    USE wrcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NRAY
    REAL(rkind):: ANGZ,ANGPH,SINP2,SINT2

    IERR=0

    IF(MDLWRI.EQ.0) THEN
       WRITE(6,*) &
            '# default values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
       WRITE(6,'(1PE12.4,0P7F9.2)') &
            RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
    ELSEIF(MDLWRI.EQ.1.OR.MDLWRI.EQ.2) THEN
       IF(ABS(RNZI).GT.1.D0) THEN
          ANGZ=0.D0
          ANGPH=0.D0
       ELSE
          ANGZ=ASIN(RNZI)*180.D0/PI
          IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
             ANGZ=0.D0
          ELSE
             ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
          ENDIF
          IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
             ANGPH=0.D0
          ELSE
             ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
          ENDIF
       ENDIF
       IF(MDLWRI.EQ.1) THEN
          WRITE(6,*) &
               '# default values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
       ELSEIF(MDLWRI.EQ.2) THEN
          WRITE(6,*) &
               '# default values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
          WRITE(6,'(1PE12.4,0P7F9.2)') &
               RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
       ENDIF
    ELSEIF(MDLWRI.EQ.11) THEN
       WRITE(6,*) &
            '# default values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
       WRITE(6,'(1PE12.4,0P7F9.2)') &
            RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
       PHII=0.D0
    ENDIF

!     --- Each ray tracing ---

    DO NRAY=1,NRAYMAX

       IF(MDLWRI.LT.100) THEN

1         WRITE(6,*) '# NRAY = ',NRAY
          IF(MDLWRI.EQ.0) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,PHII,RKR0,RNZI,RNPHII,UUI
             IF(ABS(RNZI).GT.1.D0) THEN
                ANGZ=0.D0
                ANGPH=0.D0
             ELSE
                ANGZ=ASIN(RNZI)*180.D0/PI
                IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
                   ANGZ=0.D0
                ELSE
                   ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
                ENDIF
                IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
                   ANGPH=0.D0
                ELSE
                   ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
                ENDIF
             ENDIF
          ELSEIF(MDLWRI.EQ.1) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,PHII,RKR0,ANGZ,ANGPH,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1.D0-SINT2)/(1.D0-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1.D0-SINP2)/(1.D0-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
          ELSEIF(MDLWRI.EQ.2) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,PHI,MODEW,ANGZ,ANGPH,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,PHII,MODEW,ANGZ,ANGPH,UUI
             SINP2=SIN(ANGZ*PI/180.D0)**2
             SINT2=SIN(ANGPH*PI/180.D0)**2
             RNZI  =SQRT(SINP2*(1-SINT2)/(1-SINP2*SINT2))
             RNPHII=SQRT(SINT2*(1-SINP2)/(1-SINP2*SINT2))
             IF(ANGZ.LT.0.D0) RNZI=-RNZI
             IF(ANGPH.LT.0.D0) RNPHII=-RNPHII
             WRITE(6,*) 'XX MDLWRI=2 IS NOT SUPPORTED YET.'
             GOTO 1
          ELSEIF(MDLWRI.EQ.11) THEN
             READ(5,*,ERR=1,END=9000) &
                  RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
             WRITE(6,*) &
                  '# initial values: RF,RP,ZP,RKR0,RNZ,RNPHI,UU'
             WRITE(6,'(1PE12.4,0P7F9.2)') &
                  RF,RPI,ZPI,RKR0,RNZI,RNPHII,UUI
             IF(ABS(RNZI).GT.1.D0) THEN
                ANGZ=0.D0
                ANGPH=0.D0
             ELSE
                ANGZ=ASIN(RNZI)*180.D0/PI
                IF(ABS(RNZI).GT.SQRT(1.D0-RNPHII**2)) THEN
                   ANGZ=0.D0
                ELSE
                   ANGZ=ASIN(RNZI/SQRT(1.D0-RNPHII**2))*180.D0/PI
                ENDIF
                IF(ABS(RNPHII).GT.SQRT(1.D0-RNZI**2)) THEN
                   ANGPH=0.D0
                ELSE
                   ANGPH=ASIN(RNPHII/SQRT(1.D0-RNZI**2))*180.D0/PI
                ENDIF
             ENDIF
          ENDIF
          RFIN(NRAY)=RF
          RPIN(NRAY)=RPI
          ZPIN(NRAY)=ZPI
          PHIIN(NRAY)=PHII
          RKRIN(NRAY)=RKR0
          RNZIN(NRAY)=RNZI
          RNPHIIN(NRAY)=RNPHII
          ANGZIN(NRAY)=ANGZ
          ANGPHIN(NRAY)=ANGPH
          MODEWIN(NRAY)=MODEW
          UUIN(NRAY)=UUI
       ENDIF
    END DO
    RETURN

9000 CONTINUE
    IERR=10
    RETURN
  END SUBROUTINE WRSETUP

  ! --- calculation of radial profile of power absorption ---

  SUBROUTINE wr_calc_pwr

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nrs,nray,nstp,nrs1,nrs2,ndrs,locmax
    REAL(rkind):: drho,xl,yl,zl,rhon1,rhon2,sdr,delpwr,pwrmax,dpwr,ddpwr
    INTEGER,SAVE:: nrsmax_save=0,nrrmax_save=0,nraymax_save=0,nstpmax_save=0

    !   --- allocate variables for power deposition profile ---

    IF(nstpmax.NE.nstpmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(rs_nstp_nray)) DEALLOCATE(rs_nstp_nray)
       IF(ALLOCATED(rr_nstp_nray)) DEALLOCATE(rr_nstp_nray)
       ALLOCATE(rs_nstp_nray(nstpmax,nraymax),rr_nstp_nray(nstpmax,nraymax))
    END IF
    IF(nrsmax.NE.nrsmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_nrs)) DEALLOCATE(pos_nrs)
       IF(ALLOCATED(pwr_nrs)) DEALLOCATE(pwr_nrs)
       IF(ALLOCATED(pwr_nrs_nray)) DEALLOCATE(pwr_nrs_nray)
       ALLOCATE(pos_nrs(nrsmax),pwr_nrs(nrsmax),pwr_nrs_nray(nrsmax,nraymax))
    END IF
    IF(nrrmax.NE.nrrmax_save.OR.nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_nrr)) DEALLOCATE(pos_nrr)
       IF(ALLOCATED(pwr_nrr)) DEALLOCATE(pwr_nrr)
       IF(ALLOCATED(pwr_nrr_nray)) DEALLOCATE(pwr_nrr_nray)
       ALLOCATE(pos_nrr(nrrmax),pwr_nrr(nrrmax),pwr_nrr_nray(nrrmax,nraymax))
    END IF
    IF(nraymax.NE.nraymax_save) THEN
       IF(ALLOCATED(pos_pwrmax_rs_nray)) DEALLOCATE(pos_pwrmax_rs_nray)
       IF(ALLOCATED(pwrmax_rs_nray)) DEALLOCATE(pwrmax_rs_nray)
       IF(ALLOCATED(pos_pwrmax_rr_nray)) DEALLOCATE(pos_pwrmax_rr_nray)
       IF(ALLOCATED(pwrmax_rr_nray)) DEALLOCATE(pwrmax_rr_nray)
       ALLOCATE(pos_pwrmax_rs_nray(nraymax),pwrmax_rs_nray(nraymax))
       ALLOCATE(pos_pwrmax_rr_nray(nraymax),pwrmax_rr_nray(nraymax))
    END IF
    nrsmax_save=nrsmax
    nrrmax_save=nrrmax
    nraymax_save=nraymax
    nstpmax_save=nstpmax

!     ----- CALCULATE RADIAL DEPOSITION PROFILE (Minor radius) -----

    drho=1.D0/nrsmax
    DO nrs=1,nrsmax
       pos_nrs(nrs)=(DBLE(nrs)-0.5D0)*drho
    ENDDO
    DO nray=1,nraymax
       DO nrs=1,nrsmax
          pwr_nrs_nray(nrs,nray)=0.D0
       ENDDO
    ENDDO
    DO nrs=1,nrsmax
       pwr_nrs(nrs)=0.D0
    ENDDO

    DO nray=1,nraymax
       DO nstp=0,nstpmax_nray(nray)-1
          xl=rays(1,nstp,nray)
          yl=rays(2,nstp,nray)
          zl=rays(3,nstp,nray)
          CALL pl_mag_old(xl,yl,zl,rhon1)
          IF(rhon1.LE.1.D0) THEN
             nrs1=INT(rhon1/drho)+1
             IF(nrs1.GT.nrsmax) nrs1=nrsmax
             xl=rays(1,nstp+1,nray)
             yl=rays(2,nstp+1,nray)
             zl=rays(3,nstp+1,nray)
             CALL PL_MAG_OLD(xl,yl,zl,rhon2)
             nrs2=int(rhon2/drho)+1
             IF(nrs2.GT.nrsmax) nrs2=nrsmax
             ndrs=ABS(nrs2-nrs1)
             IF(ndrs.EQ.0) THEN
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +rays(8,nstp+1,nray)
             ELSE IF(nrs1.lt.nrs2) THEN
                sdr=(rhon2-rhon1)/drho
                delpwr=rays(8,nstp+1,nray)/sdr
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +(DBLE(nrs1)-rhon1/drho)*delpwr
                DO nrs=nrs1+1,nrs2-1
                   pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray)+delpwr
                ENDDO
                pwr_nrs_nray(nrs2,nray)=pwr_nrs_nray(nrs2,nray) &
                     +(rhon2/drho-DBLE(nrs2-1))*delpwr
             ELSE
                sdr=(rhon1-rhon2)/drho
                delpwr=rays(8,nstp+1,nray)/sdr
                pwr_nrs_nray(nrs2,nray)=pwr_nrs_nray(nrs2,nray) &
                     +(DBLE(nrs2)-rhon2/drho)*delpwr
                DO nrs=nrs2+1,nrs1-1
                   pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray)+delpwr
                ENDDO
                pwr_nrs_nray(nrs1,nray)=pwr_nrs_nray(nrs1,nray) &
                     +(rhon1/drho-DBLE(nrs1-1))*delpwr
             END IF
          ENDIF
       ENDDO

       DO nrs=1,nrsmax
          pwr_nrs_nray(nrs,nray)=pwr_nrs_nray(nrs,nray) &
               /(2.D0*PI*(DBLE(nrs)-0.5D0)*drho*drho)
       ENDDO
    ENDDO

!     ----- sum of power for each ray-------

    DO nrs=1,nrsmax
       pwr_nrs(nrs)=0.D0
       DO nray=1,nraymax
          pwr_nrs(nrs)=pwr_nrs(nrs)+pwr_nrs_nray(nrs,nray)
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
                -pwr_nrs_nray(locmax-1,nray))/(2.D0*drho)
          ddpwr=(pwr_nrs_nray(locmax+1,nray) &
              -2*pwr_nrs_nray(locmax  ,nray) &
                +pwr_nrs_nray(locmax-1,nray))/drho**2
          pos_pwrmax_rs_nray(nray)=(locmax-0.5D0)/(nrsmax-1.D0)-dpwr/ddpwr
          pwrmax_rs_nray(nray)=pwrmax-dpwr**2/(2.D0*ddpwr)
       ENDIF
    END DO
   
    pwrmax=0.D0
    locmax=0
    DO nrs=1,nrsmax
       IF(pwr_nrs(nrs).GT.pwrmax) THEN
          pwrmax=pwr_nrs(nrs)
          locmax=nrs
       ENDIF
    END DO
    IF(locmax.LE.1) THEN
       pos_pwrmax_rs=pos_nrs(1)
    ELSE IF(locmax.GE.nrsmax) THEN
       pos_pwrmax_rs=pos_nrs(nrsmax)
    ELSE
       dpwr =(pwr_nrs(locmax+1) &
             -pwr_nrs(locmax-1))/(2.D0*drho)
       ddpwr=(pwr_nrs(locmax+1) &
           -2*pwr_nrs(locmax  ) &
             +pwr_nrs(locmax-1))/drho**2
       pos_pwrmax_rs=(locmax-0.5D0)/(nrsmax-1.D0)-dpwr/ddpwr
       pwrmax_rs=pwrmax-dpwr**2/(2.D0*ddpwr)
    ENDIF
   
    WRITE(6,'(A,1PE12.4,A,1PE12.4)') &
         '    PWRMAX=',pwrmax_rs,'  AT RHON =',pos_pwrmax_rs

    CALL PAGES
    CALL grd1d(0,pos_nrs,pwr_nrs_nray,nrsmax,nrsmax,nraymax, &
         '@pwr_nrs vs. pos_nrs@')
    CALL PAGEE

    RETURN
  END SUBROUTINE wr_calc_pwr

END MODULE wrcalc
