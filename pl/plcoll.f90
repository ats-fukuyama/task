! plcoll.f90

MODULE plcoll

  PRIVATE
  PUBLIC pl_set_rnuc
  PUBLIC pl_coll

CONTAINS

  SUBROUTINE pl_set_rnuc(plf)
    USE plcomm
    USE plcomm_type
    IMPLICIT NONE
    TYPE(pl_prf_type),DIMENSION(nsmax):: plf
    REAL(rkind),DIMENSION(nsmax) :: rn,rtpp,rtpr,rnu
    INTEGER:: ns

    DO ns=1,nsmax
       rn(ns)=plf(ns)%rn
       rtpp(ns)=plf(ns)%rtpp
       rtpr(ns)=plf(ns)%rtpr
    END DO
    CALL pl_coll(rn,rtpr,rtpp,rnu)
    DO ns=1,nsmax
       plf(ns)%rnuc=pnuc(ns)*rnu(ns)
    END DO
    RETURN
  END SUBROUTINE pl_set_rnuc

  SUBROUTINE pl_coll(rn,rtpr,rtpp,rnu)
    USE plcomm
    USE plcomm_type
    IMPLICIT NONE
    REAL(rkind),INTENT(IN),DIMENSION(nsmax):: rn,rtpr,rtpp
    REAL(rkind),INTENT(OUT),DIMENSION(nsmax):: rnu
    REAL(rkind):: RNE,RNTE,TE,TI,RNTI,RNZI,RLAMEE,RLAMEI,RLAMII,SN,PNN0
    REAL(rkind):: VTE,RNUEE,RNUEI,RNUEN,RNUE
    REAL(rkind):: VTI,RNUIE,RNUII,RNUIN,RNUI
    INTEGER:: ns

    SN=1.D-20 ! typical ionizatioin crosssection
    PNN0=PPN0/(PTN0*AEE) ! neutral density

    ! --- set collision frequency ---

    DO ns=1,nsmax
       RNE=0.D0
       RNTE=0.D0
       RNTI=0.D0
       RNZI=0.D0
       IF(PA(ns).LE.0.1D0) THEN  ! electron or positron
          RNE=RNE+RN(ns)
          RNTE=RNTE+RN(ns)*(RTPR(ns)+2.D0*RTPP(ns))*1.D3/3.D0   ! eV unit
       ELSE ! ion
          RNTI=RNTI+RN(ns)*(RTPR(ns)+2.D0*RTPP(ns))*1.D3/3.D0
          RNZI=RNZI+RN(ns)*PZ(ns)**2
       END IF
    END DO

    IF(RNE.GT.0.D0) THEN
       TE=RNTE/RNE
       TI=RNTI/RNE
       IF(TI.LE.0.D0) TI=0.03D0
       RLAMEE= 8.0D0+2.3D0*(LOG10(TE)-0.5D0*LOG10(RNE))
       RLAMEI= RLAMEE+0.3D0
       RLAMII=12.1D0+2.3D0*(LOG10(TI)-0.5D0*LOG10(RNE))

       DO NS=1,NSMAX
          IF(PA(NS).LE.0.1D0) THEN ! electron or positron
             VTE=SQRT(2.D0*TE*AEE/AME)
             RNUEE=RNE*RLAMEE/(1.24D-4*SQRT(TE*1.D-3)**3)
             RNUEI=RNZI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
             RNUEN=PNN0*SN*0.88D0*VTE
             SELECT CASE(model_coll)
             CASE(1)
                RNUE=RNUEE+RNUEI
             CASE DEFAULT
                RNUE=RNUEE+RNUEI+RNUEN
             END SELECT
             RNU(NS)=RNUE
          ELSE
             TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
             VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
             RNUIE=PZ(NS)**2*RNE*RLAMEI &
                  /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
             RNUII=RNZI*RLAMII &
                  /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
             RNUIN=PNN0*SN*0.88D0*VTI
             SELECT CASE(model_coll)
             CASE(1)
                RNUI=RNUIE+RNUII
             CASE DEFAULT
                RNUI=RNUIE+RNUII+RNUIN
             END SELECT
             RNU(NS)=RNUI
          ENDIF
       ENDDO
    ELSE
       DO NS=1,NSMAX
          RNU(NS)=0.D0
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE pl_coll

! ****** CALCULATE COLLISION FREQUENCY ******

  SUBROUTINE pl_coll2(rn,rtpr,rtpp,rnue,rnui,rnun)

    USE plcomm
    USE plcomm_type
    USE platmd
    IMPLICIT NONE
    REAL(rkind),INTENT(IN),DIMENSION(nsmax):: rn,rtpr,rtpp
    REAL(rkind),INTENT(OUT),DIMENSION(nsmax):: rnue,rnui,rnun
    REAL(rkind),PARAMETER:: RKEV=1.D3
    REAL(rkind):: RNSUM,RNTSUM,RTAVE,TE,VTE,SNE,SVE,SNI,PNN0,TI,VTI
    REAL(rkind):: RLAMEE,RLAMEI,RLAMII,RNUEE,RNUEI,RNUEN,RNUIE,RNUII,RNUIN
    REAL(rkind):: RNESUM,RNTESUM
    INTEGER:: NS,NSI

    RNESUM=0.D0
    RNTESUM=0.D0
    RNSUM=0.D0
    RNTSUM=0.D0
    DO NS=1,NSMAX
       IF(PA(NS).LE.0.1D0) THEN
          RNESUM=RNESUM+RN(NS)
          RNTESUM=RNTESUM+RN(NS)*(RTPR(NS)+2.D0*RTPP(NS))*RKEV/3.D0
       END IF
       RNSUM=RNSUM+RN(NS)
       RNTSUM=RNTSUM+RN(NS)*(RTPR(NS)+2.D0*RTPP(NS))*RKEV/3.D0
    END DO
    IF(RNSUM.GT.0.D0) THEN
       RTAVE=RNTSUM/RNSUM
    ELSE
       RTAVE=0.03D0
    END IF
         
    DO NS=1,NSMAX
       IF(RNESUM.GT.0.AND.RTAVE.GT.0.D0) THEN
          IF(PA(NS).LE.0.1D0) THEN ! electron or positron
             IF(RTAVE.GT.6.65D0) THEN
                RLAMEE=8.0D0 &
                      +2.3D0*(LOG10(RTAVE)-0.5D0*LOG10(RNESUM))
             ELSE
                RLAMEE=7.0D0 &
                      +2.3D0*(1.5D0*LOG10(RTAVE)-0.5D0*LOG10(RNESUM))
             END IF
          ELSE
             RLAMEE=1.D0
          END IF
          IF(RTAVE.GT.13.3D0) THEN
             RLAMEI=8.3D0 &
                   +2.3D0*(LOG10(RTAVE)-0.5D0*LOG10(RNESUM))
          ELSE
             RLAMEI=7.0D0 &
                   +2.3D0*(1.5D0*LOG10(RTAVE)-0.5D0*LOG10(RNESUM))
          END IF
          IF(NS.GT.1) THEN
             IF(RTAVE.GT.24.5D3) THEN
                RLAMII=12.1D0 &
                      +2.3D0*(LOG10(RTAVE)-0.5D0*LOG10(RNESUM))
             ELSE
                RLAMII=7.0D0 &
                      +2.3D0*(1.5D0*LOG10(RTAVE)-0.5D0*LOG10(RNESUM))
             END IF
          ELSE
             RLAMII=1.D0
          ENDIF
       ELSE
          RLAMEE= 1.D0
          RLAMEI= 1.D0
          RLAMII= 1.D0
       END IF
    END DO

    IF(MODELN.EQ.0) THEN
       TE=RNTESUM/RNESUM
       VTE=SQRT(2.D0*TE*AEE/AME)
       SNE=0.88D-20
       SVE=SNE*VTE
    ELSE
       CALL ATSIGV(TE,SVE,1)
    ENDIF
    SNI=1.D-20
    PNN0=PPN0/(PTN0*AEE)

    DO NS=1,NSMAX
       IF(PA(NS).LE.0.1D0) THEN ! electron or positron
          TE=(RTPR(NS)+2.D0*RTPP(NS))*RKEV/3.D0
          VTE=SQRT(2.D0*TE*AEE/AME)
          RNUEE=RN(NS)*RLAMEE &
               /(1.24D-4*SQRT(TE*1.D-3)**3)
          RNUEI=0.D0
          DO NSI=1,NSMAX
             IF(PA(NS).LE.0.1D0) RNUEI=RNUEI+PZ(NSI)**2*RN(NSI)
          ENDDO
          RNUEI=RNUEI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
          RNUEN=PNN0*SVE
          RNUE(NS)=RNUEE
          RNUI(NS)=RNUEI
          RNUN(NS)=RNUEN
       ELSE
          TI=(RTPR(NS)+2.D0*RTPP(NS))*RKEV/3.D0
          VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
          RNUIE=PZ(NS)**2*RN(1)*RLAMEI &
               /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
          RNUII=PZ(NS)**4*RN(NS)*RLAMII &
               /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
          RNUIN=PNN0*SNI*0.88D0*VTI
          RNUE(NS)=RNUIE
          RNUI(NS)=RNUII
          RNUN(NS)=RNUIN
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE pl_coll2

END MODULE plcoll
