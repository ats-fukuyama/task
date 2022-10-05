! plcoll.f90

MODULE plcoll

  PRIVATE
  PUBLIC pl_coll,pl_set_rnuc

CONTAINS

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
             RNUE=RNUEE+RNUEI+RNUEN
             RNU(NS)=RNUE
          ELSE
             TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
             VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
             RNUIE=PZ(NS)**2*RNE*RLAMEI &
                  /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
             RNUII=RNZI*RLAMII &
                  /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
             RNUIN=PNN0*SN*0.88D0*VTI
             RNUI=RNUIE+RNUII+RNUIN
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
END MODULE plcoll
