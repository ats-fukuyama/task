MODULE DPTENS

CONTAINS

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE DP_TENS(CW,CKPR,CKPP,NS,mag,plf,grd,CLDISP)

    USE dpcomm
    USE dptnsr1
    USE dpfpin
    USE dphotf
    USE dphotr
    USE plprof
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    TYPE(pl_mag_type),INTENT(IN):: mag
    TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
    TYPE(pl_grd_type),DIMENSION(nsmax),INTENT(IN):: grd
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
    COMPLEX(rkind):: CLDISP1(6)
    INTEGER:: ID1,ID2,IDV,IERR

      IF(plf(NS)%RN.LE.0.D0) THEN
         CALL DPTNCL(CW,NS,mag,plf,CLDISP)
!      ELSEIF(RHON_LOC.LT.RHON_MIN.OR.RHON_LOC.GT.RHON_MAX) THEN
!         ID1=MOD(MODELP(NS),100)
!         CALL DPTENS_AN(ID1,CW,CKPR,CKPP,NS,mag,plf,grd,CLDISP)
      ELSE
         ID1=MOD(MODELP(NS),100)
         ID2=MODELP(NS)/100
         IDV=MODELV(NS)
         SELECT CASE(IDV)
         CASE(0)
            CALL DPTENS_AN(ID1,CW,CKPR,CKPP,NS,mag,plf,grd,CLDISP)
         CASE(1)
            CALL DPFMFL(NS,plf,0)
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DP_HOTFI(CW,CKPR,CKPP,NS,mag,CLDISP)
            ELSE
               CALL DP_HOTF(CW,CKPR,CKPP,NS,mag,CLDISP)
            ENDIF
         CASE(2)
            CALL DPFPFL(NS,mag,IERR)
            IF(IERR.NE.0) RETURN
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DP_HOTFI(CW,CKPR,CKPP,NS,mag,CLDISP)
            ELSE
               CALL DP_HOTF(CW,CKPR,CKPP,NS,mag,CLDISP)
            ENDIF
         CASE(3)
            CALL DPFMFL(NS,plf,1)
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DP_HOTRI(CW,CKPR,CKPP,NS,mag,CLDISP)
            ELSE
               CALL DP_HOTR(CW,CKPR,CKPP,NS,mag,CLDISP)
            ENDIF
         CASE(4)
            CALL DPFPFL(NS,mag,IERR)
            IF(IERR.NE.0) RETURN
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DP_HOTRI(CW,CKPR,CKPP,NS,mag,CLDISP)
            ELSE
               CALL DP_HOTR(CW,CKPR,CKPP,NS,mag,CLDISP)
            ENDIF
         END SELECT

         IF(ID2.EQ.2) THEN
            CALL DPTNCL(CW,NS,mag,plf,CLDISP1)
         ELSEIF(ID2.EQ.3) THEN
            CALL DPTNKP(CW,CKPR,CKPP,NS,mag,plf,CLDISP1)
         ENDIF

         IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
            CLDISP(1)=DBLE(CLDISP1(1))+CI*DIMAG(CLDISP(1))
            CLDISP(2)=DBLE(CLDISP1(2))+CI*DIMAG(CLDISP(2))
            CLDISP(3)=DBLE(CLDISP1(3))+CI*DIMAG(CLDISP(3))
            CLDISP(4)=DBLE(CLDISP1(4))+CI*DIMAG(CLDISP(4))
            CLDISP(5)=CI*DIMAG(CLDISP1(5))+DBLE(CLDISP(5))
            CLDISP(6)=CI*DIMAG(CLDISP1(6))+DBLE(CLDISP(6))
         ENDIF
      ENDIF
!         WRITE(6,'(A,3I5,1PE12.4)') &
!              'NS,ID1,ID2,RN =',NS,ID1,ID2,RN(NS)
!         WRITE(6,'(A,1P6E12.4)') &
!              'CW,R,P=',CW,CKPR,CKRR
!         WRITE(6,'(A,1P6E12.4)') &
!              'CLDISP=',CLDISP(1),CLDISP(2),CLDISP(3)
!         WRITE(6,'(A,1P6E12.4)') &
!              'CLDISP=',CLDISP(4),CLDISP(5),CLDISP(6)
      RETURN
  END SUBROUTINE DP_TENS

!     ****** CALCULATE ANALYTIC DIELECTRIC TENSOR ******

  SUBROUTINE DPTENS_AN(ID1,CW,CKPR,CKPP,NS,mag,plf,grd,CLDISP)

    USE dpcomm
    USE dptnsr1
    USE dptnsr2
    USE plprof
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: ID1,NS
    TYPE(pl_mag_type),INTENT(IN):: mag
    TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
    TYPE(pl_grd_type),DIMENSION(nsmax),INTENT(IN):: grd
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)

      SELECT CASE(ID1)
      CASE(0) ! collisionless cold model
         CALL DPTNCL(CW,NS,mag,plf,CLDISP)
      CASE(1) ! collisional cold model
         CALL DPTNCC(CW,NS,mag,plf,CLDISP)
      CASE(2) ! idial MHD model
         CALL DPTNIM(CKPR,NS,mag,plf,CLDISP)
      CASE(3) ! resistive MHD model
         CALL DPTNRM(CW,NS,mag,plf,CLDISP)
      CASE(4) ! kinetic without FLR model
         CALL DPTNHP(CW,CKPR,CKPP,NS,mag,plf,CLDISP)
      CASE(5) ! kinetic with FLR model
         CALL DPTNKS(CW,CKPR,CKPP,NS,mag,plf,CLDISP)
      CASE(6) ! kinetic with FLR model
         CALL DPTNKP(CW,CKPR,CKPP,NS,mag,plf,CLDISP)
      CASE(7) ! weakly relativistic model
         CALL DPTNKR(CW,CKPR,CKPP,NS,mag,plf,CLDISP)
      CASE(11:15) ! old WM models model
         CALL DPTNFK2(CW,CKPR,CKPP,NS,mag,plf,CLDISP)
      CASE(16) ! old WM drift kinetic model
         CALL DPTNDK0(CW,CKPR,CKPP,NS,mag,plf,grd,CLDISP)
      CASE DEFAULT
         WRITE(6,*) 'XX WRONG MODELP IN DPTENS: ID1=',ID1
      END SELECT
      RETURN
  END SUBROUTINE DPTENS_AN
END MODULE DPTENS
