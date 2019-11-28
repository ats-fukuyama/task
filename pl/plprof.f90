  MODULE pllocal
    USE plcomm,ONLY: rkind,NSM

    REAL(rkind):: XPOS_LOC,YPOS_LOC,ZPOS_LOC,RHON_LOC
    REAL(rkind):: BNX,BNY,BNZ,BABS
    REAL(rkind),DIMENSION(NSM):: RN,RTPR,RTPP,RU,RUPL,RUPR,RUPP,RNUC,RZCL
    REAL(rkind),DIMENSION(NSM):: RLN,RLRPT,RLTPP,RLU,RLUPL
  END MODULE pllocal

  MODULE plprof
    USE bpsd_kinds

    TYPE pl_mag_type
       real(rkind):: BABS,BNX,BNY,BNZ,RHON
                           ! BABS: magnetic field strength [T]
                           ! BNX:  normalized X component of B (phi=0 deg)
                           ! BNY:  normalized Y component of B (phi=90 deg)
                           ! BNZ:  normalized Z component of B (vertical)
                           ! RHON: normalized minor radius
    END TYPE pl_mag_type

    TYPE pl_plf_type       ! local plasma parameter
       real(rkind):: RN,RTPR,RTPP,RU,RUPL,RNUC,RZCL
                           ! RN:   number density [10^{20}m^{-3}]
                           ! RTPR: parallel temperature [keV]
                           ! RTPP: perpendicular temperature [keV]
                           ! RU:   toroidal fluid velocity [m/s]
                           ! RUPL: poloidal fluid velocity [m/s]
                           ! RNUC: collision frequency [1/s]
                           ! RZCL: collision parameter (RNUC/OMEGA)
    END TYPE pl_plf_type

    TYPE pl_grd_type
       real(rkind):: grdn,grdtpr,grdtpp,grdu,grdupl
                           ! GRDNN:  density gradient [1/m]
                           ! GRDTPR: parallel temperature gradient [1/m]
                           ! GRDTPP: perpendicular temperature gradient [1/m]
                           ! GRDU:   toroidal fluid velocity gradient [1/m]
                           ! GRDUPL: poloidal fluid velocity gradient [1/m]
    END TYPE pl_grd_type

    INTERFACE
       SUBROUTINE GETRZ(RL,Z,PP,BR,BZ,BT,RHON)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: RL,Z,PP
         REAL(rkind),INTENT(OUT):: BR,BZ,BT,RHON
       END SUBROUTINE GETRZ
       SUBROUTINE GETQP(RHON,QL)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: RHON
         REAL(rkind),INTENT(OUT):: QL
       END SUBROUTINE GETQP
       SUBROUTINE GETRMN(RHON,RRMINL)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: RHON
         REAL(rkind),INTENT(OUT):: RRMINL
       END SUBROUTINE GETRMN
       SUBROUTINE GETRMX(RHON,RRMAXL)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: RHON
         REAL(rkind),INTENT(OUT):: RRMAXL
       END SUBROUTINE GETRMX
       SUBROUTINE GETAXS(RAXIS,ZAXIS)
         USE bpsd_kinds
         REAL(rkind),INTENT(OUT):: RAXIS,ZAXIS
       END SUBROUTINE GETAXS
       SUBROUTINE GETRSU(RSU,ZSU,NSUM,NSUMAX)
         USE bpsd_kinds
         INTEGER(ikind),INTENT(IN):: NSUM
         REAL(rkind),DIMENSION(NSUM),INTENT(OUT):: RSU,ZSU
         INTEGER(ikind),INTENT(OUT):: NSUMAX
       END SUBROUTINE GETRSU
       SUBROUTINE SPL1DF(X,DATA,XA,UDATA,NXMAX,IERR)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: X
         REAL(rkind),INTENT(OUT):: DATA
         INTEGER(ikind),INTENT(IN):: NXMAX
         REAL(rkind),DIMENSION(NXMAX),INTENT(IN):: XA
         REAL(rkind),DIMENSION(4,NXMAX),INTENT(IN):: UDATA
         INTEGER(ikind),INTENT(OUT):: IERR
       END SUBROUTINE SPL1DF
       SUBROUTINE GET_BMINMAX(rhon,BBMIN,BBMAX)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: rhon
         REAL(rkind),INTENT(OUT):: BBMIN,BBMAX
       END SUBROUTINE GET_BMINMAX
       SUBROUTINE GET_RZ(rhon,chip,R,Z)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: rhon,chip
         REAL(rkind),INTENT(OUT):: R,Z
       END SUBROUTINE GET_RZ
       SUBROUTINE GET_RZB(rhon,chip,R,Z,BR,BZ,BT,BB)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: rhon,chip
         REAL(rkind),INTENT(OUT):: R,Z,BR,BZ,BT,BB
       END SUBROUTINE GET_RZB
       SUBROUTINE GET_B(rhon,chip,BB)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: rhon,chip
         REAL(rkind),INTENT(OUT):: BB
       END SUBROUTINE GET_B
       SUBROUTINE GET_DVDRHO(rhon,dvdrho)
         USE bpsd_kinds
         REAL(rkind),INTENT(IN):: rhon
         REAL(rkind),INTENT(OUT):: dvdrho
       END SUBROUTINE GET_DVDRHO
    END INTERFACE

  CONTAINS
        
!     ****** CALCULATE LOCAL MAGNETIC FIELD (OLD) ******

    SUBROUTINE pl_mag_old(X,Y,Z,RHON)
      USE pllocal, ONLY: XPOS_LOC,YPOS_LOC,ZPOS_LOC,RHON_LOC,BNX,BNY,BNZ,BABS
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,Y,Z
      REAL(rkind),INTENT(OUT):: RHON
      TYPE(pl_mag_type):: MAG
      
      CALL pl_mag(X,Y,Z,MAG)
      RHON=MAG%RHON

      XPOS_LOC=X
      YPOS_LOC=Y
      ZPOS_LOC=Z
      RHON_LOC=RHON
      BNX=MAG%BNX
      BNY=MAG%BNY
      BNZ=MAG%BNZ
      BABS=MAG%BABS

      RETURN
    END SUBROUTINE pl_mag_old

!     ****** CALCULATE LOCAL MAGNETIC FIELD ******

    SUBROUTINE pl_mag(X,Y,Z,MAG)

      USE plcomm, ONLY: RA,RR,BB,MODELG
      USE plprof2d
      USE plload
      IMPLICIT NONE
      REAL(rkind),INTENT(in):: X,Y,Z
      TYPE(pl_mag_type),INTENT(OUT):: MAG
      REAL(rkind):: RHON
      INTEGER:: IERR
      REAL(8):: BABS,AL(3)

      REAL(8) :: RL, BR, BT, BX, BY, BZ, &
                 RSINT, RCOST

      SELECT CASE(MODELG)
      CASE(0,1)
         CALL pl_mag_rz(X,Y,Z,BR,BZ,BT,RHON)
         BX = BR
         BY = BT
      CASE(2,3,5,8)
         RL=SQRT(X**2+Y**2)
         RCOST= X/RL
         RSINT= Y/RL
         CALL pl_mag_rz(X,Y,Z,BR,BZ,BT,RHON)
         BX = BR*RCOST-BT*RSINT
         BY = BR*RSINT+BT*RCOST
      CASE(11)
         CALL PLSMAG11(X,Y,BABS,AL)
         BX = BABS*AL(1)
         BY = BABS*AL(2)
         BZ = BABS*AL(3)
      CASE(12)
         CALL pl_read_p2Dmag(X,Y,BX,BY,BZ,IERR)
      CASE(13)
         CALL PLSMAG13(X,Y,BABS,AL)
         BX = BABS*AL(1)
         BY = BABS*AL(2)
         BZ = BABS*AL(3)
      CASE(14)
         CALL pl_read_p2Dmag(X,Y,BX,BY,BZ,IERR)
      END SELECT

      MAG%BABS = SQRT(BX**2+BY**2+BZ**2)
      MAG%RHON = RHON

      IF(MAG%BABS.LE.0.D0) THEN
         MAG%BNX = 0.D0
         MAG%BNY = 1.D0
         MAG%BNZ = 0.D0
      ELSE
         MAG%BNX = BX/MAG%BABS
         MAG%BNY = BY/MAG%BABS
         MAG%BNZ = BZ/MAG%BABS
      ENDIF
      RETURN
    END SUBROUTINE pl_mag

!     ****** CALCULATE LOCAL MAGNETIC FIELD ******

    SUBROUTINE pl_mag_rz(X,Y,Z,BR,BZ,BT,RHON)

      USE plcomm, ONLY: RA,RR,BB,MODELG
      USE plprof2d
      USE plload
      IMPLICIT NONE
      REAL(rkind),INTENT(in):: X,Y,Z
      REAL(rkind),INTENT(OUT):: BR,BZ,BT,RHON

      INTEGER:: IERR
      REAL(rkind) :: BP, PP, QL, RL, RS, BABS, AL(3), &
                 RSINP, RCOSP

      SELECT CASE(MODELG)
      CASE(0)
         RS   = SQRT(X*X+Z*Z)
         RHON = RS/RA
         CALL pl_qprf(RHON,QL)
         RSINP= Z/RS
         RCOSP= X/RS
         BT   = BB
         IF(RR.EQ.0.D0) THEN
            BP=0.D0
         ELSE
            BP   = RS*BT/(RR*QL)
         END IF
         BR   =-BP*RSINP
         BZ   = BP*RCOSP

      CASE(1)
         RS =SQRT((X-RR)**2+Z*Z)
         RHON=RS/RA
         CALL pl_qprf(RHON,QL)
         RSINP= Z/RS
         RCOSP= (X-RR)/RS
         BT   = BB*RR/X
         BP   = RS*BT/(RR*QL)
         BR   =-BP*RSINP
         BZ   = BP*RCOSP

      CASE(2)
         RL=SQRT(X**2+Y**2)
         RS =SQRT((RL-RR)**2+Z**2)
         RHON=RS/RA
         IF(RS.LE.0.D0) THEN
            BT   = BB
            BR   = 0.D0
            BZ   = 0.D0
         ELSE
            CALL pl_qprf(RHON,QL)
            RSINP= Z/RS
            RCOSP= (RL-RR)/RS
            BT   = BB*RR/RL
            BP   = RS*BT/(RR*QL)
            BR   =-BP*RSINP
            BZ   = BP*RCOSP
         ENDIF

      CASE(3,5,8)
         RL=SQRT(X**2+Y**2)
         PP=0.D0
         CALL GETRZ(RL,Z,PP,BR,BZ,BT,RHON)

      CASE(11)
         CALL plsmag11(X,Y,BABS,AL)
         BR=BABS*AL(1)
         BZ=BABS*AL(2)
         BT=BABS*AL(3)

      CASE(12)
         CALL pl_read_p2Dmag(X,Y,BR,BZ,BT,IERR)

      CASE(13)
         CALL plsmag13(X,Y,BABS,AL)
         BR=BABS*AL(1)
         BZ=BABS*AL(2)
         BT=BABS*AL(3)

      CASE(14)
         CALL pl_read_p2Dmag(X,Y,BR,BZ,BT,IERR)

      END SELECT

      RETURN
    END SUBROUTINE pl_mag_rz

!     ****** CALCULATE PLASMA PROFILE ******

    SUBROUTINE pl_prof2(RHON,RN,RTPR,RTPP,RU)

      USE plcomm
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: RHON
      REAL(rkind),INTENT(OUT),DIMENSION(NSMAX):: RN,RTPR,RTPP,RU
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      INTEGER(ikind):: NS

      CALL pl_prof(RHON,PLF)
      DO NS=1,NSMAX
         RN(NS)  =PLF(NS)%RN
         RTPR(NS)=PLF(NS)%RTPR
         RTPP(NS)=PLF(NS)%RTPP
         RU(NS)  =PLF(NS)%RU
      ENDDO
      RETURN
    END SUBROUTINE pl_prof2

!     ****** CALCULATE PLASMA PROFILE ******

    SUBROUTINE pl_prof3d_old(X,Y,Z)

      USE plcomm,ONLY: NSMAX
      USE pllocal,ONLY: RN,RTPR,RTPP,RU,RUPL,RNUC,RZCL
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,Y,Z
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      INTEGER:: NS

      CALL pl_prof3d(X,Y,Z,PLF)
      DO NS=1,NSMAX
         RN(NS)=PLF(NS)%RN
         RTPR(NS)=PLF(NS)%RTPR
         RTPP(NS)=PLF(NS)%RTPP
         RU(NS)=PLF(NS)%RU
         RUPL(NS)=PLF(NS)%RUPL
         RNUC(NS)=PLF(NS)%RNUC
         RZCL(NS)=PLF(NS)%RZCL
      END DO
      RETURN
    END SUBROUTINE pl_prof3d_old

!     ****** CALCULATE PLASMA PROFILE ******

    SUBROUTINE pl_prof_old(RHON)

      USE plcomm,ONLY: NSMAX
      USE pllocal,ONLY: RN,RTPR,RTPP,RU,RUPL,RNUC,RZCL
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: RHON
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      INTEGER:: NS

      CALL pl_prof(RHON,PLF)
      DO NS=1,NSMAX
         RN(NS)=PLF(NS)%RN
         RTPR(NS)=PLF(NS)%RTPR
         RTPP(NS)=PLF(NS)%RTPP
         RU(NS)=PLF(NS)%RU
         RUPL(NS)=PLF(NS)%RUPL
         RNUC(NS)=PLF(NS)%RNUC
         RZCL(NS)=PLF(NS)%RZCL
      END DO
      RETURN
    END SUBROUTINE pl_prof_old

!     ****** CALCULATE PLASMA PROFILE ******

!     Input: RHON : 
!
!     Output: PLF(NS)%RN   : Density
!             PLF(NS)%RTPR : Parallel temperature
!             PLF(NS)%RTPP : Perpendicular temperature
!             PLF(NS)%RU   : Toroidal rotation velocity
!             PLF(NS)%RUPL : Poloidal rotation velocity
!             PLF(NS)%RNUC : collision frequency
!             PLF(NS)%RZCL : collision factor (=nu/omega) 


!     ****** CALCULATE PLASMA PROFILE in 3D ******

!     Density, temperatures, rotation and Ratio of collision
!       frequency to wave frequency evaluated at a given point, RHON
!
!     Input: X,Y,Z : X and Y : horizontal plane
!                    Z : vertical axis
!
!     Output: PLF(NS)%RN   : Density
!             PLF(NS)%RTPR : Parallel temperature
!             PLF(NS)%RTPP : Perpendicular temperature
!             PLF(NS)%RU   : Toroidal rotation velocity
!             PLF(NS)%RUPL : Poloidal rotation velocity

    SUBROUTINE pl_prof3d(X,Y,Z,PLF)

        USE plcomm,ONLY: PZ,PN,PTPR,PTPP,PU,PNS,PTS,PUS,PZCL, &
             NSMAX,MODELN,MODELG, &
             RR,RA
        USE plprof2d
        USE plload
        IMPLICIT NONE
        REAL(rkind),INTENT(IN):: X,Y,Z
        TYPE(pl_plf_type),DIMENSION(NSMAX),INTENT(OUT):: PLF
        REAL(rkind),DIMENSION(NSMAX) :: &
             RN_PL,RT_PL,RTPR_PL,RTPP_PL,RU_PL,RNUC_PL,RZCL_PL
        REAL(rkind):: RHON,FACTX,FACTY,FACTN,FACTT,FACTU
        INTEGER:: NS,NSMAXL,IERR


        SELECT CASE(MODELG)
        CASE(0)
           IF(RA.EQ.0.D0) THEN
              FACTX=1.D0
           ELSE
              FACTX=EXP(-(X*X+Z*Z)/(RA*RA))
           ENDIF
           IF(RR.EQ.0.D0) THEN
              FACTY=1.D0
           ELSE
              FACTY=EXP(-Y/RR)
           ENDIF
           FACTN=FACTX*FACTY
           FACTT=FACTX*FACTY
           FACTU=FACTX*FACTY
           DO NS=1,NSMAX
              PLF(NS)%RN  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
              PLF(NS)%RTPR=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
              PLF(NS)%RTPP=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
              PLF(NS)%RU  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
              PLF(NS)%RUPL=0.D0
              PLF(NS)%RNUC=0.D0
              PLF(NS)%RZCL=PZCL(NS)
           END DO
        CASE(1)
           IF(RA.EQ.0.D0) THEN
              FACTX=1.D0
           ELSE
              FACTX=EXP(-(X*X+Z*Z)/(RA*RA))
           ENDIF
           IF(RR.EQ.0.D0) THEN
              FACTY=1.D0
           ELSE
              FACTY=EXP(-Y/RR)
           ENDIF
           FACTN=FACTX*FACTY
           FACTT=FACTX*FACTY
           FACTU=FACTX*FACTY
           DO NS=1,NSMAX
              PLF(NS)%RN  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
              PLF(NS)%RTPR=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
              PLF(NS)%RTPP=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
              PLF(NS)%RU  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
              PLF(NS)%RUPL=0.D0
              PLF(NS)%RNUC=0.D0
              PLF(NS)%RZCL=PZCL(NS)
           END DO
        CASE(11)
           CALL PLSDEN11(X,Y,RN_PL,RTPR_PL,RTPP_PL,RZCL_PL)
           DO NS=1,NSMAX
              PLF(NS)%RN  =RN_PL(NS)
              PLF(NS)%RTPR=RTPR_PL(NS)
              PLF(NS)%RTPP=RTPP_PL(NS)
              PLF(NS)%RU  =0.D0
              PLF(NS)%RUPL=0.D0
              PLF(NS)%RNUC=0.D0
              PLF(NS)%RZCL=PZCL(NS)
           ENDDO
        CASE(12)
           CALL pl_read_p2D(X,Y,RN_PL,RT_PL,RU_PL,NSMAXL,IERR)
           DO NS=1,NSMAXL
              PLF(NS)%RN  =RN_PL(NS)
              PLF(NS)%RTPR=RT_PL(NS)
              PLF(NS)%RTPP=RT_PL(NS)
              PLF(NS)%RU  =RU_PL(NS)
              PLF(NS)%RUPL=0.D0
              PLF(NS)%RNUC=0.D0
              PLF(NS)%RZCL=PZCL(NS)
           ENDDO
        CASE(13)
           CALL PLSDEN13(X,Y,RN_PL,RTPR_PL,RTPP_PL,RZCL_PL)
           DO NS=1,NSMAX
              PLF(NS)%RN  =RN_PL(NS)
              PLF(NS)%RTPR=RTPR_PL(NS)
              PLF(NS)%RTPP=RTPP_PL(NS)
              PLF(NS)%RU  =0.D0
              PLF(NS)%RUPL=0.D0
              PLF(NS)%RNUC=0.D0
              PLF(NS)%RZCL=PZCL(NS)
           ENDDO
        CASE DEFAULT
           CALL pl_mag_old(X,Y,Z,RHON)
           CALL pl_prof(RHON,PLF)
        END SELECT

        RETURN
      END SUBROUTINE pl_prof3d
     

!     ****** CALCULATE PLASMA PROFILE ******

    SUBROUTINE pl_prof(RHON,PLF)

        USE plcomm,ONLY: PZ,PN,PTPR,PTPP,PU,PNS,PTS,PUS,PZCL, NSMAX&
             &,MODELN,RA,RB, PROFN1,PROFN2,PROFT1,PROFT2, PROFU1&
             &,PROFU2, PNITB,PTITB,PUITB,RHOITB,RHOEDG
        USE plload,ONLY: pl_read_trdata
        IMPLICIT NONE
        REAL(rkind),INTENT(IN):: RHON
        TYPE(pl_plf_type),DIMENSION(NSMAX),INTENT(OUT):: PLF
        REAL(rkind):: RHOL, FACTN, FACTT, FACTU, FACTITB, PL0, PL,&
             & FACT, FNX, DFNX, AN, BN, FTX, DFTX, AT, BT, FUX, DFUX,&
             & AU, BU, VAL, PNL, PTL
        INTEGER(ikind)  :: NS
        REAL(rkind),DIMENSION(NSMAX) :: RN_PL,RT_PL,RTPR_PL,RTPP_PL&
             &,RU_PL,RUPL_PL

        IF(RHON.LE.0.D0) THEN
           RHOL=0.D0
! Modified by AF on 2014-10-11
!        ELSEIF(RHON.GE.1.D0) THEN
        ELSEIF(RHON.GE.RB/RA) THEN
           RHOL=RB/RA
        ELSE
           RHOL=RHON
        ENDIF

        SELECT CASE(MODELN)
        CASE(0,1)
           IF(RHOL.GT.1.D0) THEN
              DO NS=1,NSMAX
                 IF(MODELN.EQ.1) THEN
                    PLF(NS)%RN  =PNS(NS)
                 ELSE
                    PLF(NS)%RN  =0.D0
                 END IF
                 PLF(NS)%RTPR=PTS(NS)
                 PLF(NS)%RTPP=PTS(NS)
                 PLF(NS)%RU  =PUS(NS)
                 PLF(NS)%RUPL=0.D0
                 PLF(NS)%RNUC=0.D0
                 PLF(NS)%RZCL=PZCL(NS)
              ENDDO
           ELSE
              DO NS=1,NSMAX
                 FACTN=(1.D0-RHOL**PROFN1(NS))**PROFN2(NS)
                 FACTT=(1.D0-RHOL**PROFT1(NS))**PROFT2(NS)
                 FACTU=(1.D0-RHOL**PROFU1(NS))**PROFU2(NS)
                 PLF(NS)%RN  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
                 PLF(NS)%RTPR=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
                 PLF(NS)%RTPP=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
                 PLF(NS)%RU  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
                 PLF(NS)%RUPL=0.D0
                 PLF(NS)%RNUC=0.D0
                 PLF(NS)%RZCL=PZCL(NS)
                 IF(RHOL.LT.RHOITB(NS)) THEN
                    FACTITB =(1.D0-(RHOL/RHOITB(NS))**4)**2
                    PLF(NS)%RN  =PLF(NS)%RN  +PNITB(NS)*FACTITB
                    PLF(NS)%RTPR=PLF(NS)%RTPR+PTITB(NS)*FACTITB
                    PLF(NS)%RTPP=PLF(NS)%RTPP+PTITB(NS)*FACTITB
                    PLF(NS)%RU  =PLF(NS)%RU  +PUITB(NS)*FACTITB
                    PLF(NS)%RUPL=0.D0
                    PLF(NS)%RNUC=0.D0
                    PLF(NS)%RZCL=PZCL(NS)
               ENDIF
            ENDDO
         ENDIF

      CASE(2)
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               PLF(NS)%RN  =PNS(NS)
               PLF(NS)%RTPR=PTS(NS)
               PLF(NS)%RTPP=PTS(NS)
               PLF(NS)%RU  =PUS(NS)
               PLF(NS)%RUPL=0.D0
               PLF(NS)%RNUC=0.D0
               PLF(NS)%RZCL=PZCL(NS)
            ENDDO
         ELSE
            CALL GETPP(0.D0,PL0)
            CALL GETPP(RHOL,PL)
            FACT=SQRT(PL/PL0)
            DO NS=1,NSMAX
               FACTU=(1.D0-RHOL**PROFU1(NS))**PROFU2(NS)
               PLF(NS)%RN  =(PN(NS)-PNS(NS))*FACT+PNS(NS)
               PLF(NS)%RTPR=(PTPR(NS)-PTS(NS))*FACT+PTS(NS)
               PLF(NS)%RTPP=(PTPP(NS)-PTS(NS))*FACT+PTS(NS)
               PLF(NS)%RU  =(PU(NS)-PUS(NS))*FACTU+PUS(NS)
               PLF(NS)%RUPL=0.D0
               PLF(NS)%RNUC=0.D0
               PLF(NS)%RZCL=PZCL(NS)
            ENDDO
         ENDIF

      CASE(3)
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               PLF(NS)%RN  =PNS(NS)
               PLF(NS)%RTPR=PTS(NS)
               PLF(NS)%RTPP=PTS(NS)
               PLF(NS)%RU  =PUS(NS)
               PLF(NS)%RUPL=0.D0
               PLF(NS)%RNUC=0.D0
               PLF(NS)%RZCL=PZCL(NS)
            ENDDO
         ELSE
            DO NS=1,NSMAX
               IF(RHOL.LE.RHOEDG) THEN
                  FACTN=(1.D0-RHOL**PROFN1(NS))**PROFN2(NS)
                  FACTT=(1.D0-RHOL**PROFT1(NS))**PROFT2(NS)
                  FACTU=(1.D0-RHOL**PROFU1(NS))**PROFU2(NS)
               ELSE
                  FNX=(1.D0-RHOEDG**PROFN1(NS))**PROFN2(NS)
                  DFNX=-PROFN1(NS)*PROFN2(NS)*RHOEDG**(PROFN1(NS)-1.D0) &
                       *(1.D0-RHOEDG**PROFN1(NS))**(PROFN2(NS)-1.D0)
                  AN= 3.D0*FNX/(1.D0-RHOEDG)**2+DFNX/(1.D0-RHOEDG)
                  BN=-2.D0*FNX/(1.D0-RHOEDG)**3-DFNX/(1.D0-RHOEDG)**2
                  FACTN=AN*(1.D0-RHOL)**2+BN*(1.D0-RHOL)**3

                  FTX=(1.D0-RHOEDG**PROFT1(NS))**PROFT2(NS)
                  DFTX=-PROFT1(NS)*PROFT2(NS)*RHOEDG**(PROFT1(NS)-1.D0) &
                       *(1.D0-RHOEDG**PROFT1(NS))**(PROFT2(NS)-1.D0)
                  AT= 3.D0*FTX/(1.D0-RHOEDG)**2+DFTX/(1.D0-RHOEDG)
                  BT=-2.D0*FTX/(1.D0-RHOEDG)**3-DFTX/(1.D0-RHOEDG)**2
                  FACTT=AT*(1.D0-RHOL)**2+BT*(1.D0-RHOL)**3

                  FUX=(1.D0-RHOEDG**PROFU1(NS))**PROFU2(NS)
                  DFUX=-PROFU1(NS)*PROFU2(NS)*RHOEDG**(PROFU1(NS)-1.D0) &
                       *(1.D0-RHOEDG**PROFU1(NS))**(PROFU2(NS)-1.D0)
                  AU= 3.D0*FUX/(1.D0-RHOEDG)**2+DFUX/(1.D0-RHOEDG)
                  BU=-2.D0*FUX/(1.D0-RHOEDG)**3-DFUX/(1.D0-RHOEDG)**2
                  FACTU=AU*(1.D0-RHOL)**2+BU*(1.D0-RHOL)**3
               ENDIF

               PLF(NS)%RN  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               PLF(NS)%RTPR=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               PLF(NS)%RTPP=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               PLF(NS)%RU  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               PLF(NS)%RUPL=0.D0
               PLF(NS)%RNUC=0.D0
               PLF(NS)%RZCL=PZCL(NS)
               IF(RHOL.LT.RHOITB(NS)) THEN
                  FACTITB =(1.D0-(RHOL/RHOITB(NS))**4)**2
                  PLF(NS)%RN  =PLF(NS)%RN  +PNITB(NS)*FACTITB
                  PLF(NS)%RTPR=PLF(NS)%RTPR+PTITB(NS)*FACTITB
                  PLF(NS)%RTPP=PLF(NS)%RTPP+PTITB(NS)*FACTITB
                  PLF(NS)%RU  =PLF(NS)%RU  +PUITB(NS)*FACTITB
                  PLF(NS)%RUPL=0.D0
                  PLF(NS)%RNUC=0.D0
                  PLF(NS)%RZCL=PZCL(NS)
               ENDIF
            ENDDO
         ENDIF

      CASE(8)
         DO NS=1,NSMAX
            CALL WMSPL_PROF(Rhol,NS,RN_PL(NS),RT_PL(NS))
         ENDDO

!----  Modification for charge neutrality after spline interpolation

         VAL=0.D0
         DO NS=2,NSMAX-1
            VAL=VAL+PZ(NS)*RN_PL(NS)
         ENDDO
         RN_PL(NSMAX)=(RN_PL(1)-VAL)/PZ(NSMAX)

!----

         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               PLF(NS)%RN=PNS(NS)
               IF (NS.EQ.1.OR.NS.GT.1) THEN
                  CALL WMSPL_PROF(1.D0,NS,RN_PL(NS),RT_PL(NS))
                  PLF(NS)%RTPR=RT_PL(NS)*1.D-3
                  PLF(NS)%RTPP=RT_PL(NS)*1.D-3
               ELSE
                  PLF(NS)%RTPR=PTS(NS)
                  PLF(NS)%RTPP=PTS(NS)
               ENDIF
               PLF(NS)%RU  =PUS(NS)
               PLF(NS)%RUPL=0.D0
               PLF(NS)%RNUC=0.D0
               PLF(NS)%RZCL=PZCL(NS)
            ENDDO
         ELSE
            DO NS=1,NSMAX
               FACTN=(1.D0-RHOL**PROFN1(NS))**PROFN2(NS)
               FACTT=(1.D0-RHOL**PROFT1(NS))**PROFT2(NS)
               FACTU=(1.D0-RHOL**PROFU1(NS))**PROFU2(NS)
               IF (NS.EQ.1.OR.NS.GT.1) THEN
                  PLF(NS)%RN  = RN_PL(NS)*1.D-20
                  PLF(NS)%RTPR= RT_PL(NS)*1.D-3
                  PLF(NS)%RTPP= RT_PL(NS)*1.D-3
               ELSE
                  PLF(NS)%RN  =((PN(NS)  -PNS(NS))*FACTN+PNS(NS))
                  PLF(NS)%RTPR=((PTPR(NS)-PTS(NS))*FACTT+PTS(NS))
                  PLF(NS)%RTPP=((PTPP(NS)-PTS(NS))*FACTT+PTS(NS))
               ENDIF
               PLF(NS)%RU  = (PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               PLF(NS)%RUPL=0.D0
               PLF(NS)%RNUC=0.D0
               PLF(NS)%RZCL=PZCL(NS)
            ENDDO
         ENDIF

      CASE(9)
         CALL pl_bpsd_get(RHOL,RN_PL,RTPR_PL,RTPP_PL,RU_PL,RUPL_PL)
         DO NS=1,NSMAX
            PLF(NS)%RN  =RN_PL(NS)
            PLF(NS)%RTPR=RTPR_PL(NS)
            PLF(NS)%RTPP=RTPP_PL(NS)
            PLF(NS)%RU  =RU_PL(NS)
            PLF(NS)%RUPL=RUPL_PL(NS)
            PLF(NS)%RNUC=0.D0
            PLF(NS)%RZCL=PZCL(NS)
         ENDDO

      CASE(21)
         DO NS=1,NSMAX
            CALL pl_read_trdata(RHOL,NS,PNL,PTL)
            PLF(NS)%RN  =PNL
            PLF(NS)%RTPR=PTL
            PLF(NS)%RTPP=PTL
            PLF(NS)%RU  =0.D0
            PLF(NS)%RUPL=0.D0
            PLF(NS)%RNUC=0.D0
            PLF(NS)%RZCL=PZCL(NS)
            IF(NS.EQ.1) WRITE(6,'(A,1P3E12.4)') 'rhol,pnl,ptl=',RHOL,PNL,PTL
         ENDDO

      END SELECT

      RETURN
    END SUBROUTINE pl_prof

    SUBROUTINE pl_grad(rhon,grd)
      USE plcomm,ONLY: ra,nsmax
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: rhon
      TYPE(pl_grd_type),DIMENSION(NSMAX),INTENT(OUT):: grd
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf1,plf2
      INTEGER:: ns
      REAL(rkind):: drhon=1.D-6

      IF(rhon.LE.0.D0.OR.rhon.GT.1.D0) THEN
         DO ns=1,nsmax
            grd(ns)%grdn  =0.D0
            grd(ns)%grdtpr=0.D0
            grd(ns)%grdtpp=0.D0
            grd(ns)%grdu  =0.D0
            grd(ns)%grdupl=0.D0
         END DO
      ELSE
         CALL pl_prof(rhon+drhon,plf1)
         CALL pl_prof(rhon-drhon,plf2)
         DO ns=1,nsmax
            grd(ns)%grdn  =(plf2(ns)%rn  -plf1(ns)%rn  )/(2.D0*rhon*RA)
            grd(ns)%grdtpr=(plf2(ns)%rtpr-plf1(ns)%rtpr)/(2.D0*rhon*RA)
            grd(ns)%grdtpp=(plf2(ns)%rtpp-plf1(ns)%rtpp)/(2.D0*rhon*RA)
            grd(ns)%grdu  =(plf2(ns)%ru  -plf1(ns)%ru  )/(2.D0*rhon*RA)
            grd(ns)%grdupl=(plf2(ns)%rupl-plf1(ns)%rupl)/(2.D0*rhon*RA)
         END DO
      END IF
      RETURN
    END SUBROUTINE pl_grad

    SUBROUTINE pl_bpsd_get(rho,rn,rtpr,rtpp,ru,rupl)
      USE plcomm
      USE bpsd
      REAL(rkind),INTENT(IN):: rho
      REAL(rkind),DIMENSION(nsmax),INTENT(OUT):: rn,rtpr,rtpp,ru,rupl
      TYPE(bpsd_plasmaf_type),save :: plasmaf
      INTEGER(ikind):: ns,ierr

      plasmaf%nrmax=1
      plasmaf%rho(1)=rho
      CALL bpsd_get_data(plasmaf,ierr)
      DO ns=1,min(nsmax,plasmaf%nsmax)
         rn(ns)  =plasmaf%data(1,ns)%density
         rtpr(ns)=plasmaf%data(1,ns)%temperature_para
         rtpp(ns)=plasmaf%data(1,ns)%temperature_perp
         ru(ns)  =plasmaf%data(1,ns)%velocity_tor
         rupl(ns)=plasmaf%data(1,ns)%velocity_pol
      ENDDO
    END SUBROUTINE pl_bpsd_get

!     ***** AVERAGE MINOR RADIUS FOR PARABOLIC Q PROFILE *****

    FUNCTION rsrhon(rhon)

      USE plcomm,ONLY: RA,MODELG
      IMPLICIT NONE
      REAL(rkind), INTENT(IN):: RHON
      REAL(rkind):: rsrhon
      REAL(rkind):: RHOL,RRMINL,RRMAXL

      SELECT CASE(MODELG)
      CASE(0:2)
         IF(RHON.LE.0.D0) THEN
           RHOL=0.D0
         ELSE
           RHOL=RHON
         ENDIF
         rsrhon=RHOL*RA
      CASE(3,5,8)
         CALL GETRMN(RHON,RRMINL)
         CALL GETRMX(RHON,RRMAXL)
         rsrhon=0.5D0*(RRMAXL-RRMINL)
      END SELECT
      RETURN
    END FUNCTION rsrhon

!     ****** CALCULATE Q PROFILE ******

    SUBROUTINE pl_qprf(RHON,QL)

      USE plcomm,ONLY: MODELG,MODELQ,Q0,QA,RHOMIN,QMIN,PI,RA,BB,RMU0, &
           RR,RIP,PROFJ
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: RHON
      REAL(rkind),INTENT(OUT):: QL
      REAL(8) :: RHOL, QSA0, QSAA

      IF(RHON.LE.0.D0) THEN
         RHOL=0.D0
      ELSE
         RHOL=RHON
      ENDIF

      SELECT CASE(MODELG)
      CASE(0:2)
         IF(MODELQ.EQ.0) THEN
            IF(RHOL.GT.1.D0) THEN
               QL = QA*RHOL**2
            ELSEIF(RHOMIN.LE.0.D0)THEN
               QL =(Q0-QA)*(1.D0-RHOL**2)+QA
            ELSE
               QSA0    =1.D0/Q0-1.D0/QMIN
               QSAA    =1.D0/QA-1.D0/QMIN
               IF(RHOL.LE.RHOMIN)THEN
                  QL =1.D0/(1.D0/Q0-QSA0*(3.D0*RHOL**2/RHOMIN**2 &
     &                                   -2.D0*RHOL**3/RHOMIN**3))
               ELSE
                  QL =1.D0/(1.D0/QMIN+3.D0*QSA0*(RHOL-RHOMIN)**2/RHOMIN**2 &
     &                    +(QSAA     -3.D0*QSA0*(1.D0-RHOMIN)**2/RHOMIN**2) &
     &                        *(RHOL-RHOMIN)**3/(1.D0-RHOMIN)**3)
               ENDIF
            ENDIF
         ELSEIF(MODELQ.EQ.1) THEN
            IF(RR.EQ.0.D0) THEN
               QA=1.D6
               Q0=1.D6
            ELSE
               QA=2.D0*PI*RA*RA*BB/(RMU0*RIP*1.D6*RR)
               Q0=QA/(1.D0+PROFJ)
            ENDIF
            IF(RHOL.GE.1.D0) THEN
               QL = QA*RHOL**2
            ELSEIF(RHOL.LE.1.D-30) THEN
               QL = Q0
            ELSE
               QL=QA*RHOL**2/(1.D0-(1.D0-RHOL**2)**(PROFJ+1.D0))
            ENDIF
         ENDIF
      CASE DEFAULT
         CALL GETQP(RHOL,QL)
      END SELECT
      RETURN
    END SUBROUTINE pl_qprf

!     ****** CALCULATE BMIN ON MAG SURFACE ******

    SUBROUTINE pl_bminmax(RHON,BMIN,BMAX)

      USE plcomm,ONLY: MODELG,BB,RR,RA
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: RHON
      REAL(rkind),INTENT(OUT):: BMIN,BMAX
      REAL(rkind) :: RS,QL,BT,BP,BTMIN,BTMAX

      SELECT CASE(MODELG)
      CASE(0,1)
         RS=RA*RHON
         BT= BB
         CALL GETQP(RHON,QL)
         BP= RS*BB/(RR*QL)
         BMIN= SQRT(BT**2+BP**2)
         BMAX= SQRT(BT**2+BP**2)
      CASE(2)
         RS=RA*RHON
         BTMIN= BB/(1+RS/RR)
         BTMAX= BB/(1-RS/RR)
         CALL pl_qprf(RHON,QL)
         BP= RS*BB/(RR*QL)
         BMIN= SQRT(BTMIN**2+BP**2)
         BMAX= SQRT(BTMAX**2+BP**2)
      CASE(3,5,8)
         CALL GET_BMINMAX(RHON,BMIN,BMAX)
      END SELECT
      RETURN
    END SUBROUTINE pl_bminmax

!     ****** CALCULATE RRMIN AND RRMAX ON MAG SURFACE ******

    SUBROUTINE pl_rrminmax(RHON,RRMINL,RRMAXL)

      USE plcomm,ONLY: MODELG,RR
      IMPLICIT NONE
      REAL(8),INTENT(IN):: RHON
      REAL(8),INTENT(OUT):: RRMINL,RRMAXL
      REAL(8)  :: RS

      SELECT CASE(MODELG)
      CASE(0)
         RRMINL=0.D0
         RRMAXL=0.D0
      CASE(1)
         RRMINL=RR
         RRMAXL=RR
      CASE(2)
         RS=rsrhon(RHON)
         RRMINL=RR-RS
         RRMAXL=RR+RS
      CASE(3,5,8)
         CALL GETRMN(RHON,RRMINL)
         CALL GETRMX(RHON,RRMAXL)
      END SELECT
      RETURN
    END SUBROUTINE pl_rrminmax

!     ****** CALCULATE RRMIN AND RRMAX ON MAG SURFACE ******

    SUBROUTINE pl_dvdrho(RHON,DVDRHO)

      USE plcomm,ONLY: MODELG,TWOPI,RR,RA,RKAP
      IMPLICIT NONE
      REAL(8),INTENT(IN):: RHON
      REAL(8),INTENT(OUT):: DVDRHO

      SELECT CASE(MODELG)
      CASE(0,1,2)
         DVDRHO=TWOPI*RR*TWOPI*RA**2*RKAP*RHON
      CASE(3,5,8)
         CALL GET_DVDRHO(RHON,DVDRHO)
      END SELECT
      RETURN
    END SUBROUTINE pl_dvdrho

!     ***** PLASMA MAGNETIC AXIS *****

    SUBROUTINE pl_axis(RAXIS,ZAXIS)

      USE plcomm,ONLY: RR,MODELG
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT) :: RAXIS,ZAXIS

      SELECT CASE(MODELG)
      CASE(0)
         RAXIS=0.D0
         ZAXIS=0.D0
      CASE(1:2)
         RAXIS=RR
         ZAXIS=0.D0
      CASE(3,5,8)
         CALL GETAXS(RAXIS,ZAXIS)
      END SELECT
      RETURN
    END SUBROUTINE pl_axis

!     ***** PLASMA BOUNDARY *****

    SUBROUTINE pl_rzsu(RSU,ZSU,NSUM,NSUMAX)

      USE plcomm,ONLY: PI,RA,RKAP,RR,MODELG
      IMPLICIT NONE
      INTEGER(ikind),INTENT(IN):: NSUM
      REAL(rkind),DIMENSION(NSUM),INTENT(OUT) :: RSU,ZSU
      INTEGER(ikind),INTENT(OUT):: NSUMAX
      REAL(rkind)     :: DTH, TH
      INTEGER(ikind)  :: NSU

      SELECT CASE(MODELG)
      CASE(0:2)
         NSUMAX=NSUM
         DTH=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            TH=(NSU-1)*DTH
            RSU(NSU)=RR+     RA*COS(TH)
            ZSU(NSU)=   RKAP*RA*SIN(TH)
         ENDDO
      CASE(3,5,8)
         CALL GETRSU(RSU,ZSU,NSUM,NSUMAX)
      END SELECT
      RETURN
    END SUBROUTINE pl_rzsu

!     ***** FILE READ FOR TASK/WM (WMXPRF) *****

    subroutine pl_wmxprf(ierr)

      USE plcomm,ONLY: NSMAX,PZ
      USE plxprf
      IMPLICIT NONE
      REAL(rkind),DIMENSION(NXPRF,NXSPC):: PRFN,PRFT
      CHARACTER(LEN=80):: TRFILE='topics-data' ! fixed name
      INTEGER(ikind):: ierr,ifno,nr,ns,irc,n,i
      REAL(rkind):: val

      ierr = 0

!----  Open profile data file and read
!----  PRFNE, PRFTE is data at the point divided equally by rho 
!----    defined by toroidal magnetic flux

      IFNO=22
      OPEN ( IFNO, FILE=TRFILE, ERR=9995 )
      READ ( IFNO, '(I3)', END=9996, ERR=9996 ) NPRFMAX
      DO N=1,NPRFMAX
         READ ( IFNO, '(13E14.7)', END=9996, ERR=9996 ) &
             PRFRHO(N), (PRFN(N,I), I=1,NXSPC), &
                        (PRFT(N,I), I=1,NXSPC)
      ENDDO

!----  Modification for charge neutrality

!      DO NR=1,NPRFMAX
!         VAL=0.D0
!         DO NS=2,NSMAX-1
!            VAL=VAL+PZ(NS)*PRFN(NR,NS)
!         ENDDO
!         PRFN(NR,NSMAX)=(PRFN(NR,1)-VAL)/PZ(NSMAX)
!      ENDDO

      DO NR=1,NPRFMAX
         VAL=0.D0
         DO NS=2,NSMAX
            VAL=VAL+PZ(NS)*PRFN(NR,NS)
         ENDDO
         PRFN(NR,1)=VAL
      ENDDO

!----  Set coefficient for spline

      DO NS=1,NSMAX
         CALL SPL1D(PRFRHO,PRFN(1,NS),DERIV,UPRFN(1,1,NS), NPRFMAX,0,IRC)
         IF (IRC.NE.0) GO TO 9997
         CALL SPL1D(PRFRHO,PRFT(1,NS),DERIV,UPRFT(1,1,NS), NPRFMAX,0,IRC)
         IF (IRC.NE.0) GO TO 9997
      ENDDO

!----  Debug write

 8000 FORMAT(' N ',3X,'PRFRHO',6X,'PRFNE',6X,'PRFNI1',5X,'PRFNI2', &
                   5X,'PRFNI3',5X,'PRFNI4')
 8010 FORMAT(' N ',3X,'PRFRHO',6X,'PRFTE',6X,'PRFTI1',5X,'PRFTI2', &
                   5X,'PRFTI3',5X,'PRFTI4')
      GO TO 9999

 9995 WRITE(6,*) '==========  PLWMXPRF FILE OPEN ERROR  =========='
      GO TO 9999
 9996 WRITE(6,*) '==========  PLWMXPRF FILE READ ERROR  =========='
      GO TO 9999
 9997 WRITE(6,*) '==========  PLWMXPRF SPL1D ERROR  =========='

 9999 CLOSE( IFNO )
      RETURN
    END SUBROUTINE pl_wmxprf

!     ***** Interpolation of profile at a given point *****

    SUBROUTINE wmspl_prof(Rhol,NS,PNL,PTL)

      USE plcomm,ONLY: PNS,PTS,modeln
      USE plxprf
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: rhol   ! Normalized radius
      INTEGER(ikind),INTENT(IN):: NS  ! Particle species
      REAL(rkind),INTENT(OUT):: PNL   ! Density at Rhol
      REAL(rkind),INTENT(OUT):: PTL   ! Temperature at Rhol
      REAL(rkind):: PPL
      INTEGER(ikind):: IERR

      IF (Rhol.GT.1.0D0) THEN
         IF(modeln.EQ.1) THEN
            PNL = PNS(NS)
         ELSE
            PNL = 0.D0
         END IF
         PTL = PTS(NS)
      ELSE
         CALL SPL1DF(Rhol,PPL,PRFRHO,UPRFN(1,1,NS),NPRFMAX,IERR)
         PNL=PPL
         CALL SPL1DF(Rhol,PPL,PRFRHO,UPRFT(1,1,NS),NPRFMAX,IERR)
         PTL=PPL
      ENDIF

      RETURN
    END SUBROUTINE wmspl_prof

!     ****** Coordinate conversion ******

    SUBROUTINE pl_getRZ(rhon,chip,R,Z)

      USE plcomm, ONLY: RA,RR,MODELG,TWOPI

      IMPLICIT NONE
      REAL(rkind),INTENT(in):: rhon,chip
      REAL(rkind),INTENT(OUT):: R,Z

      SELECT CASE(MODELG)
      CASE(0)
         R=   RA*rhon
         Z=   RA*rhon*chip/TWOPI
      CASE(1)
         R=RR+RA*rhon*COS(chip)
         Z=   RA*rhon*SIN(chip)

      CASE(2)
         R=RR+RA*rhon*COS(chip)
         Z=   RA*rhon*SIN(chip)

      CASE(3,5,8)
         CALL GET_RZ(rhon,chip,R,Z)
      END SELECT
    END SUBROUTINE pl_getRZ

!     ****** Coordinate conversion and magnetic field ******

    SUBROUTINE pl_getRZB(rhon,chip,R,Z,BR,BZ,BT)

      USE plcomm, ONLY: RA,RR,MODELG,TWOPI

      IMPLICIT NONE
      REAL(rkind),INTENT(in):: rhon,chip
      REAL(rkind),INTENT(OUT):: R,Z,BR,BZ,BT
      REAL(rkind):: RHON_DUMMY

      SELECT CASE(MODELG)
      CASE(0)
         R=   RA*rhon
         Z=   RA*rhon*chip/TWOPI
      CASE(1,2)
         R=RR+RA*rhon*COS(chip)
         Z=   RA*rhon*SIN(chip)
      CASE(3,5,8)
         CALL GET_RZ(rhon,chip,R,Z)
      END SELECT
      CALL pl_mag_rz(R,0.D0,Z,BR,BZ,BT,RHON_DUMMY)
      RETURN
    END SUBROUTINE pl_getRZB

!     ****** calculate magnetci field ******

    SUBROUTINE pl_getB(rhon,chip,B)

      USE plcomm, ONLY: RA,RR,MODELG,TWOPI

      IMPLICIT NONE
      REAL(rkind),INTENT(in):: rhon,chip
      REAL(rkind),INTENT(OUT):: B
      REAL(rkind):: R,Z,BR,BZ,BT

      CALL pl_getRZB(rhon,chip,R,Z,BR,BZ,BT)
      B=SQRT(BR*BR+BZ*BZ+BT*BT)

      RETURN
    END SUBROUTINE pl_getB

  END MODULE plprof
