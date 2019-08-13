!     $Id$

      module plprof_

      use plcomm

      type position_XYZ
         real(rkind):: X,Y,Z
      end type position_XYZ
      type position_RZP
         real(rkind):: R,Z,phi
      end type position_RZP
      type position_mag
         real(rkind):: rhon,chi,phi
      end type position_mag

      type mag_XYZ
         real(rkind):: BX,BY,BZ,Babs
      end type mag_XYZ
      type mag_RZP
         real(rkind):: BR,BZ,Bphi,Babs
      end type mag_RZP
      type mag_local
         real(rkind):: Bp,Bt,Babs
      end type mag_local
      type mag_sub
         real(kind):: Bsub1,Bsub2,Bsub3
      end type mag_sub
      type mag_sup
         real(kind):: Bsup1,Bsup2,Bsup3
      end type mag_sup

      contains

!     ****** convert position vector *****

        function convert_pos_XYZ_RZP(XYZ)
          implicit none
          type(position_XYZ),intent(in):: XYZ
          type(position_RZP),intent(out):: convert_pos_XYZ_RZP

          SELECT CASE(modelg)
          CASE(0:1) ! slab,cylinder
             convert_pos_XYZ_RZP%R=XYZ%X
             convert_pos_XYZ_RZP%Z=XYZ%Z
             convert_pos_XYZ_RZP%phi=-XYZ%Y/(2.d0*pi*RR)
          CASE DEFAULT
             convert_pos_XYZ_RZP%R=SQRT(XYZ%X**2+XYZ%Y**2)
             convert_pos_XYZ_RZP%Z=XYZ%Z
             convert_pos_XYZ_RZP%phi=ATAN2(-XYZ%Y,XYZ%X)
          END SELECT
          return
        end function convert_pos_XYZ_RZP

        function convert_pos_RZP_XYZ(RZP)
          implicit none
          type(position_RZP),intent(in):: RZP
          type(position_XYZ),intent(out):: convert_pos_RZP_XYZ

          SELECT CASE(modelg)
          CASE(0:1) ! slab,cylinder
             convert_pos_RZP_XYZ%X=RZP%R
             convert_pos_RZP_XYZ%Y=-2.d0*pi*RR*RZP%phi
             convert_pos_RZP_XYZ%Z=RZP%Z
          CASE DEFAULT
             convert_pos_RZP_XYZ%X=RZP%R*COS(-RZP%phi)
             convert_pos_RZP_XYZ%Y=RZP%R*SIN(-RZP%phi)
             convert_pos_RZP_XYZ%Z=RZP%Z
          END SELECT
          return
        end function convert_pos_RZP_XYZ

        function convert_pos_RZP_mag(RZP,posmag)
          implicit none
          type(position_RZP),intent(in):: RZP
          type(position_mag),intent(out):: convert_pos_RZP_mag
          real(rkind):: rhon,chi,phi

          SELECT CASE(modelg)
          CASE(0) ! slab
             rhon=(RZP%R-RR)/RA
             chi=ATAN2(RZP%Z,RA)
             phi=RZP%phi
          CASE(1) ! cyrinder
             rhon=SQRT((RZP%R-RR)**2+(RZP%Z/RKAP)**2)/RA
             chi=ATAN2(RZP%Z,RZP%R-RR)
             phi=RZP%phi
          CASE(2) ! torus
             rhon=SQRT((RZP%R-RR)**2+(RZP%Z/RKAP)**2)/RA
             chi=ATAN2(RZP%Z/RKAP,RZP%R-RR)
             phi=RZP%phi
          CASE(3,5)

      ELSEIF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN
         RL=SQRT(X**2+Y**2)
         PP=0.D0
         CALL GETRZ(RL,Z,PP,BR,BZ,BT,RHON)
!         WRITE(6,'(1P6E12.4)') RL,ZZ,BR,BZ,BT,RHON
         RCOST=X/RL
         RSINT=Y/RL
         BX = BR*RCOST-BT*RSINT
         BY = BR*RSINT+BT*RCOST
      ENDIF
C
      BABS = SQRT(BX**2+BY**2+BZ**2)
C
      IF(BABS.LE.0.D0) THEN
         BNX = 0.D0
         BNY = 1.D0
         BNZ = 0.D0
      ELSE
         BNX = BX/BABS
         BNY = BY/BABS
         BNZ = BZ/BABS
      ENDIF
C
      XPOS_LOC=X
      YPOS_LOC=Y
      ZPOS_LOC=Z
      RHON_LOC=RHON




----
             CALL PLQPRF(RHON,QL)
         IF(RS.LE.0.D0) QL=-QL
         BT   = BB
         BP   = RS*BT/(RR*QL)
         BX   = 0.D0
         BY   = BB
         BZ   = BP
C
      ELSEIF(MODELG.EQ.1) THEN
         RS =SQRT((X-RR)**2+Z**2)
         RHON=RS/RA
         IF(RS.LE.0.D0) THEN
            BX   = 0.D0
            BY   = BB
            BZ   = 0.D0
         ELSE
            CALL PLQPRF(RHON,QL)
            RSINT= Z/RS
            RCOST= (X-RR)/RS
            BT   = BB
            BP   = RS*BT/(RR*QL)
            BX   =-BP*RSINT
            BY   = BB
            BZ   = BP*RCOST
         ENDIF
C
      ELSEIF(MODELG.EQ.2) THEN
         RL=SQRT(X**2+Y**2)
         RS =SQRT((RL-RR)**2+Z**2)
         RHON=RS/RA
         IF(RS.LE.0.D0) THEN
            BT   = BB
            BR   = 0.D0
            BZ   = 0.D0
         ELSE
            CALL PLQPRF(RHON,QL)
            RSINP= Z/RS
            RCOSP= (RL-RR)/RS
            BT   = BB/(1.D0+RS*RCOSP/RR)
            BP   = RS*BT/(RR*QL)
            BR   =-BP*RSINP
            BZ   = BP*RCOSP
         ENDIF
         RCOST=X/RL
         RSINT=Y/RL
         BX = BR*RCOST-BT*RSINT
         BY = BR*RSINT+BT*RCOST
C
      ELSEIF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN
         RL=SQRT(X**2+Y**2)
         PP=0.D0
         CALL GETRZ(RL,Z,PP,BR,BZ,BT,RHON)
!         WRITE(6,'(1P6E12.4)') RL,ZZ,BR,BZ,BT,RHON
         RCOST=X/RL
         RSINT=Y/RL
         BX = BR*RCOST-BT*RSINT
         BY = BR*RSINT+BT*RCOST
      ENDIF
C
      BABS = SQRT(BX**2+BY**2+BZ**2)
C
      IF(BABS.LE.0.D0) THEN
         BNX = 0.D0
         BNY = 1.D0
         BNZ = 0.D0
      ELSE
         BNX = BX/BABS
         BNY = BY/BABS
         BNZ = BZ/BABS
      ENDIF
C
      XPOS_LOC=X
      YPOS_LOC=Y
      ZPOS_LOC=Z
      RHON_LOC=RHON
C
      RETURN
      END
C
!     ****** CALCULATE PLASMA PROFILE ******
C
      SUBROUTINE PLPROF(RHON)
C
      INCLUDE '../pl/plcomm.inc'
      INCLUDE '../pl/plcom2.inc'
      DIMENSION RNPL(NSM),RTPL(NSM),RUPL(NSM)
C
      IF(RHON.LE.0.D0) THEN
         RHOL=0.D0
      ELSEIF(RHON.GE.1.D0) THEN
         RHOL=1.D0
      ELSE
         RHOL=RHON
      ENDIF
C
      IF(MODELN.EQ.0) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =0.D0
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            FACTN=(1.D0-RHOL**PROFN1)**PROFN2
            FACTT=(1.D0-RHOL**PROFT1)**PROFT2
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
C
            DO NS=1,NSMAX
               RN(NS)  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               RTPR(NS)=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               RTPP(NS)=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF
C
      ELSEIF(MODELN.EQ.1) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =PNS(NS)
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            FACTN=(1.D0-RHOL**PROFN1)**PROFN2
            FACTT=(1.D0-RHOL**PROFT1)**PROFT2
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
C
            DO NS=1,NSMAX
               RN(NS)  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               RTPR(NS)=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               RTPP(NS)=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF
C
      ELSEIF(MODELN.EQ.2) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =PNS(NS)
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            PTOT=PTOT*1.D20*1.D3/AEE
            CALL GETPP(0.D0,PL0)
            CALL GETPP(PSIN,PL)
            FACT=SQRT(PTOT*PL/PL0)
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            DO NS=1,NSMAX
               RN(NS)  =PN(NS)  *FACT
               RTPR(NS)=PTPR(NS)*FACT
               RTPP(NS)=PTPP(NS)*FACT
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ENDIF
C
      ELSEIF(MODELN.EQ.3) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =PNS(NS)
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            IF(RHOL.LE.RHOEDG) THEN
               FACTN=(1.D0-RHOL**PROFN1)**PROFN2
               FACTT=(1.D0-RHOL**PROFT1)**PROFT2
               FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            ELSE
               FNX=(1.D0-RHOEDG**PROFN1)**PROFN2
               DFNX=-PROFN1*PROFN2*RHOEDG**(PROFN1-1.D0)
     &             *(1.D0-RHOEDG**PROFN1)**(PROFN2-1.D0)
               AN= 3*FNX/(1.D0-RHOEDG)**2+DFNX/(1.D0-RHOEDG)
               BN=-2*FNX/(1.D0-RHOEDG)**3-DFNX/(1.D0-RHOEDG)**2
               FACTN=AN*(1-RHOL)**2+BN*(1-RHOL)**3
C
               FTX=(1.D0-RHOEDG**PROFT1)**PROFT2
               DFTX=-PROFT1*PROFT2*RHOEDG**(PROFT1-1.D0)
     &             *(1.D0-RHOEDG**PROFT1)**(PROFT2-1.D0)
               AT= 3*FTX/(1.D0-RHOEDG)**2+DFTX/(1.D0-RHOEDG)
               BT=-2*FTX/(1.D0-RHOEDG)**3-DFTX/(1.D0-RHOEDG)**2
               FACTT=AT*(1-RHOL)**2+BT*(1-RHOL)**3
C
               FUX=(1.D0-RHOEDG**PROFU1)**PROFU2
               DFUX=-PROFU1*PROFU2*RHOEDG**(PROFU1-1.D0)
     &             *(1.D0-RHOEDG**PROFU1)**(PROFU2-1.D0)
               AU= 3*FUX/(1.D0-RHOEDG)**2+DFUX/(1.D0-RHOEDG)
               BU=-2*FUX/(1.D0-RHOEDG)**3-DFUX/(1.D0-RHOEDG)**2
               FACTU=AU*(1-RHOL)**2+BU*(1-RHOL)**3
            ENDIF
C
            DO NS=1,NSMAX
               RN(NS)  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               RTPR(NS)=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               RTPP(NS)=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF
C
      ELSEIF(MODELN.EQ.9) THEN
         CALL PLDATA_GETPL(RHOL,RNPL,RTPL,RUPL)
         DO NS=1,NSMAX
            RN(NS)  =RNPL(NS)
            RTPR(NS)=RTPL(NS)
            RTPP(NS)=RTPL(NS)
            RU(NS)  =RUPL(NS)
            RZCL(NS)=0.D0
         ENDDO
         IF(NSMAX.EQ.6) THEN
            RN(2)=RN(2)-RN(5)
            RN(4)=RN(4)-RN(6)
         ENDIF
      ENDIF
      RHON_LOC=RHON
C
      RETURN
      END
C
!     ****** CALCULATE Q PROFILE ******
C
      SUBROUTINE PLQPRF(RHON,QL)
C
      INCLUDE '../pl/plcomm.inc'
C
      IF(RHON.LE.0.D0) THEN
         RHOL=0.D0
      ELSE
         RHOL=RHON
      ENDIF
C
      IF(MODELG.LE.2) THEN
         IF(MODELQ.EQ.0) THEN
            IF(RHOL.GT.1.D0) THEN
               QL = QA*RHOL**2
            ELSEIF(RHOMIN.LE.0.D0)THEN
               QL =(Q0-QA)*(1-RHOL**2)+QA
            ELSE
               QSA0    =1/Q0-1/QMIN
               QSAA    =1/QA-1/QMIN
               IF(RHOL.LE.RHOMIN)THEN
                  QL =1/(1/Q0-QSA0*(3*RHOL**2/RHOMIN**2
     &                       -2*RHOL**3/RHOMIN**3))
               ELSE
                  QL =1/(1/QMIN+3*QSA0*(RHOL-RHOMIN)**2/RHOMIN**2
     &                  +(QSAA-3*QSA0*(1-RHOMIN)**2/RHOMIN**2)
     &                    *(RHOL-RHOMIN)**3/(1-RHOMIN)**3)
               ENDIF
            ENDIF
         ELSEIF(MODELQ.EQ.1) THEN
            QA=2.D0*PI*RA*RA*BB/(RMU0*RIP*1.D6*RR)
            Q0=QA/(1.D0+PROFJ)
            IF(RHOL.GE.1.D0) THEN
               QL = QA*RHOL**2
            ELSEIF(RHOL.LE.1.D-30) THEN
               QL = Q0
            ELSE
               QL=QA*RHOL**2/(1.D0-(1.D0-RHOL**2)**(PROFJ+1.D0))
            ENDIF
         ENDIF
      ELSE
         CALL GETQP(RHOL,QL)
      ENDIF
      RETURN
      END
C
!     ****** CALCULATE BMIN ON MAG SURFACE ******
C
      SUBROUTINE PLBMIN(RHON,BMINL)
C
      INCLUDE '../pl/plcomm.inc'
C
      IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
         RS=RSRHON(RHON)
         BMINT= BB
         CALL GETQP(RHON,QL)
         BMINP= RS*BMINT/(RR*QL)
         BMINL= SQRT(BMINT**2+BMINP**2)
      ELSEIF(MODELG.EQ.2) THEN
         RS=RSRHON(RHON)
         BMINT= BB/(1+RS/RR)
         CALL PLQPRF(RHON,QL)
         BMINP= RS*BMINT/((RR+RS)*QL)
         BMINL= SQRT(BMINT**2+BMINP**2)
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRMX(RHON,RRMAXL)
         PP=0.D0
         Z=0.D0
         CALL GETRZ(RRMAXL,Z,PP,BR,BZ,BT,RHONL)
!         write(6,'(A,1P3E12.4)') 'RHON,BR,BZ      =',RHON,BR,BZ
         BMINL = SQRT(BR**2+BT**2+BZ**2)
C - old -
!         BTL=BB*RR/RRMAXL
!         CALL PLQPRF(RHON,QL)
!         write(6,'(A,1P3E12.4)') 'RHON,RR,QL      =',RHON,RR,QL
!         write(6,'(A,1P3E12.4)') 'RHON,BTL,BT     =',RHON,BTL,BT
!         RS=RSRHON(RHON)
!         BPL=RS*BTL/(RR*QL)
!         write(6,'(A,1P3E12.4)') 'RHON,BPL,BTP    =',RHON,BPL,
!     &        SQRT(BR**2+BZ**2)
!         BMIN1=SQRT(BTL**2+BPL**2)
!         write(6,'(A,1P3E12.4)') 'RHON,BMINL,BMIN1=',RHON,BMINL,BMIN1
!         pause
      ENDIF
      RETURN
      END
C
!     ****** CALCULATE RRMIN AND RRMAX ON MAG SURFACE ******
C
      SUBROUTINE PLRRMX(RHON,RRMINL,RRMAXL)
C
      INCLUDE '../pl/plcomm.inc'
C
      IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
         RRMINL=RR
         RRMAXL=RR
      ELSEIF(MODELG.EQ.2) THEN
         RS=RSRHON(RHON)
         RRMINL=RR-RS
         RRMAXL=RR+RS
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRMN(RHON,RRMINL)
         CALL GETRMX(RHON,RRMAXL)
      ENDIF
      RETURN
      END
C
!     ***** AVERAGE MINOR RADIUS FOR PARABOLIC Q PROFILE *****
C
      FUNCTION RSRHON(RHON)
C
      INCLUDE '../pl/plcomm.inc'
C
      IF(MODELG.LT.3) THEN
         IF(RHON.LE.0.D0) THEN
           RHOL=0.D0
         ELSE
           RHOL=RHON
         ENDIF
         RSRHON=RHOL*RA
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRMN(RHON,RRMINL)
         CALL GETRMX(RHON,RRMAXL)
         RSRHON=0.5D0*(RRMAXL-RRMINL)
      ENDIF
      RETURN
      END
C
!     ***** PLASMA MAGNETIC AXIS *****
C
      SUBROUTINE PLAXIS(RAXIS,ZAXIS)
C
      INCLUDE '../pl/plcomm.inc'
C
      IF(MODELG.LT.3) THEN
         RAXIS=RR
         ZAXIS=0.D0
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETAXS(RAXIS,ZAXIS)
      ENDIF
      RETURN
      END
C
!     ***** PLASMA BOUNDARY *****
C
      SUBROUTINE PLRZSU(RSU,ZSU,NM,NSUMAX)
C
      INCLUDE '../pl/plcomm.inc'
C
      DIMENSION RSU(NM),ZSU(NM)
C
      IF(MODELG.EQ.0)THEN
         NSUMAX=5
         RSU(1)=RR-RA
         ZSU(1)=  -RA
         RSU(2)=RR+RA
         ZSU(2)=  -RA
         RSU(3)=RR+RA
         ZSU(3)=  +RA
         RSU(4)=RR-RA
         ZSU(4)=  +RA
         RSU(5)=RR-RA
         ZSU(5)=  -RA
      ELSEIF(MODELG.EQ.1.OR.MODELG.EQ.2) THEN
         NSUMAX=MIN(NM,101)
         DTH=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            TH=(NSU-1)*DTH
            RSU(NSU)=RR+     RA*COS(TH)
            ZSU(NSU)=   RKAP*RA*SIN(TH)
         ENDDO
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRSU(RSU,ZSU,NM,NSUMAX)
      ENDIF
      RETURN
      END



