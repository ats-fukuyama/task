! trprof.f90

MODULE trprof

  PRIVATE
  PUBLIC tr_prof,tr_prof_impurity,tr_prof_current

CONTAINS

!     ***********************************************************

!           SET INITIAL PROFILE

!     ***********************************************************

    SUBROUTINE tr_prof

      USE trcomm
      IMPLICIT NONE
      INTEGER:: NR, NS, NF
      REAL(8) :: PROF,qsurf,qaxis

      CALL TR_EDGE_DETERMINER(0)
      CALL TR_EDGE_SELECTOR(0)
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      DO NR=1,NRMAX
         RG(NR) = DBLE(NR)*DR
         RM(NR) =(DBLE(NR)-0.5D0)*DR
         VTOR(NR)=0.D0
         VPAR(NR)=0.D0
         VPRP(NR)=0.D0
         VPOL(NR)=0.D0
         WROT(NR)=0.D0
         DO NS=1,NSTM
            RN(NR,NS)=0.D0
            RT(NR,NS)=0.D0
            RU(NR,NS)=0.D0
         END DO
         DO NF=1,NFM
            RW(NR,NF)=0.D0
            RNF(NR,NF)=0.D0
            RTF(NR,NF)=0.D0
         END DO

         select case(MDLUF)
         case(1)
            IF(MDNI.EQ.0) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:2) = (PT(1:2)-PTS(1:2))*PROF+PTS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ELSE
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = RNU(1,NR,3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               RT(NR,1) = RTU(1,NR,1)
               RT(NR,2) = RTU(1,NR,2)
               RT(NR,3) = RTU(1,NR,3)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ENDIF

            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            SEX(NR,1:NSM) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0

            RNF(NR,1) = RNFU(1,NR)
            RNF(NR,2:NFM) = RNFU(1,NR)
            PBM(NR)   = PBMU(1,NR)
            WROT(NR)  = WROTU(1,NR)
            VTOR(NR)  = WROTU(1,NR)*RMJRHOU(1,NR)
         case(2)
            IF(MDNI.EQ.0) THEN ! ** MDNI **
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  RN(NR,1:2) = (PN(1:2)-PNS(1:2))*PROF+PNS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RN(NR,1) = RNU(1,NR,1)
                  RN(NR,2) = RNU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,3) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF +RNU(1,NRMAX,2)
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF +RNU(1,NRMAX,2)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:2) = (PT(1:2)-PTS(1:2))*PROF+PTS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               RT(NR,3) = RTU(1,NR,2)
               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:2) = (PT(1:2)-PTS(1:2))*PROF+PTS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ENDIF
            ELSE ! ** MDNI **
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  RN(NR,1:3) = (PN(1:3)-PNS(1:3))*PROF+PNS(1:3)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RN(NR,1:3) = RNU(1,NR,1:3)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF +RNU(1,NRMAX,2)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:3) = (PT(1:3)-PTS(1:3))*PROF+PTS(1:3)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1:3) = RTU(1,NR,1:3)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1:3) = RNU(1,NR,1:3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1:3) = RNU(1,NR,1:3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:3) = (PT(1:3)-PTS(1:3))*PROF+PTS(1:3)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1:3) = RTU(1,NR,1:3)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ENDIF
            ENDIF ! ** MDNI **

            PEX(NR,1)=PNBU(1,NR,1)
            PEX(NR,2)=PNBU(1,NR,2)
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0

            SEX(NR,1)=SNBU(1,NR,1)+SWLU(1,NR)/PZ(2)
            SEX(NR,2)=SNBU(1,NR,2)+SWLU(1,NR)
            SEX(NR,3)=0.D0
            SEX(NR,4)=0.D0
            PRF(NR,1:NSM)=0.D0
            RNF(NR,1)=RNFU(1,NR)
            RNF(NR,2:NFM)=0.D0
            PBM(NR)=0.D0
            WROT(NR) =WROTU(1,NR)
            VTOR(NR) =WROTU(1,NR)*RMJRHOU(1,NR)
         case(3)
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            RN(NR,1) = RNU(1,NR,1)
            RN(NR,2) = RNU(1,NR,2)
            RN(NR,3) = RNU(1,NR,3)
            RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
            RT(NR,1) = RTU(1,NR,1)
            RT(NR,2) = RTU(1,NR,2)
            RT(NR,3) = RTU(1,NR,3)
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)

            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            SEX(NR,1:NSM) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0
            RNF(NR,1:NFM) = 0.D0
            PBM(NR)=0.D0
            PRF(NR,4) = 0.D0
            WROT(NR)  = WROTU(1,NR)
            VTOR(NR)  = WROTU(1,NR)*RMJRHOU(1,NR)
         case default
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            RN(NR,1:NSM) = (PN(1:NSM)-PNS(1:NSM))*PROF+PNS(1:NSM)

            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            RT(NR,1:NSM) = (PT(1:NSM)-PTS(1:NSM))*PROF+PTS(1:NSM)

            PEX(NR,1:NSM) = 0.D0
            SEX(NR,1:NSM) = 0.D0
            PRF(NR,1:NSM) = 0.D0
            RNF(NR,1:NFM) = 0.D0
            PBM(NR)=0.D0
            WROT(NR)=0.D0
            VTOR(NR)=0.D0
         end select

         IF(MDLEQ0.EQ.1) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1)**PROFU2
            RN(NR,7) = (PN(7)-PNS(7))*PROF+PNS(7)
            RN(NR,8) = (PN(8)-PNS(8))*PROF+PNS(8)
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         ELSE
            ANNU(NR)=0.D0
         ENDIF

         RW(NR,1:NFM) = 0.D0

         SUMPBM=SUMPBM+PBM(NR)
      ENDDO
!      CALL PLDATA_SETR(RG,RM)
      CALL TR_EDGE_DETERMINER(1)
      CALL TR_EDGE_SELECTOR(1)

      qsurf=5.D0*(RA/RR)*(BB/RIPS)
      qaxis=0.5D0*qsurf
      DO nr=1,nrmax
         qpinv(nr)=1.D0/(qaxis+(qsurf-qaxis)*RG(nr)**2)
      END DO

    END SUBROUTINE tr_prof


    SUBROUTINE tr_prof_impurity

      USE trcomm
      IMPLICIT NONE
      REAL(8)   :: ANEAVE, ANI, ANZ, DILUTE, TE, TRZEC,TRZEFE
      INTEGER NR
      EXTERNAL TRZEC,TRZEFE
      
!     *** CALCULATE ANEAVE and ANC, ANFE ***

      ANEAVE=SUM(RN(1:NRMAX,1)*RM(1:NRMAX))*2.D0*DR
      SELECT CASE(MDLIMP)
      CASE(0)
         ANC (1:NRMAX)=0.D0
         ANFE(1:NRMAX)=0.D0
      CASE(1,3)
         DO NR=1,NRMAX
            ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC *1.D-2*RN(NR,1)
            ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE*1.D-2*RN(NR,1)
         END DO
      CASE(2,4)
         ANC (1:NRMAX)=PNC *RN(1:NRMAX,1)
         ANFE(1:NRMAX)=PNFE*RN(1:NRMAX,1)
      END SELECT

!     *** CALCULATE PZC,PZFE ***

      DO NR=1,NRMAX
         TE=RT(NR,1)
         PZC(NR)=TRZEC(TE)
         PZFE(NR)=TRZEFE(TE)
      ENDDO

!     *** Dilution of ION due to IMPURITY DENSITY ***

      IF(MDLUF.NE.3) THEN
         DO NR=1,NRMAX
            ANI = SUM(PZ(2:NSMAX)*RN(NR,2:NSMAX))     ! main ion charge density
            ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)   ! imputity ion charge den
            DILUTE = 1.D0-ANZ/ANI                     ! dilution factor
            IF(DILUTE.LT.0.D0) THEN
               WRITE(6,*) 'XX trprof: negative DILUTE: reduce PNC/PNFE'
               STOP
            END IF
            RN(NR,2:NSMAX) = RN(NR,2:NSMAX)*DILUTE    ! main ions diluted
         ENDDO
         PNSS(1)=PNS(1)
         PNSS(2:NSMAX)=PNS(2:NSMAX)*DILUTE
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1)=PNSA(1)
            PNSSA(2:NSMAX)=PNSA(2:NSMAX)*DILUTE
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ELSE
         PNSS(1:NSMAX)=PNS(1:NSMAX)
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1:NSMAX)=PNSA(1:NSMAX)
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ENDIF
      CALL TRZEFF

    END SUBROUTINE tr_prof_impurity

!     *** CALCULATE PROFILE OF AJ(R) ***

    SUBROUTINE tr_prof_current

      USE trcomm
      IMPLICIT NONE
      INTEGER:: NR
      REAL(8), DIMENSION(NRMAX) :: DSRHO
      REAL(8) :: FACT,FACTOR0,FACTORM,FACTORP,FCTR,PROF
      REAL(8) :: SUML

!     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***

      IF(MDLUF.EQ.1) THEN
         IF(MDLJQ.EQ.0) THEN ! *** MDLJQ ***
         NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTOR0*DR+FACTORM*RDPVRHOG(NR-1))/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO
         QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
         ELSE IF(MDLJQ.EQ.1) THEN ! *** MDLJQ ***
            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDPVRHOG(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            RDP(1:NRMAX)=RDPVRHOG(1:NRMAX)*DVRHOG(1:NRMAX)
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         NR=1
            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
            AJOH(NR)=AJ(NR)
         DO NR=2,NRMAX
            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
            AJOH(NR)=AJ(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO
         ELSE IF(MDLJQ.EQ.2) THEN ! *** MDLJQ ***
            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDPVRHOG(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            RDP(1:NRMAX)=RDPVRHOG(1:NRMAX)*DVRHOG(1:NRMAX)
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         NR=1
            AJ(NR) =1.d-6
            AJOH(NR)=AJ(NR)
         DO NR=2,NRMAX
            AJ(NR) =1.D-6
            AJOH(NR)=AJ(NR)
         ENDDO
         NR=1
            AJTOR(NR) =0.d0
         DO NR=2,NRMAX
            AJTOR(NR) =0.d0
         ENDDO
         ENDIF ! *** MDLJQ ***
      ELSEIF(MDLUF.EQ.2) THEN
         IF(MDLJQ.EQ.0) THEN  ! *** MDLJQ ***
            NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTOR0*DR+FACTORM*RDPVRHOG(NR-1))/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO

         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO
         QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))

         ELSEIF(MDLJQ.EQ.1) THEN ! *** MDLJQ ***

            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDPVRHOG(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            RDP(1:NRMAX)=RDPVRHOG(1:NRMAX)*DVRHOG(1:NRMAX)
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

            NR=1
               FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
               FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
               AJOH(NR)=AJ(NR)
            DO NR=2,NRMAX
               FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
               FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
               FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
               AJOH(NR)=AJ(NR)
            ENDDO
            NR=1
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORP=ABVRHOG(NR  )
               AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
            DO NR=2,NRMAX
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORM=ABVRHOG(NR-1)
               FACTORP=ABVRHOG(NR  )
               AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
            ENDDO
         ELSEIF(MDLJQ.EQ.2) THEN ! *** MDLJQ ***

            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDPVRHOG(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            RDP(1:NRMAX)=RDPVRHOG(1:NRMAX)*DVRHOG(1:NRMAX)
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

            NR=1
               AJ(NR) =1.D-6
               AJOH(NR)=AJ(NR)
            DO NR=2,NRMAX
               AJ(NR) =1.D-6
               AJOH(NR)=AJ(NR)
            ENDDO
            NR=1
               AJTOR(NR) =0.D0
            DO NR=2,NRMAX
               AJTOR(NR) =0.D0
            ENDDO
         ENDIF ! *** MDLJQ ***
         RIPA=ABVRHOG(NRAMAX)*RDPVRHOG(NRAMAX)*1.D-6 /(2.D0*PI*RMU0)
      ELSEIF(MDLUF.EQ.3) THEN
         DO NR=1,NRMAX
            AJOH(NR)=AJU(1,NR)
            AJ(NR)  =AJU(1,NR)
         ENDDO
         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         ENDDO
         BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX) =FACT*RDP(1:NRMAX)
         RDPVRHOG(1:NRMAX) =FACT*RDPVRHOG(1:NRMAX)
         AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)  =AJOH(1:NRMAX)
         BP(1:NRMAX)  =FACT*BP(1:NRMAX)
         QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      ELSE
         DO NR=1,NRMAX
            IF((1.D0-RM(NR)**ABS(PROFJ1)).LE.0.D0) THEN
               PROF=0.D0
            ELSE
               PROF= (1.D0-RM(NR)**ABS(PROFJ1))**ABS(PROFJ2)
            ENDIF
            AJOH(NR)= PROF
            AJ(NR)  = PROF
         ENDDO

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO

         RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)=FACT*RDP(1:NRMAX)
         RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
         AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)  =AJOH(1:NRMAX)
         BP(1:NRMAX)  =FACT*BP(1:NRMAX)
         QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      ENDIF
!      write(6,*) 'in trprof'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      pause
!     *** calculate q_axis ***
!      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))
      Q0=(20.D0*QP(1)-23.D0*QP(2)+8.D0*QP(3))/5.D0

!     calculate plasma current inside the calucated region (rho <= rhoa)
!     necessary for MDLEQB = 1 and MDLUF /= 0
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
!!         DSRHO(1:NRAMAX)=DVRHO(1:NRAMAX)/(2.D0*PI*RMJRHO(1:NRAMAX))
         DSRHO(1:NRAMAX)=DVRHO(1:NRAMAX)/(2.D0*PI*RR)
         RIPA=SUM(AJ(1:NRAMAX)*DSRHO(1:NRAMAX))*DR/1.D6
      ENDIF

!     *** THIS MODEL ASSUMES CONSTANT EZ ***

      IF(PROFJ1.LE.0.D0.OR.MDNCLS.EQ.1) THEN
         CALL TRZEFF
         CALL TRCFET
         IF(PROFJ1.GT.0.D0.AND.MDNCLS.EQ.1) GOTO 2000

         AJOH(1:NRMAX)=1.D0/ETA(1:NRMAX)
         AJ(1:NRMAX)  =1.D0/ETA(1:NRMAX)

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO
         BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)  =FACT*RDP(1:NRMAX)
         RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
         AJOH(1:NRMAX) =FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)   =AJOH(1:NRMAX)
         AJTOR(1:NRMAX)=FACT*AJTOR(1:NRMAX)
         BP(1:NRMAX)   =FACT*BP(1:NRMAX)
         QP(1:NRMAX)   =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      ENDIF
 2000 CONTINUE
      SUML=0.D0
      DO NR=1,NRMAX
         SUML=SUML+RDP(NR)*DR
         RPSI(NR)=SUML
         BPRHO(NR)=BP(NR)
      ENDDO

    END SUBROUTINE tr_prof_current

!     ***********************************************************

!           EDGE VALUE SELECTOR

!     ***********************************************************

      SUBROUTINE TR_EDGE_SELECTOR(NSW)

!        NSW = 0: store edge value; substitute rhoa value
!              1: restore original edge value

      USE TRCOMM, ONLY : MDLUF, NRAMAX, NSM, NSTM, PNSS, PNSSA, PTS, PTSA, RHOA, RN, RT, PNSSO,PTSO,PNSSAO,PTSAO
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NSW
      INTEGER(4) :: NS

      IF(RHOA.EQ.1.D0) RETURN

      IF(MDLUF.EQ.0) THEN
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS(NS)
               PTSO (NS)=PTS (NS)

               PNSS (NS)=PNSSAO(NS)
               PTS  (NS)=PTSAO (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSS (NS)=PNSSO(NS)
               PTS  (NS)=PTSO (NS)
            ENDDO
         ENDIF
      ELSE
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS (NS)
               PTSO (NS)=PTS  (NS)

               PNSS (NS)=PNSSA(NS)
               PTS  (NS)=PTSA (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSS (NS)=PNSSO(NS)
               PTS  (NS)=PTSO (NS)
            ENDDO
         ENDIF
      ENDIF
      RETURN
!
      ENTRY TR_EDGE_DETERMINER(NSW)

      IF(MDLUF.EQ.0.AND.RHOA.NE.1.D0) THEN
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSAO(NS)=PNSS(NS)
               PTSAO (NS)=PTS (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSSAO(NS)=RN(NRAMAX,NS)
               PTSAO (NS)=RT(NRAMAX,NS)
               PNSSA (NS)=PNSSAO(NS)
               PTSA  (NS)=PTSAO (NS)
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE TR_EDGE_SELECTOR
    END MODULE trprof
