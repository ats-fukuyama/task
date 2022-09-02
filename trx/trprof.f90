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
      REAL(rkind) :: PROF,qsurf,qaxis

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
      REAL(rkind)   :: ANEAVE, ANI, ANZ, DILUTE, TE, TRZEC,TRZEFE
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
      CALL TRZEFF

    END SUBROUTINE tr_prof_impurity

!     *** CALCULATE PROFILE OF AJ(R) ***

    SUBROUTINE tr_prof_current

      USE trcomm
      USE libitp
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind), DIMENSION(NRMAX) :: DSRHO
      REAL(rkind) :: FACT,FACTOR0,FACTORM,FACTORP,PROF
      REAL(rkind) :: SUML

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
         QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX) &
                      /(4.D0*PI**2*RDPVRHOG(1:NRMAX))
         Q0=(20.D0*QP(1)-23.D0*QP(2)+8.D0*QP(3))/5.D0

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

      USE TRCOMM, ONLY : NRAMAX, NSM, NSTM, PNSS, PNSSA, PTS, PTSA, RHOA, RN, RT, PNSSO,PTSO,PNSSAO,PTSAO
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSW
      INTEGER :: NS

      IF(RHOA.EQ.1.D0) RETURN

         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS(NS)
               PTSO (NS)=PTS (NS)

               PNSS (NS)=PNSSAO(NS)
               PTS  (NS)=PTSAO (NS)
            ENDDO
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

      IF(RHOA.NE.1.D0) THEN
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
