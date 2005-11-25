C
C     ***********************************************************
C
C           SET INITIAL PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRPROF
C
      INCLUDE 'trcomm.inc'
      DIMENSION DSRHO(NRM)
      COMMON /TMSLC2/ NTAMAX
C
      T     = 0.D0
      TPRE  = 0.D0
      TST   = 0.D0
      VSEC  = 0.D0
      NGR   = 0
      NGT   = 0
      NGST  = 0
      RIP   = RIPS
C
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
      ENDIF
C
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
C
         IF(MDLUF.EQ.1) THEN ! *** MDLUF ***
            IF(MDNI.EQ.0) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,2
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ELSE
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = RNU(1,NR,3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
C               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
C                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
C                  DO NS=1,3
C                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
C                  ENDDO
C               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
C                  RT(NR,1) = RTU(1,NR,1)
C                  RT(NR,2) = RTU(1,NR,2)
C                  RT(NR,3) = RTU(1,NR,3)
C               ENDIF
               RT(NR,1) = RTU(1,NR,1)
               RT(NR,2) = RTU(1,NR,2)
               RT(NR,3) = RTU(1,NR,3)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ENDIF
C
            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0
C
            PBM(NR)   = PBMU(1,NR)
            RNF(NR,1) = RNFU(1,NR)
            WROT(NR)  = WROTU(1,NR)
            VTOR(NR)  = WROTU(1,NR)*RMJRHOU(1,NR)
         ELSEIF(MDLUF.EQ.2) THEN ! *** MDLUF ***
            IF(MDNI.EQ.0) THEN !!!
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  DO NS=1,2
                     RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RN(NR,1) = RNU(1,NR,1)
                  RN(NR,2) = RNU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,3) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF
     &                    +RNU(1,NRMAX,2)
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF
     &                    +RNU(1,NRMAX,2)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,2
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               RT(NR,3) = RTU(1,NR,2)
               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,2
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ENDIF
            ELSE !!!
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  DO NS=1,3
                     RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  DO NS=1,3
                     RN(NR,NS) = RNU(1,NR,NS)
                  ENDDO
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF
     &                    +RNU(1,NRMAX,2)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,3
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  DO NS=1,3
                     RT(NR,NS) = RTU(1,NR,NS)
                  ENDDO
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               DO NS=1,3
                  RN(NR,NS) = RNU(1,NR,NS)
               ENDDO
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               DO NS=1,3
                  RN(NR,NS) = RNU(1,NR,NS)
               ENDDO
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,3
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  DO NS=1,3
                     RT(NR,NS) = RTU(1,NR,NS)
                  ENDDO
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ENDIF
            ENDIF !!!
C
            PEX(NR,1)=PNBU(1,NR,1)
            PEX(NR,2)=PNBU(1,NR,2)
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0
C
            SEX(NR,1)=SNBU(1,NR,1)+SWLU(1,NR)/PZ(2)
            SEX(NR,2)=SNBU(1,NR,2)+SWLU(1,NR)
            SEX(NR,3)=0.D0
            SEX(NR,4)=0.D0
            RNF(NR,1)=RNFU(1,NR)
            WROT(NR) =WROTU(1,NR)
            VTOR(NR) =WROTU(1,NR)*RMJRHOU(1,NR)
         ELSEIF(MDLUF.EQ.3) THEN ! *** MDLUF ***
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            RN(NR,1) = RNU(1,NR,1)
            RN(NR,2) = RNU(1,NR,2)
            RN(NR,3) = RNU(1,NR,3)
            RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
            RT(NR,1) = RTU(1,NR,1)
            RT(NR,2) = RTU(1,NR,2)
            RT(NR,3) = RTU(1,NR,3)
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                 +RTU(1,NRMAX,2)
C
c$$$            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
c$$$            DO NS=1,NSM
c$$$               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
c$$$            ENDDO
c$$$C
c$$$            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
c$$$            DO NS=1,NSM
c$$$               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
c$$$            ENDDO
C
            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0
            WROT(NR)  = WROTU(1,NR)
            VTOR(NR)  = WROTU(1,NR)*RMJRHOU(1,NR)
         ELSE ! *** MDLUF ***
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            DO NS=1,NSM
               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
            ENDDO
C
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            DO NS=1,NSM
               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
            ENDDO
C
            DO NS=1,NSM
               PEX(NR,NS) = 0.D0
            ENDDO
         ENDIF
C
         IF(MDLEQ0.EQ.1) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1)**PROFU2
            RN(NR,7) = (PN(7)-PNS(7))*PROF+PNS(7)
            RN(NR,8) = (PN(8)-PNS(8))*PROF+PNS(8)
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         ENDIF
C
         DO NF=1,NFM
            RW(NR,NF) = 0.D0
         ENDDO
C
         SUMPBM=SUMPBM+PBM(NR)
      ENDDO
      CALL PLDATA_SETR(RG,RM)
      CALL TR_EDGE_DETERMINER(1)
      CALL TR_EDGE_SELECTOR(1)
C
C     *** CALCULATE GEOMETRIC FACTOR ***
C
      CALL TRSTGF
      CALL TRGFRG
C
C     *** CALCULATE PZC,PZFE ***
C
      CALL TRZEFF
C
C     *** CALCULATE ANEAVE ***
C
      ANESUM=0.D0
      DO NR=1,NRMAX
         ANESUM=ANESUM+RN(NR,1)*RM(NR)
      ENDDO 
      ANEAVE=ANESUM*2.D0*DR
C
C     *** CALCULATE IMPURITY DENSITY
C                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***
C
      IF(MDLUF.NE.3) THEN
         DO NR=1,NRMAX
            ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC
     &               *1.D-2*RN(NR,1)
            ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE
     &               *1.D-2*RN(NR,1)
            ANI = 0.D0
            DO NS=2,NSM
               ANI=ANI+PZ(NS)*RN(NR,NS)
            ENDDO
            ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
            DILUTE = 1.D0-ANZ/ANI
            DO NS=2,NSM
               RN(NR,NS) = RN(NR,NS)*DILUTE
            ENDDO
         ENDDO
         PNSS(1)=PNS(1)
         DO NS=2,NSM
            PNSS(NS)=PNS(NS)*DILUTE
         ENDDO
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1)=PNSA(1)
            DO NS=2,NSM
               PNSSA(NS)=PNSA(NS)*DILUTE
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ELSE
         DO NS=1,NSM
            PNSS(NS)=PNS(NS)
         ENDDO
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            DO NS=1,NSM
               PNSSA(NS)=PNSA(NS)
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ENDIF
C
C     *** CALCULATE PROFILE OF AJ(R) ***
C
C     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***
C
      IF(MDLUF.EQ.1) THEN
         IF(MDLJQ.EQ.0) THEN ! *** MDLJQ ***
         NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         DO NR=1,NRMAX
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         ENDDO
         ELSE ! *** MDLJQ ***
         DO NR=1,NRMAX
            QP(NR) =QPU(1,NR)
            RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
         NR=1
            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
            AJOH(NR)=AJ(NR)
         DO NR=2,NRMAX
            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
            AJOH(NR)=AJ(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         ENDIF ! *** MDLJQ ***
      ELSEIF(MDLUF.EQ.2) THEN
         IF(MDLJQ.EQ.0) THEN  ! *** MDLJQ ***
            NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
c$$$         DO NR=1,NRMAX
c$$$            AJ(NR)=AJU(1,NR)
c$$$            AJNB(NR)=AJNBU(1,NR)
c$$$            AJOH(NR)=AJ(NR)-AJNB(NR)
c$$$            CALL TRSUMD(AJ,DVRHO,NR,RTMP)
c$$$            RDP(NR)=RMU0/(RR*DVRHOG(NR)*ABRHOG(NR))*RTMP*DR
c$$$            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
c$$$         ENDDO
C
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         DO NR=1,NRMAX
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         ENDDO
C
         ELSE ! *** MDLJQ ***
C
         DO NR=1,NRMAX
            QP(NR) =QPU(1,NR)
            RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
            NR=1
               FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
               AJOH(NR)=AJ(NR)
            DO NR=2,NRMAX
               FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
               FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
               AJOH(NR)=AJ(NR)
            ENDDO
            NR=1
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
               AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
            DO NR=2,NRMAX
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
               AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
            ENDDO
         ENDIF ! *** MDLJQ ***
         RIPA=DVRHOG(NRAMAX)*ABRHOG(NRAMAX)*RDP(NRAMAX)*1.D-6
     &       /(2.D0*PI*RMU0)
      ELSEIF(MDLUF.EQ.3) THEN
         DO NR=1,NRMAX
            AJOH(NR)=AJU(1,NR)
            AJ(NR)  =AJU(1,NR)
         ENDDO
C
         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
         ENDDO
         DO NR=1,NRMAX
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         DO NR=1,NRMAX
            AJOH(NR)=FACT*AJOH(NR)
            AJ(NR)  =AJOH(NR)
            BP(NR)  =FACT*BP(NR)
            QP(NR)  =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &              /(4.D0*PI**2*RDP(NR))
         ENDDO
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
C
         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
C
         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         DO NR=1,NRMAX
            RDP(NR)=FACT*RDP(NR)
            AJOH(NR)=FACT*AJOH(NR)
            AJ(NR)  =AJOH(NR)
            BP(NR)  =FACT*BP(NR)
            QP(NR)  =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &              /(4.D0*PI**2*RDP(NR))
         ENDDO
      ENDIF
      Q0=(4.D0*QP(1)-QP(2))/3.D0
C
C     calculate plasma current inside the calucated region (rho <= rhoa)
C     necessary for MDLEQB = 1 and MDLUF /= 0
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         VOL=0.D0
         DO NR=1,NRAMAX
            VOL=VOL+DVRHO(NR)*DR
            DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RMJRHO(NR))
         ENDDO
         CALL TRSUMD(AJ   ,DSRHO,NRAMAX,AJTSUM )
         RIPA=AJTSUM*DR/1.D6
      ENDIF
C
C     *** THIS MODEL ASSUMES CONSTANT EZ ***
C
      IF(PROFJ1.LE.0.D0.OR.MDNCLS.EQ.1) THEN
         CALL TRZEFF
         CALL TRCFET
         IF(PROFJ1.GT.0.D0.AND.MDNCLS.EQ.1) GOTO 2000
C
         DO NR=1,NRMAX
            AJOH(NR)=1.D0/ETA(NR)
            AJ(NR)  =1.D0/ETA(NR)
         ENDDO
C
         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         DO NR=1,NRMAX
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         DO NR=1,NRMAX
            RDP(NR)  =FACT*RDP(NR)
            AJOH(NR) =FACT*AJOH(NR)
            AJ(NR)   =AJOH(NR)
            AJTOR(NR)=FACT*AJTOR(NR)
            BP(NR)   =FACT*BP(NR)
            QP(NR)   =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &               /(4.D0*PI**2*RDP(NR))
         ENDDO
      ENDIF
 2000 CONTINUE
      SUM=0.D0
      DO NR=1,NRMAX
         SUM=SUM+RDP(NR)*DR
         RPSI(NR)=SUM
         BPRHO(NR)=BP(NR)
      ENDDO
C
      IF(MODELG.EQ.9) THEN
C     *** Give initial profiles to TASK/EQ ***
         CALL TREQIN(RR,RA,RKAP,RDLT,BB,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN1: IERR=',IERR
C     *** Contorol output display from TASK/EQ ***
CCC         CALL EQPARM(2,'EPSEQ=1.D-4',IERR)
CCC         CALL EQPARM(2,'NPRINT=1',IERR)
         CALL EQPARM(2,'NPRINT=0',IERR)
         CALL EQPARM(2,'NPSMAX=100',IERR)
C
         CALL TRSETG(0,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX TRPROF INITIAL EQUILIBRIUM FAILED : IERR = ',
     &                  IERR
         ENDIF
      ENDIF
C
      GRG(1)=0.0
      DO NR=1,NRMAX
         GRM(NR)  =GUCLIP(RM(NR))
         GRG(NR+1)=GUCLIP(RG(NR))
      ENDDO
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      RETURN
      END
C
C     ***********************************************************
C
C           SET GEOMETRICAL FACTOR VIA TASK/EQ
C
C     ***********************************************************
C
      SUBROUTINE TRSETG(ID,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION DUMMY(NRM)!,DSRHO(NRM)
C
      IF(IREAD.NE.0) THEN
C     *** Give initial profiles to TASK/EQ ***
         CALL TREQIN(RR,RA,RKAP,RDLT,BB,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN1: IERR=',IERR
C     *** Control output display from TASK/EQ ***
CCC         CALL EQPARM(2,'EPSEQ=1.D-4',IERR)
CCC         CALL EQPARM(2,'NPRINT=1',IERR)
         CALL EQPARM(2,'NPRINT=0',IERR)
         CALL EQPARM(2,'NPSMAX=100',IERR)
         IREAD=0
      ENDIF
C
      DO NR=1,NRMAX
         PRHO(NR)=0.D0
         TRHO(NR)=0.D0
         DO NS=1,NSM
            PRHO(NR)=PRHO(NR)+RN(NR,NS)*RT(NR,NS)*1.D14*RKEV
         ENDDO
         DO NF=1,NFM
            PRHO(NR)=PRHO(NR)+RW(NR,NF)*1.D14*RKEV
         ENDDO
         TRHO(NR)=0.D0
         DO NS=2,NSM
            TRHO(NR)=TRHO(NR)+RT(NR,NS)*RN(NR,NS)/RN(NR,1)
         ENDDO
         HJRHO(NR)=AJ(NR)
         VTRHO(NR)=0.D0
         RHOTR(NR)=RM(NR)
      ENDDO
C
      IF(ID.EQ.0) THEN
         ICONT=0
         SRIP=RIPS
      ELSE
         ICONT=1
         SRIP=0.D0
      ENDIF
C
C     *** Excute TASK/EQ ***
      CALL TREQEX (NRMAX,RM,PRHO,HJRHO,VTRHO,TRHO,SRIP,ICONT,
     &             RSA,RDPA,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX TREQEX: IERR=',IERR
         WRITE(6,*) '&TRDATA'
         WRITE(6,*) '  RR=',RR
         WRITE(6,*) '  RA=',RA
         WRITE(6,*) '  RKAP=',RKAP
         WRITE(6,*) '  RDLT=',RDLT
         WRITE(6,*) '  BB=',BB
         WRITE(6,*) '  RIP=',RIPS
         WRITE(6,*) '  PRHO=',(PRHO(NR),NR=1,NRMAX)
         WRITE(6,*) '  HJRHO=',(HJRHO(NR),NR=1,NRMAX)
         WRITE(6,*) '  VTRHO=',(VTRHO(NR),NR=1,NRMAX)
         WRITE(6,*) '  TRHO=',(TRHO(NR),NR=1,NRMAX)
         WRITE(6,*) '&END'
         RETURN
      ENDIF
C     
C     Initially(NT=0), current density as an input quantity is modified
C     in order to keep total plasma current constant.
C
      IF(ID.EQ.0) THEN
         DO NR=1,NRMAX
            AJ(NR)=HJRHO(NR)
         ENDDO
      ENDIF
C     
C     *** Provide geometric quantities on half mesh ***
      CALL TREQGET(NRMAX,RM,
     &             DUMMY,TTRHO,DVRHO,DUMMY,ABRHO,ARRHO,
     &             AR1RHO,AR2RHO,DUMMY,DUMMY,DUMMY,
     &             DUMMY,DUMMY,DUMMY,RKPRHO,
     &             IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQGET1: IERR=',IERR
C
C     *** Provide geometric quantities on grid mesh ***
      CALL TREQGET(NRMAX,RG,
     &             QRHO,TTRHOG,DVRHOG,DUMMY,ABRHOG,ARRHOG,
     &             AR1RHOG,AR2RHOG,ABB2RHOG,AIB2RHOG,ARHBRHOG,
     &             EPSRHO,RMNRHO,RMNRHO,RKPRHOG,
     &             IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TREQGET2: IERR=',IERR
C
C     *** Adjust system of unit ***
      RDPA=RDPA/(2.D0*PI)
      DO NR=1,NRMAX
         TTRHO(NR)   =TTRHO(NR)/(2.D0*PI)
         ARRHO(NR)   =ARRHO(NR)/RR**2
         TTRHOG(NR)  =TTRHOG(NR)/(2.D0*PI)
         ARRHOG(NR)  =ARRHOG(NR)/RR**2
         ABB2RHOG(NR)=ABB2RHOG(NR)*BB**2
         AIB2RHOG(NR)=AIB2RHOG(NR)/BB**2
      ENDDO
C
c$$$      DO NR=1,NRMAX
c$$$         RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
c$$$     &          /(4.D0*PI**2*QRHO(NR))
c$$$         QP(NR)=QRHO(NR)
c$$$      ENDDO
C
c$$$      DO NR=NRMAX-1,1,-1
c$$$         RDP(NR) =( DVRHOG(NR+1)*ABRHOG(NR+1)/TTRHOG(NR+1)*RDP(NR+1)
c$$$     &             -RMU0*BB*DVRHO(NR+1)*DR/TTRHO(NR+1)**2*AJ(NR+1))
c$$$     &            /(DVRHOG(NR)*ABRHOG(NR)/TTRHOG(NR))
c$$$      ENDDO
C
      NR=1
         RDP(NR) = RMU0*BB*DVRHO(NR)*DR/TTRHO(NR)**2*AJ(NR)
     &           /(DVRHOG(NR)*ABRHOG(NR)/TTRHOG(NR))
C
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
         AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
      DO NR=2,NRMAX
         RDP(NR) =( DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)*RDP(NR-1)
     &             +RMU0*BB*DVRHO(NR)*DR/TTRHO(NR)**2*AJ(NR))
     &            /(DVRHOG(NR)*ABRHOG(NR)/TTRHOG(NR))
C
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
         AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
      ENDDO
      DO NR=1,NRMAX
C         write(6,*) RM(NR),QP(NR)
         QP(NR)=2.D0*PI*BB*(RSA/SQRT(2.D0*PI))**2*RG(NR)/RDP(NR)
c$$$         DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RR)
      ENDDO
c$$$      CALL TRSUMD(AJ   ,DSRHO,NRMAX,AJTSUM)
c$$$      AJT   = AJTSUM*DR/1.D6
c$$$      write(6,*) "Parallel current =",AJT
c$$$      CALL TRSUMD(AJTOR,DSRHO,NRMAX,AJTTSUM)
c$$$      AJTT   = AJTTSUM*DR/1.D6
c$$$      write(6,*) "Toroidal current =",AJTT
c$$$      NR=NRMAX
c$$$      write(6,*) "Current from RDP =",
c$$$     &     DVRHOG(NR)*ABRHOG(NR)*RDP(NR)/(2.D0*PI*RMU0)*1.D-6
c$$$      write(6,*) "Current from RDPA=",
c$$$     &     DVRHOG(NR)*ABRHOG(NR)*RDPA/(2.D0*PI*RMU0)*1.D-6
C
      RETURN
      END
C
C     ********************************************************
C
C           RADIAL INTEGRATION ONLY FOR J CONVERGENCE
C
C     ********************************************************
C
      SUBROUTINE TRSUMJ(A,B,NMAX,SUM)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION A(NMAX),B(NMAX)
C
      SUM=0.D0
      DO N=1,NMAX
         SUM=SUM+A(N)**2*B(N)
      ENDDO
      RETURN
      END
C
C     ***********************************************************
C
C           SET GEOMETRIC FACTOR AT HALF MESH
C
C     ***********************************************************
C
      SUBROUTINE TRSTGF
C
      INCLUDE 'trcomm.inc'
C
      RKAPS=SQRT(RKAP)
      IF(MDLUF.NE.0) THEN
         IF(MODELG.EQ.2.OR.MODELG.EQ.3.OR.MODELG.EQ.9) THEN
            DO NR=1,NRMAX
               TTRHO(NR)=TTRHOU(1,NR)
               DVRHO(NR)=DVRHOU(1,NR)
               ABRHO(NR)=ABRHOU(1,NR)
               ARRHO(NR)=ARRHOU(1,NR)
               AR1RHO(NR)=AR1RHOU(1,NR)
               AR2RHO(NR)=AR2RHOU(1,NR)
               RMJRHO(NR)=RMJRHOU(1,NR)
               RMNRHO(NR)=RMNRHOU(1,NR)
               RKPRHO(NR)=RKPRHOU(1,NR)
               IF(MDPHIA.EQ.0) THEN
C     define rho_a from phi_a data
                  RHO_A=SQRT(PHIA/(PI*BB))
                  RJCB(NR)=1.D0/RHO_A
                  RHOM(NR)=RM(NR)*RHO_A
                  RHOG(NR)=RG(NR)*RHO_A
               ELSE
C     define rho_a from volume data at ideal LCFS
                  RHO_A=SQRT(VOLAU(1)/(2.D0*PI**2*RMJRHOU(1,NRMAX)))
                  RJCB(NR)=1.D0/RHO_A
                  RHOM(NR)=RM(NR)/RJCB(NR)
                  RHOG(NR)=RG(NR)/RJCB(NR)
               ENDIF
               EPSRHO(NR)=RMNRHO(NR)/RMJRHO(NR)
            ENDDO
            CALL FLUX
         ELSEIF(MODELG.EQ.5) THEN
C            CALL INITIAL_EQDSK(EPSRHO,TTRHO,DVRHO,ABRHO,ARRHO,
C     &                         AR1RHO,AR2RHO,RJCB,RMJRHO,RMNRHO,RKPRHO,
C     &                         NRMAX,NRM)
         ENDIF
      ELSE
         IF(MODELG.EQ.2.OR.MODELG.EQ.3.OR.MODELG.EQ.9) THEN
            DO NR=1,NRMAX
               BPRHO(NR)=BP(NR)
               QRHO(NR)=QP(NR)
C
               TTRHO(NR)=BB*RR
               DVRHO(NR)=2.D0*PI*RKAP*RA*RA*2.D0*PI*RR*RM(NR)
               ABRHO(NR)=1.D0/(RKAPS*RA*RR)**2
               ARRHO(NR)=1.D0/RR**2
               AR1RHO(NR)=1.D0/(RKAPS*RA)
               AR2RHO(NR)=1.D0/(RKAPS*RA)**2
               RMJRHO(NR)=RR
               RMNRHO(NR)=RA*RG(NR)
               RKPRHO(NR)=RKAP
               RJCB(NR)=1.D0/(RKAPS*RA)
               RHOM(NR)=RM(NR)/RJCB(NR)
               RHOG(NR)=RG(NR)/RJCB(NR)
               EPSRHO(NR)=RMNRHO(NR)/RMJRHO(NR)
            ENDDO
         ELSEIF(MODELG.EQ.5) THEN
C            CALL INITIAL_EQDSK(EPSRHO,TTRHO,DVRHO,ABRHO,ARRHO,
C     &                         AR1RHO,AR2RHO,RJCB,RMJRHO,RMNRHO,RKPRHO,
C     &                         NRMAX,NRM)
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           GEOMETRIC QUANTITIES AT GRID MESH
C
C     ***********************************************************
C
      SUBROUTINE TRGFRG
C
      INCLUDE 'trcomm.inc'
C
      DO NR=1,NRMAX-1
         AR1RHOG(NR)=0.5D0*(AR1RHO(NR)+AR1RHO(NR+1))
         AR2RHOG(NR)=0.5D0*(AR2RHO(NR)+AR2RHO(NR+1))
         RKPRHOG(NR)=0.5D0*(RKPRHO(NR)+RKPRHO(NR+1))
         TTRHOG (NR)=0.5D0*(TTRHO (NR)+TTRHO (NR+1))
         DVRHOG (NR)=0.5D0*(DVRHO (NR)+DVRHO (NR+1))
         ARRHOG (NR)=0.5D0*(ARRHO (NR)+ARRHO (NR+1))
         ABRHOG (NR)=0.5D0*(ABRHO (NR)+ABRHO (NR+1))
C
         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
         ARHBRHOG(NR)=AR2RHOG(NR)*AIB2RHOG(NR)
      ENDDO
      NR=NRMAX
         RGL=RG(NR)
         RML=RM(NR)
         RML1=RM(NR-1)
C
         AR1RHOG(NR)=FEDG(RGL,RML,RML1,AR1RHO(NR),AR1RHO(NR-1))
         AR2RHOG(NR)=FEDG(RGL,RML,RML1,AR2RHO(NR),AR2RHO(NR-1))
         RKPRHOG(NR)=FEDG(RGL,RML,RML1,RKPRHO(NR),RKPRHO(NR-1))
         TTRHOG (NR)=FEDG(RGL,RML,RML1,TTRHO (NR),TTRHO (NR-1))
         DVRHOG (NR)=FEDG(RGL,RML,RML1,DVRHO (NR),DVRHO (NR-1))
         ARRHOG (NR)=FEDG(RGL,RML,RML1,ARRHO (NR),ARRHO (NR-1))
         ABRHOG (NR)=FEDG(RGL,RML,RML1,ABRHO (NR),ABRHO (NR-1))
C
         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
         ARHBRHOG(NR)=AR2RHOG(NR)*AIB2RHOG(NR)
C
      RETURN
      END
C
C     ***********************************************************
C
C           EDGE VALUE SELECTOR
C
C     ***********************************************************
C
      SUBROUTINE TR_EDGE_SELECTOR(NSW)
C
      INCLUDE 'trcomm.inc'
      DIMENSION PNSSO(NSTM),PTSO(NSTM),PNSSAO(NSTM),PTSAO(NSTM)
      SAVE PNSSO,PTSO,PNSSAO,PTSAO
C
      IF(RHOA.EQ.1.D0) RETURN
C
      IF(MDLUF.EQ.0) THEN
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS(NS)
               PTSO (NS)=PTS (NS)
C
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
C
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
C     
      ENTRY TR_EDGE_DETERMINER(NSW)
C
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
C
      RETURN
      END
