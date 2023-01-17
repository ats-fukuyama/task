!     $Id: fpcalwr.f90,v 1.7 2013/01/22 04:31:15 fukuyama Exp $
!
! ************************************************************
!
!            CALCULATION OF DW (BOUNCE AVERAGED,RAY)
!
! ************************************************************

      MODULE fpcalwr

      USE fpcomm
      USE libbes,ONLY: bessjn

      contains

!--------------------------------------

      SUBROUTINE FP_CALWR(NSA)

      USE fpwrin
      USE plprof,ONLY: pl_getRZ
      IMPLICIT NONE
      INTEGER:: NSA, NSBA, NITM, NP, NTH, NRDO, NR, NAV, NS
      REAL(rkind):: FACT, DELH, ETAL, THETAL, RRAVE, RSAVE
      REAL(rkind):: X,Y,Z,RL,ZL,RXB,RYB,RZB,RRLB,RZLB
      INTEGER:: NRAY,NIT,NITMX,MINNB1,MINNB2
      REAL(rkind):: DLAMN1,RXMIN1,RYMIN1,RZMIN1,RKXMN1,RKYMN1,RKZMN1,RBMIN1
      REAL(rkind):: DLAMN2,RXMIN2,RYMIN2,RZMIN2,RKXMN2,RKYMN2,RKZMN2,RBMIN2
      COMPLEX(rkind):: CEXMN1,CEYMN1,CEZMN1,CEXMN2,CEYMN2,CEZMN2
      REAL(rkind):: RRLMN1,RRLMN2,RZLMN1,RZLMN2,DEL12,XA1,XA2,XLL2,A1,A2
      COMPLEX(rkind):: CEX,CEY,CEZ
      REAL(rkind):: RKX,RKY,RKZ,RADB,DELCR2,DELRB2,ARG
      REAL(rkind):: DWPPS,DWPTS,DWTPS,DWTTS

      REAL(rkind),DIMENSION(:,:),POINTER:: DLA !(0:NITMAXM,NRAYMAX)

! =============  CALCULATION OF DWPP AND DWPT  ===============

      IF(NRAYS_WR.EQ.0) THEN
         NRAYS=1
      ELSE
         NRAYS=NRAYS_WR
      END IF
      IF(NRAYE_WR.EQ.0.OR.NRAYE_WR.GT.NRAYMAX) THEN
         NRAYE=NRAYMAX
      ELSE
         NRAYE=NRAYE_WR
      END IF

      FACT=0.5D0

      ALLOCATE(DLA(0:NITMAXM,NRAYS:NRAYE))

      NS=NS_NSA(NSA)

      IF(MODELW(NS).EQ.1.OR.MODELW(NS).EQ.2) THEN
         DO NRDO=NRSTART,NREND
            NR=NRDO
            DO NTH=1,NTHMAX
               DELH=4.D0*ETAM(NTH,NR)/NAVMAX

               DO NAV=1,NAVMAX
                  ETAL=DELH*(NAV-0.5D0)-2.D0*ETAM(NTH,NR)
                  CALL pl_getRZ(RM(NR),ETAL,RL,ZL)

                  DO NRAY=NRAYS,NRAYE
                     NITMX=NITMAX(NRAY)
                     RF_WR=RAYIN(1,NRAY)

                     DO NIT=0,NITMX
                        RXB=RXS(NIT,NRAY)
                        RYB=RYS(NIT,NRAY)
                        RZB=RZS(NIT,NRAY)
                        RRLB=SQRT(RXB**2+RYB**2)
                        RZLB=RZB
                        DLA(NIT,NRAY)=SQRT((RRLB-RL)**2+(RZLB-ZL)**2)
!                        write(6,'(A,2I5,1P3E12.4)') &
!                          'point 2',NRAY,NIT,RRLB,RZLB,DLA(NIT,NRAY)
                     ENDDO

                     MINNB1=0
                     DO NIT=0,NITMX-1
                        IF(DLA(NIT+1,NRAY).LT.DLA(NIT,NRAY))THEN
                           MINNB1=NIT+1
                        ENDIF
                     ENDDO

                     IF(MINNB1.EQ.0) GOTO 1
                     IF(MINNB1.EQ.NITMX) GOTO 1

                     IF(DLA(MINNB1-1,NRAY).LT.DLA(MINNB1+1,NRAY))THEN
                        MINNB2=MINNB1-1
                     ELSE
                        MINNB2=MINNB1+1
                     ENDIF

!                     write(6,'(A,4I5,1P4E12.4)') &
!                    'point 3',NRAY,NIT,MINNB1,MINNB2, &
!                     DLA(MINNB1,NRAY),DLA(MINNB2,NRAY)

                     DLAMN1=   DLA(MINNB1,NRAY)
                     RXMIN1=   RXS(MINNB1,NRAY)
                     RYMIN1=   RYS(MINNB1,NRAY)
                     RZMIN1=   RZS(MINNB1,NRAY)
                     CEXMN1=  CEXS(MINNB1,NRAY)
                     CEYMN1=  CEYS(MINNB1,NRAY)
                     CEZMN1=  CEZS(MINNB1,NRAY)
                     RKXMN1=  RKXS(MINNB1,NRAY)
                     RKYMN1=  RKYS(MINNB1,NRAY)
                     RKZMN1=  RKZS(MINNB1,NRAY)
                     RBMIN1=RAYRB1(MINNB1,NRAY)

                     DLAMN2=   DLA(MINNB2,NRAY)
                     RXMIN2=   RXS(MINNB2,NRAY)
                     RYMIN2=   RYS(MINNB2,NRAY)
                     RZMIN2=   RZS(MINNB2,NRAY)
                     CEXMN2=  CEXS(MINNB2,NRAY)
                     CEYMN2=  CEYS(MINNB2,NRAY)
                     CEZMN2=  CEZS(MINNB2,NRAY)
                     RKXMN2=  RKXS(MINNB2,NRAY)
                     RKYMN2=  RKYS(MINNB2,NRAY)
                     RKZMN2=  RKZS(MINNB2,NRAY)
                     RBMIN2=RAYRB1(MINNB2,NRAY)
               
                     RRLMN1=SQRT(RXMIN1**2+RYMIN1**2)
                     RZLMN1=RZMIN1
                     RRLMN2=SQRT(RXMIN2**2+RYMIN2**2)
                     RZLMN2=RZMIN2
                     DEL12=SQRT((RRLMN1-RRLMN2)**2+(RZLMN1-RZLMN2)**2)

                     XA1=(DLAMN1**2-DLAMN2**2+DEL12**2)/(2.D0*DEL12)
                     XA2=(DLAMN2**2-DLAMN1**2+DEL12**2)/(2.D0*DEL12)             
                     XLL2=DLAMN1**2-XA1**2
                     A1=XA1/(XA1+XA2)
                     A2=XA2/(XA1+XA2)

                     CEX =A1*CEXMN2+A2*CEXMN1
                     CEY =A1*CEYMN2+A2*CEYMN1
                     CEZ =A1*CEZMN2+A2*CEZMN1
                     RKX =A1*RKXMN2+A2*RKXMN1
                     RKY =A1*RKYMN2+A2*RKYMN1
                     RKZ =A1*RKZMN2+A2*RKZMN1
                     RXB =A1*RXMIN2+A2*RXMIN1
                     RYB =A1*RYMIN2+A2*RYMIN1
                     RZB =A1*RZMIN2+A2*RZMIN1
                     RADB=A1*RBMIN2+A2*RBMIN1
                     IF(RADB.EQ.0.D0) THEN
                        DELRB2=DELY_WR**2 ! for ray tracing
                     ELSE
                        DELRB2=RADB**2   ! for beam tracing
                     END IF
                     DELCR2=XLL2
                     ARG=DELCR2/DELRB2

                     ARGB (NR,NTH,NAV,NRAY)=ARG
                     CEB(1,NR,NTH,NAV,NRAY)=CEX
                     CEB(2,NR,NTH,NAV,NRAY)=CEY
                     CEB(3,NR,NTH,NAV,NRAY)=CEZ
                     RKB(1,NR,NTH,NAV,NRAY)=RKX
                     RKB(2,NR,NTH,NAV,NRAY)=RKY
                     RKB(3,NR,NTH,NAV,NRAY)=RKZ
                     RBB(1,NR,NTH,NAV,NRAY)=RXB
                     RBB(2,NR,NTH,NAV,NRAY)=RYB
                     RBB(3,NR,NTH,NAV,NRAY)=RZB
1                    CONTINUE
            
!            IF(IDEBUG.EQ.1) THEN
!               WRITE(6,'(5I5)') NR,NTH,NAV,NRAY,NIT
!               WRITE(6,'(1P4E12.4)') RRLMN1,RRLMN2,RZLMN1,RZLMN2
!               WRITE(6,'(1P6E12.4)') CEX,CEY,CEZ
!               WRITE(6,'(1P3E12.4)') DELCR2,DELRB2,ARG
!            ENDIF

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO NRDO=NRSTART,NREND
         NR=NRDO
         DO NTH=1,NTHMAX
            IF(NTH.NE.ITL(NR).AND.NTH.NE.ITU(NR)) THEN
!               DO NP=1,NPMAX+1
               DO NP=NPSTART,NPENDWG
                  CALL FPDWAV(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP,NS), &
                              NR,NTH,DWPPS,DWPTS,DWTPS,DWTTS,NSA)
                  DWWRPP(NTH,NP,NR,NSA)=DWPPS
                  DWWRPT(NTH,NP,NR,NSA)=DWPTS
!                  IF(ABS(DWPPS).GT.1.D-12) THEN
!                     WRITE(6,'(A,3I5,1PE12.4)') &
!                          'NR,NTH,NP,DWPPS=',NR,NTH,NP,DWPPS
!                  END IF
               ENDDO
            ENDIF
         ENDDO

         IF(MODELA.EQ.1) THEN
!            DO NP=1,NPMAX+1
            DO NP=NPSTART,NPENDWG
               DO NTH=ITL(NR)+1,NTHMAX/2
                  DWWRPP(NTH,NP,NR,NSA)  =(DWWRPP(NTH,NP,NR,NSA) &
                                          +DWWRPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWWRPT(NTH,NP,NR,NSA)  =(DWWRPT(NTH,NP,NR,NSA) &
                                          +DWWRPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
                  DWWRPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWWRPP(NTH,NP,NR,NSA)
                  DWWRPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWWRPT(NTH,NP,NR,NSA)
               ENDDO
               DWWRPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0             &
                    *( DWWRPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                      +DWWRPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                      +DWWRPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                      +DWWRPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR))

               DWWRPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0             &
                    *( DWWRPT(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                      +DWWRPT(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                      +DWWRPT(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                      +DWWRPT(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR))
               DWWRPP(ITU(NR),NP,NR,NSA)  =DWWRPP(ITL(NR),NP,NR,NSA)
               DWWRPT(ITU(NR),NP,NR,NSA)  =DWWRPT(ITL(NR),NP,NR,NSA)
            ENDDO
         ENDIF
      ENDDO

! =============  CALCULATION OF DWTP AND DWTT  ===============

      IF(MODELW(NS).EQ.1.OR.MODELW(NS).EQ.2) THEN
         DO  NRDO=NRSTART,NREND
            NR=NRDO
            DO  NTH=1,NTHMAX
               DELH=4.D0*ETAG(NTH,NR)/NAVMAX
               DO NAV=1,NAVMAX
                  ETAL=DELH*(NAV-0.5D0)-2.D0*ETAG(NTH,NR)
                  CALL pl_getRZ(RM(NR),ETAL,RL,ZL)

                  DO NRAY=NRAYS,NRAYE
                     NITMX=NITMAX(NRAY)
                     RF_WR=RAYIN(1,NRAY)

                     DO NIT=0,NITMX
                        RXB=RXS(NIT,NRAY)
                        RYB=RYS(NIT,NRAY)
                        RZB=RZS(NIT,NRAY)
                        RRLB=SQRT(RXB**2+RYB**2)
                        RZLB=RZB
                        DLA(NIT,NRAY)=SQRT((RRLB-RL)**2+(RZLB-ZL)**2)
                     ENDDO

                     MINNB1=0
                     DO NIT=0,NITMX-1
                        IF(DLA(NIT+1,NRAY).LT.DLA(NIT,NRAY))THEN
                           MINNB1=NIT+1
                        ENDIF
                     ENDDO

                     IF(MINNB1.EQ.0) GOTO 2
                     IF(MINNB1.EQ.NITMX) GOTO 2

                     IF(DLA(MINNB1-1,NRAY).LT.DLA(MINNB1+1,NRAY))THEN
                        MINNB2=MINNB1-1
                     ELSE
                        MINNB2=MINNB1+1
                     ENDIF

                     DLAMN1=   DLA(MINNB1,NRAY)
                     RXMIN1=   RXS(MINNB1,NRAY)
                     RYMIN1=   RYS(MINNB1,NRAY)
                     RZMIN1=   RZS(MINNB1,NRAY)
                     CEXMN1=  CEXS(MINNB1,NRAY)
                     CEYMN1=  CEYS(MINNB1,NRAY)
                     CEZMN1=  CEZS(MINNB1,NRAY)
                     RKXMN1=  RKXS(MINNB1,NRAY)
                     RKYMN1=  RKYS(MINNB1,NRAY)
                     RKZMN1=  RKZS(MINNB1,NRAY)
                     RBMIN1=RAYRB1(MINNB1,NRAY)

                     DLAMN2=   DLA(MINNB2,NRAY)
                     RXMIN2=   RXS(MINNB2,NRAY)
                     RYMIN2=   RYS(MINNB2,NRAY)
                     RZMIN2=   RZS(MINNB2,NRAY)
                     CEXMN2=  CEXS(MINNB2,NRAY)
                     CEYMN2=  CEYS(MINNB2,NRAY)
                     CEZMN2=  CEZS(MINNB2,NRAY)
                     RKXMN2=  RKXS(MINNB2,NRAY)
                     RKYMN2=  RKYS(MINNB2,NRAY)
                     RKZMN2=  RKZS(MINNB2,NRAY)
                     RBMIN2=RAYRB1(MINNB2,NRAY)

                     RRLMN1=SQRT(RXMIN1**2+RYMIN1**2)
                     RZLMN1=RZMIN1
                     RRLMN2=SQRT(RXMIN2**2+RYMIN2**2)
                     RZLMN2=RZMIN2
                     DEL12=SQRT((RRLMN1-RRLMN2)**2+(RZLMN1-RZLMN2)**2)

                     XA1=(DLAMN1**2-DLAMN2**2+DEL12**2)/(2.D0*DEL12)
                     XA2=(DLAMN2**2-DLAMN1**2+DEL12**2)/(2.D0*DEL12)
                     XLL2=DLAMN1**2-XA1**2
                     A1=XA1/(XA1+XA2)
                     A2=XA2/(XA1+XA2)

                     CEX =A1*CEXMN2+A2*CEXMN1
                     CEY =A1*CEYMN2+A2*CEYMN1
                     CEZ =A1*CEZMN2+A2*CEZMN1
                     RKX =A1*RKXMN2+A2*RKXMN1
                     RKY =A1*RKYMN2+A2*RKYMN1
                     RKZ =A1*RKZMN2+A2*RKZMN1
                     RXB =A1*RXMIN2+A2*RXMIN1
                     RYB =A1*RYMIN2+A2*RYMIN1
                     RZB =A1*RZMIN2+A2*RZMIN1
                     RADB=A1*RBMIN2+A2*RBMIN1
                     IF(RADB.EQ.0.D0) THEN
                        DELRB2=DELY_WR**2 ! for ray tracing
                     ELSE
                        DELRB2=RADB**2   ! for beam tracing
                     END IF
                     DELCR2=XLL2
                     ARG=DELCR2/DELRB2

                     ARGB (NR,NTH,NAV,NRAY)=ARG
                     CEB(1,NR,NTH,NAV,NRAY)=CEX
                     CEB(2,NR,NTH,NAV,NRAY)=CEY
                     CEB(3,NR,NTH,NAV,NRAY)=CEZ
                     RKB(1,NR,NTH,NAV,NRAY)=RKX
                     RKB(2,NR,NTH,NAV,NRAY)=RKY
                     RKB(3,NR,NTH,NAV,NRAY)=RKZ
                     RBB(1,NR,NTH,NAV,NRAY)=RXB
                     RBB(2,NR,NTH,NAV,NRAY)=RYB
                     RBB(3,NR,NTH,NAV,NRAY)=RZB
2                    CONTINUE
                  ENDDO

!            IF(IDEBUG.EQ.1) THEN
!               WRITE(6,'(3I3)') NR,NAV,NCR
!               WRITE(6,'(1P3E12.4)') X,Y,Z
!               WRITE(6,'(1P3E12.4)') RX,RY,RZ
!               WRITE(6,'(1P3E12.4)') DELR2,DELCR2,ARG
!            ENDIF

               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DO NRDO=NRSTART,NREND
         NR=NRDO
         DO NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM 
                  CALL FPDWAV(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP,NS), &
                              NR,NTH,DWPPS,DWPTS,DWTPS,DWTTS,NSA)
                  DWWRTP(NTH,NP,NR,NSA)=DWTPS
                  DWWRTT(NTH,NP,NR,NSA)=DWTTS
               ENDDO
            ELSE
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM 
                  DWWRTP(NTH,NP,NR,NSA)=0.D0
                  DWWRTT(NTH,NP,NR,NSA)=0.D0
               ENDDO
            ENDIF
         ENDDO

         IF(MODELA.EQ.1) THEN
            DO NTH=ITL(NR)+1,NTHMAX/2
               DO NP=NPSTARTW,NPENDWM 
!               DO NP=1,NPMAX
                  DWWRTP(NTH,NP,NR,NSA)=(DWWRTP(NTH,NP,NR,NSA) &
                                        -DWWRTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWWRTT(NTH,NP,NR,NSA)=(DWWRTT(NTH,NP,NR,NSA) &
                                        +DWWRTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
                  DWWRTP(NTHMAX-NTH+2,NP,NR,NSA)=-DWWRTP(NTH,NP,NR,NSA)
                  DWWRTT(NTHMAX-NTH+2,NP,NR,NSA)= DWWRTT(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      DEALLOCATE(DLA)
      RETURN
      END SUBROUTINE FP_CALWR

!***********************************************************************
!     Calculate DW averaged over magnetic surface for (p,theta0,r)
!***********************************************************************

      SUBROUTINE FPDWAV(ETA,RSIN,RCOS,P,NR,NTH,  &
                        DWPPS,DWPTS,DWTPS,DWTTS,NSA)

      USE fpwrin
      USE plprof,ONLY: pl_getRZ
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: ETA,RSIN,RCOS,P
      INTEGER,INTENT(IN):: NR,NTH,NSA
      REAL(rkind),INTENT(OUT):: DWPPS,DWPTS,DWTPS,DWTTS
      REAL(rkind):: DELH,ETAL,THETAL,RRAVE,RSAVE,X,Y,Z,RL,ZL
      REAL(rkind):: RX,RY,RZ,RLCR,ZLCR,DELR2,DELCR2,ARG,FACTOR
      REAL(rkind):: RKX,RKY,RKZ,DWPPL,DWPTL,DWTPL,DWTTL
      REAL(rkind):: PSIN,PCOS,PSI
      COMPLEX(rkind):: CEX,CEY,CEZ
      INTEGER:: NAV,NRAY,NS,NCR

      NS=NS_NSA(NSA)
      DELH=4.D0*ETA/NAVMAX

      DWPPS=0.D0
      DWPTS=0.D0
      DWTPS=0.D0
      DWTTS=0.D0

      DO NAV=1,NAVMAX
         ETAL=DELH*(NAV-0.5D0)-2.D0*ETA
         CALL pl_getRZ(RM(NR),ETAL,RL,ZL)

         DO NRAY=NRAYS,NRAYE
            RF_WR=RAYIN(1,NRAY)

            IF(MODELW(NS).EQ.1) THEN
               DO NCR=1,NCRMAX(NR,NRAY)
                  RX=RCR(1,NCR,NR,NRAY)
                  RY=RCR(2,NCR,NR,NRAY)
                  RZ=RCR(3,NCR,NR,NRAY)
                  RLCR=SQRT(RX**2+RY**2)
                  ZLCR=RZ
                  DELR2=(RL-RLCR)**2+(ZL-ZLCR)**2
                  DELCR2=DELY_WR**2
                  ARG=DELR2/DELCR2

                  IF(ARG.LT.15.D0) THEN
                     FACTOR=EXP(-ARG)
                     CALL FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI,NSA)
                     CEX=CECR(1,NCR,NR,NRAY)*FACTOR
                     CEY=CECR(2,NCR,NR,NRAY)*FACTOR
                     CEZ=CECR(3,NCR,NR,NRAY)*FACTOR
                     RKX=RKCR(1,NCR,NR,NRAY)
                     RKY=RKCR(2,NCR,NR,NRAY)
                     RKZ=RKCR(3,NCR,NR,NRAY)
                     CALL FPDWLL(P,PSIN,PCOS,                      &
                                 CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ, &
                                 DWPPL,DWPTL,DWTPL,DWTTL,NSA)
!                     IF(NRAY.EQ.2) &
!                     WRITE(28,'(A,4I5,1P3E12.4)') &
!                          'FPDWAV:',NR,NAV,NRAY,NCR,CEX,DWPPL
!                     IF(NRAY.EQ.2) &
!                     WRITE(28,'(A,1P6E12.4)') &
!                          '------:',RX,RY,RZ,RKX,RKY,RKZ
!                     IF(NRAY.EQ.1) &
!                     WRITE(29,'(A,4I5,1P3E12.4)') &
!                          'FPDWAV:',NR,NAV,NRAY,NCR,CEX,DWPPL
!                     IF(NRAY.EQ.1) &
!                     WRITE(29,'(A,1P6E12.4)') &
!                          '------:',RX,RY,RZ,RKX,RKY,RKZ
                     DWPPS=DWPPS+DWPPL*RCOS/PCOS
                     DWPTS=DWPTS+DWPTL          /SQRT(PSI)
                     DWTPS=DWTPS+DWTPL          /SQRT(PSI)
                     DWTTS=DWTTS+DWTTL*PCOS/RCOS/PSI
!                     WRITE(6,*) NR,NTH,DWPPS
!                     IF(IDEBUG.EQ.1) THEN
!                        WRITE(6,'(3I3)') NR,NAV,NCR
!                        WRITE(6,'(1P3E12.4)') X,Y,Z
!                        WRITE(6,'(1P3E12.4)') RX,RY,RZ
!                        WRITE(6,'(1P3E12.4)') DELR2,DELCR2,ARG
!                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

            IF(MODELW(NS).EQ.2) THEN
               ARG=ARGB(NR,NTH,NAV,NRAY)
               IF(ARG.GT.0.D0.AND.ARG.LT.15.D0) THEN
                  FACTOR= EXP(-ARG)
                  CALL FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI,NSA)
                  CEX=CEB(1,NR,NTH,NAV,NRAY)*FACTOR
                  CEY=CEB(2,NR,NTH,NAV,NRAY)*FACTOR
                  CEZ=CEB(3,NR,NTH,NAV,NRAY)*FACTOR
                  RKX=RKB(1,NR,NTH,NAV,NRAY)
                  RKY=RKB(2,NR,NTH,NAV,NRAY)
                  RKZ=RKB(3,NR,NTH,NAV,NRAY)
                  RX =RBB(1,NR,NTH,NAV,NRAY)
                  RY =RBB(2,NR,NTH,NAV,NRAY)
                  RZ =RBB(3,NR,NTH,NAV,NRAY)
                  CALL FPDWLL(P,PSIN,PCOS,               &
                              CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ, &
                              DWPPL,DWPTL,DWTPL,DWTTL,NSA)

                  DWPPS=DWPPS+DWPPL*RCOS/PCOS
                  DWPTS=DWPTS+DWPTL          /SQRT(PSI)
                  DWTPS=DWTPS+DWTPL          /SQRT(PSI)
                  DWTTS=DWTTS+DWTTL*PCOS/RCOS/PSI

               ENDIF
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FPDWAV

!***********************************************************************
!     calculate local diffusion coefficient
!***********************************************************************
  
      SUBROUTINE FPDWLL(P,PSIN,PCOS,                      &
                        CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ, &
                        DWPPL,DWPTL,DWTPL,DWTTL,NSA)

      USE plprof
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: P,PSIN,PCOS,RKX,RKY,RKZ,RX,RY,RZ
      COMPLEX(rkind),INTENT(IN):: CEX,CEY,CEZ
      REAL(rkind),INTENT(OUT):: DWPPL,DWPTL,DWTPL,DWTTL
      INTEGER,INTENT(IN):: NSA
      INTEGER:: NS,NHMAX,NC,NMI,NPI
      REAL(rkind):: RHON,RW,RWC,RKPARA,RKPERP
      REAL(rkind):: U1X,U1Y,U1Z,U2X,U2Y,U2Z
      REAL(rkind):: RGAMMA,PPARA,PPERP,VPARA,VPERP
      REAL(rkind):: DWC11,DWC12,DWC21,DWC22,RKW,RGZAI
      REAL(rkind):: RJN,RJNM,RJNP,RTHETA2,A11,A12,A21,A22,DWC,EX
      COMPLEX(rkind):: CE1,CE2,CE3,CEPARA,CEPLUS,CEMINUS,CTHETA
      REAL(rkind),DIMENSION(:),POINTER::  RJ,DRJ
      TYPE(pl_mag_type):: mag

      NS=NS_NSA(NSA)
      CALL pl_mag(RX,RY,RZ,MAG)

      RW     =2.D0*PI*RF_WR*1.D6
      RWC    =AEFP(NSA)*MAG%BABS/AMFP(NSA)

      RKPARA = RKX*MAG%BNX +RKY*MAG%BNY +RKZ*MAG%BNZ
      RKPERP = SQRT(RKX*RKX+RKY*RKY+RKZ*RKZ-RKPARA*RKPARA)

      U1X = (RKX-MAG%BNX*RKPARA)/RKPERP
      U1Y = (RKY-MAG%BNY*RKPARA)/RKPERP
      U1Z = (RKZ-MAG%BNZ*RKPARA)/RKPERP

      U2X = (MAG%BNY*RKZ-MAG%BNZ*RKY)/RKPERP
      U2Y = (MAG%BNZ*RKX-MAG%BNX*RKZ)/RKPERP
      U2Z = (MAG%BNX*RKY-MAG%BNY*RKX)/RKPERP

      CE1    = CEX*U1X+CEY*U1Y+CEZ*U1Z
      CE2    = CEX*U2X+CEY*U2Y+CEZ*U2Z
      CEPARA = CEX*MAG%BNX+CEY*MAG%BNY+CEZ*MAG%BNZ

      CEPLUS =(CE1+CI*CE2)/SQRT(2.D0)
      CEMINUS=(CE1-CI*CE2)/SQRT(2.D0)

      RGAMMA =SQRT(1.D0+P*P*THETA0(NS))
      PPARA  =PTFP0(NSA)*P*PCOS
      PPERP  =PTFP0(NSA)*P*PSIN
      VPARA  =PPARA/(AMFP(NSA)*RGAMMA)
      VPERP  =PPERP/(AMFP(NSA)*RGAMMA)

      DWC11=0.D0
      DWC12=0.D0
      DWC21=0.D0
      DWC22=0.D0
      RKW  =RKPARA/RW
      RGZAI=ABS(RKPERP*PPERP/(RWC*AMFP(NSA)))
!      WRITE(6,'(A,1P6E12.4)') 'RGZAI=', &
!                       rx,RKPERP,PPERP,RWC,AMFP(NSA),RGZAI

!      WRITE(26,*) 'RKPERP,PPERP,RWC,RGZAI = ',RKPERP,PPERP,RWC,RGZAI
!      CALL BESSEL(RGZAI,RJ,NCBMAX,NJMAX+1)

      NHMAX=MAX(ABS(NCMIN(NS)-1),ABS(NCMAX(NS)+1),2)
      ALLOCATE(RJ(0:NHMAX),DRJ(0:NHMAX))
      CALL BESSJN(RGZAI,NHMAX,RJ,DRJ)

      DO NC=NCMIN(NS),NCMAX(NS)
         
         IF (NC.LT.0) THEN
             RJN=(-1)**(-NC)*RJ(-NC)
         ELSE
             RJN=RJ(NC)
         ENDIF
         NMI=NC-1
         IF (NMI.LT.0) THEN
             RJNM=(-1)**(-NMI)*RJ(-NMI)
         ELSE
             RJNM=RJ(NMI)
         ENDIF
         NPI=NC+1
         IF (NPI.LT.0) THEN
             RJNP=(-1)**(-NPI)*RJ(-NPI)
         ELSE
             RJNP=RJ(NPI)
         ENDIF
         IF (NC.EQ.0) THEN
             CTHETA=PPERP*(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0) &
                    +CEPARA*PPARA*RJN
         ELSE
             CTHETA=(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0) &
                    +CEPARA*PPARA*(RJNM+RJNP)*RKPERP      &
                    /(2*NC*RWC*AMFP(NSA))
         ENDIF
         RTHETA2=ABS(CTHETA)**2
         IF (NC.EQ.0) THEN
            A11=0
            A12=0
            A21=0
            A22=RTHETA2*RKW**2/(AMFP(NSA)**2*RGAMMA**2)
	 ELSE
            A11=RTHETA2*(1.D0-RKW*VPARA)**2
            A12=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A21=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A22=RTHETA2*RKW**2*VPERP**2
	 ENDIF        
         IF(VPARA.EQ.0.D0) THEN
            DWC=0.D0  
         ELSE
            EX=-((RGAMMA-RKPARA*PPARA/(RW*AMFP(NSA))-NC*RWC/RW) &
                 /(PPARA*DELNPR_WR/(AMFP(NSA)*VC)))**2
            IF (EX.LT.-100.D0) THEN 
                DWC=0.D0
            ELSE
                DWC=0.5D0*SQRT(PI)*AEFP(NSA)**2*EXP(EX)/PTFP0(NSA)**2 &
                     /(RW*ABS(PPARA)*DELNPR_WR/(AMFP(NSA)*VC))
!                WRITE(6,'(A,1P4E12.4)') 'RW,NC*RWC,EX,DWC=',RW,NC*RWC,EX,DWC
            ENDIF
         ENDIF
         DWC11=DWC11+DWC*A11
         DWC12=DWC12+DWC*A12
         DWC21=DWC21+DWC*A21
         DWC22=DWC22+DWC*A22
      END DO
      DEALLOCATE(RJ,DRJ)

      DWPPL=PSIN*(PSIN*DWC11+PCOS*DWC12) &
           +PCOS*(PSIN*DWC21+PCOS*DWC22)
      DWPTL=PSIN*(PCOS*DWC11-PSIN*DWC12) &
           +PCOS*(PCOS*DWC21-PSIN*DWC22)
      DWTPL=PCOS*(PSIN*DWC11+PCOS*DWC12) &
           -PSIN*(PSIN*DWC21+PCOS*DWC22)
      DWTTL=PCOS*(PCOS*DWC11-PSIN*DWC12) &
           -PSIN*(PCOS*DWC21-PSIN*DWC22)
!      IF(ABS(DWPPL).GT.1.D-12) THEN
!         WRITE(6,'(A,1PE12.4)') 'DWPPL=',DWPPL
!         WRITE(6,'(A,1P3E12.4)') 'RL,RZ=',SQRT(RX**2+RY**2),RZ
!      END IF

      RETURN
      END SUBROUTINE FPDWLL

!-------------------------------------


!
!***********************************************************************
!     Calculate PSIN, PCOS, PSI
!***********************************************************************
!
      SUBROUTINE FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI,NSA)
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,NSA
      REAL(rkind),INTENT(IN):: ETAL,RSIN,RCOS
      REAL(rkind),INTENT(OUT):: PSIN,PCOS,PSI
      REAL(rkind):: ARG

      IF(MODELA.EQ.0) THEN
         PSI=1.D0
         PSIN=RSIN
         PCOS=RCOS
      ELSE
         PSI=(1.D0+EPSRM(NR))/(1.D0+EPSRM(NR)*COS(ETAL))
         PSIN=SQRT(PSI)*RSIN
         ARG=1.D0-PSI*RSIN**2
         IF(ARG.LT.0.D0) ARG=0.D0
         IF (RCOS.GT.0.0D0) THEN
            PCOS= SQRT(ARG)
         ELSE
            PCOS=-SQRT(ARG)
         END IF
      ENDIF
      RETURN
      END SUBROUTINE FPDWRP
!---------------------------------
      END MODULE fpcalwr

