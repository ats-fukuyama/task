C     $Id$
C
C***********************************************************************
C     load ray data and calculate spline coefficient
C***********************************************************************
C
      SUBROUTINE FPLDWR(IERR)
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../dp/dpcom1.inc'
      INCLUDE '../wr/wrcom1.inc'
      INCLUDE 'fpcom2.inc'
C
      DIMENSION CFX1(0:NITM,NRAYM),CFX2(0:NITM,NRAYM),
     &          CFX3(0:NITM,NRAYM),FX4(0:NITM,NRAYM),
     &          FX5(0:NITM,NRAYM),FX6(0:NITM,NRAYM),FX7(0:NITM,NRAYM),
     &          FX8(0:NITM,NRAYM),FX9(0:NITM,NRAYM),FX0(0:NITM,NRAYM)
C      
      CALL FROPEN(21,KNAMWR,0,MODEFR,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPWRIN: FROPEN: IERR=',IERR
         RETURN
      ENDIF
C
      REWIND(21)
C
      READ(21) NRAYMX
      DO NRAY=1,NRAYMX
         READ(21) NITMAX(NRAY)
         READ(21) (RAYIN(I,NRAY),I=1,8)
         READ(21) (CEXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (CEYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (CEZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (RKXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (RKYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (RKZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (RXS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
         READ(21) (RYS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
         READ(21) (RZS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
         READ(21) (RAYRB1(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(21) (RAYRB2(NIT,NRAY),NIT=0,NITMAX(NRAY))
         DO I=0,8
            READ(21) (RAYS(I,NIT,NRAY),NIT=0,NITMAX(NRAY))
         ENDDO
      ENDDO
      CLOSE(21)
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      WRITE(6,*) 'NRAYMX=',NRAYMX
C
C     ----- Calculate spline coefficients -----
C
      DO NRAY=1,NRAYMX
         NITMX=NITMAX(NRAY)
         DO NIT=0,NITMX
            CALL PLMAG(RXS(NIT,NRAY),RYS(NIT,NRAY),RZS(NIT,NRAY),
     &                 RHON)
            PSIX(NIT,NRAY)=RHON**2
            SI(NIT,NRAY)=RAYS(0,NIT,NRAY)
         ENDDO
C
         CALL CSPL1D(SI(0,NRAY),CEXS(0,NRAY),CFX1(0,NRAY),
     &               CU1(1,0,NRAY),NITMX+1,0,IERR)
         CALL CSPL1D(SI(0,NRAY),CEYS(0,NRAY),CFX2(0,NRAY),
     &               CU2(1,0,NRAY),NITMX+1,0,IERR)
         CALL CSPL1D(SI(0,NRAY),CEZS(0,NRAY),CFX3(0,NRAY),
     &               CU3(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RKXS(0,NRAY),FX4(0,NRAY),
     &               U4(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RKYS(0,NRAY),FX5(0,NRAY),
     &               U5(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RKZS(0,NRAY),FX6(0,NRAY),
     &               U6(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RXS(0,NRAY),FX7(0,NRAY),
     &               U7(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RYS(0,NRAY),FX8(0,NRAY),
     &               U8(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RZS(0,NRAY),FX9(0,NRAY),
     &               U9(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),PSIX(0,NRAY),FX0(0,NRAY),
     &               U0(1,0,NRAY),NITMX+1,0,IERR)
      ENDDO
C
C     ----- Find crossing point -----
C
      DO NRAY=1,NRAYMX
         NITMX=NITMAX(NRAY)
C
      DO NR=1,NRMAX
         PSICR =RM(NR)**2
         PSIPRE=PSIX(0,NRAY)
         NCR=0
         DO NIT=1,NITMX
            PSIL=PSIX(NIT,NRAY)
            IF((PSIPRE-PSICR)*(PSIL-PSICR).LT.0.D0.OR.
     &          PSIL-PSICR.EQ.0.D0) THEN
               CALL FPCROS(PSICR,NIT,NRAY,SICR)
               WRITE(6,601) '# NR,NRAY,NIT,PSICR,SICR=',
     &                         NR,NRAY,NIT,PSICR,SICR
  601          FORMAT(A,3I4,1P2E15.6)
               CALL FPCREK(SICR,NRAY,
     &                     CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ)
               IF(NCR.LT.NCRM) THEN
                  NCR=NCR+1
                  CECR(1,NCR,NR,NRAY)=CEX
                  CECR(2,NCR,NR,NRAY)=CEY
                  CECR(3,NCR,NR,NRAY)=CEZ
                  RKCR(1,NCR,NR,NRAY)=RKX
                  RKCR(2,NCR,NR,NRAY)=RKY
                  RKCR(3,NCR,NR,NRAY)=RKZ
                  RCR(1,NCR,NR,NRAY)=RX
                  RCR(2,NCR,NR,NRAY)=RY
                  RCR(3,NCR,NR,NRAY)=RZ
               ENDIF
            ENDIF
            PSIPRE=PSIL
         ENDDO
         NCRMAX(NR,NRAY)=NCR
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         RHOL =RM(NR)
         CALL PLRRMX(RHOL,RRMIN(NR),RRMAX(NR))
C         WRITE(6,602) NR,RHOL,RRMIN(NR),RRMAX(NR)
C  602    FORMAT('# NR,RHOL,RRMIN,RRMAX=',I3,1P3E15.7)
      ENDDO
C
  900 IERR=0
      RETURN
      END
C
C***********************************************************************
C     Calculate crossing point by Newton's method
C***********************************************************************
C
      SUBROUTINE FPCROS(PSICR,NIT,NRAY,SICR)
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
      INCLUDE 'fpcom2.inc'
C
      Y   =PSIX(NIT-1,NRAY)-PSICR
      YNEW=PSIX(NIT  ,NRAY)-PSICR
      X   =SI(NIT-1,NRAY)
      DX  =SI(NIT,NRAY)-SI(NIT-1,NRAY)
      IF(DX.EQ.0.D0) GOTO 8000
      DYDX=(YNEW-Y)/DX
      Y   =YNEW
C
      ICOUNT=0
C
  100 CONTINUE
         ICOUNT=ICOUNT+1
C         WRITE(6,'(A,I3,1P4E15.7)') 'FPCROS:',ICOUNT,X,DX,Y,YNEW
         IF(ABS(Y).LE.EPSNWR) GOTO 7000
         IF(ICOUNT.GT.LMAXNWR) GOTO 8100
         IF(DYDX.EQ.0.D0) GOTO 8200
         DX=-Y/DYDX
         IF(DX.EQ.0.D0) GOTO 8300
         X=X+DX
         CALL SPL1DF(X,PSIL,
     &               SI(0,NRAY),U0(1,0,NRAY),NITMAX(NRAY),IERR)
         YNEW=PSIL-PSICR
         DYDX=(YNEW-Y)/DX
         Y=YNEW
      GOTO 100
C
 7000 SICR=X
      RETURN
C
 8000 WRITE(6,*) 'XX FPCROS: INVALID INITIAL VALUES: X0=X1'
      GOTO 7000
 8100 WRITE(6,*) 'XX FPCROS: ICOUNT EXCEEDS LMAXNW AT SI =',X
      GOTO 7000
 8200 WRITE(6,*) 'XX FPCROS: DYDX BECOMES 0 AT SI =',X
      GOTO 7000
 8300 WRITE(6,*) 'XX FPCROS: DX BECOMES 0 AT SI = ',X
      GOTO 7000
      END
C
C***********************************************************************
C     Calculate wave field, wave number and position at crossing point
C***********************************************************************
C
      SUBROUTINE FPCREK(SICR,NRAY,
     &                  CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ)
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
      INCLUDE 'fpcom2.inc'
C
      NITMX=NITMAX(NRAY)
      CALL CSPL1DF(SICR,CEX,SI(0,NRAY),CU1(1,0,NRAY),NITMX,IERR)
      CALL CSPL1DF(SICR,CEY,SI(0,NRAY),CU2(1,0,NRAY),NITMX,IERR)
      CALL CSPL1DF(SICR,CEZ,SI(0,NRAY),CU3(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RKX,SI(0,NRAY),U4(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RKY,SI(0,NRAY),U5(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RKZ,SI(0,NRAY),U6(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RX,SI(0,NRAY),U7(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RY,SI(0,NRAY),U8(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RZ,SI(0,NRAY),U9(1,0,NRAY),NITMX,IERR)
C
      RETURN
      END
