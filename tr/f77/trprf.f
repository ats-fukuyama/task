C     $Id$
C
C     ***********************************************************
C
C           RF HEATING AND CURRENT DRIVE (GAUSSIAN PROFILE)
C
C     ***********************************************************
C
      SUBROUTINE TRPWRF
C
      INCLUDE 'trcomm.inc'
C
      IF(PECTOT+PLHTOT+PICTOT.LE.0.D0) RETURN
C
      IF(PLHR0.LT.0.D0) THEN
         VPHLH=VC/(PLHNPR*ABS(PLHR0))
         DO NR=NRMAX,1,-1
            VTE=SQRT(ABS(RT(NR,1))*RKEV/AME)
            IF(VTE.GT.VPHLH) THEN
               IF(NR.EQ.NRMAX) THEN
                  PLHR0L=RA
               ELSE
                  VTEP=SQRT(ABS(RT(NR+1,1))*RKEV/AME)
                  FACT=(VTE-VPHLH)/(VTE-VTEP)
                  PLHR0L=(FACT*RM(NR+1)+(1.D0-FACT)*RM(NR))*RA
               ENDIF
               GOTO 6
            ENDIF
         ENDDO
         PLHR0L=0.D0
    6    CONTINUE
C         WRITE(6,*) '*** PLHR0L = ',PLHR0L
      ELSE
         PLHR0L=PLHR0
         VPHLH=VC/PLHNPR
      ENDIF
C
      SUMEC = 0.D0
      SUMLH = 0.D0
      SUMIC = 0.D0
      DO NR=1,NRMAX
         SUMEC = SUMEC
     &         +DEXP(-((RA*RM(NR)-PECR0 )/PECRW)**2)*DVRHO(NR)*DR
         SUMLH = SUMLH
     &         +DEXP(-((RA*RM(NR)-PLHR0L)/PLHRW)**2)*DVRHO(NR)*DR
         SUMIC = SUMIC
     &         +DEXP(-((RA*RM(NR)-PICR0 )/PICRW)**2)*DVRHO(NR)*DR
      ENDDO
C
      PEC0 = PECTOT*1.D6/SUMEC
      PLH0 = PLHTOT*1.D6/SUMLH
      PIC0 = PICTOT*1.D6/SUMIC
C
C      IF(ABS(PLHNPR).LE.1.D0) THEN
C         NLH=PLHR0/DR+1.D0
C         RTLH=((RM(NLH+1)-PLHR0)*RT(NLH,1)
C     &        +(PLHR0-RM(NLH))*RT(NLH+1,1))/DR
C         VTLH=SQRT(RTLH*RKEV/AME)
C         VPHLH=VTLH/ABS(PLHNPR)
C         IF(VPHLH.GT.VC) VPHLH=VC
C         IF(PLHNPR.LT.0.D0) VPHLH=-VPHLH
C         WRITE(6,*) '** PLHNPR = ',VC/VPHLH
C      ELSE
C         VPHLH=VC/PLHNPR
C      ENDIF
C
      DO NR=1,NRMAX
         PECL = PEC0*DEXP(-((RA*RM(NR)-PECR0)/PECRW)**2)
         PLHL = PLH0*DEXP(-((RA*RM(NR)-PLHR0L)/PLHRW)**2)
         PICL = PIC0*DEXP(-((RA*RM(NR)-PICR0)/PICRW)**2)
         PRF (NR,1  )=PECTOE*PECL
     &               +PLHTOE*PLHL
     &               +PICTOE*PICL
         PRFV(NR,1,1)=PECTOE*PECL
         PRFV(NR,1,2)=PLHTOE*PLHL
         PRFV(NR,1,3)=PICTOE*PICL
         PRF (NR,2  )=(1.D0-PECTOE)*PECL
     &               +(1.D0-PLHTOE)*PLHL
     &               +(1.D0-PICTOE)*PICL
         PRFV(NR,2,1)=(1.D0-PECTOE)*PECL
         PRFV(NR,2,2)=(1.D0-PLHTOE)*PLHL
         PRFV(NR,2,3)=(1.D0-PICTOE)*PICL
C
         RLNLMD=16.1D0 - 1.15D0*LOG10(RN(NR,1))
     &                 + 2.30D0*LOG10(RT(NR,1))
         VTE=SQRT(ABS(RT(NR,1))*RKEV/AME)
C
         IF(PECCD.NE.0.D0) THEN
            VPHEC=VC/(VTE*PECNPR)
            IF(PECNPR.LE.1.D0) THEN
               EFCDEC=0.D0
            ELSE
               EFCDEC=TRCDEF(VPHEC,ZEFF(NR),0.D0,EPSRHO(NR),0)
            ENDIF
         ELSE
            EFCDEC=0.D0
         ENDIF
         IF(PLHCD.NE.0.D0) THEN
            IF(ABS(VPHLH).GT.VC) THEN
               EFCDLH=0.D0
            ELSE
               EFCDLH=TRCDEF(VPHLH/VTE,ZEFF(NR),0.D0,EPSRHO(NR),0)
            ENDIF
         ELSE
            EFCDLH=0.D0
         ENDIF
         IF(PICCD.NE.0.D0) THEN
            VPHIC=VC/(VTE*PICNPR)
            IF(PICNPR.LE.1.D0) THEN
               EFCDIC=0.D0
            ELSE
               EFCDIC=TRCDEF(VPHIC,ZEFF(NR),0.D0,EPSRHO(NR),1)
            ENDIF
         ELSE
            EFCDIC=0.D0
         ENDIF
         AJRFV(NR,1)=0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD)
     &              *(PECCD*PECTOE*EFCDEC*PECL)
         AJRFV(NR,2)=0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD)
     &              *(PLHCD*PLHTOE*EFCDLH*PLHL)
         AJRFV(NR,3)=0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD)
     &              *(PICCD*PICTOE*EFCDIC*PICL)
         AJRF(NR)=0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD)
     &           *(PECCD*PECTOE*EFCDEC*PECL
     &            +PLHCD*PLHTOE*EFCDLH*PLHL
     &            +PICCD*PICTOE*EFCDIC*PICL)
      ENDDO
C
      RETURN
      END
C
C     ****** CURRENT DRIVE EFFICIENCY ******
C
C      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
C      Z  = ZEFF   : EFFECTIVE Z
C      XR = X / RR : NORMALIZED X
C      YR = Y / RR : NORMALIZED Y
C      ID : 0 : LANDAU DAMPING
C           1 : TTMP
C
      REAL*8 FUNCTION TRCDEF(WT,Z,XR,YR,ID)
C
      IMPLICIT NONE
      INTEGER ID
      REAL*8 WT,Z,XR,YR
      REAL*8 R,A,C,D,W,RM,RC,EFF0,EFF1,EFF2,EFF3,Y1,Y2,YT,ARG
C
      R=SQRT(XR*XR+YR*YR)
      IF(ID.EQ.0) THEN
         D=3.D0/Z
         C=3.83D0
         A=0.D0
         RM=1.38D0
         RC=0.389D0
      ELSE
         D=11.91D0/(0.678D0+Z)
         C=4.13D0
         A=12.3D0
         RM=2.48D0
         RC=0.0987D0
      ENDIF
      IF(WT.LE.1.D-20) THEN
         W=1.D-20
      ELSE
         W=WT
      ENDIF
      EFF0=D/W+C/Z**0.707D0+4.D0*W*W/(5.D0+Z)
      EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)
C
      Y2=(R+XR)/(1.D0+R)
      IF(Y2.LT.0.D0) Y2=0.D0
      Y1=SQRT(Y2)
      EFF2=1.D0+A*(Y1/W)**3
C
      IF(Y2.LE.1.D-20) THEN
         YT=(1.D0-Y2)*WT*WT/1.D-60
      ELSE
         YT=(1.D0-Y2)*WT*WT/Y2
      ENDIF
      IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
         ARG=(RC*YT)**RM
         IF(ARG.LE.100.D0) THEN
            EFF3=1.D0-MIN(EXP(-ARG),1.D0)
         ELSE
            EFF3=1.D0
         ENDIF
      ELSE
         EFF3=1.D0
      ENDIF
C
      TRCDEF=EFF0*EFF1*EFF2*EFF3
      RETURN
      END
