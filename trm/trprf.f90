!     ***********************************************************

!           RF HEATING AND CURRENT DRIVE (GAUSSIAN PROFILE)

!     ***********************************************************

      SUBROUTINE TRPWRF

        USE TRCOMM, ONLY : &
             AJRF, AJRFV, AME, DR, DVRHO, EPSRHO, NRMAX, PECCD, PECNPR, &
             PECR0, PECRW, PECTOE, PECIN, PICCD, PICNPR, PICR0, PICRW, &
             PICTOE, PICIN, PLHCD, PLHNPR, PLHR0, PLHRW, PLHTOE, PLHIN, &
             PRF, PRFV, RA, RKEV, RM, RN, RT, VC, ZEFF, rkind, &
             PEC_NEC,PLH_NLH,PIC_NIC,PECTOT,PLHTOT,PICTOT, &
             NECMAX,NLHMAX,NICMAX,AJRFV,PRFV
      IMPLICIT NONE
      REAL(rkind)   :: &
           EFCDEC, EFCDIC, EFCDLH, FACT, PEC0, PECL, PIC0, PICL, PLH0, &
           PLHL, PLHR0L, RLNLMD, SUMEC, SUMIC, SUMLH, VPHEC, VPHIC, VPHLH, &
           VTE, VTEP
      INTEGER:: NR,NEC,NLH,NIC
      REAL(rkind)   :: TRCDEF

      PECTOT=0.D0
      PLHTOT=0.D0
      PICTOT=0.D0
      DO NEC=1,NECMAX
         PECTOT=PECTOT+PECIN(NEC)
      END DO
      DO NLH=1,NLHMAX
         PLHTOT=PLHTOT+PLHIN(NEC)
      END DO
      DO NIC=1,NICMAX
         PICTOT=PICTOT+PICIN(NEC)
      END DO
      IF(PECTOT+PLHTOT+PICTOT.LE.0.D0) THEN
         PRF(1:NRMAX,1)=0.D0
         PRF(1:NRMAX,2)=0.D0
         PRF(1:NRMAX,3)=0.D0
         AJRF(1:NRMAX)=0.D0
         DO NEC=1,NECMAX
            PEC_NEC(NEC,1:NRMAX)=0.D0
         END DO
         DO NLH=1,NLHMAX
            PLH_NLH(NLH,1:NRMAX)=0.D0
         END DO
         DO NIC=1,NICMAX
            PIC_NIC(NIC,1:NRMAX)=0.D0
         END DO
         RETURN
      END IF

      DO NR=1,NRMAX
         PRFV(NR,1,1)=0.D0
         PRFV(NR,2,1)=0.D0
         AJRFV(NR,1)=0.D0
      END DO
      
      DO NEC=1,NECMAX
         SUMEC = 0.D0
         DO NR=1,NRMAX
            SUMEC = SUMEC + DEXP(-((RA*RM(NR)-PECR0(NEC) ) &
                                  /PECRW(NEC))**2)*DVRHO(NR)*DR
         END DO
         PEC0 = PECIN(NEC)*1.D6/SUMEC
         DO NR=1,NRMAX
            PECL = PEC0*DEXP(-((RA*RM(NR)-PECR0(NEC))/PECRW(NEC))**2)
            PRFV(NR,1,1)=PRFV(NR,1,1)+      PECTOE(NEC) *PECL
            PRFV(NR,2,1)=PRFV(NR,2,1)+(1.D0-PECTOE(NEC))*PECL

            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(NR,1))  + 2.30D0*LOG10(RT(NR,1))
            VTE=SQRT(ABS(RT(NR,1))*RKEV/AME)

            IF(PECCD(NEC).NE.0.D0) THEN
               VPHEC=VC/(VTE*PECNPR(NEC))
               IF(PECNPR(NEC).LE.1.D0) THEN
                  EFCDEC=0.D0
               ELSE
                  EFCDEC=TRCDEF(VPHEC,ZEFF(NR),0.D0,EPSRHO(NR),0)
               ENDIF
            ELSE
               EFCDEC=0.D0
            ENDIF
            AJRFV(NR,1)=AJRFV(NR,1)+0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD) &
                 *(PECCD(NEC)*PECTOE(NEC)*EFCDEC*PECL)
         END DO
      END DO

      DO NLH=1,NLHMAX
         IF(PLHR0(NLH).LT.0.D0) THEN
            VPHLH=VC/(PLHNPR(NLH)*ABS(PLHR0(NLH)))
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
6           CONTINUE
            !         WRITE(6,*) '*** PLHR0L = ',PLHR0L
         ELSE
            PLHR0L=PLHR0(NLH)
            VPHLH=VC/PLHNPR(NLH)
         ENDIF

         SUMLH = 0.D0
         DO NR=1,NRMAX
            SUMLH = SUMLH + DEXP(-((RA*RM(NR)-PLHR0L) &
                                   /PLHRW(NLH))**2)*DVRHO(NR)*DR
         ENDDO
         PLH0 = PLHIN(NLH)*1.D6/SUMLH

         DO NR=1,NRMAX
            PLHL = PLH0*DEXP(-((RA*RM(NR)-PLHR0L)/PLHRW(NLH))**2)
            PRFV(NR,1,2)=PRFV(NR,1,2)+PLHTOE(NLH)*PLHL
            PRFV(NR,2,2)=PRFV(NR,2,2)+(1.D0-PLHTOE(NLH))*PLHL

            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(NR,1))  + 2.30D0*LOG10(RT(NR,1))
            VTE=SQRT(ABS(RT(NR,1))*RKEV/AME)

            IF(PLHCD(NLH).NE.0.D0) THEN
               IF(ABS(VPHLH).GT.VC) THEN
                  EFCDLH=0.D0
               ELSE
                  EFCDLH=TRCDEF(VPHLH/VTE,ZEFF(NR),0.D0,EPSRHO(NR),0)
               ENDIF
            ELSE
               EFCDLH=0.D0
            ENDIF
            AJRFV(NR,2)=AJRFV(NR,2)+0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD) &
                 *(PLHCD(NLH)*PLHTOE(NLH)*EFCDLH*PLHL)
         END DO
      END DO
      
      DO NIC=1,NICMAX
         SUMIC = 0.D0
         DO NR=1,NRMAX
            SUMIC = SUMIC + DEXP(-((RA*RM(NR)-PICR0(NIC) ) &
                                  /PICRW(NIC))**2)*DVRHO(NR)*DR
         ENDDO
         PIC0 = PICIN(NIC)*1.D6/SUMIC

         DO NR=1,NRMAX
            PICL = PIC0*DEXP(-((RA*RM(NR)-PICR0(NIC))/PICRW(NIC))**2)
            PRFV(NR,1,3)=PRFV(NR,1,3)+PICTOE(NIC)*PICL
            PRFV(NR,2,3)=PRFV(NR,2,3)+(1.D0-PICTOE(NIC))*PICL

            RLNLMD=16.1D0 - 1.15D0*LOG10(RN(NR,1))  + 2.30D0*LOG10(RT(NR,1))
            VTE=SQRT(ABS(RT(NR,1))*RKEV/AME)

            IF(PICCD(NIC).NE.0.D0) THEN
               VPHIC=VC/(VTE*PICNPR(NIC))
               IF(PICNPR(NIC).LE.1.D0) THEN
                  EFCDIC=0.D0
               ELSE
                  EFCDIC=TRCDEF(VPHIC,ZEFF(NR),0.D0,EPSRHO(NR),1)
               ENDIF
            ELSE
               EFCDIC=0.D0
            ENDIF
            AJRFV(NR,3)=AJRFV(NR,3)+0.384D0*RT(NR,1)/(RN(NR,1)*RLNLMD) &
                       *(PICCD(NIC)*PICTOE(NIC)*EFCDIC*PICL)
         END DO
      END DO

      DO NR=1,NRMAX
         PRF(NR,1)=PRFV(NR,1,1)+PRFV(NR,1,2)+PRFV(NR,1,3)
         PRF(NR,2)=PRFV(NR,2,1)+PRFV(NR,2,2)+PRFV(NR,2,3)
         AJRF(NR)=AJRFV(NR,1)+AJRFV(NR,2)+AJRFV(NR,3)
      END DO

      RETURN
      END SUBROUTINE TRPWRF

!     ****** CURRENT DRIVE EFFICIENCY ******

!      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
!      Z  = ZEFF   : EFFECTIVE Z
!      XR = X / RR : NORMALIZED X
!      YR = Y / RR : NORMALIZED Y
!      ID : 0 : LANDAU DAMPING
!           1 : TTMP

      FUNCTION TRCDEF(WT,Z,XR,YR,ID)

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      INTEGER ID
      REAL(rkind) WT,Z,XR,YR, TRCDEF
      REAL(rkind) R,A,C,D,W,RM,RC,EFF0,EFF1,EFF2,EFF3,Y1,Y2,YT,ARG

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

      Y2=(R+XR)/(1.D0+R)
      IF(Y2.LT.0.D0) Y2=0.D0
      Y1=SQRT(Y2)
      EFF2=1.D0+A*(Y1/W)**3

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

      TRCDEF=EFF0*EFF1*EFF2*EFF3
      RETURN
      END FUNCTION TRCDEF
