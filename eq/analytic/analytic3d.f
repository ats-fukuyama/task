C
C     *** NUMERICAL vs. ANALYTIC SOLUTIONS ***
C
C     HISTORY
C
C     2002/01/16
C     2002/03/01
C     2002/03/05
C     2002/03/06
C     2002/03/09
C     2002/03/13
C     2002/03/19
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
      PARAMETER (NXM=64**2,NZM=64**2)
      DIMENSION GX(NXM),GY(NZM),GZ(NXM,NZM),GJY(NXM)
      DIMENSION RX(NXM),ZX(NZM),PSI(NXM,NZM),HJ(NXM)
      DIMENSION PSIZAX(NXM)
      CHARACTER TXT*80
C
      CALL GSOPEN
C
C     MUL : Only a multiple of 64
C
      MUL   = 64
      NXMAX = 64*MUL
      NZMAX = 64*MUL
C
C     *********************************************
C     *              GIVEN VALUES                 *
C     *********************************************
C
C     *** CONSTANTS ***
C
C        PI    : Pi
C        RMU0  : Permeability of free space
C        AMP   : Proton mass
C        AEE   : Electron charge
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AMP    = 1.6726231D-27
      AEE    = 1.60217733D-19
C
C     *** CONFIGURATION PARAMETERS ***
C
C        RR    : Plasma major radius                             (m)
C        RA    : Plasma minor radius                             (m)
C        ZE    : Plasma range of z axis                          (m)
C        EPSA  : Cross-section parameter ; -1 < EPSA < 1
C        VPHI  : Troidal rotation velocity                     (m/s)
C        TMP   : Plasma temperature                            (keV)
C        PINT  : Plasma pressure                               (MPa)
C
      RR     = 3.D0
      RA     = 1.D0
      ZE     = RA
      EPSA   = 0.D0
      VPHI   = 5.D4
      TMP    = 6.D0*1.D3
      PINT   = 1.D-1*1.D6
      RC     = RR-RA
C
C     *** CONSTANTS OBTAINED BY ANALYTIC CALCULATION ***
C
      RAXIS  = 3.16275D0
      ZAXIS  = 0.D0
C
C     *** CHANGE PARAMETER ***
C
      NCHG   = 1
      NLRC   = 0
C
C     *** ARBITRARY CONSTANT ***
C
C        HJ0   : Current density at R=RC                   (MA/m**2)
C
      HJ0    = 1.335477D-1*1.D6
C
C     *********************************************
C     *            CALCULATE CONSTANT             *
C     *********************************************
C
      OTC = (VPHI**2*AMP)/(RR**2*(TMP*AEE))
C
C     *********************************************
C     *            ANALYTIC PSI,J                 *
C     *********************************************
C
      CC1=((EPSA-1.D0)*RAXIS**2/(8.D0*RR**2))+(1.D0/(2.D0*RR**2*OTC))
     &     *(EXP(OTC*RAXIS**2/2.D0)-1.D0)
      A1=(CC1/RR**2)*(RAXIS**2-(RR-RA)**2)
      A1R=(CC1/RR**2)*(RAXIS**2-(RR+RA)**2)
      A2=-((EPSA-1.D0)*RAXIS**4)/(1.6D1*RR**4)+(1/(RR**4*OTC**2))
     &     *(1.D0+RAXIS**2*OTC/2.D0-EXP(RAXIS**2*OTC/2.D0))
      A3=((EPSA-1.D0)*(RR-RA)**4)/(1.6D1*RR**4)-(1/(RR**4*OTC**2))
     &     *(1.D0+(RR-RA)**2*OTC/2.D0-EXP((RR-RA)**2*OTC/2.D0))
      A3R=((EPSA-1.D0)*(RR+RA)**4)/(1.6D1*RR**4)-(1/(RR**4*OTC**2))
     &     *(1.D0+(RR+RA)**2*OTC/2.D0-EXP((RR+RA)**2*OTC/2.D0))
      IF (A1+A2+A3.LT.0.D0) THEN
         AA=-(A1+A2+A3)
      ELSEIF (A1+A2+A3.EQ.0.D0) THEN
         WRITE(6,*) "ERROR-LEFT"
         STOP
      ELSEIF (A1+A2+A3.GT.0.D0) THEN
         AA=A1+A2+A3
      ENDIF
C
      IF (A1R+A2+A3R.LT.0.D0) THEN
         AAR=-(A1R+A2+A3R)
      ELSEIF (A1R+A2+A3R.EQ.0.D0) THEN
         WRITE(6,*) "ERROR-RIGHT"
         STOP
      ELSEIF (A1R+A2+A3R.GT.0.D0) THEN
         AAR=A1R+A2+A3R
      ENDIF
      PP   = SQRT(RMU0*RR**4*PINT/AA)
      PPR  = SQRT(RMU0*RR**4*PINT/AAR)
      HMT = 2.D0*RMU0*PINT*RC**2
     &     *((RR**4/PP*RC)*RMU0*HJ0-EXP(RC**2*OTC/2.D0))
      HMTR = 2.D0*RMU0*PINT*RC**2
     &     *((RR**4/PPR*RC)*RMU0*HJ0-EXP(RC**2*OTC/2.D0))
      HM   = (HMT*PP)/(2.D0*RMU0*PINT*RR**2)
      HMR  = (HMTR*PPR)/(2.D0*RMU0*PINT*RR**2)
      CC2  = -CC1*PP*(RR-RA)**2/RR**2+PP*A3
      CC2R = -CC1*PPR*(RR+RA)**2/RR**2+PPR*A3R
C
      IF (HMT.LE.-0.4D0) THEN
         DD = 4.D0
      ELSEIF (HMT.GE.0.1D0) THEN
         DD = 2.D0
      ELSE
         DD = 3.D0
      ENDIF
C
C     *** CLACULATE M limit ***
C
      HMLMT=-(PP*(RR-RA)**2/RR**2)*EXP(OTC*(RR-RA)**2/2.D0)
C      WRITE(6,*) "HM >",HMLMT
C
C     *** CALCULATE RKAP around magnetic axis***
C
C        RKAP  : Plasma shape elongation
C
      HNUM = EXP(OTC*RAXIS**2/2.D0)-(1.D0-EPSA)/2.D0
      IF (NLRC.EQ.0) THEN
         DEN = ((RR**2*HM)/(RAXIS**2*PP))+(1.D0-EPSA)/2.D0
      ELSE
         DEN = ((RR**2*HMR)/(RAXIS**2*PPR))+(1.D0-EPSA)/2.D0
      ENDIF
      RKAPRA = SQRT(HNUM/DEN)
C
C     ***
C
C      DX = (2.D0*RA)/(NXMAX-1)
C      DZ = (2.D0*ZE)/(NZMAX-1)
      DX = DD/(NXMAX-1)
      DZ = DD/(NZMAX-1)
C
      ZDLTMX=0.D0
      RRMIN=3.D0
      RRMAX=3.D0
      ZZMIN=0.D0
      ZZMAX=0.D0
      PSIM=1.D0
      PSIMAX=0.D0
      DO NX=1,NXMAX
      DO NZ=1,NZMAX
C
C         IF (NX.LE.NXMAX/2) THEN
C            RX(NX)=(RR-RA)+DX*DBLE(NX-1)
C         ELSE
C            RX(NX)=(RR+DBLE(DX/2))+DX*DBLE(NX-NXMAX/2-1)
C         ENDIF
C         IF (NX.LE.NZMAX/2) THEN
C            ZX(NZ)=-DZ*(DBLE(NZMAX/2-1)+0.5D0)+DZ*DBLE(NZ-1)
C         ELSE
C            ZX(NZ)=DBLE(DZ/2)+DZ*DBLE(NZ-NZMAX/2-1)
C         ENDIF
         RX(NX) = (RR-(DD/2.D0))+DX*DBLE(NX-1)
         ZX(NZ) = -(DD/2.D0)+DZ*DBLE(NZ-1)
C
         IF (NLRC.EQ.0) THEN
         PSI(NX,NZ)=CC1*PP*RX(NX)**2/RR**2
     &        -(HM*ZX(NZ)**2/(2.D0*RR**2))
     &        +PP*((EPSA-1.D0)/4.D0)*(ZX(NZ)**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PP/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2
         ELSE
         PSI(NX,NZ)=CC1*PPR*RX(NX)**2/RR**2
     &        -(HMR*ZX(NZ)**2/(2.D0*RR**2))
     &        +PPR*((EPSA-1.D0)/4.D0)*(ZX(NZ)**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PPR/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2R
         ENDIF
         IF (PSI(NX,NZ).GT.0.D0) THEN
            IF (PSI(NX,NZ).GT.PSIMAX) PSIMAX=PSI(NX,NZ)
            IF (ZX(NZ).GT.ZDLTMX) THEN
               ZDLTMX=ZX(NZ)
               RDLTMX=RX(NX)
            ENDIF
            IF (RX(NX).LT.RRMIN) RRMIN=RX(NX)
            IF (RX(NX).GT.RRMAX) RRMAX=RX(NX)
            IF (ZX(NZ).LT.ZZMIN) ZZMIN=ZX(NZ)
            IF (ZX(NZ).GT.ZZMAX) ZZMAX=ZX(NZ)
         ENDIF
         IF (ABS(PSI(NX,NZ)).LT.PSIM) THEN
            PSIM=ABS(PSI(NX,NZ))
            PSIMIN=PSI(NX,NZ)
         ENDIF
      ENDDO
         IF (NLRC.EQ.0) THEN
         PSIZAX(NX)=CC1*PP*RX(NX)**2/RR**2
     &        -(HM*ZAXIS**2/(2.D0*RR**2))
     &        +PP*((EPSA-1.D0)/4.D0)*(ZAXIS**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PP/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2
         HJ(NX)=(HM/(RMU0*RR**2*RX(NX))+(RX(NX)*PP/(RMU0*RR**4))
     &           *EXP(RX(NX)**2*OTC/2.D0))*1.D-6
         ELSE
         PSIZAX(NX)=CC1*PPR*RX(NX)**2/RR**2
     &        -(HMR*ZAXIS**2/(2.D0*RR**2))
     &        +PPR*((EPSA-1.D0)/4.D0)*(ZAXIS**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PPR/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2R
         HJ(NX)=(HMR/(RMU0*RR**2*RX(NX))+(RX(NX)*PPR/(RMU0*RR**4))
     &           *EXP(RX(NX)**2*OTC/2.D0))*1.D-6
         ENDIF
      ENDDO
C      write(6,*) rrmin,rrmax,zzmin,zzmax,psimin,psimax
C
C     *********************************************
C     *          Calculate RKAP,RDLT              *
C     *********************************************
C
      A=ABS(RRMAX-RRMIN)
      B=ABS(ZZMAX-ZZMIN)
      RKAP=B/A
      C=0.5D0*(RRMAX-RRMIN)
      D=ABS((RRMAX-C)-RDLTMX)
      RDLT=D/C
C
C     *********************************************
C     *           DISPLAY PARAMETERS              *
C     *********************************************
C
 200  FORMAT(3H  |,8X,'PP =',1X,1F18.15,3H  |)
 210  FORMAT(3H  |,8X,'HM =',1X,1F18.15,3H  |)
 220  FORMAT(3H  |,' RKAP_AXIS =',1X,1F18.15,3H  |)
 230  FORMAT(3H  |,6X,'RKAP =',1X,1F18.15,3H  |)
 240  FORMAT(3H  |,6X,'DELT =',1X,1F18.15,3H  |)
 250  FORMAT(2H  ,' ------ CALCULATION RESULTS ------ ')
 260  FORMAT(2H  ,' --------------------------------- ')
      WRITE(6,250)
      WRITE(6,200) PP
      WRITE(6,210) HM
      WRITE(6,220) RKAPRA
      WRITE(6,230) RKAP
      WRITE(6,240) RDLT
      WRITE(6,260)
C
C     *********************************************
C     *               DRAW GRAPH                  *
C     *********************************************
C
      DO NX=1,NXMAX
         IF (MOD(NX,MUL).EQ.1) THEN
            DO NZ=1,NZMAX
               IF (MOD(NZ,MUL).EQ.1) THEN
                  GX((NX-1)/MUL+1)   = GUCLIP(RX(NX))
                  GY((NZ-1)/MUL+1)   = GUCLIP(ZX(NZ))
                  GZ((NX-1)/MUL+1,(NZ-1)/MUL+1) = GUCLIP(PSI(NX,NZ))
                  GJY((NX-1)/MUL+1)  = GUCLIP(HJ(NX))
               ENDIF
C               GX(NX)   = GUCLIP(RX(NX))
C               GY(NZ)   = GUCLIP(ZX(NZ))
C               GZ(NX,NZ) = GUCLIP(PSI(NX,NZ))
C               GJY(NX)  = GUCLIP(HJ(NX))
            ENDDO
         ENDIF
      ENDDO
C
C      CALL PAGES
C      TXT='/ANALYTIC AND NUMERICAL PSI/'
C      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GY,NXM,NXMAX,2,TXT,5.5)
C      CALL PAGEE
C
      CALL GRAPH3(15.0,15.0,8.0,GX,GY,GZ,NXM,NXMAX/MUL)
C
      CALL PAGES
      TXT='/ANALYTIC J/'
      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GJY,NXM,NXMAX/MUL,1,TXT,3.5)
      CALL PAGEE
      GOTO 9000
C
 9000 CALL GSCLOS
      STOP
      END
C
C     *********************************************
C     *               SUBROUTINES                 *
C     *********************************************
C
      SUBROUTINE GRAPH1(PXMIN,PXMAX,PYMIN,PYMAX,GX,GY,NXM,
     &                   NXMAX,NGMAX,TXT,POS)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX)
      CHARACTER*80 TXT
C
      CALL GMNMX1(GX,1,NXMAX,1,XMIN,XMAX)
      CALL GMNMX1(GY(1,1),1,NXMAX,1,YMIN,YMAX)
      DO 100 NG=2,NGMAX
         CALL GMNMX1(GY(1,NG),1,NXMAX,1,YMIN1,YMAX1)
         YMIN=MIN(YMIN,YMIN1)
         YMAX=MAX(YMAX,YMAX1)
  100 CONTINUE
      CALL GQSCAL(XMIN,XMAX,GXMIN,GXMAX,GXSCAL)
      CALL GQSCAL(YMIN,YMAX,GYMIN,GYMAX,GYSCAL)
C
C     Origin should have a scale mark
C
      IF(GXMIN*GXMAX.GT.0.0) THEN
         GXORG=GXMIN
      ELSE
         GXORG=0.0
      ENDIF
      IF(GYMIN*GYMAX.GT.0.0) THEN
         GYORG=GYMIN
      ELSE
         GYORG=0.0
      ENDIF
C
C     symbol font has better minus sign
C
      CALL SETFNT(2)
      CALL SETCHS(0.35,0.0)
      CALL GDEFIN(PXMIN,PXMAX,PYMIN,PYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,0,4)
      CALL GFRAME
C
C     NGULEN choose appropriate format of scale values
C
      CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.3,9)
      CALL GVALUE(GXORG,2*GXSCAL,0.0,0.0,NGULEN(2*GXSCAL))
      CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.3,9)
      CALL GVALUE(0.0,0.0,GYORG,2*GYSCAL,NGULEN(2*GYSCAL))
      DO 200 NG=1,NGMAX
         CALL SETLIN(-1,-1,7-MOD(NG-1,5))
         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,0)
C         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,MOD(NG-1,8))
  200 CONTINUE
      CALL SETLIN(-1,-1,7)
C
      CALL GTEXTX(POS,17.3,TXT,2)
C
      RETURN
      END
C
C
      SUBROUTINE GRAPH3(XL,YL,ZL,GX,GY,GZ,NXM,NXMAX)
C
      DIMENSION GX(NXM),GY(NXM),GZ(NXM,NXM)
      EXTERNAL R2G2B
C
      PHI=-90.1
      THETA=0.01
      RADIUS=9.0
C
      CALL GMNMX1(GX,1,NXMAX,1,XMIN,XMAX)
      CALL GMNMX1(GY,1,NXMAX,1,YMIN,YMAX)
      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NXMAX,1,ZMIN,ZMAX)
      IF (ZMIN.LT.PSIMIN) ZMIN=PSIMIN
      CALL GQSCAL(XMIN,XMAX,GXMIN,GXMAX,GXSCAL)
      CALL GQSCAL(YMIN,YMAX,GYMIN,GYMAX,GYSCAL)
      CALL GQSCAL(ZMIN,ZMAX,GZMIN,GZMAX,GZSCAL)
      OX = 0.5*(XMIN+XMAX)
      OY = 0.5*(YMIN+YMAX)
      OZ = 0.5*(ZMIN+ZMAX)
C
C     Origin should have a scale mark
C
C      IF(GXMIN*GXMAX.GT.0.0) THEN
C         GXORG=GXMIN
C      ELSE
C         GXORG=0.0
C      ENDIF
C      IF(GYMIN*GYMAX.GT.0.0) THEN
C         GYORG=GYMIN
C      ELSE
C         GYORG=0.0
C      ENDIF
C      IF(GZMIN*GZMAX.GT.0.0) THEN
C         GZORG=GZMIN
C      ELSE
C         GZORG=0.0
C      ENDIF
C
C     symbol font has better minus sign
C
      CALL PAGES
      CALL SETFNT(2)
      CALL SETCHS(0.35,0.0)
      CALL GDEFIN3D(XL,YL,ZL,GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX)
      CALL GVIEW3D(PHI,THETA,RADIUS,OX,OY,OZ)
      CALL SETLIN(0,0,7)
      CALL GFRAME
C
C     NGULEN choose appropriate format of scale values
C
      CALL GSCALE3DX(XMIN,GXSCAL,0.3,0)
      CALL GSCALE3DY(YMIN,GYSCAL,0.3,0)
      CALL GSCALE3DZ(ZMIN,GZSCAL,0.3,0)
      CALL GVALUE3DX(XMIN,GXSCAL,1,1)
      CALL GVALUE3DY(YMIN,GYSCAL,1,1)
      CALL GVALUE3DZ(ZMIN,GZSCAL,1,1)
      CALL SETLIN(-1,-1,7)
C
C     Gouraud Shading
      CALL PERS3D1(GZ,NXM,NXMAX,NXMAX,-27,R2G2B)
      CALL GAxis3D(0)
      CALL GDrawBack3D(0.3, 0.3, 0.3)
      CALL PAGEE
C
      RETURN
      END
