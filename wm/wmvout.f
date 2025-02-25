C     $Id$
C
C      ****** DRAW MAGNETIC SURFACE DATA ******
C
      SUBROUTINE WMGRMS
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
*************
*                 THIS PROGRAM - PROUT - ACCEPTS THE OUTPUT
*                 FROM THE EQUILIBRIUM PROGRAM VMEC (wout FILES)
*                 AND PRINTS/PLOTS THE APPROPRIATE DATA
*************
C
CF      OPEN(UNIT=6,FILE='OUTB',STATUS='UNKNOWN')
C
      DO MN=1,MNMAX
         IXM(MN)=NINT(XM(MN))
         IXN(MN)=NINT(XN(MN))
      ENDDO
C
      DO MN=1,NTOR0
         BMOD(MN)=1.5D0*BMOD(MN+MNMAX)-0.5D0*BMOD(MN+2*MNMAX)
         RGMOD(MN)=1.5D0*RGMOD(MN+MNMAX)-0.5D0*RGMOD(MN+2*MNMAX)
      ENDDO
C
      DO MN=1+NTOR0,MNMAX
         BMOD(MN)=0.D0
         RGMOD(MN)=0.D0
      ENDDO
C
      RJTHETA(1) = 2.D0*RJTHETA(2) - RJTHETA(3)
      RJZETA (1) = 2.D0*RJZETA (2) - RJZETA (3)
      RJTHETA(NSRMAX) = 2.D0*RJTHETA(NSRMAX-1) - RJTHETA(NSRMAX-2)
      RJZETA (NSRMAX) = 2.D0*RJZETA (NSRMAX-1) - RJZETA (NSRMAX-2)
      RIOTAS(1) = 2.D0*RIOTAS(2) - RIOTAS(3)
C
      CALL WRFCN
      CALL PLOTTER
C
      RETURN
      END
C
C     ****** RADIAL BETA PROFILE AND |B|=BSQ ROUTINE ******
C
      SUBROUTINE WRFCN
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C     ****** DETERMINE RADIAL BETA PROFILE AND |B|=BSQ ******
C
      DO J=1,NSRMAX
         LJ=(J-1)*MNMAX
         BETAP(J)=0.25D0*BMOD(MN0+LJ)**2
C
         DO MN=1,MNMAX
            BETAP(J)=0.25D0*BMOD(MN+LJ)**2+BETAP(J)
C            IF(J.LE.2) WRITE(6,*) LJ,MN,BMOD(MN+LJ),BETAP(J)
         ENDDO
      ENDDO
C
      DO J=1,NSRMAX
C         WRITE(6,'(A,1P3E12.4)') 
C     &        'PHI,PRES,BETAP=',PHI(J),PRES(J),BETAP(J)
         PHIMOD(MN0+MNMAX*(J-1))=XSQRT(ABS(PHI(J)))
         PHIMOD2(MN0+MNMAX*(J-1))=ABS(PHI(J))
      ENDDO
C
      DO J=2,NSRMAX
         BETAP(J)=PRES(J)/BETAP(J)
      ENDDO
C
      BETAP(1)=1.5D0*BETAP(2)-0.5D0*BETAP(3)
      PRES (1)=1.5D0*PRES(2)-0.5D0*PRES(3)
      TB1     =1.5D0*BETAP(NSRMAX)-0.5D0*BETAP(NSRMAX-1)
      BZCO (1)=1.5D0*BZCO(2)-0.5D0*BZCO(3)
      BPCO (1)=0.D0
C
      DO J=2,NSRMAX
         UB(J)=VP(J)/PHIPS(J)
      ENDDO
C
      DO J=2,NSRMAX
         RMASS(J)=RMASS(J)/(ABS(PHIPS(J)))**RGAM
      ENDDO
C
      RMASS(1)=1.5D0*RMASS(2)-0.5D0*RMASS(3)
      VP   (1)=1.5D0*VP(2)-0.5D0*VP(3)
      TM1     =1.5D0*RMASS(NSRMAX)-0.5D0*RMASS(NSRMAX-1)
      TV1     =1.5D0*VP(NSRMAX)-0.5D0*VP(NSRMAX-1)
      TP1     =1.5D0*PRES(NSRMAX)-0.5D0*PRES(NSRMAX-1)
      UB1     =1.5D0*UB(NSRMAX)-0.5D0*UB(NSRMAX-1)
      UB   (1)=UB(2)*TWOPI
C
      DO J=2,NSRMAX-1
         JP=J+1
         BPCO(J)=0.5D0*(BPCO(JP)+BPCO(J))
         UB(J)=0.5D0*(UB(JP)+UB(J))*TWOPI
         VP(J)=0.5D0*(VP(JP)+VP(J))
         RMASS(J)=0.5D0*(RMASS(JP)+RMASS(J))
         PRES(J)=0.5D0*(PRES(JP)+PRES(J))
         BETAP(J)=0.5D0*(BETAP(JP)+BETAP(J))
         BZCO(J)=0.5D0*(BZCO(JP)+BZCO(J))
      ENDDO
C
      UB(NSRMAX)=TWOPI*UB1
      VP(NSRMAX)=TV1
      RMASS(NSRMAX)=TM1
      BETAP(NSRMAX)=TB1
      PRES(NSRMAX)=TP1
      BPCO(NSRMAX)=2.D0*BPCO(NSRMAX)-BPCO(NSRMAX-1)
      BZCO(NSRMAX)=2.D0*BZCO(NSRMAX)-BZCO(NSRMAX-1)
C
      DO I=1,NSRMAX
         VP(I)=(TWOPI**2)*VP(I)
      ENDDO
C
C      PRINT 250, BETAP(1)
C      WRITE(3,250)BETAP(1)
C
C  200 FORMAT(//,40X,'FOURIER COEFFICIENTS X(M,N)',/)
C  210 FORMAT(//,9X,' S ',3X,6(7X,A2,I1,',',I3,')'),/)
C  220 FORMAT(1P7E15.3)
C  230 FORMAT(//,25X,'COVARIANT COMPONENTS OF B',
C     &            ' AND INTEGRATED CURRENTS',//,
C     &  9X,' S ',10X,'<BZETA>',8X,'<BTHETA>',8X,
C     &                    'RITOL',11X,'RIPOL',/)
C  240 FORMAT(//,9X,' S ',11X,'VP',12X,'DV/DPHI',10X,'RMASS',13X,
C     & 'P',12X,'BETA',/)
C  250 FORMAT(//,' BETA ON AXIS (SUM OVER MODES) = ',1PE10.3)
C
      RETURN
      END
C
C     ****** DATA PLOT ROUTINE ******
C
      SUBROUTINE PLOTTER
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      COMMON /PLTCN2/ RMIN,RMAX,ZMIN,ZMAX
C      COMMON /PLTCN6/ LXC,LYC,NOXC,NOYC
      COMMON /RZARRAY/R(NRT),Z(NRT)
C
      REAL*8 BSQ(NRT),R12(NRT),Z12(NRT),
     &       PHIM(NRT),PHIM2(NRT),RGSQRT(NRT)
      INTEGER IVAR(4),NCON(4)
      CHARACTER MCHAR*100
C 
      DATA IVAR / 0,1,2,3 / ,  NCON / 10,15,10,10 /,
     &     KPLOT / 2 /
C     &     NOXC,NOYC / 2,2 / , 
C     &     LXC,LYC / -1,0 / ,
C
C
C     KPLOT: = 0, NO PLOTS; 
C            = 1, ONLY FLUX SURFACES AND MOD-B CONTOURS;
C            = 2, ALL PLOTS EXCEPT R, Z, LAMBDA COEFFICIENTS; 
C            = 3, ALL PLOTS
C
      IF (KPLOT.EQ.0) RETURN
C
C     ****** COMPUTE R,Z SCALES ******
C
C      NHHMAX=1
C
C      IF(NTOR.NE.1) NHHMAX=4
C
      LES=1+NTHPTS*(NSRMAX-1)
C
      DO KZ=1,NHHMAX
         CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,R,Z,RMNC,ZMNS)
         IF(KZ.GT.1) GOTO 15
         RMAX=R(LES)
         RMIN=R(LES)
         ZMAX=Z(LES)
         ZMIN=Z(LES)
   15    RMAX=DMAX1(RMAX,AMAXAF(R(LES),1,NTHPTS))
         RMIN=DMIN1(RMIN,AMINAF(R(LES),1,NTHPTS))
         ZMAX=DMAX1(ZMAX,AMAXAF(Z(LES),1,NTHPTS))
         ZMIN=DMIN1(ZMIN,AMINAF(Z(LES),1,NTHPTS))
      ENDDO
c---------------------------------------------------------------------------
      print *,' rmin = ',rmin,' rmax = ',rmax
      print *,' zmin = ',zmin,' zmax = ',zmax
c---------------------------------------------------------------------------
C
      DO KZ=1,NHHMAX
C
         CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,R,Z,RMNC,ZMNS)
         CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,BSQ,Z,BMOD,ZMNS)
         CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,PHIM,Z,PHIMOD,ZMNS)
         CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,PHIM2,Z,PHIMOD2,ZMNS)
         CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,RGSQRT,Z,RGMOD,ZMNS)
C
         DO KT=1,NTHPTS
C
            R12(KT)=R(KT)
            Z12(KT)=Z(KT)
            R12(KT+NTHPTS)=0.5D0*(SQRT(2.D0)-1.D0)
     &                    *(R(KT+NTHPTS)-R(KT))
            Z12(KT+NTHPTS)=0.5D0*(SQRT(2.D0)-1.D0)
     &                    *(Z(KT+NTHPTS)-Z(KT))
C
         DO J=2,NSRMAX
C
            L=KT+NTHPTS*(J-1)
            R12(L)=0.5D0*(R(L)+R(L-NTHPTS))
            Z12(L)=0.5D0*(Z(L)+Z(L-NTHPTS))
C
         ENDDO
         ENDDO
C
         NDEG=IDNINT(360.D0*(KZ-1)/DBLE(NHHMAX))
         WRITE(MCHAR,30)'   NDEG = ',NDEG
   30    FORMAT('$',A,I7,'$')
C
C        CALL PAGES
        CALL CONTOUR(R,Z,PHIM,'$CONTOURS OF PHI$',
     &  MCHAR,NCON(1),IVAR(1))
C        CALL PAGEE
C
C        CALL PAGES
        CALL CONTOUR(R12,Z12,BSQ,'$MOD-B CONTOURS$',
     &  MCHAR,NCON(2),IVAR(2))
C        CALL PAGEE
C
C        CALL PAGES
        CALL CONTOUR(R,Z,PHIM2,'$PHI AND Q CONTOURS$',
     &  MCHAR,NCON(3),IVAR(3))
C        CALL PAGEE
C
C        CALL PAGES
        CALL CONTOUR(R12,Z12,RGSQRT,'$JACOBIAN CONTOURS$',
     &  MCHAR,NCON(4),IVAR(4))
C        CALL PAGEE
C
      ENDDO
C
      RETURN
      END
C
C     ****** COUTOUR PLOT ROUTINE ******
C
      SUBROUTINE CONTOUR(XPLOT,YPLOT,FUNC,ITOP,MTITLE,NCON,
     &                   IVAR)
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      COMMON /PLTCN2/ RMIN,RMAX,ZMIN,ZMAX
      COMMON /RZARRAY/ R(NRT),Z(NRT)
C
      DIMENSION XPLOT(*),YPLOT(*),FUNC(NTHPTS,NSD)
      DIMENSION GXPLOT(NRT),GYPLOT(NRT)
      DIMENSION GRTHET(NSD),GZTHET(NSD)
      DIMENSION GFUNC(NTHPTS,NSD),GR1(NRT),GZ1(NRT)
C
      CHARACTER ITOP*(*),MTITLE*(*)
      DIMENSION KA1(4,NTHPTS,NSD)
C
C     ****** LOCATE PHYSICAL ORIGIN ******
C      ******  ROUND-OFF AXIS ******
C
      STEP = DMAX1(RMAX-RMIN,ZMAX-ZMIN)
      EXSTEP = 2.D0 - DNINT(DLOG10(STEP))
      FAC  = 10.D0**EXSTEP
      NSTEP = IDNINT(.25D0 * (FAC*STEP)*1.05D0)
      NRAV  = IDNINT(.50D0 * (RMAX + RMIN) * FAC)
      NZAV  = IDNINT(.50D0 * (ZMAX + ZMIN) * FAC)
      STEP  = NSTEP/FAC
      XMIN  = NRAV /FAC - 2.D0*STEP
      XMAX  = NRAV /FAC + 2.D0*STEP
      YMIN  = NZAV /FAC - 2.D0*STEP
      YMAX  = NZAV /FAC + 2.D0*STEP
C
      GXMIN = GUCLIP(XMIN)
      GXMAX = GUCLIP(XMAX)
      GYMIN = GUCLIP(YMIN)
      GYMAX = GUCLIP(YMAX)
C
      GRLEN=GUCLIP(RGMAX-RGMIN)
      GZLEN=GUCLIP(ZGMAX-ZGMIN)
      IF(GRLEN.GT.GZLEN) THEN
         GPR=15.0
         GPZ=15.0*GZLEN/GRLEN
      ELSE
         GPR=15.0*GRLEN/GZLEN
         GPZ=15.0
      ENDIF
C
      DO I=1,NSRMAX*NTHPTS
         GXPLOT(I)=GUCLIP(XPLOT(I))
         GYPLOT(I)=GUCLIP(YPLOT(I))
      ENDDO
C
      CALL PAGES
C
      CALL MOVE(3.0,17.2)
      CALL TEXTX(ITOP)
      CALL TEXTX(MTITLE)
      CALL GQSCAL(GXMIN,GXMAX,GXSMIN,GXSMAX,GXSTP)
      CALL GQSCAL(GYMIN,GYMAX,GYSMIN,GYSMAX,GYSTP)
C
C      CALL GDEFIN(3.0,23.0,2.0,17.0,GXSMIN,GXMAX,GYMIN,GYMAX)
      CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ,
     &               GUCLIP(RGMIN),GUCLIP(RGMAX),
     &               GUCLIP(ZGMIN),GUCLIP(ZGMAX))
      CALL GFRAME
      CALL GSCALE(GXMIN,GXSTP,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYMIN,GYSTP,0.2,9)
      CALL GVALUE(GXSMIN,2*GXSTP,0.0,0.0,NGULEN(GXSTP))
      CALL GVALUE(0.0,0.0,GYSMIN,2*GYSTP,NGULEN(GYSTP))
C
C     ****** PLOT BOUNDARY AND MAGNETIC AXIS ******
C
      IF(IVAR.EQ.0) CALL GPLOTP(GXPLOT,GYPLOT,1,1,1,3,1,0)
      NTHRS=(NSRMAX-1)*NTHPTS
      CALL GPLOTP(GXPLOT(NTHRS+1),GYPLOT(NTHRS+1),1,NTHPTS,1,0,1,0)
C
C     ****** SETUP TO DRAW CONTOURS IN INDEX SPACE
C
      FMIN=AMINAF(FUNC,1,NSRMAX*NTHPTS)
      FMAX=AMAXAF(FUNC,1,NSRMAX*NTHPTS)
      DELF=(FMAX-FMIN)/DBLE(NCON)
C
      DO J=1,NSRMAX
      DO I=1,NTHPTS
         GFUNC(I,J)=GUCLIP(FUNC(I,J))
      ENDDO
      ENDDO
C
      DO I=1,NSRMAX*NTHPTS
         GR1(I)=GUCLIP(R(I))
         GZ1(I)=GUCLIP(Z(I))
      ENDDO
C
      GDELF=GUCLIP(DELF)
      GFMIN=GUCLIP(FMIN)
C
      CALL CONTP6(GFUNC,GR1,GZ1,NTHPTS,NTHPTS,NSRMAX,
     &            GFMIN,GDELF,NCON,0,KA1)
C
C
C     ****** SETUP TRANSFORMATION TO ORIGINAL AXES ******
C
      IF(IVAR.EQ.2) THEN
         DO KT=1,NTHPTS,2
            DO J=1,NSRMAX
               GRTHET(J)=GUCLIP(R(KT+NTHPTS*(J-1)))
               GZTHET(J)=GUCLIP(Z(KT+NTHPTS*(J-1)))
            ENDDO
            CALL GPLOTP(GRTHET,GZTHET,1,NSRMAX,1,0,1,2)
         ENDDO
      ENDIF
C
      CALL MOVE(19.0,17.0)
      CALL TEXT('MAX: ',5)
      CALL NUMBD(FMAX,'(1PE12.4)',12)
      CALL MOVE(19.0,16.5)
      CALL TEXT('MIN: ',5)
      CALL NUMBD(FMIN,'(1PE12.4)',12)
      CALL MOVE(19.0,16.0)
      CALL TEXT('STP: ',5)
      CALL NUMBD(DELF,'(1PE12.4)',12)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ****** CALCULATE R,Z VALUE ******
C
      SUBROUTINE TOTZ(NTH,NS,KZ,XM,XN,
     &                R,Z,RMNC,ZMNS)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RMNC(*),ZMNS(*),R(*),Z(*),XM(*),XN(*)
      COMMON /WMHVR1/ NSRMAX,MNMAX
      COMMON /WMHRD1/ MN0,NTOR0
      COMMON /WMHRD2/ HS,OHS,DNORM,RC
C
      PIT=1.D0/(NTH-1.D0)
c.... updated on 05/05/03 by N^2 .. start ...........
c     PIZ=RC/NHHMAX
      PIZ=1.0d0/dble(nhc*NHHMAX)
c.... updated on 05/05/03 by N^2 .. end .............
      NRTH=NS*NTH
C
      DO L=1,NRTH
          R(L)=0.D0
          Z(L)=0.D0
      ENDDO
C
      DO J=1,NS
C
         JES=NTH*(J-1)
         MES=MNMAX*(J-1)
      DO  MN=1,MNMAX
         XM0=XM(MN)*PIT
         XN0=XN(MN)*PIZ
         
C
      DO KT=1,NTH
C
         L=KT+JES
         M=MN+MES
         ARG=XM0*(KT-1)-XN0*(KZ-1)
         IARG=IDINT(ARG)
         ARG=TWOPI*(ARG-IARG)
         R(L)=R(L)+RMNC(M)*COS(ARG)
         Z(L)=Z(L)+ZMNS(M)*SIN(ARG)
C
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
C
      FUNCTION AMAXAF(R,NS,NE)
C
      REAL*8 AMAXAF,R,RMAX
      DIMENSION R(NE)
C
      RMAX=R(NS)
      DO N=NS+1,NE
           RMAX=DMAX1(RMAX,R(N))
      ENDDO
C
      AMAXAF=RMAX
C
      RETURN
      END
C
      FUNCTION AMINAF(R,NS,NE)
C
      REAL*8 AMINAF,R,RMIN
      DIMENSION R(NE)
C
      RMIN=R(NS)
      DO N=NS+1,NE
         RMIN=DMIN1(RMIN,R(N))
      ENDDO
C
      AMINAF=RMIN
C
      RETURN
      END
C
      REAL*8 FUNCTION XSQRT(X)
      REAL*8 X
      IF(X.LE.0.D0) THEN
         XSQRT=0.D0
      ELSE
         XSQRT=SQRT(X)
      ENDIF
      RETURN
      END
