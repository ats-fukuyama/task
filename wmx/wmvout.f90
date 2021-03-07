! wmvout.f90

MODULE wmvout_local
  USE wmcomm,ONLY: rkind
  REAL(rkind):: RMIN,RMAX,ZMIN,ZMAX
  REAL(rkind),ALLOCATABLE:: R(:),Z(:)
END MODULE wmvout_local

MODULE wmvout

  PRIVATE
  PUBLIC wm_grms
  
CONTAINS

!      ****** DRAW MAGNETIC SURFACE DATA ******

  SUBROUTINE wm_grms

!                 THIS PROGRAM - PROUT - ACCEPTS THE OUTPUT
!                 FROM THE EQUILIBRIUM PROGRAM VMEC (wout FILES)
!                 AND PRINTS/PLOTS THE APPROPRIATE DATA

    USE wmcomm
    USE vmcomm
    IMPLICIT NONE
    INTEGER:: MN

    DO MN=1,MNMAX
       IXM(MN)=NINT(XM(MN))
       IXN(MN)=NINT(XN(MN))
    ENDDO

    DO MN=1,NTOR0
       BMOD(MN)=1.5D0*BMOD(MN+MNMAX)-0.5D0*BMOD(MN+2*MNMAX)
       RGMOD(MN)=1.5D0*RGMOD(MN+MNMAX)-0.5D0*RGMOD(MN+2*MNMAX)
    ENDDO

    DO MN=1+NTOR0,MNMAX
       BMOD(MN)=0.D0
       RGMOD(MN)=0.D0
    ENDDO

    RJTHETA(1) = 2.D0*RJTHETA(2) - RJTHETA(3)
    RJZETA (1) = 2.D0*RJZETA (2) - RJZETA (3)
    RJTHETA(NSRMAX) = 2.D0*RJTHETA(NSRMAX-1) - RJTHETA(NSRMAX-2)
    RJZETA (NSRMAX) = 2.D0*RJZETA (NSRMAX-1) - RJZETA (NSRMAX-2)
    RIOTAS(1) = 2.D0*RIOTAS(2) - RIOTAS(3)

    CALL WRFCN
    CALL PLOTTER

    RETURN
  END SUBROUTINE wm_grms

!     ****** RADIAL BETA PROFILE AND |B|=BSQ ROUTINE ******

  SUBROUTINE WRFCN

    USE wmcomm
    USE vmcomm
    IMPLICIT NONE
    INTEGER:: J,LJ,MN,JP,I
    REAL(rkind):: TB1,TM1,TV1,TP1,UB1

!     ****** DETERMINE RADIAL BETA PROFILE AND |B|=BSQ ******

    DO J=1,NSRMAX
       LJ=(J-1)*MNMAX
       BETAP(J)=0.25D0*BMOD(MN0+LJ)**2

       DO MN=1,MNMAX
          BETAP(J)=0.25D0*BMOD(MN+LJ)**2+BETAP(J)
       ENDDO
    ENDDO

    DO J=1,NSRMAX
       PHIMOD(MN0+MNMAX*(J-1))=XSQRT(ABS(PHI(J)))
       PHIMOD2(MN0+MNMAX*(J-1))=ABS(PHI(J))
    ENDDO

    DO J=2,NSRMAX
       BETAP(J)=PRES(J)/BETAP(J)
    ENDDO

    BETAP(1)=1.5D0*BETAP(2)-0.5D0*BETAP(3)
    PRES (1)=1.5D0*PRES(2)-0.5D0*PRES(3)
    TB1     =1.5D0*BETAP(NSRMAX)-0.5D0*BETAP(NSRMAX-1)
    BZCO (1)=1.5D0*BZCO(2)-0.5D0*BZCO(3)
    BPCO (1)=0.D0

    DO J=2,NSRMAX
       UB(J)=VP(J)/PHIPS(J)
    ENDDO

    DO J=2,NSRMAX
       RMASS(J)=RMASS(J)/(ABS(PHIPS(J)))**RGAM
    ENDDO

    RMASS(1)=1.5D0*RMASS(2)-0.5D0*RMASS(3)
    VP   (1)=1.5D0*VP(2)-0.5D0*VP(3)
    TM1     =1.5D0*RMASS(NSRMAX)-0.5D0*RMASS(NSRMAX-1)
    TV1     =1.5D0*VP(NSRMAX)-0.5D0*VP(NSRMAX-1)
    TP1     =1.5D0*PRES(NSRMAX)-0.5D0*PRES(NSRMAX-1)
    UB1     =1.5D0*UB(NSRMAX)-0.5D0*UB(NSRMAX-1)
    UB   (1)=UB(2)*TWOPI

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

    UB(NSRMAX)=TWOPI*UB1
    VP(NSRMAX)=TV1
    RMASS(NSRMAX)=TM1
    BETAP(NSRMAX)=TB1
    PRES(NSRMAX)=TP1
    BPCO(NSRMAX)=2.D0*BPCO(NSRMAX)-BPCO(NSRMAX-1)
    BZCO(NSRMAX)=2.D0*BZCO(NSRMAX)-BZCO(NSRMAX-1)

    DO I=1,NSRMAX
       VP(I)=(TWOPI**2)*VP(I)
    ENDDO
    RETURN
  END SUBROUTINE WRFCN

!     ****** DATA PLOT ROUTINE ******

  SUBROUTINE PLOTTER

    USE wmcomm
    USE vmcomm
    USE wmvout_local
    IMPLICIT NONE
    REAL(rkind),ALLOCATABLE:: BSQ(:),R12(:),Z12(:),PHIM(:),PHIM2(:),RGSQRT(:)
    INTEGER,PARAMETER:: IVAR(4)=[0,1,2,3]
    INTEGER,PARAMETER:: NCON(4)=[10,15,10,10]
    INTEGER,PARAMETER:: KPLOT=2
    CHARACTER(LEN=100):: MCHAR
    INTEGER:: LES,KZ,KT,J,L,NDEG

    ALLOCATE(R(NRT),Z(NRT))
    ALLOCATE(BSQ(NRT),R12(NRT),Z12(NRT),PHIM(NRT),PHIM2(NRT),RGSQRT(NRT))

!     KPLOT: = 0, NO PLOTS; 
!            = 1, ONLY FLUX SURFACES AND MOD-B CONTOURS;
!            = 2, ALL PLOTS EXCEPT R, Z, LAMBDA COEFFICIENTS; 
!            = 3, ALL PLOTS

    IF (KPLOT.EQ.0) RETURN

!     ****** COMPUTE R,Z SCALES ******

    LES=1+NTHPTS*(NSRMAX-1)

    DO KZ=1,NHHMAX
       CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,R,Z,RMNC,ZMNS)
       IF(KZ.GT.1) GOTO 15
       RMAX=R(LES)
       RMIN=R(LES)
       ZMAX=Z(LES)
       ZMIN=Z(LES)
15     RMAX=DMAX1(RMAX,AMAXAF(R(LES),1,NTHPTS))
       RMIN=DMIN1(RMIN,AMINAF(R(LES),1,NTHPTS))
       ZMAX=DMAX1(ZMAX,AMAXAF(Z(LES),1,NTHPTS))
       ZMIN=DMIN1(ZMIN,AMINAF(Z(LES),1,NTHPTS))
    ENDDO
!---------------------------------------------------------------------------
    print *,' rmin = ',rmin,' rmax = ',rmax
    print *,' zmin = ',zmin,' zmax = ',zmax
!---------------------------------------------------------------------------

    DO KZ=1,NHHMAX

       CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,R,Z,RMNC,ZMNS)
       CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,BSQ,Z,BMOD,ZMNS)
       CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,PHIM,Z,PHIMOD,ZMNS)
       CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,PHIM2,Z,PHIMOD2,ZMNS)
       CALL TOTZ(NTHPTS,NSRMAX,KZ,XM,XN,RGSQRT,Z,RGMOD,ZMNS)

       DO KT=1,NTHPTS

          R12(KT)=R(KT)
          Z12(KT)=Z(KT)
          R12(KT+NTHPTS)=0.5D0*(SQRT(2.D0)-1.D0)*(R(KT+NTHPTS)-R(KT))
          Z12(KT+NTHPTS)=0.5D0*(SQRT(2.D0)-1.D0)*(Z(KT+NTHPTS)-Z(KT))
          DO J=2,NSRMAX
             L=KT+NTHPTS*(J-1)
             R12(L)=0.5D0*(R(L)+R(L-NTHPTS))
             Z12(L)=0.5D0*(Z(L)+Z(L-NTHPTS))
          ENDDO
       ENDDO

       NDEG=IDNINT(360.D0*(KZ-1)/DBLE(NHHMAX))
       WRITE(MCHAR,30)'   NDEG = ',NDEG
30     FORMAT('$',A,I7,'$')

       CALL CONTOUR(R,Z,PHIM,'$CONTOURS OF PHI$',MCHAR,NCON(1),IVAR(1))

       CALL CONTOUR(R12,Z12,BSQ,'$MOD-B CONTOURS$',MCHAR,NCON(2),IVAR(2))

       CALL CONTOUR(R,Z,PHIM2,'$PHI AND Q CONTOURS$',MCHAR,NCON(3),IVAR(3))

       CALL CONTOUR(R12,Z12,RGSQRT,'$JACOBIAN CONTOURS$',MCHAR,NCON(4),IVAR(4))

    ENDDO

    RETURN
  END SUBROUTINE PLOTTER

!     ****** COUTOUR PLOT ROUTINE ******

  SUBROUTINE CONTOUR(XPLOT,YPLOT,FUNC,ITOP,MTITLE,NCON,IVAR)

    USE wmcomm
    USE vmcomm
    USE wmvout_local
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XPLOT(NRT),YPLOT(NRT),FUNC(NTHPTS,NSD)
    CHARACTER(LEN=*),INTENT(IN)::  ITOP,MTITLE
    INTEGER,INTENT(IN):: NCON,IVAR
    REAL,ALLOCATABLE:: GXPLOT(:),GYPLOT(:),GRTHET(:),GZTHET(:)
    REAL,ALLOCATABLE:: GFUNC(:,:),GR1(:),GZ1(:)
    INTEGER,ALLOCATABLE:: KA1(:,:,:)
    REAL(rkind):: STEP,EXSTEP,FAC
    INTEGER:: NSTEP,NRAV,NZAV,I,NTHRS,J,KT
    REAL(rkind):: XMIN,XMAX,YMIN,YMAX,FMIN,FMAX,DELF
    REAL:: GXMIN,GXMAX,GYMIN,GYMAX,GRLEN,GZLEN,GPZ,GDELF
    REAL:: GXSMIN,GXSMAX,GYSMIN,GYSMAX,GPR,GXSTP,GYSTP,GFMIN
    EXTERNAL:: PAGES,MOVE,TEXTX,GQSCAL,GNUMBD,TEXT,NUmbD,PAGEE
    EXTERNAL:: GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,CONTP6
    
    ALLOCATE(GXPLOT(NRT),GYPLOT(NRT))
    ALLOCATE(GRTHET(NSD),GZTHET(NSD))
    ALLOCATE(GFUNC(NTHPTS,NSD),GR1(NRT),GZ1(NRT))
    ALLOCATE(KA1(4,NTHPTS,NSD))

!     ****** LOCATE PHYSICAL ORIGIN ******
!      ******  ROUND-OFF AXIS ******

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

    GXMIN = GUCLIP(XMIN)
    GXMAX = GUCLIP(XMAX)
    GYMIN = GUCLIP(YMIN)
    GYMAX = GUCLIP(YMAX)

    GRLEN=GUCLIP(RGMAX-RGMIN)
    GZLEN=GUCLIP(ZGMAX-ZGMIN)
    IF(GRLEN.GT.GZLEN) THEN
       GPR=15.0
       GPZ=15.0*GZLEN/GRLEN
    ELSE
       GPR=15.0*GRLEN/GZLEN
       GPZ=15.0
    ENDIF

    DO I=1,NSRMAX*NTHPTS
       GXPLOT(I)=GUCLIP(XPLOT(I))
       GYPLOT(I)=GUCLIP(YPLOT(I))
    ENDDO

    CALL PAGES

    CALL MOVE(3.0,17.2)
    CALL TEXTX(ITOP)
    CALL TEXTX(MTITLE)
    CALL GQSCAL(GXMIN,GXMAX,GXSMIN,GXSMAX,GXSTP)
    CALL GQSCAL(GYMIN,GYMAX,GYSMIN,GYSMAX,GYSTP)

    CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ, &
                GUCLIP(RGMIN),GUCLIP(RGMAX), &
                GUCLIP(ZGMIN),GUCLIP(ZGMAX))
    CALL GFRAME
    CALL GSCALE(GXMIN,GXSTP,0.0,0.0,0.2,9)
    CALL GSCALE(0.0,0.0,GYMIN,GYSTP,0.2,9)
    CALL GVALUE(GXSMIN,2*GXSTP,0.0,0.0,NGULEN(GXSTP))
    CALL GVALUE(0.0,0.0,GYSMIN,2*GYSTP,NGULEN(GYSTP))

!     ****** PLOT BOUNDARY AND MAGNETIC AXIS ******

    IF(IVAR.EQ.0) CALL GPLOTP(GXPLOT,GYPLOT,1,1,1,3,1,0)
    NTHRS=(NSRMAX-1)*NTHPTS
    CALL GPLOTP(GXPLOT(NTHRS+1),GYPLOT(NTHRS+1),1,NTHPTS,1,0,1,0)

!     ****** SETUP TO DRAW CONTOURS IN INDEX SPACE

    FMIN=AMINAF(FUNC,1,NSRMAX*NTHPTS)
    FMAX=AMAXAF(FUNC,1,NSRMAX*NTHPTS)
    DELF=(FMAX-FMIN)/DBLE(NCON)
    
    DO J=1,NSRMAX
       DO I=1,NTHPTS
          GFUNC(I,J)=GUCLIP(FUNC(I,J))
       ENDDO
    ENDDO

    DO I=1,NSRMAX*NTHPTS
       GR1(I)=GUCLIP(R(I))
       GZ1(I)=GUCLIP(Z(I))
    ENDDO

    GDELF=GUCLIP(DELF)
    GFMIN=GUCLIP(FMIN)

    CALL CONTP6(GFUNC,GR1,GZ1,NTHPTS,NTHPTS,NSRMAX, &
                GFMIN,GDELF,NCON,0,KA1)

!
!     ****** SETUP TRANSFORMATION TO ORIGINAL AXES ******

    IF(IVAR.EQ.2) THEN
       DO KT=1,NTHPTS,2
          DO J=1,NSRMAX
             GRTHET(J)=GUCLIP(R(KT+NTHPTS*(J-1)))
             GZTHET(J)=GUCLIP(Z(KT+NTHPTS*(J-1)))
          ENDDO
          CALL GPLOTP(GRTHET,GZTHET,1,NSRMAX,1,0,1,2)
       ENDDO
    ENDIF

    CALL MOVE(19.0,17.0)
    CALL TEXT('MAX: ',5)
    CALL NUMBD(FMAX,'(1PE12.4)',12)
    CALL MOVE(19.0,16.5)
    CALL TEXT('MIN: ',5)
    CALL NUMBD(FMIN,'(1PE12.4)',12)
    CALL MOVE(19.0,16.0)
    CALL TEXT('STP: ',5)
    CALL NUMBD(DELF,'(1PE12.4)',12)

    CALL PAGEE

    RETURN
  END SUBROUTINE CONTOUR

!     ****** CALCULATE R,Z VALUE ******

  SUBROUTINE TOTZ(NTHMAX_L,NRMAX_L,KZ,XM,XN,R,Z,RMNC,ZMNS)

    USE wmcomm
    USE vmcomm,ONLY: MNMAX
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NTHMAX_L,NRMAX_L,KZ
    REAL(rkind),INTENT(IN):: XM(MNMAX),XN(MNMAX)
    REAL(rkind),INTENT(OUT):: R(NTHMAX_L*NRMAX_L),Z(NTHMAX_L*NRMAX_L)
    REAL(rkind),INTENT(IN):: RMNC(MNMAX),ZMNS(MNMAX)
    INTEGER:: KT,NTOT,N,NR,JES,MES,MN,L,M,IARG
    REAL(rkind):: DTH,DPH,XM0,XN0,ARG

    DTH=1.D0/(NTHMAX_L-1.D0)
    DPH=1.0d0/dble(nhc*NHHMAX)
    NTOT=NRMAX_L*NTHMAX_L

    DO N=1,NTOT
       R(N)=0.D0
       Z(N)=0.D0
    ENDDO

    DO NR=1,NRMAX_L
       JES=NTHMAX_L*(NR-1)
       MES=MNMAX*(NR-1)
       DO  MN=1,MNMAX
          XM0=XM(MN)*DTH
          XN0=XN(MN)*DPH
          DO KT=1,NTHMAX_L
             L=KT+JES
             M=MN+MES
             ARG=XM0*(KT-1)-XN0*(KZ-1)
             IARG=INT(ARG)
             ARG=TWOPI*(ARG-IARG)
             R(L)=R(L)+RMNC(M)*COS(ARG)
             Z(L)=Z(L)+ZMNS(M)*SIN(ARG)
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE TOTZ

! --- max in a range of argument ---
  
  FUNCTION AMAXAF(R,NS,NE)
    USE wmcomm,ONLY: rkind
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NS,NE
    REAL(rkind),INTENT(IN):: R(NE)
    REAL(rkind):: AMAXAF
    REAL(rkind):: RMAX
    INTEGER:: N

    RMAX=R(NS)
    DO N=NS+1,NE
       RMAX=DMAX1(RMAX,R(N))
    ENDDO
    AMAXAF=RMAX

    RETURN
  END FUNCTION AMAXAF

! --- min in a range of argument ---
  
  FUNCTION AMINAF(R,NS,NE)
    USE wmcomm,ONLY: rkind
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NS,NE
    REAL(rkind),INTENT(IN):: R(NE)
    REAL(rkind):: AMINAF
    REAL(rkind):: RMIN
    INTEGER:: N

    RMIN=R(NS)
    DO N=NS+1,NE
       RMIN=DMIN1(RMIN,R(N))
    ENDDO

    AMINAF=RMIN

    RETURN
  END FUNCTION AMINAF

! --- square root for positive value ---
  
  FUNCTION XSQRT(X)
    USE wmcomm,ONLY: rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: XSQRT
    
    IF(X.LE.0.D0) THEN
       XSQRT=0.D0
    ELSE
       XSQRT=SQRT(X)
    ENDIF
    RETURN
  END FUNCTION XSQRT
END MODULE wmvout
