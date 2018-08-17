C     $Id$
C
C     ***** CALCULATE MAGNETIC AXIS AND EDGE *****
C
      SUBROUTINE EQAXIS(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      REAL(8),DIMENSION(:,:),ALLOCATABLE::  PSIRG,PSIZG,PSIRZG
      EXTERNAL PSIGD,PSIGZ0
C
      ALLOCATE(PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM))
      IERR=0
C
C     ----- calculate setup for psig(R,Z) -----
C
      CALL setup_psig
C
C     ----- calculate position of magnetic axis -----
C
      CALL find_axis
C
C      WRITE(6,*) RAXIS,ZAXIS,PSIG(RAXIS,ZAXIS)
C
      IF(MDLEQF.LT.10) THEN
         RMAX=RR+RB
         RMIN=RR-RB
         ZMAX= RKAP*RB
         ZMIN=-RKAP*RB
      ELSE
         RMAX=RGMAX
         RMIN=RGMIN
         ZMAX=ZGMAX
         ZMIN=ZGMIN
      ENDIF
C
C      write(6,'(1P4E12.4)') RMAX,RMIN,ZMAX,ZMIN
C
      IF(RAXIS.LE.RMAX.AND.
     &   RAXIS.GE.RMIN.AND.
     &   ZAXIS.LE.ZMAX.AND.
     &   ZAXIS.GE.ZMIN) THEN
         PSI0=PSIG(RAXIS,ZAXIS)
         PSIPA=-PSI0
      ELSE
         WRITE(6,'(A)') 'XX EQAXIS: AXIS OUT OF PLASMA:'
         IERR=103
         RETURN
      ENDIF
C
C     ----- calculate outer plasma surface -----
C
      IF(PSIGZ0(RR)*PSIGZ0(RMAX).GE.0.D0) THEN
         REDGE=RMAX
      ELSE
         REDGE=FBRENT(PSIGZ0,RR,RMAX,1.D-8)
      ENDIF
C      write(6,*) 'redge=',redge
C
      DEALLOCATE(PSIRG,PSIZG,PSIRZG)
      RETURN
      END
C
C     ***** INTEGRATE ALONG THE MAGNETIC FIELD LINE *****
C
      SUBROUTINE EQMAGS(RINIT,ZINIT,NMAX,XA,YA,N,IERR)
C
C     ** Input **
C       RINIT : Initial starting point for tracing
C       ZINIT : Initial starting point for tracing
C       NMAX  : Size of arraies of XA, YA
C     ** Output **
C       XA    : Length from (RINIT,ZINIT) to the current position along the field line
C       YA    : Coordinate of the current position
C       N     : Number of partitions along the magnetic surface
C       IERR  : Error indicator
C
      INCLUDE '../eq/eqcomc.inc'
C
      EXTERNAL EQDERV
      DIMENSION Y(2),DYDX(2),YOUT(2)
      DIMENSION XA(NMAX),YA(2,NMAX)
C
      NEQ=2
C
      FACT=SQRT(SQRT(2.D0))
      H=FACT*2.D0*PI*RKAP*(RINIT-RAXIS)/NMAX
      ISTEP=0
C
C      WRITE(6,'(I5,1P3E12.4)') 0,H,RAXIS,ZAXIS
C      WRITE(6,'(I5,1P3E12.4)') NMAX,FACT,PI,RKAP
C      pause
C
  100 X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
C
C      WRITE(6,'(I5,1P3E12.4)') 1,X,Y(1),Y(2)
C
      N=1
      XA(N)=X
      YA(1,N)=Y(1)
      YA(2,N)=Y(2)
C
      IMODE=0
      DO I=2,NMAX
         CALL EQDERV(X,Y,DYDX)
         CALL EQRK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
         IF(IMODE.EQ.0) THEN
            IF((YOUT(2)-ZINIT)*PSI0.GT.0.D0) IMODE=1
         ELSE
            IF((YOUT(2)-ZINIT)*PSI0.LT.0.D0) GOTO 1000
         ENDIF
         X=X+H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
C
C         WRITE(6,'(I5,1P5E12.4)') N+1,X,Y(1),Y(2),DYDX(1),DYDX(2)
C
         N=N+1
         XA(N)=X
         YA(1,N)=Y(1)
         YA(2,N)=Y(2)
      ENDDO
C
      IF(ISTEP.LE.4) THEN
         H=FACT*H
         ISTEP=ISTEP+1
         GOTO 100
      ENDIF
      WRITE(6,*) 'XX EQMAGS: NOT ENOUGH N'
      pause
      IERR=1
      RETURN
C
 1000 CONTINUE
      H=0.1D0*H
      DO I=1,11
         CALL EQDERV(X,Y,DYDX)
         CALL EQRK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
         IF((YOUT(2)-ZINIT)*PSI0.LT.0.D0) GOTO 2000
         X=X+H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
      ENDDO
      WRITE(6,*) 'XX EQMAGS: UNEXPECTED BEHAVIOR'
      IERR=2
      RETURN
C
 2000 CONTINUE
      DEL=(ZINIT-Y(2))/(YOUT(2)-Y(2))
      X=X+H*DEL
      Y(1)=Y(1)+(YOUT(1)-Y(1))*DEL
      Y(2)=Y(2)+(YOUT(2)-Y(2))*DEL
      N=N+1
      XA(N)=X
      YA(1,N)=Y(1)
      YA(2,N)=Y(2)
      IERR=0
C
      RETURN
      END
C
C     ***** DERIVATIVES *****
C
      SUBROUTINE EQDERV(X,Y,DYDX)
C
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION Y(2),DYDX(2)
C
      CALL PSIGD(Y(1),Y(2),PSIRL,PSIZL)
C
      PSID=SQRT(PSIRL**2+PSIZL**2)
C
      DYDX(1)=-PSIZL/PSID
      DYDX(2)= PSIRL/PSID
C      WRITE(6,'(1P5E12.4)') X,Y(1),Y(2),DYDX(1),DYDX(2)
      RETURN
      END
C
C     ***** SETUP PSIG *****
C
      SUBROUTINE setup_psig
C
      INCLUDE '../eq/eqcomc.inc'
C
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: PSIRG,PSIZG,PSIRZG
C
C     ----- calculate spline coef for psi(R,Z) -----
C
      ALLOCATE(PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM))

      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,UPSIRZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX setup_psig: SPL2D ERROR: IERR=',IERR
         STOP
      ENDIF
      DEALLOCATE(PSIRG,PSIZG,PSIRZG)
      RETURN
      END
C
C     ***** calculate position of magnetic axis *****
C
      SUBROUTINE find_axis
C
      INCLUDE '../eq/eqcomc.inc'
      EXTERNAL PSIGD
C
      DELT=1.D-8
      EPS=1.D-4
      ILMAX=40
      LIST=0
      RINIT=RAXIS
      ZINIT=ZAXIS
      RSAVE=RAXIS
      ZSAVE=ZAXIS
      CALL NEWTN(PSIGD,RINIT,ZINIT,RAXIS,ZAXIS,
     &           DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5,1P2E12.4)')
     &        'XX EQAXIS: NEWTN ERROR: IER=',IER,RSAVE,ZSAVE
         WRITE(6,'(A)') 'XX EQAXIS: AXIS NOT FOUND:'
         IERR=102
         RETURN
      ENDIF
      RETURN
      END
C
C     ***** calculate position of xpoint1 *****
C
      SUBROUTINE find_xpoint1
C
      INCLUDE '../eq/eqcomc.inc'
      EXTERNAL PSIGD
C
      DELT=1.D-8
      EPS=1.D-4
      ILMAX=40
      LIST=0
      RINIT=RXPNT1
      ZINIT=ZXPNT1
      RSAVE=RINIT
      ZSAVE=ZINIT
      CALL NEWTN(PSIGD,RINIT,ZINIT,RXPNT1,ZXPNT1,
     &           DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5,1P2E12.4)')
     &        'XX find_xpoint1: NEWTN ERROR: IER=',IER,RSAVE,ZSAVE
         WRITE(6,'(A)') 'XX xpint1 NOT FOUND:'
         IERR=102
         RETURN
      ENDIF
      RETURN
      END
C
C     ***** calculate position of xpoint2 *****
C
      SUBROUTINE find_xpoint2
C
      INCLUDE '../eq/eqcomc.inc'
      EXTERNAL PSIGD
C
      DELT=1.D-8
      EPS=1.D-6
      ILMAX=40
      LIST=0
      RINIT=RXPNT2
      ZINIT=ZXPNT2
      RSAVE=RINIT
      ZSAVE=ZINIT
      CALL NEWTN(PSIGD,RINIT,ZINIT,RXPNT2,ZXPNT2,
     &           DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5,1P2E12.4)')
     &        'XX find_xpoint2: NEWTN ERROR: IER=',IER,RSAVE,ZSAVE
         WRITE(6,'(A)') 'XX xpint1 NOT FOUND:'
         IERR=102
         RETURN
      ENDIF
      RETURN
      END
C
C     ***** INTEGRATE ALONG THE MAGNETIC FIELD LINE *****
C
      SUBROUTINE calc_separtrix(RINIT,ZINIT,RXP,ZXP,H,NMAX,
     &                          XA,RA,ZA,NTOT,IERR)
C
C     ** Input **
C       RINIT : Initial starting point for tracing
C       ZINIT : Initial starting point for tracing
C       RXP   : Location of xpoint
C       ZYP   : Location o f xpoint
C       H     : Step size
C       NMAX  : Size of arraies of XA, YA
C     ** Output **
C       XA    : Length along the field line from (RINIT,ZINIT)
C       YA    : Position of separatrix points
C       N     : Number of positions
C       IERR  : Error indicator
C
C      INCLUDE '../eq/eqcomc.inc'
C
      IMPLICIT NONE
      REAL(8),INTENT(IN):: RINIT,ZINIT,RXP,ZXP,H
      INTEGER,INTENT(IN):: NMAX
      REAL(8),INTENT(OUT)::XA(NMAX),RA(NMAX),ZA(NMAX)
      INTEGER,INTENT(OUT):: NTOT,IERR

      INTEGER,PARAMETER:: NEQ=2
      REAL(8):: XA1(NMAX),YA1(2,NMAX)
      REAL(8):: XA2(NMAX),YA2(2,NMAX)
      REAL(8):: X,Y(2),DYDX(2),YOUT(2),DISTANCE
      INTEGER:: N,I,N1,N2
      EXTERNAL EQDERV

      X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
      N=1
      XA1(N)=X
      YA1(1,N)=Y(1)
      YA1(2,N)=Y(2)
C
      DO I=2,NMAX
         CALL EQDERV(X,Y,DYDX)
         CALL EQRK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
         X=X+H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
         N=N+1
         XA1(N)=X
         YA1(1,N)=Y(1)
         YA1(2,N)=Y(2)
         DISTANCE=SQRT((Y(1)-RXP)**2+(Y(2)-ZXP)**2)
         IF(DISTANCE < 1.2D0*H) GOTO 1000
      ENDDO
      GOTO 9100

 1000 CONTINUE
      N1=N+1
      XA1(N1)=X+DISTANCE
      YA1(1,N1)=RXP
      YA1(2,N1)=ZXP

      X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
      N=1
      XA2(N)=X
      YA2(1,N)=Y(1)
      YA2(2,N)=Y(2)
C
      DO I=2,NMAX
         CALL EQDERV(X,Y,DYDX)
         CALL EQRK4(X,Y,DYDX,YOUT,-H,NEQ,EQDERV)
         X=X-H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
         N=N+1
         XA2(N)=X
         YA2(1,N)=Y(1)
         YA2(2,N)=Y(2)
         DISTANCE=SQRT((Y(1)-RXP)**2+(Y(2)-ZXP)**2)
         IF(DISTANCE < 1.2D0*H) GOTO 2000
      ENDDO
      GO TO 9200

 2000 CONTINUE
      N2=N+1
      XA2(N2)=X-DISTANCE
      YA2(1,N2)=RXP
      YA2(2,N2)=ZXP

      NTOT=N1+N2-1
      DO N=1,N2
         XA(N)=XA2(N2-N+1)
         RA(N)=YA2(1,N2-N+1)
         ZA(N)=YA2(2,N2-N+1)
      ENDDO
      DO N=N2+1,NTOT
         XA(N)=XA1(N-N2+1)
         RA(N)=YA1(1,N-N2+1)
         ZA(N)=YA1(2,N-N2+1)
      END DO
      IERR=0
      RETURN

 9100 WRITE(6,*) 'XX calc_separtrix: XP not found for positive angle'
      IERR=1
      RETURN
 9200 WRITE(6,*) 'XX calc_separtrix: XP not found for negative angle'
      IERR=1
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PSI(R,Z) *****
C
      FUNCTION PSIG(R,Z)
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DF(R,Z,PSIL,RG,ZG,UPSIRZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX PSIG: SPL2DF ERROR: IERR=',IERR
         WRITE(6,'(A,1P2E12.4)') '   R,Z=',R,Z
      ENDIF
      PSIG=PSIL
      RETURN
      END
C
C     ***** INTERPOLATE SUBROUTINE DPSIDR,DPSIDZ(R,Z) *****
C
      SUBROUTINE PSIGD(R,Z,DPSIDR,DPSIDZ)
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DD(R,Z,PSIL,DPSIDR,DPSIDZ,
     &            RG,ZG,UPSIRZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX PSIGD: SPL2DD ERROR: IERR=',IERR
         WRITE(6,'(A,1P2E12.4)') '   R,Z=',R,Z
      ENDIF
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PSI on ZAXIS *****
C
      FUNCTION PSIGZ0(R)
      INCLUDE '../eq/eqcomc.inc'
      REAL(8) R,PSIGZ0
      PSIGZ0=PSIG(R,ZAXIS)
      RETURN
      END
