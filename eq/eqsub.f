C     $Id$
C
C     ***** INTEGRATE ALONG THE MAGNETIC FIELD LINE *****
C
      SUBROUTINE EQMAGS(RINIT,ZINIT,NMAX,XA,YA,N,IERR)
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
C      WRITE(6,'(I5,1P3E12.4)') 0,H,RAXIS,ZAXIS
C
  100 X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
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
            IF((YOUT(2)-ZINIT)*SAXIS.GT.0.D0) IMODE=1
         ELSE
            IF((YOUT(2)-ZINIT)*SAXIS.LT.0.D0) GOTO 1000
         ENDIF
         X=X+H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
C         WRITE(6,'(I5,1P3E12.4)') N+1,X,Y(1),Y(2)
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
      IERR=1
      RETURN
C
 1000 CONTINUE
      H=0.1D0*H
      DO I=1,11
         CALL EQDERV(X,Y,DYDX)
         CALL EQRK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
         IF((YOUT(2)-ZINIT)*SAXIS.LT.0.D0) GOTO 2000
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
      CALL EQPSID(Y(1),Y(2),PSIRL,PSIZL)
C
      PSID=SQRT(PSIRL**2+PSIZL**2)
C
      DYDX(1)=-PSIZL/PSID
      DYDX(2)= PSIRL/PSID
C      WRITE(6,'(1P5E12.4)') X,Y(1),Y(2),DYDX(1),DYDX(2)
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PSI(R,ZAXIS) *****
C
      FUNCTION PSIAX(R)
C
      INCLUDE '../eq/eqcomc.inc'
C
      PSIAX=PSIG(R,ZAXIS)
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PSI(R,Z) *****
C
      FUNCTION PSIG(R,Z)
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DF(R,Z,PSIL,RG,ZG,URZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX PSIG: SPL2DF ERROR : IERR=',IERR
      PSIG=PSIL
      RETURN
      END
C
C     ***** INTERPOLATE SUBROUTINE PSI,DPSIDR,DPSIDZ(R,Z) *****
C
      SUBROUTINE  EQPSID(R,Z,DPSIDR,DPSIDZ)
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DD(R,Z,PSIL,DPSIDR,DPSIDZ,
     &            RG,ZG,URZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQPSID: SPL2DD ERROR : IERR=',IERR
      IF(IERR.NE.0) WRITE(6,'(A,1P2E12.4)') '   R,Z=',R,Z
      RETURN
      END
