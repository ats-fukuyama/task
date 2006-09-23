C     $Id$
C***********************************************************************
C  -- Iterative template routine --
C     Univ. of Tennessee and Oak Ridge National Laboratory
C     October 1, 1993
C     Details of this algorithm are described in "Templates for the 
C     Solution of Linear Systems: Building Blocks for Iterative 
C     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
C     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
C     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
C
C  Purpose
C  =======
C
C  BICGSTAB solves the linear system A*x = b using the
C  BiConjugate Gradient Stabilized iterative method with 
C  preconditioning.
C
C  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
C  For other measures, see the above reference.
C
***********************************************************************
C
      SUBROUTINE BSTABCDB(A,N,M,LA,LD,B,EPS,ITER,X,D,R,RT,P,PH,V,T,
     &                    S,SH,TEMP1,TEMP2,NM,NBSIZ,IERR)
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A,B
      COMPLEX * 16 X,D
      COMPLEX * 16 R,RT
      COMPLEX * 16 P,PH
      COMPLEX * 16 V,T
      COMPLEX * 16 S,SH
      COMPLEX * 16 RHO,RHO1,RHO2,RHO3,RHO4
      COMPLEX * 16 ALPHA,BETA,OMEGA
      COMPLEX * 16 TEMP1,TEMP2
      REAL    *  8 EPS
      REAL    *  8 RESID,TOL,RHOTOL,OMEGATOL
      REAL    *  8 BNRM2,DNORM2
C
      DIMENSION A(LA,N),B(N)
      DIMENSION X(N),D(LD,N)
      DIMENSION R(N),RT(N)
      DIMENSION P(N),PH(N)
      DIMENSION V(N),T(N)
      DIMENSION S(N),SH(N)
      DIMENSION TEMP1(LD,NBSIZ),TEMP2(NBSIZ)
      DIMENSION NM(NBSIZ)
      LOGICAL   FTFLG
C
      IERR = 0
      EPS  = 1.D-8
      ITRM = 2000
C
C     TEST THE INPUT PARAMETERS.
C
      IF(N.LT.0) THEN
         IERR = -1
      ENDIF
      IF ( IERR.NE.0 ) RETURN
C
C     GIVE THE LIMIT
C
      RESID=EPS
      TOL  =EPS
      FTFLG = .TRUE.
C
      RHOTOL= 1.D-32
      OMEGATOL= 1.D-32
C
C     MAIN ROUTINE
C
      DO I=1,N
         X(I)=(0.D0,0.D0)
      ENDDO
C
      DO J=1,N
         DO I=1,NBSIZ
            D(I,J)=(0.D0,0.D0)
         ENDDO
      ENDDO
C
      NRP=0
      DO J=1,N
         B(J)=(0.D0,0.D0)
         DO I=1,2*M+1
            A(I,J)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,J),B(J),J,LA,NRP)
      ENDDO
C
      ICF=2*NBSIZ-NBSIZ
      IDF=2*NBSIZ
      IEF=2*NBSIZ+NBSIZ
C      
      CALL CALCUDB(A,D,TEMP1,NM,LA,LD,N,M,NBSIZ,IERR)
C
      DO I=1,N
         R(I)=0.D0
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            R(I)=R(I)-A(M+1+J,I)*X(K)
         ENDDO
      ENDDO
C
      DO I=1,N
         R(I)=B(I)-R(I)
      ENDDO
C
      IF(DNORM2(N,R,1).LE.TOL) GO TO 900
C
      DO I=1,N
         RT(I)=R(I)
      ENDDO
C
      BNRM2 = DNORM2(N,B,1)
      IF(BNRM2.EQ.0.D0) BNRM2 = 1.D0
C
      ITER = 0
C
C     ******** BI CONJUGATE GRADIENT STABILIZED ITERATION ********
C
  200 ITER=ITER+1
C
      RHO=(0.D0,0.D0)
      DO I=1,N
         RHO=RHO+RT(I)*R(I)
      ENDDO
      IF(CDABS(RHO).LT.RHOTOL) GO TO 700
C
C     COMPUTE VECTOR P
C
      IF(ITER.GT.1) THEN
         BETA=(RHO/RHO1)*(ALPHA/OMEGA)
         DO I=1,N
            P(I)=R(I)+BETA*(P(I)-OMEGA*V(I))
         ENDDO
      ELSE
         DO I=1,N
            P(I)=R(I)
         ENDDO
      ENDIF
C
C     COMPUTE DIRECTION ADJUSTING VECTOR PHAT AND SCALAR ALPHA
C
C      CALL CALJACB(A,LA,N,M,PH,P,IERR)
C
      DO I=1,N
         PH(I)=P(I)
      ENDDO
      CALL CALLDUB(A,D,TEMP2,PH,LA,LD,N,M,NBSIZ,IERR)
C
      DO I=1,N
         V(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            V(I)=V(I)+A(M+1+J,I)*PH(K)
         ENDDO
      ENDDO
C
      RHO2=(0.D0,0.D0)
      DO I=1,N
         RHO2=RHO2+RT(I)*V(I)
      ENDDO
      ALPHA=RHO/RHO2
C
C     EARLY CHECK FOR TOLERANCE
C
      DO I=1,N
         R(I)=R(I)-ALPHA*V(I)
      ENDDO
C
      DO I=1,N
         S(I)=R(I)
      ENDDO
C
      IF(DNORM2(N,S,1).LE.TOL) THEN
         DO I=1,N
            X(I)=X(I)+ALPHA*PH(I)
         ENDDO
         RESID=DNORM2(N,S,1)/BNRM2
         GO TO 900
      ELSE
C
C     COMPUTE STABILIZER VECTOR SHAT AND SCALAR OMEGA
C
C         CALL CALJACB(A,LA,N,M,PH,P,IERR)
C
C
         CALL CALLDUB(A,D,TEMP2,SH,LA,LD,N,M,NBSIZ,IERR)
C
         DO I=1,N
            T(I)=(0.D0,0.D0)
            DO K=MAX(I-M,1),MIN(I+M,N)
               J=K-I
               T(I)=T(I)+A(M+1+J,I)*SH(K)
            ENDDO
         ENDDO
C
         RHO3=(0.D0,0.D0)
         RHO4=(0.D0,0.D0)
         DO I=1,N
            RHO3=RHO3+T(I)*S(I)
            RHO4=RHO4+T(I)*T(I)
         ENDDO
C
         OMEGA=RHO3/RHO4
C
         DO I=1,N
            X(I)=X(I)+ALPHA*PH(I)+OMEGA*SH(I)
         ENDDO
C
C     COMPUTE RESIDUAL R, CHECK FOR TOLERANCE
C
         DO I=1,N
            R(I)=R(I)-OMEGA*T(I)
         ENDDO
C
         IF(FTFLG) IERR = -1
C
         CALL CJUDGE(N,R,B,BNRM2,RESID,TOL,IERR)
C         WRITE(6,*) 'ITR,RESID=',ITER,RESID
C
         IF(IERR.EQ.1) GO TO 900
C
         IF(ITER.EQ.ITRM) THEN
            GO TO 800
         ENDIF
C
      ENDIF
C
      IF(CDABS(OMEGA).LT.OMEGATOL) THEN
         GO TO 700
      ELSE
         RHO1=RHO
         GO TO 200
      ENDIF
C
  800 CONTINUE
C
C     ITERATION FAILS
C
      IERR = 1
      EPS  = RESID
      WRITE(6,*) 'NO CONVERGENCE BICGSTAB'
C
      RETURN
C
  700 CONTINUE
C
C     SET BREAKDOWN FLAG
C
      IF(ABS(RHO).LT.RHOTOL) THEN
         IERR= -10
         EPS = RESID
         WRITE(6,*) 'RHO IS BELOW RHOTOL'
      ELSE IF(ABS(OMEGA).LT.OMEGATOL) THEN
         IERR=-11
         EPS = RESID
         WRITE(6,*) 'OMEGA IS BELOW OMEGATOL'
      ENDIF
C
      RETURN
C
  900 CONTINUE
C
C     ITERATION SUCCESSFUL AND RETURN
C
      DO I=1,N
         B(I)=X(I)
      ENDDO
C
      EPS=RESID
C
      IERR= 0
C
C     END OF BICGSTAB
C
      RETURN
      END
C
C     ######## CONVERGENCE JUGDGE ROUTINE ########
C
      SUBROUTINE CJUDGE(N,R,B,BNRM2,RESID,TOL,IERR)
C
      COMPLEX * 16 R(N),B(N)
      REAL    *  8 RESID,TOL
      REAL    *  8 DNORM2,BNRM2
C
      IF(IERR.EQ.-1) THEN
         BNRM2 = DNORM2(N,B,1)
         IF(BNRM2.EQ.0.D0) BNRM2=1.D0 
      ENDIF
C
      RESID = DNORM2(N,R,1)/BNRM2
C
      IERR = 0
      IF(RESID.LE.TOL) IERR = 1
C
      RETURN
      END
C
C     ######## CALCULATE NORM ########
C
      FUNCTION DNORM2(N,DX,INC)
C
      COMPLEX * 16 DX(N)
      REAL    *  8 CUTLO,CUTHI,HITEST,SUM,XMAX
      REAL    *  8 DNORM2
      REAL    *  8 ZERO,ONE
C
      DATA ZERO,ONE /0.D0,1.D0/
      DATA CUTLO,CUTHI /8.232D-11, 1.304D19/
C
      IF(N.LE.0.AND.INC.LE.0) THEN
         DNORM2=0.D0
         RETURN
      ENDIF
C
   10 ASSIGN 30 TO NEXT
C
      SUM = 0.D0
      I   = 1
      IX  = 1
C
C     BEGIN MAIN LOOP
C
   20 GO TO NEXT, (30,50,70,110)
C
   30 IF(CDABS(DX(1)).GT.CUTLO) THEN
         IX=1
         GO TO 85
      ENDIF
C
      ASSIGN 50 TO NEXT
      XMAX = 0.D0
C
C     PHASE 1.  SUM IS ZERO
C
   50 IF(DX(I).EQ.0.D0) GO TO 200
C
      IF(CDABS(DX(I)).GT.CUTLO) THEN
         IX=I
         GO TO 85
      ENDIF
C
C     PREPARE FOR PHASE 2.
C
      ASSIGN 70 TO NEXT
      GO TO 105
C
C     PREPARE FOR PHASE 4.
C
  100 CONTINUE
C
      IX = J
      ASSIGN 110 TO NEXT
      SUM=(SUM/DX(I))/DCONJG(DX(I))
C
  105 XMAX = CDABS(DX(I))
      GO TO 115
C
C     PHASE 2.  SUM IS SMALL.
C     SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF(CDABS(DX(I)).GT.CUTLO) THEN
         IX=I
         GO TO 75
      ENDIF
C
C     COMMON CODE FOR PHASES 2 AND 4.
C     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF(CDABS(DX(I)).LE.XMAX) GO TO 115
C
      SUM = 1.D0 + SUM * (XMAX / DX(I))*DCONJG(XMAX / DX(I))
      XMAX = CDABS(DX(I))
      GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)*DCONJG(DX(I)/XMAX)
C
      GO TO 200
C
C     PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
C   85 HITEST = CUTHI/FLOAT( 2*N )
   85 HITEST = CUTHI/ 2*N
C
C     PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J = IX,N
         IF(CDABS(DX(I)) .GE. HITEST) GO TO 100
         SUM = SUM + DX(I)*DCONJG(DX(I))
         I = I + INC
   95 CONTINUE
      DNORM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
C
      IX = IX + 1
      I = I + INC
      IF( IX .LE. N ) GO TO 20
C
C     END OF MAIN LOOP.
C     
C     COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNORM2 = XMAX * DSQRT(SUM)
  300 CONTINUE

      RETURN
      END
