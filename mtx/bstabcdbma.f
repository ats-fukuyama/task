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
      SUBROUTINE BSTABCDBMA(A,N,M,LA,LD,B,EPS,ITER,X,D,R,RT,P,PH,V,T,
     &                    S,SH,F,TEMP1,TEMP2,NM,NBSIZ,NFST,NEND,IERR)
C
      INCLUDE 'mpif.h'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
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
      COMPLEX * 16 F
      COMPLEX * 16 CTMP
      REAL    *  8 EPS
      REAL    *  8 RESID,TOL,RHOTOL,OMEGATOL
      REAL    *  8 BNRM2,DNRM2M
      REAL    *  8 TMP
C
      DIMENSION A(LA,NFST:NEND),B(NFST:NEND)
      DIMENSION X(NFST:NEND),D(LD,NFST:NEND)
      DIMENSION R(NFST:NEND),RT(NFST:NEND)
      DIMENSION P(NFST:NEND),PH(NFST:NEND)
      DIMENSION V(NFST:NEND),T(NFST:NEND)
      DIMENSION S(NFST:NEND),SH(NFST:NEND)
      DIMENSION TEMP1(LD,NBSIZ),TEMP2(NBSIZ)
      DIMENSION F(M*(2*M+1))
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
      DO I=ISTA-M,IEND+M
         X (I-ISTA+1)=(0.D0,0.D0)
         R (I-ISTA+1)=(0.D0,0.D0)
         RT(I-ISTA+1)=(0.D0,0.D0)
         P (I-ISTA+1)=(0.D0,0.D0)
         PH(I-ISTA+1)=(0.D0,0.D0)
         V (I-ISTA+1)=(0.D0,0.D0)
         T (I-ISTA+1)=(0.D0,0.D0)
         S (I-ISTA+1)=(0.D0,0.D0)
         SH(I-ISTA+1)=(0.D0,0.D0)
      ENDDO
C
      DO J=ISTA-M,IEND+M
         DO I=1,NBSIZ
            D(I,J-ISTA+1)=(0.D0,0.D0)
         ENDDO
      ENDDO
C
      NRP=0
      DO J=ISTA,IEND
         B(J-ISTA+1)=(0.D0,0.D0)
         DO I=1,2*M+1
            A(I,J-ISTA+1)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,J-ISTA+1),B(J-ISTA+1),J,LA,NRP)
      ENDDO
C
      ISCNF=NBSIZ*(NBST-1)+1
      IECNF=NBSIZ*NBED
C
      ICF=2*NBSIZ-NBSIZ
      IDF=2*NBSIZ
      IEF=2*NBSIZ+NBSIZ
C
      CALL SENDRCVC2(A,F,N,M,LA,NFST,NEND)
C
C      CALL CALDINV(A,D,NM,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
      CALL CALCUDBM(A,D,TEMP1,NM,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
C      CALL CALCUDBMS(A,D,TEMP1,NM,LA,LD,N,M,NBSIZ,F,NFST,NEND,IERR)
C
      DO I=ISTA,IEND
         R(I-ISTA+1)=0.D0
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            R(I-ISTA+1)=R(I-ISTA+1)-A(M+1+J,I-ISTA+1)*X(K-ISTA+1)
         ENDDO
      ENDDO
C
      DO I=ISTA,IEND
         R(I-ISTA+1)=B(I-ISTA+1)-R(I-ISTA+1)
      ENDDO
C
      IF(DNRM2M(N,R,1,NFST,NEND).LE.TOL) GO TO 900
C
      DO I=ISTA,IEND
         RT(I-ISTA+1)=R(I-ISTA+1)
      ENDDO
C
      BNRM2 = DNRM2M(N,B,1,NFST,NEND)
      IF(BNRM2.EQ.0.D0) BNRM2 = 1.D0
C
      ITER = 0
C
C     ******** BI CONJUGATE GRADIENT STABILIZED ITERATION ********
C
  200 ITER=ITER+1
C
      RHO=(0.D0,0.D0)
      DO I=ISTA,IEND
         RHO=RHO+RT(I-ISTA+1)*R(I-ISTA+1)
      ENDDO
C
      CALL MPI_REDUCE(RHO,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      RHO=CTMP
      CALL MPI_BCAST(RHO,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)
C
      IF(CDABS(RHO).LT.RHOTOL) GO TO 700
C
C     COMPUTE VECTOR P
C
      IF(ITER.GT.1) THEN
         BETA=(RHO/RHO1)*(ALPHA/OMEGA)
         DO I=ISTA,IEND
            P(I-ISTA+1)=R(I-ISTA+1)
     &                 +BETA*(P(I-ISTA+1)-OMEGA*V(I-ISTA+1))
         ENDDO
      ELSE
         DO I=ISTA,IEND
            P(I-ISTA+1)=R(I-ISTA+1)
         ENDDO
      ENDIF
C
C     COMPUTE DIRECTION ADJUSTING VECTOR PHAT AND SCALAR ALPHA
C
C
      DO I=ISTA,IEND
         PH(I-ISTA+1)=P(I-ISTA+1)
      ENDDO
C
C      CALL CALJACBM(A,D,TEMP2,PH,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
      CALL CALLDUBM(A,D,TEMP2,PH,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
C      CALL CALLDUBMS(A,D,TEMP2,PH,LA,LD,N,M,NBSIZ,F,NFST,NEND,IERR)
C
      CALL SENDRCVC(PH,N,M,NFST,NEND)
C
      DO I=ISTA,IEND
         V(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            V(I-ISTA+1)=V(I-ISTA+1)+A(M+1+J,I-ISTA+1)*PH(K-ISTA+1)
         ENDDO
      ENDDO
C
      RHO2=(0.D0,0.D0)
      DO I=ISTA,IEND
         RHO2=RHO2+RT(I-ISTA+1)*V(I-ISTA+1)
      ENDDO
C
      CALL MPI_REDUCE(RHO2,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      RHO2=CTMP
      CALL MPI_BCAST(RHO2,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)
      ALPHA=RHO/RHO2
C
C     EARLY CHECK FOR TOLERANCE
C
      DO I=ISTA,IEND
         R(I-ISTA+1)=R(I-ISTA+1)-ALPHA*V(I-ISTA+1)
      ENDDO
C
      DO I=ISTA,IEND
         S(I-ISTA+1)=R(I-ISTA+1)
      ENDDO
C
      IF(DNRM2M(N,S,1,NFST,NEND).LE.TOL) THEN
         DO I=ISTA,IEND
            X(I-ISTA+1)=X(I-ISTA+1)+ALPHA*PH(I-ISTA+1)
         ENDDO
         RESID=DNRM2M(N,S,1,NFST,NEND)/BNRM2
         GO TO 900
      ELSE
C
C     COMPUTE STABILIZER VECTOR SHAT AND SCALAR OMEGA
C
C         CALL CALJACBM(A,LA,N,M,PH,P,NFST,NEND,IERR)
C
         DO I=ISTA,IEND
            SH(I-ISTA+1)=S(I-ISTA+1)
         ENDDO
C
C         CALL CALJACBM(A,D,TEMP2,SH,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
         CALL CALLDUBM(A,D,TEMP2,SH,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
C         CALL CALLDUBMS(A,D,TEMP2,SH,LA,LD,N,M,NBSIZ,F,NFST,NEND,IERR)
C
         CALL SENDRCVC(SH,N,M,NFST,NEND)
C
         DO I=ISTA,IEND
            T(I-ISTA+1)=(0.D0,0.D0)
            DO K=MAX(I-M,1),MIN(I+M,N)
               J=K-I
               T(I-ISTA+1)=T(I-ISTA+1)
     &                    +A(M+1+J,I-ISTA+1)*SH(K-ISTA+1)
            ENDDO
         ENDDO
C
         RHO3=(0.D0,0.D0)
         RHO4=(0.D0,0.D0)
         DO I=ISTA,IEND
            RHO3=RHO3+T(I-ISTA+1)*S(I-ISTA+1)
            RHO4=RHO4+T(I-ISTA+1)*T(I-ISTA+1)
         ENDDO
C
         CALL MPI_REDUCE(RHO3,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         RHO3=CTMP
         CALL MPI_BCAST(RHO3,1,MPI_DOUBLE_COMPLEX,0,
     &        MPI_COMM_WORLD,IERR)
C
         CALL MPI_REDUCE(RHO4,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &        MPI_COMM_WORLD,IERR)
         RHO4=CTMP
         CALL MPI_BCAST(RHO4,1,MPI_DOUBLE_COMPLEX,0,
     &        MPI_COMM_WORLD,IERR)
C
         OMEGA=RHO3/RHO4
C
         DO I=ISTA,IEND
            X(I-ISTA+1)=X(I-ISTA+1)+ALPHA*PH(I-ISTA+1)
     &                             +OMEGA*SH(I-ISTA+1)
         ENDDO
C
C     COMPUTE RESIDUAL R, CHECK FOR TOLERANCE
C
         DO I=ISTA,IEND
            R(I-ISTA+1)=R(I-ISTA+1)-OMEGA*T(I-ISTA+1)
         ENDDO
C
         IF(FTFLG) IERR = -1
C
         CALL CJUDGEM(N,R,B,BNRM2,RESID,TOL,NFST,NEND,IERR)
         IF(MYRANK.EQ.0) WRITE(6,*) 'ITR,RESID=',ITER,RESID
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
      IF(MYRANK.EQ.0) WRITE(6,*) 'NO CONVERGENCE BICGSTAB'
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
         IF(MYRANK.EQ.0) WRITE(6,*) 'RHO IS BELOW RHOTOL'
      ELSE IF(ABS(OMEGA).LT.OMEGATOL) THEN
         IERR=-11
         EPS = RESID
         IF(MYRANK.EQ.0) WRITE(6,*) 'OMEGA IS BELOW OMEGATOL'
      ENDIF
C
      RETURN
C
  900 CONTINUE
C
C     ITERATION SUCCESSFUL AND RETURN
C
      DO I=ISTA,IEND
C         B(I-M-ISTA+1)=X(I-ISTA+1)
          B(I-2*NBSIZ-ISTA+1)=X(I-ISTA+1)
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
