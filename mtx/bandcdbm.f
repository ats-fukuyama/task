C     $Id$
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION IN MPI) ******
C
      SUBROUTINE BANDCDBM(A,D,X,N,L,LA,AG,LG,NG,XG,LD,F,NM,
     &                    TMPA,TMPG,TMPB,TINV1,TEMP1,TEMP2,
     &                    NBSIZ,NFST,NEND,IERR)
C
      INCLUDE 'mpif.h'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      INTEGER      ISTATUS(MPI_STATUS_SIZE)
      COMPLEX * 16 A (LA,NFST:NEND),D (NFST:NEND),X (NFST:NEND)
      COMPLEX * 16 AG(LG,NG       ),XG(NG       )
      COMPLEX * 16 TMPA (LD,NBSIZ)
      COMPLEX * 16 TMPG (LD,NBSIZ)
      COMPLEX * 16 TMPB (NBSIZ)
      COMPLEX * 16 TINV1(LD,NBSIZ)
      COMPLEX * 16 TEMP1(LD,NBSIZ)
      COMPLEX * 16 TEMP2(NBSIZ)
      COMPLEX * 16 F(6*NBSIZ*2*NBSIZ)
      INTEGER      NM(NBSIZ)
      INTEGER      IDISP(2)
C
      NRP=0
      DO J=ISTA,IEND
         D(J-ISTA+1)=(0.D0,0.D0)
         DO I=1,L
            A(I,J-ISTA+1)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,J-ISTA+1),D(J-ISTA+1),J,LA,NRP)
      ENDDO
C
      IAF=2*NBSIZ-NBSIZ
      IBF=2*NBSIZ
      ICF=2*NBSIZ+NBSIZ
C
      LMMM=  NBSIZ
      LMM =2*NBSIZ
      LMID=3*NBSIZ
      LMP =4*NBSIZ
      LMPP=5*NBSIZ
C
      LBND=3*NBSIZ
C
      DO J=1,NG
         DO I=1,2*LBND-1
            AG(I,J)=(0.D0,0.D0)
         ENDDO
         XG(J)=(0.D0,0.D0)
      ENDDO
C
      CALL SENDRCVC2(A,F,N,2*NBSIZ,LA,NFST,NEND)
      CALL SENDRCVC (D  ,N,2*NBSIZ   ,NFST,NEND)
C
      IF(MYRANK.NE.NPROCS-1) THEN
         NBEDD=NBED
      ELSE
         NBEDD=NBED-1
      ENDIF
C
C 
C     ###### NEW COEFFICIENT MATRIX FOR BAND GAUSSIAN ELIMINATION ######
C
C     ###### DOWNWORD SCANNING ######
C
      IF(MYRANK.NE.NPROCS-1) THEN
         NBSTD=NBST+1
      ENDIF
C
C     ###### SET INITIAL VALUE ######
C
      IF(MYRANK.NE.NPROCS-1) THEN
C
         IBLC=NBSIZ*(NBST-1)+1
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TINV1(I,J)=A(IBF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)
            ENDDO
         ENDDO
C     
         CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
         if(ierr.ne.0) write(6,*) '0 ierr=',ierr
C
         DO J=1,NBSIZ
C
            TMPB(J)=(0.D0,0.D0)
            DO K=1,NBSIZ
               TMPB(J)=TMPB(J)
     &                +TINV1(K,J)*D(IBLC+(K-1)-ISTA+1)
            ENDDO
C
            DO I=1,NBSIZ
               TMPA(I,J)=(0.D0,0.D0)
               TMPG(I,J)=(0.D0,0.D0)
C
               DO K=1,NBSIZ
                  TMPA(I,J)=TMPA(I,J)
     &                     -TINV1(K,J)
     &                     *A(ICF+(I-1)-(K-1),IBLC+(K-1)-ISTA+1)
               ENDDO
C
               IF(MYRANK.NE.0) THEN
                  DO K=1,NBSIZ
                     TMPG(I,J)=TMPG(I,J)
     &                        -TINV1(K,J)
     &                        *A(IAF+(I-1)-(K-1),IBLC+(K-1)-ISTA+1)
                  ENDDO
               ENDIF
C
            ENDDO
C
         ENDDO
C
         DO IBN=NBSTD,NBED
            IBLC=NBSIZ*(IBN-1)+1
C
            CALL CALNCF(A,LA,IBLC,IAF,IBF,TMPA,TINV1,
     &                  LD,NBSIZ,NFST,NEND,IERR)
C
C     ###### CALCULATE BETA ######
C
            DO J=1,NBSIZ
               TEMP2(J)=D(IBLC+(J-1)-ISTA+1)
               DO K=1,NBSIZ
                  TEMP2(J)=TEMP2(J)
     &                    -A(IAF+(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                    *TMPB(K)
               ENDDO
            ENDDO
C
            DO J=1,NBSIZ
               TMPB(J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TMPB(J)=TMPB(J)+TINV1(K,J)*TEMP2(K)
               ENDDO
            ENDDO
C
C     ###### CALCULATE ALPHA ######
C
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  TMPA(I,J)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     TMPA(I,J)=TMPA(I,J)
     &                        -TINV1(K,J)
     &                        *A(ICF+(I-1)-(K-1),IBLC+(K-1)-ISTA+1)
                  ENDDO
               ENDDO
            ENDDO
C
C     ###### CALCULATE GAMMA #####
C
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  TEMP1(I,J)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     TEMP1(I,J)=TEMP1(I,J)
     &                         -A(IAF+(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                         *TMPG(I,K)
                  ENDDO
               ENDDO
            ENDDO
C
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  TMPG(I,J)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     TMPG(I,J)=TMPG(I,J)+TINV1(K,J)*TEMP1(I,K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
C
C     ###### SET THE COEFFICIENT MATRIX FOR THE COMMON PROCESSOR 1 ######
C
C
         IF(MYRANK.EQ.0) THEN
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  AG(LMMM+(I-1)-(J-1),2*NBSIZ+J)=(0.D0,0.D0)
                  IF(I.EQ.J) THEN
                     AG(LMID+(I-1)-(J-1),2*NBSIZ+J)=(1.D0,0.D0)
                  ENDIF
                  AG(LMP +(I-1)-(J-1),2*NBSIZ+J)=-TMPA(I,J)
               ENDDO
C
               XG(J+2*NBSIZ)=TMPB(J)
            ENDDO
C
         ELSE
            IBLCGN=2*NBSIZ*(MYRANK-1)+3*NBSIZ+NBSIZ+1
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  AG(LMMM+(I-1)-(J-1),IBLCGN+(J-1))=-TMPG(I,J)
                  IF(I.EQ.J) THEN
                     AG(LMID+(I-1)-(J-1),IBLCGN+(J-1))=(1.D0,0.D0)
                  ENDIF
                  AG(LMP +(I-1)-(J-1),IBLCGN+(J-1))=-TMPA(I,J)
               ENDDO
C
               XG(IBLCGN+(J-1))=TMPB(J)
            ENDDO
C
         ENDIF
C
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
C     ###### UPWARD SCANNING ######
C
      IF(MYRANK.EQ.0) THEN
         NBSTD=NBST+1
         NBEDD=NBED-1
      ELSE
         NBSTD=NBST
         NBEDD=NBED-1
      ENDIF
C
C     ###### SET INITIAL VALUE ######
C
      IBLC=NBSIZ*(NBED-1)+1
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            TINV1(I,J)=A(IBF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)
         ENDDO
      ENDDO
C
      CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
      if(ierr.ne.0) write(6,*) '1 ierr=',ierr
C
      DO J=1,NBSIZ
C
         TMPB(J)=(0.D0,0.D0)
         DO K=1,NBSIZ
            TMPB(J)=TMPB(J)
     &             +TINV1(K,J)*D(IBLC+(K-1)-ISTA+1)
         ENDDO
C
         DO I=1,NBSIZ
C
            TMPA(I,J)=(0.D0,0.D0)
            DO K=1,NBSIZ
               TMPA(I,J)=TMPA(I,J)
     &                  -TINV1(K,J)
     &                  *A(IAF+(I-1)-(K-1),IBLC+(K-1)-ISTA+1)
            ENDDO
C
            TMPG(I,J)=(0.D0,0.D0)
            IF(MYRANK.NE.NPROCS-1) THEN
               DO K=1,NBSIZ
                  TMPG(I,J)=TMPG(I,J)
     &                     -TINV1(K,J)
     &                     *A(ICF+(I-1)-(K-1),IBLC+(K-1)-ISTA+1)
               ENDDO
            ENDIF
C
         ENDDO
C
      ENDDO
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            A(IAF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)=TMPA(I,J)
            A(IBF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)=TMPG(I,J)
         ENDDO
         D(IBLC+(J-1)-ISTA+1)=TMPB(J)
      ENDDO
C
C     ###### UPWARD SCANNING MAIN ######
C
      DO IBN=NBEDD,NBSTD,-1
         IBLC=NBSIZ*(IBN-1)+1
C
         CALL CALNCF(A,LA,IBLC,ICF,IBF,TMPA,TINV1,
     &               LD,NBSIZ,NFST,NEND,IERR)
C
C     ###### CALCULATE BETA ######
C
         DO J=1,NBSIZ
            TEMP2(J)=D(IBLC+(J-1)-ISTA+1)
            DO K=1,NBSIZ
               TEMP2(J)=TEMP2(J)
     &                 -A(ICF+(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                 *TMPB(K)
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            TMPB(J)=(0.D0,0.D0)
            DO K=1,NBSIZ
               TMPB(J)=TMPB(J)+TINV1(K,J)*TEMP2(K)
            ENDDO
         ENDDO
C
C     ###### CALCULATE ALPHA ######
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TMPA(I,J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TMPA(I,J)=TMPA(I,J)
     &                     -TINV1(K,J)
     &                     *A(IAF+(I-1)-(K-1),IBLC+(K-1)-ISTA+1)
               ENDDO
            ENDDO
         ENDDO
C
C     ###### CALCULATE GAMMA ######
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TEMP1(I,J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TEMP1(I,J)=TEMP1(I,J)
     &                      -A(ICF+(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                      *TMPG(I,K)
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TMPG(I,J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TMPG(I,J)=TMPG(I,J)+TINV1(K,J)*TEMP1(I,K)
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               A(IAF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)=TMPA(I,J)
               A(IBF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)=TMPG(I,J)
            ENDDO
            D(IBLC+(J-1)-ISTA+1)=TMPB(J)
         ENDDO
C
      ENDDO
C
C     ###### SET THE COEFFICIENT MATRIX FOR THE COMMON PROCESSOR 2 ######
C
      IF(NPROCS.NE.1) THEN
C
         IF(MYRANK.EQ.0) THEN
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  AG(LMID+(I-1)-(J-1),J      )=A(IBF+(I-1)-(J-1),J)
                  AG(LMP +(I-1)-(J-1),J      )=A(ICF+(I-1)-(J-1),J)
                  AG(LMM +(I-1)-(J-1),J+NBSIZ)=-TMPA(I,J)
                  IF(I.EQ.J) THEN
                     AG(LMID+(I-1)-(J-1),J+NBSIZ)=(1.D0,0.D0)
                  ENDIF
                  AG(LMPP+(I-1)-(J-1),J+NBSIZ)=-TMPG(I,J)
               ENDDO
C
               XG(J      )=D(J)
               XG(J+NBSIZ)=TMPB(J)
            ENDDO
C
         ELSE
            IBLCGN=2*NBSIZ*(MYRANK-1)+3*NBSIZ+1
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  AG(LMM +(I-1)-(J-1),IBLCGN+(J-1))=-TMPA(I,J)
                  IF(I.EQ.J) THEN
                     AG(LMID+(I-1)-(J-1),IBLCGN+(J-1))=(1.D0,0.D0)
                  ENDIF
                  AG(LMPP+(I-1)-(J-1),IBLCGN+(J-1))=-TMPG(I,J)
               ENDDO
C
               XG(IBLCGN+(J-1))=TMPB(J)
            ENDDO
C
         ENDIF
C 
         IF(MYRANK.NE.0) THEN
            IF(MYRANK.NE.NPROCS-1) THEN
               IGF=2*NBSIZ*(MYRANK-1)+3*NBSIZ+1
               IGL=IGF+2*NBSIZ-1
            ELSE
               IGF=2*NBSIZ*(MYRANK-1)+3*NBSIZ+1
               IGL=IGF+NBSIZ-1
            ENDIF
         ENDIF
C
C     ###### SEND AND RECV MATRIX AND VECTOR USING FOR SMALL BAND ######
C
         DO NRANK=1,NPROCS-1
            IF(MYRANK.EQ.NRANK) THEN
               IDISP(1)=IGF
               IDISP(2)=IGL
               CALL MPI_SEND(IDISP,2,MPI_INTEGER,0,1,
     &              MPI_COMM_WORLD,IERR)
            ELSE IF(MYRANK.EQ.0) THEN
               CALL MPI_RECV(IDISP,2,MPI_INTEGER,NRANK,1,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
               IGF=IDISP(1)
               IGL=IDISP(2)
            ENDIF
C
            IF(MYRANK.EQ.NRANK) THEN
               MN=0
               DO J=IGF,IGL
                  DO I=1,2*LBND-1
                     MN=MN+1
                     F(MN)=AG(I,J)
                  ENDDO
               ENDDO
               CALL MPI_SEND(F,MN,MPI_DOUBLE_COMPLEX,0,2,
     &              MPI_COMM_WORLD,IERR)
            ELSE IF(MYRANK.EQ.0) THEN
               MN=(2*LBND-1)*(IGL-IGF+1)
               CALL MPI_RECV(F,MN,MPI_DOUBLE_COMPLEX,NRANK,2,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
               MN=0
               DO J=IGF,IGL
                  DO I=1,2*LBND-1
                     MN=MN+1
                     AG(I,J)=F(MN)
                  ENDDO
               ENDDO
            ENDIF
C
            IF(MYRANK.EQ.NRANK) THEN
               CALL MPI_SEND(XG(IGF),IGL-IGF+1,MPI_DOUBLE_COMPLEX,0,3,
     &              MPI_COMM_WORLD,IERR)
            ELSE IF(MYRANK.EQ.0) THEN
               CALL MPI_RECV(XG(IGF),IGL-IGF+1,MPI_DOUBLE_COMPLEX,
     &              NRANK,3,MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
         ENDDO
C
C     ###### GAUSSIAN ELIMINATION ######
C     
         IF(MYRANK.EQ.0) THEN
            CALL BANDCDD(AG,XG,NG,2*LBND-1,LG,IERR)
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
         CALL MPI_BCAST(XG,NG,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
C
         IF(MYRANK.EQ.0) THEN
            DO I=1,NBSIZ
               X(I               -ISTA+1)=XG(I        )
               X(I+NBSIZ         -ISTA+1)=XG(I+  NBSIZ)
               X((NBED-1)*NBSIZ+I-ISTA+1)=XG(I+2*NBSIZ)
            ENDDO
C
         ELSE IF(MYRANK.NE.NPROCS-1) THEN
            IBLCG=2*NBSIZ*(MYRANK-1)+3*NBSIZ+1
            DO I=1,NBSIZ
               X((NBST-1)*NBSIZ+I-ISTA+1)=XG(IBLCG      +(I-1))
               X((NBED-1)*NBSIZ+I-ISTA+1)=XG(IBLCG+NBSIZ+(I-1))
            ENDDO
C
         ELSE
            IBLCG=2*NBSIZ*(MYRANK-1)+3*NBSIZ+1
            DO I=1,NBSIZ
               X((NBST-1)*NBSIZ+I-ISTA+1)=XG(IBLCG      +(I-1))
            ENDDO
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
C         CALL SENDRCVC(X,N,2*NBSIZ,NFST,NEND)
         CALL SENDRCVC(X,N,NBSIZ,NFST,NEND)
C
      ELSE
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TINV1(I,J)=A(IBF+(I-1)-(J-1),J)
               DO K=1,NBSIZ
                  TINV1(I,J)=TINV1(I,J)
     &                      +A(ICF+(K-1)-(J-1),J)*TMPA(I,K)
               ENDDO
            ENDDO
         ENDDO
C
         CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
         if(ierr.ne.0) write(6,*) '2 ierr=',ierr
C
         DO J=1,NBSIZ
            TEMP2(J)=D(J)
            DO K=1,NBSIZ
               TEMP2(J)=TEMP2(J)-A(ICF+(K-1)-(J-1),J)*TMPB(K)
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            X(J)=(0.D0,0.D0)
            DO K=1,NBSIZ
               X(J)=X(J)+TINV1(K,J)*TEMP2(K)
            ENDDO
         ENDDO
C
      ENDIF
C
C     ###### CALCULATE SOLUTION VECTOR ######
C
      IF(NPROCS.NE.1) THEN
C
         IF(MYRANK.EQ.0) THEN
            NBSTD=NBST+2
            NBEDD=NBED-1
         ELSE IF(MYRANK.NE.NPROCS-1) THEN
            NBSTD=NBST+1
            NBEDD=NBED-1
         ELSE
            NBSTD=NBST+1
            NBEDD=NBED
         ENDIF
C
      ELSE
C
         NBSTD=NBST+1
         NBEDD=NBED
C
      ENDIF
C
      IF(NPROCS.NE.1) THEN
C
         DO IBN=NBSTD,NBEDD
            IBLC=NBSIZ*(IBN-1)+1
C
            DO I=1,NBSIZ
               X(IBLC+(I-1)-ISTA+1)=D(IBLC+(I-1)-ISTA+1)
               DO K=1,NBSIZ
                  X(IBLC+(I-1)-ISTA+1)=X(IBLC+(I-1)-ISTA+1)
     &                 +A(IAF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *X(IBLC-NBSIZ+(K-1)-ISTA+1)
     &                 +A(IBF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *X(NBED*NBSIZ+K-ISTA+1)
               ENDDO
            ENDDO
C
         ENDDO
C
      ELSE
C
         DO IBN=NBSTD,NBEDD
            IBLC=NBSIZ*(IBN-1)+1
C
            DO I=1,NBSIZ
               X(IBLC+(I-1)-ISTA+1)=D(IBLC+(I-1)-ISTA+1)
               DO K=1,NBSIZ
                  X(IBLC+(I-1)-ISTA+1)=X(IBLC+(I-1)-ISTA+1)
     &                 +A(IAF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *X(IBLC-NBSIZ+(K-1)-ISTA+1)
               ENDDO
            ENDDO
C
         ENDDO
C
      ENDIF
C
      DO I=ISTA,IEND
         D(I-2*NBSIZ-ISTA+1)=X(I-ISTA+1)
C         write(20,*) i,x(i-ista+1)
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE NEXT BLOCK COEFFICIENT ******
C
      SUBROUTINE CALNCF(A,LA,IBLC,IAF,IBF,TEMP1,TINV1,
     &                  LD,NBSIZ,NFST,NEND,IERR)
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
C
      COMPLEX * 16 A(LA,NFST:NEND)
      COMPLEX * 16 TEMP1(LD,NBSIZ)
      COMPLEX * 16 TINV1(LD,NBSIZ)
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            TINV1(I,J)=A(IBF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)
            DO K=1,NBSIZ
               TINV1(I,J)=TINV1(I,J)
     &                   +A(IAF+(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                   *TEMP1(I,K)
            ENDDO
         ENDDO
      ENDDO
C
      CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
      if(ierr.ne.0) write(6,*) '3 ierr=',ierr
C
      RETURN
      END
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION) ******
C
      SUBROUTINE BANDCDD( A , X , N , L , LA , IERR )
C
C     COMPLEX * 16    A( LA , N ) , X( N ) , ATMP( LMAX ) , TEMP
      COMPLEX * 16    A( LA , N ) , X( N ) , TEMP
      REAL    *  8    EPS , ABS1 , ABS2
      DATA EPS/ 1.D-70 /
C
      IF( MOD(L,2) .EQ. 0 ) GOTO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1
C
      DO 30 K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO 30 I = 1 , LHMK
            LPMI = L+1-I
            DO 40 J = 2 , L
               A( J-1 , K ) = A( J , K )
   40       CONTINUE
            A( L    , K    ) = ( 0.D0 , 0.D0 )
            A( LPMI , NPMK ) = ( 0.D0 , 0.D0 )
   30 CONTINUE
C
C     DO 50 I = 1 , NM
C        IPIVOT = I
C        IP     = I+1
C        DO 60 K = IP , LH
C           IF( CDABS(A(1,K)) .GT. CDABS(A(1,IPIVOT)) ) IPIVOT=K
C  60    CONTINUE
C
C
      DO 50 I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = CDABS( A(1,IPIVOT) )
         DO 60 K = IP , LH
            ABS1 = CDABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
   60    CONTINUE
C
         IF( CDABS(A(1,IPIVOT)) .LT. EPS ) GOTO 9002
         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO 90 J = 1 , L
C              ATMP( J )          = A   ( J , I      )
               TEMP               = A   ( J , I      )
               A   ( J , I      ) = A   ( J , IPIVOT )
C              A   ( J , IPIVOT ) = ATMP( J )
               A   ( J , IPIVOT ) = TEMP
   90       CONTINUE
         END IF
C
         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP
C
         DO 120 J = 2 , L
            A( J , I ) = A( J , I ) * TEMP
  120    CONTINUE
C
         DO 130 K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            DO 140 J = 2 , L
               A( J-1 , K ) = A( J , K ) - A( J , I ) * TEMP
  140       CONTINUE
C
            A( L , K ) = ( 0.D0 , 0.D0 )
  130    CONTINUE
         IF( LH .LT. N ) LH = LH + 1
   50 CONTINUE
C
      IF( CDABS(A(1,N)) .LT. EPS ) GOTO 9002
C
      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO 160 I = 1 , NM
         K = N-I
         TEMP = ( 0.D0 , 0.D0 )
         DO 170 J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
  170    CONTINUE
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
  160 CONTINUE
C
      IERR = 0
      RETURN
C
 9000 IERR = 10000
      RETURN
 9002 IERR = 30000+I
      RETURN
C
      END
