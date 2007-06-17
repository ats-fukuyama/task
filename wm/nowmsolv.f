C     $Id$
C
C     ****** SOLVE SIMULTANEOUS EQUATIONS ******
C
      SUBROUTINE WMSOLV
C
      INCLUDE 'wmcomm.inc'
C
      COMPLEX * 16 D,X
      COMPLEX * 16 P1,P2,Q1,Q2,R1,R2
      COMPLEX * 16 P,Q,R,W,R0,E,H
      COMPLEX * 16 TMPA,TMPG
      COMPLEX * 16 TMPB
      COMPLEX * 16 TEMP1,TEMP2
      COMPLEX * 16 TINV1
      REAL    *  8 EPS
C
      PARAMETER (NBSIZM=3*MDM*NDM)
C
      DIMENSION P1(MSIZP),P2(MSIZP)
      DIMENSION Q1(MSIZP),Q2(MSIZP)
      DIMENSION R1(MSIZP),R2(MSIZP)
      DIMENSION X(MSIZP),D(NBSIZM,MSIZP)
      DIMENSION P(MSIZP),Q(MSIZP),R(MSIZP)
      DIMENSION W(MSIZP)
      DIMENSION R0(MSIZP),E(MSIZP),H(MSIZP)
      DIMENSION TEMP1(NBSIZM,NBSIZM),TEMP2(NBSIZM)
      DIMENSION TMPA (NBSIZM,NBSIZM),TMPG (NBSIZM,NBSIZM)
      DIMENSION TMPB (NBSIZM)
      DIMENSION TINV1(NBSIZM,NBSIZM)
      DIMENSION NM(NBSIZM)
      DIMENSION CFVS(MSIZP)
C
      NBSIZ=3*MDSIZ*NDSIZ
      MSIZ=NRMAX*NBSIZ
      MBND=2*NBSIZ
C
      NFRAC=(NRMAX-1)/NPROCS+1
      NBST=MYRANK*NFRAC+1
      NBED=NBST+NFRAC-1
      IF(MYRANK.EQ.NPROCS-1) NBED=NRMAX
C
      DO I=1,MSIZ
         CFVS(I)=CFV(I)
      ENDDO
C
      CALL MTXSET(NPROCS,MYRANK,NBSIZ,NRMAX)
C
C      WRITE(6,*) 'MODELM,NPROCS=',MODELM,NPROCS
C      WRITE(6,*) 'NCPUMIN,NCPUMAX=',NCPUMIN,NCPUMAX
C      WRITE(6,*) 'MYRANK,ISTA,IEND=',MYRANK,ISTA,IEND
C
C
C     ######## SELECTION OF SOLUTION METHOD ########
C
      IF(MODELM.EQ.0.AND.MYRANK.EQ.0) THEN
         CALL BANDCDNB(CEM,CFV,MSIZ,2*MBND-1,MBNDM,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XXX BANDCD PIVOT ERROR: IERR=',IERR
         ENDIF
      ENDIF
C
      IF(MODELM.EQ.1.AND.MYRANK.EQ.0) THEN
         CALL BANDCDB(CEM,CFV,X,MSIZ,2*MBND-1,MBNDM,
     &        NBSIZM,NM,TMPA,TMPG,TMPB,TINV1,TEMP1,TEMP2,
     &        NBSIZ,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XXX BANDCDB PIVOT ERROR: IERR=',IERR
         ENDIF
      ENDIF
C
      IF(MODELM.EQ.2.AND.MYRANK.EQ.0) THEN
         CALL BSTABCDB(CEM,MSIZ,MBND-1,MBNDM,NBSIZM,CFV,
     &                 EPS,ITR,X,D,R1,R2,P1,P2,Q1,Q2,H,W,
     &                 TEMP1,TEMP2,NM,NBSIZ,IERR)
         IF(ITR.GT.3) 
     &        WRITE(6,*) 'ITR=',ITR,'  EPS=',EPS
      ENDIF
C
      IF(MYRANK.EQ.0) THEN
         DO I=1,MSIZ
            CFVG(I)=CFV(I)
         ENDDO
      ENDIF
C
      RETURN
      END
