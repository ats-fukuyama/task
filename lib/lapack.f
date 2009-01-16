C     $Id$
      SUBROUTINE LAPACK_DGBTRF(M,N,KL,KU,AX,LDAB,IPIV,INFO)
C
      DIMENSION IPIV(*),AX(LDAB,*)
      CALL DGBTRF(M,N,KL,KU,AX,LDAB,IPIV,INFO)
      RETURN
      END
C
      SUBROUTINE LAPACK_DGBTRS(KCH,N,KL,KU,NRHS,AX,LDAB,IPIV,
     &                         X,LDB,INFO)
C
      CHARACTER KCH*(*)
      DIMENSION AX(LDAB,*),X(*)
      CALL DGBTRS(KCH,N,KL,KU,NRHS,AX,LDAB,IPIV,X,LDB,INFO)
      RETURN
      END
C
      SUBROUTINE LAPACK_DGBSV(N,KL,KU,NRHS,AX,LDAB,IPIV,
     &                        X,LDB,INFO)
C
      DIMENSION AX(LDAB,*),X(*)
      CALL DGBSV(N,KL,KU,NRHS,AX,LDAB,IPIV,X,LDB,INFO)
      RETURN
      END
C
      SUBROUTINE LAPACK_DSBTRD(KCH1,KCH2,MMMAX,MHMAX,
     &                         FM,MMM,DM,EM,QM,I,WORK,INFO1)
C
      CHARACTER KCH1*(*),KCH2*(*)
      DIMENSION FM(MMM,*),DM(*),EM(*),QM(1,*)
      CALL DSBTRD(KCH1,KCH2,MMMAX,MHMAX,
     &            FM,MMM,DM,EM,QM,I,WORK,INFO1)
      RETURN
      END
C
      SUBROUTINE LAPACK_DSTEBZ(KCH1,KCH2,MMMAX,VL,VU,IMIN,IMAX,EPSEG,
     &                         DM,EM,INMAX,NSPLIT,W,
     &               IBLOCK,ISPLIT,WORK,IWORK,INFO)
C
      CHARACTER KCH1*(*),KCH2*(*)
      DIMENSION DM(*),EM(*),IBLOCK(*),ISPLIT(*)
      DIMENSION WORK(*),IWORK(*),W(*)
C
      CALL DSTEBZ(KCH1,KCH2,MMMAX,VL,VU,IMIN,IMAX,EPSEG,
     &            DM,EM,INMAX,NSPLIT,W,
     &            IBLOCK,ISPLIT,WORK,IWORK,INFO)
      RETURN
      END
C
      SUBROUTINE LAPACK_DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C
      DIMENSION A(LDA,*),B(LDB,*)
      CALL DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      RETURN
      END
