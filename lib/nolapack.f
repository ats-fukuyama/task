C     $Id$
      SUBROUTINE LAPACK_DGBTRF(M,N,KL,KU,AX,LDAB,IPIV,INFO)
C
      DIMENSION IPIV(*),AX(LDAB,*)
      INFO=-1
      RETURN
      END
C
      SUBROUTINE LAPACK_DGBTRS(KCH,N,KL,KU,NRHS,AX,LDAB,IPIV,
     &                         X,LDB,INFO)
C
      CHARACTER KCH*(*)
      DIMENSION AX(LDAB,*),X(*)
      INFO=-1
      RETURN
      END
C
      SUBROUTINE LAPACK_DGBSV(N,KL,KU,NRHS,AX,LDAB,IPIV,
     &                        X,LDB,INFO)
C
      DIMENSION AX(LDAB,*),X(*)
      INFO=-1
      RETURN
      END
C
      SUBROUTINE LAPACK_DSBTRD(KCH1,KCH2,MMMAX,MHMAX,
     &                         FM,MMM,DM,EM,QM,I,WORK,INFO)
C
      CHARACTER KCH1*(*),KCH2*(*)
      DIMENSION FM(MMM,*),DM(*),EM(*),QM(1,*)
      INFO=-1
      RETURN
      END
C
      SUBROUTINE LAPACK_DSTEBZ(KCH1,KCH2,MMMAX,VL,VU,IMIN,IMAX,EPSEG,
     &                         DM,EM,INMAX,NSPLIT,W,
     &                         IBLOCK,ISPLIT,WORK,IWORK,INFO)
C
      CHARACTER KCH1*(*),KCH2*(*)
      DIMENSION DM(*),EM(*),IBLOCK(*),ISPLIT(*)
      DIMENSION WORK(*),IWORK(*),W(*)
C
      INFO=-1
      RETURN
      END
