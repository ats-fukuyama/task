C     $Id$
C
      SUBROUTINE MPINIT(nprocs1,myrank1)
      nprocs1=1
      myrank1=0
      RETURN
      END
C
      SUBROUTINE MPTERM
      RETURN
      END
C
      SUBROUTINE MPSYNC
      RETURN
      END
C
      SUBROUTINE MPSETI(NMAX,NRANK,ista,iend)
      ista = 1
      iend = NMAX
      RETURN
      END
C
      SUBROUTINE MPBCDN(dtmp,NDTMP)
      REAL*8 dtmp(NDTMP)
      RETURN
      END
C
      SUBROUTINE MPBCRN(rtmp,NRTMP)
      DIMENSION rtmp(NRTMP)
      RETURN
      END
C
      SUBROUTINE MPBCIN(itmp,NITMP)
      DIMENSION itmp(NITMP)
      RETURN
      END
C
      SUBROUTINE MPBCKN(ktmp,NKTMP)
      CHARACTER ktmp*(*)
      RETURN
      END
C
      SUBROUTINE MPBCCN(ctmp,NDTMP)
      COMPLEX*16 ctmp(NDTMP)
      RETURN
      END
C
      SUBROUTINE MPBCDA(D)
      REAL*8 D
      RETURN
      END
C
      SUBROUTINE MPBCRA(R)
      RETURN
      END
C
      SUBROUTINE MPBCIA(I)
      RETURN
      END
C
      SUBROUTINE MPBCKA(K)
      CHARACTER K*1
      RETURN
      END
C
      SUBROUTINE MPBCCA(C)
      COMPLEX*16 C
      RETURN
      END
C
      SUBROUTINE MPGTDN(dtmp,NDTMP)
      REAL*8 dtmp(NDTMP)
      RETURN
      END
C
      SUBROUTINE MPGTRN(rtmp,NRTMP)
      DIMENSION rtmp(NRTMP)
      RETURN
      END
C
      SUBROUTINE MPGTRV(GX,NV,GXTOT,NVTOT,NM)
      DIMENSION GX(NV),GXTOT(NM)
      DO N=1,NV
         GXTOT(N)=GX(N)
      ENDDO
      NVTOT=NV
      RETURN
      END
C
      SUBROUTINE MPGTCV(CX,NV,CXTOT,NVTOT,NM)
      complex*16 CX(NV),CXTOT(NM)
      DO N=1,NV
         CXTOT(N)=CX(N)
      ENDDO
      NVTOT=NV
      RETURN
      END
