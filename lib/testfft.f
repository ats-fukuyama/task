C     $Id$
C
      IMPLICIT NONE
      INTEGER NM
      PARAMETER (NM=1024)
C
      COMPLEX*16 A(NM),B(NM),C(NM),AS(NM),BS(NM)
      REAL*8 WORK(NM)
C      COMPLEX*16 WORK(NM)
      INTEGER IWORK(NM)
      INTEGER NMAX,NSAVE,N,IND,LP
C
      NSAVE=0
C
    1 WRITE(6,*) '## NMAX=?'
      READ(5,*,END=9000,ERR=1) NMAX
C
      DO N=1,NMAX
         A(N)=0.D0
      ENDDO
      A(1)=(1.D0,1.D0)
      DO N=2,NMAX
         A(N)=DCMPLX(1.D0/N,0.D0)
      ENDDO
C
      IF(NMAX.NE.NSAVE) IND=1
      NSAVE=NMAX
      LP=NINT(LOG(DBLE(N))/LOG(2.D0))
C
      DO N=1,NMAX
         AS(N)=A(N)
      ENDDO
C      CALL FFT2L(A,B,WORK,IWORK,NMAX/2,IND,0+1,LP)
      CALL FFT2L(A,B,WORK,IWORK,NMAX,IND,0)
C
      DO N=1,NMAX
         BS(N)=B(N)
      ENDDO
C      CALL FFT2L(B,C,WORK,IWORK,NMAX/2,IND,1+1,LP)
      CALL FFT2L(B,C,WORK,IWORK,NMAX,IND,1)
C
      WRITE(6,'(I5,1P6E12.4)') (N,AS(N),BS(N),C(N),N=1,NMAX)
      GOTO 1
C
 9000 STOP
      END
      
C
