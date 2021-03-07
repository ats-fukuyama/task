

MODULE libbes

  USE task_kinds,ONLY: dp
  PRIVATE
  PUBLIC BESJNX,BESYNX,BESINX,BESKNX,BESEINX,BESEKNX, &
         BESJNV,BESYNV,BESINV,BESKNV,BESEINV,BESEKNV, &
         BESSJN,LAMBDA

CONTAINS

!     ****** LIBRARY OF BESSEL FUNCTIONS *****

  FUNCTION BESJNX(N,X)
    REAL(dp):: BESJNX
    REAL(dp):: X,XX,BJ(101),ALPHA,BESJ0X,BESJ1X
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESJNX=BESJ0X(XX)
    ELSEIF(NN.EQ.1) THEN
       BESJNX=BESJ1X(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DJBESL(XX,ALPHA,NN+1,BJ,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESJNX: NCALC=',NCALC
!       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESJNX: NCALC=',NCALC
       BESJNX=BJ(NN+1)
    ELSE
       WRITE(6,*) 'XX BESJNX: ABS(N).GT.100: N=',N
       BESJNX=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESJNX=-BESJNX
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESJNX=-BESJNX
    RETURN
  END FUNCTION BESJNX

  FUNCTION BESYNX(N,X)
    REAL(dp):: BESYNX
    REAL(dp):: X,XX,BY(101),ALPHA,BESY0X,BESY1X
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESYNX=BESY0X(XX)
    ELSEIF(NN.EQ.1) THEN
       BESYNX=BESY1X(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DYBESL(XX,ALPHA,NN+1,BY,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESYNX: NCALC=',NCALC
!       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESYNX: NCALC=',NCALC
       BESYNX=BY(NN+1)
    ELSE
       WRITE(6,*) 'XX BESYNX: ABS(N).GT.100: N=',N
       BESYNX=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESYNX=-BESYNX
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESYNX=-BESYNX
    RETURN
  END FUNCTION BESYNX

  DOUBLE PRECISION FUNCTION BESINX(N,X)
    REAL*8 X,XX,BI(101),ALPHA,BESI0X,BESI1X
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESINX=BESI0X(XX)
    ELSEIF(NN.EQ.1) THEN
       BESINX=BESI1X(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DIBESL(XX,ALPHA,NN+1,1,BI,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESINX: NCALC=',NCALC
!         IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESIN: NCALC=',NCALC
       BESINX=BI(NN+1)
    ELSE
       WRITE(6,*) 'XX BESINX: ABS(N).GT.100: N=',N
       BESINX=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESINX=-BESINX
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESINX=-BESINX
    RETURN
  END FUNCTION BESINX

  DOUBLE PRECISION FUNCTION BESKNX(N,X)
    REAL*8 X,XX,BK(101),ALPHA,BESK0X,BESK1X
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESKNX=BESK0X(XX)
    ELSEIF(NN.EQ.1) THEN
       BESKNX=BESK1X(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DKBESL(XX,ALPHA,NN+1,1,BK,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESKNX: NCALC=',NCALC
!         IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESKNX: NCALC=',NCALC
       BESKNX=BK(NN+1)
    ELSE
       WRITE(6,*) 'XX BESKNX: ABS(N).GT.100: N=',N
       BESKNX=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESKNX=-BESKNX
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESKNX=-BESKNX
    RETURN
  END FUNCTION BESKNX

  DOUBLE PRECISION FUNCTION BESEINX(N,X)
    REAL*8 X,XX,BI(101),ALPHA,BESEI0X,BESEI1X
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESEINX=BESEI0X(XX)
    ELSEIF(NN.EQ.1) THEN
       BESEINX=BESEI1X(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DIBESL(XX,ALPHA,NN+1,2,BI,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESEINX: NCALC=',NCALC
!       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESEINX: NCALC=',NCALC
       BESEINX=BI(NN+1)
    ELSE
       WRITE(6,*) 'XX BESEIN: ABS(N).GT.100: N=',N
       BESEINX=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESEINX=-BESEINX
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESEINX=-BESEINX
    RETURN
  END FUNCTION BESEINX

  DOUBLE PRECISION FUNCTION BESEKNX(N,X)
    REAL*8 X,XX,BK(101),ALPHA,BESEK0X,BESEK1X
    INTEGER:: N, NN, NCALC
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESEKNX=BESEK0X(XX)
    ELSEIF(NN.EQ.1) THEN
       BESEKNX=BESEK1X(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DKBESL(XX,ALPHA,NN+1,2,BK,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESEKNX: NCALC=',NCALC
!       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESEKNX: NCALC=',NCALC
       BESEKNX=BK(NN+1)
    ELSE
       WRITE(6,*) 'XX BESEKNX: ABS(N).GT.100: N=',N
       BESEKNX=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESEKNX=-BESEKNX
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESEKNX=-BESEKNX
    RETURN
  END FUNCTION BESEKNX

  SUBROUTINE BESJNV(N,X,V,IERR)
    REAL(dp):: X,XX,V(0:N),ALPHA
    DATA ALPHA/0.D0/
    
    IF(N.LT.0) THEN
       IERR=1
       RETURN
    ENDIF
    XX=ABS(X)
         
    CALL DJBESL(XX,ALPHA,N+1,V,NCALC)
    IF(NCALC.LT.0) THEN
       WRITE(6,*) 'XX BESJNV: NCALC=',NCALC
       WRITE(6,*) 'XX BESJNV: N,X=',N,X
       IERR=100+ABS(NCALC)
    ELSEIF(NCALC.LT.N+1) THEN
       DO nn=ncalc,N
          V(nn)=0.D0
       END DO
       IERR=0
    ELSE
       IERR=0
    ENDIF
    IF(X.LT.0.D0) THEN
       DO I=1,N
          IF(MOD(I,2).EQ.1) V(I)=-V(I)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE BESJNV

  SUBROUTINE BESYNV(N,X,V,IERR)
    REAL(dp):: X,XX,V(0:N),ALPHA
    DATA ALPHA/0.D0/

    IF(N.LT.0) THEN
       IERR=1
       RETURN
    ENDIF
    XX=ABS(X)
         
    CALL DYBESL(XX,ALPHA,N+1,V,NCALC)
    IF(NCALC.LT.0) THEN
       WRITE(6,*) 'XX BESYNV: NCALC=',NCALC
       WRITE(6,*) 'XX BESYNV: N,X=',N,X
       IERR=100+ABS(NCALC)
    ELSEIF(NCALC.LT.N+1) THEN
!         WRITE(6,*) 'XX BESYNV: NCALC,N=',NCALC,N
       IERR=10+NCALC
    ELSE
       IERR=0
    ENDIF
    IF(X.LT.0.D0) THEN
       DO I=1,N+1
          IF(MOD(I,2).EQ.1) V(I)=-V(I)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE BESYNV

  SUBROUTINE BESINV(N,X,V,IERR)
    REAL(dp):: X,XX,V(0:N),ALPHA
    DATA ALPHA/0.D0/

    IF(N.LT.0) THEN
       IERR=1
       RETURN
    ENDIF
    XX=ABS(X)
         
    CALL DIBESL(XX,ALPHA,N+1,1,V,NCALC)
    IF(NCALC.LT.0) THEN
       WRITE(6,*) 'XX BESINV: NCALC=',NCALC
       WRITE(6,*) 'XX BESINV: N,X=',N,X
       IERR=100+ABS(NCALC)
    ELSEIF(NCALC.LT.N+1) THEN
!         WRITE(6,*) 'XX BESINV: NCALC,N=',NCALC,N
       IERR=10+NCALC
    ELSE
       IERR=0
    ENDIF
    IF(X.LT.0.D0) THEN
       DO I=1,N+1
          IF(MOD(I,2).EQ.1) V(I)=-V(I)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE BESINV

  SUBROUTINE BESKNV(N,X,V,IERR)
    REAL(dp):: X,XX,V(0:N),ALPHA
    DATA ALPHA/0.D0/

    IF(N.LT.0) THEN
       IERR=1
       RETURN
    ENDIF
    XX=ABS(X)
         
    CALL DKBESL(XX,ALPHA,N+1,1,V,NCALC)
    IF(NCALC.LT.0) THEN
       WRITE(6,*) 'XX BESKNV: NCALC=',NCALC
       WRITE(6,*) 'XX BESKNV: N,X=',N,X
       IERR=100+ABS(NCALC)
    ELSEIF(NCALC.LT.N+1) THEN
!         WRITE(6,*) 'XX BESKNV: NCALC,N=',NCALC,N
       IERR=10+NCALC
    ELSE
       IERR=0
    ENDIF
    IF(X.LT.0.D0) THEN
       DO I=1,N+1
          IF(MOD(I,2).EQ.1) V(I)=-V(I)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE BESKNV

  SUBROUTINE BESEINV(N,X,V,IERR)
    REAL(dp):: X,XX,V(0:N),ALPHA
    DATA ALPHA/0.D0/

    IF(N.LT.0) THEN
       IERR=1
       RETURN
    ENDIF
    XX=ABS(X)
         
    CALL DIBESL(XX,ALPHA,N+1,2,V,NCALC)
    IF(NCALC.LT.0) THEN
       WRITE(6,*) 'XX BESEINV: NCALC=',NCALC
       WRITE(6,*) 'XX BESEINV: N,X=',N,X
       IERR=100+ABS(NCALC)
    ELSEIF(NCALC.LT.N+1) THEN
!         WRITE(6,*) 'XX BESEINV: NCALC,N=',NCALC,N
       IERR=10+NCALC
    ELSE
       IERR=0
    ENDIF
    IF(X.LT.0.D0) THEN
       DO I=1,N+1
          IF(MOD(I,2).EQ.1) V(I)=-V(I)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE BESEINV

  SUBROUTINE BESEKNV(N,X,V,IERR)
    REAL(dp):: X,XX,V(0:N),ALPHA
    DATA ALPHA/0.D0/

    IF(N.LT.0) THEN
       IERR=1
       RETURN
    ENDIF
    XX=ABS(X)
         
    CALL DKBESL(XX,ALPHA,N+1,2,V,NCALC)
    IF(NCALC.LT.0) THEN
       WRITE(6,*) 'XX BESEKNV: NCALC=',NCALC
       WRITE(6,*) 'XX BESEKNV: N,X=',N,X
       IERR=100+ABS(NCALC)
    ELSEIF(NCALC.LT.N+1) THEN
!         WRITE(6,*) 'XX BESEKNV: NCALC,N=',NCALC,N
       IERR=10+NCALC
    ELSE
       IERR=0
    ENDIF
    IF(X.LT.0.D0) THEN
       DO I=1,N+1
          IF(MOD(I,2).EQ.1) V(I)=-V(I)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE BESEKNV


!     **************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND
!
!     **************************************

      SUBROUTINE BESSJN(X,NMAX,BJN,DBJN)

        IMPLICIT NONE
        REAL(dp),INTENT(IN):: X
        INTEGER,INTENT(IN):: NMAX
        REAL(dp),INTENT(OUT):: BJN(0:NMAX),DBJN(0:NMAX)
        INTEGER:: n,ierr

      CALL BESJNV(NMAX,X,BJN,IERR)
!      write(6,*) ierr,X
!      WRITE(6,'(5ES12.4)') (BJN(n),n=0,nmax)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX BESSJN: BESJNV error: IERR=',IERR
         STOP
      END IF

      DBJN(0)=BJN(1)
      IF(ABS(X).EQ.0.D0)THEN
         DBJN(1)=0.5D0
         DO N=2,NMAX
            DBJN(N)=0.D0
         ENDDO
      ELSE
         DO N=1,NMAX
            DBJN(N)=BJN(N-1)-N*BJN(N)/X
         ENDDO
      ENDIF
!      WRITE(6,'(5ES12.4)') (DBJN(n),n=0,nmax)
      RETURN
    END SUBROUTINE BESSJN

!     ******************************************************************
!
!                       ****** LAMBDA FUNCTION ******
!
!        MODIFIED BESSEL FUNCTION OF THE FIRST KIND, MULIPLIED BY EXP
!                             COMPLEX ARGUMENT
!
!     ******************************************************************

      SUBROUTINE LAMBDA(N,CX,CALAM,IERR)

      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)

      DIMENSION CALAM(0:N)
      DATA ONE/1.D0/
      DATA D55/1.D-55/

      IERR=0
      XA=ABS(CX)
      NA=ABS(N)

      IF(NA.GE.30000.OR.N.LT.0) THEN
         CALAM(0)=(0.D0,0.D0)
         IERR=1
         RETURN
      ENDIF
      IF(XA.GE.173.D0) THEN
         CALAM(0)=(0.D0,0.D0)
         IERR=2
         RETURN
      ENDIF
      IF(XA.LE.1.D-8) THEN
         IF(XA.LE.1.D-77) THEN
            CALAM(0)=1.D0
            DO I=1,NA
               CALAM(I)=0.D0
            ENDDO
         ELSE
            CALAM(0)=EXP(-CX)
            CT1=0.5D0*CX
            T2=ONE
            CT3=ONE
            DO I=1,NA
               IF(ABS(CT3).LE.1.D-77*ABS(T2/CT1)) THEN
                  CALAM(I)=0.D0
               ELSE
                  CT3=CT3*CT1/T2
                  T2=T2+ONE
                  CALAM(I)=CT3*EXP(-CX)
               ENDIF
            ENDDO
         ENDIF
      ELSE
         CZ=2.D0/CX
         IF(XA.GE.10.D0) THEN
            L=40
         ELSEIF(XA.GE.0.1D0) THEN
            L=25
         ELSE
            L=10
         ENDIF
         RX=XA
         NM=MAX(NA,INT(RX))+L
         CT3=0.D0
         CT2=1.D-75
         CS=0.D0
         DO I=1,NM
            K=NM-I
            CT1=(K+1)*CT2*CZ+CT3
            IF(K.LE.NA) CALAM(K)=CT1
            CS=CS+CT1
            IF(ABS(CS).GT.1.D55) THEN
               CT1=CT1*D55
               CT2=CT2*D55
               CS=CS*D55
               DO J=K,NA
                  CALAM(J)=CALAM(J)*D55
               ENDDO
            ENDIF
            CT3=CT2
            CT2=CT1
         ENDDO
         CS=CS+CS-CT1
         DO J=0,NA
            CALAM(J)=CALAM(J)/CS
         ENDDO
      ENDIF
      RETURN
    END SUBROUTINE LAMBDA
END MODULE libbes
