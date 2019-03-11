MODULE libbes

  PRIVATE
  PUBLIC BESJN,BESYN,BESIN,BESKN,BESEIN,BESEKN, &
         BESJNV,BESYNV,BESINV,BESKNV,BESEINV,BESEKNV, &
         BESSJN,LAMBDA

CONTAINS

!     ****** LIBRARY OF BESSEL FUNCTIONS *****

  DOUBLE PRECISION FUNCTION BESJN(N,X)
    REAL*8 X,XX,BJ(101),ALPHA,BESJ0,BESJ1
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESJN=BESJ0(XX)
    ELSEIF(NN.EQ.1) THEN
       BESJN=BESJ1(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DJBESL(XX,ALPHA,NN+1,BJ,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESJN: NCALC=',NCALC
       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESJN: NCALC=',NCALC
       BESJN=BJ(NN+1)
    ELSE
       WRITE(6,*) 'XX BESJN: ABS(N).GT.100: N=',N
       BESJN=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESJN=-BESJN
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESJN=-BESJN
    RETURN
  END FUNCTION BESJN

  DOUBLE PRECISION FUNCTION BESYN(N,X)
    REAL*8 X,XX,BY(101),ALPHA,BESY0,BESY1
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESYN=BESY0(XX)
    ELSEIF(NN.EQ.1) THEN
       BESYN=BESY1(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DYBESL(XX,ALPHA,NN+1,BY,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESYN: NCALC=',NCALC
       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESYN: NCALC=',NCALC
       BESYN=BY(NN+1)
    ELSE
       WRITE(6,*) 'XX BESYN: ABS(N).GT.100: N=',N
       BESYN=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESYN=-BESYN
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESYN=-BESYN
    RETURN
  END FUNCTION BESYN

  DOUBLE PRECISION FUNCTION BESIN(N,X)
    REAL*8 X,XX,BI(101),ALPHA,BESI0,BESI1
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESIN=BESI0(XX)
    ELSEIF(NN.EQ.1) THEN
       BESIN=BESI1(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DIBESL(XX,ALPHA,NN+1,1,BI,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESIN: NCALC=',NCALC
!         IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESIN: NCALC=',NCALC
       BESIN=BI(NN+1)
    ELSE
       WRITE(6,*) 'XX BESIN: ABS(N).GT.100: N=',N
       BESIN=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESIN=-BESIN
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESIN=-BESIN
    RETURN
  END FUNCTION BESIN

  DOUBLE PRECISION FUNCTION BESKN(N,X)
    REAL*8 X,XX,BK(101),ALPHA,BESK0,BESK1
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESKN=BESK0(XX)
    ELSEIF(NN.EQ.1) THEN
       BESKN=BESK1(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DKBESL(XX,ALPHA,NN+1,1,BK,NCALC)
       IF(NCALC.LT.0) WRITE(6,*) 'XX BESKN: NCALC=',NCALC
!         IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESKN: NCALC=',NCALC
       BESKN=BK(NN+1)
    ELSE
       WRITE(6,*) 'XX BESKN: ABS(N).GT.100: N=',N
       BESKN=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESKN=-BESKN
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESKN=-BESKN
    RETURN
  END FUNCTION BESKN

  DOUBLE PRECISION FUNCTION BESEIN(N,X)
    REAL*8 X,XX,BI(101),ALPHA,BESEI0,BESEI1
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESEIN=BESEI0(XX)
    ELSEIF(NN.EQ.1) THEN
       BESEIN=BESEI1(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DIBESL(XX,ALPHA,NN+1,2,BI,NCALC)
       IF(NCALC.LT.NN+1) WRITE(6,*) 'XX BESEIN: NCALC=',NCALC
       BESEIN=BI(NN+1)
    ELSE
       WRITE(6,*) 'XX BESEIN: ABS(N).GT.100: N=',N
       BESEIN=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESEIN=-BESEIN
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESEIN=-BESEIN
    RETURN
  END FUNCTION BESEIN

  DOUBLE PRECISION FUNCTION BESEKN(N,X)
    REAL*8 X,XX,BK(101),ALPHA,BESEK0,BESEK1
    INTEGER:: N, NN, NCALC
    DATA ALPHA/0.D0/

    NN=ABS(N)
    XX=ABS(X)
    IF(NN.EQ.0) THEN
       BESEKN=BESEK0(XX)
    ELSEIF(NN.EQ.1) THEN
       BESEKN=BESEK1(XX)
    ELSEIF(NN.LE.100) THEN
       CALL DKBESL(XX,ALPHA,NN+1,2,BK,NCALC)
       IF(NCALC.LT.NN+1) THEN
          WRITE(6,*) 'XX BESEKN: NCALC=',NCALC
          WRITE(6,'(A,I5,1P2E12.4)') 'XX NN,XX,ALPHA=',NN,XX,ALPHA
          STOP
       END IF
       BESEKN=BK(NN+1)
    ELSE
       WRITE(6,*) 'XX BESEKN: ABS(N).GT.100: N=',N
       BESEKN=0.D0
    ENDIF
    IF(N.LT.0.AND.MOD(NN,2).EQ.1) BESEKN=-BESEKN
    IF(X.LT.0.D0.AND.MOD(NN,2).EQ.1) BESEKN=-BESEKN
    RETURN
  END FUNCTION BESEKN

  SUBROUTINE BESJNV(N,X,V,IERR)
    REAL*8 X,XX,V(0:N),ALPHA
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
!         WRITE(6,*) 'XX BESJNV: NCALC,N=',NCALC,N
       IERR=10+NCALC
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
    REAL*8 X,XX,V(0:N),ALPHA
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
    REAL*8 X,XX,V(0:N),ALPHA
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
    REAL*8 X,XX,V(0:N),ALPHA
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
    REAL*8 X,XX,V(0:N),ALPHA
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
    REAL*8 X,XX,V(0:N),ALPHA
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
        REAL(kind=8),INTENT(IN):: X
        INTEGER,INTENT(IN):: NMAX
        REAL(kind=8),INTENT(OUT):: BJN(0:NMAX),DBJN(0:NMAX)
        INTEGER:: n,ierr

      CALL BESJNV(NMAX,X,BJN,IERR)
!      write(6,*) ierr,X
!      WRITE(6,'(1P5E12.4)') (BJN(n),n=0,nmax)
      IF(IERR.NE.0) RETURN

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
!      WRITE(6,'(1P5E12.4)') (DBJN(n),n=0,nmax)
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
