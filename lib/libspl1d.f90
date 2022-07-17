! libspl1d.f90

MODULE libspl1d

  PRIVATE
  PUBLIC spl1d      ! calculate real spline coefficients
  PUBLIC spl1df     ! interpolate real function by spline 
  PUBLIC spl1dd     ! interpolate real function and derivative by spline 
  PUBLIC spl1ddd    ! interpolate real function and derivatives by spline 
  PUBLIC spl1di0    ! calculate real integral for spline
  PUBLIC spl1di     ! integrate real function by spline
  PUBLIC cspl1d     ! calculate complex spline coefficients
  PUBLIC cspl1df    ! interpolate complex function by spline 
  PUBLIC cspl1dd    ! interpolate complex function and derivataive by spline 
  PUBLIC cspl1di0   ! calculate complex integral for spline
  PUBLIC cspl1di    ! integrate compled function by spline

  PUBLIC tdmsrdx,tdmprdx,tdmscdx,tdmpcdx ! tri-diagonal matrix solver
  
CONTAINS

!   ************************************************
!   **           Spline Interpolation             **
!   ************************************************

!     ****** One-Dimensional Spline Interpolation ******
!          **** Calculation of coefficients ****

      SUBROUTINE SPL1D(X,F,FX,U,NXMAX,ID,IERR)

!      INPUT : X(NXMAX)  : COORDINATES
!              F(NXMAX)  : VALUE
!              FX(NXMAX) : EDGE DERIVATIVE FOR 1<= ID <=3
!              NXMAX     : NUMBER OF VARIABLES
!              ID        : 0 : SECOND DERIVATIVES = 0 AT X(1) AND X(NXMAX)
!                          1 : DERIVATIVE FX(1) IS GIVEN
!                          2 : DERIVATIVE FX(NXMAX) IS GIVEN
!                          3 : DERIVATIVES FX(1) AND FX(NXMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NXMAX)
!      OUTPUT: U(4,NXMAX): SPLINE COEFICIENTS
!              FX(NXMAX) : ESTIMATED DERIVATIVES
!              IERR      : ERROR INDICATOR

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                  INTENT(IN) :: NXMAX, ID
      REAL(dp), DIMENSION(NXMAX),   INTENT(IN) :: X, F
      REAL(dp), DIMENSION(NXMAX), INTENT(INOUT):: FX
      INTEGER,                  INTENT(OUT):: IERR
      REAL(dp), DIMENSION(4,NXMAX), INTENT(OUT):: U
      INTEGER  :: ID1, ID2, NX
      REAL(dp)     :: DX, DX1, DX2, DX3, DXM, DXP, T11, T21, T31, T41, &
     &               V31, V32, V33, V34, V41, V42, V43, V44

      IERR=0

      ID1=MOD(ID,2)
      ID2=MOD(ID/2,2)

      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         U(1,1)=DXP
         U(2,1)=2.D0*(DXP+DXM)
         U(3,1)=DXM
      ELSE
         IF(ID1.EQ.0) THEN
            DXP=X(2)-X(1)
            U(1,1)=0.D0
            U(2,1)=2.D0*DXP
            U(3,1)=DXP
         ELSE
            U(1,1)=0.D0
            U(2,1)=1.D0
            U(3,1)=0.D0
         ENDIF
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         U(1,NX)=DXP
         U(2,NX)=2.D0*(DXP+DXM)
         U(3,NX)=DXM
      ENDDO
      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         U(1,NXMAX)=DXP
         U(2,NXMAX)=2.D0*(DXP+DXM)
         U(3,NXMAX)=DXM
      ELSE
         IF(ID2.EQ.0) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            U(1,NXMAX)=DXM
            U(2,NXMAX)=2.D0*DXM
            U(3,NXMAX)=0.D0
         ELSE
            U(1,NXMAX)=0.D0
            U(2,NXMAX)=1.D0
            U(3,NXMAX)=0.D0
         ENDIF
      ENDIF

      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         FX(1)=3.D0*(DXM*(F(2)-F(1))/DXP &
     &              +DXP*(F(NXMAX)-F(NXMAX-1))/DXM)
      ELSE
         IF(ID1.EQ.0) THEN
            FX(1)=3.D0*(F(2)-F(1))
         ENDIF
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         FX(NX)=3.D0*(DXM*(F(NX+1)-F(NX))/DXP &
     &               +DXP*(F(NX)-F(NX-1))/DXM)
      ENDDO
      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         FX(NXMAX)=3.D0*(DXM*(F(2)-F(1))/DXP &
     &                  +DXP*(F(NXMAX)-F(NXMAX-1))/DXM)
      ELSE
         IF(ID2.EQ.0) THEN
            FX(NXMAX)=3.D0*(F(NXMAX)-F(NXMAX-1))
         ENDIF
      ENDIF

      IF(ID.EQ.4) THEN
         CALL TDMPRDX(U,FX,NXMAX-1,IERR)
         IF(IERR.NE.0) GOTO 9001
         FX(NXMAX)=FX(1)
      ELSE
         CALL TDMSRDX(U,FX,NXMAX,IERR)
         IF(IERR.NE.0) GOTO 9002
      ENDIF

      NX=1
      U(1,NX)=0.D0
      U(2,NX)=0.D0
      U(3,NX)=0.D0
      U(4,NX)=0.D0
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)

         T11=  F(NX-1)
         T21= FX(NX-1)
         T31=  F(NX  )
         T41= FX(NX  )

         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         V31=-3.D0*DX2
         V32=-2.D0*DX1
         V33= 3.D0*DX2
         V34=-1.D0*DX1
         V41= 2.D0*DX3
         V42=      DX2
         V43=-2.D0*DX3
         V44=      DX2

         U(1,NX)=    T11
         U(2,NX)=            T21
         U(3,NX)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,NX)=V41*T11+V42*T21+V43*T31+V44*T41
      ENDDO
      RETURN

 9001 WRITE(6,*) 'XX SPL1D: TDMPRDX ERROR : IERR=',IERR
      IERR=1
      RETURN
 9002 WRITE(6,*) 'XX SPL1D: TDMSRDX ERROR : IERR=',IERR
      IERR=2
      RETURN
      END SUBROUTINE SPL1D

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL1DF(X0,F0,X,U,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                INTENT(IN) :: NXMAX
      REAL(dp),                   INTENT(IN) :: X0
      REAL(dp),DIMENSION(NXMAX),  INTENT(IN) :: X
      REAL(dp),DIMENSION(4,NXMAX),INTENT(IN) :: U
      REAL(dp),                   INTENT(OUT):: F0
      INTEGER,                INTENT(OUT):: IERR
      REAL(dp)     :: FS, DX
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      F0= U(1,NX) &
     &  + U(2,NX)*DX &
     &  + U(3,NX)*DX*DX &
     &  + U(4,NX)*DX*DX*DX
!      WRITE(6,'(A,I8,3ES12.4)') &
!     &     'NX,X0,X(NX-1),X(NX)=',NX,X0,X(NX-1),X(NX)
      RETURN
      END SUBROUTINE SPL1DF

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL1DD(X0,F0,DF0,X,U,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                INTENT(IN) :: NXMAX
      REAL(dp),                   INTENT(IN) :: X0
      REAL(dp),DIMENSION(NXMAX),  INTENT(IN) :: X
      REAL(dp),DIMENSION(4,NXMAX),INTENT(IN) :: U
      REAL(dp),                   INTENT(OUT):: F0, DF0
      INTEGER,                INTENT(OUT):: IERR
      REAL(dp)     :: FS, DX
      INTEGER  :: NX

      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      F0= U(1,NX) &
     &  + U(2,NX)*DX &
     &  + U(3,NX)*DX*DX &
     &  + U(4,NX)*DX*DX*DX
      DF0= U(2,NX) &
     &   + U(3,NX)*DX*2 &
     &   + U(4,NX)*DX*DX*3
      IERR=0
      RETURN
      END SUBROUTINE SPL1DD

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL1DDD(X0,F0,DF0,DDF0,X,U,NXMAX,IERR)

      USE task_kinds,ONLY: dp  
      IMPLICIT NONE
      INTEGER,                    INTENT(IN) :: NXMAX
      REAL(dp),                   INTENT(IN) :: X0
      REAL(dp),DIMENSION(NXMAX),  INTENT(IN) :: X
      REAL(dp),DIMENSION(4,NXMAX),INTENT(IN) :: U
      REAL(dp),                   INTENT(OUT):: F0, DF0, DDF0
      INTEGER,                    INTENT(OUT):: IERR
      REAL(dp) :: FS, DX
      INTEGER  :: NX

      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      F0= U(1,NX) &
     &  + U(2,NX)*DX &
     &  + U(3,NX)*DX*DX &
     &  + U(4,NX)*DX*DX*DX
      DF0= U(2,NX) &
     &   + U(3,NX)*DX*2 &
     &   + U(4,NX)*DX*DX*3
      DDF0= U(3,NX)*2 &
     &    + U(4,NX)*DX*6
      IERR=0
      RETURN
      END SUBROUTINE SPL1DDD

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of U0  ****

      SUBROUTINE SPL1DI0(X,U,U0,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                INTENT(IN) :: NXMAX
      REAL(dp),DIMENSION(NXMAX),  INTENT(IN) :: X
      REAL(dp),DIMENSION(4,NXMAX),INTENT(IN) :: U
      REAL(dp),DIMENSION(NXMAX),  INTENT(OUT):: U0
      INTEGER,                INTENT(OUT):: IERR
      REAL(dp)     :: DX
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF

      U0(1)=0.D0
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)

         U0(NX)= U0(NX-1) &
     &         + U(1,NX)*DX &
     &         + U(2,NX)*DX*DX/2.D0 &
     &         + U(3,NX)*DX*DX*DX/3.D0 &
     &         + U(4,NX)*DX*DX*DX*DX/4.D0
      ENDDO

!     WRITE(6,'(A,2I5,3ES12.4)')
!     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END SUBROUTINE SPL1DI0

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL1DI(X0,FI,X,U,U0,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                INTENT(IN) :: NXMAX
      REAL(dp),                   INTENT(IN) :: X0
      REAL(dp),DIMENSION(NXMAX) , INTENT(IN) :: X
      REAL(dp),DIMENSION(4,NXMAX),INTENT(IN) :: U
      REAL(dp),DIMENSION(NXMAX),  INTENT(IN) :: U0
      REAL(dp),                   INTENT(OUT):: FI
      INTEGER,                INTENT(OUT):: IERR
      REAL(dp)     :: DX, FS
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      FI= U0(NX-1) &
     &  + U(1,NX)*DX &
     &  + U(2,NX)*DX*DX/2.D0 &
     &  + U(3,NX)*DX*DX*DX/3.D0 &
     &  + U(4,NX)*DX*DX*DX*DX/4.D0
!      WRITE(6,'(A,2I5,3ES12.4)')
!     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END SUBROUTINE SPL1DI

!     ****** One-Dimensional Spline Interpolation ******
!          **** Calculation of coefficients ****

      SUBROUTINE CSPL1D(X,F,FX,U,NXMAX,ID,IERR)

!      INPUT : X(NXMAX)  : COORDINATES
!              F(NXMAX)  : VALUE
!              FX(NXMAX) : EDGE DERIVATIVE FOR 1<= ID <=3
!              NXMAX     : NUMBER OF VARIABLES
!              ID        : 0 : SECOND DERIVATIVES = 0 AT X(1) AND X(NXMAX)
!                          1 : DERIVATIVE FX(1) IS GIVEN
!                          2 : DERIVATIVE FX(NXMAX) IS GIVEN
!                          3 : DERIVATIVES FX(1) AND FX(NXMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NXMAX)
!      OUTPUT: U(4,NXMAX): SPLINE COEFICIENTS
!              FX(NXMAX) : ESTIMATED DERIVATIVES
!              IERR      : ERROR INDICATOR

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                   INTENT(IN)   :: NXMAX, ID
      REAL(dp),    DIMENSION(NXMAX), INTENT(IN)   :: X
      COMPLEX(dp),DIMENSION(NXMAX),  INTENT(IN)   :: F
      COMPLEX(dp),DIMENSION(NXMAX),  INTENT(INOUT):: FX
      COMPLEX(dp),DIMENSION(4,NXMAX),INTENT(INOUT):: U
      INTEGER,                   INTENT(OUT)  :: IERR
      REAL(dp)     :: DXM, DXP, DX, DX1, DX2, DX3, V31, V32, V33, V34, &
     &               V41, V42, V43, V44
      INTEGER  :: ID1, ID2, NX
      COMPLEX(dp) :: T11, T21, T31, T41

      IERR=0

      ID1=MOD(ID,2)
      ID2=MOD(ID/2,2)

      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         U(1,1)=DXP
         U(2,1)=2.D0*(DXP+DXM)
         U(3,1)=DXM
      ELSE
         IF(ID1.EQ.0) THEN
            DXP=X(2)-X(1)
            U(1,1)=0.D0
            U(2,1)=2.D0*DXP
            U(3,1)=DXP
         ELSE
            U(1,1)=0.D0
            U(2,1)=1.D0
            U(3,1)=0.D0
         ENDIF
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         U(1,NX)=DXP
         U(2,NX)=2.D0*(DXP+DXM)
         U(3,NX)=DXM
      ENDDO
      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         U(1,NXMAX)=DXP
         U(2,NXMAX)=2.D0*(DXP+DXM)
         U(3,NXMAX)=DXM
      ELSE
         IF(ID2.EQ.0) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            U(1,NXMAX)=DXM
            U(2,NXMAX)=2.D0*DXM
            U(3,NXMAX)=0.D0
         ELSE
            U(1,NXMAX)=0.D0
            U(2,NXMAX)=1.D0
            U(3,NXMAX)=0.D0
         ENDIF
      ENDIF

      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         FX(1)=3.D0*(DXM*(F(2)-F(1))/DXP &
     &              +DXP*(F(NXMAX)-F(NXMAX-1))/DXM)
      ELSE
         IF(ID1.EQ.0) THEN
            FX(1)=3.D0*(F(2)-F(1))
         ENDIF
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         FX(NX)=3.D0*(DXM*(F(NX+1)-F(NX))/DXP &
     &               +DXP*(F(NX)-F(NX-1))/DXM)
      ENDDO
      IF(ID.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         FX(NXMAX)=3.D0*(DXM*(F(2)-F(1))/DXP &
     &                  +DXP*(F(NXMAX)-F(NXMAX-1))/DXM)
      ELSE
         IF(ID2.EQ.0) THEN
            FX(NXMAX)=3.D0*(F(NXMAX)-F(NXMAX-1))
         ENDIF
      ENDIF

      IF(ID.EQ.4) THEN
         CALL TDMPCDX(U,FX,NXMAX-1,IERR)
         IF(IERR.NE.0) GOTO 9001
         FX(NXMAX)=FX(1)
      ELSE
         CALL TDMSCDX(U,FX,NXMAX,IERR)
         IF(IERR.NE.0) GOTO 9002
      ENDIF

      NX=1
      U(1,NX)=0.D0
      U(2,NX)=0.D0
      U(3,NX)=0.D0
      U(4,NX)=0.D0
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)

         T11=  F(NX-1)
         T21= FX(NX-1)
         T31=  F(NX  )
         T41= FX(NX  )

         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         V31=-3.D0*DX2
         V32=-2.D0*DX1
         V33= 3.D0*DX2
         V34=-1.D0*DX1
         V41= 2.D0*DX3
         V42=      DX2
         V43=-2.D0*DX3
         V44=      DX2

         U(1,NX)=    T11
         U(2,NX)=            T21
         U(3,NX)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,NX)=V41*T11+V42*T21+V43*T31+V44*T41
      ENDDO
      RETURN

 9001 WRITE(6,*) 'XX CSPL1D: TDMPRDX ERROR : IERR=',IERR
      IERR=1
      RETURN
 9002 WRITE(6,*) 'XX CSPL1D: TDMSRDX ERROR : IERR=',IERR
      IERR=2
      RETURN
      END SUBROUTINE CSPL1D

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL1DF(X0,F0,X,U,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: NXMAX
      REAL(dp),                  INTENT(IN) :: X0
      REAL(dp),DIMENSION(NXMAX), INTENT(IN) :: X
      COMPLEX(dp),DIMENSION(4,NXMAX),INTENT(IN):: U
      COMPLEX(dp),               INTENT(OUT):: F0
      INTEGER,               INTENT(OUT):: IERR
      REAL(dp)     :: FS, DX
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      F0= U(1,NX) &
     &  + U(2,NX)*DX &
     &  + U(3,NX)*DX*DX &
     &  + U(4,NX)*DX*DX*DX
!      WRITE(6,'(A,2I5,3ES12.4)')
!     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END SUBROUTINE CSPL1DF

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL1DD(X0,F0,DF0,X,U,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: NXMAX
      REAL(dp),                  INTENT(IN) :: X0
      REAL(dp),DIMENSION(NXMAX), INTENT(IN) :: X
      COMPLEX(dp),DIMENSION(4,NXMAX),INTENT(IN):: U
      COMPLEX(dp),               INTENT(OUT):: F0, DF0
      INTEGER,               INTENT(OUT):: IERR
      REAL(dp)     :: FS, DX
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      F0= U(1,NX) &
     &  + U(2,NX)*DX &
     &  + U(3,NX)*DX*DX &
     &  + U(4,NX)*DX*DX*DX
      DF0= U(2,NX) &
     &   + U(3,NX)*DX*2 &
     &   + U(4,NX)*DX*DX*3
      IERR=0
      RETURN
      END SUBROUTINE CSPL1DD

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of U0  ****

      SUBROUTINE CSPL1DI0(X,U,U0,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,               INTENT(IN) :: NXMAX
      REAL(dp),DIMENSION(NXMAX), INTENT(IN) :: X
      COMPLEX(dp),DIMENSION(4,NXMAX),INTENT(IN):: U
      COMPLEX(dp),DIMENSION(NXMAX), INTENT(OUT):: U0
      INTEGER,               INTENT(OUT):: IERR
      REAL(dp)     :: DX
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF

      U0(1)=0.D0
      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)

         U0(NX)= U0(NX-1) &
     &         + U(1,NX)*DX &
     &         + U(2,NX)*DX*DX/2.D0 &
     &         + U(3,NX)*DX*DX*DX/3.D0 &
     &         + U(4,NX)*DX*DX*DX*DX/4.D0
      ENDDO

!     WRITE(6,'(A,2I5,3ES12.4)')
!     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END SUBROUTINE CSPL1DI0

!     ****** One-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL1DI(X0,FI,X,U,U0,NXMAX,IERR)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                   INTENT(IN) :: NXMAX
      REAL(dp),                      INTENT(IN) :: X0
      REAL(dp),    DIMENSION(NXMAX), INTENT(IN) :: X
      COMPLEX(dp),DIMENSION(4,NXMAX),INTENT(IN) :: U
      COMPLEX(dp),DIMENSION(NXMAX),  INTENT(IN) :: U0
      COMPLEX(dp),                   INTENT(OUT):: FI
      INTEGER,                   INTENT(OUT):: IERR
      REAL(dp)     :: DX, FS
      INTEGER  :: NX


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      FS=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FS*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FS.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FS.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

      DX=X0-X(NX-1)

      FI= U0(NX-1) &
     &  + U(1,NX)*DX &
     &  + U(2,NX)*DX*DX/2.D0 &
     &  + U(3,NX)*DX*DX*DX/3.D0 &
     &  + U(4,NX)*DX*DX*DX*DX/4.D0
!      WRITE(6,'(A,2I5,3ES12.4)')
!     &     'NX,NXI,X0,X(NX-1),X(NX)=',NX,NXI,X0,X(NX-1),X(NX)
      RETURN
      END SUBROUTINE CSPL1DI


!     ***** TRI-DIAGONAL MATRIX SOLVER *****

      SUBROUTINE TDMSRDX(A,X,NMAX,IERR)

!     +++++ INPUT +++++
!        A(4,NMAX) : D : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
!        X(NMAX)   : D : RIGHT HAND SIDE VECTOR
!             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=X(N)
!        NMAX      : I : NUMBER OF EQUATIONS
!     +++++ OUTPUT ++++
!        X(NMAX)   : D : SOLUTION VECTOR
!        IERR      : I : ERROR INDICATOR (0 for naormal)
!     +++++ COMMENT +++
!         This routine requires A(4,NMAX) rather than A(3,NMAX).
!     +++++++++++++++++

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                INTENT(IN)   :: NMAX
      REAL(dp), DIMENSION(4,NMAX),INTENT(INOUT):: A
      REAL(dp), DIMENSION(NMAX),  INTENT(INOUT):: X
      INTEGER,                INTENT(OUT)  :: IERR
      REAL(dp)    :: AAN,BBN,AN,BN,CN,DN,FACT,XN
      INTEGER :: N

      AAN=0.D0
      BBN=0.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=X(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         A(1,N)=AAN
         A(2,N)=BBN
      ENDDO
      XN=BBN

      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         X(N)=AAN*X(N+1)+BBN
      ENDDO
      IERR=0
      RETURN

 9001 IERR=9001
      RETURN
      END SUBROUTINE TDMSRDX

!     ***** PERIODIC TRI-DIAGONAL MATRIX SOLVER *****

      SUBROUTINE TDMPRDX(A,X,NMAX,IERR)

!     +++++ INPUT +++++
!        A(4,NMAX) : D : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
!        XS(NMAX)  : D : RIGHT HAND SIDE VECTOR
!             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=X(N)
!        NMAX      : I : NUMBER OF EQUATIONS
!     +++++ OUTPUT ++++
!        X(NMAX)   : D : SOLUTION VECTOR
!        IERR      : I : ERROR INDICATOR (0 for naormal)
!     +++++ COMMENT +++
!        Do not include periodic points twice.
!           X(0)  =X(N)
!           X(N+1)=X(1)
!        This routine requires A(4,NMAX) rather than A(3,NMAX).
!     +++++++++++++++++

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                INTENT(IN)   :: NMAX
      REAL(dp), DIMENSION(4,NMAX),INTENT(INOUT):: A
      REAL(dp), DIMENSION(NMAX),  INTENT(INOUT):: X
      INTEGER,                INTENT(OUT)  :: IERR
      REAL(dp) :: AAN, BBN, CCN, DDN, EEN, DDNN, EENN, AN, BN, CN, DN, FACT, X1, XN
      INTEGER :: N

      AAN=0.D0
      BBN=0.D0
      CCN=1.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=X(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         CCN=-AN*CCN*FACT
         A(1,N)=AAN
         A(2,N)=BBN
         A(3,N)=CCN
      ENDDO
      IF(CCN.EQ.1.D0) GOTO 9002
      DDN=BBN/(1.D0-CCN)
      EEN=AAN/(1.D0-CCN)
      DDNN=DDN
      EENN=EEN

      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         DDN=AAN*DDN+BBN+CCN*DDNN
         EEN=AAN*EEN    +CCN*EENN
      ENDDO
      IF(EEN.EQ.1.D0) GOTO 9003
      X1=DDN/(1.D0-EEN)
      XN=DDNN+EENN*X1

      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         X(N)=AAN*X(N+1)+BBN+CCN*XN
      ENDDO
      IERR=0
      RETURN

 9001 IERR=9001
      RETURN
 9002 IERR=9002
      RETURN
 9003 IERR=9003
      RETURN
      END SUBROUTINE TDMPRDX

!     ***** TRI-DIAGONAL MATRIX SOLVER *****

      SUBROUTINE TDMSCDX(A,X,NMAX,IERR)

!     +++++ INPUT +++++
!        A(4,NMAX) : D : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
!        X(NMAX)   : D : RIGHT HAND SIDE VECTOR
!             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=X(N)
!        NMAX      : I : NUMBER OF EQUATIONS
!     +++++ OUTPUT ++++
!        X(NMAX)   : D : SOLUTION VECTOR
!        IERR      : I : ERROR INDICATOR (0 for naormal)
!     +++++++++++++++++

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                    INTENT(IN)   :: NMAX
      COMPLEX(dp), DIMENSION(4,NMAX), INTENT(INOUT):: A
      COMPLEX(dp), DIMENSION(NMAX),   INTENT(INOUT):: X
      INTEGER,                    INTENT(OUT)  :: IERR
      COMPLEX(dp) ::  AAN, BBN, AN, BN, CN, DN, FACT, XN
      INTEGER  :: N

      AAN=0.D0
      BBN=0.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=X(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         A(1,N)=AAN
         A(2,N)=BBN
      ENDDO
      XN=BBN

      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         X(N)=AAN*X(N+1)+BBN
      ENDDO
      IERR=0
      RETURN

 9001 IERR=9001
      RETURN
      END SUBROUTINE TDMSCDX

!     ***** PERIODIC TRI-DIAGONAL MATRIX SOLVER *****

      SUBROUTINE TDMPCDX(A,X,NMAX,IERR)

!     +++++ INPUT +++++
!        A(4,NMAX) : D : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
!        XS(NMAX)  : D : RIGHT HAND SIDE VECTOR
!             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=X(N)
!        NMAX      : I : NUMBER OF EQUATIONS
!     +++++ OUTPUT ++++
!        X(NMAX)   : D : SOLUTION VECTOR
!        IERR      : I : ERROR INDICATOR (0 for naormal)
!     +++++ COMMENT +++
!        X(0)  =X(N)
!        X(N+1)=X(1)
!     +++++++++++++++++

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER,                    INTENT(IN)   :: NMAX
      COMPLEX(dp), DIMENSION(4,NMAX), INTENT(INOUT):: A
      COMPLEX(dp), DIMENSION(NMAX),   INTENT(INOUT):: X
      INTEGER,                    INTENT(OUT)  :: IERR
      COMPLEX(dp):: AAN, BBN, CCN, DDN, EEN, DDNN, EENN,  &
           AN, BN, CN, DN, FACT, X1, XN
      INTEGER  :: N


      AAN=0.D0
      BBN=0.D0
      CCN=1.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=X(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         CCN=-AN*CCN*FACT
         A(1,N)=AAN
         A(2,N)=BBN
         A(3,N)=CCN
      ENDDO
      IF(CCN.EQ.1.D0) GOTO 9002
      DDN=BBN/(1.D0-CCN)
      EEN=AAN/(1.D0-CCN)
      DDNN=DDN
      EENN=EEN

      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         DDN=AAN*DDN+BBN+CCN*DDNN
         EEN=AAN*EEN    +CCN*EENN
      ENDDO
      IF(EEN.EQ.1.D0) GOTO 9003
      X1=DDN/(1.D0-EEN)
      XN=DDNN+EENN*X1

      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         X(N)=AAN*X(N+1)+BBN+CCN*XN
      ENDDO
      IERR=0
      RETURN

 9001 IERR=9001
      RETURN
 9002 IERR=9002
      RETURN
 9003 IERR=9003
      RETURN
      END SUBROUTINE TDMPCDX
END MODULE libspl1d
