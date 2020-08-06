!     $Id$
!   ************************************************
!   **           Spline Interpolation             **
!   ************************************************

!     ****** Two-Dimensional Spline Interpolation ******
!       **** Calculation of coefficients ****

      SUBROUTINE SPL2D(X,Y,F,FX,FY,FXY,U,NXM,NXMAX,NYMAX,IDX,IDY,IERR)
!      INPUT : X(NXMAX)        : COORDINATES
!              Y(NYMAX)        : COORDINATES
!              F(NXM,NYMAX)  : VALUE
!              FX(NXM,NYMAX) : EDGE DERIVATIVES FOR IDX != 0
!              FY(NXM,NYMAX) : EDGE DERIVATIVES FOR IDY != 0
!              FXY(NXM,NYMAX): CORNER DERIVATIVES FOR IDY OR IDY != 0
!              NXM       : ARRAY SIZE
!              NXMAX     : NUMBER OF VARIABLES
!              NYMAX     : NUMBER OF VARIABLES
!              IDX       : 0 : SECOND DERIVATIVES = 0 AT X(1) AND X(NXMAX)
!                          1 : DERIVATIVE FX(1) IS GIVEN
!                          2 : DERIVATIVE FX(NXMAX) IS GIVEN
!                          3 : DERIVATIVES FX(1) AXD FX(NXMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NXMAX)
!              IDY       : 0 : SECOND DERIVATIVES = 0 AT Y(1) AND Y(NYMAX)
!                          1 : DERIVATIVE FY(1) IS GIVEN
!                          2 : DERIVATIVE FY(NYMAX) IS GIVEN
!                          3 : DERIVATIVES FY(1) AXD FY(NYMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NYMAX)
!      OUTPUT: U(4,4,NXM,NYMAX): SPLINE COEFICIENTS
!              FX(NXM,NYMAX) : ESTIMATED DERIVATIVES
!              FY(NXM,NYMAX) : ESTIMATED DERIVATIVES
!              FXY(NXM,NYMAX): ESTIMATED DERIVATIVES
!              IERR          : ERROR INDICATOR

      IMPLICIT NONE

      INTEGER(4),                   INTENT(IN) :: NXM, NXMAX, NYMAX, IDX, IDY
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NXM,NYMAX),    INTENT(IN) :: F
      REAL(8), DIMENSION(NXM,NYMAX),  INTENT(INOUT):: FX, FY, FXY
      INTEGER(4),                       INTENT(OUT):: IERR
      REAL(8), DIMENSION(4,4,NXM,NYMAX),INTENT(OUT):: U

      INTEGER(4)            :: NX, NY, IDX1, IDX2, IDY1, IDY2
      REAL(8)               :: DX, DXM, DXP, DX1, DX2, DX3
      REAL(8)               :: DY, DYM, DYP, DY1, DY2, DY3
      REAL(8)               :: T11,T12,T13,T14,T21,T22,T23,T24
      REAL(8)               :: T31,T32,T33,T34,T41,T42,T43,T44
      REAL(8)               :: V13,V14,V23,V24,V31,V32,V33,V34,V41,V42,V43,V44


      REAL(8),    DIMENSION(:,:),ALLOCATABLE:: UX, UX0
      REAL(8),    DIMENSION(:,:),ALLOCATABLE:: UY, UY0
      REAL(8),    DIMENSION(:)  ,ALLOCATABLE:: BX
      REAL(8),    DIMENSION(:)  ,ALLOCATABLE:: BY

      ALLOCATE(UX(4,NXMAX),UX0(4,NXMAX),BX(NXMAX))
      ALLOCATE(UY(4,NYMAX),UY0(4,NYMAX),BY(NYMAX))

!      REAL(8),    DIMENSION(4,NXMAX):: UX, UX0
!      REAL(8),    DIMENSION(4,NYMAX):: UY, UY0
!      REAL(8),    DIMENSION(NXMAX)  :: BX
!      REAL(8),    DIMENSION(NYMAX)  :: BY

!      IF(NXMAX.GT.NMAX) GOTO 9001
!      IF(NYMAX.GT.NMAX) GOTO 9002

      IDX1=MOD(IDX,2)
      IDX2=MOD(IDX/2,2)
      IDY1=MOD(IDY,2)
      IDY2=MOD(IDY/2,2)

!     --- calculate UX ---

      IF(IDX.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         UX(1,1)=DXP
         UX(2,1)=2.D0*(DXP+DXM)
         UX(3,1)=DXM
      ELSE
         IF(IDX1.EQ.0) THEN
            DXP=X(2)-X(1)
            UX(1,1)=0.D0
            UX(2,1)=2.D0*DXP
            UX(3,1)=DXP
         ELSE
            UX(1,1)=0.D0
            UX(2,1)=1.D0
            UX(3,1)=0.D0
         ENDIF
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         UX(1,NX)=DXP
         UX(2,NX)=2.D0*(DXP+DXM)
         UX(3,NX)=DXM
      ENDDO
      IF(IDX.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         UX(1,NXMAX)=DXP
         UX(2,NXMAX)=2.D0*(DXP+DXM)
         UX(3,NXMAX)=DXM
      ELSE
         IF(IDX2.EQ.0) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            UX(1,NXMAX)=DXM
            UX(2,NXMAX)=2.D0*DXM
            UX(3,NXMAX)=0.D0
         ELSE
            UX(1,NXMAX)=0.D0
            UX(2,NXMAX)=1.D0
            UX(3,NXMAX)=0.D0
         ENDIF
      ENDIF
      DO NX=1,NXMAX
         UX0(1,NX)=UX(1,NX)
         UX0(2,NX)=UX(2,NX)
         UX0(3,NX)=UX(3,NX)
      ENDDO

!     --- calculate UY ---

      IF(IDY.EQ.4) THEN
         DYM=Y(NYMAX)-Y(NYMAX-1)
         DYP=Y(2)-Y(1)
         UY(1,1)=DYP
         UY(2,1)=2.D0*(DYP+DYM)
         UY(3,1)=DYM
      ELSE
         IF(IDY1.EQ.0) THEN
            DYP=Y(2)-Y(1)
            UY(1,1)=0.D0
            UY(2,1)=2.D0*DYP
            UY(3,1)=DYP
         ELSE
            UY(1,1)=0.D0
            UY(2,1)=1.D0
            UY(3,1)=0.D0
         ENDIF
      ENDIF
      DO NY=2,NYMAX-1
         DYM=Y(NY)-Y(NY-1)
         DYP=Y(NY+1)-Y(NY)
         UY(1,NY)=DYP
         UY(2,NY)=2.D0*(DYP+DYM)
         UY(3,NY)=DYM
      ENDDO
      IF(IDY.EQ.4) THEN
         DYM=Y(NYMAX)-Y(NYMAX-1)
         DYP=Y(2)-Y(1)
         UY(1,NYMAX)=DYP
         UY(2,NYMAX)=2.D0*(DYP+DYM)
         UY(3,NYMAX)=DYM
      ELSE
         IF(IDY2.EQ.0) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            UY(1,NYMAX)=DYM
            UY(2,NYMAX)=2.D0*DYM
            UY(3,NYMAX)=0.D0
         ELSE
            UY(1,NYMAX)=0.D0
            UY(2,NYMAX)=1.D0
            UY(3,NYMAX)=0.D0
         ENDIF
      ENDIF
      DO NY=1,NYMAX
         UY0(1,NY)=UY(1,NY)
         UY0(2,NY)=UY(2,NY)
         UY0(3,NY)=UY(3,NY)
      ENDDO

!     --- calculate FX ---

      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(F(    2,NY)-F(      1,NY))/DXP &
     &                 +DXP*(F(NXMAX,NY)-F(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(F(2,NY)-F(1,NY))
            ELSE
               BX(1)=FX(1,NY)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(F(NX+1,NY)-F(NX,  NY))/DXP &
     &                  +DXP*(F(NX,  NY)-F(NX-1,NY))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(F(    2,NY)-F(      1,NY))/DXP &
     &                     +DXP*(F(NXMAX,NY)-F(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(F(NXMAX,NY)-F(NXMAX-1,NY))
            ELSE
               BX(NXMAX)=FX(NXMAX,NY)
            ENDIF
         ENDIF

         IF(IDX.EQ.4) THEN
            CALL TDMPRDX(UX,BX,NXMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9003
            BX(NXMAX)=BX(1)
         ELSE
            CALL TDMSRDX(UX,BX,NXMAX,IERR)
            IF(IERR.NE.0) GOTO 9004
         ENDIF
         DO NX=1,NXMAX
            FX(NX,NY)=BX(NX)
         ENDDO
      ENDDO

!     --- calculate FY ---

      DO NX=1,NXMAX
         DO NY=1,NYMAX
            UY(1,NY)=UY0(1,NY)
            UY(2,NY)=UY0(2,NY)
            UY(3,NY)=UY0(3,NY)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(1)=3.D0*(DYM*(F(NX,    2)-F(NX,      1))/DYP &
     &                 +DYP*(F(NX,NYMAX)-F(NX,NYMAX-1))/DYM)
         ELSE
            IF(IDY1.EQ.0) THEN
               BY(1)=3.D0*(F(NX,2)-F(NX,1))
            ELSE
               BY(1)=FY(NX,1)
            ENDIF
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(F(NX,NY+1)-F(NX,NY  ))/DYP &
     &                  +DYP*(F(NX,NY  )-F(NX,NY-1))/DYM)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(NYMAX)=3.D0*(DYM*(F(NX,    2)-F(NX,      1))/DYP &
     &                     +DYP*(F(NX,NYMAX)-F(NX,NYMAX-1))/DYM)
         ELSE
            IF(IDY2.EQ.0) THEN
               BY(NYMAX)=3.D0*(F(NX,NYMAX)-F(NX,NYMAX-1))
            ELSE
               BY(NYMAX)=FY(NX,NYMAX)
            ENDIF
         ENDIF

         IF(IDY.EQ.4) THEN
            CALL TDMPRDX(UY,BY,NYMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9005
            BY(NYMAX)=BY(1)
         ELSE
            CALL TDMSRDX(UY,BY,NYMAX,IERR)
            IF(IERR.NE.0) GOTO 9006
         ENDIF
         DO NY=1,NYMAX
            FY(NX,NY)=BY(NY)
         ENDDO
      ENDDO

!     --- calculate FXY ---

      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(FY(    2,NY)-FY(      1,NY))/DXP &
     &                 +DXP*(FY(NXMAX,NY)-FY(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(FY(2,NY)-FY(1,NY))
            ELSE
               BX(1)=FXY(1,NY)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FY(NX+1,NY)-FY(NX,  NY))/DXP &
     &                  +DXP*(FY(NX,  NY)-FY(NX-1,NY))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(FY(    2,NY)-FY(      1,NY))/DXP &
     &                     +DXP*(FY(NXMAX,NY)-FY(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(FY(NXMAX,NY)-FY(NXMAX-1,NY))
            ELSE
               BX(NXMAX)=FXY(NXMAX,NY)
            ENDIF
         ENDIF

         IF(IDX.EQ.4) THEN
            CALL TDMPRDX(UX,BX,NXMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9007
            BX(NXMAX)=BX(1)
         ELSE
            CALL TDMSRDX(UX,BX,NXMAX,IERR)
            IF(IERR.NE.0) GOTO 9008
         ENDIF
         DO NX=1,NXMAX
            FXY(NX,NY)=BX(NX)
         ENDDO
      ENDDO

!     --- calculate spline coefficients ---

      DO NX=2,NXMAX
      DO NY=2,NYMAX
         DX=X(NX)-X(NX-1)
         DY=Y(NY)-Y(NY-1)

         DY1=1.D0/DY
         DY2=DY1*DY1
         DY3=DY2*DY1
         V13=-3.D0*DY2
         V23=-2.D0*DY1
         V33= 3.D0*DY2
         V43=-1.D0*DY1
         V14= 2.D0*DY3
         V24=      DY2
         V34=-2.D0*DY3
         V44=      DY2

         T11=  F(NX-1,NY-1)
         T12=                    FY(NX-1,NY-1)
         T13=  F(NX-1,NY-1)*V13+ FY(NX-1,NY-1)*V23 &
     &       + F(NX-1,NY  )*V33+ FY(NX-1,NY  )*V43
         T14=  F(NX-1,NY-1)*V14+ FY(NX-1,NY-1)*V24 &
     &       + F(NX-1,NY  )*V34+ FY(NX-1,NY  )*V44
         T21= FX(NX-1,NY-1)
         T22=                   FXY(NX-1,NY-1)
         T23= FX(NX-1,NY-1)*V13+FXY(NX-1,NY-1)*V23 &
     &       +FX(NX-1,NY  )*V33+FXY(NX-1,NY  )*V43
         T24= FX(NX-1,NY-1)*V14+FXY(NX-1,NY-1)*V24 &
     &       +FX(NX-1,NY  )*V34+FXY(NX-1,NY  )*V44
         T31=  F(NX,  NY-1)
         T32=                    FY(NX,  NY-1)
         T33=  F(NX,  NY-1)*V13+ FY(NX,  NY-1)*V23 &
     &       + F(NX,  NY  )*V33+ FY(NX,  NY  )*V43
         T34=  F(NX,  NY-1)*V14+ FY(NX,  NY-1)*V24 &
     &       + F(NX,  NY  )*V34+ FY(NX,  NY  )*V44
         T41= FX(NX,  NY-1)
         T42=                   FXY(NX,  NY-1)
         T43= FX(NX,  NY-1)*V13+FXY(NX,  NY-1)*V23 &
     &       +FX(NX,  NY  )*V33+FXY(NX,  NY  )*V43
         T44= FX(NX,  NY-1)*V14+FXY(NX,  NY-1)*V24 &
     &       +FX(NX,  NY  )*V34+FXY(NX,  NY  )*V44

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

         U(1,1,NX,NY)=    T11
         U(2,1,NX,NY)=            T21
         U(3,1,NX,NY)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,1,NX,NY)=V41*T11+V42*T21+V43*T31+V44*T41
         U(1,2,NX,NY)=    T12
         U(2,2,NX,NY)=            T22
         U(3,2,NX,NY)=V31*T12+V32*T22+V33*T32+V34*T42
         U(4,2,NX,NY)=V41*T12+V42*T22+V43*T32+V44*T42
         U(1,3,NX,NY)=    T13
         U(2,3,NX,NY)=            T23
         U(3,3,NX,NY)=V31*T13+V32*T23+V33*T33+V34*T43
         U(4,3,NX,NY)=V41*T13+V42*T23+V43*T33+V44*T43
         U(1,4,NX,NY)=    T14
         U(2,4,NX,NY)=            T24
         U(3,4,NX,NY)=V31*T14+V32*T24+V33*T34+V34*T44
         U(4,4,NX,NY)=V41*T14+V42*T24+V43*T34+V44*T44
      ENDDO
      ENDDO
      IERR=0
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN

! 9001 WRITE(6,*) 'XX SPL2D: NXMAX.GT.NMAX:',NXMAX,NMAX
!      IERR=1
!      RETURN
! 9002 WRITE(6,*) 'XX SPL2D: NYMAX.GT.NMAX:',NYMAX,NMAX
!      IERR=2
!      RETURN
 9003 WRITE(6,*) 'XX SPL2D: TDMPRDX ERROR: IERR=',IERR
      IERR=3
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN
 9004 WRITE(6,*) 'XX SPL2D: TDMSRDX ERROR: IERR=',IERR
      IERR=4
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN
 9005 WRITE(6,*) 'XX SPL2D: TDMPRDX ERROR: IERR=',IERR
      IERR=5
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN
 9006 WRITE(6,*) 'XX SPL2D: TDMSRDX ERROR: IERR=',IERR
      IERR=6
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN
 9007 WRITE(6,*) 'XX SPL2D: TDMPRDX ERROR: IERR=',IERR
      IERR=7
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN
 9008 WRITE(6,*) 'XX SPL2D: TDMSRDX ERROR: IERR=',IERR
      IERR=8
      DEALLOCATE(UX,UX0,BX,UY,UY0,BY)
      RETURN
      END SUBROUTINE SPL2D

!     ****** Two-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL2DF(X0,Y0,F0,X,Y,U,NXM,NXMAX,NYMAX,IERR)

      IMPLICIT NONE

      INTEGER(4),                       INTENT(IN) :: NXM, NXMAX, NYMAX
      REAL(8),                          INTENT(IN) :: X0, Y0
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(4,4,NXM,NYMAX),INTENT(IN) :: U
      INTEGER(4),                       INTENT(OUT):: IERR
      REAL(8),                          INTENT(OUT):: F0

      INTEGER(4)            :: NX, NY, I
      REAL(8)               :: DX, DY, FX, FY
      REAL(8),DIMENSION(4)  :: UF


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      FX=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FX*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
      FY=1.D0/(Y(NYMAX)-Y(1))
      NY=NINT((Y0-Y(1))*FY*(NYMAX-1))+1
      IF(NY.LT.1) THEN
         IERR=3
         NY=2
      ENDIF
      IF(NY.GT.NYMAX) THEN
         IERR=4
         NY=NYMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FX.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FX.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

 5005 IF(NY.GE.NYMAX) GOTO 5006
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 5006
         NY=NY+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NY.LE.2) GOTO 5008
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 5008
         NY=NY-1
         GOTO 5007
 5008 CONTINUE
      IF(NY.LT.2)     NY=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)

      DO I=1,4
         UF(I)=((U(I,4,NX,NY) *DY     +U(I,3,NX,NY))*DY &
                +U(I,2,NX,NY))*DY     +U(I,1,NX,NY)
      ENDDO
      F0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      IERR=0
      RETURN
      END SUBROUTINE SPL2DF

!     ****** Two-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL2DD(X0,Y0,F0,FX0,FY0,X,Y,U,NXM,NXMAX,NYMAX,IERR)

      IMPLICIT NONE

      INTEGER(4),                       INTENT(IN) :: NXM, NXMAX, NYMAX
      REAL(8),                          INTENT(IN) :: X0, Y0
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(4,4,NXM,NYMAX),INTENT(IN) :: U
      INTEGER(4),                       INTENT(OUT):: IERR
      REAL(8),                          INTENT(OUT):: F0, FX0, FY0

      INTEGER(4)  :: NX, NY, I
      REAL(8)     :: DX, DY, FX, FY
      REAL(8),DIMENSION(4):: UF


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      FX=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FX*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
      FY=1.D0/(Y(NYMAX)-Y(1))
      NY=NINT((Y0-Y(1))*FY*(NYMAX-1))+1
      IF(NY.LT.1) THEN
         IERR=3
         NY=2
      ENDIF
      IF(NY.GT.NYMAX) THEN
         IERR=4
         NY=NYMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FX.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FX.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

 5005 IF(NY.GE.NYMAX) GOTO 5006
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 5006
         NY=NY+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NY.LE.2) GOTO 5008
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 5008
         NY=NY-1
         GOTO 5007
 5008 CONTINUE
      IF(NY.LT.2)     NY=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)

      DO I=1,4
         UF(I)=((U(I,4,NX,NY) *DY +U(I,3,NX,NY))*DY &
                +U(I,2,NX,NY))*DY +U(I,1,NX,NY)
      ENDDO
      F0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)
      FX0= (3*UF(4)*DX+2*UF(3))*DX+UF(2)

      DO I=1,4
         UF(I)=(3*U(I,4,NX,NY)*DY +2*U(I,3,NX,NY))*DY &
                +U(I,2,NX,NY)
      ENDDO
      FY0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      IERR=0
      RETURN
      END SUBROUTINE SPL2DD

!     ****** Two-Dimensional Spline Interpolation ******
!       **** Calculation of coefficients ****

      SUBROUTINE CSPL2D(X,Y,F,FX,FY,FXY,U,NXM,NXMAX,NYMAX,  &
     &                 IDX,IDY,IERR)
!      INPUT : X(NXMAX)        : COORDINATES
!              Y(NYMAX)        : COORDINATES
!              F(NXM,NYMAX)  : VALUE
!              FX(NXM,NYMAX) : EDGE DERIVATIVES FOR IDX != 0
!              FY(NXM,NYMAX) : EDGE DERIVATIVES FOR IDY != 0
!              FXY(NXM,NYMAX): CORNER DERIVATIVES FOR IDY OR IDY != 0
!              NXM       : ARRAY SIZE
!              NXMAX     : NUMBER OF VARIABLES
!              NYMAX     : NUMBER OF VARIABLES
!              IDX       : 0 : SECOND DERIVATIVES = 0 AT X(1) AND X(NXMAX)
!                          1 : DERIVATIVE FX(1) IS GIVEN
!                          2 : DERIVATIVE FX(NXMAX) IS GIVEN
!                          3 : DERIVATIVES FX(1) AND FX(NXMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NXMAX)
!              IDY       : 0 : SECOND DERIVATIVES = 0 AT Y(1) AND Y(NYMAX)
!                          1 : DERIVATIVE FY(1) IS GIVEN
!                          2 : DERIVATIVE FY(NYMAX) IS GIVEN
!                          3 : DERIVATIVES FY(1) AND FY(NYMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NYMAX)
!      OUTPUT: U(4,4,NXM,NYMAX): SPLINE COEFICIENTS
!              FX(NXM,NYMAX) : ESTIMATED DERIVATIVES
!              FY(NXM,NYMAX) : ESTIMATED DERIVATIVES
!              FXY(NXM,NYMAX): ESTIMATED DERIVATIVES
!              IERR          : ERROR INDICATOR

      IMPLICIT NONE

      INTEGER(4),                          INTENT(IN) :: NXM, NXMAX, NYMAX, IDX, IDY
      REAL(8), DIMENSION(NXMAX),           INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),           INTENT(IN) :: Y
      COMPLEX(8), DIMENSION(NXM,NYMAX),    INTENT(IN) :: F
      COMPLEX(8), DIMENSION(NXM,NYMAX),    INTENT(INOUT):: FX, FY, FXY
      INTEGER(4),                          INTENT(OUT):: IERR
      COMPLEX(8), DIMENSION(4,4,NXM,NYMAX),INTENT(OUT):: U


      INTEGER(4) :: NX, NY
      REAL(8)    :: DX, DXM, DXP, DX1, DX2, DX3, DY, DYM, DYP, DY1, DY2, DY3
      INTEGER(4) :: IDX1, IDX2, IDY1, IDY2
      COMPLEX(8) :: T11, T12, T13, T14, T21, T22, T23, T24
      COMPLEX(8) :: T31, T32, T33, T34, T41, T42, T43, T44
      COMPLEX(8) :: V13, V14, V23, V24, V31, V32, V33, V34
      COMPLEX(8) :: V41, V42, V43, V44
!      INTEGER(4), PARAMETER :: NMAX = 1001
!      COMPLEX(8),    DIMENSION(4,NMAX):: UX, UY,UX0, UY0
!      COMPLEX(8),    DIMENSION(NMAX)  ::  BX
      COMPLEX(8),    DIMENSION(4,NXMAX):: UX, UX0
      COMPLEX(8),    DIMENSION(4,NYMAX):: UY, UY0
      COMPLEX(8),    DIMENSION(NXMAX)  :: BX
      COMPLEX(8),    DIMENSION(NYMAX)  :: BY

!      IF(NXMAX.GT.NMAX) GOTO 9001
!      IF(NYMAX.GT.NMAX) GOTO 9002

      IDX1=MOD(IDX,2)
      IDX2=MOD(IDX/2,2)
      IDY1=MOD(IDY,2)
      IDY2=MOD(IDY/2,2)

!     --- calculate UX ---

      IF(IDX.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         UX(1,1)=DXP
         UX(2,1)=2.D0*(DXP+DXM)
         UX(3,1)=DXM
      ELSE
         IF(IDX1.EQ.0) THEN
            DXP=X(2)-X(1)
            UX(1,1)=0.D0
            UX(2,1)=2.D0*DXP
            UX(3,1)=DXP
         ELSE
            UX(1,1)=0.D0
            UX(2,1)=1.D0
            UX(3,1)=0.D0
         ENDIF
      ENDIF
      DO NX=2,NXMAX-1
         DXM=X(NX)-X(NX-1)
         DXP=X(NX+1)-X(NX)
         UX(1,NX)=DXP
         UX(2,NX)=2.D0*(DXP+DXM)
         UX(3,NX)=DXM
      ENDDO
      IF(IDX.EQ.4) THEN
         DXM=X(NXMAX)-X(NXMAX-1)
         DXP=X(2)-X(1)
         UX(1,NXMAX)=DXP
         UX(2,NXMAX)=2.D0*(DXP+DXM)
         UX(3,NXMAX)=DXM
      ELSE
         IF(IDX2.EQ.0) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            UX(1,NXMAX)=DXM
            UX(2,NXMAX)=2.D0*DXM
            UX(3,NXMAX)=0.D0
         ELSE
            UX(1,NXMAX)=0.D0
            UX(2,NXMAX)=1.D0
            UX(3,NXMAX)=0.D0
         ENDIF
      ENDIF
      DO NX=1,NXMAX
         UX0(1,NX)=UX(1,NX)
         UX0(2,NX)=UX(2,NX)
         UX0(3,NX)=UX(3,NX)
      ENDDO

!     --- calculate UY ---

      IF(IDY.EQ.4) THEN
         DYM=Y(NYMAX)-Y(NYMAX-1)
         DYP=Y(2)-Y(1)
         UY(1,1)=DYP
         UY(2,1)=2.D0*(DYP+DYM)
         UY(3,1)=DYM
      ELSE
         IF(IDY1.EQ.0) THEN
            DYP=Y(2)-Y(1)
            UY(1,1)=0.D0
            UY(2,1)=2.D0*DYP
            UY(3,1)=DYP
         ELSE
            UY(1,1)=0.D0
            UY(2,1)=1.D0
            UY(3,1)=0.D0
         ENDIF
      ENDIF
      DO NY=2,NYMAX-1
         DYM=Y(NY)-Y(NY-1)
         DYP=Y(NY+1)-Y(NY)
         UY(1,NY)=DYP
         UY(2,NY)=2.D0*(DYP+DYM)
         UY(3,NY)=DYM
      ENDDO
      IF(IDY.EQ.4) THEN
         DYM=Y(NYMAX)-Y(NYMAX-1)
         DYP=Y(2)-Y(1)
         UY(1,NYMAX)=DYP
         UY(2,NYMAX)=2.D0*(DYP+DYM)
         UY(3,NYMAX)=DYM
      ELSE
         IF(IDY2.EQ.0) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            UY(1,NYMAX)=DYM
            UY(2,NYMAX)=2.D0*DYM
            UY(3,NYMAX)=0.D0
         ELSE
            UY(1,NYMAX)=0.D0
            UY(2,NYMAX)=1.D0
            UY(3,NYMAX)=0.D0
         ENDIF
      ENDIF
      DO NY=1,NYMAX
         UY0(1,NY)=UY(1,NY)
         UY0(2,NY)=UY(2,NY)
         UY0(3,NY)=UY(3,NY)
      ENDDO

!     --- calculate FX ---

      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(F(    2,NY)-F(      1,NY))/DXP   &
     &                 +DXP*(F(NXMAX,NY)-F(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(F(2,NY)-F(1,NY))
            ELSE
               BX(1)=FX(1,NY)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(F(NX+1,NY)-F(NX,  NY))/DXP   &
     &                  +DXP*(F(NX,  NY)-F(NX-1,NY))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(F(    2,NY)-F(      1,NY))/DXP   &
     &                     +DXP*(F(NXMAX,NY)-F(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(F(NXMAX,NY)-F(NXMAX-1,NY))
            ELSE
               BX(NXMAX)=FX(NXMAX,NY)
            ENDIF
         ENDIF

         IF(IDX.EQ.4) THEN
            CALL TDMPCDX(UX,BX,NXMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9003
            BX(NXMAX)=BX(1)
         ELSE
            CALL TDMSCDX(UX,BX,NXMAX,IERR)
            IF(IERR.NE.0) GOTO 9004
         ENDIF
         DO NX=1,NXMAX
            FX(NX,NY)=BX(NX)
         ENDDO
      ENDDO

!     --- calculate FY ---

      DO NX=1,NXMAX
         DO NY=1,NYMAX
            UY(1,NY)=UY0(1,NY)
            UY(2,NY)=UY0(2,NY)
            UY(3,NY)=UY0(3,NY)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(1)=3.D0*(DYM*(F(NX,    2)-F(NX,      1))/DYP  &
     &                 +DYP*(F(NX,NYMAX)-F(NX,NYMAX-1))/DYM)
         ELSE
            IF(IDY1.EQ.0) THEN
               BY(1)=3.D0*(F(NX,2)-F(NX,1))
            ELSE
               BY(1)=FY(NX,1)
            ENDIF
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(F(NX,NY+1)-F(NX,NY  ))/DYP  &
     &                  +DYP*(F(NX,NY  )-F(NX,NY-1))/DYM)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(NYMAX)=3.D0*(DYM*(F(NX,    2)-F(NX,      1))/DYP  &
     &                     +DYP*(F(NX,NYMAX)-F(NX,NYMAX-1))/DYM)
         ELSE
            IF(IDY2.EQ.0) THEN
               BY(NYMAX)=3.D0*(F(NX,NYMAX)-F(NX,NYMAX-1))
            ELSE
               BY(NYMAX)=FY(NX,NYMAX)
            ENDIF
         ENDIF

         IF(IDY.EQ.4) THEN
            CALL TDMPCDX(UY,BY,NYMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9005
            BY(NYMAX)=BY(1)
         ELSE
            CALL TDMSCDX(UY,BY,NYMAX,IERR)
            IF(IERR.NE.0) GOTO 9006
         ENDIF
         DO NY=1,NYMAX
            FY(NX,NY)=BY(NY)
         ENDDO
      ENDDO

!     --- calculate FXY ---

      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(FY(    2,NY)-FY(      1,NY))/DXP  &
     &                 +DXP*(FY(NXMAX,NY)-FY(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(FY(2,NY)-FY(1,NY))
            ELSE
               BX(1)=FXY(1,NY)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FY(NX+1,NY)-FY(NX,  NY))/DXP  &
     &                  +DXP*(FY(NX,  NY)-FY(NX-1,NY))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(FY(    2,NY)-FY(      1,NY))/DXP  &
     &                     +DXP*(FY(NXMAX,NY)-FY(NXMAX-1,NY))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(FY(NXMAX,NY)-FY(NXMAX-1,NY))
            ELSE
               BX(NXMAX)=FXY(NXMAX,NY)
            ENDIF
         ENDIF

         IF(IDX.EQ.4) THEN
            CALL TDMPCDX(UX,BX,NXMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9007
            BX(NXMAX)=BX(1)
         ELSE
            CALL TDMSCDX(UX,BX,NXMAX,IERR)
            IF(IERR.NE.0) GOTO 9008
         ENDIF
         DO NX=1,NXMAX
            FXY(NX,NY)=BX(NX)
         ENDDO
      ENDDO

!     --- calculate spline coefficients ---

      DO NX=2,NXMAX
      DO NY=2,NYMAX
         DX=X(NX)-X(NX-1)
         DY=Y(NY)-Y(NY-1)

         DY1=1.D0/DY
         DY2=DY1*DY1
         DY3=DY2*DY1
         V13=-3.D0*DY2
         V23=-2.D0*DY1
         V33= 3.D0*DY2
         V43=-1.D0*DY1
         V14= 2.D0*DY3
         V24=      DY2
         V34=-2.D0*DY3
         V44=      DY2

         T11=  F(NX-1,NY-1)
         T12=                    FY(NX-1,NY-1)
         T13=  F(NX-1,NY-1)*V13+ FY(NX-1,NY-1)*V23  &
     &       + F(NX-1,NY  )*V33+ FY(NX-1,NY  )*V43
         T14=  F(NX-1,NY-1)*V14+ FY(NX-1,NY-1)*V24  &
     &       + F(NX-1,NY  )*V34+ FY(NX-1,NY  )*V44
         T21= FX(NX-1,NY-1)
         T22=                   FXY(NX-1,NY-1)
         T23= FX(NX-1,NY-1)*V13+FXY(NX-1,NY-1)*V23  &
     &       +FX(NX-1,NY  )*V33+FXY(NX-1,NY  )*V43
         T24= FX(NX-1,NY-1)*V14+FXY(NX-1,NY-1)*V24  &
     &       +FX(NX-1,NY  )*V34+FXY(NX-1,NY  )*V44
         T31=  F(NX,  NY-1)
         T32=                    FY(NX,  NY-1)
         T33=  F(NX,  NY-1)*V13+ FY(NX,  NY-1)*V23  &
     &       + F(NX,  NY  )*V33+ FY(NX,  NY  )*V43
         T34=  F(NX,  NY-1)*V14+ FY(NX,  NY-1)*V24  &
     &       + F(NX,  NY  )*V34+ FY(NX,  NY  )*V44
         T41= FX(NX,  NY-1)
         T42=                   FXY(NX,  NY-1)
         T43= FX(NX,  NY-1)*V13+FXY(NX,  NY-1)*V23  &
     &       +FX(NX,  NY  )*V33+FXY(NX,  NY  )*V43
         T44= FX(NX,  NY-1)*V14+FXY(NX,  NY-1)*V24  &
     &       +FX(NX,  NY  )*V34+FXY(NX,  NY  )*V44

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

         U(1,1,NX,NY)=    T11
         U(2,1,NX,NY)=            T21
         U(3,1,NX,NY)=V31*T11+V32*T21+V33*T31+V34*T41
         U(4,1,NX,NY)=V41*T11+V42*T21+V43*T31+V44*T41
         U(1,2,NX,NY)=    T12
         U(2,2,NX,NY)=            T22
         U(3,2,NX,NY)=V31*T12+V32*T22+V33*T32+V34*T42
         U(4,2,NX,NY)=V41*T12+V42*T22+V43*T32+V44*T42
         U(1,3,NX,NY)=    T13
         U(2,3,NX,NY)=            T23
         U(3,3,NX,NY)=V31*T13+V32*T23+V33*T33+V34*T43
         U(4,3,NX,NY)=V41*T13+V42*T23+V43*T33+V44*T43
         U(1,4,NX,NY)=    T14
         U(2,4,NX,NY)=            T24
         U(3,4,NX,NY)=V31*T14+V32*T24+V33*T34+V34*T44
         U(4,4,NX,NY)=V41*T14+V42*T24+V43*T34+V44*T44
      ENDDO
      ENDDO
      IERR=0
      RETURN

! 9001 WRITE(6,*) 'XX CSPL2D: NXMAX.GT.NMAX:',NXMAX,NMAX
!      IERR=1
!      RETURN
! 9002 WRITE(6,*) 'XX CSPL2D: NYMAX.GT.NMAX:',NYMAX,NMAX
!      IERR=2
1      RETURN
 9003 WRITE(6,*) 'XX CSPL2D: TDMSCDX ERROR: IERR=',IERR
      IERR=3
      RETURN
 9004 WRITE(6,*) 'XX CSPL2D: TDMSCDX ERROR: IERR=',IERR
      IERR=4
      RETURN
 9005 WRITE(6,*) 'XX CSPL2D: TDMSCDX ERROR: IERR=',IERR
      IERR=5
      RETURN
 9006 WRITE(6,*) 'XX CSPL2D: TDMSCDX ERROR: IERR=',IERR
      IERR=6
      RETURN
 9007 WRITE(6,*) 'XX CSPL2D: TDMSCDX ERROR: IERR=',IERR
      IERR=7
      RETURN
 9008 WRITE(6,*) 'XX CSPL2D: TDMSCDX ERROR: IERR=',IERR
      IERR=8
      RETURN
      END SUBROUTINE CSPL2D

!     ****** Two-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL2DF(X0,Y0,F0,X,Y,U,NXM,NXMAX,NYMAX,IERR)

      IMPLICIT NONE

      INTEGER(4),                          INTENT(IN) :: NXM, NXMAX, NYMAX
      REAL(8),                             INTENT(IN) :: X0, Y0
      REAL(8), DIMENSION(NXMAX),           INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),           INTENT(IN) :: Y
      COMPLEX(8), DIMENSION(4,4,NXM,NYMAX),INTENT(IN) :: U
      INTEGER(4),                          INTENT(OUT):: IERR
      COMPLEX(8),                          INTENT(OUT):: F0

      INTEGER(4)            :: NX, NY
      REAL(8)               :: DX, DY, FX, FY
      COMPLEX(8)            :: F1, F2, F3, F4


      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      FX=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FX*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
      FY=1.D0/(Y(NYMAX)-Y(1))
      NY=NINT((Y0-Y(1))*FY*(NYMAX-1))+1
      IF(NY.LT.1) THEN
         IERR=3
         NY=2
      ENDIF
      IF(NY.GT.NYMAX) THEN
         IERR=4
         NY=NYMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FX.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FX.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

 5005 IF(NY.GE.NYMAX) GOTO 5006
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 5006
         NY=NY+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NY.LE.2) GOTO 5008
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 5008
         NY=NY-1
         GOTO 5007
 5008 CONTINUE
      IF(NY.LT.2)     NY=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)

!      WRITE(6,'(2I5,2ES12.4)') NX,NY,X0,Y0
!      IF(NX.EQ.33 .AND. NY.EQ.18) THEN
!         WRITE(6,'(4ES12.4)') X0,Y0,DX,DY
!         WRITE(6,'(4ES12.4)') U(1,1,NX,NY),U(1,2,NX,NY),
!     &                         U(1,3,NX,NY),U(1,4,NX,NY)
!         WRITE(6,'(4ES12.4)') U(2,1,NX,NY),U(2,2,NX,NY),
!     &                         U(2,3,NX,NY),U(2,4,NX,NY)
!         WRITE(6,'(4ES12.4)') U(3,1,NX,NY),U(3,2,NX,NY),
!     &                         U(3,3,NX,NY),U(3,4,NX,NY)
!         WRITE(6,'(4ES12.4)') U(4,1,NX,NY),U(4,2,NX,NY),
!     &                         U(4,3,NX,NY),U(4,4,NX,NY)
!         CALL GUFLSH
!      ENDIF
!
!      F0= U(1,1,NX,NY)      +U(1,2,NX,NY)*DY
!     &   +U(1,3,NX,NY)*DY*DY+U(1,4,NX,NY)*DY*DY*DY
!     &  +(U(2,1,NX,NY)      +U(2,2,NX,NY)*DY
!     &   +U(2,3,NX,NY)*DY*DY+U(2,4,NX,NY)*DY*DY*DY)*DX
!     &  +(U(3,1,NX,NY)      +U(3,2,NX,NY)*DY
!     &   +U(3,3,NX,NY)*DY*DY+U(3,4,NX,NY)*DY*DY*DY)*DX*DX
!     &  +(U(4,1,NX,NY)      +U(4,2,NX,NY)*DY
!     &   +U(4,3,NX,NY)*DY*DY+U(4,4,NX,NY)*DY*DY*DY)*DX*DX*DX

      F4= ((U(4,4,NX,NY) *DY     +U(4,3,NX,NY))*DY   &
     &     +U(4,2,NX,NY))*DY     +U(4,1,NX,NY)
      F3= ((U(3,4,NX,NY) *DY     +U(3,3,NX,NY))*DY   &
     &     +U(3,2,NX,NY))*DY     +U(3,1,NX,NY)
      F2= ((U(2,4,NX,NY) *DY     +U(2,3,NX,NY))*DY   &
     &     +U(2,2,NX,NY))*DY     +U(2,1,NX,NY)
      F1= ((U(1,4,NX,NY) *DY     +U(1,3,NX,NY))*DY   &
     &     +U(1,2,NX,NY))*DY     +U(1,1,NX,NY)
      F0= ((F4*DX+F3)*DX+F2)*DX+F1

      IERR=0
      RETURN
      END SUBROUTINE CSPL2DF

!     ****** Two-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL2DD(X0,Y0,F0,FX0,FY0,X,Y,U,NXM,NXMAX,NYMAX,IERR)

      IMPLICIT NONE

      INTEGER(4),                          INTENT(IN) :: NXM, NXMAX, NYMAX
      REAL(8),                             INTENT(IN) :: X0, Y0
      REAL(8), DIMENSION(NXMAX),           INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),           INTENT(IN) :: Y
      COMPLEX(8), DIMENSION(4,4,NXM,NYMAX),INTENT(IN) :: U
      INTEGER(4),                          INTENT(OUT):: IERR
      COMPLEX(8),                          INTENT(OUT):: F0, FX0, FY0

      INTEGER(4)  :: NX, NY
      REAL(8)     :: DX, DY, FX, FY
      COMPLEX(8)  :: F1, F2, F3, F4, FY1, FY2, FY3, FY4

      IERR=0
      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      FX=1.D0/(X(NXMAX)-X(1))
      NX=NINT((X0-X(1))*FX*(NXMAX-1))+1
      IF(NX.LT.1) THEN
         IERR=1
         NX=2
      ENDIF
      IF(NX.GT.NXMAX) THEN
         IERR=2
         NX=NXMAX
      ENDIF
      FY=1.D0/(Y(NYMAX)-Y(1))
      NY=NINT((Y0-Y(1))*FY*(NYMAX-1))+1
      IF(NY.LT.1) THEN
         IERR=3
         NY=2
      ENDIF
      IF(NY.GT.NYMAX) THEN
         IERR=4
         NY=NYMAX
      ENDIF

 5001 IF(NX.GE.NXMAX) GOTO 5002
      IF((X0-X(NX  ))*FX.LE.0.D0) GOTO 5002
         NX=NX+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NX.LE.2) GOTO 5004
      IF((X0-X(NX-1))*FX.GE.0.D0) GOTO 5004
         NX=NX-1
         GOTO 5003
 5004 CONTINUE
      IF(NX.LT.2)     NX=2

 5005 IF(NY.GE.NYMAX) GOTO 5006
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 5006
         NY=NY+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NY.LE.2) GOTO 5008
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 5008
         NY=NY-1
         GOTO 5007
 5008 CONTINUE
      IF(NY.LT.2)     NY=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)

      F4= ((U(4,4,NX,NY) *DY     +U(4,3,NX,NY))*DY  &
     &     +U(4,2,NX,NY))*DY     +U(4,1,NX,NY)
      F3= ((U(3,4,NX,NY) *DY     +U(3,3,NX,NY))*DY  &
     &     +U(3,2,NX,NY))*DY     +U(3,1,NX,NY)
      F2= ((U(2,4,NX,NY) *DY     +U(2,3,NX,NY))*DY  &
     &     +U(2,2,NX,NY))*DY     +U(2,1,NX,NY)
      F1= ((U(1,4,NX,NY) *DY     +U(1,3,NX,NY))*DY  &
     &     +U(1,2,NX,NY))*DY     +U(1,1,NX,NY)

      F0= ((F4*DX+F3)*DX+F2)*DX+F1
      FX0= (3*F4*DX+2*F3)*DX+F2

      FY4= (3*U(4,4,NX,NY) *DY     +2*U(4,3,NX,NY))*DY  &
     &       +U(4,2,NX,NY)
      FY3= (3*U(3,4,NX,NY) *DY     +2*U(3,3,NX,NY))*DY  &
     &       +U(3,2,NX,NY)
      FY2= (3*U(2,4,NX,NY) *DY     +2*U(2,3,NX,NY))*DY  &
     &       +U(2,2,NX,NY)
      FY1= (3*U(1,4,NX,NY) *DY     +2*U(1,3,NX,NY))*DY  &
     &       +U(1,2,NX,NY)
      FY0= ((FY4*DX+FY3)*DX+FY2)*DX+FY1
      IERR=0
      RETURN
      END SUBROUTINE CSPL2DD
