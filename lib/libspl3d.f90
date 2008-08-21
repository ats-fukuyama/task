!     $Id$
!
!     ****** Three-Dimensional Spline Interpolation ******
!       **** Calculation of coefficients ****

      SUBROUTINE SPL3D(X,Y,Z,F,FX,FY,FZ,FXY,FYZ,FZX,FXYZ,U, &
                       NXM,NYM,NXMAX,NYMAX,NZMAX,IDX,IDY,IDZ,IERR)
!      INPUT : X(NXMAX)        : COORDINATES
!              Y(NYMAX)        : COORDINATES
!              Z(NZMAX)        : COORDINATES
!              F(NXM,NYM,NZMAX)   : VALUE
!              FX(NXM,NYM,NZMAX)  : EDGE DERIVATIVES FOR IDX != 0
!              FY(NXM,NYM,NZMAX)  : EDGE DERIVATIVES FOR IDY != 0
!              FZ(NXM,NYM,NZMAX)  : EDGE DERIVATIVES FOR IDZ != 0
!              FXY(NXM,NYM,NZMAX) : CORNER DERIVATIVES FOR IDX OR IDY != 0
!              FYZ(NXM,NYM,NZMAX) : CORNER DERIVATIVES FOR IDY OR IDZ != 0
!              FZX(NXM,NYM,NZMAX) : CORNER DERIVATIVES FOR IDZ OR IDX != 0
!              FXYZ(NXM,NYM,NZMAX): CORNER DERIVS FOR IDX OR IDY OR IDZ != 0
!              NXM       : ARRAY SIZE
!              NYM       : ARRAY SIZE
!              NXMAX     : NUMBER OF VARIABLES
!              NYMAX     : NUMBER OF VARIABLES
!              NZMAX     : NUMBER OF VARIABLES
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
!              IDZ       : 0 : SECOND DERIVATIVES = 0 AT Z(1) AND Z(NZMAX)
!                          1 : DERIVATIVE FZ(1) IS GIVEN
!                          2 : DERIVATIVE FZ(NYMAX) IS GIVEN
!                          3 : DERIVATIVES FZ(1) AXD FZ(NYMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NZMAX)
!      OUTPUT: U(4,4,4,NXM,NYM,NZMAX): SPLINE COEFICIENTS
!              FX(NXM,NYM,NZMAX) : ESTIMATED DERIVATIVES
!              FY(NXM,NYM,NZMAX) : ESTIMATED DERIVATIVES
!              FZ(NXM,NYM,NZMAX) : ESTIMATED DERIVATIVES
!              FXY(NXM,NYM,NZMAX): ESTIMATED DERIVATIVES
!              FYZ(NXM,NYM,NZMAX): ESTIMATED DERIVATIVES
!              FZX(NXM,NYM,NZMAX): ESTIMATED DERIVATIVES
!              IERR          : ERROR INDICATOR

      IMPLICIT NONE

      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NZMAX),        INTENT(IN) :: Z
      REAL(8), DIMENSION(NXM,NYM,NZMAX),INTENT(IN) :: F
      REAL(8), DIMENSION(NXM,NYM,NZMAX),INTENT(INOUT):: FX,FY,FZ
      REAL(8), DIMENSION(NXM,NYM,NZMAX),INTENT(INOUT):: FXY,FYZ,FZX,FXYZ
      REAL(8), DIMENSION(4,4,4,NXM,NYM,NZMAX),INTENT(OUT):: U
      INTEGER(4),INTENT(IN) :: NXM,NYM,NXMAX,NYMAX,NZMAX,IDX,IDY,IDZ
      INTEGER(4),INTENT(OUT):: IERR

      INTEGER(4) :: NX,NY,NZ,I,J,K,L
      REAL(8)    :: DX,DXM,DXP,DX1,DX2,DX3
      REAL(8)    :: DY,DYM,DYP,DY1,DY2,DY3
      REAL(8)    :: DZ,DZM,DZP,DZ1,DZ2,DZ3
      INTEGER(4)    :: IDX1,IDX2,IDY1,IDY2,IDZ1,IDZ2
      REAL(8),DIMENSION(4,4):: VX,VY,VZ
      REAL(8),DIMENSION(4,4,4):: T,TT,TTT,TTTT
      REAL(8),    DIMENSION(4,NXMAX):: UX,UX0
      REAL(8),    DIMENSION(4,NYMAX):: UY,UY0
      REAL(8),    DIMENSION(4,NZMAX):: UZ,UZ0
      REAL(8),    DIMENSION(NXMAX)  :: BX
      REAL(8),    DIMENSION(NYMAX)  :: BY
      REAL(8),    DIMENSION(NZMAX)  :: BZ

      IDX1=MOD(IDX,2)
      IDX2=MOD(IDX/2,2)
      IDY1=MOD(IDY,2)
      IDY2=MOD(IDY/2,2)
      IDZ1=MOD(IDZ,2)
      IDZ2=MOD(IDZ/2,2)

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

!     --- calculate UZ ---

      IF(IDZ.EQ.4) THEN
         DZM=Z(NZMAX)-Z(NZMAX-1)
         DZP=Z(2)-Z(1)
         UZ(1,1)=DZP
         UZ(2,1)=2.D0*(DZP+DZM)
         UZ(3,1)=DZM
      ELSE
         IF(IDZ1.EQ.0) THEN
            DZP=Z(2)-Z(1)
            UZ(1,1)=0.D0
            UZ(2,1)=2.D0*DZP
            UZ(3,1)=DZP
         ELSE
            UZ(1,1)=0.D0
            UZ(2,1)=1.D0
            UZ(3,1)=0.D0
         ENDIF
      ENDIF
      DO NZ=2,NZMAX-1
         DZM=Z(NZ)-Z(NZ-1)
         DZP=Z(NZ+1)-Z(NZ)
         UZ(1,NZ)=DZP
         UZ(2,NZ)=2.D0*(DZP+DZM)
         UZ(3,NZ)=DZM
      ENDDO
      IF(IDZ.EQ.4) THEN
         DZM=Z(NZMAX)-Z(NZMAX-1)
         DZP=Z(2)-Z(1)
         UZ(1,NZMAX)=DZP
         UZ(2,NZMAX)=2.D0*(DZP+DZM)
         UZ(3,NZMAX)=DZM
      ELSE
         IF(IDZ2.EQ.0) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            UZ(1,NZMAX)=DZM
            UZ(2,NZMAX)=2.D0*DZM
            UZ(3,NZMAX)=0.D0
         ELSE
            UZ(1,NZMAX)=0.D0
            UZ(2,NZMAX)=1.D0
            UZ(3,NZMAX)=0.D0
         ENDIF
      ENDIF
      DO NZ=1,NZMAX
         UZ0(1,NZ)=UZ(1,NZ)
         UZ0(2,NZ)=UZ(2,NZ)
         UZ0(3,NZ)=UZ(3,NZ)
      ENDDO

!     --- calculate FX ---

      DO NZ=1,NZMAX
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(F(    2,NY,NZ)-F(      1,NY,NZ))/DXP &
     &                 +DXP*(F(NXMAX,NY,NZ)-F(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(F(2,NY,NZ)-F(1,NY,NZ))
            ELSE
               BX(1)=FX(1,NY,NZ)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(F(NX+1,NY,NZ)-F(NX,  NY,NZ))/DXP &
     &                  +DXP*(F(NX,  NY,NZ)-F(NX-1,NY,NZ))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(F(    2,NY,NZ)-F(      1,NY,NZ))/DXP &
     &                     +DXP*(F(NXMAX,NY,NZ)-F(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(F(NXMAX,NY,NZ)-F(NXMAX-1,NY,NZ))
            ELSE
               BX(NXMAX)=FX(NXMAX,NY,NZ)
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
            FX(NX,NY,NZ)=BX(NX)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FY ---

      DO NZ=1,NZMAX
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            UY(1,NY)=UY0(1,NY)
            UY(2,NY)=UY0(2,NY)
            UY(3,NY)=UY0(3,NY)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(1)=3.D0*(DYM*(F(NX,    2,NZ)-F(NX,      1,NZ))/DYP &
     &                 +DYP*(F(NX,NYMAX,NZ)-F(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY1.EQ.0) THEN
               BY(1)=3.D0*(F(NX,2,NZ)-F(NX,1,NZ))
            ELSE
               BY(1)=FY(NX,1,NZ)
            ENDIF
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(F(NX,NY+1,NZ)-F(NX,NY  ,NZ))/DYP &
     &                  +DYP*(F(NX,NY  ,NZ)-F(NX,NY-1,NZ))/DYM)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(NYMAX)=3.D0*(DYM*(F(NX,    2,NZ)-F(NX,      1,NZ))/DYP &
     &                     +DYP*(F(NX,NYMAX,NZ)-F(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY2.EQ.0) THEN
               BY(NYMAX)=3.D0*(F(NX,NYMAX,NZ)-F(NX,NYMAX-1,NZ))
            ELSE
               BY(NYMAX)=FY(NX,NYMAX,NZ)
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
            FY(NX,NY,NZ)=BY(NY)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FZ ---

      DO NY=1,NYMAX
      DO NX=1,NXMAX
         DO NZ=1,NZMAX
            UZ(1,NZ)=UZ0(1,NZ)
            UZ(2,NZ)=UZ0(2,NZ)
            UZ(3,NZ)=UZ0(3,NZ)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(1)=3.D0*(DZM*(F(NX,NY,    2)-F(NX,NY,      1))/DZP &
     &                 +DZP*(F(NX,NY,NZMAX)-F(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ1.EQ.0) THEN
               BZ(1)=3.D0*(F(NX,NY,2)-F(NX,NY,1))
            ELSE
               BZ(1)=FZ(NX,NY,1)
            ENDIF
         ENDIF
         DO NZ=2,NZMAX-1
            DZM=Z(NZ)-Z(NZ-1)
            DZP=Z(NZ+1)-Z(NZ)
            BZ(NZ)=3.D0*(DZM*(F(NX,NY,NZ+1)-F(NX,NY,NZ  ))/DZP &
     &                  +DZP*(F(NX,NY,NZ  )-F(NX,NY,NZ-1))/DZM)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(NZMAX)=3.D0*(DZM*(F(NX,NY,    2)-F(NX,NY,      1))/DZP &
     &                     +DZP*(F(NX,NY,NZMAX)-F(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ2.EQ.0) THEN
               BZ(NZMAX)=3.D0*(F(NX,NY,NZMAX)-F(NX,NY,NZMAX-1))
            ELSE
               BZ(NZMAX)=FZ(NX,NY,NZMAX)
            ENDIF
         ENDIF

         IF(IDZ.EQ.4) THEN
            CALL TDMPRDX(UZ,BZ,NZMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9005
            BZ(NZMAX)=BZ(1)
         ELSE
            CALL TDMSRDX(UZ,BZ,NZMAX,IERR)
            IF(IERR.NE.0) GOTO 9006
         ENDIF
         DO NZ=1,NZMAX
            FZ(NX,NY,NZ)=BZ(NZ)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FXY ---

      DO NZ=1,NZMAX
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(FY(    2,NY,NZ)-FY(      1,NY,NZ))/DXP &
     &                 +DXP*(FY(NXMAX,NY,NZ)-FY(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(FY(2,NY,NZ)-FY(1,NY,NZ))
            ELSE
               BX(1)=FXY(1,NY,NZ)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FY(NX+1,NY,NZ)-FY(NX,  NY,NZ))/DXP &
     &                  +DXP*(FY(NX,  NY,NZ)-FY(NX-1,NY,NZ))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(FY(    2,NY,NZ)-FY(      1,NY,NZ))/DXP &
     &                     +DXP*(FY(NXMAX,NY,NZ)-FY(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(FY(NXMAX,NY,NZ)-FY(NXMAX-1,NY,NZ))
            ELSE
               BX(NXMAX)=FXY(NXMAX,NY,NZ)
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
            FXY(NX,NY,NZ)=BX(NX)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FYZ ---

      DO NZ=1,NZMAX
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            UY(1,NY)=UY0(1,NY)
            UY(2,NY)=UY0(2,NY)
            UY(3,NY)=UY0(3,NY)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(1)=3.D0*(DYM*(FZ(NX,    2,NZ)-FZ(NX,      1,NZ))/DYP &
     &                 +DYP*(FZ(NX,NYMAX,NZ)-FZ(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY1.EQ.0) THEN
               BY(1)=3.D0*(FZ(NX,2,NZ)-FZ(NX,1,NZ))
            ELSE
               BY(1)=FYZ(NX,1,NZ)
            ENDIF
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(FZ(NX,NY+1,NZ)-FZ(NX,NY,  NZ))/DYP &
     &                  +DYP*(FZ(NX,NY,  NZ)-FZ(NX,NY-1,NZ))/DYM)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(NYMAX)=3.D0*(DYM*(FZ(NX,    2,NZ)-FZ(NX,      1,NZ))/DYP &
     &                     +DYP*(FZ(NX,NYMAX,NZ)-FZ(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY2.EQ.0) THEN
               BY(NYMAX)=3.D0*(FZ(NX,NYMAX,NZ)-FZ(NX,NYMAX-1,NZ))
            ELSE
               BY(NYMAX)=FYZ(NX,NYMAX,NZ)
            ENDIF
         ENDIF

         IF(IDY.EQ.4) THEN
            CALL TDMPRDX(UY,BY,NYMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9007
            BY(NYMAX)=BY(1)
         ELSE
            CALL TDMSRDX(UY,BY,NYMAX,IERR)
            IF(IERR.NE.0) GOTO 9008
         ENDIF
         DO NY=1,NYMAX
            FYZ(NX,NY,NZ)=BY(NY)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FZX ---

      DO NY=1,NYMAX
      DO NX=1,NXMAX
         DO NZ=1,NZMAX
            UZ(1,NZ)=UZ0(1,NZ)
            UZ(2,NZ)=UZ0(2,NZ)
            UZ(3,NZ)=UZ0(3,NZ)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(1)=3.D0*(DZM*(FX(NX,NY,    2)-FX(NX,NY,      1))/DZP &
     &                 +DZP*(FX(NX,NY,NZMAX)-FX(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ1.EQ.0) THEN
               BZ(1)=3.D0*(FX(NX,NY,2)-FX(NX,NY,1))
            ELSE
               BZ(1)=FZX(NX,NY,1)
            ENDIF
         ENDIF
         DO NZ=2,NZMAX-1
            DZM=Z(NZ)-Z(NZ-1)
            DZP=Z(NZ+1)-Z(NZ)
            BZ(NZ)=3.D0*(DZM*(FX(NX,NY,NZ+1)-FX(NX,NY,NZ  ))/DZP &
     &                  +DZP*(FX(NX,NY,NZ  )-FX(NX,NY,NZ-1))/DZM)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(NZMAX)=3.D0*(DZM*(FX(NX,NY,    2)-FX(NX,NY,      1))/DZP &
     &                     +DZP*(FX(NX,NY,NZMAX)-FX(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ2.EQ.0) THEN
               BZ(NZMAX)=3.D0*(FX(NX,NY,NZMAX)-FX(NX,NY,NZMAX-1))
            ELSE
               BZ(NZMAX)=FZX(NX,NY,NZMAX)
            ENDIF
         ENDIF

         IF(IDZ.EQ.4) THEN
            CALL TDMPRDX(UZ,BZ,NZMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9007
            BZ(NZMAX)=BZ(1)
         ELSE
            CALL TDMSRDX(UZ,BZ,NZMAX,IERR)
            IF(IERR.NE.0) GOTO 9008
         ENDIF
         DO NZ=1,NZMAX
            FZX(NX,NY,NZ)=BZ(NZ)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FXYZ ---
      write(6,*) '--- calculate FXYZ ---'

      DO NZ=1,NZMAX
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(FYZ(    2,NY,NZ)-FYZ(      1,NY,NZ))/DXP &
     &                 +DXP*(FYZ(NXMAX,NY,NZ)-FYZ(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(FYZ(2,NY,NZ)-FYZ(1,NY,NZ))
            ELSE
               BX(1)=FXYZ(1,NY,NZ)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FYZ(NX+1,NY,NZ)-FYZ(NX,  NY,NZ))/DXP &
     &                  +DXP*(FYZ(NX,  NY,NZ)-FYZ(NX-1,NY,NZ))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(FYZ(    2,NY,NZ)-FYZ(      1,NY,NZ))/DXP &
     &                     +DXP*(FYZ(NXMAX,NY,NZ)-FYZ(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(FYZ(NXMAX,NY,NZ)-FYZ(NXMAX-1,NY,NZ))
            ELSE
               BX(NXMAX)=FXYZ(NXMAX,NY,NZ)
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
            FXYZ(NX,NY,NZ)=BX(NX)
         ENDDO
      ENDDO
      ENDDO

!      DO NX=1,NXMAX
!      DO NY=1,NYMAX
!      DO NZ=1,NZMAX
!         WRITE(6,'(3I3,1P7E10.2)') NX,NY,NZ, &
!              F(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ), &
!              FXY(NX,NY,NZ),FYZ(NX,NY,NZ),FZX(NX,NY,NZ)
!      ENDDO
!      ENDDO
!      ENDDO

!     --- calculate spline coefficients ---

      DO NZ=2,NZMAX
         DZ=Z(NZ)-Z(NZ-1)
         DZ1=1.D0/DZ
         DZ2=DZ1*DZ1
         DZ3=DZ2*DZ1
         DO J=1,4
            DO I=1,4
               VZ(I,J)=0.D0
            ENDDO
         END DO
         VZ(1,1)= 1.D0
         VZ(2,2)= 1.D0
         VZ(3,1)=-3.D0*DZ2
         VZ(3,2)=-2.D0*DZ1
         VZ(3,3)= 3.D0*DZ2
         VZ(3,4)=-1.D0*DZ1
         VZ(4,1)= 2.D0*DZ3
         VZ(4,2)=      DZ2
         VZ(4,3)=-2.D0*DZ3
         VZ(4,4)=      DZ2

      DO NY=2,NYMAX
         DY=Y(NY)-Y(NY-1)
         DY1=1.D0/DY
         DY2=DY1*DY1
         DY3=DY2*DY1
         DO J=1,4
            DO I=1,4
               VY(I,J)=0.D0
            ENDDO
         END DO
         VY(1,1)= 1.D0
         VY(2,2)= 1.D0
         VY(3,1)=-3.D0*DY2
         VY(3,2)=-2.D0*DY1
         VY(3,3)= 3.D0*DY2
         VY(3,4)=-1.D0*DY1
         VY(4,1)= 2.D0*DY3
         VY(4,2)=      DY2
         VY(4,3)=-2.D0*DY3
         VY(4,4)=      DY2

      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)
         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         DO J=1,4
            DO I=1,4
               VX(I,J)=0.D0
            ENDDO
         END DO
         VX(1,1)= 1.D0
         VX(2,2)= 1.D0
         VX(3,1)=-3.D0*DX2
         VX(3,2)=-2.D0*DX1
         VX(3,3)= 3.D0*DX2
         VX(3,4)=-1.D0*DX1
         VX(4,1)= 2.D0*DX3
         VX(4,2)=      DX2
         VX(4,3)=-2.D0*DX3
         VX(4,4)=      DX2

         T(1,1,1)=   F(NX-1,NY-1,NZ-1)
         T(1,1,2)=  FZ(NX-1,NY-1,NZ-1)
         T(1,1,3)=   F(NX-1,NY-1,NZ)
         T(1,1,4)=  FZ(NX-1,NY-1,NZ)

         T(1,2,1)=  FY(NX-1,NY-1,NZ-1)
         T(1,2,2)= FYZ(NX-1,NY-1,NZ-1)
         T(1,2,3)=  FY(NX-1,NY-1,NZ)
         T(1,2,4)= FYZ(NX-1,NY-1,NZ)

         T(1,3,1)=   F(NX-1,NY  ,NZ-1)
         T(1,3,2)=  FZ(NX-1,NY  ,NZ-1)
         T(1,3,3)=   F(NX-1,NY  ,NZ)
         T(1,3,4)=  FZ(NX-1,NY  ,NZ)

         T(1,4,1)=  FY(NX-1,NY  ,NZ-1)
         T(1,4,2)= FYZ(NX-1,NY  ,NZ-1)
         T(1,4,3)=  FY(NX-1,NY  ,NZ)
         T(1,4,4)= FYZ(NX-1,NY  ,NZ)

         T(2,1,1)=  FX(NX-1,NY-1,NZ-1)
         T(2,1,2)= FZX(NX-1,NY-1,NZ-1)
         T(2,1,3)=  FX(NX-1,NY-1,NZ)
         T(2,1,4)= FZX(NX-1,NY-1,NZ)

         T(2,2,1)= FXY(NX-1,NY-1,NZ-1)
         T(2,2,2)=FXYZ(NX-1,NY-1,NZ-1)
         T(2,2,3)= FXY(NX-1,NY-1,NZ)
         T(2,2,4)=FXYZ(NX-1,NY-1,NZ)

         T(2,3,1)=  FX(NX-1,NY  ,NZ-1)
         T(2,3,2)= FZX(NX-1,NY  ,NZ-1)
         T(2,3,3)=  FX(NX-1,NY  ,NZ)
         T(2,3,4)= FZX(NX-1,NY  ,NZ)

         T(2,4,1)= FXY(NX-1,NY  ,NZ-1)
         T(2,4,2)=FXYZ(NX-1,NY  ,NZ-1)
         T(2,4,3)= FXY(NX-1,NY  ,NZ)
         T(2,4,4)=FXYZ(NX-1,NY  ,NZ)

         T(3,1,1)=   F(NX  ,NY-1,NZ-1)
         T(3,1,2)=  FZ(NX  ,NY-1,NZ-1)
         T(3,1,3)=   F(NX  ,NY-1,NZ)
         T(3,1,4)=  FZ(NX  ,NY-1,NZ)

         T(3,2,1)=  FY(NX  ,NY-1,NZ-1)
         T(3,2,2)= FYZ(NX  ,NY-1,NZ-1)
         T(3,2,3)=  FY(NX  ,NY-1,NZ)
         T(3,2,4)= FYZ(NX  ,NY-1,NZ)

         T(3,3,1)=   F(NX  ,NY  ,NZ-1)
         T(3,3,2)=  FZ(NX  ,NY  ,NZ-1)
         T(3,3,3)=   F(NX  ,NY  ,NZ)
         T(3,3,4)=  FZ(NX  ,NY  ,NZ)

         T(3,4,1)=  FY(NX  ,NY  ,NZ-1)
         T(3,4,2)= FYZ(NX  ,NY  ,NZ-1)
         T(3,4,3)=  FY(NX  ,NY  ,NZ)
         T(3,4,4)= FYZ(NX  ,NY  ,NZ)

         T(4,1,1)=  FX(NX  ,NY-1,NZ-1)
         T(4,1,2)= FZX(NX  ,NY-1,NZ-1)
         T(4,1,3)=  FX(NX  ,NY-1,NZ)
         T(4,1,4)= FZX(NX  ,NY-1,NZ)

         T(4,2,1)= FXY(NX  ,NY-1,NZ-1)
         T(4,2,2)=FXYZ(NX  ,NY-1,NZ-1)
         T(4,2,3)= FXY(NX  ,NY-1,NZ)
         T(4,2,4)=FXYZ(NX  ,NY-1,NZ)

         T(4,3,1)=  FX(NX  ,NY  ,NZ-1)
         T(4,3,2)= FZX(NX  ,NY  ,NZ-1)
         T(4,3,3)=  FX(NX  ,NY  ,NZ)
         T(4,3,4)= FZX(NX  ,NY  ,NZ)

         T(4,4,1)= FXY(NX  ,NY  ,NZ-1)
         T(4,4,2)=FXYZ(NX  ,NY  ,NZ-1)
         T(4,4,3)= FXY(NX  ,NY  ,NZ)
         T(4,4,4)=FXYZ(NX  ,NY  ,NZ)

         DO K=1,4
            DO J=1,4
               DO I=1,4
                  TT(I,J,K)=0.D0
                  DO L=1,4
                     TT(I,J,K)=TT(I,J,K)+VX(I,L)*T(L,J,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K=1,4
            DO J=1,4
               DO I=1,4
                  TTT(I,J,K)=0.D0
                  DO L=1,4
                     TTT(I,J,K)=TTT(I,J,K)+VY(J,L)*TT(I,L,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K=1,4
            DO J=1,4
               DO I=1,4
                  TTTT(I,J,K)=0.D0
                  DO L=1,4
                     TTTT(I,J,K)=TTTT(I,J,K)+VZ(K,L)*TTT(I,J,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K=1,4
            DO J=1,4
               DO I=1,4
                  U(I,J,K,NX,NY,NZ)=TTTT(I,J,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      IERR=0
      RETURN

! 9001 WRITE(6,*) 'XX SPL3D: NXMAX.GT.NMAX:',NXMAX,NMAX
!      IERR=1
!      RETURN
! 9002 WRITE(6,*) 'XX SPL3D: NYMAX.GT.NMAX:',NYMAX,NMAX
!      IERR=2
!      RETURN
 9003 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=3
      RETURN
 9004 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=4
      RETURN
 9005 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=5
      RETURN
 9006 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=6
      RETURN
 9007 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=7
      RETURN
 9008 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=8
      RETURN
      END SUBROUTINE SPL3D

!     ****** Three-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL3DF(X0,Y0,Z0,F0,X,Y,Z,U,NXM,NYM,NXMAX,NYMAX,NZMAX,IERR)

      IMPLICIT NONE

      REAL(8),                          INTENT(IN) :: X0, Y0, Z0
      REAL(8),                          INTENT(OUT):: F0
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NZMAX),        INTENT(IN) :: Z
      REAL(8), DIMENSION(4,4,4,NXM,NYM,NZMAX),INTENT(IN) :: U
      INTEGER(4),                       INTENT(IN) :: NXM,NYM,NXMAX,NYMAX,NZMAX
      INTEGER(4),                       INTENT(OUT):: IERR

      INTEGER(4)            :: NX, NY, NZ, I, J
      REAL(8)               :: DX, DY, DZ, FX, FY, FZ
      REAL(8),DIMENSION(4)  :: UF
      REAL(8),DIMENSION(4,4)  :: UFF

      IERR=0

      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      IF(Z(NZMAX).EQ.Z(1)) THEN
         IERR=7
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

      FZ=1.D0/(Z(NZMAX)-Z(1))
      NZ=NINT((Z0-Z(1))*FZ*(NZMAX-1))+1
      IF(NZ.LT.1) THEN
         IERR=5
         NZ=2
      ENDIF
      IF(NZ.GT.NZMAX) THEN
         IERR=6
         NZ=NZMAX
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

 6001 IF(NY.GE.NYMAX) GOTO 6002
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 6002
         NY=NY+1
         GOTO 6001
 6002 CONTINUE
 6003 IF(NY.LE.2) GOTO 6004
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 6004
         NY=NY-1
         GOTO 6003
 6004 CONTINUE
      IF(NY.LT.2)     NY=2

 7001 IF(NZ.GE.NZMAX) GOTO 7002
      IF((Z0-Z(NZ  ))*FZ.LE.0.D0) GOTO 7002
         NZ=NZ+1
         GOTO 7001
 7002 CONTINUE
 7003 IF(NZ.LE.2) GOTO 7004
      IF((Z0-Z(NZ-1))*FZ.GE.0.D0) GOTO 7004
         NZ=NZ-1
         GOTO 7003
 7004 CONTINUE
      IF(NZ.LT.2)     NZ=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)
      DZ=Z0-Z(NZ-1)

      DO J=1,4
      DO I=1,4
         UFF(I,J)=((U(I,J,4,NX,NY,NZ) *DZ     +U(I,J,3,NX,NY,NZ))*DZ &
                   +U(I,J,2,NX,NY,NZ))*DZ     +U(I,J,1,NX,NY,NZ)
      ENDDO
      ENDDO
      DO I=1,4
         UF(I)=((UFF(I,4)*DY+UFF(I,3))*DY+UFF(I,2))*DY+UFF(I,1)
      ENDDO
      F0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      IERR=0
      RETURN
      END SUBROUTINE SPL3DF

!     ****** Three-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE SPL3DD(X0,Y0,Z0,F0,FX0,FY0,FZ0, &
                        X,Y,Z,U,NXM,NYM,NXMAX,NYMAX,NZMAX,IERR)

      IMPLICIT NONE

      REAL(8),                          INTENT(IN) :: X0,Y0,Z0
      REAL(8),                          INTENT(OUT):: F0,FX0,FY0,FZ0
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NZMAX),        INTENT(IN) :: Z
      REAL(8), DIMENSION(4,4,4,NXM,NYM,NZMAX),INTENT(IN) :: U
      INTEGER(4),                       INTENT(IN) :: NXM,NYM,NXMAX,NYMAX,NZMAX
      INTEGER(4),                       INTENT(OUT):: IERR

      INTEGER(4)  :: NX,NY,NZ,I,J
      REAL(8)     :: DX,DY,DZ,FX,FY,FZ
      REAL(8),DIMENSION(4):: UF
      REAL(8),DIMENSION(4,4):: UFF

      IERR=0

      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      IF(Z(NZMAX).EQ.Z(1)) THEN
         IERR=7
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

      FZ=1.D0/(Z(NZMAX)-Z(1))
      NZ=NINT((Z0-Z(1))*FZ*(NZMAX-1))+1
      IF(NZ.LT.1) THEN
         IERR=5
         NZ=2
      ENDIF
      IF(NZ.GT.NZMAX) THEN
         IERR=6
         NZ=NZMAX
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

 6001 IF(NY.GE.NYMAX) GOTO 6002
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 6002
         NY=NY+1
         GOTO 6001
 6002 CONTINUE
 6003 IF(NY.LE.2) GOTO 6004
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 6004
         NY=NY-1
         GOTO 6003
 6004 CONTINUE
      IF(NY.LT.2)     NY=2

 7001 IF(NZ.GE.NZMAX) GOTO 7002
      IF((Z0-Z(NZ  ))*FZ.LE.0.D0) GOTO 7002
         NZ=NZ+1
         GOTO 7001
 7002 CONTINUE
 7003 IF(NZ.LE.2) GOTO 7004
      IF((Z0-Z(NZ-1))*FZ.GE.0.D0) GOTO 7004
         NZ=NZ-1
         GOTO 7003
 7004 CONTINUE
      IF(NZ.LT.2)     NZ=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)
      DZ=Z0-Z(NZ-1)

      DO J=1,4
      DO I=1,4
         UFF(I,J)=((U(I,J,4,NX,NY,NZ) *DZ     +U(I,J,3,NX,NY,NZ))*DZ &
                   +U(I,J,2,NX,NY,NZ))*DZ     +U(I,J,1,NX,NY,NZ)
      ENDDO
      ENDDO
      DO I=1,4
         UF(I)=((UFF(I,4)*DY+UFF(I,3))*DY+UFF(I,2))*DY+UFF(I,1)
      ENDDO
      F0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)
      FX0= (3*UF(4)*DX+2*UF(3))*DX+UF(2)

      DO I=1,4
         UF(I)=(3*UFF(I,4)*DY+2*UFF(I,3))*DY+UFF(I,2)
      ENDDO
      FY0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      DO J=1,4
      DO I=1,4
         UFF(I,J)=(3*U(I,J,4,NX,NY,NZ)*DZ+2*U(I,J,3,NX,NY,NZ))*DZ &
                    +U(I,J,2,NX,NY,NZ)
      ENDDO
      ENDDO
      DO I=1,4
         UF(I)=((UFF(I,4)*DY+UFF(I,3))*DY+UFF(I,2))*DY+UFF(I,1)
      ENDDO
      FZ0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      IERR=0
      RETURN
      END SUBROUTINE SPL3DD
!
!     ****** Three-Dimensional Spline Interpolation ******
!       **** Calculation of coefficients ****

      SUBROUTINE CSPL3D(X,Y,Z,F,FX,FY,FZ,FXY,FYZ,FZX,FXYZ,U, &
                       NXM,NYM,NXMAX,NYMAX,NZMAX,IDX,IDY,IDZ,IERR)
!      INPUT : X(NXMAX)        : COORDINATES
!              Y(NYMAX)        : COORDINATES
!              Z(NZMAX)        : COORDINATES
!              F(NXM,NYM,NZMAX)   : VALUE
!              FX(NXM,NYM,NZMAX)  : EDGE DERIVATIVES FOR IDX != 0
!              FY(NXM,NYM,NZMAX)  : EDGE DERIVATIVES FOR IDY != 0
!              FZ(NXM,NYM,NZMAX)  : EDGE DERIVATIVES FOR IDZ != 0
!              FXY(NXM,NYM,NZMAX) : CORNER DERIVATIVES FOR IDX OR IDY != 0
!              FYZ(NXM,NYM,NZMAX) : CORNER DERIVATIVES FOR IDY OR IDZ != 0
!              FZX(NXM,NYM,NZMAX) : CORNER DERIVATIVES FOR IDZ OR IDX != 0
!              FXYZ(NXM,NYM,NZMAX): CORNER DERIVS FOR IDX OR IDY OR IDZ != 0
!              NXM       : ARRAY SIZE
!              NYM       : ARRAY SIZE
!              NXMAX     : NUMBER OF VARIABLES
!              NYMAX     : NUMBER OF VARIABLES
!              NZMAX     : NUMBER OF VARIABLES
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
!              IDZ       : 0 : SECOND DERIVATIVES = 0 AT Z(1) AND Z(NZMAX)
!                          1 : DERIVATIVE FZ(1) IS GIVEN
!                          2 : DERIVATIVE FZ(NYMAX) IS GIVEN
!                          3 : DERIVATIVES FZ(1) AXD FZ(NYMAX) ARE GIVEN
!                          4 : PERIODIC FUNCTION: F(1)=F(NZMAX)
!      OUTPUT: U(4,4,4,NXM,NYM,NZMAX): SPLINE COEFICIENTS
!              FX(NXM,NYM,NZMAX) : ESTIMATED DERIVATIVES
!              FY(NXM,NYM,NZMAX) : ESTIMATED DERIVATIVES
!              FZ(NXM,NYM,NZMAX) : ESTIMATED DERIVATIVES
!              FXY(NXM,NYM,NZMAX): ESTIMATED DERIVATIVES
!              FYZ(NXM,NYM,NZMAX): ESTIMATED DERIVATIVES
!              FZX(NXM,NYM,NZMAX): ESTIMATED DERIVATIVES
!              IERR          : ERROR INDICATOR

      IMPLICIT NONE

      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NZMAX),        INTENT(IN) :: Z
      COMPLEX(8), DIMENSION(NXM,NYM,NZMAX),INTENT(IN) :: F
      COMPLEX(8), DIMENSION(NXM,NYM,NZMAX),INTENT(INOUT):: FX,FY,FZ
      COMPLEX(8), DIMENSION(NXM,NYM,NZMAX),INTENT(INOUT):: FXY,FYZ,FZX,FXYZ
      COMPLEX(8), DIMENSION(4,4,4,NXM,NYM,NZMAX),INTENT(OUT):: U
      INTEGER(4),INTENT(IN) :: NXM,NYM,NXMAX,NYMAX,NZMAX,IDX,IDY,IDZ
      INTEGER(4),INTENT(OUT):: IERR

      INTEGER(4) :: NX,NY,NZ,I,J,K,L
      REAL(8)    :: DX,DXM,DXP,DX1,DX2,DX3
      REAL(8)    :: DY,DYM,DYP,DY1,DY2,DY3
      REAL(8)    :: DZ,DZM,DZP,DZ1,DZ2,DZ3
      INTEGER(4)    :: IDX1,IDX2,IDY1,IDY2,IDZ1,IDZ2
      COMPLEX(8),DIMENSION(4,4):: VX,VY,VZ
      COMPLEX(8),DIMENSION(4,4,4):: T,TT,TTT,TTTT
      COMPLEX(8),    DIMENSION(4,NXMAX):: UX,UX0
      COMPLEX(8),    DIMENSION(4,NYMAX):: UY,UY0
      COMPLEX(8),    DIMENSION(4,NZMAX):: UZ,UZ0
      COMPLEX(8),    DIMENSION(NXMAX)  :: BX
      COMPLEX(8),    DIMENSION(NYMAX)  :: BY
      COMPLEX(8),    DIMENSION(NZMAX)  :: BZ

      IDX1=MOD(IDX,2)
      IDX2=MOD(IDX/2,2)
      IDY1=MOD(IDY,2)
      IDY2=MOD(IDY/2,2)
      IDZ1=MOD(IDZ,2)
      IDZ2=MOD(IDZ/2,2)

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

!     --- calculate UZ ---

      IF(IDZ.EQ.4) THEN
         DZM=Z(NZMAX)-Z(NZMAX-1)
         DZP=Z(2)-Z(1)
         UZ(1,1)=DZP
         UZ(2,1)=2.D0*(DZP+DZM)
         UZ(3,1)=DZM
      ELSE
         IF(IDZ1.EQ.0) THEN
            DZP=Z(2)-Z(1)
            UZ(1,1)=0.D0
            UZ(2,1)=2.D0*DZP
            UZ(3,1)=DZP
         ELSE
            UZ(1,1)=0.D0
            UZ(2,1)=1.D0
            UZ(3,1)=0.D0
         ENDIF
      ENDIF
      DO NZ=2,NZMAX-1
         DZM=Z(NZ)-Z(NZ-1)
         DZP=Z(NZ+1)-Z(NZ)
         UZ(1,NZ)=DZP
         UZ(2,NZ)=2.D0*(DZP+DZM)
         UZ(3,NZ)=DZM
      ENDDO
      IF(IDZ.EQ.4) THEN
         DZM=Z(NZMAX)-Z(NZMAX-1)
         DZP=Z(2)-Z(1)
         UZ(1,NZMAX)=DZP
         UZ(2,NZMAX)=2.D0*(DZP+DZM)
         UZ(3,NZMAX)=DZM
      ELSE
         IF(IDZ2.EQ.0) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            UZ(1,NZMAX)=DZM
            UZ(2,NZMAX)=2.D0*DZM
            UZ(3,NZMAX)=0.D0
         ELSE
            UZ(1,NZMAX)=0.D0
            UZ(2,NZMAX)=1.D0
            UZ(3,NZMAX)=0.D0
         ENDIF
      ENDIF
      DO NZ=1,NZMAX
         UZ0(1,NZ)=UZ(1,NZ)
         UZ0(2,NZ)=UZ(2,NZ)
         UZ0(3,NZ)=UZ(3,NZ)
      ENDDO

!     --- calculate FX ---

      DO NZ=1,NZMAX
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(F(    2,NY,NZ)-F(      1,NY,NZ))/DXP &
     &                 +DXP*(F(NXMAX,NY,NZ)-F(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(F(2,NY,NZ)-F(1,NY,NZ))
            ELSE
               BX(1)=FX(1,NY,NZ)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(F(NX+1,NY,NZ)-F(NX,  NY,NZ))/DXP &
     &                  +DXP*(F(NX,  NY,NZ)-F(NX-1,NY,NZ))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(F(    2,NY,NZ)-F(      1,NY,NZ))/DXP &
     &                     +DXP*(F(NXMAX,NY,NZ)-F(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(F(NXMAX,NY,NZ)-F(NXMAX-1,NY,NZ))
            ELSE
               BX(NXMAX)=FX(NXMAX,NY,NZ)
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
            FX(NX,NY,NZ)=BX(NX)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FY ---

      DO NZ=1,NZMAX
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            UY(1,NY)=UY0(1,NY)
            UY(2,NY)=UY0(2,NY)
            UY(3,NY)=UY0(3,NY)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(1)=3.D0*(DYM*(F(NX,    2,NZ)-F(NX,      1,NZ))/DYP &
     &                 +DYP*(F(NX,NYMAX,NZ)-F(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY1.EQ.0) THEN
               BY(1)=3.D0*(F(NX,2,NZ)-F(NX,1,NZ))
            ELSE
               BY(1)=FY(NX,1,NZ)
            ENDIF
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(F(NX,NY+1,NZ)-F(NX,NY  ,NZ))/DYP &
     &                  +DYP*(F(NX,NY  ,NZ)-F(NX,NY-1,NZ))/DYM)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(NYMAX)=3.D0*(DYM*(F(NX,    2,NZ)-F(NX,      1,NZ))/DYP &
     &                     +DYP*(F(NX,NYMAX,NZ)-F(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY2.EQ.0) THEN
               BY(NYMAX)=3.D0*(F(NX,NYMAX,NZ)-F(NX,NYMAX-1,NZ))
            ELSE
               BY(NYMAX)=FY(NX,NYMAX,NZ)
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
            FY(NX,NY,NZ)=BY(NY)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FZ ---

      DO NY=1,NYMAX
      DO NX=1,NXMAX
         DO NZ=1,NZMAX
            UZ(1,NZ)=UZ0(1,NZ)
            UZ(2,NZ)=UZ0(2,NZ)
            UZ(3,NZ)=UZ0(3,NZ)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(1)=3.D0*(DZM*(F(NX,NY,    2)-F(NX,NY,      1))/DZP &
     &                 +DZP*(F(NX,NY,NZMAX)-F(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ1.EQ.0) THEN
               BZ(1)=3.D0*(F(NX,NY,2)-F(NX,NY,1))
            ELSE
               BZ(1)=FZ(NX,NY,1)
            ENDIF
         ENDIF
         DO NZ=2,NZMAX-1
            DZM=Z(NZ)-Z(NZ-1)
            DZP=Z(NZ+1)-Z(NZ)
            BZ(NZ)=3.D0*(DZM*(F(NX,NY,NZ+1)-F(NX,NY,NZ  ))/DZP &
     &                  +DZP*(F(NX,NY,NZ  )-F(NX,NY,NZ-1))/DZM)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(NZMAX)=3.D0*(DZM*(F(NX,NY,    2)-F(NX,NY,      1))/DZP &
     &                     +DZP*(F(NX,NY,NZMAX)-F(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ2.EQ.0) THEN
               BZ(NZMAX)=3.D0*(F(NX,NY,NZMAX)-F(NX,NY,NZMAX-1))
            ELSE
               BZ(NZMAX)=FZ(NX,NY,NZMAX)
            ENDIF
         ENDIF

         IF(IDZ.EQ.4) THEN
            CALL TDMPCDX(UZ,BZ,NZMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9005
            BZ(NZMAX)=BZ(1)
         ELSE
            CALL TDMSCDX(UZ,BZ,NZMAX,IERR)
            IF(IERR.NE.0) GOTO 9006
         ENDIF
         DO NZ=1,NZMAX
            FZ(NX,NY,NZ)=BZ(NZ)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FXY ---

      DO NZ=1,NZMAX
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(FY(    2,NY,NZ)-FY(      1,NY,NZ))/DXP &
     &                 +DXP*(FY(NXMAX,NY,NZ)-FY(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(FY(2,NY,NZ)-FY(1,NY,NZ))
            ELSE
               BX(1)=FXY(1,NY,NZ)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FY(NX+1,NY,NZ)-FY(NX,  NY,NZ))/DXP &
     &                  +DXP*(FY(NX,  NY,NZ)-FY(NX-1,NY,NZ))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(FY(    2,NY,NZ)-FY(      1,NY,NZ))/DXP &
     &                     +DXP*(FY(NXMAX,NY,NZ)-FY(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(FY(NXMAX,NY,NZ)-FY(NXMAX-1,NY,NZ))
            ELSE
               BX(NXMAX)=FXY(NXMAX,NY,NZ)
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
            FXY(NX,NY,NZ)=BX(NX)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FYZ ---

      DO NZ=1,NZMAX
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            UY(1,NY)=UY0(1,NY)
            UY(2,NY)=UY0(2,NY)
            UY(3,NY)=UY0(3,NY)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(1)=3.D0*(DYM*(FZ(NX,    2,NZ)-FZ(NX,      1,NZ))/DYP &
     &                 +DYP*(FZ(NX,NYMAX,NZ)-FZ(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY1.EQ.0) THEN
               BY(1)=3.D0*(FZ(NX,2,NZ)-FZ(NX,1,NZ))
            ELSE
               BY(1)=FYZ(NX,1,NZ)
            ENDIF
         ENDIF
         DO NY=2,NYMAX-1
            DYM=Y(NY)-Y(NY-1)
            DYP=Y(NY+1)-Y(NY)
            BY(NY)=3.D0*(DYM*(FZ(NX,NY+1,NZ)-FZ(NX,NY,  NZ))/DYP &
     &                  +DYP*(FZ(NX,NY,  NZ)-FZ(NX,NY-1,NZ))/DYM)
         ENDDO
         IF(IDY.EQ.4) THEN
            DYM=Y(NYMAX)-Y(NYMAX-1)
            DYP=Y(2)-Y(1)
            BY(NYMAX)=3.D0*(DYM*(FZ(NX,    2,NZ)-FZ(NX,      1,NZ))/DYP &
     &                     +DYP*(FZ(NX,NYMAX,NZ)-FZ(NX,NYMAX-1,NZ))/DYM)
         ELSE
            IF(IDY2.EQ.0) THEN
               BY(NYMAX)=3.D0*(FZ(NX,NYMAX,NZ)-FZ(NX,NYMAX-1,NZ))
            ELSE
               BY(NYMAX)=FYZ(NX,NYMAX,NZ)
            ENDIF
         ENDIF

         IF(IDY.EQ.4) THEN
            CALL TDMPCDX(UY,BY,NYMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9007
            BY(NYMAX)=BY(1)
         ELSE
            CALL TDMSCDX(UY,BY,NYMAX,IERR)
            IF(IERR.NE.0) GOTO 9008
         ENDIF
         DO NY=1,NYMAX
            FYZ(NX,NY,NZ)=BY(NY)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FZX ---

      DO NY=1,NYMAX
      DO NX=1,NXMAX
         DO NZ=1,NZMAX
            UZ(1,NZ)=UZ0(1,NZ)
            UZ(2,NZ)=UZ0(2,NZ)
            UZ(3,NZ)=UZ0(3,NZ)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(1)=3.D0*(DZM*(FX(NX,NY,    2)-FX(NX,NY,      1))/DZP &
     &                 +DZP*(FX(NX,NY,NZMAX)-FX(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ1.EQ.0) THEN
               BZ(1)=3.D0*(FX(NX,NY,2)-FX(NX,NY,1))
            ELSE
               BZ(1)=FZX(NX,NY,1)
            ENDIF
         ENDIF
         DO NZ=2,NZMAX-1
            DZM=Z(NZ)-Z(NZ-1)
            DZP=Z(NZ+1)-Z(NZ)
            BZ(NZ)=3.D0*(DZM*(FX(NX,NY,NZ+1)-FX(NX,NY,NZ  ))/DZP &
     &                  +DZP*(FX(NX,NY,NZ  )-FX(NX,NY,NZ-1))/DZM)
         ENDDO
         IF(IDZ.EQ.4) THEN
            DZM=Z(NZMAX)-Z(NZMAX-1)
            DZP=Z(2)-Z(1)
            BZ(NZMAX)=3.D0*(DZM*(FX(NX,NY,    2)-FX(NX,NY,      1))/DZP &
     &                     +DZP*(FX(NX,NY,NZMAX)-FX(NX,NY,NZMAX-1))/DZM)
         ELSE
            IF(IDZ2.EQ.0) THEN
               BZ(NZMAX)=3.D0*(FX(NX,NY,NZMAX)-FX(NX,NY,NZMAX-1))
            ELSE
               BZ(NZMAX)=FZX(NX,NY,NZMAX)
            ENDIF
         ENDIF

         IF(IDZ.EQ.4) THEN
            CALL TDMPCDX(UZ,BZ,NZMAX-1,IERR)
            IF(IERR.NE.0) GOTO 9007
            BZ(NZMAX)=BZ(1)
         ELSE
            CALL TDMSCDX(UZ,BZ,NZMAX,IERR)
            IF(IERR.NE.0) GOTO 9008
         ENDIF
         DO NZ=1,NZMAX
            FZX(NX,NY,NZ)=BZ(NZ)
         ENDDO
      ENDDO
      ENDDO

!     --- calculate FXYZ ---
      write(6,*) '--- calculate FXYZ ---'

      DO NZ=1,NZMAX
      DO NY=1,NYMAX
         DO NX=1,NXMAX
            UX(1,NX)=UX0(1,NX)
            UX(2,NX)=UX0(2,NX)
            UX(3,NX)=UX0(3,NX)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(1)=3.D0*(DXM*(FYZ(    2,NY,NZ)-FYZ(      1,NY,NZ))/DXP &
     &                 +DXP*(FYZ(NXMAX,NY,NZ)-FYZ(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX1.EQ.0) THEN
               BX(1)=3.D0*(FYZ(2,NY,NZ)-FYZ(1,NY,NZ))
            ELSE
               BX(1)=FXYZ(1,NY,NZ)
            ENDIF
         ENDIF
         DO NX=2,NXMAX-1
            DXM=X(NX)-X(NX-1)
            DXP=X(NX+1)-X(NX)
            BX(NX)=3.D0*(DXM*(FYZ(NX+1,NY,NZ)-FYZ(NX,  NY,NZ))/DXP &
     &                  +DXP*(FYZ(NX,  NY,NZ)-FYZ(NX-1,NY,NZ))/DXM)
         ENDDO
         IF(IDX.EQ.4) THEN
            DXM=X(NXMAX)-X(NXMAX-1)
            DXP=X(2)-X(1)
            BX(NXMAX)=3.D0*(DXM*(FYZ(    2,NY,NZ)-FYZ(      1,NY,NZ))/DXP &
     &                     +DXP*(FYZ(NXMAX,NY,NZ)-FYZ(NXMAX-1,NY,NZ))/DXM)
         ELSE
            IF(IDX2.EQ.0) THEN
               BX(NXMAX)=3.D0*(FYZ(NXMAX,NY,NZ)-FYZ(NXMAX-1,NY,NZ))
            ELSE
               BX(NXMAX)=FXYZ(NXMAX,NY,NZ)
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
            FXYZ(NX,NY,NZ)=BX(NX)
         ENDDO
      ENDDO
      ENDDO

!      DO NX=1,NXMAX
!      DO NY=1,NYMAX
!      DO NZ=1,NZMAX
!         WRITE(6,'(3I3,1P7E10.2)') NX,NY,NZ, &
!              F(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ), &
!              FXY(NX,NY,NZ),FYZ(NX,NY,NZ),FZX(NX,NY,NZ)
!      ENDDO
!      ENDDO
!      ENDDO

!     --- calculate spline coefficients ---

      DO NZ=2,NZMAX
         DZ=Z(NZ)-Z(NZ-1)
         DZ1=1.D0/DZ
         DZ2=DZ1*DZ1
         DZ3=DZ2*DZ1
         DO J=1,4
            DO I=1,4
               VZ(I,J)=0.D0
            ENDDO
         END DO
         VZ(1,1)= 1.D0
         VZ(2,2)= 1.D0
         VZ(3,1)=-3.D0*DZ2
         VZ(3,2)=-2.D0*DZ1
         VZ(3,3)= 3.D0*DZ2
         VZ(3,4)=-1.D0*DZ1
         VZ(4,1)= 2.D0*DZ3
         VZ(4,2)=      DZ2
         VZ(4,3)=-2.D0*DZ3
         VZ(4,4)=      DZ2

      DO NY=2,NYMAX
         DY=Y(NY)-Y(NY-1)
         DY1=1.D0/DY
         DY2=DY1*DY1
         DY3=DY2*DY1
         DO J=1,4
            DO I=1,4
               VY(I,J)=0.D0
            ENDDO
         END DO
         VY(1,1)= 1.D0
         VY(2,2)= 1.D0
         VY(3,1)=-3.D0*DY2
         VY(3,2)=-2.D0*DY1
         VY(3,3)= 3.D0*DY2
         VY(3,4)=-1.D0*DY1
         VY(4,1)= 2.D0*DY3
         VY(4,2)=      DY2
         VY(4,3)=-2.D0*DY3
         VY(4,4)=      DY2

      DO NX=2,NXMAX
         DX=X(NX)-X(NX-1)
         DX1=1.D0/DX
         DX2=DX1*DX1
         DX3=DX2*DX1
         DO J=1,4
            DO I=1,4
               VX(I,J)=0.D0
            ENDDO
         END DO
         VX(1,1)= 1.D0
         VX(2,2)= 1.D0
         VX(3,1)=-3.D0*DX2
         VX(3,2)=-2.D0*DX1
         VX(3,3)= 3.D0*DX2
         VX(3,4)=-1.D0*DX1
         VX(4,1)= 2.D0*DX3
         VX(4,2)=      DX2
         VX(4,3)=-2.D0*DX3
         VX(4,4)=      DX2

         T(1,1,1)=   F(NX-1,NY-1,NZ-1)
         T(1,1,2)=  FZ(NX-1,NY-1,NZ-1)
         T(1,1,3)=   F(NX-1,NY-1,NZ)
         T(1,1,4)=  FZ(NX-1,NY-1,NZ)

         T(1,2,1)=  FY(NX-1,NY-1,NZ-1)
         T(1,2,2)= FYZ(NX-1,NY-1,NZ-1)
         T(1,2,3)=  FY(NX-1,NY-1,NZ)
         T(1,2,4)= FYZ(NX-1,NY-1,NZ)

         T(1,3,1)=   F(NX-1,NY  ,NZ-1)
         T(1,3,2)=  FZ(NX-1,NY  ,NZ-1)
         T(1,3,3)=   F(NX-1,NY  ,NZ)
         T(1,3,4)=  FZ(NX-1,NY  ,NZ)

         T(1,4,1)=  FY(NX-1,NY  ,NZ-1)
         T(1,4,2)= FYZ(NX-1,NY  ,NZ-1)
         T(1,4,3)=  FY(NX-1,NY  ,NZ)
         T(1,4,4)= FYZ(NX-1,NY  ,NZ)

         T(2,1,1)=  FX(NX-1,NY-1,NZ-1)
         T(2,1,2)= FZX(NX-1,NY-1,NZ-1)
         T(2,1,3)=  FX(NX-1,NY-1,NZ)
         T(2,1,4)= FZX(NX-1,NY-1,NZ)

         T(2,2,1)= FXY(NX-1,NY-1,NZ-1)
         T(2,2,2)=FXYZ(NX-1,NY-1,NZ-1)
         T(2,2,3)= FXY(NX-1,NY-1,NZ)
         T(2,2,4)=FXYZ(NX-1,NY-1,NZ)

         T(2,3,1)=  FX(NX-1,NY  ,NZ-1)
         T(2,3,2)= FZX(NX-1,NY  ,NZ-1)
         T(2,3,3)=  FX(NX-1,NY  ,NZ)
         T(2,3,4)= FZX(NX-1,NY  ,NZ)

         T(2,4,1)= FXY(NX-1,NY  ,NZ-1)
         T(2,4,2)=FXYZ(NX-1,NY  ,NZ-1)
         T(2,4,3)= FXY(NX-1,NY  ,NZ)
         T(2,4,4)=FXYZ(NX-1,NY  ,NZ)

         T(3,1,1)=   F(NX  ,NY-1,NZ-1)
         T(3,1,2)=  FZ(NX  ,NY-1,NZ-1)
         T(3,1,3)=   F(NX  ,NY-1,NZ)
         T(3,1,4)=  FZ(NX  ,NY-1,NZ)

         T(3,2,1)=  FY(NX  ,NY-1,NZ-1)
         T(3,2,2)= FYZ(NX  ,NY-1,NZ-1)
         T(3,2,3)=  FY(NX  ,NY-1,NZ)
         T(3,2,4)= FYZ(NX  ,NY-1,NZ)

         T(3,3,1)=   F(NX  ,NY  ,NZ-1)
         T(3,3,2)=  FZ(NX  ,NY  ,NZ-1)
         T(3,3,3)=   F(NX  ,NY  ,NZ)
         T(3,3,4)=  FZ(NX  ,NY  ,NZ)

         T(3,4,1)=  FY(NX  ,NY  ,NZ-1)
         T(3,4,2)= FYZ(NX  ,NY  ,NZ-1)
         T(3,4,3)=  FY(NX  ,NY  ,NZ)
         T(3,4,4)= FYZ(NX  ,NY  ,NZ)

         T(4,1,1)=  FX(NX  ,NY-1,NZ-1)
         T(4,1,2)= FZX(NX  ,NY-1,NZ-1)
         T(4,1,3)=  FX(NX  ,NY-1,NZ)
         T(4,1,4)= FZX(NX  ,NY-1,NZ)

         T(4,2,1)= FXY(NX  ,NY-1,NZ-1)
         T(4,2,2)=FXYZ(NX  ,NY-1,NZ-1)
         T(4,2,3)= FXY(NX  ,NY-1,NZ)
         T(4,2,4)=FXYZ(NX  ,NY-1,NZ)

         T(4,3,1)=  FX(NX  ,NY  ,NZ-1)
         T(4,3,2)= FZX(NX  ,NY  ,NZ-1)
         T(4,3,3)=  FX(NX  ,NY  ,NZ)
         T(4,3,4)= FZX(NX  ,NY  ,NZ)

         T(4,4,1)= FXY(NX  ,NY  ,NZ-1)
         T(4,4,2)=FXYZ(NX  ,NY  ,NZ-1)
         T(4,4,3)= FXY(NX  ,NY  ,NZ)
         T(4,4,4)=FXYZ(NX  ,NY  ,NZ)

         DO K=1,4
            DO J=1,4
               DO I=1,4
                  TT(I,J,K)=0.D0
                  DO L=1,4
                     TT(I,J,K)=TT(I,J,K)+VX(I,L)*T(L,J,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K=1,4
            DO J=1,4
               DO I=1,4
                  TTT(I,J,K)=0.D0
                  DO L=1,4
                     TTT(I,J,K)=TTT(I,J,K)+VY(J,L)*TT(I,L,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K=1,4
            DO J=1,4
               DO I=1,4
                  TTTT(I,J,K)=0.D0
                  DO L=1,4
                     TTTT(I,J,K)=TTTT(I,J,K)+VZ(K,L)*TTT(I,J,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K=1,4
            DO J=1,4
               DO I=1,4
                  U(I,J,K,NX,NY,NZ)=TTTT(I,J,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      IERR=0
      RETURN

! 9001 WRITE(6,*) 'XX SPL3D: NXMAX.GT.NMAX:',NXMAX,NMAX
!      IERR=1
!      RETURN
! 9002 WRITE(6,*) 'XX SPL3D: NYMAX.GT.NMAX:',NYMAX,NMAX
!      IERR=2
!      RETURN
 9003 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=3
      RETURN
 9004 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=4
      RETURN
 9005 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=5
      RETURN
 9006 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=6
      RETURN
 9007 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=7
      RETURN
 9008 WRITE(6,*) 'XX SPL3D: TDMSRDX ERROR: IERR=',IERR
      IERR=8
      RETURN
      END SUBROUTINE CSPL3D

!     ****** Three-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL3DF(X0,Y0,Z0,F0,X,Y,Z,U,NXM,NYM,NXMAX,NYMAX,NZMAX,IERR)

      IMPLICIT NONE

      REAL(8),                          INTENT(IN) :: X0, Y0, Z0
      COMPLEX(8),                          INTENT(OUT):: F0
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NZMAX),        INTENT(IN) :: Z
      COMPLEX(8), DIMENSION(4,4,4,NXM,NYM,NZMAX),INTENT(IN) :: U
      INTEGER(4),                       INTENT(IN) :: NXM,NYM,NXMAX,NYMAX,NZMAX
      INTEGER(4),                       INTENT(OUT):: IERR

      INTEGER(4)            :: NX, NY, NZ, I, J
      REAL(8)               :: DX, DY, DZ, FX, FY, FZ
      COMPLEX(8),DIMENSION(4)  :: UF
      COMPLEX(8),DIMENSION(4,4)  :: UFF

      IERR=0

      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      IF(Z(NZMAX).EQ.Z(1)) THEN
         IERR=7
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

      FZ=1.D0/(Z(NZMAX)-Z(1))
      NZ=NINT((Z0-Z(1))*FZ*(NZMAX-1))+1
      IF(NZ.LT.1) THEN
         IERR=5
         NZ=2
      ENDIF
      IF(NZ.GT.NZMAX) THEN
         IERR=6
         NZ=NZMAX
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

 6001 IF(NY.GE.NYMAX) GOTO 6002
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 6002
         NY=NY+1
         GOTO 6001
 6002 CONTINUE
 6003 IF(NY.LE.2) GOTO 6004
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 6004
         NY=NY-1
         GOTO 6003
 6004 CONTINUE
      IF(NY.LT.2)     NY=2

 7001 IF(NZ.GE.NZMAX) GOTO 7002
      IF((Z0-Z(NZ  ))*FZ.LE.0.D0) GOTO 7002
         NZ=NZ+1
         GOTO 7001
 7002 CONTINUE
 7003 IF(NZ.LE.2) GOTO 7004
      IF((Z0-Z(NZ-1))*FZ.GE.0.D0) GOTO 7004
         NZ=NZ-1
         GOTO 7003
 7004 CONTINUE
      IF(NZ.LT.2)     NZ=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)
      DZ=Z0-Z(NZ-1)

      DO J=1,4
      DO I=1,4
         UFF(I,J)=((U(I,J,4,NX,NY,NZ) *DZ     +U(I,J,3,NX,NY,NZ))*DZ &
                   +U(I,J,2,NX,NY,NZ))*DZ     +U(I,J,1,NX,NY,NZ)
      ENDDO
      ENDDO
      DO I=1,4
         UF(I)=((UFF(I,4)*DY+UFF(I,3))*DY+UFF(I,2))*DY+UFF(I,1)
      ENDDO
      F0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      IERR=0
      RETURN
      END SUBROUTINE CSPL3DF

!     ****** Three-Dimensional Spline Interpolation ******
!       **** Calculation of interpolation  ****

      SUBROUTINE CSPL3DD(X0,Y0,Z0,F0,FX0,FY0,FZ0, &
                        X,Y,Z,U,NXM,NYM,NXMAX,NYMAX,NZMAX,IERR)

      IMPLICIT NONE

      REAL(8),                          INTENT(IN) :: X0,Y0,Z0
      COMPLEX(8),                          INTENT(OUT):: F0,FX0,FY0,FZ0
      REAL(8), DIMENSION(NXMAX),        INTENT(IN) :: X
      REAL(8), DIMENSION(NYMAX),        INTENT(IN) :: Y
      REAL(8), DIMENSION(NZMAX),        INTENT(IN) :: Z
      COMPLEX(8), DIMENSION(4,4,4,NXM,NYM,NZMAX),INTENT(IN) :: U
      INTEGER(4),                       INTENT(IN) :: NXM,NYM,NXMAX,NYMAX,NZMAX
      INTEGER(4),                       INTENT(OUT):: IERR

      INTEGER(4)  :: NX,NY,NZ,I,J
      REAL(8)     :: DX,DY,DZ,FX,FY,FZ
      COMPLEX(8),DIMENSION(4):: UF
      COMPLEX(8),DIMENSION(4,4):: UFF

      IERR=0

      IF(X(NXMAX).EQ.X(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(Y(NYMAX).EQ.Y(1)) THEN
         IERR=8
         RETURN
      ENDIF
      IF(Z(NZMAX).EQ.Z(1)) THEN
         IERR=7
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

      FZ=1.D0/(Z(NZMAX)-Z(1))
      NZ=NINT((Z0-Z(1))*FZ*(NZMAX-1))+1
      IF(NZ.LT.1) THEN
         IERR=5
         NZ=2
      ENDIF
      IF(NZ.GT.NZMAX) THEN
         IERR=6
         NZ=NZMAX
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

 6001 IF(NY.GE.NYMAX) GOTO 6002
      IF((Y0-Y(NY  ))*FY.LE.0.D0) GOTO 6002
         NY=NY+1
         GOTO 6001
 6002 CONTINUE
 6003 IF(NY.LE.2) GOTO 6004
      IF((Y0-Y(NY-1))*FY.GE.0.D0) GOTO 6004
         NY=NY-1
         GOTO 6003
 6004 CONTINUE
      IF(NY.LT.2)     NY=2

 7001 IF(NZ.GE.NZMAX) GOTO 7002
      IF((Z0-Z(NZ  ))*FZ.LE.0.D0) GOTO 7002
         NZ=NZ+1
         GOTO 7001
 7002 CONTINUE
 7003 IF(NZ.LE.2) GOTO 7004
      IF((Z0-Z(NZ-1))*FZ.GE.0.D0) GOTO 7004
         NZ=NZ-1
         GOTO 7003
 7004 CONTINUE
      IF(NZ.LT.2)     NZ=2

      DX=X0-X(NX-1)
      DY=Y0-Y(NY-1)
      DZ=Z0-Z(NZ-1)

      DO J=1,4
      DO I=1,4
         UFF(I,J)=((U(I,J,4,NX,NY,NZ) *DZ     +U(I,J,3,NX,NY,NZ))*DZ &
                   +U(I,J,2,NX,NY,NZ))*DZ     +U(I,J,1,NX,NY,NZ)
      ENDDO
      ENDDO
      DO I=1,4
         UF(I)=((UFF(I,4)*DY+UFF(I,3))*DY+UFF(I,2))*DY+UFF(I,1)
      ENDDO
      F0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)
      FX0= (3*UF(4)*DX+2*UF(3))*DX+UF(2)

      DO I=1,4
         UF(I)=(3*UFF(I,4)*DY+2*UFF(I,3))*DY+UFF(I,2)
      ENDDO
      FY0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      DO J=1,4
      DO I=1,4
         UFF(I,J)=(3*U(I,J,4,NX,NY,NZ)*DZ+2*U(I,J,3,NX,NY,NZ))*DZ &
                    +U(I,J,2,NX,NY,NZ)
      ENDDO
      ENDDO
      DO I=1,4
         UF(I)=((UFF(I,4)*DY+UFF(I,3))*DY+UFF(I,2))*DY+UFF(I,1)
      ENDDO
      FZ0= ((UF(4)*DX+UF(3))*DX+UF(2))*DX+UF(1)

      IERR=0
      RETURN
      END SUBROUTINE CSPL3DD
