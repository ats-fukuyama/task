C     $Id$
C
      IMPLICIT NONE
      INTEGER,PARAMETER:: NXM=101,NYM=101,NZM=101,NGM=1001
      REAL(8),DIMENSION(NXM)::  X0
      REAL(8),DIMENSION(NYM)::  Y0
      REAL(8),DIMENSION(NZM)::  Z0
      REAL(8),DIMENSION(NXM,NYM,NZM):: F0,FX0,FY0,FZ0
      REAL(8),DIMENSION(4,4,4,NXM,NYM,NZM):: U
      REAL(8),DIMENSION(NXM,NYM,NZM):: FXY0,FYZ0,FZX0,FXYZ0
      REAL(4),DIMENSION(NGM):: GX
      REAL(4),DIMENSION(NGM,5):: GF
      CHARACTER STR*80
      REAL(8):: PI,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      REAL(8):: DX0,DY0,DZ0,X,Y,Z,F,FX,FY,FZ,DX,DY,DZ
      INTEGER:: NXMIN,NXMAX,NYMIN,NYMAX,NZMIN,NZMAX
      INTEGER:: NGMAX,MODE,MODEG
      INTEGER:: NX,NY,NZ,NG,IERR
      REAL(4):: GUCLIP
C
      CALL GSOPEN
C
      PI=2.D0*ASIN(1.D0)
      NXMAX=5
      NYMAX=5
      NZMAX=5
      XMIN=0.D0
      XMAX=2.D0*PI
      YMIN=0.D0
      YMAX=2.D0*PI
      ZMIN=0.D0
      ZMAX=2.D0*PI
      NGMAX=201
      MODE=0
C
    1 WRITE(6,*) '## MODE:0..3,NXMAX,NYMAX,NZMAX='
      READ(5,*,ERR=1,END=9000) MODE,NXMAX,NYMAX,NZMAX
C
      IF(MODE.LT.4) THEN
         DX0=(XMAX-XMIN)/(NXMAX-1)
         DY0=(YMAX-YMIN)/(NYMAX-1)
         DZ0=(ZMAX-ZMIN)/(NZMAX-1)
         DO NX=1,NXMAX
            X0(NX)=XMIN+DX0*(NX-1)
         ENDDO
         DO NY=1,NYMAX
            Y0(NY)=YMIN+DY0*(NY-1)
         ENDDO
         DO NZ=1,NZMAX
            Z0(NZ)=ZMIN+DZ0*(NZ-1)
         ENDDO
         DO NZ=1,NZMAX
            Z=Z0(NZ)
         DO NY=1,NYMAX
            Y=Y0(NY)
         DO NX=1,NXMAX
            X=X0(NX)
            IF(MODE.LE.1) THEN
               F=SIN(X)*SIN(Y)*SIN(Z)
            ELSE
               IF(SQRT((X-PI)**2+(Y-PI)**2+(Z-PI)**2).LE.0.5D0*PI) THEN
                  F=1.D0
               ELSE
                  F=0.D0
               ENDIF
            ENDIF
            F0(NX,NY,NZ)=F
         ENDDO
         ENDDO
         ENDDO

C         DO NX=1,NXMAX
C         DO NY=1,NYMAX
C         DO NZ=1,NZMAX
C            WRITE(6,'(3I3,1P4E12.4)') NX,NY,NZ,
C     &           X0(NX),Y0(NY),Z0(NZ),F0(NX,NY,NZ)
C         ENDDO
C         ENDDO
C         ENDDO

C
         IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
            CALL SPL3D(X0,Y0,Z0,F0,FX0,FY0,FZ0,FXY0,FYZ0,FZX0,FXYZ0,
     &                 U,NXM,NYM,NXMAX,NYMAX,NZMAX,0,0,0,IERR)
         ELSE
            CALL SPL3D(X0,Y0,Z0,F0,FX0,FY0,FZ0,FXY0,FYZ0,FZX0,FXYZ0,
     &                 U,NXM,NYM,NXMAX,NYMAX,NZMAX,4,4,4,IERR)
         ENDIF
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX SPL3D ERROR: IERR=',IERR
            GOTO 1
         ENDIF
C
         DX=(XMAX-XMIN)/(NGMAX-1)
         DY=(YMAX-YMIN)/(NGMAX-1)
         DZ=(ZMAX-ZMIN)/(NGMAX-1)
         DO NG=1,NGMAX
            GX(NG)=GUCLIP(DBLE(NG-1)/DBLE(NGMAX-1))
            X=XMIN+DX*(NG-1)
            Y=YMIN+DY*(NG-1)
            Z=ZMIN+DZ*(NG-1)
            IF(MODE.LE.1) THEN
               F=SIN(X)*SIN(Y)*SIN(Z)
            ELSE
               IF(SQRT((X-PI)**2+(Y-PI)**2+(Z-PI)**2).LE.0.5D0*PI) THEN
                  F=1.D0
               ELSE
                  F=0.D0
               ENDIF
            ENDIF
            GF(NG,1)=GUCLIP(F)
            CALL SPL3DD(X,Y,Z,F,FX,FY,FZ,X0,Y0,Z0,
     &                  U,NXM,NYM,NXMAX,NYMAX,NZMAX,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX SPL3DD ERROR: IERR=',IERR
               GOTO 1
            ENDIF
            GF(NG,2)=GUCLIP(F)
            GF(NG,3)=GUCLIP(FX)
            GF(NG,4)=GUCLIP(FY)
            GF(NG,5)=GUCLIP(FZ)
            WRITE(6,'(I3,1P6E12.4)') NG,GX(NG),GF(NG,1),GF(NG,2),
     &                                  GF(NG,3),GF(NG,4),GF(NG,5)
         ENDDO
C
         IF(MODE.LE.1) THEN
            STR='/SIN(X)/'
         ELSE
            STR='/STEP(X)/'
         ENDIF
         MODEG=0
         CALL PAGES
         CALL GRF1D(0,GX,GF,NGM,NGMAX,5,STR,MODEG)
         CALL PAGEE
      ENDIF
      GOTO 1
C
 9000 CALL GSCLOS
      STOP
      END
