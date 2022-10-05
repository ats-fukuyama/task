C     $Id$
C
      USE libspl1d
      USE libspl2d
      USE libgrf
      IMPLICIT REAL*8 (A-F,H,O-Z)
      PARAMETER (NXM=101)
      PARAMETER (NYM=101)
      DIMENSION X0(NXM),F0(NXM),FX0(NXM),UF(4,NXM)
      DIMENSION Y0(NYM),E0(NXM,NYM)
      DIMENSION EX0(NXM,NYM),EY0(NXM,NYM),EXY0(NXM,NYM)
      DIMENSION UE(4,4,NXM,NYM)
      DIMENSION U0(NXM)
      DIMENSION GX(NXM),GF(NXM,3)
      DIMENSION GY(NYM),GE(NXM,NYM,3)
      CHARACTER STR*80
C
      CALL GSOPEN
C
      PI=2.D0*ASIN(1.D0)
      XMIN=0.D0
      XMAX=2.D0*PI
      NXMAX0=11
      NXMAX=101
      YMIN=0.D0
      YMAX=2.D0*PI
      NYMAX0=11
      NYMAX=101
      MODE=1
C
    1 WRITE(6,*) '## MODE,NXMAX0,NYMAX0='
      READ(5,*,ERR=1,END=9000) MODE,NXMAX0,NYMAX0
C
      IF(MODE.LE.4) THEN
         DX0=(XMAX-XMIN)/(NXMAX0-1)
         DO NX=1,NXMAX0
            X=XMIN+DX0*(NX-1)
            IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
               F=SIN(X)
            ELSE
               IF(ABS(X-PI).LE.0.5D0*PI) THEN
                  F=1.D0
               ELSE
                  F=0.D0
               ENDIF
            ENDIF
            X0(NX)=X
            F0(NX)=F
         ENDDO
C
         IF(MODE.EQ.1.OR.MODE.EQ.3) THEN
            CALL SPL1D(X0,F0,FX0,UF,NXMAX0,0,IERR)
         ELSE
            CALL SPL1D(X0,F0,FX0,UF,NXMAX0,4,IERR)
         ENDIF
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D ERROR: IERR=',IERR
C
         CALL SPL1DI0(X0,UF,U0,NXMAX0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1DI0 ERROR: IERR=',IERR
C
         DX=(XMAX-XMIN)/(NXMAX-1)
         DO NX=1,NXMAX
            X=XMIN+DX*(NX-1)
            CALL SPL1DD(X,F,DF,X0,UF,NXMAX0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX SPL1DD ERROR: IERR=',IERR
            CALL SPL1DI(X,FI,X0,UF,U0,NXMAX0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX SPL1DI ERROR: IERR=',IERR
            GX(NX)=GUCLIP(X)
            GF(NX,1)=GUCLIP(F)
            GF(NX,2)=GUCLIP(DF)
            GF(NX,3)=GUCLIP(FI)
         ENDDO
C
         NGMAX=3
         IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
            STR='/SIN(X)/'
         ELSE
            STR='/STEP(X)/'
         ENDIF
         MODEG=0
         CALL PAGES
         CALL GRF1D(0,GX,GF,NXM,NXMAX,NGMAX,STR,MODEG)
         CALL PAGEE
      ELSE
         DX0=(XMAX-XMIN)/(NXMAX0-1)
         DY0=(YMAX-YMIN)/(NYMAX0-1)
         DO NY=1,NYMAX0
            Y=YMIN+DY0*(NY-1)
            Y0(NY)=Y
            DO NX=1,NXMAX0
               X=XMIN+DX0*(NX-1)
               IF(MODE.EQ.5.OR.MODE.EQ.6) THEN
                  F=SIN(X)*SIN(Y)
               ELSE
                  IF((X-PI)**2+(Y-PI)**2.LE.1.D0*PI) THEN
                     F=1.D0
                  ELSE
                     F=0.D0
                  ENDIF
               ENDIF
               X0(NX)=X
               E0(NX,NY)=F
            ENDDO
         ENDDO
C
         IF(MODE.EQ.5.OR.MODE.EQ.7) THEN
            CALL SPL2D(X0,Y0,E0,EX0,EY0,EXY0,UE,
     &                 NXM,NXMAX0,NYMAX0,0,0,IERR)
         ELSE
            CALL SPL2D(X0,Y0,E0,EX0,EY0,EXY0,UE,
     &                 NXM,NXMAX0,NYMAX0,4,4,IERR)
         ENDIF
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D ERROR: IERR=',IERR
C
         DX=(XMAX-XMIN)/(NXMAX-1)
         DY=(YMAX-YMIN)/(NYMAX-1)
         DO NY=1,NYMAX
            Y=YMIN+DY*(NY-1)
            GY(NY)=GUCLIP(Y)
            DO NX=1,NXMAX
               X=XMIN+DX*(NX-1)
               CALL SPL2DD(X,Y,E,EX,EY,X0,Y0,UE,
     &                    NXM,NXMAX0,NYMAX0,IERR)
               IF(IERR.NE.0) WRITE(6,*) 'XX SPL2DD ERROR: IERR=',IERR
               GX(NX)=GUCLIP(X)
               GE(NX,NY,1)=GUCLIP(E)
               GE(NX,NY,2)=GUCLIP(EX)
               GE(NX,NY,3)=GUCLIP(EY)
            ENDDO
         ENDDO
C
         IF(MODE.EQ.5.OR.MODE.EQ.6) THEN
            STR='/SIN(X)SIN(Y)/'
         ELSE
            STR='/STEP(X)STEP(Y)/'
         ENDIF
C
         MODEG=0
         CALL PAGES
         CALL GRF2DC(0,GX,GY,GE(1,1,1),NXM,NXMAX,NYMAX,STR,MODEG,0)
         CALL PAGEE
         CALL PAGES
         CALL GRF2DB(0,GX,GY,GE(1,1,1),NXM,NXMAX,NYMAX,STR)
         CALL PAGEE
         CALL PAGES
         CALL GRF2DC(0,GX,GY,GE(1,1,2),NXM,NXMAX,NYMAX,STR,MODEG,0)
         CALL PAGEE
         CALL PAGES
         CALL GRF2DC(0,GX,GY,GE(1,1,3),NXM,NXMAX,NYMAX,STR,MODEG,0)
         CALL PAGEE
      ENDIF
      GOTO 1
C
 9000 CALL GSCLOS
      STOP
      END
