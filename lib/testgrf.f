C     $Id$
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      PARAMETER (NXM=1001)
      PARAMETER (NYM=101)
      PARAMETER (NGM=10)
      DIMENSION GX(NXM),GF(NXM,NGM)
      DIMENSION GY(NYM),GZ(NXM,NYM),KA(8,NXM,NYM)
      CHARACTER STR*80
C
      CALL GSOPEN
C
      XMIN=0.D0
      XMAX=4.D0*ASIN(1.D0)
      NXMAX=101
      DX=(XMAX-XMIN)/(NXMAX-1)
      DO NX=1,NXMAX
         X=XMIN+DX*(NX-1)
         GX(NX)=GUCLIP(X)
         GF(NX,1)=GUCLIP(SIN(X))
         GF(NX,2)=GUCLIP(SIN(2*X))
         GF(NX,3)=GUCLIP(SIN(3*X))
         GF(NX,4)=GUCLIP(SIN(4*X))
         GF(NX,5)=GUCLIP(SIN(5*X))
         GF(NX,6)=GUCLIP(SIN(6*X))
      ENDDO
      NGMAX=6
      STR='/SIN(n*X)/'
      MODE=0
      CALL PAGES
      CALL GRF1D(0,GX,GF,NXM,NXMAX,NGMAX,STR,MODE)
      CALL PAGEE
C
      XMIN=0.001D0
      XMAX=1.D0
      NXMAX=101
      DX=(XMAX-XMIN)/(NXMAX-1)
      DO NX=1,NXMAX
         X=XMIN+DX*(NX-1)
         GX(NX)=GUCLIP(LOG10(X))
         GF(NX,1)=GUCLIP(LOG10(EXP(-X)))
         GF(NX,2)=GUCLIP(LOG10(EXP(-2*X)))
         GF(NX,3)=GUCLIP(LOG10(EXP(-3*X)))
         GF(NX,4)=GUCLIP(LOG10(EXP(-4*X)))
         GF(NX,5)=GUCLIP(LOG10(EXP(-5*X)))
         GF(NX,6)=GUCLIP(LOG10(EXP(-6*X)))
      ENDDO
      NGMAX=6
      STR='/EXP(n*X)/'
      MODE=3
      CALL PAGES
      CALL GRF1D(0,GX,GF,NXM,NXMAX,NGMAX,STR,MODE)
      CALL PAGEE
C
      XMIN=0.D0
      XMAX=4.D0*ASIN(1.D0)
      YMIN=0.D0
      YMAX=4.D0*ASIN(1.D0)
      NXMAX=101
      NYMAX=101
      DX=(XMAX-XMIN)/(NXMAX-1)
      DY=(YMAX-YMIN)/(NYMAX-1)
      DO NX=1,NXMAX
         X=XMIN+DX*(NX-1)
         GX(NX)=GUCLIP(X)
      ENDDO
      DO NY=1,NYMAX
         Y=YMIN+DY*(NY-1)
         GY(NY)=GUCLIP(Y)
      ENDDO
      DO NY=1,NYMAX
         Y=YMIN+DY*(NY-1)
         DO NX=1,NXMAX
            X=XMIN+DX*(NX-1)
            GZ(NX,NY)=GUCLIP(COS(X)*COS(Y))
         ENDDO
      ENDDO
      STR='/COS(X)COS(Y)/'
C
      MODE=0
      CALL PAGES
      CALL GRF2DC(0,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KA,MODE)
      CALL PAGEE
C
      CALL PAGES
      CALL GRF2DB(0,GX,GY,GZ,NXM,NXMAX,NYMAX,STR)
      CALL PAGEE
C
      CALL GSCLOS
      STOP
      END
