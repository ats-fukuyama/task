!     ******* DISPERSION RELATION *******
C
      IF(NDISP.EQ.1) THEN
         XMIN=-RA
         XMAX= RA
C         IF(NGRAPH.GT.0) THEN
C            IF(IGRAPH.EQ.0) CALL GSOPEN
C            IGRAPH=1
C         ENDIF
  100    WRITE(6,*) '## INPUT : XMIN,XMAX,KXMIN,KXMAX,NXP ?'
         READ(5,*,ERR=100,END=2) XMIN,XMAX,PKXMIN,PKXMAX,NXP
         IF(ABS(XMAX-XMIN).LE.1.D-32) GOTO 2
         DO 110 NX=1,NXP
            XAM(NX)=XMIN+(XMAX-XMIN)*(NX-0.5)/NXP
  110    CONTINUE
         CALL W1PROF
         NZ=1
         CALL W1DSPA(NALPHA)
         IF(NMODEL.LE.3) THEN
            CALL W1WKXB
         ELSE
            CALL W1WKXD
         ENDIF
         CALL W1GDSP(PKXMIN,PKXMAX)
         GOTO 100
      ENDIF
