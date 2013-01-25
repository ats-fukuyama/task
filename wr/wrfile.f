C     $Id$
C
C***********************************************************************
C     save ray data
C***********************************************************************
C
      SUBROUTINE WRSAVE
C
      USE plcomm
      INCLUDE 'wrcomm.inc'
C
      NFL=21
      CALL FWOPEN(NFL,KNAMWR,0,MODEFW,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRSAVE: FWOPEN ERROR: IERR=',IERR
         RETURN
      END IF
C
      WRITE(NFL,ERR=9) NRAYMX
      DO NRAY=1,NRAYMX
         WRITE(NFL,ERR=9) NITMAX(NRAY)
      END DO
      DO NRAY=1,NRAYMX
         WRITE(NFL,ERR=9) (RAYIN(I,NRAY),I=1,8)
         WRITE(NFL,ERR=9) (CEXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (CEYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (CEZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RKXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RKYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RKZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RAYRB1(NIT,NRAY),NIT=0,NITMAX(NRAY))
         WRITE(NFL,ERR=9) (RAYRB2(NIT,NRAY),NIT=0,NITMAX(NRAY))
         DO I=0,8
            WRITE(NFL,ERR=9) (RAYS(I,NIT,NRAY),NIT=0,NITMAX(NRAY))
         ENDDO
      ENDDO
      CLOSE(NFL)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',
     &     TRIM(KNAMWR)
      RETURN
C
    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',
     &     TRIM(KNAMWR)
      RETURN
      END
C
C***********************************************************************
C     load ray data
C***********************************************************************
C
      SUBROUTINE WRLOAD(NSTAT)
C
      USE plcomm
      INCLUDE 'wrcomm.inc'
C
      NSTAT=0
C
      NFL=21
      CALL FROPEN(NFL,KNAMWR,0,MODEFR,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRLOAX: FROPEN ERROR: IERR=',IERR
         RETURN
      END IF
C
      READ(NFL,END=8,ERR=9) NRAYMX
      DO NRAY=1,NRAYMX
         READ(NFL,END=8,ERR=9) NITMAX(NRAY)
      END DO
      DO NRAY=1,NRAYMX
         READ(NFL,END=8,ERR=9) (RAYIN(I,NRAY),I=1,8)
         READ(NFL,END=8,ERR=9) (CEXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (CEYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (CEZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RKXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RKYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RKZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB1(NIT,NRAY),NIT=0,NITMAX(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB2(NIT,NRAY),NIT=0,NITMAX(NRAY))
         DO I=0,8
            READ(NFL,END=8,ERR=9) (RAYS(I,NIT,NRAY),NIT=0,NITMAX(NRAY))
         ENDDO
      ENDDO
      CLOSE(NFL)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',
     &     TRIM(KNAMWR)
      IF(RAYRB1(1,1).EQ.0.D0.AND.RAYRB2(1,1).EQ.0.D0) THEN
         NSTAT=1
      ELSE
         NSTAT=2
      ENDIF
      RF=RAYIN(1,1)
      RPI=RAYIN(2,1)
      ZPI=RAYIN(3,1)
      PHII=RAYIN(4,1)
      RKR0=RAYIN(5,1)
      RNZI=RAYIN(6,NRAY)
      RNPHI=RAYIN(7,NRAY)
      UUI=RAYIN(8,NRAY)

      SINP2=RNZI**4 /(RNZI**2+RNPHII**2-RNZI**2*RNPHI**2)
      SINT2=RNPHI**4/(RNZI**2+RNPHII**2-RNZI**2*RNPHI**2)
      ANGZ= 180.D0/PI*ASIN(SQRT(SINP2))
      ANGPH=180.D0/PI*ASIN(SQRT(SINT2))
      IF(RNZI .LT.0.D0) ANGZ= -ANGZ
      IF(RNPHI.LT.0.D0) ANGPH=-ANGPH

      RETURN
C
    8 WRITE(6,*) 'XX WRLOAD: End of file detected: KNAMFR= ',
     &     TRIM(KNAMWR)
      RETURN
    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',
     &     TRIM(KNAMWR)
      RETURN
      END
