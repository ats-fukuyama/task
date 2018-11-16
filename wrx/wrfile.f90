! wrfile.f90

MODULE wrfile

CONTAINS

!***********************************************************************
!     save ray data
!***********************************************************************

  SUBROUTINE WRSAVE

    USE wrcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: NFL,IERR,NRAY,I,NSTP

      NFL=21
      CALL FWOPEN(NFL,KNAMWR,0,MODEFW,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRSAVE: FWOPEN ERROR: IERR=',IERR
         RETURN
      END IF

      WRITE(NFL,ERR=9) NRAYMAX
      DO NRAY=1,NRAYMAX
         WRITE(NFL,ERR=9) NSTPMAX_NRAY(NRAY)
      END DO
      WRITE(6,*,ERR=9) 'NRAYMAX=',NRAYMAX
      DO NRAY=1,NRAYMAX
         WRITE(6,*,ERR=9) 'NSTPMAX=',NSTPMAX_NRAY(NRAY)
      END DO
      DO NRAY=1,NRAYMAX
         WRITE(NFL,ERR=9) (RAYIN(I,NRAY),I=1,8)
         WRITE(NFL,ERR=9) (CEXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (CEYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (CEZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RAYRB1(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RAYRB2(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         DO I=0,8
            WRITE(NFL,ERR=9) (RAYS(I,NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         ENDDO
         WRITE(NFL,ERR=9) (BNXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BNYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BNZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BABSS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
!         DO NSTP=NSTPMAX(NRAY)-10,NSTPMAX_NRAY(NRAY)
!            WRITE(6,'(A,I5,1PE12.4)') 'NSTP,BABSS=',NSTP,BABSS(NSTP,NRAY)
!         END DO
      ENDDO
      CLOSE(NFL)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',TRIM(KNAMWR)
      RETURN

    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE WRSAVE

!***********************************************************************
!     load ray data
!***********************************************************************

  SUBROUTINE WRLOAD(NSTAT)

    USE wrcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: NSTAT
    INTEGER:: NFL,IERR,I,NRAY,NSTP
    REAL(rkind):: SINP2,SINT2,ANGZ,ANGPH

      NSTAT=0

      NFL=21
      CALL FROPEN(NFL,KNAMWR,0,MODEFR,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRLOAX: FROPEN ERROR: IERR=',IERR
         RETURN
      END IF

      READ(NFL,END=8,ERR=9) NRAYMAX
      DO NRAY=1,NRAYMAX
         READ(NFL,END=8,ERR=9) NSTPMAX_NRAY(NRAY)
      END DO
      DO NRAY=1,NRAYMAX
         READ(NFL,END=8,ERR=9) (RAYIN(I,NRAY),I=1,8)
         READ(NFL,END=8,ERR=9) (CEXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (CEYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (CEZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB1(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB2(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         DO I=0,8
            READ(NFL,END=8,ERR=9) (RAYS(I,NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         ENDDO
         READ(NFL,END=8,ERR=9) (BNXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BNYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BNZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BABSS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
      ENDDO
      CLOSE(NFL)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',TRIM(KNAMWR)
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
      RNPHII=RAYIN(7,NRAY)
      UUI=RAYIN(8,NRAY)

      SINP2=RNZI**4 /(RNZI**2+RNPHII**2-RNZI**2*RNPHII**2)
      SINT2=RNPHII**4/(RNZI**2+RNPHII**2-RNZI**2*RNPHII**2)
      ANGZ= 180.D0/PI*ASIN(SQRT(SINP2))
      ANGPH=180.D0/PI*ASIN(SQRT(SINT2))
      IF(RNZI  .LT.0.D0) ANGZ= -ANGZ
      IF(RNPHII.LT.0.D0) ANGPH=-ANGPH

      RETURN

    8 WRITE(6,*) 'XX WRLOAD: End of file detected: KNAMFR= ',TRIM(KNAMWR)
      RETURN
    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE WRLOAD
END MODULE wrfile
