C 
C     *** CAUTION! ***
C
C    -------------------------------------------------------------
C   | If you want to compile this, you must remove "C" at line 63 |
C   | and 495 on eqcalc.f.                                        |
C    -------------------------------------------------------------
C
C   ************************************************
C   **                  TRTEST                    **
C   ************************************************
C
      SUBROUTINE TREQIN(IERR)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      CHARACTER KNAMEQ1*32,KNAM*32,KPNAM*32
C
      WRITE(6,*) '## TRTEST 2002/01/31'
      CALL EQINIT
C      CALL EQPARF
      MODIFY = 2
C
      CALL EQMESH
      CALL TRPSIN(IERR)
      IF(IEER.NE.0) GOTO 1000
      CALL EQLOOP(IERR)
      IF(IERR.NE.0) GOTO 1000
      CALL EQTORZ
      CALL EQSETP
C
 1000 WRITE(6,*) 'XX TRTEST: ERROR : IERR=',IERR
C      STOP
      RETURN
      END
