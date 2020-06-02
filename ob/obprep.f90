! obprep.f90

MODULE obprep

  PRIVATE
  PUBLIC ob_prep

CONTAINS
  
!     ***** Load equilibrium data *****

  SUBROUTINE ob_prep(ierr)

    USE obcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: initeq=0

    ierr=0

    IF(initeq.EQ.0) THEN
       SELECT CASE(mdlobp)
       CASE(0) ! boozer coordinates
          SELECT CASE(modelg)
          CASE(3,5)
             WRITE(line,'(A)') 'mdleqc=1'   ! set boozer poloidal angle
             CALL eqparm(2,line,ierr)
             IF(ierr.NE.0) RETURN
             CALL eqload(modelg,knameq,ierr)
             IF(ierr.NE.0) RETURN
             WRITE(line,'(A,I5)') 'nrmax =',51
             CALL eqparm(2,line,ierr)
             IF(ierr.NE.0) RETURN
             WRITE(line,'(A,I5)') 'nthmax=',64
             CALL eqparm(2,line,ierr)
             IF(ierr.NE.0) RETURN
             WRITE(line,'(A,I5)') 'nsumax=',64
             CALL eqparm(2,line,ierr)
             IF(ierr.NE.0) RETURN
             CALL eqcalq(ierr)
             IF(ierr.NE.0) RETURN
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             INITEQ=1
          END SELECT ! modelg
       END SELECT ! mdlobp
    ENDIF
    RETURN
  END SUBROUTINE ob_prep
END MODULE obprep

