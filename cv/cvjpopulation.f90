! cvjpopulation.f90

MODULE cvjpopulation

  PRIVATE
  PUBLIC cvj_population_load

CONTAINS

!***********************************************************************
!     read cvj population data
!***********************************************************************

  SUBROUTINE cvj_population_load(ierr)

    USE cvcomm
    USE cvlib
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    CHARACTER(LEN=256):: line
    INTEGER:: ipos_comma(5)
    INTEGER:: nfl,iostat,ncountry,ipos,i,iactive

    ierr=0

    nfl=21
    CALL FROPEN(nfl,knam_cvj_population_in,1,0,'cv',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cvj_population_load: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF

    READ(nfl,'(A)') line

    DO ncountry=1,ncountry_max
       READ(NFL,'(A)',IOSTAT=iostat,END=8001,ERR=9001) line
!       WRITE(6,'(I4,A)') ncountry,TRIM(line)
       ipos=0
       iactive=1
       DO i=1,LEN_TRIM(line)
          IF(line(i:i).EQ.'"') THEN
             IF(iactive.EQ.0) THEN
                iactive=1
             ELSE
                iactive=0
             END IF
          END IF
          IF(iactive.EQ.1.AND.line(i:i).EQ.',') THEN
             ipos=ipos+1
             ipos_comma(ipos)=i
          END IF
       END DO
       country_name_ncountry(ncountry)=line(1:ipos_comma(1)-1)
       READ(line(ipos_comma(2)+1:ipos_comma(3)-1),*) &
            population_ncountry(ncountry)
       country_id_ncountry(ncountry) &
            =TRIM(line(ipos_comma(3)+1:ipos_comma(4)-1))
       region_id_ncountry(ncountry) &
            =TRIM(line(ipos_comma(4)+1:ipos_comma(5)-1))
    END DO

    CLOSE(nfl)

    ! --- population read completed ---
    WRITE(6,'(A)') '## cvj_population_load: completed'
    RETURN

8001 CONTINUE
    IERR=8001
    WRITE(6,'(A,I8)') &
         'XX cvj_population_load: ABNORMAL END OF FILE: IOSTAT=',iostat
    RETURN

9001 CONTINUE
    IERR=9001
    WRITE(6,'(A,I8)') &
         'XX cvj_population_load: FILE READ ERROR: IOSTAT=',iostat
    RETURN

  END SUBROUTINE cvj_population_load
END MODULE cvjpopulation
