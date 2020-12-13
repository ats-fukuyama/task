! cvjfile.f90

MODULE cvjfile

  PRIVATE
  PUBLIC cvj_load

CONTAINS

!***********************************************************************
!     load cvj data
!***********************************************************************

  SUBROUTINE cvj_load(ierr)

    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: ncountry,ndate

    CALL cvj_load_prep(ierr)
    IF(ierr.NE.0) RETURN

    CALL cvj_load_cases(ierr)
    IF(ierr.NE.0) RETURN

    CALL cvj_load_deaths(ierr)
    IF(ierr.NE.0) RETURN

    DO ncountry=1,ncountry_max
       ncases_new_ndate_ncountry(1,ncountry)=0
       ndeaths_new_ndate_ncountry(1,ncountry)=0
       DO ndate=2,ndate_max
          ncases_new_ndate_ncountry(ndate,ncountry) &
               =ncases_total_ndate_ncountry(ndate,ncountry) &
               -ncases_total_ndate_ncountry(ndate-1,ncountry)
          ndeaths_new_ndate_ncountry(ndate,ncountry) &
               =ndeaths_total_ndate_ncountry(ndate,ncountry) &
               -ndeaths_total_ndate_ncountry(ndate-1,ncountry)
       END DO
    END DO

    RETURN
  END SUBROUTINE cvj_load

  ! *** count cvj date and allocate data array ***
  
  SUBROUTINE cvj_load_prep(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,iostat,ncountry,ndate,i,ipos
    CHARACTER(LEN=3650):: long_line

    ierr=0

    nfl=21
    CALL FROPEN(nfl,knam_cvj_cases_in,1,0,'cvj',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cvj_load: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF

    READ(NFL,'(A)',IOSTAT=iostat,END=8001,ERR=9001) long_line
    ndate=0
    ipos=0
    DO i=1,LEN_TRIM(long_line)
       IF(long_line(i:i).EQ.',') THEN
          ndate=ndate+1
          ipos=i
       END IF
    END DO
    ndate_max=ndate  ! first column is prefecture name
    ndate_max_g=ndate_max
    WRITE(6,'(A,I8,2X,A10)') &
         '## ndate_max=',ndate_max,TRIM(long_line(ipos+1:LEN_TRIM(long_line)))

    ncountry=0
    DO
       READ(NFL,'(A)',IOSTAT=iostat,END=8002,ERR=9002) long_line
       ncountry=ncountry+1
    END DO

8002 CONTINUE
    ncountry_max=ncountry-1
    WRITE(6,'(A,I8)') &
         '## ncountry_max=',ncountry_max

    REWIND(NFL)

    IF(ALLOCATED(date_id_ndate)) DEALLOCATE(date_id_ndate)
    ALLOCATE(date_id_ndate(ndate_max))

    IF(ALLOCATED(country_id_ncountry)) THEN
       DEALLOCATE(country_id_ncountry)
       DEALLOCATE(country_name_ncountry)
       DEALLOCATE(region_id_ncountry)
       DEALLOCATE(population_ncountry)
    END IF
    ALLOCATE(country_id_ncountry(ncountry_max))
    ALLOCATE(country_name_ncountry(ncountry_max))
    ALLOCATE(region_id_ncountry(ncountry_max))
    ALLOCATE(population_ncountry(ncountry_max))

    IF(ALLOCATED(ncases_new_ndate_ncountry)) THEN
       DEALLOCATE(ncases_new_ndate_ncountry)
       DEALLOCATE(ncases_total_ndate_ncountry)
       DEALLOCATE(ndeaths_new_ndate_ncountry)
       DEALLOCATE(ndeaths_total_ndate_ncountry)
    END IF
    ALLOCATE(ncases_new_ndate_ncountry(ndate_max,ncountry_max))
    ALLOCATE(ncases_total_ndate_ncountry(ndate_max,ncountry_max))
    ALLOCATE(ndeaths_new_ndate_ncountry(ndate_max,ncountry_max))
    ALLOCATE(ndeaths_total_ndate_ncountry(ndate_max,ncountry_max))
    ncases_new_ndate_ncountry(1:ndate_max,1:ncountry_max)=0
    ncases_total_ndate_ncountry(1:ndate_max,1:ncountry_max)=0
    ndeaths_new_ndate_ncountry(1:ndate_max,1:ncountry_max)=0
    ndeaths_total_ndate_ncountry(1:ndate_max,1:ncountry_max)=0

    WRITE(6,'(A)') '## cvj_load_prep: completed'
    RETURN

8001 CONTINUE
    IERR=8001
    GO TO 8100
8100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cv_load_prep: ABNORMAL END OF FILE: IOSTAT=',iostat
    RETURN

9001 CONTINUE
    IERR=9001
    GO TO 9100
9002 CONTINUE
    IERR=9002
    GO TO 9100
9100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cv_load_prep: FILE READ ERROR: IOSTAT=',iostat
    RETURN

  END SUBROUTINE cvj_load_prep

  ! *** load cvj cases data ***

  SUBROUTINE cvj_load_cases(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,iostat,ncountry,ndate,i,ipos
    CHARACTER(LEN=3650):: long_line

    ierr=0

    nfl=21
    CALL FROPEN(nfl,knam_cvj_cases_in,1,0,'cvj',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cvj_load: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF

    READ(NFL,'(A)',IOSTAT=iostat,END=8003,ERR=9003) long_line
    ndate=0
    ipos=0
    DO i=1,LEN_TRIM(long_line)
       IF(long_line(i:i).EQ.',') THEN
          IF(ndate.GE.1) date_id_ndate(ndate)=long_line(ipos+1:i-1)
          ndate=ndate+1
          ipos=i
       END IF
    END DO
    date_id_ndate(ndate)=long_line(ipos+1:i-1)
    IF(ndate.NE.ndate_max) THEN
       WRITE(6,'(A,2I6)') &
            'XX cvj_load: ndate.NE.ndate_max:',ndate,ndate_max
       ierr=8001
       RETURN
    END IF

    DO ncountry=1,ncountry_max
       READ(NFL,'(A)',IOSTAT=iostat,END=8004,ERR=9004) long_line
       ndate=0
       ipos=0
       DO i=1,LEN_TRIM(long_line)
          IF(long_line(i:i).EQ.',') THEN
             IF(ndate.EQ.0) THEN
                country_name_ncountry(ncountry)=TRIM(long_line(ipos+1:i-1))
             ELSE
                READ(long_line(ipos+1:i-1),'(I8)') &
                     ncases_total_ndate_ncountry(ndate,ncountry)
             END IF
             ndate=ndate+1
             ipos=i
          END IF
       END DO
       READ(long_line(ipos+1:LEN_TRIM(long_line)),'(I8)') &
            ncases_total_ndate_ncountry(ndate,ncountry)
    END DO

    REWIND(nfl)

    ! --- read completed ---
    WRITE(6,'(A)') '## cvj_load_cases: completed'
    RETURN

8003 CONTINUE
    IERR=8003
    GO TO 8100
8004 CONTINUE
    IERR=8004
8100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cv_load_cases: ABNORMAL END OF FILE: IOSTAT=',iostat
    RETURN

9003 CONTINUE
    IERR=9003
    GO TO 9100
9004 CONTINUE
    IERR=9004
9100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cv_load_cases: FILE READ ERROR: IOSTAT=',iostat
    RETURN

  END SUBROUTINE cvj_load_cases

  ! *** load cvj deaths data ***

  SUBROUTINE cvj_load_deaths(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,iostat,ncountry,ndate,i,ipos
    CHARACTER(LEN=3650):: long_line

    ierr=0

    nfl=21
    CALL FROPEN(nfl,knam_cvj_deaths_in,1,0,'cvj',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cvj_load: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF

    READ(NFL,'(A)',IOSTAT=iostat,END=8003,ERR=9003) long_line
    ndate=0
    ipos=0
    DO i=1,LEN_TRIM(long_line)
       IF(long_line(i:i).EQ.',') THEN
          IF(ndate.GE.1) date_id_ndate(ndate)=long_line(ipos+1:i-1)
          ndate=ndate+1
          ipos=i
       END IF
    END DO
    date_id_ndate(ndate)=long_line(ipos+1:i-1)
    IF(ndate.NE.ndate_max) THEN
       WRITE(6,'(A,2I6)') &
            'XX cvj_load: ndate.NE.ndate_max:',ndate,ndate_max
       ierr=8001
       RETURN
    END IF

    DO ncountry=1,ncountry_max
       READ(NFL,'(A)',IOSTAT=iostat,END=8004,ERR=9004) long_line
       ndate=0
       ipos=0
       DO i=1,LEN_TRIM(long_line)
          IF(long_line(i:i).EQ.',') THEN
             IF(ndate.EQ.0) THEN
                country_name_ncountry(ncountry)=TRIM(long_line(ipos+1:i-1))
             ELSE
                READ(long_line(ipos+1:i-1),'(I8)') &
                     ndeaths_total_ndate_ncountry(ndate,ncountry)
             END IF
             ndate=ndate+1
             ipos=i
          END IF
       END DO
       READ(long_line(ipos+1:LEN_TRIM(long_line)),'(I8)') &
            ndeaths_total_ndate_ncountry(ndate,ncountry)
    END DO

    REWIND(nfl)

    ! --- read completed ---
    WRITE(6,'(A)') '## cvj_load_deaths: completed'
    RETURN

8003 CONTINUE
    IERR=8003
    GO TO 8100
8004 CONTINUE
    IERR=8004
8100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cv_load_deaths: ABNORMAL END OF FILE: IOSTAT=',iostat
    RETURN

9003 CONTINUE
    IERR=9003
    GO TO 9100
9004 CONTINUE
    IERR=9004
9100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cv_load_deaths: FILE READ ERROR: IOSTAT=',iostat
    RETURN

  END SUBROUTINE cvj_load_deaths

END MODULE cvjfile
