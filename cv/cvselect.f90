! cvselect.f90

MODULE cvselect
  PRIVATE
  PUBLIC cv_select
CONTAINS

  SUBROUTINE cv_select
    USE cvcomm
    IMPLICIT NONE
    CHARACTER(LEN=1):: KID
    INTEGER:: ierr

    CALL cv_select_read
1   CONTINUE
    WRITE(6,'(A)') '## Save cv-select? (Y or N)'
    READ(5,*,ERR=1,END=9) KID
    CALL ToUpper(KID)
    IF(KID.EQ.'N') GO TO 9
    IF(KID.NE.'Y') Go TO 1
    CALL cv_select_save(ierr)
9   CONTINUE
    RETURN
  END SUBROUTINE cv_select
    
! --- read selected counry-id from cv-select ---
  
  SUBROUTINE cv_select_read
    USE cvcomm
    USE libfio
    IMPLICIT NONE
    CHARACTER(LEN=2):: country_id
    INTEGER:: ncountry,nselect,ierr,nfl,nstat

    ierr=0

    nfl=22
    CALL FROPEN(nfl,knam_cv_select,1,0,'cv',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cvselect: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF
    nselect=0
    DO
       READ(NFL,'(A2)',END=3,ERR=2) country_id
       nselect=nselect+1
    END DO
2   CONTINUE
    WRITE(6,'(A,A,A,I6,I6)') &
         'XX cvselect: FILE READ ERROR with ',knam_cv_select, &
         ': nselect,iostat=',nselect,nstat
    ierr=8001
    RETURN

3   CONTINUE
    nselect_max=nselect
    WRITE(6,'(A,I8)') '## nselect_max=',nselect_max
    
    IF(ALLOCATED(country_id_nselect)) DEALLOCATE(country_id_nselect)
    ALLOCATE(country_id_nselect(nselect_max))  
    IF(ALLOCATED(country_name_nselect)) DEALLOCATE(country_name_nselect)
    ALLOCATE(country_name_nselect(nselect_max))  
    IF(ALLOCATED(ncountry_nselect)) DEALLOCATE(ncountry_nselect)
    ALLOCATE(ncountry_nselect(nselect_max))  

    REWIND(nfl)
    nselect=0
    DO
       READ(NFL,'(A2)',END=5,ERR=4) country_id
       nselect=nselect+1
       country_id_nselect(nselect)=country_id
    END DO
4   CONTINUE
    WRITE(6,'(A,A,A,I6)') &
         'XX cvselect: FILE READ ERROR with ',knam_cv_select,': iostat=',nstat
    ierr=8001
    RETURN

5   CONTINUE    

    DO nselect=1,nselect_max
       country_id=country_id_nselect(nselect)
       ncountry_nselect(nselect)=0
       DO ncountry=1,ncountry_max
          IF(country_id_ncountry(ncountry).EQ.country_id) THEN
             ncountry_nselect(nselect)=ncountry
             EXIT
          END IF
       END DO
       IF(ncountry_nselect(nselect).EQ.0) THEN
          WRITE(6,'(A,A2,A,I4)') &
               'XX cvselect error: unknown country_id: ',country_id, &
               '  nselect=',nselect
       ELSE
          ncountry=ncountry_nselect(nselect)
          WRITE(6,'(I4,2X,A2,2X,A)') &
               nselect,country_id,TRIM(country_name_ncountry(ncountry))
          country_name_nselect(nselect)=country_name_ncountry(ncountry)
       END IF
    END DO
       
    RETURN
  END SUBROUTINE cv_select_read

! --- save selec data to cv-select-cases and cv-select_deaths ---
  
  SUBROUTINE cv_select_save(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,nselect,ncountry,ndate,nstat
    CHARACTER(LEN=2):: country_id
    CHARACTER(LEN=256):: country_name
    CHARACTER(LEN=80):: format1,format2
    
    ierr=0

    IF(ALLOCATED(ncases_new_ndate_nselect)) THEN
       DEALLOCATE(ncases_new_ndate_nselect,ncases_total_ndate_nselect)
       DEALLOCATE(ndeaths_new_ndate_nselect,ndeaths_total_ndate_nselect)
    END IF
    ALLOCATE(ncases_new_ndate_nselect(ndate_max,nselect_max))
    ALLOCATE(ncases_total_ndate_nselect(ndate_max,nselect_max))
    ALLOCATE(ndeaths_new_ndate_nselect(ndate_max,nselect_max))
    ALLOCATE(ndeaths_total_ndate_nselect(ndate_max,nselect_max))

    DO nselect=1,nselect_max
       ncountry=ncountry_nselect(nselect)
       DO ndate=1,ndate_max
          ncases_new_ndate_nselect(ndate,nselect) &
               =ncases_new_ndate_ncountry(ndate,ncountry)
          ncases_total_ndate_nselect(ndate,nselect) &
               =ncases_total_ndate_ncountry(ndate,ncountry)
          ndeaths_new_ndate_nselect(ndate,nselect) &
               =ndeaths_new_ndate_ncountry(ndate,ncountry)
          ndeaths_total_ndate_nselect(ndate,nselect) &
               =ndeaths_total_ndate_ncountry(ndate,ncountry)
       END DO
    END DO

    ! --- write cases data ---

    nfl=21
    CALL FWOPEN(nfl,knam_csv_select_cases,1,0,'CV',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cv_select_save 1: FWOPEN ERROR: IERR=',ierr
       WRITE(6,*) '   file name=',knam_csv_select_cases
       RETURN
    END IF

    WRITE(format1,'(A,i8,A)') '(A,',ndate_max,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A,",",A2,',ndate_max,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'select cases,ID', &
         (date_id_ndate(ndate),ndate=1,ndate_max)
    DO nselect=1,nselect_max
       country_name=country_name_nselect(nselect)
       country_id=country_id_nselect(nselect)
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            TRIM(country_name),country_id, &
            (ncases_total_ndate_nselect(ndate,nselect),ndate=1,ndate_max)
    END DO
    CLOSE(nfl)

    WRITE(6,*) &
         '# CASES DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_select_cases)

    ! --- write deaths data ---

    nfl=21
    CALL FWOPEN(nfl,knam_csv_select_deaths,1,0,'CV',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cv_select_save 1: FWOPEN ERROR: IERR=',ierr
       RETURN
    END IF

    WRITE(format1,'(A,i8,A)') '(A,',ndate_max,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A,",",A2,',ndate_max,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'select deaths,ID', &
         (date_id_ndate(ndate),ndate=1,ndate_max)
    DO nselect=1,nselect_max
       country_name=country_name_nselect(nselect)
       country_id=country_id_nselect(nselect)
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            TRIM(country_name),country_id, &
            (ndeaths_total_ndate_nselect(ndate,nselect),ndate=1,ndate_max)
    END DO
    CLOSE(nfl)

    WRITE(6,*) &
         '# DEATHS DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_select_deaths)

    RETURN

9001   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_select_cases)
    ierr=9001
    RETURN
9002   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_select_cases)
    ierr=9002
    RETURN
9003   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_select_deaths)
    ierr=9001
    RETURN
9004   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_select_deaths)
    ierr=9002
    RETURN

  END SUBROUTINE cv_select_save
    
END MODULE cvselect
