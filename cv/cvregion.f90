! cvregion.f90

MODULE cvregion
  PRIVATE
  PUBLIC cv_region
CONTAINS

  SUBROUTINE cv_region
    USE cvcomm
    IMPLICIT NONE
    CHARACTER(LEN=1):: KID
    INTEGER:: ierr

    CALL cv_region_count
1   CONTINUE
    WRITE(6,'(A)') '## Save cv-region? (Y or N)'
    READ(5,*,ERR=1,END=9) KID
    CALL ToUpper(KID)
    IF(KID.EQ.'N') GO TO 9
    IF(KID.NE.'Y') Go TO 1
    CALL cv_region_save(ierr)
9   CONTINUE
    RETURN
  END SUBROUTINE cv_region
    
! --- count region_id and nregion ---
  
  SUBROUTINE cv_region_count
    USE cvcomm
    IMPLICIT NONE
    CHARACTER(LEN=4):: region_id
    INTEGER:: ncountry,nregion

    IF(ALLOCATED(nregion_ncountry)) DEALLOCATE(nregion_ncountry)
    ALLOCATE(nregion_ncountry(ncountry_max))
    DO ncountry=1,ncountry_max
       nregion_ncountry(ncountry)=0
    END DO
    nregion_max=7
    IF(ALLOCATED(region_id_nregion)) DEALLOCATE(region_id_nregion)
    ALLOCATE(region_id_nregion(nregion_max))  
    IF(ALLOCATED(ncount_nregion)) DEALLOCATE(ncount_nregion)
    ALLOCATE(ncount_nregion(nregion_max))
    IF(ALLOCATED(region_name_nregion)) DEALLOCATE(region_name_nregion)
    ALLOCATE(region_name_nregion(nregion_max))  

! --- if nregion_max is unknown, calculate nregion_max
    
    IF(nregion_max.EQ.0) THEN
       nregion_max=1
       region_id_nregion(1)=region_id_ncountry(1)
       ncount_nregion(1)=1
       nregion_ncountry(1)=1
       DO ncountry=2,ncountry_max
          region_id=region_id_ncountry(ncountry)
          DO nregion=1,nregion_max
             IF(region_id.EQ.region_id_nregion(nregion)) THEN ! region found
                ncount_nregion(nregion)=ncount_nregion(nregion)+1
                nregion_ncountry(ncountry)=nregion
                EXIT
             END IF
          END DO
          IF(nregion_ncountry(ncountry).EQ.0) THEN ! region not found: new one
             nregion_max=nregion_max+1
             nregion=nregion_max
             region_id_nregion(nregion)=region_id
             ncount_nregion(nregion)=1
             nregion_ncountry(ncountry)=nregion
          END IF
       END DO
       WRITE(6,'(A,I8)') '## nregion_max=',nregion_max
    ELSE
       region_id_nregion(1)='WPRO'
       region_id_nregion(2)='EURO'
       region_id_nregion(3)='SEAR'
       region_id_nregion(4)='EMRO'
       region_id_nregion(5)='AMRO'
       region_id_nregion(6)='AFRO'
       region_id_nregion(7)='Othe'
       DO nregion=1,nregion_max
          ncount_nregion(nregion)=0
       END DO
       DO ncountry=1,ncountry_max
          region_id=region_id_ncountry(ncountry)
          DO nregion=1,nregion_max
             IF(region_id.EQ.region_id_nregion(nregion)) THEN ! region found
                ncount_nregion(nregion)=ncount_nregion(nregion)+1
                nregion_ncountry(ncountry)=nregion
                EXIT
             END IF
          END DO
       END DO
    END IF
    
    DO nregion=1,nregion_max
       SELECT CASE(region_id_nregion(nregion))
       CASE('EMRO')
          region_name_nregion(nregion)='Eastern Mediterranean'
       CASE('EURO')
          region_name_nregion(nregion)='Europe'
       CASE('AFRO')
          region_name_nregion(nregion)='Africa'
       CASE('AMRO')
          region_name_nregion(nregion)='Americas'
       CASE('WPRO')
          region_name_nregion(nregion)='Western Pacific'
       CASE('SEAR')
          region_name_nregion(nregion)='South-East Asia'
       CASE('Othe')
          region_name_nregion(nregion)='Others'
       END SELECT
    END DO
          
    DO nregion=1,nregion_max
       WRITE(6,'(I8,2X,A4,2X,I6,2X,A)') &
            nregion,region_id_nregion(nregion),ncount_nregion(nregion), &
            TRIM(region_name_nregion(nregion))
    END DO
    RETURN
  END SUBROUTINE cv_region_count

! --- save region data to cv-region-cases and cv-region_deaths ---
  
  SUBROUTINE cv_region_save(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,nregion,ncountry,ndate,nstat,ncount
    CHARACTER(LEN=4):: region_id
    CHARACTER(LEN=256):: region_name
    CHARACTER(LEN=80):: format1,format2
    
    ierr=0

    IF(ALLOCATED(ncases_new_ndate_nregion)) THEN
       DEALLOCATE(ncases_new_ndate_nregion,ncases_total_ndate_nregion)
       DEALLOCATE(ndeaths_new_ndate_nregion,ndeaths_total_ndate_nregion)
    END IF
    ALLOCATE(ncases_new_ndate_nregion(ndate_max,nregion_max))
    ALLOCATE(ncases_total_ndate_nregion(ndate_max,nregion_max))
    ALLOCATE(ndeaths_new_ndate_nregion(ndate_max,nregion_max))
    ALLOCATE(ndeaths_total_ndate_nregion(ndate_max,nregion_max))

    ! --- sum up region data ---

    DO nregion=1,nregion_max
       DO ndate=1,ndate_max
          ncases_new_ndate_nregion(ndate,nregion)=0
          ncases_total_ndate_nregion(ndate,nregion)=0
          ndeaths_new_ndate_nregion(ndate,nregion)=0
          ndeaths_total_ndate_nregion(ndate,nregion)=0
       END DO
    END DO

    DO ncountry=1,ncountry_max
       nregion=nregion_ncountry(ncountry)
       DO ndate=1,ndate_max
          ncases_new_ndate_nregion(ndate,nregion) &
               =ncases_new_ndate_nregion(ndate,nregion) &
               +ncases_new_ndate_ncountry(ndate,ncountry)
          ncases_total_ndate_nregion(ndate,nregion) &
               =ncases_total_ndate_nregion(ndate,nregion) &
               +ncases_total_ndate_ncountry(ndate,ncountry)
          ndeaths_new_ndate_nregion(ndate,nregion) &
               =ndeaths_new_ndate_nregion(ndate,nregion) &
               +ndeaths_new_ndate_ncountry(ndate,ncountry)
          ndeaths_total_ndate_nregion(ndate,nregion) &
               =ndeaths_total_ndate_nregion(ndate,nregion) &
               +ndeaths_total_ndate_ncountry(ndate,ncountry)
       END DO
    END DO

    ! --- write cases data ---

    nfl=21
    CALL FWOPEN(nfl,knam_csv_out_region,1,0,'CV',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cv_region_save : FWOPEN ERROR: IERR=',ierr
       WRITE(6,*) '   file name=',knam_csv_out_region
       RETURN
    END IF

    ncount=(ndate_max-ndate_start_region)/ndate_step_region+1
    WRITE(format1,'(A,i8,A)') '(A,',ncount,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A,',ncount,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'region cases', &
         (date_id_ndate(ndate), &
          ndate=ndate_start_region,ndate_max,ndate_step_region)
    DO nregion=1,nregion_max
       region_name=region_name_nregion(nregion)
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            TRIM(region_name), &
            (ncases_total_ndate_nregion(ndate,nregion), &
             ndate=ndate_start_region,ndate_max,ndate_step_region)
    END DO

    WRITE(6,*) &
         '# CASES DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_out_region)

    ! --- write deaths data ---

    WRITE(format1,'(A,i8,A)') '(A,',ncount,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A,',ncount,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'region deaths', &
         (date_id_ndate(ndate), &
          ndate=ndate_start_region,ndate_max,ndate_step_region)
    DO nregion=1,nregion_max
       region_name=region_name_nregion(nregion)
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            TRIM(region_name), &
            (ndeaths_total_ndate_nregion(ndate,nregion), &
             ndate=ndate_start_region,ndate_max,ndate_step_region)
    END DO
    CLOSE(nfl)

    WRITE(6,*) &
         '# DEATHS DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_out_region)

    RETURN

9001 CONTINUE
    WRITE(6,'(A/A,I8,A)') &
         'XX cv_region_save 9001: File IO error detected:', &
         '   nstat,KNAMFR= ',nstat,TRIM(knam_csv_out_region)
    ierr=9001
    RETURN
9002 CONTINUE
    WRITE(6,'(A/A,I8,A)') &
         'XX cv_region_save 9002: File IO error detected:', &
         '   nstat,KNAMFR= ',nstat,TRIM(knam_csv_out_region)
    ierr=9002
    RETURN

  END SUBROUTINE cv_region_save
    
END MODULE cvregion
