! cvfile.f90

MODULE cvfile

  PRIVATE
  PUBLIC cv_load,cv_global

CONTAINS

!***********************************************************************
!     load WHO csv data
!***********************************************************************

  SUBROUTINE cv_load(ierr)

    USE cvcomm
    USE cvlib
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,nstat,ndata,ncountry,i,ndata_max,id1,id2,ipos,iactive
    INTEGER:: ndate,nyear,nmonth,nday,nshift
    CHARACTER(LEN=256):: line
    INTEGER:: ipos_comma(7)
    CHARACTER(LEN=256):: column_name(8)
    CHARACTER(LEN=10):: date_id
    CHARACTER(LEN=2):: country_id
    CHARACTER(LEN=80):: country_name
    CHARACTER(LEN=4):: region_id
    CHARACTER(LEN=10),ALLOCATABLE:: date_id_ndata(:)
    CHARACTER(LEN=2),ALLOCATABLE:: country_id_ndata(:)
    CHARACTER(LEN=80),ALLOCATABLE:: country_name_ndata(:)
    CHARACTER(LEN=4),ALLOCATABLE:: region_id_ndata(:)
    INTEGER:: ncases_new,ncases_total,ndeaths_new,ndeaths_total
    INTEGER,ALLOCATABLE:: &
         ncases_new_ndata(:),ncases_total_ndata(:), &
         ndeaths_new_ndata(:),ndeaths_total_ndata(:)
    INTEGER,ALLOCATABLE::&
         ncountry_ndata(:),ndate_ndata(:)
    INTEGER:: ndata_id1_id2(27,26),ncountry_id1_id2(27,26)

    ierr=0

    nfl=21
    CALL FROPEN(nfl,knam_csv_in,1,0,'cv',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cvread: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF

    READ(NFL,'(8A)',END=8001,ERR=9001) (column_name(i),i=1,8)
    ndata=0
    DO
       READ(NFL,'(A)',IOSTAT=nstat,END=10,ERR=9002) line
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
       date_id=line(1:ipos_comma(1)-1)
       country_id=line(ipos_comma(1)+1:ipos_comma(2)-1)
       country_name=line(ipos_comma(2)+1:ipos_comma(3)-1)
       region_id=line(ipos_comma(3)+1:ipos_comma(4)-1)
       READ(line(ipos_comma(4)+1:ipos_comma(5)-1),*) ncases_new
       READ(line(ipos_comma(5)+1:ipos_comma(6)-1),*) ncases_total
       READ(line(ipos_comma(6)+1:ipos_comma(7)-1),*) ndeaths_new
       READ(line(ipos_comma(7)+1:LEN_TRIM(line)),*) ndeaths_total
       ndata=ndata+1
    END DO

10  CONTINUE
    ndata_max=ndata
    WRITE(6,'(A,I8)') '## ndata_max=',ndata_max

    ALLOCATE(date_id_ndata(ndata_max))
    ALLOCATE(country_id_ndata(ndata_max))
    ALLOCATE(country_name_ndata(ndata_max))
    ALLOCATE(region_id_ndata(ndata_max))
    ALLOCATE(ncases_new_ndata(ndata_max))
    ALLOCATE(ncases_total_ndata(ndata_max))
    ALLOCATE(ndeaths_new_ndata(ndata_max))
    ALLOCATE(ndeaths_total_ndata(ndata_max))
    ALLOCATE(ncountry_ndata(ndata_max))
    ALLOCATE(ndate_ndata(ndata_max))

    REWIND(nfl)
    READ(NFL,'(8A)',IOSTAT=nstat,END=8003,ERR=9003) (column_name(i),i=1,8)
    DO ndata=1,ndata_max
       READ(NFL,'(A)',IOSTAT=nstat,END=8004,ERR=9004) line
       ipos=0
       iactive=1
       DO i=1,LEN_TRIM(line)
!          IF(ndata.EQ.24814) WRITE(6,'(I4,A1,3I4)') &
!               i,line(i:i),iactive,ipos,ipos_comma(MAX(ipos,1))
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
       date_id_ndata(ndata)=line(1:ipos_comma(1)-1)
       country_id_ndata(ndata)=line(ipos_comma(1)+1:ipos_comma(2)-1)
       IF(country_id_ndata(ndata).EQ.' ') country_id_ndata(ndata)='OT'
       IF(line(ipos_comma(2)+1:ipos_comma(2)+1).EQ.'"'.AND. &
          line(ipos_comma(3)-1:ipos_comma(3)-1).EQ.'"') THEN
          country_name_ndata(ndata)=line(ipos_comma(2)+2:ipos_comma(3)-2)
       ELSE
          country_name_ndata(ndata)=line(ipos_comma(2)+1:ipos_comma(3)-1)
       END IF
       region_id_ndata(ndata)=line(ipos_comma(3)+1:ipos_comma(4)-1)
       READ(line(ipos_comma(4)+1:ipos_comma(5)-1),*) ncases_new_ndata(ndata)
       READ(line(ipos_comma(5)+1:ipos_comma(6)-1),*) ncases_total_ndata(ndata)
       READ(line(ipos_comma(6)+1:ipos_comma(7)-1),*) ndeaths_new_ndata(ndata)
       READ(line(ipos_comma(7)+1:LEN_TRIM(line)),*) ndeaths_total_ndata(ndata)
    END DO
    CLOSE(nfl)

    ! --- initialize table to convert from country_id to ndata
    ndata_id1_id2(1:27,1:26)=0
    ncountry_id1_id2(1:27,1:26)=0
    ncountry_ndata(1:ndata_max)=0

    ! --- list up unique country_id
    ncountry=0
    DO ndata=1,ndata_max
       IF(country_id_ndata(ndata).EQ.' ') THEN
          id1=27
          id2=1
       ELSE
          CALL convert_country_id(country_id_ndata(ndata),id1,id2,ierr)
          IF(ierr.NE.0) THEN
             WRITE(6,'(A,A,I8,2X,A2,4X,I4)') &
                  'XX cvread: country_id is out of range:', &
                  'ndata,country_id,id1,id2', &
                  ndata,country_id_ndata(ndata),id1,id2
             RETURN
          END IF
       END IF

       IF(ndata_id1_id2(id1,id2).EQ.0) THEN
          ncountry=ncountry+1
          ndata_id1_id2(id1,id2)=ndata
          ncountry_id1_id2(id1,id2)=ncountry
       END IF
       ncountry_ndata(ndata)=ncountry_id1_id2(id1,id2)
    END DO
    ncountry_max=ncountry
    WRITE(6,'(A,I8)') '## ncountry_max=',ncountry_max

    ! --- make lists of country_id, country_name and region_id
    
    IF(ALLOCATED(country_id_ncountry)) THEN
       DEALLOCATE(country_id_ncountry)
       DEALLOCATE(country_name_ncountry)
       DEALLOCATE(region_id_ncountry)
    END IF
    ALLOCATE(country_id_ncountry(ncountry_max))
    ALLOCATE(country_name_ncountry(ncountry_max))
    ALLOCATE(region_id_ncountry(ncountry_max))
    
    DO ndata=1,ndata_max
       IF(ncountry_ndata(ndata).NE.0) THEN
          ncountry=ncountry_ndata(ndata)
          country_id_ncountry(ncountry)=country_id_ndata(ndata)
          country_name_ncountry(ncountry)=country_name_ndata(ndata)
          region_id_ncountry(ncountry)=region_id_ndata(ndata)
       END IF
    END DO

    DEALLOCATE(country_id_ndata,country_name_ndata,region_id_ndata)

    ! --- convert date string to ndate --- ndata=1 for 2020-01-21

    ndate_max=0
    DO ndata=1,ndata_max
       date_id=date_id_ndata(ndata)
       CALL convert_date_id_to_ndate(date_id,ndate)
       ndate_ndata(ndata)=ndate
       IF(ndate.GT.ndate_max) ndate_max=ndate
       IF(ndate.LE.0) THEN
          WRITE(6,'(A,I10,2X,A10)') &
               'XX cvread: negative ndate for ndata:',ndata,date_id
          ierr=9100
          RETURN
       END IF
    END DO
    IF(ALLOCATED(date_id_ndate)) DEALLOCATE(date_id_ndate)
    ALLOCATE(date_id_ndate(ndate_max))
    CALL convert_ndate_to_date_id(ndate_max,date_id_ndate(ndate_max))
    WRITE(6,'(A,I8,2X,A10)') &
         '## ndate_max=',ndate_max,date_id_ndate(ndate_max)
    DO ndate=1,ndate_max-1
       CALL convert_ndate_to_date_id(ndate,date_id_ndate(ndate))
    END DO

    ! --- convert (ndata) array to (ndate,ncountry) array ---

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
    
    DO ndata=1,ndata_max
       ncountry=ncountry_ndata(ndata)
       ndate=ndate_ndata(ndata)
       ncases_new_ndate_ncountry(ndate,ncountry)=ncases_new_ndata(ndata)
       ncases_total_ndate_ncountry(ndate,ncountry)=ncases_total_ndata(ndata)
       ndeaths_new_ndate_ncountry(ndate,ncountry)=ndeaths_new_ndata(ndata)
       ndeaths_total_ndate_ncountry(ndate,ncountry)=ndeaths_total_ndata(ndata)
    END DO

    ! --- correct missing data ---

    DO ncountry=1,ncountry_max
       ncases_total=ncases_total_ndate_ncountry(1,ncountry)
       ndeaths_total=ndeaths_total_ndate_ncountry(1,ncountry)
       DO ndate=2,ndate_max
          IF(ncases_total.NE.0) THEN
             IF(ncases_total_ndate_ncountry(ndate,ncountry).EQ.0) &
                ncases_total_ndate_ncountry(ndate,ncountry)=ncases_total
             IF(ndeaths_total_ndate_ncountry(ndate,ncountry).EQ.0) &
                ndeaths_total_ndate_ncountry(ndate,ncountry)=ndeaths_total
          END IF
          ncases_total=ncases_total_ndate_ncountry(1,ncountry)
          ndeaths_total=ndeaths_total_ndate_ncountry(1,ncountry)
       END DO
    END DO
       
    ! --- read completed ---
    WRITE(6,'(A)') '## cvread: completed'
    RETURN

8001 CONTINUE
    IERR=8001
    GO TO 8100
8003 CONTINUE
    IERR=8003
    GO TO 8100
8004 CONTINUE
    IERR=8004
8100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cvread: ABNORMAL END OF FILE: IOSTAT=',nstat
    RETURN

9001 CONTINUE
    IERR=9001
    GO TO 9100
9002 CONTINUE
    IERR=9002
    GO TO 9100
9003 CONTINUE
    IERR=9003
    GO TO 9100
9004 CONTINUE
    IERR=9004
9100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cvread: FILE READ ERROR: IOSTAT=',nstat
    RETURN

  END SUBROUTINE cv_load

!***********************************************************************
!     write global csv data
!***********************************************************************

  SUBROUTINE cv_global(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,ncountry,ndate,nstat,ncount
    CHARACTER(LEN=2):: country_id
    CHARACTER(LEN=256):: country_name
    CHARACTER(LEN=80):: format1,format2
    
    ierr=0

    ! --- write cases data ---
    
    nfl=21
    CALL FWOPEN(nfl,knam_csv_out_global,1,0,'CV',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cv_global: FWOPEN ERROR: IERR=',ierr
       WRITE(6,*) '   FILE NAME=',knam_csv_out_global
       RETURN
    END IF

    ncount=(ndate_max-ndate_start_global)/ndate_step_global+1
    WRITE(format1,'(A,i8,A)') '(A,',ncount,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A2,",",A,",",A4,',ncount,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'country_id,country_name,region_id', &
         (date_id_ndate(ndate), &
         ndate=ndate_start_global,ndate_max,ndate_step_global)
    DO ncountry=1,ncountry_max
       country_id=country_id_ncountry(ncountry)
       IF(country_id.EQ.'PS'.OR.country_id.EQ.'BQ') THEN
          country_name='"'//TRIM(country_name_ncountry(ncountry))//'"'
       ELSE
          country_name=country_name_ncountry(ncountry)
       END IF
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            country_id_ncountry(ncountry),TRIM(country_name), &
            region_id_ncountry(ncountry), &
            (ncases_total_ndate_ncountry(ndate,ncountry), &
            ndate=ndate_start_global,ndate_max,ndate_step_global)
    END DO

    WRITE(6,*) &
         '# CASES DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_out_global)

    ! --- write deaths data ---
    

    WRITE(format1,'(A,i8,A)') '(A,',ncount,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A2,",",A,",",A4,',ncount,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'country_id,country_name,region_id', &
         (date_id_ndate(ndate), &
         ndate=ndate_start_global,ndate_max,ndate_step_global)
    DO ncountry=1,ncountry_max
       country_id=country_id_ncountry(ncountry)
       IF(country_id.EQ.'PS'.OR.country_id.EQ.'BQ') THEN
          country_name='"'//TRIM(country_name_ncountry(ncountry))//'"'
       ELSE
          country_name=country_name_ncountry(ncountry)
       END IF
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            country_id_ncountry(ncountry),TRIM(country_name), &
            region_id_ncountry(ncountry), &
            (ndeaths_total_ndate_ncountry(ndate,ncountry), &
            ndate=ndate_start_global,ndate_max,ndate_step_global)
    END DO
    CLOSE(nfl)

    WRITE(6,*) &
         '# DEATHS DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_out_global)

    RETURN

9001 CONTINUE
    WRITE(6,'(A/A,I8,A)') &
         'XX cv_global 9001: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_out_global)
    ierr=9001
    RETURN
9002 CONTINUE
    WRITE(6,'(A,I8,A)') &
         'XX cv_write 9002: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_out_global)
    ierr=9002
    RETURN
  END SUBROUTINE cv_global

END MODULE cvfile
