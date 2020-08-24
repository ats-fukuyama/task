! cvfile.f90

MODULE cvfile

  PRIVATE
  PUBLIC cv_read,cv_write

CONTAINS

!***********************************************************************
!     read WHO csv data
!***********************************************************************

  SUBROUTINE cv_read(ierr)

    USE cvcomm
    USE cvlib
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,nstat,ndata,ncountry,i,ndata_max,id1,id2,ipos,iactive
    INTEGER:: nlist,ndate,nyear,nmonth,nday,nshift
    CHARACTER(LEN=256):: line
    INTEGER:: ipos_comma(7)
    CHARACTER(LEN=256):: column_name(8)
    CHARACTER(LEN=10):: date_id
    CHARACTER(LEN=2):: country_id
    CHARACTER(LEN=80):: country_name
    CHARACTER(LEN=4):: region_name
    CHARACTER(LEN=10),ALLOCATABLE:: date_id_ndata(:)
    CHARACTER(LEN=2),ALLOCATABLE:: country_id_ndata(:)
    CHARACTER(LEN=80),ALLOCATABLE:: country_name_ndata(:)
    CHARACTER(LEN=4),ALLOCATABLE:: region_name_ndata(:)
    INTEGER:: ncases_new,ncases_total,ndeaths_new,ndeaths_total
    INTEGER,ALLOCATABLE:: &
         ncases_new_ndata(:),ncases_total_ndata(:), &
         ndeaths_new_ndata(:),ndeaths_total_ndata(:)
    INTEGER,ALLOCATABLE::&
         nlist_ndata(:),ndate_ndata(:)
    INTEGER:: ndata_id1_id2(27,26),nlist_id1_id2(27,26)

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
       region_name=line(ipos_comma(3)+1:ipos_comma(4)-1)
       READ(line(ipos_comma(4)+1:ipos_comma(5)-1),*) ncases_new
       READ(line(ipos_comma(5)+1:ipos_comma(6)-1),*) ncases_total
       READ(line(ipos_comma(6)+1:ipos_comma(7)-1),*) ndeaths_new
       READ(line(ipos_comma(7)+1:LEN_TRIM(line)),*) ndeaths_total
       ndata=ndata+1
    END DO

10  CONTINUE
    ndata_max=ndata
    WRITE(6,'(A,I8)') '## ndata_max=',ndata_max

    IF(ALLOCATED(date_id_ndata)) THEN
       DEALLOCATE(date_id_ndata)
       DEALLOCATE(ncases_new_ndata)
       DEALLOCATE(ncases_total_ndata)
       DEALLOCATE(ndeaths_new_ndata)
       DEALLOCATE(ndeaths_total_ndata)
       DEALLOCATE(nlist_ndata)
       DEALLOCATE(ndate_ndata)
    END IF
    ALLOCATE(date_id_ndata(ndata_max))
    ALLOCATE(country_id_ndata(ndata_max))
    ALLOCATE(country_name_ndata(ndata_max))
    ALLOCATE(region_name_ndata(ndata_max))
    ALLOCATE(ncases_new_ndata(ndata_max))
    ALLOCATE(ncases_total_ndata(ndata_max))
    ALLOCATE(ndeaths_new_ndata(ndata_max))
    ALLOCATE(ndeaths_total_ndata(ndata_max))
    ALLOCATE(nlist_ndata(ndata_max))
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
       IF(line(ipos_comma(2)+1:ipos_comma(2)+1).EQ.'"'.AND. &
          line(ipos_comma(3)-1:ipos_comma(3)-1).EQ.'"') THEN
          country_name_ndata(ndata)=line(ipos_comma(2)+2:ipos_comma(3)-2)
       ELSE
          country_name_ndata(ndata)=line(ipos_comma(2)+1:ipos_comma(3)-1)
       END IF
       region_name_ndata(ndata)=line(ipos_comma(3)+1:ipos_comma(4)-1)
       READ(line(ipos_comma(4)+1:ipos_comma(5)-1),*) ncases_new_ndata(ndata)
       READ(line(ipos_comma(5)+1:ipos_comma(6)-1),*) ncases_total_ndata(ndata)
       READ(line(ipos_comma(6)+1:ipos_comma(7)-1),*) ndeaths_new_ndata(ndata)
       READ(line(ipos_comma(7)+1:LEN_TRIM(line)),*) ndeaths_total_ndata(ndata)
    END DO
    CLOSE(nfl)

    ! --- initialize table to convert from country_id to ndata
    ndata_id1_id2(1:27,1:26)=0
    nlist_id1_id2(1:27,1:26)=0
    nlist_ndata(1:ndata_max)=0

    ! --- list up unique country_id
    nlist=0
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
          nlist=nlist+1
          ndata_id1_id2(id1,id2)=ndata
          nlist_id1_id2(id1,id2)=nlist
       END IF
       nlist_ndata(ndata)=nlist_id1_id2(id1,id2)
    END DO
    nlist_max=nlist
    WRITE(6,'(A,I8)') '## nlist_max=',nlist_max

    ! --- make lists of country_id, country_name and region_name
    
    IF(ALLOCATED(country_id_nlist)) THEN
       DEALLOCATE(country_id_nlist)
       DEALLOCATE(country_name_nlist)
       DEALLOCATE(region_name_nlist)
    END IF
    ALLOCATE(country_id_nlist(nlist_max))
    ALLOCATE(country_name_nlist(nlist_max))
    ALLOCATE(region_name_nlist(nlist_max))
    
    DO ndata=1,ndata_max
       IF(nlist_ndata(ndata).NE.0) THEN
          nlist=nlist_ndata(ndata)
          country_id_nlist(nlist)=country_id_ndata(ndata)
          country_name_nlist(nlist)=country_name_ndata(ndata)
          region_name_nlist(nlist)=region_name_ndata(ndata)
       END IF
    END DO

    DEALLOCATE(country_id_ndata,country_name_ndata,region_name_ndata)

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

    ! --- convert (ndata) array to (ndate,nlist) array ---

    IF(ALLOCATED(ncases_new_ndate_nlist)) THEN
       DEALLOCATE(ncases_new_ndate_nlist)
       DEALLOCATE(ncases_total_ndate_nlist)
       DEALLOCATE(ndeaths_new_ndate_nlist)
       DEALLOCATE(ndeaths_total_ndate_nlist)
    END IF
    ALLOCATE(ncases_new_ndate_nlist(ndate_max,nlist_max))
    ALLOCATE(ncases_total_ndate_nlist(ndate_max,nlist_max))
    ALLOCATE(ndeaths_new_ndate_nlist(ndate_max,nlist_max))
    ALLOCATE(ndeaths_total_ndate_nlist(ndate_max,nlist_max))

    ncases_new_ndate_nlist(1:ndate_max,1:nlist_max)=0
    ncases_total_ndate_nlist(1:ndate_max,1:nlist_max)=0
    ndeaths_new_ndate_nlist(1:ndate_max,1:nlist_max)=0
    ndeaths_total_ndate_nlist(1:ndate_max,1:nlist_max)=0
    
    DO ndata=1,ndata_max
       nlist=nlist_ndata(ndata)
       ndate=ndate_ndata(ndata)
       ncases_new_ndate_nlist(ndate,nlist)=ncases_new_ndata(ndata)
       ncases_total_ndate_nlist(ndate,nlist)=ncases_total_ndata(ndata)
       ndeaths_new_ndate_nlist(ndate,nlist)=ndeaths_new_ndata(ndata)
       ndeaths_total_ndate_nlist(ndate,nlist)=ndeaths_total_ndata(ndata)
    END DO

    ! --- correct missing data ---

    DO nlist=1,nlist_max
       ncases_total=ncases_total_ndate_nlist(1,nlist)
       ndeaths_total=ndeaths_total_ndate_nlist(1,nlist)
       DO ndate=2,ndate_max
          IF(ncases_total.NE.0) THEN
             IF(ncases_total_ndate_nlist(ndate,nlist).EQ.0) &
                ncases_total_ndate_nlist(ndate,nlist)=ncases_total
             IF(ndeaths_total_ndate_nlist(ndate,nlist).EQ.0) &
                ndeaths_total_ndate_nlist(ndate,nlist)=ndeaths_total
          END IF
          ncases_total=ncases_total_ndate_nlist(1,nlist)
          ndeaths_total=ndeaths_total_ndate_nlist(1,nlist)
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

  END SUBROUTINE cv_read

!***********************************************************************
!     write csv data
!***********************************************************************

  SUBROUTINE cv_write(ierr)

    USE cvcomm
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nfl,nlist,ndate,nstat
    CHARACTER(LEN=2):: country_id
    CHARACTER(LEN=256):: country_name
    CHARACTER(LEN=80):: format1,format2
    
    ierr=0

    ! --- write cases data ---
    
    nfl=21
    CALL FWOPEN(nfl,knam_csv_out_cases,1,0,'CV',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX CVSAVE: FWOPEN ERROR: IERR=',ierr
       RETURN
    END IF

    WRITE(format1,'(A,i8,A)') '(A,',ndate_max,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A2,",",A,",",A4,',ndate_max,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9001) &
         'country_id,country_name,region_id', &
         (date_id_ndate(ndate),ndate=1,ndate_max)
    DO nlist=1,nlist_max
       country_id=country_id_nlist(nlist)
       IF(country_id.EQ.'PS'.OR.country_id.EQ.'BQ') THEN
          country_name='"'//TRIM(country_name_nlist(nlist))//'"'
       ELSE
          country_name=country_name_nlist(nlist)
       END IF
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9002) &
            country_id_nlist(nlist),TRIM(country_name), &
            region_name_nlist(nlist), &
            (ncases_total_ndate_nlist(ndate,nlist),ndate=1,ndate_max)
    END DO
    CLOSE(nfl)

    WRITE(6,*) &
         '# CASES DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_out_cases)

    ! --- write deaths data ---
    
    nfl=21
    CALL FWOPEN(nfl,knam_csv_out_deaths,1,0,'CV',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX CVSAVE: FWOPEN ERROR: IERR=',ierr
       RETURN
    END IF

    WRITE(format1,'(A,i8,A)') '(A,',ndate_max,'(",",A10))'
    WRITE(format2,'(A,i8,A)') '(A2,",",A,",",A4,',ndate_max,'(",",I0))'
    WRITE(nfl,format1,IOSTAT=nstat,ERR=9003) &
         'country_id,country_name,region_id', &
         (date_id_ndate(ndate),ndate=1,ndate_max)
    DO nlist=1,nlist_max
       country_id=country_id_nlist(nlist)
       IF(country_id.EQ.'PS'.OR.country_id.EQ.'BQ') THEN
          country_name='"'//TRIM(country_name_nlist(nlist))//'"'
       ELSE
          country_name=country_name_nlist(nlist)
       END IF
       WRITE(nfl,format2,IOSTAT=nstat,ERR=9004) &
            country_id_nlist(nlist),TRIM(country_name), &
            region_name_nlist(nlist), &
            (ndeaths_total_ndate_nlist(ndate,nlist),ndate=1,ndate_max)
    END DO
    CLOSE(nfl)

    WRITE(6,*) &
         '# DEATHS DATA WAS SUCCESSFULLY SAVED TO THE FILE: ', &
         TRIM(knam_csv_out_deaths)

    RETURN

9001   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_out_cases)
    ierr=9001
    RETURN
9002   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_out_cases)
    ierr=9002
    RETURN
9003   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_out_deaths)
    ierr=9001
    RETURN
9004   WRITE(6,'(A,I8,A)') &
         'XX cv_write: File IO error detected: nstat,KNAMFR= ', &
         nstat,TRIM(knam_csv_out_deaths)
    ierr=9002
    RETURN
  END SUBROUTINE cv_write

END MODULE cvfile
