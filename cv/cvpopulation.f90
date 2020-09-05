! cvpopulation.f90

MODULE cvpopulation

  PRIVATE
  PUBLIC cv_population_load

CONTAINS

!***********************************************************************
!     read population and are csv data
!***********************************************************************

  SUBROUTINE cv_population_load(ierr)

    USE cvcomm
    USE cvlib
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    CHARACTER(LEN=2),ALLOCATABLE:: country_id_ndata(:)
    CHARACTER(LEN=80),ALLOCATABLE:: country_name_ndata(:)
    REAL(dp),ALLOCATABLE:: population_ndata(:),area_ndata(:)
    CHARACTER(LEN=256):: line
    INTEGER:: ipos_comma(3)
    CHARACTER(LEN=2):: country_id
    CHARACTER(LEN=80):: country_name
    INTEGER:: nfl,nstat,ndata,ipos,iactive,i,ndata_max,ncountry
    REAL(dp):: population,area

    ierr=0

    nfl=21
    CALL FROPEN(nfl,knam_cv_population,1,0,'cv',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX cv_population_read: FROPEN ERROR: ierr=',ierr
       RETURN
    END IF

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
       country_id=line(1:ipos_comma(1)-1)
       IF(line(ipos_comma(1)+1:ipos_comma(1)+1).EQ.'"'.AND. &
          line(ipos_comma(2)-1:ipos_comma(2)-1).EQ.'"') THEN
          country_name=line(ipos_comma(1)+2:ipos_comma(2)-2)
       ELSE
          country_name=line(ipos_comma(1)+1:ipos_comma(2)-1)
       END IF
       READ(line(ipos_comma(2)+1:ipos_comma(3)-1),*) population
       READ(line(ipos_comma(3)+1:LEN_TRIM(line)),*) area
       ndata=ndata+1
    END DO

10  CONTINUE
    REWIND(nfl)

    ndata_max=ndata
    WRITE(6,'(A,I8)') '## ndata_max=',ndata_max

    ALLOCATE(country_id_ndata(ndata_max))
    ALLOCATE(country_name_ndata(ndata_max))
    ALLOCATE(population_ndata(ndata_max))
    ALLOCATE(area_ndata(ndata_max))

    DO ndata=1,ndata_max
       READ(NFL,'(A)',IOSTAT=nstat,END=8004,ERR=9004) line
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
       country_id_ndata(ndata)=line(1:ipos_comma(1)-1)
       IF(line(ipos_comma(1)+1:ipos_comma(1)+1).EQ.'"'.AND. &
          line(ipos_comma(2)-1:ipos_comma(2)-1).EQ.'"') THEN
          country_name_ndata(ndata)=line(ipos_comma(1)+2:ipos_comma(2)-2)
       ELSE
          country_name_ndata(ndata)=line(ipos_comma(1)+1:ipos_comma(2)-1)
       END IF
       READ(line(ipos_comma(2)+1:ipos_comma(3)-1),*) population_ndata(ndata)
       READ(line(ipos_comma(3)+1:LEN_TRIM(line)),*) area_ndata(ndata)
    END DO
    CLOSE(nfl)

    ! --- allocate and initialize population and area array

    IF(ALLOCATED(population_ncountry)) THEN
       DEALLOCATE(population_ncountry)
       DEALLOCATE(area_ncountry)
    END IF
    ALLOCATE(population_ncountry(ncountry_max))
    ALLOCATE(area_ncountry(ncountry_max))
    population_ncountry(1:ncountry_max)=0.D0
    area_ncountry(1:ncountry_max)=0.D0

    ! --- convert country_id to ncountry

    DO ncountry=1,ncountry_max
       country_id=country_id_ncountry(ncountry)
       DO ndata=1,ndata_max
          IF(country_id_ndata(ndata).EQ.country_id) THEN
             population_ncountry(ncountry)=population_ndata(ndata)
             area_ncountry(ncountry)=area_ndata(ndata)
             GO TO 20
          END IF
       END DO
       WRITE(6,'(A,I4,2X,A2)') &
            'XX cv_population: no population found for', &
            ncountry,country_id_ncountry(ncountry)
20     CONTINUE
    END DO

!    DO ndata=1,ndata_max
!       WRITE(22,'(I6,2X,A2,2ES12.4)') &
!            ndata,country_id_ndata(ndata), &
!            population_ndata(ndata),area_ndata(ndata)
!    END DO
           
!    DO ncountry=1,ncountry_max
!       WRITE(21,'(I6,2X,A2,2ES12.4,2X,A)') &
!            ncountry,country_id_ncountry(ncountry), &
!            population_ncountry(ncountry),area_ncountry(ncountry), &
!            TRIM(country_name_ncountry(ncountry))
!    END DO
           
    ! --- population read completed ---
    WRITE(6,'(A)') '## cv_population_read: completed'
    RETURN

8004 CONTINUE
    IERR=8004
8100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cvpopulation: ABNORMAL END OF FILE: IOSTAT=',nstat
    RETURN

9002 CONTINUE
    IERR=9002
    GO TO 9100
9004 CONTINUE
    IERR=9004
9100 CONTINUE    
    WRITE(6,'(A,I8)') 'XX cvpopulation: FILE READ ERROR: IOSTAT=',nstat
    RETURN

  END SUBROUTINE cv_population_load
END MODULE cvpopulation
