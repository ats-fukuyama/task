! cvrank.f90

MODULE cvrank

  PRIVATE
  PUBLIC cv_rank,cv_rank_exec

CONTAINS

  SUBROUTINE cv_rank
    USE cvcomm
    IMPLICIT none
    INTEGER:: id,id1,id2
    INTEGER:: nrank,ncountry,ndate
    REAL(dp):: rdata

1   CONTINUE
    WRITE(6,'(A)') '# Input ID: 0 for all, or 1-12, or -1 for end'
    READ(5,*,ERR=1,END=9000) id
    IF(id.GT.idrank_max) GO TO 1
    IF(id.LT.0) GO TO 9000

    CALL cv_rank_exec(id)
    
    IF(id.EQ.0) THEN
       id1=1
       id2=idrank_max
    ELSE
       id1=id
       id2=id
    END IF
    
    DO id=id1,id2
       SELECT CASE(id)
       CASE(1)
          WRITE(6,'(A)') 'nrank: total cases: number:'
       CASE(2)
          WRITE(6,'(A)') 'nrank: new cases: number:'
       CASE(3)
          WRITE(6,'(A)') 'nrank: average cases: number:'
       CASE(4)
          WRITE(6,'(A)') 'nrank: total deaths: number:'
       CASE(5)
          WRITE(6,'(A)') 'nrank: new deaths: number:'
       CASE(6)
          WRITE(6,'(A)') 'nrank: average deaths: number:'
       CASE(7)
          WRITE(6,'(A)') 'nrank: total cases: rate:'
       CASE(8)
          WRITE(6,'(A)') 'nrank: new cases: rate:'
       CASE(9)
          WRITE(6,'(A)') 'nrank: average cases: rate:'
       CASE(10)
          WRITE(6,'(A)') 'nrank: total deaths: rate:'
       CASE(11)
          WRITE(6,'(A)') 'nrank: new deaths: rate:'
       CASE(12)
          WRITE(6,'(A)') 'nrank: average deaths: rate:'
       END SELECT

       DO nrank=1,nrank_max
          rdata=data_nrank_idrank(nrank,id)
          ncountry=ncountry_nrank_idrank(nrank,id)
          WRITE(6,'(I3,I6,1X,A2,F12.2,I10,I8,I8,I6,2X,A20)') &
               nrank,ncountry,country_id_ncountry(ncountry),rdata, &
               ncases_total_ndate_ncountry(ndate_max,ncountry), &
               ncases_new_ndate_ncountry(ndate_max,ncountry), &
               ndeaths_total_ndate_ncountry(ndate_max,ncountry), &
               ndeaths_new_ndate_ncountry(ndate_max,ncountry), &
               TRIM(country_name_ncountry(ncountry))
       END DO

    END DO
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE cv_rank

  SUBROUTINE cv_rank_exec(id_)
    USE cvcomm
    IMPLICIT none
    INTEGER,INTENT(IN):: id_
    INTEGER:: id1,id2,id
    INTEGER:: nrank,ncountry,nrank1,nvalue1,ncountry1,ndate
    INTEGER,SAVE:: nrank_max_save=0
    REAL(dp):: rdata
    REAL(dp):: data_nrank(nrank_max)

    IF(id_.EQ.0) THEN
       id1=1
       id2=12
    ELSE
       id1=id_
       id2=id_
    END IF
    
    IF(nrank_max.NE.nrank_max_save) THEN
       IF(ALLOCATED(ncountry_nrank_idrank)) THEN
          DEALLOCATE(ncountry_nrank_idrank)
          DEALLOCATE(data_nrank_idrank)
       END IF

       ALLOCATE(ncountry_nrank_idrank(nrank_max,idrank_max))
       ALLOCATE(data_nrank_idrank(nrank_max,idrank_max))
    END IF
    
    DO id=id1,id2
       ncountry_nrank_idrank(1:nrank_max,id)=0
       data_nrank_idrank(1:nrank_max,id)=0.D0

       DO ncountry=1,ncountry_max
          IF(population_ncountry(ncountry).LE.population_min_rank) CYCLE
          SELECT CASE(id)
          CASE(1,7)
             rdata=DBLE(ncases_total_ndate_ncountry(ndate_max,ncountry))
          CASE(2,8)
             rdata=DBLE(ncases_new_ndate_ncountry(ndate_max,ncountry))
          CASE(3,9)
             rdata=0.D0
             DO ndate=ndate_max-ndays_ave+1,ndate_max
                rdata=rdata+DBLE(ncases_new_ndate_ncountry(ndate,ncountry))
             END DO
             rdata=rdata/DBLE(ndays_ave)
          CASE(4,10)
             rdata=DBLE(ndeaths_total_ndate_ncountry(ndate_max,ncountry))
          CASE(5,11)
             rdata=DBLE(ndeaths_new_ndate_ncountry(ndate_max,ncountry))
          CASE(6,12)
             rdata=0.D0
             DO ndate=ndate_max-ndays_ave+1,ndate_max
                rdata=rdata+DBLE(ndeaths_new_ndate_ncountry(ndate,ncountry))
             END DO
             rdata=rdata/DBLE(ndays_ave)
          END SELECT
          SELECT CASE(id)
          CASE(7:12)
             rdata=rdata/population_ncountry(ncountry)
          END SELECT
          IF(rdata.GE.data_nrank_idrank(nrank_max,id)) THEN
             nrank=nrank_max
10           CONTINUE
             IF(nrank.EQ.1) GO TO 20
             nrank=nrank-1
             IF(rdata.GE.data_nrank_idrank(nrank,id)) GO TO 10
             nrank=nrank+1
20           CONTINUE
             DO nrank1=nrank_max-1,nrank,-1
                data_nrank_idrank(nrank1+1,id) &
                     =data_nrank_idrank(nrank1,id)
                ncountry_nrank_idrank(nrank1+1,id) &
                     =ncountry_nrank_idrank(nrank1,id)
             END DO
             data_nrank_idrank(nrank,id)=rdata
             ncountry_nrank_idrank(nrank,id)=ncountry
          END IF
       END DO
    END DO
    
    RETURN
  END SUBROUTINE cv_rank_exec
END MODULE cvrank
