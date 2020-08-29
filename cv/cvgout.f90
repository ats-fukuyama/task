! cvgout.f90

MODULE cvgout

  PRIVATE
  PUBLIC cv_gout

CONTAINS

  SUBROUTINE cv_gout
    USE cvcomm
    IMPLICIT NONE
    CHARACTER(LEN=2):: country_id
    INTEGER:: nlist_plot,nlist,id

1   CONTINUE
    WRITE(6,'(A)') '## cvgout: INPUT country_id (LEN=2):'
    READ(5,*,ERR=1,END=9000) country_id

    CALL TOUPPER(country_id)
    nlist_plot=0
    DO nlist=1,nlist_max
       IF(country_id.EQ.country_id_nlist(nlist)) THEN
          nlist_plot=nlist
          EXIT
       END IF
    END DO
    IF(nlist_plot.EQ.0) THEN
       WRITE(6,'(A,A2)') 'XX cvgout: unknown country_id:',country_id
       GO TO 1
    END IF

    WRITE(6,'(A,A)') '   Now plotting: ',country_name_nlist(nlist_plot)

2   CONTINUE
    WRITE(6,'(A)') '## INPUT graph type id: 1:6'
    READ(5,*,ERR=2,END=1) id
    IF(id.EQ.0) GO TO 2
    IF(id.LT.1.OR.id.GT.6) THEN
       WRITE(6,'(A,I)') 'XX cvgout: id out of range [1:6] : id=',id
       GO TO 2
    END IF

    CALL cv_gsub1(nlist_plot,id)

    GO TO 2

9000 CONTINUE
    RETURN
  END SUBROUTINE cv_gout

  ! plot graph default
    
  SUBROUTINE cv_gsub1(nlist,id)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nlist,id
    REAL(dp):: xg(ndate_max),fg(ndate_max,4)
    INTEGER:: ndate,ndate_ave,ndata
    CHARACTER(LEN=80):: title

    DO ndate=1,ndate_max
       xg(ndate)=DBLE(ndate)
    END DO

    SELECT CASE(id)
    CASE(1)
       DO ndate=1,ndate_max
          fg(ndate,1)=ncases_total_ndate_nlist(ndate,nlist)
       END DO
       title='@'//country_id_nlist(nlist)//': '// &
                  TRIM(country_name_nlist(nlist))//': '// &
                  region_id_nlist(nlist)//': '// &
                  'total cases@'
       ndata=1
    CASE(2)
       DO ndate=1,ndate_max
          fg(ndate,1)=ndeaths_total_ndate_nlist(ndate,nlist)
       END DO
       title='@'//country_id_nlist(nlist)//': '// &
                  TRIM(country_name_nlist(nlist))//': '// &
                  region_id_nlist(nlist)//': '// &
                  'total deaths@'
       ndata=1
    CASE(3)
       DO ndate=1,ndate_max
          fg(ndate,1)=ncases_new_ndate_nlist(ndate,nlist)
       END DO
       title='@'//country_id_nlist(nlist)//': '// &
                  TRIM(country_name_nlist(nlist))//': '// &
                  region_id_nlist(nlist)//': '// &
                  'new cases@'
       ndata=1
    CASE(4)
       DO ndate=1,ndate_max
          fg(ndate,1)=ndeaths_new_ndate_nlist(ndate,nlist)
       END DO
       title='@'//country_id_nlist(nlist)//': '// &
                  TRIM(country_name_nlist(nlist))//': '// &
                  region_id_nlist(nlist)//': '// &
                  'new deaths@'
       ndata=1
    CASE(5)
       DO ndate=1,6
          fg(ndate,1)=ncases_new_ndate_nlist(ndate,nlist)
          fg(ndate,2)=0.D0
       END DO
       DO ndate=7,ndate_max
          fg(ndate,1)=ncases_new_ndate_nlist(ndate,nlist)
          fg(ndate,2)=0.D0
          DO ndate_ave=ndate-6,ndate
             fg(ndate,2)=fg(ndate,2)+ncases_new_ndate_nlist(ndate_ave,nlist)
          END DO
          fg(ndate,2)=fg(ndate,2)/7.D0
       END DO
       title='@'//country_id_nlist(nlist)//': '// &
                  TRIM(country_name_nlist(nlist))//': '// &
                  region_id_nlist(nlist)//': '// &
                  'new cases (7days ave)@'
       ndata=2
    CASE(6)
       DO ndate=1,6
          fg(ndate,1)=ndeaths_new_ndate_nlist(ndate,nlist)
          fg(ndate,2)=0.D0
       END DO
       DO ndate=7,ndate_max
          fg(ndate,1)=ndeaths_new_ndate_nlist(ndate,nlist)
          fg(ndate,2)=0.D0
          DO ndate_ave=ndate-6,ndate
             fg(ndate,2)=fg(ndate,2)+ndeaths_new_ndate_nlist(ndate_ave,nlist)
          END DO
          fg(ndate,2)=fg(ndate,2)/7.D0
       END DO
       title='@'//country_id_nlist(nlist)//': '// &
                  TRIM(country_name_nlist(nlist))//': '// &
                  region_id_nlist(nlist)//': '// &
                  'new deaths (7 dasy ave)@'
       ndata=2
    END SELECT

    CALL PAGES
    CALL grd1d(0,xg,fg,ndate_max,ndate_max,ndata,title)
    CALL PAGEE

    RETURN
  END SUBROUTINE cv_gsub1
END MODULE cvgout
