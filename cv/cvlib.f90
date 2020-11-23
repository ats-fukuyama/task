! cvlib.f90

MODULE cvlib

  PRIVATE
  PUBLIC convert_date_id_to_ndate,convert_ndate_to_date_id,convert_country_id

CONTAINS

  ! --- convert string date_id to integer ndate

  SUBROUTINE convert_date_id_to_ndate(date_id,ndate)

    IMPLICIT NONE
    CHARACTER(LEN=10),INTENT(IN):: date_id
    INTEGER,INTENT(OUT):: ndate
    INTEGER:: nyear,nmonth,nday,nshift
    
    READ(date_id,'(I4,1X,I2,1X,I2)') nyear,nmonth,nday
    SELECT CASE(nyear)
    CASE(2020)
!       ndate=-3    ! shift for 2020-01-04
       ndate=-2    ! shift for 2020-01-03   
       nshift=1
    CASE(2021:2023)
       ndate=365*(nyear-2020)+1-3
       nshift=0
    END SELECT
    SELECT CASE(nmonth)
    CASE(2)
       ndate=ndate+31
    CASE(3)
       ndate=ndate+31+28+nshift
    CASE(4)
       ndate=ndate+31+28+31+nshift
    CASE(5)
       ndate=ndate+31+28+31+30+nshift
    CASE(6)
       ndate=ndate+31+28+31+30+31+nshift
    CASE(7)
       ndate=ndate+31+28+31+30+31+30+nshift
    CASE(8)
       ndate=ndate+31+28+31+30+31+30+31+nshift
    CASE(9)
       ndate=ndate+31+28+31+30+31+30+31+31+nshift
    CASE(10)
       ndate=ndate+31+28+31+30+31+30+31+31+30+nshift
    CASE(11)
       ndate=ndate+31+28+31+30+31+30+31+31+30+31+nshift
    CASE(12)
       ndate=ndate+31+28+31+30+31+30+31+31+30+31+30+nshift
    END SELECT
    ndate=ndate+nday
    RETURN
  END SUBROUTINE convert_date_id_to_ndate

  ! --- convert integer ndate to string date_id

  SUBROUTINE convert_ndate_to_date_id(ndate,date_id)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: ndate
    CHARACTER(LEN=10),INTENT(OUT):: date_id
    INTEGER:: nyear,nmonth,nday,nshift
    CHARACTER(LEN=2),DIMENSION(31):: &
         kn=(/ '01','02','03','04','05','06','07','08','09','10', &
               '11','12','13','14','15','16','17','18','19','20', &
               '21','22','23','24','25','26','27','28','29','30','31' /)
    
    SELECT CASE(ndate+20)
    CASE(1:366)
       nyear=2020
       nshift=1
!       nday=ndate+3 ! start at 2020-01-04
       nday=ndate+2 ! start at 2020-01-03
    CASE(367:366+365)
       nyear=2021
       nshift=0
       nday=ndate+3-366
    CASE(367+365:366+2*365)
       nyear=2022
       nshift=0
       nday=ndate+3-366-365
    CASE(367+2*365:366+3*365)
       nyear=2023
       nshift=0
       nday=ndate+3-366-2*365
    CASE DEFAULT
       date_id='????-??-??'
       RETURN
    END SELECT

    IF(nday.LE.31+28+nshift) THEN
       SELECT CASE(nday)
       CASE(1:31)
          nmonth=1
       CASE(32:31+28+1)
          nmonth=2
          nday=nday-31
       END SELECT
    ELSE
       SELECT CASE(nday-nshift)
       CASE(32+28:31+28+31)
          nmonth=3
          nday=nday-31-28-nshift
       CASE(32+28+31:31+28+31+30)
          nmonth=4
          nday=nday-31-28-nshift-31
       CASE(32+28+31+30:31+28+31+30+31)
          nmonth=5
          nday=nday-31-28-nshift-31-30
       CASE(32+28+31+30+31:31+28+31+30+31+30)
          nmonth=6
          nday=nday-31-28-nshift-31-30-31
       CASE(32+28+31+30+31+30:31+28+31+30+31+30+31)
          nmonth=7
          nday=nday-31-28-nshift-31-30-31-30
       CASE(32+28+31+30+31+30+31:31+28+31+30+31+30+31+31)
          nmonth=8
          nday=nday-31-28-nshift-31-30-31-30-31
       CASE(32+28+31+30+31+30+31+31:31+28+31+30+31+30+31+31+30)
          nmonth=9
          nday=nday-31-28-nshift-31-30-31-30-31-31
       CASE(32+28+31+30+31+30+31+31+30: &
            31+28+31+30+31+30+31+31+30+31)
          nmonth=10
          nday=nday-31-28-nshift-31-30-31-30-31-31-30
       CASE(32+28+31+30+31+30+31+31+30+31: &
            31+28+31+30+31+30+31+31+30+31+30)
          nmonth=11
          nday=nday-31-28-nshift-31-30-31-30-31-31-30-31
       CASE(32+28+31+30+31+30+31+31+30+31+30: &
            31+28+31+30+31+30+31+31+30+31+30+31)
          nmonth=12
          nday=nday-31-28-nshift-31-30-31-30-31-31-30-31-30
       END SELECT
    END IF
    WRITE(date_id,'(I4,A1,A2,A1,A2)') nyear,'-',kn(nmonth),'-',kn(nday)
    RETURN
  END SUBROUTINE convert_ndate_to_date_id

  ! --- convert country_id to id1 and id2 ---

  SUBROUTINE convert_country_id(country_id,id1,id2,ierr)
    IMPLICIT NONE
    CHARACTER(LEN=2),INTENT(IN):: country_id
    INTEGER,INTENT(OUT):: id1,id2,ierr

    ierr=0
    
    id1=ICHAR(country_id(1:1))
    IF(id1.GE.ICHAR('a').AND.id1.LE.ICHAR('z')) &
         id1=id1-ICHAR('a')+ICHAR('A')
    IF(id1.GE.ICHAR('A').AND.id1.LE.ICHAR('Z')) THEN
       id1=id1-ICHAR('A')+1
    ELSE
       WRITE(6,'(A,A2,4X,I4)') &
            'XX cvread: country_id is out of range: country_id,id1', &
            country_id,id1
       ierr=9101
       RETURN
    END IF
    id2=ICHAR(country_id(2:2))
    IF(id2.GE.ICHAR('a').AND.id2.LE.ICHAR('z')) &
         id2=id2-ICHAR('a')+ICHAR('A')
    IF(id2.GE.ICHAR('A').AND.id2.LE.ICHAR('Z')) THEN
       id2=id2-ICHAR('A')+1
    ELSE
       WRITE(6,'(A,A2,4X,I4)') &
            'XX cvread: country_id is out of range: country_id,id2', &
            country_id,id2
       ierr=9102
       RETURN
    END IF
    RETURN
  END SUBROUTINE convert_country_id
  
END MODULE cvlib
