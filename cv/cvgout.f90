! cvgout.f90

MODULE cvgout

  PRIVATE
  PUBLIC cv_gout

CONTAINS

  SUBROUTINE cv_gout
    USE cvcomm
    IMPLICIT NONE
    CHARACTER(LEN=2):: country_id,kid
    INTEGER:: ncountry_plot,ncountry,id,idx

    id=0

1   CONTINUE
    
    WRITE(6,'(A,A)') &
         '## cvgout: INPUT country_id (LEN=2): e.g. JP for Japan,', &
         'XL for list, XX for end'
    READ(5,*,ERR=1,END=9000) country_id
    GO TO 5
    
2   CONTINUE
    WRITE(6,'(A)') &
         '## INPUT country id or graph id: XL:list XG:help XX or 9 for end'
    READ(5,*,ERR=2,END=9000) kid
    READ(kid,'(I2)',ERR=4) idx
    id=idx
    IF(id.EQ.9) GO TO 9000

3   CONTINUE
    SELECT CASE(id)
    CASE(1,10:14)
       CALL cv_gsub1(ncountry_plot,id)
    CASE(2,20:24)
       CALL cv_gsub2(ncountry_plot,id)
    CASE(3,30:34)
       CALL cv_gsub3(ncountry_plot,id)
    CASE DEFAULT
       WRITE(6,'(A,I)') 'XX cvgout: graph id out of range : id=',id
    END SELECT
    GO TO 2

4   CONTINUE
    country_id=kid
5   CONTINUE
    CALL TOUPPER(country_id)
    IF(country_id.EQ.'XX') GO TO 9000
    IF(country_id.EQ.'XL') THEN
       WRITE(6,'(25(1X,A2))') &
            (country_id_ncountry(ncountry),ncountry=1,ncountry_max)
       GO TO 2
    END IF
    IF(country_id.EQ.'XG') THEN
       WRITE(6,'((A))') &
            ' 1: 11+12+13+14', &
            '    11: total cases vs time', &
            '    12: total deaths vs time', &
            '    13: new cases vs time', &
            '    14: new deaths vs time', &
            ' 2: 21+22+23+24 (per million population) ', &
            '    21: total cases per population vs time', &
            '    22: total deaths per population vs time', &
            '    23: new cases per population vs time', &
            '    24: new deaths per population vs time', &
            ' 3: 31+32+33+34 ', &
            '    31: deaths number vs cases number (fixed range)', &
            '    32: deaths rate vs cases rate (fixed range)', &
            '    33: deaths number vs cases number (adjusted range)', &
            '    34: deaths rate vs cases rate (adjusted range)'
       GO TO 2
    END IF
    ncountry_plot=0
    DO ncountry=1,ncountry_max
       IF(country_id.EQ.country_id_ncountry(ncountry)) THEN
          ncountry_plot=ncountry
          EXIT
       END IF
    END DO
    IF(ncountry_plot.EQ.0) THEN
       WRITE(6,'(A,A2)') 'XX cvgout: unknown country_id:',country_id
       GO TO 1
    END IF

    WRITE(6,'(A,A)') &
         '   Now plotting: ',TRIM(country_name_ncountry(ncountry_plot))
    
    IF(id.EQ.0) GO TO 2
    GO TO 3
    

9000 CONTINUE
    RETURN
  END SUBROUTINE cv_gout

  ! plot graph by number
    
  SUBROUTINE cv_gsub1(ncountry,id)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: ncountry,id
    REAL(dp):: xg(ndate_max),fg(ndate_max,4)
    INTEGER:: ndata,ndate
    CHARACTER(LEN=80):: title
    INTEGER,PARAMETER:: nlmax=1
    REAL(dp):: line_rgb(3,nlmax)=(/ 1.D0,0.D0,0.D0 /)
    REAL(dp):: line_mark_size(nlmax)=(/ 0.3D0 /)
    INTEGER:: line_mark(nlmax)=(/ -1 /)
    INTEGER:: line_mark_step(nlmax)=(/ 14 /)
    INTEGER,PARAMETER:: nlmax2=2
    REAL(dp):: line_rgb2(3,nlmax2)=(/ 0.D0,0.D0,1.D0, 1.D0,0.D0,0.D0 /)
    REAL(dp):: line_mark_size2(nlmax2)=(/ 0.3D0, 0.3D0 /)
    INTEGER:: line_mark2(nlmax2)=(/ 0, -1 /)
    INTEGER:: line_mark_step2(nlmax2)=(/ 1, 14 /)

    DO ndate=1,ndate_max
       xg(ndate)=DBLE(ndate)
    END DO

    SELECT CASE(id)
    CASE(1,10)
       CALL PAGES
       CALL cv_gsub11(1,ncountry,fg,title,ndata)
       CALL grd1d(1,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub11(2,ncountry,fg,title,ndata)
       CALL grd1d(2,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub11(3,ncountry,fg,title,ndata)
       CALL grd1d(3,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb2,LINE_MARK=line_mark2, &
                  LINE_MARK_STEP=line_mark_step2, &
                  LINE_MARK_SIZE=line_mark_size2,NLMAX=nlmax2)
       CALL cv_gsub11(4,ncountry,fg,title,ndata)
       CALL grd1d(4,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb2,LINE_MARK=line_mark2, &
                  LINE_MARK_STEP=line_mark_step2, &
                  LINE_MARK_SIZE=line_mark_size2,NLMAX=nlmax2)
       CALL PAGEE
    CASE(11:12)
       CALL cv_gsub11(id-10,ncountry,fg,title,ndata)
       CALL PAGES
       CALL grd1d(0,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    CASE(13:14)
       CALL cv_gsub11(id-10,ncountry,fg,title,ndata)
       CALL PAGES
       CALL grd1d(0,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb2,LINE_MARK=line_mark2, &
                  LINE_MARK_STEP=line_mark_step2, &
                  LINE_MARK_SIZE=line_mark_size2,NLMAX=nlmax2)
       CALL PAGEE
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub1

  ! --- set up data '''

  SUBROUTINE cv_gsub11(id,ncountry,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,ncountry
    REAL(dp),INTENT(OUT):: fg(ndate_max,4)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate,ndate_ave

    SELECT CASE(id)
    CASE(1)
       DO ndate=1,ndate_max
          fg(ndate,1)=ncases_total_ndate_ncountry(ndate,ncountry)
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'total cases@'
       ndata=1
    CASE(2)
       DO ndate=1,ndate_max
          fg(ndate,1)=ndeaths_total_ndate_ncountry(ndate,ncountry)
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'total deaths@'
       ndata=1
    CASE(3)
       DO ndate=1,6
          fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry)
          fg(ndate,2)=0.D0
       END DO
       DO ndate=7,ndate_max
          fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry)
          fg(ndate,2)=0.D0
          DO ndate_ave=ndate-6,ndate
             fg(ndate,2)=fg(ndate,2) &
                  +ncases_new_ndate_ncountry(ndate_ave,ncountry)
          END DO
          fg(ndate,2)=fg(ndate,2)/7.D0
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'new cases (7days ave)@'
       ndata=2
    CASE(4)
       DO ndate=1,6
          fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry)
          fg(ndate,2)=0.D0
       END DO
       DO ndate=7,ndate_max
          fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry)
          fg(ndate,2)=0.D0
          DO ndate_ave=ndate-6,ndate
             fg(ndate,2)=fg(ndate,2) &
                  +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
          END DO
          fg(ndate,2)=fg(ndate,2)/7.D0
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'new deaths (7 days ave)@'
       ndata=2
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub11
    
  ! plot graph by rate
    
  SUBROUTINE cv_gsub2(ncountry,id)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: ncountry,id
    REAL(dp):: xg(ndate_max),fg(ndate_max,4)
    INTEGER:: ndata,ndate
    CHARACTER(LEN=80):: title
    INTEGER,PARAMETER:: nlmax=1
    REAL(dp):: line_rgb(3,nlmax)=(/ 1.D0,0.D0,0.D0 /)
    REAL(dp):: line_mark_size(nlmax)=(/ 0.3D0 /)
    INTEGER:: line_mark(nlmax)=(/ -1 /)
    INTEGER:: line_mark_step(nlmax)=(/ 14 /)
    INTEGER,PARAMETER:: nlmax2=2
    REAL(dp):: line_rgb2(3,nlmax2)=(/ 0.D0,0.D0,1.D0, 1.D0,0.D0,0.D0 /)
    REAL(dp):: line_mark_size2(nlmax2)=(/ 0.3D0, 0.3D0 /)
    INTEGER:: line_mark2(nlmax2)=(/ 0, -1 /)
    INTEGER:: line_mark_step2(nlmax2)=(/ 1, 14 /)

    DO ndate=1,ndate_max
       xg(ndate)=DBLE(ndate)
    END DO

    SELECT CASE(id)
    CASE(2,20)
       CALL PAGES
       CALL cv_gsub21(1,ncountry,fg,title,ndata)
       CALL grd1d(1,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub21(2,ncountry,fg,title,ndata)
       CALL grd1d(2,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub21(3,ncountry,fg,title,ndata)
       CALL grd1d(3,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb2,LINE_MARK=line_mark2, &
                  LINE_MARK_STEP=line_mark_step2, &
                  LINE_MARK_SIZE=line_mark_size2,NLMAX=nlmax2)
       CALL cv_gsub21(4,ncountry,fg,title,ndata)
       CALL grd1d(4,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb2,LINE_MARK=line_mark2, &
                  LINE_MARK_STEP=line_mark_step2, &
                  LINE_MARK_SIZE=line_mark_size2,NLMAX=nlmax2)
       CALL PAGEE
    CASE(21:22)
       CALL cv_gsub21(id-20,ncountry,fg,title,ndata)
       CALL PAGES
       CALL grd1d(0,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    CASE(23:24)
       CALL cv_gsub21(id-20,ncountry,fg,title,ndata)
       CALL PAGES
       CALL grd1d(0,xg,fg,ndate_max,ndate_max,ndata,title, &
                  LINE_RGB=line_rgb2,LINE_MARK=line_mark2, &
                  LINE_MARK_STEP=line_mark_step2, &
                  LINE_MARK_SIZE=line_mark_size2,NLMAX=nlmax2)
       CALL PAGEE
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub2

  ! --- set up data '''

  SUBROUTINE cv_gsub21(id,ncountry,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,ncountry
    REAL(dp),INTENT(OUT):: fg(ndate_max,4)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate,ndate_ave

    SELECT CASE(id)
    CASE(1)
       DO ndate=1,ndate_max
          fg(ndate,1)=ncases_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'total cases/pop(M)@'
       ndata=1
    CASE(2)
       DO ndate=1,ndate_max
          fg(ndate,1)=ndeaths_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'total deaths/pop(M)@'
       ndata=1
    CASE(3)
       DO ndate=1,6
          fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
          fg(ndate,2)=0.D0
       END DO
       DO ndate=7,ndate_max
          fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
          fg(ndate,2)=0.D0
          DO ndate_ave=ndate-6,ndate
             fg(ndate,2)=fg(ndate,2) &
                  +ncases_new_ndate_ncountry(ndate_ave,ncountry)
          END DO
          fg(ndate,2)=fg(ndate,2)/(7.D0*population_ncountry(ncountry))
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'new cases (7days ave)/pop(M)@'
       ndata=2
    CASE(4)
       DO ndate=1,6
          fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
          fg(ndate,2)=0.D0
       END DO
       DO ndate=7,ndate_max
          fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
          fg(ndate,2)=0.D0
          DO ndate_ave=ndate-6,ndate
             fg(ndate,2)=fg(ndate,2) &
                  +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
          END DO
          fg(ndate,2)=fg(ndate,2)/(7.D0*population_ncountry(ncountry))
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'new deaths (7 days ave)@'
       ndata=2
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub21
    
  ! plot graph by rate
    
  SUBROUTINE cv_gsub3(ncountry,id)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: ncountry,id
    REAL(dp):: xg(ndate_max),fg(ndate_max,4)
    INTEGER:: ndata,ndate
    CHARACTER(LEN=80):: title
    INTEGER,PARAMETER:: nlmax=1
    REAL(dp):: line_rgb(3,nlmax)=(/ 1.D0,0.D0,0.D0 /)
    REAL(dp):: line_mark_size(nlmax)=(/ 0.3D0 /)
    INTEGER:: line_mark(nlmax)=(/ -1 /)
    INTEGER:: line_mark_step(nlmax)=(/ 14 /)
    REAL(dp):: scale_thickness=0.035D0

    SELECT CASE(id)
    CASE(3,30)
       CALL PAGES
       CALL cv_gsub31(1,ncountry,xg,fg,title,ndata)
       CALL grd1d(1,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=2.477D0,XMAX=7.D0,FMIN=1.D0,FMAX=5.477D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub31(2,ncountry,xg,fg,title,ndata)
       CALL grd1d(2,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=0.477D0,XMAX=4.477D0,FMIN=-0.523D0,FMAX=3.477D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub31(3,ncountry,xg,fg,title,ndata)
       CALL grd1d(3,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=3.D0,FMIN=1.D0, &
                  ASPECT=0.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL cv_gsub31(4,ncountry,xg,fg,title,ndata)
       CALL grd1d(4,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=1.D0,FMIN=-1.D0, &
                  ASPECT=0.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    CASE(31)
       CALL cv_gsub31(1,ncountry,xg,fg,title,ndata)
       CALL PAGES
       CALL grd1d(1,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=2.477D0,XMAX=7.D0,FMIN=1.D0,FMAX=5.477D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    CASE(32)
       CALL cv_gsub31(2,ncountry,xg,fg,title,ndata)
       CALL PAGES
       CALL grd1d(2,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=0.477D0,XMAX=4.477D0,FMIN=-0.523D0,FMAX=3.477D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    CASE(33)
       CALL cv_gsub31(3,ncountry,xg,fg,title,ndata)
       CALL PAGES
       CALL grd1d(1,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=3.D0,FMIN=1.D0, &
                  ASPECT=0.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    CASE(34)
       CALL cv_gsub31(4,ncountry,xg,fg,title,ndata)
       CALL PAGES
       CALL grd1d(2,xg,fg,ndate_max,ndate_max,ndata,title, &
                  XMIN=1.D0,FMIN=-1.D0, &
                  ASPECT=0.D0,XSCALE_ZERO=0,FSCALE_ZERO=0,MODE_LS=3, &
                  SCALE_THICKNESS=scale_thickness, &
                  XSCALE_TYPE=0,FSCALE_TYPE=0, &
                  LINE_RGB=line_rgb,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
       CALL PAGEE
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub3

  ! --- set up data '''

  SUBROUTINE cv_gsub31(id,ncountry,xg,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,ncountry
    REAL(dp),INTENT(OUT):: xg(ndate_max),fg(ndate_max,4)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate
    REAL(dp):: rcases,rdeaths

    SELECT CASE(id)
    CASE(1,3)
       DO ndate=1,ndate_max
          rcases=DBLE(ncases_total_ndate_ncountry(ndate,ncountry))
          rdeaths=DBLE(ndeaths_total_ndate_ncountry(ndate,ncountry))
          xg(ndate)=  LOG10(MAX(0.01D0,rcases))
          fg(ndate,1)=LOG10(MAX(0.01D0,rdeaths))
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'deaths vs cases in number@'
       ndata=1
    CASE(2,4)
       DO ndate=1,ndate_max
          rcases=DBLE(ncases_total_ndate_ncountry(ndate,ncountry)) &
                     /population_ncountry(ncountry)
          rdeaths=DBLE(ndeaths_total_ndate_ncountry(ndate,ncountry)) &
                     /population_ncountry(ncountry)
          xg(ndate)=  LOG10(MAX(0.01D0,rcases))
          fg(ndate,1)=LOG10(MAX(0.01D0,rdeaths))
       END DO
       title='@'//country_id_ncountry(ncountry)//': '// &
                  TRIM(country_name_ncountry(ncountry))//': '// &
                  region_id_ncountry(ncountry)//': '// &
                  'deaths vs cases in rate@'
       ndata=1
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub31
    
END MODULE cvgout
