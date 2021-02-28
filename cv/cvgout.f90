! cvgout.f90

MODULE cvgout

  USE cvcomm,ONLY: dp
  
  PRIVATE
  PUBLIC cv_gout

  INTEGER,PARAMETER:: nlmax=24
  REAL(dp):: line_rgb(3,nlmax),line_mark_size(nlmax)
  INTEGER:: line_pat(nlmax),line_mark(nlmax),line_mark_step(nlmax)

  INTEGER,PARAMETER:: nlmax1=6
  REAL(dp),PARAMETER:: line_rgb1(3,nlmax1) &
       =RESHAPE([1.0D0,0.0D0,0.0D0, &
                 1.0D0,0.8D0,0.0D0, &
                 0.0D0,0.8D0,0.0D0, &
                 0.0D0,0.8D0,1.0D0, &
                 0.0D0,0.0D0,1.0D0, &
                 1.0D0,0.0D0,1.0D0],[3,6])
    
  INTEGER,PARAMETER:: nlmax2=2
  REAL(dp):: line_rgb2(3,nlmax2) &
       =RESHAPE( [0.D0,0.D0,1.D0, 1.D0,0.D0,0.D0],[3,2])
  INTEGER,PARAMETER:: line_pat2(nlmax2)=[ 0, 0 ]
  REAL(dp),PARAMETER:: line_mark_size2(nlmax2)=[ 0.3D0, 0.3D0 ]
  INTEGER,PARAMETER:: line_mark2(nlmax2)=[ 0, -1 ]
  INTEGER,PARAMETER:: line_mark_step2(nlmax2)=[ 1, 14 ]

CONTAINS

  SUBROUTINE cv_gout
    USE cvcomm
    USE cvparm,ONLY: cv_parm
    USE cvrank
    USE libcharx
    USE libkio
    IMPLICIT NONE
    CHARACTER(LEN=2):: kword
    CHARACTER(LEN=80):: line
    CHARACTER(LEN=2),ALLOCATABLE:: kworda(:)
    INTEGER:: mode
    CHARACTER(LEN=1):: kid
    INTEGER:: ncountry,ncountry_plot,id,idx,nword,nword_max
    INTEGER:: nl,nplot,nplot_max,idrank
    LOGICAL:: nplot_max_clear_logic
    INTEGER,PARAMETER:: nplot_m=32
    INTEGER:: ncountry_nplot(nplot_m)
    INTEGER,SAVE:: INIT=0
    EXTERNAL TOUPPER

    IF(INIT.EQ.0) THEN
       DO nl=1,nlmax
          line_rgb(1,nl)=line_rgb1(1,MOD(nl-1,nlmax1)+1)
          line_rgb(2,nl)=line_rgb1(2,MOD(nl-1,nlmax1)+1)
          line_rgb(3,nl)=line_rgb1(3,MOD(nl-1,nlmax1)+1)
          line_pat(nl)=0
          SELECT CASE((nl-1)/nlmax1+1)
          CASE(1)
             line_mark(nl)=-9
             line_mark_size(nl)=0.2D0
          CASE(2)
             line_mark(nl)=-6
             line_mark_size(nl)=0.2D0
          CASE(3)
             line_mark(nl)=-8
             line_mark_size(nl)=0.3D0
          CASE(4)
             line_mark(nl)=-3
             line_mark_size(nl)=0.2D0
          END SELECT
          line_mark_step(nl)=14
       END DO
       INIT=1
    END IF
    
    id=0
    nplot_max=0
    nplot_max_clear_logic=.TRUE.

1   CONTINUE
    
    WRITE(6,'(A,A)') &
         '## INPUT country id, graph id, ', &
         'rank:N1-6,R1-6, help:?,?H/C/N/R/G/L/D/P, end:X'
    CALL task_klin(line,kid,mode,cv_parm)
    IF(mode.EQ.2.OR.mode.EQ.3) GO TO 1  ! 2: parm input, 3:input error
    
    CALL ksplitn(line,' ,',2,nword_max,kworda)

    nplot_max_clear_logic=.TRUE.
    
    DO nword=1,nword_max
       kword=kworda(nword)
       IF(kword.EQ.'+ ') GO TO 2
       IF(kword(1:1).EQ.'-') GO TO 2
       READ(kword,'(I2)',ERR=2) idx ! positive number for graph id
!       WRITE(6,'(A,I6)') 'idx:',idx
       id=idx
       IF(id.EQ.9) GO TO 9000       ! 
       CYCLE

2      CONTINUE
       CALL TOUPPER(kword)

       SELECT CASE(kword)
       CASE('X ','XX')
          GO TO 9000
       CASE('?C')
          CALL cv_gout_help(1)
          id=-id
       CASE('?N')
          CALL cv_gout_help(2)
          id=-id
       CASE('?G')
          CALL cv_gout_help(3)
          id=-id
       CASE('?H')
          CALL cv_gout_help(4)
          id=-id
       CASE('?R')
          CALL cv_gout_help(5)
          id=-id
       CASE('?L')
          CALL cv_gout_help(6)
          id=-id
       CASE('?D')
          CALL cv_gout_help(7)
          id=-id
       CASE('?P','? ')
          WRITE(6,'(A,A)') &
               ' ncountry id  population      totC    newC    totD  newD'
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             WRITE(6,'(I3,I6,1X,A2,F12.3,I10,I8,I8,I6,2X,A20)') &
                  nplot,ncountry,country_id_ncountry(ncountry), &
                  population_ncountry(ncountry), &
                  ncases_total_ndate_ncountry(ndate_max,ncountry), &
                  ncases_new_ndate_ncountry(ndate_max,ncountry), &
                  ndeaths_total_ndate_ncountry(ndate_max,ncountry), &
                  ndeaths_new_ndate_ncountry(ndate_max,ncountry), &
                  TRIM(country_name_ncountry(ncountry))
          END DO
          id=-id
       CASE('+ ')
          nplot_max_clear_logic=.FALSE.
       CASE('- ','-1','-2','-3','-4','-5','-6','-7','-8','-9')
          SELECT CASE(kword)
          CASE('- ','-1')
             nplot_max=nplot_max-1
          CASE('-2')
             nplot_max=nplot_max-2
          CASE('-3')
             nplot_max=nplot_max-3
          CASE('-4')
             nplot_max=nplot_max-4
          CASE('-5')
             nplot_max=nplot_max-5
          CASE('-6')
             nplot_max=nplot_max-6
          CASE('-7')
             nplot_max=nplot_max-7
          CASE('-8')
             nplot_max=nplot_max-8
          CASE('-9')
             nplot_max=nplot_max-9
          END SELECT
          IF(nplot_max.LT.0) nplot_max=0
       CASE('N1','N2','N3','N4','N5','N6', &
            'R1','R2','R3','R4','R5','R6')
          SELECT CASE(kword)
          CASE('N1')
             idrank=1
          CASE('N2')
             idrank=2
          CASE('N3')
             idrank=3
          CASE('N4')
             idrank=4
          CASE('N5')
             idrank=5
          CASE('N6')
             idrank=6
          CASE('R1')
             idrank=7
          CASE('R2')
             idrank=8
          CASE('R3')
             idrank=9
          CASE('R4')
             idrank=10
          CASE('R5')
             idrank=11
          CASE('R6')
             idrank=12
          END SELECT
          CALL cv_rank_exec(idrank)
             
          nplot_max=nrank_max
          DO nplot=1,nplot_max
             ncountry=ncountry_nrank_idrank(nplot,idrank)
             ncountry_nplot(nplot)=ncountry
             WRITE(6,'(I4,2X,A2,2X,A)') &
                  nplot,country_id_ncountry(ncountry), &
                  TRIM(country_name_ncountry(ncountry))
          END DO
       CASE DEFAULT
          ncountry_plot=0
          DO ncountry=1,ncountry_max
             IF(kword.EQ.country_id_ncountry(ncountry)) THEN
                ncountry_plot=ncountry
                EXIT
             END IF
          END DO
          IF(ncountry_plot.EQ.0) THEN
             WRITE(6,'(A,A2)') 'XX cvgout: unknown country_id:',kword
          ELSE
             IF(nplot_max_clear_logic) THEN
                nplot_max=1
                nplot_max_clear_logic=.FALSE.
             ELSE
                IF(nplot_max.EQ.nplot_m) THEN
                   WRITE(6,'(A,I6)') &
                        'XX cvgout: nplot_max exceeds nplot_m:',nplot_m
                ELSE
                   nplot_max=nplot_max+1
                END IF
             END IF
             ncountry_nplot(nplot_max)=ncountry_plot
          END IF
       END SELECT
    END DO

    IF(id.LE.0) THEN
       id=-id
       GO TO 1
    END IF
    IF(nplot_max.EQ.0) GO TO 1
      
    SELECT CASE(id)
    CASE(1,10:19)
       CALL cv_gsub1(id,ncountry_nplot,nplot_max)
    CASE(2,20:29)
       CALL cv_gsub2(id,ncountry_nplot,nplot_max)
    CASE(3,30:39)
       CALL cv_gsub3(id,ncountry_nplot,nplot_max)
    CASE(4,40:49)
       CALL cv_gsub4(id,ncountry_nplot,nplot_max)
    CASE DEFAULT
       WRITE(6,'(A,I3)') 'XX cvgout: graph id out of range : id=',id
    END SELECT
    GO TO 1


9000 CONTINUE
    RETURN
  END SUBROUTINE cv_gout

  SUBROUTINE cv_gout_help(id)
    USE cvcomm
    USE cvlib
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id
    INTEGER:: ncountry,ndate
    CHARACTER(LEN=10):: date_id

    SELECT CASE(id)
    CASE(1) !XC
       WRITE(6,'(A)') '## List of country id:'
       WRITE(6,'(25(1X,A2))') &
            (country_id_ncountry(ncountry),ncountry=1,ncountry_max)
    CASE(2) !XN
       WRITE(6,'(A)') '## List of country name ane id:'
       WRITE(6,'(4(3X,A14,1X,A2))') &
            (country_name_ncountry(ncountry)//'              ', &
            country_id_ncountry(ncountry), &
            ncountry=1,ncountry_max)
    CASE(3) !XG
       WRITE(6,'(A)') '## List of graph id:'
       WRITE(6,'((A))') &
            ' 1: 10: 11+12+13+14', &
            '    11: total cases vs time', &
            '    12: total deaths vs time', &
            '    13: new cases vs time', &
            '    14: new deaths vs time', &
            '    15: 16+17+18+19', &
            '    16: total cases vs time (Log10 scale)', &
            '    17: total deaths vs time (Log10 scale)', &
            '    18: new cases vs time (Log10 scale)', &
            '    19: new deaths vs time (Log10 scale)', &
            ' 2: 20: 21+22+23+24 (per million population) ', &
            '    21: total cases per population vs time', &
            '    22: total deaths per population vs time', &
            '    23: new cases per population vs time', &
            '    24: new deaths per population vs time', &
            '    25: 26+27+28+29 (per million population) ', &
            '    26: total cases per population vs time (Log10 scale)', &
            '    27: total deaths per population vs time (Log10 scale)', &
            '    28: new cases per population vs time (Log10 scale)', &
            '    29: new deaths per population vs time (Log10 scale)', &
            ' 3: 30: 31+32+33+34 ', &
            '    31: total deaths number vs cases number (fixed range)', &
            '    32: total deaths rate vs cases rate (fixed range)', &
            '    33: total deaths number vs cases number (adjusted range)', &
            '    34: total deaths rate vs cases rate (adjusted range)', &
            ' 4  40: 41+42+43+44 ', &
            '    41: new deaths number vs cases number (fixed range)', &
            '    42: new deaths rate vs cases rate (fixed range)', &
            '    43: new deaths number vs cases number (adjusted range)', &
            '    44: new deaths rate vs cases rate (adjusted range)'
    CASE(4) !XH
       WRITE(6,'(A)')   '## Help: how to input counrty_ids and graph_id'
       WRITE(6,'(A)')   '      - input a list of country_ids or graph_id'
       WRITE(6,'(A,A)') '      - two-character country_ids should be ', &
                        'separated by space or comma.'
       WRITE(6,'(A)')   '      - the length of list is limited to one line.'
       WRITE(6,'(A,A)') '      - the order of country_ids determines ', &
                        'the line attribute.'
       WRITE(6,'(A,A)') '      - graph_id can be located on the same line ',&
                        'with a list of country_ids.'
       WRITE(6,'(A,A)') '      - if more than two graph_ids exist, ', &
                        'only the last one is effective.'
       WRITE(6,'(A,A)') '?H: this help'
       WRITE(6,'(A,A)') '?C: list of country id'
       WRITE(6,'(A,A)') '?N: list of country name'
       WRITE(6,'(A,A)') '?G: list of graph ids'
       WRITE(6,'(A,A)') '?R: list of ranking ids'
       WRITE(6,'(A,A)') '?L: list of line colors'
       WRITE(6,'(A,A)') '?D: convert ndate to data-id'
       WRITE(6,'(A,A)') '?P: list of countries to be plotted'
       WRITE(6,'(A,A)') '+ : addition of country (before any country id)'
       WRITE(6,'(A,A)') '-n: removal of countries from the last'
       WRITE(6,'(A,A)') 'ndays_ave=7                   ! ', &
            'range of day averageing'
       WRITE(6,'(A,A)') 'nrank_max=12                  ! ', &
            'number of items in ranking lists'
       WRITE(6,'(A,A)') 'ndate_min_g=1                 ! ', &
            'minimum date number in graph (1=2020-01-03) '
       WRITE(6,'(A,A)') 'ndate_max_g=ndate_max         ! ', &
            'maximum date number in graph                '
       WRITE(6,'(A,A)') 'population_min_rank=10.D0     ! ', &
            'minimum population per million in ranking list'
       WRITE(6,'(A,A)') 'cases_number_log_min=10.D0    ! ', &
            'minimum number for ncases in log'
       WRITE(6,'(A,A)') 'deaths_number_log_min=1.D0    ! ', &
            'minimum number for ndeaths in log'
       WRITE(6,'(A,A)') 'cases_rate_log_min=0.1D0      ! ', &
            'minimum rate for ncases in log'
       WRITE(6,'(A,A)') 'deaths_rate_log_min=0.01D0    ! ', &
            'minimum rate for ndeaths in log'
       WRITE(6,'(A,A)') 'ratio_new_total_log_min=0.1D0 ! ', &
            'ratio of log minimum between new and total'
       
    CASE(5) !XR
       WRITE(6,'(A)')   '## List of Ranking selection'
       WRITE(6,'(A)')   ' N1: Ranking by number of total cases'
       WRITE(6,'(A)')   ' N2: Ranking by number of new cases'
       WRITE(6,'(A)')   ' N3: Ranking by number of new cases 7 days average'
       WRITE(6,'(A)')   ' N4: Ranking by number of total deaths'
       WRITE(6,'(A)')   ' N5: Ranking by number of new deaths'
       WRITE(6,'(A)')   ' N6: Ranking by number of new deaths 7 days average'
       WRITE(6,'(A)')   ' R1: Ranking by rate of total cases'
       WRITE(6,'(A)')   ' R2: Ranking by rate of new cases'
       WRITE(6,'(A)')   ' R3: Ranking by rate of new cases 7 days average'
       WRITE(6,'(A)')   ' R4: Ranking by rate of total deaths'
       WRITE(6,'(A)')   ' R5: Ranking by rate of new deaths'
       WRITE(6,'(A)')   ' R6: Ranking by rate of new deaths 7 days average'
       WRITE(6,'(A)')   '                rate means per population in million'
    CASE(6) !XL
       WRITE(6,'(A)')   '## List of Line attribute'
       WRITE(6,'(A)')   ' 1: red        circulrar'
       WRITE(6,'(A)')   ' 2: yellow     circulrar'
       WRITE(6,'(A)')   ' 3: green      circulrar'
       WRITE(6,'(A)')   ' 4: light blue circulrar'
       WRITE(6,'(A)')   ' 5: dark blue  circulrar'
       WRITE(6,'(A)')   ' 6: purple     circulrar'
       WRITE(6,'(A)')   ' 7: red        triangle'
       WRITE(6,'(A)')   ' 8: yellow     triangle'
       WRITE(6,'(A)')   ' 9: green      triangle'
       WRITE(6,'(A)')   '10: light blue triangle'
       WRITE(6,'(A)')   '11: dark blue  triangle'
       WRITE(6,'(A)')   '12: purple     triangle'
    CASE(7) !XD
1      CONTINUE
       WRITE(6,'(A,I4,A)') '## input ndate: 1..',ndate_max,'  0: for end'
       READ(5,*,ERR=1,END=9) ndate
       IF(ndate.LE.0) GO TO 9
       IF(ndate.GT.ndate_max) GO TO 1
       CALL convert_ndate_to_date_id(ndate,date_id)
       WRITE(6,'(A,I6,A,A)') '  ndate=',ndate,'  date_id=',date_id
       GO TO 1
9      CONTINUE
    END SELECT
    RETURN
  END SUBROUTINE cv_gout_help

    ! plot graph by number
    
  SUBROUTINE cv_gsub1(id,ncountry_nplot,nplot_max)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)
    
    REAL(dp),ALLOCATABLE:: xg(:),fg(:,:)
    INTEGER:: ndata,ndate,ngid
    CHARACTER(LEN=80):: title
    EXTERNAL PAGES,PAGEE

    INTEGER:: mode_ls,id_ls,id_ltype,ncount,ncount_max

    ncount_max=ndate_max_g-ndate_min_g+1
    ALLOCATE(xg(ncount_max))
    DO ncount=1,ncount_max
       ndate=ncount+ndate_min_g-1
       xg(ncount)=DBLE(ndate)
    END DO

    SELECT CASE(id)
    CASE(1,10,11,12,13,14)
       mode_ls=0
       id_ls=0
       id_ltype=0
    CASE(15,16,17,18,19)
       mode_ls=2
       id_ls=5
       id_ltype=1
    END SELECT

    CALl PAGES
    SELECT CASE(id)
    CASE(1,10,11,15,16)
       ngid=1
       IF(id.EQ.11.OR.id.EQ.16) ngid=0
       CALL cv_gsub11(1+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size, &
                  NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
    END SELECT
    SELECT CASE(id)
    CASE(1,10,12,15,17)
       ngid=2
       IF(id.EQ.12.OR.id.EQ.17) ngid=0
       CALL cv_gsub11(2+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size, &
                  NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
    END SELECT
    SELECT CASE(id)
    CASE(1,10,13,15,18)
       ngid=3
       IF(id.EQ.13.OR.id.EQ.18) ngid=0
       CALL cv_gsub11(3+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       IF(nplot_max.EQ.1) THEN
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                     LINE_MARK_STEP=line_mark_step, &
                     LINE_MARK_SIZE=line_mark_size, &
                     NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       END IF
    END SELECT
    SELECT CASE(id)
    CASE(1,10,14,15,19)
       ngid=4
       IF(id.EQ.14.OR.id.EQ.19) ngid=0
       CALL cv_gsub11(4+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       IF(nplot_max.EQ.1) THEN
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                     LINE_MARK_STEP=line_mark_step, &
                     LINE_MARK_SIZE=line_mark_size, &
                     NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       END IF
    END SELECT
    CALL PAGEE

    RETURN
  END SUBROUTINE cv_gsub1

  ! --- set up data ---

  SUBROUTINE cv_gsub11(id,ncountry_nplot,nplot_max,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)
    REAL(dp),ALLOCATABLE,INTENT(OUT):: fg(:,:)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate,ndate_ave,nplot,ncountry,ncount,ncount_max
    REAL(dp):: rcases,rdeaths

    ncount_max=ndate_max_g-ndate_min_g+1

    SELECT CASE(id)
    CASE(1,6)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(ncount_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(1)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                fg(ncount,nplot)=ncases_total_ndate_ncountry(ndate,ncountry)
             END DO
          CASE(6)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                rcases=ncases_total_ndate_ncountry(ndate,ncountry)
                fg(ncount,nplot)=LOG10(MAX(cases_number_log_min,rcases))
             END DO
          END SELECT
       END DO
       IF(nplot_max.EQ.1) THEN
          ncountry=ncountry_nplot(1)
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'total cases@'
       ELSE
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'total cases@'
       END IF
    CASE(2,7)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(ncount_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(2)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                fg(ncount,nplot)=ndeaths_total_ndate_ncountry(ndate,ncountry)
             END DO
          CASE(7)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                rdeaths=ndeaths_total_ndate_ncountry(ndate,ncountry)
                fg(ncount,nplot)=LOG10(MAX(deaths_number_log_min,rdeaths))
             END DO
          END SELECT
       END DO
       IF(nplot_max.EQ.1) THEN
          ncountry=ncountry_nplot(1)
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'total deaths@'
       ELSE
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'total deaths@'
       END IF
    CASE(3,8)
       IF(nplot_max.EQ.1) THEN
          ndata=2
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(3)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,1)=ncases_new_ndate_ncountry(ndate,ncountry)
                   fg(ncount,2)=0.D0
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                fg(ncount,1)=ncases_new_ndate_ncountry(ndate,ncountry)
                fg(ncount,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ncount,2)=fg(ncount,2) &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ncount,2)=fg(ncount,2)/ndays_ave
             END DO
          CASE(8)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   rcases=ncases_new_ndate_ncountry(ndate,ncountry)
                   fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                         *cases_number_log_min,rcases))
                   fg(ncount,2)=LOG10(ratio_new_total_log_min &
                                         *cases_number_log_min)
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                rcases=ncases_new_ndate_ncountry(ndate,ncountry)
                fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                      *cases_number_log_min,rcases))
                rcases=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rcases=rcases &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ncount,2)=LOG10(MAX(ratio_new_total_log_min &
                                      *cases_number_log_min,rcases/ndays_ave))
             END DO
          END SELECT
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'new cases (ave)@'
       ELSE
          ndata=nplot_max
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          SELECT CASE(id)
          CASE(3)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=0.D0
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ncount,nplot)=fg(ncount,nplot) &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ncount,nplot)=fg(ncount,nplot)/ndays_ave
                END DO
             END DO
          CASE(8)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=LOG10(ratio_new_total_log_min &
                                            *cases_number_log_min)
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   rcases=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rcases=rcases &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                             *cases_number_log_min,rcases &
                                             /ndays_ave))
                END DO
             END DO
          END SELECT
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'new cases (ave)@'
       END IF
    CASE(4,9)
       IF(nplot_max.EQ.1) THEN
          ndata=2
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(4)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,1)=ndeaths_new_ndate_ncountry(ndate,ncountry)
                   fg(ncount,2)=0.D0
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                fg(ncount,1)=ndeaths_new_ndate_ncountry(ndate,ncountry)
                fg(ncount,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ncount,2)=fg(ncount,2) &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ncount,2)=fg(ncount,2)/ndays_ave
             END DO
          CASE(9)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry)
                   fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                         *deaths_number_log_min,rdeaths))
                   fg(ncount,2)=LOG10(ratio_new_total_log_min &
                                         *deaths_number_log_min)
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry)
                fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                      *deaths_number_log_min,rdeaths))
                rdeaths=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rdeaths=rdeaths &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ncount,2)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_number_log_min,rdeaths/ndays_ave))
             END DO
          END SELECT
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'new deaths (ave)@'
       ELSE
          ndata=nplot_max
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          SELECT CASE(id)
          CASE(4)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=0.D0
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ncount,nplot)=fg(ncount,nplot) &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ncount,nplot)=fg(ncount,nplot)/ndays_ave
                END DO
             END DO
          CASE(9)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=LOG10(ratio_new_total_log_min &
                                            *deaths_number_log_min)
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   rdeaths=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rdeaths=rdeaths &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ncount,nplot) &
                        =LOG10(MAX(ratio_new_total_log_min &
                                  *deaths_number_log_min,rdeaths &
                                  /ndays_ave))
                END DO
             END DO
          END SELECT
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'new deaths (ave)@'
       END IF
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub11
    
  ! plot graph by rate
    
  SUBROUTINE cv_gsub2(id,ncountry_nplot,nplot_max)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)
    
    REAL(dp),ALLOCATABLE:: xg(:),fg(:,:)
    INTEGER:: ndata,ndate,ngid
    CHARACTER(LEN=80):: title

    INTEGER:: mode_ls,id_ls,id_ltype,ncount,ncount_max
    EXTERNAL PAGES,PAGEE

    ncount_max=ndate_max_g-ndate_min_g+1
    ALLOCATE(xg(ncount_max))
    DO ncount=1,ncount_max
       ndate=ncount+ndate_min_g-1
       xg(ncount)=DBLE(ndate)
    END DO

    SELECT CASE(id)
    CASE(2,20,21,22,23,24)
       mode_ls=0
       id_ls=0
       id_ltype=0
    CASE(25,26,27,28,29)
       mode_ls=2
       id_ls=5
       id_ltype=1
    END SELECT

    CALl PAGES
    SELECT CASE(id)
    CASE(2,20,21,25,26)
       ngid=1
       IF(id.EQ.21.OR.id.EQ.26) ngid=0
       CALL cv_gsub21(1+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size, &
                  NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
    END SELECT
    SELECT CASE(id)
    CASE(2,20,22,25,27)
       ngid=2
       IF(id.EQ.22.OR.id.EQ.27) ngid=0
       CALL cv_gsub21(2+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size, &
                  NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
    END SELECT
    SELECT CASE(id)
    CASE(2,20,23,25,28)
       ngid=3
       IF(id.EQ.23.OR.id.EQ.28) ngid=0
       CALL cv_gsub21(3+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       IF(nplot_max.EQ.1) THEN
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                     LINE_MARK_STEP=line_mark_step, &
                     LINE_MARK_SIZE=line_mark_size, &
                     NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       END IF
    END SELECT
    SELECT CASE(id)
    CASE(2,20,24,25,29)
       ngid=4
       IF(id.EQ.24.OR.id.EQ.29) ngid=0
       CALL cv_gsub21(4+id_ls,ncountry_nplot,nplot_max,fg,title,ndata)
       IF(nplot_max.EQ.1) THEN
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ncount_max,ncount_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                     LINE_MARK_STEP=line_mark_step, &
                     LINE_MARK_SIZE=line_mark_size, &
                     NLMAX=nlmax,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       END IF
    END SELECT
    CALL PAGEE

    RETURN
  END SUBROUTINE cv_gsub2

  ! --- set up data ---

  SUBROUTINE cv_gsub21(id,ncountry_nplot,nplot_max,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)
    REAL(dp),ALLOCATABLE,INTENT(OUT):: fg(:,:)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate,ndate_ave,nplot,ncountry,ncount,ncount_max
    REAL(dp):: rcases,rdeaths

    ncount_max=ndate_max_g-ndate_min_g+1

    SELECT CASE(id)
    CASE(1,6)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(ncount_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(1)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                fg(ncount,nplot)=ncases_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
             END DO
          CASE(6)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                rcases=ncases_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ncount,nplot)=LOG10(MAX(cases_rate_log_min,rcases))
             END DO
          END SELECT
       END DO
       IF(nplot_max.EQ.1) THEN
          ncountry=ncountry_nplot(1)
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'total cases/pop(M)@'
       ELSE
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'total cases/pop(M)@'
       END IF
    CASE(2,7)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(ncount_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(2)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                fg(ncount,nplot)=ndeaths_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
             END DO
          CASE(7)
             DO ncount=1,ncount_max
                ndate=ndate_min_g+ncount-1
                rdeaths=ndeaths_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ncount,nplot)=LOG10(MAX(deaths_rate_log_min,rdeaths))
             END DO
          END SELECT
       END DO
       IF(nplot_max.EQ.1) THEN
          ncountry=ncountry_nplot(1)
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'total deaths/pop(M)@'
       ELSE
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'total deaths/pop(M)@'
       END IF
    CASE(3,8)
       IF(nplot_max.EQ.1) THEN
          ndata=2
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(3)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,1)=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                   fg(ncount,2)=0.D0
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                fg(ncount,1)=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ncount,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ncount,2)=fg(ncount,2) &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ncount,2)=fg(ncount,2)/ndays_ave &
                     /population_ncountry(ncountry)
             END DO
          CASE(8)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   rcases=ncases_new_ndate_ncountry(ndate,ncountry) &
                        /population_ncountry(ncountry)
                   fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                         *cases_rate_log_min,rcases))
                   fg(ncount,2)=LOG10(ratio_new_total_log_min &
                                         *cases_rate_log_min)
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                rcases=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                      *cases_rate_log_min,rcases))
                rcases=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rcases=rcases &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                rcases=rcases/population_ncountry(ncountry)
                fg(ncount,2)=LOG10(MAX(ratio_new_total_log_min &
                                      *cases_rate_log_min,rcases/ndays_ave))
             END DO
          END SELECT
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'new cases/pop(M) (ave)@'
       ELSE
          ndata=nplot_max
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          SELECT CASE(id)
          CASE(3)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=0.D0
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ncount,nplot)=fg(ncount,nplot) &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ncount,nplot)=fg(ncount,nplot)/ndays_ave &
                        /population_ncountry(ncountry)
                END DO
             END DO
          CASE(8)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=LOG10(ratio_new_total_log_min &
                                            *cases_rate_log_min)
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   rcases=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rcases=rcases &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   rcases=rcases/population_ncountry(ncountry)
                   fg(ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                             *cases_rate_log_min,rcases &
                                             /ndays_ave))
                END DO
             END DO
          END SELECT
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'new cases/pop(M) (ave)@'
       END IF
    CASE(4,9)
       IF(nplot_max.EQ.1) THEN
          ndata=2
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(4)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,1)=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                   fg(ncount,2)=0.D0
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                fg(ncount,1)=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ncount,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ncount,2)=fg(ncount,2) &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ncount,2)=fg(ncount,2)/ndays_ave &
                     /population_ncountry(ncountry)
             END DO
          CASE(9)
             IF(ndate_min_g.LE.ndays_ave-1) THEN
                DO ncount=1,ndays_ave-ndate_min_g
                   ndate=ncount+ndate_min_g-1
                   rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                   fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                         *deaths_rate_log_min,rdeaths))
                   fg(ncount,2)=LOG10(ratio_new_total_log_min &
                                         *deaths_rate_log_min)
                END DO
             END IF
             DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                ndate=ncount+ndate_min_g-1
                rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ncount,1)=LOG10(MAX(ratio_new_total_log_min &
                                      *deaths_rate_log_min,rdeaths))
                rdeaths=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rdeaths=rdeaths &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                rdeaths=rdeaths/population_ncountry(ncountry)
                fg(ncount,2)=LOG10(MAX(ratio_new_total_log_min &
                                      *deaths_rate_log_min,rdeaths/ndays_ave))
             END DO
          END SELECT
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'new deaths/pop(M) (ave)@'
       ELSE
          ndata=nplot_max
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ncount_max,ndata))
          SELECT CASE(id)
          CASE(4)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=0.D0
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   fg(ncount,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ncount,nplot)=fg(ncount,nplot) &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ncount,nplot)=fg(ncount,nplot)/ndays_ave &
                        /population_ncountry(ncountry)
                END DO
             END DO
          CASE(9)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                IF(ndate_min_g.LE.ndays_ave-1) THEN
                   DO ncount=1,ndays_ave-ndate_min_g
                      fg(ncount,nplot)=LOG10(ratio_new_total_log_min &
                                            *deaths_rate_log_min)
                   END DO
                END IF
                DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
                   ndate=ncount+ndate_min_g-1
                   rdeaths=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rdeaths=rdeaths &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   rdeaths=rdeaths/population_ncountry(ncountry)
                   fg(ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                             *deaths_rate_log_min,rdeaths &
                                             /ndays_ave))
                END DO
             END DO
          END SELECT
          title=''
          DO nplot=1,nplot_max
             ncountry=ncountry_nplot(nplot)
             title=TRIM(title)//' '//country_id_ncountry(ncountry)
          END DO
          title='@'//title(2:LEN_TRIM(title))//': '// &
               'new deaths/pop(M) (ave)@'
       END IF
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub21
    
  ! plot graph of total ratio
    
  SUBROUTINE cv_gsub3(id,ncountry_nplot,nplot_max)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)

    REAL(dp),ALLOCATABLE:: fg(:,:,:)
    INTEGER,ALLOCATABLE:: ncount_maxa(:)
    INTEGER:: ndata,ngid,ncount_max
    CHARACTER(LEN=80):: title
    EXTERNAL PAGES,PAGEE

    
    ncount_max=ndate_max_g-ndate_min_g+1

    ALLOCATE(ncount_maxa(nplot_max))
    ncount_maxa(1:nplot_max)=ncount_max

    CALL PAGES
    SELECT CASE(id)
    CASE(3,30,31)
       ngid=1
       IF(id.EQ.31) ngid=0
       CALL cv_gsub31(1,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=3.D0,XMAX=7.477D0,YMIN=1.D0,YMAX=5.477D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    SELECT CASE(id)
    CASE(3,30,32)
       ngid=2
       IF(id.EQ.32) ngid=0
       CALL cv_gsub31(2,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=1.D0,XMAX=5.D0,YMIN=-0.523D0,YMAX=3.477D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    SELECT CASE(id)
    CASE(3,30,33)
       ngid=3
       IF(id.EQ.33) ngid=0
       CALL cv_gsub31(3,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=LOG10(cases_number_log_min), &
                  YMIN=LOG10(deaths_number_log_min), &
                  ASPECT=0.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    SELECT CASE(id)
    CASE(3,30,34)
       ngid=4
       IF(id.EQ.34) ngid=0
       CALL cv_gsub31(4,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=LOG10(cases_rate_log_min), &
                  YMIN=LOG10(deaths_rate_log_min), &
                  ASPECT=0.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    CALL PAGEE

    RETURN
  END SUBROUTINE cv_gsub3

  ! --- set up data '''

  SUBROUTINE cv_gsub31(id,ncountry_nplot,nplot_max,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)
    REAL(dp),ALLOCATABLE,INTENT(OUT):: fg(:,:,:)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate,nplot,ncountry,ncount,ncount_max
    REAL(dp):: rcases,rdeaths

    ncount_max=ndate_max_g-ndate_min_g+1

    SELECT CASE(id)
    CASE(1,3)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(2,ncount_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          DO ncount=1,ncount_max
             ndate=ncount+ndate_min_g-1
             rcases=DBLE(ncases_total_ndate_ncountry(ndate,ncountry))
             rdeaths=DBLE(ndeaths_total_ndate_ncountry(ndate,ncountry))
             fg(1,ncount,nplot)=LOG10(MAX(cases_number_log_min,rcases))
             fg(2,ncount,nplot)=LOG10(MAX(deaths_number_log_min,rdeaths))
          END DO
       END DO
       title=''
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          title=TRIM(title)//' '//country_id_ncountry(ncountry)
       END DO
       title='@'//title(2:LEN_TRIM(title))//': '// &
            'total D vs C in number@'
    CASE(2,4)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(2,ncount_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          DO ncount=1,ncount_max
             ndate=ncount+ndate_min_g-1
             rcases=DBLE(ncases_total_ndate_ncountry(ndate,ncountry)) &
                  /population_ncountry(ncountry)
             rdeaths=DBLE(ndeaths_total_ndate_ncountry(ndate,ncountry)) &
                  /population_ncountry(ncountry)
             fg(1,ncount,nplot)=LOG10(MAX(cases_rate_log_min,rcases))
             fg(2,ncount,nplot)=LOG10(MAX(deaths_rate_log_min,rdeaths))
          END DO
       END DO
       title=''
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          title=TRIM(title)//' '//country_id_ncountry(ncountry)
       END DO
       title='@'//title(2:LEN_TRIM(title))//': '// &
            'total D vs C in rate@'
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub31
    
  ! plot graph of total ratio
    
  SUBROUTINE cv_gsub4(id,ncountry_nplot,nplot_max)

    USE cvcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)

    REAL(dp),ALLOCATABLE:: fg(:,:,:)
    INTEGER,ALLOCATABLE:: ncount_maxa(:)
    INTEGER:: ndata,ngid,ncount_max,mode_2d
    CHARACTER(LEN=80):: title
    EXTERNAL PAGES,PAGEE

    IF(nplot_max.EQ.1) THEN
       mode_2d=22
    ELSE
       mode_2d=21
    END IF
    
    ncount_max=ndate_max_g-ndate_min_g+1

    ALLOCATE(ncount_maxa(nplot_max))
    ncount_maxa(1:nplot_max)=ncount_max

    CALL PAGES
    SELECT CASE(id)
    CASE(4,40,41)
       ngid=1
       IF(id.EQ.41) ngid=0
       CALL cv_gsub41(1,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=0.477D0,XMAX=5.477D0,YMIN=-1.D0,YMAX=4.D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1,MODE_2D=mode_2d, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    SELECT CASE(id)
    CASE(4,40,42)
       ngid=2
       IF(id.EQ.42) ngid=0
       CALL cv_gsub41(2,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=-1.523D0,XMAX=3.477D0,YMIN=-3.D0,YMAX=2.D0, &
                  ASPECT=1.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1,MODE_2D=mode_2d, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    SELECT CASE(id)
    CASE(4,40,43)
       ngid=3
       IF(id.EQ.43) ngid=0
       CALL cv_gsub41(3,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=LOG10(ratio_new_total_log_min*cases_number_log_min), &
                  YMIN=LOG10(ratio_new_total_log_min*deaths_number_log_min), &
                  ASPECT=0.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1,MODE_2D=mode_2d, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    SELECT CASE(id)
    CASE(4,40,44)
       ngid=4
       IF(id.EQ.44) ngid=0
       CALL cv_gsub41(4,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ncount_max,ncount_maxa,ndata,title, &
                  XMIN=LOG10(ratio_new_total_log_min*cases_rate_log_min), &
                  YMIN=LOG10(ratio_new_total_log_min*deaths_rate_log_min), &
                  ASPECT=0.D0,XSCALE_ZERO=0,YSCALE_ZERO=0,MODE_LS=3, &
                  XGRID_LTYPE=1,YGRID_LTYPE=1, &
                  NOINFO=1,MODE_2D=mode_2d, &
                  LINE_RGB=line_rgb,LINE_PAT=line_pat,LINE_MARK=line_mark, &
                  LINE_MARK_STEP=line_mark_step, &
                  LINE_MARK_SIZE=line_mark_size,NLMAX=nlmax)
    END SELECT
    CALL PAGEE

    RETURN
  END SUBROUTINE cv_gsub4

  ! --- set up data '''

  SUBROUTINE cv_gsub41(id,ncountry_nplot,nplot_max,fg,title,ndata)
    USE cvcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id,nplot_max
    INTEGER,INTENT(IN):: ncountry_nplot(nplot_max)
    REAL(dp),ALLOCATABLE,INTENT(OUT):: fg(:,:,:)
    CHARACTER(LEN=80),INTENT(OUT):: title
    INTEGER,INTENT(OUT):: ndata
    INTEGER:: ndate,nplot,ncountry,ndate_ave,ncount,ncount_max
    REAL(dp):: rcases,rdeaths

    ncount_max=ndate_max_g-ndate_min_g+1
    SELECT CASE(id)
    CASE(1,3)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(2,ncount_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          IF(ndate_min_g.LE.ndays_ave-1) THEN
             DO ncount=1,ndays_ave-ndate_min_g
                ndate=ncount+ndate_min_g-1
                fg(1,ncount,nplot)=LOG10(ratio_new_total_log_min &
                                        *cases_number_log_min)
                fg(2,ncount,nplot)=LOG10(ratio_new_total_log_min &
                                       *deaths_number_log_min)
             END DO
          END IF
          DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
             ndate=ncount+ndate_min_g-1
             rcases=0.D0
             rdeaths=0.D0
             DO ndate_ave=ndate-ndays_ave+1,ndate
                rcases=rcases &
                     +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                rdeaths=rdeaths &
                     +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
             END DO
             rcases=rcases/ndays_ave
             rdeaths=rdeaths/ndays_ave
             fg(1,ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                         *cases_number_log_min,rcases))
             fg(2,ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                         *deaths_number_log_min,rdeaths))
          END DO
       END DO
       title=''
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          title=TRIM(title)//' '//country_id_ncountry(ncountry)
       END DO
       title='@'//title(2:LEN_TRIM(title))//': '// &
            'new D vs C in number@'
    CASE(2,4)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(2,ncount_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          IF(ndate_min_g.LE.ndays_ave-1) THEN
             DO ncount=1,ndays_ave-ndate_min_g
                ndate=ncount+ndate_min_g-1
                fg(1,ncount,nplot)=LOG10(ratio_new_total_log_min &
                                        *cases_rate_log_min)
                fg(2,ncount,nplot)=LOG10(ratio_new_total_log_min &
                                        *deaths_rate_log_min)
             END DO
          END IF
          DO ncount=MAX(1,ndays_ave-ndate_min_g+1),ncount_max
             ndate=ncount+ndate_min_g-1
             rcases=0.D0
             rdeaths=0.D0
             DO ndate_ave=ndate-ndays_ave+1,ndate
                rcases=rcases &
                     +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                rdeaths=rdeaths &
                     +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
             END DO
             rcases=rcases/(population_ncountry(ncountry)*ndays_ave)
             rdeaths=rdeaths/(population_ncountry(ncountry)*ndays_ave)
             fg(1,ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                         *cases_rate_log_min,rcases))
             fg(2,ncount,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                         *deaths_rate_log_min,rdeaths))
          END DO
       END DO
       title=''
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          title=TRIM(title)//' '//country_id_ncountry(ncountry)
       END DO
       title='@'//title(2:LEN_TRIM(title))//': '// &
            'new D vs C in rate@'
    END SELECT

    RETURN
  END SUBROUTINE cv_gsub41
    
END MODULE cvgout
