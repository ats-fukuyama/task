! cvgout.f90

MODULE cvgout

  PRIVATE
  PUBLIC cv_gout

  INTEGER,PARAMETER:: nlmax=12
  REAL(8),PARAMETER:: line_rgb(3,nlmax) &
       =RESHAPE([1.0D0,0.0D0,0.0D0, &
                 1.0D0,0.8D0,0.0D0, &
                 0.0D0,0.8D0,0.0D0, &
                 0.0D0,0.8D0,1.0D0, &
                 0.0D0,0.0D0,1.0D0, &
                 1.0D0,0.0D0,1.0D0, &
                 1.0D0,0.0D0,0.0D0, &
                 1.0D0,0.8D0,0.0D0, &
                 0.0D0,0.8D0,0.0D0, &
                 0.0D0,0.8D0,1.0D0, &
                 0.0D0,0.0D0,1.0D0, &
                 1.0D0,0.0D0,1.0D0],[3,12])
  INTEGER,PARAMETER:: line_pat(nlmax) &
       =[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
  REAL(8),PARAMETER:: line_mark_size(nlmax) &
       =[ 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, &
          0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0 ]
  INTEGER,PARAMETER:: line_mark(nlmax) &
       =[ -1, -1, -1, -1, -1, -1, -4, -4, -4, -4, -4, -4 ]
  INTEGER,PARAMETER:: line_mark_step(nlmax) &
       =[ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 ]
    
  INTEGER,PARAMETER:: nlmax2=2
  REAL(8):: line_rgb2(3,nlmax2) &
       =RESHAPE( [0.D0,0.D0,1.D0, 1.D0,0.D0,0.D0],[3,2])
  INTEGER,PARAMETER:: line_pat2(nlmax2)=[ 0, 0 ]
  REAL(8),PARAMETER:: line_mark_size2(nlmax2)=[ 0.3D0, 0.3D0 ]
  INTEGER,PARAMETER:: line_mark2(nlmax2)=[ 0, -1 ]
  INTEGER,PARAMETER:: line_mark_step2(nlmax2)=[ 1, 14 ]

CONTAINS

  SUBROUTINE cv_gout
    USE cvcomm
    USE libcharx
    IMPLICIT NONE
    CHARACTER(LEN=2):: country_id,kword
    CHARACTER(LEN=80):: line
    CHARACTER(LEN=2),ALLOCATABLE:: kworda(:)
    INTEGER:: ncountry,ncountry_plot,id,idx,nword,nword_max,nplot,nplot_max
    INTEGER,ALLOCATABLE:: ncountry_nword(:),ncountry_nplot(:)

    id=0
    nplot_max=0

1   CONTINUE
    
    WRITE(6,'(A,A)') &
         '## INPUT country ids or graph id, ', &
         'list of ids:XC,XG, help:XH, end:XX or 9'
    READ(5,'(A)',ERR=1,END=9000) line
    CALL ksplitn(line,' ,',2,nword_max,kworda)

    IF(ALLOCATED(ncountry_nword)) DEALLOCATE(ncountry_nword)
    ALLOCATE(ncountry_nword(nword_max))

    nplot=0
    DO nword=1,nword_max
       kword=kworda(nword)
       READ(kword,'(I2)',ERR=2) idx
       id=idx
       IF(id.EQ.9) GO TO 9000
       CYCLE

2      CONTINUE
       CALL TOUPPER(kword)

       IF(kword.EQ.'XX') GO TO 9000
       IF(kword.EQ.'XC') THEN
          CALL cv_gout_help(1)
          CYCLE
       END IF
       IF(kword.EQ.'XG') THEN
          CALL cv_gout_help(2)
          CYCLE
       END IF
       IF(kword.EQ.'XH') THEN
          CALL cv_gout_help(3)
          CYCLE
       END IF

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
          nplot=nplot+1
          ncountry_nword(nplot)=ncountry_plot
       END IF
    END DO
    IF(nplot.GE.1) THEN
       nplot_max=nplot

       IF(ALLOCATED(ncountry_nplot)) DEALLOCATE(ncountry_nplot)
       ALLOCATE(ncountry_nplot(nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nword(nplot)
          ncountry_nplot(nplot)=ncountry
          WRITE(6,'(I4,2X,A2,2X,A)') &
               nplot,country_id_ncountry(ncountry), &
                     TRIM(country_name_ncountry(ncountry))
       END DO
    END IF

    IF(id.EQ.0) GO TO 1
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
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id
    INTEGER:: ncountry

    SELECT CASE(id)
    CASE(1) !XC
       WRITE(6,'(A)') '## List of country id:'
       WRITE(6,'(25(1X,A2))') &
            (country_id_ncountry(ncountry),ncountry=1,ncountry_max)
    CASE(2) !XG
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
            ' 2: 10: 21+22+23+24 (per million population) ', &
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
    CASE(3) !XH
       WRITE(6,'(A)') '## Help: how to input counrty_ids and graph_id'
       WRITE(6,'(A)') '    - input a list of country_ids or graph_id'
       WRITE(6,'(A)') '      - two-character country_ids should be separated '
       WRITE(6,'(A)') '        by space or comma.'
       WRITE(6,'(A)') '      - the length of list is limited to one line.'
       WRITE(6,'(A)') '      - the order of country_ids determines'
       WRITE(6,'(A)') '        the line attribute.'
       WRITE(6,'(A)') '      - graph_id can be located on the same line'
       WRITE(6,'(A)') '        with a list of country_ids.'
       WRITE(6,'(A)') '      - if there are more than two graph_ids,'
       WRITE(6,'(A)') '        only the last graph_id is effective.'
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
    
    REAL(dp):: xg(ndate_max)
    REAL(dp),ALLOCATABLE:: fg(:,:)
    INTEGER:: ndata,ndate,ngid
    CHARACTER(LEN=80):: title

    INTEGER:: mode_ls,id_ls,id_ltype

    DO ndate=1,ndate_max
       xg(ndate)=DBLE(ndate)
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
       CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
       CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
    INTEGER:: ndate,ndate_ave,nplot,ncountry
    REAL(dp):: rcases,rdeaths

    SELECT CASE(id)
    CASE(1,6)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(ndate_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(1)
             DO ndate=1,ndate_max
                fg(ndate,nplot)=ncases_total_ndate_ncountry(ndate,ncountry)
             END DO
          CASE(6)
             DO ndate=1,ndate_max
                rcases=ncases_total_ndate_ncountry(ndate,ncountry)
                fg(ndate,nplot)=LOG10(MAX(cases_number_log_min,rcases))
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
       ALLOCATE(fg(ndate_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(2)
             DO ndate=1,ndate_max
                fg(ndate,nplot)=ndeaths_total_ndate_ncountry(ndate,ncountry)
             END DO
          CASE(7)
             DO ndate=1,ndate_max
                rdeaths=ndeaths_total_ndate_ncountry(ndate,ncountry)
                fg(ndate,nplot)=LOG10(MAX(deaths_number_log_min,rdeaths))
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
          ALLOCATE(fg(ndate_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(3)
             DO ndate=1,ndays_ave-1
                fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,2)=0.D0
             END DO
             DO ndate=ndays_ave,ndate_max
                fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ndate,2)=fg(ndate,2) &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ndate,2)=fg(ndate,2)/ndays_ave
             END DO
          CASE(8)
             DO ndate=1,ndays_ave-1
                rcases=ncases_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *cases_number_log_min,rcases))
                fg(ndate,2)=LOG10(ratio_new_total_log_min &
                                 *cases_number_log_min)
             END DO
             DO ndate=ndays_ave,ndate_max
                rcases=ncases_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *cases_number_log_min,rcases))
                rcases=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rcases=rcases &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ndate,2)=LOG10(MAX(ratio_new_total_log_min &
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
          ALLOCATE(fg(ndate_max,ndata))
          SELECT CASE(id)
          CASE(3)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=0.D0
                END DO
                DO ndate=ndays_ave,ndate_max
                   fg(ndate,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ndate,nplot)=fg(ndate,nplot) &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ndate,nplot)=fg(ndate,nplot)/ndays_ave
                END DO
             END DO
          CASE(8)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=LOG10(ratio_new_total_log_min &
                                        *cases_number_log_min)
                END DO
                DO ndate=ndays_ave,ndate_max
                   rcases=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rcases=rcases &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
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
          ALLOCATE(fg(ndate_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(4)
             DO ndate=1,ndays_ave-1
                fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,2)=0.D0
             END DO
             DO ndate=ndays_ave,ndate_max
                fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ndate,2)=fg(ndate,2) &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ndate,2)=fg(ndate,2)/ndays_ave
             END DO
          CASE(9)
             DO ndate=1,ndays_ave-1
                rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_number_log_min,rdeaths))
                fg(ndate,2)=LOG10(deaths_number_log_min)
             END DO
             DO ndate=ndays_ave,ndate_max
                rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_number_log_min,rdeaths))
                rdeaths=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rdeaths=rdeaths &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ndate,2)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_number_log_min,rdeaths &
                                     /ndays_ave))
             END DO
          END SELECT
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'new deaths (ave)@'
       ELSE
          ndata=nplot_max
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ndate_max,ndata))
          SELECT CASE(id)
          CASE(4)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=0.D0
                END DO
                DO ndate=ndays_ave,ndate_max
                   fg(ndate,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ndate,nplot)=fg(ndate,nplot) &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ndate,nplot)=fg(ndate,nplot)/ndays_ave
                END DO
             END DO
          CASE(9)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=LOG10(ratio_new_total_log_min &
                                        *deaths_number_log_min)
                END DO
                DO ndate=ndays_ave,ndate_max
                   rdeaths=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rdeaths=rdeaths &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ndate,nplot) &
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
    
    REAL(dp):: xg(ndate_max)
    REAL(dp),ALLOCATABLE:: fg(:,:)
    INTEGER:: ndata,ndate,ngid
    CHARACTER(LEN=80):: title

    INTEGER:: mode_ls,id_ls,id_ltype

    DO ndate=1,ndate_max
       xg(ndate)=DBLE(ndate)
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
       CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
       CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
                     NOINFO=1, &
                     LINE_RGB=line_rgb2,LINE_PAT=line_pat2, &
                     LINE_MARK=line_mark2, &
                     LINE_MARK_STEP=line_mark_step2, &
                     LINE_MARK_SIZE=line_mark_size2, &
                     NLMAX=nlmax2,MODE_LS=mode_ls,FGRID_LTYPE=id_ltype)
       ELSE
          CALL grd1d(ngid,xg,fg,ndate_max,ndate_max,ndata,title, &
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
    INTEGER:: ndate,ndate_ave,nplot,ncountry
    REAL(dp):: rcases,rdeaths

    SELECT CASE(id)
    CASE(1,6)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(ndate_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(1)
             DO ndate=1,ndate_max
                fg(ndate,nplot)=ncases_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
             END DO
          CASE(6)
             DO ndate=1,ndate_max
                rcases=ncases_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,nplot)=LOG10(MAX(cases_rate_log_min,rcases))
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
       ALLOCATE(fg(ndate_max,ndata))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          SELECT CASE(id)
          CASE(2)
             DO ndate=1,ndate_max
                fg(ndate,nplot)=ndeaths_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
             END DO
          CASE(7)
             DO ndate=1,ndate_max
                rdeaths=ndeaths_total_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,nplot)=LOG10(MAX(deaths_rate_log_min,rdeaths))
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
          ALLOCATE(fg(ndate_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(3)
             DO ndate=1,ndays_ave-1
                fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,2)=0.D0
             END DO
             DO ndate=ndays_ave,ndate_max
                fg(ndate,1)=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ndate,2)=fg(ndate,2) &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ndate,2)=fg(ndate,2)/ndays_ave &
                     /population_ncountry(ncountry)
             END DO
          CASE(8)
             DO ndate=1,ndays_ave-1
                rcases=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *cases_rate_log_min,rcases))
                fg(ndate,2)=LOG10(ratio_new_total_log_min &
                                 *cases_rate_log_min)
             END DO
             DO ndate=ndays_ave,ndate_max
                rcases=ncases_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *cases_rate_log_min,rcases))
                rcases=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rcases=rcases &
                        +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                rcases=rcases/population_ncountry(ncountry)
                fg(ndate,2)=LOG10(MAX(ratio_new_total_log_min &
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
          ALLOCATE(fg(ndate_max,ndata))
          SELECT CASE(id)
          CASE(3)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=0.D0
                END DO
                DO ndate=ndays_ave,ndate_max
                   fg(ndate,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ndate,nplot)=fg(ndate,nplot) &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ndate,nplot)=fg(ndate,nplot)/ndays_ave &
                        /population_ncountry(ncountry)
                END DO
             END DO
          CASE(8)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=LOG10(ratio_new_total_log_min &
                                        *cases_rate_log_min)
                END DO
                DO ndate=ndays_ave,ndate_max
                   rcases=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rcases=rcases &
                           +ncases_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   rcases=rcases/population_ncountry(ncountry)
                   fg(ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
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
          ALLOCATE(fg(ndate_max,ndata))
          ncountry=ncountry_nplot(1)
          SELECT CASE(id)
          CASE(4)
             DO ndate=1,ndays_ave-1
                fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,2)=0.D0
             END DO
             DO ndate=ndays_ave,ndate_max
                fg(ndate,1)=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,2)=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   fg(ndate,2)=fg(ndate,2) &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                fg(ndate,2)=fg(ndate,2)/ndays_ave &
                     /population_ncountry(ncountry)
             END DO
          CASE(9)
             DO ndate=1,ndays_ave-1
                rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_rate_log_min,rdeaths))
                fg(ndate,2)=LOG10(ratio_new_total_log_min &
                                 *deaths_rate_log_min)
             END DO
             DO ndate=ndays_ave,ndate_max
                rdeaths=ndeaths_new_ndate_ncountry(ndate,ncountry) &
                     /population_ncountry(ncountry)
                fg(ndate,1)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_rate_log_min,rdeaths))
                rdeaths=0.D0
                DO ndate_ave=ndate-ndays_ave+1,ndate
                   rdeaths=rdeaths &
                        +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                END DO
                rdeaths=rdeaths/population_ncountry(ncountry)
                fg(ndate,2)=LOG10(MAX(ratio_new_total_log_min &
                                     *deaths_rate_log_min,rdeaths &
                                     /ndays_ave))
             END DO
          END SELECT
          title='@'//country_id_ncountry(ncountry)//': '// &
               TRIM(country_name_ncountry(ncountry))//': '// &
               region_id_ncountry(ncountry)//': '// &
               'new deaths/pop(M) (ave)@'
       ELSE
          ndata=nplot_max
          IF(ALLOCATED(fg)) DEALLOCATE(fg)
          ALLOCATE(fg(ndate_max,ndata))
          SELECT CASE(id)
          CASE(4)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=0.D0
                END DO
                DO ndate=ndays_ave,ndate_max
                   fg(ndate,nplot)=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      fg(ndate,nplot)=fg(ndate,nplot) &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   fg(ndate,nplot)=fg(ndate,nplot)/ndays_ave &
                        /population_ncountry(ncountry)
                END DO
             END DO
          CASE(9)
             DO nplot=1,nplot_max
                ncountry=ncountry_nplot(nplot)
                DO ndate=1,ndays_ave-1
                   fg(ndate,nplot)=LOG10(ratio_new_total_log_min &
                                         *deaths_rate_log_min)
                END DO
                DO ndate=ndays_ave,ndate_max
                   rdeaths=0.D0
                   DO ndate_ave=ndate-ndays_ave+1,ndate
                      rdeaths=rdeaths &
                           +ndeaths_new_ndate_ncountry(ndate_ave,ncountry)
                   END DO
                   rdeaths=rdeaths/population_ncountry(ncountry)
                   fg(ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
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

    REAL(dp):: xg(ndate_max)
    REAL(dp),ALLOCATABLE:: fg(:,:,:)
    INTEGER,ALLOCATABLE:: ndate_maxa(:)
    INTEGER:: ndata,ndate,ngid
    CHARACTER(LEN=80):: title

    
    ALLOCATE(ndate_maxa(nplot_max))
    ndate_maxa(1:nplot_max)=ndate_max

    CALL PAGES
    SELECT CASE(id)
    CASE(3,30,31)
       ngid=1
       IF(id.EQ.31) ngid=0
       CALL cv_gsub31(1,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
                  XMIN=2.477D0,XMAX=7.D0,YMIN=1.D0,YMAX=5.477D0, &
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
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
                  XMIN=0.477D0,XMAX=4.477D0,YMIN=-0.523D0,YMAX=3.477D0, &
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
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
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
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
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
    INTEGER:: ndate,nplot,ncountry
    REAL(dp):: rcases,rdeaths

    SELECT CASE(id)
    CASE(1,3)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(2,ndate_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          DO ndate=1,ndate_max
             rcases=DBLE(ncases_total_ndate_ncountry(ndate,ncountry))
             rdeaths=DBLE(ndeaths_total_ndate_ncountry(ndate,ncountry))
             fg(1,ndate,nplot)=LOG10(MAX(cases_number_log_min,rcases))
             fg(2,ndate,nplot)=LOG10(MAX(deaths_number_log_min,rdeaths))
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
       ALLOCATE(fg(2,ndate_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          DO ndate=1,ndate_max
             rcases=DBLE(ncases_total_ndate_ncountry(ndate,ncountry)) &
                  /population_ncountry(ncountry)
             rdeaths=DBLE(ndeaths_total_ndate_ncountry(ndate,ncountry)) &
                  /population_ncountry(ncountry)
             fg(1,ndate,nplot)=LOG10(MAX(cases_rate_log_min,rcases))
             fg(2,ndate,nplot)=LOG10(MAX(deaths_rate_log_min,rdeaths))
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

    REAL(dp):: xg(ndate_max)
    REAL(dp),ALLOCATABLE:: fg(:,:,:)
    INTEGER,ALLOCATABLE:: ndate_maxa(:)
    INTEGER:: ndata,ndate,ngid,mode_2d
    CHARACTER(LEN=80):: title

    IF(nplot_max.EQ.1) THEN
       mode_2d=22
    ELSE
       mode_2d=21
    END IF
    
    ALLOCATE(ndate_maxa(nplot_max))
    ndate_maxa(1:nplot_max)=ndate_max

    CALL PAGES
    SELECT CASE(id)
    CASE(4,40,41)
       ngid=1
       IF(id.EQ.41) ngid=0
       CALL cv_gsub41(1,ncountry_nplot,nplot_max,fg,title,ndata)
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
                  XMIN=0.D0,XMAX=5.D0,YMIN=-1.D0,YMAX=4.D0, &
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
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
                  XMIN=-2.D0,XMAX=3.D0,YMIN=-3.D0,YMAX=2.D0, &
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
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
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
       CALL grdxy(ngid,fg,2,ndate_max,ndate_maxa,ndata,title, &
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
    INTEGER:: ndate,nplot,ncountry,ndate_ave
    REAL(dp):: rcases,rdeaths

    SELECT CASE(id)
    CASE(1,3)
       ndata=nplot_max
       IF(ALLOCATED(fg)) DEALLOCATE(fg)
       ALLOCATE(fg(2,ndate_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          DO ndate=1,ndays_ave-1
             fg(1,ndate,nplot)=LOG10(ratio_new_total_log_min &
                                    *cases_number_log_min)
             fg(2,ndate,nplot)=LOG10(ratio_new_total_log_min &
                                    *deaths_number_log_min)
          END DO
          DO ndate=ndays_ave,ndate_max
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
             fg(1,ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                        *cases_number_log_min,rcases))
             fg(2,ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
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
       ALLOCATE(fg(2,ndate_max,nplot_max))
       DO nplot=1,nplot_max
          ncountry=ncountry_nplot(nplot)
          DO ndate=1,ndays_ave-1
             fg(1,ndate,nplot)=LOG10(ratio_new_total_log_min &
                                    *cases_rate_log_min)
             fg(2,ndate,nplot)=LOG10(ratio_new_total_log_min &
                                    *deaths_rate_log_min)
          END DO
          DO ndate=ndays_ave,ndate_max
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
             fg(1,ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
                                        *cases_rate_log_min,rcases))
             fg(2,ndate,nplot)=LOG10(MAX(ratio_new_total_log_min &
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
