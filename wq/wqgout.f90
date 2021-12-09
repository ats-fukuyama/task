! wqgout.f90

MODULE wqgout

  PRIVATE
  PUBLIC wq_gout

CONTAINS

  ! *** graphic menu ***
  
  SUBROUTINE wq_gout
    USE wqcomm
    USE wqparm,ONLY: wq_parm,wq_broadcast
    USE wqview,ONLY: wq_view
    USE libchar
    USE libkio
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80):: line
    CHARACTER(LEN=1):: kid
    INTEGER:: mode

1   CONTINUE
    
    IF(nrank.EQ.0) THEN
       WRITE(6,'(A)') &
            '## wqgout menu: g:global p:profile t:time  parm  x/end'
       CALL TASK_KLIN(line,kid,mode,wq_parm)
       CALL ToUpper(kid)
    END IF
    CALL mtx_broadcast1_character(kid)
    CALL mtx_broadcast1_integer(mode)

    IF(mode.EQ.2) CALL wq_broadcast
    IF(mode.NE.1) GO TO 1

    SELECT CASE(kid)
    CASE('G')
       CALL wq_gsub_global
    CASE('P')
       CALL wq_gsub_profile
    CASE('T')
       CALL wq_gsub_time
    CASE('X')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wq_gout

  ! *** time history of global quantities ***
  
  SUBROUTINE wq_gsub_global

    USE wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: ngt
    REAL(rkind):: tn(ngt_m)

    DO ngt=1,ngt_max
       tn(ngt)=t_ngt(ngt)/period
    END DO

    CALL PAGES
    WRITE(6,'(A,I8,2ES12.4)') 'ngt_max,t,EX:',ngt_max,t_ngt(1),EX_max_ngt(1)
    CALL grd1d(1,tn,EX_max_ngt,ngt_m,ngt_max,1,'@Ex_max vs tn@')
    CALL grd1d(2,tn,EY_max_ngt,ngt_m,ngt_max,1,'@EY_max vs tn@')
    CALL grd1d(3,tn,EZ_max_ngt,ngt_m,ngt_max,1,'@EZ_max vs tn@')
    CALL grd1d(4,tn,pabs_tot_ngt,ngt_m,ngt_max,1,'@Pabs_tot vs tn@')

    CALL PAGEE
    RETURN
  END SUBROUTINE wq_gsub_global

  ! *** profile of elecric field and pabs ***
  
  SUBROUTINE wq_gsub_profile
  
    USE wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: nx,ny
    REAL(rkind):: fr(nxmax,nymax),fi(nxmax,nymax),fa(nxmax,nymax)

    CALL PAGES

    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(EX(nx,ny))
          fi(nx,ny)=AIMAG(EX(nx,ny))
          fa(nx,ny)=ABS(EX(nx,ny))
       END DO
    END DO
    CALL grd2d(5,xg,yg,fr,nxmax,nxmax,nymax,'@Real Ex@',ASPECT=0.D0)
    CALL grd2d(6,xg,yg,fi,nxmax,nxmax,nymax,'@Imag Ex@',ASPECT=0.D0)
    CALL grd2d(7,xg,yg,fa,nxmax,nxmax,nymax,'@Abs  Ex@',ASPECT=0.D0)
    
    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(EY(nx,ny))
          fi(nx,ny)=AIMAG(EY(nx,ny))
          fa(nx,ny)=ABS(EY(nx,ny))
       END DO
    END DO
    CALL grd2d( 8,xg,yg,fr,nxmax,nxmax,nymax,'@Real Ey@',ASPECT=0.D0)
    CALL grd2d( 9,xg,yg,fi,nxmax,nxmax,nymax,'@Imag Ey@',ASPECT=0.D0)
    CALL grd2d(10,xg,yg,fa,nxmax,nxmax,nymax,'@Abs  Ey@',ASPECT=0.D0)
    
    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(EZ(nx,ny))
          fi(nx,ny)=AIMAG(EZ(nx,ny))
          fa(nx,ny)=ABS(EZ(nx,ny))
       END DO
    END DO
    CALL grd2d(11,xg,yg,fr,nxmax,nxmax,nymax,'@Real Ez@',ASPECT=0.D0)
    CALL grd2d(12,xg,yg,fi,nxmax,nxmax,nymax,'@Imag Ez@',ASPECT=0.D0)
    CALL grd2d(13,xg,yg,fa,nxmax,nxmax,nymax,'@Abs  Ez@',ASPECT=0.D0)

    CALL PAGEE

    IF(ngt_max.EQ.0) RETURN
    IF(pabs_tot_ngt(ngt_max).LE.0.D0) RETURN

    CALL PAGES
    
    DO ny=1,nymax
       DO nx=1,nxmax
          fr(nx,ny)=REAL(pabs(nx,ny))
       END DO
    END DO
    CALL grd2d(5,xg,yg,fr,nxmax,nxmax,nymax,'@Pabs@',ASPECT=0.D0)

    CALL PAGEE
    RETURN
  END SUBROUTINE wq_gsub_profile

  ! *** time history of elecric field and pabs ***
  
  SUBROUTINE wq_gsub_time
  
    USE wqcomm
    USE libgrf
    IMPLICIT NONE
    INTEGER:: id,ngr,nx,ny,npage,npage_max,ng_id,ng
    CHARACTER(LEN=9):: title
    CHARACTER(LEN=16):: title_t
    REAL(rkind):: f(nxmax,nymax)

1   CONTINUE
    WRITE(6,'(A)') '## Input : 1:Exr 2:Exi 3:Exa 4:Eyr 5:Eyi 6:Eya'
    WRITE(6,'(A)') '           7:Ezr 8:Ezi 9:Eza 10:Pabs     0:end'
    READ(5,*,END=9000,ERR=1) id
    IF(id.EQ.0) GO TO 9000

    SELECT CASE(id)
    CASE(1)
       title='@Real Ex:'
    CASE(2)
       title='@Imag Ex:'
    CASE(3)
       title='@Abs  Ex:'
    CASE(4)
       title='@Real Ey:'
    CASE(5)
       title='@Imag Ey:'
    CASE(6)
       title='@Abs  Ey:'
    CASE(7)
       title='@Real Ez:'
    CASE(8)
       title='@Imag Ez:'
    CASE(9)
       title='@Abs  Ez:'
    CASE(10)
       title='@Pabs   :'
    END SELECT

    npage_max=(ngr_max-1)/16+1
    IF(npage_max.EQ.1) THEN
       SELECT CASE(ngr_max)
       CASE(1)
          ng_id=0
       CASE(2:4)
          ng_id=1
       CASE(5:9)
          ng_id=5
       CASE(10:16)
          ng_id=14
       END SELECT
    ELSE
       ng_id=14
    END IF

    DO npage=1,npage_max
       CALL PAGES
       DO ng=1,MIN(ngr_max-16*(npage-1),16)
          ngr=16*(npage-1)+ng
          WRITE(title_t,'(A,ES12.4,A)') ' t=',t_ngr(ngr),'@'
          DO ny=1,nymax
             DO nx=1,nxmax
                SELECT CASE(id)
                CASE(1)
                   f(nx,ny)=REAL(EX_save(nx,ny,ngr))
                CASE(2)
                   f(nx,ny)=AIMAG(EX_save(nx,ny,ngr))
                CASE(3)
                   f(nx,ny)=ABS(EX_save(nx,ny,ngr))
                CASE(4)
                   f(nx,ny)=REAL(EY_save(nx,ny,ngr))
                CASE(5)
                   f(nx,ny)=AIMAG(EY_save(nx,ny,ngr))
                CASE(6)
                   f(nx,ny)=ABS(EY_save(nx,ny,ngr))
                CASE(7)
                   f(nx,ny)=REAL(EZ_save(nx,ny,ngr))
                CASE(8)
                   f(nx,ny)=AIMAG(EZ_save(nx,ny,ngr))
                CASE(9)
                   f(nx,ny)=ABS(EZ_save(nx,ny,ngr))
                CASE(10)
                   f(nx,ny)=pabs_save(nx,ny,ngr)
                END SELECT
             END DO ! nx
          END DO !ny

          CALL grd2d(ng-1+ng_id,xg,yg,f,nxmax,nxmax,nymax,title//title_t, &
                     ASPECT=0.D0)
       END DO
       CALL PAGEE
    END DO
    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wq_gsub_time
END MODULE wqgout
