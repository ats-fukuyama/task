  ! sak4.f90

MODULE sak4
  PRIVATE
  PUBLIC sak_4

CONTAINS

  SUBROUTINE sak_4

    USE sakcomm
    USE saksub
    USE libgrf
    IMPLICIT NONE

    REAL(dp):: wr,wi,rk,sg
    REAL(dp):: rkmin,rkmax,sgmin,sgmax
    REAL(dp):: delrk,delsg
    REAL(dp):: delta_nw,eps_nw,rd
    INTEGER:: lmax_nw,list_nw,ierr
    INTEGER:: mode,nrkmax,nsgmax,nrk,nsg
    REAL(dp),ALLOCATABLE:: xa(:),ya(:),fa(:,:)
    REAL(dp),ALLOCATABLE:: rka(:),sga(:)
    REAL(dp),ALLOCATABLE:: wr1a(:,:),wi1a(:,:),wi2a(:,:),wra(:,:),wia(:,:)
    CHARACTER(LEN=80):: title

    mode=1  ! 1: rk 1dplot, 2: sg 1dplot, 3: rk-sg 2Dplot
    rkmin=0.D0
    rkmax=0.55D0
    sgmin=0.D0
    sgmax=1.D0
    nrkmax=101
    nsgmax=1

    delta_nw = 1.D-6
    eps_nw = 1.D-8
    lmax_nw= 40
    list_nw= 0

1   CONTINUE
    WRITE(6,'(A)') '## INPUT: mode,nrkmax,rkmin,rkmax,nsgmax,sgmin,sgmax ?'
    WRITE(6,'(A,2I6,2ES12.4,I6,2ES12.4)') &
         '## ',mode,nrkmax,rkmin,rkmax,nsgmax,sgmin,sgmax
    READ(5,*,ERR=1,END=9000) mode,nrkmax,rkmin,rkmax,nsgmax,sgmin,sgmax
    IF(mode.EQ.0) GO TO 9000

    ALLOCATE(rka(nrkmax),sga(nsgmax))
    ALLOCATE(wr1a(nrkmax,nsgmax),wi1a(nrkmax,nsgmax),wi2a(nrkmax,nsgmax))
    ALLOCATE(wra(nrkmax,nsgmax),wia(nrkmax,nsgmax))

    IF(nrkmax.EQ.1) THEN
       delrk=0.D0
    ELSE
       delrk=(rkmax-rkmin)/(nrkmax-1)
    END IF
    IF(nsgmax.EQ.1) THEN
       delsg=0.D0
    ELSE
       delsg=(sgmax-sgmin)/(nsgmax-1)
    END IF
    DO nrk=1,nrkmax
       rka(nrk)=rkmin+delrk*(nrk-1)
    END DO
    DO nsg=1,nsgmax
       sga(nsg)=sgmin+delsg*(nsg-1)
    END DO

    DO nsg=1,nsgmax
       DO nrk=1,nrkmax
          CALL cwaprx(rka(nrk),sga(nsg), &
               wr1a(nrk,nsg),wi1a(nrk,nsg),wi2a(nrk,nsg))
          CALL set_rksg(rka(nrk),sga(nsg))
          CALL newtn0(subeps1,wr1a(nrk,nsg),wi2a(nrk,nsg), &
               wra(nrk,nsg),wia(nrk,nsg),rd, &
               delta_nw,eps_nw,lmax_nw,list_nw,ierr)
       END DO
    END DO

    CALL PAGES
    SELECT CASE(mode)
    CASE(1)
       ALLOCATE(xa(nrkmax),fa(nrkmax,3))
       xa(1:nrkmax)=rka(1:nrkmax)
       fa(1:nrkmax,1)=wra(1:nrkmax,1)
       fa(1:nrkmax,2)=wr1a(1:nrkmax,1)
       title='@wr vs rk: exact,approx@'
       CALL grd1d(1,xa,fa,nrkmax,nrkmax,2,title)
       fa(1:nrkmax,1)=wia(1:nrkmax,1)
       fa(1:nrkmax,2)=wi2a(1:nrkmax,1)
       fa(1:nrkmax,3)=wi1a(1:nrkmax,1)
       title='@wi vs rk: exact,approx1,approx2@'
       CALL grd1d(2,xa,fa,nrkmax,nrkmax,3,title)
       DEALLOCATE(xa,fa)
    CASE(2)
       ALLOCATE(xa(nsgmax),fa(nsgmax,3))
       xa(1:nsgmax)=sga(1:nsgmax)
       fa(1:nsgmax,1)=wra(1,1:nsgmax)
       fa(1:nsgmax,2)=wr1a(1,1:nsgmax)
       title='@wr vs sg: exact,approx@'
       CALL grd1d(1,xa,fa,nsgmax,nsgmax,2,title)
       fa(1:nsgmax,1)=wia(1,1:nsgmax)
       fa(1:nsgmax,2)=wi2a(1,1:nsgmax)
       fa(1:nsgmax,3)=wi1a(1,1:nsgmax)
       title='@wi vs sg: exact,approx1,approx2@'
       CALL grd1d(2,xa,fa,nsgmax,nsgmax,3,title)
       DEALLOCATE(xa,fa)
    CASE(3)
       ALLOCATE(xa(nrkmax),fa(nrkmax,nsgmax))
       xa(1:nrkmax)=rka(1:nrkmax)
       DO nsg=1,nsgmax
          fa(1:nrkmax,nsg)=wra(1:nrkmax,nsg)
       END DO
       title='@wr vs rk: for various sg@'
       CALL grd1d(1,xa,fa,nrkmax,nrkmax,nsgmax,title)
       DO nsg=1,nsgmax
          fa(1:nrkmax,nsg)=wia(1:nrkmax,nsg)
       END DO
       title='@wi vs rk: for various sg@'
       CALL grd1d(2,xa,fa,nrkmax,nrkmax,nsgmax,title)
       DEALLOCATE(xa,fa)
    CASE(4)
       ALLOCATE(xa(nsgmax),fa(nsgmax,nrkmax))
       xa(1:nsgmax)=sga(1:nsgmax)
       DO nrk=1,nrkmax
          fa(1:nsgmax,nrk)=wra(nrk,1:nsgmax)
       END DO
       title='@wr vs sg: for various rk@'
       CALL grd1d(1,xa,fa,nsgmax,nsgmax,nrkmax,title)
       DO nrk=1,nrkmax
          fa(1:nsgmax,nrk)=wia(nrk,1:nsgmax)
       END DO
       title='@wi vs sg: for various rk@'
       CALL grd1d(2,xa,fa,nsgmax,nsgmax,nrkmax,title)
       DEALLOCATE(xa,fa)
    CASE(5)
       ALLOCATE(xa(nrkmax),ya(nsgmax),fa(nrkmax,nsgmax))
       xa(1:nrkmax)=rka(1:nrkmax)
       ya(1:nsgmax)=sga(1:nsgmax)
       DO nsg=1,nsgmax
          fa(1:nrkmax,nsg)=wra(1:nrkmax,nsg)
       END DO
       title='@wr(rk,sg)@'
       CALL grd2d(1,xa,ya,fa,nrkmax,nrkmax,nsgmax,title, &
                  ASPECT=1.D0)
       DO nsg=1,nsgmax
          fa(1:nrkmax,nsg)=wia(1:nrkmax,nsg)
       END DO
       title='@wi(rk,sg)@'
       CALL grd2d(2,xa,ya,fa,nrkmax,nrkmax,nsgmax,title, &
                  ASPECT=1.D0)
       DEALLOCATE(xa,ya,fa)
    END SELECT
    CALL PAGEE

    DEALLOCATE(rka,sga)
    DEALLOCATE(wr1a,wi1a,wi2a)
    DEALLOCATE(wra,wia)

    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE sak_4
END MODULE sak4

       
       
    
    
    
