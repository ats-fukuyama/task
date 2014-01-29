!C
!C Input parapmeter interface
!C
MODULE T2PARM
  
  USE T2CNST, ONLY: i0ikind,i0rkind
  
  PRIVATE
  PUBLIC T2_PARM,T2_VIEW
  
CONTAINS

  !C
  !C 
  !C

  SUBROUTINE T2_PARM(MODE,KIN,IERR)
!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR

1   CALL TASK_PARM(MODE,'T2',KIN,T2_NLIN,T2_PLST,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl T2_CHECK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE T2_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE T2_NLIN(NID,IST,IERR)

    USE T2COMM, ONLY: &
         c10rname, i0dbg, i0fnum, i0mfcs, i0wstp,&
         i0dmax,i0amax,&
!         i0tmax, d0tstp, d0tmax,&
         i0smax, i0nmax, i0lmax, i1mlvl,&
         i0pdiv_number, i1rdn2, d1rec,&
         i0pmax,d0eps,d0rmjr,d0rmnr,&
         i0m0,i0n0,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,&
         d1pa,d1pz,&
         d0qc,d0qs,d0rw,&
         dt,time_init,eps_conv, &
         ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,nconvmax, &
         idfile,idprint,idplot,idmode,idebug

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR

    NAMELIST /T2/ &
         c10rname, i0dbg, i0fnum, i0mfcs, i0wstp,&
         i0dmax,i0amax,&
!         i0tmax, d0tstp, d0tmax,&
         i0smax, i0nmax, i0lmax, i1mlvl,&
         i0pdiv_number, i1rdn2, d1rec,&
!         i0pmax,d0eps, &
         d0rmjr,d0rmnr,&
         i0m0,i0n0,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,&
         d1pa,d1pz,&
         d0qc,d0qs,d0rw, &
         dt,time_init,eps_conv, &
         ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,nconvmax, &
         idfile,idprint,idplot,idmode,idebug


    READ(NID,T2,IOSTAT=IST,ERR=9800,END=9900)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE T2_NLIN

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE T2_PLST

    IMPLICIT NONE
    WRITE(6,'(A)') '# &T2 : c10rname,i0dbg,i0fnum,i0mfcs,i0wstp,'
    WRITE(6,'(A)') '        i0dmax,i0amax,' !i0tmax,d0tstp,d0tmax,'
    WRITE(6,'(A)') '        i0smax,i0nmax,i0lmax,'
    WRITE(6,'(A)') '        i1mlvl,i0pdiv_number, i1rdn2, d1rec,'
!    WRITE(6,'(A)') '        i0pmax,d0eps,'
    WRITE(6,'(A)') '        d0rmjr,d0rmnr,'
    WRITE(6,'(A)') '        i0m0,i0n0,d0bc,'
    WRITE(6,'(A)') '        d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,'
    WRITE(6,'(A)') '        d1pa,d1pz,'
    WRITE(6,'(A)') '        d0qc,d0qs,d0rw,'
    WRITE(6,'(A)') '        dt,time_init,eps_conv,'
    WRITE(6,'(A)') '        ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,'
    WRITE(6,'(A)') '        nconvmax,idfile,idprint,idplot,idmode,idebug'
    RETURN

  END SUBROUTINE T2_PLST

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE T2_CHECK(IERR)

    USE T2COMM,ONLY: i0pdiv_number
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(i0pdiv_number <= 0) THEN
       WRITE(6,'(A,I12)') &
            'XX T2_check: INVALID i0pdiv_number: ',i0pdiv_number
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE T2_CHECK

!     ****** SHOW PARAMETERS ******

  SUBROUTINE T2_VIEW

    USE T2COMM, ONLY: &
         c10rname, &
         i0dbg, i0fnum, i0mfcs, i0wstp,&
         i0dmax,i0amax,&
!         i0tmax, d0tstp, d0tmax,&
         i0smax, i0nmax, i0lmax, i1mlvl,&
         i0pdiv_number, i1rdn2, d1rec,&
!         i0pmax,d0eps,&
         d0rmjr,d0rmnr,&
         i0m0,i0n0,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,&
         d1pa,d1pz,&
         d0qc,d0qs,d0rw,&
         dt,time_init,eps_conv, &
         ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,nconvmax, &
         idfile,idprint,idplot,idmode,idebug

    IMPLICIT NONE
    INTEGER(i0ikind):: i1

    WRITE(6,601) &
         'i0dbg        ',i0dbg , &
         'i0fnum       ',i0fnum, &
         'i0mfcs       ',i0mfcs
    WRITE(6,601) &
         'i0wstp       ',i0wstp, &
         'i0dmax       ',i0dmax, &
         'i0amax       ',i0amax
    WRITE(6,602) &
         'ntmax        ',ntmax, &
         'ntstep       ',ntstep, &
         'dt           ',dt
    WRITE(6,601) &
         'i0smax       ',i0smax, &
         'i0nmax       ',i0nmax, &
         'i0lmax       ',i0lmax
    WRITE(6,601) &
         'i0pdiv_number',i0pdiv_number

    WRITE(6,'(A)') '      i1  i1mlvl  i1rdns       d1rec'
    DO i1=1,i0lmax
       WRITE(6,'(3I8,1PE12.4)') i1,i1mlvl(i1),i1rdn2(i1),d1rec(i1)
    END DO

    WRITE(6,603) &
         'nconvmax     ',nconvmax,&
         'eps_conv     ',eps_conv

    WRITE(6,601) &     
         'idfile       ',idfile, &
         'idprint      ',idprint, &
         'idplot       ',idplot
    WRITE(6,601) &     
         'idmode       ',idmode, &
         'idebug       ',idebug


    WRITE(6,604) &
         'd0rmjr       ',d0rmjr,&
         'd0rmnr       ',d0rmnr
    WRITE(6,604) &
         'd0qc         ',d0qc, &
         'd0qs         ',d0qs,&
         'd0rw         ',d0rw
    WRITE(6,602) &
         'i0m0         ',i0m0, &
         'i0n0         ',i0n0,&
         'd0bc         ',d0bc
    WRITE(6,'(A)')
    DO i1=1,i0smax
       WRITE(6,'(I6,1P2E12.4/6X,1P6E12.4)') &
            i1,d1pa(i1),d1pz(i1),&
            d1nc(i1),d1ns(i1),d1nw(i1),d1tc(i1),d1ts(i1),d1tw(1)
    END DO
    RETURN

601 FORMAT(' ',A11,'=',I8,4X  :2X,A11,'=',I8,4X  : &
            2X,A11,'=',I8)
602 FORMAT(' ',A11,'=',I8,4X  :2X,A11,'=',I8,4X  : &
            2X,A11,'=',1PE12.4)
603 FORMAT(' ',A11,'=',I8,4X  :2X,A11,'=',1PE12.4: &
            2X,A11,'=',1PE12.4)
604 FORMAT(' ',A11,'=',1PE12.4:2X,A11,'=',1PE12.4: &
            2X,A11,'=',1PE12.4)
  END SUBROUTINE T2_VIEW

END MODULE T2PARM
