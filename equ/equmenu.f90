! equmanu.f90

!     $Id$
!
      module eqmenu_mod
      use eqinit_mod
      use eqgout_mod
      use equunit_mod
      use trunit_mod
      use equ_eqdsk
      USE equ_convert
      USE libkio
      public
      contains
!
!   ***** topics/equ menu *****
!
      subroutine eqmenu
!
      implicit none
      integer init,mstat,ierr,mode
      character knam*80,kpname*80
      character kid*1,line*80
      save init,mstat,kpname
      data init/0/
!
      if(init.eq.0) then
         mstat=0
         kpname='equparm'
         init=1
      endif
!
    1 continue
         ierr=0
         write(6,601) 
  601    format('## eq menu: R/RUN  C/CONT  P,V,F/PARM  G/GRAPH', &
                        '  S,L,E,K/FILE  Q/QUIT')
!
         call task_klin(line,kid,mode,eqparm)
      if(mode.ne.1) goto 1
!
      if(kid.eq.'R') then
         mstat=0
         call equ_prof
         call tr_prof
         call equ_calc
         mstat=1
!
      elseif(kid.eq.'C') then
         if(mstat.ge.1) then 
            call tr_exec
            call equ_calc
         else
            write(6,*) 'XX: No data for continuing calculation!'
         endif
!
      elseif(kid.eq.'P') then
         call equ_parm(0,'equ',ierr)
!
      elseif(kid.eq.'V') then
         call equ_view
!
      elseif(kid.eq.'F') then
         call equ_parm(1,'equparm',ierr)
!
      elseif(kid.eq.'G') then
         call equ_gout
!
      elseif(kid.eq.'S') then
         call equ_save
      elseif(kid.eq.'L') then
!         call ktrim(knameq,kl)
!   10    write(6,*) '#equ> input : eqdata file name : ',knameq(1:kl)
!         read(5,'(a80)',err=10,end=9000) knam
!         if(knam(1:2).ne.'/ ') knameq=knam
!
!         modelg=3
!         call eqread(0,ierr)
!         if(ierr.ne.0) goto 10
!         mstat=2
!
      elseif(kid.eq.'K') then
!         call ktrim(knameq,kl)
!   11    write(6,*) '#equ> input : eqdata file name : ',knameq(1:kl)
!         read(5,'(A80)',ERR=11,END=9000) knam
!         if(knam(1:2).ne.'/ ') knameq=knam
!
!         modelg=5
!         call eqread(0,ierr)
!         if(ierr.ne.0) goto 11
!         mstat=2
!
      elseif(kid.eq.'E') then
         CALL equ_eqdsk_write(ierr)
      elseif(kid.eq.'Z') then
         CALL equ_read_gfile_csv('coildata.csv')
      else if(kid.eq.'X'.or.kid.eq.'#') then
         continue
      elseif(kid.eq.'Q') then
         goto 9000
      else
         write(6,*) 'XX eqmenu: unknown kid'
      endif
      goto 1
!
 9000 return
      end subroutine eqmenu
!
      end module eqmenu_mod
