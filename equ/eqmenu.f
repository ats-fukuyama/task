C     $Id$
C
      module eqmenu_mod
      use eqinit_mod
      use eqgout_mod
      use equunit_mod
      use trunit_mod
      use equ_eqdsk
      public
      contains
c
c   ***** topics/equ menu *****
c
      subroutine eqmenu
c
      implicit none
      integer init,mstat,ierr,mode
      character knam*80,kpname*80
      character kid*1,line*80
      save init,mstat,kpname
      data init/0/
c
      if(init.eq.0) then
         mstat=0
         kpname='equparm'
         init=1
      endif
c
    1 continue
         ierr=0
         write(6,601) 
  601    format('## eq menu: R/RUN  C/CONT  P,V,F/PARM  G/GRAPH',
     &                  '  S,L,E,K/FILE  Q/QUIT')
c
         call task_klin(line,kid,mode,eqparm)
      if(mode.ne.1) goto 1
c
      if(kid.eq.'R') then
         mstat=0
         call equ_prof
         call tr_prof
         call equ_calc
         mstat=1
c
      elseif(kid.eq.'C') then
         if(mstat.ge.1) then 
            call tr_exec
            call equ_calc
         else
            write(6,*) 'XX: No data for continuing calculation!'
         endif
c
      elseif(kid.eq.'P') then
         call equ_parm(0,'equ',ierr)
c
      elseif(kid.eq.'V') then
         call equ_view
c
      elseif(kid.eq.'F') then
         call equ_parm(1,'equparm',ierr)
c
      elseif(kid.eq.'G') then
         call equ_gout
c
      elseif(kid.eq.'S') then
         call equ_save
      elseif(kid.eq.'L') then
c         call ktrim(knameq,kl)
c   10    write(6,*) '#equ> input : eqdata file name : ',knameq(1:kl)
c         read(5,'(a80)',err=10,end=9000) knam
c         if(knam(1:2).ne.'/ ') knameq=knam
c
c         modelg=3
c         call eqread(0,ierr)
c         if(ierr.ne.0) goto 10
c         mstat=2
c
      elseif(kid.eq.'K') then
c         call ktrim(knameq,kl)
c   11    write(6,*) '#equ> input : eqdata file name : ',knameq(1:kl)
c         read(5,'(A80)',ERR=11,END=9000) knam
c         if(knam(1:2).ne.'/ ') knameq=knam
c
c         modelg=5
c         call eqread(0,ierr)
c         if(ierr.ne.0) goto 11
c         mstat=2
c
      elseif(kid.eq.'E') then
         CALL equ_eqdsk_write(ierr)
      else if(kid.eq.'X'.or.kid.eq.'#') then
         continue
      elseif(kid.eq.'Q') then
         goto 9000
      else
         write(6,*) 'XX eqmenu: unknown kid'
      endif
      goto 1
c
 9000 return
      end subroutine eqmenu
c
      end module eqmenu_mod
