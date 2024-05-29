c@stripx   Glenn Bateman
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine stripx (nin,ntrans,nout)
c
c  This sbrtn reads input data file logical unit number nin
c  one line at a time up to 80 characters long,
c  prints out the input verbatum to output file logical unit number nout
c  then searches for the first appearance of the character ! on each line
c  and prints out everything to the left of ! to output unit number ntrans
c
c  In this way, everything to the right of ! is treated as a comment
c
c  mod dmc (29 Oct 1998) -- to make the namelist work on a DEC machine, 
c  do the following transmogrifications:
c     a) namelist names &<name> --> $<name>
c     b) character / --> $end
c
c  this is because DEC and IBM, of course, have different syntax for the
c  start and end of a namelist.
c
c
      parameter (kc = 132)
c
      character line * 132
c
      character cputype * 10            ! start dmc 29-Oct 1998 new code
c
      cputype=' '
!      call getenv('CPU',cputype)
      if(cputype(1:3).eq.'DEC') then
         idec=1
      else
         idec=0
      endif                             ! end dmc 29-Oct 1998
c
  10  read (nin,100,end=20) line
 100  format (a132)
c
c..find number of characters before spaces
c
      ilength = 0
      do j=1,kc
         if ( line(j:j) .ne. ' ' ) then ! start dmc 29-Oct 1998 changes
            if(ilength.eq.0) then
c  first non-blank character...
               if(idec.eq.1) then
c  DEC machine -- transmogrify
                  if(line(j:j).eq.'&') line(j:j)='$' ! start of namelist
                  if(line(j:j).eq.'/') line(j:j+4)='$end!' ! end of n.l.
               endif
            endif
            ilength = j
         endif                          ! end dmc 29-Oct 1998
      enddo
c
c..echo input data on output file
c
      if ( ilength .gt. 0 ) then
        write (nout,110) line(1:ilength)
      else
        write (nout,*)
      endif
c
c..ixlen = number of nontrivial characters before "!"
c
      ixlen = 0
      do j=1,kc
        if ( line(j:j) .eq. '!' ) go to 14
        if ( line(j:j) .ne. ' ' ) ixlen = j
      enddo
  14  continue
c
c..echo input data up to "!" on temporary file
c
      if ( ixlen .gt. 0 ) write (ntrans,110) line(1:ixlen)
c
 110  format (a)
c
c
      go to 10
c
  20  continue
      rewind ntrans
c
      return
      end
