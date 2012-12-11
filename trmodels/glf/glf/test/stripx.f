c@stripx   Glenn Bateman
c  rgb 06-oct-94 changed line length from 80 to 132 characters
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
      integer, parameter :: kc = 132
      integer            :: nin, ntrans, nout
      integer            :: ilength, ixlen, j
c
      character line * 132
c
  10  read (nin,100,end=20) line
 100  format (a132)
c
c..find number of characters before spaces
c
      ilength = 0
      do j=1,kc
        if ( line(j:j) .ne. ' ' ) ilength = j
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
