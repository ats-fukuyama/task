c@xtverb  .../bateman/appl/xtverb.f
c
c    This program extracts lines of text between each
c  occurance of \begin{verbatim} and \end{verbatim}.
c
c    To compile and link this program, type:
c
c  f77 -o xtverb xtverb.f 
c
c    To run the resulting controllee, type:
c
c  xtverb < inputfile > outputfile
c
c  where
c    inputfile  is the name of the LaTeX file with imbedded program
c    outputfile is the desired name for the output file.
c
c  Note:  The original version of this program searched for the
c    character strings '/begin...' and '/end...'.  However, on the
c    Sun workstations, the / in '/...' is treated as an escape
c    character.  In order to make this program work the same way
c    on all computers, I had to remove the '/' character.
c
c
      character line*80
      logical lverb
      integer lastcol, lenght, j
c
c
c  lverb = .true. only when the text is located between
c          \begin{verbatim} and \end{verbatim}
c  length = number of nonblank characters in each line (up to 80)
c  lastcol = the last column were xtverb will look for begin{verbatim}
c            or end{verbatim}
c
c
      lastcol = 20
c
      lverb = .false.
c
  10  continue
c
      read (5,100,end=90) line
 100  format (a)
c
        length = 1
      do 20 j=1,80
        if ( line(j:j) .ne. ' ' ) length = j
 20   continue
c
      do 22 j=2,lastcol
        if ( line(j:j+12) .eq. 'end{verbatim}' ) lverb = .false.
 22   continue
c
      if ( lverb ) write (6,100) line(1:length)
c
      do 24 j=2,lastcol
        if ( line(j:j+14) .eq. 'begin{verbatim}' ) lverb = .true.
 24   continue
c
      go to 10
c
  90  continue
c
      call exit(0)
      stop
      end

