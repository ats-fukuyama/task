C-----------------------------------------------------------------------
C  C9DATE  --  return a 9 character date, e.g. " 5-Mar-93"
C    do it whether on VAX/VMS or on unix
C
C    mod dmc -- Y2K -- new output format:  ddmmmyyyy e.g. 05mar1993
C    using date_and_time routine, same code works vms and unix
C
      SUBROUTINE C9DATE(ZDATE)
C
      CHARACTER*(*) ZDATE ! THE DATE (OUTPUT)
C
      CHARACTER*10 UDATE
C
      character*30 foo1,foo2
      integer ivals(8) 
C
      integer imon
      character*3 mons(12)
      data mons/'jan','feb','mar','apr','may','jun','jul','aug',
     >          'sep','oct','nov','dec'/
C------------------------------------
C
      ZDATE=' '
C
      call date_and_time(udate,foo1,foo2,ivals)
      read(udate(5:6),'(i2)') imon
      zdate=udate(7:8)//mons(imon)//udate(1:4)
C
      RETURN
      END
