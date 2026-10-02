      module eqtrn_mod
      public
      contains
c
c=======================================================================
      subroutine eqtrn
c=======================================================================
c     pre-programed coil current
c=======================================================================
      use aaa_mod
      use vac_mod
      use cnt_mod
c-----
      implicit none
      integer i,kkstep
      real*8 cvac0(icvdm)
      data kkstep/0/

      save kkstep
C      save cvac0,cvact
      save cvac0

C      namelist/eqtr/cvact
c-----------------------------------------------------------------------
      write(6,*)'  entered eqtrn '

      if(kkstep.eq.0)then
C        open(8,FILE='trparm',STATUS='old',FORM='FORMATTED')
C        rewind 8
C        read(8,eqtr,err=999,end=998)
C        close(8)
C        write(6,eqtr)
        do i=1,icvdm
          cvac0(i)=cvac(i)
        enddo
        kkstep=1
C 999   if(kkstep.eq.0)then
C         kkstep=-1
C        endif
C 998    if(kkstep.eq.0)then
C          kkstep=-1
C        endif        
      endif
c-----
C      if(kkstep.lt.0)then
C        return
C      endif
c-----
      do i=1,icvdm
        if(ivac(i).lt.0)then
          cvac(i)=cvac0(i)+cvact(i)*time
        endif
      enddo
c-----------------------------------------------------------------------
      return
      end subroutine eqtrn
     
      end module eqtrn_mod 
