! equtrn.f90

      module eqtrn_mod
      public
      contains
!
!=======================================================================
      subroutine eqtrn
!=======================================================================
!     pre-programed coil current
!=======================================================================
      use aaa_mod
      use vac_mod
      use cnt_mod
!-----
      implicit none
      integer i,kkstep
      real*8 cvac0(icvdm)
      data kkstep/0/

      save kkstep
!      save cvac0,cvact
      save cvac0

!      namelist/eqtr/cvact
!-----------------------------------------------------------------------
      write(6,*)'  entered eqtrn '

      if(kkstep.eq.0)then
!        open(8,FILE='trparm',STATUS='old',FORM='FORMATTED')
!        rewind 8
!        read(8,eqtr,err=999,end=998)
!        close(8)
!        write(6,eqtr)
        do i=1,icvdm
          cvac0(i)=cvac(i)
        enddo
        kkstep=1
! 999   if(kkstep.eq.0)then
!         kkstep=-1
!        endif
! 998    if(kkstep.eq.0)then
!          kkstep=-1
!        endif        
      endif
!-----
!      if(kkstep.lt.0)then
!        return
!      endif
!-----
      do i=1,icvdm
        if(ivac(i).lt.0)then
          cvac(i)=cvac0(i)+cvact(i)*time
        endif
      enddo
!-----------------------------------------------------------------------
      return
      end subroutine eqtrn
     
      end module eqtrn_mod 
