      module tradv_mod
	  use eqtrn_mod
      public
      contains
c
c=======================================================================
      subroutine tradv
c=======================================================================
      use aaa_mod
      use trn_mod
      use cnt_mod
      implicit none
! local variables
      integer    is,n
c=======================================================================
      do n=1,nro
      do is=0,mion
      tem(n,is)=1.05d+00*tem(n,is)
      pre(n,is)=1.05d+00*pre(n,is)
      enddo
      enddo
                  write(ft06,*)'==pre has been increased by 5%'
				  
!	  call trsvq

      time=time+1.
      call eqtrn
c-----------------------------------------------------------------------
      return
      end subroutine tradv
c
      end module tradv_mod
