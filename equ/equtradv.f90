! equtradv.f90

      module tradv_mod
	  use eqtrn_mod
      public
      contains
!
!=======================================================================
      subroutine tradv
!=======================================================================
      use aaa_mod
      use trn_mod
      use cnt_mod
      implicit none
! local variables
      integer    is,n
!=======================================================================
      do n=1,nro
      do is=0,mion
      tem(n,is)=1.05d+00*tem(n,is)
      pre(n,is)=1.05d+00*pre(n,is)
      enddo
      enddo
                  write(ft06,*)'==pre has been increased by 5%'
				  
!	  call equtrsvq

      time=time+1.
      call eqtrn
!-----------------------------------------------------------------------
      return
      end subroutine tradv
!
      end module tradv_mod
