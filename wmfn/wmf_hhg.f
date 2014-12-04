!----- calculate coefficint matrix fmd (E cylindrical +,-,para)-----

      subroutine fem_hhg(fmd,drho,fml)

      use libfem
      use wmfem_comm
      implicit none
      complex(8),dimension(4,4,4,nfcmax,nfcmax,4),intent(in) :: fmd
!      complex(8),dimension(mbmax,mbmax):: fml
      complex(8),dimension(4*8*nfcmax-1,2*8*nfcmax),intent(out):: fml
      real(8),intent(in):: drho
      integer:: mr,mc,i,j,k,nf1,nf2,inod,ml,mw
      integer::ig

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif

      fml=0d0
      !!!!!!!!!!!!!!!!atteru?
      mr=8*nfcmax  ! line interval between radial points
      mc=16*nfcmax ! diagonal column at the top line
         do nf2=1, nfcmax
         do nf1=1, nfcmax
         do j=1,4
         do i=1,4
      !!!!!!!!!!!!!!!!atteru?
            ml=8*(nf1-1)+2*(i-1)+1
            mw=mc+8*(nf2-nf1)+2*(j-i)
            do inod=1,4
               ig=2*inod-1
               fml(mw     ,ml     )=fml(mw     ,ml     ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,1,ig)*drho
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,1,ig)
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,5,ig)
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,5,ig)/drho
               fml(mw+1   ,ml     )=fml(mw+1   ,ml     )
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,2,ig)*drho**2
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,2,ig)*drho
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,6,ig)*drho
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,6,ig)
               fml(mw+mr  ,ml     )=fml(mw+mr  ,ml     )
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,3,ig)*drho
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,3,ig)
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,7,ig)
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,7,ig)/drho
               fml(mw+mr+1,ml     )=fml(mw+mr+1,ml     )
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,4,ig)*drho**2
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,4,ig)*drho
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,8,ig)*drho
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,8,ig)
               fml(mw-1   ,ml+1   )=fml(mw-1   ,ml+1   )
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,1,ig)*drho**2
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,1,ig)*drho
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,5,ig)*drho
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,5,ig)
               fml(mw     ,ml+1   )=fml(mw     ,ml+1   ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,2,ig)*drho**3 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,2,ig)*drho**2 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,6,ig)*drho**2 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,6,ig)*drho
               fml(mw+mr-1,ml+1   )=fml(mw+mr-1,ml+1   ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,3,ig)*drho**2 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,3,ig)*drho 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,7,ig)*drho 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,7,ig)
               fml(mw+mr  ,ml+1   )=fml(mw+mr  ,ml+1   ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,4,ig)*drho**3 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,4,ig)*drho**2
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,8,ig)*drho**2 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,8,ig)*drho


               fml(mw-mr  ,ml+mr  )=fml(mw-mr  ,ml+mr  ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,1,ig)*drho 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,1,ig) 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,5,ig) 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,5,ig)/drho
               fml(mw-mr+1,ml+mr  )=fml(mw-mr+1,ml+mr  ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,2,ig)*drho**2 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,2,ig)*drho 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,6,ig)*drho 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,6,ig)
               fml(mw     ,ml+mr  )=fml(mw     ,ml+mr  ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,3,ig)*drho 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,3,ig) 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,7,ig) 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,7,ig)/drho
               fml(mw+1   ,ml+mr  )=fml(mw+1   ,ml+mr  ) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,4,ig)*drho**2 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,4,ig)*drho 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,8,ig)*drho 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,8,ig)

               fml(mw-mr-1,ml+mr+1)=fml(mw-mr-1,ml+mr+1) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,1,ig)*drho**2 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,1,ig)*drho 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,5,ig)*drho 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,5,ig)
               fml(mw-mr  ,ml+mr+1)=fml(mw-mr  ,ml+mr+1) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,2,ig)*drho**3 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,2,ig)*drho**2 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,6,ig)*drho**2 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,6,ig)*drho
               fml(mw-1   ,ml+mr+1)=fml(mw-1   ,ml+mr+1) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,3,ig)*drho**2 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,3,ig)*drho 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,7,ig)*drho 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,7,ig)
               fml(mw     ,ml+mr+1)=fml(mw     ,ml+mr+1) 
     &               +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,4,ig)*drho**3 
     &               +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,4,ig)*drho**2 
     &               +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,8,ig)*drho**2 
     &               +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,8,ig)*drho
            enddo
         enddo
         enddo
         enddo
         enddo

      return
      end subroutine fem_hhg
