!---- FEM quadratic ---

!---- FEM quadratic and linear ---

      subroutine fem_ppq(fmd1,fmd2,drho,fml)

      use libfem
      use wmfem_comm
      implicit none
      complex(8),dimension(3,3,4,nfcmax,nfcmax),intent(in):: fmd1,fmd2
      complex(8),dimension(mbmax,mbmax),intent(out):: fml
      real(8),intent(in):: drho
      integer:: ic1,ic2,ip1,ip2,mu1,mu2,nfc1,nfc2,mb1,mb2
      integer,dimension(8):: ica =(/2,3,1,2,3,1,2,3/) ! column number
      integer,dimension(8):: ipa =(/1,1,1,2,2,2,3,3/) ! position number

      if(table_initialize_flag.eq.0) then
         call table_initialize
         table_initialize_flag=1
      endif

      do mu2=1,8
         do mu1=1,8
            fml(mu1,mu2)=0.d0
         end do
      end do

      do mu1=1,8
         ic1=ica(mu1)
         ip1=ipa(mu1)
      do mu2=1,8
         ic2=ica(mu2)
         ip2=ipa(mu2)
         if((ic1 == 1) .AND. (ic2 == 1)) THEN
            do nfc1=1,nfcmax
               mb1=nfcmax*(mu1-1)+nfc1
            do nfc2=1,nfcmax
               mb2=nfcmax*(mu2-1)+nfc2
               fml(mb1,mb2)=fml(mb1,mb2)
     &        +fmd1(ic1,ic2,1,nfc1,nfc2)*table_ppp(1,ip1,  ip2  )*drho
     &        +fmd1(ic1,ic2,2,nfc1,nfc2)*table_ppp(1,ip1+2,ip2  )
     &        +fmd1(ic1,ic2,3,nfc1,nfc2)*table_ppp(1,ip1,  ip2+2)
     &        +fmd1(ic1,ic2,4,nfc1,nfc2)*table_ppp(1,ip1+2,ip2+2)/drho
     &        +fmd2(ic1,ic2,1,nfc1,nfc2)*table_ppp(2,ip1,  ip2  )*drho
     &        +fmd2(ic1,ic2,2,nfc1,nfc2)*table_ppp(2,ip1+2,ip2  )
     &        +fmd2(ic1,ic2,3,nfc1,nfc2)*table_ppp(2,ip1,  ip2+2)
     &        +fmd2(ic1,ic2,4,nfc1,nfc2)*table_ppp(2,ip1+2,ip2+2)/drho
            end do
            end do

         else if(ic1 == 1) THEN
            do nfc1=1,nfcmax
               mb1=nfcmax*(mu1-1)+nfc1
            do nfc2=1,nfcmax
               mb2=nfcmax*(mu2-1)+nfc2
               fml(mb1,mb2)=fml(mb1,mb2)
     &        +fmd1(ic1,ic2,1,nfc1,nfc2)*table_ppq(1,ip1,  ip2  )*drho
     &        +fmd1(ic1,ic2,2,nfc1,nfc2)*table_ppq(1,ip1+2,ip2  )
     &        +fmd1(ic1,ic2,3,nfc1,nfc2)*table_ppq(1,ip1,  ip2+3)
     &        +fmd1(ic1,ic2,4,nfc1,nfc2)*table_ppq(1,ip1+2,ip2+3)/drho
     &        +fmd2(ic1,ic2,1,nfc1,nfc2)*table_ppq(2,ip1,  ip2  )*drho
     &        +fmd2(ic1,ic2,2,nfc1,nfc2)*table_ppq(2,ip1+2,ip2  )
     &        +fmd2(ic1,ic2,3,nfc1,nfc2)*table_ppq(2,ip1,  ip2+3)
     &        +fmd2(ic1,ic2,4,nfc1,nfc2)*table_ppq(2,ip1+2,ip2+3)/drho
            end do
            end do

         else if(ic2 == 1) THEN
            do nfc1=1,nfcmax
               mb1=nfcmax*(mu1-1)+nfc1
            do nfc2=1,nfcmax
               mb2=nfcmax*(mu2-1)+nfc2
               fml(mb1,mb2)=fml(mb1,mb2)
     &        +fmd1(ic1,ic2,1,nfc1,nfc2)*table_pqp(1,ip1,  ip2  )*drho
     &        +fmd1(ic1,ic2,2,nfc1,nfc2)*table_pqp(1,ip1+3,ip2  )
     &        +fmd1(ic1,ic2,3,nfc1,nfc2)*table_pqp(1,ip1,  ip2+2)
     &        +fmd1(ic1,ic2,4,nfc1,nfc2)*table_pqp(1,ip1+3,ip2+2)/drho
     &        +fmd2(ic1,ic2,1,nfc1,nfc2)*table_pqp(2,ip1,  ip2  )*drho
     &        +fmd2(ic1,ic2,2,nfc1,nfc2)*table_pqp(2,ip1+3,ip2  )
     &        +fmd2(ic1,ic2,3,nfc1,nfc2)*table_pqp(2,ip1,  ip2+2)
     &        +fmd2(ic1,ic2,4,nfc1,nfc2)*table_pqp(2,ip1+3,ip2+2)/drho
            end do
            end do

         else
            do nfc1=1,nfcmax
               mb1=nfcmax*(mu1-1)+nfc1
            do nfc2=1,nfcmax
               mb2=nfcmax*(mu2-1)+nfc2
               fml(mb1,mb2)=fml(mb1,mb2)
     &        +fmd1(ic1,ic2,1,nfc1,nfc2)*table_pqq(1,ip1,  ip2  )*drho
     &        +fmd1(ic1,ic2,2,nfc1,nfc2)*table_pqq(1,ip1+3,ip2  )
     &        +fmd1(ic1,ic2,3,nfc1,nfc2)*table_pqq(1,ip1,  ip2+3)
     &        +fmd1(ic1,ic2,4,nfc1,nfc2)*table_pqq(1,ip1+3,ip2+3)/drho
     &        +fmd2(ic1,ic2,1,nfc1,nfc2)*table_pqq(2,ip1,  ip2  )*drho
     &        +fmd2(ic1,ic2,2,nfc1,nfc2)*table_pqq(2,ip1+3,ip2  )
     &        +fmd2(ic1,ic2,3,nfc1,nfc2)*table_pqq(2,ip1,  ip2+3)
     &        +fmd2(ic1,ic2,4,nfc1,nfc2)*table_pqq(2,ip1+3,ip2+3)/drho
            end do
            end do
         end if
      END DO
      END DO

      return
      end subroutine fem_ppq
