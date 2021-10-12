module fpwrite

contains

  subroutine fpcsv1D(f,filename)

    USE fpcomm,ONLY: rkind
     implicit none

     character(*),intent(in) :: filename
     REAL(rkind),  intent(in) :: f(:)
     integer :: imax,i

     imax=size(f)

     open(1985,file=filename)
       do i=1,imax
         if(i==imax)then
          if ( abs(f(i))<1.d-99 ) then
            write(1985,'(e18.11)')0.d0
          else if(f(i)>0.1d99) then
            write(1985,'(e18.11)')0.1d99
          else if(f(i)<-0.1d99) then
            write(1985,'(e18.11)')-0.1d99
          else
            write(1985,'(e18.11)')f(i)
          end if           
         else
          if ( abs(f(i))<1.d-99 ) then
            write(1985,'(e18.11,",")',advance="no")0.d0
          else if(f(i)>0.1d99) then
            write(1985,'(e18.11,",")',advance="no")0.1d99
          else if(f(i)<-0.1d99) then
            write(1985,'(e18.11,",")',advance="no")-0.1d99
          else
            write(1985,'(e18.11,",")',advance="no")f(i)
          end if
         end if
       end do
     close(1985)

  end subroutine fpcsv1D

  subroutine fpcsv2D(f,filename)
    USE fpcomm,ONLY: rkind
    implicit none
    character(*),intent(in) :: filename
    REAL(rkind),  intent(in) :: f(:,:)
    integer :: imax,i,jmax,j

    imax=size(f,1)
    jmax=size(f,2)

    open(1985,file=filename)
    do i=1,imax
     do j=1,jmax
       if(j==jmax)then
        if ( abs(f(i,j))<1.d-99 ) then
          write(1985,'(e18.11)')0.d0
        else if(f(i,j)>0.1d99) then
          write(1985,'(e18.11)')0.1d99
        else if(f(i,j)<-0.1d99) then
          write(1985,'(e18.11)')-0.1d99
        else
          write(1985,'(e18.11)')f(i,j)
        end if           
       else
        if ( abs(f(i,j))<1.d-99 ) then
          write(1985,'(e18.11,",")',advance="no")0.d0
        else if(f(i,j)>0.1d99) then
          write(1985,'(e18.11,",")',advance="no")0.1d99
        else if(f(i,j)<-0.1d99) then
          write(1985,'(e18.11,",")',advance="no")-0.1d99
        else
          write(1985,'(e18.11,",")',advance="no")f(i,j)
        end if           
       end if
     end do
    end do
    close(1985)

  end subroutine fpcsv2D

  subroutine fptxt1D(f,filename)
    USE fpcomm,ONLY: rkind
    implicit none
    REAL(rkind),intent(in) :: f(:)
    character(*),intent(in) :: filename
    integer :: i, imax

    imax = size(f)

    open(100,file=filename)
    write(100,'(I5)')imax
    do i = 1, imax-1
      write(100,'(e18.11,",")',advance='no')write_for_python(f(i))
    end do
      write(100,'(e18.11)')write_for_python(f(imax))
    close(100)

  end subroutine fptxt1D

  subroutine fptxt2D(f,filename)
    USE fpcomm,ONLY: rkind
    implicit none
    REAL(rkind),intent(in) :: f(:,:)
    character(*),intent(in) :: filename
    integer :: i1,i2,imax(2)

    imax(1) = size(f,1)
    imax(2) = size(f,2)

    open(100,file=filename)
    write(100,'(I5,",")',advance='no')imax(1)
    write(100,'(I5)')imax(2)
    do i2 = 1, imax(2)
      do i1 = 1, imax(1)-1
        write(100,'(e18.11,",")',advance='no')write_for_python(f(i1,i2))
      end do
      write(100,'(e18.11)')write_for_python(f(imax(1),i2))
    end do
    close(100)

  end subroutine fptxt2D

  subroutine fptxt3D(f,filename)
    USE fpcomm,ONLY: rkind
    implicit none
    REAL(rkind),intent(in) :: f(:,:,:)
    character(*),intent(in) :: filename
    integer :: i1, i2, i3, imax(3)

    imax(1) = size(f,1) 
    imax(2) = size(f,2) 
    imax(3) = size(f,3) 

    open(100,file=filename)
    write(100,'(I5,",")',advance='no')imax(1)
    write(100,'(I5,",")',advance='no')imax(2)
    write(100,'(I5)')imax(3)
    do i3 = 1, imax(3)
      do i2 = 1, imax(2)
        do i1 = 1, imax(1)-1
          write(100,'(e18.11,",")',advance='no')write_for_python(f(i1,i2,i3))
        end do
        write(100,'(e18.11)')write_for_python(f(imax(1),i2,i3))
      end do
    end do
    close(100)

  end subroutine fptxt3D

  subroutine fptxt4D(f,filename)
    USE fpcomm,ONLY: rkind
    implicit none
    REAL(rkind),intent(in) :: f(:,:,:,:)
    character(*),intent(in) :: filename
    integer :: i1, i2, i3, i4, imax(4)

    imax(1) = size(f,1) 
    imax(2) = size(f,2) 
    imax(3) = size(f,3) 
    imax(4) = size(f,4) 

    open(100,file=filename)
    write(100,'(I5,",")',advance='no')imax(1)
    write(100,'(I5,",")',advance='no')imax(2)
    write(100,'(I5,",")',advance='no')imax(3)
    write(100,'(I5)')imax(4)
    do i4 = 1, imax(4)
      do i3 = 1, imax(3)
        do i2 = 1, imax(2)
          do i1 = 1, imax(1)-1
            write(100,'(e18.11,",")',advance='no')write_for_python(f(i1,i2,i3,i4))
          end do
          write(100,'(e18.11)')write_for_python(f(imax(1),i2,i3,i4))
        end do
      end do
    end do
    close(100)

  end subroutine fptxt4D

  subroutine fptxt5D(f,filename)
    USE fpcomm,ONLY: rkind
    implicit none
    REAL(rkind),intent(in) :: f(:,:,:,:,:)
    character(*),intent(in) :: filename
    integer :: i1, i2, i3, i4, i5, imax(5)

    imax(1) = size(f,1) 
    imax(2) = size(f,2) 
    imax(3) = size(f,3) 
    imax(4) = size(f,4) 
    imax(5) = size(f,5) 

    open(100,file=filename)
    write(100,'(I5,",")',advance='no')imax(1)
    write(100,'(I5,",")',advance='no')imax(2)
    write(100,'(I5,",")',advance='no')imax(3)
    write(100,'(I5,",")',advance='no')imax(4)
    write(100,'(I5)')imax(5)
    do i5 = 1, imax(5)
      do i4 = 1, imax(4)
        do i3 = 1, imax(3)
          do i2 = 1, imax(2)

            do i1 = 1, imax(1)-1
              write(100,'(e18.11,",")',advance='no')write_for_python(f(i1,i2,i3,i4,i5))
            end do

            write(100,'(e18.11)')write_for_python(f(imax(1),i2,i3,i4,i5))

          end do
        end do
      end do
    end do
    close(100)

  end subroutine fptxt5D

  function write_for_python(f) result(g)
    USE fpcomm,ONLY: rkind
    implicit none
    REAL(rkind),intent(in) :: f
    REAL(rkind) :: g

    if ( f /= f ) then
      g = 0.d0
    else if ( f /= 0.d0 .and. ABS(f) < 1.d-99 ) then
      g = 1.d-99
    else if ( f > 1.d99 ) then
      g = 1.d99
    else if ( f < -1.d99 ) then
      g = 1.d-99
    else 
      g = f
    end if

  end function write_for_python

!-----------------------------------------------------------------------------------------------------------
! subroutine fpwrite_diffusi(recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps)
! use fpcomm
! implicit none
! integer::nr,t=0
! REAL(rkind),dimension(nrmax,nsamax),intent(in):: recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps
! character(len=30)::filedeffe,filekeffe,filechieffe,filetempe,filedense,filegammae,filehfe, &
!                        filedeffi,filekeffi,filechieffi,filetempi,filedensi,filegammai,filehfi, &
!                        filedeff_chi,filechi_deff,fileps,filetps

! t=t+1
! ! open(497,file='source.txt')
! ! open(498,file='epower.txt')
! ! open(499,file='ipower.txt')
! ! do nr=1,nrmax
! ! write(497,'(2e15.4)')dfloat(nr)/dfloat(nrmax),tps(nr,1)
! ! write(498,'(2e15.4)')dfloat(nr)/dfloat(nrmax),rspb(nr,1)+rspf(nr,1)+rsps(nr,1)+rspl(nr,1)+rsps_cx(nr,1)
! ! write(499,'(2e15.4)')dfloat(nr)/dfloat(nrmax),rspb(nr,2)+rspf(nr,2)+rsps(nr,2)+rspl(nr,2)+rsps_cx(nr,2)
! ! end do
! ! close(497)
! ! close(498)
! ! close(499)

!     write(filedeffe,'("deffe",i2.2,".txt")')t
!     write(filechieffe,'("chieffe",i2.2,".txt")')t
!     write(filekeffe,'("keffe",i2.2,".txt")')t
!     write(filetempe,'("tempe",i2.2,".txt")')t
!     write(filedense,'("dense",i2.2,".txt")')t
!     write(filegammae,'("gammae",i2.2,".txt")')t
!     write(filehfe,'("hfe",i2.2,".txt")')t
!     open(500,file=filedeffe)
!     open(501,file=filechieffe)
!     open(502,file=filekeffe)
!     open(503,file=filetempe)
!     open(504,file=filedense)
!     open(505,file=filegammae)
!     open(506,file=filehfe)
!     write(filedeffi,'("deffi",i2.2,".txt")')t
!     write(filechieffi,'("chieffi",i2.2,".txt")')t
!     write(filekeffi,'("keffi",i2.2,".txt")')t
!     write(filetempi,'("tempi",i2.2,".txt")')t
!     write(filedensi,'("densi",i2.2,".txt")')t
!     write(filegammai,'("gammai",i2.2,".txt")')t
!     write(filehfi,'("hfi",i2.2,".txt")')t
!     open(600,file=filedeffi)
!     open(601,file=filechieffi)
!     open(602,file=filekeffi)
!     open(603,file=filetempi)
!     open(604,file=filedensi)
!     open(605,file=filegammai)
!     open(606,file=filehfi)
!     write(filedeff_chi,'("deff_per_chi",i2.2,".txt")')t
!     write(filechi_deff,'("chi_per_deff",i2.2,".txt")')t
!     open(700,file=filedeff_chi)
!     open(701,file=filechi_deff)
!     write(fileps,'("ps",i2.2,".txt")')t
!     write(filetps,'("tps",i2.2,".txt")')t
!     open(702,file=fileps)
!     open(703,file=filetps)

!     do nr=1,nrmax
!        write(500,'(1e15.4)')recv_d(nr,1)
!        write(501,'(1e15.4)')recv_chi(nr,1)
!        write(502,'(1e15.4)')recv_k(nr,1)
!        write(503,'(1e15.4)')recv_temps(nr,1)
!        write(504,'(1e15.4)')rns(nr,1)
!        write(505,'(1e15.4)')recv_gamma(nr,1)
!        write(506,'(1e15.4)')recv_hf(nr,1)
!        write(600,'(1e15.4)')recv_d(nr,2)
!        write(601,'(1e15.4)')recv_chi(nr,2)
!        write(602,'(1e15.4)')recv_k(nr,2)
!        write(603,'(1e15.4)')recv_temps(nr,2)
!        write(604,'(1e15.4)')rns(nr,2)
!        write(605,'(1e15.4)')recv_gamma(nr,2)
!        write(606,'(1e15.4)')recv_hf(nr,2)
!        write(700,'(1e15.4)')recv_d(nr,2)/(recv_chi(nr,2)+1.d-5)
!        write(701,'(1e15.4)')recv_chi(nr,2)/(recv_d(nr,2)+1.d-5)
!        write(702,'(1e15.4)')recv_ps(nr,2)
!        write(703,'(1e15.4)')tps(nr,2)
!     end do

!     close(500)
!     close(501)
!     close(502)
!     close(503)
!     close(504)
!     close(505)
!     close(506)
!     close(600)
!     close(601)
!     close(602)
!     close(603)
!     close(604)
!     close(605)
!     close(606)
!     close(700)
!     close(701)
!     close(702)
!     close(703)

! end subroutine fpwrite_diffusi
! !-----------------------------------------------------------------------------------------------------------
! subroutine fpcsv(recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps,recv_p)
!   use fpcomm
!   implicit none

!   integer::nr
!   REAL(rkind),dimension(nrmax,nsamax),intent(in):: recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps,recv_p
!   character(len=30)::filename

!   write(filename,'("pdep_data",i2,".csv")')int(pdep_exp*10)
!   open(800,file=filename)
!     do nr=1,nrmax
!        write(800,*)real(nr,8)/real(nrmax,8),',',rns(nr,1),',', recv_temps(nr,1),',',recv_d(nr,1),',',recv_chi(nr,1),',',&
!             recv_k(nr,1),',',recv_gamma(nr,1),',',recv_hf(nr,1),',',tps(nr,1),',',recv_ps(nr,1),',',recv_p(nr,1),',',&
!             rns(nr,2),',', recv_temps(nr,2),',',recv_d(nr,2),',',recv_chi(nr,2),',', &
!             recv_k(nr,2),',',recv_gamma(nr,2),',',recv_hf(nr,2),',', tps(nr,2),',',recv_ps(nr,2),',',recv_p(nr,2)
!     end do
!   close(800)
! end subroutine fpcsv
! !----------------------------------------------------------------------------

end module
