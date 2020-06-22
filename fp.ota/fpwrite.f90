module fpwrite

contains

  subroutine fpcsv1D(f,filename)

     implicit none

     character(*),intent(in) :: filename
     real(8),  intent(in) :: f(:)
     integer :: imax,i

     imax=size(f)

     open(1985,file=filename)
       do i=1,imax
         if(i==imax)then
           write(1985,'(e18.11)')f(i)
         else
           write(1985,'(e18.11,",")',advance="no")f(i)
         end if
       end do
     close(1985)

  end subroutine fpcsv1D

  subroutine fpcsv2D(f,filename)
    implicit none
    character(*),intent(in) :: filename
    real(8),  intent(in) :: f(:,:)
    integer :: imax,i,jmax,j

    imax=size(f,1)
    jmax=size(f,2)

    open(1985,file=filename)
    do j =1,jmax
     do i=1,imax
       if(i==imax)then
         write(1985,'(e18.11)')f(i,j)
       else
         write(1985,'(e18.11,",")',advance="no")f(i,j)
       end if
     end do
    end do
    close(1985)

  end subroutine fpcsv2D

!-----------------------------------------------------------------------------------------------------------
subroutine fpwrite_diffusi(recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps)
use fpcomm
implicit none
integer::nr,t=0
double precision,dimension(nrmax,nsamax),intent(in):: recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps
character(len=30)::filedeffe,filekeffe,filechieffe,filetempe,filedense,filegammae,filehfe, &
                       filedeffi,filekeffi,filechieffi,filetempi,filedensi,filegammai,filehfi, &
                       filedeff_chi,filechi_deff,fileps,filetps

t=t+1
! open(497,file='source.txt')
! open(498,file='epower.txt')
! open(499,file='ipower.txt')
! do nr=1,nrmax
! write(497,'(2e15.4)')dfloat(nr)/dfloat(nrmax),tps(nr,1)
! write(498,'(2e15.4)')dfloat(nr)/dfloat(nrmax),rspb(nr,1)+rspf(nr,1)+rsps(nr,1)+rspl(nr,1)+rsps_cx(nr,1)
! write(499,'(2e15.4)')dfloat(nr)/dfloat(nrmax),rspb(nr,2)+rspf(nr,2)+rsps(nr,2)+rspl(nr,2)+rsps_cx(nr,2)
! end do
! close(497)
! close(498)
! close(499)

    write(filedeffe,'("deffe",i2.2,".txt")')t
    write(filechieffe,'("chieffe",i2.2,".txt")')t
    write(filekeffe,'("keffe",i2.2,".txt")')t
    write(filetempe,'("tempe",i2.2,".txt")')t
    write(filedense,'("dense",i2.2,".txt")')t
    write(filegammae,'("gammae",i2.2,".txt")')t
    write(filehfe,'("hfe",i2.2,".txt")')t
    open(500,file=filedeffe)
    open(501,file=filechieffe)
    open(502,file=filekeffe)
    open(503,file=filetempe)
    open(504,file=filedense)
    open(505,file=filegammae)
    open(506,file=filehfe)
    write(filedeffi,'("deffi",i2.2,".txt")')t
    write(filechieffi,'("chieffi",i2.2,".txt")')t
    write(filekeffi,'("keffi",i2.2,".txt")')t
    write(filetempi,'("tempi",i2.2,".txt")')t
    write(filedensi,'("densi",i2.2,".txt")')t
    write(filegammai,'("gammai",i2.2,".txt")')t
    write(filehfi,'("hfi",i2.2,".txt")')t
    open(600,file=filedeffi)
    open(601,file=filechieffi)
    open(602,file=filekeffi)
    open(603,file=filetempi)
    open(604,file=filedensi)
    open(605,file=filegammai)
    open(606,file=filehfi)
    write(filedeff_chi,'("deff_per_chi",i2.2,".txt")')t
    write(filechi_deff,'("chi_per_deff",i2.2,".txt")')t
    open(700,file=filedeff_chi)
    open(701,file=filechi_deff)
    write(fileps,'("ps",i2.2,".txt")')t
    write(filetps,'("tps",i2.2,".txt")')t
    open(702,file=fileps)
    open(703,file=filetps)

    do nr=1,nrmax
       write(500,'(1e15.4)')recv_d(nr,1)
       write(501,'(1e15.4)')recv_chi(nr,1)
       write(502,'(1e15.4)')recv_k(nr,1)
       write(503,'(1e15.4)')recv_temps(nr,1)
       write(504,'(1e15.4)')rns(nr,1)
       write(505,'(1e15.4)')recv_gamma(nr,1)
       write(506,'(1e15.4)')recv_hf(nr,1)
       write(600,'(1e15.4)')recv_d(nr,2)
       write(601,'(1e15.4)')recv_chi(nr,2)
       write(602,'(1e15.4)')recv_k(nr,2)
       write(603,'(1e15.4)')recv_temps(nr,2)
       write(604,'(1e15.4)')rns(nr,2)
       write(605,'(1e15.4)')recv_gamma(nr,2)
       write(606,'(1e15.4)')recv_hf(nr,2)
       write(700,'(1e15.4)')recv_d(nr,2)/(recv_chi(nr,2)+1.d-5)
       write(701,'(1e15.4)')recv_chi(nr,2)/(recv_d(nr,2)+1.d-5)
       write(702,'(1e15.4)')recv_ps(nr,2)
       write(703,'(1e15.4)')tps(nr,2)
    end do

    close(500)
    close(501)
    close(502)
    close(503)
    close(504)
    close(505)
    close(506)
    close(600)
    close(601)
    close(602)
    close(603)
    close(604)
    close(605)
    close(606)
    close(700)
    close(701)
    close(702)
    close(703)

end subroutine fpwrite_diffusi
!-----------------------------------------------------------------------------------------------------------
subroutine fpcsv(recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps,recv_p)
  use fpcomm
  implicit none

  integer::nr
  double precision,dimension(nrmax,nsamax),intent(in):: recv_d,recv_chi,recv_k,recv_gamma,recv_hf,recv_temps,recv_ps,recv_p
  character(len=30)::filename

  write(filename,'("pdep_data",i2,".csv")')int(pdep_exp*10)
  open(800,file=filename)
    do nr=1,nrmax
       write(800,*)real(nr,8)/real(nrmax,8),',',rns(nr,1),',', recv_temps(nr,1),',',recv_d(nr,1),',',recv_chi(nr,1),',',&
            recv_k(nr,1),',',recv_gamma(nr,1),',',recv_hf(nr,1),',',tps(nr,1),',',recv_ps(nr,1),',',recv_p(nr,1),',',&
            rns(nr,2),',', recv_temps(nr,2),',',recv_d(nr,2),',',recv_chi(nr,2),',', &
            recv_k(nr,2),',',recv_gamma(nr,2),',',recv_hf(nr,2),',', tps(nr,2),',',recv_ps(nr,2),',',recv_p(nr,2)
    end do
  close(800)
end subroutine fpcsv
!----------------------------------------------------------------------------

end module
