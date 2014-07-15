!     
!==================================================c
!==================================================c
!      program depsum
!
!==================================================c
!
!                         coded by S. Murakami
!
!                       last modified 2001/08/20
!c==================================================c
!c   05/03/23    n, T limit is added.
!c==================================================c
!c
     subroutine depsum
      implicit real*8 (a-h,o-y)
!
      parameter (maxdv=10000)
      parameter	(MAXP=30)
      character *256	DBFILE
!
      dimension a1(19),a2(19),a3(19),a4(19)
!
!
      open(10,file='fit_pb4.out10',status='old')
      open(20,file='fit_pb4.out20',status='old')
      open(30,file='fit_pb4.out30',status='old')
      open(40,file='fit_pb4.out40',status='new')
      rewind(10)
      rewind(20)
      rewind(30)

!
      np=maxp-1
!--------------------------------------------c
!--------------------------------------------c
      do ir=1,np

       read(10,'(i4,19e15.5)') ii,(a1(j),j=1,19)
       read(20,'(i4,19e15.5)') ii,(a2(j),j=1,19)
       read(30,'(i4,19e15.5)') ii,(a3(j),j=1,19)
!
          a4(1)=a1(1)
       do j=2,19
          a4(j)=a1(j)+a2(j)+a3(j)
       end do
!
       write(40,'(i4,19e15.5)') ir,(a4(j),j=1,19)
!
       enddo
    end subroutine