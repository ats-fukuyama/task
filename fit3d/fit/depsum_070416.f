c     
c==================================================c
c==================================================c
c      program depsum
       subroutine depsum 
c
c==================================================c
c
c                         coded by S. Murakami
c
c                       last modified 2001/08/20
c==================================================c
c   05/03/23    n, T limit is added.
c   07/12/21    subroutine by nec
c==================================================c
c
      implicit real*8 (a-h,o-y)
c
      parameter (maxdv=10000)
      parameter	(MAXP=30)
      character *256	DBFILE
c
      dimension a1(17),a2(17),a3(17),a4(17)
c
c
!      open(10,file='fit_pb4.out10',status='unknown')
!      open(20,file='fit_pb4.out20',status='unknown')
!      open(30,file='fit_pb4.out30',status='unknown')
!      open(40,file='fit_pb4.out40',status='unknown')

c
      np=maxp-1
c       
c
      do 100 ir=1,np
c
c
c--------------------------------------------c
c--------------------------------------------c
c
c
C       read(10,7000) ii,(a1(j),j=1,17)
C       read(20,7000) ii,(a2(j),j=1,17)
C       read(30,7000) ii,(a3(j),j=1,17)
       read(60,7000) ii,(a1(j),j=1,17)
       read(70,7000) ii,(a2(j),j=1,17)
       read(80,7000) ii,(a3(j),j=1,17)
c
          a4(1)=a1(1)
       do j=2,17
          a4(j)=a1(j)+a2(j)+a3(j)
       end do
c
C       write(40,7000) ir,(a4(j),j=1,17)
       write(90,7000) ir,(a4(j),j=1,17)
c     
 100  continue
c
c     
 7000 format(i4,17e15.5)
c
c     
c      stop
      return
      end subroutine
