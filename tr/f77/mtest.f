c
      implicit none
      integer,parameter :: nrm=10
      real*8,dimension(nrm) :: datam,datam1
      real*8,dimension(nrm+1) :: datag
      integer nrmax,nr,ierr
      real*8 dr
c
      nrmax=10
      dr=1.d0/nrmax
      do nr=1,nrmax
         datam(nr)=(nr-0.5d0)*dr
      enddo
c
      call mesh_convert_mtog(datam,datag,nrmax)
      call mesh_convert_gtom(datag,datam1,nrmax)
c
      do nr=1,nrmax
         write(6,'(I5,1P3E12.4)') nr,datam(nr),datag(nr),datam1(nr)
      enddo
      stop
      end
c
c     ----- convert half mesh to origin + grid mesh -----
c
      subroutine mesh_convert_mtog(datam,datag,nrmax)
c
      implicit none
      integer nrmax
      real*8 datam(nrmax),datag(nrmax+1)
      integer nr
c
      datag(1)=(9.d0*datam(1)-datam(2))/8.d0
      do nr=2,nrmax
         datag(nr)=(datam(nr-1)+datam(nr))/2.d0
      enddo
      datag(nrmax+1)=(4.d0*datam(nrmax)-datam(nrmax-1))/3.d0
      return
      end subroutine mesh_convert_mtog
c
c     ----- convert origin + grid mesh to half mesh -----
c
      subroutine mesh_convert_gtom(datag,datam,nrmax)
c
      implicit none
      integer nrmax
      real*8 datag(nrmax+1),datam(nrmax)
      real*8 c11,c12,c21,c22,det,a11,a12,a21,a22
      integer nr,ierr
c
      c11= 9.d0/8.d0
      c12=-1.d0/8.d0
      c21= 0.5d0
      c22= 0.5d0
      det=c11*c22-c12*c21
      a11= c22/det
      a12=-c12/det
      a21=-c21/det
      a22= c11/det
      datam(1)=a11*datag(1)+a12*datag(2)
      datam(2)=a21*datag(1)+a22*datag(2)
      do nr=3,nrmax
         datam(nr)=2.d0*datag(nr)-datam(nr-1)
      enddo
      return
      end subroutine mesh_convert_gtom
