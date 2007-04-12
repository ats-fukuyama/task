!     $Id$
      module trpl_mod
      use bpsd_mod
      type(bpsd_species_type),private,save :: species
      type(bpsd_equ1D_type),private,save :: equ1D
      type(bpsd_metric1D_type),private,save :: metric1D
      type(bpsd_plasmaf_type),private,save :: plasmaf
      logical, private, save :: trpl_init_flag = .TRUE.
      public
      contains

!=======================================================================
      subroutine trpl_init
!=======================================================================
      use trcomm
! local variables
      real(8)       qpl,rgl
      integer(4)    ns,nr,ierr
      real*8 temp(nrmp,nsm,3)
!=======================================================================
      if(trpl_init_flag) then
         species%nsmax=0
         plasmaf%nsmax=0
         plasmaf%nrmax=0
         trpl_init_flag=.FALSE.
      endif

!      write(6,*) 'top of trpl_init:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of trpl_init:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      if(species%nsmax.ne.nsmax) then
         if(allocated(species%data)) then
            deallocate(species%data)
         endif
         species%nsmax=nsmax
         allocate(species%data(species%nsmax))
      endif

      do ns=1,species%nsmax
         species%data(ns)%pa=pa(ns)
         species%data(ns)%pz=pz(ns)
         species%data(ns)%pz0=pz(ns)
      enddo
      call bpsd_set_species(species,ierr)

      if((equ1D%nrmax.ne.nrmax+1)) then
         if(allocated(equ1D%s)) then
            deallocate(equ1D%s)
         endif
         if(allocated(equ1D%data)) then
            deallocate(equ1D%data)
         endif
         equ1D%nrmax=nrmax+1
         allocate(equ1D%s(equ1D%nrmax))
         allocate(equ1D%data(equ1D%nrmax))
      endif

      equ1D%time=0.d0
      equ1D%s(1)=0.d0
      do nr=1,nrmax
         equ1D%s(nr+1)=rg(nr)**2
      enddo

      if((metric1D%nrmax.ne.nrmax+1)) then
         if(allocated(metric1D%s)) then
            deallocate(metric1D%s)
         endif
         if(allocated(metric1D%data)) then
            deallocate(metric1D%data)
         endif
         metric1D%nrmax=nrmax+1
         allocate(metric1D%s(metric1D%nrmax))
         allocate(metric1D%data(metric1D%nrmax))
      endif

      metric1D%time=0.d0
      metric1D%s(1)=0.d0
      do nr=1,nrmax
         metric1D%s(nr+1)=rg(nr)**2
      enddo

      if((plasmaf%nsmax.ne.nsmax).or. &
     &   (plasmaf%nrmax.ne.nrmax+1)) then
         if(allocated(plasmaf%s)) then
            deallocate(plasmaf%s)
         endif
         if(allocated(plasmaf%data)) then
            deallocate(plasmaf%data)
         endif
         if(allocated(plasmaf%qinv)) then
            deallocate(plasmaf%qinv)
         endif
         plasmaf%nsmax=nsmax
         plasmaf%nrmax=nrmax+1
         allocate(plasmaf%s(plasmaf%nrmax))
         allocate(plasmaf%data(plasmaf%nrmax,plasmaf%nsmax))
         allocate(plasmaf%qinv(plasmaf%nrmax))
      endif

      plasmaf%time=0.d0
      plasmaf%s(1)=0.d0
      do nr=1,nrmax
         plasmaf%s(nr+1)=rg(nr)**2
      enddo
      do ns=1,nsmax
         call mesh_convert_mtog(rn(1,ns),temp(1,ns,1),nrmax)
         call mesh_convert_mtog(rt(1,ns),temp(1,ns,2),nrmax)
         call mesh_convert_mtog(ru(1,ns),temp(1,ns,3),nrmax)
      enddo
      do nr=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=temp(nr,ns,1)*1.d20
            plasmaf%data(nr,ns)%pt=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpr=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpp=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%pu=temp(nr,ns,3)
         enddo
      enddo
         QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) &
     &              /(4.D0*PI**2*RDP(1:NRMAX))
      do nr=2,plasmaf%nrmax
         plasmaf%qinv(nr)=1.d0/qp(nr-1)
      enddo
      plasmaf%qinv(1)=(plasmaf%s(3)*plasmaf%qinv(2) &
     &                -plasmaf%s(2)*plasmaf%qinv(3)) &
     &               /(plasmaf%s(3)-plasmaf%s(2))

      call bpsd_set_plasmaf(plasmaf,ierr)
      return
      end subroutine trpl_init

!=======================================================================
      subroutine trpl_set(ierr)
!=======================================================================
      use trcomm
      integer    ierr
! local variables
      integer    ns,nr
      real*8 temp(nrmp,nsm,3)
!=======================================================================

!      write(6,*) 'top of trpl_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of trpl_set: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      plasmaf%time=t
      do ns=1,nsmax
         call mesh_convert_mtog(rn(1,ns),temp(1,ns,1),nrmax)
         call mesh_convert_mtog(rt(1,ns),temp(1,ns,2),nrmax)
         call mesh_convert_mtog(ru(1,ns),temp(1,ns,3),nrmax)
      enddo
      do nr=1,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=temp(nr,ns,1)*1.d20
            plasmaf%data(nr,ns)%pt=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpr=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpp=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%pu=temp(nr,ns,3)
         enddo
      enddo

!      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) &
!     &            /(4.D0*PI**2*RDP(1:NRMAX))

      do nr=2,plasmaf%nrmax
         plasmaf%qinv(nr)=1.d0/qp(nr-1)
      enddo
      plasmaf%qinv(1)=(plasmaf%s(3)*plasmaf%qinv(2) &
     &                -plasmaf%s(2)*plasmaf%qinv(3)) &
     &               /(plasmaf%s(3)-plasmaf%s(2))

!      write(6,*) 'end of trpl_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'end of trpl_set:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      call bpsd_set_plasmaf(plasmaf,ierr)
      return
      end subroutine trpl_set

!=======================================================================
      subroutine trpl_get(ierr)
!=======================================================================
      use trcomm
      integer    ierr
! local variables
      integer    ns,nr
      real*8 temp(nrmp,nsm,3)
      real*8 tempx(nrmp,12),psita,dpsitdrho,dvdrho
      REAL(8) :: FACTOR0, FACTORM, FACTORP
!=======================================================================
!      write(6,*) 'top of trpl_get: qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of trpl_get: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause
     
      call bpsd_get_plasmaf(plasmaf,ierr)

      do ns=1,plasmaf%nsmax
         do nr=1,plasmaf%nrmax
            temp(nr,ns,1)=plasmaf%data(nr,ns)%pn*1.d-20
            temp(nr,ns,2)=plasmaf%data(nr,ns)%pt*1.D-3
            temp(nr,ns,3)=plasmaf%data(nr,ns)%pu
         enddo
      enddo
      do nr=2,plasmaf%nrmax
         qp(nr-1)=1.d0/plasmaf%qinv(nr)
      enddo

      do ns=1,nsmax
         call mesh_convert_gtom(temp(1,ns,1),rn(1,ns),nrmax)
         call mesh_convert_gtom(temp(1,ns,2),rt(1,ns),nrmax)
         call mesh_convert_gtom(temp(1,ns,3),ru(1,ns),nrmax)
      enddo

      call bpsd_get_equ1D(equ1D,ierr)
      do nr=1,equ1D%nrmax
         tempx(nr,1)=equ1D%data(nr)%psip
         tempx(nr,2)=equ1D%data(nr)%psit
         tempx(nr,3)=equ1D%data(nr)%ppp
         tempx(nr,4)=equ1D%data(nr)%piq
         tempx(nr,5)=equ1D%data(nr)%pip*rmu0/(2.d0*pi)
         tempx(nr,6)=equ1D%data(nr)%pit
      enddo
      call mesh_convert_gtom(tempx(1,1),PSITRHO,nrmax)
      call mesh_convert_gtom(tempx(1,2),PSIPRHO,nrmax)
      call mesh_convert_gtom(tempx(1,3),PPPRHO,nrmax)
      call mesh_convert_gtom(tempx(1,4),PIQRHO,nrmax)
      call mesh_convert_gtom(tempx(1,5),TTRHO,nrmax)
      call mesh_convert_gtom(tempx(1,6),PIRHO,nrmax)
      TTRHOG(1:nrmax)=tempx(2:nrmax+1,5)
      psita=equ1D%data(equ1D%nrmax)%psit

      call bpsd_get_metric1D(metric1D,ierr)
      do nr=2,metric1D%nrmax
         rgl=rg(nr-1)
         dpsitdrho=2.D0*psita*rgl
         dvdrho=metric1D%data(nr)%dvpsit*dpsitdrho
         write(26,'(I5,1P3E12.4)')  &
     &        nr,dpsitdrho,metric1D%data(nr)%dvpsit,dvdrho
         tempx(nr,1)=dvdrho
         tempx(nr,2)=metric1D%data(nr)%avegvr2/dvdrho**2
         tempx(nr,3)=metric1D%data(nr)%aver2i
         tempx(nr,4)=sqrt(metric1D%data(nr)%avegv2/dvdrho**2)
         tempx(nr,5)=metric1D%data(nr)%avegv2/dvdrho**2
         tempx(nr,6)=metric1D%data(nr)%aveb2
         tempx(nr,7)=metric1D%data(nr)%aveb2i
         tempx(nr,8)=metric1D%data(nr)%avegvr2/dvdrho**2
         tempx(nr,9)=metric1D%data(nr)%avera/metric1D%data(nr)%averr
         tempx(nr,10)=metric1D%data(nr)%averr
         tempx(nr,11)=metric1D%data(nr)%avera
         tempx(nr,12)=metric1D%data(nr)%aveelip
      enddo
      nr=1
         tempx(nr,1)=0.d0
         tempx(nr,2)=tempx(2,2)
         tempx(nr,3)=metric1D%data(nr)%aver2i
         tempx(nr,4)=tempx(2,4)
         tempx(nr,5)=tempx(2,5)
         tempx(nr,6)=metric1D%data(nr)%aveb2
         tempx(nr,7)=metric1D%data(nr)%aveb2i
         tempx(nr,8)=tempx(2,8)
         tempx(nr,9)=metric1D%data(nr)%avera/metric1D%data(nr)%averr
         tempx(nr,10)=metric1D%data(nr)%averr
         tempx(nr,11)=metric1D%data(nr)%avera
         tempx(nr,12)=metric1D%data(nr)%aveelip

      call mesh_convert_gtom0(tempx(1,1),DVRHO,nrmax)
      call mesh_convert_gtom(tempx(1,2),ABRHO,nrmax)
      call mesh_convert_gtom(tempx(1,3),ARRHO,nrmax)
      call mesh_convert_gtom(tempx(1,4),AR1RHO,nrmax)
      call mesh_convert_gtom(tempx(1,5),AR2RHO,nrmax)
      call mesh_convert_gtom(tempx(1,12),RKPRHO,nrmax)

      DVRHOG(1:nrmax)=tempx(2:nrmax+1,1)
      ABRHOG(1:nrmax)=tempx(2:nrmax+1,2)
      ARRHOG(1:nrmax)=tempx(2:nrmax+1,3)
      AR1RHOG(1:nrmax)=tempx(2:nrmax+1,4)
      AR2RHOG(1:nrmax)=tempx(2:nrmax+1,5)
      ABB2RHOG(1:nrmax)=tempx(2:nrmax+1,6)
      AIB2RHOG(1:nrmax)=tempx(2:nrmax+1,7)
      ARHBRHOG(1:nrmax)=tempx(2:nrmax+1,8)
      EPSRHO(1:nrmax)=tempx(2:nrmax+1,9)
      RMJRHO(1:nrmax)=tempx(2:nrmax+1,10)
      RMNRHO(1:nrmax)=tempx(2:nrmax+1,11)
      RKPRHOG(1:nrmax)=tempx(2:nrmax+1,12)

      RDP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) &
     &            /(4.D0*PI**2*QP(1:NRMAX))
!      BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

      do nr=1,nrmax
         write(26,'(I5,1P5E12.4)') &
     &        nr,qp(nr),TTRHO(nr),ARRHOG(nr),DVRHOG(nr),rdp(nr)
      enddo
      NR=1
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         AJOH(NR)=AJ(NR)
         write(26,'(I5,1P5E12.4)') &
     &        nr,dvrho(nr),factor0,0.d0,factorp,aj(nr)
      DO NR=2,NRMAX
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         AJOH(NR)=AJ(NR)
         write(26,'(I5,1P5E12.4)') &
     &        nr,dvrho(nr),factor0,factorm,factorp,aj(nr)
      ENDDO
!      NR=1
!         FACTOR0=RR/(RMU0*DVRHO(NR))
!         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
!         AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
!      DO NR=2,NRMAX
!         FACTOR0=RR/(RMU0*DVRHO(NR))
!         FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
!         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
!         AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
!      ENDDO

!      write(6,*) 'end of trpl_get: aj'
!      write(6,'(1P5E12.4)') (aj(nr),nr=1,nrmax)
!      write(6,*) 'end of trpl_get: qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'end of trpl_get: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      return
      end subroutine trpl_get

!     ----- convert half mesh to origin + grid mesh -----

      subroutine mesh_convert_mtog(datam,datag,nrmax)

      implicit none
      integer nrmax
      real*8 datam(nrmax),datag(nrmax+1)
      integer nr

      datag(1)=(9.d0*datam(1)-datam(2))/8.d0
      do nr=2,nrmax
         datag(nr)=(datam(nr-1)+datam(nr))/2.d0
      enddo
      datag(nrmax+1)=(4.d0*datam(nrmax)-datam(nrmax-1))/3.d0
      return
      end subroutine mesh_convert_mtog

!     ----- convert origin + grid mesh to half mesh -----

      subroutine mesh_convert_gtom(datag,datam,nrmax)

      implicit none
      integer nrmax
      real*8 datag(nrmax+1),datam(nrmax)
      real*8 c11,c12,c21,c22,det,a11,a12,a21,a22
      integer nr,ierr

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

!     ----- convert half mesh to origin + grid mesh -----

      subroutine mesh_convert_mtog0(datam,datag,nrmax)

      implicit none
      integer nrmax
      real*8 datam(nrmax),datag(nrmax+1)
      integer nr

      datag(1)=0.d0
      do nr=2,nrmax
         datag(nr)=(datam(nr-1)+datam(nr))/2.d0
      enddo
      datag(nrmax+1)=(4.d0*datam(nrmax)-datam(nrmax-1))/3.d0
      return
      end subroutine mesh_convert_mtog0

!     ----- convert origin + grid mesh to half mesh -----

      subroutine mesh_convert_gtom0(datag,datam,nrmax)

      implicit none
      integer nrmax
      real*8 datag(nrmax+1),datam(nrmax)
      real*8 c11,c12,c21,c22,det,a11,a12,a21,a22
      integer nr,ierr

      c11= 9.d0/8.d0
      c12=-1.d0/8.d0
      c21= 0.5d0
      c22= 0.5d0
      det=c11*c22-c12*c21
      a11= c22/det
      a12=-c12/det
      a21=-c21/det
      a22= c11/det
      datam(1)=0.5d0*datag(2)
      do nr=2,nrmax
         datam(nr)=2.d0*datag(nr)-datam(nr-1)
      enddo
      return
      end subroutine mesh_convert_gtom0

      end module trpl_mod
