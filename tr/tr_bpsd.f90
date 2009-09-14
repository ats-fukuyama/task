!     $Id$
      module tr_bpsd

      use bpsd
      type(bpsd_device_type),  private,save :: device
      type(bpsd_species_type), private,save :: species
      type(bpsd_equ1D_type),   private,save :: equ1D
      type(bpsd_metric1D_type),private,save :: metric1D
      type(bpsd_plasmaf_type), private,save :: plasmaf
      logical, private, save :: tr_bpsd_init_flag = .TRUE.
      public

      contains

!=======================================================================
      subroutine tr_bpsd_init
!=======================================================================
      use trcomm
! local variables
      real(8)    :: qpl,rgl
      integer(4) :: ns,nr,ierr
      real(8)    :: temp(nrmp,nsm,3)
!=======================================================================

      if(tr_bpsd_init_flag) then
         species%nsmax=0
         plasmaf%nsmax=0
         plasmaf%nrmax=0
         tr_bpsd_init_flag=.FALSE.
      endif

!      write(6,*) 'top of tr_bpsd_init:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_init:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      device%rr=RR
      device%zz=0.d0
      device%ra=RA
      device%rb=RA+0.1D0
      device%bb=BB
      device%ip=RIP
      device%elip=RKAP
      device%trig=RDLT
      call bpsd_set_data(device,ierr)

      if(species%nsmax.ne.nsmax) then
         if(associated(species%data)) then
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
      call bpsd_set_data(species,ierr)

      if((equ1D%nrmax.ne.nrmax+1)) then
         if(associated(equ1D%rho)) then
            deallocate(equ1D%rho)
         endif
         if(associated(equ1D%data)) then
            deallocate(equ1D%data)
         endif
         equ1D%nrmax=nrmax+1
         allocate(equ1D%rho(equ1D%nrmax))
         allocate(equ1D%data(equ1D%nrmax))
      endif

      equ1D%time=0.d0
      equ1D%rho(1)=0.d0
      do nr=1,nrmax
         equ1D%rho(nr+1)=rg(nr)
      enddo

      if((metric1D%nrmax.ne.nrmax+1)) then
         if(associated(metric1D%rho)) then
            deallocate(metric1D%rho)
         endif
         if(associated(metric1D%data)) then
            deallocate(metric1D%data)
         endif
         metric1D%nrmax=nrmax+1
         allocate(metric1D%rho(metric1D%nrmax))
         allocate(metric1D%data(metric1D%nrmax))
      endif

      metric1D%time=0.d0
      metric1D%rho(1)=0.d0
      do nr=1,nrmax
         metric1D%rho(nr+1)=rg(nr)
      enddo

      if((plasmaf%nsmax.ne.nsmax).or. &
     &   (plasmaf%nrmax.ne.nrmax+1)) then
         if(associated(plasmaf%rho)) then
            deallocate(plasmaf%rho)
         endif
         if(associated(plasmaf%data)) then
            deallocate(plasmaf%data)
         endif
         if(associated(plasmaf%qinv)) then
            deallocate(plasmaf%qinv)
         endif
         plasmaf%nsmax=nsmax
         plasmaf%nrmax=nrmax+1
         allocate(plasmaf%rho(plasmaf%nrmax))
         allocate(plasmaf%data(plasmaf%nrmax,plasmaf%nsmax))
         allocate(plasmaf%qinv(plasmaf%nrmax))
      endif

      plasmaf%time=0.d0
      plasmaf%rho(1)=0.d0
      do nr=1,nrmax
         plasmaf%rho(nr+1)=rg(nr)
      enddo
      do ns=1,nsmax
         call mesh_convert_mtog(rn(1:nrmax,ns),temp(1:nrmax,ns,1),nrmax)
         call mesh_convert_mtog(rt(1:nrmax,ns),temp(1:nrmax,ns,2),nrmax)
         call mesh_convert_mtog(ru(1:nrmax,ns),temp(1:nrmax,ns,3),nrmax)
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
!!      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) &
!!           &              /(4.D0*PI**2*RDP(1:NRMAX))
      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      do nr=2,plasmaf%nrmax
         plasmaf%qinv(nr)=1.d0/qp(nr-1)
      enddo
      plasmaf%qinv(1)=((plasmaf%rho(3))**2*plasmaf%qinv(2) &
     &                -(plasmaf%rho(2))**2*plasmaf%qinv(3)) &
     &               /((plasmaf%rho(3))**2-(plasmaf%rho(2))**2)

      call bpsd_set_data(plasmaf,ierr)
      return
      end subroutine tr_bpsd_init

!=======================================================================
      subroutine tr_bpsd_set(ierr)
!=======================================================================
      use trcomm
      integer(4) :: ierr
! local variables
      integer(4) :: ns,nr
      real(8)    :: temp(nrmp,nsm,3)
!=======================================================================

!      write(6,*) 'top of tr_bpsd_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_set: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      device%rr=RR
      device%zz=0.d0
      device%ra=RA
      device%rb=RA+0.1D0
      device%bb=BB
      device%ip=RIP
      device%elip=RKAP
      device%trig=RDLT
      call bpsd_set_data(device,ierr)

      plasmaf%time=t
      do ns=1,nsmax
         call mesh_convert_mtog(rn(1:nrmax,ns),temp(1:nrmax,ns,1),nrmax)
         call mesh_convert_mtog(rt(1:nrmax,ns),temp(1:nrmax,ns,2),nrmax)
         call mesh_convert_mtog(ru(1:nrmax,ns),temp(1:nrmax,ns,3),nrmax)
      enddo
      do nr=1,plasmaf%nrmax
         if(nr.eq.1) then
            plasmaf%rho(nr)=0.d0
         else
            plasmaf%rho(nr)=rg(nr-1)
         endif
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn=temp(nr,ns,1)*1.d20
            plasmaf%data(nr,ns)%pt=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpr=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%ptpp=temp(nr,ns,2)*1.D3
            plasmaf%data(nr,ns)%pu=temp(nr,ns,3)
         enddo
      enddo

!      write(6,*) 'RDP_set:',(rdp(nr),nr=1,nrmax)

      do nr=1,nrmax
!!         plasmaf%qinv(nr+1) &
!!         & =(4.D0*PI**2*RDP(nr))/(TTRHOG(nr)*ARRHOG(nr)*DVRHOG(nr))
         plasmaf%qinv(nr+1)=(4.D0*PI**2*RDPVRHOG(nr))/(TTRHOG(nr)*ARRHOG(nr))
      enddo
      plasmaf%qinv(1)=((plasmaf%rho(3))**2*plasmaf%qinv(2) &
     &                -(plasmaf%rho(2))**2*plasmaf%qinv(3)) &
     &               /((plasmaf%rho(3))**2-(plasmaf%rho(2))**2)

!      write(6,*) 'end of tr_bpsd_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'end of tr_bpsd_set:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      call bpsd_set_data(plasmaf,ierr)
      return
      end subroutine tr_bpsd_set

!=======================================================================
      subroutine tr_bpsd_get(ierr)
!=======================================================================
      use trcomm
      integer(4),intent(out) :: ierr
! local variables
      integer(4) :: ns,nr
      real(8)    :: temp(nrmp,nsm,3),tmp(nrm)
      real(8)    :: tempx(nrmp,17),psita,dpsitdrho,dvdrho
      REAL(8)    :: FACTOR0, FACTORM, FACTORP
!=======================================================================
!      write(6,*) 'top of tr_bpsd_get: qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_get: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      call bpsd_get_data(device,ierr)

      RR=device%rr
      RA=device%ra
      RB=device%rb
      BB=device%bb
      RIP=device%ip
      RKAP=device%elip
      RDLT=device%trig

      if(modelg.eq.3.or.modelg.eq.5) then
         RIPS=RIP
         RIPE=RIP
      endif

      call bpsd_get_data(plasmaf,ierr)

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
      Q0=2.d0*qp(1)-qp(2)

      do ns=1,nsmax
         call mesh_convert_gtom(temp(1:nrmax,ns,1),rn(1:nrmax,ns),nrmax)
         call mesh_convert_gtom(temp(1:nrmax,ns,2),rt(1:nrmax,ns),nrmax)
         call mesh_convert_gtom(temp(1:nrmax,ns,3),ru(1:nrmax,ns),nrmax)
      enddo

      if(modelg.eq.3.or.modelg.eq.5.or.modelg.eq.8.or.modelg.eq.9) then

      call bpsd_get_data(equ1D,ierr)

      do nr=1,equ1D%nrmax
         tempx(nr,1)=equ1D%data(nr)%psit
         tempx(nr,2)=equ1D%data(nr)%psip
         tempx(nr,3)=equ1D%data(nr)%ppp
         tempx(nr,4)=equ1D%data(nr)%piq
         tempx(nr,5)=equ1D%data(nr)%pip*rmu0/(2.d0*pi)
         tempx(nr,6)=equ1D%data(nr)%pit
      enddo
      call data_interpolate_gtom_full(tempx(1,1),PSITRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,2),PSIPRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,3),PPPRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,4),PIQRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,5),TTRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,6),PIRHO,nrmax)

      TTRHOG(1:nrmax)=tempx(2:nrmax+1,5)
      psita=equ1D%data(equ1D%nrmax)%psit

      call bpsd_get_data(metric1D,ierr)
      do nr=2,metric1D%nrmax       ! metric1D%nrmax = nrmax + 1
         rgl=rg(nr-1) ! rgl is equivalent to metric1D%rho(nr)
         dpsitdrho=2.D0*psita*rgl
         dvdrho=metric1D%data(nr)%dvpsit*dpsitdrho
         tempx(nr,1)=dvdrho
         tempx(nr,2)=metric1D%data(nr)%avegrr2
         tempx(nr,3)=metric1D%data(nr)%aver2i
         tempx(nr,4)=metric1D%data(nr)%avegr
         tempx(nr,5)=metric1D%data(nr)%avegr2
         tempx(nr,6)=metric1D%data(nr)%aveb2
         tempx(nr,7)=metric1D%data(nr)%aveb2i
         tempx(nr,8)=metric1D%data(nr)%avegvr2/dvdrho**2
         tempx(nr,9)=metric1D%data(nr)%rs/metric1D%data(nr)%rr
         tempx(nr,10)=metric1D%data(nr)%rr
         tempx(nr,11)=metric1D%data(nr)%rs
         tempx(nr,12)=metric1D%data(nr)%elip
         tempx(nr,13)=1.d0/metric1D%data(nr)%dvpsip/(2.d0*pi)
         tempx(nr,14)=metric1D%data(nr)%avegvr2
         tempx(nr,15)=metric1D%data(nr)%pvol
         tempx(nr,16)=metric1D%data(nr)%psur
         tempx(nr,17)=metric1D%data(nr)%aveb
      enddo
      nr=1 ! equivalent to "rho = 0"
         tempx(nr,1)=0.d0                     ! definition
         tempx(nr,2)=metric1D%data(nr)%avegrr2
         tempx(nr,3)=metric1D%data(nr)%aver2i
         tempx(nr,4)=metric1D%data(nr)%avegr
         tempx(nr,5)=metric1D%data(nr)%avegr2
         tempx(nr,6)=metric1D%data(nr)%aveb2
         tempx(nr,7)=metric1D%data(nr)%aveb2i
         tempx(nr,8)=tempx(2,8)               ! dummy because dvdrho=0.
         tempx(nr,9)=metric1D%data(nr)%rs/metric1D%data(nr)%rr
         tempx(nr,10)=metric1D%data(nr)%rr
         tempx(nr,11)=metric1D%data(nr)%rs
         tempx(nr,12)=metric1D%data(nr)%elip
         tempx(nr,13)=1.d0/metric1D%data(nr)%dvpsip/(2.d0*pi)
         tempx(nr,14)=metric1D%data(nr)%avegvr2

         tempx(nr,15)=metric1D%data(nr)%pvol
         tempx(nr,16)=metric1D%data(nr)%psur
         tempx(nr,17)=metric1D%data(nr)%aveb

      call data_interpolate_gtom_full(tempx(1,1), DVRHO, nrmax)
      call data_interpolate_gtom_full(tempx(1,2), ABRHO, nrmax)
      call data_interpolate_gtom_full(tempx(1,3), ARRHO, nrmax)
      call data_interpolate_gtom_full(tempx(1,4), AR1RHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,5), AR2RHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,12),RKPRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,14),ABVRHO,nrmax)
      call data_interpolate_gtom_full(tempx(1,17),ABB1RHO,nrmax)

      DVRHOG  (1:nrmax)=tempx(2:nrmax+1,1)  ! dV/drho
      ABRHOG  (1:nrmax)=tempx(2:nrmax+1,2)  ! <|grad rho|^2/R^2>; avegrr2
      ARRHOG  (1:nrmax)=tempx(2:nrmax+1,3)  ! <1/R^2>; aver2i
      AR1RHOG (1:nrmax)=tempx(2:nrmax+1,4)  ! <|grad rho|>
      AR2RHOG (1:nrmax)=tempx(2:nrmax+1,5)  ! <|grad rho|^2>
      ABB2RHOG(1:nrmax)=tempx(2:nrmax+1,6)  ! <B^2>
      AIB2RHOG(1:nrmax)=tempx(2:nrmax+1,7)  ! <1/B^2>
      ARHBRHOG(1:nrmax)=tempx(2:nrmax+1,8)  ! not used
      EPSRHO  (1:nrmax)=tempx(2:nrmax+1,9)
      RMJRHO  (1:nrmax)=tempx(2:nrmax+1,10) ! local R
      RMNRHO  (1:nrmax)=tempx(2:nrmax+1,11) ! local r
      RKPRHOG (1:nrmax)=tempx(2:nrmax+1,12) ! local kappa
      RDPVRHOG(1:nrmax)=tempx(2:nrmax+1,13) ! dpsi/dV
      ABVRHOG (1:nrmax)=tempx(2:nrmax+1,14) ! <|grad V|^2/R^2>

      PVOLRHOG(1:nrmax)=tempx(2:nrmax+1,15) ! Plasma volume
      PSURRHOG(1:nrmax)=tempx(2:nrmax+1,16) ! Plasma surface

      do nr=1,nrmax
!         RDP(nr)=TTRHOG(nr)*ARRHOG(nr)*DVRHOG(nr)/(4.D0*PI**2*QP(nr))
         RDP(nr)=DVRHOG(nr)*RDPVRHOG(nr)
      enddo
!      BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

!      write(6,*) 'RDP_get:',(rdp(nr),nr=1,nrmax)

      NR=1
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         AJOH(NR)=AJ(NR)
      DO NR=2,NRMAX
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         AJOH(NR)=AJ(NR)
      ENDDO

!      write(6,*) 'end of tr_bpsd_get: aj'
!      write(6,'(1P5E12.4)') (aj(nr),nr=1,nrmax)
!      write(6,*) 'end of tr_bpsd_get: qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'end of tr_bpsd_get: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      ! Calibration in order to keep consistency between metrics and plasma current
      if(modelg.eq.3.or.modelg.eq.5) then
         RIP  = ABVRHOG(NRMAX)*RDPVRHOG(NRMAX)/(2.D0*PI*RMU0)*1.D-6
         RIPS = RIP
         RIPE = RIP
      endif

      endif

      return
      end subroutine tr_bpsd_get

!     ----- convert half mesh to origin + grid mesh -----
!
!       Suppose that rho derivative of data be zero
!          at the axis and the boundary
!       datag(1) at rho = 0 and datag(nrmax+1) at rho = 1

      subroutine mesh_convert_mtog(datam,datag,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datam(nrmax)
      real(8),    intent(out) :: datag(nrmax+1)

      datag(1)       = (9.d0*datam(1)-datam(2))/8.d0
      datag(2:nrmax) = 0.5d0 * (datam(1:nrmax-1) + datam(2:nrmax))
      datag(nrmax+1) = (4.d0*datam(nrmax)-datam(nrmax-1))/3.d0

      return
      end subroutine mesh_convert_mtog

!     ----- convert origin + grid mesh to half mesh -----
!
!       just invert mesh_convert_gtom
      !!! CAUTION !!!=======================================================!
      !  This routine should be used only in case that one reconverts       !
      !   data converted by "mesh_convert_mtog" routine.                    !
      !  In any other case, one should use "data_interpolate_gtom" routine. !
      !=====================================================================!

      subroutine mesh_convert_gtom(datag,datam,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datag(nrmax+1)
      real(8),    intent(out) :: datam(nrmax)
      real(8) :: c11=9.d0/8.d0,c12=-1.d0/8.d0,c21=0.5d0,c22=0.5d0
      real(8) :: det,a11,a12,a21,a22

      det=c11*c22-c12*c21
      a11= c22/det
      a12=-c12/det
      a21=-c21/det
      a22= c11/det
      datam(1)=a11*datag(1)+a12*datag(2)
      datam(2)=a21*datag(1)+a22*datag(2)
      datam(3:nrmax) = 2.d0 * datag(3:nrmax) - datam(2:nrmax-1)
      return
      end subroutine mesh_convert_gtom

!     ----- convert half mesh to origin + grid mesh -----
!
!       Suppose that data be zero at the axis
!          and rho derivative of data be zero at the boundary

      subroutine mesh_convert_mtog0(datam,datag,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datam(nrmax)
      real(8),    intent(out) :: datag(nrmax+1)

      datag(1)       = 0.d0
      datag(2:nrmax) = 0.5d0 * (datam(1:nrmax-1) + datam(2:nrmax))
      datag(nrmax+1) = (4.d0*datam(nrmax)-datam(nrmax-1))/3.d0

      return
      end subroutine mesh_convert_mtog0

!     ----- convert origin + grid mesh to half mesh -----
!
!       just invert mesh_convert_mtog0
      !!! CAUTION !!!========================================================!
      !  This routine should be used only in case that one reconverts        !
      !   data converted by "mesh_convert_mtog" routine.                     !
      !  In any other case, one should use "data_interpolate_gtom0" routine. !
      !======================================================================!

      subroutine mesh_convert_gtom0(datag,datam,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datag(nrmax+1)
      real(8),    intent(out) :: datam(nrmax)
      real(8) :: c11=9.d0/8.d0,c12=-1.d0/8.d0,c21=0.5d0,c22=0.5d0
      real(8) :: det,a11,a12,a21,a22

      det=c11*c22-c12*c21
      a11= c22/det
      a12=-c12/det
      a21=-c21/det
      a22= c11/det
      datam(1)=0.5d0*datag(2)
      datam(2:nrmax) = 2.d0 * datag(2:nrmax) - datam(1:nrmax-1)
      return
      end subroutine mesh_convert_gtom0

!     ----- interpolate data on half mesh from full data on grid -----

      subroutine data_interpolate_gtom_full(datag,datam,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datag(nrmax+1)
      real(8),    intent(out) :: datam(nrmax)

      datam(1:nrmax) = 0.5d0 * (datag(1:nrmax) + datag(2:nrmax+1))

      return
      end subroutine data_interpolate_gtom_full

!     ----- interpolate data on half mesh
!              from data on grid except the axis -----

      subroutine data_interpolate_gtom(datag,datam,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datag(nrmax+1)
      real(8),    intent(out) :: datam(nrmax)

      ! linear extrapolation
      datam(1) = 1.5d0 * datag(2) - 0.5d0 * datag(3)

      datam(2:nrmax) = 0.5d0 * (datag(2:nrmax) + datag(3:nrmax+1))

      return
      end subroutine data_interpolate_gtom

      end module tr_bpsd
