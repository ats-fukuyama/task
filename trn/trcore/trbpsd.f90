  MODULE trbpsd

      use bpsd
      type(bpsd_device_type),  private,save :: device
      type(bpsd_species_type), private,save :: species
      type(bpsd_equ1D_type),   private,save :: equ1D
      type(bpsd_metric1D_type),private,save :: metric1D
      type(bpsd_plasmaf_type), private,save :: plasmaf
      logical, private, save :: tr_bpsd_init_flag = .TRUE.
      public

      CONTAINS
!=======================================================================
      SUBROUTINE tr_bpsd_init(ierr)

      USE trcomm
! local variables
      INTEGER(ikind),INTENT(out) :: ierr
      INTEGER(ikind) :: ns,nr
!      real(8)    :: temp(nrmp,nsm,3)


      if(tr_bpsd_init_flag) then
         species%nsmax = 0
         plasmaf%nsmax = 0
         plasmaf%nrmax = 0
         tr_bpsd_init_flag = .FALSE.
      endif

!      write(6,*) 'top of tr_bpsd_init:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_init:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      ! --- deivece ---
      device%rr   = RR
      device%zz   = 0.d0
      device%ra   = RA
      device%rb   = RA+0.1D0
      device%bb   = BB
      device%ip   = RIP
      device%elip = RKAP
      device%trig = RDLT
      call bpsd_set_data(device,ierr)

      ! --- species ---
      if(species%nsmax.ne.nsmax) then
         if(associated(species%data)) then
            deallocate(species%data)
         endif
         species%nsmax = nsmax
         allocate(species%data(species%nsmax))
      endif

      do ns=1,species%nsmax
         species%data(ns)%pa  = pa(ns)
         species%data(ns)%pz  = pz(ns)
         species%data(ns)%pz0 = pz(ns)
      enddo
      call bpsd_set_data(species,ierr)

      ! --- equ1D ---
      if((equ1D%nrmax.ne.nrmax)) then
         if(associated(equ1D%rho)) then
            deallocate(equ1D%rho)
         endif
         if(associated(equ1D%data)) then
            deallocate(equ1D%data)
         endif
         equ1D%nrmax = nrmax
         allocate(equ1D%rho(0:equ1D%nrmax))
         allocate(equ1D%data(0:equ1D%nrmax))
      endif

      equ1D%time   = 0.d0
      equ1D%rho(0) = 0.d0
      do nr=0,nrmax
         equ1D%rho(nr) = rg(nr)
      enddo

      ! --- metric1D ---
      if((metric1D%nrmax.ne.nrmax)) then
         if(associated(metric1D%rho)) then
            deallocate(metric1D%rho)
         endif
         if(associated(metric1D%data)) then
            deallocate(metric1D%data)
         endif
         metric1D%nrmax = nrmax
         allocate(metric1D%rho(0:metric1D%nrmax))
         allocate(metric1D%data(0:metric1D%nrmax))
      endif

      metric1D%time   = 0.d0
      metric1D%rho(0) = 0.d0
      do nr=0,nrmax
         metric1D%rho(nr) = rg(nr)
      enddo

      ! --- plasmaf ---
      if((plasmaf%nsmax.ne.nsmax).or. &
     &   (plasmaf%nrmax.ne.nrmax)) then
         if(associated(plasmaf%rho)) then
            deallocate(plasmaf%rho)
         endif
         if(associated(plasmaf%data)) then
            deallocate(plasmaf%data)
         endif
         if(associated(plasmaf%qinv)) then
            deallocate(plasmaf%qinv)
         endif
         plasmaf%nsmax = nsmax
         plasmaf%nrmax = nrmax
         allocate(plasmaf%rho(0:plasmaf%nrmax))
         allocate(plasmaf%data(0:plasmaf%nrmax,plasmaf%nsmax))
         allocate(plasmaf%qinv(0:plasmaf%nrmax))
      endif

      plasmaf%time   = 0.d0
      plasmaf%rho(0) = 0.d0
      do nr=0,nrmax
         plasmaf%rho(nr) = rg(nr)
      enddo

      do nr=0,plasmaf%nrmax
         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn   = rn(ns,nr)*1.d20
            plasmaf%data(nr,ns)%pt   = rt(ns,nr)*1.d3
            plasmaf%data(nr,ns)%ptpr = rt(ns,nr)*1.d3
            plasmaf%data(nr,ns)%ptpp = rt(ns,nr)*1.d3
            plasmaf%data(nr,ns)%pu   = ru(ns,nr)
         enddo
      enddo
!!      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) &
!!           &              /(4.D0*PI**2*RDP(1:NRMAX))
!!      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      do nr=0,plasmaf%nrmax
         plasmaf%qinv(nr)=1.d0/qp(nr)
      enddo
!      plasmaf%qinv(0)=((plasmaf%rho(2))**2*plasmaf%qinv(1)       &
!     &                -(plasmaf%rho(1))**2*plasmaf%qinv(2))      &
!     &               /((plasmaf%rho(2))**2-(plasmaf%rho(1))**2)

      call bpsd_set_data(plasmaf,ierr)
      return
      END SUBROUTINE tr_bpsd_init

!=======================================================================
      SUBROUTINE tr_bpsd_set(ierr)
!=======================================================================
      USE trcomm
      INTEGER(ikind),INTENT(out) :: ierr
! local variables
      INTEGER(ikind) :: ns,nr
!      real(8)    :: temp(nrmp,nsm,3)
!=======================================================================

!      write(6,*) 'top of tr_bpsd_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_set: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      ! --- device ---
      device%rr   = RR
      device%zz   = 0.d0
      device%ra   = RA
      device%rb   = RA+0.1D0
      device%bb   = BB
      device%ip   = RIP
      device%elip = RKAP
      device%trig = RDLT
      call bpsd_set_data(device,ierr)

      ! --- plasmaf ---
      plasmaf%time = t

      do nr=0,plasmaf%nrmax
         plasmaf%rho(nr)=rg(nr)

         do ns=1,plasmaf%nsmax
            plasmaf%data(nr,ns)%pn   = rn(ns,nr)*1.d20
            plasmaf%data(nr,ns)%pt   = rt(ns,nr)*1.D3
            plasmaf%data(nr,ns)%ptpr = rt(ns,nr)*1.D3
            plasmaf%data(nr,ns)%ptpp = rt(ns,nr)*1.D3
            plasmaf%data(nr,ns)%pu   = ru(ns,nr)
         enddo
      enddo

!      write(6,*) 'RDP_set:',(rdp(nr),nr=1,nrmax)

!      do nr=0,nrmax
!!!         plasmaf%qinv(nr+1) &
!!!         & =(4.D0*PI**2*RDP(nr))/(TTRHOG(nr)*ARRHOG(nr)*DVRHOG(nr))
!         plasmaf%qinv(nr+1)=(4.D0*PI**2*RDPVRHOG(nr+1)) &
!                             / (TTRHOG(nr+1)*ARRHOG(nr+1))
!      enddo
!      plasmaf%qinv(0)=((plasmaf%rho(2))**2*plasmaf%qinv(1) &
!     &                -(plasmaf%rho(1))**2*plasmaf%qinv(2)) &
!     &               /((plasmaf%rho(2))**2-(plasmaf%rho(1))**2)
      do nr=0, nrmax
         plasmaf%qinv(nr) = 1.d0 / qp(nr)
      end do

!      write(6,*) 'end of tr_bpsd_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'end of tr_bpsd_set:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      call bpsd_set_data(plasmaf,ierr)
      return
      END SUBROUTINE tr_bpsd_set

!=======================================================================
      SUBROUTINE tr_bpsd_get(ierr)

      USE trcomm
      INTEGER(ikind),INTENT(out) :: ierr
! local variables
      INTEGER(ikind) :: ns,nr
!      real(8)    :: temp(nrmp,nsm,3)
!      real(8)    :: tempx(nrmp,17),psita,dpsitdrho,dvdrho,rgl
      REAL(rkind)    :: psita
      REAL(rkind)    :: FACTOR0, FACTORM, FACTORP

!      write(6,*) 'top of tr_bpsd_get: qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_get: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      call bpsd_get_data(device,ierr)

      RR   = device%rr
      RA   = device%ra
      BB   = device%bb
      RIP  = device%ip
      RKAP = device%elip
      RDLT = device%trig

!      if(modelg.eq.3.or.modelg.eq.5) then
!         RIPS = RIP
!         RIPE = RIP
!      endif

      call bpsd_get_data(plasmaf,ierr)

      do ns=1,plasmaf%nsmax
         do nr=0,plasmaf%nrmax
            rn(nr,ns) = plasmaf%data(nr,ns)%pn*1.d-20
            rt(nr,ns) = plasmaf%data(nr,ns)%pt*1.D-3
            ru(nr,ns) = plasmaf%data(nr,ns)%pu
         enddo
      enddo
      do nr=0,plasmaf%nrmax
         qp(nr)=1.d0/plasmaf%qinv(nr)
      enddo
      Q0=qp(0)


      ! TASK/EQ or EQDSK output geometry, call TOPICS/EQU or TASK/EQ
      if(modelg.eq.3.or.modelg.eq.5.or.modelg.eq.8.or.modelg.eq.9) then

         equ1D%nrmax = nrmax
         do nr=0,nrmax
            equ1D%rho(nr) = rg(nr)
         enddo
         call bpsd_get_data(equ1D,ierr)

         do nr=0,equ1D%nrmax
            psitrho(nr) = equ1D%data(nr)%psit
            psiprho(nr) = equ1D%data(nr)%psip
            ppprho(nr)  = equ1D%data(nr)%ppp
            piqrho(nr)  = equ1D%data(nr)%piq
            ttrho(nr)   = equ1D%data(nr)%pip * rmu0/(2.d0*pi)
            pirho(nr)   = equ1D%data(nr)%pit
         enddo

         psita=equ1D%data(equ1D%nrmax)%psit

         metric1D%nrmax  = nrmax
         metric1D%rho(0) = 0.d0
         do nr=0,nrmax
            metric1D%rho(nr)=rg(nr)
         enddo

         call bpsd_get_data(metric1D,ierr)
         do nr=0,metric1D%nrmax
            pvolrho(nr) = metric1D%data(nr)%pvol      ! Plasma volume
            psurrho(nr) = metric1D%data(nr)%psur      ! Plasma surface
            dvrho(nr)   = metric1D%data(nr)%dvpsit * (2.d0*psita*rg(nr))
                                                      ! dV/drho
            rdpvrho(nr) = 1.d0 / (metric1D%data(nr)%dvpsip * 2.d0*pi)
                                                      ! dpsi/dV
            !             metric1D%data(nr)%aver2     ! <R^2>
            arrho(nr)   = metric1D%data(nr)%aver2i    ! <1/R^2>
            abb2rho(nr) = metric1D%data(nr)%aveb2     ! <B^2>
            aib2rho(nr) = metric1D%data(nr)%aveb2i    ! <1/B^2>
            !             metric1D%data(nr)%avegv     ! <|grad V|>
            !             metric1D%data(nr)%avegv2    ! <|grad V|^2>
            abvrho(nr)  = metric1D%data(nr)%avegvr2   ! <|grad V|^2/R^2>
            ar1rho(nr)  = metric1D%data(nr)%avegr     ! <|grad rho|>
            ar2rho(nr)  = metric1D%data(nr)%avegr2    ! <|grad rho|^2>
            abrho(nr)   = metric1D%data(nr)%avegrr2   ! <|grad rho|^2/R^2>
            !             metric1D%data(nr)%avegpp2   ! <|grad Psip|^2>
            rmjrho(nr)  = metric1D%data(nr)%rr        ! local R
            rmnrho(nr)  = metric1D%data(nr)%rs        ! local r
            rkprho(nr)  = metric1D%data(nr)%elip      ! elipticity
            !             metric1D%data(nr)%trig      ! triangularity
            abb1rho(nr) = metric1D%data(nr)%aveb      ! <B>

            if(nr /= 0) then
            arhbrho(nr) = abvrho(nr) / dvrho(nr)**2
            endif
            epsrho(nr)  = rmnrho(nr) / rmjrho(nr) ! rs/rr
         enddo
            arhbrho(0) = arhbrho(1)   ! dummy because dvrho(0)=0


! =================== now considering
! === This part should be seperated into another subroutine ??? ===

         do nr=0,nrmax
            !RDP(nr)=TTRHOG(nr)*ARRHOG(nr)*DVRHOG(nr)/(4.D0*PI**2*QP(nr))
            rdp(nr)=dvrho(nr)*rdpvrho(nr)
         enddo
         !      BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR
         !      write(6,*) 'RDP_get:',(rdp(nr),nr=1,nrmax)

!         NR=1
!         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
!         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
!         AJ(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
!         AJOH(NR)=AJ(NR)
!         DO NR=2,NRMAX
!            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
!            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
!            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
!            AJ(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
!            AJOH(NR)=AJ(NR)
!         ENDDO
         
         !      write(6,*) 'end of tr_bpsd_get: aj'
         !      write(6,'(1P5E12.4)') (aj(nr),nr=1,nrmax)
         !      write(6,*) 'end of tr_bpsd_get: qp'
         !      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
         !      write(6,*) 'end of tr_bpsd_get: rt'
         !      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
         !      pause
         
         ! TASK/EQ or EQDSK output geometry
         ! Calibration in order to keep consistency 
         !                          between metrics and plasma current
!         if(modelg.eq.3.or.modelg.eq.5) then
!            RIP  = ABVRHOG(NRMAX)*RDPVRHOG(NRMAX)/(2.D0*PI*RMU0)*1.D-6
!            RIPS = RIP
!            RIPE = RIP
!         endif

      endif

      return
      END SUBROUTINE tr_bpsd_get


!=========================================================================
!=========================================================================
!   ----- convert half mesh to origin + grid mesh -----
!
!  Suppose that rho derivative of data be zero at the axis and the boundary
!  datag(1) at rho = 0 and datag(nrmax+1) at rho = 1

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


!   ----- convert origin + grid mesh to half mesh -----
!      ***  just invert mesh_convert_gtom ***

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


!  ----- convert half mesh to origin + grid mesh -----
!
!  Suppose that data be zero at the axis
!                     and rho derivative of data be zero at the boundary
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
!       ***  just invert mesh_convert_mtog0 ***

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

! ----- interpolate data on half mesh from full data on grid -----
      subroutine data_interpolate_gtom_full(datag,datam,nrmax)

      implicit none
      integer(4), intent(in)  :: nrmax
      real(8),    intent(in)  :: datag(nrmax+1)
      real(8),    intent(out) :: datam(nrmax)

      datam(1:nrmax) = 0.5d0 * (datag(1:nrmax) + datag(2:nrmax+1))

      return
      end subroutine data_interpolate_gtom_full

! ----- interpolate data on half mesh
!                          from data on grid except the axis -----
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

  END MODULE trbpsd
