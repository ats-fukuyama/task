module torbeam_iwrap

  contains

!-------------------------------------------------------------------------------------
subroutine torbeam(equilibrium_ids,core_profiles_ids,ec_launchers_ids,waves_ids, &
     codeparam_torbeam,output_flag,output_message)
  !----------------------------
  ! TORBEAM subroutine in IMAS
  !----------------------------
  use ids_schemas
  use ids_routines
  use mod_f90_kind
  use f90_file_reader, only: file2buffer
  use codeparam_torbeam_ids
  use interpos_module

  implicit none

  double precision, parameter :: pi = 3.14159265358979

  type(ids_equilibrium), intent(in):: equilibrium_ids
  type(ids_core_profiles), intent(in) :: core_profiles_ids
  type(ids_ec_launchers), intent(in):: ec_launchers_ids
  type(ids_waves), intent(out) :: waves_ids
  type(ids_parameters_input):: codeparam_torbeam
  type(type_torbeam_codeparam_data):: codeparam_torbeam_data
  integer, intent(out):: output_flag
  character(len=:), pointer, intent(out) :: output_message

  double precision:: sgnm,rhoresult(0:19)
  integer:: npsi,npsi2,k,kp,lfd,iray,icnt,ibgout

  ! Definition of intinbeam and flotinbeam input structures
  integer:: maxint,maxflt,iii
  parameter (maxint = 50, maxflt = 50)
  integer:: intinbeam(0:maxint-1)
  double precision:: floatinbeam(0:maxflt-1)

  integer:: mmax,nmax,maxdiml,maxlen,maxvol
  parameter (mmax = 450, nmax = mmax)
  parameter (maxdiml = 1+mmax+nmax+4*mmax*nmax)
  parameter (maxlen = 2*mmax+2*nmax)
  parameter (maxvol = 100)
  double precision:: xtordeg,xpoldeg
  double precision:: eqdata(0:maxdiml-1)
  double precision:: prdata(0:maxlen-1)
  double precision:: zeffdata(0:2*mmax-1),zeffdat(0:2*mmax-1)
  double precision, dimension(2*maxvol) :: volprof

  ! Originally from topfilevals
  double precision, dimension(:), allocatable :: psi,Rarr,Zarr
  double precision, dimension(:,:), allocatable :: br,bt,bz,psi2d
  integer:: ni,nj,i,j
  double precision:: psiax,psiedge
  character(len=70) :: line

  ! Originally from profilesvals
  double precision, dimension(:), allocatable:: te,ne,Zeff

  ! Originally from antennavals
  integer:: noout, ncdroutine, nabsroutine
  integer:: nprofv,nastra,nprofcalc,ncdharm,npnts_extrap,nfreq_extrap,nrel
  integer:: npow, ncd, ianexp, ndns, nte, nshot, nrela, nmaxh
  double precision:: alpha,beta, rhostop,xzsrch
  double precision:: xrtol, xatol, xstep, xthdeg, xphdeg
  double precision:: nefac, tefac, btfac, bpfac

  integer:: ibeam,nbeam,iend,kend
  integer:: ndat,npnt,ntraj,nfreq,nlines_param
  parameter (ndat = 100000, npnt = 5000, ntraj = 10000)
  integer, allocatable :: npointsout(:)
  double precision, allocatable :: extrascal(:,:)
  double precision:: t1data(0:6*ndat-1),t1tdata(0:6*ndat-1)
  double precision:: t2data(0:5*ndat-1),t2ndata(0:3*npnt-1)
  double precision, allocatable :: profout(:,:,:),trajout(:,:,:),extradata(:,:,:)
  double precision:: proftot(3,npnt)
  double precision:: power_array(ndat),npar_array(ndat),nper_array(ndat),slength_array(ndat)
  double precision:: b0,iplasma,powout,curr,rhcd,jcdout,wcdout

  interface
     subroutine beam(intinb,floatinb,mmm,nnn,eqdat, &
          kkk,lll,prdat,rhores,ie,t1dat,t1tdat, &
          ke,t2dat,t2ndat,ic,ib,nprofv,volpr, &
          power_array,npar_array,nper_array,slength_array,powout) &
          bind(c, name="beam_")
       use mod_f90_kind
       implicit none
       integer(ikind) :: maxint,maxflt, mmax, nmax, ndat
       integer(ikind) :: maxvol, maxdim, maxlen, npnt
       parameter (maxint = 50, maxflt = 50)
       parameter (mmax = 450, nmax = mmax)
       parameter (ndat = 100000)
       parameter (maxvol = 100)
       parameter (maxdim = 1+mmax+nmax+4*mmax*nmax)
       parameter (maxlen = 2*mmax+2*nmax)
       parameter (npnt = 5000)
       integer(ikind) :: intinb(0:maxint-1)
       integer(ikind) :: ie,ke,ic,ib
       integer(ikind) :: nprofv
       integer(ikind) :: mmm,nnn,kkk,lll
       real(rkind) :: floatinb(0:maxflt-1)
       real(rkind) :: rhores(0:19)
       real(rkind), dimension(0:maxdim-1) :: eqdat
       real(rkind), dimension(0:maxlen-1) :: prdat
       real(rkind), dimension(0:6*ndat-1) :: t1dat,t1tdat
       real(rkind), dimension(0:5*ndat-1) :: t2dat
       real(rkind), dimension(0:3*npnt-1) :: t2ndat
       real(rkind), dimension(2*maxvol) :: volpr
       real(rkind), dimension(ndat):: power_array,npar_array,nper_array,slength_array
       real(rkind) :: powout
     end subroutine beam
  end interface

  common /eq/ eqdata
  common /pr/ prdata
  common /intpr/ npsi,npsi2
  common /effch/ zeffdat
  common /tmp/ psiedge,psiax
  common /datout/ curr,rhcd,jcdout,wcdout

  ! GLOBAL ERROR FLAG
  output_flag = 0
  allocate(character(200):: output_message)
  output_message = ''

  ! ... Initialize vector 'eqdata' as TORBEAM input (topfile):
  ! Psi, 1D, for normalization
  npsi = size(core_profiles_ids%profiles_1d(1)%grid%psi)
  npsi2 = npsi ! just for the common intpr not to complain, but otherwise useless
  allocate(psi(npsi))
  ! psi = core_profiles_ids%profiles_1d(1)%grid%psi
  ! To ensure consistency between the 1D and 2D psi profiles: take both from the equilibrium IDS
  call interpos(equilibrium_ids%time_slice(1)%profiles_1d%rho_tor_norm, &
       equilibrium_ids%time_slice(1)%profiles_1d%psi, &
       size(equilibrium_ids%time_slice(1)%profiles_1d%rho_tor_norm), &
       size(core_profiles_ids%profiles_1d(1)%grid%psi), &
       xout=core_profiles_ids%profiles_1d(1)%grid%rho_tor_norm, &
       yout=psi)

  ! Generate grid quantities
  ni = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1)
  nj = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2)
  allocate(Rarr(ni),Zarr(nj))
  allocate(br(ni,nj))
  allocate(bt(ni,nj))
  allocate(bz(ni,nj))
  allocate(psi2d(ni,nj))
  Rarr    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1
  Zarr    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2
  br      = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_r
  bt      = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_tor
  bz      = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_z
  psi2d   = equilibrium_ids%time_slice(1)%profiles_2d(1)%psi
  psiedge = equilibrium_ids%time_slice(1)%global_quantities%psi_boundary
  psiax   = equilibrium_ids%time_slice(1)%global_quantities%psi_axis

  ! ---------------------------------------
  ! FILL TORBEAM INTERNAL EQUILIBRIUM DATA
  ! ---------------------------------------
  eqdata = 0.
  eqdata(0) = psiedge
  do i=1,ni
     eqdata(i) = Rarr(i)
  enddo
  do j=1,nj
     eqdata(j+ni) = Zarr(j)
  enddo
  k = 0
  do j=1,nj
     do i=1,ni
        k = k+1
        eqdata(k+ni+nj) = br(i,j)
     enddo
  enddo
  k = 0
  do j=1,nj
     do i=1,ni
        k = k+1
        eqdata(k+ni+nj+ni*nj) = bt(i,j)
     enddo
  enddo
  k = 0
  do j=1,nj
     do i=1,ni
        k = k+1
        eqdata(k+ni+nj+2*ni*nj) = bz(i,j)
     enddo
  enddo
  k = 0
  do j=1,nj
     do i=1,ni
        k = k+1
        eqdata(k+ni+nj+3*ni*nj) = psi2d(i,j)
     enddo
  enddo

  b0      = -equilibrium_ids%time_slice(1)%global_quantities%magnetic_axis%b_field_tor
  iplasma = equilibrium_ids%time_slice(1)%global_quantities%ip

  ! -----------------------------------
  ! FILL TORBEAM INTERNAL PROFILE DATA
  ! -----------------------------------
  ! ... Initialize vector 'prdata' as TORBEAM input (ne.dat & Te.dat):
  prdata = 0.
  ! Psi and profiles
  allocate(te(npsi),ne(npsi),Zeff(npsi))
  te    = core_profiles_ids%profiles_1d(1)%electrons%temperature
  ne    = core_profiles_ids%profiles_1d(1)%electrons%density
  Zeff  = core_profiles_ids%profiles_1d(1)%zeff

  do i=1,npsi
     j = i-1
     prdata(j) = sqrt((psi(i)-psi(1))/(psi(npsi)-psi(1)))
     prdata(j+npsi) = ne(i)*1.d-19
  enddo

  do i=1,npsi
     j = i-1
     prdata(j+2*npsi) = sqrt((psi(i)-psi(1))/(psi(npsi)-psi(1)))
     prdata(j+2*npsi+npsi) = te(i)*1.d-3
  enddo

  zeffdat = 0.
  do i=1,npsi
     j = i-1
     zeffdat(j) = sqrt((psi(i)-psi(1))/(psi(npsi)-psi(1)))
     zeffdat(j+npsi) = Zeff(i)
  enddo

  zeffdata = 0.
  do lfd=0,2*npsi-1
     zeffdata(lfd) = zeffdat(lfd)
  enddo

  ! COUND NUMBER OF INPUT BEAMS (EVEN THOSE WITH NO POWER)
  if (associated(ec_launchers_ids%launcher)) then
     nbeam = size(ec_launchers_ids%launcher)
  else
     nbeam = 0
  endif

  ! ALLOCATIONS
  allocate(profout(nbeam,3,npnt))
  allocate(trajout(nbeam,15,ntraj))
  allocate(npointsout(nbeam))
  allocate(extrascal(nbeam,3))      ! increase the second dimension if there are more data to be saved
  allocate(extradata(nbeam,4,ndat)) ! increase the second dimension if there are more data to be saved
  profout    = 0.
  trajout    = 0.
  npointsout = 0
  extrascal  = 0.
  extradata  = 0.

  ! READ INPUT XML FILE
  call parse_torbeam_codeparam(codeparam_torbeam%parameters_value,codeparam_torbeam_data)

  ! Initialize antenna data as TORBEAM input:

  !Int parameters
  npow         = codeparam_torbeam_data%npow
  ncd          = codeparam_torbeam_data%ncd
  ncdroutine   = codeparam_torbeam_data%ncdroutine
  nprofv       = codeparam_torbeam_data%nprofv
  noout        = codeparam_torbeam_data%noout
  nrela        = codeparam_torbeam_data%nrela
  nmaxh        = codeparam_torbeam_data%nmaxh
  nabsroutine  = codeparam_torbeam_data%nabsroutine
  nastra       = codeparam_torbeam_data%nastra
  nprofcalc    = codeparam_torbeam_data%nprofcalc
  ncdharm      = codeparam_torbeam_data%ncdharm
  npnts_extrap = codeparam_torbeam_data%npnts_extrap
  nfreq_extrap = codeparam_torbeam_data%nfreq_extrap
  nrel         = codeparam_torbeam_data%nrel

  !Float parameters
  xrtol   = codeparam_torbeam_data%xrtol
  xatol   = codeparam_torbeam_data%xatol
  xstep   = codeparam_torbeam_data%xstep
  rhostop = codeparam_torbeam_data%rhostop
  xzsrch  = codeparam_torbeam_data%xzsrch

  ! DETERMINE WHETHER PSI FLUX IS MAXIMUM (1) OR MINIMUM (-1) AT THE MAGNETIC AXIS
  if (psiedge > psiax) then
     sgnm = 1.
  else
     sgnm = -1.
  endif

  ! LOOP OVER BEAMS OF THE EC_LAUNCHERS IDS
  do ibeam=1,nbeam

     ! TO START FROM CLEAN
     intinbeam   = 0
     floatinbeam = 0.

     ! ONLY DEAL WITH ACTIVE BEAMS
     if(ec_launchers_ids%launcher(ibeam)%power_launched%data(1)>0) then

        ! IT LOOKS LIKE TORBEAM NEEDS PHI = 0, OTHERWISE IT DOES NOT TREAT THE BEAM PROPERLY
        ec_launchers_ids%launcher(ibeam)%launching_position%phi = 0.

        !intinbeam
        intinbeam(0)  =  2 ! tbr
        intinbeam(1)  = 2  ! tbr
        intinbeam(2)  = ec_launchers_ids%launcher(ibeam)%mode%data(1) ! (nmod)
        intinbeam(3)  = npow
        intinbeam(4)  = ncd
        intinbeam(5)  = 2 ! tbr
        intinbeam(6)  = ncdroutine
        intinbeam(7)  = nprofv
        intinbeam(8)  = noout
        intinbeam(9)  = nrela
        intinbeam(10) = nmaxh
        intinbeam(11) = nabsroutine
        intinbeam(12) = nastra
        intinbeam(13) = nprofcalc
        intinbeam(14) = ncdharm
        intinbeam(15) = npnts_extrap
        intinbeam(16) = nfreq_extrap
        intinbeam(17) = 0 ! always (no rhot.dat shall be involved when using IMAS)
        intinbeam(18) = nrel

        !floatinbeam(17:18): obsolete --> not filled)
        !floatinbeam(6:13):  analytic --> not filled)
        !floatinbeam(26:32): analytic --> not filled)
        floatinbeam(0)  = ec_launchers_ids%launcher(ibeam)%frequency%data(1)              ! (xf)
        beta  = -ec_launchers_ids%launcher(ibeam)%steering_angle_tor%data(1)
        alpha =  ec_launchers_ids%launcher(ibeam)%steering_angle_pol%data(1)
        xpoldeg = asin(cos(beta)*sin(alpha))*180./pi
        xtordeg = atan(tan(beta)/cos(alpha))*180./pi
        floatinbeam(1)  = xtordeg
        floatinbeam(2)  = xpoldeg
        floatinbeam(3)  = 1.e2*ec_launchers_ids%launcher(ibeam)%launching_position%r(1) & ! (xxb)
             *cos(ec_launchers_ids%launcher(ibeam)%launching_position%phi(1))
        floatinbeam(4)  = 1.e2*ec_launchers_ids%launcher(ibeam)%launching_position%r(1) & ! (xyb)
             *sin(ec_launchers_ids%launcher(ibeam)%launching_position%phi(1))
        floatinbeam(5)  = 1.e2*ec_launchers_ids%launcher(ibeam)%launching_position%z(1)   ! (xzb)

        floatinbeam(14) = xrtol  ! keep
        floatinbeam(15) = xatol  ! keep
        floatinbeam(16) = xstep  ! keep
        floatinbeam(19) = -1.e2_rkind/(ec_launchers_ids%launcher(ibeam)%beam%phase%curvature%data(1,1)+1.e-4_rkind) ! (xryyb)
        floatinbeam(20) = -1.e2_rkind/(ec_launchers_ids%launcher(ibeam)%beam%phase%curvature%data(2,1)+1.e-4_rkind) ! (xrzzb)
        floatinbeam(21) = ec_launchers_ids%launcher(ibeam)%beam%spot%size%data(1,1)*1.e2             ! (xwyyb)
        floatinbeam(22) = ec_launchers_ids%launcher(ibeam)%beam%spot%size%data(2,1)*1.e2             ! (xwzzb)
        floatinbeam(23) = ec_launchers_ids%launcher(ibeam)%power_launched%data(1)*1.e-6              ! (xpw0)
        floatinbeam(24) = equilibrium_ids%time_slice(1)%boundary%geometric_axis%r*1e2                ! (xrmaj)
        floatinbeam(25) = equilibrium_ids%time_slice(1)%boundary%minor_radius*1e2                    ! (xrmin)
        floatinbeam(33) = sgnm    ! (deduced from psi_ed-psi_ax)
        floatinbeam(34) = core_profiles_ids%profiles_1d(1)%zeff(1)                                   ! (xzeff)
        floatinbeam(35) = rhostop ! keep
        floatinbeam(36) = xzsrch  ! keep

        print*,'------------------------------------------------------------'
        print*,'Input power (MW), ibeam',ec_launchers_ids%launcher(ibeam)%power_launched%data*1.e-6,ibeam

        ! CALL TORBEAM
        call beam(intinbeam,floatinbeam,ni,nj,eqdata,npsi,npsi,prdata,rhoresult, &
             iend,t1data,t1tdata,kend,t2data,t2ndata,icnt,ibgout,nprofv,volprof, &
             power_array,npar_array,nper_array,slength_array,powout)

!!$        ! For debugging purpose
!!$        do iii=0,maxint-1
!!$           write(330,*) intinbeam(iii)
!!$        enddo
!!$        do iii=0,maxflt-1
!!$           write(331,*) floatinbeam(iii)
!!$        enddo
!!$        do iii=0,maxdiml-1
!!$           write(332,*) eqdata(iii)
!!$        enddo
!!$        do iii=0,maxlen-1
!!$           write(333,*) prdata(iii)
!!$        enddo

        ! --------------------------------------------------------------------------------------------------
        ! Result structures:
        ! --------------------------------------------------------------------------------------------------
        ! Beam propagation:
        ! - t1data  = (6 variables)
        !   * R - major-radius coordinate of the central ray                  (0:iend-1)
        !   * Z - vertical coordinate of the central ray                      (iend:2*iend-1)
        !   * R - major radius of the upper peripheral ray, i.e.interaction   (2*iend:3*iend-1)
        !         of beam width with the poloidal plane above the central ray
        !   * Z - vertical coordinate of the upper peripheral ray             (3*iend:4*iend-1)
        !   * R - major radius of the lower peripheral ray                    (4*iend:5*iend-1)
        !   * Z - vertical coordinate of the lower peripheral ray             (5*iend:6*iend-1).
        ! --------------------------------------------------------------------------------------------------
        ! - t1tdata = (4 variables)
        !   * X-coordinate of the central ray
        !   * Y-coordinate of the central ray (i.e. projection of the central ray onto a horizontal plane)
        !   * X-coordinate of left and right peripheral rays
        !   * Y-coordinate of left and right peripheral rays
        !   (intersection of the beam width and the horizontal plane running throuh the central ray)
        ! --------------------------------------------------------------------------------------------------
        ! Absorption and current drive
        ! --------------------------------------------------------------------------------------------------
        ! - t2data  = OBSOLETE
        ! - t2ndata = (3 variables)
        !   * radial coordinate (rho_p or rho_t as above, 0:npnt-1)
        !   * power density in MW/m3 (npnt:2*npnt-1)
        !   * driven current density in MA/m2 (2*npnt:3*npnt-1)
        ! --------------------------------------------------------------------------------------------------

        npointsout(ibeam)    = iend
        extrascal(ibeam,1)   = ec_launchers_ids%launcher(ibeam)%power_launched%data(1)*1.e-6
        extrascal(ibeam,2)   = 1.e6*powout
        extrascal(ibeam,3)   = 1.e6*rhoresult(12)
        extradata(ibeam,1,:) = ec_launchers_ids%launcher(ibeam)%power_launched%data(1)*1.e-6*(1.-power_array(:))
        extradata(ibeam,2,:) = -npar_array(:)*b0/abs(b0)
        extradata(ibeam,3,:) = nper_array(:)
        extradata(ibeam,4,:) = 1.e-2*slength_array(:)

        ! 1D PROFILES OF RHO, DP/DV, J (CHECKED OK, DIM = NPNT = 5000)
        do kp=0,2
           do lfd=0,npnt-1
              profout(ibeam,kp+1,lfd+1) = t2ndata(kp*npnt+lfd)
           enddo
        enddo

        ! TO PREVENT GOING BEYOND THE PRE-DEFINED TRAJOUT ARRAY
        if(iend.gt.ntraj) iend = ntraj

        ! TRAJECTORY OF CENTRAL RAY AND 4 PERIPHERAL "RAYS"
        do lfd=0,iend-1
           ! 1st ray
           trajout(ibeam, 1,lfd+1) = t1data(lfd)                    ! r
           trajout(ibeam, 2,lfd+1) = t1data(iend+lfd)               ! z
           !trajout(ibeam, 3,lfd+1) = acos(t1tdata(lfd)/t1data(lfd)) ! phi
           trajout(ibeam, 3,lfd+1) = asin(t1tdata(iend+lfd)/t1data(lfd)) ! phi
           ! 2nd ray
           trajout(ibeam, 4,lfd+1) = t1data(2*iend+lfd)
           trajout(ibeam, 5,lfd+1) = t1data(3*iend+lfd)
           !trajout(ibeam, 6,lfd+1) = acos(t1tdata(lfd)/t1data(2*iend+lfd))
           trajout(ibeam, 6,lfd+1) = asin(t1tdata(iend+lfd)/t1data(2*iend+lfd))
           ! 3rd ray
           trajout(ibeam, 7,lfd+1) = t1data(4*iend+lfd)
           trajout(ibeam, 8,lfd+1) = t1data(5*iend+lfd)
           !trajout(ibeam, 9,lfd+1) = acos(t1tdata(lfd)/t1data(4*iend+lfd))
           trajout(ibeam, 9,lfd+1) = asin(t1tdata(iend+lfd)/t1data(4*iend+lfd))
           ! 4th ray
           trajout(ibeam,10,lfd+1) = sqrt(t1tdata(2*iend+lfd)**2+t1tdata(3*iend+lfd)**2)
           trajout(ibeam,11,lfd+1) = t1data(iend+lfd)
           !trajout(ibeam,12,lfd+1) = acos(t1tdata(2*iend+lfd)/trajout(ibeam,10,lfd+1))
           trajout(ibeam,12,lfd+1) = asin(t1tdata(3*iend+lfd)/trajout(ibeam,10,lfd+1))
           ! 5th ray
           trajout(ibeam,13,lfd+1) = sqrt(t1tdata(4*iend+lfd)**2+t1tdata(5*iend+lfd)**2)
           trajout(ibeam,14,lfd+1) = t1data(iend+lfd)
           !trajout(ibeam,15,lfd+1) = acos(t1tdata(4*iend+lfd)/trajout(ibeam,13,lfd+1))
           trajout(ibeam,15,lfd+1) = asin(t1tdata(5*iend+lfd)/trajout(ibeam,13,lfd+1))
        enddo

     endif ! TEST BEAM_POWER > 0

  enddo ! LOOP OVER BEAMS OF EC_LAUNCHERS IDS

  ! FOR OUTPUT PROFILES SUMMED UP OVER ALL BEAMS
  proftot = 0.
  do lfd = 1,npnt
     proftot(1,lfd) = profout(1,1,lfd)
     do ibeam = 1,nbeam
        proftot(2,lfd) = proftot(2,lfd)+profout(ibeam,2,lfd)
        proftot(3,lfd) = proftot(3,lfd)+profout(ibeam,3,lfd)
     enddo
  enddo

  ! ----------------------------
  ! SAVE RESULTS INTO WAVES IDS
  ! ----------------------------

  ! GLOBAL QUANTITIES
  waves_ids%ids_properties%homogeneous_time = 1
  allocate(waves_ids%time(1))
  waves_ids%time = equilibrium_ids%time
  allocate(waves_ids%code%name(1))
  waves_ids%code%name = 'TORBEAM'

  ! LOOP OVER BEAMS (LAUNCHERS)
  !if(nbeam.gt.10) nbeam = 10 ! MSR waiting for IMAS-3271
  allocate(waves_ids%coherent_wave(nbeam))
  do ibeam=1,nbeam
     allocate(waves_ids%coherent_wave(ibeam)%global_quantities(1))
     allocate(waves_ids%coherent_wave(ibeam)%identifier%antenna_name(1))
     waves_ids%coherent_wave(ibeam)%identifier%antenna_name                      = ec_launchers_ids%launcher(ibeam)%name
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%frequency               = ec_launchers_ids%launcher(ibeam)%frequency%data(1)
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%electrons%power_thermal = extrascal(ibeam,2)
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%power                   = extrascal(ibeam,2)
     waves_ids%coherent_wave(ibeam)%global_quantities(1)%current_tor             = extrascal(ibeam,3)
     waves_ids%coherent_wave(ibeam)%identifier%type%description                  = 'TORBEAM'
     ! if statement needed to filter empty sources for hcd2core_sources
     if(ec_launchers_ids%launcher(ibeam)%power_launched%data(1)>0) then
        waves_ids%coherent_wave(ibeam)%identifier%type%name                         = 'EC'
        waves_ids%coherent_wave(ibeam)%identifier%type%index                        = 1
        waves_ids%coherent_wave(ibeam)%wave_solver_type%index                       = 1 ! BEAM/RAY TRACING
     endif
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor_norm(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%power_density(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%electrons%power_density_thermal(npnt))
     allocate(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%current_parallel_density(npnt))
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi                        = profout(ibeam,1,1:npnt)**2*(psiedge-psiax)+psiax
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%power_density                   = 1.e6*profout(ibeam,2,1:npnt)
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%electrons%power_density_thermal = 1.e6*profout(ibeam,2,1:npnt)
     waves_ids%coherent_wave(ibeam)%profiles_1d(1)%current_parallel_density        = -1.e6*profout(ibeam,3,1:npnt)*iplasma/abs(iplasma)

     ! rho_tor and rho_tor_norm to be filled for hcd2core_sources not to complain
     call interpos(-equilibrium_ids%time_slice(1)%profiles_1d%psi, &
          equilibrium_ids%time_slice(1)%profiles_1d%rho_tor, &
          size(equilibrium_ids%time_slice(1)%profiles_1d%psi), &
          size(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi), &
          xout=-waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi, &
          yout=waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor)
     call interpos(-equilibrium_ids%time_slice(1)%profiles_1d%psi, &
          equilibrium_ids%time_slice(1)%profiles_1d%rho_tor_norm, &
          size(equilibrium_ids%time_slice(1)%profiles_1d%psi), &
          size(waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi), &
          xout=-waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%psi, &
          yout=waves_ids%coherent_wave(ibeam)%profiles_1d(1)%grid%rho_tor_norm)

     ! LOOP OVER RAYS
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(5))
     iend = min(npointsout(ibeam),ntraj)
     if(ec_launchers_ids%launcher(ibeam)%power_launched%data(1)>0) then
        do iray = 1,5 ! 5 rays (to not mix with input beams)
           allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%r(iend))
           allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%z(iend))
           allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%phi(iend))
           waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%r   = 1.e-2*trajout(ibeam,1+3*(iray-1),1:iend)
           waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%z   = 1.e-2*trajout(ibeam,2+3*(iray-1),1:iend)
           waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(iray)%position%phi = trajout(ibeam,3+3*(iray-1),1:iend)
        enddo
     endif

     ! SOME QUANTITIES ARE KNOWN ONLY ON THE CENTRAL RAY
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%electrons%power(iend))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%parallel(iend))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%perpendicular(iend))
     allocate(waves_ids%coherent_wave(ibeam)%beam_tracing%beam(1)%length(iend))
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_initial                 = extrascal(ibeam,1)*1.e6
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%electrons%power               = extradata(ibeam,1,1:iend)*1.e6
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%parallel      = extradata(ibeam,2,1:iend)
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%power_flow_norm%perpendicular = extradata(ibeam,3,1:iend)
     waves_ids%coherent_wave(ibeam)%beam_tracing(1)%beam(1)%length                        = extradata(ibeam,4,1:iend)

  enddo ! LOOP OVER BEAMS (LAUNCHERS)

  waves_ids%code%output_flag = 0 ! NO ERROR

222 format(6(1x,1pe13.5))
225 format(3(1x,1pe13.5))
999 format(1p,3e13.5,i5,e13.5)

  deallocate(psi)
  deallocate(Rarr,Zarr)
  deallocate(br)
  deallocate(bt)
  deallocate(bz)
  deallocate(psi2d)
  deallocate(te,ne,Zeff)
  deallocate(profout)
  deallocate(trajout)
  deallocate(npointsout)
  deallocate(extrascal)
  deallocate(extradata)

end subroutine torbeam

end module torbeam_iwrap
