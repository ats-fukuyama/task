module fowob
  private
  integer :: nobt_in_max
  real(8),allocatable :: penergym(:,:),penergyg(:,:)
  integer,allocatable,dimension(:,:,:,:) :: nobt_in


  public :: fow_ob_run

contains

  subroutine fow_ob_run

    use fowcomm,only : nzemax_fow=>nzemax, npsmax_fow=>npsmax,orbit_p, orbit_ze, orbit_ps,construct_orbit
    use fpcomm, only : npmax_fp=>npmax, nsamax_fp=>nsamax

    use plinit
    use obcomm
    use obinit
    use obparm
    use obprep
    use obcalc

    implicit none
    integer :: ierr,nobt,np,nze,nps,nsa,mode(3)

    ierr = 0

    call fow_ob_prep

    CALL pl_init
    CALL EQINIT
    CALL ob_init
    CALL ob_parm(1,'../fp.ota/fpparm',ierr)
    
    if ( nobt_in_max>=nobt_m ) then
      write(*,*)"ERROR at fow_ob_run : nobt_max must be less than nobt_m"
      write(*,*)"nobt_max = ",nobt_in_max
      write(*,*)"nobt_m   = ",nobt_m
      STOP
    end if

    nobt_max = nobt_in_max

    CALL ob_prep(ierr)
    CALL ob_allocate

    mode=[0,0,0]
    call fow_initial_value(penergy_in, pcangle_in, zeta_in, psipn_in, theta_in, mode)
    ! do nobt = 1, nobt_in_max
    !   write(*,*)penergy_in(nobt)
    ! end do
    WRITE(*,*)"nobt_m",nobt_m,tmax

    CALL ob_calc(ierr)

    CALL ob_deallocate
    deallocate(penergym, penergyg)

  end subroutine fow_ob_run

  subroutine fow_ob_prep

    use fowcomm
    use fpcomm
    use foworbitclassify

    implicit none

    integer :: np,nze,nps,nsa,i,mode(3)
    real(8) :: PVM,PVG

    allocate(penergym(npmax,nsamax))
    allocate(penergyg(npmax+1,nsamax))
    allocate(nobt_in(npmax,nzemax,npsmax,nsamax))

    do nsa=1,nsamax
      do np=1,npmax+1
        PVG = SQRT(1.D0+THETA0(nsa)*PG(np,nsa)**2)
        penergyg(np,nsa) = (PVG-1.D0)*PTFP0(nsa)**2/(THETA0(nsa)*AMFP(nsa))/(1.d3*aee)
        if( np/=npmax+1 ) then
          PVM = SQRT(1.D0+THETA0(nsa)*PM(np,nsa)**2)
          penergym(np,nsa) = (PVM-1.D0)*PTFP0(nsa)**2/(THETA0(nsa)*AMFP(nsa))/(1.d3*aee)
        end if
       end do
     end do

     i = 0
     mode = [0,0,0]
     do nsa = 1, nsamax
       do nps = 1, npsmax
         do nze = 1, nzemax
           do np = 1, npmax
            if ( .not.forbitten(np,nze,nps,nsa,mode) ) then
              i = i+1
              nobt_in(np,nze,nps,nsa) = i
            end if
           end do
         end do
       end do
     end do

     nobt_in_max = i
     i = 0

  end subroutine fow_ob_prep

  subroutine fow_initial_value(penergy_in, pcangle_in, zeta_in, psipn_in, theta_in, mode)

    use fowcomm
    use fpcomm
    use foworbitclassify  

    implicit none
    integer,intent(in) :: mode(3)
    real(rkind),intent(inout) :: penergy_in(:), pcangle_in(:), zeta_in(:), psipn_in(:), theta_in(:)
    integer :: i,np,nze,nps,nsa

    do nsa = 1, nsamax
      do nps = 1, npsmax+mode(3)
        do nze = 1, nzemax+mode(2)
          do np = 1, npmax+mode(1)
            if ( .not.forbitten(np,nze,nps,nsa,mode) ) then

              i = nobt_in(np,nze,nps,nsa)

              select case(mode(1))
              case(0)
                penergy_in(i) = penergym(np,nsa)
              case(1)
                penergy_in(i) = penergyg(np,nsa)
              end select

              select case(mode(2))
              case(0)
                pcangle_in(i) = zeta(nze)
              case(1)
                pcangle_in(i) = zetag(nze)
              end select

              select case(mode(3))
              case(0)
                psipn_in(i) = psim(nps)/psi0
              case(1)
                psipn_in(i) = psimg(nps)/psi0
              end select

              if ( pcangle_in(i)>0.d0 ) then
                theta_in(i) = 0.d0
              else
                theta_in(i) = pi
              end if

              zeta_in(i) = 0.d0

            end if
          end do
        end do
      end do
    end do
  
  end subroutine fow_initial_value

end module fowob
