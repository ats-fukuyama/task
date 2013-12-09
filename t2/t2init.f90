!C
!C
!C
MODULE T2INIT
  
  USE T2CNST, ONLY: i0ikind,i0rkind
  
  IMPLICIT NONE

  PUBLIC T2_INIT !C INITIALIZE INPUT PARAMETERS
  
  PRIVATE
CONTAINS

  !C
  !C
  !C

  SUBROUTINE T2_INIT
    
    USE T2COMM, ONLY: &
         c10rname, i0dbg, i0fnum, i0mfcs, i0wstp,&
         i0dmax0,i0amax0,&
         i0tmax, d0tstp, d0tmax,&
         i0spcs, i0nmax0, i0lmax, i1mlvl,&
         i0pdiv_number, i1rdn2, d1rec,&
         i0pmax,d0eps,d0rmjr,d0rmnr,&
         i0m0,i0n0,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,&
         d1pa,d1pz,&
         d0qc,d0qs,d0rw

    USE T2CNST, ONLY: i0lmaxm,i0spcsm,d0aee,d0ame,d0amp

    NAMELIST /T2/ &
         c10rname, i0dbg, i0fnum, i0mfcs, i0wstp,&
         i0dmax0,i0amax0,&
         i0tmax, d0tstp, d0tmax,&
         i0spcs, i0nmax0, i0lmax, i1mlvl,&
         i0pdiv_number, i1rdn2, d1rec,&
         i0pmax,d0eps,d0rmjr,d0rmnr,&
         i0m0,i0n0,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,&
         d1pa,d1pz,&
         d0qc,d0qs,d0rw

    c10rname = 'TEST'
    i0dbg    =  0
    i0fnum   = 10
    i0mfcs   =  1

    i0wstp   =  1  ! output timing

    i0dmax0  =  2  ! mesh dim
    i0amax0  = 32  ! gauss kyuuseki number of sample points
    i0tmax   = 10  ! 
    d0tstp   = 1.D-12  !
    d0tmax   = 1.D-12  !
    d0eps    = 1.D-4
    i0pmax   =  49     ! iteration 

    i0spcs   =  2 
    i0nmax0  =  4      ! number of nodes in a elemnt
    i0lmax   =  6      ! 
    i0pdiv_number = 6  ! 

    i1mlvl(1)=1        ! 6 x 2^0
    i1mlvl(2)=2        ! 6 x 2^1
    i1mlvl(3)=3
    i1mlvl(4)=4
    i1mlvl(5)=5        ! 6 x 2^5
    i1mlvl(6)=5        ! 6 x 2^5
    i1mlvl(  0)          = 0
    i1mlvl(  7:i0lmaxm+1) = 0

    i1rdn2(1)=  10     ! number of radial nodes in a level
    i1rdn2(2)=  10     
    i1rdn2(3)=  10
    i1rdn2(4)=  10
    i1rdn2(5)=  10
    i1rdn2(6)=  10
    i1rdn2(-1:0)     = 0
    i1rdn2(7:i0lmaxm) = 0

    d1rec(0) = 0.000D0  ! least radial point in a level
    d1rec(1) = 0.050D0
    d1rec(2) = 0.100D0
    d1rec(3) = 0.400D0
    d1rec(4) = 0.700D0
    d1rec(5) = 1.000D0
    d1rec(6) = 1.100D0
    d1rec(7:i0lmaxm) = 0.D0

    d0rmjr   =  3.0D0
    d0rmnr   =  1.0D0
    d0rw     =  1.1D0
    d0bc     =  3.0D0
    d0qc     =  1.0D0
    d0qs     =  3.0D0

!   PLASMA PARAMETER
    i0m0     =  1       ! pressure profile parameter
    i0n0     =  3       ! pressure profile parameter

!    Electron 
    d1pa(1) =  d0ame/d0amp
    d1pz(1) = -1.D0
    d1nc(1) = 1.0D0
    d1ns(1) = 2.0D-1
    d1nw(1) = 5.0D-2
    d1tc(1) = 5.0D0
    d1ts(1) = 1.0D0
    d1tw(1) = 1.0D-1

!   Deuterium
    d1pa(2) = 2.D0
    d1pz(2) = 1.D0
    d1nc(2) = 1.D0
    d1ns(2) = 2.D-1
    d1nw(2) = 5.D-2
    d1tc(2) = 5.D0
    d1ts(2) = 1.D0
    d1tw(2) = 1.D-1

!   Other plasma speces
    d1nc(3:i0spcsm) = 0.D0
    d1ns(3:i0spcsm) = 0.D0 
    d1nw(3:i0spcsm) = 0.D0
    d1tc(3:i0spcsm) = 0.D0
    d1ts(3:i0spcsm) = 0.D0
    d1tw(3:i0spcsm) = 0.D0
    d1pa(3:i0spcsm) = 0.D0
    d1pz(3:i0spcsm) = 0.D0

    RETURN
  END SUBROUTINE T2_INIT
END MODULE T2INIT
    
