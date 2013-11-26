!C 
!C NAMELIST OF CONTROL PARAMETERS IN TAST/T2
!C

&T2
c10rname     ='TEST'
!C 
i0dbg        = 4
i0fnum       = 10
i0dmax0      = 2
i0amax0      = 32
i0spcs       = 2
i0wstp       = 1
i0mfcs       = 1
i0pmax = 49
i0tmax = 10
d0tmax = 1.D-30
d0eps  = 1.D-4
d0tstp = 1.D-30

!C FOR GRID GENERATION

i0nmax0      = 4
i0lmax       = 6
i0pdiv_number= 6

i1mesh_level(1)=1
i1mesh_level(2)=2
i1mesh_level(3)=3
i1mesh_level(4)=4
i1mesh_level(5)=5
i1mesh_level(6)=5


i1rdiv_number(1)=  10
i1rdiv_number(2)=  10
i1rdiv_number(3)=  10
i1rdiv_number(4)=  10
i1rdiv_number(5)=  10
i1rdiv_number(6)=  10

d1rec_tmp(0) = 0.000D0
d1rec_tmp(1) = 0.050D0
d1rec_tmp(2) = 0.100D0
d1rec_tmp(3) = 0.400D0
d1rec_tmp(4) = 0.700D0
d1rec_tmp(5) = 1.000D0
d1rec_tmp(6) = 1.100D0

d0qc    = 1.0D0
d0qs    = 3.0D0
d0bc    = 3.0D0
d0rmjr  = 3.0D0
!d0rmjr = 1.0D6
d0rmnr  = 1.0D0
d0rw    = 1.1D0

!C PLASMA PARAMETER
i0m0 = 1
i0n0 = 3
!C ELECTRON 
d1nc_tmp(1) = 1.0D0
d1ns_tmp(1) = 2.0D-1
d1nw_tmp(1) = 5.0D-2
d1tc_tmp(1) = 5.0D0
d1ts_tmp(1) = 1.0D0
d1tw_tmp(1) = 1.0D-1
d1pa_tmp(1) = 5.446169971D-4 
d1pz_tmp(1) = -1.D0

!C DEUTERIUM
d1nc_tmp(2) = 1.D0
d1ns_tmp(2) = 2.D-1
d1nw_tmp(2) = 5.D-2
d1tc_tmp(2) = 5.D0
d1ts_tmp(2) = 1.D0
d1tw_tmp(2) = 1.D-1
d1pa_tmp(2) = 2.D0
d1pz_tmp(2) = 1.D0
&END