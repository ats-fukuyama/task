!     $Id$
!
!             ########## /TASK/WF ##########
!
!             WAVE PROPAGATION AND ABSORPTION
!
!                  HOT/COLD PLASMA
!                  3-DIMENSIONAL INHOMOGENEITY
!                  FIRST ORDER FINITE ELEMENT METHOD
!                  SCALAR AND VECTOR POTENTIAL
!
!                  CODED BY A. FUKUYAMA
!
!                  DEPARTMENT OF NUCLEAR ENGINEERING
!                  KYOTO UNIVERSITY
!                  KYOTO 606-8501, JAPAN
!
!                  V1.0   : 1996 JULY (CARTESIAN COORDINATES)
!                  V1.01  : 1996 SEPT (CYLINDRICAL COORDINATES)
!                  V1.1   : 1996 SEPT (TRANSLATIONAL THREE-DIMENSIONAL)
!                  V2.0   : 2003 JAN  (FULL THREE-DIMENSIONAL)
!                  V3.0   : 2003 JUN  (SIDE ELEMENTS)
!
!     ************************************************************
!                  V3.0T  : 2013 FEB  (GAMMA10 tamdem mirror)
!                  Modified by T. YOKOYAMA
!                  Univ. of Tsukuba (GAMMA10 team)
!                     25/02/2013-
!
!    This work is partially supported by the bidirectional collaborative 
!    research program of National Institute for Fusion Science, Japan
!    (NIFS11KUGM050) (NIFS11KUGM055)

program wfmain
  
  use wfcomm
  use libmtx
  implicit none
  
  ! --- initialize ---
  call pl_allocate_ns
  call mtx_initialize

  if(nrank.eq.0) then
     write(6,*) '## TASKX/WFX  V3.02  2010/10/09 ###'
     call gsopen
  end if

  call setaif
  call wfinit

  if (nrank.eq.0) call wfparf
  call wfparm_broadcast

  ! --- menu ---

  call wfmenu
  
  ! --- finalize ---
  
  if (nrank.eq.0) call gsclos
  if(NFOPEN.ne.0) close(26)
  call mtxc_cleanup
  call wfelm_deallocate
  call wfsid_deallocate
  call wfsrt_deallocate
  call wfmed_deallocate
  call wffld_deallocate
  call wfpwr_deallocate
  if(nrank.eq.0) call wfwin_deallocate
  call mtx_finalize
  stop
end program wfmain
