!     $Id: wfmain.f90,v 1.8 2012/03/05 06:29:02 maruyama Exp $
!
!     ########## /TASK/WF2 ##########
!
!     WAVE PROPAGATION AND ABSORPTION
!
!     FEATURES:
!          FOR AXIAL SYMMETRY PLASMA
!          2-DIMENSIONAL INHOMOGENEITY
!          COLD PLASMA DIELECTRIC TENSOR
!          FIRST ORDER FINITE ELEMENT METHOD
!          SCALAR AND VECTOR BASIS FUNCTION
!
!     CODED BY Y.MARUYAMA
!
!     DEPARTMENT OF NUCLEAR ENGINEERING
!     KYOTO UNIVERSITY
!     KYOTO 606-8501, JAPAN
!
!     V1.0   : 2011 DEC
!
!     ***************************************

program wfmain
  
  USE wfcomm
  USE plinit
  USE wfinit
  USE wfparm
  USE libmtx
  implicit none
  INTEGER:: ierr
  
  ! --- initialize ---
  call mtx_initialize
  call allocate_init
  iddiv=0

  if(nrank.eq.0) then
     write(6,*) '## TASK/WF2  V1.0  2011/12/30 ###'
     call gsopen
  end if

  call setaif
  call setaie
  CALL pl_init
  call wf_init

  if (nrank.eq.0) call wf_parm(1,'wfparm',IERR)
  call wfparm_broadcast

  ! --- menu ---
  call wfmenu

  ! --- finalize ---

  if(nrank.eq.0) call gsclos
  if(NFOPEN.ne.0) close(26)
  if(elminit.ne.0) call wfelm_deallocate
  if(sidinit.ne.0) call wfsid_deallocate
  if(srtinit.ne.0) call wfsrt_deallocate
!  if(medinit.ne.0) call wfmed_deallocate
  if(fldinit.ne.0) call wffld_deallocate
  if((nrank.eq.0).and.(wininit.ne.0)) call wfwin_deallocate

  call mtx_finalize
  stop
end program wfmain

!------------------------------------------

subroutine allocate_init
  use wfcomm
  implicit none

!  divinit=0
  elminit=0
  sidinit=0
  srtinit=0
!  medinit=0
!  medinit=0
  srfinit=0
  slvinit=0
!  fldinit=0
!  nasinit=0
  wininit=0

  return
end subroutine allocate_init
