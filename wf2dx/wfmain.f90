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
  USE wfmenu
  USE wfsub
  USE libmtx
  implicit none
  INTEGER:: ierr
  
  ! --- initialize ---

  call mtx_initialize

  if(nrank.eq.0) then
     write(6,*) '## TASK/WF2  V1.0  2011/12/30 ###'
     call gsopen
  end if

  call wf_set_aif
  call wf_set_aie
  CALL pl_init
  CALL eqinit
  call wf_init

  if (nrank.eq.0) call wf_parm(1,'wfparm',IERR)
  call wfparm_broadcast

  ! --- menu ---

  call wf_menu

  ! --- finalize ---

  if(nrank.eq.0) call gsclos
  if(NFOPEN.ne.0) close(26)

  call mtx_finalize
  stop
end program wfmain

