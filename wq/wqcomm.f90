! wqcomm.f90

MODULE wqcomm_parm

  USE bpsd_kinds
  USE bpsd_constants
  USE commpi
  
  IMPLICIT NONE
  PUBLIC

  REAL(rkind):: FREQ
  REAL(rkind):: dtfactor
  REAL(rkind):: dxfactor
  REAL(rkind):: dyfactor
  REAL(rkind):: nufactor
  REAL(rkind):: B0
  REAL(rkind):: RR
  REAL(rkind):: RA
  REAL(rkind):: q0
  REAL(rkind):: qa
  REAL(rkind):: n0
  REAL(rkind):: TMN
  INTEGER:: ntmax
  INTEGER:: nxmax
  INTEGER:: nymax
  INTEGER:: INMODE

  INTEGER:: model_pulse
  INTEGER:: model_dielectric
  INTEGER:: model_plot
  REAL(rkind):: pulse_cycle
  REAL(rkind):: dielectric_2
  REAL(rkind):: dielectric_3
  REAL(rkind):: freq_resonance
  REAL(rkind):: freq_collision
  INTEGER:: ntplot_interval
  INTEGER:: ntplot_max

END MODULE wqcomm_parm

MODULE wqcomm

  USE wqcomm_parm
  IMPLICIT NONE

  REAL(rkind),ALLOCATABLE :: &
       ne(:,:),OCE(:,:),OPE(:,:),OUH(:,:),pabs(:,:),OR(:,:),OL(:,:)
  COMPLEX(rkind),ALLOCATABLE :: &
       EX(:,:),EY(:,:),EZ(:,:),A(:,:,:,:),AA(:,:),B(:,:), &
       CD(:,:,:,:),CDplus(:,:,:,:),CDminus(:,:,:,:)
  REAL(rkind):: &
       omega,period,wavelength,dt,dx,dy,nu,omegaplus,omegaminus,domega
  INTEGER:: &
       nttot,ntplot
  REAL(rkind):: &
       ttot
  REAL(rkind),ALLOCATABLE:: &
       tplot(:)
  COMPLEX(rkind),ALLOCATABLE :: &
       EX_save(:,:,:),EY_save(:,:,:),EZ_save(:,:,:)
CONTAINS

  SUBROUTINE wq_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: nxmax_save,nymax_save
    INTEGER,SAVE:: INIT=0

    IF(INIT.EQ.0) THEN
       INIT=1
    ELSE
       IF(nxmax.EQ.nxmax_save.AND. &
          nymax.EQ.nymax_save) RETURN
       CALL wq_deallocate
    END IF

    ALLOCATE(EX(nxmax,nymax),EY(nxmax,nymax),EZ(nxmax,nymax))
    ALLOCATE(A(3,3,nxmax,nymax),AA(3,3),B(3,3))
    ALLOCATE(CD(3,3,nxmax,nymax))
    ALLOCATE(CDplus(3,3,nxmax,nymax))
    ALLOCATE(CDminus(3,3,nxmax,nymax))
    ALLOCATE(ne(nxmax,nymax),OCE(nxmax,nymax))
    ALLOCATE(OPE(nxmax,nymax),OUH(nxmax,nymax))
    ALLOCATE(pabs(nxmax,nymax),OR(nxmax,nymax),OL(nxmax,nymax))
    IF(ntplot_max.GT.0) THEN
       ALLOCATE(tplot(ntplot_max))
       ALLOCATE(EX_save(nxmax,nymax,ntplot_max))
       ALLOCATE(EY_save(nxmax,nymax,ntplot_max))
       ALLOCATE(EZ_save(nxmax,nymax,ntplot_max))
    END IF

    nxmax_save=nxmax
    nymax_save=nymax

  END SUBROUTINE wq_allocate

  SUBROUTINE wq_deallocate
    IMPLICIT NONE
   
    DEALLOCATE(EX,EY,EZ)
    DEALLOCATE(A,AA,B)
    DEALLOCATE(CD,CDplus,CDminus)
    DEALLOCATE(ne,OCE,OPE,OUH,pabs,OR,OL)
    IF(ALLOCATED(EX_save)) DEALLOCATE(EX_save,EY_save,EZ_save,tplot)

  END SUBROUTINE wq_deallocate
END MODULE wqcomm

