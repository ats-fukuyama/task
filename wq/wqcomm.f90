! wqcomm.f90

MODULE wqcomm_parm

  USE bpsd_kinds
  USE bpsd_constants
  USE commpi
  
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER:: medium_m=10   ! maximum number of medium
  INTEGER,PARAMETER:: idebug_m=99   ! size of array idebuga
  INTEGER,PARAMETER:: ngt_m=10001   ! maximum number of global value storage
  INTEGER,PARAMETER:: ngr_m=101     ! maximum number of profile storage

  INTEGER:: model_geometry
  REAL(rkind):: xnmin
  REAL(rkind):: xnmax
  REAL(rkind):: ynmin
  REAL(rkind):: ynmax

  REAL(rkind):: B0
  REAL(rkind):: RR
  REAL(rkind):: RA
  REAL(rkind):: q0
  REAL(rkind):: qa

  REAL(rkind):: freq
  REAL(rkind):: rkz
  INTEGER:: nph

  INTEGER:: model_source
  REAL(rkind):: source_position_xn
  REAL(rkind):: source_position_yn
  REAL(rkind):: source_width
  REAL(rkind):: source_thickness
  REAL(rkind):: source_angle

  INTEGER:: model_pulse
  INTEGER:: model_ramp
  REAL(rkind):: pulse_length
  REAL(rkind):: ramp_length

  INTEGER:: medium_max
  INTEGER:: id_medium(medium_m)
  REAL(rkind):: xnmin_medium(medium_m)
  REAL(rkind):: xnmax_medium(medium_m)
  REAL(rkind):: ynmin_medium(medium_m)
  REAL(rkind):: ynmax_medium(medium_m)
  REAL(rkind):: dielectric_medium(medium_m)
  REAL(rkind):: res_freq_medium(medium_m)
  REAL(rkind):: res_coll_medium(medium_m)
  REAL(rkind):: density_medium(medium_m)
  REAL(rkind):: collision_medium(medium_m)

  INTEGER:: model_solver
  INTEGER:: model_plot

  REAL(rkind):: fimplicit
  INTEGER:: ntype_mat
  REAL(rkind):: eps_mat

  REAL(rkind):: dtfactor
  REAL(rkind):: dxfactor
  REAL(rkind):: dyfactor

  INTEGER:: ntmax
  INTEGER:: ntstep
  INTEGER:: ngtstep
  INTEGER:: ngrstep

  INTEGER:: idebuga(idebug_m)

END MODULE wqcomm_parm

MODULE wqcomm

  USE wqcomm_parm
  IMPLICIT NONE

  INTEGER:: nxmax
  INTEGER:: nymax

  INTEGER:: ngt_max,ngr_max
  REAL(rkind):: omega,period,wave_length,wave_number,dt
  REAL(rkind):: omegaplus,omegaminus,domega
  INTEGER:: icount_mat,icount_mat_max
  REAL(rkind):: t_tot
  INTEGER:: nt_tot

  INTEGER,DIMENSION(:,:),ALLOCATABLE :: &
       medium_nx_ny(:,:)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &
       xg(:),yg(:),xgn(:),ygn(:)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE :: &
       ne,OCE,OPE,OUH,pabs,OR,OL
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE :: &
       EX,EY,EZ
  COMPLEX(rkind),DIMENSION(:,:,:,:),ALLOCATABLE :: &
       A,Ainv, CD,CDplus,CDminus
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE :: &
       EX_save,EY_save,EZ_save
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE :: &
       pabs_save
  REAL(rkind),DIMENSION(:),ALLOCATABLE :: &
       t_ngt,t_ngr,pabs_tot_ngt,EX_max_ngt,EY_max_ngt,EZ_max_ngt
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

    ALLOCATE(medium_nx_ny(nxmax,nymax))
    ALLOCATE(xg(nxmax),yg(nymax))
    ALLOCATE(xgn(nxmax),ygn(nymax))
    ALLOCATE(EX(nxmax,nymax),EY(nxmax,nymax),EZ(nxmax,nymax))
    ALLOCATE(A(3,3,nxmax,nymax),Ainv(3,3,nxmax,nymax))
    ALLOCATE(CD(3,3,nxmax,nymax))
    ALLOCATE(CDplus(3,3,nxmax,nymax))
    ALLOCATE(CDminus(3,3,nxmax,nymax))
    ALLOCATE(ne(nxmax,nymax),OCE(nxmax,nymax))
    ALLOCATE(OPE(nxmax,nymax),OUH(nxmax,nymax))
    ALLOCATE(pabs(nxmax,nymax),OR(nxmax,nymax),OL(nxmax,nymax))
    ALLOCATE(t_ngt(ngt_m),pabs_tot_ngt(ngt_m))
    ALLOCATE(EX_max_ngt(ngt_m),EY_max_ngt(ngt_m),EZ_max_ngt(ngt_m))
    ALLOCATE(t_ngr(ngr_m))
    ALLOCATE(EX_save(nxmax,nymax,ngr_m))
    ALLOCATE(EY_save(nxmax,nymax,ngr_m))
    ALLOCATE(EZ_save(nxmax,nymax,ngr_m))
    ALLOCATE(pabs_save(nxmax,nymax,ngr_m))

    nxmax_save=nxmax
    nymax_save=nymax

  END SUBROUTINE wq_allocate

  SUBROUTINE wq_deallocate
    IMPLICIT NONE
   
    DEALLOCATE(medium_nx_ny)
    DEALLOCATE(xg,yg)
    DEALLOCATE(xgn,ygn)
    DEALLOCATE(EX,EY,EZ)
    DEALLOCATE(A,Ainv)
    DEALLOCATE(CD,CDplus,CDminus)
    DEALLOCATE(ne,OCE,OPE,OUH,pabs,OR,OL)
    DEALLOCATE(t_ngt,EX_max_ngt,EY_max_ngt,EZ_max_ngt,pabs_tot_ngt)
    DEALLOCATE(t_ngr,EX_save,EY_save,EZ_save,pabs_save)


  END SUBROUTINE wq_deallocate
END MODULE wqcomm

