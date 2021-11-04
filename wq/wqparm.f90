! wqparm.f90

MODULE wqparm

  PRIVATE
  PUBLIC wq_parm,wq_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE wq_parm(mode,kin,ierr)
    
!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN):: kin
    INTEGER,INTENT(OUT):: ierr
    
1   CALL TASK_PARM(mode,'WQ',kin,wq_nlin,wq_plst,ierr)
    IF(ierr.NE.0) RETURN

    CALl wq_check(ierr)
    IF(mode.EQ.0.AND.ierr.NE.0) GOTO 1
    IF(ierr.NE.0) ierr=ierr+100

    RETURN
  END SUBROUTINE wq_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wq_nlin(nid,ist,ierr)

    USE wqcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nid
    INTEGER,INTENT(OUT):: ist,ierr

    NAMELIST /WQ/ &
         model_geometry,xnmin,xnmax,ynmin,ynmax, &
         B0,RR,RA,q0,qa, &
         freq,rkz,nph, &
         model_source,source_position_xn,source_position_yn, &
         source_width,source_thickness,source_angle, &
         model_pulse,pulse_length,model_ramp,ramp_length, &
         medium_max,id_medium,xnmin_medium,xnmax_medium,ynmin_medium, &
         ynmax_medium,dielectric_medium,res_freq_medium,res_coll_medium, &
         density_medium,collision_medium, &
         model_solver,model_plot, &
         fimplicit,ntype_mat,eps_mat, &
         dtfactor,dxfactor,dyfactor, &
         ntmax,ntstep,ngtstep,ngrstep,idebuga, &
         model_geometry,xnmin,xnmax,ynmin,ynmax, &
         B0,RR,RA,q0,qa, &
         freq,rkz,nph, &
         model_source,source_position_xn,source_position_yn, &
         source_width,source_thickness,source_angle, &
         model_pulse,pulse_length,model_ramp,ramp_length, &
         medium_max,id_medium,xnmin_medium,xnmax_medium,ynmin_medium, &
         ynmax_medium,dielectric_medium,res_freq_medium,res_coll_medium, &
         density_medium,collision_medium, &
         model_solver,model_plot, &
         fimplicit,ntype_mat,eps_mat, &
         dtfactor,dxfactor,dyfactor, &
         ntmax,ntstep,ngtstep,ngrstep,idebuga

    ierr=0

    READ(nid,WQ,IOSTAT=ist,ERR=9800,END=9900)
    IF(ist.NE.0) THEN
       WRITE(6,'(A,I5)') 'XX wq_nlin: READ ERROR: IOSTAT=',ist
       ierr=ist
       RETURN
    END IF
    RETURN

9800 IERR=8
    WRITE(6,'(A)') 'XX wq_nlin: READ ERROR'
    RETURN
9900 IERR=9
    WRITE(6,'(A)') 'XX wq_nlin: READ END OF FILE'
    RETURN
  END SUBROUTINE wq_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wq_plst

    WRITE(6,'(A)') &
         '# &WQ: model_geometry,xnmin,xnmax,ynmin,ynmax,', &
         '       B0,RR,RA,q0,qa,', &
         '       freq,rkz,nph,', &
         '       model_source,source_position_xn,source_position_yn,', &
         '       source_width,source_thickness,source_angle,', &
         '       model_pulse,pulse_length,model_ramp,ramp_length,', &
         '       medium_max,id_medium,', &
         '       xnmin_medium,xnmax_medium,ynmin_medium,ynmax_medium,', &
         '       dielectric_medium,res_freq_medium,res_coll_medium,', &
         '       density_medium,collision_medium,', &
         '       model_solver,model_plot,', &
         '       fimplicit,ntype_mat,eps_mat,', &
         '       dtfactor,dxfactor,dyfactor,', &
         '       ntmax,ntstep,ngtstep,ngrstep,idebuga'
    RETURN
  END SUBROUTINE wq_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE wq_check(ierr)

    USE wqcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    RETURN
  END SUBROUTINE wq_check

!     ***** BROADCAST INPUT PARAMETERS *****

  SUBROUTINE wq_broadcast

    USE wqcomm_parm
    USE libmpi
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: rdata

! --- WQ specific input parameters ---

    idata( 1)=model_geometry
    idata( 2)=nph
    idata( 3)=model_source
    idata( 4)=model_pulse
    idata( 5)=model_ramp
    idata( 6)=medium_max
    idata( 7)=model_solver
    idata( 8)=model_plot
    idata( 9)=ntype_mat
    idata(10)=ntmax
    idata(11)=ntstep
    idata(12)=ngtstep
    idata(13)=ngrstep
    CALL mtx_broadcast_integer(idata,13)
    model_geometry=idata( 1)
    nph=idata( 2)
    model_source=idata( 3)
    model_pulse=idata( 4)
    model_ramp=idata( 5)
    medium_max=idata( 6)
    model_solver=idata( 7)
    model_plot=idata( 8)
    ntype_mat=idata( 9)
    ntmax=idata(10)
    ntstep=idata(11)
    ngtstep=idata(12)
    ngrstep=idata(13)

    rdata( 1)=xnmin
    rdata( 2)=xnmax
    rdata( 3)=ynmin
    rdata( 4)=ynmax
    rdata( 5)=B0
    rdata( 6)=RR
    rdata( 7)=RA
    rdata( 8)=q0
    rdata( 9)=qa
    rdata(10)=freq
    rdata(11)=rkz
    rdata(12)=source_position_xn
    rdata(13)=source_position_yn
    rdata(14)=source_width
    rdata(15)=source_thickness
    rdata(16)=source_angle
    rdata(17)=pulse_length
    rdata(18)=ramp_length
    rdata(19)=dtfactor
    rdata(20)=dxfactor
    rdata(21)=dyfactor
    rdata(22)=fimplicit
    rdata(23)=eps_mat
    CALL mtx_broadcast_real8(rdata,23)
    xnmin=rdata( 1)
    xnmax=rdata( 2)
    ynmin=rdata( 3)
    ynmax=rdata( 4)
    B0=rdata( 5)
    RR=rdata( 6)
    RA=rdata( 7)
    q0=rdata( 8)
    qa=rdata( 9)
    freq=rdata(10)
    rkz=rdata(11)
    source_position_xn=rdata(12)
    source_position_yn=rdata(13)
    source_width=rdata(14)
    source_thickness=rdata(15)
    source_angle=rdata(16)
    pulse_length=rdata(17)
    ramp_length=rdata(18)
    dtfactor=rdata(19)
    dxfactor=rdata(20)
    dyfactor=rdata(21)
    fimplicit=rdata(22)
    eps_mat=rdata(23)

    CALL mtx_broadcast_integer(id_medium,medium_max)
    CALL mtx_broadcast_real8(xnmin_medium,medium_max)
    CALL mtx_broadcast_real8(xnmax_medium,medium_max)
    CALL mtx_broadcast_real8(ynmin_medium,medium_max)
    CALL mtx_broadcast_real8(ynmax_medium,medium_max)
    CALL mtx_broadcast_real8(dielectric_medium,medium_max)
    CALL mtx_broadcast_real8(res_freq_medium,medium_max)
    CALL mtx_broadcast_real8(res_coll_medium,medium_max)
    CALL mtx_broadcast_real8(density_medium,medium_max)
    CALL mtx_broadcast_real8(collision_medium,medium_max)
    CALL mtx_broadcast_integer(idebuga,99)
  END SUBROUTINE wq_broadcast
END MODULE wqparm
