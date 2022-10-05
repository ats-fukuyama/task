! wqview.f90

MODULE wqview

  PRIVATE
  PUBLIC wq_view

CONTAINS

!     ****** DISPLAY INPUT DATA ******

  SUBROUTINE wq_view
    USE wqcomm_parm
    IMPLICIT NONE
    INTEGER:: i,loop

    WRITE(6,'(A,I8)')     'model_geometry     =',model_geometry
    WRITE(6,'(A,ES12.4)') 'xnmin              =',xnmin
    WRITE(6,'(A,ES12.4)') 'xnmax              =',xnmax
    WRITE(6,'(A,ES12.4)') 'ynmin              =',ynmin
    WRITE(6,'(A,ES12.4)') 'ynmax              =',ynmax

    WRITE(6,'(A,ES12.4)') 'B0                 =',B0
    WRITE(6,'(A,ES12.4)') 'RR                 =',RR
    WRITE(6,'(A,ES12.4)') 'RA                 =',RA
    WRITE(6,'(A,ES12.4)') 'q0                 =',q0
    WRITE(6,'(A,ES12.4)') 'qa                 =',qa

    WRITE(6,'(A,ES12.4)') 'freq               =',freq
    WRITE(6,'(A,ES12.4)') 'rkz                =',rkz
    WRITE(6,'(A,I8)')     'nph                =',nph
    
    WRITE(6,'(A,I8)')     'model_source       =',model_source
    WRITE(6,'(A,ES12.4)') 'source_position_xn =',source_position_xn
    WRITE(6,'(A,ES12.4)') 'source_position_yn =',source_position_yn
    WRITE(6,'(A,ES12.4)') 'source_width       =',source_width
    WRITE(6,'(A,ES12.4)') 'source_thickness   =',source_thickness
    WRITE(6,'(A,ES12.4)') 'source_angle       =',source_angle

    WRITE(6,'(A,I8)')     'model_pulse        =',model_pulse
    WRITE(6,'(A,ES12.4)') 'pulse_length       =',pulse_length
    WRITE(6,'(A,I8)')     'model_ramp         =',model_ramp
    WRITE(6,'(A,ES12.4)') 'ramp_length        =',ramp_length

    WRITE(6,'(A,I8)')     'medium_max =',medium_max
    DO loop=1,(medium_max-1)/5+1
       WRITE(6,'(A,5I12)') &
            'medium             =',&
            (5*(loop-1)+i,i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5(I8,4X))') &
            'id_medium          =',&
            (id_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'xnmin_medium       =', &
            (xnmin_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'xnmax_medium       =', &
            (xnmax_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'ynmin_medium       =', &
            (ynmin_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'ynmax_medium       =', &
            (ynmax_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'dielectric_medium  =', &
            (dielectric_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'res_freq_medium    =', &
            (res_freq_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'res_coll_medium    =', &
            (res_coll_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'density_medium     =', &
            (density_medium(i),i=1,MOD(medium_max-1,5)+1)
       WRITE(6,'(A,5ES12.4)') &
            'collision_medium   =', &
            (collision_medium(i),i=1,MOD(medium_max-1,5)+1)
    END DO

    WRITE(6,'(A,I8)')     'model_solver       =',model_solver
    WRITE(6,'(A,I8)')     'model_plot         =',model_plot
    WRITE(6,'(A,ES12.4)') 'fimplicit          =',fimplicit
    WRITE(6,'(A,I8)')     'ntype_mat          =',ntype_mat
    WRITE(6,'(A,ES12.4)') 'eps_mat            =',eps_mat

    WRITE(6,'(A,ES12.4)') 'dtfactor           =',dtfactor
    WRITE(6,'(A,ES12.4)') 'dxfactor           =',dxfactor
    WRITE(6,'(A,ES12.4)') 'dyfactor           =',dyfactor

    WRITE(6,'(A,I8)')     'ntmax              =',ntmax
    WRITE(6,'(A,I8)')     'ntstep             =',ntstep
    WRITE(6,'(A,I8)')     'ngtstep            =',ngtstep
    WRITE(6,'(A,I8)')     'ngrstep            =',ngrstep

    ! --- idebuga : debug option array:
    !               idebuga(61): matrix equation coefficients
    
    DO i=1,99
       IF(idebuga(i).NE.0) WRITE(6,'(A,I2,A,I8)') &
            'idebuga(',i,')         =',idebuga(i)
    END DO
    
    RETURN
  END SUBROUTINE wq_view
END MODULE wqview
