! wqview.f90

MODULE wqview

  PRIVATE
  PUBLIC wq_view

CONTAINS

!     ****** DISPLAY INPUT DATA ******

  SUBROUTINE wq_view
    USE wqcomm_parm
    IMPLICIT NONE

    WRITE(6,'(A,ES12.4)') 'FREQ      =',FREQ
    WRITE(6,'(A,ES12.4)') 'dtfactor  =',dtfactor
    WRITE(6,'(A,ES12.4)') 'dxfactor  =',dxfactor
    WRITE(6,'(A,ES12.4)') 'dyfactor  =',dyfactor
    WRITE(6,'(A,ES12.4)') 'nufactor  =',nufactor
    WRITE(6,'(A,ES12.4)') 'B0        =',B0
    WRITE(6,'(A,ES12.4)') 'RR        =',RR
    WRITE(6,'(A,ES12.4)') 'RA        =',RA
    WRITE(6,'(A,ES12.4)') 'q0        =',q0
    WRITE(6,'(A,ES12.4)') 'qa        =',qa
    WRITE(6,'(A,ES12.4)') 'n0        =',n0
    WRITE(6,'(A,ES12.4)') 'TMN       =',TMN
    WRITE(6,'(A,I8)')     'ntmax     =',ntmax
    WRITE(6,'(A,I8)')     'nxmax     =',nxmax
    WRITE(6,'(A,I8)')     'nymax     =',nymax
    WRITE(6,'(A,I8)')     'INMODE    =',INMODE
    RETURN
  END SUBROUTINE wq_view
END MODULE wqview
