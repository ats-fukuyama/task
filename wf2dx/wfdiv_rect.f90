! wfdiv_rect.f90
  
MODULE wfdiv_rect

  USE wfcomm

  REAL(rkind):: xr_min=0.D0
  REAL(rkind):: xr_max=1.D0
  REAL(rkind):: yr_min=0.D0
  REAL(rkind):: yr_max=1.D0
  REAL(rkind):: del_xr=1.D-2
  REAL(rkind):: del_yr=1.D-2

  PRIVATE
  PUBLIC wf_div_rect_input
  PUBLIC wf_div_rect_exec
  
CONTAINS

  SUBROUTINE wf_div_rect_input
    IMPLICIT NONE

1   WRITE(6,'(A,4F10.4)') &
         '## DIV:   xr_min,xr_max,yr_min,yr_max: ',xr_min,xr_max,yr_min,yr_max
    WRITE(6,'(A)')        &
         '## INPUT: xr_min,xr_max,yr_min,yr_max ?'
    READ(5,*,ERR=1,END=9) xr_min,xr_max,yr_min,yr_max
    WRITE(6,'(A,4F10.4)') &
         '## DIV:   xr_min,xr_max,yr_min,yr_max: ',xr_min,xr_max,yr_min,yr_max
    IF(ABS(xr_max-xr_min).LE.1.D-6.OR.ABS(yr_max-yr_min).LE.1.D-6) GOTO 1
    
2   WRITE(6,'(A,2F10.4)') &
         '## DIV:   del_xr,del_yr: ',del_xr,del_yr
    WRITE(6,'(A)')        &
         '## INPUT: del_xr,del_yr ?'
    READ(5,*,ERR=2,END=1) del_xr,del_yr
    write(6,'(A,2F10.4)') &
         '## DIV:   del_xr,del_yr = ',del_xr,del_yr
    IF(ABS(del_xr).LE.1.D-6.OR.ABS(del_yr).LE.1.D-6) GOTO 2
    
    mode_mesh=1

9   CONTINUE
    RETURN
  END SUBROUTINE wf_div_rect_input
           
  SUBROUTINE wf_div_rect_exec
    USE femcomm
    IMPLICIT NONE

    INTEGER :: node,nelm
    INTEGER :: node1,node2,node3,node4
    INTEGER :: nxr,nyr
    INTEGER :: nxr_max,nyr_max
    REAL(rkind) :: dxr,dyr,xr_len,yr_len

    xr_len=xr_max-xr_min
    yr_len=yr_max-yr_min

    ! --- set node_max ---
    
    nxr_max=NINT(xr_len/del_xr)

    IF(MOD(nxr_max,2).EQ.0) nxr_max=nxr_max+1 ! nxr_max should be odd number
    dxr=DBLE(xr_len/(nxr_max-1))
    nyr_max=NINT(yr_len/del_yr)
    IF(MOD(nyr_max,2).EQ.0) nyr_max=nyr_max+1 ! nyr_max should be odd number
    dyr=DBLE(yr_len/(nyr_max-1))
    node_max=nxr_max*nyr_max

    ! --- set nelm_max ---
    nelm_max=2*(nxr_max-1)*(nyr_max-1)

    CALL fem_base_allocate

    ! --- set node coordinates ---

    node=0
    DO nyr=1,nyr_max
       DO nxr=1,nxr_max
          node=node+1
          xnode(node)=xr_min+del_xr*(nxr-1)
          ynode(node)=yr_min+del_yr*(nyr-1)
       END DO
    END DO
  
    ! --- set element ---

    nelm=0
    DO nyr=1,nyr_max-1
       DO nxr=1,nxr_max-1
          node1=nxr+nxr_max*(nyr-1)
          node2=nxr+nxr_max*(nyr-1)+1
          node3=nxr+nxr_max*nyr+1
          node4=nxr+nxr_max*nyr

          IF(nelm.GE.nelm_max) WRITE(6,'(A,6I8)') &
               'nxr,nyr,nelm,max=',nxr,nyr,nelm,nxr_max,nyr_max,nelm_max

          IF(nxr.LE.nxr_max/2) THEN
             IF(nyr.LE.nyr_max/2) then
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node1
                node_nside_nelm(2,nelm)=node3
                node_nside_nelm(3,nelm)=node4
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node1
                node_nside_nelm(2,nelm)=node2
                node_nside_nelm(3,nelm)=node3
             ELSE
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node1
                node_nside_nelm(2,nelm)=node2
                node_nside_nelm(3,nelm)=node4
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node2
                node_nside_nelm(2,nelm)=node3
                node_nside_nelm(3,nelm)=node4
             END IF
          ELSE
             IF(nyr.LE.nyr_max/2) then
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node1
                node_nside_nelm(2,nelm)=node2
                node_nside_nelm(3,nelm)=node4
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node2
                node_nside_nelm(2,nelm)=node3
                node_nside_nelm(3,nelm)=node4
             ELSE
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node1
                node_nside_nelm(2,nelm)=node3
                node_nside_nelm(3,nelm)=node4
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node1
                node_nside_nelm(2,nelm)=node2
                node_nside_nelm(3,nelm)=node3
             END IF
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE wf_div_rect_exec
END MODULE wfdiv_rect
