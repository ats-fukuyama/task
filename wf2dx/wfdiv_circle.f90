! wfdiv_circle.f90
  
MODULE wfdiv_circle

  USE wfcomm, &
       nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
       xnode=>rnode,ynode=>znode,itype_mesh=>iddiv

  REAL(rkind):: r_max=1.D0
  REAL(rkind):: x_origin=0.D0
  REAL(rkind):: y_origin=0.D0
  REAL(rkind):: del_r=1.D0

  PRIVATE
  PUBLIC wf_div_circle_input
  PUBLIC wf_div_circle_exec
  
CONTAINS

  SUBROUTINE wf_div_circle_input
    IMPLICIT NONE
    
1   WRITE(6,'(A,3F10.4)') &
         '## DIV:   r_max,x_origin,y_origin: ',r_max,x_origin,y_origin
    WRITE(6,'(A)')        &
         '## INPUT: r_max,x_origin,y_origin ? '
    READ(5,*,ERR=1,END=9) r_max,x_origin,y_origin
    WRITE(6,'(A,3F10.4)') &
         '## DIV:   r_max,x_origin,y_origin: ',r_max,x_origin,y_origin
    IF(ABS(r_max).LE.1.D-6) GOTO 1
    
2   WRITE(6,'(A,3F10.4)') '## DIV:   del_r: ',del_r
    WRITE(6,'(A)')       '## INPUT: del_r ? '
    READ(5,*,ERR=2,END=1)del_r
    IF(ABS(del_r).LE.1.D-6) goto 2
    
    itype_mesh=2
    
9   CONTINUE
    RETURN
  END SUBROUTINE wf_div_circle_input
           
  SUBROUTINE wf_div_circle_exec
    IMPLICIT NONE

    INTEGER :: node,nelm
    INTEGER :: node1,node2,node3,node4
    INTEGER :: nr,nth,nth1
    INTEGER :: nr_max
    INTEGER,ALLOCATABLE :: nthmax_nr(:)
    REAL(rkind) :: r,theta,dr

    ! --- set the number of rings ---
    
    nr_max=NINT(r_max/del_r)+1
    dr=DBLE(r_max/(nr_max-1))
    ALLOCATE(nthmax_nr(nr_max))

    ! --- set node_max & nthmax_nr---
                                     
    node=1 ! origin
    nthmax_nr(1)=1
    DO nr=2,nr_max
       nthmax_nr(nr)=(nr-1)*6
       node=node+nthmax_nr(nr)
    END DO
    node_max=node
    CALL wfelm_allocate

    ! --- set nelm_max ---

    nelm=6
    DO nr=2,nr_max-1
       nelm=nelm+nthmax_nr(nr)+nthmax_nr(nr+1)
    END DO
    nelm_max=nelm
    CALL wfelm_allocate

    ! --- set node coordinates ---
  
    node=1
    nr=1
    xnode(1)=x_origin
    ynode(1)=y_origin

    DO nr=2,nr_max
       r=del_r*DBLE(nr-1)
       DO nth=1,nthmax_nr(nr)
          node=node+1
          theta=DBLE(nth-1)*2.D0*PI/DBLE(nthmax_nr(nr))
          xnode(node)=r*COS(theta)+x_origin
          ynode(node)=r*SIN(theta)+y_origin
       END DO
    END DO
  
    ! --- set element ---
  
    nelm=0
    nr=1

    node1=1
    node2=2
    node3=3
    node4=4

    nelm=nelm+1
    node_nside_nelm(1,nelm)=node1
    node_nside_nelm(2,nelm)=node2
    node_nside_nelm(3,nelm)=node3
    nelm=nelm+1
    node_nside_nelm(1,nelm)=node1
    node_nside_nelm(2,nelm)=node3
    node_nside_nelm(3,nelm)=node4

    node1=1
    node2=4
    node3=5
    node4=6

    nelm=nelm+1
    node_nside_nelm(1,nelm)=node1
    node_nside_nelm(2,nelm)=node2
    node_nside_nelm(3,nelm)=node3
    nelm=nelm+1
    node_nside_nelm(1,nelm)=node1
    node_nside_nelm(2,nelm)=node3
    node_nside_nelm(3,nelm)=node4

    node1=1
    node2=6
    node3=7
    node4=2

    nelm=nelm+1
    node_nside_nelm(1,nelm)=node1
    node_nside_nelm(2,nelm)=node2
    node_nside_nelm(3,nelm)=node3
    nelm=nelm+1
    node_nside_nelm(1,nelm)=node1
    node_nside_nelm(2,nelm)=node3
    node_nside_nelm(3,nelm)=node4

    DO nr=2,nr_max-1

       nth=0
       nth1=0
     
1      CONTINUE
     
       nth=nth+1
       nth1=nth1+1
     
       node1=nth +node_nr(nr)
       node2=nth1+node_nr(nr+1)
       node3=nth1+node_nr(nr+1)+1
       node4=nth +node_nr(nr)  +1
       IF(nth.eq.nthmax_nr(nr)) node4=node4-nthmax_nr(nr)
     
       nelm=nelm+1
       node_nside_nelm(1,nelm)=node1
       node_nside_nelm(2,nelm)=node2
       node_nside_nelm(3,nelm)=node3
       nelm=nelm+1
       node_nside_nelm(1,nelm)=node1
       node_nside_nelm(2,nelm)=node3
       node_nside_nelm(3,nelm)=node4
     
       IF (MOD(nth,nr-1).EQ.0) THEN
          nelm=nelm+1
          node1=nth +node_nr(nr)  +1
          node2=nth1+node_nr(nr+1)+1
          node3=nth1+node_nr(nr+1)+2
          IF(nth.EQ.nthmax_nr(nr)) node1=node1-nthmax_nr(nr)
          if(nth1+1.eq.nthmax_nr(nr+1)) node3=node3-nthmax_nr(nr+1)
          node_nside_nelm(1,nelm)=node1
          node_nside_nelm(2,nelm)=node2
          node_nside_nelm(3,nelm)=node3
          nth1=nth1+1
       END IF
           
       IF(nth1.LT.nthmax_nr(nr+1)) GOTO 1
    END DO

    RETURN
  END SUBROUTINE wf_div_circle_exec

  ! ***** the number of innner-ring-nodes *****
  
  FUNCTION node_nr(nr)
    IMPLICIT NONE
    INTEGER :: nr,node_nr

    node_nr=3*(nr-1)*(nr-2)+1

    RETURN
  END FUNCTION NODE_NR

END MODULE wfdiv_circle
