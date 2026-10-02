! wfdiv.f90

MODULE wfdiv

!     ######### /TASK/WF2/WFDIV ########
!
!      MESH DATA GENERATION PROGRAM
!
!     #################################

  PRIVATE
  PUBLIC wf_div

CONTAINS
  
  SUBROUTINE wf_div

    USE libmpi
    USE libchar
    USE wfcomm
    USE wfparm
    USE wfdiv_rect
    USE wfdiv_layer
    USE wfdiv_circle
    USE wfdiv_eq
    USE wfindex
    USE wfgsub
    IMPLICIT NONE
    INTEGER   :: ne,nn,ierr
    CHARACTER :: kid*1

    ! *** div main menu ***
  
1   CONTINUE

    IF(nrank.EQ.0) THEN  
       WRITE(6,*) '## INPUT: D/DIV  G/DRAW  P,V/PARM ',&
                            'S/SAVE  L/LOAD  W/LIST  X/EXIT'
       READ(5,'(A1)',ERR=1,END=9000) kid
       CALL toupper(kid)
    END IF
    CALL mtx_barrier
    CALL mtx_broadcast_character(kid,1)

    SELECT CASE(kid)

    CASE('D')
       ! *** devider menu ***

2      CONTINUE

       kid=""
       IF(nrank.eq.0) THEN
          WRITE(6,'(A24)') '## Tyoe: R/Rect L/Layer C/Circle E/Eq X/exit'
          READ(5,'(A1)') kid
          CALL toupper(kid)
          CALL mtx_barrier
          CALL mtx_broadcast_character(kid,1)
       END IF
     
       SELECT CASE(kid)
       CASE('R')
          model_div=1
       CASE('L')
          model_div=2
       CASE('C')
          model_div=3
       CASE('E')
          model_div=4
       CASE('X')
          GO TO 9000
       CASE DEFAULT
          GO TO 2
       END SELECT

       IF(nrank.EQ.0) THEN
          SELECT CASE(model_div)
          CASE(1)
             CALL wf_div_rect_input
          CASE(2)
             CALL wf_div_layer_input
          CASE(3)
             CALL wf_div_circle_input
          CASE(4)
             CALL wf_div_eq_input
          END SELECT
       END IF
        
       CALL wfdiv_broadcast

       SELECT CASE(model_div)
       CASE(1)
          CALL wf_div_rect_exec
       CASE(2)
          CALL wf_div_layer_exec
       CASE(3)
          CALL wf_div_circle_exec
       CASE(4)
          CALL wf_div_eq_exec
       END SELECT

       CALL fem_mesh_allocate
       IF(nrank.EQ.0) WRITE(6,*) '--- WFINDEX start ---'
       CALL wf_index
       IF(nrank.EQ.0) WRITE(6,*) '--- WFFEPI start ---'
       CALL wf_set_fep
  
       DO NE=1,nelm_max
          nmed_nelm(NE)=1
       END DO
       nmed_max=1
     
       DO NN=1,node_max
          mode_node(NN)=0
       END DO

    CASE('G')
       ! *** element graph ***
       IF(nrank.EQ.0) THEN
          NWXMAX=0
          CALL wf_gr_element
       END IF

    CASE('W')
       ! *** list element data ***
       IF(nrank.EQ.0) CALL wf_list_element
       
     
    CASE('L')
       ! *** load element data ***
       !     CALL wf_load_element(ierr)
     
    CASE('P')
       if(nrank.eq.0) call WF_PARM(0,'WF',IERR)
       call wfparm_broadcast
     
    CASE('V')
       if (nrank.eq.0) call WF_VIEW
     
    CASE('S')
       !     if (nrank.eq.0) call WFWELM(0)
     
    CASE('X')
     goto 9000
  END SELECT
  GO TO 1
  
9000 CONTINUE
  RETURN
END SUBROUTINE wf_div

! *** rectangular mesh ***

subroutine SETNODX
  
  use wfcomm
  implicit none

  integer :: node,nelm
  integer :: node1,node2,node3,node4
  integer :: NX,NY
  integer :: NXMAX,NYMAX
  real(rkind) :: dx,dy,xlen,ylen

  xlen=xdiv_max-xdiv_min
  ylen=ydiv_max-ydiv_min

  ! --- set node_max ---
  nxmax=NINT(xlen/del_xdiv)
  if(MOD(nxmax,2).eq.0) nxmax=nxmax+1
  dx=DBLE(xlen/(nxmax-1))
  nymax=NINT(ylen/del_ydiv)
  if(MOD(nymax,2).eq.0) nymax=nymax+1
  dy=DBLE(ylen/(nxmax-1))
  node_max=nxmax*nymax

  ! --- set nelm_max ---
  nelm_max=2*(nxmax-1)*(nymax-1)

  ! --- set node ---
  node=0
  do ny=1,nymax
     do nx=1,nxmax
        node=node+1
        xnode(node)=xdiv_min+del_xdiv*(nx-1)
        ynode(node)=ydiv_min+del_ydiv*(ny-1)
     end do
  end do
  
  ! --- set element ---
  nelm=0
  do ny=1,nymax-1
     do nx=1,nxmax-1
        node1=nx+nxmax*(ny-1)
        node2=nx+nxmax*(ny-1)+1
        node3=nx+nxmax*ny+1
        node4=nx+nxmax*ny

        if(nx.le.nxmax/2) then
           if(ny.le.nymax/2) then
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node1
              node_nside_nelm(2,nelm)=node3
              node_nside_nelm(3,nelm)=node4         
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node1
              node_nside_nelm(2,nelm)=node2
              node_nside_nelm(3,nelm)=node3
           else
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node1
              node_nside_nelm(2,nelm)=node2
              node_nside_nelm(3,nelm)=node4         
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node2
              node_nside_nelm(2,nelm)=node3
              node_nside_nelm(3,nelm)=node4
           end if
        else
           if(ny.le.nymax/2) then
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node1
              node_nside_nelm(2,nelm)=node2
              node_nside_nelm(3,nelm)=node4
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node2
              node_nside_nelm(2,nelm)=node3
              node_nside_nelm(3,nelm)=node4
           else
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node1
              node_nside_nelm(2,nelm)=node3
              node_nside_nelm(3,nelm)=node4
              nelm=nelm+1
              node_nside_nelm(1,nelm)=node1
              node_nside_nelm(2,nelm)=node2
              node_nside_nelm(3,nelm)=node3
           end if
        end if
     end do
  end do
  return
end subroutine SETNODX

! *** circular mesh ***

subroutine SETNODC

  use wfcomm
  implicit none

  integer :: node,nelm
  integer :: node1,node2,node3,node4
  integer :: NR,NTH,NTH1
  integer :: NRMAX
  integer,dimension(:),pointer :: NTHMAX
  real(rkind) :: RRING,THETA,DR

  ! --- set the number of rings ---
                                
  NRMAX=NINT(RB/del_xdiv)+1
  DR=DBLE(RB/(NRMAX-1))
  allocate(NTHMAX(NRMAX))

  ! --- set node_max & NTHMAX---
                                     
  node=1
  NTHMAX(1)=1
  do NR=2,NRMAX
     NTHMAX(NR)=(NR-1)*6
     DO NTH=1,NTHMAX(NR)
        node=node+1
     end DO
  end do
  node_max=node

  ! --- set nelm_max ---

  nelm=0
  do NR=1,NRMAX-1
     nelm=nelm+6*(2*NR-1)
  end do
  nelm_max=nelm

 CALL fem_base_allocate

  ! --- set node ---
  NR=1
  node=1
  xnode(node)=0.d0
  ynode(node)=0.d0

  DO NR=2,NRMAX
     RRING=DR*DBLE(NR-1)
     DO NTH=1,NTHMAX(NR)
        node=node+1
        THETA=DBLE(NTH-1)*2.d0*PI/DBLE(NTHMAX(NR))
        xnode(node)=RRING*COS(THETA)
        ynode(node)=RRING*SIN(THETA)
     END DO
  END DO

  node_max=node

  IF(modelg.EQ.2) THEN
     DO node=1,node_max
        xnode(node)=xnode(node)+RR
     END DO
  END IF
  
  ! --- set element ---
  
  nelm=0
  NR=1

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

  DO NR=2,NRMAX-1

     NTH=0
     NTH1=0
     
1    continue
     
     NTH=NTH+1
     NTH1=NTH1+1
     
     node1=NTH +INNODE(NR)
     node2=NTH1+INNODE(NR+1)
     node3=NTH1+INNODE(NR+1)+1
     node4=NTH +INNODE(NR)  +1
     if(NTH.eq.NTHMAX(NR)) node4=node4-NTHMAX(NR)
     
     nelm=nelm+1
     node_nside_nelm(1,nelm)=node1
     node_nside_nelm(2,nelm)=node2
     node_nside_nelm(3,nelm)=node3
     nelm=nelm+1
     node_nside_nelm(1,nelm)=node1
     node_nside_nelm(2,nelm)=node3
     node_nside_nelm(3,nelm)=node4
     
     if (mod(NTH,NR-1).eq.0) then
        nelm=nelm+1
        node1=NTH +INNODE(NR)  +1
        node2=NTH1+INNODE(NR+1)+1
        node3=NTH1+INNODE(NR+1)+2
        if(NTH.eq.NTHMAX(NR)) node1=node1-NTHMAX(NR)
        if(NTH1+1.eq.NTHMAX(NR+1)) node3=node3-NTHMAX(NR+1)
        node_nside_nelm(1,nelm)=node1
        node_nside_nelm(2,nelm)=node2
        node_nside_nelm(3,nelm)=node3
        NTH1=NTH1+1
        
     end if
           
     if (NTH1.lt.NTHMAX(NR+1)) goto 1
  END DO

  return
end subroutine SETNODC

! ***** the number of innner-ring-nodes *****

function INNODE(NR)
  implicit none
  integer :: NR,INNODE

  INNODE=3*(NR-1)*(NR-2)+1

  return
end function INNODE


!     ****** List Element Data ******
subroutine wf_list_element

  use wfcomm
  implicit none
  integer :: I,J,NE,ISD,NSD_1,NSD_2

  write(6,100) node_max,(I,xnode(I),ynode(I),I=1,node_max)
100 format(/' ','NODE DATA',7X,'node_max=',I6/&
            ' ',8X,'xnode',10X,'ynode'/&
           (' ',I6,1P2E13.4))

  write(6,*) '----------------------------------------------'
  write(6,200) nseg_max,(I,node_nseg(1,I),node_nseg(2,I),I=1,nseg_max)
200 format(/' ','SIDE DATA',7X,'nseg_max=',I6/&
            ' ',10X,'ND1',4X,'ND2',10X,'ND1',4X,'ND2'/&
           (' ',2(I6,'(',2I6,')')))

  write(6,*) '----------------------------------------------'
  write(6,300) nelm_max,(I,(node_nside_nelm(J,I),J=1,3),I=1,nelm_max)
300 format(' ','ELEMENT NODE DATA',5X,'nelm_max=',I6/&
          (' ',2(I5,'(',3I5,')',2X)))
  
  write(6,*) '----------------------------------------------'
  write(6,*) ' ELEMENT SIDE DATA'
  do NE=1,nelm_max-1,2
     write(6,'(2(A,I4,A,2X))') 'NE=',NE,  '   NSD   ND1   ND2',&
                               'NE=',NE+1,'   NSD   ND1   ND2'
     do ISD=1,3
        NSD_1=nseg_nside_nelm(ISD,NE  )
        NSD_2=nseg_nside_nelm(ISD,NE+1)
        if(NSD_1.ge.0) then
           if(NSD_2.ge.0) then
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,node_nseg(1,NSD_1),node_nseg(2,NSD_1),&
                   NSD_2,node_nseg(1,NSD_2),node_nseg(2,NSD_2)
           else
              NSD_2=-NSD_2
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,node_nseg(1,NSD_1),node_nseg(2,NSD_1),&
                   NSD_2,node_nseg(2,NSD_2),node_nseg(1,NSD_2)
           end if
        else
           NSD_1=-NSD_1
           if(NSD_2.ge.0) then
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,node_nseg(2,NSD_1),node_nseg(1,NSD_1),&
                   NSD_2,node_nseg(1,NSD_2),node_nseg(2,NSD_2)
           else
              NSD_2=-NSD_2
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,node_nseg(2,NSD_1),node_nseg(1,NSD_1),&
                   NSD_2,node_nseg(2,NSD_2),node_nseg(1,NSD_2)
           end if
        end if
     end do
  end do

  write(6,*) '----------------------------------------------'
  write(6,500) mtx_len,(I,mode_node(I),I=1,node_max)
500 format(' ','BOUNDARY NODE DATA',4X,'mtx_len=',I4/&
          (' ',4(I5,'(',I1,') ')))
  return
end subroutine wf_list_element

!-- broadcast data --
subroutine wfdiv_broadcast

  use wfcomm
  use libmpi
  implicit none

  real(rkind),dimension(7) :: rdata

  if (nrank.eq.0) then
     rdata(1)=xdiv_min
     rdata(2)=xdiv_max
     rdata(3)=ydiv_min
     rdata(4)=ydiv_max
     rdata(5)=del_xdiv
     rdata(6)=del_ydiv
     rdata(7)=RB
  end if
  
  call mtx_barrier
  call mtx_broadcast_real8(rdata,7)

  xdiv_min=rdata(1)
  xdiv_max=rdata(2)
  ydiv_min=rdata(3)
  ydiv_max=rdata(4)
  del_xdiv=rdata(5)
  del_ydiv=rdata(6)
  RB    =rdata(7)

  call mtx_broadcast_real8(r_corner,3)
  call mtx_broadcast_real8(z_corner,3)

  return 
end subroutine wfdiv_broadcast

END MODULE wfdiv
