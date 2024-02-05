! wfdiv.f90

!     ######### /TASK/WF2/WFDIV ########
!
!      MESH DATA GENERATION PROGRAM
!
!     #################################

subroutine WFDIV

  use libmpi
  USE libchar
  use wfcomm
  use wfparm
  implicit none
  integer   :: NE,NN,IERR,mode,nlayer
  character(LEN=1) :: KID
  REAL(rkind):: th,x,y,xmin,xmax,ymin,ymax
  INTEGER:: nth,nthmax

1 continue

  if (nrank.eq.0) then  
     write(6,*) '## INPUT: D/DIV  G/DRAW  P,V/PARM ',&
                          'S/SAVE  L/LOAD  W/LIST  X/EXIT'
     read(5,'(A1)',ERR=1,END=9000) KID
     call toupper(KID)
  end if
  call mtx_barrier
  call mtx_broadcast_character(KID,1)
  
  if(KID.eq.'D') then
     NZMH = NZM/2

2    continue

     KID=""
     if (nrank.eq.0) then
        write(6,'(A)') '## TYPE: R/rect A/arc L/layer C/Cricle E/eq'
        read(5,'(A1)') KID
        call toupper(KID)
        if (KID.ne."R".AND. &
            KID.ne."A".AND. &
            KID.NE."L".AND. &
            KID.NE."C".AND. &
            KID.NE."E") goto 2
     end if
     call mtx_barrier
     call mtx_broadcast_character(KID,1)

     IF(nrank.EQ.0)THEN
        SELECT CASE(KID)
        CASE('R')
3          write(6,'(A,4F10.4)') '## DIV:   xmin,xmax,ymin,ymax = ',&
                                            BDRMIN,BDRMAX,BDZMIN,BDZMAX
           write(6,'(A)')        '## INPUT: xmin,xmax,ymin,ymax ? '
           read(5,*,ERR=3,END=2)            BDRMIN,BDRMAX,BDZMIN,BDZMAX
           write(6,'(A,4F10.4)') '## DIV:   xmin,xmax,ymin,ymax = ',&
                                            BDRMIN,BDRMAX,BDZMIN,BDZMAX
4          write(6,'(A,2F10.4)') '## DIV:   delx,dely = ',DELR,DELZ
           write(6,'(A)')        '## INPUT: delx,dely ? '
           read(5,*,ERR=4,END=2) DELR,DELZ
           write(6,'(A,2F10.4)') '## DIV:   delx,dely = ',DELR,DELZ
           write(6,'(A)')        '## INPUT: delx,dely ? '
           if(abs(DELR).le.1.d-6.or.abs(DELZ).le.1.d-6) goto 2
           
           mode_div=1
           r_corner(1)=BDRMIN
           z_corner(1)=BDZMIN
           r_corner(2)=BDRMAX
           z_corner(2)=BDZMIN
           r_corner(3)=BDRMIN
           z_corner(3)=BDZMAX
           
        CASE('A')
13          CONTINUE
           write(6,'(A,4F10.4)') &
                '## DIV:   rmin,rmax,thmin,thmax = ',&
                rmin_div,rmax_div,thmin_div,thmax_div
           write(6,'(A)') &
                '## INPUT: rmin,rmax,thmin,thmax ? '
           read(5,*,ERR=13,END=2) &
                rmin_div,rmax_div,thmin_div,thmax_div
           write(6,'(A,4F10.4)') &
                '## DIV:   rmin,rmax,thmin,yth_ax = ',&
                rmin_div,rmax_div,thmin_div,thmax_div
14         CONTINUE
           write(6,'(A,2F10.4)') &
                '## DIV:   delr,delth = ',delr_div,delth_div
           write(6,'(A)') &
                '## INPUT: delr,delth ? '
           read(5,*,ERR=14,END=2) delr_div,delth_div
           write(6,'(A,2F10.4)') &
                '## DIV:   delr,delth = ',delr_div,delth_div
           if(abs(delr_div).le.1.d-6.or.abs(delth_div).le.1.d-6) goto 2
           mode_div=2
           nthmax=NINT((thmax_div-thmin_div)/delth_div)
           th=thmin_div
           x=rmin_div*COS(th*Pi/180.D0)
           y=rmin_div*SIN(th*Pi/180.D0)
           xmin=x
           xmax=x
           ymin=y
           ymax=y
           DO nth=2,nthmax
              th=thmin_div+delth_div*(nth-1)
              x=rmin_div*COS(th*Pi/180.D0)
              y=rmin_div*SIN(th*Pi/180.D0)
              xmin=MIN(xmin,x)
              xmax=MAX(xmax,x)
              ymin=MIN(ymin,y)
              ymax=MAX(ymax,y)
              x=rmax_div*COS(th*Pi/180.D0)
              y=rmax_div*SIN(th*Pi/180.D0)
              xmin=MIN(xmin,x)
              xmax=MAX(xmax,x)
              ymin=MIN(ymin,y)
              ymax=MAX(ymax,y)
           END DO
           BDRMIN=xmin
           BDRMAX=xmax
           BDZMIN=ymin
           BDZMAX=ymax
           r_corner(1)=BDRMIN
           z_corner(1)=BDZMIN
           r_corner(2)=BDRMAX
           z_corner(2)=BDZMIN
           r_corner(3)=BDRMIN
           z_corner(3)=BDZMAX
           
        CASE('L')
31         CONTINUE
           WRITE(6,'(A,A)') '## Layers: align in X or Y?: '
           READ(5,*,ERR=31,END=2) ch_layer_mode
           CALL toupper(ch_layer_mode)
           SELECT CASE(ch_layer_mode)
           CASE('X')
              mode=1
           CASE('Y')
              mode=2
           CASE default
              WRITE(6,*) 'XX input error: X or Y'
              GO TO 31
           END SELECT
32         CONTINUE
           WRITE(6,'(A)') '## Number of layers?'
           READ(5,*,ERR=32,END=31) nlayer_max
           IF(nlayer_max.LE.0) THEN
              WRITE(6,*) 'XX Input positive number of layers:'
              GOTO 32
           END IF
           
33         CONTINUE
           WRITE(6,'(A,A1,A)') &
                '# input minimum ',ch_layer_mode,' of layer 1?'
           READ(5,*,ERR=33,END=32) posl_nlayer(1)
           DO nlayer=1,nlayer_max
              WRITE(6,'(A,A1,A,I4,A)') &
                   '# input maximum ',ch_layer_mode,' of layer ', &
                   nlayer,' and the step size?'
              READ(5,*,ERR=33,END=32) &
                   posl_nlayer(nlayer+1),thickness_nlayer(nlayer)
              IF(thickness_nlayer(nlayer).LE.0.D0) THEN
                 WRITE(6,*) 'XX wfdiv_L: wrong thickness:', &
                      thickness_nlayer(nlayer)
                 GO TO 31
              END IF
                 
           END DO
34         CONTINUE

           SELECT CASE(mode)
           CASE(1)
35            CONTINUE
              WRITE(6,'(A,3ES12.4)') &
                   '# input y minimum, y maximum, and step size:',&
                   pos_min,pos_max,step_size
              READ(5,*,ERR=35,END=34) pos_min,pos_max,step_size
           CASE(2)
36            CONTINUE
              WRITE(6,'(A,3ES12.4)') &
                   '# input x minimum and x maximum, and step size:', &
                   pos_min,pos_max,step_size
              READ(5,*,ERR=36,END=34) pos_min,pos_max,step_size
           END SELECT

           IF(step_size.LE.0.D0) THEN
              WRITE(6,*) 'XX wfdiv_L: wrong step_size:',step_size
              GO TO 34
           END IF

           DO nlayer=1,nlayer_max
              WRITE(6,'(A,I4,3ES12.3)') 'nlayer:', &
                   nlayer,posl_nlayer(nlayer),posl_nlayer(nlayer+1), &
                   thickness_nlayer(nlayer)
           END DO
           WRITE(6,'(A,3ES12.4)') 'pos_min,max,step_size=', &
                pos_min,pos_max,step_size

           mode_div=3
           SELECT CASE(mode)
           CASE(1)
              BDRMIN=posl_nlayer(1)
              BDRMAX=posl_nlayer(nlayer_max+1)
              BDZMIN=pos_min
              BDZMAX=pos_max
           CASE(2)
              BDZMIN=posl_nlayer(1)
              BDZMAX=posl_nlayer(nlayer_max+1)
              BDRMIN=pos_min
              BDRMAX=pos_max
           END SELECT
           
           r_corner(1)=BDRMIN
           z_corner(1)=BDZMIN
           r_corner(2)=BDRMAX
           z_corner(2)=BDZMIN
           r_corner(3)=BDRMIN
           z_corner(3)=BDZMAX
           
        CASE('C')
45         CONTINUE
           write(6,'(A15,F10.4)') '## DIV:   RB = ',RB
           write(6,'(A15)')       '## INPUT: RB ? '
           read(5,*,ERR=45,END=2) RB
46         CONTINUE
           write(6,'(A17,F10.4)') '## DIV:   DELR = ',DELR
           write(6,'(A17)')       '## INPUT: DELR ? '
           read(5,*,ERR=46,END=2) DELR
           if(abs(DELR).le.1.D-8) goto 2
           mode_div=4
           BDRMIN=-RB
           BDRMAX= RB
           BDZMIN=-RB
           BDZMAX= RB
           r_corner(1)=BDRMIN
           z_corner(1)=BDZMIN
           r_corner(2)=BDRMAX
           z_corner(2)=BDZMIN
           r_corner(3)=BDRMIN
           z_corner(3)=BDZMAX

        END SELECT
     END IF

     CALL wfdiv_broadcast
     
     IF(KID.EQ.'R') THEN
        CALL set_node_rect
     ELSEIF(KID.EQ.'A') THEN
        CALL set_node_arc
     ELSEIF(KID.EQ.'L') THEN
        call mtx_barrier
        call mtx_broadcast1_integer(mode)
        CALL set_node_layer(mode)
     ELSEIF(KID.EQ.'C') THEN
        CALL set_node_circle
     END IF

     if(nrank.eq.0) write(6,*) '--- WFINDX start ---'
     call WFINDX
     if(nrank.eq.0) write(6,*) '--- WFFEPI start ---'
     call WFFEPI
  
     NKMAX=1
     do NE=1,nelm_max
        KAELM(NE)=1
     end do
     NMKA(1)=0
     NMMAX=0
     
     do NN=1,node_max
        KANOD(NN)=0
     end do
    
  elseif(KID.eq.'G') then
     if (nrank.eq.0) then
        NWXMAX=0
        call wf_gr_element
     end if
     
  elseif(KID.eq.'W') then
     if (nrank.eq.0) call WFLDIV
     
  elseif(KID.eq.'L') then
     !     call WFRELM(ID)
     
  elseif(KID.eq.'P') then
     if(nrank.eq.0) call WF_PARM(0,'WF',IERR)
     call wfparm_broadcast
     
  elseif(KID.eq.'V') then
     if (nrank.eq.0) call WF_VIEW
     
  elseif(KID.eq.'S') then
     !     if (nrank.eq.0) call WFWELM(0)
     
  elseif(KID.eq.'X') then
     goto 9000
  end if
  goto 1
  
9000 continue
  return
end subroutine WFDIV

! *** rectangular mesh ***

SUBROUTINE set_node_rect
  
  use wfcomm
  implicit none

  integer :: NN,NE
  integer :: NN1,NN2,NN3,NN4
  integer :: NR,NZ
  integer :: NRMAX,NZMAX
  real(rkind) :: DR,DZ,BDRLEN,BDZLEN

  BDRLEN=BDRMAX-BDRMIN
  BDZLEN=BDZMAX-BDZMIN

  ! --- set node_max ---
  NRMAX=NINT(BDRLEN/DELR)+1
  if(MOD(NRMAX,2).eq.0) NRMAX=NRMAX+1
  DR=DBLE(BDRLEN/(NRMAX-1))
  NZMAX=NINT(BDZLEN/DELZ)+1
  if(MOD(NZMAX,2).eq.0) NZMAX=NZMAX+1
  DZ=DBLE(BDZLEN/(NZMAX-1))
  node_max=NRMAX*NZMAX

  ! --- set nelm_max ---
  nelm_max=2*(NRMAX-1)*(NZMAX-1)
  call wf_node_allocate
  call wf_elm_allocate

  ! --- set node ---
  NN=0
  do NZ=1,NZMAX
     do NR=1,NRMAX
        NN=NN+1
        RNODE(NN)=BDRMIN+DR*(NR-1)
        ZNODE(NN)=BDZMIN+DZ*(NZ-1)
     end do
  end do
  
  ! --- set element ---
  NE=0
  do NZ=1,NZMAX-1
     do NR=1,NRMAX-1
        NN1=NR+NRMAX*(NZ-1)
        NN2=NR+NRMAX*(NZ-1)+1
        NN3=NR+NRMAX*NZ+1
        NN4=NR+NRMAX*NZ

        if(NR.le.NRMAX/2) then
           if(NZ.le.NZMAX/2) then
              NE=NE+1
              NDELM(1,NE)=NN1
              NDELM(2,NE)=NN3
              NDELM(3,NE)=NN4         
              NE=NE+1
              NDELM(1,NE)=NN1
              NDELM(2,NE)=NN2
              NDELM(3,NE)=NN3
           else
              NE=NE+1
              NDELM(1,NE)=NN1
              NDELM(2,NE)=NN2
              NDELM(3,NE)=NN4         
              NE=NE+1
              NDELM(1,NE)=NN2
              NDELM(2,NE)=NN3
              NDELM(3,NE)=NN4
           end if
        else
           if(NZ.le.NZMAX/2) then
              NE=NE+1
              NDELM(1,NE)=NN1
              NDELM(2,NE)=NN2
              NDELM(3,NE)=NN4
              NE=NE+1
              NDELM(1,NE)=NN2
              NDELM(2,NE)=NN3
              NDELM(3,NE)=NN4
           else
              NE=NE+1
              NDELM(1,NE)=NN1
              NDELM(2,NE)=NN3
              NDELM(3,NE)=NN4
              NE=NE+1
              NDELM(1,NE)=NN1
              NDELM(2,NE)=NN2
              NDELM(3,NE)=NN3
           end if
        end if
     end do
  end do
  return
end subroutine set_node_rect

  
! *** arc-type mesh ***

subroutine set_node_arc
  
  use wfcomm
  implicit none

  integer :: node,nelm
  integer :: node1,node2,node3,node4
  integer :: nr,nth
  integer :: nrmax,nthmax
  real(rkind) :: rsize,thsize,dr,dth,r,th

  ! --- set node_max ---

  rsize=rmax_div-rmin_div
  thsize=thmax_div-thmin_div

  nrmax=NINT(rsize/delr_div)+1
  nthmax=NINT(thsize/delth_div)+1
  
  IF(MOD(nrmax,2).EQ.0) nrmax=nrmax+1    ! nrmax : odd number
  dr=DBLE(rsize/(nrmax-1))
  IF(MOD(nthmax,2).EQ.0) nthmax=nthmax+1 ! nthmax : odd number
  dth=DBLE(thsize/(nthmax-1))
  node_max=nrmax*nthmax

  ! --- set nelm_max ---
  nelm_max=2*(nrmax-1)*(nthmax-1)
  CALL wf_node_allocate
  CALL wf_elm_allocate

  ! --- set node ---
  node=0
  DO nth=1,nthmax
     DO nr=1,nrmax
        node=node+1
        r=rmin_div+dr*(nr-1)
        th=(thmin_div+dth*(nth-1))*Pi/180.D0
        RNODE(node)=RR+r*COS(th)
        ZNODE(node)=   r*SIN(th)
     END DO
  END DO
  
  ! --- set element ---
  nelm=0
  do nth=1,nthmax-1
     do nr=1,nrmax-1
        node1=nr+nrmax*(nth-1)
        node2=nr+nrmax*(nth-1)+1
        node3=nr+nrmax*nth+1
        node4=nr+nrmax*nth

        if(nr.LE.nrmax/2) then
           if(nth.LE.nthmax/2) then
              nelm=nelm+1
              NDELM(1,nelm)=node1
              NDELM(2,nelm)=node3
              NDELM(3,nelm)=node4         
              nelm=nelm+1
              NDELM(1,nelm)=node1
              NDELM(2,nelm)=node2
              NDELM(3,nelm)=node3
           else
              nelm=nelm+1
              NDELM(1,nelm)=node1
              NDELM(2,nelm)=node2
              NDELM(3,nelm)=node4         
              nelm=nelm+1
              NDELM(1,nelm)=node2
              NDELM(2,nelm)=node3
              NDELM(3,nelm)=node4
           end if
        else
           if(nth.le.nthmax/2) then
              nelm=nelm+1
              NDELM(1,nelm)=node1
              NDELM(2,nelm)=node2
              NDELM(3,nelm)=node4
              nelm=nelm+1
              NDELM(1,nelm)=node2
              NDELM(2,nelm)=node3
              NDELM(3,nelm)=node4
           else
              nelm=nelm+1
              NDELM(1,nelm)=node1
              NDELM(2,nelm)=node3
              NDELM(3,nelm)=node4
              nelm=nelm+1
              NDELM(1,nelm)=node1
              NDELM(2,nelm)=node2
              NDELM(3,nelm)=node3
           end if
        end if
     end do
  end do
  return
end subroutine set_node_arc

  ! *** layered mesh ***

  SUBROUTINE set_node_layer(mode)
    USE wfcomm,xnode=>RNODE,ynode=>ZNODE, &
         xnode_min=>BDRMIN,xnode_max=>BDRMAX, &
         ynode_min=>BDZMIN,ynode_max=>BDZMAX, &
         node_nside_nelm=>NDELM
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    INTEGER,ALLOCATABLE:: nposl_nlayer(:)
    REAL(rkind),ALLOCATABLE:: pos_nposl(:)
    REAL(rkind),ALLOCATABLE:: pos_npos(:)
    INTEGER:: nlayer,nelm,node,node_base,npos,npos_half,npos_max
    INTEGER:: nposl,nposl_half,nposl_base,nposl_total
    REAL(rkind):: pos,thickness,x,y

    ALLOCATE(nposl_nlayer(nlayer_max))
    
    nposl_total=1
    DO nlayer=1,nlayer_max
       nposl_nlayer(nlayer) &
            =NINT((posl_nlayer(nlayer+1)-posl_nlayer(nlayer)) &
            /thickness_nlayer(nlayer))+1
       IF(MOD(nposl_nlayer(nlayer),2).EQ.1) &
            nposl_nlayer(nlayer)=nposl_nlayer(nlayer)+1 ! nposl_nlayer: even
       nposl_total=nposl_total+nposl_nlayer(nlayer)
    END DO   ! odd nposl_total

!    WRITE(6,'(A,I6)') 'nposl_total=',nposl_total
!    DO nlayer=1,nlayer_max
!       WRITE(6,'(A,2I6)') 'nlayer,nposl:',nlayer,nposl_nlayer(nlayer)
!    END DO
              
    ALLOCATE(pos_nposl(nposl_total))
    
    nposl=1
    pos=posl_nlayer(1)
    pos_nposl(nposl)=pos
    DO nlayer=1,nlayer_max
       thickness=(posl_nlayer(nlayer+1)-posl_nlayer(nlayer)) &
            /nposl_nlayer(nlayer)
       WRITE(6,'(A,2I6,2ES12.4)') 'nlayer,nposl,pos,thick:', &
            nlayer,nposl,pos,thickness
       DO npos=1,nposl_nlayer(nlayer)
          nposl=nposl+1
          pos=pos+thickness
          pos_nposl(nposl)=pos
       END DO
       pos=posl_nlayer(nlayer+1)
       pos_nposl(nposl)=pos
    END DO

    DO nposl=1,nposl_total
       WRITE(6,'(A,I6,ES12.4)') 'nposl,pos:', &
            nposl,pos_nposl(nposl)
    END DO
                 
    npos_max=NINT((pos_max-pos_min)/step_size)
    IF(MOD(npos_max,2).EQ.0) npos_max=npos_max+1  ! odd npos_max
    WRITE(6,'(A,I6)') 'npos_max=',npos_max
    ALLOCATE(pos_npos(npos_max))
    DO npos=1,npos_max-1
       pos_npos(npos)=pos_min+(pos_max-pos_min)*(npos-1)/(npos_max-1)
    END DO
    pos_npos(npos_max)=pos_max

    DO npos=1,npos_max
       WRITE(6,'(A,I6,2ES12.4)') 'npos,pos:',npos,pos_npos(npos)
    END DO

    node_max=nposl_total*npos_max
    nelm_max=2*(nposl_total-1)*(npos_max-1)
    WRITE(6,*) 'node_max=',node_max
    WRITE(6,*) 'nelm_max=',nelm_max

    CALL wf_node_allocate
    CALL wf_elm_allocate

    SELECT CASE(mode)
    CASE(1)
       node=0
       DO nposl=1,nposl_total
          x=pos_nposl(nposl)
          DO npos=1,npos_max
             node=node+1
             xnode(node)=x
             ynode(node)=pos_npos(npos)
          END DO
       END DO
    CASE(2)
       node=0
       DO nposl=1,nposl_total
          y=pos_nposl(nposl+1)
          DO npos=1,npos_max
             node=node+1
             xnode(node)=pos_npos(npos)
             ynode(node)=y
          END DO
       END DO
    END SELECT

    WRITE(6,*) 'node_max,node:',node_max,node

    SELECT CASE(mode)
    CASE(1)
       npos_half=(npos_max-1)/2
       nposl_base=0
       nelm=0
       DO nlayer=1,nlayer_max
          nposl_half=nposl_nlayer(nlayer)/2
          DO nposl=nposl_base+1,nposl_base+nposl_half
             node_base=npos_max*(nposl-1)
             DO npos=2,npos_half
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(3,nelm)=node_base+npos_max+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(3,nelm)=node_base+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos
             END DO
             DO npos=npos_half+1,npos_max
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(3,nelm)=node_base+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos-1
             END DO
          END DO
          DO nposl=nposl_base+nposl_half+1,nposl_base+2*nposl_half
             node_base=npos_max*(nposl-1)
             DO npos=2,npos_half
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(3,nelm)=node_base+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos-1
             END DO
             DO npos=npos_half+1,npos_max
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(3,nelm)=node_base+npos_max+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(3,nelm)=node_base+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos
             END DO
          END DO
          nposl_base=nposl_base+nposl_nlayer(nlayer)
       END DO
    CASE(2)
       npos_half=(npos_max-1)/2
       nposl_base=0
       nelm=0
       DO nlayer=1,nlayer_max
          nposl_half=nposl_nlayer(nlayer)/2
          DO nposl=nposl_base+1,nposl_base+nposl_half
             node_base=npos_max*(nposl-1)
             DO npos=2,npos_half
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(2,nelm)=node_base+npos_max+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(2,nelm)=node_base+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos
             END DO
             DO npos=npos_half+1,npos_max
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(2,nelm)=node_base+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos-1
             END DO
          END DO
          DO nposl=nposl_base+nposl_half+1,nposl_base+2*nposl_half
             node_base=npos_max*(nposl-1)
             DO npos=2,npos_half
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(2,nelm)=node_base+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos
                node_nside_nelm(2,nelm)=node_base+npos_max+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos-1
             END DO
             DO npos=npos_half+1,npos_max
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(2,nelm)=node_base+npos_max+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos-1
                nelm=nelm+1
                node_nside_nelm(1,nelm)=node_base+npos-1
                node_nside_nelm(2,nelm)=node_base+npos
                node_nside_nelm(3,nelm)=node_base+npos_max+npos
             END DO
          END DO
          nposl_base=nposl_base+nposl_nlayer(nlayer)
       END DO
    END SELECT
    WRITE(6,*) 'nelm_max,nelm=',nelm_max,nelm

    SELECT CASE(mode)
    CASE(1)
       xnode_min=pos_npos(1)
       xnode_max=pos_npos(npos_max)
       ynode_min=pos_nposl(1)
       ynode_max=pos_nposl(nposl_total)
    CASE(2)
       xnode_min=pos_nposl(1)
       xnode_max=pos_nposl(nposl_total)
       ynode_min=pos_npos(1)
       ynode_max=pos_npos(npos_max)
    END SELECT
    RETURN
  END SUBROUTINE set_node_layer

! *** circular mesh ***

subroutine set_node_circle

  use wfcomm
  implicit none

  integer :: NN,NE
  integer :: NN1,NN2,NN3,NN4
  integer :: NR,NTH,NTH1
  integer :: NRMAX,INNODE
  integer,dimension(:),pointer :: NTHMAX
  real(rkind) :: RRING,THETA,DR

  ! --- set the number of rings ---
                                
  NRMAX=NINT(RB/DELR)+1
  DR=DBLE(RB/(NRMAX-1))
  allocate(NTHMAX(NRMAX))

  ! --- set node_max & NTHMAX---
                                     
  NN=1
  NTHMAX(1)=1
  do NR=2,NRMAX
     NTHMAX(NR)=(NR-1)*6
     DO NTH=1,NTHMAX(NR)
        NN=NN+1
     end DO
  end do
  node_max=NN

  ! --- set nelm_max ---

  NE=0
  do NR=1,NRMAX-1
     NE=NE+6*(2*NR-1)
  end do
  nelm_max=NE
  call wf_node_allocate
  call wf_elm_allocate

  ! --- set node ---
  NN=1
  NR=1
  RNODE(1)=0.d0 + RR
  ZNODE(1)=0.d0

  do NR=2,NRMAX
     RRING=DR*DBLE(NR-1)
     DO NTH=1,NTHMAX(NR)
        NN=NN+1
        THETA=DBLE(NTH-1)*2.d0*PI/DBLE(NTHMAX(NR))
        RNODE(NN)=RRING*cos(THETA) + RR
        ZNODE(NN)=RRING*sin(THETA)
     END DO
  END do
  
  ! --- set element ---
  
  NE=0
  NR=1

  NN1=1
  NN2=2
  NN3=3
  NN4=4

  NE=NE+1
  NDELM(1,NE)=NN1
  NDELM(2,NE)=NN2
  NDELM(3,NE)=NN3
  NE=NE+1
  NDELM(1,NE)=NN1
  NDELM(2,NE)=NN3
  NDELM(3,NE)=NN4

  NN1=1
  NN2=4
  NN3=5
  NN4=6

  NE=NE+1
  NDELM(1,NE)=NN1
  NDELM(2,NE)=NN2
  NDELM(3,NE)=NN3
  NE=NE+1
  NDELM(1,NE)=NN1
  NDELM(2,NE)=NN3
  NDELM(3,NE)=NN4

  NN1=1
  NN2=6
  NN3=7
  NN4=2

  NE=NE+1
  NDELM(1,NE)=NN1
  NDELM(2,NE)=NN2
  NDELM(3,NE)=NN3
  NE=NE+1
  NDELM(1,NE)=NN1
  NDELM(2,NE)=NN3
  NDELM(3,NE)=NN4

  DO NR=2,NRMAX-1

     NTH=0
     NTH1=0
     
1    continue
     
     NTH=NTH+1
     NTH1=NTH1+1
     
     NN1=NTH +INNODE(NR)
     NN2=NTH1+INNODE(NR+1)
     NN3=NTH1+INNODE(NR+1)+1
     NN4=NTH +INNODE(NR)  +1
     if(NTH.eq.NTHMAX(NR)) NN4=NN4-NTHMAX(NR)
     
     NE=NE+1
     NDELM(1,NE)=NN1
     NDELM(2,NE)=NN2
     NDELM(3,NE)=NN3
     NE=NE+1
     NDELM(1,NE)=NN1
     NDELM(2,NE)=NN3
     NDELM(3,NE)=NN4
     
     if (mod(NTH,NR-1).eq.0) then
        NE=NE+1
        NN1=NTH +INNODE(NR)  +1
        NN2=NTH1+INNODE(NR+1)+1
        NN3=NTH1+INNODE(NR+1)+2
        if(NTH.eq.NTHMAX(NR)) NN1=NN1-NTHMAX(NR)
        if(NTH1+1.eq.NTHMAX(NR+1)) NN3=NN3-NTHMAX(NR+1)
        NDELM(1,NE)=NN1
        NDELM(2,NE)=NN2
        NDELM(3,NE)=NN3
        NTH1=NTH1+1
        
     end if
           
     if (NTH1.lt.NTHMAX(NR+1)) goto 1
  END DO

  return
END SUBROUTINE set_node_circle
  

! ***** the number of innner-ring-nodes *****
function INNODE(NR)
  implicit none
  integer :: NR,INNODE

  INNODE=3*(NR-1)*(NR-2)+1

  return
end function INNODE


!     ****** List Element Data ******
subroutine WFLDIV

  use wfcomm
  implicit none
  integer :: I,J,NE,ISD,NSD_1,NSD_2

  write(6,100) node_max,(I,RNODE(I),ZNODE(I),I=1,node_max)
100 format(/' ','NODE DATA',7X,'node_max=',I6/&
            ' ',8X,'RNODE',10X,'ZNODE'/&
           (' ',I6,1P2E13.4))

  write(6,*) '----------------------------------------------'
  write(6,200) NSDMAX,(I,NDSID(1,I),NDSID(2,I),I=1,NSDMAX)
200 format(/' ','SIDE DATA',7X,'NSDMAX=',I6/&
            ' ',10X,'ND1',4X,'ND2',10X,'ND1',4X,'ND2'/&
           (' ',2(I6,'(',2I6,')')))

  write(6,*) '----------------------------------------------'
  write(6,300) nelm_max,(I,(NDELM(J,I),J=1,3),I=1,nelm_max)
300 format(' ','ELEMENT NODE DATA',5X,'nelm_max=',I6/&
          (' ',2(I5,'(',3I5,')',2X)))
  
  write(6,*) '----------------------------------------------'
  write(6,*) ' ELEMENT SIDE DATA'
  do NE=1,nelm_max-1,2
     write(6,'(2(A,I4,A,2X))') 'NE=',NE,  '   NSD   ND1   ND2',&
                               'NE=',NE+1,'   NSD   ND1   ND2'
     do ISD=1,3
        NSD_1=NSDELM(ISD,NE  )
        NSD_2=NSDELM(ISD,NE+1)
        if(NSD_1.ge.0) then
           if(NSD_2.ge.0) then
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,NDSID(1,NSD_1),NDSID(2,NSD_1),&
                   NSD_2,NDSID(1,NSD_2),NDSID(2,NSD_2)
           else
              NSD_2=-NSD_2
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,NDSID(1,NSD_1),NDSID(2,NSD_1),&
                   NSD_2,NDSID(2,NSD_2),NDSID(1,NSD_2)
           end if
        else
           NSD_1=-NSD_1
           if(NSD_2.ge.0) then
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,NDSID(2,NSD_1),NDSID(1,NSD_1),&
                   NSD_2,NDSID(1,NSD_2),NDSID(2,NSD_2)
           else
              NSD_2=-NSD_2
              write(6,'(7X,3I6,9X,4I6)')&
                   NSD_1,NDSID(2,NSD_1),NDSID(1,NSD_1),&
                   NSD_2,NDSID(2,NSD_2),NDSID(1,NSD_2)
           end if
        end if
     end do
  end do

  write(6,*) '----------------------------------------------'
  write(6,500) MLEN,(I,KANOD(I),I=1,node_max)
500 format(' ','BOUNDARY NODE DATA',4X,'MLEN=',I4/&
          (' ',4(I5,'(',I1,') ')))
  return
end subroutine WFLDIV

!-- broadcast data --
subroutine wfdiv_broadcast

  use wfcomm
  use libmpi
  implicit none

  real(rkind),dimension(7) :: rdata

  if (nrank.eq.0) then
     rdata(1)=BDRMIN
     rdata(2)=BDRMAX
     rdata(3)=BDZMIN
     rdata(4)=BDZMAX
     rdata(5)=DELR
     rdata(6)=DELZ
     rdata(7)=RB
  end if
  
  call mtx_barrier
  call mtx_broadcast_real8(rdata,7)

  BDRMIN=rdata(1)
  BDRMAX=rdata(2)
  BDZMIN=rdata(3)
  BDZMAX=rdata(4)
  DELR  =rdata(5)
  DELZ  =rdata(6)
  RB    =rdata(7)

  call mtx_broadcast_real8(r_corner,3)
  call mtx_broadcast_real8(z_corner,3)

  return 
end subroutine wfdiv_broadcast
