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
  character :: KID*1,ch*1,char_mode

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
        write(6,'(A)') '## TYPE: X/square C/circle L/layer A/arc'
        read(5,'(A1)') KID
        call toupper(KID)
        if (KID.ne."X".AND. &
            KID.ne."C".AND. &
            KID.NE."L".AND. &
            KID.NE."A") goto 2
     end if
     call mtx_barrier
     call mtx_broadcast_character(KID,1)

     IF(nrank.EQ.0)THEN
        SELECT CASE(KID)
        CASE('X')
3          write(6,'(A,4F10.4)') '## DIV:   BDRMIN,BDRMAX,BDZMIN,BDZMAX = ',&
                                            BDRMIN,BDRMAX,BDZMIN,BDZMAX
           write(6,'(A)')        '## INPUT: BDRMIN,BDRMAX,BDZMIN,BDZMAX ? '
           read(5,*,ERR=3,END=2)            BDRMIN,BDRMAX,BDZMIN,BDZMAX
           write(6,'(A,4F10.4)') '## DIV:   BDRMIN,BDRMAX,BDZMIN,BDZMAX = ',&
                                            BDRMIN,BDRMAX,BDZMIN,BDZMAX
4          write(6,'(A,2F10.4)') '## DIV:   DELR,DELZ = ',DELR,DELZ
           write(6,'(A)')        '## INPUT: DELR,DELZ ? '
           read(5,*,ERR=4,END=2) DELR,DELZ
           write(6,'(A,2F10.4)') '## DIV:   DELR,DELZ = ',DELR,DELZ
           write(6,'(A)')        '## INPUT: DELR,DELZ ? '
           if(abs(DELR).le.1.d-6.or.abs(DELZ).le.1.d-6) goto 2
           iddiv=1
           r_corner(1)=BDRMIN
           z_corner(1)=BDZMIN
           r_corner(2)=BDRMAX
           z_corner(2)=BDZMIN
           r_corner(3)=BDRMIN
           z_corner(3)=BDZMAX
           
        CASE('C')
           
5          write(6,'(A15,F10.4)') '## DIV:   RB = ',RB
           write(6,'(A15)')       '## INPUT: RB ? '
           read(5,*,ERR=5,END=2) RB
           BDRMIN=-RB
           BDRMAX= RB
           BDZMIN=-RB
           BDZMAX= RB
6          write(6,'(A17,F10.4)') '## DIV:   DELR = ',DELR
           write(6,'(A17)')       '## INPUT: DELR ? '
           read(5,*,ERR=6,END=2) DELR
           if(abs(DELR).le.1.D-6) goto 2
           iddiv=2
           r_corner(1)=BDRMIN
           z_corner(1)=BDZMIN
           r_corner(2)=BDRMAX
           z_corner(2)=BDZMIN
           r_corner(3)=BDRMIN
           z_corner(3)=BDZMAX

        CASE('L')
11         CONTINUE
           WRITE(6,'(A)') '## Layers: uniform in X or Y?'
           READ(5,*,ERR=11,END=2) ch
           CALL toupper(ch)
           SELECT CASE(ch)
           CASE('X')
              mode=1
              char_mode='x'
           CASE('Y')
              mode=2
              char_mode='y'
           CASE default
              WRITE(6,*) 'XX input error: X or Y'
              GO TO 11
           END SELECT
12         CONTINUE
           WRITE(6,'(A,I4,A)') '## Number of layers? (Max=',NLM,')'
           READ(5,*,ERR=12,END=11) nlayer_max
           IF(nlayer_max.LE.0.OR.nlayer_max.GT.NLM) THEN
              WRITE(6,*) 'XX Input number of layers:'
              GOTO 12
           END IF
           
13         CONTINUE
           WRITE(6,'(A,A1,A)') &
                '# input minimum ',char_mode,' of layer 1?'
           READ(5,*,ERR=13,END=12) posl_nlayer(1)
           DO nlayer=1,nlayer_max
              WRITE(6,'(A,A1,A,I4,A)') &
                   '# input maximum ',char_mode,' of layer ', &
                   nlayer,' and the step size?'
              READ(5,*,ERR=13,END=12) &
                   posl_nlayer(nlayer+1),step_size_nlayer(nlayer)
              IF(step_size_nlayer(nlayer).LE.0.D0) THEN
                 WRITE(6,*) 'XX wfdiv_L: wrong step_size:', &
                      step_size_nlayer(nlayer)
                 GO TO 11
              END IF
                 
           END DO
14         CONTINUE

           SELECT CASE(mode)
           CASE(1)
15            CONTINUE
              WRITE(6,'(A,3ES12.4)') &
                   '# input y minimum, y maximum, and step size:',&
                   pos_min,pos_max,step_size
              READ(5,*,ERR=15,END=14) pos_min,pos_max,step_size
           CASE(2)
16            CONTINUE
              WRITE(6,'(A,3ES12.4)') &
                   '# input x minimum and x maximum, and step size:', &
                   pos_min,pos_max,step_size
              READ(5,*,ERR=16,END=14) pos_min,pos_max,step_size
           END SELECT

           IF(step_size.LE.0.D0) THEN
              WRITE(6,*) 'XX wfdiv_L: wrong step_size:',step_size
              GO TO 14
           END IF

!           DO nlayer=1,nlayer_max
!              WRITE(6,'(A,I4,3ES12.3)') 'nlayer:', &
!                   nlayer,posl_nlayer(nlayer),posl_nlayer(nlayer+1), &
!                   step_size_nlayer(nlayer)
!           END DO
!           WRITE(6,'(A,3ES12.4)') 'pos_min,max,step_size=', &
!                pos_min,pos_max,step_size

           iddiv=3
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
        END SELECT
     END IF
!     WRITE(6,'(A,4ES12.4)') &
!          'BDRMIN,BDRMAX,BDZMIN,BDZMAX=',BDRMIN,BDRMAX,BDZMIN,BDZMAX

     CALL wfdiv_broadcast

     IF(KID.EQ.'X') THEN
        CALL SETNODX
     ELSEIF(KID.EQ.'C') THEN
        CALL SETNODC
     ELSEIF(KID.EQ.'L') THEN
        call mtx_barrier
        call mtx_broadcast1_integer(mode)
        CALL set_node_L(mode)
     END IF

     if(nrank.eq.0) write(6,*) '--- WFINDX start ---'
     call WFINDX
     if(nrank.eq.0) write(6,*) '--- WFFEPI start ---'
     call WFFEPI
  
     NKMAX=1
     do NE=1,NEMAX
        KAELM(NE)=1
     end do
     NMKA(1)=0
     NMMAX=0
     
!     NBMAX=0
     do NN=1,NNMAX
        KANOD(NN)=0
     end do
    
  elseif(KID.eq.'G') then
     if (nrank.eq.0) then
        NWXMAX=0
        call WFGDIV
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

subroutine SETNODX
  
  use wfcomm
  implicit none

  integer :: NN,NE
  integer :: NN1,NN2,NN3,NN4
  integer :: NR,NZ
  integer :: NRMAX,NZMAX
  real(rkind) :: DR,DZ,BDRLEN,BDZLEN

  BDRLEN=BDRMAX-BDRMIN
  BDZLEN=BDZMAX-BDZMIN

  ! --- set NNMAX ---
  NRMAX=NINT(BDRLEN/DELR)
  if(MOD(NRMAX,2).eq.0) NRMAX=NRMAX+1
  DR=DBLE(BDRLEN/(NRMAX-1))
  NZMAX=NINT(BDZLEN/DELZ)
  if(MOD(NZMAX,2).eq.0) NZMAX=NZMAX+1
  DZ=DBLE(BDZLEN/(NZMAX-1))
  NNMAX=NRMAX*NZMAX

  ! --- set NEMAX ---
  NEMAX=2*(NRMAX-1)*(NZMAX-1)
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
end subroutine SETNODX

! *** circular mesh ***

subroutine SETNODC

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

  ! --- set NNMAX & NTHMAX---
                                     
  NN=1
  NTHMAX(1)=1
  do NR=2,NRMAX
     NTHMAX(NR)=(NR-1)*6
     DO NTH=1,NTHMAX(NR)
        NN=NN+1
     end DO
  end do
  NNMAX=NN

  ! --- set NEMAX ---

  NE=0
  do NR=1,NRMAX-1
     NE=NE+6*(2*NR-1)
  end do
  NEMAX=NE
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
  end subroutine SETNODC

  ! *** layered mesh ***

  SUBROUTINE set_node_L(mode)
    USE wfcomm,nelm_max=>NEMAX,node_max=>NNMAX,xnode=>RNODE,ynode=>ZNODE, &
         xnode_min=>BDRMIN,xnode_max=>BDRMAX, &
         ynode_min=>BDZMIN,ynode_max=>BDZMAX, &
         node_nside_nelm=>NDELM
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    INTEGER,ALLOCATABLE:: nposl_nlayer(:)
    REAL(rkind),ALLOCATABLE:: pos_nposl(:)
    REAL(rkind),ALLOCATABLE:: pos_npos(:)
    INTEGER:: nlayer,nelm,node,node_base,npos,npos_half
    INTEGER:: nposl,nposl_half,nposl_base,nposl_total,npos_max
    REAL(rkind):: pos,step_size_nl,x,y

    npos_max=NINT((pos_max-pos_min)/step_size)
    
    ALLOCATE(nposl_nlayer(nlayer_max))
    
    nposl_total=1
    DO nlayer=1,nlayer_max
       nposl_nlayer(nlayer) &
            =NINT((posl_nlayer(nlayer+1)-posl_nlayer(nlayer)) &
            /step_size_nlayer(nlayer))
       IF(nposl_nlayer(nlayer).LE.1) nposl_nlayer(nlayer)=2
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
       step_size_nl=(posl_nlayer(nlayer+1)-posl_nlayer(nlayer)) &
            /nposl_nlayer(nlayer)
!       WRITE(6,'(A,2I6,2ES12.4)') 'nlayer,nposl,pos,thick:', &
!            nlayer,nposl,pos,step_size_nl
       DO npos=1,nposl_nlayer(nlayer)
          nposl=nposl+1
          pos=pos+step_size_nl
          pos_nposl(nposl)=pos
       END DO
       pos=posl_nlayer(nlayer+1)
       pos_nposl(nposl)=pos
    END DO
!       WRITE(6,'(A,2I6,2ES12.4)') 'nlayer,nposl,pos,thick:', &
!            nlayer_max+1,nposl,pos,step_size

!    DO nposl=1,nposl_total
!       WRITE(6,'(A,I6,ES12.4)') 'nposl,pos:', &
!            nposl,pos_nposl(nposl)
!    END DO
                 
    npos_max=NINT((pos_max-pos_min)/step_size)
    IF(MOD(npos_max,2).EQ.0) npos_max=npos_max+1  ! odd npos_max
 !   WRITE(6,'(A,I6)') 'npos_max=',npos_max
    ALLOCATE(pos_npos(npos_max))
    DO npos=1,npos_max-1
       pos_npos(npos)=pos_min+(pos_max-pos_min)*(npos-1)/(npos_max-1)
    END DO
    pos_npos(npos_max)=pos_max

!    DO npos=1,npos_max
!       WRITE(6,'(A,I6,2ES12.4)') 'npos,pos:',npos,pos_npos(npos)
!    END DO

    node_max=nposl_total*npos_max
    nelm_max=2*(nposl_total-1)*(npos_max-1)
    IF(nrank.EQ.0) WRITE(6,*) 'node_max=',node_max
    IF(nrank.EQ.0) WRITE(6,*) 'nelm_max=',nelm_max

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

!    WRITE(6,*) 'node_max,node:',node_max,node

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
!    WRITE(6,*) 'nelm_max,nelm=',nelm_max,nelm

    RETURN
  END SUBROUTINE set_node_L

  

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

  write(6,100) NNMAX,(I,RNODE(I),ZNODE(I),I=1,NNMAX)
100 format(/' ','NODE DATA',7X,'NNMAX=',I6/&
            ' ',8X,'RNODE',10X,'ZNODE'/&
           (' ',I6,1P2E13.4))

  write(6,*) '----------------------------------------------'
  write(6,200) NSDMAX,(I,NDSID(1,I),NDSID(2,I),I=1,NSDMAX)
200 format(/' ','SIDE DATA',7X,'NSDMAX=',I6/&
            ' ',10X,'ND1',4X,'ND2',10X,'ND1',4X,'ND2'/&
           (' ',2(I6,'(',2I6,')')))

  write(6,*) '----------------------------------------------'
  write(6,300) NEMAX,(I,(NDELM(J,I),J=1,3),I=1,NEMAX)
300 format(' ','ELEMENT NODE DATA',5X,'NEMAX=',I6/&
          (' ',2(I5,'(',3I5,')',2X)))
  
  write(6,*) '----------------------------------------------'
  write(6,*) ' ELEMENT SIDE DATA'
  do NE=1,NEMAX-1,2
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
  write(6,500) MLEN,(I,KANOD(I),I=1,NNMAX)
500 format(' ','BOUNDARY NODE DATA',4X,'MLEN=',I4/&
          (' ',4(I5,'(',I1,') ')))
  return
end subroutine WFLDIV

!-- broadcast data --
subroutine wfdiv_broadcast

  use wfcomm
  use libmpi
  implicit none

  INTEGER:: idata(1)
  real(rkind):: rdata(10)

  IF(nrank.EQ.0) THEN
     idata(1)=nlayer_max
  END IF

  call mtx_barrier
  call mtx_broadcast_integer(idata,1)
  nlayer_max=idata(1)
     

  if (nrank.eq.0) then
     rdata(1)=BDRMIN
     rdata(2)=BDRMAX
     rdata(3)=BDZMIN
     rdata(4)=BDZMAX
     rdata(5)=DELR
     rdata(6)=DELZ
     rdata(7)=RB
     rdata(8)=pos_min
     rdata(9)=pos_max
     rdata(10)=step_size
  end if
  
  call mtx_barrier
  call mtx_broadcast_real8(rdata,10)

  BDRMIN=rdata(1)
  BDRMAX=rdata(2)
  BDZMIN=rdata(3)
  BDZMAX=rdata(4)
  DELR  =rdata(5)
  DELZ  =rdata(6)
  RB    =rdata(7)
  pos_min=rdata(8)
  pos_max=rdata(9)
  step_size=rdata(10)

  call mtx_broadcast_real8(r_corner,3)
  call mtx_broadcast_real8(z_corner,3)

  call mtx_barrier
  call mtx_broadcast_real8(posl_nlayer,nlayer_max+1)
  call mtx_broadcast_real8(step_size_nlayer,nlayer_max)

  return 
end subroutine wfdiv_broadcast
