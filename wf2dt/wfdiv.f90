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
  USE wfdiv_layers
  implicit none
  integer   :: NE,NN,IERR
  character :: KID*1

1 continue

  write(6,*) '## INPUT: D/DIV  G/DRAW  P,V/PARM ',&
                       'S/SAVE  L/LOAD  W/LISTw  X/EXIT'
  read(5,'(A1)',ERR=1,END=9000) KID
  call toupper(KID)
  
  SELECT CASE(KID)
  CASE('D') ! divider 
     NZMH = NZM/2

2    CONTINUE

     KID=""
     if (nrank.eq.0) then
        write(6,'(A24)') '## TYPE: X/Rect L/layer C/Circle A/arc'
        read(5,'(A1)') KID
        call toupper(KID)
     end if

     SELECT CASE(KID)
     CASE('X')
3       write(6,'(A,4F10.4)') '## DIV:   BDRMIN,BDRMAX,BDZMIN,BDZMAX = ',&
                                         BDRMIN,BDRMAX,BDZMIN,BDZMAX
        write(6,'(A)')        '## INPUT: BDRMIN,BDRMAX,BDZMIN,BDZMAX ? '
        read(5,*,ERR=3,END=2)            BDRMIN,BDRMAX,BDZMIN,BDZMAX
        write(6,'(A,4F10.4)') '## DIV:   BDRMIN,BDRMAX,BDZMIN,BDZMAX = ',&
                                         BDRMIN,BDRMAX,BDZMIN,BDZMAX
4       write(6,'(A,2F10.4)') '## DIV:   DELR,DELZ = ',DELR,DELZ
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
           
     CASE('L')
        CALL wf_div_layers
        
     CASE DEFAULT
        GO TO 2

     END SELECT

     call wfdiv_broadcast
     SELECT CASE(KID)
     CASE('X')
        call SETNODX
     CASE('C')
        call SETNODC
     END SELECT

     if(nrank.eq.0) write(6,*) '--- WFINDX starts ---'
     call WFINDX
     if(nrank.eq.0) write(6,*) '--- WFFEPI starts ---'
     call WFFEPI
     if(nrank.eq.0) write(6,*) '--- WFFEPI ends ---'
  
     NKMAX=1
     do NE=1,NEMAX
        KAELM(NE)=1
     end do
     NMKA(1)=0
     NMMAX=0
     
     NBMAX=0
     do NN=1,NNMAX
        KANOD(NN)=0
     end do
    
  CASE('G')
     if (nrank.eq.0) then
        NWXMAX=0
        call WFGDIV
     end if
     
  CASE('W')
     if (nrank.eq.0) call WFLDIV
     
  CASE('L')
     !     call WFRELM(ID)
     
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
  go to 1
  
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
  call wfelm_allocate
  call wfelm_allocate

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
  integer,dimension(:),ALLOCATABLE :: NTHMAX
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
  call wfelm_allocate
  call wfelm_allocate

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
