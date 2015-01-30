subroutine wfmenu

  use wfcomm
  use libmpi
  implicit none
  
  integer  :: MODE
  integer,dimension(1) :: idata
  integer  :: IERR
  character:: KID*1,LINE*80

1 continue

  if(nrank.eq.0) then
     write(6,*) '## INPUT: P,V:PARM  D:DIV  A:ANT', &
                         ' W,C:WAVE G:GRAPH  S,L:FILE  Q:QUIT'
     call GUFLSH
     call WFKLIN(LINE,KID,MODE)
  end if

  call mtx_broadcast_character(KID,1)
  if(nrank.eq.0) idata(1)=MODE
  call mtx_broadcast_integer(idata,1)
  MODE=idata(1)
  call mtx_barrier
  if(MODE.ne.1) goto 1

  if     (KID.eq.'P') then
     if(nrank.eq.0) call wfparm(KID)
     call wfparm_broadcast
     goto 1
  elseif (KID.eq.'V') then
     if (nrank.eq.0) call WFVIEW
  elseif (KID.eq.'D') then
     call WFDIV
  elseif (KID.eq.'A') then
!     if (NNMAX.eq.0) call WFRELM(ID)
     call WFANT
  elseif (KID.eq.'C') then
     call WFWPRE(IERR)
  elseif (KID.eq.'W') then
     call WFWAVE
  elseif (KID.eq.'G') then
     if (nrank.eq.0) call WFGOUT
!  elseif (KID.eq.'S') then
!     call WFWFLD
!  elseif (KID.eq.'L') then
!     call WFRFLD
  elseif (KID.eq.'?') then
     if(nrank.eq.0) call WFINFO
  elseif (KID.eq.'Q') then
     goto 9000
  end  if
  KID=' '
  goto 1

9000 continue
  return
end subroutine wfmenu

! ***** INPUT KID or LINE *****
! MODE=0 : LINE INPUT
!      1 : KID INPUT
!      2 : PARM INPUT
!      3 : ERROR

subroutine WFKLIN(LINE,KID,MODE)

  implicit none
  integer,intent(out) :: MODE
  integer :: ID,I
  character,intent(inout) :: LINE*80,KID*1

  read(5,'(A80)',ERR=2,END=3) LINE

  ID=0
  do I=1,80
     if(LINE(I:I).eq.'=') ID=1
  end do
  if(ID.eq.1) then
     call WFPARL(LINE)
     MODE=2
     return
  end if
  
  KID=LINE(1:1)
  call GUCPTL(KID)
  if(KID.eq.'P'.or.&
     KID.eq.'V'.or.&
     KID.eq.'D'.or.&
     KID.eq.'A'.or.&
     KID.eq.'S'.or.&
     KID.eq.'L'.or.&
     KID.eq.'W'.or.&
     KID.eq.'G'.or.&
     KID.eq.'?'.or.&
     KID.eq.'Q') then
     MODE=1
     return
  end if
  
  KID=' '
  MODE=0
  return
  
2 write(6,*) 'XX INPUT ERROR !'
  MODE=3
  return
  
3 KID='Q'
  MODE=1
  return
end subroutine WFKLIN

!     ***** DEBUG INFORMATION ROUTINE *****

subroutine WFINFO
  
  use wfcomm
  implicit none
  integer   :: IE,IN,NN,NE,IA,IS,NSD
  real(8)   :: R,Z
  character :: KID*1
  
8001 continue

  write(6,*) '## INPUT:  E:element  N,V:node  '//&
            ' S:side  A,M:antenna  F:FEP  X:end'
  read(5,'(A1)',ERR=8001,END=9000) KID

  call GUCPTL(KID)
  
  if(KID.eq.'E') then
8002 write(6,*) '## INPUT: Element number '
     read(5,*,ERR=8002,END=8001) IE
     if(IE.eq.0) goto 8001
     write(6,'(A,5I8)') '  NE,KA,KN=',IE,KAELM(IE),&
          KNELM(1,IE),KNELM(2,IE),KNELM(3,IE)
     write(6,*)     '---------------------------------'
     write(6,'(A)') '   NN           R           Z  KA'
     write(6,*)     '---------------------------------'
     do IN=1,3
        NN=NDELM(IN,IE)
        write(6,'(I5,1P2E12.4,1X,I3)') NN,RNODE(NN),ZNODE(NN),KANOD(NN)
     end do
     goto 8002
     
  elseif(KID.eq.'N') then
8003 write(6,*) '## INPUT: Node number '
     read(5,*,ERR=8003,END=8001) NN
     if(NN.eq.0) goto 8001
     write(6,'(A,1P2E12.4,I5)') '   RNODE,ZNODE,KA =',&
          RNODE(NN),ZNODE(NN),KANOD(NN)
     write(6,*) '---------------------------'
     write(6,*) '  Sides including this node'
     write(6,*) '---------------------------'
     do NSD=1,NSDMAX
        do IN=1,2
           if(NDSID(IN,NSD).eq.NN) then
              write(6,'(A,2I8)') '   NSD,IN =',&
                   NSD,IN
           end if
        end do
     end do
     goto 8003
     
  elseif(KID.eq.'A') then
8004 write(6,*) '## INPUT: Antenna number, Segment number '
     read(5,*,ERR=8004,END=8001) IA,IS
     if(IA.eq.0) goto 8001
     write(6,*) '   RJ0,ZJ0 =',&
          RJ0(IS,IA),ZJ0(IS,IA)
     goto 8004
     
  elseif(KID.eq.'M') then
8005 write(6,*) '## INPUT: Antenna number, Segment number '
     read(5,*,ERR=8005,END=8001) IA,IS
     if(IA.eq.0) goto 8001
     write(6,*) '   JELMT,RJ,ZJ =',&
          JELMT(IS,IA),RJ(IS,IA),ZJ(IS,IA)
     goto 8005
     
  elseif(KID.eq.'F') then
8006 write(6,*) '## INPUT: NE,R,Z'
     read(5,*,ERR=8006,END=8001) NE,R,Z
     if(NE.le.0) goto 8001
     if(NE.gt.NEMAX) NE=NEMAX
     call FEP(R,Z,NE)
     write(6,*) '   NE =',NE
     if(NE.ne.0) then
        write(6,*) '  NODE=',&
             NDELM(1,NE),NDELM(2,NE),NDELM(3,NE)
     end if
     goto 8006
     
  elseif(KID.eq.'V') then
     do NN=1,NNMAX
        if(KANOD(NN).gt.0) then
           write(6,'(I5,1P3E12.4)')&
                NN,RNODE(NN),ZNODE(NN)
        end if
     end do

  elseif(KID.eq.'S') then
     do NSD=1,NSDMAX
        write(6,*) NSD,NDSID(1,NSD),NDSID(2,NSD),KASID(NSD)
     end do
     
  elseif(KID.eq.'X') then
     goto 9000
  end if
  goto 8001

9000 return
end subroutine WFINFO
