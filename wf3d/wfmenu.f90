subroutine wfmenu

  use libmpi
  use libmtx
  use wfcomm
  implicit none
  
  integer,dimension(1):: MODE
  integer  :: ID,IERR
  character:: KID*1,LINE*80

1 continue

  if(nrank.eq.0) then
     write(6,601)
601  format(' ','## INPUT: P,V:PARM  D:DIV  A:ANT', &
            '  W,C:WAVE  G:GRAPH  S,L,N:FILE  Y:G10  Q:QUIT')
     call GUFLSH
     call WFKLIN(LINE,KID,MODE)
  end if
  call mtx_barrier
  call mtx_broadcast_character(KID,1)
  call mtx_broadcast_integer(MODE,1)
  if(MODE(1).ne.1) goto 1


  if     (KID.eq.'P') then
     if(nrank.eq.0) call wfparm(KID)
     call wfparm_broadcast
     goto 1
  elseif (KID.eq.'V') then
     if (nrank.eq.0) call WFVIEW
     call mtx_barrier
     if (nrank.eq.1) call wfview
     call mtx_barrier
  elseif (KID.eq.'D') then
     call WFDIV
  elseif (KID.eq.'A') then
     if (NNMAX.eq.0) call WFRELM(ID)
     call WFANT
  elseif (KID.eq.'C') then
     call WFWPRE(IERR)
  elseif (KID.eq.'W') then
     call WFWAVE
  elseif (KID.eq.'G') then
     if (nrank.eq.0) call WFGOUT
     call mtx_barrier
  elseif (KID.eq.'N') then
     call WFNAS
  elseif (KID.eq.'S') then
     call WFWFLD
  elseif (KID.eq.'L') then
     call WFRFLD
  elseif (KID.eq.'?') then
     if(nrank.eq.0) call WFINFO
  elseif (KID.eq.'Q') then
     goto 9000
!
! ----- Add. By YOKOYAMA 25/02/2013 ----
!    KID=Y for GAMMA10 at Univ. of Tsukuba
!
  elseif (KID.eq.'Y') then
!  Loading files for calculation of magnetic field of G10
      call MAGPRP(14,-1)
!  Loading files for saving plot points of G10
      call G10BFL
      write(6,*) '## READ B-FILD LINE DATA FOR GAMMA10'
!
! ----- -----
!
  end  if
  KID=' '
  goto 1

9000 continue
  return
end subroutine wfmenu

!     ***** INPUT KID or LINE *****
!                   MODE=0: LINE INPUT
!                        1: KID INPUT
!                        2: PARM INPUT
!                        3: NEW PROMPT

subroutine WFKLIN(LINE,KID,MODE)

  implicit none
  integer,dimension(1) :: MODE
  integer   :: ID,I
  character :: LINE*80,KID*1

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
     KID.eq.'N'.or.&
     KID.eq.'S'.or.&
     KID.eq.'L'.or.&
     KID.eq.'W'.or.&
     KID.eq.'G'.or.&
     KID.eq.'?'.or.&
! ----- Add. By YOKOYAMA 25/02/2013 ----
     KID.eq.'Y'.or.&
! ----- -----
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

SUBROUTINE WFINFO
  
  use wfcomm
  implicit none
  integer   :: IE,I,IN,NN,NE,NSF,IB,IA,IS,IN1,IN2,IN3,IDEBUGS,L
  real(8)   :: RND,X,Y,Z
  character :: KID*1

8001 write(6,*) '## INPUT:  E:element  N:node  B:boundary  '//&
                       ' A,M:ant  F:FEP  D:EFINDK  X:end'
     read(5,'(A1)',ERR=8001,END=9000) KID
     ! capitalize KID
     call GUCPTL(KID)

     if(KID.eq.'E') then
8002    write(6,*) '## INPUT: Element number '
        read(5,*,ERR=8002,END=8001) IE
        if(IE.eq.0) goto 8001
        write(6,'(A,5I8)') '  NE,KN=',IE,&
                   KNELM(1,IE),KNELM(2,IE),KNELM(3,IE),KNELM(4,IE)
        write(6,'(A,5I8)') '  KA,ND=',KAELM(IE),&
                   NDELM(1,IE),NDELM(2,IE),NDELM(3,IE),NDELM(4,IE)
        do I=1,4
           IN=NDELM(I,IE)
           RND=SQRT(XND(IN)**2+YND(IN)**2)
           write(6,'(A,I5,1P4E12.4,I3,I5)') 'IN,X,Y,Z,R,KA,IM =',&
                &IN,XND(IN),YND(IN),ZND(IN),RND,KANOD(IN),IMLEN(IN)
        end do
        goto 8002
     
     elseif(KID.eq.'N') then
8003      write(6,*) '## INPUT: Node number '
          read(5,*,ERR=8003,END=8001) NN
          if(NN.eq.0) goto 8001
          write(6,'(A,1P3E12.4,2I5)') '   XND,YND,ZND,KA,IM =',&
                       XND(NN),YND(NN),ZND(NN),KANOD(NN),IMLEN(NN)
          do NE=1,NEMAX
             do IN=1,4
                if(NDELM(IN,NE).eq.NN) then
                     write(6,'(A,3I8)') '   NE,IN,KN=',&
                                            NE,IN,KNELM(IN,NE)
                end if
             end do
          end do
          do NSF=1,NSFMAX
             do IN=1,3
                if(NDSRF(IN,NSF).eq.NN) then
                   WRITE(6,'(A,6I8)') '   NSF,IN,KN,ND=',&
                                          NSF,IN,KNSRF(IN,NSF),&
                             NDSRF(1,NSF),NDSRF(2,NSF),NDSRF(3,NSF)
                end if
             end do
          end do
          goto 8003

       elseif(KID.eq.'B') then
8004    write(6,*) '## INPUT: Boundary number '
        read(5,*,ERR=8004,END=8001) IB
        if(IB.eq.0) goto 8001
        write(6,*) '   NDBDY =',NDBDY(IB)
        goto 8004
        
     elseif(KID.eq.'A') then
8005    write(6,*) '## INPUT: Antenna number, Segment number '
        read(5,*,ERR=8005,END=8001) IA,IS
        if(IA.eq.0) goto 8001
        write(6,*) '   XJ0,YJ0,ZJ0 =',&
                       XJ0(IS,IA),YJ0(IS,IA),ZJ0(IS,IA)
        goto 8005
     
     elseif(KID.eq.'M') then
8006    write(6,*) '## INPUT: Antenna number, Segment number '
        read(5,*,ERR=8006,END=8001) IA,IS
        if(IA.eq.0) goto 8001
        write(6,*) '   JELMT,XJ,YJ,ZJ =',&
                       JELMT(IS,IA),XJ(IS,IA),YJ(IS,IA),ZJ(IS,IA)
        goto 8006
     
     elseif(KID.eq.'F') then
8007    write(6,*) '## INPUT: IE,X,Y,Z'
        read(5,*,ERR=8007,END=8001) IE,X,Y,Z
        if(IE.lt.0) goto 8001
        if(IE.gt.NEMAX) IE=NEMAX
        call FEP(X,Y,Z,IE)
        write(6,*) '   IE =',IE
        if(IE.ne.0) then
           write(6,*) '  NDELM=',&
                   NDELM(1,IE),NDELM(2,IE),NDELM(3,IE),NDELM(4,IE)
        end if
        goto 8007

     elseif(KID.eq.'D') then
8008    write(6,*) '## INPUT: IE,IN1,IN2,IN3'
        read(5,*,ERR=8008,END=8001) IE,IN1,IN2,IN3
        if(IE.eq.0) goto 8001
        IDEBUGS=IDEBUG
        IDEBUG=1
        call EFINDK(IE,IN1,IN2,IN3,L)
        IDEBUG=IDEBUGS
        write(6,*) 'IE,L = ',IE,L
        goto 8008

     elseif(KID.eq.'V') then
        do NN=1,NNMAX
           if(KANOD(NN).gt.0) then
              write(6,'(I5,1P3E12.4)')&
                   NN,XND(NN),YND(NN),ZND(NN)
           end if
        end do

     elseif(KID.eq.'X') then
        goto 9000
     end if
     goto 8001

9000 return
end SUBROUTINE WFINFO
