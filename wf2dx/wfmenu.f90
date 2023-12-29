! wfmenu.f90

MODULE wfmenu

  PRIVATE
  PUBLIC wf_menu

CONTAINS

  SUBROUTINE wf_menu

  use libmpi
  use wfcomm
  USE plparm
  USE dpparm
  USE wfparm
  USE plload, ONLY: pl_load
  USE wfload, ONLY: wf_load_wg
  USE wfdiv
  USE wfant
  USE wfwave
  USE wfgout
  USE libkio
  implicit none
  
  integer  :: MODE
  integer  :: IERR
  character:: KID*1,LINE*80
  INTEGER:: IDEBUG_SAVE

1 continue

  IF(nrank.EQ.0) THEN
     WRITE(6,*) '## INPUT: P,V:parm  D:div  A:ant', &
                         ' R:run C:check G:graph  S,L:file  Q:QUIT'
     CALL TASK_KLIN(LINE,KID,MODE,WF_PARM)
  END IF
  call mtx_barrier

  call mtx_broadcast1_character(KID)
  call mtx_broadcast1_integer(MODE)
  if(MODE.ne.1) goto 1

  if     (KID.eq.'P') then
     if(nrank.eq.0) call wf_parm(0,'wf',IERR)
     call pl_broadcast
     call dp_broadcast
     call wfparm_broadcast
     goto 1
  elseif (KID.eq.'V') then
     if (nrank.eq.0) call WF_VIEW
  elseif (KID.eq.'D') then
     call wf_div
  elseif (KID.eq.'A') then
     call wf_ant
  elseif (KID.eq.'C') then
     call wf_wpre(IERR)
  elseif (KID.eq.'R') then
     call wf_wave
  elseif (KID.eq.'G') then
     if (nrank.eq.0) call wf_gout
  elseif (KID.eq.'L') then
     IDEBUG_SAVE=IDEBUG
     IDEBUG=1
     CALL pl_load(ierr)
     CALL wf_load_wg(ierr)
     IDEBUG=IDEBUG_SAVE
  elseif (KID.eq.'?') then
     if(nrank.eq.0) call WFINFO
  elseif (KID.eq.'Q') then
     goto 9000
  end  if
  KID=' '
  goto 1

9000 continue
  if(node_max_save.ne.0) call wf_node_deallocate
  if(nelm_max_save.ne.0) call wf_nelm_deallocate
  if(nseg_max_save.ne.0) call wf_nseg_deallocate
  if(nelm_max_sort_save.NE.0) call wf_sort_deallocate
  if(nelm_max_solve_save.ne.0) call wf_solve_deallocate
  if(node_max_field_save.ne.0) call wf_field_deallocate
  if(ngxmax_save.ne.0) call wf_win_deallocate
  return
end subroutine wf_menu

!     ***** DEBUG INFORMATION ROUTINE *****

subroutine WFINFO
  
  use wfcomm
  USE wfindex
  USE libchar
  implicit none
  integer   :: IE,IN,NN,NE,IA,IS,NSD
  real(rkind)   :: R,Z
  character :: KID*1
  
8001 continue

  write(6,*) '## INPUT:  E:element  N,V:node  '//&
            ' S:side  A,M:antenna  F:FEP  X:end'
  read(5,'(A1)',ERR=8001,END=9000) KID

  call toupper(KID)
  
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
        NN=node_nside_nelm(IN,IE)
        write(6,'(I5,1P2E12.4,1X,I3)') NN,xnode(NN),ynode(NN),KANOD(NN)
     end do
     goto 8002
     
  elseif(KID.eq.'N') then
8003 write(6,*) '## INPUT: Node number '
     read(5,*,ERR=8003,END=8001) NN
     if(NN.eq.0) goto 8001
     write(6,'(A,1P2E12.4,I5)') '   xnode,ynode,KA =',&
          xnode(NN),ynode(NN),KANOD(NN)
     write(6,*) '---------------------------'
     write(6,*) '  Sides including this node'
     write(6,*) '---------------------------'
     do NSD=1,nseg_max
        do IN=1,2
           if(node_nseg(IN,NSD).eq.NN) then
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
     if(NE.gt.nelm_max) NE=nelm_max
     call wf_fep(R,Z,NE)
     write(6,*) '   NE =',NE
     if(NE.ne.0) then
        write(6,*) '  NODE=',&
             node_nside_nelm(1,NE),node_nside_nelm(2,NE),node_nside_nelm(3,NE)
     end if
     goto 8006
     
  elseif(KID.eq.'V') then
     do NN=1,node_max
        if(KANOD(NN).gt.0) then
           write(6,'(I5,1P3E12.4)')&
                NN,xnode(NN),ynode(NN)
        end if
     end do

  elseif(KID.eq.'S') then
     do NSD=1,nseg_max
        write(6,*) NSD,node_nseg(1,NSD),node_nseg(2,NSD),KASID(NSD)
     end do
     
  elseif(KID.eq.'X') then
     goto 9000
  end if
  goto 8001

9000 return
end subroutine WFINFO

END MODULE wfmenu
