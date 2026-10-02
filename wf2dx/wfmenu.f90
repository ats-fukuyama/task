! wfmenu.f90

MODULE wfmenu

  PRIVATE
  PUBLIC wf_menu

CONTAINS

  SUBROUTINE wf_menu

  USE libmpi
  USE wfcomm
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
  IMPLICIT NONE
  
  INTEGER  :: MODE
  INTEGER  :: IERR
  CHARACTER(LEN=1):: KID
  CHARACTER(LEN=80):: LINE
  INTEGER:: IDEBUG_SAVE

1 CONTINUE

  IF(nrank.EQ.0) THEN
     WRITE(6,*) '## INPUT: P,V:parm  D:div  A:ant', &
                         ' R:run C:check G:graph  S,L:file  Q:QUIT'
     CALL TASK_KLIN(LINE,KID,MODE,WF_PARM)
  END IF
  CALL MTX_BARRIER

  CALL mtx_broadcast1_character(KID)
  CALL mtx_broadcast1_integer(MODE)
  IF(MODE.NE.1) GOTO 1

  IF(KID.EQ.'P') THEN ! namelist input
     IF(nrank.EQ.0) CALL wf_parm(0,'wf',IERR)
     CALL pl_broadcast
     CALL dp_broadcast
     CALL wfparm_broadcast
     GOTO 1
  ELSEIF (KID.EQ.'V') THEN
     IF (nrank.EQ.0) CALL WF_VIEW
  ELSEIF (KID.EQ.'D') THEN
     CALL wf_div
  ELSEIF (KID.EQ.'A') THEN
     CALL wf_ant
  ELSEIF (KID.EQ.'C') THEN
     CALL wf_wpre(IERR)
  ELSEIF (KID.EQ.'R') THEN
     CALL wf_wave
  ELSEIF (KID.EQ.'G') THEN
     IF (nrank.EQ.0) CALL wf_gout
  ELSEIF (KID.EQ.'L') THEN
     IDEBUG_SAVE=IDEBUG
     IDEBUG=1
     CALL pl_load(ierr)
     CALL wf_load_wg(ierr)
     IDEBUG=IDEBUG_SAVE
  ELSEIF (KID.EQ.'?') THEN
     if(nrank.EQ.0) CALL WFINFO
  ELSEIF (KID.EQ.'Q') THEN
     GOTO 9000
  end IF
  KID=' '
  goto 1

9000 CONTINUE
  IF(wf_solve_allocated) call wf_solve_deallocate
  IF(wf_win_allocated) call wf_win_deallocate
  IF(fem_mesh_allocated) call fem_mesh_deallocate
  IF(fem_base_allocated) call fem_base_deallocate

  return
end subroutine wf_menu

!     ***** DEBUG INFORMATION ROUTINE *****

SUBROUTINE WFINFO
  
  USE wfcomm
  USE wfindex
  USE libchar
  IMPLICIT NONE
  INTEGER:: nelm,nside,node,nant,np,np0,nseg
  REAL(rkind):: x,y
  CHARACTER:: KID*1
  
8001 CONTINUE

  WRITE(6,*) '## INPUT:  E:element  N,V:node  '//&
            ' S:side  A,M:antenna  F:FEP  X:end'
  READ(5,'(A1)',ERR=8001,END=9000) KID

  CALL TOUPPER(KID)
  
  IF(KID.EQ.'E') THEN
8002 CONTINUE
     WRITE(6,*) '## INPUT: Element number '
     READ(5,*,ERR=8002,END=8001) nelm
     IF(nelm.EQ.0) GOTO 8001
     WRITE(6,'(A,5I8)') '  nelm,nmed,nelm1/2/3=', &
          nelm,nmed_nelm(nelm),nelm1_nside_nelm(1,nelm), &
          nelm1_nside_nelm(2,nelm), nelm1_nside_nelm(3,nelm)
     WRITE(6,'(A)')     '------------------------------------'
     WRITE(6,'(A)') '  node           X           Y  mode'
     WRITE(6,'(A)')     '------------------------------------'
     DO nside=1,3
        node=node_nside_nelm(nside,nelm)
        WRITE(6,'(I5,1P2E12.4,1X,I6)') &
             node,xnode(node),ynode(node),mode_node(node)
     END DO
     GOTO 8002
     
  ELSEIF(KID.EQ.'N') THEN
8003 CONTINUE
     WRITE(6,*) '## INPUT: node number '
     READ(5,*,ERR=8003,END=8001) node
     IF(node.eq.0) GOTO 8001
     WRITE(6,'(A,1P2E12.4,I5)') '   xnode,ynode,mode =',&
          xnode(node),ynode(node),mode_node(node)
     WRITE(6,*) '---------------------------'
     WRITE(6,*) '  nseg including this node'
     WRITE(6,*) '---------------------------'
     DO nseg=1,nseg_max
        DO nside=1,2
           IF(node_nseg(nside,nseg).EQ.node) THEN
              WRITE(6,'(A,2I8)') '   nseg,nside =',&
                   nseg,nside
           END IF
        END DO
     END DO
     GOTO 8003
     
  ELSEIF(KID.eq.'A') then
8004 CONTINUE
     WRITE(6,*) '## INPUT: Antenna number, point number (guide) '
     READ(5,*,ERR=8004,END=8001) nant,np0
     IF(nant.EQ.0) GOTO 8001
     WRITE(6,*) '   XJ0,YJ0 =',&
          x_np0_nant(np0,nant),y_np0_nant(np0,nant)
     GOTO 8004
     
  ELSEIF(KID.EQ.'M') THEN
8005 CONTINUE
     WRITE(6,*) '## INPUT: Antenna number, point number (mesh)'
     READ(5,*,ERR=8005,END=8001) nant,np
     IF(nant.eq.0) goto 8001
     write(6,*) '   nelm,X,Y =',&
          nelm_np_nant(np,nant),x_np_nant(np,nant),y_np_nant(np,nant)
     goto 8005
     
  ELSEIF(KID.EQ.'F') THEN
8006 CONTINUE
     WRITE(6,*) '## INPUT: NE,X,Y'
     READ(5,*,ERR=8006,END=8001) nelm,X,Y
     IF(nelm.LE.0) GOTO 8001
     IF(nelm.GT.nelm_max) nelm=nelm_max
     CALL wf_fep(X,Y,nelm)
     WRITE(6,*) '   nelm =',nelm
     IF(nant.NE.0) THEN
        WRITE(6,*) '  NODE=',&
             node_nside_nelm(1,nelm),node_nside_nelm(2,nelm), &
             node_nside_nelm(3,nelm)
     END IF
     GOTO 8006
     
  ELSEIF(KID.EQ.'V') THEN
     DO node=1,node_max
        IF(mode_node(node).NE.0) then
           WRITE(6,'(I5,1P3E12.4)')&
                node,xnode(node),ynode(node)
        END IF
     END DO

  ELSEIF(KID.EQ.'S') THEN
     DO nseg=1,nseg_max
        WRITE(6,*) nseg,node_nseg(1,nseg),node_nseg(2,nseg),mode_nseg(nseg)
     END DO
     
  ELSEIF(KID.eq.'X') THEN
     GOTO 9000
  END IF
  GOTO 8001

9000 RETURN
END SUBROUTINE WFINFO

END MODULE wfmenu
