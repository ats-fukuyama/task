! wfdiv.f90

!     ######### /TASK/WF2D/WFDIV ########
!
!      MESH DATA GENERATION PROGRAM
!
!     #################################

MODULE wfdiv

  PRIVATE
  PUBLIC wf_div

CONTAINS

  SUBROUTINE wf_div

    USE wfcomm
    USE wfparm
    USE wfdiv_rect
    USE wfdiv_circle
    USE libmpi
    USE libchar
  
    IMPLICIT NONE
    INTEGER   :: NE,NN,IERR
    CHARACTER :: KID*1

1   continue

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
        write(6,'(A24)') '## TYPE: X/RECT C/CIRCLE'
        read(5,'(A1)') KID
        call toupper(KID)
        if (KID.ne."X".and.KID.ne."C") goto 2
     end if
     call mtx_barrier
     call mtx_broadcast_character(KID,1)

     IF(nrank.EQ.0) THEN
        IF(KID.EQ.'X') THEN
           CALL wf_div_rect_input
        ELSEIF(KID.EQ.'C') THEN
           CALL wf_div_circle_input
        END IF
     END IF
     CALL wfdiv_broadcast
     
     IF(KID.EQ.'X') THEN
        CALL wf_div_rect_exec
     ELSEIF(KID.EQ.'C') THEN
        CALL wf_div_circle_exec
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
     
     NBMAX=0
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
end subroutine wf_div

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

END MODULE wfdiv
