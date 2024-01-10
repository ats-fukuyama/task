! wfant.f90

MODULE wfant

  PRIVATE
  PUBLIC wf_ant
  PUBLIC wf_modant

CONTAINS
  
  SUBROUTINE wf_ant

  use libmpi
  USE libchar
  use wfcomm
  use wfparm
  USE wfgsub
  implicit none
  integer   :: IERR,nant,npoint
  character :: KID*1

  if (nrank.eq.0) WRITE(6,*) '--- wf_set_bdy start ---'
!  CALL wf_set_bdy(IERR)
  CALL wf_wpre(IERR)
  IF(IERR.NE.0) RETURN

1 continue
  if (nrank.eq.0) then
     WRITE(6,601)
601  FORMAT(' ','## INPUT: A/ANT  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',&
                          'W/LIST  X/EXIT')
     READ(5,'(A1)',ERR=1,END=9000) KID
     CALL toupper(KID)
  end if
  call mtx_barrier
  call mtx_broadcast_character(KID,1)
  
  IF    (KID.EQ.'A') THEN
     CALL wf_define_antenna

  ELSEIF(KID.EQ.'G') THEN
     if (nrank.eq.0) CALL wf_gr_antenna

  ELSEIF(KID.EQ.'P') THEN
     if (nrank.eq.0) CALL wf_parm(0,'wf',IERR)
     call wfparm_broadcast

  ELSEIF(KID.EQ.'V') THEN
     if (nrank.eq.0) CALL WF_VIEW

  ELSEIF(KID.EQ.'S') THEN
!test     if (nrank.eq.0) CALL WFWANT

  ELSEIF(KID.EQ.'L') THEN
!test     CALL WFRANT

  ELSEIF(KID.EQ.'W') THEN
     if (nrank.eq.0) then
        DO nant=1,nant_max
           WRITE(6,610) (npoint,RJ0(npoint,nant),ZJ0(npoint,nant),&
                npoint=1,JNUM0(nant))
610        FORMAT(' ',I5,2F12.5)
        ENDDO
        DO nant=1,nant_max
           WRITE(6,611) (npoint,RJ(npoint,nant),ZJ(npoint,nant), &
                JELMT(npoint,nant),npoint=1,JNUM(nant))
611        FORMAT(' ',I5,2F12.5,I10)
        ENDDO
     end if
  ELSEIF(KID.EQ.'X') THEN
     GOTO 9000
  ENDIF
  GOTO 1
  
9000 RETURN
END SUBROUTINE wf_ant

!     ****** define Antenna Data ******

SUBROUTINE wf_define_antenna

  use libmpi
  USE libchar
  use wfcomm
  implicit none
  integer   :: nant,npoint,IERR
  real(rkind)   :: DEGN,DTHETA,THETA,R,Z
  character KID*1

  DEGN=PI/180.D0

1 continue

  if (nrank.eq.0) then

     WRITE(6,*) '## ENTER NUMBER OF ANTENNAS (1-',NAM,')'
     READ(5,*,ERR=1,END=9000) nant_max
     IF(nant_max.EQ.0) goto 9000
     IF(nant_max.LT.1.OR.nant_max.GT.NAM) GOTO 1
     
     DO nant=1,nant_max
        JNUM0(nant)=0
        
2       continue
        
        WRITE(6,*) '## ANTENNA NUMBER = ',nant
        WRITE(6,*) '## TYPE: C/CIRCLE  A/ARC  P/POINTS'//&
             &            '  X/EXIT'
        READ(5,'(A1)',ERR=2,END=1) KID
        CALL toupper(KID)
        
        IF(KID.EQ.'C') THEN              ! Circular
3          WRITE(6,602) RD,npoint_max_nant(nant)
602        FORMAT(' ','## RD,npoint_max_nant(nant) = ',F10.3,I5)
           READ(5,*,ERR=3,END=2) RD(nant),npoint_max_nant(nant)
           DTHETA=2.D0*PI/(npoint_max_nant(nant)-1)
           DO npoint=1,npoint_max_nant(nant)
              THETA=DTHETA*(npoint-1)+1.D-5
              RJ0(npoint,nant)=RD(nant)*COS(THETA)
              ZJ0(npoint,nant)=RD(nant)*SIN(THETA)
           END DO
           JNUM0(nant)=npoint_max_nant(nant)
           
        ELSEIF(KID.EQ.'A') THEN           ! Arc
4          WRITE(6,603) &
                THETJ1(nant),THETJ2(nant),RD(nant),npoint_max_nant(nant)
603        FORMAT(' ','## THETJ1,THETJ2 = ',2F10.3/&
                  ' ','## RD,npoint_max_nant = ',F10.3,I5)
        READ(5,*,ERR=4,END=2) &
             THETJ1(nant),THETJ2(nant),RD(nant),npoint_max_nant(nant)

           THETA=DEGN*THETJ1(nant)
           IF(itype_mesh.EQ.1) THEN  ! rect
              RJ0(1,nant)=RD(nant)*COS(THETA)+2.0d0*RD(nant)
              ZJ0(1,nant)=RD(nant)*SIN(THETA)             
           ELSEIF(itype_mesh.EQ.2) THEN ! circle
              RJ0(1,nant)=1.5D0*RD(nant)*COS(THETA)
              ZJ0(1,nant)=1.5D0*RD(nant)*SIN(THETA)
           END IF
           DTHETA=(THETJ2(nant)-THETJ1(nant))/(npoint_max_nant(nant)-3)
           DO npoint=2,npoint_max_nant(nant)-1
              THETA=DEGN*(DTHETA*(npoint-2)+THETJ1(nant))
              RJ0(npoint,nant)=RD(nant)*COS(THETA)
              ZJ0(npoint,nant)=RD(nant)*SIN(THETA)
           END DO
           THETA=DEGN*THETJ2(nant)
           if(itype_mesh.eq.1) then
              RJ0(npoint_max_nant(nant),nant)=RD(nant)*COS(THETA)+2.0d0*RD(nant)
              ZJ0(npoint_max_nant(nant),nant)=RD(nant)*SIN(THETA)
           elseif(itype_mesh.eq.2) then
              RJ0(npoint_max_nant(nant),nant)=1.5D0*RD(nant)*COS(THETA)
              ZJ0(npoint_max_nant(nant),nant)=1.5D0*RD(nant)*SIN(THETA)
           end if
           JNUM0(nant)=npoint_max_nant(nant)
           
        ELSEIF(KID.EQ.'P') THEN    ! given points 
6          WRITE(6,*) '## NUMBER OF POINTS : npoint_max_nant(nant)=',npoint_max_nant(nant)
           READ(5,*,ERR=6,END=2) npoint_max_nant(nant)
           DO npoint=1,npoint_max_nant(nant)
7             WRITE(6,*) '## NO.',npoint,': R,Z ?'
              READ(5,*,ERR=7,END=6) R,Z
              RJ0(npoint,nant)=R
              ZJ0(npoint,nant)=Z
           END DO
           JNUM0(nant)=npoint_max_nant(nant)
           
        ELSEIF(KID.EQ.'X') THEN
           GOTO 1
        ELSE
           WRITE(6,*) 'XX UNKNOWN KID: ',KID
           GOTO 2
        ENDIF
     END DO
  end if

  SELECT CASE(MODELG)
  CASE(2)
     do nant=1,nant_max
        do npoint=1,npoint_max_nant(nant)
           RJ0(npoint,nant)=RJ0(npoint,nant)+RR
        end do
     end do
  END SELECT

  call mtx_barrier
  call wf_broadcast_antenna

  CALL wf_modant(IERR)
  IF(IERR.NE.0) GOTO 9000
  
9000 RETURN
END SUBROUTINE wf_define_antenna

! --- broadcast data ---

subroutine wf_broadcast_antenna
  
  use wfcomm
  use libmpi
  implicit none
  
  integer :: nant,npoint
  real(rkind),dimension(NJM)::ddatar,ddataz

  call mtx_broadcast1_integer(nant_max)

  call mtx_broadcast_integer(JNUM0,nant_max)
  
  do nant=1,nant_max

     if(nrank.eq.0) then
        npoint_max_nant(nant)=JNUM0(nant)
        do npoint=1,npoint_max_nant(nant)
           ddatar(npoint)=RJ0(npoint,nant)
           ddataz(npoint)=ZJ0(npoint,nant)
        end do
     end if
     
     call mtx_broadcast1_integer(npoint_max_nant(nant))
     call mtx_broadcast_real8(ddatar,npoint_max_nant(nant))
     call mtx_broadcast_real8(ddataz,npoint_max_nant(nant))

     JNUM0(nant)=npoint_max_nant(nant)
     do npoint=1,npoint_max_nant(nant)
        RJ0(npoint,nant)=ddatar(npoint)
        ZJ0(npoint,nant)=ddataz(npoint)
     end do

  end do

  return
end subroutine wf_broadcast_antenna

!     ****** MODIFY ANTENNA DATA ******

SUBROUTINE wf_modant(IERR)
  
  use wfcomm
  USE wfindex
  USE wfsub,ONLY: wf_cross
  implicit none

  integer,intent(out) :: IERR
  integer :: NE,nant,L,KN,LS,N,ID,NENEXT,NENEW,NE2
  INTEGER:: NSD
  real(rkind) :: RC,ZC

  NE=0
  DO nant=1,nant_max
     CALL wf_fep(RJ0(1,nant),ZJ0(1,nant),NE)
     IF(nrank.EQ.0) &
          WRITE(6,'(A,I5,1P2E12.4,I5)') &
          'nant,RJ0,ZJ0=',nant,RJ0(1,nant),ZJ0(1,nant),NE

     !    outside starting point

     IF(NE.EQ.0) THEN
        WRITE(6,'(A,I6)') 'JNUM0=',JNUM0(nant)
        IF(JNUM0(nant).EQ.1) GOTO 8500 ! error: one point and outside
        NE2=NE
        CALL wf_fep(RJ0(2,nant),ZJ0(2,nant),NE2)
        IF(nrank.EQ.0) &
             WRITE(6,'(A,I5,1P2E12.4,I5)') &
             'nant,RJ0,ZJ0=',nant,RJ0(2,nant),ZJ0(2,nant),NE2
        DO NSD=1,nseg_max              ! look for crossing boundary
           L =INSID(NSD)
           NE=NESID(NSD)
           KN=KNELM(L,NE)

           IF(KN.eq.0) THEN
              NE2=NE
              CALL wf_fep(RJ0(2,nant),ZJ0(2,nant),NE2)
!              IF(nrank.EQ.0) &
!                   WRITE(6,'(A,I5,1P2E12.4,I5)') &
!                   'NA,RJ0,ZJ0=',NA,RJ0(2,NA),ZJ0(2,NA),NE2
              CALL wf_cross(RJ0(1,nant),ZJ0(1,nant),&
                   &    RJ0(2,nant),ZJ0(2,nant),&
                   &    NE,L,RC,ZC,IERR)
              IF(IERR.EQ.0) THEN
                 LS=L
                 GOTO 1000   ! crossing point on boundary found
              ENDIF
           ENDIF

        ENDDO
        GOTO 8000 ! error: no crossing boundary found
        
1000    CONTINUE  ! set starting point (N=1)
        N=1
        RJ(N,nant)=RC
        ZJ(N,nant)=ZC
        JELMT(N,nant)=NE
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'nant,N,NE,R,Z=' ,&
                nant,N,NE,RJ(N,nant),ZJ(N,nant)
        ENDIF

!    inside starting point

     ELSE
        N=1
        RJ(N,nant)=RJ0(1,nant)
        ZJ(N,nant)=ZJ0(1,nant)
        JELMT(N,nant)=NE
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'nant,N,NE,R,Z=',&
                nant,N,NE,RJ(N,nant),ZJ(N,nant)
        ENDIF
        LS=0
     ENDIF
     
     DO ID=2,JNUM0(nant)
        NENEXT=NE
        CALL wf_fep(RJ0(ID,nant),ZJ0(ID,nant),NENEXT)
3000    CONTINUE
        IF(NENEXT.EQ.NE) GOTO 4500 ! next antenna node in the same element
        DO L=1,3                   ! look for crossing point
           IF(L.NE.LS) THEN
              CALL wf_cross(RJ (N ,nant),ZJ (N ,nant),&
                   &    RJ0(ID,nant),ZJ0(ID,nant),&
                   &         NE,L,RC,ZC,IERR)
              IF(IERR.EQ.0) THEN
                 LS=L
                 GOTO 4000 ! crossing point found
              ENDIF
           ENDIF
        ENDDO
        GOTO 8100 
        
4000    IF(N+1.GT.NJM) GOTO 8200 ! error: number of antenna nodes overflow 
        N=N+1
        RJ(N,nant)=RC
        ZJ(N,nant)=ZC
        JELMT(N,nant)=NE
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'nant,N,NE,R,Z=',&
                nant,N,NE,RJ(N,nant),ZJ(N,nant)
        ENDIF
        NENEW=KNELM(LS,NE) ! look for neighboring element
        IF(NENEW.GT.0) THEN ! element found
           DO L=1,3
              IF(KNELM(L,NENEW).EQ.NE) LS=L
           ENDDO
           NE=NENEW         ! set side LS and element NE
        ELSE
           IF(ID.EQ.JNUM0(nant).AND.NENEXT.LE.0) GOTO 6000 ! last point outside
           GOTO 8400 ! error: cross point to outside
        ENDIF
        GOTO 3000 ! look for next cross point 
        
4500    IF(N+1.GT.NJM) GOTO 8200 ! error: number of antenna nodes overflow 
        N=N+1
        RJ(N,nant)=RJ0(ID,nant)
        ZJ(N,nant)=ZJ0(ID,nant)
        JELMT(N,nant)=NE
        LS=0
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'nant,N,NE,R,Z=',&
                nant,N,NE,RJ(N,nant),ZJ(N,nant)
        ENDIF
     ENDDO
     
6000 JNUM(nant)=N
  ENDDO
  IERR=0
  RETURN
  
8000 IERR=8000
  JNUM(nant)=N
  if(nrank.eq.0) WRITE(6,800) IERR,nant,N
800 FORMAT(' ','## wf_modant ERROR : IERR = ',I5/&
         & ' ','                   : CANNOT FIND BOUNDARY POINT'/&
         & ' ','                   : nant,N=',2I7)
  JNUM(nant)=N
  RETURN
  
8100 IERR=8100
  if(nrank.eq.0) WRITE(6,810) IERR,nant,ID,N
810 FORMAT(' ','## wf_modant ERROR : IERR = ',I5/&
         & ' ','                   : CANNOT FIND NEXT CROSSPOINT'/&
         & ' ','                   : nant,ID,N =',3I7)
  JNUM(nant)=N
  RETURN
  
8200 if(nrank.eq.0) WRITE(6,820) nant,ID,N,NJM
820 FORMAT(' ','## MODANT ERROR : N.GT.NJM '/&
         & ' ','                : nant,ID,N,NJM = ',4I7)
  IERR=8200
  JNUM(nant)=N
  RETURN
  
8400 IERR=8400
  if(nrank.eq.0) WRITE(6,840) IERR,nant,ID,NE,N
840 FORMAT(' ','## MODANT ERROR : IERR =',I5/&
         & ' ','                : ABMORMAL END OF ANTENNA DATA '/&
         & ' ','                : nant,ID,NE,N = ',4I7)
  JNUM(nant)=N
  IERR=8400
  RETURN
  
8500 IERR=8500
  if(nrank.eq.0) WRITE(6,850) IERR,nant,NE,JNUM0(nant)
850 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         & ' ','                : nant,NE,JNUM0 = ',3I7)
  RETURN
  
END SUBROUTINE wf_modant

END MODULE wfant
