! wfant.f90

MODULE wfant

  PRIVATE
  PUBLIC wf_ant
  PUBLIC wf_modant

CONTAINS
  
  SUBROUTINE wf_ant

  use wfcomm
  use wfparm
  USE wfgsub
  USE libmpi
  USE libchar
  implicit none
  integer   :: IERR,nant,np,np0
  character :: KID*1

  if (nrank.eq.0) WRITE(6,*) '--- wf_set_bdy start ---'

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
           WRITE(6,610) (np0,x_np0_nant(np0,nant),y_np0_nant(np0,nant),&
                np0=1,np0_max_nant(nant))
610        FORMAT(' ',I5,2F12.5)
        ENDDO
        DO nant=1,nant_max
           WRITE(6,611) (np,x_np_nant(np,nant),y_np_nant(np,nant), &
                nelm_np_nant(np,nant),np=1,np_max_nant(nant))
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

     WRITE(6,*) '## ENTER NUMBER OF ANTENNAS (1-',nantm,')'
     READ(5,*,ERR=1,END=9000) nant_max
     IF(nant_max.EQ.0) goto 9000
     IF(nant_max.LT.1.OR.nant_max.GT.nantm) GOTO 1
     
     DO nant=1,nant_max
        np0_max_nant(nant)=0
        
2       continue
        
        WRITE(6,*) '## ANTENNA NUMBER = ',nant
        WRITE(6,*) '## TYPE: C/CIRCLE  A/ARC  P/POINTS'//&
             &            '  X/EXIT'
        READ(5,'(A1)',ERR=2,END=1) KID
        CALL toupper(KID)
        
        IF(KID.EQ.'C') THEN              ! Circular
3          WRITE(6,602) RD,np_max_nant(nant)
602        FORMAT(' ','## RD,np_max_nant(nant) = ',F10.3,I5)
           READ(5,*,ERR=3,END=2) RD(nant),np_max_nant(nant)
           DTHETA=2.D0*PI/(np_max_nant(nant)-1)
           DO npoint=1,np_max_nant(nant)
              THETA=DTHETA*(npoint-1)+1.D-5
              x_np0_nant(npoint,nant)=RD(nant)*COS(THETA)
              y_np0_nant(npoint,nant)=RD(nant)*SIN(THETA)
           END DO
           np0_max_nant(nant)=np_max_nant(nant)
           
        ELSEIF(KID.EQ.'A') THEN           ! Arc
4          WRITE(6,603) &
                THETJ1(nant),THETJ2(nant),RD(nant),np_max_nant(nant)
603        FORMAT(' ','## THETJ1,THETJ2 = ',2F10.3/&
                  ' ','## RD,np_max_nant = ',F10.3,I5)
        READ(5,*,ERR=4,END=2) &
             THETJ1(nant),THETJ2(nant),RD(nant),np_max_nant(nant)

           THETA=DEGN*THETJ1(nant)
           IF(mode_mesh.EQ.1) THEN  ! rect
              x_np0_nant(1,nant)=RD(nant)*COS(THETA)+2.0d0*RD(nant)
              y_np0_nant(1,nant)=RD(nant)*SIN(THETA)             
           ELSEIF(mode_mesh.EQ.2) THEN ! circle
              x_np0_nant(1,nant)=1.5D0*RD(nant)*COS(THETA)
              y_np0_nant(1,nant)=1.5D0*RD(nant)*SIN(THETA)
           END IF
           DTHETA=(THETJ2(nant)-THETJ1(nant))/(np_max_nant(nant)-3)
           DO npoint=2,np_max_nant(nant)-1
              THETA=DEGN*(DTHETA*(npoint-2)+THETJ1(nant))
              x_np0_nant(npoint,nant)=RD(nant)*COS(THETA)
              y_np0_nant(npoint,nant)=RD(nant)*SIN(THETA)
           END DO
           THETA=DEGN*THETJ2(nant)
           if(mode_mesh.eq.1) then
              x_np0_nant(np_max_nant(nant),nant)=RD(nant)*COS(THETA)+2.0d0*RD(nant)
              y_np0_nant(np_max_nant(nant),nant)=RD(nant)*SIN(THETA)
           elseif(mode_mesh.eq.2) then
              x_np0_nant(np_max_nant(nant),nant)=1.5D0*RD(nant)*COS(THETA)
              y_np0_nant(np_max_nant(nant),nant)=1.5D0*RD(nant)*SIN(THETA)
           end if
           np0_max_nant(nant)=np_max_nant(nant)
           
        ELSEIF(KID.EQ.'P') THEN    ! given points 
6          WRITE(6,*) '## NUMBER OF POINTS : np_max_nant(nant)=',np_max_nant(nant)
           READ(5,*,ERR=6,END=2) np_max_nant(nant)
           DO npoint=1,np_max_nant(nant)
7             WRITE(6,*) '## NO.',npoint,': R,Z ?'
              READ(5,*,ERR=7,END=6) R,Z
              x_np0_nant(npoint,nant)=R
              y_np0_nant(npoint,nant)=Z
           END DO
           np0_max_nant(nant)=np_max_nant(nant)
           
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
        do npoint=1,np_max_nant(nant)
           x_np0_nant(npoint,nant)=x_np0_nant(npoint,nant)+RR
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
  real(rkind),dimension(npointm)::ddatar,ddataz

  call mtx_broadcast1_integer(nant_max)

  call mtx_broadcast_integer(np0_max_nant,nant_max)
  
  do nant=1,nant_max

     if(nrank.eq.0) then
        np_max_nant(nant)=np0_max_nant(nant)
        do npoint=1,np_max_nant(nant)
           ddatar(npoint)=x_np0_nant(npoint,nant)
           ddataz(npoint)=y_np0_nant(npoint,nant)
        end do
     end if
     
     call mtx_broadcast1_integer(np_max_nant(nant))
     call mtx_broadcast_real8(ddatar,np_max_nant(nant))
     call mtx_broadcast_real8(ddataz,np_max_nant(nant))

     np0_max_nant(nant)=np_max_nant(nant)
     do npoint=1,np_max_nant(nant)
        x_np0_nant(npoint,nant)=ddatar(npoint)
        y_np0_nant(npoint,nant)=ddataz(npoint)
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
  integer :: nelm,nant,nside,nelm1,nsides,np,np0,nelm_next,nelm_new,nelm2
  INTEGER:: nseg
  real(rkind) :: xc,yc

  nelm=0
  DO nant=1,nant_max
     CALL wf_fep(x_np0_nant(1,nant),y_np0_nant(1,nant),nelm)
     IF(nrank.EQ.0) &
          WRITE(6,'(A,I5,1P2E12.4,I5)') &
          'nant,x_np0_nant,y_np0_nant=', &
           nant,x_np0_nant(1,nant),y_np0_nant(1,nant),nelm

     !    outside starting point

     IF(nelm.EQ.0) THEN
        WRITE(6,'(A,I6)') 'np0_max_nant=',np0_max_nant(nant)
        IF(np0_max_nant(nant).EQ.1) GOTO 8500 ! error: one point and outside
        nelm2=nelm
        CALL wf_fep(x_np0_nant(2,nant),y_np0_nant(2,nant),nelm2)
        IF(nrank.EQ.0) &
             WRITE(6,'(A,I5,1P2E12.4,I5)') &
             'nant,x_np0_nant,y_np0_nant=', &
              nant,x_np0_nant(2,nant),y_np0_nant(2,nant),nelm2
        DO nseg=1,nseg_max              ! look for crossing boundary
           nside=nside_nseg(1,nseg)
           nelm=nelm_nseg(1,nseg)
           nelm1=nelm1_nside_nelm(nside,nelm)

           IF(nelm1.eq.0) THEN
              nelm2=nelm
              CALL wf_fep(x_np0_nant(2,nant),y_np0_nant(2,nant),nelm2)
!              IF(nrank.EQ.0) &
!                   WRITE(6,'(A,I5,1P2E12.4,I5)') &
!                   'NA,x_np0_nant,y_np0_nant=',NA,x_np0_nant(2,NA),y_np0_nant(2,NA),NE2
              CALL wf_cross(x_np0_nant(1,nant),y_np0_nant(1,nant), &
                            x_np0_nant(2,nant),y_np0_nant(2,nant), &
                            nelm,nside,xc,yc,ierr)
              IF(IERR.EQ.0) THEN
                 nsides=nside
                 GOTO 1000   ! crossing point on boundary found
              ENDIF
           ENDIF

        ENDDO
        GOTO 8000 ! error: no crossing boundary found
        
1000    CONTINUE  ! set starting point (np=1)
        np=1
        x_np_nant(np,nant)=xc
        y_np_nant(np,nant)=yc
        nelm_np_nant(np,nant)=nelm
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'nant,N,NE,R,Z=' ,&
                nant,np,nelm,x_np_nant(np,nant),y_np_nant(np,nant)
        ENDIF

!    inside starting point

     ELSE
        np=1
        x_np_nant(np,nant)=x_np0_nant(1,nant)
        y_np_nant(np,nant)=y_np0_nant(1,nant)
        nelm_np_nant(np,nant)=nelm
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') &
                'nant,np,nelm,x,y=',&
                nant,np,nelm,x_np_nant(nelm,nant),y_np_nant(nelm,nant)
        ENDIF
        nsides=0
     ENDIF
     
     DO np0=2,np0_max_nant(nant)
        nelm_next=nelm
        CALL wf_fep(x_np0_nant(np0,nant),y_np0_nant(np0,nant),nelm_next)
3000    CONTINUE
        IF(nelm_next.EQ.nelm) GOTO 4500 ! next antenna node in the same element
        DO nside=1,3                   ! look for crossing point
           IF(nside.NE.nsides) THEN
              CALL wf_cross(x_np_nant (np,nant),y_np_nant (np,nant),&
                            x_np0_nant(np0,nant),y_np0_nant(np0,nant),&
                            nelm,nside,xc,yc,IERR)
              IF(IERR.EQ.0) THEN
                 nsides=nside
                 GOTO 4000 ! crossing point found
              ENDIF
           ENDIF
        ENDDO
        GOTO 8100 
        
4000    IF(np+1.GT.npointm) GOTO 8200 ! error: number of antenna pnt overflow 
        np=np+1
        x_np_nant(np,nant)=xc
        y_np_nant(np,nant)=yc
        nelm_np_nant(np,nant)=nelm
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') &
                'nant,np,nelm,x,z=', &
                nant,np,nelm,x_np_nant(np,nant),y_np_nant(np,nant)
        ENDIF
        nelm_new=nelm1_nside_nelm(nsides,nelm) ! look for neighboring element
        IF(nelm_new.GT.0) THEN ! element found
           DO nside=1,3
              IF(nelm1_nside_nelm(nside,nelm_new).EQ.nelm) nsides=nside
           ENDDO
           nelm=nelm_new         ! set side LS and element NE
        ELSE
           IF(np0.EQ.np0_max_nant(nant).AND. &
              nelm_next.LE.0) GOTO 6000 ! last point outside
           GOTO 8400 ! error: cross point to outside
        ENDIF
        GOTO 3000 ! look for next cross point 
        
4500    IF(np+1.GT.npointm) GOTO 8200 ! error: number of antenna point over 
        np=np+1
        x_np_nant(np,nant)=x_np0_nant(np0,nant)
        y_np_nant(np,nant)=y_np0_nant(np0,nant)
        nelm_np_nant(np,nant)=nelm
        nsides=0
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') &
                'nant,np,nelm,x,y=',&
                nant,np,nelm,x_np_nant(np,nant),y_np_nant(np,nant)
        ENDIF
     ENDDO
     
6000 np_max_nant(nant)=np
  ENDDO
  IERR=0
  RETURN
  
8000 IERR=8000
  np_max_nant(nant)=np
  if(nrank.eq.0) WRITE(6,800) IERR,nant,np
800 FORMAT(' ','## wf_modant ERROR : IERR = ',I5/&
           ' ','                   : CANNOT FIND BOUNDARY POINT'/&
           ' ','                   : nant,np=',2I7)
  np_max_nant(nant)=np
  RETURN
  
8100 IERR=8100
  if(nrank.eq.0) WRITE(6,810) IERR,nant,np0,np
810 FORMAT(' ','## wf_modant ERROR : IERR = ',I5/&
           ' ','                   : CANNOT FIND NEXT CROSSPOINT'/&
           ' ','                   : nant,np0,np =',3I7)
  np_max_nant(nant)=np
  RETURN
  
8200 if(nrank.eq.0) WRITE(6,820) nant,np0,np,npointm
820 FORMAT(' ','## MODANT ERROR : np.GT.npointm '/&
           ' ','                : nant,np0,np,npointm = ',4I7)
  IERR=8200
  np_max_nant(nant)=np
  RETURN
  
8400 IERR=8400
  if(nrank.eq.0) WRITE(6,840) IERR,nant,np0,nelm,np
840 FORMAT(' ','## MODANT ERROR : IERR =',I5/&
           ' ','                : ABMORMAL END OF ANTENNA DATA '/&
           ' ','                : nant,np0,nelm,np = ',4I7)
  np_max_nant(nant)=np
  IERR=8400
  RETURN
  
8500 IERR=8500
  if(nrank.eq.0) WRITE(6,850) IERR,nant,nelm,np0_max_nant(nant)
850 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
           ' ','                : nant,nelm,np0_max_nant = ',3I7)
  RETURN
  
END SUBROUTINE wf_modant

END MODULE wfant
