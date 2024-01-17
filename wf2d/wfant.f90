!     $Id: wfant.f90,v 1.5 2012/03/05 06:29:02 maruyama Exp $
!
!     ######### /TASK/WF2/WFANT ########
!
!      ANTENNA DATA GENERATION PROGRAM
!
!     #################################

SUBROUTINE WFANT

  use libmpi
  USE libchar
  use wfcomm
  use wfparm
  implicit none
  integer   :: IERR,NA,N
  character :: KID*1

  if (nrank.eq.0) WRITE(6,*) '--- SETBDY start ---'
  CALL SETBDY(IERR)
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
     CALL WFDEFA

  ELSEIF(KID.EQ.'G') THEN
     if (nrank.eq.0) CALL wf_gr_antenna

  ELSEIF(KID.EQ.'P') THEN
     if (nrank.eq.0) CALL WF_PARM(0,'wf',IERR)
     call wfparm_broadcast

  ELSEIF(KID.EQ.'V') THEN
     if (nrank.eq.0) CALL WF_VIEW

  ELSEIF(KID.EQ.'S') THEN
!test     if (nrank.eq.0) CALL WFWANT

  ELSEIF(KID.EQ.'L') THEN
!test     CALL WFRANT

  ELSEIF(KID.EQ.'W') THEN
     if (nrank.eq.0) then
        DO NA=1,NAMAX
           WRITE(6,610) (N,RJ0(N,NA),ZJ0(N,NA),&
                N=1,JNUM0(NA))
610        FORMAT(' ',I5,2F12.5)
        ENDDO
        DO NA=1,NAMAX
           WRITE(6,611) (N,RJ(N,NA),ZJ(N,NA),JELMT(N,NA),&
                N=1,JNUM(NA))
611        FORMAT(' ',I5,2F12.5,I10)
        ENDDO
     end if
  ELSEIF(KID.EQ.'X') THEN
     GOTO 9000
  ENDIF
  GOTO 1
  
9000 RETURN
END SUBROUTINE WFANT

!     ****** Define Antenna Data ******
SUBROUTINE WFDEFA

  use libmpi
  USE libchar
  use wfcomm
  implicit none
  integer   :: NA,NJ,IERR
  real(rkind)   :: DEGN,DTHETA,THETA,R,Z
  character KID*1

  DEGN=PI/180.D0

1 continue

  if (nrank.eq.0) then

     WRITE(6,*) '## ENTER NUMBER OF ANTENNAS (1-',NAM,')'
     READ(5,*,ERR=1,END=9000) NAMAX
     IF(NAMAX.EQ.0) goto 9000
     IF(NAMAX.LT.1.OR.NAMAX.GT.NAM) GOTO 1
     
     DO NA=1,NAMAX
        JNUM0(NA)=0
        
2       continue
        
        WRITE(6,*) '## ANTENNA NUMBER = ',NA
        WRITE(6,*) '## TYPE: C/CIRCLE  A/ARC  P/POINTS'//&
             &            '  X/EXIT'
        READ(5,'(A1)',ERR=2,END=1) KID
        CALL toupper(KID)
        
        IF(KID.EQ.'C') THEN
3          WRITE(6,602) RD,NJMAX
602        FORMAT(' ','## RD,NJMAX = ',F10.3,I5)
           READ(5,*,ERR=3,END=2) RD,NJMAX
           DTHETA=2.D0*PI/(NJMAX-1)
           DO NJ=1,NJMAX
              THETA=DTHETA*(NJ-1)+1.D-5
              RJ0(NJ,NA)=RD*COS(THETA)
              ZJ0(NJ,NA)=RD*SIN(THETA)
           END DO
           JNUM0(NA)=NJMAX
           
        ELSEIF(KID.EQ.'A') THEN
4          WRITE(6,603) THETJ1,THETJ2,RD,NJMAX
603        FORMAT(' ','## THETJ1,THETJ2 = ',2F10.3/&
                  ' ','## RD,NJMAX = ',F10.3,I5)
           READ(5,*,ERR=4,END=2) THETJ1,THETJ2,RD,NJMAX

           THETA=DEGN*THETJ1
           if(iddiv.eq.1) then
              RJ0(1,NA)=RD*COS(THETA)+2.0d0*RD
              ZJ0(1,NA)=RD*SIN(THETA)             
           elseif(iddiv.eq.2) then
              RJ0(1,NA)=1.5D0*RD*COS(THETA)
              ZJ0(1,NA)=1.5D0*RD*SIN(THETA)
           end if
           DTHETA=(THETJ2-THETJ1)/(NJMAX-3)
           DO NJ=2,NJMAX-1
              THETA=DEGN*(DTHETA*(NJ-2)+THETJ1)
              RJ0(NJ,NA)=RD*COS(THETA)
              ZJ0(NJ,NA)=RD*SIN(THETA)
           END DO
           THETA=DEGN*THETJ2
           if(iddiv.eq.1) then
              RJ0(NJMAX,NA)=RD*COS(THETA)+2.0d0*RD
              ZJ0(NJMAX,NA)=RD*SIN(THETA)
           elseif(iddiv.eq.2) then
              RJ0(NJMAX,NA)=1.5D0*RD*COS(THETA)
              ZJ0(NJMAX,NA)=1.5D0*RD*SIN(THETA)
           end if
           JNUM0(NA)=NJMAX
           
        ELSEIF(KID.EQ.'P') THEN
6          WRITE(6,*) '## NUMBER OF POINTS : NJMAX=',NJMAX
           READ(5,*,ERR=6,END=2) NJMAX
           DO NJ=1,NJMAX
7             WRITE(6,*) '## NO.',NJ,': R,Z ?'
              READ(5,*,ERR=7,END=6) R,Z
              RJ0(NJ,NA)=R
              ZJ0(NJ,NA)=Z
           END DO
           JNUM0(NA)=NJMAX
           
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
     do NA=1,NAMAX
        do NJ=1,NJMAX
           RJ0(NJ,NA)=RJ0(NJ,NA)+RR
        end do
     end do
  END SELECT

  call mtx_barrier
  call wfant_broadcast

  CALL MODANT(IERR)
  IF(IERR.NE.0) GOTO 9000
  
9000 RETURN
END SUBROUTINE WFDEFA

! --- broadcast data ---
subroutine wfant_broadcast
  
  use wfcomm
  use libmpi
  implicit none
  
  integer :: NA,NJ
  real(rkind),dimension(NJM)::ddatar,ddataz

  call mtx_broadcast1_integer(NAMAX)

  call mtx_broadcast_integer(JNUM0,NAMAX)
  
  do NA=1,NAMAX

     if(nrank.eq.0) then
        NJMAX=JNUM0(NA)
        do NJ=1,NJMAX
           ddatar(NJ)=RJ0(NJ,NA)
           ddataz(NJ)=ZJ0(NJ,NA)
        end do
     end if
     
     call mtx_broadcast1_integer(NJMAX)
     call mtx_broadcast_real8(ddatar,NJMAX)
     call mtx_broadcast_real8(ddataz,NJMAX)

     JNUM0(NA)=NJMAX
     do NJ=1,NJMAX
        RJ0(NJ,NA)=ddatar(NJ)
        ZJ0(NJ,NA)=ddataz(NJ)
     end do

  end do

  return
end subroutine wfant_broadcast
