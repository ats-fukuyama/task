!     $Id$
!
!     ######### /TASK/WF/WFANT ########
!
!      ANTENNA DATA GENERATION PROGRAM
!
!     #################################

SUBROUTINE WFANT

  use libmpi
  use libmtx
  use wfcomm
  implicit none
  integer   :: IERR,NA,N
  character :: KID*1

  if (nrank.eq.0) WRITE(6,*) '--- SETBDY start ---'
  CALL SETBDY(IERR)
  IF(IERR.NE.0) RETURN

  if (nrank.eq.0) WRITE(6,*) '--- SETSID start ---'
  CALL SETSID(IERR)
  IF(IERR.NE.0) RETURN
  
1 continue
  if (nrank.eq.0) then
     WRITE(6,601)
601  FORMAT(' ','## INPUT: A/ANT  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',&
                          'W/LIST  X/EXIT')
     READ(5,'(A1)',ERR=1,END=9000) KID
     CALL GUCPTL(KID)
  end if
  call mtx_barrier
  call mtx_broadcast_character(KID,1)
  
  IF    (KID.EQ.'A') THEN
     CALL WFDEFA

  ELSEIF(KID.EQ.'G') THEN
     if (nrank.eq.0) CALL WFPLTA

  ELSEIF(KID.EQ.'P') THEN
     if (nrank.eq.0) CALL WFPARM(KID)
     call wfparm_broadcast

  ELSEIF(KID.EQ.'V') THEN
     if (nrank.eq.0) CALL WFVIEW

  ELSEIF(KID.EQ.'S') THEN
     if (nrank.eq.0) CALL WFWANT

  ELSEIF(KID.EQ.'L') THEN
     CALL WFRANT

  ELSEIF(KID.EQ.'W') THEN
     if (nrank.eq.0) then
        DO NA=1,NAMAX
           WRITE(6,610) (N,XJ0(N,NA),YJ0(N,NA),ZJ0(N,NA),&
                N=1,JNUM0(NA))
610        FORMAT(' ',I5,3F12.5)
        ENDDO
        DO NA=1,NAMAX
           WRITE(6,611) (N,XJ(N,NA),YJ(N,NA),ZJ(N,NA),JELMT(N,NA),&
                N=1,JNUM(NA))
611        FORMAT(' ',I5,3F12.5,I10)
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
  use libmtx
  use wfcomm
  implicit none
  integer   :: NA,NJ,IERR
  real(8)   :: DEGN,DTHETA,THETA,DRD,RDL,X,Y,Z
  real(8)   :: ENDR,DFZ,DFZMIN,XPT,YPT,ZPT,RX,RY,ZTEMP
  integer:: J
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
        WRITE(6,*) '## TYPE: C/CIRCLE  A/ARC  S/SPIRAL P/POINTS'//&
             &            '  X/EXIT'
        READ(5,'(A1)',ERR=2,END=1) KID
        CALL GUCPTL(KID)
        
        IF(KID.EQ.'C') THEN
3          WRITE(6,602) RD,NJMAX,ZANT
602        FORMAT(' ','## RD,NJMAX,ZANT = ',F10.3,I5,F10.3)
           READ(5,*,ERR=3,END=2) RD,NJMAX,ZANT
           DTHETA=2.D0*PI/(NJMAX-1)
           DO NJ=1,NJMAX
              THETA=DTHETA*(NJ-1)+1.D-5
              XJ0(NJ,NA)=RD*COS(THETA)
              YJ0(NJ,NA)=RD*SIN(THETA)
              ZJ0(NJ,NA)=ZANT
           END DO
           JNUM0(NA)=NJMAX
           
        ELSEIF(KID.EQ.'A') THEN
4          WRITE(6,603) THETJ1,THETJ2,RD,NJMAX,ZANT
603        FORMAT(' ','## THETJ1,THETJ2 = ',2F10.3/&
                  ' ','## RD,NJMAX = ',F10.3,I5/&
                  ' ','## ZANT = ',F10.3)
           READ(5,*,ERR=4,END=2) THETJ1,THETJ2,RD,NJMAX,ZANT
           THETA=DEGN*THETJ1
           XJ0(1,NA)=1.5D0*RD*COS(THETA)
           YJ0(1,NA)=1.5D0*RD*SIN(THETA)
           ZJ0(1,NA)=ZANT
           DTHETA=(THETJ2-THETJ1)/(NJMAX-3)
           DO NJ=2,NJMAX-1
              THETA=DEGN*(DTHETA*(NJ-2)+THETJ1)
              Xj0(NJ,NA)=RD*COS(THETA)
              YJ0(NJ,NA)=RD*SIN(THETA)
              ZJ0(NJ,NA)=ZANT
           END DO
           THETA=DEGN*THETJ2
           XJ0(NJMAX,NA)=1.5D0*RD*COS(THETA)
           YJ0(NJMAX,NA)=1.5D0*RD*SIN(THETA)
           ZJ0(NJMAX,NA)=ZANT
           JNUM0(NA)=NJMAX
!
!
! ----- Add. By YOKOYAMA 28/02/2013 ----
!
!        Ellipsoidal Antenna
!
         ELSEIF(KID.EQ.'E') THEN
            ENDR = 0.5D0
    9       WRITE(6,605) ENDR,THETJ1,THETJ2,RD,NJMAX,ZANT
  605       FORMAT(' ','## ENDR = ',F10.3/&
                   ' ','## THETJ1,THETJ2 = ',2F10.3/&
                   ' ','## RD,NJMAX = ',F10.3,I5/&
                   ' ','## ZANT = ',F10.3)
            READ(5,*,ERR=9,END=2) ENDR,THETJ1,THETJ2,RD,NJMAX,ZANT
!
!          Find the cross-section of magnetic flux tube with the closest.
!          与えられたZ座標に，最も近いZ座標を持つ磁力管断面を探す．
!           FLZ,FLX,FLY -> 'cm' unit
!           ZPT,XPT,YPT ->  'm' unit
            DFZMIN = 1.D2
            DO J=1,NGFLIN
               ZPT = FLZ(J)/1.D2
               DFZ = ABS(ZANT-ZPT)
               IF(DFZ.LT.DFZMIN) THEN
                  DFZMIN = DFZ
                  XPT = FLX(J)/1.D2
                  YPT = FLY(J)/1.D2
                  ZTEMP = ZPT
               ENDIF
            ENDDO
!
            RX = RD/RA * XPT
            RY = RD/RA * YPT
            THETA=DEGN*THETJ1
            XJ0(1,NA)=RX*COS(THETA)+(ENDR-RX)
            YJ0(1,NA)=RY*SIN(THETA)
            ZJ0(1,NA)=ZANT
            DTHETA=(THETJ2-THETJ1)/(NJMAX-3)
            DO NJ=2,NJMAX-1
               THETA=DEGN*(DTHETA*(NJ-2)+THETJ1)
               XJ0(NJ,NA)=RX*COS(THETA)
               YJ0(NJ,NA)=RY*SIN(THETA)
               ZJ0(NJ,NA)=ZANT
            ENDDO
            THETA=DEGN*THETJ2
            XJ0(NJMAX,NA)=RX*COS(THETA)-(ENDR-RX)
            YJ0(NJMAX,NA)=RY*SIN(THETA)
            ZJ0(NJMAX,NA)=ZANT
            JNUM0(NA)=NJMAX
!
! ----- -----
!
        ELSEIF(KID.EQ.'S') THEN
5          WRITE(6,604) THETS1,THETS2,RD1,RD2,NJMAX,ZANT,ZWALL
604        FORMAT(' ','## THETS1,THETJS2 = ',2F10.3/&
                  ' ','## RD1,RD2,NJMAX = ',2F10.3,I5/&
                  ' ','## ZANT,ZWALL = ',2F10.3)
           READ(5,*,ERR=5,END=2) THETS1,THETS2,RD1,RD2,NJMAX,ZANT,ZWALL
               
           THETA=DEGN*THETS1
           XJ0(1,NA)=RD1*COS(THETA)
           YJ0(1,NA)=RD1*SIN(THETA)
           ZJ0(1,NA)=ZWALL
           DTHETA=(THETS2-THETS1)/(NJMAX-3)
           DRD   =(RD2   -RD1   )/(NJMAX-3)
           DO NJ=2,NJMAX-1
              THETA=DEGN*(DTHETA*(NJ-2)+THETS1)
              RDL  =      DRD   *(NJ-2)+RD1
              XJ0(NJ,NA)=RDL*COS(THETA)
              YJ0(NJ,NA)=RDL*SIN(THETA)
              ZJ0(NJ,NA)=ZANT
           END DO
           THETA=DEGN*THETS2
           XJ0(NJMAX,NA)=RD2*COS(THETA)
           YJ0(NJMAX,NA)=RD2*SIN(THETA)
           ZJ0(NJMAX,NA)=ZWALL
           JNUM0(NA)=NJMAX
!
! ----- Add. By YOKOYAMA 28/02/2013 ----
!
!        Antenna Along Mob-B Surface
!
         ELSEIF(KID.EQ.'Y') THEN
            CALL MBDEFA(NA)
!
!        Bar-Type Antenna
!
         ELSEIF(KID.EQ.'B') THEN
            CALL DEFBTA(NA)
!
! ----- -----
!           
        ELSEIF(KID.EQ.'P') THEN
6          WRITE(6,*) '## NUMBER OF POINTS : NJMAX=',NJMAX
           READ(5,*,ERR=6,END=2) NJMAX
           DO NJ=1,NJMAX
7             WRITE(6,*) '## NO.',NJ,': X,Y,Z ?'
              READ(5,*,ERR=7,END=6) X,Y,Z
              XJ0(NJ,NA)=X
              YJ0(NJ,NA)=Y
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

  call mtx_barrier
  call wfant_broadcast

  CALL MODANT(IERR)
  IF(IERR.NE.0) GOTO 9000
  
9000 RETURN
END SUBROUTINE WFDEFA

!------------------------------------------------

subroutine wfant_broadcast
  
  use libmpi
  use libmtx
  use wfcomm
  implicit none
  
  integer :: NA,NJ
  integer,dimension(NAM)::idata
  real(8),dimension(NJM)::ddatax,ddatay,ddataz

  if(nrank.eq.0) idata(1)=NAMAX
  call mtx_broadcast_integer(idata,1)
  NAMAX=idata(1)

  call mtx_broadcast_integer(JNUM0,NAMAX)
  
  do NA=1,NAMAX

     NJMAX=JNUM0(NA)
     if(nrank.eq.0) then
        do NJ=1,NJMAX
           ddatax(NJ)=XJ0(NJ,NA)
           ddatay(NJ)=YJ0(NJ,NA)
           ddataz(NJ)=ZJ0(NJ,NA)
        end do
     end if
     
     call mtx_broadcast_real8(ddatax,NJMAX)
     call mtx_broadcast_real8(ddatay,NJMAX)
     call mtx_broadcast_real8(ddataz,NJMAX)

     do NJ=1,NJMAX
        XJ0(NJ,NA)=ddatax(NJ)
        YJ0(NJ,NA)=ddatay(NJ)
        ZJ0(NJ,NA)=ddataz(NJ)
     end do

  end do

  return
end subroutine wfant_broadcast
