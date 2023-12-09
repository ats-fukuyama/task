!    $Id$

Module wigout
  PRIVATE
  PUBLIC wi_gout,wi_mesh,wi_gra1
 
  interface
     real function GUCLIP(X)
       USE bpsd_kinds,ONLY: rkind
       real(rkind):: X
     end function GUCLIP
     integer function NGULEN(Y)
       real:: Y
     end function NGULEN
  end interface

CONTAINS

  SUBROUTINE wi_gout

    USE wicomm,ONLY: rkind
    USE wiparm
    USE libgrf
    IMPLICIT NONE
    CHARACTER(LEN=1):: kch
    INTEGER:: IERR

    ierr=0
!    WRITE(6,'(A)') &
!         '#### WI GOUT:  R/1D  X/exit'
!    CALL TASK_KLIN(line,kid,mode,wi_parm)
!    IF(mode == 2 .OR. mode == 3) GOTO 1
!
!    ICH=ICHAR(LINE(1:1))
!    IF(ICH.GE.97.AND.ICH.LE.122) ICH=ICH-32
!    KCH=CHAR(ICH)

    kch='R'

    SELECT CASE(kch)
    CASE('R') 
       CALL wi_gra1
       CALL wi_gra2
    CASE('X') 
       GO TO 9000
    END SELECT

!    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wi_gout

!     *****  GRAPHIC   *****

  SUBROUTINE wi_gra1

    USE wicomm
    IMPLICIT NONE
    COMPLEX(rkind):: CZ
    REAL,DIMENSION(nxmax+1):: SVX,SCR,SCI,SCA,SCER,SCEI,SCEA,SPOWR 
    REAL:: SMINX,SMAXX,SMINF,SMAXF,SMINE,SMAXE,SMINP,SMAXP
    REAL:: SGMINX,SGMAXX,SGMINF,SGMAXF,SGMINE,SGMAXE,SGMINP,SGMAXP
    REAL:: SCALX,SCALF,SCALE,SCALP
    REAL(rkind):: ANB
    REAL(rkind):: R,T,S
    INTEGER:: J,JD
    EXTERNAL PAGES,SETLIN,SETCHS,SETFNT,GMNMX1,GQSCAL,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,MOVE,TEXT,NUMBI,NUMBD
    

    DO J=1,NXMAX+1 
       JD=2*(J-1) 
       SVX(J)=GUCLIP(xgrid(J-1))
       SCR(J)=GUCLIP(REAL(CFY(JD+1)))
       SCI(J)=GUCLIP(AIMAG(CFY(JD+1)))
       SCA(J)=GUCLIP(ABS(CFY(JD+1)))
       SCER(J)=GUCLIP(REAL(CFY(JD+2)))
       SCEI(J)=GUCLIP(AIMAG(CFY(JD+2)))
       SCEA(J)=GUCLIP(ABS(CFY(JD+2)))
       SPOWR(J)=GUCLIP(REAL(CPOWER(J-1)))
    END DO

    CALL PAGES
    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.)
    CALL SETFNT(32)
    CALL GMNMX1(SVX,1,NXMAX+1,1,SMINX,SMAXX)
    CALL GMNMX1(SCA,1,NXMAX+1,1,SMINF,SMAXF)
    CALL GMNMX1(SCEA,1,NXMAX+1,1,SMINE,SMAXE)
    CALL GMNMX1(SPOWR,1,NXMAX+1,1,SMINP,SMAXP)
    CALL GQSCAL(SMINX,SMAXX,SGMINX,SGMAXX,SCALX)
    CALL GQSCAL(-SMAXF,SMAXF,SGMINF,SGMAXF,SCALF)
    CALL GQSCAL(-SMAXE,SMAXE,SGMINE,SGMAXE,SCALE)
    CALL GQSCAL(SMINP,SMAXP,SGMINP,SGMAXP,SCALP)
    CALL GDEFIN(3.,12.,10.0,17.0,SMINX,SMAXX,SGMINF,SGMAXF)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALF,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALF,NGULEN(2*SCALF))
    CALL SETLIN(0,0,6)
!    CALL GPLOTP(SVX,SCR,1,NXMAX+1,1,0,0,1)
    CALL GPLOTP(SVX,SCR,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,5)
!    CALL GPLOTP(SVX,SCI,1,NXMAX+1,1,0,0,3)
    CALL GPLOTP(SVX,SCI,1,NXMAX+1,1,0,0,2)
!    CALL SETLIN(0,0,4)
!    CALL GPLOTP(SVX,SCA,1,NXMAX+1,1,0,0,5)
    CALL SETLIN(0,0,7)
    CALL GDEFIN(3.,12.,1.0,8.0,SMINX,SMAXX,SGMINE,SGMAXE)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALE,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALE,NGULEN(2*SCALE))
    CALL SETLIN(0,0,6)
!    CALL GPLOTP(SVX,SCER,1,NXMAX+1,1,0,0,1)
    CALL GPLOTP(SVX,SCER,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,5)
!    CALL GPLOTP(SVX,SCEI,1,NXMAX+1,1,0,0,3)
    CALL GPLOTP(SVX,SCEI,1,NXMAX+1,1,0,0,2)
!    CALL SETLIN(0,0,4)
!    CALL GPLOTP(SVX,SCEA,1,NXMAX+1,1,0,0,5)
    CALL SETLIN(0,0,7)
    CALL GDEFIN(16.0,25.,1.0,8.0,SMINX,SMAXX,SGMINP,SGMAXP)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALP,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALP,NGULEN(2*SCALP))
    CALL SETLIN(0,0,6)
!    CALL GPLOTP(SVX,SPOWR,1,NXMAX+1,1,0,0,1)
    CALL GPLOTP(SVX,SPOWR,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,7)

    CALL MOVE(15.5,17.0)
    CALL TEXT('  NXMAX = ',10)
    CALL NUMBI(NXMAX,'(I7)',8)
    CALL MOVE(15.5,16.5)
    CALL TEXT('  NWMAX = ',10)
    CALL NUMBI(NWMAX,'(I7)',8)
    CALL MOVE(15.5,16.0)
    CALL TEXT('  XMAX  = ',10)
    CALL NUMBD(XMAX,'(F7.1)',7)
    CALL MOVE(15.5,15.5)
    CALL TEXT('  PN0   = ',10)
    CALL NUMBD(PN0,'(F7.3)',7)
    CALL MOVE(15.5,15.0)
    CALL TEXT('  ALFA  = ',10)
    CALL NUMBD(ALFA,'(F9.5)',9)
    CALL MOVE(15.5,14.5)
    CALL TEXT('  ANY   = ',10)
    CALL NUMBD(ANY,'(F7.3)',7)
    CALL MOVE(15.5,14.0)
    CALL TEXT('  BETA  = ',10)
    CALL NUMBD(BETA,'(F7.3)',7)
    CALL MOVE(15.5,13.5)
    CALL TEXT('  PNU   = ',10)
    CALL NUMBD(PNU,'(F7.3)',7)
    CALL MOVE(15.5,12.75)

    CALL TEXT('  R     = ',10)
    R=ABS(CFY(NXMAX*2+3))**2
    CALL NUMBD(R,'(F9.5)',9)
    CALL MOVE(15.5,12.25)
    CALL TEXT('  A     = ',10)
    T=1-R
    CALL NUMBD(T,'(F9.5)',9)
    CALL MOVE(15.5,11.5)
    CALL TEXT('  P-IN  = ',10)

    IF(ALFA*xgrid(nxmax).GT.100.D0) THEN
       ANB=0.D0
    ELSE
       ANB=DEXP(-ALFA*xgrid(nxmax))
    END IF
    IF(1.D0-ANB-ANY*ANY.LT.0.D0) THEN
       S=0.D0
    ELSE
       S=T/(DSQRT(1.D0-ANB-ANY*ANY)) 
    END IF
    CALL NUMBD(S,'(ES12.4)',12)
    CALL MOVE(15.5,11.0)
    CALL TEXT('  P-ABS = ',10)
    CALL NUMBD(PTOT,'(ES12.4)',12)
    CALL MOVE(15.5,10.0)
    CALL TEXT('   CZ   = ',10)
    CZ=(1.D0+CFY(NXMAX*2+3))*DSQRT(1.D0-ANY*ANY)/(1.D0-CFY(NXMAX*2+3)) 
    S=REAL(CZ)
    CALL NUMBD(S,'(F9.5)',9)
    CALL MOVE(15.5, 9.5)
    CALL TEXT('       +i ',10)
    S=AIMAG(CZ)
    CALL NUMBD(S,'(F9.5)',9)
    CALL MOVE(6.0,17.2)
    CALL TEXT('< CEX >',7)
    CALL MOVE(6.0,8.2)
    CALL TEXT('< CEY >',7)
    CALL MOVE(19.4,8.2)
    CALL TEXT('< POWER >',9)
    CALL PAGEE
    RETURN
  END SUBROUTINE wi_gra1

!     *****  GRAPHIC   *****

  SUBROUTINE wi_gra2

    USE wicomm
    USE libspf2,ONLY: airyb
    IMPLICIT NONE
    COMPLEX(rkind):: CZ
    REAL,DIMENSION(nxmax+1):: SVX,SCR,SCI,SCA,SCER,SCEI,SCEA,SPOWR
    REAL,DIMENSION(nxmax+1):: SVA,SVAI
    REAL(rkind),DIMENSION(nxmax+1):: VA,VAI,VDAI,VBI,VDBI
    REAL:: SMINX,SMAXX,SMINF,SMAXF,SMINE,SMAXE,SMINP,SMAXP
    REAL:: SGMINX,SGMAXX,SGMINF,SGMAXF,SGMINE,SGMAXE,SGMINP,SGMAXP
    REAL:: SCALX,SCALF,SCALE,SCALP
    REAL(rkind):: ANB
    REAL(rkind):: R,T,S
    INTEGER:: J,JD
    EXTERNAL PAGES,SETLIN,SETCHS,SETFNT,GMNMX1,GQSCAL,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,MOVE,TEXT,NUMBI,NUMBD
    

    DO J=1,NXMAX+1 
       JD=2*(J-1) 
       SVX(J)=GUCLIP(xgrid(J-1))
       VA(J)=-xgrid(J-1)*alfa**(1.D0/3.D0)
       CALL airyb(VA(J),VAI(J),VDAI(J),VBI(J),VDBI(J))
       SVA(J)=GUCLIP(VA(J))
       SVAI(J)=GUCLIP(2.D0*VAI(J))
       SCR(J)=GUCLIP(REAL(CFY(JD+1)))
       SCI(J)=GUCLIP(AIMAG(CFY(JD+1)))
       SCA(J)=GUCLIP(ABS(CFY(JD+1)))
       SCER(J)=GUCLIP(REAL(CFY(JD+2)))
       SCEI(J)=GUCLIP(AIMAG(CFY(JD+2)))
       SCEA(J)=GUCLIP(ABS(CFY(JD+2)))
       SPOWR(J)=GUCLIP(REAL(CPOWER(J-1)))
    END DO

    CALL PAGES
    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.)
    CALL SETFNT(32)
    CALL GMNMX1(SVX,1,NXMAX+1,1,SMINX,SMAXX)
    CALL GMNMX1(SCA,1,NXMAX+1,1,SMINF,SMAXF)
    CALL GMNMX1(SCEA,1,NXMAX+1,1,SMINE,SMAXE)
    CALL GMNMX1(SPOWR,1,NXMAX+1,1,SMINP,SMAXP)
    CALL GQSCAL(SMINX,SMAXX,SGMINX,SGMAXX,SCALX)
    CALL GQSCAL(-SMAXF,SMAXF,SGMINF,SGMAXF,SCALF)
    CALL GQSCAL(-SMAXE,SMAXE,SGMINE,SGMAXE,SCALE)
    CALL GQSCAL(SMINP,SMAXP,SGMINP,SGMAXP,SCALP)
    CALL GDEFIN(3.,12.,10.0,17.0,SMINX,SMAXX,SGMINF,SGMAXF)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALF,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALF,NGULEN(2*SCALF))
    CALL SETLIN(0,0,6)
    CALL GPLOTP(SVX,SCR,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,4)
    CALL GPLOTP(SVX,SVAI,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(SVX,SCI,1,NXMAX+1,1,0,0,2)
    CALL SETLIN(0,0,7)
    CALL GDEFIN(3.,12.,1.0,8.0,SMINX,SMAXX,SGMINE,SGMAXE)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALE,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALE,NGULEN(2*SCALE))
    CALL SETLIN(0,0,6)
    CALL GPLOTP(SVX,SCER,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,4)
    CALL GPLOTP(SVX,SVAI,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(SVX,SCEI,1,NXMAX+1,1,0,0,2)
    CALL SETLIN(0,0,7)
    CALL GDEFIN(16.0,25.,1.0,8.0,SMINX,SMAXX,SGMINP,SGMAXP)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALP,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALP,NGULEN(2*SCALP))
    CALL SETLIN(0,0,6)
    CALL GPLOTP(SVX,SPOWR,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,7)

    CALL MOVE(15.5,17.0)
    CALL TEXT('  NXMAX = ',10)
    CALL NUMBI(NXMAX,'(I7)',8)
    CALL MOVE(15.5,16.5)
    CALL TEXT('  NWMAX = ',10)
    CALL NUMBI(NWMAX,'(I7)',8)
    CALL MOVE(15.5,16.0)
    CALL TEXT('  XMAX  = ',10)
    CALL NUMBD(XMAX,'(F7.1)',7)
    CALL MOVE(15.5,15.5)
    CALL TEXT('  PN0   = ',10)
    CALL NUMBD(PN0,'(F7.3)',7)
    CALL MOVE(15.5,15.0)
    CALL TEXT('  ALFA  = ',10)
    CALL NUMBD(ALFA,'(F9.5)',9)
    CALL MOVE(15.5,14.5)
    CALL TEXT('  ANY   = ',10)
    CALL NUMBD(ANY,'(F7.3)',7)
    CALL MOVE(15.5,14.0)
    CALL TEXT('  BETA  = ',10)
    CALL NUMBD(BETA,'(F7.3)',7)
    CALL MOVE(15.5,13.5)
    CALL TEXT('  PNU   = ',10)
    CALL NUMBD(PNU,'(F7.3)',7)
    CALL MOVE(15.5,12.75)

    CALL TEXT('  R     = ',10)
    R=ABS(CFY(NXMAX*2+3))**2
    CALL NUMBD(R,'(F9.5)',9)
    CALL MOVE(15.5,12.25)
    CALL TEXT('  A     = ',10)
    T=1-R
    CALL NUMBD(T,'(F9.5)',9)
    CALL MOVE(15.5,11.5)
    CALL TEXT('  P-IN  = ',10)

    IF(ALFA*xgrid(nxmax).GT.100.D0) THEN
       ANB=0.D0
    ELSE
       ANB=DEXP(-ALFA*xgrid(nxmax))
    END IF
    IF(1.D0-ANB-ANY*ANY.LT.0.D0) THEN
       S=0.D0
    ELSE
       S=T/(DSQRT(1.D0-ANB-ANY*ANY)) 
    END IF
    CALL NUMBD(S,'(ES12.4)',12)
    CALL MOVE(15.5,11.0)
    CALL TEXT('  P-ABS = ',10)
    CALL NUMBD(PTOT,'(ES12.4)',12)
    CALL MOVE(15.5,10.0)
    CALL TEXT('   CZ   = ',10)
    CZ=(1.D0+CFY(NXMAX*2+3))*DSQRT(1.D0-ANY*ANY)/(1.D0-CFY(NXMAX*2+3)) 
    S=REAL(CZ)
    CALL NUMBD(S,'(F9.5)',9)
    CALL MOVE(15.5, 9.5)
    CALL TEXT('       +i ',10)
    S=AIMAG(CZ)
    CALL NUMBD(S,'(F9.5)',9)
    CALL MOVE(6.0,17.2)
    CALL TEXT('< CEX >',7)
    CALL MOVE(6.0,8.2)
    CALL TEXT('< CEY >',7)
    CALL MOVE(19.4,8.2)
    CALL TEXT('< POWER >',9)
    CALL PAGEE
    RETURN
  END SUBROUTINE wi_gra2

!     *****  GRAPHIC   *****

  SUBROUTINE wi_gra3

    USE wicomm
    IMPLICIT NONE
    COMPLEX(rkind):: CZ
    REAL,DIMENSION(nxmax+1):: SVX,SCR,SCI,SCA,SCER,SCEI,SCEA,SPOWR 
    REAL:: SMINX,SMAXX,SMINF,SMAXF,SMINE,SMAXE,SMINP,SMAXP
    REAL:: SGMINX,SGMAXX,SGMINF,SGMAXF,SGMINE,SGMAXE,SGMINP,SGMAXP
    REAL:: SCALX,SCALF,SCALE,SCALP
    REAL(rkind):: ANB
    REAL(rkind):: R,T,S
    INTEGER:: J,JD
    EXTERNAL PAGES,SETLIN,SETCHS,SETFNT,GMNMX1,GQSCAL,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,GPLOTP,MOVE,TEXT,NUMBI,NUMBD
    

    DO J=1,NXMAX+1 
       SVX(J)=GUCLIP(xgrid(J-1))
    END DO

    J=1
    JD=0
    CBX(J)=(CFY(JD+3)-CFY(JD+1))/(xgrid(J)-xgrid(J-1))
    CBY(J)=(CFY(JD+4)-CFY(JD+2))/(xgrid(J)-xgrid(J-1))
    DO J=2,NXMAX
       JD=2*(J-1) 
       CBX(J)=(CFY(JD+3)-CFY(JD-1))/(xgrid(J)-xgrid(J-2))
       CBY(J)=(CFY(JD+4)-CFY(JD  ))/(xgrid(J)-xgrid(J-2))
    END DO
    J=NXMAX+1
    JD=2*(J-1) 
    CBX(J)=(CFY(JD+1)-CFY(JD-1))/(xgrid(1)-xgrid(J-2))
    CBY(J)=(CFY(JD+2)-CFY(JD  ))/(xgrid(1)-xgrid(J-2))
       
       SCR(J)=GUCLIP(REAL(CFY(JD+1)))
       SCI(J)=GUCLIP(AIMAG(CFY(JD+1)))
       SCA(J)=GUCLIP(ABS(CFY(JD+1)))
       SCER(J)=GUCLIP(REAL(CFY(JD+2)))
       SCEI(J)=GUCLIP(AIMAG(CFY(JD+2)))
       SCEA(J)=GUCLIP(ABS(CFY(JD+2)))
       SPOWR(J)=GUCLIP(REAL(CPOWER(J-1)))
    END DO

    CALL PAGES
    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.)
    CALL SETFNT(32)
    CALL GMNMX1(SVX,1,NXMAX+1,1,SMINX,SMAXX)
    CALL GMNMX1(SCA,1,NXMAX+1,1,SMINF,SMAXF)
    CALL GMNMX1(SCEA,1,NXMAX+1,1,SMINE,SMAXE)
    CALL GMNMX1(SPOWR,1,NXMAX+1,1,SMINP,SMAXP)
    CALL GQSCAL(SMINX,SMAXX,SGMINX,SGMAXX,SCALX)
    CALL GQSCAL(-SMAXF,SMAXF,SGMINF,SGMAXF,SCALF)
    CALL GQSCAL(-SMAXE,SMAXE,SGMINE,SGMAXE,SCALE)
    CALL GQSCAL(SMINP,SMAXP,SGMINP,SGMAXP,SCALP)
    CALL GDEFIN(3.,12.,10.0,17.0,SMINX,SMAXX,SGMINF,SGMAXF)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALF,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALF,NGULEN(2*SCALF))
    CALL SETLIN(0,0,6)
    CALL GPLOTP(SVX,SCR,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(SVX,SCI,1,NXMAX+1,1,0,0,2)
    CALL SETLIN(0,0,7)
    CALL GDEFIN(3.,12.,1.0,8.0,SMINX,SMAXX,SGMINE,SGMAXE)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALE,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALE,NGULEN(2*SCALE))
    CALL SETLIN(0,0,6)
    CALL GPLOTP(SVX,SCER,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,5)
    CALL GPLOTP(SVX,SCEI,1,NXMAX+1,1,0,0,2)
    CALL SETLIN(0,0,7)
    CALL GDEFIN(16.0,25.,1.0,8.0,SMINX,SMAXX,SGMINP,SGMAXP)
    CALL GFRAME
    CALL GSCALE(0.,SCALX,0.,SCALP,0.3,1)
    CALL GSCALE(0.,10000.,0.,10000.,0.,0)
    CALL GVALUE(0.,2.*SCALX,0.,0.,NGULEN(2*SCALX))
    CALL GVALUE(0.,0.,0.,2.*SCALP,NGULEN(2*SCALP))
    CALL SETLIN(0,0,6)
    CALL GPLOTP(SVX,SPOWR,1,NXMAX+1,1,0,0,0)
    CALL SETLIN(0,0,7)

    CALL MOVE(15.5,17.0)
    CALL TEXT('  NXMAX = ',10)
    CALL NUMBI(NXMAX,'(I7)',8)
    CALL MOVE(15.5,16.5)
    CALL TEXT('  NWMAX = ',10)
    CALL NUMBI(NWMAX,'(I7)',8)
    CALL MOVE(15.5,16.0)
    CALL TEXT('  XMAX  = ',10)
    CALL NUMBD(XMAX,'(F7.1)',7)
    CALL MOVE(15.5,15.5)
    CALL TEXT('  PN0   = ',10)
    CALL NUMBD(PN0,'(F7.3)',7)
    CALL MOVE(15.5,15.0)
    CALL TEXT('  ALFA  = ',10)
    CALL NUMBD(ALFA,'(F9.5)',9)
    CALL MOVE(15.5,14.5)
    CALL TEXT('  ANY   = ',10)
    CALL NUMBD(ANY,'(F7.3)',7)
    CALL MOVE(15.5,14.0)
    CALL TEXT('  BETA  = ',10)
    CALL NUMBD(BETA,'(F7.3)',7)
    CALL MOVE(15.5,13.5)
    CALL TEXT('  PNU   = ',10)
    CALL NUMBD(PNU,'(F7.3)',7)
    CALL MOVE(15.5,12.75)

    CALL TEXT('  R     = ',10)
    R=ABS(CFY(NXMAX*2+3))**2
    CALL NUMBD(R,'(F9.5)',9)
    CALL MOVE(15.5,12.25)
    CALL TEXT('  A     = ',10)
    T=1-R
    CALL NUMBD(T,'(F9.5)',9)
    CALL MOVE(15.5,11.5)
    CALL TEXT('  P-IN  = ',10)

    IF(ALFA*xgrid(nxmax).GT.100.D0) THEN
       ANB=0.D0
    ELSE
       ANB=DEXP(-ALFA*xgrid(nxmax))
    END IF
    IF(1.D0-ANB-ANY*ANY.LT.0.D0) THEN
       S=0.D0
    ELSE
       S=T/(DSQRT(1.D0-ANB-ANY*ANY)) 
    END IF
    CALL NUMBD(S,'(ES12.4)',12)
    CALL MOVE(15.5,11.0)
    CALL TEXT('  P-ABS = ',10)
    CALL NUMBD(PTOT,'(ES12.4)',12)
    CALL MOVE(15.5,10.0)
    CALL TEXT('   CZ   = ',10)
    CZ=(1.D0+CFY(NXMAX*2+3))*DSQRT(1.D0-ANY*ANY)/(1.D0-CFY(NXMAX*2+3)) 
    S=REAL(CZ)
    CALL NUMBD(S,'(F9.5)',9)
    CALL MOVE(15.5, 9.5)
    CALL TEXT('       +i ',10)
    S=AIMAG(CZ)
    CALL NUMBD(S,'(F9.5)',9)
    CALL MOVE(6.0,17.2)
    CALL TEXT('< CEX >',7)
    CALL MOVE(6.0,8.2)
    CALL TEXT('< CEY >',7)
    CALL MOVE(19.4,8.2)
    CALL TEXT('< POWER >',9)
    CALL PAGEE
    RETURN
  END SUBROUTINE wi_gra3

  SUBROUTINE wi_mesh
    USE libgrf,ONLY: GRD1D
    USE wicomm,ONLY: ikind,rkind,nxmax,xgrid
    IMPLICIT NONE
    INTEGER(ikind):: nx
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: xid
    EXTERNAL PAGES,PAGEE

    CALL PAGES
    ALLOCATE(xid(nxmax))
    DO nx=1,nxmax
       xid(nx)=dble(nx)
    END DO
    CALL GRD1D(0,xid,xgrid,nxmax,nxmax,1,'@xgrid vs nx@')
    DEALLOCATE(xid)
    CALL PAGEE
    RETURN
  END SUBROUTINE wi_mesh
END Module wigout
