Module wimgout
  PRIVATE
  PUBLIC wim_gout
 
  interface
     real(4) function GUCLIP(X)
       real(8):: X
     end function GUCLIP
     integer(4) function NGULEN(Y)
       real(4):: Y
     end function NGULEN
  end interface

CONTAINS

  SUBROUTINE wim_gout

    USE wimcomm,ONLY: rkind
    USE wimparm,ONLY: wim_parm
    USE libgrf
    IMPLICIT NONE
    INTEGER           :: ierr,mode,kid,ich
    CHARACTER         :: kch
    CHARACTER(LEN=80) :: line

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') &
         '#### WIM GOUT: 1,2,3,4 X/exit'
    CALL TASK_KLIN(line,kid,mode,wim_parm)
    IF(mode == 2 .OR. mode == 3) GOTO 1

    ICH=ICHAR(LINE(1:1))
    IF(ICH.GE.97.AND.ICH.LE.122) ICH=ICH-32
    KCH=CHAR(ICH)

    SELECT CASE(kch)
    CASE('1') ! ER,EL,EZ: scaled
       CALL WIMGRA(1,1)
    CASE('2') ! ER,EL,EZ
       CALL WIMGRA(1,2)
    CASE('3') ! EX,EY,EZ: sclaed
       CALL WIMGRA(2,1)
    CASE('4') ! EX,EY,EZ
       CALL WIMGRA(2,2)
    CASE('X')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wim_gout

!     *****  GRAPHIC  *****

  SUBROUTINE WIMGRA(ID1,ID2)

    USE wimcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: ID1  ! 1 for ER,EL,EZ
                              ! 2 for EX,EY,EZ
    INTEGER,INTENT(IN):: ID2  ! 1 scaled by total
                              ! 2 scaled by each component
    INTEGER:: NXMAX,NX
    REAL(4):: XMAX,GAMIN,GAMAX,GMIN,GMAX,SCAL,GAMIN1,GAMIN2,GAMAX1,GAMAX2, &
              GAMAX3,SCALX,GMAXN,SCALN
    REAL(4),ALLOCATABLE,DIMENSION(:):: GX,GDATA1,GDATA2,GDATA3
    REAL(rkind):: TL1,TL2,TR1,TR2,RL1,RL2,RR1,RR2,PTL1,PTL2,PTR1,PTR2, &
                  PT,PRL1,PRL2,PRR1,PRR2
    COMPLEX(rkind):: CEX,CEY,CBY,CBX,CPRIN,CPLIN,CEP,CEM
    
    XMAX=GUCLIP(ZMAX)
    NXMAX=NZMAX+1
    ALLOCATE(GX(NXMAX),GDATA1(NXMAX),GDATA2(NXMAX),GDATA3(NXMAX))

    CALL PAGES
    CALL SETMKS(0,0.1)
    CALL SETCHS(0.3,0.)
    CALL SETLIN(0,2,7)
    CALL MOVE( 2.5,17.7)
    IF(ID1.EQ.1) THEN
       CALL TEXT('ER,EL,EZ',8)
    ELSE IF(ID1.EQ.2) THEN
       CALL TEXT('EX,EY,EZ',8)
    END IF

    DO NX=1,NXMAX
       GX    (NX)=ZA(NX-1)
    END DO
    CALL GMNMX1(GX,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCALX)

    IF(ID2.EQ.1) THEN
       IF(ID1.EQ.1) THEN
          DO NX=1,NXMAX
             GDATA1(NX)=ABS(CE(3*NX)+CI*CE(3*NX+1))/SQRT(2.D0)
             GDATA2(NX)=ABS(CE(3*NX)-CI*CE(3*NX+1))/SQRT(2.D0)
             GDATA3(NX)=+ABS(CE(3*NX+2))
          END DO
       ELSE IF(ID1.EQ.2) THEN
          DO NX=1,NXMAX
             GDATA1(NX)=ABS(CE(3*NX))
             GDATA2(NX)=ABS(CE(3*NX+1))
             GDATA3(NX)=ABS(CE(3*NX+2))
          END DO
       END IF
       CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX1)
       CALL GMNMX1(GDATA2,1,NXMAX,1,GAMIN,GAMAX2)
       CALL GMNMX1(GDATA3,1,NXMAX,1,GAMIN,GAMAX3)
       GAMAX=MAX(GAMAX1,GAMAX2,GAMAX3)
       CALL GQSCAL(0.0,GAMAX,GMIN,GMAXN,SCALN)
    ENDIF
       

    DO NX=1,NXMAX
       IF(ID1.EQ.1) THEN
          CEP=(CE(3*NX)+CI*CE(3*NX+1))/SQRT(2.D0)
       ELSEIF(ID1.EQ.2) THEN
          CEP=CE(3*NX)
       ENDIF
       GDATA1(NX)= ABS(CEP)
       GDATA2(NX)=REAL(CEP)
       GDATA3(NX)=AIMAG(CEP)
    END DO

    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(0.0,GAMAX,GMIN,GMAX,SCAL)
    IF(ID2.EQ.1) THEN
       GMAX=GMAXN
       SCAL=SCALN
    END IF
    CALL GDEFIN(2.5,12.5,12.0,17.5,0.0,XMAX,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,2.*SCALX,0.,0.,0.2,9)
    CALL GSCALE(0.,0.,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,1,4)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,5)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GX, GDATA3,1,NXMAX,1, 0, 0,2)
    CALL SETLIN(0,2,7)

    DO NX=1,NXMAX
       IF(ID1.EQ.1) THEN
          CEM=(CE(3*NX)-CI*CE(3*NX+1))/SQRT(2.D0)
       ELSE IF(ID1.EQ.2) THEN
          CEM=CE(3*NX+1)
       END IF

       GDATA1(NX)= ABS(CEM)
       GDATA2(NX)=REAL(CEM)
       GDATA3(NX)=IMAG(CEM)
    END DO

    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(0.0,GAMAX,GMIN,GMAX,SCAL)
    IF(ID2.EQ.1) THEN
       GMAX=GMAXN
       SCAL=SCALN
    END IF
    CALL GDEFIN(2.5,12.5,6.5,12.,0.0,XMAX,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,2*SCALX,0.,0.,0.2,9)
    CALL GSCALE(0.,0.,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,1,4)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,5)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GX, GDATA3,1,NXMAX,1, 0, 0,2)
    CALL SETLIN(0,2,7)

    DO NX=1,NXMAX
       GDATA1(NX)= ABS(CE(3*NX+2))
       GDATA2(NX)=REAL(CE(3*NX+2))
       GDATA3(NX)=IMAG(CE(3*NX+2))
    END DO

    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(0.0,GAMAX,GMIN,GMAX,SCAL)
    IF(ID2.EQ.1) THEN
       GMAX=GMAXN
       SCAL=SCALN
    END IF
    CALL GDEFIN(2.5,12.5,1.,6.5,0.0,XMAX,-GMAX,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,2.*SCALX,0.,0.,0.2,9)
    CALL GSCALE(0.,0.,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,4*SCALX,0.,0.,NGULEN(4*SCALX))
    CALL SETLIN(0,1,4)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,5)
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,5)
    CALL GPLOTP(GX, GDATA3,1,NXMAX,1, 0, 0,2)
    CALL SETLIN(0,2,7)

    DO NX=1,NXMAX
       GDATA1(NX)=REAL(CPWR(NX-1))
    END DO

    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN,GAMAX)
    CALL GQSCAL(0.0,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,1.0,6.5,0.,XMAX,0.0,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.,0.,2.*SCAL,0.1,1)
    CALL GSCALE(0.,2.*SCALX,0.,0.,0.2,9)
    CALL GSCALE(0.,0.,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL GVALUE(0.,4*SCALX,0.,0.,NGULEN(4*SCALX))
    CALL SETLIN(0,2,6)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1, 0, 0,0)
    CALL SETLIN(0,2,7)

    DO NX=1,NXMAX
       GDATA1(NX)=GUCLIP(AKD1(NX-1))
       GDATA2(NX)=GUCLIP(AKD2(NX-1))
    END DO

    CALL GMNMX1(GDATA1,1,NXMAX,1,GAMIN1,GAMAX1)
    CALL GMNMX1(GDATA2,1,NXMAX,1,GAMIN2,GAMAX2)
    GAMIN=MIN(GAMIN1,GAMIN2)
    GAMAX=MAX(GAMAX1,GAMAX2)
    CALL GQSCAL(GAMIN,GAMAX,GMIN,GMAX,SCAL)
    CALL GDEFIN(15.0,25.0,6.5,12.0,0.,XMAX,GMIN,GMAX)
    CALL GFRAME
    CALL GSCALE(0.,0.,0.,2.*SCAL,0.1,9)
    CALL GSCALE(0.,2.*SCALX,0.,0.,0.2,9)
    CALL GSCALE(0.,0.,0.,13.*SCAL,0.1,0)
    CALL GVALUE(0.,0.,0.,4*SCAL,NGULEN(4*SCAL))
    CALL SETLIN(0,1,6)
    CALL GPLOTP(GX, GDATA1,1,NXMAX,1,1,20,0)
    CALL GPLOTP(GX, GDATA2,1,NXMAX,1,2,20,0)
    CALL SETLIN(0,2,7)

    CALL MOVE(13.5,18.0)
    CALL TEXT('  NZMAX = ',10)
    CALL NUMBI(NZMAX,'(I9)',9)
    CALL TEXT('  NWMAX = ',10)
    CALL NUMBI(NWMAX,'(I9)',9)
    CALL MOVE(13.5,17.5)
    CALL TEXT('  PN0   = ',10)
    CALL NUMBD(PN0,'(F9.5)',9)
    CALL TEXT('  ZMAX  = ',10)
    CALL NUMBD(ZMAX,'(F9.4)',9)
    CALL MOVE(13.5,17.0)
    CALL TEXT('  PB0   = ',10)
    CALL NUMBD(PB0,'(F9.5)',9)
    CALL TEXT('  DZMAX = ',10)
    CALL NUMBD(DZMAX,'(F9.4)',9)
    CALL MOVE(13.5,16.5)
    CALL TEXT('  PB1   = ',10)
    CALL NUMBD(PB1,'(F9.5)',9)
    CALL TEXT('  DZWID = ',10)
    CALL NUMBD(DZWID,'(F9.4)',9)
    CALL MOVE(13.5,16.0)
    CALL TEXT('  ANX   = ',10)
    CALL NUMBD(ANX,'(F9.5)',9)
    CALL MOVE(13.5,15.5)
    CALL TEXT('  BETA  = ',10)
    CALL NUMBD(BETA,'(F9.5)',9)

    CALL MOVE(13.5,15.0)
    CALL TEXT('  IN    = ',10)
    TL1=ABS(CER1)**2 
    CALL NUMBD(TL1,'(F7.4)',7)
    TL2=ABS(CEL1)**2 
    CALL NUMBD(TL2,'(F7.4)',7)
    TR1=ABS(CER2)**2 
    CALL NUMBD(TR1,'(F7.4)',7)
    TR2=ABS(CEL2)**2 
    CALL NUMBD(TR2,'(F7.4)',7)
    CALL MOVE(13.5,14.5)
    CALL TEXT('  REFL  = ',10)
    RL1=ABS(CE(1))**2 
    CALL NUMBD(RL1,'(F7.4)',7)
    RL2=ABS(CE(2))**2 
    CALL NUMBD(RL2,'(F7.4)',7)
    RR1=ABS(CE(NZMAX*3+6))**2 
    CALL NUMBD(RR1,'(F7.4)',7)
    RR2=ABS(CE(NZMAX*3+7))**2 
    CALL NUMBD(RR2,'(F7.4)',7)
    CALL MOVE(13.5,14.0)
    CALL TEXT('  P-IN  = ',10)
    PTL1=TL1*IMAG(CAL+CDL*CONJG(CFL))
    PTL2=TL2*IMAG(CBL+CEL*CONJG(CGL)) 
    PTR1=TR1*IMAG(CAR+CDR*CONJG(CFR))
    PTR2=TR2*IMAG(CBR+CER*CONJG(CGR)) 
    PT=PTR1+PTR2+PTL1+PTL2
    CALL NUMBD(PTL1/PT,'(F7.4)',7)
    CALL NUMBD(PTL2/PT,'(F7.4)',7)
    CALL NUMBD(PTR1/PT,'(F7.4)',7)
    CALL NUMBD(PTR2/PT,'(F7.4)',7)
    CALL MOVE(13.5,13.5)
    CALL TEXT('  P-REF = ',10)
    PRL1=RL1*IMAG(CAL+CDL*CONJG(CFL))
    PRL2=RL2*IMAG(CBL+CEL*CONJG(CGL)) 
    PRR1=RR1*IMAG(CAR+CDR*CONJG(CFR))
    PRR2=RR2*IMAG(CBR+CER*CONJG(CGR)) 
    CALL NUMBD(PRL1/PT,'(F7.4)',7)
    CALL NUMBD(PRL2/PT,'(F7.4)',7)
    CALL NUMBD(PRR1/PT,'(F7.4)',7)
    CALL NUMBD(PRR2/PT,'(F7.4)',7)
    CALL MOVE(13.5,13.0)
    CALL TEXT('  P-ABS = ',10)
    CALL NUMBD(DBLE(CPTOT)/PT,'(F7.4)',7)
    CALL TEXT('  P-IN ',7)
    CEX=     CER1+CE(1) +     CEL1+CE(2)
    CEY=CFL*(CER1+CE(1))+CGL*(CEL1+CE(2))
    CBY=CAL*(CER1-CE(1))+CBL*(CEL1-CE(2))
    CBX=CDL*(CER1-CE(1))+CEL*(CEL1-CE(2))
    CPLIN=CONJG(CEX)*CBY+CONJG(CEY)*CBX
    CALL NUMBD(IMAG(CPLIN)/PT,'(F7.4)',7)
    CEX=     CER2+CE(NZMAX*3+6) +     CEL2+CE(NZMAX*3+7)
    CEY=CFR*(CER2+CE(NZMAX*3+6))+CGR*(CEL2+CE(NZMAX*3+7))
    CBY=CAR*(CER2-CE(NZMAX*3+6))+CBR*(CEL2-CE(NZMAX*3+7))
    CBX=CDR*(CER2-CE(NZMAX*3+6))+CER*(CEL2-CE(NZMAX*3+7))
    CPRIN=CONJG(CEX)*CBY+CONJG(CEY)*CBX
    CALL NUMBD(IMAG(CPRIN)/PT,'(F7.4)',7)

    CALL PAGEE
    RETURN
  END SUBROUTINE wimgra
END Module wimgout
