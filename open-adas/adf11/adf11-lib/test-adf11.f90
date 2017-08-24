!
! Test program to read OPEN-ADAS ADF11 files (Unresolved only)
!
PROGRAM test_adf11

  USE ADF11
  USE libgrf
  IMPLICIT NONE
  INTEGER:: IERR,ND,IND,IZ,IS,ID,IT,NXMAX,NYMAX,NX,NY,MODE_2D
  INTEGER:: NXMAX_SAVE,NYMAX_SAVE
  INTEGER:: NL,NLMAX
  REAL(dp):: PN,PT,DR,DRMIN,DRMAX,XMIN,XMAX,DX,YMIN,YMAX,DY
  REAL(dp),DIMENSION(:),ALLOCATABLE:: XDATA,YDATA
  REAL(dp),DIMENSION(:,:),ALLOCATABLE:: FDATA,LINE_RGB
  REAL(4):: RGB(3)

  CALL GSOPEN
  CALL LOAD_ADF11_bin(IERR)
  IF(IERR.NE.0) THEN
     WRITE(6,*) 'XX test_adf11: LOAD_ADF11_bin: IERR=',IERR
     STOP
  END IF
  IND=1
  MODE_2D=0
  ND=1
  NXMAX_SAVE=50
  NYMAX_SAVE=50

1 CONTINUE

  WRITE(6,*) '## Input Data Number ND (0 for quit):'
  READ(5,*,ERR=1,END=9000) ND
  IF(ND.LE.0) GOTO 9000
  DRMIN=DRCOFA(1,1,1,ND)+14.0
  DRMAX=DRCOFA(1,1,1,ND)+14.0
  DO IS=1,ISMAXA(ND)
     DO IT=1,ITMAXA(ND)
        DO ID=1,IDMAXA(ND)
           DRMIN=MIN(DRMIN,DRCOFA(ID,IT,IS,ND))
           DRMAX=MAX(DRMAX,DRCOFA(ID,IT,IS,ND))
        END DO
     END DO
  END DO
  WRITE(6,'(A,1P2E12.4)') &
       'Log10_PN(n20): ',DDENSA(1,ND)-14.D0,DDENSA(IDMAXA(ND),ND)-14.D0
  WRITE(6,'(A,1P2E12.4)') &
       'Log10_PT(keV): ',DTEMPA(1,ND)-3.D0,DTEMPA(ITMAXA(ND),ND)-3.D0
  WRITE(6,'(A,1P2E12.4)') &
       'Log10_DR(n20): ',DRMIN,DRMAX

2 CONTINUE

  WRITE(6,*) '## Input Plot type IND and graph type MODE_2D (IND=0 for end)'
  WRITE(6,*) '    IND=1:DR(PN,PT), 2:DR(PT,IZ), 3:spline-1, 4:spline-2'
  WRITE(6,*) '    MODE_2D=0:1D, 1:Contour, 2:Paint, 3:Cand P, 11:Bird eye view'
  READ(5,*,ERR=2,END=1) IND,MODE_2D
  IF(IND.LE.0) GOTO 1

  SELECT CASE(IND)
  CASE(1)
3    WRITE(6,'(A,2I5)') &
          '## Input IZ (-1 for end):IZMIN,IZMAX=', &
          IS1MINA(ND)-1,IS1MAXA(ND)-1
     READ(5,*,ERR=3,END=2) IZ
     IF(IZ.LT.0) GO TO 2
     IS=IZ+1
     NXMAX=ITMAXA(ND)
     NYMAX=IDMAXA(ND)
     ALLOCATE(XDATA(NXMAX),YDATA(NYMAX),FDATA(NXMAX,NYMAX))
     XDATA(1:NXMAX)=DTEMPA(1:ITMAXA(ND),ND)-3.D0
     YDATA(1:NYMAX)=DDENSA(1:IDMAXA(ND),ND)-14.D0
     DO NY=1,NYMAX
        DO NX=1,NXMAX
           FDATA(NX,NY)=DRCOFA(NY,NX,IS,ND)+14.D0
        END DO
     END DO
     CALL PAGES
     CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T,n)@',3,0,MODE_2D)
     CALL PAGEE
     DEALLOCATE(XDATA,YDATA,FDATA)
     GO TO 3

  CASE(2)
4    WRITE(6,'(A,1P2E12.4)') &
          '## Input log_10 PN (-10 for end): PNMIN,PNMAX=', &
          DDENSA(1,ND)-14.D0,DDENSA(IDMAXA(ND),ND)-14.D0
     READ(5,*,ERR=4,END=2) PN
     IF(PN.LE.-10.D0) GOTO 2
     NXMAX=ITMAXA(ND)
     NYMAX=ISMAXA(ND)
     ALLOCATE(XDATA(NXMAX),YDATA(NYMAX),FDATA(NXMAX,NYMAX))
     XDATA(1:NXMAX)=DTEMPA(1:NXMAX,ND)-3.D0
     DO NY=1,NYMAX
        IZ=NY-2+IS1MINA(ND)
        YDATA(NY)=DBLE(IZ)
        DO NX=1,NXMAX
           PT=XDATA(NX)
           CALL CALC_ADF11(ND,IZ,PN,PT,DR,IERR)
           IF(IERR.NE.0) THEN
              WRITE(6,*) 'XX test-adf11:CALC_ADF11: IERR=',IERR
              GOTO 4
           END IF
           FDATA(NX,NY)=MAX(DR,-6.D0)
        END DO
     END DO
     CALL PAGES
     CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T,Z)@',1,0,MODE_2D)
     CALL PAGEE
     DEALLOCATE(XDATA,YDATA,FDATA)
     GO TO 4

  CASE(3)
     NXMAX=NXMAX_SAVE
     NYMAX=NYMAX_SAVE
5    WRITE(6,'(A)') &
          '## Input IZ,NXMAX,NYMAX (IZ=-1 for end):'
     WRITE(6,'(A,4I5)') &
          '      IZMIN,IZMAX,NXMAX,NYMAX=', &
                 IS1MINA(ND)-1,IS1MAXA(ND)-1,NXMAX,NYMAX
     READ(5,*,ERR=5,END=2) IZ,NXMAX,NYMAX
     IF(IZ.LT.0) GO TO 2
     IF(NXMAX.LE.1) GO TO 5
     IF(NYMAX.LE.1) GO TO 5
     NXMAX_SAVE=NXMAX
     NYMAX_SAVE=NYMAX
     XMIN=DTEMPA(         1,ND)-3.D0
     XMAX=DTEMPA(ITMAXA(ND),ND)-3.D0
     DX=(XMAX-XMIN)/(NXMAX-1)
     YMIN=DDENSA(         1,ND)-14.D0
     YMAX=DDENSA(IDMAXA(ND),ND)-14.D0
     DY=(YMAX-YMIN)/(NYMAX-1)
     write(6,'(A,1P3E12.4)') 'X: ',XMIN,XMAX,DX
     write(6,'(A,1P3E12.4)') 'Y: ',YMIN,YMAX,DY
     ALLOCATE(XDATA(NXMAX),YDATA(NYMAX),FDATA(NXMAX,NYMAX))
     DO NX=1,NXMAX
        XDATA(NX)=XMIN+DX*(NX-1)
     END DO
     DO NY=1,NYMAX
        YDATA(NY)=YMIN+DY*(NY-1)
     END DO
     DO NY=1,NYMAX
        DO NX=1,NXMAX
           PN=YDATA(NY)
           PT=XDATA(NX)
           CALL CALC_ADF11(ND,IZ,PN,PT,DR,IERR)
           IF(IERR.NE.0) GOTO 5
           FDATA(NX,NY)=DR
        END DO
     END DO
     CALL PAGES
     NLMAX=11
     ALLOCATE(LINE_RGB(3,NLMAX))
     DO NL=1,11
        CALL R2Y2W(0.1*(NL-1),RGB)
        LINE_RGB(1:3,NL)=RGB(1:3)
     END DO
     IF(MODE_2D.EQ.0) THEN
        CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T,n)@',3,0,MODE_2D)
     ELSE
        CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T,n)@',3,0,MODE_2D, &
                   NLMAX=NLMAX,LINE_RGB=LINE_RGB)
     END IF
     DEALLOCATE(LINE_RGB)
     CALL PAGEE
     DEALLOCATE(XDATA,YDATA,FDATA)
     GO TO 5

  CASE(4)
     IF(ISMAXA(ND).EQ.1) THEN
        WRITE(6,*) 'XX ISMAXA(ND)=1: no contour plot available'
        GO TO 2
     END IF
     NXMAX=NXMAX_SAVE
     NYMAX=NYMAX_SAVE
6    WRITE(6,'(A)') &
          '## Input log_10 PN (-10 for end):'
     WRITE(6,'(A,1P2E12.4,2I5)') &
          '      PNMIN,PNMAX,NXMAX,NYMAX=', &
                 DDENSA(1,ND)-14.D0,DDENSA(IDMAXA(ND),ND)-14.D0,NXMAX,NYMAX
     READ(5,*,ERR=6,END=2) PN,NXMAX,NYMAX
     IF(PN.LE.-10.D0) GOTO 2
     IF(NXMAX.LE.1) GO TO 6
     IF(NYMAX.LE.1) GO TO 6
     NXMAX_SAVE=NXMAX
     NYMAX_SAVE=NYMAX
     XMIN=DTEMPA(         1,ND)-3.D0
     XMAX=DTEMPA(ITMAXA(ND),ND)-3.D0
     DX=(XMAX-XMIN)/(NXMAX-1)
     YMIN=DBLE(IS1MINA(ND))-1.D0
     YMAX=DBLE(IS1MAXA(ND))-1.D0
     DY=(YMAX-YMIN)/(NYMAX-1)
     write(6,'(A,1P3E12.4)') 'X: ',XMIN,XMAX,DX
     write(6,'(A,1P3E12.4)') 'Y: ',YMIN,YMAX,DY
     ALLOCATE(XDATA(NXMAX),YDATA(NYMAX),FDATA(NXMAX,NYMAX))
     DO NX=1,NXMAX
        XDATA(NX)=XMIN+DX*(NX-1)
     END DO
     DO NY=1,NYMAX
        YDATA(NY)=YMIN+DY*(NY-1)
     END DO
     DO NY=1,NYMAX
        DO NX=1,NXMAX
           IZ=NINT(YDATA(NY))
           PT=XDATA(NX)
           CALL CALC_ADF11(ND,IZ,PN,PT,DR,IERR)
           IF(IERR.NE.0) GOTO 6
           FDATA(NX,NY)=DR
        END DO
     END DO
     CALL PAGES
     NLMAX=11
     ALLOCATE(LINE_RGB(3,NLMAX))
     DO NL=1,11
        CALL R2Y2W(0.1*(NL-1),RGB)
        LINE_RGB(1:3,NL)=RGB(1:3)
     END DO
     IF(MODE_2D.EQ.0) THEN
        CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T,Z)@',1,0,MODE_2D)
     ELSE
        CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T,Z)@',1,0,MODE_2D, &
                   NLMAX=NLMAX,LINE_RGB=LINE_RGB)
     END IF
     DEALLOCATE(LINE_RGB)
     CALL PAGEE
     DEALLOCATE(XDATA,YDATA,FDATA)
     GO TO 6
  END SELECT
  GOTO 2
  
9000 CONTINUE
  STOP
END PROGRAM test_adf11
