!
! Test program to read OPEN-ADAS ADF11 files (Unresolved only)
!
PROGRAM test_adf11

  USE ADF11
  USE libgrf
  IMPLICIT NONE
  INTEGER:: IERR,ND,IND,IZ,IS,ID,IT,NXMAX,NYMAX,NX,NY,MODE_2D,IZ0,IC
  INTEGER:: NXMAX_SAVE,NYMAX_SAVE,NPTMAX_SAVE
  INTEGER:: NL,NLMAX,NDI,NDR,NPTMAX,I,ND1,ND2
  INTEGER:: NZ,NZL,NZMIN,NZMAX,NZSTEP,NZCOUNT,NT,NTMAX
  REAL(dp):: PN,PT,DR,DRMIN,DRMAX,XMIN,XMAX,DX,YMIN,YMAX,DY
  REAL(dp):: PZAV,PZ2AV,PWBAV,PWLAV,DRR,DRI
  REAL(dp):: PNMIN,PNMAX,PTMIN,PTMAX,DPT,FMIN,FMAX,PTSTEP
  REAL(dp),DIMENSION(:),ALLOCATABLE:: XDATA,YDATA,PNZ
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
  NXMAX_SAVE=201
  NYMAX_SAVE=201
  NPTMAX_SAVE=5

1 CONTINUE
  WRITE(6,*) '## Input Data Number IZ0 (0 for quit):'
  READ(5,*,ERR=1,END=9000) IZ0
  IF(IZ0.LE.0) GOTO 9000
  ND=ND_TABLE(IZ0,1)
  IF(ND.EQ.0) THEN
     WRITE(6,*) 'XX No Data for IZ0=',IZ0
     GO TO 1
  END IF

2 CONTINUE
  WRITE(6,*) '## Input Data Type: 1..12: ICLASS, '
  WRITE(6,*) '   21: calc Z/Z2/PB/PL 22: PNZ vs Z, 23: Z/Z2/PB/PL: vs PT'
  WRITE(6,*) '   31: DRI-DRR vs PT, 32: DRI-DRR vs Z,  0: for end'
  READ(5,*,ERR=2,END=1) IC
  IF(IC.LE.0) GOTO 1

  SELECT CASE(IC)
  CASE(1:12)
     ND=ND_TABLE(IZ0,IC)
     IF(ND.EQ.0) GO TO 2
     DRMIN=DRCOFA(1,1,1,ND)+14.0
     DRMAX=DRCOFA(1,1,1,ND)+14.0
     DO IS=1,IZMAXA(ND)
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

3    CONTINUE

     WRITE(6,*) '## Input Plot type IND and graph type MODE_2D (IND=0 for end)'
     WRITE(6,*) '    IND=1:DR(PN,PT), 2:DR(PT,IZ), 3:spline-1, 4:spline-2'
     WRITE(6,*) '    MODE_2D=0:1D, 1:Cont, 2:Paint, 3:Cand P, 11:Bird eye view'
     READ(5,*,ERR=3,END=2) IND,MODE_2D
     IF(IND.LE.0) GO TO 2

     SELECT CASE(IND)
     CASE(1)
11      WRITE(6,'(A,2I5)') &
             '## Input IZ (0 for end):IZMIN,IZMAX=', &
             IZMINA(ND),IZMAXA(ND)
        READ(5,*,ERR=11,END=3) IZ
        IF(IZ.LE.0) GO TO 3
        NXMAX=ITMAXA(ND)
        NYMAX=IDMAXA(ND)
        ALLOCATE(XDATA(NXMAX),YDATA(NYMAX),FDATA(NXMAX,NYMAX))
        XDATA(1:NXMAX)=DTEMPA(1:ITMAXA(ND),ND)-3.D0
        YDATA(1:NYMAX)=DDENSA(1:IDMAXA(ND),ND)-14.D0
        DO NY=1,NYMAX
           DO NX=1,NXMAX
              FDATA(NX,NY)=DRCOFA(NY,NX,IZ,ND)+14.D0
           END DO
        END DO
        CALL PAGES
        CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
                   '@DR(T) for n@',3,0,MODE_2D)
        CALL PAGEE
        DEALLOCATE(XDATA,YDATA,FDATA)
        GO TO 11

     CASE(2)
12      WRITE(6,'(A,1P2E12.4)') &
             '## Input log_10 PN (-10 for end): PNMIN,PNMAX=', &
             DDENSA(1,ND)-14.D0,DDENSA(IDMAXA(ND),ND)-14.D0
        READ(5,*,ERR=12,END=3) PN
        IF(PN.LE.-10.D0) GOTO 3
        NXMAX=ITMAXA(ND)
        NYMAX=IZMAXA(ND)
        ALLOCATE(XDATA(NXMAX),YDATA(NYMAX),FDATA(NXMAX,NYMAX))
        XDATA(1:NXMAX)=DTEMPA(1:NXMAX,ND)-3.D0
        DO NY=1,NYMAX
           IZ=NY-1+IZMINA(ND)
           YDATA(NY)=DBLE(IZ)
           DO NX=1,NXMAX
              PT=XDATA(NX)
              CALL CALC_ADF11(ND,IZ,PN,PT,DR,IERR)
              IF(IERR.NE.0) THEN
                 WRITE(6,*) 'XX test-adf11:CALC_ADF11: IERR=',IERR
                 GOTO 12
              END IF
              FDATA(NX,NY)=MAX(DR,-6.D0)
           END DO
        END DO
        CALL PAGES
        CALL GRD2D(0,XDATA,YDATA,FDATA,NXMAX,NXMAX,NYMAX, &
             '@DR(T) for Z@',1,0,MODE_2D)
        CALL PAGEE
        DEALLOCATE(XDATA,YDATA,FDATA)
        GO TO 12

     CASE(3)
        NXMAX=NXMAX_SAVE
        NYMAX=NYMAX_SAVE
13      WRITE(6,'(A)') &
             '## Input IZ,NXMAX,NYMAX (IZ=0 for end):'
        WRITE(6,'(A,4I5)') &
          '      IZMIN,IZMAX,NXMAX,NYMAX=', &
                 IZMINA(ND),IZMAXA(ND),NXMAX,NYMAX
        READ(5,*,ERR=13,END=3) IZ,NXMAX,NYMAX
        IF(IZ.LE.0) GO TO 3
        IF(NXMAX.LE.1) GO TO 13
        IF(NYMAX.LE.1) GO TO 13
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
              IF(IERR.NE.0) GOTO 13
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
        GO TO 13

     CASE(4)
        IF(IZMAXA(ND).EQ.1) THEN
           WRITE(6,*) 'XX IZMAXA(ND)=1: no contour plot available'
           GO TO 3
        END IF
        NXMAX=NXMAX_SAVE
        NYMAX=NYMAX_SAVE
14      WRITE(6,'(A)') &
             '## Input log_10 PN (-10 for end):'
        WRITE(6,'(A,1P2E12.4,2I5)') &
             '      PNMIN,PNMAX,NXMAX,NYMAX=', &
             DDENSA(1,ND)-14.D0,DDENSA(IDMAXA(ND),ND)-14.D0,NXMAX,NYMAX
        READ(5,*,ERR=14,END=3) PN,NXMAX,NYMAX
        IF(PN.LE.-10.D0) GO TO 3
        IF(NXMAX.LE.1) GO TO 14
        IF(NYMAX.LE.1) GO TO 14
        NXMAX_SAVE=NXMAX
        NYMAX_SAVE=NYMAX
        XMIN=DTEMPA(         1,ND)-3.D0
        XMAX=DTEMPA(ITMAXA(ND),ND)-3.D0
        DX=(XMAX-XMIN)/(NXMAX-1)
        YMIN=DBLE(IZMINA(ND))
        YMAX=DBLE(IZMAXA(ND))
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
              IF(IERR.NE.0) GOTO 14
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
        GO TO 14
     END SELECT

  CASE(21)
     WRITE(6,'(A)') '## calc Z/Z2/PB/PL'
     NDI=ND_TABLE(IZ0,2)  ! ionization rate
     WRITE(6,'(A,1P2E12.4)') &
          'Ionize Log10_PN(n20): ', &
          DDENSA(1,NDI)-14.D0,DDENSA(IDMAXA(NDI),NDI)-14.D0
     WRITE(6,'(A,1P2E12.4)') &
          'Ionize Log10_PT(keV): ', &
          DTEMPA(1,NDI)-3.D0,DTEMPA(ITMAXA(NDI),NDI)-3.D0
     NDR=ND_TABLE(IZ0,1)  ! recombination rate
     WRITE(6,'(A,1P2E12.4)') &
          'Recomb Log10_PN(n20): ', &
          DDENSA(1,NDR)-14.D0,DDENSA(IDMAXA(NDR),NDR)-14.D0
     WRITE(6,'(A,1P2E12.4)') &
          'Recomb Log10_PT(keV): ', &
          DTEMPA(1,NDR)-3.D0, DTEMPA(ITMAXA(NDR),NDR)-3.D0

21   CONTINUE

     WRITE(6,'(A)') '## Input Log10_PN and Log10_PT (-10 for end):'
     READ(5,*,ERR=21,END=9000) PN,PT
     IF(PN.LE.-10.D0) GO TO 2
     CALL IONIZE_EQ2(IZ0,PN,PT,PZAV,PZ2AV,PWBAV,PWLAV,2,IERR)

     WRITE(6,'(A,1PE12.4)') 'Z    AV = ',PZAV
     WRITE(6,'(A,1PE12.4)') 'Z**2 AV = ',PZ2AV
     WRITE(6,'(A,1PE12.4)') 'PWB  AV = ',PWBAV
     WRITE(6,'(A,1PE12.4)') 'PWL  AV = ',PWLAV

     GOTO 21
  
  CASE(22)
     WRITE(6,'(A)') '## PNZ vs Z'
     NPTMAX=NPTMAX_SAVE
     NDI=ND_TABLE(IZ0,2)  ! ionization rate
     NDR=ND_TABLE(IZ0,1)  ! recombination rate
     PNMIN=MIN(DDENSA(1,NDI),DDENSA(1,NDR))-14.D0
     PNMAX=MAX(DDENSA(IDMAXA(NDI),NDI),DDENSA(IDMAXA(NDR),NDR))-14.D0
     PTMIN=MIN(DTEMPA(1,NDI),DTEMPA(1,NDR))-3.D0
     PTMAX=MAX(DTEMPA(ITMAXA(NDI),NDI),DTEMPA(ITMAXA(NDR),NDR))-3.D0
     WRITE(6,'(A,I5)')       '## NPTMAX= ',NPTMAX
     WRITE(6,'(A,1P2E12.4)') '## PNMIN,PNMAX= ',PNMIN,PNMAX
     WRITE(6,'(A,1P2E12.4)') '## PTMIN,PTMAX= ',PTMIN,PTMAX

22   CONTINUE
     WRITE(6,'(A)') '## Input NPTMAX,PN,PTMIN,PTMAX (NPTMAX=0 for end):'
     READ(5,*,ERR=2,END=9000) NPTMAX,PN,PTMIN,PTMAX
     IF(NPTMAX.LE.0) GO TO 2

     NPTMAX_SAVE=NPTMAX
     NYMAX=NPTMAX
     NXMAX=IZMAXA(NDI)
     ALLOCATE(XDATA(0:NXMAX),YDATA(NYMAX),FDATA(0:NXMAX,NYMAX),PNZ(0:NXMAX))
     DO NX=0,NXMAX
        XDATA(NX)=DBLE(NX)
     END DO
     DPT=(PTMAX-PTMIN)/(NYMAX-1)
     DO NY=1,NYMAX
        PT=PTMIN+DPT*(NY-1)
        CALL IONIZE_EQ1(IZ0,PN,PT,PNZ,0,IERR)
        YDATA(NY)=PT
        DO NX=0,NXMAX
           FDATA(NX,NY)=MAX(LOG10(PNZ(NX)),-10.D0)
        END DO
     END DO

     DO NY=1,NYMAX
        WRITE(6,'(A,I2,A,1PE12.4)') 'PT(',NY,')=',YDATA(NY)
     END DO

     CALL PAGES
     CALL GRD1D(0,XDATA,FDATA,NXMAX+1,NXMAX+1,NYMAX,'@PNZ(Z) for PT@',2, &
                FMIN=-8.D0)
     CALL PAGEE

     DEALLOCATE(XDATA,YDATA,FDATA,PNZ)
     GOTO 22

  CASE(23)
     WRITE(6,'(A)') '## AVE Z/Z2/PB/PL vs PT'
     NXMAX=NXMAX_SAVE
     NDI=ND_TABLE(IZ0,2)  ! ionization rate
     NDR=ND_TABLE(IZ0,1)  ! recombination rate
     PNMIN=MIN(DDENSA(1,NDI),DDENSA(1,NDR))-14.D0
     PNMAX=MAX(DDENSA(IDMAXA(NDI),NDI),DDENSA(IDMAXA(NDR),NDR))-14.D0
     PTMIN=MIN(DTEMPA(1,NDI),DTEMPA(1,NDR))-3.D0
     PTMAX=MAX(DTEMPA(ITMAXA(NDI),NDI),DTEMPA(ITMAXA(NDR),NDR))-3.D0
     WRITE(6,'(A,I5)')       '## NXMAX= ',NXMAX
     WRITE(6,'(A,1P2E12.4)') '## PNMIN,PNMAX= ',PNMIN,PNMAX
     WRITE(6,'(A,1P2E12.4)') '## PTMIN,PTMAX= ',PTMIN,PTMAX

23   CONTINUE
     WRITE(6,'(A)') '## Input NXMAX,PN,PTMIN,PTMAX (NXMAX=0 for end):'
     READ(5,*,ERR=2,END=9000) NXMAX,PN,PTMIN,PTMAX
     IF(NXMAX.LE.0) GO TO 2

     NXMAX_SAVE=NXMAX
     ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,4))
     DPT=(PTMAX-PTMIN)/(NXMAX-1)
     DO NX=1,NXMAX
        PT=PTMIN+DPT*(NX-1)
        CALL IONIZE_EQ2(IZ0,PN,PT,PZAV,PZ2AV,PWBAV,PWLAV,0,IERR)
        XDATA(NX)=PT
        FDATA(NX,1)=PZAV
        FDATA(NX,2)=PZ2AV
        FDATA(NX,3)=LOG10(PWBAV)
        FDATA(NX,4)=LOG10(PWLAV)
!        WRITE(6,'(I5,1P5E12.4)') NX,XDATA(NX),(FDATA(NX,I),I=1,4)
     END DO
     FMAX=DBLE(NINT(MAXVAL(FDATA(1:NXMAX,4))))
     FMIN=FMAX-3.D0
     CALL PAGES
     CALL GRD1D(1,XDATA,FDATA(1:NXMAX,1),NXMAX,NXMAX,1,'@PZAV vs PT@',1)
     CALL GRD1D(3,XDATA,FDATA(1:NXMAX,2),NXMAX,NXMAX,1,'@PZ2AV vs PT@',1)
     CALL GRD1D(2,XDATA,FDATA(1:NXMAX,3),NXMAX,NXMAX,1,'@PWBAV vs PT@',3, &
                FMIN=FMIN,FMAX=FMAX)
     CALL GRD1D(4,XDATA,FDATA(1:NXMAX,4),NXMAX,NXMAX,1,'@PWLAV vs PT@',3, &
                FMIN=FMIN,FMAX=FMAX)
     CALL PAGEE

     DEALLOCATE(XDATA,FDATA)
     GOTO 23

     CASE(31)
31      WRITE(6,'(A)') '## DRI-DRR vs PT for Zs'
        NTMAX=101
        PTMIN=0.5D0
        PTMAX=6.5D0
        PN=-1.D0
        NZMIN=25
        NZMAX=45
        NZSTEP=5
        WRITE(6,'(A)') &
             '## Input NTMAX,PTmin,PTmax,log(PN): NTMAX=0 for end'
        READ(5,*,ERR=31,END=2) NTMAX,PTMIN,PTMAX,PN
        IF(NTMAX.LE.0) GO TO 2
        IF(NTMAX.EQ.1) GO TO 31
        ND1=ND_TABLE(IZ0,1)
        ND2=ND_TABLE(IZ0,2)

        WRITE(6,'(A,I5,A,I5,A)') &
             '## Input NZMIN,NZMAX,NZSTEP: ', &
             IZMINA(ND1),' <= NZ <= ',IZMAXA(ND1)
        READ(5,*,ERR=31,END=2) NZMIN,NZMAX,NZSTEP
        IF(NZMAX.LT.IZMINA(ND1)) GO TO 31
        IF(NZMAX.GT.IZMAXA(ND1)) GO TO 31
        NZCOUNT=(NZMAX-NZMIN)/NZSTEP+1

        ALLOCATE(XDATA(NTMAX),FDATA(NTMAX,NZCOUNT))

        PTSTEP=(PTMAX-PTMIN)/(NTMAX-1)
        DO NT=1,NTMAX
           PT=PTMIN+PTSTEP*(NT-1)
           XDATA(NT)=PT
           DO NZ=NZMIN,NZMAX,NZSTEP
              NZL=(NZ-NZMIN)/NZSTEP+1
              CALL CALC_ADF11(ND1,NZ,PN,LOG10(PT),DRR,IERR)
              CALL CALC_ADF11(ND2,NZ,PN,LOG10(PT),DRI,IERR)
              FDATA(NT,NZL)=10.D0**DRI-10.D0**DRR
           END DO
        END DO

        CALL PAGES
        CALL GRD1D(0,XDATA,FDATA,NTMAX,NTMAX,NZCOUNT,'@DRI-DRR vs PT@',0, &
                   XMIN=0.D0,FMIN=-1.5D5,FMAX=1.5D5, &
                   XSCALE_TYPE=0,XSCALE_SIZE=1.D0, &
                   FSCALE_TYPE=0,FSCALE_SIZE=1.D0)
        CALL PAGEE
        DEALLOCATE(XDATA,FDATA)
        GO TO 31

     CASE(32)
32      WRITE(6,'(A)') '## DRI-DRR vs Z for PTs'
        ND1=ND_TABLE(IZ0,1)
        ND2=ND_TABLE(IZ0,2)
        NZMIN=1
        NZMAX=IZ0
        NTMAX=6
        PTMIN=1.D0
        PTMAX=6.D0
        PN=-1.D0
        WRITE(6,'(A,I5,A,I5,A)') &
             '## Input NZMIN,NZMAX: ', &
             IZMINA(ND1),' <= NZ <= ',IZMAXA(ND1),' : NZMIN=0 for end'
        READ(5,*,ERR=32,END=2) NZMIN,NZMAX
        IF(NZMIN.LE.0) GO TO 2
        IF(NZMAX.LT.IZMINA(ND1)) GO TO 32
        IF(NZMAX.GT.IZMAXA(ND1)) GO TO 32

        WRITE(6,'(A)') &
             '## Input NTMAX,PTmin,PTmax,log(PN): NTMAX=0 for end'
        READ(5,*,ERR=32,END=2) NTMAX,PTMIN,PTMAX,PN
        IF(NTMAX.LE.0) GO TO 2
        IF(NTMAX.EQ.1) GO TO 32
        ND1=ND_TABLE(IZ0,1)
        ND2=ND_TABLE(IZ0,2)

        NZCOUNT=NZMAX-NZMIN+1

        ALLOCATE(XDATA(NZCOUNT),FDATA(NZCOUNT,NTMAX))

        PTSTEP=(PTMAX-PTMIN)/(NTMAX-1)
        DO NZ=NZMIN,NZMAX
           XDATA(NZ)=DBLE(NZ)
           DO NT=1,NTMAX
              PT=PTMIN+PTSTEP*(NT-1)
              CALL CALC_ADF11(ND1,NZ,PN,LOG10(PT),DRR,IERR)
              CALL CALC_ADF11(ND2,NZ,PN,LOG10(PT),DRI,IERR)
              FDATA(NZ,NT)=10.D0**DRI-10.D0**DRR
           END DO
        END DO

        CALL PAGES
        CALL GRD1D(0,XDATA,FDATA,NZCOUNT,NZCOUNT,NTMAX,'@DRI-DRR vs Z@',0, &
                   XMIN=20.D0,XMAX=60.D0,FMIN=-0.5D5,FMAX=0.5D5,&
                   XSCALE_TYPE=0,XSCALE_SIZE=1.D0, &
                   FSCALE_TYPE=0,FSCALE_SIZE=1.D0)
        CALL PAGEE
        DEALLOCATE(XDATA,FDATA)
        GO TO 32
     END SELECT
  
9000 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM test_adf11
