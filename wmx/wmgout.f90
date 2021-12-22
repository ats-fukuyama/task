! wmgout.f90

MODULE wmgout

  PRIVATE
  PUBLIC wm_gout

CONTAINS

!     ****** CONTROL GRAPHICS ******

  SUBROUTINE wm_gout

    USE wmcomm
    USE wmgsub
    USE equnit_mod
    USE libchar
    IMPLICIT NONE
    CHARACTER(LEN=5):: KSTR
    CHARACTER(LEN=1):: K1,K2,K3,K4

1   WRITE(6,*) ' ## INPUT GSTR : R/EANSPB/ATMN  CPM/P/123  CP/J', &
               '  P/F/SBQ23J   R/GZ  S'
    WRITE(6,*) '                 CMP/EB/RTZsbh+-P/RIA', &
               '  G/01234  EQ  ?/HELP  X/EXIT'
    READ(5,'(A5)',ERR=1,END=900) KSTR
    K1=KSTR(1:1)

    CALL ToUpper(K1)
    IF (K1.EQ.'X') GOTO 900
    IF (K1.EQ.'G') THEN
       K2=KSTR(2:2)
       CALL ToUpper(K2)
       IF(K2.EQ.'0') NGRAPH=0
       IF(K2.EQ.'1') NGRAPH=1
       IF(K2.EQ.'2') NGRAPH=2
       IF(K2.EQ.'3') NGRAPH=3
       IF(K2.EQ.'4') NGRAPH=4
       GOTO 1
    ENDIF
    IF(K1.EQ.'?') THEN
       CALL wm_ghelp
       GOTO 1
    ENDIF

    IF((K1.EQ.'R').OR.(K1.EQ.'C').OR.(K1.EQ.'M').OR. &
       (K1.EQ.'P').OR.(K1.EQ.'S').OR.(K1.EQ.'E')) THEN
       K2=KSTR(2:2)
       CALL ToUpper(K2)
       K3=KSTR(3:3)
       CALL ToUpper(K3)
       K4=KSTR(4:4)
       CALL ToUpper(K4)
       
       IF(K1.EQ.'R') THEN
          IF(K2.EQ.'G') THEN
             CALL wm_greqg(K2,K3)
          ELSEIF(K2.EQ.'Z') THEN
!             CALL wm_gbooz
          ELSE
             CALL wm_gr1d(K2,K3)
          ENDIF
       ENDIF
       IF(K1.EQ.'P'.OR.(K1.EQ.'C'.AND.NGRAPH.EQ.0)) CALL wm_greq(K2,K3,K4)
       IF(K1.EQ.'C') CALL wm_grth(K2,K3,K4)
       IF(K1.EQ.'M') CALL wm_grmd(K2,K3,K4)
!       IF(K1.EQ.'S'.AND.(MODELG.EQ.4.OR.MODELG.EQ.6)) CALL wm_grms
       IF(K1.EQ.'E'.AND.K2.EQ.'Q'.AND.MODELG.GE.3) CALL eq_gout
    ELSE
       WRITE(6,*) '## UNDEFINED CONTROL CHARACTER: K1=',K1
    END IF
    GOTO 1

900 RETURN
  END SUBROUTINE wm_gout

!     ****** DRAW 1D GRAPH ******

  SUBROUTINE wm_gr1d(K2,K3)

    USE wmcomm
    USE wmgsub
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2,K3
    REAL,ALLOCATABLE:: GX1(:),GX2(:),GYS(:,:),GY(:,:),GY1(:,:),GY2(:,:)
    COMPLEX(rkind),ALLOCATABLE:: CF(:,:),CF1(:,:),CF2(:,:),CF3(:,:)
    REAL(rkind),ALLOCATABLE:: POWER(:,:),PF1(:,:),PF2(:,:)
    INTEGER:: NR,NX1,NX2,NG4,I,NS,NTH,NHH,MD,MD1,MD2,MDX,ND,NDX,IMAX
    REAL:: GXMIN,GXMAX,GXL
    CHARACTER(LEN=6):: KTITL(4)
    REAL:: GP(4,4)
    DATA GP/ 3.0, 10.8,  9.5, 16.5, &
             3.0, 10.8,  1.0,  8.0, &
            13.8, 21.6,  9.5, 16.5, &
            13.8, 21.6,  1.0,  8.0/
    EXTERNAL PAGES,SETCHS,MOVE,TEXT,NUMBI,DRAW,SETLIN,PAGEE

    ALLOCATE(GX1(nrmax+1),GX2(nrmax+1),GYS(nrmax+1,nsmax))
    ALLOCATE(GY(nrmax+1,3))
    ALLOCATE(GY1(nrmax+1,nthmax*nhhmax),GY2(nrmax+1,nthmax*nhhmax))
    ALLOCATE(CF(nrmax+1,3),CF1(nrmax+1,nthmax*nhhmax))
    ALLOCATE(CF2(nrmax+1,nthmax*nhhmax),CF3(nrmax+1,nthmax*nhhmax))
    ALLOCATE(POWER(nrmax+1,nsmax))
    ALLOCATE(PF1(nrmax+1,nthmax*nhhmax),PF2(nrmax+1,nthmax*nhhmax))

    GXMIN=GUCLIP(XR(1))
    GXMAX=GUCLIP(XR(NRMAX+1))

    SELECT CASE(K2)
    CASE('E','A','N','S','P') THEN
       NR=1
          GX1(NR)=GUCLIP(       XR(NR))
          GX2(NR)=GUCLIP(       XR(NR))
       DO NR=2,NRMAX+1
          GX1(NR)=GUCLIP(       XR(NR))
          GX2(NR)=GUCLIP(0.5D0*(XR(NR-1)+XR(NR)))
       ENDDO
       NX1=NRMAX+1
       NX2=NRMAX+1
       NG4=NSMAX
       KTITL(1)='Er    '
       KTITL(2)='Etheta'
       KTITL(3)='Ez    '
       KTITL(4)='Pabs  '
    CASE('A')
       NR=1
          GX1(NR)=GUCLIP(       XR(NR))
          GX2(NR)=GUCLIP(       XR(NR))
       DO NR=2,NRMAX+1
          GX1(NR)=GUCLIP(       XR(NR))
          GX2(NR)=GUCLIP(0.5D0*(XR(NR-1)+XR(NR)))
       ENDDO
       NX1=NRMAX+1
       NX2=NRMAX+1
       NG4=NSMAX
       KTITL(1)='Pabs  '
       KTITL(2)='Pabs 1'
       KTITL(3)='Pabs 2'
       KTITL(4)='Pabs 3'
    CASE('B')
       NR=1
          GX2(NR)=GUCLIP(       XR(NR))
          GX1(NR)=GUCLIP(       XR(NR))
       DO NR=2,NRMAX+1
          GX2(NR)=GUCLIP(       XR(NR))
          GX1(NR)=GUCLIP(0.5D0*(XR(NR-1)+XR(NR)))
       ENDDO
       NX1=NRMAX+1
       NX2=NRMAX+1
       NG4=1
       KTITL(1)='Br    '
       KTITL(2)='Btheta'
       KTITL(3)='Bz    '
       KTITL(4)='Jrf   '
    ENDIF


    SELECT CASE(K3)
    CASE('T','A')
1      IF(NTHMAX.EQ.1) THEN
          NTH=1
       ELSE
          WRITE(6,*) '## INPUT NTH : 1..',NTHMAX
          READ(5,*,ERR=1,END=9000) NTH
       END IF
       IF(NTH.LT.1.OR.NTH.GT.NTHMAX) THEN
          WRITE(6,*) 'XX ILLEGAL NTH'
          GOTO 1
       ENDIF
2      IF(NHHMAX.EQ.1) THEN
          NHH=1
       ELSE
          WRITE(6,*) '## INPUT NHH : 1..',NHHMAX
          READ(5,*,ERR=2,END=9000) NHH
       END IF
       IF(NHH.LT.1.OR.NHH.GT.NHHMAX) THEN
          WRITE(6,*) 'XX ILLEGAL NHH'
          GOTO 2
       ENDIF

       SELECT CASE(K2)
       CASE('E')
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CEFLD(I,NTH,NHH,NR)
             ENDDO
          ENDDO
       CASE('N')
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CEN(I,NTH,NHH,NR)
             ENDDO
          ENDDO
       CASE('S')
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CES(I,NTH,NHH,NR)
             ENDDO
          ENDDO
       CASE('P')
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CEP(I,NTH,NHH,NR)
             ENDDO
          ENDDO
       CASE('B')
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CBFLD(I,NTH,NHH,NR)
             ENDDO
          ENDDO
       END SELECT

       SELECT CASE(K3)
       CASE('A')
          DO NS=1,NSMAX
             DO NR=1,NRMAX
                POWER(NR,NS)=PABSR(NR,NS)
             ENDDO
             POWER(NRMAX+1,NS)=0.D0
          ENDDO
       CASE('T','M','N')
          DO NS=1,NSMAX
             DO NR=1,NRMAX
                POWER(NR,NS)=PABS(NTH,NHH,NR,NS)
             ENDDO
             POWER(NRMAX+1,NS)=0.D0
          ENDDO
       END SELECT
    ELSEIF(K3.EQ.'M') THEN
3      IF(NTHMAX.EQ.1) THEN
          MD=0
       ELSE
          MD1=NTH0+MDMIN
          MD2=NTH0+MDMAX-1
          WRITE(6,*) '## INPUT MD : ',MD1,'..',MD2
          READ(5,*,ERR=3,END=9000) MD
          MD=MD-NTH0
       END IF
       IF(MD.LT.MDMIN.OR.MD.GT.MDMAX) THEN
          WRITE(6,*) 'XX ILLEGAL MD'
          GOTO 3
       ENDIF
       MDX=MD-MDMIN+1
4      IF(NHHMAX.EQ.1) THEN
          ND=0
       ELSE
          WRITE(6,*) '## INPUT ND : ',NDMIN,'..',NDMAX-1
          READ(5,*,ERR=3,END=9000) ND
       END IF
       IF(ND.LT.NDMIN.OR.ND.GT.NDMAX) THEN
          WRITE(6,*) 'XX ILLEGAL ND'
          GOTO 4
       ENDIF
       NDX=ND-NDMIN+1

       IF(K2.EQ.'E') THEN
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CEFLDK(I,MDX,NDX,NR)
             ENDDO
          ENDDO
          DO NS=1,NSMAX
             DO NR=1,NRMAX
                POWER(NR,NS)=PABSK(MDX,NDX,NR,NS)
             ENDDO
             POWER(NRMAX+1,NS)=0.D0
          ENDDO
       ELSE
          DO I=1,3
             DO NR=1,NRMAX+1
                CF(NR,I)=CBFLDK(I,MDX,NDX,NR)
             ENDDO
          ENDDO
          DO NR=1,NRMAX
             POWER(NR,1)=PCUR(MDX,NDX,NR)
          ENDDO
          POWER(NRMAX+1,1)=0.D0
       ENDIF
    ELSEIF(K3.EQ.'N') THEN
7      IF(NTHMAX.EQ.1) THEN
          MD1=0
          MD2=0
       ELSE
          MD1=NTH0+MDMIN
          MD2=NTH0+MDMAX-1
          WRITE(6,*) '## INPUT MD1,MD2 : ',MD1,'..',MD2
          READ(5,*,ERR=7,END=9000) MD1,MD2
          MD1=MD1-NTH0
          MD2=MD2-NTH0
       END IF
       IF(MD1.LT.MDMIN.OR.MD1.GT.MD2.OR.MD2.GT.MDMAX) THEN
          WRITE(6,*) 'XX ILLEGAL MD1,MD2'
          GOTO 7
       ENDIF
8      IF(NHHMAX.EQ.1) THEN
          ND=0
       ELSE
          ND=NPH0
          WRITE(6,*) '## INPUT ND : ',ND
          READ(5,*,ERR=7,END=9000) ND
          ND=ND-NPH0
       END IF
       IF(ND.LT.NDMIN.OR.ND.GT.NDMAX) THEN
          WRITE(6,*) 'XX ILLEGAL ND'
          GOTO 8
       ENDIF
       NDX=ND-NDMIN+1

       IMAX=MD2-MD1+1
       IF(K2.EQ.'E') THEN
          DO I=1,IMAX
             MDX=MD1+I-1-MDMIN+1
             DO NR=1,NRMAX+1
                CF1(NR,I)=CEFLDK(1,MDX,NDX,NR)
                CF2(NR,I)=CEFLDK(2,MDX,NDX,NR)
                CF3(NR,I)=CEFLDK(3,MDX,NDX,NR)
                PF1(NR,I)=PABSK(MDX,NDX,NR,1)
                PF2(NR,I)=PABSK(MDX,NDX,NR,2)
             ENDDO
          ENDDO
       ELSE
          DO I=1,IMAX
             MDX=MD1+I-1-MDMIN+1
             DO NR=1,NRMAX+1
                CF1(NR,I)=CBFLDK(1,MDX,NDX,NR)
                CF2(NR,I)=CBFLDK(2,MDX,NDX,NR)
                CF3(NR,I)=CBFLDK(3,MDX,NDX,NR)
                PF1(NR,I)=PCUR(MDX,NDX,NR)
                PF2(NR,I)=0.D0
             ENDDO
          ENDDO
       ENDIF
    ELSEIF(K3.EQ.' ') THEN
       DO NS=1,NSMAX
          DO NR=1,NRMAX
             POWER(NR,NS)=PABSR(NR,NS)
          ENDDO
          POWER(NRMAX+1,NS)=0.D0
       ENDDO
    ENDIF

    CALL PAGES
    CALL SETCHS(0.3,0.0)

    IF(K3.EQ.'N') THEN

!      *** E/B(R) ****

       DO I=1,IMAX
          DO NR=1,NX2
             GY1(NR,I)=GUCLIP(DBLE(CF1(NR,I)))
             GY2(NR,I)=GUCLIP(DIMAG(CF1(NR,I)))
          ENDDO
       ENDDO
       CALL wm_gn1d(NX1,GX1,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,1),KTITL(1))

!      *** E/B(THETA) ****

       DO I=1,IMAX
          DO NR=1,NX1
             GY1(NR,I)=GUCLIP(DBLE(CF2(NR,I)))
             GY2(NR,I)=GUCLIP(DIMAG(CF2(NR,I)))
          ENDDO
       ENDDO
       CALL wm_gn1d(NX1,GX1,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,2),KTITL(2))

!      *** E/B(Z) ****

       DO I=1,IMAX
          DO NR=1,NX1
             GY1(NR,I)=GUCLIP(DBLE(CF3(NR,I)))
             GY2(NR,I)=GUCLIP(DIMAG(CF3(NR,I)))
          ENDDO
       ENDDO
       CALL wm_gn1d(NX1,GX1,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,3),KTITL(3))

!      *** POWER / CURRENT ***

       DO I=1,IMAX
          DO NR=1,NX2
             GY1(NR,I)=GUCLIP(PF1(NR,I))
             GY2(NR,I)=GUCLIP(PF2(NR,I))
          ENDDO
       ENDDO
       CALL wm_gn1d(NX2,GX2,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,4),KTITL(4))

    ELSEIF(K3.EQ.' ') THEN

!      *** POWER : All species ***

       DO I=1,NG4
          DO NR=1,NX2
             GYS(NR,I)=GUCLIP(POWER(NR,I))
          ENDDO
       ENDDO
       CALL wm_gsub(NX2,GX2,GXMIN,GXMAX,GYS,NG4,GP(1,1),KTITL(1))

!      *** POWER : species 1 ***

       DO NR=1,NX2
          GYS(NR,1)=GUCLIP(POWER(NR,1))
       ENDDO
       CALL wm_gsub(NX2,GX2,GXMIN,GXMAX,GYS,  1,GP(1,2),KTITL(2))

!      *** POWER : species 2 ***

       DO NR=1,NX2
          GYS(NR,1)=GUCLIP(POWER(NR,2))
       ENDDO
       CALL wm_gsub(NX2,GX2,GXMIN,GXMAX,GYS,  1,GP(1,3),KTITL(3))

!      *** POWER : species 3 ***

       IF(NRMAX.GE.3) THEN
          DO NR=1,NX2
             GYS(NR,1)=GUCLIP(POWER(NR,3))
          ENDDO
          CALL wm_gsub(NX2,GX2,GXMIN,GXMAX,GYS,  1,GP(1,4),KTITL(4))
       END IF
    ELSE

!      *** E/B(R) ****

       DO NR=1,NX2
          GY(NR,1)=GUCLIP(DBLE(CF(NR,1)))
          GY(NR,2)=GUCLIP(DIMAG(CF(NR,1)))
          GY(NR,3)=GUCLIP(ABS(CF(NR,1)))
       ENDDO
       CALL wm_gsub(NX1,GX1,GXMIN,GXMAX,GY,2,GP(1,1),KTITL(1))

!      *** E/B(THETA) ****

       DO NR=1,NX1
          GY(NR,1)=GUCLIP(DBLE(CF(NR,2)))
          GY(NR,2)=GUCLIP(DIMAG(CF(NR,2)))
          GY(NR,3)=GUCLIP(ABS(CF(NR,2)))
       ENDDO
       CALL wm_gsub(NX1,GX1,GXMIN,GXMAX,GY,2,GP(1,2),KTITL(2))

!      *** E/B(Z) ****

       DO NR=1,NX1
          GY(NR,1)=GUCLIP(DBLE(CF(NR,3)))
          GY(NR,2)=GUCLIP(DIMAG(CF(NR,3)))
          GY(NR,3)=GUCLIP(ABS(CF(NR,3)))
       ENDDO
       CALL wm_gsub(NX1,GX1,GXMIN,GXMAX,GY,2,GP(1,3),KTITL(3))

!        *** POWER ***

       DO I=1,NG4
          DO NR=1,NX2
             GY(NR,I)=GUCLIP(POWER(NR,I))
          ENDDO
       ENDDO
       CALL wm_gsub(NX2,GX2,GXMIN,GXMAX,GY,NG4,GP(1,4),KTITL(4))
    ENDIF

    CALL SETLIN(0,0,7)
    IF(K3.EQ.'A')  THEN
       CALL MOVE(GP(2,4)-2.0,GP(4,4)+0.2)
       CALL TEXT('total',5)
    END IF

    IF(K3.EQ.'M')  THEN
       CALL MOVE(GP(2,4)-2.0,GP(4,4)+0.2)
       CALL TEXT('mode ',5)
    END IF

    IF(K3.EQ.'N')  THEN
       CALL MOVE(GP(2,4)-6.0,GP(4,4)+0.2)
       CALL TEXT('md:',3)
       CALL NUMBI(NTH0+MD1,'(I3)',3)
       DO I=1,IMAX
          GXL=GP(2,4)-4.0+2.0*(I-1)/REAL(IMAX-1)
          CALL SETLIN(0,2,7-MOD(I-1,5))
          CALL MOVE(GXL,GP(4,4)+0.4)
          CALL DRAW(GXL,GP(4,4)+0.2)
       ENDDO
       CALL SETLIN(0,0,7)
       CALL MOVE(GP(2,4)-2.0,GP(4,4)+0.2)
       CALL NUMBI(NTH0+MD2,'(I3)',3)
    END IF

    CALL wm_gprm('R',K3,NTH,NHH,MD,ND)

    CALL PAGEE

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_gr1d

!     ****** DRAW LINES OF 1D GRAPH ******

  SUBROUTINE wm_gn1d(NX,GX,GXMIN,GXMAX,GY1,GY2,NY,GP,KTITL)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NX,NY
    REAL,INTENT(IN):: GX(NX),GXMIN,GXMAX,GY1(NX,NY),GY2(NX,NY),GP(4)
    CHARACTER(LEN=6),INTENT(IN):: KTITL
    INTEGER:: I
    REAL:: GYMIN,GYMAX,GYMIN1,GYMAX1,GSX,GSYMIN,GSYMAX,GSY
    REAL:: GYMIN2,GYMAX2
    EXTERNAL GMNMX1,GQSCAL,GDEFIN,GFRAME,GSCALE,GVALUE,SETLIN,GPLOTP
    EXTERNAL MOVE,TEXT

    IF(NY.EQ.0) RETURN

    GYMIN= 1.E32
    GYMAX=-1.E32
    DO I=1,NY
       CALL GMNMX1(GY1(1,I),1,NX,1,GYMIN1,GYMAX1)
       CALL GMNMX1(GY2(1,I),1,NX,1,GYMIN2,GYMAX2)
       GYMIN=MIN(GYMIN,GYMIN1,GYMIN2)
       GYMAX=MAX(GYMAX,GYMAX1,GYMAX2)
    ENDDO

!    CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)
    GSX=0.1

    IF(GYMIN.GT.0.0.AND.GYMAX.GT.0.0) THEN
       GYMIN=0.0
    ELSE IF(GYMIN.LT.0.0.AND.GYMAX.LT.0.0) THEN
       GYMAX=0.0
    ENDIF
    CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)

    IF(ABS(GSYMAX-GSYMIN).LT.1.E-15) GOTO 9000

    CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),GXMIN,GXMAX,GSYMIN,GSYMAX)
    CALL GFRAME
    CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
    CALL GVALUE(GXMIN,GSX*5,0.0,0.0,NGULEN(GSX*5))
    CALL GSCALE(0.0,0.0,0.0,GSY,0.1,9)
    CALL GVALUE(0.0,0.0,0.0,GSY*2,NGULEN(GSY*2))

    DO I=1,NY
       CALL SETLIN(0,2,7-MOD(I-1,5))
       CALL GPLOTP(GX,GY1(1,I),1,NX,1,0,0,0)
       CALL GPLOTP(GX,GY2(1,I),1,NX,1,0,0,2)
    ENDDO
    CALL SETLIN(0,2,7)

9000 CONTINUE
    CALL MOVE(GP(1),GP(4)+0.2)
    CALL TEXT(KTITL,6)

    RETURN
  END SUBROUTINE wm_gn1d

! ****** DRAW 2D POLOIDAL GRAPH ******

  SUBROUTINE wm_grth(K2,K3,K4)

    USE wmcomm
    USE wmgsub
    IMPLICIT NONE
    CHARACTER(LEN=1),INTENT(IN):: K2,K3,K4
    REAL,ALLOCATABLE:: GXR(:),GY(:,:)
    INTEGER:: NA3,NG3,NHH,NTH,NR,NX,I,NSTEP
    REAL:: GRMIN,GRMAX,GZ1MIN,GZ1MAX,GZMIN,GZMAX,GZSTEP
    COMPLEX(rkind):: CFL
    EXTERNAL PAGES,SETCHS,GMNMX1,GMNMX2,GQSCAL,GDEFIN,SETLIN
    EXTERNAL CONTQ3,MOVE,TEXT,NUMBR,PAGEE

    ALLOCATE(GXR(nrmax+1),GY(nrmax+1,nthmax*nhhmax))

    IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
       SELECT CASE(K3)
       CASE('R')
          NA3=1
          NG3=1
       CASE('T')
          NA3=1
          NG3=2
       CASE('Z')
          NA3=1
          NG3=3
       CASE('S')
          NA3=2
          NG3=1
       CASE('H')
          NA3=2
          NG3=2
       CASE('B')
          NA3=2
          NG3=3
       CASE('+')
          NA3=3
          NG3=1
       CASE('-')
          NA3=3
          NG3=2
       CASE('P')
          NA3=3
          NG3=3
       CASE DEFAULT
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRTH'
          GOTO 9000
       END SELECT
    ELSEIF(K2.EQ.'P') THEN
       SELECT CASE(K3)
       CASE('E')
          NG3=1
       CASE('D')
          NG3=2
       CASE('T')
          NG3=3
       CASE('1')
            NG3=1
       CASE('2')
          NG3=2
       CASE('3')
          NG3=3
       CASE('4')
          NG3=4
       CASE('5')
          NG3=5
       CASE('6')
          NG3=6
       CASE DEFAULT
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRTH'
          GOTO 9000
       END SELECT
    ELSE IF(K2.EQ.'J') THEN
       NG3=0
    ELSE
       WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #2 IN WMGRTH'
       GOTO 9000
    ENDIF

1   CONTINUE
    IF(NHHMAX.EQ.1) THEN
       NHH=1
    ELSE
       WRITE(6,*) '## INPUT NHH : 1..',NHHMAX
       READ(5,*,ERR=1,END=9000) NHH
    END IF
    IF(NHH.LT.1.OR.NHH.GT.NHHMAX) THEN
       WRITE(6,*) 'XX ILLEGAL NHH'
       GOTO 1
    ENDIF

    DO NTH=1,MDSIZ
       DO NR=1,NRMAX+1
          IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
             IF(K2.EQ.'E') THEN
                IF(NA3.EQ.1) THEN
                   CFL=CEN(NG3,NTH,NHH,NR)
                ELSEIF(NA3.EQ.2) THEN
                   CFL=CEB(NG3,NTH,NHH,NR)
                ELSEIF(NA3.EQ.3) THEN
                   CFL=CEP(NG3,NTH,NHH,NR)
                ENDIF
             ELSEIF(K2.EQ.'B') THEN
                CFL=CBFLD(NG3,NTH,NHH,NR)
             ENDIF

             IF(K4.EQ.'R') THEn
                GY(NR,NTH)=GUCLIP(DBLE(CFL))
             ELSEIF(K4.EQ.'I') THEN
                GY(NR,NTH)=GUCLIP(DIMAG(CFL))
             ELSEIF(K4.EQ.'A') THEN
                GY(NR,NTH)=GUCLIP(ABS(CFL))
             ELSE
                WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRTH'
                GOTO 9000
             ENDIF
             NX=NRMAX+1
          ELSEIF(K2.EQ.'P') THEN
             GY(NR,NTH)=GUCLIP(PABS(NTH,NHH,NR,NG3))
          ELSEIF(K2.EQ.'J') THEN
             GY(NR,NTH)=GUCLIP(PCUR(NTH,NHH,NR))
          ENDIF
       ENDDO
    ENDDO
    IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
       NX=NRMAX+1
    ELSEIF(K2.EQ.'P'.OR.K2.EQ.'J') THEN
       NX=NRMAX
    ENDIF

    CALL PAGES
    CALL SETCHS(0.3,0.0)

    DO I=1,NX
       GXR(I)=GUCLIP(XR(I))
    ENDDO

    CALL GMNMX1(GXR,1,NX,1,GRMIN,GRMAX)
    CALL GMNMX2(GY,nrmax+1,1,NX,1,1,MDSIZ,1,GZ1MIN,GZ1MAX)
    CALL GQSCAL(GZ1MIN,GZ1MAX,GZMIN,GZMAX,GZSTEP)
    NSTEP=NCONT
    GZSTEP=(GZMAX-GZMIN)/NSTEP

    CALL GDEFIN(2.5,17.5,2.5,17.5,-GRMAX,GRMAX,-GRMAX,GRMAX)
    CALL SETLIN(0,2,4)
    CALL wm_circle(GRMAX)
    IF(GZMIN.GT.0.0) THEN
       CALL SETLIN(0,2,6)
       CALL CONTQ3(GY,GXR,nrmax+1,NX,MDSIZ,GZMIN,GZSTEP,NSTEP,0,KACONT)
    ELSEIF(GZMAX.LT.0.0) THEN
       CALL SETLIN(0,2,5)
       CALL CONTQ3(GY,GXR,nrmax+1,NX,MDSIZ,GZMIN,GZSTEP,NSTEP,2,KACONT)
    ELSE
       CALL SETLIN(0,2,6)
       CALL CONTQ3(GY,GXR,nrmax+1,NX,MDSIZ,0.5*GZSTEP,GZSTEP,NSTEP,0,KACONT)
       CALL SETLIN(0,2,5)
       CALL CONTQ3(GY,GXR,nrmax+1,NX,MDSIZ, &
                                         -0.5*GZSTEP,-GZSTEP,NSTEP,2,KACONT)
    END IF
    CALL SETLIN(0,0,7)

    CALL MOVE(16.5,17.0)
    SELECT CASE(K2)
    CASE('E')
       CALL TEXT('E ',2)
    CASE('B')
       CALL TEXT('B ',2)
    CASE('P')
       CALL TEXT('P ',2)
    CASE('C')
       CALL TEXT('J ',2)
    END SELECT
    IF(K2.NE.'P') THEN
       SELECT CASE(K3)
       CASE('R')
          CALL TEXT('r ',2)
       CASE('T')
          CALL TEXT('theta ',6)
       CASE('Z')
          CALL TEXT('phi ',4)
       CASE('S')
          CALL TEXT('s ',2)
       CASE('H')
          CALL TEXT('h ',2)
       CASE('B')
          CALL TEXT('b ',2)
       CASE('+')
          CALL TEXT('+ ',2)
       CASE('-')
          CALL TEXT('- ',2)
       CASE('P')
          CALL TEXT('P ',2)
       END SELECT
    ELSE
       SELECT CASE(K3)
       CASE('E')
          CALL TEXT('e ',2)
       CASE('D')
          CALL TEXT('D ',2)
       CASE('T')
          CALL TEXT('T ',2)
       CASE('1')
          CALL TEXT('1 ',2)
       CASE('2')
          CALL TEXT('2 ',2)
       CASE('3')
          CALL TEXT('3 ',2)
       CASE('4')
          CALL TEXT('4 ',2)
       CASE('5')
          CALL TEXT('5 ',2)
       CASE('6')
          CALL TEXT('6 ',2)
       END SELECT
    ENDIF
    SELECT CASE(K4)
    CASE('R')
       CALL TEXT('Real ',5)
    CASE('I')
       CALL TEXT('Imag ',5)
    CASE('A')
       CALL TEXT('Abs ',4)
    END SELECT

    CALL MOVE(16.5,16.0)
    CALL TEXT('MAX :',5)
    CALL NUMBR(GZ1MAX,'(1PE12.4)',12)
    CALL MOVE(16.5,15.5)
    CALL TEXT('MIN :',5)
    CALL NUMBR(GZ1MIN,'(1PE12.4)',12)
    CALL MOVE(16.5,15.0)
    CALL TEXT('STEP:',5)
    CALL NUMBR(GZSTEP,'(1PE12.4)',12)

    CALL wm_gprm('C',K3,0,0,0,0)

    CALL PAGEE

    RETURN
9000 CONTINUE
    RETURN
  END SUBROUTINE wm_grth

!     ****** DRAW 2D MODE GRAPH ******

  SUBROUTINE wm_grmd(K2,K3,K4)

    USE wmcomm
    USE wmgsub
    IMPLICIT NONE
    CHARACTER(LEN=1):: K2,K3,K4
    REAL,ALLOCATABLE:: GY(:,:)
    INTEGER:: NG3,ND,NDX,NR,MDX,NX
    EXTERNAL PAGES,SETCHS,PAGEE

    ALLOCATE(GY(nrmax+1,nthmax+1))

    IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
       SELECT CASE(K3)
       CASE('R')
          NG3=1
       CASE('T')
          NG3=2
       CASE('Z')
          NG3=3
       CASE DEFAULT
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRMD'
          GOTO 9000
       END SELECT
    ELSEIF(K2.EQ.'P') THEN
       SELECT CASE(K3)
       CASE('E')
          NG3=1
       CASE('H')
          NG3=2
       CASE('D')
          NG3=3
       CASE('1')
          NG3=1
       CASE('2')
          NG3=2
       CASE('3')
          NG3=3
       CASE('4')
          NG3=4
       CASE('5')
          NG3=5
       CASE('6')
          NG3=6
       CASE DEFAULT
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRMD'
          GOTO 9000
       END SELECT
    ELSE IF(K2.EQ.'J') THEN
       NG3=0
    ELSE
       WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #2 IN WMGRMD'
       GOTO 9000
    ENDIF

1   CONTINUE
    IF(NHHMAX.EQ.1) THEN
       ND=0
    ELSE
       WRITE(6,*) '## INPUT ND : ',NDMIN,'..',NDMAX-1
       READ(5,*,ERR=1,END=9000) ND
    END IF
    IF(ND.LT.NDMIN.OR.ND.GT.NDMAX) THEN
       WRITE(6,*) 'XX ILLEGAL ND'
       GOTO 1
    ENDIF
    NDX=ND-NDMIN+1

    SELECT CASE(K2)
    CASE('E')
       SELECT CASE(K4)
       CASE('R')
          DO NR=1,NRMAX+1
             DO MDX=1,MDSIZ
                GY(NR,MDX+1)=GUCLIP(DBLE(CEFLDK(NG3,MDX,NDX,NR)))
             ENDDO
             GY(NR,1)=GUCLIP(DBLE(CEFLDK(NG3,MDSIZ,NDX,NR)))
          ENDDO
       CASE('I')
          DO NR=1,NRMAX+1
             DO MDX=1,MDSIZ
                GY(NR,MDX+1)=GUCLIP(DIMAG(CEFLDK(NG3,MDX,NDX,NR)))
             ENDDO
             GY(NR,1)=GUCLIP(DIMAG(CEFLDK(NG3,MDSIZ,NDX,NR)))
          ENDDO
       CASE('A')
          DO NR=1,NRMAX+1
             DO MDX=1,MDSIZ
                GY(NR,MDX+1)=GUCLIP(ABS(CEFLDK(NG3,MDX,NDX,NR)))
             ENDDO
             GY(NR,1)=GUCLIP(ABS(CEFLDK(NG3,MDSIZ,NDX,NR)))
          ENDDO
       CASE DEFAULT
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRMD'
          GOTO 9000
       END SELECT
       NX=NRMAX+1
    CASE('B')
       SELECT CASE(K4)
       CASE('R')
          DO NR=1,NRMAX+1
             DO MDX=1,MDSIZ
                GY(NR,MDX+1)=GUCLIP(DBLE(CBFLDK(NG3,MDX,NDX,NR)))
             ENDDO
             GY(NR,1)=GUCLIP(DBLE(CBFLDK(NG3,MDSIZ,NDX,NR)))
          ENDDO
       CASE('I')
          DO NR=1,NRMAX+1
             DO MDX=1,MDSIZ
                GY(NR,MDX+1)=GUCLIP(DIMAG(CBFLDK(NG3,MDX,NDX,NR)))
             ENDDO
             GY(NR,1)=GUCLIP(DIMAG(CBFLDK(NG3,MDSIZ,NDX,NR)))
          ENDDO
       CASE('A')
          DO NR=1,NRMAX+1
             DO MDX=1,MDSIZ
                GY(NR,MDX+1)=GUCLIP(ABS(CBFLDK(NG3,MDX,NDX,NR)))
             ENDDO
             GY(NR,1)=GUCLIP(ABS(CBFLDK(NG3,MDSIZ,NDX,NR)))
          ENDDO
       CASE DEFAULT
          WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRMD'
          GOTO 9000
       END SELECT
       NX=NRMAX+1
    CASE('P')
       DO NR=1,NRMAX
          DO MDX=1,MDSIZ
             GY(NR,MDX+1)=GUCLIP(PABSK(MDX,NDX,NR,NG3))
          ENDDO
          GY(NR,1)=GUCLIP(PABSK(MDSIZ,NDX,NR,NG3))
       ENDDO
       NX=NRMAX
    END SELECT

!
    CALL PAGES

    CALL SETCHS(0.3,0.0)

    CALL wm_gcon(GY,K2,K3,K4,NX)

    CALL wm_gprm('M',K3,0,0,0,0)

    CALL PAGEE

9000 CONTINUE
    RETURN
  END SUBROUTINE wm_grmd

!     ****** DRAW GRAPH OF MULTIPLE LINES ******

  SUBROUTINE wm_gcon(GZL,K2,K3,K4,NX)

    USE wmcomm
    IMPLICIT NONE
    REAL,INTENT(IN):: GZL(nrmax+1,nthmax+1)
    INTEGER,INTENT(IN):: NX
    CHARACTER(LEN=1),INTENT(IN):: K2,K3,K4
    REAL,ALLOCATABLE:: GY(:)
    INTEGER:: MD,MDX,NR,NRLCFS,NSTEP
    REAL:: GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX,GZSTEP
    REAL:: GSXMIN,GSXMAX,GSX,GSYMIN,GSYMAX,GSY,GSZMIN,GSZMAX,GSZ
    EXTERNAL GMNMX2,GQSCAL,GDEFIN,GFRAME,GSCALE,GVALUE,SETLIN,CONTQ1
    EXTERNAL MOVE,TEXT,NUMBR
    
    ALLOCATE(GY(nthmax+1))

    GXMIN=GUCLIP(XR(1))
    GXMAX=GUCLIP(XR(NRMAX+1))

    DO MD=MDMIN-1,MDMAX
       MDX=MD-MDMIN+2
       GY(MDX)=NTH0+MD
    ENDDO

    GYMIN=GY(1)
    GYMAX=GY(MDSIZ+1)

    CALL GMNMX2(GZL,nrmax+1,1,NX,1,1,MDSIZ+1,1,GZMIN,GZMAX)

    ! --- evaluate min and max only inside of plasma ---
    DO NR=1,NX
       IF(xrho(NR).LT.1.D0 ) NRLCFS=NR
    ENDDO
    NRLCFS=NRLCFS+1

    CALL GMNMX2(GZL,nrmax+1,1,NRLCFS,1,1,nthmax+1,1,GZMIN,GZMAX)

    CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)
    CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)
    CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSZ)

    NSTEP=NCONT
    GZSTEP=(GZMAX-GZMIN)/NSTEP

    CALL GDEFIN(3.0,21.6,1.0,16.5,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL GFRAME
    CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
    CALL GVALUE(GXMIN,GSX*2,0.0,0.0,1)
    CALL GSCALE(0.0,0.0,0.0,MAX(1.0,GSY),0.1,9)
    CALL GVALUE(0.0,0.0,0.0,MAX(1.0,GSY),0)

    IF (GZMIN.GT.0.0) THEN
       CALL SETLIN(0,2,6)
       CALL CONTQ1(GZL,nrmax+1,NX,nthmax+1,GZMIN,GZSTEP,NSTEP,0,0,KACONT)
    ELSE IF (GZMAX.LT.0.0) THEN
       CALL SETLIN(0,2,5)
       CALL CONTQ1(GZL,nrmax+1,NX,nthmax+1+1,GZMIN,GZSTEP,NSTEP,0,2,KACONT)
    ELSE 
       CALL SETLIN(0,2,6)
       CALL CONTQ1(GZL,nrmax+1,NX,nthmax+1, &
                       0.5*GZSTEP, GZSTEP,NSTEP,0,0,KACONT)
       CALL SETLIN(0,2,5)
       CALL CONTQ1(GZL,nrmax+1,NX,nthmax+1, &
                      -0.5*GZSTEP,-GZSTEP,NSTEP,0,2,KACONT)
    END IF
    CALL SETLIN(0,0,7)

    CALL MOVE(1.0,17.0)
    CALL TEXT('mode',4)
    CALL MOVE(4.5,17.3)
    SELECT CASE(K2)
    CASE('E')
       CALL TEXT('E ',2)
    CASE('B')
       CALL TEXT('B ',2)
    CASE('P')
       CALL TEXT('P ',2)
    CASE('C')
       CALL TEXT('J ',2)
    END SELECT
    IF(K2.NE.'P') THEN
       SELECT CASE(K3)
       CASE('R')
          CALL TEXT('  r   ',6)
       CASE('T')
          CALL TEXT('theta ',6)
       CASE('Z')
          CALL TEXT('  z   ',6)
       END SELECT
    ELSE
       SELECT CASE(K3)
       CASE('E')
          CALL TEXT('e ',2)
       CASE('D')
          CALL TEXT('D ',2)
       CASE('T')
          CALL TEXT('T ',2)
       END SELECT
    END IF
    SELECT CASE(K4)
    CASE('R')
       CALL TEXT('Real',4)
    CASE('I')
       CALL TEXT('Imag',4)
    CASE('A')
       CALL TEXT('Abs ',4)
    END SELECT
    
    CALL MOVE(10.0,17.3)
    CALL TEXT('min=',4)
    CALL NUMBR(GZMIN,'(1PE12.4)',12)
    CALL MOVE(10.0,16.8)
    CALL TEXT('max=',4)
    CALL NUMBR(GZMAX,'(1PE12.4)',12)
    CALL MOVE(16.0,17.3)
    CALL TEXT('step=',5)
    CALL NUMBR(GZSTEP,'(1PE12.4)',12)

    RETURN
  END SUBROUTINE wm_gcon

!     ****** EXPLANATION OF GRAPHIC COMMAND ******

  SUBROUTINE wm_ghelp
    IMPLICIT NONE
    WRITE(6,600)
600 FORMAT( &
         '   R: radial profile'/ &
         '    E: wave electric field'/ &
         '    B: wave magnetic field'/ &
         '     A: wave field at a poloidal angle and total Pabs'/ &
         '     T: wave field and Pabs at a poloidal angle'/ &
         '     M: wave field and Pabs with a poloidal mode number m'/ &
         '     N: wave field and Pabs for a range of m'/) 

    WRITE(6,601)
601 FORMAT( &
         '   P: poloidal cross section'/ &
         '   M: poloidal mode spectrum'/ &
         '   C: circular poloidal projection'/ &
         '    E: wave electric field'/ &
         '    B: wave magnetic field'/ &
         '     R: radial component'/ &
         '     T: poloidal component'/ &
         '     Z: toroidal component'/ &
         '     s: radial component'/ &
         '     h: perpendicular component'/ &
         '     b: parallel component'/ &
         '     +: right hand polarized component'/ &
         '     -: left hand polarized component'/ &
         '     P: parallel component'/ &
         '      R: real component'/ &
         '      I: imaginary component'/ &
         '      A: absolute value'/)

    WRITE(6,602)
602 FORMAT( &
         '   P: poloidal cross section'/ &
         '   M: poloidal mode spectrum'/ &
         '   C: circular poloidal projection'/ &
         '    P: absorbed power density'/ &
         '     n: particle species number'/ &
         '        (usually 1: electron)'/)

    WRITE(6,603)
603 FORMAT( &
         '   P: poloidal cross section'/ &
         '   C: circular poloidal projection'/ &
         '    J: currnet density'// &
         '     n: particle species number'/ &
         '        (usually 1: electron)'// &
         '   P: poloidal cross section'/ &
         '    F: various quantities'/ &
         '     S: psi'/ &
         '     B: total magnetic field'/ &
         '     Q: safety factor'/ &
         '     2: poloidal magnetic field'/ &
         '     3: toroidal magnetic field'/ &
         '     J: Jacobian'/)

    WRITE(6,604)
604 FORMAT( &
         '   R: radial profile'/ &
         '    G: Fourier components metric tensor'// &
         '   S: equilibrium profile'// &
         '   G: graphic type'/ &
         '    1: contour map'/ &
         '    2: painted map'/ &
         '    3: bird eye view 1 (only for circular for the present)'/ &
         '    4: bird eye view 2 (only for circular for the present)'// &
         '   ?: this help message'/ &
         '   X: exit'/)
    RETURN
  END SUBROUTINE wm_ghelp
END MODULE wmgout
