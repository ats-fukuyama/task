C     $Id$
C
C     ****** CONTROL GRAPHICS ******
C
      SUBROUTINE WMGOUT
C
      INCLUDE 'wmcomm.h'
C
      CHARACTER KSTR*5,K1,K2,K3,K4
C
    1 WRITE(6,*) ' ## INPUT GSTR : R/EB/ATM  RCMP/P/123  CP/J  P/G  S'
      WRITE(6,*) '                 CMP/EB/RTZ/RIA  P/F/SBQ23J  X/EXIT'
      READ(5,'(A5)',ERR=1,END=900) KSTR
      K1=KSTR(1:1)
      CALL GUCPTL(K1)
      IF (K1.EQ.'X') GOTO 900
      IF (K1.EQ.'G') THEN
	 K2=KSTR(2:2)
	 CALL GUCPTL(K2)
         IF(K2.EQ.'1') NGRAPH=1
         IF(K2.EQ.'2') NGRAPH=2
         IF(K2.EQ.'3') NGRAPH=3
         IF(K2.EQ.'4') NGRAPH=4
         GOTO 1
      ENDIF
      IF (K1.EQ.'Q') GOTO 900
C
      IF((K1.EQ.'R').OR.(K1.EQ.'C').OR.(K1.EQ.'M').OR.
     &   (K1.EQ.'P').OR.(K1.EQ.'S')) THEN
	 K2=KSTR(2:2)
	 CALL GUCPTL(K2)
         K3=KSTR(3:3)
	 CALL GUCPTL(K3)
         K4=KSTR(4:4)
	 CALL GUCPTL(K4)
         IF(K1.EQ.'R') CALL WMGR1D(K2,K3)
         IF(K1.EQ.'C') CALL WMGRTH(K2,K3,K4)
         IF(K1.EQ.'M') CALL WMGRMD(K2,K3,K4)
         IF(K1.EQ.'P') CALL WMGREQ(K2,K3,K4)
         IF(K1.EQ.'S'.AND.MODELG.GE.4) CALL WMGRMS
      ELSE
         WRITE(6,*) '## UNDEFINED CONTROL CHARACTER: K1=',K1
      END IF
      GOTO 1
C
  900 RETURN
      END
C
C     ****** DRAW 1D GRAPH ******
C
      SUBROUTINE WMGR1D(K2,K3)
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION GX1(NRM),GX2(NRM),GY(NRM,NSM),GP(4,4)
      DIMENSION GY1(NRM,MDM),GY2(NRM,MDM)
      DIMENSION CF(NRM,3),POWER(NRM,NSM)
      DIMENSION CF1(NRM,MDM),CF2(NRM,MDM),CF3(NRM,MDM)
      DIMENSION PF1(NRM,MDM),PF2(NRM,MDM)
C
      CHARACTER*6 KTITL(4)
      CHARACTER K2,K3
C
      DATA GP/ 3.0, 10.8,  9.5, 16.5,
     &         3.0, 10.8,  1.0,  8.0,
     &        13.8, 21.6,  9.5, 16.5,
     &        13.8, 21.6,  1.0,  8.0/
C
      GXMIN=GCLIP(XR(1))
      GXMAX=GCLIP(XR(NRMAX+1))
C
      IF(K2.EQ.'E') THEN
         NR=1
            GX1(NR)=GCLIP(       XR(NR))
            GX2(NR)=GCLIP(       XR(NR))
         DO 10 NR=2,NRMAX+1
            GX1(NR)=GCLIP(       XR(NR))
            GX2(NR)=GCLIP(0.5D0*(XR(NR-1)+XR(NR)))
   10    CONTINUE
         NX1=NRMAX+1
         NX2=NRMAX+1
         NG4=NSMAX
         KTITL(1)='Er    '
         KTITL(2)='Etheta'
         KTITL(3)='Ez    '
         KTITL(4)='Pabs  '
      ELSE
         NR=1
            GX2(NR)=GCLIP(       XR(NR))
            GX1(NR)=GCLIP(       XR(NR))
         DO 20 NR=2,NRMAX+1
            GX2(NR)=GCLIP(       XR(NR))
            GX1(NR)=GCLIP(0.5D0*(XR(NR-1)+XR(NR)))
   20    CONTINUE
         NX1=NRMAX+1
         NX2=NRMAX+1
         NG4=1
         KTITL(1)='Br    '
         KTITL(2)='Btheta'
         KTITL(3)='Bz    '
         KTITL(4)='Jrf   '
      ENDIF
C
      IF(K3.EQ.'T'.OR.K3.EQ.'A') THEN
    1    IF(NTHMAX.EQ.1) THEN
            NTH=1
         ELSE
            WRITE(6,*) '## INPUT NTH : 1..',NTHMAX
            READ(5,*,ERR=1,END=9000) NTH
         END IF
         IF(NTH.LT.1.OR.NTH.GT.NTHMAX) THEN
            WRITE(6,*) 'XX ILLEGAL NTH'
            GOTO 1
         ENDIF
    2    IF(NPHMAX.EQ.1) THEN
            NPH=1
         ELSE
            WRITE(6,*) '## INPUT NPH : 1..',NPHMAX
            READ(5,*,ERR=2,END=9000) NPH
         END IF
         IF(NPH.LT.1.OR.NPH.GT.NPHMAX) THEN
            WRITE(6,*) 'XX ILLEGAL NPH'
            GOTO 2
         ENDIF
C
         IF(K2.EQ.'E') THEN
            DO 100 I=1,3
            DO 100 NR=1,NRMAX+1
               CF(NR,I)=CEFLD(I,NTH,NPH,NR)
  100       CONTINUE
            IF(K3.EQ.'A') THEN
               DO NS=1,NSMAX
                  DO NR=1,NRMAX
                     POWER(NR,NS)=PABSR(NR,NS)
                  ENDDO
                  POWER(NRMAX+1,NS)=0.D0
               ENDDO
            ELSE
               DO NS=1,NSMAX
                  DO NR=1,NRMAX
                     POWER(NR,NS)=PABS(NTH,NPH,NR,NS)
                  ENDDO
                  POWER(NRMAX+1,NS)=0.D0
               ENDDO
            ENDIF
         ELSE
            DO 200 I=1,3
            DO 200 NR=1,NRMAX+1
               CF(NR,I)=CBFLD(I,NTH,NPH,NR)
  200       CONTINUE
            IF(K3.EQ.'A') THEN
               DO 210 NR=1,NRMAX
C                  POWER(NR,1)=PCURR(NR)
                  POWER(NR,1)=QPS(NR)
  210          CONTINUE
               POWER(NRMAX+1,1)=0.D0
            ELSE
               DO 220 NR=1,NRMAX
                  POWER(NR,1)=PCUR(NTH,NPH,NR)
  220          CONTINUE
               POWER(NRMAX+1,1)=0.D0
            ENDIF
         ENDIF
      ELSEIF(K3.EQ.'M') THEN
    3    IF(NTHMAX.EQ.1) THEN
            MD=0
         ELSE
            WRITE(6,*) '## INPUT MD : ',MDMIN,'..',MDMAX-1
            READ(5,*,ERR=3,END=9000) MD
         END IF
         IF(MD.LT.MDMIN.OR.MD.GT.MDMAX) THEN
            WRITE(6,*) 'XX ILLEGAL MD'
            GOTO 3
         ENDIF
         MDX=MD-MDMIN+1
    4    IF(NPHMAX.EQ.1) THEN
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
C
         IF(K2.EQ.'E') THEN
            DO 300 I=1,3
            DO 300 NR=1,NRMAX+1
               CF(NR,I)=CEFLDK(I,MDX,NDX,NR)
  300       CONTINUE
            DO NS=1,NSMAX
               DO NR=1,NRMAX
                  POWER(NR,NS)=PABSK(MDX,NDX,NR,NS)
               ENDDO
               POWER(NRMAX+1,NS)=0.D0
            ENDDO
         ELSE
            DO 400 I=1,3
            DO 400 NR=1,NRMAX+1
               CF(NR,I)=CBFLDK(I,MDX,NDX,NR)
  400       CONTINUE
            DO 410 NR=1,NRMAX
               POWER(NR,1)=PCUR(MDX,NDX,NR)
  410       CONTINUE
               POWER(NRMAX+1,1)=0.D0
         ENDIF
      ELSEIF(K3.EQ.'N') THEN
    7    IF(NTHMAX.EQ.1) THEN
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
    8    IF(NPHMAX.EQ.1) THEN
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
C
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
      ENDIF
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
      IF(K3.NE.'N') THEN
C
C        *** E/B(R) ****
C
         DO NR=1,NX2
            GY(NR,1)=GCLIP(DBLE(CF(NR,1)))
            GY(NR,2)=GCLIP(DIMAG(CF(NR,1)))
            GY(NR,3)=GCLIP(ABS(CF(NR,1)))
         ENDDO
         CALL WMGSUB(NX1,GX1,GXMIN,GXMAX,GY,2,GP(1,1),KTITL(1))
C
C        *** E/B(THETA) ****
C
         DO NR=1,NX1
            GY(NR,1)=GCLIP(DBLE(CF(NR,2)))
            GY(NR,2)=GCLIP(DIMAG(CF(NR,2)))
            GY(NR,3)=GCLIP(ABS(CF(NR,2)))
         ENDDO
         CALL WMGSUB(NX1,GX1,GXMIN,GXMAX,GY,2,GP(1,2),KTITL(2))
C
C        *** E/B(Z) ****
C
         DO NR=1,NX1
            GY(NR,1)=GCLIP(DBLE(CF(NR,3)))
            GY(NR,2)=GCLIP(DIMAG(CF(NR,3)))
            GY(NR,3)=GCLIP(ABS(CF(NR,3)))
         ENDDO
         CALL WMGSUB(NX1,GX1,GXMIN,GXMAX,GY,2,GP(1,3),KTITL(3))
C
C        *** POWER / CURRENT ***
C
         DO I=1,NG4
         DO NR=1,NX2
            GY(NR,I)=GCLIP(POWER(NR,I))
         ENDDO
         ENDDO
         CALL WMGSUB(NX2,GX2,GXMIN,GXMAX,GY,NG4,GP(1,4),KTITL(4))
      ELSE
C
C        *** E/B(R) ****
C
         DO I=1,IMAX
         DO NR=1,NX2
            GY1(NR,I)=GCLIP(DBLE(CF1(NR,I)))
            GY2(NR,I)=GCLIP(DIMAG(CF1(NR,I)))
         ENDDO
         ENDDO
         CALL WMGN1D(NX1,GX1,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,1),KTITL(1))
C
C        *** E/B(THETA) ****
C
         DO I=1,IMAX
         DO NR=1,NX1
            GY1(NR,I)=GCLIP(DBLE(CF2(NR,I)))
            GY2(NR,I)=GCLIP(DIMAG(CF2(NR,I)))
         ENDDO
         ENDDO
         CALL WMGN1D(NX1,GX1,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,2),KTITL(2))
C
C        *** E/B(Z) ****
C
         DO I=1,IMAX
         DO NR=1,NX1
            GY1(NR,I)=GCLIP(DBLE(CF3(NR,I)))
            GY2(NR,I)=GCLIP(DIMAG(CF3(NR,I)))
         ENDDO
         ENDDO
         CALL WMGN1D(NX1,GX1,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,3),KTITL(3))
C
C        *** POWER / CURRENT ***
C
         DO I=1,IMAX
         DO NR=1,NX2
            GY1(NR,I)=GCLIP(PF1(NR,I))
            GY2(NR,I)=GCLIP(PF2(NR,I))
         ENDDO
         ENDDO
         CALL WMGN1D(NX2,GX2,GXMIN,GXMAX,GY1,GY2,IMAX,GP(1,4),KTITL(4))
      ENDIF
C
      CALL SETLIN(0,0,7)
      IF(K3.EQ.'A')  THEN
        CALL MOVE(GP(2,4)-2.0,GP(4,4)+0.2)
        CALL TEXT('total',5)
      END IF
C
      IF(K3.EQ.'M')  THEN
        CALL MOVE(GP(2,4)-2.0,GP(4,4)+0.2)
        CALL TEXT('mode ',5)
      END IF
C
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
C
      CALL WMGPRM('R',K3,NTH,NPH,MD,ND)
C
      CALL PAGEE
C
 9000 RETURN
      END
C
C     ****** DRAW LINES OF 1D GRAPH ******
C
      SUBROUTINE WMGSUB(NX,GX,GXMIN,GXMAX,GY,NY,GP,KTITL)
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION GX(NRM),GY(NRM,NSM),GP(4)
      DIMENSION ICL(6),IPAT(6)
      CHARACTER KTITL*6
      DATA ICL/7,6,5,4,3,2/,IPAT/0,2,4,6,3,1/
C
      IF(NY.EQ.0) RETURN
C
      GYMIN= 1.E32
      GYMAX=-1.E32
      DO 100 I=1,NY
         CALL GMNMX1(GY(1,I),1,NX,1,GYMIN1,GYMAX1)
         GYMIN=MIN(GYMIN,GYMIN1)
         GYMAX=MAX(GYMAX,GYMAX1)
  100 CONTINUE
C
C     CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)
      GSX=0.1
C
      IF(GYMIN.GT.0.0.AND.GYMAX.GT.0.0) THEN
         GYMIN=0.0
      ELSE IF(GYMIN.LT.0.0.AND.GYMAX.LT.0.0) THEN
         GYMAX=0.0
      ENDIF
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)
C
      IF(ABS(GSYMAX-GSYMIN).LT.1.E-15) GOTO 9000
C
      CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),GXMIN,GXMAX,GSYMIN,GSYMAX)
      CALL GFRAME
      CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
      CALL GVALUE(GXMIN,GSX*5,0.0,0.0,NGVLEN(GSX*5))
      CALL GSCALE(0.0,0.0,0.0,GSY,0.1,9)
      CALL GVALUE(0.0,0.0,0.0,GSY*2,NGVLEN(GSY*2))
C
      DO 200 I=1,NY
         CALL SETLIN(0,2,ICL(I))
         CALL GPLOTP(GX,GY(1,I),1,NX,1,0,0,IPAT(I))
  200 CONTINUE
      CALL SETLIN(0,2,7)
C
 9000 CALL MOVE(GP(1),GP(4)+0.2)
      CALL TEXT(KTITL,6)
C
      RETURN
      END
C
C     ****** DRAW LINES OF 1D GRAPH ******
C
      SUBROUTINE WMGN1D(NX,GX,GXMIN,GXMAX,GY1,GY2,NY,GP,KTITL)
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION GX(NRM),GY1(NRM,NDM),GY2(NRM,NDM),GP(4)
      CHARACTER KTITL*6
C
      IF(NY.EQ.0) RETURN
C
      GYMIN= 1.E32
      GYMAX=-1.E32
      DO I=1,NY
         CALL GMNMX1(GY1(1,I),1,NX,1,GYMIN1,GYMAX1)
         CALL GMNMX1(GY2(1,I),1,NX,1,GYMIN2,GYMAX2)
         GYMIN=MIN(GYMIN,GYMIN1,GYMIN2)
         GYMAX=MAX(GYMAX,GYMAX1,GYMAX2)
      ENDDO
C
C     CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)
      GSX=0.1
C
      IF(GYMIN.GT.0.0.AND.GYMAX.GT.0.0) THEN
         GYMIN=0.0
      ELSE IF(GYMIN.LT.0.0.AND.GYMAX.LT.0.0) THEN
         GYMAX=0.0
      ENDIF
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)
C
      IF(ABS(GSYMAX-GSYMIN).LT.1.E-15) GOTO 9000
C
      CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),GXMIN,GXMAX,GSYMIN,GSYMAX)
      CALL GFRAME
      CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
      CALL GVALUE(GXMIN,GSX*5,0.0,0.0,NGVLEN(GSX*5))
      CALL GSCALE(0.0,0.0,0.0,GSY,0.1,9)
      CALL GVALUE(0.0,0.0,0.0,GSY*2,NGVLEN(GSY*2))
C
      DO 200 I=1,NY
         CALL SETLIN(0,2,7-MOD(I-1,5))
         CALL GPLOTP(GX,GY1(1,I),1,NX,1,0,0,0)
         CALL GPLOTP(GX,GY2(1,I),1,NX,1,0,0,2)
  200 CONTINUE
      CALL SETLIN(0,2,7)
C
 9000 CALL MOVE(GP(1),GP(4)+0.2)
      CALL TEXT(KTITL,6)
C
      RETURN
      END
C
C     ****** DRAW 2D POLOIDAL GRAPH ******
C
      SUBROUTINE WMGRTH(K2,K3,K4)
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION GY(NRM,MDM),GXR(NRM)
C
      CHARACTER K2,K3,K4
C
      IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
         IF(K3.EQ.'R') THEN
            NG3=1
         ELSEIF(K3.EQ.'T') THEN
            NG3=2
         ELSEIF(K3.EQ.'Z') THEN
            NG3=3
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRTH'
            GOTO 9000
         ENDIF
      ELSEIF(K2.EQ.'P') THEN
         IF(K3.EQ.'E') THEN
            NG3=1
         ELSEIF(K3.EQ.'D') THEN
            NG3=2
         ELSEIF(K3.EQ.'T') THEN
            NG3=3
         ELSEIF(K3.EQ.'1') THEN
            NG3=1
         ELSEIF(K3.EQ.'2') THEN
            NG3=2
         ELSEIF(K3.EQ.'3') THEN
            NG3=3
         ELSEIF(K3.EQ.'4') THEN
            NG3=4
         ELSEIF(K3.EQ.'5') THEN
            NG3=5
         ELSEIF(K3.EQ.'6') THEN
            NG3=6
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRTH'
            GOTO 9000
         ENDIF
      ELSE IF(K2.EQ.'J') THEN
         NG3=0
      ELSE
         WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #2 IN WMGRTH'
         GOTO 9000
      ENDIF
C
    1 IF(NPHMAX.EQ.1) THEN
         NPH=1
      ELSE
         WRITE(6,*) '## INPUT NPH : 1..',NPHMAX
         READ(5,*,ERR=1,END=9000) NPH
      END IF
      IF(NPH.LT.1.OR.NPH.GT.NPHMAX) THEN
         WRITE(6,*) 'XX ILLEGAL NPH'
         GOTO 1
      ENDIF
C
      IF(K2.EQ.'E') THEN
         IF(K4.EQ.'R') THEn
            DO 100 NTH=1,MDSIZ
            DO 100 NR=1,NRMAX+1
               GY(NR,NTH)=GCLIP(DBLE(CEFLD(NG3,NTH,NPH,NR)))
  100       CONTINUE
         ELSEIF(K4.EQ.'I') THEN
            DO 110 NTH=1,MDSIZ
            DO 110 NR=1,NRMAX+1
               GY(NR,NTH)=GCLIP(DIMAG(CEFLD(NG3,NTH,NPH,NR)))
  110       CONTINUE
         ELSEIF(K4.EQ.'A') THEN
            DO 120 NTH=1,MDSIZ
            DO 120 NR=1,NRMAX+1
               GY(NR,NTH)=GCLIP(ABS(CEFLD(NG3,NTH,NPH,NR)))
  120       CONTINUE
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRTH'
            GOTO 9000
         ENDIF
         NX=NRMAX+1
      ELSEIF(K2.EQ.'B') THEN
         IF(K4.EQ.'R') THEN
            DO 200 NTH=1,MDSIZ
            DO 200 NR=1,NRMAX+1
               GY(NR,NTH)=GCLIP(DBLE(CBFLD(NG3,NTH,NPH,NR)))
  200       CONTINUE
         ELSEIF(K4.EQ.'I') THEN
            DO 210 NTH=1,MDSIZ
            DO 210 NR=1,NRMAX+1
               GY(NR,NTH)=GCLIP(DIMAG(CBFLD(NG3,NTH,NPH,NR)))
  210       CONTINUE
         ELSEIF(K4.EQ.'A') THEN
            DO 220 NTH=1,MDSIZ
            DO 220 NR=1,NRMAX+1
               GY(NR,NTH)=GCLIP(ABS(CBFLD(NG3,NTH,NPH,NR)))
  220       CONTINUE
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRTH'
            GOTO 9000
         ENDIF
         NX=NRMAX+1
      ELSEIF(K2.EQ.'P') THEN
         DO 300 NTH=1,MDSIZ
         DO 300 NR=1,NRMAX
            GY(NR,NTH)=GCLIP(PABS(NTH,NPH,NR,NG3))
  300    CONTINUE
         NX=NRMAX
      ELSE IF (K2.EQ.'J') THEN
          DO 310 NTH=1,MDSIZ
          DO 310 NR=1,NRMAX
             GY(NR,NTH)=GCLIP(PCUR(NTH,NPH,NR))
  310     CONTINUE
          NX=NRMAX
      END IF
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
C
      DO 400 I=1,NX
        GXR(I)=GCLIP(XR(I))
  400 CONTINUE
C
      CALL GMNMX1(GXR,1,NX,1,GRMIN,GRMAX)
      CALL GMNMX2(GY,NRM,1,NX,1,1,MDSIZ,1,GZ1MIN,GZ1MAX)
      CALL GQSCAL(GZ1MIN,GZ1MAX,GZMIN,GZMAX,GZSTEP)
      NSTEP=10
      GZSTEP=(GZMAX-GZMIN)/NSTEP
C
      CALL GDEFIN(2.5,17.5,2.5,17.5,-GRMAX,GRMAX,-GRMAX,GRMAX)
      CALL SETLIN(0,2,4)
      CALL CIRCLE(GRMAX)
      IF(GZMIN.GT.0.0) THEN
         CALL SETLIN(0,2,6)
         CALL CONTQ3(GY,GXR,NRM,NX,MDSIZ,
     &               GZMIN,GZSTEP,NSTEP,0,KACONT)
      ELSEIF(GZMAX.LT.0.0) THEN
         CALL SETLIN(0,2,5)
         CALL CONTQ3(GY,GXR,NRM,NX,MDSIZ,
     &               GZMIN,GZSTEP,NSTEP,2,KACONT)
      ELSE
         CALL SETLIN(0,2,6)
         CALL CONTQ3(GY,GXR,NRM,NX,MDSIZ,
     &               0.5*GZSTEP,GZSTEP,NSTEP,0,KACONT)
         CALL SETLIN(0,2,5)
         CALL CONTQ3(GY,GXR,NRM,NX,MDSIZ,
     &              -0.5*GZSTEP,-GZSTEP,NSTEP,2,KACONT)
      END IF
      CALL SETLIN(0,0,7)
C
      CALL MOVE(16.5,17.0)
      IF(K2.EQ.'E') CALL TEXT('E ',2)
      IF(K2.EQ.'B') CALL TEXT('B ',2)
      IF(K2.EQ.'P') CALL TEXT('P ',2)
      IF(K2.EQ.'C') CALL TEXT('J ',2)
      IF(K3.EQ.'R') CALL TEXT('r ',2)
      IF(K3.EQ.'T') CALL TEXT('theta ',6)
      IF(K3.EQ.'Z') CALL TEXT('phi ',4)
      IF(K3.EQ.'E') CALL TEXT('e ',2)
      IF(K3.EQ.'D') CALL TEXT('D ',2)
      IF(K3.EQ.'T') CALL TEXT('T ',2)
      IF(K3.EQ.'1') CALL TEXT('e ',2)
      IF(K3.EQ.'2') CALL TEXT('D ',2)
      IF(K3.EQ.'3') CALL TEXT('T ',2)
      IF(K3.EQ.'4') CALL TEXT('a ',2)
      IF(K3.EQ.'5') CALL TEXT('He',2)
      IF(K3.EQ.'6') CALL TEXT('Be',2)
      IF(K4.EQ.'R') CALL TEXT('Real ',5)
      IF(K4.EQ.'I') CALL TEXT('Imag ',5)
      IF(K4.EQ.'A') CALL TEXT('Abs ',4)
      CALL MOVE(16.5,16.0)
      CALL TEXT('MAX :',5)
      CALL NUMBR(GZ1MAX,'(1PE12.4)',12)
      CALL MOVE(16.5,15.5)
      CALL TEXT('MIN :',5)
      CALL NUMBR(GZ1MIN,'(1PE12.4)',12)
      CALL MOVE(16.5,15.0)
      CALL TEXT('STEP:',5)
      CALL NUMBR(GZSTEP,'(1PE12.4)',12)
C
      CALL WMGPRM('C',K3,0,0,0,0)
C
      CALL PAGEE
C
      RETURN
 9000 END
C
C     ****** DRAW 2D MODE GRAPH ******
C
      SUBROUTINE WMGRMD(K2,K3,K4)
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION GY(NRM,MDM+1)
C
      CHARACTER K2,K3,K4
C
      IF(K2.EQ.'E'.OR.K2.EQ.'B') THEN
         IF(K3.EQ.'R') THEN
            NG3=1
         ELSEIF(K3.EQ.'T') THEN
            NG3=2
         ELSEIF(K3.EQ.'Z') THEN
            NG3=3
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRMD'
            GOTO 9000
         ENDIF
      ELSEIF(K2.EQ.'P') THEN
         IF(K3.EQ.'E') THEN
            NG3=1
         ELSEIF(K3.EQ.'H') THEN
            NG3=2
         ELSEIF(K3.EQ.'D') THEN
            NG3=3
         ELSEIF(K3.EQ.'1') THEN
            NG3=1
         ELSEIF(K3.EQ.'2') THEN
            NG3=2
         ELSEIF(K3.EQ.'3') THEN
            NG3=3
         ELSEIF(K3.EQ.'4') THEN
            NG3=4
         ELSEIF(K3.EQ.'5') THEN
            NG3=5
         ELSEIF(K3.EQ.'6') THEN
            NG3=6
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #3 IN WMGRMD'
            GOTO 9000
         ENDIF
      ELSE IF(K2.EQ.'J') THEN
         NG3=0
      ELSE
         WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #2 IN WMGRMD'
         GOTO 9000
      ENDIF
C
    1 IF(NPHMAX.EQ.1) THEN
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
C
      IF (K2.EQ.'E') THEN
         IF (K4.EQ.'R') THEN
            DO 100 NR=1,NRMAX+1
               DO 101 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(DBLE(CEFLDK(NG3,MDX,NDX,NR)))
  101          CONTINUE
               GY(NR,1)=GCLIP(DBLE(CEFLDK(NG3,MDSIZ,NDX,NR)))
  100       CONTINUE
         ELSE IF (K4.EQ.'I') THEN
            DO 110 NR=1,NRMAX+1
               DO 111 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(DIMAG(CEFLDK(NG3,MDX,NDX,NR)))
  111          CONTINUE
               GY(NR,1)=GCLIP(DIMAG(CEFLDK(NG3,MDSIZ,NDX,NR)))
  110       CONTINUE
         ELSE IF (K4.EQ.'A') THEN
            DO 120 NR=1,NRMAX+1
               DO 121 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(ABS(CEFLDK(NG3,MDX,NDX,NR)))
  121          CONTINUE
               GY(NR,1)=GCLIP(ABS(CEFLDK(NG3,MDSIZ,NDX,NR)))
  120       CONTINUE
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRMD'
            GOTO 9000
         ENDIF
         NX=NRMAX+1
      ELSE IF (K2.EQ.'B') THEN
         IF (K4.EQ.'R') THEN
            DO 105 NR=1,NRMAX+1
               DO 106 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(DBLE(CBFLDK(NG3,MDX,NDX,NR)))
  106          CONTINUE
               GY(NR,1)=GCLIP(DBLE(CBFLDK(NG3,MDSIZ,NDX,NR)))
  105       CONTINUE
         ELSE IF (K4.EQ.'I') THEN
            DO 115 NR=1,NRMAX+1
               DO 116 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(DIMAG(CBFLDK(NG3,MDX,NDX,NR)))
  116          CONTINUE
               GY(NR,1)=GCLIP(DIMAG(CBFLDK(NG3,MDSIZ,NDX,NR)))
  115       CONTINUE
         ELSE IF (K4.EQ.'A') THEN
            DO 125 NR=1,NRMAX+1
               DO 126 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(ABS(CBFLDK(NG3,MDX,NDX,NR)))
  126          CONTINUE
               GY(NR,1)=GCLIP(ABS(CBFLDK(NG3,MDSIZ,NDX,NR)))
  125       CONTINUE
         ELSE
            WRITE(6,*) 'XX UNDEFINED CONTROL CHAR #4 IN WMGRMD'
            GOTO 9000
         ENDIF
         NX=NRMAX+1
      ELSE IF (K2.EQ.'P') THEN
            DO 300 NR=1,NRMAX
               DO 290 MDX=1,MDSIZ
                  GY(NR,MDX+1)=GCLIP(PABSK(MDX,NDX,NR,NG3))
  290          CONTINUE
               GY(NR,1)=GCLIP(PABSK(MDSIZ,NDX,NR,NG3))
  300       CONTINUE
            NX=NRMAX
      END IF
C
C
      CALL PAGES
C
      CALL SETCHS(0.3,0.0)
C
      CALL WMGCON(GY,K2,K3,K4,NX)
C
      CALL WMGPRM('M',K3,0,0,0,0)
C
      CALL PAGEE
C
 9000 RETURN
      END
C
C     ****** DRAW GRAPH OF MULTIPLE LINES ******
C
      SUBROUTINE WMGCON(GZL,K2,K3,K4,NX)
C
      INCLUDE 'wmcomm.h'
C
      DIMENSION GY(MDM+1),GZL(NRM,MDM+1)
      CHARACTER K2,K3,K4
C
      GXMIN=GCLIP(XR(1))
      GXMAX=GCLIP(XR(NRMAX+1))
C
      DO 200 MD=MDMIN-1,MDMAX
         MDX=MD-MDMIN+2
         GY(MDX)=NTH0+MD
  200 CONTINUE
C
      GYMIN=GY(1)
      GYMAX=GY(MDSIZ+1)
C
      CALL GMNMX2(GZL,NRM,1,NX,1,1,MDSIZ+1,1,GZMIN,GZMAX) 
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSZ)
C
      NSTEP=10
      GZSTEP=(GZMAX-GZMIN)/NSTEP
C
      CALL GDEFIN(3.0,21.6,1.0,16.5,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL GFRAME
      CALL GSCALE(GXMIN,GSX,0.0,0.0,0.1,9)
      CALL GVALUE(GXMIN,GSX*2,0.0,0.0,1)
      CALL GSCALE(0.0,0.0,0.0,MAX(1.0,GSY),0.1,9)
      CALL GVALUE(0.0,0.0,0.0,MAX(1.0,GSY),0)
C
      IF (GZMIN.GT.0.0) THEN
         CALL SETLIN(0,2,6)
         CALL CONTQ1(GZL,NRM,NX,MDSIZ+1,GZMIN,GZSTEP,NSTEP,0,0,KACONT)
      ELSE IF (GZMAX.LT.0.0) THEN
         CALL SETLIN(0,2,5)
         CALL CONTQ1(GZL,NRM,NX,MDSIZ+1,GZMIN,GZSTEP,NSTEP,0,2,KACONT)
      ELSE 
         CALL SETLIN(0,2,6)
         CALL CONTQ1(GZL,NRM,NX,MDSIZ+1, 0.5*GZSTEP,GZSTEP,NSTEP,
     &             0,0,KACONT)
         CALL SETLIN(0,2,5)
         CALL CONTQ1(GZL,NRM,NX,MDSIZ+1,-0.5*GZSTEP,-GZSTEP,NSTEP,
     &             0,2,KACONT)
      END IF
      CALL SETLIN(0,0,7)
C
      CALL MOVE(1.0,17.0)
      CALL TEXT('mode',4)
      CALL MOVE(4.5,17.3)
      IF(K2.EQ.'E') CALL TEXT('E ',2)
      IF(K2.EQ.'B') CALL TEXT('B ',2)
      IF(K2.EQ.'P') CALL TEXT('P ',2)
      IF(K2.EQ.'C') CALL TEXT('J ',2)
      IF(K3.EQ.'R') CALL TEXT('  r   ',6)
      IF(K3.EQ.'T') CALL TEXT('theta ',6)
      IF(K3.EQ.'Z') CALL TEXT('  z   ',6)
      IF(K3.EQ.'E') CALL TEXT('e ',2)
      IF(K3.EQ.'D') CALL TEXT('D ',2)
      IF(K3.EQ.'T') CALL TEXT('T ',2)
      IF(K4.EQ.'R') CALL TEXT('Real',4)
      IF(K4.EQ.'I') CALL TEXT('Imag',4)
      IF(K4.EQ.'A') CALL TEXT('Abs ',4)
C
      CALL MOVE(10.0,17.3)
      CALL TEXT('min=',4)
      CALL NUMBR(GZMIN,'(1PE12.4)',12)
      CALL MOVE(10.0,16.8)
      CALL TEXT('max=',4)
      CALL NUMBR(GZMAX,'(1PE12.4)',12)
      CALL MOVE(16.0,17.3)
      CALL TEXT('step=',5)
      CALL NUMBR(GZSTEP,'(1PE12.4)',12)
C
      RETURN
      END
C
C
C     ****** WRITE PARAMETERS ******
C
      SUBROUTINE WMGPRM(K1,K3,NTH,NPH,MD,ND)
C
      INCLUDE 'wmcomm.h'
C
      CHARACTER K1,K3
      REAL*4 XPOS,YPOS,DY
C
      CALL SETCHS(0.3,0.0)
      XPOS=22.0
      YPOS=17.5
      DY=0.4
C
      CALL MOVE(XPOS,YPOS)
C
      IF(MODELP.EQ.0) THEN
        CALL TEXT('Wave Guide',10)
      ELSE IF (MODELP.EQ.1) THEN
         CALL TEXT('MHD',3) 
      ELSE IF (MODELP.EQ.2) THEN
         CALL TEXT('COLD',4) 
      ELSE IF (MODELP.EQ.3) THEN
         CALL TEXT('HOT1',4) 
      ELSE IF (MODELP.EQ.4) THEN
         CALL TEXT('HOT2',4) 
      ELSE IF (MODELP.EQ.5) THEN
         CALL TEXT('HOT3',4) 
      END IF
C
C     ****** DEVICE PARAMETERS ******
C
      CALL MOVE(XPOS,YPOS-DY)
      CALL TEXT('BB: ',4)
      CALL NUMBD(BB,'(F8.4)',8)
C
      CALL MOVE(XPOS,YPOS-2*DY)
      CALL TEXT('RR: ',4)
      CALL NUMBD(RR,'(F8.4)',8)
C
      CALL MOVE(XPOS,YPOS-3*DY)
      CALL TEXT('Q0: ',4)
      CALL NUMBD(Q0,'(F8.4)',8)
C
      CALL MOVE(XPOS,YPOS-4*DY)
      CALL TEXT('QA: ',4)
      CALL NUMBD(QA,'(F8.4)',8)
C
C     ****** PLASMA PARAMETERS ******
C
      YPOS=YPOS-5*DY
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PN: ',4)
      DO 100 NS=1,NSMAX
         CALL MOVE(XPOS,YPOS-NS*DY)
         CALL NUMBD(PN(NS),'(F12.4)',12)
  100 CONTINUE
C
      YPOS=YPOS-7*DY
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PTPR: ',6)
      DO 200 NS=1,NSMAX
         CALL MOVE(XPOS,YPOS-NS*DY)
         CALL NUMBD(PTPR(NS),'(F12.4)',12)
  200 CONTINUE
C
C     ****** WAVE PARAMETERS ******
C
      YPOS= YPOS-8*DY
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RF:',3)
      CALL NUMBD(DBLE(CRF),'(F9.4)',9)
C
      CALL MOVE(XPOS,YPOS-DY)
      CALL TEXT('RFI:',4)
      CALL MOVE(XPOS,YPOS-2*DY)
      CALL NUMBD(DIMAG(CRF),'(1PD12.4)',12)
C
      CALL MOVE(XPOS,YPOS-3*DY)
      CALL TEXT('NPH0:',5)
      CALL NUMBI(NPH0,'(I7)',7)
C
      CALL MOVE(XPOS,YPOS-4*DY)
      CALL TEXT('NTH0:',5)
      CALL NUMBI(NTH0,'(I7)',7)
C
      CALL MOVE(XPOS,YPOS-5*DY)
      CALL TEXT('NRMAX:',6)
      CALL NUMBI(NRMAX,'(I6)',6)
C
      CALL MOVE(XPOS,YPOS-6*DY)
      CALL TEXT('NTHMAX:',7)
      CALL NUMBI(NTHMAX,'(I5)',5)
C
      CALL MOVE(XPOS,YPOS-7*DY)
      CALL TEXT('NPHMAX:',7)
      CALL NUMBI(NPHMAX,'(I5)',5)
C
C     ****** COMPUTATION RESULTS ******
C
      YPOS=YPOS-9*DY
      CALL MOVE(XPOS,YPOS)
      IF(K1.EQ.'R'.AND.K3.EQ.'M') THEN
         MDX=MD-MDMIN+1
         NDX=ND-NDMIN+1
         CALL TEXT('Pabs(MD):',9)
         DO 300 NS=1,NSMAX
            CALL MOVE(XPOS,YPOS-NS*DY)
            CALL NUMBD(PABSKT(MDX,NDX,NS),'(1PD12.4)',12)
  300    CONTINUE
      ELSE
         CALL TEXT('Pabs:',5)
         DO 400 NS=1,NSMAX
            CALL MOVE(XPOS,YPOS-NS*DY)
            CALL NUMBD(PABST(NS),'(1PD12.4)',12)
  400    CONTINUE
      ENDIF
C
      YPOS=YPOS-7*DY
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('CUR:',4)
      CALL MOVE(XPOS,YPOS-DY)
      CALL NUMBD(PCURTT,'(1PE12.4)',12)
C
      YPOS=YPOS-2*DY
      CALL MOVE(XPOS,YPOS)
      IF(K1.EQ.'R'.AND.K3.EQ.'M') THEN
         NDX=ND-NDMIN+1
         MDX=MD-MDMIN+1
         CALL TEXT('Imp(MD):',8)
         CALL MOVE(XPOS,YPOS-DY)
         CALL NUMBD(DBLE(CRADKT(MDX,NDX)),'(1PE12.4)',12)
         CALL MOVE(XPOS,YPOS-2*DY)
         CALL NUMBD(DIMAG(CRADKT(MDX,NDX)),'(1PE12.4)',12)
      ELSE
         CALL TEXT('Imp:',4)
         CALL MOVE(XPOS,YPOS-DY)
         CALL NUMBD(DBLE(CRADTT),'(1PE12.4)',12)
         CALL MOVE(XPOS,YPOS-2*DY)
         CALL NUMBD(DIMAG(CRADTT),'(1PE12.4)',12)
      ENDIF
C
C     ****** THETA/MODE ******
C
      YPOS=YPOS-4*DY
      CALL MOVE(XPOS,YPOS)
      IF(K1.EQ.'R') THEN
         IF(K3.EQ.'M') THEN
            CALL TEXT('MD:',3)
            CALL NUMBI(MD,'(I7)',7)
            CALL MOVE(XPOS,YPOS-DY)
            CALL TEXT('ND:',3)
            CALL NUMBI(ND,'(I7)',7)
         ELSE
            CALL TEXT('Th:',3)
            CALL NUMBD(DBLE(NTH-1)/DBLE(NTHMAX)*360.D0,'(F9.4)',9)
            CALL MOVE(XPOS,YPOS-DY)
            CALL TEXT('Ph:',3)
            CALL NUMBD(DBLE(NPH-1)/DBLE(NPHMAX)*360.D0,'(F9.4)',9)
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     ****** AVOID REAL*4 UNDERFLOW ******
C
      REAL FUNCTION GCLIP(D)
C
      REAL*8 D
C
      IF(ABS(D).GT.1.D-30) then
         GCLIP=SNGL(D)
      ELSE
         GCLIP=0.0
      ENDIF
      RETURN
      END
C
C     *****************************
C
C     OPTIMUM NUM LENGTH FOR GVALUE
C
C     *****************************
C
      FUNCTION NGVLEN(GSTEP)
C
      NGX = -INT(LOG10(DBLE(GSTEP*0.11)))
      IF(NGX.LT.-5)THEN
         NGX=-1
      ELSE
         IF(NGX.LT.0) NGX=0
         IF(NGX.GT.5) NGX=-1
      ENDIF
      NGVLEN=NGX
      RETURN
      END
