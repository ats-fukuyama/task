C     $Id$
C
C     ********** F.E.M. DIVIDER ( FIRST ORDER ) **********
C
      SUBROUTINE WFDIV
C
      USE libchar
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
    1 WRITE(6,601) 
  601 FORMAT(1H ,'## INPUT: D/DIV  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',
     &           'W/LIST  X/EXIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL toupper(KID)
C
      IF(KID.EQ.'D') THEN
C
    2    WRITE(6,602) 
  602    FORMAT(1H ,'## TYPE:  X,R,P/POLY  M,V/MIRROR  ',
     &              'C/CIRCLE  T/TOKAMAK  H/HELICAL')
         READ(5,'(A1)',ERR=2,END=1) KID
         CALL toupper(KID)
C
    3    CONTINUE
         IF(KID.EQ.'X') THEN
            MODELV=0
            ID1DIV=0
            CALL DFNODX(IERR)
         ELSEIF(KID.EQ.'R') THEN
            MODELV=1
            ID1DIV=0
            CALL DFNODR(IERR)
         ELSEIF(KID.EQ.'P') THEN
            MODELV=2
            ID1DIV=1
            CALL DFNODP(IERR)
         ELSEIF(KID.EQ.'M') THEN
            MODELV=3
            ID1DIV=1
            CALL DFNODM(IERR)
         ELSEIF(KID.EQ.'V') THEN
            MODELV=4
            ID1DIV=1
            CALL DFNODV(IERR)
         ELSEIF(KID.EQ.'C') THEN
            MODELV=5
            ID1DIV=1
            CALL DFNODC(IERR)
         ELSEIF(KID.EQ.'T') THEN
            MODELV=6
            ID1DIV=1
            CALL DFNODC(IERR)
         ELSEIF(KID.EQ.'H') THEN
            MODELV=7
            ID1DIV=1
            CALL DFNODC(IERR)
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID: ',KID
            IERR=1
         ENDIF
            IF(IERR.NE.0) THEN
               WRITE(6,'(A,I5)') 'XX WFDIV: IERR=',IERR
               GOTO 2
            ENDIF
C
         CALL SETNOD(IERR)
            IF(IERR.NE.0) GOTO 3
         CALL SETELM(IERR)
            IF(IERR.NE.0) GOTO 3
C
         WRITE(6,600) NNOD,NELM
  600    FORMAT(1H ,'NNOD  =',I6,2X,'NELM  =',I6)
C
         CALL WFGINI
         NDMAX=0
         NBWMAX=0
         NBPMAX=0
         NMMAX=1
         NVWMAX=0
         NVPMAX=0
C
      ELSEIF(KID.EQ.'G') THEN
         CALL WFGDIV
      ELSEIF(KID.EQ.'W') THEN
         CALL WFLDIV
      ELSEIF(KID.EQ.'L') THEN
         CALL WFRELM
         CALL WFRZON
      ELSEIF(KID.EQ.'P') THEN
         CALL WFPARM
      ELSEIF(KID.EQ.'V') THEN
         CALL WFVIEW
      ELSEIF(KID.EQ.'S') THEN
         CALL WFWELM
      ELSEIF(KID.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ****** Define 2-D Node Array (RECTANGULAR) ******
C
      SUBROUTINE DFNODX(IERR)
C
      INCLUDE 'wfcomm.inc'
C
    1 WRITE(6,601) BXMIN,BXMAX,BYMIN,BYMAX
  601 FORMAT(1H ,'## DIV:   BXMIN,BXMAX,BYMIN,BYMAX = ',4F10.4/
     &       1H ,'## INPUT: BXMIN,BXMAX,BYMIX,BYMAX ? ')
      READ(5,*,ERR=1,END=9100) BXMIN,BXMAX,BYMIN,BYMAX
C
    7 WRITE(6,607) DELX,DELY
  607 FORMAT(1H ,'## DIV:   DELX,DELY = ',2F10.4/
     &       1H ,'## INPUT: DELX,DELY ? ')
      READ(5,*,ERR=7,END=1) DELX,DELY
      IF(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 1
C
      XWEIGH=0.5D0*(BXMIN+BXMAX)
      YWEIGH=0.5D0*(BYMIN+BYMAX)
C
      NYMAX=NINT((BYMAX-BYMIN)/DELY)+1
      NYMAX=NYMAX+MOD(NYMAX+1,2)
         IF(NYMAX.GT.NYM) GOTO 9200
      DY=(BYMAX-BYMIN)/DBLE(NYMAX-1)
C
      DO 10 NY=1,NYMAX
         Y=BYMIN+DY*(NY-1)
C
         XL(NY)=BXMIN
         XR(NY)=BXMAX
         NXMAX=NINT((XR(NY)-XL(NY))/DELX)+1
         NXMAX=NXMAX+MOD(NXMAX+1,2)
            IF(NXMAX.GT.NXM) GOTO 9300
         DX=(XR(NY)-XL(NY))/DBLE(NXMAX-1)
C
         NXA(NY)=NXMAX
      DO 10 NX=1,NXMAX
         X=DX*(NX-1)+XL(NY)
         XDA(NX,NY)=X
         YDA(NX,NY)=Y
         XDW(NX,NY)=XWEIGH
         YDW(NX,NY)=YWEIGH
   10 CONTINUE
C
      IERR=0
      RETURN
C
 9100 IERR=1
      RETURN
C
 9200 WRITE(6,602) NYMAX,NYM
  602 FORMAT(1H ,'DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9300 WRITE(6,603) NY,NXMAX,NXM
  603 FORMAT(1H ,'DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Define 2-D Node Array (RECTANGULARS) ******
C
      SUBROUTINE DFNODR(IERR)
C
      INCLUDE 'wfcomm.inc'
      PARAMETER (NXQM=11,NYQM=11)
      DIMENSION XRECT(NXQM,NYQM),YRECT(NYQM)
      DIMENSION DELXR(NXQM,NYQM),DELYR(NYQM)
      DIMENSION NXQMAX(NYQM)
C     
    1 WRITE(6,601)
  601 FORMAT(1H ,'## DIV:   NYQMAX : number of rectangulars in y + 1')
      READ(5,*,ERR=1,END=9100) NYQMAX
      IF(NYQMAX.LE.1) THEN
         WRITE(6,*) 'XX NYQMAX should be greater than 1'
         GOTO 1
      ENDIF
C     
    2 WRITE(6,602)
  602 FORMAT(1H ,'## DIV:   YRECT(1) : lower boundary')
      READ(5,*,ERR=2,END=1) YRECT(1)
C
      DO NYQ=1,NYQMAX-1
C
    3    WRITE(6,603) NYQ
  603    FORMAT(1H ,'## DIV: NYQ= ',I2,
     &              ' : NXQMAX : number of rectangulars in x + 1 ')
         READ(5,*,ERR=3,END=2) NXQMAX(NYQ)
         IF(NXQMAX(NYQ).LE.1) THEN
            WRITE(6,*) 'XX NXQMAX should be greater than 1'
            GOTO 1
         ENDIF
C
         DO NXQ=1,NXQMAX(NYQ)
    4       WRITE(6,604) NXQ,NYQ
  604       FORMAT(1H ,'## DIV:   XRECT(',I2,',',I2,')')
            READ(5,*,ERR=4,END=3) XRECT(NXQ,NYQ)
         ENDDO
C
    5    WRITE(6,605) NYQ+1
  605    FORMAT(1H ,'## DIV:   YRECT(',I2,')')
         READ(5,*,ERR=5,END=3) YRECT(NYQ+1)
      ENDDO
      NXQMAX(NYQMAX)=NXQMAX(NYQMAX-1)
      DO NXQ=1,NXQMAX(NYQMAX)
         XRECT(NXQ,NYQMAX)=XRECT(NXQ,NYQMAX-1)
      ENDDO
C
      BYMIN=YRECT(1)
      BYMAX=YRECT(NYQMAX)
      BXMIN=XRECT(1,1)
      BXMAX=XRECT(NXQMAX(1),1)
      DO NYQ=2,NYQMAX
         BXMIN=MIN(BXMIN,XRECT(1,NYQ))
         BXMAX=MAX(BXMAX,XRECT(NXQMAX(NYQ),NYQ))
      ENDDO
C
      NYMAX=1
      DO NYQ=1,NYQMAX-1
    6    WRITE(6,606) NYQ,YRECT(NYQ),YRECT(NYQ+1),DELY
  606    FORMAT(1H ,'## DIV:   NYQ,YLOWER,YUPPER,DELY = ',I4,3F10.4/
     &          1H ,'## INPUT: DELY ? ')
         READ(5,*,ERR=6,END=1) DELY
         IF(ABS(DELY).LE.1.D-6) GOTO 6
         DYR=YRECT(NYQ+1)-YRECT(NYQ)
         NYMAXL=MAX(NINT(DYR/DELY),1)
         IF(NYMAX+NYMAXL.GT.NYM) GOTO 9200
         DY=DYR/DBLE(NYMAXL)
         DELYR(NYQ)=DY
C
         NXMAX=1
         DO NXQ=1,NXQMAX(NYQ)-1
    7       WRITE(6,607) NYQ,NXQ,XRECT(NXQ,NYQ),XRECT(NXQ+1,NYQ),DELX
  607       FORMAT(1H ,'## DIV:   NYQ,NXQ,XLEFT,YRIGHT,DELX = ',
     &                 2I4,3F10.4/
     &             1H ,'## INPUT: DELX ? ')
            READ(5,*,ERR=7,END=6) DELX
            IF(ABS(DELX).LE.1.D-6) GOTO 7
            DXR=XRECT(NXQ+1,NYQ)-XRECT(NXQ,NYQ)
            NXMAXL=MAX(NINT(DXR/DELX),1)
            IF(NXMAX+NXMAXL.GT.NXM) GOTO 9300
            DX=DXR/DBLE(NXMAXL)
            DELXR(NXQ,NYQ)=DX
         ENDDO
      ENDDO
C
      NYMAX=0
      MODE=0
      DO NYQ=1,NYQMAX-1
         NYMIN=NYMAX+1
         DYR=YRECT(NYQ+1)-YRECT(NYQ)
         DY=DELYR(NYQ)
         IF(MODE.EQ.0) THEN
            YMIN=YRECT(NYQ)
         ELSE
            YMIN=YRECT(NYQ)+DY
         ENDIF
         NYMAXL=MAX(NINT(DYR/DY),1)
         YWEIGH1=YRECT(NYQ)+DY*DBLE(NYMAXL/2)
         IF(NXQMAX(NYQ).GE.NXQMAX(NYQ+1)) THEN
            IF(MODE.EQ.0) THEN
               NYMAX=NYMIN+NYMAXL
            ELSE
               NYMAX=NYMIN+NYMAXL-1
            ENDIF
            MODE=1
            IF(NYQ.EQ.NYQMAX-1) THEN
               YWEIGH2=YWEIGH1
            ELSE
               DYR2=YRECT(NYQ+2)-YRECT(NYQ+1)
               DY2=DELYR(NYQ+1)
               NYMAXL2=MAX(NINT(DYR2/DY2),1)
               YWEIGH2=YRECT(NYQ+1)+DY*DBLE(NYMAXL2/2)
            ENDIF
         ELSE
            IF(MODE.EQ.0) THEN
               NYMAX=NYMIN+NYMAXL-1
            ELSE
               NYMAX=NYMIN+NYMAXL-2
            ENDIF
            MODE=0
         ENDIF
C
         NXMAX=0
         DO NXQ=1,NXQMAX(NYQ)-1
            NXMIN=NXMAX+1
            DXR=XRECT(NXQ+1,NYQ)-XRECT(NXQ,NYQ)
            DX=DELXR(NXQ,NYQ)
            NXMAXL=MAX(NINT(DXR/DX),1)
            IF(NXQ.EQ.NXQMAX(NYQ)-1) THEN
               NXMAX=NXMIN+NXMAXL
            ELSE
               NXMAX=NXMIN+NXMAXL-1
            ENDIF
            XWEIGH=XRECT(NXQ,NYQ)+DX*DBLE(NXMAXL/2)
C
            DO NX=NXMIN,NXMAX
               X=XRECT(NXQ,NYQ)+DX*(NX-NXMIN)
               DO NY=NYMIN,NYMAX
                  Y=YMIN+DY*(NY-NYMIN)
                  XDA(NX,NY)=X
                  YDA(NX,NY)=Y
                  IF(NY.EQ.NYMIN+NYMAXL) THEN
                     YWEIGH=YWEIGH2
                  ELSE
                     YWEIGH=YWEIGH1
                  ENDIF
                  XDW(NX,NY)=XWEIGH
                  YDW(NX,NY)=YWEIGH
C                  WRITE(6,'(4I4,2F8.4)') NYQ,NXQ,NX,NY,X,Y
               ENDDO
            ENDDO
         ENDDO
         DO NY=NYMIN,NYMAX
            NXA(NY)=NXMAX
         ENDDO
      ENDDO
C
      IERR=0
      RETURN
C
 9100 IERR=1
      RETURN
C
 9200 WRITE(6,692) NYMAX,NYM
  692 FORMAT(1H ,'DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9300 WRITE(6,693) NY,NXMAX,NXM
  693 FORMAT(1H ,'DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Define 2-D Node Array (OBLIQUE) ******
C
      SUBROUTINE DFNODP(IERR)
C
      INCLUDE 'wfcomm.inc'
C      DATA EPSL/1.D-12/
C
      IERR=0
    1 WRITE(6,610) BXMIN,BXMAX,BYMIN,BYMAX
  610 FORMAT(1H ,'## DIV:   BXMIN,BXMAX,BYMIN,BYMAX = ',4F10.4/
     &       1H ,'## INPUT: BXMIN,BXMAX,BYMIX,BYMAX ? ')
      READ(5,*,ERR=1,END=9100) BXMIN,BXMAX,BYMIN,BYMAX
C
    7 WRITE(6,607) DELX,DELY
  607 FORMAT(1H ,'## DIV:   DELX,DELY = ',2F10.4/
     &       1H ,'## INPUT: DELX,DELY ? ')
      READ(5,*,ERR=7,END=1) DELX,DELY
      IF(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 1
C
      XWEIGH=0.5D0*       BXMAX
      YWEIGH=0.5D0*(BYMIN+BYMAX)
C
      NYMAX=NINT((BYMAX-BYMIN)/DELY)+1
      DY=(BYMAX-BYMIN)/DBLE(NYMAX-1)
C
      DO 30 NY=1,NYMAX
         Y=BYMIN+DY*(NY-1)
         XMAX=BXMIN+(BXMAX-BXMIN)*(NY-1)/DBLE(NYMAX-1)
         NXMAX=NINT(XMAX/DELX)+1
         IF(NXMAX.EQ.1) NXMAX=2
         DX=XMAX/DBLE(NXMAX-1)
         DO NX=1,NXMAX
            XDA(NX,NY)=DX*(NX-1)
            YDA(NX,NY)=Y
            XDW(NX,NY)=XWEIGH
            YDW(NX,NY)=YWEIGH
         ENDDO
         IF(NY.EQ.NYMAX) THEN
            NXA(NY)=NXMAX
         ELSE
            XMAX1=BXMIN+(BXMAX-BXMIN)*NY/DBLE(NYMAX-1)
            NXMAX1=NINT((XMAX1-XMAX)/DELX)
            DX=(XMAX1-XMAX)/DBLE(NXMAX1)
            NXA(NY)=NXMAX+NXMAX1-1
            DO NX=NXMAX+1,NXMAX+NXMAX1-1
               XDA(NX,NY)=XMAX+DX*(NX-NXMAX)
               YDA(NX,NY)=Y+DY*(NX-NXMAX)/DBLE(NXMAX1)
               XDW(NX,NY)=XWEIGH
               YDW(NX,NY)=YWEIGH
            ENDDO
         ENDIF
         
   30 CONTINUE
      RETURN
C
 9100 IERR=1
      RETURN
      END
C
C     ****** Define 2-D Node Array (MIRROR) ******
C
      SUBROUTINE DFNODM(IERR)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION PSIH(NXM)
      EXTERNAL BOUNDM
      DATA EPS/1.D-8/
C
    1 WRITE(6,608) BXMIN,BXMAX,BYMIN,BYMAX
  608 FORMAT(1H ,'## DIV:   BXMIN,BXMAX,BYMIN,BYMAX = ',4F10.4/
     &       1H ,'## INPUT: BXMIN,BXMAX,BYMIX,BYMAX ? ')
      READ(5,*,ERR=1,END=9100) BXMIN,BXMAX,BYMIN,BYMAX
C
    7 WRITE(6,607) DELX,DELY
  607 FORMAT(1H ,'## DIV:   DELX,DELY = ',2F10.4/
     &       1H ,'## INPUT: DELX,DELY ? ')
      READ(5,*,ERR=7,END=1) DELX,DELY
      IF(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 1
C
      IERR=0
C
      XWEIGH=0.5D0*(BXMIN+BXMAX)
      YWEIGH=0.5D0*(BYMIN+BYMAX)
C
      NYMAX=NINT((BYMAX-BYMIN)/DELY)+1
      NYMAX=NYMAX+MOD(NYMAX+1,2)
         IF(NYMAX.GT.NYM) GOTO 9200
      DY=(BYMAX-BYMIN)/DBLE(NYMAX-1)
C
      NXMAX=NINT((BXMAX-BXMIN)/DELX)+1
      NXMAX=NXMAX+MOD(NXMAX+1,2)
         IF(NXMAX.GT.NXM) GOTO 9300
      DX=(BXMAX-BXMIN)/DBLE(NXMAX-1)
C
      NY=1
      Y=0.D0
      PSIH(1)=0.D0
      DO 10 NX=2,NXMAX
         X=DX*(NX-1)
         CALL WFSPSI(X,Y,PSI)
         PSIH(NX)=PSI
   10 CONTINUE
C
      DO 30 NY=1,NYMAX
         X=0.D0
         Y=BYMIN+DY*(NY-1)
         NXA(NY)=NXMAX
         XDA(1,NY)=X
         YDA(1,NY)=Y
         XDW(1,NY)=XWEIGH
         YDW(1,NY)=YWEIGH
         DO 20 NX=2,NXMAX
            YBP=Y
            PSIBP=PSIH(NX)
            XS=XDA(NX-1,NY)
            XE=BXMAX
            DX=DELX*0.3D0
            ILL=0
            CALL FRGFLS(XS,XE,DX,X,BOUNDM,EPS,ILL)
C
            WRITE(6,'(A,2I5,1P3E12.4)') 'NY,NX,XS,XE,X =',NY,NX,XS,XE,X
            WRITE(6,*) 'ILL =',ILL
C
            IF(ILL.GT.0)   GOTO 9000
            XDA(NX,NY)=X
            YDA(NX,NY)=Y
            XDW(NX,NY)=XWEIGH
            YDW(NX,NY)=YWEIGH
   20    CONTINUE
   30 CONTINUE
      RETURN
C
 9000 IERR=1
      RETURN
C
 9100 IERR=1
      RETURN
C
 9200 WRITE(6,602) NYMAX,NYM
  602 FORMAT(1H ,'DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=2
      RETURN
C
 9300 WRITE(6,603) NY,NXMAX,NXM
  603 FORMAT(1H ,'DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=3
      END
C
C
C     ****** Definititon of Boundary for Fixed Y ******
C
      FUNCTION BOUNDM(X)
C
      INCLUDE 'wfcomm.inc'
C
      X1=X
      YBP1=YBP
      CALL WFSPSI(X1,YBP1,PSI)
      BOUNDM=PSI-PSIBP
      RETURN
      END
C
C     ****** Define 2-D Node Array (MIRROR) ******
C
      SUBROUTINE DFNODV(IERR)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION PSIH(NXM)
      EXTERNAL BOUNDV
      DATA EPS/1.D-8/
C
      IERR=0
    1 WRITE(6,610) BXMIN,BXMAX,BYMIN,BYMAX
  610 FORMAT(1H ,'## DIV:   BXMIN,BXMAX,BYMIN,BYMAX = ',4F10.4/
     &       1H ,'## INPUT: BXMIN,BXMAX,BYMIX,BYMAX ? ')
      READ(5,*,ERR=1,END=9100) BXMIN,BXMAX,BYMIN,BYMAX
C
    7 WRITE(6,607) DELX,DELY
  607 FORMAT(1H ,'## DIV:   DELX,DELY = ',2F10.4/
     &       1H ,'## INPUT: DELX,DELY ? ')
      READ(5,*,ERR=7,END=1) DELX,DELY
      IF(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 1
C
      XWEIGH=0.5D0*(BXMIN+BXMAX)
      YWEIGH=0.5D0*(BYMIN+BYMAX)
C
      NYMAX=NINT((BYMAX-BYMIN)/DELY)+1
      NYMAX=NYMAX+MOD(NYMAX+1,2)
         IF(NYMAX.GT.NYM) GOTO 9200
      DY=(BYMAX-BYMIN)/DBLE(NYMAX-1)
C
      NXMAX=NINT((BXMAX-BXMIN)/DELX)+1
      NXMAX=NXMAX+MOD(NXMAX+1,2)
         IF(NXMAX.GT.NXM) GOTO 9300
      DX=(BXMAX-BXMIN)/DBLE(NXMAX-1)
C
      Y=0.D0
      DO 10 NX=1,NXMAX
         X=BXMIN+DX*(NX-1)
         CALL WFSPSI(X,Y,PSI)
         IF(X.LT.0.D0) THEN
            PSI=-SQRT(PSI)
         ELSEIF(X.GT.0.D0) THEN
            PSI= SQRT(PSI)
         ELSE
            PSI=0.D0
         ENDIF
         PSIH(NX)=PSI
C         WRITE(6,*) NX,X,PSI
   10 CONTINUE
C
      DO 30 NY=1,NYMAX
         Y=BYMIN+DY*(NY-1)
         NXA(NY)=NXMAX
         DO 20 NX=1,NXMAX
            YBP=Y
            PSIBP=PSIH(NX)
            IF(NX.EQ.1) THEN
               XS=BXMIN
            ELSE
               XS=XDA(NX-1,NY)
            ENDIF
            XE=BXMAX
            DX=DELX*0.3D0
            ILL=0
            CALL FRGFLS(XS,XE,DX,X,BOUNDV,EPS,ILL)
C
            IF(ILL.GT.0) THEN
            WRITE(6,*) 'NX,NY,YBP,PSIBP=',NX,NY,YBP,PSIBP
            WRITE(6,*) 'XS,XE,DX=',XS,XE,DX
            WRITE(6,*) 'XDA(NX,NY) =',XDA(NX,NY)
            WRITE(6,*) 'ILL =',ILL
            ENDIF
C
            IF(ILL.GT.0)   GOTO 9000
            XDA(NX,NY)=X
            YDA(NX,NY)=Y
            XDW(NX,NY)=XWEIGH
            YDW(NX,NY)=YWEIGH
   20    CONTINUE
   30 CONTINUE
      RETURN
C
 9000 IERR=1
      RETURN
C
 9100 IERR=1
      RETURN
C
 9200 WRITE(6,602) NYMAX,NYM
  602 FORMAT(1H ,'DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9300 WRITE(6,603) NY,NXMAX,NXM
  603 FORMAT(1H ,'DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=1
      END
C
C
C     ****** Definititon of Boundary for Fixed Y ******
C
      FUNCTION BOUNDV(X)
C
      INCLUDE 'wfcomm.inc'
C
      X1=X
      YBP1=YBP
      CALL WFSPSI(X1,YBP1,PSI)
      IF(X.LT.0.D0) THEN
         PSI=-SQRT(PSI)
      ELSEIF(X.GT.0.D0) THEN
         PSI= SQRT(PSI)
      ELSE
         PSI=0.D0
      ENDIF
      BOUNDV=PSI-PSIBP
      RETURN
      END
C
C
C     ****** Define 2-D Node Array (CIRCULAR) ******
C
      SUBROUTINE DFNODC(IERR)
C
      INCLUDE 'wfcomm.inc'
      PARAMETER (NYMH=NYM/2)
      DIMENSION XL1(NYMH),XR1(NYMH)
      DIMENSION NXA1(NYMH),NXR1(NYMH),NXL1(NYMH)
      DIMENSION XDA1(NXM,NYMH),YDA1(NXM,NYMH)
      DIMENSION XU(NXM),YU(NXM),DYU(NXM),NYU(NXM)
      EXTERNAL BOUNDX,BOUNDY
      DATA EPS/1.D-8/
C
    1 CONTINUE
      IF(MODELV.EQ.5) THEN
         WRITE(6,601) RB
  601    FORMAT(1H ,'## DIV:   RB = ',F10.4/
     &          1H ,'## INPUT: RB ? ')
         READ(5,*,ERR=1,END=9100) RB
         BXMIN=-RB
         BXMAX= RB
         BYMIN=-RB
         BYMAX= RB
      ELSEIF(MODELV.EQ.6) THEN
         WRITE(6,605) RB,BKAP,BDEL
  605    FORMAT(1H ,'## DIV:   RB,BKAP,BDEL = ',3F10.4/
     &          1H ,'## INPUT: RB,BKAP,BDEL ? ')
         READ(5,*,ERR=1,END=9100) RB,BKAP,BDEL
         BXMIN=-RB
         BXMAX= RB
         BYMIN=-RB*RKAP
         BYMAX= RB*RKAP
      ELSEIF(MODELV.EQ.7) THEN
         WRITE(6,606) RB,BXMIN,BXMAX
  606    FORMAT(1H ,'## DIV:   RB,BXMIN,BXMAX = ',3F10.4/
     &          1H ,'## INPUT: RB,BXMIN,BXMAX ? ')
         READ(5,*,ERR=1,END=9100) RB,BXMIN,BXMAX
         BYMIN=-RB
         BYMAX= RB
      ENDIF
C
    7 WRITE(6,607) DELX,DELY
  607 FORMAT(1H ,'## DIV:   DELX,DELY = ',2F10.4/
     &       1H ,'## INPUT: DELX,DELY ? ')
      READ(5,*,ERR=7,END=1) DELX,DELY
      IF(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 1
C
      XWEIGH=0.5D0*(BXMIN+BXMAX)
      YWEIGH=0.5D0*(BYMIN+BYMAX)
C
      NYMAX=NINT(BYMAX/DELY)+1
         IF(NYMAX.GT.NYM) GOTO 9200
      IND=0
      XL1(1)=BXMIN
      XR1(1)=BXMAX
C
      DO 10 NYDO=1,NYMAX
         NY=NYDO
         YF=DELY*DBLE(NY-1)
         IF(NY.EQ.1) THEN
            XS=0.5D0*(BXMIN+BXMAX)
         ELSE
            XS=0.5D0*(XL1(NY-1)+XR1(NY-1))
         ENDIF
         XE=BXMIN-1.D-5
         DX=-DELX*0.3D0
         ILL=0
         CALL FRGFLS(XS,XE,DX,XL1(NY),BOUNDX,EPS,ILL)
         IF(ILL.EQ.902) THEN
            NY=NY-1
            GOTO 20
         ENDIF
         IF(ILL.GT.0)   GOTO 9000
         XE=BXMAX+1.D-5
         DX=DELX*0.3D0
         ILL=0
         CALL FRGFLS(XS,XE,DX,XR1(NY),BOUNDX,EPS,ILL)
         IF(ILL.EQ.902) THEN
            NY=NY-1
            GOTO 20
         ENDIF
         IF(ILL.GT.0)   GOTO 9000
         IF(NY.GT.1) THEN
            DXL=ABS(XL1(NY)-XL1(NY-1))/DELX
            DXR=ABS(XR1(NY)-XR1(NY-1))/DELX
            IF(DXL**2+DXR**2.LT.1.D0) THEN
               NYMID=NY
            ELSE
               IND=1
            ENDIF
         ENDIF
   10 CONTINUE
   20 CONTINUE
C
      NYMAX=MIN(NYMAX,NY)
C
      DO 1000 NY=1,NYMID
         NXA1(NY)=NINT((XR1(NY)-XL1(NY))/DELX)+1
            IF(NXA1(NY).GT.NXM) GOTO 9300
         IF (NXA1(NY).EQ.1) THEN
            NXL1(NY)=1
            NXR1(NY)=0
            XDA1(1,NY)=XL1(NY)
            YDA1(1,NY)=DELY*DBLE(NY-1)
            IND=0
            GOTO 1010
         ENDIF
         IF(MOD(NXA1(NY),2).EQ.0) NXA1(NY)=NXA1(NY)+1
         DX=(XR1(NY)-XL1(NY))/DBLE(NXA1(NY)-1)
         NXL1(NY)=1
         NXR1(NY)=1
         DO 1000 NX=1,NXA1(NY)
            XDA1(NX,NY)=DX*DBLE(NX-1)+XL1(NY)
            YDA1(NX,NY)=DELY*DBLE(NY-1)
 1000    CONTINUE
 1010 CONTINUE
C
      IF(IND.GE.1.AND.NYMID.LT.NYMAX) THEN
         NLMID=NXL1(NYMID)
         NRMID=NXR1(NYMID)
         NXMID=NXA1(NYMID)
         YMID =DELY*DBLE(NYMID-1)
         DO 4100 NX=NLMID+1,NXMID-NRMID
            XU(NX)=XDA1(NX,NYMID)
            XF=XU(NX)
            YS=YMID
            YE=BYMAX+1.D-5
            DY=DELY
            ILL=0
            CALL FRGFLS(YS,YE,DY,YU(NX),BOUNDY,EPS,ILL)
            IF(ILL.GT.0) GOTO 9400
            NYU(NX)=MAX(1,NINT((YU(NX)-YMID)/DELY))
            DYU(NX)=(YU(NX)-YMID)/DBLE(NYU(NX))
            NYMAX=MAX(NYMAX,NYU(NX)+NYMID)
 4100    CONTINUE
C
         DO 4300 NY=NYMID+1,NYMAX
            NX=0
            INDC=0
            NXL1(NY)=0
            NXR1(NY)=0
            NYUL=NY-NYMID
            DO 4200 I=NLMID+1,NXMID-NRMID
               IF(NYU(I).EQ.NYUL) THEN
                  NX=NX+1
                  XDA1(NX,NY)=XDA1(I,NYMID)
                  YDA1(NX,NY)=YU(I)
                  IF(INDC.EQ.0) THEN
                     NXL1(NY)=NXL1(NY)+1
                  ELSE
                     NXR1(NY)=NXR1(NY)+1
                  ENDIF
               ELSE IF(NYU(I).GT.NYUL) THEN
                  NX=NX+1
                  XDA1(NX,NY)=XDA1(I,NYMID)
                  YDA1(NX,NY)=DYU(I)*DBLE(NYUL)+YMID
                  INDC=1
               ENDIF
 4200       CONTINUE
            IF(NX.EQ.0) GOTO 4305
            NXA1(NY)=NX
 4300    CONTINUE
         GOTO 4310
 4305    NYMAX=NY-1
 4310    CONTINUE
      ELSE
         NXL1(NYMAX)=NXA1(NYMAX)
         NXR1(NYMAX)=0
      ENDIF
C
      DO 5000 NY=1,NYMAX
         NYA=NYMAX-1+NY
         NYB=NYMAX+1-NY
         XL(NYA)=XL1(NY)
         XL(NYB)=XL1(NY)
         XR(NYA)=XR1(NY)
         XR(NYB)=XR1(NY)
         NXA(NYA)=NXA1(NY)
         NXA(NYB)=NXA1(NY)
      DO 5000 NX=1,NXA1(NY)
         XDA(NX,NYA)= XDA1(NX,NY)
         XDA(NX,NYB)= XDA1(NX,NY)
         YDA(NX,NYA)= YDA1(NX,NY)
         YDA(NX,NYB)=-YDA1(NX,NY)
         XDW(NX,NYA)= XWEIGH
         XDW(NX,NYB)= XWEIGH
         YDW(NX,NYA)= YWEIGH
         YDW(NX,NYB)=-YWEIGH
 5000 CONTINUE
      NYMAX=2*NYMAX-1
C
      IERR=0
      RETURN
C
 9000 WRITE(6,611) NX,NY,YF,XS,XE,XL(NY),XR(NY),ILL
  611 FORMAT(1H ,'DFNODE : FRGFLS ERROR : NX,NY,YF,XS,XE,XL,XR,ILL'/
     &       1H ,2I5,1P5E15.5,I5)
C
 9100 IERR=1
      RETURN
C
 9200 WRITE(6,612) NYMAX,NYM
  612 FORMAT(1H ,'DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9300 WRITE(6,613) NY,NXA(NY),NXM
  613 FORMAT(1H ,'DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=1
      RETURN
C
 9400 WRITE(6,614) NX,NY,XF,YS,YE,YU(NX),ILL
  614 FORMAT(1H ,'NODE1 : FRGFLS ERROR : NX,NY,XF,YS,YE,YU,ILL'/
     &       1H ,2I5,1P4E15.5,I5)
      IERR=1
      RETURN
      END
C
C     ****** Definititon of Boundary for Fixed Y ******
C
      FUNCTION BOUNDX(X)
C
      INCLUDE 'wfcomm.inc'
C
      BOUNDX=BOUNDF(X,YF)
      RETURN
      END
C
C     ****** Definititon of Boundary for Fixed X ******
C
      FUNCTION BOUNDY(Y)
C
      INCLUDE 'wfcomm.inc'
C
      BOUNDY=BOUNDF(XF,Y)
      RETURN
      END
C
C     ****** Definititon of Boundary ******
C
      FUNCTION BOUNDF(X,Y)
C
      INCLUDE 'wfcomm.inc'
C
      IF(MODELV.EQ.0) THEN
         BOUNDF=(X-BXMIN)*(BXMAX-X)*(Y-BYMIN)*(BYMAX-Y)
      ELSEIF(MODELV.EQ.3) THEN
         BOUNDF=(X-BXMIN)*(BXMAX-X)*(Y-BYMIN)*(BYMAX-Y)
      ELSEIF(MODELV.EQ.5) THEN
         BOUNDF=(RB*RB-X*X-Y*Y)
      ELSEIF(MODELV.EQ.6) THEN
         BOUNDF=(RB*RB-X*X-Y*Y/(BKAP*BKAP))
      ELSEIF(MODELV.EQ.7) THEN
         BOUNDF=(RB*RB-X*X-Y*Y)*(X-BXMIN)*(BXMAX-X)
      ENDIF
      RETURN
      END
C
C     ****** Set Node Array ******
C
      SUBROUTINE SETNOD(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      IN=0
      INMAX=NNODM
      DO 15 NY=1,NYMAX
         DO 10 NX=1,NXA(NY)
            IF(IN.GE.INMAX) GOTO 9000
            IN=IN+1
            NDA(NX,NY)=IN
            XD(IN)= XDA(NX,NY)
            YD(IN)= YDA(NX,NY)
C            WRITE(6,'(A,I5,1P2E12.4)') 'IN,XD,YD=',IN,XD(IN),YD(IN)
   10    CONTINUE
   15 CONTINUE
C
      NNOD=IN
      IERR=0
      RETURN
C
 9000 WRITE(6,601) IN,INMAX
  601 FORMAT(1H ,'SETNOD : NNOD EXCEEDS NNODM : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Set Element Array ******
C
      SUBROUTINE SETELM(IERR)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION ICK(NNODM)
C
      IE=0
      SEPS=1.D-12
C
      DO 2000 NYDO=1,NYMAX-1
         NY=NYDO
         NX1MAX=NXA(NY)
         NX2MAX=NXA(NY+1)
         NX1=1
         NX2=1
C
         XAI=XDA(NX1  ,NY  )
         XAJ=XDA(NX2  ,NY+1)
         IF(XAI.LT.XAJ) THEN
   10       CONTINUE
            IF(ID1DIV.EQ.1) THEN
               IF(IE.GE.NELMM) GOTO 9001
               IE=IE+1
               IELM(1,IE)=NDA(NX1  ,NY  )
               IELM(2,IE)=NDA(NX1+1,NY  )
               IELM(3,IE)=NDA(NX2  ,NY+1)
               CALL WFSELM(IE,S)
               IF(S.LE.SEPS) IE=IE-1
            ENDIF
            NX1=NX1+1
            IF(NX1.GT.NX1MAX) GOTO 9003
            XAI=XDA(NX1  ,NY  )
            XAJ=XDA(NX2  ,NY+1)
            IF(XAI.LT.XAJ) GOTO 10
         ELSE IF(XAI.GT.XAJ) THEN
   11       CONTINUE
            IF(ID1DIV.EQ.1) THEN
               IF(IE.GE.NELMM) GOTO 9001
               IE=IE+1
               IELM(1,IE)=NDA(NX1  ,NY  )
               IELM(2,IE)=NDA(NX2+1,NY+1)
               IELM(3,IE)=NDA(NX2  ,NY+1)
               CALL WFSELM(IE,S)
               IF(S.LE.SEPS) IE=IE-1
            ENDIF
            NX2=NX2+1
            IF(NX2.GT.NX2MAX) GOTO 9003
            XAI=XDA(NX1  ,NY  )
            XAJ=XDA(NX2  ,NY+1)
            IF(XAI.GT.XAJ) GOTO 11
         ENDIF
C
   20    CONTINUE
         XWEIGH=XDW(NX1,NY)
         YWEIGH=YDW(NX1,NY)
C
         X=XDA(NX1,  NY   )
         Y=YDA(NX1,  NY   )
         R=SQRT((X-XWEIGH)**2+(Y-YWEIGH)**2)
         FACTOR=1.D0-R*1.D-12
         XAI=XWEIGH+(X-XWEIGH)*FACTOR
         YAI=YWEIGH+(Y-YWEIGH)*FACTOR
         X=XDA(NX1+1,NY   )
         Y=YDA(NX1+1,NY   )
         R=SQRT((X-XWEIGH)**2+(Y-YWEIGH)**2)
         FACTOR=1.D0-R*1.D-12
         XBI=XWEIGH+(X-XWEIGH)*FACTOR
         YBI=YWEIGH+(Y-YWEIGH)*FACTOR
         X=XDA(NX2,  NY+1)
         Y=YDA(NX2,  NY+1)
         R=SQRT((X-XWEIGH)**2+(Y-YWEIGH)**2)
         FACTOR=1.D0-R*1.D-12
         XAJ=XWEIGH+(X-XWEIGH)*FACTOR
         YAJ=YWEIGH+(Y-YWEIGH)*FACTOR
         X=XDA(NX2+1,NY+1)
         Y=YDA(NX2+1,NY+1)
         R=SQRT((X-XWEIGH)**2+(Y-YWEIGH)**2)
         FACTOR=1.D0-R*1.D-12
         XBJ=XWEIGH+(X-XWEIGH)*FACTOR
         YBJ=YWEIGH+(Y-YWEIGH)*FACTOR
C
         VA=(XAI-XBJ)**2+(YAI-YBJ)**2
         VB=(XBI-XAJ)**2+(YBI-YAJ)**2
         IF(IE.GE.NELMM) GOTO 9001
         IE=IE+1
         IF(VA.GT.VB) THEN
            IELM(1,IE)=NDA(NX1  ,NY  )
            IELM(2,IE)=NDA(NX1+1,NY   )
            IELM(3,IE)=NDA(NX2  ,NY+1)
            CALL WFSELM(IE,S)
            IF(S.LE.SEPS) IE=IE-1
            NX1=NX1+1
         ELSE
            IELM(1,IE)=NDA(NX1  ,NY  )
            IELM(2,IE)=NDA(NX2+1,NY+1)
            IELM(3,IE)=NDA(NX2  ,NY+1)
            CALL WFSELM(IE,S)
            IF(S.LE.SEPS) IE=IE-1
            NX2=NX2+1
         ENDIF
C
         IF(NX1.LT.NX1MAX) THEN
            IF(NX2.LT.NX2MAX) THEN
               GOTO 20
            ELSE
   21          CONTINUE
               XAI=XDA(NX1, NY   )
               XAJ=XDA(NX2, NY+1 )
               IF(ID1DIV.EQ.1.OR.XAI.LT.XAJ) THEN
                  IF(IE.GE.NELMM) GOTO 9001
                  IE=IE+1
                  IELM(1,IE)=NDA(NX1  ,NY  )
                  IELM(2,IE)=NDA(NX1+1,NY  )
                  IELM(3,IE)=NDA(NX2  ,NY+1)
                  CALL WFSELM(IE,S)
                  IF(S.LE.SEPS) IE=IE-1
                  NX1=NX1+1
                  IF(NX1.LT.NX1MAX) GOTO 21
               ENDIF
            ENDIF
         ELSE
            IF(NX2.LT.NX2MAX) THEN
   22          CONTINUE
               XAI=XDA(NX1,  NY   )
               XAJ=XDA(NX2, NY+1 )
               IF(ID1DIV.EQ.1.OR.XAJ.LT.XAI) THEN
                  IF(IE.GE.NELMM) GOTO 9001
                  IE=IE+1
                  IELM(1,IE)=NDA(NX1   ,NY  )
                  IELM(2,IE)=NDA(NX2+1,NY+1)
                  IELM(3,IE)=NDA(NX2  ,NY+1)
                  CALL WFSELM(IE,S)
                  IF(S.LE.SEPS) IE=IE-1
                  NX2=NX2+1
                  IF(NX2.LT.NX2MAX) GOTO 22
               ENDIF
            ENDIF
         ENDIF
 2000 CONTINUE
C
      NELM=IE
C
C     LOOK FOR UNREFERENCED NODE
C
      DO IN=1,NNOD
         ICK(IN)=0
      ENDDO
      DO IE=1,NELM
         DO J=1,3
            IN=IELM(J,IE)
            ICK(IN)=ICK(IN)+1
         ENDDO
      ENDDO
      IERR=0
      DO IN=1,NNOD
         IF(ICK(IN).EQ.0) THEN
            WRITE(6,*) IN,XD(IN),YD(IN)
            IERR=2
         ENDIF
      ENDDO
      IF(IERR.NE.0) GOTO 9002
C
      RETURN
C
 9001 NELM=IE
      WRITE(6,601) NELM,NELMM
  601 FORMAT(1H ,'SETELM : NELM EXCEEDS NELMM : ',2I8)
      WRITE(6,'(A,3I5,1P2E12.4)') 'NY,NX1,NX2,XAI,XAJ=',
     &                             NY,NX1,NX2,XAI,XAJ
      IERR=1
      RETURN
C
 9002 WRITE(6,602) 
  602 FORMAT(1H ,'SETELM : UNREFERENCED NODE EXISTS')
      IERR=2
      RETURN
C
 9003 WRITE(6,603) 
  603 FORMAT(1H ,'SETELM : FIRST LINE OF ELEMENTS ERROR')
      IERR=3
      RETURN
      END
C
C     ****** List Element Data ******
C
      SUBROUTINE WFLDIV
C
      INCLUDE 'wfcomm.inc'
C
      WRITE(6,100) NNOD,(I,XD(I),YD(I),I=1,NNOD)
  100 FORMAT(1H0,'NODE DATA',7X,'NNOD=',I6/
     &       1H ,9X,'XD',11X,'YD',17X,'XD',11X,'YD'/
     &       (1H ,2(I6,1P2E13.4)))
C
      WRITE(6,200) NELM,(I,(IELM(J,I),J=1,3),
     &                    I=1,NELM)
  200 FORMAT(1H ,'ELEMENT DATA',5X,'NELM=',I6/
     &      (1H ,3(I5,'(',3I5,')',2X)))
C
      RETURN
      END
C
C     ****** Draw Element Data ******
C
      SUBROUTINE WFGDIV
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGVEW(1,1,GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &                GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL PAGES
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,7)
      IF(NDRAWD.EQ.0) THEN
         CALL WFGBDY
      ELSE
         CALL WFGELM
      ENDIF
C
      CALL PAGEE
      RETURN
      END
C
C     ****** Draw Element Paramters ******
C
      SUBROUTINE WFPRME
C
      INCLUDE 'wfcomm.inc'
C
      GXMIN=20.
      GYMAX=17.
      CALL SETCHS(0.3,0.)
      GDY=1.5*0.3
C
      GXL=GXMIN
      GYL=GYMAX
      CALL MOVE(GXL,GYL)
      CALL TEXT('NNOD =',6)
      CALL NUMBI(NNOD,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('ELMT =',6)
      CALL NUMBI(NELM,'(I5)',5)
      RETURN
      END
