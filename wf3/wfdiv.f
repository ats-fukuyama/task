C     $Id$
C
C     ********** F.E.M. DIVIDER ( FIRST ORDER ) **********
C
      SUBROUTINE WFDIV
C
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
    1 WRITE(6,601) 
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'D') THEN
C
    2    WRITE(6,602) 
         READ(5,'(A1)',ERR=2,END=1) KID
         CALL GUCPTL(KID)
         IF(KID.EQ.'X') THEN
    3       WRITE(6,603) BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
            READ(5,*,ERR=3,END=2) BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
            IDDIV=0
         ELSEIF(KID.EQ.'C') THEN
    4       WRITE(6,604) RB,BZMIN,BZMAX
            READ(5,*,ERR=4,END=2) RB,BZMIN,BZMAX
            BXMIN=-RB
            BXMAX= RB
            BYMIN=-RB
            BYMAX= RB
            IDDIV=1
         ELSEIF(KID.EQ.'A') THEN
    5       WRITE(6,605) RB,RBAX,BZMIN,BZMAX
            READ(5,*,ERR=5,END=2) RB,RBAX,BZMIN,BZMAX
            BXMIN=-RB
            BXMAX= RB
            BYMIN=-RB
            BYMAX= RB
            IDDIV=2
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID: ',KID
         ENDIF
    7    WRITE(6,607) DELX,DELY,DELZ
         READ(5,*,ERR=7,END=2) DELX,DELY,DELZ
         IF(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 2
C
         CALL DFNODE(IERR)
            IF(IERR.NE.0) GOTO 7
         CALL SETNOD(IERR)
            IF(IERR.NE.0) GOTO 7
         IF(IDDIV.NE.2) THEN
            CALL SETELM(IERR)
         ELSE
            CALL SETELMX(IERR)
         ENDIF
            IF(IERR.NE.0) GOTO 7
         WRITE(6,*) '--- WFINDX start ---'
         CALL WFINDX
         WRITE(6,*) '--- WFFEPI start ---'
         CALL WFFEPI
C
         NKMAX=1
         DO NE=1,NEMAX
            KAELM(NE)=1
         ENDDO
         NMKA(1)=0
         NMMAX=0
C
         NBMAX=0
         DO NN=1,NNMAX
            KANOD(NN)=0
         ENDDO
C
      ELSEIF(KID.EQ.'G') THEN
         CALL WFGDIV
         IF(NDRAWD.EQ.-1) THEN
            CALL WFGNAS(0)
            CALL WFGNAS(1)
            CALL WFGNAS(2)
            CALL WFGNAS(3)
            CALL WFGNAS(4)
            CALL WFGNAS(5)
         ENDIF
      ELSEIF(KID.EQ.'W') THEN
         CALL WFLDIV
      ELSEIF(KID.EQ.'L') THEN
         CALL WFRELM(ID)
      ELSEIF(KID.EQ.'P') THEN
         CALL WFPARM(KID)
      ELSEIF(KID.EQ.'V') THEN
         CALL WFVIEW
      ELSEIF(KID.EQ.'S') THEN
         CALL WFWELM(0)
      ELSEIF(KID.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
  601 FORMAT('## INPUT: D/DIV  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',
     &       'W/LIST  X/EXIT')
  602 FORMAT('## TYPE:  X/RECT  C/CIRCLE  A/COAXIAL')
  603 FORMAT('## DIV:   BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX = '/
     &       '          ',6F10.4/
     &       '## INPUT: BXMIN,BXMAX,BYMAX,BYMAX,BZMIN,BZMAX ? ')
  604 FORMAT('## DIV:   RB,BZMIN,BZMAX = ',3F10.4/
     &       '## INPUT: RB,BZMIN,BZMAX ? ')
  605 FORMAT('## DIV:   RB,RBAX,BZMIN,BZMAX = ',4F10.4/
     &       '## INPUT: RB,RBAX,BZMIN,BZMAX ? ')
  607 FORMAT('## DIV:   DELX,DELY,DELZ = ',3F10.4/
     &       '## INPUT: DELX,DELY,DELZ ? ')
      END
C
C     ****** Definititon of Boundary ******
C
      FUNCTION BOUNDF(X,Y)
C
      INCLUDE 'wfcomm.inc'
C
      IF(IDDIV.EQ.0) THEN
         BOUNDF=(X-BXMIN)*(BXMAX-X)*(Y-BYMIN)*(BYMAX-Y)
      ELSEIF(IDDIV.EQ.1) THEN
         BOUNDF=(RB*RB-X*X-Y*Y)
      ENDIF
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
C     ****** Define 2-D Node Array ******
C
      SUBROUTINE DFNODE(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      IF(IDDIV.EQ.0) THEN
         CALL DFNODX(IERR)
      ELSEIF(IDDIV.EQ.1) THEN
         CALL DFNODC(IERR)
      ELSEIF(IDDIV.EQ.2) THEN
         CALL DFNODCX(IERR)
      ENDIF
      RETURN
      END
C
C     ****** Define 2-D Node Array (RECTANGULAR) ******
C
      SUBROUTINE DFNODX(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      NYMAX=NINT((BYMAX-BYMIN)/DELY)+1
      NYMAX=NYMAX+MOD(NYMAX+1,2)
         IF(NYMAX.GT.NYM) GOTO 9200
      DY=(BYMAX-BYMIN)/DBLE(NYMAX-1)
C
      DO NY=1,NYMAX
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
      DO NX=1,NXMAX
         X=DX*(NX-1)+XL(NY)
         XREL=DX*DBLE(NX-(NXMAX+1)/2)
         YREL=DY*DBLE(NY-(NYMAX+1)/2)
         R=SQRT(XREL**2+YREL**2)
         FACTOR=1.D0-R*1.D-6
         XNDA(NX,NY)=X*FACTOR
         YNDA(NX,NY)=Y*FACTOR
      ENDDO
      ENDDO
C
      IERR=0
      RETURN
C
 9200 WRITE(6,602) NYMAX,NYM
  602 FORMAT(' ','DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9300 WRITE(6,603) NY,NXMAX,NXM
  603 FORMAT(' ','DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Define 2-D Node Array (CIRCULAR) ******
C
      SUBROUTINE DFNODC(IERR)
C
      INCLUDE 'wfcomm.inc'
      PARAMETER (NYMH=NYM/2)
      DIMENSION XL1(NYMH),XR1(NYMH)
      DIMENSION NXA1(NYMH),NXR1(NYMH),NXL1(NYMH)
      DIMENSION XNDA1(NXM,NYMH),YNDA1(NXM,NYMH)
      DIMENSION XU(NXM),YU(NXM),DYU(NXM),NYU(NXM)
      EXTERNAL BOUNDX,BOUNDY
      DATA EPS/1.D-8/
C
      NYMAX=NINT(BYMAX/DELY)+1
         IF(NYMAX.GT.NYM) GOTO 9200
      IND=0
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
      DO NY=1,NYMID
         NXA1(NY)=NINT((XR1(NY)-XL1(NY))/DELX)+1
            IF(NXA1(NY).GT.NXM) GOTO 9300
         IF (NXA1(NY).EQ.1) THEN
            NXL1(NY)=1
            NXR1(NY)=0
            XNDA1(1,NY)=XL1(NY)
            YNDA1(1,NY)=DELY*DBLE(NY-1)
            IND=0
            GOTO 1010
         ENDIF
         IF(MOD(NXA1(NY),2).EQ.0) NXA1(NY)=NXA1(NY)+1
         DX=(XR1(NY)-XL1(NY))/DBLE(NXA1(NY)-1)
         NXL1(NY)=1
         NXR1(NY)=1
         DO NX=1,NXA1(NY)
            XNDA1(NX,NY)=DX*DBLE(NX-1)+XL1(NY)
            YNDA1(NX,NY)=DELY*DBLE(NY-1)
         ENDDO
      ENDDO
 1010 CONTINUE
C
      IF(IND.GE.1.AND.NYMID.LT.NYMAX) THEN
         NLMID=NXL1(NYMID)
         NRMID=NXR1(NYMID)
         NXMID=NXA1(NYMID)
         YMID =DELY*DBLE(NYMID-1)
         DO 4100 NX=NLMID+1,NXMID-NRMID
            XU(NX)=XNDA1(NX,NYMID)
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
                  XNDA1(NX,NY)=XNDA1(I,NYMID)
                  YNDA1(NX,NY)=YU(I)
                  IF(INDC.EQ.0) THEN
                     NXL1(NY)=NXL1(NY)+1
                  ELSE
                     NXR1(NY)=NXR1(NY)+1
                  ENDIF
               ELSE IF(NYU(I).GT.NYUL) THEN
                  NX=NX+1
                  XNDA1(NX,NY)=XNDA1(I,NYMID)
                  YNDA1(NX,NY)=DYU(I)*DBLE(NYUL)+YMID
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
      DO NY=1,NYMAX
         NYA=NYMAX+1-NY
         NYB=NYMAX-1+NY
         XL(NYA)=XL1(NY)
         XL(NYB)=XL1(NY)
         XR(NYA)=XR1(NY)
         XR(NYB)=XR1(NY)
         NXA(NYA)=NXA1(NY)
         NXA(NYB)=NXA1(NY)
      DO NX=1,NXA1(NY)
         XNDA(NX,NYA)= XNDA1(NX,NY)
         XNDA(NX,NYB)= XNDA1(NX,NY)
         YNDA(NX,NYA)=-YNDA1(NX,NY)
         YNDA(NX,NYB)= YNDA1(NX,NY)
      ENDDO
      ENDDO
      NYMAX=2*NYMAX-1
C
      IERR=0
      RETURN
C
 9000 WRITE(6,601) NX,NY,YF,XS,XE,XL(NY),XR(NY),ILL
  601 FORMAT(' ','DFNODE : FRGFLS ERROR : NX,NY,YF,XS,XE,XL,XR,ILL'/
     &       ' ',2I5,1P5E15.5,I5)
      IERR=1
      RETURN
C
 9200 WRITE(6,602) NYMAX,NYM
  602 FORMAT(' ','DFNODE : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9300 WRITE(6,603) NY,NXA(NY),NXM
  603 FORMAT(' ','DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
      IERR=1
      RETURN
C
 9400 WRITE(6,604) NX,NY,XF,YS,YE,YU(NX),ILL
  604 FORMAT(' ','NODE1 : FRGFLS ERROR : NX,NY,XF,YS,YE,YU,ILL'/
     &       ' ',2I5,1P4E15.5,I5)
      IERR=1
      RETURN
      END
C
C     ****** Define 2-D Node Array (CIRCULAR COAXIAL) ******
C
      SUBROUTINE DFNODCX(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      NYMAX=NINT((RB-RBAX)/DELY)+1
         IF(NYMAX.GT.NYM) GOTO 9100
         IF(NYMAX.LE.1) GOTO 9100
      DY=(RB-RBAX)/(NYMAX-1)
C
      DO NY=1,NYMAX
         R=RBAX+DY*(NY-1)
         IF(R.LE.0.D0) THEN
            NXA(NY)=1
            XNDA(1,NY)=0.D0
            YNDA(1,NY)=0.D0
         ELSE
            NXMAX=NINT((2.D0*PI*R)/DELX)+1
            IF(NXMAX.GT.NXM) GOTO 9200
            NXMAX=4*((NXMAX-1)/4+1)
            IF(NXMAX.EQ.4) NXMAX=8
            DX=2.D0*PI/NXMAX
            NXA(NY)=NXMAX
            DO NX=1,NXMAX
               TH=DX*(NX-1)
               XNDA(NX,NY)=R*COS(TH)
               YNDA(NX,NY)=R*SIN(TH)
            ENDDO
         ENDIF
C         WRITE(6,'(I5,1PE12.4,I5,1PE12.4)') NY,R,NXMAX,DX
      ENDDO
C
      IERR=0
      RETURN
C
 9100 WRITE(6,601) NYMAX,NYM
  601 FORMAT(' ','DFNODCX : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN
C
 9200 WRITE(6,602) NY,NXMAX,NXM
  602 FORMAT(' ','DFNODCX : NXMAX EXCEEDS NXM AT NX=',I8,' : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Set Node Array ******
C
      SUBROUTINE SETNOD(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      NZMAX=NINT((BZMAX-BZMIN)/DELZ)+1
      IF(NZMAX.GT.NZM) GOTO 9100
      DZ=(BZMAX-BZMIN)/(NZMAX-1)
      DO NZ=1,NZMAX
         ZNDA(NZ)=(NZ-1)*DZ+BZMIN
      ENDDO
C
      IN=0
      INMAX=NNM
      DO NY=1,NYMAX
         DO NX=1,NXA(NY)
            DO NZ=1,NZMAX
               IF(IN.GE.INMAX) GOTO 9000
               IN=IN+1
               NDA(NX,NY,NZ)=IN
               FACTOR=1.D0
C
               XND(IN)= XNDA(NX,NY)*FACTOR
               YND(IN)= YNDA(NX,NY)*FACTOR
               ZND(IN)= ZNDA(NZ)
            ENDDO
         ENDDO
      ENDDO
C
      NNMAX=IN
      IERR=0
      RETURN
C
 9000 WRITE(6,601) IN,INMAX
  601 FORMAT(' ','SETNOD : NNMAX EXCEEDS NNM : ',2I8)
      IERR=1
      RETURN
C
 9100 WRITE(6,602) NZMAX,NZM
  602 FORMAT(' ','SETNOD : NZMAX EXCEEDS NZM : ',2I8)
      IERR=2
      RETURN
      END
C
C     ****** Set Element Array ******
C
      SUBROUTINE SETELM(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      IE=0
C
      DO 2000 NYDO=1,NYMAX-1
         NY=NYDO
         NXMAX=NXA(NY)
         NX1MAX=NXA(NY+1)
         NX=1
         NX1=1
C
         XAI=XNDA(NX,  NY   )
   18    XBJ=XNDA(NX1+1,NY+1)
         IF(XBJ.LT.XAI) THEN
            NX1=NX1+1
            GOTO 18
         ENDIF
         XAJ=XNDA(NX1,  NY+1)
   19    XBI=XNDA(NX+1,NY   )
         IF(XBI.LT.XAJ) THEN
            NX=NX+1
            GOTO 19
         ENDIF
   20    XAI=XNDA(NX,  NY   )
         XBI=XNDA(NX+1,NY   )
         XAJ=XNDA(NX1,  NY+1)
         XBJ=XNDA(NX1+1,NY+1)
         YAI=YNDA(NX,  NY   )
         YBI=YNDA(NX+1,NY   )
         YAJ=YNDA(NX1,  NY+1)
         YBJ=YNDA(NX1+1,NY+1)
         VA=(XAI-XBJ)**2+(YAI-YBJ)**2
         VB=(XBI-XAJ)**2+(YBI-YAJ)**2
         IF(IE.GE.NEM) GOTO 9000
         IF(VA.GT.VB) THEN
            DO NZ=1,NZMAX-1
               CALL SETELL(NX,NX1,NY,NZ,IE,0,IERR)
               IF(IERR.NE.0) GOTO 9000
            ENDDO
            NX=NX+1
         ELSE
            DO NZ=1,NZMAX-1
               CALL SETELL(NX,NX1,NY,NZ,IE,1,IERR)
               IF(IERR.NE.0) GOTO 9000
            ENDDO
            NX1=NX1+1
         ENDIF
         IF(IE.GE.NEM) GOTO 9000
         IF(NX.LT.NXMAX) THEN
            IF(NX1.LT.NX1MAX) THEN
               GOTO 20
            ELSE
   21          CONTINUE
               XAI=XNDA(NX, NY  )
               XAJ=XNDA(NX1,NY+1)
               IF(XAI.LE.XAJ) THEN
                  DO NZ=1,NZMAX-1
                     CALL SETELL(NX,NX1,NY,NZ,IE,0,IERR)
                     IF(IERR.NE.0) GOTO 9000
                  ENDDO
                  NX=NX+1
                  IF(NX.LT.NXMAX) GOTO 21
               ENDIF
            ENDIF
         ELSE
            IF(NX1.LT.NX1MAX) THEN
   22          CONTINUE
               XAI=XNDA(NX,  NY   )
               XAJ=XNDA(NX1, NY+1)
               IF(XAJ.LE.XAI) THEN
                  DO NZ=1,NZMAX-1
                     CALL SETELL(NX,NX1,NY,NZ,IE,1,IERR)
                     IF(IERR.NE.0) GOTO 9000
                  ENDDO
                  NX1=NX1+1
                  IF(NX1.LT.NX1MAX) GOTO 22
               ENDIF
            ENDIF
         ENDIF
 2000 CONTINUE
C
      NEMAX=IE
      IERR=0
      RETURN
C
 9000 WRITE(6,601) IE,NEM
  601 FORMAT(' ','SETELM : NEMAX EXCEEDS NEM : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Set Element Array SUB******
C
      SUBROUTINE SETELL(NX,NX1,NY,NZ,IE,ID,IERR)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION NEL(6),NID(4,3,2)
      DATA NID/1,2,3,6, 1,4,5,6, 1,5,2,6, 2,3,1,6, 2,4,5,6, 2,1,4,6/
C
      IF(IE+2.GE.NEM) GOTO 9000
C
      IF(ID.EQ.0) THEN
         IF(NX+NX1.GE.(NXA(NY)+1)/2+(NXA(NY+1)+1)/2) THEN
            IF(NY.GE.(NYMAX+1)/2) THEN
               NEL(1)=NDA(NX   ,NY  ,NZ  )
               NEL(2)=NDA(NX+1 ,NY  ,NZ  )
               NEL(3)=NDA(NX1  ,NY+1,NZ  )
               NEL(4)=NDA(NX   ,NY  ,NZ+1)
               NEL(5)=NDA(NX+1 ,NY  ,NZ+1)
               NEL(6)=NDA(NX1  ,NY+1,NZ+1)
               IEL=1
            ELSE
               NEL(1)=NDA(NX1  ,NY+1,NZ  )
               NEL(2)=NDA(NX   ,NY  ,NZ  )
               NEL(3)=NDA(NX+1 ,NY  ,NZ  )
               NEL(4)=NDA(NX1  ,NY+1,NZ+1)
               NEL(5)=NDA(NX   ,NY  ,NZ+1)
               NEL(6)=NDA(NX+1 ,NY  ,NZ+1)
               IEL=1
            ENDIF
         ELSE
            IF(NY.GE.(NYMAX+1)/2) THEN
               NEL(1)=NDA(NX   ,NY  ,NZ  )
               NEL(2)=NDA(NX+1 ,NY  ,NZ  )
               NEL(3)=NDA(NX1  ,NY+1,NZ  )
               NEL(4)=NDA(NX   ,NY  ,NZ+1)
               NEL(5)=NDA(NX+1 ,NY  ,NZ+1)
               NEL(6)=NDA(NX1  ,NY+1,NZ+1)
               IEL=2
            ELSE
               NEL(1)=NDA(NX+1 ,NY  ,NZ  )
               NEL(2)=NDA(NX1  ,NY+1,NZ  )
               NEL(3)=NDA(NX   ,NY  ,NZ  )
               NEL(4)=NDA(NX+1 ,NY  ,NZ+1)
               NEL(5)=NDA(NX1  ,NY+1,NZ+1)
               NEL(6)=NDA(NX   ,NY  ,NZ+1)
               IEL=2
            ENDIF
         ENDIF
      ELSE
         IF(NX+NX1.GE.(NXA(NY)+1)/2+(NXA(NY+1)+1)/2) THEN
            IF(NY.GE.(NYMAX+1)/2) THEN
               NEL(1)=NDA(NX1  ,NY+1,NZ  )
               NEL(2)=NDA(NX   ,NY  ,NZ  )
               NEL(3)=NDA(NX1+1,NY+1,NZ  )
               NEL(4)=NDA(NX1  ,NY+1,NZ+1)
               NEL(5)=NDA(NX   ,NY  ,NZ+1)
               NEL(6)=NDA(NX1+1,NY+1,NZ+1)
               IEL=2
            ELSE
               NEL(1)=NDA(NX1+1,NY+1,NZ  )
               NEL(2)=NDA(NX1  ,NY+1,NZ  )
               NEL(3)=NDA(NX   ,NY  ,NZ  )
               NEL(4)=NDA(NX1+1,NY+1,NZ+1)
               NEL(5)=NDA(NX1  ,NY+1,NZ+1)
               NEL(6)=NDA(NX   ,NY  ,NZ+1)
               IEL=2
            ENDIF
         ELSE
            IF(NY.GE.(NYMAX+1)/2) THEN
               NEL(1)=NDA(NX   ,NY  ,NZ  )
               NEL(2)=NDA(NX1+1,NY+1,NZ  )
               NEL(3)=NDA(NX1  ,NY+1,NZ  )
               NEL(4)=NDA(NX   ,NY  ,NZ+1)
               NEL(5)=NDA(NX1+1,NY+1,NZ+1)
               NEL(6)=NDA(NX1  ,NY+1,NZ+1)
               IEL=1
            ELSE
               NEL(1)=NDA(NX1+1,NY+1,NZ  )
               NEL(2)=NDA(NX1  ,NY+1,NZ  )
               NEL(3)=NDA(NX   ,NY  ,NZ  )
               NEL(4)=NDA(NX1+1,NY+1,NZ+1)
               NEL(5)=NDA(NX1  ,NY+1,NZ+1)
               NEL(6)=NDA(NX   ,NY  ,NZ+1)
               IEL=1
            ENDIF
         ENDIF
      ENDIF
      IE=IE+1
      NDELM(1,IE)=NEL(NID(1,1,IEL))
      NDELM(2,IE)=NEL(NID(2,1,IEL))
      NDELM(3,IE)=NEL(NID(3,1,IEL))
      NDELM(4,IE)=NEL(NID(4,1,IEL))
      IE=IE+1
      NDELM(1,IE)=NEL(NID(1,2,IEL))
      NDELM(2,IE)=NEL(NID(2,2,IEL))
      NDELM(3,IE)=NEL(NID(3,2,IEL))
      NDELM(4,IE)=NEL(NID(4,2,IEL))
      IE=IE+1
      NDELM(1,IE)=NEL(NID(1,3,IEL))
      NDELM(2,IE)=NEL(NID(2,3,IEL))
      NDELM(3,IE)=NEL(NID(3,3,IEL))
      NDELM(4,IE)=NEL(NID(4,3,IEL))
C
      IERR=0
      RETURN
C
 9000 IERR=9000
      RETURN
      END
C
C     ****** Set Element Array ******
C
      SUBROUTINE SETELMX(IERR)
C
      INCLUDE 'wfcomm.inc'
C
      IE=0
C
      DO NYDO=1,NYMAX-1
         NY=NYDO
         NX0MAX=NXA(NY)
         NX1MAX=NXA(NY+1)
C
         IF(NX0MAX.EQ.1) THEN
            DO NX10=1,NX1MAX
               NX11=NX10+1
               IF(NX11.GT.NX1MAX) NX11=1
               DO NZ=1,NZMAX-1
                  CALL SETELLX(1,NX10,NX11,NY,NZ,IE,0,IERR)
                  IF(IERR.NE.0) GOTO 9000
               ENDDO
            ENDDO
         ELSE
C
            NX00=1
            NX10=1
C
  100       NX01=NX00+1
            IF(NX01.GT.NX0MAX) NX01=1
            NX11=NX10+1
            IF(NX11.GT.NX1MAX) NX11=1
C            
            RL00=(XNDA(NX11,NY+1)-XNDA(NX00,NY))**2
     &          +(YNDA(NX11,NY+1)-YNDA(NX00,NY))**2
            RL01=(XNDA(NX10,NY+1)-XNDA(NX01,NY))**2
     &          +(YNDA(NX10,NY+1)-YNDA(NX01,NY))**2
C
C            WRITE(6,'(5I5,1P,2E12.4)') 
C     &           NY,NX00,NX10,NX0MAX,NX1MAX,RL00,RL01
C            IF(NX10.GT.NX1MAX) STOP
C
            IF(RL00.LT.RL01) THEN
               DO NZ=1,NZMAX-1
                  CALL SETELLX(NX00,NX10,NX11,NY,NZ,IE,0,IERR)
                  IF(IERR.NE.0) GOTO 9000
               ENDDO
               IF(NX10.EQ.NX1MAX) THEN
                  DO NZ=1,NZMAX-1
                     CALL SETELLX(NX00,NX01,NX11,NY,NZ,IE,1,IERR)
                     IF(IERR.NE.0) GOTO 9000
                  ENDDO
                  GOTO 200
               ELSE
                  NX10=NX10+1
               ENDIF
            ELSEIF(RL00.GE.RL01) THEN
               DO NZ=1,NZMAX-1
                  CALL SETELLX(NX00,NX01,NX10,NY,NZ,IE,1,IERR)
                  IF(IERR.NE.0) GOTO 9000
               ENDDO
               IF(NX00.EQ.NX0MAX) THEN
                  DO NZ=1,NZMAX-1
                     CALL SETELLX(NX01,NX10,NX11,NY,NZ,IE,0,IERR)
                     IF(IERR.NE.0) GOTO 9000
                  ENDDO
                  GOTO 200
               ELSE
                  NX00=NX00+1
               ENDIF
            ENDIF
C            IF(NX01.EQ.1.AND.NX11.EQ.1) THEN
C               IF(IC.EQ.1) GOTO 200
C               IC=IC+1
C               GOTO 100
C            ENDIF
             GOTO 100
C
  200       CONTINUE
         ENDIF
C
C      WRITE(6,*) NY,IE
      ENDDO
C
      NEMAX=IE
      IERR=0
      RETURN
C
 9000 WRITE(6,601) IE,NEM
  601 FORMAT(' ','SETELMX : NEMAX EXCEEDS NEM : ',2I8)
      IERR=1
      RETURN
      END
C
C     ****** Set Element Array SUB******
C
      SUBROUTINE SETELLX(NX0,NX1,NX2,NY,NZ,IE,ID,IERR)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION NEL(6),NID(4,3,4)
      DATA NID/1,2,3,6, 1,4,5,6, 1,5,2,6, 1,2,3,5, 1,6,4,5, 1,3,6,5,
     &         1,2,3,6, 2,4,5,6, 1,4,2,6, 1,2,3,6, 1,4,5,6, 1,5,2,6/
C
      IF(IE+2.GE.NEM) GOTO 9000
C
      IF(ID.EQ.0) THEN
         NEL(1)=NDA(NX0,NY  ,NZ  )
         NEL(2)=NDA(NX1,NY+1,NZ  )
         NEL(3)=NDA(NX2,NY+1,NZ  )
         NEL(4)=NDA(NX0,NY  ,NZ+1)
         NEL(5)=NDA(NX1,NY+1,NZ+1)
         NEL(6)=NDA(NX2,NY+1,NZ+1)
         IEL=1
      ELSE
         NEL(1)=NDA(NX1,NY  ,NZ  )
         NEL(2)=NDA(NX0,NY  ,NZ  )
         NEL(3)=NDA(NX2,NY+1,NZ  )
         NEL(4)=NDA(NX1,NY  ,NZ+1)
         NEL(5)=NDA(NX0,NY  ,NZ+1)
         NEL(6)=NDA(NX2,NY+1,NZ+1)
         IEL=3
      ENDIF
C
      XC=XND(NEL(1))+XND(NEL(2))+XND(NEL(3))
      YC=YND(NEL(1))+YND(NEL(2))+YND(NEL(3))
      IF(XC*YC.LT.0.D0) IEL=IEL+1
C
      IE=IE+1
      NDELM(1,IE)=NEL(NID(1,1,IEL))
      NDELM(2,IE)=NEL(NID(2,1,IEL))
      NDELM(3,IE)=NEL(NID(3,1,IEL))
      NDELM(4,IE)=NEL(NID(4,1,IEL))
      IE=IE+1
      NDELM(1,IE)=NEL(NID(1,2,IEL))
      NDELM(2,IE)=NEL(NID(2,2,IEL))
      NDELM(3,IE)=NEL(NID(3,2,IEL))
      NDELM(4,IE)=NEL(NID(4,2,IEL))
      IE=IE+1
      NDELM(1,IE)=NEL(NID(1,3,IEL))
      NDELM(2,IE)=NEL(NID(2,3,IEL))
      NDELM(3,IE)=NEL(NID(3,3,IEL))
      NDELM(4,IE)=NEL(NID(4,3,IEL))
C
      IERR=0
      RETURN
C
 9000 IERR=9000
      RETURN
      END
C
C     ****** List Element Data ******
C
      SUBROUTINE WFLDIV
C
      INCLUDE 'wfcomm.inc'
C
      WRITE(6,100) NNMAX,(I,XND(I),YND(I),ZND(I),I=1,NNMAX)
  100 FORMAT(/' ','NODE DATA',7X,'NNMAX=',I6/
     &       ' ',8X,'XND',10X,'YND',10X,'ZND'/
     &       (' ',I6,1P3E13.4))
C
      WRITE(6,200) NEMAX,(I,(NDELM(J,I),J=1,4),
     &                    I=1,NEMAX)
  200 FORMAT(/' ','ELEMENT DATA',5X,'NEMAX=',I6/
     &      (' ',2(I5,'(',4I5,')',2X)))
C
      WRITE(6,400) NBMAX,(NDBDY(I),I=1,NBMAX)
  400 FORMAT(' ','BOUNDARY DATA',4X,'NBMAX=',I6/
     &      (' ',10(I6)))
C
      WRITE(6,500) MBND,MLEN,(I,KANOD(I),IMLEN(I),I=1,NNMAX)
  500 FORMAT(' ','NODE AND MATRIX',2X,'MBND=',I6,4X,'MLEN=',I4/
     &      (' ',4(I5,'(',2I5,')')))
      RETURN
      END
