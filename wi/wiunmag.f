C
C     *****  INTEGRAL EQUATION  *****
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /FIELD/  CFY(NLEN)
      COMMON /PLASMA/ CWP(0:NXMAX),CWE(0:NXMAX)
      COMMON /KERNEL/ CU(2,-NXMAX:NXMAX)
      COMMON /MATRIX/ CK(NWID,NLEN),CSO(NLEN)
      COMMON /POWER/  CPOWER(0:NXMAX)
      COMMON /ARRYDS/ D0(0:1,0:1),D1(0:1,0:1),
     &                D2(0:1,0:1),D3(0:1,0:1,0:1)
      REAL GT1,GT2
C
      CALL GSOPEN

      NX=400
      NW=200
      XMAX=200.D0
      PN0=2.0D0
      ALFA=-0.005D0
      CFYN=(1.D0,0.D0)
      AKY=0.2D0
      BETA=0.1D0
      CALL INITDS
C
  100 CONTINUE
      CALL GUTIME(GT1)
      WRITE(6,*) '## INPUT NX,NW,XMAX,PN0,ALFA,KY,BETA'
      READ(5,*,END=9000) NX,NW,XMAX,PN0,ALFA,AKY,BETA
      IF(NX.EQ.0) GOTO 9000

      CALL SUBFW(NX,XMAX,PN0,ALFA,NW,AKY,BETA)
      CALL GUTIME(GT2)
      WRITE(6,601) 'SUBFW ',GT2-GT1
      CALL SUBCK2(NX,NW,XMAX,AKY,BETA)
      CALL SUBINI(NX,CFYN,AKY,BETA)
      CALL GUTIME(GT2)
      WRITE(6,601) 'SUBCK2',GT2-GT1
      IF(NW.EQ.NX) THEN
         CALL INVMCD(CK,NX*2+3,IERR)
            IF(IERR.NE.0) GOTO 9900
         CALL SUBFY(NX)
      ELSE
         CALL BANDCD(CK,CSO,2*NX+3,4*NW+3,NWID,IERR)
            IF(IERR.NE.0) GOTO 9900
         CALL SUBFYW(NX)
      ENDIF
      CALL GUTIME(GT2)
      WRITE(6,601) 'SOLVER',GT2-GT1
      CALL SUBPOW(NX,NW,XMAX,AKY,BETA,PTOT)
      CALL GUTIME(GT2)
      WRITE(6,601) 'POWER ',GT2-GT1
      CALL GRA(NX,NW,XMAX,PN0,ALFA,AKY,BETA,PTOT)
      GOTO 100
C
 9000 CALL GSCLOS
      STOP
 9900 WRITE(6,*) 'IERR = ',IERR
      STOP
C
  601 FORMAT(1H ,'## END OF ',A6,' ##  CPU TIME = ',F8.3,' SEC')
      END
C
C     *****  INITIALIZE D0,D1,D2,D3  *****
C
      SUBROUTINE INITDS
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ARRYDS/ D0(0:1,0:1),D1(0:1,0:1),
     &                D2(0:1,0:1),D3(0:1,0:1,0:1)
C
      D0(0,0)=1.D0/3.D0
      D0(0,1)=1.D0/6.D0
      D0(1,0)=1.D0/6.D0
      D0(1,1)=1.D0/3.D0
      D1(0,0)=-0.5D0
      D1(0,1)=-0.5D0
      D1(1,0)=0.5D0
      D1(1,1)=0.5D0
      D2(0,0)=1.D0
      D2(0,1)=-1.D0
      D2(1,0)=-1.D0
      D2(1,1)=1.D0
      DO 1 L=0,1
        DO 2 M=0,1
          DO 3 N=0,1
            D3(L,M,N)=1.D0/12.D0
    3     CONTINUE
    2   CONTINUE
    1 CONTINUE
      D3(0,0,0)=1.D0/4.D0
      D3(1,1,1)=1.D0/4.D0
      RETURN
      END
C
C     *****  SET KERNEL FUNCTION  ***** 
C
      SUBROUTINE SUBFW(NX,XMAX,PN0,ALFA,NW,AKY,BETA)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /PLASMA/ CWP(0:NXMAX),CWE(0:NXMAX)
      COMMON /KERNEL/ CU(2,-NXMAX:NXMAX)
      COMMON /COEF/ G1,G2,G3,G4,G5,N1
C
      DX=XMAX/NX 
      RKY=AKY*BETA
      DO 20 J=1,2
         N1=J
         DO 10 I=0,NW
            X=I*DX
            CALL EUL(X,ALFA,RKY,CS,5,L)
            CU(J, I)=CS
            CU(J,-I)=CS
   10    CONTINUE
   20 CONTINUE
      DO 30 I=0,NX 
         CWE(I)=DEXP(ALFA*I*DX)
         CWP(I)=PN0
   30 CONTINUE
      RETURN
      END
C
C     *****  CALCULATION OF COEFFICIENT MATRIX  ***** 
C
      SUBROUTINE SUBCK2(NX,NW,XMAX,AKY,BETA)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /PLASMA/ CWP(0:NXMAX),CWE(0:NXMAX)
      COMMON /KERNEL/ CU(2,-NXMAX:NXMAX)
      COMMON /MATRIX/ CK(NWID,NLEN),CSO(NLEN)
      COMMON /ARRYDS/ D0(0:1,0:1),D1(0:1,0:1),
     &                D2(0:1,0:1),D3(0:1,0:1,0:1)
C
      DATA CAI/(0.D0,1.D0)/
C
      RKY=AKY*BETA
      RKY2=RKY**2
      DX=XMAX/DBLE(NX)
      DX2=DX**2
      NDUB=2*NX
      BETA2=BETA*BETA
      DKY=AKY*AKY
      CIKY=CAI*AKY/BETA
      CBB=CAI/(DSQRT(1.D0-AKY*AKY)*BETA)
C
      IF(NW.EQ.NX) THEN
         NBAND=0
         NWDUB=NDUB
         NWDDUB=NDUB
      ELSE
         NBAND=1
         NWDUB=2*NW
         NWDDUB=4*NW
      ENDIF
C
      DO 41 I=1,NDUB+3
         DO 42 J=1,NWDDUB+3
            CK(J,I)=(0.,0.)
   42    CONTINUE
   41 CONTINUE
      DO 43 MM=0,NX-1
         DO 44 I=MM,MM+1
            ID=2*I
            DO 45 J=MM,MM+1
               JD=2*J
               IF(NW.NE.NX) JD=2*NW+1+2*J-2*I
               CK(JD+1      ,ID+1)=CK(JD+1      ,ID+1)
     &                              +(DKY-1.D0)*DX*D0(I-MM,J-MM)
               CK(JD+2      ,ID+1)=CK(JD+2      ,ID+1)
     &                              +CIKY*D1(J-MM,I-MM)
               CK(JD+1-NBAND,ID+2)=CK(JD+1-NBAND,ID+2)
     &                              -CIKY*D1(I-MM,J-MM)
               CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2)
     &                              +D2(I-MM,J-MM)/(DX*BETA2)
     &                              -DX*D0(I-MM,J-MM)
   45       CONTINUE
   44    CONTINUE
   43 CONTINUE
      CK(NWDUB+3      ,NDUB+2)=-CBB
      CK(NWDUB+2-NBAND,NDUB+3)=1.D0
      CK(NWDUB+3-NBAND,NDUB+3)=-1.D0
C
      DO 51 MM=0,NX-1
        NS=MM-NW+1
        NE=MM+NW-1
        IF(NS.LE.0) NS=0
        IF(NE.GE.NX-1) NE=NX-1
        DO 52 NN=NS,NE
          DO 53 I=MM,MM+1
            ID=2*I
            DO 54 J=NN,NN+1
              JD=2*J
              IF(NW.NE.NX) JD=2*NW+1+2*J-2*I
              DO 55 KK=MM,MM+1
                DO 56 KD=NN,NN+1
                  CK(JD+1,ID+1)=CK(JD+1,ID+1)
     &                          -CWP(KD)*CWE(KK)*CWE(KD)
     &                          *(DX2*RKY2*CU(1,KK-KD)
     &                          *D0(I-MM,KK-MM)*D0(J-NN,KD-NN)
     &                          +(CU(1,KK-KD)-CAI*CU(2,KK-KD))
     &                          *D1(I-MM,KK-MM)*D1(J-NN,KD-NN))
                  CK(JD+2,ID+1)=CK(JD+2,ID+1)
     &                          -CWP(KD)*CWE(KK)*CWE(KD)
     &                          *(-DX*RKY*CU(2,KK-KD)
     &                          *D0(I-MM,KK-MM)*D1(J-NN,KD-NN))
                  CK(JD+1-NBAND,ID+2)=CK(JD+1-NBAND,ID+2)
     &                          -CWP(KD)*CWE(KK)*CWE(KD)
     &                          *(DX*RKY*CU(2,KK-KD)
     &                          *D1(I-MM,KK-MM)*D0(J-NN,KD-NN))
                  CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2)
     &                          -CWP(KD)*CWE(KK)*CWE(KD)
     &                          *(RKY2*DX2*(CU(1,KK-KD)-CAI*CU(2,KK-KD))
     &                          *D0(I-MM,KK-MM)*D0(J-NN,KD-NN)
     &                          +CU(1,KK-KD)
     &                          *D1(I-MM,KK-MM)*D1(J-NN,KD-NN))
   56           CONTINUE
   55         CONTINUE 
   54       CONTINUE
   53     CONTINUE
   52   CONTINUE 
   51 CONTINUE 
C
      DO 61 MM=0,NX-1
        DO 62 I=MM,MM+1
          ID=2*I
          DO 63 J=MM,MM+1
            JD=2*J
            IF(NW.NE.NX) JD=2*NW+1+2*J-2*I
            DO 64 KS=MM,MM+1
              CK(JD+1      ,ID+1)=CK(JD+1      ,ID+1)
     &                           +CWP(KS)*CWE(KS)*CWE(KS)*DX
     &                           *D3(I-MM,J-MM,KS-MM)
              CK(JD+2-NBAND,ID+2)=CK(JD+2-NBAND,ID+2)
     &                           +CWP(KS)*CWE(KS)*CWE(KS)*DX
     &                           *D3(I-MM,J-MM,KS-MM)
   64       CONTINUE
   63     CONTINUE
   62   CONTINUE
   61 CONTINUE
C
      DO 100 IO=1,NWDDUB+3
         IF(NW.NE.NX) THEN
            IOB=2*NW+4-IO
         ELSE
            IOB=2
         ENDIF
         IF(IOB.GE.1) CK(IOB,IO)=(0.D0,0.D0)
         CK(IO,2)=(0.D0,0.D0)
  100 CONTINUE
      IF(NW.NE.NX) THEN
         I2=2*NW+2
      ELSE
         I2=2
      ENDIF
      CK(I2,2)=(1.D0,0.D0)
      RETURN
      END
C
C     *****  CALCULATION OF RHS VECTOR  *****   
C
      SUBROUTINE SUBINI(NX,CFYN,AKY,BETA)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /MATRIX/ CK(NWID,NLEN),CSO(NLEN)
      DATA CAI/(0.D0,1.D0)/
C
      CBB=CAI/(DSQRT(1.D0-AKY*AKY)*BETA)
C
      DO 10 I=1,NX*2+1
        CSO(I)=(0.D0,0.D0)
   10 CONTINUE
      CSO(NX*2+2)=-CBB*CFYN
      CSO(NX*2+3)=     CFYN
      RETURN
      END
C
C     *****  SET FIELD (FULL MATRIX)  ***** 
C
      SUBROUTINE SUBFY(NX)
C
      IMPLICIT COMPLEX*16(C),REAL(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /FIELD/  CFY(NLEN)
      COMMON /MATRIX/ CK(NWID,NLEN),CSO(NLEN)
C
      DO 10 I=1,NX*2+3
         CFY(I)=(0.D0,0.D0)
         DO 20 J=1,NX*2+3
            CFY(I)=CFY(I)+CK(J,I)*CSO(J)
   20    CONTINUE
   10 CONTINUE
      RETURN
      END
C
C     *****  SET FIELD (BAND MATRIX)  ***** 
C
      SUBROUTINE SUBFYW(NX)
C
      IMPLICIT COMPLEX*16(C),REAL(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /FIELD/  CFY(NLEN)
      COMMON /MATRIX/ CK(NWID,NLEN),CSO(NLEN)
C
      DO 10 I=1,NX*2+3
         CFY(I)=CSO(I)
   10 CONTINUE
      RETURN
      END
C
C     *****  ABSORBED POWER  *****
C
      SUBROUTINE SUBPOW(NX,NW,XMAX,AKY,BETA,PTOT)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /FIELD/  CFY(NLEN)
      COMMON /PLASMA/ CWP(0:NXMAX),CWE(0:NXMAX)
      COMMON /KERNEL/ CU(2,-NXMAX:NXMAX)
      COMMON /POWER/  CPOWER(0:NXMAX)
      COMMON /ARRYDS/ D0(0:1,0:1),D1(0:1,0:1),
     &                D2(0:1,0:1),D3(0:1,0:1,0:1)
      DATA CAI/(0.D0,1.D0)/
C
      RKY=AKY*BETA
      RKY2=RKY**2
      DX=XMAX/DBLE(NX)
      DX2=DX**2
C
      DO 10 N=0,NX
        CPOWER(N)=(0.D0,0.D0)
   10 CONTINUE
      PTOT=0.D0
C
      DO 20 MM=0,NX-1
        NS=MM-NW+1
        NE=MM+NW-1
        AD=1.D0/(2.D0*DX) 
        BD=1.D0/(2.D0*DX) 
        IF(MM.EQ.0) AD=1.D0/DX
        IF(MM.EQ.NX-1) BD=1.D0/DX
        IF(NS.LE.0) NS=0
        IF(NE.GE.NX-1) NE=NX-1
        DO 21 NN=NS,NE
          DO 22 I=MM,MM+1
            ID=2*I
            DO 23 J=NN,NN+1
              JD=2*J
              DO 24 KK=MM,MM+1
                DO 25 KD=NN,NN+1
                  CP1=DX2*RKY2*CU(1,KK-KD)
     &                *D0(I-MM,KK-MM)*D0(J-NN,KD-NN)
     &                +(CU(1,KK-KD)-CAI*CU(2,KK-KD))
     &                *D1(I-MM,KK-MM)*D1(J-NN,KD-NN)
                  CP2=-DX*RKY*CU(2,KK-KD)
     &                *D0(I-MM,KK-MM)*D1(J-NN,KD-NN)
                  CP3= DX*RKY*CU(2,KK-KD)
     &                *D1(I-MM,KK-MM)*D0(J-NN,KD-NN)
                  CP4=DX2*RKY2*(CU(1,KK-KD)-CAI*CU(2,KK-KD))
     &                *D0(I-MM,KK-MM)*D0(J-NN,KD-NN)
     &                +CU(1,KK-KD)*D1(I-MM,KK-MM)*D1(J-NN,KD-NN)
                  CPA=CWP(KD)*CWE(KK)*CWE(KD)
     &                *(CONJG(CFY(ID+1))*(CP1*CFY(JD+1)+CP2*CFY(JD+2)) 
     &                 +CONJG(CFY(ID+2))*(CP3*CFY(JD+1)+CP4*CFY(JD+2)))
                  CPOWER(MM  )=CPOWER(MM  )-CAI*AD*CPA
                  CPOWER(MM+1)=CPOWER(MM+1)-CAI*BD*CPA
                  PTOT=PTOT-REAL(CAI*CPA)
   25           CONTINUE
   24         CONTINUE
   23       CONTINUE
   22     CONTINUE
   21   CONTINUE
   20 CONTINUE
      RETURN
      END
C
C     *****  GRAPHIC   *****
C
      SUBROUTINE GRA(NX,NW,XMAX,PN0,ALFA,AKY,BETA,PTOT)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-R,T-Z),REAL*4(S)
      PARAMETER (NXMAX=401,NWMAX=200,NLEN=NXMAX*2+3,NWID=4*NWMAX+3)
      COMMON /FIELD/  CFY(NLEN)
      COMMON /POWER/  CPOWER(0:NXMAX)
      DIMENSION SVX(NXMAX),SCR(NXMAX),SCI(NXMAX),SCA(NXMAX)
      DIMENSION SCER(NXMAX),SCEI(NXMAX),SCEA(NXMAX),SPOWR(NXMAX) 
C
      DX=XMAX/DBLE(NX)
      DO 20 J=1,NX+1 
         JD=2*(J-1) 
         SVX(J)=(J-1)*XMAX/NX 
         SCR(J)=REAL(CFY(JD+1))
         SCI(J)=AIMAG(CFY(JD+1))
         SCA(J)=ABS(CFY(JD+1))
         SCER(J)=REAL(CFY(JD+2))
         SCEI(J)=AIMAG(CFY(JD+2))
         SCEA(J)=ABS(CFY(JD+2))
         SPOWR(J)=REAL(CPOWER(J-1))
   20 CONTINUE
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(SVX,1,NX+1,1,SMINX,SMAXX)
      CALL GMNMX1(SCA,1,NX+1,1,SMINF,SMAXF)
      CALL GMNMX1(SCEA,1,NX+1,1,SMINE,SMAXE)
      CALL GMNMX1(SPOWR,1,NX+1,1,SMINP,SMAXP)
      CALL GQSCAL(SMINX,SMAXX,SGMINX,SGMAXX,SCALX)
      CALL GQSCAL(-SMAXF,SMAXF,SGMINF,SGMAXF,SCALF)
      CALL GQSCAL(-SMAXE,SMAXE,SGMINE,SGMAXE,SCALE)
      CALL GQSCAL(SMINP,SMAXP,SGMINP,SGMAXP,SCALP)
      CALL GDEFIN(3.,12.,11.0,18.0,SMINX,SMAXX,SGMINF,SGMAXF)
      CALL GFRAME
      CALL GSCALE(0.,SCALX,0.,SCALF,0.3,1)
      CALL GSCALE(0.,10000.,0.,10000.,0.,0)
      CALL GVALUE(0.,2.*SCALX,0.,0.,0)
      CALL GVALUE(0.,0.,0.,2.*SCALF,-2)
      CALL SETLIN(0,0,6)
      CALL GPLOTP(SVX,SCR,1,NX+1,1,0,0,1)
      CALL SETLIN(0,0,5)
      CALL GPLOTP(SVX,SCI,1,NX+1,1,0,0,3)
      CALL SETLIN(0,0,4)
      CALL GPLOTP(SVX,SCA,1,NX+1,1,0,0,5)
      CALL SETLIN(0,0,7)
      CALL GDEFIN(3.,12.,2.0,9.0,SMINX,SMAXX,SGMINE,SGMAXE)
      CALL GFRAME
      CALL GSCALE(0.,SCALX,0.,SCALE,0.3,1)
      CALL GSCALE(0.,10000.,0.,10000.,0.,0)
      CALL GVALUE(0.,2.*SCALX,0.,0.,0)
      CALL GVALUE(0.,0.,0.,2.*SCALE,-2)
      CALL SETLIN(0,0,6)
      CALL GPLOTP(SVX,SCER,1,NX+1,1,0,0,1)
      CALL SETLIN(0,0,5)
      CALL GPLOTP(SVX,SCEI,1,NX+1,1,0,0,3)
      CALL SETLIN(0,0,4)
      CALL GPLOTP(SVX,SCEA,1,NX+1,1,0,0,5)
      CALL SETLIN(0,0,7)
      CALL GDEFIN(16.0,25.,2.0,9.0,SMINX,SMAXX,SGMINP,SGMAXP)
      CALL GFRAME
      CALL GSCALE(0.,SCALX,0.,SCALP,0.3,1)
      CALL GSCALE(0.,10000.,0.,10000.,0.,0)
      CALL GVALUE(0.,2.*SCALX,0.,0.,0)
      CALL GVALUE(0.,0.,0.,2.*SCALP,-2)
      CALL SETLIN(0,0,6)
      CALL GPLOTP(SVX,SPOWR,1,NX+1,1,0,0,1)
      CALL SETLIN(0,0,7)
C
      CALL MOVE(15.5,18.0)
      CALL TEXT('  NX    = ',10)
      CALL NUMBI(NX,'(I7)',8)
      CALL MOVE(15.5,17.5)
      CALL TEXT('  NW    = ',10)
      CALL NUMBI(NW,'(I7)',8)
      CALL MOVE(15.5,17.0)
      CALL TEXT('  XMAX  = ',10)
      CALL NUMBD(XMAX,'(F7.3)',7)
      CALL MOVE(15.5,16.5)
      CALL TEXT('  PN0   = ',10)
      CALL NUMBD(PN0,'(F7.3)',7)
      CALL MOVE(15.5,16.0)
      CALL TEXT('  ALFA  = ',10)
      CALL NUMBD(ALFA,'(F9.5)',9)
      CALL MOVE(15.5,15.5)
      CALL TEXT('  AKY   = ',10)
      CALL NUMBD(AKY,'(F7.3)',7)
      CALL MOVE(15.5,15.0)
      CALL TEXT('  BETA  = ',10)
      CALL NUMBD(BETA,'(F7.3)',7)
      CALL MOVE(15.5,14.0)
      CALL TEXT('  R     = ',10)
      R=ABS(CFY(NX*2+3))**2 
      CALL NUMBD(R,'(F9.5)',9)
      CALL MOVE(15.5,13.5)
      CALL TEXT('  T     = ',10)
      T=1-R
      CALL NUMBD(T,'(F9.5)',9)
      CALL MOVE(15.5,12.5)
      CALL TEXT('  P-IN  = ',10)
      S=T/(DSQRT(1.D0-AKY*AKY)*BETA) 
      CALL NUMBR(S,'(F9.5)',9)
      CALL MOVE(15.5,12.0)
      CALL TEXT('  P-ABS = ',10)
      CALL NUMBD(PTOT,'(F9.5)',9)
      CALL MOVE(15.5,11.0)
      CALL TEXT('   CZ   = ',10)
      CZ=(1.D0+CFY(NX*2+3))*DSQRT(1.D0-AKY*AKY)/(1.D0-CFY(NX*2+3)) 
      S=REAL(CZ)
      CALL NUMBR(S,'(F9.5)',9)
      CALL MOVE(15.5,10.5)
      CALL TEXT('       +i ',10)
      S=IMAG(CZ)
      CALL NUMBR(S,'(F9.5)',9)
      CALL MOVE(6.0,18.1)
      CALL TEXT('< CEX >',7)
      CALL MOVE(6.0,9.1)
      CALL TEXT('< CEY >',7)
      CALL MOVE(19.4,9.1)
      CALL TEXT('< POWER >',9)
      CALL PAGEE
      RETURN
      END
C
C     *****  EULER TRANSFOMATION  *****
C
      SUBROUTINE EUL(X,ALFA,RKY,CS,M,L)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)      
      DIMENSION A(500),B(500)
      COMMON /COEF/ G1,G2,G3,G4,G5,N1
      DATA HP/1.5707 96326 79489 66192D0/
C
      G2=HP
      G3=X
      G4=ALFA
      G5=RKY
      H0=0.5D0
      EPS=1.D-5
      ILST=0
C
      SR1=0.D0
      SI1=0.D0
      IF(M.NE.0) THEN 
         DO 10 K=M-1,0,-1
            G1=DBLE(K)
            CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS,ILST)
            SR1=SR1+SR
            SI1=SI1+SI
   10    CONTINUE
      ENDIF
C
      G1=DBLE(M)
      CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS,ILST)
      A(1)=SR
      SR2=0.5D0*SR
      B(1)=SI
      SI2=0.5D0*SI
      PARITY=-1.D0
      L=0
C
   30 CONTINUE
      L=L+1
      IF(L.GE.200) GOTO 9000
      G1=DBLE(M+L)
      CALL DEFTC2(SR,SI,ESR,ESI,H0,EPS,ILST)
      A(L+1)=SR*PARITY
      B(L+1)=SI*PARITY
      DO 40 K=L,1,-1    
         A(K)=A(K+1)-A(K)
         B(K)=B(K+1)-B(K)
   40 CONTINUE
      SKR=A(1)*PARITY*0.5D0**(L+1)
      SR2=SR2+SKR
      SKI=B(1)*PARITY*0.5D0**(L+1)
      SI2=SI2+SKI
      PARITY=-PARITY
      IF(DABS(SKR).GT.EPS.OR.DABS(SKI).GT.EPS) GOTO 30
C
      SR=(SR1+SR2)/DSQRT(4.D0*G2)
      SI=(SI1+SI2)/DSQRT(4.D0*G2)
      CS=CMPLX(SR,SI)
      RETURN
C
 9000 WRITE(6,*) ' ## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
      RETURN
      END
C
C     *****  REAL PART  *****
C
      REAL*8 FUNCTION FUNR(X,XM,XP)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /COEF/ G1,G2,G3,G4,G5,N1
C
      Y1=XM
      IF(IDINT(G1).EQ.0) THEN 
         T=0.5D0*G2*XP
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            IF(N1.EQ.1) THEN
               AN2=1.D0
            ELSEIF(N1.EQ.2) THEN
               AN2=T
            ENDIF
               FUNR=AN2*0.5*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR=0.D0
         ENDIF 
      ELSE
         T=G2*(X+2.D0*G1)
         T2=T*T
         YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
            IF(N1.EQ.1) THEN
               AN2=1.D0
            ELSEIF(N1.EQ.2) THEN
               AN2=T
            ENDIF
               FUNR=AN2*G2*DEXP(YY)*DCOS(T)
         ELSE
            FUNR=0.D0
         ENDIF
      ENDIF
      RETURN
      END
C
C     *****  IMAG PART  *****
C
      REAL*8 FUNCTION FUNI(X,XM,XP)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /COEF/ G1,G2,G3,G4,G5,N1
C
      Y1=X
      Y2=XM
      T=G2*(XP+2.D0*G1)
      T2=T*T
      YY=-0.5D0*(G3*G3/T2+(G4*G4+G5*G5)*T2)
      IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
         IF(N1.EQ.1) THEN
            AN2=1.D0
         ELSEIF(N1.EQ.2) THEN
            AN2=T
         ENDIF
         FUNI=AN2*G2*DEXP(YY)*DSIN(T)
      ELSE
         FUNI=0.D0
      ENDIF
      RETURN
      END
C
C     *****  DOUBLE EXPONENTIAL FORMULA  *****
C
      SUBROUTINE DEFTC2(CSR,CSI,ESR,ESI,H0,EPS,ILST)
C
C        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
C                    (-1.D0, +1.D0)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER KL*1
C
      DATA HP/1.5707 96326 79489 66192D0/
C
      EPS1=EPS**0.75
      H=H0
      X=0.D0
      CSR=HP*H*FUNR(X,1.D0-X,1.D0+X)
      CSRP=0.D0
      CSI=HP*H*FUNI(X,1.D0-X,1.D0+X)
      CSIP=0.D0
      N=0
      NP=0
      NM=0
      NMIN=1
C
    5 IND=0
      ATPR=1.D0
      ATPI=1.D0
      ATMR=1.D0
      ATMI=1.D0
      ND=2
      EPSI=MAX(EPS1*H,2.D-17)
      IF(N.EQ.0) ND=1
C
   10 N=N+ND
      HN=DBLE(N)*H
      HC=HP*H*COSH(HN)
      IF(IND.NE.1) THEN
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         NP=NP+1
         ATR=ATPR
         ATPR=ABS(CTR)
         ATI=ATPI
         ATPI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATPR)**2+(ATI+ATPI)**2) 
     &           .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
                  IF(IND.EQ.-1) GO TO 100
                  IND=1
            ENDIF
         ENDIF
      ENDIF
C
      IF(IND.NE.-1) THEN
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CTR=HC*FUNR(X,XM,XP)*CC*CC
         CSR=CSR+CTR
         CTI=HC*FUNI(X,XM,XP)*CC*CC
         CSI=CSI+CTI
         NM=NM+1
         ATR=ATMR
         ATMR=ABS(CTR)
         ATI=ATMI
         ATMI=ABS(CTI)
         IF(N.GE.NMIN) THEN
            IF(SQRT((ATR+ATMR)**2+(ATI+ATMI)**2) 
     &           .LT.EPSI*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) THEN 
                  IF(IND.EQ.1) GO TO 100
                  IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10
C
  100 CONTINUE
      ESR=ABS(CSR-CSRP)
      CSRP=CSR
      ESI=ABS(CSI-CSIP)
      CSIP=CSI
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) THEN
            WRITE(6,601) H,NP,NM,CSR
            WRITE(6,604) CSI
         ENDIF
         IF(H.LT.H0) THEN
            WRITE(6,602) H,NP,NM,CSR,ESR
            WRITE(6,605) CSI,ESI
         ENDIF
      ENDIF
      IF(SQRT(ESR*ESR+ESI*ESI)
     &  .LT.EPS1*MAX(SQRT(CSR*CSR+CSI*CSI),1.D0)) GOTO 200
C
      IF(N.GT.1000) THEN
  110    WRITE(6,603)
         READ(5,501,END=120) KL
         IF(KL.EQ.'Q') GOTO 200
         IF(KL.NE.'C') GOTO 110
      ENDIF
  120 CONTINUE
C
      H=0.5D0*H
      CSR=0.5D0*CSR
      CSI=0.5D0*CSI
      NMIN=N/2
      N=-1
      GOTO 5
C
  200 RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  604 FORMAT(1H ,13X,16X,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  605 FORMAT(1H ,13X,16X,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'# SUB DEFTC2 # C or CR : CONTINUE / Q : QUIT')
      END
