C     $Id$
C
C     ***** TASK/DP MENU *****
C
      SUBROUTINE DPMENU
C
      INCLUDE 'dpcomm.inc'
C
      PARAMETER (NXGM=101)
      DIMENSION GX(NXGM),GY(NXGM,6)
C
      DIMENSION CD4(6),CD5(6),CD6(6),CD7(6)
      DIMENSION CLDISP(9)
      CHARACTER KID*1,LINE*80
      EXTERNAL DPPARM
C
    1 CONTINUE
         WRITE(6,*) '## DP MENU: P,V/PARM  ',
     &       'D1,D2,D3/DISP  F/ROOT  T,S/TEST  Q/QUIT'
C
         CALL TASK_KLIN(LINE,KID,MODE,DPPARM)
      IF(MODE.NE.1) GOTO 1
C
      IF(KID.EQ.'P') THEN
         CALL DPPARM(0,'DP',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PLVIEW
         CALL DPVIEW
      ELSEIF(KID.EQ.'D') THEN
         READ(LINE(2:),*,ERR=1,END=1) NID
         IF(NID.EQ.1) THEN
            CALL DPGRP1
         ELSEIF(NID.EQ.2) THEN
            CALL DPCONT
         ELSEIF(NID.EQ.3) THEN
            CALL DPCONTX
         ELSE
            WRITE(6,*) 'XX WRMENU: unknown NID'
         ENDIF
      ELSEIF(KID.EQ.'F') THEN
         CALL DPROOT
      ELSEIF(KID.EQ.'T') THEN
         MODELPS=MODELP(1)
         RFS=RF0
         RKZS=RKZ0
         RKXS=RKX0
         RL=RR
 1001    WRITE(6,*) '# INPUT: RL,RF0,RKZ0,RKX0 ='
         READ(5,*,ERR=1001,END=1002) RL,RF0,RKZ0,RKX0
         IF(RF0.LE.0.D0) GOTO 1002
         CW=2.D0*PI*DCMPLX(RF0,RFI0)*1.D6
         CKPR=MAX(RKZ0,1.D-4)
         CKPP=RKX0
         CALL PLMAG(RL,0.D0,0.D0,RHON)
         CALL PLPROF(RHON)
         MODELP(1)=1
         CALL DPTENS(CW,CKPR,CKPP,1,CD4)
         MODELP(1)=4
         CALL DPTENS(CW,CKPR,CKPP,1,CD5)
         MODELP(1)=5
         CALL DPTENS(CW,CKPR,CKPP,1,CD6)
         MODELP(1)=16
         CALL DPTENS(CW,CKPR,CKPP,1,CD7)
         WRITE(6,602) 
  602    FORMAT(6X,'MODELP=1',13X,'MODELP=4',
     &         13X,'MODELP=5',13X,'MODELP=16')
         WRITE(6,603) (CD4(I),CD5(I),CD6(I),CD7(I),I=1,6)
  603    FORMAT((1PE9.2,1P7E10.2))
         GOTO 1001
C
 1002    MODELP(1)=MODELPS
         RF0=RFS
         RKZ0=RKZS
         RKX0=RKXS
         GOTO 1
      ELSEIF(KID.EQ.'K') THEN
         NS=1
         MM=1
         RF0=0.1D0
         RKPR0=0.1D0
         MM=1
         R1=0.5D0*RA
         DR=0.05D0*RA
         NCH1=1
         NCH2=1
 3001    WRITE(6,*) '## INPUT NS,RF0,RKPR0,MM,R1,DR: '
         READ(5,*,ERR=3001,END=1) NS,RF0,RKPR0,MM,R1,DR
         IF(RF0.EQ.0.D0) GOTO 1
         CW=2.D0*PI*DCMPLX(RF0,RFI0)*1.D6
         CKPR=MAX(RKPR0,1.D-4)
         CALL DPFMFLR(NS,R1,DR,0)
         write(6,*) pmax,npmax,nthmax
         CALL DPFPGRA(PMAX,NPMAX,NTHMAX,PM,THM,FP(1,1,2),'F2')
         NR=2
         CALL DPDKDTR(CW,CKPR,NS,NR,NCH1,NCH2,MM,CLDISP)
         WRITE(6,*) (CLDISP(I),I=1,9)
         GOTO 3001
      ELSEIF(KID.EQ.'S') THEN
         XMIN=-50.D0
         XMAX= 50.D0
         NXGMAX=101
         Q=2.5D0
 2001    WRITE(6,*) '# XMIN,XMAX,NXGMAX,Q='
         READ(5,*,ERR=2001,END=1) XMIN,XMAX,NXGMAX,Q
         DX=(XMAX-XMIN)/(NXGMAX-1)
         DO NX=1,NXGMAX
            X=XMIN+DX*(NX-1)
            GX(NX)=GUCLIP(X)
            CZ=X
            CFZ=CFQZ(Q,CZ)
            GY(NX,1)=GUCLIP(DBLE(CFZ))
            GY(NX,2)=GUCLIP(DIMAG(CFZ))
            CFZ=CFQZ_Z(Q,CZ)
            GY(NX,3)=GUCLIP(DBLE(CFZ))
            GY(NX,4)=GUCLIP(DIMAG(CFZ))
C            CFZ=CFQZ_EXP(Q,CZ)
C            WRITE(6,'(1P3E12.4)') X,CFZ
C            GY(NX,5)=GUCLIP(DBLE(CFZ))
C            GY(NX,6)=GUCLIP(DIMAG(CFZ))
            WRITE(6,'(1P5E12.4)') GX(NX),GY(NX,1),GY(NX,2),GY(NX,3),
     &                                   GY(NX,4)
         ENDDO
         CALL DPGTMP(GX,GY,NXGM,NXGMAX,4)
         GOTO 2001
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ###############################################################
C
C            SUBROUTINE DPGTMP
C
C     ###############################################################
C
      SUBROUTINE DPGTMP(GX,GY,NXM,NXMAX,NGMAX)
C
      DIMENSION GX(NXMAX),GY(NXM,NGMAX)
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY(1,1),1,NXMAX,1,GYMIN,GYMAX)
      DO NG=2,NGMAX
         CALL GMNMX1(GY(1,NG),1,NXMAX,1,GYMIN1,GYMAX1)
         GYMIN=MIN(GYMIN,GYMIN1)
         GYMAX=MAX(GYMAX,GYMAX1)
      ENDDO
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL SETLIN(0,0,5)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GXSMN
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GYSMN
      ENDIF
      CALL GDEFIN(4.,17.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(GXORG,GSCALX,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYORG,GSCALY,0.2,9)
      CALL GVALUE(GXORG,2*GSCALX,0.0,0.0,2)
      CALL GVALUE(0.0,0.0,GYORG,2*GSCALY,-2)
      DO NG=1,NGMAX
         CALL SETLIN(0,0,7-MOD(NG-1,5))
         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,MOD(NG-1,5)+1)
      ENDDO
      CALL PAGEE
      RETURN
      END
