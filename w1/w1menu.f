C
C     ***** TASK/W1 MENU *****
C
      SUBROUTINE W1MENU
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1PRM6/ IELEC(ISM)
      COMMON /W1ANT1/ AJYL(IAM),AJYH(IAM),AJZL(IAM),AJZH(IAM),
     &                ALYL(IAM),ALYH(IAM),APYL(IAM),APYH(IAM)
      COMMON /W1BND1/ CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1ZDAT/ CJ1(NZPM),CJ2(NZPM),CJ3(NZPM),CJ4(NZPM),
     &                ZA(NZPM),AKZ(NZPM),NZANT1(IAM),NZANT2(IAM),NANT
      COMMON /W1CTRL/ DRF,DRKZ,DXFACT,DXWDTH,APRFPN,APRFTR,APRFTP,
     &                NPRINT,NFILE,NGRAPH,NLOOP,NSYM,NFLR,NALPHA,
     &                NSYS,NDISP,NXABS
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
      COMMON /W1PRM5/ EPSH
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,NCDTYP
      CHARACTER KID*1,LINE*80
      EXTERNAL W1PARM

    1 CONTINUE
         IERR=0
         WRITE(6,'(A)')
     &        '## W1 MENU: P,V/PARM  D/DISP  R/RUN  G/GRAPH  Q/QUIT'

         CALL TASK_KLIN(LINE,KID,MODE,W1PARM)
               IF(MODE.NE.1) GOTO 1

         SELECT CASE(KID)
         CASE('P')
            CALL W1PARM(0,'W1',IERR)
         CASE('V')
            CALL W1VIEW
         CASE('D')
            CALL W1GDSP
         CASE('R')
            CALL W1EXEC
         CASE('Q')
            GO TO 9000
         END SELECT
      GO TO 1
 9000 CONTINUE
      RETURN
      END
