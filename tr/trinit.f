C     $Id$
C
C     ***********************************************************
C
C           INITIALIZE CONSTANTS AND DEFAULT VALUES
C
C     ***********************************************************
C
      SUBROUTINE TRINIT
C
      INCLUDE 'trcomm.inc'
C
      NT=0
C
      KUFDEV='jet'
      KUFDCG='19649'
      TIME_INT=0.D0
C
      VOID    = 0.D0
C
      PI      = ASIN(1.D0)*2.D0
      AME     = 9.1093897D-31
      AMM     = 1.6726231D-27
      AEE     = 1.60217733D-19
      VC      = 2.99792458D8
      AMYU0   = 4.D0*PI*1.D-7
      AEPS0   = 1.D0/(VC*VC*AMYU0)
      RKEV    = AEE*1.D3
C
      RR      = 3.0D0
      RA      = 1.2D0
      RKAP    = 1.5D0
      RKAPS   = SQRT(RKAP)
      RDLT    = 0.0D0
      BB      = 3.D0
      RIPS    = 3.D0
      RIPE    = 3.D0
      RIPSS   = 3.D0
      RHOA    = 1.D0
C
      PA(1)   = AME/AMM
      PZ(1)   =-1.D0
      PN(1)   = 0.5D0
      PT(1)   = 1.5D0
      PTS(1)  = 0.05D0
      PNS(1)  = 0.05D0
C
      PA(2)   = 2.D0
      PZ(2)   = 1.D0
      PN(2)   = 0.5D0-2.D-7
      PT(2)   = 1.5D0
      PTS(2)  = 0.05D0
      PNS(2)  = 0.05D0-2.D-8
C
      PA(3)   = 3.D0
      PZ(3)   = 1.D0
      PN(3)   = 1.D-7
      PT(3)   = 1.5D0
      PTS(3)  = 0.05D0
      PNS(3)  = 1.D-8
C
      PA(4)   = 4.D0
      PZ(4)   = 2.D0
      PN(4)   = 1.D-7
      PT(4)   = 1.5D0
      PTS(4)  = 0.05D0
      PNS(4)  = 1.D-8
C
      PA(5)   = 12.D0
      PZ(5)   = 2.D0
      PN(5)   = VOID
      PT(5)   = 0.D0
      PTS(5)  = 0.D0
      PNS(5)  = VOID
C
      PA(6)   = 12.D0
      PZ(6)   = 4.D0
      PN(6)   = VOID
      PT(6)   = 0.D0
      PTS(6)  = 0.D0
      PNS(6)  = VOID
C
      PA(7)   = 2.D0
      PZ(7)   = 0.D0
      PN(7)   = 1.D-15
      PT(7)   = 0.D0
      PTS(7)  = 0.D0
      PNS(7)  = 2.D-4
C
      PA(8)   = 2.D0
      PZ(8)   = 0.D0
      PN(8)   = 1.D-15
      PT(8)   = 0.D0
      PTS(8)  = 0.D0
      PNS(8)  = 1.D-15
C
C      PNC     = 1.D0
C      PNFE    = 1.D-2
      PNC     = 0.D0
      PNFE    = 0.D0
      PNNU    = 0.D0
      PNNUS   = 0.D0
C
      PROFN1 = 2.D0
      PROFN2 = 0.5D0
      PROFT1 = 2.D0
      PROFT2 = 1.D0
      PROFU1 =12.D0
      PROFU2 = 1.D0
      PROFJ1 =-2.D0
      PROFJ2 = 1.D0
C
      ALP(1) = 1.0D0
      ALP(2) = 0.D0
      ALP(3) = 0.D0
C
      AV0    = 0.5D0
      AD0    = 0.5D0
C
      CNC    = 1.D0
      CDW(1) = 0.04D0
      CDW(2) = 0.04D0
      CDW(3) = 0.04D0
      CDW(4) = 0.04D0
      CDW(5) = 0.04D0
      CDW(6) = 0.04D0
      CDW(7) = 0.04D0
      CDW(8) = 0.04D0
C
C        *****  0.GE.MDLKAI.LT.10 : CONSTANT COEFFICIENT MODEL *****
C        ***** 10.GE.MDLKAI.LT.20 : DRIFT WAVE (+ITG +ETG) MODEL *****
C        ***** 20.GE.MDLKAI.LT.30 : REBU-LALLA MODEL *****
C        ***** 30.GE.MDLKAI.LT.40 : CURRENT-DIFFUSIVITY DRIVEN MODEL *****
C        ***** 40.GE.MDLKAI.LT.60 : DRIFT WAVE BALLOONING MODEL *****
C        *****       MDLKAI.GE.60 : ITG(/TEM, ETG) MODEL ETC *****
C
C           *****  MDLKAI.EQ. 0   : CONSTANT *****
C           *****  MDLKAI.EQ. 1   : CONSTANT/(1-A*rho^2) *****
C           *****  MDLKAI.EQ. 2   : CONSTANT*(dTi/drho)^B/(1-A*rho^2) *****
C           *****  MDLKAI.EQ. 3   : CONSTANT*(dTi/drho)^B*Ti^C *****
C                                                                  
C           *****  MDLKAI.EQ. 10  : etac=1 *****
C           *****  MDLKAI.EQ. 11  : etac=1 1/(1+exp) *****
C           *****  MDLKAI.EQ. 12  : etac=1 1/(1+exp) *q *****
C           *****  MDLKAI.EQ. 13  : etac=1 1/(1+exp) *(1+q^2) *****
C           *****  MDLKAI.EQ. 14  : etac=1+2.5*(Ln/RR-0.2) 1/(1+exp) *****
C           *****  MDLKAI.EQ. 15  : etac=1 1/(1+exp) func(q,eps,Ln) *****
C                                                                  
C           *****  MDLKAI.EQ. 20  : Rebu-Lalla model *****
C                                                                  
C           *****  MDLKAI.EQ. 30  : CDBM 1/(1+s) *****
C           *****  MDLKAI.EQ. 31  : CDBM F(s,alpha,kappaq) *****
C           *****  MDLKAI.EQ. 32  : CDBM F(s,alpha,kappaq)/(1+WE1^2) *****
C           *****  MDLKAI.EQ. 33  : CDBM F(s,0,kappaq) *****
C           *****  MDLKAI.EQ. 34  : CDBM F(s,0,kappaq)/(1+WE1^2) *****
C           *****  MDLKAI.EQ. 35  : CDBM (s-alpha)^2/(1+s^2.5) *****
C           *****  MDLKAI.EQ. 36  : CDBM (s-alpha)^2/(1+s^2.5)/(1+WE1^2) *****
C           *****  MDLKAI.EQ. 37  : CDBM s^2/(1+s^2.5) *****
C           *****  MDLKAI.EQ. 38  : CDBM s^2/(1+s^2.5)/(1+WE1^2) *****
C           *****  MDLKAI.EQ. 39  : CDBM F2(s,alpha,kappaq,a/R) *****
C           *****  MDLKAI.EQ. 40  : CDBM F3(s,alpha,kappaq,a/R)/(1+WS1^2) *****
C
C           *****  MDLKAI.EQ. 60  : GLF23 model *****
C           *****  MDLKAI.EQ. 61  : GLF23 (numerical stability enhanced version) *****
C           *****  MDLKAI.EQ. 62  : IFS/PPPL model *****
C           *****  MDLKAI.EQ. 63  : Weiland model *****
C           *****  MDLKAI.EQ. 64  : Bohm/Gyro-Bohm model *****
C
      MDLKAI = 31
      MDLETA = 3
      MDLAD  = 3
      MDLAVK = 3
      MDLJBS = 5
      MDLKNC = 1
C
C     MDLWLD : mode selector
C            0    : using effective transport coefficients
C            else : using transport coefficients' vectors
      MDLWLD=0
C
C     MDDW : mode selector for anomalous particle transport coefficient.
C            you must not modify this parameter.
C            0    : if MDDW=0 from start to finish when you choose
C                   a certain transport model (MDLKAI),
C                   you could control a ratio of anomalous particle
C                   transport to total particle transport to manipulate
C                   the factor of AD0.
C            else : this is because you chose MDLKAI=60, 61, or 63
C                   which assign the transport models that can calculate
C                   an anomalous particle transport coefficient
C                   on their own.
      MDDW=0
C
      DT     = 0.01D0 
      NRMAX  = 50
      NTMAX  = 100
      NTSTEP = 10
      NGRSTP = 100
      NGTSTP = 2
      NGPST  = 4
      TSST   = 1.D9
C
C     *** Convergence Parameter ***
      EPSLTR = 0.001D0
C      EPSLTR = 1.D99
      LMAXTR = 10
C
C     *** Semi-Empirical Parameter for Anomalous Transport ***
      CHP    = 0.D0
      CK0    = 12.D0 ! for electron
      CK1    = 12.D0 ! for ions
      CWEB   = 0.D0  ! for omega ExB
      CALF   = 1.D0  ! for s-alpha
      CKALFA = 0.D0
      CKBETA = 0.D0
      CKGUMA = 0.D0
C
C     *** Sawtooth ***
      TPRST  = 0.1D0
      MDLST  = 0
C
      MDLNF  = 0
      IZERO  = 3
C
C     *** NBI ***
      PNBTOT = 0.D0
      PNBR0  = 0.D0
      PNBRW  = 0.5D0
      PNBENG = 80.D0
      PNBRTG = 3.D0
      MDLNB  = 1
C
C     *** ECH ***
      PECTOT = 0.D0
      PECR0  = 0.D0
      PECRW  = 0.2D0
      PECTOE = 1.D0
      PECNPR = 0.D0
      MDLEC  = 0
C
C     *** LH ***
      PLHTOT = 0.D0
      PLHR0  = 0.D0
      PLHRW  = 0.2D0
      PLHTOE = 1.D0
      PLHNPR = 2.D0
      MDLLH  = 0
C
C     *** ICH ***
      PICTOT = 0.D0
      PICR0  = 0.D0
      PICRW  = 0.5D0
      PICTOE = 0.5D0
      PICNPR = 2.D0
      MDLIC  = 0
C
C     *** Current Drive ***
      PNBCD  = 1.D0
      PECCD  = 0.D0
      PLHCD  = 1.D0
      PICCD  = 0.D0
      PBSCD  = 1.D0
      MDLCD  = 0
C
C     *** Pelet ***
      PELTOT = 0.D0
      PELR0  = 0.D0
      PELRW  = 0.5D0
      PELRAD = 0.D0
      PELVEL = 0.D0
      PELTIM = -10.D0
      MDLPEL = 1
C
      PELPAT(1) = 1.0D0
      PELPAT(2) = 1.0D0
      PELPAT(3) = 0.0D0
      PELPAT(4) = 0.0D0
C
      KFNLOG='trf.log'
C
C     *****
C
C     *** TR or TR/EQ ***
C        0 : TR
C        3 : TR/EQ
      MODELG=0
C     iteration interval for TASK/EQ
      NTEQIT=10
C
C     *** UFILE ***
C        0 : not used
C        1 : time evolution
C        2 : steady state
C        3 : compared with TOPICS
      MDLUF=0
C
C        0 : NSMAX=2, ne=ni
C        1 : calculate nimp and zeff profiles from NE, ZIMP and NM1
C        2 : calculate nimp and ni profiles from NE, ZIMP and ZEFFR
C        3 : calculate zeff and ni profiles from NE, ZIMP and NIMP
      MDNI=0
C
C     Initial profile switch
      MODEP=3
C
C        0 : create AJ(NR) profile from experimental Q profile
C        1 : create QP(NR) profile from experimental CURTOT profile
      MDLJQ=0
C
C     *** Error Indicator for UFILE Reader ***
      MDALL=0
      MDQ=0
      MDCUR=0
      MDTT=0
      MDAREA=0
      MDRMJ=0
      MDRMN=0
      MDGR1=0
      MDGR2=0
C
C     *** Eqs. Selection Parameter ***
      MDLEQB=1  ! 0/1 for B_theta
      MDLEQN=0  ! 0/1 for density
      MDLEQT=1  ! 0/1 for heat
      MDLEQU=0  ! 0/1 for rotation
      MDLEQZ=0  ! 0/1 for impurity
      MDLEQ0=0  ! 0/1 for neutral
      MDLEQE=0  ! 0/1 for electron density
C
      MDLEOI=0  ! 0/1/2 for electron only or bulk ion only if NSMAX=1
C
      NSMAX=2   ! the number of e, D, T and He
      NSZMAX=0  ! the number of impurities
      NSNMAX=2  ! the number of neutrals, 0 or 2 fixed
C      NFMAX=0   ! the number of fast particles
      NSCMAX=NSMAX+NSZMAX ! the number of charged particles
      NSTMAX=NSMAX+NSZMAX+NSNMAX ! the number of all particles
C
C     *** NCLASS SWITCH ***
C        0    : off
C        else : on
      MDNCLS=0
C
C     *** MODERATE TIME EVOLUTION FOR ANOMALOUS TRANSPORT COEFFICIENTS ***
C        0    : off
C        else : multiplyer for TAUK (which is the required time of
C               averaging magnetic surface)
      MDTC=0
C
      CNB=1.D0
C
      RETURN
      END
C
C     ***********************************************************
C
C           PARAMETER INPUT
C
C     ***********************************************************
C
      SUBROUTINE TRPARM(MODE,KIN,IERR)
C
C     MODE=0 : standard namelinst input
C     MODE=1 : namelist file input
C     MODE=2 : namelist line input
C
C     IERR=0 : normal end
C     IERR=1 : namelist standard input error
C     IERR=2 : namelist file does not exist
C     IERR=3 : namelist file open error
C     IERR=4 : namelist file read error
C     IERR=5 : namelist file abormal end of file
C     IERR=6 : namelist line input error
C     IERR=7 : unknown MODE
C     IERR=10X : input parameter out of range
C
      INCLUDE 'trcomm.inc'
C
      EXTERNAL TRNLIN,TRPLST
      CHARACTER KIN*(*)
C
    1 CALL TASK_PARM(MODE,'TR',KIN,TRNLIN,TRPLST,IERR)
      IF(IERR.NE.0) RETURN
C
      CALL TRCHEK(IERR)
      NTMAX_SAVE=NTMAX
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100
C
      RETURN
      END
C
C     ****** INPUT NAMELIST ******
C
      SUBROUTINE TRNLIN(NID,IST,IERR)
C
      INCLUDE 'trcomm.inc'
C
      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RIPSS,RHOA,
     &              PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              PROFJ1,PROFJ2,ALP,AD0,AV0,CNC,CDW,CWEB,CALF,
     &              MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,
     &              DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST,
     &              EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA,
     &              TPRST,
     &              MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDLUF,MDNCLS,MDLWLD,
     &              PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB,
     &              PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC,
     &              PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH,
     &              PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC,
     &              PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD,
     &              PELTOT,PELR0,PELRW,PELRAD,PELVEL,MDLPEL,
     &              PELTIM,PELPAT,KFNLOG,
     &              MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,
     &              MDLEOI,NSMAX,NSZMAX,NSNMAX,
     &              KUFDEV,KUFDCG,TIME_INT,MODEP,MDNI,MDLJQ,MDTC,CNB
C
      IF(NID.LT.0) THEN
         WRITE(-NID,TR,IOSTAT=IST,ERR=9800)
      ELSE
         READ(NID,TR,IOSTAT=IST,ERR=9800,END=9900)
         NTMAX_SAVE=NTMAX
      ENDIF
      IERR=0
      RETURN
C
 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE TRPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &TR : RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RIPSS,RHOA'/
     &       ' ',8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/
     &       ' ',8X,'PNC,PNFE,PNNU,PNNUS'/
     &       ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/
     &       ' ',8X,'PROFJ1,PROFJ2,ALP'/
     &       ' ',8X,'CK0,CK1,CNC,CDW,CWEB,CALF,CKALFA,CKBETA'/
     &       ' ',8X,'KFNLOG,MDLKNC'/
     &       ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/
     &       ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/
     &       ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/
     &       ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,PNBCD,MDLNB'/
     &       ' ',8X,'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/
     &       ' ',8X,'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/
     &       ' ',8X,'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/
     &       ' ',8X,'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/
     &       ' ',8X,'PELTIM,PELPAT,MODELG,NTEQIT,MDLUF,MDNCLS,MDLWLD'/
     &       ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/
     &       ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDEV,KUFDCG'/
     &       ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDTC,CNB')
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE TRCHEK(IERR)
C
      INCLUDE 'trcomm.inc'
C
      IERR=0
C
      IF(NRMAX.LT.1.OR.NRMAX.GT.NRM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX,NRM =',NRMAX,NRM
         NRMAX=NRM
         IERR=1
      ENDIF
C
      IF(NTMAX.LT.1.OR.NTMAX/NGTSTP.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM*NGTSTP
         IERR=1
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           VIEW INPUT PARAMETER
C
C     ***********************************************************
C
      SUBROUTINE TRVIEW(ID)
C
      INCLUDE 'trcomm.inc'
C
      WRITE(6,*) '** TRANSPORT **'
      WRITE(6,602) 'MDLEQB',MDLEQB,
     &             'MDLEQN',MDLEQN,
     &             'MDLEQT',MDLEQT,
     &             'MDLEQU',MDLEQU
      WRITE(6,602) 'MDLEQZ',MDLEQZ,
     &             'MDLEQ0',MDLEQ0,
     &             'MDLEQE',MDLEQE,
     &             'MDLEOI',MDLEOI
      WRITE(6,602) 'NSMAX ',NSMAX,
     &             'NSZMAX',NSZMAX,
     &             'NSNMAX',NSNMAX
      WRITE(6,601) 'RR    ',RR,
     &             'RA    ',RA,
     &             'RKAP  ',RKAP,
     &             'RDLT  ',RDLT
      WRITE(6,601) 'RIPS  ',RIPS,
     &             'RIPE  ',RIPE,
     &             'RIPSS ',RIPSS,
     &             'BB    ',BB
C
      WRITE(6,611)
  611 FORMAT(' ','NS',2X,'PA           PZ    PN(E20)  PNS(E20) ',
     &                   'PT(KEV)  PTS(KEV)  PELPAT')
      DO NS=1,NSTMAX
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS),
     &                PELPAT(NS)
  612    FORMAT(' ',I2,1PD12.4,0P,F6.1,5F9.4)
      ENDDO
C
      WRITE(6,601) 'PNC   ',PNC,
     &             'PNFE  ',PNFE,
     &             'PNNU  ',PNNU,
     &             'PNNUS ',PNNUS
C
      WRITE(6,601) 'PROFN1',PROFN1,
     &             'PROFT1',PROFT1,
     &             'PROFU1',PROFU1,
     &             'PROFJ1',PROFJ1
C
      WRITE(6,601) 'PROFN2',PROFN2,
     &             'PROFT2',PROFT2,
     &             'PROFU2',PROFU2,
     &             'PROFJ2',PROFJ2
C
      WRITE(6,601) 'ALP(1)',ALP(1),
     &             'ALP(2)',ALP(2),
     &             'ALP(3)',ALP(3),
     &             'PBSCD ',PBSCD
C
      WRITE(6,602) 'MDLKAI',MDLKAI,
     &             'MDLETA',MDLETA,
     &             'MDLAD ',MDLAD,
     &             'MDLAVK',MDLAVK
C
      WRITE(6,602) 'MDLJBS',MDLJBS,
     &             'MDLKNC',MDLKNC,
     &             'MDNCLS',MDNCLS,
     &             'MDLWLD',MDLWLD
C
      WRITE(6,604) 'MDLUF ',MDLUF,
     &             'KUFDEV',KUFDEV,
     &             'KUFDCG',KUFDCG,
     &             'MDNI  ',MDNI
      WRITE(6,605) 'MDLJQ ',MDLJQ,
     &             'MDTC  ',MDTC,
     &             'RHOA  ',RHOA
C
      WRITE(6,602) 'MODELG',MODELG,
     &             'NTEQIT',NTEQIT
C
      WRITE(6,601) 'CK0   ',CK0,
     &             'CK1   ',CK1,
     &             'CNC   ',CNC,
     &             'AD0   ',AD0
C
      WRITE(6,601) 'CHP   ',CHP,
     &             'CWEB  ',CWEB,
     &             'CALF  ',CALF
C
      IF((MDLKAI.GE.1.AND.MDLKAI.LT.10).OR.ID.EQ.1)
     &   WRITE(6,601) 'CKALFA',CKALFA,
     &                'CKBETA',CKBETA,
     &                'CKGUMA',CKGUMA
C
      IF((MDLKAI.GE.10.AND.MDLKAI.LT.20).OR.ID.EQ.1) 
     &   WRITE(6,613) CDW(1),CDW(2),CDW(3),CDW(4),
     &                CDW(5),CDW(6),CDW(7),CDW(8)
  613 FORMAT(' ','    AKDW(E) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       ' ','    AKDW(D) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       ' ','    AKDW(T) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       ' ','    AKDW(A) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW')
C
      WRITE(6,601) 'DT    ',DT,
     &             'EPSLTR',EPSLTR,
     &             'TSST  ',TSST,
     &             'TPRST ',TPRST
      WRITE(6,602) 'LMAXTR',LMAXTR,
     &             'NRMAX ',NRMAX,
     &             'NTMAX ',NTMAX,
     &             'NTSTEP',NTSTEP
      WRITE(6,602) 'NGRSTP',NGRSTP,
     &             'NGTSTP',NGTSTP,
     &             'NGPST ',NGPST,
     &             'IZERO ',IZERO
C
      WRITE(6,602) 'MDLST ',MDLST,
     &             'MDLCD ',MDLCD,
     &             'MDLNF ',MDLNF
C
      IF((PNBTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PNBTOT',PNBTOT,
     &                'PNBR0 ',PNBR0,
     &                'PNBRW ',PNBRW,
     &                'PNBENG',PNBENG
         WRITE(6,603) 'MDLNB ',MDLNB,
     &                'PNBRTG',PNBRTG,
     &                'PNBCD ',PNBCD
      ENDIF
C
      IF((PECTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PECTOT',PECTOT,
     &                'PECR0 ',PECR0,
     &                'PECRW ',PECRW,
     &                'PECTOE',PECTOE
         WRITE(6,603) 'MDLEC ',MDLEC,
     &                'PECNPR',PECNPR,
     &                'PECCD ',PECCD
      ENDIF
C
      IF((PLHTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PLHTOT',PLHTOT,
     &                'PLHR0 ',PLHR0,
     &                'PLHRW ',PLHRW,
     &                'PLHTOE',PLHTOE
         WRITE(6,603) 'MDLLH ',MDLLH,
     &                'PLHNPR',PLHNPR,
     &                'PLHCD ',PLHCD
      ENDIF
C
      IF((PICTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PICTOT',PICTOT,
     &                'PICR0 ',PICR0,
     &                'PICRW ',PICRW,
     &                'PICTOE',PICTOE
         WRITE(6,603) 'MDLIC ',MDLIC,
     &                'PICNPR',PICNPR,
     &                'PICCD ',PICCD
      ENDIF
C
      IF((PELTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PELTOT',PELTOT,
     &                'PELR0 ',PELR0,
     &                'PELRW ',PELRW
         WRITE(6,603) 'MDLPEL',MDLPEL,
     &                'PELRAD',PELRAD,
     &                'PELVEL',PELVEL,
     &                'PELTIM',PELTIM
      ENDIF
C
      WRITE(6,601) 'CNB   ',CNB 
      RETURN
C
  601 FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X   :2X,A6,'=',I7)
  603 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1X,A6,4X:
     &        2X,A6,'=',1X,A6,4X:2X,A6,'=',I7)
  605 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
      END
C
C     ***********************************************************
C
C           SET INITIAL PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRPROF
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC2/ NTAMAX
C
C     ZEFF=1
C
C      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
C      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
C      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
C      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
C      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
C      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
      FACTJ   = 1.D0
C
      T     = 0.D0
      TPRE  = 0.D0
      TST   = 0.D0
      VSEC  = 0.D0
      NGR   = 0
      NGT   = 0
      NGST  = 0
      RIP   = RIPS
C
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
         RR=RRU(1)
         RA=RAU(1)
         RKAP=RKAPU(1)
         BB=BBU(1)
      ENDIF
      RKAPS=SQRT(RKAP)
C
      CALL TR_EDGE_DETERMINER(0)
      CALL TR_EDGE_SELECTOR(0)
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      DO NR=1,NRMAX
         RG(NR) = DBLE(NR)*DR
         RM(NR) =(DBLE(NR)-0.5D0)*DR
C
         IF(MDLUF.EQ.1) THEN
            IF(MDNI.EQ.0) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,2
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ELSE
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = RNU(1,NR,3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
C               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
C                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
C                  DO NS=1,3
C                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
C                  ENDDO
C               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
C                  RT(NR,1) = RTU(1,NR,1)
C                  RT(NR,2) = RTU(1,NR,2)
C                  RT(NR,3) = RTU(1,NR,3)
C               ENDIF
               RT(NR,1) = RTU(1,NR,1)
               RT(NR,2) = RTU(1,NR,2)
               RT(NR,3) = RTU(1,NR,3)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ENDIF
C
            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0
C
            PBM(NR)  = PBMU(1,NR)
            RNFS(NR) = RNFU(1,NR)
         ELSEIF(MDLUF.EQ.2) THEN
            IF(MDNI.EQ.0) THEN !!!
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  DO NS=1,2
                     RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RN(NR,1) = RNU(1,NR,1)
                  RN(NR,2) = RNU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,3) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF
     &                    +RNU(1,NRMAX,2)
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF
     &                    +RNU(1,NRMAX,2)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,2
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               RT(NR,3) = RTU(1,NR,2)
               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,2
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ENDIF
            ELSE !!!
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  DO NS=1,3
                     RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  DO NS=1,3
                     RN(NR,NS) = RNU(1,NR,NS)
                  ENDDO
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF
     &                    +RNU(1,NRMAX,2)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,3
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  DO NS=1,3
                     RT(NR,NS) = RTU(1,NR,NS)
                  ENDDO
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               DO NS=1,3
                  RN(NR,NS) = RNU(1,NR,NS)
               ENDDO
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               DO NS=1,3
                  RN(NR,NS) = RNU(1,NR,NS)
               ENDDO
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
C
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  DO NS=1,3
                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
                  ENDDO
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  DO NS=1,3
                     RT(NR,NS) = RTU(1,NR,NS)
                  ENDDO
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF
     &                    +RTU(1,NRMAX,2)
            ENDIF
            ENDIF !!!
C
            PEX(NR,1)=PNBU(1,NR,1)
            PEX(NR,2)=PNBU(1,NR,2)
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0
C
            BOGUS=0.D0
            SEX(NR,1)=SNBU(1,NR,1)*BOGUS
            SEX(NR,2)=SNBU(1,NR,2)*BOGUS
            SEX(NR,3)=0.D0
            SEX(NR,4)=0.D0
         ELSEIF(MDLUF.EQ.3) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            DO NS=1,NSM
               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
            ENDDO
C
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            DO NS=1,NSM
               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
            ENDDO
C
            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0
         ELSE
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            DO NS=1,NSM
               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
            ENDDO
C
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            DO NS=1,NSM
               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
            ENDDO
C
            DO NS=1,NSM
               PEX(NR,NS) = 0.D0
            ENDDO
         ENDIF
C
C     *****
C
C         PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
C         DO NS=1,NSM
C            RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
C         ENDDO
C
C         PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
C         DO NS=1,NSM
C            RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
C         ENDDO
C
C     *****
C     
         IF(MDLEQ0.EQ.1) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1)**PROFU2
            RN(NR,7) = (PN(7)-PNS(7))*PROF+PNS(7)
            RN(NR,8) = (PN(8)-PNS(8))*PROF+PNS(8)
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         ENDIF
C
         DO NF=1,NFM
           RW(NR,NF) = 0.D0
         ENDDO
      ENDDO
      CALL PLDATA_SETR(RG,RM)
      CALL TR_EDGE_DETERMINER(1)
      CALL TR_EDGE_SELECTOR(1)
C
C     *** CALCULATE GEOMETRIC FACTOR ***
C
      CALL TRSTGF
      CALL TRGFRG
C
C     *** CALCULATE PZC,PZFE ***
C
      CALL TRZEFF
C
C     *** CALCULATE ANEAVE ***
C
      ANESUM=0.D0
      DO NR=1,NRMAX
         ANESUM=ANESUM+RN(NR,1)*RM(NR)
      ENDDO 
      ANEAVE=ANESUM*2.D0*DR
C
C     *** CALCULATE IMPURITY DENSITY
C                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***
C
      IF(MDLUF.NE.3) THEN
         DO NR=1,NRMAX
            ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC
     &               *1.D-2*RN(NR,1)
            ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE
     &               *1.D-2*RN(NR,1)
            ANI = 0.D0
            DO NS=2,NSM
               ANI=ANI+PZ(NS)*RN(NR,NS)
            ENDDO
            ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
            DILUTE = 1.D0-ANZ/ANI
            DO NS=2,NSM
               RN(NR,NS) = RN(NR,NS)*DILUTE
            ENDDO
         ENDDO
         PNSS(1)=PNS(1)
         DO NS=2,NSM
            PNSS(NS)=PNS(NS)*DILUTE
         ENDDO
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1)=PNSA(1)
            DO NS=2,NSM
               PNSSA(NS)=PNSA(NS)*DILUTE
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ELSE
         DO NS=1,NSM
            PNSS(NS)=PNS(NS)
         ENDDO
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            DO NS=1,NSM
               PNSSA(NS)=PNSA(NS)
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ENDIF
C
C     *** CALCULATE PROFILE OF AJ(R) ***
C
C     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***
C
      IF(MDLUF.EQ.1) THEN
         IF(MDLJQ.EQ.0) THEN
         NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         DO NR=1,NRMAX
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         ENDDO
         ELSEIF(MDLJQ.EQ.1) THEN
         DO NR=1,NRMAX
            QP(NR) =QPU(1,NR)
            RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
         NR=1
            FACTOR0=TTRHO(NR)**2/(AMYU0*BB*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
            AJOH(NR)=AJ(NR)
         DO NR=2,NRMAX
            FACTOR0=TTRHO(NR)**2/(AMYU0*BB*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
            AJOH(NR)=AJ(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         ENDIF
C
         RIP   = RIPU(1)
         RIPS  = RIPU(1)
         RIPSS = RIPU(1)
         RIPE  = RIPU(1)
      ELSEIF(MDLUF.EQ.2) THEN
         IF(MDLJQ.EQ.0) THEN
            NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         DO NR=1,NRMAX
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*RDP(NR))
         ENDDO
C
         ELSE
C
         DO NR=1,NRMAX
            QP(NR) =QPU(1,NR)
            RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
            NR=1
               FACTOR0=TTRHO(NR)**2/(AMYU0*BB*DVRHO(NR))
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
               AJOH(NR)=AJ(NR)
            DO NR=2,NRMAX
               FACTOR0=TTRHO(NR)**2/(AMYU0*BB*DVRHO(NR))
               FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
               AJOH(NR)=AJ(NR)
            ENDDO
            NR=1
               FACTOR0=RR/(AMYU0*DVRHO(NR))
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
               AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
            DO NR=2,NRMAX
               FACTOR0=RR/(AMYU0*DVRHO(NR))
               FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
               AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
            ENDDO
         ENDIF
C
C         RIP   = 2.D0*PI*RA*RKAPS*BP(NRMAX)/AMYU0*1.D-6
         RIP   = RIPS
         RIPSS = RIP
         RIPE  = RIP
      ELSEIF(MDLUF.EQ.3) THEN
         DO NR=1,NRMAX
            AJOH(NR)=AJU(1,NR)
            AJ(NR)  =AJU(1,NR)
         ENDDO
C
         NR=1
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
         DO NR=2,NRMAX
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
         ENDDO
         DO NR=1,NRMAX
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
         RDPS=2.D0*PI*AMYU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         DO NR=1,NRMAX
            AJOH(NR)=FACT*AJOH(NR)
            AJ(NR)  =AJOH(NR)
            BP(NR)  =FACT*BP(NR)
            QP(NR)  =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &              /(4.D0*PI**2*RDP(NR))
         ENDDO
      ELSE
         DO NR=1,NRMAX
            IF((1.D0-RM(NR)**ABS(PROFJ1)).LE.0.D0) THEN
               PROF=0.D0    
            ELSE             
               PROF= (1.D0-RM(NR)**ABS(PROFJ1))**ABS(PROFJ2)
            ENDIF             
            AJOH(NR)= PROF
            AJ(NR)  = PROF
         ENDDO
C
         NR=1
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
C
         RDPS=2.D0*PI*AMYU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         DO NR=1,NRMAX
            RDP(NR)=FACT*RDP(NR)
            AJOH(NR)=FACT*AJOH(NR)
            AJ(NR)  =AJOH(NR)
            BP(NR)  =FACT*BP(NR)
            QP(NR)  =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &              /(4.D0*PI**2*RDP(NR))
         ENDDO
      ENDIF
C      Q0=(4.D0*QP(1)-QP(2))/3.D0
C
C     *** THIS MODEL ASSUMES CONSTANT EZ ***
C
      IF(PROFJ1.LE.0.D0.OR.MDNCLS.EQ.1) THEN
         CALL TRZEFF
         DO NR=1,NRMAX
C
C        ****** CLASSICAL RESISTIVITY ******
C
            ANE=RN(NR,1)
            TEL =ABS(RT(NR,1))
            ZEFFL=ZEFF(NR)
C
            COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &           /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
            TAUE = COEF*SQRT(AME)*(TEL*RKEV)**1.5D0/SQRT(2.D0)
C
            ETA(NR) = AME/(ANE*1.D20*AEE*AEE*TAUE)
     &              *(0.29D0+0.46D0/(1.08D0+ZEFFL))
C
C        ****** NEOCLASSICAL RESISTIVITY ******
C
            IF(MDLETA.EQ.1) THEN
               EPS=RA*RM(NR)/RR
               EPSS=SQRT(EPS)**3
               IF(NR.EQ.1) THEN
                  QL= 0.25D0*(3.D0*Q0+QP(NR))
               ELSE
                  QL= 0.5D0*(QP(NR-1)+QP(NR))
               ENDIF
               VTE=SQRT(TEL*RKEV/AME)
               RNUE=QL*RR/(TAUE*VTE*EPSS)
               RK33E=RK33/(1.D0+RA33*SQRT(ABS(RNUE))+RB33*RNUE)
     &                   /(1.D0+RC33*RNUE*EPSS)
C
               FT     = 1.D0-SQRT(EPS)*RK33E
               ETA(NR)= ETA(NR)/FT
C     
C        ****** NEOCLASSICAL RESISTIVITY PART II ******
C
            ELSEIF(MDLETA.EQ.2) THEN
               EPS=RM(NR)*RA/RR
               EPSS=SQRT(EPS)**3
               IF(NR.EQ.1) THEN
                  Q0=(4.D0*QP(1)-QP(2))/3.D0
                  QL= 0.25D0*(3.D0*Q0+QP(NR))
                  ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))
               ELSE
                  QL= 0.5D0*(QP(NR-1)+QP(NR))
                  ZEFFL=0.5D0*(ZEFF(NR-1)+ZEFF(NR))
               ENDIF
C
C               VTE=1.33D+7*DSQRT(TEL)
               VTE=SQRT(ABS(TEL)*RKEV/AME)
               FT=1.D0-(1.D0-EPS)**2
     &         /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
               rLnLam=15.2D0-0.5D0*DLOG(ANE)+DLOG(TEL)
               TAUE=6.D0*PI*SQRT(2.D0*PI)*AEPS0**2*DSQRT(AME)
     &             *(TEL*RKEV)**1.5D0/(ANE*1.D20*AEE**4*rLnLam)
               RNUSE=RR*QL/(VTE*TAUE*EPSS)
               PHI=FT/(1.D0+(0.58D0+0.2D0*ZEFFL)*RNUSE)                
               ETAS=1.65D-9*rLnLam/(ABS(TEL)**1.5D0)
               CH=0.56D0*(3.D0-ZEFFL)/((3.D0+ZEFFL)*ZEFFL)
C
               ETA(NR)=ETAS*ZEFFL*(1.D0+0.27D0*(ZEFFL-1.D0))
     &                           /((1.D0-PHI)*(1.D0-CH*PHI)
     &                           *(1.D0+0.47D0*(ZEFFL-1.D0)))
C
C        ****** NEOCLASSICAL RESISTIVITY BY O. SAUTER  ******
C
            ELSEIF(MDLETA.EQ.3) THEN
               EPS=RM(NR)*RA/RR
               EPSS=SQRT(EPS)**3
               IF(NR.EQ.1) THEN
                  Q0=(4.D0*QP(1)-QP(2))/3.D0
                  QL= 0.25D0*(3.D0*Q0+QP(NR))
                  ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))
               ELSE
                  QL= 0.5D0*(QP(NR-1)+QP(NR))
                  ZEFFL=0.5D0*(ZEFF(NR-1)+ZEFF(NR))
               ENDIF
               rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/(TEL*1.D3))
               RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
               SGMSPTZ=1.9012D4*(TEL*1.D3)**1.5D0/(ZEFFL*RNZ*rLnLame)
               FT=1.D0-(1.D0-EPS)**2.D0
     &          /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
               RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame
     &             /((TEL*1.D3)**2*EPSS)
               F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)
     &                +0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5D0)
               ETA(NR)=1.D0/(SGMSPTZ*F33(F33TEFF,ZEFFL))
            ENDIF
         ENDDO
         IF(PROFJ1.GT.0.D0.AND.MDNCLS.EQ.1) GOTO 2000
C
         DO NR=1,NRMAX
            AJOH(NR)=1.D0/ETA(NR)
            AJ(NR)  =1.D0/ETA(NR)
         ENDDO
C
         NR=1
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
         DO NR=2,NRMAX
            FACTOR0=AMYU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
         ENDDO
         NR=1
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(AMYU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         DO NR=1,NRMAX
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
C
         RDPS=2.D0*PI*AMYU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         DO NR=1,NRMAX
            RDP(NR)  =FACT*RDP(NR)
            AJOH(NR) =FACT*AJOH(NR)
            AJ(NR)   =AJOH(NR)
            AJTOR(NR)=FACT*AJTOR(NR)
            BP(NR)   =FACT*BP(NR)
            QP(NR)   =TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &               /(4.D0*PI**2*RDP(NR))
         ENDDO
         Q0=(4.D0*QP(1)-QP(2))/3.D0
      ENDIF
 2000 CONTINUE
      SUM=0.D0
      DO NR=1,NRMAX
         SUM=SUM+RDP(NR)*DR
         RPSI(NR)=SUM
         BPRHO(NR)=BP(NR)
         QRHO(NR)=QP(NR)
      ENDDO
C
      IF(MODELG.EQ.3) THEN
         DO NR=1,NRMAX
            RHOTR(NR)=RM(NR)
            AJ(NR)   =AJOH(NR)
            HJRHO(NR)=AJ(NR)*1.D-6
         ENDDO
         CALL TREQIN(RR,RA,RKAP,RDLT,BB,RIP,
     &               NRMAX,RHOTR,HJRHO,QRHO,MDLUF,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TREQIN1: IERR=',IERR
C
         DO NR=1,NRMAX
C            WRITE(6,'(A,I3,1P5E12.4)') 'NR,R,AJ,E,Q,QP=',
C     &           NR,RHOTR(NR),AJ(NR),ETA(NR),QRHO(NR),QP(NR)
C            WRITE(6,'(A,I3,1P5E12.4)') 'NR,A,K,RDLT,B,I=',
C     &           NR,RA,RKAP,RDLT,BB,RIP
            QP(NR)=-QRHO(NR+1)
         ENDDO
         Q0=(4.D0*QP(1)-QP(2))/3.D0
      ENDIF
C
      DRIP  = (RIPSS-RIPS)/1.D1
 1000 RIP   = RIPSS
      IF(MODELG.EQ.3) THEN
C         CALL TRCONV(L,IERR)
C         WRITE(6,*) "L=",L
         CALL TRSETG
      ENDIF
C
      GRG(1)=0.0
      DO NR=1,NRMAX
         GRM(NR)  =GUCLIP(RM(NR))
         GRG(NR+1)=GUCLIP(RG(NR))
      ENDDO
C
      IF(MODELG.NE.0) THEN
         RIPSS=RIPSS-DRIP
C         write(6,'(A,1P4E12.5)') "RIP,RIPSS,RIPS,RIPE= ",RIP,RIPSS,RIPS
C     &        ,RIPE
         IF(DRIP.NE.0.AND.ABS(RIPSS-RIPS).LT.1.D-10) THEN
            GOTO 1000
         ENDIF
      ENDIF
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      RETURN
      END
C
C     ***********************************************************
C
C           SET GEOMETRICAL FACTOR VIA TASK/EQ
C
C     ***********************************************************
C
      SUBROUTINE TRSETG
C
      INCLUDE 'trcomm.inc'
C
      DO NR=1,NRMAX
         PRHO(NR)=0.D0
         TRHO(NR)=0.D0
         DO NS=1,NSM
            PRHO(NR)=PRHO(NR)+RN(NR,NS)*RT(NR,NS)*1.D14*RKEV
         ENDDO
         DO NF=1,NFM
            PRHO(NR)=PRHO(NR)+RW(NR,NF)*1.D14*RKEV
         ENDDO
         TRHO(NR)=0.D0
         DO NS=2,NSM
            TRHO(NR)=TRHO(NR)+RT(NR,NS)*RN(NR,NS)/RN(NR,1)
         ENDDO
         HJRHO(NR)=AJ(NR)*1.D-6 !*FACTJ
         VTRHO(NR)=0.D0
C         WRITE(6,'(A,I5,1P4E12.4)')
C     &           'NR,P/HJ/T/VTRHO=',NR,PRHO(NR),
C     &           HJRHO(NR),TRHO(NR),VTRHO(NR)
      ENDDO
C      PAUSE
C
C      DO NR=1,NRMAX
C         WRITE(6,'(A,2I5,1P4E12.4)')
C     &           'NR,I/RM/J/V/T=',NR,NRMAX,RIP,
C     &           HJRHO(NR),VTRHO(NR),TRHO(NR)
C      ENDDO
C
      CALL TREQEX(RIP,NRMAX,PRHO,HJRHO,VTRHO,TRHO,
     &            QRHO,TTRHO,DVRHO,DSRHO,
     &            ABRHO,ARRHO,AR1RHO,AR2RHO,
     &            EPSRHO,MDLUF,IERR)
C
C      DO NR=1,NRMAX
C         WRITE(6,'(A,I5,1P4E12.4)')
C     &           'NR,Q/TT/DV/AB=',NR,QRHO(NR),
C     &           TTRHO(NR),DVRHO(NR),ABRHO(NR)
C         WRITE(6,'(A,I5,1P4E12.4)')
C     &           'NR,Q/TT/AB/EP=',NR,QRHO(NR),
C     &           TTRHO(NR),ABRHO(NR),EPSRHO(NR)
C      ENDDO
C      PAUSE
C
C      DO NR=1,NRMAX
C         BPNR=RA*RHOTR(NR)*TTRHO(NR)*ARRHO(NR)/QRHO(NR)
C         write(6,*) NR,BPNR,BP(NR)
C         BP(NR)=BPNR
C      ENDDO
C
      DO NR=1,NRMAX
         BPRHO(NR)=RKAPS*RA*RG(NR)*BB/(RR*QRHO(NR))
      ENDDO
C
      NR=1
         FACTOR0=RR*RR/(AMYU0*DVRHO(NR))
         FACTOR2=DVRHO(NR  )*ABRHO(NR  )
         FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)
         FACTORP=0.5D0*(FACTOR2+FACTOR3)
         AJ(NR)= FACTOR0*FACTORP*BPRHO(NR)/DR/AR1RHO(NR) !/RJCB(NR)
C         BPRHO(NR)= AJ(NR)*DR*AR1RHO(NR)/(FACTOR0*FACTORP)
      DO NR=2,NRMAX-1
         FACTOR0=RR*RR/(AMYU0*DVRHO(NR))
         FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)
         FACTOR2=DVRHO(NR  )*ABRHO(NR  )
         FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)
         FACTORM=0.5D0*(FACTOR1+FACTOR2)
         FACTORP=0.5D0*(FACTOR2+FACTOR3)
         AJ(NR)= FACTOR0*(FACTORP*BPRHO(NR)-FACTORM*BPRHO(NR-1))
     &          /DR/AR1RHO(NR)!/RJCB(NR)
C         BPRHO(NR)=(AJ(NR)*DR*AR1RHO(NR)/FACTOR0+FACTORM*BPRHO(NR-1))
C     &            /FACTORP
      ENDDO
      NR=NRMAX
         FACTOR0=RR*RR/(AMYU0*DVRHO(NR))
         FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)
         FACTOR2=DVRHO(NR  )*ABRHO(NR  )
         FACTORM=0.5D0*(FACTOR1+FACTOR2)
         FACTORP=(3.D0*FACTOR2-FACTOR1)/2.D0
         AJ(NR)= FACTOR0*(FACTORP*BPRHO(NR)-FACTORM*BPRHO(NR-1))
     &          /DR/AR1RHO(NR)!/RJCB(NR)
C         BPRHO(NR)=(AJ(NR)*DR*AR1RHO(NR)/FACTOR0+FACTORM*BPRHO(NR-1))
C     &            /FACTORP
C
c$$$      CALL TRSUMD(AJ  ,DSRHO,NRMAX,AJTSUM)
c$$$      AJT=AJTSUM*DR/1.D6
c$$$C      write(6,*) RIP,AJT,AJT/RIP
c$$$      FACTJ=RIP/AJT
C     
      DO NR=1,NRMAX
         BP(NR)=BPRHO(NR)!*FACTJ
         AJ(NR)=AJ(NR)   !*FACTJ
C     HJRHO(NR)=AJ(NR)*1.D-6
      ENDDO
c$$$      BPSEQ=BP(NRMAX)
c$$$      RIPEQ=RIP
C
c$$$      DO NR=1,NRMAX
c$$$         QL=RKAPS*RA*RG(NR)*BB/(RR*BP(NR))
c$$$         FACTQ(NR)=QRHO(NR)/QL
c$$$      ENDDO
C
      CALL TRGFRG
C
      RETURN
      END
C
C     ***********************************************************
C
C           CONVERGENCE TEST
C
C     ***********************************************************
C
      SUBROUTINE TRCONV(L,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION AJOLD(NRM),AJDLT(NRM)
C
      IERR=0
      L=0
      DO NR=1,NRMAX
         AJOLD(NR)=0.D0
      ENDDO
 100  L=L+1
      IF (L.GT.10) THEN
         WRITE(6,*) 'XX ITERATION IS TOO MUCH! (OVER 10)'
         IERR=1
         RETURN
      ENDIF
      CALL TRSETG
      DO NR=1,NRMAX
         AJDLT(NR)=AJ(NR)-AJOLD(NR)
      ENDDO
      CALL TRSUMJ(AJDLT,RHOTR,NRMAX,SUMJDLT)
      CALL TRSUMJ(AJ   ,RHOTR,NRMAX,SUMJNOW)
      CONV=SQRT((SUMJDLT/DBLE(NRMAX))/(SUMJNOW/DBLE(NRMAX)))
      IF(CONV.GT.1.D-5) THEN
         DO NR=1,NRMAX
            AJOLD(NR)=AJ(NR)
         ENDDO
         GOTO 100
      ENDIF
C
      RETURN
      END
C
C     ********************************************************
C
C           RADIAL INTEGRATION ONLY FOR J CONVERGENCE
C
C     ********************************************************
C
      SUBROUTINE TRSUMJ(A,B,NMAX,SUM)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION A(NMAX),B(NMAX)
C
      SUM=0.D0
      DO 100 N=1,NMAX
         SUM=SUM+A(N)**2*B(N)
  100 CONTINUE
      RETURN
      END
C
C     ***********************************************************
C
C           MODEL SELECTOR
C
C     ***********************************************************
C
      SUBROUTINE TR_EQS_SELECT(INIT)
C
      INCLUDE 'trcomm.inc'
      COMMON /TRINS1/ INS
      SAVE NSSMAX
C
      IF(NSMAX.EQ.1.AND.MDLEOI.EQ.0) MDLEOI=1
      IF(INIT.EQ.0) THEN
         INS=0
      ELSE
         IF(NSSMAX.NE.NSMAX) THEN
            INS=0
         ENDIF
      ENDIF
      NSSMAX=NSMAX
      IF((MDLUF.NE.0.AND.MDNI.NE.0).AND.(NSMAX.EQ.1.OR.NSMAX.EQ.2)) THEN
         IF(NSMAX.EQ.1) INS=1
         NSMAX=3
         PA(3)=12.D0
C         PZ(3)=6.D0
      ENDIF
C
      IF(MDLEQT.EQ.0) THEN
         NEQMAX=MDLEQB+(MDLEQN+MDLEQT+MDLEQU)*NSMAX
     &         +MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ELSEIF(MDLEQT.EQ.1) THEN
         NEQMAX=MDLEQB+NSMAX+(MDLEQT+MDLEQU)*NSMAX
     &         +MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ENDIF
C
      IF(MDLEQN.EQ.0.AND.MDLEQE.EQ.1) THEN
         WRITE(6,*) "ERROR! : MDLEQE can be 1 when MDLEQN is 1."
         STOP
      ENDIF
C
      DO NEQ=1,NEQM
         NSS(NEQ)=-1
         NSV(NEQ)=-1
         NNS(NEQ)=0
         NST(NEQ)=0
      ENDDO
      NEQ=0
      IF(MDLEQB.EQ.1) THEN
         NEQ=NEQ+1
         NSS(NEQ)=0
         NSV(NEQ)=0
      ENDIF
      IND=0
      INDH=0
      INDHD=0
      DO NS=1,NSM
         IF(MDLEQN.EQ.1.OR.MDLEQT.EQ.1) THEN
            IF(MDLEQN.EQ.1) IND=1
            CALL TR_TABLE(NS,NEQ,1,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
         IF(MDLEQT.EQ.1) THEN
            CALL TR_TABLE(NS,NEQ,2,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
         IF(MDLEQU.EQ.1) THEN
            CALL TR_TABLE(NS,NEQ,3,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
      ENDDO
      IF(MDLEQZ.EQ.1) THEN
         IF(NSZMAX.EQ.1) THEN
            NEQ=NEQ+1
            NSS(NEQ  )=5
            NSV(NEQ  )=1
         ELSEIF(NSZMAX.EQ.2) THEN
            NEQ=NEQ+2
            NSS(NEQ-1)=5
            NSV(NEQ-1)=1
            NSS(NEQ  )=6
            NSV(NEQ  )=1
         ENDIF
      ENDIF
      IF(MDLEQ0.EQ.1) THEN
         NEQ=NEQ+2
         NSS(NEQ-1)=7
         NSV(NEQ-1)=1
         NSS(NEQ  )=8
         NSV(NEQ  )=1
      ENDIF
C
      IF(INS.NE.0) THEN
         NEQI=0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.NE.MDLEOI.OR.NSVN.NE.2) THEN
               NEQI=NEQI+1
               NNS(NEQI)=NEQ
            ENDIF
         ENDDO
      ENDIF
C
      NNSC=0
      DO NEQ=1,NEQMAX
         NNSN=NNS(NEQ)
         IF(NNSN.NE.0) THEN
            NNSC=NNSC+1
         ELSE
            GOTO 1000
         ENDIF
      ENDDO
 1000 CONTINUE
      NNSMAX=NNSC
      NEQRMAX=NEQMAX-NNSMAX
C
      NEQS=1
      NEQT=1
      DO NEQ=1,NEQMAX
         NNSN=NNS(NEQ)
         DO NEQ1=1,NEQMAX
            IF(NEQ1.GE.NEQS) THEN
               IF(NEQ1.LT.NNSN.OR.NNSN.EQ.0) THEN
                  NST(NEQ1)=NEQT
                  NEQT=NEQT+1
                  IF(NEQT.GT.NEQRMAX) GOTO 1200
               ELSEIF(NNSN.EQ.NEQ1) THEN
                  NST(NEQ1)=0
                  NEQS=NNSN+1
                  GOTO 1100
               ENDIF
            ENDIF
         ENDDO
 1100    CONTINUE
      ENDDO
 1200 CONTINUE
C
      DO NS=1,NSTMAX
         IF(PZ(NS).NE.0.D0) THEN
            AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
         ELSE
            AMZ(NS)=0.D0
C     I don't know whether this representation is true or not.
         ENDIF
      ENDDO
C
      WRITE(6,600) 'NEQ','NSS','NSV','NNS','NST'
      DO NEQ=1,NEQMAX
         WRITE(6,610) NEQ,NSS(NEQ),NSV(NEQ),NNS(NEQ),NST(NEQ)
      ENDDO
 600  FORMAT(' ',5(' ',A))
 610  FORMAT(' ',5I4)
C
      IF(MDNCLS.EQ.0) THEN
         IF(MDLKAI.EQ.63) THEN
            IF(MDLWLD.EQ.0) THEN
               MDDIAG=0
            ELSE
               MDDIAG=2
            ENDIF
         ELSE
            MDDIAG=0
         ENDIF
      ELSE
         IF(MDLKAI.EQ.63) THEN
            IF(MDLWLD.EQ.0) THEN
               MDDIAG=1
            ELSE
               MDDIAG=3
            ENDIF
         ELSE
            MDDIAG=1
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           SORTER AS MAIN PART OF MODEL SELECTOR
C
C     ***********************************************************
C
      SUBROUTINE TR_TABLE(NS,NEQ,NSW,IND,INDH,INDHD)
C
      INCLUDE 'trcomm.inc'
C
      REM=AME/AMM
      IF(NS.LE.NSMAX) THEN
         IF(ABS(PA(NS)-REM).LE.1.D-10) THEN
            NEQ=NEQ+1
            NSS(NEQ)=1
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1.AND.MDLEQE.EQ.0) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 100
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 100           CONTINUE
            ELSEIF(NSW.EQ.2.AND.IND.EQ.1) THEN
               IF(MDLEQE.EQ.0) THEN
                  IF(NSS(1).EQ.0) THEN
                     NNS(1)=2
                  ELSEIF(NSS(1).EQ.1) THEN
                     NNS(1)=1
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF(PA(NS).EQ.1.D0.OR.PA(NS).EQ.2.D0) THEN
            IF(PA(NS).EQ.1) INDH=1
            NEQ=NEQ+1
            IF(INDH.EQ.1.AND.PA(NS).EQ.2.D0) THEN
               NSS(NEQ)=3
               INDHD=1
            ELSE
               NSS(NEQ)=2
            ENDIF
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 200
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 200           CONTINUE
            ENDIF
         ELSEIF(PA(NS).EQ.3.D0) THEN
            NEQ=NEQ+1
            IF(INDHD.EQ.1) THEN
               NSS(NEQ)=4
            ELSE
               NSS(NEQ)=3
            ENDIF
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 300
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 300           CONTINUE
            ENDIF
         ELSEIF(PA(NS).EQ.4.D0) THEN
            NEQ=NEQ+1
            NSS(NEQ)=4
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 400
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 400           CONTINUE
            ENDIF
         ELSEIF(PA(NS).EQ.12.D0.AND.NSMAX.EQ.3) THEN
            NEQ=NEQ+1
            NSS(NEQ)=3
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 500
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 500           CONTINUE
            ENDIF
         ELSEIF(PA(NS).EQ.0.D0) THEN
            IND=-1
         ENDIF
      ENDIF
C     
      RETURN
      END
C
C     ***********************************************************
C
C           SET GEOMETRIC FACTOR AT HALF MESH
C
C     ***********************************************************
C
      SUBROUTINE TRSTGF
C
      INCLUDE 'trcomm.inc'
C
      IF(MDLUF.NE.0) THEN
         DO NR=1,NRMAX
            EPSRHO(NR)=RA*RG(NR)/RR
C
            TTRHO(NR)=TTRHOU(1,NR)
            DVRHO(NR)=DVRHOU(1,NR)
            DSRHO(NR)=DSRHOU(1,NR)
            ABRHO(NR)=ABRHOU(1,NR)
            ARRHO(NR)=ARRHOU(1,NR)
            AR1RHO(NR)=AR1RHOU(1,NR)
            AR2RHO(NR)=AR2RHOU(1,NR)
            RJCB(NR)=1.D0/(RKAPS*RA)
C            RJCB(NR)=AR1RHOU(1,NR)
            RMJRHO(NR)=RMJRHOU(1,NR)
            RMNRHO(NR)=RMNRHOU(1,NR)
            EKAPPA(NR)=RKAP
         ENDDO
      ELSE
         DO NR=1,NRMAX
            EPSRHO(NR)=RA*RG(NR)/RR
            BPRHO(NR)=BP(NR)
            QRHO(NR)=QP(NR)
            TTRHO(NR)=BB*RR
            DVRHO(NR)=2.D0*PI*RKAP*RA*RA*2.D0*PI*RR*RM(NR)
            DSRHO(NR)=2.D0*PI*RKAP*RA*RA*RM(NR)
            ABRHO(NR)=1.D0/(RKAPS*RA*RR)**2
            ARRHO(NR)=1.D0/RR**2
            AR1RHO(NR)=1.D0/(RKAPS*RA)
            AR2RHO(NR)=1.D0/(RKAPS*RA)**2
            RJCB(NR)=1.D0/(RKAPS*RA)
C
            EKAPPA(NR)=RKAP
            RMJRHO(NR)=RR
            RMNRHO(NR)=RA*RM(NR)
         ENDDO
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           GEOMETRIC FACTOR AT GRID MESH
C
C     ***********************************************************
C
      SUBROUTINE TRGFRG
C
      INCLUDE 'trcomm.inc'
C
      DO NR=1,NRMAX-1
         AR1RHOG(NR)=0.5D0*(AR1RHO(NR)+AR1RHO(NR+1))
         AR2RHOG(NR)=0.5D0*(AR2RHO(NR)+AR2RHO(NR+1))
         RMJRHOG(NR)=0.5D0*(RMJRHO(NR)+RMJRHO(NR+1))
         RMNRHOG(NR)=0.5D0*(RMNRHO(NR)+RMNRHO(NR+1))
         TTRHOG (NR)=0.5D0*(TTRHO (NR)+TTRHO (NR+1))
         DVRHOG (NR)=0.5D0*(DVRHO (NR)+DVRHO (NR+1))
         ARRHOG (NR)=0.5D0*(ARRHO (NR)+ARRHO (NR+1))
         ABRHOG (NR)=0.5D0*(ABRHO (NR)+ABRHO (NR+1))
      ENDDO
      NR=NRMAX
         RGL=RG(NR)
         RML=RM(NR)
         RML1=RM(NR-1)
         AR1RHOG(NR)=FEDG(RGL,RML,RML1,AR1RHO(NR),AR1RHO(NR-1))
         AR2RHOG(NR)=FEDG(RGL,RML,RML1,AR2RHO(NR),AR2RHO(NR-1))
         RMJRHOG(NR)=FEDG(RGL,RML,RML1,RMJRHO(NR),RMJRHO(NR-1))
         RMNRHOG(NR)=FEDG(RGL,RML,RML1,RMNRHO(NR),RMNRHO(NR-1))
         TTRHOG (NR)=FEDG(RGL,RML,RML1,TTRHO (NR),TTRHO (NR-1))
         DVRHOG (NR)=FEDG(RGL,RML,RML1,DVRHO (NR),DVRHO (NR-1))
         ARRHOG (NR)=FEDG(RGL,RML,RML1,ARRHO (NR),ARRHO (NR-1))
         ABRHOG (NR)=FEDG(RGL,RML,RML1,ABRHO (NR),ABRHO (NR-1))
C
      RETURN
      END
C
C     ***********************************************************
C
C           EDGE VALUE SELECTOR
C
C     ***********************************************************
C
      SUBROUTINE TR_EDGE_SELECTOR(NSW)
C
      INCLUDE 'trcomm.inc'
      DIMENSION PNSSO(NSTM),PTSO(NSTM),PNSSAO(NSTM),PTSAO(NSTM)
      SAVE PNSSO,PTSO,PNSSAO,PTSAO
C
      IF(RHOA.EQ.1.D0) RETURN
C
      IF(MDLUF.EQ.0) THEN
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS(NS)
               PTSO (NS)=PTS (NS)
C
               PNSS (NS)=PNSSAO(NS)
               PTS  (NS)=PTSAO (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSS (NS)=PNSSO(NS)
               PTS  (NS)=PTSO (NS)
            ENDDO
         ENDIF
      ELSE
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS (NS)
               PTSO (NS)=PTS  (NS)
C
               PNSS (NS)=PNSSA(NS)
               PTS  (NS)=PTSA (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSS (NS)=PNSSO(NS)
               PTS  (NS)=PTSO (NS)
            ENDDO
         ENDIF
      ENDIF
      GOTO 9000
C
      ENTRY TR_EDGE_DETERMINER(NSW)
C
      IF(MDLUF.EQ.0.AND.RHOA.NE.1.D0) THEN
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSAO(NS)=PNSS(NS)
               PTSAO (NS)=PTS (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSSAO(NS)=RN(NRAMAX,NS)
               PTSAO (NS)=RT(NRAMAX,NS)
               PNSSA (NS)=PNSSAO(NS)
               PTSA  (NS)=PTSAO (NS)
            ENDDO
         ENDIF
      ELSE
         RETURN
      ENDIF
C
 9000 RETURN
      END

