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
      INCLUDE 'trcomm.h'
C
      PI      = ASIN(1.D0)*2.D0
      AME     = 9.1093897D-31
      AMM     = 1.6605402D-27
      AEE     = 1.60217733D-19
      VC      = 2.99792458D8
      AMYU0   = 4.D0*PI*1.D-7
      AEPS0   = 1.D0/(VC*VC*AMYU0)
      RKEV    = AEE*1.D3
C
      RR      = 3.0D0     
      RA      = 1.2D0      
      RKAP    = 1.5D0
      RDLT    = 0.0D0
      BB      = 3.D0        
      RIPS    = 3.D0         
      RIPE    = 3.D0         
      RIPSS   = 3.D0
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
      PROFJ1 = -2.D0
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
      MDLKAI = 31
      MDLETA = 1
C      MDLAD  = 1
      MDLAD  = 0
      MDLAVK = 0
      MDLJBS = 3
      MDLKNC = 1
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
      EPSLTR = 0.001D0
C      EPSLTR = 1.D99
      LMAXTR = 10
C
      CHP    = 0.D0
      CK0    = 12.D0
      CKALFA = 0.D0
      CKBETA = 0.D0
      CKGUMA = 0.D0
C
      TPRST  = 0.1D0
      MDLST  = 0
C
      MDLNF  = 0
      IZERO  = 3
C
      PNBTOT = 0.D0
      PNBR0  = 0.D0
      PNBRW  = 0.5D0
      PNBENG = 80.D0
      PNBRTG = 3.D0
      MDLNB  = 1
C
      PECTOT = 0.D0
      PECR0  = 0.D0
      PECRW  = 0.2D0
      PECTOE = 1.D0
      PECNPR = 0.D0
      MDLEC  = 0
C
      PLHTOT = 0.D0
      PLHR0  = 0.D0
      PLHRW  = 0.2D0
      PLHTOE = 1.D0
      PLHNPR = 2.D0
      MDLLH  = 0
C
      PICTOT = 0.D0
      PICR0  = 0.D0
      PICRW  = 0.5D0
      PICTOE = 0.5D0
      PICNPR = 2.D0
      MDLIC  = 0
C
      PNBCD  = 1.D0
      PECCD  = 0.D0
      PLHCD  = 1.D0
      PICCD  = 0.D0
      PBSCD  = 1.D0
      MDLCD  = 0
C
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
      MODELG=0
      NTEQIT=0
      MDLUF=0
C
      MDALL=0
      MDQ=0
      MDCUR=0
      MDTT=0
      MDVOL=0
      MDAREA=0
      MDRMJ=0
      MDGR1=0
      MDGR2=0
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
      SUBROUTINE TRPARM(KID)
C
      INCLUDE 'trcomm.h'
C
      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RIPSS,
     &              PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              PROFJ1,PROFJ2,ALP,AD0,AV0,CNC,CDW,
     &              MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,
     &              DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST,
     &              EPSLTR,LMAXTR,CHP,CK0,CKALFA,CKBETA,CKGUMA,TPRST,
     &              MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDLUF,
     &              PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB,
     &              PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC,
     &              PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH,
     &              PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC,
     &              PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD,
     &              PELTOT,PELR0,PELRW,PELRAD,PELVEL,MDLPEL,
     &              PELTIM,PELPAT,KFNLOG
C
      LOGICAL LEX
      CHARACTER KPNAME*32,KLINE*70,KNAME*80,KID*1
C
      MODE=0
    1    CONTINUE
         WRITE(6,*) '# INPUT &TR :'
         READ(5,TR,ERR=2,END=3)
         KID=' '
         GOTO 4
C
    2    CALL TRPLST
      GOTO 1
C
    3 KID='Q'
    4 GOTO 3000
C
      ENTRY TRPARL(KLINE)
C
      MODE=1
      KNAME=' &TR '//KLINE//' &END'
      WRITE(33,'(A80)') KNAME
      REWIND(33)
      READ(33,TR,ERR=8,END=8)
      WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
      GOTO 9
    8 CALL TRPLST
    9 REWIND(33)
      GOTO 3000
C
      ENTRY TRPARF
C
      MODE=2
      KPNAME='trparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(.NOT.LEX) RETURN
C
      OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
      READ(25,TR,ERR=9800,END=9900)
      WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
C
 3000 IERR=0
C
      IF(NRMAX.LT.1.OR.NRMAX.GT.NRM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX,NRM =',NRMAX,NRM
         NRMAX=NRM
         IERR=1
      ENDIF
C
      IF(NTMAX.LT.1.OR.NTMAX.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM
         IERR=1
      ENDIF
C
      IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1
C
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      RETURN
C
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE TRPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &TR : RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RIPSS'/
     &       ' ',8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/
     &       ' ',8X,'PNC,PNFE,PNNU,PNNUS'/
     &       ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/
     &       ' ',8X,'PROFJ1,PROFJ2,ALP'/
     &       ' ',8X,'CK0,CNC,CDW,CKALFA,CKBETA,KFNLOG,MDLKNC'/
     &       ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/
     &       ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/
     &       ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/
     &       ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,PNBCD,MDLNB'/
     &       ' ',8X,'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/
     &       ' ',8X,'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/
     &       ' ',8X,'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/
     &       ' ',8X,'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/
     &       ' ',8X,'PELTIM,PELPAT,MODELG,NTEQIT,MDLUF')
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
      INCLUDE 'trcomm.h'
C
      WRITE(6,*) '** TRANSPORT **'
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
      DO NS=1,NSM
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
      WRITE(6,602) 'MDLJBS',MDLJBS,
     &             'MDLKNC',MDLKNC
C
      WRITE(6,601) 'CK0   ',CK0,
     &             'CNC   ',CNC,
     &             'AD0   ',AD0,
     &             'CHP   ',CHP
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
     &             'MDLNF ',MDLNF,
     &             'MODELG',MODELG
C
      WRITE(6,602) 'NTEQIT',NTEQIT,
     &             'MDLUF ',MDLUF
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
      RETURN
C
  601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  603 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
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
      INCLUDE 'trcomm.h'
      PARAMETER(NURM=NRM+1)
      DIMENSION RAD(NURM),FQ(NURM),FNE(NURM),FTE(NURM),FTI(NURM)
      DIMENSION FDVOL(NURM),FPBE(NURM),FPBI(NURM),FCUR(NURM),DERIV(NURM)
      DIMENSION UQPH(4,NURM),UNEPH(4,NURM),UTEPH(4,NURM),UTIPH(4,NURM)
      DIMENSION UDVOLPH(4,NURM),UPBEPH(4,NURM),UPBIPH(4,NURM)
      DIMENSION UCURPH(4,NURM)
      CHARACTER KFILE*10
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
C
      MODEP=2
      IF(MDLUF.EQ.2.AND.(MODEP.EQ.1.OR.MODEP.EQ.2)) THEN
         WRITE(6,*) "*****"
         KFILE='NE.PHI'
         CALL UFREAD(KFILE,RAD,FNE,NUFMAX,MDNE,IERR)
         CALL SPL1D(RAD,FNE,DERIV,UNEPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FNE: IERR=',IERR
         RADNOW=RAD(NRMAX+1)
         CALL SPL1DF(RADNOW,PNEPH,RAD,UNEPH,NUFMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF FNE: IERR=',IERR
         PNS(1)=PNEPH*1.D-20
         PNS(2)=PNEPH*1.D-20-2.D-8
         PNS(3)=1.D-8
         PNS(4)=1.D-8
C
         KFILE='TE.PHI'
         CALL UFREAD(KFILE,RAD,FTE,NUFMAX,MDTE,IERR)
         CALL SPL1D(RAD,FTE,DERIV,UTEPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FTE: IERR=',IERR
         RADNOW=RAD(NRMAX+1)
         CALL SPL1DF(RADNOW,PTEPH,RAD,UTEPH,NUFMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF FTE: IERR=',IERR
         PTS(1)=PTEPH*1.D-3
C     
         KFILE='TI.PHI'
         CALL UFREAD(KFILE,RAD,FTI,NUFMAX,MDTE,IERR)
         CALL SPL1D(RAD,FTI,DERIV,UTIPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FTI: IERR=',IERR
         RADNOW=RAD(NRMAX+1)
         CALL SPL1DF(RADNOW,PTIPH,RAD,UTIPH,NUFMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF FTI: IERR=',IERR
         PTS(2)=PTIPH*1.D-3
         PTS(3)=PTIPH*1.D-3
         PTS(4)=PTIPH*1.D-3
C
         KFILE='DVOL.PHI'
         CALL UFREAD(KFILE,RAD,FDVOL,NUFMAX,MDDVOL,IERR)
         CALL SPL1D(RAD,FDVOL,DERIV,UDVOLPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FDVOL: IERR=',IERR
C
         KFILE='PBE.PHI'
         CALL UFREAD(KFILE,RAD,FPBE,NUFMAX,MDPBE,IERR)
         CALL SPL1D(RAD,FPBE,DERIV,UPBEPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FPBE: IERR=',IERR
C
         KFILE='PBI.PHI'
         CALL UFREAD(KFILE,RAD,FPBI,NUFMAX,MDPBI,IERR)
         CALL SPL1D(RAD,FPBI,DERIV,UPBIPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FPBI: IERR=',IERR
      ENDIF
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
      DR = 1.D0/DBLE(NRMAX)
      RKAPX=(RKAP-1.D0)/(RKAP+1.D0)
      FKAP=0.5D0*(RKAP+1.D0)
     &     *(1.D0+RKAPX/4.D0+RKAPX*RKAPX/64.D0)
C      FKAP=RKAP
C
      DO NR=1,NRMAX
         RG(NR)  = DBLE(NR*DR)
         RM(NR)  = DBLE(NR-0.5D0)*DR
C
         IF(MDLUF.NE.2.OR.MODEP.NE.2) THEN
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
         ELSE
            RMNOW=RM(NR)
            CALL SPL1DF(RMNOW,PNE,RAD,UNEPH,NUFMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF PNE: IERR=',IERR
            RN(NR,1) = PNE*1.D-20
            RN(NR,2) = PNE*1.D-20-2.D-8
            RN(NR,3) = 1.D-8
            RN(NR,4) = 1.D-8
C
            RMNOW=RM(NR)
            CALL SPL1DF(RMNOW,PDVOL,RAD,UDVOLPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRPROF: SPL1DF PDVOL: IERR=',IERR
C
            RMNOW=RM(NR)
            CALL SPL1DF(RMNOW,PTE,RAD,UTEPH,NUFMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF PTE: IERR=',IERR
            RT(NR,1) = PTE*1.D-3
            CALL SPL1DF(RMNOW,PTI,RAD,UTIPH,NUFMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF PTI: IERR=',IERR
            DO NS=2,NSM
               RT(NR,NS) = PTI*1.D-3
            ENDDO
C
            RMNOW=RM(NR)
            CALL SPL1DF(RMNOW,PBE,RAD,UPBEPH,NUFMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF PBE: IERR=',IERR
            PEX(NR,1)=PBE
            CALL SPL1DF(RMNOW,PBI,RAD,UPBIPH,NUFMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF PBI: IERR=',IERR
            PEX(NR,2)=PBI
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0
         ENDIF
         NS=1
C
         PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1)**PROFU2
         ANNU(NR)= (PNNU-PNNUS)*PROF+PNNUS
C
         DO NF=1,NFM
           RW(NR,NF) = 0.D0
         ENDDO
      ENDDO
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
      DO NR=1,NRMAX
         ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC
     &            *1.D-2*RN(NR,1)
         ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE
     &            *1.D-2*RN(NR,1)
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
C
C     *** CALCULATE PROFILE OF AJ(R) ***
C
C     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***
C
      IF(MDLUF.NE.2) THEN
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
      BP(1)=(RM(1)*RA*DR*AMYU0*AJ(1))/RG(1)
      DO NR=2,NRMAX
         BP(NR)=(RG(NR-1)*BP(NR-1)+RM(NR)*RA*DR*AMYU0*AJ(NR))/RG(NR)
      ENDDO
C
      BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*FKAP)
      FACT=BPS/BP(NRMAX)
      DO NR=1,NRMAX
         AJOH(NR)=FACT*AJOH(NR)
         AJ(NR)  =AJOH(NR)
         BP(NR)  =FACT*BP(NR)
         QP(NR)  =FKAP*RA*RG(NR)*BB/(RR*BP(NR))
c$$$         write(6,*) QP(NR),BP(NR),AJ(NR)
      ENDDO
      ELSE
         KFILE='Q.PHI'
         CALL UFREAD(KFILE,RAD,FQ,NUFMAX,MDQ,IERR)
         CALL SPL1D(RAD,FQ,DERIV,UQPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FQ: IERR=',IERR
         KFILE='CUR.PHI'
         CALL UFREAD(KFILE,RAD,FCUR,NUFMAX,MDCUR,IERR)
         CALL SPL1D(RAD,FCUR,DERIV,UCURPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1D FCUR: IERR=',IERR
         DO NR=1,NRMAX
            RGNOW=RG(NR)
            CALL SPL1DF(RGNOW,PQ,RAD,UQPH,NUFMAX,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRPROF: SPL1DF PQ: IERR=',IERR
            CALL SPL1DF(RGNOW,PCUR,RAD,UCURPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRPROF: SPL1DF PCUR: IERR=',IERR
            QP(NR)=PQ
            BP(NR)=FKAP*RA*RG(NR)*BB/(RR*QP(NR))
         ENDDO
c$$$         AJ(1)=BP(1)*RG(1)/(RM(1)*RA*DR*AMYU0)
c$$$         AJOH(1)=AJ(1)
c$$$         DO NR=2,NRMAX
c$$$            AJ(NR)=(RG(NR)*BP(NR)-RG(NR-1)*BP(NR-1))
c$$$     &            /(RM(NR)*RA*DR*AMYU0)
c$$$            AJOH(NR)=AJ(NR)
c$$$         ENDDO
         DO NR=1,NRMAX
            RGNOW=RG(NR)
            CALL SPL1DF(RGNOW,PCUR,RAD,UCURPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRPROF: SPL1DF PCUR: IERR=',IERR
            AJ(NR)=PCUR
            AJOH(NR)=PCUR
         ENDDO
         RIP=2.D0*PI*RA*FKAP*BP(NRMAX)/AMYU0*1.D-6
         RIPS=RIP
         RIPSS=RIP
         RIPE=RIP
c$$$         DO NR=1,NRMAX
c$$$            write(6,*) NR,AJ(NR)
c$$$         ENDDO
      ENDIF
C      Q0=(4.D0*QP(1)-QP(2))/3.D0
C      STOP
C
C     *** THIS MODEL ASSUMES CONSTANT EZ ***
C
      IF(PROFJ1.LE.0.D0) THEN
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
               VTE=SQRT(ABS(TE)*RKEV/AME)
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
            ENDIF 
            AJOH(NR)=1.D0/ETA(NR)
            AJ(NR)  =1.D0/ETA(NR)
         ENDDO
C
         BP(1)=(RM(1)*RA*DR*AMYU0*AJ(1))/RG(1)
         DO NR=2,NRMAX
            BP(NR)=(RG(NR-1)*BP(NR-1)+RM(NR)*RA*DR*AMYU0*AJ(NR))/RG(NR)
         ENDDO
C
         BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*FKAP)
         FACT=BPS/BP(NRMAX)
         DO NR=1,NRMAX
            AJOH(NR)=FACT*AJOH(NR)
            AJ(NR)  =AJOH(NR)
            BP(NR)  =FACT*BP(NR)
            QP(NR)  =FKAP*RA*RG(NR)*BB/(RR*BP(NR))
         ENDDO
         Q0=(4.D0*QP(1)-QP(2))/3.D0
      ENDIF
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
            QP(NR)=-QRHO(NR)
C            write(6,*) NR,QP(NR)
         ENDDO
C         PAUSE
C         STOP
         Q0=(4.D0*QP(1)-QP(2))/3.D0
      ENDIF
C
      DRIP  = (RIPSS-RIPS)/1.D1
 1000 RIP   = RIPSS
      IF(MODELG.EQ.3) THEN
         CALL TRCONV(L,IERR)
         write(6,*) "L=",L
      ELSE
         CALL TRSETG
      ENDIF
C
      IF(MODELG.EQ.0) THEN
         GRG(1)=0.0
         DO NR=1,NRMAX
            GRM(NR)  =GCLIP(RA*RM(NR))
            GRG(NR+1)=GCLIP(RA*RG(NR))
         ENDDO
      ELSE
         GRG(1)=0.0
         DO NR=1,NRMAX
            GRM(NR)  =GCLIP(RM(NR))
            GRG(NR+1)=GCLIP(RG(NR))
         ENDDO
C
         RIPSS=RIPSS-DRIP
C         write(6,'(A,1P4E12.5)') "RIP,RIPSS,RIPS,RIPE= ",RIP,RIPSS,RIPS
C     &        ,RIPE
         IF(DRIP.NE.0.AND.ABS(RIPSS-RIPS).LT.1.D-10) THEN
C            write(6,*) RIPSS-RIPS
            GOTO 1000
         ENDIF
C         RIP=RIPSS
C         write(6,*) RIP
C         STOP
      ENDIF
      RETURN
      END
C
C     ***********************************************************
C
C           SET GEOMETRICAL FACTOR
C
C     ***********************************************************
C
      SUBROUTINE TRSETG
C
      INCLUDE 'trcomm.h'
      PARAMETER(NURM=NRM+1)
      DIMENSION DERIV(NURM+2)
      DIMENSION ASR(NURM)
      DIMENSION FFQ(NURM),FFTT(NURM),FFVOL(NURM),FFDVOL(NURM)
      DIMENSION FFAREA(NURM),FFRMJ(NURM),FFGR1(NURM),FFGR2(NURM)
      DIMENSION UQEQ(4,NURM),UTTEQ(4,NURM)
      DIMENSION UVOLEQ(4,NURM),UDVOLEQ(4,NURM),UAREAEQ(4,NURM)
      DIMENSION URMJEQ(4,NURM),UGR1EQ(4,NURM),UGR2EQ(4,NURM)
      CHARACTER KFILE*10
C
      IF(MODELG.EQ.3) THEN
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
            HJRHO(NR)=AJ(NR)*1.D-6
            VTRHO(NR)=0.D0
C            WRITE(6,'(A,I5,1P4E12.4)')
C     &           'NR,P/HJ/T/VTRHO=',NR,PRHO(NR),
C     &           HJRHO(NR),TRHO(NR),VTRHO(NR)
         ENDDO
C         PAUSE
C
         DO NR=1,NRMAX
C            WRITE(6,'(A,2I5,1P4E12.4)')
C     &           'NR,I/RM/J/V/T=',NR,NRMAX,RIP,
C     &           HJRHO(NR),VTRHO(NR),TRHO(NR)
         ENDDO
C
         CALL TREQEX(RIP,NRMAX,PRHO,HJRHO,VTRHO,TRHO,
     &               QRHO,TTRHO,DVRHO,DSRHO,
     &               ABRHO,ARRHO,AR1RHO,AR2RHO,
     &               EPSRHO,MDLUF,IERR)
C
         DO NR=1,NRMAX
C            WRITE(6,'(A,I5,1P4E12.4)')
C     &           'NR,Q/TT/DV/AB=',NR,QRHO(NR),
C     &           TTRHO(NR),DVRHO(NR),ABRHO(NR)
C            WRITE(6,'(A,I5,1P4E12.4)')
C     &           'NR,Q/TT/AB/EP=',NR,QRHO(NR),
C     &           TTRHO(NR),ABRHO(NR),EPSRHO(NR)
         ENDDO
C         PAUSE
C
C         DO NR=1,NRMAX
C            BPNR=RA*RHOTR(NR)*TTRHO(NR)*ARRHO(NR)/QRHO(NR)
C            write(6,*) NR,BPNR,BP(NR)
C            BP(NR)=BPNR
C         ENDDO
         NR=1
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
            FACTORP=0.5D0*(FACTOR2+FACTOR3)
            AJ(NR)= FACTOR0*FACTORP*BP(NR)/DR
            BPRHO(NR)= AJ(NR)*DR/(FACTOR0*FACTORP)
         DO NR=2,NRMAX-1
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
            FACTORM=0.5D0*(FACTOR1+FACTOR2)
            FACTORP=0.5D0*(FACTOR2+FACTOR3)
            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR
            BPRHO(NR)=(AJ(NR)*DR/FACTOR0+FACTORM*BPRHO(NR-1))/FACTORP
         ENDDO
         NR=NRMAX
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTORM=0.5D0*(FACTOR1+FACTOR2)
            FACTORP=(3.D0*FACTOR2-FACTOR1)/2.D0
            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR
            BPRHO(NR)=(AJ(NR)*DR/FACTOR0+FACTORM*BPRHO(NR-1))/FACTORP
C            write(6,*) "AJ(NR)=",AJ(NR)
c$$$         DO NR=2,NRMAX-1
c$$$C            write(6,*) TTRHO(NR),ARRHO(NR),DVRHO(NR)
c$$$C            write(6,*) ABRHO(NR)
c$$$            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
c$$$            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
c$$$            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
c$$$            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
c$$$            FACTORM=0.5D0*(FACTOR1+FACTOR2)
c$$$            FACTORP=0.5D0*(FACTOR2+FACTOR3)
c$$$            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR
c$$$C           write(6,'(A,1P5E12.4)') 'F0,FP,FM,-,AJ='
c$$$C     &           ,FACTOR0,FACTORP*BP(NR),FACTORM*BP(NR-1)
c$$$C     &           ,FACTORP*BP(NR)-FACTORM*BP(NR-1),AJ(NR)
c$$$         ENDDO
C
         DO NR=1,NRMAX
            BP(NR)=BPRHO(NR)
            QP(NR)=QRHO(NR)
            HJRHO(NR)=AJ(NR)*1.D-6
C            write(6,*) NR,HJRHO(NR)
         ENDDO
      ELSE
         IF(MDLUF.NE.2) THEN
            DO NR=1,NRMAX
               BPRHO(NR)=BP(NR)
               QRHO(NR)=QP(NR)
               TTRHO(NR)=BB*RR
               DVRHO(NR)=2.D0*PI*RKAP*RA*RA*2.D0*PI*RR*RM(NR)
               DSRHO(NR)=2.D0*PI*FKAP*RA*RA*RM(NR)
               ABRHO(NR)=1.D0/(RA*RR)**2
               ARRHO(NR)=1.D0/RR**2
               AR1RHO(NR)=1.D0/RA
               AR2RHO(NR)=1.D0/RA**2
C               EPSRHO(NR)=RA*RM(NR)/RR
               EPSRHO(NR)=RA*RG(NR)/RR
            ENDDO
         ELSE
C
         KFILE='VOL.PHI'
         CALL UFREAD(KFILE,ASR,FFVOL,NUFMAX,MDVOL,IERR)
         CALL SPL1D(ASR,FFVOL,DERIV,UVOLEQ,NUFMAX,0,IERR)
C         DO NUF=1,NUFMAX
C            RHOFDV=ASR(NUF)
C            CALL SPL1DD(RHOFDV,VPL,DVPL,ASR,UVOLEQ,NUFMAX,IERR)
C            FFDVOL(NUF)=DVPL
C         ENDDO
C         CALL SPL1D(ASR,FFDVOL,DERIV,UDVOLEQ,NUFMAX,0,IERR)
C
         KFILE='AREA.PHI'
         CALL UFREAD(KFILE,ASR,FFAREA,NUFMAX,MDAREA,IERR)
         CALL SPL1D(ASR,FFAREA,DERIV,UAREAEQ,NUFMAX,0,IERR)
C
         KFILE='RMJ.PHI'
         CALL UFREAD(KFILE,ASR,FFRMJ,NUFMAX,MDRMJ,IERR)
         CALL SPL1D(ASR,FFRMJ,DERIV,URMJEQ,NUFMAX,0,IERR)
C
         KFILE='GR1.PHI'
         CALL UFREAD(KFILE,ASR,FFGR1,NUFMAX,MDGR1,IERR)
         CALL SPL1D(ASR,FFGR1,DERIV,UGR1EQ,NUFMAX,0,IERR)
C     
         KFILE='GR2.PHI'
         CALL UFREAD(KFILE,ASR,FFGR2,NUFMAX,MDGR2,IERR)
         CALL SPL1D(ASR,FFGR2,DERIV,UGR2EQ,NUFMAX,0,IERR)
C
c$$$         KFILE='TT.PHI'
c$$$         CALL UFREAD(KFILE,ASR,FFTT,NUFMAX,MDTT,IERR)
c$$$         CALL SPL1D(ASR,FFTT,DERIV,UTTEQ,NUFMAX,0,IERR)
         MDTT=1
         CALL UFTTRHO(FFTT)
         CALL SPL1D(ASR,FFTT,DERIV,UTTEQ,NUFMAX,0,IERR)
C
C     *****
C
         DO NR=1,NRMAX
            BPRHO(NR)=BP(NR)
            QRHO(NR)=QP(NR)
            EPSRHO(NR)=RA*RG(NR)/RR
C
            RADNOW=RM(NR)
C
            CALL SPL1DF(RADNOW,TTL,ASR,UTTEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF TTL2: IERR=',IERR
            TTRHO(NR)=TTL
C
            CALL SPL1DD(RADNOW,VPL,DVPL,ASR,UVOLEQ,NUFMAX,IERR)
C            CALL SPL1DF(RADNOW,VPL,ASR,UDVOLEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF VPL2: IERR=',IERR
            DVRHO(NR)=DVPL
C
            CALL SPL1DD(RADNOW,SPL,DSPL,ASR,UAREAEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF SPL2: IERR=',IERR
            DSRHO(NR)=DSPL
C
            CALL SPL1DF(RADNOW,AVRRL,ASR,URMJEQ,NUFMAX,IERR)
            CALL SPL1DF(RADNOW,AVR2L,ASR,UGR2EQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF AVBRL2: IERR=',IERR
            ABRHO(NR)=AVRRL*AVR2L
C
            CALL SPL1DF(RADNOW,AVRRL,ASR,URMJEQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF AVRRL2: IERR=',IERR
            ARRHO(NR)=1.D0/(AVRRL**2)
C
            CALL SPL1DF(RADNOW,AVR1L,ASR,UGR1EQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF AVR1L2: IERR=',IERR
            AR1RHO(NR)=AVR1L
C
            CALL SPL1DF(RADNOW,AVR2L,ASR,UGR2EQ,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRINIT: SPL1DF AVR2L2: IERR=',IERR
            AR2RHO(NR)=AVR2L
            ENDDO
         ENDIF
      ENDIF
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
      INCLUDE 'trcomm.h'
      DIMENSION AJOLD(NRM),AJDLT(NRM)
C
      IERR=0
      L=0
      DO NR=1,NRMAX
         AJOLD(NR)=0.D0
      ENDDO
 200  L=L+1
      IF (L.GT.50) THEN
         WRITE(6,*) 'XX ITERATION IS TOO MUCH! (OVER 50)'
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
         GOTO 200
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
