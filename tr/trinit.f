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
      AME     = 9.10953D-31
      AMM     = 1.67290D-27
      AEE     = 1.60219D-19
      VC      = 2.99792458D8
      AMYU0   = 4.D0*PI*1.D-7
      AEPS0   = 1.D0/(VC*VC*AMYU0)
      RKEV    = AEE*1.D3
C
      RR      = 3.0D0     
      RA      = 1.2D0      
      RKAP    = 1.5D0
      BB      = 3.D0        
      RIPS    = 3.D0         
      RIPE    = 3.D0         
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
      NAMELIST /TR/ RR,RA,RKAP,BB,RIPS,RIPE,
     &              PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              PROFJ1,PROFJ2,ALP,AD0,AV0,CNC,CDW,
     &              MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,
     &              DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST,
     &              EPSLTR,LMAXTR,CHP,CK0,CKALFA,CKBETA,CKGUMA,TPRST,
     &              MDLST,MDLNF,IZERO,
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
  601 FORMAT(1H ,'# &TR : RR,RA,RKAP,BB,RIPS,RIPE'/
     &       1H ,8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/
     &       1H ,8X,'PNC,PNFE,PNNU,PNNUS'/
     &       1H ,8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/
     &       1H ,8X,'PROFJ1,PROFJ2,ALP'/
     &       1H ,8X,'CK0,CNC,CDW,CKALFA,CKBETA,KFNLOG,MDLKNC'/
     &       1H ,8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/
     &       1H ,8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/
     &       1H ,8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/
     &       1H ,8X,'PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,PNBCD,MDLNB'/
     &       1H ,8X,'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/
     &       1H ,8X,'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/
     &       1H ,8X,'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/
     &       1H ,8X,'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/
     &       1H ,8X,'PELTIM,PELPAT')
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
     &             'BB    ',BB
      WRITE(6,601) 'RIPS  ',RIPS,
     &             'RIPE  ',RIPE
C
      WRITE(6,611)
  611 FORMAT(1H ,'NS',2X,'PA           PZ    PN(E20)  PNS(E20) ',
     &                   'PT(KEV)  PTS(KEV)  PELPAT')
      DO NS=1,NSM
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS),
     &                PELPAT(NS)
  612    FORMAT(1H ,I2,1PD12.4,0PF6.1,0P5F9.4)
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
  613 FORMAT(1H ,'    AKDW(E) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       1H ,'    AKDW(D) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       1H ,'    AKDW(T) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       1H ,'    AKDW(A) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW')
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
      RETURN
C
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  603 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',1PE11.3:
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
C
C     ZEFF=1
C
      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
      T     = 0.D0
      TPRE  = 0.D0
      TST   = 0.D0
      VSEC  = 0.D0
      NGR   = 0
      NGT   = 0
      NGST  = 0
C
      DR = RA/DBLE(NRMAX)
      RKAPX=(RKAP-1.D0)/(RKAP+1.D0)
      FKAP=0.5D0*(RKAP+1.D0)
     &     *(1.D0+RKAPX/4.D0+RKAPX*RKAPX/64.D0)
C
      DO 100 NR=1,NRMAX
         RG(NR)  = NR*DR
         RM(NR)  = (NR-0.5D0)*DR
C
         PROF   = (1.D0-(ALP(1)*RM(NR)/RA)**PROFN1)**PROFN2
         DO 10 NS=1,NSM
            RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
   10    CONTINUE
C
         PROF   = (1.D0-(ALP(1)*RM(NR)/RA)**PROFT1)**PROFT2
         DO 20 NS=1,NSM
            RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
   20    CONTINUE
C
         PROF   = (1.D0-(ALP(1)*RM(NR)/RA)**PROFU1)**PROFU2
         ANNU(NR)= (PNNU-PNNUS)*PROF+PNNUS
C
         DO 30 NF=1,NFM
           RW(NR,NF) = 0.D0
   30    CONTINUE
  100 CONTINUE
C
C     *** CALCULATE PZC,PZFE ***
C
      CALL TRZEFF
C
C     *** CALCULATE ANEAVE ***
C
      ANESUM=0.D0
      DO 200 NR=1,NRMAX
         ANESUM=ANESUM+RN(NR,1)*RM(NR)
  200 CONTINUE
      ANEAVE=ANESUM*2.D0*PI*DR/(PI*RA*RA)
C
C     *** CALCULATE IMPURITY DENSITY
C                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***
C
      DO 300 NR=1,NRMAX
         ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC
     &            *1.D-2*RN(NR,1)
         ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE
     &            *1.D-2*RN(NR,1)
         ANI = 0.D0
         DO 210 NS=2,NSM
            ANI=ANI+PZ(NS)*RN(NR,NS)
  210    CONTINUE
         ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DILUTE = 1.D0-ANZ/ANI
         DO 220 NS=2,NSM
            RN(NR,NS) = RN(NR,NS)*DILUTE
  220    CONTINUE
  300 CONTINUE
      PNSS(1)=PNS(1)
      DO 310 NS=2,NSM
         PNSS(NS)=PNS(NS)*DILUTE
  310 CONTINUE
C
C     *** CALCULATE Q(R) ***
C
C
C     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***
C
      DO 600 NR=1,NRMAX
         IF((1.D0-(RM(NR)/RA)**ABS(PROFJ1)).LE.0.D0) THEN
            PROF=0.D0    
         ELSE             
            PROF= (1.D0-(RM(NR)/RA)**PROFJ1)**ABS(PROFJ2)
         ENDIF             
         AJOH(NR)= PROF
  600 CONTINUE
C
         BP(1)=(RM(1)*DR*AJOH(1))/RG(1)
         DO 700 NR=2,NRMAX
            BP(NR)=(RG(NR-1)*BP(NR-1)+RM(NR)*DR*AJOH(NR))/RG(NR)
  700    CONTINUE
C
      BPS= AMYU0*RIPS*1.D6/(2.D0*PI*RA*FKAP)
      FACT=BPS/BP(NRMAX)
C
      DO 800 NR=1,NRMAX
         BP(NR)=FACT*BP(NR)
         QP(NR)=FKAP*RG(NR)*BB/(RR*BP(NR))
  800 CONTINUE
C
      IF(PROFJ1.LE.0.D0) THEN
C
C     *** THIS MODEL ASSUMES CONSTANT EZ ***
C
         CALL TRZEFF
         DO 810 NR=1,NRMAX
C
C        ****** CLASSICAL RESISTIVITY ******
C
         ANE=RN(NR,1)
         TEL =ABS(RT(NR,1))
         ZEFFL=ZEFF(NR)
C
         COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &         /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(TEL*RKEV)**1.5D0/SQRT(2.D0)
C
         ETA(NR) = AME/(ANE*1.D20*AEE*AEE*TAUE)
     &             *(0.29D0+0.46D0/(1.08D0+ZEFFL))
C
C        ****** NEOCLASSICAL RESISTIVITY ******
C
         IF(MDLETA.EQ.1) THEN
            EPS=RM(NR)/RR
            EPSS=SQRT(EPS)**3
            IF(NR.EQ.1) THEN
               QL= 0.25D0*(3.D0*Q0+QP(NR))
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
            ENDIF
            VTE=SQRT(TEL*RKEV/AME)
            RNUE=QL*RR/(TAUE*VTE*EPSS)
            RK33E=RK33/(1.D0+RA33*SQRT(ABS(RNUE))+RB33*RNUE)
     &                /(1.D0+RC33*RNUE*EPSS)
C
            FT     = 1.D0-SQRT(EPS)*RK33E
            ETA(NR)= ETA(NR)/FT
C
C        ****** NEOCLASSICAL RESISTIVITY PART II ******
C
         ELSEIF(MDLETA.EQ.2) THEN
            EPS=RM(NR)/RR
            EPSS=SQRT(EPS)**3
            IF(NR.EQ.1) THEN
               Q0=(4.D0*QP(1)-QP(2))/3.D0
               QL= 0.25D0*(3.D0*Q0+QP(NR))
               ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))
C               ZEFFL=ZEFF(NR)
            ELSE
               QL= 0.5D0*(QP(NR-1)+QP(NR))
               ZEFFL=0.5D0*(ZEFF(NR-1)+ZEFF(NR))
            ENDIF
C
         VTE=1.33D+7*DSQRT(TEL)
         FT=1.D0-(1.D0-EPS)**2.D0
     &         /(DSQRT(1.D0-EPS**2.D0)*(1.D0+1.46D0*DSQRT(EPS)))
         rLnLam=15.2D0-DLOG(ANE)/2+DLOG(TEL)
         TAUE=6.D0*PI*SQRT(2*PI)*AEPS0**2*DSQRT(AME)
     &             *(TEL*RKEV)**1.5D0/(ANE*1.D20*AEE**4*rLnLam)
         RNUSE=RR*QL/(VTE*TAUE*EPSS)
         PHI=FT/(1.D0+(0.58D0+0.2D0*ZEFFL)*RNUSE)                
         ETAS=1.65D-9*rLnLam/ABS(TEL)**1.5D0
         CH=0.56D0*(3.D0-ZEFFL)/((3.D0+ZEFFL)*ZEFFL)

         ETA(NR)=ETAS*ZEFFL*(1.D0+0.27D0*(ZEFFL-1.D0))
     &            /(1.D0-PHI)*(1.D0-CH*PHI)*(1.D0+0.47D0*(ZEFFL-1.D0))
         ENDIF 
  810    CONTINUE
C
         BP(1)=(RM(1)*DR/ETA(1))/RG(1)
         DO 820 NR=2,NRMAX
            BP(NR)=(RG(NR-1)*BP(NR-1)+RM(NR)*DR/ETA(NR))/RG(NR)
  820    CONTINUE
C
         BPS= AMYU0*RIPS*1.D6/(2.D0*PI*RA*FKAP)
         FACT=BPS/BP(NRMAX)
C
         DO 830 NR=1,NRMAX
            BP(NR)=FACT*BP(NR)
            QP(NR)=FKAP*RG(NR)*BB/(RR*BP(NR))
  830    CONTINUE
C
      ENDIF
C
      Q0=(4.D0*QP(1)-QP(2))/3.D0
C
C      WRITE(6,691) 'COEF',COEF
C      WRITE(6,691) 'RKEV',RKEV
C      WRITE(6,691) 'RTE ',(RT(NR,1),NR=1,NRMAX)
C      WRITE(6,691) 'ZEF ',(ZEFF(NR),NR=1,NRMAX)
C      WRITE(6,691) 'ETA ',(ETA(NR),NR=1,NRMAX)
C      WRITE(6,691) 'BP ',(BP(NR),NR=1,NRMAX)
C      WRITE(6,691) 'QP ',(QP(NR),NR=1,NRMAX)
  691 FORMAT(1H ,A3/
     &      (1H ,1P5E12.4))
C
      GRG(1)=0.0
      DO 900 NR=1,NRMAX
         GRM(NR)  =GCLIP(RM(NR))
         GRG(NR+1)=GCLIP(RG(NR))
  900 CONTINUE
C
      RETURN
      END
