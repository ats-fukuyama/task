C     $Id: wminit.f,v 1.53 2013/01/25 12:54:10 fukuyama Exp $
C
C     ****** INITIALIZE CONSTANTS & PARAMETERS ******
C
      SUBROUTINE WMINIT
C
      INCLUDE 'wmcomm.inc'
C
C     **** ALPHA PARTICLE PARAMETERS ****
C
C        PNA   : Alpha density at center               (1.0E20/Mm*3)
C        PNAL  : Density scale length                            (m)
C        PTA   : Effective temperature                         (keV)
C
      PNA  = 0.02D0
      PNAL = 0.5D0
      PTA  = 3.5D3
C
C     **** ZEFF PARAMETERS ****
C
C     ZEFF  : Effective Z (sum n Z^2 / sum n Z)
C
      ZEFF  = 2.D0
C
C     *** WAVE PARAMETER ***
C
C     RF    : Wave frequency                            (MHz)
C     RFI   : Wave growth rate                          (MHz)
C     RD    : Antenna minor radius                      (m)
C     BETAJ : Antenna current profile parameter
C     NTH0  : Central value of poloidal mode number
C     NPH0  : Central value of toroidal mode number
C     NHC   : Number of helical coils
C     PRFIN : Input Power (0 for given antenna current) (W)
C
      RF     = 50.D0
      RFI    = 0.D0
      CRF    = DCMPLX(RF,RFI)
      RD     = 1.1D0
      BETAJ  = 0.D0
      NTH0   = 0
      NPH0   = 8
      NHC    = 10
      PRFIN  = 0.D0
C
C     (Multi-mode case)
C
C
C     *** ANTENNA PARAMETERS ***
C
C        NAMAX : Number of antennae
C        AJ    : Antenna current density                       (A/m)
C        AEWGT : Waveguide electric field (poloidal)           (V/m)
C        AEWGZ : Waveguide electric field (toroidal)           (V/m)
C        APH   : Antenna phase                              (degree)
C        THJ1  : Start poloidal angle of antenna            (degree)
C        THJ2  : End poloidal angle of antenna              (degree)
C        PHJ1  : Start toroidal angle of antenna            (degree)
C        PHJ2  : End toroidal angle of antenna              (degree)
C        ANTANG: Antenna angle: 0 for vertical antenna or perp WG
C
      NAMAX  = 1
      DO NA=1,NAM
         AJ(NA)    = 1.D0
         AEWGT(NA) = 0.D0
         AEWGZ(NA) = 0.D0
         APH(NA)   = 0.D0
         THJ1(NA)  =-45.D0
         THJ2(NA)  = 45.D0
         PHJ1(NA)  = 0.D0
         PHJ2(NA)  = 0.D0
         ANTANG(NA)= 0.D0
      ENDDO
      MWGMAX = 0
C
C     *** MESH PARAMETERS ***
C
C        NRMAX  : Number of radial mesh points
C        NTHMAX : Number of poloidal mesh points
C        NHHMAX : Number of helically coupled toroidal modes
C                 =1 : axisymmetric calculation
C                 >1 : helical calculation
C        NPHMAX : Number of toroidal modes
C                 =1 : single toroidal mode calculation
C                 >1 : multi toroidal mode calculation (-NPHMAX/2+1..NPHMAX/2)
C
      NRMAX   = 50
      NTHMAX  = 1
      NHHMAX  = 1
      NPHMAX  = 1

      NSUMAX  = 64
      NSWMAX  = 64
C
C     NPH0L  : Current central value of toroidal mode number
C     PFRACL : Fraction of power for each toroidal mode number
C
      DO NPH = 1, NPHMAX
         NPH0L(NPH) = 8
         PFRACL(NPH)= 0.D0
      ENDDO
C
      DO NS=1,NSM
         DO NR=1,NRMAX
            MODELPR(NR,NS)=0
            MODELVR(NR,NS)=0
         ENDDO
      ENDDO
C
C     *** CONTROL PARAMETERS ***
C
C        NPRINT: Control print output
C                   0: No print out
C                   1: Minimum print out (without input data)
C                   2: Minimum print out (with input data)
C                   3: Standard print out
C                   4: More print out
C        NGRAPH: Control graphic output
C                   0: No graphic out
C                   1: Standard graphic out (2D: Coutour)
C                   2: Standard graphic out (2D: Paint)
C                   3: Standard graphic out (2D: Bird's eye)
C        MODELJ: Control antenna current model
C                   0: Loop antenna
C                   1: Waveguide
C                   2: Poloidal current
C                   3: Toroidal current
C                  2X: Vacuum eigen mode, poloidal current
C                  3X: Vacuum eigen mode, toroidal current
C
C        MODELA: Control alpha particle contribution
C                   0: No alpha effect
C                   1: Precession of alpha particles
C                   2: Precession of electrons
C                   3: Precession of both alpha particles and electrons
C                   4: Calculate alpha particle density using slowing down
C        MODELM: Control matrix solver
C                   0: BANDCD
C                   1: BANDCDB
C                   2: BSTABCDB
C                   8: BANDCDM
C                   9: BANDCDBM
C                  10: BSTABCDBM
C        MODELW: Control writing a data of absorped power
C                   0: Not writting
C                   1: Writting
C        MDMAX_IPS: Number of poloidal modes for dielectric tensor calculation
C        NDMAX_IPS: Number of toroidal modes for dielectric tensor calculation
C
      NPRINT = 2
      NGRAPH = 1
      MODELJ = 0
      MODELA = 0
      MODELM = 0
      MODELW = 0
      MDMAX_IPS=64
      NDMAX_IPS=8
C
C     *** EIGEN VALUE PARAMETERS ***
C
C        FRMIN : Minimum real part of frequency in amplitude scan
C        FRMAX : Maximum real part of frequency in amplitude scan
C        FIMIN : Minimum imag part of frequency in amplitude scan
C        FIMAX : Maximum imag part of frequency in amplitude scan
C        FI0   : Imag part of frequency in 1D amplitude scan
C
C        NGFMAX: Number of real freq mesh in 1D amplitude scan
C        NGXMAX: Number of real freq mesh in 2D amplitude scan
C        NGYMAX: Number of imag freq mesh in 2D amplitude scan
C
C        SCMIN : Minimum value in parameter scan
C        SCMAX : Maximum value in parameter scan
C        NSCMAX: Number of mesh in parameter scan
C
C        LISTEG: Listing in parameter scan
C
C        FRINI : Initial real part of frequency in Newton method
C        FIINI : Initial imag part of frequency in Newton method
C
C        DLTNW : Step size in evaluating derivatives in Newton method
C        EPSNW : Convergence criterion in Newton method
C        LMAXNW: Maximum iteration count in Newton method
C        LISTNW: Listing in Newton method
C        MODENW: Type of Newton method
C
C        NCONT : Number of contour lines
C        ILN1  : Line type of lower contours
C        IBL1  : Line boldness of lower contours
C        ICL1  : Line color of lower contours
C        ILN2  : Line type of higher contours
C        IBL2  : Line boldness of higher contours
C        ICL2  : Line color of higher contours
C
      FRMIN = 0.1D0
      FRMAX = 1.D0
      FIMIN =-0.1D0
      FIMAX = 0.1D0
      FI0   = 0.D0
C
      FRINI = DBLE(CRF)
      FIINI = DIMAG(CRF)
C
      NGFMAX= 11
      NGXMAX= 11
      NGYMAX= 11
C
      SCMIN = 0.1D0
      SCMAX = 1.D0
      NSCMAX= 11
C
      LISTEG= 1
C
      DLTNW = 1.D-6
      EPSNW = 1.D-6
      LMAXNW= 10
      LISTNW= 1
      MODENW= 0

      NCONT=30
C
C     *** ALFVEN FREQUENCY PARAMETERS ***
C
C        WAEMIN : Minimum frequency in Alfven frequency scan
C        WAEMAX : Maximum frequency in Alfven frequency scan
C
      WAEMIN = 0.001D0
      WAEMAX = 0.200D0
      KNAMEQ_SAVE=' '
      KNAMTR_SAVE=' '
C
      KNAMWM=' '

C    seki for vmec field strength 
      B0_FACT=1.0

      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE WMPARM(MODE,KIN,IERR)
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
      USE libkio
      INCLUDE 'wmcomm.inc'
C
      EXTERNAL WMNLIN,WMPLST
      CHARACTER KIN*(*)
C
    1 CALL TASK_PARM(MODE,'WM',KIN,WMNLIN,WMPLST,IERR)
      IF(IERR.NE.0) RETURN
C
      CALl WMCHEK(IERR)
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100
C
      RETURN
      END
C
C     ****** INPUT NAMELIST ******
C
      SUBROUTINE WMNLIN(NID,IST,IERR)
C
      INCLUDE 'wmcomm.inc'
C
      NAMELIST /WM/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,
     &              PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PU,PUS,NSMAX,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              PNA,PNAL,PTA,ZEFF,NDISP1,NDISP2,
     &              RF,RFI,RD,BETAJ,AJ,AEWGT,AEWGZ,
     &              APH,THJ1,THJ2,PHJ1,PHJ2,NAMAX,
     &              NRMAX,NTHMAX,NHHMAX,NTH0,NPH0,NHC,
     &              NPRINT,NGRAPH,MODELG,MODELJ,MODELP,MODELN,MODELA,
     &              MODELQ,MODELM,MODELW,MODELV,
     &              MWGMAX,MODEFR,MODEFW,ANTANG,
     &              FRMIN,FRMAX,FIMIN,FIMAX,FI0,FRINI,FIINI,
     &              NGFMAX,NGXMAX,NGYMAX,SCMIN,SCMAX,NSCMAX,LISTEG,
     &              DLTNW,EPSNW,LMAXNW,LISTNW,MODENW,NCONT,
     &              RHOMIN,QMIN,RHOEDG,
     &              RHOITB,PNITB,PTITB,PUITB,
     &              KNAMEQ,KNAMTR,KNAMWM,KNAMFP,KNAMFO,KNAMPF,
     &              WAEMIN,WAEMAX,PRFIN,MODELPR,MODELVR,
!     &              NSUMAX,NSWMAX,NPHMAX,NPH0L,PFRACL
! seki
     &              NSUMAX,NSWMAX,NPHMAX,NPH0L,PFRACL,
     &              MDMAX_IPS,NDMAX_IPS,
     &              B0_FACT
C
      RF=DREAL(CRF)
      RFI=DIMAG(CRF)
      READ(NID,WM,IOSTAT=IST,ERR=9800,END=9900)
      CRF=DCMPLX(RF,RFI)
      IF(NHHMAX.EQ.1) THEN
         NHHMAX2=1
      ELSE
         NHHMAX2=2*NHHMAX
      ENDIF
      IF(NTHMAX.EQ.1) THEN
         NTHMAX2=1
      ELSE
         NTHMAX2=2*NTHMAX
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
      SUBROUTINE WMPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &WM : BB,RR,RA,RB,Q0,QA,RKAP,RDLT,'/
     &       9X,'PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,'/
     &       9X,'PROFN1,PROFN2,PROFT1,PROFT2,ZEFF,'/
     &       9X,'NSMAX,PNA,PNAL,PTA,RF,RFI,RD,BETAJ,'/
     &       9X,'AJ,AEWGT,AEWGZ,APH,THJ1,THJ2,PHJ1,PHJ2,NAMAX,MWGMAX,'/
     &       9X,'NRMAX,NTHMAX,NHHMAX,NTH0,NPH0,NHC,'/
     &       9X,'MODELG,MODELJ,MODELP,MODELA,MODELN,'/
     &       9X,'MODELQ,MODELM,MODELW,,'/
     &       9X,'KNAMEQ,KNAMTR,KNAMPF,MODEFR,MODEFW,'/
     &       9X,'NPRINT,NGRAPH,PRFIN,MODELPR,MODELVR,'/
     &       9X,'FRMIN,FRMAX,FIMIN,FIMAX,FI0,'/
     &       9X,'FRINI,FIINI,NGFMAX,NGXMAX,NGYMAX,'/
     &       9X,'SCMIN,SCMAX,NSCMAX,LISTEG,ANTANG,'/
     &       9X,'DLTNW,EPSNW,LMAXNW,LISTNW,MODENW,NCONT,'/
     &       9X,'RHOMIN,QMIN,PU,PUS,PROFU1,PROFU2'/
     &       9X,'RHOITB,PNITB,PTITB,PUITB'/
     &       9X,'WAEMIN,WAEMAX,KNAMWM,KNAMFP,KNAMFO'/
!     &       9X,'NSUMAX,NSWMAX,NPHMAX,NPH0L,PFRACL')
!    seki
     &       9X,'NSUMAX,NSWMAX,NPHMAX,NPH0L,PFRACL'/
     &       9X,'MDMAX_IPS,NDMAX_IPS'/
     &       9X,'B0_FACT')
!    seki
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE WMCHEK(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      IERR=0
C
      IF(NSMAX.LT.0.OR.NSMAX.GT.NSM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NSMAX'
         WRITE(6,*) '                  NSMAX,NSM =',NSMAX,NSM
         IERR=1
      ENDIF
C
      IF((NAMAX.LT.1).OR.(NAMAX.GT.NAM)) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NAMAX'
         WRITE(6,*) '                  NAMAX,NAM =',NAMAX,NAM         
         IERR=1
      ENDIF
C
      IF(NTHMAX.LT.1.OR.NTHMAX.GT.MDM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : TOO LARGE NTHMAX'
         WRITE(6,*) '                  NTHMAX,MDM =',NTHMAX,MDM
         IERR=1
      ELSE
         MDP=NINT(LOG(DBLE(NTHMAX))/LOG(2.D0))
         IF(2**MDP.NE.NTHMAX) THEN
            WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTHMAX'
            WRITE(6,*) '                  NTHMAX,MDM =',NTHMAX,MDM
            IERR=1
         ENDIF
      ENDIF
C
      IF(NHHMAX.LT.1.OR.NHHMAX.GT.NDM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NHHMAX'
         WRITE(6,*) '                  NHHMAX,NDM =',NHHMAX,NDM
         IERR=1
      ELSE
         NDP=NINT(LOG(DBLE(NHHMAX))/LOG(2.D0))
         IF(2**NDP.NE.NHHMAX) THEN
            WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NHHMAX'
            WRITE(6,*) '                  NHHMAX,NDM =',NHHMAX,NDM
            IERR=1
         ENDIF
      ENDIF
C
      IF(NRMAX.LT.1.OR.NRMAX+1.GT.NRM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX,NRM =',NRMAX,NRM
         IERR=1
      ENDIF
C
      IF(MODELJ.LT.0) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL MODELJ'
         WRITE(6,*) '                  MODELJ =',MODELJ
         IERR=1
      ENDIF
C
      RETURN
      END
C
C     ****** DISPLAY INPUT DATA ******
C
      SUBROUTINE WMVIEW
C
      INCLUDE 'wmcomm.inc'
C
      IF(NPRINT.LT.2) RETURN
C
      IF(MODELG.EQ.0) THEN
         WRITE(6,*) '## MODELG=0: UNIFORM ##'
      ELSE IF(MODELG.EQ.1) THEN 
         WRITE(6,*) '## MODELG=1: CYLINDRICAL ##'
      ELSE IF(MODELG.EQ.2) THEN 
         WRITE(6,*) '## MODELG=2: TOROIDAL ##'
      ELSE IF(MODELG.EQ.3) THEN 
         WRITE(6,*) '## MODELG=3: TASK/EQ ##'
      ELSE IF(MODELG.EQ.4) THEN 
         WRITE(6,*) '## MODELG=4: VMEC ##'
      ELSE IF(MODELG.EQ.5) THEN 
         WRITE(6,*) '## MODELG=5: EQDSK ##'
      ELSE IF(MODELG.EQ.6) THEN 
         WRITE(6,*) '## MODELG=6: BOOZER ##'
      END IF
C
      IF(MODELJ.EQ.0) THEN
         WRITE(6,*) '## MODELJ=0: Loop antenna ##'
      ELSE IF(MODELJ.EQ.1) THEN 
         WRITE(6,*) '## MODELJ=1: Waveguide ##'
      ELSE IF(MODELJ.EQ.2) THEN 
         WRITE(6,*) '## MODELJ=2: POLOIDAL MODE ##'
      ELSE IF(MODELJ.EQ.3) THEN 
         WRITE(6,*) '## MODELJ=3: TOROIDAL MODE ##'
      ELSE
         WRITE(6,*) '## MODELJ=',MODELJ,': VACUUM EIGEN MODE ##'
      END IF
C
      IF(MODELN.EQ.7) THEN
         WRITE(6,*) '## MODELN=7: READ PROFILE DATA : WMDPRF ##'
      ELSEIF(MODELN.EQ.8) THEN
         WRITE(6,*) '## MODELN=8: READ PROFILE DATA : WMXPRF ##'
      ELSEIF(MODELN.EQ.9) THEN
         WRITE(6,*) '## MODELN=9: READ PROFILE DATA : TRDATA ##'
      ENDIF
C
      IF(MOD(MODELA,4).EQ.1) THEN
         WRITE(6,*) '## MODELA=1: ALPHA PARTICLE EFFECT ##'
      ELSE IF(MOD(MODELA,4).EQ.2) THEN
         WRITE(6,*) '## MODELA=2: ELECTRON BETA EFFECT ##'
      ELSE IF(MOD(MODELA,4).EQ.3) THEN
         WRITE(6,*) '## MODELA=3: ALPHA PARTICLE AND ELECTRON ',
     &              'BETA EFFECTS ##'
      ENDIF
      IF(MODELA.GE.4) THEN
         WRITE(6,*) '## MODELA=4: ALPHA PARTICLE DENSITY CALCULATED ##'
      ENDIF
C
      IF(MODELM.EQ.0) THEN
         WRITE(6,*) '## MODELM=0: BANDCD ##'
      ELSE IF(MODELM.EQ.1) THEN
         WRITE(6,*) '## MODELM=1: BANDCDB'
      ELSE IF(MODELM.EQ.2) THEN
         WRITE(6,*) '## MODELM=2: BSTABCDB ##'
      ELSE IF(MODELM.EQ.8) THEN
         WRITE(6,*) '## MODELM=8: BANDCDM ##'
      ELSE IF(MODELM.EQ.9) THEN
         WRITE(6,*) '## MODELM=9: BANDCDBM ##'
      ELSE IF(MODELM.EQ.10) THEN
         WRITE(6,*) '## MODELM=10: BSTABCDBM ##'
      ENDIF
C
      RF =DBLE(CRF)
      RFI=DIMAG(CRF)
      WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    ,
     &             'RA    ',RA    ,'RB    ',RB
      WRITE(6,601) 'Q0    ',Q0    ,'QA    ',QA    ,
     &             'RKAP  ',RKAP  ,'RDLT  ',RDLT
      WRITE(6,601) 'ZEFF  ',ZEFF  ,'PNA   ',PNA   ,
     &             'PNAL  ',PNAL  ,'PTA   ',PTA
      WRITE(6,601) 'RHOMIN',RHOMIN,'QMIN  ',QMIN  ,
     &             'PRFIN ',PRFIN
      WRITE(6,601) 'RF    ',RF    ,'RFI   ',RFI   ,
     &             'RD    ',RD    ,'BETAJ ',BETAJ
      WRITE(6,602) 'NRMAX ',NRMAX ,'NTHMAX',NTHMAX,
     &             'NHHMAX',NHHMAX
      WRITE(6,602) 'NTH0  ',NTH0  ,'NPH0  ',NPH0  ,
     &             'NHC   ',NHC   ,'MWGMAX',MWGMAX
      WRITE(6,602) 'MODELG',MODELG,'MODELJ',MODELJ,
     &             'MODELN',MODELN,'MODELA',MODELA
      WRITE(6,602) 'MODELM',MODELM,'MODELQ',MODELQ,
     &             'MODEFR',MODEFR,'MODEFW',MODEFW
      WRITE(6,603) 'NDMAX_IPS   ',NDMAX_IPS,
     &             'MDMAX_IPS   ',MDMAX_IPS
C
      WRITE(6,692)
      DO NS=1,NSMAX
         WRITE(6,611) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),
     &                   PTPR(NS),PTPP(NS),PTS(NS)
      ENDDO
      DO NS=1,NSMAX
         WRITE(6,612) NS,MODELP(NS),MODELV(NS),NDISP1(NS),NDISP2(NS),
     &          PZCL(NS),PU(NS),PUS(NS),PNITB(NS),PTITB(NS),PUITB(NS)
      ENDDO

      DO NS=1,NSMAX
         WRITE(6,613) NS,PROFN1(NS),PROFN2(NS),
     &                   PROFT1(NS),PROFT2(NS),
     &                   PROFU1(NS),PROFU2(NS),RHOITB(NS)
      END DO
C
      WRITE(6,693)
      DO NA=1,NAMAX
         WRITE(6,610) NA,AJ(NA),APH(NA),THJ1(NA),THJ2(NA),
     &                                  PHJ1(NA),PHJ2(NA)
         WRITE(6,614)    AEWGT(NA),AEWGZ(NA),ANTANG(NA)
      ENDDO
      RETURN
C
  601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  603 FORMAT(' ',A12,'=',I7,4X  :2X,A12,'=',I7,4X  :
     &        2X,A12,'=',I7)
  610 FORMAT(' ',I1,6(1PE11.3))
  611 FORMAT(' ',I1,7(1PE11.3))
  612 FORMAT(' ',I1,3I3,I2,6(1PE11.3))
  613 FORMAT(' ',I1,7(1PE11.3))
  614 FORMAT(' ',1X,6(1PE11.3))
  692 FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS',
     &                      8X,'PTPR',7X,'PTPP',7X,'PTS'/
     &       ' ','  MP MV ND1 ND2',2X,'PZCL',7X,'PU',9X,'PUS',
     &                      8X,'PNITB',6X,'PTITB',6X,'PUITB')
  693 FORMAT(' ','NA    AJ',9X,'APH',8X,'THJ1',7X,'THJ2',
     &                      7X,'PHJ1',7X,'PHJ2')
      END
C
C     ****** BROADCAST PARAMETERS ******
C
      SUBROUTINE WMPRBC
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION IPARA(25),DPARA(21)
C
      IF(NRANK.EQ.0) THEN
         RF=DBLE(CRF)
         RFI=DIMAG(CRF)
         IPARA(1) =NSMAX
         IPARA(2) =NAMAX
         IPARA(3) =NRMAX
         IPARA(4) =NTHMAX
         IPARA(5) =NHHMAX
         IPARA(6) =NTH0
         IPARA(7) =NPH0
         IPARA(8) =NHC
         IPARA(9) =NPRINT
         IPARA(10)=NGRAPH
         IPARA(11)=MODELG
         IPARA(12)=MODELJ
         IPARA(13)=MODELN
         IPARA(14)=MODELA
         IPARA(15)=MODELK
         IPARA(16)=MODELM
         IPARA(17)=MODELW
         IPARA(18)=LISTEG
         IPARA(19)=LMAXNW
         IPARA(20)=LISTNW
         IPARA(21)=MODENW
         IPARA(22)=NCONT
         IPARA(23)=NPHMAX
         IPARA(24)=NDMAX_IPS
         IPARA(25)=MDMAX_IPS
C
         DPARA(1) =BB
         DPARA(2) =RR
         DPARA(3) =RA
         DPARA(4) =RB
         DPARA(5) =Q0
         DPARA(6) =QA
         DPARA(7) =RKAP
         DPARA(8) =RDLT
         DPARA(9) =ZEFF
         DPARA(10)=PNA
         DPARA(11)=PNAL
         DPARA(12)=PTA
         DPARA(13)=RF
         DPARA(14)=RFI
         DPARA(15)=RD
         DPARA(16)=BETAJ
         DPARA(17)=DLTNW
         DPARA(18)=EPSNW
         DPARA(19)=RHOMIN
         DPARA(20)=QMIN
         DPARA(21)=PRFIN
      ENDIF
C
      CALL mtx_broadcast_integer(IPARA,25)
      CALL mtx_broadcast_real8(DPARA,21)
C
      IF(NRANK.NE.0) THEN
         NSMAX =IPARA(1) 
         NAMAX =IPARA(2) 
         NRMAX =IPARA(3) 
         NTHMAX=IPARA(4) 
         NHHMAX=IPARA(5) 
         NTH0  =IPARA(6) 
         NPH0  =IPARA(7) 
         NHC   =IPARA(8) 
         NPRINT=IPARA(9) 
         NGRAPH=IPARA(10)
         MODELG=IPARA(11)
         MODELJ=IPARA(12)
         MODELN=IPARA(13)
         MODELA=IPARA(14)
         MODELK=IPARA(15)
         MODELM=IPARA(16)
         MODELW=IPARA(17)
         LISTEG=IPARA(18)
         LMAXNW=IPARA(19)
         LISTNW=IPARA(20)
         MODENW=IPARA(21)
         NCONT =IPARA(22)
         NPHMAX=IPARA(23)
         NDMAX_IPS=IPARA(24)
         MDMAX_IPS=IPARA(25)
C       
         BB    =DPARA(1) 
         RR    =DPARA(2) 
         RA    =DPARA(3) 
         RB    =DPARA(4) 
         Q0    =DPARA(5) 
         QA    =DPARA(6) 
         RKAP  =DPARA(7) 
         RDLT  =DPARA(8) 
         ZEFF  =DPARA(9)
         PNA   =DPARA(10)
         PNAL  =DPARA(11)
         PTA   =DPARA(12)
         RF    =DPARA(13)
         RFI   =DPARA(14)
         RD    =DPARA(15)
         BETAJ =DPARA(16)
         DLTNW =DPARA(17)
         EPSNW =DPARA(18)
         RHOMIN=DPARA(19)
         QMIN  =DPARA(20)
         PRFIN =DPARA(21)
         CRF=DCMPLX(RF,RFI)
      ENDIF
C
      CALL mtx_broadcast_real8(PA,NSMAX)
      CALL mtx_broadcast_real8(PZ,NSMAX)
      CALL mtx_broadcast_real8(PN,NSMAX)
      CALL mtx_broadcast_real8(PNS,NSMAX)
      CALL mtx_broadcast_real8(PZCL,NSMAX)
      CALL mtx_broadcast_real8(PTPR,NSMAX)
      CALL mtx_broadcast_real8(PTPP,NSMAX)
      CALL mtx_broadcast_real8(PTS,NSMAX)
      CALL mtx_broadcast_real8(PU,NSMAX)
      CALL mtx_broadcast_real8(PUS,NSMAX)
      CALL mtx_broadcast_real8(PNITB,NSMAX)
      CALL mtx_broadcast_real8(PTITB,NSMAX)
      CALL mtx_broadcast_real8(PUITB,NSMAX)
      CALL mtx_broadcast_real8(PROFN1,NSMAX)
      CALL mtx_broadcast_real8(PROFN2,NSMAX)
      CALL mtx_broadcast_real8(PROFT1,NSMAX)
      CALL mtx_broadcast_real8(PROFT2,NSMAX)
      CALL mtx_broadcast_real8(PROFU1,NSMAX)
      CALL mtx_broadcast_real8(PROFU2,NSMAX)
      CALL mtx_broadcast_integer(MODELP,NSMAX)
      CALL mtx_broadcast_real8(AJ,NAMAX)
      CALL mtx_broadcast_real8(AEWGT,NAMAX)
      CALL mtx_broadcast_real8(AEWGZ,NAMAX)
      CALL mtx_broadcast_real8(APH,NAMAX)
      CALL mtx_broadcast_real8(THJ1,NAMAX)
      CALL mtx_broadcast_real8(THJ2,NAMAX)
      CALL mtx_broadcast_real8(PHJ1,NAMAX)
      CALL mtx_broadcast_real8(PHJ2,NAMAX)
      CALL mtx_broadcast_real8(ANTANG,NAMAX)
      CALL mtx_broadcast_character(KNAMEQ,80)
      CALL mtx_broadcast_character(KNAMTR,80)
      CALL mtx_broadcast_character(KNAMPF,80)
C
      RETURN
      END
