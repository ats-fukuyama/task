! fpparm.f90

module fpparm

  PRIVATE
  PUBLIC fp_parm,fp_broadcast,fp_view

contains

!     ***********************
!          PARAMETER INPUT
!     ***********************

  SUBROUTINE fp_parm(mode,kin,ierr)

!     mode=0 : standard namelinst input
!     mode=1 : namelist file input
!     mode=2 : namelist line input

!     ierr=0 : normal end
!     ierr=1 : namelist standard input error
!     ierr=2 : namelist file does not exist
!     ierr=3 : namelist file open error
!     ierr=4 : namelist file read error
!     ierr=5 : namelist file abormal end of file
!     ierr=6 : namelist line input error
!     ierr=7 : unknown MODE
!     ierr=10X : input parameter out of range

      USE plcomm,ONLY: MODEL_PROF,NSMAX, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
      IMPLICIT NONE
      INTEGER,INTENT(IN):: mode
      CHARACTER(LEN=*),INTENT(IN)::  kin
      INTEGER,INTENT(OUT):: ierr
      INTEGER:: NS
      EXTERNAL task_parm

    1 CALL task_parm(mode,'FP',kin,fp_nlin,fp_plst,ierr)
      IF(ierr.NE.0) RETURN

      IF(MODEL_PROF.EQ.0) THEN
         DO NS=1,NSMAX
            PROFN1(NS)=PROFN1(1)
            PROFN2(NS)=PROFN2(1)
            PROFT1(NS)=PROFT1(1)
            PROFT2(NS)=PROFT2(1)
            PROFU1(NS)=PROFU1(1)
            PROFU2(NS)=PROFu2(1)
         END DO
      END IF

      CALL fp_check(ierr)
      IF(mode.EQ.0.AND.ierr.NE.0) GO TO 1
      IF(ierr.NE.0) ierr=ierr+100

      RETURN
  END SUBROUTINE fp_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE fp_nlin(nid,ist,ierr)

      use fpcomm_parm
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: nid
      INTEGER,INTENT(OUT) :: ist,ierr

      NAMELIST /FP/ &
           NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW, &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX, &
           NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
           PNITB,PTITB,PUITB,PZCL, &
           KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF, &
           KNAMFO,KNAMTR,KNAMEQ2,KID_NS,ID_NS, &
           NSAMAX,NSBMAX,NS_NSA,NS_NSB, &
           LMAX_WR,NRAYS_WR,NRAYE_WR,NCMIN,NCMAX,NBEAMMAX,NSSPB,NSSPF, &
           NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX, &
           NTMAX,NTSTEP_COEF,NTSTEP_COLL, &
           NTG1STEP,NTG1MIN,NTG1MAX, &
           NTG2STEP,NTG2MIN,NTG2MAX, &
           MODELE,MODELA,MODELC,MODELR,MODELS,MODELW, &
           MODELD,MODELD_RDEP,MODELD_PDEP,MODELD_EDGE, &
           MODELD_PINCH,MODELD_BOUNDARY,MODELD_CDBM, &
           MODEL_LOSS,MODEL_SYNCH,MODEL_NBI,MODEL_WAVE, &
           IMTX,MODEL_KSP,MODEL_PC,LMAXFP,LMAXE, &
           NGLINE,NGRAPH,LLMAX,LLMAX_NF,IDBGFP, &
           MODEL_DISRUPT,MODEL_Connor_fp,MODEL_BS,MODEL_jfp, &
           MODEL_LNL,MODEL_RE_pmax,MODELD_n_RE,MODEL_IMPURITY, &
           MODEL_SINK,N_IMPU,MODEL_DELTA_F, &
           N_partition_r,N_partition_s,N_partition_p, &
           PMAX,PMAX_BB,EMAX, &
           R1,DELR1,RMIN,RMAX,E0,ZEFF, &
           PABS_LH,PABS_FW,PABS_EC,PABS_wr,PABS_WM,RF_WM, &
           FACT_WM,FACT_WR,FACT_NRAY,DELNPR_WR,DELNPR_WM,EPS_WR,DELY_WR, &
           DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC, &
           DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW, &
           CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH, &
           SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG, &
           SPFTOT,SPFR0,SPFRW,SPFENG, &
           DRR0,DRRS,FACTOR_CDBM,DRR_EDGE,RHO_EDGE, &
           FACTOR_DRR_EDGE,FACTOR_PINCH,deltaB_B,TLOSS, &
           DELT,RIMPL,EPSFP,EPSM,EPSE,EPSDE,H0DE, &
           PGMAX,RGMAX,RGMIN, &
           T0_quench,tau_quench,tau_mgi, &
           time_quench_start,RJPROF1,RJPROF2, &
           v_RE,target_zeff,SPITOT, MODEL_EX_READ_Tn, MODEL_DELTA_F_NI_RATIO, MODEL_DELTA_F_CN,
           FACT_BULK, time_exp_offset, MODEL_BULK_CONST, RN_NEU0, MODEL_CX_LOSS, RN_NEUS, &
           EG_NAME_TMS, EG_NAME_CX, SV_FILE_NAME_H, SV_FILE_NAME_D, NSA_F1, NTH_F1, NR_F1, &
           OUTPUT_TXT_F1, OUTPUT_TXT_DELTA_F, OUTPUT_TXT_HEAT_PROF, OUTPUT_TXT_BEAM_WIDTH, &
           OUTPUT_TXT_BEAM_DENS, NI_RATIO, MODEL_NF_CS, NF_IDMAX, &
           DH_RATIO, HHe_RATIO, DT_RATIO, given_zeff, OUTPUT_NFID, SHOT_NUMBER, EG_PATH,&
           MODELS_full, MODELS_bt, MODELS_tt

      READ(nid,FP,IOSTAT=ist,ERR=9800,END=9900)

      ierr=0
      RETURN

 9800 ierr=8
      RETURN
 9900 ierr=9
      RETURN
  END SUBROUTINE fp_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE fp_plst

      WRITE(6,*) '&FP : NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW,'
      WRITE(6,*) '      RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'
      WRITE(6,*) '      PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
      WRITE(6,*) '      RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX,'
      WRITE(6,*) '      NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,'
      WRITE(6,*) '      PNITB,PTITB,PUITB,PZCL,'
      WRITE(6,*) '      KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,'
      WRITE(6,*) '      KNAMFO,KNAMTR,KNAMEQ2,KID_NS,ID_NS,'
      WRITE(6,*) '      NSAMAX,NSBMAX,NS_NSA,NS_NSB,'
      WRITE(6,*) '      LMAX_WR,NRAYS_WR,NRAYE_WR,'
      WRITE(6,*) '      NCMIN,NCMAX,NBEAMMAX,NSSPB,NSSPF,'
      WRITE(6,*) '      NPMAX,NTHMAX,NRMAX,NAVMAX,NP2MAX,'
      WRITE(6,*) '      NTMAX,NTSTEP_COEF,NTSTEP_COLL,'
      WRITE(6,*) '      NTG1STEP,NTG1MIN,NTG1MAX,'
      WRITE(6,*) '      NTG2STEP,NTG2MIN,NTG2MAX,'
      WRITE(6,*) '      MODELE,MODELA,MODELC,MODELR,MODELS,MODELW,'
      WRITE(6,*) '      MODELD,MODELD_RDEP,MODELD_PDEP,MODELD_EDGE,'
      WRITE(6,*) '      MODELD_BOUNDARY,MODELD_CDBM,MODELD_PINCH,'
      WRITE(6,*) '      MODEL_LOSS,MODEL_SYNCH,MODEL_NBI,MODEL_WAVE,'
      WRITE(6,*) '      IMTX,MODEL_KSP,MODEL_PC,LMAXFP,LMAXE,'
      WRITE(6,*) '      NGLINE,NGRAPH,LLMAX,LLMAX_NF,IDBGFP,'
      WRITE(6,*) '      MODEL_DISRUPT,MODEL_Connor_fp,MODEL_BS,MODEL_jfp,'
      WRITE(6,*) '      MODEL_LNL,MODEL_RE_pmax,MODELD_n_RE,MODEL_IMPURITY,'
      WRITE(6,*) '      MODEL_SINK,N_IMPU,MODEL_DELTA_F'
      WRITE(6,*) '      N_partition_r,N_partition_s,N_partition_p,'
      WRITE(6,*) '      PMAX,PMAX_BB,EMAX'
      WRITE(6,*) '      R1,DELR1,RMIN,RMAX,E0,ZEFF,'
      WRITE(6,*) '      PABS_LH,PABS_FW,PABS_EC,PABS_WR,PABS_WM,RF_WM,'
      WRITE(6,*) '      FACT_WM,FACT_WR,DELNPR_WR,DELNPR_WM,EPS_WR,DELY_WR,'
      WRITE(6,*) '      FACT_NRAY,Y0_WM,DELY_WM,'
      WRITE(6,*) '      DEC,PEC1,PEC2,PEC3,PEC4,RFEC,DELYEC,'
      WRITE(6,*) '      DLH,PLH1,PLH2,RLH,DFW,PFW1,PFW2,RFW,'
      WRITE(6,*) '      CEWR,CEWTH,CEWPH,RKWR,RKWTH,RKWPH,'
      WRITE(6,*) '      SPBTOT,SPBR0,SPBRW,SPBENG,SPBANG,SPBPANG,'
      WRITE(6,*) '      SPFTOT,SPFR0,SPFRW,SPFENG,'
      WRITE(6,*) '      DRR0,DRRS,FACTOR_CDBM,DRR_EDGE,RHO_EDGE,'
      WRITE(6,*) '      FACTOR_DRR_EDGE,FACTOR_PINCH,deltaB_B,TLOSS,'
      WRITE(6,*) '      DELT,RIMPL,EPSFP,EPSM,EPSE,EPSDE,H0DE,'
      WRITE(6,*) '      PGMAX,RGMAX,RGMIN,'
      WRITE(6,*) '      T0_quench,tau_quench,tau_mgi,'
      WRITE(6,*) '      time_quench_start,RJPROF1,RJPROF2,'
      WRITE(6,*) '      v_RE,target_zeff,SPITOT,MODEL_EX_READ_Tn, MODEL_EX_READ_DH_RATIO, FACT_BULK'
      WRITE(6,*) '      time_exp_offset, MODEL_BULK_CONST, RN_NEU0, MODEL_CX_LOSS, RN_NEUS'
      WRITE(6,*) '      EG_NAME_TMS, EG_NAME_CX, SV_FILE_NAME_H, SV_FILE_NAME_D, NSA_F1, NTH_F1, NR_F1'
      WRITE(6,*) '      OUTPUT_TXT_F1, OUTPUT_TXT_DELTA_F, OUTPUT_TXT_HEAT_PROF, OUTPUT_TXT_BEAM_WIDTH'
      WRITE(6,*) '      OUTPUT_TXT_BEAM_DENS, NI_RATIO'

      RETURN
  END SUBROUTINE fp_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE fp_check(ierr)

      USE fpcomm
      IMPLICIT NONE
      integer,intent(out):: ierr

      ierr=0

      RETURN
  END SUBROUTINE fp_check

!     ***** BROADCAST INPUT PARAMETERS *****

  SUBROUTINE fp_broadcast

      USE fpcomm_parm
      USE libmpi
      USE libmtx
      IMPLICIT NONE
      INTEGER,DIMENSION(99):: idata
      real(8),DIMENSION(99):: rdata
      complex(8),DIMENSION(3):: cdata
      INTEGER:: NS

!----- PL input parameters -----     

      idata( 1)=NSMAX
      idata( 2)=MODELG
      idata( 3)=MODELN
      idata( 4)=MODELQ
      idata( 5)=IDEBUG
      idata( 6)=MODEFR
      idata( 7)=MODEFW
      idata( 8)=MODEL_PROF
      idata( 9)=MODEL_PROF

      CALL mtx_broadcast_integer(idata,7)
      NSMAX =idata( 1)
      MODELG=idata( 2)
      MODELN=idata( 3)
      MODELQ=idata( 4)
      IDEBUG=idata( 5)
      MODEFR=idata( 6)
      MODEFW=idata( 7)

      rdata( 1)=RR
      rdata( 2)=RA
      rdata( 3)=RB
      rdata( 4)=RKAP
      rdata( 5)=RDLT
      rdata( 6)=BB
      rdata( 7)=Q0
      rdata( 8)=QA
      rdata( 9)=RIP
      rdata(10)=PROFJ
      rdata(11)=RHOMIN
      rdata(12)=QMIN
      rdata(13)=RHOEDG
      rdata(14)=RHOGMN
      rdata(15)=RHOGMX

      CALL mtx_broadcast_real8(rdata,15)
      RR    =rdata( 1)
      RA    =rdata( 2)
      RB    =rdata( 3)
      RKAP  =rdata( 4)
      RDLT  =rdata( 5)
      BB    =rdata( 6)
      Q0    =rdata( 7)
      QA    =rdata( 8)
      RIP   =rdata( 9)
      PROFJ =rdata(10)
      RHOMIN=rdata(11)
      QMIN  =rdata(12)
      RHOEDG=rdata(13)
      RHOGMN=rdata(14)
      RHOGMX=rdata(15)

      CALL mtx_broadcast_integer(NPA,NSMAX)
      CALL mtx_broadcast_integer(ID_NS,NSMAX)
      CALL mtx_broadcast_real8(PA,NSMAX)
      CALL mtx_broadcast_real8(PZ,NSMAX)
      CALL mtx_broadcast_real8(PN,NSMAX)
      CALL mtx_broadcast_real8(PNS,NSMAX)
      CALL mtx_broadcast_real8(PTPR,NSMAX)
      CALL mtx_broadcast_real8(PTPP,NSMAX)
      CALL mtx_broadcast_real8(PTS,NSMAX)
      CALL mtx_broadcast_real8(PU,NSMAX)
      CALL mtx_broadcast_real8(PUS,NSMAX)
      CALL mtx_broadcast_real8(RHOITB,NSMAX)
      CALL mtx_broadcast_real8(PNITB,NSMAX)
      CALL mtx_broadcast_real8(PTITB,NSMAX)
      CALL mtx_broadcast_real8(PUITB,NSMAX)
      CALL mtx_broadcast_real8(PROFN1,NSMAX)
      CALL mtx_broadcast_real8(PROFN2,NSMAX)
      CALL mtx_broadcast_real8(PROFT1,NSMAX)
      CALL mtx_broadcast_real8(PROFT2,NSMAX)
      CALL mtx_broadcast_real8(PROFU1,NSMAX)
      CALL mtx_broadcast_real8(PROFU2,NSMAX)
      CALL mtx_broadcast_real8(PZCL,NSMAX)

      CALL mtx_broadcast_character(KNAMEQ,80)
      CALL mtx_broadcast_character(KNAMWR,80)
      CALL mtx_broadcast_character(KNAMFP,80)
      CALL mtx_broadcast_character(KNAMWM,80)
      CALL mtx_broadcast_character(KNAMPF,80)
      CALL mtx_broadcast_character(KNAMFO,80)
      CALL mtx_broadcast_character(KNAMTR,80)
      CALL mtx_broadcast_character(KNAMEQ2,80)
      CALL mtx_broadcast_character(EG_NAME_TMS,80)
      CALL mtx_broadcast_character(EG_NAME_CX,80)
      CALL mtx_broadcast_character(SV_FILE_NAME_H,80)
      CALL mtx_broadcast_character(SV_FILE_NAME_D,80)

      DO NS=1,NSMAX
         CALL mtx_broadcast_character(KID_NS(NS),2)
      END DO

!----- FP input parameters -----

      idata( 1)=NSAMAX
      idata( 2)=NSBMAX
      idata( 3)=LMAX_WR
      idata( 4)=NBEAMMAX
      idata( 5)=NSSPF
      idata( 6)=NPMAX
      idata( 7)=NTHMAX
      idata( 8)=NRMAX
      idata( 9)=NAVMAX
      idata(10)=NP2MAX

      idata(11)=NTMAX
      idata(12)=NTSTEP_COEF
      idata(13)=NTSTEP_COLL
      idata(14)=NTG1STEP
      idata(15)=NTG1MIN
      idata(16)=NTG1MAX
      idata(17)=NTG2STEP
      idata(18)=NTG2MIN
      idata(19)=NTG2MAX

      idata(20)=MODELE
      idata(21)=MODELA
!      idata(22)=MODELC
      idata(22)=0.D0
      idata(23)=MODELR
      idata(24)=MODELS
      idata(25)=MODELD
      idata(26)=MODELD_RDEP
      idata(27)=MODELD_PDEP
      idata(28)=MODELD_EDGE
      idata(29)=MODELD_PINCH
      idata(30)=MODELD_BOUNDARY
      idata(31)=MODELD_CDBM
      idata(32)=MODEL_LOSS
      idata(33)=MODEL_SYNCH
      idata(34)=MODEL_NBI
      idata(35)=MODEL_WAVE
      idata(36)=IMTX
      idata(37)=MODEL_KSP
      idata(38)=MODEL_PC

      idata(39)=LMAXFP
      idata(40)=LMAXE
      idata(41)=NGLINE
      idata(42)=NGRAPH
      idata(43)=LLMAX
      idata(44)=LLMAX_NF
      idata(45)=IDBGFP

      idata(46)=MODEL_DISRUPT
      idata(47)=MODEL_Connor_FP
      idata(48)=MODEL_BS
      idata(49)=MODEL_jfp
      idata(50)=MODEL_LNL
      idata(51)=MODEL_RE_pmax
      idata(52)=MODELD_n_RE
      idata(53)=MODEL_IMPURITY
      idata(54)=MODEL_SINK
      idata(55)=n_impu
      idata(56)=N_partition_r
      idata(57)=N_partition_s
      idata(58)=N_partition_p
      idata(59)=MODEL_EX_READ_Tn
      idata(60)=MODEL_BULK_CONST
      idata(61)=MODEL_CX_LOSS
      idata(62)=NSA_F1
      idata(63)=NTH_F1
      idata(64)=NR_F1
      idata(65)=MODEL_EX_READ_DH_RATIO
      idata(66)=OUTPUT_TXT_F1
      idata(67)=OUTPUT_TXT_DELTA_F
      idata(68)=OUTPUT_TXT_HEAT_PROF
      idata(69)=OUTPUT_TXT_BEAM_WIDTH
      idata(70)=OUTPUT_TXT_BEAM_DENS
      idata(71)=NRAYS_WR
      idata(72)=NRAYE_WR

      CALL mtx_broadcast_integer(idata,72)
      NSAMAX         =idata( 1)
      NSBMAX         =idata( 2)
      LMAX_WR        =idata( 3)
      NBEAMMAX       =idata( 4)
      NSSPF          =idata( 5)
      NPMAX          =idata( 6)
      NTHMAX         =idata( 7)
      NRMAX          =idata( 8)
      NAVMAX         =idata( 9)
      NP2MAX         =idata(10)

      NTMAX          =idata(11)
      NTSTEP_COEF    =idata(12)
      NTSTEP_COLL    =idata(13)
      NTG1STEP       =idata(14)
      NTG1MIN        =idata(15)
      NTG1MAX        =idata(16)
      NTG2STEP       =idata(17)
      NTG2MIN        =idata(18)
      NTG2MAX        =idata(19)

      MODELE         =idata(20)
      MODELA         =idata(21)
!      MODELC         =idata(22)
      MODELR         =idata(23)
      MODELS         =idata(24)
      MODELD         =idata(25)
      MODELD_RDEP    =idata(26)
      MODELD_PDEP    =idata(27)
      MODELD_EDGE    =idata(28)
      MODELD_PINCH   =idata(29)
      MODELD_BOUNDARY=idata(30)
      MODELD_CDBM    =idata(31)
      MODEL_LOSS     =idata(32)
      MODEL_SYNCH    =idata(33)
      MODEL_NBI      =idata(34)
      MODEL_WAVE     =idata(35)
      IMTX           =idata(36)
      MODEL_KSP      =idata(37)
      MODEL_PC       =idata(38)

      LMAXFP         =idata(39)
      LMAXE          =idata(40)
      NGLINE         =idata(41)
      NGRAPH         =idata(42)
      LLMAX          =idata(43)
      LLMAX_NF       =idata(44)
      IDBGFP         =idata(45)
      MODEL_DISRUPT  =idata(46)
      MODEL_Connor_FP=idata(47)
      MODEL_BS       =idata(48)
      MODEL_jfp      =idata(49)
      MODEL_LNL      =idata(50)
      MODEL_RE_pmax  =idata(51)
      MODELD_n_RE    =idata(52)
      MODEL_IMPURITY =idata(53)
      MODEL_SINK     =idata(54)
      n_impu         =idata(55)
      N_partition_r  =idata(56)
      N_partition_s  =idata(57)
      N_partition_p  =idata(58)
      MODEL_EX_READ_Tn  =idata(59)
      MODEL_BULK_CONST =idata(60)
      MODEL_CX_LOSS  =idata(61)
      NSA_F1 = idata(62)
      NTH_F1 = idata(63)
      NR_F1 = idata(64)
      MODEL_EX_READ_DH_RATIO=idata(65)
      OUTPUT_TXT_F1=idata(66)
      OUTPUT_TXT_DELTA_F=idata(67)
      OUTPUT_TXT_HEAT_PROF=idata(68)
      OUTPUT_TXT_BEAM_WIDTH=idata(69)
      OUTPUT_TXT_BEAM_DENS=idata(70)
      NRAYS_WR=idata(71)
      NRAYE_WR=idata(72)

      CALL mtx_broadcast_integer(NS_NSA,NSAMAX)
      CALL mtx_broadcast_integer(NS_NSB,NSBMAX)
      CALL mtx_broadcast_integer(NCMIN,NSAMAX)
      CALL mtx_broadcast_integer(NCMAX,NSAMAX)
      CALL mtx_broadcast_integer(NSSPB,NBEAMMAX)
      CALL mtx_broadcast_integer(MODELW,NSMAX)
      CALL mtx_broadcast_integer(MODELC,NSMAX)
      CALL mtx_broadcast_integer(MODEL_DELTA_F,NSMAX)

      rdata( 1)=R1
      rdata( 2)=DELR1
      rdata( 3)=RMIN
      rdata( 4)=RMAX
      rdata( 5)=E0
      rdata( 6)=ZEFF
      rdata( 7)=PABS_LH
      rdata( 8)=PABS_FW
      rdata( 9)=PABS_EC
      rdata(10)=PABS_WR

      rdata(11)=PABS_WM
      rdata(12)=RF_WM
      rdata(13)=FACT_WM
      rdata(14)=FACT_WR
      rdata(15)=DELNPR_WR
      rdata(16)=DELNPR_WM
      rdata(17)=EPS_WR
      rdata(18)=DELY_WR
      rdata(19)=Y0_WM
      rdata(20)=DELY_WM

      rdata(21)=DEC
      rdata(22)=PEC1
      rdata(23)=PEC2
      rdata(24)=PEC3
      rdata(25)=PEC4
      rdata(26)=RFEC
      rdata(27)=DELYEC
      rdata(28)=DLH
      rdata(29)=PLH1
      rdata(30)=PLH2

      rdata(31)=RLH
      rdata(32)=DFW
      rdata(33)=PFW1
      rdata(34)=PFW2
      rdata(35)=RFW
      rdata(36)=RKWR
      rdata(37)=RKWTH
      rdata(38)=RKWPH
      rdata(39)=SPFTOT
      rdata(40)=SPFR0

      rdata(41)=SPFRW
      rdata(42)=SPFENG
      rdata(43)=DRR0
      rdata(44)=DRRS
      rdata(45)=FACTOR_CDBM
      rdata(46)=DRR_EDGE
      rdata(47)=RHO_EDGE
      rdata(48)=FACTOR_DRR_EDGE
      rdata(49)=FACTOR_PINCH
      rdata(50)=DELTAB_B

      rdata(51)=DELT
      rdata(52)=RIMPL
      rdata(53)=EPSFP
      rdata(54)=EPSM
      rdata(55)=EPSE
      rdata(56)=EPSDE
      rdata(57)=H0DE
      rdata(58)=PGMAX
      rdata(59)=RGMAX
      rdata(60)=RGMIN

      rdata(61)=T0_quench
      rdata(62)=tau_quench
      rdata(63)=tau_mgi
      rdata(64)=time_quench_start
      rdata(65)=RJPROF1
      rdata(66)=RJPROF2
      rdata(67)=v_RE
      rdata(68)=target_zeff
      rdata(69)=SPITOT
      rdata(70)=FACT_BULK
      rdata(71)=time_exp_offset
      rdata(72)=RN_NEU0
      rdata(73)=RN_NEUS

      CALL mtx_broadcast_real8(rdata,73)

      R1               =rdata( 1)
      DELR1            =rdata( 2)
      RMIN             =rdata( 3)
      RMAX             =rdata( 4)
      E0               =rdata( 5)
      ZEFF             =rdata( 6)
      PABS_LH          =rdata( 7)
      PABS_FW          =rdata( 8)
      PABS_EC          =rdata( 9)
      PABS_WR          =rdata(10)

      PABS_WM          =rdata(11)
      RF_WM            =rdata(12)
      FACT_WM          =rdata(13)
      FACT_WR          =rdata(14)
      DELNPR_WR        =rdata(15)
      DELNPR_WM        =rdata(16)
      EPS_WR           =rdata(17)
      DELY_WR          =rdata(18)
      Y0_WM            =rdata(19)
      DELY_WM          =rdata(20)

      DEC              =rdata(21)
      PEC1             =rdata(22)
      PEC2             =rdata(23)
      PEC3             =rdata(24)
      PEC4             =rdata(25)
      RFEC             =rdata(26)
      DELYEC           =rdata(27)
      DLH              =rdata(28)
      PLH1             =rdata(29)
      PLH2             =rdata(30)

      RLH              =rdata(31)
      DFW              =rdata(32)
      PFW1             =rdata(33)
      PFW2             =rdata(34)
      RFW              =rdata(35)
      RKWR             =rdata(36)
      RKWTH            =rdata(37)
      RKWPH            =rdata(38)
      SPFTOT           =rdata(39)
      SPFR0            =rdata(40)

      SPFRW            =rdata(41)
      SPFENG           =rdata(42)
      DRR0             =rdata(43)
      DRRS             =rdata(44)
      FACTOR_CDBM      =rdata(45)
      DRR_EDGE         =rdata(46)
      RHO_EDGE         =rdata(47)
      FACTOR_DRR_EDGE  =rdata(48)
      FACTOR_PINCH     =rdata(49)
      DELTAB_B         =rdata(50)

      DELT             =rdata(51)
      RIMPL            =rdata(52)
      EPSFP            =rdata(53)
      EPSM             =rdata(54)
      EPSE             =rdata(55)
      EPSDE            =rdata(56)
      H0DE             =rdata(57)
      PGMAX            =rdata(58)
      RGMAX            =rdata(59)
      RGMIN            =rdata(60)

      T0_quench        =rdata(61)
      tau_quench       =rdata(62)
      tau_mgi          =rdata(63)
      time_quench_start=rdata(64)
      RJPROF1          =rdata(65)
      RJPROF2          =rdata(66)
      v_RE             =rdata(67)
      target_zeff      =rdata(68)
      SPITOT           =rdata(69)
      FACT_BULK        =rdata(70)

      time_exp_offset  =rdata(71)
      RN_NEU0          =rdata(72)
      RN_NEUS          =rdata(73)

      CALL mtx_broadcast_real8(pmax,NSMAX)
      CALL mtx_broadcast_real8(pmax_bb,NSMAX)
      CALL mtx_broadcast_real8(Emax,NSMAX)
      CALL mtx_broadcast_real8(TLOSS,NSMAX)
      CALL mtx_broadcast_real8(NI_RATIO,NSMAX)

      CALL mtx_broadcast_real8(SPBTOT,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBR0 ,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBRW ,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBENG,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBANG,NBEAMMAX)
      CALL mtx_broadcast_real8(SPBPANG,NBEAMMAX)

      CALL mtx_broadcast_real8(FACT_NRAY,NRAYM)

      cdata (1)=CEWR
      cdata (2)=CEWTH
      cdata (3)=CEWPH
      CALL mtx_broadcast_complex8(cdata,3)
      CEWR =cdata(1)
      CEWTH=cdata(2)
      CEWPH=cdata(3)

      RETURN
  END SUBROUTINE fp_broadcast

!     ***********************
!          PARAMETER VIEW
!     ***********************

  SUBROUTINE fp_view

      use fpcomm_parm
      IMPLICIT NONE
      integer:: nsa,nsb,ns,NBEAM

      WRITE(6,600) 'E0      ',E0      ,'ZEFF    ',ZEFF
      IF(NRMAX.EQ.1) THEN
         WRITE(6,600) &
              'R1      ',R1      ,'DELR1   ',DELR1   ,'DELT    ',DELT 
      ELSE
         WRITE(6,600) &
              'RMIN    ',RMIN    ,'RMAX    ',RMAX    ,'DELT    ',DELT 
      ENDIF

      WRITE(*,*) "----- ANALYZED SPECIES ----- "
      WRITE(*,*) "-----    TEST SPECIES    ----- "
      DO nsa=1,nsamax
         WRITE(6,603) 'ns_nsa  ',ns_nsa(nsa)
      END DO
      WRITE(*,*) "-----    BACKGROUND SPECIES    ----- "
      DO nsb=1,nsbmax
         WRITE(6,603) 'ns_nsb  ',ns_nsb(nsb)
      END DO
      WRITE(*,*) "----- Maximum normalized momentum and loss time ----- "
      DO nsb=1,nsbmax
         WRITE(6,600) 'pmax    ',pmax(nsb), 'pmax_bb ',pmax_bb(nsb), &
                      'tloss   ',tloss(nsb)
      END DO

      WRITE(*,*) "----- PARAMETERS OF HEATINGS -----"
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         WRITE(6,'(A,I2,A,I2)') 'NSA = ',NSA,'  NS  = ',NS
         
         IF(MODELW(NS).EQ.0) THEN
            
            WRITE(6,600) 'PABS_LH ',PABS_LH ,'PABS_FW ',PABS_FW , &
                         'DABS_EC ',PABS_EC
            WRITE(6,600) 'DEC     ',DEC     ,'DELYEC  ',DELYEC
            WRITE(6,600) 'PEC1    ',PEC1    ,'PEC2    ',PEC2
            WRITE(6,600) 'PEC3    ',PEC3    ,'PEC4    ',PEC4
            WRITE(6,600) 'DLH     ',DLH     ,'RLH     ',RLH
            WRITE(6,600) 'PLH1    ',PLH1    ,'PLH2    ',PLH2
            WRITE(6,600) 'DFW     ',DFW     ,'RFW     ',RFW
            WRITE(6,600) 'PFW1    ',PFW1    ,'PFW2    ',PFW2
            
         ELSEIF(MODELW(NS).EQ.1) THEN
            WRITE(6,600) 'PABS_WR ',PABS_WR ,'DELNPR_R',DELNPR_WR, &
                         'DELY_WR ',DELY_WR
            WRITE(6,602) 'EPS_WR  ',EPS_WR  ,'LMAX_WR ',LMAX_WR
            WRITE(6,603) 'NRAYS_WR',NRAYS_WR,'NRAYE_WR',NRAYE_WR

         ELSEIF(MODELW(NS).EQ.2) THEN
            WRITE(6,600) 'PABS_WR ',PABS_WR ,'DELNPR_R',DELNPR_WR, &
                         'DELY_WR ',DELY_WR
            WRITE(6,601) 'FACT_WR ',FACT_WR ,'EPS_WR  ',EPS_WR   , &
                         'LMAX_WR ',LMAX_WR
            
         ELSEIF(MODELW(NS).EQ.3) THEN
            WRITE(6,600) 'PABS_WM ',PABS_WM ,'RF_WM   ',RF_WM, &
                         'DELNPR_M',DELNPR_WM
            WRITE(6,600) 'Y0_WM   ',Y0_WM      ,'DELY_WM ',DELY_WM
            WRITE(6,600) 'CEWR/R  ',DBLE(CEWR) ,'CEWR/I  ',DIMAG(CEWR)
            WRITE(6,600) 'CEWTH/R ',DBLE(CEWTH),'CEWTH/I ',DIMAG(CEWTH)
            WRITE(6,600) 'CEWPH/R ',DBLE(CEWPH),'CEWPH/I ',DIMAG(CEWPH)
            WRITE(6,600) 'RKWR    ',RKWR       ,'RKWTH   ',RKWTH   , &
                         'RKWPH   ',RKWPH
            
         ELSEIF(MODELW(NS).EQ.4) THEN
            WRITE(6,600) 'PABS_WM  ',PABS_WM   ,'DELNPR_M',DELNPR_WM, &
                         'FACT_WM  ',FACT_WM
         ENDIF
         
         IF(TLOSS(NS).NE.0.D0) THEN
            WRITE(6,600) 'TLOSS   ',TLOSS(NS)
         ENDIF

         IF(MODELW(NS).EQ.0)THEN
            WRITE(6,*) 'GIVEN WAVE DIFFUSION COEFFICIENTS'
         ELSE IF(MODELW(NS).EQ.1)THEN
            WRITE(6,*) 'RAY TRACING WAVE DATA'
         ELSE IF(MODELW(NS).EQ.2)THEN
            WRITE(6,*) 'BEAM TRACING WAVE DATA'
         ELSE IF(MODELW(NS).EQ.3)THEN
            WRITE(6,*) 'GIVEN WAVE AMPLITUDE'
         ELSE IF(MODELW(NS).EQ.4)THEN
            WRITE(6,*) 'FULL WAVE DATA'
         ELSE
            WRITE(6,*) 'XX UNKNOWN MODELW: MODELW =',MODELW(NSA)
         END IF
      
      END DO

      DO NBEAM=1,NBEAMMAX
         IF(SPBTOT(NBEAM).NE.0.D0) THEN
            WRITE(6,*) "NBI NUMBER NBEAM & NSSPB=",NBEAM,NSSPB(NBEAM)
            WRITE(6,600) 'SPBTOT   ',SPBTOT(NBEAM), &
                         'SPBENG   ',SPBENG(NBEAM), &
                         'SPBANG   ',SPBANG(NBEAM)
            WRITE(6,600) 'SPBPANG   ',SPBPANG(NBEAM), &
                         'SPBR0    ',SPBR0(NBEAM), &
                         'SPBRW    ',SPBRW(NBEAM)
         ENDIF
      END DO

      IF(MODELS.eq.0)THEN
         WRITE(6,*) 'No fusion reaction'
      ELSEIF(MODELS.eq.1)THEN
         IF(SPFTOT.NE.0.D0) THEN
            WRITE(6,*) 'Given profile fusion reaction'
            WRITE(6,*) "NSSPF=",NSSPF
            WRITE(6,600) 'SPFTOT   ',SPFTOT, &
                         'SPFENG   ',SPFENG
            WRITE(6,600) 'SPFR0    ',SPFR0, &
                         'SPFRW    ',SPFRW
         ENDIF
      ELSE
         WRITE(6,*) 'Self-consistent fusion reaction'
      END IF

      WRITE(6,604)'DRR0            ',DRR0            , &
                  'DRRS            ',DRRS
      WRITE(6,604)'FACTOR_CDBM     ',FACTOR_CDBM     , &
                  'DELTAB_B        ',DELTAB_B
      WRITE(6,604)'DRR_EDGE        ',DRR_EDGE        , &
                  'RHO_EDGE        ',RHO_EDGE
      WRITE(6,604)'FACTOR_DRR_EDGE ',FACTOR_DRR_EDGE , &
                  'FACTOR_PINCH    ',FACTOR_PINCH


      WRITE(6,*) "----- OTHER PARAMETERS -----"
      WRITE(6,600) 'RIMPL   ',RIMPL   ,'EPSM    ',EPSM    ,'EPSFP   ',EPSFP
      WRITE(6,600) 'EPSE    ',EPSE    ,'EPSDE   ',EPSDE   ,'H0DE    ',H0DE    
      WRITE(6,600) 'PGMAX   ',PGMAX   ,'RGMAX   ',RGMAX   ,'RGMIN   ',RGMIN
      WRITE(6,603) 'LMAXE   ',LMAXE   ,'LLMAX   ',LLMAX   ,'NGLINE  ',NGLINE
      WRITE(6,603) 'IDBGFP  ',IDBGFP  ,'NGRAPH  ',NGRAPH  ,'LLMAX_NF',LLMAX_NF

      WRITE(6,603) 'NPMAX   ',NPMAX   ,'NTHMAX  ',NTHMAX  ,'NRMAX   ',NRMAX
      WRITE(6,603) 'NAVMAX  ',NAVMAX  ,'NP2MAX  ',NP2MAX  ,'NTMAX   ',NTMAX   
      WRITE(6,603) 'NTG1STEP',NTG1STEP,'NTG1MIN ',NTG1MIN ,'NTG1MAX ',NTG1MAX
      WRITE(6,603) 'NTG2STEP',NTG2STEP,'NTG2MIN ',NTG2MIN ,'NTG2MAX ',NTG2MAX
      WRITE(6,606) 'NTSTEP_COEF     ',NTSTEP_COEF, &
                   'NTSTEP_COLL     ',NTSTEP_COLL
      WRITE(6,603) 'MODELE  ',MODELE  ,'MODELA  ',MODELA
      WRITE(6,603) 'MODELR  ',MODELR  ,'MODELS  ',MODELS  ,'MODELD  ',MODELD
      WRITE(6,606) 'MODELD_RDEP     ',MODELD_RDEP    , &
                   'MODELD_PDEP     ',MODELD_PDEP
      WRITE(6,606) 'MODELD_EDGE     ',MODELD_EDGE    , &
                   'MODELD_PINCH    ',MODELD_PINCH
      WRITE(6,606) 'MODELD_BOUNDARY ',MODELD_BOUNDARY, &
                   'MODELD_CDBM     ',MODELD_CDBM
      WRITE(6,606) 'MODEL_LOSS      ',MODEL_LOSS     , &
                   'MODEL_SYNCH     ',MODEL_SYNCH
      WRITE(6,606) 'MODEL_NBI       ',MODEL_NBI      , &
                   'MODEL_WAVE      ',MODEL_WAVE
      WRITE(6,603) 'IMTX    ',IMTX    , &
                   'LMAXFP  ',LMAXFP  , &
                   'LMAXE   ',LMAXE
      WRITE(6,606) 'MODEL_KSP       ',MODEL_KSP      , &
                   'MODEL_PC        ',MODEL_PC
      WRITE(6,603) 'NGLINE  ',NGLINE  , &
                   'NGRAPH  ',NGRAPH
      WRITE(6,603) 'LLMAX   ',LLMAX   , &
                   'LLMAX_NF',LLMAX_NF, &
                   'IDBGFP  ',IDBGFP
      WRITE(6,606) 'MODEL_DISRUPT   ',MODEL_DISRUPT   , &
                   'MODEL_Connor_fp ',MODEL_Connor_fp
      WRITE(6,606) 'MODEL_BS        ',MODEL_BS        , &
                   'MODEL_jfp       ',MODEL_jfp
      WRITE(6,606) 'MODEL_LNL       ',MODEL_LNL       , &
                   'MODEL_RE_pmax   ',MODEL_RE_pmax
      WRITE(6,606) 'MODELD_n_RE     ',MODELD_n_RE     , &
                   'MODEL_IMPURITY  ',MODEL_IMPURITY
      WRITE(6,606) 'MODEL_LNL       ',MODEL_LNL       , &
                   'MODEL_RE_pmax   ',MODEL_RE_pmax
      WRITE(6,606) 'MODEL_SINK      ',MODEL_SINK      , &
                   'N_IMPU          ',N_IMPU
      WRITE(6,604) 'T0_quench       ',T0_quench       , &
                   'tau_quench      ',tau_quench
      WRITE(6,604) 'tau_mgi         ',tau_mgi         , &
                   'time_quench_star',time_quench_start
      WRITE(6,600) 'RJPROF1 ',RJPROF1 , &
                   'RJPROF2 ',RJPROF2 , &
                   'v_RE    ',v_RE
      WRITE(6,604) 'target_zeff     ',target_zeff     , &
                   'SPITOT          ',SPITOT

      WRITE(6,*) "-------- PLASMA MODELS --------"

      IF(MODELE.EQ.0)THEN
         WRITE(6,*) 'FIXED ELECTRIC FIELD'
      ELSE IF(MODELE.EQ.1)THEN
         WRITE(6,*) 'CONSISTENT ELECTRIC FIELD'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELE: MODELE =',MODELE
      ENDIF

      IF(MODELR.EQ.0)THEN
         WRITE(6,*) 'NONRELATIVISTIC'
      ELSE IF(MODELR.EQ.1)THEN
         WRITE(6,*) 'RELATIVISTIC'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELR: MODELR =',MODELR
      ENDIF

      IF(MODELD.EQ.0)THEN
         WRITE(6,*) 'WITHOUT RADIAL TRANPORT'
      ELSE IF(MODELD.EQ.1)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (const. for r,p,th, without pinch)'
      ELSE IF(MODELD.EQ.2)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (const. for r,p,th, with pinch)'
      ELSE IF(MODELD.EQ.3)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (p dependence without pinch)'
      ELSE IF(MODELD.EQ.4)THEN
         WRITE(6,*) 'WITH RADIAL TRANSPORT (p dependence with pinch)'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELD: MODELD =',MODELD
      ENDIF

      IF(MODELA.EQ.0)THEN
         WRITE(6,*) 'NOT BOUNCE AVERAGED'
      ELSE IF(MODELA.EQ.1)THEN
         WRITE(6,*) 'BOUNCE AVERAGED'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELA: MODELA =',MODELA
      ENDIF

      DO NSB=1,NSBMAX
         NS=NS_NSB(NSB)
         IF(MODELC(NS).EQ.0)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': LINEAR COLLISION OPERATOR & CONST. T'
         ELSE IF(MODELC(NS).eq.1)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': LINEAR COLLISION OPERATOR & VARIABLE. T'
         ELSE IF(MODELC(NS).EQ.2)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': NONLINEAR COLLISION OPERATOR FOR LIKE PARTILCES & CONST. T'
         ELSE IF(MODELC(NS).eq.3)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': NONLINEAR COLLISION OPERATOR FOR LIKE PARTILCES & VAR. T'
         ELSE IF(MODELC(NS).eq.4)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': NONLINEAR COLLISION OPERATOR'
         ELSE IF(MODELC(NS).EQ.5)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': NONLINEAR COLLISION OPERATOR'
         ELSE IF(MODELC(NS).EQ.-1)THEN
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': LINEAR COLLISION OPERATOR WITH ION SCATTERING'
         ELSE
            WRITE(6,'(A,I5,A)') 'NSB=',NSB, &
                 ': XX UNKNOWN MODELC: MODELC =',MODELC(NS)
         END IF
      END DO

      IF(MODELG.EQ.2)THEN
         WRITE(6,*) 'GIVEN PLASMA GEOMETRY'
      ELSE IF(MODELG.EQ.3)THEN
         WRITE(6,*) 'MHD EQUILIBRIUM FROM TASK/EQ'
      ELSE IF(MODELG.EQ.5)THEN
         WRITE(6,*) 'MHD EQUILIBRIUM FROM EFIT'
      ELSE IF(MODELG.EQ.8)THEN
         WRITE(6,*) 'MHD EQUILIBRIUM FROM TOPICS/EQU'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELG: MODELG =',MODELG
      END IF

      WRITE(6,*) "-------- MPI CONFIGURATION --------"
      WRITE(6,'(A,I4)') "THE NUMBER MPI PROCESS   =", nsize
      WRITE(6,'(A,I4)') "PARTITION NUMBER FOR NSA =", N_partition_s
      WRITE(6,'(A,I4)') "PARTITION NUMBER FOR NR  =", N_partition_r
      WRITE(6,'(A,I4)') "PARTITION NUMBER FOR NP  =", N_partition_p

      RETURN

  600 FORMAT(' ',A8,'=',1PE12.4:3X,A8,'=',1PE12.4:3X,A8,'=',1PE12.4)
  601 FORMAT(' ',A8,'=',1PE12.4:3X,A8,'=',1PE12.4:3X,A8,'=',I8)
  602 FORMAT(' ',A8,'=',1PE12.4:3X,A8,'=',I8,4X  :3X,A8,'=',I8)
  603 FORMAT(' ',A8,'=',I8,4X  :3X,A8,'=',I8,4X  :3X,A8,'=',I8)
  604 FORMAT(' ',A16,'=',1PE12.4:3X,A16,'=',1PE12.4)
  605 FORMAT(' ',A16,'=',1PE12.4:3X,A16,'=',I8)
  606 FORMAT(' ',A16,'=',I8,4X  :3X,A16,'=',I8)
  END SUBROUTINE fp_view

END module fpparm
