!     ***********************************************************

!           MAIN ROUTINE FOR TRANSPORT CALCULATION

!     ***********************************************************

      SUBROUTINE TRLOOP

      USE TRCOMM, ONLY : &
     &   ABRHOG, AJ, AJU, AMM, ANC, ANFE, ANNU, AR1RHOG, AX, AY, AZ, BB, &
     &   DR, DT, DVRHO, DVRHOG, EPSLTR, LDAB, LMAXTR, MDLEQB, MDLEQE, &
     &   MDLEQN, MDLPCK, MDLUF, MDTC, MLM, MODELG, NEQMAX, NFM, NGR, NGRSTP, &
     &   NGT, NGTSTP, NRAMAX, NRM, NRMAX, NROMAX, NSM, NSMAX, NSS, NST, NSV, &
     &   NT, NTEQIT, NTMAX, NTSTEP, NTUM, &
     &   PA, PZ, PZC, PZFE, Q0, QP, RDP, RG, RHOA, RIP, RIPE, RIPS, RIPU, &
     &   RMU0, RN, RR, RT, RU, RW, T, TPRST, TST, TTRHO, TTRHOG, &
     &   VLOOP, VSEC, X, XV, Y, YV, Z, ZV ,NEQMAXM, DIPDT, MDLEOI, &
     &   MDLEQN, MDLEQT,PTS,PT,TIME_INT,PECU,nstm,PEX, &
     &   MDLKAI, CK0, CK1, PI,EPS0,AEE,RTM,AMZ,PIN,ra,AKNC,AKDW,AR2RHO,er,KUFDCG
      USE TRCOM1, ONLY : TMU, TMU1, NTAMAX, NTXMAX, NTXMAX1
      use tr_bpsd,only: tr_bpsd_set, tr_bpsd_get
      use trunit
      use equnit_mod
      use equunit_mod


      use task3d_tentative_flags, only : FIT3D_tentative_flags, T3D_POWERBALANCE_CALC_flg, flg_teti2ufile, flg_teti_comparison, T3D_POWERBALANCE_CALC_SHOTLOOP_flg
      use T3D_FILE_IO      
      
      IMPLICIT NONE
      INTEGER(4):: IERR
      REAL(8)   :: FCTR
      
      real(8) :: cpu_time_start, cpu_time_end

      call cpu_time(cpu_time_start)

      
      PT(1)=RT(1,1)*1.05
      PT(2)=RT(1,2)*1.05
      PTS(1)=RT(NRMAX,1)*0.95
      PTS(2)=RT(NRMAX,2)*0.95

      if(nt /= 0) then 
        print*,'Sorry, continuous calculation routine(NT.eq.1) is imcomplete.'
        print*, 'NT =',nt
        stop
      endif
   

!     温度分布比較のため初回温度分布（実験測定値を想定）及び各輸送係数係数の保存ルーチン
!     パワーバランス解析時には本ルーチンで算出した各種値を使用．
!     変更＠20121101 by Arimitsu Wakasa
!======================================================>>>>@20121101v01
      if(T3D_POWERBALANCE_CALC_flg == 0 .and. flg_teti_comparison == 1) call CAMPARISON_TETI_T0
      CALL TRCALC(IERR)
      CALL TRGLOB

      if(T3D_POWERBALANCE_CALC_flg==1) then
        call POWERBALANCE_CALC
        if (T3D_POWERBALANCE_CALC_SHOTLOOP_flg ==1) return
      endif
      if(flg_OUTPUT_FILE_WP_EVO==1) call t3d_write_file_evo
!======================================================<<<<@20121101v01

      write(6,'(a8,$)') ": nt ="
      
      do nt=0,ntmax

      write(6,'(i3,$)') nt
      
      if(FIT3D_tentative_flags==1) call FIT3D_EXEC

!===============================================================================
!===============================================================================

!      CALL TREVAL(NT,IERR)
!      IF(IERR.NE.0) exit


      CALL trexec(DT,IERR)
      IF(IERR.NE.0) exit

!     /* Sawtooth Oscillation */
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))
      IF(Q0.LT.1.D0) TST=TST+DT

      IF(TST+0.5D0*DT.GT.TPRST) THEN
         print *, '+++ start sub TRSAWT'
         CALL TRSAWT
         print *, '+++ end sub TRSAWT'
         TST=0.D0
      ENDIF

!     *** SET GEOMETRY VIA TASK/EQ ***

      IF(NTEQIT.NE.0) THEN
         IF(MOD(NT,NTEQIT).EQ.0) THEN
            if(modelg.eq.8) THEN
               call equ_calc
               call tr_bpsd_get(IERR)
               if(ierr.ne.0) return
            endif
            if(modelg.eq.9) THEN
               call eq_calc
               call tr_bpsd_get(IERR)
               if(ierr.ne.0) return
            endif
         ENDIF
      ENDIF

      if(mod(nt-1,ngtstp)==0)then
        if(flg_OUTPUT_FILE_WP_EVO==1) call t3d_write_file_evo
      endif

!     *** READING DATA FROM UFILES FOR NEXT STEP ***
      
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) CALL TR_UFREAD
!      IF(MDLUF.EQ.2.AND.MODEP.EQ.3) CALL TR_UFREAD_S
      IF(MDLUF.EQ.2) then
        CALL TR_UFREAD_S
      endif
    enddo
    

      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         RIPS=RIP
         RIPE=RIP
      ELSE
         RIPS=RIPE
      ENDIF
 

!     各種ファイル出力処理ルーチン
!     出力ファイル生成の有無は, task3d_tentative_flags.f90 内の各フラグ変数にて決定．
!     1: ON 0: OFF
!           - flg_teti2ufiles     :: 最終温度分布をUFILE形式で出力するか否か．
!           - flg_teti_comparison :: 最終温度分布をと初期温度分との誤差評価を出力するか否か．．
!     出力ファイル名のメモリ解放処理
!     出力ファイル関連フラグ初期化
!     変更＠20121101 by Arimitsu Wakasa
!======================================================>>>>@20121101v02
      if(flg_teti2ufile == 1) call teti2ufiles
      if(T3D_POWERBALANCE_CALC_flg == 0 .and. flg_teti_comparison == 1) then
        call CAMPARISON_TETI_T0
        CAMPARISON_TETI_T0_file_flg=0
      endif
      if(T3D_POWERBALANCE_CALC_FLG == 0 .and. flg_OUTPUT_FILE_TRCOEFLOG_LASTEST==1)then
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_LASTEST)
        TRCOEFLOG_LASTEST_file_flg=0      
      endif
      if(T3D_POWERBALANCE_CALC_FLG == 0 .and. flg_OUTPUT_FILE_TRCOEFLOG_EVO==1)then
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_EVO)
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010)
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025)
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050)
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075)
        deallocate(STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090)
        TRCOEFLOG_EVO_file_flg=0
      endif
      if(flg_OUTPUT_FILE_WP_EVO==1) then
        deallocate(STR_OUTPUT_FILE_WP_EVO)
        FILE_WP_EVO_file_flg=0  
      endif
!======================================================>>>>@20121101v02
      call cpu_time(cpu_time_end)
      print*
      print*,'========================================================'
      print*,'CPU TIME @tr_loop', cpu_time_end - cpu_time_start, 'sec'
      print*
 
 
 
      RETURN
      END SUBROUTINE TRLOOP