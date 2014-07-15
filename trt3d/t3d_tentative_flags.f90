!
!	task3d_tentative_flags.f90
!	エネルギーバランス式＋WOUT対応版@20111205
!
!	Created by WAKASA Arimitsu on 11/12/05.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!
module task3d_tentative_flags

    implicit none
    public

!======================================================================================================================!
!    Files入出力関連
!    File名決定関連変数群
!======================================================================================================================!
!    変数名：flg_OUTPUT_FILE_TRCOEFLOG_LASTEST
!    変数名：flg_OUTPUT_FILE_TRCOEFLOG_EVO
!    変数名：flg_OUTPUT_FILE_ERF_TT0
!    変数名：flg_OUTPUT_FILE_WP_EVO
!    タイプ：整数型
!    上記各出力の有無の決定．
!    1: ON 0: OFF

     integer(4), parameter :: flg_OUTPUT_FILE_TRCOEFLOG_LASTEST=1
     integer(4), parameter :: flg_OUTPUT_FILE_TRCOEFLOG_EVO=1
     integer(4), parameter :: flg_OUTPUT_FILE_ERF_TT0=1
     integer(4), parameter :: flg_OUTPUT_FILE_WP_EVO=1
     
!    変数名：flg_OUTPUT_FILE_NAME_LHDSHOT
!    変数名：flg_OUTPUT_FILE_NAME_SHOTTIME
!    変数名：flg_OUTPUT_FILE_NAME_TBMODEL
!    変数名：flg_OUTPUT_FILE_NAME_CKeCKi
!    変数名：flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i
!    変数名：flg_OUTPUT_FILE_NAME_muemui
!    タイプ：整数型
!    上記各出力ファイルに対してファイル名に記載する情報の有無の決定．
!    1: ON 0: OFF

     integer(4), parameter :: flg_OUTPUT_FILE_NAME_LHDSHOT=1
     integer(4), parameter :: flg_OUTPUT_FILE_NAME_SHOTTIME=1
     integer(4), parameter :: flg_OUTPUT_FILE_NAME_TBMODEL=1
     integer(4), parameter :: flg_OUTPUT_FILE_NAME_CKeCKi=0
     integer(4), parameter :: flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i=1
     integer(4), parameter :: flg_OUTPUT_FILE_NAME_muemui=1

!======================================================================================================================!   
!   特殊ループ制御
!   現在未完成：ブラッシュアップが必要です．
!   shot loop & TB_C_const servey loop
!======================================================================================================================!
!   変数名 : flg_SHOT_LOOP_CONTROL
!   変数名 : flg_TB_C_LOOP_CONTROL
!    タイプ：整数型
!    上記各出力ファイルに対してファイル名に記載する情報の有無の決定
!    1: ON 0: OFF

        integer(4), parameter :: flg_SHOT_LOOP_CONTROL=0
        integer(4), parameter :: flg_TB_C_LOOP_CONTROL=0

!======================================================================================================================!   
!   変数名：PIN_EQ
!   タイプ：実数型配列
!   タグ：加熱, エネルギー等分配
!   加熱項詳細保存用配列（エネルギー等分配項）
!   注意：現状配列要素数は一時的なものです．配列破壊に注意してください．

    real(8), dimension(10,10,500) :: PIN_EQ=0.0d0

!   変数名：FIT3D_ONOFF_flg
!   タイプ：整数型
!   タグ：加熱, FIT3D
!   FIT3Dによる加熱パワー分布計算を行うか否か
!   1: ON 0: OFF
    integer, parameter :: FIT3D_tentative_flags=0

    !   変数名：FIT3DtoUFILE_flg
    !   タイプ：整数型
    !   タブ：加熱，FIT3D, UFILE
    !   FIT3Dでの計算結果を定常解析用にUFILE仕様で書きだすか否か
    !   1: ON 0: OFF
        integer:: fit3dtoufile_flg=0

        !   FIT3DtoUFILE_flg から派生．
        !   変数名：nTPolynomialtoUFILE_flg
        !   タイプ：整数型
        !   タブ：加熱，FIT3D, UFILE, 密度, 温度, 多項式
        !   「FIT3Dでの計算結果を定常解析用にUFILE仕様で書きだす際に使用する温度密度分布」を多項式形式で読み込んでUfile化するかどうか．
        !   1: ON 0: OFF
            integer:: ntPolynomial2UFILE_flg=0

    
!   変数名：FIT3D_CALC_TIMING
!   タイプ：整数型配列
!   タグ：加熱, FIT3D, FILEIO
!   FIT3Dでの計算をどのタイミング（nt=??）で計算するかを決定（上限1000とする．）
!   初期値を-1として，-1が現出以降にはFIT3Dの計算はさせない．
    integer,save,dimension(1000) :: FIT3D_CALC_TIMING=-1
    
!======================================================================================================================!
!   変数名：T3D_POWERBALANCE_CALC_flg
!   タイプ：整数型
!   タグ：パワーバランス, 定常
!   Ufile等で読み込んだ温度密度分布を元に，必要ならば加熱パワー吸収分布を計算の後，
!   定常状態（パワーバランス）を想定して，各熱輸送フラックス，熱輸送係数の算出を行う．
!   フラグ有効時には，TRでの時間発展計算は行われない．
!   1:ON 0:OFF
    integer,parameter :: T3D_POWERBALANCE_CALC_flg=0

    !   変数名：T3D_POWERBALANCE_CALC_SHOTLOOP_flg
    !   タイプ：整数型
    !   タグ：パワーバランス, 定常, ショットループ
    !       - Ufile等で読み込んだ温度密度分布を元に，必要ならば加熱パワー吸収分布を計算の後，
    !       - 定常状態（パワーバランス）を想定して，各熱輸送フラックス，熱輸送係数の算出を行う．
    !       - フラグ有効時には，TRでの時間発展計算は行われない．
    !   以上の動作後，他のショットのループを行うためtrmenuへ回帰するかどうかの判定フラグ．
    !   1:ON 0:OFF
        integer,parameter :: T3D_POWERBALANCE_CALC_SHOTLOOP_flg=0

!======================================================================================================================!
!   変数名：need_boz_Translate_flg
!   タイプ：整数型
!   タグ：NEWBOZ, 座標系
!   BOZデータの生成時に右手系・左手系など物理量の正負が異なる場合があったので，
!   それなの整合性を取るための変換ルーチンを作成した．
!   フラグ有効時には，BOZデータ読み込み後，左手系から右手系などの変換を行う．
!   2012.12.14以降，変換の必要のないBOZデータが作成されているはずである，
!   1:ON 0:OFF
    integer, parameter ::need_boz_Translate_flg=0

!======================================================================================================================!    
!   変数名：flg_teti2ufile
!   タイプ：整数型
!   タグ：UFILE
!   計算結果の温度分布Te&TiをUfile形式で書きだす．               
!   FIT３Dへのデータ受け渡し用に制作                            
!   現在(2011.12.06以降）はUFILEを介さずFIT3Dを実行しているので未使用）
!   1:ON 0:OFF    
    integer, parameter :: flg_teti2ufile=0

!======================================================================================================================!
!   変数名：flg_teti_comparison
!   タイプ：整数型
!   タグ：温度, 実験
!   Initialの温度分布と最終温度分布との誤差評価を行い結果を出力する．          
!   主に実験データ（初期条件）との差異を評価（輸送モデルやコンスタントファクタ等）するために使用．
!   1:ON 0:OFF    
    integer, parameter :: flg_teti_comparison=1

!======================================================================================================================!
!   CKe, CKi                :: 電子イオン其々について複数個のコンスタントファクタを有するする際に使用する変数
!   gb_mu, gb_mu_e, gb_mu_i :: 温度勾配項のべき乗を決定する変数
!   gB+gradTモデル（Case (209)）等で使用．
!   タイプ：実数型
!   タグ：GBモデル, 乱流輸送
!======================================================================================================================!
    real(8),save :: CKe0, CKe1_5, CKi0, CKi1_5
    real(8),save :: gb_mu,gb_mu_e,gb_mu_i
    real(8),save,dimension(100) :: Collective_CKe0,Collective_CKe1_5  
!======================================================================================================================!


    contains
    
    subroutine FIT3D_CALC_TIMING_initialize
    FIT3D_CALC_TIMING(1:10)=(/0,1,2,5,7,10,15,30,50,100&
&       /)
    end subroutine
    
    subroutine CKe0CKe1_SETTING
    Collective_CKe0(1:11)=(/8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0,16d0,17d0,18d0/)
    Collective_CKe1_5(1:11)=(/0.d0,0.5d0,1.d0,1.5d0,2.d0,2.5d0,3.d0,3.5d0,4.d0,4.5d0,5.d0/)
    end subroutine

end module

