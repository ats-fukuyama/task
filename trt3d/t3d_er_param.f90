!
!	t3d_er_param.f90
!	task.20110214
!
!	Created by WAKASA Arimitsu on 11/02/23.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!
    module t3d_er_param
!-------------------------------------------------------------------------!   
!    ER_CALC_METHOD, is ER calculation method: 10 or 20 or 30
!    
!    10: ER_EXP. 実験結果測定値で固定
!    20: ER_AMBIPOLAR, 両極性条件で決定
!    30: ER_DYNAMICS, 電場の発展方程式を解くことで電場算出
!        31: Euler Explicit,  
!        32: Euler Implicit,
!        33: Crank-Nicolson, 
!        34: balanced
!        35: TEST routine
!
!
!    ER_CALC_CONDITION, is calculation condition: 0 or n or -1 or -2
!
!     0: 初回のみ電場算出，以降，この値で固定し計算
!     n: DTに関して，n step 毎に電場の計算
!    -1: DTに関して，すべてのステップにて電場の計算
!    -2: 指定されたタイミングで電場の計算（指定変数：ER_CALC_TIMING
!
!    DEA_SETTING_CONDITION, is HOW does DEA set: 0 or 1
!
!       0: DEAは径方向に一定．本ファイル（t3d_er_param.f90）で指定されたDEA(1:nrmax)を 使用する．
!       1: DEAはERモジュール内( t3d_er.f90）にて与える．
!
!    ER_CALC_INITIAL_VALUE, is an initial value of Er: 0 or 1 or 2
!        
!        Using for ER_dynamics only.
!        0: constant value defined in t3d_er.f90
!        1: determined by Amipolar condition @ initial nT profiles
!        2: Experimental (or other) data 
!
!
!    ER_CALC_BOUNDARY_CONDITION, is select use or not-use the boundary condition.
!    
!     0: not use the boundary condition.
!     1: use the boundary condition.
!
!
!    CNT_TREXE_ITRATION: 
!
!       trexe.f90 内にて行列解放ループのループ回数保持用変数
!       電場算出ルーチン，S.TRERAD の呼び出し判断に使用する
!       -1: trloop全体を通して初回時
!        0: TRERAD呼び出し有効
!       1〜: TRERAD呼び出し無効（及び，TREXE内Iteration回数保持用）
!
!       flg_CNT_TREXE_ITERATION:
!           CNT_TREXE_ITRATIONに付随するフラグ用変数
!           0:TRITERATION内にてntが同一のままである
!           1:TRITERATION内にてntが同一である最後の周回（ERの計算を行うタイミング(nt)を次へ更新する）
!-------------------------------------------------------------------------!   
!    dt_er: ER＿DYNAMICS適用時に使用．発展方程式を解く際のタイムステップ幅
!    
!    dea: ER＿DYNAMICS適用時に使用．電場の拡散係数
!    
!    er_gamma_factor: ER＿DYNAMICS適用時に使用．発展方程式に適用されるファクタ
!-------------------------------------------------------------------------! 
        implicit none
        public
        
!------------------  PARAMETERS  -----------------------------------------!
        integer(4), parameter :: ER_CALC_METHOD = 20
        integer(4), parameter :: ER_CALC_CONDITION = -2
        integer(4), parameter :: ER_CALC_INITIAL_VALUE = 1
        integer(4), parameter :: ER_CALC_BOUNDARY_CONDITION = 1
        integer(4), parameter :: DEA_SETTING_CONDITION = 0
        
        integer(4), parameter :: TEMPORARY_SETTING = 1
!        integer(4), parameter :: TEMPORARY_SETTING_20110909 = 0
!-------------------------------------------------------------------------!

        integer(4), save ::CNT_TREXE_ITERATION =-1
        integer(4), save ::flg_CNT_TREXE_ITERATION =0
        integer(4),save,dimension(1000) ::ER_CALC_TIMING=-1
        real(8),save:: dt_er
        real(8),save:: er_gamma_factor
        real(8),save:: er_boundary_value
        real(8),allocatable,save:: dea(:)

    contains
        subroutine t3d_er_calc_initialize
            ER_CALC_TIMING(1:51)=(/0,&
&	1	,	2	,	3	,	4	,	5	,	6	,	7	,	8	,	9	,	10	,	&
&	11	,	12	,	13	,	14	,	15	,	16	,	17	,	18	,	19	,	20	,	&
&	21	,	22	,	23	,	24	,	25	,	26	,	27	,	28	,	29	,	30	,	&
&	31	,	32	,	33	,	34	,	35	,	36	,	37	,	38	,	39	,	40	,	&
&	41	,	42	,	43	,	44	,	45	,	46	,	47	,	48	,	49	,	50		&
&       /)

!&	51	,	52	,	53	,	54	,	55	,	56	,	57	,	58	,	59	,	60	,	&
!&	61	,	62	,	63	,	64	,	65	,	66	,	67	,	68	,	69	,	70	,	&
!&	71	,	72	,	73	,	74	,	75	,	76	,	77	,	78	,	79	,	80	,	&
!&	81	,	82	,	83	,	84	,	85	,	86	,	87	,	88	,	89	,	90	,	&
!&	91	,	92	,	93	,	94	,	95	,	96	,	97	,	98	,	99	,	100	,	&
!&	101	,	102	,	103	,	104	,	105	,	106	,	107	,	108	,	109	,	110	,	&
!&	111	,	112	,	113	,	114	,	115	,	116	,	117	,	118	,	119	,	120	,	&
!&	121	,	122	,	123	,	124	,	125	,	126	,	127	,	128	,	129	,	130	,	&
!&	131	,	132	,	133	,	134	,	135	,	136	,	137	,	138	,	139	,	140	,	&
!&	141	,	142	,	143	,	144	,	145	,	146	,	147	,	148	,	149	,	150	,	&
!&	151	,	152	,	153	,	154	,	155	,	156	,	157	,	158	,	159	,	160	,	&
!&	161	,	162	,	163	,	164	,	165	,	166	,	167	,	168	,	169	,	170	,	&
!&	171	,	172	,	173	,	174	,	175	,	176	,	177	,	178	,	179	,	180	,	&
!&	181	,	182	,	183	,	184	,	185	,	186	,	187	,	188	,	189	,	190	,	&
!&	191	,	192	,	193	,	194	,	195	,	196	,	197	,	198	,	199	,	200	,	&
!&	201	,	202	,	203	,	204	,	205	,	206	,	207	,	208	,	209	,	210	,	&
!&	211	,	212	,	213	,	214	,	215	,	216	,	217	,	218	,	219	,	220	,	&
!&	221	,	222	,	223	,	224	,	225	,	226	,	227	,	228	,	229	,	230	,	&
!&	231	,	232	,	233	,	234	,	235	,	236	,	237	,	238	,	239	,	240	,	&
!&	241	,	242	,	243	,	244	,	245	,	246	,	247	,	248	,	249	,	250	    &    
!&       /)
!            ER_CALC_TIMING(1:35)=(/ &
!            & 0,1,2,3,5,7,10,15,20,25, &
!            & 30,40,50,75,100,125,150,175,200,250,&                    
!            & 300,350,400,450,500,550,600,650,700,750,&
!            & 800,850,900,950,1000 &
!            &/)
        end subroutine
    
        subroutine t3d_er_param_initialize
            use trcomm, only:nrmax
            
            if(allocated(dea)) then
                continue
            else
                allocate(dea(nrmax))
            endif
            
            dea=1.0d5
            dt_er = 1.0d-6
            er_gamma_factor= 1.0d-0

            er_boundary_value=13.0d3
            end subroutine

    end module t3d_er_param