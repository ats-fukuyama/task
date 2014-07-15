!
!	t3d_parasurvey_loop.f90
!	task3D_modified
!
!	Created by WAKASA Arimitsu on 11/08/28.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!

subroutine PARAMETER_SURVEY_LOOP
      use trunit
      USE TRCOMM, ONLY : &
           & MDLUF, MDLXP, NT, NTMAX, NTMAX_SAVE, ALLOCATE_TRCOMM,NRMAX,&
           & CK0,CK1,MDLKAI,KUFDCG,BB,TIME_INT,RR,ra
      use T3D_FILE_IO
      use task3d_tentative_flags   
      use t3d_er_param
      
      implicit none
      EXTERNAL TRPARM
      INTEGER(4)       :: IERR
            
      integer(4):: cnt_Cloop
      integer(4):: cnt_gb_mu_loop
      integer(4)::  cnt_SHOT_loop

      integer(4),parameter :: cnt_Cloop_MAX=1
      integer(4),parameter :: cnt_gb_mu_loop_MAX=1
      real(8),dimension(1:cnt_Cloop_MAX) :: CK_loop=0.d0,CK0_loop=0.d0,CK1_loop=0.d0
      real(8),dimension(1:cnt_Cloop_MAX) :: CKe0_loop=0.d0,CKe1_5_loop=0.d0,CKi0_loop=0.d0,CKi1_5_loop=0.d0

      real(8),dimension(1:cnt_gb_mu_loop_MAX) :: gb_mu_loop
      real(8),dimension(1:cnt_gb_mu_loop_MAX) :: gb_mu_e_loop
      real(8),dimension(1:cnt_gb_mu_loop_MAX) :: gb_mu_i_loop
      integer(4)::cnt_MAX_SHOT_loop=0
      real(8) :: CKe0_tmp
      
      integer(4):: cnt_CKe0,cnt_CKe1_5,cnt_CKi1_5,cnt_CKi0,cnt_collective_CK=0,cnt_collective_CK_loop=0

     write(6,*) "################################################################"
     write(6,*) "##                                                            ##"   
     write(6,*) "##                                                            ##"   
     write(6,*) "##   ----------------   SHOT LOOP START    ----------------   ##"
     write(6,*) "##                                                            ##"   
     write(6,*) "##                                                            ##"   
     write(6,*) "################################################################"
     write(6,*)
         
!    ##################################################
!    ########                                  ########
!    ########    PARAMETER 4 TRPARM SETTING    ########
!    ########                                  ########
!    ##################################################
#define ONE_SHOT_LOOP_gBmodel_C_CONST
!#define TWO_SHOT_LOOP_gBmodel_C_CONST
!#define THREE_SHOT_LOOP_gBmodel_C_CONST
!#define C_OPTIMIZATION_LOOP_C_CONST
!#define ALLSHOT_LOOP_C_CONST
!#define C_OPTIMIZATION_LOOP
!#define C_OPTIMIZATION_LOOP_C_CONST
!#define POWERBALANCE_SHOTLOOP
!  #
!  # 1: HIGHTI_EXP_C_CONST
!  # 2: ALLSHOT_LOOP_C_CONST
!  # 3: C_OPTIMIZATION_LOOP
!  # 4: ALLSHOT_LOOP
!  #
!本ルーチンにおけるShotループはあくまで簡易版である．
!すなわち，traparmの再読み込みを省略し，一部のデータのみを本ループ内で上書きすることで代用している．
!基本的には最初に読み込んだtrparmのデータを使用していることに注意せよ．
#ifdef  POWERBALANCE_SHOTLOOP
        cnt_MAX_SHOT_loop=14
#endif
#ifdef  ONE_SHOT_LOOP_gBmodel_C_CONST
        cnt_MAX_SHOT_loop=1
#endif
#ifdef  TWO_SHOT_LOOP_gBmodel_C_CONST
        cnt_MAX_SHOT_loop=2
#endif
#ifdef  THREE_SHOT_LOOP_gBmodel_C_CONST
        cnt_MAX_SHOT_loop=3
#endif
#ifdef  C_OPTIMIZATION_LOOP
        cnt_MAX_SHOT_loop=3
#endif
#ifdef  C_OPTIMIZATION_LOOP_C_CONST
        cnt_MAX_SHOT_loop=3
#endif
#ifdef HIGHTI_EXP_C_CONST
        cnt_MAX_SHOT_loop=8
#ENDIF
#IFDEF ALLSHOT_LOOP_C_CONST
        cnt_MAX_SHOT_loop=14
#endif
#ifdef OLDCHOSEN_SHOT         
        cnt_MAX_SHOT_loop=15
#endif         

    do cnt_SHOT_loop=1,cnt_MAX_SHOT_loop

#ifdef ONE_SHOT_LOOP_gBmodel_C_CONST
         select case(cnt_SHOT_loop)
         case(1)
!        ========================
         KUFDCG='109697'
         TIME_INT=4.54d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=0.0d0
         CK1=0.0d0
!        ========================
!         KUFDCG='109695'
!         TIME_INT=3.74d0
!         RR=3.6d0
!         ra=0.554d0
!         BB=2.75d0
!         MDLKAI=219
!         CK0=1.0d0
!         CK1=1.0d0         
!        ========================
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select
        if(MDLKAI==209 .or. MDLKAI==218 .or. MDLKAI==219)then
            CKe0    = 13.8d0
            CKe1_5  = 2.68d0
            CKi0    = 0.d0
            CKi1_5  = 7.07
            gb_mu_e =1.5d0
            gb_mu_i =1.5d0
        endif                     
#endif
#ifdef TWO_SHOT_LOOP_gBmodel_C_CONST
         select case(cnt_SHOT_loop)
         case(1)
         KUFDCG='109125'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=219
         CK0=1.0d0
         CK1=1.0d0
         case(2)
         KUFDCG='109133'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=219
         CK0=1.0d0
         CK1=1.0d0             
!        ========================
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select
        if(MDLKAI==209 .or. MDLKAI==218 .or. MDLKAI==219)then
            CKe0    = 13.8d0
            CKe1_5  = 2.68d0
            CKi0    = 0.d0
            CKi1_5  = 7.07
            gb_mu_e =1.5d0
            gb_mu_i =1.5d0
        endif                     
#endif
#ifdef THREE_SHOT_LOOP_gBmodel_C_CONST
         select case(cnt_SHOT_loop)
         case(1)
         KUFDCG='109082'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=219
         CK0=1.0d0
         CK1=1.0d0
         case(2)
         KUFDCG='109129'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=219
         CK0=1.0d0
         CK1=1.0d0
         case(3)
         KUFDCG='109135'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=219
         CK0=1.0d0
         CK1=1.0d0                     
!        ========================
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select
        if(MDLKAI==209 .or. MDLKAI==218 .or. MDLKAI==219)then
            CKe0    = 13.8d0
            CKe1_5  = 2.68d0
            CKi0    = 0.d0
            CKi1_5  = 7.07
            gb_mu_e =1.5d0
            gb_mu_i =1.5d0
        endif                     
#endif
#ifdef  C_OPTIMIZATION_LOOP
         select case(cnt_SHOT_loop)
         case(1)
         KUFDCG='109125'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(2)
         KUFDCG='109134'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(3)
         KUFDCG='109695'
         TIME_INT=3.74d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select          
#endif
#ifdef  C_OPTIMIZATION_LOOP_C_CONST
!       minimization @ 109125
        Collective_CKe0(1)=13.68738107d0
        Collective_CKe1_5(1)=2.682695795d0
        Collective_CKi0(1)=0d0
        Collective_CKi1_5(1)=7.429971446d0

         select case(cnt_SHOT_loop)
         case(1)
         KUFDCG='109125'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(2)
         KUFDCG='109134'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(3)
         KUFDCG='109695'
         TIME_INT=3.74d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select          
#endif
#ifdef HIGHTI_EXP_C_CONST
         select case(cnt_SHOT_loop)
         case(1)
         KUFDCG='109081'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(2)
         KUFDCG='109082'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(3)
         KUFDCG='109125'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(4)
         KUFDCG='109129'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(5)
         KUFDCG='109131'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(6)
         KUFDCG='109133'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(7)
         KUFDCG='109134'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case(8)
         KUFDCG='109135'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
        Collective_CKe0(1)=1.60133d01
        Collective_CKe1_5(1)=2.68270d00
        Collective_CKi0(1)=0.0d0
        Collective_CKi1_5(1)=7.42997d00
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select         
#ENDIF
#IFDEF ALLSHOT_LOOP_C_CONST
!       3SHOT minimize
!        Collective_CKe0(1)=7.19685673d0
!        Collective_CKe1_5(1)=3.72759372d0
!        Collective_CKi0(1)=0d0
!        Collective_CKi1_5(1)=6.095068271d0
!       minimization @ 109125
!        Collective_CKe0(1)=13.68738107d0
!        Collective_CKe1_5(1)=2.682695795d0
!        Collective_CKi0(1)=0d0
!        Collective_CKi1_5(1)=7.429971446d0
!       minimization @ 109134
!        Collective_CKe0(1)=10d0
!        Collective_CKe1_5(1)=3.72759372d0
!        Collective_CKi0(1)=0d0
!        Collective_CKi1_5(1)=6.095068271d0
!       minimization @ 109695
!        Collective_CKe0(1)=2.682695795d0
!        Collective_CKe1_5(1)=3.72759372d0
!        Collective_CKi0(1)=0d0
!        Collective_CKi1_5(1)=7.429971446d0


         select case(cnt_SHOT_loop)
         case(1)
         KUFDCG='109081'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(2)
         KUFDCG='109082'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(3)
         KUFDCG='109125'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(4)
         KUFDCG='109129'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(5)
         KUFDCG='109131'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(6)
         KUFDCG='109133'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(7)
         KUFDCG='109134'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(8)
         KUFDCG='109135'
         TIME_INT=4.24d0
         RR=3.6d0
         ra=0.554d0
         BB=2.85d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case(9)
         KUFDCG='109695'
         TIME_INT=3.74d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         ! s088343
         case(10)
         KUFDCG='088343'
         TIME_INT=1.84d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         ! s011369
         case(11)
         KUFDCG='011369'
         TIME_INT=2.2d0
         RR=3.60d0
         ra=0.554d0
         BB=1.52d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         ! s16727         
         case(12)
         KUFDCG='016727'
         TIME_INT=1.02d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.5d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         ! s32940       
         case(13)
         KUFDCG='032940'
         TIME_INT=2.003d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.52d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         ! s081079
         case(14)
         KUFDCG='081079'
         TIME_INT=0.866d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.5d0
         MDLKAI=209
         CK0=1.0d0
         CK1=1.0d0
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select          
        if(MDLKAI==209 .or. MDLKAI==218 .or. MDLKAI==219)then
            CKe0    = 13.8d0
            CKe1_5  = 2.68d0
            CKi0    = 0.d0
            CKi1_5  = 7.07
            gb_mu_e =1.5d0
            gb_mu_i =1.5d0
        endif         
         
#endif
#ifdef OLDCHOSEN_SHOT
         select case(cnt_SHOT_loop)
         ! s088343
         case(1)
         KUFDCG='088343'
         TIME_INT=1.84d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=4
         CK0=0.251d0
         CK1=0.251d0
         case(2)
         KUFDCG='088343'
         TIME_INT=1.84d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=210
         CK0=5.62d0
         CK1=5.62d0
         case(3)
         KUFDCG='088343'
         TIME_INT=1.84d0
         RR=3.6d0
         ra=0.554d0
         BB=2.75d0
         MDLKAI=200
!         CK0=25d0
!         CK1=25d0
         CK0=21.2
         CK1=21.2
         ! s11369
         case(4)
         KUFDCG='011369'
         TIME_INT=2.2d0
         RR=3.60d0
         ra=0.554d0
         BB=1.52d0
         MDLKAI=4
         CK0=0.251d0
         CK1=0.251d0
         case(5)
         KUFDCG='011369'
         TIME_INT=2.2d0
         RR=3.60d0
         ra=0.554d0
         BB=1.52d0
         MDLKAI=210
         CK0=5.62d0
         CK1=5.62d0
         case(6)
         KUFDCG='011369'
         TIME_INT=2.2d0
         RR=3.60d0
         ra=0.554d0
         BB=1.52d0
         MDLKAI=200
         CK0=21.2d0
         CK1=21.2d0
         ! s16727         
         case(7)
         KUFDCG='016727'
         TIME_INT=1.02d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.5d0
         MDLKAI=4
         CK0=0.251d0
         CK1=0.251d0
         case(8)
         KUFDCG='016727'
         TIME_INT=1.02d0
         RR=3.75d0
         ra=0.5634d0 
         BB=1.5d0
         MDLKAI=210
         CK0=5.62d0
         CK1=5.62d0
!         CK0=1.8d0
!         CK1=1.8d0
         case(9)
         KUFDCG='016727'
         TIME_INT=1.02d0
         RR=3.75d0
         ra=0.5634d0 
         BB=1.5d0
         MDLKAI=200
         CK0=21.2d0
         CK1=21.2d0       
         ! s32940       
         case(10)
         KUFDCG='032940'
         TIME_INT=2.003d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.52d0
         MDLKAI=4
         CK0=0.251d0
         CK1=0.251d0
         case(11)
         KUFDCG='032940'
         TIME_INT=2.003d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.52d0
         MDLKAI=210
         CK0=5.62d0
         CK1=5.62d0
         case(12)
         KUFDCG='032940'
         TIME_INT=2.003d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.52d0
         MDLKAI=200
         CK0=21.2d0
         CK1=21.2d0
!         CK0=25d0
!         CK1=25d0
         ! s81079      
         case(13)
         KUFDCG='081079'
         TIME_INT=0.866d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.5d0
         MDLKAI=4
         CK0=0.251d0
         CK1=0.251d0
         case(14)
         KUFDCG='081079'
         TIME_INT=0.866d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.5d0
         MDLKAI=210
         CK0=5.62d0
         CK1=5.62d0
!         CK0=1.8d0
!         CK1=1.8d0
         case(15)
         KUFDCG='081079'
         TIME_INT=0.866d0
         RR=3.75d0
         ra=0.5634d0
         BB=1.5d0
         MDLKAI=200
!         CK0=25d0
!         CK1=25d0
         CK0=21.2d0
         CK1=21.2d0
         case default
            print*,"Error in the SHOT loop @",__FILE__,__LINE__
            stop
         end select          
#endif         

!    call CKe0CKe1_SETTING
!    do cnt_CKe0=1,11    
!        CKe0=Collective_CKe0(cnt_CKe0)
!    do cnt_CKe1_5=1,11
!        CKe1_5=Collective_CKe1_5(cnt_CKe1_5)

    flg_ERF_Difference_Coptimize=0
    flg_t3d_er_calc_timing =0
    CNT_TREXE_ITERATION =-1

    call T3D_FILENAME_MAKER

    CALL tr_prof(ierr)
    CALL TRLOOP
            

    enddo
!    enddo
!    enddo
end subroutine PARAMETER_SURVEY_LOOP
