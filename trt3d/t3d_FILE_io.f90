!
!	FILE_WRITE_FLGS.f90
!	task3D_modified
!
!	Created by WAKASA Arimitsu on 11/08/28.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!

!module FILE_WRITE_FLGS
module T3D_FILE_IO
use TRCOMM,only:NRMAX,CK0,CK1,MDLKAI,KUFDCG,TIME_INT,PIN,PEX,ER,RT,RN,AVKDW,AVKNC,AVK,AVDW,AVNC,AV,AKDW,AKNC,AK,ADDW,ADNC,AD,RG,T
use task3d_tentative_flags


    implicit none
    integer(4),parameter:: NRMAX4M_FILE_WRITE_FLGS=500
    integer(4),save :: TRCOEF_01

    integer(4),save:: TRCOEFLOG_EVO_file_flg=0
    integer(4),save:: CAMPARISON_TETI_T0_file_flg=0
    integer(4),save:: FILE_WP_EVO_file_flg=0  
    integer(4),save:: TRCOEFLOG_LASTEST_file_flg=0

    real(8),dimension(1:NRMAX4M_FILE_WRITE_FLGS, 0:10),save :: pflx_save, hflx_save, dd_save, vv_save, chi0_save,chi_save, kv_save ! for check the NC parameters

!    character(38) ,save:: FILENAME_CK0CK1mu
!    character(51) ,save:: FILENAME_CK0CK1
!    character(15) ,save:: FILENAME_gbmu
!    character(6),save:: FILENAME_shotnum
    integer(4),save :: flg_t3d_er_calc_timing
    integer(4),save :: flg_ERF_Difference_Coptimize
    real(8),save  ::pre_AKDWIL=0.d0,pre_AKDWEL=0.d0

    real(8),dimension(1:NRMAX4M_FILE_WRITE_FLGS),save :: AKDWe_GB_0=0.d0,AKDWe_GB_1_5=0.d0,AKDWi_GB_0=0.d0,AKDWi_GB_1_5=0.d0
    real(8) :: tmpAKDWe_GB_0=0.d0,tmpAKDWe_GB_1_5=0.d0,tmpAKDWi_GB_0=0.d0,tmpAKDWi_GB_1_5=0.d0
    
!    character(81) :: FILENAME_LOG_LASTEST
!    character(22) :: FILENAME_LOG_LASTEST_endpart

    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_LASTEST
    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_EVO
    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010
    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025
    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050
    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075
    character(len=:),allocatable :: STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090
    character(len=:),allocatable :: STR_OUTPUT_FILE_ERF_TT0
    character(len=:),allocatable :: STR_OUTPUT_FILE_WP_EVO
    character(len=:),allocatable :: STR_OUTPUT_FILE_POWERBALANCE
    
    
    character(7),save :: STR_OUTPUT_FILE_NAME_LHDSHOT="sxxxxxx"
    character(7),save :: STR_OUTPUT_FILE_NAME_SHOTTIME="txx.xxx"
    character(8),save :: STR_OUTPUT_FILE_NAME_TBMODEL="xxxxxxxx"
    character(9), dimension(2), save :: STR_OUTPUT_FILE_NAME_CKeCKi="CKsxx.xxx"
    character(9), dimension(4), save :: STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i="CKsxx.xxx"
    character(7), dimension(2), save  :: STR_OUTPUT_FILE_NAME_mu="musx.xx"

contains

SUBROUTINE T3D_FILENAME_MAKER

!     ************* MAKE LHD SHOT NUMBER *********************    
    if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)then
        if(len_trim(KUFDCG)==6)then
            STR_OUTPUT_FILE_NAME_LHDSHOT="s"//KUFDCG
            print'(a7)',STR_OUTPUT_FILE_NAME_LHDSHOT
        else
            print"(a60)",  "ERROR: The amount of digits of LHDSHOT_NO.(KUFDCG) is not 6."
            print"(a16,a)","       KUFDCG = ",KUFDCG
            print"(a14,a,a16,i)",     "       STOP @ ",__FILE__,"          line =",__LINE__
            stop
        endif
    endif
    
!     ************* MAKE LHD SHOT TIME ************************
    if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)then
        if((TIME_INT/100.d0)>=1.d0)then
            print*,"Sorry, this version is not allow TIME_INIT >= 100."
            print*,"TASK3D is stopped. Check TIME_INIT value."
            stop
        else if ((TIME_INT/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_SHOTTIME,'(a1f6.3)') "t",TIME_INT
        else
            write(STR_OUTPUT_FILE_NAME_SHOTTIME,'(a2f5.3)') "t0",TIME_INT
        endif
        print'(a7)',STR_OUTPUT_FILE_NAME_SHOTTIME
    endif
!     ************* MAKE TB MODEL NAME **************************    
    if(flg_OUTPUT_FILE_NAME_TBMODEL==1)then
        select case(MDLKAI)
        case(4)
            STR_OUTPUT_FILE_NAME_TBMODEL="edgeBohm"
        case(200)
            STR_OUTPUT_FILE_NAME_TBMODEL="gyroBohm"
        case(209)
            STR_OUTPUT_FILE_NAME_TBMODEL="gB+gradT"
        case(210)
            STR_OUTPUT_FILE_NAME_TBMODEL="Alcator_"
        case(218)
            STR_OUTPUT_FILE_NAME_TBMODEL="gB+grdTi"
        case(219)
            STR_OUTPUT_FILE_NAME_TBMODEL="A_gB+grT"
        case default
            print"(a40,i)","this version is not including MDLKAI = ",MDLKAI
            print*,"TASK3D is stopped. Check MDLKAI value."
            print"(a14,a,a16,i)",     "       STOP @ ",__FILE__,"          line =",__LINE__
            stop
        end select
        print'(a8)',STR_OUTPUT_FILE_NAME_TBMODEL
    endif
!     ************* MAKE CKe and CKi **************************
    if(flg_OUTPUT_FILE_NAME_CKeCKi==1)then
        if((CK0/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_CKeCKi(1),'(a3f6.3)') "CKe",CK0
        else
            write(STR_OUTPUT_FILE_NAME_CKeCKi(1),'(a4f5.3)') "CKe0",CK0
        endif

        if((CK1/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_CKeCKi(2),'(a3f6.3)') "CKi",CK1
        else
            write(STR_OUTPUT_FILE_NAME_CKeCKi(2),'(a4f5.3)') "CKi0",CK1
        endif
        print'(a9)',STR_OUTPUT_FILE_NAME_CKeCKi(1)
        print'(a9)',STR_OUTPUT_FILE_NAME_CKeCKi(2)
    endif
!     ***********************************************************    
    ! 現バージョンでは電子イオン両熱拡散係数に関してそれぞれ2つずつのコンスタントファクタ（C1e,C2e,C1i,C2i）を想定している．@20121110
    if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)then
        if((CKe0/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1),'(a3f6.3)') "C1e",CKe0
        else
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1),'(a4f5.3)') "C1e0",CKe0
        endif
        if((CKe1_5/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2),'(a3f6.3)') "C2e",CKe1_5
        else
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2),'(a4f5.3)') "C2e0",CKe1_5
        endif
        if((CKi0/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3),'(a3f6.3)') "C1i",CKi0
        else
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3),'(a4f5.3)') "C1i0",CKi0
        endif
        if((CKi1_5/10.d0)>=1.d0)then
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4),'(a3f6.3)') "C2i",CKi1_5
        else
            write(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4),'(a4f5.3)') "C2i0",CKi1_5
        endif
        print'(a9)',STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1)
        print'(a9)',STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2)
        print'(a9)',STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3)
        print'(a9)',STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4)
    endif
!     ***********************************************************    
    if(flg_OUTPUT_FILE_NAME_muemui==1)then
        write(STR_OUTPUT_FILE_NAME_mu(1),'(a3f4.2)') "mue",gb_mu_e
        write(STR_OUTPUT_FILE_NAME_mu(2),'(a3f4.2)') "mui",gb_mu_i
        print'(a7)',STR_OUTPUT_FILE_NAME_mu(1)
        print'(a7)',STR_OUTPUT_FILE_NAME_mu(2)
    endif
!     ***********************************************************
    return

END SUBROUTINE T3D_FILENAME_MAKER
    
!     ***********************************************************
!           Te & Ti >> Ufile format 
!     ***********************************************************

SUBROUTINE teti2ufiles

    USE TRCOMM, ONLY : nrmax,rt,rg

    integer(4) :: i
    integer(4) :: cnt_nr_forUFILE

    open(761,file='Handover_TE_UFILE.txt',status='replace')
    open(762,file='Handover_TI_UFILE.txt',status='replace')

    cnt_nr_forUFILE=(nrmax+1)/5
    do i=0,cnt_nr_forUFILE-1
        if(i==0) then
            write(761,'(5e14.5e3)') 0.d0,rg(i*5+1),rg(i*5+2),rg(i*5+3),rg(i*5+4)
        else
            write(761,'(5e14.5e3)') rg(i*5+0),rg(i*5+1),rg(i*5+2),rg(i*5+3),rg(i*5+4)
        endif
    enddo
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(761,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0)
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(761,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1)
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(761,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1),rg(cnt_nr_forUFILE*5+2)
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(761,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1),rg(cnt_nr_forUFILE*5+2),rg(cnt_nr_forUFILE*5+3)

    write(761,'(5e14.5e3)') 1.01d0
      
    do i=0,cnt_nr_forUFILE-1
        if(i==0) then
            write(761,'(5e14.5e3)') (rt(1,1)+(rt(1,1)-rt(2,1)))*1d3,rt(i*5+1,1)*1d3,rt(i*5+2,1)*1d3,rt(i*5+3,1)*1d3,rt(i*5+4,1)*1d3
        else
            write(761,'(5e14.5e3)') rt(i*5+0,1)*1d3,rt(i*5+1,1)*1d3,rt(i*5+2,1)*1d3,rt(i*5+3,1)*1d3,rt(i*5+4,1)*1d3
        endif      
    enddo
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(761,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,1)*1d3
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(761,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,1)*1d3,rt(cnt_nr_forUFILE*5+1,1)*1d3
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(761,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,1)*1d3,rt(cnt_nr_forUFILE*5+1,1)*1d3,rt(cnt_nr_forUFILE*5+2,1)*1d3
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(761,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,1)*1d3,rt(cnt_nr_forUFILE*5+1,1)*1d3,rt(cnt_nr_forUFILE*5+2,1)*1d3,rt(cnt_nr_forUFILE*5+3,1)*1d3
                  
    do i=0,cnt_nr_forUFILE-1
        if(i==0) then
            write(762,'(5e14.5e3)') 0.d0,rg(i*5+1),rg(i*5+2),rg(i*5+3),rg(i*5+4)
        else
            write(762,'(5e14.5e3)') rg(i*5+0),rg(i*5+1),rg(i*5+2),rg(i*5+3),rg(i*5+4)
        endif
    enddo
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(762,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0)
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(762,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1)
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(762,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1),rg(cnt_nr_forUFILE*5+2)
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(762,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1),rg(cnt_nr_forUFILE*5+2),rg(cnt_nr_forUFILE*5+3)

    write(762,'(5e14.5e3)') 1.01d0
      
    do i=0,cnt_nr_forUFILE-1
        if(i==0) then
            write(762,'(5e14.5e3)') (rt(1,2)+(rt(1,2)-rt(2,2)))*1d3,rt(i*5+1,2)*1d3,rt(i*5+2,2)*1d3,rt(i*5+3,2)*1d3,rt(i*5+4,2)*1d3
        else
            write(762,'(5e14.5e3)') rt(i*5+0,2)*1d3,rt(i*5+1,2)*1d3,rt(i*5+2,2)*1d3,rt(i*5+3,2)*1d3,rt(i*5+4,2)*1d3
        endif          
    enddo

    if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(762,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,2)*1d3
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(762,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,2)*1d3,rt(cnt_nr_forUFILE*5+1,2)*1d3
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(762,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,2)*1d3,rt(cnt_nr_forUFILE*5+1,2)*1d3,rt(cnt_nr_forUFILE*5+2,2)*1d3
    if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(762,'(1e14.5e3)') rt(cnt_nr_forUFILE*5+0,2)*1d3,rt(cnt_nr_forUFILE*5+1,2)*1d3,rt(cnt_nr_forUFILE*5+2,2)*1d3,rt(cnt_nr_forUFILE*5+3,2)*1d3

    close(761)
    close(762)

    return
    
end SUBROUTINE teti2ufiles

SUBROUTINE CAMPARISON_TETI_T0
!===========================================================================      
!=======   calculate the difference between initial T and lastest T.  ======      
!=======     for determing the Constant Factor of TB models.         ======      
!===========================================================================

! NOW We assume only one ion.(H only plasma)
    USE TRCOMM, ONLY : nrmax,nsmax,rg,rt,MDLKAI,CK0,CK1
!    USE T3D_FILE_IO
!    integer(4),save :: cnt_teti_comparison=0
    integer(4) :: i,j,nr,ns
    integer(4) :: cnt_len_filename=0
   
    real(8),allocatable,save :: TT0_save(:,:)
    real(8), dimension(nrmax,NSMAX) :: Err_TT0andLastTT
    real(8), dimension(nsmax) ::  SUM_Err_TT0andLastTT
    real(8) :: TOTAL_SUM_Err_TT0andLastTT    

    if(CAMPARISON_TETI_T0_file_flg == 0)then
        allocate(TT0_save(nrmax,nsmax))
        do i=1,nrmax
            do j=1,NSMAX
                TT0_save(i,j)=RT(i,j)
            enddo
        enddo
        
        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_LHDSHOT)+1
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_SHOTTIME)+1
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_TBMODEL)+1
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_CKeCKi(1))*2)+1
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1))*4)+1
        if(flg_OUTPUT_FILE_NAME_muemui==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_mu(1))*2)+1
    
        cnt_len_filename=cnt_len_filename+len("ERF_Difference_TT0andLastestT.txt")
    
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_ERF_TT0)

        STR_OUTPUT_FILE_ERF_TT0=""

        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//STR_OUTPUT_FILE_NAME_LHDSHOT//"_"
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//STR_OUTPUT_FILE_NAME_SHOTTIME//"_"
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//STR_OUTPUT_FILE_NAME_TBMODEL//"_"
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//STR_OUTPUT_FILE_NAME_CKeCKi(1)//STR_OUTPUT_FILE_NAME_CKeCKi(2)//"_"
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4)//"_"
        if(flg_OUTPUT_FILE_NAME_muemui==1)      STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//STR_OUTPUT_FILE_NAME_mu(1)//STR_OUTPUT_FILE_NAME_mu(2)//"_"

        STR_OUTPUT_FILE_ERF_TT0=STR_OUTPUT_FILE_ERF_TT0//"ERF_Difference_TT0andLastestT.txt"        
        CAMPARISON_TETI_T0_file_flg=1
        return
    endif    
    
    SUM_Err_TT0andLastTT=0.d0
    TOTAL_SUM_Err_TT0andLastTT=0.d0
    do nr=1,nrmax
        do ns=1,2
            Err_TT0andLastTT(nr,ns)=((TT0_save(nr,ns)-RT(nr,ns))/TT0_save(nr,ns))**2
            SUM_Err_TT0andLastTT(ns) = SUM_Err_TT0andLastTT(ns) + Err_TT0andLastTT(nr,ns)
        enddo
    enddo
    TOTAL_SUM_Err_TT0andLastTT=SUM_Err_TT0andLastTT(1)+SUM_Err_TT0andLastTT(2)

    open(780, file=STR_OUTPUT_FILE_ERF_TT0, status='replace')
    if(MDLKAI==209 .or. MDLKAI==218 .or. MDLKAI==219)then
    write(780,'(a1,a13,a3,i5)') '#','TB_MODEL', ':', MDLKAI
    write(780,'(a1,a13,a3,e20.10)') '#', 'mue = ',':',gB_mu_e
    write(780,'(a1,a13,a3,e20.10)') '#', 'mui = ',':',gB_mu_i
    write(780,'(a1,a13,a3,e20.10)') '#', 'CKe_GB_0.0 = ',':',CKe0
    write(780,'(a1,a13,a3,e20.10)') '#', 'CKe_GB_1.5 = ',':',CKe1_5
    write(780,'(a1,a13,a3,e20.10)') '#', 'CKi_GB_0.0 = ',':',CKi0
    write(780,'(a1,a13,a3,e20.10)') '#', 'CKi_GB_1.5 = ',':',CKi1_5
    else
    write(780,'(a1,a10,a3,i5)') '#','TB_MODEL', ':', MDLKAI
    write(780,'(a1,a10,a3,e20.10)') '#', 'CK0 = ',':',CK0
    write(780,'(a1,a10,a3,e20.10)') '#', 'CK1 = ',':',CK1
    write(780,'(a1,a10,a3,e20.10)') '#', 'mue = ',':',gB_mu_e
    write(780,'(a1,a10,a3,e20.10)') '#', 'mui = ',':',gB_mu_i
    endif
    
    write(780,'(a1,a5,11a20)') '#','nr','rg(nr)','TT0_ele','TT0_ion','TTsim_ele','TTsim_ion','ERF_ele','ERF_ion','SUM_ERF_ele','SUM_ERF_ion','SUM_ERF_TOT','av_sqrt_SEF'

    do nr=1,nrmax
       write(780,'(i6,11e20.10)') nr, rg(nr),TT0_save(nr,1),TT0_save(nr,2),rt(nr,1),rt(nr,2),Err_TT0andLastTT(nr,1),Err_TT0andLastTT(nr,2),SUM_Err_TT0andLastTT(1),SUM_Err_TT0andLastTT(2),TOTAL_SUM_Err_TT0andLastTT,sqrt(TOTAL_SUM_Err_TT0andLastTT/(1.d0*nrmax))
    enddo

    close(780)
    deallocate(TT0_save)
    deallocate(STR_OUTPUT_FILE_ERF_TT0)    
    
END SUBROUTINE CAMPARISON_TETI_T0

subroutine t3d_write_file_evo
    use trcomm, only : rn, rt, nrmax, ak, av, ad, rm, er, t,WPT,TAUE1,Q0,POHT,PNBT,PRFST,PNFT,PEXST
    implicit none
    real(8) :: PINT
    integer(4) :: cnt_len_filename=0    

    PINT=POHT+PNBT+PRFST+PNFT+PEXST

    if(FILE_WP_EVO_file_flg==0)then
        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_LHDSHOT)+1
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_SHOTTIME)+1
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_TBMODEL)+1
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_CKeCKi(1))*2)+1
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1))*4)+1
        if(flg_OUTPUT_FILE_NAME_muemui==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_mu(1))*2)+1
    
        cnt_len_filename=cnt_len_filename+len("tr_out_WP_evo.txt")
    
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_WP_EVO)

        STR_OUTPUT_FILE_WP_EVO=""

        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//STR_OUTPUT_FILE_NAME_LHDSHOT//"_"
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//STR_OUTPUT_FILE_NAME_SHOTTIME//"_"
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//STR_OUTPUT_FILE_NAME_TBMODEL//"_"
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//STR_OUTPUT_FILE_NAME_CKeCKi(1)//STR_OUTPUT_FILE_NAME_CKeCKi(2)//"_"
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4)//"_"
        if(flg_OUTPUT_FILE_NAME_muemui==1)      STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//STR_OUTPUT_FILE_NAME_mu(1)//STR_OUTPUT_FILE_NAME_mu(2)//"_"

        STR_OUTPUT_FILE_WP_EVO=STR_OUTPUT_FILE_WP_EVO//"tr_out_WP_evo.txt"
        
        open(154,file=STR_OUTPUT_FILE_WP_EVO,status='replace')
        write(154,'(a,a21,21a22)') '#','1:t','2:WPT','3:TAUE1','4:Q0','5:PINT','6:POHT','7:PNBT','8:PRFST','9:PNFT','10:PEXST','11:Te_r000','12:Ti_r000','13:Te_r010','14:Ti_r010','15:Te_r025','16:Ti_r025','17:Te_r050','18:Ti_r050','19:Te_r075','20:Ti_r075','21:Te_r090','22:Ti_r090'
        FILE_WP_EVO_file_flg=1
        close(154)
    endif

    open(154,file=STR_OUTPUT_FILE_WP_EVO,position='append')
    write(154,'(22e22.14)') t,WPT,TAUE1,Q0,PINT,POHT,PNBT,PRFST,PNFT,PEXST,rt(1,1),rt(1,2),rt(int((nrmax*0.1)),1),rt(int((nrmax*0.1)),2),rt(int((nrmax*0.25)),1),rt(int((nrmax*0.25)),2),rt(int((nrmax*0.5)),1),rt(int((nrmax*0.5)),2),rt(int((nrmax*0.75)),1),rt(int((nrmax*0.75)),2),rt(int((nrmax*0.9)),1),rt(int((nrmax*0.9)),2)
    close(154)

end subroutine t3d_write_file_evo

subroutine OUTPUT_FILE_TRCOEFLOG_LASTEST
    integer(4) :: cnt_len_filename=0
    integer(4) :: nr
    
    if(TRCOEFLOG_LASTEST_file_flg==0)then
        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_LHDSHOT)+1
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_SHOTTIME)+1
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_TBMODEL)+1
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_CKeCKi(1))*2)+1
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1))*4)+1
        if(flg_OUTPUT_FILE_NAME_muemui==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_mu(1))*2)+1
    
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_lastest.txt")
    
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_LASTEST)

        STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=""

        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//STR_OUTPUT_FILE_NAME_LHDSHOT//"_"
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//STR_OUTPUT_FILE_NAME_SHOTTIME//"_"
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//STR_OUTPUT_FILE_NAME_TBMODEL//"_"
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//STR_OUTPUT_FILE_NAME_CKeCKi(1)//STR_OUTPUT_FILE_NAME_CKeCKi(2)//"_"
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4)//"_"
        if(flg_OUTPUT_FILE_NAME_muemui==1)      STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//STR_OUTPUT_FILE_NAME_mu(1)//STR_OUTPUT_FILE_NAME_mu(2)//"_"

        STR_OUTPUT_FILE_TRCOEFLOG_LASTEST=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST//"TRCOEFLOG_lastest.txt"
        TRCOEFLOG_LASTEST_file_flg=1
    endif

    if(TRCOEF_01==1)then
            open(717,file=STR_OUTPUT_FILE_TRCOEFLOG_LASTEST,status='replace')
            if(MDLKAI==209 .or. MDLKAI==218 .or. MDLKAI==219)then
            write(717,'(a6,53a20)')'#1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)',&
                &   '11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)', &
                &   '21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)',&
                &   '31:rt(nr,2)','32:er(nr)','33:PEX(nr,1)','34:PEX(nr,2)','35:PIN(nr,1)','36:PIN(nr,2)',&
                &   '37:pflx(nr,1)','38:pflx(nr,2)','39:hflx(nr,1)','40:hflx(nr,2)','41:dd(nr,1)','42:dd(nr,2)','43:vv(nr,1)','44:vv(nr,2)','45:chi0(nr,1)','46:chi0(nr,2)','47:chi(nr,1)','48:chi(nr,2)','49:kv(nr,1)','50:kv(nr,2)',&
                &   '51:AKDWeGB0','52:AKDWeGB1.5','53:AKDWiGB0','54:AKDWeGB1.5'
            do nr=1,nrmax
                write(717,'(i6,53e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)/1.d3, PEX(nr,1)/1.d6,PEX(nr,2)/1.d6,PIN(nr,1)/1.d6,PIN(nr,2)/1.d6,&
                &   pflx_save(nr,0),pflx_save(nr,1), hflx_save(nr,0),hflx_save(nr,1), dd_save(nr,0), dd_save(nr,1), vv_save(nr,0), vv_save(nr,1), chi0_save(nr,0), chi0_save(nr,1),chi_save(nr,0),chi_save(nr,1), kv_save(nr,0), kv_save(nr,1),&
                &   AKDWe_GB_0(nr),AKDWe_GB_1_5(nr),AKDWi_GB_0(nr),AKDWi_GB_1_5(nr)
            enddo
            else
            write(717,'(a6,49a20)')'#1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)',&
                &   '11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)', &
                &   '21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)',&
                &   '31:rt(nr,2)','32:er(nr)','33:PEX(nr,1)','34:PEX(nr,2)','35:PIN(nr,1)','36:PIN(nr,2)',&
                &   '37:pflx(nr,1)','38:pflx(nr,2)','39:hflx(nr,1)','40:hflx(nr,2)','41:dd(nr,1)','42:dd(nr,2)','43:vv(nr,1)','44:vv(nr,2)','45:chi0(nr,1)','46:chi0(nr,2)','47:chi(nr,1)','48:chi(nr,2)','49:kv(nr,1)','50:kv(nr,2)'
            do nr=1,nrmax
                write(717,'(i6,49e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)/1.d3, PEX(nr,1)/1.d6,PEX(nr,2)/1.d6,PIN(nr,1)/1.d6,PIN(nr,2)/1.d6,&
                &   pflx_save(nr,0),pflx_save(nr,1), hflx_save(nr,0),hflx_save(nr,1), dd_save(nr,0), dd_save(nr,1), vv_save(nr,0), vv_save(nr,1), chi0_save(nr,0), chi0_save(nr,1),chi_save(nr,0),chi_save(nr,1), kv_save(nr,0), kv_save(nr,1)
            enddo
            endif

            close(717)
      endif
end subroutine   

subroutine OUTPUT_FILE_TRCOEFLOG_EVO

    integer(4) :: cnt_len_filename=0,cnt_use_filename
    integer(4) :: nr
    
    if(TRCOEFLOG_EVO_file_flg==0 )then    
        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_LHDSHOT)+1
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_SHOTTIME)+1
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_TBMODEL)+1
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_CKeCKi(1))*2)+1
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1))*4)+1
        if(flg_OUTPUT_FILE_NAME_muemui==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_mu(1))*2)+1
    
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_EVO.txt")
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_EVO)
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_EVO_rho010.txt")
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010)
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_EVO_rho025.txt")
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025)
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_EVO_rho050.txt")
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050)
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_EVO_rho075.txt")
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075)
        cnt_len_filename=cnt_len_filename+len("TRCOEFLOG_EVO_rho090.txt")
        allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090)

        STR_OUTPUT_FILE_TRCOEFLOG_EVO=""
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010=""
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025=""
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050=""
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075=""
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090=""

        if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//STR_OUTPUT_FILE_NAME_LHDSHOT//"_"
        if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//STR_OUTPUT_FILE_NAME_SHOTTIME//"_"
        if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//STR_OUTPUT_FILE_NAME_TBMODEL//"_"
        if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//STR_OUTPUT_FILE_NAME_CKeCKi(1)//STR_OUTPUT_FILE_NAME_CKeCKi(2)//"_"
        if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4)//"_"
        if(flg_OUTPUT_FILE_NAME_muemui==1)      STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//STR_OUTPUT_FILE_NAME_mu(1)//STR_OUTPUT_FILE_NAME_mu(2)//"_"

        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010=STR_OUTPUT_FILE_TRCOEFLOG_EVO//"TRCOEFLOG_EVO_rho010.txt"
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025=STR_OUTPUT_FILE_TRCOEFLOG_EVO//"TRCOEFLOG_EVO_rho025.txt"
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050=STR_OUTPUT_FILE_TRCOEFLOG_EVO//"TRCOEFLOG_EVO_rho050.txt"
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075=STR_OUTPUT_FILE_TRCOEFLOG_EVO//"TRCOEFLOG_EVO_rho075.txt"
        STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090=STR_OUTPUT_FILE_TRCOEFLOG_EVO//"TRCOEFLOG_EVO_rho090.txt"
        STR_OUTPUT_FILE_TRCOEFLOG_EVO=STR_OUTPUT_FILE_TRCOEFLOG_EVO//"TRCOEFLOG_EVO.txt"

        open(711,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010,status='replace')
        open(712,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025,status='replace')
        open(713,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050,status='replace')
        open(714,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075,status='replace')
        open(715,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090,status='replace')
        open(716,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO,status='replace')
        write(711,'(a6,31a20)')'# 1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)','11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)','21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)','31:rt(nr,2)','32:er(nr)'
        write(712,'(a6,31a20)')'# 1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)','11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)','21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)','31:rt(nr,2)','32:er(nr)'
        write(713,'(a6,31a20)')'# 1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)','11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)','21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)','31:rt(nr,2)','32:er(nr)'
        write(714,'(a6,31a20)')'# 1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)','11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)','21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)','31:rt(nr,2)','32:er(nr)'
        write(715,'(a6,31a20)')'# 1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)','11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)','21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)','31:rt(nr,2)','32:er(nr)'
        close(711)
        close(712)
        close(713)
        close(714)
        close(715)
        close(716)

        TRCOEFLOG_EVO_file_flg=1
    endif
      
    if(TRCOEF_01==1)then
        open(716,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO,position='APPEND')
        write(716,'(a6,49a20)')'# 1:nr','2:time','3:rg(nr)','4:AD(nr,1)','5:ADNC(nr,1)','6:ADDW(nr,1)','7:AK(nr,1)','8:AKNC(nr,1)','9:AKDW(nr,1)','10:AV(nr,1)','11:AVNC(nr,1)','12:AVDW(nr,1)','13:AVK(nr,1)','14:AVKNC(nr,1)','15:AVKDW(nr,1)','16:AD(nr,2)','17:ADNC(nr,2)','18:ADDW(nr,2)','19:AK(nr,2)','20:AKNC(nr,2)','21:AKDW(nr,2)','22:AV(nr,2)','23:AVNC(nr,2)','24:AVDW(nr,2)','25:AVK(nr,2)','26:AVKNC(nr,2)','27:AVKDW(nr,2)','28:rn(nr,1)','29:rn(nr,2)','30:rt(nr,1)','31:rt(nr,2)','32:er(nr)','33:PEX(nr,1)','34:PEX(nr,2)','35:PIN(nr,1)','36:PIN(nr,2)',&
                            &   '37:pflx(nr,1)','38:pflx(nr,2)','39:hflx(nr,1)','40:hflx(nr,2)','41:dd(nr,1)','42:dd(nr,2)','43:vv(nr,1)','44:vv(nr,2)','45:chi0(nr,1)','46:chi0(nr,2)','47:chi(nr,1)','48:chi(nr,2)','49:kv(nr,1)','50:kv(nr,2)'
        do nr=1,nrmax
            write(716,'(i6,49e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr),pex(nr,1),pex(nr,2),pin(nr,1),pin(nr,2),&
                            &   pflx_save(nr,0),pflx_save(nr,1), hflx_save(nr,0),hflx_save(nr,1), dd_save(nr,0), dd_save(nr,1), vv_save(nr,0), vv_save(nr,1), chi0_save(nr,0), chi0_save(nr,1),chi_save(nr,0),chi_save(nr,1), kv_save(nr,0), kv_save(nr,1)
        enddo
        write(716,*)
        close(716)

        open(711,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho010,position='APPEND')
        open(712,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho025,position='APPEND')
        open(713,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho050,position='APPEND')
        open(714,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho075,position='APPEND')
        open(715,file=STR_OUTPUT_FILE_TRCOEFLOG_EVO_rho090,position='APPEND')
        nr=int((nrmax*0.1))
        write(711,'(i6,31e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)
        nr=int((nrmax*0.25))
        write(712,'(i6,31e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)
        nr=int((nrmax*0.50))
        write(713,'(i6,31e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)
        nr=int((nrmax*0.75))
        write(714,'(i6,31e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)
        nr=int((nrmax*0.90))
        write(715,'(i6,31e20.10)')nr,T,rg(nr),AD(nr,1),ADNC(nr,1),ADDW(nr,1),AK(nr,1),AKNC(nr,1),AKDW(nr,1),AV(nr,1),AVNC(nr,1),AVDW(nr,1),AVK(nr,1),AVKNC(nr,1),AVKDW(nr,1),AD(nr,2),ADNC(nr,2),ADDW(nr,2),AK(nr,2),AKNC(nr,2),AKDW(nr,2),AV(nr,2),AVNC(nr,2),AVDW(nr,2),AVK(nr,2),AVKNC(nr,2),AVKDW(nr,2),rn(nr,1),rn(nr,2),rt(nr,1),rt(nr,2),er(nr)
        close(711)
        close(712)
        close(713)
        close(714)
        close(715)
    endif      
      
end subroutine     

end module
