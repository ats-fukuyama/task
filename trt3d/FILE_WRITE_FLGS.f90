!
!	FILE_WRITE_FLGS.f90
!	task3D_modified
!
!	Created by WAKASA Arimitsu on 11/08/28.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!

module FILE_WRITE_FLGS
use TRCOMM,only:NRMAX
    implicit none
    integer(4),parameter:: NRMAX4M_FILE_WRITE_FLGS=500
    integer(4) :: TRCOEF_01
    
    
    integer(4) ,save ::PEX_flg=0
    real(8),dimension(1:NRMAX4M_FILE_WRITE_FLGS, 0:10),save :: pflx_save, hflx_save, dd_save, vv_save, chi0_save,chi_save, kv_save ! for check the NC parameters
    real(8),save :: gb_mu,gb_mu_e,gb_mu_i
    character(38) ,save:: FILENAME_CK0CK1mu
!    character(23) ,save:: FILENAME_CK0CK1
!    character(59) ,save:: FILENAME_CK0CK1mu
    character(51) ,save:: FILENAME_CK0CK1
    character(15) ,save:: FILENAME_gbmu
    character(6),save:: FILENAME_shotnum
    integer(4),save :: flg_t3d_er_calc_timing
    integer(4),save :: flg_ERF_Difference_Coptimize
    real(8),save  ::pre_AKDWIL=0.d0,pre_AKDWEL=0.d0
    real(8),save :: CKe0, CKe1_5, CKi0, CKi1_5
    real(8),dimension(1:NRMAX4M_FILE_WRITE_FLGS),save :: AKDWe_GB_0=0.d0,AKDWe_GB_1_5=0.d0,AKDWi_GB_0=0.d0,AKDWi_GB_1_5=0.d0
    real(8) :: tmpAKDWe_GB_0=0.d0,tmpAKDWe_GB_1_5=0.d0,tmpAKDWi_GB_0=0.d0,tmpAKDWi_GB_1_5=0.d0
    
    character(81) :: FILENAME_LOG_LASTEST
    character(22) :: FILENAME_LOG_LASTEST_endpart  
    
end module