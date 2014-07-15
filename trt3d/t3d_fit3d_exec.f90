!     ***********************************************************

!           FIT 3D calculation 

!     ***********************************************************

    SUBROUTINE FIT3D_EXEC

    USE TRCOMM, ONLY : nrmax,rt,rn,nt,rg,PEX,KUFDCG,T,BB,RR
    use task3d_tentative_flags

    integer(4) :: nr, i
    integer(4) :: cnt_nr_forUFILE    
!       FIT3DとTASK3Dを連携する拡張に使用した変数群
    real(8), save, dimension(1000) :: pex_nbi_e, pex_nbi_i,pex_save_e,pex_save_i
    real(8),dimension(nrmax) :: rteforc,rtiforc,rneforc,rgforc
    integer(4),save::flg_FIT3D_calc_timing_initialize=0 !0:未初期化, 1〜:初期化済み＆何度目のFIT3D計算タイミングかの記録

!================================================================================
!==========  FIT3D calculation for PEX                                   =========
!==========     NBI加熱吸収パワー分布をFIT3Dで計算するかの判定，及び，             =========
!==========     FIT3Dでの計算タイミングの判定，FIT3Dへの受け渡し                 =========
!==========     各パラメータは現在(2011,12,06)                               =========
!==========     >> module task3d_tentative_flagsに記述                    =========
!======================================================================@2011,12,06   

    if(FIT3D_tentative_flags==1)then
!       FIT3Dで加熱吸収パワーの計算を行うか否かの判定．
!       FIT3D_tentative_falgs=0: unuse FIT3D and use UFILE for PEX
!       FIT3D_tentative_falgs=1: use FIT3D and unuse UFILE for PEX
        if(flg_FIT3D_calc_timing_initialize==0)then
            do nr=1,nrmax
                pex_save_e(nr)=pex(nr,1)
                pex_save_i(nr)=pex(nr,2)
            enddo        
            call FIT3D_CALC_TIMING_initialize
            flg_FIT3D_calc_timing_initialize=1
 !===============================================================================
 !==========    初回計算時にFIT3Dの計算結果をUfile形式で書きだす.               =========
 !==========    Ufileを用いた定常解析に使用.                                 =========
 !==========    竹田くんの卒論用に制作＠20121220．                           =========
 !=====================================ここから======竹田くんの卒論用に制作＠20121220
            if(fit3dtoufile_flg==1) then
                do nr=1,nrmax
                    rgforc(nr)=rg(nr)
                    rteforc(nr)=RT(nr,1)
                    rtiforc(nr)=RT(nr,2)
                    rneforc(nr)=RN(nr,1)*1d1
                enddo
                call tr_fit3d_connect(KUFDCG,T,nrmax,BB,RR,rgforc,rteforc,rtiforc,rneforc,pex_nbi_e,pex_nbi_i)
                do nr=1,nrmax
                    if(pex_nbi_e(nr) <0.d0) pex_nbi_e(nr)=0.d0
                    if(pex_nbi_i(nr) <0.d0) pex_nbi_i(nr)=0.d0
                enddo
            
                do nr=1,nrmax
                    pex(nr,1)=pex_nbi_e(nr)
                    pex(nr,2)=pex_nbi_i(nr)
                enddo        
            
                open(761,file='Handover_PEXE_UFILE.txt',status='replace')
                open(762,file='Handover_PEXI_UFILE.txt',status='replace')
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

                write(761,'(5e14.5e3)') 4.24
      
                do i=0,cnt_nr_forUFILE-1
                    if(i==0) then
                        if((pex(1,1)+(pex(1,1)-pex(2,1)))<0)then
                            write(761,'(5e14.5e3)') 0d0,pex(i*5+1,1),pex(i*5+2,1),pex(i*5+3,1),pex(i*5+4,1)
                        else
                            write(761,'(5e14.5e3)') (pex(1,1)+(pex(1,1)-pex(2,1))),pex(i*5+1,1),pex(i*5+2,1),pex(i*5+3,1),pex(i*5+4,1)
                        endif
                    else
                        write(761,'(5e14.5e3)') pex(i*5+0,1),pex(i*5+1,1),pex(i*5+2,1),pex(i*5+3,1),pex(i*5+4,1)
                    endif      
                enddo
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(761,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,1)
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(761,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,1),pex(cnt_nr_forUFILE*5+1,1)
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(761,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,1),pex(cnt_nr_forUFILE*5+1,1),pex(cnt_nr_forUFILE*5+2,1)
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(761,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,1),pex(cnt_nr_forUFILE*5+1,1),pex(cnt_nr_forUFILE*5+2,1),pex(cnt_nr_forUFILE*5+3,1)
                close(761)
            
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

                write(762,'(5e14.5e3)') 4.24
      
                do i=0,cnt_nr_forUFILE-1
                    if(i==0) then
                        if((pex(1,2)+(pex(1,2)-pex(2,2)))<0)then
                            write(762,'(5e14.5e3)') 0d0,pex(i*5+1,2),pex(i*5+2,2),pex(i*5+3,2),pex(i*5+4,2)
                        else
                            write(762,'(5e14.5e3)') (pex(1,2)+(pex(1,2)-pex(2,2))),pex(i*5+1,2),pex(i*5+2,2),pex(i*5+3,2),pex(i*5+4,2)
                        endif
                    else
                        write(762,'(5e14.5e3)') pex(i*5+0,2),pex(i*5+1,2),pex(i*5+2,2),pex(i*5+3,2),pex(i*5+4,2)
                    endif          
                enddo
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(762,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,2)
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(762,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,2),pex(cnt_nr_forUFILE*5+1,2)
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(762,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,2),pex(cnt_nr_forUFILE*5+1,2),pex(cnt_nr_forUFILE*5+2,2)
                if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(762,'(1e14.5e3)') pex(cnt_nr_forUFILE*5+0,2),pex(cnt_nr_forUFILE*5+1,2),pex(cnt_nr_forUFILE*5+2,2),pex(cnt_nr_forUFILE*5+3,2)
                close(762)

                if(ntPolynomial2UFILE_flg==1) then
                    do nr=1,nrmax
                        rt(nr,1)=rteforc(nr)
                        rt(nr,2)=rtiforc(nr)
                        rn(nr,1)=rneforc(nr)
                    enddo
                
                    open(761,file='Handover_TE_Poly2UFILE.txt',status='replace')
                    open(762,file='Handover_TI_Poly2UFILE.txt',status='replace')
                    open(763,file='Handover_nn_Poly2UFILE.txt',status='replace')

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

                    write(761,'(5e14.5e3)') 4.24d0
      
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
                    close(761)
                  
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

                    write(762,'(5e14.5e3)') 4.24d0
      
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
                    close(762)
                
                    do i=0,cnt_nr_forUFILE-1
                        if(i==0) then
                            write(763,'(5e14.5e3)') 0.d0,rg(i*5+1),rg(i*5+2),rg(i*5+3),rg(i*5+4)
                        else
                            write(763,'(5e14.5e3)') rg(i*5+0),rg(i*5+1),rg(i*5+2),rg(i*5+3),rg(i*5+4)
                        endif
                    enddo
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(763,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0)
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(763,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1)
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(763,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1),rg(cnt_nr_forUFILE*5+2)
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(763,'(1e14.5e3)') rg(cnt_nr_forUFILE*5+0),rg(cnt_nr_forUFILE*5+1),rg(cnt_nr_forUFILE*5+2),rg(cnt_nr_forUFILE*5+3)

                    write(763,'(5e14.5e3)') 4.24d0
      
                    do i=0,cnt_nr_forUFILE-1
                        if(i==0) then
                            write(763,'(5e14.5e3)') (rn(1,1)+(rn(1,1)-rn(2,1)))*1d20,rn(i*5+1,1)*1d20,rn(i*5+2,1)*1d20,rn(i*5+3,1)*1d20,rn(i*5+4,1)*1d20
                        else
                            write(763,'(5e14.5e3)') rn(i*5+0,1)*1d20,rn(i*5+1,1)*1d20,rn(i*5+2,1)*1d20,rn(i*5+3,1)*1d20,rn(i*5+4,1)*1d20
                        endif          
                    enddo
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==1)write(763,'(1e14.5e3)') rn(cnt_nr_forUFILE*5+0,1)*1d20
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==2)write(763,'(1e14.5e3)') rn(cnt_nr_forUFILE*5+0,1)*1d20,rn(cnt_nr_forUFILE*5+1,1)*1d20
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==3)write(763,'(1e14.5e3)') rn(cnt_nr_forUFILE*5+0,1)*1d20,rn(cnt_nr_forUFILE*5+1,1)*1d20,rn(cnt_nr_forUFILE*5+2,1)*1d20
                    if (((nrmax+1)-(cnt_nr_forUFILE*5))==4)write(763,'(1e14.5e3)') rn(cnt_nr_forUFILE*5+0,1)*1d20,rn(cnt_nr_forUFILE*5+1,1)*1d20,rn(cnt_nr_forUFILE*5+2,1)*1d20,rn(cnt_nr_forUFILE*5+3,1)*1d20            
                    close(763)
                endif

                print*,'Completed FIT3D to UFILE. and STOPPED.'           
                stop
      
            endif
        endif
!===============================================================================
!=====================================ここまで======竹田くんの卒論用に制作＠20121220

        if(nt == FIT3D_CALC_TIMING(flg_FIT3D_calc_timing_initialize))then 
            do nr=1,nrmax
                rgforc(nr)=rg(nr)
                rteforc(nr)=RT(nr,1)
                rtiforc(nr)=RT(nr,2)
                rneforc(nr)=RN(nr,1)*1d1
            enddo
            call tr_fit3d_connect(KUFDCG,T,nrmax,BB,RR,rgforc,rteforc,rtiforc,rneforc,pex_nbi_e,pex_nbi_i)
            do nr=1,nrmax
                if(pex_nbi_e(nr) <0.d0) pex_nbi_e(nr)=0.d0
                if(pex_nbi_i(nr) <0.d0) pex_nbi_i(nr)=0.d0
            enddo
            do nr=1,nrmax
                pex_save_e(nr)=pex_nbi_e(nr)
                pex_save_i(nr)=pex_nbi_i(nr)
                pex(nr,1)=pex_save_e(nr)
                pex(nr,2)=pex_save_i(nr)        
            enddo
            flg_FIT3D_calc_timing_initialize=flg_FIT3D_calc_timing_initialize+1
        else
            do nr=1,nrmax
                pex(nr,1)=pex_save_e(nr)
                pex(nr,2)=pex_save_i(nr)
            enddo
        endif
    endif
    
!===============================================================================
!===============================================================================
!======================================================================@2011,12,06   

      RETURN
      END SUBROUTINE FIT3D_EXEC
    