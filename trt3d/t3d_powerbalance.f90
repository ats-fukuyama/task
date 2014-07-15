!     ***********************************************************

!           Power Balance calculation 

!     ***********************************************************

SUBROUTINE POWERBALANCE_CALC
    
    USE TRCOMM, ONLY : nrmax,nsmax,er,RN,RT,AEE,PI,EPS0,RTM,AMZ,DVRHO,DR,PIN,RA,BB,RG,RR
    use task3d_tentative_flags
    use ncdn_mod, only : ncdn_analysis
    use T3D_FILE_IO
    
    real(8),dimension(nrmax,nsmax)::VolumeIntegral_PIN_wi_EQ,chi_tot_enrgyblance,drTemp,QQNC,drho_dtem
    real(8), dimension(1:300,0:10):: ddd,vvv,chi0,chi,kv,pflx,hflx
    real(8),dimension(10,10,500)::TAU_EQ
    real(8)::COEF,COULOG
    real(8)::BETA_AVRGV,TOTAL_V
    integer(4) :: nr, ns, ns1
    integer(4) :: cnt_len_filename=0     

    ddd=0.d0
    vvv=0.d0
    chi=0.d0
    kv=0.d0
    pflx=0.d0
    hflx=0.d0

	print*
    print*,'============================================'        
    print*,'====> Power Balance Calculating route.'        
    print*,'============================================'        
    print*
    
    chi=0.d0
    chi0=0.d0
    
    COEF = AEE**4*1.D20/(3.D0*SQRT(2.D0*PI)*PI*EPS0**2)
    do nr=1,nrmax
        DO NS=1,NSMAX
            DO NS1=1,NSMAX
                IF(NS1.NE.NS) THEN
                    TAU_EQ(ns,ns1,nr)=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1)) &
            &                  *COULOG(NS,NS1,RN(NR,1),RT(NR,NS)) *RN(NR,NS1)
                ENDIF
            ENDDO
        ENDDO
    enddo
    ! TAU_EQ(2,2,nr)を電子イオン間緩和時間として再計算
    do nr=1,nrmax
        TAU_EQ(2,2,nr)=(3.d0*sqrt(2.d0)*PI**(1.5d0)*(8.85418710d-12)**2.d0*(9.10938188d-31)*(1.672621d-27)) &
            &   /(rn(nr,1)*1.d20*(1.6021d-19**4.d0)*COULOG(1,2,RN(NR,1),RT(NR,1))) &
            &   *(((rt(nr,2)*1e3*1.6021*1.d-19)/1.672621d-27)+((rt(nr,1)*1d3*1.6021*1.d-19)/9.1093d-31))**(1.5d0)
    enddo
    
    do nr=1,nrmax
        PIN_EQ(1,2,nr)=1.5d0*rn(nr,2)*1.d20*(-rt(nr,1)+rt(nr,2))*1.60217646*1d-16/TAU_EQ(2,2,nr)
        PIN_EQ(2,1,nr)=1.5d0*rn(nr,1)*1.d20*(-rt(nr,2)+rt(nr,1))*1.60217646*1d-16/TAU_EQ(2,2,nr)
    enddo

    call ncdn_analysis(pflx, hflx, ddd, vvv,chi0, chi, kv)
    
    VolumeIntegral_PIN_wi_EQ=0.d0
    VolumeIntegral_PIN_wi_EQ=0.d0

    nr=1
    VolumeIntegral_PIN_wi_EQ(nr,1)=DVRHO(nr)*DR*(PIN(nr,1)+PIN_EQ(1,2,nr))
    VolumeIntegral_PIN_wi_EQ(nr,2)=DVRHO(nr)*DR*(PIN(nr,2)+PIN_EQ(2,1,nr))
    do nr=1,nrmax
        QQNC(nr,1)=DVRHO(nr)/RA*hflx(nr,0)*1d20*1d3*1.60217646d-19
        QQNC(nr,2)=DVRHO(nr)/RA*hflx(nr,1)*1d20*1d3*1.60217646d-19
    enddo
    do nr=2,nrmax
        VolumeIntegral_PIN_wi_EQ(nr,1)=VolumeIntegral_PIN_wi_EQ(nr-1,1)+DVRHO(nr)*DR*(PIN(nr,1)+PIN_EQ(1,2,nr))
        VolumeIntegral_PIN_wi_EQ(nr,2)=VolumeIntegral_PIN_wi_EQ(nr-1,2)+DVRHO(nr)*DR*(PIN(nr,2)+PIN_EQ(2,1,nr))
    enddo 
    
    
    nr=1
    drTemp(nr,1)=((rt(nr+1,1)-rt(nr,1))/dr)/ra
    drTemp(nr,2)=((rt(nr+1,2)-rt(nr,2))/dr)/ra
    do nr=2,nrmax-1
        drTemp(nr,1)=((rt(nr+1,1)-rt(nr-1,1))/(2.d0*dr))/ra
        drTemp(nr,2)=((rt(nr+1,2)-rt(nr-1,2))/(2.d0*dr))/ra
    enddo
    nr=nrmax
    drTemp(nr,1)=((rt(nr,1)-rt(nr-1,1))/dr)/ra
    drTemp(nr,2)=((rt(nr,2)-rt(nr-1,2))/dr)/ra

    nr=1
    drho_dtem(nr,1)=((rt(nr+1,1)-rt(nr,1))/dr)
    drho_dtem(nr,2)=((rt(nr+1,2)-rt(nr,2))/dr)
    do nr=2,nrmax-1
        drho_dtem(nr,1)=((rt(nr+1,1)-rt(nr-1,1))/(2.d0*dr))
        drho_dtem(nr,2)=((rt(nr+1,2)-rt(nr-1,2))/(2.d0*dr))
    enddo
    nr=nrmax
    drho_dtem(nr,1)=((rt(nr,1)-rt(nr-1,1))/dr)
    drho_dtem(nr,2)=((rt(nr,2)-rt(nr-1,2))/dr)            
    
    do nr=1,nrmax
        chi_tot_enrgyblance(nr,1)=VolumeIntegral_PIN_wi_EQ(nr,1)/((-1.d0*drTemp(nr,1)*1.60217646*1d-16)*(DVRHO(nr)/ra)*rn(nr,1)*1.d20)
        chi_tot_enrgyblance(nr,2)=VolumeIntegral_PIN_wi_EQ(nr,2)/((-1.d0*drTemp(nr,2)*1.60217646*1d-16)*(DVRHO(nr)/ra)*rn(nr,2)*1.d20)
    enddo
    
    if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_LHDSHOT)+1
    if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_SHOTTIME)+1
    if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     cnt_len_filename=cnt_len_filename+len_trim(STR_OUTPUT_FILE_NAME_TBMODEL)+1
    if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_CKeCKi(1))*2)+1
    if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1))*4)+1
    if(flg_OUTPUT_FILE_NAME_muemui==1)      cnt_len_filename=cnt_len_filename+(len_trim(STR_OUTPUT_FILE_NAME_mu(1))*2)+1
    
    cnt_len_filename=cnt_len_filename+len("PowerBlanceCalc_results.txt")
    
    allocate(character(cnt_len_filename)::STR_OUTPUT_FILE_POWERBALANCE)

    STR_OUTPUT_FILE_POWERBALANCE=""

    if(flg_OUTPUT_FILE_NAME_LHDSHOT==1)     STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//STR_OUTPUT_FILE_NAME_LHDSHOT//"_"
    if(flg_OUTPUT_FILE_NAME_SHOTTIME==1)    STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//STR_OUTPUT_FILE_NAME_SHOTTIME//"_"
    if(flg_OUTPUT_FILE_NAME_TBMODEL==1)     STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//STR_OUTPUT_FILE_NAME_TBMODEL//"_"
    if(flg_OUTPUT_FILE_NAME_CKeCKi==1)      STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//STR_OUTPUT_FILE_NAME_CKeCKi(1)//STR_OUTPUT_FILE_NAME_CKeCKi(2)//"_"
    if(flg_OUTPUT_FILE_NAME_C1eC2eC1iC2i==1)STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(1)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(2)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(3)//STR_OUTPUT_FILE_NAME_C1eC2eC1iC2i(4)//"_"
    if(flg_OUTPUT_FILE_NAME_muemui==1)      STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//STR_OUTPUT_FILE_NAME_mu(1)//STR_OUTPUT_FILE_NAME_mu(2)//"_"

    STR_OUTPUT_FILE_POWERBALANCE=STR_OUTPUT_FILE_POWERBALANCE//"PowerBlanceCalc_results.txt"    
    print*,STR_OUTPUT_FILE_POWERBALANCE

    open(890,file=STR_OUTPUT_FILE_POWERBALANCE,status='replace')
    write(890,'(a6,58a28)')'1:nr','2:rg(nr)','3:rt(nr,1)','4:rt(nr,2)','5:dr_dtem(nr,1)','6:dr_dtem(nr,2)','7:VolIntPINwiPinEQ_e[MW]','8:VolIntPINwiPinEQ_i[MW]', &
     &        '9:chi_tot_Powblnc_e','10:chi_tot_Powblnc_i','11:D3_NC(nr,1)','12:D3_NC(nr,2)','13:D3-1.5D2(nr,1)','14:D3-1.5D2(nr,2)',&
     &        '15:chi_tot_eff_ele','16:chi_tot_eff_ion',&
     &        '17:chi_ano_eff_ele','18:chi_ano_eff_ion', &  
     &        '19:PIN(nr,1)','20:PIN(nr,2)','21:PIN_EQ(1,2,nr)','22:PIN_EQ(2,1,nr)','23:TAU_EQ(1,2,nr)','24:TAU_EQ(2,1,nr)','25:TAU_EQ(2,2,nr)',&
     &        '26:qNC*Surface(nr,1)[MW]','27:qNC*Surface(nr,2)[MW]','28:RVRHO(nt)','29:2*PI*R2*PI*rr','30:rn(nr,1)','31:hflx(nr,1)[MW/m2]','32:hflx(nr,2)[MW/m2]',&
     &        '33:Te-Ti[J]','34:PINwiEQ(nr,1)','35:PINwiEQ(nr,2)','36:1/TAU(2,2,nr)','37:n*nu*(Te-Ti)','38:((Te/me)+(Ti/mi))^(1.5)',&
     &        '39:Pe[Pa]','40:Pi[Pa]','41:BB0[T]','42:---','43:beta(nr)',&
     &        '44:beta*dV','45:(beta*dV)/V',&
     &        '46:er(nr)[keV]','47:Gamma_e[10^19/m2]','48:Gamma_i[10^19/m2]',&
     &        '49:drho_dtem(nr,1)','50:drho_dtem(nr,1)',&
     &        '51:DV(nr)','52:4PI^2rR[m^2]','53:Qe/DVRHO','54:Qi/DVRHO','55:Qe/4PI^2rR','56:Qi/4PI^2rR',&
     &        '57:pflx(nr,0)','58:pflx(nr,1)'
    BETA_AVRGV=0.d0
    TOTAL_V=0.d0
    do nr=1,nrmax
        TOTAL_V=TOTAL_V+(DVRHO(nr)*dr)
        BETA_AVRGV=(BETA_AVRGV+(2.d0*4.0d0*PI*1d-7)*(rt(nr,1)+rt(nr,2))*1d3*0.5*1.603*1d-19*rn(nr,1)*1d20/(BB*BB)*DVRHO(nr)*DR)
        write(890,'(i6,58e28.13)')nr,rg(nr),rt(nr,1),rt(nr,2),drtemp(nr,1),drtemp(nr,2),VolumeIntegral_PIN_wi_EQ(nr,1)/1d6,VolumeIntegral_PIN_wi_EQ(nr,2)/1d6,&
        &   chi_tot_enrgyblance(nr,1),chi_tot_enrgyblance(nr,2),chi0(nr,0),chi0(nr,1),chi(nr,0),chi(nr,1),&
        &   -(VolumeIntegral_PIN_wi_EQ(nr,1))/(rn(nr,1)*1d20*drtemp(nr,1)*1d3*1.603*1d-19)/(2.d0*PI*RR*2.d0*PI*rg(nr)*ra),-(VolumeIntegral_PIN_wi_EQ(nr,2))/(rn(nr,1)*1d20*drtemp(nr,2)*1d3*1.603*1d-19)/(2.d0*PI*RR*2.d0*PI*rg(nr)*ra),&
        &   -(VolumeIntegral_PIN_wi_EQ(nr,1))/(rn(nr,1)*1d20*drtemp(nr,1)*1d3*1.603*1d-19)/(2.d0*PI*RR*2.d0*PI*rg(nr)*ra)-chi(nr,0),-(VolumeIntegral_PIN_wi_EQ(nr,2))/(rn(nr,1)*1d20*drtemp(nr,2)*1d3*1.603*1d-19)/(2.d0*PI*RR*2.d0*PI*rg(nr)*ra)-chi(nr,1),&
        &   pin(nr,1),pin(nr,2),pin_eq(1,2,nr),pin_eq(2,1,nr),TAU_EQ(1,2,nr),TAU_EQ(2,1,nr),TAU_EQ(2,2,nr),QQNC(nr,1)/1d6,QQNC(nr,2)/1d6,&
        &   DVRHO(nr)/RA,2.d0*PI*RR*2.d0*PI*rg(nr)*RA,rn(nr,1),hflx(nr,0)*1d20*1d3*1.60217646d-19/1d6,hflx(nr,1)*1d20*1d3*1.60217646d-19/1d6,&
        &   ((rt(nr,1)-rt(nr,2))*1d3*1.6021d-19),(PIN(nr,1)+PIN_EQ(1,2,nr)),(PIN(nr,2)+PIN_EQ(2,1,nr)),(1.d0/TAU_EQ(2,2,nr)),rn(nr,1)*(1.d0/TAU_EQ(2,2,nr))*((rt(nr,1)-rt(nr,2))*1d3*1.6021d-19),1.d0/((((rt(nr,2)*1e3*1.6021*1.d-19)/1.672621d-27)+((rt(nr,1)*1d3*1.6021*1.d-19)/9.1093d-31))**(1.5d0)),&
        &   rn(nr,1)*1d20*rt(nr,1)*1d3*1.603*1d-19,rn(nr,1)*1d20*rt(nr,2)*1d3*1.603*1d-19,BB,0.d0,(2.d0*4.0d0*PI*1d-7)*(rt(nr,1)+rt(nr,2))*1d3*0.5*1.603*1d-19*rn(nr,1)*1d20/(BB*BB),&
        &   BETA_AVRGV,BETA_AVRGV/TOTAL_V,&
        &   er(nr)/1e3,pflx(nr,0)*1.d1,pflx(nr,1)*1.d1,&
        &   drho_dtem(nr,1),drho_dtem(nr,2),&
        &   dvrho(nr),4.d0*PI*PI*rg(nr)*ra*RR,VolumeIntegral_PIN_wi_EQ(nr,1)/1d6/dvrho(nr),VolumeIntegral_PIN_wi_EQ(nr,2)/1d6/dvrho(nr),VolumeIntegral_PIN_wi_EQ(nr,1)/1d6/(4.d0*PI*PI*rg(nr)*ra*RR),VolumeIntegral_PIN_wi_EQ(nr,2)/1d6/(4.d0*PI*PI*rg(nr)*ra*RR),&
        &   pflx(nr,0),pflx(nr,1)
    enddo
    close(890)
    deallocate(STR_OUTPUT_FILE_POWERBALANCE)
    
    if (T3D_POWERBALANCE_CALC_SHOTLOOP_flg ==1) return

    stop
  
    RETURN
END SUBROUTINE POWERBALANCE_CALC