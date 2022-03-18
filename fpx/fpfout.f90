! fpfout.f90

MODULE fpfout

  USE fpcomm,ONLY: rkind
  
  PRIVATE
  PUBLIC fp_fout,fp_fout_tot

  REAL(rkind),allocatable,dimension(:,:,:) :: &
       density, p_flux, Deff, temperature, Sr, Sr_Dp, Sr_Dth, Sr_Dr, Sr_F
  REAL(rkind),allocatable,dimension(:,:,:,:,:) :: &
       fI, source

CONTAINS

  SUBROUTINE fp_fout
    use fpcomm
    RETURN
  END SUBROUTINE fp_fout
    
  SUBROUTINE fp_fout_tot(nt)
    use fowcomm
    use fpcomm
    use fpwrite
    use foworbit
    use fowdistribution
    use fowevalNC
    use fowobclassify
    implicit none
    integer,intent(in) :: nt
    INTEGER,SAVE:: nthmax_save=0
    INTEGER,SAVE:: npmax_save=0
    INTEGER,SAVE:: nrmax_save=0
    INTEGER,SAVE:: nsamax_save=0
    INTEGER,SAVE:: ntmax_save=0
   
    REAL(rkind) :: r_, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0
    REAL(rkind) :: dpsipdr0, MJ2keV
    REAL(rkind),dimension(nrmax) :: rmg
    REAL(rkind),dimension(npmax,nsamax) :: momm
    REAL(rkind),dimension(npmax+1,nsamax) :: momg
    REAL(rkind),dimension(nthmax,npmax,nrmax,nsamax) :: taup, r0
    REAL(rkind),dimension(nrmax,nsamax) :: Drhmrhm,Drw,Drwav,nud_mono,Dbanana
    REAL(rkind),dimension(nrmax,nsamax) :: Dnewba,Dnewpla
    REAL(rkind),dimension(nrmax,nsamax) :: jaceffect,heatfow,heatnewba,heatrw
    integer :: nth,np,nr,nsa

    IF(nthmax.NE.nthmax_save.OR. &
       npmax.NE.npmax_save.OR. &
       nrmax.NE.nrmax_save.OR. &
       nsamax.NE.nsamax_save.OR. &
       ntmax.NE.ntmax_save) THEN
       IF(ALLOCATED(density)) THEN
          DEALLOCATE(density,p_flux,Sr,Sr_Dp,Sr_Dth,Sr_Dr,Sr_F,Deff)
          DEALLOCATE(fI,source,temperature)
       END IF
       ALLOCATE(p_flux (             nrmax,nsamax,ntmax))
       ALLOCATE(Sr     (             nrmax,nsamax,ntmax))
       ALLOCATE(Sr_Dp  (             nrmax,nsamax,ntmax))
       ALLOCATE(Sr_Dth (             nrmax,nsamax,ntmax))
       ALLOCATE(Sr_Dr  (             nrmax,nsamax,ntmax))
       ALLOCATE(Sr_F   (             nrmax,nsamax,ntmax))
       ALLOCATE(Deff   (             nrmax,nsamax,ntmax))
       ALLOCATE(fI     (nthmax,npmax,nrmax,nsamax,ntmax))
       ALLOCATE(source (nthmax,npmax,nrmax,nsamax,ntmax))
       ALLOCATE(temperature(         nrmax,nsamax,ntmax))
    END IF

    MJ2keV = 1.d-3/aee*1.d6

    if ( nt == 1 ) then
       do nr = 1, nrmax
          rmg(nr) = rg(nr+1)
       end do

       do nsa = 1, nsamax
          do np = 1, npmax
             momm(np,nsa) = pm(np,nsa)*ptfp0(nsa)
             momg(np,nsa) = pg(np,nsa)*ptfp0(nsa)
          end do
          momg(npmax+1,nsa) = pg(npmax+1,nsa)*ptfp0(nsa)
       end do
    
       do nsa = 1, nsamax
          do nr = 1, nrmax
             do np = 1, npmax
                do nth = 1, nthmax
                   call mean_ra_quantities(orbit_m(nth,np,nr,nsa), &
                        r_, psip0, costh0, sinth0, B0, F0, &
                        dBdr0, dFdr0, dpsipdr0)
                   r0(nth,np,nr,nsa) = r_
                   taup(nth,np,nr,nsa) &
                        =orbit_m(nth,np,nr,nsa) &
                        %time(orbit_m(nth,np,nr,nsa)%nstp_max)
                end do
             end do
          end do
       end do

       call fptxt1D(psimg,"dat/psimg.txt")
       call fptxt1D(Fpsig,"dat/Fpsig.txt")
       call fptxt1D(Boutg,"dat/Boutg.txt")
       call fptxt1D(Bing,"dat/Bing.txt")
       call fptxt1D(psim,"dat/psim.txt")
       call fptxt1D(Fpsi,"dat/Fpsi.txt")
       call fptxt1D(Bout,"dat/Bout.txt")
       call fptxt1D(Bin,"dat/Bin.txt")
       call fptxt1D(rm,"dat/rm.txt")
       call fptxt1D(rg,"dat/rg.txt")
       call fptxt1D(rmg,"dat/rmg.txt")

       call fptxt1D(safety_factor,"dat/safety_factor.txt")
    
       call fptxt2D(pm,"dat/pm.txt")
       call fptxt2D(pg,"dat/pg.txt")
       call fptxt2D(momm,"dat/momentum.txt")
       call fptxt2D(momg,"dat/momentumg.txt")
    
       call fptxt4D(thetam,"dat/thetam.txt")
       call fptxt4D(thetamg,"dat/thetamg.txt")
       call fptxt4D(thetam_pg,"dat/thetam_pg.txt")
       call fptxt4D(thetam_rg,"dat/thetam_rg.txt")
       call fptxt4D(JI,"dat/Jacobian.txt")
       call fptxt4D(JIR,"dat/JIR.txt")
       call fptxt4D(taup,"dat/taup.txt")
       call fptxt4D(r0,"dat/r_mean.txt")
    
       call fptxt4D(Dppfow,"dat/Dpp.txt")
       call fptxt4D(Dptfow,"dat/Dpt.txt")
       call fptxt4D(Dprfow,"dat/Dpr.txt")
    
       call fptxt4D(Dtpfow,"dat/Dtp.txt")
       call fptxt4D(Dttfow,"dat/Dtt.txt")
       call fptxt4D(Dtrfow,"dat/Dtr.txt")
    
       call fptxt4D(Drpfow,"dat/Drp.txt")
       call fptxt4D(Drtfow,"dat/Drt.txt")
       call fptxt4D(Drrfow,"dat/Drr.txt") 
      
       call fptxt4D(Fppfow,"dat/Fpp.txt")
       call fptxt4D(Fthfow,"dat/Fth.txt")
       call fptxt4D(Frrfow,"dat/Frr.txt")

       call fptxt4D(Dcpp,"dat/Dpp_fp.txt")
       call fptxt4D(Dcpt,"dat/Dpt_fp.txt")

       call fptxt4D(Dctp,"dat/Dtp_fp.txt")
       call fptxt4D(Dctt,"dat/Dtt_fp.txt")

       call fptxt4D(Fcpp,"dat/Fpp_fp.txt")
       call fptxt4D(Fcth,"dat/Fth_fp.txt")

       CALL fow_evaluate_NC(Drhmrhm,Drw,Drwav,Dbanana,nud_mono, &
                            heatrw,heatrwav)
       
       CALL fptxt2D(Drhmrhm,"dat/Drhmrhm.txt")
       CALL fptxt2D(heatfow, "dat/heatfow.txt")![2022/2/19] edited by anzai

       CALL fptxt2D(Drw,"dat/Drw.txt")
       CALL fptxt2D(Drwav,"dat/Drwav.txt")
!       CALL fptxt2D(nud_mono,"dat/nud_mono.txt")
!       CALL fptxt2D(Dbanana,"dat/Dbanana.txt")
       CALL fptxt2D(Dnewba, "dat/Dnewba.txt") ![2022/1/31] edited by anzai
       CALL fptxt2D(Dnewpla, "dat/Dnewpla.txt") ![2022/1/31] edited by anzai
!       CALL fptxt2D(jaceffect, "dat/jaceffect.txt")![2022/2/5] edited by anzai
       CALL fptxt2D(heatnewba, "dat/heatnewba.txt")![2022/2/19] edited by anzai
       CALL fptxt2D(heatrw, "dat/heatrw.txt")![2022/2/23] edited by anzai
       CALL fptxt2D(heatrwav, "dat/heatrwav.txt")![2022/2/23] edited by anzai

       CALL fow_ob_classify

       call fptxt1D(theta_m,"dat/thetam_obclass.txt")
       call fptxt1D(xi,"dat/xi_obclass.txt")
       call fptxt1D(rm,"dat/rm.txt")
       call fptxt3D(beta_D,"dat/beta_D_obclass.txt")
       call fptxt3D(beta_stag,"dat/beta_stag_obclass.txt")
       call fptxt3D(beta_pinch,"dat/beta_pinch_obclass.txt")
       call fptxt3D(theta_pinch,"dat/thetam_pinch_obclass.txt")
       call fptxt2D(xi_Xtype_boundary, &
            "dat/xi_Xtype_boundary_obclass.txt")
    end if

    do nsa = 1, nsamax
       do nr = 1, nrmax
          do np = 1, npmax
             do nth = 1, nthmax
                fI(nth,np,nr,nsa,nt)     = fnsp(nth,np,nr,nsa)
                source(nth,np,nr,nsa,nt) = SPP(nth,np,nr,nsa)
             end do
          end do
       end do
    end do

    call moment_0th_order_COM(density(:,:,nt), fnsp)
    call moment_2nd_order_COM(temperature(:,:,nt), fnsp)

    do nsa = 1, nsamax
       do nr = 1, nrmax
          temperature(nr,nsa,nt) &
               = temperature(nr,nsa,nt)*MJ2keV/(1.5d0*density(nr,nsa,nt)*1.d20)
       end do
    end do

    call particle_flux_element(Sr(:,:,nt), Sr_Dp(:,:,nt), Sr_Dth(:,:,nt), &
         Sr_Dr(:,:,nt), Sr_F(:,:,nt))
    call effective_diffusion_coefficient(Deff(:,:,nt))

    if ( nt == ntmax ) then
       call fptxt5D(fI,"dat/f.txt")
       call fptxt5D(source,"dat/source.txt")
       call fptxt3D(density,"dat/density.txt")
       call fptxt3D(Deff,"dat/Deff.txt")
       call fptxt3D(p_flux,"dat/p_flux.txt")
       call fptxt3D(Sr,"dat/Sr.txt")
       call fptxt3D(Sr_Dp,"dat/Sr_Dp.txt")
       call fptxt3D(Sr_Dth,"dat/Sr_Dth.txt")
       call fptxt3D(Sr_Dr,"dat/Sr_Dr.txt")
       call fptxt3D(Sr_F,"dat/Sr_F.txt")
       call fptxt3D(temperature,"dat/temperature.txt")
    end if

  end subroutine fp_fout_tot
END MODULE fpfout
