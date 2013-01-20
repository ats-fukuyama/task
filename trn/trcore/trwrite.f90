MODULE trwrite

  USE trcomm, ONLY: ikind,rkind

  PUBLIC

CONTAINS

  SUBROUTINE tr_write_open(ierr)
    USE trcomm,ONLY: unitid,mdlwrt,kwpnam,kwtnam
    INTEGER(ikind),INTENT(OUT) :: ierr

    unitid(1) = 20 ! csv ouput of radial profiles
    unitid(2) = 21 ! csv ouput of time evolutions

    IF(mdlwrt /= 0)THEN
       OPEN(unitid(1),FILE=kwpnam,IOSTAT=ierr)
       OPEN(unitid(2),FILE=kwtnam,IOSTAT=ierr)
    END IF

    IF(ierr /= 0)THEN
       WRITE(6,*) 'XX tr_setup: File open error. IERR=',ierr
    END IF

    RETURN
  END SUBROUTINE tr_write_open

! ************************************************************************

  SUBROUTINE tr_write_close(ierr)
    USE trcomm, ONLY: unitid
    INTEGER(ikind),INTENT(OUT) :: ierr

    CLOSE(unitid(1),IOSTAT=ierr)
    CLOSE(unitid(2),IOSTAT=ierr)

    RETURN
  END SUBROUTINE tr_write_close

! ************************************************************************

  SUBROUTINE tr_writep_csv
    ! *** The csv output routine for radial profile ***
    ! *** The unit of the variables correspond to the list in trcomm.f90
    USE trcomm,ONLY: unitid,mdluf,nrmax,ntmax,neqmax,nsamax,  &
         ns_nsa,nsa_neq,nva_neq,kidns,idnsa,t,rhog,rhom,&
         rn,ru,rt,dtr,vtr,dtr_tb,vtr_tb,dtr_nc,vtr_nc,  &
         fluxtb,fluxnc,eta,qp,jtot,joh,jbs_nc,jex_nc,   &
         jcd_nb,jcd_ec,jcd_lh,jcd_ic,ptot,poh,pnb,pec,  &
         pibw,pic,plh,pnf,prl,pwl,snb,spl,swl,          &
         vtor,vpol,vpar,vprp,wrot,                      &
         pvolrho,psurrho,dvrho,rdpvrho,arrho,abb2rho,   &
         aib2rho,abvrho,ar1rho,ar2rho,abrho,rmjrho,     &
         rmnrho,rkprho,abb1rho,epsrho,bp,er,ezoh,z_eff, &
         mshear,mcurv,vexbp,dvexbpdr,wexbp,v_se,alpha,  &
         wp_t,wp_th,wp_inc,rw,taue1,taue2,taue3,        &
         taue89,taue98,h89,h98y2,betan,                 &
         pin_t,poh_t,pnb_t,pec_t,pibw_t,pic_t,          &
         plh_t,pnf_t,prl_t,pout_t,stdrt,offrt

    USE trufin,ONLY: rtug,rnug,qpug

    CHARACTER(LEN=15)  :: kwnsa
    INTEGER(ikind)     :: nr, ns, nsa, nwr, nwrmax

    CHARACTER(LEN=15),DIMENSION(100) :: kwcsv
    REAL(rkind),      DIMENSION(100) :: wcsv


    DO nr = 0, nrmax
       nwr = 1

       ! ---  substitution  ---
       wcsv(nwr) = rhog(nr) ; kwcsv(nwr) = 'RHOG' ; nwr = nwr + 1
       wcsv(nwr) = rhom(nr) ; kwcsv(nwr) = 'RHOM' ; nwr = nwr + 1
!
       DO nsa = 1, nsamax
          ns = ns_nsa(nsa)
          kwnsa = 'RN'//TRIM(kidns(ns))
          IF(idnsa(nsa)==2) kwnsa = TRIM(kwnsa)//'_f'
          wcsv(nwr) = rn(nsa,nr) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       END DO
!!$       DO nsa = 1, nsamax
!!$          ns = ns_nsa(nsa)
!!$          kwnsa = 'RU'//TRIM(kidns(ns))
!!$          IF(idnsa(nsa)==2) kwnsa = TRIM(kwnsa)//'_f'
!!$          wcsv(nwr) = ru(nsa,nr) ; kwscv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$       END DO
       DO nsa = 1, nsamax
          ns = ns_nsa(nsa)
          kwnsa = 'RT'//TRIM(kidns(ns))
          IF(idnsa(nsa)==2) kwnsa = TRIM(kwnsa)//'_f'
          wcsv(nwr) = rt(nsa,nr) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       END DO

       IF(mdluf > 0)THEN ! experimental data
          DO nsa = 1, nsamax
             ns = ns_nsa(nsa)
             kwnsa = 'RN_exp'//TRIM(kidns(ns))
             IF(idnsa(nsa)==2) kwnsa = TRIM(kwnsa)//'_f'
             wcsv(nwr) = rnug(nsa,nr) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
          END DO
!!$          DO nsa = 1, nsamax
!!$             ns = ns_nsa(nsa)
!!$             kwnsa = 'RU_exp'//TRIM(kidns(ns))
!!$             IF(idnsa(nsa)==2) kwnsa = TRIM(kwnsa)//'_f'
!!$             wcsv(nwr) = ru(nsa,nr) ; kwscv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$          END DO
          DO nsa = 1, nsamax
             ns = ns_nsa(nsa)
             kwnsa = 'RT_exp'//TRIM(kidns(ns))
             IF(idnsa(nsa)==2) kwnsa = TRIM(kwnsa)//'_f'
             wcsv(nwr) = rtug(nsa,nr) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
          END DO
       END IF

       DO nsa = 1, nsamax
          ns = ns_nsa(nsa)
          kwnsa = 'D'//TRIM(kidns(ns))
          IF(idnsa(nsa)==2 .OR. idnsa(nsa)==0) CYCLE
          wcsv(nwr) = dtr(1+3*nsa-2,1+3*nsa-2,nr) ; kwcsv(nwr) = TRIM(kwnsa)
          nwr = nwr + 1
       END DO
!!$       DO nsa = 1, nsamax
!!$          WRITE(knsa,'(I1)') nsa
!!$          kwnsa = 'D_u'//TRIM(knsa)
!!$          IF(idnsa(nsa)==2 .OR. idnsa(nsa)==0) CYCLE
!!$          wcsv(nwr) = dtr(1+3*nsa-1,1+3*nsa-1,nr) ; kwcsv(nwr) = TRIM(kwnsa)
!!$          nwr = nwr + 1
!!$       END DO
       DO nsa = 1, nsamax
          ns = ns_nsa(nsa)
          kwnsa = 'chi'//TRIM(kidns(ns))
          IF(idnsa(nsa)==2 .OR. idnsa(nsa)==0) CYCLE
          wcsv(nwr) = dtr(1+3*nsa  ,1+3*nsa  ,nr) ; kwcsv(nwr) = TRIM(kwnsa)
          nwr = nwr + 1
       END DO

       wcsv(nwr) = qp(nr)     ; kwcsv(nwr) = 'qp'    ; nwr = nwr + 1
       IF(mdluf > 0)THEN
          wcsv(nwr) = qpug(nr)  ; kwcsv(nwr) = 'qp_exp'  ; nwr = nwr + 1
       END IF
       wcsv(nwr) = jtot(nr)   ; kwcsv(nwr) = 'j_tot' ; nwr = nwr + 1
       wcsv(nwr) = joh(nr)    ; kwcsv(nwr) = 'j_oh'  ; nwr = nwr + 1
       wcsv(nwr) = jbs_nc(nr) ; kwcsv(nwr) = 'j_bs'  ; nwr = nwr + 1
       wcsv(nwr) = jcd_nb(nr) ; kwcsv(nwr) = 'j_nb'  ; nwr = nwr + 1
       wcsv(nwr) = jcd_ec(nr) ; kwcsv(nwr) = 'j_ec'  ; nwr = nwr + 1
       wcsv(nwr) = jcd_ic(nr) ; kwcsv(nwr) = 'j_ic'  ; nwr = nwr + 1
       wcsv(nwr) = jcd_lh(nr) ; kwcsv(nwr) = 'j_lh'  ; nwr = nwr + 1

       wcsv(nwr) = SUM(poh(:,nr))  ; kwcsv(nwr) = 'P_oh'  ; nwr = nwr + 1
       wcsv(nwr) = SUM(pnb(:,nr))  ; kwcsv(nwr) = 'P_nb'  ; nwr = nwr + 1
       wcsv(nwr) = SUM(pec(:,nr))  ; kwcsv(nwr) = 'P_ec'  ; nwr = nwr + 1
       wcsv(nwr) = SUM(pic(:,nr))  ; kwcsv(nwr) = 'P_ic'  ; nwr = nwr + 1
       wcsv(nwr) = SUM(plh(:,nr))  ; kwcsv(nwr) = 'P_lh'  ; nwr = nwr + 1
       wcsv(nwr) = SUM(prl(:,nr))  ; kwcsv(nwr) = 'P_rl'  ; nwr = nwr + 1
       wcsv(nwr) = SUM(pnf(:,nr))  ; kwcsv(nwr) = 'P_nf'  ; nwr = nwr + 1


       ! ---  write into the file in csv format  ---
       nwrmax = nwr - 1
       IF(nr == 0)THEN ! label
          DO nwr = 1, nwrmax
             kwcsv(nwr) = ADJUSTR(kwcsv(nwr))
             IF(nwr == 1) kwcsv(nwr)(1:1) = '#' ! the comment out identifier
             WRITE(unitid(1),'(A15,A1)',ADVANCE='NO') kwcsv(nwr),','
          END DO
          WRITE(unitid(1),*) ! line break
       END IF
       DO nwr = 1, nwrmax
          WRITE(unitid(1),'(PE15.6,A1)',ADVANCE='NO') wcsv(nwr),','
       END DO
       WRITE(unitid(1),*) ! line break

    END DO

    RETURN
  END SUBROUTINE tr_writep_csv


  SUBROUTINE tr_writet_csv
    ! *** The csv output routine for radial profile ***
    USE trcomm,ONLY: unitid,mdluf,nrmax,ntmax,neqmax,nsamax,    &
         ns_nsa,nsa_neq,nva_neq,kidns,idnsa,t,rhog,rhom,t,&
         rn,ru,rt,dtr,vtr,dtr_tb,vtr_tb,dtr_nc,vtr_nc,    &
         fluxtb,fluxnc,eta,qp,rip,jtot,joh,jbs_nc,jex_nc, &
         jcd_nb,jcd_ec,jcd_lh,jcd_ic,ptot,poh,pnb,pec,    &
         pibw,pic,plh,pnf,prl,pwl,snb,spl,swl,            &
         vtor,vpol,vpar,vprp,wrot,                        &
         pvolrho,psurrho,dvrho,rdpvrho,arrho,abb2rho,     &
         aib2rho,abvrho,ar1rho,ar2rho,abrho,rmjrho,       &
         rmnrho,rkprho,abb1rho,epsrho,bp,er,ezoh,z_eff,   &
         mshear,mcurv,vexbp,dvexbpdr,wexbp,v_se,alpha,    &
         ws_t,wp_t,wp_th,wp_inc,wpu_inc,rw,               &
         taue1,taue2,taue3,taue89,taue98,h89,h98y2,betan, &
         pin_t,poh_t,pnb_t,pec_t,pibw_t,pic_t,            &
         plh_t,pnf_t,prl_t,pout_t,stdrt,offrt

    USE trufin,ONLY: rtug,rnug,qpug,wthug,wtotug

    CHARACTER(LEN=15)  :: kwnsa
    INTEGER(ikind)     :: nr0,nr02,nr09, ns, nsa, nwr, nwrmax

    CHARACTER(LEN=15),DIMENSION(100) :: kwcsv
    REAL(rkind),      DIMENSION(100) :: wcsv

    nr0  = 0
    nr02 = INT(0.2*nrmax)
    nr09 = INT(0.9*nrmax)

    nwr = 1

    ! ---  substitution  ---
    wcsv(nwr) = t ; kwcsv(nwr) = 'TIME' ; nwr = nwr + 1

    DO nsa = 1, nsamax
       IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
       ns = ns_nsa(nsa)
       kwnsa = 'RN(0)'//TRIM(kidns(nsa))
       wcsv(nwr) = rn(nsa,nr0)  ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       kwnsa = 'RN(0.2)'//TRIM(kidns(nsa))
       wcsv(nwr) = rn(nsa,nr02) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       kwnsa = 'RN(0.9)'//TRIM(kidns(nsa))
       wcsv(nwr) = rn(nsa,nr09) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
    END DO
!!$    DO nsa = 1, nsamax
!!$       IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
!!$       ns = ns_nsa(nsa)
!!$       kwnsa = 'RU(0)'//TRIM(kidns(nsa))
!!$       wcsv(nwr) = ru(nsa,nr0)  ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$       kwnsa = 'RU(0.2)'//TRIM(kidns(nsa))
!!$       wcsv(nwr) = ru(nsa,nr02) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$       kwnsa = 'RU(0.9)'//TRIM(kidns(nsa))
!!$       wcsv(nwr) = ru(nsa,nr09) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$    END DO
    DO nsa = 1, nsamax
       IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
       ns = ns_nsa(nsa)
       kwnsa = 'RT(0)'//TRIM(kidns(nsa))
       wcsv(nwr) = rt(nsa,nr0)  ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       kwnsa = 'RT(0.2)'//TRIM(kidns(nsa))
       wcsv(nwr) = rt(nsa,nr02) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       kwnsa = 'RT(0.9)'//TRIM(kidns(nsa))
       wcsv(nwr) = rt(nsa,nr09) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
    END DO

    IF(mdluf > 0)THEN
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
          ns = ns_nsa(nsa)
          kwnsa = 'RN_exp(0)'//TRIM(kidns(nsa))
          wcsv(nwr) = rnug(nsa,nr0)  ; kwcsv(nwr) = TRIM(kwnsa); nwr = nwr + 1
          kwnsa = 'RN_exp(0.2)'//TRIM(kidns(nsa))
          wcsv(nwr) = rnug(nsa,nr02) ; kwcsv(nwr) = TRIM(kwnsa); nwr = nwr + 1
          kwnsa = 'RN_exp(0.9)'//TRIM(kidns(nsa))
          wcsv(nwr) = rnug(nsa,nr09) ; kwcsv(nwr) = TRIM(kwnsa); nwr = nwr + 1
       END DO
!!$    DO nsa = 1, nsamax
!!$       IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
!!$       ns = ns_nsa(nsa)
!!$       kwnsa = 'RU(0)'//TRIM(kidns(nsa))
!!$       wcsv(nwr) = ru(nsa,nr0)  ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$       kwnsa = 'RU(0.2)'//TRIM(kidns(nsa))
!!$       wcsv(nwr) = ru(nsa,nr02) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$       kwnsa = 'RU(0.9)'//TRIM(kidns(nsa))
!!$       wcsv(nwr) = ru(nsa,nr09) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
!!$    END DO
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
          ns = ns_nsa(nsa)
          kwnsa = 'RT_exp(0)'//TRIM(kidns(nsa))
          wcsv(nwr) = rtug(nsa,nr0)  ; kwcsv(nwr) = TRIM(kwnsa); nwr = nwr + 1
          kwnsa = 'RT_exp(0.2)'//TRIM(kidns(nsa))
          wcsv(nwr) = rtug(nsa,nr02) ; kwcsv(nwr) = TRIM(kwnsa); nwr = nwr + 1
          kwnsa = 'RT_exp(0.9)'//TRIM(kidns(nsa))
          wcsv(nwr) = rtug(nsa,nr09) ; kwcsv(nwr) = TRIM(kwnsa); nwr = nwr + 1
       END DO
    END IF


    wcsv(nwr) = qp(0)     ; kwcsv(nwr) = 'qp(0)' ; nwr = nwr + 1
    wcsv(nwr) = qp(nrmax) ; kwcsv(nwr) = 'qp(a)' ; nwr = nwr + 1
    IF(mdluf > 0)THEN
       wcsv(nwr) = qpug(0)     ; kwcsv(nwr) = 'qp_exp(0)' ; nwr = nwr + 1
       wcsv(nwr) = qpug(nrmax) ; kwcsv(nwr) = 'qp_exp(a)' ; nwr = nwr + 1
    END IF
    wcsv(nwr) = rip       ; kwcsv(nwr) = 'Ip'    ; nwr = nwr + 1

    wcsv(nwr) = wp_t   ; kwcsv(nwr) = 'Wp_tot'  ; nwr = nwr + 1
    wcsv(nwr) = wp_th  ; kwcsv(nwr) = 'Wp_th'   ; nwr = nwr + 1
    wcsv(nwr) = wp_inc ; kwcsv(nwr) = 'Wp_inc'  ; nwr = nwr + 1
    IF(mdluf > 0)THEN
       wcsv(nwr) = wtotug  ; kwcsv(nwr) = 'Wp_tot(exp)'  ; nwr = nwr + 1
       wcsv(nwr) = wthug   ; kwcsv(nwr) = 'Wp_th(exp)'   ; nwr = nwr + 1
       wcsv(nwr) = wpu_inc ; kwcsv(nwr) = 'Wp_inc(exp)'  ; nwr = nwr + 1
    END IF
    wcsv(nwr) = taue1  ; kwcsv(nwr) = 'tau_e1'  ; nwr = nwr + 1
    wcsv(nwr) = taue2  ; kwcsv(nwr) = 'tau_e2'  ; nwr = nwr + 1
    wcsv(nwr) = taue3  ; kwcsv(nwr) = 'tau_e3'  ; nwr = nwr + 1
    wcsv(nwr) = taue89 ; kwcsv(nwr) = 'tau_e89' ; nwr = nwr + 1
    wcsv(nwr) = taue98 ; kwcsv(nwr) = 'tau_e98' ; nwr = nwr + 1
    wcsv(nwr) = betan  ; kwcsv(nwr) = 'beta_n'  ; nwr = nwr + 1

    wcsv(nwr) = pin_t  ; kwcsv(nwr) = 'Pin'  ; nwr = nwr + 1
    wcsv(nwr) = poh_t  ; kwcsv(nwr) = 'Poh'  ; nwr = nwr + 1
    wcsv(nwr) = pnb_t  ; kwcsv(nwr) = 'Pnb'  ; nwr = nwr + 1
    wcsv(nwr) = pec_t  ; kwcsv(nwr) = 'Pec'  ; nwr = nwr + 1
    wcsv(nwr) = pic_t  ; kwcsv(nwr) = 'Pic'  ; nwr = nwr + 1
    wcsv(nwr) = plh_t  ; kwcsv(nwr) = 'Plh'  ; nwr = nwr + 1
    wcsv(nwr) = pnf_t  ; kwcsv(nwr) = 'Pnf'  ; nwr = nwr + 1
    wcsv(nwr) = prl_t  ; kwcsv(nwr) = 'Prl'  ; nwr = nwr + 1
    wcsv(nwr) = pout_t ; kwcsv(nwr) = 'Pout' ; nwr = nwr + 1

    ! comparison with experimental data
    IF(mdluf > 0)THEN
       wcsv(nwr) = rw      ; kwcsv(nwr) = 'RW'          ; nwr = nwr + 1
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==0 .OR. idnsa(nsa)==2) CYCLE
          ns = ns_nsa(nsa)
          kwnsa = 'STD_T'//TRIM(kidns(nsa))
          wcsv(nwr) = stdrt(nsa) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
          kwnsa = 'OFF_T'//TRIM(kidns(nsa))
          wcsv(nwr) = offrt(nsa) ; kwcsv(nwr) = TRIM(kwnsa) ; nwr = nwr + 1
       END DO
    END IF


    ! ---  write into the file in csv format  ---
    nwrmax = nwr - 1
    IF(t == 0.d0)THEN ! label
       DO nwr = 1, nwrmax
          kwcsv(nwr) = ADJUSTR(kwcsv(nwr))
          IF(nwr == 1) kwcsv(nwr)(1:1) = '#' ! the comment out identifier
          WRITE(unitid(2),'(A15,A1)',ADVANCE='NO') kwcsv(nwr),','
       END DO
       WRITE(unitid(2),*) ! line break
    END IF
    DO nwr = 1, nwrmax
       WRITE(unitid(2),'(PE15.6,A1)',ADVANCE='NO') wcsv(nwr),','
    END DO
    WRITE(unitid(2),*) ! line break

    RETURN
  END SUBROUTINE tr_writet_csv


  SUBROUTINE tr_write3d_csv

    RETURN
  END SUBROUTINE tr_write3d_csv

END MODULE trwrite
