! fowexec.f90
! $Id: fpexec.f90,v 1.29 2013/02/08 07:36:24 nuga Exp $
! [2022/3/26]
! **************************
!    EXECUTE TIME ADVANCE
! **************************
! made by ota /  modified by anzai
!

module fowexec

  use fpcomm
  use fowcomm

contains

  subroutine fow_exec(nsa,IERR,its)
  !=========================================
  ! Main routine for time advance
  !=========================================
    use libmpi
    use libmtx
    use fpmpi

    implicit none
    real(8),dimension(nmend-nmstart+1):: Bm_L
    real(8),dimension(nthmax):: sendbuf_p, recvbuf_p
    real(8),dimension(nthmax*(npend-npstart+1)):: sendbuf_r, recvbuf_r
    double precision :: begin_time, end_time
    integer:: nsa, np, nth, nr, nl, nm, nn, ns
    integer:: nthS, nlL
    integer:: IERR,its,i,j,ll1
    integer:: imtxstart1,imtxend1
    integer :: nm_pnc, nm_D, nm_Xstg, nl_lost, nm_Ostg

    call cpu_time(begin_time)

    ns=ns_nsa(nsa)
    
    !**** allocate initial distribution
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          f(nth,np,nr) = fnsp(nth,np,nr,nsa)
        end do
      end do
    end do

    call mtx_set_communicator(comm_nrnp) !3D

    !**** Set up matrix solver 
    call mtx_setup(imtxsize,imtxstart1,imtxend1,imtxwidth)
    if(imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend) then
        write(6,*) 'XX fp_exec: '
        write(6,*) '   imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend'
        write(6,*) '   imtxstart1,imtxstart = ',imtxstart1,imtxstart
        write(6,*) '   imtxend1,imtxend     = ',imtxend1,imtxend
        stop
    end if

    !**** Set up weight array

    call fowweight(nsa,IERR)

    !**** Set up index array nma

    call SET_FM_NMA(nsa,fnsm)

    do nm = nmstart, nmend
        nlmax(nm)=0
        Bm(nm)=0.D0
        do nl = 1, nlmaxm
          LL(nm,nl)=0
          AL(nm,nl)=0.D0
        end do
    end do


    !**** Calculate matrix coefficients in a row 
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          nm = nma(nth,np,nr)
          call fowsetm(nth, np, nr, nsa, nlmax(nm))
        end do
      end do
    end do

    !**** internal boundary condition
    do nr = 1, nrmax
      do np = 1, npmax

        if ( theta_pnc(np,nr,nsa) /= NO_pinch_orbit ) then
          if ( pz(ns) >= 0.d0 ) then
            nm_D = nma(nth_pnc(nsa)-1,np,nr)
            nm_Xstg = nma(nth_stg(nsa),np,nr_pnc_point(np,nr,nsa))
          else
            nm_D = nma(nth_pnc(nsa),np,nr)
            nm_Xstg = nma(nth_stg(nsa)-1,np,nr_pnc_point(np,nr,nsa))
          end if
          call IBC_pinch(np,nr,nsa,nlmax(nm_D),nlmax(nm_Xstg))  
        end if


        if ( pz(ns) >= 0.d0 ) then
          nm_Xstg = nma(nth_stg(nsa),np,nr)
          if ( nr_rhom_pinch(np,nr,nsa) <= nrmax ) then
            nm_pnc = nma(nth_pnc(nsa)-1,np,nr_rhom_pinch(np,nr,nsa))
            call IBC_X_stagnation(np,nr,nsa,nlmax(nm_Xstg),nlmax(nm_pnc))
          else
            nl_lost = -1
            call IBC_X_stagnation(np,nr,nsa,nlmax(nm_Xstg),nl_lost)
          end if

        else
          nm_Xstg = nma(nth_stg(nsa)-1,np,nr)
          if ( nr_rhom_pinch(np,nr,nsa) <= nrmax ) then
            nm_pnc = nma(nth_pnc(nsa),np,nr_rhom_pinch(np,nr,nsa))
            call IBC_X_stagnation(np,nr,nsa,nlmax(nm_Xstg),nlmax(nm_pnc))
          else
            nl_lost = -1
            call IBC_X_stagnation(np,nr,nsa,nlmax(nm_Xstg),nl_lost)
          end if

        end if

        if ( pz(ns) >= 0.d0 ) then
          nm_Ostg = nma(nth_pnc(nsa)-1,np,nr)
        else
          nm_Ostg = nma(nth_pnc(nsa),np,nr)
        end if
        call IBC_O_stagnation(np,nr,nsa,nlmax(nm_Ostg))

      end do
    end do


    !**** Diagonal term 
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          nm = nma(nth,np,nr)
          Bm(nm) = (1.d0+(1.d0-rimpl)*delt*dl(nm))*fm(nm) &
                    +delt*spp(nth,np,nr,nsa)
          if ( nm.ge.imtxstart.and.nm.le.imtxend ) then
            call mtx_set_matrix(nm, nm, 1.d0-rimpl*delt*dl(nm))
            call mtx_set_vector(nm, fm(nm))
          end if
        end do
      end do
    end do

    !**** Off diagonal term 
    do nm = nmstart,nmend ! LHS
      if(nm.GE.imtxstart.AND.nm.LE.imtxend) then
        do nl=1,nlmax(nm)
            if(LL(nm,nl).NE.0) then
              call mtx_set_matrix(nm,LL(nm,nl),-RIMPL*DELT*AL(NM,NL))
            end if
        end do
      end if
    end do

    !**** Source vector: contribution from off-diagonal term 
    do nm=nmstart,nmend ! RHS
      do nl=1,nlmax(nm)
        nn=LL(nm,nl)
        if(nn.NE.0) then
          if(nn.ge.nmSTART-nthMAX.and.NN.le.nmend+NTHMAX)then
            Bm(nm)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,nl)*FM(nn)
          elseif(nn.lt.nmSTART-nthMAX)then
            Bm(nm)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,nl)*FM_shadow_m(nn)
          else
            Bm(nm)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,nl)*FM_shadow_p(nn)
          end if
        end if
      end do
      if(nm.GE.imtxstart.AND.nm.LE.imtxend) then
         call mtx_set_source(nm,Bm(nm))
      end if
    end do


    !**** Solve matrix equation 
    call mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC)
    !** ncom is nessesary for MUMPS not PETSc

    ierr=0

    !**** Get solution vector 
    call mtx_get_vector(Bm_L)

    do nr=nrstart, NRend
        do np=npstart, NPend
          do nth=1,NTHMAX
              nm=NMA(nth,np,nr)
              if(ABS(Bm_L(nm-nmstart+1)).LT.1.D-100) then
                Fns0(nth,np,nr,nsa)=0.D0
              else
                Fns0(nth,np,nr,nsa)=Bm_L(nm-nmstart+1)
              end if
          end do
        end do
    end do
    !** shadow requires to communicate

    call mtx_set_communicator(comm_np)
    do nr=nrstart, NRend
        call shadow_comm_np(nr,nsa)
    end do
    call mtx_set_communicator(comm_nr)
    call shadow_comm_nr(nsa)
    call mtx_set_communicator(comm_nrnp) !3D


    !**** Clean up matrix solver
    call mtx_cleanup

    call mtx_reset_communicator

    call cpu_time(end_time)
    write(6,'(A,I0,A,ES10.3,A)')'fowexec time (nsa=',nsa,') : ' &
                                   ,end_time-begin_time,'[sec]'

    return
  end subroutine fow_exec

!==================================================
!
!==================================================

  subroutine SET_FM_nma(nsa,func_in)

    implicit none
    integer:: nth, np, nr, nsa, nm, nrS, NPS
    double precision,dimension(nthmax,npstartw:npendwm&
                    ,nrstartw:nrendwm,nsastart:nsaend), &
                    intent(in):: func_in

    if(nrstart.eq.1)then
        nrs=1
    else
        nrs=nrstart-1
    end if

    do nr=nrstartw,nrendwm
        do np=npstartw,npendwm
          do nth=1,nthmax
              nm = nth+nthmax*(np-1)+npmax*nthmax*(nr-1)
              nma(nth,np,nr) = nm
          end do
        end do
    end do

    do nr=nrstart,nrend
        do np=npstartw,npendwm
          do nth=1,nthmax
              nm = nma(nth,np,nr)
              FM(nm) = func_in(nth,np,nr,nsa)
          end do
        end do
    end do 

    nr = nrstartw

    if(nr.ne.nrstart)then
        do np = npstartw, npendwm
          do nth=1,nthmax
              nm = nma(nth,np,nr)
              FM_shadow_m(nm)=func_in(nth,np,nr,nsa)
          end do
        end do
    end if

    nr = nrendwm

    if(nr.ne.nrend)then
        do np = npstartw, npendwm
          do nth = 1, nthmax
              nm = nma(nth,np,nr)
              FM_shadow_p(nm) = func_in(nth,np,nr,nsa)
          end do
        end do
    end if


  end subroutine SET_FM_nmA

  subroutine fowweight(nsa,IERR) 
  !------------------------------------------------
  !        CALCULATION OF weight
  ! proposed by Chang and Cooper [30] in Karney
  !------------------------------------------------

    use fpcomm
    use fowcomm

    implicit none
    real(8) :: dfdth, fvel, dfdp
    real(8) ::epswt=1.d-70
    double precision :: dfdrm, width_p, width_r, width_t, dth
    integer :: nrl, nrr, npl, npr
    integer :: nsa, np, nth, nr, nthl, nthr, ns
    integer :: ierr
    !**** calculation of weigthing (including off-diagonal terms)


    ns = ns_nsa(nsa)

    do nr = nrstart, nrend
      do np = npstart, npendwg
        do nth = 1, nthmax
            dfdth = 0.d0
            dfdrm  = 0.d0
            if ( np /= 1 ) then
              nthl = min(nth+1,nthmax)
              nthr = max(nth-1,1)
              nrl  = min(nr+1,nrmax)
              nrr  = max(nr-1,1)

              if ( np == npmax+1 ) then
                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                        / (2.d0*pg(np,ns) &
                        * delthm(nth,np-1,nr,nsa)*f(nth,np-1,nr))

                  dfdrm  = (f(nth,np-1,nrl)-f(nth,np-1,nrr)) &
                         / (2.d0*delr*f(nth,np-1,nr))
                end if
              else

                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  if ( abs(f(nth,np,nr)) > epswt ) then
                  
                    dth = (delthm(nth,np,nr,nsa) & 
                        + delthm(nth,np-1,nr,nsa))*0.5d0
                    if ( dth > 0.d0 ) then
                      dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                            / (4.d0*pg(np,ns)*dth*f(nth,np-1,nr)) &
                            + (f(nthl,np  ,nr)-f(nthr,np,nr)) &
                            / (4.d0*pg(np,ns)*dth*f(nth,np,nr))
                    end if

                    dfdrm  = (f(nth,np-1,nrl)-f(nthr,np-1,nrr)) &
                           / (4.d0*delr*f(nth,np-1,nr)) &
                           + (f(nthl,np  ,nrl)-f(nthr,np,nrr)) &
                           / (4.d0*delr*f(nth,np,nr))

                  else
                    dth = ( delthm(nth,np,nr,nsa) & 
                        + delthm(nth,np-1,nr,nsa) )*0.5d0
                    if ( dth > 0.d0 ) then
                      dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                            / (2.d0*pg(np,ns) & 
                            * delthm_pg(nth,np,nr,nsa)*f(nth,np-1,nr)) 
                    end if
                    
                    dfdrm  = (f(nth,np-1,nrl)-f(nth,np-1,nrr)) &
                           / (2.d0*delr*f(nth,np-1,nr)) 
                  end if
                else
                  if ( abs(f(nth,np,nr)) > epswt ) then
                    dth = (delthm(nth,np,nr,nsa) & 
                        + delthm(nth,np-1,nr,nsa))*0.5d0
                    if ( dth >= 0.d0 ) then
                      dfdth = (f(nthl,np,nr)-f(nthr,np,nr)) &
                            / (2.d0*pg(np,ns)*dth*f(nth,np  ,nr))
                    end if

                    dfdrm  = (f(nth,np,nrl)-f(nth,np,nrr)) &
                           / (2.d0*delr*f(nth,np,nr))
                  end if
                end if
              end if
            end if
            fvel = Fpp_j(nth,np,nr,nsa) & 
                 - Dpt_j(nth,np,nr,nsa)*dfdth & 
                 - Dpr_j(nth,np,nr,nsa)*dfdrm
            weighp(nth,np,nr,nsa) = fowwegh(-delp(ns)*fvel, & 
                                      Dpp_j(nth,np,nr,nsa))

        end do
      end do
    end do

    do nr = nrstart, nrend
      do np = npstart, npend
        do nth = 1, nthmax+1
          dfdp = 0.d0          
          npl = min(npmax,np+1)
          npr = max(1, np-1)
          width_p = dble(npl-npr)

          dfdrm = 0.d0
          nrl = min(nrmax,nr+1)
          nrr = max(1,nr-1)
          width_r = dble(nrl-nrr)
          
          if ( nth == 1 ) then
            if ( abs(f(nth,np,nr)) > epswt ) then
              dfdp = (f(nth,npl,nr)-f(nth,npr,nr)) &
                   / (width_p*delp(ns)*f(nth,np,nr))
              dfdrm = (f(nth,np,nrl)-f(nth,np,nrr)) &
                    / (width_r*delr   *f(nth,np,nr))  
            end if
          else if ( nth == nthmax+1 ) then
            if ( abs(f(nthmax,np,nr)) > epswt ) then
              dfdp = (f(nthmax,npl,nr)-f(nthmax,npr,nr)) &
                   / (width_p*delp(ns)*f(nthmax,np,nr))
              dfdrm = (f(nthmax,np,nrl)-f(nthmax,np,nrr)) & 
                    / (width_r*delr*f(nthmax,np,nr))  
            end if
          else
            if ( abs(f(nth,np,nr)) > epswt  & 
                   .and. abs(f(nth-1,np,nr)) > epswt ) then
              dfdp = ((f(nth-1,npl,nr)-f(nth-1,npr,nr)) & 
                   / (width_p*delp(ns)*f(nth-1,np,nr)) &
                   + (f(nth  ,npl,nr)-f(nth  ,npr,nr)) & 
                   / (width_p*delp(ns)*f(nth  ,np,nr)))/2.d0

              dfdrm = ((f(nth-1,np,nrl)-f(nth-1,np,nrr)) &
                    / (width_r*delr*f(nth-1,np,nr)) &
                    + (f(nth,np,nrl)-f(nth,np,nrr)) & 
                    / (width_r*delr*f(nth  ,np,nrr)))

            else if ( abs(f(nth,np,nr)) > epswt & 
                        .and. abs(f(nth-1,np,nr)) <= epswt ) then
              dfdp = (f(nth,npl,nr)-f(nth,npr,nr)) &
                   / (width_p*delp(ns)*f(nth,np,nr))
              dfdrm = (f(nth  ,np,nrl)-f(nth,np,nrr)) &
                    / (width_r*delr*f(nth  ,np,nr))

            else if (abs(f(nth,np,nr)) <= epswt  &
                       .and. abs(f(nth-1,np,nr)) > epswt ) then
              dfdp = (f(nth-1,npl,nr)-f(nth-1,npr,nr)) & 
                   / (width_p*delp(ns)*f(nth-1,np,nr))
              dfdrm = (f(nth-1,np,nrl)-f(nth-1,np,nrr)) & 
                    / (width_r*delr*f(nth-1,np,nr))
            end if
          end if

          fvel = Fth_j(nth,np,nr,nsa)-Dtp_j(nth,np,nr,nsa) & 
               * dfdp-Dtr_j(nth,np,nr,nsa)*dfdrm
          if ( nth <= nthmax ) then
            weight(nth,np,nr,nsa) = fowwegh(-delthm(nth,np,nr,nsa) & 
                                  * fvel,Dtt_j(nth,np,nr,nsa))
          else
            weight(nth,np,nr,nsa) = fowwegh(-delthm(nth-1,np,nr,nsa) & 
                                  * fvel,Dtt_j(nth,np,nr,nsa))
          end if
        end do
      end do
    end do

    do nr = nrstart, nrendwg
      do np = npstart, npend
        do nth = 1, nthmax
          dfdp = 0.d0          
          npl = min(npmax,np+1)
          npr = max(1, np-1)
          width_p = dble(npl-npr)

          dfdth = 0.d0
          nthl = min(nth+1,nthmax)
          nthr = max(nth-1,1)
          width_t = dble(nthl-nthr)

          if ( nr == 1 ) then
            if ( abs(f(nth,np,1)) > epswt ) then
              dfdp = (f(nth,npl,1)-f(nth,npr,1)) & 
                   / (width_p*delp(ns)*f(nth,np,1))
              if ( delthm_rg(nth,np,nr,nsa) /= 0.d0 ) then
                dfdth = (f(nthl,np,1)-f(nthr,np,1)) & 
                      / (width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,1))
              else
                dfdth = 0.d0
              end if
            end if
          else if ( nr == nrmax+1 ) then
            if ( abs(f(nth,np,nrmax)) > epswt ) then
              dfdp = (f(nth,npl,nrmax)-f(nth,npr,nrmax)) & 
                   / (width_p*delp(ns)*f(nth,np,nrmax))
              if ( delthm_rg(nth,np,nr,nsa) /= 0.d0 ) then
                dfdth= (f(nthl,np,nrmax)-f(nthr,np,nrmax)) & 
                     /(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nrmax))
              end if
            end if
          else
            if ( abs(f(nth,np,nr)) > epswt &
                   .and. abs(f(nth,np,nr-1)) > epswt ) then
              dfdp = ((f(nth,npl,nr-1)-f(nth,npr,nr-1)) & 
                   / (width_p*delp(ns)*f(nth,np,nr-1)) &
                   + (f(nth,npl,nr)-f(nth,npr,nr)) & 
                   / (width_p*delp(ns)*f(nth,np,nr)))/2.d0

              if ( nr == 2 ) then
                if ( delthm_rg(nth,np,nr,nsa) > 0.d0 ) then
                  dfdth = (f(nthl,np,nr)-f(nthr,np,nr)) &
                        / (width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nr))
                end if
              else
                if ( delthm_rg(nth,np,nr-1,nsa) > 0.d0  &
                       .and. delthm_rg(nth,np,nr  ,nsa) > 0.d0 ) then
                  if ( delthm_rg(nth,np,nr-1,nsa) > 0.d0 &
                         .and. delthm_rg(nth,np,nr  ,nsa) > 0.d0 ) then
                    dfdth = ((f(nthl,np,nr-1)-f(nthr,np,nr-1)) &
                          / (width_t*delthm_rg(nth,np,nr-1,nsa) & 
                          * f(nth,np,nr)) &
                          + (f(nthl,np,nr)-f(nthr,np,nr)) & 
                          / (width_t*delthm_rg(nth,np,nr,nsa) & 
                          * f(nth,np,nr)))/2.d0                  
                  end if
                end if
              end if

            else if ( abs(f(nth,np,nr)) > epswt  & 
                        .and. abs(f(nth,np,nr-1)) <= epswt ) then
              dfdp = (f(nth,npl,nr)-f(nth,npr,nr)) & 
                   / (width_p*delp(ns)*f(nth,np,nr))
              if ( delthm_rg(nth,np,nr,nsa) > 0.d0 ) then
                dfdth = (f(nthl,np,nr)-f(nthr,np,nr)) & 
                      / (width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nr))    
              end if

            else if ( abs(f(nth,np,nr)) <= epswt  & 
                        .and. abs(f(nth,np,nr-1)) > epswt ) then
              dfdp = (f(nth,npl,nr-1)-f(nth,npr,nr-1)) & 
                   / (width_p*delp(ns)*f(nth,np,nr-1))
              if ( delthm_rg(nth,np,nr-1,nsa) > 0.d0  ) then
                dfdth = (f(nthl,np,nr-1)-f(nthr,np,nr-1)) & 
                      / (width_t*delthm_rg(nth,np,nr-1,nsa)*f(nth,np,nr-1))
              end if

            end if
          end if
          fvel = Frr_j(nth,np,nr,nsa)-Drp_j(nth,np,nr,nsa) & 
               * dfdp-Drt_j(nth,np,nr,nsa)*dfdth
          weighr(nth,np,nr,nsa) = fowwegh(-delr*fvel,Drr_j(nth,np,nr,nsa))
        end do
      end do
    end do

    return
  end subroutine fowweight

  function fowwegh(X,Y)
  !-------------------------------------------------
  ! weighting function for convection effect
  !-------------------------------------------------

    implicit none
    real(8):: X, Y, Z
    real(8):: fowwegh

    if(abs(Y).LT.1.D-70) then
        if(X.GT.0.D0) then
          fowwegh=0.D0
        elseif(X.LT.0.D0) then
          fowwegh=1.D0
        else
          fowwegh=0.5D0
        end if
    else
        Z=X/Y
        if(abs(Z).LT.1.D-5)then
          fowwegh=0.5D0-Z/12.D0+Z**3/720.D0
        else if(Z.GE.100.D0)then
          fowwegh=1.D0/Z
        else if(Z.LE.-100.D0)then
          fowwegh=1.D0/Z+1.D0
        else
          fowwegh=1.D0/Z-1.D0/(EXP(Z)-1.D0)
        end if
    end if
    return
  end function fowwegh
  
  subroutine fowsetm(nth,np,nr,nsa,nl)
  !------------------------------------------
  !     Calculation of matrix coefficients
  !------------------------------------------
    implicit none
    integer,intent(in):: nth,np,nr,nsa
    integer,intent(out):: nl
    double precision:: pl, sl, rl
    double precision :: deltath
    double precision :: DIVD(3,3), DIVF(3), del(3), Ffow(2,3)
    double precision :: D_term, F_term
    integer :: si, sj, sk, sign_to_index(-1:1)
    integer :: alpha, beta, gama, loc(4), loc_pnc(4)
    integer :: nm
    integer :: ierr, ns
    integer :: boundary_flag


    ns=ns_nsa(nsa)

    nl=0
    nm=nma(nth,np,nr)

    Ffow(1,1) = Fpp_j(nth,np,nr,nsa)
    Ffow(2,1) = Fpp_j(nth,np+1,nr,nsa)

    Ffow(1,2) = Fth_j(nth,np,nr,nsa)
    Ffow(2,2) = Fth_j(nth+1,np,nr,nsa)

    Ffow(1,3) = Frr_j(nth,np,nr,nsa)
    Ffow(2,3) = Frr_j(nth,np,nr+1,nsa)


    !**** discretized (div(d/dX))_Y
    pl = pm(np,nsa)
    sl = sin( thetam(nth,np,nr,nsa) )
    rl = rm(nr)
    deltath= delthm(nth,np,nr,nsa)

    DIVD(1,1) = 1.d0/(delp(nsa)**2)
    DIVD(1,2) = 1.d0/(delp(nsa) * deltath * 2.d0)
    DIVD(1,3) = 1.d0/(delp(nsa) * delr * 2.d0)
    DIVD(2,1) = 1.d0/(delp(nsa) * deltath * 2.d0)
    DIVD(2,2) = 1.d0/(deltath**2 )
    DIVD(2,3) = 1.d0/(deltath * delr * 2.d0)
    DIVD(3,1) = 1.d0/(delp(nsa)  * delr * 2.d0)
    DIVD(3,2) = 1.d0/(deltath * delr * 2.d0)
    DIVD(3,3) = 1.d0/(delr**2)
    DIVF(1)   = 1.d0/(delp(nsa))
    DIVF(2)   = 1.d0/(deltath)
    DIVF(3)   = 1.d0/(delr)

    sign_to_index(-1) = 1
    sign_to_index(1)  = 2
    sign_to_index(0)  = 0
    loc = [nth,np,nr,nsa]

    !**** term of f(nth, np, nr)
    DL(nm) = 0.d0
    do alpha = 1, 3
      do si = -1, 1, 2
        D_term = -1.d0*Dfow(alpha,alpha,sign_to_index(si)-1,0,loc) & 
               * DIVD(alpha,alpha)
        F_term = -1.d0*si*Ffow(sign_to_index(si),alpha) & 
               * w(si,alpha,0,sign_to_index(si)-1,0,loc)*DIVF(alpha)
        dl(nm) = dl(nm) + (D_term + F_term)
      end do
    end do

    if ( nth == nth_pnc(nsa) .and. & 
           theta_pnc(np,nr,nsa) /= NO_PinCH_ORBIT ) then

      D_term = Dfow(2,2,0,0,loc)*DIVD(2,2)*IBCflux_ratio(np,nr,nsa)
      F_term = -1.d0*Ffow(1,2)*w(-1,2,0,0,0,loc) & 
             * DIVF(2)*IBCflux_ratio(np,nr,nsa)
      dl(nm) = dl(nm) + (D_term + F_term)

    end if

    !**** terms of f(i+si, j, k) 
    do alpha = 1, 3
      do si = -1, 1, 2
        boundary_flag = check_external_boundary(alpha,0,si,0,loc)

        if ( boundary_flag == 1  ) cycle
        D_term = Dfow(alpha,alpha,sign_to_index(si)-1,0,loc) & 
               * DIVD(alpha,alpha)
        F_term = -1.d0*si*Ffow(sign_to_index(si),alpha) & 
               * w(-1*si,alpha,0,sign_to_index(si)-1,0,loc)*DIVF(alpha)

        do beta = 1, 3
          if ( alpha == beta ) cycle
          do sj = -1, 1, 2
            D_term = D_term + si * sj & 
                   * Dfow(beta,alpha,sign_to_index(sj)-1,0,loc) &
                   * w(sj,beta,alpha,sign_to_index(sj)-1,si,loc) & 
                   * DIVD(beta,alpha)
          end do
        end do

        nl = nl+1
        ll(nm,nl) = get_nma(alpha,0,si,0,loc)
        al(nm,nl) = al(nm,nl) + (D_term + F_term)

        if ( abs(al(nm,nl)) < 1.d-70 ) then
          ll(nm,nl) = 0
          al(nm,nl) = 0.d0
          nl = nl-1
        end if

      end do
    end do

    !**** terms of f(i+si, j+sj, k) 
    do alpha = 1, 3
      do beta = 1, 3
        if ( alpha >= beta) cycle
        do si = -1, 1, 2
          do sj = -1, 1, 2

            boundary_flag = check_external_boundary(alpha,beta,si,sj,loc)
            if ( boundary_flag == 1 ) cycle
            
            D_term = si*sj*(Dfow(alpha,beta,sign_to_index(si)-1,0,loc) & 
                   * DIVD(alpha,beta) &
                   * w(-1*si,alpha,beta,sign_to_index(si)-1,sj,loc) &
                   + Dfow(beta,alpha,sign_to_index(sj)-1,0,loc) & 
                   * DIVD(beta,alpha) &
                   * w(-1*sj,beta,alpha,sign_to_index(sj)-1,si,loc))

            nl = nl+1
            ll(nm,nl) = get_nma(alpha,beta,si,sj,loc)
            al(nm,nl) = al(nm,nl) + D_term

            if ( abs(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    

          end do
        end do
      end do
    end do

    spp(nth,np,nr,nsa) = sppb(nth,np,nr,nsa)

    return
  end subroutine fowsetm

  function Dfow(alpha,beta,si,sj,loc)
  !-------------------------------------------
  !
  !------------------------------------------

    implicit none
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    double precision :: Dfow
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    if ( alpha == 1 ) then
      if ( beta == 1 ) then
        Dfow = Dpp_j(nth,np+si,nr,nsa)
      else if ( beta == 2 ) then
        Dfow = Dpt_j(nth+sj,np+si,nr,nsa)
      else if ( beta == 3 ) then
        Dfow = Dpr_j(nth,np+si,nr+sj,nsa)
      end if

    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        Dfow = Dtp_j(nth+si,np+sj,nr,nsa)
      else if ( beta == 2 ) then
        Dfow = Dtt_j(nth+si,np,nr,nsa)
      else if ( beta == 3 ) then
        Dfow = Dtr_j(nth+si,np,nr+sj,nsa)
      end if

    else if ( alpha == 3 ) then
      if ( beta == 1 ) then
        Dfow = Drp_j(nth,np+sj,nr+si,nsa)
      else if ( beta == 2 ) then
        Dfow = Drt_j(nth+sj,np,nr+si,nsa)
      else if ( beta == 3 ) then
        Dfow = Drr_j(nth,np,nr+si,nsa)
      end if

    end if

  end function Dfow

  function w(sign,alpha,beta,si,sj,loc)
  !----------------------------------------------
  ! weight function for alpha, beta, gamma
  !----------------------------------------------

    implicit none
    integer,intent(in) :: sign,alpha,beta,si,sj,loc(4)
    double precision :: w
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)
    
    if ( alpha == 1 ) then
      if ( beta == 2 ) then
        w = weighp(nth+sj, np+si, nr, nsa)
      else if ( beta == 3 ) then
        w = weighp(nth, np+si, nr+sj, nsa)
      else
        w = weighp(nth, np+si, nr, nsa)
      end if

    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        w = weight(nth+si, np+sj, nr, nsa)
      else if ( beta == 3 ) then
        w = weight(nth+si, np, nr+sj, nsa)
      else 
        w = weight(nth+si, np, nr, nsa)
      end if

    else
      if ( beta == 1 ) then
        w = weighr(nth, np+sj, nr+si, nsa)
      else if ( beta == 2 ) then
        w = weighr(nth+sj, np, nr+si, nsa)
      else
        w = weighr(nth, np, nr+si, nsa)
      end if
    end if

    if ( sign < 0 ) w = 1.d0-w

  end function w

  function check_external_boundary(alpha,beta,si,sj,loc) result(flag)
  !-------------------------------------------------------
  ! flag == 0 -> inside external boundary
  !         1 -> outside external boundary
  !-------------------------------------------------------

    implicit none
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: flag
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    flag = 0

    !**** check external boundary
    if ( alpha == 1 ) then
      if ( beta == 2 ) then
        if ( (np+si <= 0 .or. npmax+1 <= np+si)&
            .or. (nth+sj <= 0 .or. nthmax+1 <= nth+sj) )&
            flag = 1

      else if ( beta == 3 ) then
        if ( (np+si <= 0 .or. npmax+1 <= np+si)&
            .or. (nr+sj <= 0 .or. nrmax+1 <= nr+sj) )&
          flag = 1

      else !** beta == 0
        if ( np+si <= 0 .or. npmax+1 <= np+si ) flag = 1

      end if

    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        if ( (np+sj <= 0 .or. npmax+1 <= np+sj)&
            .or. (nth+si <= 0 .or. nthmax+1 <= nth+si) )&
          flag = 1

      else if ( beta == 3 ) then
        if ( (nth+si <= 0 .or. nthmax+1 <= nth+si)&
            .or. (nr+sj <= 0 .or. nrmax+1 <= nr+sj) )&
          flag = 1

      else !** beta == 0
        if ( nth+si <= 0 .or. nthmax+1 <= nth+si ) flag = 1

      end if

    else
      if ( beta == 1 ) then
        if ( (np+sj <= 0 .or. npmax+1 <= np+sj)&
            .or. (nr+si <= 0 .or. nrmax+1 <= nr+si) )&
          flag = 1

      else if ( beta == 2 ) then
          if ( (nth+sj <= 0 .or. nthmax+1 <= nth+sj)&
            .or. (nr+si <= 0 .or. nrmax+1 <= nr+si) )&
          flag = 1

      else !** beta == 0
        if ( nr+si <= 0 .or. nrmax+1 <= nr+si ) flag = 1

      end if
    end if

    !** check internal boundary
   
  end function check_external_boundary

  function get_nma(alpha,beta,si,sj,loc) result(n)
  !-----------------------------------------------
  !
  !-----------------------------------------------
 
    implicit none
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: n
    integer :: nth, np, nr, nsa

    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    if ( alpha == 1 ) then
      if ( beta == 2 ) then
        n = nma(nth+sj,np+si,nr)
      else if ( beta == 3 ) then
        n = nma(nth,np+si,nr+sj)
      else !** beta == 0
        n = nma(nth,np+si,nr)
      end if
    
    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        n = nma(nth+si,np+sj,nr)
      else if ( beta == 3 ) then
        n = nma(nth+si,np,nr+sj)
      else !** beta == 0
        n = nma(nth+si,np,nr)
      end if
    
    else if ( alpha == 3 ) then
      if ( beta == 1 ) then
        n = nma(nth,np+sj,nr+si)
      else if ( beta == 2 ) then
        n = nma(nth+sj,np,nr+si)
      else !** beta == 0
        n = nma(nth,np,nr+si)
      end if

    end if

  end function get_nma

!==========================================================
!
!=========================================================

  subroutine IBC_pinch(np,nr,nsa,nl_D,nl_Xstg)
  !------------------------------------------------
  !
  !-----------------------------------------------

    implicit none
    integer,intent(in) :: np,nr,nsa
    integer,intent(inout) :: nl_D, nl_Xstg
    double precision :: S_pnc_to_trp, S_pnc_to_pas
    double precision:: pl, sl, rl
    double precision :: deltath, ratio, fact1, fact2
    double precision :: DIVD(3), DIVF, Fthm, Fthp, Dtx(3)
    double precision :: D_term, D_tp, D_tt, D_tr, F_term 
    double precision :: dfnspdp, dfnspdth, dfnspdr, fnspg
    integer:: ierr, ns, nth, nm_D, nm_Xstg, nl
    integer :: loc_pnc(4), loc_Xstg(4), loc_D(4), loc(4)

    if ( theta_pnc(np,nr,nsa) == NO_PinCH_ORBIT ) then
      return
    end if

    ns = ns_nsa(nsa)

    nth = nth_pnc(nsa)

    pl = pm(np,nsa)
    sl = Sin( thetam(nth-1,np,nr,nsa) )
    rl = rm(nr)
    deltath= delthm(nth-1,np,nr,nsa)

    Dtx(1)  = Dtp_j(nth,np,nr,nsa)/(delp(nsa)*deltath*2.d0)
    Dtx(2)  = Dtt_j(nth,np,nr,nsa)/(deltath**2)
    Dtx(3)  = Dtr_j(nth,np,nr,nsa)/(deltath*delr*2.d0)
    Fthp    = Fth_j(nth,np,nr,nsa)/deltath !**

    loc_pnc = [nth,np,nr,nsa]
    loc_D   = [nth,np,nr,nsa]
    loc_Xstg= [nth_stg(nsa),np,nr_pnc_point(np,nr,nsa),nsa]

    if ( pz(ns) >= 0.d0 ) then
      !** calcurate sign of S_pnc_to_trp
      if ( np == npmax ) then
        dfnspdp = (f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                / ( delp(ns))
      else if ( np == 1 ) then
        dfnspdp = (f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                / (delp(ns))
      else
        dfnspdp = ( f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                / ( 2.d0*delp(ns) )
      end if

      if ( nr == nrmax ) then
        dfnspdr = (f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                / ( delr )  
      else if ( nr == 1 ) then
        dfnspdr = (f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                / ( delr )
      else
        dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                / ( 2.d0*delr )
      end if

      dfnspdth = ( fnsp(nth,np,nr,nsa)-fnsp(nth-1,np,nr,nsa) )&
               / (0.5d0*pl*(delthm(nth,np,nr,nsa) & 
               + delthm(nth-1,np,nr,nsa)))

      fnspg = f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0)

      S_pnc_to_trp = -1.d0*(Dtpfow(nth,np,nr,nsa)*dfnspdp  &
                   + Dttfow(nth,np,nr,nsa)*dfnspdth &
                   + Dtrfow(nth,np,nr,nsa)*dfnspdr)  &
                   + Fthfow(nth,np,nr,nsa)*fnspg

      if ( S_pnc_to_trp >= 0.d0 ) then
        return
      end if


      fact1 = (1.d0-IBCflux_ratio(np,nr,nsa))/JI(nth-1,np,nr,nsa)
      fact2 = IBCflux_ratio(np,nr,nsa)/JI(nth_stg(nsa), &
                              np,nr_pnc_point(np,nr,nsa),nsa)
      nm_D    = nma(nth-1,np,nr)
      nm_Xstg = nma(nth_stg(nsa), np, nr_pnc_point(np,nr,nsa))
      dl(nm_D) = dl(nm_D) + fact1*( Dtx(2)+Fthp*weight(nth,np,nr,nsa) )

      do nl = 1, nl_D
        if ( ll(nm_D,nl) == nma_boundary(nth,np+1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(1) &
                      * (1.d0-weight(nth,np+1,nr,nsa))*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np+1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(1)*weight(nth,np+1,nr,nsa)*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np-1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(1) & 
                      * (1.d0-weight(nth,np-1,nr,nsa))*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np-1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(1)*weight(nth,np-1,nr,nsa)*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) -( Dtx(2) &
                      - Fthp*(1.d0-weight(nth,np,nr,nsa)) )*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np,nr+1) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(3) & 
                      * (1.d0-weight(nth,np,nr+1,nsa))*fact1
          
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np,nr+1) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(3)*weight(nth,np,nr+1,nsa)*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np,nr-1) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(3) & 
                      * (1.d0-weight(nth,np,nr-1,nsa))*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np,nr-1) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(3)*weight(nth,np,nr-1,nsa)*fact1
  
        end if
      end do
  
      nl_Xstg = nl_Xstg+1 
      ll(nm_Xstg,nl_Xstg) = nma(nth-1,np,nr)
      al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + fact2 & 
                          * (Dtx(2)+Fthp*weight(nth,np,nr,nsa))
  
      if ( nma_boundary(nth,np+1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np+1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(1) & 
                            * (1.d0-weight(nth,np+1,nr,nsa))*fact2  
      end if
  
      if ( nma_boundary(nth-1,np+1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np+1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(1) & 
                            * weight(nth,np+1,nr,nsa)*fact2  
      end if
  
      if ( nma_boundary(nth,np-1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np-1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(1) & 
                            * (1.d0-weight(nth,np-1,nr,nsa))*fact2  
      end if
  
      if (  nma_boundary(nth-1,np-1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np-1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(1) & 
                            * weight(nth,np-1,nr,nsa)*fact2  
      end if
  
      if ( nma_boundary(nth,np,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) -( Dtx(2) &
                            + Fthp*(1.d0-weight(nth,np,nr,nsa)) )*fact2  
      end if
  
      if ( nma_boundary(nth,np,nr+1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np,nr+1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(3) & 
                            * (1.d0-weight(nth,np,nr+1,nsa))*fact2  
      end if
  
      if ( nma_boundary(nth-1,np,nr+1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np,nr+1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(3) & 
                            * weight(nth,np,nr+1,nsa)*fact2  
      end if
  
      if ( nma_boundary(nth,np,nr-1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np,nr-1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(3) & 
                            * (1.d0-weight(nth,np,nr-1,nsa))*fact2  
      end if
  
      if ( nma_boundary(nth-1,np,nr-1)  /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np,nr-1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(3) & 
                            * weight(nth,np,nr-1,nsa)*fact2  
      end if
  
    else !** PZ(ns) < 0
      if ( np == npmax ) then
        dfnspdp = (f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) &
                - f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                / ( delp(ns) )
      else if ( np == 1 ) then
        dfnspdp = (f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                / ( delp(ns) )
      else
        dfnspdp = (f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                / ( 2.d0*delp(ns) )
      end if

      if ( nr == nrmax ) then
        dfnspdr = (f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                / ( delr )  
      else if ( nr == 1 ) then
        dfnspdr = (f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                / ( delr )
      else
        dfnspdr = (f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                / ( 2.d0*delr )
      end if

      dfnspdth = (fnsp(nth,np,nr,nsa)-fnsp(nth-1,np,nr,nsa))&
               / (0.5d0*pl*(delthm(nth,np,nr,nsa) &
               + delthm(nth-1,np,nr,nsa)) )

      fnspg   = f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0)

      S_pnc_to_pas = -1.d0*(Dtpfow(nth,np,nr,nsa)*dfnspdp  &
                   + Dttfow(nth,np,nr,nsa)*dfnspdth &
                   + Dtrfow(nth,np,nr,nsa)*dfnspdr)  &
                   + Fthfow(nth,np,nr,nsa)*fnspg

      if ( S_pnc_to_pas <= 0.d0 ) then
        return
      end if

      fact1 = -1.d0*IBCflux_ratio(np,nr,nsa)/JI(nth,np,nr,nsa)
      fact2 = -1.d0*(1.d0-IBCflux_ratio(np,nr,nsa)) &
            / JI(nth_stg(nsa)-1,np,nr_pnc_point(np,nr,nsa),nsa)
      nm_D    = nma(nth,np,nr)
      nm_Xstg = nma(nth_stg(nsa)-1, np, nr_pnc_point(np,nr,nsa))

      dl(nm_D) = dl(nm_D) + fact1*( -1.d0 * Dtx(2) &
               + Fthp*(1.d0-weight(nth,np,nr,nsa)) )

      do nl = 1, nl_D
        if ( ll(nm_D,nl) == nma_boundary(nth,np+1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(1) & 
                      * (1.d0-weight(nth,np+1,nr,nsa))*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np+1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(1) & 
                      * weight(nth,np+1,nr,nsa)*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np-1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(1) & 
                      * (1.d0-weight(nth,np-1,nr,nsa))*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np-1,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(1) & 
                      * weight(nth,np-1,nr,nsa)*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np,nr) ) then
          al(nm_D,nl) = al(nm_D,nl) +( Dtx(2) &
                      + Fthp*weight(nth,np,nr,nsa) )*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np,nr+1) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(3) & 
                      * (1.d0-weight(nth,np,nr+1,nsa))*fact1
          
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np,nr+1) ) then
          al(nm_D,nl) = al(nm_D,nl) - Dtx(3)*weight(nth,np,nr+1,nsa)*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth,np,nr-1) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(3) & 
                      * (1.d0-weight(nth,np,nr-1,nsa))*fact1
  
        else if ( ll(nm_D,nl) == nma_boundary(nth-1,np,nr-1) ) then
          al(nm_D,nl) = al(nm_D,nl) + Dtx(3)*weight(nth,np,nr-1,nsa)*fact1
  
        end if
      end do
  
      nl_Xstg = nl_Xstg+1 
      ll(nm_Xstg,nl_Xstg) = nm_D
      al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + fact2 & 
                          * (-1.d0*Dtx(2) + Fthp &
                          * (1.d0-weight(nth,np,nr,nsa)) )
  
      if ( nma_boundary(nth,np+1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np+1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(1) & 
                            * (1.d0-weight(nth,np+1,nr,nsa))*fact2  
      end if
  
      if ( nma_boundary(nth-1,np+1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np+1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(1) & 
                            * weight(nth,np+1,nr,nsa)*fact2  
      end if
  
      if ( nma_boundary(nth,np-1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np-1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(1) & 
                            * (1.d0-weight(nth,np-1,nr,nsa))*fact2  
      end if
  
      if (  nma_boundary(nth-1,np-1,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np-1,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(1) & 
                            * weight(nth,np-1,nr,nsa)*fact2  
      end if
  
      if ( nma_boundary(nth,np,nr) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np,nr)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) -( Dtx(2) &
                            - Fthp*weight(nth,np,nr,nsa) )*fact2  
      end if
  
      if ( nma_boundary(nth,np,nr+1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np,nr+1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(3) & 
                            * (1.d0-weight(nth,np,nr+1,nsa))*fact2  
      end if
  
      if ( nma_boundary(nth-1,np,nr+1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np,nr+1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) - Dtx(3) & 
                            * weight(nth,np,nr+1,nsa)*fact2  
      end if
  
      if ( nma_boundary(nth,np,nr-1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth,np,nr-1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(3) & 
                            * (1.d0-weight(nth,np,nr-1,nsa))*fact2  
      end if
  
      if ( nma_boundary(nth-1,np,nr-1) /= -1 ) then
        nl_Xstg = nl_Xstg+1
        ll(nm_Xstg,nl_Xstg) = nma(nth-1,np,nr-1)
        al(nm_Xstg,nl_Xstg) = al(nm_Xstg,nl_Xstg) + Dtx(3) & 
                            * weight(nth,np,nr-1,nsa)*fact2  
      end if  

    end if

  end subroutine IBC_pinch

  subroutine IBC_X_stagnation(np,nr,nsa,nl_Xstg,nl_pnc)
  !----------------------------------------------------
  !
  !----------------------------------------------------

    implicit none
    integer,intent(in) :: np,nr,nsa
    integer,intent(inout) :: nl_Xstg, nl_pnc
    double precision:: pl, sl, rl
    double precision :: deltath, S_x, fact1, fact2
    double precision :: DIVD(3), DIVF, Fthm, Fthp, Dtx(3)
    double precision :: D_tp, D_tt, D_tr
    double precision :: dfnspdp, dfnspdth, dfnspdr, fnspg
    integer:: ierr, ns, nth, nm_Xstg, nm_pnc, nl
    integer :: loc_Xstg(4), loc_pnc(4), loc(4)

    ns = ns_nsa(nsa)

    if ( pz(ns) >= 0.d0 ) then
      
      nth = nth_stg(nsa)
      loc_Xstg = [nth,np,nr,nsa]
      if ( nl_pnc /= -1 ) then
        loc_pnc = [nth_pnc(nsa),np,nr_rhom_pinch(np,nr,nsa),nsa]
      end if

      pl = pm(np,nsa)
      sl = Sin( thetam(nth,np,nr,nsa) )
      rl = rm(nr)
      deltath= delthm(nth,np,nr,nsa)
  
      Dtx(1)  = Dtp_j(nth,np,nr,nsa)/(delp(nsa)*deltath*2.d0)
      Dtx(2)  = Dtt_j(nth,np,nr,nsa)/(deltath**2)
      Dtx(3)  = Dtr_j(nth,np,nr,nsa)/(deltath*delr*2.d0)
      Fthp    = Fth_j(nth,np,nr,nsa)/deltath !**  

      if ( np == npmax ) then
        dfnspdp = (f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) &
                - f_grid(loc_Xstg, 0.d0, -0.5d0, 0.5d0) )&
                / (delp(ns))
      else if ( np == 1 ) then
        dfnspdp = (f_grid(loc_Xstg, 0.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) )&
                / (delp(ns))
      else
        dfnspdp = (f_grid(loc_Xstg, 0.d0, 1.5d0, 0.5d0) &
                - f_grid(loc_Xstg, 0.d0, -0.5d0, 0.5d0) )&
                / (2.d0*delp(ns))
      end if

      if ( nr == nrmax ) then
        dfnspdr = (f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) & 
                - f_grid(loc_Xstg, 0.d0, 0.5d0, -0.5d0) )&
                / (delr)  
      else if ( nr == 1 ) then
        dfnspdr = (f_grid(loc_Xstg, 0.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) )&
                / (delr)
      else
        dfnspdr = (f_grid(loc_Xstg, 0.d0, 0.5d0, 1.5d0) &
                - f_grid(loc_Xstg, 0.d0, 0.5d0, -0.5d0) )&
                / (2.d0 * delr)
      end if

      dfnspdth = (fnsp(nth,np,nr,nsa)-fnsp(nth-1,np,nr,nsa) )&
               / (0.5d0*pl*(delthm(nth,np,nr,nsa) & 
               + delthm(nth-1,np,nr,nsa)) )

      fnspg   = f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0)

      S_x = -1.d0*(Dtp_j(nth,np,nr,nsa)*dfnspdp  &
          + Dttfow(nth,np,nr,nsa)*dfnspdth &
          + Dtrfow(nth,np,nr,nsa)*dfnspdr)  &
          + Fthfow(nth,np,nr,nsa)*fnspg

      nm_Xstg = nma(nth, np, nr)

      if ( S_x <= 0.d0 .and. nl_pnc == -1 ) then !** orbital loss
        return

      else if ( S_x > 0.d0 .and. nl_pnc == -1 ) then

        fact1 = -1.d0/JI(nth,np,nr,nsa)

        dl(nm_Xstg) = dl(nm_Xstg) + fact1*( -1.d0*Dtx(2) &
                    + Fthp*(1.d0-weight(nth,np,nr,nsa)) )

        do nl = 1, nl_Xstg
          if ( ll(nm_Xstg,nl) == nma_boundary(nth,np+1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(1) & 
                           * (1.d0-weight(nth,np+1,nr,nsa))*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth-1,np+1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(1) &
                           * weight(nth,np+1,nr,nsa)*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np-1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(1) & 
                           * (1.d0-weight(nth,np-1,nr,nsa))*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth-1,np-1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(1) & 
                           * weight(nth,np-1,nr,nsa)*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) +( Dtx(2) &
                           + Fthp*weight(nth,np,nr,nsa) )*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np,nr+1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(3) & 
                           * (1.d0-weight(nth,np,nr+1,nsa))*fact1
            
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth-1,np,nr+1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(3) & 
                           * weight(nth,np,nr+1,nsa)*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np,nr-1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(3) & 
                           * (1.d0-weight(nth,np,nr-1,nsa))*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth-1,np,nr-1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(3) & 
                           * weight(nth,np,nr-1,nsa)*fact1
    
          end if
        end do
    
        return

      else

        nm_pnc=nma(loc_pnc(1),loc_pnc(2),loc_pnc(3))

        fact2 = 1.d0/JI(loc_pnc(1),loc_pnc(2),loc_pnc(3),loc_pnc(4))
        nl_pnc = nl_pnc+1 
        ll(nm_pnc,nl_pnc) = nm_Xstg
        al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + fact2 & 
                          * (-1.d0*Dtx(2)+Fthp & 
                          * (1.d0-weight(nth,np,nr,nsa)) )
    
        if ( nma_boundary(nth,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth,np+1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(1) & 
                            * (1.d0-weight(nth,np+1,nr,nsa))*fact2  
        end if
    
        if ( nma_boundary(nth-1,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth-1,np+1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(1) & 
                            * weight(nth,np+1,nr,nsa)*fact2  
        end if
    
        if ( nma_boundary(nth,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth,np-1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(1) & 
                            * (1.d0-weight(nth,np-1,nr,nsa))*fact2  
        end if
    
        if (  nma_boundary(nth-1,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth-1,np-1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(1) & 
                            * weight(nth,np-1,nr,nsa)*fact2  
        end if
    
        if ( nma_boundary(nth,np,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth,np,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) -( Dtx(2) &
                            - Fthp*weight(nth,np,nr,nsa) )*fact2  
        end if
    
        if ( nma_boundary(nth,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth,np,nr+1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(3) & 
                            * (1.d0-weight(nth,np,nr+1,nsa))*fact2  
        end if
    
        if ( nma_boundary(nth-1,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth-1,np,nr+1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(3) & 
                            * weight(nth,np,nr+1,nsa)*fact2  
        end if
    
        if ( nma_boundary(nth,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth,np,nr-1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(3) & 
                            * (1.d0-weight(nth,np,nr-1,nsa))*fact2  
        end if
    
        if ( nma_boundary(nth-1,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(nm_pnc,nl_pnc) = nma(nth-1,np,nr-1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(3) & 
                            * weight(nth,np,nr-1,nsa)*fact2  
        end if  
  
      end if

    else !** pz(ns) < 0.d0
       
      nth = nth_stg(nsa)-1
      loc_Xstg = [nth,np,nr,nsa]
      if ( nl_pnc /= -1 ) then
        loc_pnc = [nth_pnc(nsa)-1,np,nr_rhom_pinch(np,nr,nsa),nsa]
      end if

      nm_Xstg = nma(nth, np, nr)

      pl = pm(np,nsa)
      sl = Sin( thetam(nth,np,nr,nsa) )
      rl = rm(nr)
      deltath= delthm(nth,np,nr,nsa)
  
      Dtx(1)  = Dtp_j(nth+1,np,nr,nsa)/(delp(nsa)*deltath*2.d0)
      Dtx(2)  = Dtt_j(nth+1,np,nr,nsa)/(deltath**2)
      Dtx(3)  = Dtr_j(nth+1,np,nr,nsa)/(deltath*delr*2.d0)
      Fthp    = Fth_j(nth+1,np,nr,nsa)/deltath !**  

      if ( np == npmax ) then
        dfnspdp = (f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) & 
                - f_grid(loc_Xstg, 1.d0, -0.5d0, 0.5d0) )&
                / (delp(ns))
      else if ( np == 1 ) then
        dfnspdp = (f_grid(loc_Xstg, 1.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) )&
                / (delp(ns))
      else
        dfnspdp = (f_grid(loc_Xstg, 1.d0, 1.5d0, 0.5d0) & 
                - f_grid(loc_Xstg, 1.d0, -0.5d0, 0.5d0) )&
                / (2.d0*delp(ns))
      end if

      if ( nr == nrmax ) then
        dfnspdr = (f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) & 
                - f_grid(loc_Xstg, 1.d0, 0.5d0, -0.5d0) )&
                / (delr)
  
      else if ( nr == 1 ) then
        dfnspdr = (f_grid(loc_Xstg, 1.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) )&
                / (delr)
      else
        dfnspdr = (f_grid(loc_Xstg, 1.d0, 0.5d0, 1.5d0) & 
                - f_grid(loc_Xstg, 1.d0, 0.5d0, -0.5d0) )&
                / (2.d0*delr)
      end if

      dfnspdth = (fnsp(nth+1,np,nr,nsa)-fnsp(nth,np,nr,nsa) )&
               / (0.5d0*pl*(delthm(nth+1,np,nr,nsa) & 
               + delthm(nth,np,nr,nsa)) )

      fnspg   = f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0)

      S_x = -1.d0*(Dtpfow(nth+1,np,nr,nsa)*dfnspdp  &
          + Dttfow(nth+1,np,nr,nsa)*dfnspdth &
          + Dtrfow(nth+1,np,nr,nsa)*dfnspdr) &
          + Fthfow(nth+1,np,nr,nsa)*fnspg

      if ( S_x >= 0.d0 .and. nl_pnc == -1 ) then !** orbital loss
        return

      else if ( S_x < 0.d0 .and. nl_pnc == -1 ) then

        fact1 = 1.d0/JI(nth,np,nr,nsa)
        nm_Xstg = nma(nth, np, nr)

        dl(nm_Xstg) = dl(nm_Xstg) + fact1*(Dtx(2) & 
                    + Fthp*weight(nth+1,np,nr,nsa) )
        do nl = 1, nl_Xstg
          if ( ll(nm_Xstg,nl) == nma_boundary(nth+1,np+1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(1) & 
                           * (1.d0-weight(nth+1,np+1,nr,nsa))*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np+1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(1) & 
                           * weight(nth+1,np+1,nr,nsa)*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth+1,np-1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(1) & 
                           * (1.d0-weight(nth+1,np-1,nr,nsa))*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np-1,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(1) & 
                           * weight(nth+1,np-1,nr,nsa)*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth+1,np,nr) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) -( Dtx(2) &
                           - Fthp*(1.d0-weight(nth+1,np,nr,nsa)) )*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth+1,np,nr+1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(3) &
                           * (1.d0-weight(nth+1,np,nr+1,nsa))*fact1
            
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np,nr+1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) - Dtx(3) & 
                           * weight(nth+1,np,nr+1,nsa)*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth+1,np,nr-1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(3) & 
                           * (1.d0-weight(nth+1,np,nr-1,nsa))*fact1
    
          else if ( ll(nm_Xstg,nl) == nma_boundary(nth,np,nr-1) ) then
            al(nm_Xstg,nl) = al(nm_Xstg,nl) + Dtx(3) & 
                           * weight(nth+1,np,nr-1,nsa)*fact1
    
          end if
        end do
      
        return

      else

        nm_pnc = nma(loc_pnc(1),loc_pnc(2),loc_pnc(3))

        fact1 = 1.d0/JI(loc_pnc(1),loc_pnc(2),loc_pnc(3),loc_pnc(4))

        nl_pnc = nl_pnc+1 
        ll(nm_pnc,nl_pnc) = nm_Xstg
        al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + fact1 & 
                          * (Dtx(2)+Fthp*weight(nth+1,np,nr,nsa) )

        if ( nma_boundary(nth+1,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth+1,np+1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(1) & 
                            * (1.d0-weight(nth+1,np+1,nr,nsa))*fact1
  
        else if ( nma_boundary(nth,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth,np+1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(1) & 
                            * weight(nth+1,np+1,nr,nsa)*fact1
  
        else if ( nma_boundary(nth+1,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth+1,np-1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(1) & 
                            * (1.d0-weight(nth+1,np-1,nr,nsa))*fact1
  
        else if ( nma_boundary(nth,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth,np-1,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(1) & 
                            * weight(nth+1,np-1,nr,nsa)*fact1
  
        else if ( nma_boundary(nth+1,np,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth+1,np,nr)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) -( Dtx(2) &
                            - Fthp*(1.d0-weight(nth+1,np,nr,nsa)) )*fact1
  
        else if ( nma_boundary(nth+1,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth+1,np,nr+1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(3) & 
                            * (1.d0-weight(nth+1,np,nr+1,nsa))*fact1
          
        else if ( nma_boundary(nth,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth,np,nr+1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) - Dtx(3) &
                            * weight(nth+1,np,nr+1,nsa)*fact1
  
        else if ( nma_boundary(nth+1,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth+1,np,nr-1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(3) & 
                            * (1.d0-weight(nth+1,np,nr-1,nsa))*fact1
  
        else if ( nma_boundary(nth,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(nm_pnc,nl_pnc) = nma(nth,np,nr-1)
          al(nm_pnc,nl_pnc) = al(nm_pnc,nl_pnc) + Dtx(3) & 
                            * weight(nth+1,np,nr-1,nsa)*fact1
  
        end if
      end if
    end if

  end subroutine IBC_X_stagnation

  subroutine IBC_O_stagnation(np,nr,nsa,nl_Ostg)
  !----------------------------------------------
  !
  !---------------------------------------------

    implicit none
    integer,intent(in) :: np,nr,nsa
    integer,intent(inout) :: nl_Ostg
    real(rkind) :: fact
    double precision :: DIVF, Fthm, Fthp, Dtx(3),pl,sl,deltath,rl
    integer :: ns, nl, nm_Ostg,nth

    ns = ns_nsa(nsa)

    if ( pz(ns) >= 0.d0 ) then
      nth = nth_stg(nsa)-1
      pl = pm(np,nsa)
      sl = Sin( thetam(nth,np,nr,nsa) )
      rl = rm(nr)
      deltath= delthm(nth,np,nr,nsa)
      nm_Ostg = nma(nth,np,nr)
  
      Dtx(1)  = Dtp_j(nth+1,np,nr,nsa)/(delp(nsa)*deltath*2.d0)
      Dtx(2)  = Dtt_j(nth+1,np,nr,nsa)/(deltath**2)
      Dtx(3)  = Dtr_j(nth+1,np,nr,nsa)/(deltath*delr*2.d0)
      Fthp    = Fth_j(nth+1,np,nr,nsa)/deltath !**

      fact = 1.d0/JI(nth,np,nr,nsa)

      dl(nm_Ostg) = dl(nm_Ostg) + (Dtx(2)+Fthp &
                  * weight(nth+1,np,nr,nsa))*fact

      do nl = 1, nl_Ostg
        if ( ll(nm_Ostg,nl) == nma_boundary(nth+1,np+1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(1) & 
                         * (1-weight(nth+1,np+1,nr,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np+1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(1) & 
                         * weight(nth+1,np+1,nr,nsa)*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth+1,np-1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(1) & 
                         * (1.d0-weight(nth+1,np-1,nr,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np+1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(1) & 
                         * weight(nth+1,np-1,nr,nsa)*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth+1,np,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - (Dtx(2) & 
                         + Fthp*(1.d0-weight(nth+1,np,nr,nsa)))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth+1,np,nr+1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(3) & 
                         * (1.d0-weight(nth+1,np,nr+1,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np,nr+1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(3) & 
                         * weight(nth+1,np,nr+1,nsa)*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth+1,np,nr-1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(3) &
                         * (1.d0-weight(nth+1,np,nr-1,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np,nr-1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(3) &
                         * weight(nth+1,np,nr-1,nsa)*fact
        end if
      end do

    else !** pz(ns) < 0.d0
      nth = nth_stg(nsa)
      pl = pm(np,nsa)
      sl = Sin( thetam(nth,np,nr,nsa) )
      rl = rm(nr)
      deltath= delthm(nth,np,nr,nsa)
      nm_Ostg = nma(nth,np,nr)
  
      Dtx(1)  = Dtp_j(nth,np,nr,nsa)/(delp(nsa)*deltath*2.d0)
      Dtx(2)  = Dtt_j(nth,np,nr,nsa)/(deltath**2)
      Dtx(3)  = Dtr_j(nth,np,nr,nsa)/(deltath*delr*2.d0)
      Fthp    = Fth_j(nth,np,nr,nsa)/deltath !**

      fact = -1.d0/JI(nth,np,nr,nsa)

      dl(nm_Ostg) = dl(nm_Ostg) + fact*( -1.d0*Dtx(2) &
                  + Fthp*(1.d0-weight(nth,np,nr,nsa)) )
      do nl = 1, nl_Ostg
        if ( ll(nm_Ostg,nl) == nma_boundary(nth,np+1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(1) & 
                         * (1.d0-weight(nth,np+1,nr,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth-1,np+1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(1) & 
                         * weight(nth,np+1,nr,nsa)*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np-1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(1) & 
                         * (1.d0-weight(nth,np-1,nr,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth-1,np-1,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(1) & 
                         * weight(nth,np-1,nr,nsa)*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np,nr) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) +(Dtx(2) & 
                         + Fthp*weight(nth,np,nr,nsa) )*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np,nr+1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(3) & 
                         * (1.d0-weight(nth,np,nr+1,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth-1,np,nr+1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) - Dtx(3) & 
                         * weight(nth,np,nr+1,nsa)*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth,np,nr-1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(3) &
                         * (1.d0-weight(nth,np,nr-1,nsa))*fact
        else if ( ll(nm_Ostg,nl) == nma_boundary(nth-1,np,nr-1) ) then
          al(nm_Ostg,nl) = al(nm_Ostg,nl) + Dtx(3) & 
                         * weight(nth,np,nr-1,nsa)*fact
        end if
      end do
    end if

  end subroutine IBC_O_stagnation

  function f_grid(loc,dth,dp,dr)
  !----------------------------------------------------
  ! dx = -1/2, 1/2, 3/2    -> half integer grid point
  ! dx = 0, 1              -> integer grid point
  !----------------------------------------------------

    implicit none
    integer,intent(in) :: loc(4)
    real(rkind) :: f_grid
    real(rkind),intent(in) :: dth, dp, dr
    real(rkind) :: f_plus, f_minus
    integer :: nth, np, nr, nsa
    integer :: ith, ip, ir

    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    if ( int(2.d0*dth) == 2*INT(dth) ) then
      ith = nth+ int(dth)
      ip  = np + int(dp-0.5d0)
      ir  = nr + int(dr-0.5d0)

      f_plus  = (1.d0-weight(ith, ip, ir, nsa)) * fnsp(ith,ip,ir,nsa)
      f_minus = weight(ith, ip, ir, nsa) * fnsp(ith-1,ip,ir,nsa)

    else if ( int(2.d0*dp) == 2*INT(dp) ) then
      ith = nth + int(dth-0.5d0)
      ip  = np  + int(dp)
      ir  = nr  + int(dr -0.5d0)

      f_plus  = (1.d0-weighp(ith, ip, ir, nsa)) * fnsp(ith,ip,ir,nsa)
      f_minus = weighp(ith, ip, ir, nsa) * fnsp(ith,ip-1,ir,nsa)

    else if ( int(2.d0*dr) == 2*INT(dr) ) then
      ith = nth + int(dth-0.5d0)
      ip  = np  + int(dp -0.5d0)
      ir  = nr  + int(dr)

      f_plus  = (1.d0-weighr(ith, ip, ir, nsa)) * fnsp(ith,ip,ir,nsa)
      f_minus = weighr(ith, ip, ir, nsa) * fnsp(ith,ip,ir-1,nsa)

    end if

    f_grid = f_plus + f_minus

  end function f_grid

  function nma_boundary(nth,np,nr) result(n)
  !-------------------------------------------------
  !
  !-------------------------------------------------

    implicit none
    logical :: nth_,np_,nr_
    integer :: n
    integer,intent(in) :: nth,np,nr

    nth_ = 1 <= nth .and. nth <= nthmax 
    np_  = 1 <= np  .and. np  <= npmax 
    nr_  = 1 <= nr  .and. nr  <= nrmax 

    if ( nth_ .and. np_ .and. nr_ ) then
      n = nma(nth,np,nr)
    else
      n = -1
    end if

  end function nma_boundary


end module fowexec