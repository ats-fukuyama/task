!     $Id: fpexec.f90,v 1.29 2013/02/08 07:36:24 nuga Exp $

!     **************************
!        EXECUTE TIME ADVANCE
!     **************************

MODULE fowexec

  use fpcomm
  use fowcomm

contains

  subroutine fow_exec(NSA,IERR,its)

    USE libmpi
    USE libmtx
    USE fpmpi
    IMPLICIT NONE
    integer:: NSA, NP, NTH, NR, NL, NM, NN, NS
    integer:: NTHS, NLL
    integer:: IERR,its,i,j,ll1
    integer:: imtxstart1,imtxend1
    real(8),dimension(nmend-nmstart+1):: BM_L
    real(8),dimension(nthmax):: sendbuf_p, recvbuf_p
    real(8),dimension(nthmax*(npend-npstart+1)):: sendbuf_r, recvbuf_r

    double precision :: begin_time, end_time

    call cpu_time(begin_time)

    NS=NS_NSA(NSA)

    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          f(nth,np,nr) = fnsp(nth,np,nr,nsa)
        end do
      end do
    end do

    CALL mtx_set_communicator(comm_nrnp) !3D

    !     ----- Set up matrix solver -----
    CALL mtx_setup(imtxsize,imtxstart1,imtxend1,imtxwidth)
    IF(imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend) THEN
        WRITE(6,*) 'XX fp_exec: '
        WRITE(6,*) '   imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend'
        WRITE(6,*) '   imtxstart1,imtxstart = ',imtxstart1,imtxstart
        WRITE(6,*) '   imtxend1,imtxend     = ',imtxend1,imtxend
        STOP
    ENDIF

    !     ----- Set up weight array -----

    CALL fowweight(NSA,IERR)

    !     ----- Set up index array NMA -----
    !               NM: line number of the coefficient matrix
    !               NL: 

    CALL SET_FM_NMA(NSA,FNSM)

    DO NM=NMSTART,NMEND
        NLMAX(NM)=0
        BM(NM)=0.D0
        DO NL=1,NLMAXM
          LL(NM,NL)=0
          AL(NM,NL)=0.D0
        ENDDO
    ENDDO


    !     ----- Calculate matrix coefficients in a row -----
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          nm = nma(nth,np,nr)
          call fowsetm(nth, np, nr, nsa, nlmax(nm))
        end do
      end do
    end do

    !     ----- Diagonal term -----
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          nm = nma(nth,np,nr)
          bm(nm) = (1.d0+(1.d0-rimpl)*delt*dl(nm))*fm(nm) &
                    +delt*spp(nth,np,nr,nsa)
          if(nm.ge.imtxstart.and.nm.le.imtxend) then
            call mtx_set_matrix(nm, nm, 1.d0-rimpl*delt*dl(nm))
            call mtx_set_vector(nm, fm(nm))
          ENDIF
        end do
      end do
    end do

    !     ----- Off diagonal term -----

    DO NM=NMSTART,NMEND ! LHS
      IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
        DO NL=1,NLMAX(NM)
            IF(LL(NM,NL).NE.0) THEN
              CALL mtx_set_matrix(nm,LL(NM,NL),-RIMPL*DELT*AL(NM,NL))
            ENDIF
        ENDDO
      ENDIF
    ENDDO
close(21)
close(901)
    !     ----- Source vector: contribution from off-diagonal term -----

    DO NM=NMSTART,NMEND ! RHS
      DO NL=1,NLMAX(NM)
        NN=LL(NM,NL)
        IF(NN.NE.0) THEN
          IF(NN.ge.NMSTART-NTHMAX.and.NN.le.NMEND+NTHMAX)THEN
            BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM(NN)
          ELSEIF(NN.lt.NMSTART-NTHMAX)THEN
            BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM_shadow_m(NN)
          ELSE
            BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM_shadow_p(NN)
          END IF
        ENDIF
      ENDDO
      IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
         CALL mtx_set_source(nm,BM(NM))
      ENDIF
    ENDDO


    !     ----- Solve matrix equation -----

    CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) ! ncom is nessesary for MUMPS not PETSc
    ierr=0

    !     ----- Get solution vector -----

    CALL mtx_get_vector(BM_L)

    DO NR=NRSTART, NREND
        DO NP=NPSTART, NPEND
          DO NTH=1,NTHMAX
              NM=NMA(NTH,NP,NR)
              IF(ABS(BM_L(NM-NMSTART+1)).LT.1.D-100) THEN
                FNS0(NTH,NP,NR,NSA)=0.D0
              ELSE
                FNS0(NTH,NP,NR,NSA)=BM_L(NM-NMSTART+1)
              END IF
          END DO
        END DO
    END DO
    !     shadow requires to communicate
    CALL mtx_set_communicator(comm_np)
    DO NR=NRSTART, NREND
        CALL shadow_comm_np(NR,NSA)
    END DO
    CALL mtx_set_communicator(comm_nr)
    CALL shadow_comm_nr(NSA)
    CALL mtx_set_communicator(comm_nrnp) !3D


    !     ----- Clean up matrix solver -----
    CALL mtx_cleanup

    CALL mtx_reset_communicator

    call cpu_time(end_time)
    write(6,'(A,I0,A,ES10.3,A)')'fowexec time (nsa=',nsa,') : ',end_time-begin_time,'[sec]'

    RETURN
  END subroutine fow_exec

  !     ---------------------------------

  SUBROUTINE SET_FM_NMA(NSA,func_in)

    IMPLICIT NONE
    integer:: NTH, NP, NR, NSA, NM, NRS, NPS
    double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND), &
          intent(IN):: func_in

    IF(NRSTART.eq.1)THEN
        NRS=1
    ELSE
        NRS=NRSTART-1
    END IF

    DO NR=NRSTARTW,NRENDWM
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              NM=NTH+NTHMAX*(NP-1)+NPMAX*NTHMAX*(NR-1)
              NMA(NTH,NP,NR)=NM
          END DO
        END DO
    END DO

    DO NR=NRSTART,NREND
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              NM=NMA(NTH,NP,NR)
              FM(NM)=func_in(NTH,NP,NR,NSA)
          ENDDO
        ENDDO
    ENDDO
    NR=NRSTARTW
    IF(NR.ne.NRSTART)THEN
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              NM=NMA(NTH,NP,NR)
              FM_shadow_m(NM)=func_in(NTH,NP,NR,NSA)
          ENDDO
        ENDDO
    END IF
    NR=NRENDWM
    IF(NR.ne.NREND)THEN
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              NM=NMA(NTH,NP,NR)
              FM_shadow_p(NM)=func_in(NTH,NP,NR,NSA)
          ENDDO
        ENDDO
    END IF


  END SUBROUTINE SET_FM_NMA

  !
  !     ***************************
  !        CALCULATION OF WEIGHT
  !     ***************************
  !
  subroutine fowweight(NSA,IERR) ! proposed by Chang and Cooper [30] in Karney

    use fpcomm
    use fowcomm

    implicit none
    integer:: nsa, np, nth, nr, nthl, nthr, ns
    integer:: ierr
    real(8):: dfdth, fvel, dfdp
    
    integer :: nrl, nrr, npl, npr
    double precision :: dfdr, delps_l, width_p, width_r, width_t

    !     +++++ calculation of weigthing (including off-diagonal terms) +++++

    real(8)::epswt=1.d-70

    ns = ns_nsa(nsa)

    do nr = nrstart, nrend
      delps_l = psimg(nr+1)-psimg(nr)
      do np = npstart, npendwg
        do nth = 1, nthmax
            dfdth = 0.d0
            dfdr  = 0.d0
            if ( np /= 1 ) then
              nthl = min(nth+1,nthmax)
              nthr = max(nth-1,1)
              nrl  = min(nr+1,nrmax)
              nrr  = max(nr-1,1)

              if ( np == npmax+1 ) then
                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                          /(2.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)*f(nth,np-1,nr))

                  dfdr  = (f(nth,np-1,nrl)-f(nth,np-1,nrr)) &
                          /(2.d0*delps_l*f(nth,np-1,nr))
                end if
              else

                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  if ( abs(f(nth,np,nr)) > epswt ) then
                    dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                            /(4.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)*f(nth,np-1,nr)) &
                          + (f(nthl,np  ,nr)-f(nthr,np  ,nr)) &
                            /(4.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)*f(nth,np  ,nr))

                    dfdr  = (f(nth,np-1,nrl)-f(nthr,np-1,nrr)) &
                            /(4.d0*delps_l*f(nth,np-1,nr)) &
                          + (f(nthl,np  ,nrl)-f(nthr,np  ,nrr)) &
                            /(4.d0*delps_l*f(nth,np  ,nr))
                  else
                    dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                          /(2.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)*f(nth,np-1,nr)) 

                    dfdr  = (f(nth,np-1,nrl)-f(nth,np-1,nrr)) &
                          /(2.d0*delps_l*f(nth,np-1,nr)) 
                  end if
                else
                  if ( abs(f(nth,np,nr)) > epswt ) then
                    dfdth = (f(nthl,np  ,nr)-f(nthr,np  ,nr)) &
                          /(2.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)*f(nth,np  ,nr))

                    dfdr  = (f(nth,np  ,nrl)-f(nth,np  ,nrr)) &
                          /(2.d0*delps_l*f(nth,np  ,nr))
                  end if
                end if
              end if
            end if
            fvel = Fppfow(nth,np,nr,nsa)-Dptfow(nth,np,nr,nsa)*dfdth-Dprfow(nth,np,nr,nsa)*dfdr
            weighp(nth,np,nr,nsa) = fowwegh(-delp(ns)*fvel,dppfow(nth,np,nr,nsa))

        end do
      end do
    end do

    do nr = nrstart, nrend
      delps_l = psimg(nr+1)-psimg(nr)
      do np = npstart, npend
        do nth = 1, nthmax+1
          dfdp = 0.d0          
          npl = min(npmax,np+1)
          npr = max(1, np-1)
          width_p = dble(npl-npr)

          dfdr = 0.d0
          nrl = min(nrmax,nr+1)
          nrr = max(1,nr-1)
          width_r = dble(nrl-nrr)
          
          if ( nth == 1 ) then
            if ( abs(f(nth,np,nr)) > epswt ) then
              dfdp = (f(nth,npl,nr)-f(nth,npr,nr))/(width_p*delp(ns)*f(nth,np,nr))
              dfdr = (f(nth,np,nrl)-f(nth,np,nrr))/(width_r*delps_l   *f(nth,np,nr))  
            end if
          else if ( nth == nthmax+1 ) then
            if ( abs(f(nthmax,np,nr)) > epswt ) then
              dfdp = (f(nthmax,npl,nr)-f(nthmax,npr,nr))/(width_p*delp(ns)*f(nthmax,np,nr))
              dfdr = (f(nthmax,np,nrl)-f(nthmax,np,nrr))/(width_r*delps_l*f(nthmax,np,nr))  
            end if
          else
            if ( abs(f(nth,np,nr)) > epswt .and. abs(f(nth-1,np,nr)) > epswt ) then
              dfdp = ( &
                (f(nth-1,npl,nr)-f(nth-1,npr,nr))/(width_p*delp(ns)*f(nth-1,np,nr))+&
                (f(nth  ,npl,nr)-f(nth  ,npr,nr))/(width_p*delp(ns)*f(nth  ,np,nr)) &
              )/2.d0

              dfdr = ( &
                (f(nth-1,np,nrl)-f(nth-1,np,nrr))/(width_r*delps_l*f(nth-1,np,nr))+&
                (f(nth  ,np,nrl)-f(nth  ,np,nrr))/(width_r*delps_l*f(nth  ,np,nr)) &
              )
            else if ( abs(f(nth,np,nr)) > epswt .and. abs(f(nth-1,np,nr)) <= epswt ) then
              dfdp = (f(nth  ,npl,nr)-f(nth  ,npr,nr))/(width_p*delp(ns)*f(nth  ,np,nr))
              dfdr = (f(nth  ,np,nrl)-f(nth  ,np,nrr))/(width_r*delps_l*f(nth  ,np,nr))
            else if ( abs(f(nth,np,nr)) <= epswt .and. abs(f(nth-1,np,nr)) > epswt ) then
              dfdp = (f(nth-1,npl,nr)-f(nth-1,npr,nr))/(width_p*delp(ns)*f(nth-1,np,nr))
              dfdr = (f(nth-1,np,nrl)-f(nth-1,np,nrr))/(width_r*delps_l*f(nth-1,np,nr))
            end if
          end if
          fvel = Fthfow(nth,np,nr,nsa)-Dtpfow(nth,np,nr,nsa)*dfdp-Dtrfow(nth,np,nr,nsa)*dfdr
          if ( nth <= nthmax ) then
            weight(nth,np,nr,nsa) = fowwegh(-delthm(nth,np,nr,nsa)*fvel,Dttfow(nth,np,nr,nsa))
          else
            weight(nth,np,nr,nsa) = fowwegh(-delthm(nth-1,np,nr,nsa)*fvel,Dttfow(nth,np,nr,nsa))
          end if
        end do
      end do
    end do

    do nr = nrstart, nrendwg
      if ( nr == nrmax+1 ) then
        delps_l = psimg(nr)-psimg(nr-1)
      else
        delps_l = psimg(nr+1)-psimg(nr)
      end if
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
            if ( abs(f(nth,np,nr)) > epswt ) then
              dfdp = (f(nth,npl,1)-f(nth,npr,1))/(width_p*delp(ns)*f(nth,np,1))
              dfdth= (f(nthl,np,1)-f(nthr,np,1))/(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,1))  
            end if
          else if ( nr == nrmax+1 ) then
            if ( abs(f(nth,np,nrmax)) > epswt ) then
              dfdp = (f(nth,npl,nrmax)-f(nth,npr,nrmax))/(width_p*delp(ns)*f(nth,np,nrmax))
              dfdth= (f(nthl,np,nrmax)-f(nthr,np,nrmax))/(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nrmax))
            end if
          else
            if ( abs(f(nth,np,nr)) > epswt .and. abs(f(nth,np,nr-1)) > epswt ) then
              dfdp = ( &
                (f(nth,npl,nr-1)-f(nth,npr,nr-1))/(width_p*delp(ns)*f(nth,np,nr-1))+&
                (f(nth,npl,nr  )-f(nth,npr,nr  ))/(width_p*delp(ns)*f(nth,np,nr  )) &
              )/2.d0

              if ( nr == 2 ) then
                dfdth= (f(nthl,np,nr)-f(nthr,np,nr))/(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nr))
              else
                dfdth= ( &
                  (f(nthl,np,nr-1)-f(nthr,np,nr-1))/(width_t*delthm_rg(nth,np,nr-1,nsa)*f(nth,np,nr))+&
                  (f(nthl,np,nr  )-f(nthr,np,nr  ))/(width_t*delthm_rg(nth,np,nr  ,nsa)*f(nth,np,nr)) &
                )/2.d0
              end if

            else if ( abs(f(nth,np,nr)) > epswt .and. abs(f(nth,np,nr-1)) <= epswt ) then
              dfdp = (f(nth,npl,nr  )-f(nth,npr,nr  ))/(width_p*delp(ns)*f(nth,np,nr))

              dfdth= (f(nthl,np,nr  )-f(nthr,np,nr  ))/(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nr))

            else if ( abs(f(nth,np,nr)) <= epswt .and. abs(f(nth,np,nr-1)) > epswt ) then
              dfdp = (f(nth,npl,nr-1)-f(nth,npr,nr-1))/(width_p*delp(ns)*f(nth,np,nr-1))

              dfdth= (f(nthl,np,nr-1)-f(nthr,np,nr-1))/(width_t*delthm_rg(nth,np,nr-1,nsa)*f(nth,np,nr-1))
            end if
          end if
          fvel = Frrfow(nth,np,nr,nsa)-Drpfow(nth,np,nr,nsa)*dfdp-Drtfow(nth,np,nr,nsa)*dfdth
          weighr(nth,np,nr,nsa) = fowwegh(-delps_l*fvel,Drrfow(nth,np,nr,nsa))
        end do
      end do
    end do

    return
  end subroutine fowweight

  ! ************************************************
  !     WEIGHTING FUNCTION FOR CONVECTION EFFECT
  ! ************************************************

  FUNCTION fowwegh(X,Y)

    IMPLICIT NONE
    real(8):: X, Y, Z
    real(8):: fowwegh

    IF(ABS(Y).LT.1.D-70) THEN
        IF(X.GT.0.D0) THEN
          fowwegh=0.D0
        ELSEIF(X.LT.0.D0) THEN
          fowwegh=1.D0
        ELSE
          fowwegh=0.5D0
        end if
    ELSE
        Z=X/Y
        IF(ABS(Z).LT.1.D-5)THEN
          fowwegh=0.5D0-Z/12.D0+Z**3/720.D0
        ELSE IF(Z.GE.100.D0)THEN
          fowwegh=1.D0/Z
        ELSE IF(Z.LE.-100.D0)THEN
          fowwegh=1.D0/Z+1.D0
        ELSE
          fowwegh=1.D0/Z-1.D0/(EXP(Z)-1.D0)
        END IF
    ENDIF
    RETURN
  END FUNCTION fowwegh

  ! ******************************************
  !     Calculation of matrix coefficients
  ! ******************************************

  subroutine fowsetm(NTH,NP,NR,NSA,NL)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NTH,NP,NR,NSA
    INTEGER,INTENT(OUT):: NL
    INTEGER:: NM

    integer:: IERR, NS
    real(8):: DFDTH, FVEL, DVEL, DFDP, DFDB
    real(8):: DPPM, DPPP, DPTM, DPTP, DTPM, DTPP, DTTM, DTTP
    ! real(8):: WPM, WPP, WTM, WTP, VPM, VPP, VTM, VTP
    real(8):: WTB, VTB, WRBM, VRBM, WRBP, VRBP
    real(8):: DIVDPP, DIVDPT, DIVDTP, DIVDTT, DIVFPP, DIVFTH
    real(8):: DRRM,DRRP,WRM,WRP,VRM,VRP,DIVDRR,DIVFRR
    DOUBLE PRECISION:: WPBM, VPBM, WPBP, VPBP
    DOUBLE PRECISION:: PL, SL, RL
    ! fow extension
    double precision :: DIVDPR, DIVDTR, DIVDRP, DIVDRT
    double precision :: deltap, deltath, deltaps
    double precision :: DIVD(3,3), DIVF(3), del(3), Ffow(2,3)
    double precision :: D_term, F_term
    integer :: si, sj, sk, sign_to_index(-1:1), alpha, beta, gama, loc(4), loc_pnc(4)
    integer :: boundary_flag
              ! out of external boundary -> boundary_flag = 1
    type(Xstg_as_pnc_point) :: xap
    integer :: i, nrpo, nthpo


    NS=NS_NSA(NSA)

    NL=0
    NM=NMA(NTH,NP,NR)

    Ffow(1,1) = Fppfow(nth,np,nr,nsa)*pg(np,nsa)**2
    Ffow(2,1) = Fppfow(nth,np+1,nr,nsa)*pg(np+1,nsa)**2

    Ffow(1,2) = Fthfow(nth,np,nr,nsa)*sin(thetamg(nth,np,nr,nsa))
    Ffow(2,2) = Fthfow(nth+1,np,nr,nsa)*sin(thetamg(nth+1,np,nr,nsa))

    Ffow(1,3) = Frrfow(nth,np,nr,nsa)*psimg(nr)
    Ffow(2,3) = Frrfow(nth,np,nr+1,nsa)*psimg(nr+1)


    ! discretized (div(d/dX))_Y
    PL = pm(np,nsa)
    SL = sin(thetam(nth,np,nr,nsa))
    RL = psim(nr)
    deltap = delp(nsa)
    deltath= thetamg(nth+1,np,nr,nsa)-thetamg(nth,np,nr,nsa)
    deltaps=psimg(nr+1)-psimg(nr)

    DIVD(1,1) = 1.d0/(pl**2 * deltap**2)
    DIVD(1,2) = 1.d0/(pl**2 * deltap * deltath * 2.d0)
    DIVD(1,3) = 1.d0/(pl**2 * deltap * deltaps * 2.d0)
    DIVD(2,1) = 1.d0/(pl * sl * deltap * deltath * 2.d0)
    DIVD(2,2) = 1.d0/(pl * sl * deltath**2 )
    DIVD(2,3) = 1.d0/(pl * sl * deltath * deltaps * 2.d0)
    DIVD(3,1) = 1.d0/(rl * deltap  * deltaps * 2.d0)
    DIVD(3,2) = 1.d0/(rl * deltath * deltaps * 2.d0)
    DIVD(3,3) = 1.d0/(rl * deltaps**2)
    DIVF(1)   = 1.d0/(pl**2 * deltap)
    DIVF(2)   = 1.d0/(pl * sl * deltath)
    DIVF(3)   = 1.d0/(rl * deltaps)

    sign_to_index(-1) = 1
    sign_to_index(1)  = 2
    sign_to_index(0)  = 0
    loc = [nth,np,nr,nsa]

    ! term of f(nth, np, nr)
    DL(NM) = 0.d0
    do alpha = 1, 3
      do si = -1, 1, 2
        D_term = -1.d0*Dfow(alpha,alpha,sign_to_index(si)-1,0,loc)*DIVD(alpha,alpha)
        F_term = -1.d0*si*Ffow(sign_to_index(si),alpha)*w(si,alpha,0,sign_to_index(si)-1,0,loc)*DIVF(alpha)
        dl(nm) = dl(nm) + (D_term + F_term)/JI(nth,np,nr,nsa)
      end do
    end do
    if ( nth == nth_pnc(nsa) .and. theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
      D_term = Dfow(2,2,0,0,loc)*DIVD(2,2)*IBCflux_ratio(np,nr,nsa)
      F_term = -1.d0*Ffow(1,2)*w(-1,2,0,0,0,loc)*DIVF(2)*IBCflux_ratio(np,nr,nsa)
      dl(nm) = dl(nm) + (D_term + F_term)/JI(nth,np,nr,nsa)

    end if

    ! terms of f(i+si, j, k) 
    do alpha = 1, 3
      do si = -1, 1, 2
        boundary_flag = check_external_boundary(alpha,0,si,0,loc)
        if ( boundary_flag == 1  ) cycle

        D_term = Dfow(alpha,alpha,sign_to_index(si)-1,0,loc)*DIVD(alpha,alpha)
        F_term = -1.d0*si*Ffow(sign_to_index(si),alpha)*w(-1*si,alpha,0,sign_to_index(si)-1,0,loc)*DIVF(alpha)

        do beta = 1, 3
          if ( alpha == beta ) cycle
          do sj = -1, 1, 2
            D_term = D_term + si * sj*Dfow(beta,alpha,sign_to_index(sj)-1,0,loc) &
                    *w(sj,beta,alpha,sign_to_index(sj)-1,si,loc)*DIVD(beta,alpha)
          end do
        end do

        nl = nl+1
        ll(nm,nl) = get_nma(alpha,0,si,0,loc)
        al(nm,nl) = al(nm,nl) + (D_term + F_term)/JI(nth,np,nr,nsa)

        if ( nth == nth_pnc(nsa) .and. theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
          if ( alpha == 1 ) then
            if ( si == 1 ) then
              D_term = Dfow(2,1,0,0,loc)*w(-1,2,1,0,1,loc)*DIVD(2,1)
              F_term = 0.d0
            else
              D_term = -1.d0*Dfow(2,1,0,0,loc)*w(-1,2,1,0,-1,loc)*DIVD(2,1)
              F_term = 0.d0
            end if
          else if ( alpha == 2 ) then
            if ( si == -1 ) then
              D_term = -1.d0*Dfow(2,2,0,0,loc)*DIVD(2,2)
              F_term = -1.d0*Ffow(1,2)*DIVF(2)*w(1,2,0,0,0,loc)
            end if
          else
            if ( si == 1 ) then
              D_term = Dfow(2,3,0,0,loc)*w(-1,2,3,0,1,loc)*DIVD(2,3)
              F_term = 0.d0
            else
              D_term = Dfow(2,3,0,0,loc)*w(-1,2,3,0,-1,loc)*DIVD(2,3)
              F_term = 0.d0
            end if           
          end if
          al(nm,nl) = al(nm,nl) + IBCflux_ratio(np,nr,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
        end if

        if ( ABS(al(nm,nl)) < 1.d-70 ) then
          ll(nm,nl) = 0
          al(nm,nl) = 0.d0
          nl = nl-1
        end if

      end do
    end do

    ! terms of f(i+si, j+sj, k) 
    do alpha = 1, 3
      do beta = 1, 3
        if ( alpha >= beta) cycle
        do si = -1, 1, 2
          do sj = -1, 1, 2

            boundary_flag = check_external_boundary(alpha,beta,si,sj,loc)
            if ( boundary_flag == 1 ) cycle
            
            D_term = si*sj*( &
              Dfow(alpha,beta,sign_to_index(si)-1,0,loc)*DIVD(alpha,beta) &
              *w(-1*si,alpha,beta,sign_to_index(si)-1,sj,loc) &
              +Dfow(beta,alpha,sign_to_index(sj)-1,0,loc)*DIVD(beta,alpha) &
              *w(-1*sj,beta,alpha,sign_to_index(sj)-1,si,loc) &
            )

            nl = nl+1
            ll(nm,nl) = get_nma(alpha,beta,si,sj,loc)
            al(nm,nl) = al(nm,nl) + D_term/JI(nth,np,nr,nsa)

            if ( nth == nth_pnc(nsa) .and. theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              if ( alpha == 1 .and. beta == 2) then
                if ( si == 1 .and. sj == -1 ) then
                  D_term = Dfow(2,1,0,0,loc)*w(1,2,1,0,1,loc)*DIVD(2,1)
                else if ( si == -1 .and. sj == -1 ) then
                  D_term = -1.d0*Dfow(2,1,0,0,loc)*w(1,2,1,0,-1,loc)*DIVD(2,1)
                end if

              else if ( alpha == 2 .and. beta == 3 ) then
                if ( si == -1 .and. sj == 1 ) then
                  D_term = Dfow(2,3,0,0,loc)*w(1,2,3,0,1,loc)*DIVD(2,3)
                else if ( si == -1 .and. sj == -1 ) then
                  D_term = -1.d0*Dfow(2,3,0,0,loc)*w(1,2,3,0,-1,loc)+DIVD(2,3)
                end if
              end if
              al(nm,nl) = al(nm,nl) + IBCflux_ratio(np,nr,nsa)*D_term/JI(nth,np,nr,nsa)
            end if

            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    

          end do
        end do
      end do
    end do

    if ( aefp(nsa) > 0.d0 .and. nth == nth_stg(nsa) .or. aefp(nsa) < 0.d0 .and. nth == nth_stg(nsa)-1 ) then
      xap = Xstg_as_pncp(np,nr,nsa)

      if ( xap%number >= 1 ) then
        do i = 1, xap%number
          nthpo = nth_pnc(nsa)
          nrpo  = xap%nr(i)
          loc_pnc = [nthpo, np, nrpo, nsa]

          Ffow(1,1) = Fppfow(nthpo,np,nrpo,nsa)*pg(np,nsa)**2
          Ffow(2,1) = Fppfow(nthpo,np+1,nrpo,nsa)*pg(np+1,nsa)**2      
          Ffow(1,2) = Fthfow(nthpo,np,nrpo,nsa)*sin(thetamg(nthpo,np,nrpo,nsa))
          Ffow(2,2) = Fthfow(nthpo+1,np,nrpo,nsa)*sin(thetamg(nthpo+1,np,nrpo,nsa))
          Ffow(1,3) = Frrfow(nthpo,np,nrpo,nsa)*psimg(nrpo)
          Ffow(2,3) = Frrfow(nthpo,np,nrpo+1,nsa)*psimg(nrpo+1)
          PL = pm(np,nsa)
          SL = sin(thetam(nthpo,np,nrpo,nsa))
          RL = psim(nrpo)
          deltap = delp(nsa)
          deltath= thetamg(nthpo+1,np,nrpo,nsa)-thetamg(nthpo,np,nrpo,nsa)
          deltaps=psimg(nrpo+1)-psimg(nrpo)      
          DIVD(1,1) = 1.d0/(pl**2 * deltap**2)
          DIVD(1,2) = 1.d0/(pl**2 * deltap * deltath * 2.d0)
          DIVD(1,3) = 1.d0/(pl**2 * deltap * deltaps * 2.d0)
          DIVD(2,1) = 1.d0/(pl * sl * deltap * deltath * 2.d0)
          DIVD(2,2) = 1.d0/(pl * sl * deltath**2 )
          DIVD(2,3) = 1.d0/(pl * sl * deltath * deltaps * 2.d0)
          DIVD(3,1) = 1.d0/(rl * deltap  * deltaps * 2.d0)
          DIVD(3,2) = 1.d0/(rl * deltath * deltaps * 2.d0)
          DIVD(3,3) = 1.d0/(rl * deltaps**2)
          DIVF(1)   = 1.d0/(pl**2 * deltap)
          DIVF(2)   = 1.d0/(pl * sl * deltath)
          DIVF(3)   = 1.d0/(rl * deltaps)      

          D_term = Dfow(2,2,0,0,loc_pnc)*DIVD(2,2)
          F_term = -1.d0*Ffow(1,2)*w(-1,2,0,0,0,loc_pnc)*DIVF(2)
          dl(nm) = dl(nm) - IBCflux_ratio(np,nrpo,nsa) * ( D_term + F_term )/JI(nth,np,nr,nsa)

          boundary_flag = check_external_boundary(1,0,1,0,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = Dfow(2,1,0,0,loc_pnc)*w(-1,2,1,0,1,loc_pnc)*DIVD(2,1)
            F_term = 0.d0
            nl = nl+1
            ll(nm,nl) = nma(nthpo,np+1,nrpo)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(1,0,-1,0,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = -1.d0*Dfow(2,1,0,0,loc_pnc)*w(-1,2,1,0,-1,loc_pnc)*DIVD(2,1)
            F_term = 0.d0
            nl = nl+1
            ll(nm,nl) = nma(nthpo,np-1,nrpo)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(2,0,-1,0,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = -1.d0*Dfow(2,2,0,0,loc_pnc)*DIVD(2,2)
            F_term = -1.d0*Ffow(1,2)*DIVF(2)*w(1,2,0,0,0,loc_pnc)
            nl = nl+1
            ll(nm,nl) = nma(nthpo-1,np,nrpo)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(3,0,1,0,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = Dfow(2,3,0,0,loc_pnc)*w(-1,2,3,0,1,loc_pnc)*DIVD(2,3)
            F_term = 0.d0
            nl = nl+1
            ll(nm,nl) = nma(nthpo,np,nrpo+1)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(3,0,-1,0,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = Dfow(2,3,0,0,loc_pnc)*w(-1,2,3,0,-1,loc_pnc)*DIVD(2,3)
            F_term = 0.d0
            nl = nl+1
            ll(nm,nl) = nma(nthpo,np,nrpo-1)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(1,2,1,-1,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = Dfow(2,1,0,0,loc_pnc)*w(1,2,1,0,1,loc_pnc)*DIVD(2,1)
            nl = nl+1
            ll(nm,nl) = nma(nthpo-1,np+1,nrpo)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*( D_term + F_term )/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(1,2,-1,-1,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = -1.d0*Dfow(2,1,0,0,loc_pnc)*w(1,2,1,0,-1,loc_pnc)*DIVD(2,1)
            nl = nl+1
            ll(nm,nl) = nma(nthpo-1,np-1,nrpo)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*D_term/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(2,3,-1,1,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = Dfow(2,3,0,0,loc_pnc)*w(1,2,3,0,1,loc_pnc)*DIVD(2,3)
            nl = nl+1
            ll(nm,nl) = nma(nthpo-1,np,nrpo+1)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*D_term/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

          boundary_flag = check_external_boundary(2,3,-1,-1,loc_pnc)
          if ( boundary_flag == 0 ) then
            D_term = -1.d0*Dfow(2,3,0,0,loc_pnc)*w(1,2,3,0,-1,loc_pnc)+DIVD(2,3)
            nl = nl+1
            ll(nm,nl) = nma(nthpo-1,np,nrpo-1)
            al(nm,nl) = al(nm,nl) - IBCflux_ratio(np,nrpo,nsa)*D_term/JI(nth,np,nr,nsa)
            if ( ABS(al(nm,nl)) < 1.d-70 ) then
              ll(nm,nl) = 0
              al(nm,nl) = 0.d0
              nl = nl-1
            end if    
          end if

        end do  
      end if
      
    end if



    SPP(NTH,NP,NR,NSA) &
            =( SPPB(NTH,NP,NR,NSA) )

    RETURN
  END subroutine fowsetm

  function Dfow(alpha,beta,si,sj,loc)
    implicit none
    double precision :: Dfow
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    if ( alpha == 1 ) then
      if ( beta == 1 ) then
        Dfow = Dppfow(nth,np+si,nr,nsa)*pg(np+si,nsa)**2
      else if ( beta == 2 ) then
        Dfow = Dptfow(nth+sj,np+si,nr,nsa)*pg(np+si,nsa)
      else if ( beta == 3 ) then
        Dfow = Dprfow(nth,np+si,nr+sj,nsa)*pg(np+si,nsa)**2
      end if

    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        Dfow = Dtpfow(nth+si,np+sj,nr,nsa)*sin(thetamg(nth+si,np+sj,nr,nsa))
      else if ( beta == 2 ) then
        Dfow = Dttfow(nth+si,np,nr,nsa)*sin(thetamg(nth+si,np+sj,nr,nsa))/pm(np,nsa)
      else if ( beta == 3 ) then
        Dfow = Dtrfow(nth+si,np,nr+sj,nsa)*sin(thetamg(nth+si,np,nr+sj,nsa))
      end if

    else if ( alpha == 3 ) then
      if ( beta == 1 ) then
        Dfow = Drpfow(nth,np+sj,nr+si,nsa)*psimg(nr+si)
      else if ( beta == 2 ) then
        Dfow = Drtfow(nth+sj,np,nr+si,nsa)*psimg(nr+si)/pm(np,nsa)
      else if ( beta == 3 ) then
        Dfow = Drrfow(nth,np,nr+si,nsa)*psimg(nr+si)
      end if

    end if

  end function

  function w(sign,alpha,beta,si,sj,loc)
    ! weight function for alpha, beta, gamma
    implicit none
    double precision :: w
    integer,intent(in) :: sign,alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)
    
    if ( alpha == 1 ) then
      if ( beta == 2 ) then
        w = WEIGHP(nth+sj, np+si, nr, nsa)
      else if ( beta == 3 ) then
        w = WEIGHP(nth, np+si, nr+sj, nsa)
      else
        w = WEIGHP(nth, np+si, nr, nsa)
      end if

    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        w = WEIGHT(nth+si, np+sj, nr, nsa)
      else if ( beta == 3 ) then
        w = WEIGHT(nth+si, np, nr+sj, nsa)
      else 
        w = WEIGHT(nth+si, np, nr, nsa)
      end if

    else
      if ( beta == 1 ) then
        w = WEIGHR(nth, np+sj, nr+si, nsa)
      else if ( beta == 2 ) then
        w = WEIGHR(nth+sj, np, nr+si, nsa)
      else
        w = WEIGHR(nth, np, nr+si, nsa)
      end if
    end if

    if ( sign < 0 ) w = 1.d0-w

  end function w

  function check_external_boundary(alpha,beta,si,sj,loc) result(flag)
    ! flag == 0 -> inside external boundary
    !         1 -> outside external boundary
    implicit none
    integer :: flag
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    flag = 0

    ! check external boundary
    if ( alpha == 1 ) then
      if ( beta == 2 ) then
        if ( (np+si <= 0 .or. npmax+1 <= np+si)&
            .or. (nth+sj <= 0 .or. nthmax+1 <= nth+sj) )&
            flag = 1

      else if ( beta == 3 ) then
        if ( (np+si <= 0 .or. npmax+1 <= np+si)&
            .or. (nr+sj <= 0 .or. nrmax+1 <= nr+sj) )&
          flag = 1

      else ! beta == 0
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

      else ! beta == 0
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

      else ! beta == 0
        if ( nr+si <= 0 .or. nrmax+1 <= nr+si ) flag = 1

      end if
    end if

    ! check internal boundary
    

  end function check_external_boundary

  function get_nma(alpha,beta,si,sj,loc) result(n)
    implicit none
    integer :: n
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
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
      else ! beta == 0
        n = nma(nth,np+si,nr)
      end if
    
    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        n = nma(nth+si,np+sj,nr)
      else if ( beta == 3 ) then
        n = nma(nth+si,np,nr+sj)
      else ! beta == 0
        n = nma(nth+si,np,nr)
      end if
    
    else if ( alpha == 3 ) then
      if ( beta == 1 ) then
        n = nma(nth,np+sj,nr+si)
      else if ( beta == 2 ) then
        n = nma(nth+sj,np,nr+si)
      else ! beta == 0
        n = nma(nth,np,nr+si)
      end if

    end if

  end function

END MODULE fowexec
